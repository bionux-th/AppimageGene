//! 文本格式导入：**FASTA / GenBank → `.dna` 字节**。
//!
//! `.dna` 是 SnapGene 的原生格式，但真实工作流里序列大多来自 NCBI 的
//! FASTA / GenBank 下载。SnapGene 的「Open Files…」本身就能直接打开这些文件
//! （打开对话框的 `DNA`/`RNA`/`Protein` 过滤器包含它们），导入后成为
//! **未命名文档**，第一次保存时才写成 `.dna`。这里实现同一条路径。
//!
//! 支持子集（够覆盖下载文件）：
//!
//! - **FASTA**：多记录、标题行 `>名字 描述`、序列跨行；空白与非字母字符跳过。
//! - **GenBank**：`LOCUS`（名称/长度/拓扑/分子类型）、`DEFINITION`、`FEATURES`
//!   的位置表达式 `n`、`n..m`、`<n..m`、`n..>m`、`complement(...)`、`join(...)`、
//!   `order(...)`（含 `complement(join(...))` 嵌套）；限定词 `/label`、`/gene`、
//!   `/product`、`/note`、`/gene_synonym`。多记录文件逐条返回。
//!
//! 不做的事（诚实标注，见 `docs/snapgene-gap-analysis-*.zh-CN.md`）：
//! `/translation=` 不写入（SnapGene 由序列推翻译）、REFERENCE/COMMENT 不入库、
//! EMBL/AB1/SCF 未实现、单个特征的重叠外显子按原样分段但不再合并。

use crate::write::{Feature, NewFile};
use hx_core::{MoleculeFlags, MoleculeKind, Topology};

/// 导入的一条记录（一条 FASTA 记录 / 一个 GenBank LOCUS）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Imported {
    /// 记录名（FASTA 标题或 LOCUS 名）。
    pub name: String,
    /// 描述（FASTA 标题的剩余部分 / GenBank DEFINITION）。
    pub description: Option<String>,
    pub kind: MoleculeKind,
    pub circular: bool,
    /// 单链（GenBank `ss-DNA`）。
    pub single_stranded: bool,
    pub sequence: Vec<u8>,
    pub features: Vec<ImportedFeature>,
}

/// 导入的一条特征：`segments` 为 1-based 闭区间，按基因组顺序排列。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImportedFeature {
    pub name: String,
    /// SnapGene 特征类型（未知的 GenBank key 归入 `misc_feature`）。
    pub kind: String,
    pub segments: Vec<(usize, usize)>,
    /// 反向链（GenBank `complement`）。
    pub reverse: bool,
}

/// 可识别的文本格式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportFormat {
    Fasta,
    GenBank,
}

/// 按内容嗅探格式（不依赖扩展名：下载的文件常被改名）。
pub fn detect_format(bytes: &[u8]) -> Option<ImportFormat> {
    let text = std::str::from_utf8(bytes).ok()?;
    let head = text.trim_start();
    if head.starts_with('>') {
        return Some(ImportFormat::Fasta);
    }
    // GenBank：以 `LOCUS` 开头（允许前导空白/BOM）
    let head = head.trim_start_matches('\u{feff}');
    if head.starts_with("LOCUS") || head.contains("\nLOCUS ") || head.contains("\nORIGIN") {
        return Some(ImportFormat::GenBank);
    }
    None
}

/// 导入全部记录（FASTA 多记录 / GenBank 多 LOCUS）。
pub fn parse(bytes: &[u8]) -> Result<Vec<Imported>, String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "不是 UTF-8 文本文件".to_string())?;
    match detect_format(bytes) {
        Some(ImportFormat::Fasta) => Ok(parse_fasta(text)),
        Some(ImportFormat::GenBank) => parse_genbank(text),
        None => Err("unrecognized sequence file (expected FASTA or GenBank)".to_string()),
    }
}

// ---------------------------------------------------------------- FASTA

/// 解析 FASTA：`>` 起新记录，其余行拼序列。
pub fn parse_fasta(text: &str) -> Vec<Imported> {
    let mut out: Vec<Imported> = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(rest) = line.strip_prefix('>') {
            let mut it = rest.splitn(2, char::is_whitespace);
            let name = it.next().unwrap_or("").trim().to_string();
            let desc = it
                .next()
                .map(|d| d.trim().to_string())
                .filter(|d| !d.is_empty());
            out.push(Imported {
                name: if name.is_empty() {
                    "sequence".into()
                } else {
                    name
                },
                description: desc,
                kind: MoleculeKind::Dna,
                circular: false,
                single_stranded: false,
                sequence: Vec::new(),
                features: Vec::new(),
            });
            continue;
        }
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        let Some(rec) = out.last_mut() else {
            continue;
        };
        for b in line.bytes() {
            if b.is_ascii_alphabetic() || b == b'*' {
                rec.sequence.push(b.to_ascii_uppercase());
            }
        }
    }
    for rec in &mut out {
        rec.kind = guess_kind(&rec.sequence);
        if rec.kind == MoleculeKind::Protein {
            // 蛋白序列里的终止符不写入（SnapGene 的蛋白序列不含 `*`）
            rec.sequence.retain(|b| *b != b'*');
        }
    }
    out
}

/// 按字母表猜分子类型：出现 E/F/I/L/P/Q 等蛋白专有字母 ⇒ 蛋白；出现 U ⇒ RNA；
/// 只有 A/C/G/T/N 及 IUPAC 简并码 ⇒ DNA。
pub fn guess_kind(seq: &[u8]) -> MoleculeKind {
    let mut has_u = false;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'U' => has_u = true,
            b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' | b'Z' | b'J' | b'O' | b'X' | b'*' => {
                return MoleculeKind::Protein
            }
            _ => {}
        }
    }
    if has_u {
        MoleculeKind::Rna
    } else {
        MoleculeKind::Dna
    }
}

// ---------------------------------------------------------------- GenBank

/// 解析 GenBank（可含多条 `LOCUS`…`//` 记录）。
pub fn parse_genbank(text: &str) -> Result<Vec<Imported>, String> {
    let mut out = Vec::new();
    // 记录按 `//` 切分（`//` 只在行首出现时才是记录结束）
    let mut cur = String::new();
    let mut records: Vec<String> = Vec::new();
    for line in text.lines() {
        cur.push_str(line);
        cur.push('\n');
        if line.trim_end() == "//" {
            records.push(std::mem::take(&mut cur));
        }
    }
    if cur.trim().len() > 0 {
        records.push(cur);
    }
    for rec in records {
        if !rec.contains("LOCUS") {
            continue;
        }
        if let Some(r) = parse_genbank_record(&rec) {
            out.push(r);
        }
    }
    if out.is_empty() {
        return Err("GenBank 文件里没有可导入的 LOCUS 记录".to_string());
    }
    Ok(out)
}

/// 单条 GenBank 记录。
fn parse_genbank_record(rec: &str) -> Option<Imported> {
    let mut name = String::new();
    let mut kind = MoleculeKind::Dna;
    let mut circular = false;
    let mut single_stranded = false;
    let mut description: Vec<String> = Vec::new();
    let mut sequence: Vec<u8> = Vec::new();
    let mut features: Vec<ImportedFeature> = Vec::new();

    /// 解析状态：头部 / FEATURES / ORIGIN。
    #[derive(PartialEq)]
    enum Mode {
        Head,
        Definition,
        Features,
        Origin,
    }
    let mut mode = Mode::Head;
    // 当前正在收集的特征
    let mut feat_key = String::new();
    let mut feat_loc = String::new();
    let mut feat_quals: Vec<(String, String)> = Vec::new(); // 值保留引号，flush 时去引号
    let mut open_qual: Option<usize> = None; // 引号未闭合的限定词下标

    fn flush(
        features: &mut Vec<ImportedFeature>,
        key: &mut String,
        loc: &mut String,
        quals: &mut Vec<(String, String)>,
    ) {
        if !key.is_empty() && !loc.is_empty() {
            let (segments, reverse) = parse_location(&loc.replace(char::is_whitespace, ""));
            if let Some(f) = build_feature(key, quals, &segments, reverse) {
                features.push(f);
            }
        }
        key.clear();
        loc.clear();
        quals.clear();
    }

    for line in rec.lines() {
        let (tag, rest) = split_gb_line(line);
        match tag {
            "LOCUS" => {
                mode = Mode::Head;
                let mut it = rest.split_whitespace();
                name = it.next().unwrap_or("sequence").to_string();
                for tok in rest.split_whitespace() {
                    let t = tok.to_ascii_lowercase();
                    if t == "circular" {
                        circular = true;
                    }
                    if t.starts_with("ss-") || t == "ssdna" {
                        single_stranded = true;
                    }
                    kind = match t.as_str() {
                        "dna" | "ss-dna" | "ds-dna" | "ssdna" | "dsdna" => MoleculeKind::Dna,
                        "rna" | "mrna" | "trna" | "rrna" | "snrna" | "ncrna" | "ss-rna" => {
                            MoleculeKind::Rna
                        }
                        "protein" | "aa" | "peptide" => MoleculeKind::Protein,
                        _ => kind,
                    };
                }
            }
            "DEFINITION" => {
                mode = Mode::Definition;
                if !rest.trim().is_empty() {
                    description.push(rest.trim().to_string());
                }
            }
            "ORIGIN" => {
                mode = Mode::Origin;
            }
            "FEATURES" => {
                mode = Mode::Features;
            }
            _ => match mode {
                Mode::Definition => {
                    let t = line.trim();
                    // 下一段（如 ACCESSION）开始 ⇒ 描述结束
                    if !t.is_empty() && line.starts_with("  ") && !line.starts_with("   ") {
                        if !is_gb_section(t) {
                            description.push(t.to_string());
                        }
                    }
                }
                Mode::Features => {
                    let t = line.trim_end();
                    if t.starts_with(' ') && !t.trim().is_empty() && tag.is_empty() {
                        let t = t.trim();
                        if let Some(q) = t.strip_prefix('/') {
                            let (k, v) = match q.split_once('=') {
                                Some((k, v)) => (k.to_string(), v.trim().to_string()),
                                None => (q.to_string(), String::new()),
                            };
                            open_qual = (v.starts_with('"') && !v.ends_with('"'))
                                .then_some(feat_quals.len());
                            feat_quals.push((k, v));
                        } else if let Some(i) = open_qual {
                            let v = &mut feat_quals[i].1;
                            v.push(' ');
                            v.push_str(t);
                            if t.ends_with('"') {
                                open_qual = None;
                            }
                        } else {
                            // 位置表达式的续行（`join(...)` 跨行）
                            feat_loc.push_str(t);
                        }
                    }
                }
                Mode::Origin => {
                    for b in line.bytes() {
                        if b.is_ascii_alphabetic() {
                            sequence.push(b.to_ascii_uppercase());
                        }
                    }
                }
                Mode::Head => {}
            },
        }

        // FEATURES 段里 `tag` 非空且缩进 5 空格 = 新的特征键行
        if mode == Mode::Features && !tag.is_empty() && tag != "FEATURES" {
            flush(&mut features, &mut feat_key, &mut feat_loc, &mut feat_quals);
            open_qual = None;
            feat_key = tag.to_string();
            feat_loc = rest.trim().to_string();
        }
    }
    flush(&mut features, &mut feat_key, &mut feat_loc, &mut feat_quals);

    if sequence.is_empty() {
        return None;
    }
    Some(Imported {
        name: if name.is_empty() {
            "sequence".into()
        } else {
            name
        },
        description: if description.is_empty() {
            None
        } else {
            Some(description.join(" "))
        },
        kind,
        circular,
        single_stranded,
        sequence,
        features,
    })
}

/// GenBank 的段落关键词（用于判断 DEFINITION 续行何时结束）。
fn is_gb_section(t: &str) -> bool {
    matches!(
        t.split_whitespace().next().unwrap_or(""),
        "ACCESSION"
            | "VERSION"
            | "KEYWORDS"
            | "SOURCE"
            | "REFERENCE"
            | "AUTHORS"
            | "TITLE"
            | "JOURNAL"
            | "PUBMED"
            | "COMMENT"
            | "PRIMARY"
            | "FEATURES"
            | "ORIGIN"
            | "BASE"
            | "CONTIG"
            | "DBSOURCE"
    )
}

/// 拆一行 GenBank：`tag` = 行首关键词（去掉两端空白），`rest` = 其后内容。
fn split_gb_line(line: &str) -> (&str, &str) {
    let trimmed = line.trim_start();
    let lead = line.len() - trimmed.len();
    // 前 12 列是关键词区（`FEATURES`/`ORIGIN` 等），之后可能是续行数据
    let (tag, rest) = if lead <= 12 && !trimmed.is_empty() {
        match trimmed.split_once(char::is_whitespace) {
            Some((k, r)) => (k, r),
            None => (trimmed, ""),
        }
    } else {
        ("", trimmed)
    };
    (tag, rest)
}

fn unquote(s: &str) -> String {
    s.trim().trim_matches('"').trim().to_string()
}

/// 由 GenBank key + 限定词构造特征。
fn build_feature(
    key: &str,
    quals: &[(String, String)],
    segments: &[(usize, usize)],
    reverse: bool,
) -> Option<ImportedFeature> {
    if segments.is_empty() {
        return None;
    }
    let q = |name: &str| -> Option<String> {
        quals
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| unquote(v))
            .filter(|v| !v.is_empty())
    };
    let name = q("label")
        .or_else(|| q("gene"))
        .or_else(|| q("product"))
        .or_else(|| q("note"))
        .or_else(|| q("gene_synonym"))
        .unwrap_or_else(|| key.to_string())
        .trim()
        .to_string();
    let name = if name.is_empty() {
        key.to_string()
    } else {
        name
    };
    Some(ImportedFeature {
        name,
        kind: snapgene_type(key).to_string(),
        segments: segments.to_vec(),
        reverse,
    })
}

/// GenBank key → SnapGene 特征类型。未知的键归 `misc_feature`
/// （SnapGene 的特征类型词表见其 `standardCommonFeatures.ftrs`）。
pub fn snapgene_type(key: &str) -> &'static str {
    const KNOWN: [&str; 44] = [
        "CDS",
        "promoter",
        "terminator",
        "misc_feature",
        "rep_origin",
        "primer_bind",
        "protein_bind",
        "RBS",
        "polyA_signal",
        "polyA_site",
        "enhancer",
        "intron",
        "exon",
        "gene",
        "LTR",
        "misc_RNA",
        "ncRNA",
        "mRNA",
        "rRNA",
        "tRNA",
        "5'UTR",
        "3'UTR",
        "5'clip",
        "3'clip",
        "regulatory",
        "repeat_region",
        "mobile_element",
        "misc_recomb",
        "misc_signal",
        "sig_peptide",
        "propeptide",
        "mat_peptide",
        "transit_peptide",
        "oriT",
        "source",
        "variation",
        "unsure",
        "conserved",
        "STS",
        "telomere",
        "centromere",
        "D-loop",
        "attenuator",
        "protein",
    ];
    KNOWN
        .iter()
        .find(|k| k.eq_ignore_ascii_case(key))
        .copied()
        .unwrap_or("misc_feature")
}

/// SnapGene 的特征类型默认色（取自其 `standardCommonFeatures.ftrs` 与样例文件；
/// 未知类型用 `misc_feature` 的灰蓝）。
pub fn feature_color(kind: &str) -> &'static str {
    match kind {
        "CDS" => "#ccffcc",
        "promoter" | "terminator" | "enhancer" | "regulatory" | "source" => "#ffffff",
        "misc_feature" => "#99ccff",
        "rep_origin" => "#ffff00",
        "primer_bind" => "#a020f0",
        "protein_bind" => "#31849b",
        "repeat_region" | "RBS" | "mRNA" | "polyA_signal" | "polyA_site" | "5'UTR"
        | "3'UTR" | "oriT" => "#a6acb3",
        "misc_RNA" => "#00ccff",
        "ncRNA" => "#8fbc8f",
        "gene" => "#ff7f50",
        "tRNA" => "#228b22",
        "exon" | "intron" => "#708090",
        "LTR" | "misc_recomb" | "mobile_element" | "misc_signal" => "#ffe4c4",
        "sig_peptide" | "propeptide" | "mat_peptide" | "transit_peptide" => "#993366",
        _ => "#a6acb3",
    }
}

// ---- 位置表达式 ----------------------------------------------------------

/// 解析 GenBank 位置表达式 → （1-based 闭区间列表，是否反向）。
/// 支持 `complement(...)`、`join(...)`、`order(...)` 与 `<`/`>` 模糊端点。
pub fn parse_location(s: &str) -> (Vec<(usize, usize)>, bool) {
    let s = s.trim();
    if let Some(inner) = strip_call(s, "complement") {
        let (segs, _) = parse_location(inner);
        return (segs, true);
    }
    for call in ["join", "order"] {
        if let Some(inner) = strip_call(s, call) {
            let mut segs = Vec::new();
            for part in split_top(inner) {
                let (mut sub, rev) = parse_location(&part);
                if rev {
                    sub.reverse();
                }
                segs.extend(sub);
            }
            segs.sort_by_key(|(a, _)| *a);
            return (segs, false);
        }
    }
    // 简单位置：`n`、`n..m`、`<n..m`、`n..>m`
    let part = s.trim_start_matches(['<', '>']);
    let (a, b) = match part.split_once("..") {
        Some((a, b)) => (
            a.trim_start_matches(['<', '>']).trim(),
            b.trim_start_matches(['<', '>']).trim(),
        ),
        None => (part, part),
    };
    let (Ok(a), Ok(b)) = (a.parse::<usize>(), b.parse::<usize>()) else {
        return (Vec::new(), false);
    };
    if a == 0 || b == 0 {
        return (Vec::new(), false);
    }
    if a <= b {
        (vec![(a, b)], false)
    } else {
        (vec![(b, a)], false)
    }
}

/// `name(...)` → 括号内内容（要求整串就是这个调用）。
fn strip_call<'a>(s: &'a str, name: &str) -> Option<&'a str> {
    let s = s.trim();
    if s.len() <= name.len() + 2 || !s[..name.len()].eq_ignore_ascii_case(name) {
        return None;
    }
    let rest = s[name.len()..].trim_start();
    let inner = rest.strip_prefix('(')?.strip_suffix(')')?;
    Some(inner)
}

/// 按顶层逗号切分（忽略括号内的逗号）。
fn split_top(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut cur = String::new();
    for c in s.chars() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth == 0 => {
                out.push(std::mem::take(&mut cur));
                continue;
            }
            _ => {}
        }
        cur.push(c);
    }
    if !cur.trim().is_empty() {
        out.push(cur);
    }
    out
}

// ---- 写 .dna ------------------------------------------------------------

/// 导入记录 → `.dna` 字节（未命名文档；第一次保存时才落地）。
///
/// 拓扑/链数按导入信息设置；特征写进 block 10，颜色按类型给出 SnapGene 的默认色。
pub fn to_dna_bytes(rec: &Imported) -> Vec<u8> {
    let len = rec.sequence.len();
    let mut flags = MoleculeFlags::from_byte(0);
    if rec.kind.is_nucleotide() {
        flags = flags
            .with_topology(if rec.circular {
                Topology::Circular
            } else {
                Topology::Linear
            })
            .with_strandedness(if rec.single_stranded {
                hx_core::Strandedness::Single
            } else {
                hx_core::Strandedness::Double
            });
    }
    let feats: Vec<Feature> = rec
        .features
        .iter()
        .filter_map(|f| {
            let segs: Vec<(usize, usize)> = f
                .segments
                .iter()
                .map(|(a, b)| {
                    // 越界坐标夹到序列内（下载文件里偶见）
                    (a.saturating_sub(1).min(len) + 1, (*b).min(len))
                })
                .filter(|(a, b)| a <= b)
                .collect();
            let (a, b) = *segs.first()?;
            Some(Feature {
                start: a,
                end: b,
                name: f.name.clone(),
                kind: f.kind.clone(),
                color: feature_color(&f.kind).to_string(),
                directionality: if !rec.kind.features_have_direction() {
                    0
                } else if f.reverse {
                    2
                } else {
                    1
                },
                extra_segments: segs.into_iter().skip(1).collect(),
            })
        })
        .collect();
    let nf = NewFile {
        kind: rec.kind,
        flags,
        sequence: rec.sequence.clone(),
        name: Some(rec.name.clone()),
        description: rec.description.clone(),
        created_by: Some("AppimageGene (imported)".to_string()),
        features: feats,
        extra_blocks: Vec::new(),
    };
    nf.to_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA: &str = ">seq1 test record\nACGTacgt\nAAAA\n>seq2\nMKV*L\n";

    #[test]
    fn fasta_two_records() {
        let recs = parse(FASTA.as_bytes()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "seq1");
        assert_eq!(recs[0].description.as_deref(), Some("test record"));
        assert_eq!(recs[0].sequence, b"ACGTACGTAAAA");
        assert_eq!(recs[0].kind, MoleculeKind::Dna);
        assert_eq!(recs[1].kind, MoleculeKind::Protein);
        // 蛋白终止符不写入
        assert_eq!(recs[1].sequence, b"MKVL");
    }

    #[test]
    fn fasta_rna_detected_by_u() {
        let recs = parse(b">r\nAUGCUU\n").unwrap();
        assert_eq!(recs[0].kind, MoleculeKind::Rna);
    }

    #[test]
    fn guess_kind_alphabet() {
        assert_eq!(guess_kind(b"ACGTN"), MoleculeKind::Dna);
        assert_eq!(guess_kind(b"ACGU"), MoleculeKind::Rna);
        assert_eq!(guess_kind(b"MPEPTIDE"), MoleculeKind::Protein);
    }

    #[test]
    fn location_simple_and_complement() {
        assert_eq!(parse_location("123..456"), (vec![(123, 456)], false));
        assert_eq!(parse_location("complement(10..20)"), (vec![(10, 20)], true));
        assert_eq!(parse_location("42"), (vec![(42, 42)], false));
        assert_eq!(parse_location("<1..>99"), (vec![(1, 99)], false));
    }

    #[test]
    fn location_join_and_nested_complement() {
        assert_eq!(
            parse_location("join(1..10,20..30)"),
            (vec![(1, 10), (20, 30)], false)
        );
        assert_eq!(
            parse_location("complement(join(20..30,1..10))"),
            (vec![(1, 10), (20, 30)], true)
        );
    }

    const GB: &str = "\
LOCUS       TestSeq                 20 bp    DNA     circular SYN 01-JAN-2020
DEFINITION  A test record for import.
ACCESSION   ABC123
FEATURES             Location/Qualifiers
     source          1..20
                     /organism=\"synthetic\"
     CDS             complement(join(5..9,12..16))
                     /gene=\"gfp\"
                     /product=\"green fluorescent protein\"
     promoter        1..4
                     /label=\"T7 promoter\"
ORIGIN
        1 atgcatgcat gcatgcatgc
//
";

    #[test]
    fn genbank_record_header() {
        let recs = parse(GB.as_bytes()).unwrap();
        assert_eq!(recs.len(), 1);
        let r = &recs[0];
        assert_eq!(r.name, "TestSeq");
        assert!(r.circular);
        assert_eq!(r.kind, MoleculeKind::Dna);
        assert_eq!(r.sequence, b"ATGCATGCATGCATGCATGC");
        assert_eq!(r.features.len(), 3);
        let cds = r.features.iter().find(|f| f.kind == "CDS").unwrap();
        assert!(cds.reverse);
        assert_eq!(cds.name, "gfp");
        assert_eq!(cds.segments, vec![(5, 9), (12, 16)]);
        let prom = r.features.iter().find(|f| f.kind == "promoter").unwrap();
        assert_eq!(prom.name, "T7 promoter");
    }

    #[test]
    fn genbank_to_dna_round_trips_through_parser() {
        let recs = parse(GB.as_bytes()).unwrap();
        let bytes = to_dna_bytes(&recs[0]);
        let file = crate::DnaFile::parse(&bytes).unwrap();
        assert_eq!(file.header.type_of_sequence, 1);
        let xml = file.block(10).map(|b| b.payload.clone()).unwrap_or_default();
        let feats = crate::annot::parse_features(&xml);
        assert_eq!(feats.len(), 3);
        let cds = feats.iter().find(|f| f.kind == "CDS").unwrap();
        assert_eq!(cds.directionality, 2);
        assert_eq!(cds.segments.len(), 2);
        assert_eq!(cds.segments[0].range.start, 4);
        assert_eq!(cds.segments[1].range.end, 16);
        assert_eq!(feats.iter().filter(|f| f.kind == "promoter").count(), 1);
    }

    #[test]
    fn genbank_multiple_loci() {
        let two = format!("{GB}{GB}");
        let recs = parse(two.as_bytes()).unwrap();
        assert_eq!(recs.len(), 2);
    }

    #[test]
    fn feature_colors_come_from_type_table() {
        assert_eq!(feature_color("CDS"), "#ccffcc");
        assert_eq!(feature_color("promoter"), "#ffffff");
        assert_eq!(feature_color("whatever"), "#a6acb3");
    }
}
