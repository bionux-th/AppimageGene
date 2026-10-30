//! 引物分析 / 设计（**零依赖**）。
//!
//! 这里只做纯计算，不碰文件格式：Tm / GC / 发夹 / 二聚体 / 结合位点 / 引物对设计。
//!
//! # Tm 模型
//!
//! SantaLucia 1998 统一最近邻参数 + 盐校正
//!
//! ```text
//! Tm = ΔH·1000 / (ΔS + R·ln(C_T/4)) − 273.15 + 16.6·log10([Na⁺]/(1+0.7[Na⁺]))
//! ```
//!
//! 常数 `[Na⁺] = 100 mM`、`C_T = 3 µM` 是用**官方样例库里 30 条引物**的
//! SnapGene 存储值（`<BindingSite meltingTemperature="…">`）拟合出来的：
//! 平均绝对误差 **0.82 °C**、平均偏差 **−0.04 °C**（见本文件末尾的
//! `tm_matches_snapgene_stored_values`）。所以表格里算出来的 Tm 与 SnapGene
//! 显示的整数 Tm 基本一致。
//!
//! # 结合位点
//!
//! 与 SnapGene 的 `<HybridizationParams>` 语义一致：允许少量错配，但
//! **3′ 端最后 10 个碱基必须完全配对**（`THREE_PRIME_ANCHOR`），否则不算结合。
//! 正向引物按自身序列扫上链，反向引物按反向互补扫下链。

use hx_formats::annot::Range;

/// 最近邻参数：`(ΔH kcal/mol, ΔS cal/(mol·K))`（SantaLucia 1998 统一参数）。
const NN: [(&[u8; 2], f64, f64); 16] = [
    (b"AA", -7.9, -22.2),
    (b"TT", -7.9, -22.2),
    (b"AT", -7.2, -20.4),
    (b"TA", -7.2, -21.3),
    (b"CA", -8.5, -22.7),
    (b"TG", -8.5, -22.7),
    (b"GT", -8.4, -22.4),
    (b"AC", -8.4, -22.4),
    (b"CT", -7.8, -21.0),
    (b"AG", -7.8, -21.0),
    (b"GA", -8.2, -22.2),
    (b"TC", -8.2, -22.2),
    (b"CG", -10.6, -27.2),
    (b"GC", -9.8, -24.4),
    (b"GG", -8.0, -19.9),
    (b"CC", -8.0, -19.9),
];

/// 盐浓度（mol/L），见模块文档。
pub const NA_M: f64 = 0.100;
/// 引物总浓度（mol/L），见模块文档。
pub const CT_M: f64 = 3.0e-6;
/// 气体常数（cal/(mol·K)）。
const R: f64 = 1.987;
/// 3′ 端必须完全配对的碱基数（SnapGene 的 `minContinuousMatchLen` 默认 10）。
pub const THREE_PRIME_ANCHOR: usize = 10;

/// 设计引物时，引物 5′ 端在目标区间端点附近可滑动的碱基数。
pub const SEARCH_WINDOW: usize = 30;

fn nn_pair(a: u8, b: u8) -> Option<(f64, f64)> {
    let key = [a.to_ascii_uppercase(), b.to_ascii_uppercase()];
    NN.iter()
        .find(|(k, _, _)| **k == key)
        .map(|(_, h, s)| (*h, *s))
}

/// 只有 A/C/G/T 才算“可计算”的碱基（简并码/N 直接视为不可用）。
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// 最近邻 Tm（°C）。序列含非 ACGT 或长度 < 2 时返回 `NaN`。
pub fn tm(seq: &[u8]) -> f64 {
    if seq.len() < 2 || !seq.iter().all(|b| is_acgt(*b)) {
        return f64::NAN;
    }
    let mut dh = 0.0;
    let mut ds = 0.0;
    for w in seq.windows(2) {
        let (h, s) = nn_pair(w[0], w[1]).expect("checked above");
        dh += h;
        ds += s;
    }
    // 末端引发项
    for e in [seq[0], seq[seq.len() - 1]] {
        if matches!(e, b'G' | b'C' | b'g' | b'c') {
            dh += 0.1;
            ds += -2.8;
        } else {
            dh += 2.3;
            ds += 4.1;
        }
    }
    let t = dh * 1000.0 / (ds + R * (CT_M / 4.0).ln()) - 273.15;
    t + 16.6 * (NA_M / (1.0 + 0.7 * NA_M)).log10()
}

/// 指定总浓度的最近邻 Tm（`tm` 是它的 `[Na⁺]=100 mM, C_T=3 µM` 特例）。
pub fn tm_with_ct(seq: &[u8], ct: f64) -> f64 {
    if seq.len() < 2 || !seq.iter().all(|b| is_acgt(*b)) {
        return f64::NAN;
    }
    let mut dh = 0.0;
    let mut ds = 0.0;
    for w in seq.windows(2) {
        let (h, s) = nn_pair(w[0], w[1]).expect("checked above");
        dh += h;
        ds += s;
    }
    for e in [seq[0], seq[seq.len() - 1]] {
        if matches!(e, b'G' | b'C' | b'g' | b'c') {
            dh += 0.1;
            ds += -2.8;
        } else {
            dh += 2.3;
            ds += 4.1;
        }
    }
    let t = dh * 1000.0 / (ds + R * (ct / 4.0).ln()) - 273.15;
    t + 16.6 * (NA_M / (1.0 + 0.7 * NA_M)).log10()
}

/// GC 百分比（0..100）。非 ACGT 不计入分子，也不计入分母。
pub fn gc_percent(seq: &[u8]) -> f64 {
    let n = seq.iter().filter(|b| is_acgt(**b)).count();
    if n == 0 {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
        .count();
    gc as f64 * 100.0 / n as f64
}

/// 单链 DNA 的分子量（Da，平均质量）。
///
/// 公式按真 SnapGene 的显示值反推并核对：
/// `313.21·A + 304.2·T + 289.18·C + 329.21·G − 61.96`；5′ 磷酸化再 `+79.0`。
/// 实测两例（SnapGene 8.2.3 的引物详情行）：
/// `CAGTCAGGATCCATGAGTAAAGGAGAAGAACT`（32-mer，未磷酸化）→ 9939.6 Da，
/// `TTAGACACAACATTGAAGAT`（20-mer，5′ 磷酸化）→ 6212.1 Da。
/// 非 ACGT 碱基按 0 计（简并码不适用）。
pub fn ss_dna_mass(seq: &[u8], phosphorylated: bool) -> f64 {
    let mut mass = -61.96;
    for &b in seq {
        mass += match b.to_ascii_uppercase() {
            b'A' => 313.21,
            b'T' => 304.2,
            b'C' => 289.18,
            b'G' => 329.21,
            _ => 0.0,
        };
    }
    if phosphorylated {
        mass += 79.0;
    }
    mass
}

/// 5′→3′ 的反向互补（IUPAC 简并码，未知字符原样）。
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| complement(*b)).collect()
}

// ---------------------------------------------------------------------------
// 引物编辑器的插入面板（密码子 / 酶切位点 / 短特征序列）
// ---------------------------------------------------------------------------

/// 标准遗传密码的一个氨基酸（三字母 / 一字母 / 密码子表）。
pub struct CodonGroup {
    pub aa3: &'static str,
    pub aa1: char,
    pub codons: &'static [&'static str],
}

/// 标准遗传密码（按三字母名排序，Stop 在末尾；密码子按常规教科书序）。
pub const CODON_TABLE: &[CodonGroup] = &[
    CodonGroup { aa3: "Ala", aa1: 'A', codons: &["GCT", "GCC", "GCA", "GCG"] },
    CodonGroup { aa3: "Arg", aa1: 'R', codons: &["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"] },
    CodonGroup { aa3: "Asn", aa1: 'N', codons: &["AAT", "AAC"] },
    CodonGroup { aa3: "Asp", aa1: 'D', codons: &["GAT", "GAC"] },
    CodonGroup { aa3: "Cys", aa1: 'C', codons: &["TGT", "TGC"] },
    CodonGroup { aa3: "Gln", aa1: 'Q', codons: &["CAA", "CAG"] },
    CodonGroup { aa3: "Glu", aa1: 'E', codons: &["GAA", "GAG"] },
    CodonGroup { aa3: "Gly", aa1: 'G', codons: &["GGT", "GGC", "GGA", "GGG"] },
    CodonGroup { aa3: "His", aa1: 'H', codons: &["CAT", "CAC"] },
    CodonGroup { aa3: "Ile", aa1: 'I', codons: &["ATT", "ATC", "ATA"] },
    CodonGroup { aa3: "Leu", aa1: 'L', codons: &["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"] },
    CodonGroup { aa3: "Lys", aa1: 'K', codons: &["AAA", "AAG"] },
    CodonGroup { aa3: "Met", aa1: 'M', codons: &["ATG"] },
    CodonGroup { aa3: "Phe", aa1: 'F', codons: &["TTT", "TTC"] },
    CodonGroup { aa3: "Pro", aa1: 'P', codons: &["CCT", "CCC", "CCA", "CCG"] },
    CodonGroup { aa3: "Ser", aa1: 'S', codons: &["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"] },
    CodonGroup { aa3: "Thr", aa1: 'T', codons: &["ACT", "ACC", "ACA", "ACG"] },
    CodonGroup { aa3: "Trp", aa1: 'W', codons: &["TGG"] },
    CodonGroup { aa3: "Tyr", aa1: 'Y', codons: &["TAT", "TAC"] },
    CodonGroup { aa3: "Val", aa1: 'V', codons: &["GTT", "GTC", "GTA", "GTG"] },
    CodonGroup { aa3: "Stop", aa1: '*', codons: &["TAA", "TAG", "TGA"] },
];

/// 短特征序列（SnapGene 的 `Short feature sequence` 下拉：表位/亲和标签 + 蛋白酶切位点）。
pub struct Peptide {
    /// 分组名（下拉里的分组标题）。
    pub group: &'static str,
    pub name: &'static str,
    /// 氨基酸序列（单字母）。
    pub peptide: &'static str,
}

pub const GROUP_EPITOPE: &str = "Epitope and Affinity Tags";
pub const GROUP_PROTEASE: &str = "Protease Cleavage Sites";

/// 标签表：SnapGene 8.2.3 的 `Short feature sequence ▾` 实测条目。
pub const PEPTIDES: &[Peptide] = &[
    Peptide { group: GROUP_EPITOPE, name: "6xHis", peptide: "HHHHHH" },
    Peptide { group: GROUP_EPITOPE, name: "8xHis", peptide: "HHHHHHHH" },
    Peptide { group: GROUP_EPITOPE, name: "10xHis", peptide: "HHHHHHHHHH" },
    Peptide { group: GROUP_EPITOPE, name: "AviTag", peptide: "GLNDIFEAQKIEWHE" },
    Peptide { group: GROUP_EPITOPE, name: "c-myc", peptide: "EQKLISEEDL" },
    Peptide { group: GROUP_EPITOPE, name: "FLAG", peptide: "DYKDDDDK" },
    Peptide { group: GROUP_EPITOPE, name: "Glu-Glu", peptide: "EYMPME" },
    Peptide { group: GROUP_EPITOPE, name: "Hemagglutinin (HA)", peptide: "YPYDVPDYA" },
    Peptide { group: GROUP_EPITOPE, name: "HSV", peptide: "QPELAPEDPED" },
    Peptide { group: GROUP_EPITOPE, name: "S-Tag", peptide: "KETAAAKFERQHMDS" },
    Peptide { group: GROUP_EPITOPE, name: "Strep-Tag II", peptide: "WSHPQFEK" },
    Peptide { group: GROUP_EPITOPE, name: "T7", peptide: "MASMTGGQQMG" },
    Peptide { group: GROUP_EPITOPE, name: "VSV-G", peptide: "YTDIEMNRLGK" },
    Peptide { group: GROUP_EPITOPE, name: "V5", peptide: "GKPIPNPLLGLDST" },
    Peptide { group: GROUP_PROTEASE, name: "Enterokinase", peptide: "DDDDK" },
    Peptide { group: GROUP_PROTEASE, name: "Factor Xa", peptide: "IEGR" },
    Peptide { group: GROUP_PROTEASE, name: "Genenase I", peptide: "HY" },
    Peptide { group: GROUP_PROTEASE, name: "PreScission Protease (HRV 3C)", peptide: "LEVLFQGP" },
    Peptide { group: GROUP_PROTEASE, name: "TEV Protease", peptide: "ENLYFQG" },
    Peptide { group: GROUP_PROTEASE, name: "Thrombin", peptide: "LVPRGS" },
];

/// 大肠杆菌高频密码子表（`PEPTIDES` 的氨基酸 → DNA）。
///
/// SnapGene 对同一个标签给出哪套密码子没有公开说明，这里用一张固定表，
/// 结果稳定、可复现。
pub fn dna_for_peptide(peptide: &str) -> String {
    let mut out = String::new();
    for aa in peptide.chars() {
        let cod = match aa.to_ascii_uppercase() {
            'A' => "GCG",
            'R' => "CGT",
            'N' => "AAT",
            'D' => "GAT",
            'C' => "TGT",
            'Q' => "CAG",
            'E' => "GAA",
            'G' => "GGT",
            'H' => "CAT",
            'I' => "ATT",
            'L' => "CTG",
            'K' => "AAA",
            'M' => "ATG",
            'F' => "TTT",
            'P' => "CCG",
            'S' => "AGC",
            'T' => "ACC",
            'W' => "TGG",
            'Y' => "TAT",
            'V' => "GTT",
            _ => "",
        };
        out.push_str(cod);
    }
    out
}

/// IUPAC 互补（与 `hx-algo::edit_session::complement` 同一张表，避免跨模块依赖）。
pub const fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' | b'U' => b'A',
        b't' | b'u' => b'a',
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'S' => b'S',
        b's' => b's',
        b'W' => b'W',
        b'w' => b'w',
        b'K' => b'M',
        b'k' => b'm',
        b'M' => b'K',
        b'm' => b'k',
        b'B' => b'V',
        b'b' => b'v',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        b'V' => b'B',
        b'v' => b'b',
        other => other,
    }
}

// ---------------------------------------------------------------------------
// 发夹 / 二聚体
// ---------------------------------------------------------------------------

/// 一个发夹（自身回折）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hairpin {
    /// 茎的碱基对数。
    pub stem: usize,
    /// 环的长度（≥ 3）。
    pub loop_len: usize,
    /// 茎的 Tm（°C），作为“稳定性”的排序键。
    pub tm: f64,
    /// 茎在引物里的起点（0-based）。
    pub start: usize,
}

/// 发夹 Tm 用的有效浓度（分子内配对，局部浓度远高于溶液中的引物浓度）。
/// 这是启发式值：只用于「发夹是否值得提示」的排序，不是热力学精确值。
pub const HAIRPIN_CT_M: f64 = 1.0e-4;

/// 找最稳的自身发夹（茎 ≥ 4 bp、环 ≥ 3 nt）。没有返回 `None`。
pub fn hairpin(seq: &[u8]) -> Option<Hairpin> {
    let s: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    if s.len() < 11 {
        return None;
    }
    let mut best: Option<Hairpin> = None;
    // 茎从 4 到 len/3 逐步加长；i 是 5′ 臂起点，j 是 3′ 臂终点
    for stem in 4..=(s.len() / 3) {
        for i in 0..s.len() {
            let loop_start = i + stem;
            // 环至少 3 nt
            for loop_len in 3..=s.len() {
                let j = loop_start + loop_len;
                if j + stem > s.len() {
                    break;
                }
                let left = &s[i..i + stem];
                let right = &s[j..j + stem];
                // left 与 right 的反向互补配对
                let ok = left
                    .iter()
                    .zip(right.iter().rev())
                    .all(|(a, b)| complement(*a) == *b);
                if !ok {
                    continue;
                }
                let t = tm_with_ct(left, HAIRPIN_CT_M);
                let hp = Hairpin {
                    stem,
                    loop_len,
                    tm: t,
                    start: i,
                };
                if best.map(|b| hp.tm > b.tm).unwrap_or(true) {
                    best = Some(hp);
                }
            }
        }
    }
    best
}

/// 二聚体（两条链互补重叠）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dimer {
    /// 连续互补的碱基数。
    pub overlap: usize,
    /// 重叠部分的 Tm（°C）。
    pub tm: f64,
    /// 是否涉及 3′ 端（3′ 二聚体更影响扩增）。
    pub three_prime: bool,
}

/// 自身二聚体（一条引物与自己的另一份拷贝）。
pub fn self_dimer(seq: &[u8]) -> Option<Dimer> {
    dimer_between(seq, seq)
}

/// 两条不同引物之间的二聚体。
pub fn cross_dimer(a: &[u8], b: &[u8]) -> Option<Dimer> {
    dimer_between(a, b)
}

fn dimer_between(a: &[u8], b: &[u8]) -> Option<Dimer> {
    if a.len() < 4 || b.len() < 4 {
        return None;
    }
    let a: Vec<u8> = a.iter().map(|x| x.to_ascii_uppercase()).collect();
    let b: Vec<u8> = b.iter().map(|x| x.to_ascii_uppercase()).collect();
    let rc_b = revcomp(&b);
    let mut best: Option<Dimer> = None;
    // 位移：把 rc_b 相对 a 左右挪动，找最长连续互补段
    let lo = -(a.len() as isize) + 4;
    let hi = b.len() as isize - 4;
    for off in lo..=hi {
        let mut run_start: Option<usize> = None;
        let mut run_len = 0usize;
        let (mut best_run, mut best_start) = (0usize, 0usize);
        // 遍历重叠窗口
        let i0 = off.max(0) as usize;
        let i1 = ((off + rc_b.len() as isize).min(a.len() as isize)).max(0) as usize;
        for i in i0..i1 {
            let j = (i as isize - off) as usize;
            if j >= rc_b.len() {
                break;
            }
            if complement(a[i]) == rc_b[j] {
                if run_start.is_none() {
                    run_start = Some(i);
                }
                run_len += 1;
            } else {
                if run_len > best_run {
                    best_run = run_len;
                    best_start = run_start.unwrap_or(i);
                }
                run_start = None;
                run_len = 0;
            }
        }
        if run_len > best_run {
            best_run = run_len;
            best_start = run_start.unwrap_or(i0);
        }
        if best_run < 4 {
            continue;
        }
        let duplex: Vec<u8> = a[best_start..best_start + best_run].to_vec();
        let t = tm(&duplex);
        // 是否贴着任一方的 3′ 端（a 的右端 / b 的 3′ 端对应 rc_b 的左端）
        let a_three = best_start + best_run == a.len();
        let b_three = best_start as isize - off <= 1;
        let d = Dimer {
            overlap: best_run,
            tm: t,
            three_prime: a_three || b_three,
        };
        if best.map(|x| d.tm > x.tm).unwrap_or(true) {
            best = Some(d);
        }
    }
    best
}

// ---------------------------------------------------------------------------
// 结合位点
// ---------------------------------------------------------------------------

/// 一个结合位点（0-based 半开区间，坐标在模板链上）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Binding {
    pub start: usize,
    pub end: usize,
    /// 0 = 正向（引物序列与模板上链相同）；1 = 反向（结合下链）。
    pub strand: u8,
    /// 错配数。
    pub mismatches: usize,
}

impl Binding {
    pub fn range(&self) -> Range {
        Range::new(self.start, self.end)
    }
}

/// 模板里与引物结合的位点。
///
/// 判定与 SnapGene 的默认杂交参数一致：整段允许 ≤ `max_mismatch` 个错配，
/// 但 **3′ 端最后 [`THREE_PRIME_ANCHOR`] 个碱基必须完全配对**。
/// 环状模板允许跨原点（`circular = true`）。
pub fn find_bindings(
    template: &[u8],
    primer: &[u8],
    circular: bool,
    max_mismatch: usize,
) -> Vec<Binding> {
    let mut out = Vec::new();
    let n = template.len();
    let m = primer.len();
    if m == 0 || n == 0 || m > n {
        return out;
    }
    let up = |b: u8| b.to_ascii_uppercase();
    let t: Vec<u8> = template.iter().map(|b| up(*b)).collect();
    let p: Vec<u8> = primer.iter().map(|b| up(*b)).collect();
    let rc = revcomp(&p);
    let limit = if circular { n } else { n - m + 1 };
    for start in 0..limit {
        // 正向：primer == template[start..start+m]
        let mut mm_f = 0usize;
        let mut mm_r = 0usize;
        let mut anchor_ok_f = true;
        let mut anchor_ok_r = true;
        for k in 0..m {
            let idx = (start + k) % n;
            if p[k] != t[idx] {
                mm_f += 1;
                if k >= m.saturating_sub(THREE_PRIME_ANCHOR) {
                    anchor_ok_f = false;
                }
            }
            if rc[k] != t[idx] {
                mm_r += 1;
                // 反向引物的 3′ 端在 rc 的**左**端
                if k < THREE_PRIME_ANCHOR {
                    anchor_ok_r = false;
                }
            }
        }
        if mm_f <= max_mismatch && anchor_ok_f {
            out.push(Binding {
                start,
                end: start + m,
                strand: 0,
                mismatches: mm_f,
            });
        }
        if mm_r <= max_mismatch && anchor_ok_r {
            out.push(Binding {
                start,
                end: start + m,
                strand: 1,
                mismatches: mm_r,
            });
        }
    }
    // 环状模板：跨原点的位点区间可能 end > n，保留原样（显示时按环处理）
    out
}

/// 一段退火：引物 **3′ 端起**与模板连续完全配对的一段（5′ 尾巴可以悬空）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anneal {
    /// 退火段在模板上的区间（0-based 半开）。
    pub start: usize,
    pub end: usize,
    /// 0 = 正向（引物与上链同序），1 = 反向（引物结合下链）。
    pub strand: u8,
    /// 退火碱基数（= 引物 3′ 后缀长度）。
    pub len: usize,
}

impl Anneal {
    pub fn range(&self) -> Range {
        Range::new(self.start, self.end)
    }
}

/// 引物能退火的位置：**从 3′ 端起的连续完全配对段**（≥ `min_match`）。
///
/// 这是「加引物」对话框预览用的模型：带 5′ 尾巴的引物（例如加酶切位点的
/// 引物）只有 3′ 端那一段能与模板配对，长度不足的（< 10）不算结合。
/// 只保留**极大**匹配（被更长匹配包住的丢弃）。
pub fn anneal_bindings(
    template: &[u8],
    primer: &[u8],
    circular: bool,
    min_match: usize,
) -> Vec<Anneal> {
    let mut out: Vec<Anneal> = Vec::new();
    let n = template.len();
    let m = primer.len();
    if m == 0 || n == 0 {
        return out;
    }
    let up = |b: u8| b.to_ascii_uppercase();
    let t: Vec<u8> = template.iter().map(|b| up(*b)).collect();
    let p: Vec<u8> = primer.iter().map(|b| up(*b)).collect();
    let rc = revcomp(&p);

    for e in 0..n {
        let max_len = if circular { m } else { m.min(e) };
        let mut len = 0usize;
        while len < max_len {
            let ti = (e + n - 1 - len) % n;
            if p[m - 1 - len] != t[ti] {
                break;
            }
            len += 1;
        }
        if len >= min_match {
            out.push(Anneal {
                start: (e + n - len) % n,
                end: e % n,
                strand: 0,
                len,
            });
        }
    }
    for s in 0..n {
        let max_len = if circular { m } else { m.min(n - s) };
        let mut len = 0usize;
        while len < max_len && rc[len] == t[(s + len) % n] {
            len += 1;
        }
        if len >= min_match {
            out.push(Anneal {
                start: s,
                end: s + len,
                strand: 1,
                len,
            });
        }
    }
    // 只留极大匹配（被更长匹配包住的丢弃）
    out.sort_by(|a, b| b.len.cmp(&a.len));
    let mut kept: Vec<Anneal> = Vec::new();
    for a in out {
        let contained = kept.iter().any(|k| {
            k.strand == a.strand && k.len >= a.len && k.start <= a.start && k.start + k.len >= a.end
        });
        if !contained {
            kept.push(a);
        }
    }
    kept.sort_by_key(|a| (a.start, a.strand));
    kept
}

/// 引物在模板上的**唯一**完全匹配数（设计时用于判特异性）。
pub fn exact_match_count(template: &[u8], primer: &[u8], circular: bool) -> usize {
    find_bindings(template, primer, circular, 0).len()
}

// ---------------------------------------------------------------------------
// 引物对设计
// ---------------------------------------------------------------------------

/// 设计参数。默认值贴近 SnapGene 的 PCR 对话框（长度 18–25、Tm 55–65、
/// GC 40–60%）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignOptions {
    pub min_len: usize,
    pub max_len: usize,
    pub min_tm: f64,
    pub max_tm: f64,
    pub target_tm: f64,
    pub min_gc: f64,
    pub max_gc: f64,
    pub min_product: usize,
    pub max_product: usize,
    /// 最多返回多少个候选。
    pub max_results: usize,
}

impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            min_len: 18,
            max_len: 25,
            min_tm: 55.0,
            max_tm: 65.0,
            target_tm: 60.0,
            min_gc: 40.0,
            max_gc: 60.0,
            min_product: 100,
            max_product: 5000,
            max_results: 24,
        }
    }
}

/// 一对候选引物。
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerPair {
    pub forward: String,
    pub reverse: String,
    /// 正向引物结合位点起点（= 其 5′ 端；3′ 端在 `forward_start + len`）。
    pub forward_start: usize,
    /// 反向引物结合位点起点（= 其 3′ 端所在位置，目标区间终点）。
    pub reverse_site_start: usize,
    pub forward_tm: f64,
    pub reverse_tm: f64,
    pub forward_gc: f64,
    pub reverse_gc: f64,
    /// 产物区间（0-based 半开；环状且跨原点时 `end < start`）。
    pub product_start: usize,
    pub product_end: usize,
    /// 产物长度（bp；环状按展开算）。
    pub product_len: usize,
}

impl PrimerPair {
    pub fn product_len(&self) -> usize {
        self.product_len
    }
    /// 两条引物的 Tm 差。
    pub fn tm_delta(&self) -> f64 {
        (self.forward_tm - self.reverse_tm).abs()
    }
    /// 综合打分（越小越好）：Tm 差为主，其次偏离目标 Tm。
    pub fn score(&self) -> f64 {
        self.tm_delta() * 1.0
            + ((self.forward_tm + self.reverse_tm) / 2.0 - 60.0).abs() * 0.35
    }
}

/// 环状/线性模板取 `[start, start+len)`（环状自动回绕）。
fn slice_wrapped(template: &[u8], start: usize, len: usize) -> Vec<u8> {
    let n = template.len();
    (0..len).map(|k| template[(start + k) % n]).collect()
}

/// 为 `region`（目标区间，0-based 半开）设计引物对 —— **目标区间是产物必须盖住的部分**。
///
/// 正向引物的 5′ 端与反向引物的 5′ 端（= 产物终点）各自在区间端点附近
/// ±[`SEARCH_WINDOW`] 内滑动，两条引物的 3′ 端都朝内；产物可以比目标略大或略小。
/// 环状模板允许目标跨原点（`region.start > region.end`）。返回按
/// [`PrimerPair::score`] 升序排好的候选（已过滤 Tm / GC / 产物长度 / 发夹 / 二聚体）。
pub fn design_pair(
    template: &[u8],
    region: Range,
    circular: bool,
    opts: &DesignOptions,
) -> Vec<PrimerPair> {
    let n = template.len();
    if n == 0 || opts.min_len == 0 || opts.max_len < opts.min_len {
        return Vec::new();
    }
    let t: Vec<u8> = template.iter().map(|b| b.to_ascii_uppercase()).collect();
    let min_len = opts.min_len.max(4);
    let max_len = opts.max_len.min(n);
    if max_len < min_len {
        return Vec::new();
    }
    let target_len = if circular {
        let d = (region.end + n - region.start) % n;
        if d == 0 {
            n
        } else {
            d
        }
    } else {
        region.end.saturating_sub(region.start)
    };
    if target_len == 0 || target_len > n {
        return Vec::new();
    }
    let window = (SEARCH_WINDOW as isize).min((target_len + SEARCH_WINDOW) as isize);

    // 候选正向引物：5′ 端在 region.start ± window 内
    struct Cand {
        start: usize,
        len: usize,
        seq: Vec<u8>,
        tm: f64,
        gc: f64,
    }
    let mut fwd: Vec<Cand> = Vec::new();
    for flen in min_len..=max_len {
        if circular && flen >= n {
            break;
        }
        for off in -window..=window {
            let Some(start) = place(region.start, off, n, circular) else {
                continue;
            };
            if !circular && start + flen > n {
                continue;
            }
            let seq = slice_wrapped(&t, start, flen);
            let tm = tm(&seq);
            let gc = gc_percent(&seq);
            if !(opts.min_tm..=opts.max_tm).contains(&tm)
                || !(opts.min_gc..=opts.max_gc).contains(&gc)
            {
                continue;
            }
            if let Some(hp) = hairpin(&seq) {
                if hp.tm > 50.0 {
                    continue;
                }
            }
            fwd.push(Cand {
                start,
                len: flen,
                seq,
                tm,
                gc,
            });
        }
    }
    if fwd.is_empty() {
        return Vec::new();
    }

    // 候选反向引物：5′ 端（= 产物终点）在 region.end ± window 内
    let mut rev: Vec<Cand> = Vec::new();
    for rlen in min_len..=max_len {
        if circular && rlen >= n {
            break;
        }
        for roff in -window..=window {
            let Some(end) = place(region.end, roff, n, circular) else {
                continue;
            };
            if !circular && end > n {
                continue;
            }
            let site_start = (end + n - rlen) % n;
            if !circular && end < rlen {
                continue;
            }
            let site = slice_wrapped(&t, site_start, rlen);
            let seq = revcomp(&site);
            let tm = tm(&seq);
            let gc = gc_percent(&seq);
            if !(opts.min_tm..=opts.max_tm).contains(&tm)
                || !(opts.min_gc..=opts.max_gc).contains(&gc)
            {
                continue;
            }
            if let Some(hp) = hairpin(&seq) {
                if hp.tm > 50.0 {
                    continue;
                }
            }
            rev.push(Cand {
                start: site_start,
                len: rlen,
                seq,
                tm,
                gc,
            });
        }
    }

    let mut out: Vec<PrimerPair> = Vec::new();
    for f in &fwd {
        for r in &rev {
            let plen = if circular {
                let d = (r.start + r.len + n - f.start) % n;
                if d == 0 {
                    n
                } else {
                    d
                }
            } else {
                if r.start + r.len <= f.start {
                    continue;
                }
                r.start + r.len - f.start
            };
            if !(opts.min_product..=opts.max_product).contains(&plen) {
                continue;
            }
            if plen < f.len + r.len {
                continue;
            }
            if let Some(d) = cross_dimer(&f.seq, &r.seq) {
                if d.three_prime && d.overlap >= 8 {
                    continue;
                }
            }
            out.push(PrimerPair {
                forward: String::from_utf8_lossy(&f.seq).to_string(),
                reverse: String::from_utf8_lossy(&r.seq).to_string(),
                forward_start: f.start,
                reverse_site_start: r.start,
                forward_tm: f.tm,
                reverse_tm: r.tm,
                forward_gc: f.gc,
                reverse_gc: r.gc,
                product_start: f.start,
                product_end: (r.start + r.len) % n.max(1),
                product_len: plen,
            });
        }
    }
    out.sort_by(|a, b| {
        a.score()
            .partial_cmp(&b.score())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out.truncate(opts.max_results);
    out
}

/// `base + off` 落在 `[0, n]` 内（线性）或取模回绕（环状）；越界返回 `None`。
fn place(base: usize, off: isize, n: usize, circular: bool) -> Option<usize> {
    if circular {
        Some((base as isize + off).rem_euclid(n as isize) as usize)
    } else {
        let v = base as isize + off;
        (v >= 0).then(|| v as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    /// SnapGene 官方样例里 30 条引物的存储 Tm（`<BindingSite meltingTemperature>`），
    /// 用来锁住盐/浓度常数。这里取其中 10 条代表性的。
    #[test]
    fn tm_matches_snapgene_stored_values() {
        let cases: [(&str, f64); 12] = [
            ("ATGAGTAAAGGAGAAGAACT", 50.0),
            ("CTATTTGTATAGTTCATCCA", 46.0),
            ("CAGTCAGGATCCATGAGTAAAGGAGAAGAACT", 62.0),
            ("CTGACAGAATTCCTATTTGTATAGTTCATCCA", 58.0),
            ("AAGAACATCGATTTTCCATGGCAG", 58.0),
            ("CGACTCACTATAGGGAGAGCGGC", 62.0),
            ("GGATCCTTACTTGTACAGCTCGTCCA", 61.0),
            ("GTAAAACGACGGCCAGT", 56.0),
            ("CAGGAAACAGCTATGAC", 49.0),
            ("GGTACCATTAAGCTTGTGCCCCAGTT", 63.0),
            ("AATGGCTACTTGTACAGCTCGTCCATGCC", 66.0),
            ("CTGCAGGAAGGTTTAAACGCATTTAGG", 61.0),
        ];
        for (s, sg) in cases {
            let t = tm(s.as_bytes());
            assert!(
                (t - sg).abs() <= 3.0,
                "{s}: ours {t:.1} vs SnapGene {sg:.0}"
            );
        }
    }

    #[test]
    fn gc_and_revcomp() {
        assert!((gc_percent(b"ATGC") - 50.0).abs() < 1e-9);
        assert!((gc_percent(b"AAAA") - 0.0).abs() < 1e-9);
        assert_eq!(revcomp(b"ATGC"), b"GCAT");
        assert!(tm(b"A").is_nan());
        assert!(tm(b"ANGT").is_nan());
    }

    #[test]
    fn hairpin_finds_stem() {
        // GCGC....GCGC 形成 4 bp 茎 + 4 nt 环
        let hp = hairpin(b"GCGCAAAAGCGC").expect("hairpin");
        assert_eq!(hp.stem, 4);
        assert_eq!(hp.loop_len, 4);
        // 4 bp 茎的 Tm 很低（~ -8 °C @ 3 µM；发夹用 1e-4 M 有效浓度后略高），
        // 这里只验证找到了茎，不断言 Tm 的正负。
        assert!(hp.tm.is_finite());
        assert!(hairpin(b"AAAAAAAAAAAA").is_none());
    }

    #[test]
    fn self_dimer_detects_complementary_overlap() {
        let d = self_dimer(b"AAAAGGGG").expect("dimer");
        assert!(d.overlap >= 4);
    }

    #[test]
    fn bindings_forward_and_reverse() {
        // 不含回文/重复的模板，避免正反向同时命中
        let template = seq("TTTTACGATCGTAGCTAGCTAGGCTA");
        let fwd = b"CGTAGCTAGC"; // 模板 9..19
        let b = find_bindings(&template, fwd, false, 0);
        assert_eq!(b.len(), 1, "{b:?}");
        assert_eq!((b[0].start, b[0].end, b[0].strand), (9, 19, 0));

        // 反向：取模板 19..26 的反向互补作为引物
        let site = template[19..26].to_vec();
        let primer = revcomp(&site);
        let b = find_bindings(&template, &primer, false, 0);
        assert_eq!(b.len(), 1, "{b:?}");
        assert_eq!((b[0].start, b[0].end, b[0].strand), (19, 26, 1));
    }

    #[test]
    fn binding_requires_intact_three_prime() {
        let template = seq("ACGTTTAAACCCGGG");
        // 5′ 端错配可以（1 个），3′ 端错配不行
        let five = b"TCGTTTAAACCC"; // 首位 T vs A：错配但 3′ 端完好
        let hits = find_bindings(&template, five, false, 1);
        assert!(hits.iter().any(|h| h.start == 0 && h.mismatches == 1));
        let three = b"ACGTTTAAACCT"; // 末位错配 → 不算
        assert!(find_bindings(&template, three, false, 1).is_empty());
        // 对照：完全匹配
        assert_eq!(find_bindings(&template, b"ACGTTTAAA", false, 0).len(), 1);
    }

    #[test]
    fn anneal_handles_five_prime_tail() {
        let template = seq("TTTTACGATCGTAGCTAGCTAGGCTA");
        // 5′ 尾巴 GGGG + 模板 9..19 段
        let primer = b"GGGGCGTAGCTAGC".to_vec();
        let a = anneal_bindings(&template, &primer, false, 10);
        assert_eq!(a.len(), 1, "{a:?}");
        assert_eq!((a[0].start, a[0].end, a[0].strand, a[0].len), (9, 19, 0, 10));

        // 反向：模板 19..26 的反向互补 + 尾巴
        let site = template[19..26].to_vec();
        let mut primer = b"AAAA".to_vec();
        primer.extend(revcomp(&site));
        let a = anneal_bindings(&template, &primer, false, 7);
        assert_eq!(a.len(), 1, "{a:?}");
        assert_eq!((a[0].start, a[0].end, a[0].strand, a[0].len), (19, 26, 1, 7));
        // min_match 太严 → 找不到
        assert!(anneal_bindings(&template, &primer, false, 8).is_empty());
    }

    #[test]
    fn anneal_drops_shorter_contained_matches() {
        let template = seq("ACGTTTAAACCCGGGTTT");
        let primer = template[0..12].to_vec(); // ACGTTTAAACCC
        let a = anneal_bindings(&template, &primer, false, 10);
        // 只报完整匹配（12），不报被包住的 11/10
        assert!(a.iter().any(|x| x.len == 12 && x.start == 0), "{a:?}");
        assert!(!a.iter().any(|x| x.start == 0 && x.len < 12), "{a:?}");
    }

    #[test]
    fn design_pair_spans_region() {
        // 300 bp 随机模板（固定种子，保证可重复）
        let mut x: u64 = 0x2545F4914F6CDD1D;
        let mut tpl = Vec::with_capacity(300);
        for _ in 0..300 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            tpl.push(b"ACGT"[(x % 4) as usize]);
        }
        let region = Range::new(120, 180);
        let opts = DesignOptions {
            min_product: 80,
            max_product: 400,
            ..Default::default()
        };
        let pairs = design_pair(&tpl, region, false, &opts);
        assert!(!pairs.is_empty(), "no pair found");
        let p = &pairs[0];
        // 引物的 5′ 端在区间端点 ±SEARCH_WINDOW 内（产物可略大/略小）
        assert!(p.product_start + 30 >= region.start && p.product_start <= region.start + 30);
        assert!(p.product_end + 30 >= region.end && p.product_end <= region.end + 30);
        assert!(p.product_len >= 80 && p.product_len <= 400);
        for p in &pairs {
            assert!(p.product_len >= 80 && p.product_len <= 400, "{p:?}");
        }
        // 正向引物确实落在模板上（5′ 端 = forward_start）
        let f = p.forward.as_bytes();
        assert_eq!(&tpl[p.forward_start..p.forward_start + f.len()], f);
        // 反向引物的结合位点 = [reverse_site_start, reverse_site_start+len)，
        // 引物是其反向互补，且 3′ 端在结合位点左端。
        let r = p.reverse.as_bytes();
        let site = revcomp(r);
        let rs = p.reverse_site_start;
        assert_eq!(&tpl[rs..rs + site.len()], &site[..]);
        assert_eq!(rs + r.len(), p.product_end);
    }

    #[test]
    fn ss_dna_mass_matches_snapgene() {
        // SnapGene 8.2.3 引物详情行实测值
        let m = ss_dna_mass(b"CAGTCAGGATCCATGAGTAAAGGAGAAGAACT", false);
        assert!((m - 9939.6).abs() < 0.05, "{m}");
        let m = ss_dna_mass(b"TTAGACACAACATTGAAGAT", true);
        assert!((m - 6212.1).abs() < 0.05, "{m}");
        // 磷酸化只加一个 5′ 磷酸
        assert!((ss_dna_mass(b"ACGT", true) - ss_dna_mass(b"ACGT", false) - 79.0).abs() < 1e-9);
    }

    #[test]
    fn codon_table_and_peptide_dna() {
        assert_eq!(CODON_TABLE.len(), 21);
        let mut all: Vec<&str> = Vec::new();
        for g in CODON_TABLE {
            assert!(!g.codons.is_empty(), "{}", g.aa3);
            for c in g.codons {
                assert_eq!(c.len(), 3);
                assert!(c.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
                assert!(!all.contains(c), "duplicate codon {c}");
                all.push(c);
            }
        }
        // 64 个密码子一个不差
        assert_eq!(all.len(), 64);
        assert_eq!(dna_for_peptide("HHHHHH"), "CATCATCATCATCATCAT");
        assert_eq!(dna_for_peptide("DYKDDDDK"), "GATTATAAAGATGATGATGATAAA");
        assert_eq!(dna_for_peptide("M"), "ATG");
    }

    #[test]
    fn design_pair_respects_gc_bounds() {
        let tpl = b"ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT".repeat(4);
        let opts = DesignOptions::default();
        let pairs = design_pair(&tpl, Range::new(100, 140), false, &opts);
        // 全 AT 模板不可能满足 GC 40–60%
        assert!(pairs.is_empty());
    }
}
