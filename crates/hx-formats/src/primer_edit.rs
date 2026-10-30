//! 引物（block 5）编辑：对 `<Primers>` XML 做**局部手术**。
//!
//! 与特征编辑（[`crate::feature_edit`]）同一策略：只在必要处改字符串，
//! 其余字节原样保留（`HybridizationParams`、`PrimerAnnotation`、未知属性等）。
//!
//! SnapGene 的引物元素形状（实测，8.2.2）：
//!
//! ```xml
//! <Primer recentID="0" name="GFP_FOR" sequence="CAGTCAGGATCCATG..."
//!         description="&lt;html&gt;…&lt;/html&gt;" dateAdded="2023-05-19T00:04:55Z">
//!   <BindingSite location="0-19" boundStrand="0"
//!                annealedBases="ATGAGTAAAGGAGAAGAACT" meltingTemperature="50">
//!     <Component bases="CAGTCAGGATCC"/>
//!     <Component hybridizedRange="0-19" bases="ATGAGTAAAGGAGAAGAACT"/>
//!   </BindingSite>
//!   <BindingSite simplified="1" location="0-19" boundStrand="0"
//!                annealedBases="ATGAGTAAAGGAGAAGAACT" meltingTemperature="50">…</BindingSite>
//! </Primer>
//! ```
//!
//! `location` 是 **0-based 闭区间**（与 [`crate::annot::parse_primers`] 的
//! `Range { start, end }` 相差 1：`end` 为半开）。

use crate::annot::Range;

/// 一个结合位点的写侧模型。
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerSiteSpec {
    /// 模板上的结合区间（0-based 半开；环状可 `start > end` 表示跨原点）。
    pub range: Range,
    /// 0 = 正向（结合上链），1 = 反向（结合下链）。
    pub strand: u8,
    /// 引物里实际退火的那一段（3′ 端，5′→3′）。
    pub annealed: String,
    /// 退火段 Tm（°C）。
    pub tm: f64,
}

/// 一个引物的写侧模型。
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerSpec {
    pub name: String,
    /// 完整引物序列（5′→3′，可含 5′ 尾巴）。
    pub sequence: String,
    /// 描述（纯文本；空表示不写）。
    pub description: String,
    /// 原样保留的 `dateAdded`（编辑时沿用；空表示写当前 UTC 时间）。
    pub date_added: String,
    /// 5′ 磷酸化（`phosphorylated="1"`）。
    pub phosphorylated: bool,
    pub sites: Vec<PrimerSiteSpec>,
}

impl PrimerSpec {
    pub fn new(name: impl Into<String>, sequence: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            sequence: sequence.into(),
            description: String::new(),
            date_added: String::new(),
            phosphorylated: false,
            sites: Vec::new(),
        }
    }
}

/// 引物编辑（可重放的原子操作）。
#[derive(Debug, Clone, PartialEq)]
pub enum PrimerEdit {
    Recalculate { params: crate::annot::HybridizationParams, sites: Vec<(String,Vec<PrimerSiteSpec>)> },
    Hybridization(crate::annot::HybridizationParams),
    /// 新增一个引物（追加到 `</Primers>` 前）。
    Add(Box<PrimerSpec>),
    /// 按名字删除（匹配 `<Primer … name="…">`）。
    Remove { name: String },
    /// 改名（保留序列/位点/描述）。
    Rename { name: String, new_name: String },
    /// 整条替换（序列/位点/描述；`spec.name` 用于定位旧条目）。
    Replace { name: String, spec: Box<PrimerSpec> },
}

impl PrimerEdit {
    pub fn label(&self) -> String {
        match self {
            PrimerEdit::Recalculate {..} => "recalculate primer binding sites".into(),
            PrimerEdit::Hybridization(_) => "edit hybridization parameters".into(),
            PrimerEdit::Add(p) => format!("add primer {}", p.name),
            PrimerEdit::Remove { name } => format!("remove primer {name}"),
            PrimerEdit::Rename { name, new_name } => format!("rename primer {name} → {new_name}"),
            PrimerEdit::Replace { name, .. } => format!("edit primer {name}"),
        }
    }
}

// ---------------------------------------------------------------------------
// 字符串手术
// ---------------------------------------------------------------------------

/// XML 属性转义（与项目写回契约一致：只用四个实体 + `&#10;`）。
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\n' => out.push_str("&#10;"),
            c => out.push(c),
        }
    }
    out
}

fn attr_in<'a>(tag: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key}=\"");
    let i = tag.find(&pat)? + pat.len();
    let j = i + tag[i..].find('"')?;
    Some(&tag[i..j])
}

/// 找到名为 `name` 的 `<Primer …>…</Primer>` 的字节区间；找不到返回 `None`。
fn primer_span(text: &str, name: &str) -> Option<(usize, usize)> {
    let mut from = 0usize;
    while let Some(rel) = text[from..].find("<Primer ") {
        let start = from + rel;
        let open_end = start + text[start..].find('>')? + 1;
        let open = &text[start..open_end];
        if attr_in(open, "name") == Some(name) {
            let end = match text[open_end..].find("</Primer>") {
                Some(i) => open_end + i + "</Primer>".len(),
                None => open_end,
            };
            return Some((start, end));
        }
        from = open_end;
    }
    None
}

fn next_valid_id(text: &str) -> u32 {
    text.find("<Primers")
        .and_then(|i| {
            let tag = &text[i..i + text[i..].find('>').map(|j| j + 1).unwrap_or(0)];
            attr_in(tag, "nextValidID").and_then(|v| v.parse().ok())
        })
        .unwrap_or(1)
}

fn set_next_valid_id(text: &mut String, id: u32) {
    let Some(i) = text.find("<Primers") else {
        return;
    };
    let Some(j) = text[i..].find('>') else {
        return;
    };
    let tag = &text[i..i + j];
    if let Some(pos) = tag.find("nextValidID=\"") {
        let vstart = i + pos + "nextValidID=\"".len();
        let vend = vstart + text[vstart..].find('"').unwrap_or(0);
        text.replace_range(vstart..vend, &id.to_string());
    } else {
        text.insert_str(i + "<Primers".len(), &format!(" nextValidID=\"{id}\""));
    }
}

/// UTC 时间戳 `YYYY-MM-DDTHH:MM:SSZ`（不引依赖；用 Hinnant 的 civil_from_days）。
fn utc_timestamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3_600,
        (rem % 3_600) / 60,
        rem % 60
    )
}

/// 生成一条 `<Primer>`；`recent_id` 只用于 SnapGene 的“最近使用”排序。
fn build_primer(spec: &PrimerSpec, recent_id: u32) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "<Primer recentID=\"{recent_id}\" name=\"{}\" sequence=\"{}\"",
        escape(&spec.name),
        escape(&spec.sequence)
    ));
    if !spec.description.is_empty() {
        // SnapGene 存 HTML；我们只写纯文本（属性转义一次），读回时原样展示。
        out.push_str(&format!(" description=\"{}\"", escape(&spec.description)));
    }
    let date = if spec.date_added.is_empty() {
        utc_timestamp()
    } else {
        spec.date_added.clone()
    };
    out.push_str(&format!(" dateAdded=\"{}\">", escape(&date)));
    // 5′ 磷酸化：SnapGene 用 `<Primer … phosphorylated="1">`（见二进制里的
    // `…useAdjustedBindingSites.phosphorylated=` 与 `phosphorylated.PrimerAnnotation`
    // 字符串）。位置在 dateAdded 之后。
    if spec.phosphorylated {
        out.insert_str(out.rfind('>').unwrap_or(out.len()), " phosphorylated=\"1\"");
    }

    for site in &spec.sites {
        // 0-based 闭区间；跨原点时 `end` 是回绕后的（半开）终点，
        // 因此两种情形都写 `end - 1`（实测 SnapGene：`location="2824-4"`）。
        let loc = format!("{}-{}", site.range.start, site.range.end.saturating_sub(1));
        let tm = if site.tm.is_finite() {
            format!("{:.0}", site.tm)
        } else {
            "0".to_string()
        };
        // 不配对的那段按 SnapGene 的写法排在退火段的外侧：5′ 尾巴（退火段是
        // `sequence` 后缀）在前，序列末尾多打的碱基（退火段是前缀）在后。
        let (pre, post) = if site.annealed.is_empty() {
            ("", "")
        } else if spec.sequence.ends_with(&site.annealed) {
            (&spec.sequence[..spec.sequence.len() - site.annealed.len()], "")
        } else if spec.sequence.starts_with(&site.annealed) {
            ("", &spec.sequence[site.annealed.len()..])
        } else {
            ("", "")
        };
        let comps = {
            let mut c = String::new();
            if !pre.is_empty() {
                c.push_str(&format!("<Component bases=\"{}\"/>", escape(pre)));
            }
            c.push_str(&format!(
                "<Component hybridizedRange=\"{loc}\" bases=\"{}\"/>",
                escape(&site.annealed)
            ));
            if !post.is_empty() {
                c.push_str(&format!("<Component bases=\"{}\"/>", escape(post)));
            }
            c
        };
        for simplified in [false, true] {
            out.push_str(&format!(
                "<BindingSite{} location=\"{}\" boundStrand=\"{}\" annealedBases=\"{}\" meltingTemperature=\"{}\">",
                if simplified { " simplified=\"1\"" } else { "" },
                loc,
                site.strand,
                escape(&site.annealed),
                tm
            ));
            out.push_str(&comps);
            out.push_str("</BindingSite>");
        }
    }
    out.push_str("</Primer>");
    out
}

/// 施加一批引物编辑；XML 不是合法 UTF-8 时原样返回。
pub fn apply_primer_edits(xml: &[u8], edits: &[PrimerEdit]) -> Vec<u8> {
    let Ok(mut text) = String::from_utf8(xml.to_vec()) else {
        return xml.to_vec();
    };
    for edit in edits {
        match edit {
            PrimerEdit::Recalculate { params, sites } => {
                text=String::from_utf8(apply_primer_edits(text.as_bytes(),&[PrimerEdit::Hybridization(*params)])).unwrap_or(text);
                for(name,sites)in sites {
                    if let Some((a,b))=primer_span(&text,name){
                        let mut node=text[a..b].to_owned();
                        while let Some(start)=node.find("<BindingSite"){
                            let Some(gt)=node[start..].find('>').map(|x|start+x+1)else{break};
                            let end=if node[start..gt].trim_end_matches('>').trim_end().ends_with('/') {gt}else{match node[gt..].find("</BindingSite>"){Some(i)=>gt+i+14,None=>break}};
                            node.replace_range(start..end,"");
                        }
                        let spec=PrimerSpec{sites:sites.clone(),..PrimerSpec::new(name,"")};
                        let rendered=build_primer(&spec,0);
                        if let Some(open)=rendered.find('>'){
                            if let Some(end)=rendered.rfind("</Primer>"){
                                if let Some(at)=node.rfind("</Primer>"){node.insert_str(at,&rendered[open+1..end]);}
                            }
                        }
                        text.replace_range(a..b,&node);
                    }
                }
            }
            PrimerEdit::Hybridization(p) => {
                let tag = format!("<HybridizationParams minContinuousMatchLen=\"{}\" allowMismatch=\"{}\" minMeltingTemperature=\"{}\" showAdditionalFivePrimeMatches=\"{}\" minimumFivePrimeAnnealing=\"{}\"/>", p.min_continuous_match_len, u8::from(p.allow_mismatch), p.min_melting_temperature, u8::from(p.show_additional_five_prime_matches), p.minimum_five_prime_annealing);
                if let Some(start) = text.find("<HybridizationParams") {
                    if let Some(end) = text[start..].find("/>") { text.replace_range(start..start + end + 2, &tag); }
                } else if let Some(start) = text.find("<Primers") {
                    if let Some(at) = text[start..].find('>') { text.insert_str(start + at + 1, &tag); }
                }
            }
            PrimerEdit::Add(spec) => {
                if !text.contains("<Primers") {
                    text = "<?xml version=\"1.0\"?><Primers nextValidID=\"1\"></Primers>".into();
                }
                let id = next_valid_id(&text);
                let node = build_primer(spec, id);
                match text.rfind("</Primers>") {
                    Some(i) => text.insert_str(i, &node),
                    None => text.push_str(&node),
                }
                set_next_valid_id(&mut text, id + 1);
            }
            PrimerEdit::Remove { name } => {
                if let Some((a, b)) = primer_span(&text, name) {
                    text.replace_range(a..b, "");
                }
            }
            PrimerEdit::Rename { name, new_name } => {
                if let Some((a, b)) = primer_span(&text, name) {
                    let open_end = a + text[a..].find('>').unwrap_or(0);
                    let open = text[a..open_end].to_string();
                    if let Some(pos) = open.find("name=\"") {
                        let vs = a + pos + "name=\"".len();
                        let ve = vs + text[vs..].find('"').unwrap_or(0);
                        text.replace_range(vs..ve, &escape(new_name));
                    }
                    let _ = b;
                }
            }
            PrimerEdit::Replace { name, spec } => {
                if let Some((a, b)) = primer_span(&text, name) {
                    let old = text[a..b].to_string();
                    let open_end = a + text[a..].find('>').unwrap_or(0) + 1;
                    // 保留原 recentID（若存在）
                    let open = &text[a..open_end];
                    let rid = attr_in(open, "recentID")
                        .and_then(|v| v.parse::<u32>().ok())
                        .unwrap_or(next_valid_id(&text));
                    let mut node = build_primer(spec, rid);
                    // 序列没变时保留原有的 `<PrimerAnnotation …/>`（5′ 尾巴上的
                    // 酶位点 / 标签）；序列一变这些位置就不可信了，直接丢弃。
                    if attr_in(&old, "sequence").unwrap_or("") == spec.sequence {
                        let body_start = old.find('>').map(|i| i + 1).unwrap_or(0);
                        let body_end = old.rfind("</Primer>").unwrap_or(old.len());
                        if body_end > body_start {
                            let body = &old[body_start..body_end];
                            let mut keep = String::new();
                            let mut at = 0usize;
                            while let Some(rel) = body[at..].find("<PrimerAnnotation") {
                                let s0 = at + rel;
                                let e0 = body[s0..]
                                    .find("/>")
                                    .map(|i| s0 + i + 2)
                                    .or_else(|| {
                                        body[s0..].find("</PrimerAnnotation>").map(|i| {
                                            s0 + i + "</PrimerAnnotation>".len()
                                        })
                                    })
                                    .unwrap_or(body.len());
                                keep.push_str(&body[s0..e0]);
                                at = e0;
                            }
                            if !keep.is_empty() {
                                if let Some(i) = node.rfind("</Primer>") {
                                    node.insert_str(i, &keep);
                                }
                            }
                        }
                    }
                    // 名称以 spec 为准（允许改名）
                    text.replace_range(a..b, &node);
                }
            }
        }
    }
    text.into_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annot::parse_primers;

    fn empty() -> Vec<u8> {
        b"<?xml version=\"1.0\"?><Primers nextValidID=\"1\"><HybridizationParams minContinuousMatchLen=\"10\"/></Primers>".to_vec()
    }

    fn spec() -> PrimerSpec {
        PrimerSpec {
            name: "F1".into(),
            sequence: "CAGTCAGGATCCATGAGTAAAGGAGAAGAACT".into(),
            description: String::new(),
            date_added: String::new(),
            phosphorylated: false,
            sites: vec![PrimerSiteSpec {
                range: Range::new(0, 20),
                strand: 0,
                annealed: "ATGAGTAAAGGAGAAGAACT".into(),
                tm: 50.3,
            }],
        }
    }

    #[test]
    fn add_round_trips_through_parser() {
        let xml = apply_primer_edits(&empty(), &[PrimerEdit::Add(Box::new(spec()))]);
        let primers = parse_primers(&xml);
        assert_eq!(primers.len(), 1);
        assert_eq!(primers[0].name, "F1");
        assert_eq!(primers[0].sequence, "CAGTCAGGATCCATGAGTAAAGGAGAAGAACT");
        assert_eq!(primers[0].sites.len(), 1);
        assert_eq!(primers[0].sites[0].range, Range::new(0, 20));
        assert_eq!(primers[0].sites[0].strand, 0);
        assert_eq!(primers[0].sites[0].annealed_bases, 20);
        assert_eq!(primers[0].sites[0].melting_temperature, Some(50.0));
        let text = String::from_utf8(xml).unwrap();
        assert!(text.contains("nextValidID=\"2\""), "{text}");
        assert!(text.contains("<HybridizationParams"), "{text}");
        // 简化版 + 详细版各一份
        assert_eq!(text.matches("<BindingSite").count(), 2);
    }

    #[test]
    fn add_then_remove_and_rename() {
        let mut xml = apply_primer_edits(&empty(), &[PrimerEdit::Add(Box::new(spec()))]);
        xml = apply_primer_edits(&xml, &[PrimerEdit::Rename { name: "F1".into(), new_name: "GFP_F".into() }]);
        let p = parse_primers(&xml);
        assert_eq!(p[0].name, "GFP_F");
        assert_eq!(p[0].sites.len(), 1, "rename must keep body");
        xml = apply_primer_edits(&xml, &[PrimerEdit::Remove { name: "GFP_F".into() }]);
        assert!(parse_primers(&xml).is_empty());
    }

    #[test]
    fn three_prime_overhang_goes_after_hybridized_component() {
        // pET28b-Rep1.dna 的 Primer 2：反向引物，退火段是 `sequence` 前缀，
        // 末尾多出来的 CCC 是**3′ 端**的额外碱基。
        let mut s = spec();
        s.name = "Primer 2".into();
        s.sequence = "CTAAAGGGAGCCCCCGATTTAGAGCTTGACGCCC".into();
        s.sites = vec![PrimerSiteSpec {
            range: Range::new(148, 179),
            strand: 1,
            annealed: "CTAAAGGGAGCCCCCGATTTAGAGCTTGACG".into(),
            tm: 66.0,
        }];
        let xml = apply_primer_edits(&empty(), &[PrimerEdit::Add(Box::new(s))]);
        let text = String::from_utf8(xml.clone()).unwrap();
        assert!(
            text.contains(
                "<Component hybridizedRange=\"148-178\" \
                 bases=\"CTAAAGGGAGCCCCCGATTTAGAGCTTGACG\"/><Component bases=\"CCC\"/>"
            ),
            "{text}"
        );
        let p = parse_primers(&xml);
        assert_eq!(p[0].sites[0].annealed_bases, 31);
        assert_eq!(p[0].sequence, "CTAAAGGGAGCCCCCGATTTAGAGCTTGACGCCC");
    }

    #[test]
    fn remove_leaves_other_primers() {
        let mut a = spec();
        a.name = "A".into();
        let mut b = spec();
        b.name = "B".into();
        b.sites[0].strand = 1;
        let xml = apply_primer_edits(
            &empty(),
            &[PrimerEdit::Add(Box::new(a)), PrimerEdit::Add(Box::new(b))],
        );
        assert_eq!(parse_primers(&xml).len(), 2);
        let xml = apply_primer_edits(&xml, &[PrimerEdit::Remove { name: "A".into() }]);
        let p = parse_primers(&xml);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].name, "B");
    }

    #[test]
    fn escaping_in_name_and_description() {
        let mut s = spec();
        s.name = "a<b&\"c".into();
        s.description = "x < y & z".into();
        let xml = apply_primer_edits(&empty(), &[PrimerEdit::Add(Box::new(s))]);
        let text = String::from_utf8(xml.clone()).unwrap();
        assert!(text.contains("name=\"a&lt;b&amp;&quot;c\""), "{text}");
        assert!(text.contains("description=\"x &lt; y &amp; z\""), "{text}");
        assert_eq!(parse_primers(&xml)[0].name, "a<b&\"c");
    }

    #[test]
    fn replace_keeps_one_entry_and_updates_sequence() {
        let xml = apply_primer_edits(&empty(), &[PrimerEdit::Add(Box::new(spec()))]);
        let mut s = spec();
        s.sequence = "TTTTTTTTTTTTTTTTTTTT".into();
        s.sites[0].annealed = s.sequence.clone();
        s.sites[0].tm = 42.0;
        let xml = apply_primer_edits(&xml, &[PrimerEdit::Replace { name: "F1".into(), spec: Box::new(s) }]);
        let p = parse_primers(&xml);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].sequence, "TTTTTTTTTTTTTTTTTTTT");
        assert_eq!(p[0].sites[0].annealed_bases, 20);
        assert_eq!(p[0].sites[0].melting_temperature, Some(42.0));
    }
}
