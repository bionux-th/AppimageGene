//! 特征（block 10）与引物（block 5）XML 的**只读**解析。
//!
//! 这两块的 XML 由 SnapGene 单行写出、结构规整，这里用极简扫描器解析，
//! 不引入 XML 依赖；未识别的部分一律忽略（只读展示，不参与写回）。
//!
//! 坐标约定（与 spec §7.1/§7.4 一致）：
//! - 特征 `<Segment range="A-B">` 为 **1-based 闭区间** ⇒ 0-based 半开 `[A-1, B)`；
//! - 引物 `<BindingSite location="A-B">` 为 **0-based 闭区间** ⇒ `[A, B+1)`；
//! - `A > B` 表示环状分子上跨原点的区间（原样保留，不拆开）。

/// 0-based 半开区间；`start > end` 表示跨原点。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Range {
    pub start: usize,
    pub end: usize,
}

impl Range {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub fn len(self) -> usize {
        if self.start <= self.end {
            self.end - self.start
        } else {
            // 跨原点：调用方需提供序列长度才能算长度；这里只报“非负差值”
            self.end
        }
    }
}

/// 特征的一个片段（多段特征有多个 `<Segment>`）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureSegment {
    pub range: Range,
    /// `#rrggbb`。
    pub color: String,
    /// `type="standard|translated|..."` 原文。
    pub seg_type: String,
    pub translated: bool,
}

/// 一个特征的展示信息。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureAnnotation {
    pub name: String,
    /// `CDS`/`promoter`/`misc_feature`…
    pub kind: String,
    /// 元素上的 `subType`（SnapGene 只在 `Site` 等类型上写，如 `active`/`acetylation`）。
    ///
    /// oracle：SG 的特征工具卡标题旁芯片写 `Site · active` = `type` + `subType`。
    pub sub_type: String,
    /// `/note` 限定符的纯文本（`<br>` 已换成换行；工具卡的 `NOTE` 段）。
    pub note: String,
    /// 文件里的 `visible="0"`（SnapGene 的「隐藏特征」）：图谱与序列视图不画它，
    /// Features 面板里默认不勾选（oracle：THYX_MYCTU.prot 的 source/gene/Protein）。
    pub visible: bool,
    /// 0=无向、1=正向、2=反向、3=双向。
    pub directionality: i32,
    /// 1/2/3 或 -1/-2/-3；绝对值表示首个完整密码子的偏移 + 1。
    pub reading_frame: i32,
    pub translate_first_codon_as_met: bool,
    /// 正向翻译延续的最远终点（文件的 maxRunOn，1-based 闭终点）。
    pub max_run_on: Option<usize>,
    pub segments: Vec<FeatureSegment>,
}

impl FeatureAnnotation {
    /// 编码链方向；不能用它判断是否绘制箭头。
    pub fn is_forward(&self) -> bool {
        self.directionality != 2
    }

    pub fn has_forward_arrow(&self) -> bool {
        matches!(self.directionality, 1 | 3)
    }

    pub fn has_reverse_arrow(&self) -> bool {
        matches!(self.directionality, 2 | 3)
    }
}

/// 引物结合位点。
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerSite {
    pub range: Range,
    /// 0=正/上链，1=反/下链。
    pub strand: u8,
    /// 退火碱基数（SnapGene 写的是序列字符串，这里取长度；缺省为 0）。
    pub annealed_bases: usize,
    /// 文件里存的 Tm（`meltingTemperature`，整数 °C）；缺省为 `None`。
    pub melting_temperature: Option<f64>,
}

/// 一个引物的展示信息。
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerAnnotation {
    pub name: String,
    /// `<Primer sequence="...">` 的序列（可能为空）。
    pub sequence: String,
    /// `description` 原文（SnapGene 存 HTML，可能为空）。
    pub description: String,
    /// `<Primer dateAdded="…">` 原文（UTC，`YYYY-MM-DDTHH:MM:SSZ`；可能为空）。
    pub date_added: String,
    /// `phosphorylated="1"` = 5′ 磷酸化修饰（Primers 面板的 `MODIFICATION` 行）。
    pub phosphorylated: bool,
    /// `boundStrand=0` 的结合位点（简化版优先；无简化版时取详细版首段）。
    pub sites: Vec<PrimerSite>,
}

/// Parameters which govern SnapGene's additional-primer-binding display.
/// Missing attributes use the application defaults (not zero); these values
/// are read-only and the original XML remains byte-preserved on write.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HybridizationParams {
    pub min_continuous_match_len: usize,
    pub allow_mismatch: bool,
    pub min_melting_temperature: f64,
    pub show_additional_five_prime_matches: bool,
    pub minimum_five_prime_annealing: usize,
}

impl Default for HybridizationParams {
    fn default() -> Self {
        Self {
            min_continuous_match_len: 10,
            allow_mismatch: true,
            min_melting_temperature: 40.0,
            show_additional_five_prime_matches: true,
            minimum_five_prime_annealing: 15,
        }
    }
}

// ---------------------------------------------------------------------------
// 极简 XML 扫描
// ---------------------------------------------------------------------------

fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    s.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
}

/// 从 `text[from..]` 找下一个 `<name`（含边界检查）标签，返回 `(属性表, 标签结束位置)`。
pub(crate) fn next_tag<'a>(text: &'a str, name: &str, from: usize) -> Option<(Vec<(&'a str, String)>, usize)> {
    let pat = format!("<{name}");
    let mut at = from;
    loop {
        let i = text.get(at..)?.find(&pat)? + at;
        let after = i + pat.len();
        // 保证是完整标签名（`<Feature` 不匹配 `<Features`）
        match text.as_bytes().get(after) {
            Some(b' ' | b'>' | b'\t' | b'\n' | b'/') => {
                let mut attrs = Vec::new();
                let mut j = after;
                let bytes = text.as_bytes();
                while j < bytes.len() && bytes[j] != b'>' {
                    // 跳到属性名
                    while j < bytes.len() && (bytes[j] as char).is_whitespace() {
                        j += 1;
                    }
                    let name_start = j;
                    while j < bytes.len()
                        && bytes[j] != b'='
                        && bytes[j] != b'>'
                        && !(bytes[j] as char).is_whitespace()
                    {
                        j += 1;
                    }
                    if j >= bytes.len() || bytes[j] != b'=' {
                        // 无值属性/自闭合：结束
                        if j < bytes.len() && bytes[j] == b'>' {
                            break;
                        }
                        if j >= bytes.len() {
                            break;
                        }
                        j += 1;
                        continue;
                    }
                    let key = &text[name_start..j];
                    j += 1;
                    if bytes.get(j) != Some(&b'"') {
                        // 非引号值：跳过
                        while j < bytes.len() && bytes[j] != b'"' && bytes[j] != b'>' {
                            j += 1;
                        }
                        attrs.push((key, String::new()));
                        continue;
                    }
                    j += 1;
                    let v_start = j;
                    while j < bytes.len() && bytes[j] != b'"' {
                        j += 1;
                    }
                    attrs.push((key, unescape(&text[v_start..j])));
                    j += 1;
                }
                return Some((attrs, j));
            }
            _ => {
                at = after;
            }
        }
    }
}

fn attr<'a>(attrs: &'a [(&str, String)], key: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(k, _)| *k == key)
        .map(|(_, v)| v.as_str())
}

fn parse_usize(v: &str) -> Option<usize> {
    v.trim().parse().ok()
}

/// 特征片段区间：`"A-B"`（1-based 闭）→ 0-based 半开。
pub fn parse_feature_range(s: &str) -> Option<Range> {
    let (a, b) = s.split_once('-')?;
    let a = parse_usize(a)?;
    let b = parse_usize(b)?;
    if a == 0 {
        return None;
    }
    Some(Range::new(a - 1, b))
}

/// 引物结合位点区间：`"A-B"`（0-based 闭）→ 0-based 半开。
/// `annealedBases` 既可能是数字（旧式/测试文件），也可能是序列字符串（SnapGene 现状）。
fn parse_annealed(v: &str) -> usize {
    v.parse().unwrap_or_else(|_| v.len())
}

fn parse_tm(v: &str) -> Option<f64> {
    v.parse().ok()
}

pub fn parse_primer_range(s: &str) -> Option<Range> {
    let (a, b) = s.split_once('-')?;
    let a = parse_usize(a)?;
    let b = parse_usize(b)?;
    Some(Range::new(a, b + 1))
}

/// 解析 block 10（特征 XML）。无法识别的元素被跳过。
pub fn parse_features(xml: &[u8]) -> Vec<FeatureAnnotation> {
    let Ok(text) = std::str::from_utf8(xml) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut at = 0usize;
    while let Some((attrs, tag_end)) = next_tag(text, "Feature", at) {
        let name = attr(&attrs, "name").unwrap_or("").to_string();
        let kind = attr(&attrs, "type").unwrap_or("misc_feature").to_string();
        let sub_type = attr(&attrs, "subType").unwrap_or("").to_string();
        // `visible="0"` 是 SnapGene 的「特征被隐藏」（导入的 UniProt 文件里 source/gene
        // /Protein 就是这种）；缺失 = 可见。
        let visible = attr(&attrs, "visible") != Some("0");
        let directionality: i32 = attr(&attrs, "directionality")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        let reading_frame = attr(&attrs, "readingFrame")
            .and_then(|v| v.parse::<i32>().ok())
            .filter(|v| matches!(*v, -3..=-1 | 1..=3))
            .unwrap_or(1);
        let translate_first_codon_as_met = attr(&attrs, "translateFirstCodonAsMet") == Some("1");
        let max_run_on = attr(&attrs, "maxRunOn").and_then(parse_usize);
        // 元素体：到 `</Feature>` 为止（若无则到文本结尾）
        let close = text[tag_end..]
            .find("</Feature>")
            .map(|i| tag_end + i)
            .unwrap_or(text.len());
        let body = &text[tag_end..close];
        let mut segments = Vec::new();
        let mut seg_at = 0usize;
        while let Some((sattrs, s_end)) = next_tag(body, "Segment", seg_at) {
            if let Some(range) = attr(&sattrs, "range").and_then(parse_feature_range) {
                segments.push(FeatureSegment {
                    range,
                    color: attr(&sattrs, "color").unwrap_or("#c0c0c0").to_string(),
                    seg_type: attr(&sattrs, "type").unwrap_or("standard").to_string(),
                    translated: attr(&sattrs, "translated")
                        .map(|v| v == "1")
                        .unwrap_or(false),
                });
            }
            seg_at = s_end;
        }
        out.push(FeatureAnnotation {
            name,
            kind,
            sub_type,
            note: crate::feature_edit::note_of(body),
            visible,
            directionality,
            reading_frame,
            translate_first_codon_as_met,
            max_run_on,
            segments,
        });
        at = close.max(tag_end);
    }
    out
}

/// 解析 block 5（引物 XML）。优先取 `simplified="1"` 的简化版结合位点；
/// 没有简化版时退回详细版（`<Component hybridizedRange>` 的 0-based 闭区间）。
pub fn parse_primers(xml: &[u8]) -> Vec<PrimerAnnotation> {
    let Ok(text) = std::str::from_utf8(xml) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut at = 0usize;
    while let Some((attrs, tag_end)) = next_tag(text, "Primer", at) {
        let name = attr(&attrs, "name").unwrap_or("").to_string();
        let sequence = attr(&attrs, "sequence").unwrap_or("").to_string();
        let close = text[tag_end..]
            .find("</Primer>")
            .map(|i| tag_end + i)
            .unwrap_or(text.len());
        let body = &text[tag_end..close];

        let mut simple: Vec<PrimerSite> = Vec::new();
        let mut detailed: Vec<PrimerSite> = Vec::new();
        let mut b_at = 0usize;
        while let Some((battrs, b_end)) = next_tag(body, "BindingSite", b_at) {
            if attr(&battrs, "simplified")
                .map(|v| v == "1")
                .unwrap_or(false)
            {
                if let Some(range) = attr(&battrs, "location").and_then(parse_primer_range) {
                    simple.push(PrimerSite {
                        range,
                        strand: attr(&battrs, "boundStrand")
                            .and_then(parse_usize)
                            .unwrap_or(0) as u8,
                        annealed_bases: attr(&battrs, "annealedBases")
                            .map(parse_annealed)
                            .unwrap_or(0),
                        melting_temperature: attr(&battrs, "meltingTemperature").and_then(parse_tm),
                    });
                }
            } else if let Some(range) = attr(&battrs, "location").and_then(parse_primer_range) {
                detailed.push(PrimerSite {
                    range,
                    strand: attr(&battrs, "boundStrand")
                        .and_then(parse_usize)
                        .unwrap_or(0) as u8,
                    annealed_bases: attr(&battrs, "annealedBases")
                        .map(parse_annealed)
                        .unwrap_or(0),
                    melting_temperature: attr(&battrs, "meltingTemperature").and_then(parse_tm),
                });
            }
            b_at = b_end;
        }
        if simple.is_empty() {
            // 详细版：<Component hybridizedRange="A-B" .../>（0-based 闭）
            let mut b_at = 0usize;
            while let Some((cattrs, c_end)) = next_tag(body, "Component", b_at) {
                if let Some(range) = attr(&cattrs, "hybridizedRange").and_then(parse_primer_range) {
                    detailed.push(PrimerSite {
                        range,
                        strand: 0,
                        annealed_bases: attr(&cattrs, "bases").map(|b| b.len()).unwrap_or(0),
                        melting_temperature: None,
                    });
                }
                b_at = c_end;
            }
        }
        let sites = if simple.is_empty() { detailed } else { simple };
        out.push(PrimerAnnotation {
            name,
            sequence,
            description: attr(&attrs, "description").unwrap_or("").to_string(),
            date_added: attr(&attrs, "dateAdded").unwrap_or("").to_string(),
            phosphorylated: attr(&attrs, "phosphorylated")
                .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
                .unwrap_or(false),
            sites,
        });
        at = close.max(tag_end);
    }
    out
}

/// Parse the optional `<HybridizationParams>` tag with SnapGene defaults.
pub fn parse_hybridization_params(xml: &[u8]) -> HybridizationParams {
    let Ok(text) = std::str::from_utf8(xml) else { return HybridizationParams::default() };
    let Some((attrs, _)) = next_tag(text, "HybridizationParams", 0) else {
        return HybridizationParams::default();
    };
    let defaults = HybridizationParams::default();
    HybridizationParams {
        min_continuous_match_len: attr(&attrs, "minContinuousMatchLen")
            .and_then(parse_usize).unwrap_or(defaults.min_continuous_match_len),
        allow_mismatch: attr(&attrs, "allowMismatch")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(defaults.allow_mismatch),
        min_melting_temperature: attr(&attrs, "minMeltingTemperature")
            .and_then(parse_tm).unwrap_or(defaults.min_melting_temperature),
        show_additional_five_prime_matches: attr(&attrs, "showAdditionalFivePrimeMatches")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(defaults.show_additional_five_prime_matches),
        minimum_five_prime_annealing: attr(&attrs, "minimumFivePrimeAnnealing")
            .and_then(parse_usize).unwrap_or(defaults.minimum_five_prime_annealing),
    }
}

/// 从 block 13（酶显示偏好）载荷里取显示集合名（off 18 起 ASCII、NUL 填充）。
pub fn enzyme_set_name(payload: &[u8]) -> Option<String> {
    if payload.len() < 18 {
        return None;
    }
    let rest = &payload[18..];
    let end = rest.iter().position(|b| *b == 0).unwrap_or(rest.len());
    let s = std::str::from_utf8(&rest[..end]).ok()?.trim();
    (!s.is_empty()).then(|| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const FEATS: &str = r##"<?xml version="1.0"?><Features nextValidID="2"><Feature recentID="0" name="GFP" directionality="1" type="CDS" allowSegmentOverlaps="0"><Segment range="1-717" color="#05fd14" type="standard" translated="1"/><Q name="product"><V text="&lt;html&gt;x&lt;/html&gt;"/></Q></Feature><Feature recentID="1" name="ori" directionality="2" type="rep_origin"><Segment range="700-100" color="#ff0000" type="standard" translated="0"/></Feature></Features>"##;

    #[test]
    fn features_basic() {
        let f = parse_features(FEATS.as_bytes());
        assert_eq!(f.len(), 2);
        assert_eq!(f[0].name, "GFP");
        assert_eq!(f[0].kind, "CDS");
        assert_eq!(f[0].directionality, 1);
        assert!(f[0].is_forward());
        assert_eq!(f[0].segments.len(), 1);
        assert_eq!(f[0].segments[0].range, Range::new(0, 717));
        assert_eq!(f[0].segments[0].color, "#05fd14");
        assert!(f[0].segments[0].translated);
        // 跨原点 + 反向
        assert_eq!(f[1].segments[0].range, Range::new(699, 100));
        assert!(!f[1].is_forward());
    }

    #[test]
    fn features_multi_segment_and_entities() {
        let xml = r##"<Features><Feature name="A &amp; B" type="misc_feature" directionality="3"><Segment range="10-20" color="#00ff00" type="standard" translated="0"/><Segment range="30-40" color="#0000ff" type="standard" translated="0"/></Feature></Features>"##;
        let f = parse_features(xml.as_bytes());
        assert_eq!(f[0].name, "A & B");
        assert_eq!(f[0].segments.len(), 2);
        assert_eq!(f[0].segments[1].range, Range::new(29, 40));
        assert!(f[0].is_forward()); // 双向按正向
    }

    #[test]
    fn no_features_tag() {
        assert!(parse_features(b"<Primers></Primers>").is_empty());
        assert!(parse_features(b"not xml").is_empty());
    }

    #[test]
    fn omitted_directionality_is_undirected_and_translation_settings_are_preserved() {
        let xml = br#"<Features><Feature name="RBS" type="RBS"><Segment range="1-6"/></Feature><Feature name="lac operator"><Segment range="7-12"/></Feature><Feature name="lacI" directionality="2" readingFrame="-2" translateFirstCodonAsMet="1"><Segment range="13-21" translated="1"/></Feature><Feature name="both" directionality="3"><Segment range="22-30"/></Feature></Features>"#;
        let f = parse_features(xml);
        for feature in &f[..2] {
            assert_eq!(feature.directionality, 0);
            assert!(!feature.has_forward_arrow());
            assert!(!feature.has_reverse_arrow());
            assert!(!feature.translate_first_codon_as_met);
        }
        assert_eq!(f[2].reading_frame, -2);
        assert!(f[2].translate_first_codon_as_met);
        assert!(!f[2].has_forward_arrow());
        assert!(f[2].has_reverse_arrow());
        assert!(f[3].has_forward_arrow() && f[3].has_reverse_arrow());
    }

    #[test]
    fn primers_simplified_preferred() {
        let xml = r#"<?xml version="1.0"?><Primers nextValidID="1"><HybridizationParams minContinuousMatchLen="10"/><Primer name="F1" sequence="ATGC"><BindingSite simplified="0" location="100-119" boundStrand="0" annealedBases="20" meltingTemperature="60"/><BindingSite simplified="1" location="100-119" boundStrand="0" annealedBases="20"/></Primer><Primer name="R1" sequence="TTTT"><BindingSite simplified="1" location="200-219" boundStrand="1" annealedBases="20"/></Primer></Primers>"#;
        let p = parse_primers(xml.as_bytes());
        assert_eq!(p.len(), 2);
        assert_eq!(p[0].name, "F1");
        assert_eq!(p[0].sites.len(), 1, "同时存在时应只取简化版");
        assert_eq!(p[0].sites[0].range, Range::new(100, 120));
        assert_eq!(p[0].sites[0].strand, 0);
        assert_eq!(p[1].sites[0].strand, 1);
        assert_eq!(p[1].sequence, "TTTT");
    }

    #[test]
    fn hybridization_parameters_parse_saved_snapgene_values_and_defaults() {
        let p = parse_hybridization_params(br#"<Primers><HybridizationParams minContinuousMatchLen="10" allowMismatch="1" minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" minimumFivePrimeAnnealing="15"/></Primers>"#);
        assert_eq!(p, HybridizationParams::default());
        let p = parse_hybridization_params(br#"<Primers><HybridizationParams minContinuousMatchLen="8" allowMismatch="0" minMeltingTemperature="35" showAdditionalFivePrimeMatches="0" minimumFivePrimeAnnealing="12"/></Primers>"#);
        assert_eq!(p.min_continuous_match_len, 8);
        assert!(!p.allow_mismatch);
        assert_eq!(p.min_melting_temperature, 35.0);
        assert!(!p.show_additional_five_prime_matches);
        assert_eq!(p.minimum_five_prime_annealing, 12);
        assert_eq!(parse_hybridization_params(b"<Primers/>"), HybridizationParams::default());
    }

    #[test]
    fn primers_detailed_fallback() {
        // 只有详细版：取 <Component hybridizedRange>
        let xml = r#"<Primers><Primer name="P" sequence="GGG"><BindingSite location="50-69" boundStrand="0" annealedBases="20"><Component hybridizedRange="50-59" bases="GGGGGGGGGG"/><Component hybridizedRange="60-69" bases="GGGGGGGGGG"/></BindingSite><BindingSite simplified="1" location="50-69" boundStrand="0" annealedBases="20"/></Primer></Primers>"#;
        let p = parse_primers(xml.as_bytes());
        assert_eq!(p[0].sites.len(), 1, "有简化版时走简化版");
        assert_eq!(p[0].sites[0].range, Range::new(50, 70));

        let only_detailed = r#"<Primers><Primer name="P" sequence="GGG"><Component hybridizedRange="50-59" bases="GGGGGGGGGG"/></Primer></Primers>"#;
        let p2 = parse_primers(only_detailed.as_bytes());
        assert_eq!(p2[0].sites.len(), 1);
        assert_eq!(p2[0].sites[0].range, Range::new(50, 60));
    }

    #[test]
    fn range_parsers() {
        assert_eq!(parse_feature_range("1-1"), Some(Range::new(0, 1)));
        assert_eq!(parse_feature_range("717-717"), Some(Range::new(716, 717)));
        assert_eq!(parse_feature_range("0-5"), None);
        assert_eq!(parse_feature_range("x"), None);
        assert_eq!(parse_primer_range("0-0"), Some(Range::new(0, 1)));
        assert_eq!(parse_primer_range("10-19"), Some(Range::new(10, 20)));
    }

    #[test]
    fn enzyme_set_name_from_block13() {
        let mut p = vec![0u8; 18];
        p[17] = 3;
        p.extend_from_slice(b"Unique 6+ Cutters\0\0\0");
        assert_eq!(enzyme_set_name(&p).as_deref(), Some("Unique 6+ Cutters"));
        assert_eq!(enzyme_set_name(&[0u8; 10]), None);
    }
}
