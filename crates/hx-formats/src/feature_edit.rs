//! 特征编辑（`block 10` 特征 XML 的**局部手术**）。
//!
//! 与序列编辑同一条原则：**只改被编辑的那个 `<Feature>` 元素**，其余字节原样保留
//! （未编辑文件读→写逐字节一致）。因此这里不做「整体重新序列化」，
//! 而是在文本上做定点的元素替换/删除/属性改写。
//!
//! 编辑以 `recentID` 定位元素（SnapGene 自己的身份标识；全部样例的 `block 10`
//! 都有 `recentID`，见 `docs/new-file-and-features-2026-10-18.zh-CN.md`）。
//!
//! 属性词汇表（203 样例实测出现）：
//! `recentID name type directionality allowSegmentOverlaps consecutiveTranslationNumbering
//!   readingFrame swappedSegmentNumbering translationMW hitsStopCodon subType cleavageArrows
//!   visible prioritize maxRunOn maxFusedRunOn consecutiveNumberingStartsFrom detectionMode
//!   isFavorite geneticCode translateFirstCodonAsMet`
//! `<Segment>`：`range color type translated name translationNumberingStartsFrom`。
//! 本模块只改写**已知字段**，未知属性一律原样保留。

/// SnapGene「类型」下拉的 SO 特征类型表（顺序、内容取自真 SnapGene 8.2.2 截图；
/// 见 `docs/new-file-and-features-2026-10-18.zh-CN.md` 图 5）。
pub const FEATURE_TYPES: &[&str] = &[
    "-10_signal",
    "-35_signal",
    "3'UTR",
    "5'UTR",
    "assembly_gap",
    "attenuator",
    "C_region",
    "CAAT_signal",
    "CDS",
    "centromere",
    "conflict",
    "D-loop",
    "D_segment",
    "enhancer",
    "exon",
    "gap",
    "GC_signal",
    "gene",
    "iDNA",
    "intron",
    "J_segment",
    "LTR",
    "mat_peptide",
    "misc_binding",
    "misc_difference",
    "misc_feature",
    "misc_recomb",
    "misc_RNA",
    "misc_signal",
    "misc_structure",
    "mobile_element",
    "modified_base",
    "mRNA",
    "N_region",
    "ncRNA",
    "old_sequence",
    "operon",
    "oriT",
    "polyA_signal",
    "polyA_site",
    "precursor_RNA",
    "prim_transcript",
    "primer_bind",
    "promoter",
    "propeptide",
    "protein",
    "protein_bind",
    "RBS",
    "regulatory",
    "rep_origin",
    "repeat_region",
    "rRNA",
    "S_region",
    "sig_peptide",
    "source",
    "stem_loop",
    "STS",
    "TATA_signal",
    "telomere",
    "terminator",
    "tmRNA",
    "transit_peptide",
    "tRNA",
    "unsure",
    "V_region",
    "V_segment",
    "variation",
];

/// 「无颜色」的取值（SnapGene 的 `color="noColor"`）。
pub const NO_COLOR: &str = "noColor";

/// 颜色选择器调色板（SnapGene 8 的颜色弹出网格；值为样例实测出现过的颜色）。
pub const PALETTE: &[&str] = &[
    "#ff0000", "#ff6666", "#ff9900", "#ffcc00", "#ffff00", "#e1ff2f", //
    "#05fd14", "#00ff00", "#00cc00", "#00ffff", "#00ccff", "#99ccff", //
    "#3366ff", "#0000ff", "#a020f0", "#ff00ff", "#ff7f50", "#cc99b2", //
    "#993366", "#a6acb3", "#ffffff", "#cccccc", "#808080", "#000000", //
];

/// 每种特征类型的默认颜色（样例统计出的 SnapGene 内建默认值）。
pub fn default_color(kind: &str) -> &'static str {
    match kind {
        "CDS" => "#ccffcc",
        "misc_feature" => "#99ccff",
        "primer_bind" => "#a020f0",
        "rep_origin" => "#ffff00",
        "protein_bind" => "#31849b",
        "misc_RNA" => "#00ccff",
        "rRNA" => "#00ff00",
        "tRNA" => "#228b22",
        "Region" => "#cc99b2",
        "Site" => NO_COLOR,
        "promoter" | "terminator" | "regulatory" | "source" => "#ffffff",
        _ => "#a6acb3",
    }
}

/// 特征的一个片段（1-based 闭区间，`start > end` = 跨原点）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatSeg {
    pub start: usize,
    pub end: usize,
    pub color: String,
    /// 该片段自己的 `translated="1"`。
    ///
    /// 真 SnapGene 里它是**按片段**存的：CDS 的 gap 片段（`type="gap"`）就没有它，
    /// 所以不能拿特征级开关去覆盖每个片段（那会把 gap 片段也改成翻译片段）。
    pub translated: bool,
}

/// 可编辑的特征（对话框的数据模型）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatSpec {
    pub recent_id: u32,
    pub name: String,
    pub kind: String,
    /// 0 = 无向、1 = 正向、2 = 反向、3 = 双向。
    pub directionality: i32,
    /// `translated="1"`：在序列视图里翻译此特征。
    pub translate: bool,
    /// `prioritize="1"`：在图中优先显示此特征。
    pub prioritize: bool,
    /// `visible="0"` = 隐藏。
    pub visible: bool,
    /// `/note` 限定符（去掉富文本外壳后的纯文本）。
    pub note: String,
    /// 其它限定符（对话框只读展示，写回时原样保留）。
    pub other_quals: Vec<(String, String)>,
    pub segments: Vec<FeatSeg>,
}

impl FeatSpec {
    /// 新建特征（SnapGene 的 `Add Feature`）：1 个片段、类型 `misc_feature`、
    /// 颜色取该类型的默认色。
    pub fn new(start: usize, end: usize) -> Self {
        Self {
            recent_id: 0,
            name: String::new(),
            kind: "misc_feature".to_string(),
            directionality: 1,
            translate: false,
            prioritize: false,
            visible: true,
            note: String::new(),
            other_quals: Vec::new(),
            segments: vec![FeatSeg {
                start,
                end,
                color: default_color("misc_feature").to_string(),
                translated: false,
            }],
        }
    }

    /// 特征级「翻译此特征」开关（SnapGene 的复选框）：写回时**每个片段一起**变。
    ///
    /// 读取时是「任一开心片段处于翻译状态」；编辑时统一设定所有片段。
    pub fn set_translate(&mut self, on: bool) {
        self.translate = on;
        for s in &mut self.segments {
            s.translated = on;
        }
    }

    /// 片段总长（跨原点片段按「到序列末尾 + 从头开始」无法在不知序列长度时计算，
    /// 这里返回 `start..=end` 的朴素跨度；调用方需要精确长度时用 [`Self::span_len`]）。
    pub fn span_len(&self, seq_len: usize) -> usize {
        self.segments
            .iter()
            .map(|s| {
                if s.start <= s.end {
                    s.end - s.start + 1
                } else {
                    (seq_len - s.start + 1) + s.end
                }
            })
            .sum()
    }

    /// 首个片段（对话框的「位置」列与图谱判定用）。
    pub fn first_seg(&self) -> Option<&FeatSeg> {
        self.segments.first()
    }
}

/// 一次特征编辑。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureEdit {
    /// 添加（`recent_id` 由写回时分配）。
    Add(FeatSpec),
    /// 修改（按 `recent_id` 定位；只改写有变化的字段）。
    Update(FeatSpec),
    /// 删除。
    Remove { recent_id: u32 },
    /// 显示 / 隐藏（`Show or Hide Feature`，只动 `visible` 属性）。
    SetVisible { recent_id: u32, visible: bool },
}

impl FeatureEdit {
    /// 状态栏/undo 用的人读描述。
    pub fn label(&self) -> String {
        match self {
            FeatureEdit::Add(s) => {
                if s.name.is_empty() {
                    format!("add {}", s.kind)
                } else {
                    format!("add feature {}", s.name)
                }
            }
            FeatureEdit::Update(s) => {
                if s.name.is_empty() {
                    format!("edit {}", s.kind)
                } else {
                    format!("edit feature {}", s.name)
                }
            }
            FeatureEdit::Remove { recent_id } => format!("remove feature #{recent_id}"),
            FeatureEdit::SetVisible { visible, .. } => {
                if *visible {
                    "show feature".to_string()
                } else {
                    "hide feature".to_string()
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 解析
// ---------------------------------------------------------------------------

/// XML 属性/文本值转义。
///
/// 严格跟真 SnapGene 对齐（全语料统计）：它只用 `&amp; &lt; &gt; &quot;` 与数字
/// 实体 `&#10;`（换行）；从不写 `&apos;`/`&#39;`——`'` 一律**原样**写在双引号属性里
/// （如 `name="3' LTR (ΔU3)"`）。
fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\n' => out.push_str("&#10;"),
            _ => out.push(c),
        }
    }
    out
}

/// 反转义实体。
pub fn unesc(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    s.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&#39;", "'")
        .replace("&#10;", "\n")
        .replace("&amp;", "&")
}

/// 只读展示用的富文本 → 纯文本：在 [`plain_text`] 之上把行标签换成换行。
///
/// oracle：SG 的特征工具卡里 `/evidence=…` 是**另起一行**的，而文件里的 note 是
/// `<html><body>N6-acetyllysine.<br>/evidence=…</body></html>`——`<br>` 必须变成
/// 换行，否则两行会粘在一起。这里只服务**显示**，写回路径仍用 [`plain_text`]
/// （它把 `<br>` 整个丢掉，保持与 SnapGene 对话框一致的纯文本）。
pub fn plain_text_lines(html: &str) -> String {
    if !html.contains('<') {
        return html.to_string();
    }
    let mut s = String::with_capacity(html.len() + 8);
    let mut i = 0usize;
    while i < html.len() {
        if html.as_bytes()[i] == b'<' {
            // 只在**确实有闭合 `>`** 时按标签处理（没闭合的 `<` 当普通字符）。
            let Some(rel) = html[i..].find('>') else {
                s.push_str(&html[i..]);
                break;
            };
            let end = i + rel;
            let tag = html[i + 1..end].trim_end_matches('/').trim().to_ascii_lowercase();
            if matches!(tag.as_str(), "br" | "p" | "div" | "li" | "tr" | "hr") {
                s.push('\n');
            } else {
                s.push_str(&html[i..=end]);
            }
            i = end + 1;
        } else {
            let ch = html[i..].chars().next().unwrap();
            s.push(ch);
            i += ch.len_utf8();
        }
    }
    plain_text(&s)
}

/// 去掉 SnapGene 富文本外壳（`<html><body>…</body></html>`）与 HTML 标签，
/// 给对话框显示纯文本。
pub fn plain_text(html: &str) -> String {
    let mut s = html.to_string();
    if let Some(i) = s.find("<body>") {
        s = s[i + "<body>".len()..].to_string();
        if let Some(j) = s.rfind("</body>") {
            s.truncate(j);
        }
    }
    let mut out = String::with_capacity(s.len());
    let mut depth = 0usize;
    for c in s.chars() {
        match c {
            '<' => depth += 1,
            '>' => depth = depth.saturating_sub(1),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out.replace("&nbsp;", " ")
        .replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
}

/// 一个 `<Feature>` 元素在文本中的字节范围与属性。
#[derive(Debug, Clone)]
struct Elem {
    start: usize,
    open_end: usize,
    end: usize,
    attrs: Vec<(String, String)>,
    inner: String,
}

impl Elem {
    fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

/// 解析开放标签 `attrs`（`key="value"` 序列）。
fn parse_attrs(s: &str) -> Vec<(String, String)> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < b.len() {
        while i < b.len() && (b[i] as char).is_whitespace() {
            i += 1;
        }
        let ks = i;
        while i < b.len() && b[i] != b'=' && !(b[i] as char).is_whitespace() {
            i += 1;
        }
        if ks == i {
            break;
        }
        let key = s[ks..i].to_string();
        while i < b.len() && (b[i] as char).is_whitespace() {
            i += 1;
        }
        if i < b.len() && b[i] == b'=' {
            i += 1;
            while i < b.len() && (b[i] as char).is_whitespace() {
                i += 1;
            }
            let (qs, qe) = if i < b.len() && (b[i] == b'"' || b[i] == b'\'') {
                let q = b[i];
                let s2 = i + 1;
                let mut j = s2;
                while j < b.len() && b[j] != q {
                    j += 1;
                }
                let r = (s2, j);
                i = (j + 1).min(b.len());
                r
            } else {
                let s2 = i;
                while i < b.len() && !(b[i] as char).is_whitespace() {
                    i += 1;
                }
                (s2, i)
            };
            out.push((key, s[qs..qe].to_string()));
        } else {
            out.push((key, String::new()));
        }
    }
    out
}

/// 扫描所有 `<Feature …>…</Feature>` 元素。
fn scan(text: &str) -> Vec<Elem> {
    let mut out = Vec::new();
    let mut at = 0usize;
    while let Some(rel) = text[at..].find("<Feature ") {
        let start = at + rel;
        let Some(gt) = text[start..].find('>') else {
            break;
        };
        let open_end = start + gt + 1;
        let open = &text[start + "<Feature".len()..open_end - 1];
        let attrs = parse_attrs(open);
        let self_closing = text[start..open_end].ends_with("/>");
        let (inner, end) = if self_closing {
            (String::new(), open_end)
        } else {
            match text[open_end..].find("</Feature>") {
                Some(r) => (
                    text[open_end..open_end + r].to_string(),
                    open_end + r + "</Feature>".len(),
                ),
                None => (text[open_end..].to_string(), text.len()),
            }
        };
        out.push(Elem {
            start,
            open_end,
            end,
            attrs,
            inner,
        });
        at = end.max(open_end);
    }
    out
}

/// 解析 `block 10` 的可编辑特征列表（按文件中出现顺序）。
pub fn parse_specs(xml: &[u8]) -> Vec<FeatSpec> {
    let Ok(text) = std::str::from_utf8(xml) else {
        return Vec::new();
    };
    scan(text).iter().filter_map(spec_from_elem).collect()
}

/// 从一个 `<Feature>` 元素解析出可编辑模型。
fn spec_from_elem(e: &Elem) -> Option<FeatSpec> {
    let recent_id: u32 = e.attr("recentID").and_then(|v| v.parse().ok())?;
    let segments = parse_segs(&e.inner);
    // 限定符：`<Q name="…"><V text="…|int="…"/></Q>`
    let mut note = String::new();
    let mut other = Vec::new();
    let mut qat = 0usize;
    while let Some(rel) = e.inner[qat..].find("<Q name=\"") {
        let s = qat + rel;
        let name_start = s + "<Q name=\"".len();
        let Some(qend) = e.inner[name_start..].find('"') else {
            break;
        };
        let qname = e.inner[name_start..name_start + qend].to_string();
        let close = e.inner[s..]
            .find("</Q>")
            .map(|i| s + i + "</Q>".len())
            .unwrap_or(e.inner.len());
        let val = q_value(&e.inner[s..close]);
        if qname == "note" {
            note = plain_text(&val);
        } else {
            other.push((qname, plain_text(&val)));
        }
        qat = close;
    }
    let dir = e
        .attr("directionality")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let translate = segments_meta(&e.inner).0;
    Some(FeatSpec {
        recent_id,
        name: unesc(e.attr("name").unwrap_or("")),
        kind: unesc(e.attr("type").unwrap_or("misc_feature")),
        directionality: dir,
        translate,
        prioritize: e.attr("prioritize") == Some("1"),
        visible: e.attr("visible") != Some("0"),
        note,
        other_quals: other,
        segments,
    })
}

/// 一个 `<Q>` 元素的取值（`<V text="…">` 优先，其次 `int="…"`）。
/// 从一个 `<Feature>` 元素内文里取 `/note` 限定符，转成**多行纯文本**（工具卡用）。
///
/// 与 [`parse_specs`] 里的 note 同源，但换行走 [`plain_text_lines`]（`<br>` → `\n`）。
pub fn note_of(inner: &str) -> String {
    let Some(rel) = inner.find("<Q name=\"note\">") else {
        return String::new();
    };
    let s = rel;
    let close = inner[s..]
        .find("</Q>")
        .map(|i| s + i + "</Q>".len())
        .unwrap_or(inner.len());
    plain_text_lines(&q_value(&inner[s..close]))
}

fn q_value(q: &str) -> String {
    let Some(i) = q.find("<V ") else {
        return String::new();
    };
    let v = &q[i..];
    if let Some(j) = v.find("text=\"") {
        let rest = &v[j + 6..];
        let e = rest.find('"').unwrap_or(rest.len());
        return unesc(&rest[..e]);
    }
    if let Some(j) = v.find("int=\"") {
        let rest = &v[j + 5..];
        let e = rest.find('"').unwrap_or(rest.len());
        return rest[..e].to_string();
    }
    String::new()
}

/// `(是否有 translated="1", note 原文, note 纯文本)` 之类 的内文查询。
fn segments_meta(inner: &str) -> (bool, Option<String>, String) {
    let translate = inner
        .split("<Segment ")
        .skip(1)
        .any(|s| s.split('>').next().unwrap_or("").contains("translated=\"1\""));
    let mut raw = None;
    let mut text = String::new();
    let mut at = 0usize;
    while let Some(rel) = inner[at..].find("<Q name=\"note\">") {
        let s = at + rel;
        let close = inner[s..]
            .find("</Q>")
            .map(|i| s + i + "</Q>".len())
            .unwrap_or(inner.len());
        raw = Some(inner[s..close].to_string());
        text = plain_text(&q_value(&inner[s..close]));
        at = close;
    }
    (translate, raw, text)
}


// ---------------------------------------------------------------------------
// 写回
// ---------------------------------------------------------------------------

/// 把属性表写成开放标签内文（`k="v" k2="v2"`，单空格分隔）。
fn render_attrs(attrs: &[(String, String)]) -> String {
    let mut s = String::new();
    for (k, v) in attrs {
        s.push(' ');
        s.push_str(k);
        s.push_str("=\"");
        s.push_str(&esc(v));
        s.push('"');
    }
    s
}

fn set_attr(attrs: &mut Vec<(String, String)>, key: &str, val: Option<&str>) {
    match val {
        Some(v) => match attrs.iter_mut().find(|(k, _)| k == key) {
            Some(e) => e.1 = v.to_string(),
            None => attrs.push((key.to_string(), v.to_string())),
        },
        None => attrs.retain(|(k, _)| k != key),
    }
}

/// 取原元素里的 `<Segment …/>` 原文标签（保留原属性，避免无谓改写）。
fn segment_tags(inner: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut at = 0usize;
    while let Some(rel) = inner[at..].find("<Segment ") {
        let s = at + rel;
        let Some(gt) = inner[s..].find('>') else {
            break;
        };
        out.push(inner[s..s + gt + 1].to_string());
        at = s + gt + 1;
    }
    out
}

/// 解析元素内文的片段列表。
fn parse_segs(inner: &str) -> Vec<FeatSeg> {
    let mut out = Vec::new();
    for tag in segment_tags(inner) {
        let a = parse_attrs(&tag["<Segment".len()..]);
        let get = |k: &str| {
            a.iter()
                .find(|(kk, _)| kk == k)
                .map(|(_, v)| v.clone())
                .unwrap_or_default()
        };
        if let Some((sa, sb)) = get("range").split_once('-') {
            if let (Ok(sa), Ok(sb)) = (sa.parse::<usize>(), sb.parse::<usize>()) {
                out.push(FeatSeg {
                    start: sa,
                    end: sb,
                    color: get("color"),
                    translated: get("translated") == "1",
                });
            }
        }
    }
    out
}

/// 改写一个 `<Segment>` 标签的 `translated` 属性。
///
/// SnapGene **只写 `translated="1"`**（6530 个特征 / 7606 个片段里 `translated="0"`
/// 出现 0 次）——默认（不翻译）时该属性整个省略，所以关闭翻译 = 删掉属性。
fn set_translated(tag: &str, translate: bool) -> String {
    if let Some(i) = tag.find("translated=\"") {
        let s = i + "translated=\"".len();
        if let Some(e) = tag[s..].find('"') {
            if translate {
                let mut out = String::with_capacity(tag.len());
                out.push_str(&tag[..s]);
                out.push('1');
                out.push_str(&tag[s + e..]);
                return out;
            }
            // 删掉整个属性（连同前面的空白；`s+e` 是闭合引号处，跳过它）
            let mut start = i;
            while start > 0 && tag.as_bytes()[start - 1].is_ascii_whitespace() {
                start -= 1;
            }
            let mut out = String::with_capacity(tag.len());
            out.push_str(&tag[..start]);
            out.push_str(&tag[(s + e + 1).min(tag.len())..]);
            return out;
        }
    }
    if !translate {
        return tag.to_string();
    }
    // 没有该属性：插到 `/>` 前
    match tag.rfind("/>") {
        Some(i) => format!("{} translated=\"1\"/>", tag[..i].trim_end()),
        None => tag.to_string(),
    }
}
/// 片段序列化；`orig` = 原元素（区间与颜色未变时原样保留原文，只按需改 `translated`）。
fn render_segments(spec: &FeatSpec, orig: Option<&Elem>) -> String {
    if let Some(o) = orig {
        let raw = segment_tags(&o.inner);
        let old = parse_segs(&o.inner);
        let same = !raw.is_empty()
            && old.len() == spec.segments.len()
            && old.iter().zip(&spec.segments).all(|(a, b)| {
                a.start == b.start && a.end == b.end && a.color == b.color
            });
        if same {
            // 逐片段对齐 `translated`：原文没动过的片段**逐字节保留**
            return raw
                .iter()
                .zip(&spec.segments)
                .map(|(t, seg)| {
                    let cur = t.contains("translated=\"1\"");
                    if cur == seg.translated {
                        t.clone()
                    } else {
                        set_translated(t, seg.translated)
                    }
                })
                .collect::<Vec<_>>()
                .concat();
        }
    }
    spec.segments
        .iter()
        .map(|seg| {
            // 与 SnapGene 一致：不翻译时省略 `translated`，只在翻译时写 `translated="1"`
            format!(
                "<Segment range=\"{}-{}\" color=\"{}\" type=\"standard\"{}/>",
                seg.start,
                seg.end,
                esc(&seg.color),
                if seg.translated { " translated=\"1\"" } else { "" }
            )
        })
        .collect()
}

/// `/note` 限定符元素（SnapGene 写纯文本，不套富文本外壳）。
fn render_note_elem(note: &str) -> String {
    format!("<Q name=\"note\"><V text=\"{}\"/></Q>", esc(note))
}

/// 原元素里除 `<Segment>`/`<Q>` 之外的残余内嵌内容（罕见，原样保留）。
fn residual(inner: &str) -> String {
    let mut rest = String::new();
    let mut i = 0usize;
    while i < inner.len() {
        let next_seg = inner[i..].find("<Segment ").map(|r| i + r);
        let next_q = inner[i..].find("<Q name=\"").map(|r| i + r);
        let next = match (next_seg, next_q) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => {
                rest.push_str(&inner[i..]);
                break;
            }
        };
        rest.push_str(&inner[i..next]);
        if Some(next) == next_seg {
            i = inner[next..]
                .find("/>")
                .map(|r| next + r + 2)
                .unwrap_or(inner.len());
        } else {
            i = inner[next..]
                .find("</Q>")
                .map(|r| next + r + 4)
                .unwrap_or(inner.len());
        }
    }
    if rest.trim().is_empty() {
        String::new()
    } else {
        rest
    }
}

/// 渲染**一条**限定符元素（新值）。
///
/// 数字且原文用 `int=` 形式（或属于惯用 `int` 的限定符）时保持 `int=`，否则用 `text=`。
fn render_qual(name: &str, val: &str, orig: Option<&str>) -> String {
    let numeric = val.parse::<i64>().is_ok()
        && match orig {
            Some(e) => e.contains("<V int=\""),
            None => matches!(
                name,
                "codon_start" | "transl_table" | "citation" | "calculated_mol_wt"
            ),
        };
    format!(
        "<Q name=\"{}\"><V {}=\"{}\"/></Q>",
        esc(name),
        if numeric { "int" } else { "text" },
        esc(val)
    )
}

/// 渲染限定符：**按原文顺序**逐个判定，值未变 ⇒ 原元素逐字节保留（可能带 `<html>` 外壳
/// 或 `int=` 形式）；值被改过 ⇒ 重写；原文没有而 spec 有的（新加）⇒ 追加在末尾。
///
/// 顺序很要紧：真 SnapGene 的 `note` 并不总在最后（6530 个特征里 1470 个 note 在末尾、
/// 其余夹在中间），所以只能跟着原文走。
fn render_quals(spec: &FeatSpec, orig: Option<&Elem>) -> String {
    let mut out = String::new();
    let Some(o) = orig else {
        for (n, v) in &spec.other_quals {
            out.push_str(&render_qual(n, v, None));
        }
        if !spec.note.is_empty() {
            out.push_str(&render_note_elem(&spec.note));
        }
        return out;
    };
    let mut seen: Vec<String> = Vec::new();
    let mut at = 0usize;
    while let Some(rel) = o.inner[at..].find("<Q name=\"") {
        let s = at + rel;
        let name_start = s + "<Q name=\"".len();
        let Some(qend) = o.inner[name_start..].find('"') else {
            break;
        };
        let qname = o.inner[name_start..name_start + qend].to_string();
        let close = o.inner[s..]
            .find("</Q>")
            .map(|i| s + i + "</Q>".len())
            .unwrap_or(o.inner.len());
        let elem = &o.inner[s..close];
        at = close;
        if qname == "note" {
            let same = plain_text(&q_value(elem)) == spec.note;
            if same {
                // 含富文本外壳；`<html><body><br></body></html>`（纯文本为空）也原样保留
                out.push_str(elem);
            } else if !spec.note.is_empty() {
                out.push_str(&render_note_elem(&spec.note));
            }
            // 原文有、纯文本不同且新值为空 ⇒ 用户清空了 note，丢弃
        } else if let Some((n, v)) = spec.other_quals.iter().find(|(n, _)| *n == qname) {
            seen.push(n.clone());
            if plain_text(&q_value(elem)) == *v {
                out.push_str(elem);
            } else {
                out.push_str(&render_qual(n, v, Some(elem)));
            }
        }
        // 原文有、spec 里没有 ⇒ 用户删掉了
    }
    for (n, v) in &spec.other_quals {
        if !seen.contains(n) {
            out.push_str(&render_qual(n, v, None));
        }
    }
    if !o.inner.contains("<Q name=\"note\"") && !spec.note.is_empty() {
        out.push_str(&render_note_elem(&spec.note));
    }
    out
}

/// 把一个特征渲染成 XML 元素；`orig` = 原元素（保留未触及的属性与限定符）。
fn render_feature(id: u32, spec: &FeatSpec, orig: Option<&Elem>, raw: Option<&str>) -> String {
    // 快速通路：可编辑字段与原文完全一致 ⇒ 逐字节返回原文
    if let (Some(o), Some(t)) = (orig, raw) {
        if let Some(orig_spec) = spec_from_elem(o) {
            if orig_spec.name == spec.name
                && orig_spec.kind == spec.kind
                && orig_spec.directionality == spec.directionality
                && orig_spec.translate == spec.translate
                && orig_spec.prioritize == spec.prioritize
                && orig_spec.visible == spec.visible
                && orig_spec.note == spec.note
                && orig_spec.other_quals == spec.other_quals
                && orig_spec.segments == spec.segments
            {
                return t.to_string();
            }
        }
    }
    let mut attrs: Vec<(String, String)> = orig.map(|e| e.attrs.clone()).unwrap_or_default();
    set_attr(&mut attrs, "recentID", Some(&id.to_string()));
    set_attr(&mut attrs, "name", Some(&spec.name));
    set_attr(&mut attrs, "type", Some(&spec.kind));
    if spec.directionality == 0 {
        set_attr(&mut attrs, "directionality", None);
    } else {
        set_attr(
            &mut attrs,
            "directionality",
            Some(&spec.directionality.to_string()),
        );
    }
    set_attr(&mut attrs, "allowSegmentOverlaps", Some("0"));
    set_attr(
        &mut attrs,
        "prioritize",
        if spec.prioritize { Some("1") } else { None },
    );
    set_attr(
        &mut attrs,
        "visible",
        if spec.visible { None } else { Some("0") },
    );
    if orig.is_none()
        && !attrs
            .iter()
            .any(|(k, _)| k == "consecutiveTranslationNumbering")
    {
        attrs.push(("consecutiveTranslationNumbering".into(), "1".into()));
    }
    let mut inner = render_segments(spec, orig);
    inner.push_str(&render_quals(spec, orig));
    if let Some(o) = orig {
        inner.push_str(&residual(&o.inner));
    }
    format!("<Feature{}>{}</Feature>", render_attrs(&attrs), inner)
}

/// 在 `block 10` 载荷上顺序施加特征编辑，返回新的载荷。
pub fn apply_feature_edits(xml: &[u8], edits: &[FeatureEdit]) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(xml) else {
        return xml.to_vec();
    };
    let mut cur = text.to_string();
    for e in edits {
        cur = apply_one(&cur, e);
    }
    cur.into_bytes()
}

fn next_valid_id(text: &str) -> u32 {
    let Some(i) = text.find("nextValidID=\"") else {
        return 0;
    };
    let s = i + "nextValidID=\"".len();
    text[s..]
        .split('"')
        .next()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

fn apply_one(text: &str, edit: &FeatureEdit) -> String {
    let elems = scan(text);
    match edit {
        FeatureEdit::Add(spec) => {
            let max_id = elems
                .iter()
                .filter_map(|e| e.attr("recentID").and_then(|v| v.parse::<u32>().ok()))
                .max()
                .unwrap_or(0);
            let id = next_valid_id(text).max(max_id + 1);
            let elem = render_feature(id, spec, None, None);
            // 插到最后一个 </Feature> 之后（保持文件顺序：新特征在末尾）
            let insert_at = match text.rfind("</Feature>") {
                Some(i) => i + "</Feature>".len(),
                None => match text.rfind("</Features>") {
                    Some(i) => i,
                    None => text.len(),
                },
            };
            let mut out = String::with_capacity(text.len() + elem.len() + 8);
            out.push_str(&text[..insert_at]);
            out.push_str(&elem);
            out.push_str(&text[insert_at..]);
            // 更新 nextValidID
            if let Some(i) = out.find("nextValidID=\"") {
                let s = i + "nextValidID=\"".len();
                if let Some(e) = out[s..].find('"') {
                    let new_id = (id + 1).to_string();
                    out.replace_range(s..s + e, &new_id);
                }
            }
            out
        }
        FeatureEdit::Update(spec) => {
            let Some(e) = elems.iter().find(|e| {
                e.attr("recentID").and_then(|v| v.parse::<u32>().ok()) == Some(spec.recent_id)
            }) else {
                return text.to_string();
            };
            let elem = render_feature(spec.recent_id, spec, Some(e), Some(&text[e.start..e.end]));
            let mut out = String::with_capacity(text.len() + elem.len());
            out.push_str(&text[..e.start]);
            out.push_str(&elem);
            out.push_str(&text[e.end..]);
            out
        }
        FeatureEdit::Remove { recent_id } => {
            let Some(e) = elems.iter().find(|e| {
                e.attr("recentID").and_then(|v| v.parse::<u32>().ok()) == Some(*recent_id)
            }) else {
                return text.to_string();
            };
            let mut out = String::with_capacity(text.len());
            out.push_str(&text[..e.start]);
            out.push_str(&text[e.end..]);
            out
        }
        FeatureEdit::SetVisible { recent_id, visible } => {
            let Some(e) = elems.iter().find(|e| {
                e.attr("recentID").and_then(|v| v.parse::<u32>().ok()) == Some(*recent_id)
            }) else {
                return text.to_string();
            };
            let mut attrs = e.attrs.clone();
            set_attr(
                &mut attrs,
                "visible",
                if *visible { None } else { Some("0") },
            );
            let mut out = String::with_capacity(text.len() + 16);
            out.push_str(&text[..e.start]);
            out.push_str("<Feature");
            out.push_str(&render_attrs(&attrs));
            out.push('>');
            out.push_str(&text[e.open_end..]);
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const XML: &str = "<?xml version=\"1.0\"?><Features nextValidID=\"2\">\
<Feature recentID=\"0\" name=\"A\" type=\"CDS\" allowSegmentOverlaps=\"0\" consecutiveTranslationNumbering=\"1\">\
<Segment range=\"3-8\" color=\"#ccffcc\" type=\"standard\" translated=\"1\"/>\
<Q name=\"note\"><V text=\"hello\"/></Q></Feature>\
<Feature recentID=\"1\" name=\"B\" directionality=\"2\" type=\"promoter\" allowSegmentOverlaps=\"0\">\
<Segment range=\"12-20\" color=\"#ffffff\" type=\"standard\"/></Feature></Features>";

    /// 工具卡的 note：`<br>` 要变成换行（oracle：SG 卡片里 `/evidence=…` 另起一行），
    /// 但对话框用的 `plain_text` 仍然把 `<br>` 整个丢掉（写回路径不能变）。
    #[test]
    fn note_of_keeps_line_breaks_from_br() {
        let inner = "<Segment range=\"76-76\" color=\"noColor\" type=\"standard\"/>\
<Q name=\"note\"><V text=\"&lt;html&gt;&lt;body&gt;N6-acetyllysine.&lt;br&gt;/evidence=ECO:0000269.&lt;/body&gt;&lt;/html&gt;\"/></Q>";
        assert_eq!(note_of(inner), "N6-acetyllysine.\n/evidence=ECO:0000269.");
        // 同一个 note 走对话框路径仍然粘成一行（`plain_text` 不变）
        assert_eq!(plain_text("N6-acetyllysine.<br>/evidence=ECO:0000269."), "N6-acetyllysine./evidence=ECO:0000269.");
        // 没有 note 限定符 ⇒ 空串（卡片不画 NOTE 段）
        assert_eq!(note_of("<Segment range=\"1-2\" color=\"#fff\" type=\"standard\"/>"), "");
    }

    #[test]
    fn parse_reads_ids_segments_note_translate() {
        let s = parse_specs(XML.as_bytes());
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].recent_id, 0);
        assert_eq!(s[0].name, "A");
        assert_eq!(s[0].kind, "CDS");
        assert!(s[0].translate);
        assert_eq!(s[0].note, "hello");
        assert_eq!(s[0].segments[0].start, 3);
        assert_eq!(s[0].segments[0].end, 8);
        assert_eq!(s[1].directionality, 2);
        assert!(!s[1].translate);
    }

    #[test]
    fn update_unchanged_is_byte_identical() {
        let specs = parse_specs(XML.as_bytes());
        for sp in specs {
            let out = apply_feature_edits(XML.as_bytes(), &[FeatureEdit::Update(sp)]);
            assert_eq!(String::from_utf8(out).unwrap(), XML);
        }
    }

    #[test]
    fn add_assigns_next_id_and_bumps_header() {
        let mut spec = FeatSpec::new(30, 40);
        spec.name = "new one".into();
        spec.kind = "misc_feature".into();
        let out = apply_feature_edits(XML.as_bytes(), &[FeatureEdit::Add(spec)]);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("nextValidID=\"3\""), "{s}");
        assert!(s.contains("<Feature recentID=\"2\" name=\"new one\" type=\"misc_feature\""), "{s}");
        assert!(s.contains("<Segment range=\"30-40\""), "{s}");
        // 已有元素逐字节保留
        assert!(s.contains("<Feature recentID=\"1\" name=\"B\" directionality=\"2\" type=\"promoter\" allowSegmentOverlaps=\"0\">"));
        let again = parse_specs(s.as_bytes());
        assert_eq!(again.len(), 3);
        assert_eq!(again[2].name, "new one");
    }

    #[test]
    fn update_changes_only_touched_fields() {
        let mut sp = parse_specs(XML.as_bytes())[1].clone();
        sp.name = "B2".into();
        sp.directionality = 0;
        let out = apply_feature_edits(XML.as_bytes(), &[FeatureEdit::Update(sp)]);
        let s = String::from_utf8(out).unwrap();
        // 名字改了、directionality 掉了，其余原样
        assert!(s.contains("<Feature recentID=\"1\" name=\"B2\" type=\"promoter\" allowSegmentOverlaps=\"0\">"), "{s}");
        assert!(s.contains("<Feature recentID=\"0\" name=\"A\" type=\"CDS\" allowSegmentOverlaps=\"0\" consecutiveTranslationNumbering=\"1\">"));
    }

    #[test]
    fn remove_and_visibility() {
        let s = apply_feature_edits(XML.as_bytes(), &[FeatureEdit::Remove { recent_id: 0 }]);
        let s = String::from_utf8(s).unwrap();
        assert!(!s.contains("name=\"A\""));
        assert!(s.contains("name=\"B\""));

        let s = apply_feature_edits(
            XML.as_bytes(),
            &[FeatureEdit::SetVisible { recent_id: 1, visible: false }],
        );
        let s = String::from_utf8(s).unwrap();
        assert!(s.contains("<Feature recentID=\"1\" name=\"B\" directionality=\"2\" type=\"promoter\" allowSegmentOverlaps=\"0\" visible=\"0\">"), "{s}");
        // 再显示回来 → 与原 XML 逐字节一致
        let s2 = apply_feature_edits(
            s.as_bytes(),
            &[FeatureEdit::SetVisible { recent_id: 1, visible: true }],
        );
        assert_eq!(String::from_utf8(s2).unwrap(), XML);
    }

    #[test]
    fn plain_text_strips_rich_wrapper() {
        assert_eq!(
            plain_text("<html><body>direct repeat</body></html>"),
            "direct repeat"
        );
        assert_eq!(
            plain_text("<html><body>x<i>i</i>y</body></html>"),
            "xiy"
        );
    }

    /// 多段 + 富文本限定符的样本（字段顺序、`int=`/富文本外壳都按真 SnapGene 来）。
    const XML_RICH: &str = "<?xml version=\"1.0\"?><Features nextValidID=\"2\">\
<Feature recentID=\"0\" name=\"gfp\" directionality=\"1\" type=\"CDS\" allowSegmentOverlaps=\"0\" consecutiveTranslationNumbering=\"1\">\
<Segment range=\"1-3\" color=\"#ccffcc\" type=\"standard\" translated=\"1\"/>\
<Segment range=\"4-6\" color=\"#ccffcc\" type=\"standard\" translated=\"1\"/>\
<Q name=\"codon_start\"><V int=\"1\"/></Q>\
<Q name=\"product\"><V text=\"&lt;html&gt;&lt;body&gt;enhanced GFP&lt;/body&gt;&lt;/html&gt;\"/></Q>\
<Q name=\"transl_table\"><V int=\"1\"/></Q>\
<Q name=\"note\"><V text=\"&lt;html&gt;&lt;body&gt;codon optimized&lt;/body&gt;&lt;/html&gt;\"/></Q></Feature>\
<Feature recentID=\"1\" name=\"ori\" type=\"rep_origin\" allowSegmentOverlaps=\"0\">\
<Segment range=\"10-20\" color=\"#ffff00\" type=\"standard\"/></Feature></Features>";

    #[test]
    fn segment_edit_keeps_other_qualifiers_byte_exact() {
        let mut sp = parse_specs(XML_RICH.as_bytes())[0].clone();
        assert_eq!(sp.segments.len(), 2);
        sp.segments.pop();
        let out = String::from_utf8(apply_feature_edits(
            XML_RICH.as_bytes(),
            &[FeatureEdit::Update(sp)],
        ))
        .unwrap();
        assert!(
            out.contains("<Segment range=\"1-3\" color=\"#ccffcc\" type=\"standard\" translated=\"1\"/>"),
            "{out}"
        );
        assert!(!out.contains("range=\"4-6\""), "{out}");
        // 未动的限定符逐字节保留（富文本外壳、int 形式都在）
        assert!(out.contains(
            "<Q name=\"product\"><V text=\"&lt;html&gt;&lt;body&gt;enhanced GFP&lt;/body&gt;&lt;/html&gt;\"/></Q>"
        ), "{out}");
        assert!(out.contains("<Q name=\"codon_start\"><V int=\"1\"/></Q>"), "{out}");
        assert!(out.contains(
            "<Q name=\"note\"><V text=\"&lt;html&gt;&lt;body&gt;codon optimized&lt;/body&gt;&lt;/html&gt;\"/></Q>"
        ), "{out}");
        // 另一个特征完全没动
        assert!(out.ends_with(
            "<Feature recentID=\"1\" name=\"ori\" type=\"rep_origin\" allowSegmentOverlaps=\"0\">\
<Segment range=\"10-20\" color=\"#ffff00\" type=\"standard\"/></Feature></Features>"
        ), "{out}");
    }

    #[test]
    fn editing_other_qualifier_rewrites_value() {
        let mut sp = parse_specs(XML_RICH.as_bytes())[0].clone();
        assert_eq!(sp.other_quals[1].0, "product");
        assert_eq!(sp.other_quals[1].1, "enhanced GFP");
        sp.other_quals[1].1 = "GFPmut3".into();
        let out = String::from_utf8(apply_feature_edits(
            XML_RICH.as_bytes(),
            &[FeatureEdit::Update(sp)],
        ))
        .unwrap();
        assert!(
            out.contains("<Q name=\"product\"><V text=\"GFPmut3\"/></Q>"),
            "{out}"
        );
        // 未改的 int 限定符仍是 int 形式、逐字节不变
        assert!(out.contains("<Q name=\"transl_table\"><V int=\"1\"/></Q>"), "{out}");
    }

    #[test]
    fn removing_other_qualifier_drops_it() {
        let mut sp = parse_specs(XML_RICH.as_bytes())[0].clone();
        sp.other_quals.retain(|(n, _)| n != "product");
        let out = String::from_utf8(apply_feature_edits(
            XML_RICH.as_bytes(),
            &[FeatureEdit::Update(sp)],
        ))
        .unwrap();
        assert!(!out.contains("name=\"product\""), "{out}");
        assert!(out.contains("name=\"codon_start\""), "{out}");
        assert!(out.contains("name=\"transl_table\""), "{out}");
    }

    #[test]
    fn new_feature_writes_translated_only_when_translated() {
        let mut sp = FeatSpec::new(1, 3);
        sp.name = "x".into();
        let out = String::from_utf8(apply_feature_edits(
            XML.as_bytes(),
            &[FeatureEdit::Add(sp.clone())],
        ))
        .unwrap();
        assert!(
            out.contains("<Segment range=\"1-3\" color=\"#99ccff\" type=\"standard\"/>"),
            "{out}"
        );
        sp.set_translate(true);
        let out2 = String::from_utf8(apply_feature_edits(XML.as_bytes(), &[FeatureEdit::Add(sp)]))
            .unwrap();
        assert!(
            out2.contains("<Segment range=\"1-3\" color=\"#99ccff\" type=\"standard\" translated=\"1\"/>"),
            "{out2}"
        );
    }

    #[test]
    fn turning_translation_off_removes_attribute() {
        let mut sp = parse_specs(XML_RICH.as_bytes())[0].clone();
        sp.set_translate(false);
        let out = String::from_utf8(apply_feature_edits(
            XML_RICH.as_bytes(),
            &[FeatureEdit::Update(sp)],
        ))
        .unwrap();
        assert!(
            out.contains("<Segment range=\"1-3\" color=\"#ccffcc\" type=\"standard\"/>"),
            "{out}"
        );
        assert!(
            out.contains("<Segment range=\"4-6\" color=\"#ccffcc\" type=\"standard\"/>"),
            "{out}"
        );
    }
}
