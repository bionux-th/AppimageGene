//! `block 7`（历史树）写出：明文 XML。
//!
//! 依据（oracle T3/T8）：
//!
//! - 历史面板显示**只由 block 7 驱动**；block 11 仅服务“回到历史点”；
//! - 明文 XML 树（非 XZ）被真 SnapGene 正常读取——随安装包分发的
//!   `Azurite BFP.dna` / `mTagBFP2.dna` / `pHELLSGATE 8.dna` / `pUC57-T.dna`
//!   四例均是明文树；
//! - 结构：`<HistoryTree><Node …>…`，**根节点 = 当前状态**，子节点 = 更早状态；
//! - 非叶节点必须有 `<InputSummary/>`（缺了会让真 SnapGene 崩溃；空元素合法）；
//! - 实测块序：block 7 紧跟 block 6 之后（`9,0,3,[11…],8,10,5,6,7,17,…`）。
//!
//! 我们写出的树：根 = 编辑后状态（`operation` 取 SnapGene 词汇表中的
//! `replace`/`insertFragment`），子 = 编辑前状态（`operation="invalid"`，
//! 携带当时的特征/引物 XML），从而在真 SnapGene 的 History 面板里看到
//! “本轮编辑”这一条记录。

use crate::dna::DnaFile;
use crate::error::FormatError;
use crate::util::xml_escape;
use crate::write::Edit;
use hx_core::{MoleculeFlags, MoleculeKind, SeqRecord, Strandedness};

/// A display node. Parent indices retain branching and duplicate names/IDs.
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryEntry {
    pub attrs: Vec<(String, String)>,
    pub summaries: Vec<InputSummary>,
    pub parent: Option<usize>,
    pub features: Vec<crate::annot::FeatureAnnotation>,
    pub primers: Vec<crate::annot::PrimerAnnotation>,
}

impl HistoryEntry {
    pub fn attr(&self, key: &str) -> &str {
        self.attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str()).unwrap_or("")
    }

    pub fn label(&self) -> &str {
        if self.attr("useCustomMapLabel") == "1" && !self.attr("customMapLabel").is_empty() {
            self.attr("customMapLabel")
        } else {
            self.attr("name").trim_end_matches(".dna").trim_end_matches(".rna").trim_end_matches(".prot")
        }
    }
}

/// Block 7 is either UTF-8 XML or an XZ container (both occur in SG samples).
pub fn decode_payload(payload: &[u8]) -> Result<String, String> {
    let bytes = if payload.starts_with(b"\xfd7zXZ\0") {
        crate::xz::decompress(payload).map_err(|e| format!("Cannot decode history: {e}"))?
    } else { payload.to_vec() };
    String::from_utf8(bytes).map_err(|e| format!("History is not UTF-8: {e}"))
}

/// Read node hierarchy and all input summaries, ignoring annotation subtrees.
/// Quote-aware tokenization prevents attribute text from becoming markup.
pub fn parse_entries(xml: &str) -> Result<Vec<HistoryEntry>, String> {
    let mut entries = Vec::new();
    let mut stack = Vec::new();
    let mut pos = 0;
    if !xml.contains("<HistoryTree") { return Err("Missing HistoryTree element".into()); }
    while let Some(rel) = xml[pos..].find('<') {
        let start = pos + rel;
        if xml[start..].starts_with("<!--") {
            pos = start + xml[start..].find("-->").ok_or("Unclosed history comment")? + 3;
            continue;
        }
        let mut quote = None;
        let end = xml.as_bytes()[start..].iter().position(|&b| {
            if let Some(q) = quote { if b == q { quote = None; } false }
            else if b == b'\'' || b == b'"' { quote = Some(b); false }
            else { b == b'>' }
        }).map(|i| start + i).ok_or("Unclosed history tag")?;
        let tag = &xml[start..=end];
        let name = tag[1..].trim_start().split(|c: char| c.is_whitespace() || c == '/' || c == '>').next().unwrap_or("");
        if tag.starts_with("</Node") {
            stack.pop().ok_or("Unexpected history node closing tag")?;
        } else if name == "Node" {
            let attrs = crate::annot::next_tag(tag, "Node", 0).ok_or("Invalid history node")?.0
                .into_iter().map(|(k, v)| (k.to_owned(), v)).collect();
            let idx = entries.len();
            entries.push(HistoryEntry { attrs, summaries: Vec::new(), parent: stack.last().copied(), features: Vec::new(), primers: Vec::new() });
            if !tag.trim_end_matches('>').trim_end().ends_with('/') {
                if stack.len() >= 256 { return Err("History nesting exceeds 256 states".into()); }
                stack.push(idx);
            }
        } else if name == "Features" {
            if let Some(&idx) = stack.last() {
                if let Some(close) = xml[end + 1..].find("</Features>") {
                    entries[idx].features = crate::annot::parse_features(xml[start..end + 1 + close + 11].as_bytes());
                }
            }
        } else if name == "Primers" {
            if let Some(&idx) = stack.last() {
                if let Some(close) = xml[end + 1..].find("</Primers>") {
                    entries[idx].primers = crate::annot::parse_primers(xml[start..end + 1 + close + 10].as_bytes());
                }
            }
        } else if name == "InputSummary" {
            if let Some(&idx) = stack.last() {
                let attrs = crate::annot::next_tag(tag, "InputSummary", 0).ok_or("Invalid history summary")?.0
                    .into_iter().map(|(k, v)| (k.to_owned(), v)).collect();
                entries[idx].summaries.push(InputSummary { attrs });
            }
        }
        pos = end + 1;
    }
    if !stack.is_empty() { return Err("Unclosed history node".into()); }
    if entries.iter().filter(|n| n.parent.is_none()).count() != 1 { return Err("History must have one current state".into()); }
    Ok(entries)
}

/// English operation captions; unknown SG operations remain visible, never discarded.
pub fn operation_label(op: &str) -> String {
    let known = match op {
        "invalid" | "" => return String::new(),
        "insertFragment" => "Insert fragment",
        "insertFragments" => "Insert fragments",
        "goldenGateAssembly" => "Golden Gate assembly",
        "amplifyFragment" => "PCR amplification",
        "newFileFromSelection" => "New file from selection",
        "convertToSingleStranded" => "Convert to single-stranded",
        "convertToDoubleStranded" => "Convert to double-stranded",
        "flip" => "Reverse complement",
        "changeMethylation" => "Change methylation",
        "primerDirectedMutagenesis" => "Primer-directed mutagenesis",
        "replace" => "Replace sequence",
        "delete" => "Delete sequence",
        "insert" | "insertBases" => "Insert sequence",
        other => return other.to_owned(),
    };
    known.into()
}

/// `block 8` 里的末端修饰与粘性末端。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqProperties {
    pub upstream_stickiness: i64,
    pub downstream_stickiness: i64,
    pub upstream_modification: String,
    pub downstream_modification: String,
}

impl Default for SeqProperties {
    fn default() -> Self {
        Self {
            upstream_stickiness: 0,
            downstream_stickiness: 0,
            upstream_modification: "Unmodified".into(),
            downstream_modification: "Unmodified".into(),
        }
    }
}

/// `<InputSummary …/>` 的属性表（顺序保留）。
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct InputSummary {
    pub attrs: Vec<(String, String)>,
}

impl InputSummary {
    /// 空元素 `<InputSummary/>`。
    pub fn empty() -> Self {
        Self::default()
    }

    /// `manipulation="…" val1="…" val2="…"`。
    pub fn range(manipulation: &str, val1: usize, val2: usize) -> Self {
        Self {
            attrs: vec![
                ("manipulation".into(), manipulation.into()),
                ("val1".into(), val1.to_string()),
                ("val2".into(), val2.to_string()),
            ],
        }
    }

    fn render(&self, out: &mut String) {
        if self.attrs.is_empty() {
            out.push_str("<InputSummary/>");
            return;
        }
        out.push_str("<InputSummary");
        for (k, v) in &self.attrs {
            out.push_str(&format!(" {k}=\"{}\"", xml_escape(v)));
        }
        out.push_str("/>");
    }
}

/// 历史树节点。
#[derive(Debug, Clone)]
pub struct Node {
    pub name: String,
    pub kind: MoleculeKind,
    pub seq_len: usize,
    pub strandedness: Strandedness,
    pub id: u32,
    pub circular: bool,
    /// 可“回到”标记；实测 1/2（语义留给激活版实验），我们不写 block 11 时通常留空。
    pub resurrectable: Option<u8>,
    /// SnapGene 操作名（起始状态的父节点用 `invalid`）。
    pub operation: String,
    pub properties: Option<SeqProperties>,
    pub custom_map_label: Option<String>,
    pub input_summary: Option<InputSummary>,
    /// 原样嵌入的子 XML（`<Features>…</Features><Primers>…</Primers>`）。
    pub body: String,
    pub children: Vec<Node>,
    /// 原样嵌入的完整子元素（如从旧树移植来的父节点 `<Node>…</Node>`）。
    pub raw_children: Vec<String>,
}

impl Node {
    fn render(&self, out: &mut String) {
        out.push_str(&format!("<Node name=\"{}\"", xml_escape(&self.name)));
        out.push_str(&format!(" type=\"{}\"", kind_name(self.kind)));
        out.push_str(&format!(" seqLen=\"{}\"", self.seq_len));
        out.push_str(&format!(
            " strandedness=\"{}\"",
            match self.strandedness {
                Strandedness::Double => "double",
                Strandedness::Single => "single",
            }
        ));
        out.push_str(&format!(" ID=\"{}\" circular=\"{}\"", self.id, self.circular as u8));
        if let Some(r) = self.resurrectable {
            out.push_str(&format!(" resurrectable=\"{r}\""));
        }
        if let Some(p) = &self.properties {
            out.push_str(&format!(
                " upstreamModification=\"{}\" downstreamModification=\"{}\"",
                xml_escape(&p.upstream_modification),
                xml_escape(&p.downstream_modification)
            ));
            if p.upstream_stickiness != 0 {
                out.push_str(&format!(" upstreamStickiness=\"{}\"", p.upstream_stickiness));
            }
            if p.downstream_stickiness != 0 {
                out.push_str(&format!(
                    " downstreamStickiness=\"{}\"",
                    p.downstream_stickiness
                ));
            }
        }
        if let Some(l) = &self.custom_map_label {
            out.push_str(&format!(
                " useCustomMapLabel=\"1\" customMapLabel=\"{}\"",
                xml_escape(l)
            ));
        }
        out.push_str(&format!(" operation=\"{}\"", self.operation));

        let has_body = !self.body.is_empty()
            || !self.children.is_empty()
            || !self.raw_children.is_empty()
            || self.input_summary.is_some();
        if !has_body {
            out.push_str("/>");
            return;
        }
        out.push('>');
        match &self.input_summary {
            Some(is) => is.render(out),
            // 非叶节点必须有 InputSummary（缺了会让真 SnapGene 崩溃）
            None => out.push_str("<InputSummary/>"),
        }
        out.push_str(&self.body);
        for c in &self.children {
            c.render(out);
        }
        for r in &self.raw_children {
            out.push_str(r);
        }
        out.push_str("</Node>");
    }
}

/// `<HistoryTree>` 包装，输出 block 7 载荷字节（明文 XML，UTF-8）。
pub fn to_block_payload(root: &Node) -> Vec<u8> {
    let mut s = String::with_capacity(512);
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?><HistoryTree>");
    root.render(&mut s);
    s.push_str("</HistoryTree>");
    s.into_bytes()
}

fn kind_name(kind: MoleculeKind) -> &'static str {
    match kind {
        MoleculeKind::Dna => "DNA",
        MoleculeKind::Rna => "RNA",
        MoleculeKind::Protein => "Protein",
    }
}

// ---------------------------------------------------------------------------
// 从编辑前后文件构造历史树
// ---------------------------------------------------------------------------

/// 一个节点的基础属性（供根/父两个状态共用）。
fn base_node(
    name: &str,
    rec: &SeqRecord,
    seq_len: usize,
    id: u32,
    operation: &str,
    properties: Option<SeqProperties>,
    custom_map_label: Option<String>,
) -> Node {
    Node {
        name: name.to_string(),
        kind: rec.kind,
        seq_len,
        strandedness: rec.flags.strandedness(),
        id,
        circular: rec.flags.has(MoleculeFlags::CIRCULAR),
        resurrectable: None,
        operation: operation.to_string(),
        properties,
        custom_map_label,
        input_summary: None,
        body: String::new(),
        children: Vec::new(),
        raw_children: Vec::new(),
    }
}

/// 构造“本轮编辑”的历史树：根 = 编辑后状态，子 = 编辑前状态。
///
/// - `new_file`：`rewrite_after_edit` 的产物（含调整后的块 10/5，用于根节点 body）；
/// - `src`：编辑前的文件（块 10/5 用于父节点 body，块 6/8 提供标签与末端修饰）；
/// - `document_name` / `source_name`：根/父节点显示名（通常是两个文件名）；
/// - `parent_resurrectable`：父节点是否标记 `resurrectable`（我们未写 block 11，
///   默认 `None`；oracle 实验可对比写法）。
pub fn edit_history(
    new_file: &DnaFile,
    src: &DnaFile,
    edits: &[Edit],
    document_name: &str,
    source_name: &str,
    parent_resurrectable: Option<u8>,
) -> Result<Node, FormatError> {
    if let Some(b) = src.block(7) {
        let text = decode_payload(&b.payload).map_err(FormatError::BadHeader)?;
        parse_entries(&text).map_err(FormatError::BadHeader)?;
    }
    let graft = existing_root_xml(src);
    if src.header.molecule_kind().is_none() {
        return Err(FormatError::BadHeader("未知分子类型".into()));
    }
    let old_rec = src
        .sequence()
        .ok_or_else(|| FormatError::BadHeader("缺少序列块".into()))?;
    let new_rec = new_file
        .sequence()
        .ok_or_else(|| FormatError::BadHeader("缺少序列块".into()))?;

    let notes = src
        .block(6)
        .and_then(|b| std::str::from_utf8(&b.payload).ok());
    let label = notes.and_then(|x| {
        let use_label = tag_text(x, "UseCustomMapLabel").unwrap_or_default() == "1";
        if use_label {
            tag_text(x, "CustomMapLabel").filter(|l| !l.is_empty())
        } else {
            None
        }
    });
    let properties = src
        .block(8)
        .and_then(|b| std::str::from_utf8(&b.payload).ok())
        .map(parse_properties)
        .or_else(|| Some(SeqProperties::default()));

    let (operation, summary) = summarize_edits(edits);

    // 若旧树可解析（明文），把它整棵挂到新根之下，谱系得以保留；
    // 新根 ID 取旧树最大 ID + 1（block 11 的 node_index 引用这些 ID）。
    let new_id = graft
        .as_deref()
        .and_then(max_node_id)
        .map(|m| m + 1)
        .unwrap_or(1);
    let mut root = base_node(
        document_name,
        &new_rec,
        new_rec.len(),
        new_id,
        &operation,
        properties.clone(),
        label.clone(),
    );
    root.input_summary = Some(summary);
    root.body = state_body(&new_file.blocks);

    match &graft {
        Some(old_root) => {
            let parent = graft_parent(
                old_root,
                source_name,
                old_rec.len(),
                old_rec.flags,
                parent_resurrectable,
            );
            root.raw_children.push(parent);
        }
        None => {
            let mut parent = base_node(
                source_name,
                &old_rec,
                old_rec.len(),
                0,
                "invalid",
                properties,
                label,
            );
            parent.resurrectable = parent_resurrectable;
            parent.body = state_body(&src.blocks);
            root.children.push(parent);
        }
    }
    Ok(root)
}

/// 由编辑序列归纳 SnapGene 的 `operation` 与 `<InputSummary>` 属性。
///
/// 词表取自 13 个自带历史树的样例（见 spec §8）：单次区间替换 → `replace`，
/// 单次插入 → `insertFragment` + `manipulation="insert"`，其余（删除/多次编辑）
/// 归并为一段 `replace`。val1/val2 为**编辑前坐标系**的 0-based 区间。
fn summarize_edits(edits: &[Edit]) -> (String, InputSummary) {
    let span = edit_span(edits);
    match edits {
        [Edit::Delete { start, end }] => (
            "delete".into(), InputSummary::range("remove", *start, *end),
        ),
        [Edit::Insert { at, seq }] => (
            "insert".to_string(),
            InputSummary::range("insert", *at, at + seq.len()),
        ),
        _ => ("replace".to_string(), InputSummary::range("replace", span.0, span.1)),
    }
}

/// 编辑影响范围（编辑前坐标，0-based 半开）；插入视为零宽点。
fn edit_span(edits: &[Edit]) -> (usize, usize) {
    let (mut lo, mut hi) = (usize::MAX, 0usize);
    for e in edits {
        let (s, en) = match e {
            Edit::Delete { start, end } | Edit::Replace { start, end, .. } => (*start, *end),
            Edit::Insert { at, .. } => (*at, *at),
        };
        lo = lo.min(s);
        hi = hi.max(en);
    }
    if lo == usize::MAX {
        (0, 0)
    } else {
        (lo, hi)
    }
}

/// 取出 `<Features>…</Features>` 与 `<Primers>…</Primers>`（历史树内无 `<?xml?>` 前缀、
/// 无属性，实测形态见 `2. mCherry BsaI.dna`）。
fn state_body(blocks: &[crate::block::Block]) -> String {
    let mut out = String::new();
    for (t, tag) in [(10u8, "Features"), (5u8, "Primers")] {
        let Some(b) = blocks.iter().find(|b| b.block_type == t) else {
            continue;
        };
        let Ok(text) = std::str::from_utf8(&b.payload) else {
            continue;
        };
        let Some(start) = text.find(&format!("<{tag}")) else {
            continue;
        };
        let Some(open_end) = text[start..].find('>').map(|i| start + i + 1) else {
            continue;
        };
        let Some(close) = text.rfind(&format!("</{tag}>")) else {
            continue;
        };
        if close < open_end {
            continue;
        }
        out.push_str(&format!("<{tag}>"));
        out.push_str(&text[open_end..close]);
        out.push_str(&format!("</{tag}>"));
    }
    out
}

/// Extract the existing root, decoding XZ first and accepting self-closing leaves.
pub fn existing_root_xml(src: &DnaFile) -> Option<String> {
    let b = src.block(7)?;
    let text = decode_payload(&b.payload).ok()?;
    parse_entries(&text).ok()?;
    let t = text.trim_start();
    if !t.starts_with("<?xml") && !t.starts_with("<HistoryTree") {
        return None; // XZ 或未知编码
    }
    let (_, tag_end) = crate::annot::next_tag(&text, "Node", 0)?;
    let start = text[..tag_end].rfind("<Node")?;
    let end = text.rfind("</Node>").map(|i| i + "</Node>".len())
        .or_else(|| text[start..].find("/>").map(|i| start + i + 2))?;
    if end <= start {
        return None;
    }
    Some(text[start..end].to_string())
}

/// 树中出现过的最大 `ID` 属性值（`recentID` 这类不会误匹配：要求空格前置）。
pub fn max_node_id(xml: &str) -> Option<u32> {
    let mut best: Option<u32> = None;
    let mut pos = 0usize;
    while let Some(i) = xml[pos..].find(" ID=\"") {
        let at = pos + i + " ID=\"".len();
        let end = at + xml[at..].find('"')?;
        if let Ok(v) = xml[at..end].parse::<u32>() {
            best = Some(best.map_or(v, |b: u32| b.max(v)));
        }
        pos = end;
    }
    best
}

/// 移植旧树根节点：保留其内部内容与子节点，只把顶层标量属性更新为编辑前状态。
fn graft_parent(
    old_root: &str,
    name: &str,
    seq_len: usize,
    flags: MoleculeFlags,
    resurrectable: Option<u8>,
) -> String {
    let Some(gt) = old_root.find('>') else {
        return old_root.to_string();
    };
    let self_closing = old_root[..gt].trim_end().ends_with('/');
    let mut open = old_root[..gt].trim_end().trim_end_matches('/').to_string();
    open = patch_attr(&open, "name", name);
    open = patch_attr(&open, "seqLen", &seq_len.to_string());
    open = patch_attr(
        &open,
        "strandedness",
        match flags.strandedness() {
            Strandedness::Double => "double",
            Strandedness::Single => "single",
        },
    );
    open = patch_attr(&open, "circular", if flags.has(MoleculeFlags::CIRCULAR) { "1" } else { "0" });
    if let Some(r) = resurrectable {
        open = patch_attr(&open, "resurrectable", &r.to_string());
    }
    format!("{open}{}{}", if self_closing { "/" } else { "" }, &old_root[gt..])
}

/// 在开标签文本（不含结尾 `>`）里替换/插入 `key="value"`。
fn patch_attr(open_tag: &str, key: &str, value: &str) -> String {
    let needle = format!(" {key}=\"");
    if let Some(at) = open_tag.find(&needle) {
        let vstart = at + needle.len();
        if let Some(rel) = open_tag[vstart..].find('"') {
            let vend = vstart + rel;
            let mut s = String::with_capacity(open_tag.len() + value.len());
            s.push_str(&open_tag[..vstart]);
            s.push_str(&xml_escape(value));
            s.push_str(&open_tag[vend..]);
            return s;
        }
    }
    format!("{open_tag} {key}=\"{}\"", xml_escape(value))
}

/// 取 `<tag>文本</tag>` 的内容（首个匹配）。
fn tag_text(xml: &str, tag: &str) -> Option<String> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let s = xml.find(&open)? + open.len();
    let e = xml[s..].find(&close)? + s;
    Some(xml[s..e].to_string())
}

fn parse_properties(xml: &str) -> SeqProperties {
    let int = |tag: &str| -> i64 {
        tag_text(xml, tag)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(0)
    };
    SeqProperties {
        upstream_stickiness: int("UpstreamStickiness"),
        downstream_stickiness: int("DownstreamStickiness"),
        upstream_modification: tag_text(xml, "UpstreamModification")
            .unwrap_or_else(|| "Unmodified".into()),
        downstream_modification: tag_text(xml, "DownstreamModification")
            .unwrap_or_else(|| "Unmodified".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::write::{rewrite_after_edit, NewFile};

    #[test]
    fn parses_branches_summaries_and_escaped_names() {
        let xml = r#"<HistoryTree><Node name="A &amp; B" ID="3"><InputSummary val1="7"/><InputSummary name1="BamHI"/><Node name="first"/><Node name="second"><Node name="old"/></Node></Node></HistoryTree>"#;
        let n = parse_entries(xml).unwrap();
        assert_eq!(n.len(), 4);
        assert_eq!(n[0].label(), "A & B");
        assert_eq!(n[0].summaries.len(), 2);
        assert_eq!(n.iter().map(|n| n.parent).collect::<Vec<_>>(), vec![None, Some(0), Some(0), Some(2)]);
        assert!(parse_entries("<HistoryTree><Node name=\"x\">").is_err());
    }

    #[test]
    fn self_closing_history_is_grafted_not_replaced() {
        let mut src = demo_file();
        src.blocks.push(crate::Block::new(7, br#"<HistoryTree><Node name="ancestor" ID="17" seqLen="31" operation="invalid"/></HistoryTree>"#.to_vec()));
        let edits = [Edit::Insert { at: 2, seq: b"AA".to_vec() }];
        let new = rewrite_after_edit(&src, &edits).unwrap();
        let root = edit_history(&new, &src, &edits, "new", "old", None).unwrap();
        assert_eq!(root.id, 18);
        let text = decode_payload(&to_block_payload(&root)).unwrap();
        let nodes = parse_entries(&text).unwrap();
        assert_eq!(nodes.len(), 2);
        assert_eq!(nodes[1].attr("ID"), "17");
    }

    #[test]
    fn compressed_sample_history_survives_edit() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile/RNA files/Homo sapiens mitochondrial tRNA-Ala.rna");
        if !path.exists() { return; }
        let src = DnaFile::parse(&std::fs::read(path).unwrap()).unwrap();
        let xml = decode_payload(&src.block(7).unwrap().payload).unwrap();
        assert_eq!(parse_entries(&xml).unwrap().len(), 2);
        let edits = [Edit::Insert { at: 2, seq: b"AA".to_vec() }];
        let new = rewrite_after_edit(&src, &edits).unwrap();
        let root = edit_history(&new, &src, &edits, "new", "old", None).unwrap();
        assert_eq!(parse_entries(&decode_payload(&to_block_payload(&root)).unwrap()).unwrap().len(), 3);
    }

    fn demo_file() -> DnaFile {
        let bytes = NewFile::dna(b"ATGCGTACGTTAGCCTAGCTAGCTTAAAGGCC".to_vec())
            .with_name("demo")
            .to_bytes();
        DnaFile::parse(&bytes).unwrap()
    }

    #[test]
    fn edit_history_shape() {
        let src = demo_file();
        let edits = [Edit::Delete { start: 5, end: 10 }];
        let out = rewrite_after_edit(&src, &edits).unwrap();
        let node = edit_history(&out, &src, &edits, "out.dna", "src.dna", None).unwrap();
        let xml = String::from_utf8(to_block_payload(&node)).unwrap();
        assert!(xml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?><HistoryTree><Node"));
        assert!(xml.contains("name=\"out.dna\""), "{xml}");
        assert!(xml.contains("seqLen=\"27\""), "{xml}"); // 32 − 5
        assert!(xml.contains("operation=\"delete\""), "{xml}");
        assert!(xml.contains("<InputSummary manipulation=\"remove\" val1=\"5\" val2=\"10\"/>"), "{xml}");
        assert!(xml.contains("name=\"src.dna\" type=\"DNA\" seqLen=\"32\""), "{xml}");
        assert!(xml.contains("operation=\"invalid\""), "{xml}");
        assert!(xml.contains("customMapLabel=\"demo\""), "{xml}");
        assert!(xml.ends_with("</Node></Node></HistoryTree>"), "{xml}");
    }

    #[test]
    fn insert_uses_insert_operation() {
        let src = demo_file();
        let edits = [Edit::Insert {
            at: 3,
            seq: b"GGGG".to_vec(),
        }];
        let out = rewrite_after_edit(&src, &edits).unwrap();
        let node = edit_history(&out, &src, &edits, "o.dna", "s.dna", Some(1)).unwrap();
        let xml = String::from_utf8(to_block_payload(&node)).unwrap();
        assert!(xml.contains("operation=\"insert\""), "{xml}");
        assert!(xml.contains("<InputSummary manipulation=\"insert\" val1=\"3\" val2=\"7\"/>"), "{xml}");
        assert!(xml.contains("resurrectable=\"1\""), "{xml}");
    }

    #[test]
    fn chained_edits_preserve_lineage() {
        let src = demo_file();
        let e1 = [Edit::Delete { start: 5, end: 10 }];
        let f1 = rewrite_after_edit(&src, &e1).unwrap();
        let n1 = edit_history(&f1, &src, &e1, "v1.dna", "src.dna", None).unwrap();
        // 模拟“第一次编辑已保存”：结果文件带 block 7
        let mut f1_saved = f1.clone();
        f1_saved
            .blocks
            .push(crate::block::Block::new(7, to_block_payload(&n1)));

        let e2 = [Edit::Insert {
            at: 0,
            seq: b"TT".to_vec(),
        }];
        let f2 = rewrite_after_edit(&f1_saved, &e2).unwrap();
        let n2 = edit_history(&f2, &f1_saved, &e2, "v2.dna", "v1.dna", None).unwrap();
        let xml = String::from_utf8(to_block_payload(&n2)).unwrap();
        assert_eq!(xml.matches("<Node ").count(), 3, "{xml}");
        assert!(xml.contains("ID=\"2\""), "{xml}"); // 新根 ID = 旧树最大 ID + 1
        assert!(xml.contains("name=\"v2.dna\""), "{xml}");
        assert!(xml.contains("name=\"v1.dna\""), "{xml}");
        assert!(xml.contains("name=\"src.dna\""), "{xml}");
    }

    #[test]
    fn non_leaf_node_always_has_input_summary() {
        let node = Node {
            name: "x".into(),
            kind: MoleculeKind::Dna,
            seq_len: 1,
            strandedness: Strandedness::Double,
            id: 1,
            circular: false,
            resurrectable: None,
            operation: "replace".into(),
            properties: None,
            custom_map_label: None,
            input_summary: None,
            body: String::new(),
            raw_children: Vec::new(),
            children: vec![Node {
                name: "y".into(),
                kind: MoleculeKind::Dna,
                seq_len: 1,
                strandedness: Strandedness::Double,
                id: 0,
                circular: false,
                resurrectable: Some(1),
                operation: "invalid".into(),
                properties: None,
                custom_map_label: None,
                input_summary: None,
                body: String::new(),
                raw_children: Vec::new(),
                children: Vec::new(),
            }],
        };
        let xml = String::from_utf8(to_block_payload(&node)).unwrap();
        assert!(xml.contains("<InputSummary/><Node "), "{xml}");
    }
}
