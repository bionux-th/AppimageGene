//! 写入方向：新建文件与编辑写回（依据 `docs/SnapGene_dna_format_spec.md` §12）。
//!
//! 两条路径：
//!
//! 1. [`NewFile`] —— 从零构造一份 `.dna`（oracle 实验 T1 已证：真 SnapGene
//!    接受的最小集合是 block 9 + 0，这里按 SnapGene 习惯写 9/0/8/10/5/6）；
//! 2. [`rewrite_after_edit`] —— 序列编辑后的写回：重建序列块，丢弃派生缓存
//!    （block 3/7/11/34，T2/T3 已由 oracle 裁决），其余块原样保留，并按编辑
//!    平移/裁剪块 5（引物结合位点）与块 10（特征）中的区间。
//!
//! 设计原则与读入侧一致：**不解释的字节原样保留**，只动必须动的部分。

use crate::block::Block;
use crate::dna::DnaFile;
use crate::error::FormatError;
use crate::util::{html_wrap, snapgene_date, snapgene_utc_attr, utc_now, uuid_v4, xml_escape};
use hx_core::{MoleculeFlags, MoleculeKind, SeqRecord};

/// 序列变更后必须丢弃的派生块。
///
/// - 3：酶位点缓存（T2：真 SnapGene 原样信任，旧缓存会导致错误显示）；
/// - 7/11：历史树与节点快照（T3：旧树会与新序列自相矛盾；丢弃后显示
///   `No recorded history`）——P1 之后由历史写回模块重新生成；
/// - 34：RNA 二级结构（随序列变化失效）。
pub const DROPPED_ON_SEQUENCE_CHANGE: &[u8] = &[3, 7, 11, 34];

// ---------------------------------------------------------------------------
// 新建文件
// ---------------------------------------------------------------------------

/// 待写入的特征（`block 10` 的一条 `<Feature>`，单片段）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    /// 1-based 闭区间起点。
    pub start: usize,
    /// 1-based 闭区间终点。
    pub end: usize,
    pub name: String,
    /// SnapGene 特征类型（`CDS`/`promoter`/`misc_feature`…）。
    pub kind: String,
    /// 颜色 `#rrggbb`。
    pub color: String,
    /// 方向：1 = 正向，2 = 反向，0 = 无方向。
    pub directionality: i32,
    /// 除 `start..end` 之外的**附加片段**（1-based 闭区间）。
    ///
    /// `start`/`end` 仍是首个片段（新建特征 = 唯一片段）；导入 `join(...)` 位置
    /// 时把其余外显子放这里，写出为同名的多个 `<Segment>`。
    pub extra_segments: Vec<(usize, usize)>,
}

impl Feature {
    pub fn new(name: &str, start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            name: name.to_string(),
            kind: "misc_feature".into(),
            color: "#ff0000".into(),
            directionality: 1,
            extra_segments: Vec::new(),
        }
    }

    /// 属性顺序/取值对齐真 SnapGene 实测（6530 个特征）：`recentID name type
    /// directionality allowSegmentOverlaps consecutiveTranslationNumbering`；
    /// `directionality="0"` **省略属性**（全语料 0 次出现），`translated` 只在
    /// 片段上、且只在为真时写。
    fn to_xml(&self, index: usize) -> String {
        let dir = match self.directionality {
            1 => " directionality=\"1\"",
            2 => " directionality=\"2\"",
            _ => "",
        };
        let mut segs = format!(
            "<Segment range=\"{}-{}\" color=\"{}\" type=\"standard\"/>",
            self.start,
            self.end,
            xml_escape(&self.color)
        );
        for (a, b) in &self.extra_segments {
            segs.push_str(&format!(
                "<Segment range=\"{a}-{b}\" color=\"{}\" type=\"standard\"/>",
                xml_escape(&self.color)
            ));
        }
        format!(
            "<Feature recentID=\"{index}\" name=\"{name}\" type=\"{kind}\"{dir} \
             allowSegmentOverlaps=\"0\" consecutiveTranslationNumbering=\"1\">{segs}</Feature>",
            name = xml_escape(&self.name),
            kind = xml_escape(&self.kind),
        )
    }
}

/// 新建一份文件的描述。
#[derive(Debug, Clone)]
pub struct NewFile {
    pub kind: MoleculeKind,
    pub flags: MoleculeFlags,
    pub sequence: Vec<u8>,
    /// 图谱标签（`CustomMapLabel`）；`None` 时不写该元素。
    pub name: Option<String>,
    /// 描述（纯文本，会包装成 SnapGene 富文本 HTML）。
    pub description: Option<String>,
    /// `CreatedBy` 字段。
    pub created_by: Option<String>,
    pub features: Vec<Feature>,
    /// 附带的原始块（如复制来的 block 13/14），按给定顺序追加在标准块之后。
    pub extra_blocks: Vec<Block>,
}

impl NewFile {
    /// 新建 DNA 文件：线性、双链、无甲基化（`0x02`）。
    pub fn dna(sequence: Vec<u8>) -> Self {
        Self::of(MoleculeKind::Dna, sequence)
    }

    /// 按分子类型新建文件。
    ///
    /// 缺省链性/拓扑与真 SnapGene 的 New File 面板一致：DNA = 双链线性
    /// （`0x02`）、RNA = **单链**线性（`0x00`）、蛋白无链性/拓扑（`0x00`）。
    pub fn of(kind: MoleculeKind, sequence: Vec<u8>) -> Self {
        let flags = match kind {
            MoleculeKind::Dna => MoleculeFlags::from_byte(MoleculeFlags::DOUBLE_STRANDED),
            MoleculeKind::Rna | MoleculeKind::Protein => MoleculeFlags::from_byte(0),
        };
        Self {
            kind,
            flags,
            sequence,
            name: None,
            description: None,
            created_by: None,
            features: Vec::new(),
            extra_blocks: Vec::new(),
        }
    }

    pub fn with_name(mut self, name: &str) -> Self {
        self.name = Some(name.to_string());
        self
    }

    pub fn with_description(mut self, desc: &str) -> Self {
        self.description = Some(desc.to_string());
        self
    }

    pub fn with_flags(mut self, flags: MoleculeFlags) -> Self {
        self.flags = flags;
        self
    }

    pub fn with_features(mut self, features: Vec<Feature>) -> Self {
        self.features = features;
        self
    }

    /// 块顺序取 203 样例中最常见者（79 份）：`9, 0, 3, 8, 10, 5, 6, 13, 14`，
    /// 去掉本实现暂不生成的 3/13/14。
    pub fn to_bytes(&self) -> Vec<u8> {
        let now = utc_now();
        let mut blocks: Vec<Block> = Vec::new();
        blocks.push(Block::new(
            9,
            crate::dna::Header::new(self.kind.to_type_of_sequence()).raw,
        ));
        let rec = SeqRecord::new(self.kind, self.flags, self.sequence.clone());
        blocks.push(Block::new(
            self.kind.sequence_block_type(),
            rec.to_block_payload(),
        ));
        blocks.push(Block::new(8, additional_properties_xml().as_bytes().to_vec()));
        blocks.push(Block::new(10, self.features_xml().into_bytes()));
        blocks.push(Block::new(5, empty_primers_xml().as_bytes().to_vec()));
        blocks.push(Block::new(6, self.notes_xml(&now).into_bytes()));
        blocks.extend(self.extra_blocks.iter().cloned());
        let mut out = Vec::new();
        crate::block::write_blocks(&mut out, &blocks);
        out
    }

    fn features_xml(&self) -> String {
        let mut s = format!(
            "<?xml version=\"1.0\"?><Features nextValidID=\"{}\">",
            self.features.len()
        );
        for (i, f) in self.features.iter().enumerate() {
            s.push_str(&f.to_xml(i));
        }
        s.push_str("</Features>");
        s
    }

    fn notes_xml(&self, now: &(i64, u32, u32, u32, u32, u32)) -> String {
        let date = snapgene_date(*now);
        let utc = snapgene_utc_attr(*now);
        let mut s = String::from("<Notes>");
        s.push_str(&format!("<UUID>{}</UUID>", uuid_v4()));
        s.push_str("<Type>Synthetic</Type>");
        s.push_str("<ConfirmedExperimentally>0</ConfirmedExperimentally>");
        if let Some(name) = &self.name {
            s.push_str(&format!(
                "<CustomMapLabel>{}</CustomMapLabel><UseCustomMapLabel>1</UseCustomMapLabel>",
                xml_escape(name)
            ));
        }
        if let Some(desc) = &self.description {
            s.push_str(&format!("<Description>{}</Description>", html_wrap(desc)));
        }
        s.push_str(&format!("<Created UTC=\"{utc}\">{date}</Created>"));
        s.push_str(&format!("<LastModified UTC=\"{utc}\">{date}</LastModified>"));
        if let Some(by) = &self.created_by {
            s.push_str(&format!("<CreatedBy>{}</CreatedBy>", xml_escape(by)));
        }
        s.push_str("<SequenceClass>UNA</SequenceClass>");
        s.push_str("<TransformedInto>Unspecified</TransformedInto>");
        s.push_str("</Notes>");
        s
    }
}

/// 空引物表模板（含 `HybridizationParams` 默认值，取自样例实测）。
pub fn empty_primers_xml() -> &'static str {
    "<?xml version=\"1.0\"?><Primers nextValidID=\"0\"><HybridizationParams \
     minContinuousMatchLen=\"10\" allowMismatch=\"1\" minMeltingTemperature=\"40\" \
     showAdditionalFivePrimeMatches=\"1\" minimumFivePrimeAnnealing=\"15\"/></Primers>"
}

/// Serialize mapped annotations without reducing them to the minimal `Feature`
/// creation model. Only fields represented by `FeatureAnnotation` are written;
/// this is not a lossless serializer for arbitrary source feature XML.
pub fn feature_annotations_xml(features:&[crate::annot::FeatureAnnotation],len:usize)->Vec<u8> {
    let esc=|s:&str|xml_escape(s).replace('\n',"&#10;").replace('\r',"&#13;");
    let mut body=String::new();
    let mut count=0;
    for f in features {
        let segments=f.segments.iter().filter(|s|s.range.start<len && s.range.end<=len && s.range.start!=s.range.end).collect::<Vec<_>>();
        if segments.is_empty() {continue;}
        body.push_str(&format!("<Feature recentID=\"{count}\" name=\"{}\" type=\"{}\"",esc(&f.name),esc(&f.kind)));
        if !f.sub_type.is_empty() {body.push_str(&format!(" subType=\"{}\"",esc(&f.sub_type)));}
        if !f.visible {body.push_str(" visible=\"0\"");}
        if matches!(f.directionality,1..=3) {body.push_str(&format!(" directionality=\"{}\"",f.directionality));}
        if matches!(f.reading_frame,-3..=-1 | 1..=3) {body.push_str(&format!(" readingFrame=\"{}\"",f.reading_frame));}
        if f.translate_first_codon_as_met {body.push_str(" translateFirstCodonAsMet=\"1\"");}
        if let Some(end)=f.max_run_on.filter(|end|*end>0 && *end<=len) {body.push_str(&format!(" maxRunOn=\"{end}\""));}
        body.push_str(" allowSegmentOverlaps=\"0\" consecutiveTranslationNumbering=\"1\">");
        for s in segments {
            body.push_str(&format!("<Segment range=\"{}-{}\" color=\"{}\" type=\"{}\"{}/>",
                s.range.start+1,s.range.end,esc(&s.color),esc(&s.seg_type),if s.translated {" translated=\"1\""}else{""}));
        }
        if !f.note.is_empty() {body.push_str(&format!("<Q name=\"note\"><V text=\"{}\"/></Q>",esc(&f.note)));}
        body.push_str("</Feature>");
        count+=1;
    }
    format!("<?xml version=\"1.0\"?><Features nextValidID=\"{count}\">{body}</Features>").into_bytes()
}

/// 附加序列属性（`block 8`）：未修饰末端（质粒/新文件的常规值）。
pub fn additional_properties_xml() -> &'static str {
    "<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness>\
     <DownstreamStickiness>0</DownstreamStickiness>\
     <UpstreamModification>Unmodified</UpstreamModification>\
     <DownstreamModification>Unmodified</DownstreamModification>\
     </AdditionalSequenceProperties>"
}

// ---------------------------------------------------------------------------
// 编辑
// ---------------------------------------------------------------------------

/// 一次序列编辑（坐标一律 **0-based 半开区间**，`Insert.at` 为插入点）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Edit {
    Delete { start: usize, end: usize },
    Insert { at: usize, seq: Vec<u8> },
    Replace { start: usize, end: usize, seq: Vec<u8> },
}

impl Edit {
    /// 该编辑“删除”的碱基数。
    fn removed(&self) -> usize {
        match self {
            Edit::Delete { start, end } => end - start,
            Edit::Replace { start, end, .. } => end - start,
            Edit::Insert { .. } => 0,
        }
    }

    /// 该编辑“插入”的碱基数。
    fn added(&self) -> usize {
        match self {
            Edit::Insert { seq, .. } | Edit::Replace { seq, .. } => seq.len(),
            Edit::Delete { .. } => 0,
        }
    }

    /// 在序列上施加编辑（坐标以当前序列为准）。
    pub fn apply(&self, seq: &[u8]) -> Result<Vec<u8>, FormatError> {
        let bad = |msg: &str| FormatError::BadEdit(msg.to_string());
        let check = |start: usize, end: usize| -> Result<(), FormatError> {
            if start > end || end > seq.len() {
                Err(bad(&format!(
                    "编辑区间 {start}..{end} 越界（序列长 {}）",
                    seq.len()
                )))
            } else {
                Ok(())
            }
        };
        let mut out = Vec::with_capacity(seq.len());
        match self {
            Edit::Delete { start, end } => {
                check(*start, *end)?;
                out.extend_from_slice(&seq[..*start]);
                out.extend_from_slice(&seq[*end..]);
            }
            Edit::Insert { at, seq: ins } => {
                if *at > seq.len() {
                    return Err(bad(&format!("插入点 {at} 越界（序列长 {}）", seq.len())));
                }
                out.extend_from_slice(&seq[..*at]);
                out.extend_from_slice(ins);
                out.extend_from_slice(&seq[*at..]);
            }
            Edit::Replace { start, end, seq: ins } => {
                check(*start, *end)?;
                out.extend_from_slice(&seq[..*start]);
                out.extend_from_slice(ins);
                out.extend_from_slice(&seq[*end..]);
            }
        }
        Ok(out)
    }
}

/// 在序列上顺序施加多次编辑。
pub fn apply_edits(seq: &[u8], edits: &[Edit]) -> Result<Vec<u8>, FormatError> {
    let mut cur = seq.to_vec();
    for e in edits {
        cur = e.apply(&cur)?;
    }
    Ok(cur)
}

/// 编辑后的整份文件写回（不落盘，返回新的 [`DnaFile`]）。
///
/// 契约（规格 §12.2，T2/T3 已由 oracle 裁决）：
/// - 序列块按新序列重建（保留 flags）；
/// - 丢弃 [`DROPPED_ON_SEQUENCE_CHANGE`]；
/// - 块 5/10 的区间按编辑平移/裁剪，完全落入删除区的元素被移除；
/// - 块 6 的 `<LastModified>` 更新为当前时间；
/// - 其余块（13/14/28/16/17/20/2/未知）原样保留、顺序不变。
pub fn rewrite_after_edit(file: &DnaFile, edits: &[Edit]) -> Result<DnaFile, FormatError> {
    let kind = file.header.molecule_kind().ok_or_else(|| {
        FormatError::BadHeader(format!("未知分子类型 {}", file.header.type_of_sequence))
    })?;
    let seq_block_type = kind.sequence_block_type();
    let old = file
        .block(seq_block_type)
        .ok_or_else(|| FormatError::BadHeader("缺少序列块".into()))?;
    let old_rec = SeqRecord::from_block_payload(kind, &old.payload)
        .ok_or_else(|| FormatError::BadHeader("序列块载荷为空".into()))?;
    let new_seq = apply_edits(&old_rec.sequence, edits)?;
    let changed = new_seq != old_rec.sequence;

    let old_len = old_rec.sequence.len();
    let mut blocks: Vec<Block> = Vec::with_capacity(file.blocks.len());
    for b in &file.blocks {
        if b.block_type == seq_block_type {
            let rec = SeqRecord::new(kind, old_rec.flags, new_seq.clone());
            blocks.push(Block::new(seq_block_type, rec.to_block_payload()));
            continue;
        }
        if changed && DROPPED_ON_SEQUENCE_CHANGE.contains(&b.block_type) {
            continue;
        }
        if changed && b.block_type == 5 {
            blocks.push(Block::new(5, adjust_primer_xml(&b.payload, edits, old_len)));
            continue;
        }
        if changed && b.block_type == 10 {
            blocks.push(Block::new(10, adjust_feature_xml(&b.payload, edits, old_len)));
            continue;
        }
        if b.block_type == 6 {
            blocks.push(Block::new(6, touch_last_modified(&b.payload)));
            continue;
        }
        blocks.push(b.clone());
    }
    Ok(DnaFile {
        header: file.header.clone(),
        blocks,
    })
}

/// 在已重写的文件里插入 **block 7 历史树**（根 = 编辑后状态，子 = 编辑前状态），
/// 位置紧跟 block 6（实测块序 `…, 6, 7, 17, …`）。
///
/// - `new_file`：[`rewrite_after_edit`] 的产物（就地插入）；
/// - `src`：编辑前的文件（提供父节点 body 与旧树谱系）；
/// - `document_name` / `source_name`：根/父节点显示名（通常是保存后与打开时的文件名）；
/// - `parent_resurrectable`：父节点 `resurrectable` 标记（无 block 11 快照时传 `None`，
///   T10 实测该标记不影响显示）。
///
/// 返回树中节点数。`hx-cli edit` 与 UI 保存（`hx-algo::EditSession`）共用此实现，
/// 显示效果已由 oracle T10 裁决（`docs/oracle_T1_T4_results.zh-CN.md` §6）。
pub fn insert_history_block(
    new_file: &mut DnaFile,
    src: &DnaFile,
    edits: &[Edit],
    document_name: &str,
    source_name: &str,
    parent_resurrectable: Option<u8>,
) -> Result<usize, FormatError> {
    let node = crate::history::edit_history(
        new_file,
        src,
        edits,
        document_name,
        source_name,
        parent_resurrectable,
    )?;
    let payload = crate::history::to_block_payload(&node);
    let nodes = payload.windows(6).filter(|w| *w == b"<Node ").count();
    new_file.blocks.retain(|b| b.block_type != 7);
    let pos = new_file
        .blocks
        .iter()
        .position(|b| b.block_type == 6)
        .map(|i| i + 1)
        .unwrap_or(new_file.blocks.len());
    new_file.blocks.insert(pos, Block::new(7, payload));
    Ok(nodes)
}

/// 更新 `block 6` 的 `<LastModified>`（无则插入在 `</Notes>` 前）。
pub fn touch_last_modified(payload: &[u8]) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(payload) else {
        return payload.to_vec();
    };
    let now = utc_now();
    let elem = format!(
        "<LastModified UTC=\"{}\">{}</LastModified>",
        snapgene_utc_attr(now),
        snapgene_date(now)
    );
    if let Some(start) = text.find("<LastModified") {
        if let Some(rel_end) = text[start..].find("</LastModified>") {
            let end = start + rel_end + "</LastModified>".len();
            let mut out = String::with_capacity(text.len());
            out.push_str(&text[..start]);
            out.push_str(&elem);
            out.push_str(&text[end..]);
            return out.into_bytes();
        }
    }
    match text.find("</Notes>") {
        Some(pos) => {
            let mut out = String::with_capacity(text.len() + elem.len());
            out.push_str(&text[..pos]);
            out.push_str(&elem);
            out.push_str(&text[pos..]);
            out.into_bytes()
        }
        None => payload.to_vec(),
    }
}

// ---------------------------------------------------------------------------
// 区间变换（特征 / 引物结合位点）
// ---------------------------------------------------------------------------

/// 0-based 半开区间。
type Iv = (usize, usize);

/// 对一个 0-based 半开区间施加编辑；返回 `None` 表示区间被完全删除。
fn apply_edit_to_iv(iv: Iv, edit: &Edit) -> Option<Iv> {
    let (s, e) = iv;
    match *edit {
        Edit::Delete { start: ds, end: de } => {
            if de <= s {
                // 删除区在区间左侧 → 整体左移
                let d = de - ds;
                Some((s - d, e - d))
            } else if ds >= e {
                // 删除区在区间右侧 → 不动
                Some(iv)
            } else {
                // 有重叠：按重叠量裁剪两端
                let cut_before = s.min(de).saturating_sub(ds);
                let cut_upto = e.min(de).saturating_sub(ds);
                let (s2, e2) = (s - cut_before, e - cut_upto);
                if e2 <= s2 {
                    None
                } else {
                    Some((s2, e2))
                }
            }
        }
        Edit::Insert { at, ref seq } => {
            let n = seq.len();
            if s >= at {
                Some((s + n, e + n))
            } else if e > at {
                Some((s, e + n))
            } else {
                Some(iv)
            }
        }
        Edit::Replace { start, end, ref seq } => {
            // 先删后插：删除 [start,end) 后，插入点就是 `start`（新坐标系）。
            let (s0, _) = iv;
            let (s, e) = apply_edit_to_iv(iv, &Edit::Delete { start, end })?;
            let n = seq.len();
            if s0 <= start && e > start {
                // 区间包住了替换起点（它不是被插到后面去的）→ 起点不动，只补上插入长度
                Some((s, e + n))
            } else if s >= start {
                Some((s + n, e + n))
            } else {
                Some((s, e))
            }
        }
    }
}

/// 区间片段：`tail` 表示终点跟随序列末端（环状跨原点区间的尾段）。
#[derive(Debug, Clone, Copy)]
struct Piece {
    s: usize,
    e: usize,
    tail: bool,
}

/// 区间坐标的基准。
///
/// SnapGene 的两种写法混用，实测（pJET1.2.dna / pET28b-SprB-*.dna）：
/// 特征 `<Segment range="A-B">` 是 **1-based 闭区间**；
/// 引物 `<BindingSite location>` 与 `<Component hybridizedRange>` 是
/// **0-based 闭区间**（0 是合法起点，不能减一）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Base {
    One,
    Zero,
}

impl Base {
    /// 闭区间起点 → 0-based 半开起点。
    fn start(self, a: usize) -> usize {
        match self {
            Base::One => a.saturating_sub(1),
            Base::Zero => a,
        }
    }

    /// 闭区间终点 → 0-based 半开终点。
    fn end(self, b: usize) -> usize {
        match self {
            Base::One => b,
            Base::Zero => b + 1,
        }
    }

    /// 0-based 半开片段 → 该基准下的闭区间文本。
    fn text(self, s: usize, e: usize) -> String {
        match self {
            Base::One => format!("{}-{}", s + 1, e),
            Base::Zero => format!("{}-{}", s, e.saturating_sub(1)),
        }
    }
}

/// 单个 `<Segment>`/`<BindingSite>` 变换后的片段集合。
///
/// 输入 `(a, b)` 为闭区间（基准见 [`Base`]）；`a > b` 表示环状跨原点
/// （覆盖 `[start, len) ∪ [0, end)`）。输出为若干 0-based 半开片段。
fn transform_range(a: usize, b: usize, base: Base, edits: &[Edit], old_len: usize) -> Vec<Piece> {
    let mut len = old_len;
    let mut parts: Vec<Piece> = if a <= b {
        vec![Piece {
            s: base.start(a),
            e: base.end(b),
            tail: false,
        }]
    } else {
        vec![
            Piece {
                s: base.start(a),
                e: len,
                tail: true,
            },
            Piece {
                s: 0,
                e: base.end(b),
                tail: false,
            },
        ]
    };
    for ed in edits {
        let mut next = Vec::with_capacity(parts.len());
        for p in parts {
            if let Some((s2, e2)) = apply_edit_to_iv((p.s, p.e), ed) {
                next.push(Piece {
                    s: s2,
                    e: e2,
                    tail: p.tail,
                });
            }
        }
        parts = next;
        len = (len + ed.added()).saturating_sub(ed.removed());
        for p in &mut parts {
            if p.tail {
                p.e = len; // 尾段终点＝序列末端
            }
        }
        parts.retain(|p| p.e > p.s);
        if parts.is_empty() {
            break;
        }
    }
    parts
}

/// 片段集合 → 一个或多个闭区间文本（基准见 [`Base`]）。
///
/// 若恰好构成跨原点形态（首段始于 0、尾段止于序列末端），合并为 `A-B`（`A > B`）。
fn ranges_to_text(parts: &[Piece], len: usize, base: Base) -> Vec<String> {
    if parts.len() == 2 {
        for (i, t) in parts.iter().enumerate() {
            // 尾段止于序列末端 + 首段始于 0 → 合并为跨原点形式（A > B）
            if t.tail && t.e == len && parts[1 - i].s == 0 {
                return vec![base.text(t.s, parts[1 - i].e)];
            }
        }
    }
    parts.iter().map(|p| base.text(p.s, p.e)).collect()
}

/// 在 XML 文本里查 `attr="A-B"` 属性值（仅匹配纯数字区间）。
fn find_range_attr(text: &str, attr: &str) -> Option<(usize, usize, usize, usize)> {
    // 返回 (属性值起点, 属性值终点, A, B)
    let key = format!("{attr}=\"");
    let pos = text.find(&key)?;
    let vstart = pos + key.len();
    let vend = text[vstart..].find('"')? + vstart;
    let value = &text[vstart..vend];
    let (a, b) = value.split_once('-')?;
    let a: usize = a.parse().ok()?;
    let b: usize = b.parse().ok()?;
    Some((vstart, vend, a, b))
}

/// 块 10（特征 XML）：按编辑变换 `<Segment range="A-B">`。
///
/// 特征内所有片段都被删除时，整个 `<Feature>…</Feature>` 移除；部分保留时
/// 只重写幸存片段的 `range`。`<Feature …/>` 自闭合形式原样保留。
pub fn adjust_feature_xml(payload: &[u8], edits: &[Edit], old_len: usize) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(payload) else {
        return payload.to_vec();
    };
    let mut out = String::with_capacity(text.len());
    let mut pos = 0usize;
    while let Some(rel) = text[pos..].find("<Feature ") {
        let fstart = pos + rel;
        out.push_str(&text[pos..fstart]);
        // 找开标签结束
        let Some(gt) = text[fstart..].find('>') else {
            out.push_str(&text[fstart..]);
            return out.into_bytes();
        };
        let open_end = fstart + gt + 1;
        let self_closing = text[fstart..open_end].ends_with("/>");
        if self_closing {
            out.push_str(&text[fstart..open_end]);
            pos = open_end;
            continue;
        }
        let Some(rel_close) = text[open_end..].find("</Feature>") else {
            out.push_str(&text[fstart..]);
            return out.into_bytes();
        };
        let close_start = open_end + rel_close;
        let fend = close_start + "</Feature>".len();
        let inner = &text[open_end..close_start];

        // 逐个 <Segment …/> 处理
        let new_len = seq_len_after(edits, old_len);
        let mut new_inner = String::with_capacity(inner.len());
        let mut seg_pos = 0usize;
        let mut kept = 0usize;
        while let Some(rel_seg) = inner[seg_pos..].find("<Segment ") {
            let sstart = seg_pos + rel_seg;
            new_inner.push_str(&inner[seg_pos..sstart]);
            let Some(sgt) = inner[sstart..].find("/>") else {
                new_inner.push_str(&inner[sstart..]);
                seg_pos = inner.len();
                break;
            };
            let send = sstart + sgt + 2;
            let tag = &inner[sstart..send];
            if let Some((vs, ve, a, b)) = find_range_attr(tag, "range") {
                let parts = transform_range(a, b, Base::One, edits, old_len);
                let texts = ranges_to_text(&parts, new_len, Base::One);
                for (i, t) in texts.iter().enumerate() {
                    if i > 0 {
                        new_inner.push(' ');
                    }
                    new_inner.push_str(&tag[..vs]);
                    new_inner.push_str(t);
                    new_inner.push_str(&tag[ve..]);
                }
                if !texts.is_empty() {
                    kept += 1;
                }
            } else {
                new_inner.push_str(tag);
                kept += 1;
            }
            seg_pos = send;
        }
        new_inner.push_str(&inner[seg_pos..]);

        if kept > 0 || !inner.contains("<Segment ") {
            out.push_str(&text[fstart..open_end]);
            out.push_str(&new_inner);
            out.push_str("</Feature>");
        }
        pos = fend;
    }
    out.push_str(&text[pos..]);
    out.into_bytes()
}

/// 块 5（引物 XML）：变换 `<BindingSite location="A-B">` 与
/// `<Component hybridizedRange="A-B">`。
///
/// 坐标是 **0-based 闭区间**（见 [`Base`]），`0` 是合法起点。
/// 某个结合位点整体落在删除区里时只丢这**一个**位点；`<Primer>` 本身与
/// 它的其它位点、全部属性都保留（早期实现会把整个引物一起丢掉，
/// 于是“删一个碱基 → 第一个引物消失”）。
pub fn adjust_primer_xml(payload: &[u8], edits: &[Edit], old_len: usize) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(payload) else {
        return payload.to_vec();
    };
    let new_len = seq_len_after(edits, old_len);
    let mut out = String::with_capacity(text.len());
    let mut pos = 0usize;
    while let Some(rel) = text[pos..].find("<Primer ") {
        let pstart = pos + rel;
        out.push_str(&text[pos..pstart]);
        let Some(gt) = text[pstart..].find('>') else {
            out.push_str(&text[pstart..]);
            return out.into_bytes();
        };
        let open_end = pstart + gt + 1;
        if text[pstart..open_end].ends_with("/>") {
            // 没有结合位点的引物（自闭合）→ 原样保留
            out.push_str(&text[pstart..open_end]);
            pos = open_end;
            continue;
        }
        let Some(rel_close) = text[open_end..].find("</Primer>") else {
            out.push_str(&text[pstart..]);
            return out.into_bytes();
        };
        let pend = open_end + rel_close + "</Primer>".len();
        out.push_str(&adjust_binding_sites(
            &text[pstart..pend],
            edits,
            old_len,
            new_len,
        ));
        pos = pend;
    }
    out.push_str(&text[pos..]);
    out.into_bytes()
}

/// 逐 `<BindingSite>` 变换一个 `<Primer>…</Primer>` 元素的内部。
fn adjust_binding_sites(primer: &str, edits: &[Edit], old_len: usize, new_len: usize) -> String {
    let mut out = String::with_capacity(primer.len());
    let mut pos = 0usize;
    while let Some(rel) = primer[pos..].find("<BindingSite") {
        let bstart = pos + rel;
        out.push_str(&primer[pos..bstart]);
        let Some(gt) = primer[bstart..].find('>') else {
            out.push_str(&primer[bstart..]);
            return out;
        };
        let open_end = bstart + gt + 1;
        let bend = if primer[bstart..open_end].ends_with("/>") {
            open_end
        } else if let Some(rc) = primer[open_end..].find("</BindingSite>") {
            open_end + rc + "</BindingSite>".len()
        } else {
            out.push_str(&primer[bstart..]);
            return out;
        };
        let site = &primer[bstart..bend];
        match rewrite_ranges_in(site, edits, old_len, new_len, &["location", "hybridizedRange"]) {
            Some(new_site) => out.push_str(&new_site),
            None => {} // 该结合位点已不在序列上 → 只丢这一个位点
        }
        pos = bend;
    }
    out.push_str(&primer[pos..]);
    out
}

/// 重写一段 XML 里全部 `attr="A-B"` 属性；某区间整体落在删除区里时返回 `None`。
///
/// 非区间值（解析不了）原样保留，不当成“被删除”。
fn rewrite_ranges_in(
    text: &str,
    edits: &[Edit],
    old_len: usize,
    new_len: usize,
    attrs: &[&str],
) -> Option<String> {
    let mut out = String::with_capacity(text.len());
    let mut pos = 0usize;
    loop {
        // 找到最近的、属于 attrs 的属性名
        let mut best: Option<(usize, &str)> = None;
        for a in attrs {
            let key = format!("{a}=\"");
            if let Some(i) = text[pos..].find(&key) {
                let at = pos + i;
                if best.map_or(true, |(b, _)| at < b) {
                    best = Some((at, a));
                }
            }
        }
        let Some((at, attr)) = best else {
            out.push_str(&text[pos..]);
            return Some(out);
        };
        let vstart = at + attr.len() + 2;
        let Some(q) = text[vstart..].find('"') else {
            out.push_str(&text[pos..]);
            return Some(out);
        };
        let vend = vstart + q;
        let parsed = text[vstart..vend]
            .split_once('-')
            .and_then(|(a, b)| Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?)));
        let Some((a, b)) = parsed else {
            // 不是区间 → 原样保留
            out.push_str(&text[pos..vend + 1]);
            pos = vend + 1;
            continue;
        };
        let parts = transform_range(a, b, Base::Zero, edits, old_len);
        let texts = ranges_to_text(&parts, new_len, Base::Zero);
        let first = texts.first()?; // 空 = 该区间落在删除区 → 调用方丢弃这个元素
        out.push_str(&text[pos..vstart]);
        out.push_str(first);
        out.push('"');
        pos = vend + 1;
    }
}

/// 施加编辑后的序列长度。
fn seq_len_after(edits: &[Edit], old_len: usize) -> usize {
    let mut len = old_len;
    for e in edits {
        len = (len + e.added()).saturating_sub(e.removed());
    }
    len
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::DnaFile;

    fn ps(v: &[Piece]) -> Vec<(usize, usize)> {
        v.iter().map(|p| (p.s, p.e)).collect()
    }

    #[test]
    fn delete_shifts_and_clips() {
        // 删除 [2,4)（0-based）
        let e = Edit::Delete { start: 2, end: 4 };
        assert_eq!(apply_edit_to_iv((0, 2), &e), Some((0, 2))); // 之前不动
        assert_eq!(apply_edit_to_iv((4, 6), &e), Some((2, 4))); // 之后平移
        assert_eq!(apply_edit_to_iv((1, 5), &e), Some((1, 3))); // 跨越 → 缩短
        assert_eq!(apply_edit_to_iv((2, 4), &e), None); // 完全删除
    }

    #[test]
    fn insert_shifts_and_extends() {
        let e = Edit::Insert {
            at: 3,
            seq: b"AA".to_vec(),
        };
        assert_eq!(apply_edit_to_iv((3, 5), &e), Some((5, 7))); // 起点在插入点 → 平移
        assert_eq!(apply_edit_to_iv((5, 6), &e), Some((7, 8))); // 之后平移
        assert_eq!(apply_edit_to_iv((1, 3), &e), Some((1, 3))); // 前沿不扩
        assert_eq!(apply_edit_to_iv((1, 5), &e), Some((1, 7))); // 内部插入 → 延长
    }

    #[test]
    fn replace_is_delete_then_insert() {
        let e = Edit::Replace {
            start: 2,
            end: 5,
            seq: b"TTTT".to_vec(),
        };
        assert_eq!(apply_edit_to_iv((0, 2), &e), Some((0, 2)));
        assert_eq!(apply_edit_to_iv((5, 8), &e), Some((6, 9)));
        assert_eq!(apply_edit_to_iv((1, 6), &e), Some((1, 7))); // 跨越 → 净 +1
        assert_eq!(apply_edit_to_iv((2, 5), &e), None); // 完全被替换掉
    }

    #[test]
    fn wrap_range_tracks_sequence_end() {
        // 100 bp 环状：range=90-10（1-based），删除末端 5 bp → 仍是跨原点形式
        let e = Edit::Delete { start: 95, end: 100 };
        let parts = transform_range(90, 10, Base::One, &[e], 100);
        assert_eq!(ps(&parts), vec![(89, 95), (0, 10)]);
        let texts = ranges_to_text(&parts, 95, Base::One);
        assert_eq!(texts, vec!["90-10".to_string()]);
    }

    #[test]
    fn wrap_range_shifts_with_edit_before_origin() {
        // 删除 0-based 4..6：尾段前移、首段缩短，仍是合法 wrap（88-8）
        let e = Edit::Delete { start: 4, end: 6 };
        let parts = transform_range(90, 10, Base::One, &[e], 100);
        assert_eq!(ps(&parts), vec![(87, 98), (0, 8)]);
        let texts = ranges_to_text(&parts, 98, Base::One);
        assert_eq!(texts, vec!["88-8".to_string()]);
    }

    #[test]
    fn wrap_range_degrades_to_plain() {
        // 首段被整体删除 → 退化为普通区间
        let parts = transform_range(90, 10, Base::One, &[Edit::Delete { start: 0, end: 10 }], 100);
        assert_eq!(ps(&parts), vec![(79, 90)]);
        assert_eq!(ranges_to_text(&parts, 90, Base::One), vec!["80-90".to_string()]);
    }

    #[test]
    fn feature_xml_adjusts_and_drops() {
        let xml = b"<?xml version=\"1.0\"?><Features nextValidID=\"2\">\
            <Feature recentID=\"0\" name=\"A\" type=\"CDS\"><Segment range=\"10-20\" color=\"#ff0000\" type=\"standard\"/></Feature>\
            <Feature recentID=\"1\" name=\"B\" type=\"misc_feature\"><Segment range=\"30-40\" color=\"#00ff00\" type=\"standard\"/></Feature>\
            </Features>";
        // 删除 0-based 15..35：A(10-20) 缩短为 10-15；B(30-40) 去掉被删段后变为 16-20
        let out = adjust_feature_xml(xml, &[Edit::Delete { start: 15, end: 35 }], 60);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("range=\"10-15\""), "{s}");
        assert!(s.contains("range=\"16-20\""), "{s}");

        // 整个特征落在删除区内 → 该 <Feature> 被移除
        let out2 = adjust_feature_xml(xml, &[Edit::Delete { start: 9, end: 20 }], 60);
        let s2 = String::from_utf8(out2).unwrap();
        assert!(!s2.contains("name=\"A\""), "{s2}");
        assert!(s2.contains("name=\"B\""), "{s2}");
    }

    #[test]
    fn primer_xml_adjusts_and_drops() {
        let xml = b"<?xml version=\"1.0\"?><Primers nextValidID=\"2\">\
            <Primer recentID=\"0\" name=\"P1\" sequence=\"AAA\"><BindingSite location=\"5-7\" boundStrand=\"1\">\
            <Component hybridizedRange=\"5-7\" bases=\"AAA\"/></BindingSite><BindingSite simplified=\"1\" location=\"5-7\" boundStrand=\"1\"/></Primer>\
            <Primer recentID=\"1\" name=\"P2\" sequence=\"CCC\"><BindingSite location=\"20-22\" boundStrand=\"1\"/></Primer>\
            </Primers>";
        let out = adjust_primer_xml(xml, &[Edit::Insert { at: 0, seq: b"GG".to_vec() }], 30);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("location=\"7-9\""), "{s}");
        assert!(s.contains("hybridizedRange=\"7-9\""), "{s}");
        assert!(s.contains("location=\"22-24\""), "{s}");
    }

    /// 引物 `location` 是 0-based 闭区间：`0` 是合法起点，不能下溢。
    /// （真实文件里 `location="0-36"` 就是“从第 1 个碱基开始”。）
    #[test]
    fn zero_based_range_starts_at_origin() {
        let e = Edit::Delete { start: 3, end: 4 };
        let parts = transform_range(0, 36, Base::Zero, std::slice::from_ref(&e), 100);
        assert_eq!(ps(&parts), vec![(0, 36)]);
        assert_eq!(ranges_to_text(&parts, 99, Base::Zero), vec!["0-35".to_string()]);
        // 同一起点在 1-based 特征里是 range="1-37"
        let feature = transform_range(1, 37, Base::One, &[e], 100);
        assert_eq!(ranges_to_text(&feature, 99, Base::One), vec!["1-36".to_string()]);
    }

    #[test]
    fn zero_based_wrap_range_round_trips() {
        // 0-based 闭区间 90-4（环状）→ 片段 [90,100) ∪ [0,5)
        let parts = transform_range(90, 4, Base::Zero, &[], 100);
        assert_eq!(ps(&parts), vec![(90, 100), (0, 5)]);
        assert_eq!(ranges_to_text(&parts, 100, Base::Zero), vec!["90-4".to_string()]);
    }

    /// 回归：删掉一个碱基时，**只有真正被删掉的结合位点消失**，
    /// 引物本身和其它位点必须保留（此前第一个引物会整个消失）。
    #[test]
    fn primer_keeps_other_sites_when_one_is_deleted() {
        let xml = b"<?xml version=\"1.0\"?><Primers nextValidID=\"2\">\
            <Primer recentID=\"0\" name=\"Primer 1\" sequence=\"AAA\">\
            <BindingSite location=\"0-36\" boundStrand=\"0\"><Component hybridizedRange=\"0-36\" bases=\"AAA\"/></BindingSite>\
            <BindingSite location=\"469-497\" boundStrand=\"0\"><Component hybridizedRange=\"469-497\" bases=\"AAA\"/></BindingSite>\
            <BindingSite simplified=\"1\" location=\"929-956\" boundStrand=\"0\"><Component hybridizedRange=\"929-956\" bases=\"AAA\"/></BindingSite>\
            </Primer>\
            <Primer recentID=\"1\" name=\"Primer 2\" sequence=\"CCC\">\
            <BindingSite location=\"148-178\" boundStrand=\"0\"><Component hybridizedRange=\"148-178\" bases=\"CCC\"/></BindingSite>\
            </Primer></Primers>";
        // 删掉 0-based 3..4（第 4 个碱基）
        let out = adjust_primer_xml(xml, &[Edit::Delete { start: 3, end: 4 }], 1000);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("name=\"Primer 1\""), "引物 1 被整体丢了：{s}");
        assert!(s.contains("location=\"0-35\""), "{s}");
        assert!(s.contains("location=\"468-496\""), "{s}");
        assert!(s.contains("location=\"928-955\""), "{s}");
        assert!(s.contains("name=\"Primer 2\""), "{s}");
        assert!(s.contains("location=\"147-177\""), "{s}");

        // 位点被整体删掉时：只丢那一个位点，引物与兄弟位点仍在
        let out2 = adjust_primer_xml(xml, &[Edit::Delete { start: 469, end: 498 }], 1000);
        let s2 = String::from_utf8(out2).unwrap();
        assert!(s2.contains("name=\"Primer 1\""), "{s2}");
        assert!(!s2.contains("location=\"469-497\""), "{s2}");
        assert!(!s2.contains("location=\"468-496\""), "{s2}");
        assert!(s2.contains("location=\"0-36\""), "{s2}");
        assert!(s2.contains("location=\"900-927\""), "{s2}");
    }

    #[test]
    fn new_file_round_trips() {
        let seq = b"ATGCGTACGTTAGCCTAGCTAGCTTAAAGGCC".to_vec();
        let f = NewFile::dna(seq.clone())
            .with_name("demo")
            .with_features(vec![Feature::new("ORF1", 1, 12)]);
        let bytes = f.to_bytes();
        let parsed = DnaFile::parse(&bytes).unwrap();
        let rec = parsed.sequence().unwrap();
        assert_eq!(rec.sequence, seq);
        assert_eq!(parsed.header.type_of_sequence, 1);
        assert!(parsed.block(6).unwrap().payload.starts_with(b"<Notes>"));
        let feats = String::from_utf8_lossy(&parsed.block(10).unwrap().payload).to_string();
        assert!(feats.contains("name=\"ORF1\""), "{feats}");
        assert!(feats.contains("range=\"1-12\""), "{feats}");
    }

    /// 三种分子类型的**头部字段与序列块类型**必须与真 SnapGene 一致
    /// （8.2.2 oracle 实测：`.prot` = type 2 / block 21 / flags 0x00；
    /// `.rna` = type 7 / block 32 / flags 0x00（New File 默认单链线性）；
    /// `.dna` = type 1 / block 0 / flags 0x02（双链线性））。
    #[test]
    fn new_file_kinds_match_snapgene_headers() {
        use hx_core::MoleculeKind;
        for (kind, tos, block, flags) in [
            (MoleculeKind::Dna, 1u16, 0u8, 0x02u8),
            (MoleculeKind::Rna, 7, 32, 0x00),
            (MoleculeKind::Protein, 2, 21, 0x00),
        ] {
            let seq = b"MKTXAYIAKGG".to_vec();
            let bytes = NewFile::of(kind, seq.clone()).to_bytes();
            let parsed = DnaFile::parse(&bytes).unwrap();
            assert_eq!(
                parsed.header.type_of_sequence, tos,
                "{kind:?}: 头部 type_of_sequence"
            );
            assert_eq!(parsed.header.molecule_kind(), Some(kind));
            let sb = parsed
                .blocks
                .iter()
                .find(|b| b.block_type == block)
                .unwrap_or_else(|| panic!("{kind:?}: 没有 block {block}"));
            let rec = parsed.sequence().unwrap();
            assert_eq!(rec.sequence, seq, "{kind:?}: 序列");
            assert_eq!(rec.flags.to_byte(), flags, "{kind:?}: flags");
            assert_eq!(sb.payload[0], flags, "{kind:?}: 载荷首字节 = flags");
            // 往返：读回再写出逐字节一致
            assert_eq!(parsed.to_bytes(), bytes, "{kind:?}: 读→写往返");
            let again = crate::write::rewrite_after_edit(&parsed, &[]).unwrap();
            assert_eq!(again.to_bytes(), bytes, "{kind:?}: 空编辑写回");
        }
    }

    #[test]
    fn edit_drops_derived_blocks() {
        let seq = b"ATGCGTACGTTAGCCTAGCTAGCTTAAAGGCC".to_vec();
        let base = NewFile::dna(seq).to_bytes();
        let mut file = DnaFile::parse(&base).unwrap();
        // 手工塞入派生块，验证写回丢弃
        for t in DROPPED_ON_SEQUENCE_CHANGE {
            file.blocks.push(Block::new(*t, vec![1, 2, 3]));
        }
        let out = rewrite_after_edit(&file, &[Edit::Delete { start: 0, end: 3 }]).unwrap();
        for t in DROPPED_ON_SEQUENCE_CHANGE {
            assert!(out.block(*t).is_none(), "block {t} 应被丢弃");
        }
        let rec = out.sequence().unwrap();
        assert_eq!(rec.sequence, b"CGTACGTTAGCCTAGCTAGCTTAAAGGCC");
        // LastModified 已刷新为今天
        let notes = out.block(6).unwrap().payload.clone();
        let notes = String::from_utf8(notes).unwrap();
        let today = snapgene_date(utc_now());
        assert!(notes.contains(&format!(">{today}</LastModified>")), "{notes}");
    }

    #[test]
    fn noop_edit_keeps_everything() {
        let base = NewFile::dna(b"ATGCGT".to_vec()).to_bytes();
        let file = DnaFile::parse(&base).unwrap();
        let out = rewrite_after_edit(&file, &[]).unwrap();
        assert_eq!(out.blocks.len(), file.blocks.len());
    }
}
