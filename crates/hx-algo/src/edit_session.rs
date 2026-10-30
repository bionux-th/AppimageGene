//! # 编辑会话（工作副本 + 编辑历史 + undo/redo + 保存字节）
//!
//! 设计要点（与 P1 写回契约一致，规格 §12.2）：
//!
//! - 会话以**打开（或最近一次保存）时的文件**为基准；所有编辑以
//!   [`hx_formats::write::Edit`] 列表按序记录；
//! - 工作状态（序列 / 特征 / 引物）= 基准 + 编辑列表**重放**得到：
//!   序列用 [`apply_edits`]；特征/引物用块 5/10 XML 的区间变换后**重新解析**
//!   ——与写回走同一条代码路径，显示与落盘逐字节一致；
//! - undo/redo 只挪动编辑列表两端，不保存快照；
//! - 保存时：[`EditSession::save_bytes`] 直接给出待写字节（无编辑 ⇒ 原样字节，
//!   保证未编辑文件读→写逐字节一致）。
//! - 线性双链 DNA 的底链独有末端碱基保留在原生工作副本；公开序列、注释和
//!   编辑使用顶链显示坐标，在重放/保存边界统一映射，不改动核心 SeqRecord。

use hx_core::{MoleculeFlags, MoleculeKind, SeqRecord};
use hx_formats::annot::{parse_features, parse_primers, FeatureAnnotation, PrimerAnnotation};
use hx_formats::feature_edit::{apply_feature_edits, parse_specs, FeatureEdit, FeatSpec};
use hx_formats::primer_edit::{apply_primer_edits, PrimerEdit};
use hx_formats::write::{adjust_feature_xml, adjust_primer_xml, apply_edits, rewrite_after_edit, Edit};
use hx_formats::DnaFile;

/// 会话里的一次操作（按时间顺序重放；序列编辑与特征编辑共用一条列表，
/// 因为两者互相影响：特征坐标随序列编辑平移，特征编辑也改工作 XML）。
#[derive(Debug, Clone, PartialEq)]
pub enum Op {
    /// 序列编辑（[`Edit`]）。
    Seq(Edit),
    /// 特征编辑（[`FeatureEdit`]）。
    Feat(FeatureEdit),
    /// 引物编辑（[`PrimerEdit`]）。
    Primer(PrimerEdit),
    Notes {tag:String,value:String},
}

impl Op {
    /// 人读描述（状态栏 / undo 提示）。
    pub fn label(&self) -> String {
        match self {
            Op::Seq(e) => edit_label(e),
            Op::Feat(f) => f.label(),
            Op::Primer(p) => p.label(),
            Op::Notes {tag,..}=>format!("Edit {tag}"),
        }
    }

    /// 是否为“注释”编辑（特征或引物）——序列编辑之外的都算。
    fn is_annotation(&self) -> bool {
        matches!(self, Op::Feat(_) | Op::Primer(_))
    }
}

/// 一个编辑会话（UI 无关）。
pub struct EditSession {
    /// 基准文件（保存/写回的输入）。
    file: DnaFile,
    /// 基准字节（无编辑时保存逐字节原样）。
    bytes: Vec<u8>,
    kind: MoleculeKind,
    flags: MoleculeFlags,
    /// 基准序列。
    base_seq: Vec<u8>,
    /// 基准块 5/10 的 XML 载荷。
    base_features_xml: Option<Vec<u8>>,
    base_primers_xml: Option<Vec<u8>>,
    /// 原生坐标的工作副本（基准 + `ops` 重放）；公开的 sequence/features/primers
    /// 使用顶链坐标。sequence() 借用此串，不再复制一份完整顶链。
    stored_sequence: Vec<u8>,
    stored_features: Vec<FeatureAnnotation>,
    stored_primers: Vec<PrimerAnnotation>,
    ends: hx_formats::dna::EndProperties,
    features: Vec<FeatureAnnotation>,
    primers: Vec<PrimerAnnotation>,
    /// 工作块 10 载荷（含特征编辑）。
    features_xml: Option<Vec<u8>>,
    /// 工作块 5 载荷（含引物编辑）。
    primers_xml: Option<Vec<u8>>,
    /// 工作块 10 的可编辑模型（对话框用）。
    specs: Vec<FeatSpec>,
    ops: Vec<Op>,
    undone: Vec<Op>,
}

impl EditSession {
    /// 从已解析文件（及其原始字节）建立会话。
    pub fn new(file: DnaFile, bytes: Vec<u8>) -> Result<Self, String> {
        let kind = file
            .header
            .molecule_kind()
            .ok_or_else(|| "unknown molecule type".to_string())?;
        let seq_type = kind.sequence_block_type();
        let block = file
            .block(seq_type)
            .ok_or_else(|| "missing sequence block".to_string())?;
        let rec = SeqRecord::from_block_payload(kind, &block.payload)
            .ok_or_else(|| "empty sequence block".to_string())?;
        let base_features_xml = file.block(10).map(|b| b.payload.clone());
        // 块 5 可能整块缺失（手工构造的最小文件）；给一份空骨架，
        // 这样「加引物」能落地（无引物编辑时保存仍走原字节，不受影响）。
        let base_primers_xml = Some(
            file.block(5)
                .map(|b| b.payload.clone())
                .unwrap_or_else(|| hx_formats::write::empty_primers_xml().as_bytes().to_vec()),
        );
        let ends = if kind == MoleculeKind::Dna && !rec.flags.topology().is_circular()
            && rec.flags.strandedness() == hx_core::Strandedness::Double
        { file.end_properties() } else { Default::default() };
        let mut s = Self {
            base_seq: rec.sequence.clone(),
            flags: rec.flags,
            kind,
            file,
            bytes,
            base_features_xml,
            base_primers_xml,
            stored_sequence: rec.sequence,
            stored_features: Vec::new(),
            stored_primers: Vec::new(),
            ends,
            features: Vec::new(),
            primers: Vec::new(),
            features_xml: None,
            primers_xml: None,
            specs: Vec::new(),
            ops: Vec::new(),
            undone: Vec::new(),
        };
        s.rebuild()?;
        Ok(s)
    }

    // ---- 只读访问 ----

    pub fn file(&self) -> &DnaFile {
        &self.file
    }

    /// 基准字节（未编辑时即保存字节）。
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub fn kind(&self) -> MoleculeKind {
        self.kind
    }

    pub fn flags(&self) -> MoleculeFlags {
        self.flags
    }

    pub fn circular(&self) -> bool {
        self.flags.topology().is_circular()
    }

    /// 当前顶链序列（工作副本；底链独有的末端标记不计入编辑坐标）。
    pub fn sequence(&self) -> &[u8] {
        &self.stored_sequence[self.ends.top_range(self.stored_sequence.len())]
    }

    /// 克隆 / 原生文件构造使用存储坐标，不可把已裁掉标记的顶链再附上末端属性。
    pub fn stored_sequence(&self) -> &[u8] { &self.stored_sequence }
    pub fn stored_features(&self) -> &[FeatureAnnotation] { &self.stored_features }
    pub fn stored_primers(&self) -> &[PrimerAnnotation] { &self.stored_primers }
    pub fn end_properties(&self) -> &hx_formats::dna::EndProperties { &self.ends }

    pub fn features(&self) -> &[FeatureAnnotation] {
        &self.features
    }

    pub fn primers(&self) -> &[PrimerAnnotation] {
        &self.primers
    }

    pub fn len(&self) -> usize {
        self.sequence().len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence().is_empty()
    }

    /// 基准顶链长度（编辑前的 UI 坐标基准）。
    pub fn base_len(&self) -> usize {
        self.ends.top_range(self.base_seq.len()).len()
    }

    /// 序列编辑条数（不含特征编辑）——`hx-cli edit` 的等价物。
    pub fn seq_edits(&self) -> Vec<Edit> {
        self.ops
            .iter()
            .filter_map(|o| match o {
                Op::Seq(e) => Some(e.clone()),
                Op::Feat(_) | Op::Primer(_) | Op::Notes{..} => None,
            })
            .collect()
    }

    fn storage_seq_edits(&self) -> Vec<Edit> {
        let offset = self.ends.top_range(self.base_seq.len()).start;
        self.seq_edits().iter().map(|e| offset_edit(e, offset)).collect()
    }

    /// 会话操作总数（序列 + 特征）。
    pub fn op_count(&self) -> usize {
        self.ops.len()
    }

    /// 特征编辑条数。
    pub fn feature_edit_count(&self) -> usize {
        self.ops
            .iter()
            .filter(|o| matches!(o, Op::Feat(_)))
            .count()
    }

    /// 引物编辑条数。
    pub fn primer_edit_count(&self) -> usize {
        self.ops
            .iter()
            .filter(|o| matches!(o, Op::Primer(_)))
            .count()
    }

    /// 注释编辑（特征 + 引物）条数。
    pub fn annotation_edit_count(&self) -> usize {
        self.ops.iter().filter(|o| o.is_annotation()).count()
    }

    /// 原生存储坐标的工作块 10 载荷（含特征编辑；无块 10 时为 `None`）。
    pub fn features_xml(&self) -> Option<&[u8]> {
        self.features_xml.as_deref()
    }

    pub fn primers_xml(&self) -> Option<&[u8]> { self.primers_xml.as_deref() }

    /// 工作块 10 的可编辑模型（对话框用；索引与 [`Self::features`] 一一对应）。
    pub fn feature_specs(&self) -> &[FeatSpec] {
        &self.specs
    }

    /// 某个特征的当前模型（按 `recentID`）。
    pub fn feature_spec(&self, recent_id: u32) -> Option<&FeatSpec> {
        self.specs.iter().find(|s| s.recent_id == recent_id)
    }

    pub fn is_dirty(&self) -> bool {
        !self.ops.is_empty()
    }

    pub fn can_undo(&self) -> bool {
        !self.ops.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.undone.is_empty()
    }

    /// 块数（基准文件）。
    pub fn block_count(&self) -> usize {
        self.file.blocks.len()
    }

    // ---- 编辑 ----

    /// 施加一条序列编辑（坐标为**当前**序列坐标，0-based 半开）；失败则不改变状态。
    pub fn apply(&mut self, edit: Edit) -> Result<(), String> {
        self.push_op(Op::Seq(edit))
    }

    /// 施加一条特征编辑。
    pub fn apply_feature(&mut self, edit: FeatureEdit) -> Result<(), String> {
        self.push_op(Op::Feat(edit))
    }

    /// 施加一条引物编辑。
    pub fn apply_primer(&mut self, edit: PrimerEdit) -> Result<(), String> {
        self.push_op(Op::Primer(edit))
    }

    pub fn recalculate_hybridization(&mut self, params:hx_formats::annot::HybridizationParams)->Result<(),String>{
        if self.ends.upstream_stickiness!=0||self.ends.downstream_stickiness!=0{return Err("Recalculating binding sites on sticky-ended sources is not yet supported.".into());}
        let sites=self.primers.iter().map(|p|{
            let bindings=crate::hybridization::search(self.sequence(),p.sequence.as_bytes(),self.circular(),params);
            let sites=bindings.into_iter().map(|b|{
                let n=b.site.annealed_bases.min(p.sequence.len());
                hx_formats::primer_edit::PrimerSiteSpec{range:b.site.range,strand:b.site.strand,annealed:if b.five_prime{p.sequence[..n].to_owned()}else{p.sequence[p.sequence.len()-n..].to_owned()},tm:b.site.melting_temperature.unwrap_or(0.0)}
            }).collect();(p.name.clone(),sites)
        }).collect();
        self.push_op(Op::Primer(PrimerEdit::Recalculate{params,sites}))
    }

    pub fn notes_xml(&self)->Vec<u8> {
        let mut xml=self.file.block(6).map(|b|b.payload.clone()).unwrap_or_else(||b"<Notes></Notes>".to_vec());
        for op in &self.ops {if let Op::Notes{tag,value}=op {xml=hx_formats::notes::set_field(&xml,tag,value);}}
        xml
    }

    pub fn set_note(&mut self,tag:&str,value:&str)->Result<(),String> {
        if tag.is_empty() || !tag.bytes().all(|b|b.is_ascii_alphanumeric()||b==b'_') {return Err("Invalid Notes field".into());}
        self.push_op(Op::Notes{tag:tag.into(),value:value.into()})
    }

    fn push_op(&mut self, op: Op) -> Result<(), String> {
        self.ops.push(op);
        match self.rebuild() {
            Ok(()) => {
                self.undone.clear();
                Ok(())
            }
            Err(e) => {
                self.ops.pop();
                Err(e)
            }
        }
    }

    /// 撤销一步；返回被撤销编辑的描述（无历史返回 `None`）。
    pub fn undo(&mut self) -> Option<String> {
        let e = self.ops.pop()?;
        let label = e.label();
        self.undone.push(e);
        // 该编辑此前已成功重放过，这里不可能失败
        let _ = self.rebuild();
        Some(label)
    }

    /// 重做一步；返回被重做编辑的描述。
    pub fn redo(&mut self) -> Option<String> {
        let e = self.undone.pop()?;
        let label = e.label();
        self.ops.push(e);
        let _ = self.rebuild();
        Some(label)
    }

    /// 待保存字节：无编辑 ⇒ 原样字节；有编辑 ⇒ 按写回契约重建整份文件（不含历史树）。
    ///
    /// 序列编辑走 [`rewrite_after_edit`]（丢弃块 3/7/11/34 等缓存）；
    /// 特征编辑不改序列，只把**工作块 10** 覆盖回去，并丢弃历史树
    /// （SnapGene 的历史节点语义未解，不伪造）。
    pub fn save_bytes(&self) -> Result<Vec<u8>, String> {
        if self.ops.is_empty() {
            return Ok(self.bytes.clone());
        }
        let f = self.build_file()?;
        Ok(f.to_bytes())
    }

    /// 待保存字节（含 **block 7 历史树**，与 `hx-cli edit` 同一条实现）。
    ///
    /// `doc_name` = 保存后的文件名（新根显示名），`src_name` = 打开时的文件名（父节点）。
    /// 返回 `(字节, 历史树节点数)`。序列操作逐项记录；纯注释编辑保留原历史树。
    pub fn save_bytes_with_history(
        &self,
        doc_name: &str,
        src_name: &str,
        parent_resurrectable: Option<u8>,
    ) -> Result<(Vec<u8>, usize), String> {
        if self.ops.is_empty() {
            return Ok((self.bytes.clone(), 0));
        }
        let seq_edits = self.storage_seq_edits();
        let mut f = self.build_file()?;
        if seq_edits.is_empty() {
            // Annotation/notes changes must not destroy molecular lineage.
            // Do not invent an unsupported SG molecular operation for them.
            for b in self.file.blocks.iter().filter(|b| matches!(b.block_type, 7 | 11 | 29 | 30)) {
                if !f.blocks.iter().any(|x| x.block_type == b.block_type && x.payload == b.payload) {
                    f.blocks.push(b.clone());
                }
            }
            return Ok((f.to_bytes(), 0));
        }
        // Record each accepted sequence operation in chronological order. Undo
        // removes operations from this list, so it cannot leave phantom nodes.
        let mut replay = EditSession::new(self.file.clone(), self.bytes.clone())?;
        let mut previous = self.file.clone();
        let mut nodes = 0;
        for op in &self.ops {
            let offset = replay.ends.top_range(replay.stored_sequence.len()).start;
            replay.push_op(op.clone())?;
            if let Op::Seq(edit) = op {
                let mut next = replay.build_file()?;
                nodes = hx_formats::write::insert_history_block(
                    &mut next, &previous, &[offset_edit(edit, offset)], doc_name,
                    if nodes == 0 { src_name } else { doc_name }, parent_resurrectable,
                ).map_err(|e| e.to_string())?;
                previous = next;
            }
        }
        f.blocks.retain(|b| b.block_type != 7);
        if let Some(history) = previous.block(7) {
            let pos = f.blocks.iter().position(|b| b.block_type == 6).map(|i| i + 1).unwrap_or(f.blocks.len());
            f.blocks.insert(pos, history.clone());
        }
        Ok((f.to_bytes(), nodes))
    }

    /// 组装待写文件：序列编辑 → [`rewrite_after_edit`]；块 10 → 工作载荷。
    fn build_file(&self) -> Result<DnaFile, String> {
        let seq_edits = self.storage_seq_edits();
        let mut f = rewrite_after_edit(&self.file, &seq_edits).map_err(|e| e.to_string())?;
        if self.feature_edit_count() > 0 {
            let Some(xml) = &self.features_xml else {
                return Ok(f);
            };
            match f.blocks.iter_mut().find(|b| b.block_type == 10) {
                Some(b) => b.payload = xml.clone(),
                None => {
                    // 原文件没有块 10：插在块 8 之后（没有则插在序列块之后）
                    let at = f
                        .blocks
                        .iter()
                        .position(|b| b.block_type == 8)
                        .map(|i| i + 1)
                        .unwrap_or(1);
                    f.blocks.insert(at, hx_formats::block::Block::new(10, xml.clone()));
                }
            }
            // 历史树（块 7/11）承载的是**旧**特征集，纯特征编辑时丢弃
            if seq_edits.is_empty() {
                f.blocks.retain(|b| b.block_type != 7 && b.block_type != 11);
            }
        }
        if self.primer_edit_count() > 0 {
            let Some(xml) = &self.primers_xml else {
                return Ok(f);
            };
            match f.blocks.iter_mut().find(|b| b.block_type == 5) {
                Some(b) => b.payload = xml.clone(),
                None => {
                    // 原文件没有块 5：按块号顺序插在第一个 > 5 的块之前
                    let at = f
                        .blocks
                        .iter()
                        .position(|b| b.block_type > 5)
                        .unwrap_or(f.blocks.len());
                    f.blocks.insert(at, hx_formats::block::Block::new(5, xml.clone()));
                }
            }
            if seq_edits.is_empty() {
                f.blocks.retain(|b| b.block_type != 7 && b.block_type != 11);
            }
        }
        if self.ops.iter().any(|op|matches!(op,Op::Notes{..})) {
            let notes=self.notes_xml();
            if let Some(b)=f.blocks.iter_mut().find(|b|b.block_type==6){b.payload=notes;}
            else{f.blocks.push(hx_formats::Block::new(6,notes));}
        }
        Ok(f)
    }

    /// 保存成功后以新文件为新基准（清空历史）。
    pub fn commit_saved(&mut self, bytes: Vec<u8>) -> Result<(), String> {
        let file = DnaFile::parse(&bytes).map_err(|e| format!("parse saved file: {e}"))?;
        let mut fresh = EditSession::new(file, bytes)?;
        // 路径等外部状态由调用方维护，这里只接替文件与工作状态
        std::mem::swap(self, &mut fresh);
        Ok(())
    }

    /// 用操作列表重放工作状态（序列 → 特征/引物 XML，按时间顺序）。
    fn rebuild(&mut self) -> Result<(), String> {
        let mut seq = self.base_seq.clone();
        let mut fx = self.base_features_xml.clone();
        let mut px = self.base_primers_xml.clone();
        for op in &self.ops {
            match op {
                Op::Seq(e) => {
                    let old_len = seq.len();
                    let top = self.ends.top_range(old_len);
                    // 先在顶链坐标验证，不能让越界 UI 编辑落到隐藏标记上。
                    let valid = match e {
                        Edit::Insert { at, .. } => *at <= top.len(),
                        Edit::Delete { start, end } | Edit::Replace { start, end, .. } => start <= end && *end <= top.len(),
                    };
                    if !valid { return Err(format!("Edit lies outside the displayed sequence ({} bases).", top.len())); }
                    validate_terminal_edit(e, &self.ends, top.len())?;
                    let native = offset_edit(e, top.start);
                    let one = std::slice::from_ref(&native);
                    seq = apply_edits(&seq, one).map_err(|e| e.to_string())?;
                    if let Some(x) = &fx {
                        fx = Some(adjust_feature_xml(x, one, old_len));
                    }
                    if let Some(x) = &px {
                        px = Some(adjust_primer_xml(x, one, old_len));
                    }
                }
                Op::Feat(f) => {
                    if let Some(x) = &fx {
                        let top = self.ends.top_range(seq.len());
                        let native = storage_feature_edit(f, x, seq.len(), top)?;
                        fx = Some(apply_feature_edits(x, std::slice::from_ref(&native)));
                    }
                }
                Op::Primer(p) => {
                    if let Some(x) = &px {
                        let native = storage_primer_edit(p, x, self.ends.top_range(seq.len()))?;
                        px = Some(apply_primer_edits(x, std::slice::from_ref(&native)));
                    }
                }
                Op::Notes{..}=>{},
            }
        }
        let top = self.ends.top_range(seq.len());
        let display_fx = fx.as_ref().map(|x| project_annotation_xml(x, seq.len(), top.clone(), false));
        let display_px = px.as_ref().map(|x| project_annotation_xml(x, seq.len(), top.clone(), true));
        self.stored_sequence = seq;
        self.stored_features = fx.as_deref().map(parse_features).unwrap_or_default();
        self.stored_primers = px.as_deref().map(parse_primers).unwrap_or_default();
        self.features = display_fx.as_deref().map(parse_features).unwrap_or_default();
        self.primers = display_px.as_deref().map(parse_primers).unwrap_or_default();
        self.specs = display_fx.as_deref().map(parse_specs).unwrap_or_default();
        self.features_xml = fx;
        self.primers_xml = px;
        Ok(())
    }
}

fn offset_edit(edit: &Edit, offset: usize) -> Edit {
    match edit {
        Edit::Insert { at, seq } => Edit::Insert { at: at + offset, seq: seq.clone() },
        Edit::Delete { start, end } => Edit::Delete { start: start + offset, end: end + offset },
        Edit::Replace { start, end, seq } => Edit::Replace { start: start + offset, end: end + offset, seq: seq.clone() },
    }
}

/// 外侧标记及突出端长度不因普通双链编辑改变。突出端本身的化学编辑暂不支持，
/// 明确拒绝而不是静默保留错误的 stickiness。
fn validate_terminal_edit(edit: &Edit, ends: &hx_formats::dna::EndProperties, len: usize) -> Result<(), String> {
    if ends.upstream_stickiness == 0 && ends.downstream_stickiness == 0 { return Ok(()); }
    let left = ends.upstream_stickiness.max(0) as usize;
    let right = len.saturating_sub(ends.downstream_stickiness.min(0).unsigned_abs() as usize);
    let (a, b, added) = match edit {
        Edit::Insert { at, seq } => (*at, *at, seq.len()),
        Edit::Delete { start, end } => (*start, *end, 0),
        Edit::Replace { start, end, seq } => (*start, *end, seq.len()),
    };
    if a < left || b > right || right.saturating_sub(left).saturating_sub(b - a) + added == 0 {
        return Err("Editing sticky overhangs or removing the entire double-stranded region is not supported. Blunt the fragment first.".into());
    }
    Ok(())
}

/// 只构造显示用 XML；原生工作 XML 始终保留标记区和未知属性。
fn project_annotation_xml(xml: &[u8], len: usize, top: std::ops::Range<usize>, primer: bool) -> Vec<u8> {
    if top.start == 0 && top.end == len { return xml.to_vec(); }
    let edits = [Edit::Delete { start: top.end, end: len }, Edit::Delete { start: 0, end: top.start }];
    if primer { adjust_primer_xml(xml, &edits, len) } else { adjust_feature_xml(xml, &edits, len) }
}

fn storage_feature_edit(
    edit: &FeatureEdit, xml: &[u8], len: usize, top: std::ops::Range<usize>,
) -> Result<FeatureEdit, String> {
    if top.start == 0 && top.end == len { return Ok(edit.clone()); }
    let mut spec = match edit {
        FeatureEdit::Add(spec) | FeatureEdit::Update(spec) => spec.clone(),
        _ => return Ok(edit.clone()),
    };
    for seg in &spec.segments {
        if seg.start == 0 || seg.start > top.len() || seg.end == 0 || seg.end > top.len() {
            return Err("Feature coordinates lie outside the displayed sequence.".into());
        }
    }
    let original = parse_specs(xml).into_iter().find(|s| s.recent_id == spec.recent_id);
    let displayed = parse_specs(&project_annotation_xml(xml, len, top.clone(), false))
        .into_iter().find(|s| s.recent_id == spec.recent_id);
    for seg in &mut spec.segments { seg.start += top.start; seg.end += top.start; }
    if matches!(edit, FeatureEdit::Update(_)) {
        if let (Some(original), Some(displayed)) = (original, displayed) {
            // 元数据编辑不得把被显示窗口裁掉的部分写没。
            let unchanged = match edit {
                FeatureEdit::Update(s) => s.segments.len() == displayed.segments.len()
                    && s.segments.iter().zip(&displayed.segments).all(|(a, b)| a.start == b.start && a.end == b.end),
                _ => false,
            };
            if unchanged {
                let input = match edit { FeatureEdit::Update(s) => s, _ => unreachable!() };
                let mut raw_segments = original.segments;
                for raw in &mut raw_segments {
                    let a = raw.start.saturating_sub(1).max(top.start);
                    let b = raw.end.min(top.end);
                    if a >= b { continue; }
                    if let Some(visible) = input.segments.iter().find(|s| {
                        s.start == a - top.start + 1 && s.end == b - top.start
                    }) {
                        raw.color = visible.color.clone();
                        raw.translated = visible.translated;
                    }
                }
                spec.segments = raw_segments;
            }
        }
        Ok(FeatureEdit::Update(spec))
    } else { Ok(FeatureEdit::Add(spec)) }
}

fn storage_primer_edit(edit: &PrimerEdit, xml: &[u8], top: std::ops::Range<usize>) -> Result<PrimerEdit, String> {
    if let PrimerEdit::Replace { name, .. } = edit {
        if parse_primers(xml).iter().find(|p| p.name == *name).is_some_and(|p| {
            p.sites.iter().any(|s| s.range.start < top.start || s.range.end > top.end)
        }) {
            return Err("Editing a primer that extends into bottom-strand end markers is not supported. Its original binding sites have been preserved.".into());
        }
    }
    let mut native = edit.clone();
    let spec = match &mut native {
        PrimerEdit::Add(spec) | PrimerEdit::Replace { spec, .. } => spec,
        _ => return Ok(native),
    };
    for site in &mut spec.sites {
        if site.range.start >= top.len() || site.range.end > top.len() {
            return Err("Primer coordinates lie outside the displayed sequence.".into());
        }
        site.range.start += top.start; site.range.end += top.start;
    }
    Ok(native)
}

/// 序列编辑的人读描述（状态栏用）。
pub fn edit_label(e: &Edit) -> String {
    match e {
        Edit::Delete { start, end } => format!("delete {}..{}", start + 1, end),
        Edit::Insert { at, seq } => format!("insert {} bp at {}", seq.len(), at + 1),
        Edit::Replace { start, end, seq } => {
            format!("replace {}..{} with {} bp", start + 1, end, seq.len())
        }
    }
}

/// 解析用户粘贴的序列文本：只保留 ASCII 字母，其余（空白/数字/符号）丢弃。
pub fn parse_seq_text(s: &str) -> Vec<u8> {
    s.bytes().filter(|b| b.is_ascii_alphabetic()).collect()
}

/// IUPAC 互补表（含简并码；未知字符原样）。
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
        b'V' => b'B',
        b'v' => b'b',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

/// 反向互补。
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| complement(*b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_notes_edits_are_undoable_and_preserve_unrelated_blocks_and_sequence() {
        let new=hx_formats::write::NewFile::dna(b"ACGT".repeat(25));
        let original=new.to_bytes();let file=DnaFile::parse(&original).unwrap();
        let mut session=EditSession::new(file.clone(),original.clone()).unwrap();
        session.set_note("Description","new & description").unwrap();
        assert!(session.is_dirty());assert!(session.can_undo());
        assert_eq!(hx_formats::notes::field(&session.notes_xml(),"Description").as_deref(),Some("new & description"));
        let saved=DnaFile::parse(&session.save_bytes().unwrap()).unwrap();
        for b in file.blocks.iter().filter(|b|b.block_type!=6) {assert_eq!(saved.block(b.block_type).unwrap().payload,b.payload);}
        session.undo().unwrap();assert!(!session.is_dirty());assert_eq!(session.save_bytes().unwrap(),original);
        session.redo().unwrap();assert!(session.is_dirty());
        session.apply(Edit::Insert{at:10,seq:b"TT".to_vec()}).unwrap();
        let saved=DnaFile::parse(&session.save_bytes().unwrap()).unwrap();
        assert_eq!(saved.sequence().unwrap().len(),102);
        assert_eq!(hx_formats::notes::field(&saved.block(6).unwrap().payload,"Description").as_deref(),Some("new & description"));
    }

    fn sticky_sample() -> EditSession {
        let mut file = DnaFile::parse(&hx_formats::write::NewFile::dna(b"TGGGGAATTCCCCAA".to_vec()).to_bytes()).unwrap();
        for block in &mut file.blocks {
            match block.block_type {
                8 => block.payload = b"<AdditionalSequenceProperties extra=\"keep\"><UpstreamStickiness>-1</UpstreamStickiness><DownstreamStickiness>2</DownstreamStickiness><UpstreamModification>Phosphorylated</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>".to_vec(),
                10 => block.payload = b"<Features nextValidID=\"2\"><Feature recentID=\"0\" name=\"visible\" type=\"misc_feature\" extra=\"keep\"><Segment range=\"2-5\" color=\"#ff0000\"/></Feature><Feature recentID=\"1\" name=\"marker\" type=\"misc_feature\"><Segment range=\"1-1\" color=\"#00ff00\"/></Feature></Features>".to_vec(),
                5 => block.payload = b"<Primers nextValidID=\"1\"><Primer recentID=\"0\" name=\"p\" sequence=\"GGGG\"><BindingSite location=\"1-4\" boundStrand=\"0\" annealedBases=\"GGGG\"><Component hybridizedRange=\"1-4\" bases=\"GGGG\"/></BindingSite></Primer></Primers>".to_vec(),
                _ => {},
            }
        }
        let bytes = file.to_bytes();
        EditSession::new(file, bytes).unwrap()
    }

    #[test]
    fn sticky_editor_uses_top_coordinates_without_discarding_native_markers() {
        let mut s = sticky_sample();
        assert_eq!(s.len(), 12, "editor must not count the three bottom-strand markers");
        assert_eq!(s.sequence(), b"GGGGAATTCCCC");
        assert_eq!(s.features().len(), 1);
        assert_eq!(s.features()[0].segments[0].range, hx_formats::annot::Range::new(0, 4));
        assert_eq!(s.feature_spec(0).unwrap().segments[0].start, 1);
        assert_eq!(s.primers()[0].sites[0].range, hx_formats::annot::Range::new(0, 4));
        assert_eq!(s.save_bytes().unwrap(), s.bytes());
        s.apply(Edit::Insert { at: 0, seq: b"AC".to_vec() }).unwrap();
        assert_eq!(s.sequence(), b"ACGGGGAATTCCCC");
        let saved = DnaFile::parse(&s.save_bytes().unwrap()).unwrap();
        assert_eq!(saved.sequence().unwrap().sequence, b"TACGGGGAATTCCCCAA");
        assert_eq!(saved.block(8).unwrap().payload, s.file().block(8).unwrap().payload);
        assert_eq!(s.features()[0].segments[0].range, hx_formats::annot::Range::new(2, 6));
        assert_eq!(s.primers()[0].sites[0].range, hx_formats::annot::Range::new(2, 6));
        s.undo().unwrap();
        assert_eq!(s.sequence(), b"GGGGAATTCCCC");
        assert_eq!(s.save_bytes().unwrap(), s.bytes());
    }

    #[test]
    fn sticky_editor_append_delete_replace_save_and_commit_use_display_coordinates() {
        let mut s = sticky_sample();
        s.apply(Edit::Insert { at: s.len(), seq: b"GG".to_vec() }).unwrap();
        s.apply(Edit::Replace { start: 4, end: 8, seq: b"AT".to_vec() }).unwrap();
        s.apply(Edit::Delete { start: 0, end: 1 }).unwrap();
        assert_eq!(s.sequence(), b"GGGATCCCCGG");
        assert_eq!(s.stored_sequence(), b"TGGGATCCCCGGAA");
        assert_eq!(s.stored_features().len(), 2, "hidden marker annotation must survive edits");
        let (bytes, nodes) = s.save_bytes_with_history("edited", "original", None).unwrap();
        assert!(nodes > 0);
        let mut loaded = EditSession::new(DnaFile::parse(&bytes).unwrap(), bytes.clone()).unwrap();
        assert_eq!(loaded.sequence(), s.sequence());
        assert_eq!(loaded.stored_sequence(), s.stored_sequence());
        loaded.commit_saved(bytes).unwrap();
        loaded.apply(Edit::Insert { at: 0, seq: b"C".to_vec() }).unwrap();
        assert_eq!(loaded.stored_sequence(), b"TCGGGATCCCCGGAA");
        assert_eq!(loaded.base_len(), 11);
        s.undo().unwrap(); s.undo().unwrap(); s.undo().unwrap();
        assert_eq!(s.save_bytes().unwrap(), s.bytes());
        s.redo().unwrap(); s.redo().unwrap(); s.redo().unwrap();
        assert_eq!(s.sequence(), b"GGGATCCCCGG");
    }

    #[test]
    fn sticky_feature_and_primer_edits_are_mapped_back_to_native_coordinates() {
        use hx_formats::annot::Range;
        use hx_formats::primer_edit::{PrimerSiteSpec, PrimerSpec};
        let mut s = sticky_sample();
        let mut spec = s.feature_spec(0).unwrap().clone();
        spec.name = "renamed".into();
        s.apply_feature(FeatureEdit::Update(spec)).unwrap();
        let mut added = FeatSpec::new(5, 8); added.name = "new".into();
        s.apply_feature(FeatureEdit::Add(added)).unwrap();
        let mut primer = PrimerSpec::new("new primer", "AATT");
        primer.sites.push(PrimerSiteSpec { range: Range::new(4, 8), strand: 0, annealed: "AATT".into(), tm: 40.0 });
        s.apply_primer(PrimerEdit::Add(Box::new(primer))).unwrap();
        s.apply(Edit::Insert { at: 0, seq: b"G".to_vec() }).unwrap();
        let saved = DnaFile::parse(&s.save_bytes().unwrap()).unwrap();
        let features = parse_features(&saved.block(10).unwrap().payload);
        let primer = parse_primers(&saved.block(5).unwrap().payload);
        assert_eq!(s.features()[0].name, "renamed");
        assert_eq!(s.features()[1].segments[0].range, Range::new(5, 9));
        assert_eq!(features[2].segments[0].range, Range::new(6, 10));
        assert_eq!(s.primers()[1].sites[0].range, Range::new(5, 9));
        assert_eq!(primer[1].sites[0].range, Range::new(6, 10));
        assert!(String::from_utf8_lossy(&saved.block(10).unwrap().payload).contains("extra=\"keep\""));
        assert_eq!(saved.block(8).unwrap().payload, s.file().block(8).unwrap().payload);
    }

    #[test]
    fn sticky_editor_rejects_out_of_view_and_overhang_edits_transactionally() {
        let mut s = sticky_sample();
        for edit in [
            Edit::Insert { at: 13, seq: b"A".to_vec() },
            Edit::Delete { start: 0, end: 13 },
            Edit::Replace { start: 0, end: 12, seq: Vec::new() },
        ] {
            assert!(s.apply(edit).is_err());
            assert!(!s.is_dirty());
            assert_eq!(s.sequence(), b"GGGGAATTCCCC");
        }
        let bytes = hx_formats::write::NewFile::dna(b"AATTGGGGTT".to_vec()).to_bytes();
        let mut file = DnaFile::parse(&bytes).unwrap();
        file.blocks.iter_mut().find(|b| b.block_type == 8).unwrap().payload =
            b"<AdditionalSequenceProperties><UpstreamStickiness>4</UpstreamStickiness><DownstreamStickiness>-2</DownstreamStickiness></AdditionalSequenceProperties>".to_vec();
        let mut s = EditSession::new(file.clone(), file.to_bytes()).unwrap();
        for edit in [Edit::Delete { start: 0, end: 1 }, Edit::Insert { at: 10, seq: b"A".to_vec() }] {
            assert!(s.apply(edit).unwrap_err().contains("sticky overhangs"));
            assert!(!s.is_dirty());
        }
        s.apply(Edit::Insert { at: 4, seq: b"AC".to_vec() }).unwrap();
        assert_eq!(s.sequence(), b"AATTACGGGGTT");
        assert_eq!((s.end_properties().upstream_stickiness, s.end_properties().downstream_stickiness), (4, -2));
    }

    #[test]
    fn sticky_feature_rename_preserves_segments_hidden_in_marker_regions() {
        let mut s = sticky_sample();
        let mut bytes = s.save_bytes().unwrap();
        let mut file = DnaFile::parse(&bytes).unwrap();
        file.blocks.iter_mut().find(|b| b.block_type == 10).unwrap().payload =
            b"<Features nextValidID=\"1\"><Feature recentID=\"0\" name=\"mixed\" type=\"misc_feature\"><Segment range=\"1-1\" color=\"#00ff00\"/><Segment range=\"2-5\" color=\"#ff0000\"/></Feature></Features>".to_vec();
        bytes = file.to_bytes(); s = EditSession::new(file, bytes).unwrap();
        assert_eq!(s.features()[0].segments.len(), 1);
        let mut spec = s.feature_spec(0).unwrap().clone(); spec.name = "renamed".into();
        s.apply_feature(FeatureEdit::Update(spec)).unwrap();
        assert_eq!(s.stored_features()[0].segments.len(), 2);
        assert_eq!(s.stored_features()[0].segments[0].range, hx_formats::annot::Range::new(0, 1));
    }

    #[test]
    fn all_terminal_polarities_use_the_same_display_and_storage_contract() {
        for (up, down, expected) in [(-2, 2, "GGAATTCC"), (-2, -2, "GGAATTCCAA"), (2, 2, "TTGGAATTCC"), (2, -2, "TTGGAATTCCAA")] {
            let mut s = sticky_sample();
            let mut file = s.file().clone();
            file.blocks.iter_mut().find(|b| b.block_type == 0).unwrap().payload =
                [vec![2], b"TTGGAATTCCAA".to_vec()].concat();
            file.blocks.iter_mut().find(|b| b.block_type == 8).unwrap().payload = format!(
                "<AdditionalSequenceProperties><UpstreamStickiness>{up}</UpstreamStickiness><DownstreamStickiness>{down}</DownstreamStickiness></AdditionalSequenceProperties>"
            ).into_bytes();
            let bytes = file.to_bytes(); s = EditSession::new(file, bytes.clone()).unwrap();
            assert_eq!(s.sequence(), expected.as_bytes());
            assert_eq!(s.save_bytes().unwrap(), bytes);
            let at = 4;
            s.apply(Edit::Insert { at, seq: b"AC".to_vec() }).unwrap();
            let mut top = expected.as_bytes().to_vec(); top.splice(at..at, *b"AC");
            assert_eq!(s.sequence(), top);
            let bytes = s.save_bytes().unwrap();
            let loaded = EditSession::new(DnaFile::parse(&bytes).unwrap(), bytes).unwrap();
            assert_eq!(loaded.sequence(), s.sequence());
            assert_eq!(loaded.end_properties(), s.end_properties());
        }
    }

    #[test]
    fn terminal_metadata_does_not_trim_circular_or_single_stranded_molecules() {
        for (kind, flags) in [(MoleculeKind::Dna, 3), (MoleculeKind::Dna, 0), (MoleculeKind::Rna, 0), (MoleculeKind::Protein, 0)] {
            let mut file = DnaFile::parse(&hx_formats::write::NewFile::of(kind, b"ACGTACGT".to_vec())
                .with_flags(MoleculeFlags::from_byte(flags)).to_bytes()).unwrap();
            file.blocks.iter_mut().find(|b| b.block_type == 8).unwrap().payload =
                b"<AdditionalSequenceProperties><UpstreamStickiness>-2</UpstreamStickiness><DownstreamStickiness>2</DownstreamStickiness></AdditionalSequenceProperties>".to_vec();
            let bytes = file.to_bytes();
            let s = EditSession::new(file, bytes.clone()).unwrap();
            assert_eq!(s.sequence(), b"ACGTACGT");
            assert_eq!(s.save_bytes().unwrap(), bytes);
        }
    }

    #[test]
    fn replacing_a_clipped_primer_cannot_silently_discard_hidden_binding_bases() {
        use hx_formats::primer_edit::PrimerSpec;
        let s = sticky_sample();
        let mut file = s.file().clone();
        file.blocks.iter_mut().find(|b| b.block_type == 5).unwrap().payload =
            b"<Primers nextValidID=\"1\"><Primer name=\"p\" sequence=\"TGGG\"><BindingSite location=\"0-3\" boundStrand=\"0\" annealedBases=\"TGGG\"/></Primer></Primers>".to_vec();
        let bytes = file.to_bytes();
        let mut s = EditSession::new(file, bytes.clone()).unwrap();
        assert_eq!(s.primers()[0].sites[0].range, hx_formats::annot::Range::new(0, 3));
        assert!(s.apply_primer(PrimerEdit::Replace { name: "p".into(), spec: Box::new(PrimerSpec::new("p", "TGGG")) }).is_err());
        assert_eq!(s.save_bytes().unwrap(), bytes);
        s.apply_primer(PrimerEdit::Rename { name: "p".into(), new_name: "renamed".into() }).unwrap();
        assert_eq!(s.stored_primers()[0].sites[0].range, hx_formats::annot::Range::new(0, 4));
    }

    /// 用 `NewFile` 造骨架，再把 block 10/5 的载荷替换为受控 XML。
    fn sample() -> EditSession {
        // 22 bp 环状质粒；特征 3..9（1-based 闭 = 4-9）；引物结合位点 10..19（0-based 闭）
        let seq = b"ACGTACGTACGGGGCCCCTTTT".to_vec();
        let bytes0 = hx_formats::write::NewFile::dna(seq).to_bytes();
        let mut f = DnaFile::parse(&bytes0).unwrap();
        let feats = "<?xml version=\"1.0\"?><Features nextValidID=\"1\"><Feature recentID=\"0\" name=\"f1\" directionality=\"1\" type=\"CDS\"><Segment range=\"3-8\" color=\"#ff0000\" type=\"standard\" translated=\"0\"/></Feature></Features>";
        let primers = "<?xml version=\"1.0\"?><Primers nextValidID=\"1\"><Primer name=\"p1\" sequence=\"AAAA\"><BindingSite simplified=\"1\" location=\"10-19\" boundStrand=\"0\" annealedBases=\"10\"/></Primer></Primers>";
        for b in f.blocks.iter_mut() {
            match b.block_type {
                10 => b.payload = feats.as_bytes().to_vec(),
                5 => b.payload = primers.as_bytes().to_vec(),
                _ => {}
            }
        }
        let bytes = f.to_bytes();
        let parsed = DnaFile::parse(&bytes).unwrap();
        EditSession::new(parsed, bytes).unwrap()
    }

    #[test]
    fn new_session_matches_base() {
        let s = sample();
        assert_eq!(s.len(), 22);
        assert_eq!(s.features().len(), 1);
        assert_eq!(s.features()[0].segments[0].range, hx_formats::annot::Range::new(2, 8));
        assert_eq!(s.primers().len(), 1);
        assert_eq!(s.primers()[0].sites[0].range, hx_formats::annot::Range::new(10, 20));
        assert!(!s.is_dirty());
    }

    #[test]
    fn unedited_save_is_byte_identical() {
        let s = sample();
        assert_eq!(s.save_bytes().unwrap(), s.bytes());
    }

    #[test]
    fn delete_shifts_features_and_primers() {
        let mut s = sample();
        s.apply(Edit::Delete { start: 0, end: 2 }).unwrap();
        assert_eq!(s.sequence(), b"GTACGTACGGGGCCCCTTTT");
        assert_eq!(s.features()[0].segments[0].range, hx_formats::annot::Range::new(0, 6));
        assert_eq!(s.primers()[0].sites[0].range, hx_formats::annot::Range::new(8, 18));
        // 保存后重新解析，特征/引物位置一致
        let bytes = s.save_bytes().unwrap();
        let f = DnaFile::parse(&bytes).unwrap();
        let f_ann = parse_features(&f.block(10).unwrap().payload);
        assert_eq!(f_ann[0].segments[0].range, hx_formats::annot::Range::new(0, 6));
        let p_ann = parse_primers(&f.block(5).unwrap().payload);
        assert_eq!(p_ann[0].sites[0].range, hx_formats::annot::Range::new(8, 18));
        // 序列块也已重建
        assert_eq!(f.sequence().unwrap().sequence, s.sequence());
    }

    /// 回归：结合位点从 0 开始（`location="0-4"`，0-based 闭区间）的引物
    /// 在编辑后必须还在——早期实现 `0 - 1` 下溢，会把整个 `<Primer>` 丢掉。
    #[test]
    fn primer_starting_at_zero_survives_edits() {
        let seq = b"ACGTACGTACGGGGCCCCTTTT".to_vec();
        let bytes0 = hx_formats::write::NewFile::dna(seq).to_bytes();
        let mut f = DnaFile::parse(&bytes0).unwrap();
        let primers = "<?xml version=\"1.0\"?><Primers nextValidID=\"2\">\
            <Primer recentID=\"0\" name=\"p1\" sequence=\"AAAA\">\
            <BindingSite location=\"0-4\" boundStrand=\"0\" annealedBases=\"ACGTA\">\
            <Component hybridizedRange=\"0-4\" bases=\"ACGTA\"/></BindingSite>\
            <BindingSite location=\"10-19\" boundStrand=\"0\" annealedBases=\"GGGGCCCCTT\">\
            <Component hybridizedRange=\"10-19\" bases=\"GGGGCCCCTT\"/></BindingSite>\
            <BindingSite simplified=\"1\" location=\"0-4\" boundStrand=\"0\" annealedBases=\"ACGTA\"/>\
            <BindingSite simplified=\"1\" location=\"10-19\" boundStrand=\"0\" annealedBases=\"GGGGCCCCTT\"/></Primer>\
            <Primer recentID=\"1\" name=\"p2\" sequence=\"CCCC\">\
            <BindingSite location=\"12-17\" boundStrand=\"0\" annealedBases=\"GGCCCC\">\
            <Component hybridizedRange=\"12-17\" bases=\"GGCCCC\"/></BindingSite>\
            <BindingSite simplified=\"1\" location=\"12-17\" boundStrand=\"0\" annealedBases=\"GGCCCC\"/></Primer>\
            </Primers>";
        for b in f.blocks.iter_mut() {
            if b.block_type == 5 {
                b.payload = primers.as_bytes().to_vec();
            }
        }
        let bytes = f.to_bytes();
        let parsed = DnaFile::parse(&bytes).unwrap();
        let mut s = EditSession::new(parsed, bytes).unwrap();
        assert_eq!(s.primers().len(), 2);

        // 删掉第 2 个碱基（0-based 1..2）
        s.apply(Edit::Delete { start: 1, end: 2 }).unwrap();
        assert_eq!(s.primers().len(), 2, "引物被整体丢掉了");
        let p1 = &s.primers()[0];
        assert_eq!(p1.name, "p1");
        assert_eq!(p1.sites.len(), 2);
        assert_eq!(p1.sites[0].range, hx_formats::annot::Range::new(0, 4));
        assert_eq!(p1.sites[1].range, hx_formats::annot::Range::new(9, 19));

        // 存盘后再解析，仍然是两个引物、两个位点
        let out = s.save_bytes().unwrap();
        let f2 = DnaFile::parse(&out).unwrap();
        let p_ann = parse_primers(&f2.block(5).unwrap().payload);
        assert_eq!(p_ann.len(), 2);
        assert_eq!(p_ann[0].sites.len(), 2);
        assert_eq!(p_ann[0].sites[0].range, hx_formats::annot::Range::new(0, 4));
        assert_eq!(p_ann[1].sites[0].range, hx_formats::annot::Range::new(11, 17));
    }

    #[test]
    fn add_primer_when_block5_missing() {
        use hx_formats::primer_edit::{PrimerEdit, PrimerSiteSpec, PrimerSpec};
        // 造一个没有块 5 的最小文件
        let seq = b"ACGTACGTACGGGGCCCCTTTT".to_vec();
        let bytes0 = hx_formats::write::NewFile::dna(seq).to_bytes();
        let mut f = DnaFile::parse(&bytes0).unwrap();
        f.blocks.retain(|b| b.block_type != 5);
        let bytes = f.to_bytes();
        let mut s = EditSession::new(DnaFile::parse(&bytes).unwrap(), bytes).unwrap();
        assert!(s.primers().is_empty());
        // 未编辑：逐字节原样
        assert!(!s.is_dirty());
        s.apply_primer(PrimerEdit::Add(Box::new(PrimerSpec {
            name: "P".into(),
            sequence: "ACGTACGTAC".into(),
            description: String::new(),
            date_added: String::new(),
            phosphorylated: false,
            sites: vec![PrimerSiteSpec {
                range: hx_formats::annot::Range::new(0, 10),
                strand: 0,
                annealed: "ACGTACGTAC".into(),
                tm: 40.0,
            }],
        })))
        .unwrap();
        assert_eq!(s.primers().len(), 1);
        let out = DnaFile::parse(&s.save_bytes().unwrap()).unwrap();
        let block5 = out.block(5).expect("block 5 inserted");
        assert_eq!(parse_primers(&block5.payload).len(), 1);
        // 块 5 落在块 8 之前（写回时按块号插在第一个大于 5 的块前面）
        let types: Vec<u8> = out.blocks.iter().map(|b| b.block_type).collect();
        let i5 = types.iter().position(|t| *t == 5).unwrap();
        let i8 = types.iter().position(|t| *t == 8).unwrap();
        assert!(i5 < i8, "{types:?}");
    }

    #[test]
    fn add_and_remove_primer_writes_block5() {
        use hx_formats::primer_edit::{PrimerEdit, PrimerSiteSpec, PrimerSpec};
        let mut s = sample();
        let spec = PrimerSpec {
            name: "new".into(),
            sequence: "AAAACCCCGGGGTTTT".into(),
            description: String::new(),
            date_added: String::new(),
            phosphorylated: false,
            sites: vec![PrimerSiteSpec {
                range: hx_formats::annot::Range::new(0, 16),
                strand: 0,
                annealed: "AAAACCCCGGGGTTTT".into(),
                tm: 55.0,
            }],
        };
        s.apply_primer(PrimerEdit::Add(Box::new(spec))).unwrap();
        assert_eq!(s.primers().len(), 2);
        assert_eq!(s.primer_edit_count(), 1);
        // 纯引物编辑不写历史树
        let (bytes, nodes) = s
            .save_bytes_with_history("x.dna", "x.dna", None)
            .unwrap();
        assert_eq!(nodes, 0);
        let f = DnaFile::parse(&bytes).unwrap();
        let p = parse_primers(&f.block(5).unwrap().payload);
        assert_eq!(p.len(), 2);
        assert_eq!(p[1].name, "new");
        assert_eq!(p[1].sites[0].range, hx_formats::annot::Range::new(0, 16));
        // undo 回到一个引物
        s.undo().unwrap();
        assert_eq!(s.primers().len(), 1);
        s.redo().unwrap();
        assert_eq!(s.primers().len(), 2);
        // 删除后落盘只剩原引物
        s.apply_primer(PrimerEdit::Remove { name: "new".into() }).unwrap();
        assert_eq!(s.primers().len(), 1);
        let f = DnaFile::parse(&s.save_bytes().unwrap()).unwrap();
        assert_eq!(parse_primers(&f.block(5).unwrap().payload).len(), 1);
    }

    #[test]
    fn deleting_whole_feature_removes_it() {
        let mut s = sample();
        s.apply(Edit::Delete { start: 2, end: 9 }).unwrap();
        assert!(s.features().is_empty(), "{:?}", s.features());
    }

    #[test]
    fn insert_and_replace() {
        let mut s = sample();
        s.apply(Edit::Insert { at: 4, seq: b"TTTT".to_vec() }).unwrap();
        assert_eq!(s.sequence(), b"ACGTTTTTACGTACGGGGCCCCTTTT".to_vec());
        s.apply(Edit::Replace { start: 0, end: 4, seq: b"GG".to_vec() }).unwrap();
        assert_eq!(&s.sequence()[..6], b"GGTTTT");
        assert_eq!(s.op_count(), 2);
    }

    #[test]
    fn undo_redo_round_trip() {
        let mut s = sample();
        let before = s.sequence().to_vec();
        s.apply(Edit::Delete { start: 5, end: 8 }).unwrap();
        let after = s.sequence().to_vec();
        assert_ne!(before, after);
        s.undo().unwrap();
        assert_eq!(s.sequence(), before.as_slice());
        assert!(s.can_redo());
        s.redo().unwrap();
        assert_eq!(s.sequence(), after.as_slice());
        assert!(!s.can_redo());
    }

    #[test]
    fn undo_across_multiple_edits() {
        let mut s = sample();
        let base = s.sequence().to_vec();
        s.apply(Edit::Delete { start: 0, end: 3 }).unwrap();
        s.apply(Edit::Insert { at: 0, seq: b"CC".to_vec() }).unwrap();
        let two = s.sequence().to_vec();
        s.undo().unwrap();
        assert_ne!(s.sequence(), two.as_slice());
        s.undo().unwrap();
        assert_eq!(s.sequence(), base.as_slice());
        assert!(!s.can_undo());
        assert!(s.undo().is_none());
    }

    #[test]
    fn invalid_edit_is_rejected_and_state_kept() {
        let mut s = sample();
        let before = s.sequence().to_vec();
        assert!(s.apply(Edit::Delete { start: 5, end: 999 }).is_err());
        assert_eq!(s.sequence(), before.as_slice());
        assert!(!s.is_dirty());
    }

    #[test]
    fn save_with_history_plants_block7_after_block6() {
        let mut s = sample();
        s.apply(Edit::Delete { start: 0, end: 2 }).unwrap();
        let (bytes, nodes) = s
            .save_bytes_with_history("saved.dna", "src.dna", None)
            .unwrap();
        assert_eq!(nodes, 2, "新根 + 父节点");
        let f = DnaFile::parse(&bytes).unwrap();
        let types: Vec<u8> = f.blocks.iter().map(|b| b.block_type).collect();
        let i6 = types.iter().position(|t| *t == 6).expect("block 6");
        assert_eq!(types[i6 + 1], 7, "block 7 紧跟 block 6：{types:?}");
        assert!(!types.contains(&3), "block 3 缓存已丢弃：{types:?}");
        let payload = &f.block(7).unwrap().payload;
        let text = std::str::from_utf8(payload).unwrap();
        assert!(text.contains("saved.dna") && text.contains("src.dna"), "{text}");
        // 序列/特征与无历史版本一致
        assert_eq!(f.sequence().unwrap().sequence, s.sequence());
        // 未编辑 ⇒ 不写历史，原样字节
        let fresh = sample();
        let (b2, n2) = fresh.save_bytes_with_history("a.dna", "b.dna", None).unwrap();
        assert_eq!((b2, n2), (fresh.bytes().to_vec(), 0));
    }

    #[test]
    fn commit_saved_resets_history() {
        let mut s = sample();
        s.apply(Edit::Delete { start: 0, end: 2 }).unwrap();
        let bytes = s.save_bytes().unwrap();
        s.commit_saved(bytes.clone()).unwrap();
        assert!(!s.is_dirty());
        assert!(!s.can_undo());
        assert_eq!(s.save_bytes().unwrap(), bytes, "保存后再保存应逐字节一致");
        assert_eq!(s.len(), 20);
    }

    #[test]
    fn seq_text_and_revcomp() {
        assert_eq!(parse_seq_text("acgt ACGT\n123 tt"), b"acgtACGTtt".to_vec());
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());
        assert_eq!(reverse_complement(b"GaTt"), b"aAtC".to_vec());
    }

    // ---- 特征编辑 ----

    #[test]
    fn add_feature_shows_up_and_saves() {
        let mut s = sample();
        let mut sp = hx_formats::feature_edit::FeatSpec::new(10, 15);
        sp.name = "new".into();
        sp.kind = "promoter".into();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Add(sp)).unwrap();
        assert!(s.is_dirty());
        assert_eq!(s.features().len(), 2);
        assert_eq!(s.feature_specs().len(), 2);
        assert_eq!(s.feature_specs()[1].name, "new");
        // 无序列编辑时，序列字节与基准一致，只有块 10 变了
        let b = s.save_bytes().unwrap();
        let f = DnaFile::parse(&b).unwrap();
        let specs = hx_formats::feature_edit::parse_specs(&f.block(10).unwrap().payload);
        assert_eq!(specs.len(), 2);
        assert_eq!(specs[1].segments[0].start, 10);
        assert_eq!(f.sequence().unwrap().sequence, s.sequence());
        // 保存后回到干净状态
        s.commit_saved(b).unwrap();
        assert!(!s.is_dirty());
        assert_eq!(s.features().len(), 2);
    }

    #[test]
    fn edit_and_remove_feature() {
        let mut s = sample();
        let mut sp = s.feature_specs()[0].clone();
        sp.name = "renamed".into();
        sp.note = "a note".into();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Update(sp)).unwrap();
        assert_eq!(s.feature_specs()[0].name, "renamed");
        assert_eq!(s.feature_specs()[0].note, "a note");
        let b = s.save_bytes().unwrap();
        let f = DnaFile::parse(&b).unwrap();
        let x = String::from_utf8_lossy(&f.block(10).unwrap().payload).to_string();
        assert!(x.contains("name=\"renamed\""), "{x}");
        assert!(x.contains("<Q name=\"note\"><V text=\"a note\"/></Q>"), "{x}");

        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Remove { recent_id: 0 })
            .unwrap();
        assert!(s.features().is_empty());
        // undo 两步回到原样
        s.undo().unwrap();
        s.undo().unwrap();
        assert_eq!(s.features().len(), 1);
        assert!(!s.is_dirty());
        assert_eq!(s.save_bytes().unwrap(), s.bytes());
    }

    #[test]
    fn feature_edit_after_sequence_edit_uses_current_coords() {
        let mut s = sample();
        // 删掉开头 2 bp：特征 3-8（1-based）→ 1-6
        s.apply(Edit::Delete { start: 0, end: 2 }).unwrap();
        assert_eq!(s.feature_specs()[0].segments[0].start, 1);
        assert_eq!(s.feature_specs()[0].segments[0].end, 6);
        // 再改名为 CDS2 并新加一个特征（用当前坐标）
        let mut sp = s.feature_specs()[0].clone();
        sp.kind = "gene".into();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Update(sp)).unwrap();
        let mut ns = hx_formats::feature_edit::FeatSpec::new(1, 3);
        ns.name = "tail".into();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Add(ns)).unwrap();
        let b = s.save_bytes().unwrap();
        let f = DnaFile::parse(&b).unwrap();
        let specs = hx_formats::feature_edit::parse_specs(&f.block(10).unwrap().payload);
        assert_eq!(specs.len(), 2);
        assert_eq!(specs[0].kind, "gene");
        assert_eq!(specs[0].segments[0].start, 1);
        assert_eq!(specs[1].name, "tail");
        assert_eq!(specs[1].segments[0].end, 3);
        // 序列编辑后的纯特征编辑**不**伪造历史树
        let (hb, nodes) = s.save_bytes_with_history("a.dna", "b.dna", None).unwrap();
        assert_eq!(nodes, 2);
        assert!(DnaFile::parse(&hb).unwrap().block(7).is_some());
    }

    #[test]
    fn pure_feature_edit_drops_history_and_keeps_sequence() {
        let mut s = sample();
        let mut sp = s.feature_specs()[0].clone();
        sp.name = "x".into();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::Update(sp)).unwrap();
        let (b, nodes) = s.save_bytes_with_history("a.dna", "b.dna", None).unwrap();
        assert_eq!(nodes, 0);
        let f = DnaFile::parse(&b).unwrap();
        assert!(f.block(7).is_none());
        // 块集合（除块 6 的时间戳外）保持不变
        let before: Vec<u8> = s.file().blocks.iter().map(|x| x.block_type).collect();
        let after: Vec<u8> = f.blocks.iter().map(|x| x.block_type).collect();
        assert_eq!(before, after, "纯特征编辑不该增删块");
        assert_eq!(f.sequence().unwrap().sequence, s.sequence());
    }

    #[test]
    fn annotation_save_preserves_existing_lineage() {
        let mut s = sample();
        s.apply(Edit::Insert { at: 2, seq: b"AA".to_vec() }).unwrap();
        let (bytes, _) = s.save_bytes_with_history("x.dna", "original.dna", None).unwrap();
        s.commit_saved(bytes).unwrap();
        let history = s.file().block(7).unwrap().payload.clone();
        let mut sp = s.feature_specs()[0].clone();
        sp.name = "renamed".into();
        s.apply_feature(FeatureEdit::Update(sp)).unwrap();
        s.set_note("Description", "new notes").unwrap();
        let (bytes, count) = s.save_bytes_with_history("x.dna", "x.dna", None).unwrap();
        assert_eq!(count, 0);
        let saved = DnaFile::parse(&bytes).unwrap();
        assert_eq!(saved.block(7).unwrap().payload, history);
        assert_eq!(saved.blocks_of(7).count(), 1);
    }

    #[test]
    fn each_sequence_operation_tracks_undo_redo_and_save() {
        let mut s = sample();
        s.apply(Edit::Insert { at: 2, seq: b"AA".to_vec() }).unwrap();
        s.apply(Edit::Delete { start: 5, end: 7 }).unwrap();
        assert_eq!(s.save_bytes_with_history("x.dna", "src.dna", None).unwrap().1, 3);
        s.undo().unwrap();
        assert_eq!(s.save_bytes_with_history("x.dna", "src.dna", None).unwrap().1, 2);
        s.redo().unwrap();
        let (bytes, count) = s.save_bytes_with_history("x.dna", "src.dna", None).unwrap();
        assert_eq!(count, 3);
        s.commit_saved(bytes.clone()).unwrap();
        assert_eq!(s.save_bytes_with_history("x.dna", "x.dna", None).unwrap().0, bytes);
        s.apply(Edit::Insert { at: 0, seq: b"A".to_vec() }).unwrap();
        assert_eq!(s.save_bytes_with_history("x.dna", "x.dna", None).unwrap().1, 4);
    }

    #[test]
    fn hybridization_recalculation_is_atomic_undoable_and_roundtrips(){
        let bytes=hx_formats::write::NewFile::dna(b"AACAGGTT".to_vec()).to_bytes();
        let mut s=EditSession::new(DnaFile::parse(&bytes).unwrap(),bytes).unwrap();
        s.apply_primer(PrimerEdit::Add(Box::new(hx_formats::primer_edit::PrimerSpec::new("p","AACCGGTT")))).unwrap();
        let before=s.save_bytes().unwrap();
        let params=hx_formats::annot::HybridizationParams{min_continuous_match_len:8,allow_mismatch:true,min_melting_temperature:-100.0,show_additional_five_prime_matches:false,..Default::default()};
        s.recalculate_hybridization(params).unwrap();assert!(!s.primers()[0].sites.is_empty());
        let bytes=s.save_bytes_with_history("x","x",None).unwrap().0;
        let reopened=EditSession::new(DnaFile::parse(&bytes).unwrap(),bytes).unwrap();assert_eq!(reopened.primers(),s.primers());
        s.undo();assert_eq!(s.save_bytes().unwrap(),before);s.redo();assert!(!s.primers()[0].sites.is_empty());
    }

    #[test]
    fn hide_feature_is_reversible_byte_identically() {
        let mut s = sample();
        let before = s.bytes().to_vec();
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::SetVisible {
            recent_id: 0,
            visible: false,
        })
        .unwrap();
        assert!(!s.feature_specs()[0].visible);
        let hidden = s.save_bytes().unwrap();
        assert!(String::from_utf8_lossy(&hidden).contains("visible=\"0\""));
        s.apply_feature(hx_formats::feature_edit::FeatureEdit::SetVisible {
            recent_id: 0,
            visible: true,
        })
        .unwrap();
        assert_eq!(s.save_bytes().unwrap(), before);
    }
}
