//! 打开的文件：编辑会话（`hx-algo`）+ 酶位点缓存 + 展示信息。

use std::path::{Path, PathBuf};

use hx_algo::EditSession;
use hx_core::MoleculeKind;
use hx_enzymes::ENZYMES;
use hx_formats::annot::{enzyme_set_name, FeatureAnnotation, PrimerAnnotation};
use hx_formats::feature_edit::{FeatureEdit, FeatSpec};
use hx_formats::primer_edit::PrimerEdit;
use hx_formats::write::Edit;
use hx_formats::DnaFile;

/// 只能**原样保存**、不能编辑的原生 SnapGene 文档（DNA/蛋白比对、凝胶…）。
///
/// 这类文件用同一套 TLV 容器（块结构已解，能逐字节往返），但内容属于 SnapGene 的
/// 比对编辑器 / 凝胶编辑器，我们没有对应视图。做法：保留原始字节，
/// 视图页签只留一个「概览」页，画布上说明情况，保存时原样写回。
#[derive(Debug, Clone)]
pub struct NativeOnly {
    pub dirty: bool,
    saved_bytes: Vec<u8>,
    undo: Vec<Vec<u8>>,
    redo: Vec<Vec<u8>>,
    /// 类别名（`DNA alignment` / `Gel`…）。
    pub class: &'static str,
    /// 头部 `type_of_sequence`。
    pub type_of_sequence: u16,
    /// 原始字节（保存时原样写出）。
    pub bytes: Vec<u8>,
    /// 块摘要（`类型·字节数`），用于概览面板。
    pub blocks: Vec<(u8, usize)>,
}

/// 一个已打开的文档（含未保存编辑）。
pub struct Document {
    pub path: Option<PathBuf>,
    /// 编辑会话：基准文件 + 工作副本 + undo/redo。
    pub session: EditSession,
    /// block 13 里记录的酶显示集合名（若有）。
    pub enzyme_set_name: Option<String>,
    /// 未保存文档的显示名（导入的 FASTA/GenBank 记录名）。
    pub name_hint: Option<String>,
    /// 只能原样保存的原生文档（比对/凝胶）；`Some` 时 `session` 是空占位。
    pub native_only: Option<NativeOnly>,
    /// `Secondary Structure` 页的解析/布局缓存（按块长度 + 编辑计数失效）。
    pub rna_cache: Option<RnaStructCache>,
    /// 本地折叠任务：没有 block 34（或序列被改过）时算结构 + 次优结构。
    /// 文件里已经存了 block 34 时不再重算 —— SnapGene 的着色**就是**按块里的
    /// 矩阵画的（见 `rna_structure`），重算会让同一个文件两边颜色对不上。
    fold_job: Option<FoldJob>,
}

/// block 34（XZ 压缩 JSON）的解析结果 + 布局缓存。
///
/// 16S 的块解压 + 布局都在毫秒级，但 `Secondary Structure` 页每帧都要用，
/// 所以按文档缓存；用户改动序列后编辑计数变化，缓存自动重算。
pub struct RnaStructCache {
    /// `(block 34 载荷长度, 编辑计数)`。
    key: (usize, usize),
    /// 主序列（原样；显示成 `U` 的映射在绘制层做）。
    seq: Vec<u8>,
    pub structure: Result<hx_formats::rna_struct::RnaStructure, String>,
    /// 界面上的结构列表（MFE + 次优），与 `RnaStructure::structures()` 同序。
    pub structures: Option<Vec<(Vec<(usize, usize)>, f64)>>,
    /// 当前画出来的那个结构：`(在 `structures()` 里的下标, 布局)`。
    pub layout: Option<(usize, hx_render::RnaLayout)>,
    /// 次优结构列表的缩略图布局（与 `structures()` 同序）。
    pub thumbs: Option<Vec<hx_render::RnaLayout>>,
    /// 每个结构的 naview 坐标（内嵌 ViennaRNA 的 `RNAplot` 布局）。
    /// `None` = 还没算；`Some(vec![])` = 布局失败，回退到内置布局。
    naview: Option<Vec<Vec<(f32, f32)>>>,
}

impl RnaStructCache {
    /// 结构列表（MFE 在前，其后为次优结构）。
    pub fn ensure_structures(&mut self, _len: usize) {
        if self.structures.is_some() {
            return;
        }
        self.structures = Some(match &self.structure {
            Ok(st) => st.structures(),
            Err(_) => Vec::new(),
        });
    }

    /// 取（必要时用内嵌 ViennaRNA 算）全部结构的 naview 坐标。
    fn ensure_naview(&mut self) {
        if self.naview.is_some() {
            return;
        }
        let Some(list) = self.structures.as_ref() else {
            self.naview = Some(Vec::new());
            return;
        };
        let len = self.seq.len();
        let dbs: Vec<String> = match &self.structure {
            Ok(st) => list
                .iter()
                .map(|(pairs, _)| st.dot_bracket(pairs, len))
                .collect(),
            Err(_) => Vec::new(),
        };
        self.naview = Some(
            dbs.iter()
                .map(|db| hx_vienna::naview(db).unwrap_or_default())
                .collect(),
        );
    }

    /// 取得（必要时计算）某个结构在 100% 坐标下的布局。
    pub fn ensure_layout(&mut self, idx: usize, len: usize) {
        if self.layout.as_ref().is_some_and(|(i, _)| *i == idx) {
            return;
        }
        self.ensure_structures(len);
        self.ensure_naview();
        if let Some(pos) = self.naview.as_ref().and_then(|all| all.get(idx)) {
            if pos.len() == len {
                self.layout = Some((idx, hx_render::RnaLayout::from_positions(len, pos)));
                return;
            }
        }
        let Some((pairs, _)) = self
            .structures
            .as_ref()
            .and_then(|v| v.get(idx))
            .map(|(p, e)| (p.clone(), *e))
        else {
            return;
        };
        self.layout = Some((idx, hx_render::RnaLayout::new(len, &pairs)));
    }

    /// 取得（必要时计算）全部次优结构的缩略图布局。
    pub fn ensure_thumbs(&mut self, len: usize) {
        if self.thumbs.is_some() {
            return;
        }
        self.ensure_structures(len);
        self.ensure_naview();
        let naview = self.naview.clone().unwrap_or_default();
        let list = self.structures.as_ref().map(|v| {
            v.iter()
                .enumerate()
                .map(|(i, (pairs, _))| match naview.get(i) {
                    Some(pos) if pos.len() == len => {
                        hx_render::RnaLayout::from_positions(len, pos)
                    }
                    _ => hx_render::RnaLayout::new(len, pairs),
                })
                .collect()
        });
        self.thumbs = Some(list.unwrap_or_default());
    }
}

/// `Secondary Structure` 页拿到的缓存状态。
pub enum RnaCache<'a> {
    /// 已就绪（`Ok` 来自 block 34，或本地折叠的结果）。
    Ready(&'a mut RnaStructCache),
    /// 本地折叠正在后台线程里跑。
    Folding,
    /// 没有可显示的结构（空序列）。
    None,
}

/// 一个本地折叠任务（没有 block 34，或序列被改过时用）。
///
/// 16S 这种 1500 nt 的序列，MFE + 配分函数在 release 下要 8 秒左右，
/// 绝不能放在 UI 线程上，所以一律后台算，界面先显示 `Folding…`。
/// 次优结构枚举的能量带宽（单位 10 cal/mol，200 = MFE + 2.0 kcal/mol）。
///
/// SnapGene 的列表规律（tRNA 72 nt 与 445566 365 nt 两个样本实测一致）：
/// 第一条是 MFE，其后最多 9 条次优，能量下界到 MFE + ~2.0 kcal/mol ——
/// 正是 `vrna_subopt(delta = 200)` 的前 9 条。
const SUBOPT_DELTA: i32 = 200;
/// 列表里保留的次优结构条数（不含 MFE）。
const SUBOPT_MAX: usize = 9;
/// 超过这个长度就不再枚举次优：Wuchty 枚举是结构数级别的开销
/// （1541 nt 的 16S 在 delta=1.0 时已经要 67 s / 57 万条），不值得为侧栏列表
/// 卡住后台线程。这类文件只显示 MFE；SnapGene 打开后会自己补算。
const SUBOPT_MAX_LEN: usize = 600;

/// 折叠结果 + 次优结构（配对表, 能量）。
type FoldResult = (hx_vienna::Folded, Vec<(Vec<(usize, usize)>, f64)>);

struct FoldJob {
    /// 与 `RnaStructCache::key` 同一套失效键。
    key: (usize, usize),
    seq: Vec<u8>,
    handle: std::thread::JoinHandle<Result<FoldResult, String>>,
}

impl FoldJob {
    fn spawn(key: (usize, usize), seq: Vec<u8>) -> Self {
        let fold_seq: String = seq
            .iter()
            .map(|&b| {
                let c = b.to_ascii_uppercase();
                if c == b'T' {
                    'U'
                } else {
                    c as char
                }
            })
            .collect();
        let handle = std::thread::spawn(move || {
            let folded = hx_vienna::fold(&fold_seq, 37.0, true, false)?;
            let subopt = if fold_seq.len() <= SUBOPT_MAX_LEN {
                // Wuchty 的枚举会把 MFE 自己也算一条（SnapGene 的列表不含它），
                // 所以先剔掉与 MFE 配对表相同的条目。
                let mfe = folded.pairs();
                hx_vienna::subopt(&fold_seq, 37.0, true, false, SUBOPT_DELTA)
                    .unwrap_or_default()
                    .into_iter()
                    .filter(|(pairs, _)| *pairs != mfe)
                    .take(SUBOPT_MAX)
                    .collect()
            } else {
                Vec::new()
            };
            Ok((folded, subopt))
        });
        Self { key, seq, handle }
    }
}

/// 把 ViennaRNA 的折叠结果包装成与 block 34 同构的结构。
fn structure_from_fold((f, subopt): FoldResult) -> hx_formats::rna_struct::RnaStructure {
    const R: f64 = 0.0019872041;
    let rt = R * (37.0 + 273.15);
    let freq = (-(f.mfe - f.ensemble) / rt).exp();
    let suboptimal = subopt
        .into_iter()
        .map(|(pairs, energy)| hx_formats::rna_struct::Suboptimal { pairs, energy })
        .collect();
    hx_formats::rna_struct::RnaStructure::from_fold(
        f.pairs(),
        f.probs.clone(),
        suboptimal,
        hx_formats::rna_struct::StructStats {
            free_energy_ensemble: f.ensemble,
            free_energy_optimal: f.mfe,
            freq_optimal: freq,
            ensemble_diversity: f.diversity,
        },
        37.0,
        true,
        false,
    )
}

impl Document {
    pub fn from_path(path: &Path) -> Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        let mut doc = Self::from_bytes(bytes)?;
        doc.path = Some(path.to_path_buf());
        Ok(doc)
    }

    /// 由导入结果建立未命名文档（FASTA/GenBank；`hint` = 记录名）。
    pub fn from_imported(bytes: Vec<u8>, hint: String) -> Result<Self, String> {
        let mut doc = Self::from_bytes(bytes)?;
        doc.name_hint = Some(hint);
        Ok(doc)
    }

    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self, String> {
        let file = DnaFile::parse(&bytes).map_err(|e| format!("parse .dna: {e}"))?;
        // 比对 / 凝胶等：没有可显示的序列，保留字节 + 类别名
        if MoleculeKind::from_type_of_sequence(file.header.type_of_sequence).is_none() {
            let type_of_sequence = file.header.type_of_sequence;
            let class = hx_core::document_class_name(type_of_sequence);
            let blocks: Vec<(u8, usize)> = file
                .blocks
                .iter()
                .map(|b| (b.block_type, b.payload.len()))
                .collect();
            let empty = hx_formats::write::NewFile::dna(Vec::new()).to_bytes();
            let empty_file = DnaFile::parse(&empty).map_err(|e| format!("parse .dna: {e}"))?;
            let session = EditSession::new(empty_file, empty)?;
            return Ok(Self {
                path: None,
                session,
                enzyme_set_name: None,
                name_hint: None,
                native_only: Some(NativeOnly {
                    dirty: false,
                    saved_bytes: bytes.clone(),
                    undo: Vec::new(),
                    redo: Vec::new(),
                    class,
                    type_of_sequence,
                    bytes,
                    blocks,
                }),
                rna_cache: None,
                fold_job: None,
            });
        }
        let enzyme_set_name = file.block(13).and_then(|b| enzyme_set_name(&b.payload));
        let session = EditSession::new(file, bytes)?;
        Ok(Self {
            path: None,
            session,
            enzyme_set_name,
            name_hint: None,
            native_only: None,
            rna_cache: None,
            fold_job: None,
        })
    }

    // ---- 转发给编辑会话 ----

    pub fn len(&self) -> usize {
        self.session.len()
    }

    pub fn sequence(&self) -> &[u8] {
        self.session.sequence()
    }

    pub fn features(&self) -> &[FeatureAnnotation] {
        self.session.features()
    }

    pub fn primers(&self) -> &[PrimerAnnotation] {
        self.session.primers()
    }

    pub fn circular(&self) -> bool {
        self.session.circular()
    }

    pub fn block_count(&self) -> usize {
        self.session.block_count()
    }

    /// `block 8` 的末端属性（粘性末端 / 修饰 / 原始首碱基）。
    pub fn end_properties(&self) -> hx_formats::dna::EndProperties {
        self.session.end_properties().clone()
    }

    /// 克隆入口始终使用原生存储并集；主编辑器公开的是顶链坐标。
    pub fn cloning_source(&self, label: impl Into<String>) -> hx_algo::cloning::Source {
        hx_algo::cloning::Source::new(label, self.session.stored_sequence().to_vec(), self.circular())
            .with_features(self.session.stored_features().to_vec())
            .with_primers(self.session.stored_primers().to_vec())
            .with_ends(self.end_properties())
    }

    /// `block 6` 里用户自定义的图谱标签（商品化载体的显示名）。
    pub fn custom_map_label(&self) -> Option<String> {
        let notes=self.session.notes_xml();
        if hx_formats::notes::field(&notes,"UseCustomMapLabel").as_deref()!=Some("1"){return None;}
        hx_formats::notes::field(&notes,"CustomMapLabel").filter(|s|!s.trim().is_empty())
    }

    /// 某个块类型的载荷（用于展示「文件里存了什么」；编辑后仍以基准文件为准）。
    pub fn block_payload(&self, ty: u8) -> Option<&[u8]> {
        if ty == 5 { return self.session.primers_xml(); }
        self.session.file().block(ty).map(|b| b.payload.as_slice())
    }

    /// `Secondary Structure` 页的数据：文件里有 block 34 就直接读（与 SnapGene
    /// 完全一致）；没有、或者用户改过序列，就用内嵌的 ViennaRNA 本地折叠。
    ///
    /// **块里的 `probabilityMatrix` 就是 SnapGene 屏幕上的着色依据**（用一份
    /// 全 0.05 的假矩阵在 SnapGene 里验证过：配对碱基全变蓝、未配对碱基变红），
    /// 所以这条路径不能拿现算的概率去替换它 —— 否则同一个文件在两边颜色会不同。
    ///
    /// 折叠在后台线程，第一次调用返回 `Folding`，界面轮询到完成后再取结果。
    pub fn rna_structure(&mut self) -> Result<RnaCache<'_>, String> {
        let payload = self.block_payload(34).map(|p| p.to_vec());
        let edits = self.edit_count();
        let key = (payload.as_ref().map_or(0, |p| p.len()), edits);
        let seq_now = self.session.sequence().to_vec();

        // 缓存命中就直接复用（读块的结果与本地折叠的结果都算）。
        // 少了这一步，无 block 34 的文件在折叠完成后 `fold_job` 已经取空，
        // 下一帧又会 `map_or(true, …)` 判成「过期」→ 重新起折叠任务，
        // 界面在 `Folding…` 和结果之间逐帧来回切，看起来就是不断闪烁。
        if self
            .rna_cache
            .as_ref()
            .is_some_and(|c| c.key == key && c.seq == seq_now)
        {
            return Ok(RnaCache::Ready(self.rna_cache.as_mut().expect("cache")));
        }

        // 原样未改的文件：block 34 就是 SnapGene 自己算的那份，直接读。
        if edits == 0 {
            if let Some(p) = payload {
                if self.rna_cache.as_ref().map_or(true, |c| c.key != key) {
                    let structure = hx_formats::rna_struct::RnaStructure::parse_payload(&p)
                        .map_err(|e| e.to_string());
                    self.rna_cache = Some(RnaStructCache {
                        key,
                        seq: self.session.sequence().to_vec(),
                        structure,
                        structures: None,
                        layout: None,
                        thumbs: None,
                        naview: None,
                    });
                }
                return Ok(RnaCache::Ready(self.rna_cache.as_mut().expect("cache")));
            }
        }

        let seq = self.session.sequence().to_vec();
        if seq.is_empty() {
            return Ok(RnaCache::None);
        }

        let stale = self
            .fold_job
            .as_ref()
            .map_or(true, |j| j.key != key || j.seq != seq);
        if stale {
            self.fold_job = Some(FoldJob::spawn(key, seq));
            return Ok(RnaCache::Folding);
        }
        if self.fold_job.as_ref().is_some_and(|j| !j.handle.is_finished()) {
            return Ok(RnaCache::Folding);
        }

        let job = self.fold_job.take().expect("fold job");
        let folded = job
            .handle
            .join()
            .map_err(|_| "本地折叠线程异常退出".to_string())?;
        let structure = folded.map(structure_from_fold);
        self.rna_cache = Some(RnaStructCache {
            key,
            seq: job.seq,
            structure,
            structures: None,
            layout: None,
            thumbs: None,
            naview: None,
        });
        Ok(RnaCache::Ready(self.rna_cache.as_mut().expect("cache")))
    }

    pub fn is_dirty(&self) -> bool {
        self.session.is_dirty() || self.native_only.as_ref().is_some_and(|n|n.dirty)
    }

    pub fn can_undo(&self) -> bool {
        self.native_only.as_ref().map_or_else(||self.session.can_undo(),|n|!n.undo.is_empty())
    }

    pub fn can_redo(&self) -> bool {
        self.native_only.as_ref().map_or_else(||self.session.can_redo(),|n|!n.redo.is_empty())
    }

    /// Commit one accepted native edit, shared by gel controls and history routing.
    pub fn edit_native(&mut self, bytes: Vec<u8>) -> Result<(),String> {
        let file=DnaFile::parse(&bytes).map_err(|e|e.to_string())?;
        let n=self.native_only.as_mut().ok_or("Not a native document")?;
        if file.header.type_of_sequence!=n.type_of_sequence{return Err("Native document type changed".into());}
        if bytes==n.bytes{return Ok(());}
        n.undo.push(std::mem::replace(&mut n.bytes,bytes));n.redo.clear();
        n.blocks=file.blocks.iter().map(|b|(b.block_type,b.payload.len())).collect();n.dirty=n.bytes!=n.saved_bytes;Ok(())
    }

    pub fn edit_count(&self) -> usize {
        self.session.op_count()
    }

    // ---- 特征编辑 ----

    /// 工作块 10 的可编辑模型（对话框用）。
    pub fn feature_specs(&self) -> &[FeatSpec] {
        self.session.feature_specs()
    }

    /// 某个特征（按 `recentID`）的当前模型。
    pub fn feature_spec(&self, recent_id: u32) -> Option<&FeatSpec> {
        self.session.feature_spec(recent_id)
    }

    /// 施加一条特征编辑。
    pub fn apply_feature(&mut self, edit: FeatureEdit) -> Result<(), String> {
        self.session.apply_feature(edit)
    }

    /// 施加一条引物编辑。
    pub fn apply_primer(&mut self, edit: PrimerEdit) -> Result<(), String> {
        self.session.apply_primer(edit)
    }

    /// 是否为 DNA 分子（酶切显示只对 DNA 有意义）。
    pub fn is_dna(&self) -> bool {
        self.session.kind() == MoleculeKind::Dna
    }

    /// 分子类型（DNA / RNA / 蛋白）。
    pub fn kind(&self) -> MoleculeKind {
        self.session.kind()
    }

    pub fn is_protein(&self) -> bool {
        self.kind() == MoleculeKind::Protein
    }

    pub fn is_rna(&self) -> bool {
        self.kind() == MoleculeKind::Rna
    }

    /// 是否为单链分子（RNA 与单链 DNA 都按单链渲染；蛋白等价于单链）。
    pub fn is_single_stranded(&self) -> bool {
        !self.is_dna() || self.session.flags().strandedness() == hx_core::Strandedness::Single
    }

    /// 长度单位（`bp` / `bases` / `aa`）。
    pub fn unit(&self) -> &'static str {
        self.kind().unit()
    }

    /// 画布右上角的摘要，如 `pET160 GW D-TOPO.dna (5835 bp)`、`THYX_MYCTU.prot (250 aa)`。
    pub fn summary(&self) -> String {
        if let Some(class) = self.native_class() {
            return format!("{} ({class})", self.display_name());
        }
        format!(
            "{} ({} {})",
            self.path
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| self.display_name()),
            crate::app::nfmt_plain(self.len()),
            self.unit()
        )
    }

    /// 标题栏用（SnapGene：`名字 (Topology / N unit)`；蛋白没有拓扑）。
    pub fn title_summary(&self) -> String {
        if self.kind().is_nucleotide() {
            format!(
                "{} ({} / {} {})",
                self.display_name(),
                if self.circular() { "Circular" } else { "Linear" },
                crate::app::nfmt_plain(self.len()),
                self.unit()
            )
        } else {
            format!(
                "{} ({} {})",
                self.display_name(),
                crate::app::nfmt_plain(self.len()),
                self.unit()
            )
        }
    }

    pub fn topo_label(&self) -> &'static str {
        if self.circular() {
            "circular"
        } else {
            "linear"
        }
    }

    pub fn display_name(&self) -> String {
        match &self.path {
            Some(p) => p
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("untitled")
                .to_string(),
            // 导入的 FASTA/GenBank 记录名（尚未保存成 .dna）
            None => self
                .name_hint
                .clone()
                .unwrap_or_else(|| "untitled".to_string()),
        }
    }

    /// Native map title is independent of the filename (Save As must not
    /// silently replace an explicitly enabled CustomMapLabel).
    pub fn map_name(&self) -> String {
        self.custom_map_label().unwrap_or_else(||{
            if self.path.is_none() {
                if let Some(name)=&self.name_hint {
                    let path=Path::new(name);
                    if path.extension().is_some_and(|e|matches!(e.to_str(),Some("dna"|"rna"|"prot"))) {
                        return path.file_stem().unwrap_or_default().to_string_lossy().into_owned();
                    }
                }
            }
            self.display_name()
        })
    }

    /// 文档页签用的名字：**带扩展名**的文件名（SnapGene 实测页签写作
    /// `16S rRNA.rna` / `THYX_MYCTU.prot`；`display_name()` 是去掉扩展名的名字，
    /// 用于图谱/序列标题与窗口标题）。
    pub fn file_label(&self) -> String {
        match &self.path {
            Some(p) => p
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| self.display_name()),
            None => self.display_name(),
        }
    }

    /// 只能原样保存的原生文档类别（比对/凝胶…）。
    pub fn native_class(&self) -> Option<&'static str> {
        self.native_only.as_ref().map(|n| n.class)
    }

    /// 原生只读文档的原始字节数（`None` = 普通可编辑文档）。
    pub fn native_bytes(&self) -> Option<usize> {
        self.native_only.as_ref().map(|n| n.bytes.len())
    }

    // ---- 导出 ----

    /// 图谱 SVG（与屏幕一致，优先启用的原生图谱名）。
    pub fn export_map_svg(&self, enzymes: &[usize], style: &hx_render::svg::Style) -> hx_algo::export::Export {
        hx_algo::export::map_svg(&self.session, enzymes, style, &self.map_name())
    }

    /// 序列视图 SVG。
    pub fn export_sequence_svg(
        &self,
        enzymes: &[usize],
        style: &hx_render::svg::Style,
    ) -> hx_algo::export::Export {
        hx_algo::export::sequence_svg(&self.session, enzymes, style, &self.display_name())
    }

    // ---- 编辑 ----

    pub fn apply(&mut self, edit: Edit) -> Result<(), String> {
        self.session.apply(edit)
    }

    /// 保存前把本地折叠出来的结构写进 block 34。
    ///
    /// SnapGene 的 `Secondary Structure` 页**只显示块里存的结构**：文件没有
    /// block 34 时它会一直转圈（实测：同一段 72 nt tRNA，带块的文件秒开，
    /// 去掉块后 20 秒仍在转、进程 0% CPU）。所以只要本地折叠结果与当前序列
    /// 对得上，就把这块补上——这样我们写出的 `.rna` 在 SnapGene 里也能直接
    /// 看到结构。
    ///
    /// 原文件本来就带 block 34 时**不碰**：那是 SnapGene 自己算的，次优结构
    /// 列表比我们全。没有配分函数（`probabilityMatrix` 为空）时也不写——
    /// 半个块会让 SnapGene 显示假的置信度颜色。
    fn inject_rna_structure(&self, bytes: Vec<u8>) -> Vec<u8> {
        let Some(cache) = self.rna_cache.as_ref() else {
            return bytes;
        };
        let Ok(st) = &cache.structure else {
            return bytes;
        };
        if st.probs.is_empty() || st.stats.is_none() || cache.seq != self.session.sequence() {
            return bytes;
        }
        let Ok(mut file) = DnaFile::parse(&bytes) else {
            return bytes;
        };
        if file.header.molecule_kind() != Some(MoleculeKind::Rna) || file.block(34).is_some() {
            return bytes;
        }
        file.blocks
            .push(hx_formats::block::Block::new(34, st.to_payload()));
        file.to_bytes()
    }

    pub fn undo(&mut self) -> Option<String> {
        if let Some(n)=&mut self.native_only{let bytes=n.undo.pop()?;n.redo.push(std::mem::replace(&mut n.bytes,bytes));n.dirty=n.bytes!=n.saved_bytes;n.blocks=DnaFile::parse(&n.bytes).ok()?.blocks.iter().map(|b|(b.block_type,b.payload.len())).collect();return Some("Edit Gel".into());}
        self.session.undo()
    }

    pub fn redo(&mut self) -> Option<String> {
        if let Some(n)=&mut self.native_only{let bytes=n.redo.pop()?;n.undo.push(std::mem::replace(&mut n.bytes,bytes));n.dirty=n.bytes!=n.saved_bytes;n.blocks=DnaFile::parse(&n.bytes).ok()?.blocks.iter().map(|b|(b.block_type,b.payload.len())).collect();return Some("Edit Gel".into());}
        self.session.redo()
    }

    /// 保存到 `path`：未编辑时逐字节原样写出；有编辑时按写回契约重建，
    /// 并写入 block 7 历史树（与 `hx-cli edit` 同实现，T10 oracle 已验证显示）。
    ///
    /// 返回 `(字节数, 历史树节点数)`；保存成功后以新文件为新基准。
    pub fn save_to(&mut self, path: &Path) -> Result<(usize, usize), String> {
        // 比对/凝胶：字节原样写回（我们没有可编辑的视图）
        if let Some(n) = &mut self.native_only {
            std::fs::write(path, &n.bytes).map_err(|e| format!("write {}: {e}", path.display()))?;
            self.path = Some(path.to_path_buf());
            n.dirty=false;
            n.saved_bytes=n.bytes.clone();
            return Ok((n.bytes.len(), 0));
        }
        let doc_name = path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| path.display().to_string());
        let src_name = self
            .path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| doc_name.clone());
        let (bytes, nodes) = self
            .session
            .save_bytes_with_history(&doc_name, &src_name, None)?;
        let mut bytes = self.inject_rna_structure(bytes);
        // Default names follow the path; explicit custom names remain independent.
        // Do this in the outgoing copy, so failed saves cannot change the session.
        if self.custom_map_label().is_none() && (self.path.as_deref()!=Some(path)) {
            let mut file=DnaFile::parse(&bytes).map_err(|e|e.to_string())?;
            let notes=file.block(6).map(|b|b.payload.clone()).unwrap_or_else(||b"<Notes></Notes>".to_vec());
            let stem=path.file_stem().map(|s|s.to_string_lossy().into_owned()).unwrap_or_else(||"Untitled".into());
            let notes=hx_formats::notes::set_field(&notes,"CustomMapLabel",&stem);
            let notes=hx_formats::notes::set_field(&notes,"UseCustomMapLabel","0");
            if let Some(b)=file.blocks.iter_mut().find(|b|b.block_type==6){b.payload=notes;}
            else{file.blocks.push(hx_formats::Block::new(6,notes));}
            bytes=file.to_bytes();
        }
        std::fs::write(path, &bytes).map_err(|e| format!("write {}: {e}", path.display()))?;
        self.session.commit_saved(bytes)?;
        self.path = Some(path.to_path_buf());
        Ok((self.session.bytes().len(), nodes))
    }
}

#[cfg(test)]
mod map_name_tests {
    use super::*;

    #[test]
    fn newly_created_rna_title_follows_first_save_and_save_as_in_native_notes() {
        let mut state=crate::newfile::NewFileState::new(MoleculeKind::Rna);state.seq="ACGU".repeat(20);
        let mut doc=crate::newfile::build_doc(&state).unwrap();
        let dir=PathBuf::from("/tmp/opencode").join(format!("rna-default-title-{}",std::process::id()));std::fs::create_dir_all(&dir).unwrap();
        for name in ["445566.rna","renamed.rna"] {
            let path=dir.join(name);assert!(!path.exists());doc.save_to(&path).unwrap();
            let reopened=Document::from_path(&path).unwrap();let stem=path.file_stem().unwrap().to_str().unwrap();
            assert_eq!(reopened.map_name(),stem);
            let notes=reopened.session.notes_xml();
            assert_eq!(hx_formats::notes::field(&notes,"CustomMapLabel").as_deref(),Some(stem));
            assert_eq!(hx_formats::notes::field(&notes,"UseCustomMapLabel").as_deref(),Some("0"));
            std::fs::remove_file(path).unwrap();
        }
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn create_dialog_name_is_not_an_explicit_custom_title() {
        let mut state=crate::newfile::NewFileState::new(MoleculeKind::Rna);state.seq="ACGU".repeat(20);state.name="445566.rna".into();
        let doc=Document::from_imported(state.build_bytes(),state.name.clone()).unwrap();
        assert_eq!(doc.map_name(),"445566");assert_eq!(doc.custom_map_label(),None);
    }

    #[test]
    fn disabled_or_missing_custom_map_label_falls_back_to_filename_in_map_and_export() {
        for enabled in [false,true] {
            let mut new=hx_formats::write::NewFile::of(MoleculeKind::Rna,b"ACGT".repeat(20));
            new.name=Some("Untitled".into());
            let mut file=DnaFile::parse(&new.to_bytes()).unwrap();
            if !enabled {
                let b=file.blocks.iter_mut().find(|b|b.block_type==6).unwrap();
                b.payload=String::from_utf8(b.payload.clone()).unwrap().replace("<UseCustomMapLabel>1</UseCustomMapLabel>","<UseCustomMapLabel>0</UseCustomMapLabel>").into_bytes();
            }
            let mut doc=Document::from_bytes(file.to_bytes()).unwrap();doc.path=Some(PathBuf::from("445566.rna"));
            let name=if enabled {"Untitled"}else{"445566"};assert_eq!(doc.map_name(),name);
            let svg=doc.export_map_svg(&[],&hx_render::svg::Style::map(900.0,900.0)).svg;
            assert!(svg.contains(&format!(">{name}<")));assert_eq!(doc.file_label(),"445566.rna");
        }
        let mut doc=Document::from_bytes(hx_formats::write::NewFile::of(MoleculeKind::Rna,b"ACGT".to_vec()).to_bytes()).unwrap();
        doc.path=Some(PathBuf::from("without-custom-name.rna"));assert_eq!(doc.map_name(),"without-custom-name");
    }

    #[test]
    fn save_as_keeps_the_native_custom_title_instead_of_renaming_it() {
        let mut new=hx_formats::write::NewFile::of(MoleculeKind::Rna,b"ACGT".repeat(20));new.name=Some("Untitled".into());
        let bytes=new.to_bytes();let mut doc=Document::from_bytes(bytes.clone()).unwrap();
        let dir=PathBuf::from("/tmp/opencode").join(format!("rna-name-regression-{}",std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();let path=dir.join("445566.rna");
        assert!(!path.exists(),"test must not overwrite an existing file");doc.save_to(&path).unwrap();
        let reopened=Document::from_path(&path).unwrap();
        assert_eq!(reopened.file_label(),"445566.rna");assert_eq!(reopened.map_name(),"Untitled");
        assert_eq!(std::fs::read(&path).unwrap(),bytes);
        std::fs::remove_file(path).unwrap();std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn rna_custom_map_title_is_not_the_filename() {
        let mut file=hx_formats::write::NewFile::of(MoleculeKind::Rna,[b"ACGT".repeat(91),vec![b'A']].concat());
        file.name=Some("Untitled".into());
        let bytes=file.to_bytes();let mut doc=Document::from_bytes(bytes.clone()).unwrap();
        doc.path=Some(PathBuf::from("445566.rna"));
        assert_eq!(doc.file_label(),"445566.rna");
        assert_eq!(doc.display_name(),"445566");
        let ctx=eframe::egui::Context::default();
        let mut app=crate::HelixApp::preview(&ctx,doc);
        app.map_layout=hx_render::MapLayout::Circular;
        let mut output=ctx.run_ui(eframe::egui::RawInput{screen_rect:Some(eframe::egui::Rect::from_min_size(eframe::egui::Pos2::ZERO,eframe::egui::vec2(1800.0,1100.0))),..Default::default()},|ui|{
            crate::view_map::draw(ui,&mut app);
        });
        let text=output.shapes.iter().filter_map(|s|if let eframe::egui::Shape::Text(t)=&s.shape{Some(t.galley.job.text.as_str())}else{None}).collect::<Vec<_>>().join("\n");
        output.textures_delta.clear();
        assert!(text.contains("Untitled"),"map must honor the enabled native label, not the basename: {text}");
        for mark in ["50","100","150","200","250","300","350"] {assert!(text.lines().any(|s|s==mark),"RNA ruler is missing {mark}: {text}");}
        assert!(!text.lines().any(|s|s=="20" || s=="40"),"RNA map should use the sparse numbering in the supplied SG reference");
        assert!(!text.contains("445566"));output.textures_delta.clear();
        assert_eq!(app.doc.as_ref().unwrap().session.save_bytes().unwrap(),bytes,"display must not rewrite native metadata");
    }
}

#[cfg(test)]
mod sticky_document_tests {
    use super::*;
    use hx_algo::cloning::{ligate, product_bytes, Source};
    use hx_formats::dna::EndProperties;

    fn ligation_bytes() -> Vec<u8> {
        let left = Source::new("left", b"TGGGGAATT".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: -1, downstream_stickiness: 4, ..Default::default()
        });
        let right = Source::new("right", b"AATTCCCCAA".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: 4, downstream_stickiness: 2, ..Default::default()
        });
        product_bytes(&ligate(&[left, right], false).unwrap(), "Ligation")
    }

    #[test]
    fn ligation_document_titles_editing_and_cloning_share_correct_coordinates() {
        let bytes = ligation_bytes();
        let mut doc = Document::from_bytes(bytes.clone()).unwrap();
        doc.name_hint = Some("Ligation.dna".into());
        assert_eq!(doc.len(), 12);
        assert!(doc.title_summary().contains("12 bp"));
        assert!(doc.summary().contains("12 bp"));
        assert_eq!(doc.sequence(), b"GGGGAATTCCCC");
        assert_eq!(doc.session.save_bytes().unwrap(), bytes);
        let source = doc.cloning_source("again");
        assert_eq!(source.seq, b"TGGGGAATTCCCCAA");
        assert_eq!(source.top_seq(), doc.sequence());
        assert_eq!(source.ends, doc.end_properties());
        doc.apply(Edit::Insert { at: 0, seq: b"AC".to_vec() }).unwrap();
        doc.apply(Edit::Insert { at: doc.len(), seq: b"GT".to_vec() }).unwrap();
        assert_eq!(doc.len(), 16);
        let saved = doc.session.save_bytes().unwrap();
        let native = DnaFile::parse(&saved).unwrap();
        assert_eq!(native.sequence().unwrap().sequence, b"TACGGGGAATTCCCCGTAA");
        assert_eq!(native.end_properties(), source.ends);
        let loaded = Document::from_bytes(saved).unwrap();
        assert_eq!(loaded.sequence(), doc.sequence());
        assert_eq!(loaded.cloning_source("again").seq, doc.session.stored_sequence());
        doc.undo().unwrap(); doc.undo().unwrap();
        assert_eq!(doc.session.save_bytes().unwrap(), bytes);
    }
}

/// 酶位点命中缓存：`hits[i]` 与 `selected[i]` 一一对应（0-based 起点列表）。
#[derive(Default, Clone)]
pub struct SiteHits {
    pub selected: Vec<usize>,
    pub hits: Vec<Vec<usize>>,
}

impl SiteHits {
    /// 重新计算（`selected` = [`ENZYMES`] 下标）。
    pub fn recompute(&mut self, seq: &[u8], circular: bool) {
        self.hits.clear();
        for &idx in &self.selected {
            let e = match ENZYMES.get(idx) {
                Some(e) => e,
                None => {
                    self.hits.push(Vec::new());
                    continue;
                }
            };
            self.hits.push(
                hx_enzymes::sets::positions(idx, seq, circular)
                    .into_iter()
                    .map(|p| p as usize)
                    .collect(),
            );
        }
    }

    /// 非 DNA 分子：只保留选择、清空命中。
    pub fn clear_hits(&mut self) {
        self.hits = vec![Vec::new(); self.selected.len()];
    }

    pub fn total(&self) -> usize {
        self.hits.iter().map(|v| v.len()).sum()
    }

    /// 命中是否为唯一（该酶在序列中只有一处）。
    pub fn is_unique(&self, i: usize) -> bool {
        self.hits.get(i).map(|v| v.len() == 1).unwrap_or(false)
    }
}

/// 常用酶（SnapGene 风格默认集，名称大小写不敏感）。
/// 常用酶集（名称表在 [`hx_enzymes::COMMON`]）。
pub const COMMON_ENZYMES: &[&str] = hx_enzymes::COMMON;

/// 由名称表解析为 [`ENZYMES`] 下标（跳过未知名称）。
pub fn indices_of(names: &[&str]) -> Vec<usize> {
    let mut out = Vec::new();
    for n in names {
        if let Some(e) = hx_enzymes::by_name(n) {
            if let Some(i) = ENZYMES.iter().position(|x| std::ptr::eq(x, e)) {
                out.push(i);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}
