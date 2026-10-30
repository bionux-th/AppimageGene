//! 分子克隆引擎（SnapGene `Actions ▸ … Cloning` 的计算部分）。
//!
//! 各方法分别解释片段边界、双链末端或重组位点，再投影产物特征：
//!
//! - [`Source::slice`]：从一条来源序列取一段（支持环状跨原点、反向互补），
//!   并把落在这一段上的特征一起搬到片段坐标系；
//! - [`assemble`]：把若干片段拼成产物，自动识别相邻片段之间的重叠并去重，
//!   同时把特征平移到产物坐标。
//! - [`restriction_cuts`] / [`digest_fragment`] / [`ligate`]：双链酶切、补平、
//!   末端相容性和线性产物外侧末端；不得用无条件 assemble 代替连接校验。
//!
//! 方法的实际等价覆盖见 `docs/cloning-equivalence.md`；共用拼接原语不代表
//! 各方法的化学语义或 SnapGene 工作流相同。
//!
//! UI 侧（`hx-app/src/clonedlg.rs`）只负责收集「用哪条序列、哪一段、什么方向」。

use hx_core::MoleculeFlags;
use hx_formats::annot::{FeatureAnnotation, FeatureSegment, Range};
use hx_formats::write::NewFile;

mod digest;
pub mod golden_gate;
pub mod homology;
pub mod gateway;
pub use digest::{blunt_ends, digest_fragment, ligate, ligate_reaction, restriction_product, restriction_cuts, RestrictionCut};

/// 克隆方式。对应 SnapGene `Actions` 菜单里的每一条命令。
///
/// 变体名字直接决定对话框的页签集合、默认产物名和重叠长度；带数字的
/// Gateway 变体表示一次反应里参与重组的片段数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    Restriction,
    RestrictionMulti,
    Ligation,
    Gibson,
    GibsonMulti,
    GibsonAssemble,
    InFusion,
    InFusionMulti,
    NeBuilder,
    NeBuilderMulti,
    NeBuilderAssemble,
    GoldenGate,
    GoldenGateMulti,
    TaCloning,
    GcCloning,
    TopoTa,
    TopoBlunt,
    TopoDirectional,
    GatewayBp,
    GatewayLr(usize),
    GatewayBpLr(usize),
}

impl Method {
    /// 对话框标题（不含方法名后的冒号与产物名）。
    pub fn title(self) -> String {
        match self {
            Method::Restriction | Method::RestrictionMulti => {
                "Restriction and Insertion Cloning".into()
            }
            Method::Ligation => "Ligate Fragments".into(),
            Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble => {
                "Gibson Assembly®".into()
            }
            Method::InFusion | Method::InFusionMulti => "In-Fusion® Cloning".into(),
            Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble => {
                "NEBuilder® HiFi DNA Assembly".into()
            }
            Method::GoldenGate | Method::GoldenGateMulti => "Golden Gate Assembly".into(),
            Method::TaCloning => "TA Cloning".into(),
            Method::GcCloning => "GC Cloning".into(),
            Method::TopoTa => "TOPO® TA Cloning".into(),
            Method::TopoBlunt => "Blunt TOPO® Cloning".into(),
            Method::TopoDirectional => "Directional TOPO® Cloning".into(),
            Method::GatewayBp => "Gateway® BP Cloning".into(),
            Method::GatewayLr(n) => match n {
                1 => "Gateway® LR Cloning: One Fragment".into(),
                2 => "Gateway® LR Cloning: Two Fragments".into(),
                3 => "Gateway® LR Cloning: Three Fragments".into(),
                _ => "Gateway® LR Cloning: Four Fragments".into(),
            },
            Method::GatewayBpLr(n) => match n {
                1 => "Gateway® BP+LR Cloning: One Fragment".into(),
                2 => "Gateway® BP+LR Cloning: Two Fragments".into(),
                3 => "Gateway® BP+LR Cloning: Three Fragments".into(),
                _ => "Gateway® BP+LR Cloning: Four Fragments".into(),
            },
        }
    }

    /// 对话框里 Insert Fragment 这一动作的名字（标题用）。
    pub fn action(self) -> &'static str {
        match self {
            Method::Restriction => "Insert Fragment",
            Method::RestrictionMulti => "Insert Multiple Fragments",
            Method::Ligation => "Ligate Fragments",
            Method::Gibson | Method::InFusion | Method::NeBuilder | Method::GoldenGate => {
                "Insert Fragment"
            }
            Method::GibsonMulti
            | Method::InFusionMulti
            | Method::NeBuilderMulti
            | Method::GoldenGateMulti => "Insert Multiple Fragments",
            Method::GibsonAssemble | Method::NeBuilderAssemble => "Assemble Multiple Fragments",
            Method::TaCloning => "TA Cloning",
            Method::GcCloning => "GC Cloning",
            Method::TopoTa => "TOPO® TA Cloning",
            Method::TopoBlunt => "Blunt TOPO® Cloning",
            Method::TopoDirectional => "Directional TOPO® Cloning",
            Method::GatewayBp => "BP Cloning",
            Method::GatewayLr(_) => "LR Cloning",
            Method::GatewayBpLr(_) => "BP + LR Cloning",
        }
    }

    /// 产物默认文件名（SnapGene 的命名规则）。
    pub fn product_name(self) -> &'static str {
        match self {
            Method::Restriction
            | Method::RestrictionMulti
            | Method::TaCloning
            | Method::GcCloning
            | Method::TopoTa
            | Method::TopoBlunt
            | Method::TopoDirectional => "Cloned.dna",
            Method::Ligation => "Ligation.dna",
            Method::Gibson
            | Method::GibsonMulti
            | Method::GibsonAssemble
            | Method::InFusion
            | Method::InFusionMulti
            | Method::NeBuilder
            | Method::NeBuilderMulti
            | Method::NeBuilderAssemble
            | Method::GoldenGate
            | Method::GoldenGateMulti => "Assembled.dna",
            Method::GatewayBp => "Entry Clone.dna",
            Method::GatewayLr(_) | Method::GatewayBpLr(_) => "Expression Clone.dna",
        }
    }

    /// 该方法一次要选几个插入片段（载体另算）。
    pub fn fragment_count(self) -> usize {
        match self {
            Method::GatewayLr(n) | Method::GatewayBpLr(n) => n.clamp(1, 4),
            Method::RestrictionMulti
            | Method::GibsonMulti
            | Method::GibsonAssemble
            | Method::InFusionMulti
            | Method::NeBuilderMulti
            | Method::NeBuilderAssemble
            | Method::GoldenGateMulti => 2,
            Method::Ligation => 2,
            _ => 1,
        }
    }

    pub fn is_restriction(self) -> bool {
        matches!(self, Method::Restriction | Method::RestrictionMulti)
    }

    pub fn is_ligation(self) -> bool {
        matches!(self, Method::Ligation)
    }

    pub fn is_gateway(self) -> bool {
        matches!(
            self,
            Method::GatewayBp | Method::GatewayLr(_) | Method::GatewayBpLr(_)
        )
    }

    pub fn is_bp(self) -> bool {
        matches!(self, Method::GatewayBp | Method::GatewayBpLr(_))
    }

    pub fn is_lr(self) -> bool {
        matches!(self, Method::GatewayLr(_) | Method::GatewayBpLr(_))
    }

    /// TA/GC/TOPO 用粘性末端而不是同源臂。
    pub fn is_ta_topo(self) -> bool {
        matches!(
            self,
            Method::TaCloning
                | Method::GcCloning
                | Method::TopoTa
                | Method::TopoBlunt
                | Method::TopoDirectional
        )
    }

    pub fn is_multiple(self) -> bool {
        match self {
            Method::GatewayLr(n) | Method::GatewayBpLr(n) => n > 1,
            Method::RestrictionMulti
            | Method::GibsonMulti
            | Method::GibsonAssemble
            | Method::InFusionMulti
            | Method::NeBuilderMulti
            | Method::NeBuilderAssemble
            | Method::GoldenGateMulti => true,
            _ => false,
        }
    }

    /// 该方法历史上的同源臂上限（bp），仅作为参考值保留。
    ///
    /// 同源组装现在走 `homology::assemble_circular`，上限是片段长度而不是
    /// 固定的 40/20 bp，因此这里不再参与 Gibson/In-Fusion 的实际判定。
    pub fn max_overlap(self) -> usize {
        match self {
            Method::Restriction | Method::RestrictionMulti | Method::Ligation => 0,
            Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble => 40,
            Method::InFusion | Method::InFusionMulti => 20,
            Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble => 40,
            Method::GoldenGate | Method::GoldenGateMulti => 8,
            Method::TaCloning | Method::GcCloning => 3,
            Method::TopoTa | Method::TopoBlunt | Method::TopoDirectional => 3,
            Method::GatewayBp | Method::GatewayLr(_) | Method::GatewayBpLr(_) => 25,
        }
    }

    /// 该方法要求的最短同源臂（bp）。限制性酶切克隆不靠同源臂。
    pub fn min_overlap(self) -> usize {
        match self {
            Method::Restriction | Method::RestrictionMulti | Method::Ligation => 0,
            Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble => 15,
            Method::InFusion | Method::InFusionMulti => 12,
            Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble => 15,
            Method::GoldenGate | Method::GoldenGateMulti => 4,
            Method::TaCloning | Method::GcCloning => 1,
            Method::TopoTa | Method::TopoBlunt | Method::TopoDirectional => 0,
            Method::GatewayBp | Method::GatewayLr(_) | Method::GatewayBpLr(_) => 25,
        }
    }

    /// 产物是否环状（TA/TOPO/限制性/Gateway 都成环；线性连接看两端）。
    pub fn circular_product(self) -> bool {
        !matches!(self, Method::Ligation)
    }
}

/// 一条来源序列（打开的文档）。
#[derive(Debug, Clone)]
pub struct Source {
    pub label: String,
    pub seq: Vec<u8>,
    pub circular: bool,
    pub features: Vec<FeatureAnnotation>,
    /// 来源文件里的引物（克隆对话框的 Primers 子页签显示，不参与拼接）。
    pub primers: Vec<hx_formats::annot::PrimerAnnotation>,
    /// `block 8` 的末端属性（粘性末端 / 修饰 / 原始首碱基）。
    pub ends: hx_formats::dna::EndProperties,
}

/// 从来源取出的一段（已是线性、已按方向排好）。
#[derive(Debug, Clone)]
pub struct Slice {
    /// 片段标签（`probe.dna` / `Fragment 1`）。
    pub label: String,
    pub seq: Vec<u8>,
    /// 特征坐标已换算到片段坐标系（0-based 半开）。
    pub features: Vec<FeatureAnnotation>,
    pub reversed: bool,
}

impl Slice {
    /// 从片段里再取一段（不跨原点；用于 Gateway 的位点裁切）。
    pub fn sub(&self, start: usize, end: usize) -> Slice {
        let (a, b) = (start.min(self.seq.len()), end.min(self.seq.len()));
        let mut features = Vec::new();
        for f in &self.features {
            let mut segs = Vec::new();
            for seg in &f.segments {
                let lo = seg.range.start.max(a);
                let hi = seg.range.end.min(b);
                if hi > lo {
                    segs.push(FeatureSegment {
                        range: Range::new(lo - a, hi - a),
                        color: seg.color.clone(),
                        ..seg.clone()
                    });
                }
            }
            if segs.is_empty() {
                continue;
            }
            let mut nf = f.clone();
            nf.max_run_on=f.max_run_on.filter(|end|*end>a && *end<=b).map(|end|end-a);
            let first = segs.remove(0);
            nf.segments = {
                let mut v = vec![first];
                v.extend(segs);
                v
            };
            features.push(nf);
        }
        Slice {
            label: self.label.clone(),
            seq: self.seq.get(a..b).unwrap_or(&[]).to_vec(),
            features,
            reversed: self.reversed,
        }
    }

    /// 去掉末尾 `n` 个碱基（Gateway 环化时避免核心序列重复）。
    pub fn truncate_end(&mut self, n: usize) {
        let keep = self.seq.len().saturating_sub(n);
        self.seq.truncate(keep);
        for f in &mut self.features {
            f.max_run_on=f.max_run_on.filter(|end|*end<=keep);
            f.segments.retain_mut(|s| {
                if s.range.start >= keep {
                    return false;
                }
                s.range.end = s.range.end.min(keep);
                s.range.end > s.range.start
            });
        }
    }

    /// 在 5′ 端加引物尾巴，特征坐标整体后移。
    pub fn prepend(&mut self, bases: &[u8]) {
        if bases.is_empty() {
            return;
        }
        let n = bases.len();
        self.seq.splice(0..0, bases.iter().copied());
        for f in &mut self.features {
            f.max_run_on=f.max_run_on.map(|end|end+n);
            for seg in &mut f.segments {
                seg.range = Range::new(seg.range.start + n, seg.range.end + n);
            }
        }
    }

    /// 在 3′ 端加引物尾巴（已是正向链序列）。
    pub fn append(&mut self, bases: &[u8]) {
        self.seq.extend_from_slice(bases);
    }
}

/// Gateway 重组位点的核心序列。attB/attP/attL/attR 共享同一段 25 bp 核心，
/// 差别只在两侧翼序列；模拟克隆时用核心做同源臂即可得到正确大小的产物。
pub const ATT_B1: &[u8] = b"ACAAGTTTGTACAAAAAAGCAGGCT";
pub const ATT_B2: &[u8] = b"ACCCAGCTTTCTTGTACAAAGTGGT";

/// 真实 attL1（pENTR/pDONR221 的 BP 产物，100 bp，正向）。
/// 前 95 bp 来自供体 attP1，末 5 bp 来自片段 attB1（`ATT_B1[20..]`）。
pub const ATT_L1: &[u8] = b"CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCT";
/// 真实 attL2（100 bp，产物里通常画成反向箭头）。
/// 前 22 bp 来自片段 attB2（`ATT_B2[..22]`），其余 78 bp 来自供体 attP2。
pub const ATT_L2: &[u8] = b"ACCCAGCTTTCTTGTACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG";

/// BP 的 attL1：供体 attP1 起点前 1 bp 起的 95 bp + attB1 的末 5 bp。
fn bp_attl1(donor: &[u8], attp1_start: usize) -> Vec<u8> {
    let start = attp1_start.saturating_sub(1);
    let end = (start + 95).min(donor.len());
    let mut v = donor[start..end].to_vec();
    v.extend_from_slice(&ATT_B1[ATT_B1.len() - 5..]);
    v
}

/// BP 的 attL2：attB2 的前 22 bp + 供体 attP2 末端的 78 bp（含 attP2 后 1 bp）。
fn bp_attl2(donor: &[u8], attp2_end: usize) -> Vec<u8> {
    let mut v = ATT_B2[..22].to_vec();
    let end = (attp2_end + 1).min(donor.len());
    let start = end.saturating_sub(78);
    v.extend_from_slice(&donor[start..end]);
    v
}

/// 拼接处。
#[derive(Debug, Clone, Copy)]
pub struct Junction {
    /// 拼接点在产物上的 0-based 坐标。
    pub at: usize,
    /// 被去掉的重复碱基数（同源臂长度）。
    pub overlap: usize,
    /// 真实重叠序列（用于界面提示）。
    pub seq_start: usize,
    pub seq_end: usize,
}

/// 产物。
#[derive(Debug, Clone)]
pub struct Product {
    pub seq: Vec<u8>,
    pub features: Vec<FeatureAnnotation>,
    pub junctions: Vec<Junction>,
    /// 每个来源片段在产物上的范围（0-based 半开），用于预览里高亮插入片段。
    pub spans: Vec<(usize, usize)>,
    pub circular: bool,
    /// 线性产物的外侧末端；序列沿用 Source 的标记碱基约定。
    pub ends: hx_formats::dna::EndProperties,
}

impl Product {
    pub fn display_len(&self) -> usize {
        self.seq.len()
            .saturating_sub(self.ends.upstream_stickiness.min(0).unsigned_abs() as usize)
            .saturating_sub(self.ends.downstream_stickiness.max(0) as usize)
    }
}

impl Source {
    pub fn new(label: impl Into<String>, seq: Vec<u8>, circular: bool) -> Self {
        Self {
            label: label.into(),
            seq,
            circular,
            features: Vec::new(),
            primers: Vec::new(),
            ends: hx_formats::dna::EndProperties::default(),
        }
    }

    pub fn with_ends(mut self, ends: hx_formats::dna::EndProperties) -> Self {
        self.ends = ends;
        self
    }

    /// 存储序列里「顶链」的区间（去掉两端属于底链突出端的标记碱基）。
    pub fn top_range(&self) -> (usize, usize) {
        let range = self.ends.top_range(self.seq.len());
        (range.start, range.end)
    }

    /// 顶链（显示用序列）。粘性末端的标记碱基不属于顶链。
    pub fn top_seq(&self) -> &[u8] {
        let (a, b) = self.top_range();
        &self.seq[a..b]
    }

    /// 界面显示的长度（顶链长度，等于 SnapGene 的 `N bp`）。
    pub fn display_len(&self) -> usize {
        let (a, b) = self.top_range();
        b - a
    }

    /// 是否是带粘性末端的线性来源（TA/TOPO 载体或 PCR 产物）。
    pub fn has_sticky_ends(&self) -> bool {
        self.ends.upstream_stickiness != 0 || self.ends.downstream_stickiness != 0
    }

    /// 反向使用（SnapGene 的 `←`）：反向互补序列与特征，两端粘性末端互换。
    pub fn flipped(&self) -> Source {
        let mut seq = self.seq.clone();
        revcomp(&mut seq);
        let mut features = self.features.clone();
        for f in &mut features {
            reverse_feature(f, seq.len());
        }
        let ends = hx_formats::dna::EndProperties {
            upstream_stickiness: self.ends.downstream_stickiness,
            downstream_stickiness: self.ends.upstream_stickiness,
            upstream_modification: self.ends.downstream_modification.clone(),
            downstream_modification: self.ends.upstream_modification.clone(),
            original_first_base: None,
        };
        Source {
            label: self.label.clone(),
            seq,
            circular: self.circular,
            features,
            primers: self.primers.clone(),
            ends,
        }
    }

    /// 产物旋转量：让存储序列第 `OriginalFirstBase` 个碱基落在产物位置 1。
    ///
    /// `concat_start` 是拼接串在**存储坐标**里的起点（顶链起点），`total` 是
    /// 拼接串长度（= 产物长度）。
    pub fn rotation_offset(&self, concat_start: usize, total: usize) -> usize {
        let Some(ofb) = self.ends.original_first_base else {
            return 0;
        };
        (ofb as i64 - concat_start as i64).rem_euclid(total.max(1) as i64) as usize
    }

    pub fn with_features(mut self, features: Vec<FeatureAnnotation>) -> Self {
        self.features = features;
        self
    }

    pub fn with_primers(
        mut self,
        primers: Vec<hx_formats::annot::PrimerAnnotation>,
    ) -> Self {
        self.primers = primers;
        self
    }

    /// 取 `[start, end)`：环状时 `start > end` 表示跨原点；`reversed` 取反向互补。
    pub fn slice(&self, start: usize, end: usize, reversed: bool) -> Slice {
        let n = self.seq.len();
        if n == 0 {
            return Slice {
                label: self.label.clone(),
                seq: Vec::new(),
                features: Vec::new(),
                reversed,
            };
        }
        // 归一化成「若干线性区段」（跨原点拆两段）
        let pieces: Vec<(usize, usize)> = if self.circular && start == end {
            // 环状序列在同一个切点打开 = 整圈
            let at = start % n;
            vec![(at, n), (0, at)]
        } else if start < end {
            vec![(start.min(n), end.min(n))]
        } else if start == end {
            Vec::new()
        } else if self.circular {
            vec![(start.min(n), n), (0, end.min(n))]
        } else {
            Vec::new()
        };
        let mut seq: Vec<u8> = Vec::new();
        for (a, b) in &pieces {
            seq.extend_from_slice(&self.seq[*a..*b]);
        }
        let mut feats = self.project_features(&pieces);
        if reversed {
            revcomp(&mut seq);
            for f in &mut feats {
                reverse_feature(f, seq.len());
            }
        }
        Slice {
            label: self.label.clone(),
            seq,
            features: feats,
            reversed,
        }
    }

    /// 把特征投影到给定区段（区段按来源顺序排列，拼接后即片段坐标）。
    fn project_features(&self, pieces: &[(usize, usize)]) -> Vec<FeatureAnnotation> {
        let mut out = Vec::new();
        let source_len = self.seq.len();
        for f in &self.features {
            let mut segs: Vec<FeatureSegment> = Vec::new();
            for seg in &f.segments {
                // 片段自己的范围也允许跨原点（`start > end`）
                let ranges: Vec<(usize, usize)> = if seg.range.start <= seg.range.end {
                    vec![(seg.range.start, seg.range.end)]
                } else {
                    vec![
                        (seg.range.start, source_len),
                        (0, seg.range.end),
                    ]
                };
                for (rs, re) in ranges {
                    // 逐区段求交集，并换算到片段坐标
                    let mut offset = 0usize;
                    for (ps, pe) in pieces {
                        let a = rs.max(*ps);
                        let b = re.min(*pe);
                        if b > a {
                            segs.push(FeatureSegment {
                                range: Range::new(offset + a - ps, offset + b - ps),
                                color: seg.color.clone(),
                                ..seg.clone()
                            });
                        }
                        offset += pe - ps;
                    }
                }
            }
            if segs.is_empty() {
                continue;
            }
            let mut nf = f.clone();
            // 首段为主段，其余按出现的先后作为附加片段
            nf.max_run_on=f.max_run_on.and_then(|end|{
                let mut offset=0;
                for &(a,b) in pieces {
                    if end>a && end<=b {return Some(offset+end-a);}
                    offset+=b-a;
                }
                None
            });
            let first = segs.remove(0);
            nf.segments = {
                let mut v = vec![first];
                v.extend(segs);
                v
            };
            out.push(nf);
        }
        out
    }
}

/// 从一条来源里按酶切位点取出被切开的线性片段（切成两段以上时由调用方挑）。
pub fn cut_positions(seq: &[u8], circular: bool, enzyme: &hx_enzymes::Enzyme) -> Vec<usize> {
    let mut out: Vec<usize> = restriction_cuts(seq, circular, enzyme)
        .into_iter().map(|c| c.top as usize).collect();
    out.sort_unstable();
    out.dedup();
    out
}

/// 突出端长度；不能用来给两侧切点同时加偏移。补平须处理具体片段。
pub fn blunt_shift(enzyme: &hx_enzymes::Enzyme) -> usize {
    (enzyme.site_len() as i64 - enzyme.reverse_cut as i64 - enzyme.cut as i64)
        .unsigned_abs() as usize
}

/// `a` 的后缀与 `b` 的前缀的最长公共长度（`min ≤ k ≤ max` 才算数）。
///
/// 上限同时受两段长度限制：不允许把某一段整个吃掉。
pub fn overlap_len(a: &[u8], b: &[u8], min: usize, max: usize) -> usize {
    let hi = a.len().min(b.len()).min(max);
    let lo = min.max(1);
    if hi < lo {
        return 0;
    }
    let mut best = 0usize;
    for k in (lo..=hi).rev() {
        if a[a.len() - k..].eq_ignore_ascii_case(&b[..k]) {
            best = k;
            break;
        }
    }
    best
}

/// 在序列（大小写不敏感）里找一个位点。
pub fn find_motif(seq: &[u8], motif: &[u8]) -> Option<usize> {
    if motif.is_empty() || seq.len() < motif.len() {
        return None;
    }
    seq.windows(motif.len())
        .position(|w| w.eq_ignore_ascii_case(motif))
}

/// 在一条来源里按特征名找 Gateway 位点，并带方向（1=正向、2=反向）。
/// 方向用于判断插入片段是否要反向互补。
fn site_feature_dir(src: &Source, names: &[&str]) -> Option<(usize, usize, i32)> {
    src.features.iter().find_map(|f| {
        let n = f.name.to_ascii_lowercase();
        if names.iter().any(|want| n == want.to_ascii_lowercase()) {
            f.segments
                .first()
                .map(|seg| (seg.range.start, seg.range.end, f.directionality))
        } else {
            None
        }
    })
}

/// 按 motif 找一对位点，兼顾跨原点与反向互补。返回 (start, end, reversed)，
/// 区间包含两个位点本身。
fn motif_pair(src: &Source, m1: &[u8], m2: &[u8]) -> Option<(usize, usize, bool)> {
    let scan = |seq: &[u8]| -> Option<(usize, usize, usize, usize)> {
        let a = find_motif(seq, m1)?;
        let b = find_motif(seq, m2)?;
        Some((a, a + m1.len(), b, b + m2.len()))
    };
    if let Some((a1, a2, b1, b2)) = scan(&src.seq) {
        return if a1 <= b1 {
            Some((a1, b2, false))
        } else {
            // 位点跨原点：取 B 位点开头到 A 位点结尾，再反向互补。
            Some((b1, a2, true))
        };
    }
    let rc = revcomp_str(&src.seq);
    let (a1, a2, b1, b2) = scan(&rc)?;
    let n = src.seq.len();
    // rc 坐标 → 原序列坐标（整段反向）。
    let (s, e) = if a1 <= b1 {
        (n - b2, n - a1)
    } else {
        (n - a2, n - b1)
    };
    Some((s, e, true))
}

/// Gateway BP/LR 重组。
///
/// SnapGene 的产物布局（用 pDONR221 + 650 bp attB 片段实测）：
///   BP: donor[0..attP1.start-1] + attL1(100) + gene + attL2(100) + donor[attP2.end+1..]
///   LR: dest[0..attR1.start] + attB1(25) + gene + attB2(25) + dest[attR2.end..]
/// 这样 BP 产物是 3144 bp、LR 产物是 6403 bp，和 SnapGene 的尺寸/位点位置一致。
pub fn gateway_recombine(
    inserts: &[&Source],
    vector: &Source,
    insert_label: &str,
    bp: bool,
) -> Result<Product, String> {
    if inserts.is_empty() {
        return Err("Please choose a source file.".into());
    }
    if !bp && inserts.len()>1 {return gateway::lr(inserts,vector);}
    if bp && inserts.len()!=1 {return Err("Each BP reaction requires one attB Fragment and one matched Donor Vector.".into());}
    if bp && !(vector.features.iter().any(|f|f.name.eq_ignore_ascii_case("attP1"))
        && vector.features.iter().any(|f|f.name.eq_ignore_ascii_case("attP2"))) {
        return gateway::bp(inserts[0],vector);
    }
    // 不能仅凭注释名字无条件切片，更不能把 attB 载体当成 attP/attR 载体。
    for source in std::iter::once(vector).chain(inserts.iter().copied()) {
        for feature in source.features.iter().filter(|f| is_att_name(&f.name)) {
            if feature.segments.is_empty() || feature.segments.iter().any(|s| {
                s.range.start >= source.seq.len() || s.range.end > source.seq.len() || s.range.start >= s.range.end
            }) {
                return Err(format!("The att-site coordinates in {} are invalid.", source.label));
            }
        }
    }
    let mut slices: Vec<Slice> = Vec::new();

    if bp {
        let has_p = vector
            .features
            .iter()
            .any(|f| f.name.eq_ignore_ascii_case("attP1"))
            && vector
                .features
                .iter()
                .any(|f| f.name.eq_ignore_ascii_case("attP2"));
        if !has_p {
            return Err("Please specify a suitable Donor Vector with attP1 and attP2 sites.".into());
        }
        let (p1, p2) = vector
            .features
            .iter()
            .find(|f| f.name.eq_ignore_ascii_case("attP1"))
            .zip(
                vector
                    .features
                    .iter()
                    .find(|f| f.name.eq_ignore_ascii_case("attP2")),
            )
            .map(|(a, b)| {
                (
                    (
                        a.segments[0].range.start,
                        a.segments[0].range.end,
                        a.directionality,
                    ),
                    (
                        b.segments[0].range.start,
                        b.segments[0].range.end,
                        b.directionality,
                    ),
                )
            })
            .ok_or_else(|| "Please specify a suitable Donor Vector.".to_string())?;
        let (s1, _e1, _d1) = p1;
        let (_s2, e2, _d2) = p2;
        if s1 == 0 || e2 >= vector.seq.len() {
            return Err("The Donor Vector att sites are too close to the ends.".into());
        }
        // 供体骨架前半段：attP1 起点前 1 bp 之前。
        slices.push(vector.slice(0, s1 - 1, false));
        // attL1 是**杂合位点**：95 bp 供体 attP1 左翼 + 5 bp 片段 attB1 右端。
        slices.push(Slice {
            label: "attL1".into(),
            seq: bp_attl1(&vector.seq, s1),
            features: Vec::new(),
            reversed: false,
        });

        for ins in inserts {
            let (b1, b2) =
                site_pair(ins, &["attB1"], &["attB2"], ATT_B1, ATT_B2).ok_or_else(|| {
                    format!(
                        "The {insert_label} must contain two different att sites that face each other."
                    )
                })?;
            slices.push(gene_slice(ins, b1, b2)?);
        }

        // attL2 也是杂合位点：22 bp 片段 attB2 左端 + 78 bp 供体 attP2 右翼。
        slices.push(Slice {
            label: "attL2".into(),
            seq: bp_attl2(&vector.seq, e2),
            features: Vec::new(),
            reversed: false,
        });
        // 供体骨架后半段：attP2 末 1 bp 之后。
        slices.push(vector.slice(e2 + 1, vector.seq.len(), false));

        let mut product = assemble(&slices, true, 0, 0);
        product.features.retain(|f| !is_att_name(&f.name));
        let pre = s1 - 1;
        let gene_len: usize = product.seq.len() - pre - 100 - 100 - (vector.seq.len() - e2 - 1);
        product.features.push(make_site_feature("attL1", pre, 100, 1));
        product
            .features
            .push(make_site_feature("attL2", pre + 100 + gene_len, 100, 2));
        product.junctions.clear();
        return Ok(product);
    }

    // LR：目的载体（attR1/attR2）
    let has_r = vector
        .features
        .iter()
        .any(|f| f.name.eq_ignore_ascii_case("attR1"))
        && vector
            .features
            .iter()
            .any(|f| f.name.eq_ignore_ascii_case("attR2"));
    if !has_r {
        return Err("Please specify a suitable Destination Vector with attR1 and attR2 sites.".into());
    }
    let (r1s, r2e) = vector
        .features
        .iter()
        .find(|f| f.name.eq_ignore_ascii_case("attR1"))
        .zip(
            vector
                .features
                .iter()
                .find(|f| f.name.eq_ignore_ascii_case("attR2")),
        )
        .map(|(a, b)| {
            (
                a.segments[0].range.start,
                b.segments[0].range.end,
            )
        })
        .ok_or_else(|| "Please specify a suitable Destination Vector.".to_string())?;
    if r2e >= vector.seq.len() {
        return Err("The Destination Vector att sites are too close to the ends.".into());
    }
    // 目的载体前半段（attR1 之前）。
    slices.push(vector.slice(0, r1s, false));
    let mut att_offsets: Vec<(usize, usize)> = Vec::new(); // (start, gene_len)
    let mut off = r1s;
    for ins in inserts {
        let (l1, l2) = match (
            site_feature_dir(ins, &["attL1"]),
            site_feature_dir(ins, &["attL2"]),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => return Err(format!(
                "The {insert_label} must contain two different att sites that face each other."
            )),
        };
        if !matches!((l1.2, l2.2), (1, 2) | (2, 1)) {
            return Err("The Entry Clone must contain attL1 and attL2 sites that face each other.".into());
        }
        let (l1s, l1e, d1) = l1;
        let (l2s, l2e, _d2) = l2;
        if l1e < l1s || l2e < l2s {
            return Err("The Entry Clone att sites are invalid.".into());
        }
        if d1 == 1 {
            // attB1 = attL1 靠近基因的末 25 bp；gene = attL1 末到 attL2 头。
            slices.push(ins.slice(l1e.saturating_sub(25), l1e, false));
            slices.push(gene_slice(ins, (l1s, l1e, 1), (l2s, l2e, 2))?);
            slices.push(ins.slice(l2s, l2s + 25, true));
            let gene_len = if l1e <= l2s {
                l2s - l1e
            } else {
                ins.seq.len() - l1e + l2s
            };
            att_offsets.push((off, gene_len));
            off += 25 + gene_len + 25;
        } else {
            // 反向 Entry Clone：gene = attL2 末到 attL1 头，整体反向互补。
            slices.push(ins.slice(l2e, l2e + 25, false));
            slices.push(ins.slice(l2e, l1s, true));
            slices.push(ins.slice(l1s.saturating_sub(25), l1s, true));
            let gene_len = if l2e <= l1s {
                l1s - l2e
            } else {
                ins.seq.len() - l2e + l1s
            };
            att_offsets.push((off, gene_len));
            off += 25 + gene_len + 25;
        }
    }
    // 目的载体后半段（attR2 之后）。
    slices.push(vector.slice(r2e, vector.seq.len(), false));

    let mut product = assemble(&slices, true, 0, 0);
    product.features.retain(|f| !is_att_name(&f.name));
    for (start, gene_len) in &att_offsets {
        product.features.push(make_site_feature("attB1", *start, 25, 1));
        product.features.push(make_site_feature("attB2", *start + 25 + *gene_len, 25, 2));
    }
    product.junctions.clear();
    Ok(product)
}

fn is_att_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    n.starts_with("attb")
        || n.starts_with("attp")
        || n.starts_with("attl")
        || n.starts_with("attr")
}

fn make_site_feature(name: &str, start: usize, len: usize, dir: i32) -> FeatureAnnotation {
    FeatureAnnotation {
        name: name.into(),
        kind: "protein_bind".into(),
        sub_type: String::new(),
        note: String::new(),
        visible: true,
        directionality: dir,
        reading_frame: 0,
        translate_first_codon_as_met: false,
        max_run_on: None,
        segments: vec![FeatureSegment {
            range: Range::new(start, start + len),
            color: "#3f8f9f".into(),
            seg_type: "standard".into(),
            translated: false,
        }],
    }
}

/// 找一对位点：优先用特征；没有特征时用 25 bp 核心 motif（支持反向链）。
fn site_pair(
    src: &Source,
    names1: &[&str],
    names2: &[&str],
    m1: &[u8],
    m2: &[u8],
) -> Option<((usize, usize, i32), (usize, usize, i32))> {
    if let (Some(a), Some(b)) = (
        site_feature_dir(src, names1),
        site_feature_dir(src, names2),
    ) {
        return Some((a, b));
    }
    if let (Some(a), Some(b)) = (find_motif(&src.seq, m1), find_motif(&src.seq, m2)) {
        return Some(((a, a + m1.len(), 1), (b, b + m2.len(), 1)));
    }
    let mut rc = src.seq.clone();
    revcomp(&mut rc);
    let a = find_motif(&rc, m1)?;
    let b = find_motif(&rc, m2)?;
    let n = src.seq.len();
    let ma = n.saturating_sub(a + m1.len());
    let mb = n.saturating_sub(b + m2.len());
    Some(((ma, ma + m1.len(), 2), (mb, mb + m2.len(), 2)))
}

/// 两个位点之间的基因片段；反向位点取反向互补。
fn gene_slice(
    src: &Source,
    b1: (usize, usize, i32),
    b2: (usize, usize, i32),
) -> Result<Slice, String> {
    let (s1, e1, d1) = b1;
    let (s2, e2, _d2) = b2;
    if d1 == 2 {
        Ok(src.slice(e2, s1, true))
    } else {
        Ok(src.slice(e1, s2, false))
    }
}

/// 把片段按顺序拼成产物：相邻片段自动识别重叠（≥ `min_overlap`）并从**前一
/// 片段末尾**去掉重复部分。`circular` 只标拓扑，不影响拼接。
pub fn assemble(
    slices: &[Slice],
    circular: bool,
    min_overlap: usize,
    max_overlap: usize,
) -> Product {
    let mut seq: Vec<u8> = Vec::new();
    let mut features: Vec<FeatureAnnotation> = Vec::new();
    let mut junctions: Vec<Junction> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::new();

    // 每段的有效范围（去掉被后一段覆盖的尾巴）
    let mut trimmed: Vec<(usize, usize)> = Vec::new(); // (start, end) in slice coords
    for i in 0..slices.len() {
        let s = 0usize;
        let mut e = slices[i].seq.len();
        if i + 1 < slices.len() && min_overlap > 0 {
            let ov = overlap_len(
                &slices[i].seq,
                &slices[i + 1].seq,
                min_overlap,
                max_overlap,
            );
            if ov > 0 {
                e = e.saturating_sub(ov);
            }
        }
        trimmed.push((s, e));
    }

    // 环状产物还要去掉「最后一段的尾巴」与「第一段开头」之间的重复同源臂。
    // 线性拼接的相邻片段已经在上面的循环里去重，但首尾相接处不会再经过一遍。
    let mut circular_junction: Option<Junction> = None;
    if circular && slices.len() > 1 && min_overlap > 0 {
        let first = &slices[0];
        let last_i = slices.len() - 1;
        let (_, last_e) = trimmed[last_i];
        let last = &slices[last_i];
        let ov = overlap_len(
            &last.seq[..last_e],
            &first.seq[..trimmed[0].1],
            min_overlap,
            max_overlap,
        );
        if ov > 0 {
            trimmed[last_i].1 = last_e.saturating_sub(ov);
            let total: usize = trimmed.iter().map(|(a, b)| b - a).sum();
            circular_junction = Some(Junction {
                at: total,
                overlap: ov,
                seq_start: total,
                seq_end: ov,
            });
        }
    }

    for (i, sl) in slices.iter().enumerate() {
        let (s, e) = trimmed[i];
        let base = seq.len();
        if i + 1 < slices.len() {
            let ov = slices[i].seq.len() - e;
            if ov > 0 {
                junctions.push(Junction {
                    at: base + (e - s),
                    overlap: ov,
                    seq_start: base + (e - s),
                    seq_end: base + (e - s) + ov,
                });
            } else {
                junctions.push(Junction {
                    at: base + (e - s),
                    overlap: 0,
                    seq_start: base + (e - s),
                    seq_end: base + (e - s),
                });
            }
        }
        spans.push((base, base + (e - s)));
        seq.extend_from_slice(&sl.seq[s..e]);
        for f in &sl.features {
            let mut segs: Vec<FeatureSegment> = Vec::new();
            for seg in &f.segments {
                let a = seg.range.start.max(s);
                let b = seg.range.end.min(e);
                if b > a {
                    segs.push(FeatureSegment {
                        range: Range::new(base + a - s, base + b - s),
                        color: seg.color.clone(),
                        ..seg.clone()
                    });
                }
            }
            if segs.is_empty() {
                continue;
            }
            let mut nf = f.clone();
            let first = segs.remove(0);
            nf.segments = {
                let mut v = vec![first];
                v.extend(segs);
                v
            };
            nf.max_run_on=f.max_run_on.filter(|end|*end>s && *end<=e).map(|end|base+end-s);
            features.push(nf);
        }
    }

    if let Some(j) = circular_junction {
        junctions.push(j);
    }
    Product {
        seq,
        features,
        junctions,
        spans,
        circular,
        ends: Default::default(),
    }
}

/// TA / GC / TOPO 克隆的选项。
#[derive(Debug, Clone, Copy)]
pub struct TaTopoOptions {
    /// 给片段补 3′ A 突出端（SnapGene 的 `Add 3' A overhangs to Insert`）。
    pub add_a: bool,
    /// 是否要求粘性末端互补（TA/GC/定向 TOPO）；平末端 TOPO 不要求。
    pub require_compatible: bool,
}

impl Default for TaTopoOptions {
    fn default() -> Self {
        Self {
            add_a: false,
            require_compatible: true,
        }
    }
}

/// 提示文案（与 SnapGene 一致）。
pub const TA_NEEDS_A: &str = "TA cloning requires 3' A overhang in Fragment.";
pub const TA_INCOMPATIBLE: &str = "Please make the ends compatible for TOPO® TA cloning.";

/// IUPAC 碱基 → 位集。
fn iupac_bases(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

/// 两段 IUPAC 序列是否可能相同（逐位交集非空）。
pub fn iupac_compatible(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| iupac_bases(*x) & iupac_bases(*y) != 0)
}

/// 把 `[a, b)` 旋转 `r`（左移）到长度 `n` 的坐标里，跨原点时拆成两段。
fn rotate_range(a: usize, b: usize, r: usize, n: usize) -> Vec<(usize, usize)> {
    if a >= b || n == 0 {
        return Vec::new();
    }
    let r = r % n;
    let na = (a + n - r) % n;
    let nb = (b + n - r) % n;
    if na < nb {
        vec![(na, nb)]
    } else {
        vec![(na, n), (0, nb)]
    }
}

/// 把一个来源的特征投影到产物：先裁到顶链区间，再平移、旋转。
fn project_rotated(
    src: &Source,
    (lo, hi): (usize, usize),
    offset: usize,
    rot: usize,
    total: usize,
) -> Vec<FeatureAnnotation> {
    let mut out = Vec::new();
    for f in &src.features {
        let mut segs: Vec<FeatureSegment> = Vec::new();
        for seg in &f.segments {
            let a = seg.range.start.max(lo);
            let b = seg.range.end.min(hi);
            if b <= a {
                continue;
            }
            let (a, b) = (offset + a - lo, offset + b - lo);
            for (ra, rb) in rotate_range(a, b, rot, total) {
                segs.push(FeatureSegment {
                    range: Range::new(ra, rb),
                    color: seg.color.clone(),
                    ..seg.clone()
                });
            }
        }
        if segs.is_empty() {
            continue;
        }
        let mut nf = f.clone();
        let first = segs.remove(0);
        nf.segments = {
            let mut v = vec![first];
            v.extend(segs);
            v
        };
        nf.max_run_on=f.max_run_on.filter(|end|*end>lo && *end<=hi)
            .map(|end|(offset+end-lo-1+total-rot%total)%total+1);
        out.push(nf);
    }
    out
}

/// 给平末端片段补 3′ A 突出端。
///
/// 存储序列前置底链 A 尾的互补标记 T，末尾追加顶链 A 尾。
/// 特征随前置标记平移一位；已有粘性末端时原样返回。
pub fn add_a_overhangs(src: &Source) -> Source {
    add_three_prime_overhangs(src, b'A')
}

/// Add an explicit 3′ A or G base to each strand, retaining the native union.
pub fn add_three_prime_overhangs(src: &Source, base: u8) -> Source {
    if src.ends.upstream_stickiness != 0 || src.ends.downstream_stickiness != 0 {
        return src.clone();
    }
    let mut out = src.clone();
    let mut seq = Vec::with_capacity(src.seq.len() + 2);
    seq.push(crate::complement(base));
    seq.extend_from_slice(&src.seq);
    seq.push(base);
    out.seq = seq;
    for feature in &mut out.features {
        feature.max_run_on=feature.max_run_on.map(|end|end+1);
        for segment in &mut feature.segments {
            segment.range.start += 1;
            segment.range.end += 1;
        }
    }
    out.ends.upstream_stickiness = -1;
    out.ends.downstream_stickiness = -1;
    out.ends.original_first_base = None;
    out
}

/// Method-specific end rules. The legacy TA options remain available for its
/// established native-file oracle; other TOPO branches must not bypass checks.
pub fn end_cloning_product(vector: &Source, insert: &Source, method: Method, add_tail: bool) -> Result<Product, String> {
    if vector.circular || insert.circular {
        return Err("Please choose linear Vector and Fragment sequences.".into());
    }
    for source in [vector, insert] {
        let n = source.seq.len();
        let top = source.ends.top_range(n);
        let bottom = source.ends.bottom_range(n);
        if top.start.max(bottom.start) >= top.end.min(bottom.end)
            || source.seq.iter().any(|b| hx_enzymes::iupac_bits(*b).is_none())
        {
            return Err("The source must contain a valid double-stranded DNA region.".into());
        }
    }
    let mut working = insert.clone();
    match method {
        Method::TaCloning | Method::TopoTa | Method::GcCloning => {
            let gc = method == Method::GcCloning;
            if add_tail { working = add_three_prime_overhangs(&working, if gc { b'G' } else { b'A' }); }
            if vector.ends.upstream_stickiness != -1 || vector.ends.downstream_stickiness != -1 {
                return Err(if gc { "GC cloning requires a Vector with single-base 3' C overhangs." } else { "TA cloning requires a Vector with single-base 3' T overhangs." }.into());
            }
            if working.ends.upstream_stickiness != -1 || working.ends.downstream_stickiness != -1 {
                return Err(if gc { "GC cloning requires 3' G overhangs in the Fragment." } else { TA_NEEDS_A }.into());
            }
            let (vbase, ibase) = if gc { (b'C', b'G') } else { (b'T', b'A') };
            let matches = |source: &Source, base: u8| source.seq.first().zip(source.seq.last()).is_some_and(|(&a, &b)|
                iupac_compatible(&[a], &[crate::complement(base)]) && iupac_compatible(&[b], &[base]));
            if !matches(vector, vbase) || !matches(&working, ibase) {
                return Err(if gc { "Please make the ends compatible for GC cloning." } else { TA_INCOMPATIBLE }.into());
            }
            ta_topo_product(vector, &working, TaTopoOptions { add_a: false, require_compatible: true })
        }
        Method::TopoBlunt => {
            if vector.has_sticky_ends() || insert.has_sticky_ends() {
                return Err("Blunt TOPO cloning requires blunt ends in both the Vector and Fragment.".into());
            }
            ta_topo_product(vector, insert, TaTopoOptions { add_a: false, require_compatible: false })
        }
        Method::TopoDirectional => {
            if insert.has_sticky_ends() {
                return Err("Directional TOPO cloning requires a blunt-ended Fragment with a 5' CACC sequence.".into());
            }
            if vector.ends.upstream_stickiness != 0 || vector.ends.downstream_stickiness != 4
                || !vector.seq.ends_with(b"CACC") || !insert.seq.starts_with(b"CACC")
            {
                return Err("Directional TOPO cloning requires the CACC directional junction in the Vector and Fragment.".into());
            }
            ta_topo_product(vector, insert, TaTopoOptions { add_a: false, require_compatible: false })
        }
        _ => Err("Select a TA, GC or TOPO cloning method.".into()),
    }
}

/// TA / GC / TOPO 克隆：把线性化载体和带粘性末端的片段拼成环状产物。
///
/// SnapGene 的模型（用 pCR®2.1-TOPO® + 651 bp Taq PCR 产物实测）：
///
/// - 产物顶链 = 载体顶链 + 片段顶链（`4581 = 3930 + 651`）；
/// - 两个接头的粘性末端必须互补：载体存储序列末位与片段存储序列首位
///   IUPAC 相容（`T` vs `Y`），环化处片段末位与载体首位相容（`R` vs `A`）；
/// - 产物按载体的 `OriginalFirstBase` 旋转，回到原始质粒编号。
pub fn ta_topo_product(
    vector: &Source,
    insert: &Source,
    opt: TaTopoOptions,
) -> Result<Product, String> {
    let (vstart, vend) = vector.top_range();
    if vend <= vstart {
        return Err("Please choose a source file.".into());
    }

    // 算法与示意图共用同一个加 A 后的工作副本。
    let ins = if opt.add_a {
        add_a_overhangs(insert)
    } else {
        insert.clone()
    };
    let (iup, idown) = (ins.ends.upstream_stickiness, ins.ends.downstream_stickiness);
    let (istart, iend) = ins.top_range();
    if iend <= istart {
        return Err("Please choose a source file.".into());
    }

    // 粘性末端互补检查（TA/GC/定向 TOPO）。
    if opt.require_compatible {
        if idown >= 0 {
            return Err(TA_NEEDS_A.into());
        }
        // 载体 → 片段
        let k = (-vector.ends.downstream_stickiness).max(-iup).max(1) as usize;
        if vector.seq.len() >= k
            && ins.seq.len() >= k
            && !iupac_compatible(&vector.seq[vector.seq.len() - k..], &ins.seq[..k])
        {
            return Err(TA_INCOMPATIBLE.into());
        }
        // 片段 → 载体（环化）
        let k2 = (-ins.ends.downstream_stickiness)
            .max(-vector.ends.upstream_stickiness)
            .max(1) as usize;
        if ins.seq.len() >= k2
            && vector.seq.len() >= k2
            && !iupac_compatible(&ins.seq[ins.seq.len() - k2..], &vector.seq[..k2])
        {
            return Err(TA_INCOMPATIBLE.into());
        }
    }

    let vt = &vector.seq[vstart..vend];
    let it = &ins.seq[istart..iend];
    let mut cat = Vec::with_capacity(vt.len() + it.len());
    cat.extend_from_slice(vt);
    cat.extend_from_slice(it);
    // 环化处：片段末位的简并碱基按载体首位解析（oracle：R → A）。
    if let (Some(&first), Some(&last)) = (cat.first(), cat.last()) {
        if iupac_bases(last).count_ones() > 1
            && vector.ends.upstream_stickiness < 0 && ins.ends.downstream_stickiness < 0
            && iupac_bases(first).count_ones() == 1
            && iupac_bases(last) & iupac_bases(first) != 0
        {
            let n = cat.len();
            cat[n - 1] = first;
        }
    }
    let total = cat.len();
    let rot = vector.rotation_offset(vstart, total);

    let mut seq = Vec::with_capacity(total);
    seq.extend_from_slice(&cat[rot..]);
    seq.extend_from_slice(&cat[..rot]);

    let mut features =
        project_rotated(vector, (vstart, vend), 0, rot, total);
    features.extend(project_rotated(&ins, (istart, iend), vt.len(), rot, total));

    let insert_start = (vt.len() + total - rot % total) % total;
    let insert_end = (vt.len() + it.len() + total - rot % total) % total;
    let spans = if insert_start < insert_end {
        vec![(insert_start, insert_end)]
    } else {
        vec![(insert_start, total), (0, insert_end)]
    };

    Ok(Product {
        seq,
        features,
        junctions: vec![Junction {
            at: insert_start,
            overlap: vector.ends.downstream_stickiness.unsigned_abs() as usize,
            seq_start: insert_start,
            seq_end: insert_start,
        }],
        spans,
        circular: true,
        ends: Default::default(),
    })
}

/// 产物 → `.dna` 字节（New File 骨架 + 特征）。
pub fn product_bytes(prod: &Product, name: &str) -> Vec<u8> {
    let flags = MoleculeFlags::from_byte(
        MoleculeFlags::DOUBLE_STRANDED
            | if prod.circular {
                MoleculeFlags::CIRCULAR
            } else {
                0
            },
    );
    let bytes = NewFile::dna(prod.seq.clone())
        .with_flags(flags)
        .with_name(name)
        .to_bytes();
    let mut file = hx_formats::DnaFile::parse(&bytes).expect("NewFile generates valid DNA");
    if let Some(block)=file.blocks.iter_mut().find(|b|b.block_type==10) {
        block.payload=hx_formats::write::feature_annotations_xml(&prod.features,prod.seq.len());
    }
    if let Some(block) = file.blocks.iter_mut().find(|b| b.block_type == 8) {
        let e = &prod.ends;
        let original = e.original_first_base.map(|base| format!("<OriginalFirstBase>{base}</OriginalFirstBase>"))
            .unwrap_or_default();
        let xml = format!(
            "<AdditionalSequenceProperties><UpstreamStickiness>{}</UpstreamStickiness>\
             <DownstreamStickiness>{}</DownstreamStickiness>\
             <UpstreamModification>{}</UpstreamModification>\
             <DownstreamModification>{}</DownstreamModification>{original}</AdditionalSequenceProperties>",
            e.upstream_stickiness, e.downstream_stickiness,
            hx_formats::util::xml_escape(&e.upstream_modification),
            hx_formats::util::xml_escape(&e.downstream_modification),
        );
        *block = hx_formats::Block::new(8, xml.into_bytes());
    }
    file.to_bytes()
}

fn revcomp(seq: &mut Vec<u8>) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        };
    }
}

fn revcomp_str(s: &[u8]) -> Vec<u8> {
    let mut v = s.to_vec();
    revcomp(&mut v);
    v
}

/// 反向互补后的特征：坐标翻转、方向翻转。
fn reverse_feature(f: &mut FeatureAnnotation, len: usize) {
    f.max_run_on=f.max_run_on.filter(|end|*end>0 && *end<=len).map(|end|len-end+1);
    for seg in &mut f.segments {
        let (a, b) = (seg.range.start, seg.range.end);
        seg.range = Range::new(len.saturating_sub(b), len.saturating_sub(a));
    }
    // 多段特征的段顺序也要反过来（SnapGene 用段顺序表示 5′→3′）
    f.segments.reverse();
    f.directionality = match f.directionality {
        1 => 2,
        2 => 1,
        other => other,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    fn src(seq: &str, circular: bool) -> Source {
        Source::new("test.dna", seq.as_bytes().to_vec(), circular)
    }

    #[test]
    fn slice_takes_region_and_reverse_complement() {
        let s = src("AAAACCCCGGGGTTTT", false);
        assert_eq!(s.slice(4, 8, false).seq, b"CCCC");
        assert_eq!(s.slice(4, 8, true).seq, b"GGGG");
    }

    #[test]
    fn circular_slice_wraps_origin() {
        let s = src("AAAACCCCGGGGTTTT", true);
        assert_eq!(s.slice(14, 2, false).seq, b"TTAA");
    }

    #[test]
    fn ta_topo_concatenates_top_strands_and_rotates() {
        // 载体：上游 3′ 突出端标记 A + 顶链 GGGCCCTT（下游 3′ T 突出端）。
        let mut vector = src("AGGGCCCTT", false);
        vector.ends.upstream_stickiness = -1;
        vector.ends.downstream_stickiness = -1;
        vector.ends.original_first_base = Some(7);
        // 片段：Taq 产物 Y + ACAAGGT + R。
        let mut insert = src("YACAAGGTR", false);
        insert.ends.upstream_stickiness = -1;
        insert.ends.downstream_stickiness = -1;

        let p = ta_topo_product(&vector, &insert, TaTopoOptions::default()).unwrap();
        // 顶链 8 + 8；环化处的 R 按载体首位（这里是 G）解析。
        assert_eq!(p.seq.len(), 16);
        let unrot = b"GGGCCCTTACAAGGTG";
        let rot = 7 - 1;
        let expect: Vec<u8> = unrot[rot..].iter().chain(unrot[..rot].iter()).copied().collect();
        assert_eq!(p.seq, expect);
        assert!(p.circular);
    }

    #[test]
    fn ta_topo_rejects_blunt_insert() {
        let mut vector = src("AGGGCCCTT", false);
        vector.ends.upstream_stickiness = -1;
        vector.ends.downstream_stickiness = -1;
        let insert = src("ACAAGGT", false);
        assert_eq!(
            ta_topo_product(&vector, &insert, TaTopoOptions::default()).unwrap_err(),
            TA_NEEDS_A
        );
        // 补 3′ A 之后可以克隆。
        let a = add_a_overhangs(&insert);
        assert_eq!(a.seq, b"TACAAGGTA");
        assert_eq!(a.display_len(), 8);
        let p = ta_topo_product(
            &vector,
            &insert,
            TaTopoOptions {
                add_a: true,
                require_compatible: true,
            },
        )
        .unwrap();
        assert_eq!(p.seq.len(), 16);
    }

    #[test]
    fn overlap_detection() {
        assert_eq!(overlap_len(b"AAAACCCC", b"CCCCGGGG", 3, 40), 4);
        assert_eq!(overlap_len(b"AAAACCCC", b"CCCCGGGG", 5, 40), 0);
        assert_eq!(overlap_len(b"ACGT", b"ACGT", 4, 40), 4);
    }

    #[test]
    fn assemble_trims_overlap_and_keeps_one_copy() {
        let a = src("GGGGAAAACCCC", false).slice(0, 12, false);
        let b = src("CCCCGGGGTTTT", false).slice(0, 12, false);
        // 重叠 CCCC(4) → 产物去掉 a 末尾的 CCCC
        let p = assemble(&[a, b], true, 4, 40);
        assert_eq!(p.seq, b"GGGGAAAACCCCGGGGTTTT");
        assert_eq!(p.junctions[0].overlap, 4);
    }

    #[test]
    fn features_follow_slices_and_flip() {
        let mut s = src("AAAACCCCGGGGTTTT", false);
        s.features = vec![FeatureAnnotation {
            name: "cds".into(),
            kind: "CDS".into(),
            sub_type: String::new(),
            note: String::new(),
            visible: true,
            directionality: 1,
            reading_frame: 1,
            translate_first_codon_as_met: false,
            max_run_on: None,
            segments: vec![FeatureSegment {
                range: Range::new(4, 8),
                color: "#ff0000".into(),
                seg_type: "standard".into(),
                translated: false,
            }],
        }];
        let tailed = add_a_overhangs(&s);
        assert_eq!(tailed.features[0].segments[0].range, Range::new(5, 9));
        assert_eq!(s.features[0].segments[0].range, Range::new(4, 8));
        assert_eq!(add_a_overhangs(&tailed).seq, tailed.seq);
        let fwd = s.slice(2, 10, false);
        assert_eq!(fwd.features[0].segments[0].range, Range::new(2, 6));
        let rev = s.slice(2, 10, true);
        assert_eq!(rev.features[0].segments[0].range, Range::new(2, 6));
        assert_eq!(rev.features[0].directionality, 2);
    }

    #[test]
    fn circular_full_circle_slice_when_start_equals_end() {
        let s = src("AAAACCCCGGGGTTTT", true);
        assert_eq!(s.slice(5, 5, false).seq.len(), 16);
        assert_eq!(s.slice(5, 5, true).seq.len(), 16);
    }

    #[test]
    fn overlap_is_capped() {
        let s = "A".repeat(100);
        let a = src(&s, false).slice(0, 100, false);
        let b = src(&s, false).slice(0, 100, false);
        let p = assemble(&[a, b], true, 15, 40);
        assert_eq!(p.junctions[0].overlap, 40);
        // 环状首尾也同源：两个 100 bp 片段各被裁 40 bp，产物 120 bp。
        assert_eq!(p.seq.len(), 120);
    }

    #[test]
    fn product_bytes_parse_round_trip() {
        let a = src("ATGAAACCCGGGTTT", false).slice(0, 15, false);
        let p = assemble(&[a], false, 0, 0);
        let bytes = product_bytes(&p, "Cloned");
        let f = hx_formats::dna::DnaFile::parse(&bytes).unwrap();
        assert_eq!(f.sequence().unwrap().sequence, p.seq);
    }

    #[test]
    fn cloned_feature_metadata_survives_native_save_and_reload() {
        let xml=br##"<Features><Feature name="CDS &amp; tag" type="CDS" subType="active" visible="0" directionality="3" readingFrame="-2" translateFirstCodonAsMet="1" maxRunOn="36"><Segment range="4-12" color="#aabbcc" type="standard" translated="1"/><Segment range="13-15" color="#112233" type="gap"/><Q name="note"><V text="keep &lt;tag&gt; &amp; both colors"/></Q></Feature></Features>"##;
        let features=hx_formats::annot::parse_features(xml);
        let p=Product{seq:b"ATG".repeat(20),features:features.clone(),junctions:vec![],spans:vec![],circular:true,ends:Default::default()};
        let bytes=product_bytes(&p,"metadata-check");
        let file=hx_formats::DnaFile::parse(&bytes).unwrap();
        let loaded=hx_formats::annot::parse_features(&file.blocks.iter().find(|b|b.block_type==10).unwrap().payload);
        assert_eq!(loaded,features,"cloning must not reduce annotations to name/type/range");
        assert_eq!(file.sequence().unwrap().sequence,p.seq);
        assert_eq!(file.end_properties(),p.ends);
    }

    #[test]
    fn cloned_run_on_endpoint_uses_the_same_coordinate_transforms_as_features() {
        let mut feature=make_site_feature("coding",3,12,1);
        feature.max_run_on=Some(18);
        feature.reading_frame=2;feature.segments[0].translated=true;
        let source=src("ACGTACGTACGTACGTACGTACGT",false).with_features(vec![feature]);
        let mut slice=source.slice(3,21,false);
        assert_eq!(slice.features[0].max_run_on,Some(15));
        slice.prepend(b"GGG");assert_eq!(slice.features[0].max_run_on,Some(18));
        let assembled=assemble(&[src("CCC",false).slice(0,3,false),slice.clone()],false,0,0);
        assert_eq!(assembled.features[0].max_run_on,Some(21));
        let reverse=source.flipped();assert_eq!(reverse.features[0].max_run_on,Some(7));
        assert_eq!(reverse.flipped().features[0],source.features[0]);
        let ta=add_a_overhangs(&source);assert_eq!(ta.features[0].max_run_on,Some(19));
        slice.truncate_end(6);assert_eq!(slice.features[0].max_run_on,None);
    }

    /// pDONR221 的 attP1/attP2 正向序列（232 bp）。
    const ATTP1: &[u8] = b"AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA";
    const ATTP2: &[u8] = b"TACAGGTCACTAATACCATCTAAGTAGTTGATTCATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTT";

    #[test]
    fn gateway_bp_uses_true_hybrid_attl_sites() {
        // 供体：4 bp 前缀 + attP1 + 6 bp 间隔 + attP2 + "GTTT"。
        // 前缀末位 C 与 pDONR221 attP1 前 1 bp 相同；末位 G 正是 pDONR221
        // attP2 之后的那 1 bp，attL2 的 78 bp 供体段含它。
        let mut donor_seq = b"GGGC".to_vec();
        donor_seq.extend_from_slice(ATTP1);
        donor_seq.extend_from_slice(b"CCCCCC");
        donor_seq.extend_from_slice(ATTP2);
        donor_seq.extend_from_slice(b"GTTT");
        assert_eq!(donor_seq.len(), 478);
        let donor = Source::new("donor.dna", donor_seq, true).with_features(vec![
            make_site_feature("attP1", 4, ATTP1.len(), 1),
            make_site_feature("attP2", 4 + ATTP1.len() + 6, ATTP2.len(), 2),
        ]);

        let mut insert_seq = ATT_B1.to_vec();
        insert_seq.extend_from_slice(b"AAA");
        insert_seq.extend_from_slice(ATT_B2);
        let insert = Source::new("attB.dna", insert_seq, false);
        let p = gateway_recombine(&[&insert], &donor, "attB Fragment", true).unwrap();

        let mut want = b"GGG".to_vec();
        want.extend_from_slice(ATT_L1);
        want.extend_from_slice(b"AAA");
        want.extend_from_slice(ATT_L2);
        want.extend_from_slice(b"TTT");
        assert_eq!(p.seq, want);
        assert_eq!(p.seq.len(), 209);

        let l1 = p.features.iter().find(|f| f.name == "attL1").unwrap();
        assert_eq!(l1.directionality, 1);
        assert_eq!(l1.segments[0].range, Range::new(3, 103));
        let l2 = p.features.iter().find(|f| f.name == "attL2").unwrap();
        assert_eq!(l2.directionality, 2);
        assert_eq!(l2.segments[0].range, Range::new(106, 206));
    }

    #[test]
    fn gateway_replaces_the_cassette_between_sites() {
        let mut insert_seq = ATT_B1.to_vec();
        insert_seq.extend_from_slice(b"AAA");
        insert_seq.extend_from_slice(ATT_B2);
        let mut vector_seq = b"GGG".to_vec();
        vector_seq.extend_from_slice(ATT_B1);
        vector_seq.extend_from_slice(b"CCCCCC");
        vector_seq.extend_from_slice(ATT_B2);
        vector_seq.extend_from_slice(b"TTT");
        let insert = Source::new("attB.dna", insert_seq.clone(), false);
        let vector = Source::new("donor.dna", vector_seq, true);
        assert!(gateway_recombine(&[&insert], &vector, "attB Fragment", true).is_err());
        assert!(gateway_recombine(&[&insert], &vector, "Entry Clone", false).unwrap_err().contains("attR1"));
    }
}
