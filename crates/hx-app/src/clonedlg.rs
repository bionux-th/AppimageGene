//! 「分子克隆」对话框（SnapGene `Actions ▸ … Cloning`）。
//!
//! 布局与页签照真 SnapGene 8.2.3 实拍（见 `docs/cloning-2026-10-22.zh-CN.md`）：
//!
//! - 限制性/Gibson/Golden Gate/In-Fusion/NEBuilder：`Vector | Fragment | Product`
//! - TA/GC/TOPO：`Fragment | Linearized Vector | Product`
//! - Gateway BP：`attB Fragment | Donor Vector | Entry Clone`
//! - Gateway LR：`Entry Clone(s) | Destination Vector | Expression Clone`
//! - Gateway BP+LR：前五者合并
//! - 线性连接：`Fragments (2) | Product`
//!
//! 左栏是预览（Map/Sequence/Enzymes/Primers 子页签）、PCR 引物面板和底部的
//! 拼接示意图；右栏是来源/酶/方向等参数、黄色提示、产物命名与 Clone/Cancel。
//!
//! 计算全部在 [`hx_algo::cloning`]；本文件只收参数、画图、给提示。

use std::path::{Path, PathBuf};

use eframe::egui::{self, Color32, Pos2, Rect, RichText, Stroke, Vec2};
use hx_algo::cloning::{self, Method, Product, Slice, Source};

use crate::theme;

mod golden_gate;
mod output;
#[cfg(test)]
mod success_matrix;
#[cfg(test)]
mod multisite_tests;

// ---------------------------------------------------------------------------
// 数据
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linearize {
    Pcr,
    Restriction,
    Direct,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoldenGateMode { PcrDigest, Digest, OtherEnzymes, Direct }

/// 一侧（载体 / 片段 / 供体 / 目的地）的选择。
#[derive(Debug, Clone)]
pub struct Pick {
    /// 文档下标（`CloningDlg::docs`）。
    pub doc: Option<usize>,
    /// 商品化标准载体的名字（选中后按需从 SnapGene 资源目录读序列）。
    pub std_name: Option<String>,
    /// 标准载体页签里选的是「商品化」而不是自定义文件。
    pub use_std: bool,
    pub std_src: Option<Source>,
    /// 反向使用（SnapGene 的 `→` / `←`）。
    pub reversed: bool,
    /// 限制性酶切克隆的酶（`hx_enzymes::ENZYMES` 下标）。
    pub enzyme: Option<usize>,
    pub enzyme2: Option<usize>,
    pub blunt: bool,
    /// 已选双链切口；None 时按骨架最长/插入片段最短的完整酶切片段初始化。
    pub cut_pair: Option<(cloning::RestrictionCut, cloning::RestrictionCut)>,
    /// 线性化方式。
    pub lin: Linearize,
    /// 线性化用的限制性酶。
    pub lin_enzyme: Option<usize>,
    /// 只用到序列的一段（0-based 半开）；`None` = 整条。
    pub region: Option<(usize, usize)>,
    /// `Use as a template for PCR`（否则 `Use directly`）。
    pub use_pcr: bool,
    /// 引物 5′ 尾巴（正向链）。
    pub primer1: String,
    pub primer2: String,
    pub phos1: bool,
    pub phos2: bool,
    /// `Add 3' A overhangs to Insert`（TA/GC/TOPO 的平末端片段补救）。
    pub add_a: bool,
    pub gg_mode: GoldenGateMode,
}

impl Default for Pick {
    fn default() -> Self {
        Self {
            doc: None,
            std_name: None,
            use_std: false,
            std_src: None,
            reversed: false,
            enzyme: None,
            enzyme2: None,
            blunt: false,
            cut_pair: None,
            lin: Linearize::Pcr,
            lin_enzyme: None,
            region: None,
            use_pcr: true,
            primer1: String::new(),
            primer2: String::new(),
            phos1: false,
            phos2: false,
            add_a: false,
            gg_mode: GoldenGateMode::PcrDigest,
        }
    }
}

pub struct DocSnap {
    pub name: String,
    pub src: Source,
    /// 打开时的 `.dna` 字节（含引物/特征）。克隆预览复用主视图时需要完整的
    /// `Document`；`Source` 只有序列/特征/引物，没有块级信息。
    pub bytes: Option<Vec<u8>>,
}

impl DocSnap {
    pub fn new(name: impl Into<String>, src: Source) -> Self {
        Self {
            name: name.into(),
            src,
            bytes: None,
        }
    }

    pub fn with_bytes(mut self, bytes: Vec<u8>) -> Self {
        self.bytes = Some(bytes);
        self
    }
}

#[derive(Debug, Clone)]
pub struct StdEntry {
    pub name: String,
    pub path: PathBuf,
}

/// 预览区的小页签。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreviewTab {
    Map,
    Sequence,
    Enzymes,
    Primers,
}

/// 当前页签的语义。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TabKind {
    Vector,
    Fragment,
    Fragments,
    LinearizedVector,
    DonorVector,
    DestinationVector,
    Product,
}

pub struct CloningDlg {
    pub method: Method,
    /// 打开时使用的克隆默认参数（SnapGene `[actions]`）。
    pub settings: crate::settings::CloningSettings,
    pub docs: Vec<DocSnap>,
    pub vector: Pick,
    /// BP+LR 时的 Destination Vector（BP 的 Donor 用 `vector`）。
    pub destination: Pick,
    /// MultiSite BP+LR has one independently selected donor for each PCR fragment.
    pub gateway_donors: Vec<Pick>,
    pub fragments: Vec<Pick>,
    /// 多片段时正在编辑的片段下标。
    pub frag_index: usize,
    pub tab: usize,
    pub preview_tab: PreviewTab,
    pub polymerase: usize,
    /// 全局 Golden Gate Type IIS 酶，影响两侧处理和产物计算。
    pub gg_enzyme: usize,
    /// Golden Gate Settings 勾选的酶；`gg_enzyme` 是这里的下标。
    pub gg_enzymes: Vec<&'static str>,
    pub gg_adjust_open: bool,
    pub gg_overhangs: Vec<String>,
    pub product_name: String,
    pub create_product: bool,
    /// SnapGene 的 `and close this window`：取消勾选时克隆完保留对话框。
    pub close_after_clone: bool,
    pub strain: usize,
    pub hint: String,
    pub error: Option<String>,
    pub preview: Option<Product>,
    preview_key: String,
    /// 酶选择下拉的过滤文本。
    pub enzyme_filter: String,
    /// 商品化载体列表。
    pub std_donor: Vec<StdEntry>,
    pub std_dest: Vec<StdEntry>,
    pub std_ta: Vec<StdEntry>,
    /// `Use attB1 and attB2 adapter primers`。
    pub adapter_primers: bool,
    /// 打开后的重绘帧数：字体图集可能因新字号扩容，需要再排一两帧。
    frames_open: u8,
    /// 克隆对话框里独立的主视图（Map/Sequence/Enzymes/Primers）。
    pub preview_app: Option<Box<crate::HelixApp>>,
    /// 当前预览文档的签名：文档/页签/参数变化时才重建。
    preview_doc_key: String,
    /// 主视图正在显示的文档（`take` 出来交给 `preview_app`）。
    preview_doc: Option<crate::doc::Document>,
}

pub const POLYMERASES: &[&str] = &[
    "Creates Blunt Ends",
    "T4 DNA Polymerase",
    "DNA Polymerase I, Large (Klenow) Fragment",
];

pub const TA_POLYMERASES: &[&str] = &["Creates 3' A or G Overhangs", "Taq DNA Polymerase"];

pub const STRAINS: &[&str] = &[
    "Unspecified",
    "DH5alpha",
    "TOP10",
    "Stbl3",
    "JM109",
    "XL1-Blue",
    "BL21(DE3)",
    "ccdB Survival",
    "DB3.1",
];

// Gateway 商品化载体的显示名。序列在打开对话框时从 SnapGene 资源目录按文件名读取。
const DONOR_FILES: &[&str] = &[
    "pDONR221 P5-P2",
    "pDONR P2r-P3",
    "pDONR P4-P1r",
    "pDONR Zeo",
    "pDONR201",
    "pDONR207",
    "pDONR221",
    "pDONR221 P1-P4",
    "pDONR221 P1-P5r",
    "pDONR221 P3-P2",
    "pDONR221 P4r-P3r",
    "pDONR221 P5-P4",
    "pDONR222",
    "pDONR223",
];

const DEST_FILES: &[&str] = &[
    "pAd BLOCK-iT-DEST",
    "pAd CMV V5-DEST",
    "pAd PL-DEST",
    "pBAD-DEST49",
    "pcDNA-DEST40",
    "pcDNA-DEST47",
    "pcDNA-DEST53",
    "pcDNA3.1 nV5-DEST",
    "pcDNA3.2 V5-DEST",
    "pcDNA6.2 C-EmGFP-DEST",
    "pcDNA6.2 C-YFP-DEST",
    "pcDNA6.2 N-EmGFP-DEST",
    "pcDNA6.2 N-YFP-DEST",
    "pcDNA6.2 V5-DEST",
    "pDEST R4-R3",
    "pDEST R4-R3 Vector II",
    "pDEST8",
    "pDEST10",
    "pDEST14",
    "pDEST15",
    "pDEST17",
    "pDEST20",
    "pDEST22",
    "pDEST24",
];

/// 菜单「Donor Vectors ▸」列表。
pub const DONOR_NAMES: &[&str] = &[
    "pDONR221 P5-P2",
    "pDONR221 P2r-P3",
    "pDONR221 P4-P1r",
    "pDONR Zeo",
    "pDONR201",
    "pDONR207",
    "pDONR221",
    "pDONR221 P1-P4",
    "pDONR221 P1-P5r",
    "pDONR221 P3-P2",
    "pDONR221 P4r-P3r",
    "pDONR221 P5-P4",
    "pDONR222",
    "pDONR223",
];

/// 菜单「Destination Vectors ▸」列表。
pub const DEST_NAMES: &[&str] = &[
    "pAd BLOCK-iT-DEST",
    "pAd CMV V5-DEST",
    "pAd PL-DEST",
    "pBAD-DEST49",
    "pcDNA-DEST40",
    "pcDNA-DEST47",
    "pcDNA-DEST53",
    "pcDNA3.1 nV5-DEST",
    "pcDNA3.2 V5-DEST",
    "pcDNA6.2 C-EmGFP-DEST",
    "pcDNA6.2 N-EmGFP-DEST",
    "pcDNA6.2 V5-DEST",
    "pDEST R4-R3",
    "pDEST8",
    "pDEST10",
    "pDEST14",
    "pDEST15",
    "pDEST17",
    "pDEST20",
    "pDEST22",
    "pDEST24",
];

pub const TA_VECTOR_NAMES: &[&str] = &[
    "Allele TA Vector (linearized)",
    "pDrive (linearized)",
    "pETBlue-1 AccepTor",
    "pGEM-T Easy (linearized)",
    "pGEM-T (linearized)",
    "pLPS-T (linearized)",
    "pLUG (linearized)",
    "pPrime (linearized)",
    "pSC-A-amp kan (linearized)",
];

pub const GC_VECTOR_NAMES: &[&str] = &[
    "pGC Blue (linearized)",
    "pSMART GC HK (linearized)",
    "pSMART GC LK (linearized)",
];

pub const TOPO_TA_VECTOR_NAMES: &[&str] = &[
    "pBAD Thio-TOPO (linearized)",
    "pBAD-TOPO (linearized)",
    "pBlue-TOPO (linearized)",
    "pcDNA3.1 CT-GFP-TOPO (linearized)",
    "pcDNA3.1 NT-GFP-TOPO (linearized)",
    "pcDNA3.1 V5-His-TOPO (linearized)",
    "pcDNA3.3-TOPO (linearized)",
    "pcDNA3.4-TOPO (linearized)",
    "pcDNA4 HisMax-TOPO (linearized)",
    "pcDNA5 FRT TO-TOPO (linearized)",
    "pcDNA6.2 C-EmGFP-GW TOPO (linearized)",
];

pub const TOPO_BLUNT_VECTOR_NAMES: &[&str] = &[
    "pCR4Blunt-TOPO (linearized)",
    "pCR-Blunt II-TOPO (linearized)",
    "pCR-XL-2-TOPO (linearized)",
    "pLPS-B (linearized)",
    "pTOP Blunt V2 (linearized)",
];

pub const TOPO_DIR_VECTOR_NAMES: &[&str] = &[
    "pBAD202 D-TOPO",
    "pcDNA3.1D V5-His-TOPO",
    "pcDNA3.2 V5 GW D-TOPO",
    "pcDNA6.2 V5 GW D-TOPO",
    "pEF5 FRT V5-D-TOPO",
];


impl CloningDlg {
    pub fn new(method: Method, docs: Vec<DocSnap>, start: Option<usize>) -> Self {
        Self::new_with_settings(method, docs, start, crate::settings::CloningSettings::default())
    }

    /// 用持久化的 SnapGene `[actions]` 默认参数打开克隆对话框。
    pub fn new_with_settings(
        method: Method,
        docs: Vec<DocSnap>,
        start: Option<usize>,
        settings: crate::settings::CloningSettings,
    ) -> Self {
        Self::new_with_settings_and_gg(
            method,
            docs,
            start,
            settings,
            &crate::ggsettings::GgSettings::default(),
        )
    }

    /// 同时使用 Golden Gate Settings 的酶清单打开对话框。
    pub fn new_with_settings_and_gg(
        method: Method,
        docs: Vec<DocSnap>,
        start: Option<usize>,
        settings: crate::settings::CloningSettings,
        gg: &crate::ggsettings::GgSettings,
    ) -> Self {
        let n = method.fragment_count().max(1);
        // Golden Gate Settings 里勾选的酶决定下拉清单；默认酶不在清单里时
        // 退回第一个可用酶，避免对话框指向一个不存在的下标。
        let gg_enzymes = gg.enzyme_names();
        let wanted_gg_enzyme = if settings.golden_gate_default_enzyme.trim().is_empty() {
            gg.default_enzyme_name().to_string()
        } else {
            settings.golden_gate_default_enzyme.clone()
        };
        let default_gg_enzyme = gg_enzymes
            .iter()
            .position(|n| *n == wanted_gg_enzyme)
            .unwrap_or(0);
        // 载体（BP 的供体 / LR 的目的载体）不能默认成当前文档，否则
        // 片段和载体是同一个文件，产物会退化。
        let non_start: Vec<usize> = (0..docs.len())
            .filter(|j| Some(*j) != start)
            .collect();
        // 片段优先选带对应位点的文件（attB / attL），没有就退回当前文档。
        let has_sites = |i: usize, names: &[&str]| {
            names.iter().all(|n| {
                docs[i]
                    .src
                    .features
                    .iter()
                    .any(|f| f.name.eq_ignore_ascii_case(n))
            })
        };
        let name_has = |i: usize, words: &[&str]| {
            let name = docs[i].name.to_ascii_lowercase();
            words.iter().any(|w| name.contains(w))
        };
        let any_site = |i: usize, names: &[&str]| {
            docs[i]
                .src
                .features
                .iter()
                .any(|f| names.iter().any(|n| f.name.eq_ignore_ascii_case(n)))
        };
        let pure_attb = |i: usize| {
            has_sites(i, &["attB1", "attB2"])
                && !any_site(i, &["attL1", "attL2"])
                && !any_site(i, &["attR1", "attR2"])
        };
        let pure_attl = |i: usize| has_sites(i, &["attL1", "attL2"]);
        let is_bplr = matches!(method, Method::GatewayBpLr(_));
        let is_gateway = method.is_bp() || method.is_lr();
        // 供体/目的载体可能在当前页签（用户最后打开的文件），所以这里搜**全部**
        // 文档，而不是只看 non_start；否则 active=供体时会把片段自己当供体。
        let donor_doc = (0..docs.len())
            .find(|&i| has_sites(i, &["attP1", "attP2"]));
        let dest_doc = (0..docs.len())
            .find(|&i| has_sites(i, &["attR1", "attR2"]));
        let (frag_doc, vec_doc) = if is_gateway {
            // BP 找 attB 片段，LR 找 Entry Clone，BP+LR 也优先 attB 片段。
            let frag = if method.is_bp() {
                non_start
                    .iter()
                    .copied()
                    .find(|&i| name_has(i, &["attb"]))
                    .or_else(|| start.filter(|&i| pure_attb(i)))
                    .or_else(|| non_start.iter().copied().find(|&i| pure_attb(i)))
                    .or(start)
                    .or_else(|| docs.first().map(|_| 0))
            } else {
                start
                    .filter(|&i| pure_attl(i))
                    .or_else(|| non_start.iter().copied().find(|&i| name_has(i, &["entry"])))
                    .or_else(|| non_start.iter().copied().find(|&i| pure_attl(i)))
                    .or(start)
                    .or_else(|| docs.first().map(|_| 0))
            };
            let vec = if method.is_bp() {
                donor_doc
            } else {
                dest_doc
            };
            (frag, vec)
        } else if method.is_ta_topo() {
            // TA/TOPO：当前文档是 PCR 片段，载体用标准 T 载体或另一个文件。
            let frag = start
                .or_else(|| non_start.first().copied())
                .or_else(|| docs.first().map(|_| 0));
            let vec = non_start
                .iter()
                .copied()
                .find(|&i| name_has(i, &["t-vector", "pmd19", "topo", "pcr", "togo"]))
                .or_else(|| non_start.first().copied());
            (frag, vec)
        } else {
            // Gibson/In-Fusion/NEBuilder/酶切/连接：当前文档是载体，另一个是片段。
            let vec = start
                .or_else(|| non_start.first().copied())
                .or_else(|| docs.first().map(|_| 0));
            let frag = non_start.first().copied().or(vec);
            (frag, vec)
        };
        let dest_pick = if is_bplr {
            dest_doc.or(vec_doc)
        } else if docs.len() > 1 {
            vec_doc.or(Some(1))
        } else {
            None
        };
        let mut fragments = Vec::with_capacity(n);
        for i in 0..n {
            let doc = if i == 0 {
                frag_doc
            } else {
                non_start
                    .iter()
                    .copied()
                    .find(|&j| Some(j) != vec_doc && Some(j) != frag_doc)
                    .or(vec_doc)
                    .or(frag_doc)
            };
            fragments.push(Pick {
                doc: if matches!(method, Method::GoldenGate | Method::GoldenGateMulti
                    | Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble
                    | Method::InFusion | Method::InFusionMulti | Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble)
                    { None } else { doc },
                ..Default::default()
            });
        }
        // TA/GC/TOPO：SnapGene 默认 `Use directly as the Fragment`
        // （片段已经是 PCR 产物，不再当作模板）。
        if method.is_ta_topo() {
            for p in fragments.iter_mut() {
                p.use_pcr = false;
            }
        }
        // Gateway BP：片段本身已经带 attB1/attB2 时，SnapGene 默认“直接使用”
        // 而不是再做一轮 PCR（oracle：attb.dna 打开 BP 对话框时 direct 选中）。
        if method.is_bp() {
            if let Some(p) = fragments.first_mut() {
                if let Some(i) = p.doc {
                    let s = &docs[i].src.seq;
                    let fwd = has_motif(s, cloning::ATT_B1) && has_motif(s, cloning::ATT_B2);
                    let rc = revcomp_bases(s);
                    let rev = has_motif(&rc, cloning::ATT_B1) && has_motif(&rc, cloning::ATT_B2);
                    p.use_pcr = !(fwd || rev);
                }
            }
        }
        let mut dlg = Self {
            method,
            settings,
            tab: 0,
            preview_tab: PreviewTab::Map,
            vector: Pick {
                // TA/GC/TOPO 的载体通常是商品化 T 载体：没有另一个文档时不要
                // 回退到当前片段本身，交给用户在右侧选标准载体。
                doc: if method.is_ta_topo() || method.is_gateway() {
                    vec_doc
                } else {
                    vec_doc.or(if docs.is_empty() { None } else { Some(0) })
                },
                ..Default::default()
            },
            destination: Pick {
                doc: dest_pick,
                ..Default::default()
            },
            gateway_donors: vec![Pick::default();n],
            fragments,
            frag_index: 0,
            polymerase: 0,
            gg_enzyme: default_gg_enzyme,
            gg_enzymes,
            gg_adjust_open: false,
            gg_overhangs: vec![String::new(); n + 1],
            product_name: method.product_name().to_string(),
            create_product: true,
            close_after_clone: true,
            strain: 0,
            hint: String::new(),
            error: None,
            preview: None,
            preview_key: String::new(),
            enzyme_filter: String::new(),
            std_donor: Vec::new(),
            std_dest: Vec::new(),
            std_ta: Vec::new(),
            adapter_primers: false,
            frames_open: 0,
            preview_app: None,
            preview_doc_key: String::new(),
            preview_doc: None,
            docs,
        };
        if matches!(method,Method::GatewayBpLr(n) if n>1) {
            for (i,&expected) in cloning::gateway::layout(n).unwrap().iter().enumerate() {
                dlg.gateway_donors[i].doc=dlg.docs.iter().position(|d|cloning::gateway::donor_matches(&d.src,expected));
                dlg.fragments[i].doc=None;
            }
        }
        if golden_gate::active(&dlg) {
            let enzyme = dlg.gg_enzymes.get(dlg.gg_enzyme).and_then(|n| hx_enzymes::by_name(n));
            for pick in std::iter::once(&mut dlg.vector).chain(dlg.fragments.iter_mut()) {
                if pick.doc.and_then(|i| dlg.docs.get(i)).zip(enzyme).is_some_and(|(d, e)|
                    !cloning::golden_gate::fragments(&d.src, e).is_empty()) {
                    pick.gg_mode = GoldenGateMode::Digest;
                }
            }
        }
        dlg.select_gateway_lr_default_sources();
        dlg.load_std_lists();
        dlg.refresh();
        dlg
    }

    pub fn title(&self) -> String {
        match self.method {
            Method::Restriction | Method::RestrictionMulti => {
                format!("{}: {}", self.method.title(), self.method.action())
            }
            Method::Gibson
            | Method::GibsonMulti
            | Method::GibsonAssemble
            | Method::InFusion
            | Method::InFusionMulti
            | Method::NeBuilder
            | Method::NeBuilderMulti
            | Method::NeBuilderAssemble
            | Method::GoldenGate
            | Method::GoldenGateMulti => {
                format!("{}: {}", self.method.title(), self.method.action())
            }
            _ => self.method.title(),
        }
    }

    /// 页签标题（按 SnapGene 的方法差异）。
    pub fn tabs(&self) -> Vec<String> {
        match self.method {
            Method::GatewayBp => vec![
                "attB Fragment".into(),
                "Donor Vector".into(),
                "Entry Clone".into(),
            ],
            Method::GatewayLr(n) => vec![
                if n > 1 {
                    "Entry Clones".into()
                } else {
                    "Entry Clone".into()
                },
                "Destination Vector".into(),
                "Expression Clone".into(),
            ],
            Method::GatewayBpLr(n) => vec![
                if n>1 {"attB Fragments".into()}else{"attB Fragment".into()},
                if n>1 {"Donor Vectors".into()}else{"Donor Vector".into()},
                if n>1 {"Entry Clones".into()}else{"Entry Clone".into()},
                "Destination Vector".into(),
                "Expression Clone".into(),
            ],
            Method::TaCloning
            | Method::GcCloning
            | Method::TopoTa
            | Method::TopoBlunt
            | Method::TopoDirectional => vec![
                "Fragment".into(),
                "Linearized Vector".into(),
                "Product".into(),
            ],
            Method::Ligation => vec![format!("Fragments ({})", self.fragments.len()), "Product".into()],
            _ if self.method.is_multiple() => vec![
                "Vector".into(),
                format!("Fragments ({})", self.fragments.len()),
                "Product".into(),
            ],
            _ => vec!["Vector".into(), "Fragment".into(), "Product".into()],
        }
    }

    fn tab_kind(&self, idx: usize) -> TabKind {
        match self.method {
            Method::GatewayBp => match idx {
                0 => TabKind::Fragment,   // attB
                1 => TabKind::DonorVector,
                _ => TabKind::Product,
            },
            Method::GatewayLr(_) => match idx {
                0 => TabKind::Fragments,
                1 => TabKind::DestinationVector,
                _ => TabKind::Product,
            },
            Method::GatewayBpLr(_) => match idx {
                0 => TabKind::Fragment,
                1 => TabKind::DonorVector,
                2 => TabKind::Product,
                3 => TabKind::DestinationVector,
                _ => TabKind::Product,
            },
            Method::TaCloning
            | Method::GcCloning
            | Method::TopoTa
            | Method::TopoBlunt
            | Method::TopoDirectional => match idx {
                0 => TabKind::Fragment,
                1 => TabKind::LinearizedVector,
                _ => TabKind::Product,
            },
            Method::Ligation => match idx {
                0 => TabKind::Fragments,
                _ => TabKind::Product,
            },
            _ if self.method.is_multiple() => match idx {
                0 => TabKind::Vector,
                1 => TabKind::Fragments,
                _ => TabKind::Product,
            },
            _ => match idx {
                0 => TabKind::Vector,
                1 => TabKind::Fragment,
                _ => TabKind::Product,
            },
        }
    }

    fn active_pick(&self) -> &Pick {
        match self.tab_kind(self.tab) {
            TabKind::Fragment | TabKind::Fragments => {
                self.fragments.get(self.frag_index).unwrap_or(&self.fragments[0])
            }
            TabKind::DestinationVector if matches!(self.method, Method::GatewayBpLr(_)) => {
                &self.destination
            }
            TabKind::DonorVector if matches!(self.method,Method::GatewayBpLr(n) if n>1) => {
                &self.gateway_donors[self.frag_index.min(self.gateway_donors.len()-1)]
            }
            _ => &self.vector,
        }
    }

    fn select_gateway_lr_default_sources(&mut self) {
        if !matches!(self.method,Method::GatewayLr(n) if n>1){return;}
        let mut used=Vec::new();
        for (i,&(a,ad,b,bd)) in cloning::gateway::layout(self.fragments.len()).unwrap().iter().enumerate() {
            self.fragments[i].doc=self.docs.iter().enumerate().find(|(j,d)| {
                if used.contains(j){return false;}
                let Ok(sites)=cloning::gateway::sites(&d.src) else{return false};
                let wanted_a=if ad==1{cloning::gateway::Kind::L}else{cloning::gateway::Kind::R};
                let wanted_b=if bd==2{cloning::gateway::Kind::L}else{cloning::gateway::Kind::R};
                sites.len()==2 && (sites[0].number,sites[0].direction,sites[0].kind)==(a,ad,wanted_a)
                    && (sites[1].number,sites[1].direction,sites[1].kind)==(b,bd,wanted_b)
            }).map(|(j,_)|j);
            if let Some(j)=self.fragments[i].doc{used.push(j);}
            self.fragments[i].use_pcr=false;
        }
    }

    fn pick_source<'a>(&'a self, pick: &'a Pick) -> Option<&'a Source> {
        if pick.use_std {
            return pick.std_src.as_ref();
        }
        pick.doc.and_then(|i| self.docs.get(i)).map(|d| &d.src)
    }

    fn source_name(&self, pick: &Pick) -> String {
        if let Some(n) = &pick.std_name {
            return n.clone();
        }
        if pick.use_std {
            return "< Please choose >".to_string();
        }
        pick.doc
            .and_then(|i| self.docs.get(i))
            .map(|d| d.name.clone())
            .unwrap_or_else(|| "< Please choose >".to_string())
    }

    /// 当前方法页签在左侧主视图里应显示的文档签名。
    fn view_doc_key(&self) -> String {
        match self.tab_kind(self.tab) {
            TabKind::Product => format!(
                "product|{}|{}|{}",
                self.preview_key,
                self.product_name,
                self.preview.as_ref().map(|p| p.seq.len()).unwrap_or(0)
            ),
            _ => {
                let pick = self.active_pick();
                let name = self.source_name(pick);
                let n = self.pick_source(pick).map(|s| s.seq.len()).unwrap_or(0);
                let std = pick.std_name.clone().unwrap_or_default();
                format!("source|{}|{}|{:?}|{}|{}|{}|{}", self.tab, self.frag_index,
                    pick.doc, pick.use_std, name, n, std)
            }
        }
    }

    /// 为左侧主视图构造一个完整的 `Document`。打开的文件优先用原始 `.dna`
    /// 字节（保留引物、块级信息）；商品化载体/合成来源则从 `Source` 生成。
    fn build_view_doc(&self) -> Option<crate::doc::Document> {
        match self.tab_kind(self.tab) {
            TabKind::Product => {
                let entry = if matches!(self.method, Method::GatewayBpLr(_)) && self.tab == 2 {
                    self.build_gateway_entry(self.frag_index).ok()
                } else { None };
                let p = if matches!(self.method, Method::GatewayBpLr(_)) && self.tab == 2 {
                    entry.as_ref()?
                } else { self.preview.as_ref()? };
                // 编辑会话统一处理原生末端 → 顶链显示坐标；预览不得提前丢弃末端。
                let name = if entry.is_some() { "Entry Clone" } else { &self.product_name };
                let bytes = output::product_bytes(self,p,name);
                let mut doc = crate::doc::Document::from_bytes(bytes).ok()?;
                doc.name_hint = Some(stem_label(name));
                Some(doc)
            }
            _ => {
                let pick = self.active_pick();
                let src = self.pick_source(pick)?;
                if let Some(i) = pick.doc.filter(|_| !pick.use_std) {
                    let label = self.docs.get(i).map(|d| d.name.clone());
                    if let Some(bytes) = self.docs.get(i).and_then(|d| d.bytes.clone()) {
                        if let Ok(mut doc) = crate::doc::Document::from_bytes(bytes) {
                            doc.name_hint = label.as_deref().map(stem_label);
                            return Some(doc);
                        }
                    }
                }
                let name = self.source_name(pick);
                let bytes = source_bytes(src, &name);
                let mut doc = crate::doc::Document::from_bytes(bytes).ok()?;
                doc.name_hint = Some(stem_label(&name));
                Some(doc)
            }
        }
    }

    /// 确保预览主视图存在，并在页签/来源/产物变化时换掉它的文档。
    fn ensure_preview(&mut self, ctx: &egui::Context) {
        let key = self.view_doc_key();
        let changed = key != self.preview_doc_key;
        if changed {
            self.preview_doc_key = key;
            self.preview_doc = self.build_view_doc();
            // 新页签没有可显示的来源时，丢掉旧的预览，让中央画布回到占位提示。
            if self.preview_doc.is_none() {
                self.preview_app = None;
            }
        }
        if self.preview_app.is_none() {
            if let Some(doc) = self.preview_doc.take() {
                self.preview_app = Some(Box::new(crate::HelixApp::preview(ctx, doc)));
            }
        } else if changed {
            if let Some(app) = self.preview_app.as_mut() {
                app.doc = self.preview_doc.take();
                app.needs_fit = true;
                let preset = app.enzyme_preset;
                app.apply_preset(preset);
                app.selection = None;
                app.cursor_bp = None;
                app.hover_bp = None;
            }
        }
        if let Some(app) = self.preview_app.as_mut() {
            app.tab = match self.preview_tab {
                PreviewTab::Map => crate::app::Tab::Map,
                PreviewTab::Sequence => crate::app::Tab::Sequence,
                PreviewTab::Enzymes => crate::app::Tab::Enzymes,
                PreviewTab::Primers => crate::app::Tab::Primers,
            };
        }
    }

    /// 两种来源都选好、能拼出产物时返回产物。
    pub fn build(&self) -> Result<Product, String> {
        match self.method {
            Method::Restriction | Method::RestrictionMulti => self.build_restriction(),
            Method::Ligation => self.build_ligation(),
            Method::GatewayBp => self.build_gateway_bp(),
            Method::GatewayLr(_) => self.build_gateway_lr(),
            Method::GatewayBpLr(_) => self.build_gateway_bp_lr(),
            m if m.is_ta_topo() => self.build_ta_topo(),
            Method::GoldenGate | Method::GoldenGateMulti => golden_gate::build(self),
            _ => self.build_overlap(),
        }
    }

    fn build_gateway_bp(&self) -> Result<Product, String> {
        self.build_gateway_entry(0)
    }

    fn build_gateway_entry(&self,index:usize) -> Result<Product,String> {
        let pick=self.fragments.get(index).ok_or("Please choose an attB Fragment.")?;
        let insert = self
            .pick_source(pick)
            .ok_or_else(|| "Please choose a source file.".to_string())?;
        let donor = self
            .pick_source(if matches!(self.method,Method::GatewayBpLr(n) if n>1) {
                self.gateway_donors.get(index).ok_or("Please choose a Donor Vector for every Fragment.")?
            }else{&self.vector})
            .ok_or_else(|| "Please specify a suitable Donor Vector.".to_string())?;
        let mut working = if pick.reversed { insert.flipped() } else { insert.clone() };
        if pick.use_pcr {
            if pick.primer1.trim().is_empty() || pick.primer2.trim().is_empty() {
                return Err("Please specify two primers.".into());
            }
            let region=pick.region.or_else(|| (!insert.circular).then(|| insert.top_range()))
                .ok_or("Please select the attB Fragment region to amplify.")?;
            let slice = cloning::homology::amplify(insert,region,pick.reversed,
                pick.primer1.trim().to_ascii_uppercase().as_bytes(),pick.primer2.trim().to_ascii_uppercase().as_bytes())?;
            working = Source::new(working.label, slice.seq, false).with_features(slice.features);
        }
        if matches!(self.method,Method::GatewayBpLr(n) if n>1) {cloning::gateway::bp(&working,donor)}
        else {cloning::gateway_recombine(&[&working], donor, "attB Fragment", true)}
    }

    fn build_gateway_lr(&self) -> Result<Product, String> {
        let mut inserts: Vec<&Source> = Vec::new();
        let oriented: Vec<_> = self.fragments.iter().map(|p| {
            let source = self.pick_source(p).ok_or_else(|| "Please choose an Entry Clone.".to_string())?;
            Ok(if p.reversed { source.flipped() } else { source.clone() })
        }).collect::<Result<_, String>>()?;
        for s in &oriented {
            inserts.push(s);
        }
        let dest = self
            .pick_source(&self.vector)
            .ok_or_else(|| "Please choose a Destination Vector.".to_string())?;
        cloning::gateway_recombine(&inserts, dest, "Entry Clone", false)
    }

    fn build_gateway_bp_lr(&self) -> Result<Product, String> {
        if self.fragments.len()>1 {
            let entries=(0..self.fragments.len()).map(|i|self.build_gateway_entry(i)
                .map(|p|Source::new(format!("Entry Clone {}",i+1),p.seq,p.circular).with_features(p.features)))
                .collect::<Result<Vec<_>,String>>()?;
            let dest=self.pick_source(&self.destination).ok_or("Please choose a Destination Vector.")?;
            return cloning::gateway::lr(&entries.iter().collect::<Vec<_>>(),dest);
        }
        let bp = self.build_gateway_bp()?;
        let entry = Source::new("Entry Clone", bp.seq, bp.circular).with_features(bp.features);
        let oriented = self.fragments.iter().skip(1).map(|p| {
            let s = self.pick_source(p).ok_or_else(|| "Please choose an Entry Clone.".to_string())?;
            Ok(if p.reversed { s.flipped() } else { s.clone() })
        }).collect::<Result<Vec<_>, String>>()?;
        let first = &entry;
        let mut all: Vec<&Source> = vec![first];
        all.extend(oriented.iter());
        let dest = self
            .pick_source(&self.destination)
            .ok_or_else(|| "Please choose a Destination Vector.".to_string())?;
        cloning::gateway_recombine(&all, dest, "Entry Clone", false)
    }

    fn build_restriction(&self) -> Result<Product, String> {
        let v = self
            .pick_source(&self.vector)
            .ok_or_else(|| "Please choose a source file.".to_string())?;
        let mut fragments = vec![restriction_fragment(v, &self.vector, true, true)?];
        let multi = self.method == Method::RestrictionMulti;
        let count = if multi { self.fragments.len() } else { 1 };
        for pick in self.fragments.iter().take(count) {
            let f = self
                .pick_source(pick)
                .ok_or_else(|| "Please choose a source file.".to_string())?;
            fragments.push(restriction_fragment(f, pick, false, true)?);
        }
        cloning::restriction_product(&fragments)
    }

    fn build_ligation(&self) -> Result<Product, String> {
        let mut fragments = Vec::new();
        for pick in &self.fragments {
            let s = self
                .pick_source(pick)
                .ok_or_else(|| "Please choose a source file.".to_string())?;
            fragments.push(restriction_fragment(s, pick, false, false)?);
        }
        cloning::ligate_reaction(&fragments, false)
    }

    fn build_ta_topo(&self) -> Result<Product, String> {
        // SnapGene 的 TA/GC/TOPO 需要一条**线性**载体（商品化 T 载体本身带
        // 3′ 突出端）；环状文件不能用。
        let v = self
            .pick_source(&self.vector)
            .ok_or_else(|| "Please choose a linear sequence.".to_string())?;
        if v.circular {
            return Err("Please choose a linear sequence.".into());
        }
        let f = self
            .pick_source(&self.fragments[0])
            .ok_or_else(|| "Please choose a source file.".to_string())?;
        let pick = &self.fragments[0];
        let insert = end_cloning_fragment(f, pick, self.method)?;
        cloning::end_cloning_product(v, &insert, self.method, pick.add_a)
    }

    fn build_overlap(&self) -> Result<Product, String> {
        let v = self
            .pick_source(&self.vector)
            .ok_or_else(|| "Please choose a source file.".to_string())?;
        let vs = homology_vector(v, &self.vector)?;
        let mut slices = vec![vs];
        let count = if self.method.is_multiple() {
            self.fragments.len()
        } else {
            1
        };
        for pick in self.fragments.iter().take(count) {
            let f = self
                .pick_source(pick)
                .ok_or_else(|| "Please choose a source file.".to_string())?;
            slices.push(homology_fragment(f, pick)?);
        }
        cloning::homology::assemble_circular(&slices, self.min_overlap())
    }

    fn vector_slice(&self, v: &Source) -> Slice {
        if self.overlap_method() {
            return homology_vector(v, &self.vector).unwrap_or_else(|_| v.slice(0, 0, false).sub(0, 0));
        }
        let n = v.seq.len();
        if self.vector.lin == Linearize::Restriction {
            if let Some(e) = self.vector.lin_enzyme.and_then(|i| hx_enzymes::ENZYMES.get(i)) {
                let cuts = cloning::cut_positions(&v.seq, v.circular, e);
                if let Some(&c) = cuts.first() {
                    return v.slice(c, c, self.vector.reversed);
                }
            }
        }
        let (a, b) = self.vector.region.unwrap_or((0, n));
        apply_primers(v.slice(a, b, self.vector.reversed), &self.vector)
    }

    fn fragment_slice(&self, f: &Source, pick: &Pick) -> Slice {
        if self.overlap_method() {
            return homology_fragment(f, pick).unwrap_or_else(|_| f.slice(0, 0, false).sub(0, 0));
        }
        let n = f.seq.len();
        let (a, b) = pick.region.unwrap_or((0, n));
        let mut s = f.slice(a, b, pick.reversed);
        if pick.use_pcr {
            s = apply_primers(s, pick);
        }
        s
    }

    fn refresh(&mut self) {
        let key = self.state_key();
        if key == self.preview_key && self.preview.is_some() {
            return;
        }
        self.preview_key = key;
        match self.build() {
            Ok(p) => {
                self.preview = Some(p);
                self.error = None;
                self.hint = self.success_hint();
            }
            Err(e) => {
                self.preview = None;
                self.error = Some(e.clone());
                self.hint = e;
            }
        }
    }

    fn success_hint(&self) -> String {
        if golden_gate::active(self) {
            return self.preview.as_ref().map(|p| format!("Ready to assemble\nProduct: {} bp", p.display_len()))
                .unwrap_or_default();
        }
        if !self.primers_ready() {
            return "Please specify two primers.".into();
        }
        // TA/GC/TOPO：能克隆时三页都显示绿色就绪提示；未就绪时提示框空白
        // （错误信息由 `refresh` 写进 `hint`）。
        if self.method.is_ta_topo() {
            if let Some(p) = &self.preview {
                return format!("Ready to clone\nProduct: {} bp", p.display_len());
            }
            return String::new();
        }
        if self.method.is_gateway() {
            if let Some(p) = &self.preview {
                return format!("Ready to clone\nProduct: {} bp", p.display_len());
            }
        }
        if self.method == Method::Ligation {
            if let Some(p) = &self.preview {
                return format!("Ready to ligate\nProduct: {} bp · linear", p.display_len());
            }
        }
        if self.method.is_restriction() {
            let v = enzyme_names(&self.vector);
            let f = enzyme_names(&self.fragments[0]);
            if let Some(p) = &self.preview {
                return format!("Ready to clone\nProduct: {} bp", p.display_len());
            }
            return format!("Vector: {v}; Fragment: {f}.");
        }
        if let Some(p) = &self.preview {
            if let Some(j) = p.junctions.first() {
                if j.overlap >= self.min_overlap().max(1) {
                    return format!("Overlap at the junction: {} bp.", j.overlap);
                }
            }
            return format!("Product: {} bp.", p.display_len());
        }
        "Please choose a source file.".into()
    }

    /// 提示框里实际显示的文案。TA/GC/TOPO 的 Product 页把「加 A」提示换成
    /// 通用的「末端不相容」，Linearized Vector 页则在没选载体时保持空白
    /// （照 SnapGene 8.2.3 实拍）。
    fn hint_text(&self) -> String {
        if matches!(self.method, Method::GatewayBpLr(_)) && self.tab == 2 {
            return match self.build_gateway_entry(self.frag_index) {
                Ok(p) => format!("Entry Clone: {} bp", p.display_len()),
                Err(e) => e,
            };
        }
        if self.method.is_ta_topo() {
            let kind = self.tab_kind(self.tab);
            if matches!(kind, TabKind::Product)
                && self.error.as_deref() == Some(cloning::TA_NEEDS_A)
            {
                return cloning::TA_INCOMPATIBLE.to_string();
            }
            if matches!(kind, TabKind::LinearizedVector)
                && self.error.as_deref() == Some("Please choose a linear sequence.")
            {
                return String::new();
            }
        }
        self.hint.clone()
    }

    fn overlap_method(&self) -> bool {
        matches!(
            self.method,
            Method::Gibson
                | Method::GibsonMulti
                | Method::GibsonAssemble
                | Method::InFusion
                | Method::InFusionMulti
                | Method::NeBuilder
                | Method::NeBuilderMulti
                | Method::NeBuilderAssemble
        )
    }

    /// 设计 PCR 引物时使用的同源臂长度（SnapGene `[actions]`）。
    ///
    /// Gibson / NEBuilder HiFi 读 `gibsonAssembly` / `hifiAssembly` 的
    /// `useExactLength` + `exactOverlapLength`；In-Fusion 单片段用
    /// `overlapLength`、多片段用 `multipleOverlapLength`。
    fn homology_overlap(&self) -> usize {
        match self.method {
            Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble => {
                self.settings.gibson.primer_overlap()
            }
            Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble => {
                self.settings.hifi.primer_overlap()
            }
            Method::InFusion => self.settings.in_fusion_overlap.max(1),
            Method::InFusionMulti => self.settings.in_fusion_multiple_overlap.max(1),
            _ => self.method.min_overlap(),
        }
    }

    /// 同源拼接要求的最小重叠长度。
    fn min_overlap(&self) -> usize {
        match self.method {
            Method::Gibson | Method::GibsonMulti | Method::GibsonAssemble => {
                self.settings.gibson.minimum_overlap()
            }
            Method::NeBuilder | Method::NeBuilderMulti | Method::NeBuilderAssemble => {
                self.settings.hifi.minimum_overlap()
            }
            Method::InFusion | Method::InFusionMulti => self.homology_overlap(),
            _ => self.method.min_overlap(),
        }
    }

    /// Golden Gate PCR 引物在 Type IIS 位点上游加的保护碱基数。
    fn gg_upstream_bases(&self) -> usize {
        self.settings.golden_gate_upstream_bases.clamp(0, 10)
    }

    fn fragment_uses_pcr(&self, kind: TabKind) -> bool {
        // 限制性克隆/连接全部用酶切，不显示 PCR 引物面板。
        if self.method.is_restriction() || self.method.is_ligation() {
            return false;
        }
        match kind {
            TabKind::Vector => self.overlap_method(),
            TabKind::Fragment => {
                self.overlap_method() || self.method.is_ta_topo() || self.method.is_bp()
            }
            TabKind::Fragments => !self.method.is_lr() && !self.method.is_ta_topo(),
            _ => false,
        }
    }

    /// 当前页签是否应显示引物编辑面板。
    fn tab_shows_primers(&self) -> bool {
        let kind = self.tab_kind(self.tab);
        if golden_gate::active(self) {
            return matches!(kind, TabKind::Vector | TabKind::Fragment | TabKind::Fragments)
                && self.active_pick().gg_mode == GoldenGateMode::PcrDigest;
        }
        if !self.fragment_uses_pcr(kind) {
            return false;
        }
        match kind {
            TabKind::Vector => self.vector.lin == Linearize::Pcr,
            TabKind::Fragment => self.fragments.get(self.frag_index).is_some_and(|p| p.use_pcr),
            TabKind::Fragments => self
                .fragments
                .get(self.frag_index)
                .is_some_and(|p| p.use_pcr),
            _ => false,
        }
    }

    fn primers_ready(&self) -> bool {
        if golden_gate::active(self) {
            return std::iter::once(&self.vector).chain(self.fragments.iter()).all(|p| {
                p.gg_mode != GoldenGateMode::PcrDigest || (!p.primer1.trim().is_empty() && !p.primer2.trim().is_empty())
            });
        }
        if self.overlap_method() && self.vector.lin == Linearize::Pcr {
            if self.vector.primer1.trim().is_empty()
                || self.vector.primer2.trim().is_empty()
            {
                return false;
            }
        }
        if self.method.is_ta_topo() || self.method.is_bp() {
            for p in &self.fragments {
                if p.use_pcr
                    && (p.primer1.trim().is_empty() || p.primer2.trim().is_empty())
                {
                    return false;
                }
            }
        }
        if self.overlap_method() && !self.method.is_lr() {
            let count = if self.method.is_multiple() {
                self.fragments.len()
            } else {
                1
            };
            for p in self.fragments.iter().take(count) {
                if p.use_pcr
                    && (p.primer1.trim().is_empty() || p.primer2.trim().is_empty())
                {
                    return false;
                }
            }
        }
        true
    }

    fn state_key(&self) -> String {
        let mut s = format!("{:?}|{}|{}|{}|", self.method, self.tab, self.gg_enzyme,self.frag_index);
        for p in std::iter::once(&self.vector)
            .chain(std::iter::once(&self.destination))
            .chain(self.fragments.iter()).chain(self.gateway_donors.iter())
        {
            s.push_str(&format!(
                "{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}/{:?}|",
                p.doc,
                p.use_std,
                p.std_name,
                p.reversed,
                p.enzyme,
                p.enzyme2,
                p.blunt,
                p.lin,
                p.lin_enzyme,
                p.region,
                p.use_pcr,
                p.add_a,
                p.primer1, p.primer2, p.phos1, p.phos2
            ));
            s.push_str(&format!("{:?}|", p.cut_pair));
            s.push_str(&format!("{:?}|", p.gg_mode));
        }
        s
    }

    // -----------------------------------------------------------------------
    // 商品化载体
    // -----------------------------------------------------------------------

    fn load_std_lists(&mut self) {
        let Some(base) = sg_resources_dir() else {
            return;
        };
        let plasmids = base.join("Plasmids");
        if self.method.is_bp() || self.method == Method::GatewayLr(1) {
            self.std_donor = list_std(&plasmids.join("Gateway Donor Vectors"), if self.fragments.len()>1 {&[]}else{DONOR_FILES});
        }
        if self.method.is_lr() {
            self.std_dest = list_std(&plasmids.join("Gateway Destination Vectors"), DEST_FILES);
        }
        let ta_dir = match self.method {
            Method::TaCloning => "TA Cloning Vectors",
            Method::GcCloning => "GC Cloning Vectors",
            Method::TopoTa => "TOPO TA Vectors",
            Method::TopoBlunt => "TOPO Blunt Vectors",
            Method::TopoDirectional => "TOPO Directional Vectors",
            _ => "",
        };
        if !ta_dir.is_empty() {
            self.std_ta = list_std(&plasmids.join(ta_dir), &[]);
        }
        // SG `topoTaCloning\vector` 记住上次使用的 TOPO TA 载体。只有在没有
        // 从打开的文档里自动选到载体时才套用，避免覆盖用户当前的来源。
        if self.method == Method::TopoTa
            && self.vector.doc.is_none()
            && !self.vector.use_std
        {
            let wanted = normalize_vector_name(&self.settings.topo_ta_vector);
            if let Some(entry) = self
                .std_ta
                .iter()
                .find(|e| normalize_vector_name(&e.name) == wanted)
            {
                let entry = entry.clone();
                select_std(&mut self.vector, &entry);
            }
        }
    }

}

/// 比较 SG 显示名（`pCR®2.1-TOPO®`）与资源文件名（`pCR2.1-TOPO (linearized).dna`）。
fn normalize_vector_name(name: &str) -> String {
    name.chars()
        .filter(|c| !matches!(c, '®' | '™' | '©'))
        .collect::<String>()
        .to_ascii_lowercase()
        .replace(" (linearized)", "")
        .replace(".dna", "")
        .trim()
        .to_string()
}

/// 读取标准载体的序列，挂到 pick 上。
fn select_std(pick: &mut Pick, entry: &StdEntry) -> bool {
    if let Some(src) = load_std_source(&entry.path) {
        pick.std_name = Some(entry.name.clone());
        pick.std_src = Some(src);
        pick.doc = None;
        pick.use_std = true;
        true
    } else {
        false
    }
}

// ---------------------------------------------------------------------------
// 引擎侧小工具
// ---------------------------------------------------------------------------

fn homology_vector(src: &Source, pick: &Pick) -> Result<Slice, String> {
    match pick.lin {
        Linearize::Direct => {
            if src.circular { return Err("Please choose a linearized Vector, or linearize it by PCR or restriction digestion.".into()); }
            let source = if pick.reversed { src.flipped() } else { src.clone() };
            let (a, b) = source.top_range();
            Ok(source.slice(a, b, false))
        }
        Linearize::Pcr => cloning::homology::amplify(src, pick.region.ok_or("Please select the Vector region to amplify.")?,
            pick.reversed, pick.primer1.trim().to_ascii_uppercase().as_bytes(), pick.primer2.trim().to_ascii_uppercase().as_bytes()),
        Linearize::Restriction => {
            let mut digest_pick = pick.clone();
            if let Some(e) = pick.lin_enzyme { digest_pick.enzyme = Some(e); }
            let f = restriction_fragment(src, &digest_pick, true, true)?;
            let (a, b) = f.top_range();
            Ok(f.slice(a, b, false))
        }
    }
}

fn end_cloning_fragment(source: &Source, pick: &Pick, method: Method) -> Result<Source, String> {
    let mut working = if pick.use_pcr {
        let region=pick.region.or_else(|| (!source.circular).then(|| source.top_range()))
            .ok_or("Please select the Fragment region to amplify.")?;
        let slice = cloning::homology::amplify(source, region,
            pick.reversed, pick.primer1.trim().to_ascii_uppercase().as_bytes(), pick.primer2.trim().to_ascii_uppercase().as_bytes())?;
        let pcr = Source::new(source.label.clone(), slice.seq, false).with_features(slice.features);
        if matches!(method, Method::TaCloning | Method::TopoTa) { cloning::add_a_overhangs(&pcr) }
        else if method == Method::GcCloning { cloning::add_three_prime_overhangs(&pcr, b'G') }
        else { pcr }
    } else if let Some((a, b)) = pick.region {
        if a >= source.seq.len() || b > source.seq.len() || (!source.circular && a >= b) {
            return Err("Please select a valid Fragment region.".into());
        }
        if a == 0 && b == source.seq.len() {
            if pick.reversed { source.flipped() } else { source.clone() }
        } else {
            let slice = source.slice(a, b, pick.reversed);
            Source::new(source.label.clone(), slice.seq, false).with_features(slice.features)
        }
    } else if pick.reversed { source.flipped() } else { source.clone() };
    if pick.add_a { working = cloning::add_three_prime_overhangs(&working, if method == Method::GcCloning { b'G' } else { b'A' }); }
    Ok(working)
}

fn homology_fragment(src: &Source, pick: &Pick) -> Result<Slice, String> {
    if pick.use_pcr {
        cloning::homology::amplify(src, pick.region.ok_or("Please select the Fragment region to amplify.")?,
            pick.reversed, pick.primer1.trim().to_ascii_uppercase().as_bytes(), pick.primer2.trim().to_ascii_uppercase().as_bytes())
    } else {
        if src.circular && pick.region.is_none() {
            return Err("Please select a linear Fragment region; a circular source cannot be inserted directly.".into());
        }
        let (a, b) = pick.region.unwrap_or(src.top_range());
        if a >= src.seq.len() || b > src.seq.len() || (!src.circular && a >= b) {
            return Err("Please select a valid Fragment region.".into());
        }
        Ok(src.slice(a, b, pick.reversed))
    }
}

fn apply_primers(mut s: Slice, pick: &Pick) -> Slice {
    let p1 = pick.primer1.trim().to_ascii_uppercase();
    let p2 = pick.primer2.trim().to_ascii_uppercase();
    if !p1.is_empty() {
        s.prepend(p1.as_bytes());
    }
    if !p2.is_empty() {
        let mut tail = p2.into_bytes();
        for b in tail.iter_mut() {
            *b = match *b {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                other => other,
            };
        }
        tail.reverse();
        s.append(&tail);
    }
    s
}

fn double_cuts_of(src: &Source, pick: &Pick) -> Vec<cloning::RestrictionCut> {
    let mut out = Vec::new();
    for e in [pick.enzyme, pick.enzyme2].into_iter().flatten() {
        let Some(en) = hx_enzymes::ENZYMES.get(e) else {
            continue;
        };
        out.extend(cloning::restriction_cuts(&src.seq, src.circular, en));
    }
    out.sort_by_key(|c| (c.top, c.bottom));
    out.dedup();
    out
}

fn digest_choices(src: &Source, pick: &Pick) -> Vec<cloning::RestrictionCut> {
    let mut cuts = double_cuts_of(src, pick);
    if !src.circular {
        let n = src.seq.len() as i64;
        cuts.push(cloning::RestrictionCut {
            top: src.ends.upstream_stickiness.min(0).unsigned_abs() as i64,
            bottom: src.ends.upstream_stickiness.max(0) as i64,
        });
        cuts.push(cloning::RestrictionCut {
            top: n - src.ends.downstream_stickiness.max(0) as i64,
            bottom: n - src.ends.downstream_stickiness.min(0).unsigned_abs() as i64,
        });
    }
    cuts.sort_by_key(|c| (c.top, c.bottom));
    cuts.dedup();
    cuts
}

fn selected_digest_pair(src: &Source, pick: &Pick, vector: bool)
    -> Option<(cloning::RestrictionCut, cloning::RestrictionCut)>
{
    let cuts = digest_choices(src, pick);
    if let Some((a, b)) = pick.cut_pair {
        if cuts.contains(&a) && cuts.contains(&b) { return Some((a, b)); }
    }
    // 两个及以上实际切口时默认选择它们之间的片段，而不是更短的线性末端碎屑。
    let actual = double_cuts_of(src, pick);
    let defaults = if actual.len() >= 2 { &actual } else { &cuts };
    let mut pairs: Vec<_> = defaults.windows(2).map(|w| (w[0], w[1])).collect();
    if src.circular && !cuts.is_empty() { pairs.push((*cuts.last().unwrap(), cuts[0])); }
    pairs.into_iter().filter_map(|(a, b)| {
        cloning::digest_fragment(src, a, b, pick.blunt, false).ok()
            .map(|f| ((a, b), f.display_len()))
    }).min_by_key(|(_, len)| if vector { usize::MAX - len } else { *len })
        .map(|(pair, _)| pair)
}

fn restriction_fragment(src: &Source, pick: &Pick, vector: bool, required: bool) -> Result<Source, String> {
    let kind = if vector { "Vector" } else { "Fragment" };
    let enzymes: Vec<_> = [pick.enzyme, pick.enzyme2].into_iter().flatten().collect();
    if enzymes.is_empty() {
        // SG 8.2.3 可直接使用已线性化来源的现有末端；只有环状来源必须选酶。
        if required && src.circular { return Err(format!("Select the enzyme(s) used to cut the {kind}.")); }
        if src.circular { return Err("A circular source must be linearized before ligation.".into()); }
        let mut fragment = if let Some((a, b)) = pick.region {
            if a >= b || b > src.seq.len() { return Err("Select a valid fragment region.".into()); }
            let slice = src.slice(a, b, false);
            let mut ends = src.ends.clone();
            if a != 0 { ends.upstream_stickiness = 0; ends.upstream_modification = "Unmodified".into(); }
            if b != src.seq.len() { ends.downstream_stickiness = 0; ends.downstream_modification = "Unmodified".into(); }
            ends.original_first_base = None;
            Source::new(src.label.clone(), slice.seq, false).with_features(slice.features).with_ends(ends)
        } else { src.clone() };
        if pick.blunt { fragment = cloning::blunt_ends(&fragment)?; }
        return Ok(if pick.reversed { fragment.flipped() } else { fragment });
    }
    for index in enzymes {
        let enzyme = hx_enzymes::ENZYMES.get(index).ok_or("Select a valid restriction enzyme.")?;
        if cloning::restriction_cuts(&src.seq, src.circular, enzyme).is_empty() {
            return Err(format!("The {kind} ({}) has no complete {} cut.", src.label, enzyme.name));
        }
    }
    let (a, b) = selected_digest_pair(src, pick, vector)
        .ok_or("Select two cuts that define a double-stranded fragment.")?;
    cloning::digest_fragment(src, a, b, pick.blunt, pick.reversed)
}

fn enzyme_names(pick: &Pick) -> String {
    let mut v = Vec::new();
    for e in [pick.enzyme, pick.enzyme2].into_iter().flatten() {
        if let Some(en) = hx_enzymes::ENZYMES.get(e) {
            v.push(en.name.to_string());
        }
    }
    v.join(" + ")
}

fn sg_resources_dir() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("HELIX_SNAP_SGENE_RESOURCES") {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
    }
    for p in [
        "/opt/gslbiotech/snapgene/resources",
        "/Applications/SnapGene.app/Contents/Resources",
    ] {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
    }
    None
}

/// Existing local templates only: no redistribution, download or fabricated DNA.
pub(crate) fn vector_library_entries(category: &str) -> Vec<(String, PathBuf)> {
    let Some(base) = sg_resources_dir() else { return Vec::new(); };
    list_std(&base.join("Plasmids").join(category), &[]).into_iter()
        .map(|e| (e.name, e.path)).collect()
}

fn list_std(dir: &Path, preferred: &[&str]) -> Vec<StdEntry> {
    let mut out = Vec::new();
    if preferred.is_empty() {
        if let Ok(rd) = std::fs::read_dir(dir) {
            let mut v: Vec<StdEntry> = rd
                .filter_map(|e| e.ok())
                .filter(|e| {
                    e.path()
                        .extension()
                        .map(|x| x.eq_ignore_ascii_case("dna"))
                        .unwrap_or(false)
                })
                .map(|e| {
                    let path = e.path();
                    // SnapGene 的商品化载体列表显示 `CustomMapLabel`（如
                    // `pCR®2.1-TOPO®`），文件名只在没有标签时兜底。
                    let label = std::fs::read(&path)
                        .ok()
                        .and_then(|b| hx_formats::DnaFile::parse(&b).ok())
                        .and_then(|f| f.custom_map_label());
                    let name = label.unwrap_or_else(|| {
                        path.file_stem()
                            .map(|s| s.to_string_lossy().to_string())
                            .unwrap_or_default()
                    });
                    StdEntry { name, path }
                })
                .collect();
            v.sort_by(|a, b| a.name.cmp(&b.name));
            out = v;
        }
        return out;
    }
    for name in preferred {
        let path = dir.join(format!("{name}.dna"));
        if path.is_file() {
            out.push(StdEntry {
                name: (*name).to_string(),
                path,
            });
        }
    }
    out
}

fn load_std_source(path: &Path) -> Option<Source> {
    let bytes = std::fs::read(path).ok()?;
    let doc = crate::doc::Document::from_bytes(bytes).ok()?;
    Some(doc.cloning_source(
            path.file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "Vector".into()),
    ))
}

/// 打开文件并作为来源追加到对话框（拖放用）。
pub fn doc_from_path(path: &Path) -> Option<DocSnap> {
    let bytes = std::fs::read(path).ok()?;
    let doc = crate::doc::Document::from_bytes(bytes.clone()).ok()?;
    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "Sequence".into());
    Some(
        DocSnap::new(
            name,
            doc.cloning_source(
                path.file_stem()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| "Sequence".into()),
            ),
        )
        .with_bytes(bytes),
    )
}

// ---------------------------------------------------------------------------
// 窗口
// ---------------------------------------------------------------------------

const SIDE_W: f32 = 330.0;

pub enum Outcome {
    /// 生成产物（`.dna` 字节 + 建议文件名 + 是否关窗）。
    Clone(Vec<u8>, String, bool),
    Cancel,
}

pub fn draw(dlg: &mut CloningDlg, ctx: &egui::Context) -> Option<Outcome> {
    dlg.refresh();
    if dlg.frames_open < 2 {
        dlg.frames_open = dlg.frames_open.saturating_add(1);
        ctx.request_repaint();
    }

    // SnapGene 的克隆对话框是一个**独立的原生窗口**：有自己的标题栏、
    // 可以移动/缩放。这里用 immediate viewport 把它挂到主窗口的同一帧里。
    let viewport_id = egui::ViewportId::from_hash_of("cloning_dialog");
    let builder = egui::ViewportBuilder::default()
        .with_title({let title=dlg.method.title();let action=dlg.method.action();if dlg.title()==title{crate::i18n::text(ctx,&title).to_owned()}else{format!("{}: {}",crate::i18n::text(ctx,&title),crate::i18n::text(ctx,action))}})
        .with_inner_size([1220.0, 840.0])
        .with_min_inner_size([920.0, 620.0])
        .with_resizable(true);
    let mut out: Option<Outcome> = None;
    ctx.show_viewport_immediate(viewport_id, builder, |ui, _class| {
        let vctx = ui.ctx().clone();
        if vctx.input(|i| i.viewport().close_requested()) {
            out = Some(Outcome::Cancel);
            vctx.send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }
        out = draw_dialog(dlg, &vctx, ui);
        if out.is_none() && vctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            out = Some(Outcome::Cancel);
        }
        vctx.request_repaint();
    });
    out
}

/// 原生窗口内部布局：顶部方法页签 / 右侧参数栏 / 中间主视图 / 底部引物或示意图。
fn draw_dialog(
    dlg: &mut CloningDlg,
    ctx: &egui::Context,
    ui: &mut egui::Ui,
) -> Option<Outcome> {
    if golden_gate::active(dlg) { return golden_gate::draw_dialog(ui, ctx, dlg); }
    let mut out: Option<Outcome> = None;
    let tabs = dlg.tabs();
    let gg = false; // Golden Gate has already returned through its dedicated pane layout.
    let side_w = SIDE_W;

    // ---- 顶部方法页签（attB Fragment / Donor Vector / Entry Clone）----
    egui::Panel::top("cloning_method_tabs")
        .exact_size(36.0)
        .frame(
            egui::Frame::new()
                .fill(Color32::from_rgb(0xf2, 0xf2, 0xf2))
                .inner_margin(egui::Margin::symmetric(8, 3)),
        )
        .show_inside(ui, |ui| {
            let widths: Vec<f32> = tabs
                .iter()
                .map(|t| (t.chars().count() as f32 * 7.4 + 30.0).max(86.0))
                .collect();
            let total: f32 = widths.iter().sum::<f32>() + 2.0 * tabs.len() as f32;
            ui.horizontal(|ui| {
                ui.add_space(((ui.available_width() - total) * 0.5).max(0.0));
                for (i, label) in tabs.iter().enumerate() {
                    if method_tab(ui, dlg.tab == i, label, widths[i]) {
                        dlg.tab = i;
                        if i + 1 == tabs.len() {
                            dlg.frag_index = 0;
                        }
                    }
                }
            });
        });

    // 产物就绪时 SnapGene 的提示框是浅绿色，否则浅黄色（TA/GC/TOPO 三页
    // 都一样）。
    let ready = dlg.preview.is_some() && dlg.primers_ready();
    let green = ready;
    let (hint_fill, hint_stroke) = if green {
        (
            Color32::from_rgb(0xd1, 0xfa, 0xd1),
            Color32::from_rgb(0x8e, 0x8e, 0x8e),
        )
    } else {
        (
            Color32::from_rgb(0xff, 0xff, 0x85),
            Color32::from_rgb(0x8e, 0x8e, 0x8e),
        )
    };

    // ---- 右侧参数栏（SG 的右侧是固定宽度的白色面板）----
    egui::Panel::right("cloning_side")
        .exact_size(side_w)
        .frame(
            egui::Frame::new()
                .fill(Color32::WHITE)
                .inner_margin(egui::Margin::symmetric(9, 7)),
        )
        .show_inside(ui, |ui| {
            // 底部：黄色提示 + Create product + 产物名 + Clone/Cancel
            egui::Panel::bottom("cloning_actions")
                .exact_size(if gg { 278.0 } else { 234.0 })
                .frame(
                    egui::Frame::new()
                        .fill(Color32::WHITE)
                        .inner_margin(egui::Margin {
                            left: 0,
                            right: 0,
                            top: 0,
                            bottom: 10,
                        }),
                )
                .show_inside(ui, |ui| {
                    egui::Frame::NONE
                        .fill(hint_fill)
                        .stroke(Stroke::new(1.0, hint_stroke))
                        .inner_margin(egui::Margin::symmetric(8, 6))
                        .show(ui, |ui| {
                            ui.set_width(side_w - 34.0);
                            ui.set_height(if gg { 170.0 } else { 104.0 });
                            let text = crate::i18n::message(ui.ctx(),&dlg.hint_text());
                            let needs_a = (dlg.error.as_deref() == Some(cloning::TA_NEEDS_A)
                                || (dlg.method == Method::GcCloning && dlg.error.as_deref() == Some("GC cloning requires 3' G overhangs in the Fragment.")))
                                && matches!(dlg.tab_kind(dlg.tab), TabKind::Fragment);
                            if needs_a {
                                ui.add_space(12.0);
                                ui.label(
                                    RichText::new(&text)
                                        .size(11.5)
                                        .color(Color32::from_rgb(0x33, 0x33, 0x22)),
                                );
                                ui.add_space(7.0);
                                if ui
                                    .button(
                                        RichText::new(crate::i18n::text(ui.ctx(),if dlg.method == Method::GcCloning { "Add 3' G overhangs to Insert" } else { "Add 3' A overhangs to Insert" }))
                                            .size(11.5),
                                    )
                                    .clicked()
                                {
                                    if let Some(p) = dlg.fragments.first_mut() {
                                        p.add_a = true;
                                    }
                                }
                            } else if green {
                                ui.centered_and_justified(|ui| {
                                    ui.label(
                                        RichText::new(&text)
                                            .size(11.5)
                                            .color(Color32::from_rgb(0x33, 0x33, 0x22)),
                                    );
                                });
                            } else if gg {
                                ui.centered_and_justified(|ui| { ui.label(RichText::new(&text).size(11.5)); });
                            } else if !text.is_empty() {
                                // 未就绪的提示靠上显示（SnapGene 的黄色框里文字
                                // 贴着上沿，下面留白）。
                                ui.add_space(14.0);
                                ui.label(
                                    RichText::new(&text)
                                        .size(11.5)
                                        .color(Color32::from_rgb(0x33, 0x33, 0x22)),
                                );
                            }
                        });
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Create product:")).size(11.0));
                        // SnapGene 允许取消勾选后只做模拟、不生成产物文档。HX 目前
                        // 总是创建产物，因此这里保持勾选但明确禁用，避免假控件。
                        if !gg {
                            ui.add_enabled(
                                false,
                                egui::Checkbox::new(&mut dlg.create_product, ""),
                            )
                            .on_disabled_hover_text(
                                crate::i18n::text(ui.ctx(),"Always creates the product; skipping product creation is not implemented yet."),
                            );
                        }
                        ui.with_layout(
                            egui::Layout::right_to_left(egui::Align::Center),
                            |ui| {
                                ui.checkbox(
                                    &mut dlg.close_after_clone,
                                    RichText::new(crate::i18n::text(ui.ctx(),"and close this window")).size(11.0),
                                );
                            },
                        );
                    });
                    ui.add_space(3.0);
                    ui.add(
                        egui::TextEdit::singleline(&mut dlg.product_name)
                            .desired_width(ui.available_width()),
                    );
                    ui.add_space(6.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let label = match dlg.method {
                            Method::Gibson
                            | Method::GibsonMulti
                            | Method::GibsonAssemble
                            | Method::InFusion
                            | Method::InFusionMulti
                            | Method::NeBuilder
                            | Method::NeBuilderMulti
                            | Method::NeBuilderAssemble
                            | Method::GoldenGate
                            | Method::GoldenGateMulti => "Assemble",
                            Method::Ligation => "Ligate",
                            _ => "Clone",
                        };
                        // 右对齐布局里先加的控件在最右边，所以先放 Cancel、
                        // 再放 Clone，视觉上才是 SnapGene 的 `Clone | Cancel`。
                        if ui
                            .add_sized(egui::vec2(78.0, 26.0), egui::Button::new(crate::i18n::text(ui.ctx(),"Cancel")))
                            .clicked()
                        {
                            out = Some(Outcome::Cancel);
                        }
                        let clone = ui.add_enabled_ui(ready, |ui| {
                            ui.add_sized(egui::vec2(78.0, 26.0), egui::Button::new(crate::i18n::text(ui.ctx(),label)))
                        });
                        if clone.inner.clicked() {
                            if let Some(p) = &dlg.preview {
                                let bytes = output::product_bytes(dlg,p, &dlg.product_name);
                                out = Some(Outcome::Clone(
                                    bytes,
                                    dlg.product_name.clone(),
                                    dlg.close_after_clone,
                                ));
                            }
                        }
                    });
                });

            // 上半部分可滚动（参数很多，窗口缩小后不能把按钮挤走）。
            egui::CentralPanel::default()
                .frame(egui::Frame::new().fill(Color32::WHITE))
                .show_inside(ui, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("cloning_side_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            ui.set_width(ui.available_width());
                            side_panel(ui, dlg);
                        });
                });
        });

    // Restriction's SG Vector/Fragment pages retain the reaction scheme too.
    let kind = dlg.tab_kind(dlg.tab);
    let primer_h: f32 = if show_primer_panel(dlg)
        && matches!(
            kind,
            TabKind::Vector | TabKind::Fragment | TabKind::Fragments
        )
    {
        if gg { 230.0 } else { 148.0 }
    } else {
        0.0
    };
    let schematic_h: f32 = if matches!(kind, TabKind::Product) || dlg.method.is_ta_topo() || dlg.method.is_restriction() {
        if let Method::GatewayBpLr(n)=dlg.method { if n>1 {((n+1) as f32*92.0+24.0).min(490.0)}else{236.0} }
        else if gg {236.0}else{150.0}
    } else {
        0.0
    };
    // SnapGene 的 attB 页在“直接使用”时，底部画 attB1/attB2 位点细节，
    // 而不是引物面板或产物示意图。
    let attb_h: f32 = if shows_attb_detail(dlg) { 150.0 } else { 0.0 };
    let bottom_h = primer_h.max(schematic_h).max(attb_h);
    if bottom_h > 0.0 {
        egui::Panel::bottom("cloning_bottom")
            .exact_size(bottom_h)
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(0xfa, 0xfa, 0xfa))
                    .inner_margin(egui::Margin::symmetric(8, 5)),
            )
            .show_inside(ui, |ui| {
                if primer_h > 0.0 {
                    if gg { golden_gate::primer_panel(ui, dlg); } else { primer_panel(ui, dlg); }
                } else if attb_h > 0.0 {
                    let rect = ui.max_rect();
                    paint_attb_detail(ui, rect, dlg);
                } else {
                    let rect = ui.max_rect();
                    if gg { golden_gate::draw_join(ui.painter(), rect, false, ready); }
                    else { paint_schematic(ui, rect, dlg); }
                }
            });
    }

    // ---- 中央：完整复用主窗口的 DNA 视图（SnapGene 把它放在一张白色卡片里）----
    egui::CentralPanel::default()
        .frame(
            egui::Frame::new()
                .fill(Color32::from_rgb(0xef, 0xef, 0xef))
                .inner_margin(egui::Margin::symmetric(8, 6)),
        )
        .show_inside(ui, |ui| {
            egui::Frame::new()
                .fill(theme::SEQ_BG)
                .stroke(Stroke::new(0.8, theme::BORDER))
                .show(ui, |ui| {
                    ui.set_min_size(ui.available_size());
                    dlg.ensure_preview(ctx);
                    let overlap = dlg.overlap_method();
                    if let Some(app) = dlg.preview_app.as_mut() {
                        app.doc_canvas(ui, crate::app::CanvasMode::Cloning);
                        if gg && dlg.tab < 2 {
                            if let Some((a, b)) = app.selection {
                                let pick = if dlg.tab == 0 { &mut dlg.vector } else { &mut dlg.fragments[dlg.frag_index] };
                                let range = if dlg.tab == 0 { (b, a) } else { (a, b) };
                                if pick.region != Some(range) { pick.region = Some(range); pick.cut_pair = None; }
                            }
                        }
                        dlg.preview_tab = match app.tab {
                            crate::app::Tab::Sequence => PreviewTab::Sequence,
                            crate::app::Tab::Enzymes => PreviewTab::Enzymes,
                            crate::app::Tab::Primers => PreviewTab::Primers,
                            _ => PreviewTab::Map,
                        };
                        if overlap && dlg.tab < 2 {
                            let selected = app.selection.or_else(|| app.cursor_bp.map(|at| (at, at)));
                            if let Some((a, b)) = selected {
                                let circular = app.doc.as_ref().is_some_and(|d| d.circular());
                                let pick = if dlg.tab == 0 { &mut dlg.vector } else { &mut dlg.fragments[dlg.frag_index] };
                                pick.region = Some(if dlg.tab == 0 && circular { (b, a) } else { (a, b) });
                            }
                        }
                    } else {
                        ui.centered_and_justified(|ui| {
                            paint_empty_canvas_hint(ui, dlg);
                        });
                    }
                });
        });

    out
}

/// 中央画布没有来源文件时的占位文字（照 SnapGene 的 Linearized Vector 页）。
fn paint_empty_canvas_hint(ui: &mut egui::Ui, dlg: &CloningDlg) {
    if golden_gate::active(dlg) && matches!(dlg.tab_kind(dlg.tab), TabKind::Fragment | TabKind::Fragments) {
        ui.vertical_centered(|ui| {
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"To choose the source file for the Fragment, either:")).strong());
            ui.label(crate::i18n::text(ui.ctx(),"• Use the menu at the right."));
            ui.label(crate::i18n::text(ui.ctx(),"• Drag the source file into this box."));
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"• Click here.")).color(Color32::from_rgb(0x2a, 0x5f, 0xb8)));
        });
        return;
    }
    if matches!(dlg.tab_kind(dlg.tab), TabKind::LinearizedVector) {
        ui.vertical_centered(|ui| {
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"To choose the Linearized Vector, either:")).size(12.5).strong());
            ui.add_space(8.0);
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"• Use the menu at the right.")).size(12.0));
            ui.add_space(3.0);
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"• Drag the source file into this box.")).size(12.0));
            ui.add_space(3.0);
            ui.label(
                RichText::new(crate::i18n::text(ui.ctx(),"• Click here."))
                    .size(12.0)
                    .color(Color32::from_rgb(0x2a, 0x5f, 0xb8)),
            );
        });
        return;
    }
    ui.label(
        RichText::new(crate::i18n::text(ui.ctx(),"Please choose a source file."))
            .size(13.0)
            .color(theme::TEXT_DIM),
    );
}

/// BP 的 attB 片段页在“直接使用”时是否要画 attB1/attB2 位点细节。
fn shows_attb_detail(dlg: &CloningDlg) -> bool {
    matches!(dlg.tab_kind(dlg.tab), TabKind::Fragment)
        && dlg.method.is_bp()
        && dlg.fragments.get(dlg.frag_index).is_some_and(|p| !p.use_pcr)
}

/// attB 位点细节（SnapGene 的 attB Fragment 页底部）：
///
/// ```text
///        ...CAAGTTT GTACAAAAAAGCAG ...   (To Be Inserted)   ...CAGCTTT CTTGTACAAAGTGG ...
///        ...GTTCAAA CATGTTTTTTCGTC ...                       ...GTCGAAA GAACATGTTTCACC ...
///         attB1 ->                                              <- attB2
/// ```
///
/// 两条链上下对齐，重组核心区用粉色带标出，中间写 `(To Be Inserted)`。
fn paint_attb_detail(ui: &mut egui::Ui, rect: Rect, dlg: &CloningDlg) {
    use hx_algo::cloning::{ATT_B1, ATT_B2};
    let Some(src) = current_source(dlg) else {
        return;
    };
    let source=if dlg.fragments.get(dlg.frag_index).is_some_and(|p|p.reversed) {src.flipped()}else{src.clone()};
    let seq=&source.seq;
    let find = |pat: &[u8]| seq.windows(pat.len()).position(|w| w == pat);
    let (b1_flank,b1_core,b2_flank,b2_core,left_label,right_label)=if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {
        let Ok(sites)=cloning::gateway::sites(&source) else{return};
        let sites=sites.iter().filter(|s|s.kind==cloning::gateway::Kind::B).collect::<Vec<_>>();
        if sites.len()!=2{return;}
        let (a,b)=(sites[0],sites[1]);
        let label=|s:&cloning::gateway::Site|if s.direction==1 {format!("attB{} →",s.number)}else{format!("← attB{}",s.number)};
        (&seq[a.start..a.cut-7],&seq[a.cut-7..a.cut+7],&seq[b.start..b.cut-7],&seq[b.cut-7..b.cut+7],label(a),label(b))
    } else {
        if find(&ATT_B1[8..22]).is_none() || find(&ATT_B2[10..24]).is_none() {return;}
        (&ATT_B1[1..8],&ATT_B1[8..22],&ATT_B2[3..10],&ATT_B2[10..24],"attB1 →".into(),"← attB2".into())
    };
    let comp = |s: &[u8]| -> String {
        s.iter()
            .map(|b| match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                c => *c as char,
            })
            .collect()
    };
    let txt = |s: &[u8]| String::from_utf8_lossy(s).to_string();

    let plain = Color32::from_rgb(0x8a, 0x8a, 0x8a);
    let ink = theme::TEXT;
    let pink = Color32::from_rgb(0xfd, 0x96, 0x96);
    let font = egui::FontId::monospace(11.0);
    let painter = ui.painter();

    // 每条链的两行 + 下方位点名。
    let top_y = rect.top() + 38.0;
    let bot_y = top_y + 17.0;

    let segs = |flank: &[u8], core: &[u8]| -> Vec<(String, Color32)> {
        vec![
            ("...".to_string(), plain),
            (txt(flank), plain),
            (txt(core), ink),
            ("...".to_string(), plain),
        ]
    };
    let left_segs = segs(b1_flank, b1_core);
    let right_segs = segs(b2_flank, b2_core);
    let widths = |segs: &[(String, Color32)]| -> Vec<f32> {
        segs.iter()
            .map(|(t, c)| {
                painter
                    .layout_no_wrap(t.clone(), font.clone(), *c)
                    .size()
                    .x
            })
            .collect()
    };
    let lw = widths(&left_segs);
    let rw = widths(&right_segs);
    let left_x = rect.left() + 26.0;
    let right_end = rect.right() - 26.0;
    let right_x = right_end - rw.iter().sum::<f32>();
    let core1_x0 = left_x + lw[0] + lw[1];
    let core2_x0 = right_x + rw[0] + rw[1];
    let core2_x1 = core2_x0 + rw[2];

    // 粉色带：从左核心区起点到右核心区终点（= 将被换入的片段）。
    let band = Rect::from_min_max(
        Pos2::new(core1_x0 - 1.0, top_y - 3.0),
        Pos2::new(core2_x1 + 1.0, bot_y + 15.0),
    );
    painter.rect_filled(band, 0.0, pink);

    // 两条链的序列。
    let draw_segs = |y: f32, x0: f32, segs: &[(String, Color32)], ws: &[f32]| {
        let mut x = x0;
        for ((t, c), w) in segs.iter().zip(ws) {
            painter.text(
                Pos2::new(x, y),
                egui::Align2::LEFT_TOP,
                t,
                font.clone(),
                *c,
            );
            x += w;
        }
    };
    draw_segs(top_y, left_x, &left_segs, &lw);
    draw_segs(top_y, right_x, &right_segs, &rw);

    // 下链用互补碱基（不反向），与 SnapGene 一样逐位对齐。
    let left_bot = format!(
        "{}...",
        comp(&[b1_flank, b1_core].concat())
    );
    let right_bot = format!(
        "{}...",
        comp(&[b2_flank, b2_core].concat())
    );
    painter.text(
        Pos2::new(left_x, bot_y),
        egui::Align2::LEFT_TOP,
        format!("...{left_bot}"),
        font.clone(),
        plain,
    );
    painter.text(
        Pos2::new(right_x, bot_y),
        egui::Align2::LEFT_TOP,
        format!("...{right_bot}"),
        font.clone(),
        plain,
    );

    // 中间提示。
    painter.text(
        Pos2::new(band.center().x, top_y + 8.0),
        egui::Align2::CENTER_CENTER,
        crate::i18n::text(ui.ctx(),"(To Be Inserted)"),
        egui::FontId::proportional(11.0),
        ink,
    );

    // 位点名。
    let label_y = bot_y + 24.0;
    painter.text(
        Pos2::new(core1_x0 + 16.0, label_y),
        egui::Align2::LEFT_TOP,
        left_label,
        egui::FontId::proportional(11.5),
        ink,
    );
    painter.text(
        Pos2::new(core2_x1, label_y),
        egui::Align2::RIGHT_TOP,
        right_label,
        egui::FontId::proportional(11.5),
        ink,
    );
}

fn show_primer_panel(dlg: &CloningDlg) -> bool {
    dlg.tab_shows_primers()
}

fn method_tab(ui: &mut egui::Ui, active: bool, label: &str, width: f32) -> bool {
    let label=crate::i18n::text(ui.ctx(),label);
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(width, 30.0), egui::Sense::click());
    let fill = if active {
        Color32::WHITE
    } else if resp.hovered() {
        Color32::from_rgb(0xf8, 0xf8, 0xf8)
    } else {
        Color32::from_rgb(0xf4, 0xf4, 0xf4)
    };
    let p = ui.painter();
    p.rect(
        rect,
        0.0,
        fill,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    if active {
        p.line_segment(
            [rect.left_bottom(), rect.right_bottom()],
            Stroke::new(2.0, theme::ACCENT),
        );
    }
    p.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        label,
        egui::FontId::proportional(12.0),
        if active { theme::TEXT } else { theme::TEXT_DIM },
    );
    resp.clicked()
}

fn tab_button(ui: &mut egui::Ui, active: bool, label: &str, width: f32) -> bool {
    let label=crate::i18n::text(ui.ctx(),label);
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(width, 25.0), egui::Sense::click());
    let fill = if active {
        Color32::WHITE
    } else {
        Color32::from_rgb(0xf4, 0xf4, 0xf4)
    };
    ui.painter().rect(
        rect,
        0.0,
        fill,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    ui.painter().text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        label,
        egui::FontId::proportional(12.0),
        theme::TEXT,
    );
    resp.clicked()
}

fn small_tab(ui: &mut egui::Ui, active: bool, label: &str) -> bool {
    let label=crate::i18n::text(ui.ctx(),label);
    let w = (label.chars().count() as f32 * 7.0 + 20.0).max(54.0);
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(w, 24.0), egui::Sense::click());
    let fill = if active {
        Color32::WHITE
    } else {
        Color32::from_rgb(0xf2, 0xf2, 0xf2)
    };
    ui.painter().rect(
        rect,
        0.0,
        fill,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    ui.painter().text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        label,
        egui::FontId::proportional(11.5),
        theme::TEXT,
    );
    resp.clicked()
}

fn number_badge(ui: &mut egui::Ui, active: bool, n: usize) -> bool {
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(21.0, 20.0), egui::Sense::click());
    let fill = if active {
        theme::ACCENT
    } else {
        Color32::WHITE
    };
    ui.painter().rect(
        rect,
        1.0,
        fill,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    ui.painter().text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        n.to_string(),
        egui::FontId::proportional(11.0),
        if active { Color32::WHITE } else { theme::TEXT },
    );
    resp.clicked()
}

// ---------------------------------------------------------------------------
// 左栏：预览
// ---------------------------------------------------------------------------

fn preview_panel(ui: &mut egui::Ui, size: Vec2, dlg: &mut CloningDlg) {
    let (rect, _) = ui.allocate_exact_size(size, egui::Sense::hover());
    ui.painter().rect_stroke(
        rect,
        0.0,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect.shrink(6.0)));
    child.horizontal(|ui| {
        for (tab, label) in [
            (PreviewTab::Map, "Map"),
            (PreviewTab::Sequence, "Sequence"),
            (PreviewTab::Enzymes, "Enzymes"),
            (PreviewTab::Primers, "Primers"),
        ] {
            if small_tab(ui, dlg.preview_tab == tab, label) {
                dlg.preview_tab = tab;
            }
        }
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Unique 6+ Cutters | Nonredundant")).size(9.5).color(theme::TEXT_DIM));
            ui.add_space(4.0);
            egui::ComboBox::from_id_salt("cloning_display")
                .selected_text(crate::i18n::text(ui.ctx(),"Display"))
                .width(74.0)
                .show_ui(ui, |ui| {
                    ui.label(crate::i18n::text(ui.ctx(),"Map / Sequence"));
                });
            ui.label(RichText::new("⇄").size(12.0).color(theme::TEXT_DIM));
            ui.label(RichText::new("◉").size(14.0).color(theme::ACCENT));
        });
    });
    let content = Rect::from_min_max(
        Pos2::new(rect.left() + 6.0, rect.top() + 34.0),
        Pos2::new(rect.right() - 6.0, rect.bottom() - 6.0),
    );
    match dlg.preview_tab {
        PreviewTab::Map => paint_preview(ui, content, dlg),
        PreviewTab::Sequence => paint_sequence(ui, content, dlg),
        PreviewTab::Enzymes => paint_enzymes(ui, content, dlg),
        PreviewTab::Primers => paint_primers(ui, content, dlg),
    }
}

fn drop_zone_panel(ui: &mut egui::Ui, size: Vec2, dlg: &CloningDlg) {
    let (rect, _) = ui.allocate_exact_size(size, egui::Sense::hover());
    let inner = rect.shrink2(egui::vec2(44.0, 30.0));
    let painter = ui.painter_at(rect);
    painter.rect_filled(inner, 0.0, Color32::from_rgb(0xe9, 0xe9, 0xe9));
    painter.rect_stroke(
        inner,
        0.0,
        Stroke::new(0.8, Color32::from_rgb(0xc8, 0xc8, 0xc8)),
        egui::StrokeKind::Inside,
    );
    let title = if dlg.tab_kind(dlg.tab) == TabKind::DestinationVector {
        "To choose the Destination Vector, either:"
    } else {
        "To choose the Donor Vector, either:"
    };
    let cy = inner.center().y - 34.0;
    painter.text(
        Pos2::new(inner.center().x, cy),
        egui::Align2::CENTER_CENTER,
        crate::i18n::text(ui.ctx(),title),
        egui::FontId::proportional(12.0),
        theme::TEXT,
    );
    for (i, line) in [
        "• Use the menu at the right.",
        "• Drag the file into this box.",
        "• Click here.",
    ]
    .iter()
    .enumerate()
    {
        let color = if i == 2 {
            theme::ACCENT
        } else {
            theme::TEXT_DIM
        };
        painter.text(
            Pos2::new(inner.center().x - 92.0, cy + 24.0 + i as f32 * 18.0),
            egui::Align2::LEFT_CENTER,
            crate::i18n::text(ui.ctx(),line),
            egui::FontId::proportional(11.5),
            color,
        );
    }
}

fn current_source<'a>(dlg: &'a CloningDlg) -> Option<&'a Source> {
    match dlg.tab_kind(dlg.tab) {
        TabKind::Product => None,
        TabKind::Fragment | TabKind::Fragments => dlg
            .fragments
            .get(dlg.frag_index)
            .and_then(|p| dlg.pick_source(p)),
        TabKind::DestinationVector if matches!(dlg.method, Method::GatewayBpLr(_)) => {
            dlg.pick_source(&dlg.destination)
        }
        _ => dlg.pick_source(&dlg.vector),
    }
}

/// `Source` → `.dna` 字节（无原始文件时的回退；保留序列/特征/拓扑）。
fn source_bytes(src: &Source, name: &str) -> Vec<u8> {
    let prod = Product {
        seq: src.seq.clone(),
        features: src.features.clone(),
        junctions: Vec::new(),
        spans: Vec::new(),
        circular: src.circular,
        ends: src.ends.clone(),
    };
    cloning::product_bytes(&prod, name)
}

/// 预览文档的显示名：与 SnapGene 一样用去掉扩展名的文件名。
fn stem_label(name: &str) -> String {
    Path::new(name)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(name)
        .to_string()
}

fn paint_sequence(ui: &mut egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let seq = match dlg.tab_kind(dlg.tab) {
        TabKind::Product => dlg.preview.as_ref().map(|p| p.seq.clone()).unwrap_or_default(),
        _ => current_source(dlg).map(|s| s.seq.clone()).unwrap_or_default(),
    };
    if seq.is_empty() {
        ui.painter().text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            crate::i18n::text(ui.ctx(),"Please choose a source file."),
            egui::FontId::proportional(12.0),
            theme::TEXT_DIM,
        );
        return;
    }
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect));
    egui::ScrollArea::vertical()
        .id_salt("cloning_seq")
        .show(&mut child, |ui| {
            let text: String = seq
                .chunks(60)
                .map(|c| String::from_utf8_lossy(c).to_string())
                .collect::<Vec<_>>()
                .join("\n");
            ui.label(RichText::new(text).monospace().size(10.5));
        });
}

fn paint_enzymes(ui: &mut egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let Some(src) = current_source(dlg) else {
        ui.painter().text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "Please choose a source file.",
            egui::FontId::proportional(12.0),
            theme::TEXT_DIM,
        );
        return;
    };
    let mut rows: Vec<(String, usize)> = Vec::new();
    for e in hx_enzymes::ENZYMES.iter() {
        let n = hx_enzymes::find_site_positions(&src.seq, e.site.as_bytes(), src.circular).len();
        if n > 0 {
            rows.push((e.name.to_string(), n));
        }
    }
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect));
    egui::ScrollArea::vertical()
        .id_salt("cloning_enz")
        .show(&mut child, |ui| {
            for (name, n) in rows {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(name).size(11.0).strong());
                    ui.label(
                        RichText::new(if n == 1 {
                            "1 site".to_string()
                        } else {
                            format!("{n} sites")
                        })
                        .size(10.5)
                        .color(theme::TEXT_DIM),
                    );
                });
            }
        });
}

fn paint_primers(ui: &mut egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let Some(src) = current_source(dlg) else {
        ui.painter().text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "Please choose a source file.",
            egui::FontId::proportional(12.0),
            theme::TEXT_DIM,
        );
        return;
    };
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect));
    egui::ScrollArea::vertical()
        .id_salt("cloning_primers")
        .show(&mut child, |ui| {
            if src.primers.is_empty() {
                ui.label(
                    RichText::new(crate::i18n::text(ui.ctx(),"No primers in this sequence."))
                        .size(11.0)
                        .color(theme::TEXT_DIM),
                );
                return;
            }
            for p in &src.primers {
                let r = p.sites.first().map(|s| s.range);
                let text = match r {
                    Some(r) => format!(
                        "{}  ({}..{})",
                        p.name,
                        r.start + 1,
                        r.end.max(r.start + 1)
                    ),
                    None => p.name.clone(),
                };
                ui.label(RichText::new(text).size(11.0).color(theme::PRIMER));
            }
        });
}

fn paint_preview(ui: &egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let painter = ui.painter_at(rect);
    let (seq, circular, features, highlight, primers): (
        Vec<u8>,
        bool,
        Vec<(usize, usize, Color32)>,
        Option<(usize, usize)>,
        Vec<(usize, usize)>,
    ) = match dlg.tab_kind(dlg.tab) {
        TabKind::Product => match &dlg.preview {
            Some(p) => (
                p.seq.clone(),
                p.circular,
                p.features
                    .iter()
                    .flat_map(|f| {
                        f.segments.iter().map(move |s| {
                            (
                                s.range.start,
                                s.range.end,
                                parse_color(&s.color).unwrap_or(theme::FEATURE_FALLBACK),
                            )
                        })
                    })
                    .collect(),
                if dlg.method.is_gateway() {
                    p.spans.first().copied()
                } else {
                    p.spans
                        .get(1)
                        .copied()
                        .or_else(|| p.spans.first().copied())
                },
                Vec::new(),
            ),
            None => (Vec::new(), false, Vec::new(), None, Vec::new()),
        },
        _ => match current_source(dlg) {
            Some(s) => (
                s.seq.clone(),
                s.circular,
                s.features
                    .iter()
                    .flat_map(|f| {
                        f.segments.iter().map(move |g| {
                            (
                                g.range.start,
                                g.range.end,
                                parse_color(&g.color).unwrap_or(theme::FEATURE_FALLBACK),
                            )
                        })
                    })
                    .collect(),
                dlg.active_pick().region,
                s.primers
                    .iter()
                    .filter_map(|p| p.sites.first().map(|s| (s.range.start, s.range.end)))
                    .collect(),
            ),
            None => (Vec::new(), false, Vec::new(), None, Vec::new()),
        },
    };
    if seq.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "Please choose a source file.",
            egui::FontId::proportional(12.0),
            theme::TEXT_DIM,
        );
        return;
    }
    let n = seq.len();
    let highlight = highlight.filter(|(a, b)| b > a && *b <= n);
    if circular {
        let c = rect.center();
        let r = rect.width().min(rect.height()) * 0.40;
        painter.circle_stroke(c, r, Stroke::new(2.0, theme::BACKBONE));
        painter.circle_stroke(c, r * 0.92, Stroke::new(0.6, theme::RULER));
        for (a, b, col) in features {
            let (a, b) = (a.min(n), b.min(n));
            if b <= a || b - a > n / 2 {
                continue;
            }
            let t0 = a as f32 / n as f32 * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
            let t1 = b as f32 / n as f32 * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
            let steps = ((b - a) / 8).clamp(2, 48);
            let pts: Vec<Pos2> = (0..=steps)
                .map(|i| {
                    let t = t0 + (t1 - t0) * i as f32 / steps as f32;
                    Pos2::new(
                        c.x + r * 0.92 * t.cos(),
                        c.y + r * 0.92 * t.sin(),
                    )
                })
                .collect();
            painter.add(egui::Shape::line(pts, Stroke::new(5.0, col)));
        }
        for (a, _b) in primers {
            if a >= n {
                continue;
            }
            let t = a as f32 / n as f32 * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
            let p0 = Pos2::new(c.x + r * 1.03 * t.cos(), c.y + r * 1.03 * t.sin());
            let p1 = Pos2::new(c.x + r * 1.14 * t.cos(), c.y + r * 1.14 * t.sin());
            painter.line_segment([p0, p1], Stroke::new(1.0, theme::PRIMER));
        }
        if let Some((a, b)) = highlight {
            let t0 = a as f32 / n as f32 * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
            let t1 = b as f32 / n as f32 * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
            let steps = 24;
            let pts: Vec<Pos2> = (0..=steps)
                .map(|i| {
                    let t = t0 + (t1 - t0) * i as f32 / steps as f32;
                    Pos2::new(
                        c.x + r * 0.99 * t.cos(),
                        c.y + r * 0.99 * t.sin(),
                    )
                })
                .collect();
            painter.add(egui::Shape::line(
                pts,
                Stroke::new(7.0, Color32::from_rgb(0xff, 0x8a, 0x8a)),
            ));
        }
        painter.text(
            c,
            egui::Align2::CENTER_CENTER,
            format!("{n} bp"),
            egui::FontId::proportional(12.0),
            theme::TEXT,
        );
    } else {
        let y = rect.center().y;
        painter.line_segment(
            [Pos2::new(rect.left(), y), Pos2::new(rect.right(), y)],
            Stroke::new(1.6, theme::BACKBONE),
        );
        let x_of = |p: usize| rect.left() + rect.width() * (p as f32 / n.max(1) as f32);
        for (a, b, col) in features {
            let (a, b) = (a.min(n), b.min(n));
            if b <= a {
                continue;
            }
            painter.line_segment(
                [Pos2::new(x_of(a), y - 5.0), Pos2::new(x_of(b), y - 5.0)],
                Stroke::new(4.0, col),
            );
        }
        for (a, _b) in primers {
            if a >= n {
                continue;
            }
            painter.line_segment(
                [Pos2::new(x_of(a), y - 14.0), Pos2::new(x_of(a), y - 2.0)],
                Stroke::new(1.0, theme::PRIMER),
            );
        }
        if let Some((a, b)) = highlight {
            painter.line_segment(
                [Pos2::new(x_of(a), y + 5.0), Pos2::new(x_of(b), y + 5.0)],
                Stroke::new(5.0, Color32::from_rgb(0xff, 0x8a, 0x8a)),
            );
        }
        painter.text(
            Pos2::new(rect.left() + 4.0, rect.top() + 10.0),
            egui::Align2::LEFT_TOP,
            format!("{n} bp"),
            egui::FontId::proportional(11.0),
            theme::TEXT_DIM,
        );
    }
}

// ---------------------------------------------------------------------------
// 左栏：引物
// ---------------------------------------------------------------------------

fn primer_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    let source_kind = match dlg.tab_kind(dlg.tab) {
        TabKind::Fragment | TabKind::Fragments | TabKind::LinearizedVector => {
            if dlg.method == Method::GatewayBp || dlg.method == Method::GatewayBpLr(1) {
                "attB Fragment"
            } else {
                "Fragment"
            }
        }
        TabKind::Vector => "Vector",
        _ => "Fragment",
    };
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let gg = golden_gate::active(dlg);
    let pick = match dlg.tab_kind(dlg.tab) {
        TabKind::Fragment | TabKind::Fragments | TabKind::LinearizedVector => {
            &mut dlg.fragments[idx]
        }
        _ => &mut dlg.vector,
    };
    let _ = idx;
    egui::Frame::NONE
        .fill(Color32::from_rgb(0xf7, 0xf7, 0xf7))
        .stroke(Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::symmetric(8, 6))
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            primer_row(
                ui,
                "Primer 1:",
                &format!("{source_kind}.FOR"),
                &mut pick.primer1,
                &mut pick.phos1,
                "p1",
            );
            ui.add_space(if gg { 16.0 } else { 2.0 });
            primer_row(
                ui,
                "Primer 2:",
                &format!("{source_kind}.REV"),
                &mut pick.primer2,
                &mut pick.phos2,
                "p2",
            );
        });
}

fn primer_row(
    ui: &mut egui::Ui,
    label: &str,
    name: &str,
    tail: &mut String,
    phos: &mut bool,
    id: &str,
) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),label)).size(11.0));
        egui::ComboBox::from_id_salt(format!("primer_{id}"))
            .selected_text(format!("<{name}>"))
            .width(150.0)
            .show_ui(ui, |ui| {
                ui.label(format!("<{name}>"));
            });
        ui.add_space(6.0);
        ui.label(RichText::new("5′").size(11.0));
        ui.add(
            egui::TextEdit::singleline(tail)
                .desired_width((ui.available_width() - 80.0).max(120.0)),
        );
        ui.label(RichText::new("3′").size(11.0));
    });
    ui.horizontal(|ui| {
        ui.checkbox(phos, crate::i18n::text(ui.ctx(),"5′ Phosphorylated"));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.button(crate::i18n::text(ui.ctx(),"Reverse Complement")).clicked() {
                let mut b = tail.as_bytes().to_vec();
                for x in b.iter_mut() {
                    *x = match *x {
                        b'A' => b'T',
                        b'T' => b'A',
                        b'G' => b'C',
                        b'C' => b'G',
                        other => other,
                    };
                }
                b.reverse();
                *tail = String::from_utf8_lossy(&b).to_string();
            }
        });
    });
}

// ---------------------------------------------------------------------------
// 右栏：参数面板
// ---------------------------------------------------------------------------

fn side_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    if golden_gate::active(dlg) {
        match dlg.tab_kind(dlg.tab) {
            TabKind::Vector => golden_gate::panel(ui, dlg, true),
            TabKind::Fragment | TabKind::Fragments => golden_gate::panel(ui, dlg, false),
            _ => product_panel(ui, dlg),
        }
        return;
    }
    // 调用方（原生窗口右侧栏）已经套了 ScrollArea；这里只排版。
    ui.set_width((SIDE_W - 30.0).max(180.0));
    match dlg.tab_kind(dlg.tab) {
        TabKind::Vector => vector_panel(ui, dlg),
        TabKind::Fragment => fragment_panel(ui, dlg),
        TabKind::Fragments => fragments_panel(ui, dlg),
        TabKind::LinearizedVector => linearized_vector_panel(ui, dlg),
        TabKind::DonorVector => gateway_vector_panel(ui, dlg, false),
        TabKind::DestinationVector => gateway_vector_panel(ui, dlg, true),
        TabKind::Product => product_panel(ui, dlg),
    }
}

fn source_combo(ui: &mut egui::Ui, docs: &[DocSnap], pick: &mut Pick) {
    let label = if let Some(n) = &pick.std_name {
        n.clone()
    } else {
        pick.doc
            .and_then(|i| docs.get(i))
            .map(|d| d.name.clone())
            .unwrap_or_else(|| "< Please choose >".to_string())
    };
    egui::ComboBox::from_id_salt(format!("src_{:?}_{}", pick.doc, label))
        .selected_text(if pick.doc.is_none()&&pick.std_name.is_none(){crate::i18n::text(ui.ctx(),"< Please choose >").to_string()}else{label})
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in docs.iter().enumerate() {
                if ui.selectable_label(pick.doc == Some(i), &d.name).clicked() {
                    pick.cut_pair = None;
                    pick.doc = Some(i);
                    pick.std_name = None;
                    pick.std_src = None;
                }
            }
        });
}

fn vector_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Vector:")).size(11.5));
    let idx = dlg.vector.doc;
    let label = idx
        .and_then(|i| dlg.docs.get(i))
        .map(|d| d.name.clone())
        .unwrap_or_else(|| "< Please choose >".to_string());
    egui::ComboBox::from_id_salt("vector_source")
        .selected_text(if idx.is_none(){crate::i18n::text(ui.ctx(),"< Please choose >").to_string()}else{label})
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in dlg.docs.iter().enumerate() {
                if ui.selectable_value(&mut dlg.vector.doc, Some(i), &d.name).clicked() {
                    dlg.vector.cut_pair = None;
                    dlg.vector.std_name = None;
                    dlg.vector.std_src = None;
                }
            }
        });

    ui.add_space(6.0);
    if matches!(dlg.method, Method::GoldenGate | Method::GoldenGateMulti) {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Golden Gate enzyme:")).size(10.5));
        let names = dlg.gg_enzymes.clone();
        let selected = names
            .get(dlg.gg_enzyme.min(names.len().saturating_sub(1)))
            .copied()
            .unwrap_or("BsaI");
        egui::ComboBox::from_id_salt("gg_enzyme")
            .selected_text(selected)
            .width(SIDE_W - 34.0)
            .show_ui(ui, |ui| {
                for (i, e) in names.iter().enumerate() {
                    ui.selectable_value(&mut dlg.gg_enzyme, i, *e);
                }
            });
        ui.add_space(4.0);
    }
    if dlg.method.is_restriction() {
        restriction_controls(ui, dlg, true);
        return;
    }
    ui.radio_value(&mut dlg.vector.lin, Linearize::Pcr, crate::i18n::text(ui.ctx(),"Linearize by PCR"));
    ui.radio_value(
        &mut dlg.vector.lin,
        Linearize::Restriction,
        crate::i18n::text(ui.ctx(),"Linearize with restriction enzymes"),
    );
    if dlg.vector.lin == Linearize::Restriction {
        enzyme_combo(ui, &mut dlg.vector, 0, &dlg.enzyme_filter.clone(), true, 0);
    }
    ui.radio_value(
        &mut dlg.vector.lin,
        Linearize::Direct,
        crate::i18n::text(ui.ctx(),"Use directly as the linearized Vector"),
    );
    ui.add_space(5.0);
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Polymerase for making linearized Vector:")).size(11.0));
    egui::ComboBox::from_id_salt("pol_v")
        .selected_text(POLYMERASES[dlg.polymerase.min(POLYMERASES.len() - 1)])
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, p) in POLYMERASES.iter().enumerate() {
                ui.selectable_value(&mut dlg.polymerase, i, *p);
            }
        });
    region_control(ui, dlg, true);
    orientation(ui, dlg, true);
    source_stats(ui, dlg, true);
}

fn fragment_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {gateway_fragment_badges(ui,dlg,"attB Fragment:");}
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let src_label = if dlg.method.is_bp() {
        "Source of attB Fragment:"
    } else {
        "Source of Fragment:"
    };
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),src_label)).size(11.0));
    let label = dlg
        .fragments
        .get(idx)
        .and_then(|p| p.std_name.clone().or_else(|| p.doc.and_then(|i| dlg.docs.get(i)).map(|d| d.name.clone())))
        .unwrap_or_else(|| "< Please choose >".to_string());
    egui::ComboBox::from_id_salt("fragment_source")
        .selected_text(if dlg.fragments[idx].doc.is_none()&&dlg.fragments[idx].std_name.is_none(){crate::i18n::text(ui.ctx(),"< Please choose >").to_string()}else{label})
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in dlg.docs.iter().enumerate() {
                if ui
                    .selectable_value(&mut dlg.fragments[idx].doc, Some(i), &d.name)
                    .clicked()
                {
                    dlg.fragments[idx].cut_pair = None;
                    dlg.fragments[idx].std_name = None;
                    dlg.fragments[idx].std_src = None;
                }
            }
        });

    // 限制性克隆/连接：片段也用酶切，不显示 PCR/聚合酶。
    if dlg.method.is_restriction() || dlg.method.is_ligation() {
        ui.add_space(5.0);
        restriction_controls(ui, dlg, false);
        source_stats(ui, dlg, false);
        return;
    }

    ui.add_space(5.0);
    ui.radio_value(&mut dlg.fragments[idx].use_pcr, true, crate::i18n::text(ui.ctx(),"Use as a template for PCR"));
    ui.radio_value(
        &mut dlg.fragments[idx].use_pcr,
        false,
        crate::i18n::text(ui.ctx(),if dlg.method.is_bp() {
            "Use directly as the attB Fragment"
        } else {
            "Use directly as the Fragment"
        }),
    );
    ui.add_space(4.0);
    if dlg.method.is_bp() && !dlg.fragments[idx].use_pcr {
        // SnapGene 的 BP 直接模式：选完来源后右侧没有别的参数。
        return;
    }
    // TA/GC/TOPO：只有“当作 PCR 模板”时才需要聚合酶；直接使用时右侧
    // 只剩底部的片段方向（SnapGene 8.2.3 实拍）。
    if dlg.method.is_ta_topo() {
        if dlg.fragments[idx].use_pcr {
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Polymerase for making Fragment:")).size(11.0));
            let polymerases = if matches!(dlg.method, Method::TopoBlunt | Method::TopoDirectional) { POLYMERASES } else { TA_POLYMERASES };
            egui::ComboBox::from_id_salt("pol_ta")
                .selected_text(polymerases[dlg.polymerase.min(polymerases.len() - 1)])
                .width(SIDE_W - 34.0)
                .show_ui(ui, |ui| {
                    for (i, p) in polymerases.iter().enumerate() {
                        ui.selectable_value(&mut dlg.polymerase, i, *p);
                    }
                });
            ui.add_space(4.0);
        }
        ui.add_space(4.0);
        let pad = (ui.available_height() - 46.0).max(16.0);
        ui.add_space(pad);
        orientation(ui, dlg, false);
        return;
    }
    if matches!(dlg.method, Method::TaCloning | Method::GcCloning) {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Polymerase for making Fragment:")).size(11.0));
        egui::ComboBox::from_id_salt("pol_f")
            .selected_text(TA_POLYMERASES[dlg.polymerase.min(TA_POLYMERASES.len() - 1)])
            .width(SIDE_W - 34.0)
            .show_ui(ui, |ui| {
                for (i, p) in TA_POLYMERASES.iter().enumerate() {
                    ui.selectable_value(&mut dlg.polymerase, i, *p);
                }
            });
    } else if dlg.method.is_bp() {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Polymerase for making attB Fragment:")).size(10.5));
        egui::ComboBox::from_id_salt("pol_bp")
            .selected_text(POLYMERASES[dlg.polymerase.min(POLYMERASES.len() - 1)])
            .width(SIDE_W - 34.0)
            .show_ui(ui, |ui| {
                for (i, p) in POLYMERASES.iter().enumerate() {
                    ui.selectable_value(&mut dlg.polymerase, i, *p);
                }
            });
        // SnapGene 会按 attB 适配引物重新设计 Gateway 片段引物；HX 尚未实现
        // 这套适配逻辑，因此明确禁用，避免勾选后没有任何效果。
        let adapter_label = if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {
            "Use matching numbered attB adapter primers"
        } else {
            "Use attB1 and attB2 adapter primers"
        };
        ui.add_enabled(
            false,
            egui::Checkbox::new(&mut dlg.adapter_primers, crate::i18n::text(ui.ctx(),adapter_label)),
        )
        .on_disabled_hover_text(crate::i18n::text(ui.ctx(),"attB adapter primer design is not implemented yet."));
    } else if !matches!(dlg.method, Method::TopoTa | Method::TopoBlunt | Method::TopoDirectional) {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Polymerase for making Fragment:")).size(11.0));
        egui::ComboBox::from_id_salt("pol_f2")
            .selected_text(POLYMERASES[dlg.polymerase.min(POLYMERASES.len() - 1)])
            .width(SIDE_W - 34.0)
            .show_ui(ui, |ui| {
                for (i, p) in POLYMERASES.iter().enumerate() {
                    ui.selectable_value(&mut dlg.polymerase, i, *p);
                }
            });
    }
    region_control(ui, dlg, false);
    orientation(ui, dlg, false);
    source_stats(ui, dlg, false);
}

fn fragments_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    let n = dlg.fragments.len();
    if dlg.method.is_lr() {
        // SnapGene 的 LR 对话框不对 Entry Clone 做 PCR：多个 Entry Clone
        // 用 “Entry Clone: 1 2” 数字小标签切换，下面只显示当前一个下拉框。
        let n = dlg.fragments.len();
        if n > 1 {
            ui.horizontal(|ui| {
                ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Entry Clone:")).size(11.0));
                for i in 0..n {
                    if number_badge(ui, dlg.frag_index == i, i + 1) {
                        dlg.frag_index = i;
                    }
                }
            });
            ui.add_space(4.0);
        }
        let idx = dlg.frag_index.min(n.saturating_sub(1));
        let label = if n > 1 {
            crate::ui_format!(ui.ctx(),"Entry Clone {}:", idx + 1)
        } else {
            "Entry Clone:".to_string()
        };
        entry_clone_combo(ui, dlg, idx, &label);
        return;
    }
    ui.horizontal(|ui| {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Fragments:")).size(11.5));
        if ui.small_button("−").clicked() && n > 2 {
            dlg.fragments.pop();
            dlg.frag_index = dlg.frag_index.min(dlg.fragments.len() - 1);
        }
        if ui.small_button("+").clicked() && n < 8 {
            let doc = dlg.docs.get(1).or_else(|| dlg.docs.first()).map(|_| 1.min(dlg.docs.len().saturating_sub(1)));
            dlg.fragments.push(Pick { doc, ..Default::default() });
            dlg.frag_index = dlg.fragments.len() - 1;
        }
        for i in 0..n {
            let label = format!("{}", i + 1);
            if ui
                .selectable_label(dlg.frag_index == i, label)
                .clicked()
            {
                dlg.frag_index = i;
            }
        }
    });
    ui.add_space(4.0);
    let idx = dlg.frag_index.min(n.saturating_sub(1));
    let flabel = if dlg.method.is_lr() {
        if dlg.fragments.len() > 1 {
            format!("Entry Clone {}:", idx + 1)
        } else {
            "Entry Clone:".to_string()
        }
    } else {
        crate::ui_format!(ui.ctx(),"Source of Fragment {}:", idx + 1)
    };
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),&flabel)).size(11.0));
    let label = dlg.fragments[idx]
        .doc
        .and_then(|i| dlg.docs.get(i))
        .map(|d| d.name.clone())
        .unwrap_or_else(|| "< Please choose >".to_string());
    egui::ComboBox::from_id_salt("multi_source")
        .selected_text(label)
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in dlg.docs.iter().enumerate() {
                if ui.selectable_value(&mut dlg.fragments[idx].doc, Some(i), &d.name).clicked() {
                    dlg.fragments[idx].cut_pair = None;
                }
            }
        });
    if dlg.method.is_restriction() || dlg.method.is_ligation() {
        restriction_controls(ui, dlg, false);
    } else {
        ui.add_space(5.0);
        ui.radio_value(&mut dlg.fragments[idx].use_pcr, true, crate::i18n::text(ui.ctx(),"Use as a template for PCR"));
        ui.radio_value(&mut dlg.fragments[idx].use_pcr, false, crate::i18n::text(ui.ctx(),"Use directly as the Fragment"));
        region_control(ui, dlg, false);
        orientation(ui, dlg, false);
    }
}

fn entry_clone_combo(ui: &mut egui::Ui, dlg: &mut CloningDlg, idx: usize, label: &str) {
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),label)).size(11.0));
    let text = dlg.fragments[idx]
        .doc
        .and_then(|i| dlg.docs.get(i))
        .map(|d| d.name.clone())
        .unwrap_or_else(|| "< Please choose >".to_string());
    egui::ComboBox::from_id_salt(("entry_clone", idx))
        .selected_text(text)
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in dlg.docs.iter().enumerate() {
                ui.selectable_value(&mut dlg.fragments[idx].doc, Some(i), &d.name);
            }
        });
    ui.add_space(7.0);
}

fn linearized_vector_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    // SnapGene 的「Linearized Vector」页（8.2.3 实拍）：
    //
    //   ◉ Custom Linearized Vector:
    //   [ < Please choose a linear sequence >        ]
    //   ○ Commercial Linearized Vector:
    //   ┌───────────────────────────────────────────┐
    //   │ pBAD-TOPO®                                │
    //   │ ...                                       │
    //   └───────────────────────────────────────────┘
    //
    // 两个单选永远可见，下拉框与商品化列表也一直在（不随单选切换隐藏）。
    let mut custom = !dlg.vector.use_std;
    if ui
        .radio_value(
            &mut custom,
            true,
            RichText::new(crate::i18n::text(ui.ctx(),"Custom Linearized Vector:")).size(11.5),
        )
        .clicked()
    {
        dlg.vector.use_std = false;
        dlg.vector.std_name = None;
        dlg.vector.std_src = None;
    }
    let label = if dlg.vector.use_std {
        dlg.vector
            .std_name
            .clone()
            .unwrap_or_else(|| "< Please choose a linear sequence >".to_string())
    } else {
        dlg.vector
            .doc
            .and_then(|i| dlg.docs.get(i))
            .map(|d| d.name.clone())
            .unwrap_or_else(|| "< Please choose a linear sequence >".to_string())
    };
    egui::ComboBox::from_id_salt("linvec_source")
        .selected_text(label)
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, d) in dlg.docs.iter().enumerate() {
                let on = !dlg.vector.use_std && dlg.vector.doc == Some(i);
                if ui.selectable_label(on, &d.name).clicked() {
                    dlg.vector.doc = Some(i);
                    dlg.vector.use_std = false;
                    dlg.vector.std_name = None;
                    dlg.vector.std_src = None;
                }
            }
        });

    ui.add_space(7.0);
    let mut commercial = dlg.vector.use_std;
    if ui
        .radio_value(
            &mut commercial,
            true,
            RichText::new(crate::i18n::text(ui.ctx(),"Commercial Linearized Vector:")).size(11.5),
        )
        .clicked()
    {
        dlg.vector.use_std = true;
        dlg.vector.doc = None;
    }

    let h = (ui.available_height() - 30.0).clamp(220.0, 620.0);
    egui::Frame::NONE
        .stroke(Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::same(3))
        .show(ui, |ui| {
            ui.set_min_height(h);
            egui::ScrollArea::vertical()
                .id_salt("std_linvec")
                .max_height(h)
                .show(ui, |ui| {
                    if dlg.std_ta.is_empty() {
                        ui.label(
                            RichText::new(crate::i18n::text(ui.ctx(),"No standard vectors found."))
                                .size(10.5)
                                .color(theme::TEXT_DIM),
                        );
                    }
                    for e in &dlg.std_ta {
                        let on = dlg.vector.use_std
                            && dlg.vector.std_name.as_deref() == Some(e.name.as_str());
                        if ui.selectable_label(on, &e.name).clicked() {
                            let name = e.name.clone();
                            if !select_std(&mut dlg.vector, e) {
                                dlg.hint = format!("Could not read {name}.");
                            }
                        }
                    }
                });
        });
}

fn gateway_vector_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg, dest: bool) {
    if !dest && matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {
        gateway_fragment_badges(ui,dlg,"Donor for Fragment:");
        let i=dlg.frag_index.min(dlg.gateway_donors.len()-1);
        if let Some(&(a,ad,b,bd))=cloning::gateway::layout(dlg.fragments.len()).and_then(|v|v.get(i)) {
            ui.label(crate::ui_format!(ui.ctx(),"Required donor: attP{} {} / attP{} {}",a,if ad==1 {"→"}else{"←"},b,if bd==1 {"→"}else{"←"}));
        }
        if let Some(error)=std_chooser(ui,&dlg.docs,&mut dlg.gateway_donors[i],&dlg.std_donor,"Donor Vector") {dlg.hint=error;}
        return;
    }
    let title = if dest { "Destination Vector" } else { "Donor Vector" };
    let bp_lr = matches!(dlg.method, Method::GatewayBpLr(_));
    let hint = if dest && bp_lr {
        std_chooser(
            ui,
            &dlg.docs,
            &mut dlg.destination,
            &dlg.std_dest,
            title,
        )
    } else {
        let list = if dest {
            dlg.std_dest.clone()
        } else {
            dlg.std_donor.clone()
        };
        std_chooser(ui, &dlg.docs, &mut dlg.vector, &list, title)
    };
    if let Some(h) = hint {
        dlg.hint = h;
    }
}

/// 自定义/商品化来源二选一。返回 `Some(hint)` 表示读取商品化载体失败。
fn std_chooser(
    ui: &mut egui::Ui,
    docs: &[DocSnap],
    pick: &mut Pick,
    list: &[StdEntry],
    title: &str,
) -> Option<String> {
    let mut custom = !pick.use_std;
    let mut err = None;
    ui.horizontal(|ui| {
        ui.radio_value(&mut custom, true, crate::ui_format!(ui.ctx(),"Custom {}:",crate::i18n::text(ui.ctx(),title)));
    });
    if custom {
        pick.use_std = false;
        pick.std_name = None;
        pick.std_src = None;
        source_combo(ui, docs, pick);
    }
    ui.add_space(4.0);
    let mut standard = !custom;
    if ui
        .radio_value(&mut standard, true, crate::ui_format!(ui.ctx(),"Standard {}:",crate::i18n::text(ui.ctx(),title)))
        .clicked()
    {
        custom = false;
    }
    // SnapGene 始终显示标准载体列表（即使当前选的是自定义来源）。
    egui::Frame::NONE
        .stroke(Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::same(3))
        .show(ui, |ui| {
            ui.set_min_height(230.0);
            egui::ScrollArea::vertical()
                .id_salt(format!("std_{title}"))
                .max_height(230.0)
                .show(ui, |ui| {
                    if list.is_empty() {
                        ui.label(
                            RichText::new(crate::i18n::text(ui.ctx(),"No standard vectors found."))
                                .size(10.5)
                                .color(theme::TEXT_DIM),
                        );
                    }
                    for e in list {
                        let on = pick.use_std && pick.std_name.as_deref() == Some(e.name.as_str());
                        if ui.selectable_label(on, &e.name).clicked() {
                            custom = false;
                            pick.use_std = true;
                            if !select_std(pick, e) {
                                err = Some(format!("Could not read {}.", e.name));
                            }
                        }
                    }
                });
        });
    if !custom {
        pick.use_std = true;
    }
    let _ = standard;
    err
}

fn restriction_controls(ui: &mut egui::Ui, dlg: &mut CloningDlg, is_vector: bool) {
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    ui.horizontal(|ui| {
        ui.label(crate::i18n::text(ui.ctx(),"Enzyme filter:"));
        ui.add(
            egui::TextEdit::singleline(&mut dlg.enzyme_filter)
                .hint_text(crate::i18n::text(ui.ctx(),"e.g. EcoRI"))
                .desired_width(150.0),
        );
    });
    let filter = dlg.enzyme_filter.clone();
    if is_vector {
        ui.horizontal(|ui| {
            ui.label(crate::i18n::text(ui.ctx(),"Cut"));
            enzyme_combo(ui, &mut dlg.vector, 0, &filter, true, 0);
            ui.label(crate::i18n::text(ui.ctx(),"with"));
            enzyme_combo(ui, &mut dlg.vector, 1, &filter, true, 0);
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut dlg.vector.blunt, crate::i18n::text(ui.ctx(),"and then blunt"));
        });
    } else {
        ui.horizontal(|ui| {
            ui.label(crate::i18n::text(ui.ctx(),"Cut"));
            enzyme_combo(ui, &mut dlg.fragments[idx], 0, &filter, false, idx + 1);
            ui.label(crate::i18n::text(ui.ctx(),"with"));
            enzyme_combo(ui, &mut dlg.fragments[idx], 1, &filter, false, idx + 1);
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut dlg.fragments[idx].blunt, crate::i18n::text(ui.ctx(),"and then blunt"));
        });
    }
    ui.add_space(5.0);
    ui.label(crate::i18n::text(ui.ctx(),"Site of insertion:"));
    let mut choices = Vec::new();
    let mut selected = None;
    let mut names = String::new();
    let mut n = 1i64;
    let mut circular = false;
    let mut empty = "Choose a source file.".to_string();
    {
        let pick = if is_vector {
            &dlg.vector
        } else {
            &dlg.fragments[idx]
        };
        if let Some(doc) = dlg.pick_source(pick) {
            names = enzyme_names(pick);
            if pick.enzyme.is_none() && pick.enzyme2.is_none() {
                empty = if dlg.method.is_ligation() {
                    "Using the source's existing ends. Choose an enzyme to digest it first.".to_string()
                } else { "Choose an enzyme to see the insertion sites.".to_string() };
            } else {
                empty = format!("No {} site in {}.", names, dlg.source_name(pick));
            }
            if !double_cuts_of(doc, pick).is_empty() {
                choices = digest_choices(doc, pick);
                selected = selected_digest_pair(doc, pick, is_vector);
                n = doc.seq.len().max(1) as i64;
                circular = doc.circular;
            }
        }
    }
    egui::Frame::NONE
        .stroke(Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::same(4))
        .show(ui, |ui| {
            ui.set_min_height(170.0);
            ui.set_width(SIDE_W - 32.0);
            if choices.is_empty() {
                ui.label(RichText::new(&empty).size(10.5).color(theme::TEXT_DIM));
            } else {
                if let Some((mut left, mut right)) = selected {
                    for (label, value) in [("From:", &mut left), ("To:", &mut right)] {
                        ui.horizontal(|ui| {
                            ui.label(label);
                            egui::ComboBox::from_id_salt(("digest_boundary", is_vector, idx, label))
                                .selected_text(format!("{} / {}", value.top,
                                    if circular { value.bottom.rem_euclid(n) } else { value.bottom }))
                                .width(180.0)
                                .show_ui(ui, |ui| {
                                    for cut in &choices {
                                        ui.selectable_value(value, *cut, format!("{names}   {} / {}", cut.top,
                                            if circular { cut.bottom.rem_euclid(n) } else { cut.bottom }));
                                    }
                                });
                        });
                    }
                    let pick = if is_vector { &mut dlg.vector } else { &mut dlg.fragments[idx] };
                    pick.cut_pair = Some((left, right));
                    ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Top / bottom strand cut coordinates")).size(10.5).color(theme::TEXT_DIM));
                } else {
                    ui.label(RichText::new(crate::i18n::text(ui.ctx(),"No double-stranded fragment between these cuts.")).size(10.5).color(theme::TEXT_DIM));
                }
            }
        });
    if is_vector {
        orientation(ui, dlg, true);
    } else {
        orientation(ui, dlg, false);
    }
}

fn region_control(ui: &mut egui::Ui, dlg: &mut CloningDlg, is_vector: bool) {
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let n = {
        let pick = if is_vector { &dlg.vector } else { &dlg.fragments[idx] };
        dlg.pick_source(pick).map(|s| s.seq.len()).unwrap_or(0)
    };
    if n == 0 {
        return;
    }
    let pick = if is_vector { &mut dlg.vector } else { &mut dlg.fragments[idx] };
    let (mut a, mut b) = match pick.region {
        Some((a, b)) => (a + 1, b),
        None => (1, n),
    };
    let mut changed = false;
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label(crate::i18n::text(ui.ctx(),"Region:"));
        changed |= ui
            .add(egui::DragValue::new(&mut a).range(1..=n).speed(1.0))
            .changed();
        ui.label("..");
        changed |= ui
            .add(egui::DragValue::new(&mut b).range(1..=n).speed(1.0))
            .changed();
    });
    if changed {
        let (a, b) = (a.min(n), b.min(n));
        pick.region = if a <= 1 && b >= n {
            None
        } else if a <= b {
            Some((a - 1, b))
        } else {
            Some((b - 1, a))
        };
    }
}

fn orientation(ui: &mut egui::Ui, dlg: &mut CloningDlg, is_vector: bool) {
    ui.add_space(5.0);
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let pick = if is_vector {
        &mut dlg.vector
    } else {
        &mut dlg.fragments[idx]
    };
    ui.horizontal(|ui| {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Orientation:")).size(10.5));
        if ui.button(RichText::new("→").size(12.0)).clicked() {
            pick.reversed = false;
        }
        if ui.button(RichText::new("←").size(12.0)).clicked() {
            pick.reversed = true;
        }
    });
}

fn source_stats(ui: &mut egui::Ui, dlg: &CloningDlg, is_vector: bool) {
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let pick = if is_vector {
        &dlg.vector
    } else {
        &dlg.fragments[idx]
    };
    if let Some(src) = dlg.pick_source(pick) {
        ui.add_space(5.0);
        let mut text = format!("{} bp", src.seq.len());
        if let Some((a, b)) = pick.region {
            let n = b.saturating_sub(a);
            text.push_str(&format!(" · region {n} bp"));
        }
        if src.circular {
            text.push_str(" · circular");
        }
        ui.label(RichText::new(text).size(10.5).color(theme::TEXT_DIM));
    }
}

fn enzyme_combo(
    ui: &mut egui::Ui,
    pick: &mut Pick,
    row: usize,
    filter: &str,
    is_vector: bool,
    salt: usize,
) {
    let sel = [pick.enzyme, pick.enzyme2][row]
        .and_then(|i| hx_enzymes::ENZYMES.get(i))
        .map(|e| e.name.to_string())
        .unwrap_or_else(|| "<none>".to_string());
    egui::ComboBox::from_id_salt(format!(
        "enz_{}_{}_{}",
        if is_vector { "v" } else { "f" },
        salt,
        row
    ))
    .selected_text(sel)
    .width(92.0)
    .show_ui(ui, |ui| {
        let needle = filter.to_ascii_lowercase();
        egui::ScrollArea::vertical()
            .max_height(240.0)
            .id_salt(format!("enzsc{is_vector}{salt}{row}"))
            .show(ui, |ui| {
                if ui
                    .selectable_label([pick.enzyme, pick.enzyme2][row].is_none(), "<none>")
                    .clicked()
                {
                    if row == 0 {
                        pick.enzyme = None;
                    } else {
                        pick.enzyme2 = None;
                    }
                    pick.cut_pair = None;
                }
                for (i, e) in hx_enzymes::ENZYMES.iter().enumerate() {
                    if !needle.is_empty() && !e.name.to_ascii_lowercase().contains(&needle) {
                        continue;
                    }
                    let on = [pick.enzyme, pick.enzyme2][row] == Some(i);
                    if ui.selectable_label(on, e.name).clicked() {
                        pick.cut_pair = None;
                        if row == 0 {
                            pick.enzyme = Some(i);
                        } else {
                            pick.enzyme2 = Some(i);
                        }
                    }
                }
            });
    });
}

fn product_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) && dlg.tab==2 {
        gateway_fragment_badges(ui,dlg,"Entry Clone:");
        match dlg.build_gateway_entry(dlg.frag_index) {
            Ok(entry)=>{ui.label(crate::ui_format!(ui.ctx(),"Entry Clone {} · {} bp",dlg.frag_index+1,entry.display_len()));},
            Err(error)=>{ui.label(error);},
        }
        ui.label(crate::i18n::text(ui.ctx(),"This intermediate is generated by the matching attB Fragment and Donor Vector BP reaction."));
        return;
    }
    if golden_gate::active(dlg) {
        if ui.add_sized([ui.available_width(), 26.0], egui::Button::new(crate::i18n::text(ui.ctx(),"Choose PCR Primers…"))).clicked() {
            if let Err(error) = golden_gate::design(dlg) { dlg.hint = error; }
            dlg.preview_key.clear();
        }
        if ui.add_sized([ui.available_width(), 26.0], egui::Button::new(crate::i18n::text(ui.ctx(),"Adjust Overhangs…"))).clicked() { dlg.gg_adjust_open = true; }
        ui.add_space(24.0);
    }
    ui.add_space(2.0);
    // 限制性克隆/连接不需要设计 PCR 引物。
    if !(golden_gate::active(dlg) || dlg.method.is_restriction() || dlg.method.is_ligation() || matches!(dlg.method, Method::GatewayLr(_))) {
        let btn = if dlg.method.is_bp() {
            "Choose attB PCR Primers…"
        } else if dlg.method.is_ta_topo() {
            "Choose PCR Primers…"
        } else {
            "Choose Overlapping PCR Primers…"
        };
        if ui
            .add_sized(egui::vec2(SIDE_W - 34.0, 22.0), egui::Button::new(btn))
            .clicked()
        {
            design_primers(dlg);
        }
    }
    ui.add_space(5.0);
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Bacterial Transformation Strain")).size(10.5));
    egui::ComboBox::from_id_salt("strain")
        .selected_text(STRAINS[dlg.strain.min(STRAINS.len() - 1)])
        .width(SIDE_W - 34.0)
        .show_ui(ui, |ui| {
            for (i, s) in STRAINS.iter().enumerate() {
                ui.selectable_value(&mut dlg.strain, i, *s);
            }
        });
    ui.add_space(7.0);
    let sections = product_sections(dlg);
    egui::Frame::NONE
        .fill(Color32::from_rgb(0xf7, 0xf7, 0xf7))
        .stroke(Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::same(6))
        .show(ui, |ui| {
            ui.set_width(if golden_gate::active(dlg) { ui.available_width() } else { SIDE_W - 46.0 });
            ui.set_min_height(250.0);
            for (i, row) in sections.iter().enumerate() {
                if i > 0 {
                    ui.add_space(4.0);
                }
                ui.horizontal(|ui| {
                    ui.label(RichText::new(&row.title).size(11.5).strong());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if let Some(n) = row.size {
                            ui.label(
                                RichText::new(format!("{n} bp"))
                                    .size(10.5)
                                    .color(theme::TEXT_DIM),
                            );
                        }
                    });
                });
                if let Some(source) = &row.source {
                    ui.label(
                        RichText::new(source)
                            .size(10.5)
                            .color(theme::TEXT_DIM),
                    );
                }
                if let Some(file) = &row.file {
                    // SnapGene 的 `Make file:` 会在克隆时把中间产物写成文件。
                    // HX 尚未实现写文件的位置/覆盖策略，因此这里明确禁用，
                    // 不显示一个可勾选但不产生文件的假控件。
                    ui.horizontal(|ui| {
                        let mut on = false;
                        ui.add_enabled(
                            false,
                            egui::Checkbox::new(&mut on, RichText::new(crate::i18n::text(ui.ctx(),"Make file:")).size(10.5)),
                        )
                        .on_disabled_hover_text(
                            crate::i18n::text(ui.ctx(),"Writing intermediate files is not implemented yet."),
                        );
                        ui.label(RichText::new(file).size(10.5).color(theme::TEXT_DIM));
                    });
                }
            }
            // Gateway 的右侧卡片只有来源信息（产物大小在绿色提示框里）；
            // TA/GC/TOPO 也一样（SnapGene 的 Product 页没有产物摘要行）。
            if !dlg.method.is_gateway() && !dlg.method.is_ta_topo() {
                ui.separator();
                match &dlg.preview {
                    Some(p) => {
                        ui.label(
                            RichText::new(crate::ui_format!(ui.ctx(),"Product: {} bp", p.display_len()))
                                .size(11.0)
                                .strong(),
                        );
                        for (i, j) in p.junctions.iter().enumerate() {
                            ui.label(
                                RichText::new(format!(
                                    "Junction {}: {} bp at {}",
                                    i + 1,
                                    j.overlap,
                                    if p.circular { j.at % p.seq.len().max(1) } else {
                                        j.at.saturating_sub(p.ends.upstream_stickiness.min(0).unsigned_abs() as usize)
                                    }
                                ))
                                .size(10.0)
                                .color(theme::TEXT_DIM),
                            );
                        }
                    }
                    None => {
                        ui.label(
                            RichText::new(crate::i18n::text(ui.ctx(),"No product yet — choose both sources."))
                                .size(10.5)
                                .color(theme::TEXT_DIM),
                        );
                    }
                }
            }
        });
}

/// 右侧「产物组成」卡片的一行，照 SnapGene 的排版：
/// 标题 + 右对齐长度，来源一行，自定义文件才有 `Make file:`。
struct CardRow {
    title: String,
    /// 来源行。商品化载体只写载体名（SnapGene 的 `pDONR™221`），
    /// 自定义文件写 `Source: 文件名`。
    source: Option<String>,
    /// `Make file:` 的建议文件名；商品化载体没有这一行。
    file: Option<String>,
    size: Option<usize>,
}

fn card_row(
    dlg: &CloningDlg,
    title: &str,
    pick: &Pick,
    file: &str,
    size: Option<usize>,
) -> CardRow {
    let name = dlg.source_name(pick);
    let (source, file) = if pick.use_std && pick.std_name.is_some() {
        (name, None)
    } else {
        (format!("Source: {name}"), Some(file.to_string()))
    };
    CardRow {
        title: title.into(),
        source: Some(source),
        file,
        size,
    }
}

fn product_sections(dlg: &CloningDlg) -> Vec<CardRow> {
    let size = |p: &Pick| dlg.pick_source(p).map(|s| s.display_len());
    match dlg.method {
        Method::Restriction | Method::RestrictionMulti => {
            let processed_size = |pick: &Pick, vector: bool| {
                dlg.pick_source(pick).and_then(|src| restriction_fragment(src, pick, vector, true).ok())
                    .map(|f| f.display_len())
            };
            let mut rows = vec![card_row(dlg, "Linearized Vector", &dlg.vector,
                "Linearized Vector.dna", processed_size(&dlg.vector, true))];
            for (i, pick) in dlg.fragments.iter().enumerate() {
                let title = if dlg.method.is_multiple() { format!("Fragment {}", i + 1) } else { "Fragment".into() };
                rows.push(card_row(dlg, &title, pick, &format!("{title}.dna"), processed_size(pick, false)));
            }
            rows
        }
        Method::Ligation => dlg.fragments.iter().enumerate().map(|(i, pick)| {
            let size = dlg.pick_source(pick).and_then(|src| restriction_fragment(src, pick, false, false).ok())
                .map(|f| f.display_len());
            let title = format!("Fragment {}", i + 1);
            card_row(dlg, &title, pick, &format!("{title}.dna"), size)
        }).collect(),
        Method::GatewayBp => vec![
            card_row(
                dlg,
                "attB Fragment",
                &dlg.fragments[0],
                "attB Fragment.dna",
                size(&dlg.fragments[0]),
            ),
            card_row(
                dlg,
                "Donor Vector",
                &dlg.vector,
                "Donor Vector.dna",
                size(&dlg.vector),
            ),
        ],
        Method::GatewayLr(_) => {
            let mut out = Vec::new();
            for (i, p) in dlg.fragments.iter().enumerate() {
                out.push(card_row(
                    dlg,
                    &format!("Entry Clone {}", i + 1),
                    p,
                    &format!("Entry Clone {}.dna", i + 1),
                    size(p),
                ));
            }
            out.push(card_row(
                dlg,
                "Destination Vector",
                &dlg.vector,
                "Destination Vector.dna",
                size(&dlg.vector),
            ));
            out
        }
        Method::GatewayBpLr(n) if n>1 => {
            let mut rows=Vec::new();
            for (i,pick) in dlg.fragments.iter().enumerate() {
                rows.push(card_row(dlg,&format!("attB Fragment {}",i+1),pick,&format!("attB Fragment {}.dna",i+1),size(pick)));
                let donor=&dlg.gateway_donors[i];
                rows.push(card_row(dlg,&format!("Donor Vector {}",i+1),donor,&format!("Donor Vector {}.dna",i+1),size(donor)));
                rows.push(CardRow{title:format!("Entry Clone {}",i+1),source:None,file:None,size:dlg.build_gateway_entry(i).ok().map(|p|p.display_len())});
            }
            rows.push(card_row(dlg,"Destination Vector",&dlg.destination,"Destination Vector.dna",size(&dlg.destination)));
            rows
        }
        Method::GatewayBpLr(_) => vec![
            card_row(
                dlg,
                "attB Fragment",
                &dlg.fragments[0],
                "attB Fragment.dna",
                size(&dlg.fragments[0]),
            ),
            card_row(
                dlg,
                "Donor Vector",
                &dlg.vector,
                "Donor Vector.dna",
                size(&dlg.vector),
            ),
            CardRow {
                title: "Entry Clone".into(),
                source: None,
                file: None,
                size: dlg.build_gateway_bp().ok().map(|p| p.display_len()),
            },
            card_row(
                dlg,
                "Destination Vector",
                &dlg.destination,
                "Destination Vector.dna",
                size(&dlg.destination),
            ),
        ],
        m if m.is_ta_topo() => {
            // SnapGene 的 TA/GC/TOPO Product 页：Fragment 在前、载体在后，
            // 长度用显示长度（顶链），不是文件里存的标记碱基长度。
            let top = |p: &Pick| dlg.pick_source(p).map(|s| s.display_len());
            let fragment_size=dlg.pick_source(&dlg.fragments[0])
                .and_then(|s|end_cloning_fragment(s,&dlg.fragments[0],dlg.method).ok())
                .map(|s|s.display_len());
            vec![
                card_row(
                    dlg,
                    "Fragment",
                    &dlg.fragments[0],
                    "Fragment.dna",
                    fragment_size,
                ),
                card_row(
                    dlg,
                    "Linearized Vector",
                    &dlg.vector,
                    "Linearized Vector.dna",
                    top(&dlg.vector),
                ),
            ]
        }
        _ => vec![
            card_row(
                dlg,
                "Linearized Vector",
                &dlg.vector,
                "Linearized Vector.dna",
                size(&dlg.vector),
            ),
            card_row(
                dlg,
                "Fragment",
                &dlg.fragments[0],
                "Fragment.dna",
                size(&dlg.fragments[0]),
            ),
        ],
    }
}

/// 自动设计 PCR 引物：片段两端加上载体两端同源臂。
fn design_primers(dlg: &mut CloningDlg) {
    if golden_gate::active(dlg) {
        if let Err(error) = golden_gate::design(dlg) { dlg.hint = error; }
        dlg.preview_key.clear();
        return;
    }
    if dlg.overlap_method() {
        if let Err(error) = design_homology_primers(dlg) { dlg.hint = error; return; }
        dlg.preview_key.clear(); dlg.refresh();
        return;
    }
    if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {
        let original=dlg.fragments.clone();let selected=dlg.frag_index;
        for i in 0..dlg.fragments.len() {
            if !dlg.fragments[i].use_pcr {continue;}
            dlg.frag_index=i;
            if let Err(error)=design_end_primers(dlg) {dlg.fragments=original;dlg.frag_index=selected;dlg.hint=error;return;}
        }
        dlg.frag_index=selected;
    } else if let Err(error)=design_end_primers(dlg) { dlg.hint=error;return; }
    dlg.preview_key.clear();
    dlg.refresh();
}

fn design_end_primers(dlg:&mut CloningDlg)->Result<(),String> {
    let idx=if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {dlg.frag_index.min(dlg.fragments.len()-1)}else{0};
    let pick=&dlg.fragments[idx];
    let f=dlg.pick_source(pick).ok_or("Please choose a Fragment source.")?;
    let (a,b)=pick.region.or_else(|| (!f.circular).then(|| f.top_range()))
        .ok_or("Please select the Fragment region to amplify.")?;
    if a>=f.seq.len() || b>f.seq.len() || (!f.circular && a>=b) { return Err("Please select a valid amplified region.".into()); }
    let body=f.slice(a,b,pick.reversed);
    if body.seq.len()<40 { return Err("Select at least 40 bases for automatic primer design.".into()); }
    let mut fwd=body.seq[..20].to_vec();
    let mut rev=hx_algo::reverse_complement(&body.seq[body.seq.len()-20..]);
    if dlg.method.is_bp() {
        let (left,reverse_right)=cloning::gateway::primer_sites(if matches!(dlg.method,Method::GatewayBpLr(_)) {dlg.fragments.len()}else{1},idx)?;
        if !body.seq.starts_with(&left) { fwd=[left,fwd].concat(); }
        if !body.seq.ends_with(&hx_algo::reverse_complement(&reverse_right)) { rev=[reverse_right,rev].concat(); }
    } else if dlg.method==Method::TopoDirectional && !body.seq.starts_with(b"CACC") {
        fwd=[b"CACC".to_vec(),fwd].concat();
    }
    cloning::homology::amplify(f,(a,b),pick.reversed,&fwd,&rev)?;
    dlg.fragments[idx].primer1=String::from_utf8_lossy(&fwd).into_owned();
    dlg.fragments[idx].primer2=String::from_utf8_lossy(&rev).into_owned();
    dlg.fragments[idx].use_pcr=true;
    dlg.fragments[idx].region=Some((a,b));
    Ok(())
}

fn gateway_fragment_badges(ui:&mut egui::Ui,dlg:&mut CloningDlg,label:&str) {
    ui.horizontal(|ui|{
        ui.label(crate::i18n::text(ui.ctx(),label));
        for i in 0..dlg.fragments.len() {if number_badge(ui,dlg.frag_index==i,i+1){dlg.frag_index=i;}}
    });ui.add_space(6.0);
}

fn design_homology_primers(dlg: &mut CloningDlg) -> Result<(), String> {
    let count = if dlg.method.is_multiple() { dlg.fragments.len() } else { 1 };
    let picks: Vec<_> = std::iter::once(&dlg.vector).chain(dlg.fragments.iter().take(count)).collect();
    let mut bodies = Vec::new();
    let mut pcr = Vec::new();
    for (i, pick) in picks.iter().enumerate() {
        let source = dlg.pick_source(pick).ok_or("Please choose all source files.")?;
        let amplify = if i == 0 { pick.lin == Linearize::Pcr } else { pick.use_pcr };
        let body = if amplify {
            let (a, b) = pick.region.ok_or("Select the region to amplify on every PCR template.")?;
            if a >= source.seq.len() || b > source.seq.len() || (!source.circular && a >= b) {
                return Err("Select a valid amplified region.".into());
            }
            source.slice(a, b, pick.reversed)
        } else if i == 0 { homology_vector(source, pick)? } else { homology_fragment(source, pick)? };
        bodies.push(body); pcr.push(amplify);
    }
    // 同源臂长度来自设置；不能再用固定的 40 bp 门槛拒绝本来可行的短片段。
    // 同时要保证两端各有 20 bp 退火区和一段非重叠核心，否则设计出的引物
    // 会互相覆盖或让拼接阶段报「overlaps consume all of Fragment」。
    let configured_arm = dlg.homology_overlap();
    let shortest = bodies.iter().map(|b| b.seq.len()).min().unwrap_or(0);
    let minimum = (configured_arm * 2 + 1).max(20);
    if shortest < minimum {
        return Err(format!(
            "Select at least {minimum} bases for automatic primer design with {configured_arm} bp overlaps."
        ));
    }
    let arm = configured_arm.min(shortest);
    let mut designs = Vec::new();
    for i in 0..bodies.len() {
        if !pcr[i] { continue; }
        let prev = (i + bodies.len() - 1) % bodies.len();
        let next = (i + 1) % bodies.len();
        // One tail creates each junction: add the previous fragment's suffix to
        // this forward primer, except when the previous fragment alone is PCR.
        let f_tail = &bodies[prev].seq[bodies[prev].seq.len() - arm..];
        let r_tail = if !pcr[next] { &bodies[next].seq[..arm] } else { &[] };
        let forward = [f_tail, &bodies[i].seq[..20]].concat();
        let reverse = hx_algo::reverse_complement(&[&bodies[i].seq[bodies[i].seq.len() - 20..], r_tail].concat());
        designs.push((i, String::from_utf8_lossy(&forward).into_owned(), String::from_utf8_lossy(&reverse).into_owned()));
    }
    // Validate the whole graph before changing any visible fields.
    let mut check = bodies.clone();
    for (i, f, r) in &designs {
        let temp = Source::new("PCR", bodies[*i].seq.clone(), false).with_features(bodies[*i].features.clone());
        check[*i] = cloning::homology::amplify(&temp, (0, temp.seq.len()), false, f.as_bytes(), r.as_bytes())?;
    }
    cloning::homology::assemble_circular(&check, dlg.min_overlap())?;
    for (i, f, r) in designs {
        let pick = if i == 0 { &mut dlg.vector } else { &mut dlg.fragments[i - 1] };
        pick.primer1 = f; pick.primer2 = r;
    }
    Ok(())
}

fn revcomp_bases(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&x| match x.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            o => o,
        })
        .collect()
}

fn has_motif(seq: &[u8], motif: &[u8]) -> bool {
    !motif.is_empty()
        && seq
            .windows(motif.len())
            .any(|w| w.eq_ignore_ascii_case(motif))
}

// ---------------------------------------------------------------------------
// 底部示意图
// ---------------------------------------------------------------------------

mod reaction_scheme;

fn paint_schematic(ui: &egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let painter = ui.painter_at(rect);
    painter.rect_stroke(
        rect,
        0.0,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    if dlg.method.is_ta_topo() {
        paint_ta_topo_schematic(&painter, rect, dlg);
        return;
    }
    // BP 的专用重组几何已有 SG 实拍校准；保留它，未就绪时整体淡化，
    // 避免为了修复通用图而回退已校准的 BP 布局。
    if dlg.method == Method::GatewayBp {
        let ready = dlg.preview.is_some() && dlg.primers_ready();
        let mut p = painter.clone();
        p.set_opacity(if ready { 1.0 } else { 0.38 });
        paint_gateway_bp_schematic(&p, rect);
        if !ready {
            painter.text(rect.left_top() + Vec2::new(8.0, 4.0), egui::Align2::LEFT_TOP,
                "BP reaction scheme · product not generated", egui::FontId::proportional(10.0), theme::TEXT_DIM);
        }
        return;
    }
    reaction_scheme::draw(ui, rect, dlg);
}

/// Gateway BP 的拼接示意图，照 SnapGene 的 Entry Clone 页：
///
/// ```text
///   attB1 ->        <- attB2
///   [======= 红色片段 =======]
///       X                 X
///   ╭───────────────────────╮
///   │ attP1 ->    <- attP2  │   →   ╭────────────────╮
///   ╰───────────────────────╯       │ attL1 -> <- attL2│
///                                   ╰────────────────╯
/// ```
///
/// 供体和产物都画成圆角胶囊（不是圆环）：SnapGene 把载体画成“一段
/// 双链”，重组位点落在顶边上；attB 片段悬在上方，用 X 表示交换。
fn paint_gateway_bp_schematic(painter: &egui::Painter, rect: Rect) {
    let red = Color32::from_rgb(0xe8, 0x20, 0x2a);
    let ink = Color32::from_rgb(0x1a, 0x1a, 0x1a);
    let font = egui::FontId::proportional(10.0);
    let frag_cx = rect.left() + rect.width() * 0.36;
    let prod_cx = rect.left() + rect.width() * 0.63;
    let mid_y = rect.center().y + 7.0;

    let pill_w = 130.0;
    let donor = Rect::from_center_size(
        Pos2::new(frag_cx, mid_y + 6.0),
        egui::vec2(pill_w, 46.0),
    );
    let product = Rect::from_center_size(
        Pos2::new(prod_cx, mid_y + 5.0),
        egui::vec2(pill_w, 52.0),
    );

    // ---- 上方：线性 attB 片段 ----
    let bar_w = 68.0;
    let bar_y = donor.top() - 21.0;
    let bar_l = frag_cx - bar_w * 0.5;
    let bar_r = frag_cx + bar_w * 0.5;
    painter.rect_filled(
        Rect::from_min_max(
            Pos2::new(bar_l, bar_y - 2.5),
            Pos2::new(bar_r, bar_y + 2.5),
        ),
        0.0,
        red,
    );
    // 两端的空心小方块 + 到供体顶边的交换叉。
    for sx in [bar_l, bar_r] {
        painter.rect_stroke(
            Rect::from_center_size(Pos2::new(sx, bar_y), egui::vec2(7.0, 7.0)),
            0.0,
            Stroke::new(1.0, ink),
            egui::StrokeKind::Inside,
        );
        let (y0, y1) = (bar_y + 4.0, donor.top());
        painter.line_segment(
            [Pos2::new(sx - 5.0, y0), Pos2::new(sx + 5.0, y1)],
            Stroke::new(1.0, ink),
        );
        painter.line_segment(
            [Pos2::new(sx + 5.0, y0), Pos2::new(sx - 5.0, y1)],
            Stroke::new(1.0, ink),
        );
    }
    painter.text(
        Pos2::new(bar_l + 1.0, bar_y - 5.0),
        egui::Align2::LEFT_BOTTOM,
        "attB1 →",
        font.clone(),
        ink,
    );
    painter.text(
        Pos2::new(bar_r - 1.0, bar_y - 5.0),
        egui::Align2::RIGHT_BOTTOM,
        "← attB2",
        font.clone(),
        ink,
    );

    // ---- 左：供体（圆角胶囊），attP 位点在顶边上 ----
    painter.rect_stroke(
        donor,
        23.0,
        Stroke::new(2.6, ink),
        egui::StrokeKind::Inside,
    );
    painter.line_segment(
        [
            Pos2::new(bar_l, donor.top()),
            Pos2::new(bar_r, donor.top()),
        ],
        Stroke::new(3.0, ink),
    );
    painter.text(
        Pos2::new(donor.left() + 22.0, donor.top() + 7.0),
        egui::Align2::LEFT_TOP,
        "attP1 →",
        font.clone(),
        ink,
    );
    painter.text(
        Pos2::new(donor.right() - 22.0, donor.top() + 7.0),
        egui::Align2::RIGHT_TOP,
        "← attP2",
        font.clone(),
        ink,
    );

    // ---- 右：Entry Clone（胶囊 + 顶边上的 attL 红带） ----
    painter.rect_stroke(
        product,
        26.0,
        Stroke::new(2.6, ink),
        egui::StrokeKind::Inside,
    );
    let lw = 40.0;
    let l1 = prod_cx - lw;
    let l2 = prod_cx + lw;
    painter.line_segment(
        [Pos2::new(l1, product.top()), Pos2::new(l2, product.top())],
        Stroke::new(5.0, red),
    );
    painter.text(
        Pos2::new(l1 + 1.0, product.top() - 4.0),
        egui::Align2::LEFT_BOTTOM,
        "attL1 →",
        font.clone(),
        ink,
    );
    painter.text(
        Pos2::new(l2 - 1.0, product.top() - 4.0),
        egui::Align2::RIGHT_BOTTOM,
        "← attL2",
        font.clone(),
        ink,
    );

    // ---- 中间的反应箭头 ----
    let ax0 = donor.right() + 14.0;
    let ax1 = product.left() - 14.0;
    let ay = mid_y + 2.0;
    painter.line_segment(
        [Pos2::new(ax0, ay), Pos2::new(ax1 - 1.0, ay)],
        Stroke::new(1.8, ink),
    );
    let head = 8.0;
    painter.add(egui::Shape::convex_polygon(
        vec![
            Pos2::new(ax1, ay),
            Pos2::new(ax1 - head, ay - head * 0.62),
            Pos2::new(ax1 - head, ay + head * 0.62),
        ],
        ink,
        Stroke::NONE,
    ));
}

/// TA / GC / TOPO 的拼接示意图，照 SnapGene 8.2.3 实拍：
///
/// ```text
///   ▬▬▬▬▬▬                     ╭─ ─ ─ ─ ─╮
///  ┌─────────┐   +    ╭─────────────────╮   →   ╭───────────────╮
///  │Fragment │        │ Linearized Vec. │       │    Product    │
///  │ 650 bp  │  → ←   ╰─────────────────╯       ╰───────────────╯
///  └─────────┘          ..GCCCTT   AGGGC..      ..GCCCTTACAAGT…
///  ACAAGT...AGTGGT      ..CGGGA    TTCCG..      ..CGGAATGTTCA…
///  TGTTCA...TCACCA
/// ```
///
/// - 片段：实心圆角块（有 3′ A 突出端时蓝色，平末端灰色）+ 红色提示条 +
///   两个方向按钮 + 粉色末端盒；
/// - 载体：顶部开口的胶囊（开口 = 线性化切开的位置），下面画两端序列；
/// - 产物：闭合胶囊 + 顶部粉色插入段；两端相容时末端盒合并成一条。
fn paint_ta_topo_schematic(painter: &egui::Painter, rect: Rect, dlg: &CloningDlg) {
    let ink = theme::TEXT;
    let bg = Color32::from_rgb(0xfa, 0xfa, 0xfa);
    let gray_line = Color32::from_rgb(0xa0, 0xa0, 0xa0);
    let dark_line = Color32::from_rgb(0x2b, 0x2b, 0x2b);
    let blue = Color32::from_rgb(0x2f, 0x6f, 0xc4);
    let red_bar = Color32::from_rgb(0xe0, 0x1c, 0x1c);
    let pink_line = Color32::from_rgb(0xff, 0x8f, 0x8f);
    let pink_box = Color32::from_rgb(0xff, 0x9c, 0x9c);
    let pink_faint = Color32::from_rgb(0xff, 0xc6, 0xc6);
    let gray_box = Color32::from_rgb(0xc8, 0xc8, 0xc8);
    let gray_faint = Color32::from_rgb(0xdd, 0xdd, 0xdd);

    let frag_pick = dlg.fragments.first();
    let fsrc = frag_pick.and_then(|p| dlg.pick_source(p));
    let vsrc = dlg.pick_source(&dlg.vector);
    let ready = dlg.preview.is_some() && dlg.primers_ready();
    let fragment = fsrc.zip(frag_pick).and_then(|(s, pick)| end_cloning_fragment(s, pick, dlg.method).ok());
    // 有载体参与时片段条变蓝（SnapGene 用它表示「这个片段已进入反应」）；
    // 还没选载体时是灰色。
    let active = vsrc.is_some();

    // ---- 版面 ----
    let ends_w = 118.0;
    let label_w = 74.0;
    let cap_w = 196.0;
    let prod_w = 244.0;
    let gap1 = 44.0;
    let gap2 = 48.0;
    let total = ends_w + gap1 + cap_w + gap2 + prod_w;
    let x0 = rect.left() + ((rect.width() - total) * 0.5).max(8.0);
    let top = rect.top() + 5.0;
    let cap_top = top + 10.0;
    let cap_h = 44.0;
    let box_top = top + 76.0;
    let box_h = 34.0;

    // ---- 片段 ----
    let frag_x = x0;
    let cx = frag_x + ends_w * 0.5;
    painter.rect_filled(
        Rect::from_min_size(Pos2::new(cx - 46.0, top), Vec2::new(92.0, 4.0)),
        0.0,
        red_bar,
    );
    let label = Rect::from_min_size(
        Pos2::new(cx - label_w * 0.5, top + 10.0),
        Vec2::new(label_w, 30.0),
    );
    painter.rect_filled(
        label,
        4.0,
        if active {
            blue
        } else {
            Color32::from_rgb(0x8a, 0x8a, 0x8a)
        },
    );
    let white = Color32::WHITE;
    bold_text(
        painter,
        Pos2::new(cx, label.top() + 6.0),
        egui::Align2::CENTER_TOP,
        "Fragment",
        10.0,
        white,
    );
    painter.text(
        Pos2::new(cx, label.top() + 18.0),
        egui::Align2::CENTER_TOP,
        format!("{} bp", fragment.as_ref().map(Source::display_len).unwrap_or(0)),
        egui::FontId::proportional(9.0),
        white,
    );

    // 方向按钮（→ 蓝底 = 当前正向）。
    let reversed = frag_pick.is_some_and(|p| p.reversed);
    for (i, back) in [false, true].into_iter().enumerate() {
        let r = Rect::from_min_size(
            Pos2::new(cx - 24.0 + i as f32 * 27.0, top + 48.0),
            Vec2::new(22.0, 17.0),
        );
        let on = back == reversed;
        painter.rect(
            r,
            2.0,
            if on { blue } else { Color32::from_rgb(0xe8, 0xe8, 0xe8) },
            Stroke::new(0.8, Color32::from_rgb(0x8a, 0x8a, 0x8a)),
            egui::StrokeKind::Inside,
        );
        let (a, b) = if back {
            (
                Pos2::new(r.right() - 4.0, r.center().y),
                Pos2::new(r.left() + 4.0, r.center().y),
            )
        } else {
            (
                Pos2::new(r.left() + 4.0, r.center().y),
                Pos2::new(r.right() - 4.0, r.center().y),
            )
        };
        painter.line_segment([a, b], Stroke::new(1.4, if on { white } else { ink }));
        let head = 3.5;
        let tip = b;
        let dir = if back { 1.0 } else { -1.0 };
        painter.add(egui::Shape::convex_polygon(
            vec![
                tip,
                Pos2::new(tip.x - dir * head, tip.y - head),
                Pos2::new(tip.x - dir * head, tip.y + head),
            ],
            if on { white } else { ink },
            Stroke::NONE,
        ));
    }

    let frag_rows = fragment.as_ref().map(ends_both);
    if let Some(rows) = &frag_rows {
        draw_ends_box(
            painter,
            Rect::from_min_size(Pos2::new(frag_x, box_top), Vec2::new(ends_w, box_h)),
            rows,
            pink_box,
            ink,
        );
    }

    // ---- 加号 ----
    painter.text(
        Pos2::new(frag_x + ends_w + gap1 * 0.5, cap_top + cap_h * 0.5),
        egui::Align2::CENTER_CENTER,
        "+",
        egui::FontId::proportional(19.0),
        ink,
    );

    // ---- 线性化载体 ----
    let vx = frag_x + ends_w + gap1;
    let vcap = Rect::from_min_size(Pos2::new(vx, cap_top), Vec2::new(cap_w, cap_h));
    let (vcol, vwidth) = if active {
        (dark_line, 4.2)
    } else {
        (gray_line, 3.4)
    };
    capsule_open(painter, vcap, 0.66, Stroke::new(vwidth, vcol));
    bold_text(
        painter,
        Pos2::new(vcap.center().x, vcap.top() + 9.0),
        egui::Align2::CENTER_TOP,
        "Linearized Vector",
        10.0,
        ink,
    );
    if let Some(s) = vsrc {
        painter.text(
            Pos2::new(vcap.center().x, vcap.top() + 22.0),
            egui::Align2::CENTER_TOP,
            format!("{} bp", s.display_len()),
            egui::FontId::proportional(9.0),
            Color32::from_rgb(0x66, 0x66, 0x66),
        );
        let r = ends_right(s);
        draw_ends_box(
            painter,
            Rect::from_min_size(Pos2::new(vx, box_top), Vec2::new(60.0, box_h)),
            &r,
            gray_box,
            ink,
        );
        let l = ends_left(s);
        draw_ends_box(
            painter,
            Rect::from_min_size(
                Pos2::new(vx + cap_w - 60.0, box_top),
                Vec2::new(60.0, box_h),
            ),
            &l,
            gray_box,
            ink,
        );
    }

    // ---- 箭头 ----
    let ax0 = vx + cap_w + 7.0;
    let ax1 = vx + cap_w + gap2 - 7.0;
    let ay = cap_top + cap_h * 0.5;
    painter.line_segment(
        [Pos2::new(ax0, ay), Pos2::new(ax1 - 1.0, ay)],
        Stroke::new(1.8, ink),
    );
    painter.add(egui::Shape::convex_polygon(
        vec![
            Pos2::new(ax1, ay),
            Pos2::new(ax1 - 9.0, ay - 5.5),
            Pos2::new(ax1 - 9.0, ay + 5.5),
        ],
        ink,
        Stroke::NONE,
    ));

    // ---- 产物 ----
    let px = vx + cap_w + gap2;
    let pcap = Rect::from_min_size(Pos2::new(px, cap_top), Vec2::new(prod_w, cap_h));
    let bar_w = prod_w * 0.34;
    if ready {
        // 就绪：开口胶囊，红色插入段正好补上缺口。
        capsule_open(painter, pcap, bar_w / prod_w, Stroke::new(4.2, dark_line));
        painter.line_segment(
            [
                Pos2::new(pcap.center().x - bar_w * 0.5, pcap.top()),
                Pos2::new(pcap.center().x + bar_w * 0.5, pcap.top()),
            ],
            Stroke::new(4.6, red_bar),
        );
    } else {
        painter.rect(
            pcap,
            cap_h * 0.5,
            bg,
            Stroke::new(3.4, gray_line),
            egui::StrokeKind::Inside,
        );
        // 未就绪：粉色插入段浮在胶囊上方。
        painter.line_segment(
            [
                Pos2::new(pcap.center().x - bar_w * 0.5, pcap.top() - 4.0),
                Pos2::new(pcap.center().x + bar_w * 0.5, pcap.top() - 4.0),
            ],
            Stroke::new(4.0, pink_line),
        );
    }
    bold_text(
        painter,
        Pos2::new(pcap.center().x, pcap.top() + 9.0),
        egui::Align2::CENTER_TOP,
        "Product",
        10.0,
        if ready { ink } else { Color32::from_rgb(0x33, 0x33, 0x33) },
    );
    if let Some(p) = &dlg.preview {
        painter.text(
            Pos2::new(pcap.center().x, pcap.top() + 22.0),
            egui::Align2::CENTER_TOP,
            format!("{} bp", p.display_len()),
            egui::FontId::proportional(9.0),
            Color32::from_rgb(0x66, 0x66, 0x66),
        );
    }

    // 产物末端盒：就绪时三段连成一条（交界处按载体首碱基解析）；否则三块
    // 淡色分开。
    let (Some(v), Some(rows)) = (vsrc, &frag_rows) else {
        return;
    };
    let right = ends_right(v);
    let left = ends_left(v);
    let font = egui::FontId::monospace(9.5);
    if ready {
        let f = fragment.as_ref().expect("fragment rows imply a source");
        let mid = frag_mid(f);
        let j = junction_base(v);
        let (r1, r2) = split_ends(&left);
        let box_rect = Rect::from_min_size(Pos2::new(px, box_top), Vec2::new(prod_w, box_h));
        painter.rect_filled(box_rect, 0.0, gray_box);
        let jr = format!("{}{}", j as char, r1);
        let jr2 = r2;
        // 按真实字宽连接三段，避免固定比例让序列在交界处重叠。
        for (row, parts) in [
            [right.0.clone(), mid.0.clone(), format!("{jr}..")],
            [right.1.clone(), mid.1.clone(), format!("{jr2}..")],
        ].into_iter().enumerate() {
            let y = box_top + row as f32 * 17.0;
            let mut x = px + 3.0;
            for (index, text) in parts.into_iter().enumerate() {
                let galley = painter.layout_no_wrap(text, font.clone(), ink);
                let width = galley.size().x;
                if index == 1 {
                    painter.rect_filled(
                        Rect::from_min_size(Pos2::new(x, y), Vec2::new(width, 17.0)),
                        0.0,
                        pink_box,
                    );
                }
                painter.galley(Pos2::new(x, y + 2.0), galley, ink);
                x += width;
            }
        }
    } else {
        let w = prod_w / 3.0;
        for (i, (rows, fill)) in [
            (right, gray_faint),
            (rows.clone(), pink_faint),
            (left, gray_faint),
        ]
        .into_iter()
        .enumerate()
        {
            draw_ends_box(
                painter,
                Rect::from_min_size(
                    Pos2::new(px + w * i as f32, box_top),
                    Vec2::new(w - 2.0, box_h),
                ),
                &rows,
                fill,
                ink,
            );
        }
    }
}

/// 画一遍文字再叠半像素，模拟 SnapGene 的粗体（egui 默认没有粗体字族）。
fn bold_text(
    painter: &egui::Painter,
    pos: Pos2,
    align: egui::Align2,
    text: impl ToString,
    size: f32,
    color: Color32,
) {
    let t = text.to_string();
    let font = egui::FontId::proportional(size);
    painter.text(pos, align, &t, font.clone(), color);
    painter.text(pos + Vec2::new(0.55, 0.0), align, &t, font, color);
}

/// 末端盒：两行等宽序列。
fn draw_ends_box(
    painter: &egui::Painter,
    rect: Rect,
    rows: &(String, String),
    fill: Color32,
    ink: Color32,
) {
    painter.rect_filled(rect, 0.0, fill);
    let font = egui::FontId::monospace(9.5);
    painter.text(
        Pos2::new(rect.left() + 3.0, rect.top() + 2.0),
        egui::Align2::LEFT_TOP,
        &rows.0,
        font.clone(),
        ink,
    );
    painter.text(
        Pos2::new(rect.left() + 3.0, rect.top() + 17.0),
        egui::Align2::LEFT_TOP,
        &rows.1,
        font,
        ink,
    );
}

/// 顶链首/末碱基（不含粘性末端的标记碱基）。
fn top_ends(src: &Source, n: usize) -> (Vec<u8>, Vec<u8>) {
    let (a, b) = src.top_range();
    let head_end = (a + n).min(b);
    let tail_start = b.saturating_sub(n).max(a);
    (
        src.seq[a..head_end].to_vec(),
        src.seq[tail_start..b].to_vec(),
    )
}

fn complement_str(s: &[u8]) -> String {
    s.iter()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'S' => 'S',
            b'W' => 'W',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            b'N' => 'N',
            other => *other as char,
        })
        .collect()
}

fn as_str(b: &[u8]) -> String {
    String::from_utf8_lossy(b).to_string()
}

/// 片段两端：顶链 6 列、底链按 3′ 突出端少几列，标记碱基补在首/末列
/// （SnapGene 的 `ACAAG...AGTGGTR` / `RTGTTC...CACCA`）。
fn ends_both(src: &Source) -> (String, String) {
    let n = src.seq.len();
    let up = src.ends.upstream_stickiness.min(0).unsigned_abs() as usize;
    let dn = src.ends.downstream_stickiness.min(0).unsigned_abs() as usize;
    let head_n = 6usize.saturating_sub(up).max(2);
    let tail_n = 6usize;
    let (a, b) = src.top_range();
    let head = &src.seq[a..(a + head_n).min(b)];
    let tail = &src.seq[b.saturating_sub(tail_n).max(a)..b];
    let bt = tail_n.saturating_sub(dn).max(1);
    let bottom_end = b.saturating_sub(dn).max(a);
    let bottom_tail = &src.seq[bottom_end.saturating_sub(bt).max(a)..bottom_end];
    let row1 = format!("{}...{}", as_str(head), as_str(tail));
    let mut row2 = String::new();
    if up > 0 {
        row2.push_str(&complement_str(&src.seq[..up.min(n)]));
    }
    row2.push_str(&format!(
        "{}...{}",
        complement_str(head),
        complement_str(bottom_tail)
    ));
    (row1, row2)
}

/// 产物盒中间的粉色段：片段顶链前 5 列 + 省略号 + 末 5 列。
fn frag_mid(src: &Source) -> (String, String) {
    let (a, b) = src.top_range();
    let head = &src.seq[a..(a + 5).min(b)];
    let overhang = src.ends.downstream_stickiness.min(0).unsigned_abs() as usize;
    let paired_end = b.saturating_sub(overhang).max(a);
    let tail = &src.seq[paired_end.saturating_sub(5).max(a)..paired_end];
    let marker_len = src.ends.upstream_stickiness.min(0).unsigned_abs() as usize;
    let mut bottom_head = complement_str(&src.seq[..marker_len.min(src.seq.len())]);
    bottom_head.push_str(&complement_str(head));
    (
        format!("{}...{}", as_str(head), as_str(tail)),
        format!("{}...{}", bottom_head, complement_str(tail)),
    )
}

/// 载体下游端：`..GCCCTT` / `..CGGGAA`。
fn ends_right(src: &Source) -> (String, String) {
    let (_, tail) = top_ends(src, 6);
    let t = as_str(&tail);
    let overhang = src.ends.downstream_stickiness.min(0).unsigned_abs() as usize;
    let paired = &tail[..tail.len().saturating_sub(overhang)];
    (format!("..{t}"), format!("..{}", complement_str(paired)))
}

/// 载体上游端：`AGGGC..` / `TCCCG..`。
fn ends_left(src: &Source) -> (String, String) {
    let (head, _) = top_ends(src, 5);
    let t = as_str(&head);
    let marker_len = src.ends.upstream_stickiness.min(0).unsigned_abs() as usize;
    let mut bottom = complement_str(&src.seq[..marker_len.min(src.seq.len())]);
    bottom.push_str(&complement_str(&head));
    (format!("{t}.."), format!("{bottom}.."))
}

/// 去掉 `..` 的纯序列部分（产物盒右侧要接上交界碱基）。
fn split_ends(rows: &(String, String)) -> (String, String) {
    let f = |s: &str| s.trim_matches('.').to_string();
    (f(&rows.0), f(&rows.1))
}

/// 环化交界处解析出来的碱基：载体顶链第一个碱基。
fn junction_base(src: &Source) -> u8 {
    let (a, b) = src.top_range();
    src.seq.get(a).copied().filter(|_| a < b).unwrap_or(b'A')
}

#[cfg(test)]
mod cloning_settings_tests {
    use super::*;

    #[test]
    fn defaults_match_snapgene_actions_section() {
        let s = crate::settings::CloningSettings::default();
        assert_eq!(s.gibson.primer_overlap(), 40);
        assert_eq!(s.hifi.primer_overlap(), 40);
        assert_eq!(s.gibson.minimum_overlap(), 15);
        assert_eq!(s.in_fusion_overlap, 15);
        assert_eq!(s.in_fusion_multiple_overlap, 20);
        assert_eq!(s.golden_gate_upstream_bases, 3);
        assert_eq!(s.topo_ta_vector, "pCR®2.1-TOPO®");
        for (method, overlap) in [
            (Method::Gibson, 40),
            (Method::NeBuilder, 40),
            (Method::InFusion, 15),
            (Method::InFusionMulti, 20),
        ] {
            let dlg = CloningDlg::new(method, Vec::new(), None);
            assert_eq!(dlg.homology_overlap(), overlap, "{method:?}");
        }
    }

    #[test]
    fn custom_overlap_settings_drive_assembly_and_primer_tails() {
        let mut settings = crate::settings::CloningSettings::default();
        settings.gibson.exact_length = 25;
        settings.gibson.min_length = 12;
        let dlg = CloningDlg::new_with_settings(Method::Gibson, Vec::new(), None, settings);
        assert_eq!(dlg.homology_overlap(), 25);
        assert_eq!(dlg.min_overlap(), 12);
    }

    #[test]
    fn homology_primer_minimum_follows_the_configured_overlap() {
        let vector = Source::new("vector", b"ATGCCGTAGCTAGGCTAACGTACGA".to_vec(), true);
        let insert = Source::new("insert", b"CGATCGTTAACCGGTTAAGCTTCAG".to_vec(), false);
        let docs = || {
            vec![
                DocSnap::new("vector", vector.clone()),
                DocSnap::new("insert", insert.clone()),
            ]
        };
        let mut settings = crate::settings::CloningSettings::default();
        settings.gibson.exact_length = 10;
        settings.gibson.min_length = 10;
        let mut dlg = CloningDlg::new_with_settings(Method::Gibson, docs(), Some(0), settings);
        dlg.vector.region = Some((0, 24));
        dlg.fragments[0].doc = Some(1);
        dlg.fragments[0].use_pcr = true;
        dlg.fragments[0].region = Some((0, 24));
        design_homology_primers(&mut dlg).unwrap();
        // Each junction gets one 10 bp tail on the forward primer; the reverse
        // primer only carries the 20 bp annealing region when both sides are PCR.
        assert_eq!(dlg.vector.primer1.len(), 30);
        assert_eq!(dlg.vector.primer2.len(), 20);

        // The default 40 bp overlap cannot fit into a 24 bp body; the message
        // must name the settings-derived minimum, not a fixed 40 bases.
        let mut strict = CloningDlg::new(Method::Gibson, docs(), Some(0));
        strict.vector.region = Some((0, 24));
        strict.fragments[0].doc = Some(1);
        strict.fragments[0].use_pcr = true;
        strict.fragments[0].region = Some((0, 24));
        let error = design_homology_primers(&mut strict).unwrap_err();
        assert!(error.contains("81 bases") && error.contains("40 bp"), "{error}");
        assert!(!error.contains("40 bases"), "{error}");
    }

    #[test]
    fn golden_gate_upstream_bases_and_default_enzyme_come_from_settings() {
        let mut settings = crate::settings::CloningSettings::default();
        settings.golden_gate_upstream_bases = 5;
        settings.golden_gate_default_enzyme = "SapI".into();
        let docs = vec![
            DocSnap::new("vector", Source::new("vector", b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".to_vec(), true)),
            DocSnap::new("insert", Source::new("insert", b"CCCCAAAATTTTGGGGCCCCAAAATTTTGGGG".to_vec(), false)),
        ];
        let mut dlg = CloningDlg::new_with_settings(Method::GoldenGate, docs, Some(0), settings);
        assert_eq!(dlg.gg_enzyme, dlg.gg_enzymes.iter().position(|n| *n == "SapI").unwrap());
        dlg.fragments[0].doc = Some(1);
        dlg.vector.region = Some((0, 32));
        dlg.fragments[0].region = Some((0, 32));
        dlg.gg_overhangs = vec!["AAT".into(), "GCT".into()];
        golden_gate::design(&mut dlg).unwrap();
        assert!(dlg.vector.primer1.starts_with("GCGCGGCTCTTC"), "{}", dlg.vector.primer1);
    }

    #[test]
    fn golden_gate_settings_choose_the_dialog_enzyme_list() {
        let mut gg = crate::ggsettings::GgSettings::default();
        let names = gg.enzyme_names();
        assert!(names.contains(&"BsaI"), "{names:?}");
        assert!(names.contains(&"BpiI"), "isoschizomer names resolve: {names:?}");
        assert!(!names.contains(&"AarI"), "unchecked standard enzyme: {names:?}");

        for (name, on) in gg.standard.iter_mut() {
            if name == "SapI" || name == "BbsI" {
                *on = false;
            }
        }
        for (name, on) in gg.additional.iter_mut() {
            if name == "FokI" {
                *on = true;
            }
        }
        gg.additional.push(("UnknownTestEnzyme".into(), true));
        let names = gg.enzyme_names();
        assert!(!names.contains(&"SapI"), "{names:?}");
        assert!(!names.contains(&"BbsI"), "{names:?}");
        assert!(names.contains(&"FokI"), "{names:?}");
        assert!(!names.contains(&"UnknownTestEnzyme"), "unknown enzymes are skipped: {names:?}");

        let dlg = CloningDlg::new_with_settings_and_gg(
            Method::GoldenGate,
            Vec::new(),
            None,
            crate::settings::CloningSettings::default(),
            &gg,
        );
        assert_eq!(dlg.gg_enzymes, names);
        assert!(dlg.gg_enzyme < dlg.gg_enzymes.len());
    }

    #[test]
    fn empty_golden_gate_settings_fall_back_to_default_enzymes() {
        let mut gg = crate::ggsettings::GgSettings::default();
        for (_, on) in gg.standard.iter_mut().chain(gg.additional.iter_mut()) {
            *on = false;
        }
        let names = gg.enzyme_names();
        assert_eq!(names, crate::ggsettings::DEFAULTS.to_vec());
        assert!(names.iter().all(|n| hx_enzymes::by_name(n).is_some()));
    }

    #[test]
    fn unchecked_default_enzyme_falls_back_to_first_enabled() {
        let mut gg = crate::ggsettings::GgSettings::default();
        for (name, on) in gg.standard.iter_mut() {
            if name == "BsaI" {
                *on = false;
            }
        }
        let settings = crate::settings::CloningSettings {
            golden_gate_default_enzyme: "BsaI".into(),
            ..Default::default()
        };
        let dlg = CloningDlg::new_with_settings_and_gg(
            Method::GoldenGate,
            Vec::new(),
            None,
            settings,
            &gg,
        );
        assert!(!dlg.gg_enzymes.contains(&"BsaI"));
        assert_eq!(dlg.gg_enzyme, 0);
        assert!(hx_enzymes::by_name(dlg.gg_enzymes[0]).is_some());
    }

    #[test]
    fn topo_vector_setting_matches_the_resource_file_name() {
        assert_eq!(
            normalize_vector_name("pCR®2.1-TOPO®"),
            normalize_vector_name("pCR2.1-TOPO (linearized).dna")
        );
    }
}

#[cfg(test)]
mod reaction_scheme_tests {
    use super::*;

    #[test]
    fn bp_lr_diagram_must_include_entry_clone_and_two_reaction_stages() {
        let dlg = CloningDlg::new(Method::GatewayBpLr(1), Vec::new(), None);
        let scheme = reaction_scheme::model(&dlg);
        assert_eq!(scheme.rows.len(), 2);
        assert_eq!(scheme.rows[0].output.label, "Entry Clone");
        assert_eq!(scheme.rows[1].inputs[0].label, "Entry Clone");
        assert_eq!(scheme.rows[1].inputs[1].label, "Destination Vector");
    }

    #[test]
    fn multiple_assembly_diagram_must_not_drop_extra_inputs() {
        let mut dlg = CloningDlg::new(Method::GibsonMulti, Vec::new(), None);
        dlg.fragments.push(Pick::default());
        let scheme = reaction_scheme::model(&dlg);
        assert!(scheme.rows[0].inputs.iter().any(|p| p.label == "Fragment 3"), "every selected fragment must appear");
    }
}

#[cfg(test)]
mod restriction_tests {
    use super::*;

    fn enzyme(name: &str) -> usize {
        hx_enzymes::ENZYMES.iter().position(|e| e.name == name).unwrap()
    }

    fn restriction_dialog(insert: &str, insert_enzyme: &str) -> CloningDlg {
        let docs = vec![
            DocSnap::new("vector", Source::new("vector", b"AAAAGAATTCCCCCCCCCCCCCCCCCGGATCCGGGG".to_vec(), true)),
            DocSnap::new("insert", Source::new("insert", insert.as_bytes().to_vec(), false)),
        ];
        let mut dlg = CloningDlg::new(Method::Restriction, docs, Some(1));
        dlg.vector.doc = Some(0);
        dlg.fragments[0].doc = Some(1);
        dlg.vector.enzyme = Some(enzyme("EcoRI"));
        dlg.vector.enzyme2 = Some(enzyme("BamHI"));
        dlg.fragments[0].enzyme = Some(enzyme(insert_enzyme));
        dlg.fragments[0].enzyme2 = Some(enzyme("BamHI"));
        let cuts = double_cuts_of(&dlg.docs[0].src, &dlg.vector);
        dlg.vector.cut_pair = Some((cuts[1], cuts[0]));
        dlg
    }

    #[test]
    fn restriction_build_checks_both_joins_and_explicit_cut_selection() {
        let mut dlg = restriction_dialog("TTGAATTCAAACCCGGATCCAA", "EcoRI");
        let p = dlg.build_restriction().unwrap();
        assert!(p.circular);
        assert_eq!(p.seq, b"AAAAGAATTCAAACCCGGATCCGGGG");
        assert_eq!(p.junctions.len(), 2);
        assert!(p.junctions.iter().all(|j| j.overlap == 4));
        dlg.fragments[0].reversed = true;
        assert!(dlg.build_restriction().is_err());
        dlg.fragments[0].reversed = false;
        let cuts = double_cuts_of(&dlg.docs[0].src, &dlg.vector);
        dlg.vector.cut_pair = Some((cuts[0], cuts[1]));
        assert!(dlg.build_restriction().is_err());
    }

    #[test]
    fn restriction_incompatible_ends_require_actual_blunting_not_coordinate_shift() {
        let mut dlg = restriction_dialog("TTCTCGAGAAACCCGGATCCAA", "XhoI");
        assert!(dlg.build_restriction().unwrap_err().contains("not compatible"));
        dlg.vector.blunt = true;
        assert!(dlg.build_restriction().is_err());
        dlg.fragments[0].blunt = true;
        let p = dlg.build_restriction().unwrap();
        assert!(p.circular);
        assert!(p.junctions.iter().all(|j| j.overlap == 0));
    }

    #[test]
    fn linear_ligation_does_not_infer_circular_product_from_circular_source() {
        let docs = vec![
            DocSnap::new("a", Source::new("a", b"GGGGAATT".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
                downstream_stickiness: 4, downstream_modification: "FivePrimePhosphorylated".into(), ..Default::default()
            })),
            DocSnap::new("b", Source::new("b", b"AATTCCCC".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
                upstream_stickiness: 4, ..Default::default()
            })),
        ];
        let mut dlg = CloningDlg::new(Method::Ligation, docs, Some(0));
        dlg.fragments[0].doc = Some(0);
        dlg.fragments[1].doc = Some(1);
        // 旧逻辑会错误地把此处的引物字段直接加到连接片段。
        dlg.fragments[0].primer1 = "AAAA".into();
        let p = dlg.build_ligation().unwrap();
        assert!(!p.circular);
        assert_eq!(p.seq, b"GGGGAATTCCCC");
        dlg.docs[0].src.circular = true;
        assert!(dlg.build_ligation().is_err());
    }

    #[test]
    fn linear_ligation_rejects_unphosphorylated_junctions_like_sg() {
        let docs = vec![DocSnap::new("unmodified", Source::new("unmodified", b"ACGT".to_vec(), false))];
        let mut dlg = CloningDlg::new(Method::Ligation, docs, Some(0));
        for pick in &mut dlg.fragments { pick.doc = Some(0); }
        dlg.refresh();
        assert!(dlg.preview.is_none());
        assert!(dlg.error.as_deref().unwrap().contains("at least one DNA end must be 5' phosphorylated"));
        dlg.docs[0].src.ends.upstream_modification = "FivePrimePhosphorylated".into();
        // 文档快照不可变；测试显式失效对应新来源快照。
        dlg.preview_key.clear(); dlg.refresh();
        assert!(dlg.preview.is_some());
    }

    #[test]
    fn restriction_can_use_a_prelinearized_source_without_selecting_an_enzyme() {
        let s = Source::new("linear", b"AATTGGGG".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
            upstream_stickiness: 4, downstream_modification: "FivePrimePhosphorylated".into(), ..Default::default()
        });
        let f = restriction_fragment(&s, &Pick::default(), true, true).unwrap();
        assert_eq!(f.seq, s.seq);
        assert_eq!(f.ends, s.ends);
    }

    #[test]
    fn single_linear_cut_produces_a_terminal_fragment_not_an_empty_circle() {
        let src = Source::new("linear", b"AAAAGAATTCTTTT".to_vec(), false);
        let pick = Pick { enzyme: Some(enzyme("EcoRI")), ..Default::default() };
        let fragment = restriction_fragment(&src, &pick, false, true).unwrap();
        assert!(!fragment.seq.is_empty());
        assert!(!fragment.circular);
        let bad = Pick { enzyme: Some(enzyme("BsaI")), ..Default::default() };
        assert!(restriction_fragment(&src, &bad, false, true).is_err());
    }

    #[test]
    fn selecting_different_cut_boundaries_invalidates_product_preview() {
        let mut dlg = restriction_dialog("TTGAATTCAAACCCGGATCCAA", "EcoRI");
        dlg.refresh();
        assert!(dlg.preview.is_some());
        let before = dlg.state_key();
        let (a, b) = dlg.vector.cut_pair.unwrap();
        dlg.vector.cut_pair = Some((b, a));
        assert_ne!(dlg.state_key(), before);
        dlg.refresh();
        assert!(dlg.preview.is_none());
        assert!(dlg.error.as_deref().unwrap().contains("not compatible"));
    }

    #[test]
    fn ligation_includes_every_fragment_exposed_by_the_fragment_count_control() {
        let docs = vec![DocSnap::new("blunt", Source::new("blunt", b"ACGT".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
            upstream_modification: "FivePrimePhosphorylated".into(), ..Default::default()
        }))];
        let mut dlg = CloningDlg::new(Method::Ligation, docs, Some(0));
        for pick in &mut dlg.fragments { pick.doc = Some(0); }
        dlg.fragments.push(Pick { doc: Some(0), ..Default::default() });
        assert_eq!(dlg.build_ligation().unwrap().seq, b"ACGTACGTACGT");
        assert_eq!(product_sections(&dlg).len(), 3);
    }

    #[test]
    fn linear_product_preview_counts_top_strand_but_saved_file_retains_outer_ends() {
        let docs = vec![
            DocSnap::new("left", Source::new("left", b"TGGGGAATT".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
                upstream_stickiness: -1, downstream_stickiness: 4, downstream_modification: "FivePrimePhosphorylated".into(), ..Default::default()
            })),
            DocSnap::new("right", Source::new("right", b"AATTCCCCAA".to_vec(), false).with_ends(hx_formats::dna::EndProperties {
                upstream_stickiness: 4, downstream_stickiness: 2, ..Default::default()
            })),
        ];
        let mut dlg = CloningDlg::new(Method::Ligation, docs, Some(0));
        dlg.fragments[0].doc = Some(0);
        dlg.fragments[1].doc = Some(1);
        dlg.tab = 1;
        dlg.refresh();
        let p = dlg.preview.as_ref().unwrap();
        assert_eq!(p.display_len(), 12);
        assert_eq!(dlg.build_view_doc().unwrap().sequence(), b"GGGGAATTCCCC");
        let file = hx_formats::DnaFile::parse(&cloning::product_bytes(p, "ligation")).unwrap();
        assert_eq!(file.sequence().unwrap().sequence, b"TGGGGAATTCCCCAA");
        assert_eq!(file.end_properties(), p.ends);
    }
}

#[cfg(test)]
mod ta_topo_display_tests {
    use super::*;

    #[test]
    fn sticky_end_rows_exclude_markers_from_top_and_unpaired_bottom() {
        let mut vector = Source::new("vector", b"AAGGGCGCCCTT".to_vec(), false);
        vector.ends.upstream_stickiness = -1;
        vector.ends.downstream_stickiness = -1;
        assert_eq!(ends_left(&vector), ("AGGGC..".into(), "TTCCCG..".into()));
        assert_eq!(ends_right(&vector), ("..GCCCTT".into(), "..CGGGA".into()));

        let insert = Source::new("insert", b"ACAAGTTTAGTGGT".to_vec(), false);
        let tailed = cloning::add_a_overhangs(&insert);
        assert_eq!(ends_both(&tailed), ("ACAAG...GTGGTA".into(), "ATGTTC...CACCA".into()));
        assert_eq!(frag_mid(&tailed), ("ACAAG...GTGGT".into(), "ATGTTC...CACCA".into()));
        assert_eq!(tailed.display_len(), insert.display_len() + 1);
    }

    #[test]
    fn source_preview_uses_top_strand_length() {
        let mut source = Source::new("vector", b"AAGGGCGCCCTT".to_vec(), false);
        source.ends.upstream_stickiness = -1;
        source.ends.downstream_stickiness = -1;
        source.ends.original_first_base = Some(7);
        let bytes = source_bytes(&source, "vector");
        let file = hx_formats::DnaFile::parse(&bytes).unwrap();
        assert_eq!(file.sequence().unwrap().sequence, source.seq);
        let doc = crate::doc::Document::from_bytes(bytes).unwrap();
        assert_eq!(doc.sequence(), source.top_seq());
        assert_eq!(doc.end_properties(), source.ends);
    }

    #[test]
    fn installed_topo_vector_keeps_3930_bp_display_and_3931_native_bases() {
        let path = Path::new("/opt/gslbiotech/snapgene/resources/Plasmids/TOPO TA Vectors/pCR2.1-TOPO (linearized).dna");
        if !path.is_file() {
            eprintln!("Skipping installed-vector regression: SnapGene resources not present");
            return;
        }
        let bytes = std::fs::read(path).unwrap();
        let doc = crate::doc::Document::from_bytes(bytes.clone()).unwrap();
        assert_eq!(doc.len(), 3930);
        assert_eq!(doc.session.stored_sequence().len(), 3931);
        assert_eq!(doc.session.save_bytes().unwrap(), bytes);
        let source = load_std_source(path).unwrap();
        assert_eq!(source.seq.len(), 3931);
        assert_eq!(source.display_len(), 3930);
        assert_eq!(source.ends.original_first_base, Some(3637));
        let preview = crate::doc::Document::from_bytes(source_bytes(&source, "vector")).unwrap();
        assert_eq!(preview.sequence(), doc.sequence());
        assert_eq!(preview.end_properties(), source.ends);
        let snap = doc_from_path(path).unwrap();
        assert_eq!(snap.src.seq, source.seq);
        assert_eq!(snap.src.features.len(), source.features.len());
        assert_eq!(snap.src.ends, source.ends);
    }
}

/// 顶部开口的胶囊：从左开口沿左侧绕到底部再到右开口（= 线性化切开的环）。
/// `gap` 是开口占宽度的比例。
fn capsule_open(painter: &egui::Painter, rect: Rect, gap: f32, stroke: Stroke) {
    let r = rect.height() * 0.5;
    let g = (rect.width() * gap.clamp(0.05, 0.9)) * 0.5;
    let gl = rect.center().x - g;
    let gr = rect.center().x + g;
    let mut pts: Vec<Pos2> = Vec::new();
    pts.push(Pos2::new(gl, rect.top()));
    pts.push(Pos2::new(rect.left() + r, rect.top()));
    // 左上角：-90° → -180°
    for i in 0..=8 {
        let t = -std::f32::consts::FRAC_PI_2 - std::f32::consts::FRAC_PI_2 * i as f32 / 8.0;
        pts.push(Pos2::new(
            rect.left() + r + r * t.cos(),
            rect.top() + r + r * t.sin(),
        ));
    }
    pts.push(Pos2::new(rect.left(), rect.bottom() - r));
    // 左下角：180° → 90°
    for i in 0..=8 {
        let t = std::f32::consts::PI - std::f32::consts::FRAC_PI_2 * i as f32 / 8.0;
        pts.push(Pos2::new(
            rect.left() + r + r * t.cos(),
            rect.bottom() - r + r * t.sin(),
        ));
    }
    pts.push(Pos2::new(rect.right() - r, rect.bottom()));
    // 右下角：90° → 0°
    for i in 0..=8 {
        let t = std::f32::consts::FRAC_PI_2 - std::f32::consts::FRAC_PI_2 * i as f32 / 8.0;
        pts.push(Pos2::new(
            rect.right() - r + r * t.cos(),
            rect.bottom() - r + r * t.sin(),
        ));
    }
    pts.push(Pos2::new(rect.right(), rect.top() + r));
    // 右上角：0° → -90°
    for i in 0..=8 {
        let t = -std::f32::consts::FRAC_PI_2 * i as f32 / 8.0;
        pts.push(Pos2::new(
            rect.right() - r + r * t.cos(),
            rect.top() + r + r * t.sin(),
        ));
    }
    pts.push(Pos2::new(gr, rect.top()));
    painter.add(egui::Shape::line(pts, stroke));
}

// ---------------------------------------------------------------------------
// 小工具
// ---------------------------------------------------------------------------

fn parse_color(s: &str) -> Option<Color32> {
    let h = s.strip_prefix('#')?;
    if h.len() != 6 {
        return None;
    }
    let r = u8::from_str_radix(&h[0..2], 16).ok()?;
    let g = u8::from_str_radix(&h[2..4], 16).ok()?;
    let b = u8::from_str_radix(&h[4..6], 16).ok()?;
    Some(Color32::from_rgb(r, g, b))
}
