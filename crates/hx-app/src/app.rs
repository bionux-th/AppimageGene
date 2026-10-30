//! 应用外壳：对齐 SnapGene 的浅色界面 —— 菜单栏 / 图标工具栏 / 文档页签 / 视图页签 /
//! 左侧文件面板 + 工具条 / 画布 / 状态栏。

use std::path::PathBuf;

use eframe::egui;

use hx_algo::primer;
use hx_formats::feature_edit::FeatureEdit;
use hx_formats::primer_edit::PrimerEdit;
use hx_formats::write::Edit;

use crate::doc::{indices_of, Document, SiteHits, COMMON_ENZYMES};
use crate::{icons, theme};

/// 主区视图页签。
///
/// 页签集合随分子类型变化（oracle 实测）：
/// - DNA：`Map | Sequence | Enzymes | Features | Primers | History`
/// - RNA：`Map | Sequence | Features | Secondary Structure | History`
/// - 蛋白：`Map | Sequence | Properties | Features | History`
/// 文档画布的宿主：主窗口或克隆对话框（独立的原生窗口）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanvasMode {
    Main,
    Cloning,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Map,
    Sequence,
    Enzymes,
    Features,
    Primers,
    History,
    /// RNA 二级结构（SnapGene 用 ViennaRNA 计算；我们只显示已存数据的说明）。
    SecondaryStructure,
    /// 蛋白属性（氨基酸组成表）。
    Properties,
    /// 不能编辑的原生 SnapGene 文档（比对/凝胶）的**概览**页。
    NativeOverview,
}

impl Tab {
    /// DNA 的页签集合（历史名 `ALL`；DNA 是默认分子类型）。
    pub const ALL: [Tab; 6] = [
        Tab::Map,
        Tab::Sequence,
        Tab::Enzymes,
        Tab::Features,
        Tab::Primers,
        Tab::History,
    ];
    pub const RNA: [Tab; 5] = [
        Tab::Map,
        Tab::Sequence,
        Tab::Features,
        Tab::SecondaryStructure,
        Tab::History,
    ];
    pub const PROTEIN: [Tab; 5] = [
        Tab::Map,
        Tab::Sequence,
        Tab::Properties,
        Tab::Features,
        Tab::History,
    ];
    /// 不能编辑的原生文档（比对/凝胶）：只有概览页。
    pub const NATIVE: [Tab; 1] = [Tab::NativeOverview];

    /// 该分子类型可用的视图页签。
    pub fn list(kind: hx_core::MoleculeKind) -> &'static [Tab] {
        match kind {
            hx_core::MoleculeKind::Dna => &Self::ALL,
            hx_core::MoleculeKind::Rna => &Self::RNA,
            hx_core::MoleculeKind::Protein => &Self::PROTEIN,
        }
    }

    /// 按文档决定页签集合（原生比对/凝胶文档只有概览页）。
    pub fn list_for(doc: Option<&crate::doc::Document>) -> &'static [Tab] {
        match doc {
            Some(d) if d.native_class().is_some() => &Self::NATIVE,
            Some(d) => Self::list(d.kind()),
            None => &Self::ALL,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Tab::Map => "Map",
            Tab::Sequence => "Sequence",
            Tab::Enzymes => "Enzymes",
            Tab::Features => "Features",
            Tab::Primers => "Primers",
            Tab::History => "History",
            Tab::SecondaryStructure => "Secondary Structure",
            Tab::Properties => "Properties",
            Tab::NativeOverview => "Overview",
        }
    }

    /// 下一个视图页签（Window → Show Next Tab）。
    pub fn next(self) -> Tab {
        let all = Tab::ALL;
        let i = all.iter().position(|t| *t == self).unwrap_or(0);
        all[(i + 1) % all.len()]
    }

    /// 上一个视图页签（Window → Show Previous Tab）。
    pub fn prev(self) -> Tab {
        let all = Tab::ALL;
        let i = all.iter().position(|t| *t == self).unwrap_or(0);
        all[(i + all.len() - 1) % all.len()]
    }
}

/// 「打开」小窗口的类型。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpenKind {
    File,
    Folder,
}

/// 表格视图状态：SnapGene 的列头右键菜单（`ColumnHeaderContextMenu`：排序 / 隐藏列）。
#[derive(Debug, Clone, Default)]
pub struct TableState {
    /// `(列下标, 升序)`；`None` = 该表的默认顺序。
    pub sort: Option<(usize, bool)>,
    /// 被隐藏的列下标。
    pub hidden: std::collections::HashSet<usize>,
}

impl TableState {
    pub fn visible(&self, col: usize) -> bool {
        !self.hidden.contains(&col)
    }
}

/// 酶显示预设（状态栏右侧显示当前集合，对应 SnapGene 的 "Unique 6+ Cutters | Nonredundant"）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnzymePreset {
    Unique6Nonredundant,
    Unique6,
    Nonredundant,
    SixPlus,
    All6,
    Common,
    None,
    Custom,
}

impl EnzymePreset {
    /// 菜单顺序。
    pub const ALL: [EnzymePreset; 7] = [
        EnzymePreset::Unique6Nonredundant,
        EnzymePreset::Unique6,
        EnzymePreset::Nonredundant,
        EnzymePreset::SixPlus,
        EnzymePreset::All6,
        EnzymePreset::Common,
        EnzymePreset::None,
    ];

    pub fn label(self) -> &'static str {
        match self {
            EnzymePreset::Unique6Nonredundant => "Unique 6+ Cutters | Nonredundant",
            EnzymePreset::Unique6 => "Unique 6+ Cutters",
            EnzymePreset::Nonredundant => "Nonredundant",
            EnzymePreset::SixPlus => "6+ Cutters",
            EnzymePreset::All6 => "All 6-Cutters",
            EnzymePreset::Common => "Common Enzymes",
            EnzymePreset::None => "No Enzymes",
            EnzymePreset::Custom => "Custom",
        }
    }
}

/// 单个文档标签页的**视图状态**。
///
/// SnapGene 每个文档标签页都有自己的视图状态（当前视图页签、缩放、选段、显示集合…）；
/// 我们同样按标签页保存：`HelixApp` 上的同名字段始终表示**活动标签页**的值，
/// 切换标签页时用 [`DocViewState::take`] / [`DocViewState::apply`] 存储与恢复。
/// （这样既有多标签页语义，又不必把上百处 `self.zoom` 之类的访问改写成 `self.cur().zoom`。）
pub struct DocViewState {
    pub enzyme_parent: hx_enzymes::sets::Parent,
    pub enzyme_filter: hx_enzymes::sets::Filter,
    pub translation_frames: u8,
    pub orf_min_bp: usize,
    pub sequence_color_mode: u8,
    /// 该文档当前停留的视图页签（Map / Sequence / …）。
    pub tab: Tab,
    pub zoom: f32,
    pub pan_bp: f64,
    pub rotation: f32,
    pub seq_zoom: f32,
    pub seq_pan_x: f32,
    pub seq_scroll_y: f32,
    pub selection: Option<(usize, usize)>,
    pub cursor_bp: Option<usize>,
    pub hover_bp: Option<usize>,
    pub selected_feature: Option<usize>,
    pub selected_residue: Option<(usize, usize)>,
    pub sel_anchor: Option<usize>,
    /// 蛋白 `Properties` 页：展开的明细行与选中的 pH（按文档保存）。
    pub props_expand: [bool; 3],
    pub props_ph: f64,
    pub hidden_features: std::collections::HashSet<usize>,
    pub hidden_enzymes: std::collections::HashSet<usize>,
    pub highlighted_enzymes: std::collections::HashSet<usize>,
    pub selected_enzymes: std::collections::HashSet<usize>,
    pub selected_features: std::collections::HashSet<usize>,
    pub selected_primer: Option<usize>,
    pub sites: SiteHits,
    pub sites_dirty: bool,
    pub enzyme_preset: EnzymePreset,
    pub enzyme_search: String,
    pub needs_fit: bool,
}

impl Default for DocViewState {
    fn default() -> Self {
        Self {
            enzyme_parent: hx_enzymes::sets::Parent::Nonredundant,
            enzyme_filter: hx_enzymes::sets::Filter::UniqueSixPlus,
            translation_frames: 63,
            orf_min_bp: hx_algo::orf::DEFAULT_MIN_BP,
            sequence_color_mode: 0,
            tab: Tab::Map,
            zoom: 1.0,
            pan_bp: 0.0,
            rotation: 0.0,
            seq_zoom: 1.0,
            seq_pan_x: 0.0,
            seq_scroll_y: 0.0,
            selection: None,
            cursor_bp: None,
            hover_bp: None,
            selected_feature: None,
            selected_residue: None,
            sel_anchor: None,
            props_expand: [false; 3],
            props_ph: 7.0,
            hidden_features: Default::default(),
            hidden_enzymes: Default::default(),
            highlighted_enzymes: Default::default(),
            selected_enzymes: Default::default(),
            selected_features: Default::default(),
            selected_primer: None,
            sites: SiteHits::default(),
            sites_dirty: true,
            enzyme_preset: EnzymePreset::Unique6Nonredundant,
            enzyme_search: String::new(),
            needs_fit: true,
        }
    }
}

impl DocViewState {
    /// 把活动标签页的视图状态从 `app` 上摘下来（字段被置为默认值）。
    fn take(app: &mut HelixApp) -> Self {
        Self {
            enzyme_parent: app.enzyme_parent,
            enzyme_filter: app.enzyme_filter,
            translation_frames: app.translation_frames,
            orf_min_bp: app.orf_min_bp,
            sequence_color_mode: app.sequence_color_mode,
            tab: app.tab,
            zoom: app.zoom,
            pan_bp: app.pan_bp,
            rotation: app.rotation,
            seq_zoom: app.seq_zoom,
            seq_pan_x: app.seq_pan_x,
            seq_scroll_y: app.seq.scroll_y,
            selection: app.selection.take(),
            cursor_bp: app.cursor_bp.take(),
            hover_bp: app.hover_bp.take(),
            selected_feature: app.selected_feature.take(),
            selected_residue: app.selected_residue.take(),
            sel_anchor: app.sel_anchor.take(),
            props_expand: app.props_expand,
            props_ph: app.props_ph,
            hidden_features: std::mem::take(&mut app.hidden_features),
            hidden_enzymes: std::mem::take(&mut app.hidden_enzymes),
            highlighted_enzymes: std::mem::take(&mut app.highlighted_enzymes),
            selected_enzymes: std::mem::take(&mut app.selected_enzymes),
            selected_features: std::mem::take(&mut app.selected_features),
            selected_primer: app.selected_primer.take(),
            sites: std::mem::take(&mut app.sites),
            sites_dirty: app.sites_dirty,
            enzyme_preset: app.enzyme_preset,
            enzyme_search: std::mem::take(&mut app.enzyme_search),
            needs_fit: app.needs_fit,
        }
    }

    /// 把该状态铺回 `app` 的字段上（成为活动标签页的视图状态）。
    fn apply(self, app: &mut HelixApp) {
        app.enzyme_parent = self.enzyme_parent;
        app.enzyme_filter = self.enzyme_filter;
        app.translation_frames = self.translation_frames;
        app.orf_min_bp = self.orf_min_bp;
        app.sequence_color_mode = self.sequence_color_mode;
        app.tab = self.tab;
        app.zoom = self.zoom;
        app.pan_bp = self.pan_bp;
        app.rotation = self.rotation;
        app.seq_zoom = self.seq_zoom;
        app.seq_pan_x = self.seq_pan_x;
        app.seq.scroll_y = self.seq_scroll_y;
        app.selection = self.selection;
        app.cursor_bp = self.cursor_bp;
        app.hover_bp = self.hover_bp;
        app.selected_feature = self.selected_feature;
        app.selected_residue = self.selected_residue;
        app.sel_anchor = self.sel_anchor;
        app.props_expand = self.props_expand;
        app.props_ph = self.props_ph;
        app.hidden_features = self.hidden_features;
        app.hidden_enzymes = self.hidden_enzymes;
        app.highlighted_enzymes = self.highlighted_enzymes;
        app.selected_enzymes = self.selected_enzymes;
        app.selected_features = self.selected_features;
        app.selected_primer = self.selected_primer;
        app.sites = self.sites;
        app.sites_dirty = self.sites_dirty;
        app.enzyme_preset = self.enzyme_preset;
        app.enzyme_search = self.enzyme_search;
        app.needs_fit = self.needs_fit;
    }
}

/// 一个打开的文档标签页。
///
/// 不变式：**活动标签页的文档在 `HelixApp::doc`**（所以这里 `doc` 为 `None`），
/// 其余标签页的文档保存在自己的条目里。`tabs` 为空 ⟺ `HelixApp::doc` 为空。
pub struct TabDoc {
    pub doc: Option<Document>,
    /// 只读文本 / PDF 阅读页签（与分子文档互斥）。
    pub viewer: Option<crate::viewer::ViewerDoc>,
    pub view: DocViewState,
    /// 另存为/新建时用来显示名字（有路径时以路径为准）。
    pub untitled_name: Option<String>,
}

impl TabDoc {
    /// 标签页显示名（活动标签页由调用方传入其文档/阅读器）。
    pub fn title(
        &self,
        active_doc: Option<&Document>,
        active_viewer: Option<&crate::viewer::ViewerDoc>,
    ) -> String {
        // `untitled_name` 优先：New File 面板的预览页显示 `New Sequence`（SnapGene 实测）
        if let Some(n) = &self.untitled_name {
            return n.clone();
        }
        if let Some(v) = active_viewer.or(self.viewer.as_ref()) {
            let t = v.title();
            return if v.is_dirty() { format!("{t} *") } else { t };
        }
        let doc = active_doc.or(self.doc.as_ref());
        match doc {
            Some(d) => {
                let n = d.file_label();
                if d.is_dirty() {
                    format!("{n} *")
                } else {
                    n
                }
            }
            None => self
                .untitled_name
                .clone()
                .unwrap_or_else(|| "Untitled".to_string()),
        }
    }
}

pub struct HelixApp {
    pub doc: Option<Document>,
    /// 活动的只读文本 / PDF 阅读页签（与 `doc` 互斥）。
    pub viewer: Option<crate::viewer::ViewerDoc>,
    /// 打开的文档标签页（SnapGene 的选项卡栏）。
    pub tabs: Vec<TabDoc>,
    /// 活动标签页下标（`tabs` 为空时无意义）。
    pub active: usize,
    pub status: String,
    /// 侧栏"打开"路径输入。
    pub open_path: String,
    /// 侧栏"另存为"路径输入。
    pub save_path: String,
    /// 本次会话打开过的文件（侧栏 Other Files）。
    pub recent: Vec<PathBuf>,
    /// 主区页签。
    pub tab: Tab,
    /// 环状图谱的显示原点（弧度）。用户拖拽不会改变它；“居中特征”等
    /// 明确的定位操作可以程序化地改变原点。
    pub rotation: f32,
    /// 图谱缩放（环状：半径倍数；线性：px/bp；序列视图用 `seq_zoom`）。
    pub zoom: f32,
    /// 线性视图最左可见碱基。
    pub pan_bp: f64,
    /// 酶位点开关（SnapGene 的 Display → Enzyme Sites）。
    pub show_sites: bool,
    /// 刻度开关（图谱）。
    pub show_ruler: bool,
    /// 图谱显示布局（Auto/Circular/Linear；SnapGene 的 Display → 环/直线）。
    pub map_layout: hx_render::MapLayout,
    /// 顶部工具栏开关（View → Top Toolbar）。
    pub show_toolbar: bool,
    /// 左侧「文件」侧栏是否展开（SnapGene 侧栏的 « / » 折叠）。
    pub sidebar_open: bool,
    pub sidebar_width:f32,
    /// 侧栏文件过滤框（SnapGene 的「搜索和过滤文件」）。
    pub file_filter: String,
    /// 侧栏**文件夹树**过滤框（SnapGene 打开文件夹后的 `Search and filter files`）。
    pub tree_filter: String,
    /// 侧栏文件树当前文件夹（持久化）。
    pub folder: Option<std::path::PathBuf>,
    /// 「Open Folder…」的路径输入框。
    pub folder_input: String,
    /// 已展开的子文件夹。
    pub expanded: std::collections::HashSet<std::path::PathBuf>,
    /// 目录列表缓存（key = 目录；打开文件/刷新时清空）。
    pub dir_cache: std::collections::HashMap<std::path::PathBuf, Vec<crate::settings::Entry>>,
    /// 目录读取错误提示。
    pub folder_error: Option<String>,
    /// 侧栏是否显示文件扩展名。
    /// 左侧工具条开关（View → Side Toolbar）。
    pub show_side_tools: bool,
    /// 当前画布是否在克隆对话框里：SnapGene 的对话框不画线性图谱底部的
    /// 平移/缩放条，分子名也放在轴下方而不是窗口底部。
    pub cloning_canvas: bool,
    /// 线性图谱两端 `Start`/`End` 标注开关（View → Show "Start" and "End"）。
    pub show_end_labels: bool,
    /// 「View Options…」窗口开关（等价于 Display 下拉）。
    pub show_view_options: bool,
    /// 「Enter Sequence Text…」窗口（SnapGene 的 Paste 直接插入；我们额外提供文本框）。
    pub show_paste_dialog: bool,
    pub molecular_action: Option<crate::molecular_actions::Dialog>,
    pub pending_delete: Option<(usize, usize)>,
    pub insert_unknown: bool,
    pub insert_unknown_count: usize,
    pub insert_reverse: bool,
    /// 「New DNA/RNA File / New Protein File」对话框（SnapGene 8 的 New File 面板）。
    pub new_file: Option<crate::newfile::NewFileState>,
    /// 「Add Feature… / Edit Feature…」对话框（SnapGene 8 的特征面板）。
    pub feature_dlg: Option<crate::featuredlg::FeatureDlg>,
    /// 「Add Primer… / Edit Primer…」对话框（`Primers ▸ Add Primer…`）。
    pub primer_dlg: Option<crate::primerdlg::PrimerDlg>,
    /// 有选区时先问 `Make a primer from the selection?`（0-based 半开区间）。
    pub primer_choice: Option<(usize, usize)>,
    /// Primers 面板的搜索框文案（`Find primers`）。
    pub primer_search: String,
    /// Primers 面板里展开详情的行下标。
    pub primer_expanded: std::collections::BTreeSet<usize>,
    /// Primers 面板里勾选的引物下标（批量操作）。
    pub primer_checked: std::collections::BTreeSet<usize>,
    /// 「Design Primers…」对话框（SnapGene 的 `Actions ▸ PCR…`）。
    pub design_dlg: Option<crate::primerdlg::DesignDlg>,
    /// 分子克隆对话框（`Actions ▸ … Cloning ▸ Insert Fragment…`）。
    pub cloning: Option<crate::clonedlg::CloningDlg>,
    /// 「Golden Gate Settings…」窗口状态（跨启动持久化）。
    pub gg_settings: crate::ggsettings::GgSettings,
    /// Golden Gate Settings 窗口是否打开。
    pub gg_settings_open: bool,
    /// 打开设置窗口前的快照；Cancel / × 时恢复。
    pub gg_settings_backup: Option<crate::ggsettings::GgSettings>,
    /// 克隆默认参数（SnapGene `[actions]`，跨启动持久化）。
    pub cloning_settings: crate::settings::CloningSettings,
    /// 「Analyze Selected Primers…」结果窗口（`(标签, 值)` 行）。
    pub primer_analysis: Option<Vec<(String, String)>>,
    /// 「Feature Color…」调色板（作用对象 = `selected_feature`）。
    pub feature_color_popup: bool,
    /// 打开文件/文件夹的小窗口（无系统对话框；SnapGene 用系统对话框）。
    pub show_open_dialog: Option<OpenKind>,
    /// 被隐藏的特征下标（SnapGene 的 Features 视图「Toggle visibility」）。
    pub hidden_features: std::collections::HashSet<usize>,
    /// 被隐藏的酶（Enzymes 视图最左的勾选框；索引 = `sites.selected` 下标）。
    pub hidden_enzymes: std::collections::HashSet<usize>,
    /// 被高亮的酶位点（Enzymes ▸ Actions ▸ Highlight Selected Enzyme Sites）。
    pub highlighted_enzymes: std::collections::HashSet<usize>,
    /// Enzymes 视图里被选中的行（Actions 菜单的作用对象）。
    pub selected_enzymes: std::collections::HashSet<usize>,
    /// Features 视图里被选中的行（Actions 菜单的作用对象）。
    pub selected_features: std::collections::HashSet<usize>,
    /// Primers 视图里被选中的行（`Edit Primer…` 等菜单项的作用对象）。
    pub selected_primer: Option<usize>,
    /// 酶位点显示方式：`true` = 线（每行一条全长轨道 + 刻度），`false` = 数字（SnapGene 的
    /// 「Show Sites as Lines / Numbers」）。
    pub sites_as_lines: bool,
    /// 序列视图行宽：`None` = 换行到窗口（自适应），`Some(n)` = 固定每行 n 个碱基
    /// （SnapGene 的「换行到窗口 / 使用固定线宽」）。
    pub seq_line_bp: Option<usize>,
    /// 固定线宽下序列内容横向平移（像素）。
    pub seq_pan_x: f32,
    /// Enzymes 视图里已展开的酶（行内明细）。
    pub expanded_enzymes: std::collections::HashSet<usize>,
    /// Features 视图里已展开的特征（行内明细）。
    pub expanded_features: std::collections::HashSet<usize>,
    /// Features 面板的查找框。
    pub feature_search: String,
    /// 表状态（排序 / 隐藏列）：Features / Enzymes / Primers。
    pub tbl_features: TableState,
    pub tbl_enzymes: TableState,
    pub tbl_primers: TableState,
    /// 统一的注释颜色开关（Map/Sequence 的特征、引物、翻译和 ORF 颜色）。
    pub show_colors: bool,
    /// 序列视图：显示互补链（SnapGene 默认显示）。
    pub show_complement: bool,
    /// Map/Sequence：显示额外的翻译 / ORF；特征内置翻译由 show_features 控制。
    pub show_translation: bool,
    pub translation_frames: u8,
    pub orf_min_bp: usize,
    pub hybridization_draft: Option<hx_formats::annot::HybridizationParams>,
    pub show_translation_options: bool,
    pub enzyme_preferences_open: bool,
    pub preferences: crate::preferences::Preferences,
    pub rename_file: Option<crate::rename::Rename>,
    enzyme_manager: crate::enzyme_manager::Manager,
    pub(crate) gel_view: Option<Box<crate::gel_view::GelView>>,
    restoring_workspace: bool,
    persistence_enabled: bool,
    pub sequence_color_mode: u8,
    /// 氨基酸用单字母缩写（SnapGene 的 `Use 1-Letter Amino Acid Codes`；默认三字母）。
    pub aa_one_letter: bool,
    /// Map/Sequence：统一显示特征箭头、箭头内名称、外置标签和概览。
    pub show_features: bool,
    /// Map/Sequence：显示引物结合位点。
    pub show_primers: bool,
    /// RNA `Secondary Structure` 页的缩放（`1.0` = 100%；`<= 0` = 还没适配过，按 Zoom to fit）。
    pub struct_zoom: f32,
    /// 结构视图的平移（屏幕像素）。
    pub struct_pan: egui::Vec2,
    /// 结构视图当前显示的条目（`0` = MFE，其余为次优结构）。
    pub struct_sel: usize,
    /// 结构视图里点中的碱基（0-based）；它的配对碱基一起高亮。
    pub struct_sel_base: Option<usize>,
    /// 结构视图：显示每 10 个碱基的序号标签。
    pub struct_show_index: bool,
    /// 上一次算出的 Zoom to fit 比例（组合框显示用）。
    pub struct_fit: f32,
    /// **蛋白**左侧工具条的三段特征开关（SnapGene 的 `Show "Region" features` /
    /// `Show "Site" and "Het" features` / `Show "Bond" features`；核苷酸用 `show_features`）。
    pub show_regions: bool,
    pub show_sites_het: bool,
    pub show_bonds: bool,
    /// 蛋白图谱条下的特征起止数字（SnapGene 的 `Show translation numbering`）。
    pub show_trans_numbering: bool,
    pub selected_feature: Option<usize>,
    /// 蛋白 `Properties` 页展开的明细行（0 = 分子量、1 = 消光系数、2 = 吸光度）。
    pub props_expand: [bool; 3],
    /// 蛋白 `Properties` 页「Charge at pH」下拉选中的 pH。
    pub props_ph: f64,
    /// 序列视图排版（`line_bp` 按窗口宽度自适应）。
    pub seq: hx_render::SeqView,
    pub seq_zoom: f32,
    pub cursor_bp: Option<usize>,
    /// 鼠标悬停处的碱基（0-based）。
    pub hover_bp: Option<usize>,
    /// 选段（当前序列坐标，0-based 半开；环状分子 s > e 表示跨原点）。
    pub selection: Option<(usize, usize)>,
    /// 被点击的翻译残基（原始特征下标、密码子下标），用于区分重复特征。
    pub digest_selection: Option<crate::gel_view::DigestSelection>,
    pub selected_residue: Option<(usize, usize)>,
    /// 拖选锚点。
    pub sel_anchor: Option<usize>,
    /// 编辑输入框（插入/替换用的序列文本）。
    pub paste: String,
    /// 酶位点缓存。
    pub sites: SiteHits,
    pub sites_dirty: bool,
    pub enzyme_preset: EnzymePreset,
    pub enzyme_search: String,
    /// 中央面板可用尺寸（每帧更新）。
    pub enzyme_parent: hx_enzymes::sets::Parent,
    pub enzyme_filter: hx_enzymes::sets::Filter,
    pub last_size: Option<egui::Vec2>,
    /// 上一帧中央区（画布）的矩形：给 New File 面板之类「贴文档区左上角」的窗口定位。
    pub central_rect: Option<egui::Rect>,
    /// 打开新文件后待触发一次「适配窗口」（需等首帧尺寸）。
    pub needs_fit: bool,
    /// 没有原生保存对话框时的内部「Save As」窗口。
    pub show_save_dialog: bool,
    /// 内部 Save As 窗口里的路径。
    pub save_as_path: String,
    pub show_about: bool,
    pub show_shortcuts: bool,
    pub show_info: bool,
    /// 待写入剪贴板的文本（`ui()` 末统一提交，`Context` 写剪贴板需要它）。
    pub pending_copy: Option<String>,
    pending_paste_reverse:Option<bool>,
    last_title: Option<String>,
    /// 正在拖动的文档页签（下标）与拖动时的插入位置（0..=tabs.len()）。
    tab_drag: Option<usize>,
    tab_drag_at: Option<usize>,
}

/// 蛋白特征的三个显示类别（SnapGene 蛋白视图左侧工具条的三段开关）。
///
/// oracle `/tmp/opencode/sg/hvp-*.png`：`Show "Region" features` /
/// `Show "Site" and "Het" features` / `Show "Bond" features`。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FeatGroup {
    /// 其余一切（`Region` / `Conflict` / `Mature chain` / …）= 条。
    Region,
    /// `Site` + `Het`（画在字母上方的细刻度）。
    SiteHet,
    /// `Bond`（二硫键一类的连接特征）。
    Bond,
}

/// 按 `type` 属性判定蛋白特征的类别（SnapGene 的分类，大小写不敏感）。
pub fn feat_group(kind_str: &str) -> FeatGroup {
    if kind_str.eq_ignore_ascii_case("site") || kind_str.eq_ignore_ascii_case("het") {
        FeatGroup::SiteHet
    } else if kind_str.eq_ignore_ascii_case("bond") {
        FeatGroup::Bond
    } else {
        FeatGroup::Region
    }
}

/// 左侧竖排工具条按钮的动作（按钮表由 `HelixApp::palette_items` 按视图组装）。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PalAction {
    /// DNA 图谱/序列：显示酶切位点。
    Sites,
    /// 核苷酸：显示特征。
    Features,
    /// DNA：显示引物。
    Primers,
    /// 显示翻译与 ORF。
    Translation,
    /// 显示特征颜色。
    Colors,
    /// 蛋白：`Region` 特征。
    Regions,
    /// 蛋白：`Site` + `Het` 特征。
    SitesHet,
    /// 蛋白：`Bond` 特征。
    Bonds,
    /// 蛋白图谱：条下起止数字。
    TransNumbering,
}

impl PalAction {
    fn apply(self, app: &mut HelixApp) {
        let slot = match self {
            PalAction::Sites => &mut app.show_sites,
            PalAction::Features => &mut app.show_features,
            PalAction::Primers => &mut app.show_primers,
            PalAction::Translation => &mut app.show_translation,
            PalAction::Colors => &mut app.show_colors,
            PalAction::Regions => &mut app.show_regions,
            PalAction::SitesHet => &mut app.show_sites_het,
            PalAction::Bonds => &mut app.show_bonds,
            PalAction::TransNumbering => &mut app.show_trans_numbering,
        };
        *slot = !*slot;
    }
}

impl HelixApp {
    /// 所有字段的默认值（主窗口与克隆预览窗口共用）。
    fn blank() -> Self {
        Self {
            doc: None,
            viewer: None,
            tabs: Vec::new(),
            active: 0,
            status: "Drop a .dna file here, or use File → Open File…".to_string(),
            open_path: String::new(),
            save_path: String::new(),
            recent: Vec::new(),
            tab: Tab::Map,
            rotation: 0.0,
            zoom: 1.0,
            pan_bp: 0.0,
            show_sites: true,
            show_ruler: true,
            show_colors: true,
            show_complement: true,
            show_translation: true,
            translation_frames: 63,
            orf_min_bp: hx_algo::orf::DEFAULT_MIN_BP,
            hybridization_draft: None,
            show_translation_options: false,
            enzyme_preferences_open: false,
            preferences: crate::preferences::Preferences::default(),
            rename_file: None,
            enzyme_manager: Default::default(),
            gel_view: None,
            restoring_workspace: false,
            persistence_enabled: false,
            sequence_color_mode: 0,
            aa_one_letter: false,
            show_features: true,
            show_primers: true,
            show_regions: true,
            show_sites_het: true,
            show_bonds: true,
            show_trans_numbering: true,
            struct_zoom: 0.0,
            struct_pan: egui::Vec2::ZERO,
            struct_sel: 0,
            struct_sel_base: None,
            struct_show_index: true,
            struct_fit: 1.0,
            selected_feature: None,
            props_expand: [false; 3],
            props_ph: 7.0,
            seq: hx_render::SeqView::default(),
            seq_zoom: 1.0,
            cursor_bp: None,
            hover_bp: None,
            selection: None,
            digest_selection: None,
            selected_residue: None,
            sel_anchor: None,
            paste: String::new(),
            sites: SiteHits::default(),
            sites_dirty: true,
            enzyme_preset: EnzymePreset::Unique6Nonredundant,
            enzyme_search: String::new(),
            enzyme_parent: hx_enzymes::sets::Parent::Nonredundant,
            enzyme_filter: hx_enzymes::sets::Filter::UniqueSixPlus,
            map_layout: hx_render::MapLayout::Auto,
            show_toolbar: true,
            sidebar_open: true,
            sidebar_width:crate::chrome::SIDEBAR_W,
            file_filter: String::new(),
            tree_filter: String::new(),
            folder: None,
            folder_input: String::new(),
            expanded: std::collections::HashSet::new(),
            dir_cache: std::collections::HashMap::new(),
            folder_error: None,
            show_side_tools: true,
            cloning_canvas: false,
            show_end_labels: true,
            show_view_options: false,
            show_paste_dialog: false,
            molecular_action: None,
            pending_delete: None,
            insert_unknown: false,
            insert_unknown_count: 1000,
            insert_reverse: false,
            new_file: None,
            feature_dlg: None,
            primer_dlg: None,
            primer_choice: None,
            primer_search: String::new(),
            primer_expanded: std::collections::BTreeSet::new(),
            primer_checked: std::collections::BTreeSet::new(),
            design_dlg: None,
            cloning: None,
            gg_settings: crate::ggsettings::GgSettings::default(),
            gg_settings_open: false,
            gg_settings_backup: None,
            cloning_settings: crate::settings::CloningSettings::default(),
            primer_analysis: None,
            feature_color_popup: false,
            show_open_dialog: None,
            hidden_features: Default::default(),
            hidden_enzymes: Default::default(),
            highlighted_enzymes: Default::default(),
            selected_enzymes: Default::default(),
            selected_features: Default::default(),
            selected_primer: None,
            sites_as_lines: false,
            seq_line_bp: None,
            seq_pan_x: 0.0,
            expanded_enzymes: Default::default(),
            expanded_features: Default::default(),
            feature_search: String::new(),
            tbl_features: TableState::default(),
            tbl_enzymes: TableState::default(),
            tbl_primers: TableState::default(),
            last_size: None,
            central_rect: None,
            needs_fit: false,
            show_save_dialog: false,
            save_as_path: String::new(),
            show_about: false,
            show_shortcuts: false,
            show_info: false,
            pending_copy: None,
            pending_paste_reverse:None,
            last_title: None,
            tab_drag: None,
            tab_drag_at: None,
        }
    }

    /// 克隆对话框里的独立文档视图。它有自己的 doc / 页签 / 缩放 / 选区，
    /// 主窗口的 `HelixApp` 完全不受影响。
    pub fn preview(ctx: &egui::Context, doc: crate::doc::Document) -> Self {
        theme::apply(ctx);
        let mut app = Self::blank();
        app.doc = Some(doc);
        app.sidebar_open = false;
        app.show_toolbar = false;
        app.status = "Cloning preview".into();
        // 预览实例不走主窗口的 `ui()`，酶位点要在这里自己算一遍，
        // 否则克隆对话框的 Enzymes 页和地图上的酶标签是空的。
        let preset = app.enzyme_preset;
        app.apply_preset(preset);
        app
    }

    pub fn new(cc: &eframe::CreationContext<'_>, paths: Vec<PathBuf>) -> Self {
        theme::apply(&cc.egui_ctx);
        let mut app = Self::blank();
        // 持久化设置：文件夹 + 最近文件（先载入，再打开命令行给出的文件）
        app.persistence_enabled=true;
        app.restoring_workspace=true;
        let st = crate::settings::Settings::load();
        app.preferences = crate::preferences::Preferences::load();
        app.preferences.language.install(&cc.egui_ctx);
        app.seq_zoom = app.preferences.font_scale;
        app.sites_as_lines = app.preferences.enzyme_lines;
        if !app.preferences.default_directory.trim().is_empty() {
            app.open_path=app.preferences.default_directory.clone();
            app.save_path=app.preferences.default_directory.clone();
        }
        if app.preferences.fullscreen { cc.egui_ctx.send_viewport_cmd(egui::ViewportCommand::Fullscreen(true)); }
        if let Some(f) = &st.folder {
            app.folder_input = f.display().to_string();
        }
        app.folder = st.folder;
        // Restoring the root path alone does not restore its tree state.
        app.restore_folder_tree();
        for p in st.recent.into_iter().rev() {
            if !app.recent.contains(&p) {
                app.recent.insert(0, p);
            }
        }
        // 显示/工具栏状态与 SnapGene 的 settings.ini 一样跨启动保留。
        if let Some(names) = &st.display {
            app.apply_display_preferences(names);
        }
        if let Some(v) = st.sidebar_open {
            app.sidebar_open = v;
        }
        if let Some(v) = st.show_toolbar {
            app.show_toolbar = v;
        }
        if let Some(v) = st.show_side_tools {
            app.show_side_tools = v;
        }
        // 克隆参数与 Golden Gate 酶集合（SG `[actions]` / 设置窗口）。
        app.cloning_settings = st.cloning;
        if let Some(gg) = st.gg {
            app.gg_settings = gg;
        }
        let (restore,restore_active)=if paths.is_empty()&&app.preferences.restore_workspace{crate::preferences::Preferences::read_workspace()}else{(Vec::new(),0)};
        let restoring=!restore.is_empty();
        for p in paths.into_iter().chain(restore) {
            app.open(&p);
            let preset = app.enzyme_preset;
            app.apply_preset(preset);
        }
        if restoring && restore_active<app.tabs.len(){app.activate_tab(restore_active);}
        app.restoring_workspace=false;
        app.save_settings();
        app
    }

    /// 当前文档是否为蛋白（菜单与左侧工具条按分子类型切换）。
    pub fn is_protein(&self) -> bool {
        self.doc.as_ref().is_some_and(|d| d.is_protein())
    }

    /// 当前文档是否为 DNA。
    pub fn is_dna(&self) -> bool {
        self.doc.as_ref().is_some_and(|d| d.is_dna())
    }

    pub fn displayed_orfs(&self,seq:&[u8],circular:bool)->Vec<hx_algo::orf::Orf>{
        let index=if self.is_dna(){0}else{1};let defaults=&self.preferences.molecule_defaults[index];
        let starts=if defaults.alternative_starts{vec![*b"ATG",*b"GTG",*b"TTG"]}else{vec![*b"ATG"]};
        hx_algo::orf::find_orfs_with_options(seq,circular,self.orf_min_bp,&starts,defaults.require_start,defaults.allow_end_orfs)
    }

    /// 当前文档里该**特征类别**是否允许显示（不含 `hidden_features` 的判定）。
    ///
    /// 核苷酸只有一个 `Features` 开关；蛋白按 SnapGene 的
    /// `Region` / `Site`+`Het` / `Bond` 三段分开（oracle `/tmp/opencode/sg/hvp-*.png`）。
    pub fn feature_type_shown(&self, f: &hx_formats::annot::FeatureAnnotation) -> bool {
        match self.doc.as_ref().map(|d| d.kind()) {
            Some(hx_core::MoleculeKind::Protein) => match feat_group(&f.kind) {
                FeatGroup::Region => self.show_regions,
                FeatGroup::SiteHet => self.show_sites_het,
                FeatGroup::Bond => self.show_bonds,
            },
            _ => self.show_features,
        }
    }

    /// 特征是否可以画：类别开关 + `hidden_features`（文件 `visible="0"` 也落在后者里）。
    pub fn feature_visible(&self, idx: usize, f: &hx_formats::annot::FeatureAnnotation) -> bool {
        !self.hidden_features.contains(&idx) && self.feature_type_shown(f)
    }

    /// 把当前文件夹设为 `dir`（并清缓存、持久化）。
    ///
    /// **只有**「Open Folder…」/ 路径框回车这类**显式**动作才调用它：打开文件、切标签页
    /// 都不会换侧栏根目录（用户裁定 2026-10-18：打开过的文件夹就是左侧区固定的根）。
    pub fn set_folder(&mut self, dir: Option<std::path::PathBuf>) {
        self.folder = dir.clone();
        self.expanded.clear();
        self.dir_cache.clear();
        self.folder_error = None;
        self.tree_filter.clear();
        if let Some(d) = &dir {
            self.folder_input = d.display().to_string();
        }
        self.save_settings();
        // 新根目录立刻展开到当前文档（SnapGene 实测：打开文档后树展开到它的路径）
        if let Some(p) = self.doc.as_ref().and_then(|d| d.path.clone()) {
            self.expand_ancestors(&p);
        }
    }

    /// 把「当前文档的祖先文件夹」逐个加入展开集合（从根一路到它的父目录）。
    ///
    /// SnapGene 实测截图 26：根下只有 `4. Sample Files` 自动展开，因为打开的文档在里面。
    /// 文档不在根目录里时什么都不做。
    pub fn expand_ancestors(&mut self, path: &std::path::Path) {
        let Some(root) = self.folder.clone() else {
            return;
        };
        let Some(mut cur) = path.parent().map(|p| p.to_path_buf()) else {
            return;
        };
        if !cur.starts_with(&root) {
            return;
        }
        loop {
            let done = cur == root;
            self.expanded.insert(cur.clone());
            if done {
                break;
            }
            match cur.parent() {
                Some(p) if p.starts_with(&root) => cur = p.to_path_buf(),
                _ => break,
            }
        }
    }

    /// 刷新当前文件夹的目录列表缓存。
    pub fn refresh_folder(&mut self) {
        self.dir_cache.clear();
        self.folder_error = None;
        if let Some(d) = self.folder.clone() {
            match crate::settings::list_dir(&d, 2000) {
                Ok(list) => {
                    for e in &list {
                        let _ = e;
                    }
                    self.dir_cache.insert(d, list);
                }
                Err(err) => self.folder_error = Some(err),
            }
        }
    }

    /// 取某目录的列表（带缓存）。
    pub fn dir_listing(&mut self, dir: &std::path::Path) -> &Vec<crate::settings::Entry> {
        if !self.dir_cache.contains_key(dir) {
            match crate::settings::list_dir(dir, 2000) {
                Ok(list) => {
                    self.dir_cache.insert(dir.to_path_buf(), list);
                }
                Err(err) => {
                    self.folder_error = Some(err);
                    self.dir_cache.insert(dir.to_path_buf(), Vec::new());
                }
            }
        }
        &self.dir_cache[dir]
    }

    /// 写回设置（文件夹 / 最近文件 / 显示开关 / 工具栏与侧栏可见性）。
    pub fn save_settings(&self) {
        if !self.persistence_enabled||self.restoring_workspace{return;}
        let _ = self.preferences.save();
        let (files,active)=self.workspace_snapshot();
        let _=crate::preferences::Preferences::save_workspace(&files,active);
        let mut cloning = self.cloning_settings.clone();
        cloning.golden_gate_default_enzyme =
            self.gg_settings.default_enzyme_name().to_string();
        let st = crate::settings::Settings {
            folder: self.folder.clone(),
            recent: self.recent.clone(),
            display: Some(self.display_preferences()),
            sidebar_open: Some(self.sidebar_open),
            show_toolbar: Some(self.show_toolbar),
            show_side_tools: Some(self.show_side_tools),
            cloning,
            gg: Some(self.gg_settings.clone()),
        };
        st.save();
    }

    fn restore_folder_tree(&mut self){
        if let Some(folder)=self.folder.clone(){
            self.folder_input=folder.display().to_string();
            self.expanded.insert(folder.clone());
            let _=self.dir_listing(&folder);
        }
    }

    fn workspace_snapshot(&self)->(Vec<PathBuf>,usize){
        let mut paths=Vec::new();let mut active=0;
        for(i,t)in self.tabs.iter().enumerate(){
            let path=if i==self.active{self.doc.as_ref().and_then(|d|d.path.clone()).or_else(||self.viewer.as_ref().map(|v|v.path().to_path_buf()))}
                else{t.doc.as_ref().and_then(|d|d.path.clone()).or_else(||t.viewer.as_ref().map(|v|v.path().to_path_buf()))};
            if let Some(path)=path{if i==self.active{active=paths.len();}paths.push(path);}
        }
        (paths,active)
    }

    fn relocate_open_file(&mut self,old:&std::path::Path,new:&std::path::Path){
        if let Some(doc)=self.doc.as_mut(){if doc.path.as_deref()==Some(old){doc.path=Some(new.to_path_buf());}}
        fn update_viewer(v:&mut crate::viewer::ViewerDoc,old:&std::path::Path,new:&std::path::Path){if v.path()==old{match v{crate::viewer::ViewerDoc::Text(d)=>d.path=new.to_path_buf(),crate::viewer::ViewerDoc::Pdf(d)=>d.path=new.to_path_buf()}}}
        if let Some(v)=self.viewer.as_mut(){update_viewer(v,old,new);}
        for tab in &mut self.tabs{if let Some(doc)=tab.doc.as_mut(){if doc.path.as_deref()==Some(old){doc.path=Some(new.to_path_buf());}}if let Some(v)=tab.viewer.as_mut(){update_viewer(v,old,new);}}
        for path in &mut self.recent{if path==old{*path=new.to_path_buf();}}
        if self.open_path==old.display().to_string(){self.open_path=new.display().to_string();}
        if self.save_path==old.display().to_string(){self.save_path=new.display().to_string();}
        self.dir_cache.clear();self.expand_ancestors(new);self.save_settings();
    }

    /// 当前启用的显示开关（稳定名字；只列 `true`，默认值由 [`Self::blank`] 决定）。
    ///
    /// SnapGene 同样把 `show\regions` / `show\translations` / `show\bonds`
    /// 写进 `settings.ini` 并在下次启动恢复；这里覆盖 HX 现有的全部显示开关。
    pub fn display_preferences(&self) -> Vec<String> {
        [
            ("sites", self.show_sites),
            ("features", self.show_features),
            ("primers", self.show_primers),
            ("translation", self.show_translation),
            ("colors", self.show_colors),
            ("ruler", self.show_ruler),
            ("complement", self.show_complement),
            ("aa_one_letter", self.aa_one_letter),
            ("regions", self.show_regions),
            ("sites_het", self.show_sites_het),
            ("bonds", self.show_bonds),
            ("trans_numbering", self.show_trans_numbering),
            ("end_labels", self.show_end_labels),
        ]
        .into_iter()
        .filter_map(|(name, on)| on.then(|| name.to_string()))
        .collect()
    }

    /// 恢复显示开关。列表是完整快照：未列出的开关一律关闭。
    pub fn apply_display_preferences(&mut self, names: &[String]) {
        let on = |name: &str| names.iter().any(|n| n == name);
        self.show_sites = on("sites");
        self.show_features = on("features");
        self.show_primers = on("primers");
        self.show_translation = on("translation");
        self.show_colors = on("colors");
        self.show_ruler = on("ruler");
        self.show_complement = on("complement");
        self.aa_one_letter = on("aa_one_letter");
        self.show_regions = on("regions");
        self.show_sites_het = on("sites_het");
        self.show_bonds = on("bonds");
        self.show_trans_numbering = on("trans_numbering");
        self.show_end_labels = on("end_labels");
    }

    /// 打开一个文件到**新的标签页**（SnapGene：同一窗口内多标签页；已打开则切过去）。
    ///
    /// 内容嗅探：SnapGene 原生 TLV（`.dna/.rna/.prot/.dnaln/.praln/.gel`，以及
    /// 库里那种 `.ftrs`）走原生路径；FASTA / GenBank 文本按导入处理，
    /// 成为**未命名文档**（首次保存才落成 `.dna`），与 SnapGene 的「Open Files…」一致。
    pub fn open(&mut self, path: &std::path::Path) {
        if let Some(i) = self.tab_of_path(path) {
            self.activate_tab(i);
            return;
        }
        // `.txt` / `.pdf`：只读阅读视图，绝不进入序列导入/编辑流程。
        if crate::viewer::kind_of_path(path).is_some() {
            self.open_viewer(path);
            return;
        }
        // FASTA / GenBank：导入
        if let Ok(bytes) = std::fs::read(path) {
            let native = bytes.starts_with(b"\x09") && bytes.len() > 8 && &bytes[5..13] == b"SnapGene";
            if !native && hx_formats::import::detect_format(&bytes).is_some() {
                self.import_bytes(bytes, Some(path.to_path_buf()));
                return;
            }
        }
        match Document::from_path(path) {
            Ok(doc) => self.push_doc(doc, Some(path.to_path_buf())),
            Err(e) => {
                self.status = format!("Open failed: {e}");
            }
        }
    }

    /// 打开 `.txt` / `.pdf` 为**只读**页签（不创建 `Document`，不参与保存/编辑）。
    pub fn open_viewer(&mut self, path: &std::path::Path) {
        if let Some(i) = self.tab_of_path(path) {
            self.activate_tab(i);
            return;
        }
        match crate::viewer::ViewerDoc::open(path) {
            Ok(Some(v)) => self.push_viewer(v, path.to_path_buf()),
            Ok(None) => {
                // 扩展名在调用与读取之间变了：退回普通打开流程。
                if let Ok(doc) = Document::from_path(path) {
                    self.push_doc(doc, Some(path.to_path_buf()));
                } else {
                    self.status = format!("Open failed: {} is not a viewable file", path.display());
                }
            }
            Err(e) => self.status = format!("Open failed: {e}"),
        }
    }

    /// 把只读阅读器推进标签页并激活。
    pub fn push_viewer(&mut self, viewer: crate::viewer::ViewerDoc, path: PathBuf) {
        let title = viewer.title();
        let summary = viewer.summary();
        // 当前标签页的文档/阅读器存回它自己的条目
        if !self.tabs.is_empty() {
            let v = DocViewState::take(self);
            if let Some(d) = self.doc.take() {
                self.tabs[self.active].doc = Some(d);
            }
            if let Some(vw) = self.viewer.take() {
                self.tabs[self.active].viewer = Some(vw);
            }
            self.tabs[self.active].view = v;
        } else {
            DocViewState::default().apply(self);
        }
        self.tabs.push(TabDoc {
            doc: None,
            viewer: None,
            view: DocViewState::default(),
            untitled_name: None,
        });
        self.active = self.tabs.len() - 1;
        self.doc = None;
        self.viewer = Some(viewer);
        DocViewState::default().apply(self);
        self.save_path = path.display().to_string();
        self.open_path = path.display().to_string();
        self.expand_ancestors(&path);
        self.recent.retain(|r| r != &path);
        self.recent.insert(0, path.clone());
        self.recent.truncate(12);
        self.save_settings();
        self.status = format!("Opened {title}: {summary}");
    }

    /// 导入 FASTA / GenBank 字节；每条记录一个标签页（SnapGene 也是各占一个文档）。
    /// 导入的文档没有路径：第一次 Save/Save As 时才写成 `.dna`。
    pub fn import_bytes(&mut self, bytes: Vec<u8>, src: Option<PathBuf>) {
        let recs = match hx_formats::import::parse(&bytes) {
            Ok(r) => r,
            Err(e) => {
                self.status = format!("Import failed: {e}");
                return;
            }
        };
        let stem = src
            .as_ref()
            .and_then(|p| p.file_stem())
            .map(|s| s.to_string_lossy().to_string());
        let mut opened: Vec<String> = Vec::new();
        for rec in &recs {
            let dna = hx_formats::import::to_dna_bytes(rec);
            // 未命名文档的显示名：记录名（单记录且文件名更有信息量时用文件名）
            let hint = if recs.len() == 1 {
                stem.clone().unwrap_or_else(|| rec.name.clone())
            } else {
                rec.name.clone()
            };
            match Document::from_imported(dna, hint) {
                Ok(doc) => {
                    let feats = doc.features().len();
                    let len = doc.len();
                    let unit = doc.unit();
                    self.push_doc(doc, None);
                    opened.push(format!("{} ({} {unit}, {feats} features)", rec.name, nfmt(len)));
                }
                Err(e) => self.status = format!("Import failed: {e}"),
            }
        }
        if !opened.is_empty() {
            self.status = format!("Imported {}", opened.join("; "));
        }
    }

    /// 该路径是否已经在某个标签页里（返回标签页下标）。
    pub fn tab_of_path(&self, path: &std::path::Path) -> Option<usize> {
        if self
            .doc
            .as_ref()
            .and_then(|d| d.path.as_deref())
            .is_some_and(|p| p == path)
        {
            return Some(self.active);
        }
        if self.viewer.as_ref().is_some_and(|v| v.path() == path) {
            return Some(self.active);
        }
        self.tabs.iter().position(|t| {
            t.doc.as_ref().and_then(|d| d.path.as_deref()) == Some(path)
                || t.viewer.as_ref().is_some_and(|v| v.path() == path)
        })
    }

    /// 把文档推进标签页并激活（新开 / 拖放 / Open File… / 命令行都用它）。
    pub fn push_doc(&mut self, mut doc: Document, path: Option<PathBuf>) {
        let use_defaults = path.is_none() && doc.native_class().is_none();
        if use_defaults && !self.preferences.default_author.trim().is_empty() {
            let _ = doc.session.set_note("Author", &self.preferences.default_author);
        }
        if use_defaults && doc.kind().is_nucleotide() {
            let params=if doc.len()>self.preferences.primer_large_threshold {self.preferences.primer_large}else{self.preferences.primer_small};
            let _=doc.session.apply_primer(hx_formats::primer_edit::PrimerEdit::Hybridization(params));
        }
        let n = doc.len();
        let name = doc.display_name();
        let unit = doc.unit();
        let feats = doc.features().len();
        let primers = doc.primers().len();
        // 文件里标了 `visible="0"` 的特征（SnapGene 的隐藏特征）默认不画
        let invisible: Vec<usize> = doc
            .features()
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.visible)
            .map(|(i, _)| i)
            .collect();
        // 当前文档存回它自己的标签页
        if !self.tabs.is_empty() {
            let v = DocViewState::take(self);
            if let Some(d) = self.doc.take() {
                self.tabs[self.active].doc = Some(d);
            }
            if let Some(vw) = self.viewer.take() {
                self.tabs[self.active].viewer = Some(vw);
            }
            self.tabs[self.active].view = v;
        } else {
            // 空状态：把字段上的残留状态清零
            DocViewState::default().apply(self);
        }
        // 新标签页继承当前的显示集合（Custom 时保留用户选的酶，预设则由 `apply_preset` 重新筛选）
        let mut view = DocViewState::default();
        let preset_to_inherit = self.enzyme_preset;
        view.enzyme_preset = preset_to_inherit;
        self.tabs.push(TabDoc {
            doc: None,
            viewer: None,
            view,
            untitled_name: None,
        });
        self.active = self.tabs.len() - 1;
        self.doc = Some(doc);
        self.viewer = None;
        DocViewState::default().apply(self);
        if use_defaults {
            let index=match self.doc.as_ref().unwrap().kind(){hx_core::MoleculeKind::Dna=>0,hx_core::MoleculeKind::Rna=>1,hx_core::MoleculeKind::Protein=>2};
            let d=&self.preferences.molecule_defaults[index];
            self.show_features=d.features;self.show_primers=d.primers;self.show_translation=d.translations;
            self.show_regions=d.regions;self.show_sites_het=d.sites;self.show_bonds=d.bonds;
            self.orf_min_bp=d.min_orf_aa.saturating_mul(3);self.aa_one_letter=d.one_letter;
            self.seq_line_bp=if d.line_width==0{None}else{Some(d.line_width)};self.struct_show_index=d.structure_index;
            self.seq_zoom=self.preferences.font_scale;self.show_info=self.preferences.info_imported;
        }
        for i in invisible {
            self.hidden_features.insert(i);
        }
        self.enzyme_preset = preset_to_inherit;
        if self.enzyme_preset != EnzymePreset::Custom {
            let preset = self.enzyme_preset;
            self.apply_preset(preset);
        }
        self.needs_fit = true;
        if let Some(p) = &path {
            self.save_path = p.display().to_string();
            self.open_path = p.display().to_string();
            // 侧栏根目录**不跟随**打开的文件（用户裁定 2026-10-18：打开过的文件夹
            // 就是左侧区固定的根）；只把树的展开状态推进到该文档，让它显出来。
            self.expand_ancestors(p);
            self.recent.retain(|r| r != p);
            self.recent.insert(0, p.clone());
            self.recent.truncate(12);
            self.save_settings();
        } else if let Some(p) = self.default_save_path() {
            // 未保存文档（克隆产物、导入记录、新建序列）：默认落点是侧栏
            // 根目录 + 文档名，不能沿用上一个文档的 `save_path`。
            self.save_path = p.display().to_string();
        }
        // 原生只读文档（比对/凝胶）没有序列，写类别与字节数更有信息量
        let native = self
            .doc
            .as_ref()
            .and_then(|d| d.native_class().map(|c| (c, d.native_bytes().unwrap_or(0))));
        self.status = if let Some((class, bytes)) = native {
            format!("Opened {name}: SnapGene {class} document, {} bytes", nfmt(bytes))
        } else if unit == "bp" {
            format!("Opened {name}: {n} bp, {feats} features, {primers} primers")
        } else {
            format!("Opened {name}: {n} {unit}, {feats} features")
        };
    }

    /// 无路径的新文档（File ▸ New DNA or RNA File…）。
    pub fn new_untitled(&mut self) {
        self.open_new_file_dialog(hx_core::MoleculeKind::Dna);
    }

    // ---- New DNA/RNA/Protein File 对话框 ----

    /// 打开 New File 面板（SnapGene 8 的 New File pane）。已开着时只切换分子类型。
    pub fn open_new_file_dialog(&mut self, kind: hx_core::MoleculeKind) {
        match self.new_file.as_mut() {
            Some(st) => st.set_kind(kind),
            None => self.new_file = Some(crate::newfile::NewFileState::new(kind)),
        }
        self.sync_new_file_preview();
    }

    /// 预览文档与对话框状态同步（内容变了才重建；SnapGene 在主窗口里实时预览）。
    fn sync_new_file_preview(&mut self) {
        let Some(st) = self.new_file.as_ref() else {
            return;
        };
        let sig = st.signature();
        let first = st.preview_tab.is_none();
        if !first && st.sig == sig {
            return;
        }
        let doc = match crate::newfile::build_doc(st) {
            Ok(d) => d,
            Err(e) => {
                self.status = format!("New file failed: {e}");
                return;
            }
        };
        let tab = st.preview_tab;
        match tab {
            None => {
                self.push_doc(doc, None);
                let i = self.active;
                // 预览期的页签名是 `New Sequence`（SnapGene 实测）
                self.tabs[i].untitled_name = Some("New Sequence".to_string());
                let (n, unit) = self
                    .doc
                    .as_ref()
                    .map(|d| (d.len(), d.unit()))
                    .unwrap_or((0, "bp"));
                self.status = format!("New file: {n} {unit} — not saved yet");
                if let Some(st) = self.new_file.as_mut() {
                    st.preview_tab = Some(i);
                    st.sig = sig;
                }
            }
            Some(i) => {
                if i >= self.tabs.len() {
                    // 预览页被用户关掉了
                    if let Some(st) = self.new_file.as_mut() {
                        st.preview_tab = None;
                    }
                    return;
                }
                if self.active != i {
                    self.activate_tab(i);
                }
                self.doc = Some(doc);
                DocViewState::default().apply(self);
                self.needs_fit = true;
                let (n, unit, topo) = self
                    .doc
                    .as_ref()
                    .map(|d| {
                        (
                            d.len(),
                            d.unit(),
                            if d.circular() { "circular" } else { "linear" },
                        )
                    })
                    .unwrap_or((0, "bp", "linear"));
                self.status = format!("New file: {n} {unit}, {topo} — not saved yet");
                if let Some(st) = self.new_file.as_mut() {
                    st.sig = sig;
                }
            }
        }
    }

    /// 该标签页是不是 New File 的预览页。
    fn is_preview_tab(&self, i: usize) -> bool {
        self.tabs
            .get(i)
            .is_some_and(|t| t.untitled_name.as_deref() == Some("New Sequence"))
    }

    /// 对话框「Cancel」：丢弃预览文档。
    pub fn new_file_cancel(&mut self) {
        let Some(st) = self.new_file.take() else {
            return;
        };
        if let Some(i) = st.preview_tab {
            if i < self.tabs.len() && self.is_preview_tab(i) {
                let n = self.tabs.len();
                self.close_tab(i);
                let _ = n;
            }
        }
        self.status = "New file cancelled".to_string();
    }

    /// 对话框「Create」：预览文档变成真正的（未保存）文档，名字带扩展名。
    pub fn new_file_create(&mut self) {
        let Some(st) = self.new_file.take() else {
            return;
        };
        let mut name = st.name.trim().to_string();
        if name.is_empty() {
            name = crate::newfile::default_name(st.kind);
        } else if !name.contains('.') {
            // 没写扩展名：按分子类型补（SnapGene 用文件类型决定扩展名）
            name = format!("{name}.{}", crate::newfile::default_name(st.kind).rsplit('.').next().unwrap_or("dna"));
        }
        let doc = match Document::from_imported(st.build_bytes(), name.clone()) {
            Ok(d) => d,
            Err(e) => {
                self.status = format!("Create failed: {e}");
                return;
            }
        };
        let (n, unit, feats) = (doc.len(), doc.unit(), doc.features().len());
        match st.preview_tab {
            Some(i) if i < self.tabs.len() => {
                if self.active != i {
                    self.activate_tab(i);
                }
                self.doc = Some(doc);
                self.tabs[i].untitled_name = None;
                DocViewState::default().apply(self);
                self.needs_fit = true;
            }
            _ => self.push_doc(doc, None),
        }
        // 未保存文档：首次 Save/Save As 的默认落点是侧栏根目录 + 这个文件名
        let dir = self
            .folder
            .clone()
            .unwrap_or_else(|| std::env::current_dir().unwrap_or_default());
        self.save_path = dir.join(&name).display().to_string();
        self.status = format!("Created {name}: {n} {unit}, {feats} features — not saved yet");
    }

    // ---- Add/Edit Feature 对话框（SnapGene 的 Features 面板）----

    /// 「Add Feature…」：有选段就用选段（1-based 闭），否则用光标处的单碱基。
    pub fn open_add_feature_dialog(&mut self, translated: bool) {
        let Some(doc) = self.doc.as_ref() else {
            self.status = "No file open — nothing to annotate".into();
            return;
        };
        if doc.native_class().is_some() {
            self.status = "This document type cannot be edited".into();
            return;
        }
        let len = doc.len();
        if len == 0 {
            self.status = "Cannot add a feature to an empty sequence".into();
            return;
        }
        let (s, e) = match self.selection {
            Some((a, b)) if b > a => (a.min(len) + 1, b.min(len)),
            _ => {
                let at = self.cursor_bp.unwrap_or(0).min(len - 1);
                (at + 1, at + 1)
            }
        };
        let name = doc.display_name();
        let mut dlg = crate::featuredlg::FeatureDlg::add(&name, len, s, e);
        if translated {
            dlg.spec.translate = true;
        }
        self.feature_dlg = Some(dlg);
    }

    /// 「Edit Feature…」：编辑 Features 面板里选中的特征。
    pub fn open_edit_feature_dialog(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            self.status = "No file open".into();
            return;
        };
        let Some(i) = self.selected_feature else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let Some(spec) = doc.feature_specs().get(i).cloned() else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let (name, len) = (doc.display_name(), doc.len());
        self.feature_dlg = Some(crate::featuredlg::FeatureDlg::edit(&name, len, i, spec));
    }

    /// 「Duplicate Feature…」：以选中特征的属性打开 Add 对话框（得到一个新的 recentID）。
    pub fn open_duplicate_feature_dialog(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            self.status = "No file open".into();
            return;
        };
        let Some(i) = self.selected_feature else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let Some(spec) = doc.feature_specs().get(i).cloned() else {
            return;
        };
        let (name, len) = (doc.display_name(), doc.len());
        self.feature_dlg = Some(crate::featuredlg::FeatureDlg::add_spec(&name, len, spec));
    }

    /// 对话框 OK。
    fn feature_dlg_apply(&mut self, edit: FeatureEdit) {
        let (label, is_add) = match &edit {
            FeatureEdit::Add(s) => (
                format!(
                    "Added feature {}{}",
                    if s.name.is_empty() { s.kind.clone() } else { s.name.clone() },
                    s.first_seg()
                        .map(|g| format!(" ({}..{})", g.start, g.end))
                        .unwrap_or_default()
                ),
                true,
            ),
            FeatureEdit::Update(s) => (
                format!(
                    "Edited feature {}",
                    if s.name.is_empty() { s.kind.clone() } else { s.name.clone() }
                ),
                false,
            ),
            FeatureEdit::Remove { .. } => ("Removed feature".to_string(), false),
            FeatureEdit::SetVisible { visible, .. } => (
                if *visible { "Feature shown".to_string() } else { "Feature hidden".to_string() },
                false,
            ),
        };
        let Some(doc) = self.doc.as_mut() else {
            return;
        };
        match doc.apply_feature(edit) {
            Ok(()) => {
                self.sites_dirty = false;
                if is_add {
                    // 新特征在末尾；让它成为当前选中项
                    let n = doc.features().len();
                    self.selected_feature = Some(n.saturating_sub(1));
                    self.hidden_features.clear();
                }
                let feats = doc.features().len();
                self.status = format!("{label} — {feats} features");
                self.clamp_feature_selection();
            }
            Err(e) => self.status = format!("Feature edit failed: {e}"),
        }
        self.feature_color_popup = false;
    }

    /// 「Remove Feature」：删除选中特征。
    pub fn remove_selected_feature(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let Some(i) = self.selected_feature else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let Some(spec) = doc.feature_specs().get(i) else {
            return;
        };
        let (id, name) = (
            spec.recent_id,
            if spec.name.is_empty() { spec.kind.clone() } else { spec.name.clone() },
        );
        self.feature_dlg_apply(FeatureEdit::Remove { recent_id: id });
        self.selected_feature = None;
        self.selected_features.clear();
        self.hidden_features.clear();
        self.status = format!("Removed feature {name}");
    }

    /// 「Show or Hide Feature」：翻转选中特征的 `visible`（写回文件）。
    pub fn toggle_selected_feature_visible(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let Some(i) = self.selected_feature else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let Some(spec) = doc.feature_specs().get(i) else {
            return;
        };
        let edit = FeatureEdit::SetVisible {
            recent_id: spec.recent_id,
            visible: !spec.visible,
        };
        self.feature_dlg_apply(edit);
    }

    /// 「Feature Color…」：把选中特征的所有片段改成 `color`。
    pub fn set_selected_feature_color(&mut self, color: &str) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let Some(i) = self.selected_feature else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let Some(mut spec) = doc.feature_specs().get(i).cloned() else {
            return;
        };
        for s in &mut spec.segments {
            s.color = color.to_string();
        }
        self.feature_dlg_apply(FeatureEdit::Update(spec));
        self.feature_color_popup = false;
    }

    /// 特征对话框窗口（每帧）。
    fn feature_dlg_window(&mut self, ctx: &egui::Context) {
        let (outcome, unchanged) = {
            let Some(dlg) = self.feature_dlg.as_mut() else {
                return;
            };
            // 打开后原样点 OK：不动文档（避免把「没改」记成一次编辑、把文件弄脏）
            let unchanged = dlg.mode == crate::featuredlg::Mode::Edit && !dlg.changed();
            (crate::featuredlg::draw(dlg, ctx), unchanged)
        };
        match outcome {
            Some(crate::featuredlg::Outcome::Apply(e)) => {
                let e = *e;
                self.feature_dlg = None;
                if unchanged {
                    self.status = "No change".into();
                } else {
                    self.feature_dlg_apply(e);
                }
            }
            Some(crate::featuredlg::Outcome::Cancel) => {
                self.feature_dlg = None;
                self.status = "Feature dialog cancelled".into();
            }
            None => {}
        }
    }

    // ---- 引物（Add / Edit / Duplicate / Remove / Analyze / Design） ----

    /// 当前文档的模板序列与拓扑（引物对话框预览用）。
    fn primer_template(&self) -> (Vec<u8>, bool) {
        match self.doc.as_ref() {
            Some(d) => (d.sequence().to_vec(), d.circular()),
            None => (Vec::new(), false),
        }
    }

    /// 该名字在文档里已存在时，追加 ` 2`、` 3`…（SnapGene 的 Duplicate 也会改名）。
    fn unique_primer_name(&self, base: &str) -> String {
        let Some(doc) = self.doc.as_ref() else {
            return base.to_string();
        };
        let existing: Vec<&str> = doc.primers().iter().map(|p| p.name.as_str()).collect();
        if !existing.contains(&base) {
            return base.to_string();
        }
        for i in 2..1000 {
            let cand = format!("{base} {i}");
            if !existing.contains(&cand.as_str()) {
                return cand;
            }
        }
        base.to_string()
    }

    /// `Primers ▸ Add Primer…`（Ctrl+R）。
    pub fn open_add_primer_dialog(&mut self) {
        match self.doc.as_ref() {
            None => self.status = "No document open".into(),
            Some(d) if !d.is_dna() => {
                self.status = "Primers can only be added to DNA sequences".into()
            }
            Some(d) => {
                // SnapGene：序列视图里有选区时先问 `Make a primer from the selection?`
                let sel = self.selection.filter(|(a, b)| b > a && *b <= d.len());
                match sel {
                    Some((a, b)) => self.primer_choice = Some((a, b)),
                    None => {
                        let name = d.display_name();
                        self.primer_dlg = Some(crate::primerdlg::PrimerDlg::add(&name));
                    }
                }
            }
        }
    }

    /// `Make a primer from the selection?` → Top/Bottom Strand。
    fn primer_from_selection(&mut self, strand: u8) {
        let Some((a, b)) = self.primer_choice.take() else {
            return;
        };
        let (template, _) = self.primer_template();
        let doc_name = self
            .doc
            .as_ref()
            .map(|d| d.display_name())
            .unwrap_or_default();
        let name = self.unique_primer_name("Primer 1");
        if let Some(dlg) =
            crate::primerdlg::PrimerDlg::from_selection(&doc_name, &template, a, b, strand, &name)
        {
            self.primer_dlg = Some(dlg);
        }
    }

    /// 「Make a primer from the selection?」对话框。
    fn primer_choice_window(&mut self, ctx: &egui::Context) {
        let Some((a, b)) = self.primer_choice else {
            return;
        };
        let seq: String = self
            .doc
            .as_ref()
            .map(|d| {
                d.sequence()[a..b.min(d.len())]
                    .iter()
                    .map(|x| *x as char)
                    .collect()
            })
            .unwrap_or_default();
        let mut choice: Option<u8> = None;
        let mut cancel = false;
        let mut open = true;
        egui::Window::new("primer_choice")
            .id(egui::Id::unique("primer_choice"))
            .title_bar(false)
            .resizable(false)
            .collapsible(false)
            .open(&mut open)
            .fixed_size(egui::vec2(460.0, 0.0))
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -40.0))
            .frame(
                egui::Frame::window(&ctx.style_of(egui::Theme::Light))
                    .fill(egui::Color32::WHITE)
                    .stroke(egui::Stroke::new(0.8, theme::BORDER))
                    .inner_margin(egui::Margin::symmetric(14, 12)),
            )
            .show(ctx, |ui| {
                ui.set_width(460.0);
                ui.label(
                    egui::RichText::new(crate::i18n::text(ui.ctx(),"Make a primer from the selection?"))
                        .size(13.0)
                        .strong(),
                );
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new(format!("{} .. {}", a + 1, b))
                            .monospace()
                            .size(11.5)
                            .color(theme::TEXT_DIM),
                    );
                    ui.label(egui::RichText::new(&seq).monospace().size(11.5));
                });
                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    if ui.button(crate::i18n::text(ui.ctx(),"Top Strand")).clicked() {
                        choice = Some(0);
                    }
                    if ui.button(crate::i18n::text(ui.ctx(),"Bottom Strand")).clicked() {
                        choice = Some(1);
                    }
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() {
                            cancel = true;
                        }
                    });
                });
            });
        if let Some(strand) = choice {
            self.primer_from_selection(strand);
        } else if cancel || !open {
            self.primer_choice = None;
        }
    }

    /// 选中的引物下标（未选中时：只有一个引物就用它）。
    fn active_primer_index(&self) -> Option<usize> {
        let doc = self.doc.as_ref()?;
        let n = doc.primers().len();
        if n == 0 {
            return None;
        }
        self.selected_primer.filter(|i| *i < n).or(if n == 1 { Some(0) } else { None })
    }

    /// `Primers ▸ Edit Primer…`。
    pub fn open_edit_primer_dialog(&mut self) {
        let Some(i) = self.active_primer_index() else {
            self.status = "No primer selected — pick one in the Primers tab".into();
            return;
        };
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let p = doc.primers()[i].clone();
        let (template, circular) = self.primer_template();
        let doc_name = doc.display_name();
        let mut dlg = crate::primerdlg::PrimerDlg::edit(&p, &template, circular);
        dlg.doc_name = doc_name;
        self.primer_dlg = Some(dlg);
    }

    /// `Primers ▸ Duplicate Primer…`（Ctrl+Shift+R）。
    pub fn duplicate_selected_primer(&mut self) {
        let Some(i) = self.active_primer_index() else {
            self.status = "No primer selected — pick one in the Primers tab".into();
            return;
        };
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let template = doc.sequence().to_vec();
        let mut spec = crate::primerdlg::spec_from_annotation(&doc.primers()[i], &template);
        spec.name = self.unique_primer_name(&format!("{} copy", spec.name));
        let name = spec.name.clone();
        self.primer_dlg_apply(PrimerEdit::Add(Box::new(spec)));
        self.status = format!("Duplicated primer as {name}");
    }

    /// `Primers ▸ Remove Primer`。
    pub fn remove_selected_primer(&mut self) {
        let Some(i) = self.active_primer_index() else {
            self.status = "No primer selected — pick one in the Primers tab".into();
            return;
        };
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let name = doc.primers()[i].name.clone();
        self.selected_primer = None;
        self.primer_dlg_apply(PrimerEdit::Remove { name: name.clone() });
        self.status = format!("Removed primer {name}");
    }

    /// 引物编辑落地。
    fn primer_dlg_apply(&mut self, edit: PrimerEdit) {
        let label = edit.label();
        let is_add = matches!(edit, PrimerEdit::Add(_));
        let Some(doc) = self.doc.as_mut() else {
            return;
        };
        match doc.apply_primer(edit) {
            Ok(()) => {
                let n = doc.primers().len();
                if is_add {
                    self.selected_primer = Some(n.saturating_sub(1));
                }
                self.status = format!("{label} — {n} primers");
            }
            Err(e) => self.status = format!("Primer edit failed: {e}"),
        }
    }

    /// `Primers ▸ Analyze Selected Primers…`：Tm / GC / 发夹 / 二聚体。
    pub fn analyze_selected_primer(&mut self) {
        let Some(i) = self.active_primer_index() else {
            self.status = "No primer selected — pick one in the Primers tab".into();
            return;
        };
        let (template, circular) = self.primer_template();
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let p = doc.primers()[i].clone();
        let others: Vec<(String, String)> = doc
            .primers()
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, q)| (q.name.clone(), q.sequence.clone()))
            .collect();
        let sites: Vec<primer::Anneal> =
            primer::anneal_bindings(&template, p.sequence.as_bytes(), circular, 10);
        let annealed = match sites.iter().max_by_key(|a| a.len) {
            Some(a) if a.len <= p.sequence.len() => {
                p.sequence[p.sequence.len() - a.len..].to_string()
            }
            _ => p.sequence.clone(),
        };
        let mut rows: Vec<(String, String)> = Vec::new();
        rows.push(("Name".into(), p.name.clone()));
        rows.push(("Sequence".into(), p.sequence.clone()));
        rows.push((
            "Length".into(),
            format!("{}-mer", p.sequence.len().max(annealed.len())),
        ));
        rows.push((
            "GC content".into(),
            format!("{:.0}%", primer::gc_percent(p.sequence.as_bytes())),
        ));
        rows.push((
            format!("Tm ({}-mer annealed)", annealed.len()),
            format!("{:.1} °C", primer::tm(annealed.as_bytes())),
        ));
        rows.push((
            "Tm (full primer)".into(),
            format!("{:.1} °C", primer::tm(p.sequence.as_bytes())),
        ));
        if let Some(s) = sites.iter().max_by_key(|a| a.len) {
            rows.push((
                "Annealing site".into(),
                format!(
                    "{}..{} ({})",
                    s.start + 1,
                    s.end.max(s.start + s.len),
                    if s.strand == 0 { "forward" } else { "reverse" }
                ),
            ));
        } else {
            rows.push(("Annealing site".into(), "none".into()));
        }
        match primer::hairpin(p.sequence.as_bytes()) {
            Some(h) => rows.push((
                "Hairpin".into(),
                format!("{}-bp stem, {}-nt loop, Tm {:.1} °C", h.stem, h.loop_len, h.tm),
            )),
            None => rows.push(("Hairpin".into(), "none detected".into())),
        }
        match primer::self_dimer(p.sequence.as_bytes()) {
            Some(d) => rows.push((
                "Self-dimer".into(),
                format!(
                    "{}-bp overlap, Tm {:.1} °C{}",
                    d.overlap,
                    d.tm,
                    if d.three_prime { ", 3′ end" } else { "" }
                ),
            )),
            None => rows.push(("Self-dimer".into(), "none detected".into())),
        }
        let mut cross: Vec<(String, primer::Dimer)> = others
            .iter()
            .filter_map(|(name, seq)| {
                primer::cross_dimer(p.sequence.as_bytes(), seq.as_bytes())
                    .map(|d| (name.clone(), d))
            })
            .collect();
        cross.sort_by(|a, b| {
            b.1.tm
                .partial_cmp(&a.1.tm)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (name, d) in cross.into_iter().take(3) {
            rows.push((
                format!("Dimer with {name}"),
                format!("{}-bp overlap, Tm {:.1} °C", d.overlap, d.tm),
            ));
        }
        self.primer_analysis = Some(rows);
    }

    /// `Actions ▸ PCR…`：引物对设计（长度/Tm/GC/产物长度约束）。
    pub fn open_design_dialog(&mut self) {
        match self.doc.as_ref() {
            None => self.status = "No document open".into(),
            Some(d) if !d.is_dna() => {
                self.status = "Primer design needs a DNA sequence".into()
            }
            Some(d) => {
                let len = d.len();
                self.design_dlg = Some(crate::primerdlg::DesignDlg::new(len, self.selection));
            }
        }
    }

    /// 把一对设计好的引物加入文档。
    fn add_primer_pair(&mut self, pair: &primer::PrimerPair) {
        let Some(len) = self.doc.as_ref().map(|d| d.len()) else {
            return;
        };
        let f_name = self.unique_primer_name("Fwd");
        let r_name = self.unique_primer_name("Rev");
        let (f, r) = crate::primerdlg::pair_specs(pair, len, &f_name, &r_name);
        self.primer_dlg_apply(PrimerEdit::Add(Box::new(f)));
        self.primer_dlg_apply(PrimerEdit::Add(Box::new(r)));
        self.selected_primer = None;
        self.status = format!(
            "Added {f_name} / {r_name} — {} bp product",
            pair.product_len()
        );
    }

    /// 「Add/Edit Primer」窗口。
    fn primer_dlg_window(&mut self, ctx: &egui::Context) {
        if self.primer_dlg.is_none() {
            return;
        }
        let (template, circular) = self.primer_template();
        let Some(dlg) = self.primer_dlg.as_mut() else {
            return;
        };
        match crate::primerdlg::draw(dlg, ctx, &template, circular) {
            Some(crate::primerdlg::Outcome::Apply(e)) => {
                self.primer_dlg = None;
                self.primer_dlg_apply(*e);
            }
            Some(crate::primerdlg::Outcome::Cancel) => {
                self.primer_dlg = None;
            }
            None => {}
        }
    }

    /// 「Design Primers」窗口。
    fn design_dlg_window(&mut self, ctx: &egui::Context) {
        if self.design_dlg.is_none() {
            return;
        }
        let (template, circular) = self.primer_template();
        let Some(dlg) = self.design_dlg.as_mut() else {
            return;
        };
        match crate::primerdlg::draw_design(dlg, ctx, &template, circular) {
            Some(crate::primerdlg::DesignOutcome::Add(pair)) => {
                self.design_dlg = None;
                self.add_primer_pair(&pair);
            }
            Some(crate::primerdlg::DesignOutcome::Cancel) => {
                self.design_dlg = None;
            }
            None => {}
        }
    }

    /// 「Analyze Selected Primers…」结果窗口。
    fn primer_analysis_window(&mut self, ctx: &egui::Context) {
        let Some(rows) = self.primer_analysis.as_ref() else {
            return;
        };
        let mut open = true;
        egui::Window::new(crate::i18n::text(ctx,"Primer Analysis"))
            .id(egui::Id::unique("primer_analysis"))
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                egui::Grid::new("primer_analysis_grid")
                    .num_columns(2)
                    .spacing([14.0, 4.0])
                    .show(ui, |ui| {
                        for (k, v) in rows {
                            ui.label(egui::RichText::new(crate::i18n::message(ui.ctx(),&k)).color(theme::TEXT_DIM).size(11.5));
                            let value=if matches!(k.as_str(),"Name"|"Sequence"){v.clone()}else{crate::i18n::message(ui.ctx(),v)};
                            ui.label(egui::RichText::new(value).monospace().size(11.5));
                            ui.end_row();
                        }
                    });
            });
        if !open {
            self.primer_analysis = None;
        }
    }

    /// 「Feature Color…」的小调色板窗口。
    fn feature_color_window(&mut self, ctx: &egui::Context) {
        if !self.feature_color_popup {
            return;
        }
        let mut open = true;
        let mut picked: Option<String> = None;
        egui::Window::new(crate::i18n::text(ctx,"Feature Color"))
            .id(egui::Id::unique("feature_color"))
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.label(
                    egui::RichText::new(crate::i18n::text(ui.ctx(),"Color for the selected feature"))
                        .size(11.0)
                        .color(theme::TEXT_DIM),
                );
                for chunk in hx_formats::feature_edit::PALETTE.chunks(6) {
                    ui.horizontal(|ui| {
                        for hex in chunk {
                            let c = theme::hex(hex, egui::Color32::GRAY);
                            let (rect, resp) =
                                ui.allocate_exact_size(egui::vec2(24.0, 20.0), egui::Sense::click());
                            let sw = egui::Rect::from_center_size(
                                rect.center(),
                                egui::vec2(22.0, 18.0),
                            );
                            ui.painter().rect_filled(sw, 2.0, c);
                            ui.painter().rect_stroke(
                                sw,
                                2.0,
                                egui::Stroke::new(1.0, theme::BORDER),
                                egui::StrokeKind::Inside,
                            );
                            if resp.clicked() {
                                picked = Some((*hex).to_string());
                            }
                        }
                    });
                }
                if ui.button(crate::i18n::text(ui.ctx(),"No Color")).clicked() {
                    picked = Some(hx_formats::feature_edit::NO_COLOR.to_string());
                }
            });
        if let Some(c) = picked {
            self.set_selected_feature_color(&c);
        }
        if !open {
            self.feature_color_popup = false;
        }
    }

    /// 激活某个标签页（保存当前标签页状态 → 恢复目标标签页状态）。
    pub fn activate_tab(&mut self, idx: usize) {
        if self.tabs.is_empty() || idx >= self.tabs.len() || idx == self.active {
            return;
        }
        let v = DocViewState::take(self);
        if let Some(d) = self.doc.take() {
            self.tabs[self.active].doc = Some(d);
        }
        if let Some(vw) = self.viewer.take() {
            self.tabs[self.active].viewer = Some(vw);
        }
        self.tabs[self.active].view = v;
        self.active = idx;
        self.doc = self.tabs[idx].doc.take();
        self.viewer = self.tabs[idx].viewer.take();
        let v = std::mem::take(&mut self.tabs[idx].view);
        v.apply(self);
        self.sync_path_fields();
        if let Some(vw) = &self.viewer {
            self.status = format!("Switched to {}", vw.title());
        } else if let Some(d) = &self.doc {
            let (name, n, feats) = (d.display_name(), d.len(), d.features().len());
            // 原生只读文档：状态栏左段已经写了类别与字节数，这里不再重复（否则挤爆状态栏）
            self.status = if d.native_class().is_some() {
                format!("Switched to {name}")
            } else {
                format!("Switched to {name} ({n} {unit}, {feats} features)", unit = d.unit())
            };
        }
    }

    /// 下一个 / 上一个标签页（`delta` = ±1；Ctrl+Tab / Ctrl+Shift+Tab）。
    pub fn next_tab(&mut self, delta: isize) {
        let n = self.tabs.len();
        if n < 2 {
            return;
        }
        let i = (self.active as isize + delta).rem_euclid(n as isize) as usize;
        self.activate_tab(i);
    }

    /// 拖动重排页签。`insert` 是**插入边界**（0..=tabs.len()，位于该下标页签之前）。
    /// 活动页签的文档在 `self.doc` / `self.viewer` 里，移动后只需修正 `active`。
    pub fn move_tab(&mut self, from: usize, insert: usize) {
        let n = self.tabs.len();
        if from >= n {
            return;
        }
        let mut to = insert.min(n);
        if to > from {
            to -= 1;
        }
        if to == from {
            return;
        }
        let item = self.tabs.remove(from);
        self.tabs.insert(to, item);
        self.active = if self.active == from {
            to
        } else if from < self.active && to >= self.active {
            self.active - 1
        } else if from > self.active && to <= self.active {
            self.active + 1
        } else {
            self.active
        };
        self.status = format!("Moved tab to position {}", to + 1);
    }

    /// 关闭某个标签页。
    pub fn close_tab(&mut self, idx: usize) {
        if idx >= self.tabs.len() {
            return;
        }
        let was_active = idx == self.active;
        let closed = if was_active {
            self.doc.take()
        } else {
            self.tabs[idx].doc.take()
        };
        let closed_viewer = if was_active {
            self.viewer.take()
        } else {
            self.tabs[idx].viewer.take()
        };
        let name = closed
            .as_ref()
            .map(|d| d.display_name())
            .or_else(|| closed_viewer.as_ref().map(|v| v.title()))
            .unwrap_or_else(|| "document".into());
        let dirty = closed.as_ref().is_some_and(|d| d.is_dirty())
            || closed_viewer.as_ref().is_some_and(|v| v.is_dirty());
        drop(closed);
        drop(closed_viewer);
        self.tabs.remove(idx);
        self.status = if dirty {
            format!("Closed {name} (unsaved edits discarded)")
        } else {
            format!("Closed {name}")
        };
        if self.tabs.is_empty() {
            self.active = 0;
            self.doc = None;
            DocViewState::default().apply(self);
            return;
        }
        if was_active {
            let next = idx.min(self.tabs.len() - 1);
            self.active = next;
            self.doc = self.tabs[next].doc.take();
            self.viewer = self.tabs[next].viewer.take();
            let v = std::mem::take(&mut self.tabs[next].view);
            v.apply(self);
            self.sync_path_fields();
        } else if idx < self.active {
            self.active -= 1;
        }
    }

    /// 让 `save_path` / `open_path` 跟随活动文档（切标签页后调用）。
    ///
    /// 侧栏根目录**不跟随**（用户裁定 2026-10-18）；只展开到该文档，让它在树里可见。
    fn sync_path_fields(&mut self) {
        // 切换文档时关闭内部 Save As，避免把新标签页写到上一个文档的路径。
        self.show_save_dialog = false;
        if let Some(p) = self.viewer.as_ref().map(|v| v.path().to_path_buf()) {
            // 只读文件没有保存路径的概念，但侧栏/路径框应指向它。
            self.save_path = p.display().to_string();
            self.open_path = p.display().to_string();
            self.expand_ancestors(&p);
            return;
        }
        if let Some(p) = self.doc.as_ref().and_then(|d| d.path.clone()) {
            self.save_path = p.display().to_string();
            self.open_path = p.display().to_string();
            self.expand_ancestors(&p);
        } else if let Some(p) = self.default_save_path() {
            // 未保存文档（克隆产物、导入记录、新建序列）也必须有自己的落点。
            // 否则切到另一个未保存标签页后，Ctrl+S / Save As 会沿用上一个
            // 文档的路径，甚至覆盖源文件。
            self.save_path = p.display().to_string();
        }
    }

    /// 未保存文档的默认保存路径：侧栏根目录 + 文档名（补上分子类型扩展名）。
    fn default_save_path(&self) -> Option<PathBuf> {
        let doc = self.doc.as_ref()?;
        let dir = self
            .folder
            .clone()
            .unwrap_or_else(|| std::env::current_dir().unwrap_or_default());
        let mut name = doc.display_name();
        let has_molecule_ext = std::path::Path::new(&name)
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| matches!(e.to_ascii_lowercase().as_str(), "dna" | "rna" | "prot"));
        if !has_molecule_ext {
            name.push_str(match doc.kind() {
                hx_core::MoleculeKind::Dna => ".dna",
                hx_core::MoleculeKind::Rna => ".rna",
                hx_core::MoleculeKind::Protein => ".prot",
            });
        }
        Some(dir.join(name))
    }

    /// 关闭当前文档（SnapGene 的 File → Close Tab）。
    pub fn close_doc(&mut self) {
        if self.tabs.is_empty() {
            return;
        }
        self.close_tab(self.active);
    }

    pub fn recompute_sites(&mut self) {
        let Some(doc) = self.doc.take() else {
            return;
        };
        let mut hits = std::mem::take(&mut self.sites);
        if doc.is_dna() {
            hits.recompute(doc.sequence(), doc.circular());
        } else {
            hits.clear_hits();
        }
        self.sites = hits;
        self.doc = Some(doc);
        self.sites_dirty = false;
    }

    // ---- 酶集合 ----

    /// 设定显示集合（预设）；命中数会随后由 `sites_dirty` 触发重算。
    pub fn apply_preset(&mut self, preset: EnzymePreset) {
        self.enzyme_preset = preset;
        use hx_enzymes::sets::{Parent, Filter};
        match preset {
            EnzymePreset::Unique6Nonredundant => { self.enzyme_parent = Parent::Nonredundant; self.enzyme_filter = Filter::UniqueSixPlus; }
            EnzymePreset::Unique6 => { self.enzyme_parent = Parent::All; self.enzyme_filter = Filter::UniqueSixPlus; }
            EnzymePreset::Nonredundant => { self.enzyme_parent = Parent::Nonredundant; self.enzyme_filter = Filter::All; }
            EnzymePreset::SixPlus => { self.enzyme_parent = Parent::All; self.enzyme_filter = Filter::SixPlus; }
            _ => {}
        }
        let Some(doc) = self.doc.as_ref() else {
            self.sites.selected.clear();
            self.sites_dirty = true;
            return;
        };
        if !doc.is_dna() {
            self.sites.selected.clear();
            self.sites_dirty = true;
            return;
        }
        let seq = doc.sequence().to_vec();
        let circular = doc.circular();
        self.sites.selected = match preset {
            EnzymePreset::None => Vec::new(),
            EnzymePreset::Common => indices_of(COMMON_ENZYMES),
            EnzymePreset::All6 => hx_enzymes::six_cutter_indices(),
            EnzymePreset::SixPlus | EnzymePreset::Unique6 | EnzymePreset::Nonredundant | EnzymePreset::Unique6Nonredundant => hx_enzymes::sets::select_ordered(self.enzyme_parent, self.enzyme_filter, &seq, circular, &self.preferences.suppliers),
            EnzymePreset::Custom => return,
        };
        self.sites_dirty = true;
    }

    /// 让图谱在给定尺寸下「适配窗口」（线性：整条序列铺满；环状：默认半径）。
    pub fn fit(&mut self) {
        let (Some(doc), Some(size)) = (self.doc.as_ref(), self.last_size) else {
            return;
        };
        let len = doc.len().max(1) as f32;
        if self.map_layout.is_ring(doc.circular(), doc.len()) {
            self.zoom = 1.0;
        } else {
            self.zoom = ((size.x - 90.0) / len).max(0.02);
            self.pan_bp = 0.0;
        }
        self.seq_zoom = 1.0;
    }

    pub fn zoom_in(&mut self) {
        self.zoom_by(1.25);
    }

    pub fn zoom_out(&mut self) {
        self.zoom_by(1.0 / 1.25);
    }

    pub fn zoom_by(&mut self, f: f32) {
        if self.tab == Tab::Sequence {
            self.seq_zoom = (self.seq_zoom * f).clamp(0.4, 3.0);
        } else {
            self.zoom = (self.zoom * f).clamp(0.05, 40.0);
        }
    }

    /// 程序化调整环状图谱的显示原点。普通拖拽路径不会调用此方法。
    pub fn rotate_by(&mut self, d: f32) {
        self.rotation = hx_render::norm_angle(self.rotation + d);
    }

    /// 视图复位（缩放/显示原点/平移/滚动）。
    pub fn reset_view(&mut self) {
        self.zoom = 1.0;
        self.seq_zoom = 1.0;
        self.rotation = 0.0;
        self.pan_bp = 0.0;
        self.seq.scroll_y = 0.0;
        self.status = "View reset".into();
    }

    /// 把某个特征带到视野中央（环状：程序化调整显示原点；线性：平移到可见；序列：滚到该行）。
    pub fn center_on_feature(&mut self, idx: usize) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let Some(f) = doc.features().get(idx) else {
            return;
        };
        let len = doc.len().max(1);
        let Some(seg) = f.segments.first() else {
            return;
        };
        let mid = ((seg.range.start + seg.range.end) / 2).min(len - 1);
        self.selected_feature = Some(idx);
        if doc.circular() {
            let frac = mid as f32 / len as f32;
            self.rotation = hx_render::norm_angle(-frac * std::f32::consts::TAU);
        } else {
            if let Some(size) = self.last_size {
                self.zoom = ((size.x - 90.0) / len as f32)
                    .clamp(0.05, 40.0)
                    .max(self.zoom);
                let visible_bp = (size.x - 90.0) / self.zoom.max(0.001) as f32;
                self.pan_bp = (mid as f32 - visible_bp * 0.5).max(0.0) as f64;
            }
        }
        // 序列视图：滚到该行
        let line_h = self.seq.line_h.max(1.0);
        let row = mid / self.seq.line_bp.max(1);
        self.seq.scroll_y = (row as f32 * line_h - line_h * 1.5).max(0.0);
        self.status = format!(
            "Centered on {} ({}..{})",
            f.name,
            nfmt(seg.range.start + 1),
            nfmt(seg.range.end)
        );
    }

    /// 全选序列。
    pub fn select_all(&mut self) {
        let Some(len) = self.doc.as_ref().map(|d| d.len()) else {
            return;
        };
        if len > 0 {
            self.selection = Some((0, len));
            self.selected_residue = None;
            self.status = format!("Selected whole sequence ({} bp)", nfmt(len));
        }
    }

    /// 添加酶（按名称；未知名返回 false）。
    pub fn add_enzyme(&mut self, name: &str) -> bool {
        let Some(e) = hx_enzymes::by_name(name.trim()) else {
            return false;
        };
        let Some(i) = hx_enzymes::ENZYMES.iter().position(|x| std::ptr::eq(x, e)) else {
            return false;
        };
        if self.sites.selected.contains(&i) {
            return true;
        }
        self.sites.selected.push(i);
        self.sites.selected.sort_unstable();
        self.sites_dirty = true;
        self.enzyme_preset = EnzymePreset::Custom;
        true
    }

    pub fn remove_enzyme(&mut self, idx: usize) {
        self.sites.selected.retain(|&i| i != idx);
        self.sites_dirty = true;
        self.enzyme_preset = EnzymePreset::Custom;
    }

    // ---- 编辑 ----

    /// 可读取/复制的区间，包括跨原点选段。
    pub fn selected_range(&self) -> Option<(usize, usize)> {
        let doc = self.doc.as_ref()?;
        let len = doc.len();
        if len == 0 {
            return None;
        }
        if let Some((s, e)) = self.selection {
            return (s < len && e <= len && s != e && (s < e || doc.circular()))
                .then_some((s, e));
        }
        let c = self.cursor_bp?;
        (c < len).then(|| (c, (c + 1).min(len)))
    }

    /// 现有编辑引擎接受的连续区间；跨原点选段不能退回旧光标处编辑。
    pub fn effective_range(&self) -> Option<(usize, usize)> {
        self.selected_range().filter(|(s, e)| s < e)
    }

    fn apply_edit(&mut self, edit: Edit, label: String, cursor_after: Option<usize>) {
        let Some(doc) = self.doc.as_mut() else {
            return;
        };
        match doc.apply(edit) {
            Ok(()) => {
                let len = doc.len();
                let unit = doc.unit().to_string();
                self.selection = None;
                self.selected_residue = None;
                // 光标可以停在 `len`（序列末尾之后）= 追加位置
                self.cursor_bp = match cursor_after {
                    Some(c) => Some(c.min(len)),
                    None => self.cursor_bp.map(|c| c.min(len)),
                };
                self.sites_dirty = true;
                self.status = format!("{label} — now {} {unit}", nfmt(len));
            }
            Err(e) => self.status = format!("Edit failed: {e}"),
        }
    }

    /// 历史/状态文案里带单位的数量：`3 bases` / `1 base` / `2 residues`。
    ///
    /// （`unit()` 给的是 `bp` / `bases` / `aa`，写进句子里不成话，所以单独造词。）
    fn unit_count(n: usize, protein: bool) -> String {
        let (one, many) = if protein {
            ("residue", "residues")
        } else {
            ("base", "bases")
        };
        format!("{n} {}", if n == 1 { one } else { many })
    }

    /// 序列视图里的**直接输入**（SnapGene 官方说明：单击放置光标 = 插入，
    /// 选中一段 = 替换）。
    ///
    /// 只接受字母：核酸只收 IUPAC 碱基（大小写照输入保留），蛋白收 A–Z。
    /// 位置：有选区 → 替换选区；只有光标 → 在该碱基**之前**插入（与绘制一致）。
    pub fn type_residues(&mut self, text: &str) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        if doc.native_class().is_some() {
            return;
        }
        let protein = doc.kind() == hx_core::MoleculeKind::Protein;
        let len = doc.len();
        let letters: Vec<u8> = text
            .bytes()
            .filter(|b| {
                if protein {
                    b.is_ascii_alphabetic()
                } else {
                    is_iupac_base(*b)
                }
            })
            .collect();
        if letters.is_empty() || len == 0 {
            return;
        }
        if let Some((s, e)) = self.selection.filter(|(s, e)| s < e) {
            let label = format!(
                "Replaced {}..{} with {}",
                nfmt(s + 1),
                nfmt(e),
                Self::unit_count(letters.len(), protein)
            );
            let after = s + letters.len();
            self.apply_edit(
                Edit::Replace {
                    start: s,
                    end: e,
                    seq: letters,
                },
                label,
                Some(after),
            );
        } else if let Some(c) = self.cursor_bp {
            let at = c.min(len);
            let label = format!(
                "Inserted {} at {}",
                Self::unit_count(letters.len(), protein),
                nfmt(at + 1)
            );
            let after = at + letters.len();
            self.apply_edit(Edit::Insert { at, seq: letters }, label, Some(after));
        }
    }

    /// Backspace：删掉光标**前**一个碱基（SnapGene 的退格），有选区则删选区。
    pub fn do_backspace(&mut self) {
        if let Some((s, e)) = self.selection.filter(|(s, e)| s < e) {
            let label = format!("Deleted {}..{}", nfmt(s + 1), nfmt(e));
            self.apply_edit(Edit::Delete { start: s, end: e }, label, Some(s));
            return;
        }
        match self.cursor_bp {
            Some(c) if c > 0 => {
                let at = c - 1;
                let protein = self
                    .doc
                    .as_ref()
                    .is_some_and(|d| d.kind() == hx_core::MoleculeKind::Protein);
                let label = format!(
                    "Deleted {} at {}",
                    if protein { "residue" } else { "base" },
                    nfmt(c)
                );
                self.apply_edit(
                    Edit::Delete {
                        start: at,
                        end: c,
                    },
                    label,
                    Some(at),
                );
            }
            Some(_) => self.status = "Already at the start of the sequence".into(),
            None => self.status = "Place the cursor in the sequence first".into(),
        }
    }

    /// 光标左右移动（`Shift` = 扩选）。
    pub fn move_caret(&mut self, delta: isize, extend: bool) {
        let Some(len) = self.doc.as_ref().map(|d| d.len()) else {
            return;
        };
        let cur = match (self.selection, self.cursor_bp) {
            (Some((s, e)), _) if extend => Some(if delta > 0 { e } else { s }),
            _ => self.cursor_bp.or(Some(0)),
        };
        let Some(cur) = cur else { return };
        let next = (cur as isize + delta).clamp(0, len as isize) as usize;
        if extend {
            let anchor = match self.selection {
                Some((s, e)) => {
                    if delta > 0 {
                        s
                    } else {
                        e
                    }
                }
                None => cur,
            };
            let (a, b) = if anchor <= next {
                (anchor, next)
            } else {
                (next, anchor)
            };
            self.selection = (a < b).then_some((a, b));
        } else {
            self.selection = None;
        }
        self.selected_residue = None;
        self.cursor_bp = Some(next.min(len));
    }

    /// 光标跳到序列首/尾。
    pub fn move_caret_home(&mut self, end: bool, extend: bool) {
        let Some(len) = self.doc.as_ref().map(|d| d.len()) else {
            return;
        };
        let target = if end { len } else { 0 };
        if extend {
            let anchor = match self.selection {
                Some((s, e)) => {
                    if end {
                        s
                    } else {
                        e
                    }
                }
                None => self.cursor_bp.unwrap_or(0),
            };
            let (a, b) = if anchor <= target {
                (anchor, target)
            } else {
                (target, anchor)
            };
            self.selection = (a < b).then_some((a, b));
        } else {
            self.selection = None;
        }
        self.selected_residue = None;
        self.cursor_bp = Some(target.min(len));
    }

    pub fn do_delete(&mut self) {
        let range = self.selection.filter(|(s, e)| s < e).or_else(|| {
            self.cursor_bp.filter(|&c| self.doc.as_ref().is_some_and(|d| c < d.len())).map(|c| (c, c + 1))
        });
        self.pending_delete = range;
    }

    fn delete_immediately(&mut self) {
        if let Some((s, e)) = self.selection.filter(|(s, e)| s < e) {
            let label = format!("Deleted {}..{}", nfmt(s + 1), nfmt(e));
            self.apply_edit(Edit::Delete { start: s, end: e }, label, Some(s));
            return;
        }
        let Some(c) = self.cursor_bp else {
            self.status = "Nothing selected — click or drag in the sequence first".into();
            return;
        };
        let len = self.doc.as_ref().map(|d| d.len()).unwrap_or(0);
        if c < len {
            let label = format!("Deleted base {}", nfmt(c + 1));
            self.apply_edit(
                Edit::Delete {
                    start: c,
                    end: c + 1,
                },
                label,
                Some(c),
            );
        }
    }

    /// 选段的碱基文本（复制/剪切用）。
    pub fn selected_bases(&self) -> Option<String> {
        let (s, e) = self.selected_range()?;
        let doc = self.doc.as_ref()?;
        let seq = doc.sequence();
        if s < e {
            Some(String::from_utf8_lossy(seq.get(s..e)?).to_string())
        } else {
            let mut bases = seq.get(s..)?.to_vec();
            bases.extend_from_slice(seq.get(..e)?);
            Some(String::from_utf8_lossy(&bases).to_string())
        }
    }

    pub fn do_copy(&mut self) {
        self.do_copy_form(CopyForm::TopStrand);
    }

    /// 选段按指定形式取文本（不落剪贴板）。
    ///
    /// SnapGene 的 `Copy` 下拉有五种取法（oracle `:99` 实测，菜单与剪贴板逐项核对）：
    /// 上链、下链 5′→3′（反向互补）、下链 3′→5′（互补）、氨基酸 1 字母 / 3 字母。
    pub fn copy_text_for(&self, form: CopyForm) -> Option<String> {
        let (s, e) = self.selected_range()?;
        let doc = self.doc.as_ref()?;
        let seq = doc.sequence();
        if seq.is_empty() {
            return None;
        }
        // 跨原点选段（环状）按「从 s 绕回 e」拼成连续片段
        let mut sub: Vec<u8> = Vec::with_capacity(e.saturating_sub(s));
        if s < e {
            sub.extend_from_slice(seq.get(s..e)?);
        } else {
            sub.extend_from_slice(seq.get(s..)?);
            sub.extend_from_slice(seq.get(..e)?);
        }
        let text = match form {
            CopyForm::TopStrand => String::from_utf8_lossy(&sub).to_string(),
            CopyForm::BottomStrand5to3 => {
                String::from_utf8_lossy(&hx_algo::reverse_complement(&sub)).to_string()
            }
            CopyForm::BottomStrand3to5 => {
                let comp: Vec<u8> = sub.iter().map(|b| hx_algo::complement(*b)).collect();
                String::from_utf8_lossy(&comp).to_string()
            }
            CopyForm::AminoAcids1 | CopyForm::AminoAcids3 => {
                let three = form == CopyForm::AminoAcids3;
                if doc.kind() == hx_core::MoleculeKind::Protein {
                    // 蛋白文件里选段本身就是残基：1 字母就是大写字母本身
                    if three {
                        let codes: Vec<&'static str> = sub
                            .iter()
                            .map(|b| hx_algo::three_letter(*b))
                            .collect();
                        join_codes(&codes, true)
                    } else {
                        sub.iter()
                            .map(|b| b.to_ascii_uppercase() as char)
                            .collect()
                    }
                } else {
                    // 核酸：从选段第一个碱基起按整密码子翻译（SnapGene 实测：
                    // 选 1..11 的 `TGGCGAATGGG` → 1 字母 `WRM`、3 字母 `Trp Arg Met`）
                    let up: Vec<u8> = sub
                        .iter()
                        .map(|b| match b.to_ascii_uppercase() {
                            b'U' => b'T',
                            c => c,
                        })
                        .collect();
                    let codes: Vec<&'static str> = up
                        .chunks_exact(3)
                        .map(hx_algo::codon_aa)
                        .map(|a| if three { a } else { hx_algo::one_letter(a) })
                        .collect();
                    join_codes(&codes, three)
                }
            }
        };
        Some(text)
    }

    /// SnapGene 的 `Copy ▸ …`：复制并把结果写到系统剪贴板。
    pub fn do_copy_form(&mut self, form: CopyForm) {
        match self.copy_text_for(form) {
            Some(s) if !s.is_empty() => {
                // 3 字母代码之间用空格分隔：残基个数 = 空格数 + 1
                let n = if form == CopyForm::AminoAcids3 {
                    s.chars().filter(|c| *c == ' ').count() + 1
                } else {
                    s.chars().count()
                };
                self.copy_text(s);
                self.status = if form.is_amino() {
                    format!(
                        "Copied {} {} to clipboard",
                        n,
                        if form == CopyForm::AminoAcids3 {
                            "amino acid(s), 3-letter codes"
                        } else {
                            "amino acid(s)"
                        }
                    )
                } else {
                    format!("Copied {} base(s) to clipboard", nfmt(n))
                };
            }
            Some(_) => self.status = "Nothing selected — drag in the sequence first".into(),
            None => self.status = "Nothing selected — drag in the sequence first".into(),
        }
    }

    pub fn do_cut(&mut self) {
        if self.effective_range().is_none() {
            self.status = "Select a contiguous range to cut".into();
            return;
        }
        match self.selected_bases() {
            Some(s) => {
                let n = s.len();
                self.copy_text(s);
                let (a, b) = self.effective_range().unwrap_or((0, 0));
                let label = format!("Cut {}..{}", nfmt(a + 1), nfmt(b));
                self.apply_edit(Edit::Delete { start: a, end: b }, label, Some(a));
                self.status = format!(
                    "Cut {} base(s) to clipboard — now {} bp",
                    nfmt(n),
                    nfmt(self.doc.as_ref().map(|d| d.len()).unwrap_or(0))
                );
            }
            None => self.status = "Nothing selected — click or drag in the sequence first".into(),
        }
    }

    /// 请求写入系统剪贴板（实际提交在 `ui()` 末，那里拿得到 `Context`）。
    fn copy_text(&mut self, text: String) {
        self.pending_copy = Some(text);
    }

    /// 剪贴板文本 → 插入（SnapGene 的 Edit → Paste）。
    pub fn do_paste_text(&mut self, text: &str) {
        if self.doc.is_none() {
            self.status = "No file open".into();
            return;
        }
        let seq = hx_algo::parse_seq_text(text);
        if seq.len() < 4 {
            self.paste = text.trim().to_string();
            self.status =
                "Clipboard does not look like a sequence — pasted into the edit box instead".into();
            return;
        }
        let at = self
            .effective_range()
            .map(|(s, _)| s)
            .or(self.cursor_bp)
            .unwrap_or(0);
        let label = format!("Pasted {} bp at {}", nfmt(seq.len()), nfmt(at + 1));
        let after = at + seq.len();
        self.apply_edit(Edit::Insert { at, seq }, label, Some(after));
    }

    /// 把剪贴板/文本缓冲插入光标处（Edit → Paste）。
    fn request_paste(&mut self,ctx:&egui::Context,reverse:bool) {
        self.pending_paste_reverse=Some(reverse);
        ctx.send_viewport_cmd(egui::ViewportCommand::RequestPaste);
    }

    fn paste_sequence(&mut self,text:&str,reverse:bool) {
        self.paste=text.to_string();
        if reverse {
            if self.is_protein(){self.status="Reverse-complement paste requires DNA or RNA".into();return;}
            let bases=hx_algo::parse_seq_text(text);
            let rc=hx_algo::reverse_complement(&bases);self.do_paste_text(&String::from_utf8_lossy(&rc));
        }else{self.do_paste_text(text);}
    }

    pub fn paste_buffer_into_seq(&mut self) {
        let t = self.paste.clone();
        if t.trim().is_empty() {
            self.status =
                "Clipboard is empty (Ctrl+V fills the buffer, or use Edit → Paste Sequence Text…)"
                    .into();
            return;
        }
        self.do_paste_text(&t);
    }

    /// 选区大小写转换（Edit → Make Uppercase / Lowercase）。
    pub fn do_change_case(&mut self, upper: bool) {
        let Some((s, e)) = self.effective_range() else {
            return;
        };
        let Some(doc) = self.doc.as_ref() else { return };
        let seq = doc.sequence();
        let mut sub: Vec<u8> = seq[s.min(seq.len())..e.min(seq.len())].to_vec();
        if sub.is_empty() {
            return;
        }
        for b in sub.iter_mut() {
            *b = if upper {
                b.to_ascii_uppercase()
            } else {
                b.to_ascii_lowercase()
            };
        }
        let text = String::from_utf8_lossy(&sub).to_string();
        let before = self.selection;
        self.selection = Some((s, e));
        self.selected_residue = None;
        self.do_replace_with(&text);
        if self.selection.is_some() {
            self.selection = before;
            self.selected_residue = None;
        }
        self.status = format!(
            "Selection {} (0-based {} .. {})",
            if upper { "uppercased" } else { "lowercased" },
            s,
            e
        );
    }

    /// 用文本替换当前选区（大小写转换等内部使用）。
    fn do_replace_with(&mut self, text: &str) {
        let Some((s, e)) = self.effective_range() else {
            return;
        };
        let seq = hx_algo::parse_seq_text(text);
        let label = format!("Replaced {}..{} ({} bp)", nfmt(s + 1), nfmt(e), seq.len());
        let after = s + seq.len();
        self.apply_edit(
            Edit::Replace {
                start: s,
                end: e,
                seq,
            },
            label,
            Some(after),
        );
    }

    pub fn do_insert(&mut self) {
        let Some(at) = self.effective_range().map(|(s, _)| s) else {
            self.status = "Click in the sequence to place the cursor".into();
            return;
        };
        let seq = hx_algo::parse_seq_text(&self.paste);
        if seq.is_empty() {
            self.status = "Nothing to insert — paste a sequence into the text field".into();
            return;
        }
        let label = format!("Inserted {} bp at {}", seq.len(), nfmt(at + 1));
        let after = at + seq.len();
        self.apply_edit(Edit::Insert { at, seq }, label, Some(after));
    }

    pub fn do_replace(&mut self) {
        let Some((s, e)) = self.effective_range() else {
            self.status = "Nothing selected — click or drag in the sequence first".into();
            return;
        };
        let seq = hx_algo::parse_seq_text(&self.paste);
        let label = if seq.is_empty() {
            format!(
                "Replaced {}..{} with nothing (deleted)",
                nfmt(s + 1),
                nfmt(e)
            )
        } else {
            format!(
                "Replaced {}..{} with {} bp",
                nfmt(s + 1),
                nfmt(e),
                seq.len()
            )
        };
        let after = s + seq.len();
        self.apply_edit(
            Edit::Replace {
                start: s,
                end: e,
                seq,
            },
            label,
            Some(after),
        );
    }

    pub fn do_revcomp(&mut self) {
        let Some((s, e)) = self.effective_range() else {
            self.status = "Nothing selected".into();
            return;
        };
        let slice = self
            .doc
            .as_ref()
            .and_then(|d| d.sequence().get(s..e).map(|x| x.to_vec()));
        let Some(slice) = slice else {
            return;
        };
        let rc = hx_algo::reverse_complement(&slice);
        let label = format!("Reverse-complemented {}..{}", nfmt(s + 1), nfmt(e));
        self.apply_edit(
            Edit::Replace {
                start: s,
                end: e,
                seq: rc,
            },
            label,
            Some(e),
        );
    }

    /// 把当前选中特征的第一个片段设为选段（便于整段删除/替换）。
    pub fn select_feature_range(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let Some(f) = self.selected_feature.and_then(|i| doc.features().get(i)) else {
            self.status = "No feature selected — pick one in the Features tab".into();
            return;
        };
        let len = doc.len();
        if let Some(seg) = f.segments.first() {
            if let Some((s, e)) = hx_render::split_wrap(seg.range.start, seg.range.end, len).first()
            {
                self.selection = Some((*s, *e));
                self.selected_residue = None;
                self.status = format!(
                    "Selected {} ({} bp){}",
                    f.name,
                    e - s,
                    if seg.range.start > seg.range.end {
                        " — wrapped feature: first segment only"
                    } else {
                        ""
                    }
                );
                return;
            }
        }
        self.status = "Feature has no segments".into();
    }

    pub fn do_undo(&mut self) {
        let msg = self.doc.as_mut().and_then(|d| d.undo());
        if let Some(l) = msg {
            self.selected_residue = None;
            self.sites_dirty = true;
            // 特征编辑撤销后，选中项/隐藏项可能已越界
            self.clamp_feature_selection();
            self.status = format!("Undo: {l}");
        }
    }

    pub fn do_redo(&mut self) {
        let msg = self.doc.as_mut().and_then(|d| d.redo());
        if let Some(l) = msg {
            self.selected_residue = None;
            self.sites_dirty = true;
            self.clamp_feature_selection();
            self.status = format!("Redo: {l}");
        }
    }

    /// 特征数变化后把视图层的特征下标（选中 / 隐藏 / 展开）收进合法范围。
    fn clamp_feature_selection(&mut self) {
        let n = self.doc.as_ref().map(|d| d.features().len()).unwrap_or(0);
        if self.selected_feature.is_some_and(|i| i >= n) {
            self.selected_feature = if n == 0 { None } else { Some(n - 1) };
        }
        self.selected_features.retain(|i| *i < n);
        self.hidden_features.retain(|i| *i < n);
        self.expanded_features.retain(|i| *i < n);
    }

    /// 保存（`save_as` = 强制另存到新路径）。
    pub fn do_save(&mut self, save_as: bool) {
        // 文本阅读页签：直接写回原文件（Ctrl+S）。Save As 走系统原生对话框。
        if self.doc.is_none() {
            let Some(v) = self.viewer.as_mut() else {
                return;
            };
            if !v.can_save() {
                self.status =
                    "This file is a read-only preview; there is nothing to save.".into();
                return;
            }
            let target = if save_as {
                let current = v.path().to_path_buf();
                match crate::sysopen::pick_save_path(crate::i18n::language_text(self.preferences.language,"Save As"), &current) {
                    Some(p) => Some(p),
                    None => {
                        self.status = "Save cancelled".into();
                        return;
                    }
                }
            } else {
                None
            };
            let title = v.title();
            let result = v.save_to(target.as_deref());
            match result {
                Ok(n) => {
                    self.sync_path_fields();
                    self.status = format!("Saved {title} ({n} bytes)");
                }
                Err(e) => self.status = format!("Save failed: {e}"),
            }
            return;
        }
        // Save As（或首次保存还没路径的文件）：先用**系统原生**「另存为」对话框
        // 选路径；系统上没有 zenity/kdialog 时打开内部 Save As 窗口。两条路都
        // 不会在用户未确认新路径时静默覆盖当前文件。
        let needs_path = save_as || self.doc.as_ref().is_some_and(|d| d.path.is_none());
        if needs_path {
            if crate::sysopen::has_native_dialog() {
                let current = self
                    .doc
                    .as_ref()
                    .and_then(|d| d.path.clone())
                    .or_else(|| {
                        let t = self.save_path.trim();
                        (!t.is_empty()).then(|| std::path::PathBuf::from(t))
                    })
                    .or_else(|| self.default_save_path());
                let default = current.unwrap_or_else(|| std::path::PathBuf::from("untitled.dna"));
                match crate::sysopen::pick_save_path(crate::i18n::language_text(self.preferences.language,"Save As"), &default) {
                    Some(p) => self.save_path = p.display().to_string(),
                    None => {
                        self.status = "Save cancelled".into();
                        return;
                    }
                }
            } else {
                self.begin_save_dialog();
                return;
            }
        }
        self.write_save_target(save_as);
    }

    /// 没有原生对话框时，打开内部 Save As 窗口并预填当前/默认路径。
    fn begin_save_dialog(&mut self) {
        let default = self
            .doc
            .as_ref()
            .and_then(|d| d.path.clone())
            .or_else(|| {
                let t = self.save_path.trim();
                (!t.is_empty()).then(|| PathBuf::from(t))
            })
            .or_else(|| self.default_save_path())
            .unwrap_or_else(|| PathBuf::from("untitled.dna"));
        self.save_as_path = default.display().to_string();
        self.show_save_dialog = true;
    }

    /// 把当前文档写到 `save_path`（`save_as` 或文档还没有路径时）。
    fn write_save_target(&mut self, save_as: bool) {
        let target: Option<PathBuf> = {
            let Some(doc) = self.doc.as_ref() else {
                return;
            };
            if save_as || doc.path.is_none() {
                let t = self.save_path.trim();
                if t.is_empty() {
                    self.begin_save_dialog();
                    self.status = "Type a path to save to (Save As)".into();
                    return;
                }
                Some(PathBuf::from(t))
            } else {
                doc.path.clone()
            }
        };
        let Some(path) = target else { return };
        let res = self.doc.as_mut().map(|d| d.save_to(&path));
        match res {
            Some(Ok((n, nodes))) => {
                self.save_path = path.display().to_string();
                self.status = if nodes > 0 {
                    format!(
                        "Saved {} ({n} bytes, history tree with {nodes} node(s))",
                        path.display()
                    )
                } else {
                    format!("Saved {} ({n} bytes, no changes)", path.display())
                };
            }
            Some(Err(e)) => self.status = format!("Save failed: {e}"),
            None => {}
        }
    }

    /// 导出当前视图为 SVG（同目录 `<名字>.map.svg` / `<名字>.seq.svg`）。
    ///
    /// 几何/文本生成在 `hx_algo::export` + `hx_render::svg`（与 `hx-cli export` 同路径）。
    pub fn do_export_svg(&mut self) {
        let target: Option<PathBuf> = {
            let Some(doc) = self.doc.as_ref() else {
                return;
            };
            let base = doc.path.clone().or_else(|| {
                let t = self.save_path.trim();
                if t.is_empty() {
                    None
                } else {
                    Some(PathBuf::from(t))
                }
            });
            let Some(base) = base else {
                self.status = "Export SVG: open or save a file first".into();
                return;
            };
            let seq_view = self.tab == Tab::Sequence;
            Some(base.with_extension(if seq_view { "seq.svg" } else { "map.svg" }))
        };
        let Some(out) = target else { return };
        let enzymes = self.sites.selected.clone();
        let seq_view = self.tab == Tab::Sequence;
        let ring = match self.doc.as_ref() {
            Some(doc) => self.map_layout.is_ring(doc.circular(), doc.len()),
            None => true,
        };
        let res = self.doc.as_ref().map(|doc| {
            if seq_view {
                doc.export_sequence_svg(&enzymes, &hx_algo::export::default_sequence_style())
            } else {
                doc.export_map_svg(&enzymes, &hx_algo::export::default_map_style(ring))
            }
        });
        match res {
            Some(ex) => match std::fs::write(&out, ex.svg.as_bytes()) {
                Ok(()) => {
                    self.status = format!(
                        "Exported {} ({} bytes, {} feature band(s), {} site(s))",
                        out.display(),
                        ex.svg.len(),
                        ex.bands,
                        ex.sites
                    );
                }
                Err(e) => self.status = format!("Export failed: {e}"),
            },
            None => {}
        }
    }

    /// 结构视图导出用的文件名（`<文档名>.structure.svg` / `.structure.txt`）。
    fn struct_export_path(&self, ext: &str) -> Option<PathBuf> {
        let doc = self.doc.as_ref()?;
        let base = doc.path.clone().or_else(|| {
            let t = self.save_path.trim();
            (!t.is_empty()).then(|| PathBuf::from(t))
        })?;
        Some(base.with_extension(ext))
    }

    /// 当前显示的结构（下标 + 配对表 + 布局）的只读快照。
    ///
    /// 返回 `(len, seq, pairs, svg)`：`svg` 需要 `hx_render::rna_svg`，这里一次算完，
    /// 免得调用方同时持有文档与布局两层借用。
    fn structure_snapshot(&mut self) -> Option<(usize, Vec<u8>, Vec<(usize, usize)>, String)> {
        let idx = self.struct_sel;
        let show_index = self.struct_show_index;
        let doc = self.doc.as_mut()?;
        let len = doc.len();
        let seq = doc.sequence().to_vec();
        let cache = match doc.rna_structure() {
            Ok(crate::doc::RnaCache::Ready(c)) => c,
            Ok(crate::doc::RnaCache::Folding) => {
                self.status = "Folding secondary structure…".into();
                return None;
            }
            Ok(crate::doc::RnaCache::None) => {
                self.status = "No secondary structure in this file".into();
                return None;
            }
            Err(e) => {
                self.status = format!("Secondary structure block: {e}");
                return None;
            }
        };
        cache.ensure_structures(len);
        let n = cache.structures.as_ref().map_or(0, |v| v.len());
        if n == 0 {
            return None;
        }
        let idx = idx.min(n - 1);
        cache.ensure_layout(idx, len);
        let (_, layout) = cache.layout.as_ref()?;
        let st = cache.structure.as_ref().ok()?;
        let pairs = cache.structures.as_ref()?.get(idx)?.0.clone();
        let partner = crate::view_struct::partner_map(len, &pairs);
        let colors = crate::view_struct::confidence_colors(st, len, &partner);
        let rgb: Vec<[u8; 3]> = colors.iter().map(|c| [c.r(), c.g(), c.b()]).collect();
        let svg = hx_render::rna_svg::to_svg(layout, &seq, &rgb, &pairs, show_index);
        Some((len, seq, pairs, svg))
    }

    /// `Secondary Structure ▸ Export as… ▸ Image…`：与屏幕同款的 SVG。
    pub fn do_export_structure_svg(&mut self) {
        let Some(out) = self.struct_export_path("structure.svg") else {
            self.status = "Export: open or save a file first".into();
            return;
        };
        let Some((len, _seq, pairs, svg)) = self.structure_snapshot() else {
            if self.status.is_empty() {
                self.status = "Nothing to export".into();
            }
            return;
        };
        match std::fs::write(&out, svg.as_bytes()) {
            Ok(()) => {
                self.status = format!(
                    "Exported {} ({} bytes, {len} bases, {} pair(s))",
                    out.display(),
                    svg.len(),
                    pairs.len()
                )
            }
            Err(e) => self.status = format!("Export failed: {e}"),
        }
    }

    /// 当前显示结构的 dot-bracket 文本（`Export as… ▸ Dot-Bracket…` 与复制按钮共用）。
    pub fn structure_dot_bracket(&mut self) -> Option<(String, usize, usize, f64)> {
        let idx = self.struct_sel;
        let doc = self.doc.as_mut()?;
        let len = doc.len();
        let cache = match doc.rna_structure() {
            Ok(crate::doc::RnaCache::Ready(c)) => c,
            _ => return None,
        };
        cache.ensure_structures(len);
        let st = cache.structure.as_ref().ok()?;
        let list = cache.structures.as_ref()?;
        let (pairs, energy) = list.get(idx.min(list.len().saturating_sub(1)))?.clone();
        Some((st.dot_bracket(&pairs, len), len, idx, energy))
    }

    pub fn do_export_dot_bracket(&mut self) {
        let Some(out) = self.struct_export_path("structure.txt") else {
            self.status = "Export: open or save a file first".into();
            return;
        };
        let Some((db, len, idx, energy)) = self.structure_dot_bracket() else {
            self.status = "No secondary structure to export".into();
            return;
        };
        let seq = self
            .doc
            .as_ref()
            .map(|d| String::from_utf8_lossy(d.sequence()).to_string())
            .unwrap_or_default();
        let name = self
            .doc
            .as_ref()
            .map(|d| d.display_name())
            .unwrap_or_else(|| "structure".into());
        let text = format!(
            ">{name} (structure {}, {energy:.2} kcal/mol)\n{seq}\n{db}\n",
            idx + 1
        );
        match std::fs::write(&out, text.as_bytes()) {
            Ok(()) => {
                self.status = format!(
                    "Exported {} ({len} bases, {} pair(s) in dot-bracket form)",
                    out.display(),
                    db.matches('(').count()
                )
            }
            Err(e) => self.status = format!("Export failed: {e}"),
        }
    }

    fn do_copy_dot_bracket(&mut self) {
        match self.structure_dot_bracket() {
            Some((db, len, _, _)) => {
                self.copy_text(db);
                self.status = format!("Copied dot-bracket structure ({len} bases)");
            }
            None => self.status = "No secondary structure to copy".into(),
        }
    }

    // ---- 界面：菜单 / 工具栏 / 页签 / 状态栏 / 工具条 ----

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        // Compact desktop menu spacing from the supplied expanded/collapsed views.
        ui.spacing_mut().item_spacing.x = 4.0;
        ui.spacing_mut().button_padding=egui::vec2(4.0,0.0);
        ui.style_mut().override_font_id=Some(egui::FontId::proportional(12.5));
        // 菜单结构对齐 SnapGene（`resources/translations/SnapGene_zh.qm` 抽取的
        // `*MenuState` 上下文）：条目名/顺序/分隔一致，未实现的动作置灰。
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button(crate::i18n::text(ui.ctx(),"File"), |ui| self.file_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Edit"), |ui| self.edit_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"View"), |ui| self.view_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Enzymes"), |ui| self.enzymes_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Features"), |ui| self.features_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Primers"), |ui| self.primers_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Actions"), |ui| self.actions_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Tools"), |ui| self.tools_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Window"), |ui| self.window_menu(ui));
            ui.menu_button(crate::i18n::text(ui.ctx(),"Help"), |ui| self.help_menu(ui));
        });
    }

    fn file_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        let has_viewer = self.viewer.is_some();
        let can_save = has_doc || self.viewer.as_ref().is_some_and(|v| v.can_save());
        let tabs = self.tabs.len();
        let any_dirty = self.any_dirty();
        let native = crate::sysopen::has_native_dialog();
        if item(ui, "New DNA or RNA File…", Some("Ctrl+N"), true) {
            self.open_new_file_dialog(hx_core::MoleculeKind::Dna);
            ui.close();
        }
        if item(ui, "New Protein File…", Some("Ctrl+Alt+N"), true) {
            self.open_new_file_dialog(hx_core::MoleculeKind::Protein);
            ui.close();
        }
        item(ui, "New File from Selection…", Some("Ctrl+Shift+N"), false);
        ui.separator();
        if item(ui, "Open Files…", Some("Ctrl+O"), true) {
            self.open_files_with_dialog();
            ui.close();
        }
        if item(ui, "Open Folder…", None, true) {
            self.open_folder_with_dialog();
            ui.close();
        }
        if native && item(ui, "Open File with System Dialog…", None, true) {
            self.open_file_system_dialog();
            ui.close();
        }
        ui.menu_button(crate::i18n::text(ui.ctx(),"Open Recent"), |ui| {
            let recent = self.recent.clone();
            if recent.is_empty() {
                let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"(empty)")));
            }
            for p in &recent {
                let label = p
                    .file_name()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| p.display().to_string());
                if ui
                    .button(label)
                    .on_hover_text(p.display().to_string())
                    .clicked()
                {
                    self.open(p);
                    ui.close();
                }
            }
            if !recent.is_empty() {
                ui.separator();
                if ui.button(crate::i18n::text(ui.ctx(),"Clear Menu")).clicked() {
                    self.recent.clear();
                    self.save_settings();
                    ui.close();
                }
            }
        });
        ui.separator();
        if item(ui, "Close Tab", Some("Ctrl+W"), tabs > 0) {
            self.close_doc();
            ui.close();
        }
        if item(ui, "Close All Tabs", None, tabs > 0) {
            self.close_all_tabs();
            ui.close();
        }
        if item(ui, "Close All Windows", Some("Ctrl+Shift+W"), tabs > 0) {
            self.close_all_tabs();
            ui.close();
        }
        ui.separator();
        if item(ui, "Save", Some("Ctrl+S"), can_save) {
            self.do_save(false);
            ui.close();
        }
        if item(ui, "Save As…", Some("Ctrl+Shift+S"), can_save) {
            self.do_save(true);
            ui.close();
        }
        if item(ui, "Save All", Some("Ctrl+Alt+S"), any_dirty) {
            self.save_all();
            ui.close();
        }
        ui.separator();
        if item(ui, "Show in File Manager", None, has_doc || has_viewer) {
            self.reveal_current();
            ui.close();
        }
        if item(ui, "Show Folder in File Manager", None, self.folder.is_some()) {
            if let Some(f) = self.folder.clone() {
                self.reveal_path(&f);
            }
            ui.close();
        }
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Export"), |ui| {
            if item(ui, "Map…", Some("Ctrl+E"), has_doc) {
                self.tab = Tab::Map;
                self.do_export_svg();
                ui.close();
            }
            if item(ui, "Sequence…", None, has_doc) {
                self.tab = Tab::Sequence;
                self.do_export_svg();
                ui.close();
            }
            item(ui, "History…", None, false);
            item(ui, "DNA…", None, false);
            item(ui, "Properties…", None, false);
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"Import"), |ui| {
            if item(ui, "Sequence File (FASTA / GenBank)…", None, true) {
                self.import_with_dialog();
                ui.close();
            }
            if item(ui, "FASTA File…", None, true) {
                self.import_with_dialog();
                ui.close();
            }
            if item(ui, "GenBank File…", None, true) {
                self.import_with_dialog();
                ui.close();
            }
            ui.separator();
            item(ui, "From Addgene Plasmids…", None, false);
            item(ui, "From NCBI Sequences…", None, false);
            item(ui, "From UniProt Sequences…", None, false);
            item(ui, "Vector NTI Database…", None, false);
            item(ui, "Clone Manager Primer Collection…", None, false);
            ui.separator();
            item(ui, "Import Primers…", None, false);
        });
        item(ui, "Batch Convert File Format…", None, false);
        item(ui, "Trim File History…", None, false);
        ui.separator();
        item(ui, "Page Setup…", Some("Ctrl+Shift+P"), false);
        item(ui, "Print…", Some("Ctrl+P"), false);
        ui.separator();
        if item(ui, "Exit", Some("Ctrl+Q"), true) {
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
        }
    }

    /// File ▸ Open Files…：有原生对话框就用它（可多选），否则用内置路径输入小窗口。
    pub fn open_files_with_dialog(&mut self) {
        if crate::sysopen::has_native_dialog() {
            self.open_files_system_dialog();
        } else {
            self.show_open_dialog = Some(OpenKind::File);
        }
    }

    /// 兼容旧名（单文件对话框）。
    pub fn open_file_with_dialog(&mut self) {
        self.open_files_with_dialog();
    }

    /// 系统原生对话框**多选**打开（SnapGene 的 Open Files… 能一次开多个）。
    pub fn open_files_system_dialog(&mut self) {
        let start = self
            .doc
            .as_ref()
            .and_then(|d| d.path.as_ref())
            .and_then(|p| p.parent().map(|p| p.to_path_buf()))
            .or_else(|| self.folder.clone());
        let picked = crate::sysopen::pick_files(crate::i18n::language_text(self.preferences.language,"Open Files"), start.as_deref());
        if picked.is_empty() {
            return;
        }
        for p in picked {
            self.open(&p);
        }
    }

    /// File ▸ Import ▸ …：系统原生对话框选一个序列文本文件并导入。
    pub fn import_with_dialog(&mut self) {
        let start = self.folder.clone().or_else(|| {
            self.doc
                .as_ref()
                .and_then(|d| d.path.as_ref())
                .and_then(|p| p.parent().map(|p| p.to_path_buf()))
        });
        let Some(p) = crate::sysopen::pick_file(crate::i18n::language_text(self.preferences.language,"Import Sequence File"), start.as_deref()) else {
            return;
        };
        match std::fs::read(&p) {
            Ok(bytes) => self.import_bytes(bytes, Some(p)),
            Err(e) => self.status = format!("Import failed: {e}"),
        }
    }

    /// 无条件是弹出**系统原生**文件对话框（`zenity`/`kdialog`）。
    pub fn open_file_system_dialog(&mut self) {
        let start = self
            .doc
            .as_ref()
            .and_then(|d| d.path.as_ref())
            .and_then(|p| p.parent().map(|p| p.to_path_buf()))
            .or_else(|| self.folder.clone());
        match crate::sysopen::pick_file(crate::i18n::language_text(self.preferences.language,"Open File"), start.as_deref()) {
            Some(p) => self.open(&p),
            None => {
                self.show_open_dialog = Some(OpenKind::File);
                self.status = "No system file dialog available — type a path instead".into();
            }
        }
    }

    /// File ▸ Open Folder…：同上。
    pub fn open_folder_with_dialog(&mut self) {
        if !crate::sysopen::has_native_dialog() {
            self.show_open_dialog = Some(OpenKind::Folder);
            return;
        }
        let start = self.folder.clone();
        match crate::sysopen::pick_folder(crate::i18n::language_text(self.preferences.language,"Open Folder"), start.as_deref()) {
            Some(p) => self.set_folder(Some(p)),
            None => {
                self.show_open_dialog = Some(OpenKind::Folder);
                self.status = "No system file dialog available — type a path instead".into();
            }
        }
    }

    /// 在系统文件管理器里显示活动文档（SnapGene 的 Show in Finder / Explorer）。
    pub fn reveal_current(&mut self) {
        let path = self
            .doc
            .as_ref()
            .and_then(|d| d.path.clone())
            .or_else(|| self.viewer.as_ref().map(|v| v.path().to_path_buf()))
            .or_else(|| self.folder.clone());
        if let Some(p) = path {
            self.reveal_path(&p);
        }
    }

    /// 在系统文件管理器里显示任意路径。
    pub fn reveal_path(&mut self, path: &std::path::Path) {
        match crate::sysopen::reveal(path) {
            Ok(()) => self.status = format!("Showed {} in file manager", path.display()),
            Err(e) => self.status = format!("File manager failed: {e}"),
        }
    }

    /// 是否有任何标签页有未保存改动（File ▸ Save All 的可用性）。
    pub fn any_dirty(&self) -> bool {
        self.viewer.as_ref().is_some_and(|v| v.is_dirty())
            || self.doc.as_ref().is_some_and(|d| d.is_dirty())
            || self
                .tabs
                .iter()
                .any(|t| {
                    t.doc.as_ref().is_some_and(|d| d.is_dirty())
                        || t.viewer.as_ref().is_some_and(|v| v.is_dirty())
                })
    }

    /// File ▸ Save All：把所有有改动的标签页写回各自路径。
    pub fn save_all(&mut self) {
        let mut saved = 0usize;
        let mut failed = 0usize;
        if let Some(v) = self.viewer.as_mut() {
            if v.is_dirty() {
                match v.save_to(None) {
                    Ok(_) => saved += 1,
                    Err(_) => failed += 1,
                }
            }
        }
        if let Some(d) = self.doc.as_mut() {
            if d.is_dirty() {
                match d.path.clone() {
                    Some(p) => match d.save_to(&p) {
                        Ok(_) => saved += 1,
                        Err(_) => failed += 1,
                    },
                    None => failed += 1,
                }
            }
        }
        for t in self.tabs.iter_mut() {
            if let Some(v) = t.viewer.as_mut() {
                if v.is_dirty() {
                    match v.save_to(None) {
                        Ok(_) => saved += 1,
                        Err(_) => failed += 1,
                    }
                }
            }
            if let Some(d) = t.doc.as_mut() {
                if d.is_dirty() {
                    match d.path.clone() {
                        Some(p) => match d.save_to(&p) {
                            Ok(_) => saved += 1,
                            Err(_) => failed += 1,
                        },
                        None => failed += 1,
                    }
                }
            }
        }
        self.status = if failed > 0 {
            format!("Save All: saved {saved} document(s), {failed} without a path")
        } else {
            format!("Save All: saved {saved} document(s)")
        };
    }

    /// File ▸ Close All Tabs。
    pub fn close_all_tabs(&mut self) {
        let n = self.tabs.len();
        let dirty = self.any_dirty();
        self.doc = None;
        self.viewer = None;
        self.tabs.clear();
        self.active = 0;
        DocViewState::default().apply(self);
        self.status = if dirty {
            format!("Closed {n} tab(s) (unsaved edits discarded)")
        } else {
            format!("Closed {n} tab(s)")
        };
    }

    fn edit_menu(&mut self, ui: &mut egui::Ui) {
        let (can_undo, can_redo) = self
            .doc
            .as_ref()
            .map(|d| (d.can_undo(), d.can_redo()))
            .unwrap_or((false, false));
        let has_doc = self.doc.is_some();
        let has_sel = self.effective_range().is_some();
        if item(ui, "Undo", Some("Ctrl+Z"), can_undo) {
            self.do_undo();
            ui.close();
        }
        if item(ui, "Redo", Some("Ctrl+Shift+Z"), can_redo) {
            self.do_redo();
            ui.close();
        }
        ui.separator();
        if item(ui, "Cut", Some("Ctrl+X"), has_sel) {
            self.do_cut();
            ui.close();
        }
        if item(ui, "Copy", Some("Ctrl+C"), self.selected_range().is_some()) {
            self.do_copy();
            ui.close();
        }
        let copy_has_sel = self.selected_range().is_some();
        if let Some(form) = copy_menu(ui, self, copy_has_sel) {
            self.do_copy_form(form);
            ui.close();
        }
        if item(ui, "Paste", Some("Ctrl+V"), has_doc) {
            self.request_paste(ui.ctx(),false);
            ui.close();
        }
        if item(ui,"Paste Reverse Complement",Some("Ctrl+Shift+V"),has_doc && !self.is_protein()) {self.request_paste(ui.ctx(),true);ui.close();}
        if item(ui, "Delete", Some("Del"), has_sel) {
            self.do_delete();
            ui.close();
        }
        if item(ui, "Enter Sequence Text…", None, has_doc) {
            self.show_paste_dialog = true;
            ui.close();
        }
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Find"), |ui| {
            item(ui, "Find DNA Sequence", Some("Ctrl+F"), false);
            item(ui, "Find Enzyme", None, false);
            item(ui, "Find Feature", None, false);
            item(ui, "Find Primer Data", None, false);
        });
        item(ui, "Go To…", Some("Ctrl+G"), false);
        if item(ui, "Select Range…", None, false) {}
        if item(ui, "Invert Selection", None, has_sel) {
            if let Some((s, e)) = self.effective_range() {
                if let Some(doc) = self.doc.as_ref() {
                    let len = doc.len();
                    if s > 0 && e < len {
                        // 反转：把选段之外的部分作为新选择
                        self.selection = Some((e, len));
                        self.selected_residue = None;
                        self.status = format!("Inverted selection → {} .. {}", e + 1, len);
                    } else if s == 0 && e < len {
                        self.selection = Some((e, len));
                        self.selected_residue = None;
                    } else if s > 0 && e == len {
                        self.selection = Some((0, s));
                        self.selected_residue = None;
                    }
                }
            }
            ui.close();
        }
        ui.separator();
        if item(
            ui,
            if self.aa_one_letter {
                "Use 3-Letter Amino Acid Codes"
            } else {
                "Use 1-Letter Amino Acid Codes"
            },
            None,
            true,
        ) {
            self.aa_one_letter = !self.aa_one_letter;
            ui.close();
        }
        if item(ui, "Make Uppercase", None, has_sel) {
            self.do_change_case(true);
            ui.close();
        }
        if item(ui, "Make Lowercase", None, has_sel) {
            self.do_change_case(false);
            ui.close();
        }
        ui.separator();
        item(ui, "Edit DNA Ends…", None, false);
        ui.menu_button(crate::i18n::text(ui.ctx(),"Topology"), |ui| {
            item(ui, "Topology (linear to circular)…", None, false);
            item(ui, "Topology (circular to linear)…", None, false);
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"Strandedness"), |ui| {
            item(ui, "Strandedness (ds to ss)…", None, false);
            item(ui, "Strandedness (ss to ds)…", None, false);
        });
        item(ui, "Methylation…", None, false);
        ui.separator();
        if item(ui, "Preferences…", Some("Ctrl+,"), true) {
            self.preferences.page = 0;
            self.enzyme_preferences_open = true;
            ui.close();
        }
    }

    fn view_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        if item(ui, "View Options…", None, has_doc) {
            self.show_view_options = true;
            ui.close();
        }
        let end_label = if self.show_end_labels {
            "Hide \"Start\" and \"End\""
        } else {
            "Show \"Start\" and \"End\""
        };
        if item(ui, end_label, None, true) {
            self.show_end_labels = !self.show_end_labels;
            ui.close();
        }
        item(ui, "Set Origin…", None, false);
        item(ui, "Set Translation Numbering…", None, false);
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Font Size"), |ui| {
            if item(ui, "Increase", Some("Ctrl++"), true) {
                self.zoom_in();
                ui.close();
            }
            if item(ui, "Decrease", Some("Ctrl+-"), true) {
                self.zoom_out();
                ui.close();
            }
        });
        let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Horizontal Scrolling")));
        ui.separator();
        if item_checked(ui, "Side Panel", self.sidebar_open) {
            self.sidebar_open = !self.sidebar_open;
        }
        if item_checked(ui, "Side Toolbar", self.show_side_tools) {
            self.show_side_tools = !self.show_side_tools;
        }
        if item_checked(ui, "Top Toolbar", self.show_toolbar) {
            self.show_toolbar = !self.show_toolbar;
        }
        ui.separator();
        if item_checked(ui, "Enzymes", self.show_sites) {
            self.show_sites = !self.show_sites;
        }
        if item_checked(ui, "Features", self.show_features) {
            self.show_features = !self.show_features;
        }
        if item_checked(ui, "Primers", self.show_primers) {
            self.show_primers = !self.show_primers;
        }
        if item_checked(ui, "Translations and ORFs", self.show_translation) {
            self.show_translation = !self.show_translation;
        }
        if item_checked(ui, "Colors", self.show_colors) {
            self.show_colors = !self.show_colors;
        }
        if item_checked(ui, "Ruler", self.show_ruler) {
            self.show_ruler = !self.show_ruler;
        }
        // 蛋白的三段特征开关（SnapGene 的 View 菜单同样按分子类型换内容）
        if self.is_protein() {
            if item_checked(ui, r#""Region" features"#, self.show_regions) {
                self.show_regions = !self.show_regions;
            }
            if item_checked(ui, r#""Site" and "Het" features"#, self.show_sites_het) {
                self.show_sites_het = !self.show_sites_het;
            }
            if item_checked(ui, r#""Bond" features"#, self.show_bonds) {
                self.show_bonds = !self.show_bonds;
            }
            if item_checked(ui, "Translation numbering", self.show_trans_numbering) {
                self.show_trans_numbering = !self.show_trans_numbering;
            }
        }
        ui.separator();
        if item(ui, "Fit to Window", Some("Ctrl+0"), true) {
            self.fit();
            ui.close();
        }
        item(ui, "Flip / Reverse Complement", None, false);
        item(ui, "Focus on Region", None, false);
    }

    fn enzymes_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        if item(ui, "Choose Enzymes…", Some("Ctrl+E"), has_doc) {
            self.tab = Tab::Enzymes;
            ui.close();
        }
        ui.menu_button(crate::i18n::text(ui.ctx(),"Show or Hide Enzymes"), |ui| {
            if item(ui, "Show All Enzymes", None, has_doc) {
                self.show_sites = true;
                self.status = "Showing all enzymes".into();
                ui.close();
            }
            if item(ui, "Hide All Enzymes", None, has_doc) {
                self.show_sites = false;
                self.status = "Hiding all enzymes".into();
                ui.close();
            }
        });
        item(ui, "Highlight Selected Enzyme Sites", None, false);
        item(ui, "Clear Highlighting for All Enzyme Sites", None, false);
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Use Enzyme Set"), |ui| {
            let mut chosen: Option<EnzymePreset> = None;
            for p in EnzymePreset::ALL {
                if item_checked(ui, p.label(), self.enzyme_preset == p) {
                    chosen = Some(p);
                }
            }
            if let Some(p) = chosen {
                self.apply_preset(p);
                ui.close();
            }
        });
        ui.separator();
        item(ui, "Remove Selected Enzyme Site…", None, false);
        item(ui, "Manage Enzyme Sets…", None, false);
        if item(ui, "Enzyme Database…", None, true){self.enzyme_manager.open=true;ui.close();}
        ui.separator();
        if item(ui, "Restriction Enzymes", None, has_doc) {
            self.tab = Tab::Enzymes;
            ui.close();
        }
        item(ui, "Noncutters", None, false);
        item(ui, "Enzymes Preferences", None, false);
    }

    /// Features 菜单（顺序与文案对齐真 SnapGene 8.2.2，见 8.2.2 截图「Features」展开）。
    fn features_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        let editable = self
            .doc
            .as_ref()
            .is_some_and(|d| d.native_class().is_none() && d.len() > 0);
        let has_feat = editable && self.selected_feature.is_some();

        if item(ui, "Add Feature…", Some("Ctrl+T"), editable) {
            self.open_add_feature_dialog(false);
            ui.close();
        }
        if item(ui, "Add Translated Feature…", None, editable) {
            self.open_add_feature_dialog(true);
            ui.close();
        }
        if item(ui, "Edit Feature…", None, has_feat) {
            self.open_edit_feature_dialog();
            ui.close();
        }
        if item(ui, "Duplicate Feature…", None, has_feat) {
            self.open_duplicate_feature_dialog();
            ui.close();
        }
        if item(ui, "Remove Feature", None, has_feat) {
            self.remove_selected_feature();
            ui.close();
        }
        ui.separator();
        // 片段级操作在 Edit Feature 对话框里做（SnapGene 也在那里）
        item(ui, "Create Feature Segment…", None, false);
        item(ui, "Delete Feature Segment", None, false);
        item(ui, "Merge Feature Segments…", None, false);
        item(ui, "Add Cleavage Site…", None, false);
        item(ui, "Remove Cleavage Site", None, false);
        ui.separator();
        if item(ui, "Feature Color…", None, has_feat) {
            self.feature_color_popup = true;
            ui.close();
        }
        if item(ui, "Show or Hide Feature", None, has_feat) {
            self.toggle_selected_feature_visible();
            ui.close();
        }
        let has_types = self
            .doc
            .as_ref()
            .map(|d| d.features().iter().any(|f| f.kind == "Site"))
            .unwrap_or(false);
        item(ui, "Show or Hide Feature Types…", None, has_types);
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Import Features"), |ui| {
            item(ui, "Import Features from a SnapGene File…", None, false);
            item(ui, "Import Features from a BED File…", None, false);
            item(ui, "Import Features from a GTF File…", None, false);
            item(ui, "Import Features from a GFF3 File…", None, false);
        });
        if item(ui, "Export Feature Data…", None, has_doc) {
            self.export_feature_data();
            ui.close();
        }
        ui.separator();
        item(ui, "Detect Common Features…", None, false);
        item(ui, "Add to Common Features", None, false);
        item(ui, "Browse Common Features", None, false);
        ui.separator();
        item(ui, "Add Feature Type…", None, false);
        item(ui, "Manage Feature Types…", None, false);
        ui.separator();
        ui.separator();
        if item(ui, "Select Feature Range", None, has_feat) {
            self.select_feature_range();
            ui.close();
        }
        if item(ui, "Center on Feature", None, has_feat) {
            if let Some(i) = self.selected_feature {
                self.center_on_feature(i);
            }
            ui.close();
        }
    }

    fn primers_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.as_ref().map(|d| d.is_dna()).unwrap_or(false);
        let has_sel = has_doc && self.active_primer_index().is_some();
        if item(ui, "Add Primer…", Some("Ctrl+R"), has_doc) {
            self.open_add_primer_dialog();
            ui.close();
        }
        if item(ui, "Edit Primer…", None, has_sel) {
            self.open_edit_primer_dialog();
            ui.close();
        }
        if item(ui, "Duplicate Primer…", Some("Ctrl+Shift+R"), has_sel) {
            self.duplicate_selected_primer();
            ui.close();
        }
        if item(ui, "Remove Primer", None, has_sel) {
            self.remove_selected_primer();
            ui.close();
        }
        // 菜单顺序按真 SnapGene 8.2.2 实测（`sg-primers-menu2-c.png`）：
        // Add / Edit / Duplicate / Remove ─ Primer Color / Hybridization Parameters
        // ─ Show or Hide ─ Import / Export ×2 / Analyze。
        ui.separator();
        item(ui, "Primer Color…", None, false);
        if item(ui, "Hybridization Parameters…", None, self.is_dna()) {
            self.hybridization_draft = self.doc.as_ref().map(|d| hx_formats::annot::parse_hybridization_params(d.block_payload(5).unwrap_or_default()));
            ui.close();
        }
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Show or Hide Primers"), |ui| {
            if item(ui, "Show All Primers", None, has_doc) {
                self.show_primers = true;
                self.status = "Showing all primers".into();
                ui.close();
            }
            if item(ui, "Hide All Primers", None, has_doc) {
                self.show_primers = false;
                self.status = "Hiding all primers".into();
                ui.close();
            }
        });
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Import Primers"), |ui| {
            item(ui, "Import Primers from a List…", None, false);
            item(ui, "Import Primers from a SnapGene File…", None, false);
        });
        item(ui, "Export Primer Data…", None, false);
        item(ui, "Export Primers to ssDNA…", None, false);
        // SnapGene 里这一项是付费锁；我们真算（Tm/GC/发夹/二聚体），所以可用。
        if item(ui, "Analyze Selected Primers…", None, has_sel) {
            self.analyze_selected_primer();
            ui.close();
        }
    }

    fn actions_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Restriction and Insertion Cloning"), |ui| {
            if item(ui, "Insert Fragment…", Some("Ctrl+Alt+I"), has_doc) {
                self.open_cloning(hx_algo::cloning::Method::Restriction);
                ui.close();
            }
            if item(ui, "Insert Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::RestrictionMulti);
                ui.close();
            }
            item(ui, "Delete Restriction Fragment…", None, false);
            item(ui, "Destroy Restriction Site…", None, false);
            ui.separator();
        });
        if item(ui, "Linear Ligation…", None, has_doc) {
            self.open_cloning(hx_algo::cloning::Method::Ligation);
            ui.close();
        }
        ui.separator();
        if item(ui, "PCR…", Some("Ctrl+D"), has_doc) {
            self.open_design_dialog();
            ui.close();
        }
        item(ui, "Overlap Extension PCR…", None, false);
        item(ui, "Mutagenesis…", Some("Ctrl+J"), false);
        ui.separator();
        ui.menu_button(crate::i18n::text(ui.ctx(),"Gateway® Cloning"), |ui| {
            if item(ui, "BP Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GatewayBp);
                ui.close();
            }
            for n in 1..=4 {
                let label = match n {
                    1 => "LR Cloning: One Fragment…",
                    2 => "LR Cloning: Two Fragments…",
                    3 => "LR Cloning: Three Fragments…",
                    _ => "LR Cloning: Four Fragments…",
                };
                if item(ui, label, None, has_doc) {
                    self.open_cloning(hx_algo::cloning::Method::GatewayLr(n));
                    ui.close();
                }
            }
            for n in 1..=4 {
                let label = match n {
                    1 => "BP + LR Cloning: One Fragment…",
                    2 => "BP + LR Cloning: Two Fragments…",
                    3 => "BP + LR Cloning: Three Fragments…",
                    _ => "BP + LR Cloning: Four Fragments…",
                };
                if item(ui, label, None, has_doc) {
                    self.open_cloning(hx_algo::cloning::Method::GatewayBpLr(n));
                    ui.close();
                }
            }
            ui.separator();
            ui.menu_button(crate::i18n::text(ui.ctx(),"Donor Vectors"), |ui| {
                self.vector_library_menu(ui, "Gateway Donor Vectors");
            });
            ui.menu_button(crate::i18n::text(ui.ctx(),"Destination Vectors"), |ui| {
                self.vector_library_menu(ui, "Gateway Destination Vectors");
            });
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"Gibson Assembly®"), |ui| {
            if item(ui, "Insert Fragment…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::Gibson);
                ui.close();
            }
            if item(ui, "Insert Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GibsonMulti);
                ui.close();
            }
            if item(ui, "Assemble Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GibsonAssemble);
                ui.close();
            }
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"Golden Gate Assembly"), |ui| {
            if item(ui, "Insert Fragment…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GoldenGate);
                ui.close();
            }
            if item(ui, "Insert Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GoldenGateMulti);
                ui.close();
            }
            if item(ui, "Golden Gate Settings…", None, true) {
                self.gg_settings_backup = Some(self.gg_settings.clone());
                self.gg_settings_open = true;
                ui.close();
            }
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"In-Fusion® Cloning"), |ui| {
            if item(ui, "Insert Fragment…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::InFusion);
                ui.close();
            }
            if item(ui, "Insert Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::InFusionMulti);
                ui.close();
            }
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"NEBuilder® HiFi DNA Assembly"), |ui| {
            if item(ui, "Insert Fragment…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::NeBuilder);
                ui.close();
            }
            if item(ui, "Insert Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::NeBuilderMulti);
                ui.close();
            }
            if item(ui, "Assemble Multiple Fragments…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::NeBuilderAssemble);
                ui.close();
            }
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"TA or GC Cloning"), |ui| {
            if item(ui, "TA Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::TaCloning);
                ui.close();
            }
            if item(ui, "GC Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::GcCloning);
                ui.close();
            }
            ui.separator();
            ui.menu_button(crate::i18n::text(ui.ctx(),"TA Vectors"), |ui| {
                self.vector_library_menu(ui, "TA Cloning Vectors");
            });
            ui.menu_button(crate::i18n::text(ui.ctx(),"GC Vectors"), |ui| {
                self.vector_library_menu(ui, "GC Cloning Vectors");
            });
            ui.separator();
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"TOPO® Cloning"), |ui| {
            if item(ui, "TOPO® TA Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::TopoTa);
                ui.close();
            }
            if item(ui, "Blunt TOPO® Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::TopoBlunt);
                ui.close();
            }
            if item(ui, "Directional TOPO® Cloning…", None, has_doc) {
                self.open_cloning(hx_algo::cloning::Method::TopoDirectional);
                ui.close();
            }
            ui.separator();
            ui.menu_button(crate::i18n::text(ui.ctx(),"TOPO® TA Vectors"), |ui| {
                self.vector_library_menu(ui, "TOPO TA Vectors");
            });
            ui.menu_button(crate::i18n::text(ui.ctx(),"Blunt TOPO® Vectors"), |ui| {
                self.vector_library_menu(ui, "TOPO Blunt Vectors");
            });
            ui.menu_button(crate::i18n::text(ui.ctx(),"Directional TOPO® Vectors"), |ui| {
                self.vector_library_menu(ui, "TOPO Directional Vectors");
            });
            ui.separator();
        });
        ui.separator();
        for (label, action) in [("Linearize…", crate::molecular_actions::Action::Linearize), ("Convert to Single-Stranded…", crate::molecular_actions::Action::SingleStrand)] {
            if item(ui, label, None, action.allowed(self.doc.as_ref())) { self.open_molecular_action(action); ui.close(); }
        }
        ui.separator();
        for (label, shortcut, action) in [("Anneal Oligos…", None, crate::molecular_actions::Action::Anneal), ("Make RNA…", None, crate::molecular_actions::Action::Rna), ("Make Protein…", Some("Ctrl+K"), crate::molecular_actions::Action::Protein)] {
            if item(ui, label, shortcut, action.allowed(self.doc.as_ref())) { self.open_molecular_action(action); ui.close(); }
        }
    }

    /// 打开克隆对话框：快照当前所有打开的文档。
    fn open_cloning(&mut self, method: hx_algo::cloning::Method) {
        let mut docs: Vec<crate::clonedlg::DocSnap> = Vec::new();
        let mut active_pos = None;
        for (i, tab) in self.tabs.iter().enumerate() {
            let d = if i == self.active {
                self.doc.as_ref()
            } else {
                tab.doc.as_ref()
            };
            let Some(d) = d else { continue };
            let name = tab.title(Some(d), None);
            let name = name.trim_end_matches(" *").to_string();
            if i == self.active {
                active_pos = Some(docs.len());
            }
            let src = d.cloning_source(name.clone());
            let mut snap = crate::clonedlg::DocSnap::new(name, src);
            if let Ok(bytes) = d.session.save_bytes() {
                snap = snap.with_bytes(bytes);
            }
            docs.push(snap);
        }
        self.cloning = Some(crate::clonedlg::CloningDlg::new_with_settings_and_gg(
            method,
            docs,
            active_pos,
            self.cloning_settings.clone(),
            &self.gg_settings,
        ));
    }

    /// 克隆对话框窗口（每帧）。
    fn cloning_dlg_window(&mut self, ctx: &egui::Context) {
        if self.cloning.is_none() {
            return;
        }
        let Some(dlg) = self.cloning.as_mut() else {
            return;
        };
        match crate::clonedlg::draw(dlg, ctx) {
            Some(crate::clonedlg::Outcome::Clone(bytes, name, close)) => {
                if close {
                    self.cloning = None;
                }
                self.adopt_cloned_product(bytes, name);
                // 新产物带来新的特征标签/字形，字体图集可能在同一帧扩容；
                // 再排一帧，避免侧栏和图标停留在旧 UV（显示成同一个字形）。
                ctx.request_repaint();
            }
            Some(crate::clonedlg::Outcome::Cancel) => {
                self.cloning = None;
                self.status = "Cloning cancelled".into();
                ctx.request_repaint();
            }
            None => {}
        }
    }

    /// 克隆成功后把产物接进主窗口。
    ///
    /// 产物是**未保存文档**：必须同时把 `save_path` 指向「根目录 + 产物名」。
    /// 否则随后按 Ctrl+S / Save As 会沿用上一个文档的路径（常见是载体），
    /// 在没有原生对话框的环境里会直接覆盖源文件。
    fn adopt_cloned_product(&mut self, bytes: Vec<u8>, name: String) {
        match crate::doc::Document::from_bytes(bytes) {
            Ok(mut doc) => {
                doc.name_hint = Some(name.clone());
                self.push_doc(doc, None);
                let i = self.active;
                self.tabs[i].untitled_name = Some(name.clone());
                let dir = self
                    .folder
                    .clone()
                    .unwrap_or_else(|| std::env::current_dir().unwrap_or_default());
                self.save_path = dir.join(&name).display().to_string();
                self.status = format!("Created {name} — not saved yet");
            }
            Err(e) => self.status = format!("Clone failed: {e}"),
        }
    }

    fn open_molecular_action(&mut self, action: crate::molecular_actions::Action) {
        match crate::molecular_actions::Dialog::new(action, self.doc.as_ref(), self.selection) {
            Ok(dlg) => self.molecular_action = Some(dlg),
            Err(e) => self.status = e,
        }
    }

    fn tools_menu(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        item(ui, "Assemble Contigs…", None, false);
        item(ui, "Align Open Sequences…", None, false);
        ui.separator();
        item(ui, "Show Codon Usage Tables", None, false);
        item(ui, "Genetic Codes", None, false);
        item(ui, "Letter Codes", None, false);
        item(ui, "Show DNA Calculations", None, false);
        ui.separator();
        if item(ui, "Recompute Enzyme Sites", None, has_doc) {
            self.sites_dirty = true;
            ui.close();
        }
        if item(ui, "Block Summary (status bar)", None, has_doc) {
            self.dump_blocks();
            ui.close();
        }
    }

    fn window_menu(&mut self, ui: &mut egui::Ui) {
        let mut want_activate: Option<usize> = None;
        if item(ui, "Minimize", None, true) {
            ui.ctx()
                .send_viewport_cmd(egui::ViewportCommand::Minimized(true));
            ui.close();
        }
        if item(ui, "Zoom", None, true) {
            let max = ui.ctx().input(|i| i.viewport().maximized.unwrap_or(false));
            ui.ctx()
                .send_viewport_cmd(egui::ViewportCommand::Maximized(!max));
            ui.close();
        }
        if item(ui, "Enter Full Screen", Some("F11"), true) {
            let fs = ui.ctx().input(|i| i.viewport().fullscreen.unwrap_or(false));
            ui.ctx()
                .send_viewport_cmd(egui::ViewportCommand::Fullscreen(!fs));
            ui.close();
        }
        ui.separator();
        if item(ui, "Show Next Tab", Some("Ctrl+Tab"), self.tabs.len() > 1) {
            self.next_tab(1);
            ui.close();
        }
        if item(ui, "Show Previous Tab", Some("Ctrl+Shift+Tab"), self.tabs.len() > 1) {
            self.next_tab(-1);
            ui.close();
        }
        if item(ui, "Move Tab to New Window", None, false) {}
        if item(ui, "Move All Tabs to New Windows", None, false) {}
        if item(ui, "Merge All Windows", None, false) {}
        // ---- 底部：已打开文档列表（oracle：Window 菜单底部列出全部标签页，带类型色圆点）----
        if !self.tabs.is_empty() {
            ui.separator();
            let rows: Vec<(String, Option<hx_core::MoleculeKind>)> = (0..self.tabs.len())
                .map(|i| {
                    let d = if i == self.active {
                        self.doc.as_ref()
                    } else {
                        self.tabs[i].doc.as_ref()
                    };
                    let vw = if i == self.active {
                        self.viewer.as_ref()
                    } else {
                        self.tabs[i].viewer.as_ref()
                    };
                    (
                        self.tabs[i].title(d, vw),
                        d.map(|d| d.kind()),
                    )
                })
                .collect();
            let active = self.active;
            for (i, (title, kind)) in rows.iter().enumerate() {
                let (r, g, b) = kind
                    .map(|k| k.icon_rgb())
                    .unwrap_or((0x9a, 0x9a, 0x9a));
                let text = egui::RichText::new(format!("● {title}"))
                    .color(egui::Color32::from_rgb(r, g, b));
                if ui
                    .selectable_label(i == active, text)
                    .on_hover_text(format!(
                        "{} document",
                        kind.map(|k| k.display_name()).unwrap_or("Read-only")
                    ))
                    .clicked()
                {
                    want_activate = Some(i);
                    ui.close();
                }
            }
        }
        ui.separator();
        if item(ui, "Fit to Window", Some("Ctrl+0"), true) {
            self.fit();
            ui.close();
        }
        if item(ui, "Reset View", None, true) {
            self.reset_view();
            ui.close();
        }
        if let Some(i) = want_activate {
            self.activate_tab(i);
        }
    }

    fn help_menu(&mut self, ui: &mut egui::Ui) {
        if item(ui, "User Guide", None, false) {}
        if item(ui, "Keyboard Shortcuts", None, true) {
            self.show_shortcuts = true;
            ui.close();
        }
        if item(ui, "Symbols and Abbreviations", None, false) {}
        if item(ui, "Search Tips", None, false) {}
        ui.separator();
        item(ui, "Check for Updates…", None, false);
        item(ui, "Report a Bug…", None, false);
        item(ui, "Online Resources", None, false);
        ui.separator();
        if item(ui, "About AppimageGene", None, true) {
            self.show_about = true;
            ui.close();
        }
    }

    fn save_button_menu(&mut self,ui:&mut egui::Ui) {
        ui.set_min_width(280.0);
        let has_doc=self.doc.is_some();
        if item(ui,"Save",Some("Ctrl+S"),has_doc){self.do_save(false);ui.close();}
        if item(ui,"Save As…",Some("Ctrl+Shift+S"),has_doc){self.do_save(true);ui.close();}
        ui.separator();
        let molecule=self.doc.as_ref().map(|d|d.kind().display_name()).unwrap_or("DNA");
        if item(ui,&format!("Export {molecule}…"),None,has_doc){self.export_sequence_data(false);ui.close();}
        if item(ui,&format!("Export Selected {molecule}…"),None,has_doc && self.effective_range().is_some()){self.export_sequence_data(true);ui.close();}
        item(ui,"Export Selected Protein…",None,false);
        if item(ui,"Export Map…",None,has_doc){self.export_map_from_menu();ui.close();}
        item(ui,"Export History…",None,false);ui.separator();
        if item(ui,"Export Feature Data…",None,has_doc){self.export_annotation_data(false);ui.close();}
        if item(ui,"Export Primer Data…",None,has_doc){self.export_annotation_data(true);ui.close();}
    }

    fn paste_button_menu(&mut self,ui:&mut egui::Ui) {
        ui.set_min_width(290.0);
        if item(ui,"Paste",Some("Ctrl+V"),self.doc.is_some()){self.request_paste(ui.ctx(),false);ui.close();}
        if item(ui,"Paste Reverse Complement",Some("Ctrl+Shift+V"),self.doc.is_some() && !self.is_protein()){
            self.request_paste(ui.ctx(),true);ui.close();
        }
    }

    fn exported_sequence_text(&self,selected:bool)->Option<String> {
        let doc=self.doc.as_ref()?;
        let bases=if selected {let (a,b)=self.selected_range()?;
            if a<=b {doc.sequence()[a..b].to_vec()}else{[&doc.sequence()[a..],&doc.sequence()[..b]].concat()}}
            else{doc.sequence().to_vec()};
        Some(format!(">{}{}\n{}\n",doc.display_name().replace(['\n','\r']," "),if selected {" selected"}else{""},
            bases.chunks(70).map(|chunk|chunk.iter().map(|b|doc.kind().display_base(*b) as char).collect::<String>()).collect::<Vec<_>>().join("\n")))
    }

    fn export_sequence_data(&mut self,selected:bool) {
        let Some(doc)=self.doc.as_ref() else{return};
        let index=match doc.kind(){hx_core::MoleculeKind::Dna=>0,hx_core::MoleculeKind::Rna=>1,hx_core::MoleculeKind::Protein=>2};
        let Some(fasta)=self.exported_sequence_text(selected) else{return};
        let (text,extension)=match self.preferences.export_formats[index] {
            1=>{
                let text=if selected {
                    let Some((a,b))=self.selected_range() else{return};
                    let part=doc.cloning_source("selection").slice(a,b,false);
                    let mut file=hx_formats::DnaFile::parse(&hx_formats::write::NewFile::of(doc.kind(),part.seq).with_flags(doc.session.flags().with(hx_core::MoleculeFlags::CIRCULAR,false)).to_bytes()).unwrap();
                    let payload=hx_formats::write::feature_annotations_xml(&part.features,file.sequence().unwrap().len());
                    if let Some(block)=file.blocks.iter_mut().find(|b|b.block_type==10){block.payload=payload;}
                    let selected=crate::doc::Document::from_bytes(file.to_bytes()).unwrap();
                    crate::preferences::genbank_export(&selected,&self.preferences)
                }else{crate::preferences::genbank_export(doc,&self.preferences)};
                (text,if index==2{"gp"}else{"gb"})
            }
            2=>(fasta.lines().skip(1).collect::<Vec<_>>().join("\n")+"\n","txt"),
            _=>(fasta,"fasta"),
        };
        let text=if self.preferences.export_crlf{text.replace('\n',"\r\n")}else{text};
        self.save_export_text("Export Sequence",extension,text);
    }

    fn export_annotation_data(&mut self,primers:bool) {
        let Some(doc)=self.doc.as_ref() else{return};let mut text="Name\tStart\tEnd\tType\n".to_string();
        if primers {for p in doc.primers(){for s in &p.sites {text.push_str(&format!("{}\t{}\t{}\tPrimer\n",p.name.replace(['\t','\n','\r']," "),s.range.start+1,s.range.end));}}}
        else{for f in doc.features(){for s in &f.segments{text.push_str(&format!("{}\t{}\t{}\t{}\n",f.name.replace(['\t','\n','\r']," "),s.range.start+1,s.range.end,f.kind.replace(['\t','\n','\r']," ")));}}}
        self.save_export_text(if primers {"Export Primer Data (TSV)"}else{"Export Feature Data (TSV)"},if primers {"primers.tsv"}else{"features.tsv"},text);
    }

    fn save_export_text(&mut self,title:&str,extension:&str,text:String) {
        let Some(doc)=self.doc.as_ref() else{return};
        let default=doc.path.clone().unwrap_or_else(||PathBuf::from(&self.save_path)).with_extension(extension);
        let Some(path)=crate::sysopen::pick_save_path(title,&default) else{return};
        match std::fs::write(&path,text.as_bytes()){Ok(())=>self.status=format!("Exported {}",path.display()),Err(e)=>self.status=format!("Export failed: {e}")}
    }

    fn export_map_from_menu(&mut self) {
        let Some(doc)=self.doc.as_ref() else{return};
        let ring=self.map_layout.is_ring(doc.circular(),doc.len());
        let svg=doc.export_map_svg(&self.sites.selected,&hx_algo::export::default_map_style(ring)).svg;
        self.save_export_text("Export Map (SVG)","map.svg",svg);
    }

    /// SnapGene 式工具栏：图标 + 文字（Save / Undo / Redo / Cut / Copy / Paste / Print）。
    fn toolbar(&mut self, ui: &mut egui::Ui) {
        let has_doc = self.doc.is_some();
        let has_sel = self.effective_range().is_some();
        let (can_undo, can_redo) = self
            .doc
            .as_ref()
            .map(|d| (d.can_undo(), d.can_redo()))
            .unwrap_or((false, false));
        ui.spacing_mut().item_spacing.x = 2.0;
        ui.horizontal(|ui| {
            let (save_main,_)=icons::tool_split_menu(ui,icons::Icon::Save,"Save",has_doc,|ui|self.save_button_menu(ui));
            if save_main
                .on_hover_text(crate::i18n::text(ui.ctx(),"Save (Ctrl+S)"))
                .clicked()
            {
                self.do_save(false);
            }
            if icons::tool(ui, icons::Icon::Undo, "Undo", can_undo)
                .on_hover_text(crate::i18n::text(ui.ctx(),"Undo (Ctrl+Z)"))
                .clicked()
            {
                self.do_undo();
            }
            if icons::tool(ui, icons::Icon::Redo, "Redo", can_redo)
                .on_hover_text(crate::i18n::text(ui.ctx(),"Redo (Ctrl+Shift+Z)"))
                .clicked()
            {
                self.do_redo();
            }
            ui.add_space(4.0);
            ui.separator();
            ui.add_space(4.0);
            if icons::tool(ui, icons::Icon::Cut, "Cut", has_sel)
                .on_hover_text(crate::i18n::text(ui.ctx(),"Cut selection (Ctrl+X)"))
                .clicked()
            {
                self.do_cut();
            }
            let copy_has_sel = self.selected_range().is_some();
            let (copy_main, _) = icons::tool_split_menu(
                ui,
                icons::Icon::Copy,
                "Copy",
                copy_has_sel,
                |ui| {
                    if let Some(form) = copy_menu(ui, self, copy_has_sel) {
                        self.do_copy_form(form);
                        ui.close();
                    }
                },
            );
            if copy_main.on_hover_text(crate::i18n::text(ui.ctx(),"Copy top strand bases (Ctrl+C)")).clicked() {
                self.do_copy();
            }
            let (paste_main,_)=icons::tool_split_menu(ui,icons::Icon::Paste,"Paste",has_doc,|ui|self.paste_button_menu(ui));
            if paste_main
                .on_hover_text(crate::i18n::text(ui.ctx(),"Paste clipboard at cursor (Ctrl+V)"))
                .clicked()
            {
                self.request_paste(ui.ctx(),false);
            }
            ui.add_space(4.0);
            ui.separator();
            ui.add_space(4.0);
            let _ = icons::tool(ui, icons::Icon::Print, "Print", false)
                .on_hover_text(crate::i18n::text(ui.ctx(),"Not implemented yet"));
            // 右端：`名字 (N bp)` + `Info`（SnapGene 实测：这两样在工具栏右端，
            // 不在文档页签行里——页签行要留给页签本身）
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add_space(8.0);
                if icons::tool(ui,icons::Icon::Info,"Info",has_doc)
                    .on_hover_text(crate::i18n::text(ui.ctx(),"File information (metadata, blocks, features)"))
                    .clicked()
                {
                    self.show_info = !self.show_info;
                }
                if let Some(d) = self.doc.as_ref() {
                    ui.label(
                        egui::RichText::new(d.summary())
                            .size(11.5)
                            .color(theme::TEXT_DIM),
                    );
                }
            });
        });
    }

    /// 文档页签行（SnapGene 顶部：`[tab 名字 ×] [tab2 ×] … [+] ... 名字 (N bp) [Info]`）。
    fn doc_tabs(&mut self, ui: &mut egui::Ui) -> Option<egui::Rect> {
        let mut active_rect=None;
        // 活动标签页的文档在 `self.doc`，其余在 `self.tabs[i].doc`。
        let titles: Vec<String> = (0..self.tabs.len())
            .map(|i| {
                let d = if i == self.active {
                    self.doc.as_ref()
                } else {
                    self.tabs[i].doc.as_ref()
                };
                let vw = if i == self.active {
                    self.viewer.as_ref()
                } else {
                    self.tabs[i].viewer.as_ref()
                };
                self.tabs[i].title(d, vw)
            })
            .collect();
        // 文档摘要 `名字 (N bp)` 现在画在工具栏右端（见 [`HelixApp::toolbar`]），页签行留给页签。
        // 每个标签页的分子类型（用于图标着色）
        // `None` = 原生但非序列的文档（比对/凝胶），页签画灰页而不是 DNA 双螺旋
        let kinds: Vec<Option<hx_core::MoleculeKind>> = (0..self.tabs.len())
            .map(|i| {
                let d = if i == self.active {
                    self.doc.as_ref()
                } else {
                    self.tabs[i].doc.as_ref()
                };
                d.and_then(|d| d.native_class().is_none().then(|| d.kind()))
            })
            .collect();
        let active = self.active;
        let mut want_activate: Option<usize> = None;
        let mut want_close: Option<usize> = None;
        let mut want_open = false;
        let mut want_folder = false;
        let mut want_new: Option<hx_core::MoleculeKind> = None;
        let mut want_drop: Option<(usize, usize)> = None;

        // ---- 页签宽度预算 ----
        // SnapGene 实测（14 个文档）：页签**等比压窄**、名字用 `…` 截断，而不是溢出到
        // 窗口外（我们原来会把后面的页签挤出去，用户点不到）。页签行的右端只有 `+` 与 `⋮`
        // （`名字 (N bp)` 与 `Info` 在工具栏右端，见 [`HelixApp::toolbar`]）。
        const MIN_TAB_W: f32 = 62.0;
        // 左边距 8 + `+` 28 + 间距 6 + `⋮` 28 + 右边距 10
        let reserve = 8.0 + 28.0 + 6.0 + 28.0 + 10.0;
        let avail = (ui.available_width() - reserve).max(MIN_TAB_W);
        // Standard-width tabs, compressed evenly only when the row is crowded.
        let widths=vec![crate::chrome::document_tab_width(avail+reserve,titles.len());titles.len()];

        ui.horizontal_top(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            ui.spacing_mut().item_spacing.y=0.0;
            if !titles.is_empty() {
            egui::ScrollArea::horizontal().id_salt("document_tab_strip").scroll_bar_visibility(egui::scroll_area::ScrollBarVisibility::AlwaysHidden).max_width(avail).auto_shrink([true,true]).show(ui,|ui|{
            ui.horizontal_top(|ui|{
            ui.spacing_mut().item_spacing.x=0.0;
            let mut tab_rects: Vec<egui::Rect> = Vec::with_capacity(titles.len());
            for (i, title) in titles.iter().enumerate() {
                let is_active = i == active;
                let (resp, close) =
                    tab_button_with_close(ui, is_active, title, i, kinds[i], widths[i]);
                tab_rects.push(resp.rect);
                if is_active {active_rect=Some(resp.rect.intersect(ui.clip_rect()));}
                // 拖动页签：记录起点，并在每一帧按指针位置计算插入点。
                if resp.drag_started() {
                    self.tab_drag = Some(i);
                    self.tab_drag_at = Some(i);
                }
                if resp.dragged() {
                    ui.ctx().set_cursor_icon(egui::CursorIcon::Grabbing);
                    if let Some(x) = ui.input(|inp| inp.pointer.interact_pos()).map(|p| p.x) {
                        self.tab_drag_at = Some(
                            tab_rects.iter().filter(|r| x > r.center().x).count(),
                        );
                    }
                }
                if resp.drag_stopped() {
                    if let Some(from) = self.tab_drag.take() {
                        want_drop = Some((from, self.tab_drag_at.take().unwrap_or(i)));
                    }
                }
                if close.clicked() {
                    want_close = Some(i);
                } else if resp.clicked() {
                    want_activate = Some(i);
                }
                resp.on_hover_text(if is_active {
                    format!("{title} — click × to close (Ctrl+W), drag to reorder")
                } else {
                    format!("{title} — click to switch (Ctrl+Tab), drag to reorder")
                });
                close.on_hover_text(crate::i18n::text(ui.ctx(),"Close tab (Ctrl+W)"));
            }
            // 拖动时在插入位置画一条蓝色竖线。
            if let (Some(_), Some(at)) = (self.tab_drag, self.tab_drag_at) {
                if let (Some(first), Some(last)) = (tab_rects.first(), tab_rects.last()) {
                    let x = if at == 0 {
                        first.left()
                    } else if at >= tab_rects.len() {
                        last.right()
                    } else {
                        tab_rects[at].left()
                    };
                    ui.painter().line_segment(
                        [
                            egui::pos2(x, first.top() + 3.0),
                            egui::pos2(x, first.bottom() - 3.0),
                        ],
                        egui::Stroke::new(2.0, theme::ACCENT),
                    );
                }
            }
            // 指针已经松开但没收到 drag_stopped（例如在窗口外释放）：清掉拖动状态。
            if self.tab_drag.is_some() && !ui.input(|inp| inp.pointer.any_down()) {
                self.tab_drag = None;
                self.tab_drag_at = None;
            }
            });
            });
            }
            // SnapGene 8 的 `+`：既能新建序列，也能打开文件/文件夹
            icons::plus_menu(ui, |ui| {
                if ui.button(crate::i18n::text(ui.ctx(),"New DNA or RNA File…")).clicked() {
                    want_new = Some(hx_core::MoleculeKind::Dna);
                    ui.close();
                }
                if ui.button(crate::i18n::text(ui.ctx(),"New Protein File…")).clicked() {
                    want_new = Some(hx_core::MoleculeKind::Protein);
                    ui.close();
                }
                ui.separator();
                if ui.button(crate::i18n::text(ui.ctx(),"Open Files…")).clicked() {
                    want_open = true;
                    ui.close();
                }
                if ui.button(crate::i18n::text(ui.ctx(),"Open Folder…")).clicked() {
                    want_folder = true;
                    ui.close();
                }
            })
            .0
            .on_hover_text(crate::i18n::text(ui.ctx(),"New sequence or open files/folder"));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add_space(6.0);
                // SnapGene 的 `⋮`：列出所有已打开文档（点击切换）
                icons::more_menu(ui, |ui| {
                    for (i, title) in titles.iter().enumerate() {
                        let (r, g, b) = kinds[i]
                            .map(|k| k.icon_rgb())
                            .unwrap_or((0x9a, 0x9a, 0x9a));
                        let text = egui::RichText::new(format!("● {title}"))
                            .color(egui::Color32::from_rgb(r, g, b));
                        if ui
                            .selectable_label(i == active, text)
                            .on_hover_text(format!(
                                "{} document",
                                kinds[i].map(|k| k.display_name()).unwrap_or("SnapGene")
                            ))
                            .clicked()
                        {
                            want_activate = Some(i);
                            ui.close();
                        }
                    }
                })
                .0
                .on_hover_text(crate::i18n::text(ui.ctx(),"All open documents"));
            });
        });
        if let Some((from, at)) = want_drop {
            self.move_tab(from, at);
        }
        if let Some(i) = want_activate {
            self.activate_tab(i);
        }
        if let Some(i) = want_close {
            self.close_tab(i);
        }
        if let Some(kind) = want_new {
            self.open_new_file_dialog(kind);
        }
        if want_folder {
            self.open_folder_with_dialog();
        }
        if want_open {
            // SnapGene 的 `+`：走系统原生文件对话框（KDE 上是 Dolphin 那套）
            self.open_files_with_dialog();
        }
        active_rect
    }

    /// Navigation and document tabs occupy one full-width row. Collapsing the
    /// project panel removes it entirely; it must not leave an 18px vertical rail.
    fn navigation_row(&mut self,ui:&mut egui::Ui) {
        let row=ui.max_rect();
        let mut active=None;
        ui.horizontal_top(|ui|{
            ui.spacing_mut().item_spacing.x=0.0;
            ui.spacing_mut().item_spacing.y=0.0;
            let width=crate::chrome::navigation_width(self.sidebar_open,self.sidebar_width);
            ui.allocate_ui_with_layout(egui::vec2(width,crate::chrome::DOC_TABS_H),egui::Layout::left_to_right(egui::Align::Center),|ui|{
                ui.add_space(2.0);
                if icons::navigation(ui,if self.sidebar_open {icons::Icon::Sidebar}else{icons::Icon::SidebarClosed}).on_hover_text(crate::i18n::text(ui.ctx(),if self.sidebar_open {"Hide Files panel"}else{"Show Files panel"})).clicked(){self.sidebar_open=!self.sidebar_open;}
                if self.sidebar_open {ui.with_layout(egui::Layout::right_to_left(egui::Align::Center),|ui|{
                    // 与下面侧栏 `⋮` 共用同一右内缩，两个图标的中心列对齐。
                    ui.add_space(crate::chrome::NAV_RIGHT_INSET);
                    if icons::navigation(ui,icons::Icon::Home).on_hover_text(crate::i18n::text(ui.ctx(),"Open Folder…")).clicked(){self.open_folder_with_dialog();}
                    if icons::navigation(ui,icons::Icon::Settings).on_hover_text(crate::i18n::text(ui.ctx(),"View Options…")).clicked(){self.show_view_options=true;}
                });}else{
                    if icons::navigation(ui,icons::Icon::Settings).on_hover_text(crate::i18n::text(ui.ctx(),"View Options…")).clicked(){self.show_view_options=true;}
                    if icons::navigation(ui,icons::Icon::Home).on_hover_text(crate::i18n::text(ui.ctx(),"Open Folder…")).clicked(){self.open_folder_with_dialog();}
                }
            });
            active=self.doc_tabs(ui);
        });
        // Paint after the tab scroll area, inside the row clip. A baseline at
        // rect.bottom() was clipped/overpainted by the toolbar's white panel.
        let y=row.top()+crate::chrome::DOC_TABS_H-0.5;
        let mut outline=Vec::new();
        // SG 实测（`oracle-2026-10-18-many-tabs-snapgene.png`，y=102）：这条
        // 蓝线横跨整行，包括左侧导航图标区；只在活动文档下方断开。
        let start=row.left();
        for (a,b) in crate::chrome::tab_separator_segments(start,row.right(),active.map(|r|(r.left(),r.right()))) {
            outline.push([egui::pos2(a,y),egui::pos2(b,y)]);
        }
        if let Some(rect)=active {
            // Keep both edges a half pixel INSIDE the scroll clip, including
            // the first tab; drawing on its left boundary lost half the stroke.
            for x in [rect.left()+0.5,rect.right()-0.5] {
                outline.push([egui::pos2(x,row.top()+0.5),egui::pos2(x,y)]);
            }
        }
        // SG keeps a vertical rule on the navigation/tab boundary for the whole
        // row, joined to the sidebar edge below — not just a tab border that
        // disappears when the first tab is active. Blue while the Files panel is
        // open, neutral grey when it is collapsed. Painted before the blue
        // outline so the horizontal separator still wins at the junction.
        let nav_right=row.left()+crate::chrome::navigation_width(self.sidebar_open,self.sidebar_width);
        // The left SidePanel's one-pixel stroke is centred half a pixel inside
        // its response edge (the response is inset from the nominal 230 px
        // navigation width). Match that stroke centre rather than the nominal
        // width, otherwise the boundary makes a visible one-pixel jog at y=56.
        let boundary_x=nav_right-1.5;
        let boundary=egui::Rect::from_min_max(
            egui::pos2(boundary_x-0.5,row.top()+2.0),
            egui::pos2(boundary_x+0.5,row.top()+crate::chrome::DOC_TABS_H-0.5));
        let boundary_color=if self.sidebar_open {crate::sidebar::EDGE}
            else{egui::Color32::from_rgb(0xd9,0xd9,0xd9)};
        ui.ctx().data_mut(|d|d.insert_temp(egui::Id::unique("document_row_nav_boundary"),(boundary,boundary_color)));
        ui.ctx().data_mut(|d|d.insert_temp(egui::Id::unique("document_row_outline"),outline));
        ui.ctx().data_mut(|d|d.insert_temp(egui::Id::unique("document_row_active_gap"),
            active.map(|r|egui::Rect::from_min_max(egui::pos2(r.left()+1.0,y-0.5),egui::pos2(r.right()-1.0,y+0.5)))));
    }

    /// 视图页签行：按分子类型决定（DNA：Map/Sequence/Enzymes/Features/Primers/History；
    /// RNA：Map/Sequence/Features/Secondary Structure/History；
    /// 蛋白：Map/Sequence/Properties/Features/History）。
    fn view_tabs(&mut self, ui: &mut egui::Ui, tabs: &[Tab]) {
        ui.spacing_mut().item_spacing.x = 0.0;
        if !tabs.contains(&self.tab) {
            // 换文档/换分子类型后停在当前文档没有的视图上 → 回到该文档的第一个页签
            self.tab = tabs.first().copied().unwrap_or(Tab::Map);
        }
        ui.horizontal(|ui| {
            for t in tabs {
                let t = *t;
                if tab_button(ui, self.tab == t, crate::i18n::text(ui.ctx(),t.label())).clicked() {
                    self.tab = t;
                }
            }
        });
    }

    /// 视图页签悬浮在画布左上角，与右上角缩放/布局工具条处于同一层。
    /// 这样切换 Map、Sequence 或数据面板时不会再占用中央内容的固定顶栏高度。
    fn view_tabs_floating(&mut self, ctx: &egui::Context, rect: egui::Rect, tabs: &[Tab], owner: egui::Id) {
        egui::Area::new(owner.with("view_tabs_floating"))
            .order(egui::Order::Background)
            .fixed_pos(rect.left_top() + egui::vec2(6.0, 6.0))
            .movable(false)
            .show(ctx, |ui| {
                Self::canvas_tool_frame().show(ui, |ui| self.view_tabs(ui, tabs));
            });
    }

    /// 「Open File… / Open Folder…」小窗口（我们无系统文件对话框，用路径输入替代）。
    /// New DNA/RNA/Protein File 面板（SnapGene 8 的 New File pane）。
    fn new_file_window(&mut self, ctx: &egui::Context) {
        // SnapGene 的 New File 面板贴着文档区左上角（不是屏幕居中）
        let rect = self.central_rect;
        let done = match self.new_file.as_mut() {
            Some(st) => crate::newfile::draw(st, ctx, rect),
            None => return,
        };
        match done {
            Some(true) => self.new_file_create(),
            Some(false) => self.new_file_cancel(),
            None => self.sync_new_file_preview(),
        }
    }

    /// 内部「Save As」窗口：只在没有 zenity/kdialog 的环境里使用。
    fn save_dialog_window(&mut self, ctx: &egui::Context) {
        if !self.show_save_dialog {
            return;
        }
        let mut win_open = true;
        let mut cancel = false;
        let mut go = false;
        egui::Window::new(crate::i18n::text(ctx,"Save As"))
            .id(egui::Id::unique("save_dialog"))
            .open(&mut win_open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.set_min_width(420.0);
                let te = ui.add(
                    egui::TextEdit::singleline(&mut self.save_as_path)
                        .desired_width(400.0)
                        .hint_text(crate::i18n::text(ui.ctx(),"/path/to/file.dna")),
                );
                te.request_focus();
                if ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    go = true;
                }
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new(crate::i18n::text(ui.ctx(),"Include the .dna / .rna / .prot extension."))
                            .size(10.5)
                            .color(theme::TEXT_DIM),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() {
                            cancel = true;
                        }
                        if ui.button(crate::i18n::text(ui.ctx(),"Save")).clicked() {
                            go = true;
                        }
                    });
                });
            });
        if go {
            self.confirm_save_dialog();
        } else if cancel || !win_open {
            self.show_save_dialog = false;
            self.status = "Save cancelled".into();
        }
    }

    /// 内部 Save As 窗口点 Save / 回车后的落盘动作。
    fn confirm_save_dialog(&mut self) {
        let t = self.save_as_path.trim().to_string();
        if t.is_empty() {
            self.status = "Type a path to save to (Save As)".into();
            return;
        }
        self.save_path = t;
        self.show_save_dialog = false;
        self.write_save_target(true);
    }

    fn vector_library_menu(&mut self, ui: &mut egui::Ui, category: &str) {
        let entries = crate::clonedlg::vector_library_entries(category);
        if entries.is_empty() {
            ui.label(crate::i18n::text(ui.ctx(),"No local vector templates found."));
            ui.label(crate::i18n::text(ui.ctx(),"Open a suitable .dna file, or configure HELIX_SNAP_SGENE_RESOURCES."));
        }
        for (name, path) in entries {
            if ui.button(&name).on_hover_text(crate::i18n::text(ui.ctx(),"Open an independent template copy; the library file is not overwritten.")).clicked() {
                self.open_vector_template(&path, &name);
                ui.close();
            }
        }
        ui.separator();
        if ui.button(crate::i18n::text(ui.ctx(),"Open Vector File…")).clicked() {
            self.show_open_dialog = Some(OpenKind::File);
            ui.close();
        }
    }

    fn open_vector_template(&mut self, path: &std::path::Path, name: &str) {
        let result = std::fs::read(path).map_err(|e| e.to_string())
            .and_then(|bytes| Document::from_imported(bytes, format!("{name}.dna")));
        match result {
            Ok(doc) => {
                self.push_doc(doc, None);
                self.status = format!("Opened vector template {name}; save to a new file before editing the library.");
            }
            Err(e) => self.status = format!("Cannot open vector template {name}: {e}"),
        }
    }

    fn open_dialog_window(&mut self, ctx: &egui::Context) {
        let Some(kind) = self.show_open_dialog else {
            return;
        };
        let (title, hint) = match kind {
            OpenKind::File => ("Open File", "/path/to/file.dna"),
            OpenKind::Folder => ("Open Folder", "/path/to/folder"),
        };
        let mut win_open = true;
        let mut cancel = false;
        let mut go = false;
        let mut use_system = false;
        egui::Window::new(crate::i18n::text(ctx,title))
            .id(egui::Id::unique("open_dialog"))
            .open(&mut win_open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.set_min_width(420.0);
                let field = match kind {
                    OpenKind::File => &mut self.open_path,
                    OpenKind::Folder => &mut self.folder_input,
                };
                let te = ui.add(
                    egui::TextEdit::singleline(field)
                        .desired_width(400.0)
                        .hint_text(hint),
                );
                te.request_focus();
                if ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    go = true;
                }
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    if crate::sysopen::has_native_dialog() {
                        let tool = crate::sysopen::dialog_tool_name().unwrap_or("");
                        if ui
                            .button(crate::i18n::text(ui.ctx(),"Use System Dialog…"))
                            .on_hover_text(format!("zenity/kdialog ({tool})"))
                            .clicked()
                        {
                            use_system = true;
                        }
                    }
                    ui.label(
                        egui::RichText::new(crate::i18n::text(ui.ctx(),"Tip: you can also drop a file onto the window."))
                            .size(10.5)
                            .color(theme::TEXT_DIM),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() {
                            cancel = true;
                        }
                        if ui.button(crate::i18n::text(ui.ctx(),"Open")).clicked() {
                            go = true;
                        }
                    });
                });
            });
        if use_system {
            match kind {
                OpenKind::File => self.open_file_system_dialog(),
                OpenKind::Folder => self.open_folder_with_dialog(),
            }
            cancel = true;
        }
        if go {
            match kind {
                OpenKind::File => {
                    let t = self.open_path.trim().to_string();
                    if !t.is_empty() {
                        self.open(&PathBuf::from(t));
                    }
                }
                OpenKind::Folder => {
                    let t = self.folder_input.trim().to_string();
                    if !t.is_empty() {
                        let pb = PathBuf::from(&t);
                        if pb.is_dir() {
                            self.set_folder(Some(pb));
                            self.refresh_folder();
                        } else {
                            self.folder_error = Some(format!("not a folder: {t}"));
                        }
                    }
                }
            }
            cancel = true;
        }
        if !win_open || cancel {
            self.show_open_dialog = None;
        }
    }

    /// 导出当前酶组（SnapGene 的 Export Enzyme Set…）：写成 `<文档>.enzymes.txt`。
    pub fn export_enzyme_set(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            self.status = "No file open".into();
            return;
        };
        let Some(base) = doc.path.clone() else {
            self.status = "Export enzyme set: save the file first".into();
            return;
        };
        let mut out = String::new();
        out.push_str(&format!(
            "# AppimageGene enzyme set: {}\n# name\trecognition\tcut\treverse_cut\toverhang\tsites\n",
            self.enzyme_preset.label()
        ));
        for (k, &idx) in self.sites.selected.iter().enumerate() {
            let e = &hx_enzymes::ENZYMES[idx];
            let n = self.sites.hits.get(k).map(|h| h.len()).unwrap_or(0);
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                e.name,
                e.site,
                e.cut,
                e.reverse_cut,
                match e.overhang {
                    hx_enzymes::Overhang::Blunt => "blunt",
                    hx_enzymes::Overhang::FivePrime => "5'",
                    hx_enzymes::Overhang::ThreePrime => "3'",
                    hx_enzymes::Overhang::Iis => "IIS",
                },
                n
            ));
        }
        let target = base.with_extension("enzymes.txt");
        match std::fs::write(&target, out) {
            Ok(()) => {
                self.status = format!(
                    "Exported {} enzyme(s) → {}",
                    self.sites.selected.len(),
                    target.display()
                );
            }
            Err(e) => self.status = format!("Export enzyme set failed: {e}"),
        }
    }

    /// 导出特征数据（SnapGene 的 Export Feature Data…）：写成 `<文档>.features.txt`。
    pub fn export_feature_data(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            self.status = "No file open".into();
            return;
        };
        let Some(base) = doc.path.clone() else {
            self.status = "Export feature data: save the file first".into();
            return;
        };
        let len = doc.len();
        let mut out = String::from("# name\ttype\tstart\tend\tlength\tstrand\tcolor\n");
        for f in doc.features() {
            let seg = f.segments.first();
            let (start, end, color) = match seg {
                Some(s) => (s.range.start, s.range.end, s.color.clone()),
                None => (0, 0, String::new()),
            };
            let span = hx_render::span_len(start, end, len);
            let strand = match f.directionality {
                1 => "forward",
                2 => "reverse",
                3 => "both",
                _ => "none",
            };
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                f.name,
                f.kind,
                start + 1,
                end,
                span,
                strand,
                color
            ));
        }
        let target = base.with_extension("features.txt");
        match std::fs::write(&target, out) {
            Ok(()) => {
                self.status = format!(
                    "Exported {} feature(s) → {}",
                    doc.features().len(),
                    target.display()
                );
            }
            Err(e) => self.status = format!("Export feature data failed: {e}"),
        }
    }

    /// 「Enter Sequence Text…」：粘贴框 + Insert/Replace（原顶部工具条的输入框）。
    fn paste_dialog_window(&mut self, ctx: &egui::Context) {
        if let Some((s, e)) = self.pending_delete {
            let unit = self.doc.as_ref().map(|d| if d.kind() == hx_core::MoleculeKind::Protein { "aa" } else if d.session.flags().strandedness() == hx_core::Strandedness::Single { "bases" } else { "bp" }).unwrap_or("bp");
            egui::Modal::new(egui::Id::unique("delete_sequence_confirmation")).show(ctx, |ui| {
                ui.set_min_width(260.0);
                ui.strong(crate::ui_format!(ui.ctx(),"Delete {} {}?", nfmt(e - s),crate::i18n::text(ui.ctx(),unit)));
                ui.add_space(20.0);
                ui.horizontal(|ui| {
                    if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() { self.pending_delete = None; }
                    if ui.button(crate::i18n::text(ui.ctx(),"OK")).clicked() {
                        self.pending_delete = None;
                        self.apply_edit(Edit::Delete { start: s, end: e }, format!("Deleted {}..{}", s + 1, e), Some(s));
                    }
                });
                if ui.input(|i| i.key_pressed(egui::Key::Escape)) { self.pending_delete = None; }
            });
        }
        if !self.show_paste_dialog {
            return;
        }
        self.sequence_insert_modal(ctx);
    }

    fn sequence_insert_modal(&mut self, ctx: &egui::Context) {
        let Some(doc) = self.doc.as_ref() else { self.show_paste_dialog = false; return; };
        let protein = doc.kind() == hx_core::MoleculeKind::Protein;
        let rna = doc.kind() == hx_core::MoleculeKind::Rna;
        let unit = if protein { "aa" } else if doc.session.flags().strandedness() == hx_core::Strandedness::Single { "bases" } else { "bp" };
        let title = crate::ui_format!(ctx,"Insert {} into {}", crate::i18n::text(ctx,if protein { "Amino Acids" } else { "Bases" }), doc.map_name());
        let Some((s, e)) = self.effective_range() else { self.show_paste_dialog = false; return; };
        egui::Modal::new(egui::Id::unique("insert_sequence_modal")).show(ctx, |ui| {
            ui.set_width(460.0);
            ui.heading(title);
            ui.separator();
            ui.strong(if s < e { crate::ui_format!(ctx,"Replace {} {} .. {}:", crate::i18n::text(ctx,if protein { "amino acids" } else { "bases" }), s + 1, e) }
                else { crate::ui_format!(ctx,"Insert after {} {}:", crate::i18n::text(ctx,if protein { "amino acid" } else { "base" }), s) });
            ui.horizontal(|ui| {
                ui.radio_value(&mut self.insert_unknown, false, crate::i18n::text(ctx,"Insert the following sequence:"));
                ui.label(format!("{} {}", hx_algo::parse_seq_text(&self.paste).len(),crate::i18n::text(ctx,unit)));
                if ui.small_button(crate::i18n::text(ctx,"Clear")).clicked() { self.paste.clear(); }
            });
            ui.add_enabled(!self.insert_unknown, egui::TextEdit::multiline(&mut self.paste).font(egui::TextStyle::Monospace).desired_width(450.0).desired_rows(10));
            if !protein { ui.checkbox(&mut self.insert_reverse, crate::i18n::text(ui.ctx(),"Reverse complement")); }
            ui.horizontal(|ui| {
                ui.radio_value(&mut self.insert_unknown, true, crate::i18n::text(ui.ctx(),"Insert"));
                ui.add(egui::DragValue::new(&mut self.insert_unknown_count).range(1..=1_000_000));
                ui.label(if protein { "X's" } else { "N's" });
            });
            let mut seq = if self.insert_unknown { vec![if protein { b'X' } else { b'N' }; self.insert_unknown_count] } else { hx_algo::parse_seq_text(&self.paste) };
            let valid = seq.iter().all(|b| if protein { b"ACDEFGHIKLMNPQRSTVWYBXZJUO*".contains(&b.to_ascii_uppercase()) } else { b"ACGTURYSWKMBDHVN".contains(&b.to_ascii_uppercase()) });
            if !protein && self.insert_reverse { seq = hx_algo::reverse_complement(&seq); }
            if rna { for b in &mut seq { if *b == b'T' { *b = b'U'; } if *b == b't' { *b = b'u'; } } }
            egui::Frame::new().fill(egui::Color32::from_rgb(210, 250, 212)).inner_margin(10.0).show(ui, |ui| {
                ui.set_min_width(440.0);
            if s < e { ui.label(crate::ui_format!(ui.ctx(),"{} {} will be deleted", e - s,crate::i18n::text(ui.ctx(),unit))); }
            ui.label(crate::ui_format!(ui.ctx(),"{} {} will be inserted", seq.len(),crate::i18n::text(ui.ctx(),unit)));
            });
            if !valid { ui.colored_label(egui::Color32::RED, crate::i18n::text(ui.ctx(),"Sequence contains symbols not valid for this molecule.")); }
            ui.horizontal(|ui| {
                if ui.add_enabled(valid && !seq.is_empty(), egui::Button::new(crate::i18n::text(ui.ctx(),"Insert"))).clicked() {
                    let len = seq.len();
                    let edit = if s < e { Edit::Replace { start: s, end: e, seq } } else { Edit::Insert { at: s, seq } };
                    self.apply_edit(edit, format!("Inserted {len} {unit}"), Some(s + len));
                    self.selection = Some((s, s + len));
                    self.show_paste_dialog = false;
                }
                if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() { self.show_paste_dialog = false; }
            });
            if ui.input(|i| i.key_pressed(egui::Key::Escape)) { self.show_paste_dialog = false; }
        });
    }

    #[allow(dead_code)]
    fn legacy_paste_dialog_window(&mut self, ctx: &egui::Context) {
        let mut open = true;
        egui::Window::new(crate::i18n::text(ctx,"Enter Sequence Text"))
            .id(egui::Id::unique("paste_dialog"))
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.set_min_width(340.0);
                ui.add(
                    egui::TextEdit::multiline(&mut self.paste)
                        .desired_rows(3)
                        .desired_width(330.0)
                        .hint_text(crate::i18n::text(ui.ctx(),"type/paste bases (whitespace ignored)")),
                );
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new(crate::i18n::text(ui.ctx(),"Insert at cursor / replace selection"))
                            .size(11.0)
                            .color(theme::TEXT_DIM),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button(crate::i18n::text(ui.ctx(),"Insert")).clicked() {
                            self.do_insert();
                        }
                        if ui.button(crate::i18n::text(ui.ctx(),"Replace")).clicked() {
                            self.do_replace();
                        }
                    });
                });
            });
        if !open {
            self.show_paste_dialog = false;
        }
    }

    /// SnapGene 的序列视图右键菜单（SequenceContextMenu 的子集）。
    pub fn seq_context_menu(&mut self, ui: &mut egui::Ui) {
        let has_sel = self.effective_range().is_some();
        if item(ui, "Cut", Some("Ctrl+X"), has_sel) {
            self.do_cut();
            ui.close();
        }
        if item(ui, "Copy", Some("Ctrl+C"), has_sel) {
            self.do_copy();
            ui.close();
        }
        if item(ui, "Paste", Some("Ctrl+V"), has_sel) {
            self.paste_buffer_into_seq();
            ui.close();
        }
        if item(ui, "Delete", Some("Del"), has_sel) {
            self.do_delete();
            ui.close();
        }
        ui.separator();
        if item(ui, "Copy Reverse Complement", None, has_sel) {
            if let Some((s, e)) = self.effective_range() {
                if let Some(doc) = self.doc.as_ref() {
                    let seq = doc.sequence();
                    let mut sub =
                        hx_algo::reverse_complement(&seq[s.min(seq.len())..e.min(seq.len())]);
                    sub.reverse();
                    let text = String::from_utf8_lossy(&sub).to_string();
                    self.copy_text(text);
                }
            }
            ui.close();
        }
        if item(ui, "Make Uppercase", None, has_sel) {
            self.do_change_case(true);
            ui.close();
        }
        if item(ui, "Make Lowercase", None, has_sel) {
            self.do_change_case(false);
            ui.close();
        }
        ui.separator();
        // 光标落在酶位点上时：选中识别序列
        let site_sel = self.cursor_bp.and_then(|bp| {
            let doc = self.doc.as_ref()?;
            let len = doc.len();
            for (k, &eidx) in self.sites.selected.iter().enumerate() {
                let e = hx_enzymes::ENZYMES.get(eidx)?;
                for &pos in self.sites.hits.get(k).map(|h| h.as_slice()).unwrap_or(&[]) {
                    let pos = pos as usize;
                    if pos == bp {
                        let end = (pos + e.site_len()).min(len);
                        return Some((pos, end));
                    }
                }
            }
            None
        });
        if item(ui, "Select Recognition Sequence", None, site_sel.is_some()) {
            if let Some((a, b)) = site_sel {
                self.selection = Some((a, b));
                self.selected_residue = None;
                self.status = format!("Selected {} .. {}", nfmt(a + 1), nfmt(b));
            }
            ui.close();
        }
        if item(ui, "Edit DNA Ends…", None, false) {}
        ui.separator();
        if item(ui, "Find DNA Sequence", Some("Ctrl+F"), false) {}
        if item(ui, "Export Map…", None, self.doc.is_some()) {
            self.tab = Tab::Map;
            self.do_export_svg();
            ui.close();
        }
    }

    /// 「View Options…」（≡ Display 下拉的内容，作为可停靠窗口呈现）。
    fn view_options_window(&mut self, ctx: &egui::Context) {
        if !self.show_view_options {
            return;
        }
        let mut open = true;
        egui::Window::new(crate::i18n::text(ctx,"View Options"))
            .id(egui::Id::unique("view_options"))
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.set_min_width(280.0);
                self.display_menu(ui);
            });
        if !open {
            self.show_view_options = false;
        }
    }

    /// 光标/悬停位置上的特征或酶（状态栏中段，SnapGene 的 `KanR - 560 - 1,375 (816 bp)`）。
    fn hover_info(&self) -> Option<String> {
        let doc = self.doc.as_ref()?;
        let bp = self.hover_bp.or(self.cursor_bp)?;
        if bp >= doc.len() {
            return None;
        }
        let mut best: Option<(usize, String, usize, usize)> = None;
        for f in doc.features() {
            for s in &f.segments {
                let (a, b) = (s.range.start, s.range.end);
                if a <= bp && bp < b {
                    let span = b - a;
                    if best.as_ref().map(|(sp, ..)| span < *sp).unwrap_or(true) {
                        best = Some((span, f.name.clone(), a, b));
                    }
                }
            }
        }
        if let Some((_, name, a, b)) = best {
            return Some(format!(
                "{name} - {} - {} ({} {})",
                nfmt(a + 1),
                nfmt(b),
                nfmt(b - a),
                doc.unit()
            ));
        }
        for (k, &eidx) in self.sites.selected.iter().enumerate() {
            let Some(e) = hx_enzymes::ENZYMES.get(eidx) else {
                continue;
            };
            for &pos in self.sites.hits.get(k).map(|h| h.as_slice()).unwrap_or(&[]) {
                if pos == bp {
                    return Some(e.name.to_string());
                }
            }
        }
        None
    }

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.add_space(6.0);
            if let Some(v) = self.viewer.as_ref() {
                ui.label(
                    egui::RichText::new(format!("{} · {}", v.title(), crate::i18n::viewer_summary(ui.ctx(),&v.summary())))
                        .size(11.5)
                        .color(theme::TEXT_DIM),
                );
                ui.separator();
            } else if let Some(doc) = self.doc.as_ref() {
                let head = if let Some(class) = doc.native_class() {
                    // 比对/凝胶：没有序列可显示，写类别与字节数
                    let size = doc
                        .native_only
                        .as_ref()
                        .map(|n| n.bytes.len())
                        .unwrap_or(0);
                    format!(
                        "{} · SnapGene {class} document · {} bytes",
                        doc.display_name(),
                        nfmt(size)
                    )
                } else if doc.kind().is_nucleotide() {
                    format!(
                        "{} · {} {} · {}",
                        doc.display_name(),
                        nfmt(doc.len()),
                        doc.unit(),
                        doc.topo_label()
                    )
                } else {
                    format!(
                        "{} · {} {}",
                        doc.display_name(),
                        nfmt(doc.len()),
                        doc.unit()
                    )
                };
                ui.label(
                    egui::RichText::new(head)
                        .size(11.5)
                        .color(theme::TEXT_DIM),
                );
                ui.separator();
                // SnapGene 底部在选中特征时写
                //   蛋白  `Selected:  domain interfaces [active] (p3 .. p414 / 39 aa)`
                //   核苷酸 `Selected:  ori (867 .. 1455 = 589 bp)   [54% GC]`
                // （oracle `/tmp/opencode/sg/st-prot.png`、`dn1-status-z.png`：
                //  核苷酸用 `=` 分隔长度并追加选中区的 GC 含量，蛋白用 `/` 且无 GC）。
                if let Some(f) = self
                    .selected_feature
                    .and_then(|i| doc.features().get(i))
                {
                    let segs: Vec<(usize, usize)> = f
                        .segments
                        .iter()
                        .map(|s| (s.range.start, s.range.end))
                        .collect();
                    let (a, b, n, _) = hx_render::feature_span(&segs, doc.len());
                    let pp = doc.kind().position_prefix();
                    let text = if doc.kind().is_nucleotide() {
                        format!(
                            "Selected: {}  ({} .. {} = {} {})   [{}% GC]",
                            f.name,
                            nfmt(a + 1),
                            nfmt(b),
                            nfmt(n),
                            doc.unit(),
                            gc_percent(doc.sequence(), a, b, doc.len())
                        )
                    } else {
                        format!(
                            "Selected: {}  ({}{} .. {}{} / {} {})",
                            f.name,
                            pp,
                            nfmt(a + 1),
                            pp,
                            nfmt(b),
                            nfmt(n),
                            doc.unit()
                        )
                    };
                    ui.label(egui::RichText::new(text).size(11.5).color(theme::TEXT));
                    ui.separator();
                }
                let pos = match (self.selection, self.cursor_bp, self.hover_bp) {
                    (Some((s, e)), _, _) => crate::ui_format!(ui.ctx(),
                        "Selection {}..{} ({} {})",
                        nfmt(s + 1),
                        nfmt(e),
                        nfmt(hx_render::span_len(s, e, doc.len())),
                        crate::i18n::text(ui.ctx(),doc.unit())
                    ),
                    (None, Some(c), _) => crate::ui_format!(ui.ctx(),"Cursor:{}", nfmt(c + 1)),
                    (None, None, Some(h)) => crate::ui_format!(ui.ctx(),"Position {}", nfmt(h + 1)),
                    _ => String::new(),
                };
                if !pos.is_empty() {
                    ui.label(egui::RichText::new(pos).size(11.5).color(theme::TEXT));
                }
                // SnapGene 在这里显示光标/悬停位置上的特征（`Name - start - end (N bp)`）
                if let Some(info) = self.hover_info() {
                    ui.separator();
                    ui.label(egui::RichText::new(info).size(11.5).color(theme::TEXT));
                }
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add_space(6.0);
                let native = self
                    .doc
                    .as_ref()
                    .is_some_and(|d| d.native_class().is_some());
                let preset = match self.doc.as_ref().map(|d| d.kind()) {
                    Some(hx_core::MoleculeKind::Dna) if !native => crate::ui_format!(ui.ctx(),
                        "{} · {} site(s)",
                        crate::i18n::text(ui.ctx(),self.enzyme_preset.label()),
                        nfmt(self.sites.total())
                    ),
                    Some(hx_core::MoleculeKind::Dna) => crate::i18n::text(ui.ctx(),"no sequence view").to_string(),
                    Some(hx_core::MoleculeKind::Rna) => crate::i18n::text(ui.ctx(),"RNA — enzymes n/a").to_string(),
                    Some(hx_core::MoleculeKind::Protein) => crate::i18n::text(ui.ctx(),"Protein — enzymes n/a").to_string(),
                    None => String::new(),
                };
                ui.label(
                    egui::RichText::new(preset)
                        .size(11.5)
                        .color(theme::TEXT_DIM),
                );
                ui.separator();
                ui.label(
                    egui::RichText::new(crate::i18n::message(ui.ctx(),&self.status))
                        .size(11.5)
                        .color(theme::ACCENT),
                );
            });
        });
    }

    /// 左侧竖向显示工具条（SnapGene 顺序：酶、特征、引物、翻译/ORF、颜色）。
    fn tool_palette(&mut self, ctx: &egui::Context, rect: egui::Rect) {
        ctx.data_mut(|d|d.insert_temp(egui::Id::unique("enzyme_render_preferences"),self.preferences.clone()));
        egui::Area::new(egui::Id::unique("canvas_palette"))
            .order(egui::Order::Background)
            .fixed_pos(rect.left_top() + egui::vec2(6.0, 62.0))
            .movable(false)
            .show(ctx, |ui| {
                egui::Frame::new().show(ui, |ui| {
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 3.0;
                        let has_doc = self.doc.is_some();
                        let hover_id = egui::Id::unique("palette_hover_row");
                        let pointer = ctx.pointer_hover_pos();
                        let hover_row = ctx.data(|d| d.get_temp::<(usize, egui::Rect)>(hover_id)).filter(|(_, r)| pointer.is_some_and(|p| r.contains(p))).map(|(i, _)| i);
                        let mut next_hover = None;
                        for (index, (icon, active, tip, act)) in self.palette_items().into_iter().enumerate() {
                            let row = ui.horizontal(|ui| {
                            ui.spacing_mut().item_spacing.x = 2.0;
                            let response = Self::canvas_tool_frame().inner_margin(egui::Margin::same(2)).show(ui, |ui| icons::palette(ui, icon, active && icon!=icons::Icon::Colors,
                                has_doc && (icon!=icons::Icon::Translation || !self.is_protein()))
                                .on_hover_text(crate::i18n::text(ui.ctx(),tip))).inner;
                            if response.clicked() {
                                act.apply(self);
                            }
                            // The arrow has its own hit target: opening options
                            // must not toggle the adjacent visibility button.
                            if has_doc && (response.hovered() || hover_row == Some(index) || egui::Popup::is_any_open(ctx)) && matches!(icon, icons::Icon::Enzymes | icons::Icon::Features | icons::Icon::Primers | icons::Icon::Translation | icons::Icon::Colors) {
                                let (arrow_rect, arrow) = ui.allocate_exact_size(egui::vec2(16.0, 30.0), egui::Sense::click());
                                ui.painter().rect_filled(arrow_rect, 3.0, egui::Color32::WHITE);
                                ui.painter().rect_stroke(arrow_rect, 3.0, egui::Stroke::new(0.7, theme::BORDER), egui::StrokeKind::Inside);
                                let c = arrow_rect.center();
                                ui.painter().add(egui::Shape::line(vec![c + egui::vec2(-3.0,-1.5), c + egui::vec2(0.0,1.5), c + egui::vec2(3.0,-1.5)], egui::Stroke::new(1.0, theme::TEXT)));
                                egui::Popup::menu(&arrow).show(|ui| { self.palette_options(ui, icon); });
                            }
                            });
                            let hit = row.response.rect.expand2(egui::vec2(22.0, 1.0));
                            if pointer.is_some_and(|p| hit.contains(p)) { next_hover = Some((index, hit)); }
                        }
                        if let Some(value) = next_hover { ctx.data_mut(|d| d.insert_temp(hover_id, value)); ctx.request_repaint(); }
                    });
                });
            });
    }

    /// 左侧竖排工具条的按钮表（**按视图与分子类型**决定内容）。
    ///
    /// oracle 实测（`/tmp/opencode/sg/hvp-*.png`、`hvd-*.png`、`hv-rna-*`，鼠标悬停取 tooltip）：
    /// - DNA 图谱/序列：`Show enzymes` / `Show features` / `Show primers` /
    ///   `Show translations and ORFs` / `Show colors`（图谱另有 `Show locations`、`Show alignments`，
    ///   V1 暂无对应功能）；
    /// - **RNA**：`Show features` / `Show translations and ORFs` / `Show colors`（没有酶和引物）；
    /// - **蛋白**：`Show "Region" features` / `Show "Site" and "Het" features` /
    ///   `Show "Bond" features` / `Show colors`（图谱再加 `Show translation numbering`）。
    fn palette_items(&self) -> Vec<(icons::Icon, bool, &'static str, PalAction)> {
        use icons::Icon;
        let protein = self.is_protein();
        let dna = self.is_dna();
        let mut items: Vec<(Icon, bool, &'static str, PalAction)> = Vec::new();
        if protein {
            items.push((Icon::Region, self.show_regions, r#"Show "Region" features"#, PalAction::Regions));
            items.push((Icon::Sites, self.show_sites_het, r#"Show "Site" and "Het" features"#, PalAction::SitesHet));
            items.push((Icon::Bonds, self.show_bonds, r#"Show "Bond" features"#, PalAction::Bonds));
        } else {
            if dna {
                items.push((Icon::Enzymes, self.show_sites, "Show enzymes", PalAction::Sites));
            }
            items.push((Icon::Features, self.show_features, "Show features", PalAction::Features));
            if dna {
                items.push((Icon::Primers, self.show_primers, "Show primers", PalAction::Primers));
            }
        }
        // 蛋白没有「翻译与 ORF」这一档（oracle：SG 蛋白视图左侧 4 个按钮 + 图谱的 `(123)`）。
        if !protein {
            items.push((
                Icon::Translation,
                self.show_translation,
                "Show translations and ORFs",
                PalAction::Translation,
            ));
        }
        items.push((Icon::Colors, self.show_colors, "Show colors", PalAction::Colors));
        if protein && self.tab == Tab::Map {
            items.push((
                Icon::Numbers,
                self.show_trans_numbering,
                "Show translation numbering",
                PalAction::TransNumbering,
            ));
        }
        items
    }

    fn palette_options(&mut self, ui: &mut egui::Ui, icon: icons::Icon) {
        use icons::Icon;
        match icon {
            Icon::Enzymes => {
                self.sg_enzyme_options(ui);
            }
            Icon::Primers => {
                if ui.button(crate::i18n::text(ui.ctx(),"Hybridization Parameters…")).clicked() {
                    self.hybridization_draft = self.doc.as_ref().map(|d| hx_formats::annot::parse_hybridization_params(d.block_payload(5).unwrap_or_default()));
                    ui.close();
                }
            }
            Icon::Features => { ui.checkbox(&mut self.show_features, crate::i18n::text(ui.ctx(),"Show features")); if ui.button(crate::i18n::text(ui.ctx(),"Choose Features…")).clicked() { self.tab = Tab::Features; ui.close(); } }
            Icon::Translation => {
                if ui.selectable_label(self.translation_frames == 63, crate::i18n::text(ui.ctx(),"ORFs Only")).clicked() { self.translation_frames = 63; self.show_translation = true; ui.close(); }
                ui.separator();
                for (label, mask) in [("Frame +1 Only",1),("Frame +2 Only",2),("Frame +3 Only",4),("Frame -1 Only",8),("Frame -2 Only",16),("Frame -3 Only",32),("Top 3 Frames",7),("Bottom 3 Frames",56),("All 6 Frames",63)] {
                    if ui.selectable_label(self.translation_frames == mask, crate::i18n::text(ui.ctx(),label)).clicked() { self.translation_frames = mask; self.show_translation = true; ui.close(); }
                }
                ui.separator();
                if ui.button(crate::i18n::text(ui.ctx(),"Translation Options…")).clicked() { self.show_translation_options = true; ui.close(); }
            }
            Icon::Colors => {
                ui.checkbox(&mut self.show_colors, crate::i18n::text(ui.ctx(),"Show annotation colors"));
                ui.separator();
                for (label, mode) in [("Neutral sequence colors",0),("GC Content",1),("Purine vs Pyrimidine",2),("Nucleotide Colors",3)] {
                    if ui.selectable_label(self.sequence_color_mode == mode, crate::i18n::text(ui.ctx(),label)).clicked() { self.sequence_color_mode = mode; self.show_colors = true; ui.close(); }
                }
                ui.label(crate::i18n::text(ui.ctx(),"Sequence schemes apply in Sequence view. GC mode colors G/C versus A/T; it is not a sliding-window plot."));
            }
            _ => {}
        }
    }

    fn apply_enzyme_options(&mut self) {
        if let Some(doc) = self.doc.as_ref() {
            self.sites.selected = hx_enzymes::sets::select_ordered(self.enzyme_parent, self.enzyme_filter, doc.sequence(), doc.circular(), &self.preferences.suppliers);
            self.enzyme_preset = EnzymePreset::Custom;
            self.sites_dirty = true; self.show_sites = true;
        }
    }

    fn sg_enzyme_options(&mut self, ui: &mut egui::Ui) {
        use hx_enzymes::sets::{Parent, Filter};
        ui.set_min_width(320.0);
        let mut parents=vec![("All Commercial".to_owned(),Parent::All),("Nonredundant Commercial".to_owned(),Parent::Nonredundant)];
        parents.extend(self.preferences.suppliers.iter().take(self.preferences.supplier_menus).map(|&id|(hx_enzymes::sets::SUPPLIER_NAMES[id as usize].to_owned(),Parent::Supplier(id))));
        for (label, parent) in parents {
            let count = hx_enzymes::sets::parent_indices_ordered(parent, &self.preferences.suppliers).len();
            let shown=if matches!(parent,Parent::All|Parent::Nonredundant){crate::i18n::text(ui.ctx(),&label)}else{&label};
            if ui.selectable_label(self.enzyme_parent == parent && self.enzyme_filter == Filter::All, format!("{shown} ({count})")).clicked() {
                self.enzyme_parent = parent; self.enzyme_filter = Filter::All; self.apply_enzyme_options(); ui.close();
            }
        }
        ui.separator();
        for (label, filter) in [("Unique Cutters",Filter::Unique),("Unique & Dual Cutters",Filter::UniqueDual),("6+ Cutters",Filter::SixPlus),("Unique 6+ Cutters",Filter::UniqueSixPlus),("Type IIS Enzymes",Filter::TypeIis),("Golden Gate Enzymes",Filter::GoldenGate),("Nicking Endonucleases",Filter::Nicking)] {
            let count = self.doc.as_ref().map(|d| hx_enzymes::sets::select_ordered(self.enzyme_parent,filter,d.sequence(),d.circular(),&self.preferences.suppliers).len()).unwrap_or(0);
            if ui.selectable_label(self.enzyme_filter == filter, format!("{} ({count})",crate::i18n::text(ui.ctx(),label))).clicked() { self.enzyme_filter = filter; self.apply_enzyme_options(); ui.close(); }
        }
        ui.separator();
        ui.label(crate::ui_format!(ui.ctx(),"Parental enzyme set: {}", crate::i18n::text(ui.ctx(),match self.enzyme_parent { hx_enzymes::sets::Parent::All => "All Commercial", hx_enzymes::sets::Parent::Nonredundant => "Nonredundant Commercial", hx_enzymes::sets::Parent::Neb => "New England Biolabs", hx_enzymes::sets::Parent::Thermo => "Thermo Fisher (Fermentas)",hx_enzymes::sets::Parent::Supplier(id)=>hx_enzymes::sets::SUPPLIER_NAMES[id as usize] })));
        if ui.button(crate::i18n::text(ui.ctx(),"Choose Enzymes…                 Ctrl+E")).clicked() { self.tab = Tab::Enzymes; ui.close(); }
        if ui.button(crate::i18n::text(ui.ctx(),"Save Enzyme Set…")).clicked() { self.export_enzyme_set(); ui.close(); }
        if ui.button(crate::i18n::text(ui.ctx(),"Enzymes Preferences…")).clicked() { self.preferences.page = 9; self.enzyme_preferences_open = true; ui.close(); }
    }

    fn display_options_dialogs(&mut self, ctx: &egui::Context) {
        self.enzyme_manager.show(ctx);
        if let Some(mut dialog)=self.rename_file.take(){
            match dialog.show(ctx){
                Some(false)=>{},
                Some(true)=>match crate::rename::rename_file(&dialog.path,&dialog.name){
                    Ok(target)=>{self.relocate_open_file(&dialog.path,&target);self.status=format!("Renamed to {}",dialog.name);},
                    Err(e)=>{dialog.error=e;self.rename_file=Some(dialog);}
                },
                None=>self.rename_file=Some(dialog),
            }
        }
        if self.enzyme_preferences_open {
            let old = self.preferences.clone();
            self.enzyme_preferences_open = self.preferences.show(ctx);
            if old != self.preferences {
                if old.suppliers != self.preferences.suppliers || old.parent != self.preferences.parent {
                    self.enzyme_parent = self.preferences.parent;
                    self.apply_enzyme_options();
                }
                self.sites_as_lines = self.preferences.enzyme_lines;
                if let Err(e) = self.preferences.save() { self.status = format!("Cannot save preferences: {e}"); }
            }
        }
        if let Some(mut draft) = self.hybridization_draft {
            let mut close = false;
            egui::Modal::new(egui::Id::unique("hybridization_options")).show(ctx, |ui| {
                ui.set_width(520.0); ui.heading(crate::i18n::text(ui.ctx(),"Hybridization Parameters"));
                ui.strong(crate::i18n::text(ui.ctx(),"Hybridization parameters control how many binding sites are displayed."));
                ui.separator();
                ui.horizontal(|ui| { ui.label(crate::i18n::text(ui.ctx(),"Match at least")); ui.add(egui::DragValue::new(&mut draft.min_continuous_match_len).range(1..=100)); ui.label(crate::i18n::text(ui.ctx(),"bases at the 3′ end of the primer")); });
                ui.checkbox(&mut draft.allow_mismatch, crate::i18n::text(ui.ctx(),"Allow a single isolated mismatch"));
                ui.horizontal(|ui| { ui.label(crate::i18n::text(ui.ctx(),"Minimum melting temperature (°C):")); ui.add(egui::DragValue::new(&mut draft.min_melting_temperature).range(-100.0..=150.0)); });
                ui.separator();
                ui.checkbox(&mut draft.show_additional_five_prime_matches, crate::i18n::text(ui.ctx(),"Show additional binding sites matching the 5′ end"));
                ui.horizontal(|ui| { ui.label(crate::i18n::text(ui.ctx(),"Minimum 5′ match length:")); ui.add(egui::DragValue::new(&mut draft.minimum_five_prime_annealing).range(1..=100)); });
                ui.small(crate::i18n::text(ui.ctx(),"Recalculates ungapped binding sites. Tm uses the existing nearest-neighbour screening estimate; mismatch/bulge thermodynamics are not modeled."));
                ui.horizontal(|ui| {
                    if ui.button(crate::i18n::text(ui.ctx(),"OK")).clicked() {
                        if let Some(doc) = self.doc.as_mut() {
                            match doc.session.recalculate_hybridization(draft) {
                                Ok(()) => self.status = "Hybridization parameters updated".into(), Err(e) => self.status = e,
                            }
                        }
                        close = true;
                    }
                    if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() { close = true; }
                });
                if ui.input(|i| i.key_pressed(egui::Key::Escape)) { close = true; }
            });
            self.hybridization_draft = if close { None } else { Some(draft) };
        }
        if self.show_translation_options {
            egui::Modal::new(egui::Id::unique("translation_options")).show(ctx, |ui| {
                ui.heading(crate::i18n::text(ui.ctx(),"Translation Options"));
                ui.horizontal(|ui| { ui.label(crate::i18n::text(ui.ctx(),"Minimum ORF length (bases):")); ui.add(egui::DragValue::new(&mut self.orf_min_bp).range(3..=100000)); });
                ui.checkbox(&mut self.aa_one_letter, crate::i18n::text(ui.ctx(),"Use 1-letter amino acid codes"));
                ui.label(crate::i18n::text(ui.ctx(),"Frame filters currently apply to ORFs, not continuous six-frame translation."));
                if ui.button(crate::i18n::text(ui.ctx(),"Close")).clicked() { self.show_translation_options = false; }
            });
        }
    }

    fn canvas_tool_frame() -> egui::Frame {
        egui::Frame::new()
            .fill(egui::Color32::WHITE)
            .stroke(egui::Stroke::new(0.7, theme::BORDER))
            .corner_radius(4.0)
            .inner_margin(egui::Margin::same(4))
            .shadow(egui::epaint::Shadow {
                offset: [0, 1],
                blur: 6,
                spread: 0,
                color: egui::Color32::from_black_alpha(30),
            })
    }

    /// 画布右上角悬浮的布局与缩放控制；显示开关统一放在左侧。
    fn canvas_controls(&mut self, ctx: &egui::Context, rect: egui::Rect) {
        if self.tab == Tab::SecondaryStructure {
            self.struct_controls(ctx, rect);
            return;
        }
        let pos = egui::Pos2::new(rect.right() - 6.0, rect.top() + 6.0);
        egui::Area::new(egui::Id::unique("canvas_controls"))
            .order(egui::Order::Background)
            .fixed_pos(pos)
            .pivot(egui::Align2::RIGHT_TOP)
            .movable(false)
            .show(ctx, |ui| {
                Self::canvas_tool_frame().show(ui, |ui| {
                    ui.spacing_mut().item_spacing.x = 2.0;
                    ui.horizontal(|ui| {
                        if self.tab==Tab::Map && icons::palette(ui,icons::Icon::Fit,false,true).on_hover_text(crate::i18n::text(ui.ctx(),"Zoom to fit")).clicked(){self.fit();}
                        // 图谱拓扑与序列行宽设置；显示开关统一放在左侧。
                        // 蛋白图谱没有环状布局（SnapGene 也不给这两个按钮）。
                        if self.tab == Tab::Map && !self.is_protein() {
                            for (l, icon, tip) in [
                                (
                                    hx_render::MapLayout::Circular,
                                    icons::Icon::Circular,
                                    "Show as circular",
                                ),
                                (
                                    hx_render::MapLayout::Linear,
                                    icons::Icon::Linear,
                                    "Show as linear",
                                ),
                            ] {
                                let ring=self.doc.as_ref().is_some_and(|d|self.map_layout.is_ring(d.circular(),d.len()));
                                let on = if l==hx_render::MapLayout::Circular {ring}else{!ring};
                                if icons::palette(ui, icon, on, true)
                                    .on_hover_text(crate::i18n::text(ui.ctx(),tip))
                                    .clicked()
                                {
                                    self.map_layout = l;
                                    self.fit();
                                }
                            }
                        }
                        if self.tab == Tab::Sequence {
                            // SnapGene 的两个行宽按钮：换行到窗口 / 使用固定线宽 ▾
                            let fit_mode = self.seq_line_bp.is_none();
                            if icons::palette(ui, icons::Icon::Fit, fit_mode, true)
                                .on_hover_text(crate::i18n::text(ui.ctx(),"Wrap to window"))
                                .clicked()
                            {
                                self.seq_line_bp = None;
                                self.seq_pan_x = 0.0;
                            }
                            let tip = match self.seq_line_bp {
                                None => "Line width: wrap to window".to_owned(),
                                Some(n) => format!("Line width: {n} bp"),
                            };
                            let line_width =
                                icons::palette(ui, icons::Icon::LineWidth, !fit_mode, true)
                                    .on_hover_text(crate::i18n::message(ui.ctx(),&tip));
                            egui::Popup::menu(&line_width).show(|ui| {
                                for n in [30usize, 50, 75, 100, 150, 300, 500] {
                                    if ui
                                        .selectable_label(
                                            self.seq_line_bp == Some(n),
                                            format!("{n} bp"),
                                        )
                                        .clicked()
                                    {
                                        self.seq_line_bp = Some(n);
                                        self.seq_pan_x = 0.0;
                                        ui.close();
                                    }
                                }
                                if ui
                                    .selectable_label(self.seq_line_bp.is_none(), "Wrap to window")
                                    .clicked()
                                {
                                    self.seq_line_bp = None;
                                    self.seq_pan_x = 0.0;
                                    ui.close();
                                }
                            });
                            ui.separator();
                        }
                        if self.tab!=Tab::Map && icons::palette(ui, icons::Icon::ZoomIn, false, true)
                            .on_hover_text(crate::i18n::text(ui.ctx(),"Zoom in"))
                            .clicked()
                        {
                            self.zoom_in();
                        }
                        if self.tab!=Tab::Map && icons::palette(ui, icons::Icon::ZoomOut, false, true)
                            .on_hover_text(crate::i18n::text(ui.ctx(),"Zoom out"))
                            .clicked()
                        {
                            self.zoom_out();
                        }
                        let display=struct_button(ui,"Display",theme::TEXT);
                        egui::Popup::menu(&display).show(|ui|{
                            self.display_menu(ui);
                            ui.separator();
                            if ui.button(crate::i18n::text(ui.ctx(),"Zoom in")).clicked(){self.zoom_in();ui.close();}
                            if ui.button(crate::i18n::text(ui.ctx(),"Zoom out")).clicked(){self.zoom_out();ui.close();}
                            if ui.button(crate::i18n::text(ui.ctx(),"Export current view… (SVG)")).clicked(){self.do_export_svg();ui.close();}
                        });
                    });
                });
            });
    }

    fn canvas_search_button(&mut self,ctx:&egui::Context,rect:egui::Rect) {
        egui::Area::new(egui::Id::unique("canvas_search"))
            .order(egui::Order::Background).fixed_pos(egui::pos2(rect.left()+6.0,rect.bottom()-42.0))
            .movable(false).show(ctx,|ui|{
                Self::canvas_tool_frame().show(ui,|ui|{
                    icons::palette(ui,icons::Icon::Search,false,false).on_hover_text(crate::i18n::text(ui.ctx(),"Find Sequence — not implemented yet"));
                });
            });
    }

    /// `Secondary Structure` 页的右上工具条（SnapGene：`Export as… ▾` / 缩放组合框 /
    /// Zoom to fit / `Display ▾`）。
    fn struct_controls(&mut self, ctx: &egui::Context, rect: egui::Rect) {
        let pos = egui::Pos2::new(rect.right() - 10.0, rect.top() + 6.0);
        let mut export_img = false;
        let mut export_dot = false;
        let mut reset_all = false;
        egui::Area::new(egui::Id::unique("struct_controls"))
            .order(egui::Order::Background)
            .fixed_pos(pos)
            .pivot(egui::Align2::RIGHT_TOP)
            .movable(false)
            .show(ctx, |ui| {
                Self::canvas_tool_frame().show(ui, |ui| {
                    ui.spacing_mut().item_spacing.x = 6.0;
                    ui.horizontal(|ui| {
                        let r = struct_button(ui, "Export as…", theme::ACCENT);
                        egui::Popup::menu(&r).show(|ui| {
                            if ui.button(crate::i18n::text(ui.ctx(),"Image…")).clicked() {
                                export_img = true;
                                ui.close();
                            }
                            if ui.button(crate::i18n::text(ui.ctx(),"Dot-Bracket…")).clicked() {
                                export_dot = true;
                                ui.close();
                            }
                        });
                        ui.separator();
                        // 缩放组合框：SnapGene 显示四位有效数字（`68.46%`）
                        let shown = if self.struct_zoom > 0.0 {
                            self.struct_zoom
                        } else {
                            self.struct_fit
                        };
                        ui.spacing_mut().item_spacing.x = 4.0;
                        egui::ComboBox::from_id_salt("struct_zoom")
                            .width(126.0)
                            .selected_text(format!("{:.2}%", shown * 100.0))
                            .show_ui(ui, |ui| {
                                let cur = if self.struct_zoom > 0.0 {
                                    self.struct_zoom
                                } else {
                                    self.struct_fit
                                };
                                let mut want: Option<f32> = None;
                                if ui.button(crate::i18n::text(ui.ctx(),"Zoom to fit")).clicked() {
                                    want = Some(0.0);
                                    ui.close();
                                }
                                for p in [0.25, 0.5, 0.75, 1.0] {
                                    if ui
                                        .selectable_label(
                                            self.struct_zoom > 0.0
                                                && (self.struct_zoom - p).abs() < 1e-4,
                                            format!("{:.0}%", p * 100.0),
                                        )
                                        .clicked()
                                    {
                                        want = Some(p);
                                        ui.close();
                                    }
                                }
                                ui.separator();
                                if ui.button(crate::i18n::text(ui.ctx(),"Increase base size")).clicked() {
                                    want = Some(cur * 1.15);
                                    ui.close();
                                }
                                if ui.button(crate::i18n::text(ui.ctx(),"Decrease base size")).clicked() {
                                    want = Some(cur / 1.15);
                                    ui.close();
                                }
                                if ui.button(crate::i18n::text(ui.ctx(),"Reset base size")).clicked() {
                                    want = Some(0.0);
                                    ui.close();
                                }
                                ui.separator();
                                if ui.button(crate::i18n::text(ui.ctx(),"Reset current display")).clicked() {
                                    want = Some(0.0);
                                    ui.close();
                                }
                                if ui.button(crate::i18n::text(ui.ctx(),"Reset all display")).clicked() {
                                    want = Some(0.0);
                                    reset_all = true;
                                    ui.close();
                                }
                                if let Some(z) = want {
                                    self.struct_zoom = z;
                                    self.struct_pan = egui::Vec2::ZERO;
                                }
                            });
                        if icons::palette(ui, icons::Icon::Fit, false, true)
                            .on_hover_text(crate::i18n::text(ui.ctx(),"Zoom to fit"))
                            .clicked()
                        {
                            self.struct_zoom = 0.0;
                            self.struct_pan = egui::Vec2::ZERO;
                        }
                        ui.separator();
                        let r = struct_button(ui, "Display", theme::TEXT);
                        egui::Popup::menu(&r).show(|ui| {
                            ui.checkbox(
                                &mut self.struct_show_index,
                                "Show base index labels on structure",
                            );
                        });
                        if icons::palette(ui, icons::Icon::Copy, false, true)
                            .on_hover_text(crate::i18n::text(ui.ctx(),"Copy structure as dot-bracket"))
                            .clicked()
                        {
                            self.do_copy_dot_bracket();
                        }
                    });
                });
            });
        if reset_all {
            self.struct_sel = 0;
            self.struct_sel_base = None;
        }
        if export_img {
            self.do_export_structure_svg();
        }
        if export_dot {
            self.do_export_dot_bracket();
        }
    }

    fn display_menu(&mut self, ui: &mut egui::Ui) {
        // 蛋白的三段特征开关（SnapGene 的 Display 菜单按分子类型换内容）
        if self.is_protein() {
            ui.checkbox(&mut self.show_regions, crate::i18n::text(ui.ctx(),"Show \"Region\" features"));
            ui.checkbox(&mut self.show_sites_het, crate::i18n::text(ui.ctx(),"Show \"Site\" and \"Het\" features"));
            ui.checkbox(&mut self.show_bonds, crate::i18n::text(ui.ctx(),"Show \"Bond\" features"));
            ui.checkbox(&mut self.show_colors, crate::i18n::text(ui.ctx(),"Colors"));
            ui.checkbox(&mut self.show_trans_numbering, crate::i18n::text(ui.ctx(),"Translation numbering"));
            ui.separator();
            ui.checkbox(&mut self.show_translation, crate::i18n::text(ui.ctx(),"Translations and ORFs"));
            ui.checkbox(&mut self.show_ruler, crate::i18n::text(ui.ctx(),"Ruler"));
        } else {
            if self.is_dna() {
                ui.checkbox(&mut self.show_sites, crate::i18n::text(ui.ctx(),"Enzyme sites"));
            }
            ui.checkbox(&mut self.show_features, crate::i18n::text(ui.ctx(),"Features and feature labels"));
            if self.is_dna() {
                ui.checkbox(&mut self.show_primers, crate::i18n::text(ui.ctx(),"Primers"));
            }
            if self.is_dna(){ui.checkbox(&mut self.show_translation,crate::i18n::text(ui.ctx(),"Translations and ORFs"));}
            else{let mut available=false;ui.add_enabled(false,egui::Checkbox::new(&mut available,crate::i18n::text(ui.ctx(),"Translations and ORFs")));}
            ui.checkbox(&mut self.show_colors, crate::i18n::text(ui.ctx(),"Colors"));
            ui.checkbox(&mut self.show_ruler, crate::i18n::text(ui.ctx(),"Ruler"));
        }
        ui.separator();
        ui.label(crate::i18n::text(ui.ctx(),"Map layout"));
        for l in [
            hx_render::MapLayout::Auto,
            hx_render::MapLayout::Circular,
            hx_render::MapLayout::Linear,
        ] {
            // 布局切换会改变缩放语义（线性铺满 vs 环状默认半径），切换后重新适配
            if ui
                .selectable_value(&mut self.map_layout, l, crate::i18n::text(ui.ctx(),l.label()))
                .clicked()
            {
                self.fit();
            }
        }
        ui.separator();
        ui.checkbox(&mut self.show_complement, crate::i18n::text(ui.ctx(),"Complement strand (Sequence)"));
        ui.checkbox(&mut self.aa_one_letter, crate::i18n::text(ui.ctx(),"1-letter amino acids"));
        ui.separator();
        if ui.button(crate::i18n::text(ui.ctx(),"Feature list")).clicked() {
            self.tab = Tab::Features;
            ui.close();
        }
        if ui.add_enabled(self.is_dna(),egui::Button::new(crate::i18n::text(ui.ctx(),"Enzyme list"))).clicked() {
            self.tab = Tab::Enzymes;
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"History")).clicked() {
            self.tab = Tab::History;
            ui.close();
        }
    }

    /// epoch 秒 → `YYYY-MM-DD HH:MM UTC`（Howard Hinnant 的 civil_from_days；不查时区、
    /// 不引依赖）。只用于 About 里显示二进制构建时间，UTC 足够比对先后。
    fn utc_stamp(secs: u64) -> String {
        let (days, rem) = ((secs / 86_400) as i64, secs % 86_400);
        let z = days + 719_468; // 0000-03-01 起的天数
        let era = z.div_euclid(146_097);
        let doe = z.rem_euclid(146_097);
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let d = doy - (153 * mp + 2) / 5 + 1;
        let m = if mp < 10 { mp + 3 } else { mp - 9 };
        let y = yoe + era * 400 + i64::from(m <= 2);
        format!(
            "{y:04}-{m:02}-{d:02} {:02}:{:02} UTC",
            rem / 3600,
            (rem % 3600) / 60
        )
    }

    fn about_window(&mut self, ctx: &egui::Context) {
        if !self.show_about {
            return;
        }
        let mut open = true;
        egui::Window::new(crate::i18n::text(ctx,"About AppimageGene"))
            .collapsible(false)
            .resizable(false)
            .default_width(560.0)
            .open(&mut open)
            .show(ctx, |ui| {
                ui.label(egui::RichText::new("AppimageGene").strong().size(15.0));
                ui.label(crate::i18n::text(ui.ctx(),"A Rust re-implementation of the SnapGene plasmid viewer/editor."));
                ui.add_space(4.0);
                ui.label(crate::ui_format!(ui.ctx(),"Version {}", env!("CARGO_PKG_VERSION")));
                ui.label(crate::i18n::text(ui.ctx(),"Formats: .dna (round-trip verified)"));
                // 二级结构用的 ViennaRNA 是**编进二进制**的，不依赖系统安装。
                ui.label(crate::ui_format!(ui.ctx(),"ViennaRNA {} (bundled)", hx_vienna::VERSION));
                ui.label(
                    egui::RichText::new(crate::i18n::text(ui.ctx(),"RNA folding & layout · no Python required"))
                        .size(11.0)
                        .color(theme::TEXT_DIM),
                );
                // 二进制自身的构建时间/路径：跑 release 时最容易踩的坑是**忘了重新 build**，
                // 于是「刚做的功能没出现」。把出处写在脸上，一眼看出是不是旧二进制。
                if let Ok(exe) = std::env::current_exe() {
                    let built = std::fs::metadata(&exe).and_then(|m| m.modified()).ok();
                    let (stamp, age) = match built {
                        Some(t) => {
                            let epoch = t
                                .duration_since(std::time::UNIX_EPOCH)
                                .map(|d| d.as_secs())
                                .unwrap_or(0);
                            let age = std::time::SystemTime::now()
                                .duration_since(t)
                                .map(|d| {
                                    let s = d.as_secs();
                                    match s {
                                        0..=90 => format!("{s} s ago"),
                                        91..=5400 => format!("{} min ago", s / 60),
                                        5401..=172_800 => format!("{} h ago", s / 3600),
                                        _ => format!("{} d ago", s / 86_400),
                                    }
                                })
                                .unwrap_or_else(|_| "in the future".to_string());
                            (Self::utc_stamp(epoch), age)
                        }
                        None => ("unknown".to_string(), "unknown".to_string()),
                    };
                    ui.add_space(6.0);
                    // 只显示路径尾部两级：整条路径太长，About 窗口里会被裁掉，
                    // 而「是不是 target/release 下那个二进制」看尾部就够。
                    let tail = {
                        let mut it = exe.components().rev();
                        match (it.next(), it.next()) {
                            (Some(file), Some(dir)) => {
                                format!("…/{}/{}", dir.as_os_str().to_string_lossy(), file.as_os_str().to_string_lossy())
                            }
                            _ => exe.display().to_string(),
                        }
                    };
                    ui.label(
                    egui::RichText::new(crate::ui_format!(ui.ctx(),"Binary: {}",tail))
                            .size(11.0)
                            .color(theme::TEXT_DIM),
                    );
                    ui.label(
                    egui::RichText::new(crate::ui_format!(ui.ctx(),"Built: {} ({})",stamp,age))
                            .size(11.0)
                            .color(theme::TEXT_DIM),
                    );
                }
            });
        self.show_about = open;
    }

    fn shortcuts_window(&mut self, ctx: &egui::Context) {
        if !self.show_shortcuts {
            return;
        }
        let mut open = true;
        egui::Window::new(crate::i18n::text(ctx,"Keyboard Shortcuts"))
            .collapsible(false)
            .resizable(false)
            .open(&mut open)
            .show(ctx, |ui| {
                egui::Grid::new("shortcuts").num_columns(2).show(ui, |ui| {
                    for (k, v) in [
                        ("Ctrl+O", "Open (focus path box)"),
                        ("Ctrl+S", "Save"),
                        ("Ctrl+Shift+S", "Save As…"),
                        ("Ctrl+E", "Export current view as SVG"),
                        ("Ctrl+Z / Ctrl+Shift+Z", "Undo / Redo"),
                        ("Ctrl+X / Ctrl+C / Ctrl+V", "Cut / Copy / Paste at cursor"),
                        ("Del", "Delete selection"),
                        ("Ctrl+A", "Select whole sequence"),
                        ("Esc", "Clear selection"),
                        ("Ctrl+0", "Fit to window"),
                        ("Ctrl+ + / Ctrl+ -", "Zoom in / out"),
                        ("Drag", "Sequence: select · Linear map: pan"),
                        ("Ctrl+wheel", "Zoom"),
                    ] {
                        ui.label(egui::RichText::new(k).strong());
                        ui.label(crate::i18n::text(ui.ctx(),v));
                        ui.end_row();
                    }
                });
            });
        self.show_shortcuts = open;
    }

    fn info_panel(&mut self, ui:&mut egui::Ui) {
        if !self.show_info {return;}
        if self.viewer.is_some() {
            // 阅读器没有分子元数据可编辑。
            self.show_info = false;
            return;
        }
        let Some(doc)=self.doc.as_mut() else{self.show_info=false;return;};
        egui::Panel::right("native_info_panel").default_size(490.0).min_size(340.0).max_size(650.0)
            .frame(egui::Frame::new().fill(egui::Color32::from_gray(236)))
            .show(ui,|ui|{if doc.native_class().is_some(){ui.label(crate::i18n::text(ui.ctx(),"Metadata editing is unavailable for this native-only document."));}
                else if let Err(e)=crate::info_panel::draw(ui,doc){self.status=format!("Info edit failed: {e}");}});
    }

    fn dump_blocks(&mut self) {
        let Some(doc) = self.doc.as_ref() else {
            return;
        };
        let mut parts: Vec<String> = Vec::new();
        for b in &doc.session.file().blocks {
            parts.push(format!("{}·{}B", b.block_type, b.payload.len()));
        }
        self.status = format!("{} block(s): {}", parts.len(), parts.join(" "));
    }

    /// 剪贴板写入请求（`ui()` 末统一提交给 `Context`）。
    fn take_pending_copy(&mut self) -> Option<String> {
        self.pending_copy.take()
    }
    /// 文档画布：主窗口与克隆对话框共用的一套 Map / Sequence / Enzymes / Primers 视图。
    ///
    /// 克隆对话框的预览文档放在一个独立的 `HelixApp` 里（见
    /// [`HelixApp::preview`]），因此这里的 view 状态互不污染；`mode` 只影响
    /// 页签集合和“无文档”提示。
    pub fn doc_canvas(&mut self, ui: &mut egui::Ui, mode: CanvasMode) {
        self.last_size = Some(ui.available_size());
        self.cloning_canvas = mode == CanvasMode::Cloning;
        let central_rect = ui.max_rect();
        self.central_rect = Some(central_rect);
        if self.doc.is_some() && self.sites_dirty {
            self.recompute_sites();
        }
        if self.doc.is_none() {
            if mode == CanvasMode::Cloning {
                ui.centered_and_justified(|ui| {
                    ui.label(crate::i18n::text(ui.ctx(),"Please choose a source file."));
                });
            }
            return;
        }
        let tabs: &[Tab] = match mode {
            CanvasMode::Main => Tab::list_for(self.doc.as_ref()),
            CanvasMode::Cloning => &[Tab::Map, Tab::Sequence, Tab::Enzymes, Tab::Primers],
        };
        if !tabs.contains(&self.tab) {
            self.tab = tabs.first().copied().unwrap_or(Tab::Map);
        }
        let canvas = match self.tab {
            Tab::Map => Some(crate::view_map::draw(ui, self)),
            Tab::Sequence => {
                // RNA / 蛋白 / 单链 DNA 走「单链序列视图」（SnapGene 也是这样分开的）
                if self.doc.as_ref().is_some_and(|d| d.is_single_stranded()) {
                    Some(crate::view_flat::draw(ui, self))
                } else {
                    Some(crate::view_seq::draw(ui, self))
                }
            }
            Tab::SecondaryStructure => Some(crate::view_struct::draw(ui, self)),
            Tab::Properties => {
                crate::panels::draw_properties(ui, self);
                None
            }
            Tab::NativeOverview => {
                crate::panels::draw_native_overview(ui, self);
                None
            }
            Tab::Enzymes | Tab::Features | Tab::Primers | Tab::History => {
                crate::panels::draw(ui, self);
                None
            }
        };
        let ctx = ui.ctx().clone();
        // New File 面板盖在文档区左上角（SnapGene 实测）：它打开时视图页签被挡住；
        // 克隆对话框没有这些浮动面板，始终显示四个页签。
        if mode == CanvasMode::Cloning
            || (self.new_file.is_none() && self.feature_dlg.is_none())
        {
            self.view_tabs_floating(&ctx, central_rect, tabs, ui.id());
        }
        if let Some(r) = canvas {
            // 对话框打开时左侧工具条/视图页签让位（SnapGene 的对话框盖住它们）
            if self.show_side_tools
                && self.tab != Tab::SecondaryStructure
                && self.feature_dlg.is_none()
                && self.new_file.is_none()
            {
                self.tool_palette(&ctx, r);
            }
            self.canvas_controls(&ctx, r);
            if mode==CanvasMode::Main {self.canvas_search_button(&ctx,r);}
        }
    }

}

/// 6+ 切点、非 IIS 酶的库下标。
pub fn cutters_six_plus() -> Vec<usize> {
    hx_enzymes::ENZYMES
        .iter()
        .enumerate()
        .filter(|(_, e)| e.site_len() >= 6 && !e.is_type_iis())
        .map(|(i, _)| i)
        .collect()
}

/// 全库中识别位点**与切点**去重（同切点只留名字最短者，等长时偏好常用酶——SnapGene 在
/// `ACTAGT` 上显示 SpeI 而不是 AhlI）——SnapGene "Nonredundant" 的近似。
///
/// 键是 `位点#切点偏移`：`GGGCCC` 上的 ApaI（cut 5）与 PspOMI（cut 1）切在不同位置，
/// SnapGene 两个都列（4675 / 4671），因此不能只按位点去重。
pub fn nonredundant_cutters() -> Vec<usize> {
    use std::collections::{HashMap, HashSet};
    let common: HashSet<usize> = hx_enzymes::common_indices().into_iter().collect();
    let preferred: HashSet<usize> = hx_enzymes::preferred_indices().into_iter().collect();
    // 排序键：SnapGene 偏好名 > 短名 > 常用集内 > 表序（表序即 `best` 的插入顺序）
    let rank = |j: usize| {
        (
            !preferred.contains(&j),
            hx_enzymes::ENZYMES[j].name.len(),
            !common.contains(&j),
        )
    };
    let mut best: HashMap<(String, i16), usize> = HashMap::new();
    for i in cutters_six_plus() {
        let e = &hx_enzymes::ENZYMES[i];
        let key = (e.site.to_ascii_uppercase(), e.cut);
        match best.get(&key).copied() {
            Some(j) if rank(j) <= rank(i) => {}
            _ => {
                best.insert(key, i);
            }
        }
    }
    let mut v: Vec<usize> = best.into_values().collect();
    v.sort_unstable();
    v
}

/// 在序列中**只切一刀**的 6+ 切点酶（`redundant` = 再去掉同识别位点的重复酶）。
pub fn unique_cutters(seq: &[u8], circular: bool, redundant: bool) -> Vec<usize> {
    let candidates = if redundant {
        nonredundant_cutters()
    } else {
        cutters_six_plus()
    };
    let mut out = Vec::new();
    for i in candidates {
        let e = &hx_enzymes::ENZYMES[i];
        if hx_enzymes::find_site_positions(seq, e.site.as_bytes(), circular).len() == 1 {
            out.push(i);
        }
    }
    out.sort_unstable();
    out
}

impl eframe::App for HelixApp {
    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        // eframe::App::save is not guaranteed without its optional persistence
        // feature. Our own settings must also be written on normal shutdown.
        self.save_settings();
    }
    /// eframe 每 30 秒及退出时调用：把显示/工具栏/侧栏状态写回设置文件，
    /// 与 SnapGene 重启后恢复显示开关的行为一致。
    fn save(&mut self, _storage: &mut dyn eframe::Storage) {
        self.save_settings();
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        // 文件拖放
        let dropped: Vec<PathBuf> = ui.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .map(|f| f.path().to_path_buf())
                .collect()
        });
        if let Some(p) = dropped.into_iter().next() {
            self.open(&p);
        }

        // 快捷键（文本框获得焦点时不拦截）。
        // 用 `consume_key` 按**事件自身**的修饰键匹配：先具体后一般
        // （Ctrl+Shift+Z 先于 Ctrl+Z —— matches_logically 允许模式忽略多余的 Shift）。
        let dialog_open = self.feature_dlg.is_some()
            || self.rename_file.is_some()
            || self.enzyme_preferences_open
            || self.hybridization_draft.is_some() || self.show_translation_options
            || self.molecular_action.is_some()
            || self.show_paste_dialog || self.pending_delete.is_some()
            || self.new_file.is_some()
            || self.primer_dlg.is_some()
            || self.primer_choice.is_some()
            || self.design_dlg.is_some()
            || self.cloning.is_some();
        if !dialog_open&&!ui.ctx().egui_wants_keyboard_input()&&ui.input_mut(|i|i.consume_key(egui::Modifiers::NONE,egui::Key::F2)){
            let path=self.doc.as_ref().and_then(|d|d.path.clone()).or_else(||self.viewer.as_ref().map(|v|v.path().to_path_buf()));
            if let Some(path)=path{self.rename_file=Some(crate::rename::Rename::new(path));}
        }
        let (save_as, save, export, open, fit, new_file, new_alt) = ui.input_mut(|i| {
            let cmd_shift = egui::Modifiers::COMMAND | egui::Modifiers::SHIFT;
            // Ctrl+N 与 Ctrl+Alt+N 共用一个匹配（`matches_logically` 会忽略多余的修饰键），
            // 用**事件自身的 alt**区分：Ctrl+N = DNA/RNA，Ctrl+Alt+N = 蛋白。
            let alt = i.modifiers.alt;
            (
                i.consume_key(cmd_shift, egui::Key::S),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::S),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::E),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::O),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::Num0),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::N),
                alt,
            )
        });
        // Ctrl+T = Add Feature…（文本框获得焦点或对话框打开时不拦截）
        let add_feature = if dialog_open {
            false
        } else {
            ui.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::T))
        };
        if add_feature {
            self.open_add_feature_dialog(false);
        }
        // Ctrl+Alt+I = Restriction and Insertion Cloning ▸ Insert Fragment…
        let insert_fragment = if dialog_open {
            false
        } else {
            ui.input_mut(|i| {
                i.consume_key(
                    egui::Modifiers::COMMAND | egui::Modifiers::ALT,
                    egui::Key::I,
                )
            })
        };
        if insert_fragment {
            self.open_cloning(hx_algo::cloning::Method::Restriction);
        }
        // Ctrl+R = Add Primer…、Ctrl+Shift+R = Duplicate Primer（SnapGene 的快捷键）。
        // `consume_key` 的逻辑匹配允许模式忽略事件里多余的修饰键，所以**先**匹配
        // 更具体的 Ctrl+Shift+R，否则它会被 Ctrl+R 那条吃掉。
        let (dup_primer, add_primer) = if dialog_open {
            (false, false)
        } else {
            ui.input_mut(|i| {
                (
                    i.consume_key(
                        egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                        egui::Key::R,
                    ),
                    i.consume_key(egui::Modifiers::COMMAND, egui::Key::R),
                )
            })
        };
        if dup_primer {
            self.duplicate_selected_primer();
        } else if add_primer {
            self.open_add_primer_dialog();
        }
        if new_file {
            let kind = if new_alt {
                hx_core::MoleculeKind::Protein
            } else {
                hx_core::MoleculeKind::Dna
            };
            self.open_new_file_dialog(kind);
        }
        if save_as {
            self.do_save(true);
        }
        if save {
            self.do_save(false);
        }
        if export {
            self.do_export_svg();
        }
        if open {
            self.open_file_with_dialog();
        }
        if fit {
            // 阅读页签的 Ctrl+0 是“适合宽度 / 默认字号”，不是分子视图的 fit。
            if let Some(v) = self.viewer.as_mut() {
                v.fit();
            } else {
                self.fit();
            }
        }
        // 阅读器缩放：Ctrl+`+` / Ctrl+`-`。必须放在
        // `egui_wants_keyboard_input` 判断**之前**，否则文本框获得焦点时
        // （正是编辑 .txt 的时候）快捷键会被吞掉。
        if self.viewer.is_some() {
            let (viewer_in, viewer_out) = ui.input_mut(|i| {
                let cmd = egui::Modifiers::COMMAND;
                (
                    i.consume_key(cmd, egui::Key::Plus) || i.consume_key(cmd, egui::Key::Equals),
                    i.consume_key(cmd, egui::Key::Minus),
                )
            });
            if let Some(v) = self.viewer.as_mut() {
                if viewer_in {
                    v.zoom_in();
                }
                if viewer_out {
                    v.zoom_out();
                }
            }
        }

        // 文档标签页：Ctrl+Tab / Ctrl+Shift+Tab（SnapGene 的 Next/Previous tab）、
        // Ctrl+W（关闭标签页）。`consume_key` 的匹配**允许忽略多余 Shift**，
        // 所以先消费更具体的 Ctrl+Shift+Tab，再消费 Ctrl+Tab。
        let (prev_tab, next_tab, close_tab) = ui.input_mut(|i| {
            (
                i.consume_key(
                    egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                    egui::Key::Tab,
                ),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::Tab),
                i.consume_key(egui::Modifiers::COMMAND, egui::Key::W),
            )
        });
        if next_tab {
            self.next_tab(1);
        }
        if prev_tab {
            self.next_tab(-1);
        }
        if close_tab {
            self.close_doc();
        }

        // 其余快捷键在文本框获得焦点时让位（字段内 Ctrl+Z/Delete/Esc 归 egui 处理）
        if !ui.ctx().egui_wants_keyboard_input() && !dialog_open {
            // SnapGene 的 Copy 下拉快捷键：Ctrl+Shift+C = 下链 5′→3′、Ctrl+Alt+C = 1 字母氨基酸。
            // `consume_key` 会忽略多余的 Shift，所以**先**消费这两条，再消费 Ctrl+C；
            // 另外 egui-winit 把 Ctrl+Shift+C 也翻译成 `Event::Copy`，下面用标志拦住它。
            let (copy_bottom, copy_aa1) = ui.input_mut(|i| {
                let cmd_shift = egui::Modifiers::COMMAND | egui::Modifiers::SHIFT;
                let cmd_alt = egui::Modifiers::COMMAND | egui::Modifiers::ALT;
                let bottom = i.consume_key(cmd_shift, egui::Key::C);
                let aa1 = i.consume_key(cmd_alt, egui::Key::C);
                (bottom, aa1)
            });

            // 剪贴板：Ctrl+C / Ctrl+X / Ctrl+V（Paste 由 egui 的 Paste 事件送达）
            let (copy, cut, paste_text) = ui.input_mut(|i| {
                let cmd = egui::Modifiers::COMMAND;
                let copy = !copy_bottom
                    && !copy_aa1
                    && (i.consume_key(cmd, egui::Key::C)
                        || i.events.iter().any(|e| matches!(e, egui::Event::Copy)));
                let cut = i.consume_key(cmd, egui::Key::X)
                    || i.events.iter().any(|e| matches!(e, egui::Event::Cut));
                let paste = i.events.iter().find_map(|e| match e {
                    egui::Event::Paste(s) => Some(s.clone()),
                    _ => None,
                });
                (copy, cut, paste)
            });
            if copy {
                self.do_copy();
            }
            if cut {
                self.do_cut();
            }
            let reverse_shortcut=ui.input_mut(|i|i.consume_key(egui::Modifiers::COMMAND|egui::Modifiers::SHIFT,egui::Key::V));
            let _=ui.input_mut(|i|i.consume_key(egui::Modifiers::COMMAND,egui::Key::V));
            if let Some(text) = paste_text {
                let reverse=self.pending_paste_reverse.take().unwrap_or(reverse_shortcut);
                if !text.is_empty() {
                    self.paste_sequence(&text,reverse);
                }
            }else if reverse_shortcut {self.request_paste(ui.ctx(),true);
            }
            if copy_bottom {
                self.do_copy_form(CopyForm::BottomStrand5to3);
            }
            if copy_aa1 {
                self.do_copy_form(CopyForm::AminoAcids1);
            }

            let (redo, undo, del, backspace, esc, sel_all, zoom_in, zoom_out) =
                ui.input_mut(|i| {
                    let cmd = egui::Modifiers::COMMAND;
                    let cmd_shift = egui::Modifiers::COMMAND | egui::Modifiers::SHIFT;
                    let redo =
                        i.consume_key(cmd_shift, egui::Key::Z) || i.consume_key(cmd, egui::Key::Y);
                    let undo = i.consume_key(cmd, egui::Key::Z);
                    let del = i.consume_key(egui::Modifiers::NONE, egui::Key::Delete);
                    let backspace = i.consume_key(egui::Modifiers::NONE, egui::Key::Backspace);
                    let esc = i.consume_key(egui::Modifiers::NONE, egui::Key::Escape);
                    let sel_all = i.consume_key(cmd, egui::Key::A);
                    let zoom_in = i.consume_key(cmd, egui::Key::Plus)
                        || i.consume_key(cmd, egui::Key::Equals);
                    let zoom_out = i.consume_key(cmd, egui::Key::Minus);
                    (redo, undo, del, backspace, esc, sel_all, zoom_in, zoom_out)
                });
            if redo {
                self.do_redo();
            }
            if undo {
                self.do_undo();
            }
            if del {
                self.do_delete();
            }
            if backspace {
                if self.selection.is_some_and(|(s, e)| s < e) { self.do_delete(); }
                else if let Some(c) = self.cursor_bp.filter(|c| *c > 0) { self.pending_delete = Some((c - 1, c)); }
            }
            if esc {
                self.selection = None;
                self.selected_residue = None;
                self.status = "Selection cleared".into();
            }
            if sel_all {
                self.select_all();
            }
            if zoom_in {
                self.zoom_in();
            }
            if zoom_out {
                self.zoom_out();
            }

            // ---- 序列视图里的光标与直接输入 ----
            // SnapGene：单击放置光标后输入 = 插入；先选中一段再输入 = 替换。
            // 只有序列页签、且光标/选区存在时才接管键盘，其它视图的按键不受影响。
            if self.tab == Tab::Sequence
                && self.doc.as_ref().is_some_and(|d| d.native_class().is_none())
            {
                let (left, right, left_sel, right_sel, home, end, home_sel, end_sel) =
                    ui.input_mut(|i| {
                        let sh = egui::Modifiers::SHIFT;
                        let none = egui::Modifiers::NONE;
                        // Shift 先消费：`consume_key` 会忽略多余的 Shift
                        let ls = i.consume_key(sh, egui::Key::ArrowLeft);
                        let rs = i.consume_key(sh, egui::Key::ArrowRight);
                        let hs = i.consume_key(sh, egui::Key::Home);
                        let es = i.consume_key(sh, egui::Key::End);
                        let l = i.consume_key(none, egui::Key::ArrowLeft);
                        let r = i.consume_key(none, egui::Key::ArrowRight);
                        let h = i.consume_key(none, egui::Key::Home);
                        let e = i.consume_key(none, egui::Key::End);
                        (l, r, ls, rs, h, e, hs, es)
                    });
                if left_sel {
                    self.move_caret(-1, true);
                } else if right_sel {
                    self.move_caret(1, true);
                } else if left {
                    self.move_caret(-1, false);
                } else if right {
                    self.move_caret(1, false);
                }
                if home_sel {
                    self.move_caret_home(false, true);
                } else if end_sel {
                    self.move_caret_home(true, true);
                } else if home {
                    self.move_caret_home(false, false);
                } else if end {
                    self.move_caret_home(true, false);
                }

                if self.selection.is_some() || self.cursor_bp.is_some() {
                    let typed = ui.input_mut(|i| {
                        let mut s = String::new();
                        i.events.retain(|e| match e {
                            egui::Event::Text(t) => {
                                s.push_str(t);
                                false
                            }
                            _ => true,
                        });
                        s
                    });
                    if !typed.is_empty() {
                        self.paste = typed;
                        self.insert_unknown = false;
                        self.insert_reverse = false;
                        self.show_paste_dialog = true;
                    }
                }
            }
        }

        if self.needs_fit && self.last_size.is_some() {
            self.fit();
            self.needs_fit = false;
        }
        if self.sites_dirty {
            self.recompute_sites();
        }

        // 标题栏（SnapGene 实测：`文件名 (Circular / N bp)`、`文件.prot (N aa)`；
        // `*` 紧跟在名字后面 = 有未保存编辑；文件名带扩展名）
        let title = if let Some(v) = &self.viewer {
            let star = if v.is_dirty() { "*" } else { "" };
            format!("{}{star}  ({}) — AppimageGene", v.title(), crate::i18n::viewer_summary(ui.ctx(),&v.summary()))
        } else {
            match &self.doc {
            Some(d) => {
                let name = d.file_label();
                let star = if d.is_dirty() { "*" } else { "" };
                let body = if let Some(class) = d.native_class() {
                    // 比对/凝胶：没有可显示的序列，正文写类别
                    crate::i18n::text(ui.ctx(),class).to_string()
                } else if d.kind() == hx_core::MoleculeKind::Protein {
                    format!("{} aa", nfmt(d.len()))
                } else {
                    let topo = crate::i18n::text(ui.ctx(),if d.circular() { "Circular" } else { "Linear" });
                    format!("{topo} / {} {}", nfmt(d.len()), d.unit())
                };
                format!("{name}{star}  ({body}) — AppimageGene")
            }
            None => "AppimageGene".to_string(),
            }
        };
        if self.last_title.as_deref() != Some(title.as_str()) {
            ui.ctx()
                .send_viewport_cmd(egui::ViewportCommand::Title(title.clone()));
            self.last_title = Some(title);
        }

        // ---- 区域划分（对齐 SnapGene）----
        // Full-width menu → navigation/tabs → project panel alongside toolbar/canvas.
        egui::Panel::top("menubar")
            .exact_size(crate::chrome::MENU_H)
            .frame(egui::Frame::new().fill(theme::PANEL))
            .show(ui, |_ui| {});
        egui::Panel::top("doctabs")
            .exact_size(crate::chrome::DOC_TABS_H)
            .frame(egui::Frame::new().fill(egui::Color32::WHITE))
            .show(ui,|ui|self.navigation_row(ui));
        if self.sidebar_open {
            let sidebar=egui::Panel::left("files")
                .resizable(true)
                .default_size(self.sidebar_width)
                .min_size(200.0)
                .max_size(460.0)
                // SG folder view 侧栏是**纯白**底 + 右缘 1 px 亮蓝线（`sg-sb.png`）
                .frame(egui::Frame::new().fill(crate::sidebar::BG))
                .show(ui, |ui| {
                    crate::sidebar::draw(ui, self);
                    let r = ui.max_rect();
                    ui.painter().line_segment(
                        [
                            egui::pos2(r.right() - 0.5, r.top()),
                            egui::pos2(r.right() - 0.5, r.bottom()),
                        ],
                        egui::Stroke::new(1.0, crate::sidebar::EDGE),
                    );
                });
            self.sidebar_width=sidebar.response.rect.width();
        }
        if self.show_toolbar && self.viewer.is_none() {
            egui::Panel::top("toolbar")
                .exact_size(crate::chrome::TOOLBAR_H)
                .frame(egui::Frame::new().fill(egui::Color32::WHITE).inner_margin(egui::Margin::symmetric(2,3)))
                .show(ui, |ui| self.toolbar(ui));
        }
        egui::Panel::bottom("status")
            .exact_size(24.0)
            .frame(egui::Frame::new().fill(egui::Color32::WHITE))
            .show(ui, |ui| {
                let quiet=self.selection.is_none() && self.selected_feature.is_none() && self.cursor_bp.is_none() && self.hover_bp.is_none()
                    && (self.status.is_empty() || self.status.starts_with("Opened ") || self.status.starts_with("Created ") || self.status=="Selection cleared");
                if !quiet {self.status_bar(ui);}
            });
        self.info_panel(ui);
        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(match &self.viewer {
                // 文本编辑器是白纸；PDF 用中灰底衬托白页。
                Some(crate::viewer::ViewerDoc::Text(_)) => egui::Color32::WHITE,
                Some(crate::viewer::ViewerDoc::Pdf(_)) => egui::Color32::from_gray(232),
                None => match self.tab {
                    // Supplied SG reference: Map and Sequence have a light-gray canvas.
                    Tab::Sequence|Tab::Map => self.preferences.canvas_color(),
                    _ => theme::CANVAS,
                },
            }))
            .show(ui, |ui| {
                if let Some(v) = self.viewer.as_mut() {
                    self.central_rect = Some(ui.max_rect());
                    self.last_size = Some(ui.available_size());
                    v.draw(ui);
                } else {
                    self.doc_canvas(ui, CanvasMode::Main);
                }
            });

        if let Some(lines)=ui.ctx().data(|d|d.get_temp::<Vec<[egui::Pos2;2]>>(egui::Id::unique("document_row_outline"))) {
            let painter=ui.painter().with_clip_rect(ui.ctx().content_rect());
            if let Some(Some(gap))=ui.ctx().data(|d|d.get_temp::<Option<egui::Rect>>(egui::Id::unique("document_row_active_gap"))) {
                painter.rect_filled(gap,0.0,egui::Color32::WHITE);
            }
            if let Some((rect,color))=ui.ctx().data(|d|d.get_temp::<(egui::Rect,egui::Color32)>(egui::Id::unique("document_row_nav_boundary"))) {
                painter.rect_filled(rect,0.0,color);
            }
            for line in lines {painter.line_segment(line,egui::Stroke::new(1.0,crate::sidebar::EDGE));}
        }

        // 菜单栏覆盖层：必须在所有 Panel 之后绘制，egui 0.37 的 Popup 才会用全局坐标定位。
        egui::Area::new(egui::Id::unique("menubar_overlay"))
            .order(egui::Order::Foreground)
            .fixed_pos(egui::pos2(0.0, 0.0))
            .show(ui.ctx(), |ui| {
                egui::Frame::new()
                    .fill(theme::PANEL)
                    .inner_margin(egui::Margin::symmetric(6, 0))
                    .show(ui, |ui| {
                        ui.set_width(ui.ctx().content_rect().width() - 12.0);
                        ui.set_height(crate::chrome::MENU_H);
                        self.menu_bar(ui);
                    });
            });

        // 浮动窗口
        self.view_options_window(ui.ctx());
        self.new_file_window(ui.ctx());
        self.feature_dlg_window(ui.ctx());
        self.primer_choice_window(ui.ctx());
        self.primer_dlg_window(ui.ctx());
        self.design_dlg_window(ui.ctx());
        self.cloning_dlg_window(ui.ctx());
        if self.gg_settings_open {
            if let Some(result) = crate::ggsettings::draw(&mut self.gg_settings, ui.ctx()) {
                self.gg_settings_open = false;
                match result {
                    crate::ggsettings::GgClose::Ok => {
                        // 与 SG 一样，OK 后默认值/勾选清单跨启动保留。
                        self.cloning_settings.golden_gate_default_enzyme =
                            self.gg_settings.default_enzyme_name().to_string();
                        self.save_settings();
                    }
                    crate::ggsettings::GgClose::Cancel => {
                        if let Some(old) = self.gg_settings_backup.take() {
                            self.gg_settings = old;
                        }
                    }
                }
                self.gg_settings_backup = None;
            }
        }
        self.primer_analysis_window(ui.ctx());
        self.feature_color_window(ui.ctx());
        self.paste_dialog_window(ui.ctx());
        self.display_options_dialogs(ui.ctx());
        if !dialog_open && !ui.ctx().egui_wants_keyboard_input() && ui.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::K)) {
            self.open_molecular_action(crate::molecular_actions::Action::Protein);
        }
        if let Some(result) = self.molecular_action.as_mut().and_then(|d| d.show(ui.ctx())) {
            self.molecular_action = None;
            if let Ok((bytes, name)) = result { self.adopt_cloned_product(bytes, name); }
        }
        self.open_dialog_window(ui.ctx());
        self.save_dialog_window(ui.ctx());
        self.about_window(ui.ctx());
        self.shortcuts_window(ui.ctx());

        // 剪贴板写入（需要 `Context`）
        if let Some(text) = self.take_pending_copy() {
            ui.ctx().copy_text(text);
        }
    }
}

/// 画布工具条上的白底文字按钮（`Export as… ▾` / `Display ▾`）。
fn struct_button(ui: &mut egui::Ui, label: &str, color: egui::Color32) -> egui::Response {
    let label=crate::i18n::text(ui.ctx(),label);
    let font = egui::FontId::proportional(13.0);
    let w = ui
        .painter()
        .layout_no_wrap(label.to_string(), font.clone(), color)
        .size()
        .x;
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(w + 30.0, 26.0), egui::Sense::click());
    let bg = if resp.hovered() {
        theme::ACCENT_BG
    } else {
        egui::Color32::WHITE
    };
    let p = ui.painter();
    p.rect(
        rect,
        4.0,
        bg,
        egui::Stroke::new(0.7, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    p.text(
        egui::pos2(rect.left() + 9.0 + w * 0.5, rect.center().y),
        egui::Align2::CENTER_CENTER,
        label,
        font,
        color,
    );
    // 下拉小三角（SnapGene 的 ▾；系统 Noto Sans 缺 U+25BE，直接画）
    let cx = rect.left() + 9.0 + w + 7.0;
    let cy = rect.center().y + 1.0;
    p.add(egui::Shape::convex_polygon(
        vec![
            egui::pos2(cx - 3.5, cy - 2.0),
            egui::pos2(cx + 3.5, cy - 2.0),
            egui::pos2(cx, cy + 2.5),
        ],
        color,
        egui::Stroke::NONE,
    ));
    resp
}

/// 视图页签按钮（统一白底，活动页签用底部蓝线标识）。
pub fn tab_button(ui: &mut egui::Ui, active: bool, label: &str) -> egui::Response {
    let font = egui::FontId::proportional(crate::chrome::VIEW_TAB_FONT);
    let text_w = ui
        .painter()
        .layout_no_wrap(label.to_string(), font.clone(), theme::TEXT)
        .size()
        .x;
    // oracle 实测：SnapGene 的视图页签条高 30 px、页签内边距比我们原来宽
    let size = egui::vec2(text_w + 18.0, crate::chrome::VIEW_TAB_H);
    let (rect, resp) = ui.allocate_exact_size(size, egui::Sense::click());
    let bg = if active {
        egui::Color32::WHITE
    } else if resp.hovered() {
        theme::ACCENT_BG
    } else {
        egui::Color32::WHITE
    };
    let p = ui.painter();
    p.rect_filled(rect, egui::CornerRadius::ZERO, bg);
    if active {
        p.line_segment(
            [rect.left_bottom(), rect.right_bottom()],
            egui::Stroke::new(2.0, theme::ACCENT),
        );
    }
    p.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        label,
        font,
        theme::TEXT,
    );
    resp
}

/// 文档标签页按钮：`[ 名字 × ]`。
///
/// 关闭区用 [`egui::Ui::interact`] 单独做命中判定（而不是按坐标猜），
/// 返回 `(整块响应, 关闭按钮响应)`。宽度由调用方算好（见 [`HelixApp::doc_tabs`]）：
/// SnapGene 实测文档多时页签**等比压窄**并把名字用 `…` 截断成 `01-pD…`，
/// 只有活动页签画 ×。
pub fn tab_button_with_close(
    ui: &mut egui::Ui,
    active: bool,
    label: &str,
    salt: usize,
    kind: Option<hx_core::MoleculeKind>,
    width: f32,
) -> (egui::Response, egui::Response) {
    let font = crate::fonts::document_title_font(ui.ctx(),active);
    // 左内边距 5 + 图标 14 + 间隙 5 + 文字 + 右内边距 3 + （活动页签才有的）关闭按钮 16
    let close_w = if active { 16.0 } else { 3.0 };
    let size = egui::vec2(width, crate::chrome::DOC_TABS_H);
    let (rect, resp) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
    let close_rect = egui::Rect::from_min_max(
        egui::pos2(rect.right() - 16.0, rect.top() + 4.0),
        egui::pos2(rect.right() - 2.0, rect.bottom() - 4.0),
    );
    // 只有活动页签才有 ×（oracle 实测：非活动页签不画关闭按钮）
    let close = ui.interact(
        if active {
            close_rect
        } else {
            egui::Rect::from_min_size(rect.right_top(), egui::Vec2::ZERO)
        },
        ui.scope_id().with(("tab_close", salt)),
        egui::Sense::click(),
    );
    let bg = if active {
        theme::CANVAS
    } else if resp.hovered() || close.hovered() {
        theme::ACCENT_BG
    } else {
        theme::CANVAS
    };
    {
        let p = ui.painter();
        // Active document tab joins the white toolbar underneath, with NO blue
        // lower edge. Cover the row baseline at the bottom pixel as well.
        p.rect_filled(if active {rect.expand2(egui::vec2(0.0,1.0))}else{rect}, egui::CornerRadius::ZERO, bg);
        let edge=egui::Stroke::new(0.8,if active {crate::sidebar::EDGE}else{theme::BORDER});
        p.line_segment([rect.left_top(),rect.left_bottom()],edge);
        p.line_segment([rect.right_top(),rect.right_bottom()],edge);
        if !active {p.line_segment([rect.left_bottom(),rect.right_bottom()],egui::Stroke::new(0.8,crate::sidebar::EDGE));}
        // 分子类型图标
        let icon = egui::Rect::from_center_size(
            egui::pos2(rect.left() + 5.0 + 7.0, rect.center().y),
            egui::vec2(14.0, 14.0),
        );
        match kind {
            // Supplied SG reference: outline on inactive tabs, filled on active.
            Some(kind) => theme::kind_icon(p, icon, kind, active),
            None => crate::branding::file_icon(p,icon,label,false),
        }
        let text_rect = egui::Rect::from_min_max(
            egui::pos2(rect.left() + 5.0 + 14.0 + 5.0, rect.top()),
            egui::pos2(rect.right() - close_w - 2.0, rect.bottom()),
        );
        let text_col = if active { theme::TEXT } else { theme::TEXT_DIM };
        let galley = truncate_galley(p, label, font, text_col, text_rect.width());
        p.galley(
            egui::pos2(
                text_rect.left(),
                text_rect.center().y - galley.size().y * 0.5,
            ),
            galley,
            text_col,
        );
        // 关闭小叉（悬停时加深底色）
        if active {
            if close.hovered() {
                p.rect_filled(close_rect, 2.0, theme::with_alpha(theme::TEXT_DIM, 40));
            }
            let c = close_rect.center();
            let col = if close.hovered() {
                theme::TEXT
            } else {
                theme::TEXT_DIM
            };
            let s = 4.0;
            let stroke = egui::Stroke::new(1.2, col);
            p.line_segment([c + egui::vec2(-s, -s), c + egui::vec2(s, s)], stroke);
            p.line_segment([c + egui::vec2(-s, s), c + egui::vec2(s, -s)], stroke);
        }
    }
    (resp, close)
}

/// 单行文字；放不下就用 `…` 截断（SnapGene 的文档页签就是这样压名字的）。
fn truncate_galley(
    p: &egui::Painter,
    text: &str,
    font: egui::FontId,
    color: egui::Color32,
    max_w: f32,
) -> std::sync::Arc<egui::Galley> {
    let mut job = egui::text::LayoutJob::default();
    job.wrap.max_width = max_w.max(6.0);
    job.wrap.max_rows = 1;
    job.wrap.break_anywhere = true;
    job.wrap.overflow_character = Some('…');
    job.append(
        text,
        0.0,
        egui::TextFormat {
            font_id: font,
            color,
            ..Default::default()
        },
    );
    p.layout_job(job)
}

/// 千分位（统一实现见 [`hx_render::thousands`]）。
/// 菜单里下链两项的预览：SnapGene 取复制文本的**前 6 个字符**再加省略号
/// （oracle 实测：选 11 bp 显示 `5' → 3'  (CCCATT…)`）。
fn copy_preview(text: &str) -> String {
    let mut head: String = text.chars().take(6).collect();
    if text.chars().count() > 6 {
        head.push('…');
    }
    head
}

/// `Copy` 下拉菜单（工具栏分体按钮与 Edit 菜单共用同一份结构）。
///
/// 逐项照 SnapGene 8.2.3 的菜单原文与快捷键；`Copy Map` 只在图谱视图可用，
/// 序列页签里灰显（与截图一致）。
fn copy_menu(ui: &mut egui::Ui, app: &mut HelixApp, has_sel: bool) -> Option<CopyForm> {
    let mut chosen = None;
    if item(ui, "Copy Top Strand Bases", Some("Ctrl+C"), has_sel) {
        chosen = Some(CopyForm::TopStrand);
    }
    // 子菜单标题里的预览要现算（选段变了就变）
    let preview53 = app.copy_text_for(CopyForm::BottomStrand5to3);
    let preview35 = app.copy_text_for(CopyForm::BottomStrand3to5);
    let label53 = match &preview53 {
        Some(t) if !t.is_empty() => format!("5′ → 3′  ({})", copy_preview(t)),
        _ => "5′ → 3′".to_string(),
    };
    let label35 = match &preview35 {
        Some(t) if !t.is_empty() => format!("3′ → 5′  ({})", copy_preview(t)),
        _ => "3′ → 5′".to_string(),
    };
    ui.menu_button(crate::i18n::text(ui.ctx(),"Copy Bottom Strand Bases"), |ui| {
        if item(ui, &label53, Some("Ctrl+Shift+C"), has_sel) {
            chosen = Some(CopyForm::BottomStrand5to3);
            ui.close();
        }
        if item(ui, &label35, None, has_sel) {
            chosen = Some(CopyForm::BottomStrand3to5);
            ui.close();
        }
    });
    ui.menu_button(crate::i18n::text(ui.ctx(),"Copy Amino Acids"), |ui| {
        if item(ui, "Copy 1-Letter Amino Acids", Some("Ctrl+Alt+C"), has_sel) {
            chosen = Some(CopyForm::AminoAcids1);
            ui.close();
        }
        if item(ui, "Copy 3-Letter Amino Acids", None, has_sel) {
            chosen = Some(CopyForm::AminoAcids3);
            ui.close();
        }
    });
    // `Copy Map`：图谱视图才有意义（SG 在 Sequence 页签里同样灰显）
    let _ = item(ui, "Copy Map", None, false);
    chosen
}

/// IUPAC 核酸字母（含 `U`；大小写都收）。直接输入时用它过滤按键。
fn is_iupac_base(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B'
            | b'D' | b'H' | b'V' | b'N'
    )
}

/// 菜单项：`enabled=false` ⇒ 置灰（SnapGene 菜单结构里我们尚未实现的动作）。
fn item(ui: &mut egui::Ui, label: &str, shortcut: Option<&str>, enabled: bool) -> bool {
    let mut button=egui::Button::new(crate::i18n::text(ui.ctx(),label));
    if let Some(shortcut)=shortcut {button=button.shortcut_text(shortcut);}
    ui.add_enabled(enabled,button).clicked()
}

/// 勾选式菜单项（SnapGene 的显示开关项）。
fn item_checked(ui: &mut egui::Ui, label: &str, checked: bool) -> bool {
    let label=crate::i18n::text(ui.ctx(),label);
    let text = if checked {
        format!("✓ {label}")
    } else {
        format!("   {label}")
    };
    ui.selectable_label(checked, text).clicked()
}

pub fn nfmt(n: usize) -> String {
    hx_render::thousands(n)
}

/// 不带千位分隔的数字（SnapGene 的表格、刻度与坐标都不分组；分组只出现在图谱标签里）。
pub fn nfmt_plain(n: usize) -> String {
    n.to_string()
}

/// SnapGene 的 `Copy` 下拉：同一条选段的五种取法。
///
/// 菜单原文与快捷键（oracle `:99` 的 SnapGene Viewer 8.2.3 实测）：
/// `Copy Top Strand Bases  Ctrl+C` /
/// `Copy Bottom Strand Bases ▸ [5'→3' (CCCATT…)  Ctrl+Shift+C, 3'→5' (ACCGCT…)]` /
/// `Copy Amino Acids ▸ [Copy 1-Letter Amino Acids  Ctrl+Alt+C, Copy 3-Letter Amino Acids]` /
/// `Copy Map`。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CopyForm {
    TopStrand,
    BottomStrand5to3,
    BottomStrand3to5,
    AminoAcids1,
    AminoAcids3,
}

impl CopyForm {
    pub fn is_amino(self) -> bool {
        matches!(self, Self::AminoAcids1 | Self::AminoAcids3)
    }
}

/// 三字母用空格分隔（SnapGene 序列视图和 `Copy 3-Letter Amino Acids` 都是这个写法），
/// 单字母连写（实测：`TGGCGAATGGG` → `WRM`）。
fn join_codes(codes: &[&'static str], three: bool) -> String {
    if three {
        codes.join(" ")
    } else {
        codes.concat()
    }
}

/// 强调色（旧的 `sidebar`/`view_*` 引用名）。
pub fn accent() -> egui::Color32 {
    theme::ACCENT
}

/// 选中区间的 GC 含量（整数百分比，四舍五入）。
///
/// oracle `/tmp/opencode/sg/dn1-status-z.png`：SG 的 DNA 状态栏在 `Selected:` 后写
/// `(867 .. 1455 = 589 bp)   [54% GC]` —— 长度用 `=` 分隔、后面方括号里是选中区的
/// GC 百分比（pUC19 的 ori 段实测 54%）。区间按 0-based 半开、可跨原点。
pub fn gc_percent(seq: &[u8], a: usize, b: usize, len: usize) -> usize {
    if len == 0 {
        return 0;
    }
    let mut gc = 0usize;
    let mut total = 0usize;
    for i in a..b {
        let Some(&base) = seq.get(i % len) else {
            break;
        };
        total += 1;
        if matches!(base.to_ascii_uppercase(), b'G' | b'C') {
            gc += 1;
        }
    }
    if total == 0 {
        return 0;
    }
    // 四舍五入（SG 的 54% 来自 589 bp 里 320 个 G/C = 54.3%）。
    (gc * 200 + total) / (total * 2)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_and_preview_apps_cannot_persist_over_real_workspace(){
        let app=super::HelixApp::blank();assert!(!app.persistence_enabled);
        let bytes=hx_formats::write::NewFile::dna(b"ACGT".to_vec()).to_bytes();
        let doc=crate::doc::Document::from_bytes(bytes).unwrap();
        let preview=super::HelixApp::preview(&eframe::egui::Context::default(),doc);
        assert!(!preview.persistence_enabled);preview.save_settings();app.save_settings();
    }
    #[test]
    fn restored_folder_is_expanded_and_loaded_without_open_tabs(){
        let mut app=super::HelixApp::blank();
        let folder=std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        app.folder=Some(folder.clone());app.restore_folder_tree();
        assert!(app.expanded.contains(&folder));assert!(app.dir_cache.contains_key(&folder));
        assert_eq!(app.folder_input,folder.display().to_string());
    }
    #[test]
    fn workspace_active_index_excludes_unsaved_tabs(){
        let mut app=super::HelixApp::blank();app.restoring_workspace=true;
        let bytes=hx_formats::write::NewFile::dna(b"ACGT".repeat(3)).to_bytes();
        for path in [Some("first.dna"),None,Some("last.dna")]{let mut doc=crate::doc::Document::from_bytes(bytes.clone()).unwrap();doc.path=path.map(std::path::PathBuf::from);app.push_doc(doc,path.map(std::path::PathBuf::from));}
        let(files,active)=app.workspace_snapshot();assert_eq!(files.len(),2);assert_eq!(active,1);
    }
    #[test]
    fn imported_preferences_apply_without_changing_existing_native_files() {
        let bytes=hx_formats::write::NewFile::dna(b"ACGT".repeat(50)).to_bytes();
        let mut app=super::HelixApp::blank();
        app.preferences.default_author="Test author".into();
        app.preferences.molecule_defaults[0].line_width=100;
        app.preferences.molecule_defaults[0].min_orf_aa=30;
        app.preferences.primer_small.min_melting_temperature=52.0;
        let doc=crate::doc::Document::from_bytes(bytes.clone()).unwrap();
        app.push_doc(doc,None);
        assert_eq!(app.seq_line_bp,Some(100));assert_eq!(app.orf_min_bp,90);
        let doc=app.doc.as_ref().unwrap();
        assert_eq!(hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap()).min_melting_temperature,52.0);
        assert_eq!(hx_formats::notes::field(&doc.session.notes_xml(),"Author").as_deref(),Some("Test author"));
        app.push_doc(crate::doc::Document::from_bytes(bytes.clone()).unwrap(),Some(std::path::PathBuf::from("existing.dna")));
        assert_eq!(app.doc.as_ref().unwrap().session.save_bytes().unwrap(),bytes);
    }
    #[test]
    fn hybridization_parameters_are_live_undoable_and_saved() {
        let bytes = hx_formats::write::NewFile::dna(b"ACGT".repeat(30)).to_bytes();
        let mut doc = crate::doc::Document::from_bytes(bytes).unwrap();
        let old = hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap_or_default());
        let mut params = old; params.min_melting_temperature = 55.0; params.min_continuous_match_len = 14;
        doc.session.apply_primer(hx_formats::primer_edit::PrimerEdit::Hybridization(params)).unwrap();
        assert_eq!(hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap()), params);
        let (bytes, _) = doc.session.save_bytes_with_history("x", "x", None).unwrap();
        let reopened = crate::doc::Document::from_bytes(bytes).unwrap();
        assert_eq!(hx_formats::annot::parse_hybridization_params(reopened.block_payload(5).unwrap()), params);
        doc.undo();
        assert_eq!(hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap()), old);
    }
    #[test]
    fn vector_template_opens_as_independent_document_without_changing_source() {
        let path = std::path::PathBuf::from("/tmp/opencode").join(format!("hx-vector-template-{}.dna", std::process::id()));
        let bytes = hx_formats::write::NewFile::dna(b"ACGTACGT".to_vec()).to_bytes();
        assert!(!path.exists());
        std::fs::write(&path, &bytes).unwrap();
        let mut app = super::HelixApp::blank();
        app.open_vector_template(&path, "Test vector");
        let doc = app.doc.as_ref().unwrap();
        assert!(doc.path.is_none(), "template must require Save As");
        assert_eq!(doc.sequence(), b"ACGTACGT");
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        std::fs::remove_file(path).unwrap();
    }
    use super::{feat_group, FeatGroup, HelixApp};

    #[test]
    fn view_options_is_above_canvas_buttons_even_when_they_overlap() {
        use eframe::egui;
        let ctx=egui::Context::default();crate::theme::apply(&ctx);
        let mut app=HelixApp::blank();
        let doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(100)).to_bytes()).unwrap();
        app.push_doc(doc,None);app.sidebar_open=false;app.show_view_options=true;
        let screen=egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(1800.0,1100.0));
        for _ in 0..3 {
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(screen),..Default::default()},|ui|{
                app.doc_canvas(ui,super::CanvasMode::Main);
                // Deliberately cover the floating view tabs and tool palette.
                let window=egui::Window::new("View Options").id(egui::Id::unique("view_options"))
                    .fixed_pos(egui::pos2(4.0,4.0)).fixed_size(egui::vec2(400.0,300.0)).show(ui,|ui|{ui.set_min_size(egui::vec2(400.0,300.0));ui.label(crate::i18n::text(ui.ctx(),"Options content"));});
                assert!(window.unwrap().response.rect.contains(egui::pos2(20.0,90.0)),"test point must lie inside the actual window");
            });output.textures_delta.clear();
        }
        let front=ctx.layer_id_at(egui::pos2(80.0,22.0)).expect("overlapping window must be visible");
        assert_ne!(front.id,egui::Id::unique("view_tabs_floating"),"canvas tabs cover the options window");
        let front=ctx.layer_id_at(egui::pos2(20.0,90.0)).unwrap();
        assert_ne!(front.id,egui::Id::unique("canvas_palette"),"canvas palette covers the options window");
    }

    #[test]
    fn paste_reverse_complement_uses_the_current_clipboard_event_and_is_undoable() {
        use eframe::{App,egui};
        let ctx=egui::Context::default();let mut app=HelixApp::blank();
        let doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"AAAA".to_vec()).to_bytes()).unwrap();
        app.push_doc(doc,None);app.cursor_bp=Some(2);app.paste="CCCC".into();
        app.request_paste(&ctx,true);
        let mut frame=eframe::Frame::_new_kittest();
        let mut output=ctx.run_ui(egui::RawInput{events:vec![egui::Event::Paste("AACG".into())],..Default::default()},|ui|app.ui(ui,&mut frame));output.textures_delta.clear();
        assert_eq!(app.doc.as_ref().unwrap().sequence(),b"AACGTTAA");
        app.do_undo();assert_eq!(app.doc.as_ref().unwrap().sequence(),b"AAAA");
    }

    #[test]
    fn save_and_paste_menu_contents_match_the_supplied_workflow() {
        use eframe::egui;
        let ctx=egui::Context::default();let mut app=HelixApp::blank();
        let doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGTACGT".to_vec()).to_bytes()).unwrap();app.push_doc(doc,None);
        for paste in [false,true] {
            let mut output=ctx.run_ui(egui::RawInput::default(),|ui|if paste {app.paste_button_menu(ui)}else{app.save_button_menu(ui)});
            let text=output.shapes.iter().filter_map(|s|if let egui::Shape::Text(t)=&s.shape{Some(t.galley.job.text.as_str())}else{None}).collect::<Vec<_>>().join("\n");output.textures_delta.clear();
            if paste {assert!(text.contains("Paste Reverse Complement"));assert!(text.contains("Ctrl+Shift+V"));assert!(!text.contains("Sequence Text"));}
            else{for label in ["Save As","Export DNA","Export Selected DNA","Export Selected Protein","Export Map","Export History","Export Feature Data","Export Primer Data"] {assert!(text.contains(label),"missing {label}: {text}");}}
        }
    }

    #[test]
    fn sequence_exports_have_correct_selected_bases_and_rna_alphabet() {
        let doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::of(hx_core::MoleculeKind::Rna,b"ACGTTGCA".to_vec()).to_bytes()).unwrap();
        let mut app=HelixApp::blank();app.push_doc(doc,None);app.selection=Some((2,6));
        assert!(app.exported_sequence_text(false).unwrap().ends_with("ACGUUGCA\n"));
        assert!(app.exported_sequence_text(true).unwrap().ends_with("GUUG\n"));
    }

    #[test]
    fn info_is_a_right_docked_panel_and_updates_map_label_through_undo() {
        use eframe::{App,egui};
        let ctx=egui::Context::default();let mut app=HelixApp::blank();
        let doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(200)).to_bytes()).unwrap();app.push_doc(doc,None);
        let mut frame=eframe::Frame::_new_kittest();app.show_info=true;
        for _ in 0..3 {
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(2000.0,1125.0))),..Default::default()},|ui|app.ui(ui,&mut frame));
            output.textures_delta.clear();
        }
        assert!(app.central_rect.unwrap().right()<1600.0);
        let d=app.doc.as_mut().unwrap();d.session.set_note("CustomMapLabel","Manual construct").unwrap();d.session.set_note("UseCustomMapLabel","1").unwrap();
        assert_eq!(d.map_name(),"Manual construct");d.undo();assert_ne!(d.map_name(),"Manual construct");d.redo();assert_eq!(d.map_name(),"Manual construct");
    }

    #[test]
    fn main_chrome_collapses_completely_and_restores_sidebar_with_native_map_title() {
        use eframe::{App,egui};
        let ctx=egui::Context::default();crate::theme::apply(&ctx);
        let mut new=hx_formats::write::NewFile::of(hx_core::MoleculeKind::Rna,b"ACGT".repeat(91));new.name=Some("Untitled".into());
        let mut doc=crate::doc::Document::from_bytes(new.to_bytes()).unwrap();doc.path=Some("445566.rna".into());
        let mut app=HelixApp::blank();app.push_doc(doc,None);app.map_layout=hx_render::MapLayout::Circular;
        let mut frame=eframe::Frame::_new_kittest();
        let screen=egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(2000.0,1125.0));
        let mut render=|app:&mut HelixApp,events:Vec<egui::Event>| {
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(screen),events,..Default::default()},|ui|app.ui(ui,&mut frame));
            output.textures_delta.clear();output
        };
        for _ in 0..3 {render(&mut app,vec![]);}
        let open=app.central_rect.unwrap();
        assert!((open.left()-crate::chrome::SIDEBAR_W).abs()<3.0,"{open:?}");
        assert!(open.top()>=crate::chrome::MENU_H+crate::chrome::DOC_TABS_H+crate::chrome::TOOLBAR_H);
        for expected in [false,true] {
            let pos=egui::pos2(16.0,crate::chrome::MENU_H+18.0);
            for pressed in [true,false] {render(&mut app,vec![egui::Event::PointerMoved(pos),egui::Event::PointerButton{pos,button:egui::PointerButton::Primary,pressed,modifiers:egui::Modifiers::NONE}]);}
            let output=render(&mut app,vec![]);assert_eq!(app.sidebar_open,expected);
            let canvas=app.central_rect.unwrap();
            assert!((canvas.left()-if expected {app.sidebar_width}else{0.0}).abs()<3.0,"{canvas:?}");
            let text=output.shapes.iter().filter_map(|s|if let egui::Shape::Text(t)=&s.shape{Some(t.galley.job.text.as_str())}else{None}).collect::<Vec<_>>().join("\n");
            assert!(text.contains("Untitled") && text.contains("445566.rna"));
            assert!(text.contains("Display"));assert!(!text.contains("Export"),"Export must not be a stray top-toolbar button");
        }
        let output=render(&mut app,vec![]);
        let display=output.shapes.iter().find_map(|s|if let egui::Shape::Text(t)=&s.shape {
            (t.galley.job.text=="Display").then(||t.pos+t.galley.size()*0.5)
        }else{None}).expect("Display menu trigger must be present");
        for pressed in [true,false] {render(&mut app,vec![egui::Event::PointerMoved(display),egui::Event::PointerButton{pos:display,button:egui::PointerButton::Primary,pressed,modifiers:egui::Modifiers::NONE}]);}
        let output=render(&mut app,vec![]);
        let text=output.shapes.iter().filter_map(|s|if let egui::Shape::Text(t)=&s.shape{Some(t.galley.job.text.as_str())}else{None}).collect::<Vec<_>>().join("\n");
        assert!(text.contains("Export current view") && text.contains("Features and feature labels"),"Display must open a functional menu: {text}");
    }

    #[test]
    fn native_document_tab_uses_standard_width_and_contains_close_hit_target() {
        use eframe::egui;
        let ctx=egui::Context::default();let mut targets=Vec::new();
        let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(1800.0,1100.0))),..Default::default()},|ui|{
            let (tab,close)=super::tab_button_with_close(ui,true,"445566.rna",0,Some(hx_core::MoleculeKind::Rna),230.0);
            assert_eq!(tab.rect.size(),egui::vec2(230.0,crate::chrome::DOC_TABS_H));
            assert!(tab.rect.contains_rect(close.rect));targets.push(close.rect.center());
        });output.textures_delta.clear();assert_eq!(targets.len(),1);
    }

    #[test]
    fn active_document_tab_has_no_blue_lower_edge_but_view_tabs_keep_the_underline() {
        use eframe::egui;
        let ctx=egui::Context::default();
        let mut rect=egui::Rect::NOTHING;
        let mut output=ctx.run_ui(egui::RawInput::default(),|ui|{
            rect=super::tab_button_with_close(ui,true,"active.dna",0,Some(hx_core::MoleculeKind::Dna),196.0).0.rect;
        });output.textures_delta.clear();
        for shape in output.shapes {if let egui::Shape::LineSegment{points,stroke}=shape.shape {
            assert!(!(points[0].y-points[1].y).abs().lt(&0.1) || (points[0].y-rect.bottom()).abs()>1.0 || stroke.color!=crate::sidebar::EDGE,
                "active document has a blue lower edge");
        }}
    }

    #[test]
    fn cloned_product_save_target_never_reuses_the_source_path() {
        let dir = std::env::temp_dir().join(format!(
            "hx-clone-save-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let source = dir.join("source-vector.dna");
        std::fs::write(&source, b"SOURCE").unwrap();
        let mut app = HelixApp::blank();
        app.folder = Some(dir.clone());
        app.save_path = source.display().to_string();
        let bytes = hx_formats::write::NewFile::dna(b"ACGT".repeat(40)).to_bytes();
        app.adopt_cloned_product(bytes, "Assembled.dna".into());
        let target = std::path::PathBuf::from(&app.save_path);
        assert_eq!(
            target,
            dir.join("Assembled.dna"),
            "Save/Save As after Clone must default to the product name"
        );
        assert_eq!(app.doc.as_ref().unwrap().display_name(), "Assembled.dna");
        assert_eq!(app.tabs.len(), 1);
        // 真正写到默认落点：产物可重新解析，源载体一个字节都没动。
        let (written, _) = app.doc.as_mut().unwrap().save_to(&target).unwrap();
        assert!(written > 0);
        assert_eq!(std::fs::read(&source).unwrap(), b"SOURCE");
        let reloaded = crate::doc::Document::from_bytes(std::fs::read(&target).unwrap()).unwrap();
        assert_eq!(reloaded.sequence(), app.doc.as_ref().unwrap().sequence());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unsaved_tabs_get_their_own_save_target() {
        let dir = std::env::temp_dir().join(format!(
            "hx-unsaved-tabs-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let mut app = HelixApp::blank();
        app.folder = Some(dir.clone());
        for name in ["First.dna", "Second.dna"] {
            let mut doc =
                crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(10)).to_bytes())
                    .unwrap();
            doc.name_hint = Some(name.to_string());
            app.push_doc(doc, None);
            assert_eq!(
                app.save_path,
                dir.join(name).display().to_string(),
                "each unsaved tab must own its default save target"
            );
        }
        app.activate_tab(0);
        assert_eq!(
            app.save_path,
            dir.join("First.dna").display().to_string(),
            "switching back must restore the first tab's target, not the second's"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn internal_save_as_dialog_requires_a_confirmed_path() {
        let dir = std::env::temp_dir().join(format!(
            "hx-save-dialog-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let mut app = HelixApp::blank();
        app.folder = Some(dir.clone());
        let mut doc =
            crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(10)).to_bytes())
                .unwrap();
        doc.name_hint = Some("Product.dna".into());
        app.push_doc(doc, None);
        app.begin_save_dialog();
        assert!(app.show_save_dialog);
        assert_eq!(app.save_as_path, dir.join("Product.dna").display().to_string());
        app.save_as_path = dir.join("Renamed.dna").display().to_string();
        app.confirm_save_dialog();
        assert!(!app.show_save_dialog);
        assert!(dir.join("Renamed.dna").is_file(), "{}", app.status);
        assert!(!dir.join("Product.dna").exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn display_preferences_roundtrip_restores_toggles() {
        let mut a=HelixApp::blank();
        a.show_sites=false;a.show_primers=false;a.show_colors=false;a.aa_one_letter=true;
        let names=a.display_preferences();
        let mut b=HelixApp::blank();
        b.apply_display_preferences(&names);
        assert_eq!(b.display_preferences(),names);
        assert!(!b.show_sites&&!b.show_primers&&!b.show_colors);
        assert!(b.aa_one_letter&&b.show_features);
        b.apply_display_preferences(&["features".into(),"unknown-future-key".into()]);
        assert_eq!(b.display_preferences(),vec!["features".to_string()]);
    }

    #[test]
    fn collapsed_sidebar_still_keeps_the_full_width_tab_row_separator() {
        use eframe::{egui,App};
        let ctx=egui::Context::default();crate::theme::apply(&ctx);
        let mut app=HelixApp::blank();app.sidebar_open=false;
        for i in 0..2 {
            let mut doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(40)).to_bytes()).unwrap();
            doc.name_hint=Some(format!("collapsed-{i}.dna"));app.push_doc(doc,None);
        }
        app.activate_tab(1);let mut frame=eframe::Frame::_new_kittest();
        for _ in 0..3 {
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(2000.0,1125.0))),..Default::default()},|ui|app.ui(ui,&mut frame));
            output.textures_delta.clear();
            let bottom=crate::chrome::MENU_H+crate::chrome::DOC_TABS_H-0.5;
            let lines=output.shapes.iter().filter_map(|s|if let egui::Shape::LineSegment{points,stroke}=&s.shape {
                (stroke.color==crate::sidebar::EDGE).then_some(*points)
            }else{None}).collect::<Vec<_>>();
            let spans:Vec<(f32,f32)>=lines.iter()
                .filter(|p|(p[0].y-bottom).abs()<0.1 && (p[1].y-bottom).abs()<0.1)
                .map(|p|(p[0].x.min(p[1].x),p[0].x.max(p[1].x))).collect();
            assert!(spans.iter().any(|(a,_)|*a<1.0),"collapsed sidebar must keep the left end of the separator: {spans:?}");
            assert!(spans.iter().any(|(_,b)|*b>1900.0),"separator must reach the right edge: {spans:?}");
            assert!(!spans.iter().any(|(a,b)|*a<269.0&&*b>271.0),"separator must not cross the active tab: {spans:?}");
        }
    }

    #[test]
    fn document_row_keeps_full_blue_separator_with_first_or_middle_tab_selected() {
        use eframe::{egui,App};
        for selected in [0,1] {
            let ctx=egui::Context::default();crate::theme::apply(&ctx);
            let mut app=HelixApp::blank();
            for i in 0..3 {
                let mut doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(40)).to_bytes()).unwrap();
                doc.name_hint=Some(format!("document-{i}.dna"));app.push_doc(doc,None);
            }
            app.activate_tab(selected);let mut frame=eframe::Frame::_new_kittest();
            for _ in 0..3 {
                let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(2000.0,1125.0))),..Default::default()},|ui|app.ui(ui,&mut frame));
                output.textures_delta.clear();
                let x=app.sidebar_width+selected as f32*196.0;
                let bottom=crate::chrome::MENU_H+crate::chrome::DOC_TABS_H-0.5;
                let lines=output.shapes.iter().filter_map(|s|if let egui::Shape::LineSegment{points,stroke}=&s.shape {
                    (stroke.color==crate::sidebar::EDGE).then_some(*points)
                }else{None}).collect::<Vec<_>>();
                assert!(lines.iter().any(|p|(p[0].y-bottom).abs()<0.1 && (p[1].y-bottom).abs()<0.1 && p[1].x>1900.0),"separator must continue through unused tab-row space: expected y={bottom}, {lines:?}");
                assert!(lines.iter().any(|p|(p[0].x-(x+0.5)).abs()<1.0 && (p[1].x-p[0].x).abs()<0.1 && p[1].y>=bottom-0.1),"active left edge missing: selected {selected}");
                assert!(!lines.iter().any(|p|(p[0].y-bottom).abs()<0.1 && (p[1].y-bottom).abs()<0.1 && p[0].x<x+8.0 && p[1].x>x+180.0),"separator must not span the active tab");
            }
        }
    }

    /// SG keeps a vertical rule on the navigation/tab boundary for the whole
    /// tab row, so it meets the sidebar edge below instead of vanishing with
    /// the first tab's own border. It is blue while Files is open and grey
    /// when the panel is collapsed.
    #[test]
    fn document_row_keeps_a_navigation_boundary_joined_to_the_separator() {
        use eframe::{egui,App};
        for open in [true,false] {
            let ctx=egui::Context::default();crate::theme::apply(&ctx);
            let mut app=HelixApp::blank();app.sidebar_open=open;
            for i in 0..2 {
                let mut doc=crate::doc::Document::from_bytes(hx_formats::write::NewFile::dna(b"ACGT".repeat(40)).to_bytes()).unwrap();
                doc.name_hint=Some(format!("boundary-{i}.dna"));app.push_doc(doc,None);
            }
            app.activate_tab(1);let mut frame=eframe::Frame::_new_kittest();
            let expected=if open {crate::sidebar::EDGE}else{egui::Color32::from_rgb(0xd9,0xd9,0xd9)};
            let nav=crate::chrome::navigation_width(open,app.sidebar_width)-1.5;
            let mut found=false;
            for _ in 0..3 {
                let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(2000.0,1125.0))),..Default::default()},|ui|app.ui(ui,&mut frame));
                output.textures_delta.clear();
                for s in &output.shapes {
                    let egui::Shape::Rect(r)=&s.shape else{continue};
                    if r.fill!=expected||r.rect.width()>1.5 {continue;}
                    // 1 px wide, starts 2 px below the menu and reaches the
                    // horizontal blue separator row.
                    assert!((r.rect.left()-nav).abs()<1.0,"boundary x {} != nav {nav}",r.rect.left());
                    assert!(r.rect.top()<=crate::chrome::MENU_H+2.5,"boundary must start at the row top: {:?}",r.rect);
                    assert!(r.rect.bottom()>=crate::chrome::MENU_H+crate::chrome::DOC_TABS_H-1.0,"boundary must meet the separator: {:?}",r.rect);
                    found=true;
                }
            }
            assert!(found,"missing navigation boundary (sidebar_open={open})");
        }
    }

    /// The `+` is painted, not typeset: both arms share one centre and the
    /// cross uses a regular-weight stem and stays centred in the row.
    #[test]
    fn new_document_plus_is_a_centred_regular_cross() {
        use eframe::egui;
        let ctx=egui::Context::default();crate::theme::apply(&ctx);
        let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(400.0,120.0))),..Default::default()},|ui|{
            egui::CentralPanel::default().show(ui,|ui|{
                crate::icons::plus_menu(ui,|_ui|{});
            });
        });
        output.textures_delta.clear();
        let rects:Vec<egui::Rect>=output.shapes.iter().filter_map(|s|if let egui::Shape::Rect(r)=&s.shape {
            (r.fill==crate::theme::TEXT).then_some(r.rect)
        }else{None}).collect();
        let bars:Vec<&egui::Rect>=rects.iter().filter(|r|r.width().min(r.height())>=1.5).collect();
        assert_eq!(bars.len(),2,"expected two plus bars, got {rects:?}");
        assert!((bars[0].center()-bars[1].center()).length()<0.01,"bars must share a centre: {bars:?}");
        let (a,b)=(bars[0],bars[1]);
        let (h,v)=if a.width()>a.height() {(a,b)}else{(b,a)};
        assert!((h.width()-10.0).abs()<0.1 && (h.height()-1.5).abs()<0.1,"horizontal bar {h:?}");
        assert!((v.height()-10.0).abs()<0.1 && (v.width()-1.5).abs()<0.1,"vertical bar {v:?}");
    }

    #[test]
    fn empty_document_row_has_no_placeholder_and_plus_starts_after_navigation() {
        use eframe::egui;
        for open in [true,false] {
            let ctx=egui::Context::default();crate::theme::apply(&ctx);
            let mut app=HelixApp::blank();app.sidebar_open=open;
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(1800.0,1100.0))),..Default::default()},|ui|{
                app.navigation_row(ui);
            });
            output.textures_delta.clear();
            assert!(!output.shapes.iter().any(|s|matches!(&s.shape,egui::Shape::Text(t) if t.galley.job.text.contains("No file open"))));
            let bar=output.shapes.iter().find_map(|s|match &s.shape {
                egui::Shape::Rect(r) if r.fill==crate::theme::TEXT && (r.rect.width()-10.0).abs()<0.01=>Some(r.rect),
                _=>None,
            }).expect("new/open plus remains visible without documents");
            let nav=crate::chrome::navigation_width(open,app.sidebar_width);
            assert!(bar.center().x>=nav && bar.center().x<nav+30.0,"plus should immediately follow navigation: {bar:?}, nav={nav}");
        }
    }

    #[test]
    fn empty_main_canvas_is_blank_but_cloning_keeps_its_source_prompt() {
        use eframe::egui;
        for mode in [super::CanvasMode::Main,super::CanvasMode::Cloning] {
            let ctx=egui::Context::default();crate::theme::apply(&ctx);
            let mut app=HelixApp::blank();
            let mut output=ctx.run_ui(egui::RawInput::default(),|ui|{
                app.doc_canvas(ui,mode);
            });
            output.textures_delta.clear();
            let text=output.shapes.iter().filter_map(|s|match &s.shape {
                egui::Shape::Text(t)=>Some(t.galley.job.text.as_str()),
                _=>None,
            }).collect::<Vec<_>>();
            if mode==super::CanvasMode::Main {
                assert!(text.is_empty(),"empty main canvas should not show placeholder text: {text:?}");
            } else {
                assert!(text.contains(&"Please choose a source file."));
            }
        }
    }

    #[test]
    fn document_title_uses_one_text_shape_and_a_matching_real_bold_face() {
        use eframe::egui;
        for active in [true,false] {
            let ctx=egui::Context::default();crate::theme::apply(&ctx);
            let installed=crate::fonts::install(&ctx);
            let mut output=ctx.run_ui(egui::RawInput::default(),|ui|{
                super::tab_button_with_close(ui,active,"Title 123",0,None,196.0);
            });
            output.textures_delta.clear();
            let titles=output.shapes.iter().filter_map(|s|match &s.shape {
                egui::Shape::Text(t) if t.galley.job.text=="Title 123"=>Some(t),
                _=>None,
            }).collect::<Vec<_>>();
            assert_eq!(titles.len(),1,"do not overpaint text to fake bold");
            let font=&titles[0].galley.job.sections[0].format.font_id;
            assert_eq!(font.size,13.0);
            assert_eq!(*font,crate::fonts::document_title_font(&ctx,active));
            if active && installed==Some("/usr/share/fonts/noto/NotoSans-Regular.ttf")
                && std::path::Path::new("/usr/share/fonts/noto/NotoSans-Bold.ttf").is_file() {
                assert_eq!(font.family,egui::FontFamily::Name("hb-sans-bold".into()));
            }
        }
    }

    /// 蛋白特征的类别判定（左侧工具条三段开关的后端依据）。
    #[test]
    fn feat_group_classifies_protein_feature_types() {
        // oracle：KPYK1_ECOLI.prot 里的 `type` 取值
        assert_eq!(feat_group("Region"), FeatGroup::Region);
        assert_eq!(feat_group("Conflict"), FeatGroup::Region);
        assert_eq!(feat_group("Site"), FeatGroup::SiteHet);
        assert_eq!(feat_group("site"), FeatGroup::SiteHet); // 大小写不敏感
        assert_eq!(feat_group("Het"), FeatGroup::SiteHet);
        assert_eq!(feat_group("Bond"), FeatGroup::Bond);
    }

    /// 构建时间戳的日历换算：跨闰日、跨世纪都要对（About 里用来识别旧二进制）。
    #[test]
    fn utc_stamp_matches_known_instants() {
        assert_eq!(HelixApp::utc_stamp(0), "1970-01-01 00:00 UTC");
        assert_eq!(HelixApp::utc_stamp(951_782_400), "2000-02-29 00:00 UTC");
        assert_eq!(HelixApp::utc_stamp(951_782_460), "2000-02-29 00:01 UTC");
        // 2025-10-18 13:05:24 UTC（`date -u -d @1760792724`）
        assert_eq!(HelixApp::utc_stamp(1_760_792_724), "2025-10-18 13:05 UTC");
    }

    /// 拖动重排页签：顺序变化，但活动页签仍然指向同一个文档。
    #[test]
    fn tab_reorder_keeps_the_active_document() {
        let mut app = HelixApp::blank();
        for n in [10usize, 20, 30, 40] {
            let doc = crate::doc::Document::from_bytes(
                hx_formats::write::NewFile::dna(vec![b'A'; n]).to_bytes(),
            )
            .unwrap();
            app.push_doc(doc, None);
        }
        app.activate_tab(2);
        let lens = |app: &HelixApp| -> Vec<usize> {
            (0..app.tabs.len())
                .map(|i| {
                    if i == app.active {
                        app.doc.as_ref().unwrap().len()
                    } else {
                        app.tabs[i].doc.as_ref().unwrap().len()
                    }
                })
                .collect()
        };
        assert_eq!(lens(&app), vec![10, 20, 30, 40]);
        // 把第一个页签拖到第四个之前。
        app.move_tab(0, 3);
        assert_eq!(lens(&app), vec![20, 30, 10, 40]);
        assert_eq!(app.active, 1);
        assert_eq!(app.doc.as_ref().unwrap().len(), 30);
        // 拖动活动页签本身。
        app.move_tab(1, 0);
        assert_eq!(lens(&app), vec![30, 20, 10, 40]);
        assert_eq!(app.active, 0);
        assert_eq!(app.doc.as_ref().unwrap().len(), 30);
    }
}
