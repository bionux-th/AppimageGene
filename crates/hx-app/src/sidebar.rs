//! 左侧「Files」侧栏（对齐 SnapGene 的 **folder view**）。
//!
//! 几何全部来自 **folder view 的 1280×1024 oracle**（`docs/app-shots/156-snapgene-sidebar-folder-view.png`，
//! 用户真实配置）：白底、右缘 1 px 亮蓝线、标题 16.5 px、行距 28 px、行文字 15 px、
//! 缩进 17 px、文件图标 15 px 圆环；选中行是**两侧内缩 16 px 的圆角蓝底 `#0066EC`**。
//!
//! **没有打开文件夹时**（空状态）：
//!
//! ```text
//! [◧]                                        ⋯   ← 折叠 / 溢出菜单
//! [ Open File… ]                                 ← 蓝底主按钮（系统原生对话框）
//! FOLDER
//! [ Open Folder… ]                          ↻    ← 系统原生文件夹对话框 + 刷新
//! [ /path/to/folder                        ]     ← 路径框，回车跳转
//! Tip: you can also drop a file onto the window.
//! No folder opened
//! ```
//!
//! **打开文件夹之后**（用户裁定 2026-10-18：`Open File…` / `Open Folder…` 不再出现）：
//!
//! ```text
//! [◧]                                     ↻        ← 头行（折叠 / 刷新）
//! Sample project                              ⋯    ← 根文件夹名 + 溢出菜单（SG 的 ⋮）
//! [ 🔍 Search and filter files             ]       ← 过滤框（26 px、斜体灰提示）
//!   ▾ 📁 4. Sample Files                           ← 子文件夹：三角 + **文件夹图标**
//!     ▾ 📁 DNA files                               ← 展开内容**就地**缩进一层
//!       ◯ pUC19.dna                                ← 文件名**带扩展名**；当前文档蓝底白字
//!   ▤ 输出的酶组.txt                                ← 其它文件：灰页图标
//! ───────────────────────────────────────
//! OTHER OPEN FILES
//!   ◯ pIB2-SEC13-mEGFP.dna                         ← **只列根文件夹之外**的打开文档
//! ```
//!
//! 四条与 SnapGene 对齐的行为：
//!
//! 0. **不隐藏扩展名**、文件夹画 `📁`、图标按类型画彩色圆环（原生 `.dna/.rna/.prot` 是
//!    圆环 + 双螺旋，FASTA/GenBank 是亮蓝圆环 + 折角页）—— 用户 2026-10-19 裁定。
//!    `OTHER OPEN FILES` 只列不在根文件夹里的文档（SG 把文件夹内的都画在树上）。
//!
//!
//! 1. **打开文件夹后不再显示** `Open File…` / `FOLDER` / `Open Folder…` / 路径框 /
//!    拖放提示 —— 那一段只属于空状态；打开文件夹后顶上换成文件树的过滤框。
//!    换文件夹走 `⋯ ▸ Open Folder…` 或 `File ▸ Open Folder…`。
//! 2. **根目录固定**：打开文件、切标签页都**不会**改侧栏根目录，只有显式
//!    `Open Folder…` / 路径框回车才换（见 [`crate::app::HelixApp::set_folder`]）。
//! 3. **就地展开**：展开的子文件夹，内容紧跟它自己（深一层缩进）。旧实现把待展开的
//!    目录压进一个 `Vec` 当栈、边生成边渲染，展开目录的内容会掉到**整个列表末尾**，
//!    看起来像挂在别的文件夹下面（用户 2026-10-18 的反馈）。
//!
//! 打开文档时会自动展开到它的祖先文件夹（SnapGene 实测截图 26 的行为）。
//! 打开文件/文件夹一律走**系统原生对话框**（KDE 上是 Dolphin 那套，见 [`crate::sysopen`]），
//! 内置路径小窗口只在系统没有对话框工具时兜底。

use std::path::PathBuf;

use eframe::egui::{self, Color32, Rect, Vec2};

use crate::app::HelixApp;
use crate::theme;

/// 行高（SnapGene 的列表行）。
///
/// oracle 实测（`sg-sb.png`，1280×1024，用户真实配置的 **folder view**）：行距 **28 px**、
/// 行文字 ≈ **15 px**（中文同名文件 `输出的酶组.txt` 在 SG 宽 84 px、我们 13 px 字号时 69 px）。
const ROW_H: f32 = 28.0;
/// 列表行字号（SG folder view 实测 ≈15 px）。
const ROW_FONT: f32 = 13.0;
/// 根文件夹标题字号：SG 实测大写字母 `F` 高 12 px（≈16.5 px 字号）。
const TITLE_FONT: f32 = 15.0;
/// 段标题（`OTHER OPEN FILES`）字号：SG 实测 cap 高 9 px（≈12 px 粗体）。
const SECTION_FONT: f32 = 12.0;
/// 文档/文件图标尺寸（SG 实测 15×15）。
const ICON: f32 = 15.0;
/// 每一层缩进（SG 实测：depth0 三角 x≈22、depth1 三角 x≈39 → 17 px）。
const INDENT: f32 = 17.0;
/// 展开三角占的列宽（文件行留空，保证同级名字对齐）。
const EXPAND_W: f32 = 20.0;
/// 文件图标占的列宽（文件夹行留空，理由同上）。
const ICON_SLOT: f32 = 16.0;
/// 行左内边距（SG：depth0 的三角中心在 22 px 处）。
const ROW_PAD: f32 = 14.0;
/// 选中行的圆角蓝底从面板两侧各内缩多少（SG 实测 16 px）。
const SEL_INSET: f32 = 16.0;
/// 侧栏底色：SG folder view 是**纯白**（`sg-sb.png` 取样 255）。
pub const BG: Color32 = Color32::WHITE;
/// 侧栏右缘 1 px 亮蓝分隔线（SG 实测 `(43,189,255)`）。
pub const EDGE: Color32 = Color32::from_rgb(0x2b, 0xbd, 0xff);
/// 树/`OTHER OPEN FILES` 选中行底色（SG 实测 `(0,102,236)`）。
pub const SEL_BLUE: Color32 = Color32::from_rgb(0x00, 0x66, 0xec);
/// 搜索框边框（SG 实测 `#A9A9A9`）。
const SEARCH_BORDER: Color32 = Color32::from_rgb(0xa9, 0xa9, 0xa9);
/// 搜索框提示字 / 放大镜灰（SG 斜体占位符）。
const HINT_GRAY: Color32 = Color32::from_rgb(0x8c, 0x8c, 0x8c);
/// 文件夹图标描边（SG 实测 `#333333`）。
const FOLDER_INK: Color32 = Color32::from_rgb(0x33, 0x33, 0x33);
/// 序列文本格式（FASTA/GenBank…）图标色（SG 实测 `(43,189,255)`）。
const TEXT_SEQ: Color32 = Color32::from_rgb(0x2b, 0xbd, 0xff);
/// `OTHER OPEN FILES` 上方的分隔线（SG 实测 `#D6D6D6`）。
const RULE: Color32 = Color32::from_rgb(0xd6, 0xd6, 0xd6);

/// `OTHER OPEN FILES` 的入选规则：**不在当前根文件夹里**的单独打开文档。
///
/// oracle `sg-sb.png`：根文件夹 `Fw__R1-4_sequences` 里的 4 个打开文档都显示在树上，
/// 这一段只列出文件夹外的 `pUC19.dna`。
pub fn is_outside_root(path: &std::path::Path, root: Option<&std::path::Path>) -> bool {
    match root {
        Some(r) => !path.starts_with(r),
        None => true,
    }
}

/// 侧栏显示名。**不再去扩展名**：SnapGene 的 folder view 里显示带后缀的完整文件名
/// （`sg-sb.png`：`2333.dna` / `pET28b-Rep1.dna` / `输出的酶组.txt`），用户 2026-10-19
/// 明确要求「不希望我们软件隐藏文件后缀名」。
pub fn display_name(name: String) -> String {
    name
}

/// 单行截断（SG 的 `Screenshot 20…9.41.16 AM.png`）：放不下时**在中间**省略，
/// 保住文件名的开头和**扩展名**（用户明确要求后缀可见）。
fn ellipsize(
    p: &egui::Painter,
    text: &str,
    font: egui::FontId,
    col: Color32,
    max_w: f32,
) -> std::sync::Arc<egui::Galley> {
    let layout = |s: &str| p.layout_no_wrap(s.to_string(), font.clone(), col);
    let full = layout(text);
    if full.size().x <= max_w {
        return full;
    }
    let chars: Vec<char> = text.chars().collect();
    let n = chars.len();
    let (mut lo, mut hi) = (0usize, n);
    let mut best: String = "…".to_string();
    while lo < hi {
        let mid = (lo + hi).div_ceil(2);
        let head = (mid * 55 / 100).min(n);
        let tail = (mid - head).min(n - head);
        let s: String = chars[..head].iter().collect::<String>()
            + "…"
            + &chars[n - tail..].iter().collect::<String>();
        if layout(&s).size().x <= max_w {
            best = s;
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    layout(&best)
}

/// 原生文档图标（`.dna`/`.rna`/`.prot`）：**彩色圆环 + 双螺旋**。
///
/// oracle `sg-sb.png` 8 倍放大：SG folder view 里 `.dna` = 红环 `#DB5550` + 粉红双螺旋，
/// `.rna` = 绿环，蛋白 = 青环；选中行（蓝底）整套转白。
fn doc_badge(p: &egui::Painter, rect: Rect, kind: hx_core::MoleculeKind, white: bool) {
    crate::branding::document_icon(p,rect,kind,white);
}

#[allow(dead_code)]
fn legacy_doc_badge(p: &egui::Painter, rect: Rect, kind: hx_core::MoleculeKind, white: bool) {
    if kind==hx_core::MoleculeKind::Rna {
        if !white {crate::theme::kind_icon(p,rect,kind,true);}else{
            // Selected native RNA keeps its single-strand glyph, never a DNA helix.
            let points=(0..=14).map(|i|{let t=i as f32/14.0;egui::pos2(rect.center().x+(t*std::f32::consts::TAU).sin()*rect.width()*0.20,
                rect.top()+2.0+t*(rect.height()-4.0))}).collect();
            p.add(egui::Shape::line(points,egui::Stroke::new(1.4,Color32::WHITE)));
        }
        return;
    }
    let (r, g, b) = kind.icon_rgb();
    let col = if white {
        Color32::WHITE
    } else {
        Color32::from_rgb(r, g, b)
    };
    let c = rect.center();
    let rad = rect.width() * 0.5 - 0.6;
    // Native documents use the sequence ribbon outline, not a topology circle.
    let outline=vec![egui::pos2(rect.left()+3.0,rect.top()+0.5),egui::pos2(rect.right()-3.0,rect.top()+0.5),
        egui::pos2(rect.right()-4.5,c.y),egui::pos2(rect.right()-1.0,rect.bottom()-0.5),
        egui::pos2(rect.left()+3.0,rect.bottom()-0.5),egui::pos2(rect.left()+4.5,c.y),egui::pos2(rect.left()+1.0,rect.top()+0.5)];
    p.add(egui::Shape::closed_line(outline,egui::Stroke::new(1.1,col)));
    // 双螺旋：两条反相正弦，竖直方向
    let amp = rect.width() * 0.17;
    let steps = 12;
    let (y0, y1) = (c.y - rad * 0.78, c.y + rad * 0.78);
    for phase in [0.0f32, std::f32::consts::PI] {
        let pts: Vec<egui::Pos2> = (0..=steps)
            .map(|i| {
                let t = i as f32 / steps as f32;
                egui::pos2(
                    c.x + amp * (t * 9.0 + phase).sin(),
                    y0 + t * (y1 - y0),
                )
            })
            .collect();
        p.add(egui::Shape::line(pts, egui::Stroke::new(1.1, col)));
    }
}

/// 序列文本格式（FASTA / GenBank / `.fa` / `.gbk` …）：**亮蓝圆环 + 折角页**。
/// oracle `sg-sb.png`：SG 把这类文件与 SnapGene 原生文档区分开（不是红/绿螺旋）。
fn text_seq_badge(p: &egui::Painter, rect: Rect, white: bool) {
    let col = if white {
        Color32::WHITE
    } else {
        TEXT_SEQ
    };
    let c = rect.center();
    let rad = rect.width() * 0.5 - 0.6;
    p.circle_stroke(c, rad, egui::Stroke::new(1.3, col));
    // 折角页（与其它文件的灰页同一个形状，只是染色）；略向左，别和圆环右缘糊在一起
    let page = Rect::from_center_size(
        egui::Pos2::new(c.x - 1.0, c.y),
        Vec2::new(rect.width() * 0.44, rect.height() * 0.58),
    );
    p.rect_stroke(
        page,
        1.0,
        egui::Stroke::new(1.0, col),
        egui::StrokeKind::Inside,
    );
    p.line_segment(
        [
            egui::pos2(page.right() - page.width() * 0.45, page.top()),
            egui::pos2(page.right(), page.top() + page.height() * 0.35),
        ],
        egui::Stroke::new(1.0, col),
    );
}

/// `.pdf` 图标：SG 是**灰页 + 压在左下角的红色小块**（oracle `ic-pdf.png` 8 倍放大：
/// 16×16 灰页、折角在右上，红色 `#EB2F3B` 方块 10×9 盖住左半，里面有白色音符）。
fn pdf_icon(p: &egui::Painter, rect: Rect) {
    file_icon(p, rect);
    let red = Color32::from_rgb(0xeb, 0x2f, 0x3b);
    let badge = Rect::from_min_size(
        egui::Pos2::new(rect.left() + rect.width() * 0.16, rect.top() + rect.height() * 0.30),
        Vec2::new(rect.width() * 0.58, rect.height() * 0.55),
    );
    p.rect_filled(badge, 0.8, red);
    // 白色音符：实心符头 + 符干
    let head = egui::Pos2::new(badge.center().x - 0.6, badge.center().y + 1.8);
    p.circle_filled(head, rect.width() * 0.09, Color32::WHITE);
    p.line_segment(
        [
            egui::Pos2::new(head.x + rect.width() * 0.08, head.y),
            egui::Pos2::new(head.x + rect.width() * 0.08, badge.top() + rect.height() * 0.10),
        ],
        egui::Stroke::new(1.0, Color32::WHITE),
    );
}

/// 其它文件（`.pdf`/`.png`/`.txt`…）：灰页图标（手绘，避免字体缺字）。
pub(crate) fn file_icon(p: &egui::Painter, rect: Rect) {
    let page = rect.shrink2(Vec2::new(1.5, 0.5));
    let bg = Color32::from_rgb(0xf2, 0xf3, 0xf5);
    p.rect_filled(page, 1.5, bg);
    p.rect_stroke(
        page,
        1.5,
        egui::Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    // Unknown/unsupported files are deliberately blank, not text documents.
}


/// 整行列表项：返回 (响应, 行矩形)。（占满剩余宽度，便于整行高亮）
fn list_row(ui: &mut egui::Ui) -> (egui::Response, Rect) {
    let w = ui.available_width();
    let (rect, resp) = ui.allocate_exact_size(Vec2::new(w, ROW_H), egui::Sense::click());
    (resp, rect)
}

/// 画一行的底与内容（先 allocate 再画，保证 hover 反馈）。
///
/// `kind` = 原生分子文档类型（画类型色圆环 + 双螺旋）；`text_seq` = FASTA/GenBank 这类
/// 序列文本（亮蓝圆环 + 折角页）；两者都没有 ⇒ 灰页图标。
///
/// 选中样式对齐 SG folder view（`sg-sb.png`）：**两侧各内缩 16 px 的圆角蓝底**
/// `#0066EC`，文字与图标转白（不是整行通栏）。
#[allow(clippy::too_many_arguments)]
fn paint_row(
    ui: &egui::Ui,
    rect: Rect,
    selected: bool,
    hovered: bool,
    kind: Option<hx_core::MoleculeKind>,
    text_seq: bool,
    label: &str,
    dim: bool,
) {
    let p = ui.painter();
    let body = Rect::from_min_max(
        egui::Pos2::new(rect.left() + SEL_INSET, rect.top() + 1.0),
        egui::Pos2::new(rect.right() - SEL_INSET, rect.bottom() - 1.0),
    );
    if selected {
        p.rect_filled(body, 5.0, SEL_BLUE);
    } else if hovered {
        p.rect_filled(body, 5.0, theme::ACCENT_BG);
    }
    // SG：图标列与树里的文件行对齐（x=36，文本 57）
    let mut x = rect.left() + SEL_INSET + 20.0;
    let r = Rect::from_center_size(
        egui::Pos2::new(x + ICON * 0.5, rect.center().y),
        Vec2::splat(ICON),
    );
    let _=(kind,text_seq);
    crate::branding::file_icon(p,r,label,selected);
    x += ICON + 6.0;
    let text_col = if selected {
        Color32::WHITE
    } else if dim {
        theme::TEXT_DIM
    } else {
        theme::TEXT
    };
    let avail = (rect.right() - SEL_INSET - 4.0 - x).max(20.0);
    let galley = ellipsize(
        p,
        label,
        egui::FontId::proportional(ROW_FONT),
        text_col,
        avail,
    );
    p.galley(
        egui::Pos2::new(x, rect.center().y - galley.size().y * 0.5),
        galley,
        text_col,
    );
}

/// 段标题（SnapGene 的 `OTHER OPEN FILES`：12 px 粗体灰字大写）。
fn section_label(ui: &mut egui::Ui, text: &str) {
    let text=crate::i18n::text(ui.ctx(),text);
    ui.horizontal(|ui| {
        ui.add_space(SEL_INSET);
        ui.label(
            egui::RichText::new(text)
                .size(SECTION_FONT)
                .color(Color32::from_rgb(0x4d, 0x4d, 0x4d))
                .strong(),
        );
    });
}

/// 侧栏搜索框（SG 实测：26 px 高、圆角 4、1 px `#A9A9A9` 边、左侧放大镜、
/// 斜体灰提示 `Search and filter files` / `filter open files`）。
fn search_box(ui: &mut egui::Ui, text: &mut String, hint: &str) -> egui::Response {
    ui.horizontal(|ui| {
        ui.add_space(SEL_INSET);
        let w = ui.available_width() - SEL_INSET;
        let (rect, _) = ui.allocate_exact_size(Vec2::new(w, 26.0), egui::Sense::hover());
        let p = ui.painter();
        p.rect_filled(rect, 4.0, Color32::WHITE);
        p.rect_stroke(
            rect,
            4.0,
            egui::Stroke::new(1.0, SEARCH_BORDER),
            egui::StrokeKind::Inside,
        );
        // 放大镜
        let c = egui::Pos2::new(rect.left() + 14.0, rect.center().y - 0.5);
        p.circle_stroke(c, 4.0, egui::Stroke::new(1.2, HINT_GRAY));
        p.line_segment(
            [c + Vec2::new(2.9, 2.9), c + Vec2::new(6.0, 6.0)],
            egui::Stroke::new(1.3, HINT_GRAY),
        );
        let inner = Rect::from_min_max(
            egui::Pos2::new(rect.left() + 25.0, rect.top() + 2.0),
            egui::Pos2::new(rect.right() - 6.0, rect.bottom() - 2.0),
        );
        let mut child = ui.new_child(egui::UiBuilder::new().max_rect(inner));
        child.add_sized(
            inner.size(),
            egui::TextEdit::singleline(text).frame(egui::Frame::NONE).hint_text(
                egui::RichText::new(hint)
                    .italics()
                    .size(14.0)
                    .color(HINT_GRAY),
            ),
        )
    })
    .inner
}

/// 刷新图标（手绘圆弧 + 箭头，避免默认字体缺 `↻` 字形）。
fn refresh_icon(p: &egui::Painter, rect: Rect, col: Color32) {
    use std::f32::consts::TAU;
    let c = rect.center();
    let r = rect.width() * 0.30;
    let pts: Vec<egui::Pos2> = (0..=24)
        .map(|i| {
            let t = i as f32 / 24.0;
            let a = -0.7 + t * TAU * 0.80;
            egui::pos2(c.x + r * a.cos(), c.y + r * a.sin())
        })
        .collect();
    let stroke = egui::Stroke::new(1.3, col);
    // 末端小箭头
    if let (Some(&a), Some(&b)) = (pts.get(pts.len() - 2), pts.last()) {
        let dir = (b - a).normalized();
        let n = egui::vec2(-dir.y, dir.x);
        p.line_segment([b, b - dir * 3.4 + n * 2.2], stroke);
        p.line_segment([b, b - dir * 3.4 - n * 2.2], stroke);
    }
    p.add(egui::Shape::line(pts, stroke));
}

/// 侧栏路径框回车：跳到该文件夹（SnapGene 的 FOLDER 路径框行为）。
fn navigate_to_input(app: &mut HelixApp) {
    let t = app.folder_input.trim().to_string();
    if t.is_empty() {
        return;
    }
    let pb = PathBuf::from(&t);
    if pb.is_dir() {
        app.set_folder(Some(pb));
        app.refresh_folder();
    } else {
        app.folder_error = Some(format!("not a folder: {t}"));
    }
}

// ---------------------------------------------------------------------------
// 文件树（扁平化 + 就地展开）
// ---------------------------------------------------------------------------

/// 文件树的一行（文件夹或文件），`depth` 从根的子项 0 起算。
#[derive(Debug, Clone, PartialEq)]
pub struct TreeRow {
    pub path: PathBuf,
    /// 该行显示的文本。
    pub label: String,
    pub depth: usize,
    /// 文件夹行（可展开）。
    pub is_dir: bool,
    /// 文件夹行：当前是否展开。
    pub open: bool,
    /// 文件行：是否是我们能打开的序列文件。
    pub is_seq: bool,
    /// 文件行：是否能以只读阅读视图打开（`.txt` / `.pdf`）。
    pub is_view: bool,
    /// 文件行：**原生**分子文档类型（`.dna`/`.rna`/`.prot`，画类型色圆环 + 双螺旋）。
    pub kind: Option<hx_core::MoleculeKind>,
    /// 文件行：序列**文本**格式（`.fa`/`.gbk`/…，画亮蓝圆环 + 折角页）。
    pub text_seq: bool,
}

/// 树的最大展开深度（防御性上限，避免符号链接环）。
const MAX_TREE_DEPTH: usize = 12;

/// 递归展开目录树 —— **就地**：展开的子文件夹内容紧跟它自己。
///
/// `list` 是取目录列表的闭包（真实使用时是带缓存的 `HelixApp::dir_listing`，
/// 单元测试里是内存里的假目录）。`filter` 非空时：只保留名字命中的文件，以及
/// 子树里有命中的文件夹（这些文件夹自动算作展开，不改动 `expanded`）。
pub fn flatten_tree(
    list: &mut dyn FnMut(&std::path::Path) -> Vec<crate::settings::Entry>,
    expanded: &std::collections::HashSet<PathBuf>,
    dir: &std::path::Path,
    depth: usize,
    filter: &str,
    out: &mut Vec<TreeRow>,
) {
    if depth > MAX_TREE_DEPTH {
        return;
    }
    for e in list(dir) {
        // `list_dir` 截断时的「… (只显示前 N 项)」占位行（path 为空）
        if e.path.as_os_str().is_empty() {
            continue;
        }
        let hit = filter.is_empty() || e.name.to_lowercase().contains(filter);
        if e.is_dir {
            if filter.is_empty() {
                let open = expanded.contains(&e.path);
                out.push(TreeRow {
                    path: e.path.clone(),
                    label: display_name(e.name.clone()),
                    depth,
                    is_dir: true,
                    open,
                    is_seq: false,
                    is_view: false,
                    kind: None,
                    text_seq: false,
                });
                if open {
                    flatten_tree(list, expanded, &e.path, depth + 1, filter, out);
                }
            } else if hit {
                // 文件夹名自己命中：连它的内容一起显示（子文件夹仍按展开状态）
                out.push(TreeRow {
                    path: e.path.clone(),
                    label: display_name(e.name.clone()),
                    depth,
                    is_dir: true,
                    open: true,
                    is_seq: false,
                    is_view: false,
                    kind: None,
                    text_seq: false,
                });
                flatten_tree(list, expanded, &e.path, depth + 1, "", out);
            } else {
                // 只看子树里有没有命中项：先递归收集，有命中才把这个文件夹插回它们前面
                let mark = out.len();
                flatten_tree(list, expanded, &e.path, depth + 1, filter, out);
                if out.len() == mark {
                    continue;
                }
                out.insert(
                    mark,
                    TreeRow {
                        path: e.path.clone(),
                        label: display_name(e.name.clone()),
                        depth,
                        is_dir: true,
                        open: true,
                        is_seq: false,
                        is_view: false,
                        kind: None,
                        text_seq: false,
                    },
                );
            }
        } else if hit {
            out.push(TreeRow {
                path: e.path.clone(),
                label: display_name(e.name.clone()),
                depth,
                is_dir: false,
                open: false,
                is_seq: e.is_seq,
                is_view: e.is_view,
                kind: if crate::settings::is_text_seq(&e.ext) {
                    None
                } else {
                    crate::settings::kind_of_ext(&e.ext)
                },
                text_seq: crate::settings::is_text_seq(&e.ext),
            });
        }
    }
}

/// 文件夹图标：SG 的 `📁` 轮廓（带左上角的页签），16×12、`#333` 1.2 px。
fn folder_icon(p: &egui::Painter, rect: Rect, col: Color32) {
    let (x0, x1) = (rect.left(), rect.right());
    let (y0, y1) = (rect.top(), rect.bottom());
    let pts = vec![
        egui::pos2(x0, y0 + 3.0),
        egui::pos2(x0, y0 + 1.2),
        egui::pos2(x0 + 0.9, y0 + 0.6),
        egui::pos2(x0 + 6.0, y0 + 0.6),
        egui::pos2(x0 + 8.0, y0 + 2.8),
        egui::pos2(x1 - 0.9, y0 + 2.8),
        egui::pos2(x1, y0 + 3.8),
        egui::pos2(x1, y1 - 0.8),
        egui::pos2(x1 - 0.9, y1),
        egui::pos2(x0 + 0.9, y1),
        egui::pos2(x0, y1 - 0.8),
        egui::pos2(x0, y0 + 3.0),
    ];
    p.add(egui::Shape::line(pts, egui::Stroke::new(1.2, col)));
}

/// 画一行树节点：`[缩进] [展开三角] [图标列] 名字`（SG folder view 的几何）。
fn paint_tree_row(
    ui: &egui::Ui,
    rect: Rect,
    row: &TreeRow,
    hovered: bool,
    is_current: bool,
) {
    let p = ui.painter();
    let body = Rect::from_min_max(
        egui::Pos2::new(rect.left() + SEL_INSET, rect.top() + 1.0),
        egui::Pos2::new(rect.right() - SEL_INSET, rect.bottom() - 1.0),
    );
    if is_current {
        p.rect_filled(body, 5.0, SEL_BLUE);
    } else if hovered {
        p.rect_filled(body, 5.0, theme::ACCENT_BG);
    }
    let x0 = rect.left() + ROW_PAD + row.depth as f32 * INDENT;
    let cy = rect.center().y;
    // 展开三角（文件夹才有；文件行留空，保证同级名字对齐）
    if row.is_dir {
        let c = egui::Pos2::new(x0 + EXPAND_W * 0.5 - 4.0, cy);
        let pts = if row.open {
            vec![
                c + Vec2::new(-3.2, -1.5),
                c + Vec2::new(3.2, -1.5),
                c + Vec2::new(0.0, 2.4),
            ]
        } else {
            vec![
                c + Vec2::new(-1.5, -3.2),
                c + Vec2::new(2.4, 0.0),
                c + Vec2::new(-1.5, 3.2),
            ]
        };
        let ink = if is_current {
            Color32::WHITE
        } else {
            Color32::from_rgb(0x5a, 0x5a, 0x5a)
        };
        p.add(egui::Shape::convex_polygon(pts, ink, egui::Stroke::NONE));
        // 文件夹图标
        let fr = Rect::from_center_size(
            egui::Pos2::new(x0 + EXPAND_W + ICON_SLOT * 0.5, cy),
            Vec2::new(ICON + 1.0, ICON * 0.80),
        );
        folder_icon(
            p,
            fr,
            if is_current {
                Color32::WHITE
            } else {
                FOLDER_INK
            },
        );
    } else {
        // 原生文档 / 序列文本 / 其它文件
        let cx = x0 + EXPAND_W + ICON_SLOT * 0.5;
        let r = Rect::from_center_size(egui::Pos2::new(cx, cy), Vec2::splat(ICON));
        crate::branding::file_icon(p,r,&row.label,is_current);
    }
    let text_col = if is_current {
        Color32::WHITE
    } else if !row.is_dir && !row.is_seq && !row.is_view {
        theme::TEXT_DIM
    } else {
        theme::TEXT
    };
    // 名字从图标列右边 6 px 起（SG：depth0 文本 x=56 = 14 + 20 + 16 + 6）
    let text_x = x0 + EXPAND_W + ICON_SLOT + 6.0;
    let avail = (rect.right() - SEL_INSET - 4.0 - text_x).max(20.0);
    let galley = ellipsize(
        p,
        &row.label,
        egui::FontId::proportional(ROW_FONT),
        text_col,
        avail,
    );
    p.galley(
        egui::Pos2::new(text_x, cy - galley.size().y * 0.5),
        galley,
        text_col,
    );
}

/// 侧栏右上角的 `...` 溢出菜单（SG 里对应根文件夹标题右侧的 `⋮` 按钮）。
fn overflow_menu(ui: &mut egui::Ui, app: &mut HelixApp) {
    crate::icons::more_menu(ui, |ui| {
        if ui.button(crate::i18n::text(ui.ctx(),"Open Files…    Ctrl+O")).clicked() {
            app.open_files_with_dialog();
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Open Folder…")).clicked() {
            app.open_folder_with_dialog();
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Refresh Folder")).clicked() {
            app.refresh_folder();
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Show Folder in File Manager")).clicked() {
            if let Some(f) = app.folder.clone() {
                app.reveal_path(&f);
            }
            ui.close();
        }
        if app.folder.is_some() && ui.button(crate::i18n::text(ui.ctx(),"Close Folder")).clicked() {
            app.set_folder(None);
            app.folder_input.clear();
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Close Tab    Ctrl+W")).clicked() {
            app.close_doc();
            ui.close();
        }
        ui.separator();
        if ui.button(crate::i18n::text(ui.ctx(),"Save    Ctrl+S")).clicked() {
            app.do_save(false);
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Save As…    Ctrl+Shift+S")).clicked() {
            app.do_save(true);
            ui.close();
        }
        if ui.button(crate::i18n::text(ui.ctx(),"Export View as SVG    Ctrl+E")).clicked() {
            app.do_export_svg();
            ui.close();
        }
    });
}

pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) {
    ui.spacing_mut().item_spacing.y = 0.0;

    // Navigation buttons live in the full-width document-tab row.
    if app.folder.is_none() {
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.add_space(6.0);
        // 手绘「面板」图标：外框 + 左侧实心条（SnapGene 的折叠按钮）
        let (irect, iresp) = ui.allocate_exact_size(Vec2::new(16.0, 16.0), egui::Sense::click());
        {
            let p = ui.painter();
            let r = irect.shrink(1.5);
            p.rect_stroke(
                r,
                2.0,
                egui::Stroke::new(1.1, theme::TEXT_DIM),
                egui::StrokeKind::Inside,
            );
            let mut left = r;
            left.set_right(r.left() + r.width() * 0.34);
            p.rect_filled(left.shrink(0.5), 1.5, theme::TEXT_DIM);
        }
        if iresp.on_hover_text(crate::i18n::text(ui.ctx(),"Hide Files panel")).clicked() {
            app.sidebar_open = false;
        }
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(crate::chrome::NAV_RIGHT_INSET);
            // 打开文件夹后，`...` 菜单跟着**根文件夹标题**走（SG 的 `⋮` 就在标题右侧）；
            // 空状态没有标题行，菜单留在头行。
            if app.folder.is_none() {
                overflow_menu(ui, app);
            }
            ui.add_space(2.0);
            let (r, resp) = ui.allocate_exact_size(Vec2::new(20.0, 20.0), egui::Sense::click());
            let col = if resp.hovered() {
                theme::TEXT
            } else {
                theme::TEXT_DIM
            };
            refresh_icon(ui.painter(), r, col);
            if resp.clicked() {
                app.refresh_folder();
            }
            resp.on_hover_text(crate::i18n::text(ui.ctx(),"Refresh folder"));
        });
    });
    ui.add_space(4.0);

    }
    // ---- 根文件夹标题行 + 过滤框（SG folder view：标题在上、搜索框在下）----
    if let Some(root) = app.folder.clone() {
        let root_name = root
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| root.display().to_string());
        ui.horizontal(|ui| {
            ui.add_space(SEL_INSET);
            ui.label(
                egui::RichText::new(root_name)
                    .size(TITLE_FONT)
                    .color(theme::TEXT),
            )
            .on_hover_text(root.display().to_string());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add_space(crate::chrome::NAV_RIGHT_INSET);
                overflow_menu(ui, app);
            });
        });
        ui.add_space(6.0);
        search_box(ui, &mut app.tree_filter, "Search and filter files");
        ui.add_space(2.0);
    }

    // ---- 空状态：`Open File…` / FOLDER 段 / 路径框 / 拖放提示 ----
    // SnapGene 实测（截图 26/27）里这一段在打开文件夹后仍然在，但用户裁定
    // （2026-10-18）只要打开过文件夹就不再出现 —— 见本文件头部说明。
    if app.folder.is_none() {
        ui.horizontal(|ui| {
            ui.add_space(8.0);
            let btn = egui::Button::new(
                egui::RichText::new(crate::i18n::text(ui.ctx(),"Open File…"))
                    .color(Color32::WHITE)
                    .size(12.5),
            )
            .fill(theme::ACCENT)
            .min_size(Vec2::new(136.0, 24.0));
            if ui.add(btn).clicked() {
                app.open_files_with_dialog();
            }
        });
        ui.add_space(8.0);

        section_label(ui, "FOLDER");
        ui.add_space(3.0);
        ui.horizontal(|ui| {
            ui.add_space(8.0);
            let b = egui::Button::new(egui::RichText::new(crate::i18n::text(ui.ctx(),"Open Folder…")).size(12.0))
                .min_size(Vec2::new(120.0, 22.0));
            if ui.add(b).clicked() {
                app.open_folder_with_dialog();
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add_space(8.0);
                let (r, resp) = ui.allocate_exact_size(Vec2::new(20.0, 20.0), egui::Sense::click());
                let col = if resp.hovered() {
                    theme::TEXT
                } else {
                    theme::TEXT_DIM
                };
                refresh_icon(ui.painter(), r, col);
                if resp.clicked() {
                    app.refresh_folder();
                }
                resp.on_hover_text(crate::i18n::text(ui.ctx(),"Refresh folder"));
            });
        });
        ui.add_space(3.0);
        // 路径框（SnapGene 的 FOLDER 路径输入）：回车跳转
        ui.horizontal(|ui| {
            ui.add_space(8.0);
            let w = ui.available_width() - 8.0;
            let te = ui.add_sized(
                Vec2::new(w, 22.0),
                egui::TextEdit::singleline(&mut app.folder_input).hint_text(crate::i18n::text(ui.ctx(),"/path/to/folder")),
            );
            if te.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                navigate_to_input(app);
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.add_space(10.0);
            ui.label(
                egui::RichText::new(crate::i18n::text(ui.ctx(),"Tip: you can also drop a file onto the window."))
                    .size(10.5)
                    .color(theme::TEXT_DIM),
            );
        });
    }

    // 文件树占上方，`OTHER OPEN FILES` 段固定在下面（SnapGene 的侧栏分区）
    const OTHER_H: f32 = 148.0;
    let tree_h = (ui.available_height() - OTHER_H).max(80.0);

    let current = app
        .doc
        .as_ref()
        .and_then(|d| d.path.as_ref().map(|p| p.display().to_string()))
        .or_else(|| {
            app.viewer
                .as_ref()
                .map(|v| v.path().display().to_string())
        });
    let mut to_open: Option<PathBuf> = None;
    let mut toggle: Option<PathBuf> = None;
    // 右键菜单的动作（在闭包外执行，避免同时借用 app 的其它字段）
    let mut ctx_open: Option<PathBuf> = None;
    let mut ctx_reveal: Option<PathBuf> = None;
    let mut ctx_copy_path: Option<String> = None;

    if let Some(root) = app.folder.clone() {
        // 扁平化（就地展开）→ 渲染
        let filter = app.tree_filter.trim().to_lowercase();
        let expanded = app.expanded.clone();
        let mut rows: Vec<TreeRow> = Vec::new();
        flatten_tree(
            &mut |d| app.dir_listing(d).clone(),
            &expanded,
            &root,
            0,
            &filter,
            &mut rows,
        );

        egui::ScrollArea::vertical()
            .id_salt("folder_tree")
            .auto_shrink([false, true])
            .max_height(tree_h)
            .show(ui, |ui| {
                for row in &rows {
                    let (resp, rect) = list_row(ui);
                    let is_cur = !row.is_dir
                        && current.as_deref() == Some(row.path.display().to_string().as_str());
                    paint_tree_row(ui, rect, row, resp.hovered(), is_cur);
                    if resp.clicked() {
                        if row.is_dir {
                            toggle = Some(row.path.clone());
                        } else if is_cur && resp.interact_pointer_pos().is_some_and(|p|p.x>rect.left()+50.0) {
                            app.rename_file=Some(crate::rename::Rename::new(row.path.clone()));
                        } else if row.is_seq || row.is_view {
                            to_open = Some(row.path.clone());
                        }
                    }
                    if !row.is_dir {
                        // 右键：SnapGene 的文件上下文菜单（Reveal in File Browser / 打开）
                        let path = row.path.clone();
                        let is_openable = row.is_seq || row.is_view;
                        resp.context_menu(|ui| {
                            if ui.button(crate::i18n::text(ui.ctx(),"Rename…    F2")).clicked(){app.rename_file=Some(crate::rename::Rename::new(path.clone()));ui.close();}
                            if ui
                                .add_enabled(is_openable, egui::Button::new(crate::i18n::text(ui.ctx(),"Open")))
                                .clicked()
                            {
                                ctx_open = Some(path.clone());
                                ui.close();
                            }
                            if ui.button(crate::i18n::text(ui.ctx(),"Show in File Manager")).clicked() {
                                ctx_reveal = Some(path.clone());
                                ui.close();
                            }
                            if ui.button(crate::i18n::text(ui.ctx(),"Copy Path")).clicked() {
                                ctx_copy_path = Some(path.display().to_string());
                                ui.close();
                            }
                        });
                    }
                }
            });
    } else {
        // ---- 空状态（SnapGene 截图：No folder opened）----
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.add_space(10.0);
            ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"No folder opened")).strong().size(12.5));
        });
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.add_space(10.0);
            ui.label(
                egui::RichText::new(crate::i18n::text(ui.ctx(),"Project files and subfolders will be shown here."))
                    .size(11.5)
                    .color(theme::TEXT_DIM),
            );
        });
    }

    if let Some(p) = toggle {
        if !app.expanded.insert(p.clone()) {
            app.expanded.remove(&p);
        }
    }
    if let Some(p) = to_open {
        app.open(&p);
    }
    if let Some(p) = ctx_open {
        app.open(&p);
    }
    if let Some(p) = ctx_reveal {
        app.reveal_path(&p);
    }
    if let Some(p) = ctx_copy_path {
        app.pending_copy = Some(p);
        app.status = "Copied path to clipboard".into();
    }
    if let Some(e) = app.folder_error.clone() {
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.add_space(10.0);
            ui.label(egui::RichText::new(e).size(10.5).color(theme::SITE));
        });
    }

    // ---- OTHER OPEN FILES 段 ----
    // SG 的语义：**只列「不在当前文件夹里」的单独打开文件**（`sg-sb.png`：文件夹里
    // 的 pET28b-Rep1.gbk / ECO 16S rRNA.rna 等都在树上，这里只有 pUC19.dna）。
    ui.add_space(8.0);
    {
        let r = ui.available_rect_before_wrap();
        let y = r.top() + 0.5;
        ui.painter().line_segment(
            [
                egui::Pos2::new(r.left() + SEL_INSET, y),
                egui::Pos2::new(r.right() - SEL_INSET, y),
            ],
            egui::Stroke::new(1.0, RULE),
        );
    }
    ui.add_space(8.0);
    section_label(ui, "OTHER OPEN FILES");
    ui.add_space(4.0);
    // SG：**folder view 下这一段没有过滤框**（`sg-sb.png` 只有一个蓝色条目），
    // 只有没打开文件夹的 doc view 才有 `filter open files`（旧证据 26/27）。
    if app.folder.is_none() {
        search_box(ui, &mut app.file_filter, "filter open files");
        ui.add_space(2.0);
    }

    let filter = app.file_filter.trim().to_ascii_lowercase();
    let root = app.folder.clone();
    let open_tabs: Vec<(usize, String, String)> = app
        .tabs
        .iter()
        .enumerate()
        .filter_map(|(i, t)| {
            let d = if i == app.active {
                app.doc.as_ref()
            } else {
                t.doc.as_ref()
            };
            let vw = if i == app.active {
                app.viewer.as_ref()
            } else {
                t.viewer.as_ref()
            };
            let path = d
                .and_then(|d| d.path.as_deref())
                .or_else(|| vw.map(|v| v.path()))?;
            // 只在「不在根文件夹下的文件」里列（在根文件夹内的都已经显示在树上）
            if !is_outside_root(path, root.as_deref()) {
                return None;
            }
            let name = path.file_name().map(|s| s.to_string_lossy().to_string())?;
            Some((i, name, path.display().to_string()))
        })
        .collect();
    let mut switch: Option<usize> = None;
    let list_h = (ui.available_height() - 6.0).max(40.0);
    egui::ScrollArea::vertical()
        .id_salt("other_open")
        .auto_shrink([false, false])
        .max_height(list_h)
        .show(ui, |ui| {
            for (i, name, full) in &open_tabs {
                if !filter.is_empty() && !name.to_ascii_lowercase().contains(&filter) {
                    continue;
                }
                let is_cur = current.as_deref() == Some(full.as_str());
                let (resp, rect) = list_row(ui);
                let ext = std::path::Path::new(full)
                    .extension()
                    .map(|e| e.to_string_lossy().to_ascii_lowercase())
                    .unwrap_or_default();
                let text_seq = crate::settings::is_text_seq(&ext);
                let kind = if text_seq {
                    None
                } else {
                    crate::settings::kind_of_ext(&ext)
                };
                let is_seq = crate::settings::is_openable(&ext);
                let can_open = is_seq || crate::settings::is_viewable(&ext);
                paint_row(
                    ui,
                    rect,
                    is_cur,
                    resp.hovered(),
                    kind,
                    text_seq,
                    name,
                    !can_open,
                );
                if resp.clicked() {
                    if is_cur{app.rename_file=Some(crate::rename::Rename::new(std::path::PathBuf::from(full)));}else{switch = Some(*i);}
                }
            }
        });
    if let Some(i) = switch {
        app.activate_tab(i);
    }

    // ---- 未保存提示（SnapGene 在页签标 *，这里仅在脏时补一行）----
    if let Some(d) = app.doc.as_ref() {
        if d.is_dirty() {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.add_space(10.0);
                ui.label(
                    egui::RichText::new(format!(
                        "{} edit(s) not saved — Ctrl+S",
                        crate::app::nfmt(d.edit_count())
                    ))
                    .size(10.5)
                    .color(theme::SITE),
                );
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::settings::Entry;

    /// 假目录列表：`dirs` 里的路径是文件夹，`files` 里的路径是文件。
    fn fake(
        root: &std::path::Path,
        dirs: &[&str],
        files: &[&str],
    ) -> impl FnMut(&std::path::Path) -> Vec<Entry> {
        let root = root.to_path_buf();
        let dirs: Vec<PathBuf> = dirs.iter().map(|d| root.join(d)).collect();
        let files: Vec<PathBuf> = files.iter().map(|f| root.join(f)).collect();
        move |dir: &std::path::Path| {
            let mut out: Vec<Entry> = Vec::new();
            for p in &dirs {
                if p.parent() == Some(dir) {
                    out.push(Entry {
                        name: p.file_name().unwrap().to_string_lossy().to_string(),
                        path: p.clone(),
                        is_dir: true,
                        is_seq: false,
                        is_view: false,
                        ext: String::new(),
                    });
                }
            }
            for p in &files {
                if p.parent() == Some(dir) {
                    out.push(Entry {
                        name: p.file_name().unwrap().to_string_lossy().to_string(),
                        path: p.clone(),
                        is_dir: false,
                        is_seq: p.extension().is_some_and(|e| e == "dna"),
                        is_view: p.extension().is_some_and(|e| e == "txt" || e == "pdf"),
                        ext: p
                            .extension()
                            .map(|e| e.to_string_lossy().to_ascii_lowercase())
                            .unwrap_or_default(),
                    });
                }
            }
            out
        }
    }

    /// `OTHER OPEN FILES` 只列根文件夹之外的文件（用户在文件夹里打开的文件不重复列）。
    #[test]
    fn other_open_files_excludes_folder_members() {
        let root = std::path::PathBuf::from("/r/Fw__R1-4_sequences");
        assert!(!is_outside_root(&root.join("2333.dna"), Some(&root)));
        assert!(!is_outside_root(
            &root.join("新建文件夹/RNA files/ECO 16S rRNA.rna"),
            Some(&root)
        ));
        assert!(is_outside_root(
            std::path::Path::new("/r/pUC19.dna"),
            Some(&root)
        ));
        // 没打开文件夹时，所有打开的文档都算「其它」
        assert!(is_outside_root(
            std::path::Path::new("/r/pUC19.dna"),
            None
        ));
    }

    /// 展开的子文件夹内容必须**紧跟它自己**，不能掉到列表末尾（用户反馈的回归）。
    #[test]
    fn tree_expansion_is_in_place() {
        let root = std::path::PathBuf::from("/r");
        let mut list = fake(
            &root,
            &["a", "b", "a/inner"],
            &["zz.dna", "a/inner/x.dna", "b/y.dna"],
        );
        let mut expanded = std::collections::HashSet::new();
        expanded.insert(root.join("a"));
        expanded.insert(root.join("a/inner"));
        let mut rows = Vec::new();
        flatten_tree(&mut list, &expanded, &root, 0, "", &mut rows);
        let got: Vec<(String, usize)> = rows.iter().map(|r| (r.label.clone(), r.depth)).collect();
        // a → inner → x(.dna) 必须连在一起，b 在 a 子树之后
        // （序列文件显示文档名，扩展名被去掉，见 `display_name`）
        assert_eq!(
            got,
            vec![
                ("a".to_string(), 0),
                ("inner".to_string(), 1),
                ("x.dna".to_string(), 2),
                ("b".to_string(), 0),
                ("zz.dna".to_string(), 0),
            ]
        );
        // 折叠的 b 不展开 → 内容不出现
        assert!(!rows.iter().any(|r| r.label == "y"));
    }

    /// 过滤：命中文件才留；子树有命中的文件夹自动展开；无命中的文件夹整枝丢掉。
    #[test]
    fn tree_filter_keeps_hits_and_their_folders() {
        let root = std::path::PathBuf::from("/r");
        let mut list = fake(
            &root,
            &["a", "b", "a/inner"],
            &["zz.dna", "a/inner/x.dna", "b/y.txt"],
        );
        let expanded = std::collections::HashSet::new();
        let mut rows = Vec::new();
        flatten_tree(&mut list, &expanded, &root, 0, "x.", &mut rows);
        let got: Vec<(String, usize, bool)> = rows
            .iter()
            .map(|r| (r.label.clone(), r.depth, r.is_dir))
            .collect();
        assert_eq!(
            got,
            vec![
                ("a".to_string(), 0, true),
                ("inner".to_string(), 1, true),
                ("x.dna".to_string(), 2, false),
            ],
            "只留下命中项与它的祖先文件夹，b 整枝丢掉"
        );
        assert!(rows[0].open && rows[1].open);
    }

    /// 过滤命中**文件夹名**时，把它的内容一并显示（子文件夹仍按展开状态折叠）。
    #[test]
    fn tree_filter_folder_hit_shows_its_contents() {
        let root = std::path::PathBuf::from("/r");
        let mut list = fake(&root, &["vec", "vec/inner"], &["notes.txt", "vec/x.dna"]);
        let expanded = std::collections::HashSet::new();
        let mut rows = Vec::new();
        flatten_tree(&mut list, &expanded, &root, 0, "vec", &mut rows);
        let got: Vec<(String, usize)> = rows.iter().map(|r| (r.label.clone(), r.depth)).collect();
        assert_eq!(
            got,
            vec![
                ("vec".to_string(), 0),
                ("inner".to_string(), 1),
                ("x.dna".to_string(), 1)
            ],
            "文件夹 vec 命中 → 它的内容一并显示；未展开的 inner 不再深入，notes.txt 被过滤掉"
        );
    }
}
