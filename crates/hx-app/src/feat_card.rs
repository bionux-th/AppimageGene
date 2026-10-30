//! 序列视图的**特征工具卡**（悬停特征/位点时浮出的卡片）。
//!
//! 版式按真 SnapGene 逐像素实测（oracle 证据图）：
//!
//! - `/tmp/opencode/sg/tt-t1-card.png`：蛋白序列视图悬停 `Region` 条（Mature chain）→
//!   标题 `Mature chain` + 灰芯片 `Region` + 标题行**最右**一个特征色圆角色块 /
//!   区间行 `p1 .. p470  ·  470 aa` / 虚线分隔 / `NOTE` / `Pyruvate kinase I. /id=…`。
//! - `/tmp/opencode/sg/tt-v1-full.png`：悬停多片段 `Site` 标签 → 标题 + 芯片
//!   `Site · active` / `p3 .. p414` / 虚线 / `25 SEGMENTS` / 三个段行
//!   （浅灰圆底数字 + `p3 .. p4`、`p60`、`p62`）/ `+22 Segments` 灰芯片 +
//!   `ⓘ View all segments in features view` 淡蓝行 / 虚线 / `NOTE` + 正文。
//! - `/tmp/opencode/sg/tt-x1-full.png`：单残基 `Site` → 区间只写 `p76`（不写 `..`、不写长度）。
//! - `/tmp/opencode/sg/tt-y1c.png`：单片段 `Region`（`p4 .. p9 · 6 aa`）。
//!
//! 实测尺寸：卡片内边距 ~15 px、标题 ~12.5 px 粗体、芯片底 `#F2F2F2`、正文近黑、
//! 分隔是**虚线** `#C8C8C8`、段行的圆底 `#F5F5F5` 直径 16 px、段行间距 28 px。
//!
//! 刻意差异（做不了或没实测到的，勿当 bug）：
//! 1. SG 卡片外面那圈 8 px 纯黑圆角框是 **X11 窗口管理器给的 tooltip 窗口边框**
//!    （四边等宽、纯黑、圆角），不是 SG 自己画的；这里用 1 px `#D6D6D6` + 阴影代替。
//! 2. 段号在 SG 是浅灰圆底数字（不是圈码），这里照画圆底 + 数字。
//! 3. `+N Segments` / `View all segments…` 在 SG 里可点，这里只作展示。

use eframe::egui;
use egui::Color32;

/// 卡片底色（oracle 实测 `#FFFFFF`）。
const CARD_BG: Color32 = Color32::WHITE;
/// 卡片边框（窗口边框去掉后，取 SG 内部那条 1 px 浅灰 `#D6D6D6`）。
const CARD_EDGE: Color32 = Color32::from_rgb(0xd6, 0xd6, 0xd6);
/// 类型芯片底（oracle 实测 `#EBEBEB`；`tt-v1.png` 横剖 y=254：x 381..462 全是 235）。
const CHIP_BG: Color32 = Color32::from_rgb(0xeb, 0xeb, 0xeb);
/// `+N Segments` 芯片底（oracle 实测 `#F5F5F5`）。
const SEG_CHIP_BG: Color32 = Color32::from_rgb(0xf5, 0xf5, 0xf5);
/// 段号圆底（oracle 实测 `#F5F5F5`）。
const DOT_BG: Color32 = Color32::from_rgb(0xf5, 0xf5, 0xf5);
/// 正文（oracle 实测近黑）。
const TEXT: Color32 = Color32::from_rgb(0x21, 0x21, 0x21);
/// 次级文字（段号、`NOTE`/`N SEGMENTS` 标题、`+N Segments`）。
const DIM: Color32 = Color32::from_rgb(0x56, 0x56, 0x56);
/// 虚线分隔（oracle 实测 `#C8C8C8`）。
const RULE: Color32 = Color32::from_rgb(0xc8, 0xc8, 0xc8);
/// `View all segments…` 行的底（oracle 实测 `#E6F4FF`）。
const VIEW_BG: Color32 = Color32::from_rgb(0xe6, 0xf4, 0xff);
/// 信息圆点（oracle 实测 SG 蓝）。
const INFO: Color32 = Color32::from_rgb(0x0a, 0x84, 0xff);

/// 卡片宽度上限（SG 的卡片最宽约 415 px）。
const MAX_W: f32 = 400.0;

/// 一张工具卡要显示的内容。
pub struct Card<'a> {
    pub name: &'a str,
    /// 特征类型（`Region`/`Site`/`CDS`…）。
    pub kind: &'a str,
    /// 元素上的 `subType`（多为空；`Site` 才有）。
    pub sub_type: &'a str,
    /// `/note` 纯文本（可空）。
    pub note: &'a str,
    /// 0-based 半开分段（整个分子坐标，已按升序）。
    pub segs: &'a [(usize, usize)],
    /// 坐标前缀（蛋白 `p`、核苷酸空串）。
    pub prefix: &'a str,
    /// 单位（`aa`/`bp`）。
    pub unit: &'a str,
    /// 标题行右侧的色块（`None` = 不画，如 `Site`）。
    pub swatch: Option<Color32>,
    /// 方向箭头（DNA 序列视图才有；`Some(true)` = 反向、`Some(false)` = 正向）。
    /// 蛋白视图恒为 `None`（oracle：蛋白卡片的标题行只有色块）。
    pub direction: Option<bool>,
}

/// 卡片的 egui 外框（白底 + 1 px 灰边 + 柔和投影）。
pub fn frame() -> egui::Frame {
    egui::Frame::new()
        .fill(CARD_BG)
        .stroke(egui::Stroke::new(1.0, CARD_EDGE))
        .corner_radius(4.0)
        .inner_margin(egui::Margin::symmetric(14, 10))
        .shadow(egui::epaint::Shadow {
            offset: [1, 2],
            blur: 8,
            spread: 0,
            color: Color32::from_black_alpha(48),
        })
}

/// 画一张工具卡。
pub fn show(ui: &mut egui::Ui, c: &Card<'_>) {
    ui.set_max_width(MAX_W);
    ui.spacing_mut().item_spacing = egui::vec2(6.0, 3.0);

    // ---- 标题行：名字 + 类型芯片（+ 右端色块）----
    let swatch_w = if c.swatch.is_some() { 16.0 } else { 0.0 };
    let arrow_w = if c.direction.is_some() { 16.0 } else { 0.0 };
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 6.0;
        ui.label(
            egui::RichText::new(c.name)
                .size(12.5)
                .color(TEXT)
                .strong(),
        );
        let chip = if c.sub_type.is_empty() {
            c.kind.to_string()
        } else {
            format!("{} · {}", c.kind, c.sub_type)
        };
        egui::Frame::new()
            .fill(CHIP_BG)
            .corner_radius(3.0)
            .inner_margin(egui::Margin::symmetric(5, 1))
            .show(ui, |ui| {
                ui.label(egui::RichText::new(chip).size(11.0).color(DIM));
            });
        // 标题行右端：**先方向箭头、再色块**（oracle `dn4c.png` 的 `ori` 卡片：
        // `←` 在黄块左边）。两者都右对齐到卡片内缘。
        if c.swatch.is_some() || c.direction.is_some() {
            let gap = (ui.available_width() - swatch_w - arrow_w - 2.0).max(0.0);
            ui.add_space(gap);
        }
        if let Some(reverse) = c.direction {
            let (a, _) = ui.allocate_exact_size(egui::vec2(16.0, 16.0), egui::Sense::hover());
            let dir = if reverse { -1.0 } else { 1.0 };
            let cy = a.center().y;
            let tip = egui::pos2(a.center().x + 5.0 * dir, cy);
            let tail = egui::pos2(a.center().x - 5.0 * dir, cy);
            let stroke = egui::Stroke::new(1.4, TEXT);
            ui.painter().line_segment([tail, tip], stroke);
            for dy in [-3.5f32, 3.5] {
                ui.painter().line_segment(
                    [tip, egui::pos2(tip.x - 5.0 * dir, cy + dy)],
                    stroke,
                );
            }
        }
        if let Some(color) = c.swatch {
            let (r, _) = ui.allocate_exact_size(egui::vec2(16.0, 16.0), egui::Sense::hover());
            ui.painter().rect_filled(r, 3.0, color);
            ui.painter().rect_stroke(
                r,
                3.0,
                egui::Stroke::new(0.8, CARD_EDGE),
                egui::StrokeKind::Inside,
            );
        }
    });

    // ---- 区间行 ----
    // oracle：单片段且多于 1 个残基写 `p1 .. p470  ·  470 aa`；**单残基只写 `p76`**；
    // 多片段只写包络 `p3 .. p414`（长度会误导，SG 不写）。
    if let Some(text) = range_text(c) {
        ui.add_space(1.0);
        ui.label(egui::RichText::new(text).size(12.0).color(TEXT));
    }

    // ---- 段列表（多片段才有）----
    if c.segs.len() > 1 {
        ui.add_space(4.0);
        rule(ui);
        ui.add_space(4.0);
        ui.label(
                egui::RichText::new(crate::ui_format!(ui.ctx(),"{} SEGMENTS", c.segs.len()))
                .size(10.5)
                .color(DIM),
        );
        ui.add_space(2.0);
        for (i, (s, e)) in c.segs.iter().take(3).enumerate() {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 8.0;
                let (r, _) =
                    ui.allocate_exact_size(egui::vec2(16.0, 16.0), egui::Sense::hover());
                ui.painter().circle_filled(r.center(), 8.0, DOT_BG);
                ui.painter().text(
                    r.center(),
                    egui::Align2::CENTER_CENTER,
                    format!("{}", i + 1),
                    egui::FontId::proportional(10.0),
                    DIM,
                );
                ui.label(
                    egui::RichText::new(seg_text(c.prefix, *s, *e))
                        .size(12.0)
                        .color(TEXT),
                );
            });
            ui.add_space(4.0);
        }
        if c.segs.len() > 3 {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 8.0;
                egui::Frame::new()
                    .fill(SEG_CHIP_BG)
                    .corner_radius(3.0)
                    .inner_margin(egui::Margin::symmetric(5, 1))
                    .show(ui, |ui| {
                        ui.label(
                    egui::RichText::new(crate::ui_format!(ui.ctx(),"+{} Segments", c.segs.len() - 3))
                                .size(10.5)
                                .color(DIM),
                        );
                    });
                // `ⓘ View all segments in features view`：SG 里是一个淡蓝行 + 蓝色信息点。
                // 我们没有对应动作（Features 面板就是页签），只画出来、不接交互。
                let text = crate::i18n::text(ui.ctx(),"View all segments in features view");
                let font = egui::FontId::proportional(11.5);
                let tw = ui
                    .painter()
                    .layout_no_wrap(text.to_owned(), font.clone(), DIM)
                    .size()
                    .x;
                let (r, _) = ui.allocate_exact_size(
                    egui::vec2(tw + 30.0, 19.0),
                    egui::Sense::hover(),
                );
                ui.painter().rect_filled(r, 3.0, VIEW_BG);
                ui.painter().circle_filled(
                    egui::pos2(r.left() + 10.0, r.center().y),
                    7.0,
                    INFO,
                );
                ui.painter().text(
                    egui::pos2(r.left() + 10.0, r.center().y),
                    egui::Align2::CENTER_CENTER,
                    "i",
                    egui::FontId::proportional(9.5),
                    Color32::WHITE,
                );
                ui.painter().text(
                    egui::pos2(r.left() + 21.0, r.center().y),
                    egui::Align2::LEFT_CENTER,
                    text,
                    font,
                    TEXT,
                );
            });
        }
        if c.note.is_empty() {
            return;
        }
        ui.add_space(2.0);
    }

    // ---- note ----
    if !c.note.is_empty() {
        ui.add_space(4.0);
        rule(ui);
        ui.add_space(4.0);
        ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"NOTE")).size(10.5).color(DIM));
        ui.add_space(1.0);
        ui.label(egui::RichText::new(c.note).size(12.0).color(TEXT));
    }
}

/// 卡片里的虚线分隔（SG 用 1 px 虚线，不是实线）。
fn rule(ui: &mut egui::Ui) {
    let (r, _) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), 1.0),
        egui::Sense::hover(),
    );
    let mut x = r.left();
    while x < r.right() {
        let x1 = (x + 2.0).min(r.right());
        ui.painter().line_segment(
            [egui::pos2(x, r.center().y), egui::pos2(x1, r.center().y)],
            egui::Stroke::new(1.0, RULE),
        );
        x += 4.0;
    }
}

/// 区间行文案。
fn range_text(c: &Card<'_>) -> Option<String> {
    let (s, e) = *c.segs.first()?;
    if c.segs.len() == 1 {
        return Some(if e - s > 1 {
            format!(
                "{p}{a} .. {p}{b}  ·  {n} {u}",
                p = c.prefix,
                a = s + 1,
                b = e,
                n = e - s,
                u = c.unit
            )
        } else {
            // 单个残基：只有起点（oracle `tt-x1-full.png` 的 `p76`）。
            format!("{}{}", c.prefix, s + 1)
        });
    }
    let lo = c.segs.iter().map(|s| s.0).min().unwrap_or(s);
    let hi = c.segs.iter().map(|s| s.1).max().unwrap_or(e);
    Some(format!("{}{} .. {}{}", c.prefix, lo + 1, c.prefix, hi))
}

/// 段行文案：单残基段只写起点（oracle `2 p60`）。
fn seg_text(prefix: &str, s: usize, e: usize) -> String {
    if e - s > 1 {
        format!("{p}{} .. {p}{}", s + 1, e, p = prefix)
    } else {
        format!("{}{}", prefix, s + 1)
    }
}
