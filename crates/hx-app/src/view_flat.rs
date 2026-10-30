//! 单链序列视图：**RNA / 蛋白 / 单链 DNA**（双链 DNA 走 `view_seq`）。
//!
//! 版式按真 SnapGene **1:1 像素实测**重写（蛋白 `/tmp/opencode/sg/or2-prot-seq2.png`、
//! RNA `/tmp/opencode/sg/i0-rna-seq.png`；量法见 `docs/protein-rna-seq-parity-*.zh-CN.md`）：
//!
//! ```text
//!      ┊                    ┊                            ← 位点轴（蛋白 `type="Site"`；多车道）
//!   ───┴────────────────────┴─────────────────────       ← 车道线 + 每片段小刻度条 + 立杆
//!   M K T K I V C T I G P K T E S E E M L A K M L        ← 字母：等宽、**连续无组间隙**
//!   ─┴────┴────┴────┴────┴────┴───                       ← 标尺线（刻度**朝上**：每 5 短、每 10 长）
//!   1    10    20    30    40    50    60                 ← 刻度数字在**线的下方**
//!   [═══ Mature chain ═══════════════════════]            ← 特征条 17 px，放得下时名字写在条内
//!   [Beta-strand]  [Helical region]                       ← 放不下时排到条下行（细引线连回条上）
//! ```
//!
//! 与 DNA 双链视图（`view_seq`）的区别：只有一条链，没有互补链、酶切轨道、翻译轨道与 ORF。
//! 大小写按文件原样（SnapGene 保留小写）；**行末不再单独写坐标**——标尺自己的数字
//! （10/20/…）就含行末那一个；分子第一行的第一个残基额外写 `1`。
//! 蛋白刻度写纯数字（SnapGene 序列视图没有 `p` 前缀，`p` 只出现在状态栏等处）。
//!
//! oracle 实测数值（z = 1）：字母 cap 10–11 px、字距 11.0；标尺线在字母带上沿 +16（蛋白）
//! / +21（RNA，字母带含降部）；首条车道在字母带上沿 +38；条高 17、车道间距 29；
//! 条下标签行距 16；位点轴车道间距 14、最下轴线在字母上方 11、小刻度条在轴线下 7。
//! 背景 `#F2F2F2`、字母近黑、标尺线灰、条描边中性深灰 `#434343`、条下标签近黑。

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke};

use crate::app::{feat_group, nfmt_plain, FeatGroup, HelixApp};
use crate::{feat_card, theme};

// ---- 版式常量（z = 1 时与 SnapGene 1:1 实测一致）----
/// 每个残基的水平步距（蛋白/RNA 实测都是 11 px）。
const BASE_W: f32 = 11.0;
/// 字母字号（Hack cap/em ≈ 0.727 ⇒ cap ≈ 10.5 px，对齐实测的 10–11 px）。
const LETTER_FONT: f32 = 14.5;
/// 标尺刻度数字字号（实测数字高 8–9 px ⇒ cap/em 0.7 时约 13 pt）。
const RULER_FONT: f32 = 13.0;
/// 特征条内/条下标签字号。
const SMALL_FONT: f32 = 10.5;
/// 位点标签字号。
const SITE_FONT: f32 = 11.5;
/// 标尺线相对字母带上沿的下移（蛋白；RNA 的小写字母有降部，线要更低）。
const RULER_DY_PROT: f32 = 16.0;
const RULER_DY_NUCL: f32 = 21.0;
/// 刻度长度（朝上）：每 5 短、每 10 长。
const TICK_MIN_H: f32 = 4.0;
const TICK_MAJ_H: f32 = 9.0;
/// 标尺数字行高。
const RULER_NUM_H: f32 = 12.0;
/// 字母带上沿 → 第一条车道（蛋白/RNA 实测都是 38）。
const BARS_TOP: f32 = 38.0;
/// 特征条高度（实测 17）。
const BAR_H: f32 = 17.0;
/// 车道间距（实测 29 = 17 + 12）。
const LANE_DY: f32 = 29.0;
/// 条下标签行距（实测 16）。
const LABEL_DY: f32 = 16.0;
/// 条底 → 第一条标签行。
const LABEL_GAP: f32 = 6.0;
/// 位点轴车道间距（实测 14）。
const SITE_LANE_DY: f32 = 14.0;
/// 位点小刻度条在其车道线下方的距离（实测 7）。
const SITE_MARK_DY: f32 = 7.0;
/// 位点小刻度条厚度。
const SITE_MARK_H: f32 = 2.0;
/// 最下一条轴线与字母带上方的距离（实测 11）。
const AXIS_GAP: f32 = 11.0;
/// 轴线上方标签行：底边距最上方轴线 13，行高 15。
const SITE_LABEL_UP: f32 = 13.0;
const SITE_LABEL_H: f32 = 15.0;
/// 行间空隙（上一行内容底 → 下一行的轴区/字母带）。
const ROW_GAP_PROT_SITE: f32 = 42.0;
const ROW_GAP: f32 = 34.0;
/// 视图顶部留白：视图页签与右上角控件是**浮在画布上**的
/// （`app.rs::view_tabs_floating`），所以第一行的轴区/字母带要整体下移让开它们。
const TOP_PAD: f32 = 44.0;
/// 底部概览条高度。
const BOTTOM_H: f32 = 30.0;
/// 序列左侧留白（或实测：蛋白与 RNA 的字母都从视图左缘 +130 px 起，给
/// 悬浮竖排工具条和 `5'` 留位）。
const PAD_LEFT_PROT: f32 = 128.0;
const PAD_LEFT_NUCL: f32 = 128.0;
/// 序列区右侧留白：oracle 实测 SnapGene 的序列视图**换行到窗口**时右端留 150 px
/// （1230 px 窗口里每行 60 个单位 = 内容宽 940 − 左 130 − 右 150），
/// 这段空位放行末的绝对坐标。行宽的自适应按它算（用户 2026-10-18 指摘 2K 桌面下太窄）。
const PAD_RIGHT: f32 = 150.0;
/// `5'` 标记右边缘到第一个字母的距离（实测 ≈45）。
const FIVE_PRIME_DX: f32 = 45.0;
/// 条下标签最多排几行（排不下就重叠，不再加高行）。
const MAX_LABEL_ROWS: usize = 5;
/// 字母的中心线相对字母带上沿的下移（cap 10.5 px ⇒ 中心在 cap 带中间）。
const LETTER_CENTER_DY: f32 = 5.5;

/// 一条特征在单链视图里的排版输入。
struct FlatFeat {
    /// 文档内特征下标（用于选中态）。
    idx: usize,
    name: String,
    /// 特征类型原文（工具卡芯片用）。
    kind: String,
    /// 元素上的 `subType`（`Site` 才有）。
    sub_type: String,
    /// `/note` 纯文本（工具卡的 `NOTE` 段）。
    note: String,
    /// `type="Site"`：蛋白里 SnapGene 画到字母**上方**的细轴（不画条）。
    site: bool,
    color: Color32,
    /// 是否正向（蛋白恒为 true：无方向）。
    forward: bool,
    /// 回绕后的分段（0-based 半开，整个分子坐标）。
    segs: Vec<(usize, usize)>,
}

/// 本行里一条特征的排版结果。
struct RowFeat {
    /// 行内 0-based 半开分段（空 = 本行不画）。
    spans: Vec<(usize, usize)>,
    /// 车道号（Site 与普通条各自编号）；`usize::MAX` = 本行不画。
    lane: usize,
}

/// 条下标签的排布（x 是**行内**像素，0 = `content_left`）。
struct LabelPlan {
    row: usize,
    x: f32,
    /// 引线所在的 x（条的中心）。
    leader_x: f32,
    /// 引线起点 y（条的底边，行内像素）。
    leader_y: f32,
    text: String,
    /// 特征下标（点标签 = 选中该特征）。
    idx: usize,
}

/// 位点轴上的标签（轴区最上面一行）。
struct SiteLabel {
    /// 中心 x（行内像素）。
    cx: f32,
    text: String,
    idx: usize,
}

/// 一行的完整排版。
struct RowPlan {
    feats: Vec<RowFeat>,
    site_lanes: usize,
    bar_lanes: usize,
    labels: Vec<LabelPlan>,
    /// Site 标签（中心 x、文字、特征下标），全部排在轴区最上面一行。
    site_labels: Vec<SiteLabel>,
    /// 字母带**上方**的高度（位点轴 + 轴标签）。
    axis_h: f32,
    /// 字母带上沿到本行内容底的高度。
    body_h: f32,
    /// 整行高度（= axis_h + body_h + 行间空隙）。
    height: f32,
}

pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) -> egui::Rect {
    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    if rect.width() < 160.0 || rect.height() < 120.0 {
        return rect;
    }
    // 序列页整块底色 `#F2F2F2`（图谱页才是纯白）
    painter.rect_filled(rect, 0.0, theme::SEQ_BG);
    let Some(doc) = app.doc.as_ref() else {
        return rect;
    };
    let kind = doc.kind();
    let len = doc.len();
    if len == 0 {
        return rect;
    }
    let seq: Vec<u8> = doc.sequence().to_vec();
    let directional = kind.features_have_direction();
    let nucleotide = kind.is_nucleotide();

    // ---- 特征（隐藏的跳过；蛋白无方向 ⇒ 无箭头）----
    let mut feats: Vec<FlatFeat> = Vec::new();
    for (fi, f) in doc.features().iter().enumerate() {
        // 类别开关（核苷酸 = `Features`；蛋白 = `Region` / `Site`+`Het` / `Bond`）
        // 与文件里 `visible="0"` 的隐藏特征都在这一个判据里。
        if !app.feature_visible(fi, f) {
            continue;
        }
        let color = theme::hex(
            f.segments.first().map(|s| s.color.as_str()).unwrap_or(""),
            theme::FEATURE_FALLBACK,
        );
        let forward = !directional || !f.has_reverse_arrow();
        let segs = f
            .segments
            .iter()
            .flat_map(|s| hx_render::split_wrap(s.range.start, s.range.end, len))
            .collect();
        // 蛋白的 `Site` / `Het` / `Bond`：SnapGene 画在字母上方的细刻度轴（oracle：
        // KPYK1_ECOLI 的 `Site` / `domain interfaces [active]`），不画条。
        // 核苷酸没有这种画法（核苷酸的 Site 是普通特征条）。
        let site = !nucleotide && !matches!(feat_group(&f.kind), FeatGroup::Region);
        feats.push(FlatFeat {
            idx: fi,
            name: f.name.clone(),
            kind: f.kind.clone(),
            sub_type: f.sub_type.clone(),
            note: f.note.clone(),
            site,
            color,
            forward,
            segs,
        });
    }

    // ---- 版式 ----
    let z = app.seq_zoom.clamp(0.45, 2.6);
    let base_w = BASE_W * z;
    let ruler_dy = if nucleotide {
        RULER_DY_NUCL
    } else {
        RULER_DY_PROT
    } * z;
    let pad_left = if nucleotide {
        PAD_LEFT_NUCL
    } else {
        PAD_LEFT_PROT
    } * z;
    let content_left = rect.left() + pad_left;
    // 行宽：`None` = 换行到窗口（自适应，取 10 的整数倍）；`Some(n)` = 固定每行 n 个单位。
    // oracle：SnapGene 的序列视图就是按窗口宽度换行（1230 px 窗口 → 每行 60 个单位），
    // 2K 桌面上必须铺满视图区（用户 2026-10-18 指摘我们固定 60 太窄）。
    let avail = (rect.right() - PAD_RIGHT * z - 8.0 - content_left).max(40.0);
    let row_units = fit_row_units(avail, base_w, app.seq_line_bp);
    let row_px = row_units as f32 * base_w;
    let content_right = (content_left + row_px).min(rect.right() - 8.0);
    let content_w = (content_right - content_left).max(40.0);

    let mono = FontId::monospace(LETTER_FONT * z);
    let ruler_font = FontId::proportional(RULER_FONT * z);
    let small = FontId::proportional(SMALL_FONT * z);
    let site_font = FontId::proportional(SITE_FONT * z);
    let measure = {
        let p = painter.clone();
        move |text: &str, size: f32| -> f32 {
            p.layout_no_wrap(text.to_owned(), FontId::proportional(size), theme::TEXT)
                .size()
                .x
        }
    };
    // 文字标签的命中区（点标签 = 选中该特征；oracle 图五：SG 里点 `Site`/特征名
    // 会直接高亮对应残基）。`x` 是**行内**像素，`anchor_y` 是屏幕 y。
    let text_hit = |p: Pos2, anchor_y: f32, x: f32, text: &str, size: f32, align: Align2| -> bool {
        let tw = measure(text, size);
        let x0 = content_left + x;
        let (x0, x1) = if align == Align2::CENTER_BOTTOM {
            (x0 - tw * 0.5, x0 + tw * 0.5)
        } else {
            (x0, x0 + tw)
        };
        let (y0, y1) = if align == Align2::CENTER_BOTTOM {
            (anchor_y - size * 1.05, anchor_y + size * 0.3)
        } else {
            (anchor_y - size * 0.25, anchor_y + size * 1.0)
        };
        p.x >= x0 - 3.0 && p.x <= x1 + 3.0 && p.y >= y0 && p.y <= y1
    };

    // 顶部让开两条浮动工具条（视图页签 + 右上角控件），与 DNA 序列视图同一套留白；
    // 核苷酸还要在这条留白里塞下顶部固定标尺，所以行区整体再下移 `62 + 22z`。
    let top_inset = if nucleotide {
        62.0 + 22.0 * z
    } else {
        TOP_PAD + 8.0 * z
    };
    let content_top = rect.top() + top_inset;
    let content_bottom = rect.bottom() - BOTTOM_H - 6.0;
    let view_h = (content_bottom - content_top).max(40.0);

    let n_rows = len.div_ceil(row_units);
    let mut plans: Vec<RowPlan> = Vec::with_capacity(n_rows);
    for row in 0..n_rows {
        let start = row * row_units;
        let end = (start + row_units).min(len);
        plans.push(plan_row(&feats, start, end, base_w, z, content_w, &measure));
    }
    let total_h: f32 = plans.iter().map(|p| p.height).sum();
    // 每行字母带上沿（行内坐标）
    let mut row_letters: Vec<f32> = Vec::with_capacity(n_rows);
    {
        let mut y = 0.0f32;
        for p in &plans {
            y += p.axis_h;
            row_letters.push(y);
            y += p.body_h + (p.height - p.axis_h - p.body_h);
        }
    }

    // ---- 滚动 / 缩放 ----
    let max_scroll = (total_h + TOP_PAD - view_h).max(0.0);
    if resp.hovered() {
        let (scroll, ctrl) = ui.input(|i| {
            (
                i.smooth_scroll_delta.y,
                i.modifiers.ctrl || i.modifiers.command,
            )
        });
        if ctrl {
            app.seq_zoom = (app.seq_zoom * (1.0 + scroll * 0.0015)).clamp(0.45, 2.6);
        } else if scroll.abs() > 0.1 {
            app.seq.scroll_y = (app.seq.scroll_y - scroll).clamp(0.0, max_scroll);
        }
    }
    app.seq.scroll_y = app.seq.scroll_y.clamp(0.0, max_scroll);
    let scroll = app.seq.scroll_y;

    // 行内 x（残基下标 → 像素）
    let gx = |unit: usize| -> f32 { content_left + unit as f32 * base_w };
    let center_x = |unit: usize| -> f32 { content_left + (unit as f32 + 0.5) * base_w };

    // 命中测试：屏幕坐标 → 残基下标（只有字母带里才定位）
    let unit_at = |p: Pos2| -> Option<usize> {
        let cy = p.y - content_top + scroll;
        if cy < 0.0 || cy > total_h {
            return None;
        }
        let mut row = 0usize;
        while row + 1 < n_rows && cy > row_letters[row] + plans[row].body_h {
            row += 1;
        }
        let ly = content_top + row_letters[row] - scroll;
        if p.y < ly - 6.0 || p.y > ly + ruler_dy + RULER_NUM_H * z {
            return None;
        }
        if p.x < content_left - 4.0 || p.x > content_right + 4.0 {
            return None;
        }
        let col = ((p.x - content_left) / base_w).floor().max(0.0) as usize;
        let bp = row * row_units + col.min(row_units.saturating_sub(1));
        (bp < len).then_some(bp)
    };

    // ---- 交互：点击定位、拖拽选段、点击特征条/标签选中特征 ----
    // 命中测试与绘制用**同一坐标系**：`ly` 是字母带上沿的屏幕 y（绘制里的 `ly`）。
    let feature_at = |p: Pos2| -> Option<usize> {
        for row in 0..n_rows {
            let plan = &plans[row];
            let ly = content_top + row_letters[row] - scroll;
            // 行内容的下沿：条下标签在**绘制**时从 `bars_y0 + bar_lanes*LANE_DY` 起排
            // （比 `body_h` 里记的 `bars_h` 多一个车道间距），命中测试要覆盖到那里。
            let labels_h = plan
                .labels
                .iter()
                .map(|l| (l.row + 1) as f32 * LABEL_DY)
                .fold(0.0f32, f32::max)
                * z;
            let row_bottom = ly
                + BARS_TOP * z
                + plan.bar_lanes.max(1) as f32 * LANE_DY * z
                + LABEL_GAP * z
                + labels_h
                + SMALL_FONT * z;
            if p.y < ly - plan.axis_h || p.y > row_bottom {
                continue;
            }
            let start = row * row_units;
            let end = (start + row_units).min(len);
            // 普通条
            let bars_y0 = ly + BARS_TOP * z;
            for (fi, f) in feats.iter().enumerate() {
                let lane = plan.feats[fi].lane;
                if f.site || lane == usize::MAX {
                    continue;
                }
                let by = bars_y0 + lane as f32 * LANE_DY * z;
                if p.y < by - 2.0 || p.y > by + BAR_H * z + 2.0 {
                    continue;
                }
                for (a, b) in visible_spans(f, start, end) {
                    let x0 = gx(a);
                    let x1 = gx(b);
                    if p.x >= x0 - 1.0 && p.x <= x1 + 1.0 {
                        return Some(f.idx);
                    }
                }
            }
            // 位点轴上的小刻度条
            for (fi, f) in feats.iter().enumerate() {
                let lane = plan.feats[fi].lane;
                if !f.site || lane == usize::MAX {
                    continue;
                }
                let line_y = ly - AXIS_GAP * z - lane as f32 * SITE_LANE_DY * z;
                let my = line_y + SITE_MARK_DY * z;
                if p.y < line_y - 3.0 || p.y > my + SITE_MARK_H * z + 3.0 {
                    continue;
                }
                for (a, b) in visible_spans(f, start, end) {
                    if p.x >= gx(a) - 2.0 && p.x <= gx(b) + 2.0 {
                        return Some(f.idx);
                    }
                }
            }
            // 位点标签（轴区最上一行）
            if plan.site_lanes > 0 {
                let label_y = ly
                    - AXIS_GAP * z
                    - (plan.site_lanes - 1) as f32 * SITE_LANE_DY * z
                    - SITE_LABEL_UP * z;
                for lab in &plan.site_labels {
                    if text_hit(
                        p,
                        label_y,
                        lab.cx,
                        &lab.text,
                        SITE_FONT * z,
                        Align2::CENTER_BOTTOM,
                    ) {
                        return Some(lab.idx);
                    }
                }
            }
            // 条下标签（细引线下方那几行字）
            let bars_bottom = bars_y0 + plan.bar_lanes.max(1) as f32 * LANE_DY * z;
            for lab in &plan.labels {
                let ly2 = bars_bottom + LABEL_GAP * z + lab.row as f32 * LABEL_DY * z;
                if text_hit(p, ly2, lab.x, &lab.text, SMALL_FONT * z, Align2::LEFT_TOP) {
                    return Some(lab.idx);
                }
            }
        }
        None
    };

    if resp.clicked() {
        // 单击 = 放弃文本控件焦点（否则侧栏搜索框持有键盘，输入进不了序列）
        ui.ctx().memory_mut(|m| m.stop_text_input());
        if let Some(p) = resp.interact_pointer_pos() {
            if let Some(fi) = feature_at(p) {
                app.selected_feature = Some(fi);
                app.cursor_bp = None;
                app.selection = None;
                app.sel_anchor = None;
            } else {
                app.cursor_bp = unit_at(p);
                app.selection = None;
                app.selected_feature = None;
                app.sel_anchor = None;
            }
        }
    }
    if resp.drag_started() {
        app.sel_anchor = resp.interact_pointer_pos().and_then(unit_at);
        app.selected_feature = None;
    }
    if resp.dragged() {
        if let (Some(a), Some(p)) = (app.sel_anchor, resp.interact_pointer_pos()) {
            if let Some(b) = unit_at(p) {
                let (s, e) = if a <= b { (a, b + 1) } else { (b, a + 1) };
                app.selection = Some((s, e));
                app.cursor_bp = Some(b);
            }
        }
    }
    if resp.drag_stopped() {
        app.sel_anchor = None;
    }

    // ---- 绘制 ----
    // 行区裁剪的顶边：核苷酸的行**不能画到顶部固定标尺上面**去（oracle：SnapGene 的
    // 序列永远不会出现在标尺上方；用户 2026-10-18 指摘 图三 的滚动残影）。
    // 蛋白视图没有顶部标尺，行一直可以画到画布上沿（悬浮工具条是浮在上面的）。
    let rows_top = if nucleotide && app.show_ruler {
        content_top - 6.0 * z + 1.0
    } else {
        rect.top() + 2.0
    };
    // 原始 painter 留给顶部标尺：`with_clip_rect` 是**求交**，在行区裁剪的 painter 上
    // 再裁标尺带会得到空区（标尺画不出来）。
    let painter_all = painter.clone();
    let clip = Rect::from_min_max(
        Pos2::new(rect.left() + 1.0, rows_top),
        Pos2::new(rect.right(), content_bottom),
    );
    let painter = painter.with_clip_rect(clip);
    let sel = app.selection;
    let cur = app.cursor_bp;
    // 悬停的特征（点位点括线/标签 → 该特征）：oracle 实测 SnapGene 悬停位点括线时
    // 括线变红、标签底色变浅蓝并弹工具卡；这里取悬停特征用于这两处着色。
    let hover_feat = resp.hover_pos().and_then(feature_at);

    // ---- 顶部固定标尺（**只有核苷酸**）----
    // oracle：SnapGene 的 DNA/RNA 序列视图在字母上方都有一条蓝色刻度标尺
    // （刻度每 5 短、每 10 长，数字 10/20/… 右对齐在刻度上，右端一个小方块抓手）；
    // **蛋白视图没有**这条标尺，蛋白的坐标只写在每行字母下方。
    // 几何与 `view_seq.rs` 的 DNA 序列视图一致（`ruler_h = 62 + 22z`）。
    if nucleotide && app.show_ruler {
        // 标尺自己的一块裁剪区（行区裁剪会把它切掉；用未裁剪的 painter 求交）
        let painter = painter_all.with_clip_rect(Rect::from_min_max(
            Pos2::new(rect.left() + 1.0, rect.top() + 2.0),
            Pos2::new(rect.right(), rows_top),
        ));
        let mut first = 0usize;
        while first + 1 < n_rows && row_letters[first] + plans[first].body_h < scroll {
            first += 1;
        }
        let start = first * row_units;
        let end = (start + row_units).min(len);
        let y = content_top - 6.0 * z;
        let row_right = gx(end - start).min(rect.right() - 8.0);
        painter.rect_filled(
            Rect::from_min_max(
                Pos2::new(rect.left() + 1.0, rect.top() + 2.0),
                Pos2::new(rect.right(), y + 1.0),
            ),
            0.0,
            theme::SEQ_BG,
        );
        painter.line_segment(
            [Pos2::new(content_left, y), Pos2::new(row_right, y)],
            Stroke::new(1.0, theme::ACCENT),
        );
        for i in 0..(end - start) {
            let unit = start + i + 1;
            let major = unit % 10 == 0 || (start == 0 && i == 0);
            if !(major || unit % 5 == 0) {
                continue;
            }
            let x = center_x(i);
            let h = if major { TICK_MAJ_H } else { TICK_MIN_H } * z;
            painter.line_segment(
                [Pos2::new(x, y - h), Pos2::new(x, y + 1.0 * z)],
                Stroke::new(1.0, theme::ACCENT),
            );
            if unit % 10 == 0 {
                painter.text(
                    Pos2::new(x - 2.0 * z, y - 2.0 * z),
                    Align2::RIGHT_BOTTOM,
                    nfmt_plain(unit),
                    ruler_font.clone(),
                    theme::RULER_TEXT,
                );
            }
        }
        // 右端小方块（SnapGene 标尺的抓手）
        painter.rect_filled(
            Rect::from_min_max(
                Pos2::new(row_right - 5.0 * z, y - 9.0 * z),
                Pos2::new(row_right, y + 1.0 * z),
            ),
            1.0,
            theme::ACCENT,
        );
    }

    for row in 0..n_rows {
        let plan = &plans[row];
        let start = row * row_units;
        let end = (start + row_units).min(len);
        let ly = content_top + row_letters[row] - scroll;
        if ly - plan.axis_h > content_bottom + 40.0 || ly + plan.body_h < content_top - 40.0 {
            continue;
        }
        let row_units_here = end - start;
        let row_left = content_left;
        let row_right = gx(row_units_here);

        // ---- 位点轴（字母上方）----
        if plan.site_lanes > 0 {
            let label_y = ly
                - AXIS_GAP * z
                - (plan.site_lanes - 1) as f32 * SITE_LANE_DY * z
                - SITE_LABEL_UP * z;
            for lane in 0..plan.site_lanes {
                let line_y = ly - AXIS_GAP * z - lane as f32 * SITE_LANE_DY * z;
                let mark_y = line_y + SITE_MARK_DY * z;
                for (fi, f) in feats.iter().enumerate() {
                    let p = &plan.feats[fi];
                    if !f.site || p.lane != lane || p.spans.is_empty() {
                        continue;
                    }
                    // 选中的位点：SG 把轴线与刻度都画成选择蓝（oracle `/tmp/opencode/sg/hv-site-sel*.png`）；
                    // 悬停的位点：轴线与刻度变红（oracle 实测 `(253,1,1)`）。
                    let is_sel = app.selected_feature == Some(f.idx);
                    let is_hover = !is_sel && hover_feat == Some(f.idx);
                    let axis_color = if is_sel {
                        theme::FEAT_SEL_EDGE
                    } else if is_hover {
                        theme::SITE_HOVER
                    } else {
                        theme::SITE_AXIS
                    };
                    let axis_w = if is_sel || is_hover { 1.8 * z } else { 1.0 * z };
                    // 车道线：从本行第一个片段中心到本行最后一个片段中心；
                    // 特征在本行之后还有片段 ⇒ 拉到行末（oracle 的写法）。
                    let x0 = center_x(p.spans[0].0);
                    let continues = f.segs.iter().any(|(a, _)| *a >= end);
                    let x1 = if continues {
                        row_right
                    } else {
                        center_x(p.spans[p.spans.len() - 1].1 - 1)
                    };
                    painter.line_segment(
                        [
                            Pos2::new(x0, line_y),
                            Pos2::new(x1.max(x0 + 6.0 * z), line_y),
                        ],
                        Stroke::new(axis_w, axis_color),
                    );
                    // 每个片段：轴线下方的小刻度条 + 连回轴线的立杆
                    for (a, b) in &p.spans {
                        let mx0 = gx(*a);
                        let mx1 = (gx(*b)).max(mx0 + 4.0 * z);
                        let mid = (mx0 + mx1) * 0.5;
                        painter.rect_filled(
                            Rect::from_min_max(
                                Pos2::new(mx0, mark_y),
                                Pos2::new(mx1, mark_y + SITE_MARK_H * z),
                            ),
                            0.0,
                            axis_color,
                        );
                        painter.line_segment(
                            [Pos2::new(mid, line_y), Pos2::new(mid, mark_y)],
                            Stroke::new(axis_w, axis_color),
                        );
                    }
                }
            }
            for lab in &plan.site_labels {
                let x = (content_left + lab.cx).clamp(rect.left() + 4.0, rect.right() - 4.0);
                let selected = app.selected_feature == Some(lab.idx);
                let hovered = !selected && hover_feat == Some(lab.idx);
                // 引线：标签下沿 → 本特征所在的括线（oracle 实测 `#7F7F7F` 细线；
                // 选中时同色变蓝、悬停时变红）。
                let lane = feats
                    .iter()
                    .position(|f| f.idx == lab.idx)
                    .map(|i| plan.feats[i].lane)
                    .filter(|l| *l != usize::MAX);
                if let Some(lane) = lane {
                    let line_y = ly - AXIS_GAP * z - lane as f32 * SITE_LANE_DY * z;
                    let col = if selected {
                        theme::FEAT_SEL_EDGE
                    } else if hovered {
                        theme::SITE_HOVER
                    } else {
                        theme::SITE_LEADER
                    };
                    painter.line_segment(
                        [
                            Pos2::new(x, label_y + 2.0 * z),
                            Pos2::new(x, line_y),
                        ],
                        Stroke::new(0.8, col),
                    );
                }
                let tw = measure(&lab.text, SITE_FONT * z);
                let box_rect = Rect::from_min_max(
                    Pos2::new(x - tw * 0.5 - 5.0 * z, label_y - SITE_FONT * z * 0.92 - 3.0 * z),
                    Pos2::new(x + tw * 0.5 + 5.0 * z, label_y + 3.0 * z),
                );
                if selected {
                    // oracle：选中标签 = 蓝底白字圆角块
                    painter.rect_filled(box_rect, 3.0 * z, theme::FEAT_SEL_EDGE);
                } else if hovered {
                    // oracle `tt-v1-lab.png` / `tt-x1-lab.png`：悬停标签 = **红字 + 红圆角
                    // 描边**，底色不变（此前记的「浅蓝底 + 蓝边」是误读：`#E6F4FF` 那块
                    // 是卡片里 `View all segments…` 按钮的底）。
                    painter.rect_stroke(
                        box_rect,
                        3.0 * z,
                        Stroke::new(1.2, theme::SITE_HOVER),
                        egui::StrokeKind::Inside,
                    );
                }
                painter.text(
                    Pos2::new(x, label_y),
                    Align2::CENTER_BOTTOM,
                    &lab.text,
                    site_font.clone(),
                    if selected {
                        Color32::WHITE
                    } else if hovered {
                        theme::SITE_HOVER
                    } else {
                        theme::TEXT
                    },
                );
            }
        }

        // ---- 选区 / 光标（画在字母带里）----
        let letters_top = ly - 2.0 * z;
        let letters_bottom = ly + ruler_dy;
        if let Some((s0, e0)) = sel {
            for (a, b) in hx_render::split_wrap(s0, e0, len) {
                let a = a.max(start);
                let b = b.min(end);
                if a >= b {
                    continue;
                }
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(gx(a - start), letters_top),
                        Pos2::new(gx(b - start), letters_bottom),
                    ),
                    0.0,
                    theme::with_alpha(theme::SELECTION, 60),
                );
            }
        }
        if let Some(c) = cur {
            if c >= start && c <= end {
                let x = gx(c.min(end) - start);
                painter.line_segment(
                    [
                        Pos2::new(x, letters_top - 2.0 * z),
                        Pos2::new(x, letters_bottom + 2.0 * z),
                    ],
                    Stroke::new(1.4 * z, theme::CURSOR),
                );
            }
        }

        // ---- 每行标尺（字母下方）----
        // oracle 实测两种写法完全不同：
        //   蛋白：线 + **向下**的短刻度（每 5 单位 2 px、行首/行末 11 px），
        //         数字写在**线下方**、右对齐在各十进位刻度前（首行补 `1`）。
        //   核苷酸：线 + **向上**的密刻度（每单位 2 px、每 5 单位 4 px、
        //         行首 8 px、行末 6 px），**不写数字**——数字在顶部固定标尺里。
        let ruler_y = ly + ruler_dy;
        painter.line_segment(
            [Pos2::new(row_left, ruler_y), Pos2::new(row_right, ruler_y)],
            Stroke::new(1.0, theme::SEQ_RULER),
        );
        if nucleotide {
            for i in 0..=row_units_here {
                let unit = start + i + 1;
                let h = if i == 0 {
                    8.0
                } else if i == row_units_here {
                    6.0
                } else if unit % 5 == 0 {
                    4.0
                } else {
                    2.0
                };
                let x = if i == row_units_here {
                    row_right
                } else {
                    center_x(i)
                };
                painter.line_segment(
                    [Pos2::new(x, ruler_y), Pos2::new(x, ruler_y - h * z)],
                    Stroke::new(1.0, theme::SEQ_RULER),
                );
            }
        } else {
            for i in 0..row_units_here {
                let unit = start + i + 1;
                let h = if i == 0 || i + 1 == row_units_here {
                    TICK_MAJ_H + 2.0
                } else if unit % 5 == 0 {
                    TICK_MIN_H * 0.6
                } else {
                    continue;
                };
                let x = center_x(i);
                painter.line_segment(
                    [Pos2::new(x, ruler_y), Pos2::new(x, ruler_y + h * z)],
                    Stroke::new(1.0, theme::SEQ_RULER),
                );
            }
            // 行末刻度（本行最后一个残基中心）
            let xe = center_x(row_units_here - 1);
            painter.line_segment(
                [
                    Pos2::new(xe, ruler_y),
                    Pos2::new(xe, ruler_y + (TICK_MAJ_H + 2.0) * z),
                ],
                Stroke::new(1.0, theme::SEQ_RULER),
            );
            // 数字：十进位 + 首行的 `1` + 非十进位行末
            let mut nums: Vec<usize> = (start / 10..=end / 10)
                .map(|g| g * 10)
                .filter(|n| *n > start && *n <= end)
                .collect();
            if start == 0 {
                nums.insert(0, 1);
            }
            if end % 10 != 0 {
                nums.push(end);
            }
            for n in nums {
                let x = center_x(n - 1 - start);
                painter.text(
                    Pos2::new(x - 4.0 * z, ruler_y + 2.0 * z),
                    Align2::RIGHT_TOP,
                    nfmt_plain(n),
                    ruler_font.clone(),
                    theme::TEXT,
                );
            }
        }
        // 核苷酸：行末绝对坐标写在**右侧空白区**（oracle：右对齐在视图右缘 −64）
        if nucleotide {
            painter.text(
                Pos2::new(rect.right() - 64.0 * z, ruler_y + 1.0 * z),
                Align2::RIGHT_TOP,
                nfmt_plain(end),
                ruler_font.clone(),
                theme::TEXT,
            );
        }

        // ---- 字母（等宽、连续）----
        // RNA 的 `.rna` 文件里存的是 DNA 字母，显示时按 SnapGene 转成 `u`。
        // 选中特征的残基先铺一层蓝底（oracle 实测 `#B5D5FF`，整段覆盖字母带；
        // SG 的选中效果 = 残基蓝底 + 标签蓝底白字 + 条描边变蓝）。
        //
        // **Site 不给残基加底**：oracle 实测 SG 选中 Site 只把括线轴与标签变蓝，字母底色不动
        // （用户 2026-10-19 裁定「改成与 SG 完全一致」，此前我们加了底，属刻意差异，已撤销）。
        if let Some(sfi) = app.selected_feature {
            if let Some(f) = feats.iter().find(|f| f.idx == sfi && !f.site) {
                for (a, b) in visible_spans(f, start, end) {
                    painter.rect_filled(
                        Rect::from_min_max(
                            Pos2::new(gx(a), ly - 6.0 * z),
                            Pos2::new(gx(b), ly + 13.0 * z),
                        ),
                        0.0,
                        theme::FEAT_SEL_BG,
                    );
                }
            }
        }
        let letter_y = ly + LETTER_CENTER_DY * z;
        for (i, &b) in seq[start..end].iter().enumerate() {
            painter.text(
                Pos2::new(center_x(i), letter_y),
                Align2::CENTER_CENTER,
                (kind.display_base(b) as char).to_string(),
                mono.clone(),
                theme::SEQ_LINE,
            );
        }

        // ---- 行首 `5'`（只在分子起点）----
        if nucleotide && row == 0 {
            painter.text(
                Pos2::new(content_left - FIVE_PRIME_DX * z, letter_y),
                Align2::RIGHT_CENTER,
                "5'",
                mono.clone(),
                theme::SEQ_LINE,
            );
        }

        // ---- 特征条 ----
        let bars_y0 = ly + BARS_TOP * z;
        for (fi, f) in feats.iter().enumerate() {
            let p = &plan.feats[fi];
            if f.site || p.lane == usize::MAX || p.spans.is_empty() {
                continue;
            }
            let by = bars_y0 + p.lane as f32 * LANE_DY * z;
            // 画笔关掉 ⇒ 所有特征条画成 SnapGene 那个近黑色（条内文字转白）
            let color = if app.show_colors {
                f.color
            } else {
                theme::FEATURE_NO_COLOR
            };
            let hovered = app.selected_feature != Some(f.idx) && hover_feat == Some(f.idx);
            for (a, b) in &p.spans {
                let r =
                    Rect::from_min_max(Pos2::new(gx(*a), by), Pos2::new(gx(*b), by + BAR_H * z));
                // 悬停条 = **红色描边 + 条内文字变红**（oracle `tt-t1-600-615.png`：
                // 红框 `(253,1,1)` 包住整条、条名也转红；底色不变，不是我们此前以为的蓝底）。
                draw_box(
                    &painter,
                    r,
                    color,
                    app.selected_feature == Some(f.idx),
                    hovered,
                );
                // 名字：放得下就写在条内（深色字）
                let tw = measure(&f.name, SMALL_FONT * z);
                if tw + 8.0 * z <= r.width() {
                    painter.text(
                        r.center(),
                        Align2::CENTER_CENTER,
                        &f.name,
                        small.clone(),
                        if hovered {
                            theme::SITE_HOVER
                        } else {
                            ink_on(color)
                        },
                    );
                }
                // 续行指示：条延伸到本行末而分子还没结束 ⇒ 灰色 `⟩`
                // （oracle：只有**核苷酸**画这个箭头——RNA 的 rRNA 条右端有，
                //  蛋白的 Mature chain / PRK09206 拉到行末却没有）
                if directional && *b == row_units_here && end < len {
                    let c = Pos2::new(r.right() + 10.0 * z, r.center().y);
                    draw_chevron(&painter, c, BAR_H * z * 0.42, z, f.forward);
                }
            }
        }

        // ---- 条下标签（细引线连回条上）----
        let bars_bottom = bars_y0 + plan.bar_lanes.max(1) as f32 * LANE_DY * z;
        for lab in &plan.labels {
            let ly2 = bars_bottom + LABEL_GAP * z + lab.row as f32 * LABEL_DY * z;
            let lx = content_left + lab.x;
            let leader_x = content_left + lab.leader_x;
            let selected = app.selected_feature == Some(lab.idx);
            // oracle `tt-y1c.png`：悬停条或条下标签 ⇒ 标签**红字 + 红圆角框**、
            // 引线也红（与位点标签同一套悬停色）。
            let hovered = !selected && hover_feat == Some(lab.idx);
            // 引线：条底 → 标签上沿（正下方的标签不画）
            if (lx + 1.0 - leader_x).abs() > 4.0 * z || lab.row > 0 {
                painter.line_segment(
                    [
                        Pos2::new(leader_x, ly + lab.leader_y),
                        Pos2::new(leader_x, ly2 + 1.0 * z),
                    ],
                    Stroke::new(
                        if selected {
                            1.4 * z
                        } else if hovered {
                            1.0 * z
                        } else {
                            0.8
                        },
                        if selected {
                            theme::FEAT_SEL_EDGE
                        } else if hovered {
                            theme::SITE_HOVER
                        } else {
                            theme::with_alpha(theme::SITE_AXIS, 150)
                        },
                    ),
                );
            }
            if selected {
                // oracle：选中标签 = 蓝底白字圆角块
                let tw = measure(&lab.text, SMALL_FONT * z);
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(lx - 5.0 * z, ly2 - 2.5 * z),
                        Pos2::new(lx + tw + 5.0 * z, ly2 + SMALL_FONT * z + 3.0 * z),
                    ),
                    3.0 * z,
                    theme::FEAT_SEL_EDGE,
                );
            } else if hovered {
                let tw = measure(&lab.text, SMALL_FONT * z);
                painter.rect_stroke(
                    Rect::from_min_max(
                        Pos2::new(lx - 5.0 * z, ly2 - 2.5 * z),
                        Pos2::new(lx + tw + 5.0 * z, ly2 + SMALL_FONT * z + 3.0 * z),
                    ),
                    3.0 * z,
                    Stroke::new(1.2, theme::SITE_HOVER),
                    egui::StrokeKind::Inside,
                );
            }
            painter.text(
                Pos2::new(lx, ly2),
                Align2::LEFT_TOP,
                &lab.text,
                small.clone(),
                if selected {
                    Color32::WHITE
                } else if hovered {
                    theme::SITE_HOVER
                } else {
                    theme::TEXT
                },
            );
        }
    }

    // ---- 底部：全序列概览条 + 视口指示 ----
    draw_bottom(&painter, rect, app, &feats, len, total_h, view_h, z);

    // ---- 悬停工具卡 ----
    // oracle（`/tmp/opencode/sg/tt-*.png`）：SG 在序列视图里悬停**特征条**或
    // **位点标签**都会弹白底卡片（标题 + 类型芯片 + 区间 + 段列表 + NOTE）。
    // 悬停字母（无特征）不弹；RNA/DNA 用同一版式（单位换 `bp`、前缀留空）。
    if let Some(fi) = hover_feat {
        if let Some(f) = feats.iter().find(|f| f.idx == fi) {
            let card = feat_card::Card {
                name: &f.name,
                kind: &f.kind,
                sub_type: &f.sub_type,
                note: &f.note,
                segs: &f.segs,
                prefix: if nucleotide { "" } else { "p" },
                // 单位随分子类型（DNA `bp`、RNA `bases`、蛋白 `aa`，与状态栏同一个来源）。
                unit: doc.unit(),
                // oracle：`Region` 卡片的标题行右侧有特征色块；`Site` 没有。
                swatch: (!f.site).then_some(f.color),
                // 蛋白没有方向箭头（oracle `tt-t1-card.png` 标题行只有色块）。
                direction: None,
            };
            egui::Tooltip::for_enabled(&resp)
                .at_pointer()
                .gap(14.0)
                .show(|ui| {
                    let frame = feat_card::frame();
                    frame.show(ui, |ui| {
                        ui.style_mut().override_font_id =
                            Some(egui::FontId::proportional(11.5));
                        feat_card::show(ui, &card);
                    });
                });
        }
    }

    rect
}

/// 序列视图的行宽（单位数）：`None` = 换行到窗口（自适应，取 10 的整数倍）。
///
/// oracle：SnapGene 的序列视图本来就按窗口宽度换行（1230 px 窗口里每行 60 个单位 =
/// 内容宽 940 − 左 130 − 右 150），2K 桌面上必须铺满视图区（用户 2026-10-18 指摘）。
fn fit_row_units(avail: f32, base_w: f32, fixed: Option<usize>) -> usize {
    match fixed {
        Some(n) => n.clamp(10, 2000),
        None => (((avail / base_w).floor() as usize) / 10 * 10).clamp(10, 2000),
    }
}

/// 本行的排版：分段、车道、标签行与行高。
#[allow(clippy::too_many_arguments)]
fn plan_row(
    feats: &[FlatFeat],
    start: usize,
    end: usize,
    base_w: f32,
    z: f32,
    content_w: f32,
    measure: &dyn Fn(&str, f32) -> f32,
) -> RowPlan {
    let mut row: Vec<RowFeat> = feats
        .iter()
        .map(|f| RowFeat {
            spans: visible_spans(f, start, end),
            lane: usize::MAX,
        })
        .collect();

    // 车道贪心：Site 与普通条各自编号（与已放置的条重叠就换下一层）。
    //
    // Site 的**占位**是它在这一行的整条轴线，而不是各个小片段：SnapGene 把
    // 延续到后续行的 Site 轴线一直拉到行末（oracle：行 1 的 `Site` 与
    // `domain interfaces [active]` 片段互不重叠，却仍各占一条车道）。
    let site_spans: Vec<Option<(usize, usize)>> = (0..feats.len())
        .map(|i| {
            let f = &feats[i];
            let spans = &row[i].spans;
            if spans.is_empty() {
                return None;
            }
            let a = spans[0].0;
            let continues = f.segs.iter().any(|(s, _)| *s >= end);
            let b = if continues {
                end - start
            } else {
                spans[spans.len() - 1].1
            };
            Some((a, b.max(a + 1)))
        })
        .collect();
    let mut assign = |is_site: bool| -> usize {
        let mut lanes: Vec<Vec<(usize, usize)>> = Vec::new();
        for (i, f) in feats.iter().enumerate() {
            if f.site != is_site || row[i].spans.is_empty() {
                continue;
            }
            let place: Vec<(usize, usize)> = if is_site {
                site_spans[i].into_iter().collect()
            } else {
                row[i].spans.clone()
            };
            let mut lv = 0usize;
            'search: loop {
                if lv >= lanes.len() {
                    lanes.push(Vec::new());
                }
                for (a, b) in &place {
                    if lanes[lv].iter().any(|(c, d)| a < d && c < b) {
                        lv += 1;
                        continue 'search;
                    }
                }
                break;
            }
            lanes[lv].extend(place.iter().copied());
            row[i].lane = lv;
        }
        lanes.len()
    };
    let bar_lanes = assign(false);
    let site_lanes = assign(true);

    // ---- Site 标签：轴区最上面一行，中心对齐各自的「车道线」中点，重叠时向右让位 ----
    let mut site_labels: Vec<SiteLabel> = Vec::new();
    if site_lanes > 0 {
        let mut items: Vec<(f32, String, usize)> = Vec::new();
        for (i, f) in feats.iter().enumerate() {
            if !f.site || row[i].spans.is_empty() {
                continue;
            }
            let spans = &row[i].spans;
            let cx0 = (spans[0].0 as f32 + 0.5) * base_w;
            let continues = f.segs.iter().any(|(a, _)| *a >= end);
            let cx1 = if continues {
                content_w
            } else {
                (spans[spans.len() - 1].1 as f32 - 0.5) * base_w
            };
            items.push(((cx0 + cx1) * 0.5, f.name.clone(), f.idx));
        }
        items.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut cursor = 0.0f32;
        for (cx, text, idx) in items {
            let tw = measure(&text, SITE_FONT * z);
            let mut x0 = cx - tw * 0.5;
            if x0 < cursor {
                x0 = cursor;
            }
            if x0 + tw > content_w {
                x0 = (content_w - tw).max(0.0);
            }
            cursor = x0 + tw + 10.0 * z;
            site_labels.push(SiteLabel {
                cx: x0 + tw * 0.5,
                text,
                idx,
            });
        }
    }

    // ---- 条下标签：名字放不进条里就排到条下方（贪心，从左到右）----
    let small = SMALL_FONT * z;
    let mut labels: Vec<LabelPlan> = Vec::new();
    let mut label_rows = 0usize;
    let mut cands: Vec<(f32, f32, String, usize)> = Vec::new();
    for (i, f) in feats.iter().enumerate() {
        if f.site || row[i].spans.is_empty() {
            continue;
        }
        let widest = row[i]
            .spans
            .iter()
            .map(|(a, b)| (*b - *a) as f32 * base_w)
            .fold(0.0f32, f32::max);
        if measure(&f.name, small) + 8.0 * z <= widest {
            continue; // 写进条里就够了
        }
        let (a, b) = row[i].spans[0];
        let leader_y = (BARS_TOP + row[i].lane as f32 * LANE_DY) * z + BAR_H * z;
        cands.push((
            (a + b) as f32 * 0.5 * base_w,
            leader_y,
            f.name.clone(),
            f.idx,
        ));
    }
    cands.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut rows: Vec<Vec<(f32, f32)>> = vec![Vec::new(); MAX_LABEL_ROWS];
    for (cx, leader_y, text, idx) in cands {
        let tw = measure(&text, small);
        let mut placed = None;
        for (ri, spans) in rows.iter_mut().enumerate() {
            let x = (cx - tw * 0.5).clamp(0.0, (content_w - tw).max(0.0));
            if spans
                .iter()
                .all(|(a, b)| x > *b + 6.0 * z || x + tw < *a - 6.0 * z)
            {
                spans.push((x, x + tw));
                placed = Some((ri, x));
                break;
            }
        }
        let Some((ri, x)) = placed else { continue };
        label_rows = label_rows.max(ri + 1);
        labels.push(LabelPlan {
            row: ri,
            x,
            leader_x: cx,
            leader_y,
            text,
            idx,
        });
    }

    let axis_h = if site_lanes > 0 {
        (AXIS_GAP + (site_lanes - 1) as f32 * SITE_LANE_DY + SITE_LABEL_UP + SITE_LABEL_H + 2.0) * z
    } else {
        0.0
    };
    let bars_h = if bar_lanes > 0 {
        (bar_lanes - 1) as f32 * LANE_DY * z + BAR_H * z
    } else {
        0.0
    };
    let labels_h = if label_rows > 0 {
        (LABEL_GAP + label_rows as f32 * LABEL_DY) * z
    } else {
        0.0
    };
    let body_h = BARS_TOP * z + bars_h + labels_h;
    let row_gap = if site_lanes > 0 {
        ROW_GAP_PROT_SITE
    } else {
        ROW_GAP
    } * z;
    RowPlan {
        feats: row,
        site_lanes,
        bar_lanes,
        labels,
        site_labels,
        axis_h,
        body_h,
        height: axis_h + body_h + row_gap,
    }
}

/// 特征在本行内的可见分段（**行内** 0-based 半开坐标）。
fn visible_spans(f: &FlatFeat, start: usize, end: usize) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    for (a, b) in &f.segs {
        let s = (*a).max(start);
        let e = (*b).min(end);
        if s < e {
            v.push((s - start, e - start));
        }
    }
    v
}

/// 条底色的墨色（浅底黑字 / 深底白字）。
fn ink_on(fill: Color32) -> Color32 {
    if theme::luma(fill) > 0.55 {
        theme::TEXT
    } else {
        Color32::WHITE
    }
}

/// 一条特征条：填充 + 细描边（SnapGene 的条就是描边矩形，没有内部分界线）。
fn draw_box(p: &egui::Painter, r: Rect, color: Color32, selected: bool, hovered: bool) {
    p.rect_filled(r, 0.0, color);
    let edge = if selected {
        Stroke::new(1.6, theme::ACCENT)
    } else if hovered {
        Stroke::new(1.6, theme::SITE_HOVER)
    } else {
        Stroke::new(1.0, theme::BAR_EDGE)
    };
    p.rect_stroke(r, 0.0, edge, egui::StrokeKind::Inside);
}

/// 续行指示（SnapGene 在条尾外侧画的灰色 `⟩`）。
fn draw_chevron(p: &egui::Painter, c: Pos2, h: f32, z: f32, forward: bool) {
    let w = h * 0.62;
    let dir = if forward { 1.0 } else { -1.0 };
    let stroke = Stroke::new(1.5 * z.max(0.6), theme::CHEVRON);
    let tip = Pos2::new(c.x + w * dir, c.y);
    let back = Pos2::new(c.x - w * 0.35 * dir, c.y);
    p.add(egui::Shape::line(
        vec![
            Pos2::new(back.x, c.y - h),
            tip,
            Pos2::new(back.x, c.y + h),
        ],
        stroke,
    ));
}

/// 底部概览条：整条序列按比例映射，画特征色块 + 视口框。
#[allow(clippy::too_many_arguments)]
fn draw_bottom(
    p: &egui::Painter,
    rect: Rect,
    app: &HelixApp,
    feats: &[FlatFeat],
    len: usize,
    total_h: f32,
    view_h: f32,
    z: f32,
) {
    let bar = Rect::from_min_max(
        Pos2::new(rect.left() + 24.0, rect.bottom() - BOTTOM_H + 6.0),
        Pos2::new(rect.right() - 24.0, rect.bottom() - BOTTOM_H + 6.0 + 13.0 * z),
    );
    p.rect_filled(bar, 2.0, theme::PANEL_ALT);
    p.rect_stroke(
        bar,
        2.0,
        Stroke::new(0.8, theme::BORDER),
        egui::StrokeKind::Middle,
    );
    let x_of = |u: usize| -> f32 { bar.left() + (u as f32 / len.max(1) as f32) * bar.width() };
    for f in feats {
        for (a, b) in &f.segs {
            let x0 = x_of(*a);
            let x1 = x_of(*b).max(x0 + 1.0);
            p.rect_filled(
                Rect::from_min_max(
                    Pos2::new(x0, bar.top() + 2.0),
                    Pos2::new(x1, bar.bottom() - 2.0),
                ),
                1.0,
                f.color,
            );
        }
    }
    // 视口指示
    if total_h > view_h {
        let top_u = (app.seq.scroll_y / total_h * len as f32) as usize;
        let bot_u = ((app.seq.scroll_y + view_h) / total_h * len as f32) as usize;
        let x0 = x_of(top_u.min(len));
        let x1 = x_of(bot_u.min(len)).max(x0 + 3.0);
        p.rect_stroke(
            Rect::from_min_max(Pos2::new(x0, bar.top()), Pos2::new(x1, bar.bottom())),
            2.0,
            Stroke::new(1.2, theme::ACCENT),
            egui::StrokeKind::Middle,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 行宽自适应：`None` 按可用宽度取 10 的整数倍，`Some(n)` 用固定值。
    #[test]
    fn fit_row_units_follows_window_width() {
        // oracle：SG 在 1230 px 窗口（内容宽 940 = 左 130 + 60 个单位 + 右 150）
        assert_eq!(fit_row_units(660.0, 11.0, None), 60);
        // 零头向下取整到 10 的倍数，别溢出行末坐标
        assert_eq!(fit_row_units(702.9, 11.0, None), 60);
        assert_eq!(fit_row_units(1105.0, 11.0, None), 100);
        // 太窄也要至少画 10 个
        assert_eq!(fit_row_units(12.0, 11.0, None), 10);
        // 缩放后步距变大，同一宽度能放的更少
        assert_eq!(fit_row_units(660.0, 22.0, None), 30);
        // 固定行宽优先（SnapGene 的 Line width 下拉）
        assert_eq!(fit_row_units(660.0, 11.0, Some(120)), 120);
        assert_eq!(fit_row_units(660.0, 11.0, Some(3)), 10);
    }

    fn feat(idx: usize, segs: &[(usize, usize)]) -> FlatFeat {
        FlatFeat {
            idx,
            name: format!("f{idx}"),
            kind: "misc_feature".into(),
            sub_type: String::new(),
            note: String::new(),
            site: false,
            color: Color32::RED,
            forward: true,
            segs: segs.to_vec(),
        }
    }

    fn site_feat(idx: usize, segs: &[(usize, usize)]) -> FlatFeat {
        let mut f = feat(idx, segs);
        f.site = true;
        f
    }

    /// 测试用的「字宽」估计：每字符 6 px。
    fn measure(text: &str, _size: f32) -> f32 {
        text.chars().count() as f32 * 6.0
    }

    fn plan(fs: &[FlatFeat], start: usize, end: usize) -> RowPlan {
        plan_row(fs, start, end, 11.0, 1.0, 660.0, &measure)
    }

    #[test]
    fn visible_spans_are_row_local() {
        let f = feat(0, &[(0, 1541)]);
        assert_eq!(visible_spans(&f, 60, 120), vec![(0, 60)]);
        // 第二行从 0 开始（行内坐标）
        assert_eq!(visible_spans(&f, 120, 180), vec![(0, 60)]);
    }

    #[test]
    fn overlapping_features_get_separate_lanes() {
        let fs = vec![
            feat(0, &[(10, 40)]),
            feat(1, &[(20, 50)]),
            feat(2, &[(60, 80)]),
        ];
        let p = plan(&fs, 0, 100);
        assert_eq!(p.bar_lanes, 2);
        assert_eq!(p.feats[0].lane, 0);
        assert_eq!(p.feats[1].lane, 1);
        // 与 0/1 都不重叠 ⇒ 回到第 0 层
        assert_eq!(p.feats[2].lane, 0);
    }

    #[test]
    fn sites_get_their_own_lanes_and_axis_height() {
        let fs = vec![feat(0, &[(10, 40)]), site_feat(1, &[(20, 30)])];
        let p = plan(&fs, 0, 100);
        // 普通条 1 车道、位点 1 车道；位点不占条的车道
        assert_eq!(p.bar_lanes, 1);
        assert_eq!(p.site_lanes, 1);
        assert_eq!(p.feats[1].lane, 0);
        assert!(p.axis_h > 0.0);
        assert_eq!(p.site_labels.len(), 1);
    }

    #[test]
    fn invisible_features_have_no_lane() {
        let fs = vec![feat(0, &[(10, 20)])];
        let p = plan(&fs, 100, 200);
        assert_eq!(p.feats[0].lane, usize::MAX);
        assert_eq!(p.bar_lanes, 0);
        assert!(p.body_h > 0.0);
    }

    #[test]
    fn narrow_feature_gets_a_label_row() {
        // 1 个残基 = 11 px 宽，名字 `f0`（2 字符 = 12 px + 8 边距）放不进 ⇒ 排标签
        let fs = vec![feat(0, &[(10, 11)])];
        let p = plan(&fs, 0, 100);
        assert_eq!(p.labels.len(), 1);
        // 行高里含一条标签行（38 + 17 + 6 + 16）
        assert!((p.body_h - (38.0 + 17.0 + 6.0 + 16.0)).abs() < 1e-3);
        // 宽条：名字放得进条里 ⇒ 不占标签行
        let fs2 = vec![feat(0, &[(0, 60)])];
        let p2 = plan(&fs2, 0, 100);
        assert!(p2.labels.is_empty());
        assert!((p2.body_h - (38.0 + 17.0)).abs() < 1e-3);
    }
}
