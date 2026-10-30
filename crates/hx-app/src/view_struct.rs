//! RNA 文档的 **`Secondary Structure`** 页（对齐 SnapGene 8.2.3 的实测版面）。
//!
//! 版面（oracle `/tmp/opencode/sg/r2-tRNA-struct.png`、`r20-eco-struct.png`）：
//!
//! - 画布 `#F2F2F2`；左上 `Confidence:` 四色图例；右侧白色**次优结构列表**
//!   （缩略图 + `kcal/mol`，选中项蓝底）；底部五格统计条（第一格 `#CBCBCB`，其余 `#DFDFDF`）。
//! - 碱基画成**彩色空心环**（颜色 = 该碱基的置信度）+ 黑色粗体字母；骨架线 `#252525`
//!   画在环上面、只在环与环之间露出来（环里没有线穿过字母）。
//! - 配对的两碱基之间有一条**浅灰细横档**（氢键，实测 100% 约 1.2 px、`#a8a8a8`）：
//!   它把茎画成梯子，也是缩小后仍能读出配对关系的依据。
//! - 每 10 个碱基一个序号标签 + 细引线；5′/3′ 是深色圆角小标签。
//!
//! 数据完全来自 block 34（`hx_formats::rna_struct`），几何来自 `hx_render::RnaLayout`。

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

use hx_formats::rna_struct::{RnaStructure, StructStats};
use hx_render::rna_layout::RnaLayout;
use hx_render::{RING_RADIUS, RING_STROKE};

use crate::app::HelixApp;
use crate::theme;

/// 显示字母：`.rna` 文件里存的是 DNA 字母，SnapGene 的 RNA 视图显示 `U`。
#[inline]
fn nt(b: u8) -> char {
    // 一律大写：SnapGene 的结构页不保留序列里的大小写（`.rna` 里的软屏蔽
    // 小写只影响序列视图的配色，不影响结构页）。
    match b.to_ascii_uppercase() {
        b'T' => 'U',
        other => other as char,
    }
}

/// 画布底色（SnapGene 的二级结构页）。
const BG: Color32 = Color32::from_rgb(0xf2, 0xf2, 0xf2);
/// 置信度四档（`Confidence:` 图例与碱基环同色）。
pub const CONF_COLORS: [Color32; 4] = [
    Color32::from_rgb(0xf8, 0x00, 0x00),
    Color32::from_rgb(0xff, 0xa8, 0x00),
    Color32::from_rgb(0x00, 0xb4, 0xf2),
    Color32::from_rgb(0x00, 0x37, 0xff),
];
pub const CONF_LABELS: [&str; 4] = [
    "90% and greater",
    "70% - 89%",
    "50% - 69%",
    "Less than 50%",
];
/// 碱基字母 / 骨架 / 序号文字。
const INK: Color32 = Color32::from_rgb(0x25, 0x25, 0x25);
/// 配对横档（氢键）的浅灰。
const RUNG: Color32 = Color32::from_rgb(0xa8, 0xa8, 0xa8);
/// 5′/3′ 标签底色。
const BADGE: Color32 = Color32::from_rgb(0x3a, 0x3a, 0x3a);
/// 次优结构卡片（选中）。
const CARD_ON: Color32 = Color32::from_rgb(0x00, 0x66, 0xec);
/// 统计条（第一格更深：SnapGene 用它标出「当前显示的结构」）。
const STATS_BG: Color32 = Color32::from_rgb(0xdf, 0xdf, 0xdf);
const STATS_FIRST: Color32 = Color32::from_rgb(0xcb, 0xcb, 0xcb);

/// 视图顶部留白（与视图页签、右上工具条同一行）。
const TOP: f32 = 44.0;
/// 右侧次优结构列的宽度（含滚动条）。
const THUMB_W: f32 = 113.0;
/// 底部统计条高度。
const STATS_H: f32 = 67.0;

/// 置信度 → 图例颜色。
pub fn conf_color(p: f64) -> Color32 {
    if p >= 0.90 {
        CONF_COLORS[0]
    } else if p >= 0.70 {
        CONF_COLORS[1]
    } else if p >= 0.50 {
        CONF_COLORS[2]
    } else {
        CONF_COLORS[3]
    }
}

/// 某碱基的置信度：配对碱基用「这一对」的概率；不配对的用 `1 - Σp`。
pub(crate) fn base_conf(st: &RnaStructure, i: usize, partner: Option<usize>) -> f64 {
    match partner {
        Some(j) => {
            let p = st.pair_prob(i, j);
            if p > 0.0 {
                p
            } else {
                st.paired_prob(i)
            }
        }
        None => st.unpaired_prob(i),
    }
}

/// 逐碱基配色（显示结构 + 配对表）。
pub(crate) fn confidence_colors(st: &RnaStructure, len: usize, partner: &[Option<usize>]) -> Vec<Color32> {
    (0..len)
        .map(|i| conf_color(base_conf(st, i, partner[i])))
        .collect()
}

/// 配对表：`partner[i] = Some(j)`。
pub(crate) fn partner_map(len: usize, pairs: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut v = vec![None; len];
    for &(i, j) in pairs {
        if i < len && j < len {
            v[i] = Some(j);
            v[j] = Some(i);
        }
    }
    v
}

/// 屏幕坐标版的「最空方向」（规则在 `hx_render::rna_layout::free_dir`，
/// 与 SVG 导出共用）。
fn free_dir(scr: &[hx_render::rna_layout::Pt], i: usize, away: Vec2, dist: f32) -> Vec2 {
    let d = hx_render::rna_layout::free_dir(
        scr,
        i,
        hx_render::rna_layout::Pt::new(away.x, away.y),
        dist,
    );
    Vec2::new(d.x, d.y)
}

fn bold_text(p: &egui::Painter, pos: Pos2, text: &str, font: FontId, color: Color32, bold: f32) {
    if bold > 0.0 {
        p.text(
            pos + Vec2::new(bold, 0.0),
            Align2::CENTER_CENTER,
            text,
            font.clone(),
            color,
        );
    }
    p.text(pos, Align2::CENTER_CENTER, text, font, color);
}

/// 画一页（返回画布矩形，供右上工具条定位）。
pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) -> egui::Rect {
    let (resp, p) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    if rect.width() < 200.0 || rect.height() < 180.0 {
        return rect;
    }
    p.rect_filled(rect, 0.0, Color32::WHITE);

    // 版面：画布 | 右列（次优结构）
    //       统计条（只到右列左边）
    let canvas = Rect::from_min_max(
        rect.min + Vec2::new(0.0, TOP),
        Pos2::new(rect.right() - THUMB_W, rect.bottom() - STATS_H),
    );
    let thumbs = Rect::from_min_max(
        Pos2::new(rect.right() - THUMB_W, rect.top() + TOP),
        rect.right_bottom(),
    );
    let stats = Rect::from_min_max(
        Pos2::new(rect.left(), rect.bottom() - STATS_H),
        Pos2::new(rect.right() - THUMB_W, rect.bottom()),
    );
    p.rect_filled(canvas, 0.0, BG);
    p.rect_filled(thumbs, 0.0, Color32::WHITE);
    p.rect_filled(stats, 0.0, STATS_BG);

    // ---- 界面状态（`doc` 借用期间不能再写 app） ----
    let mut zoom = app.struct_zoom;
    let mut pan = app.struct_pan;
    let mut sel = app.struct_sel;
    let mut sel_base = app.struct_sel_base;
    let show_index = app.struct_show_index;

    let Some(doc) = app.doc.as_mut() else { return rect };
    let len = doc.len();
    let seq: Vec<u8> = doc.sequence().to_vec();
    // 结构（以及置信度着色）直接来自块里的 `probabilityMatrix`，不会再有
    // 「先画结构、等概率算完再变色」的中间态。
    let cache = match doc.rna_structure() {
        Ok(crate::doc::RnaCache::Ready(c)) => c,
        Ok(crate::doc::RnaCache::Folding) => {
            // 本地折叠在后台线程里（大序列要几秒），算完再画。
            centered_note(ui, canvas, "Folding secondary structure…");
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_millis(120));
            return rect;
        }
        Ok(crate::doc::RnaCache::None) => {
            centered_note(ui, canvas, "No secondary structure is stored in this file.");
            return rect;
        }
        Err(e) => {
            centered_note(ui, canvas, &format!("Cannot read the secondary structure block:\n{e}"));
            return rect;
        }
    };
    if let Err(e) = &cache.structure {
        centered_note(ui, canvas, &format!("Cannot decode the secondary structure:\n{e}"));
        return rect;
    }
    // `structures()` 会克隆配对表；缓存起来只在数据变化时算一次。
    cache.ensure_structures(len);
    let stats_data: Option<StructStats> = cache.structure.as_ref().ok().and_then(|s| s.stats);
    let count = cache.structures.as_ref().map_or(0, |v| v.len());
    if count == 0 {
        centered_note(ui, canvas, "The secondary structure block is empty.");
        return rect;
    }
    if sel >= count {
        sel = 0;
    }
    cache.ensure_layout(sel, len);
    cache.ensure_thumbs(len);
    let structures = cache.structures.as_ref().expect("structures");
    let Some((_, layout)) = cache.layout.as_ref() else {
        return rect;
    };
    let thumbs_layout: &[RnaLayout] = cache.thumbs.as_deref().unwrap_or(&[]);
    let st: &RnaStructure = match &cache.structure {
        Ok(s) => s,
        Err(_) => return rect,
    };
    let partner = partner_map(len, &structures[sel].0);
    let colors = confidence_colors(st, len, &partner);

    // ---- 适配（Zoom to fit） ----
    let (bmin, bmax) = layout.bounds();
    let bw = (bmax.x - bmin.x).max(1.0);
    let bh = (bmax.y - bmin.y).max(1.0);
    // SnapGene 的 Zoom to fit 只按**高度**算，结构占画布高度的 72%（实测：
    // tRNA 68.46%、16S 7.51%，两者高度比都是 0.72），宽度只做溢出保护。
    let fit = (0.72 * canvas.height() / bh)
        .min(0.98 * canvas.width() / bw)
        .clamp(0.002, 4.0);
    if zoom <= 0.0 {
        zoom = fit;
        pan = Vec2::ZERO;
    }
    // ---- 输入：滚轮缩放（以指针为中心）、拖拽平移、点碱基 ----
    let pointer = ui.input(|i| i.pointer.hover_pos());
    if resp.hovered() {
        if let Some(q) = pointer.filter(|q| canvas.contains(*q)) {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll.abs() > 0.01 {
                let f = (scroll * 0.0028).exp();
                let nz = (zoom * f).clamp(fit * 0.25, 6.0);
                let center = canvas.center();
                pan = q - center - (q - center - pan) * (nz / zoom);
                zoom = nz;
            }
        }
    }
    if resp.dragged_by(egui::PointerButton::Primary) {
        let d = resp.drag_delta();
        if d != Vec2::ZERO {
            pan += d;
        }
    }
    let scale = zoom;
    let c0 = hx_render::rna_layout::Pt::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let center = canvas.center();
    let scr: Vec<Pos2> = layout
        .pos
        .iter()
        .map(|q| {
            Pos2::new(
                center.x + pan.x + (q.x - c0.x) * scale,
                center.y + pan.y + (q.y - c0.y) * scale,
            )
        })
        .collect();
    let scr_pt: Vec<hx_render::rna_layout::Pt> = scr
        .iter()
        .map(|q| hx_render::rna_layout::Pt::new(q.x, q.y))
        .collect();
    let radius = RING_RADIUS * scale;
    let hit_lim = radius.max(9.0);
    let hit = |q: Pos2| -> Option<usize> {
        let mut best: Option<(f32, usize)> = None;
        for (i, s) in scr.iter().enumerate() {
            let d = s.distance(q);
            if d <= hit_lim && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, i));
            }
        }
        best.map(|(_, i)| i)
    };
    let hover = pointer
        .filter(|q| canvas.contains(*q))
        .and_then(hit);
    if resp.clicked() {
        if let Some(q) = pointer.filter(|q| canvas.contains(*q)) {
            sel_base = match hit(q) {
                Some(i) if sel_base == Some(i) => None,
                other => other,
            };
        }
    }

    // ---- 画结构 ----
    let cp = p.with_clip_rect(canvas.intersect(rect));
    // 环 + 字母
    let letter_font = FontId::proportional((hx_render::rna_layout::LETTER_SIZE * scale).max(4.0));
    for i in 0..len {
        let col = colors[i];
        if radius < 2.2 {
            // 缩得太小：直接画实心小点（SnapGene 在 7.5% 也是这样）
            cp.circle_filled(scr[i], radius.max(0.9), col);
        } else {
            // SnapGene 的碱基圈**不填色**：骨架线只在环之间露出来。
            cp.circle_stroke(scr[i], radius, Stroke::new((RING_STROKE * scale).max(0.5), col));
        }
        if radius >= 5.0 {
            let ch = nt(seq.get(i).copied().unwrap_or(b'?')).to_string();
            bold_text(
                &cp,
                scr[i] + Vec2::new(0.0, radius * 0.04),
                &ch,
                letter_font.clone(),
                INK,
                (scale * 0.55).max(0.2),
            );
        }
    }
    // 骨架 + 配对横档：画在环**之上**，两端压在环的描边上。
    // 环是空心的，若按「环外缘」截断，描边会把缝隙吃光、连线整个消失；
    // SnapGene 也是让连线盖住一点环边（但绝不穿过中间的字母）。
    let link = (radius - 1.0).max(radius * 0.5);
    let back_w = (hx_render::rna_layout::BACKBONE_WIDTH * scale).max(0.6);
    for i in 0..len.saturating_sub(1) {
        let (a, b) = (scr[i], scr[i + 1]);
        let d = b - a;
        let l = d.length();
        if l <= 2.0 * link + 0.2 {
            continue;
        }
        let u = d / l;
        cp.line_segment([a + u * link, b - u * link], Stroke::new(back_w, INK));
    }
    // 配对横档（氢键）：实测 100% 约 1.2 px、浅灰 `#a8a8a8`。
    let rung_w = (1.2 * scale).max(0.5);
    for &(i, j) in &structures[sel].0 {
        if i >= len || j >= len {
            continue;
        }
        let (a, b) = (scr[i], scr[j]);
        let d = b - a;
        let l = d.length();
        if l <= 2.0 * link + 0.2 {
            continue;
        }
        let u = d / l;
        cp.line_segment(
            [a + u * link, b - u * link],
            Stroke::new(rung_w, RUNG),
        );
    }
    // 选中 / 悬停的碱基（连同配对的另一半）
    for (i, strong) in [(sel_base, true), (hover, false)] {
        let Some(i) = i else { continue };
        if i >= len {
            continue;
        }
        for k in [Some(i), partner[i]] {
            let Some(k) = k else { continue };
            let r = radius.max(2.0);
            cp.circle_stroke(
                scr[k],
                r + 2.5,
                Stroke::new(if strong { 2.0 } else { 1.2 }, theme::ACCENT),
            );
        }
    }

    // ---- 每 10 个碱基的序号 ----
    //
    // SnapGene 的序号是**稀疏**的：缩到能看见整个结构时（大 rRNA 的 fit 只有
    // 个位数百分比）一个都不画，中等缩放也只保留不互相重叠的那些；否则一百多个
    // 数字会把结构盖成一团（用户实测 16S 在 37% 就是这样）。这里按「最小屏幕
    // 间距」取舍：先到先得，后来的挤不下就跳过。
    if show_index && scale >= 0.45 {
        let font = FontId::proportional(11.5);
        let mut placed: Vec<Pos2> = Vec::new();
        for (i, s) in scr.iter().enumerate() {
            if i % 10 != 9 {
                continue;
            }
            let li = layout.loop_of.get(i).copied().unwrap_or(0);
            let lc = layout
                .loops
                .get(li)
                .map(|l| Pos2::new(l.center.x, l.center.y))
                .unwrap_or(*s);
            let d = free_dir(&scr_pt, i, (*s - lc).normalized(), radius + 22.0);
            let gap = radius + 7.0;
            let tip = *s + d * gap;
            let anchor = *s + d * (gap + 20.0);
            if placed.iter().any(|q| q.distance(anchor) < 40.0) {
                continue;
            }
            placed.push(anchor);
            cp.line_segment([tip, tip + d * 11.0], Stroke::new(0.8, Color32::from_gray(0x8a)));
            let align = if d.x > 0.25 {
                Align2::LEFT_CENTER
            } else if d.x < -0.25 {
                Align2::RIGHT_CENTER
            } else {
                Align2::CENTER_CENTER
            };
            cp.text(anchor, align, format!("{}", i + 1), font.clone(), INK);
        }
    }

    // ---- 5′ / 3′ 标签 ----
    if len >= 2 && scale >= 0.35 {
        let badge_font = FontId::proportional(11.5);
        // 两个末端常常挨在一起（tRNA 的接受臂顶端）：沿半径、沿骨架、远离重心都可能
        // 指进结构里，所以直接找「周围最空」的方向放标签。
        let centroid = scr.iter().fold(Vec2::ZERO, |a, b| a + b.to_vec2()) / len.max(1) as f32;
        for (idx, text, away) in [
            (0usize, "5'", (scr[0].to_vec2() - centroid).normalized()),
            (len - 1, "3'", (scr[len - 1].to_vec2() - centroid).normalized()),
        ] {
            let d = free_dir(&scr_pt, idx, away, radius + 24.0);
            let at = scr[idx] + d * (radius + 24.0);
            let box_w = 22.0;
            let box_h = 16.0;
            let r = Rect::from_center_size(at, Vec2::new(box_w, box_h));
            cp.rect_filled(r, 3.5, BADGE);
            bold_text(
                &cp,
                r.center(),
                text,
                badge_font.clone(),
                Color32::WHITE,
                0.3,
            );
        }
    }

    // ---- 左上 Confidence 图例 ----
    {
        let x = canvas.left() + 10.0;
        let mut y = canvas.top() + 18.0;
        let font = FontId::proportional(12.5);
        p.text(
            Pos2::new(x, y),
            Align2::LEFT_CENTER,
            crate::i18n::text(ui.ctx(),"Confidence:"),
            font.clone(),
            INK,
        );
        y += 20.0;
        for (c, label) in CONF_COLORS.iter().zip(CONF_LABELS) {
            p.circle_filled(Pos2::new(x + 5.0, y), 4.5, *c);
            p.text(
                Pos2::new(x + 17.0, y),
                Align2::LEFT_CENTER,
                label,
                font.clone(),
                INK,
            );
            y += 20.0;
        }
    }

    // ---- 左下两个小按钮（Display / Zoom to fit） ----
    {
        let bw = 30.0;
        let top_y = stats.top() - 66.0;
        let r1 = Rect::from_min_size(Pos2::new(canvas.left() + 8.0, top_y), Vec2::splat(bw));
        let r2 = Rect::from_min_size(
            Pos2::new(canvas.left() + 8.0, top_y + bw + 4.0),
            Vec2::splat(bw),
        );
        for (r, id, icon) in [
            (r1, "struct_display_btn", crate::icons::Icon::LineWidth),
            (r2, "struct_fit_btn", crate::icons::Icon::Search),
        ] {
            p.rect_filled(r, 4.0, Color32::WHITE);
            p.rect_stroke(
                r,
                4.0,
                Stroke::new(0.7, theme::BORDER),
                egui::StrokeKind::Inside,
            );
            crate::icons::paint(&p, icon, r.shrink(6.0), theme::TEXT, Color32::WHITE);
            let btn = ui.interact(r, ui.scope_id().with(id), Sense::click());
            if btn.clicked() {
                if id == "struct_fit_btn" {
                    zoom = fit;
                    pan = Vec2::ZERO;
                }
            }
            if id == "struct_fit_btn" {
                btn.on_hover_text(crate::i18n::text(ui.ctx(),"Zoom to fit"));
            } else {
                btn.on_hover_text(crate::i18n::text(ui.ctx(),"Display options"));
            }
        }
    }

    // ---- 右侧次优结构列表 ----
    let mut clicked_card: Option<usize> = None;
    {
        let mut child = ui.new_child(
            egui::UiBuilder::new()
                .max_rect(thumbs.shrink2(Vec2::new(6.0, 6.0)))
                .layout(egui::Layout::top_down(egui::Align::Center)),
        );
        child.set_clip_rect(thumbs.shrink2(Vec2::new(6.0, 6.0)));
        egui::ScrollArea::vertical()
            .id_salt("struct_thumbs")
            .auto_shrink([false, false])
            .show(&mut child, |ui| {
                for (i, (pairs, energy)) in structures.iter().enumerate() {
                    let on = i == sel;
                    let (card, cresp) =
                        ui.allocate_exact_size(Vec2::new(88.0, 112.0), Sense::click());
                    if cresp.clicked() {
                        clicked_card = Some(i);
                    }
                    let cp = ui.painter();
                    if on {
                        cp.rect_filled(card, 6.0, CARD_ON);
                    } else {
                        cp.rect_filled(card, 6.0, Color32::WHITE);
                        cp.rect_stroke(
                            card,
                            6.0,
                            Stroke::new(0.8, Color32::from_rgb(0xe2, 0xe2, 0xe2)),
                            egui::StrokeKind::Inside,
                        );
                    }
                    // 序号圆牌
                    let badge_c = Pos2::new(card.left() + 15.0, card.top() + 14.0);
                    cp.circle_filled(
                        badge_c,
                        9.0,
                        if on {
                            Color32::WHITE
                        } else {
                            Color32::from_rgb(0xec, 0xec, 0xec)
                        },
                    );
                    bold_text(
                        cp,
                        badge_c,
                        &(i + 1).to_string(),
                        FontId::proportional(11.0),
                        if on { CARD_ON } else { INK },
                        0.3,
                    );
                    // 缩略图
                    if let Some(tl) = thumbs_layout.get(i) {
                        let tcolors = confidence_colors(st, len, &partner_map(len, pairs));
                        let box_ = Rect::from_min_max(
                            Pos2::new(card.left() + 6.0, card.top() + 26.0),
                            Pos2::new(card.right() - 6.0, card.bottom() - 26.0),
                        );
                        draw_thumb(cp, tl, &tcolors, box_, on);
                    }
                    // 能量
                    p.text(
                        Pos2::new(card.center().x, card.bottom() - 13.0),
                        Align2::CENTER_CENTER,
                        format!("{energy:.2} kcal/mol"),
                        FontId::proportional(10.8),
                        if on { Color32::WHITE } else { INK },
                    );
                }
            });
        let _ = ui;
    }
    if let Some(i) = clicked_card {
        sel = i;
    }

    // ---- 底部统计条 ----
    {
        let labels = [
            "Free energy of shown structure",
            "Free energy of optimal structure",
            "Free energy of ensemble",
            "Freq. optimal in ensemble",
            "Ensemble diversity",
        ];
        let km = |v: Option<f64>| match v {
            Some(v) => format!("{v:.2} kcal/mol"),
            None => "-- kcal/mol".to_string(),
        };
        let cells = [
            km(Some(structures[sel].1)),
            km(stats_data.map(|s| s.free_energy_optimal)),
            km(stats_data.map(|s| s.free_energy_ensemble)),
            match stats_data {
                Some(s) => format!("{:.2}%", s.freq_optimal * 100.0),
                None => "--%".to_string(),
            },
            match stats_data {
                Some(s) => format!("{:.2}", s.ensemble_diversity),
                None => "--".to_string(),
            },
        ];
        let w = stats.width() / 5.0;
        for (i, (label, value)) in labels.iter().zip(cells.iter()).enumerate() {
            let r = Rect::from_min_size(
                Pos2::new(stats.left() + i as f32 * w, stats.top()),
                Vec2::new(w, stats.height()),
            );
            p.rect_filled(
                r,
                0.0,
                if i == 0 { STATS_FIRST } else { STATS_BG },
            );
            if i > 0 {
                p.line_segment(
                    [r.left_top(), r.left_bottom()],
                    Stroke::new(0.8, Color32::from_rgb(0xca, 0xca, 0xca)),
                );
            }
            let clip = p.with_clip_rect(r.shrink2(Vec2::new(4.0, 0.0)));
            clip.text(
                Pos2::new(r.left() + 10.0, r.top() + 14.0),
                Align2::LEFT_CENTER,
                *label,
                FontId::proportional(12.0),
                Color32::from_rgb(0x22, 0x22, 0x22),
            );
            let vpos = Pos2::new(r.left() + 10.0, r.top() + 38.0);
            let vfont = FontId::proportional(15.0);
            clip.text(
                vpos + Vec2::new(0.4, 0.0),
                Align2::LEFT_CENTER,
                value,
                vfont.clone(),
                Color32::from_rgb(0x11, 0x11, 0x11),
            );
            clip.text(
                vpos,
                Align2::LEFT_CENTER,
                value,
                vfont,
                Color32::from_rgb(0x11, 0x11, 0x11),
            );
        }
    }

    // tooltip（最后一个，避免与其他交互抢）
    if let Some(i) = hover {
        let base = nt(seq.get(i).copied().unwrap_or(b'?')).to_string();
        let conf = base_conf(st, i, partner[i]) * 100.0;
        let mut s = format!("Base {} ({base})", i + 1);
        match partner[i] {
            Some(j) => {
                let nj = nt(seq.get(j).copied().unwrap_or(b'?')).to_string();
                s.push_str(&format!(" · paired with {} ({nj}) · {conf:.0}%", j + 1));
            }
            None => s.push_str(&format!(" · unpaired · {conf:.0}%")),
        }
        egui::Tooltip::for_enabled(&resp)
            .at_pointer()
            .gap(14.0)
            .show(|ui| {
                ui.style_mut().override_font_id = Some(FontId::proportional(11.5));
                ui.label(s);
            });
    }

    // ---- 交回 app ----
    app.struct_zoom = zoom;
    app.struct_pan = pan;
    app.struct_sel = sel;
    app.struct_sel_base = sel_base;
    app.struct_fit = fit;
    rect
}

/// 缩略图：骨架一条折线 + 小圆点（大结构抽稀，避免每帧上千个形状）。
fn draw_thumb(
    p: &egui::Painter,
    layout: &RnaLayout,
    colors: &[Color32],
    area: Rect,
    on_card: bool,
) {
    if layout.pos.is_empty() {
        return;
    }
    let (bmin, bmax) = layout.bounds();
    let bw = (bmax.x - bmin.x).max(1.0);
    let bh = (bmax.y - bmin.y).max(1.0);
    let s = (area.width() / bw).min(area.height() / bh);
    let c0 = Vec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let to = |q: &hx_render::rna_layout::Pt| {
        Pos2::new(
            area.center().x + (q.x - c0.x) * s,
            area.center().y + (q.y - c0.y) * s,
        )
    };
    let pts: Vec<Pos2> = layout.pos.iter().map(to).collect();
    let line = if on_card {
        Color32::from_rgba_unmultiplied(255, 255, 255, 150)
    } else {
        Color32::from_rgb(0x9a, 0x9a, 0x9a)
    };
    p.add(egui::Shape::line(pts.clone(), Stroke::new(0.6, line)));
    let step = (pts.len() / 400).max(1);
    let r = if pts.len() > 400 { 0.9 } else { 1.5 };
    for (i, q) in pts.iter().enumerate() {
        if i % step != 0 {
            continue;
        }
        let c = colors.get(i).copied().unwrap_or(INK);
        p.circle_filled(*q, r, c);
    }
}

/// 画布中央的一行说明（没有结构数据时）。
fn centered_note(ui: &mut egui::Ui, rect: Rect, text: &str) {
    let mut child = ui.new_child(
        egui::UiBuilder::new()
            .max_rect(rect)
            .layout(egui::Layout::top_down(egui::Align::Center)),
    );
    child.add_space(rect.height() * 0.35);
    child.label(
        egui::RichText::new(text)
            .size(13.0)
            .color(theme::TEXT_DIM),
    );
}
