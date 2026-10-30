//! 配色与颜色工具：对齐 SnapGene 的**浅色**界面（浅灰面板 + 白画布 + 近黑文字 + 蓝色选中）。
//!
//! 特征/酶的颜色仍来自文件（`#rrggbb`）；这里只提供界面骨架色与"在浅底上可读"的调整工具。

use eframe::egui::{Color32, Stroke};

/// 解析 `#rrggbb` / `rrggbb`；失败返回 `fallback`。
pub fn hex(s: &str, fallback: Color32) -> Color32 {
    let t = s.trim().trim_start_matches('#');
    if t.len() < 6 || !t[..6].bytes().all(|b| b.is_ascii_hexdigit()) {
        return fallback;
    }
    let v = u32::from_str_radix(&t[..6], 16).unwrap_or(0);
    Color32::from_rgb((v >> 16) as u8, (v >> 8) as u8, v as u8)
}

/// 线性混合（`t` = 0 → `a`，`t` = 1 → `b`）。
pub fn mix(a: Color32, b: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let f = |x: u8, y: u8| -> u8 {
        (x as f32 + (y as f32 - x as f32) * t)
            .round()
            .clamp(0.0, 255.0) as u8
    };
    Color32::from_rgb(f(a.r(), b.r()), f(a.g(), b.g()), f(a.b(), b.b()))
}

/// 向黑色压暗。
pub fn darken(c: Color32, t: f32) -> Color32 {
    mix(c, Color32::BLACK, t)
}

/// 感知亮度（0..1）。
pub fn luma(c: Color32) -> f32 {
    (0.299 * c.r() as f32 + 0.587 * c.g() as f32 + 0.114 * c.b() as f32) / 255.0
}

/// 在浅色（白）底上保证可读：过亮的颜色按亮度压暗（黄色特征 → 深橄榄/褐）。
pub fn readable_on_white(c: Color32) -> Color32 {
    let l = luma(c);
    if l > 0.62 {
        darken(c, (l - 0.62) * 1.7)
    } else {
        c
    }
}

/// 特征带填充（SnapGene 风格：实色但略淡，便于压在带上的文字辨认）。
pub fn band_fill(c: Color32) -> Color32 {
    mix(c, Color32::WHITE, 0.12)
}

pub fn with_alpha(c: Color32, a: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a)
}

/// 碱基颜色（A 绿 / C 蓝 / G 橙 / T 红 / 其它灰）——仅在用户开启"按碱基着色"时使用。
pub fn base_letter_color(b: u8) -> Color32 {
    match b.to_ascii_uppercase() {
        b'A' => Color32::from_rgb(0x1f, 0x7a, 0x1f),
        b'C' => Color32::from_rgb(0x1b, 0x4f, 0xb8),
        b'G' => Color32::from_rgb(0xc8, 0x74, 0x00),
        b'T' | b'U' => Color32::from_rgb(0xc0, 0x28, 0x28),
        _ => Color32::from_rgb(0x70, 0x70, 0x70),
    }
}

// ---- 界面骨架色（浅色，对齐 SnapGene）----
/// 画布（图谱）底色。
pub const CANVAS: Color32 = Color32::WHITE;
/// 序列视图（DNA/RNA/蛋白的 `Sequence` 页）底色。
///
/// oracle 实测：SnapGene 的**图谱**页是纯白，而**序列**页整块是 `#F2F2F2`
/// （`/tmp/opencode/oracle/out/sg_dna_seq_c.png`、`/tmp/opencode/sg/or2-prot-seq2.png`、
/// `i0-rna-seq.png` 三张在同一位置取样都是 242）。
pub const SEQ_BG: Color32 = Color32::from_rgb(0xf2, 0xf2, 0xf2);
/// 面板底色（工具栏/侧栏/状态栏）。
pub const PANEL: Color32 = Color32::from_rgb(0xf0, 0xf0, 0xf0);
/// 面板上的浅色条（页签底、输入框底）。
pub const PANEL_ALT: Color32 = Color32::from_rgb(0xf7, 0xf7, 0xf7);
/// 分隔线/边框。
pub const BORDER: Color32 = Color32::from_rgb(0xd0, 0xd0, 0xd0);
/// 主文字。
pub const TEXT: Color32 = Color32::from_rgb(0x1d, 0x1d, 0x1d);
/// 次要文字（尺寸/坐标/提示）。
pub const TEXT_DIM: Color32 = Color32::from_rgb(0x63, 0x69, 0x71);
/// SnapGene 蓝（选中/链接/活动页签）。
pub const ACCENT: Color32 = Color32::from_rgb(0x1a, 0x73, 0xe8);
/// 活动页签/按钮的淡蓝底。
pub const ACCENT_BG: Color32 = Color32::from_rgb(0xcf, 0xe2, 0xfa);
/// 列表选中行底色。
pub const ROW_SEL: Color32 = Color32::from_rgb(0xcc, 0xdd, 0xf5);
/// 斑马纹。
pub const ROW_ALT: Color32 = Color32::from_rgb(0xf7, 0xf8, 0xfa);

/// 蛋白 `Properties` 页的隔行底色（真 SnapGene 实测 `#FFEFC3`）。
pub const STRIPE: Color32 = Color32::from_rgb(0xff, 0xef, 0xc3);
/// 表格表头下的细分隔线（真 SnapGene 实测 `#7F7F7F`）。
pub const RULE: Color32 = Color32::from_rgb(0x7f, 0x7f, 0x7f);

/// 骨架（环/轴）：近黑。
pub const BACKBONE: Color32 = Color32::from_rgb(0x1a, 0x1a, 0x1a);
/// 骨架细线（互补链）。
pub const BACKBONE_THIN: Color32 = Color32::from_rgb(0x9a, 0x9a, 0x9a);
/// 刻度线/尺线。
pub const RULER: Color32 = Color32::from_rgb(0xb8, 0xbc, 0xc2);
/// 刻度数字（蓝灰，SnapGene 序列视图刻度是蓝灰）。
pub const RULER_TEXT: Color32 = Color32::from_rgb(0x46, 0x5a, 0x78);
/// 主刻度。
pub const TICK_MAJOR: Color32 = Color32::from_rgb(0x5a, 0x60, 0x68);
/// 次刻度。
pub const TICK_MINOR: Color32 = Color32::from_rgb(0xb4, 0xb8, 0xbe);
/// 酶位点（非唯一）。
pub const SITE: Color32 = Color32::from_rgb(0x62, 0x62, 0x62);
/// 酶位点（唯一，SnapGene 用粗体）。
pub const SITE_UNIQUE: Color32 = Color32::from_rgb(0x24, 0x24, 0x24);
/// 无颜色特征的回退色。
pub const FEATURE_FALLBACK: Color32 = Color32::from_rgb(0xa8, 0xac, 0xb4);

/// `Show colors`（左侧竖列的画笔按钮）关掉时特征条的填充色。
///
/// oracle 实测：SnapGene 关掉画笔后，图谱里所有特征条（含 `Region` 小条）都画成
/// 这个近黑色 `(37,37,37)`，条内文字转白；`FEATURE_FALLBACK` 是**配色解析失败**时
/// 的兜底灰，两者不是一回事。
pub const FEATURE_NO_COLOR: Color32 = Color32::from_rgb(0x25, 0x25, 0x25);
/// 选中（高亮/光标/选段）。
pub const SELECTION: Color32 = Color32::from_rgb(0x2f, 0x7f, 0xd8);
/// 序列光标。
pub const CURSOR: Color32 = Color32::from_rgb(0x1a, 0x73, 0xe8);
/// 标签白色晕圈底（压在有色带上的文字）。
pub const HALO: Color32 = Color32::WHITE;
/// 序列视图行底。
pub const SEQ_LINE: Color32 = Color32::from_rgb(0x1a, 0x1a, 0x1a);

/// 链下单碱基刻度线（序列视图里链与特征条之间那条深色刻度线）。
pub const TICK_DARK: Color32 = Color32::from_rgb(0x3c, 0x3c, 0x3c);

/// 特征条续行指示（`⟩`）的灰色。
pub const CHEVRON: Color32 = Color32::from_rgb(0x8a, 0x8a, 0x8a);

/// 特征条描边的压暗目标色。
pub const BLACK_EDGE: Color32 = Color32::from_rgb(0x20, 0x20, 0x20);
/// 单链序列视图里特征条的描边色（oracle 实测：与条底色无关的中性深灰
/// `#434343`，RNA 的绿色条也是同一个灰）。
pub const BAR_EDGE: Color32 = Color32::from_rgb(0x43, 0x43, 0x43);
/// 单链序列视图的标尺线灰（oracle：蛋白 `#7F7F7F`、RNA `#666666`，取后者偏中）。
pub const SEQ_RULER: Color32 = Color32::from_rgb(0x6e, 0x6e, 0x6e);
/// 位点轴（蛋白 `type="Site"`）的线色：近黑。
pub const SITE_AXIS: Color32 = Color32::from_rgb(0x1a, 0x1a, 0x1a);
/// 鼠标悬停在位点括线上时的线色（oracle 实测 `(253,1,1)`）。
pub const SITE_HOVER: Color32 = Color32::from_rgb(0xfd, 0x01, 0x01);
/// 悬停时位点标签的底色 / 描边（oracle 实测 `(230,244,255)` / `(0,142,255)`）。
pub const SITE_HOVER_BG: Color32 = Color32::from_rgb(0xe6, 0xf4, 0xff);
pub const SITE_HOVER_EDGE: Color32 = Color32::from_rgb(0x00, 0x8e, 0xff);
/// 位点标签连回自己括线的细引线（oracle 实测未选中时 `#7F7F7F`）。
pub const SITE_LEADER: Color32 = Color32::from_rgb(0x7f, 0x7f, 0x7f);
/// 选中特征在字母带里的底（oracle 实测 `#B5D5FF`：SG 点中特征后其残基整段盖这个蓝底）。
/// 引物（SnapGene 的引物紫，序列视图 / 图谱 / 面板统一用它）。
pub const PRIMER: Color32 = Color32::from_rgb(0xa0, 0x20, 0xf0);

pub const FEAT_SEL_BG: Color32 = Color32::from_rgb(0xb5, 0xd5, 0xff);
/// 选中特征的标签底色（oracle 实测 `#3875D7`：标签变蓝底白字）。
pub const FEAT_SEL_EDGE: Color32 = Color32::from_rgb(0x38, 0x75, 0xd7);
/// 互补链碱基。
pub const BASE_COMPLEMENT: Color32 = Color32::from_rgb(0x55, 0x55, 0x55);

/// egui 视觉样式（浅色）。
/// 分子类型小图标（SnapGene 的页签/侧栏图标：圆角方框 + 类型字形）。
///
/// 实测：**活动项**是实心色块 + 白色字形，**非活动项**是同一形状的描边版；
/// DNA 画双螺旋（两条交叉折线）、RNA 画单链（一条折线）、蛋白画链状（波浪线）。
/// `rect` 是图标方框（约 13×13）。
pub fn kind_icon(
    p: &eframe::egui::Painter,
    rect: eframe::egui::Rect,
    kind: hx_core::MoleculeKind,
    active: bool,
) {
    let _=active;
    crate::branding::document_icon(p,rect,kind,false);
}

#[allow(dead_code)]
fn legacy_kind_icon(p:&eframe::egui::Painter,rect:eframe::egui::Rect,kind:hx_core::MoleculeKind,active:bool){
    use eframe::egui::{Align2, CornerRadius, Pos2, Stroke};
    let (r, g, b) = kind.icon_rgb();
    let color = Color32::from_rgb(r, g, b);
    let text_color = if active {
        Color32::WHITE
    } else {
        color
    };
    let radius = CornerRadius::same((rect.width() * 0.24) as u8);
    if active {
        // 实心底 + 稍亮的边（活动页签的图标更醒目）
        p.rect_filled(rect, radius, color);
    } else {
        p.rect_stroke(
            rect,
            radius,
            Stroke::new(1.0, with_alpha(color, 170)),
            eframe::egui::StrokeKind::Inside,
        );
    }
    // 字形区域：图标框内部（留 22% 边距）
    let inner = rect.shrink(rect.width() * 0.22);
    let (x0, x1) = (inner.left(), inner.right());
    let (y0, y1) = (inner.top(), inner.bottom());
    let mid_y = inner.center().y;
    let stroke = Stroke::new((inner.width() * 0.14).max(0.9), text_color);
    match kind {
        hx_core::MoleculeKind::Dna => {
            // 双螺旋：两条反相正弦（SnapGene 的图标就是小方框里一条双螺旋；
            // 之前的两段折线在小尺寸下会读成 `≫`）
            let steps = 14;
            let amp = inner.width() * 0.30;
            for phase in [0.0f32, std::f32::consts::PI] {
                let pts: Vec<Pos2> = (0..=steps)
                    .map(|i| {
                        let t = i as f32 / steps as f32;
                        let y = y0 + t * (y1 - y0);
                        let x = inner.center().x + amp * (t * 9.42 + phase).sin();
                        Pos2::new(x, y)
                    })
                    .collect();
                p.add(eframe::egui::Shape::line(
                    pts,
                    Stroke::new(stroke.width * 0.9, text_color),
                ));
            }
        }
        hx_core::MoleculeKind::Rna => {
            // 单链：一条正弦
            let steps = 14;
            let amp = inner.width() * 0.30;
            let pts: Vec<Pos2> = (0..=steps)
                .map(|i| {
                    let t = i as f32 / steps as f32;
                    let y = y0 + t * (y1 - y0);
                    let x = inner.center().x + amp * (t * std::f32::consts::TAU).sin();
                    Pos2::new(x, y)
                })
                .collect();
            p.add(eframe::egui::Shape::line(
                pts,
                Stroke::new(stroke.width * 0.9, text_color),
            ));
        }
        hx_core::MoleculeKind::Protein => {
            // 链：三个小圆点 + 连线
            let dots = [
                Pos2::new(x0, y1),
                Pos2::new(inner.center().x, mid_y),
                Pos2::new(x1, y0),
            ];
            p.add(eframe::egui::Shape::line(
                vec![dots[0], dots[1], dots[2]],
                Stroke::new(stroke.width * 0.85, text_color),
            ));
            for d in dots {
                p.circle_filled(d, (inner.width() * 0.16).max(0.9), text_color);
            }
        }
    }
    let _ = Align2::CENTER_CENTER;
}

pub fn visuals() -> eframe::egui::Visuals {
    use eframe::egui::style::Widgets;
    let mut v = eframe::egui::Visuals::light();
    v.panel_fill = PANEL;
    v.window_fill = PANEL_ALT;
    v.window_stroke = Stroke::new(1.0, BORDER);
    v.extreme_bg_color = Color32::WHITE;
    v.faint_bg_color = ROW_ALT;
    v.override_text_color = Some(TEXT);
    v.hyperlink_color = ACCENT;
    v.selection.bg_fill = ACCENT_BG;
    v.selection.stroke = Stroke::new(1.0, ACCENT);
    v.widgets = Widgets {
        noninteractive: eframe::egui::style::WidgetVisuals {
            bg_fill: PANEL,
            weak_bg_fill: PANEL,
            bg_stroke: Stroke::new(1.0, BORDER),
            fg_stroke: Stroke::new(1.0, TEXT),
            corner_radius: 3.into(),
            expansion: 0.0,
        },
        inactive: eframe::egui::style::WidgetVisuals {
            bg_fill: Color32::from_rgb(0xf2, 0xf2, 0xf2),
            weak_bg_fill: Color32::from_rgb(0xf2, 0xf2, 0xf2),
            bg_stroke: Stroke::new(1.0, Color32::from_rgb(0xc4, 0xc4, 0xc4)),
            fg_stroke: Stroke::new(1.0, TEXT),
            corner_radius: 3.into(),
            expansion: 0.0,
        },
        hovered: eframe::egui::style::WidgetVisuals {
            bg_fill: ACCENT_BG,
            weak_bg_fill: ACCENT_BG,
            bg_stroke: Stroke::new(1.0, ACCENT),
            fg_stroke: Stroke::new(1.0, TEXT),
            corner_radius: 3.into(),
            expansion: 0.0,
        },
        active: eframe::egui::style::WidgetVisuals {
            bg_fill: mix(ACCENT_BG, ACCENT, 0.25),
            weak_bg_fill: mix(ACCENT_BG, ACCENT, 0.25),
            bg_stroke: Stroke::new(1.0, ACCENT),
            fg_stroke: Stroke::new(1.0, TEXT),
            corner_radius: 3.into(),
            expansion: 0.0,
        },
        open: eframe::egui::style::WidgetVisuals {
            bg_fill: Color32::from_rgb(0xe8, 0xe8, 0xe8),
            weak_bg_fill: Color32::from_rgb(0xe8, 0xe8, 0xe8),
            bg_stroke: Stroke::new(1.0, BORDER),
            fg_stroke: Stroke::new(1.0, TEXT),
            corner_radius: 3.into(),
            expansion: 0.0,
        },
    };
    v
}

/// UI 正文字号（菜单栏、工具栏标签、正文、面板、对话框）。
///
/// oracle 实测（同一 Xvfb、1:1，`/tmp/opencode/sg/ink.py` 量文字墨迹高度）：
/// SnapGene 菜单/页签/序列字母的 cap 高都是 **10 px**，我们换字体前是 9 / 7 / 7。
/// Noto Sans 的 cap/em ≈ 0.714 → 10 px cap 对应字号 ≈ 14 px。
pub const FONT_UI: f32 = 14.0;
/// 视图页签（SnapGene 页签文字 cap 10 px，与正文同号）。
pub const FONT_TAB: f32 = 14.0;
/// 次要说明文字（面板里的小标题、单位、注释）。
pub const FONT_SMALL: f32 = 12.0;
/// 极小字号（上标、刻度标注等）。
pub const FONT_MICRO: f32 = 11.0;

/// 把浅色样式装到上下文（`eframe::CreationContext` 里调用一次）。
pub fn apply(ctx: &eframe::egui::Context) {
    // 字体先装：SnapGene 用系统 sans（本机 Noto Sans），与 egui 自带的 Ubuntu-Light
    // 在字高/字重上差别明显（见 `crate::fonts`）。
    crate::fonts::install(ctx);
    let visuals = visuals();
    ctx.all_styles_mut(|s| {
        use eframe::egui::{FontFamily, FontId, TextStyle};

        s.visuals = visuals.clone();
        // 字号统一走 `TextStyle`，这样 egui 自带的控件（菜单、按钮、输入框、
        // 对话框）会跟着一起变，不必逐处改。
        s.text_styles.insert(
            TextStyle::Body,
            FontId::new(FONT_UI, FontFamily::Proportional),
        );
        s.text_styles.insert(
            TextStyle::Button,
            FontId::new(FONT_UI, FontFamily::Proportional),
        );
        s.text_styles.insert(
            TextStyle::Small,
            FontId::new(FONT_SMALL, FontFamily::Proportional),
        );
        s.text_styles.insert(
            TextStyle::Heading,
            FontId::new(FONT_UI + 4.0, FontFamily::Proportional),
        );
        s.spacing.item_spacing = eframe::egui::vec2(8.0, 5.0);
        s.spacing.button_padding = eframe::egui::vec2(10.0, 5.0);
        s.spacing.menu_margin = eframe::egui::Margin::same(8);
        // SnapGene 的菜单栏行高 ≈30 px、工具栏按钮 ≈52 px（实测），字号变大后
        // 命中区也要跟着长，否则控件会挤在一起。
        s.spacing.interact_size.y = 24.0;
    });
}
