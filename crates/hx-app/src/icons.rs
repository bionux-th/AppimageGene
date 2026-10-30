//! 工具栏/工具面板图标：用 egui 图元**手绘**（不依赖字体里的图标字形，任何平台一致）。

use eframe::egui::{
    self, Align2, Color32, FontId, Pos2, Rect, Response, Sense, Shape, Stroke, StrokeKind, Vec2,
};

use crate::theme;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Icon {
    Save,
    Undo,
    Redo,
    Cut,
    Copy,
    Paste,
    Print,
    Info,
    Sidebar,
    SidebarClosed,
    Settings,
    Home,
    Fit,
    LineWidth,
    ZoomIn,
    ZoomOut,
    RotateCcw,
    RotateCw,
    Ruler,
    Labels,
    Tag,
    Enzymes,
    Search,
    Folder,
    List,
    /// 视图显示开关（SnapGene 的 ViewOptionsToggle）。
    Colors,
    Features,
    Primers,
    Translation,
    /// 图谱拓扑（Show as circular / linear）。
    Circular,
    Linear,
    /// 蛋白 `Show "Region" features`：空心矩形（无方向）。
    Region,
    /// 蛋白 `Show "Site" and "Het" features`：中轴 + 两端下垂的刻度。
    Sites,
    /// 蛋白 `Show "Bond" features`：开口朝上的弧（二硫键）。
    Bonds,
    /// 蛋白图谱 `Show translation numbering`：`(123)`。
    Numbers,
}

/// 顶部工具栏按钮：图标在上、文字在下（SnapGene 样式）。
///
/// oracle 实测（`/tmp/opencode/sg/cmp-toolbar.png`，同一放大幅度并排对照）：
/// SnapGene 的按钮约 54×50、图标框 ≈23 px、图标线宽 ≈2 px、标签文字 cap 高 8 px；
/// 我们原来是 48×40 / 20 px / 1.4 px / 5 px —— 图标和标签都偏小。
pub fn tool(ui: &mut egui::Ui, icon: Icon, label: &str, enabled: bool) -> Response {
    let label=crate::i18n::text(ui.ctx(),label);
    let size = Vec2::new(42.0, 44.0);
    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, resp) = ui.allocate_exact_size(size, sense);
    let hovered = enabled && resp.hovered();
    let bg = if hovered {
        theme::ACCENT_BG
    } else {
        ui.visuals().panel_fill
    };
    let fg = if enabled {
        theme::TEXT
    } else {
        Color32::from_rgb(0xb0, 0xb0, 0xb0)
    };
    if hovered {
        ui.painter().rect_filled(rect, 3.0, bg);
        ui.painter().rect_stroke(
            rect,
            3.0,
            Stroke::new(1.0, theme::ACCENT),
            StrokeKind::Inside,
        );
    }
    let icon_rect = Rect::from_center_size(
        Pos2::new(rect.center().x, rect.top() + 17.0),
        Vec2::splat(21.0),
    );
    paint(ui.painter(), icon, icon_rect, fg, bg);
    ui.painter().text(
        Pos2::new(rect.center().x, rect.bottom() - 9.0),
        Align2::CENTER_CENTER,
        label,
        // oracle 实测：SnapGene 工具栏标签 cap 高 8 px（略大于 FONT_SMALL）
        FontId::proportional(11.5),
        fg,
    );
    resp
}

/// 分体按钮（SnapGene 工具栏的 `Copy`）：主体 + 右下角小三角。
///
/// 主体点下去是 [`Icon`] 的默认动作（Copy = 复制上链），小三角打开下拉菜单。
/// 三角必须用 [`egui::menu::MenuButton`] 承载：`Popup::menu` 里再嵌 `menu_button`
/// 时子菜单拿不到菜单栈（`UiKind::Menu`），悬停能展开但点不中项；`MenuButton`
/// 才有完整的 `MenuState`（菜单栏里的 Edit ▸ Copy 下拉就是这么工作的）。
pub fn tool_split_menu<R>(
    ui: &mut egui::Ui,
    icon: Icon,
    label: &str,
    enabled: bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> (Response, Option<egui::InnerResponse<R>>) {
    let label=crate::i18n::text(ui.ctx(),label);
    let size = Vec2::new(44.0, 44.0);
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let split = rect.right() - 10.0;
    let main_rect = Rect::from_min_max(rect.min, Pos2::new(split, rect.bottom()));
    let arrow_rect = Rect::from_min_max(Pos2::new(split, rect.top()), rect.max);

    let hovered = enabled && ui.rect_contains_pointer(rect);
    let bg = if hovered {
        theme::ACCENT_BG
    } else {
        ui.visuals().panel_fill
    };
    let fg = if enabled {
        theme::TEXT
    } else {
        Color32::from_rgb(0xb0, 0xb0, 0xb0)
    };
    if hovered {
        ui.painter().rect_filled(rect, 3.0, bg);
        ui.painter().rect_stroke(
            rect,
            3.0,
            Stroke::new(1.0, theme::ACCENT),
            StrokeKind::Inside,
        );
    }
    let icon_rect = Rect::from_center_size(
        Pos2::new(main_rect.center().x, rect.top() + 17.0),
        Vec2::splat(21.0),
    );
    paint(ui.painter(), icon, icon_rect, fg, bg);
    ui.painter().text(
        Pos2::new(main_rect.center().x, rect.bottom() - 9.0),
        Align2::CENTER_CENTER,
        label,
        FontId::proportional(11.5),
        fg,
    );
    // 小三角：底边两角 + 顶点，实心
    let c = Pos2::new(rect.right() - 6.0, rect.top() + 17.0);
    ui.painter().add(Shape::convex_polygon(
        vec![
            Pos2::new(c.x - 3.2, c.y - 1.6),
            Pos2::new(c.x + 3.2, c.y - 1.6),
            Pos2::new(c.x, c.y + 1.8),
        ],
        fg,
        Stroke::NONE,
    ));

    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let main = ui.interact(main_rect, ui.scope_id().with(("split", label, 0)), sense);

    // 三角永远可点：菜单里没有可用的复制项时，各项自己置灰（SG 也是这样）
    let mut inner = None;
    ui.scope_builder(
        egui::UiBuilder::new()
            .max_rect(arrow_rect)
            .layout(egui::Layout::centered_and_justified(
                egui::Direction::TopDown,
            )),
        |ui| {
            let btn = egui::Button::new("")
                .frame(false)
                .min_size(arrow_rect.size());
            let response=ui.add(btn);
            let config=egui::menu::MenuConfig::default();
            let popup=egui::Popup::menu(&response).anchor(rect)
                .frame(egui::Frame::new().fill(Color32::WHITE).stroke(Stroke::new(0.8,Color32::from_gray(150))).corner_radius(0).inner_margin(egui::Margin::same(4)))
                .info(egui::UiStackInfo::new(egui::UiKind::Menu).with_tag_value(egui::menu::MenuConfig::MENU_CONFIG_TAG,config))
                .show(|ui|{
                    ui.visuals_mut().window_fill=Color32::WHITE;
                    ui.visuals_mut().window_corner_radius=egui::CornerRadius::ZERO;
                    add_contents(ui)
                });
            inner = popup;
        },
    );
    (main, inner)
}

/// 显示开关按钮（只有图标；`active` = 蓝底白图标）。
pub fn palette(ui: &mut egui::Ui, icon: Icon, active: bool, enabled: bool) -> Response {
    let size = Vec2::splat(28.0);
    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, resp) = ui.allocate_exact_size(size, sense);
    let hovered = enabled && resp.hovered();
    let fg = if !enabled {
        Color32::from_rgb(0xb4, 0xb4, 0xb4)
    } else if active {
        Color32::WHITE
    } else {
        theme::TEXT
    };
    let bg = if active {
        theme::ACCENT
    } else if hovered {
        theme::ACCENT_BG
    } else {
        Color32::WHITE
    };
    if bg != Color32::TRANSPARENT {
        ui.painter().rect_filled(rect, 4.0, bg);
    }
    let icon_rect = Rect::from_center_size(rect.center(), Vec2::splat(20.0));
    paint(ui.painter(), icon, icon_rect, fg, bg);
    resp
}

pub fn navigation(ui:&mut egui::Ui,icon:Icon)->Response {
    let (rect,response)=ui.allocate_exact_size(Vec2::splat(24.0),Sense::click());
    if response.hovered(){ui.painter().rect_filled(rect,2.0,theme::ACCENT_BG);}
    paint(ui.painter(),icon,rect.shrink(3.0),theme::TEXT,Color32::WHITE);response
}

/// Avoid relying on U+22EE support in the system font: the popup trigger and
/// its three dots share a measured hit rectangle.
pub fn more_menu<R>(ui:&mut egui::Ui,contents:impl FnOnce(&mut egui::Ui)->R)
    ->(Response,Option<egui::InnerResponse<R>>) {
    let result=egui::menu::MenuButton::from_button(egui::Button::new("").frame(false).min_size(Vec2::new(24.0,26.0))).ui(ui,contents);
    let c=result.0.rect.center();
    for dy in [-4.0,0.0,4.0] {ui.painter().circle_filled(c+Vec2::new(0.0,dy),1.0,theme::TEXT);}
    result
}

/// New-document `+` trigger. The glyph is painted rather than typeset: a font
/// `+` can sit above the row centre depending on the font. Equal 10 pt arms
/// and a normal 1.5 pt stem keep this cross centred without looking bold.
pub fn plus_menu<R>(ui:&mut egui::Ui,contents:impl FnOnce(&mut egui::Ui)->R)
    ->(Response,Option<egui::InnerResponse<R>>) {
    let result=egui::menu::MenuButton::from_button(
        egui::Button::new("").frame(false)
            .min_size(Vec2::new(30.0,crate::chrome::DOC_TABS_H)))
        .ui(ui,contents);
    let c=result.0.rect.center();
    let arm=5.0;
    let thick=1.5;
    let p=ui.painter();
    p.rect_filled(Rect::from_center_size(c,Vec2::new(arm*2.0,thick)),0.0,theme::TEXT);
    p.rect_filled(Rect::from_center_size(c,Vec2::new(thick,arm*2.0)),0.0,theme::TEXT);
    result
}

/// 画一个图标：`rect` 为图标框（约 20×20），`fg` 前景、`bg` 背景（部分图标需要挖空/遮挡）。
pub fn paint(p: &egui::Painter, icon: Icon, rect: Rect, fg: Color32, bg: Color32) {
    let lp = |x: f32, y: f32| -> Pos2 {
        Pos2::new(
            rect.left() + rect.width() * x / 22.0,
            rect.top() + rect.height() * y / 22.0,
        )
    };
    let u = |v: f32| rect.width() * v / 22.0;
    let boxed = |x0: f32, y0: f32, x1: f32, y1: f32| Rect::from_min_max(lp(x0, y0), lp(x1, y1));
    let s = |w: f32| Stroke::new(u(w), fg);
    match icon {
        Icon::Save => {
            p.rect_stroke(boxed(4.0, 3.0, 18.0, 19.0), 2.0, s(1.5), StrokeKind::Inside);
            p.rect_filled(boxed(7.5, 3.0, 14.5, 8.0), 0.0, fg);
            p.rect_stroke(
                boxed(7.5, 12.0, 14.5, 19.0),
                0.0,
                s(1.3),
                StrokeKind::Inside,
            );
        }
        Icon::Info=>{
            p.circle_stroke(lp(11.0,11.0),u(7.2),s(1.2));
            p.circle_filled(lp(11.0,7.0),u(0.9),fg);
            p.line_segment([lp(11.0,10.0),lp(11.0,15.0)],s(1.3));
            p.line_segment([lp(9.0,15.0),lp(13.0,15.0)],s(1.1));
        }
        Icon::Sidebar|Icon::SidebarClosed=>{
            // Same visual height/centre as Settings and Home. Only the outer
            // panel, its divider and a chevron: no nested box or filled block.
            p.rect_stroke(boxed(2.0,2.0,20.0,20.0),1.0,s(1.3),StrokeKind::Inside);
            p.line_segment([lp(7.0,2.5),lp(7.0,19.5)],s(1.1));
            let pts=if icon==Icon::Sidebar {vec![lp(15.0,6.0),lp(11.0,11.0),lp(15.0,16.0)]}
                else{vec![lp(11.0,6.0),lp(15.0,11.0),lp(11.0,16.0)]};
            p.add(Shape::line(pts,s(1.2)));
        }
        Icon::Home=>{
            p.add(Shape::line(vec![lp(2.5,10.0),lp(11.0,3.0),lp(19.5,10.0)],s(1.2)));
            p.add(Shape::line(vec![lp(5.0,9.0),lp(5.0,19.0),lp(9.0,19.0),lp(9.0,13.0),lp(13.0,13.0),lp(13.0,19.0),lp(17.0,19.0),lp(17.0,9.0)],s(1.2)));
        }
        Icon::Settings=>{
            let points=(0..32).map(|i|{let a=i as f32*std::f32::consts::TAU/32.0;let r=if i%4<2 {8.7}else{6.5};lp(11.0+a.sin()*r,11.0+a.cos()*r)}).collect();
            p.add(Shape::closed_line(points,s(1.15)));p.circle_stroke(lp(11.0,11.0),u(3.0),s(1.1));
        }
        Icon::Undo | Icon::RotateCcw => {
            let c = lp(11.5, 12.5);
            let r = u(6.0);
            let pts = arc(c, r, 200.0, 340.0, 18);
            p.add(Shape::line(pts.clone(), s(1.5)));
            // 箭头（200° 端，指向切线方向）
            let at=pts[0];p.add(Shape::line(vec![at+Vec2::new(u(0.6),-u(3.5)),at,at+Vec2::new(u(3.5),u(0.6))],s(1.3)));
        }
        Icon::Redo | Icon::RotateCw => {
            let c = lp(10.5, 12.5);
            let r = u(6.0);
            let pts = arc(c, r, 200.0, 340.0, 18);
            // 镜像：先镜像点再画
            let m = |q: Pos2| Pos2::new(2.0 * c.x - q.x, q.y);
            p.add(Shape::line(pts.iter().copied().map(m).collect(), s(1.5)));
            let at=pts[0];p.add(Shape::line(vec![m(at+Vec2::new(u(0.6),-u(3.5))),m(at),m(at+Vec2::new(u(3.5),u(0.6)))],s(1.3)));
        }
        Icon::Cut => {
            p.circle_stroke(lp(6.0, 5.5), u(2.1), s(1.3));
            p.circle_stroke(lp(6.0, 16.5), u(2.1), s(1.3));
            p.line_segment([lp(8.0, 6.6), lp(18.0, 16.0)], s(1.3));
            p.line_segment([lp(8.0, 15.4), lp(18.0, 6.0)], s(1.3));
        }
        Icon::Copy => {
            p.rect_stroke(boxed(3.5, 3.5, 13.0, 15.5), 1.5, s(1.3), StrokeKind::Inside);
            let front = boxed(8.0, 7.5, 18.5, 19.0);
            p.rect_filled(front, 1.5, bg);
            p.rect_stroke(front, 1.5, s(1.3), StrokeKind::Inside);
        }
        Icon::Paste => {
            p.rect_stroke(boxed(5.0, 4.5, 17.0, 20.0), 1.5, s(1.4), StrokeKind::Inside);
            p.rect_filled(boxed(8.5, 2.5, 13.5, 6.5), 1.0, fg);
            for i in 0..3 {
                let y = 10.0 + i as f32 * 3.0;
                p.line_segment([lp(8.0, y), lp(14.0, y)], s(1.1));
            }
        }
        Icon::Print => {
            p.rect_stroke(boxed(7.0, 2.5, 15.0, 7.0), 0.5, s(1.2), StrokeKind::Inside);
            p.rect_stroke(boxed(3.0, 7.0, 19.0, 15.0), 1.0, s(1.4), StrokeKind::Inside);
            let paper = boxed(6.0, 13.0, 16.0, 20.0);
            p.rect_filled(paper, 0.0, bg);
            p.rect_stroke(paper, 0.0, s(1.3), StrokeKind::Inside);
        }
        Icon::Fit => {
            p.rect_stroke(boxed(6.0,4.0,16.0,18.0),1.5,s(1.2),StrokeKind::Inside);
            p.line_segment([lp(2.0,11.0),lp(20.0,11.0)],s(1.0));
            p.line_segment([lp(2.0,9.0),lp(2.0,13.0)],s(1.0));
            p.line_segment([lp(20.0,9.0),lp(20.0,13.0)],s(1.0));
        }
        Icon::LineWidth => {
            // 固定行宽：两端边界、水平双向箭头和下拉提示，全部由图元绘制。
            for x in [4.0, 17.0] {
                p.line_segment([lp(x, 3.0), lp(x, 14.0)], s(1.3));
            }
            p.line_segment([lp(4.5, 8.5), lp(16.5, 8.5)], s(1.3));
            p.add(Shape::line(
                vec![lp(7.5, 5.5), lp(4.5, 8.5), lp(7.5, 11.5)],
                s(1.3),
            ));
            p.add(Shape::line(
                vec![lp(13.5, 5.5), lp(16.5, 8.5), lp(13.5, 11.5)],
                s(1.3),
            ));
            p.add(Shape::convex_polygon(
                vec![lp(8.0, 17.0), lp(13.0, 17.0), lp(10.5, 20.0)],
                fg,
                Stroke::NONE,
            ));
        }
        Icon::ZoomIn | Icon::ZoomOut | Icon::Search => {
            p.circle_stroke(lp(9.5, 9.5), u(5.0), s(1.4));
            p.line_segment([lp(13.3, 13.3), lp(19.0, 19.0)], s(1.6));
            match icon {
                Icon::ZoomIn => {
                    p.line_segment([lp(6.5, 9.5), lp(12.5, 9.5)], s(1.3));
                    p.line_segment([lp(9.5, 6.5), lp(9.5, 12.5)], s(1.3));
                }
                Icon::ZoomOut => {
                    p.line_segment([lp(6.5, 9.5), lp(12.5, 9.5)], s(1.3));
                }
                _ => {}
            }
        }
        Icon::Ruler => {
            // 横向标尺：条 + 长短刻度
            p.rect_stroke(boxed(2.5, 7.5, 19.5, 14.5), 1.0, s(1.3), StrokeKind::Inside);
            let mut i = 0;
            while i < 5 {
                let x = 5.0 + i as f32 * 3.5;
                let y1 = if i % 2 == 0 { 12.2 } else { 10.6 };
                p.line_segment([lp(x, 8.1), lp(x, y1)], s(1.1));
                i += 1;
            }
        }
        Icon::Tag => {
            // 标签牌（左侧小孔 + 横线）
            let pts = vec![
                lp(3.5, 7.5),
                lp(15.0, 7.5),
                lp(19.5, 11.5),
                lp(15.0, 15.5),
                lp(3.5, 15.5),
            ];
            p.add(Shape::convex_polygon(pts, bg, s(1.3)));
            p.circle_filled(lp(6.4, 11.5), u(1.1), fg);
            p.line_segment([lp(9.2, 11.5), lp(14.2, 11.5)], s(1.1));
        }
        Icon::Colors => {
            // 画笔：斜向笔杆、笔颈与刷毛。
            p.add(Shape::closed_line(
                vec![lp(9.0, 12.0), lp(17.0, 3.0), lp(19.0, 5.0), lp(11.0, 14.0)],
                s(1.3),
            ));
            p.line_segment([lp(12.0, 9.0), lp(15.0, 12.0)], s(1.1));
            p.add(Shape::closed_line(
                vec![
                    lp(9.0, 12.0),
                    lp(11.0, 14.0),
                    lp(9.0, 18.0),
                    lp(3.0, 19.0),
                    lp(6.0, 16.0),
                    lp(6.0, 14.0),
                ],
                s(1.2),
            ));
        }
        Icon::Features => {
            // 单体空心特征箭头，与图谱中的特征符号一致。
            p.add(Shape::closed_line(
                vec![
                    lp(2.5, 5.5),
                    lp(14.0, 5.5),
                    lp(20.0, 11.0),
                    lp(14.0, 16.5),
                    lp(2.5, 16.5),
                ],
                s(1.4),
            ));
        }
        Icon::Primers => {
            // 上下链 + 相反方向的两条短引物。
            p.line_segment([lp(2.5, 9.0), lp(19.5, 9.0)], s(1.1));
            p.line_segment([lp(2.5, 13.0), lp(19.5, 13.0)], s(1.1));
            p.line_segment([lp(7.0, 5.0), lp(18.5, 5.0)], s(1.3));
            p.add(Shape::line(
                vec![lp(15.5, 2.5), lp(18.5, 5.0), lp(15.5, 7.5)],
                s(1.3),
            ));
            p.line_segment([lp(3.5, 17.0), lp(15.0, 17.0)], s(1.3));
            p.add(Shape::line(
                vec![lp(6.5, 14.5), lp(3.5, 17.0), lp(6.5, 19.5)],
                s(1.3),
            ));
        }
        Icon::Translation => {
            // 两条相反方向的阅读框箭头。
            p.line_segment([lp(3.0, 6.0), lp(19.0, 6.0)], s(1.4));
            p.add(Shape::line(
                vec![lp(15.0, 2.5), lp(19.0, 6.0), lp(15.0, 9.5)],
                s(1.4),
            ));
            p.line_segment([lp(3.0, 16.0), lp(19.0, 16.0)], s(1.4));
            p.add(Shape::line(
                vec![lp(7.0, 12.5), lp(3.0, 16.0), lp(7.0, 19.5)],
                s(1.4),
            ));
        }
        Icon::Circular => {
            p.circle_stroke(lp(11.0, 11.0), u(6.5), s(1.4));
            p.circle_stroke(lp(11.0,11.0),u(4.6),s(0.8));
        }
        Icon::Region => {
            // 空心矩形（oracle：SG 蛋白视图的第一颗按钮就是一个无方向矩形）。
            p.rect_stroke(boxed(3.5, 6.0, 18.5, 16.0), 0.0, s(1.4), StrokeKind::Inside);
        }
        Icon::Sites => {
            // 中轴竖线 + 横杆 + 两端下垂的腿（oracle：SG 的 Site 刻度符号）。
            p.line_segment([lp(11.0, 3.5), lp(11.0, 12.0)], s(1.4));
            p.line_segment([lp(5.5, 12.0), lp(16.5, 12.0)], s(1.4));
            p.add(Shape::line(
                vec![lp(5.5, 12.0), lp(5.5, 17.0), lp(3.0, 17.0)],
                s(1.4),
            ));
            p.add(Shape::line(
                vec![lp(16.5, 12.0), lp(16.5, 17.0), lp(19.0, 17.0)],
                s(1.4),
            ));
        }
        Icon::Bonds => {
            // 开口朝上的半圆（oracle：SG 的 Bond 特征图标）。
            let pts: Vec<Pos2> = (0..=18)
                .map(|i| {
                    let t = std::f32::consts::PI * i as f32 / 18.0;
                    lp(11.0 - 6.5 * t.cos(), 5.0 + 6.5 * t.sin())
                })
                .collect();
            p.add(Shape::line(pts, s(1.4)));
        }
        Icon::Numbers => {
            p.text(
                lp(11.0, 11.0),
                Align2::CENTER_CENTER,
                "(123)",
                FontId::proportional(u(8.5)),
                fg,
            );
        }
        Icon::Linear => {
            p.line_segment([lp(3.5, 9.5), lp(18.5, 9.5)], s(1.1));
            p.line_segment([lp(3.5, 12.5), lp(18.5, 12.5)], s(0.8));
            for x in [5.0,17.0] {p.line_segment([lp(x,8.0),lp(x,14.0)],s(1.0));}
        }
        Icon::Enzymes => {
            // 酶名及切点标记，避免与工具栏中的 Cut 剪刀混淆。
            p.text(
                lp(11.0, 7.5),
                Align2::CENTER_CENTER,
                "SalI",
                FontId::proportional(u(8.0)),
                fg,
            );
            p.add(Shape::line(
                vec![lp(6.0, 12.0), lp(6.0, 17.0), lp(17.0, 17.0), lp(17.0, 21.0)],
                s(1.3),
            ));
        }
        Icon::Labels => {
            for i in 0..3 {
                let y = 6.0 + i as f32 * 5.0;
                p.rect_filled(boxed(3.5, y - 1.4, 7.0, y + 1.4), 0.0, fg);
                p.line_segment([lp(8.5, y), lp(18.5, y)], s(1.2));
            }
        }
        Icon::Folder => {
            let pts = vec![
                lp(2.5, 5.0),
                lp(9.0, 5.0),
                lp(11.0, 7.5),
                lp(19.5, 7.5),
                lp(19.5, 18.5),
                lp(2.5, 18.5),
            ];
            p.add(Shape::convex_polygon(pts, bg, s(1.4)));
        }
        Icon::List => {
            for i in 0..3 {
                let y = 6.0 + i as f32 * 5.0;
                p.circle_filled(lp(4.2, y), u(1.1), fg);
                p.line_segment([lp(7.0, y), lp(18.0, y)], s(1.2));
            }
        }
    }
}

/// 圆弧采样点（角度单位为度，屏幕坐标顺时针：270° = 正上方）。
fn arc(c: Pos2, r: f32, a0: f32, a1: f32, n: usize) -> Vec<Pos2> {
    (0..=n)
        .map(|i| {
            let a = (a0 + (a1 - a0) * i as f32 / n as f32).to_radians();
            Pos2::new(c.x + r * a.cos(), c.y + r * a.sin())
        })
        .collect()
}
