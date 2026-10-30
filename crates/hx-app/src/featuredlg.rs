//! 「Add Feature… / Edit Feature…」对话框（SnapGene 8.x 的特征面板）。
//!
//! 英文文案取自真 SnapGene 8.2.2 二进制字符串（`sg-strings.txt`），
//! 字段与排布依据用户提供的 8.2.2 截图（见 `docs/new-file-and-features-2026-10-18.zh-CN.md`
//! 图 4/图 5）：
//!
//! ```text
//!  ┌ Add Feature for Untitled ───────────────────────────────────────────┐
//!  │ Feature Name:  [_______________]                                    │
//!  │ Type:  [misc_feature ▾]  [→][←][↔][·]                               │
//!  │ [ ] Translate this feature in Sequence view   [Feature Translation Options...] │
//!  │ ┌ Segments ───────────────────────────────────────────────────────┐ │
//!  │ │ 28 base pairs · 1 segment    [Split Feature…][Merge Segments][Delete Segment] │
//!  │ │ #  Name        Location     Size     Color                      │ │
//!  │ │ 1  misc_feature 28..55       28 bp    ■                          │ │
//!  │ └─────────────────────────────────────────────────────────────────┘ │
//!  │ [ ] Prioritize display of this feature in maps                      │
//!  │ ┌ Qualifiers ─────────────────────────────────────────────────────┐ │
//!  │ │ [ /note ▾ ] [+] [−]                                            │ │
//!  │ │ [ free text                                                    ] │ │
//!  │ └─────────────────────────────────────────────────────────────────┘ │
//!  │                                          [Cancel]  [OK]             │
//!  └─────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! 位置一律是 **1-based 闭区间**（与 `block 10` 的 `<Segment range="A-B">` 同构）；
//! `A > B` 表示环状分子上跨原点。

use eframe::egui::{self, RichText};
use hx_formats::feature_edit::{
    default_color, FeatSpec, FeatureEdit, FEATURE_TYPES, NO_COLOR, PALETTE,
};

use crate::theme;

/// 对话框模式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Add,
    Edit,
}

/// 对话框结果。
pub enum Outcome {
    /// 用户按下 OK（携带要施加的编辑）。
    Apply(Box<FeatureEdit>),
    /// 用户按下 Cancel / Esc / 关闭。
    Cancel,
}

/// 「Add/Edit Feature」对话框状态。存在即表示对话框打开。
pub struct FeatureDlg {
    pub mode: Mode,
    /// Edit 模式：被编辑特征在文档里的下标（Features 面板行号）。
    pub index: Option<usize>,
    /// 文档名（标题 `Add Feature for %1`）。
    pub doc_name: String,
    /// 序列长度（片段大小与跨原点片段）。
    pub seq_len: usize,
    /// 工作模型。
    pub spec: FeatSpec,
    /// 打开时的模型（判断是否有实际改动）。
    pub orig: FeatSpec,
    /// 每个片段的「位置」输入缓冲。
    loc_buf: Vec<String>,
    /// 片段表格里选中的行（Merge/Delete Segment 的作用对象）。
    pub sel_seg: usize,
    /// 颜色调色板展开的片段下标。
    pub color_popup: Option<usize>,
    /// 「限定符」下拉的选项：`0` = note，其余 = `other_quals[i-1]`。
    pub sel_qual: usize,
    /// 「添加限定符 ▾」菜单是否展开。
    pub add_qual_menu: bool,
    /// 错误提示（位置格式等）。
    pub err: Option<String>,
}

impl FeatureDlg {
    /// Add 模式：`start..=end` 为 1-based 闭区间。
    pub fn add(doc_name: &str, seq_len: usize, start: usize, end: usize) -> Self {
        let mut spec = FeatSpec::new(start.max(1), end.max(start));
        // SnapGene 的新特征名默认 = 类型名
        spec.name = spec.kind.clone();
        let orig = spec.clone();
        Self {
            mode: Mode::Add,
            index: None,
            doc_name: doc_name.to_string(),
            seq_len,
            spec,
            orig,
            loc_buf: Vec::new(),
            sel_seg: 0,
            color_popup: None,
            sel_qual: 0,
            add_qual_menu: false,
            err: None,
        }
        .synced()
    }

    /// Edit 模式：编辑已有特征。
    pub fn edit(doc_name: &str, seq_len: usize, index: usize, spec: FeatSpec) -> Self {
        Self {
            mode: Mode::Edit,
            index: Some(index),
            doc_name: doc_name.to_string(),
            seq_len,
            orig: spec.clone(),
            spec,
            loc_buf: Vec::new(),
            sel_seg: 0,
            color_popup: None,
            sel_qual: 0,
            add_qual_menu: false,
            err: None,
        }
        .synced()
    }

    /// Add 模式（带初始模型）：`Duplicate Feature…` 用它复制一个已有特征。
    pub fn add_spec(doc_name: &str, seq_len: usize, spec: FeatSpec) -> Self {
        Self {
            mode: Mode::Add,
            index: None,
            doc_name: doc_name.to_string(),
            seq_len,
            orig: spec.clone(),
            spec,
            loc_buf: Vec::new(),
            sel_seg: 0,
            color_popup: None,
            sel_qual: 0,
            add_qual_menu: false,
            err: None,
        }
        .synced()
    }

    /// 同步「位置」输入缓冲。
    fn synced(mut self) -> Self {
        self.sync_loc();
        self
    }

    fn sync_loc(&mut self) {
        self.loc_buf = self
            .spec
            .segments
            .iter()
            .map(|s| format!("{}..{}", s.start, s.end))
            .collect();
        if self.sel_seg >= self.spec.segments.len() {
            self.sel_seg = self.spec.segments.len().saturating_sub(1);
        }
    }

    /// 标题：`Add Feature for %1` / `Edit Feature in %1`（SnapGene 原文）。
    pub fn title(&self) -> String {
        match self.mode {
            Mode::Add => format!("Add Feature for {}", self.doc_name),
            Mode::Edit => format!("Edit Feature in {}", self.doc_name),
        }
    }

    /// 是否有实际改动。
    pub fn changed(&self) -> bool {
        self.spec != self.orig
    }

    /// 生成要施加的编辑（OK 时调用）。
    pub fn build_edit(&self) -> FeatureEdit {
        match self.mode {
            Mode::Add => FeatureEdit::Add(self.spec.clone()),
            Mode::Edit => FeatureEdit::Update(self.spec.clone()),
        }
    }

    /// 片段总长（bp / aa）。
    pub fn total_len(&self) -> usize {
        self.spec.span_len(self.seq_len)
    }
}

/// 「位置」文本解析：`a..b` / `a-b` / `a b` / 单个数（= a..a）。1-based 闭区间。
pub fn parse_location(s: &str) -> Option<(usize, usize)> {
    let t = s.trim();
    if t.is_empty() {
        return None;
    }
    let norm = t.replace("..", "-").replace(' ', "-").replace(',', "");
    let parts: Vec<&str> = norm.split('-').filter(|p| !p.is_empty()).collect();
    match parts.len() {
        1 => {
            let a: usize = parts[0].parse().ok()?;
            Some((a, a))
        }
        2 => {
            let a: usize = parts[0].parse().ok()?;
            let b: usize = parts[1].parse().ok()?;
            Some((a, b))
        }
        _ => None,
    }
}

/// 片段长度（1-based 闭；`a > b` 走环状换算）。
fn seg_len(a: usize, b: usize, seq_len: usize) -> usize {
    if a <= b {
        b - a + 1
    } else {
        seq_len.saturating_sub(a) + 1 + b.min(seq_len)
    }
}

// ---------------------------------------------------------------------------
// UI
// ---------------------------------------------------------------------------

const DLG_W: f32 = 620.0;

/// 画对话框；返回 `Some(..)` 表示用户做了决定。
pub fn draw(dlg: &mut FeatureDlg, ctx: &egui::Context) -> Option<Outcome> {
    let mut out: Option<Outcome> = None;
    let mut open = true;
    egui::Window::new("feature_dialog")
        .id(egui::Id::unique("feature_dialog"))
        .title_bar(false)
        .resizable(false)
        .collapsible(false)
        .open(&mut open)
        .fixed_size(egui::vec2(DLG_W, 0.0))
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -30.0))
        .frame(
            egui::Frame::window(&ctx.style_of(egui::Theme::Light))
                .fill(egui::Color32::WHITE)
                .stroke(egui::Stroke::new(0.8, theme::BORDER))
                .inner_margin(egui::Margin::symmetric(12, 10)),
        )
        .show(ctx, |ui| {
            ui.set_width(DLG_W);
            // ---- 标题条 ----
            ui.horizontal(|ui| {
                let template=match dlg.mode{Mode::Add=>"Add Feature for {}",Mode::Edit=>"Edit Feature in {}"};
                ui.label(RichText::new(crate::i18n::format_values(ui.ctx(),template,&[dlg.doc_name.clone()])).size(13.5).strong());
            });
            ui.separator();
            ui.add_space(4.0);

            // ---- Feature Name ----
            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ui.ctx(),"Feature Name:"));
                ui.add_space(4.0);
                ui.add(
                    egui::TextEdit::singleline(&mut dlg.spec.name).desired_width(DLG_W - 200.0),
                );
            });
            ui.add_space(4.0);

            // ---- Type + 链向箭头 ----
            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ui.ctx(),"Type:"));
                ui.add_space(28.0);
                let mut picked: Option<String> = None;
                egui::ComboBox::from_id_salt("feat_type")
                    .selected_text(dlg.spec.kind.clone())
                    .width(220.0)
                    .show_ui(ui, |ui| {
                        egui::ScrollArea::vertical()
                            .max_height(280.0)
                            .id_salt("feat_type_scroll")
                            .show(ui, |ui| {
                                for t in FEATURE_TYPES {
                                    if ui.selectable_label(*t == dlg.spec.kind, *t).clicked() {
                                        picked = Some((*t).to_string());
                                    }
                                }
                            });
                    });
                if let Some(t) = picked {
                    dlg.spec.kind = t.clone();
                    // SnapGene 实测：名字还等于旧类型名时随类型联动
                    if dlg.spec.name == dlg.orig.kind || dlg.spec.name.is_empty() {
                        dlg.spec.name = t.clone();
                    }
                    // 该类型的内建默认色
                    let c = default_color(&t).to_string();
                    for s in &mut dlg.spec.segments {
                        s.color = c.clone();
                    }
                }
                ui.add_space(6.0);
                // 链向：4 个小钮（SnapGene 的 4 态箭头；用画笔画三角，避免缺字形）
                for (val, kind, tip) in [
                    (1u8, 0u8, "Forward"),
                    (2, 1, "Reverse"),
                    (3, 2, "Both"),
                    (0, 3, "No direction"),
                ] {
                    let on = dlg.spec.directionality == val as i32;
                    let (rect, resp) =
                        ui.allocate_exact_size(egui::vec2(26.0, 22.0), egui::Sense::click());
                    let p = ui.painter();
                    let bg = if on { theme::ACCENT } else { egui::Color32::WHITE };
                    p.rect_filled(rect, 3.0, bg);
                    p.rect_stroke(
                        rect.shrink(0.5),
                        3.0,
                        egui::Stroke::new(1.0, theme::BORDER),
                        egui::StrokeKind::Inside,
                    );
                    let ink = if on { egui::Color32::WHITE } else { theme::TEXT };
                    let c = rect.center();
                    let tri = |p: &egui::Painter, x: f32, dir: f32| {
                        p.add(egui::Shape::convex_polygon(
                            vec![
                                egui::pos2(x - 4.5 * dir, c.y - 5.0),
                                egui::pos2(x - 4.5 * dir, c.y + 5.0),
                                egui::pos2(x + 4.5 * dir, c.y),
                            ],
                            ink,
                            egui::Stroke::NONE,
                        ));
                    };
                    match kind {
                        0 => tri(p, c.x, 1.0),
                        1 => tri(p, c.x, -1.0),
                        2 => {
                            tri(p, c.x - 2.5, 1.0);
                            tri(p, c.x + 2.5, -1.0);
                        }
                        _ => {
                            p.line_segment(
                                [egui::pos2(c.x - 5.0, c.y), egui::pos2(c.x + 5.0, c.y)],
                                egui::Stroke::new(1.6, ink),
                            );
                        }
                    }
                    if resp.on_hover_text(crate::i18n::text(ui.ctx(),tip)).clicked() {
                        dlg.spec.directionality = val as i32;
                    }
                }
            });
            ui.add_space(4.0);

            // ---- Translate ----
            ui.horizontal(|ui| {
                // 复选框是**特征级**开关：勾选/取消会把所有片段的 `translated` 一起设定
                let mut tr = dlg.spec.translate;
                if ui
                    .checkbox(&mut tr, crate::i18n::text(ui.ctx(),"Translate this feature in Sequence view"))
                    .changed()
                {
                    dlg.spec.set_translate(tr);
                }
                ui.add_space(10.0);
                ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Feature Translation Options...")))
                    .on_hover_text(crate::i18n::text(ui.ctx(),"Translation options (genetic code, frame) are not implemented yet"));
            });
            ui.add_space(6.0);

            // ---- Segments ----
            egui::Frame::new()
                .stroke(egui::Stroke::new(1.0, theme::BORDER))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(8, 6))
                .show(ui, |ui| {
                    ui.set_width(DLG_W - 26.0);
                    ui.horizontal(|ui| {
                        let unit = "bp";
                        ui.label(
                            RichText::new(format!("{} {unit}", crate::app::nfmt(dlg.total_len())))
                                .size(11.5),
                        );
                        ui.label(RichText::new("·").size(11.5).color(theme::TEXT_DIM));
                        let n = dlg.spec.segments.len();
                        ui.label(
                            RichText::new(if crate::i18n::Language::get(ui.ctx())==crate::i18n::Language::English{format!("{n} segment{}",if n==1{""}else{"s"})}else{crate::ui_format!(ui.ctx(),"{} SEGMENTS",n)})
                                .size(11.5),
                        );
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            let del = ui
                                .add_enabled(
                                    dlg.spec.segments.len() > 1,
                                    egui::Button::new(crate::i18n::text(ui.ctx(),"Delete Segment")),
                                )
                                .clicked();
                            let merge = ui
                                .add_enabled(
                                    dlg.spec.segments.len() > 1,
                                    egui::Button::new(crate::i18n::text(ui.ctx(),"Merge Segments")),
                                )
                                .clicked();
                            let split = ui
                                .add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Split Feature...")))
                                .on_hover_text(crate::i18n::text(ui.ctx(),"Splitting a feature needs a chosen split point (not implemented yet)"));
                            let _ = split;
                            if del {
                                dlg.err = None;
                                let i = dlg.sel_seg.min(dlg.spec.segments.len() - 1);
                                dlg.spec.segments.remove(i);
                                dlg.sync_loc();
                            }
                            if merge {
                                dlg.err = None;
                                let (i, j) = (dlg.sel_seg, dlg.sel_seg + 1);
                                if j < dlg.spec.segments.len() {
                                    let a = dlg.spec.segments[i].clone();
                                    let b = dlg.spec.segments[j].clone();
                                    // 只有相邻/重叠才能合并
                                    let adjacent = if a.start <= a.end {
                                        a.end + 1 >= b.start || b.start <= a.end
                                    } else {
                                        false
                                    };
                                    if adjacent {
                                        // 合并后只要有一段原本是翻译片段，就保持翻译状态
                                        dlg.spec.segments[i].translated |= b.translated;
                                        dlg.spec.segments[i].end =
                                            dlg.spec.segments[i].end.max(b.end);
                                        dlg.spec.segments.remove(j);
                                        dlg.sync_loc();
                                    } else {
                                        dlg.err = Some(
                                            "Only adjacent or overlapping segments can be merged"
                                                .into(),
                                        );
                                    }
                                }
                            }
                        });
                    });
                    ui.add_space(4.0);
                    // 片段表： # / Name / Location / Size / Color
                    let mut new_sel = dlg.sel_seg;
                    ui.horizontal(|ui| {
                        for (w, h) in [
                            (28.0, "#"),
                            (150.0, "Name"),
                            (130.0, "Location"),
                            (72.0, "Size"),
                            (0.0, "Color"),
                        ] {
                            if w > 0.0 {
                                ui.allocate_ui_with_layout(
                                    egui::vec2(w, 16.0),
                                    egui::Layout::left_to_right(egui::Align::Center),
                                    |ui| {
                                        ui.label(
                                            RichText::new(h).size(11.0).color(theme::TEXT_DIM),
                                        );
                                    },
                                );
                            } else {
                                ui.label(RichText::new(h).size(11.0).color(theme::TEXT_DIM));
                            }
                        }
                    });
                    // 逐行
                    let n = dlg.spec.segments.len();
                    for i in 0..n {
                        ui.horizontal(|ui| {
                            let sel = i == dlg.sel_seg;
                            let resp = ui
                                .allocate_ui_with_layout(
                                    egui::vec2(28.0, 22.0),
                                    egui::Layout::left_to_right(egui::Align::Center),
                                    |ui| ui.selectable_label(sel, format!("{}", i + 1)),
                                )
                                .inner;
                            if resp.clicked() {
                                new_sel = i;
                            }
                            ui.allocate_ui_with_layout(
                                egui::vec2(150.0, 22.0),
                                egui::Layout::left_to_right(egui::Align::Center),
                                |ui| {
                                    ui.label(
                                        RichText::new(if dlg.spec.name.is_empty() {
                                            dlg.spec.kind.clone()
                                        } else {
                                            dlg.spec.name.clone()
                                        })
                                        .size(11.5),
                                    );
                                },
                            );
                            // 位置（可编辑）
                            let loc = ui.add_sized(
                                egui::vec2(130.0, 22.0),
                                egui::TextEdit::singleline(&mut dlg.loc_buf[i])
                                    .font(egui::TextStyle::Monospace),
                            );
                            if loc.gained_focus() || loc.changed() || loc.lost_focus() {
                                new_sel = i;
                                match parse_location(&dlg.loc_buf[i]) {
                                    Some((a, b)) if a >= 1 => {
                                        dlg.spec.segments[i].start = a;
                                        dlg.spec.segments[i].end = b;
                                        dlg.err = None;
                                    }
                                    _ => {
                                        dlg.err = Some(format!(
                                            "Invalid location \u{201c}{}\u{201d} — use start..end (1-based)",
                                            dlg.loc_buf[i]
                                        ));
                                    }
                                }
                            }
                            let seg = dlg.spec.segments[i].clone();
                            ui.allocate_ui_with_layout(
                                egui::vec2(72.0, 22.0),
                                egui::Layout::left_to_right(egui::Align::Center),
                                |ui| {
                                    ui.label(
                                        RichText::new(format!(
                                            "{} bp",
                                            crate::app::nfmt(seg_len(
                                                seg.start, seg.end, dlg.seq_len
                                            ))
                                        ))
                                        .size(11.0)
                                        .color(theme::TEXT_DIM),
                                    );
                                },
                            );
                            // 颜色
                            let (rect, cresp) =
                                ui.allocate_exact_size(egui::vec2(40.0, 20.0), egui::Sense::click());
                            let p = ui.painter();
                            let sw = egui::Rect::from_center_size(rect.center(), egui::vec2(26.0, 14.0));
                            if seg.color == NO_COLOR || seg.color.is_empty() {
                                p.rect_filled(sw, 2.0, egui::Color32::WHITE);
                                p.rect_stroke(
                                    sw,
                                    2.0,
                                    egui::Stroke::new(1.0, theme::BORDER),
                                    egui::StrokeKind::Inside,
                                );
                                p.line_segment(
                                    [sw.left_bottom(), sw.right_top()],
                                    egui::Stroke::new(1.4, egui::Color32::from_rgb(0xd0, 0, 0)),
                                );
                            } else {
                                let c = theme::hex(&seg.color, theme::FEATURE_FALLBACK);
                                p.rect_filled(sw, 2.0, c);
                                p.rect_stroke(
                                    sw,
                                    2.0,
                                    egui::Stroke::new(1.0, theme::BORDER),
                                    egui::StrokeKind::Inside,
                                );
                            }
                            if cresp.clicked() {
                                dlg.color_popup = if dlg.color_popup == Some(i) { None } else { Some(i) };
                            }
                        });
                    }
                    if let Some(i) = dlg.color_popup {
                        if i < dlg.spec.segments.len() {
                            if let Some(c) = palette_row(ui, dlg.spec.segments[i].color.as_str()) {
                                dlg.spec.segments[i].color = c;
                                dlg.color_popup = None;
                            }
                        }
                    }
                    dlg.sel_seg = new_sel.min(n.saturating_sub(1));
                    if let Some(e) = &dlg.err {
                        ui.add_space(2.0);
                        ui.label(RichText::new(e.clone()).size(11.0).color(egui::Color32::from_rgb(0xb0, 0x20, 0x20)));
                    }
                });
            ui.add_space(4.0);

            // ---- Prioritize ----
            ui.add(egui::Checkbox::new(
                &mut dlg.spec.prioritize,
                "Prioritize display of this feature in maps",
            ));
            ui.add_space(6.0);

            // ---- Qualifiers ----
            egui::Frame::new()
                .stroke(egui::Stroke::new(1.0, theme::BORDER))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(8, 6))
                .show(ui, |ui| {
                    ui.set_width(DLG_W - 26.0);
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Qualifiers")).size(11.0).color(theme::TEXT_DIM));
                        ui.add_space(8.0);
                        // 下拉：note + 其它限定符
                        let cur = if dlg.sel_qual == 0 {
                            "/note".to_string()
                        } else {
                            format!(
                                "/{}",
                                dlg.spec
                                    .other_quals
                                    .get(dlg.sel_qual - 1)
                                    .map(|(n, _)| n.as_str())
                                    .unwrap_or("?")
                            )
                        };
                        let mut pick: Option<usize> = None;
                        egui::ComboBox::from_id_salt("feat_qual")
                            .selected_text(cur)
                            .width(150.0)
                            .show_ui(ui, |ui| {
                                if ui.selectable_label(dlg.sel_qual == 0, "/note").clicked() {
                                    pick = Some(0);
                                }
                                for (i, (n, _)) in dlg.spec.other_quals.iter().enumerate() {
                                    if ui
                                        .selectable_label(dlg.sel_qual == i + 1, format!("/{n}"))
                                        .clicked()
                                    {
                                        pick = Some(i + 1);
                                    }
                                }
                            });
                        if let Some(p) = pick {
                            dlg.sel_qual = p;
                        }
                        ui.menu_button("+", |ui| {
                            for name in [
                                "label",
                                "gene",
                                "product",
                                "function",
                                "EC_number",
                                "db_xref",
                                "note",
                            ] {
                                if ui.button(format!("/{name}")).clicked() {
                                    dlg.spec.other_quals.push((name.to_string(), String::new()));
                                    dlg.sel_qual = dlg.spec.other_quals.len();
                                    ui.close();
                                }
                            }
                        })
                        .response
                        .on_hover_text(crate::i18n::text(ui.ctx(),"Add qualifier"));
                        if ui
                            .add_enabled(
                                dlg.sel_qual > 0,
                                egui::Button::new("−"),
                            )
                            .on_hover_text(crate::i18n::text(ui.ctx(),"Remove qualifier"))
                            .clicked()
                        {
                            let i = dlg.sel_qual - 1;
                            if i < dlg.spec.other_quals.len() {
                                dlg.spec.other_quals.remove(i);
                            }
                            dlg.sel_qual = 0;
                        }
                    });
                    ui.add_space(3.0);
                    // 值编辑框
                    if dlg.sel_qual == 0 {
                        ui.add(
                            egui::TextEdit::multiline(&mut dlg.spec.note)
                                .desired_width(DLG_W - 46.0)
                                .desired_rows(3),
                        );
                    } else {
                        let i = dlg.sel_qual - 1;
                        if i < dlg.spec.other_quals.len() {
                            let mut v = dlg.spec.other_quals[i].1.clone();
                            if ui
                                .add(
                                    egui::TextEdit::multiline(&mut v)
                                        .desired_width(DLG_W - 46.0)
                                        .desired_rows(3),
                                )
                                .changed()
                            {
                                dlg.spec.other_quals[i].1 = v;
                            }
                        }
                    }
                    // 其它限定符一览（只读提示）
                    if !dlg.spec.other_quals.is_empty() {
                        ui.add_space(3.0);
                        let list = dlg
                            .spec
                            .other_quals
                            .iter()
                            .map(|(n, v)| {
                                format!(
                                    "/{n}={}",
                                    if v.is_empty() { "…" } else { v.as_str() }
                                )
                            })
                            .collect::<Vec<_>>()
                            .join("   ");
                        ui.label(RichText::new(list).size(10.5).color(theme::TEXT_DIM));
                    }
                });
            ui.add_space(10.0);

            // ---- OK / Cancel ----
            ui.horizontal(|ui| {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let ok = ui.add(
                        egui::Button::new(RichText::new(crate::i18n::text(ui.ctx(),"OK")).color(egui::Color32::WHITE))
                            .fill(theme::ACCENT)
                            .min_size(egui::vec2(80.0, 24.0)),
                    );
                    if ok.clicked() {
                        out = Some(Outcome::Apply(Box::new(dlg.build_edit())));
                    }
                    if ui
                        .add(egui::Button::new(crate::i18n::text(ui.ctx(),"Cancel")).min_size(egui::vec2(80.0, 24.0)))
                        .clicked()
                    {
                        out = Some(Outcome::Cancel);
                    }
                });
            });
            if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
                out = Some(Outcome::Cancel);
            }
        });
    if !open {
        out = Some(Outcome::Cancel);
    }
    // 片段编辑（合并 / 删除 / 改区间）后把特征级「翻译」开关重新对齐片段
    dlg.spec.translate = dlg.spec.segments.iter().any(|s| s.translated);
    out
}

/// 颜色调色板一行（SnapGene 的 6 列网格 + No Color）；返回被选中的颜色。
fn palette_row(ui: &mut egui::Ui, cur: &str) -> Option<String> {
    let mut picked = None;
    ui.add_space(2.0);
    egui::Frame::new()
        .fill(egui::Color32::from_rgb(0xf7, 0xf7, 0xf7))
        .stroke(egui::Stroke::new(1.0, theme::BORDER))
        .corner_radius(3.0)
        .inner_margin(egui::Margin::symmetric(6, 6))
        .show(ui, |ui| {
            for chunk in PALETTE.chunks(6) {
                ui.horizontal(|ui| {
                    for hex in chunk {
                        let c = theme::hex(hex, egui::Color32::GRAY);
                        let (rect, resp) =
                            ui.allocate_exact_size(egui::vec2(22.0, 18.0), egui::Sense::click());
                        let sw = egui::Rect::from_center_size(rect.center(), egui::vec2(20.0, 16.0));
                        ui.painter().rect_filled(sw, 2.0, c);
                        let w = if cur.eq_ignore_ascii_case(hex) { 2.0 } else { 1.0 };
                        ui.painter().rect_stroke(
                            sw,
                            2.0,
                            egui::Stroke::new(w, theme::BORDER),
                            egui::StrokeKind::Inside,
                        );
                        if resp.clicked() {
                            picked = Some((*hex).to_string());
                        }
                    }
                });
            }
            ui.horizontal(|ui| {
                let c = theme::hex(NO_COLOR, egui::Color32::WHITE);
                let (rect, resp) =
                    ui.allocate_exact_size(egui::vec2(60.0, 18.0), egui::Sense::click());
                let sw = egui::Rect::from_center_size(rect.center(), egui::vec2(58.0, 16.0));
                ui.painter().rect_filled(sw, 2.0, c);
                ui.painter().rect_stroke(
                    sw,
                    2.0,
                    egui::Stroke::new(1.0, theme::BORDER),
                    egui::StrokeKind::Inside,
                );
                ui.painter().line_segment(
                    [sw.left_bottom(), sw.right_top()],
                    egui::Stroke::new(1.4, egui::Color32::from_rgb(0xd0, 0, 0)),
                );
                if resp.clicked() {
                    picked = Some(NO_COLOR.to_string());
                }
                ui.label(RichText::new(crate::i18n::text(ui.ctx(),"No Color")).size(11.0));
            });
        });
    picked
}
