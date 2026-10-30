//! 主区列表视图（SnapGene 的 Enzymes / Features / Primers / History 页签）。
//!
//! 这些页签在 SnapGene 里是**主区**文档视图（不是侧栏），所以放在这里而不是 `sidebar.rs`。

use eframe::egui::{self, Align2, Color32};

use crate::app::{nfmt, nfmt_plain, EnzymePreset, HelixApp, Tab};
use crate::theme;

/// 数据页共用的内容边距；顶部让出悬浮页签，再留出一段阅读间距。
pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) {
    let rect = ui.available_rect_before_wrap();
    // 数据页底色是 SnapGene 的浅灰（实测 #F2F2F2），不是画布的白色
    ui.painter().rect_filled(rect, 0.0, theme::PANEL);
    let content = egui::Rect::from_min_max(
        rect.min + egui::vec2(24.0, 58.0),
        rect.max - egui::vec2(24.0, 16.0),
    );
    ui.scope_builder(egui::UiBuilder::new().max_rect(content), |ui| {
        ui.set_clip_rect(ui.clip_rect().intersect(content));
        ui.spacing_mut().item_spacing.y = 8.0;
        match app.tab {
            Tab::Enzymes => enzymes(ui, app),
            Tab::Features => features(ui, app),
            Tab::Primers => primers(ui, app),
            Tab::History => history(ui, app),
            _ => {}
        }
    });
}

/// **不能编辑的原生 SnapGene 文档**（`DNA alignment` / `Protein alignment` / `Gel`）的概览页。
///
/// 这类文件与 `.dna` 用同一套 TLV 容器（我们的解析器能完整读、能逐字节写回），
/// 但内容属于 SnapGene 的比对编辑器 / 凝胶编辑器。页面如实说明：
/// 类别、块结构、保存行为，以及为什么没有序列视图。
pub fn draw_native_overview(ui: &mut egui::Ui, app: &mut HelixApp) {
    let rect = ui.available_rect_before_wrap();
    let content = egui::Rect::from_min_max(
        rect.min + egui::vec2(24.0, 58.0),
        rect.max - egui::vec2(24.0, 16.0),
    );
    let Some(doc) = app.doc.as_ref() else {
        return;
    };
    let Some(n) = doc.native_only.as_ref() else {
        return;
    };
    let name = doc
        .path
        .as_ref()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| doc.display_name());
    let class = n.class;
    let tos = n.type_of_sequence;
    let size = n.bytes.len();
    let blocks = n.blocks.clone();
    if tos==8{
        let bytes=n.bytes.clone();
        let mut available=Vec::new();
        for t in &app.tabs{if let Some(d)=t.doc.as_ref().filter(|d|d.is_dna()&&d.native_class().is_none()){if let Ok((bytes,_))=d.session.save_bytes_with_history(&d.map_name(),&d.map_name(),None){available.push(crate::gel_view::Source{name:d.display_name(),bytes});}}}
        for path in app.recent.iter().filter(|p|p.extension().is_some_and(|s|s.eq_ignore_ascii_case("dna"))).take(12){
            if let Ok(d)=crate::doc::Document::from_path(path){if d.is_dna()&&d.native_class().is_none(){let bytes=d.session.bytes().to_vec();if !available.iter().any(|s|s.bytes==bytes){available.push(crate::gel_view::Source{name:format!("{} (recent)",d.display_name()),bytes});}}}
        }
        if app.gel_view.as_ref().is_none_or(|g|g.key!=bytes){
            match crate::gel_view::GelView::parse(&bytes){Ok(mut g)=>{if let Some(old)=&app.gel_view{g.selected=old.selected.min(g.lanes.len()-1);}app.gel_view=Some(Box::new(g));},Err(e)=>{ui.colored_label(Color32::RED,e);return;}}
        }
        ui.scope_builder(egui::UiBuilder::new().max_rect(content),|ui|{if let Some(g)=app.gel_view.as_mut(){g.available=available;g.show(ui);}});
        if let Some(g)=app.gel_view.as_mut(){if g.changed{match g.serialize(){Ok(bytes)=>{if let Some(d)=app.doc.as_mut(){match d.edit_native(bytes.clone()){Ok(())=>{g.key=bytes;g.changed=false;},Err(e)=>g.error=e}}},Err(e)=>g.error=e}}}
        return;
    }
    if matches!(tos,3|4|8){
        ui.scope_builder(egui::UiBuilder::new().max_rect(content),|ui|{if let Err(e)=crate::native_view::show(ui,&n.bytes){ui.colored_label(Color32::RED,e);}});
        return;
    }
    ui.scope_builder(egui::UiBuilder::new().max_rect(content), |ui| {
        ui.set_clip_rect(ui.clip_rect().intersect(content));
        ui.spacing_mut().item_spacing.y = 8.0;
        ui.heading(class);
        ui.label(egui::RichText::new(name).size(12.0).color(theme::TEXT_DIM));
        ui.add_space(6.0);
        ui.label(
            egui::RichText::new(format!(
                "This is a SnapGene {class} document (header type_of_sequence = {tos}).\n\
                 AppimageGene reads and writes the container faithfully, but the \
                 {class} editor itself is not implemented, so there is no sequence \
                 or map view for it.",
            ))
            .size(12.5),
        );
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new(
                crate::i18n::text(ui.ctx(),"Saving writes the original bytes back unchanged (no edits are applied)."),
            )
            .size(12.0)
            .color(theme::TEXT_DIM),
        );
        ui.add_space(10.0);
        ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"Container")).strong());
        egui::Grid::new("native_blocks")
            .num_columns(3)
            .striped(true)
            .spacing([14.0, 4.0])
            .show(ui, |ui| {
                ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"Field")).strong());
                ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"Value")).strong());
                ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"Notes")).strong());
                ui.end_row();
                ui.label(crate::i18n::text(ui.ctx(),"File size"));
                ui.label(crate::ui_format!(ui.ctx(),"{} bytes", nfmt_plain(size)));
                ui.label("");
                ui.end_row();
                ui.label(crate::i18n::text(ui.ctx(),"Blocks"));
                ui.label(format!("{}", blocks.len()));
                ui.label(crate::i18n::text(ui.ctx(),"type · payload bytes"));
                ui.end_row();
                for (i, (ty, len)) in blocks.iter().enumerate() {
                    ui.label(format!("#{i}"));
                    ui.label(format!("{ty} · {len}"));
                    ui.label(block_note(*ty));
                    ui.end_row();
                }
            });
    });
}

/// 常见块类型的一句话说明（仅用于概览页）。
fn block_note(ty: u8) -> &'static str {
    match ty {
        9 => "file header (SnapGene magic + type/version)",
        0 => "DNA sequence",
        21 => "protein sequence",
        32 => "RNA sequence",
        16 => "trace / multi-sequence container",
        10 => "features XML",
        5 => "primers XML",
        6 => "notes XML",
        3 => "enzyme site cache",
        13 => "enzyme display preferences",
        24 => "alignment settings (JSON)",
        23 => "attachments",
        33 => "gel summary",
        _ => "",
    }
}

const COLUMN_GAP: f32 = 12.0;

/// 属性的单元格样式：普通文字，或消光系数（`M⁻¹ cm⁻¹` 用手工上标）。
#[derive(Clone, Copy, PartialEq, Eq)]
enum CellKind {
    Text,
    Extinction,
}

/// **蛋白文档的 `Properties` 页签**（SnapGene 的 `Properties and Amino Acid Composition`）。
///
/// 版式按 oracle 实测（真 SnapGene 8.2.2，`/tmp/opencode/oracle/out/kpyk1_props*.png`、
/// `x2_props.png`、`stop_props.png`、`kpyk1_props_sel.png`）：
///
/// ```text
///                                                       Properties   ← 页名在页签行右端
///                       Whole Protein     Selection
/// Length                    470 aa                18 aa
/// ▸ Molecular Weight     50,729.56 Da        1,996.35 Da   ← 展开：average / monoisotopic
/// ▸ Extinction Coefficient (280 nm)  7450 M⁻¹ cm⁻¹  0 M⁻¹ cm⁻¹ ← 展开：fully reduced / fully oxidized
/// ▸ Absorbance (280 nm, 0.1%)   0.15           0.00          ← 展开：fully reduced / fully oxidized
/// Isoelectric Point (pI)    5.86               4.10
/// Charge at pH [7.0 ▾]     -8.52              -2.77
///   Amino Acid     Number  Percent     Number  Percent   ← 两对列**总是**存在（无选中留白）
/// A  Ala  Alanine    41      8.72          1     5.56      ← 偶数行浅黄底（实测 #FFEFC3）
/// ```
///
/// 非标准残基的实测行为：
/// - 未知残基（`X`/`B`/`Z`/`J`/`U`/`O` 等）：数值照算，**分子量前加 `~`**，
///   组成表不单列它们（百分比仍按全长算）。
/// - 终止符 `*`：`Length` 的单位变成 `codons`，其余数值一律显示 `—`。
///
/// 数字来自 [`hx_core::protein`]：分子量、消光系数与吸光度与 SnapGene 逐位对齐；
/// pI 与净电荷用标准 ExPASy 模型，与 SnapGene 有零点几的差别（见该模块文档）。
pub fn draw_properties(ui: &mut egui::Ui, app: &mut HelixApp) {
    use hx_core::protein::Composition;

    /// 属性行 / 展开明细行 / 氨基酸表行的行高（oracle 实测 22.5 / 24 / 21）。
    const ROW_H: f32 = 22.5;
    const SUB_H: f32 = 24.0;
    const AA_H: f32 = 21.0;
    /// 数值不可用时的占位（SnapGene 用长破折号）。
    const DASH: &str = "—";

    let rect = ui.available_rect_before_wrap();
    // 数据页底色是 SnapGene 的浅灰（实测 #F2F2F2），不是画布的白色
    ui.painter().rect_filled(rect, 0.0, theme::PANEL);
    // 页名：与悬浮页签同一行的右端
    ui.painter().text(
        egui::pos2(rect.right() - 16.0, rect.top() + 20.0),
        Align2::RIGHT_CENTER,
        crate::i18n::text(ui.ctx(),"Properties"),
        egui::FontId::proportional(15.0),
        theme::TEXT,
    );
    let content = egui::Rect::from_min_max(
        rect.min + egui::vec2(24.0, 58.0),
        rect.max - egui::vec2(16.0, 12.0),
    );
    if content.width() < 320.0 || content.height() < 140.0 {
        return;
    }
    // 取数据（`Composition` 是 `Copy`，取完就放开对 `app` 的借用——后面还要改它的状态）
    let taken = app.doc.as_ref().map(|doc| {
        let seq = doc.sequence();
        let whole = Composition::of(seq);
        let sel = app
            .selection
            .filter(|(s, e)| e > s && *s < seq.len() && *e <= seq.len())
            .map(|(s, e)| Composition::of(&seq[s..e]));
        (whole, sel, doc.unit())
    });
    let Some((whole, sel, unit)) = taken else {
        return;
    };

    let x = content.left();
    let val_r = x + 282.0;
    let sel_r = x + 494.0;
    let font = egui::FontId::proportional(12.5);
    let sub_font = egui::FontId::proportional(12.0);
    let p = ui.painter().clone();

    // ---- 列表头 ----
    let mut y = content.top() + 8.0;
    p.text(
        egui::pos2(x + 251.0, y),
        Align2::CENTER_CENTER,
        crate::i18n::text(ui.ctx(),"Whole Protein"),
        font.clone(),
        theme::TEXT,
    );
    p.text(
        egui::pos2(x + 402.0, y),
        Align2::LEFT_CENTER,
        crate::i18n::text(ui.ctx(),"Selection"),
        font.clone(),
        theme::TEXT,
    );
    y += 27.0;

    // ---- 属性行（`*` → 一律 `—`；未知残基让分子量带 `~`）----
    let len_of = |c: &Composition| {
        if c.has_stops() {
            format!("{} codons", nfmt_plain(c.len()))
        } else {
            format!("{} {}", nfmt_plain(c.len()), unit)
        }
    };
    let mass_of = |c: &Composition| {
        if c.has_stops() {
            DASH.to_string()
        } else if c.has_unknown() {
            format!("~{} Da", mass_str(hx_core::protein::average_mass(c)))
        } else {
            format!("{} Da", mass_str(hx_core::protein::average_mass(c)))
        }
    };
    let mono_of = |c: &Composition| {
        if c.has_stops() {
            DASH.to_string()
        } else if c.has_unknown() {
            format!("~{} Da", mass_str(hx_core::protein::monoisotopic_mass(c)))
        } else {
            format!("{} Da", mass_str(hx_core::protein::monoisotopic_mass(c)))
        }
    };
    let ext_of = |c: &Composition, oxidized: bool| {
        if c.has_stops() {
            DASH.to_string()
        } else {
            let (r, o) = hx_core::protein::extinction_280(c);
            nfmt(if oxidized { o } else { r } as usize)
        }
    };
    let abs_of = |c: &Composition, oxidized: bool| {
        if c.has_stops() {
            DASH.to_string()
        } else {
            let (r, o) = hx_core::protein::absorbance_280(c);
            format!("{:.2}", if oxidized { o } else { r })
        }
    };
    let pi_of = |c: &Composition| {
        if c.has_stops() {
            DASH.to_string()
        } else {
            format!("{:.2}", hx_core::protein::isoelectric_point(c))
        }
    };
    let charge_of = |c: &Composition, ph: f64| {
        if c.has_stops() {
            DASH.to_string()
        } else {
            format!("{:.2}", hx_core::protein::charge(c, ph))
        }
    };
    // 「选中区段」列的值：没有选中区段时留空。
    let sv = |f: &dyn Fn(&Composition) -> String| sel.map(|c| f(&c)).unwrap_or_default();

    prop_line(
        &p, x, val_r, sel_r, 0.0, y, CellKind::Text,
        "Length", &len_of(&whole), &sv(&len_of), &font,
    );
    y += ROW_H;

    prop_expandable(
        ui, &p, x, content.right() - x, val_r, sel_r, y, true, CellKind::Text,
        "Molecular Weight",
        &mass_of(&whole),
        &sv(&mass_of),
        &mut app.props_expand[0],
        &font,
    );
    y += ROW_H;
    if app.props_expand[0] {
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Text, "average",
            &mass_of(&whole), &sv(&mass_of), &sub_font);
        y += SUB_H;
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Text, "monoisotopic",
            &mono_of(&whole), &sv(&mono_of), &sub_font);
        y += SUB_H;
    }

    prop_expandable(
        ui, &p, x, content.right() - x, val_r, sel_r, y, true, CellKind::Extinction,
        "Extinction Coefficient (280 nm)",
        &ext_of(&whole, false),
        &sv(&|c| ext_of(c, false)),
        &mut app.props_expand[1],
        &font,
    );
    y += ROW_H;
    if app.props_expand[1] {
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Extinction, "fully reduced",
            &ext_of(&whole, false), &sv(&|c| ext_of(c, false)), &sub_font);
        y += SUB_H;
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Extinction, "fully oxidized",
            &ext_of(&whole, true), &sv(&|c| ext_of(c, true)), &sub_font);
        y += SUB_H;
    }

    prop_expandable(
        ui, &p, x, content.right() - x, val_r, sel_r, y, true, CellKind::Text,
        "Absorbance (280 nm, 0.1%)",
        &abs_of(&whole, false),
        &sv(&|c| abs_of(c, false)),
        &mut app.props_expand[2],
        &font,
    );
    y += ROW_H;
    if app.props_expand[2] {
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Text, "fully reduced",
            &abs_of(&whole, false), &sv(&|c| abs_of(c, false)), &sub_font);
        y += SUB_H;
        prop_line(&p, x, val_r, sel_r, 20.0, y, CellKind::Text, "fully oxidized",
            &abs_of(&whole, true), &sv(&|c| abs_of(c, true)), &sub_font);
        y += SUB_H;
    }

    prop_line(
        &p, x, val_r, sel_r, 0.0, y, CellKind::Text,
        "Isoelectric Point (pI)",
        &pi_of(&whole),
        &sv(&pi_of),
        &font,
    );
    y += ROW_H;

    // Charge at pH：标签后跟一个 0.5 步长的下拉（SnapGene 的 1.0…14.0）
    prop_line(&p, x, val_r, sel_r, 0.0, y, CellKind::Text, "Charge at pH", "", "", &font);
    let combo = egui::Rect::from_min_size(egui::pos2(x + 84.0, y - 11.0), egui::vec2(74.0, 22.0));
    ui.scope_builder(egui::UiBuilder::new().max_rect(combo), |ui| {
        egui::ComboBox::from_id_salt("props_ph")
            .width(56.0)
            .selected_text(format!("{:.1}", app.props_ph))
            .show_ui(ui, |ui| {
                let mut v = 1.0f64;
                while v <= 14.0001 {
                    if ui
                        .selectable_label((app.props_ph - v).abs() < 0.05, format!("{v:.1}"))
                        .clicked()
                    {
                        app.props_ph = v;
                        ui.close();
                    }
                    v += 0.5;
                }
            });
    });
    let ph = app.props_ph;
    p.text(
        egui::pos2(val_r, y),
        Align2::RIGHT_CENTER,
        charge_of(&whole, ph),
        font.clone(),
        theme::TEXT,
    );
    if let Some(c) = sel {
        p.text(
            egui::pos2(sel_r, y),
            Align2::RIGHT_CENTER,
            charge_of(&c, ph),
            font.clone(),
            theme::TEXT,
        );
    }
    y += ROW_H + 10.0;

    // ---- 氨基酸组成表 ----
    // 两对 Number / Percent 列**总是**存在；没有选中区段时右两列留白（oracle 实测）。
    p.text(egui::pos2(x + 119.0, y), Align2::CENTER_CENTER, crate::i18n::text(ui.ctx(),"Amino Acid"), font.clone(), theme::TEXT);
    p.text(egui::pos2(x + 244.0, y), Align2::RIGHT_CENTER, crate::i18n::text(ui.ctx(),"Number"), font.clone(), theme::TEXT);
    p.text(egui::pos2(x + 319.0, y), Align2::RIGHT_CENTER, crate::i18n::text(ui.ctx(),"Percent"), font.clone(), theme::TEXT);
    p.text(egui::pos2(x + 429.0, y), Align2::RIGHT_CENTER, crate::i18n::text(ui.ctx(),"Number"), font.clone(), theme::TEXT);
    p.text(egui::pos2(x + 494.0, y), Align2::RIGHT_CENTER, crate::i18n::text(ui.ctx(),"Percent"), font.clone(), theme::TEXT);
    let sep = y + 13.0;
    p.line_segment(
        [egui::pos2(x, sep), egui::pos2(content.right(), sep)],
        egui::Stroke::new(1.0, theme::RULE),
    );
    // 20 种标准氨基酸（按单字母码排序；计数为 0 的也列出来，与 SnapGene 一致）。
    // 未知残基 / 终止符不单列，百分比仍按全长算。
    let rows: Vec<(char, &str, &str, usize, usize)> = hx_core::protein::STANDARD_AA
        .iter()
        .enumerate()
        .map(|(i, (c, tri, name))| (*c, *tri, *name, whole.counts[i], sel.map_or(0, |s| s.counts[i])))
        .collect();
    let n = rows.len();
    let total = whole.len().max(1) as f64;
    let sel_total = sel.map(|s| s.len()).unwrap_or(0).max(1) as f64;
    let has_sel = sel.is_some();
    let area = egui::Rect::from_min_max(egui::pos2(x, sep + 1.0), content.max);
    ui.scope_builder(egui::UiBuilder::new().max_rect(area), |ui| {
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                let (r, _) = ui.allocate_exact_size(
                    egui::vec2(area.width(), n as f32 * AA_H + 4.0),
                    egui::Sense::hover(),
                );
                let p = ui.painter();
                for (i, (c, tri, name, cnt, scnt)) in rows.iter().enumerate() {
                    let cy = r.top() + i as f32 * AA_H + AA_H / 2.0;
                    if i % 2 == 0 {
                        p.rect_filled(
                            egui::Rect::from_min_max(
                                egui::pos2(r.left(), cy - AA_H / 2.0),
                                egui::pos2(r.right(), cy + AA_H / 2.0),
                            ),
                            0.0,
                            theme::STRIPE,
                        );
                    }
                    let txt = |px: f32, align: Align2, s: &str| {
                        p.text(egui::pos2(px, cy), align, s, font.clone(), theme::TEXT);
                    };
                    txt(r.left(), Align2::LEFT_CENTER, &c.to_string());
                    txt(r.left() + 28.0, Align2::LEFT_CENTER, tri);
                    txt(r.left() + 73.0, Align2::LEFT_CENTER, name);
                    txt(r.left() + 244.0, Align2::RIGHT_CENTER, &nfmt_plain(*cnt));
                    txt(r.left() + 319.0, Align2::RIGHT_CENTER, &format!("{:.2}", *cnt as f64 / total * 100.0));
                    if has_sel {
                        txt(r.left() + 429.0, Align2::RIGHT_CENTER, &nfmt_plain(*scnt));
                        txt(r.left() + 494.0, Align2::RIGHT_CENTER, &format!("{:.2}", *scnt as f64 / sel_total * 100.0));
                    }
                }
            });
    });
}

/// `Properties` 页的一行：左侧名字、右侧值（整列右对齐）；`indent` 用于展开的明细行。
#[allow(clippy::too_many_arguments)]
fn prop_line(
    p: &egui::Painter,
    x: f32,
    val_r: f32,
    sel_r: f32,
    indent: f32,
    y: f32,
    kind: CellKind,
    name: &str,
    whole: &str,
    sel: &str,
    font: &egui::FontId,
) {
    p.text(
        egui::pos2(x + indent, y),
        Align2::LEFT_CENTER,
        crate::i18n::text(p.ctx(),name),
        font.clone(),
        theme::TEXT,
    );
    if !whole.is_empty() {
        cell(p, val_r, y, kind, whole, font);
    }
    if !sel.is_empty() {
        cell(p, sel_r, y, kind, sel, font);
    }
}

/// 画一个值单元格：普通文字，或带手工上标的消光系数。
fn cell(p: &egui::Painter, right: f32, y: f32, kind: CellKind, text: &str, font: &egui::FontId) {
    let numeric = !text.is_empty()
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || b == b',');
    if kind == CellKind::Extinction && numeric {
        let v: u32 = text.replace(',', "").parse().unwrap_or(0);
        paint_extinction(p, right, y, v, font);
    } else {
        p.text(
            egui::pos2(right, y),
            Align2::RIGHT_CENTER,
            text,
            font.clone(),
            theme::TEXT,
        );
    }
}

/// 可展开的属性行：左边一个三角（`▸` 收起 / `▾` 展开），点整行切换。
#[allow(clippy::too_many_arguments)]
fn prop_expandable(
    ui: &mut egui::Ui,
    p: &egui::Painter,
    x: f32,
    width: f32,
    val_r: f32,
    sel_r: f32,
    y: f32,
    open_now: bool,
    kind: CellKind,
    name: &str,
    whole: &str,
    sel: &str,
    open: &mut bool,
    font: &egui::FontId,
) {
    let hit = egui::Rect::from_min_max(
        egui::pos2(x - 16.0, y - 11.0),
        egui::pos2(x + width, y + 11.0),
    );
    if ui
        .interact(hit, ui.make_persistent_id(("prop_row", name)), egui::Sense::click())
        .clicked()
    {
        *open = !*open;
    }
    prop_triangle(p, egui::pos2(x - 10.0, y), open_now);
    prop_line(p, x, val_r, sel_r, 0.0, y, kind, name, whole, sel, font);
}

/// 展开三角（`▸` 收起 / `▾` 展开）。
fn prop_triangle(p: &egui::Painter, c: egui::Pos2, open: bool) {
    let s = 4.6;
    let pts = if open {
        vec![
            c + egui::vec2(-s, -s * 0.55),
            c + egui::vec2(s, -s * 0.55),
            c + egui::vec2(0.0, s * 0.8),
        ]
    } else {
        vec![
            c + egui::vec2(-s * 0.55, -s),
            c + egui::vec2(-s * 0.55, s),
            c + egui::vec2(s * 0.8, 0.0),
        ]
    };
    p.add(egui::Shape::convex_polygon(
        pts,
        theme::TEXT,
        egui::Stroke::NONE,
    ));
}

/// 右对齐画 `26,930 M⁻¹ cm⁻¹`：上标用更小字号抬高绘制（egui 默认字体缺 `⁻¹`）。
fn paint_extinction(p: &egui::Painter, right: f32, y: f32, v: u32, font: &egui::FontId) {
    let main = format!("{} M", nfmt(v as usize));
    let sup = egui::FontId::proportional(font.size - 3.5);
    let w = |s: &str, f: &egui::FontId| {
        p.layout_no_wrap(s.to_string(), f.clone(), theme::TEXT).size().x
    };
    let (w_main, w_sup, w_cm) = (w(&main, font), w("-1", &sup), w(" cm", font));
    let mut x = right - (w_main + 2.0 * w_sup + w_cm);
    p.text(egui::pos2(x, y), Align2::LEFT_CENTER, &main, font.clone(), theme::TEXT);
    x += w_main;
    p.text(
        egui::pos2(x, y - font.size * 0.33),
        Align2::LEFT_CENTER,
        "-1",
        sup.clone(),
        theme::TEXT,
    );
    x += w_sup;
    p.text(egui::pos2(x, y), Align2::LEFT_CENTER, " cm", font.clone(), theme::TEXT);
    x += w_cm;
    p.text(
        egui::pos2(x, y - font.size * 0.33),
        Align2::LEFT_CENTER,
        "-1",
        sup,
        theme::TEXT,
    );
}

/// `50,729.56` 这样的质量：整数部分按千位分组，保留两位小数。
fn mass_str(v: f64) -> String {
    let s = format!("{v:.2}");
    match s.split_once('.') {
        Some((int, frac)) => format!("{}.{frac}", nfmt(int.parse::<usize>().unwrap_or(0))),
        None => s,
    }
}

/// 控制列保持紧凑，剩余宽度分给主要数据列；窄窗口保留列宽并允许横向滚动。
fn column_widths<const N: usize>(
    available: f32,
    state: &crate::app::TableState,
    mut widths: [f32; N],
    stretch: usize,
) -> [f32; N] {
    let visible = (0..N).filter(|c| state.visible(*c)).count();
    let used: f32 = widths
        .iter()
        .enumerate()
        .filter(|(c, _)| state.visible(*c))
        .map(|(_, w)| *w)
        .sum::<f32>()
        + visible.saturating_sub(1) as f32 * COLUMN_GAP;
    let stretch = if state.visible(stretch) {
        Some(stretch)
    } else {
        (0..N).rev().find(|c| state.visible(*c))
    };
    if let Some(c) = stretch {
        widths[c] += (available - used).max(0.0);
    }
    widths
}

/// 白底表格用浅色横线划分行，避免紧密的交替色块影响长表阅读。
fn end_table_row(ui: &mut egui::Ui) {
    ui.end_row();
    let y = ui.cursor().top();
    ui.painter().line_segment(
        [
            egui::pos2(ui.max_rect().left(), y),
            egui::pos2(ui.max_rect().right(), y),
        ],
        egui::Stroke::new(0.7, theme::BORDER),
    );
}

fn detail_cell(ui: &mut egui::Ui, width: f32, text: String) {
    let galley = ui
        .painter()
        .layout(text, egui::FontId::monospace(11.0), theme::TEXT_DIM, width);
    // Grid 会垂直居中单元格，先量出换行后的实际高度，避免子 UI 从行中间向下溢出。
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(width, galley.size().y + 8.0),
        egui::Sense::hover(),
    );
    ui.painter()
        .galley(rect.min + egui::vec2(0.0, 4.0), galley, theme::TEXT_DIM);
}

/// 页头：粗体标题 + 右侧浅色说明。
fn page_header(ui: &mut egui::Ui, title: &str, note: &str) {
    ui.add_space(8.0);
    ui.horizontal(|ui| {
        ui.add_space(10.0);
        ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),title)).strong().size(13.5));
        ui.separator();
        ui.label(egui::RichText::new(note).size(11.5).color(theme::TEXT_DIM));
    });
    ui.add_space(6.0);
}

/// 表头（SnapGene 的 `ColumnHeaderContextMenu`）：可点列头排序、右键隐藏列；
/// 返回本次请求的排序 `(列下标, 升序)`。
fn table_header(
    ui: &mut egui::Ui,
    state: &mut crate::app::TableState,
    cols: &[(&str, bool)],
    widths: &[f32],
) -> Option<(usize, bool)> {
    let mut want: Option<(usize, bool)> = None;
    let mut hide: Option<usize> = None;
    for (ci, (name, sortable)) in cols.iter().enumerate() {
        if !state.visible(ci) {
            continue;
        }
        let arrow = match state.sort {
            Some((c, asc)) if c == ci => {
                if asc {
                    " ^"
                } else {
                    " v"
                }
            }
            _ => "",
        };
        // 用无边框按钮（与 Label 外观一致但可靠响应点击）
        let resp = ui
            .allocate_ui_with_layout(
                egui::vec2(widths[ci], 30.0),
                egui::Layout::left_to_right(egui::Align::Center),
                |ui| {
                    ui.set_min_size(egui::vec2(widths[ci], 30.0));
                    ui.add(
                        egui::Button::new(
                            egui::RichText::new(format!("{}{arrow}",crate::i18n::text(ui.ctx(),name)))
                                .strong()
                                .size(12.0),
                        )
                        .frame(false),
                    )
                },
            )
            .inner;
        if *sortable {
            // 左键点列头：切换该列升/降序（SnapGene 的列头排序）
            if resp.clicked() {
                let next = match state.sort {
                    Some((c, true)) if c == ci => (ci, false),
                    _ => (ci, true),
                };
                want = Some(next);
            }
            resp.context_menu(|ui| {
                if ui.button(crate::i18n::text(ui.ctx(),"Sort ascending")).clicked() {
                    want = Some((ci, true));
                    ui.close();
                }
                if ui.button(crate::i18n::text(ui.ctx(),"Sort descending")).clicked() {
                    want = Some((ci, false));
                    ui.close();
                }
                ui.separator();
                if ui.button(crate::i18n::text(ui.ctx(),"Hide column")).clicked() {
                    hide = Some(ci);
                    ui.close();
                }
            });
        }
    }
    end_table_row(ui);
    if let Some(h) = hide {
        state.hidden.insert(h);
    }
    if want.is_some() {
        state.sort = want;
    }
    want
}

/// 页头右侧的「Columns ▾」：勾选列可见性（隐藏列后从这里恢复）。
fn columns_button(ui: &mut egui::Ui, state: &mut crate::app::TableState, cols: &[&str]) {
    ui.menu_button(crate::i18n::text(ui.ctx(),"Columns"), |ui| {
        for (ci, name) in cols.iter().enumerate() {
            if name.is_empty() {
                continue;
            }
            let mut on = state.visible(ci);
            if ui.checkbox(&mut on, crate::i18n::text(ui.ctx(),name)).changed() {
                if on {
                    state.hidden.remove(&ci);
                } else {
                    state.hidden.insert(ci);
                }
            }
        }
    });
}

// ---------------------------------------------------------------------------
// Enzymes
// ---------------------------------------------------------------------------

pub fn enzymes(ui: &mut egui::Ui, app: &mut HelixApp) {
    let n_sel = app.sites.selected.len();
    let mut tbl = app.tbl_enzymes.clone();
    let len = app.doc.as_ref().map(|d| d.len()).unwrap_or(0);
    let circular = app.doc.as_ref().map(|d| d.circular()).unwrap_or(false);
    let mut chosen: Option<EnzymePreset> = None;
    let mut remove: Option<usize> = None;
    let mut add: Option<String> = None;
    let mut add_buf = app.enzyme_search.clone();

    // ---- 顶栏（对齐 SnapGene）：找到酶 / 添加酶 ▾ / 操作 ▾ …… N 酶 / 集合下拉 / 排序 ----
    ui.horizontal(|ui| {
        ui.add(
            egui::TextEdit::singleline(&mut app.enzyme_search)
                .desired_width(120.0)
                .hint_text(crate::i18n::text(ui.ctx(),"Find enzyme")),
        );
        ui.menu_button(crate::i18n::text(ui.ctx(),"Add enzyme"), |ui| {
            ui.set_min_width(220.0);
            ui.label(
                egui::RichText::new(crate::i18n::text(ui.ctx(),"Type a name and press Add"))
                    .size(10.5)
                    .color(theme::TEXT_DIM),
            );
            let te = ui.add(
                egui::TextEdit::singleline(&mut add_buf)
                    .desired_width(200.0)
                    .hint_text(crate::i18n::text(ui.ctx(),"e.g. EcoRI")),
            );
            let enter = te.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            ui.horizontal(|ui| {
                if ui.button(crate::i18n::text(ui.ctx(),"Add")).clicked() || enter {
                    add = Some(add_buf.clone());
                    ui.close();
                }
                if ui.button(crate::i18n::text(ui.ctx(),"Clear all")).clicked() {
                    chosen = Some(EnzymePreset::None);
                    ui.close();
                }
            });
        });
        ui.menu_button(crate::i18n::text(ui.ctx(),"Actions"), |ui| {
            // 选中行（若没有选中，则对全部已选酶生效）
            let n = app.sites.selected.len();
            let targets: Vec<usize> = if app.selected_enzymes.is_empty() {
                (0..n).collect()
            } else {
                app.selected_enzymes
                    .iter()
                    .copied()
                    .filter(|k| *k < n)
                    .collect()
            };
            let has = !targets.is_empty();
            if ui
                .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Highlight Selected Enzyme Sites")))
                .clicked()
            {
                for k in &targets {
                    app.highlighted_enzymes.insert(*k);
                }
                ui.close();
            }
            if ui
                .add_enabled(
                    !app.highlighted_enzymes.is_empty(),
                    egui::Button::new(crate::i18n::text(ui.ctx(),"Clear Highlighting for All Enzyme Sites")),
                )
                .clicked()
            {
                app.highlighted_enzymes.clear();
                ui.close();
            }
            ui.separator();
            if ui
                .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Show Selected Enzymes")))
                .clicked()
            {
                for k in &targets {
                    app.hidden_enzymes.remove(k);
                }
                ui.close();
            }
            if ui
                .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Hide Selected Enzymes")))
                .clicked()
            {
                for k in &targets {
                    app.hidden_enzymes.insert(*k);
                }
                ui.close();
            }
            if ui
                .add_enabled(
                    has,
                    egui::Button::new(crate::i18n::text(ui.ctx(),"Remove Selected from Chosen Enzymes")),
                )
                .clicked()
            {
                let mut idxs: Vec<usize> = targets
                    .iter()
                    .filter_map(|k| app.sites.selected.get(*k).copied())
                    .collect();
                idxs.sort_unstable();
                idxs.reverse();
                for i in idxs {
                    app.remove_enzyme(i);
                }
                app.selected_enzymes.clear();
                ui.close();
            }
            let _ = ui.add_enabled(
                false,
                egui::Button::new(crate::i18n::text(ui.ctx(),"Remove Selected Enzyme Site…  (mutagenesis, not implemented)")),
            );
            ui.separator();
            // 按选择范围隐藏（SnapGene: Hide Enzymes That Cut Inside/Outside the Selection）
            let sel = app.effective_range();
            if ui
                .add_enabled(
                    sel.is_some(),
                    egui::Button::new(crate::i18n::text(ui.ctx(),"Hide Enzymes That Cut Outside the Selection")),
                )
                .clicked()
            {
                if let Some((s0, e0)) = sel {
                    for k in 0..app.sites.selected.len() {
                        let inside = app
                            .sites
                            .hits
                            .get(k)
                            .map(|h| h.iter().all(|p| (*p as usize) >= s0 && (*p as usize) < e0))
                            .unwrap_or(false);
                        if !inside {
                            app.hidden_enzymes.insert(k);
                        }
                    }
                }
                ui.close();
            }
            if ui
                .add_enabled(
                    sel.is_some(),
                    egui::Button::new(crate::i18n::text(ui.ctx(),"Hide Enzymes That Cut Inside the Selection")),
                )
                .clicked()
            {
                if let Some((s0, e0)) = sel {
                    for k in 0..app.sites.selected.len() {
                        let inside = app
                            .sites
                            .hits
                            .get(k)
                            .map(|h| h.iter().any(|p| (*p as usize) >= s0 && (*p as usize) < e0))
                            .unwrap_or(false);
                        if inside {
                            app.hidden_enzymes.insert(k);
                        }
                    }
                }
                ui.close();
            }
            ui.separator();
            if ui.button(crate::i18n::text(ui.ctx(),"Export Enzyme Set…")).clicked() {
                app.export_enzyme_set();
                ui.close();
            }
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Save Enzyme Set…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Manage Enzyme Sets…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Enzyme Database")));
            ui.separator();
            if ui.button(crate::i18n::text(ui.ctx(),"Show All Enzymes")).clicked() {
                app.hidden_enzymes.clear();
                app.show_sites = true;
                ui.close();
            }
            if ui.button(crate::i18n::text(ui.ctx(),"Hide All Enzymes")).clicked() {
                for k in 0..app.sites.selected.len() {
                    app.hidden_enzymes.insert(k);
                }
                ui.close();
            }
            ui.separator();
            if ui
                .selectable_label(!app.sites_as_lines, crate::i18n::text(ui.ctx(),"Show Sites as Numbers"))
                .clicked()
            {
                app.sites_as_lines = false;
                ui.close();
            }
            if ui
                .selectable_label(app.sites_as_lines, crate::i18n::text(ui.ctx(),"Show Sites as Lines"))
                .clicked()
            {
                app.sites_as_lines = true;
                ui.close();
            }
            ui.separator();
            columns_button(ui, &mut tbl, &["", "", "Enzyme", "Sites", "Location", ""]);
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(8.0);
            // 酶组下拉
            let mut chosen_set: Option<EnzymePreset> = None;
            egui::ComboBox::from_id_salt("enzyme_set")
                .selected_text(app.enzyme_preset.label())
                .width(160.0)
                .show_ui(ui, |ui| {
                    for p in EnzymePreset::ALL {
                        if ui
                            .selectable_label(app.enzyme_preset == p, p.label())
                            .clicked()
                        {
                            chosen_set = Some(p);
                            ui.close();
                        }
                    }
                });
            if let Some(p) = chosen_set {
                chosen = Some(p);
            }
            ui.label(
            egui::RichText::new(crate::ui_format!(ui.ctx(),"{} enzymes", nfmt(n_sel)))
                    .size(11.5)
                    .color(theme::TEXT_DIM),
            );
        });
    });

    // ---- 无切点内切酶横幅（SnapGene: N of M enzymes in the chosen set do not cut）----
    let noncutters = app
        .sites
        .selected
        .iter()
        .enumerate()
        .filter(|(k, _)| app.sites.hits.get(*k).map(|h| h.is_empty()).unwrap_or(true))
        .count();
    if noncutters > 0 {
        ui.add_space(4.0);
        let (rect, _) =
            ui.allocate_exact_size(egui::vec2(ui.available_width(), 26.0), egui::Sense::hover());
        let p = ui.painter();
        p.rect_filled(rect, 0.0, Color32::from_rgb(0xe8, 0xf0, 0xfe));
        p.rect_stroke(
            rect,
            0.0,
            egui::Stroke::new(0.8, Color32::from_rgb(0xc5, 0xd8, 0xf7)),
            egui::StrokeKind::Inside,
        );
        let info = egui::Rect::from_min_size(
            egui::Pos2::new(rect.left() + 6.0, rect.center().y - 7.0),
            egui::vec2(14.0, 14.0),
        );
        p.circle_filled(info.center(), 6.5, theme::ACCENT);
        p.text(
            info.center(),
            Align2::CENTER_CENTER,
            "i",
            egui::FontId::proportional(10.0),
            Color32::WHITE,
        );
        p.text(
            egui::Pos2::new(rect.left() + 28.0, rect.center().y),
            Align2::LEFT_CENTER,
            format!(
                "{} of {} enzymes in the chosen set do not cut",
                nfmt(noncutters),
                nfmt(n_sel)
            ),
            egui::FontId::proportional(12.0),
            theme::TEXT,
        );
    }

    if n_sel == 0 {
        ui.add_space(10.0);
        ui.horizontal(|ui| {
            ui.add_space(10.0);
            ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"No enzymes selected.")).color(theme::TEXT_DIM));
        });
        app.tbl_enzymes = tbl;
        if let Some(p) = chosen {
            app.apply_preset(p);
        }
        return;
    }

    let selected = app.sites.selected.clone();
    let hits = app.sites.hits.clone();
    let mut rows: Vec<(usize, usize)> = selected
        .iter()
        .enumerate()
        .map(|(k, &idx)| (k, idx))
        .collect();
    let needle = app.enzyme_search.trim().to_ascii_lowercase();
    rows.retain(|(_, idx)| {
        hx_enzymes::ENZYMES[*idx]
            .name
            .to_ascii_lowercase()
            .contains(&needle)
    });
    if let Some((ci, asc)) = tbl.sort {
        rows.sort_by(|a, b| {
            let ea = &hx_enzymes::ENZYMES[a.1];
            let eb = &hx_enzymes::ENZYMES[b.1];
            let ka = match ci {
                2 => ea.name.to_string(),
                3 => format!("{:08}", hits.get(a.0).map(|h| h.len()).unwrap_or(0)),
                _ => format!(
                    "{:010}",
                    hits.get(a.0)
                        .and_then(|h| h.first())
                        .map(|p| hx_enzymes::display_pos(*p as usize, ea.cut, len, false))
                        .unwrap_or(0)
                ),
            };
            let kb = match ci {
                2 => eb.name.to_string(),
                3 => format!("{:08}", hits.get(b.0).map(|h| h.len()).unwrap_or(0)),
                _ => format!(
                    "{:010}",
                    hits.get(b.0)
                        .and_then(|h| h.first())
                        .map(|p| hx_enzymes::display_pos(*p as usize, eb.cut, len, false))
                        .unwrap_or(0)
                ),
            };
            if asc {
                ka.cmp(&kb)
            } else {
                kb.cmp(&ka)
            }
        });
    }

    let ecols: [(&str, bool); 6] = [
        ("", false),
        ("", false),
        ("Enzyme", true),
        ("Sites", true),
        ("Location", true),
        ("", false),
    ];
    let nvis = (0..ecols.len()).filter(|c| tbl.visible(*c)).count().max(1);
    let as_lines = app.sites_as_lines;
    // 表格总宽（Grid 单元格里的 available_width 是剩余宽度，量不出整表宽度）
    let panel_w = (ui.available_width() - 12.0).max(0.0);
    let widths = column_widths(panel_w, &tbl, [24.0, 24.0, 140.0, 72.0, 160.0, 24.0], 4);
    let mut expand_toggle: Option<usize> = None;
    let mut vis_toggle: Option<usize> = None;
    let mut pick_enzyme: Option<usize> = None;
    ui.add_space(4.0);
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new("enzyme_table")
                .num_columns(nvis)
                .striped(false)
                .spacing([COLUMN_GAP, 0.0])
                .min_col_width(24.0)
                .min_row_height(32.0)
                .show(ui, |ui| {
                    table_header(ui, &mut tbl, &ecols, &widths);
                    for (k, idx) in &rows {
                        let (k, idx) = (*k, *idx);
                        let e = &hx_enzymes::ENZYMES[idx];
                        let h = hits.get(k).map(|v| v.as_slice()).unwrap_or(&[]);
                        let hidden = app.hidden_enzymes.contains(&k);
                        // 可见性勾选框（SnapGene 最左列）
                        if tbl.visible(0) {
                            let mut vis = !hidden;
                            if ui.checkbox(&mut vis, "").changed() {
                                vis_toggle = Some(k);
                            }
                        }
                        // 展开箭头
                        let expanded = app.expanded_enzymes.contains(&k);
                        if tbl.visible(1) {
                            if ui
                                .add(
                                    egui::Button::new(if expanded { "v" } else { ">" })
                                        .frame(false),
                                )
                                .clicked()
                            {
                                expand_toggle = Some(k);
                            }
                        }
                        if tbl.visible(2) {
                            let is_picked = app.selected_enzymes.contains(&k);
                            let highlighted = app.highlighted_enzymes.contains(&k);
                            let col = if highlighted {
                                theme::ACCENT
                            } else {
                                theme::TEXT
                            };
                            if ui
                                .selectable_label(
                                    is_picked,
                                    egui::RichText::new(e.name).strong().color(col),
                                )
                                .on_hover_text(enzyme_tip(ui.ctx(),e, h.len()))
                                .clicked()
                            {
                                pick_enzyme = Some(k);
                            }
                        }
                        if tbl.visible(3) {
                            ui.label(
                                egui::RichText::new(nfmt(h.len()))
                                    .monospace()
                                    .color(if hidden { theme::TEXT_DIM } else { theme::TEXT }),
                            );
                        }
                        if tbl.visible(4) {
                            if as_lines {
                                // SnapGene 的「以线显示位点」：每行一条全长轨道 + 刻度
                                let (rect, _) = ui.allocate_exact_size(
                                    egui::vec2(widths[4], 16.0),
                                    egui::Sense::hover(),
                                );
                                let track = egui::Rect::from_min_max(
                                    egui::Pos2::new(rect.left(), rect.center().y),
                                    egui::Pos2::new(rect.right(), rect.center().y + 1.0),
                                );
                                ui.painter().rect_filled(
                                    track,
                                    0.0,
                                    theme::with_alpha(theme::BORDER, 120),
                                );
                                if len > 0 {
                                    for &site in h {
                                        let x =
                                            rect.left() + (site as f32 / len as f32) * rect.width();
                                        ui.painter().line_segment(
                                            [
                                                egui::Pos2::new(x, rect.center().y - 6.0),
                                                egui::Pos2::new(x, rect.center().y + 6.0),
                                            ],
                                            egui::Stroke::new(1.2, theme::TEXT),
                                        );
                                    }
                                }
                            } else {
                                // 数字：空格分隔、不带到千位分隔（SnapGene 写法）
                                let mut pos: Vec<String> = h
                                    .iter()
                                    .take(14)
                                    .map(|p| {
                                        hx_enzymes::display_pos(*p as usize, e.cut, len, circular)
                                            .to_string()
                                    })
                                    .collect();
                                if h.len() > 14 {
                                    pos.push("…".into());
                                }
                                ui.add(
                                    egui::Label::new(
                                        egui::RichText::new(pos.join(" "))
                                            .monospace()
                                            .size(11.5)
                                            .color(if hidden {
                                                theme::TEXT_DIM
                                            } else {
                                                theme::TEXT
                                            }),
                                    )
                                    .halign(egui::Align::Min),
                                );
                            }
                        }
                        if tbl.visible(5) && ui.small_button("×").on_hover_text(crate::i18n::text(ui.ctx(),"Remove")).clicked()
                        {
                            remove = Some(idx);
                        }
                        end_table_row(ui);
                        if expanded {
                            // 展开明细（SnapGene 的行内展开：识别序列 / 切点 / 突出端）
                            if tbl.visible(0) {
                                ui.label("");
                            }
                            if tbl.visible(1) {
                                ui.label("");
                            }
                            if tbl.visible(2) {
                                detail_cell(
                                    ui,
                                    widths[2],
                                    format!(
                                        "Recognition site: {}   ·   cut {} / {}   ·   {}",
                                        e.site,
                                        e.cut,
                                        e.reverse_cut,
                                        overhang_label(e.overhang)
                                    ),
                                );
                            }
                            if tbl.visible(3) {
                                ui.label("");
                            }
                            if tbl.visible(4) {
                                ui.label(
                                    egui::RichText::new(if h.is_empty() {
                                        crate::i18n::text(ui.ctx(),"does not cut this sequence").to_string()
                                    } else {
                                        crate::ui_format!(ui.ctx(),"{} site(s) in this sequence", h.len())
                                    })
                                    .size(11.0)
                                    .color(theme::TEXT_DIM),
                                );
                            }
                            if tbl.visible(5) {
                                ui.label("");
                            }
                            end_table_row(ui);
                        }
                    }
                });
        });
    app.tbl_enzymes = tbl;
    if let Some(k) = pick_enzyme {
        if !app.selected_enzymes.insert(k) {
            app.selected_enzymes.remove(&k);
        }
    }
    if let Some(k) = vis_toggle {
        if !app.hidden_enzymes.insert(k) {
            app.hidden_enzymes.remove(&k);
        }
    }
    if let Some(k) = expand_toggle {
        if !app.expanded_enzymes.insert(k) {
            app.expanded_enzymes.remove(&k);
        }
    }
    if let Some(idx) = remove {
        app.remove_enzyme(idx);
    }
    if let Some(p) = chosen {
        app.apply_preset(p);
    }
    if let Some(name) = add {
        let name = name.trim().to_string();
        if !name.is_empty() && app.add_enzyme(&name) {
            app.status = format!("Added {name}");
        }
    }
}

/// 酶 tooltip（SnapGene 的 EnzymeToolTipAttributes）。
fn enzyme_tip(ctx:&egui::Context,e: &hx_enzymes::Enzyme, n: usize) -> String {
    [e.name.to_string(),crate::ui_format!(ctx,"Recognition site: {}",e.site),crate::ui_format!(ctx,"Cut: {} / {}",e.cut,e.reverse_cut),crate::ui_format!(ctx,"Overhang: {}",crate::i18n::text(ctx,overhang_label(e.overhang))),crate::ui_format!(ctx,"Sites in this sequence: {}",n)].join("\n")
}

fn overhang_label(o: hx_enzymes::Overhang) -> &'static str {
    match o {
        hx_enzymes::Overhang::Blunt => "Blunt",
        hx_enzymes::Overhang::FivePrime => "5' overhang",
        hx_enzymes::Overhang::ThreePrime => "3' overhang",
        hx_enzymes::Overhang::Iis => "IIS",
    }
}

// ---------------------------------------------------------------------------
// Features
// ---------------------------------------------------------------------------

pub fn features(ui: &mut egui::Ui, app: &mut HelixApp) {
    /// 一行特征（绘制前抽取的只读快照，避免借用冲突）。
    struct Row {
        idx: usize,
        name: String,
        kind: String,
        /// 首段起点（0-based）。
        start: usize,
        /// 末段终点（0-based 半开；跨原点时 `< start`）。
        end: usize,
        /// 各段长度之和（SnapGene 的 Size 列）。
        len: usize,
        directionality: i32,
        color: egui::Color32,
        n_seg: usize,
        segs: Vec<(usize, usize)>,
        /// 提示气泡里的限定符小节（`PRODUCT`/`NOTE`…，按文件顺序，名字转大写）。
        sections: Vec<(String, String)>,
    }

    /// 悬浮提示（对齐 SnapGene 的悬浮卡：名称 + 类型 → 位置 + 长度 → 各片段 → 限定符小节）。
    ///
    /// 与 SnapGene 的已知差异：真 SnapGene 对 CDS 只写 `GENETIC CODE / Standard` 与
    /// 计算出来的 `219 amino acids = 25.7 kDa`，我们按原文写限定符小节；
    /// `/translation`（整条氨基酸序列）SnapGene 不显示，我们也不显示。
    fn feature_tip(r: &Row, prefix: &str, unit: &str) -> String {
        let mut t = format!("{}   {}", r.name, r.kind);
        t.push_str(&format!(
            "\n{}{} .. {}{}   {} {}",
            prefix,
            nfmt_plain(r.start + 1),
            prefix,
            nfmt_plain(r.end),
            nfmt_plain(r.len),
            unit
        ));
        if r.n_seg > 1 {
            t.push_str(&format!("\n{} SEGMENTS", r.n_seg));
            for (i, (a, b)) in r.segs.iter().take(3).enumerate() {
                // SnapGene 的单残基片段只写一个位置（`2  p60`），区间才写 `p3 .. p4`
                let loc = if *b == *a + 1 {
                    format!("{}{}", prefix, nfmt_plain(*a + 1))
                } else {
                    format!(
                        "{}{} .. {}{}",
                        prefix,
                        nfmt_plain(*a + 1),
                        prefix,
                        nfmt_plain(*b)
                    )
                };
                t.push_str(&format!("\n{}   {}", i + 1, loc));
            }
            if r.n_seg > 3 {
                t.push_str(&format!("\n+{} more segments", r.n_seg - 3));
            }
        }
        // 限定符小节（SnapGene 把 `NOTE` 放最后）
        for (k, v) in &r.sections {
            t.push_str(&format!("\n{}\n{}", k, v));
        }
        t
    }

    let len = app.doc.as_ref().map(|d| d.len()).unwrap_or(0);
    let unit = app.doc.as_ref().map(|d| d.unit()).unwrap_or("bp");
    // 提示气泡/状态栏的坐标前缀：蛋白写 `p3 .. p414`，核苷酸不写（oracle）
    let doc_prefix = app
        .doc
        .as_ref()
        .map(|d| d.kind().position_prefix())
        .unwrap_or("");
    let rows: Vec<Row> = app
        .doc
        .as_ref()
        .map(|doc| {
            let specs = doc.feature_specs();
            doc.features()
                .iter()
                .enumerate()
                .map(|(idx, f)| {
                    let segs: Vec<(usize, usize)> = f
                        .segments
                        .iter()
                        .map(|s| (s.range.start, s.range.end))
                        .collect();
                    // oracle（SnapGene 8.2.2）：Location = 首段起点 .. 末段终点、
                    // Size = 各段长度之和、跨原点照原样写（无 `*`）
                    let (start, end, total, _wrapped) = hx_render::feature_span(&segs, len);
                    let color = f
                        .segments
                        .first()
                        .map(|s| theme::hex(&s.color, theme::FEATURE_FALLBACK))
                        .unwrap_or(theme::FEATURE_FALLBACK);
                    // 提示气泡里的限定符小节：块 10 的两个视图按文档顺序一一对应
                    let mut sections: Vec<(String, String)> = specs
                        .get(idx)
                        .filter(|s| s.name == f.name)
                        .map(|s| {
                            s.other_quals
                                .iter()
                                .filter(|(n, _)| n != "translation")
                                .map(|(n, v)| (n.replace('_', " ").to_ascii_uppercase(), v.clone()))
                                .collect()
                        })
                        .unwrap_or_default();
                    // SnapGene 把 NOTE 排在最后（oracle：GENETIC CODE / PRODUCT / NOTE）
                    if let Some(note) = specs
                        .get(idx)
                        .filter(|s| s.name == f.name)
                        .map(|s| s.note.trim())
                        .filter(|n| !n.is_empty())
                    {
                        sections.push(("NOTE".to_string(), note.to_string()));
                    }
                    Row {
                        idx,
                        name: f.name.clone(),
                        kind: f.kind.clone(),
                        start,
                        end,
                        len: total,
                        directionality: f.directionality,
                        color,
                        n_seg: segs.len(),
                        segs,
                        sections,
                    }
                })
                .collect()
        })
        .unwrap_or_default();

    let mut tbl = app.tbl_features.clone();
    // ---- 顶栏：找到特征 / 操作 ▾ …… N 特征 ----
    ui.horizontal(|ui| {
        ui.add(
            egui::TextEdit::singleline(&mut app.feature_search)
                .desired_width(140.0)
                .hint_text(crate::i18n::text(ui.ctx(),"Find features")),
        );
        ui.menu_button(crate::i18n::text(ui.ctx(),"Actions"), |ui| {
            let targets: Vec<usize> = if app.selected_features.is_empty() {
                (0..rows.len()).collect()
            } else {
                app.selected_features.iter().copied().collect()
            };
            let has = !targets.is_empty();
            let one = app.selected_feature;
            if ui.button(crate::i18n::text(ui.ctx(),"Add Feature…")).clicked() {
                app.open_add_feature_dialog(false);
                ui.close();
            }
            if ui.button(crate::i18n::text(ui.ctx(),"Add Translated Feature…")).clicked() {
                app.open_add_feature_dialog(true);
                ui.close();
            }
            if ui
                .add_enabled(one.is_some(), egui::Button::new(crate::i18n::text(ui.ctx(),"Edit Feature…")))
                .clicked()
            {
                app.open_edit_feature_dialog();
                ui.close();
            }
            if ui
                .add_enabled(one.is_some(), egui::Button::new(crate::i18n::text(ui.ctx(),"Remove Feature")))
                .clicked()
            {
                app.remove_selected_feature();
                ui.close();
            }
            if ui
                .add_enabled(one.is_some(), egui::Button::new(crate::i18n::text(ui.ctx(),"Duplicate Feature…")))
                .clicked()
            {
                app.open_duplicate_feature_dialog();
                ui.close();
            }
            if ui
                .add_enabled(one.is_some(), egui::Button::new(crate::i18n::text(ui.ctx(),"Feature Color…")))
                .clicked()
            {
                app.feature_color_popup = true;
                ui.close();
            }
            if ui
                .add_enabled(one.is_some(), egui::Button::new(crate::i18n::text(ui.ctx(),"Show or Hide Feature")))
                .clicked()
            {
                app.toggle_selected_feature_visible();
                ui.close();
            }
            ui.separator();
            if ui
                .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Show Selected Features")))
                .clicked()
            {
                for i in &targets {
                    app.hidden_features.remove(i);
                }
                ui.close();
            }
            if ui
                .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Hide Selected Features")))
                .clicked()
            {
                for i in &targets {
                    app.hidden_features.insert(*i);
                }
                ui.close();
            }
            if ui.button(crate::i18n::text(ui.ctx(),"Show All Features")).clicked() {
                app.show_features = true;
                app.hidden_features.clear();
                ui.close();
            }
            if ui.button(crate::i18n::text(ui.ctx(),"Hide All Features")).clicked() {
                for i in 0..rows.len() {
                    app.hidden_features.insert(i);
                }
                ui.close();
            }
            ui.separator();
            if ui.button(crate::i18n::text(ui.ctx(),"Export Feature Data…")).clicked() {
                app.export_feature_data();
                ui.close();
            }
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Import Features…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Detect Common Features…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Browse Common Features")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Manage Feature Types…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Add to Common Features…")));
            let _ = ui.add_enabled(false, egui::Button::new(crate::i18n::text(ui.ctx(),"Show Codon Frequencies")));
            ui.separator();
            if ui.button(crate::i18n::text(ui.ctx(),"Select Feature Range")).clicked() {
                app.select_feature_range();
                ui.close();
            }
            if ui.button(crate::i18n::text(ui.ctx(),"Center on Feature")).clicked() {
                if let Some(i) = app.selected_feature {
                    app.center_on_feature(i);
                    app.tab = Tab::Map;
                }
                ui.close();
            }
            ui.separator();
            columns_button(
                ui,
                &mut tbl,
                &[
                    "",
                    "",
                    "Feature",
                    "Location",
                    "Size",
                    "Color",
                    "Direction",
                    "Type",
                ],
            );
            let _ = &mut tbl;
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(8.0);
            // oracle：未选中写 `8 Features`，选中若干写 `1 of 68 Features selected`
            let text = match app.selected_features.len() {
                0 => format!("{} Features", nfmt(rows.len())),
                k => format!("{} of {} Features selected", nfmt(k), nfmt(rows.len())),
            };
            ui.label(
                egui::RichText::new(text)
                    .size(11.5)
                    .color(theme::TEXT_DIM),
            );
        });
    });

    let needle = app.feature_search.trim().to_ascii_lowercase();
    let mut clicked: Option<usize> = None;
    let mut toggle_hidden: Option<usize> = None;
    let mut expand_toggle: Option<usize> = None;
    let mut pick_feature: Option<usize> = None;
    let selected = app.selected_feature;

    // SnapGene 的 Features 视图默认**按特征名升序**（列头「特征 ↑」）
    let mut order: Vec<usize> = (0..rows.len()).collect();
    match tbl.sort {
        None => order.sort_by(|a, b| {
            rows[*a]
                .name
                .to_ascii_lowercase()
                .cmp(&rows[*b].name.to_ascii_lowercase())
                .then_with(|| rows[*a].idx.cmp(&rows[*b].idx))
        }),
        Some((ci, asc)) => order.sort_by(|a, b| {
            let (ra, rb) = (&rows[*a], &rows[*b]);
            let ka = match ci {
                2 => ra.name.to_ascii_lowercase(),
                3 => format!("{:010}", ra.start),
                4 => format!("{:010}", ra.len),
                5 => format!(
                    "#{:02x}{:02x}{:02x}",
                    ra.color.r(),
                    ra.color.g(),
                    ra.color.b()
                ),
                6 => format!("{}", ra.directionality),
                _ => ra.kind.to_ascii_lowercase(),
            };
            let kb = match ci {
                2 => rb.name.to_ascii_lowercase(),
                3 => format!("{:010}", rb.start),
                4 => format!("{:010}", rb.len),
                5 => format!(
                    "#{:02x}{:02x}{:02x}",
                    rb.color.r(),
                    rb.color.g(),
                    rb.color.b()
                ),
                6 => format!("{}", rb.directionality),
                _ => rb.kind.to_ascii_lowercase(),
            };
            if asc {
                ka.cmp(&kb)
            } else {
                kb.cmp(&ka)
            }
        }),
    }

    let fcols: [(&str, bool); 8] = [
        ("", false),
        ("", false),
        ("Feature", true),
        ("Location", true),
        ("Size", true),
        ("Color", true),
        ("Direction", true),
        ("Type", true),
    ];
    let nvis = (0..fcols.len()).filter(|c| tbl.visible(*c)).count().max(1);
    let panel_w = (ui.available_width() - 12.0).max(0.0);
    let widths = column_widths(
        panel_w,
        &tbl,
        [
            24.0,
            24.0,
            (panel_w * 0.18).clamp(180.0, 340.0),
            130.0,
            100.0,
            80.0,
            90.0,
            160.0,
        ],
        7,
    );
    ui.add_space(4.0);
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new("feature_table")
                .num_columns(nvis)
                .striped(false)
                .spacing([COLUMN_GAP, 0.0])
                .min_col_width(24.0)
                .min_row_height(42.0)
                .show(ui, |ui| {
                    table_header(ui, &mut tbl, &fcols, &widths);
                    for oi in order.clone() {
                        let r = &rows[oi];
                        if !needle.is_empty() {
                            let hay = format!("{} {}", r.name, r.kind).to_ascii_lowercase();
                            if !hay.contains(&needle) {
                                continue;
                            }
                        }
                        let is_sel = selected == Some(r.idx);
                        let hidden = app.hidden_features.contains(&r.idx);
                        // 可见性勾选框
                        if tbl.visible(0) {
                            let mut vis = !hidden;
                            if ui.checkbox(&mut vis, "").changed() {
                                toggle_hidden = Some(r.idx);
                            }
                        }
                        // 展开箭头
                        let expanded = app.expanded_features.contains(&r.idx);
                        if tbl.visible(1)
                            && ui
                                .add(
                                    egui::Button::new(if expanded { "v" } else { ">" })
                                        .frame(false),
                                )
                                .clicked()
                        {
                            expand_toggle = Some(r.idx);
                        }
                        if tbl.visible(2) {
                            let label = r.name.clone();
                            let text = egui::RichText::new(label).strong().color(if is_sel {
                                theme::ACCENT
                            } else if hidden {
                                theme::TEXT_DIM
                            } else {
                                theme::TEXT
                            });
                            let text = if is_sel { text.strong() } else { text };
                            let tip = feature_tip(&r, doc_prefix, unit);
                            // oracle：单击行只是**选中**（停留在 Features 页，底部写
                            // `Selected: …`），双击才跳到图谱/序列看这个特征
                            let resp = ui.selectable_label(is_sel, text).on_hover_text(tip);
                            if resp.clicked() {
                                pick_feature = Some(r.idx);
                            }
                            if resp.double_clicked() {
                                clicked = Some(r.idx);
                            }
                        }
                        if tbl.visible(3) {
                            // 面板里只有 `A .. B`（cross-origin 照原样写，SnapGene 不画 `*`）
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} .. {}",
                                    nfmt_plain(r.start + 1),
                                    nfmt_plain(r.end)
                                ))
                                .monospace()
                                .size(11.5),
                            );
                        }
                        if tbl.visible(4) {
                            ui.label(
                                egui::RichText::new(format!("{} {}", nfmt_plain(r.len), unit))
                                    .monospace()
                                    .size(11.5),
                            );
                        }
                        if tbl.visible(5) {
                            // 颜色块（SnapGene 的颜色列）
                            let (rect, _) = ui
                                .allocate_exact_size(egui::vec2(38.0, 24.0), egui::Sense::hover());
                            let sw =
                                egui::Rect::from_center_size(rect.center(), egui::vec2(22.0, 22.0));
                            ui.painter().rect_filled(sw, 1.5, r.color);
                            ui.painter().rect_stroke(
                                sw,
                                1.5,
                                egui::Stroke::new(0.8, theme::BORDER),
                                egui::StrokeKind::Inside,
                            );
                        }
                        if tbl.visible(6) {
                            // 方向：→ / ← / ├─（手绘，避免字形缺失）
                            let (rect, _) = ui
                                .allocate_exact_size(egui::vec2(48.0, 24.0), egui::Sense::hover());
                            let c = rect.center();
                            let p = ui.painter();
                            let col = theme::TEXT;
                            p.line_segment(
                                [c - egui::vec2(7.0, 0.0), c + egui::vec2(7.0, 0.0)],
                                egui::Stroke::new(1.2, col),
                            );
                            if r.directionality == 0 {
                                for x in [-7.0, 7.0] {
                                    p.line_segment(
                                        [c + egui::vec2(x, -3.5), c + egui::vec2(x, 3.5)],
                                        egui::Stroke::new(1.2, col),
                                    );
                                }
                            } else {
                                for (dir, enabled) in [
                                    (1.0, r.directionality == 1 || r.directionality == 3),
                                    (-1.0, r.directionality == 2 || r.directionality == 3),
                                ] {
                                    if enabled {
                                        p.add(egui::Shape::line(
                                            vec![
                                                c + egui::vec2(2.5 * dir, -3.5),
                                                c + egui::vec2(7.0 * dir, 0.0),
                                                c + egui::vec2(2.5 * dir, 3.5),
                                            ],
                                            egui::Stroke::new(1.2, col),
                                        ));
                                    }
                                }
                            }
                        }
                        if tbl.visible(7) {
                            ui.label(egui::RichText::new(&r.kind).size(11.5));
                        }
                        end_table_row(ui);
                        if expanded {
                            // 行内明细：每个片段一行（Location / Size 用与主行同一套列）
                            for (i, (a, b)) in r.segs.iter().enumerate() {
                                if tbl.visible(0) {
                                    ui.label("");
                                }
                                if tbl.visible(1) {
                                    ui.label("");
                                }
                                if tbl.visible(2) {
                                    detail_cell(ui, widths[2], format!("segment {}", i + 1));
                                }
                                if tbl.visible(3) {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "{} .. {}",
                                            nfmt_plain(*a + 1),
                                            nfmt_plain(*b)
                                        ))
                                        .monospace()
                                        .size(11.5),
                                    );
                                }
                                if tbl.visible(4) {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "{} {}",
                                            nfmt_plain(hx_render::span_len(*a, *b, len)),
                                            unit
                                        ))
                                        .monospace()
                                        .size(11.5),
                                    );
                                }
                                for c in 5..8 {
                                    if tbl.visible(c) {
                                        ui.label("");
                                    }
                                }
                                end_table_row(ui);
                            }
                        }
                    }
                });
        });
    app.tbl_features = tbl;
    if let Some(idx) = pick_feature {
        if !app.selected_features.insert(idx) {
            app.selected_features.remove(&idx);
        }
        // 单击即成为「当前特征」（Edit/Remove/颜色… 都作用于它），底部状态栏因此能写
        // `Selected: <名> (p3 .. p414 / 39 aa)`；不切换页签（SnapGene 行为）
        app.selected_feature = Some(idx);
        app.status.clear();
    }
    if let Some(idx) = expand_toggle {
        if !app.expanded_features.insert(idx) {
            app.expanded_features.remove(&idx);
        }
    }
    if let Some(idx) = toggle_hidden {
        if !app.hidden_features.insert(idx) {
            app.hidden_features.remove(&idx);
        }
        app.status = if app.hidden_features.contains(&idx) {
            "Feature hidden in views".into()
        } else {
            "Feature shown again".into()
        };
    }
    if let Some(idx) = clicked {
        app.center_on_feature(idx);
        app.tab = Tab::Map;
    }
}

// ---------------------------------------------------------------------------
// Primers
// ---------------------------------------------------------------------------

/// 引物面板。
///
/// 版式按真 SnapGene 8.2.3 的 Primers 页实测（1280×1024、默认缩放）：
///   * 页底色 `#F2F2F2`；表格白底，外框与行/列分隔线 1 px `#D6D6D6`；
///   * 工具条：`Find primers` 搜索框（白底、圆角 4、边框 `#A9A9A9`）、
///     `Actions ▾`（蓝字 `#0066EC`）；右端 `N Primers` + 排序图标；
///   * 表头 32 px：`[▢] [▸] Primer ⇧ | Length ⇅ | Color ⇅ | Binding site ⇅ |
///     Direction ⇅ | Tm ⇅`，列宽 399/73/99/100/99/166；
///   * 数据行 45 px（44 + 1 px 分隔线）：展开三角、复选框、粗体名字
///     `#212121`、`N-mer`、色块（无颜色时留空）、`a .. b`、方向箭头、`nn°C`；
///   * 展开后是 `#EBEBEB` 的详情块：`SEQUENCE`（+ `X% GC / N Da`）、
///     `DATE ADDED`、`MODIFICATION`（仅 5′ 磷酸化时）、`Show Structure` 链接。
pub fn primers(ui: &mut egui::Ui, app: &mut HelixApp) {
    struct Row {
        idx: usize,
        name: String,
        sequence: String,
        site: String,
        /// 供排序用的绑定位置（补零便于字符串比较）。
        sort_site: String,
        mer: usize,
        annealed: usize,
        forward: bool,
        strand: &'static str,
        tm: Option<f64>,
        tm_text: String,
        gc: f64,
        mass: f64,
        date_added: String,
        phosphorylated: bool,
    }

    let rows: Vec<Row> = app
        .doc
        .as_ref()
        .map(|doc| {
            doc.primers()
                .iter()
                .enumerate()
                .map(|(idx, p)| {
                    let site = match p.sites.first() {
                        Some(s) => format!(
                            "{} .. {}",
                            nfmt_plain(s.range.start + 1),
                            nfmt_plain(s.range.end)
                        ),
                        None => "—".into(),
                    };
                    let seq = if p.sequence.is_empty() {
                        "—".to_string()
                    } else {
                        p.sequence.clone()
                    };
                    let site_start = p.sites.first().map(|s| s.range.start).unwrap_or(0);
                    // SnapGene 的 Length 列 = 整条引物长度（含 5′ 尾巴），
                    // 退火碱基数只出现在详情与「结合位点」里。
                    let mer = p.sequence.len();
                    let annealed = p.sites.first().map(|s| s.annealed_bases).unwrap_or(0);
                    // 退火段 = 序列的 3′ 后缀
                    let ann_seq: String = if annealed > 0 && annealed <= p.sequence.len() {
                        p.sequence[p.sequence.len() - annealed..].to_string()
                    } else {
                        p.sequence.clone()
                    };
                    let tm = p
                        .sites
                        .first()
                        .and_then(|s| s.melting_temperature)
                        .or_else(|| {
                            let t = hx_algo::primer::tm(ann_seq.as_bytes());
                            t.is_finite().then_some(t)
                        });
                    let forward = p.sites.first().map(|s| s.strand != 1).unwrap_or(true);
                    Row {
                        idx,
                        name: p.name.clone(),
                        sequence: seq,
                        site,
                        sort_site: format!("{site_start:010}"),
                        mer,
                        annealed,
                        forward,
                        strand: if forward { "forward" } else { "reverse" },
                        tm,
                        tm_text: tm.map(|t| format!("{t:.0}°C")).unwrap_or_else(|| "—".into()),
                        gc: hx_algo::primer::gc_percent(p.sequence.as_bytes()),
                        mass: hx_algo::primer::ss_dna_mass(
                            p.sequence.as_bytes(),
                            p.phosphorylated,
                        ),
                        date_added: p.date_added.clone(),
                        phosphorylated: p.phosphorylated,
                    }
                })
                .collect()
        })
        .unwrap_or_default();

    // ---- 工具条 ----
    ui.add_space(10.0);
    ui.horizontal(|ui| {
        ui.add_space(10.0);
        egui::Frame::new()
            .fill(egui::Color32::WHITE)
            .stroke(egui::Stroke::new(1.0, theme::hex("#A9A9A9", theme::BORDER)))
            .corner_radius(4.0)
            .inner_margin(egui::Margin::symmetric(6, 2))
            .show(ui, |ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut app.primer_search)
                        .frame(egui::Frame::NONE)
                        .hint_text(egui::RichText::new(crate::i18n::text(ui.ctx(),"Find primers")).size(12.0))
                        .desired_width(120.0)
                        .text_color(theme::TEXT),
                );
            });
        ui.add_space(12.0);
        ui.menu_button(
            egui::RichText::new(crate::i18n::text(ui.ctx(),"Actions ▾")).color(theme::ACCENT).size(12.5),
            |ui| {
                ui.set_min_width(180.0);
                if ui.button(crate::i18n::text(ui.ctx(),"Add Primer…")).clicked() {
                    app.open_add_primer_dialog();
                    ui.close();
                }
                let has = app.selected_primer.is_some();
                if ui
                    .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Edit Primer…")))
                    .clicked()
                {
                    app.open_edit_primer_dialog();
                    ui.close();
                }
                if ui
                    .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Duplicate Primer")))
                    .clicked()
                {
                    app.duplicate_selected_primer();
                    ui.close();
                }
                if ui
                    .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Remove Primer")))
                    .clicked()
                {
                    app.remove_selected_primer();
                    ui.close();
                }
                ui.separator();
                if ui
                    .add_enabled(has, egui::Button::new(crate::i18n::text(ui.ctx(),"Analyze Selected Primers…")))
                    .clicked()
                {
                    app.analyze_selected_primer();
                    ui.close();
                }
                ui.separator();
                ui.checkbox(&mut app.show_primers, crate::i18n::text(ui.ctx(),"Show or Hide Primers"));
            },
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(12.0);
            // 排序方向切换（SnapGene 的 ↑↓ 图标）
            let asc = match app.tbl_primers.sort {
                Some((_, a)) => a,
                None => true,
            };
            let (rect, resp) = ui.allocate_exact_size(egui::vec2(22.0, 20.0), egui::Sense::click());
            let p = ui.painter();
            let c = rect.center();
            let up = if asc { theme::ACCENT } else { theme::hex("#4F94EE", theme::ACCENT) };
            let down = if asc { theme::hex("#4F94EE", theme::ACCENT) } else { theme::ACCENT };
            p.add(egui::Shape::convex_polygon(
                vec![
                    egui::Pos2::new(c.x - 5.0, c.y - 6.0),
                    egui::Pos2::new(c.x + 5.0, c.y - 6.0),
                    egui::Pos2::new(c.x, c.y - 1.0),
                ],
                up,
                egui::Stroke::NONE,
            ));
            p.add(egui::Shape::convex_polygon(
                vec![
                    egui::Pos2::new(c.x - 5.0, c.y + 6.0),
                    egui::Pos2::new(c.x + 5.0, c.y + 6.0),
                    egui::Pos2::new(c.x, c.y + 1.0),
                ],
                down,
                egui::Stroke::NONE,
            ));
            if resp.clicked() {
                let cur = app.tbl_primers.sort.map(|(c, _)| c).unwrap_or(0);
                app.tbl_primers.sort = Some((cur, !asc));
            }
            ui.label(
                egui::RichText::new(format!(
                    "{} Primer{}",
                    rows.len(),
                    if rows.len() == 1 { "" } else { "s" }
                ))
                .size(12.0)
                .color(theme::TEXT_DIM),
            );
        });
    });
    ui.add_space(10.0);

    // ---- 过滤 / 排序 ----
    let needle = app.primer_search.trim().to_ascii_uppercase();
    let mut sorted: Vec<&Row> = rows
        .iter()
        .filter(|r| {
            needle.is_empty()
                || r.name.to_ascii_uppercase().contains(&needle)
                || r.sequence.to_ascii_uppercase().contains(&needle)
        })
        .collect();
    if let Some((ci, asc)) = app.tbl_primers.sort {
        sorted.sort_by(|a, b| {
            let key = |r: &Row| -> String {
                match ci {
                    0 => r.name.to_ascii_lowercase(),
                    1 => format!("{:010}", r.mer),
                    2 => String::new(),
                    3 => r.sort_site.clone(),
                    4 => r.strand.to_string(),
                    5 => format!("{:08.2}", r.tm.unwrap_or(f64::MIN)),
                    _ => String::new(),
                }
            };
            let (ka, kb) = (key(a), key(b));
            if asc {
                ka.cmp(&kb)
            } else {
                kb.cmp(&ka)
            }
        });
    }

    if rows.is_empty() {
        ui.horizontal(|ui| {
            ui.add_space(10.0);
            ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"No primers in this file.")).color(theme::TEXT_DIM));
        });
        return;
    }

    // ---- 列宽（SnapGene 实测：399 / 73 / 99 / 100 / 99 / 剩余）----
    let total = (ui.available_width() - 12.0).max(640.0);
    let w_len = 73.0;
    let w_col = 99.0;
    let w_site = 100.0;
    let w_dir = 99.0;
    let w_tm = 120.0;
    let w_primer = (total - w_len - w_col - w_site - w_dir - w_tm).max(180.0);
    let cols = [w_primer, w_len, w_col, w_site, w_dir, w_tm];

    let headers = ["Primer", "Length", "Color", "Binding site", "Direction", "Tm"];
    let mut clicked: Option<usize> = None;
    let mut check_toggle: Option<usize> = None;
    let mut expand_toggle: Option<usize> = None;
    let selected = app.selected_primer;

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.add_space(2.0);
            let table_w = total.min(ui.available_width());
            egui::Frame::new()
                .fill(egui::Color32::WHITE)
                .stroke(egui::Stroke::new(1.0, theme::hex("#D6D6D6", theme::BORDER)))
                .show(ui, |ui| {
                    ui.set_width(table_w);
                    ui.spacing_mut().item_spacing = egui::vec2(0.0, 0.0);

                    // 表头
                    let (hrect, hresp) = ui.allocate_exact_size(
                        egui::vec2(table_w, 31.0),
                        egui::Sense::click(),
                    );
                    let p = ui.painter();
                    let mut x = hrect.left();
                    let mut hsort: Option<usize> = None;
                    for (ci, name) in headers.iter().enumerate() {
                        let w = cols[ci];
                        let cell = egui::Rect::from_min_size(
                            egui::Pos2::new(x, hrect.top()),
                            egui::vec2(w, hrect.height()),
                        );
                        if ci > 0 {
                            p.line_segment(
                                [
                                    egui::Pos2::new(cell.left(), cell.top() + 2.0),
                                    egui::Pos2::new(cell.left(), cell.bottom() - 2.0),
                                ],
                                egui::Stroke::new(1.0, theme::hex("#D6D6D6", theme::BORDER)),
                            );
                        }
                        // 表头文字 + 排序箭头
                        let active = app.tbl_primers.sort.map(|(c, _)| c) == Some(ci);
                        let asc = app.tbl_primers.sort.map(|(_, a)| a).unwrap_or(true);
                        let glyph = if active {
                            if asc {
                                " \u{2191}"
                            } else {
                                " \u{2193}"
                            }
                        } else {
                            " \u{21C5}"
                        };
                        let text_x = if ci == 0 { cell.left() + 49.0 } else { cell.left() + 10.0 };
                        p.text(
                            egui::Pos2::new(text_x, cell.center().y),
                            egui::Align2::LEFT_CENTER,
                            crate::i18n::text(ui.ctx(),name),
                            egui::FontId::proportional(12.5),
                            theme::hex("#545454", theme::TEXT),
                        );
                        let gw = p
                            .layout_no_wrap(
                                name.to_string(),
                                egui::FontId::proportional(12.5),
                                theme::TEXT,
                            )
                            .size()
                            .x;
                        p.text(
                            egui::Pos2::new(text_x + gw + 4.0, cell.center().y),
                            egui::Align2::LEFT_CENTER,
                            glyph,
                            egui::FontId::proportional(10.5),
                            if active {
                                theme::ACCENT
                            } else {
                                theme::hex("#9A9A9A", theme::TEXT_DIM)
                            },
                        );
                        if hresp.clicked() && hresp.interact_pointer_pos().map(|q| cell.contains(q)).unwrap_or(false)
                        {
                            hsort = Some(ci);
                        }
                        x += w;
                    }
                    // 全选复选框（SnapGene 的表头也有一个，位置与行内复选框同列）
                    let all_checked = !rows.is_empty()
                        && rows.iter().all(|r| app.primer_checked.contains(&r.idx));
                    let hcb = egui::Rect::from_center_size(
                        egui::Pos2::new(hrect.left() + 33.0, hrect.center().y),
                        egui::vec2(14.0, 14.0),
                    );
                    if hresp.clicked()
                        && hresp
                            .interact_pointer_pos()
                            .map(|q| hcb.expand(4.0).contains(q))
                            .unwrap_or(false)
                    {
                        if all_checked {
                            app.primer_checked.clear();
                        } else {
                            app.primer_checked.extend(rows.iter().map(|r| r.idx));
                        }
                        hsort = None;
                    }
                    if let Some(ci) = hsort {
                        let cur = app.tbl_primers.sort;
                        let asc = match cur {
                            Some((c, a)) if c == ci => !a,
                            _ => true,
                        };
                        app.tbl_primers.sort = Some((ci, asc));
                    }
                    p.rect_filled(hcb, 3.0, egui::Color32::WHITE);
                    p.rect_stroke(
                        hcb,
                        3.0,
                        egui::Stroke::new(1.2, theme::hex("#B0BEC5", theme::BORDER)),
                        egui::StrokeKind::Inside,
                    );
                    if all_checked {
                        p.rect_filled(hcb.shrink(3.0), 1.5, theme::ACCENT);
                    }
                    p.line_segment(
                        [
                            egui::Pos2::new(hrect.left(), hrect.bottom() - 0.5),
                            egui::Pos2::new(hrect.right(), hrect.bottom() - 0.5),
                        ],
                        egui::Stroke::new(1.0, theme::hex("#D6D6D6", theme::BORDER)),
                    );

                    // 数据行
                    for r in sorted {
                        let expanded = app.primer_expanded.contains(&r.idx);
                        let is_sel = selected == Some(r.idx);
                        let row_h = 44.0;
                        let (rect, resp) = ui.allocate_exact_size(
                            egui::vec2(table_w, row_h),
                            egui::Sense::click(),
                        );
                        let p = ui.painter();
                        if resp.clicked() {
                            clicked = Some(r.idx);
                        }
                        if resp.double_clicked() {
                            clicked = Some(r.idx);
                            app.selected_primer = Some(r.idx);
                            app.open_edit_primer_dialog();
                        }
                        // 列分隔线
                        let mut x = rect.left();
                        for (ci, w) in cols.iter().enumerate() {
                            if ci > 0 {
                                p.line_segment(
                                    [
                                        egui::Pos2::new(x, rect.top()),
                                        egui::Pos2::new(x, rect.bottom()),
                                    ],
                                    egui::Stroke::new(1.0, theme::hex("#D6D6D6", theme::BORDER)),
                                );
                            }
                            x += w;
                        }
                        // 展开三角
                        let chev = egui::Rect::from_min_size(
                            egui::Pos2::new(rect.left() + 6.0, rect.center().y - 7.0),
                            egui::vec2(18.0, 14.0),
                        );
                        if resp.clicked()
                            && resp
                                .interact_pointer_pos()
                                .map(|q| chev.expand(3.0).contains(q))
                                .unwrap_or(false)
                        {
                            expand_toggle = Some(r.idx);
                            clicked = None;
                        }
                        let cc = chev.center();
                        let tri = if expanded {
                            vec![
                                egui::Pos2::new(cc.x - 4.5, cc.y - 1.5),
                                egui::Pos2::new(cc.x + 4.5, cc.y - 1.5),
                                egui::Pos2::new(cc.x, cc.y + 3.5),
                            ]
                        } else {
                            vec![
                                egui::Pos2::new(cc.x - 1.5, cc.y - 4.5),
                                egui::Pos2::new(cc.x + 3.5, cc.y),
                                egui::Pos2::new(cc.x - 1.5, cc.y + 4.5),
                            ]
                        };
                        p.add(egui::Shape::convex_polygon(
                            tri,
                            theme::hex("#212121", theme::TEXT),
                            egui::Stroke::NONE,
                        ));
                        // 复选框
                        let cb = egui::Rect::from_center_size(
                            egui::Pos2::new(rect.left() + 33.0, rect.center().y),
                            egui::vec2(14.0, 14.0),
                        );
                        let marked = app.primer_checked.contains(&r.idx);
                        p.rect_filled(cb, 3.0, egui::Color32::WHITE);
                        p.rect_stroke(
                            cb,
                            3.0,
                            egui::Stroke::new(1.2, theme::hex("#B0BEC5", theme::BORDER)),
                            egui::StrokeKind::Inside,
                        );
                        if marked {
                            let inner = cb.shrink(3.0);
                            p.rect_filled(inner, 1.5, theme::ACCENT);
                        }
                        if resp.clicked()
                            && resp
                                .interact_pointer_pos()
                                .map(|q| cb.expand(4.0).contains(q))
                                .unwrap_or(false)
                        {
                            check_toggle = Some(r.idx);
                        }

                        // 各列内容
                        let name_x = rect.left() + 49.0;
                        p.text(
                            egui::Pos2::new(name_x, rect.center().y),
                            egui::Align2::LEFT_CENTER,
                            &r.name,
                            egui::FontId::proportional(13.0),
                            theme::hex("#212121", theme::TEXT),
                        );
                        let cx = |ci: usize| -> f32 {
                            rect.left() + cols[..ci].iter().sum::<f32>() + 10.0
                        };
                        p.text(
                            egui::Pos2::new(cx(1), rect.center().y),
                            egui::Align2::LEFT_CENTER,
                            format!("{}-mer", r.mer),
                            egui::FontId::proportional(12.5),
                            theme::hex("#323232", theme::TEXT),
                        );
                        // Color：没有显式颜色时整格留空（与 SnapGene 一致）
                        let _ = cx(2);
                        p.text(
                            egui::Pos2::new(cx(3), rect.center().y),
                            egui::Align2::LEFT_CENTER,
                            &r.site,
                            egui::FontId::proportional(12.5),
                            theme::hex("#323232", theme::TEXT),
                        );
                        // 方向箭头：正向右上、反向左下（SnapGene 的斜箭头）
                        let dc = egui::Pos2::new(cx(4) + 12.0, rect.center().y);
                        let arrow = if r.forward {
                            [
                                egui::Pos2::new(dc.x - 8.0, dc.y + 4.0),
                                egui::Pos2::new(dc.x + 7.0, dc.y - 4.0),
                            ]
                        } else {
                            [
                                egui::Pos2::new(dc.x + 8.0, dc.y - 4.0),
                                egui::Pos2::new(dc.x - 7.0, dc.y + 4.0),
                            ]
                        };
                        p.add(egui::Shape::line_segment(
                            arrow,
                            egui::Stroke::new(1.6, theme::hex("#212121", theme::TEXT)),
                        ));
                        let (tip, a, b) = if r.forward {
                            (
                                arrow[1],
                                egui::Pos2::new(arrow[1].x - 5.0, arrow[1].y - 1.0),
                                egui::Pos2::new(arrow[1].x - 1.0, arrow[1].y + 5.0),
                            )
                        } else {
                            (
                                arrow[1],
                                egui::Pos2::new(arrow[1].x + 5.0, arrow[1].y + 1.0),
                                egui::Pos2::new(arrow[1].x + 1.0, arrow[1].y - 5.0),
                            )
                        };
                        p.add(egui::Shape::convex_polygon(
                            vec![tip, a, b],
                            theme::hex("#212121", theme::TEXT),
                            egui::Stroke::NONE,
                        ));
                        p.text(
                            egui::Pos2::new(cx(5), rect.center().y),
                            egui::Align2::LEFT_CENTER,
                            &r.tm_text,
                            egui::FontId::proportional(12.5),
                            theme::hex("#323232", theme::TEXT),
                        );
                        // 选中：整行 1 px 蓝框（SnapGene 的焦点框）
                        if is_sel {
                            p.rect_stroke(
                                rect.shrink(0.5),
                                0.0,
                                egui::Stroke::new(1.0, theme::ACCENT),
                                egui::StrokeKind::Inside,
                            );
                        }
                        p.line_segment(
                            [
                                egui::Pos2::new(rect.left(), rect.bottom() - 0.5),
                                egui::Pos2::new(rect.right(), rect.bottom() - 0.5),
                            ],
                            egui::Stroke::new(1.0, theme::hex("#D6D6D6", theme::BORDER)),
                        );

                        // ---- 展开详情 ----
                        if expanded {
                            let detail_rows = 1
                                + 1
                                + usize::from(r.phosphorylated)
                                + 1
                                + usize::from(!r.date_added.is_empty());
                            let dh = 16.0 + detail_rows as f32 * 29.0 + 22.0;
                            let (drect, _) = ui.allocate_exact_size(
                                egui::vec2(table_w, dh),
                                egui::Sense::hover(),
                            );
                            let p = ui.painter();
                            p.rect_filled(drect, 0.0, theme::hex("#EBEBEB", theme::PANEL));
                            let label_x = drect.left() + 82.0;
                            let value_x = drect.left() + 196.0;
                            let mut y = drect.top() + 16.0;
                            let label_font = egui::FontId::proportional(10.5);
                            let value_font = egui::FontId::proportional(12.5);
                            let mono = egui::FontId::monospace(12.5);
                            let lc = theme::hex("#6E6E6E", theme::TEXT_DIM);
                            let vc = theme::hex("#2F2F2F", theme::TEXT);
                            // SEQUENCE
                            p.text(
                                egui::Pos2::new(label_x, y + 7.0),
                                egui::Align2::LEFT_CENTER,
                                crate::i18n::text(ui.ctx(),"SEQUENCE"),
                                label_font.clone(),
                                lc,
                            );
                            p.text(
                                egui::Pos2::new(value_x, y + 7.0),
                                egui::Align2::LEFT_CENTER,
                                &r.sequence,
                                mono.clone(),
                                vc,
                            );
                            y += 29.0;
                            p.text(
                                egui::Pos2::new(value_x, y + 7.0),
                                egui::Align2::LEFT_CENTER,
                                format!(
                                    "{:.0}% GC / {} Da",
                                    r.gc,
                                    nfmt_mass(r.mass)
                                ),
                                value_font.clone(),
                                theme::hex("#6E6E6E", theme::TEXT_DIM),
                            );
                            y += 29.0;
                            if !r.date_added.is_empty() {
                                p.text(
                                    egui::Pos2::new(label_x, y + 7.0),
                                    egui::Align2::LEFT_CENTER,
                                    crate::i18n::text(ui.ctx(),"DATE ADDED"),
                                    label_font.clone(),
                                    lc,
                                );
                                p.text(
                                    egui::Pos2::new(value_x, y + 7.0),
                                    egui::Align2::LEFT_CENTER,
                                    fmt_date_added(&r.date_added),
                                    value_font.clone(),
                                    vc,
                                );
                                y += 29.0;
                            }
                            if r.phosphorylated {
                                p.text(
                                    egui::Pos2::new(label_x, y + 7.0),
                                    egui::Align2::LEFT_CENTER,
                                    crate::i18n::text(ui.ctx(),"MODIFICATION"),
                                    label_font.clone(),
                                    lc,
                                );
                                p.text(
                                    egui::Pos2::new(value_x, y + 7.0),
                                    egui::Align2::LEFT_CENTER,
                                    "5' phosphorylated",
                                    value_font.clone(),
                                    vc,
                                );
                                y += 29.0;
                            }
                            p.text(
                                egui::Pos2::new(value_x, y + 7.0),
                                egui::Align2::LEFT_CENTER,
                                crate::i18n::text(ui.ctx(),"Show Structure"),
                                value_font.clone(),
                                theme::ACCENT,
                            );
                            let _ = (r.annealed, r.tm);
                        }
                    }
                });
            ui.add_space(10.0);
        });

    if let Some(i) = check_toggle {
        if !app.primer_checked.remove(&i) {
            app.primer_checked.insert(i);
        }
    }
    if let Some(i) = expand_toggle {
        if !app.primer_expanded.remove(&i) {
            app.primer_expanded.insert(i);
        }
    }
    if let Some(i) = clicked {
        app.selected_primer = Some(i);
        // 选中引物时把它的结合位点设为当前选段（两个视图都会高亮）
        let info = app.doc.as_ref().and_then(|d| {
            d.primers()
                .get(i)
                .map(|p| (p.name.clone(), p.sites.first().map(|s| s.range)))
        });
        if let Some((name, site)) = info {
            if let Some(range) = site {
                app.selection = Some((range.start, range.end));
            }
            let n = app.doc.as_ref().map(|d| d.primers().len()).unwrap_or(0);
            app.status = format!("Selected primer {name}  —  {n} primers");
        }
    }
}

/// `dateAdded`（`YYYY-MM-DDTHH:MM:SSZ`）→ 面板里的日期文案。
///
/// SnapGene 显示的是**本地化长日期**（中文界面下是「2026年10月19日星期一」）；
/// 我们的界面是英文，因此用 `Weekday, Month D, YYYY` 的英文写法。
fn fmt_date_added(s: &str) -> String {
    let date = s.split('T').next().unwrap_or(s);
    let mut it = date.split('-');
    let (Some(y), Some(m), Some(d)) = (it.next(), it.next(), it.next()) else {
        return date.to_string();
    };
    let (Ok(y), Ok(m), Ok(d)) = (y.parse::<i32>(), m.parse::<u32>(), d.parse::<u32>()) else {
        return date.to_string();
    };
    let months = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ];
    let weekday = weekday_name(y, m, d);
    let mon = months.get((m as usize).saturating_sub(1)).copied().unwrap_or("");
    format!("{weekday}, {mon} {d}, {y}")
}

/// Zeller 同余（Gregorian）取星期名。
fn weekday_name(y: i32, m: u32, d: u32) -> &'static str {
    let (mut yy, mut mm) = (y, m);
    if mm < 3 {
        mm += 12;
        yy -= 1;
    }
    let k = yy % 100;
    let j = yy / 100;
    let h = (d as i32 + (13 * (mm as i32 + 1)) / 5 + k + k / 4 + j / 4 + 5 * j).rem_euclid(7);
    // h: 0=Saturday, 1=Sunday, 2=Monday …
    match h {
        0 => "Saturday",
        1 => "Sunday",
        2 => "Monday",
        3 => "Tuesday",
        4 => "Wednesday",
        5 => "Thursday",
        _ => "Friday",
    }
}

/// 分子量文案：与 SnapGene 一致，**5 位以上才加千位分隔**（`9939.6` / `10,812.1`）。
fn nfmt_mass(v: f64) -> String {
    let s = format!("{v:.1}");
    let (int, frac) = s.split_once('.').unwrap_or((s.as_str(), "0"));
    if int.len() >= 5 {
        let neg = int.starts_with('-');
        let digits: Vec<char> = int.trim_start_matches('-').chars().collect();
        let mut grouped = String::new();
        for (i, c) in digits.iter().enumerate() {
            if i > 0 && (digits.len() - i) % 3 == 0 {
                grouped.push(',');
            }
            grouped.push(*c);
        }
        format!("{}{grouped}.{frac}", if neg { "-" } else { "" })
    } else {
        s
    }
}

// ---------------------------------------------------------------------------
// History
// ---------------------------------------------------------------------------

/// 从 block 7 载荷里扫出节点（极简：只取 `name`/`seqLen`/`operation` 与嵌套深度）。
fn scan_history(xml: &str) -> Vec<(usize, String, usize, String)> {
    fn attr(tag: &str, key: &str) -> Option<String> {
        let pat = format!("{key}=\"");
        let i = tag.find(&pat)? + pat.len();
        let rest = &tag[i..];
        let j = rest.find('"')?;
        Some(rest[..j].to_string())
    }

    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut idx = 0usize;
    while idx < xml.len() {
        let open = xml[idx..].find("<Node ").map(|i| idx + i);
        let close = xml[idx..].find("</Node>").map(|i| idx + i);
        let open_first = match (open, close) {
            (Some(o), Some(c)) => o < c,
            (Some(_), None) => true,
            _ => false,
        };
        if open_first {
            let o = open.unwrap_or(idx);
            let end = xml[o..].find('>').map(|e| o + e + 1).unwrap_or(xml.len());
            let tag = &xml[o..end];
            let name = attr(tag, "name").unwrap_or_else(|| "(unnamed)".into());
            let seq_len = attr(tag, "seqLen")
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(0);
            let op = attr(tag, "operation").unwrap_or_default();
            out.push((depth, name, seq_len, op));
            depth += 1;
            idx = end;
        } else if let Some(c) = close {
            depth = depth.saturating_sub(1);
            idx = c + "</Node>".len();
        } else {
            break;
        }
    }
    out
}

/// `operation` 属性 → 显示名（SnapGene 的历史条目文案）。
fn op_label(op: &str) -> String {
    match op {
        "insertFragment" => "Insert fragment".into(),
        "replace" => "Replace".into(),
        "insert" | "insertBases" => "Insert".into(),
        "delete" => "Delete".into(),
        "" => String::new(),
        other => other.to_string(),
    }
}

pub fn history(ui: &mut egui::Ui, app: &mut HelixApp) {
    crate::history_view::show(ui, app);
}

#[allow(dead_code)]
fn legacy_history(ui: &mut egui::Ui, app: &mut HelixApp) {
    let xml: Option<String> = app.doc.as_ref().and_then(|doc| {
        doc.session
            .file()
            .blocks
            .iter()
            .find(|b| b.block_type == 7)
            .map(|b| String::from_utf8_lossy(&b.payload).to_string())
    });
    let nodes = xml.as_deref().map(scan_history).unwrap_or_default();
    page_header(ui, "History", &format!("{} entries", nodes.len()));

    if nodes.is_empty() {
        ui.horizontal(|ui| {
            ui.add_space(10.0);
            ui.label(
                egui::RichText::new(crate::i18n::text(ui.ctx(),"No recorded history."))
                    .size(12.5)
                    .color(theme::TEXT_DIM),
            );
        });
        return;
    }

    let mut hstate = crate::app::TableState::default();
    let widths = column_widths(
        (ui.available_width() - 12.0).max(0.0),
        &hstate,
        [260.0, 120.0, 240.0],
        0,
    );
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new("history_table")
                .num_columns(3)
                .striped(false)
                .spacing([COLUMN_GAP, 0.0])
                .min_col_width(24.0)
                .min_row_height(36.0)
                .show(ui, |ui| {
                    table_header(
                        ui,
                        &mut hstate,
                        &[("Name", false), ("Length", false), ("Operation", false)],
                        &widths,
                    );
                    for (depth, name, len, op) in &nodes {
                        ui.label(format!("{}{}", "    ".repeat(*depth), name));
                        ui.label(
                            egui::RichText::new(crate::ui_format!(ui.ctx(),"{} bp", nfmt(*len)))
                                .monospace()
                                .color(theme::TEXT_DIM),
                        );
                        ui.label(crate::i18n::text(ui.ctx(),&op_label(op)));
                        end_table_row(ui);
                    }
                });
        });
}
