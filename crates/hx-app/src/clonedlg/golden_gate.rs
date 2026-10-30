//! Golden Gate workflow and layout, based on the supplied and live SG 8.2.3 pages.
//! All modes prepare physical double-stranded fragments before assembly.
use super::*;

pub(super) fn active(dlg: &CloningDlg) -> bool {
    matches!(dlg.method, Method::GoldenGate | Method::GoldenGateMulti)
}

fn enzyme(dlg: &CloningDlg) -> Result<&'static hx_enzymes::Enzyme, String> {
    dlg.gg_enzymes.get(dlg.gg_enzyme).and_then(|n| hx_enzymes::by_name(n))
        .ok_or_else(|| "The selected Golden Gate enzyme is not available in the enzyme database.".into())
}

fn pcr_source(src: &Source, pick: &Pick) -> Result<Source, String> {
    let (a, b) = pick.region.ok_or("Please select the region to be amplified.")?;
    if a >= src.seq.len() || b > src.seq.len() || (!src.circular && a >= b) {
        return Err("Please select a valid amplified region.".into());
    }
    let slice = src.slice(a, b, pick.reversed);
    if slice.seq.len() < 12 { return Err("The amplified region is too short for PCR primers.".into()); }
    let forward = pick.primer1.trim().to_ascii_uppercase().into_bytes();
    let reverse = pick.primer2.trim().to_ascii_uppercase().into_bytes();
    if forward.is_empty() || reverse.is_empty() { return Err("Please specify two PCR primers.".into()); }
    let rc = hx_algo::reverse_complement(&slice.seq);
    let anneal = |primer: &[u8], template: &[u8]| (12..=primer.len().min(template.len())).rev()
        .find(|&n| primer[primer.len() - n..].eq_ignore_ascii_case(&template[..n]));
    let f = anneal(&forward, &slice.seq).ok_or("Primer 1 does not anneal to the selected region.")?;
    let r = anneal(&reverse, &rc).ok_or("Primer 2 does not anneal to the selected region.")?;
    let mut pcr = slice;
    pcr.prepend(&forward[..forward.len() - f]);
    pcr.append(&hx_algo::reverse_complement(&reverse[..reverse.len() - r]));
    Ok(Source::new(src.label.clone(), pcr.seq, false).with_features(pcr.features))
}

pub(super) fn prepared(dlg: &CloningDlg, pick: &Pick, vector: bool) -> Result<Source, String> {
    let source = dlg.pick_source(pick).ok_or(if vector { "Please specify a linearized Vector." } else { "Please choose a source file for the Fragment." })?;
    let e = enzyme(dlg)?;
    match pick.gg_mode {
        GoldenGateMode::Direct => {
            if source.circular { return Err("Please choose a linear sequence to use directly.".into()); }
            Ok(if pick.reversed { source.flipped() } else { source.clone() })
        }
        GoldenGateMode::OtherEnzymes => restriction_fragment(source, pick, vector, true),
        GoldenGateMode::PcrDigest | GoldenGateMode::Digest => {
            let pcr;
            if vector && pick.gg_mode == GoldenGateMode::PcrDigest && pick.region.is_none() {
                return Err("Please specify a linearized Vector.".into());
            }
            let src = if pick.gg_mode == GoldenGateMode::PcrDigest { pcr = pcr_source(source, pick)?; &pcr } else { source };
            let candidates = cloning::golden_gate::fragments(src, e);
            if candidates.is_empty() { return Err(format!("No usable fragment bounded by two {} cuts. Include two flanking sites and remove internal sites.", e.name)); }
            let chosen = if let Some((a, b)) = pick.cut_pair {
                candidates.iter().find(|(x, y, _)| *x == a && *y == b)
                    .ok_or("Please select a fragment from the digest list.")?
            } else if candidates.len() == 1 { &candidates[0] }
            else { return Err("More than one usable digest fragment. Please select one explicitly.".into()); };
            let fragment = chosen.2.clone();
            Ok(if pick.gg_mode == GoldenGateMode::Digest && pick.reversed { fragment.flipped() } else { fragment })
        }
    }
}

pub(super) fn build(dlg: &CloningDlg) -> Result<Product, String> {
    let mut parts = vec![prepared(dlg, &dlg.vector, true)?];
    for pick in dlg.fragments.iter().take(if dlg.method.is_multiple() { dlg.fragments.len() } else { 1 }) {
        parts.push(prepared(dlg, pick, false)?);
    }
    cloning::golden_gate::assemble(&parts, enzyme(dlg)?)
}

pub(super) fn header(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Golden Gate Enzyme:")).strong());
        let old = dlg.gg_enzyme;
        let names = dlg.gg_enzymes.clone();
        let selected = names
            .get(dlg.gg_enzyme.min(names.len().saturating_sub(1)))
            .copied()
            .unwrap_or("BsaI");
        egui::ComboBox::from_id_salt("gg_global_enzyme")
            .selected_text(selected)
            .width(140.0).show_ui(ui, |ui| {
                for (i, name) in names.iter().enumerate() { ui.selectable_value(&mut dlg.gg_enzyme, i, *name); }
            });
        if old != dlg.gg_enzyme {
            dlg.vector.cut_pair = None;
            for p in &mut dlg.fragments { p.cut_pair = None; }
        }
        if let Ok(e) = enzyme(dlg) {
            let extra = (e.cut as i32 - e.site_len() as i32).max(0) as usize;
            let top = format!("{}{}▼", e.site.chars().map(|c| c.to_string()).collect::<Vec<_>>().join(" "), " N".repeat(extra));
            let bottom = format!("{}{} (N){}▲", e.site.bytes().map(|b| hx_algo::complement(b) as char).map(|c| c.to_string()).collect::<Vec<_>>().join(" "),
                " N".repeat(extra), e.site_len() as i32 - e.reverse_cut as i32 - e.cut as i32);
            ui.label(RichText::new(format!("{top}\n{bottom}")).monospace().size(11.0));
        }
    });
}

pub(super) fn panel(ui: &mut egui::Ui, dlg: &mut CloningDlg, vector: bool) {
    if !vector && dlg.method.is_multiple() {
        ui.horizontal(|ui| {
            ui.label(crate::i18n::text(ui.ctx(),"Fragment:"));
            for i in 0..dlg.fragments.len() { ui.selectable_value(&mut dlg.frag_index, i, (i + 1).to_string()); }
            if ui.small_button("+").clicked() && dlg.fragments.len() < 8 { dlg.fragments.push(Pick::default()); }
            if ui.small_button("−").clicked() && dlg.fragments.len() > 2 { dlg.fragments.pop(); dlg.frag_index = 0; }
        });
    }
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let pick = if vector { &mut dlg.vector } else { &mut dlg.fragments[idx] };
    ui.label(crate::i18n::text(ui.ctx(),if vector { "Vector:" } else { "Source of Fragment:" }));
    let old_doc = pick.doc;
    egui::ComboBox::from_id_salt(("gg_source", vector, idx))
        .selected_text(pick.doc.and_then(|i| dlg.docs.get(i)).map(|d| d.name.clone()).unwrap_or_else(|| crate::i18n::text(ui.ctx(),"< Please choose >").into()))
        .width(ui.available_width()).show_ui(ui, |ui| {
            ui.selectable_value(&mut pick.doc, None, crate::i18n::text(ui.ctx(),"< Please choose >"));
            for (i, doc) in dlg.docs.iter().enumerate() { ui.selectable_value(&mut pick.doc, Some(i), &doc.name); }
        });
    if old_doc != pick.doc {
        pick.region = None; pick.cut_pair = None;
        pick.gg_mode = pick.doc.and_then(|i| dlg.docs.get(i))
            .zip(dlg.gg_enzymes.get(dlg.gg_enzyme).and_then(|name| hx_enzymes::by_name(name)))
            .filter(|(d, e)| !cloning::golden_gate::fragments(&d.src, e).is_empty())
            .map(|_| GoldenGateMode::Digest).unwrap_or(GoldenGateMode::PcrDigest);
    }
    let source = pick.doc.and_then(|i| dlg.docs.get(i)).map(|d| &d.src);
    let e = dlg.gg_enzymes.get(dlg.gg_enzyme).and_then(|n| hx_enzymes::by_name(n));
    let name = e.map(|e| e.name).unwrap_or("enzyme");
    let has_source = source.is_some();
    let can_digest = source.zip(e).is_some_and(|(s, e)| !cloning::golden_gate::fragments(s, e).is_empty());
    let can_direct = source.is_some_and(|s| !s.circular);
    let old = pick.gg_mode;
    ui.add_space(9.0);
    ui.add_enabled_ui(has_source, |ui| { ui.radio_value(&mut pick.gg_mode, GoldenGateMode::PcrDigest, crate::ui_format!(ui.ctx(),"Run PCR and then digest with {}",name)); });
    ui.add_enabled_ui(can_digest, |ui| { ui.radio_value(&mut pick.gg_mode, GoldenGateMode::Digest, crate::ui_format!(ui.ctx(),"Digest with {}",name)); });
    if vector { ui.add_enabled_ui(has_source, |ui| { ui.radio_value(&mut pick.gg_mode, GoldenGateMode::OtherEnzymes, crate::i18n::text(ui.ctx(),"Digest with other enzymes")); }); }
    ui.add_enabled_ui(can_direct, |ui| { ui.radio_value(&mut pick.gg_mode, GoldenGateMode::Direct,
        crate::i18n::text(ui.ctx(),if vector { "Use directly as the linearized Vector" } else { "Use directly as the Fragment" })); });
    if old != pick.gg_mode { pick.cut_pair = None; }
    let mode = pick.gg_mode;
    let direct_size = if mode == GoldenGateMode::Direct { source.map(Source::display_len) } else { None };
    ui.add_space(10.0);
    if mode == GoldenGateMode::OtherEnzymes {
        let filter = String::new();
        ui.horizontal(|ui| { enzyme_combo(ui, pick, 0, &filter, vector, idx); enzyme_combo(ui, pick, 1, &filter, vector, idx); });
    }
    let table_h = (ui.available_height() - 185.0).max(120.0);
    egui::Frame::NONE.stroke(Stroke::new(1.0, theme::BORDER)).show(ui, |ui| {
        ui.set_width(ui.available_width());
        ui.horizontal(|ui| {
            ui.add_sized([ui.available_width() * 0.65, 24.0], egui::Label::new(if mode == GoldenGateMode::PcrDigest { "Amplified Region" } else { "Fragment" }));
            ui.label(crate::i18n::text(ui.ctx(),"Fragment Size (bp)"));
        });
        ui.separator();
        egui::ScrollArea::vertical().max_height(table_h).min_scrolled_height(table_h).show(ui, |ui| {
            if let Some(source) = source {
                match mode {
                    GoldenGateMode::PcrDigest => {
                        if let Some((a, b)) = pick.region {
                            ui.label(format!("{} .. {}       {}", a + 1, b, source.slice(a, b, false).seq.len()));
                        }
                        ui.add_space((table_h * 0.23).max(10.0));
                        ui.label(crate::i18n::text(ui.ctx(),"To choose PCR primers automatically:"));
                        ui.add_space(10.0);
                        ui.label(crate::i18n::text(ui.ctx(),"1. In the Vector tab, select the site of insertion or the region to be replaced."));
                        ui.add_space(10.0);
                        ui.label(crate::i18n::text(ui.ctx(),"2. In the Fragment tab, select the region to be amplified and choose the orientation."));
                        ui.add_space(10.0);
                        ui.label(crate::i18n::text(ui.ctx(),"3. In the Product tab, press ‘Choose PCR Primers…’."));
                    }
                    GoldenGateMode::Digest => if let Some(e) = e {
                        for (a, b, fragment) in cloning::golden_gate::fragments(source, e) {
                            let pair = Some((a, b));
                            ui.selectable_value(&mut pick.cut_pair, pair, format!("{} – {}          {}", a.top, b.top, fragment.display_len()));
                        }
                    },
            GoldenGateMode::Direct => { ui.label(crate::ui_format!(ui.ctx(),"Start – End          {}", source.display_len())); },
                    GoldenGateMode::OtherEnzymes => {
                        let cuts = digest_choices(source, pick);
                        if let Some((mut a, mut b)) = selected_digest_pair(source, pick, vector) {
                            for (label, value) in [("From", &mut a), ("To", &mut b)] {
                                egui::ComboBox::from_id_salt(("gg_other_cut", vector, idx, label))
                                    .selected_text(format!("{}: {} / {}", crate::i18n::text(ui.ctx(),label),value.top, value.bottom))
                                    .show_ui(ui, |ui| { for cut in &cuts { ui.selectable_value(value, *cut, format!("{} / {}", cut.top, cut.bottom)); } });
                            }
                            pick.cut_pair = Some((a, b));
                        } else { ui.label(crate::i18n::text(ui.ctx(),"Choose restriction enzymes to see digest fragments.")); }
                    }
                }
            }
        });
    });
    if mode == GoldenGateMode::PcrDigest && has_source {
        // Same coordinates as the map selection, exposed here for keyboard access.
        region_control(ui, dlg, vector);
    }
    orientation(ui, dlg, vector);
    let (rect, _) = ui.allocate_exact_size(Vec2::new(ui.available_width(), 105.0), egui::Sense::hover());
    if vector { draw_vector(ui.painter(), rect, false); }
    else { draw_join(ui.painter(), rect, true, direct_size.is_some()); }
    let pick = if vector { &dlg.vector } else { &dlg.fragments[idx] };
    if let Ok(fragment) = prepared(dlg, pick, vector) {
        let up = fragment.ends.upstream_stickiness.max(0) as usize;
        let down = fragment.ends.downstream_stickiness.max(0) as usize;
        let left = &fragment.seq[..up.min(fragment.seq.len())];
        let right = &fragment.seq[fragment.seq.len().saturating_sub(down)..];
        for (at, text) in [(Pos2::new(rect.left() + rect.width() * 0.28, rect.center().y - 15.0), left),
            (Pos2::new(rect.left() + rect.width() * 0.72, rect.center().y - 15.0), right)] {
            ui.painter().text(at, egui::Align2::CENTER_CENTER, String::from_utf8_lossy(text), egui::FontId::monospace(12.0), theme::TEXT);
        }
    }
}

fn draw_vector(painter: &egui::Painter, rect: Rect, ready: bool) {
    let body = Rect::from_center_size(rect.center(), Vec2::new(rect.width() * 0.88, rect.height() * 0.65));
    capsule_open(painter, body, 0.42, Stroke::new(16.0, Color32::from_gray(if ready { 90 } else { 220 })));
}

/// Explicit panes: unlike nested panels, the PCR editors and action buttons may
/// not be clipped by the source table or a second CentralPanel.
pub(super) fn draw_dialog(ui: &mut egui::Ui, ctx: &egui::Context, dlg: &mut CloningDlg) -> Option<Outcome> {
    let all = ui.max_rect();
    ui.painter().rect_filled(all, 0.0, Color32::WHITE);
    let gap = 9.0;
    let menu_h = 25.0;
    let header_h = 102.0;
    let tab_h = 34.0;
    let side_w = (all.width() * 0.275).clamp(350.0, 480.0);
    let bottom_h = (all.height() * 0.235).clamp(215.0, 285.0);
    let top = all.top() + header_h + tab_h;
    let split_x = all.right() - side_w;
    let split_y = all.bottom() - bottom_h;
    let left_top = Rect::from_min_max(Pos2::new(all.left() + gap, top), Pos2::new(split_x - gap, split_y - gap));
    let right_top = Rect::from_min_max(Pos2::new(split_x, top), Pos2::new(all.right() - gap, split_y - gap));
    let left_bottom = Rect::from_min_max(Pos2::new(all.left() + gap, split_y), Pos2::new(split_x - gap, all.bottom() - gap));
    let right_bottom = Rect::from_min_max(Pos2::new(split_x, split_y), Pos2::new(all.right() - gap, all.bottom() - gap));
    let mut menu = ui.new_child(egui::UiBuilder::new().id_salt("gg_menu").max_rect(
        Rect::from_min_size(all.left_top(), Vec2::new(all.width(), menu_h))));
    egui::MenuBar::new().ui(&mut menu, |ui| {
        ui.menu_button(crate::i18n::text(ui.ctx(),"File"), |ui| {
            if ui.button(crate::i18n::text(ui.ctx(),"Close this cloning window")).clicked() { ctx.send_viewport_cmd(egui::ViewportCommand::Close); }
        });
        for title in ["Edit", "View", "Enzymes", "Features", "Primers", "Actions", "Tools", "Window", "Help"] {
            ui.menu_button(crate::i18n::text(ui.ctx(),title), |ui| {
                if title == "View" {
                    if let Some(app) = dlg.preview_app.as_mut() {
                        for (label, tab) in [("Map", crate::app::Tab::Map), ("Sequence", crate::app::Tab::Sequence),
                            ("Enzymes", crate::app::Tab::Enzymes), ("Primers", crate::app::Tab::Primers)] {
                            if ui.button(crate::i18n::text(ui.ctx(),label)).clicked() { app.tab = tab; ui.close(); }
                        }
                    }
                } else if title == "Enzymes" {
                    ui.label(crate::i18n::text(ui.ctx(),"Golden Gate Enzymes"));
                    let names = dlg.gg_enzymes.clone();
                    for (i, name) in names.iter().enumerate() { ui.selectable_value(&mut dlg.gg_enzyme, i, *name); }
                } else if title == "Primers" {
                    if ui.button(crate::i18n::text(ui.ctx(),"Choose PCR Primers…")).clicked() {
                        if let Err(error) = design(dlg) { dlg.hint = error; }
                        dlg.preview_key.clear(); ui.close();
                    }
                } else if title == "Actions" {
                    if ui.button(crate::i18n::text(ui.ctx(),"Adjust Overhangs…")).clicked() { dlg.gg_adjust_open = true; ui.close(); }
                } else if title == "Help" {
                    ui.label(crate::i18n::text(ui.ctx(),"Golden Gate uses Type IIS digestion and compatible cohesive ends."));
                    ui.label(crate::i18n::text(ui.ctx(),"Edit the source sequence in the main sequence window."));
                } else {
                    ui.add_enabled(false, egui::Label::new(crate::i18n::text(ui.ctx(),"Source editing is available in the main sequence window.")));
                }
            });
        }
    });
    let mut head = ui.new_child(egui::UiBuilder::new().id_salt("gg_header").max_rect(
        Rect::from_min_size(all.left_top() + Vec2::new(gap, menu_h + 8.0), Vec2::new(all.width() - 2.0 * gap, header_h - menu_h - 8.0))));
    header(&mut head, dlg);
    let mut tabs = ui.new_child(egui::UiBuilder::new().id_salt("gg_tabs").max_rect(
        Rect::from_min_size(Pos2::new(all.left(), all.top() + header_h), Vec2::new(all.width(), tab_h))));
    tabs.horizontal(|ui| {
        ui.add_space((all.width() - 270.0) * 0.5);
        for (i, label) in dlg.tabs().iter().enumerate() {
            if method_tab(ui, dlg.tab == i, label, 86.0) { dlg.tab = i; }
        }
    });
    let painter = ui.painter().clone();
    for rect in [left_top, right_top, left_bottom] {
        painter.rect_filled(rect, 0.0, Color32::from_gray(0xf2));
        painter.rect_stroke(rect, 0.0, Stroke::new(1.0, theme::BORDER), egui::StrokeKind::Inside);
    }
    let mut source_ui = ui.new_child(egui::UiBuilder::new().id_salt("gg_source_pane").max_rect(right_top.shrink(7.0)));
    egui::ScrollArea::vertical().auto_shrink([false, false]).show(&mut source_ui, |ui| {
        ui.set_min_width(right_top.width() - 30.0);
        if dlg.tab == 2 { product(ui, dlg); }
        else { panel(ui, dlg, dlg.tab == 0); }
    });
    let ready = dlg.preview.is_some() && dlg.primers_ready();
    let actions_h = 106.0;
    let hint = Rect::from_min_max(right_bottom.left_top(), Pos2::new(right_bottom.right(), right_bottom.bottom() - actions_h));
    painter.rect_filled(hint, 0.0, if ready { Color32::from_rgb(0xd1, 0xfa, 0xd1) } else { Color32::from_rgb(0xff, 0xff, 0x85) });
    painter.rect_stroke(hint, 0.0, Stroke::new(1.0, theme::BORDER), egui::StrokeKind::Inside);
    let mut hint_ui = ui.new_child(egui::UiBuilder::new().id_salt("gg_hint").max_rect(hint.shrink(10.0)));
    hint_ui.centered_and_justified(|ui| { ui.label(crate::i18n::message(ui.ctx(),&dlg.hint_text())); });
    let mut result = None;
    let mut actions = ui.new_child(egui::UiBuilder::new().id_salt("gg_actions").max_rect(
        Rect::from_min_max(Pos2::new(right_bottom.left(), hint.bottom() + 8.0), right_bottom.right_bottom())));
    actions.set_clip_rect(right_bottom);
    actions.horizontal(|ui| {
        ui.label(crate::i18n::text(ui.ctx(),"Create product:"));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| { ui.checkbox(&mut dlg.close_after_clone, crate::i18n::text(ui.ctx(),"and close this window")); });
    });
    actions.add(egui::TextEdit::singleline(&mut dlg.product_name).desired_width(actions.available_width()));
    actions.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
        if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() { result = Some(Outcome::Cancel); }
        if ui.add_enabled(ready, egui::Button::new(crate::i18n::text(ui.ctx(),"Assemble"))).clicked() {
            if let Some(p) = &dlg.preview { result = Some(Outcome::Clone(super::output::product_bytes(dlg,p, &dlg.product_name), dlg.product_name.clone(), dlg.close_after_clone)); }
        }
    });
    let mut bottom = ui.new_child(egui::UiBuilder::new().id_salt("gg_bottom").max_rect(left_bottom.shrink(9.0)));
    bottom.set_clip_rect(left_bottom);
    if dlg.tab == 2 {
        let rect = bottom.max_rect();
        draw_join(bottom.painter(), rect, false, ready);
        let x = rect.center().x;
        let vector_y = rect.center().y + rect.height() * 0.35;
        let fragment_y = rect.center().y - rect.height() * 0.08;
        for (y, vector) in [(vector_y, true), (fragment_y, false)] {
            for (dx, reversed, label) in [(-65.0, false, "→"), (35.0, true, "←")] {
                if bottom.put(Rect::from_min_size(Pos2::new(x + dx, y), Vec2::new(30.0, 21.0)), egui::Button::new(label)).clicked() {
                    if vector { dlg.vector.reversed = reversed; } else { dlg.fragments[dlg.frag_index].reversed = reversed; }
                }
            }
        }
    }
    else { primer_panel(&mut bottom, dlg); }
    let mut canvas = ui.new_child(egui::UiBuilder::new().id_salt("gg_canvas").max_rect(left_top.shrink(5.0)));
    let mut selected_region = None;
    dlg.ensure_preview(ctx);
    if dlg.tab == 2 && dlg.preview_app.is_none() {
        // Keep the same preview tabs/toolbar in the not-ready product state as
        // SG. This empty view is never offered as a generated product.
        let bytes = hx_formats::write::NewFile::dna(Vec::new()).to_bytes();
        if let Ok(doc) = crate::doc::Document::from_bytes(bytes) {
            dlg.preview_app = Some(Box::new(crate::HelixApp::preview(ctx, doc)));
        }
    }
    if let Some(app) = dlg.preview_app.as_mut() {
        let chosen: Vec<_> = hx_enzymes::ENZYMES.iter().enumerate().filter(|(_, e)| dlg.gg_enzymes.contains(&e.name)).map(|(i, _)| i).collect();
        if app.sites.selected != chosen {
            app.sites.selected = chosen;
            if let Some(doc) = app.doc.as_ref() { app.sites.recompute(doc.sequence(), doc.circular()); }
            app.sites_dirty = false;
        }
        app.doc_canvas(&mut canvas, crate::app::CanvasMode::Cloning);
        if dlg.tab < 2 {
            if let Some((a, b)) = app.selection {
                selected_region = Some((a, b, app.doc.as_ref().is_some_and(|d| d.circular())));
            }
        }
    } else if dlg.tab == 1 {
        canvas.centered_and_justified(|ui| { paint_empty_canvas_hint(ui, dlg); });
    }
    if let Some((a, b, circular)) = selected_region {
        let pick = if dlg.tab == 0 { &mut dlg.vector } else { &mut dlg.fragments[dlg.frag_index] };
        // A vector selection identifies the region to replace only on a circle;
        // a linear source has no complementary arc and must not be wrapped.
        pick.region = Some(if dlg.tab == 0 && circular { (b, a) } else { (a, b) });
    }
    adjust_window(ctx, dlg);
    result
}

fn product(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    if ui.add_sized([ui.available_width(), 25.0], egui::Button::new(crate::i18n::text(ui.ctx(),"Choose PCR Primers…"))).clicked() {
        if let Err(error) = design(dlg) {
            dlg.hint = error;
        }
        dlg.preview_key.clear();
    }
    let has_regions = dlg.vector.region.is_some() && dlg.fragments.iter().all(|p| p.region.is_some() || p.gg_mode != GoldenGateMode::PcrDigest);
    if ui.add_enabled_ui(has_regions, |ui| ui.add_sized([ui.available_width(), 25.0], egui::Button::new(crate::i18n::text(ui.ctx(),"Adjust Overhangs…")))).inner.clicked() {
        dlg.gg_adjust_open = true;
    }
    ui.add_space(25.0);
    ui.label(crate::i18n::text(ui.ctx(),"Bacterial Transformation Strain"));
    egui::ComboBox::from_id_salt("gg_strain").selected_text(STRAINS[dlg.strain.min(STRAINS.len() - 1)])
        .width(ui.available_width()).show_ui(ui, |ui| {
            for (i, name) in STRAINS.iter().enumerate() { ui.selectable_value(&mut dlg.strain, i, *name); }
        });
    ui.add_space(22.0);
    let width = ui.available_width();
    egui::Frame::NONE.fill(Color32::from_gray(0xef)).stroke(Stroke::new(1.0, theme::BORDER)).inner_margin(10).show(ui, |ui| {
        ui.set_width(width - 22.0);
        let mut picks = vec![("Linearized Vector".to_string(), &dlg.vector, true)];
        picks.extend(dlg.fragments.iter().enumerate().map(|(i, p)|
            (if dlg.method.is_multiple() { format!("Fragment {}", i + 1) } else { "Fragment".into() }, p, false)));
        for (title, pick, vector) in picks {
            ui.horizontal(|ui| {
                ui.label(&title);
                if let Ok(source) = prepared(dlg, pick, vector) { ui.label(crate::ui_format!(ui.ctx(),"{} bp", source.display_len())); }
            });
            ui.separator();
                    ui.label(crate::ui_format!(ui.ctx(),"Source: {}", pick.doc.and_then(|i| dlg.docs.get(i)).map(|d| d.name.as_str()).unwrap_or("")));
            ui.add_enabled_ui(false, |ui| {
                let mut make = false; ui.checkbox(&mut make, crate::i18n::text(ui.ctx(),"Make file:"));
                let mut name = format!("{title}.dna"); ui.add(egui::TextEdit::singleline(&mut name).desired_width(ui.available_width()));
            }).response.on_hover_text(crate::i18n::text(ui.ctx(),"Intermediate-file creation is not implemented yet; the final product can be generated."));
            let (rect, _) = ui.allocate_exact_size(Vec2::new(ui.available_width(), 62.0), egui::Sense::hover());
            let c = if vector { Color32::from_gray(220) } else { Color32::from_rgb(255, 191, 195) };
            let x = rect.left() + 5.0; let y = rect.center().y;
            for (dx, dy, col) in [(0.0, 7.0, c), (58.0, -7.0, if vector { Color32::from_rgb(255, 191, 195) } else { Color32::from_gray(220) })] {
                ui.painter().line_segment([Pos2::new(x + dx, y + dy), Pos2::new(x + dx + 40.0, y + dy)], Stroke::new(10.0, col));
            }
            ui.add_space(12.0);
        }
        ui.set_min_height(320.0);
    });
}

pub(super) fn primer_panel(ui: &mut egui::Ui, dlg: &mut CloningDlg) {
    let vector = dlg.tab == 0;
    let idx = dlg.frag_index.min(dlg.fragments.len().saturating_sub(1));
    let pick = if vector { &mut dlg.vector } else { &mut dlg.fragments[idx] };
    let enabled = pick.doc.is_some() && pick.gg_mode == GoldenGateMode::PcrDigest;
    let prefix = if vector { "Vector" } else { "Fragment" };
    for (i, (primer, phosphorylated)) in [(&mut pick.primer1, &mut pick.phos1), (&mut pick.primer2, &mut pick.phos2)].into_iter().enumerate() {
        ui.add_enabled_ui(enabled, |ui| {
            ui.horizontal(|ui| {
                ui.label(crate::ui_format!(ui.ctx(),"Primer {}:", i + 1));
                egui::ComboBox::from_id_salt(("gg_primer", vector, idx, i))
                    .selected_text(if enabled { format!("<{prefix}.{}>", if i == 0 { "FOR" } else { "REV" }) } else { String::new() })
                    .width(200.0).show_ui(ui, |ui| { ui.label(crate::i18n::text(ui.ctx(),"Enter the complete primer sequence below.")); });
            });
            ui.add_space(7.0);
            ui.horizontal(|ui| {
                ui.add_space(62.0); ui.label(RichText::new("5′").italics());
                ui.add(egui::TextEdit::singleline(primer).desired_width((ui.available_width() - 35.0).max(80.0)));
                ui.label(RichText::new("3′").italics());
            });
            ui.horizontal(|ui| {
                ui.add_space(75.0); ui.checkbox(phosphorylated, crate::i18n::text(ui.ctx(),"5′ Phosphorylated"));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button(crate::i18n::text(ui.ctx(),"Reverse Complement")).clicked() { *primer = String::from_utf8_lossy(&hx_algo::reverse_complement(primer.as_bytes())).into_owned(); }
                });
            });
        });
        ui.add_space(14.0);
    }
}

pub(super) fn design(dlg: &mut CloningDlg) -> Result<(), String> {
    let e = enzyme(dlg)?;
    let width = (e.site_len() as i32 - e.reverse_cut as i32 - e.cut as i32) as usize;
    if width == 0 || width > 8 { return Err("Unsupported overhang length for PCR design.".into()); }
    let count = if dlg.method.is_multiple() { dlg.fragments.len() } else { 1 };
    if dlg.gg_overhangs.len() != count + 1 { dlg.gg_overhangs.resize(count + 1, String::new()); }
    if dlg.gg_overhangs.iter().any(|s| s.len() != width || !s.bytes().all(|b| b"ACGT".contains(&b.to_ascii_uppercase()))) {
        dlg.gg_adjust_open = true;
        return Err(format!("Please specify {width}-base junction overhangs in Adjust Overhangs before designing PCR primers."));
    }
    let mut designs = Vec::new();
    for (i, pick) in std::iter::once(&dlg.vector).chain(dlg.fragments.iter().take(count)).enumerate() {
        if pick.gg_mode != GoldenGateMode::PcrDigest { continue; }
        let source = dlg.pick_source(pick).ok_or("Please choose both sources.")?;
        let (a, b) = pick.region.ok_or("Please select the amplified region on each PCR source.")?;
        let body = source.slice(a, b, pick.reversed);
        if body.seq.len() < 24 { return Err("The amplified region is too short for PCR primers.".into()); }
        let next = (i + 1) % (count + 1);
        let spacer = "A".repeat((e.cut as i32 - e.site_len() as i32).max(0) as usize);
        // SG `goldenGateCloning\upstreamBases`：Type IIS 位点上游的保护碱基。
        // 用 GC 交替，既避免形成新的识别位点，也给酶足够的切割空间。
        let clamp: String = (0..dlg.gg_upstream_bases())
            .map(|j| if j % 2 == 0 { 'G' } else { 'C' })
            .collect();
        let left = format!("{clamp}{}{spacer}{}", e.site, dlg.gg_overhangs[i].to_ascii_uppercase());
        let right = format!("{clamp}{}{spacer}{}", e.site,
            String::from_utf8_lossy(&hx_algo::reverse_complement(dlg.gg_overhangs[next].as_bytes())));
        designs.push((i, format!("{left}{}", String::from_utf8_lossy(&body.seq[..20])),
            format!("{right}{}", String::from_utf8_lossy(&hx_algo::reverse_complement(&body.seq[body.seq.len() - 20..])))));
    }
    for (i, forward, reverse) in designs {
        let pick = if i == 0 { &mut dlg.vector } else { &mut dlg.fragments[i - 1] };
        pick.primer1 = forward; pick.primer2 = reverse; pick.cut_pair = None;
    }
    Ok(())
}

pub(super) fn adjust_window(ctx: &egui::Context, dlg: &mut CloningDlg) {
    if !dlg.gg_adjust_open { return; }
    let count = dlg.fragments.len() + 1; dlg.gg_overhangs.resize(count, String::new());
    let mut open = true;
    egui::Window::new(crate::i18n::text(ctx,"Adjust Overhangs")).open(&mut open).resizable(false).show(ctx, |ui| {
        ui.label(crate::i18n::text(ui.ctx(),"Junction sequences are written in the top-strand direction."));
        for i in 0..count {
            ui.horizontal(|ui| { ui.label(if i == 0 { crate::i18n::text(ctx,"Vector left end").to_string() } else { crate::ui_format!(ctx,"Junction {}",i) });
                ui.add(egui::TextEdit::singleline(&mut dlg.gg_overhangs[i]).desired_width(100.0)); });
        }
    });
    dlg.gg_adjust_open = open;
}

pub(super) fn draw_join(painter: &egui::Painter, rect: Rect, fragment_only: bool, ready: bool) {
    let gray = Color32::from_gray(if ready { 90 } else { 220 });
    let pink = Color32::from_rgb(255, if ready { 110 } else { 191 }, if ready { 125 } else { 195 });
    let body = Rect::from_center_size(rect.center(), Vec2::new(rect.width() * 0.50, rect.height() * 0.55));
    if !fragment_only { capsule_open(painter, body, 0.45, Stroke::new(15.0, gray)); }
    let y = if fragment_only { rect.center().y } else { body.top() };
    let x = rect.center().x; let w = rect.width() * 0.25;
    painter.line_segment([Pos2::new(x - w * 0.5, y - 4.0), Pos2::new(x + w * 0.5, y - 4.0)], Stroke::new(9.0, pink));
    painter.line_segment([Pos2::new(x - w * 0.5 + 12.0, y + 5.0), Pos2::new(x + w * 0.5 + 12.0, y + 5.0)], Stroke::new(9.0, pink));
    painter.text(Pos2::new(x, y - 14.0), egui::Align2::CENTER_BOTTOM, crate::i18n::text(painter.ctx(),"Fragment"), egui::FontId::proportional(11.0), theme::TEXT);
    if !fragment_only { painter.text(Pos2::new(x, body.bottom() + 15.0), egui::Align2::CENTER_TOP, crate::i18n::text(painter.ctx(),"Vector"), egui::FontId::proportional(11.0), theme::TEXT); }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> CloningDlg {
        let docs = vec![DocSnap::new("vector", Source::new("vector", b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".to_vec(), true)),
            DocSnap::new("insert", Source::new("insert", b"CCCCAAAATTTTGGGGCCCCAAAATTTTGGGG".to_vec(), false))];
        let mut dlg = CloningDlg::new(Method::GoldenGate, docs, Some(0));
        dlg.fragments[0].doc = Some(1);
        dlg.vector.region = Some((0, 32));
        dlg.fragments[0].region = Some((0, 32));
        dlg.gg_overhangs = vec!["AATG".into(), "GCTT".into()];
        dlg
    }

    #[test]
    fn defaults_follow_sg_not_generic_overlap_workflow() {
        let mut dlg = fixture();
        assert_eq!(dlg.vector.gg_mode, GoldenGateMode::PcrDigest);
        assert!(!dlg.overlap_method());
        assert!(dlg.tab_shows_primers());
        dlg.fragments[0].doc = None;
        dlg.tab = 1;
        assert!(dlg.tab_shows_primers(), "disabled PCR editors remain visible when no source is selected");
        let fresh = CloningDlg::new(Method::GoldenGate, Vec::new(), None);
        assert!(fresh.fragments[0].doc.is_none());
    }

    #[test]
    fn selected_enzyme_is_part_of_the_product_cache_key() {
        let mut dlg = fixture(); let key = dlg.state_key(); dlg.gg_enzyme = 1;
        assert_ne!(key, dlg.state_key());
    }

    #[test]
    fn pcr_primers_create_real_type_iis_cuts_and_seamless_product() {
        let mut dlg = fixture();
        assert!(build(&dlg).is_err());
        design(&mut dlg).unwrap();
        let v = prepared(&dlg, &dlg.vector, true).unwrap();
        let f = prepared(&dlg, &dlg.fragments[0], false).unwrap();
        assert_eq!(v.ends.upstream_stickiness, 4); assert_eq!(f.ends.downstream_stickiness, 4);
        let product = build(&dlg).unwrap();
        assert_eq!(product.seq.len(), 72);
        assert!(product.circular);
        assert!(!String::from_utf8_lossy(&product.seq).contains("GGTCTC"));
        dlg.fragments[0].primer1 = "ACGT".into();
        assert!(build(&dlg).is_err());
    }

    #[test]
    fn direct_circular_input_and_missing_flanking_sites_are_not_assembled() {
        let mut dlg = fixture(); dlg.vector.gg_mode = GoldenGateMode::Direct;
        assert!(prepared(&dlg, &dlg.vector, true).unwrap_err().contains("linear"));
        dlg.vector.gg_mode = GoldenGateMode::Digest;
        assert!(prepared(&dlg, &dlg.vector, true).is_err());
    }

    #[test]
    fn digest_mode_uses_real_flanking_sites_and_sg_default_detection() {
        let src = Source::new("flanked", b"GGTCTCAAATGCCCCCCGCTTTGAGACC".to_vec(), false);
        let dlg = CloningDlg::new(Method::GoldenGate, vec![DocSnap::new("flanked", src)], Some(0));
        assert_eq!(dlg.vector.gg_mode, GoldenGateMode::Digest);
        assert!(!dlg.fragments[0].doc.is_some());
        let fragment = prepared(&dlg, &dlg.vector, true).unwrap();
        assert_eq!(fragment.seq, b"AATGCCCCCCGCTT");
        assert_eq!(fragment.display_len(), 10);
    }

    #[test]
    fn sap_i_uses_three_base_ends_and_global_enzyme_changes_the_result() {
        let mut dlg = fixture();
        dlg.gg_enzyme = dlg.gg_enzymes.iter().position(|n| *n == "SapI").unwrap();
        dlg.gg_overhangs = vec!["AAT".into(), "GCT".into()];
        design(&mut dlg).unwrap();
        assert_eq!(prepared(&dlg, &dlg.vector, true).unwrap().ends.upstream_stickiness, 3);
        assert_eq!(build(&dlg).unwrap().seq.len(), 70);
        dlg.gg_enzyme = 0;
        assert!(build(&dlg).is_err());
    }

    #[test]
    fn multiple_fragment_pcr_design_prepares_every_join() {
        let mut dlg = fixture(); dlg.method = Method::GoldenGateMulti;
        dlg.fragments.push(dlg.fragments[0].clone());
        dlg.gg_overhangs = vec!["AATG".into(), "GCTT".into(), "CGAC".into()];
        design(&mut dlg).unwrap();
        assert_eq!(build(&dlg).unwrap().seq.len(), 108);
        assert!(dlg.fragments.iter().all(|p| !p.primer1.is_empty() && !p.primer2.is_empty()));
    }

    #[test]
    fn every_golden_gate_tab_renders_at_large_and_minimum_sizes() {
        for size in [Vec2::new(1800.0, 1100.0), Vec2::new(1220.0, 840.0)] {
            for tab in 0..3 {
                let mut dlg = fixture(); dlg.fragments[0].doc = None; dlg.tab = tab;
                dlg.refresh();
                let ctx = egui::Context::default();
                let input = egui::RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, size)), ..Default::default() };
                let mut output = ctx.run_ui(input, |ui| {
                    let ctx = ui.ctx().clone();
                    let result = draw_dialog(ui, &ctx, &mut dlg);
                    assert!(result.is_none());
                });
                assert!(!output.shapes.is_empty());
                output.textures_delta.clear();
                assert!(dlg.preview.is_none(), "a rendered not-ready product is not an assembled molecule");
            }
        }
    }
}
