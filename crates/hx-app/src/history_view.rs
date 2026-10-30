//! Molecular lineage canvas: current product above its source molecules.
use eframe::egui;
use hx_formats::history::{HistoryEntry, operation_label};
use crate::{app::HelixApp, theme};

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn history_thumbnail_keeps_features_after_the_twelfth_entry() {
        let features = (0..17).map(|i| format!(r##"<Feature name="feature{i}" directionality="1"><Segment range="{}-{}" color="#993366"/></Feature>"##, i * 100 + 1, i * 100 + 80)).collect::<String>();
        let nodes = hx_formats::history::parse_entries(&format!("<HistoryTree><Node seqLen=\"2000\"><Features>{features}</Features></Node></HistoryTree>")).unwrap();
        assert_eq!(thumbnail_features(&nodes[0]).len(), 17);
        assert!(thumbnail_features(&nodes[0]).iter().any(|f| f.name == "feature16"));
    }

    #[test]
    fn history_nodes_keep_their_own_annotations() {
        let xml = r##"<HistoryTree><Node seqLen="90"><Features><Feature name="current"><Segment range="1-9" color="#ffffff"/></Feature></Features><Node seqLen="100"><Features><Feature name="ancestor"><Segment range="20-30" color="#993366"/></Feature></Features><Primers><Primer name="old primer" sequence="ACGT"/></Primers></Node></Node></HistoryTree>"##;
        let n = hx_formats::history::parse_entries(xml).unwrap();
        assert_eq!(n[0].features[0].name, "current");
        assert_eq!(n[1].features[0].name, "ancestor");
        assert!(n[0].primers.is_empty());
        assert_eq!(n[1].primers[0].name, "old primer");
    }
    #[test]
    fn captions_use_one_based_ranges_and_insertion_boundaries() {
        for (op, kind, caption, detail) in [
            ("delete", "DNA", "DELETE\nBASES", "Remove\n526 .. 532"),
            ("replace", "DNA", "INSERT\nBASES", "Replace\n526 .. 532"),
            ("insert", "Protein", "INSERT\nAMINO\nACIDS", "Insert at\n525"),
        ] {
            let xml = format!(r#"<HistoryTree><Node type="{kind}" operation="{op}"><InputSummary val1="525" val2="532"/></Node></HistoryTree>"#);
            let nodes = hx_formats::history::parse_entries(&xml).unwrap();
            assert_eq!(edit_caption(&nodes[0]), (caption.into(), detail.into()));
        }
    }
}

fn thumbnail_features(node: &HistoryEntry) -> Vec<&hx_formats::annot::FeatureAnnotation> {
    node.features.iter().collect()
}

fn edit_caption(node: &HistoryEntry) -> (String, String) {
    let protein = node.attr("type") == "Protein";
    let noun = if protein { "AMINO\nACIDS" } else { "BASES" };
    let summary = node.summaries.first();
    let get = |key: &str| summary.and_then(|s| s.attrs.iter().find(|(k, _)| k == key)).map(|(_, v)| v.as_str()).unwrap_or("");
    let a = get("val1").parse::<usize>().unwrap_or(0);
    let b = get("val2").parse::<usize>().unwrap_or(a);
    match node.attr("operation") {
        "delete" => (format!("DELETE\n{noun}"), format!("Remove\n{} .. {b}", a + 1)),
        "replace" => (format!("INSERT\n{noun}"), format!("Replace\n{} .. {b}", a + 1)),
        "insert" | "insertBases" | "insertAminoAcids" => (format!("INSERT\n{noun}"), format!("Insert at\n{a}")),
        "insertFragment" if get("manipulation") == "insert" => ("INSERT\nFRAGMENT".into(), format!("Insert\n{a} .. {b}")),
        op => (operation_label(op).to_uppercase(), String::new()),
    }
}

fn draw_molecule(painter: &egui::Painter, p: egui::Pos2, node: &HistoryEntry, current: bool, change: Option<&HistoryEntry>) {
    let len = node.attr("seqLen").parse::<usize>().unwrap_or(1).max(1) as f32;
    let circular = node.attr("circular") == "1";
    let point = |at: f32, radius: f32| {
        let angle = at / len * std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
        p + egui::vec2(angle.cos(), angle.sin()) * radius
    };
    let stroke = egui::Stroke::new(if node.attr("strandedness") == "single" { 1.4 } else { 4.0 }, theme::TEXT);
    if circular { painter.circle_stroke(p, 88.0, stroke); }
    else { painter.line_segment([p - egui::vec2(94.0, 0.0), p + egui::vec2(94.0, 0.0)], stroke); }
    if current {
        let step = (len / 4.0).max(1.0);
        for i in 1..=4 {
            let at = step * i as f32;
            let (a, b, text) = if circular { (point(at, 83.0), point(at, 88.0), point(at, 73.0)) }
                else { let x = p.x - 94.0 + 188.0 * at / len; (egui::pos2(x, p.y), egui::pos2(x, p.y + 5.0), egui::pos2(x, p.y + 12.0)) };
            painter.line_segment([a, b], egui::Stroke::new(1.0, theme::TEXT));
            painter.text(text, egui::Align2::CENTER_CENTER, format!("{}", at as usize), egui::FontId::proportional(8.0), theme::TEXT);
        }
    }
    for (i, feature) in thumbnail_features(node).into_iter().enumerate() {
        let color = feature.segments.first().and_then(|s| egui::Color32::from_hex(&s.color).ok()).unwrap_or(egui::Color32::from_rgb(160, 200, 160));
        for seg in &feature.segments {
            if circular {
                let start = seg.range.start as f32;
                let span = if seg.range.end >= seg.range.start { (seg.range.end - seg.range.start) as f32 } else { len - start + seg.range.end as f32 };
                // Shared coordinates and duplicate annotations must not push
                // unrelated features onto arbitrary alternating tracks.
                let radius = 64.0;
                let width = 11.0;
                let arc_px = span / len * std::f32::consts::TAU * radius;
                let head = (9.0_f32).min(arc_px * 0.45) / (std::f32::consts::TAU * radius) * len;
                let forward = feature.has_forward_arrow();
                let reverse = feature.has_reverse_arrow();
                let lo = start + if reverse { head } else { 0.0 };
                let hi = start + span - if forward { head } else { 0.0 };
                // Tessellate annular bands as convex quads; a concave polygon
                // passed to convex_polygon would incorrectly fill the centre.
                let steps = ((arc_px / 3.0).ceil() as usize).clamp(1, 180);
                for j in 0..steps {
                    let a = lo + (hi - lo) * j as f32 / steps as f32;
                    let b = lo + (hi - lo) * (j + 1) as f32 / steps as f32;
                    let mut mesh = egui::Mesh::default();
                    for pos in [point(a, radius - width * 0.5), point(a, radius + width * 0.5), point(b, radius + width * 0.5), point(b, radius - width * 0.5)] { mesh.colored_vertex(pos, color); }
                    mesh.add_triangle(0, 1, 2); mesh.add_triangle(0, 2, 3);
                    painter.add(egui::Shape::mesh(mesh));
                }
                for (enabled, at, base) in [(forward, start + span, hi), (reverse, start, lo)] {
                    if enabled {
                        let mut mesh = egui::Mesh::default();
                        for pos in [point(at, radius), point(base, radius + width * 0.85), point(base, radius - width * 0.85)] { mesh.colored_vertex(pos, color); }
                        mesh.add_triangle(0, 1, 2);
                        painter.add(egui::Shape::mesh(mesh));
                    }
                }
                // SG thumbnails suppress labels which do not fit; the former
                // unconditional labels piled every tiny cloning tag at 12 o'clock.
                let label = &feature.name;
                let galley = painter.layout_no_wrap(label.clone(), egui::FontId::proportional(9.0), theme::TEXT);
                let duplicate = node.features[..i].iter().any(|f| f.name == feature.name && f.segments == feature.segments);
                if !duplicate && label != "misc_feature" && (arc_px >= galley.size().x + 5.0 || span / len > 0.018 || label == "RBS") {
                    let mid = start + span * 0.5;
                    let mut angle = mid / len * std::f32::consts::TAU;
                    if angle.cos() < 0.0 { angle += std::f32::consts::PI; }
                    let anchor = point(mid, radius - if label == "RBS" { 25.0 } else { 15.0 });
                    painter.add(egui::epaint::TextShape::new(anchor - galley.size() * 0.5, galley, theme::TEXT).with_angle_and_anchor(angle, egui::Align2::CENTER_CENTER));
                }
            } else {
                let y = p.y + 18.0 + (i % 2) as f32 * 9.0;
                let rect = egui::Rect::from_min_max(egui::pos2(p.x - 94.0 + 188.0 * seg.range.start as f32 / len, y), egui::pos2(p.x - 94.0 + 188.0 * seg.range.end as f32 / len, y + 7.0));
                painter.rect_filled(rect, 0.0, color);
            }
        }
    }
    let source = change.or_else(|| current.then_some(node));
    if circular && !current {
        let mut labels = Vec::new();
        for primer in &node.primers {
            for site in &primer.sites {
                let at = (site.range.start + site.range.end) as f32 * 0.5;
                labels.push((at, primer.name.as_str()));
            }
        }
        labels.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut previous_y = [f32::NEG_INFINITY; 2];
        for (at, name) in labels {
            let root = point(at, 90.0);
            let right = root.x >= p.x;
            let side = usize::from(right);
            let y = (root.y - 20.0).max(previous_y[side] + 11.0);
            previous_y[side] = y;
            let end = egui::pos2(p.x + if right { 116.0 } else { -116.0 }, y);
            let color = egui::Color32::from_rgb(185, 65, 245);
            painter.line_segment([root, end], egui::Stroke::new(0.6, color));
            painter.text(end, if right { egui::Align2::LEFT_CENTER } else { egui::Align2::RIGHT_CENTER }, name, egui::FontId::proportional(8.0), color);
        }
    }
    if let Some(source) = source {
        if let Some(summary) = source.summaries.first() {
            for (k, v) in &summary.attrs {
                if k != "val1" && k != "val2" { continue; }
                let Ok(at) = v.parse::<f32>() else { continue; };
                let color = if current { egui::Color32::from_rgb(100, 195, 245) } else { theme::TEXT };
                let (a, b, label) = if circular { (point(at, 81.0), point(at, 93.0), point(at, 73.0)) }
                    else { let x = p.x - 94.0 + 188.0 * at / len; (egui::pos2(x, p.y - 7.0), egui::pos2(x, p.y + 7.0), egui::pos2(x, p.y + 10.0)) };
                painter.line_segment([a, b], egui::Stroke::new(1.2, color));
                if !current {
                    let insertion = matches!(source.attr("operation"), "insert" | "insertBases" | "insertAminoAcids");
                    if insertion && k == "val2" { continue; }
                    let shown = if !insertion && k == "val1" { at as usize + 1 } else { at as usize };
                    painter.text(label, if k == "val1" && !insertion { egui::Align2::RIGHT_TOP } else { egui::Align2::LEFT_TOP }, format!("{shown}"), egui::FontId::proportional(8.0), theme::TEXT);
                }
            }
        }
    }
}

pub fn show(ui: &mut egui::Ui, app: &mut HelixApp) {
    let Some(doc) = app.doc.as_ref() else { return; };
    let name = doc.path.as_ref().and_then(|p| p.file_name()).map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| doc.map_name());
    let result = doc.session.save_bytes_with_history(&name, &name, None)
        .and_then(|(bytes, _)| hx_formats::DnaFile::parse(&bytes).map_err(|e| e.to_string()))
        .and_then(|file| file.block(7).map(|b| hx_formats::history::decode_payload(&b.payload)
            .and_then(|xml| hx_formats::history::parse_entries(&xml))).unwrap_or_else(|| Ok(Vec::new())));
    let nodes = match result {
        Ok(nodes) => nodes,
        Err(e) => { ui.colored_label(theme::TEXT_DIM, crate::ui_format!(ui.ctx(),"History could not be read: {}",e)); return; }
    };
    if doc.session.is_dirty() { ui.label(egui::RichText::new(crate::i18n::text(ui.ctx(),"Unsaved changes")).small().color(theme::TEXT_DIM)); }
    if nodes.is_empty() { ui.label(crate::i18n::text(ui.ctx(),"No recorded molecular history in this file.")); return; }
    egui::ScrollArea::both().id_salt("molecular_history_canvas").show(ui, |ui| {
        let mut positions = vec![egui::Pos2::ZERO; nodes.len()];
        let circular = nodes.iter().any(|n| n.attr("circular") == "1");
        let row_height = if circular { 340.0 } else { 165.0 };
        let top = if circular { 100.0 } else { 25.0 };
        let mut leaf = 0usize;
        let mut max_depth = 0usize;
        fn layout(i: usize, depth: usize, nodes: &[HistoryEntry], positions: &mut [egui::Pos2], leaf: &mut usize, max_depth: &mut usize, top: f32, row_height: f32) {
            *max_depth = (*max_depth).max(depth);
            let children: Vec<_> = nodes.iter().enumerate().filter(|(_, n)| n.parent == Some(i)).map(|(j, _)| j).collect();
            if children.is_empty() {
                positions[i] = egui::pos2(110.0 + *leaf as f32 * 270.0, top + depth as f32 * row_height);
                *leaf += 1;
            } else {
                for &j in &children { layout(j, depth + 1, nodes, positions, leaf, max_depth, top, row_height); }
                let x = children.iter().map(|&j| positions[j].x).sum::<f32>() / children.len() as f32;
                positions[i] = egui::pos2(x, top + depth as f32 * row_height);
            }
        }
        layout(0, 0, &nodes, &mut positions, &mut leaf, &mut max_depth, top, row_height);
        let (rect, _) = ui.allocate_exact_size(egui::vec2((leaf as f32 * 270.0 + 20.0).max(ui.available_width()), top + 150.0 + max_depth as f32 * row_height), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        for (i, node) in nodes.iter().enumerate() {
            if let Some(parent) = node.parent {
                let a = positions[i] + rect.min.to_vec2() - egui::vec2(0.0, if circular { 94.0 } else { 22.0 });
                let b = positions[parent] + rect.min.to_vec2() + egui::vec2(0.0, if circular { 145.0 } else { 66.0 });
                let bend = egui::pos2(a.x, a.y);
                painter.line_segment([a, bend], egui::Stroke::new(1.5, theme::TEXT_DIM));
                let junction = egui::pos2(b.x, bend.y);
                painter.line_segment([bend, junction], egui::Stroke::new(1.5, theme::TEXT_DIM));
                if nodes.iter().position(|n| n.parent == Some(parent)) == Some(i) {
                    painter.arrow(junction, b - junction, egui::Stroke::new(1.5, theme::TEXT_DIM));
                    let (caption, detail) = edit_caption(&nodes[parent]);
                    let mid = (junction.y + b.y) * 0.5;
                    painter.text(egui::pos2(b.x - 7.0, mid), egui::Align2::RIGHT_CENTER, crate::i18n::history_caption(ui.ctx(),&caption), egui::FontId::proportional(12.0), theme::TEXT);
                    painter.text(egui::pos2(b.x + 7.0, mid), egui::Align2::LEFT_CENTER, detail, egui::FontId::proportional(11.0), theme::TEXT);
                }
            }
        }
        for (i, node) in nodes.iter().enumerate() {
            let p = positions[i] + rect.min.to_vec2();
            let color = if i == 0 { theme::TEXT } else { egui::Color32::from_rgb(0, 104, 240) };
            draw_molecule(&painter, p, node, i == 0, node.parent.map(|j| &nodes[j]));
            let label_y = if circular { 116.0 } else { 38.0 };
            painter.text(p + egui::vec2(0.0, label_y), egui::Align2::CENTER_CENTER, node.label(), egui::FontId::proportional(12.0), color);
            let unit = match node.attr("type") { "Protein" => "aa", _ if node.attr("strandedness") == "single" => "bases", _ => "bp" };
            painter.text(p + egui::vec2(0.0, label_y + 14.0), egui::Align2::CENTER_CENTER, format!("{} {}", node.attr("seqLen"),crate::i18n::text(ui.ctx(),unit)), egui::FontId::proportional(11.0), color);
            ui.interact(egui::Rect::from_center_size(p + egui::vec2(0.0, 18.0), egui::vec2(280.0, 90.0)), ui.scope_id().with(("history_node", i)), egui::Sense::hover()).on_hover_ui(|ui| {
                ui.strong(node.attr("name"));
                for (k, v) in &node.attrs { ui.label(format!("{k}: {v}")); }
                for summary in &node.summaries { for (k, v) in &summary.attrs { ui.label(format!("{k}: {v}")); } }
                if i > 0 { ui.label(crate::i18n::text(ui.ctx(),"Historical-state restoration is not yet supported.")); }
            });
        }
    });
}
