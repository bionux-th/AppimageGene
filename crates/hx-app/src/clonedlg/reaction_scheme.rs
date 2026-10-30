//! 反应图不是“按页签顺序摆图标”：每行有明确的共同反应物和一个产物。
//! 模型与绘制分开，防止 Destination 等输入被误画成上一步的产物。
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MoleculeShape { Linear, Circular, Open }

#[derive(Debug, Clone)]
pub(super) struct Molecule {
    pub label: String,
    pub shape: MoleculeShape,
    pub size: Option<usize>,
    pub available: bool,
    pub sites: String,
    pub inserts: usize,
    pub output: bool,
}

pub(super) struct Reaction {
    pub inputs: Vec<Molecule>,
    pub output: Molecule,
    pub mechanism: &'static str,
}

pub(super) struct Scheme { pub rows: Vec<Reaction> }

fn input(dlg: &CloningDlg, pick: &Pick, label: &str, fallback: MoleculeShape, sites: &str, inserts: usize) -> Molecule {
    let source = dlg.pick_source(pick);
    Molecule {
        label: label.into(),
        shape: source.map(|s| if s.circular { MoleculeShape::Circular } else { MoleculeShape::Linear }).unwrap_or(fallback),
        size: source.map(Source::display_len), available: source.is_some(),
        sites: sites.into(), inserts, output: false,
    }
}

fn output(label: &str, shape: MoleculeShape, product: Option<&Product>, sites: &str, inserts: usize) -> Molecule {
    Molecule { label: label.into(), shape, size: product.map(Product::display_len), available: product.is_some(),
        sites: sites.into(), inserts, output: true }
}

pub(super) fn model(dlg: &CloningDlg) -> Scheme {
    let ready = if dlg.primers_ready() { dlg.preview.as_ref() } else { None };
    let mut rows = Vec::new();
    match dlg.method {
        Method::GatewayBp | Method::GatewayBpLr(_) => {
            if matches!(dlg.method,Method::GatewayBpLr(n) if n>1) {
                let mut lr_inputs=Vec::new();
                for (i,pick) in dlg.fragments.iter().enumerate() {
                    let bp=dlg.build_gateway_entry(i).ok();
                    let entry=output(&format!("Entry Clone {}",i+1),MoleculeShape::Circular,bp.as_ref(),"numbered attL / attR",1);
                    rows.push(Reaction{inputs:vec![input(dlg,pick,&format!("attB Fragment {}",i+1),MoleculeShape::Linear,"numbered attB",1),
                        input(dlg,&dlg.gateway_donors[i],&format!("Donor Vector {}",i+1),MoleculeShape::Circular,"matching attP",0)],
                        output:entry.clone(),mechanism:"BP"});
                    lr_inputs.push(Molecule{output:false,..entry});
                }
                lr_inputs.push(input(dlg,&dlg.destination,"Destination Vector",MoleculeShape::Circular,"attR sites",0));
                rows.push(Reaction{inputs:lr_inputs,output:output("Expression Clone",MoleculeShape::Circular,ready,"numbered attB scars",dlg.fragments.len()),mechanism:"LR"});
                return Scheme{rows};
            }
            let bp = dlg.build_gateway_bp().ok();
            let bp = if dlg.primers_ready() { bp.as_ref() } else { None };
            let entry = output("Entry Clone", MoleculeShape::Circular, bp, "attL1 / attL2", 1);
            rows.push(Reaction {
                inputs: vec![input(dlg, &dlg.fragments[0], "attB Fragment", MoleculeShape::Linear, "attB1 / attB2", 1),
                    input(dlg, &dlg.vector, "Donor Vector", MoleculeShape::Circular, "attP1 / attP2", 0)],
                output: entry.clone(), mechanism: "BP",
            });
            if matches!(dlg.method, Method::GatewayBpLr(_)) {
                let mut inputs = vec![Molecule { output: false, ..entry }];
                for (i, pick) in dlg.fragments.iter().enumerate().skip(1) {
                    inputs.push(input(dlg, pick, &format!("Entry Clone {}", i + 1), MoleculeShape::Circular, "att sites", 1));
                }
                inputs.push(input(dlg, &dlg.destination, "Destination Vector", MoleculeShape::Circular, "attR sites", 0));
                rows.push(Reaction { inputs,
                    output: output("Expression Clone", MoleculeShape::Circular, ready, "attB sites", dlg.fragments.len()), mechanism: "LR" });
            }
        }
        Method::GatewayLr(_) => {
            let mut inputs = Vec::new();
            for (i, pick) in dlg.fragments.iter().enumerate() {
                let label = if dlg.fragments.len() == 1 { "Entry Clone".into() } else { format!("Entry Clone {}", i + 1) };
                inputs.push(input(dlg, pick, &label, MoleculeShape::Circular,
                    if dlg.fragments.len()>1 {"numbered attL / attR"}else{"attL sites"}, 1));
            }
            inputs.push(input(dlg, &dlg.vector, "Destination Vector", MoleculeShape::Circular, "attR sites", 0));
            rows.push(Reaction { inputs,
                output: output("Expression Clone", MoleculeShape::Circular, ready, "attB sites", dlg.fragments.len()), mechanism: "LR" });
        }
        _ => {
            let ligation = dlg.method.is_ligation();
            let assembly = matches!(dlg.method, Method::GibsonAssemble | Method::NeBuilderAssemble);
            let mut inputs = Vec::new();
            if !ligation {
                let mut vector = input(dlg, &dlg.vector, if assembly { "Fragment 1" } else { "Linearized Vector" }, MoleculeShape::Open, "", 0);
                vector.shape = MoleculeShape::Open;
                if let Some(source) = dlg.pick_source(&dlg.vector) {
                    if dlg.method.is_restriction() {
                        let prepared = restriction_fragment(source, &dlg.vector, true, true).ok();
                        vector.size = prepared.as_ref().map(Source::display_len);
                        vector.available = prepared.is_some();
                    } else {
                        vector.size = Some(dlg.vector_slice(source).seq.len());
                    }
                }
                inputs.push(vector);
            }
            let count = if dlg.method.is_multiple() || ligation { dlg.fragments.len() } else { 1 };
            for (i, pick) in dlg.fragments.iter().take(count).enumerate() {
                let label = if count == 1 { "Fragment".into() } else { format!("Fragment {}", i + 1 + usize::from(assembly)) };
                let mut molecule = input(dlg, pick, &label, MoleculeShape::Linear, "", 1);
                molecule.shape = MoleculeShape::Linear;
                if let Some(source) = dlg.pick_source(pick) {
                    if dlg.method.is_restriction() || ligation {
                        let prepared = restriction_fragment(source, pick, false, dlg.method.is_restriction()).ok();
                        molecule.size = prepared.as_ref().map(Source::display_len);
                        molecule.available = prepared.is_some();
                    } else { molecule.size = Some(dlg.fragment_slice(source, pick).seq.len()); }
                }
                inputs.push(molecule);
            }
            rows.push(Reaction { inputs,
                output: output("Product", if ligation { MoleculeShape::Linear } else { MoleculeShape::Circular }, ready, "", count),
                mechanism: if ligation { "Ligate" } else if dlg.method.is_restriction() { "Clone" } else { "Assemble" },
            });
        }
    }
    Scheme { rows }
}

pub(super) fn draw(ui: &egui::Ui, rect: Rect, dlg: &CloningDlg) {
    let scheme = model(dlg);
    let painter = ui.painter_at(rect);
    let title = crate::i18n::text(ui.ctx(),"Reaction scheme · not to scale");
    painter.text(rect.left_top() + Vec2::new(8.0, 4.0), egui::Align2::LEFT_TOP, title,
        egui::FontId::proportional(10.0), theme::TEXT_DIM);
    let area = Rect::from_min_max(rect.left_top() + Vec2::new(6.0, 20.0), rect.right_bottom() - Vec2::new(6.0, 3.0));
    let height = area.height() / scheme.rows.len().max(1) as f32;
    for (index, row) in scheme.rows.iter().enumerate() {
        let row_rect = Rect::from_min_size(area.left_top() + Vec2::new(0.0, height * index as f32), Vec2::new(area.width(), height));
        let slots = row.inputs.len() + 1;
        let slot = row_rect.width() / slots as f32;
        let y = row_rect.top() + row_rect.height() * 0.54;
        let scale = (slot / 140.0).min(1.0);
        for (i, molecule) in row.inputs.iter().chain(std::iter::once(&row.output)).enumerate() {
            let x = row_rect.left() + slot * (i as f32 + 0.5);
            draw_molecule(&painter, molecule, Pos2::new(x, y), slot, scale);
            if i > 0 {
                let at = Pos2::new(row_rect.left() + slot * i as f32, y);
                let is_output = i == row.inputs.len();
                painter.text(at, egui::Align2::CENTER_CENTER, if is_output { "→" } else { "+" },
                    egui::FontId::proportional(18.0), theme::TEXT);
                if is_output {
                    painter.text(at - Vec2::new(0.0, 19.0), egui::Align2::CENTER_BOTTOM, row.mechanism,
                        egui::FontId::proportional(9.0), theme::TEXT_DIM);
                }
            }
        }
        if index + 1 < scheme.rows.len() {
            painter.line_segment([Pos2::new(area.left(), row_rect.bottom()), Pos2::new(area.right(), row_rect.bottom())], Stroke::new(0.5, theme::BORDER));
        }
    }
}

fn draw_molecule(painter: &egui::Painter, molecule: &Molecule, center: Pos2, slot: f32, scale: f32) {
    let ink = if molecule.available { theme::BACKBONE } else { Color32::from_gray(0xb5) };
    let color = if molecule.available { Color32::from_rgb(0xf4, 0x67, 0x72) } else { Color32::from_gray(0xcc) };
    let width = (slot * 0.65).min(135.0);
    let half_height = 14.0 * scale.max(0.6);
    let body = Rect::from_center_size(center, Vec2::new(width, half_height * 2.0));
    painter.text(Pos2::new(center.x, center.y - half_height - 6.0), egui::Align2::CENTER_BOTTOM,
        crate::i18n::text(painter.ctx(),&molecule.label), egui::FontId::proportional((10.5 * scale).max(8.0)), theme::TEXT);
    match molecule.shape {
        MoleculeShape::Linear => {
            for dy in [-1.5, 1.5] {
                painter.line_segment([Pos2::new(body.left(), center.y + dy), Pos2::new(body.right(), center.y + dy)], Stroke::new(1.0, ink));
            }
        }
        MoleculeShape::Circular => {
            // 一条闭合骨架及其中的插入段；不再画成“两个互不相关的同心环”。
            painter.rect_stroke(body, half_height as u8, Stroke::new(2.0, ink), egui::StrokeKind::Inside);
        }
        MoleculeShape::Open => capsule_open(painter, body, 0.40, Stroke::new(2.0, ink)),
    }
    if molecule.inserts > 0 {
        let y = if molecule.shape == MoleculeShape::Circular { body.top() + 1.0 } else { center.y };
        let span = width * 0.45;
        for i in 0..molecule.inserts {
            let w = span / molecule.inserts as f32;
            let tint = if i % 2 == 0 { color } else if molecule.available { Color32::from_rgb(0x7d, 0xb2, 0xe8) } else { color };
            painter.line_segment([Pos2::new(center.x - span * 0.5 + w * i as f32, y),
                Pos2::new(center.x - span * 0.5 + w * (i + 1) as f32 - 1.0, y)], Stroke::new(4.0, tint));
        }
    }
    let size = molecule.size.map(|n| crate::ui_format!(painter.ctx(),"{} bp",n)).unwrap_or_else(|| {
        crate::i18n::text(painter.ctx(),if molecule.output { "Not generated" } else { "Not selected" }).into()
    });
    painter.text(Pos2::new(center.x, center.y + half_height + 5.0), egui::Align2::CENTER_TOP, size,
        egui::FontId::proportional(9.0), theme::TEXT_DIM);
    if !molecule.sites.is_empty() {
        painter.text(Pos2::new(center.x, center.y + half_height + 18.0), egui::Align2::CENTER_TOP, &molecule.sites,
            egui::FontId::proportional(8.5), theme::TEXT_DIM);
    }
}
