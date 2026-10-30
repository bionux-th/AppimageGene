//! Docked native Notes editor. Only fields with a known write-back contract are
//! editable; unsupported references/attachments are visibly read-only.
use eframe::egui::{self,Color32,RichText};
use crate::{doc::Document,theme};

fn value(xml:&[u8],tag:&str)->String {hx_formats::notes::field(xml,tag).unwrap_or_default()}

fn text_field(ui:&mut egui::Ui,doc_id:&str,tag:&str,title:&str,current:&str,multiline:bool,edits:&mut Vec<(String,String)>) {
    let id=egui::Id::new(("notes-draft",doc_id,tag));
    let mut draft=current.to_string();
    ui.label(RichText::new(crate::i18n::text(ui.ctx(),title)).size(12.5).color(theme::TEXT_DIM));
    let edit=if multiline {egui::TextEdit::multiline(&mut draft).desired_rows(4)}else{egui::TextEdit::singleline(&mut draft)};
    let response=ui.add_sized([ui.available_width(),if multiline {92.0}else{30.0}],edit.id(id).font(egui::FontId::proportional(13.0)));
    // Apply during editing, not on focus loss: Save/Ctrl+S must include the
    // current field, and Undo must not be overwritten by a stale UI draft.
    if response.changed() && draft!=current {edits.push((tag.to_string(),draft));}
}

fn section(ui:&mut egui::Ui,contents:impl FnOnce(&mut egui::Ui)) {
    egui::Frame::new().fill(Color32::WHITE).inner_margin(egui::Margin::same(16)).show(ui,|ui|{
        ui.spacing_mut().item_spacing.y=8.0;ui.set_width(ui.available_width());contents(ui);
    });ui.add_space(8.0);
}

pub fn draw(ui:&mut egui::Ui,doc:&mut Document)->Result<(),String> {
    let notes=doc.session.notes_xml();
    let id=value(&notes,"UUID");let identity=if id.is_empty(){doc.file_label()}else{id};
    let mut edits=Vec::new();
    egui::ScrollArea::vertical().id_salt("native_info_scroll").auto_shrink([false,false]).show(ui,|ui|{
        section(ui,|ui|{
            ui.columns(2,|cols|{
                for (col,tag,title) in [(0,"Created","Created"),(1,"LastModified","Modified")] {
                    let translated=crate::i18n::text(cols[col].ctx(),title);cols[col].label(RichText::new(translated).size(12.5).color(theme::TEXT_DIM));
                    cols[col].label(value(&notes,tag));
                }
            });
            let mut confirmed=value(&notes,"ConfirmedExperimentally")=="1";
            egui::Frame::new().fill(Color32::from_gray(236)).inner_margin(egui::Margin::same(8)).show(ui,|ui|{
                ui.set_width(ui.available_width());
                if ui.checkbox(&mut confirmed,crate::i18n::text(ui.ctx(),"Confirmed experimentally")).changed(){edits.push(("ConfirmedExperimentally".into(),if confirmed {"1"}else{"0"}.into()));}
            });
        });
        section(ui,|ui|{
            let mut kind=value(&notes,"Type");if kind.is_empty(){kind="Synthetic".into();}
            egui::ComboBox::from_id_salt("info_sequence_type").selected_text(format!("{} {}",crate::i18n::text(ui.ctx(),&kind),crate::i18n::text(ui.ctx(),doc.kind().display_name()))).show_ui(ui,|ui|{
                for k in ["Synthetic","Natural"] {if ui.selectable_value(&mut kind,k.to_string(),crate::i18n::text(ui.ctx(),k)).changed(){edits.push(("Type".into(),kind.clone()));}}
            });
            text_field(ui,&identity,"Organism","Laboratory Host",&value(&notes,"Organism"),false,&mut edits);
            text_field(ui,&identity,"TransformedInto","Bacterial Transformation Strain",&value(&notes,"TransformedInto"),false,&mut edits);
            if doc.is_dna(){ui.label(RichText::new(crate::i18n::text(ui.ctx(),"Methylation:")).size(12.5).color(theme::TEXT_DIM));
                let flags=doc.session.flags();
                ui.horizontal(|ui|{for (name,bit) in [("Dam",hx_core::MoleculeFlags::DAM_METHYLATED),("Dcm",hx_core::MoleculeFlags::DCM_METHYLATED),("EcoKI",hx_core::MoleculeFlags::ECOKI_METHYLATED)] {
                    ui.label(format!("{name}{}",if flags.has(bit){"⁺"}else{"⁻"}));
                }});
                ui.add_enabled(false,egui::Button::new(crate::i18n::text(ui.ctx(),"Change…"))).on_hover_text(crate::i18n::text(ui.ctx(),"Methylation editing is not implemented"));}
        });
        section(ui,|ui|{
            // Keep unknown/HTML description formatting intact unless the user
            // explicitly changes it; plain empty/simple values are editable.
            let description=value(&notes,"Description");
            if description.contains('<') {ui.label(crate::i18n::text(ui.ctx(),"Description"));ui.add_enabled_ui(false,|ui|{let mut html=description.clone();ui.add(egui::TextEdit::multiline(&mut html).desired_rows(4));});ui.small(crate::i18n::text(ui.ctx(),"Rich-text description retained; editing is not yet supported."));}
            else{text_field(ui,&identity,"Description","Description",&description,true,&mut edits);}
            ui.columns(2,|cols|{
                text_field(&mut cols[0],&identity,"AccessionNumber","Accession Number",&value(&notes,"AccessionNumber"),false,&mut edits);
                let title=crate::i18n::text(cols[1].ctx(),"Code Number");cols[1].label(RichText::new(title).size(12.5).color(theme::TEXT_DIM));
                let tip=crate::i18n::text(cols[1].ctx(),"Native Code Number field mapping has not been verified");cols[1].add_enabled_ui(false,|ui|{let mut code=value(&notes,"CodeNumber");ui.add_sized([ui.available_width(),30.0],egui::TextEdit::singleline(&mut code));}).response.on_hover_text(tip);
            });
            text_field(ui,&identity,"CreatedBy","Sequence Author",&value(&notes,"CreatedBy"),false,&mut edits);
            text_field(ui,&identity,"Comments","Comments",&value(&notes,"Comments"),true,&mut edits);
        });
        section(ui,|ui|{
            ui.collapsing("Map Label",|ui|{
            let mut custom=doc.custom_map_label().is_some();
            if ui.checkbox(&mut custom,crate::i18n::text(ui.ctx(),"Use a custom map label")).changed(){
                if custom {edits.push(("CustomMapLabel".into(),doc.map_name()));}
                edits.push(("UseCustomMapLabel".into(),if custom {"1"}else{"0"}.into()));
            }
            if custom {text_field(ui,&identity,"CustomMapLabel","Custom Map Label",&doc.map_name(),false,&mut edits);}
            else{ui.label(crate::ui_format!(ui.ctx(),"Use file name: {}",doc.display_name()));}
            });
        });
        section(ui,|ui|{
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),"References")).color(theme::ACCENT));
            ui.add_enabled_ui(false,|ui|{let mut refs=value(&notes,"References");ui.add_sized([ui.available_width(),180.0],egui::TextEdit::multiline(&mut refs).hint_text(crate::i18n::text(ui.ctx(),"Reference editing is not implemented.")));});
            ui.label(crate::i18n::text(ui.ctx(),"Embedded Files"));
            ui.add_enabled(false,egui::Button::new(crate::i18n::text(ui.ctx(),"Drag & drop or Browse to embed files"))).on_hover_text(crate::i18n::text(ui.ctx(),"Embedded-file editing is not implemented"));
        });
    });
    for (tag,value) in edits {doc.session.set_note(&tag,&value)?;}
    Ok(())
}
