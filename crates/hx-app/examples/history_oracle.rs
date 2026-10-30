//! Read-only native History visual fixture on a private display.
fn main() -> eframe::Result<()> {
    let paths = std::env::args_os().skip(1).map(std::path::PathBuf::from).collect();
    let options = eframe::NativeOptions { viewport: eframe::egui::ViewportBuilder::default()
        .with_inner_size([1800.0, 1100.0]).with_title("HX History oracle"), ..Default::default() };
    eframe::run_native("HX History oracle", options, Box::new(move |cc| {
        let mut app = hx_app::HelixApp::new(cc, paths);
        app.tab = hx_app::Tab::History;
        if std::env::var_os("AG_DENSE_ENZYMES").is_some(){
            app.tab=hx_app::Tab::Sequence;
            app.enzyme_parent=hx_enzymes::sets::Parent::All;
            app.enzyme_filter=hx_enzymes::sets::Filter::All;
            if let Some(doc)=app.doc.as_ref(){app.sites.selected=hx_enzymes::sets::select(app.enzyme_parent,app.enzyme_filter,doc.sequence(),doc.circular());app.sites_dirty=true;}
            app.show_sites=true;
        }
        if std::env::var_os("HX_PREFERENCES").is_some() {
            app.enzyme_preferences_open = true;
            app.preferences.page = 9;
            if let Ok(page)=std::env::var("HX_PREFERENCES_PAGE"){app.preferences.page=page.parse().unwrap_or(9);}
        }
        if std::env::var_os("HX_ENZYME_MENU").is_some() {
            app.tab = hx_app::Tab::Map;
            app.enzyme_parent = hx_enzymes::sets::Parent::All;
            app.enzyme_filter = hx_enzymes::sets::Filter::All;
        }
        if std::env::var_os("HX_HYBRIDIZATION_DIALOG").is_some() {
            app.tab = hx_app::Tab::Map;
            app.hybridization_draft = Some(hx_formats::annot::HybridizationParams::default());
        }
        if let Ok(action) = std::env::var("HX_MOLECULAR_ACTION") {
            use hx_app::molecular_actions::{Action, Dialog};
            let action = match action.as_str() { "linearize" => Action::Linearize, "single" => Action::SingleStrand, "anneal" => Action::Anneal, "rna" => Action::Rna, _ => Action::Protein };
            let mut dialog = Dialog::new(action, app.doc.as_ref(), None).unwrap();
            if action == Action::Anneal { dialog.oligo1 = "AAACGT".into(); dialog.oligo2 = "TTACGT".into(); dialog.offset = 2; }
            app.molecular_action = Some(dialog);
        }
        if std::env::var_os("HX_HISTORY_EDITS").is_some() {
            if let Some(doc) = app.doc.as_mut() {
                doc.session.apply(hx_formats::write::Edit::Delete { start: 5, end: 12 }).unwrap();
                doc.session.apply(hx_formats::write::Edit::Replace { start: 20, end: 30, seq: b"AAAA".to_vec() }).unwrap();
                doc.session.apply(hx_formats::write::Edit::Insert { at: 40, seq: b"GGGGG".to_vec() }).unwrap();
            }
        }
        if std::env::var_os("HX_HISTORY_DIALOG").is_some() {
            app.tab = hx_app::Tab::Sequence;
            app.selection = Some((20, 30));
            app.paste = "AATTCCGTA".into();
            app.show_paste_dialog = true;
        }
        Ok(Box::new(app))
    }))
}
