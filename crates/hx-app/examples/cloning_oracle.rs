//! Reproducible HX-only visual fixture; no SG process or user configuration.
//! Run with a private X display and XDG_CONFIG_HOME, passing the synthetic
//! EcoRI/BamHI vector.dna and insert.dna as arguments. Clone opens a new tab;
//! no file is overwritten automatically.
use hx_app::clonedlg::{CloningDlg,doc_from_path};
use hx_algo::cloning::Method;

fn main()->eframe::Result<()> {
    let paths=std::env::args_os().skip(1).map(std::path::PathBuf::from).collect::<Vec<_>>();
    assert_eq!(paths.len(),2,"pass synthetic vector.dna and insert.dna");
    let docs=paths.iter().map(|p|doc_from_path(p).expect("native cloning source")).collect::<Vec<_>>();
    let mut dlg=CloningDlg::new(Method::Restriction,docs,Some(0));
    let enzyme=|name|hx_enzymes::ENZYMES.iter().position(|e|e.name==name).unwrap();
    dlg.vector.doc=Some(0);dlg.vector.enzyme=Some(enzyme("EcoRI"));dlg.vector.enzyme2=Some(enzyme("BamHI"));
    dlg.fragments[0].doc=Some(1);dlg.fragments[0].enzyme=Some(enzyme("EcoRI"));dlg.fragments[0].enzyme2=Some(enzyme("BamHI"));
    let options=eframe::NativeOptions{viewport:eframe::egui::ViewportBuilder::default()
        .with_inner_size([1800.0,1100.0]).with_title("HX cloning oracle"),..Default::default()};
    eframe::run_native("HX cloning oracle",options,Box::new(move|cc|{
        let mut app=hx_app::HelixApp::new(cc,paths);app.cloning=Some(dlg);Ok(Box::new(app))
    }))
}
