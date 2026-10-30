//! Optional native SG-saved oracle check. No GUI or SG process is launched.
//! HELIX_SG_CLONING_ORACLES points at the controlled synthetic oracle folder,
//! never at a user's general documents or SG configuration directory.
use std::path::PathBuf;

#[test]
fn sg_saved_cloning_products_open_and_save_byte_exact_in_hx() {
    let Some(root)=std::env::var_os("HELIX_SG_CLONING_ORACLES").map(PathBuf::from) else {
        eprintln!("Skipped native SG files: set HELIX_SG_CLONING_ORACLES; sequence-only oracle tests still run.");return;
    };
    for (name,stored,top,circular) in [("restriction-product.dna",952,952,true),
        ("restriction-blunt-product.dna",960,960,true),("ligation-product.dna",166,163,false)] {
        let bytes=std::fs::read(root.join(name)).unwrap();
        let file=hx_formats::DnaFile::parse(&bytes).unwrap();
        assert_eq!(file.sequence().unwrap().sequence.len(),stored);
        let doc=hx_app::doc::Document::from_bytes(bytes.clone()).unwrap();
        assert_eq!(doc.len(),top);assert_eq!(doc.circular(),circular);
        assert_eq!(doc.session.save_bytes().unwrap(),bytes,"{name}: retain SG history and unknown blocks");
        if !circular {
            let ends=doc.end_properties();
            assert_eq!((ends.upstream_stickiness,ends.downstream_stickiness),(-1,2));
        }
        eprintln!("SG → HX byte-exact: {name}, {top} bp, {} blocks",file.blocks.len());
    }
}
