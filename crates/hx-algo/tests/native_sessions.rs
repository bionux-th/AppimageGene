//! 原生样例通过编辑会话读取/未编辑保存，必须保留全部字节及末端信息。
use std::path::{Path, PathBuf};
use hx_algo::EditSession;
use hx_formats::DnaFile;

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap().flatten() {
        let path = entry.path();
        if path.is_dir() { collect(&path, out); }
        else if matches!(path.extension().and_then(|x| x.to_str()), Some("dna" | "rna" | "prot")) { out.push(path); }
    }
}

#[test]
fn native_edit_sessions_keep_unedited_files_byte_exact_and_project_linear_ends() {
    let root = PathBuf::from(std::env::var("HELIX_SAMPLES_DIR")
        .unwrap_or_else(|_| "/opt/gslbiotech/snapgene/resources".into()));
    if !root.is_dir() {
        eprintln!("Skipping native session regression: set HELIX_SAMPLES_DIR to native samples");
        return;
    }
    let mut files = Vec::new(); collect(&root, &mut files); files.sort();
    let mut checked = 0;
    let mut sticky = 0;
    for path in files {
        let bytes = std::fs::read(&path).unwrap();
        let file = DnaFile::parse(&bytes).unwrap();
        // Alignments and gels are retained by the app's NativeOnly path, not EditSession.
        let Some(rec) = file.sequence() else { continue; };
        let s = EditSession::new(file.clone(), bytes.clone()).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert_eq!(s.save_bytes().unwrap(), bytes, "{}", path.display());
        assert_eq!(s.stored_sequence(), rec.sequence, "{}", path.display());
        let range = s.end_properties().top_range(rec.sequence.len());
        assert_eq!(s.sequence(), &rec.sequence[range.clone()], "{}", path.display());
        assert_eq!(s.len(), range.len(), "{}", path.display());
        if s.len() != rec.sequence.len() { sticky += 1; }
        checked += 1;
    }
    assert!(checked > 0, "No native sequence files found in {}", root.display());
    eprintln!("Checked {checked} native edit sessions, including {sticky} with bottom-strand end markers");
}
