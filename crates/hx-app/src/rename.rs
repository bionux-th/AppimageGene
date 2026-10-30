//! File rename, never overwrite an existing target or change file content.
use std::path::{Path,PathBuf};
use eframe::egui;

pub struct Rename { pub path:PathBuf, pub name:String, pub error:String }
impl Rename {
    pub fn new(path:PathBuf)->Self{Self{name:path.file_name().unwrap_or_default().to_string_lossy().into_owned(),path,error:String::new()}}
    pub fn show(&mut self,ctx:&egui::Context)->Option<bool>{
        let mut outcome=None;
        egui::Modal::new(egui::Id::unique("rename_file")).show(ctx,|ui|{
            ui.set_width(420.0);ui.heading(crate::i18n::text(ctx,"Rename File"));ui.label(crate::i18n::text(ctx,"File contents and extension are preserved. An existing file will not be overwritten."));
            ui.text_edit_singleline(&mut self.name);
            if !self.error.is_empty(){ui.colored_label(egui::Color32::RED,crate::i18n::message(ctx,&self.error));}
            ui.horizontal(|ui|{if ui.button(crate::i18n::text(ctx,"Rename")).clicked(){outcome=Some(true);}if ui.button(crate::i18n::text(ctx,"Cancel")).clicked(){outcome=Some(false);}});
            if ui.input(|i|i.key_pressed(egui::Key::Escape)){outcome=Some(false);}
        });outcome
    }
}
pub fn rename_file(path:&Path,name:&str)->Result<PathBuf,String>{
    if name.trim().is_empty()||name!=name.trim()||name.contains(['/', '\\','\n','\r'])||matches!(name,"."|".."){return Err("Enter a file name without directory separators or leading/trailing spaces.".into());}
    let target=path.with_file_name(name);
    if target==path{return Ok(target);}
    if path.extension()!=target.extension(){return Err("Keep the original file extension.".into());}
    if !path.is_file(){return Err("The source file no longer exists.".into());}
    // Hard-link creation is atomic/no-clobber on the same filesystem. Remove
    // the old name only after the new name exists; no byte copy or overwrite.
    std::fs::hard_link(path,&target).map_err(|e|format!("Cannot rename: {e}"))?;
    if let Err(e)=std::fs::remove_file(path){let _=std::fs::remove_file(&target);return Err(format!("Cannot remove old file name: {e}"));}
    Ok(target)
}
#[cfg(test)]mod tests{
    #[test]fn rename_preserves_bytes_and_rejects_overwrite(){
        let dir=std::path::PathBuf::from("/tmp/opencode").join(format!("ag-rename-{}",std::process::id()));std::fs::create_dir_all(&dir).unwrap();let a=dir.join("a.txt");let b=dir.join("b.txt");assert!(!a.exists()&&!b.exists());std::fs::write(&a,b"unchanged").unwrap();std::fs::write(&b,b"existing").unwrap();assert!(super::rename_file(&a,"b.txt").is_err());assert!(super::rename_file(&a,"a.dna").is_err());assert!(super::rename_file(&a,"../bad.txt").is_err());let c=super::rename_file(&a,"c.txt").unwrap();assert_eq!(std::fs::read(&c).unwrap(),b"unchanged");assert_eq!(std::fs::read(&b).unwrap(),b"existing");assert!(!a.exists());std::fs::remove_file(c).unwrap();std::fs::remove_file(b).unwrap();std::fs::remove_dir(dir).unwrap();
    }
}
