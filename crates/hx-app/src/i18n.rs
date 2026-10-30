//! UI-only translations. Never apply to sequence data, native XML or user names.
use eframe::egui;
#[derive(Clone,Copy,Debug,Default,PartialEq,Eq)]
pub enum Language{#[default] English,SimplifiedChinese}
impl Language{
    pub fn label(self)->&'static str{match self{Self::English=>"English",Self::SimplifiedChinese=>"简体中文"}}
    pub fn get(ctx:&egui::Context)->Self{ctx.data(|d|d.get_temp::<Language>(egui::Id::new("ag-ui-language"))).unwrap_or_default()}
    pub fn code(self)->&'static str{match self{Self::English=>"en",Self::SimplifiedChinese=>"zh-CN"}}
    pub fn decode(code:&str)->Self{if code=="zh-CN"{Self::SimplifiedChinese}else{Self::English}}
    pub fn install(self,ctx:&egui::Context){ctx.data_mut(|d|d.insert_temp(egui::Id::new("ag-ui-language"),self));}
}
pub fn text<'a>(ctx:&egui::Context,en:&'a str)->&'a str{
    if ctx.data(|d|d.get_temp::<Language>(egui::Id::new("ag-ui-language"))).unwrap_or_default()==Language::English{return en;}
    translate(en)
}
pub fn language_text(language:Language,en:&str)->&str{if language==Language::English{en}else{translate(en)}}
pub fn history_caption(ctx:&egui::Context,en:&str)->String{
    if Language::get(ctx)==Language::English{return en.into();}
    en.lines().map(|line|{let lower=line.to_ascii_lowercase();let title=lower.chars().next().map(|c|format!("{}{}",c.to_ascii_uppercase(),&lower[c.len_utf8()..])).unwrap_or_default();let translated=text(ctx,&title);if translated!=title{translated.to_owned()}else{text(ctx,line).to_owned()}}).collect::<Vec<_>>().join("\n")
}
fn catalog()->&'static std::collections::HashMap<&'static str,&'static str>{
    static CATALOG:std::sync::OnceLock<std::collections::HashMap<&'static str,&'static str>>=std::sync::OnceLock::new();
    CATALOG.get_or_init(||include_str!("zh_cn.tsv").lines().filter(|l|!l.starts_with('#')).filter_map(|l|l.split_once('\t')).collect())
}
/// Translate a UI-owned formatted message without inspecting interpolated data.
/// Callers pass the English template, never an arbitrary document string.
#[macro_export]
macro_rules! ui_format {
    ($ctx:expr,$en:literal,$($arg:expr),+ $(,)?)=>{{
        let values=[$(format!("{}",$arg)),+];
        $crate::i18n::format_values($ctx,$en,&values)
    }};
}
pub fn format_values(ctx:&egui::Context,en:&str,values:&[String])->String{
    let template=text(ctx,en);let mut result=String::new();let mut rest=template;
    for value in values{if let Some((before,after))=rest.split_once("{}"){result.push_str(before);result.push_str(value);rest=after;}}
    result.push_str(rest);result
}
/// Status/error strings produced by AG, not names or sequence annotations.
/// Known prefixes are localized without altering the trailing diagnostic/path.
pub fn message(ctx:&egui::Context,en:&str)->String{
    if Language::get(ctx)==Language::English{return en.to_owned();}
    if en.contains('\n'){return en.lines().map(|l|message(ctx,l)).collect::<Vec<_>>().join("\n");}
    let exact=text(ctx,en);if exact!=en{return exact.to_owned();}
    for prefix in ["Open failed: ","Import failed: ","Save failed: ","Export failed: ","Edit failed: ","Clone failed: ","Feature edit failed: ","Primer edit failed: ","Info edit failed: ","Cannot save preferences: ","Undo: ","Redo: ","Renamed to ","Imported ","Exported ","Selected: "]{
        if let Some(value)=en.strip_prefix(prefix){return format!("{}{}",text(ctx,prefix),text(ctx,value));}
    }
    if let Some(value)=en.strip_prefix("Product: ").and_then(|s|s.strip_suffix(" bp")){return format_values(ctx,"Product: {} bp",&[value.into()]);}
    // Only UI-owned templates participate; never call this on document fields.
    let mut templates=catalog().iter().filter(|(k,_)|k.contains("{}")&&!k.starts_with("{}")).collect::<Vec<_>>();
    templates.sort_by(|a,b|b.0.len().cmp(&a.0.len()).then(a.0.cmp(b.0)));
    for (template,_) in templates{if let Some(values)=template_values(template,en){return format_values(ctx,template,&values);}}
    en.to_owned()
}
fn template_values(template:&str,value:&str)->Option<Vec<String>>{
    let parts=template.split("{}").collect::<Vec<_>>();let mut rest=value.strip_prefix(parts[0])?;let mut values=Vec::new();
    for (i,part)in parts.iter().enumerate().skip(1){if i==parts.len()-1{let last=rest.strip_suffix(part)?;values.push(last.to_string());return Some(values);}if part.is_empty(){return None;}let at=rest.find(part)?;values.push(rest[..at].to_owned());rest=&rest[at+part.len()..];}None
}
pub fn viewer_summary(ctx:&egui::Context,en:&str)->String{
    if Language::get(ctx)==Language::English{return en.into();}
    en.replace("Text · ","文本 · ").replace(" line(s)"," 行").replace(" · modified"," · 已修改").replace(" · read-only preview"," · 只读预览").replace("PDF · page ","PDF · 第 ").replace(" of "," / ")
}
fn translate(en:&str)->&str{match en{
    "File"=>"文件","Edit"=>"编辑","View"=>"视图","Enzymes"=>"酶","Features"=>"特征","Primers"=>"引物","Actions"=>"操作","Tools"=>"工具","Window"=>"窗口","Help"=>"帮助",
    "Map"=>"图谱","Sequence"=>"序列","Overview"=>"概览","History"=>"历史","Properties"=>"属性","Secondary Structure"=>"二级结构","Display"=>"显示",
    "Save"=>"保存","Undo"=>"撤销","Redo"=>"重做","Cut"=>"剪切","Copy"=>"复制","Paste"=>"粘贴","Print"=>"打印","Clear"=>"清除","Cancel"=>"取消","Close"=>"关闭","Load"=>"加载","None"=>"无",
    "Open File…"|"Open File..."=>"打开文件…","Open Folder…"|"Open Folder..."=>"打开文件夹…","Open Recent"=>"最近打开","Save As…"|"Save As..."=>"另存为…","Save All"=>"全部保存","Close All"=>"全部关闭","Quit"=>"退出","Exit"=>"退出","Preferences…"|"Preferences..."=>"偏好设置…","Preferences"=>"偏好设置","Select All"=>"全选","Find"=>"查找","Export"=>"导出","Import"=>"导入",
    "General"=>"常规","Files & Tabs"=>"文件与标签页","Collections"=>"集合","Proteins"=>"蛋白质","Agarose Gels"=>"琼脂糖凝胶","Traces & Alignment"=>"峰图与比对",
    "Open windows full screen by default"=>"默认全屏打开窗口","Show the Info Panel by default for imported files"=>"导入文件时默认显示信息面板","Light Background Color:"=>"浅色背景：","Sequence font scale at startup (AppimageGene):"=>"启动时的序列字体大小：","Small"=>"小","Normal"=>"正常","Large"=>"大",
    "Snow"=>"雪白","Mist"=>"薄雾","Overcast"=>"阴天","Frost"=>"霜蓝","Sunny"=>"暖阳",
    "AppimageGene has no automatic update service, usage telemetry, account service or third-party account integration. No network requests are made by this window."=>"AppimageGene 没有自动更新、使用情况遥测、账户服务或第三方账户集成。此窗口不会发起网络请求。",
    "Project and File Tabs"=>"项目与文件标签页","When Starting AppimageGene:"=>"启动 AppimageGene 时：","Start with an empty workspace"=>"打开空白工作区","Continue where you left off"=>"恢复上次工作区","Default location for opening and saving files:"=>"打开和保存文件的默认位置：","Default Author for new files:"=>"新文件的默认作者：",
    "Restores saved file tabs and active tab. Unsaved changes still require saving before exit."=>"恢复已保存的文件标签页及活动标签页。退出前仍须保存未保存的修改。","Text files remain editable text documents. PDF files remain read-only."=>"文本文件保持可编辑；PDF 文件保持只读。",
    "FOLDER"=>"文件夹","OTHER OPEN FILES"=>"其他打开的文件","No folder opened"=>"尚未打开文件夹","Project files and subfolders will be shown here."=>"项目文件和子文件夹将显示在此处。","filter open files"=>"筛选打开的文件",
    "Molecular Weight Marker"=>"分子量标记","Embedded sequences"=>"内嵌序列","Open / recent DNA files"=>"打开或最近使用的 DNA 文件","Browse DNA File…"=>"浏览 DNA 文件…","Apply Enzymes to All Lanes"=>"将酶应用于所有泳道","cut with"=>"使用以下酶切割","lanes"=>"条泳道","% agarose"=>"% 琼脂糖","Move:"=>"移动：","DNA file:"=>"DNA 文件：","Stored PCR product"=>"已保存的 PCR 产物","Select a sample lane to inspect its embedded DNA."=>"选择样品泳道以查看其内嵌 DNA。",
    "Edits are saved to the native gel; migration is approximate."=>"修改可保存至原生凝胶文件；迁移位置为近似模拟。","Approximate migration · uncut circular lanes retain topology labels"=>"近似迁移模拟 · 未切环状 DNA 保留拓扑标记",
    _=>catalog().get(en).copied().unwrap_or(en),
}}
#[cfg(test)]mod catalog_tests{
    #[test]fn translations_have_unique_keys_and_preserve_format_slots(){
        let mut keys=std::collections::HashSet::new();
        for(line_no,line)in include_str!("zh_cn.tsv").lines().enumerate(){if line.starts_with('#')||line.is_empty(){continue;}let(en,zh)=line.split_once('\t').unwrap_or_else(||panic!("invalid translation line {}",line_no+1));assert!(keys.insert(en),"duplicate translation: {en}");assert!(!zh.is_empty());assert_eq!(en.matches("{}").count(),zh.matches("{}").count(),"format arguments: {en}");}
    }
    #[test]fn formatted_translation_does_not_translate_interpolated_names(){let ctx=eframe::egui::Context::default();super::Language::SimplifiedChinese.install(&ctx);assert_eq!(super::format_values(&ctx,"except at {} ends",&["DNA".into()]),"DNA 末端除外");}
    #[test]fn every_literal_menu_command_has_a_translation(){
        let source=include_str!("app.rs");let source=source.split("#[cfg(test)]").next().unwrap();
        let mut missing=Vec::new();
        for prefix in ["item(ui, \"","item_checked(ui, \""]{for rest in source.split(prefix).skip(1){let Some((label,_))=rest.split_once('"')else{continue};if !["DNA…","PCR…"].contains(&label)&&!label.contains('\\')&&label.chars().any(|c|c.is_ascii_alphabetic())&&super::translate(label)==label{missing.push(label);}}}
        assert!(missing.is_empty(),"Untranslated menu entries: {missing:?}");
    }
    #[test]fn chinese_status_keeps_diagnostic_paths_untouched(){let ctx=eframe::egui::Context::default();super::Language::SimplifiedChinese.install(&ctx);assert_eq!(super::message(&ctx,"Open failed: /data/My DNA/EcoRI.dna"),"打开失败：/data/My DNA/EcoRI.dna");super::Language::English.install(&ctx);assert_eq!(super::message(&ctx,"Open failed: /data/My DNA/EcoRI.dna"),"Open failed: /data/My DNA/EcoRI.dna");}
    #[test]fn dynamic_status_translation_retains_names_and_numbers(){let ctx=eframe::egui::Context::default();super::Language::SimplifiedChinese.install(&ctx);assert_eq!(super::message(&ctx,"Saved Map.dna (103 bytes)"),"已保存 Map.dna（103 字节）");assert_eq!(super::message(&ctx,"Ready to clone\nProduct: 219 bp"),"可以克隆\n产物：219 碱基对");}
    #[test]fn all_explicit_literal_translation_calls_are_catalogued(){
        fn scan(directory:&std::path::Path,missing:&mut Vec<String>){for entry in std::fs::read_dir(directory).unwrap(){let path=entry.unwrap().path();if path.is_dir(){scan(&path,missing);continue;}if path.extension().and_then(|s|s.to_str())!=Some("rs")||path.file_name().unwrap()=="i18n.rs"{continue;}let source=std::fs::read_to_string(&path).unwrap();
            for prefix in ["i18n::text(ctx,\"","i18n::text(ui.ctx(),\"","i18n::text(&ctx,\""]{for rest in source.split(prefix).skip(1){let Some((label,_))=rest.split_once('"')else{continue};if label.contains('\\'){continue;}if super::translate(label)==label&&!super::catalog().contains_key(label)&&label.chars().any(|c|c.is_ascii_alphabetic()){missing.push(format!("{}: {label}",path.display()));}}}
        }}
        let mut missing=Vec::new();scan(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),&mut missing);assert!(missing.is_empty(),"Missing translation entries: {missing:?}");
    }
}
#[cfg(test)]mod tests{use super::*;#[test]fn language_is_context_local_and_leaves_user_content_intact(){let a=egui::Context::default();let b=egui::Context::default();Language::SimplifiedChinese.install(&a);assert_eq!(text(&a,"Map"),"图谱");assert_eq!(text(&b,"Map"),"Map");assert_eq!(text(&a,"pET28b-Rep1 / EcoRI"),"pET28b-Rep1 / EcoRI");Language::English.install(&a);assert_eq!(text(&a,"Map"),"Map");}}
