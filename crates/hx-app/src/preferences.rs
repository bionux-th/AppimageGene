//! Application preferences. Only implemented behavior is editable.
use eframe::egui;

pub const PAGES: &[&str] = &["General","Files & Tabs","Export","Collections","DNA","RNA","Proteins","Agarose Gels","Traces & Alignment","Enzymes","Primers"];

#[derive(Clone, Debug, PartialEq)]
pub struct Preferences {
    pub language: crate::i18n::Language,
    pub molecule_defaults: [MoleculeDefaults; 3],
    pub default_directory: String,
    pub default_author: String,
    pub info_imported: bool,
    pub restore_workspace: bool,
    pub map_internal_only: [bool;3],
    pub sequence_minimap: [bool;3],
    pub background: u8,
    pub export_formats: [u8;3],
    pub export_crlf: bool,
    pub locus_limit: usize,
    pub export_labels: bool,
    pub primer_small: hx_formats::annot::HybridizationParams,
    pub primer_large: hx_formats::annot::HybridizationParams,
    pub primer_large_threshold: usize,
    pub page: usize,
    pub italic_enzymes: bool,
    pub bold_unique: bool,
    pub green_blunt: bool,
    pub enzyme_lines: bool,
    pub suppliers: Vec<u8>,
    pub supplier_menus: usize,
    pub selected_supplier: usize,
    pub parent: hx_enzymes::sets::Parent,
    pub fullscreen: bool,
    pub font_scale: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeDefaults {
    pub features: bool, pub primers: bool, pub translations: bool,
    pub regions: bool, pub sites: bool, pub bonds: bool,
    pub min_orf_aa: usize, pub one_letter: bool, pub line_width: usize,
    pub structure_index: bool,
    pub require_start: bool,
    pub allow_end_orfs: bool,
    pub alternative_starts: bool,
}
impl Default for MoleculeDefaults {
    fn default()->Self {Self {features:true,primers:true,translations:false,regions:true,sites:true,bonds:true,min_orf_aa:75,one_letter:true,line_width:0,structure_index:true,require_start:true,allow_end_orfs:true,alternative_starts:false}}
}

impl Default for Preferences {
    fn default() -> Self { Self { language:Default::default(), page: 0, italic_enzymes: false, bold_unique: true, green_blunt: false,
        enzyme_lines: false, suppliers: hx_enzymes::sets::DEFAULT_SUPPLIER_ORDER.to_vec(), supplier_menus: 2,
        selected_supplier: 0, parent: hx_enzymes::sets::Parent::Nonredundant, fullscreen: false, font_scale: 1.0,
        molecule_defaults:std::array::from_fn(|_|MoleculeDefaults::default()),default_directory:String::new(),default_author:String::new(),info_imported:true,restore_workspace:true,map_internal_only:[false;3],sequence_minimap:[true;3],background:2,
        export_formats:[0;3],export_crlf:false,locus_limit:16,export_labels:true,primer_small:Default::default(),primer_large:hx_formats::annot::HybridizationParams{min_continuous_match_len:15,allow_mismatch:false,..Default::default()},primer_large_threshold:1_000_000 } }
}

impl Preferences {
    pub fn canvas_color(&self)->egui::Color32{match self.background{0=>egui::Color32::WHITE,1=>egui::Color32::from_rgb(247,247,247),2=>egui::Color32::from_rgb(242,242,242),3=>egui::Color32::from_rgb(239,246,250),_=>egui::Color32::from_rgb(255,251,236)}}
    pub fn read_workspace()->(Vec<std::path::PathBuf>,usize) {
        let Some(path)=Self::path().map(|p|p.with_file_name("workspace"))else{return(Default::default(),0)};
        Self::decode_workspace(&std::fs::read_to_string(path).unwrap_or_default())
    }
    pub fn decode_workspace(text:&str)->(Vec<std::path::PathBuf>,usize){
        let mut files=Vec::new();let mut active=0;
        for line in text.lines(){if let Some(s)=line.strip_prefix("active="){active=s.parse().unwrap_or(0);}else if let Some(s)=line.strip_prefix("file="){let p=std::path::PathBuf::from(s);if p.is_file()&&!files.contains(&p){files.push(p);}}}
        (files,active)
    }
    pub fn save_workspace(files:&[std::path::PathBuf],active:usize)->Result<(),String>{
        let path=Self::path().ok_or("No workspace path")?.with_file_name("workspace");
        std::fs::create_dir_all(path.parent().unwrap()).map_err(|e|e.to_string())?;
        let mut text=format!("active={active}\n");for p in files{let s=p.to_string_lossy();if !s.contains(['\n','\r']){text.push_str(&format!("file={s}\n"));}}
        let temp=path.with_extension("pending");std::fs::write(&temp,text).map_err(|e|e.to_string())?;std::fs::rename(temp,path).map_err(|e|e.to_string())
    }
    fn path() -> Option<std::path::PathBuf> { crate::settings::Settings::path().map(|p| p.with_file_name("preferences")) }
    pub fn load() -> Self {
        Self::path().and_then(|p| std::fs::read_to_string(p).ok()).map(|s| Self::decode(&s)).unwrap_or_default()
    }
    pub fn decode(s: &str) -> Self {
        let mut p = Self::default();
        for l in s.lines() {
            let Some((k,v))=l.split_once('=') else {continue};
            match k {
                "language"=>p.language=crate::i18n::Language::decode(v),
                "restore_workspace"=>p.restore_workspace=v=="true",
                "background"=>p.background=v.parse::<u8>().unwrap_or(2).min(4),
                "internal_labels"|"minimap"=>{let dst=if k=="internal_labels"{&mut p.map_internal_only}else{&mut p.sequence_minimap};for(i,value)in v.split(',').take(3).enumerate(){dst[i]=value=="true";}},
                "export_formats"=>{for (i,v) in v.split(',').take(3).enumerate(){p.export_formats[i]=v.parse::<u8>().unwrap_or(0).min(2);}},
                "export_crlf"=>p.export_crlf=v=="true","locus_limit"=>p.locus_limit=v.parse::<usize>().unwrap_or(16).clamp(1,100),"export_labels"=>p.export_labels=v=="true",
                "primer_threshold"=>p.primer_large_threshold=v.parse::<usize>().unwrap_or(1_000_000).max(1),
                "primer_small"|"primer_large"=>{let f:Vec<_>=v.split(',').collect();if f.len()==5{let target=if k=="primer_small"{&mut p.primer_small}else{&mut p.primer_large};target.min_continuous_match_len=f[0].parse::<usize>().unwrap_or(10).clamp(1,100);target.min_melting_temperature=f[1].parse::<f64>().ok().filter(|n|n.is_finite()).unwrap_or(40.0).clamp(-100.0,150.0);target.allow_mismatch=f[2]=="true";target.show_additional_five_prime_matches=f[3]=="true";target.minimum_five_prime_annealing=f[4].parse::<usize>().unwrap_or(15).clamp(1,100);}},
                "directory"=>p.default_directory=v.to_owned(), "author"=>p.default_author=v.to_owned(),
                "info_imported"=>p.info_imported=v=="true",
                "italic"=>p.italic_enzymes=v=="true", "bold_unique"=>p.bold_unique=v=="true", "green_blunt"=>p.green_blunt=v=="true",
                "enzyme_lines"=>p.enzyme_lines=v=="true", "fullscreen"=>p.fullscreen=v=="true",
                "font_scale"=>p.font_scale=v.parse::<f32>().ok().filter(|n|n.is_finite()).unwrap_or(1.0).clamp(0.8,1.5),
                "supplier_menus"=>p.supplier_menus=v.parse::<usize>().unwrap_or(2).min(26),
                "suppliers"=> { let list:Vec<u8>=v.split(',').filter_map(|s|s.parse().ok()).collect(); let set:std::collections::HashSet<_>=list.iter().copied().collect(); if list.len()==26 && set.len()==26 && list.iter().all(|v|*v<26) {p.suppliers=list;} },
                "parent"=>p.parent=match v {"all"=>hx_enzymes::sets::Parent::All,"neb"=>hx_enzymes::sets::Parent::Neb,"thermo"=>hx_enzymes::sets::Parent::Thermo,_=>hx_enzymes::sets::Parent::Nonredundant},
                _=>{if let Some(index)=k.strip_prefix("orf.").and_then(|s|s.parse::<usize>().ok()).filter(|&n|n<3){let f:Vec<_>=v.split(',').collect();if f.len()==3{p.molecule_defaults[index].require_start=f[0]=="true";p.molecule_defaults[index].allow_end_orfs=f[1]=="true";p.molecule_defaults[index].alternative_starts=f[2]=="true";}}
                    if let Some(index)=k.strip_prefix("molecule.").and_then(|s|s.parse::<usize>().ok()).filter(|&n|n<3){
                    let fields:Vec<_>=v.split(',').collect();if fields.len()==10 {let d=&mut p.molecule_defaults[index];d.features=fields[0]=="true";d.primers=fields[1]=="true";d.translations=fields[2]=="true";d.regions=fields[3]=="true";d.sites=fields[4]=="true";d.bonds=fields[5]=="true";d.min_orf_aa=fields[6].parse().unwrap_or(75).clamp(1,100000);d.one_letter=fields[7]=="true";d.line_width=fields[8].parse::<usize>().unwrap_or(0).min(10000);d.structure_index=fields[9]=="true";}
                }}
            }
        } p
    }
    pub fn encode(&self)->String { format!("italic={}\nbold_unique={}\ngreen_blunt={}\nenzyme_lines={}\nfullscreen={}\nfont_scale={}\nsupplier_menus={}\nsuppliers={}\nparent={}\n",self.italic_enzymes,self.bold_unique,self.green_blunt,self.enzyme_lines,self.fullscreen,self.font_scale,self.supplier_menus,self.suppliers.iter().map(|v|v.to_string()).collect::<Vec<_>>().join(","),match self.parent {hx_enzymes::sets::Parent::All=>"all",hx_enzymes::sets::Parent::Neb=>"neb",hx_enzymes::sets::Parent::Thermo=>"thermo",_=>"nonredundant"}) }
    pub fn save(&self)->Result<(),String> {
        let path=Self::path().ok_or("No preferences path")?;
        std::fs::create_dir_all(path.parent().unwrap()).map_err(|e|e.to_string())?;
        std::fs::write(path,self.encode_all()).map_err(|e|e.to_string())
    }
    pub fn encode_all(&self)->String {
        let mut s=self.encode();
        s.push_str(&format!("language={}\n",self.language.code()));
        s.push_str(&format!("background={}\n",self.background));
        s.push_str(&format!("restore_workspace={}\ninternal_labels={},{},{}\nminimap={},{},{}\n",self.restore_workspace,self.map_internal_only[0],self.map_internal_only[1],self.map_internal_only[2],self.sequence_minimap[0],self.sequence_minimap[1],self.sequence_minimap[2]));
        s.push_str(&format!("export_formats={},{},{}\nexport_crlf={}\nlocus_limit={}\nexport_labels={}\nprimer_threshold={}\n",self.export_formats[0],self.export_formats[1],self.export_formats[2],self.export_crlf,self.locus_limit,self.export_labels,self.primer_large_threshold));
        for (name,p) in [("primer_small",self.primer_small),("primer_large",self.primer_large)]{s.push_str(&format!("{name}={},{},{},{},{}\n",p.min_continuous_match_len,p.min_melting_temperature,p.allow_mismatch,p.show_additional_five_prime_matches,p.minimum_five_prime_annealing));}
        s.push_str(&format!("directory={}\nauthor={}\ninfo_imported={}\n",self.default_directory.replace(['\n','\r'],""),self.default_author.replace(['\n','\r'],""),self.info_imported));
        for (i,d) in self.molecule_defaults.iter().enumerate(){s.push_str(&format!("molecule.{i}={},{},{},{},{},{},{},{},{},{}\n",d.features,d.primers,d.translations,d.regions,d.sites,d.bonds,d.min_orf_aa,d.one_letter,d.line_width,d.structure_index));}
        for(i,d)in self.molecule_defaults.iter().enumerate(){s.push_str(&format!("orf.{i}={},{},{}\n",d.require_start,d.allow_end_orfs,d.alternative_starts));}
        s
    }
    pub fn show(&mut self,ctx:&egui::Context)->bool {
        let tr=|s|crate::i18n::text(ctx,s);
        self.language.install(ctx);
        let mut open=true;
        egui::Window::new(format!("AppimageGene {}: {}",tr("Preferences"),tr(PAGES[self.page.min(PAGES.len()-1)])))
            .id(egui::Id::unique("application_preferences")).open(&mut open).default_size([1120.0,730.0]).min_size([900.0,650.0]).collapsible(false).show(ctx,|ui| {
            ui.horizontal_top(|ui| {
                egui::Frame::new().fill(egui::Color32::WHITE).stroke(egui::Stroke::new(0.7,egui::Color32::GRAY)).show(ui,|ui|{ui.vertical(|ui| {ui.set_width(145.0);ui.set_min_height(690.0);ui.spacing_mut().item_spacing.y=0.0;for (i,name) in PAGES.iter().enumerate(){let(rect,response)=ui.allocate_exact_size(egui::vec2(145.0,20.0),egui::Sense::click());if self.page==i{ui.painter().rect_filled(rect,0.0,egui::Color32::from_rgb(58,173,222));}ui.painter().text(rect.left_center()+egui::vec2(3.0,0.0),egui::Align2::LEFT_CENTER,tr(name),egui::FontId::proportional(14.0),if self.page==i{egui::Color32::WHITE}else{crate::theme::TEXT});if response.clicked(){self.page=i;}}});});
                ui.separator();
                ui.vertical(|ui| {ui.set_min_width(690.0);ui.set_max_width(960.0);
                    match self.page {
                        0=> {
                            ui.horizontal(|ui|{ui.label("Language / 语言：");egui::ComboBox::from_id_salt("interface_language").selected_text(match self.language{crate::i18n::Language::English=>"English",crate::i18n::Language::SimplifiedChinese=>"简体中文"}).show_ui(ui,|ui|{ui.selectable_value(&mut self.language,crate::i18n::Language::English,"English");ui.selectable_value(&mut self.language,crate::i18n::Language::SimplifiedChinese,"简体中文");});});
                            self.language.install(ctx);
                            ui.separator();
                            ui.checkbox(&mut self.fullscreen,tr("Open windows full screen by default"));
                            ui.checkbox(&mut self.info_imported,tr("Show the Info Panel by default for imported files"));
                            ui.horizontal(|ui|{ui.label(tr("Light Background Color:"));egui::ComboBox::from_id_salt("light_background").selected_text(tr(["Snow","Mist","Overcast","Frost","Sunny"][self.background as usize])).show_ui(ui,|ui|{for(i,name)in ["Snow","Mist","Overcast","Frost","Sunny"].iter().enumerate(){ui.selectable_value(&mut self.background,i as u8,tr(name));}});});
                            ui.add_space(24.0);
                            ui.label(tr("Sequence font scale at startup (AppimageGene):"));
                            egui::ComboBox::from_id_salt("preference_font").selected_text(tr(if self.font_scale<1.0{"Small"}else if self.font_scale>1.1{"Large"}else{"Normal"})).show_ui(ui,|ui|{for (name,value) in [("Small",0.9),("Normal",1.0),("Large",1.2)]{ui.selectable_value(&mut self.font_scale,value,tr(name));}});
                            ui.separator();
                            ui.label(tr("AppimageGene has no automatic update service, usage telemetry, account service or third-party account integration. No network requests are made by this window."));
                        }
                        1=>{
                            ui.strong(tr("Project and File Tabs"));
                            ui.label(tr("When Starting AppimageGene:"));
                            ui.radio_value(&mut self.restore_workspace,false,tr("Start with an empty workspace"));
                            ui.radio_value(&mut self.restore_workspace,true,tr("Continue where you left off"));
                            ui.small(tr("Restores saved file tabs and active tab. Unsaved changes still require saving before exit."));
                            ui.add_space(20.0);ui.separator();
                            ui.label(tr("Default location for opening and saving files:"));ui.text_edit_singleline(&mut self.default_directory);
                            ui.add_space(16.0);ui.label(tr("Default Author for new files:"));ui.text_edit_singleline(&mut self.default_author);
                            ui.separator();ui.label(tr("Text files remain editable text documents. PDF files remain read-only."));
                        }
                        4..=6=>self.molecule_page(ui,self.page-4),
                        2=>{
                            for (i,name) in ["DNA","RNA","protein"].iter().enumerate(){ui.horizontal(|ui|{ui.label(crate::ui_format!(ctx,"Default export format for {} files:",tr(name)));egui::ComboBox::from_id_salt(("export_format",i)).selected_text(tr(["FASTA","GenBank / GenPept","Plain Text"][self.export_formats[i] as usize])).show_ui(ui,|ui|{for (v,n) in ["FASTA","GenBank / GenPept","Plain Text"].iter().enumerate(){ui.selectable_value(&mut self.export_formats[i],v as u8,tr(n));}});});}
                            ui.horizontal(|ui|{ui.label(tr("Line breaks in exported files:"));egui::ComboBox::from_id_salt("export_linebreak").selected_text(tr(if self.export_crlf{"Windows Format (CRLF)"}else{"macOS / Linux Format (LF)"})).show_ui(ui,|ui|{ui.selectable_value(&mut self.export_crlf,false,tr("macOS / Linux Format (LF)"));ui.selectable_value(&mut self.export_crlf,true,tr("Windows Format (CRLF)"));});});
                            ui.add_space(30.0);ui.label(tr("GenBank Customization:"));ui.horizontal(|ui|{ui.label(tr("Limit LOCUS name length:"));ui.add(egui::DragValue::new(&mut self.locus_limit).range(1..=100));ui.label(tr("characters"));});ui.checkbox(&mut self.export_labels,tr("Export feature names to /label qualifiers"));
                        }
                        10=>{
                            ui.strong(tr("These hybridization parameters are for newly entered or imported sequences."));ui.separator();
                            ui.label(tr("For small sequences:"));hybridization_group(ui,"small",&mut self.primer_small);
                            ui.add_space(14.0);
                            ui.horizontal(|ui|{ui.label(tr("For sequences larger than"));ui.add(egui::DragValue::new(&mut self.primer_large_threshold).range(1..=1_000_000_000));ui.label(tr("bp:"));});
                            hybridization_group(ui,"large",&mut self.primer_large);
                            ui.add_space(15.0);ui.separator();ui.add_space(10.0);ui.label(tr("For small and large sequences:"));
                            ui.horizontal(|ui|{binding_diagram(ui,2);ui.vertical(|ui|{ui.checkbox(&mut self.primer_small.show_additional_five_prime_matches,tr("Show additional binding sites that do not match at the 3′ end of the primer"));ui.horizontal(|ui|{ui.label(tr("but match at least"));ui.add(egui::DragValue::new(&mut self.primer_small.minimum_five_prime_annealing).range(1..=100));ui.label(tr("bases at the 5′ end"));});});});
                            self.primer_large.show_additional_five_prime_matches=self.primer_small.show_additional_five_prime_matches;self.primer_large.minimum_five_prime_annealing=self.primer_small.minimum_five_prime_annealing;
                            ui.add_space(25.0);ui.separator();ui.add_space(20.0);ui.strong(tr("Primer-template duplex display"));
                            ui.label(tr("Ungapped duplex search is implemented. Bulges, complex loops and minimum-free-energy structure selection are not yet modeled."));
                        }
                        9=>self.enzymes(ui),
                        _=>{
                            ui.strong(tr("This preference page requires a dedicated data model."));
                            ui.label(tr("Collection management, gel simulation and trace rendering are not implemented by this preferences window. Existing native files are preserved without modification."));
                        }
                    }
                });
            });
        });open
    }
    fn enzymes(&mut self,ui:&mut egui::Ui) {
        let ctx=ui.ctx().clone();let tr=|s|crate::i18n::text(&ctx,s);
        ui.checkbox(&mut self.italic_enzymes,tr("Show restriction enzyme names in italics"));
        ui.checkbox(&mut self.bold_unique,tr("Show unique cutters in bold"));
        ui.checkbox(&mut self.green_blunt,tr("Show blunt cutters in green"));
        ui.add_space(20.0);
        ui.horizontal(|ui|{ui.label(tr("Default tab for Enzymes view:"));ui.radio_value(&mut self.enzyme_lines,false,tr("Numbers"));ui.radio_value(&mut self.enzyme_lines,true,tr("Lines"));});
        ui.add_space(16.0);ui.separator();ui.add_space(14.0);
        ui.label(tr("Enzyme suppliers in order of preference:"));
        egui::Frame::new().fill(egui::Color32::WHITE).stroke(egui::Stroke::new(1.0,egui::Color32::LIGHT_GRAY)).show(ui,|ui|{
            egui::ScrollArea::vertical().id_salt("preferred_suppliers").max_height(210.0).show(ui,|ui|{
                ui.set_min_width(680.0);
                for (row,&id) in self.suppliers.iter().enumerate(){let count=hx_enzymes::catalog::CATALOG.iter().filter(|e|e.suppliers.contains(&id)).count();if ui.selectable_label(self.selected_supplier==row,format!("{} ({count})",hx_enzymes::sets::SUPPLIER_NAMES[id as usize])).clicked(){self.selected_supplier=row;}}
            });
        });
        ui.horizontal(|ui|{
            if ui.add_enabled(self.selected_supplier>0,egui::Button::new(tr("Move Up"))).clicked(){self.suppliers.swap(self.selected_supplier,self.selected_supplier-1);self.selected_supplier-=1;}
            if ui.add_enabled(self.selected_supplier>0,egui::Button::new(tr("Move to Top"))).clicked(){let s=self.suppliers.remove(self.selected_supplier);self.suppliers.insert(0,s);self.selected_supplier=0;}
            if ui.add_enabled(self.selected_supplier+1<self.suppliers.len(),egui::Button::new(tr("Move Down"))).clicked(){self.suppliers.swap(self.selected_supplier,self.selected_supplier+1);self.selected_supplier+=1;}
        });
        ui.horizontal(|ui|{ui.label(tr("List the top"));egui::ComboBox::from_id_salt("num_supplier_menus").selected_text(self.supplier_menus.to_string()).show_ui(ui,|ui|{for n in 0..=26 {ui.selectable_value(&mut self.supplier_menus,n,n.to_string());}});ui.label(tr("suppliers in the program menus"));});
        ui.add_space(16.0);ui.separator();ui.add_space(14.0);
        ui.horizontal_wrapped(|ui| {
            ui.label(tr("Use"));egui::ComboBox::from_id_salt("parent_preference").selected_text(tr(match self.parent{hx_enzymes::sets::Parent::All=>"All Commercial",hx_enzymes::sets::Parent::Nonredundant=>"Nonredundant Commercial",hx_enzymes::sets::Parent::Neb=>"New England Biolabs",hx_enzymes::sets::Parent::Thermo=>"Thermo Fisher (Fermentas)",_=>"Supplier"})).show_ui(ui,|ui|{
                for (label,parent) in [("All Commercial",hx_enzymes::sets::Parent::All),("Nonredundant Commercial",hx_enzymes::sets::Parent::Nonredundant),("New England Biolabs",hx_enzymes::sets::Parent::Neb),("Thermo Fisher (Fermentas)",hx_enzymes::sets::Parent::Thermo)]{ui.selectable_value(&mut self.parent,parent,tr(label));}
            });ui.label(tr("as the parental enzyme set for Unique Cutters, Unique & Dual Cutters, 6+ Cutters, Unique 6+ Cutters, Type IIS Enzymes, and Nicking Endonucleases"));
        });
    }
    fn molecule_page(&mut self,ui:&mut egui::Ui,index:usize){
        let ctx=ui.ctx().clone();let tr=|s|crate::i18n::text(&ctx,s);
        let kind=["DNA","RNA","protein"][index];
        ui.strong(crate::ui_format!(&ctx,"These settings are for newly entered or imported {} sequences.",tr(kind)));
        ui.add_space(16.0);ui.separator();ui.add_space(20.0);
        let d=&mut self.molecule_defaults[index];
        if index==2 {ui.checkbox(&mut d.regions,tr("Show \"Region\" features"));ui.checkbox(&mut d.sites,tr("Show \"Site\" and \"Het\" features"));ui.checkbox(&mut d.bonds,tr("Show \"Bond\" features"));}
        else {ui.checkbox(&mut d.features,tr("Show features"));if index==0 {ui.checkbox(&mut d.primers,tr("Show primers"));}ui.checkbox(&mut d.translations,tr("Show translations and ORFs"));ui.horizontal(|ui|{ui.label(tr("ORFs: Minimum Length:"));ui.add(egui::DragValue::new(&mut d.min_orf_aa).range(1..=100000));ui.label(tr("amino acids"));});}
        if index<2{ui.horizontal(|ui|{ui.checkbox(&mut d.require_start,tr("Require a start codon:"));egui::ComboBox::from_id_salt(("start_codons",index)).selected_text(if d.alternative_starts{"ATG, GTG, TTG"}else if index==1{"AUG"}else{"ATG"}).show_ui(ui,|ui|{ui.selectable_value(&mut d.alternative_starts,false,if index==1{"AUG"}else{"ATG"});ui.selectable_value(&mut d.alternative_starts,true,if index==1{"AUG, GUG, UUG"}else{"ATG, GTG, TTG"});});ui.checkbox(&mut d.allow_end_orfs,crate::ui_format!(&ctx,"except at {} ends",if index==1{"RNA"}else{"DNA"}));});}
        ui.add_space(35.0);
        ui.columns(2,|cols|{
            cols[0].label(tr("Map View:"));
            egui::ComboBox::from_id_salt(("label_policy",index)).selected_text(tr(if self.map_internal_only[index]{"Internal Labels Only"}else{"Internal and External Labels"})).show_ui(&mut cols[0],|ui|{ui.selectable_value(&mut self.map_internal_only[index],false,tr("Internal and External Labels"));ui.selectable_value(&mut self.map_internal_only[index],true,tr("Internal Labels Only"));});
            cols[1].label(tr("Sequence View:"));cols[1].checkbox(&mut d.one_letter,tr("Use 1-letter amino acid codes"));
            cols[1].checkbox(&mut self.sequence_minimap[index],tr("Show minimap"));
            let mut fixed=d.line_width>0;
            if cols[1].checkbox(&mut fixed,tr("Use fixed line width")).changed(){d.line_width=if fixed {60}else{0};}
            if fixed {cols[1].horizontal(|ui|{ui.add(egui::DragValue::new(&mut d.line_width).range(1..=10000));ui.label(tr(if index==2{"amino acids per line"}else{"bases per line"}));});}
        });
        if index==1 {ui.add_space(25.0);ui.label(tr("Structure View:"));ui.checkbox(&mut d.structure_index,tr("Show base index labels on structure"));}
    }
}

#[cfg(test)] mod tests {
    #[test]fn language_roundtrip_and_legacy_default(){let mut p=super::Preferences::default();assert_eq!(p.language,crate::i18n::Language::English);p.language=crate::i18n::Language::SimplifiedChinese;assert_eq!(super::Preferences::decode(&p.encode_all()),p);assert_eq!(super::Preferences::decode("language=unknown").language,crate::i18n::Language::English);}
    #[test]fn chinese_molecule_preferences_render_translated_labels(){let ctx=egui::Context::default();crate::i18n::Language::SimplifiedChinese.install(&ctx);let mut preferences=Preferences::default();let mut output=ctx.run_ui(egui::RawInput::default(),|ui|preferences.molecule_page(ui,0));output.textures_delta.clear();let labels=output.shapes.iter().filter_map(|s|match &s.shape{egui::Shape::Text(t)=>Some(t.galley.job.text.as_str()),_=>None}).collect::<Vec<_>>();assert!(labels.contains(&"显示特征"));assert!(!labels.contains(&"Show features"));}
    use super::*;
    #[test]fn all_preference_pages_render_in_both_languages(){
        for language in [crate::i18n::Language::English,crate::i18n::Language::SimplifiedChinese]{let ctx=egui::Context::default();let mut p=Preferences::default();p.language=language;for page in 0..PAGES.len(){p.page=page;for _ in 0..2{let mut out=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(1800.0,1100.0))),..Default::default()},|_|{p.show(&ctx);});out.textures_delta.clear();assert!(!out.shapes.is_empty());}}}
    }
    #[test] fn preferences_roundtrip_and_invalid_supplier_order() {
        let mut p=Preferences::default();p.suppliers.swap(0,1);p.italic_enzymes=true;p.supplier_menus=5;
        assert_eq!(Preferences::decode(&p.encode_all()),p);
        assert_eq!(Preferences::decode("suppliers=11,11\nfont_scale=NaN").suppliers,Preferences::default().suppliers);
    }
    #[test] fn all_editable_defaults_survive_restart() {
        let mut p=Preferences::default();p.default_author="Example author".into();p.default_directory="/tmp/opencode".into();p.export_formats=[1,2,1];p.export_crlf=true;p.primer_large_threshold=5000;p.primer_small.min_melting_temperature=51.0;p.molecule_defaults[1].line_width=120;p.molecule_defaults[2].bonds=false;
        assert_eq!(Preferences::decode(&p.encode_all()),p);
    }
    #[test] fn workspace_and_display_defaults_are_serialized() {
        let mut p=Preferences::default();p.restore_workspace=false;p.background=4;p.map_internal_only=[true,false,true];p.sequence_minimap=[false,true,false];
        assert_eq!(Preferences::decode(&p.encode_all()),p);
        assert_ne!(p.canvas_color(),Preferences::default().canvas_color());
    }
    #[test] fn workspace_ignores_missing_and_duplicate_files() {
        let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
        let text=format!("active=1\nfile={}\nfile={}\nfile=/does/not/exist\n",path.display(),path.display());
        let (files,active)=Preferences::decode_workspace(&text);assert_eq!(files,vec![path]);assert_eq!(active,1);
    }
    #[test] fn genbank_sequence_export_roundtrips() {
        let bytes=hx_formats::write::NewFile::dna(b"ACGT".repeat(25)).with_name("test").to_bytes();
        let doc=crate::doc::Document::from_bytes(bytes).unwrap();
        let text=genbank_export(&doc,&Preferences::default());
        let imported=hx_formats::import::parse(text.as_bytes()).unwrap();
        assert_eq!(imported.len(),1);
        let file=hx_formats::DnaFile::parse(&hx_formats::import::to_dna_bytes(&imported[0])).unwrap();
        assert_eq!(file.sequence().unwrap().sequence.iter().map(u8::to_ascii_uppercase).collect::<Vec<_>>(),doc.sequence());
    }
}

pub fn genbank_export(doc:&crate::doc::Document,prefs:&Preferences)->String {
    let name=doc.map_name().chars().map(|c|if c.is_ascii_alphanumeric()||c=='_'{c}else{'_'}).take(prefs.locus_limit).collect::<String>();
    let protein=doc.kind()==hx_core::MoleculeKind::Protein;
    let mut out=format!("LOCUS       {name:<16} {} {} {} {}\nDEFINITION  {}\nFEATURES             Location/Qualifiers\n",doc.len(),if protein{"aa"}else{"bp"},doc.kind().display_name(),if doc.circular(){"circular"}else{"linear"},doc.map_name().replace(['\n','\r']," "));
    for f in doc.features(){
        let spans=f.segments.iter().map(|s|format!("{}..{}",s.range.start+1,s.range.end)).collect::<Vec<_>>();
        if spans.is_empty(){continue;}
        let mut loc=if spans.len()==1{spans[0].clone()}else{format!("join({})",spans.join(","))};
        if f.has_reverse_arrow()&&!f.has_forward_arrow(){loc=format!("complement({loc})");}
        let kind=if f.kind.bytes().all(|b|b.is_ascii_alphanumeric()||b==b'_'){f.kind.as_str()}else{"misc_feature"};
        out.push_str(&format!("     {kind:<16}{loc}\n"));
        if prefs.export_labels{out.push_str(&format!("                     /label=\"{}\"\n",f.name.replace('"',"\"\"").replace(['\n','\r']," ")));}
    }
    out.push_str("ORIGIN\n");
    for (row,chunk) in doc.sequence().chunks(60).enumerate(){out.push_str(&format!("{:>9} {}\n",row*60+1,chunk.chunks(10).map(|c|c.iter().map(|&b|(doc.kind().display_base(b) as char).to_ascii_lowercase()).collect::<String>()).collect::<Vec<_>>().join(" ")));}
    out.push_str("//\n");out
}

fn binding_diagram(ui:&mut egui::Ui,mode:u8){
    let(r,_)=ui.allocate_exact_size(egui::vec2(100.0,52.0),egui::Sense::hover());let p=ui.painter();
    p.rect_filled(r,0.0,egui::Color32::WHITE);p.rect_stroke(r,0.0,egui::Stroke::new(0.7,egui::Color32::GRAY),egui::StrokeKind::Inside);
    let y=r.center().y+4.0;p.line_segment([egui::pos2(r.left()+12.0,y),egui::pos2(r.right()-12.0,y)],egui::Stroke::new(1.5,egui::Color32::GRAY));
    for i in 0..21{let x=r.left()+15.0+i as f32*3.3;p.line_segment([egui::pos2(x,y-3.0),egui::pos2(x,y+3.0)],egui::Stroke::new(0.5,egui::Color32::GRAY));}
    let purple=egui::Color32::from_rgb(175,45,245);let a=egui::pos2(r.left()+30.0,y-9.0);let b=egui::pos2(r.right()-25.0,y-9.0);
    p.line_segment([a,b],egui::Stroke::new(1.5,purple));p.arrow(b-egui::vec2(6.0,0.0),egui::vec2(6.0,0.0),egui::Stroke::new(1.5,purple));
    if mode==1{p.line_segment([a+egui::vec2(10.0,0.0),a+egui::vec2(14.0,-5.0)],egui::Stroke::new(1.3,purple));}
    if mode==2{p.line_segment([b,b+egui::vec2(6.0,-6.0)],egui::Stroke::new(1.3,purple));}
}
fn hybridization_group(ui:&mut egui::Ui,id:&str,p:&mut hx_formats::annot::HybridizationParams){
    ui.push_id(id,|ui|{
        ui.horizontal(|ui|{binding_diagram(ui,0);ui.vertical(|ui|{ui.horizontal(|ui|{ui.label(crate::i18n::text(ui.ctx(),"Match at least"));ui.add(egui::DragValue::new(&mut p.min_continuous_match_len).range(1..=100));ui.label(crate::i18n::text(ui.ctx(),"bases at the 3′ end of the primer"));});ui.checkbox(&mut p.allow_mismatch,crate::i18n::text(ui.ctx(),"Allowing a single isolated mismatch"));});});
        ui.horizontal(|ui|{binding_diagram(ui,1);ui.label(crate::i18n::text(ui.ctx(),"Require a melting temperature (Tm) of at least"));ui.add(egui::DragValue::new(&mut p.min_melting_temperature).range(-100.0..=150.0));ui.label("°C");});
    });
}
