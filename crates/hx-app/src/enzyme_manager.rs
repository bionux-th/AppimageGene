//! Restriction enzyme browser backed by the actual catalog, not a mock list.
use eframe::egui;
#[derive(Default)]pub struct Manager{pub open:bool,pub selected:usize,pub query:String,pub supplier:usize,pub version:usize}
impl Manager{
    pub fn show(&mut self,ctx:&egui::Context){
        if !self.open{return;}let mut open=true;
        egui::Window::new(crate::i18n::text(ctx,"Restriction Enzymes")).open(&mut open).default_size([760.0,650.0]).show(ctx,|ui|{
            ui.horizontal_top(|ui|{
                ui.vertical(|ui|{ui.set_width(145.0);ui.label(crate::ui_format!(ctx,"{} enzymes",hx_enzymes::catalog::CATALOG.len()));ui.text_edit_singleline(&mut self.query);egui::ScrollArea::vertical().id_salt("enzyme_manager_names").max_height(540.0).show(ui,|ui|{for(i,e)in hx_enzymes::catalog::CATALOG.iter().enumerate(){if !e.name.to_ascii_lowercase().contains(&self.query.to_ascii_lowercase()){continue;}if ui.selectable_label(self.selected==i,e.name).clicked(){self.selected=i;self.supplier=0;}}});});
                ui.separator();ui.vertical(|ui|{
                    let e=&hx_enzymes::catalog::CATALOG[self.selected];ui.heading(e.name);
                    egui::Frame::new().fill(egui::Color32::from_rgb(255,239,195)).inner_margin(12.0).show(ui,|ui|{
                        ui.monospace(format!("5′  {}  3′",e.site));
                        for pattern in e.site.split_whitespace(){ui.monospace(format!("3′  {}  5′",String::from_utf8_lossy(&hx_algo::reverse_complement(pattern.as_bytes()).into_iter().rev().collect::<Vec<_>>())));}
                        for(a,b)in e.cuts{ui.label(crate::ui_format!(ctx,"Top cut: {}    Bottom cut: {}",a.map(|v|v.to_string()).unwrap_or_else(||crate::i18n::text(ctx,"uncut").into()),b.map(|v|v.to_string()).unwrap_or_else(||crate::i18n::text(ctx,"uncut").into())));}
                    });
                    if let Some(simple)=hx_enzymes::by_name(e.name){ui.label(crate::ui_format!(ctx,"Overhang type: {}",crate::i18n::text(ctx,if hx_enzymes::sets::nicking(self.selected){"Single-strand nick"}else{simple.overhang.label()})));}
                    ui.add_space(16.0);let supplier=e.suppliers.get(self.supplier).copied().unwrap_or(0);
                    egui::ComboBox::from_id_salt("enzyme_supplier").selected_text(hx_enzymes::sets::SUPPLIER_NAMES[supplier as usize]).show_ui(ui,|ui|{for(i,&id)in e.suppliers.iter().enumerate(){ui.selectable_value(&mut self.supplier,i,hx_enzymes::sets::SUPPLIER_NAMES[id as usize]);}});
                    let detail=&hx_enzymes::details::DETAILS[self.selected];
                    let versions:Vec<_>=detail.versions.iter().filter(|v|v.supplier==supplier).collect();self.version=self.version.min(versions.len().saturating_sub(1));
                    egui::ComboBox::from_id_salt("enzyme_version").selected_text(versions.get(self.version).map(|v|crate::ui_format!(ctx,"Variant {}",v.variant)).unwrap_or_default()).show_ui(ui,|ui|{for(i,v)in versions.iter().enumerate(){ui.selectable_value(&mut self.version,i,format!("{}{}",crate::ui_format!(ctx,"Variant {}",v.variant),if v.discontinued{crate::i18n::text(ctx," (discontinued)")}else{""}));}});
                    if let Some(v)=versions.get(self.version){ui.label(crate::ui_format!(ctx,"Recommended buffer: {}",v.buffer));ui.label(crate::ui_format!(ctx,"Buffer activities: {}",v.activities));ui.label(crate::ui_format!(ctx,"Incubation temperature: {} °C",v.temperature));}
                    ui.add_space(18.0);ui.label(crate::i18n::text(ctx,"Same recognition sequence:"));
                    let same:Vec<_>=hx_enzymes::catalog::CATALOG.iter().filter(|other|other.site==e.site).collect();
                    egui::ComboBox::from_id_salt("same_recognition").selected_text(crate::ui_format!(ctx,"{} enzymes",same.len())).show_ui(ui,|ui|{for other in &same{if ui.button(other.name).clicked(){self.selected=hx_enzymes::catalog::CATALOG.iter().position(|x|x.name==other.name).unwrap();ui.close();}}});
                    ui.label(crate::i18n::text(ctx,"Same recognition and cutting sites:"));for other in same.iter().filter(|other|other.cuts==e.cuts){ui.label(other.name);}
                    ui.add_space(20.0);ui.strong(crate::i18n::text(ctx,"Methylation sensitivity"));
                    for(name,key)in [("Dam","damSeqs"),("Dcm","dcmSeqs"),("EcoKI","ecoKISeqs")]{let value=detail.methylation.iter().find(|(k,_)|*k==key);ui.label(format!("{name}: {}",crate::i18n::text(ctx,if value.is_some(){"Context-dependent (rules available)"}else{"No effect recorded"})));if let Some((_,rules))=value{ui.collapsing(crate::ui_format!(ctx,"{} context rules",name),|ui|{ui.label(*rules);});}}
                    if detail.variable_ends{ui.label(crate::i18n::text(ctx,"Sticky ends from different recognition sites may not be compatible."));}
                });
            });
        });self.open=open;
    }
}
