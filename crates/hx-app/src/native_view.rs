//! Real stored alignment and gel contents, with original-byte save preservation.
use eframe::egui;
pub fn show(ui:&mut egui::Ui,bytes:&[u8])->Result<(),String>{
    let file=hx_formats::DnaFile::parse(bytes).map_err(|e|e.to_string())?;
    if matches!(file.header.type_of_sequence,3|4){
        let a=hx_formats::alignment::parse(&file)?;
        ui.heading(crate::i18n::text(ui.ctx(),if a.protein{"Protein Alignment"}else{"DNA Alignment"}));ui.label(crate::ui_format!(ui.ctx(),"{} · {} sequences · {} columns",a.method,a.rows.len(),a.rows[0].sequence.len()));
        ui.label(crate::i18n::text(ui.ctx(),"Stored alignment · read-only · saving preserves the original file"));
        egui::ScrollArea::both().id_salt("native_alignment").show(ui,|ui|{
            let width=a.rows[0].sequence.len();
            for start in (0..width).step_by(80){let end=(start+80).min(width);ui.monospace(format!("{:22} {} .. {}","",start+1,end));
                for row in &a.rows{ui.horizontal(|ui|{ui.add_sized([175.0,20.0],egui::Label::new(&row.name));let mut job=egui::text::LayoutJob::default();for col in start..end{let b=row.sequence[col];let same=a.rows.iter().all(|r|r.sequence[col]==b);let mut f=egui::TextFormat::simple(egui::FontId::monospace(13.0),if b==b'-'{egui::Color32::GRAY}else{crate::theme::TEXT});f.background=if same{egui::Color32::from_rgb(207,233,250)}else if a.protein{egui::Color32::from_rgb(248,229,210)}else{crate::theme::base_letter_color(b).gamma_multiply(0.15)};job.append(&(b as char).to_string(),0.0,f);}ui.label(job);});}
                let consensus=(start..end).map(|col|{let b=a.rows[0].sequence[col];if a.rows.iter().all(|r|r.sequence[col]==b){'*'}else{' '}}).collect::<String>();ui.monospace(format!("{:22} {}","Conservation",consensus));ui.add_space(15.0);
            }
        });return Ok(());
    }
    if file.header.type_of_sequence==8{return gel(ui,&file);}
    Err("Unsupported native document type".into())
}
fn attrs(tag:&str,key:&str)->String{let pat=format!("{key}=\"");tag.find(&pat).and_then(|s|{let r=&tag[s+pat.len()..];r.find('"').map(|e|r[..e].to_owned())}).unwrap_or_default()}
fn tags<'a>(s:&'a str,name:&str)->Vec<&'a str>{let pat=format!("<{name} ");let mut out=Vec::new();let mut at=0;while let Some(i)=s[at..].find(&pat){let start=at+i;let Some(end)=s[start..].find('>')else{break};out.push(&s[start..start+end+1]);at=start+end+1;}out}
fn gel(ui:&mut egui::Ui,file:&hx_formats::DnaFile)->Result<(),String>{
    let summary=std::str::from_utf8(&file.block(33).ok_or("Missing gel summary")?.payload).map_err(|e|e.to_string())?;
    ui.heading(crate::i18n::text(ui.ctx(),"Agarose Gel"));ui.label(crate::ui_format!(ui.ctx(),"Agarose: {}% · Stored lane definitions",attrs(summary,"percentAgarose")));ui.label(crate::i18n::text(ui.ctx(),"Sequence-derived band sizes; migration rendering is approximate, not a validated electrophoresis prediction. Read-only."));
    let inputs=tags(summary,"InputFile");let mut lanes:Vec<(String,Vec<usize>)>=Vec::new();
    let markers=tags(summary,"Fragment").iter().filter_map(|t|attrs(t,"length").parse().ok()).collect();lanes.push(("Marker".into(),markers));
    for lane in tags(summary,"TemplateLane"){
        let id=attrs(lane,"segId");let input=inputs.iter().find(|t|attrs(t,"segId")==id);let name=input.map(|t|attrs(t,"name")).unwrap_or_else(||format!("Lane {}",attrs(lane,"number")));
        let block=file.blocks.iter().find(|b|b.block_type==23&&b.payload.len()>=4&&u32::from_be_bytes(b.payload[..4].try_into().unwrap()).to_string()==id).ok_or("Missing gel embedded sequence")?;
        let source=hx_formats::DnaFile::parse(&block.payload[4..]).map_err(|e|e.to_string())?;let rec=source.sequence().ok_or("Missing gel sequence")?;
        let from=summary.find(lane).unwrap();let end=summary[from..].find("</TemplateLane>").map(|e|from+e).unwrap_or(summary.len());let body=&summary[from..end];
        let sizes=if let Some(pcr)=tags(body,"AmplificationInputs").first(){let range=attrs(pcr,"twoPrimerSelectionRange");let pair:Vec<_>=range.split('-').filter_map(|v|v.parse::<usize>().ok()).collect();if pair.len()==2{vec![if pair[1]>=pair[0]{pair[1]-pair[0]+1}else{rec.len()-pair[0]+pair[1]+1}]}else{return Err("Unsupported PCR lane definition".into());}}
        else{let mut cuts=Vec::new();if let Some(digest)=tags(body,"InputDigestions").first(){for key in ["name1","name2"]{let name=attrs(digest,key);if let Some(e)=hx_enzymes::by_name(&name){cuts.extend(hx_algo::cloning::restriction_cuts(&rec.sequence,rec.flags.topology().is_circular(),e).iter().map(|c|c.top));}}}cuts.sort();cuts.dedup();if cuts.is_empty(){vec![rec.len()]}else{let mut sizes:Vec<_>=cuts.windows(2).map(|w|(w[1]-w[0])as usize).collect();if rec.flags.topology().is_circular(){sizes.push((rec.len()as i64-cuts.last().unwrap()+cuts[0])as usize);}else{sizes.push(cuts[0]as usize);sizes.push((rec.len()as i64-cuts.last().unwrap())as usize);}sizes}};
        lanes.push((name,sizes));
    }
    let(r,_)=ui.allocate_exact_size(egui::vec2(ui.available_width().max(400.0),520.0),egui::Sense::hover());let p=ui.painter();p.rect_filled(r,2.0,egui::Color32::from_rgb(24,27,34));let width=r.width()/lanes.len().max(1)as f32;
    for(i,(name,sizes))in lanes.iter().enumerate(){let x=r.left()+width*(i as f32+0.5);p.text(egui::pos2(x,r.top()+18.0),egui::Align2::CENTER_CENTER,name,egui::FontId::proportional(12.0),egui::Color32::WHITE);for &size in sizes{if size==0{continue;}let y=r.top()+55.0+(1.0-((size as f32).log10()-2.0)/2.2).clamp(0.0,1.0)*420.0;p.line_segment([egui::pos2(x-width*0.3,y),egui::pos2(x+width*0.3,y)],egui::Stroke::new(3.0,egui::Color32::from_rgb(150,240,180)));p.text(egui::pos2(x,y+5.0),egui::Align2::CENTER_TOP,format!("{size} bp"),egui::FontId::proportional(10.0),egui::Color32::WHITE);}}
    Ok(())
}

#[cfg(test)]mod tests{
    #[test]fn native_sample_views_render_without_overwriting_sources(){
        let root=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile");
        for path in ["DNA Alignment/Ape COX1 Alignment.dnaln","DNA Alignment/Fluorescent Protein Gene Alignment.dnaln","Protein Alignment/Fluorescent Protein Alignment.praln","Gel Files/Agarose Gel.gel"]{
            let path=root.join(path);if !path.exists(){continue;}let bytes=std::fs::read(&path).unwrap();
            let context=eframe::egui::Context::default();let mut result=None;
            let mut output=context.run_ui(eframe::egui::RawInput{screen_rect:Some(eframe::egui::Rect::from_min_size(eframe::egui::Pos2::ZERO,eframe::egui::vec2(1800.0,1100.0))),..Default::default()},|ui|{result=Some(super::show(ui,&bytes));});output.textures_delta.clear();
            result.unwrap().unwrap();assert_eq!(std::fs::read(path).unwrap(),bytes);
        }
    }
}
