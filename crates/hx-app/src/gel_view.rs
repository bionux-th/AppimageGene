//! Split gel workspace. Embedded source viewers are independent of main tabs.
use eframe::egui::{self,Color32,Pos2,Rect,Stroke,Vec2};
use crate::{app::{HelixApp,CanvasMode},doc::Document};

#[derive(Clone)]
pub struct Lane {pub number:usize,pub name:String,pub sizes:Vec<usize>,pub treatment:String,pub supercoiled:bool,pub bytes:Option<Vec<u8>>,pub enzymes:[String;4],pub marker:Option<usize>,pub pcr:Option<String>}
#[derive(Clone)]pub struct Marker{pub supplier:String,pub name:String,pub sizes:Vec<usize>,pub definition:String}
#[derive(Clone)]pub struct Source{pub name:String,pub bytes:Vec<u8>}
/// Source coordinates are top-strand, zero-based gaps; start > end wraps the origin.
#[derive(Clone,Debug,PartialEq,Eq)]
pub struct Fragment{pub start:usize,pub end:usize,pub length:usize,pub left:String,pub right:String}
#[derive(Clone,Debug,PartialEq,Eq)]
pub struct DigestSelection{pub range:(usize,usize),pub top_start:i64,pub top_end:i64,pub bottom_start:i64,pub bottom_end:i64}
fn digest_selection(lane:&Lane,f:&Fragment)->Option<DigestSelection>{
    if lane.pcr.is_some(){return None;}
    let file=hx_formats::DnaFile::parse(lane.bytes.as_ref()?).ok()?;let rec=file.sequence()?;
    let mut ends=std::collections::BTreeMap::new();
    for name in lane.enzymes.iter().filter(|n|!n.is_empty()){let enzyme=hx_enzymes::by_name(name)?;for cut in hx_algo::cloning::restriction_cuts(&rec.sequence,rec.flags.topology().is_circular(),enzyme){
        // Two enzymes cutting the top at the same gap but the bottom differently
        // do not define one unambiguous pair of duplex ends.
        let entry=ends.entry(cut.top as usize).or_insert(Some(cut.bottom-cut.top));if *entry!=Some(cut.bottom-cut.top){*entry=None;}
    }}
    let start_offset=match ends.get(&f.start){Some(v)=>(*v)?,None=>0};let end_offset=match ends.get(&f.end){Some(v)=>(*v)?,None=>0};
    if ends.is_empty(){return None;}
    let top_start=f.start as i64;let top_end=top_start+f.length as i64;
    Some(DigestSelection{range:if f.length==rec.len(){(0,rec.len())}else{(f.start,f.end)},top_start,top_end,bottom_start:top_start+start_offset,bottom_end:top_end+end_offset})
}
impl Fragment{fn description(&self)->String{format!("{} bp = {} — {}",self.length,self.left,self.right)}}
fn lane_fragments(lane:&Lane)->Result<Vec<Fragment>,String>{
    let Some(bytes)=&lane.bytes else{return Ok(Vec::new())};
    let file=hx_formats::DnaFile::parse(bytes).map_err(|e|e.to_string())?;let rec=file.sequence().ok_or("Missing DNA")?;let n=rec.len();if n==0{return Ok(Vec::new());}
    if let Some(pcr)=&lane.pcr{let pair=attr(pcr,"twoPrimerSelectionRange").split('-').map(str::parse::<usize>).collect::<Result<Vec<_>,_>>().map_err(|_|"Invalid PCR range")?;
        if pair.len()!=2||pair.iter().any(|p|*p==0||*p>n){return Err("Invalid PCR range".into());}let start=pair[0]-1;let end=pair[1];return Ok(vec![Fragment{start,end,length:if end>start{end-start}else{n-start+end},left:format!("{} ({})",attr(pcr,"primer1"),pair[0]),right:format!("{} ({})",attr(pcr,"primer2"),pair[1])}]);}
    let circular=rec.flags.topology().is_circular();let mut cuts=std::collections::BTreeMap::<usize,Vec<String>>::new();
    for name in lane.enzymes.iter().filter(|s|!s.is_empty()){
        let i=hx_enzymes::catalog::CATALOG.iter().position(|e|e.name==*name).ok_or_else(||format!("Unknown enzyme {name}"))?;
        if hx_enzymes::sets::nicking(i)||hx_enzymes::catalog::CATALOG[i].cuts.len()!=1{return Err(format!("{name}: nicking/multiple-cut digestion is not supported by this simulator"));}
        for cut in hx_algo::cloning::restriction_cuts(&rec.sequence,circular,&hx_enzymes::ENZYMES[i]){let labels=cuts.entry(cut.top as usize).or_default();let label=format!("{name} ({})",cut.top);if !labels.contains(&label){labels.push(label);}}
    }
    if cuts.is_empty(){return Ok(vec![Fragment{start:0,end:n,length:n,left:"Start (1)".into(),right:format!("End ({n})")}]);}
    if !circular{cuts.entry(0).or_insert_with(||vec!["Start (1)".into()]);cuts.entry(n).or_insert_with(||vec![format!("End ({n})")]);}
    let cuts=cuts.into_iter().map(|(p,names)|(p,names.join(" / "))).collect::<Vec<_>>();let mut fragments=Vec::new();
    for w in cuts.windows(2){fragments.push(Fragment{start:w[0].0,end:w[1].0,length:w[1].0-w[0].0,left:w[0].1.clone(),right:w[1].1.clone()});}
    if circular{let last=cuts.last().unwrap();let first=&cuts[0];fragments.push(Fragment{start:last.0,end:first.0,length:n-last.0+first.0,left:last.1.clone(),right:first.1.clone()});}
    fragments.retain(|f|f.length>0);fragments.sort_by(|a,b|b.length.cmp(&a.length).then(a.start.cmp(&b.start)));Ok(fragments)
}
pub struct GelView {pub key:Vec<u8>,pub lanes:Vec<Lane>,pub selected:usize,pub band:Option<usize>,pub preview:Option<Box<HelixApp>>,pub max_lanes:usize,pub percent:f32,pub error:String,pub markers:Vec<Marker>,pub sources:Vec<Source>,pub available:Vec<Source>,pub changed:bool,pub browse_path:String,pub show_browse:bool}
fn empty_lane(number:usize)->Lane{Lane{number,name:"Empty".into(),sizes:Vec::new(),treatment:String::new(),supercoiled:false,bytes:None,enzymes:Default::default(),marker:None,pcr:None}}
fn attr(tag:&str,key:&str)->String{let p=format!("{key}=\"");tag.find(&p).and_then(|i|{let r=&tag[i+p.len()..];r.find('"').map(|j|r[..j].replace("&quot;","\"").replace("&apos;","'").replace("&lt;","<").replace("&gt;",">").replace("&amp;","&"))}).unwrap_or_default()}
fn tags<'a>(text:&'a str,name:&str)->Vec<&'a str>{let mut out=Vec::new();let mut pos=0;let pat=format!("<{name} ");while let Some(i)=text[pos..].find(&pat){let at=pos+i;let Some(end)=text[at..].find('>')else{break};out.push(&text[at..at+end+1]);pos=at+end+1;}out}
fn popup_style(ui:&mut egui::Ui){
    ui.spacing_mut().item_spacing.y=0.0;ui.spacing_mut().interact_size.y=20.0;ui.spacing_mut().button_padding=Vec2::new(6.0,1.0);
    ui.spacing_mut().scroll=egui::style::ScrollStyle::solid();ui.spacing_mut().scroll.bar_width=8.0;
    let visuals=ui.visuals_mut();for w in [&mut visuals.widgets.inactive,&mut visuals.widgets.hovered,&mut visuals.widgets.active]{w.bg_stroke=Stroke::NONE;w.corner_radius=egui::CornerRadius::ZERO;}
    visuals.widgets.inactive.bg_fill=Color32::TRANSPARENT;visuals.widgets.inactive.weak_bg_fill=Color32::TRANSPARENT;visuals.widgets.hovered.weak_bg_fill=Color32::from_rgb(215,233,248);
}
fn menu_row(ui:&mut egui::Ui,text:&str,selected:bool)->egui::Response{
    let (rect,response)=ui.allocate_exact_size(Vec2::new(ui.available_width(),20.0),egui::Sense::click());
    if selected||response.hovered()||response.has_focus(){ui.painter().rect_filled(rect,0.0,Color32::from_rgb(215,233,248));}
    ui.painter().with_clip_rect(rect).text(rect.left_center()+Vec2::new(6.0,0.0),egui::Align2::LEFT_CENTER,text,egui::TextStyle::Button.resolve(ui.style()),crate::theme::TEXT);
    response.widget_info(||egui::WidgetInfo::selected(egui::WidgetType::Button,ui.is_enabled(),selected,text));response
}
fn approximate_band_fraction(n:usize,percent:f32)->f32{
    // Keep sub-100 bp products separate. This remains an explicitly approximate
    // display curve, not SG's electrophoresis model. Do not clamp small products
    // onto the 100 bp position (which hid 62/34 bp as one band).
    let distance=(16000.0/n.max(1)as f32).ln()/(16000.0_f32/10.0).ln();
    distance.max(0.0).powf(percent.max(0.1).sqrt())
}
impl GelView {
    pub fn parse(bytes:&[u8])->Result<Self,String>{
        let file=hx_formats::DnaFile::parse(bytes).map_err(|e|e.to_string())?;
        let xml=std::str::from_utf8(&file.block(33).ok_or("Missing gel summary")?.payload).map_err(|e|e.to_string())?;
        let inputs=tags(xml,"InputFile");
        let mut markers=Vec::new();
        // Embedded definitions win over catalog entries with the same identity.
        for database in [xml,include_str!("gel_markers.xml")]{for supplier in tags(database,"Supplier"){
            let from=database.find(supplier).unwrap();let end=database[from..].find("</Supplier>").map(|i|from+i).unwrap_or(database.len());let body=&database[from..end];
            for ladder in tags(body,"Ladder"){let from=body.find(ladder).unwrap();let end=body[from..].find("</Ladder>").map(|i|from+i).unwrap_or(body.len());let supplier=attr(supplier,"name");let name=attr(ladder,"name");if markers.iter().any(|m:&Marker|m.supplier==supplier&&m.name==name){continue;}markers.push(Marker{supplier,name,sizes:tags(&body[from..end],"Fragment").iter().filter_map(|t|attr(t,"length").parse().ok()).collect(),definition:body[from..(end+9).min(body.len())].to_owned()});}
        }}
        let mut sources=Vec::new();
        let mut lanes=Vec::new();
        for tag in tags(xml,"LadderLane"){
            let number=attr(tag,"number").parse().map_err(|_|"Invalid marker lane")?;
            let marker=markers.iter().position(|m|m.name==attr(tag,"ladder")&&m.supplier==attr(tag,"supplier")).ok_or("Unknown lane marker")?;
            let mut lane=empty_lane(number);lane.name=markers[marker].name.clone();lane.sizes=markers[marker].sizes.clone();lane.marker=Some(marker);lanes.push(lane);
        }
        for tag in tags(xml,"TemplateLane"){
            let id=attr(tag,"segId");let number=attr(tag,"number").parse().map_err(|_|"Invalid lane number")?;
            let name=inputs.iter().find(|t|attr(t,"segId")==id).map(|t|attr(t,"name")).unwrap_or_else(||format!("Lane {number}"));
            let b=file.blocks.iter().find(|b|b.block_type==23&&b.payload.len()>=4&&u32::from_be_bytes(b.payload[..4].try_into().unwrap()).to_string()==id).ok_or("Missing embedded lane source")?;
            let native=b.payload[4..].to_vec();let source=hx_formats::DnaFile::parse(&native).map_err(|e|e.to_string())?;let rec=source.sequence().ok_or("Lane source has no DNA")?;
            let start=xml.find(tag).ok_or("Missing lane")?;let end=xml[start..].find("</TemplateLane>").map(|i|start+i).unwrap_or(xml.len());let body=&xml[start..end];
            let mut treatment=String::new();let mut supercoiled=false;
            let sizes=if let Some(pcr)=tags(body,"AmplificationInputs").first(){
                treatment=format!("{} + {}",attr(pcr,"primer1"),attr(pcr,"primer2"));
                let s=attr(pcr,"twoPrimerSelectionRange");let pair=s.split('-').map(str::parse::<usize>).collect::<Result<Vec<_>,_>>().map_err(|_|"Invalid PCR range")?;
                if pair.len()!=2||pair[0]==0||pair.iter().any(|n|*n>rec.len()){return Err("Invalid PCR lane coordinates".into());}
                vec![if pair[1]>=pair[0]{pair[1]-pair[0]+1}else{rec.len()-pair[0]+pair[1]+1}]
            }else{
                let mut cuts=Vec::new();let mut names=Vec::new();
                if let Some(digest)=tags(body,"InputDigestions").first(){for key in ["name1","name2","name3","name4"]{let name=attr(digest,key);if name.is_empty(){continue;}let e=hx_enzymes::by_name(&name).ok_or_else(||format!("Unknown enzyme {name}"))?;cuts.extend(hx_algo::cloning::restriction_cuts(&rec.sequence,rec.flags.topology().is_circular(),e).iter().map(|c|c.top));names.push(name);}}
                treatment=names.join(" + ");cuts.sort_unstable();cuts.dedup();
                if cuts.is_empty(){supercoiled=names.is_empty()&&rec.flags.topology().is_circular();vec![rec.len()]}
                else{let mut sizes:Vec<_>=cuts.windows(2).map(|w|(w[1]-w[0])as usize).collect();if rec.flags.topology().is_circular(){sizes.push((rec.len()as i64-cuts.last().unwrap()+cuts[0])as usize);}else{sizes.push(cuts[0]as usize);sizes.push((rec.len()as i64-cuts.last().unwrap())as usize);}sizes.retain(|s|*s>0);sizes.sort_unstable_by(|a,b|b.cmp(a));sizes}
            };
            let enzymes=std::array::from_fn(|i|tags(body,"InputDigestions").first().map(|t|attr(t,&format!("name{}",i+1))).unwrap_or_default());
            let pcr=tags(body,"AmplificationInputs").first().map(|s|s.to_string());
            sources.push(Source{name:name.clone(),bytes:native.clone()});
            lanes.push(Lane{number,name,sizes,treatment,supercoiled,bytes:Some(native),enzymes,marker:None,pcr});
            let lane=lanes.last_mut().unwrap();lane.sizes=lane_fragments(lane)?.iter().map(|f|f.length).collect();
        }
        let max_lanes=tags(xml,"Lanes").first().and_then(|t|attr(t,"maxLanes").parse().ok()).unwrap_or(10).max(lanes.iter().map(|l|l.number+1).max().unwrap_or(1));
        let percent=attr(xml,"percentAgarose").parse().unwrap_or(1.0);
        lanes.sort_by_key(|l|l.number);
        for number in 0..max_lanes{if !lanes.iter().any(|l|l.number==number){lanes.push(empty_lane(number));}}
        lanes.sort_by_key(|l|l.number);
        Ok(Self{key:bytes.to_vec(),lanes,selected:if max_lanes>1{1}else{0},band:None,preview:None,max_lanes,percent,error:String::new(),markers,sources,available:Vec::new(),changed:false,browse_path:String::new(),show_browse:false})
    }
    pub fn select(&mut self,ctx:&egui::Context,index:usize){
        if index>=self.lanes.len(){return;}
        self.selected=index;self.band=None;self.preview=None;
        if let Some(bytes)=self.lanes[index].bytes.as_ref(){match Document::from_bytes(bytes.clone()){
            Ok(mut doc)=>{doc.name_hint=Some(self.lanes[index].name.clone());let mut app=HelixApp::preview(ctx,doc);app.show_side_tools=true;self.preview=Some(Box::new(app));},Err(e)=>self.error=e,
        }}
    }
    pub fn select_fragment(&mut self,ctx:&egui::Context,index:usize,band:Option<usize>){
        if self.selected!=index||self.preview.is_none(){self.select(ctx,index);}
        self.band=band;
        if let Some(preview)=self.preview.as_mut(){preview.selection=None;preview.digest_selection=None;preview.cursor_bp=None;
            if let Some(fragment)=band.and_then(|j|lane_fragments(&self.lanes[index]).ok()?.get(j).cloned()){
                let n=preview.doc.as_ref().map(|d|d.len()).unwrap_or(0);
                preview.selection=Some(if fragment.length==n{(0,n)}else{(fragment.start,fragment.end)});
                preview.digest_selection=digest_selection(&self.lanes[index],&fragment);
                preview.status=format!("Selected: {}",fragment.description());
            }
        }
    }
    pub fn resize(&mut self,count:usize)->Result<(),String>{
        let count=count.clamp(1,40);
        if self.lanes.iter().skip(count).any(|l|l.bytes.is_some()||l.marker.is_some()){return Err("Clear populated trailing lanes before reducing the lane count.".into());}
        self.lanes.truncate(count);while self.lanes.len()<count{self.lanes.push(empty_lane(self.lanes.len()));}self.max_lanes=count;if self.selected>=count{self.selected=count-1;self.preview=None;self.band=None;}self.changed=true;Ok(())
    }
    pub fn recalculate(&mut self,index:usize)->Result<(),String>{
        let lane=self.lanes.get_mut(index).ok_or("Invalid lane")?;
        if lane.marker.is_some(){return Ok(());}let Some(bytes)=lane.bytes.as_ref()else{lane.sizes.clear();return Ok(())};
        let file=hx_formats::DnaFile::parse(bytes).map_err(|e|e.to_string())?;let rec=file.sequence().ok_or("Missing DNA")?;
        let mut candidate=lane.clone();candidate.pcr=None;let fragments=lane_fragments(&candidate)?;
        lane.pcr=None;lane.treatment=lane.enzymes.iter().filter(|s|!s.is_empty()).cloned().collect::<Vec<_>>().join(" + ");lane.supercoiled=rec.flags.topology().is_circular()&&fragments.len()==1&&fragments[0].left=="Start (1)";
        lane.sizes=fragments.iter().map(|f|f.length).collect();self.changed=true;Ok(())
    }
    pub fn serialize(&self)->Result<Vec<u8>,String>{
        use hx_formats::util::xml_escape as esc;
        let mut file=hx_formats::DnaFile::parse(&self.key).map_err(|e|e.to_string())?;
        let old=std::str::from_utf8(&file.block(33).ok_or("Missing gel summary")?.payload).map_err(|e|e.to_string())?;
        let mut ladder=old.find("<LadderDatabase>").and_then(|s|old.find("</LadderDatabase>").map(|e|old[s..e].to_owned())).unwrap_or_else(||"<LadderDatabase>".into());
        let mut added=std::collections::HashSet::new();
        for i in self.lanes.iter().filter_map(|l|l.marker){if !added.insert(i){continue;}let m=&self.markers[i];
            let present=tags(&ladder,"Supplier").iter().any(|s|{if attr(s,"name")!=m.supplier{return false;}let start=ladder.find(s).unwrap();let end=ladder[start..].find("</Supplier>").map(|e|start+e).unwrap_or(ladder.len());tags(&ladder[start..end],"Ladder").iter().any(|l|attr(l,"name")==m.name)});
            if !present{ladder.push_str(&format!("<Supplier name=\"{}\">{}</Supplier>",esc(&m.supplier),m.definition));}
        }ladder.push_str("</LadderDatabase>");
        let time=attr(old,"electrophoresisTime");let mut xml=format!("<GelSummary electrophoresisTime=\"{}\" percentAgarose=\"{}\"><InputFiles>",esc(&time),self.percent);
        let mut embedded=Vec::new();
        for lane in &self.lanes{if let Some(bytes)=&lane.bytes{let id=lane.number+100;let f=hx_formats::DnaFile::parse(bytes).map_err(|e|e.to_string())?;let rec=f.sequence().ok_or("Missing DNA")?;xml.push_str(&format!("<InputFile name=\"{}\" sequenceLength=\"{}\" topology=\"{}\" segId=\"{id}\"/>",esc(&lane.name),rec.len(),if rec.flags.topology().is_circular(){"circular"}else{"linear"}));let mut payload=(id as u32).to_be_bytes().to_vec();payload.extend(bytes);embedded.push(hx_formats::Block::new(23,payload));}}
        xml.push_str(&format!("</InputFiles>{ladder}<Lanes maxLanes=\"{}\">",self.max_lanes));
        for lane in &self.lanes{if let Some(i)=lane.marker{let m=&self.markers[i];xml.push_str(&format!("<LadderLane number=\"{}\" supplier=\"{}\" ladder=\"{}\"/>",lane.number,esc(&m.supplier),esc(&m.name)));}else if lane.bytes.is_some(){xml.push_str(&format!("<TemplateLane number=\"{}\" segId=\"{}\">",lane.number,lane.number+100));if let Some(pcr)=&lane.pcr{xml.push_str(pcr);}else{xml.push_str("<InputDigestions");for(i,name)in lane.enzymes.iter().enumerate(){if !name.is_empty(){xml.push_str(&format!(" name{}=\"{}\"",i+1,esc(name)));}}xml.push_str("/>");}xml.push_str("</TemplateLane>");}}
        xml.push_str("</Lanes></GelSummary>");file.blocks.retain(|b|b.block_type!=23);let block=file.blocks.iter_mut().find(|b|b.block_type==33).unwrap();block.payload=xml.into_bytes();file.blocks.extend(embedded);Ok(file.to_bytes())
    }
    pub fn show(&mut self,ui:&mut egui::Ui){
        let context=ui.ctx().clone();let tr=|s|crate::i18n::text(&context,s);
        if self.preview.is_none()&&self.lanes.get(self.selected).is_some_and(|l|l.bytes.is_some()){self.select_fragment(ui.ctx(),self.selected,self.band);}
        let all=ui.available_rect_before_wrap();let side=(all.width()*0.22).clamp(260.0,340.0);let mut choice=None;
        let left=Rect::from_min_max(all.min,Pos2::new(all.left()+side,all.bottom()));let right=Rect::from_min_max(Pos2::new(left.right()+12.0,all.top()),all.max);
        ui.scope_builder(egui::UiBuilder::new().max_rect(left),|ui|{
            ui.spacing_mut().item_spacing.y=1.0;
            ui.spacing_mut().interact_size.y=18.0;
            ui.spacing_mut().button_padding=Vec2::new(3.0,1.0);
            ui.spacing_mut().scroll=egui::style::ScrollStyle::solid();
            ui.spacing_mut().scroll.bar_width=9.0;
            ui.visuals_mut().widgets.inactive.bg_fill=Color32::from_gray(210);
            ui.visuals_mut().widgets.hovered.bg_fill=Color32::from_gray(185);
            ui.visuals_mut().widgets.active.bg_fill=Color32::from_gray(160);
            let (r,response)=ui.allocate_exact_size(Vec2::new(side,300.0),egui::Sense::click());let p=ui.painter();let gel=Rect::from_min_max(r.min+Vec2::new(43.0,25.0),r.max-Vec2::new(10.0,6.0));
            p.rect_filled(gel,0.0,Color32::from_rgb(46,46,46));let lane_width=gel.width()/self.max_lanes.max(1)as f32;
            let band_y=|n:usize|gel.top()+25.0+approximate_band_fraction(n,self.percent)*(gel.height()-42.0);
            p.text(Pos2::new(r.left()+17.0,r.top()+12.0),egui::Align2::CENTER_CENTER,"kb",egui::FontId::proportional(10.0),crate::theme::TEXT);
            for n in [10000,8000,6000,5000,4000,3000,2000,1500,1000,500]{let y=band_y(n);p.text(Pos2::new(gel.left()-6.0,y),egui::Align2::RIGHT_CENTER,format!("{:.1}",n as f32/1000.0),egui::FontId::proportional(9.0),crate::theme::TEXT);}
            for number in 0..self.max_lanes{let x=gel.left()+(number as f32+0.5)*lane_width;p.text(Pos2::new(x,r.top()+12.0),egui::Align2::CENTER_CENTER,if self.lanes[number].marker.is_some(){"MW".into()}else{number.to_string()},egui::FontId::proportional(10.0),if self.lanes.get(self.selected).is_some_and(|l|l.number==number){crate::theme::ACCENT}else{crate::theme::TEXT});p.rect_stroke(Rect::from_center_size(Pos2::new(x,gel.top()+6.0),Vec2::new(lane_width*0.75,4.0)),0.0,Stroke::new(0.7,Color32::GRAY),egui::StrokeKind::Inside);if response.clicked()&&response.interact_pointer_pos().is_some_and(|pos|(pos.x-x).abs()<lane_width*0.5&&pos.y<gel.top()+15.0){choice=Some((number,None));}}
            let band_painter=p.with_clip_rect(gel);
            for(i,lane)in self.lanes.iter().enumerate(){let x=gel.left()+(lane.number as f32+0.5)*lane_width;for(j,&size)in lane.sizes.iter().enumerate(){let y=band_y(size);let band=Rect::from_center_size(Pos2::new(x,y),Vec2::new(lane_width*0.72,if lane.marker.is_some(){1.2}else{2.0}));let selected=i==self.selected&&self.band==Some(j);band_painter.rect_filled(band,0.0,if selected{Color32::from_rgb(75,170,250)}else if lane.marker.is_some(){Color32::GRAY}else{Color32::WHITE});if response.clicked()&&response.interact_pointer_pos().is_some_and(|pos|gel.contains(pos)&&band.expand(3.0).contains(pos)){choice=Some((i,Some(j)));}}}
            ui.small(crate::i18n::text(ui.ctx(),"Approximate migration · uncut circular lanes retain topology labels"));
            ui.horizontal(|ui|{let mut count=self.max_lanes;if ui.add(egui::DragValue::new(&mut count).range(1..=40)).changed(){if let Err(e)=self.resize(count){self.error=e;}}ui.label(tr("lanes"));if ui.add(egui::DragValue::new(&mut self.percent).range(0.1..=5.0).speed(0.1)).changed(){self.changed=true;}ui.label(tr("% agarose"));});ui.separator();
            let tree_height=(left.bottom()-ui.cursor().top()-36.0).max(100.0);
            egui::ScrollArea::vertical().id_salt("gel_lane_tree").max_height(tree_height).min_scrolled_height(tree_height).show(ui,|ui|{
                popup_style(ui);ui.spacing_mut().interact_size.y=18.0;
                for(i,lane)in self.lanes.iter().enumerate(){
                    let id=ui.id().with(("lane-expanded",i));let mut expanded=ui.ctx().data_mut(|d|d.get_temp::<bool>(id)).unwrap_or(lane.marker.is_none());
                    ui.horizontal(|ui|{ui.spacing_mut().item_spacing.x=2.0;
                        if lane.bytes.is_some()||lane.marker.is_some(){let (r,response)=ui.allocate_exact_size(Vec2::new(14.0,18.0),egui::Sense::click());let c=r.center();let points=if expanded{vec![c+Vec2::new(-3.5,-2.0),c+Vec2::new(3.5,-2.0),c+Vec2::new(0.0,2.5)]}else{vec![c+Vec2::new(-2.0,-3.5),c+Vec2::new(-2.0,3.5),c+Vec2::new(2.5,0.0)]};ui.painter().add(egui::Shape::convex_polygon(points,Color32::GRAY,Stroke::NONE));if response.clicked(){expanded=!expanded;ui.ctx().data_mut(|d|d.insert_temp(id,expanded));}}else{ui.add_space(16.0);}
                        let label=if lane.marker.is_some(){format!("MW: {}",lane.name)}else if lane.bytes.is_none(){format!("{}:",lane.number)}else{format!("{}: {}",lane.number,lane.name)};
                        if menu_row(ui,&label,self.selected==i&&self.band.is_none()).clicked(){choice=Some((i,None));}
                    });
                    if expanded{let fragments=lane_fragments(lane).unwrap_or_default();
                        ui.scope(|ui|{ui.spacing_mut().indent=28.0;ui.indent(("lane-details",i),|ui|{
                            if !lane.treatment.is_empty(){ui.label(egui::RichText::new(&lane.treatment).size(12.0));}
                            for(j,size)in lane.sizes.iter().enumerate(){let label=format!("{}.  {} {}{}",j+1,size,tr("bp"),if lane.supercoiled{tr("  (supercoiled)")}else{""});
                                let response=menu_row(ui,&label,self.selected==i&&self.band==Some(j));if response.clicked(){choice=Some((i,Some(j)));}
                                if let Some(fragment)=fragments.get(j){response.on_hover_text(fragment.description());}
                            }
                        });});
                    }
                }
            });
            ui.add_space((left.bottom()-ui.cursor().top()-30.0).max(0.0));ui.separator();ui.horizontal(|ui|{
                if ui.add_enabled(self.lanes[self.selected].bytes.is_some()||self.lanes[self.selected].marker.is_some(),egui::Button::new(tr("Copy"))).clicked(){
                    let lane=&self.lanes[self.selected];let text=self.band.and_then(|j|lane_fragments(lane).ok()?.get(j).map(Fragment::description)).unwrap_or_else(||format!("{}\n{}\n{}",lane.name,lane.treatment,lane.sizes.iter().map(|n|format!("{n} bp")).collect::<Vec<_>>().join("\n")));ui.ctx().copy_text(text);
                }
                if ui.button(tr("Clear")).clicked(){self.lanes[self.selected]=empty_lane(self.selected);self.preview=None;self.band=None;self.changed=true;}
                ui.add_space(12.0);ui.label(tr("Move:"));
                for(up,label)in [(true,"↑"),(false,"↓")]{if ui.add_enabled(if up{self.selected>0}else{self.selected+1<self.lanes.len()},egui::Button::new(label)).clicked(){let next=if up{self.selected-1}else{self.selected+1};self.lanes.swap(self.selected,next);for(i,l)in self.lanes.iter_mut().enumerate(){l.number=i;}self.selected=next;self.changed=true;}}
            });
        });
        if let Some((i,band))=choice.take(){self.select_fragment(ui.ctx(),i,band);}
        ui.painter().line_segment([Pos2::new(left.right()+5.0,all.top()),Pos2::new(left.right()+5.0,all.bottom())],Stroke::new(1.0,crate::theme::BORDER));
        ui.scope_builder(egui::UiBuilder::new().max_rect(right),|ui|{
            ui.spacing_mut().interact_size.y=22.0;
            ui.spacing_mut().button_padding=Vec2::new(5.0,2.0);
            let mut replacement=None;let selected=self.selected;
            ui.horizontal(|ui|{ui.label(crate::ui_format!(&context,"Lane {}:",selected));egui::ComboBox::from_id_salt("gel_lane_content").selected_text(if self.lanes[selected].bytes.is_none()&&self.lanes[selected].marker.is_none(){tr("Empty")}else{&self.lanes[selected].name}).width(230.0).show_ui(ui,|ui|{
                popup_style(ui);
                if menu_row(ui,tr("Molecular Weight Marker"),self.lanes[selected].marker.is_some()).clicked(){if let Some(m)=self.markers.first(){let mut l=empty_lane(selected);l.marker=Some(0);l.name=m.name.clone();l.sizes=m.sizes.clone();replacement=Some(l);}ui.close();}
                ui.separator();ui.strong(tr("Embedded sequences"));
                for source in &self.sources{if menu_row(ui,&source.name,self.lanes[selected].bytes.as_ref()==Some(&source.bytes)).clicked(){let mut l=empty_lane(selected);l.name=source.name.clone();l.bytes=Some(source.bytes.clone());replacement=Some(l);ui.close();}}
                ui.separator();ui.strong(tr("Open / recent DNA files"));
                for source in &self.available{if menu_row(ui,&source.name,self.lanes[selected].bytes.as_ref()==Some(&source.bytes)).clicked(){let mut l=empty_lane(selected);l.name=source.name.clone();l.bytes=Some(source.bytes.clone());replacement=Some(l);ui.close();}}
                ui.separator();if menu_row(ui,tr("Browse DNA File…"),false).clicked(){self.show_browse=true;ui.close();}if menu_row(ui,tr("Clear"),false).clicked(){replacement=Some(empty_lane(selected));ui.close();}
            });if self.lanes[selected].bytes.is_some(){ui.label(tr(if self.lanes[selected].pcr.is_some(){"Stored PCR product"}else{"cut with"}));}});
            if let Some(lane)=replacement{self.lanes[selected]=lane;self.changed=true;if self.lanes[selected].bytes.is_some(){if let Err(e)=self.recalculate(selected){self.error=e;}}self.select(ui.ctx(),selected);}
            if let Some(marker)=self.lanes[selected].marker{
                let mut chosen=marker;egui::ComboBox::from_id_salt("gel_marker_choice").selected_text(&self.markers[marker].name).width(310.0).height(560.0).show_ui(ui,|ui|{popup_style(ui);let mut indices=(0..self.markers.len()).collect::<Vec<_>>();indices.sort_by_key(|&i|(&self.markers[i].supplier,&self.markers[i].name));let mut supplier=String::new();for i in indices{let m=&self.markers[i];if supplier!=m.supplier{ui.separator();ui.strong(&m.supplier);supplier=m.supplier.clone();}if menu_row(ui,&m.name,chosen==i).clicked(){chosen=i;ui.close();}}});
                if chosen!=marker{self.lanes[selected].marker=Some(chosen);self.lanes[selected].name=self.markers[chosen].name.clone();self.lanes[selected].sizes=self.markers[chosen].sizes.clone();self.changed=true;}
            }else if self.lanes[selected].bytes.is_some(){
                let mut names=self.lanes[selected].enzymes.clone();let mut apply_all=false;
                ui.horizontal(|ui|{for(i,name)in names.iter_mut().enumerate(){if i>0{ui.label("+");}egui::ComboBox::from_id_salt(("gel_digest_enzyme",i)).selected_text(if name.is_empty(){""}else{name.as_str()}).width(95.0).height(340.0).show_ui(ui,|ui|{
                    popup_style(ui);
                    if menu_row(ui,tr("None"),name.is_empty()).clicked(){name.clear();ui.close();}for e in hx_enzymes::ENZYMES{if menu_row(ui,e.name,*name==e.name).clicked(){*name=e.name.to_owned();ui.close();}}
                });}apply_all=ui.button(tr("Apply Enzymes to All Lanes")).clicked();});
                if names!=self.lanes[selected].enzymes||apply_all{let targets=if apply_all{(0..self.lanes.len()).filter(|&i|self.lanes[i].bytes.is_some()).collect::<Vec<_>>()}else{vec![selected]};let old=self.lanes.clone();let changed=self.changed;let mut failure=None;for i in targets{self.lanes[i].enzymes=names.clone();if let Err(e)=self.recalculate(i){failure=Some(e);break;}}if let Some(e)=failure{self.lanes=old;self.changed=changed;self.error=e;}else{self.select_fragment(ui.ctx(),selected,None);}}
            }
            if self.show_browse{ui.horizontal(|ui|{ui.label(tr("DNA file:"));ui.text_edit_singleline(&mut self.browse_path);if ui.button(tr("Load")).clicked(){match Document::from_path(std::path::Path::new(&self.browse_path)){Ok(doc)if doc.is_dna()&&doc.native_class().is_none()=>{match doc.session.save_bytes_with_history(&doc.map_name(),&doc.map_name(),None){Ok((bytes,_))=>{let name=doc.map_name();self.sources.push(Source{name:name.clone(),bytes:bytes.clone()});let mut l=empty_lane(selected);l.name=name;l.bytes=Some(bytes);self.lanes[selected]=l;self.changed=true;if let Err(e)=self.recalculate(selected){self.error=e;}self.select(ui.ctx(),selected);self.show_browse=false;},Err(e)=>self.error=e}},Ok(_)=>self.error="Choose a DNA sequence file".into(),Err(e)=>self.error=e}}if ui.button(tr("Cancel")).clicked(){self.show_browse=false;}});}
            if !self.error.is_empty(){ui.colored_label(Color32::RED,tr(&self.error));}
            ui.small(crate::i18n::text(ui.ctx(),"Edits are saved to the native gel; migration is approximate."));ui.separator();
            let summary=self.band.and_then(|j|lane_fragments(&self.lanes[selected]).ok()?.get(j).map(Fragment::description));
            if let Some(preview)=self.preview.as_mut(){let size=ui.available_size()-Vec2::new(0.0,24.0);if preview.needs_fit{preview.last_size=Some(size);preview.fit();preview.needs_fit=false;}ui.allocate_ui(size,|ui|{ui.push_id("gel_embedded_source",|ui|preview.doc_canvas(ui,CanvasMode::Cloning));});if let Some(summary)=summary{ui.separator();ui.label(crate::ui_format!(&context,"Selected: {}",summary));}}else{ui.label(tr("Select a sample lane to inspect its embedded DNA."));}
        });
        if let Some((i,band))=choice{self.select_fragment(ui.ctx(),i,band);}
    }
}

#[cfg(test)]mod tests{
    /// Entirely synthetic sequences; no installed application/sample dependency.
    fn synthetic_gel()->Vec<u8>{
        let mut file=hx_formats::DnaFile::parse(&hx_formats::write::NewFile::dna(Vec::new()).to_bytes()).unwrap();file.header=hx_formats::dna::Header::new(8);
        let mut xml=String::from("<GelSummary electrophoresisTime=\"00060000\" percentAgarose=\"1\"><InputFiles>");
        for(number,len)in [(1,4505),(2,2686),(3,2746),(4,3387),(5,2174)]{
            let mut seq=vec![b'A';len];if number==3{seq[99..105].copy_from_slice(b"GAATTC");seq[150..156].copy_from_slice(b"AAGCTT");}
            if number==5{for at in [99,1477,2077,2181usize].into_iter().filter(|at|*at+6<=len){seq[at..at+6].copy_from_slice(b"CAGCTG");}}
            let native=hx_formats::write::NewFile::dna(seq).with_flags(hx_core::MoleculeFlags::ECOLI_PLASMID_CIRCULAR).to_bytes();let mut payload=(number as u32).to_be_bytes().to_vec();payload.extend(native);file.blocks.push(hx_formats::Block::new(23,payload));xml.push_str(&format!("<InputFile name=\"Synthetic {number}\" segId=\"{number}\"/>"));
        }
        xml.push_str("</InputFiles><LadderDatabase><Supplier name=\"Synthetic\"><Ladder name=\"Synthetic ladder\"><Fragment length=\"1000\" mass=\"100\"/></Ladder></Supplier></LadderDatabase><Lanes maxLanes=\"10\"><LadderLane number=\"0\" supplier=\"Synthetic\" ladder=\"Synthetic ladder\"/>");
        for number in 1..=5{xml.push_str(&format!("<TemplateLane number=\"{number}\" segId=\"{number}\">"));match number{3=>xml.push_str("<InputDigestions name1=\"EcoRI\" name2=\"HindIII\"/>"),4=>xml.push_str("<AmplificationInputs primer1=\"Fw\" primer2=\"Rv\" twoPrimerSelectionRange=\"101-904\"/>"),_=>xml.push_str("<InputDigestions name1=\"\"/>")};xml.push_str("</TemplateLane>");}xml.push_str("</Lanes></GelSummary>");file.blocks.push(hx_formats::Block::new(33,xml.into_bytes()));file.to_bytes()
    }
    #[test]fn small_fragments_must_not_collapse_onto_one_bottom_band(){
        for percent in [1.0,3.0]{let y=[123,62,34].map(|n|super::approximate_band_fraction(n,percent));assert!(y[0]<y[1]&&y[1]<y[2],"distinct fragment lengths require distinct migration positions: {y:?}");}
        let mut seq=vec![b'A';219];seq[86..90].copy_from_slice(b"GCGC");seq[209..213].copy_from_slice(b"GCGC");seq[26..33].copy_from_slice(b"TCCAGGA");let mut native=hx_formats::write::NewFile::dna(seq);native.flags=hx_core::MoleculeFlags::ECOLI_PLASMID_CIRCULAR;
        let mut lane=super::empty_lane(1);lane.bytes=Some(native.to_bytes());lane.enzymes=["AspLEI".into(),"PfoI".into(),String::new(),String::new()];assert_eq!(super::lane_fragments(&lane).unwrap().iter().map(|f|f.length).collect::<Vec<_>>(),vec![123,62,34]);
    }
    #[test]fn embedded_tabs_can_return_to_map_after_every_view(){
        use eframe::egui;
        let original=synthetic_gel();
        let ctx=egui::Context::default();let mut gel=crate::app::HelixApp::preview(&ctx,crate::doc::Document::from_bytes(original.clone()).unwrap());gel.tab=crate::app::Tab::NativeOverview;crate::theme::apply(&ctx);
        let mut time=0.0;
        let mut frame=|gel:&mut crate::app::HelixApp,events:Vec<egui::Event>|{time+=0.02;let mut out=ctx.run_ui(egui::RawInput{screen_rect:Some(egui::Rect::from_min_size(egui::Pos2::ZERO,egui::vec2(1800.0,1100.0))),time:Some(time),events,..Default::default()},|ui|gel.doc_canvas(ui,crate::app::CanvasMode::Main));out.textures_delta.clear();out};
        let mut output=frame(&mut gel,vec![]);for _ in 0..3{output=frame(&mut gel,vec![]);}
        for label in ["Sequence","Map","Enzymes","Map","Primers","Map"]{
            let pos=output.shapes.iter().find_map(|s|match &s.shape{egui::Shape::Text(t)if t.galley.job.text==label=>Some(t.pos+t.galley.size()*0.5),_=>None}).unwrap_or_else(||panic!("missing tab {label}"));
            frame(&mut gel,vec![egui::Event::PointerMoved(pos)]);
            frame(&mut gel,vec![egui::Event::PointerButton{pos,button:egui::PointerButton::Primary,pressed:true,modifiers:Default::default()}]);
            output=frame(&mut gel,vec![egui::Event::PointerButton{pos,button:egui::PointerButton::Primary,pressed:false,modifiers:Default::default()}]);
            assert_eq!(gel.gel_view.as_ref().unwrap().preview.as_ref().unwrap().tab.label(),label,"clicking {label} must switch the embedded view");
            output=frame(&mut gel,vec![]);
        }
        assert_eq!(gel.gel_view.as_ref().unwrap().key,original);assert!(!gel.doc.as_ref().unwrap().is_dirty());
    }
    #[test]fn fragment_selection_tracks_cut_coordinates_and_wraps_without_editing_source(){
        let mut seq=vec![b'A';219];seq[88..94].copy_from_slice(b"GAATTC");seq[211..217].copy_from_slice(b"GAATTC");
        let mut native=hx_formats::write::NewFile::dna(seq);native.flags=hx_core::MoleculeFlags::ECOLI_PLASMID_CIRCULAR;let bytes=native.to_bytes();let mut lane=super::empty_lane(1);lane.bytes=Some(bytes.clone());lane.enzymes[0]="EcoRI".into();
        let fragments=super::lane_fragments(&lane).unwrap();assert_eq!(fragments.iter().map(|f|f.length).collect::<Vec<_>>(),vec![123,96]);assert_eq!((fragments[1].start,fragments[1].end),(212,89));assert_eq!(fragments[1].description(),"96 bp = EcoRI (212) — EcoRI (89)");
        let ends=super::digest_selection(&lane,&fragments[1]).unwrap();assert_eq!((ends.top_start,ends.top_end,ends.bottom_start,ends.bottom_end),(212,308,216,312));
        let original=synthetic_gel();let mut gel=super::GelView::parse(&original).unwrap();gel.lanes[1]=lane;gel.recalculate(1).unwrap();gel.changed=false;
        let ctx=eframe::egui::Context::default();gel.select_fragment(&ctx,1,Some(1));assert_eq!(gel.preview.as_ref().unwrap().selection,Some((212,89)));assert_eq!(gel.preview.as_ref().unwrap().doc.as_ref().unwrap().session.bytes(),bytes);assert_eq!(gel.key,original);assert!(!gel.changed);
        assert_eq!(gel.preview.as_ref().unwrap().digest_selection,Some(ends));
        gel.select_fragment(&ctx,1,Some(0));assert_eq!(gel.preview.as_ref().unwrap().selection,Some((89,212)));gel.select_fragment(&ctx,1,None);assert!(gel.preview.as_ref().unwrap().selection.is_none());
        gel.lanes[1].enzymes=["EcoRI".into(),"Nt.BspQI".into(),String::new(),String::new()];let sizes=gel.lanes[1].sizes.clone();assert!(gel.recalculate(1).is_err());assert_eq!(gel.lanes[1].sizes,sizes);assert!(!gel.changed);
    }
    #[test]fn linear_and_pcr_fragments_keep_actual_endpoints(){
        let mut seq=vec![b'A';120];seq[19..25].copy_from_slice(b"GAATTC");let mut native=hx_formats::write::NewFile::dna(seq);native.flags=hx_core::MoleculeFlags::ECOLI_PLASMID_LINEAR;
        let mut lane=super::empty_lane(1);lane.bytes=Some(native.to_bytes());lane.enzymes[0]="EcoRI".into();
        let f=super::lane_fragments(&lane).unwrap();assert_eq!(f.iter().map(|f|f.length).collect::<Vec<_>>(),vec![100,20]);assert_eq!(f[0].description(),"100 bp = EcoRI (20) — End (120)");assert_eq!(f[1].left,"Start (1)");
        lane.pcr=Some("<AmplificationInputs primer1=\"Fw\" primer2=\"Rv\" twoPrimerSelectionRange=\"11-90\"/>".into());let f=super::lane_fragments(&lane).unwrap();assert_eq!((f[0].start,f[0].end,f[0].length),(10,90,80));assert_eq!(f[0].left,"Fw (11)");
        assert!(super::digest_selection(&lane,&f[0]).is_none());
    }
    #[test]fn duplex_end_geometry_distinguishes_five_prime_three_prime_and_blunt(){
        for(name,site,offset)in [("EcoRI","GAATTC",4),("PstI","CTGCAG",-4),("EcoRV","GATATC",0)]{
            let mut seq=vec![b'A';100];seq[20..26].copy_from_slice(site.as_bytes());let mut native=hx_formats::write::NewFile::dna(seq);native.flags=hx_core::MoleculeFlags::ECOLI_PLASMID_LINEAR;
            let mut lane=super::empty_lane(1);lane.bytes=Some(native.to_bytes());lane.enzymes[0]=name.into();let fragments=super::lane_fragments(&lane).unwrap();
            let right=fragments.iter().find(|f|f.start>0).unwrap();let d=super::digest_selection(&lane,right).unwrap();assert_eq!(d.bottom_start-d.top_start,offset,"{name}");assert_eq!(d.bottom_end,d.top_end);
        }
    }
    #[test]fn catalog_markers_save_with_band_facts_and_native_history(){
        let original=synthetic_gel();let mut g=super::GelView::parse(&original).unwrap();
        assert!(g.markers.len()>=530);assert!(g.markers.iter().any(|m|m.supplier=="Jena Bioscience"));
        let i=g.markers.iter().position(|m|m.supplier=="Accuris"&&m.name=="SmartCheck™ 1 kb DNA Ladder").unwrap();
        assert_eq!(g.markers[i].sizes,vec![10000,8000,6000,5000,4000,3000,2500,2000,1500,1000,750,500,250]);
        assert!(g.markers[i].definition.contains("mass=\"100\""));
        let mut doc=crate::doc::Document::from_bytes(original.clone()).unwrap();
        let mut versions=vec![original.clone()];
        g.lanes[7]=g.lanes[0].clone();g.lanes[7].number=7;g.lanes[7].marker=Some(i);g.lanes[7].name=g.markers[i].name.clone();g.lanes[7].sizes=g.markers[i].sizes.clone();
        versions.push(g.serialize().unwrap());g.percent=1.7;versions.push(g.serialize().unwrap());g.resize(12).unwrap();versions.push(g.serialize().unwrap());
        g.lanes[8]=g.lanes[3].clone();g.lanes[8].number=8;g.lanes[8].enzymes=Default::default();g.recalculate(8).unwrap();versions.push(g.serialize().unwrap());
        for bytes in versions.iter().skip(1){doc.edit_native(bytes.clone()).unwrap();assert!(doc.can_undo());assert!(doc.is_dirty());}
        for bytes in versions.iter().rev().skip(1){doc.undo().unwrap();assert_eq!(&doc.native_only.as_ref().unwrap().bytes,bytes);}
        assert!(!doc.is_dirty());assert!(!doc.can_undo());assert!(doc.can_redo());
        for bytes in versions.iter().skip(1){doc.redo().unwrap();assert_eq!(&doc.native_only.as_ref().unwrap().bytes,bytes);}
        let saved=std::path::PathBuf::from(format!("/tmp/opencode/ag-gel-history-{}.gel",std::process::id()));doc.save_to(&saved).unwrap();assert!(!doc.is_dirty());
        doc.undo().unwrap();assert!(doc.is_dirty());doc.redo().unwrap();assert!(!doc.is_dirty());
        let restored=super::GelView::parse(&std::fs::read(&saved).unwrap()).unwrap();let m=&restored.markers[restored.lanes[7].marker.unwrap()];assert_eq!(m.sizes,g.markers[i].sizes);assert_eq!(m.definition,g.markers[i].definition);assert_eq!(restored.lanes[4].sizes,vec![804]);
        doc.undo().unwrap();doc.edit_native(versions[1].clone()).unwrap();assert!(!doc.can_redo());let before=doc.native_only.as_ref().unwrap().bytes.clone();doc.edit_native(before.clone()).unwrap();assert_eq!(doc.native_only.as_ref().unwrap().bytes,before);
        std::fs::remove_file(saved).unwrap();
    }
    #[test]fn embedded_source_switching_does_not_mutate_gel(){let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile/Gel Files/Agarose Gel.gel");if !path.exists(){return;}let bytes=std::fs::read(path).unwrap();let mut g=super::GelView::parse(&bytes).unwrap();let ctx=eframe::egui::Context::default();for i in [1,3,4,5]{g.select(&ctx,i);assert_eq!(g.preview.as_ref().unwrap().doc.as_ref().unwrap().len(),match i{1=>4505,3=>2746,4=>3387,_=>2174});assert_eq!(g.key,bytes);} }
    #[test]fn example_lane_definitions_and_fragments_match_stored_design(){let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile/Gel Files/Agarose Gel.gel");if !path.exists(){return;}let b=std::fs::read(path).unwrap();let g=super::GelView::parse(&b).unwrap();assert_eq!(g.max_lanes,10);assert_eq!(g.lanes.len(),10);assert_eq!(g.lanes[1].sizes,vec![4505]);assert!(g.lanes[1].supercoiled);assert_eq!(g.lanes[3].sizes,vec![2695,51]);assert_eq!(g.lanes[4].sizes,vec![804]);assert_eq!(g.lanes[5].sizes,vec![1378,600,104,92]);}
    #[test]fn edited_lanes_marker_and_enzymes_roundtrip(){
        let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile/Gel Files/Agarose Gel.gel");if !path.exists(){return;}let bytes=std::fs::read(path).unwrap();let mut g=super::GelView::parse(&bytes).unwrap();g.resize(12).unwrap();g.percent=1.5;
        let mut lane=g.lanes[3].clone();lane.number=10;lane.enzymes=["EcoRI".into(),"HindIII".into(),String::new(),String::new()];g.lanes[10]=lane;g.recalculate(10).unwrap();assert_eq!(g.lanes[10].sizes,vec![2695,51]);
        let mut marker=g.lanes[0].clone();marker.number=11;g.lanes[11]=marker;assert!(g.resize(10).is_err());
        let edited=g.serialize().unwrap();let restored=super::GelView::parse(&edited).unwrap();assert_eq!(restored.max_lanes,12);assert_eq!(restored.percent,1.5);assert_eq!(restored.lanes[10].enzymes[0],"EcoRI");assert_eq!(restored.lanes[10].sizes,vec![2695,51]);assert_eq!(restored.lanes[11].marker,Some(0));assert_eq!(restored.lanes[4].sizes,vec![804]);
        let mut doc=crate::doc::Document::from_bytes(edited.clone()).unwrap();doc.native_only.as_mut().unwrap().dirty=true;assert!(doc.is_dirty());let p=std::path::PathBuf::from("/tmp/opencode").join(format!("ag-gel-save-{}.gel",std::process::id()));assert!(!p.exists());doc.save_to(&p).unwrap();assert!(!doc.is_dirty());assert_eq!(std::fs::read(&p).unwrap(),edited);std::fs::remove_file(p).unwrap();
    }
}
