//! AppimageGene identity, rendered without external image/font dependencies.
use eframe::egui::{self,Color32,Pos2,Rect,Stroke};

pub fn app_icon()->egui::IconData{
    let size=128;let mut rgba=vec![0;size*size*4];
    for y in 0..size {for x in 0..size {
        let dx=(x as f32-64.0)/1.12;let dy=(y as f32-64.0)/1.06;
        let u=64.0+dx*0.848+dy*0.530;let v=64.0-dx*0.530+dy*0.848;
        let inside=(53.0..75.0).contains(&u)&&(6.0..14.0).contains(&v)
            ||(61.0..67.0).contains(&u)&&(14.0..27.0).contains(&v)
            ||(46.0..82.0).contains(&u)&&(26.0..34.0).contains(&v)
            ||(47.0..81.0).contains(&u)&&(34.0..72.0).contains(&v)
            ||(72.0..87.0).contains(&v)&&(u-64.0).abs()<14.0-(v-72.0)*0.65
            ||(60.0..68.0).contains(&u)&&(87.0..106.0).contains(&v)
            ||(106.0..128.0).contains(&v)&&(u-64.0).abs()<4.0-(v-106.0)*3.0/22.0
            ||(78.0..88.0).contains(&u)&&(37.0..58.0).contains(&v);
        if inside{let i=(y*size+x)*4;rgba[i..i+4].copy_from_slice(&[40,186,255,255]);}
    }}egui::IconData{rgba,width:size as u32,height:size as u32}
}

#[cfg(test)]mod tests{
    #[test]fn app_icon_is_nonempty_rgba_with_transparent_background(){let icon=super::app_icon();assert_eq!(icon.rgba.len(),(icon.width*icon.height*4)as usize);assert!(icon.rgba.chunks_exact(4).any(|p|p[3]==255));assert!(icon.rgba.chunks_exact(4).any(|p|p[3]==0));}
    #[test]fn supported_formats_have_distinct_icons_and_unknown_stays_blank(){use super::*;for(n,k)in [("readme.TXT",FileKind::Text),("paper.pdf",FileKind::Pdf),("gene.fa",FileKind::Fasta),("gene.frna",FileKind::Fasta),("vector.gbk",FileKind::Genbank),("protein.gp",FileKind::Genbank),("x.dnaln",FileKind::Alignment),("x.gel",FileKind::Gel),("x.ftrs",FileKind::Features),("x.zip",FileKind::Unknown)]{assert_eq!(file_kind(n),k);}}
}

pub fn document_icon(p:&egui::Painter,r:Rect,kind:hx_core::MoleculeKind,white:bool){
    let color=if white{Color32::WHITE}else{let(a,b,c)=kind.icon_rgb();Color32::from_rgb(a,b,c)};
    let pt=|x:f32,y:f32|Pos2::new(r.left()+x*r.width(),r.top()+y*r.height());
    p.add(egui::Shape::closed_line(vec![pt(0.15,0.05),pt(0.68,0.05),pt(0.88,0.25),pt(0.88,0.95),pt(0.15,0.95)],Stroke::new(1.1,color)));
    p.add(egui::Shape::line(vec![pt(0.68,0.05),pt(0.68,0.25),pt(0.88,0.25)],Stroke::new(1.0,color)));
    match kind{
        hx_core::MoleculeKind::Dna=>{
            for phase in [0.0,std::f32::consts::PI]{p.add(egui::Shape::line((0..=20).map(|i|{let t=i as f32/20.0;pt(0.51+(t*8.0+phase).sin()*0.18,0.30+t*0.52)}).collect(),Stroke::new(1.0,color)));}
        },
        hx_core::MoleculeKind::Rna=>{p.add(egui::Shape::line((0..=20).map(|i|{let t=i as f32/20.0;pt(0.51+(t*10.0).sin()*0.18,0.30+t*0.52)}).collect(),Stroke::new(1.4,color)));},
        hx_core::MoleculeKind::Protein=>{
            let points=vec![pt(0.35,0.36),pt(0.64,0.34),pt(0.71,0.58),pt(0.38,0.60),pt(0.33,0.81),pt(0.65,0.79)];
            p.add(egui::Shape::line(points.clone(),Stroke::new(1.0,color)));for at in points{p.circle_filled(at,r.width()*0.065,color);}
        }
    }
}

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum FileKind{Dna,Rna,Protein,Text,Pdf,Fasta,Genbank,Alignment,Gel,Features,Sequence,Unknown}
pub fn file_kind(name:&str)->FileKind{
    let ext=std::path::Path::new(name).extension().map(|s|s.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
    match ext.as_str(){
        "dna"=>FileKind::Dna,"rna"=>FileKind::Rna,"prot"=>FileKind::Protein,
        "txt"=>FileKind::Text,"pdf"=>FileKind::Pdf,
        "fasta"|"fa"|"fas"|"fna"|"ffn"|"faa"|"frna"=>FileKind::Fasta,
        "gb"|"gbk"|"genbank"|"gp"=>FileKind::Genbank,
        "dnaln"|"praln"=>FileKind::Alignment,"gel"=>FileKind::Gel,"ftrs"=>FileKind::Features,
        "seq"|"aa"|"prt"|"xdna"|"snapgene"=>FileKind::Sequence,_=>FileKind::Unknown,
    }
}

/// Shared format glyphs for tree rows and document tabs. No font symbols.
pub fn file_icon(p:&egui::Painter,r:Rect,name:&str,white:bool){
    let kind=file_kind(name);
    match kind{FileKind::Dna=>return document_icon(p,r,hx_core::MoleculeKind::Dna,white),FileKind::Rna=>return document_icon(p,r,hx_core::MoleculeKind::Rna,white),FileKind::Protein=>return document_icon(p,r,hx_core::MoleculeKind::Protein,white),_=>{}}
    let color=if white{Color32::WHITE}else{match kind{FileKind::Pdf=>Color32::from_rgb(215,63,72),FileKind::Fasta=>Color32::from_rgb(43,153,112),FileKind::Genbank=>Color32::from_rgb(73,109,192),FileKind::Text=>Color32::from_rgb(110,122,139),FileKind::Alignment=>Color32::from_rgb(119,91,180),FileKind::Gel=>Color32::from_rgb(75,120,155),FileKind::Features=>Color32::from_rgb(201,133,41),_=>Color32::from_rgb(150,158,170)}};
    let pt=|x:f32,y:f32|Pos2::new(r.left()+x*r.width(),r.top()+y*r.height());
    let line=|a:(f32,f32),b:(f32,f32)|{p.line_segment([pt(a.0,a.1),pt(b.0,b.1)],Stroke::new(1.0,color));};
    p.add(egui::Shape::closed_line(vec![pt(0.12,0.04),pt(0.67,0.04),pt(0.90,0.26),pt(0.90,0.96),pt(0.12,0.96)],Stroke::new(1.05,color)));
    line((0.67,0.04),(0.67,0.26));line((0.67,0.26),(0.90,0.26));
    match kind{
        FileKind::Text=>{for y in [0.40,0.54,0.68,0.82]{line((0.28,y),(if y>0.8{0.59}else{0.75},y));}},
        FileKind::Pdf=>{p.add(egui::Shape::line(vec![pt(0.27,0.79),pt(0.56,0.36),pt(0.44,0.72),pt(0.76,0.66),pt(0.35,0.62),pt(0.27,0.79)],Stroke::new(1.2,color)));},
        FileKind::Fasta=>{line((0.27,0.38),(0.41,0.47));line((0.41,0.47),(0.27,0.56));for y in [0.64,0.79]{line((0.27,y),(0.76,y));}},
        FileKind::Genbank=>{for y in [0.4,0.6,0.8]{line((0.26,y),(0.75,y));p.rect_filled(Rect::from_min_max(pt(0.36,y-0.05),pt(0.62,y+0.05)),0.0,color);}},
        FileKind::Alignment=>{for y in [0.41,0.57,0.73]{for x in [0.28,0.44,0.64]{line((x,y),(x+0.09,y));}}},
        FileKind::Gel=>{for x in [0.32,0.53,0.73]{line((x,0.38),(x,0.83));}for (x,y)in [(0.32,0.47),(0.32,0.68),(0.53,0.55),(0.73,0.72)]{line((x-0.07,y),(x+0.07,y));}},
        FileKind::Features=>{p.add(egui::Shape::closed_line(vec![pt(0.27,0.43),pt(0.61,0.43),pt(0.77,0.62),pt(0.61,0.80),pt(0.27,0.80)],Stroke::new(1.1,color)));},
        FileKind::Sequence=>{for y in [0.43,0.63,0.81]{line((0.29,y),(0.73,y));}},_=>{}
    }
}
