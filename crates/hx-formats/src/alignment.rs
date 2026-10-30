//! Native block-24 metadata and block-25 stored multiple alignment.
use crate::{DnaFile,block::parse_blocks,json};
#[derive(Clone,Debug)]pub struct Row{pub id:u32,pub name:String,pub sequence:Vec<u8>}
#[derive(Clone,Debug)]pub struct Alignment{pub protein:bool,pub rows:Vec<Row>,pub method:String}
pub fn parse(file:&DnaFile)->Result<Alignment,String>{
    let meta=json::parse(&file.block(24).ok_or("Missing alignment metadata")?.payload).map_err(|e|e.to_string())?;
    let protein=meta.get("seqType").and_then(|s|s.as_str())==Some("prot");
    let names=meta.get("sequences").and_then(|s|s.as_array()).ok_or("Missing sequence names")?;
    let method=meta.get("alignments").and_then(|s|s.as_array()).and_then(|s|s.first()).and_then(|s|s.get("name")).and_then(|s|s.as_str()).unwrap_or("Stored alignment").to_owned();
    let payload=&file.block(25).ok_or("No stored alignment result")?.payload;
    if payload.len()<4{return Err("Truncated alignment result".into());}
    let blocks=parse_blocks(&payload[4..]).map_err(|e|e.to_string())?;
    let mut rows=Vec::new();
    for b in blocks.iter().filter(|b|b.block_type==16){
        if b.payload.len()<4{return Err("Truncated alignment row".into());}
        let id=u32::from_be_bytes(b.payload[..4].try_into().unwrap());
        let nested=parse_blocks(&b.payload[4..]).map_err(|e|e.to_string())?;
        let seq=nested.iter().find(|b|b.block_type==if protein{21}else{0}).ok_or("Missing aligned sequence")?;
        let sequence=seq.payload.get(1..).ok_or("Truncated aligned sequence")?.to_vec();
        let name=names.iter().find(|n|n.get("ID").and_then(|s|s.as_usize())==Some(id as usize)).and_then(|n|n.get("name")).and_then(|s|s.as_str()).unwrap_or("Unnamed").to_owned();
        rows.push(Row{id,name,sequence});
    }
    if rows.is_empty(){return Err("Empty alignment".into());}
    let width=rows[0].sequence.len();if rows.iter().any(|r|r.sequence.len()!=width){return Err("Unequal stored alignment row lengths".into());}
    Ok(Alignment{protein,rows,method})
}
#[cfg(test)]mod tests{#[test]fn stored_samples_have_aligned_rows(){let root=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/TemplateFile");for(path,n)in [("DNA Alignment/Ape COX1 Alignment.dnaln",5),("Protein Alignment/Fluorescent Protein Alignment.praln",4)]{let p=root.join(path);if !p.exists(){continue;}let bytes=std::fs::read(p).unwrap();let f=crate::DnaFile::parse(&bytes).unwrap();let a=super::parse(&f).unwrap();assert_eq!(a.rows.len(),n);assert!(a.rows[0].sequence.len()>200);assert_eq!(f.to_bytes(),bytes);}}}
