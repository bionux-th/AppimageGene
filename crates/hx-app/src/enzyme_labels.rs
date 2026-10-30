//! Deterministic collision-free label lanes above each sequence row.
#[derive(Clone,Debug)]
pub struct Label {pub site:usize,pub x:f32,pub left:f32,pub width:f32,pub lane:usize}
pub fn lane_offsets(labels:&[Label],sites:&[(String,usize,bool)])->Vec<usize>{let n=labels.iter().map(|l|l.lane+1).max().unwrap_or(0);let mut heights=vec![1;n];for l in labels{heights[l.lane]=heights[l.lane].max(sites[l.site].0.lines().count());}let mut offsets=vec![0];for h in heights{offsets.push(offsets.last().unwrap()+h);}offsets}

#[cfg(test)]mod side_label_tests{
    #[test]fn text_is_beside_not_centred_on_leader(){
        let(labels,_)=super::layout(&[(0,100.0,35.0),(1,5.0,30.0)],0.0,200.0);
        for label in labels{assert!(label.left>=label.x+3.0||label.left+label.width<=label.x-3.0);}
    }
}
#[derive(Clone,Debug,PartialEq,Eq)]
pub struct SiteGeometry{pub start:usize,pub pattern:String,pub cuts:Vec<(Option<i64>,Option<i64>)>}
pub fn geometry(index:usize,start:usize,seq:&[u8],circular:bool)->SiteGeometry{
    let e=&hx_enzymes::catalog::CATALOG[index];
    for pattern in e.site.split_whitespace(){
        let matches=|p:&[u8]|p.iter().enumerate().all(|(i,&b)|seq.get(if circular{(start+i)%seq.len().max(1)}else{start+i}).is_some_and(|&s|hx_enzymes::base_matches(b,s)));
        let forward=matches(pattern.as_bytes());let rc=hx_algo::reverse_complement(pattern.as_bytes());
        if forward||matches(&rc){let n=pattern.len()as i64;let cuts=e.cuts.iter().map(|(a,b)|if forward{(a.map(|v|start as i64+v as i64),b.map(|v|start as i64+v as i64))}else{(b.map(|v|start as i64+n-v as i64),a.map(|v|start as i64+n-v as i64))}).collect();return SiteGeometry{start,pattern:if forward{pattern.to_owned()}else{String::from_utf8_lossy(&rc).into_owned()},cuts};}
    }SiteGeometry{start,pattern:e.site.split_whitespace().next().unwrap_or("").into(),cuts:Vec::new()}
}
pub fn layout(items:&[(usize,f32,f32)],left:f32,right:f32)->(Vec<Label>,usize){
    let mut ordered=items.to_vec();ordered.sort_by(|a,b|a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut occupied:Vec<Vec<(f32,f32)>>=Vec::new();let mut labels=Vec::new();
    for(site,x,width)in ordered{let width=width.min((right-left).max(1.0));
        // Prefer text immediately left of its leader; near the left edge use
        // the right side. Never centre the word across the vertical stem.
        let a=if x-width-4.0>=left{x-width-4.0}else{(x+4.0).min((right-width).max(left))};let b=a+width;
        let lane=occupied.iter().position(|v|v.iter().all(|&(c,d)|b+5.0<=c||a>=d+5.0)).unwrap_or(occupied.len());if lane==occupied.len(){occupied.push(Vec::new());}occupied[lane].push((a,b));labels.push(Label{site,x,left:a,width,lane});
    }(labels,occupied.len())
}
#[cfg(test)]mod tests{#[test]fn dense_labels_are_retained_and_do_not_overlap(){let input=(0..150).map(|i|(i,50.0+(i%10)as f32,40.0)).collect::<Vec<_>>();let(out,n)=super::layout(&input,0.0,200.0);assert_eq!(out.len(),150);assert!(n>1);for(a,i)in out.iter().enumerate(){for j in &out[a+1..]{if i.lane==j.lane{assert!(i.left+i.width<=j.left||j.left+j.width<=i.left);}}}}}

#[cfg(test)]mod geometry_tests{#[test]fn nonpalindromic_reverse_cut_is_mirrored(){let i=hx_enzymes::catalog::CATALOG.iter().position(|e|e.name=="BsaI").unwrap();let g=super::geometry(i,0,b"GAGACC",false);assert_eq!(g.pattern,"GAGACC");assert_eq!(g.cuts,vec![(Some(-5),Some(-1))]);}#[test]fn identical_cut_geometry_groups_isoschizomers(){let find=|n|hx_enzymes::catalog::CATALOG.iter().position(|e|e.name==n).unwrap();assert_eq!(super::geometry(find("XhoI"),0,b"CTCGAG",false),super::geometry(find("PaeR7I"),0,b"CTCGAG",false));}}
