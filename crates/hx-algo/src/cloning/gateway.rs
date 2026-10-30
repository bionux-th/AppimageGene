//! Numbered MultiSite Gateway recombination.
//! Site-specificity cores and L/R polarity drive crossover, never generic overlap.
//! Reference: MultiSite Gateway Pro manual, pp. vii, 6–10 and 13–18;
//! local site annotations/sequence from publicly shipped donor vectors.
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind { B, P, L, R }

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub kind: Kind,
    pub number: u8,
    pub start: usize,
    pub end: usize,
    pub direction: i32,
    pub cut: usize,
}

fn name(name: &str) -> Option<(Kind,u8)> {
    let name=name.trim().to_ascii_lowercase();let tail=name.strip_prefix("att")?;
    let kind=match tail.as_bytes().first()? { b'b'=>Kind::B,b'p'=>Kind::P,b'l'=>Kind::L,b'r'=>Kind::R,_=>return None };
    let digits=tail[1..].trim_end_matches('r');let number=digits.parse().ok()?;
    (1..=5).contains(&number).then_some((kind,number))
}

fn core(number:u8)->Option<&'static [u8]> {
    match number {1=>Some(b"TTTGTACAAAAAAG"),2=>Some(b"TTTGTACAAGAAAG"),3=>Some(b"TTTGTATAATAAAG"),
        4=>Some(b"TTTGTATAGAAAAG"),5=>Some(b"TTTGTATACAAAAG"),_=>None}
}

/// Canonical short attB sequence (5′→3′ in the site arrow direction).
/// 1/2 retain this project's established 25-base oracle convention; 3/4/5 use
/// the manufacturer's short functional sequence, without optional GGGG clamp.
pub fn att_b(number:u8)->Option<Vec<u8>> {
    match number {
        1=>Some(ATT_B1.to_vec()),2=>Some(revcomp_str(ATT_B2)),
        3=>Some(b"ACAACTTTGTATAATAAAGTTG".to_vec()),
        4=>Some(b"ACAACTTTGTATAGAAAAGTTG".to_vec()),
        5=>Some(b"ACAACTTTGTATACAAAAGTTG".to_vec()),_=>None,
    }
}

pub fn sites(src:&Source)->Result<Vec<Site>,String> {
    let mut result=Vec::new();
    for feature in &src.features {
        let Some((kind,number))=name(&feature.name) else {continue};
        if feature.segments.len()!=1 {return Err(format!("{} must have one contiguous range for {}.",src.label,feature.name));}
        let r=feature.segments[0].range;
        if r.start>=r.end || r.end>src.seq.len() || !matches!(feature.directionality,1|2) {
            return Err(format!("Invalid {} coordinates or direction in {}.",feature.name,src.label));
        }
        let seq=if feature.directionality==2 {revcomp_str(&src.seq[r.start..r.end])} else {src.seq[r.start..r.end].to_vec()};
        let motif=core(number).unwrap();
        let matches:Vec<_>=seq.windows(motif.len()).enumerate().filter(|(_,s)|s.eq_ignore_ascii_case(motif)).map(|(i,_)|i).collect();
        if matches.len()!=1 {return Err(format!("The sequence of {} in {} does not contain its unique numbered specificity core.",feature.name,src.label));}
        let offset=matches[0]+7;
        let cut=if feature.directionality==2 {r.end-offset} else {r.start+offset};
        result.push(Site{kind,number,start:r.start,end:r.end,direction:feature.directionality,cut});
    }
    // Unannotated PCR products: recognize complete short attB cores, not arbitrary
    // same-name features or the 14-base core alone in a longer attL/P/R site.
    for number in 1..=5 {
        let motif=att_b(number).unwrap();
        for reverse in [false,true] {
            let pattern=if reverse {revcomp_str(&motif)} else {motif.clone()};
            for pos in hx_enzymes::find_site_positions(&src.seq,&pattern,false) {
                let start=pos as usize;let end=start+pattern.len();
                if result.iter().any(|s|start<s.end && end>s.start) {continue;}
                let anchor=find_motif(&motif,core(number).unwrap()).unwrap()+7;
                result.push(Site{kind:Kind::B,number,start,end,direction:if reverse {2}else{1},
                    cut:if reverse {end-anchor}else{start+anchor}});
            }
        }
    }
    result.sort_by_key(|s|s.start);Ok(result)
}

fn pair(src:&Source,kinds:&[Kind])->Result<(Site,Site),String> {
    let sites:Vec<_>=sites(src)?.into_iter().filter(|s|kinds.contains(&s.kind)).collect();
    if sites.len()!=2 || sites[0].number==sites[1].number || sites[0].end>sites[1].start {
        return Err(format!("{} must contain exactly two distinct, non-overlapping compatible att sites.",src.label));
    }
    Ok((sites[0].clone(),sites[1].clone()))
}

fn complementary(a:&Site,b:&Site)->bool {
    a.number==b.number && a.direction==b.direction
        && matches!((a.kind,b.kind),(Kind::L,Kind::R)|(Kind::R,Kind::L))
}

/// Product from an ordered list of 2–4 entry clones and a destination cassette.
/// L/R can occur on either side (MultiSite Pro and R4/R3 layouts).
pub fn lr(inserts:&[&Source],destination:&Source)->Result<Product,String> {
    if !(2..=4).contains(&inserts.len()) {return Err("MultiSite LR requires two to four Entry Clones.".into());}
    if !destination.circular {return Err("MultiSite LR currently requires a circular Destination Vector; a linear source cannot be silently circularized.".into());}
    let (dl,dr)=pair(destination,&[Kind::R])?;
    let mut pairs=Vec::new();
    for source in inserts {pairs.push(pair(source,&[Kind::L,Kind::R])?);}
    if !complementary(&dl,&pairs[0].0) || !complementary(&pairs.last().unwrap().1,&dr) {
        return Err("The outer Entry Clone sites do not match the Destination Vector attR sites and directions.".into());
    }
    let mut used=vec![dl.number];
    for i in 0..pairs.len()-1 {
        if !complementary(&pairs[i].1,&pairs[i+1].0) {
            return Err(format!("Entry Clones {} and {} require matching-number attL/attR sites in the same direction at their junction.",i+1,i+2));
        }
        let n=pairs[i].1.number;
        if used.contains(&n) {return Err("Repeated att-site specificity makes the MultiSite reaction ambiguous.".into());}
        used.push(n);
    }
    if used.contains(&dr.number) {return Err("Repeated outer att-site specificity makes the reaction ambiguous.".into());}
    let mut slices=vec![destination.slice(0,dl.cut,false)];
    for (source,(left,right)) in inserts.iter().zip(&pairs) {slices.push(source.slice(left.cut,right.cut,false));}
    slices.push(destination.slice(dr.cut,destination.seq.len(),false));
    let mut product=assemble(&slices,true,0,0);product.features.retain(|f|!is_att_name(&f.name));
    let mut offset=dl.cut;
    add_junction_feature(&mut product,&dl,&pairs[0].0,offset)?;
    for i in 0..pairs.len()-1 {
        offset+=pairs[i].1.cut-pairs[i].0.cut;
        add_junction_feature(&mut product,&pairs[i].1,&pairs[i+1].0,offset)?;
    }
    offset+=pairs.last().unwrap().1.cut-pairs.last().unwrap().0.cut;
    add_junction_feature(&mut product,&pairs.last().unwrap().1,&dr,offset)?;
    product.junctions.clear();Ok(product)
}

fn add_junction_feature(product:&mut Product,left:&Site,right:&Site,at:usize)->Result<(),String> {
    let a=left.cut-left.start;let b=right.end-right.cut;
    let start=at.checked_sub(a).ok_or("Invalid recombinant att boundary.")?;
    if at+b>product.seq.len() {return Err("Invalid recombinant att boundary.".into());}
    product.features.push(make_site_feature(&format!("attB{}",left.number),start,a+b,left.direction));Ok(())
}

/// General numbered BP: match the two attB endpoints to the donor's attP
/// endpoints. Reverse-specific sites produce attR, not incorrectly named attL.
pub fn bp(insert:&Source,donor:&Source)->Result<Product,String> {
    if !donor.circular {return Err("Numbered BP currently requires a circular Donor Vector; a linear source cannot be silently circularized.".into());}
    let (pl,pr)=pair(donor,&[Kind::P])?;
    let (bl,br)=pair(insert,&[Kind::B])?;
    if pl.number!=bl.number || pr.number!=br.number || pl.direction!=bl.direction || pr.direction!=br.direction {
        return Err("The attB pair and directions must match the selected Donor Vector's attP pair.".into());
    }
    let slices=vec![donor.slice(0,pl.cut,false),insert.slice(bl.cut,br.cut,false),donor.slice(pr.cut,donor.seq.len(),false)];
    let mut product=assemble(&slices,true,0,0);product.features.retain(|f|!is_att_name(&f.name));
    let mut left=pl.start;let mut right=pl.cut+bl.end-bl.cut;
    if pl.direction==1 && left>0 {left-=1;}
    if pl.direction==2 && right<product.seq.len() {right+=1;}
    product.features.push(make_site_feature(&format!("att{}{}",if pl.direction==1 {'L'}else{'R'},pl.number),left,right-left,pl.direction));
    let at=pl.cut+br.cut-bl.cut;
    let mut left=at-(br.cut-br.start);let mut right=at+pr.end-pr.cut;
    if pr.direction==2 && right<product.seq.len() {right+=1;}
    if pr.direction==1 && left>0 {left-=1;}
    product.features.push(make_site_feature(&format!("att{}{}",if pr.direction==2 {'L'}else{'R'},pr.number),left,right-left,pr.direction));
    product.junctions.clear();Ok(product)
}

/// Standard MultiSite Pro donor arrangement (manual figures for 2/3/4 fragments).
pub fn layout(count:usize)->Option<&'static [(u8,i32,u8,i32)]> {
    match count {
        1=>Some(&[(1,1,2,2)]),2=>Some(&[(1,1,5,1),(5,1,2,2)]),
        3=>Some(&[(1,1,4,2),(4,2,3,1),(3,1,2,2)]),
        4=>Some(&[(1,1,5,1),(5,1,4,2),(4,2,3,1),(3,1,2,2)]),_=>None,
    }
}

pub fn donor_matches(source:&Source,expected:(u8,i32,u8,i32))->bool {
    pair(source,&[Kind::P]).is_ok_and(|(a,b)|(a.number,a.direction,b.number,b.direction)==expected)
}

pub fn primer_sites(count:usize,index:usize)->Result<(Vec<u8>,Vec<u8>),String> {
    let &(a,ad,b,bd)=layout(count).and_then(|s|s.get(index)).ok_or("Select one to four Gateway fragments.")?;
    let a=att_b(a).unwrap();let b=att_b(b).unwrap();
    Ok((if ad==2 {revcomp_str(&a)}else{a},if bd==2 {b}else{revcomp_str(&b)}))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn annotated(kind:Kind,n:u8,dir:i32)->Vec<u8> {
        let prefix=match kind {Kind::L=>b"GGGGGGGGGG".as_slice(),Kind::R=>b"ACAACT".as_slice(),_=>b"AAAAAA".as_slice()};
        let suffix=match kind {Kind::L=>b"TTGG".as_slice(),Kind::R=>b"CCCCCCCCCC".as_slice(),_=>b"GGGG".as_slice()};
        let seq=[prefix,core(n).unwrap(),suffix].concat();
        if dir==2 {revcomp_str(&seq)}else{seq}
    }

    fn entry(left:(Kind,u8,i32),right:(Kind,u8,i32),body:&[u8])->Source {
        let a=annotated(left.0,left.1,left.2);let b=annotated(right.0,right.1,right.2);
        let mut source=Source::new("controlled Entry",[b"AAAA".to_vec(),a.clone(),body.to_vec(),b.clone(),b"TTTT".to_vec()].concat(),true);
        let letter=|k|match k {Kind::L=>'L',Kind::R=>'R',Kind::B=>'B',Kind::P=>'P'};
        source.features=vec![make_site_feature(&format!("att{}{}",letter(left.0),left.1),4,a.len(),left.2),
            make_site_feature("payload",4+a.len(),body.len(),1),
            make_site_feature(&format!("att{}{}",letter(right.0),right.1),4+a.len()+body.len(),b.len(),right.2)];
        source
    }

    #[test]
    fn all_pro_layouts_recombine_exactly_and_preserve_payloads_without_installed_resources() {
        for n in 2..=4 {
            let specs=layout(n).unwrap();
            let entries:Vec<_>=specs.iter().enumerate().map(|(i,&(a,ad,b,bd))| {
                entry((if ad==1 {Kind::L}else{Kind::R},a,ad),
                    (if bd==2 {Kind::L}else{Kind::R},b,bd),
                    format!("ATG{}TAA",["ACGT","GTCA","TGAC","CAGT"][i].repeat(15)).as_bytes())
            }).collect();
            let mut destination=entry((Kind::R,1,1),(Kind::R,2,2),&b"C".repeat(80));
            destination.label="destination".into();
            let pairs:Vec<_>=entries.iter().map(|s|pair(s,&[Kind::L,Kind::R]).unwrap()).collect();
            let (dl,dr)=pair(&destination,&[Kind::R]).unwrap();
            // Independently assemble the expected crossover pieces from the
            // documented core offsets, without calling the assembly primitive.
            let mut expected=destination.seq[..dl.cut].to_vec();
            for (s,(a,b)) in entries.iter().zip(&pairs) {expected.extend_from_slice(&s.seq[a.cut..b.cut]);}
            expected.extend_from_slice(&destination.seq[dr.cut..]);
            let product=lr(&entries.iter().collect::<Vec<_>>(),&destination).unwrap();
            assert_eq!(product.seq,expected);assert!(product.circular);
            assert_eq!(product.features.iter().filter(|f|f.name=="payload").count(),n);
            assert_eq!(product.features.iter().filter(|f|f.name.starts_with("attB")).count(),n+1);
            let mut broken=entries.clone();broken[1]=broken[1].flipped();
            assert!(lr(&broken.iter().collect::<Vec<_>>(),&destination).is_err());
            let repeated=vec![&entries[0];n];assert!(lr(&repeated,&destination).is_err());
        }
    }

    #[test]
    fn matching_number_but_wrong_site_family_or_direction_is_rejected() {
        let dest=entry((Kind::R,1,1),(Kind::R,2,2),b"CCCC");
        let a=entry((Kind::L,1,1),(Kind::R,5,1),b"AAAA");
        let b=entry((Kind::R,5,1),(Kind::L,2,2),b"GGGG");
        assert!(lr(&[&a,&b],&dest).unwrap_err().contains("matching-number"));
        let b=entry((Kind::L,5,2),(Kind::L,2,2),b"GGGG");
        assert!(lr(&[&a,&b],&dest).is_err());
    }

    #[test]
    fn r4_r3_destination_uses_a_different_numbered_chain_without_special_casing() {
        let dest=entry((Kind::R,4,1),(Kind::R,3,2),b"CCCC");
        let a=entry((Kind::L,4,1),(Kind::R,1,1),b"AAAACCCC");
        let b=entry((Kind::L,1,1),(Kind::L,2,2),b"GGGGTTTT");
        let c=entry((Kind::R,2,2),(Kind::L,3,2),b"ATGCGCAT");
        let product=lr(&[&a,&b,&c],&dest).unwrap();
        let names=product.features.iter().filter(|f|f.name.starts_with("attB")).map(|f|f.name.as_str()).collect::<Vec<_>>();
        assert_eq!(names,vec!["attB4","attB1","attB2","attB3"]);
    }

    #[test]
    fn unrelated_unicode_feature_names_are_ignored_without_panicking() {
        let mut source=entry((Kind::L,1,1),(Kind::R,5,1),b"AAAA");
        source.features.push(make_site_feature("attβ1",0,3,1));
        assert_eq!(sites(&source).unwrap().len(),2);
    }

    #[test]
    fn linear_destination_is_not_silently_circularized_by_multisite_lr() {
        let mut dest=entry((Kind::R,1,1),(Kind::R,2,2),b"CCCC");dest.circular=false;
        let a=entry((Kind::L,1,1),(Kind::R,5,1),b"AAAA");
        let b=entry((Kind::L,5,1),(Kind::L,2,2),b"GGGG");
        assert!(lr(&[&a,&b],&dest).is_err(),"linear destination topology is not supported by this circular-product workflow");
    }

    #[test]
    fn linear_donor_is_not_silently_circularized_by_numbered_bp() {
        let mut donor=entry((Kind::P,1,1),(Kind::P,5,1),b"CCCC");donor.circular=false;
        let insert=entry((Kind::B,1,1),(Kind::B,5,1),b"AAAA");
        assert!(bp(&insert,&donor).is_err(),"linear donor topology is not supported by this circular-product workflow");
    }
}
