//! Uses the manufacturer's MultiSite Pro layouts and installed public donor
//! sequences, not a claim that a local test constitutes an SG output oracle.
use hx_algo::cloning::{self,gateway,Source};
use hx_formats::{DnaFile,annot::{parse_features,Range}};
use std::path::Path;

fn load(name:&str)->Option<Source> {
    let path=Path::new("/opt/gslbiotech/snapgene/resources/Plasmids").join(name);
    if !path.is_file(){return None;}
    let file=DnaFile::parse(&std::fs::read(path).unwrap()).unwrap();let rec=file.sequence().unwrap();
    Some(Source::new(name,rec.sequence,rec.flags.topology().is_circular())
        .with_features(file.block(10).map(|b|parse_features(&b.payload)).unwrap_or_default()).with_ends(file.end_properties()))
}

fn donor_names(n:usize)->Vec<&'static str> {
    match n {2=>vec!["pDONR221 P1-P5r.dna","pDONR221 P5-P2.dna"],
        3=>vec!["pDONR221 P1-P4.dna","pDONR221 P4r-P3r.dna","pDONR221 P3-P2.dna"],
        _=>vec!["pDONR221 P1-P5r.dna","pDONR221 P5-P4.dna","pDONR221 P4r-P3r.dna","pDONR221 P3-P2.dna"]}
}

fn entry_set(n:usize)->Option<Vec<Source>> {
    donor_names(n).iter().enumerate().map(|(i,name)| {
        let donor=load(&format!("Gateway Donor Vectors/{name}"))?;
        assert!(gateway::donor_matches(&donor,gateway::layout(n).unwrap()[i]));
        let (f,r)=gateway::primer_sites(n,i).unwrap();
        let gene=format!("ATG{}TAA",["ACGTGAC","GTCACGA","TCGAGAC","GACTGTC"][i].repeat(12)).into_bytes();
        let pcr=Source::new(format!("fragment {}",i+1),[f,gene,cloning::Source::new("reverse",r,false).flipped().seq].concat(),false);
        let product=gateway::bp(&pcr,&donor).unwrap();
        let source=Source::new(format!("entry {}",i+1),product.seq,true).with_features(product.features);
        let sites=gateway::sites(&source).unwrap();
        assert_eq!(sites.len(),2);
        Some(source)
    }).collect()
}

#[test]
fn two_three_four_fragment_lr_succeeds_with_real_numbered_donor_sequences() {
    let Some(dest)=load("Gateway Destination Vectors/pDEST17.dna") else {eprintln!("Skipping installed Gateway vectors: resources absent");return;};
    for n in 2..=4 {
        let entries=entry_set(n).unwrap();
        let refs=entries.iter().collect::<Vec<_>>();
        let product=gateway::lr(&refs,&dest).unwrap();
        assert!(product.circular);assert!(product.seq.len()>3000);
        let att:Vec<_>=product.features.iter().filter(|f|f.name.starts_with("attB")).collect();
        assert_eq!(att.len(),n+1);
        let numbers:Vec<_>=att.iter().map(|f|f.name.clone()).collect();
        assert_eq!(numbers,match n {2=>vec!["attB1","attB5","attB2"],3=>vec!["attB1","attB4","attB3","attB2"],_=>vec!["attB1","attB5","attB4","attB3","attB2"]});
        for f in &product.features {for s in &f.segments {assert!(s.range.start<s.range.end && s.range.end<=product.seq.len());}}
        let file=DnaFile::parse(&cloning::product_bytes(&product,"MultiSite")).unwrap();
        assert_eq!(file.sequence().unwrap().sequence,product.seq);
        assert_eq!(parse_features(&file.block(10).unwrap().payload).len(),product.features.len());
        let wrong:Vec<_>=entries.iter().rev().collect();assert!(gateway::lr(&wrong,&dest).is_err());
        let incomplete:Vec<_>=entries.iter().take(n-1).collect();assert!(gateway::lr(&incomplete,&dest).is_err());
        // Reverse the whole coordinate representation, then restore it before
        // use: site identities and crossover mapping must round-trip.
        let restored:Vec<_>=entries.iter().map(|e|e.flipped().flipped()).collect();
        assert_eq!(gateway::lr(&restored.iter().collect::<Vec<_>>(),&dest).unwrap().seq,product.seq);
        let cores=gateway::sites(&Source::new("product",product.seq.clone(),true).with_features(product.features.clone())).unwrap();
        assert_eq!(cores.len(),n+1);
        assert!(cores.iter().all(|s|s.kind==gateway::Kind::B));
        if let Ok(dir)=std::env::var("HELIX_MULTISITE_OUTPUT_DIR") {
            let path=Path::new(&dir).join(format!("{n}-fragments"));std::fs::create_dir_all(&path).unwrap();
            std::fs::write(path.join("product.dna"),file.to_bytes()).unwrap();
            for (i,entry) in entries.iter().enumerate() {
                let p=cloning::assemble(&[entry.slice(0,entry.seq.len(),false)],true,0,0);
                std::fs::write(path.join(format!("entry-{}.dna",i+1)),cloning::product_bytes(&p,"Entry")).unwrap();
            }
        }
    }
}

#[test]
fn malformed_numbered_sites_are_rejected_instead_of_panicking() {
    let mut source=Source::new("bad",b"A".repeat(200),true);
    let file=hx_formats::write::NewFile::dna(b"A".repeat(200)).to_bytes();
    let mut f=DnaFile::parse(&file).unwrap();
    f.blocks.iter_mut().find(|b|b.block_type==10).unwrap().payload=b"<Features><Feature name=\"attL1\" directionality=\"1\"><Segment range=\"1-100\"/></Feature></Features>".to_vec();
    source.features=parse_features(&f.block(10).unwrap().payload);
    assert!(gateway::sites(&source).is_err());
    source.features[0].segments[0].range=Range::new(0,500);assert!(gateway::sites(&source).is_err());
    assert!(gateway::att_b(0).is_none());assert!(gateway::att_b(6).is_none());
}

#[test]
fn generalized_bp_preserves_the_established_single_fragment_oracle_sequence() {
    let Some(donor)=load("Gateway Donor Vectors/pDONR221.dna") else{return;};
    let insert=Source::new("attB",[cloning::ATT_B1.to_vec(),b"ACGT".repeat(30),cloning::ATT_B2.to_vec()].concat(),false);
    let legacy=cloning::gateway_recombine(&[&insert],&donor,"attB Fragment",true).unwrap();
    let numbered=gateway::bp(&insert,&donor).unwrap();
    assert_eq!(numbered.seq,legacy.seq);
    for name in ["attL1","attL2"] {
        let a=numbered.features.iter().find(|f|f.name==name).unwrap();
        let b=legacy.features.iter().find(|f|f.name==name).unwrap();
        assert_eq!(a.segments[0].range,b.segments[0].range);
    }
}
