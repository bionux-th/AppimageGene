//! Native product output shared by the Product preview and Clone action.
//! Binding sites are recomputed on the final product, never copied with source
//! coordinates. This uses HX's exact 3′ suffix matcher, not SG's mismatch/Tm oracle.
use super::*;
use hx_formats::primer_edit::{PrimerEdit,PrimerSiteSpec,PrimerSpec};

fn bind(spec:&mut PrimerSpec,product:&Product) {
    let top=product.ends.top_range(product.seq.len());
    let hits=hx_algo::primer::anneal_bindings(&product.seq[top.clone()],spec.sequence.as_bytes(),product.circular,10);
    // PCR tails may be removed by digestion. Keep the longest surviving exact
    // suffix rather than publishing every incidental short match as a main site.
    let longest=hits.iter().map(|h|h.len).max().unwrap_or(0);
    spec.sites=hits.into_iter().filter(|h|h.len==longest).map(|h|{
        let annealed=spec.sequence[spec.sequence.len()-h.len..].to_string();
        let end=if h.end==0 {top.len()}else if h.end>top.len() {h.end%top.len()}else{h.end};
        PrimerSiteSpec{range:hx_formats::annot::Range::new(top.start+h.start,top.start+end),strand:h.strand,
            tm:hx_algo::primer::tm(annealed.as_bytes()),annealed}
    }).collect();
}

pub(super) fn product_bytes(dlg:&CloningDlg,product:&Product,name:&str)->Vec<u8> {
    let mut specs:Vec<PrimerSpec>=Vec::new();
    let mut add=|mut spec:PrimerSpec,designed:bool| {
        if spec.sequence.is_empty() || !spec.sequence.is_ascii() {return;}
        if specs.iter().any(|p|p.name==spec.name && p.sequence.eq_ignore_ascii_case(&spec.sequence)) {return;}
        bind(&mut spec,product);
        // Source primers without a surviving exact site are not retained;
        // designed PCR primers remain listed even if subsequent digestion has
        // removed the binding suffix. No binding coordinates are invented.
        if spec.sites.is_empty() && !designed {return;}
        let base=spec.name.clone();let mut suffix=2;
        while specs.iter().any(|p|p.name==spec.name) {spec.name=format!("{base} ({suffix})");suffix+=1;}
        specs.push(spec);
    };
    let mut picks=Vec::new();
    if dlg.method!=Method::Ligation {picks.push(&dlg.vector);}
    picks.extend(&dlg.fragments);
    if matches!(dlg.method,Method::GatewayBpLr(_)) {picks.push(&dlg.destination);}
    picks.extend(&dlg.gateway_donors);
    for pick in &picks {
        if let Some(source)=dlg.pick_source(pick) {
            for primer in &source.primers {
                let mut spec=PrimerSpec::new(&primer.name,&primer.sequence);
                spec.description=primer.description.clone();spec.date_added=primer.date_added.clone();
                spec.phosphorylated=primer.phosphorylated;
                add(spec,false);
            }
        }
    }
    let mut number=1;
    for pick in std::iter::once(&dlg.vector).chain(dlg.fragments.iter()) {
        let vector=std::ptr::eq(pick,&dlg.vector);
        let uses_pcr=if golden_gate::active(dlg) {pick.gg_mode==GoldenGateMode::PcrDigest}
            else if dlg.overlap_method() {if vector {pick.lin==Linearize::Pcr}else{pick.use_pcr}}
            else {!vector && (dlg.method.is_ta_topo() || dlg.method.is_bp()) && pick.use_pcr};
        if !uses_pcr || dlg.pick_source(pick).is_none() {continue;}
        for (sequence,phosphorylated) in [(&pick.primer1,pick.phos1),(&pick.primer2,pick.phos2)] {
            let mut spec=PrimerSpec::new(format!("Primer {number}"),sequence.trim().to_ascii_uppercase());
            spec.phosphorylated=phosphorylated;add(spec,true);number+=1;
        }
    }
    let bytes=cloning::product_bytes(product,name);
    if specs.is_empty() {return bytes;}
    let mut file=hx_formats::DnaFile::parse(&bytes).expect("valid cloning product");
    let edits=specs.into_iter().map(|s|PrimerEdit::Add(Box::new(s))).collect::<Vec<_>>();
    if let Some(block)=file.blocks.iter_mut().find(|b|b.block_type==5) {
        block.payload=hx_formats::primer_edit::apply_primer_edits(&block.payload,&edits);
    }
    file.to_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_primers_are_rebound_on_product_and_removed_sites_are_not_copied() {
        let seq=b"TTGACCGTATGCTAGCCTGAGTACCGATGCTAGTGACTGACCTAGTCAGTCA".to_vec();
        let primer=hx_formats::annot::PrimerAnnotation{name:"retained".into(),sequence:String::from_utf8(seq[..22].to_vec()).unwrap(),
            description:"source primer".into(),date_added:"2026-10-28T00:00:00Z".into(),phosphorylated:true,
            sites:vec![hx_formats::annot::PrimerSite{range:hx_formats::annot::Range::new(900,922),strand:0,annealed_bases:22,melting_temperature:Some(60.0)}]};
        let mut deleted=primer.clone();deleted.name="discarded".into();deleted.sequence="AAAAAAAAAAAAAAAAAAAA".into();
        let source=Source::new("input",seq.clone(),false).with_primers(vec![primer,deleted]);
        let dlg=CloningDlg::new(Method::Ligation,vec![DocSnap::new("input",source)],Some(0));
        let product=Product{seq:[b"GGG".to_vec(),seq].concat(),features:vec![],spans:vec![],junctions:vec![],circular:false,ends:Default::default()};
        let bytes=product_bytes(&dlg,&product,"rebound");
        let doc=crate::doc::Document::from_bytes(bytes.clone()).unwrap();
        assert_eq!(doc.primers().len(),1);
        let p=&doc.primers()[0];assert_eq!(p.name,"retained");assert!(p.phosphorylated);
        assert_eq!(p.sites[0].range,hx_formats::annot::Range::new(3,25));
        assert_eq!(doc.session.save_bytes().unwrap(),bytes);
    }

    #[test]
    fn binding_at_origin_uses_a_valid_closed_range_without_bottom_end_markers() {
        let sequence="ACGTCAGTGCTAGTCCATGCTA";
        let product=Product{seq:[sequence.as_bytes().to_vec(),b"TTGACCGATCGT".to_vec()].concat(),features:vec![],spans:vec![],junctions:vec![],circular:true,ends:Default::default()};
        let mut primer=PrimerSpec::new("origin",sequence);bind(&mut primer,&product);
        assert!(primer.sites.iter().any(|s|s.range==hx_formats::annot::Range::new(0,sequence.len())));
        let end_sequence=String::from_utf8(product.seq[product.seq.len()-22..].to_vec()).unwrap();
        let mut end_primer=PrimerSpec::new("ends at origin",end_sequence);bind(&mut end_primer,&product);
        assert!(end_primer.sites.iter().any(|s|s.range.end==product.seq.len()));
        let mut linear=product;linear.circular=false;linear.seq.insert(0,b'T');linear.ends.upstream_stickiness=-1;
        bind(&mut primer,&linear);
        assert!(primer.sites.iter().all(|s|s.range.start>=1));
        assert_eq!(primer.sites[0].range,hx_formats::annot::Range::new(1,sequence.len()+1));
    }
}
