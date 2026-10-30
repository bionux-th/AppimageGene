//! Exercises the actual dialog build/preview/native-document path, not only a
//! low-level concatenation primitive. Local success is NOT SG equivalence.
use super::*;

fn fasta(path:&str)->Vec<u8>{
    std::fs::read_to_string(path).unwrap().lines().filter(|l|!l.starts_with('>')).flat_map(|l|l.trim().bytes()).collect()
}

fn source(seq: Vec<u8>, circular: bool) -> Source { Source::new("controlled", seq, circular) }
fn doc(seq: Vec<u8>, circular: bool) -> DocSnap { DocSnap::new("controlled", source(seq, circular)) }
fn enzyme(name: &str) -> usize { hx_enzymes::ENZYMES.iter().position(|e| e.name == name).unwrap() }

fn persisted(dlg: &mut CloningDlg) -> Product {
    dlg.preview_key.clear(); dlg.refresh();
    let product = dlg.preview.clone().unwrap_or_else(|| panic!("{:?}: {:?}", dlg.method, dlg.error));
    assert!(dlg.primers_ready(), "{:?} must be actionable, not only a preview", dlg.method);
    assert!(product.circular, "{:?} must construct a plasmid", dlg.method);
    assert!(!product.seq.is_empty());
    let bytes = output::product_bytes(dlg,&product, "controlled-product");
    let loaded = crate::doc::Document::from_bytes(bytes.clone()).unwrap();
    assert!(loaded.circular());
    assert_eq!(loaded.sequence(), product.seq);
    assert_eq!(loaded.session.save_bytes().unwrap(), bytes);
    assert!(product.features.iter().flat_map(|f| &f.segments).all(|s| s.range.end <= product.seq.len()));
    // Render the ready Product page too. This checks application wiring and
    // layout execution, not physical mouse interaction or pixel parity with SG.
    let old_tab=dlg.tab;
    dlg.tab=dlg.tabs().len()-1;
    let ctx=egui::Context::default();
    let input=egui::RawInput{screen_rect:Some(Rect::from_min_size(Pos2::ZERO,Vec2::new(1800.0,1100.0))),..Default::default()};
    let mut frame=ctx.run_ui(input,|ui|{
        let ctx=ui.ctx().clone();
        assert!(draw_dialog(dlg,&ctx,ui).is_none());
    });
    assert!(!frame.shapes.is_empty());frame.textures_delta.clear();
    dlg.tab=old_tab;
    if let Ok(root) = std::env::var("HELIX_CLONING_SMOKE_DIR") {
        let name = format!("{:?}", dlg.method).replace(['(', ')'], "");
        let variant = if (dlg.overlap_method() && dlg.vector.lin == Linearize::Pcr)
            || (golden_gate::active(dlg) && dlg.vector.gg_mode == GoldenGateMode::PcrDigest)
            || ((dlg.method.is_ta_topo() || dlg.method.is_bp()) && dlg.fragments[0].use_pcr) { "pcr" } else { "direct" };
        let dir = Path::new(&root).join(name).join(variant).join(product.seq.len().to_string());
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("product.dna"), &bytes).unwrap();
        let mut details = format!("Method: {:?}\nProduct: {} bp, circular\nEvidence: LOCAL BUILD/SAVE/RELOAD; NOT SG equivalence\n\n", dlg.method, product.display_len());
        for (i, snapshot) in dlg.docs.iter().enumerate() {
            std::fs::write(dir.join(format!("source-{i}.dna")), snapshot.bytes.clone().unwrap_or_else(|| source_bytes(&snapshot.src, &snapshot.name))).unwrap();
        }
        for (label, pick) in std::iter::once(("Vector".into(), &dlg.vector))
            .chain(dlg.fragments.iter().enumerate().map(|(i, p)| (format!("Fragment {}", i + 1), p)))
            .chain(std::iter::once(("Destination".into(), &dlg.destination)))
        {
            details.push_str(&format!("{label}: source={:?}, region={:?}, reverse={}, linearize={:?}, PCR={}, GoldenGate={:?}\nPrimer1={}\nPrimer2={}\nCuts={:?}\n", pick.doc, pick.region, pick.reversed, pick.lin, pick.use_pcr, pick.gg_mode, pick.primer1, pick.primer2, pick.cut_pair));
        }
        details.push_str(&format!("Golden Gate junctions: {:?}\n", dlg.gg_overhangs));
        std::fs::write(dir.join("recipe.txt"), details).unwrap();
    }
    product
}

fn saved(file:&hx_formats::DnaFile)->Vec<u8>{
    file.sequence().unwrap().sequence
}

#[test]
fn practical_saved_snapgene_oracles_rebuild_and_reload_identical_products(){
    // Same synthetic input recipe that was actually run through SG 8.2.3 GUI;
    // compare exact native product sequence, topology and terminal properties.
    let vseq=format!("{}GAATTC{}GGATCC{}","A".repeat(100),"C".repeat(30),"ACCG".repeat(170));
    let iseq=format!("{}GAATTC{}GGATCC{}","T".repeat(20),"ACCA".repeat(40),"T".repeat(20));
    let v=Source::new("synthetic-vector",vseq.into_bytes(),true);
    let i=Source::new("synthetic-insert",iseq.into_bytes(),false);
    let eco=hx_enzymes::by_name("EcoRI").unwrap();let bam=hx_enzymes::by_name("BamHI").unwrap();
    let cuts=|s:&Source,e:&hx_enzymes::Enzyme,circular| cloning::restriction_cuts(&s.seq,circular,e)[0];
    let vc1=cuts(&v,eco,true);let vc2=cuts(&v,bam,true);
    let ic1=cuts(&i,eco,false);let ic2=cuts(&i,bam,false);
    let sticky_oracle=fasta(concat!(env!("CARGO_MANIFEST_DIR"),"/../../crates/hx-algo/tests/fixtures/cloning/restriction-ecori-bamhi.fasta"));
    for (fill,len) in [(false,952),(true,960)]{
        let fragments=vec![
            cloning::digest_fragment(&v,vc2,vc1,fill,false).unwrap(),
            cloning::digest_fragment(&i,ic1,ic2,fill,false).unwrap(),
        ];
        let product=cloning::restriction_product(&fragments).unwrap();
        let expected=if fill {
            [sticky_oracle[..105].to_vec(),b"AATT".to_vec(),sticky_oracle[105..271].to_vec(),b"GATC".to_vec(),sticky_oracle[271..].to_vec()].concat()
        } else {sticky_oracle.clone()};
        assert_eq!(product.seq,expected,"fill={fill}");
        assert!(product.circular);assert_eq!(product.display_len(),len);
        let bytes=cloning::product_bytes(&product,"Cloned");
        let file=hx_formats::DnaFile::parse(&bytes).unwrap();
        assert!(file.sequence().unwrap().flags.topology()==hx_core::Topology::Circular);
        assert_eq!(saved(&file),product.seq);
        assert_eq!(file.end_properties(),product.ends);
    }

    let mut left=Source::new("left",[b"T".to_vec(),b"G".repeat(79),b"AATT".to_vec()].concat(),false).with_ends(hx_formats::dna::EndProperties{
        upstream_stickiness:-1,downstream_stickiness:4,upstream_modification:"Unmodified".into(),downstream_modification:"FivePrimePhosphorylated".into(),original_first_base:None,
    });
    let right=Source::new("right",[b"AATT".to_vec(),b"C".repeat(80),b"AA".to_vec()].concat(),false).with_ends(hx_formats::dna::EndProperties{
        upstream_stickiness:4,downstream_stickiness:2,upstream_modification:"Unmodified".into(),downstream_modification:"Unmodified".into(),original_first_base:None,
    });
    let product=cloning::ligate_reaction(&[left.clone(),right.clone()],false).unwrap();
    let expected=fasta(concat!(env!("CARGO_MANIFEST_DIR"),"/../../crates/hx-algo/tests/fixtures/cloning/ligation-phosphorylated.fasta"));
    assert_eq!(product.seq,expected);assert!(!product.circular);assert_eq!(product.display_len(),163);
    assert_eq!(product.seq.len(),166);
    let bytes=cloning::product_bytes(&product,"Ligation");let file=hx_formats::DnaFile::parse(&bytes).unwrap();
    assert_eq!(saved(&file),expected);assert_eq!(file.end_properties(),product.ends);
    left.ends.downstream_modification="Unmodified".into();
    assert!(cloning::ligate_reaction(&[left,right],false).unwrap_err().contains("at least one DNA end"));
}

#[test]
fn practical_topo_ta_installed_vector_cloning_saves_and_reloads(){
    let root=Path::new("/opt/gslbiotech/snapgene/resources/Plasmids");
    let vp=root.join("TOPO TA Vectors/pCR2.1-TOPO (linearized).dna");
    if !vp.is_file(){eprintln!("Skipping installed TOPO TA practice: vector resource missing");return;}
    let vector=doc_from_path(&vp).unwrap();
    assert_eq!(vector.src.display_len(),3930);
    let mut insert_seq=b"ACGTTGCAAGTCCTGACTGAGTACCATGGCATCGTACGAT".repeat(16);
    insert_seq.extend_from_slice(b"ACGTTGCAAGTCCTGACTGA");
    assert_eq!(insert_seq.len(),660);
    let mut dlg=CloningDlg::new(Method::TopoTa,vec![vector,doc(insert_seq,false)],Some(1));
    dlg.vector.doc=Some(0);dlg.fragments[0].doc=Some(1);dlg.fragments[0].use_pcr=false;dlg.fragments[0].add_a=true;
    let product=persisted(&mut dlg);
    // The specific SG-saved 4581 bp example remains a separate legacy oracle;
    // this practical run deliberately uses a 660 bp insert and different input.
    // Installed vector adds exactly 3930 top-strand bases; insert is 660 bp
    // plus two terminal A bases. This runs the dialog through native save/reload.
    assert_eq!(product.seq.len(),4591);
    let bytes=cloning::product_bytes(&product,"Practical-TOPO-TA");
    let file=hx_formats::DnaFile::parse(&bytes).unwrap();
    assert_eq!(saved(&file),product.seq);assert_eq!(file.sequence().unwrap().flags.topology(),hx_core::Topology::Circular);
    assert_eq!(file.end_properties(),product.ends);
}

#[test]
fn existing_4581bp_topo_ta_sg_oracle_is_available_in_repository_fixtures(){
    let expected=fasta(concat!(env!("CARGO_MANIFEST_DIR"),"/tests/fixtures/cloning/snapgene-cloned-4581.fasta"));
    assert_eq!(expected.len(),4581);
    let xml=include_str!("../../tests/fixtures/cloning/snapgene-cloned-4581-primers.xml");
    let primers=hx_formats::annot::parse_primers(xml.as_bytes());
    assert_eq!(primers.iter().map(|p|p.name.as_str()).collect::<Vec<_>>(),vec!["Primer 1","Primer 2"]);
    assert_eq!(primers[0].sites.iter().filter(|s|s.melting_temperature.unwrap_or(0.0)>=40.0).count(),1);
    assert_eq!(primers[1].sites.iter().filter(|s|s.melting_temperature.unwrap_or(0.0)>=40.0).count(),1);
}

#[test]
fn homology_family_success_matrix_direct_single_multi_and_pcr_design() {
    let a = "ACGTTGCAAGTCCTGACTGA";
    let b = "GTACCATGGCATCGTACGAT";
    let c = "CATGAGTCACGTTCAGTGAC";
    for method in [Method::Gibson, Method::InFusion, Method::NeBuilder,
        Method::GibsonMulti, Method::InFusionMulti, Method::NeBuilderMulti,
        Method::GibsonAssemble, Method::NeBuilderAssemble]
    {
        let multi = method.is_multiple();
        let vector = format!("{b}{}{a}", "G".repeat(60));
        let first = format!("{a}{}{}", "C".repeat(60), if multi { c } else { b });
        let last = format!("{c}{}{b}", "T".repeat(60));
        let mut dlg = CloningDlg::new(method, vec![doc(vector.into_bytes(), false), doc(first.into_bytes(), false), doc(last.into_bytes(), false)], Some(0));
        dlg.vector.doc = Some(0); dlg.vector.lin = Linearize::Direct;
        dlg.fragments[0].doc = Some(1); dlg.fragments[0].use_pcr = false;
        if multi { dlg.fragments[1].doc = Some(2); dlg.fragments[1].use_pcr = false; }
        let p = persisted(&mut dlg);
        assert_eq!(p.seq.len(), if multi { 240 } else { 160 }, "{:?}", method);
        dlg.fragments[0].reversed = true;
        assert!(dlg.build().is_err(), "wrong direction must not silently assemble");
    }
    // SG 默认 Gibson 同源臂是 40 bp；片段本身必须留出非重叠核心，
    // 所以这里用 100 bp 片段（40 bp 的合成片段会被正确地拒绝）。
    for method in [Method::Gibson, Method::InFusion, Method::NeBuilder] {
        let mut dlg = CloningDlg::new(method, vec![
            doc(b"AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTAT".to_vec(), true),
            doc(b"TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA".to_vec(), false)], Some(0));
        dlg.vector.doc = Some(0); dlg.vector.region = Some((0, 100));
        dlg.fragments[0].doc = Some(1); dlg.fragments[0].region = Some((0, 100));
        design_homology_primers(&mut dlg).unwrap();
        let product=persisted(&mut dlg);
        assert_eq!(product.seq.len(), 200, "{:?}", method);
        let bytes=output::product_bytes(&dlg,&product,"PCR cloning");
        let loaded=crate::doc::Document::from_bytes(bytes).unwrap();
        assert_eq!(loaded.primers().len(),4,"{method:?} must save both PCR pairs");
        for (primer,expected) in loaded.primers().iter().zip([
            &dlg.vector.primer1,&dlg.vector.primer2,&dlg.fragments[0].primer1,&dlg.fragments[0].primer2]) {
            assert_eq!(&primer.sequence,expected);
            assert!(!primer.sites.is_empty());
        }
        dlg.preview=Some(product);dlg.tab=dlg.tabs().len()-1;
        let preview=dlg.build_view_doc().unwrap();
        let stable=|doc:&crate::doc::Document|doc.primers().iter().map(|p|
            (p.name.clone(),p.sequence.clone(),p.phosphorylated,p.sites.clone())).collect::<Vec<_>>();
        assert_eq!(stable(&preview),stable(&loaded),"preview and Clone must share primer sequence/sites (timestamps may differ)");
    }
}

#[test]
fn restriction_family_builds_plasmids_from_real_cut_sites() {
    let v = b"AAAAGAATTCCCCCCCCCCCCCCCCCGGATCCGGGG".to_vec();
    let i = b"TTGAATTCAAACCCGGATCCAA".to_vec();
    let mut dlg = CloningDlg::new(Method::Restriction, vec![doc(v, true), doc(i, false)], Some(0));
    dlg.vector.doc = Some(0); dlg.vector.enzyme = Some(enzyme("EcoRI")); dlg.vector.enzyme2 = Some(enzyme("BamHI"));
    let cuts = double_cuts_of(&dlg.docs[0].src, &dlg.vector); dlg.vector.cut_pair = Some((cuts[1], cuts[0]));
    dlg.fragments[0].doc = Some(1); dlg.fragments[0].enzyme = Some(enzyme("EcoRI")); dlg.fragments[0].enzyme2 = Some(enzyme("BamHI"));
    assert_eq!(persisted(&mut dlg).seq, b"AAAAGAATTCAAACCCGGATCCGGGG");
    for tab in [0,1,2] {
        dlg.tab=tab;
        let ctx=egui::Context::default();crate::theme::apply(&ctx);
        let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(Rect::from_min_size(Pos2::ZERO,Vec2::new(1800.0,1100.0))),..Default::default()},|ui|{
            draw_dialog(&mut dlg,&ctx,ui);
        });
        output.textures_delta.clear();
        let bottom_text=output.shapes.iter().filter_map(|s|match &s.shape {
            egui::Shape::Text(t) if t.pos.y>900.0=>Some(t.galley.job.text.as_str()),_=>None,
        }).collect::<Vec<_>>();
        assert!(bottom_text.iter().any(|t|t.contains("Product")),"restriction reaction scheme missing on tab {tab}: {bottom_text:?}");
    }
}

#[test]
fn restriction_multiple_joins_every_fragment_including_the_closing_pair() {
    use hx_formats::dna::EndProperties;
    let part = |label: &str, left: &str, body: &str, right: &str| {
        DocSnap::new(label,Source::new(label,format!("{left}{body}{right}").into_bytes(),false).with_ends(EndProperties {
            upstream_stickiness:4,downstream_stickiness:4,
            upstream_modification:"FivePrimePhosphorylated".into(),downstream_modification:"FivePrimePhosphorylated".into(),
            ..Default::default()
        }))
    };
    let docs=vec![part("vector","AATG","GGGGGGGG","GCTT"),part("insert1","GCTT","CCCCCCCC","CGAC"),part("insert2","CGAC","TTTTTTTT","AATG")];
    let mut dlg=CloningDlg::new(Method::RestrictionMulti,docs,Some(0));
    dlg.vector.doc=Some(0);dlg.fragments[0].doc=Some(1);dlg.fragments[1].doc=Some(2);
    assert_eq!(persisted(&mut dlg).seq,b"AATGGGGGGGGGGCTTCCCCCCCCCGACTTTTTTTT");
    dlg.docs[2].src.seq[0]=b'T';
    assert!(dlg.build().is_err());
}

#[test]
fn golden_gate_builds_single_and_multiple_insert_plasmids() {
    for method in [Method::GoldenGate, Method::GoldenGateMulti] {
        let mut dlg = CloningDlg::new(method, vec![doc(b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".to_vec(), true),
            doc(b"CCCCAAAATTTTGGGGCCCCAAAATTTTGGGG".to_vec(), false)], Some(0));
        dlg.vector.doc = Some(0); dlg.vector.region = Some((0, 32));
        for p in &mut dlg.fragments { p.doc = Some(1); p.region = Some((0, 32)); }
        dlg.gg_overhangs = if method.is_multiple() { vec!["AATG".into(), "GCTT".into(), "CGAC".into()] }
            else { vec!["AATG".into(), "GCTT".into()] };
        golden_gate::design(&mut dlg).unwrap();
        assert_eq!(persisted(&mut dlg).seq.len(), if method.is_multiple() { 108 } else { 72 });
    }
}

#[test]
fn ta_gc_and_all_topo_modes_construct_method_specific_plasmids() {
    use hx_formats::dna::EndProperties;
    for method in [Method::TaCloning, Method::GcCloning, Method::TopoTa, Method::TopoBlunt, Method::TopoDirectional] {
        let vector = match method {
            Method::GcCloning => source(b"GAAAACCCCGGGGTTTTC".to_vec(), false).with_ends(EndProperties {
                upstream_stickiness: -1, downstream_stickiness: -1, ..Default::default()
            }),
            Method::TopoBlunt => source(b"AAAACCCCGGGGTTTT".to_vec(), false),
            Method::TopoDirectional => source(b"AAAACCCCGGGGTTTTCACC".to_vec(), false).with_ends(EndProperties {
                downstream_stickiness: 4, ..Default::default()
            }),
            _ => source(b"AAAAACCCCGGGGTTTTT".to_vec(), false).with_ends(EndProperties {
                upstream_stickiness: -1, downstream_stickiness: -1, ..Default::default()
            }),
        };
        let insert = source(if method == Method::TopoDirectional { b"CACCATGGGGCCCCTTTT".to_vec() } else { b"ATGGGGCCCCTTTT".to_vec() }, false);
        let mut dlg = CloningDlg::new(method, vec![DocSnap::new("vector", vector), DocSnap::new("insert", insert)], Some(1));
        dlg.vector.doc = Some(0); dlg.fragments[0].doc = Some(1); dlg.fragments[0].use_pcr = false;
        dlg.fragments[0].add_a = matches!(method, Method::TaCloning | Method::GcCloning | Method::TopoTa);
        let p = persisted(&mut dlg);
        assert_eq!(p.seq.len(), if method == Method::TopoDirectional { 34 } else if method == Method::TopoBlunt { 30 } else { 32 });
        if method == Method::TopoBlunt {
            dlg.docs[1].src = cloning::add_a_overhangs(&dlg.docs[1].src);
            assert!(dlg.build().unwrap_err().contains("blunt ends"));
        }
    }
}

#[test]
fn every_method_rejects_missing_inputs_instead_of_publishing_a_product() {
    for method in [Method::Restriction,Method::RestrictionMulti,Method::Ligation,
        Method::Gibson,Method::GibsonMulti,Method::GibsonAssemble,
        Method::InFusion,Method::InFusionMulti,Method::NeBuilder,Method::NeBuilderMulti,Method::NeBuilderAssemble,
        Method::GoldenGate,Method::GoldenGateMulti,Method::TaCloning,Method::GcCloning,
        Method::TopoTa,Method::TopoBlunt,Method::TopoDirectional,Method::GatewayBp,
        Method::GatewayLr(1),Method::GatewayBpLr(1),Method::GatewayLr(2),Method::GatewayBpLr(2)]
    {
        let mut dlg=CloningDlg::new(method,Vec::new(),None);dlg.refresh();
        assert!(dlg.preview.is_none(),"{:?}",method);
        assert!(dlg.error.is_some(),"{:?}",method);
    }
}

#[test]
fn method_specific_ends_are_not_interchangeable() {
    use hx_formats::dna::EndProperties;
    let blunt=source(b"ATGGGGCCCCTTTT".to_vec(),false);
    let tvector=source(b"AAAAACCCCGGGGTTTTT".to_vec(),false).with_ends(EndProperties {
        upstream_stickiness:-1,downstream_stickiness:-1,..Default::default()
    });
    let atailed=cloning::add_a_overhangs(&blunt);
    assert!(cloning::end_cloning_product(&tvector,&atailed,Method::TopoTa,false).is_ok());
    assert!(cloning::end_cloning_product(&tvector,&atailed,Method::GcCloning,false).is_err());
    assert!(cloning::end_cloning_product(&tvector,&atailed,Method::TopoBlunt,false).is_err());
    assert!(cloning::end_cloning_product(&tvector,&atailed,Method::TopoDirectional,false).is_err());
    let mut invalid=atailed.clone();invalid.ends.upstream_stickiness=0;
    assert!(cloning::end_cloning_product(&tvector,&invalid,Method::TopoTa,false).is_err());
}

#[test]
fn invalid_multisite_gateway_never_masquerades_as_success() {
    let insert=source([cloning::ATT_B1.to_vec(),b"ACGT".repeat(25),cloning::ATT_B2.to_vec()].concat(),false);
    let vector=source(b"A".repeat(500),true).with_features(vec![site("attR1",30,100,1),site("attR2",300,100,2)]);
    let error=cloning::gateway_recombine(&[&insert,&insert],&vector,"Entry Clone",false).unwrap_err();
    assert!(error.contains("att") || error.contains("specificity"));
}

fn site(name: &str, start: usize, len: usize, dir: i32) -> hx_formats::annot::FeatureAnnotation {
    hx_formats::annot::FeatureAnnotation {
        name: name.into(), kind: "protein_bind".into(), sub_type: String::new(), note: String::new(), visible: true,
        directionality: dir, reading_frame: 0, translate_first_codon_as_met: false, max_run_on: None,
        segments: vec![hx_formats::annot::FeatureSegment { range: hx_formats::annot::Range::new(start, start + len),
            color: "#338899".into(), seg_type: "standard".into(), translated: false }],
    }
}

#[test]
fn gateway_single_bp_lr_and_bp_lr_preview_and_save_products() {
    let p1 = b"AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA";
    let p2 = b"TACAGGTCACTAATACCATCTAAGTAGTTGATTCATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTT";
    let donor_seq = [b"GGGC".to_vec(), p1.to_vec(), b"CCCCCC".to_vec(), p2.to_vec(), b"GTTT".to_vec()].concat();
    let donor = source(donor_seq, true).with_features(vec![site("attP1",4,232,1),site("attP2",242,232,2)]);
    let insert = source([cloning::ATT_B1.to_vec(), b"ACGT".repeat(25), cloning::ATT_B2.to_vec()].concat(),false);
    // A controlled annotation-level fixture, not a claim of native-site sequence
    // fidelity. Installed-vector tests below use real attR/attP sequence sources.
    let dest = source(b"C".repeat(400), true).with_features(vec![site("attR1",30,100,1),site("attR2",250,100,2)]);
    let docs = vec![DocSnap::new("attb",insert),DocSnap::new("donor",donor),DocSnap::new("destination",dest)];
    let mut bp = CloningDlg::new(Method::GatewayBpLr(1), docs, Some(0));
    bp.fragments[0].doc = Some(0); bp.fragments[0].use_pcr = false;
    bp.vector.doc = Some(1); bp.destination.doc = Some(2);
    let entry = bp.build_gateway_bp().unwrap();
    assert_eq!(persisted(&mut bp).seq.len(),230);
    bp.tab = 2;
    let intermediate = bp.build_view_doc().unwrap(); assert_eq!(intermediate.len(),entry.seq.len());
    assert_eq!(intermediate.sequence(),entry.seq);
    bp.method = Method::GatewayBp; assert_eq!(persisted(&mut bp).seq,entry.seq);
    let entry = source(entry.seq,true).with_features(entry.features);
    bp.docs.push(DocSnap::new("entry",entry)); bp.method = Method::GatewayLr(1);
    bp.vector.doc = Some(2); bp.fragments[0].doc = Some(3);
    assert_eq!(persisted(&mut bp).seq.len(),230);
}

#[test]
fn installed_vectors_construct_real_backbone_products_when_resources_are_available() {
    let root = Path::new("/opt/gslbiotech/snapgene/resources/Plasmids");
    if !root.is_dir() { eprintln!("Skipping installed-vector success matrix: SnapGene resources not available"); return; }
    let gene = b"ACGTTGCAAGTCCTGACTGAGTACCATGGCATCGTACGAT".repeat(3);
    let cases = [
        (Method::TaCloning, "TA Cloning Vectors/pGEM-T (linearized).dna"),
        (Method::TopoTa, "TOPO TA Vectors/pCR2.1-TOPO (linearized).dna"),
        (Method::GcCloning, "GC Cloning Vectors/pGC Blue (linearized).dna"),
        (Method::TopoBlunt, "TOPO Blunt Vectors/pCR-Blunt II-TOPO (linearized).dna"),
        (Method::TopoDirectional, "TOPO Directional Vectors/pENTR D-TOPO.dna"),
    ];
    for (method, path) in cases {
        let path = root.join(path);
        if !path.is_file() { eprintln!("Skipping optional vector {}", path.display()); continue; }
        let vector = doc_from_path(&path).unwrap();
        assert!(vector.src.seq.len() > 1500);
        let seq = if method == Method::TopoDirectional { [b"CACC".to_vec(),gene.clone()].concat() } else { gene.clone() };
        let mut dlg = CloningDlg::new(method, vec![vector,doc(seq,false)],Some(1));
        dlg.vector.doc = Some(0); dlg.fragments[0].doc = Some(1); dlg.fragments[0].use_pcr = false;
        dlg.fragments[0].add_a = matches!(method,Method::TaCloning|Method::TopoTa|Method::GcCloning);
        let product = persisted(&mut dlg);
        assert!(product.seq.len() > 1600);
    }
    let backbone_path = root.join("Gateway Destination Vectors/pET160-DEST.dna");
    if backbone_path.is_file() {
        for method in [Method::Gibson,Method::InFusion,Method::NeBuilder] {
            let vector = doc_from_path(&backbone_path).unwrap(); let n = vector.src.seq.len();
            let mut dlg = CloningDlg::new(method,vec![vector,doc(gene.clone(),false)],Some(0));
            dlg.vector.doc = Some(0);dlg.vector.region = Some((0,n));
            dlg.fragments[0].doc = Some(1);dlg.fragments[0].region = Some((0,gene.len()));
            design_homology_primers(&mut dlg).unwrap();
            assert_eq!(persisted(&mut dlg).seq.len(),n+gene.len());
        }
    }
    let gg_backbone=root.join("GC Cloning Vectors/pGC Blue (linearized).dna");
    if gg_backbone.is_file() {
        let mut vector=doc_from_path(&gg_backbone).unwrap();
        vector.src=cloning::blunt_ends(&vector.src).unwrap();
        vector.bytes=None;
        let n=vector.src.seq.len();
        let mut dlg=CloningDlg::new(Method::GoldenGate,vec![vector,doc(gene.clone(),false)],Some(0));
        dlg.vector.doc=Some(0);dlg.vector.region=Some((0,n));dlg.vector.gg_mode=GoldenGateMode::PcrDigest;
        dlg.fragments[0].doc=Some(1);dlg.fragments[0].region=Some((0,gene.len()));
        dlg.gg_overhangs=vec!["AATG".into(),"GCTT".into()];
        golden_gate::design(&mut dlg).unwrap();
        assert_eq!(persisted(&mut dlg).seq.len(),n+gene.len()+8);
    }
    let donor = root.join("Gateway Donor Vectors/pDONR221.dna");
    let dest = root.join("Gateway Destination Vectors/pDEST17.dna");
    if donor.is_file() && dest.is_file() {
        let insert = doc([cloning::ATT_B1.to_vec(),gene,cloning::ATT_B2.to_vec()].concat(),false);
        let docs = vec![insert,doc_from_path(&donor).unwrap(),doc_from_path(&dest).unwrap()];
        let mut dlg = CloningDlg::new(Method::GatewayBpLr(1),docs,Some(0));
        dlg.vector.doc=Some(1);dlg.destination.doc=Some(2);dlg.fragments[0].doc=Some(0);dlg.fragments[0].use_pcr=false;
        assert!(persisted(&mut dlg).seq.len()>2000);
    }
}

#[test]
fn ta_gc_topo_and_gateway_pcr_buttons_produce_valid_primers_not_vector_homology_tails() {
    use hx_formats::dna::EndProperties;
    let gene=b"ACGTTGCAAGTCCTGACTGAGTACCATGGCATCGTACGAT".repeat(3);
    for method in [Method::TaCloning,Method::GcCloning,Method::TopoTa,Method::TopoBlunt,Method::TopoDirectional] {
        let vector=match method {
            Method::GcCloning=>source(b"GAAAACCCCGGGGTTTTC".to_vec(),false).with_ends(EndProperties{upstream_stickiness:-1,downstream_stickiness:-1,..Default::default()}),
            Method::TopoBlunt=>source(b"AAAACCCCGGGGTTTT".to_vec(),false),
            Method::TopoDirectional=>source(b"AAAACCCCGGGGTTTTCACC".to_vec(),false).with_ends(EndProperties{downstream_stickiness:4,..Default::default()}),
            _=>source(b"AAAAACCCCGGGGTTTTT".to_vec(),false).with_ends(EndProperties{upstream_stickiness:-1,downstream_stickiness:-1,..Default::default()}),
        };
        let mut dlg=CloningDlg::new(method,vec![DocSnap::new("vector",vector),doc(gene.clone(),false)],Some(1));
        dlg.vector.doc=Some(0);dlg.fragments[0].doc=Some(1);dlg.fragments[0].region=Some((0,gene.len()));
        design_end_primers(&mut dlg).unwrap();
        let source=end_cloning_fragment(&dlg.docs[1].src,&dlg.fragments[0],method).unwrap();
        assert!(source.seq.windows(gene.len()).any(|w| w==gene));
        assert!(persisted(&mut dlg).seq.len()>gene.len());
    }
    let mut dlg=CloningDlg::new(Method::GatewayBp,vec![doc(gene.clone(),false)],Some(0));
    dlg.fragments[0].doc=Some(0);dlg.fragments[0].region=Some((0,gene.len()));
    design_end_primers(&mut dlg).unwrap();
    let pick=&dlg.fragments[0];
    let pcr=cloning::homology::amplify(&dlg.docs[0].src,(0,gene.len()),false,pick.primer1.as_bytes(),pick.primer2.as_bytes()).unwrap();
    assert!(pcr.seq.starts_with(cloning::ATT_B1));assert!(pcr.seq.ends_with(cloning::ATT_B2));
    assert_eq!(pcr.seq.len(),gene.len()+50);
}
