use super::*;

#[test]
fn selected_gateway_source_identity_invalidates_the_visible_document() {
    let docs=vec![DocSnap::new("same name",Source::new("first",b"A".repeat(80),false)),
        DocSnap::new("same name",Source::new("second",b"C".repeat(80),false))];
    let mut dlg=CloningDlg::new(Method::GatewayBpLr(2),docs,None);
    dlg.tab=0;dlg.fragments[0].doc=Some(0);
    let ctx=egui::Context::default();dlg.ensure_preview(&ctx);
    let key=dlg.view_doc_key();
    assert_eq!(dlg.build_view_doc().unwrap().sequence(),b"A".repeat(80));
    dlg.fragments[0].doc=Some(1);
    assert_ne!(dlg.view_doc_key(),key,"equal names and lengths must not keep the previous source canvas");
    dlg.ensure_preview(&ctx);
    assert_eq!(dlg.preview_app.as_ref().unwrap().doc.as_ref().unwrap().sequence(),b"C".repeat(80));
    dlg.tab=1;dlg.gateway_donors[0].doc=Some(0);dlg.gateway_donors[1].doc=Some(1);
    let key=dlg.view_doc_key();dlg.frag_index=1;
    assert_ne!(dlg.view_doc_key(),key,"switching donors must refresh the canvas even with equal names and lengths");
    dlg.ensure_preview(&ctx);
    assert_eq!(dlg.preview_app.as_ref().unwrap().doc.as_ref().unwrap().sequence(),b"C".repeat(80));
}

#[test]
fn standard_gateway_source_is_rendered_instead_of_a_retained_custom_file() {
    let source=Source::new("custom",b"A".repeat(80),false);
    let bytes=source_bytes(&source,"custom");
    let mut dlg=CloningDlg::new(Method::GatewayBpLr(2),vec![DocSnap::new("custom",source).with_bytes(bytes)],None);
    dlg.tab=1;
    dlg.gateway_donors[0]=Pick{doc:Some(0),use_std:true,std_name:Some("standard".into()),
        std_src:Some(Source::new("standard",b"C".repeat(80),true)),..Default::default()};
    assert_eq!(dlg.build_view_doc().unwrap().sequence(),b"C".repeat(80));
}

#[test]
fn gateway_source_mode_and_separate_primer_fields_have_distinct_cache_keys() {
    let mut dlg=CloningDlg::new(Method::GatewayBpLr(2),vec![],None);
    let key=dlg.state_key();dlg.gateway_donors[0].use_std=true;
    assert_ne!(dlg.state_key(),key,"a source-mode toggle must invalidate a ready product");
    dlg.fragments[0].primer1="AC".into();dlg.fragments[0].primer2="GT".into();
    let key=dlg.state_key();dlg.fragments[0].primer1="ACG".into();dlg.fragments[0].primer2="T".into();
    assert_ne!(dlg.state_key(),key,"different primer pairs must not collide by concatenation");
}

#[test]
fn direct_multisite_attb_details_follow_the_selected_fragment_and_numbered_sites() {
    for n in 2..=4 {
        let docs=(0..n).map(|i| {
            let (left,right)=cloning::gateway::primer_sites(n,i).unwrap();
            DocSnap::new(format!("attB {}",i+1),Source::new("PCR",
                [left,b"ACGT".repeat(20),hx_algo::reverse_complement(&right)].concat(),false))
        }).collect();
        let mut dlg=CloningDlg::new(Method::GatewayBpLr(n),docs,None);dlg.tab=0;
        for i in 0..n {dlg.fragments[i].doc=Some(i);dlg.fragments[i].use_pcr=false;}
        dlg.fragments[0].use_pcr=true;dlg.frag_index=1;
        assert!(shows_attb_detail(&dlg),"direct mode belongs to the selected fragment, not always the first");
        dlg.fragments[0].use_pcr=false;
        for i in 0..n {
            dlg.frag_index=i;
            let ctx=egui::Context::default();
            let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(Rect::from_min_size(Pos2::ZERO,Vec2::new(1800.0,1100.0))),..Default::default()},|ui|{
                paint_attb_detail(ui,ui.max_rect(),&dlg);
            });
            let text=output.shapes.iter().filter_map(|s|match &s.shape {
                egui::epaint::Shape::Text(t)=>Some(t.galley.job.text.as_str()),_=>None
            }).collect::<Vec<_>>().join("\n");
            let &(a,ad,b,bd)=cloning::gateway::layout(n).unwrap().get(i).unwrap();
            let label=|number,direction|if direction==1 {format!("attB{number} →")}else{format!("← attB{number}")};
            assert!(text.contains(&label(a,ad)),"missing actual left site/direction: {text}");
            assert!(text.contains(&label(b,bd)),"missing actual right site/direction: {text}");
            output.textures_delta.clear();
        }
    }
}

fn load(relative:&str)->Option<DocSnap> {
    doc_from_path(&Path::new("/opt/gslbiotech/snapgene/resources/Plasmids").join(relative))
}

#[test]
fn multisite_two_to_four_bp_lr_and_lr_use_every_selected_source() {
    let Some(destination)=load("Gateway Destination Vectors/pDEST17.dna") else {eprintln!("Skipping real-vector MultiSite UI test: resources absent");return;};
    for n in 2..=4 {
        let names=match n {2=>vec!["pDONR221 P1-P5r.dna","pDONR221 P5-P2.dna"],
            3=>vec!["pDONR221 P1-P4.dna","pDONR221 P4r-P3r.dna","pDONR221 P3-P2.dna"],
            _=>vec!["pDONR221 P1-P5r.dna","pDONR221 P5-P4.dna","pDONR221 P4r-P3r.dna","pDONR221 P3-P2.dna"]};
        let mut docs=vec![DocSnap::new("destination",destination.src.clone()).with_bytes(destination.bytes.clone().unwrap())];
        for name in names {docs.push(load(&format!("Gateway Donor Vectors/{name}")).unwrap());}
        for i in 0..n {docs.push(DocSnap::new(format!("gene {}",i+1),Source::new(format!("gene {}",i+1),
            format!("ATG{}TAA",["ACGTCGA","TGACTGC","CGTAGTC","GACAGTC"][i].repeat(15)).into_bytes(),false)));}
        let mut dlg=CloningDlg::new(Method::GatewayBpLr(n),docs,None);
        dlg.destination.doc=Some(0);
        for i in 0..n {
            dlg.gateway_donors[i].doc=Some(i+1);
            dlg.fragments[i].doc=Some(n+1+i);
            dlg.fragments[i].region=Some((0,dlg.docs[n+1+i].src.seq.len()));
            dlg.frag_index=i;design_end_primers(&mut dlg).unwrap();
        }
        dlg.refresh();let product=dlg.preview.clone().unwrap_or_else(||panic!("{:?}",dlg.error));
        assert!(dlg.primers_ready());assert!(product.circular);assert_eq!(product.features.iter().filter(|f|f.name.starts_with("attB")).count(),n+1);
        // Source-mode changes must also invalidate an already successful product,
        // not only produce a different string key in an empty dialog.
        dlg.gateway_donors[0].use_std=true;dlg.refresh();
        assert!(dlg.preview.is_none());assert!(dlg.error.is_some());
        dlg.gateway_donors[0].use_std=false;dlg.refresh();
        assert_eq!(dlg.preview.as_ref().unwrap().seq,product.seq);
        let mut entries=Vec::new();
        for i in 0..n {
            dlg.frag_index=i;dlg.tab=2;
            let entry=dlg.build_gateway_entry(i).unwrap();
            let preview=dlg.build_view_doc().unwrap();assert_eq!(preview.sequence(),entry.seq);
            entries.push(Source::new(format!("Entry {}",i+1),entry.seq,true).with_features(entry.features));
        }
        let scheme=reaction_scheme::model(&dlg);assert_eq!(scheme.rows.len(),n+1);
        assert_eq!(scheme.rows.last().unwrap().inputs.len(),n+1);
        for entry in entries {
            dlg.docs.push(DocSnap::new("generated Entry",entry));
        }
        let mut independent_docs=Vec::new();
        independent_docs.push(DocSnap::new("destination",destination.src.clone()));
        independent_docs.extend(dlg.docs.iter().rev().take(n).map(|d|DocSnap::new(&d.name,d.src.clone())));
        let auto=CloningDlg::new(Method::GatewayLr(n),independent_docs,None);
        assert!(auto.fragments.iter().all(|p|p.doc.is_some()));
        assert_eq!(auto.build().unwrap().seq,product.seq,"default sources are ordered by actual att identity, not file ordering");
        dlg.method=Method::GatewayLr(n);dlg.vector.doc=Some(0);dlg.tab=2;
        let offset=dlg.docs.len()-n;
        for i in 0..n {dlg.fragments[i].doc=Some(offset+i);dlg.fragments[i].use_pcr=false;}
        let lr=dlg.build().unwrap();assert_eq!(lr.seq,product.seq);
        let bytes=cloning::product_bytes(&lr,"MultiSite Expression");
        let doc=crate::doc::Document::from_bytes(bytes.clone()).unwrap();assert_eq!(doc.sequence(),product.seq);assert_eq!(doc.session.save_bytes().unwrap(),bytes);
        for method in [Method::GatewayLr(n),Method::GatewayBpLr(n)] {
            dlg.method=method;
            if method.is_bp() {
                for i in 0..n {dlg.fragments[i].doc=Some(n+1+i);dlg.fragments[i].use_pcr=true;}
            }
            for tab in 0..dlg.tabs().len() {
                dlg.tab=tab;dlg.refresh();
                let ctx=egui::Context::default();
                let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(Rect::from_min_size(Pos2::ZERO,Vec2::new(1800.0,1100.0))),..Default::default()},|ui|{
                    let ctx=ui.ctx().clone();assert!(draw_dialog(&mut dlg,&ctx,ui).is_none());
                });
                assert!(!output.shapes.is_empty());output.textures_delta.clear();
            }
        }
        // Incomplete automatic primer design is transactional across all PCR
        // fragments, so a later missing source cannot half-overwrite earlier fields.
        let saved=dlg.fragments.clone();let selected=dlg.frag_index;
        dlg.method=Method::GatewayBpLr(n);dlg.fragments[n-1].doc=None;
        dlg.fragments[0].primer1="ACGT".into();dlg.fragments[0].primer2="TGCA".into();
        let before=dlg.fragments[0].primer1.clone();design_primers(&mut dlg);
        assert_eq!(dlg.fragments[0].primer1,before);assert_eq!(dlg.fragments[0].primer2,"TGCA");assert_eq!(dlg.frag_index,selected);
        dlg.fragments=saved;
        dlg.method=Method::GatewayLr(n);dlg.vector.doc=Some(0);
        for i in 0..n {dlg.fragments[i].doc=Some(offset+i);dlg.fragments[i].use_pcr=false;}
        let key=dlg.state_key();dlg.gateway_donors[0].doc=None;assert_ne!(key,dlg.state_key());
        // Reverse selected input direction must be checked, not silently ignored.
        dlg.fragments[0].reversed=true;assert!(dlg.build().is_err());dlg.fragments[0].reversed=false;
        dlg.fragments[1].doc=None;assert!(dlg.build().is_err());
        if let Ok(root)=std::env::var("HELIX_MULTISITE_OUTPUT_DIR") {
            let dir=Path::new(&root).join(format!("{n}-fragments"));std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("hx-product.dna"),bytes).unwrap();
            for (i,snapshot) in dlg.docs.iter().enumerate() {
                std::fs::write(dir.join(format!("source-{i}.dna")),snapshot.bytes.clone().unwrap_or_else(||source_bytes(&snapshot.src,&snapshot.name))).unwrap();
            }
        }
    }
}
