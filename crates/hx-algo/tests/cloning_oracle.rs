//! SnapGene 8.2.3 GUI 对照，2026-10-25。受控合成输入，无许可证/用户序列。
//! 证据步骤、输出摘要和覆盖边界见 tests/fixtures/cloning/README.md。
use hx_algo::cloning::{digest_fragment, ligate_reaction, product_bytes, restriction_cuts, restriction_product, Source};
use hx_formats::dna::EndProperties;

fn sequence(file: &str) -> Vec<u8> {
    file.lines().filter(|l| !l.starts_with('>')).flat_map(|l| l.trim().bytes()).collect()
}

#[test]
fn restriction_ecori_bamhi_matches_sg_sequence_without_rotation_normalization() {
    let vseq = format!("{}GAATTC{}GGATCC{}", "A".repeat(100), "C".repeat(30), "ACCG".repeat(170));
    let iseq = format!("{}GAATTC{}GGATCC{}", "T".repeat(20), "ACCA".repeat(40), "T".repeat(20));
    let v = Source::new("vector", vseq.into_bytes(), true);
    let i = Source::new("insert", iseq.into_bytes(), false);
    let ecori = hx_enzymes::by_name("EcoRI").unwrap();
    let bamhi = hx_enzymes::by_name("BamHI").unwrap();
    let vc1 = restriction_cuts(&v.seq, true, ecori)[0];
    let vc2 = restriction_cuts(&v.seq, true, bamhi)[0];
    let ic1 = restriction_cuts(&i.seq, false, ecori)[0];
    let ic2 = restriction_cuts(&i.seq, false, bamhi)[0];
    let vector = digest_fragment(&v, vc2, vc1, false, false).unwrap();
    let insert = digest_fragment(&i, ic1, ic2, false, false).unwrap();
    assert_eq!((vector.display_len(), insert.display_len()), (786, 166));
    let p = restriction_product(&[vector, insert]).unwrap();
    assert_eq!(p.seq, sequence(include_str!("fixtures/cloning/restriction-ecori-bamhi.fasta")));
    assert_eq!(p.seq.len(), 952);
    assert!(p.circular);
    let f = hx_formats::DnaFile::parse(&product_bytes(&p, "Cloned")).unwrap();
    assert_eq!(f.end_properties().upstream_modification, "FivePrimePhosphorylated");
    assert_eq!(f.end_properties().downstream_modification, "FivePrimePhosphorylated");
    // 同轮 SG 对照：两侧的 EcoRI/BamHI 5′ 突出端全部填平。
    let vector = digest_fragment(&v, vc2, vc1, true, false).unwrap();
    let insert = digest_fragment(&i, ic1, ic2, true, false).unwrap();
    assert_eq!((vector.display_len(), insert.display_len()), (790, 170));
    let blunt = restriction_product(&[vector, insert]).unwrap();
    let expected = [p.seq[..105].to_vec(), b"AATT".to_vec(), p.seq[105..271].to_vec(), b"GATC".to_vec(), p.seq[271..].to_vec()].concat();
    assert_eq!(blunt.seq, expected);
    assert_eq!(blunt.seq.len(), 960);
}

#[test]
fn linear_ligation_matches_sg_phosphorylation_gate_and_native_sequence() {
    let mut left = Source::new("left", format!("T{}AATT", "G".repeat(79)).into_bytes(), false).with_ends(EndProperties {
        upstream_stickiness: -1, downstream_stickiness: 4, ..Default::default()
    });
    let right = Source::new("right", format!("AATT{}AA", "C".repeat(80)).into_bytes(), false).with_ends(EndProperties {
        upstream_stickiness: 4, downstream_stickiness: 2, ..Default::default()
    });
    assert_eq!(ligate_reaction(&[left.clone(), right.clone()], false).unwrap_err(),
        "For ligation, at least one DNA end must be 5' phosphorylated at the junction between Fragments 1 and 2.");
    left.ends.downstream_modification = "FivePrimePhosphorylated".into();
    let p = ligate_reaction(&[left, right], false).unwrap();
    assert_eq!(p.seq, sequence(include_str!("fixtures/cloning/ligation-phosphorylated.fasta")));
    assert_eq!(p.display_len(), 163);
    assert_eq!(p.seq.len(), 166);
    assert!(!p.circular);
    assert_eq!((p.ends.upstream_stickiness, p.ends.downstream_stickiness), (-1, 2));
    assert_eq!((&*p.ends.upstream_modification, &*p.ends.downstream_modification), ("Unmodified", "Unmodified"));
    let file = hx_formats::DnaFile::parse(&product_bytes(&p, "Ligation")).unwrap();
    assert_eq!(file.sequence().unwrap().sequence, p.seq);
    assert_eq!(file.end_properties(), p.ends);
}
