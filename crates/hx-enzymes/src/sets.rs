//! Commercial parent sets and sequence-dependent filters; counts are derived.
use crate::catalog::CATALOG;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Parent { All, Nonredundant, Neb, Thermo, Supplier(u8) }
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Filter { All, Unique, UniqueDual, SixPlus, UniqueSixPlus, TypeIis, GoldenGate, Nicking }
pub const DEFAULT_SUPPLIER_ORDER: &[u8] = &[11,19,16,23,24,20,3,15,0,18,14,22,10,4,9,12,6,7,8,5,25,1,2,21,13,17];
pub const SUPPLIER_NAMES: &[&str] = &["Alator Biosciences","American Allied Biochemical","CinnaGen","Enzynomics","EURx","GeNei","Jena Bioscience","MilliporeSigma","Minotech Biotechnology","Molecular Biology Resources - CHIMERx","MP Biomedicals","New England Biolabs","Nippon Gene","NZYTech","Promega","Roche","SibEnzyme","Stratagene","TaKaRa Bio (Clontech)","Thermo Fisher (Fermentas)","Thermo Fisher (Invitrogen)","Thermo Fisher (United States Biochemical Corp.)","Toyobo","USB","Vivantis","Yeasen"];

/// SG 8.2.3 default selected Golden Gate list (0x4e9f2e), distinct
/// from the 14-entry candidate list (0x4795eb).
pub const GOLDEN_GATE_DEFAULT: &[&str] = &["BbsI","BfuAI","BsaI","BsmBI","BspQI","BtgZI","Esp3I","PaqCI","SapI"];

pub fn nicking(i: usize) -> bool { CATALOG[i].cuts.iter().any(|(a,b)| a.is_none() || b.is_none()) }
pub fn type_iis(i: usize) -> bool {
    // SG 8.2.3 filter 0x1211800: one cut pair, neither cut strictly
    // inside the first recognition sequence. Boundary cuts count here.
    let e = &CATALOG[i];
    let len = e.site.split_whitespace().next().unwrap_or("").len() as i16;
    !nicking(i) && e.cuts.len() == 1 && e.cuts.iter().all(|(a,b)|
        [a,b].iter().all(|v| v.is_some_and(|v| v <= 0 || v >= len)))
}
pub fn parent_indices(parent: Parent) -> Vec<usize> {
    parent_indices_ordered(parent, DEFAULT_SUPPLIER_ORDER)
}
pub fn parent_indices_ordered(parent: Parent, order: &[u8]) -> Vec<usize> {
    let mut indices: Vec<_> = (0..CATALOG.len()).filter(|&i| match parent {
        Parent::Neb => CATALOG[i].suppliers.contains(&11), Parent::Thermo => CATALOG[i].suppliers.contains(&19),
        Parent::Supplier(id) => CATALOG[i].suppliers.contains(&id),
        _ => true,
    }).collect();
    if parent == Parent::Nonredundant {
        // SG supplier-batch selection, not one globally ranked winner per
        // equivalence group. All eligible members of one supplier survive.
        // Order is the SG default preferredSupplierOrder, recorded in docs.
        let mut state = vec![0u8; CATALOG.len()];
        for &supplier in order {
            let batch: Vec<_> = indices.iter().copied().filter(|&i| CATALOG[i].suppliers.contains(&supplier)).collect();
            for &i in &batch { if state[i] == 0 { state[i] = 1; } }
            for i in batch {
                if state[i] != 1 { continue; }
                state[i] = 2;
                for name in CATALOG[i].equivalents {
                    if let Some(j) = CATALOG.iter().position(|e| e.name == *name) {
                        if state[j] == 0 { state[j] = 3; }
                    }
                }
            }
        }
        indices.retain(|&i| state[i] == 2);
    }
    indices
}
pub fn positions(i: usize, seq: &[u8], circular: bool) -> Vec<u32> {
    let mut out = Vec::new();
    for pattern in CATALOG[i].site.split_whitespace() {
        out.extend(crate::find_site_positions(seq, pattern.as_bytes(), circular));
        let rc: Vec<_> = pattern.bytes().rev().map(|b| match b { b'A'=>b'T',b'T'=>b'A',b'C'=>b'G',b'G'=>b'C',b'R'=>b'Y',b'Y'=>b'R',b'K'=>b'M',b'M'=>b'K',b'B'=>b'V',b'V'=>b'B',b'D'=>b'H',b'H'=>b'D',x=>x }).collect();
        if rc != pattern.as_bytes() { out.extend(crate::find_site_positions(seq, &rc, circular)); }
    }
    out.sort_unstable();out.dedup();out
}
pub fn select(parent: Parent, filter: Filter, seq: &[u8], circular: bool) -> Vec<usize> {
    select_ordered(parent,filter,seq,circular,DEFAULT_SUPPLIER_ORDER)
}
pub fn select_ordered(parent: Parent, filter: Filter, seq: &[u8], circular: bool, order: &[u8]) -> Vec<usize> {
    if filter == Filter::GoldenGate {
        return CATALOG.iter().enumerate().filter(|(_,e)| GOLDEN_GATE_DEFAULT.contains(&e.name)).map(|(i,_)| i).collect();
    }
    parent_indices_ordered(parent,order).into_iter().filter(|&i| {
        let specific = CATALOG[i].site.split_whitespace().next().unwrap_or("").bytes().filter(|b| *b != b'N').count();
        match filter {
            Filter::All => true,
            Filter::SixPlus => specific >= 6,
            Filter::Unique => !nicking(i) && positions(i,seq,circular).len() * CATALOG[i].cuts.len() == 1,
            Filter::UniqueDual => !nicking(i) && matches!(positions(i,seq,circular).len() * CATALOG[i].cuts.len(),1|2),
            Filter::UniqueSixPlus => specific >= 6 && !nicking(i) && positions(i,seq,circular).len() * CATALOG[i].cuts.len() == 1,
            Filter::TypeIis => type_iis(i),
            Filter::Nicking => nicking(i),
            Filter::GoldenGate => unreachable!(),
        }
    }).collect()
}

#[cfg(test)] mod tests {
    use super::*;
    use crate::ENZYMES;
    #[test] fn sg_supplier_batches_and_filter_semantics() {
        let parent = parent_indices(Parent::Nonredundant);
        assert_eq!(parent.len(),303);
        let names: Vec<_> = parent.iter().map(|&i| CATALOG[i].name).collect();
        assert!(names.contains(&"Nb.BsmI"));
        assert!(!names.contains(&"Nb.Mva1269I"));
        // Same supplier isoschizomers are intentionally retained.
        assert!(names.contains(&"BfuAI") && names.contains(&"BspMI"));
        assert_eq!(select(Parent::Nonredundant,Filter::SixPlus,b"",false).len(),209);
        assert_eq!(select(Parent::Nonredundant,Filter::TypeIis,b"",false).len(),70);
        assert_eq!(select(Parent::Nonredundant,Filter::Nicking,b"",false).len(),12);
        let find = |name| CATALOG.iter().position(|e| e.name==name).unwrap();
        assert!(type_iis(find("DpnII")), "SG boundary-cut filter includes DpnII");
        assert!(!type_iis(find("BcgI")), "multiple cuts excluded");
        assert!(!type_iis(find("Nt.BspQI")), "nicking excluded");
        let gg: Vec<_> = select(Parent::Nonredundant,Filter::GoldenGate,b"",false).iter().map(|&i| CATALOG[i].name).collect();
        assert_eq!(gg, GOLDEN_GATE_DEFAULT);
        assert!(!gg.contains(&"BsmAI"));
        assert!(gg.contains(&"PaqCI") && gg.contains(&"BfuAI"));
    }
    #[test] fn unique_means_double_strand_cuts_not_recognition_sites() {
        let seq = b"AAAAAAAAAAAAGCTCTTCAAAAAAAAAAAATGCAAAAAAACGAAAAAAAAAAAA";
        let hits = select(Parent::All, Filter::Unique, seq, false);
        assert!(!hits.iter().any(|&i| nicking(i)));
        let multi = CATALOG.iter().position(|e| e.name == "BcgI").unwrap();
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(positions(multi, seq, false).len(),1);
        assert!(!select(Parent::All,Filter::Unique,seq,false).contains(&multi));
        assert!(select(Parent::All,Filter::UniqueDual,seq,false).contains(&multi));
    }
    #[test] fn pet28b_rep1_same_input_counts_when_user_fixture_is_available() {
        let Some(path) = std::env::var_os("HELIX_ENZYME_ORACLE_DNA") else { return; };
        let bytes = std::fs::read(path).unwrap(); let mut at = 0;
        while at + 5 <= bytes.len() {
            let ty = bytes[at]; let n = u32::from_be_bytes(bytes[at+1..at+5].try_into().unwrap()) as usize;
            if ty == 0 {
                let seq = &bytes[at+6..at+5+n]; let circular = bytes[at+5]&1 != 0;
                assert_eq!(seq.len(),9325);
                assert_eq!(select(Parent::Nonredundant,Filter::Unique,seq,circular).len(),54);
                assert_eq!(select(Parent::Nonredundant,Filter::UniqueDual,seq,circular).len(),72);
                assert_eq!(select(Parent::Nonredundant,Filter::UniqueSixPlus,seq,circular).len(),54);
                return;
            }
            at += 5+n;
        }
        panic!("Oracle is missing its DNA sequence block");
    }
    #[test] fn commercial_catalog_identity_and_supplier_counts() {
        assert_eq!(CATALOG.len(),686);assert_eq!(ENZYMES.len(),686);
        for (a,b) in CATALOG.iter().zip(ENZYMES) {assert_eq!(a.name,b.name);}
        assert_eq!(parent_indices(Parent::Neb).len(),258);
        assert_eq!(parent_indices(Parent::Thermo).len(),209);
    }
}
