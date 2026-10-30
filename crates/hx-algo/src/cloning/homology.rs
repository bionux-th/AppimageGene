//! Homology assembly: actual PCR products and a required homology at every join.
use super::{assemble, overlap_len, Product, Slice, Source};

/// Simulate an explicitly selected PCR region from complete 5′→3′ primers.
/// Only the part preceding the annealing suffix is added as a tail.
pub fn amplify(source: &Source, region: (usize, usize), reversed: bool,
    forward: &[u8], reverse: &[u8]) -> Result<Slice, String>
{
    let (a, b) = region;
    if source.seq.is_empty() || a >= source.seq.len() || b > source.seq.len()
        || (!source.circular && a >= b)
    { return Err("Please select a valid region to amplify.".into()); }
    let mut region = source.slice(a, b, reversed);
    if forward.is_empty() || reverse.is_empty() { return Err("Please specify two PCR primers.".into()); }
    if !forward.iter().chain(reverse).all(|b| hx_enzymes::iupac_bits(*b).is_some()) {
        return Err("PCR primers must contain only IUPAC DNA bases.".into());
    }
    let anneal = |primer: &[u8], template: &[u8]| (12..=primer.len().min(template.len())).rev()
        .find(|&n| primer[primer.len() - n..].eq_ignore_ascii_case(&template[..n]));
    let f = anneal(forward, &region.seq).ok_or("Primer 1 does not anneal to the beginning of the selected region.")?;
    let r = anneal(reverse, &crate::reverse_complement(&region.seq))
        .ok_or("Primer 2 does not anneal to the end of the selected region.")?;
    region.prepend(&forward[..forward.len() - f]);
    region.append(&crate::reverse_complement(&reverse[..reverse.len() - r]));
    Ok(region)
}

/// Every neighboring pair, including the closing pair, must share an overlap.
/// The overlap upper bound is the fragment length, not an arbitrary 20/40 bp cap.
pub fn assemble_circular(parts: &[Slice], min_overlap: usize) -> Result<Product, String> {
    if parts.len() < 2 { return Err("Please choose at least two assembly fragments.".into()); }
    let mut overlaps = Vec::new();
    for i in 0..parts.len() {
        let j = (i + 1) % parts.len();
        let a = &parts[i].seq; let b = &parts[j].seq;
        if a.is_empty() || b.is_empty() { return Err("An assembly fragment is empty.".into()); }
        let max = a.len().min(b.len()).saturating_sub(1);
        let overlap = overlap_len(a, b, min_overlap, max);
        if overlap == 0 {
            return Err(format!("No suitable overlap at the junction between Fragments {} and {} (at least {} bp required).", i + 1, j + 1, min_overlap));
        }
        overlaps.push(overlap);
    }
    // Two overlapping joins must not silently remove an entire fragment.
    for i in 0..parts.len() {
        let left = overlaps[(i + parts.len() - 1) % parts.len()];
        if left + overlaps[i] >= parts[i].seq.len() {
            return Err(format!("The overlaps consume all of Fragment {}. Select a region with a non-overlapping core.", i + 1));
        }
    }
    // Use exactly the validated overlaps. Re-running generic overlap discovery
    // after trimming can choose a different closing overlap on repeat-rich DNA.
    let trimmed:Vec<_>=parts.iter().zip(&overlaps).map(|(part,&overlap)| part.sub(0,part.seq.len()-overlap)).collect();
    let mut product=assemble(&trimmed,true,0,0);
    for (junction,&overlap) in product.junctions.iter_mut().zip(&overlaps) {
        junction.overlap=overlap;
        junction.seq_end=junction.seq_start+overlap;
    }
    product.junctions.push(super::Junction{at:product.seq.len(),overlap:*overlaps.last().unwrap(),
        seq_start:product.seq.len(),seq_end:*overlaps.last().unwrap()});
    Ok(product)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(s: &str) -> Slice { Source::new("part", s.as_bytes().to_vec(), false).slice(0, s.len(), false) }

    #[test]
    fn both_internal_and_circular_joins_are_required() {
        let a = "ACGTTGCAAGTCCTGACTGA";
        let b = "GTACCATGGCATCGTACGAT";
        let vector = part(&format!("{b}{}{}", "G".repeat(40), a));
        let insert = part(&format!("{a}{}{}", "C".repeat(40), b));
        let p = assemble_circular(&[vector.clone(), insert], 15).unwrap();
        assert_eq!(p.seq.len(), 120);
        assert!(p.circular);
        let bad = part(&format!("{a}{}{}", "C".repeat(40), "T".repeat(20)));
        assert!(assemble_circular(&[vector, bad], 15).unwrap_err().contains("Fragments 2 and 1"));
    }

    #[test]
    fn full_primers_do_not_duplicate_their_annealed_bases() {
        let source = Source::new("PCR", b"ACGTTGCAAGTCCTGACTGAGGGGGTACCATGGCATCGTACGAT".to_vec(), false);
        let forward = b"TTTTACGTTGCAAGTCCTGACTGA";
        let reverse = crate::reverse_complement(b"GTACCATGGCATCGTACGATCCCC");
        let p = amplify(&source, (0, source.seq.len()), false, forward, &reverse).unwrap();
        assert_eq!(p.seq, [b"TTTT".to_vec(), source.seq.clone(), b"CCCC".to_vec()].concat());
        assert!(amplify(&source, (0, source.seq.len()), false, b"AAAA", &reverse).is_err());
    }

    #[test]
    fn long_homology_is_not_capped_at_forty_bases() {
        let a="ACGTTGCAAGTCCTGACTGAGTACCATGGCATCGTACGATCGTACGCTAGT";
        let b="GTCAGATCGTACGTGACCATGCTAGCTAGGTCAGTACGATCGATTCGAGTCA";
        let first=part(&format!("{a}{}{b}","C".repeat(80)));
        let second=part(&format!("{b}{}{a}","G".repeat(80)));
        let product=assemble_circular(&[first,second],15).unwrap();
        assert_eq!(product.seq.len(),a.len()+b.len()+160);
        assert_eq!(product.junctions.len(),2);
        assert_eq!(product.junctions[0].overlap,b.len());
        assert_eq!(product.junctions[1].overlap,a.len());
    }
}
