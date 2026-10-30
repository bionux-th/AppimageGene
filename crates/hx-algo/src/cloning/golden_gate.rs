//! Golden Gate uses Type IIS double-stranded cuts and matching cohesive ends,
//! never a generic homology-overlap assembly.
use super::{digest_fragment, restriction_cuts, restriction_product, Product, RestrictionCut, Source};

/// Complete digest fragments which exclude the selected enzyme's recognition
/// sequence. Boundaries remain explicit, so repeated sites aren't silently ignored.
pub fn fragments(src: &Source, enzyme: &hx_enzymes::Enzyme) -> Vec<(RestrictionCut, RestrictionCut, Source)> {
    let cuts = restriction_cuts(&src.seq, src.circular, enzyme);
    let mut pairs: Vec<_> = cuts.windows(2).map(|w| (w[0], w[1])).collect();
    if src.circular && !cuts.is_empty() { pairs.push((*cuts.last().unwrap(), cuts[0])); }
    pairs.into_iter().filter_map(|(a, b)| {
        let fragment = digest_fragment(src, a, b, false, false).ok()?;
        let rc = super::revcomp_str(enzyme.site.as_bytes());
        let contains = !hx_enzymes::find_site_positions(&fragment.seq, enzyme.site.as_bytes(), false).is_empty()
            || !hx_enzymes::find_site_positions(&fragment.seq, &rc, false).is_empty();
        (!contains).then_some((a, b, fragment))
    }).collect()
}

/// Assemble prepared fragments in the user's specified order and orientation.
/// Circular closing ends are checked too. The ligase gate requires phosphorylation.
pub fn assemble(prepared: &[Source], enzyme: &hx_enzymes::Enzyme) -> Result<Product, String> {
    if prepared.len() < 2 { return Err("Please choose a Vector and at least one Fragment.".into()); }
    let width = enzyme.site_len() as i32 - enzyme.reverse_cut as i32 - enzyme.cut as i32;
    if width <= 0 { return Err("Select a Type IIS enzyme that creates 5' overhangs.".into()); }
    for (i, source) in prepared.iter().enumerate() {
        if source.circular || source.ends.upstream_stickiness != width || source.ends.downstream_stickiness != width {
            return Err(format!("{} must have {width}-base 5' overhangs for {}.",
                if i == 0 { "The linearized Vector".into() } else { format!("Fragment {i}") }, enzyme.name));
        }
    }
    let product = restriction_product(prepared)?;
    if !restriction_cuts(&product.seq, true, enzyme).is_empty() {
        return Err(format!("The assembled product contains a {} site and would be cut again. Choose different junctions or remove the internal site.", enzyme.name));
    }
    Ok(product)
}

#[cfg(test)]
mod tests {
    use super::*;
    use hx_formats::dna::EndProperties;

    fn part(label: &str, seq: &str, ends: (&str, &str)) -> Source {
        Source::new(label, format!("{}{}{}", ends.0, seq, ends.1).into_bytes(), false).with_ends(EndProperties {
            upstream_stickiness: 4, downstream_stickiness: 4,
            upstream_modification: "FivePrimePhosphorylated".into(), downstream_modification: "FivePrimePhosphorylated".into(),
            ..Default::default()
        })
    }

    #[test]
    fn cohesive_ends_are_not_homology_overlaps() {
        let e = hx_enzymes::by_name("BsaI").unwrap();
        let vector = part("vector", "GGGGGG", ("AATG", "GCTT"));
        let insert = part("insert", "CCCCCC", ("GCTT", "AATG"));
        let product = assemble(&[vector.clone(), insert.clone()], e).unwrap();
        assert_eq!(product.seq, b"AATGGGGGGGGCTTCCCCCC");
        assert!(product.circular);
        assert!(assemble(&[vector, insert.flipped()], e).is_err());
    }

    #[test]
    fn digest_discards_recognition_sites_and_preserves_real_overhangs() {
        let src = Source::new("insert", b"GGTCTCAGCTTCCCCCCAATGTGAGACC".to_vec(), false);
        let e = hx_enzymes::by_name("BsaI").unwrap();
        let parts = fragments(&src, e);
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0].2.seq, b"GCTTCCCCCCAATG");
        assert_eq!(parts[0].2.top_seq(), b"GCTTCCCCCC");
        assert_eq!(parts[0].2.ends.upstream_stickiness, 4);
        assert_eq!(parts[0].2.ends.downstream_stickiness, 4);
    }

    #[test]
    fn a_residual_site_is_not_a_successful_golden_gate_product() {
        let e = hx_enzymes::by_name("BsaI").unwrap();
        assert!(assemble(&[part("vector", "GGGG", ("AATG", "GCTT")),
            part("insert", "CCGGTCTCAAAAAACC", ("GCTT", "AATG"))], e).unwrap_err().contains("cut again"));
    }
}
