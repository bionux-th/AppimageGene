//! 双链酶切与末端感知连接。坐标都是两条链共用的顶链方向、0-based 间隙。
use super::{assemble, revcomp_str, rotate_range, Product, Source};
use hx_formats::dna::EndProperties;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestrictionCut {
    pub top: i64,
    /// 环状时保留相对 top 的有符号距离，不单独取模，否则跨原点会翻转极性。
    pub bottom: i64,
}

impl RestrictionCut {
    pub fn blunt(at: usize) -> Self {
        Self { top: at as i64, bottom: at as i64 }
    }

    pub fn stickiness(self) -> i32 {
        (self.bottom - self.top) as i32
    }
}

/// 数据库 reverse_cut 是互补链 5′ 端量起的偏移。
/// 正向位点: top=p+cut, bottom=p+site_len-reverse_cut。
/// 反向位点: top=p+reverse_cut, bottom=p+site_len-cut。
pub fn restriction_cuts(seq: &[u8], circular: bool, enzyme: &hx_enzymes::Enzyme) -> Vec<RestrictionCut> {
    if let Some(i) = hx_enzymes::ENZYMES.iter().position(|e| e.name == enzyme.name) {
        if hx_enzymes::sets::nicking(i) || hx_enzymes::catalog::CATALOG[i].cuts.len() != 1 { return Vec::new(); }
    }
    let n = seq.len() as i64;
    if n == 0 { return Vec::new(); }
    let pat = enzyme.site.as_bytes();
    let rc = revcomp_str(pat);
    let mut cuts = Vec::new();
    for reverse in [false, true] {
        if reverse && rc.eq_ignore_ascii_case(pat) { continue; }
        let pattern = if reverse { rc.as_slice() } else { pat };
        for p in hx_enzymes::find_site_positions(seq, pattern, circular) {
            let (a, b) = if reverse {
                (enzyme.reverse_cut as i64, enzyme.site_len() as i64 - enzyme.cut as i64)
            } else {
                (enzyme.cut as i64, enzyme.site_len() as i64 - enzyme.reverse_cut as i64)
            };
            let top = p as i64 + a;
            let bottom = p as i64 + b;
            if circular {
                // 不支持突出端绕完整圈或更长的极短底物。
                if (bottom - top).abs() >= n { continue; }
                let normalized = top.rem_euclid(n);
                cuts.push(RestrictionCut { top: normalized, bottom: normalized + bottom - top });
            } else if (0..=n).contains(&top) && (0..=n).contains(&bottom) {
                cuts.push(RestrictionCut { top, bottom });
            }
        }
    }
    cuts.sort_by_key(|c| (c.top, c.bottom));
    cuts.dedup();
    cuts
}

/// 取两个双链切口之间的一段。相同环状切口表示打开整圈；线性不允许回绕。
/// 5′ 突出端填平，3′ 突出端切除；不能通过统一平移切点模拟。
pub fn digest_fragment(
    src: &Source, left: RestrictionCut, right: RestrictionCut, blunt: bool, reversed: bool,
) -> Result<Source, String> {
    let n = src.seq.len() as i64;
    if n == 0 { return Err("The source sequence is empty.".into()); }
    for cut in [left, right] {
        if !(0..=n).contains(&cut.top)
            || (cut.bottom - cut.top).abs() >= n
            || (!src.circular && !(0..=n).contains(&cut.bottom))
        {
            return Err("The selected cuts lie outside the source sequence.".into());
        }
        if !src.circular {
            let (a, b) = src.top_range();
            let bottom_start = src.ends.upstream_stickiness.max(0) as i64;
            let bottom_end = n - src.ends.downstream_stickiness.min(0).unsigned_abs() as i64;
            if !(a as i64..=b as i64).contains(&cut.top)
                || !(bottom_start..=bottom_end).contains(&cut.bottom)
            {
                return Err("The selected cut extends beyond a source strand.".into());
            }
        }
    }
    let mut rt = right.top;
    let mut rb = right.bottom;
    if src.circular && rt <= left.top { rt += n; rb += n; }
    if rt <= left.top || rb <= left.bottom {
        return Err("The selected cuts do not define a double-stranded fragment.".into());
    }
    if left.top.max(left.bottom) >= rt.min(rb) {
        return Err("The selected cuts leave no double-stranded region.".into());
    }
    let (lt, mut lb) = (left.top, left.bottom);
    if blunt {
        // 左侧顶链的 5′ 端保持不动；右侧填补 5′ 悬垂或切除 3′ 悬垂后
        // 顶链止于原底链切口（两侧分别处理，绝非所有切点加同一偏移）。
        lb = lt;
        rt = rb;
    }
    let a = lt.min(lb);
    let b = rt.max(rb);
    if (!src.circular && (a < 0 || b > n)) || b <= a {
        return Err("The selected cuts lie outside the source sequence.".into());
    }
    let mut seq = Vec::with_capacity((b - a) as usize);
    let mut pieces = Vec::new();
    let mut at = a;
    while at < b {
        let start = at.rem_euclid(n) as usize;
        let len = (b - at).min(n - start as i64) as usize;
        seq.extend_from_slice(&src.seq[start..start + len]);
        pieces.push((start, start + len));
        at += len as i64;
    }
    let ends = EndProperties {
        upstream_stickiness: (lb - lt) as i32,
        downstream_stickiness: (rb - rt) as i32,
        upstream_modification: if !src.circular && left.top == src.top_range().0 as i64 {
            src.ends.upstream_modification.clone()
        } else { "FivePrimePhosphorylated".into() },
        downstream_modification: if !src.circular && right.top == src.top_range().1 as i64 {
            src.ends.downstream_modification.clone()
        } else { "FivePrimePhosphorylated".into() },
        // 原序列位置 0 仍在片段中时记下其存储坐标，限制性克隆保留载体原点。
        original_first_base: if src.circular {
            let original = (-a).rem_euclid(n);
            (original < b - a).then_some(original as usize)
        } else { None },
    };
    let fragment = Source::new(src.label.clone(), seq, false)
        .with_features(src.project_features(&pieces)).with_ends(ends);
    Ok(if reversed { fragment.flipped() } else { fragment })
}

fn end_bases(src: &Source, left: bool) -> Result<&[u8], String> {
    let e = if left { src.ends.upstream_stickiness } else { src.ends.downstream_stickiness };
    let len = e.unsigned_abs() as usize;
    if len >= src.seq.len() && len != 0 {
        return Err(format!("The ends of {} do not define a double-stranded fragment.", src.label));
    }
    Ok(if left { &src.seq[..len] } else { &src.seq[src.seq.len() - len..] })
}

/// 聚合酶填补 5′ 突出端、切除 3′ 突出端；特征随裁切一起移动。
pub fn blunt_ends(src: &Source) -> Result<Source, String> {
    end_bases(src, true)?; end_bases(src, false)?;
    let (top_start, top_end) = src.top_range();
    let bottom_start = src.ends.upstream_stickiness.max(0) as usize;
    let bottom_end = src.seq.len().saturating_sub(src.ends.downstream_stickiness.min(0).unsigned_abs() as usize);
    if top_start.max(bottom_start) >= top_end.min(bottom_end) {
        return Err("The ends leave no double-stranded region.".into());
    }
    let a = src.ends.upstream_stickiness.min(0).unsigned_abs() as usize;
    let b = src.seq.len() - src.ends.downstream_stickiness.min(0).unsigned_abs() as usize;
    if a >= b { return Err("Blunting leaves no double-stranded fragment.".into()); }
    let slice = src.slice(a, b, false);
    let ends = EndProperties {
        upstream_modification: src.ends.upstream_modification.clone(),
        downstream_modification: src.ends.downstream_modification.clone(),
        ..Default::default()
    };
    Ok(Source::new(src.label.clone(), slice.seq, false).with_features(slice.features).with_ends(ends))
}

fn compatible(left: &Source, right: &Source) -> Result<usize, String> {
    let a = end_bases(left, false)?;
    let b = end_bases(right, true)?;
    if left.ends.downstream_stickiness != right.ends.upstream_stickiness
        || a.len() != b.len()
        || !a.iter().zip(b).all(|(&a, &b)| {
            // 两端均存为顶链方向：须同向相等，不再取一次反向互补。
            let a = hx_enzymes::iupac_bits(a).unwrap_or(0);
            let b = hx_enzymes::iupac_bits(b).unwrap_or(0);
            a != 0 && a == b
        })
    {
        return Err(format!("The ends of {} and {} are not compatible.", left.label, right.label));
    }
    Ok(a.len())
}

/// 序列构造原语，不判断磷酸化；真实连接模拟须调用 ligate_reaction。
/// 连接各片段的顶链，只在外侧保留底链突出端标记。
pub fn ligate(fragments: &[Source], circular: bool) -> Result<Product, String> {
    if fragments.is_empty() { return Err("Please choose a source file.".into()); }
    for f in fragments {
        if f.circular { return Err("A ligation fragment must be linearized first.".into()); }
        if f.top_seq().is_empty() { return Err("The fragment is empty.".into()); }
        end_bases(f, true)?; end_bases(f, false)?;
        let bottom_start = f.ends.upstream_stickiness.max(0) as usize;
        let bottom_end = f.seq.len().saturating_sub(f.ends.downstream_stickiness.min(0).unsigned_abs() as usize);
        let (a, b) = f.top_range();
        if a.max(bottom_start) >= b.min(bottom_end) {
            return Err("The ends leave no double-stranded region.".into());
        }
    }
    let mut overlaps = Vec::new();
    for pair in fragments.windows(2) { overlaps.push(compatible(&pair[0], &pair[1])?); }
    let closing_overlap = if circular { Some(compatible(fragments.last().unwrap(), &fragments[0])?) } else { None };
    let last = fragments.len() - 1;
    let slices: Vec<_> = fragments.iter().enumerate().map(|(i, f)| {
        let (a, b) = f.top_range();
        f.slice(if !circular && i == 0 { 0 } else { a },
            if !circular && i == last { f.seq.len() } else { b }, false)
    }).collect();
    let mut product = assemble(&slices, circular, 0, 0);
    // 底链悬垂区虽然不在本片段的顶链区间，连接后却落在相邻片段的
    // 顶链上。不能因去掉内部标记而丢掉这些碱基上的注释。
    product.features.clear();
    let total = product.seq.len() as i64;
    for (i, fragment) in fragments.iter().enumerate() {
        let top_start = fragment.top_range().0;
        let shift = product.spans[i].0 as i64
            - if !circular && i == 0 { 0 } else { top_start as i64 };
        for mut feature in fragment.project_features(&[(0, fragment.seq.len())]) {
            feature.max_run_on=feature.max_run_on.and_then(|end|{
                let mapped=end as i64-1+shift;
                if circular {Some(mapped.rem_euclid(total) as usize+1)}
                else {(0..total).contains(&mapped).then_some((mapped+1) as usize)}
            });
            let mut segments = Vec::new();
            for segment in &feature.segments {
                let a = segment.range.start as i64 + shift;
                let b = segment.range.end as i64 + shift;
                let ranges = if circular {
                    let start = a.rem_euclid(total) as usize;
                    let end = start + (b - a).min(total) as usize;
                    if end <= total as usize { vec![(start, end)] }
                    else { vec![(start, total as usize), (0, end - total as usize)] }
                } else { vec![(a.clamp(0, total) as usize, b.clamp(0, total) as usize)] };
                for (a, b) in ranges {
                    if a < b {
                        let mut segment = segment.clone();
                        segment.range = hx_formats::annot::Range::new(a, b);
                        // 打开整圈时悬垂区在存储串里出现两次，连接后落到
                        // 同一产物位置。合并同一特征中重叠的重复投影。
                        let mut index = 0;
                        while index < segments.len() {
                            let old: &hx_formats::annot::FeatureSegment = &segments[index];
                            if old.color == segment.color && old.seg_type == segment.seg_type
                                && old.translated == segment.translated
                                && old.range.start < segment.range.end && segment.range.start < old.range.end
                            {
                                let old = segments.remove(index);
                                segment.range.start = segment.range.start.min(old.range.start);
                                segment.range.end = segment.range.end.max(old.range.end);
                                index = 0;
                            } else { index += 1; }
                        }
                        segments.push(segment);
                    }
                }
            }
            if !segments.is_empty() {
                feature.segments = segments;
                product.features.push(feature);
            }
        }
    }
    for (junction, overlap) in product.junctions.iter_mut().zip(overlaps) {
        junction.overlap = overlap;
    }
    if let Some(overlap) = closing_overlap {
        product.junctions.push(super::Junction {
            at: product.seq.len(), overlap, seq_start: product.seq.len(), seq_end: 0,
        });
    }
    if !circular {
        product.ends = EndProperties {
            upstream_stickiness: fragments[0].ends.upstream_stickiness,
            downstream_stickiness: fragments[last].ends.downstream_stickiness,
            upstream_modification: fragments[0].ends.upstream_modification.clone(),
            downstream_modification: fragments[last].ends.downstream_modification.clone(),
            original_first_base: None,
        };
    }
    Ok(product)
}

/// SG 8.2.3 运行对照：每个接头至少一个参与连接的 5′ 端必须磷酸化。
/// UpstreamModification 标记顶链左端，DownstreamModification 标记底链右端。
/// 未识别的修饰不当作磷酸化；TOPO 等共价连接机制不走普通 ligase 路径。
pub fn ligate_reaction(fragments: &[Source], circular: bool) -> Result<Product, String> {
    let product = ligate(fragments, circular)?;
    let joins = fragments.len().saturating_sub(1) + usize::from(circular);
    for i in 0..joins {
        let j = (i + 1) % fragments.len();
        if fragments[i].ends.downstream_modification != "FivePrimePhosphorylated"
            && fragments[j].ends.upstream_modification != "FivePrimePhosphorylated"
        {
            return Err(format!("For ligation, at least one DNA end must be 5' phosphorylated at the junction between Fragments {} and {}.", i + 1, j + 1));
        }
    }
    Ok(product)
}

/// 限制性插入克隆：连接后尽可能保留骨架的原始位置 1（SG 实测）。
pub fn restriction_product(fragments: &[Source]) -> Result<Product, String> {
    let mut product = ligate_reaction(fragments, true)?;
    if let Some(vector) = fragments.first() {
        if let Some(first) = vector.ends.original_first_base {
            let (a, b) = vector.top_range();
            if (a..b).contains(&first) {
                rotate_product(&mut product, first - a);
            }
        }
    }
    // SG 保存的环状限制性克隆产物保留此修饰对（无外侧突出端）。
    product.ends.upstream_modification = "FivePrimePhosphorylated".into();
    product.ends.downstream_modification = "FivePrimePhosphorylated".into();
    Ok(product)
}

fn rotate_product(product: &mut Product, amount: usize) {
    let n = product.seq.len();
    if n == 0 || amount % n == 0 { return; }
    let amount = amount % n;
    product.seq.rotate_left(amount);
    for feature in &mut product.features {
        feature.max_run_on=feature.max_run_on.map(|end|(end-1+n-amount)%n+1);
        feature.segments = feature.segments.iter().flat_map(|segment| {
            rotate_range(segment.range.start, segment.range.end, amount, n).into_iter().map(|(a, b)| {
                let mut segment = segment.clone();
                segment.range = hx_formats::annot::Range::new(a, b);
                segment
            })
        }).collect();
    }
    product.spans = product.spans.iter().flat_map(|&(a, b)| rotate_range(a, b, amount, n)).collect();
    for junction in &mut product.junctions {
        junction.at = (junction.at + n - amount) % n;
        junction.seq_start = (junction.seq_start + n - amount) % n;
        junction.seq_end = (junction.seq_end + n - amount) % n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cuts(seq: &str, circular: bool, enzyme: &str) -> Vec<RestrictionCut> {
        restriction_cuts(seq.as_bytes(), circular, hx_enzymes::by_name(enzyme).unwrap())
    }

    #[test]
    fn common_cut_coordinates_are_not_same_origin_offsets() {
        assert_eq!(cuts("AAGAATTCTT", false, "EcoRI"), vec![RestrictionCut { top: 3, bottom: 7 }]);
        assert_eq!(cuts("AAGACGTCTT", false, "AatII"), vec![RestrictionCut { top: 7, bottom: 3 }]);
        assert_eq!(cuts("AACCCGGGTT", false, "SmaI"), vec![RestrictionCut::blunt(5)]);
        for e in hx_enzymes::ENZYMES {
            let d = e.site_len() as i64 - e.reverse_cut as i64 - e.cut as i64;
            match e.overhang {
                hx_enzymes::Overhang::FivePrime => assert!(d > 0, "{}", e.name),
                hx_enzymes::Overhang::ThreePrime => assert!(d < 0, "{}", e.name),
                hx_enzymes::Overhang::Blunt => assert_eq!(d, 0, "{}", e.name),
                hx_enzymes::Overhang::Iis => {},
            }
        }
    }

    #[test]
    fn nonpalindromic_equal_offsets_still_match_reverse_sites() {
        let enzyme = hx_enzymes::by_name("AciI").unwrap();
        assert_eq!(enzyme.cut, enzyme.reverse_cut);
        assert_eq!(cuts("AAGCGGTT", false, "AciI"), vec![RestrictionCut { top: 3, bottom: 5 }]);
    }

    #[test]
    fn type_iis_both_orientations_and_incomplete_linear_sites() {
        assert_eq!(cuts("GGTCTCAACGTTTT", false, "BsaI"), vec![RestrictionCut { top: 7, bottom: 11 }]);
        assert_eq!(cuts("AAAACGTTGAGACC", false, "BsaI"), vec![RestrictionCut { top: 3, bottom: 7 }]);
        assert!(cuts("GGTCTCA", false, "BsaI").is_empty());
        assert!(cuts("GAGACC", false, "BsaI").is_empty());
        assert_eq!(cuts("GAGACCAAAAAAAA", true, "BsaI"), vec![RestrictionCut { top: 9, bottom: 13 }]);
    }

    #[test]
    fn circular_cuts_keep_signed_overhang_across_origin() {
        assert_eq!(cuts("AATTCGGG", true, "EcoRI"), vec![RestrictionCut { top: 0, bottom: 4 }]);
        assert_eq!(cuts("TTCGGGAA", true, "EcoRI"), vec![RestrictionCut { top: 6, bottom: 10 }]);
        let s = Source::new("ring", b"TTCGGGAA".to_vec(), true);
        let c = cuts("TTCGGGAA", true, "EcoRI")[0];
        let f = digest_fragment(&s, c, c, false, false).unwrap();
        assert_eq!(f.seq, b"AATTCGGGAATT");
        assert_eq!(f.top_seq(), b"AATTCGGG");
        assert_eq!((f.ends.upstream_stickiness, f.ends.downstream_stickiness), (4, 4));
        assert_eq!(ligate(&[f], true).unwrap().seq, b"AATTCGGG");
        let s = Source::new("three prime ring", b"GACGTCAAAA".to_vec(), true);
        let c = cuts("GACGTCAAAA", true, "AatII")[0];
        let f = digest_fragment(&s, c, c, false, false).unwrap();
        assert_eq!(f.ends.upstream_stickiness, -4);
        assert_eq!(f.ends.downstream_stickiness, -4);
        assert_eq!(ligate(&[f], true).unwrap().seq, b"CAAAAGACGT");
    }

    #[test]
    fn cut_then_rejoin_restores_both_five_and_three_prime_sequences() {
        for (seq, enzyme) in [("AAAAGAATTCTTTT", "EcoRI"), ("AAAAGACGTCTTTT", "AatII")] {
            let src = Source::new("linear", seq.as_bytes().to_vec(), false);
            let c = cuts(seq, false, enzyme)[0];
            let a = digest_fragment(&src, RestrictionCut::blunt(0), c, false, false).unwrap();
            let b = digest_fragment(&src, c, RestrictionCut::blunt(seq.len()), false, false).unwrap();
            assert_eq!(ligate(&[a.clone(), b.clone()], false).unwrap().seq, seq.as_bytes());
            let p = ligate(&[b.flipped(), a.flipped()], false).unwrap();
            assert_eq!(p.seq, revcomp_str(seq.as_bytes()));
            assert_eq!(a.flipped().ends.upstream_stickiness, c.stickiness());
        }
    }

    #[test]
    fn incompatible_sequence_polarity_and_blunt_sticky_joins_are_rejected() {
        let a = Source::new("A", b"GGGGAATT".to_vec(), false).with_ends(EndProperties {
            downstream_stickiness: 4, ..Default::default()
        });
        let b = Source::new("B", b"AATTCCCC".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: 4, ..Default::default()
        });
        assert_eq!(ligate(&[a.clone(), b.clone()], false).unwrap().seq, b"GGGGAATTCCCC");
        let mut bad = b.clone(); bad.seq[0] = b'C';
        assert!(ligate(&[a.clone(), bad], false).is_err());
        let mut bad = b.clone(); bad.ends.upstream_stickiness = -4;
        assert!(ligate(&[a.clone(), bad], false).is_err());
        let mut bad = b.clone(); bad.ends.upstream_stickiness = 0;
        assert!(ligate(&[a, bad], false).is_err());
        assert!(ligate(&[b], true).is_err());
        assert!(ligate(&[Source::new("ring", b"ACGT".to_vec(), true)], false).is_err());
    }

    #[test]
    fn blunting_fills_five_prime_and_removes_three_prime_tails() {
        for (seq, enzyme, expected_left, expected_right) in [
            ("AAAAGAATTCTTTT", "EcoRI", "AAAAGAATT", "AATTCTTTT"),
            ("AAAAGACGTCTTTT", "AatII", "AAAAG", "CTTTT"),
        ] {
            let src = Source::new("linear", seq.as_bytes().to_vec(), false);
            let c = cuts(seq, false, enzyme)[0];
            for (left, right, expected) in [
                (RestrictionCut::blunt(0), c, expected_left),
                (c, RestrictionCut::blunt(seq.len()), expected_right),
            ] {
                let sticky = digest_fragment(&src, left, right, false, false).unwrap();
                let blunt = digest_fragment(&src, left, right, true, false).unwrap();
                assert_eq!(blunt.seq, expected.as_bytes());
                assert_eq!((blunt.ends.upstream_stickiness, blunt.ends.downstream_stickiness), (0, 0));
                assert_eq!(blunt_ends(&sticky).unwrap().seq, blunt.seq);
            }
        }
    }

    #[test]
    fn linear_product_retains_outer_markers_and_native_end_metadata() {
        let a = Source::new("A", b"TGGGG".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: -1, upstream_modification: "Phosphorylated".into(),
            ..Default::default()
        });
        let b = Source::new("B", b"CCCCAA".to_vec(), false).with_ends(EndProperties {
            downstream_stickiness: 2, ..Default::default()
        });
        let p = ligate(&[a, b], false).unwrap();
        assert_eq!(p.seq, b"TGGGGCCCCAA");
        let f = hx_formats::DnaFile::parse(&super::super::product_bytes(&p, "joined")).unwrap();
        assert_eq!(f.end_properties(), p.ends);
        assert!(!p.circular);
    }

    #[test]
    fn circular_equal_boundaries_rotate_and_linear_boundaries_never_wrap() {
        let s = Source::new("ring", b"AAAACCCC".to_vec(), true);
        assert_eq!(s.slice(4, 4, false).seq, b"CCCCAAAA");
        let s = Source::new("line", b"AAAACCCC".to_vec(), false);
        assert!(s.slice(6, 2, false).seq.is_empty());
        assert!(digest_fragment(&s, RestrictionCut::blunt(6), RestrictionCut::blunt(2), false, false).is_err());
    }

    #[test]
    fn features_follow_marker_removal_reverse_and_ligation() {
        use hx_formats::annot::{FeatureAnnotation, FeatureSegment, Range};
        let mut s = Source::new("annotated", b"AAAAGACGTCTTTT".to_vec(), false);
        s.features.push(FeatureAnnotation {
            name: "test".into(), kind: "CDS".into(), sub_type: String::new(), note: String::new(),
            visible: true, directionality: 1, reading_frame: 1,
            translate_first_codon_as_met: false, max_run_on: None,
            segments: vec![FeatureSegment {
                range: Range::new(9, 13), color: "#ff0000".into(), seg_type: "standard".into(), translated: false,
            }],
        });
        let c = cuts("AAAAGACGTCTTTT", false, "AatII")[0];
        let f = digest_fragment(&s, c, RestrictionCut::blunt(s.seq.len()), false, false).unwrap();
        assert_eq!(f.top_range(), (4, 9));
        assert_eq!(f.features[0].segments[0].range, Range::new(4, 8));
        let p = ligate(&[f], false).unwrap();
        assert_eq!(p.features[0].segments[0].range, Range::new(4, 8));
        let blunt = blunt_ends(&Source::new("f", p.seq, false).with_features(p.features).with_ends(p.ends)).unwrap();
        assert_eq!(blunt.features[0].segments[0].range, Range::new(0, 4));
        let reverse = blunt.flipped();
        assert_eq!(reverse.features[0].segments[0].range, Range::new(1, 5));
        assert_eq!(reverse.features[0].directionality, 2);
    }

    #[test]
    fn asymmetric_ends_reverse_without_changing_three_prime_polarity() {
        let s = Source::new("asymmetric", b"TACGTCC".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: -1, downstream_stickiness: 2,
            upstream_modification: "Phosphorylated".into(), ..Default::default()
        });
        let f = s.flipped();
        assert_eq!(f.seq, b"GGACGTA");
        assert_eq!((f.ends.upstream_stickiness, f.ends.downstream_stickiness), (2, -1));
        assert_eq!(f.top_seq(), b"GGACGTA");
        assert_eq!(f.ends.downstream_modification, "Phosphorylated");
        assert_eq!(f.flipped().seq, s.seq);
        assert_eq!(f.flipped().ends, s.ends);
    }

    #[test]
    fn annotations_on_bottom_strand_markers_survive_internal_join() {
        use hx_formats::annot::{FeatureAnnotation, FeatureSegment, Range};
        let mut a = Source::new("left", b"GGGGAATT".to_vec(), false).with_ends(EndProperties {
            downstream_stickiness: 4, ..Default::default()
        });
        a.features.push(FeatureAnnotation {
            name: "sticky annotation".into(), kind: "misc_feature".into(), sub_type: String::new(), note: String::new(),
            visible: true, directionality: 2, reading_frame: 1,
            translate_first_codon_as_met: false, max_run_on: None,
            segments: vec![FeatureSegment {
                range: Range::new(4, 8), color: "#ff0000".into(), seg_type: "standard".into(), translated: false,
            }],
        });
        let b = Source::new("right", b"AATTCCCC".to_vec(), false).with_ends(EndProperties {
            upstream_stickiness: 4, ..Default::default()
        });
        let p = ligate(&[a, b], false).unwrap();
        assert_eq!(p.seq, b"GGGGAATTCCCC");
        assert_eq!(p.features[0].segments[0].range, Range::new(4, 8));
        assert_eq!(p.features[0].directionality, 2);
    }

    #[test]
    fn overlapping_overhangs_cannot_be_rescued_by_blunting() {
        let s = Source::new("short", b"AATTCCGG".to_vec(), false);
        let left = RestrictionCut { top: 0, bottom: 4 };
        let right = RestrictionCut { top: 3, bottom: 7 };
        for blunt in [false, true] {
            assert!(digest_fragment(&s, left, right, blunt, false).is_err());
        }
        let s = s.with_ends(EndProperties {
            upstream_stickiness: 4, downstream_stickiness: -4, ..Default::default()
        });
        assert!(ligate(&[s.clone()], false).is_err());
        assert!(blunt_ends(&s).is_err());
    }

    #[test]
    fn blunting_does_not_silently_remove_terminal_chemical_modifications() {
        let s = Source::new("modified", b"AAAAGAATTCTTTT".to_vec(), false).with_ends(EndProperties {
            upstream_modification: "Phosphorylated".into(), ..Default::default()
        });
        let c = cuts("AAAAGAATTCTTTT", false, "EcoRI")[0];
        let f = digest_fragment(&s, RestrictionCut::blunt(0), c, true, false).unwrap();
        assert_eq!(f.ends.upstream_modification, "Phosphorylated");
        assert_eq!(blunt_ends(&s).unwrap().ends.upstream_modification, "Phosphorylated");
    }

    #[test]
    fn circular_annotation_does_not_duplicate_the_repeated_sticky_region() {
        use hx_formats::annot::{FeatureAnnotation, FeatureSegment, Range};
        let mut s = Source::new("ring", b"TTCGGGAA".to_vec(), true);
        s.features.push(FeatureAnnotation {
            name: "across cut".into(), kind: "misc_feature".into(), sub_type: String::new(), note: String::new(),
            visible: true, directionality: 1, reading_frame: 1,
            translate_first_codon_as_met: false, max_run_on: None,
            segments: vec![FeatureSegment {
                range: Range::new(0, 4), color: "#ff0000".into(), seg_type: "standard".into(), translated: false,
            }],
        });
        let c = cuts("TTCGGGAA", true, "EcoRI")[0];
        let f = digest_fragment(&s, c, c, false, false).unwrap();
        let p = ligate(&[f], true).unwrap();
        assert_eq!(p.features.len(), 1);
        assert_eq!(p.features[0].segments.len(), 1);
        assert_eq!(p.features[0].segments[0].range, Range::new(2, 6));
        let file = hx_formats::DnaFile::parse(&super::super::product_bytes(&p, "circular")).unwrap();
        let features = hx_formats::annot::parse_features(&file.block(10).unwrap().payload);
        assert_eq!(features[0].segments[0].range, Range::new(2, 6));
    }

    #[test]
    fn reaction_checks_every_internal_and_closing_phosphate_pair() {
        let plain = Source::new("blunt", b"ACGT".to_vec(), false);
        let mut a = plain.clone();
        let b = plain.clone();
        assert!(ligate_reaction(&[a.clone(), b.clone()], false).is_err());
        a.ends.downstream_modification = "FivePrimePhosphorylated".into();
        assert!(ligate_reaction(&[a.clone(), b.clone()], false).is_ok());
        assert!(ligate_reaction(&[a.clone(), b.clone()], true).unwrap_err().contains("Fragments 2 and 1"));
        // 把第一个片段反向后磷酸基移到外侧，不再能救助内部接头。
        assert!(ligate_reaction(&[a.flipped(), b.clone()], false).is_err());
        a.ends.upstream_modification = "FivePrimePhosphorylated".into();
        assert!(ligate_reaction(&[a.clone(), b.clone()], true).is_ok());
        assert!(ligate_reaction(&[a, b.clone(), b], false).unwrap_err().contains("Fragments 2 and 3"));
    }

    #[test]
    fn restoring_vector_origin_rotates_features_spans_and_junctions_together() {
        use hx_formats::annot::{FeatureAnnotation, FeatureSegment, Range};
        let v = Source::new("vector", b"AAAAGAATTCCCCCCCCCCCCCCCCCGGATCCGGGG".to_vec(), true);
        let i = Source::new("insert", b"TTGAATTCAAACCCGGATCCAA".to_vec(), false);
        let cut = |s: &Source, name: &str| restriction_cuts(&s.seq, s.circular, hx_enzymes::by_name(name).unwrap())[0];
        let mut backbone = digest_fragment(&v, cut(&v, "BamHI"), cut(&v, "EcoRI"), false, false).unwrap();
        let first = backbone.ends.original_first_base.unwrap();
        backbone.features.push(FeatureAnnotation {
            name: "origin".into(), kind: "misc_feature".into(), sub_type: String::new(), note: String::new(),
            visible: true, directionality: 1, reading_frame: 1,
            translate_first_codon_as_met: false, max_run_on: None,
            segments: vec![FeatureSegment { range: Range::new(first, first + 3), color: "#ff0000".into(), seg_type: "standard".into(), translated: false }],
        });
        let insert = digest_fragment(&i, cut(&i, "EcoRI"), cut(&i, "BamHI"), false, false).unwrap();
        let p = restriction_product(&[backbone, insert]).unwrap();
        assert_eq!(p.seq, b"AAAAGAATTCAAACCCGGATCCGGGG");
        assert_eq!(p.features[0].segments[0].range, Range::new(0, 3));
        assert_eq!(p.junctions.iter().map(|j| j.at).collect::<Vec<_>>(), vec![5, 17]);
        assert_eq!(p.spans.iter().map(|(a, b)| b - a).sum::<usize>(), p.seq.len());
    }
}
