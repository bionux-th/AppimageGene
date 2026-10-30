//! 引物「尾巴」（不与模板配对的额外碱基）落在 `sequence` 的哪一端。
//!
//! SnapGene 的 `<Primer sequence="…">` 永远写成引物 5′→3′，而
//! `<BindingSite annealedBases="…">` 只给退火段的**长度**，没有说明它是
//! `sequence` 的前缀还是后缀 —— 两种写法在真文件里都存在：
//!
//! * 常见的 5′ 尾巴（如 `pET28b-SprB-CDT368.dna` 的反向引物
//!   `GCTCTA…CAGGGT`，退火段是后缀）——尾巴在 `sequence` 前缀；
//! * 在序列框**末尾**多打的碱基（如 `pET28b-Rep1.dna` 的 Primer 2
//!   `CTAAAG…ACGCCC`，退火段是前缀）——尾巴在 `sequence` 后缀，
//!   在序列视图里画在结合框的**左侧**（反向引物的 3′ 端）。
//!
//! 判据只能是模板：退火段必须等于结合区间上按引物方向读出的序列。

use hx_formats::annot::{PrimerSite, Range};

/// 结合区间上的模板碱基，按**引物方向**给出（`strand=1` 时反向互补）。
///
/// 环状文件里 `range.start > range.end` 表示跨原点。
pub fn site_bases(template: &[u8], range: Range, strand: u8) -> Vec<u8> {
    let n = template.len();
    if n == 0 {
        return Vec::new();
    }
    let start = range.start.min(n - 1);
    let count = if range.end >= range.start {
        range.end - range.start
    } else {
        n - range.start + range.end
    };
    let mut out: Vec<u8> = (0..count.min(n))
        .map(|k| template[(start + k) % n].to_ascii_uppercase())
        .collect();
    if strand == 1 {
        out = hx_algo::primer::revcomp(&out);
    }
    out
}

/// `true` = 退火段是 `sequence` 的后缀，也就是常见的 **5′ 尾巴**写法。
///
/// 判不出来时按 5′ 尾巴处理（SnapGene 里绝大多数引物都是这种）。
pub fn annealed_is_suffix(sequence: &str, site: &PrimerSite, template: &[u8]) -> bool {
    let seq: Vec<u8> = sequence.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let ann = site.annealed_bases.min(seq.len());
    if ann == 0 || ann == seq.len() {
        return true;
    }
    let want = site_bases(template, site.range, site.strand);
    if want.len() < ann {
        return true;
    }
    let w = &want[..ann];
    let pre = seq[..ann].eq_ignore_ascii_case(w);
    let suf = seq[seq.len() - ann..].eq_ignore_ascii_case(w);
    match (pre, suf) {
        (true, false) => false,
        (false, true) => true,
        _ => {
            // A permitted internal mismatch must not turn a prefix-bound site
            // into a suffix-bound site and move its unpaired tail to 5′.
            let mismatches=|part:&[u8]|part.iter().zip(w).filter(|(a,b)|!a.eq_ignore_ascii_case(b)).count();
            let pre_mm=mismatches(&seq[..ann]);let suf_mm=mismatches(&seq[seq.len()-ann..]);
            !(pre_mm<suf_mm)
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(range: Range, strand: u8, annealed: usize) -> PrimerSite {
        PrimerSite {
            range,
            strand,
            annealed_bases: annealed,
            melting_temperature: None,
        }
    }

    #[test]
    fn five_prime_tail_is_prefix() {
        // pET28b 上的正向引物片段：退火段 = 模板，尾巴在 5′ 端。
        let template = b"TTTTTTGGGACGCGCCCTGTAGCGGCGCATTAAGCTTTTTT";
        let s = site(Range::new(6, 33), 0, 27);
        assert!(annealed_is_suffix(
            "TGGCGAATGGGACGCGCCCTGTAGCGGCGCATTAAGC",
            &s,
            template
        ));
    }

    #[test]
    fn three_prime_tail_is_suffix() {
        // pET28b-Rep1.dna 的 Primer 2：反向引物，多出来的 CCC 在 sequence 末尾。
        let template = b"AAAACGTCAAGCTCTAAATCGGGGGCTCCCTTTAGAAAA";
        let s = site(Range::new(4, 35), 1, 31);
        assert!(!annealed_is_suffix(
            "CTAAAGGGAGCCCCCGATTTAGAGCTTGACGCCC",
            &s,
            template
        ));
    }
}
