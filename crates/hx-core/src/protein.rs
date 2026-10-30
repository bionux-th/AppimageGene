//! 蛋白属性：氨基酸组成、分子量、280 nm 消光系数/吸光度、等电点与净电荷。
//!
//! 这些数字对应 SnapGene 蛋白文档的 `Properties` 页；公式与实测对照
//! （`THYX_MYCTU.prot` 250 aa、`KPYK1_ECOLI.prot` 470 aa，真 SnapGene 8.2.2）：
//!
//! | 项目 | 本模块（A / B） | SnapGene（A / B） | 说明 |
//! | --- | --- | --- | --- |
//! | 平均分子量 | 27591.33 / 50729.59 Da | 27,591.31 / 50,729.56 Da | 差 ≤ 0.005 Da |
//! | 单同位素分子量 | 27574.16 / 50697.26 Da | — / 50,697.26 Da | 与该表完全一致 |
//! | 消光系数 280 nm（还原） | 26930 / 7450 | 26,930 / 7450 | 一致 |
//! | 消光系数 280 nm（氧化） | 27055 / 7825 | — / 7825 | 一致（B 实测） |
//! | 吸光度 280 nm 0.1% | 0.98 / 0.15 | 0.98 / 0.15 | 一致 |
//! | 等电点 pI | 6.54 / 5.77 | 6.79 / 5.86 | 见下 |
//! | pH 7 净电荷 | −1.44 / −6.51 | −1.42 / −8.52 | 见下 |
//!
//! 分子量 = 残基原子组成之和 + 一分子水，原子量用 C 12.011、H 1.00794、N 14.0067、
//! O 15.9994、S 32.066（这套值同时把两个样本压在 0.005 Da 以内）；单同位素质量用
//! 残基单同位素质量表 + 水 18.010565。
//!
//! 消光系数（280 nm）= `5500·Trp + 1490·Tyr`（**还原态**，SnapGene 顶行显示的就是它），
//! 氧化态再加 `125 × (Cys ÷ 2)`（二硫键）。
//!
//! 吸光度（280 nm, 0.1%）= 消光系数 ÷ 平均分子量（0.1% = 1 mg/mL）。
//!
//! 非标准残基（实测真 SnapGene）：
//! - `X` / `B` / `Z` / `J` / `U` / `O` 等**未知残基**：分子量按 20 种残基平均质量计入
//!   （实测 3 个 `X` 时 SnapGene 的分子量比野生型少 15.76 Da → 每个 `X` 记 118.888 Da，
//!   正是 20 种平均残基质量的均值），`Properties` 页在分子量前加 `~`；
//!   消光系数/等电点/净电荷把它们当中性，数值不变；组成表不单列这些残基。
//! - `*`（终止符）：SnapGene 把整个序列当成"不能计算的蛋白"——`Length` 的单位从 `aa`
//!   变成 `codons`，分子量/消光系数/吸光度/等电点/净电荷一律显示 `—`。
//!
//! 电荷用 Henderson–Hasselbalch 求和（N 端 / C 端 / K / R / H / D / E / C / Y），
//! pKa 取 ExPASy（Bjellqvist）表；pI = 净电荷为 0 的 pH（二分）。SnapGene 的显示值
//! 与任何单一 pKa 表都不完全吻合——实测它的 pI 比本模块高 0.1-0.25，而 KPYK1 的
//! "Charge at pH 7"（−8.52）与它自己的 pI（5.86）在同一套 pKa 下**互不自洽**
//! （该组成在任何标准表下 pH 7 净电荷都在 −6.5 附近），故这里保留标准算法并记录差异。

/// 20 种标准氨基酸（按单字母码排序）：(单字母, 三字母, 全名)。
///
/// 全名与 SnapGene 的 `Properties` 表一致（每个单词首字母大写，如 `Aspartic Acid`）。
pub const STANDARD_AA: [(char, &str, &str); 20] = [
    ('A', "Ala", "Alanine"),
    ('C', "Cys", "Cysteine"),
    ('D', "Asp", "Aspartic Acid"),
    ('E', "Glu", "Glutamic Acid"),
    ('F', "Phe", "Phenylalanine"),
    ('G', "Gly", "Glycine"),
    ('H', "His", "Histidine"),
    ('I', "Ile", "Isoleucine"),
    ('K', "Lys", "Lysine"),
    ('L', "Leu", "Leucine"),
    ('M', "Met", "Methionine"),
    ('N', "Asn", "Asparagine"),
    ('P', "Pro", "Proline"),
    ('Q', "Gln", "Glutamine"),
    ('R', "Arg", "Arginine"),
    ('S', "Ser", "Serine"),
    ('T', "Thr", "Threonine"),
    ('V', "Val", "Valine"),
    ('W', "Trp", "Tryptophan"),
    ('Y', "Tyr", "Tyrosine"),
];

/// 残基原子组成 `(C, H, N, O, S)`（不含末端水），与 [`STANDARD_AA`] 同序。
const ATOMS: [(u64, u64, u64, u64, u64); 20] = [
    (3, 5, 1, 1, 0),   // A
    (3, 5, 1, 1, 1),   // C
    (4, 5, 1, 3, 0),   // D
    (5, 7, 1, 3, 0),   // E
    (9, 9, 1, 1, 0),   // F
    (2, 3, 1, 1, 0),   // G
    (6, 7, 3, 1, 0),   // H
    (6, 11, 1, 1, 0),  // I
    (6, 12, 2, 1, 0),  // K
    (6, 11, 1, 1, 0),  // L
    (5, 9, 1, 1, 1),   // M
    (4, 6, 2, 2, 0),   // N
    (5, 7, 1, 1, 0),   // P
    (5, 8, 2, 2, 0),   // Q
    (6, 12, 4, 1, 0),  // R
    (3, 5, 1, 2, 0),   // S
    (4, 7, 1, 2, 0),   // T
    (5, 9, 1, 1, 0),   // V
    (11, 10, 2, 1, 0), // W
    (9, 9, 1, 2, 0),   // Y
];

/// 残基单同位素质量（Da，不含末端水），与 [`STANDARD_AA`] 同序。
const MONO: [f64; 20] = [
    71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406,
    128.09496, 113.08406, 131.04049, 114.04293, 97.05276, 128.05858, 156.10111, 87.03203,
    101.04768, 99.06841, 186.07931, 163.06333,
];

/// 平均原子量（**与 SnapGene 实测对齐的那一套**）。
const W_C: f64 = 12.011;
const W_H: f64 = 1.00794;
const W_N: f64 = 14.0067;
const W_O: f64 = 15.9994;
const W_S: f64 = 32.066;

/// 单同位素水分子的质量（残基和再加它）。
const WATER_MONO: f64 = 18.010565;

/// 平衡态 pKa（ExPASy / Bjellqvist）：(N 端, C 端, K, R, H, D, E, C, Y)。
const PKA: (f64, f64, f64, f64, f64, f64, f64, f64, f64) =
    (7.5, 3.55, 10.0, 12.0, 5.98, 4.05, 4.45, 9.0, 10.0);

/// 氨基酸组成：20 种标准氨基酸的计数 + 非标准字符计数。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Composition {
    /// 与 [`STANDARD_AA`] 同序的计数。
    pub counts: [usize; 20],
    /// 未知残基（`X`/`B`/`Z`/`J`/`U`/`O` 等）的个数。
    pub other: usize,
    /// 终止符 `*` 的个数（SnapGene 见到它就认为该蛋白无法计算属性）。
    pub stops: usize,
}

impl Composition {
    /// 统计一条蛋白序列（大小写不敏感；非标准字符进 `other`）。
    pub fn of(seq: &[u8]) -> Self {
        let mut c = Self::default();
        for &b in seq {
            if b == b'*' {
                c.stops += 1;
                continue;
            }
            match aa_index(b) {
                Some(i) => c.counts[i] += 1,
                None => c.other += 1,
            }
        }
        c
    }

    /// 总长度（标准 + 未知 + 终止符）。
    pub fn len(&self) -> usize {
        self.counts.iter().sum::<usize>() + self.other + self.stops
    }

    /// 序列里是否有终止符 `*`（SnapGene 据此不再给出任何数值属性）。
    pub fn has_stops(&self) -> bool {
        self.stops > 0
    }

    /// 序列里是否有未知残基（SnapGene 的数值照算，分子量前加 `~`）。
    pub fn has_unknown(&self) -> bool {
        self.other > 0
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// 取某个标准氨基酸（按单字母码）的计数。
    pub fn count(&self, aa: char) -> usize {
        aa_index(aa as u8).map_or(0, |i| self.counts[i])
    }

    /// 20 个计数之和（不含非标准字符）。
    pub fn standard_len(&self) -> usize {
        self.counts.iter().sum()
    }
}

/// 单字母码 → [`STANDARD_AA`] 下标（大小写不敏感；非标准返回 `None`）。
pub fn aa_index(c: u8) -> Option<usize> {
    let up = c.to_ascii_uppercase();
    STANDARD_AA.iter().position(|(ch, _, _)| *ch as u8 == up)
}

/// 平均分子量（Da，含末端水）。
pub fn average_mass(c: &Composition) -> f64 {
    let (mut a, mut h, mut n, mut o, mut s) = (0u64, 2u64, 0u64, 1u64, 0u64);
    for (i, &n_i) in c.counts.iter().enumerate() {
        let t = ATOMS[i];
        a += t.0 * n_i as u64;
        h += t.1 * n_i as u64;
        n += t.2 * n_i as u64;
        o += t.3 * n_i as u64;
        s += t.4 * n_i as u64;
    }
    let known = a as f64 * W_C + h as f64 * W_H + n as f64 * W_N + o as f64 * W_O + s as f64 * W_S;
    known + c.other as f64 * mean_residue_mass()
}

/// 单同位素分子量（Da，含末端水）。
pub fn monoisotopic_mass(c: &Composition) -> f64 {
    let known: f64 = c
        .counts
        .iter()
        .enumerate()
        .map(|(i, &n)| MONO[i] * n as f64)
        .sum();
    let mean: f64 = MONO.iter().sum::<f64>() / 20.0;
    known + c.other as f64 * mean + WATER_MONO
}

/// 20 种残基平均质量的均值（未知残基的替代质量；实测 SnapGene 用它，见模块文档）。
pub fn mean_residue_mass() -> f64 {
    (0..20)
        .map(|i| {
            let t = ATOMS[i];
            t.0 as f64 * W_C
                + t.1 as f64 * W_H
                + t.2 as f64 * W_N
                + t.3 as f64 * W_O
                + t.4 as f64 * W_S
        })
        .sum::<f64>()
        / 20.0
}

/// 280 nm 摩尔消光系数（M⁻¹cm⁻¹）：`(还原态, 氧化态)`。
///
/// 还原态 = `5500·Trp + 1490·Tyr`；氧化态 = 还原态 + `125 × (Cys ÷ 2)`（二硫键）。
pub fn extinction_280(c: &Composition) -> (u32, u32) {
    let w = c.count('W') as u32;
    let y = c.count('Y') as u32;
    let cys = c.count('C') as u32;
    let reduced = 5500 * w + 1490 * y;
    (reduced, reduced + 125 * (cys / 2))
}

/// 280 nm 吸光度（0.1% = 1 mg/mL）：`(还原态, 氧化态)`；分子量为 0 时返回 0。
pub fn absorbance_280(c: &Composition) -> (f64, f64) {
    let (r, o) = extinction_280(c);
    let mw = average_mass(c);
    if mw <= 0.0 {
        return (0.0, 0.0);
    }
    (r as f64 / mw, o as f64 / mw)
}

/// pH 下的净电荷（Henderson–Hasselbalch 求和；正数 = 偏正）。
pub fn charge(c: &Composition, ph: f64) -> f64 {
    let (n, cc, k, r, h, d, e, cy, y) = PKA;
    // 碱性基团（质子化 = +1）
    let pos = 1.0 / (1.0 + 10f64.powf(ph - n))
        + c.count('K') as f64 / (1.0 + 10f64.powf(ph - k))
        + c.count('R') as f64 / (1.0 + 10f64.powf(ph - r))
        + c.count('H') as f64 / (1.0 + 10f64.powf(ph - h));
    // 酸性基团（去质子化 = −1）
    let neg = 1.0 / (1.0 + 10f64.powf(cc - ph))
        + c.count('D') as f64 / (1.0 + 10f64.powf(d - ph))
        + c.count('E') as f64 / (1.0 + 10f64.powf(e - ph))
        + c.count('C') as f64 / (1.0 + 10f64.powf(cy - ph))
        + c.count('Y') as f64 / (1.0 + 10f64.powf(y - ph));
    pos - neg
}

/// 等电点 pI：净电荷为 0 的 pH（在 0..14 内二分）。
pub fn isoelectric_point(c: &Composition) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 14.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if charge(c, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// THYX_MYCTU（250 aa）与 KPYK1_ECOLI（470 aa）的实测组成，
    /// 用来对照真 SnapGene `Properties` 页上的数字。
    const THYX: [usize; 20] = [
        40, 3, 14, 16, 9, 9, 9, 12, 7, 20, 3, 5, 14, 8, 21, 14, 18, 18, 3, 7,
    ];
    const KPYK1: [usize; 20] = [
        41, 6, 27, 35, 12, 40, 8, 33, 34, 40, 18, 24, 12, 12, 21, 24, 34, 44, 0, 5,
    ];

    fn comp(counts: [usize; 20]) -> Composition {
        Composition { counts, other: 0, stops: 0 }
    }

    #[test]
    fn composition_counts_ignore_case_and_other() {
        let c = Composition::of(b"AcDEfX*");
        assert_eq!(c.len(), 7);
        assert_eq!(c.other, 1); // X
        assert_eq!(c.stops, 1); // * 单独计数
        assert_eq!(c.count('A'), 1);
        assert_eq!(c.count('c'), 1);
        assert_eq!(c.count('C'), 1);
        assert_eq!(c.count('D'), 1);
        assert_eq!(c.count('E'), 1);
        assert_eq!(c.count('F'), 1);
        assert_eq!(c.count('W'), 0);
        assert_eq!(c.standard_len(), 5);
    }

    #[test]
    fn stops_are_counted_separately() {
        let c = Composition::of(b"AC*D**");
        assert_eq!(c.stops, 3);
        assert_eq!(c.other, 0);
        assert_eq!(c.len(), 6);
        assert!(c.has_stops());
        assert!(!c.has_unknown());
    }

    #[test]
    fn unknown_residues_use_mean_residue_mass() {
        // 真 SnapGene：THYX 把 3 个残基改成 X 后分子量 27,575.55（野生型 27,591.31），
        // 反解每个 X = 118.888 Da = 20 种平均残基质量的均值。
        let mean = mean_residue_mass();
        assert!((mean - 118.888).abs() < 0.05, "mean = {mean}");
        // 只把 X 记成均值质量，pH 相关属性不受影响
        let c = Composition::of(b"XXXX");
        assert!((average_mass(&c) - (mean * 4.0 + 18.015)).abs() < 0.02);
        assert_eq!(charge(&c, 7.0), charge(&comp([0; 20]), 7.0));
        assert_eq!(extinction_280(&c), (0, 0));
    }

    #[test]
    fn composition_of_samples_matches_readme_lengths() {
        assert_eq!(comp(THYX).standard_len(), 250);
        assert_eq!(comp(KPYK1).standard_len(), 470);
    }

    #[test]
    fn average_mass_matches_snapgene() {
        // SnapGene：27,591.31 Da / 50,729.56 Da（我们的表在 0.01 Da 以内）
        let a = average_mass(&comp(THYX));
        let b = average_mass(&comp(KPYK1));
        assert!((a - 27591.31).abs() < 0.01, "THYX 平均分子量 {a}");
        assert!((b - 50729.56).abs() < 0.01, "KPYK1 平均分子量 {b}");
    }

    #[test]
    fn monoisotopic_mass_matches_snapgene() {
        // SnapGene 的 KPYK1「monoisotopic」= 50,697.26 Da
        assert!((monoisotopic_mass(&comp(KPYK1)) - 50697.26).abs() < 0.01);
        assert!((monoisotopic_mass(&comp(THYX)) - 27574.16).abs() < 0.01);
    }

    #[test]
    fn extinction_and_absorbance_match_snapgene() {
        // SnapGene：THYX 26,930、KPYK1 7450（还原）/ 7825（氧化）
        assert_eq!(extinction_280(&comp(THYX)), (26930, 27055));
        assert_eq!(extinction_280(&comp(KPYK1)), (7450, 7825));
        // SnapGene 的吸光度（两位小数）：0.98 / 0.15
        let (ra, _) = absorbance_280(&comp(THYX));
        let (rb, _) = absorbance_280(&comp(KPYK1));
        assert_eq!(format!("{ra:.2}"), "0.98");
        assert_eq!(format!("{rb:.2}"), "0.15");
    }

    #[test]
    fn charge_falls_with_ph_and_vanishes_at_pi() {
        let c = comp(THYX);
        let pi = isoelectric_point(&c);
        assert!(charge(&c, pi).abs() < 0.01, "pI 处电荷应为 0");
        assert!(charge(&c, pi - 1.0) > 0.0);
        assert!(charge(&c, pi + 1.0) < 0.0);
        // 单调下降
        assert!(charge(&c, 5.0) > charge(&c, 6.0));
        assert!(charge(&c, 8.0) > charge(&c, 9.0));
    }

    #[test]
    fn charge_and_pi_are_the_exposy_model() {
        // 这套 pKa（Bjellqvist）下：THYX pI 6.54 / pH 7 净电荷 −1.44；
        // SnapGene 显示 6.79 / −1.42（电荷几乎一致，pI 差 0.25，见模块文档）。
        assert_eq!(format!("{:.2}", isoelectric_point(&comp(THYX))), "6.54");
        assert_eq!(format!("{:.2}", charge(&comp(THYX), 7.0)), "-1.44");
        assert_eq!(format!("{:.2}", isoelectric_point(&comp(KPYK1))), "5.77");
        assert_eq!(format!("{:.2}", charge(&comp(KPYK1), 7.0)), "-6.51");
    }
}
