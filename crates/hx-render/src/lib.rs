//! # hx-render — 图谱与序列视图的**纯几何**层
//!
//! 与 egui 解耦：这里只做「碱基坐标 → 屏幕坐标」的换算、特征带状形状、
//! 刻度与标签避让、序列视图排版；绘制原语（线段/多边形/文本）由 `hx-app`
//! 或 SVG/PNG 导出器负责。屏幕坐标系 y 向下（与 egui 一致）。
//!
//! 约定：碱基坐标一律 **0-based 半开**（`[start, end)`），序列长度 `len`。
//! 环状特征可**跨原点**（`start > end`，此时表示 `[start, len) ∪ [0, end)`）。

#![forbid(unsafe_code)]

pub mod map;
pub mod circular_labels;
pub mod rna_layout;
pub mod rna_svg;
pub mod seq;
pub mod svg;

pub use map::{
    assign_label_levels, assign_label_levels_strict, ruler_ticks, ruler_ticks_subdiv, MapLayout,
    MapView, Topology,
};
pub use rna_layout::{RnaLayout, RING_RADIUS, RING_STROKE};
pub use seq::SeqView;

/// 千分位整数（`1234567` → `1,234,567`）；刻度坐标与长度标注统一用它。
/// SnapGene 的位点/长度显示：**只有 5 位及以上（≥ 10000）才加千位分隔**，
/// 4 位数照原样写（oracle：`BglII (5600)`、`SphI (2697)`、标题 `2737 bp`；
/// 而 pAd 写 `(25,803) SpeI`、`34,864 bp`）。
pub fn thousands(n: usize) -> String {
    if n < 10_000 {
        return n.to_string();
    }
    let digits = n.to_string();
    let bytes = digits.as_bytes();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

/// 角度归一化到 `[0, 2π)`。
#[inline]
pub fn norm_angle(a: f32) -> f32 {
    let t = std::f32::consts::TAU;
    let mut r = a % t;
    if r < 0.0 {
        r += t;
    }
    r
}

/// 把（可能跨原点的）区间拆成不跨原点的若干半开区间。
///
/// `start == end` 视为**整圈**（返回 `[(0, len)]`）；`start < end` 原样返回。
pub fn split_wrap(start: usize, end: usize, len: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }
    let start = start.min(len);
    let end = end.min(len);
    if start < end {
        vec![(start, end)]
    } else if start > end {
        vec![(start, len), (0, end)]
    } else if start == 0 {
        // `0..0`：整圈
        vec![(0, len)]
    } else {
        // `start == end > 0`：环状上视为跨原点整圈（与 SnapGene `A > A` 的退化写法一致）
        vec![(start, len), (0, end)]
    }
}

/// 区间长度（环状语义：`start > end` 时＝跨原点长度）。
pub fn span_len(start: usize, end: usize, len: usize) -> usize {
    if start <= end {
        end - start
    } else {
        len - start + end
    }
}

/// 多段特征在 **Features 面板 / 底部状态栏** 里显示为一行时的取值
/// （oracle 实测，SnapGene 8.2.2）。返回 `(start, end, len, wrapped)`：
///
/// - `start` = **首段**起点、`end` = **末段**终点（0-based 半开）：`domain interfaces [active]`
///   的 25 段跨 `3 .. 414`，SnapGene 就写 `3 .. 414`（不是首段的 `3 .. 4`）；
/// - 跨原点照原样写（`CmR` 写 `8886 .. 219`，**不加 `*` 或其他记号**），此时 `end < start`，
///   `wrapped` 只留给提示文本用；
/// - `len` = **各段长度之和**（同一特征 SnapGene 的 Size 列写 `39 aa`，而单段只有 2 aa）。
pub fn feature_span(segs: &[(usize, usize)], seqlen: usize) -> (usize, usize, usize, bool) {
    let start = segs.first().map(|s| s.0).unwrap_or(0);
    let end = segs.last().map(|s| s.1).unwrap_or(0);
    let len = segs.iter().map(|(a, b)| span_len(*a, *b, seqlen)).sum();
    let wrapped = segs.iter().any(|(a, b)| a > b);
    (start, end, len, wrapped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_non_wrapping() {
        assert_eq!(split_wrap(3, 8, 100), vec![(3, 8)]);
        assert_eq!(split_wrap(0, 0, 100), vec![(0, 100)]);
        assert_eq!(split_wrap(5, 5, 100), vec![(5, 100), (0, 5)]);
        assert_eq!(split_wrap(0, 0, 0), Vec::<(usize, usize)>::new());
    }

    #[test]
    fn split_wrapping() {
        assert_eq!(split_wrap(90, 10, 100), vec![(90, 100), (0, 10)]);
        assert_eq!(span_len(90, 10, 100), 20);
        assert_eq!(span_len(10, 90, 100), 80);
    }

    #[test]
    fn angle_normalization() {
        use std::f32::consts::{PI, TAU};
        assert!((norm_angle(-PI / 2.0) - (3.0 * PI / 2.0)).abs() < 1e-5);
        assert!((norm_angle(TAU + 1.0) - 1.0).abs() < 1e-5);
    }

    /// oracle（SnapGene 8.2.2 Features 面板）：Location = 首段起点 .. 末段终点，
    /// Size = 各段长度之和，跨原点不加记号。
    #[test]
    fn feature_span_matches_snapgene_rows() {
        // KPYK1_ECOLI.prot `domain interfaces [active]`：25 段、3 .. 414、39 aa
        let segs = [(2, 4), (59, 60), (61, 62)];
        assert_eq!(feature_span(&segs, 1000), (2, 62, 4, false));
        // 单段（pCas9.dna `Cas9` → 2225 .. 6331 / 4107 bp）
        assert_eq!(feature_span(&[(2224, 6331)], 9326), (2224, 6331, 4107, false));
        // 跨原点（pCas9.dna `CmR` → 8886 .. 219 / 660 bp，无 `*`）
        assert_eq!(feature_span(&[(8885, 219)], 9326), (8885, 219, 660, true));
        // 空特征不 panic
        assert_eq!(feature_span(&[], 100), (0, 0, 0, false));
    }
}
