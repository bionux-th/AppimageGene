//! # hx-enzymes — 酶识别位点匹配
//!
//! 匹配语义与 SnapGene 对齐（由 `docs/SnapGene_dna_format_spec.md` §5 的
//! block 3 缓存全量比对裁定）：
//!
//! - IUPAC 感知（`N`、`R`、`Y`、`S`、`W`、`K`、`M`、`B`、`D`、`H`、`V`）；
//! - **包含语义**：模式字符的碱基集必须覆盖序列位置的全部可能碱基
//!   （序列 `N` 只能匹配模式 `N`；见 [`base_matches`]）；
//! - 允许重叠（逐碱基滑动，不做跳跃）；
//! - 环状序列：起点 ∈ `[0, L)`，包含跨原点命中；线性序列不跨末端截断。
//!
//! 酶数据库（[`data::ENZYMES`]，584 条）由 `research/enzymes-migrate.py` 从旧
//! HelixBench 资产 `assets/enzymes-data.js` 生成；[`sites`] 是 block 3 的规范
//! 站点表（473 条，[`block3`] 编解码用）；本模块提供与数据库无关的底层匹配原语。

#![forbid(unsafe_code)]

pub mod block3;
pub mod data;
pub mod catalog;
pub mod details;
pub mod sets;
pub mod sites;

pub use data::{
    by_name, common_indices, preferred_indices, six_cutter_indices, Enzyme, Overhang, COMMON,
    ENZYMES, PREFERRED,
};

/// SnapGene 显示用的酶位点编号（0-based，直接作为标签里的数字）。
///
/// 由 oracle 实拍裁定（pAd 34,864 bp 的 Enzymes 面板与图谱标签）：数字 =
/// **识别序列起点（0-based）+ 正向链切点偏移 `cut`**，等价于「切点前那个碱基」的
/// 1-based 坐标；`cut == 1` 的酶（EcoRI 等）与「起点 1-based」重合，`cut > 1`
/// 的酶（MreI 2、PmeI 4、AsiSI 5）才能区分开。
///
/// 环状序列越界时按长度回绕；线性序列夹到 `[0, len]`。
pub fn display_pos(start: usize, cut: i16, len: usize, circular: bool) -> usize {
    let p = start as i64 + cut as i64;
    if circular && len > 0 {
        p.rem_euclid(len as i64) as usize
    } else {
        p.clamp(0, len as i64) as usize
    }
}

/// IUPAC 碱基 → 位集合（bit0=A, bit1=C, bit2=G, bit3=T）。未知字符返回 `None`。
#[inline]
pub fn iupac_bits(base: u8) -> Option<u8> {
    Some(match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,             // A|G
        b'Y' => 0b1010,             // C|T
        b'S' => 0b0110,             // C|G
        b'W' => 0b1001,             // A|T
        b'K' => 0b1100,             // G|T
        b'M' => 0b0011,             // A|C
        b'B' => 0b1110,             // C|G|T
        b'D' => 0b1101,             // A|G|T
        b'H' => 0b1011,             // A|C|T
        b'V' => 0b0111,             // A|C|G
        b'N' => 0b1111,
        _ => return None,
    })
}

/// 单个位点字符是否匹配序列字符（两者均按 IUPAC 解释）。
///
/// 采用**包含语义**：仅当模式字符的可能碱基集**覆盖**序列字符的全部可能碱基时
/// 才算命中（`seq_set ⊆ pattern_set`）。对具体碱基（A/C/G/T）退化为通常的
/// “模式包含该碱基”；对序列中的歧义码则更严格，例如序列 `N` 只能匹配模式 `N`。
///
/// 该语义由 203 个自带样例的 block 3 缓存逐字节比对裁定（N/D 碱基样例，
/// 见 `tests/block3_rebuild.rs`）：宽松的“交集非空”语义会多出 7–8 处命中。
#[inline]
pub fn base_matches(pattern: u8, seq: u8) -> bool {
    match (iupac_bits(pattern), iupac_bits(seq)) {
        (Some(p), Some(s)) => p & s == s,
        _ => false,
    }
}

/// 在序列中查找位点的 0-based 起点（升序，允许重叠）。
///
/// `circular = true` 时按“起点 ∈ [0, L)”计跨原点命中（与 SnapGene
/// block 3 计数/位置流一致）。
pub fn find_site_positions(seq: &[u8], site: &[u8], circular: bool) -> Vec<u32> {
    let l = seq.len();
    let m = site.len();
    if m == 0 || l < m {
        return Vec::new();
    }
    let mut out = Vec::new();
    if circular {
        // 环状：在 seq+seq 上匹配长度为 m 的窗口，起点 < L 者计入
        let doubled_len = 2 * l;
        for start in 0..l {
            let mut ok = true;
            for k in 0..m {
                let idx = (start + k) % doubled_len;
                let s = seq[idx % l];
                if !base_matches(site[k], s) {
                    ok = false;
                    break;
                }
            }
            if ok {
                out.push(start as u32);
            }
        }
    } else {
        if m > l {
            return out;
        }
        for start in 0..=(l - m) {
            let mut ok = true;
            for k in 0..m {
                if !base_matches(site[k], seq[start + k]) {
                    ok = false;
                    break;
                }
            }
            if ok {
                out.push(start as u32);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_overlapping() {
        assert_eq!(find_site_positions(b"GATCGGATC", b"GATC", false), vec![0, 5]);
        assert_eq!(find_site_positions(b"AAAA", b"AA", false), vec![0, 1, 2]);
    }

    #[test]
    fn circular_wraps_origin() {
        // ATAT 环状上的 TA：1 与 3（3 为跨原点命中）
        assert_eq!(find_site_positions(b"ATAT", b"TA", false), vec![1]);
        assert_eq!(find_site_positions(b"ATAT", b"TA", true), vec![1, 3]);
    }

    #[test]
    fn iupac_aware() {
        // N 位点：任意碱基；与 block 3 验证一致的简单样例
        assert_eq!(find_site_positions(b"ACGT", b"NNNN", false), vec![0]);
        // ASST = A S S T（S = C|G）：AACGT 的 ACGT 段在偏移 1 命中
        assert_eq!(find_site_positions(b"ATCTT", b"ASST", false), vec![]);
        assert_eq!(find_site_positions(b"AACGT", b"ASST", false), vec![1]);
    }

    #[test]
    fn ambiguous_bases_require_subset() {
        // 序列 N（未知碱基）只能匹配模式 N —— 由 203 个样例的 block 3 缓存裁定
        assert!(!base_matches(b'G', b'N'));
        assert!(!base_matches(b'R', b'N'));
        assert!(base_matches(b'N', b'N'));
        // 序列 D = A|G|T：模式须覆盖 {A,G,T}
        assert!(!base_matches(b'H', b'D')); // H = A|C|T，缺 G
        assert!(!base_matches(b'B', b'D')); // B = C|G|T，缺 A
        assert!(base_matches(b'D', b'D'));
        assert!(base_matches(b'N', b'D'));
        // 具体碱基仍按“模式包含该碱基”
        assert!(base_matches(b'R', b'A'));
        assert!(!base_matches(b'Y', b'A'));
        // 端到端：GAATTC 中含一个 N 的位置不计命中
        assert_eq!(find_site_positions(b"GAATTC", b"GAAYTC", false), vec![0]);
        assert_eq!(find_site_positions(b"GAANTC", b"GAAYTC", false), vec![]);
        assert_eq!(find_site_positions(b"GAANTC", b"GAANTC", false), vec![0]);
    }

    #[test]
    fn no_match_when_site_longer() {
        assert!(find_site_positions(b"ACGT", b"ACGTA", false).is_empty());
        assert!(find_site_positions(b"", b"ACG", true).is_empty());
    }
}
