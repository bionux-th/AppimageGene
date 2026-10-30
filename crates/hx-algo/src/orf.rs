//! 开放阅读框架（ORF）查找：6 个阅读框里从 `ATG` 到终止密码子的区段。
//!
//! 语义对齐 SnapGene 的两条显示规则（其翻译库 `SequenceView` 上下文）：
//! - “ORFs are allowed to begin at DNA ends that lack start codons” ⇒ **线性**序列里，
//!   某个阅读框从位置 0 起、在长度阈值内没有起始密码子但有终止密码子，也算一个 ORF；
//! - 环状分子允许 ORF **跨原点**（起点 ∈ `[0, L)`，终点可回绕）。
//!
//! 只做显示用途：不处理可变起始密码子、不解析内含子；同一阅读框内互相嵌套的 ORF
//! 只保留最外层（否则环状图谱上会叠成一团）。

/// 显示阈值为 75 个密码子（225 bp）；有真实终止密码子时将其计入筛选长度。
/// `Orf` 的坐标和 `aa_len` 仍只描述编码区，不把终止密码子计为氨基酸。
pub const DEFAULT_MIN_BP: usize = 75 * 3;

/// 一个 ORF（0-based 半开区间；`end < start` 表示跨原点）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Orf {
    pub start: usize,
    /// 终点（不含终止密码子）。
    pub end: usize,
    /// 1 = 上链（正向），-1 = 下链（反向）。
    pub strand: i8,
    /// 是否由起始密码子 `ATG` 起（否则是“缺起始密码子的端部 ORF”）。
    pub has_start: bool,
}

impl Orf {
    /// 氨基酸个数（不含终止密码子）。
    pub fn aa_len(&self, len: usize) -> usize {
        let span = if self.end >= self.start {
            self.end - self.start
        } else {
            len - self.start + self.end
        };
        span / 3
    }
}

#[inline]
fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

/// 在编码序列上按 3 碱基步长找 ORF，返回 `(start, end, has_start)`（编码坐标，bp）。
#[cfg(test)]
fn orfs_in_frame(code: &[u8], min_bp: usize, allow_end_orf: bool) -> Vec<(usize, usize, bool)> {
    orfs_in_frame_options(code,min_bp,allow_end_orf,&[ *b"ATG" ],true)
}
fn orfs_in_frame_options(code:&[u8],min_bp:usize,allow_end_orf:bool,starts:&[[u8;3]],require_start:bool)->Vec<(usize,usize,bool)>{
    let mut out = Vec::new();
    let n = code.len() / 3;
    if n == 0 {
        return out;
    }
    // 端部 ORF（仅线性序列）：从位置 0 到第一个终止，且其间无 ATG
    if allow_end_orf {
        if let Some(stop) = (0..n).find(|i| is_stop(&code[i * 3..i * 3 + 3])) {
            let has_atg = (0..=stop).any(|i| starts.iter().any(|s|code[i*3..i*3+3]==*s));
            if stop > 0 && !has_atg && (stop + 1) * 3 >= min_bp {
                out.push((0, stop * 3, false));
            }
        }
    }
    for k in 0..n {
        let has_start=starts.iter().any(|s|code[k*3..k*3+3]==*s);
        if if require_start{!has_start}else{k>0&&!is_stop(&code[(k-1)*3..k*3])} {
            continue;
        }
        let mut stop = n;
        for j in k + 1..n {
            if is_stop(&code[j * 3..j * 3 + 3]) {
                stop = j;
                break;
            }
        }
        let (s, e) = (k * 3, stop * 3);
        // SnapGene 的最短 ORF 长度包括真实末端终止密码子。
        // 若阅读框延伸至序列末端，则没有这额外的 3 bp。
        let screening_bp = e - s + if stop < n { 3 } else { 0 };
        if e > s && screening_bp >= min_bp {
            out.push((s, e, has_start));
        }
    }
    // 同一框架内去掉被别的 ORF 完全包含者
    let mut keep: Vec<(usize, usize, bool)> = Vec::new();
    for (s, e, has) in out {
        let nested = keep
            .iter()
            .any(|(s2, e2, _)| *s2 <= s && e <= *e2 && (*s2, *e2) != (s, e));
        if !nested {
            keep.retain(|(s2, e2, _)| !(s <= *s2 && *e2 <= e));
            keep.push((s, e, has));
        }
    }
    keep
}

/// 查找 6 个阅读框的 ORF。`seq` 为上链序列（大小写不敏感）；`circular` = 允许跨原点。
///
/// `min_bp` 为最短区段长度（包含真实末端终止密码子）；返回坐标不含终止密码子。
///
/// 返回按 `(起点, 终点, 链)` 排序的列表。
pub fn find_orfs(seq: &[u8], circular: bool, min_bp: usize) -> Vec<Orf> {
    find_orfs_with_options(seq,circular,min_bp,&[*b"ATG"],true,true)
}
pub fn find_orfs_with_options(seq:&[u8],circular:bool,min_bp:usize,starts:&[[u8;3]],require_start:bool,allow_ends:bool)->Vec<Orf>{
    let len = seq.len();
    if len < 6 {
        return Vec::new();
    }
    let min_bp = min_bp.max(3);
    let l2 = if circular { 2 * len } else { len };
    let fwd: Vec<u8> = (0..l2).map(|i| {let b=seq[i % len].to_ascii_uppercase();if b==b'U'{b'T'}else{b}}).collect();
    let rev: Vec<u8> = fwd
        .iter()
        .rev()
        .map(|b| crate::edit_session::complement(*b))
        .collect();
    let mut out: Vec<Orf> = Vec::new();
    for strand in [1i8, -1] {
        let code = if strand == 1 { &fwd } else { &rev };
        for off in 0..3usize {
            for (cs, ce, has_start) in orfs_in_frame_options(&code[off..], min_bp, !circular&&allow_ends,starts,require_start) {
                let (s, e) = (off + cs, off + ce);
                // 环状命中只从第一圈起算，且不能把超过一圈的区间取模成短 ORF。
                if circular && (s >= len || e - s > len) {
                    continue;
                }
                let (gs, ge) = if strand == 1 {
                    (s, e)
                } else {
                    (l2 - e, l2 - s)
                };
                if circular {
                    let ge = ge % len;
                    let gs = gs % len;
                    // 环状：只保留起点在 [0, L)、且长度不超过整圈的命中
                    let span = if ge >= gs { ge - gs } else { len - gs + ge };
                    if span == 0 || span > len {
                        continue;
                    }
                    out.push(Orf {
                        start: gs,
                        end: ge,
                        strand,
                        has_start,
                    });
                } else {
                    if gs >= len || ge > len {
                        continue;
                    }
                    out.push(Orf {
                        start: gs,
                        end: ge,
                        strand,
                        has_start,
                    });
                }
            }
        }
    }
    out.sort_by_key(|o| (o.start, o.end, o.strand));
    out.dedup_by_key(|o| (o.start, o.end, o.strand));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_simple_orf() {
        // 3 + ATG GCT GCT TTC GCT TAA
        let seq = b"CCCATGGCTGCTTTCGCTTAA".to_vec();
        let orfs = find_orfs(&seq, false, 12);
        let o = orfs.iter().find(|o| o.has_start && o.strand == 1).unwrap();
        assert_eq!((o.start, o.end), (3, 18));
        assert_eq!(o.aa_len(seq.len()), 5);
    }

    #[test]
    fn respects_min_length() {
        let seq = b"CCCATGGCTGCTTTCGCTTAA".to_vec();
        assert!(find_orfs(&seq, false, 19).is_empty());
    }

    #[test]
    fn reverse_strand_orf() {
        // 上链插入 ATGGCTGCTTAA 的反向互补 = TTAAGCAGCCAT（终止子在上链的低位一侧）
        let seq = b"GGGTTAAGCAGCCATCCC".to_vec();
        let orfs = find_orfs(&seq, false, 9);
        let o = orfs
            .iter()
            .find(|o| o.strand == -1 && o.has_start)
            .expect("reverse ORF");
        // 编码 ATG GCT GCT 后的终止子位于基因组 3..6，故 ORF 本体 = 6..15
        assert_eq!((o.start, o.end), (6, 15));
        assert_eq!(o.aa_len(seq.len()), 3);
    }

    #[test]
    fn end_orf_without_start_codon() {
        let seq = b"GCTGCTTTCGCTTAA".to_vec();
        let orfs = find_orfs(&seq, false, 9);
        assert_eq!(orfs.len(), 1);
        assert!(!orfs[0].has_start);
        assert_eq!((orfs[0].start, orfs[0].end), (0, 12));
    }

    #[test]
    fn circular_wraps_origin() {
        // 环状：终止子在原点右侧（TAA）、ATG 在中段 ⇒ ORF 跨原点（end < start）
        let seq = b"TAAATGGGGGGGGGGGGG".to_vec();
        assert_eq!(seq.len(), 18);
        let orfs = find_orfs(&seq, true, 12);
        let o = orfs
            .iter()
            .find(|o| o.strand == 1 && o.start == 3)
            .expect("wrap ORF");
        assert_eq!(o.end, 0, "{orfs:?}");
        assert_eq!(o.aa_len(seq.len()), 5);
    }

    #[test]
    fn nested_orfs_collapse() {
        // 同一框架里 ATG…ATG…TAA：只保留最外层
        let seq = b"ATGCCCATGCCCTAAGGG".to_vec();
        let orfs = find_orfs(&seq, false, 3);
        let fwd: Vec<&Orf> = orfs.iter().filter(|o| o.strand == 1).collect();
        assert_eq!(fwd.len(), 1, "{orfs:?}");
        assert_eq!((fwd[0].start, fwd[0].end), (0, 12));
    }

    #[test]
    fn default_length_counts_the_terminal_stop_codon() {
        let accepted = [
            b"ATG".as_slice(),
            b"GCT".repeat(73).as_slice(),
            b"TAA".as_slice(),
        ]
        .concat();
        let rejected = [
            b"ATG".as_slice(),
            b"GCT".repeat(72).as_slice(),
            b"TAA".as_slice(),
        ]
        .concat();
        assert!(find_orfs(&accepted, false, DEFAULT_MIN_BP)
            .iter()
            .any(|o| o.strand == 1 && o.start == 0));
        assert!(!find_orfs(&rejected, false, DEFAULT_MIN_BP)
            .iter()
            .any(|o| o.strand == 1 && o.start == 0));
        let o = find_orfs(&accepted, false, DEFAULT_MIN_BP)
            .into_iter()
            .find(|o| o.strand == 1 && o.start == 0)
            .unwrap();
        assert_eq!(o.end, 222);
        assert_eq!(o.aa_len(accepted.len()), 74);
    }

    #[test]
    fn terminal_stop_threshold_works_on_reverse_strand_and_across_origin() {
        let coding = [
            b"ATG".as_slice(),
            b"GCT".repeat(73).as_slice(),
            b"TGA".as_slice(),
        ]
        .concat();
        let rev: Vec<u8> = coding
            .iter()
            .rev()
            .map(|b| crate::edit_session::complement(*b))
            .collect();
        let reverse = find_orfs(&rev, false, DEFAULT_MIN_BP);
        assert!(
            reverse.iter().any(|o| o.strand == -1
                && (o.start, o.end) == (3, 225)
                && o.aa_len(rev.len()) == 74),
            "{reverse:?}"
        );

        // 避免整圈与空区间使用相同端点；将真实终止子和编码区分置于原点两侧。
        let mut ring = [coding.as_slice(), b"CCC"].concat();
        ring.rotate_left(100);
        let forward = find_orfs(&ring, true, DEFAULT_MIN_BP);
        assert!(
            forward.iter().any(|o| o.strand == 1
                && (o.start, o.end) == (128, 122)
                && o.aa_len(ring.len()) == 74),
            "{forward:?}"
        );
        let ring_rev: Vec<u8> = ring
            .iter()
            .rev()
            .map(|b| crate::edit_session::complement(*b))
            .collect();
        let reverse = find_orfs(&ring_rev, true, DEFAULT_MIN_BP);
        assert!(
            reverse.iter().any(|o| o.strand == -1
                && (o.start, o.end) == (106, 100)
                && o.aa_len(ring.len()) == 74),
            "{reverse:?}"
        );
    }

    #[test]
    fn threshold_does_not_credit_a_missing_stop_and_handles_startless_linear_end() {
        assert!(find_orfs(b"TAATAA", false, 3).is_empty());
        let short = [b"ATG".as_slice(), b"GCT".repeat(73).as_slice()].concat();
        assert!(!find_orfs(&short, false, DEFAULT_MIN_BP)
            .iter()
            .any(|o| o.strand == 1 && o.start == 0));
        let long = [short.as_slice(), b"GCT"].concat();
        assert!(find_orfs(&long, false, DEFAULT_MIN_BP)
            .iter()
            .any(|o| o.strand == 1 && (o.start, o.end) == (0, 225)));
        let startless = [b"GCT".repeat(74).as_slice(), b"TAG".as_slice()].concat();
        assert!(find_orfs(&startless, false, DEFAULT_MIN_BP)
            .iter()
            .any(|o| o.strand == 1 && !o.has_start && (o.start, o.end) == (0, 222)));
    }

    #[test]
    fn linear_end_coordinates_do_not_wrap_and_circular_spans_do_not_exceed_one_turn() {
        let fwd = find_orfs(b"ATGGCTGCT", false, 9);
        assert!(fwd
            .iter()
            .any(|o| o.strand == 1 && (o.start, o.end) == (0, 9)));
        let rev = find_orfs(b"AGCAGCCAT", false, 9);
        assert!(rev
            .iter()
            .any(|o| o.strand == -1 && (o.start, o.end) == (0, 9)));
        assert!(find_orfs(b"ATGCCCCCCC", true, 3).is_empty());
    }
}
