//! `block 3`（酶位点分析缓存）编解码。
//!
//! 布局（由 `docs/SnapGene_dna_format_spec.md` §5 全量解码，203/203 样例验证）：
//!
//! ```text
//! 01                       // 1B 版本标记
//! <u32 BE 4763>            // 站点表字节长
//! <4763B>                  // 473 条逗号分隔 IUPAC 位点（[`sites::BLOCK3_SITES`]）
//! <473 × u32 BE>           // 各位点命中计数
//! <473 × u16 BE>           // 同上（16 位副本，与 u32 表逐值相同）
//! <7 × u32 BE>             // 位置流分桶边界的**文件绝对偏移**
//! <N × u16 LE>             // 命中起点位置流（按站点序分组，组内升序）
//! ```
//!
//! 计数语义 = **环状感知的正向命中数**（与 [`crate::find_site_positions`] 一致）：
//! 线性只算序列内完整匹配；环状算起点 ∈ `[0, L)` 的全部匹配（含跨原点命中，
//! 存原始起点、不取模）。
//!
//! 7 个偏移的语义（自载荷首字节起算 7634 处为位置流起点）：
//! `[流起点, @184, @221, @250, @392, @450, 流结束]`，对应
//! `[`[`sites::BLOCK3_BUCKETS`]`]` 的 6 个分桶边界。

#![forbid(unsafe_code)]

use crate::find_site_positions;
use crate::sites::{BLOCK3_BUCKETS, BLOCK3_SITES, BLOCK3_SITE_COUNT, BLOCK3_SITES_LEN};

/// 载荷头部固定长度（版本 1B + 长度 4B + 站点表 + 两份计数表 + 7 偏移）。
pub const PAYLOAD_HEADER_LEN: usize = 1 + 4 + BLOCK3_SITES_LEN + BLOCK3_SITE_COUNT * 4 + BLOCK3_SITE_COUNT * 2 + 28;

/// 位置流起点相对载荷首字节的偏移（7634）。
pub const POSITION_STREAM_OFFSET: usize = PAYLOAD_HEADER_LEN;

/// 解码结果。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block3Data {
    /// 站点表原文（逗号分隔）。
    pub site_text: String,
    /// 各站点计数（u32 表；解码时会校验 u16 副本与之一致）。
    pub counts: Vec<u32>,
    /// 7 个文件绝对偏移。
    pub offsets: [u32; 7],
    /// 各站点命中起点（与 `counts` 等长；组内升序）。
    pub positions: Vec<Vec<u32>>,
}

/// 解码 block 3 载荷。遇到结构性错误返回 `Err(说明)`。
pub fn decode(payload: &[u8]) -> Result<Block3Data, String> {
    if payload.len() < 5 {
        return Err(format!("载荷过短：{} 字节", payload.len()));
    }
    if payload[0] != 1 {
        return Err(format!("版本标记不是 01：{:#x}", payload[0]));
    }
    let site_len = u32::from_be_bytes(payload[1..5].try_into().unwrap()) as usize;
    let sites_end = 5 + site_len;
    if sites_end > payload.len() {
        return Err(format!("站点表越界：需 {} 字节", sites_end));
    }
    let site_text = std::str::from_utf8(&payload[5..sites_end])
        .map_err(|e| format!("站点表非 UTF-8：{e}"))?
        .to_string();

    let n = if site_text.is_empty() {
        0
    } else {
        site_text.split(',').count()
    };
    let c32_end = sites_end + n * 4;
    let c16_end = c32_end + n * 2;
    let off_end = c16_end + 28;
    if off_end > payload.len() {
        return Err(format!("计数/偏移表越界：需 {} 字节", off_end));
    }
    let counts: Vec<u32> = (0..n)
        .map(|i| u32::from_be_bytes(payload[sites_end + i * 4..sites_end + i * 4 + 4].try_into().unwrap()))
        .collect();
    let mut counts16 = Vec::with_capacity(n);
    for i in 0..n {
        let v = u16::from_be_bytes(payload[c32_end + i * 2..c32_end + i * 2 + 2].try_into().unwrap());
        counts16.push(v as u32);
    }
    if counts != counts16 {
        return Err("u32 计数表与 u16 副本不一致".into());
    }
    let mut offsets = [0u32; 7];
    for (i, o) in offsets.iter_mut().enumerate() {
        *o = u32::from_be_bytes(payload[c16_end + i * 4..c16_end + i * 4 + 4].try_into().unwrap());
    }

    let stream = &payload[off_end..];
    if stream.len() % 2 != 0 {
        return Err(format!("位置流长度 {} 不是偶数", stream.len()));
    }
    let total: usize = counts.iter().map(|c| *c as usize).sum();
    if stream.len() / 2 != total {
        return Err(format!(
            "位置流条目 {} 与计数和 {} 不符",
            stream.len() / 2,
            total
        ));
    }
    let mut positions = Vec::with_capacity(n);
    let mut cur = 0usize;
    for c in &counts {
        let n = *c as usize;
        let mut v = Vec::with_capacity(n);
        for k in 0..n {
            let at = (cur + k) * 2;
            v.push(u16::from_le_bytes(stream[at..at + 2].try_into().unwrap()) as u32);
        }
        cur += n;
        positions.push(v);
    }
    Ok(Block3Data {
        site_text,
        counts,
        offsets,
        positions,
    })
}

/// 用规范站点表（[`sites::BLOCK3_SITES`]）为序列重建 block 3 载荷。
///
/// `payload_offset` = 本载荷首字节（版本标记）在**文件中的绝对偏移**——
/// 7 个偏移字段写的是文件绝对偏移，因此调用方必须在最终块序确定后传入。
pub fn encode(seq: &[u8], circular: bool, payload_offset: u64) -> Vec<u8> {
    encode_with_sites(BLOCK3_SITES, seq, circular, payload_offset)
}

/// 与 [`encode`] 相同，但可指定站点表（用于回放任意历史文件的表）。
pub fn encode_with_sites(site_text: &str, seq: &[u8], circular: bool, payload_offset: u64) -> Vec<u8> {
    let sites: Vec<&str> = site_text.split(',').collect();
    let n = sites.len();
    let mut counts: Vec<u32> = Vec::with_capacity(n);
    let mut positions: Vec<Vec<u32>> = Vec::with_capacity(n);
    for s in &sites {
        let pos = find_site_positions(seq, s.as_bytes(), circular);
        counts.push(pos.len() as u32);
        positions.push(pos);
    }
    encode_tables(site_text, &counts, &positions, payload_offset)
}

/// 由计数与位置表直接编码（供需要自行计算/缓存命中的调用方使用）。
///
/// `positions` 必须与 `counts` 等长，且 `positions[i].len() == counts[i]`。
pub fn encode_tables(
    site_text: &str,
    counts: &[u32],
    positions: &[Vec<u32>],
    payload_offset: u64,
) -> Vec<u8> {
    assert_eq!(counts.len(), positions.len(), "计数表与位置表长度不一致");
    for (i, (c, p)) in counts.iter().zip(positions).enumerate() {
        assert_eq!(*c as usize, p.len(), "第 {i} 个站点的计数与位置数不符");
    }
    let n = counts.len();
    let header_len = 1 + 4 + site_text.len() + n * 4 + n * 2 + 28;
    let total_hits: u64 = counts.iter().map(|c| *c as u64).sum();

    let mut out = Vec::with_capacity(header_len + 2 * total_hits as usize);
    out.push(1);
    out.extend_from_slice(&(site_text.len() as u32).to_be_bytes());
    out.extend_from_slice(site_text.as_bytes());

    for c in counts {
        out.extend_from_slice(&c.to_be_bytes()); // u32 BE 表
    }
    for c in counts {
        let v = u16::try_from(*c).unwrap_or_else(|_| panic!("计数 {c} 超过 u16 上限（block 3 格式限制）"));
        out.extend_from_slice(&v.to_be_bytes()); // u16 BE 副本
    }

    // 7 个文件绝对偏移：流起点 + 5 个分桶边界（2 字节/命中）+ 流结束
    let stream_start = payload_offset + header_len as u64;
    let mut offsets = [stream_start as u32; 7];
    if n == BLOCK3_SITE_COUNT {
        let mut acc: u64 = 0;
        for (i, (lo, hi)) in BLOCK3_BUCKETS.iter().enumerate() {
            acc += counts[*lo..*hi].iter().map(|c| *c as u64).sum::<u64>();
            offsets[i + 1] = (stream_start + 2 * acc) as u32;
        }
    } else {
        // 非规范站点表长度：各边界退化为流结束（仅用于回放/测试）
        let end = (stream_start + 2 * total_hits) as u32;
        for o in offsets.iter_mut().skip(1) {
            *o = end;
        }
    }
    for o in offsets {
        out.extend_from_slice(&o.to_be_bytes());
    }

    for p in positions {
        for start in p {
            let v = u16::try_from(*start).unwrap_or_else(|_| {
                panic!("命中起点 {start} 超过 u16 上限（block 3 格式限制）")
            });
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    debug_assert_eq!(out.len() as u64, header_len as u64 + 2 * total_hits);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_small() {
        // 小序列上做一次编码→解码回环（把载荷放在文件偏移 1000 处）
        let seq = b"GAATTCGGATCCGAATTC";
        let payload = encode(seq, false, 1000);
        let d = decode(&payload).unwrap();
        assert_eq!(d.site_text, BLOCK3_SITES);
        assert_eq!(d.counts.len(), BLOCK3_SITE_COUNT);
        assert_eq!(d.positions.len(), BLOCK3_SITE_COUNT);
        assert_eq!(d.offsets[0], 1000 + POSITION_STREAM_OFFSET as u32);
        // EcoRI 位点 GAATTC 有两处（0、12）
        let i = crate::sites::sites().position(|s| s == "GAATTC").unwrap();
        assert_eq!(d.counts[i], 2);
        assert_eq!(d.positions[i], vec![0, 12]);
        // 位置流结束 = 流起点 + 2·Σ计数
        let total: u32 = d.counts.iter().sum();
        assert_eq!(d.offsets[6], d.offsets[0] + 2 * total);
    }

    #[test]
    fn circular_counts_wrap_origin() {
        // ATAT 环状：TA 命中 1、3（3 跨原点）；线性只有 1
        let lin = encode_with_sites("TA", b"ATAT", false, 0);
        let cir = encode_with_sites("TA", b"ATAT", true, 0);
        assert_eq!(decode(&lin).unwrap().counts, vec![1]);
        assert_eq!(decode(&cir).unwrap().counts, vec![1 + 1]);
        assert_eq!(decode(&cir).unwrap().positions[0], vec![1, 3]);
    }
}
