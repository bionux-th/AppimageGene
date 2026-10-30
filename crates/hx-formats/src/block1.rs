//! Block 1 — 压缩 DNA 序列（历史节点快照的序列存储；顶层罕见）。
//!
//! 结构（详见 `docs/SnapGene_dna_format_spec.md` §4）：
//!
//! ```text
//! off 0  4B  cl          // “cl 之后所有字节”的总长（BE u32）
//! off 4  4B  ul          // 最终解码字符数
//! off 8  1B  stamp       // 写入方标记（SnapGene 解析时忽略；保真沿用原值）
//! off 9  4B  chunks      // section 数量
//! off 13 4B  lowercases  // 小写区间对数
//! off 17 1B  marker      // 第一个 section 标记
//! off 18 4B  count       // 第一个 section 字符数
//! off 22 …   第一个 section 数据；其后每个 section 自带 5B 帧（marker+count）
//! 末尾   lowercases × 8B = (start:u32 BE, end:u32 BE) 闭区间
//! ```
//!
//! Section 标记：`0x01` 纯 A/C/G/T（2-bit 打包，`00=G, 01=A, 10=T, 11=C`，
//! MSB 在前）、`0x02` IUPAC（4-bit 半字节，高半字节在前）、`0x03` 连续 N。
//!
//! 解析 → `to_payload()` 必须与原字节一致（保真）；`from_sequence()` 用于
//! 新写出（如历史节点）。

use crate::error::FormatError;

pub const SECTION_PLAIN: u8 = 0x01;
pub const SECTION_IUPAC: u8 = 0x02;
pub const SECTION_N_RUN: u8 = 0x03;

/// 2-bit 值 → 碱基：`00=G, 01=A, 10=T, 11=C`。
const PLAIN_BASES: [u8; 4] = [b'G', b'A', b'T', b'C'];

/// 4-bit 半字节 → 碱基（0x0-0x3 与 2-bit 同序；0xF 按 N 处理）。
const IUPAC_BASES: [u8; 16] = [
    b'G', b'A', b'T', b'C', b'N', b'B', b'D', b'H', b'K', b'M', b'R', b'S', b'V', b'W', b'Y', b'N',
];

/// 防御性上限：section 数与小写区间数。
const MAX_SECTIONS: u32 = 1_000_000;
const MAX_LOWERCASE_RANGES: u32 = 1_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Section {
    pub marker: u8,
    pub count: u32,
    /// 原始打包数据（解析时原样保留；`SECTION_N_RUN` 恒为空）。
    pub data: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompressedSeq {
    pub stamp: u8,
    pub uncompressed_len: u32,
    pub sections: Vec<Section>,
    /// 小写区间（闭区间，作用于拼装后的序列）。
    pub lowercases: Vec<(u32, u32)>,
}

impl CompressedSeq {
    /// 解析 block 1 载荷。
    pub fn parse(payload: &[u8]) -> Result<Self, FormatError> {
        let need = |o: usize, n: usize| -> Result<(), FormatError> {
            match o.checked_add(n) {
                Some(end) if end <= payload.len() => Ok(()),
                _ => Err(FormatError::Truncated {
                    offset: o,
                    needed: n,
                    avail: payload.len().saturating_sub(o),
                }),
            }
        };

        need(0, 22)?;
        let cl = u32::from_be_bytes(payload[0..4].try_into().unwrap());
        if cl as usize != payload.len() - 4 {
            return Err(FormatError::BadCompressedSeq(format!(
                "cl={cl} 与载荷长 {} 不符",
                payload.len()
            )));
        }
        let uncompressed_len = u32::from_be_bytes(payload[4..8].try_into().unwrap());
        let stamp = payload[8];
        let chunks = u32::from_be_bytes(payload[9..13].try_into().unwrap());
        let lowercase_count = u32::from_be_bytes(payload[13..17].try_into().unwrap());
        if chunks > MAX_SECTIONS || lowercase_count > MAX_LOWERCASE_RANGES {
            return Err(FormatError::BadCompressedSeq(format!(
                "chunks={chunks} lowercases={lowercase_count} 超上限"
            )));
        }

        let mut o = 17usize;
        let mut sections = Vec::with_capacity(chunks as usize);
        for _ in 0..chunks {
            need(o, 5)?;
            let marker = payload[o];
            let count = u32::from_be_bytes(payload[o + 1..o + 5].try_into().unwrap());
            o += 5;
            let n = section_data_len(marker, count)?;
            need(o, n)?;
            sections.push(Section {
                marker,
                count,
                data: payload[o..o + n].to_vec(),
            });
            o += n;
        }

        let mut lowercases = Vec::with_capacity(lowercase_count as usize);
        for _ in 0..lowercase_count {
            need(o, 8)?;
            let start = u32::from_be_bytes(payload[o..o + 4].try_into().unwrap());
            let end = u32::from_be_bytes(payload[o + 4..o + 8].try_into().unwrap());
            lowercases.push((start, end));
            o += 8;
        }
        if o != payload.len() {
            return Err(FormatError::BadCompressedSeq(format!(
                "解析到 {o}，载荷共 {} 字节",
                payload.len()
            )));
        }

        let total: u64 = sections.iter().map(|s| s.count as u64).sum();
        if total != uncompressed_len as u64 {
            return Err(FormatError::BadCompressedSeq(format!(
                "section 字符和 {total} != ul {uncompressed_len}"
            )));
        }
        Ok(Self {
            stamp,
            uncompressed_len,
            sections,
            lowercases,
        })
    }

    /// 按解析结果原样序列化（与输入字节一致）。
    pub fn to_payload(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(22 + self.sections.iter().map(|s| 5 + s.data.len()).sum::<usize>());
        out.extend_from_slice(&[0, 0, 0, 0]); // cl 占位
        out.extend_from_slice(&self.uncompressed_len.to_be_bytes());
        out.push(self.stamp);
        out.extend_from_slice(&(self.sections.len() as u32).to_be_bytes());
        out.extend_from_slice(&(self.lowercases.len() as u32).to_be_bytes());
        for s in &self.sections {
            out.push(s.marker);
            out.extend_from_slice(&s.count.to_be_bytes());
            out.extend_from_slice(&s.data);
        }
        for &(a, b) in &self.lowercases {
            out.extend_from_slice(&a.to_be_bytes());
            out.extend_from_slice(&b.to_be_bytes());
        }
        let cl = (out.len() - 4) as u32;
        out[0..4].copy_from_slice(&cl.to_be_bytes());
        out
    }

    /// 拼装为 ASCII 序列（含小写还原）。
    pub fn to_sequence(&self) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::with_capacity(self.uncompressed_len as usize);
        for s in &self.sections {
            let start = out.len();
            let want = s.count as usize;
            match s.marker {
                SECTION_PLAIN => {
                    'outer: for &byte in &s.data {
                        for sh in [6u32, 4, 2, 0] {
                            if out.len() - start >= want {
                                break 'outer;
                            }
                            out.push(PLAIN_BASES[((byte >> sh) & 3) as usize]);
                        }
                    }
                }
                SECTION_IUPAC => {
                    'outer: for &byte in &s.data {
                        for nib in [byte >> 4, byte & 0x0f] {
                            if out.len() - start >= want {
                                break 'outer;
                            }
                            out.push(IUPAC_BASES[(nib & 0x0f) as usize]);
                        }
                    }
                }
                SECTION_N_RUN => {
                    out.resize(start + want, b'N');
                }
                _ => {}
            }
        }
        for &(s, e) in &self.lowercases {
            let (s, e) = (s as usize, e as usize);
            if s <= e && e < out.len() {
                for b in &mut out[s..=e] {
                    *b = b.to_ascii_lowercase();
                }
            }
        }
        out
    }

    /// 从 ASCII 序列构造（自动分组 0x01/0x02/0x03 section 并归并小写区间）。
    ///
    /// `stamp` 为写入方标记字节；SnapGene 解析时忽略，保真写回时沿用原值，
    /// 新写出暂用 `0x1f`（待 P1 oracle 实验确认）。
    pub fn from_sequence(seq: &[u8], stamp: u8) -> Self {
        const CLASS_PLAIN: u8 = 0;
        const CLASS_N: u8 = 1;
        const CLASS_IUPAC: u8 = 2;
        fn classify(up: u8) -> u8 {
            match up {
                b'A' | b'C' | b'G' | b'T' => CLASS_PLAIN,
                b'N' => CLASS_N,
                _ => CLASS_IUPAC,
            }
        }

        // 小写区间（整段扫描，跨 section 连续小写会归并）
        let mut lowercases: Vec<(u32, u32)> = Vec::new();
        let mut k = 0usize;
        while k < seq.len() {
            if seq[k].is_ascii_lowercase() {
                let s = k;
                while k < seq.len() && seq[k].is_ascii_lowercase() {
                    k += 1;
                }
                lowercases.push((s as u32, (k - 1) as u32));
            } else {
                k += 1;
            }
        }

        let mut sections: Vec<Section> = Vec::new();
        let mut i = 0usize;
        while i < seq.len() {
            let class = classify(seq[i].to_ascii_uppercase());
            let mut j = i;
            while j < seq.len() && classify(seq[j].to_ascii_uppercase()) == class {
                j += 1;
            }
            let run = &seq[i..j];
            let count = run.len() as u32;
            match class {
                CLASS_PLAIN => {
                    let mut data = Vec::with_capacity(run.len().div_ceil(4));
                    let mut acc: u32 = 0;
                    let mut n: u32 = 0;
                    for &b in run {
                        let bits = match b.to_ascii_uppercase() {
                            b'G' => 0,
                            b'A' => 1,
                            b'T' => 2,
                            _ => 3, // C
                        };
                        acc = (acc << 2) | bits;
                        n += 1;
                        if n == 4 {
                            data.push(acc as u8);
                            acc = 0;
                            n = 0;
                        }
                    }
                    if n > 0 {
                        data.push(((acc << (2 * (4 - n))) & 0xff) as u8);
                    }
                    sections.push(Section {
                        marker: SECTION_PLAIN,
                        count,
                        data,
                    });
                }
                CLASS_IUPAC => {
                    let mut data = Vec::with_capacity(run.len().div_ceil(2));
                    let mut high: Option<u8> = None;
                    for &b in run {
                        let nib = iupac_nibble(b.to_ascii_uppercase());
                        match high.take() {
                            None => high = Some(nib),
                            Some(h) => data.push((h << 4) | nib),
                        }
                    }
                    if let Some(h) = high {
                        data.push(h << 4); // 低半字节补 0（待写回验证）
                    }
                    sections.push(Section {
                        marker: SECTION_IUPAC,
                        count,
                        data,
                    });
                }
                _ => {
                    sections.push(Section {
                        marker: SECTION_N_RUN,
                        count,
                        data: Vec::new(),
                    });
                }
            }
            i = j;
        }

        Self {
            stamp,
            uncompressed_len: seq.len() as u32,
            sections,
            lowercases,
        }
    }
}

fn section_data_len(marker: u8, count: u32) -> Result<usize, FormatError> {
    Ok(match marker {
        SECTION_PLAIN => {
            if count > i32::MAX as u32 {
                return Err(FormatError::BadCompressedSeq(format!("count={count} 过大")));
            }
            count.div_ceil(4) as usize
        }
        SECTION_IUPAC => count.div_ceil(2) as usize,
        SECTION_N_RUN => 0,
        other => {
            return Err(FormatError::BadCompressedSeq(format!(
                "未知 section 标记 0x{other:02x}"
            )))
        }
    })
}

fn iupac_nibble(up: u8) -> u8 {
    match up {
        b'G' => 0x0,
        b'A' => 0x1,
        b'T' => 0x2,
        b'C' => 0x3,
        b'N' => 0x4,
        b'B' => 0x5,
        b'D' => 0x6,
        b'H' => 0x7,
        b'K' => 0x8,
        b'M' => 0x9,
        b'R' => 0xa,
        b'S' => 0xb,
        b'V' => 0xc,
        b'W' => 0xd,
        b'Y' => 0xe,
        _ => 0x4, // 未知字符按 N
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_gatc_msb_first() {
        let c = CompressedSeq::from_sequence(b"GATC", 0);
        assert_eq!(c.sections.len(), 1);
        assert_eq!(c.sections[0].data, vec![0b00_01_10_11]);
    }

    #[test]
    fn synthetic_roundtrip_mixed() {
        let seq: &[u8] = b"ACGTNNRYatcg";
        let c = CompressedSeq::from_sequence(seq, 0x1f);
        let payload = c.to_payload();
        let back = CompressedSeq::parse(&payload).unwrap();
        assert_eq!(back.to_sequence(), seq);
        assert_eq!(back.uncompressed_len, seq.len() as u32);
        assert_eq!(back.lowercases, vec![(8, 11)]);
        let markers: Vec<u8> = back.sections.iter().map(|s| s.marker).collect();
        assert_eq!(
            markers,
            vec![SECTION_PLAIN, SECTION_N_RUN, SECTION_IUPAC, SECTION_PLAIN]
        );
        // 再序列化必须字节一致
        assert_eq!(back.to_payload(), payload);
    }

    #[test]
    fn partial_last_byte_padding() {
        // 5 个碱基：最后一个字节仅高 2 bit 有效
        let c = CompressedSeq::from_sequence(b"GGGGG", 0);
        assert_eq!(c.sections[0].data, vec![0x00, 0x00]);
        let back = CompressedSeq::parse(&c.to_payload()).unwrap();
        assert_eq!(back.to_sequence(), b"GGGGG");
    }

    #[test]
    fn rejects_bad_cl() {
        let mut payload = CompressedSeq::from_sequence(b"ACGT", 0).to_payload();
        payload[0] = 0xff;
        assert!(matches!(
            CompressedSeq::parse(&payload),
            Err(FormatError::BadCompressedSeq(_))
        ));
    }
}
