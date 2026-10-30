//! 通用 TLV 块层。
//!
//! ```text
//! [block_type: u8][length: u32 BE][payload: length bytes]
//! ```
//!
//! 203/203 样例从偏移 0 起连续解析到 EOF（无填充/对齐）。载荷不经解释地
//! 原样保存，是读→写字节保真的基础。

use crate::error::FormatError;

/// 单块载荷防御性上限（512 MiB）。
pub const MAX_BLOCK_LEN: u32 = 512 * 1024 * 1024;

/// 一个 TLV 块。`payload` 为原始字节。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub block_type: u8,
    /// 读入时的块起始偏移（仅用于诊断；写出时忽略并重新计算）。
    pub offset: u32,
    pub payload: Vec<u8>,
}

impl Block {
    pub fn new(block_type: u8, payload: Vec<u8>) -> Self {
        Self {
            block_type,
            offset: 0,
            payload,
        }
    }
}

/// 解析整份文件为块序列。
pub fn parse_blocks(bytes: &[u8]) -> Result<Vec<Block>, FormatError> {
    let mut blocks = Vec::new();
    let mut o = 0usize;
    while o < bytes.len() {
        if o + 5 > bytes.len() {
            return Err(FormatError::Truncated {
                offset: o,
                needed: 5,
                avail: bytes.len() - o,
            });
        }
        let block_type = bytes[o];
        let len = u32::from_be_bytes(bytes[o + 1..o + 5].try_into().unwrap());
        if len > MAX_BLOCK_LEN {
            return Err(FormatError::BlockTooLarge { offset: o, len });
        }
        let start = o + 5;
        let end = start + len as usize;
        if end > bytes.len() {
            return Err(FormatError::Truncated {
                offset: start,
                needed: len as usize,
                avail: bytes.len() - start,
            });
        }
        blocks.push(Block {
            block_type,
            offset: o as u32,
            payload: bytes[start..end].to_vec(),
        });
        o = end;
    }
    Ok(blocks)
}

/// 序列化块到输出缓冲。
pub fn write_blocks(out: &mut Vec<u8>, blocks: &[Block]) {
    out.reserve(blocks.iter().map(|b| 5 + b.payload.len()).sum::<usize>());
    for b in blocks {
        out.push(b.block_type);
        out.extend_from_slice(&(b.payload.len() as u32).to_be_bytes());
        out.extend_from_slice(&b.payload);
    }
}

/// 块类型名称（8.2.2 实测规格；未知类型返回 `"unknown"`）。
pub fn block_type_name(t: u8) -> &'static str {
    match t {
        0 => "DNA sequence",
        1 => "compressed DNA sequence",
        2 => "opaque metadata (rare)",
        3 => "enzyme site cache",
        5 => "primers XML",
        6 => "notes XML",
        7 => "history tree XML",
        8 => "additional sequence properties XML",
        9 => "file header",
        10 => "features XML",
        11 => "history node snapshot",
        13 => "enzyme display preferences",
        14 => "custom enzyme sets XML",
        16 => "trace container",
        17 => "alignable sequences XML",
        18 => "ZTR chromatogram data",
        20 => "strand colors XML",
        21 => "protein sequence",
        23 => "attachments",
        27 => "trace alignment (BGZF BAM)",
        28 => "enzyme visibilities XML",
        29 => "history manipulation XML (LZMA)",
        30 => "history node content (nested)",
        32 => "RNA sequence",
        34 => "RNA secondary structure (LZMA JSON)",
        35 => "compressed empty JSON array (LZMA)",
        _ => "unknown",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_blocks() {
        let blocks = vec![
            Block::new(9, b"SnapGene 14B".to_vec()),
            Block::new(0, vec![0x1f, b'A', b'C', b'G', b'T']),
            Block::new(3, vec![1, 0, 0, 0x12, 0x9b]),
        ];
        let mut buf = Vec::new();
        write_blocks(&mut buf, &blocks);
        let back = parse_blocks(&buf).unwrap();
        assert_eq!(back.len(), 3);
        assert_eq!(back[0].block_type, 9);
        assert_eq!(back[1].payload, vec![0x1f, b'A', b'C', b'G', b'T']);
        assert_eq!(back[2].payload, vec![1, 0, 0, 0x12, 0x9b]);
        assert_eq!(back[2].offset as usize, 5 + 12 + 5 + 5); // "SnapGene 14B" = 12 字节
    }

    #[test]
    fn truncated_is_error() {
        let mut buf = Vec::new();
        write_blocks(&mut buf, &[Block::new(0, vec![1, 2, 3])]);
        buf.truncate(buf.len() - 1);
        assert!(matches!(
            parse_blocks(&buf),
            Err(FormatError::Truncated { .. })
        ));
    }
}
