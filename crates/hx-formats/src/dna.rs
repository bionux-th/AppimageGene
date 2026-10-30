//! `.dna` / `.rna` / `.prot` 顶层模型：头部块（9）+ 块序列。
//!
//! 读取保持全部载荷原样；`to_bytes()` 重建的字节流与输入一致（未编辑时）。

use crate::block::{self, Block};
use crate::error::FormatError;
use hx_core::{MoleculeKind, SeqRecord};

pub const BLOCK_HEADER: u8 = 9;
pub const BLOCK_DNA: u8 = 0;
pub const BLOCK_PROTEIN: u8 = 21;
pub const BLOCK_RNA: u8 = 32;

const HEADER_MAGIC: &[u8; 8] = b"SnapGene";

/// 头部信息（block 9 载荷的 14 字节）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// 1 = DNA，2 = 蛋白，7 = RNA。
    pub type_of_sequence: u16,
    /// 写出格式版本（8.2.2 恒为 15）。
    pub export_version: u16,
    /// 导入格式版本（8.2.2 恒为 19）。
    pub import_version: u16,
    /// 原始 14 字节载荷（写回字节保真）。
    pub raw: Vec<u8>,
}

impl Header {
    /// 按 8.2.2 规范构造：`SnapGene` + tos + 15 + 19。
    pub fn new(type_of_sequence: u16) -> Self {
        let mut raw = Vec::with_capacity(14);
        raw.extend_from_slice(HEADER_MAGIC);
        raw.extend_from_slice(&type_of_sequence.to_be_bytes());
        raw.extend_from_slice(&15u16.to_be_bytes());
        raw.extend_from_slice(&19u16.to_be_bytes());
        Self {
            type_of_sequence,
            export_version: 15,
            import_version: 19,
            raw,
        }
    }

    pub fn parse(payload: &[u8]) -> Result<Self, FormatError> {
        if payload.len() != 14 {
            return Err(FormatError::BadHeader(format!(
                "载荷 {} 字节（应为 14）",
                payload.len()
            )));
        }
        if &payload[..8] != &HEADER_MAGIC[..] {
            return Err(FormatError::BadHeader("魔术字不是 SnapGene".into()));
        }
        let u16_at = |i: usize| u16::from_be_bytes([payload[i], payload[i + 1]]);
        Ok(Self {
            type_of_sequence: u16_at(8),
            export_version: u16_at(10),
            import_version: u16_at(12),
            raw: payload.to_vec(),
        })
    }

    pub fn molecule_kind(&self) -> Option<MoleculeKind> {
        MoleculeKind::from_type_of_sequence(self.type_of_sequence)
    }

    pub fn to_block(&self) -> Block {
        Block::new(BLOCK_HEADER, self.raw.clone())
    }
}

/// `block 8`（`<AdditionalSequenceProperties>`）里与克隆粘性末端相关的字段。
///
/// SnapGene 的约定（实测 8.2.3 的 TOPO/TA 载体与 PCR 产物）：
///
/// - `UpstreamStickiness` / `DownstreamStickiness`：**负数 = 3′ 突出端，
///   正数 = 5′ 突出端，0 = 平末端**；数值是突出碱基数。
/// - 存储的序列是「标记碱基 + 顶链」：上游 3′ 突出端属于底链，其互补碱基
///   作为标记放在序列开头；下游 5′ 突出端属于底链，其互补碱基作为标记放在
///   序列末尾。因此顶链是 `seq[max(0,-up) .. len-max(0,down)]`。
/// - `OriginalFirstBase`：环状原始序列里对应「存储序列第 `ofb` 个碱基」的
///   1-based 位置；产物按它旋转回原始编号。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EndProperties {
    pub upstream_stickiness: i32,
    pub downstream_stickiness: i32,
    pub upstream_modification: String,
    pub downstream_modification: String,
    pub original_first_base: Option<usize>,
}

impl Default for EndProperties {
    fn default() -> Self {
        Self {
            upstream_stickiness: 0,
            downstream_stickiness: 0,
            upstream_modification: "Unmodified".into(),
            downstream_modification: "Unmodified".into(),
            original_first_base: None,
        }
    }
}

impl EndProperties {
    /// 存储并集中的顶链区间；只适用于线性双链 DNA。
    pub fn top_range(&self, stored_len: usize) -> std::ops::Range<usize> {
        let start = (self.upstream_stickiness.min(0).unsigned_abs() as usize).min(stored_len);
        let end = stored_len.saturating_sub(self.downstream_stickiness.max(0) as usize).max(start);
        start..end
    }

    /// 同一存储并集中的底链区间（沿顶链方向，从左向右）。
    pub fn bottom_range(&self, stored_len: usize) -> std::ops::Range<usize> {
        let start = (self.upstream_stickiness.max(0) as usize).min(stored_len);
        let end = stored_len.saturating_sub(self.downstream_stickiness.min(0).unsigned_abs() as usize).max(start);
        start..end
    }
}

/// 一份 SnapGene 文件：头部 + 其余块（不含头部块）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaFile {
    pub header: Header,
    pub blocks: Vec<Block>,
}

impl DnaFile {
    /// 解析整份文件。首个块必须是 block 9。
    pub fn parse(bytes: &[u8]) -> Result<Self, FormatError> {
        let mut all = block::parse_blocks(bytes)?;
        if all.is_empty() || all[0].block_type != BLOCK_HEADER {
            return Err(FormatError::MissingHeader);
        }
        let header = Header::parse(&all[0].payload)?;
        all.remove(0);
        Ok(Self {
            header,
            blocks: all,
        })
    }

    /// 序列化为完整文件字节。
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        block::write_blocks(&mut out, std::slice::from_ref(&self.header.to_block()));
        block::write_blocks(&mut out, &self.blocks);
        out
    }

    /// 取首个指定类型的块。
    pub fn block(&self, t: u8) -> Option<&Block> {
        self.blocks.iter().find(|b| b.block_type == t)
    }

    /// 取所有指定类型的块（如 block 11 历史节点）。
    pub fn blocks_of(&self, t: u8) -> impl Iterator<Item = &Block> + '_ {
        self.blocks.iter().filter(move |b| b.block_type == t)
    }

    /// `block 8` 的末端属性（粘性末端 / 修饰 / 原始首碱基）。
    ///
    /// 没有 block 8 时返回默认值（平末端、未修饰、无原始首碱基）。
    pub fn end_properties(&self) -> EndProperties {
        let Some(b) = self.block(8) else {
            return EndProperties::default();
        };
        let xml = String::from_utf8_lossy(&b.payload);
        let text = |tag: &str| -> Option<String> {
            let open = format!("<{tag}>");
            let close = format!("</{tag}>");
            let s = xml.find(&open)? + open.len();
            let e = xml[s..].find(&close)? + s;
            Some(xml[s..e].to_string())
        };
        let int = |tag: &str| text(tag).and_then(|v| v.trim().parse::<i32>().ok());
        EndProperties {
            upstream_stickiness: int("UpstreamStickiness").unwrap_or(0),
            downstream_stickiness: int("DownstreamStickiness").unwrap_or(0),
            upstream_modification: text("UpstreamModification")
                .unwrap_or_else(|| "Unmodified".into()),
            downstream_modification: text("DownstreamModification")
                .unwrap_or_else(|| "Unmodified".into()),
            original_first_base: text("OriginalFirstBase")
                .and_then(|v| v.trim().parse::<usize>().ok()),
        }
    }

    /// `block 6`（`<Notes>`）里的图谱标签（`UseCustomMapLabel=1` 时的
    /// `CustomMapLabel`）。SnapGene 的商品化载体列表用这个名字而不是文件名。
    pub fn custom_map_label(&self) -> Option<String> {
        let xml = String::from_utf8_lossy(&self.block(6)?.payload);
        let enabled=xml.split_once("<UseCustomMapLabel>")?.1.split_once("</UseCustomMapLabel>")?.0.trim();
        if enabled!="1" {return None;}
        let open = "<CustomMapLabel>";
        let s = xml.find(open)? + open.len();
        let e = xml[s..].find("</CustomMapLabel>")? + s;
        let label = xml[s..e].trim();
        (!label.is_empty()).then(|| label.replace("&lt;","<").replace("&gt;",">")
            .replace("&quot;","\"").replace("&apos;","'").replace("&amp;","&"))
    }

    /// 主序列记录（按分子类型取 block 0/32/21）。
    pub fn sequence(&self) -> Option<SeqRecord> {
        let kind = self.header.molecule_kind()?;
        let b = self.block(kind.sequence_block_type())?;
        SeqRecord::from_block_payload(kind, &b.payload)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hx_core::{Strandedness, Topology};

    #[test]
    fn map_label_honors_native_enable_flag_and_xml_entities() {
        let mut new=crate::write::NewFile::of(hx_core::MoleculeKind::Rna,b"ACGT".to_vec());
        new.name=Some("RNA & <test>".into());
        let mut file=DnaFile::parse(&new.to_bytes()).unwrap();
        assert_eq!(file.custom_map_label().as_deref(),Some("RNA & <test>"));
        let notes=&mut file.blocks.iter_mut().find(|b|b.block_type==6).unwrap().payload;
        *notes=String::from_utf8(notes.clone()).unwrap().replace("<UseCustomMapLabel>1</UseCustomMapLabel>","<UseCustomMapLabel>0</UseCustomMapLabel>").into_bytes();
        assert_eq!(file.custom_map_label(),None);
    }

    #[test]
    fn terminal_ranges_distinguish_real_strands_from_marker_bases() {
        let ends = EndProperties { upstream_stickiness: -1, downstream_stickiness: 2, ..Default::default() };
        assert_eq!(ends.top_range(15), 1..13);
        assert_eq!(ends.bottom_range(15), 0..15);
        let ends = EndProperties { upstream_stickiness: 4, downstream_stickiness: -2, ..Default::default() };
        assert_eq!(ends.top_range(15), 0..15);
        assert_eq!(ends.bottom_range(15), 4..13);
        let malformed = EndProperties { upstream_stickiness: i32::MIN, downstream_stickiness: i32::MAX, ..Default::default() };
        assert_eq!(malformed.top_range(15), 15..15);
        assert_eq!(malformed.bottom_range(15), 0..15);
        assert_eq!(malformed.top_range(0), 0..0);
    }

    fn minimal_file() -> Vec<u8> {
        let mut out = Vec::new();
        block::write_blocks(&mut out, &[Header::new(1).to_block()]);
        let mut seq = vec![0x1f];
        seq.extend_from_slice(b"ACGTACGT");
        block::write_blocks(&mut out, &[Block::new(BLOCK_DNA, seq)]);
        out
    }

    #[test]
    fn parse_and_roundtrip_minimal() {
        let bytes = minimal_file();
        let f = DnaFile::parse(&bytes).unwrap();
        assert_eq!(f.header.type_of_sequence, 1);
        assert_eq!(f.header.export_version, 15);
        assert_eq!(f.header.molecule_kind(), Some(MoleculeKind::Dna));
        assert_eq!(f.blocks.len(), 1);
        let rec = f.sequence().unwrap();
        assert_eq!(rec.kind, MoleculeKind::Dna);
        assert_eq!(rec.flags.topology(), Topology::Circular);
        assert_eq!(rec.flags.strandedness(), Strandedness::Double);
        assert_eq!(rec.sequence, b"ACGTACGT");
        assert_eq!(f.to_bytes(), bytes);
    }

    #[test]
    fn missing_header_is_error() {
        let mut out = Vec::new();
        block::write_blocks(&mut out, &[Block::new(BLOCK_DNA, vec![0x00])]);
        assert_eq!(DnaFile::parse(&out), Err(FormatError::MissingHeader));
    }

    #[test]
    fn empty_file_is_error() {
        assert_eq!(DnaFile::parse(&[]), Err(FormatError::MissingHeader));
    }
}
