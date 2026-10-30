//! 分子类型与序列记录。

use crate::flags::MoleculeFlags;

/// 分子类型。
///
/// 与 `.dna` 头部 `type_of_sequence` 的映射（实测）：`1 = DNA`、`2 = 蛋白`、`7 = RNA`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MoleculeKind {
    Dna,
    Rna,
    Protein,
}

impl MoleculeKind {
    pub const fn from_type_of_sequence(v: u16) -> Option<Self> {
        match v {
            1 => Some(Self::Dna),
            2 => Some(Self::Protein),
            7 => Some(Self::Rna),
            _ => None,
        }
    }

    pub const fn to_type_of_sequence(self) -> u16 {
        match self {
            Self::Dna => 1,
            Self::Protein => 2,
            Self::Rna => 7,
        }
    }

    pub const fn display_name(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Protein => "Protein",
        }
    }

    /// 序列块类型：DNA = 0，蛋白 = 21，RNA = 32。
    pub const fn sequence_block_type(self) -> u8 {
        match self {
            Self::Dna => 0,
            Self::Rna => 32,
            Self::Protein => 21,
        }
    }

    /// 长度单位（SnapGene 界面用词，oracle 实测）：
    /// DNA 写 `bp`、RNA 写 `bases`、蛋白写 `aa`。
    pub const fn unit(self) -> &'static str {
        match self {
            Self::Dna => "bp",
            Self::Rna => "bases",
            Self::Protein => "aa",
        }
    }

    /// 是否为核苷酸序列（DNA/RNA）。酶切、互补链、翻译只对核苷酸有意义。
    pub const fn is_nucleotide(self) -> bool {
        !matches!(self, Self::Protein)
    }

    /// 特征是否**有方向**（画箭头）：DNA/RNA 是，蛋白不是（SnapGene 画无方向方块）。
    pub const fn features_have_direction(self) -> bool {
        self.is_nucleotide()
    }

    /// 序列视图每行多少个单位（oracle 实测：蛋白与 RNA 都是 60 —— 蛋白刻度尺
    /// 第一行写 `1 10 … 60`、第二行从 `70` 起，RNA 同样每行 60 个碱基）。
    pub const fn default_row_units(self) -> usize {
        60
    }

    /// 序列行标尺的坐标前缀（蛋白写 `p10`/`p20`…，核苷酸直接写数字）。
    pub const fn position_prefix(self) -> &'static str {
        match self {
            Self::Protein => "p",
            _ => "",
        }
    }

    /// 文档/页签/侧栏图标的分子类型颜色（自 oracle 截图取样：
    /// DNA `#DB5550`、RNA `#6CC780`、蛋白 `#6CC0DC`）。
    pub const fn icon_rgb(self) -> (u8, u8, u8) {
        match self {
            Self::Dna => (0xdb, 0x55, 0x50),
            Self::Rna => (0x6c, 0xc7, 0x80),
            Self::Protein => (0x6c, 0xc0, 0xdc),
        }
    }

    /// 显示用碱基/残基字母。`.rna` 文件里存的是 DNA 字母（`t`），
    /// SnapGene 的 RNA 视图显示为 `u`（oracle：`aaa u u g a a g a g u u u g`），
    /// 所以 RNA 要把 T/t 映射成 U/u；其它分子类型原样返回。
    pub const fn display_base(self, b: u8) -> u8 {
        match self {
            Self::Rna => match b {
                b't' => b'u',
                b'T' => b'U',
                other => other,
            },
            _ => b,
        }
    }
}

/// 一条序列记录：分子类型 + 属性字节 + ASCII 序列（保留原始大小写）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    pub kind: MoleculeKind,
    pub flags: MoleculeFlags,
    /// ASCII 序列字节，原样保留（`N`/IUPAC/小写不在此层规范）。
    pub sequence: Vec<u8>,
}

impl SeqRecord {
    pub fn new(kind: MoleculeKind, flags: MoleculeFlags, sequence: Vec<u8>) -> Self {
        Self {
            kind,
            flags,
            sequence,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// 按序列块载荷格式解码（[`Self::to_block_payload`] 的逆操作）。
    pub fn from_block_payload(kind: MoleculeKind, payload: &[u8]) -> Option<Self> {
        let (&flags, sequence) = payload.split_first()?;
        Some(Self {
            kind,
            flags: MoleculeFlags::from_byte(flags),
            sequence: sequence.to_vec(),
        })
    }

    /// 按序列块载荷格式编码：`[flags: u8][ASCII 序列]`。
    pub fn to_block_payload(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(1 + self.sequence.len());
        out.push(self.flags.to_byte());
        out.extend_from_slice(&self.sequence);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn type_of_sequence_mapping() {
        assert_eq!(MoleculeKind::from_type_of_sequence(1), Some(MoleculeKind::Dna));
        assert_eq!(MoleculeKind::from_type_of_sequence(2), Some(MoleculeKind::Protein));
        assert_eq!(MoleculeKind::from_type_of_sequence(7), Some(MoleculeKind::Rna));
        assert_eq!(MoleculeKind::from_type_of_sequence(3), None);
        assert_eq!(MoleculeKind::Rna.to_type_of_sequence(), 7);
        assert_eq!(MoleculeKind::Rna.sequence_block_type(), 32);
    }

    /// `.rna` 文件存的是 DNA 字母，RNA 视图要显示成 `u`（oracle 实测）。
    #[test]
    fn rna_display_base_maps_t_to_u() {
        assert_eq!(MoleculeKind::Rna.display_base(b't'), b'u');
        assert_eq!(MoleculeKind::Rna.display_base(b'T'), b'U');
        assert_eq!(MoleculeKind::Rna.display_base(b'a'), b'a');
        assert_eq!(MoleculeKind::Dna.display_base(b't'), b't');
        assert_eq!(MoleculeKind::Protein.display_base(b'T'), b'T');
    }
}

/// SnapGene **文档类别**名（头部 `type_of_sequence`）。
///
/// 实测（样例 + oracle）：`1` = DNA、`2` = 蛋白、`3` = DNA 比对、
/// `4` = 蛋白比对、`7` = RNA、`8` = 凝胶；其余值返回 `"unknown"`。
/// 类别 3/4/8 是 SnapGene 的比对编辑器 / 凝胶编辑器文档，本实现只能原样保存。
pub const fn document_class_name(v: u16) -> &'static str {
    match v {
        1 => "DNA",
        2 => "Protein",
        3 => "DNA alignment",
        4 => "Protein alignment",
        7 => "RNA",
        8 => "Gel",
        _ => "unknown",
    }
}

/// 该类别是否是我们能编辑/显示的分子（DNA/RNA/蛋白）。
pub const fn document_class_editable(v: u16) -> bool {
    matches!(v, 1 | 2 | 7)
}
