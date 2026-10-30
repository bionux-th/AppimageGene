//! SnapGene 序列属性位（`block 0`/`21`/`32` 载荷首字节）。
//!
//! 实测位义（8.2.2 样例 203 文件全量验证）：
//!
//! | bit | 掩码 | 含义 |
//! |---|---|---|
//! | 0 | 0x01 | 环状（0 = 线性） |
//! | 1 | 0x02 | 双链（0 = 单链） |
//! | 2 | 0x04 | Dam 甲基化 |
//! | 3 | 0x08 | Dcm 甲基化 |
//! | 4 | 0x10 | EcoKI 甲基化 |

/// 分子拓扑。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Topology {
    Linear,
    Circular,
}

impl Topology {
    pub const fn is_circular(self) -> bool {
        matches!(self, Topology::Circular)
    }
}

/// 链性。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strandedness {
    Single,
    Double,
}

/// 序列属性字节。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub struct MoleculeFlags(pub u8);

impl MoleculeFlags {
    pub const CIRCULAR: u8 = 0x01;
    pub const DOUBLE_STRANDED: u8 = 0x02;
    pub const DAM_METHYLATED: u8 = 0x04;
    pub const DCM_METHYLATED: u8 = 0x08;
    pub const ECOKI_METHYLATED: u8 = 0x10;

    /// 大肠杆菌常规制备质粒：环状 + 双链 + 三种甲基化（0x1f）。
    pub const ECOLI_PLASMID_CIRCULAR: Self = Self(0x1f);
    /// 大肠杆菌常规制备质粒的线性化形式（0x1e）。
    pub const ECOLI_PLASMID_LINEAR: Self = Self(0x1e);

    #[inline]
    pub const fn from_byte(b: u8) -> Self {
        Self(b)
    }

    #[inline]
    pub const fn to_byte(self) -> u8 {
        self.0
    }

    #[inline]
    pub const fn has(self, bit: u8) -> bool {
        self.0 & bit != 0
    }

    #[inline]
    pub const fn with(self, bit: u8, on: bool) -> Self {
        if on {
            Self(self.0 | bit)
        } else {
            Self(self.0 & !bit)
        }
    }

    #[inline]
    pub const fn topology(self) -> Topology {
        if self.has(Self::CIRCULAR) {
            Topology::Circular
        } else {
            Topology::Linear
        }
    }

    #[inline]
    pub const fn strandedness(self) -> Strandedness {
        if self.has(Self::DOUBLE_STRANDED) {
            Strandedness::Double
        } else {
            Strandedness::Single
        }
    }

    #[inline]
    pub const fn with_topology(self, t: Topology) -> Self {
        self.with(Self::CIRCULAR, matches!(t, Topology::Circular))
    }

    #[inline]
    pub const fn with_strandedness(self, s: Strandedness) -> Self {
        self.with(Self::DOUBLE_STRANDED, matches!(s, Strandedness::Double))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typical_plasmid_flags() {
        let f = MoleculeFlags::ECOLI_PLASMID_CIRCULAR;
        assert_eq!(f.to_byte(), 0x1f);
        assert_eq!(f.topology(), Topology::Circular);
        assert_eq!(f.strandedness(), Strandedness::Double);
        assert!(f.has(MoleculeFlags::DAM_METHYLATED));
        assert!(f.has(MoleculeFlags::DCM_METHYLATED));
        assert!(f.has(MoleculeFlags::ECOKI_METHYLATED));
    }

    #[test]
    fn single_stranded_ssdna_flag_zero() {
        // ssDNA.dna 实测 0x00：线性 + 单链 + 无甲基化
        let f = MoleculeFlags::from_byte(0x00);
        assert_eq!(f.topology(), Topology::Linear);
        assert_eq!(f.strandedness(), Strandedness::Single);
        assert!(!f.has(MoleculeFlags::DAM_METHYLATED));
    }

    #[test]
    fn setters() {
        let f = MoleculeFlags::from_byte(0x02).with_topology(Topology::Circular);
        assert_eq!(f.to_byte(), 0x03);
        assert_eq!(f.strandedness(), Strandedness::Double);
        let g = f.with_strandedness(Strandedness::Single);
        assert_eq!(g.to_byte(), 0x01);
    }
}
