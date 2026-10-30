//! # hx-core — HelixBench 核心域模型
//!
//! 与 UI 和文件格式无关的最小共享类型：
//! - 分子类型（DNA/RNA/蛋白）与 `type_of_sequence` 映射；
//! - SnapGene 序列属性位（block 0/21/32 首字节）的语义封装；
//! - 序列记录 `SeqRecord`。
//!
//! 位定义来自对 SnapGene 8.2.2 自带 203 个 `.dna` 样例的实测，
//! 详见 `docs/SnapGene_dna_format_spec.md` §3。

#![forbid(unsafe_code)]

pub mod flags;
pub mod molecule;
pub mod protein;

pub use flags::{MoleculeFlags, Strandedness, Topology};
pub use molecule::{document_class_editable, document_class_name, MoleculeKind, SeqRecord};
pub use protein::Composition;
