//! # hx-formats — SnapGene 文件格式层
//!
//! 依据 `docs/SnapGene_dna_format_spec.md`（8.2.2 实测版）实现：
//!
//! - [`block`]：通用 TLV 块层（`[type: u8][len: u32 BE][payload]`）；
//! - [`dna`]：`.dna`/`.rna`/`.prot` 顶层模型（头部块 9 + 块序列）；
//! - [`block1`]：压缩 DNA 序列块（2-bit/4-bit/N 段 + 小写区间）。
//!
//! 设计原则：**载荷一律原样保存**，未编辑文件的读→写必须字节一致
//! （黄金测试见 `tests/samples.rs`）。

#![forbid(unsafe_code)]

pub mod annot;
pub mod block;
pub mod block1;
pub mod dna;
pub mod error;
pub mod feature_edit;
pub mod history;
pub mod import;
pub mod json;
pub mod notes;
pub mod primer_edit;
pub mod rna_struct;
pub mod util;
pub mod write;
pub mod xz;

pub use block::Block;
pub use dna::{DnaFile, Header};
pub use error::FormatError;
pub mod alignment;
