//! # hx-algo — 编辑 / 历史 / 工具算法（占位）
//!
//! 规划内容（V1）：
//! - 编辑事务模型与 undo/redo（与 block 7/11 历史写出一致化）；
//! - block 3 酶位点缓存重建（依赖 [`hx_enzymes`] 匹配 + 计数/位置流布局）；
//! - 序列变换（反向互补、大小写、环状旋转）等。

#![forbid(unsafe_code)]

pub mod cloning;
pub mod edit_session;
pub mod export;
pub mod orf;
pub mod primer;
pub mod translate;

pub use edit_session::{complement, edit_label, parse_seq_text, reverse_complement, EditSession};
pub use orf::{find_orfs, Orf};
pub use primer::{design_pair, find_bindings, gc_percent, tm, DesignOptions, PrimerPair};
pub use translate::{codon_aa, one_letter, three_letter, translate, translate_orf};
pub mod hybridization;
