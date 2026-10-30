//! 格式层错误类型（无外部依赖）。

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FormatError {
    /// 数据在 `offset` 处截断：需要 `needed` 字节，仅剩 `avail`。
    Truncated {
        offset: usize,
        needed: usize,
        avail: usize,
    },
    /// 块声明长度超出防御性上限。
    BlockTooLarge { offset: usize, len: u32 },
    /// 缺少 SnapGene 头部块（block 9 必须是首个块）。
    MissingHeader,
    /// 头部块内容不符（长度/魔术字）。
    BadHeader(String),
    /// block 1 压缩序列描述符/内容非法。
    BadCompressedSeq(String),
    /// block 34 二级结构块非法（XZ 解压或 JSON）。
    BadStructure(String),
    /// 编辑操作非法（越界、区间次序颠倒等）。
    BadEdit(String),
}

impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FormatError::Truncated {
                offset,
                needed,
                avail,
            } => write!(
                f,
                "数据截断：偏移 {offset} 处需要 {needed} 字节，仅剩 {avail}"
            ),
            FormatError::BlockTooLarge { offset, len } => {
                write!(f, "块长度异常：偏移 {offset} 处声明 {len} 字节")
            }
            FormatError::MissingHeader => write!(f, "缺少 SnapGene 头部块（block 9）"),
            FormatError::BadHeader(msg) => write!(f, "头部块格式不符：{msg}"),
            FormatError::BadCompressedSeq(msg) => write!(f, "block 1 压缩序列非法：{msg}"),
            FormatError::BadStructure(msg) => write!(f, "block 34 二级结构非法：{msg}"),
            FormatError::BadEdit(msg) => write!(f, "编辑非法：{msg}"),
        }
    }
}

impl std::error::Error for FormatError {}
