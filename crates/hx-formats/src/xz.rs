//! XZ 容器 + LZMA2 解压（纯 Rust，零依赖）。
//!
//! SnapGene 把两类块压成 `.xz`：
//! - block 34（RNA 二级结构 JSON）——**必须**能读，否则 RNA 的
//!   `Secondary Structure` 视图只能显示“未实现”；
//! - block 29（历史树，LZMA 压缩）——可选。
//!
//! 实测样本（8.2.3）用的是标准 XZ 容器、单个块、单一 `--lzma2=dict=8MiB`
//! 过滤器、CRC64 校验（`xz --list -vv` 可见）。这里实现解码所需的完整子集：
//! 流头/块头/索引/流尾解析 + 标准 LZMA2 解码器（含状态、重复距离、
//! 距离槽位、字面量上下文、长度解码器）。
//!
//! 不校验 CRC（只做长度校验）：我们只读字节，校验失败也应由上层决定；
//! 读坏的数据会先撞上长度/索引错误。

use std::fmt;

/// XZ 解码错误。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum XzError {
    /// 输入不是 XZ（魔术字不符）。
    NotXz,
    /// 数据在 `offset` 处截断。
    Truncated { offset: usize, needed: usize },
    /// 结构不合法（块头/过滤器/索引）。
    Bad(String),
    /// LZMA2 数据不合法。
    BadLzma2(String),
}

impl fmt::Display for XzError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            XzError::NotXz => write!(f, "不是 XZ 流（魔术字不符）"),
            XzError::Truncated { offset, needed } => {
                write!(f, "XZ 数据截断：偏移 {offset} 处需要 {needed} 字节")
            }
            XzError::Bad(msg) => write!(f, "XZ 结构非法：{msg}"),
            XzError::BadLzma2(msg) => write!(f, "LZMA2 数据非法：{msg}"),
        }
    }
}

impl std::error::Error for XzError {}

const MAGIC: [u8; 6] = [0xFD, b'7', b'z', b'X', b'Z', 0x00];
const FOOTER_MAGIC: [u8; 2] = [b'Y', b'Z'];

/// 把一段完整的 XZ 流解压为原始字节。
pub fn decompress(src: &[u8]) -> Result<Vec<u8>, XzError> {
    let mut r = Reader { buf: src, pos: 0 };
    if r.take(6)? != MAGIC {
        return Err(XzError::NotXz);
    }
    // 流标志：第一个字节必须为 0，第二个字节低 4 位是校验类型。
    let flags0 = r.u8()?;
    let flags1 = r.u8()?;
    if flags0 != 0 || flags1 & 0xF0 != 0 {
        return Err(XzError::Bad("流标志保留位非 0".into()));
    }
    let check = Check::from_id(flags1 & 0x0F);
    r.skip(4)?; // 流头 CRC32（不校验）

    let mut out = Vec::new();
    loop {
        let hdr_start = r.pos;
        let header_size_byte = r.u8()?;
        if header_size_byte == 0 {
            // 索引指示符
            parse_index(&mut r)?;
            break;
        }
        let header_size = (header_size_byte as usize + 1) * 4;
        let block_flags = r.u8()?;
        if block_flags & 0x3C != 0 {
            return Err(XzError::Bad("块标志保留位非 0".into()));
        }
        let n_filters = (block_flags & 0x03) as usize + 1;
        let has_comp_size = block_flags & 0x40 != 0;
        let has_uncomp_size = block_flags & 0x80 != 0;
        let comp_size = if has_comp_size {
            Some(read_varint(&mut r)?)
        } else {
            None
        };
        let uncomp_size = if has_uncomp_size {
            Some(read_varint(&mut r)?)
        } else {
            None
        };
        let mut filters = Vec::with_capacity(n_filters);
        for _ in 0..n_filters {
            let id = read_varint(&mut r)?;
            let props_len = read_varint(&mut r)? as usize;
            let props = r.take(props_len)?.to_vec();
            filters.push((id, props));
        }
        // 头部填充 + CRC32（4 字节）。
        let consumed = r.pos - hdr_start;
        if consumed + 4 > header_size {
            return Err(XzError::Bad("块头长度不自洽".into()));
        }
        r.skip(header_size - consumed)?;

        if filters.len() != 1 || filters[0].0 != 0x21 {
            return Err(XzError::Bad(format!(
                "只支持单个 LZMA2 过滤器，实际 {:?}",
                filters.iter().map(|(id, _)| *id).collect::<Vec<_>>()
            )));
        }
        let body_start = r.pos;
        let before = out.len();
        lzma2_decode(&mut r, &mut out, comp_size)?;
        let produced = out.len() - before;
        if let Some(expect) = uncomp_size {
            if produced as u64 != expect {
                return Err(XzError::Bad(format!(
                    "块解压长度不符：声明 {expect}，实际 {produced}"
                )));
            }
        }
        if let Some(cs) = comp_size {
            if (r.pos - body_start) as u64 != cs {
                return Err(XzError::Bad(format!(
                    "块压缩长度不符：声明 {cs}，实际 {}",
                    r.pos - body_start
                )));
            }
        }
        // 块填充到 4 字节边界，然后跳过校验值。
        let pad = (4 - (r.pos % 4)) % 4;
        r.skip(pad)?;
        r.skip(check.size())?;
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Check {
    None,
    Crc32,
    Crc64,
    Sha256,
    /// 未知/保留的校验类型：按 fmt 里给出的长度处理。
    Other(usize),
}

impl Check {
    fn from_id(id: u8) -> Self {
        match id {
            0 => Check::None,
            1 => Check::Crc32,
            4 => Check::Crc64,
            10 => Check::Sha256,
            other => Check::Other((other as usize + 1).next_power_of_two()),
        }
    }
    fn size(self) -> usize {
        match self {
            Check::None => 0,
            Check::Crc32 => 4,
            Check::Crc64 => 8,
            Check::Sha256 => 32,
            Check::Other(n) => n,
        }
    }
}

fn parse_index(r: &mut Reader<'_>) -> Result<(), XzError> {
    let count = read_varint(r)?;
    for _ in 0..count {
        let _unpadded = read_varint(r)?;
        let _uncomp = read_varint(r)?;
    }
    // 索引填充到 4 字节边界，然后 CRC32+填充。
    let pad = (4 - (r.pos % 4)) % 4;
    r.skip(pad)?;
    r.skip(4)?;
    // 流尾：CRC32(4) 反向大小(4) 流标志(2) "YZ"。
    let tail = r.take(12)?;
    if tail[10..12] != FOOTER_MAGIC {
        return Err(XzError::Bad("流尾魔术字不符".into()));
    }
    if tail[8] != 0 {
        return Err(XzError::Bad("流尾标志非法".into()));
    }
    Ok(())
}

fn read_varint(r: &mut Reader<'_>) -> Result<u64, XzError> {
    let mut value = 0u64;
    for i in 0..9 {
        let b = r.u8()?;
        value |= ((b & 0x7F) as u64) << (7 * i);
        if b & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(XzError::Bad("varint 超过 9 字节".into()))
}

struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn u8(&mut self) -> Result<u8, XzError> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or(XzError::Truncated { offset: self.pos, needed: 1 })?;
        self.pos += 1;
        Ok(b)
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8], XzError> {
        let end = self.pos.checked_add(n).ok_or(XzError::Truncated { offset: self.pos, needed: n })?;
        if end > self.buf.len() {
            return Err(XzError::Truncated { offset: self.pos, needed: n });
        }
        let s = &self.buf[self.pos..end];
        self.pos = end;
        Ok(s)
    }
    fn skip(&mut self, n: usize) -> Result<(), XzError> {
        self.take(n).map(|_| ())
    }
}

// ---------------------------------------------------------------------------
// LZMA2
// ---------------------------------------------------------------------------

/// 概率模型总位数（LZMA 固定 11）。
const PROB_BITS: u32 = 11;
const PROB_INIT: u16 = 1 << (PROB_BITS - 1);
const MOVE_BITS: u32 = 5;

const NUM_STATES: usize = 12;
const NUM_POS_BITS_MAX: usize = 4;
const NUM_POS_STATES: usize = 1 << NUM_POS_BITS_MAX;
const NUM_LEN_TO_POS_STATES: usize = 4;
const NUM_ALIGN_BITS: usize = 4;
const START_POS_MODEL_INDEX: u32 = 4;
const END_POS_MODEL_INDEX: u32 = 14;
const NUM_FULL_DISTANCES: usize = 1 << 7;
const MATCH_MIN_LEN: u32 = 2;

/// LZMA2 解码：从 `r` 读取（`comp_size` 给定时限制在该长度内），写进 `out`。
fn lzma2_decode(
    r: &mut Reader<'_>,
    out: &mut Vec<u8>,
    comp_size: Option<u64>,
) -> Result<(), XzError> {
    let start = r.pos;
    let limit = match comp_size {
        Some(n) => start + n as usize,
        None => r.buf.len(),
    };
    let mut dec = LzmaDecoder::new();
    // 字典起点：遇到“重置字典”时前移，防止匹配跨过重置点。
    let mut dict_start = 0usize;
    loop {
        if r.pos >= limit {
            return Err(XzError::BadLzma2("缺少结束标记（0x00）".into()));
        }
        let control = r.u8()?;
        if control == 0 {
            break;
        }
        if control < 0x80 {
            // 未压缩块：0x01 = 重置字典，0x02 = 不重置。
            if control != 1 && control != 2 {
                return Err(XzError::BadLzma2(format!("未压缩块控制字节 0x{control:02X}")));
            }
            let hi = r.u8()? as usize;
            let lo = r.u8()? as usize;
            let size = ((hi << 8) | lo) + 1;
            if control == 1 {
                // 重置字典：旧数据不能再被引用；LZMA 状态同时归零。
                dict_start = out.len();
            }
            // 规范：未压缩块之后 LZMA 状态（state/reps）必须重置；概率模型保留。
            dec.reset_state();
            let bytes = r.take(size)?;
            out.extend_from_slice(bytes);
            continue;
        }
        let reset = (control >> 5) & 0x03;
        let unpack_hi = (control & 0x1F) as usize;
        let b1 = r.u8()? as usize;
        let b2 = r.u8()? as usize;
        let unpack_size = ((unpack_hi << 16) | (b1 << 8) | b2) + 1;
        let c1 = r.u8()? as usize;
        let c2 = r.u8()? as usize;
        let comp_len = ((c1 << 8) | c2) + 1;
        // 实测顺序：控制字节 → 解压长度(2) → 压缩长度(2) → [新属性字节]。
        if reset >= 2 {
            let props = r.u8()?;
            if props > (4 * 5 + 4) * 9 + 8 {
                return Err(XzError::BadLzma2(format!("LZMA 属性字节 0x{props:02X}")));
            }
            dec.set_props(props);
            dec.reset_probs();
        }
        if reset >= 1 {
            dec.reset_state();
        }
        if reset == 3 {
            dict_start = out.len();
        }
        let payload = r.take(comp_len)?;
        dec.decode_chunk(payload, unpack_size, out, dict_start)?;
        if r.pos > limit {
            return Err(XzError::BadLzma2("块压缩数据越过声明的压缩长度".into()));
        }
    }
    Ok(())
}

struct LzmaDecoder {
    lc: u32,
    lp: u32,
    pb: u32,
    state: usize,
    reps: [u32; 4],
    is_match: [u16; NUM_STATES * NUM_POS_STATES],
    is_rep: [u16; NUM_STATES],
    is_rep_g0: [u16; NUM_STATES],
    is_rep_g1: [u16; NUM_STATES],
    is_rep_g2: [u16; NUM_STATES],
    is_rep0_long: [u16; NUM_STATES * NUM_POS_STATES],
    pos_slot: [[u16; 64]; NUM_LEN_TO_POS_STATES],
    pos_decoders: [u16; NUM_FULL_DISTANCES],
    align: [u16; 1 << NUM_ALIGN_BITS],
    len: LenDecoder,
    rep_len: LenDecoder,
    literals: Vec<u16>,
}

impl LzmaDecoder {
    fn new() -> Self {
        let mut d = Self {
            lc: 3,
            lp: 0,
            pb: 2,
            state: 0,
            reps: [0; 4],
            is_match: [PROB_INIT; NUM_STATES * NUM_POS_STATES],
            is_rep: [PROB_INIT; NUM_STATES],
            is_rep_g0: [PROB_INIT; NUM_STATES],
            is_rep_g1: [PROB_INIT; NUM_STATES],
            is_rep_g2: [PROB_INIT; NUM_STATES],
            is_rep0_long: [PROB_INIT; NUM_STATES * NUM_POS_STATES],
            pos_slot: [[PROB_INIT; 64]; NUM_LEN_TO_POS_STATES],
            pos_decoders: [PROB_INIT; NUM_FULL_DISTANCES],
            align: [PROB_INIT; 1 << NUM_ALIGN_BITS],
            len: LenDecoder::new(),
            rep_len: LenDecoder::new(),
            literals: vec![PROB_INIT; 0x300],
        };
        d.reset_probs();
        d
    }

    fn set_props(&mut self, props: u8) {
        let mut p = props as u32;
        self.lc = p % 9;
        p /= 9;
        self.lp = p % 5;
        p /= 5;
        self.pb = p;
    }

    fn reset_state(&mut self) {
        self.state = 0;
        self.reps = [0; 4];
    }

    fn reset_probs(&mut self) {
        self.is_match.fill(PROB_INIT);
        self.is_rep.fill(PROB_INIT);
        self.is_rep_g0.fill(PROB_INIT);
        self.is_rep_g1.fill(PROB_INIT);
        self.is_rep_g2.fill(PROB_INIT);
        self.is_rep0_long.fill(PROB_INIT);
        for row in self.pos_slot.iter_mut() {
            row.fill(PROB_INIT);
        }
        self.pos_decoders.fill(PROB_INIT);
        self.align.fill(PROB_INIT);
        self.len.reset();
        self.rep_len.reset();
        let n = 0x300usize << (self.lc + self.lp);
        self.literals.clear();
        self.literals.resize(n, PROB_INIT);
    }

    fn decode_chunk(
        &mut self,
        data: &[u8],
        unpack_size: usize,
        out: &mut Vec<u8>,
        dict_start: usize,
    ) -> Result<(), XzError> {
        let mut rc = RangeDecoder::new(data)?;
        let target = out.len() + unpack_size;
        // 概率数组尺寸随 lc+lp 变化，确保够用。
        let need = 0x300usize << (self.lc + self.lp);
        if self.literals.len() != need {
            self.literals.clear();
            self.literals.resize(need, PROB_INIT);
        }
        while out.len() < target {
            let pos = out.len();
            let pos_state = pos & ((1 << self.pb) - 1);
            let is_match = rc.decode_bit(&mut self.is_match[(self.state << NUM_POS_BITS_MAX) + pos_state]);
            if is_match == 0 {
                // 字面量
                let prev = if pos == 0 { 0u8 } else { out[pos - 1] } as usize;
                let lit_state = ((pos & ((1 << self.lp) - 1)) << self.lc)
                    + ((prev >> (8 - self.lc)) & ((1 << self.lc) - 1));
                let base = 0x300 * lit_state;
                let mut symbol: u32 = 1;
                if self.state < 7 {
                    while symbol < 0x100 {
                        let bit = rc.decode_bit(&mut self.literals[base + symbol as usize]);
                        symbol = (symbol << 1) | bit;
                    }
                } else {
                    let dist = self.reps[0] as usize + 1;
                    if dist > pos - dict_start {
                        return Err(XzError::BadLzma2("字面量匹配字节距离越界".into()));
                    }
                    let mut match_byte = out[pos - dist] as u32;
                    let mut offs: u32 = 0x100;
                    while symbol < 0x100 {
                        match_byte <<= 1;
                        let match_bit = match_byte & offs;
                        let bit = rc.decode_bit(
                            &mut self.literals[base + (offs + match_bit + symbol) as usize],
                        );
                        symbol = (symbol << 1) | bit;
                        // offs 表示「此前每一位都与 matchByte 相同」；一旦失配就永久归零。
                        // 注意不能写成 offs &= !(match_byte ^ symbol)：解码器的 symbol
                        // 未左移到 bit8，与编码器的 symbol 对齐方式不同。
                        if bit == 0 {
                            offs &= !match_bit;
                        } else {
                            offs &= match_bit;
                        }
                    }
                }
                out.push((symbol & 0xFF) as u8);
                self.state = if self.state < 4 {
                    0
                } else if self.state < 10 {
                    self.state - 3
                } else {
                    self.state - 6
                };
                continue;
            }

            let len;
            if rc.decode_bit(&mut self.is_rep[self.state]) != 0 {
                if rc.decode_bit(&mut self.is_rep_g0[self.state]) == 0 {
                    if rc.decode_bit(
                        &mut self.is_rep0_long[(self.state << NUM_POS_BITS_MAX) + pos_state],
                    ) == 0
                    {
                        // 短重复：复制 1 字节
                        self.state = if self.state < 7 { 9 } else { 11 };
                        let dist = self.reps[0] as usize + 1;
                        if dist > pos - dict_start {
                            return Err(XzError::BadLzma2("短重复距离越界".into()));
                        }
                        let b = out[pos - dist];
                        out.push(b);
                        continue;
                    }
                } else {
                    let dist;
                    if rc.decode_bit(&mut self.is_rep_g1[self.state]) == 0 {
                        dist = self.reps[1];
                    } else {
                        if rc.decode_bit(&mut self.is_rep_g2[self.state]) == 0 {
                            dist = self.reps[2];
                        } else {
                            dist = self.reps[3];
                            self.reps[3] = self.reps[2];
                        }
                        self.reps[2] = self.reps[1];
                    }
                    self.reps[1] = self.reps[0];
                    self.reps[0] = dist;
                }
                len = self.decode_len(&mut rc, pos_state, true) + MATCH_MIN_LEN;
                self.state = if self.state < 7 { 8 } else { 11 };
            } else {
                self.reps[3] = self.reps[2];
                self.reps[2] = self.reps[1];
                self.reps[1] = self.reps[0];
                let l = self.decode_len(&mut rc, pos_state, false);
                self.state = if self.state < 7 { 7 } else { 10 };
                self.reps[0] = self.decode_distance(&mut rc, l)?;
                len = l + MATCH_MIN_LEN;
            }
            let dist = self.reps[0] as usize + 1;
            if dist > pos - dict_start {
                return Err(XzError::BadLzma2(format!(
                    "匹配距离 {dist} 越过字典起点（已输出 {}）",
                    pos - dict_start
                )));
            }
            if pos + len as usize > target {
                return Err(XzError::BadLzma2("匹配越过块解压长度".into()));
            }
            for _ in 0..len {
                let b = out[out.len() - dist];
                out.push(b);
            }
        }
        Ok(())
    }

    fn decode_len(&mut self, rc: &mut RangeDecoder<'_>, pos_state: usize, rep: bool) -> u32 {
        let d = if rep { &mut self.rep_len } else { &mut self.len };
        let c0 = rc.decode_bit(&mut d.choice);
        let v = if c0 == 0 {
            rc.decode_bit_tree(&mut d.low[pos_state], 3)
        } else {
            let c1 = rc.decode_bit(&mut d.choice2);
            if c1 == 0 {
                8 + rc.decode_bit_tree(&mut d.mid[pos_state], 3)
            } else {
                16 + rc.decode_bit_tree(&mut d.high, 8)
            }
        };
        v
    }

    fn decode_distance(&mut self, rc: &mut RangeDecoder<'_>, len: u32) -> Result<u32, XzError> {
        let len_state = (len as usize).min(NUM_LEN_TO_POS_STATES - 1);
        let pos_slot = rc.decode_bit_tree(&mut self.pos_slot[len_state], 6);
        if pos_slot < START_POS_MODEL_INDEX {
            return Ok(pos_slot);
        }
        let num_direct = (pos_slot >> 1) - 1;
        let mut dist = (2 | (pos_slot & 1)) << num_direct;
        if pos_slot < END_POS_MODEL_INDEX {
            // 7-Zip SDK 把这段逆序位树平铺在 `SpecPos` 起始的 128 项里：
            // 首地址 = 距离基数 `dist`，树内下标从 1 开始（与 SDK 的 distance++ 一致）。
            let idx = dist as usize;
            dist += rc.decode_reverse_bit_tree(&mut self.pos_decoders[idx..], num_direct);
        } else {
            dist += rc.decode_direct_bits(num_direct - 4)? << 4;
            dist += rc.decode_reverse_bit_tree(&mut self.align, NUM_ALIGN_BITS as u32);
        }
        Ok(dist)
    }
}

#[derive(Clone)]
struct LenDecoder {
    choice: u16,
    choice2: u16,
    low: [[u16; 8]; NUM_POS_STATES],
    mid: [[u16; 8]; NUM_POS_STATES],
    high: [u16; 256],
}

impl LenDecoder {
    fn new() -> Self {
        Self {
            choice: PROB_INIT,
            choice2: PROB_INIT,
            low: [[PROB_INIT; 8]; NUM_POS_STATES],
            mid: [[PROB_INIT; 8]; NUM_POS_STATES],
            high: [PROB_INIT; 256],
        }
    }
    fn reset(&mut self) {
        self.choice = PROB_INIT;
        self.choice2 = PROB_INIT;
        for row in self.low.iter_mut() {
            row.fill(PROB_INIT);
        }
        for row in self.mid.iter_mut() {
            row.fill(PROB_INIT);
        }
        self.high.fill(PROB_INIT);
    }
}

struct RangeDecoder<'a> {
    data: &'a [u8],
    pos: usize,
    range: u32,
    code: u32,
}

impl<'a> RangeDecoder<'a> {
    fn new(data: &'a [u8]) -> Result<Self, XzError> {
        if data.len() < 5 || data[0] != 0 {
            return Err(XzError::BadLzma2("range coder 头部非法".into()));
        }
        let code = u32::from_be_bytes([data[1], data[2], data[3], data[4]]);
        Ok(Self {
            data,
            pos: 5,
            range: 0xFFFF_FFFF,
            code,
        })
    }

    fn next(&mut self) -> u8 {
        let b = self.data.get(self.pos).copied().unwrap_or(0);
        self.pos += 1;
        b
    }

    fn normalize(&mut self) {
        if self.range < (1 << 24) {
            self.range <<= 8;
            self.code = (self.code << 8) | self.next() as u32;
        }
    }

    fn decode_bit(&mut self, prob: &mut u16) -> u32 {
        let bound = (self.range >> PROB_BITS) * (*prob as u32);
        let bit;
        if self.code < bound {
            self.range = bound;
            let delta = (1u32 << PROB_BITS) - *prob as u32;
            *prob += (delta >> MOVE_BITS) as u16;
            bit = 0;
        } else {
            self.range -= bound;
            self.code -= bound;
            *prob -= *prob >> MOVE_BITS;
            bit = 1;
        }
        self.normalize();
        bit
    }

    fn decode_bit_tree(&mut self, probs: &mut [u16], num_bits: u32) -> u32 {
        let mut m: u32 = 1;
        for _ in 0..num_bits {
            let bit = self.decode_bit(&mut probs[m as usize]);
            m = (m << 1) | bit;
        }
        m - (1 << num_bits)
    }

    fn decode_reverse_bit_tree(&mut self, probs: &mut [u16], num_bits: u32) -> u32 {
        let mut m: u32 = 1;
        let mut symbol = 0u32;
        for i in 0..num_bits {
            let bit = self.decode_bit(&mut probs[m as usize]);
            m = (m << 1) | bit;
            symbol |= bit << i;
        }
        symbol
    }

    fn decode_direct_bits(&mut self, num_bits: u32) -> Result<u32, XzError> {
        let mut res = 0u32;
        for _ in 0..num_bits {
            self.range >>= 1;
            self.code = self.code.wrapping_sub(self.range);
            let t = 0u32.wrapping_sub(self.code >> 31);
            self.code = self.code.wrapping_add(self.range & t);
            res = (res << 1).wrapping_add(t.wrapping_add(1));
            self.normalize();
        }
        Ok(res)
    }
}

// ---------------------------------------------------------------------------
// 写出：只做「原样存」的 XZ 流
// ---------------------------------------------------------------------------

/// 把一段字节封成 XZ 流（LZMA2 非压缩块 + CRC64 校验）。
///
/// 用途只有一处：`block 34`（RNA 二级结构）必须是 XZ 压缩的 JSON，而 SnapGene
/// 自己的流用的是 CRC64。压不压缩对 SnapGene 没有影响——它读的是标准 XZ 流，
/// 所以这里直接写 LZMA2 的 `0x01/0x02`（原样拷贝）块，省掉一个 LZMA 匹配器，
/// 也不会引入任何外部依赖。几 KB 的结构 JSON 体积翻几倍也无所谓。
pub fn compress_store(src: &[u8]) -> Vec<u8> {
    // LZMA2 的「原样块」最大 64 KiB。
    const CHUNK: usize = 1 << 16;
    let mut out = Vec::with_capacity(src.len() + 64);

    // 流头：magic + 流标志（0x04 = CRC64）+ 标志的 CRC32。
    out.extend_from_slice(&MAGIC);
    let stream_flags = [0u8, 0x04];
    out.extend_from_slice(&stream_flags);
    out.extend_from_slice(&crc32(&stream_flags).to_le_bytes());

    // 块头（12 字节）：1 个 LZMA2 过滤器，字典 16 MiB，不带块内尺寸。
    let mut hdr = vec![0x02u8, 0x00, 0x21, 0x01, 0x18, 0, 0, 0];
    hdr.extend_from_slice(&crc32(&hdr).to_le_bytes());
    out.extend_from_slice(&hdr);

    // LZMA2 数据：第一个块 0x01（重置字典），其后 0x02；最后是 0x00 结束标记。
    let mut data = Vec::with_capacity(src.len() + src.len() / CHUNK * 3 + 1);
    for (i, chunk) in src.chunks(CHUNK).enumerate() {
        data.push(if i == 0 { 0x01 } else { 0x02 });
        let n = chunk.len() - 1;
        data.push((n >> 8) as u8);
        data.push((n & 0xFF) as u8);
        data.extend_from_slice(chunk);
    }
    data.push(0x00);

    // 未填充大小 = 块头 + 数据 + 校验（不含填充，XZ 规范的定义）。
    let unpadded = hdr.len() + data.len() + 8;
    out.extend_from_slice(&data);
    out.extend(std::iter::repeat_n(0u8, (4 - data.len() % 4) % 4));
    out.extend_from_slice(&crc64(src).to_le_bytes());

    // 索引：1 条记录（未填充大小、解压后大小）。
    let mut idx = vec![0x00u8, 0x01];
    put_varint(unpadded, &mut idx);
    put_varint(src.len(), &mut idx);
    while idx.len() % 4 != 0 {
        idx.push(0);
    }
    idx.extend_from_slice(&crc32(&idx).to_le_bytes());
    out.extend_from_slice(&idx);

    // 流尾：索引大小（以 4 字节为单位 − 1）+ 流标志 + 二者的 CRC32。
    let mut foot = Vec::with_capacity(6);
    foot.extend_from_slice(&(((idx.len() / 4) - 1) as u32).to_le_bytes());
    foot.extend_from_slice(&stream_flags);
    out.extend_from_slice(&crc32(&foot).to_le_bytes());
    out.extend_from_slice(&foot);
    out.extend_from_slice(&FOOTER_MAGIC);
    out
}

/// XZ 的可变长整数：小端 7 位一组，最高位是「还有后续」。
fn put_varint(mut v: usize, out: &mut Vec<u8>) {
    loop {
        let b = (v & 0x7F) as u8;
        v >>= 7;
        if v == 0 {
            out.push(b);
            return;
        }
        out.push(b | 0x80);
    }
}

/// CRC32（IEEE，与 zlib / XZ 头尾用的一致）。
fn crc32(data: &[u8]) -> u32 {
    let mut c = 0xFFFF_FFFFu32;
    for &b in data {
        c ^= b as u32;
        for _ in 0..8 {
            let m = (c & 1).wrapping_neg();
            c = (c >> 1) ^ (0xEDB8_8320 & m);
        }
    }
    !c
}

/// CRC64（ECMA-182，XZ 的 `Check ID = 4`）。
fn crc64(data: &[u8]) -> u64 {
    let mut c = 0xFFFF_FFFF_FFFF_FFFFu64;
    for &b in data {
        c ^= b as u64;
        for _ in 0..8 {
            let m = (c & 1).wrapping_neg();
            c = (c >> 1) ^ (0xC96C_5795_D787_0F42 & m);
        }
    }
    !c
}

#[cfg(test)]
mod write_tests {
    use super::*;

    #[test]
    fn store_roundtrip_small_and_large() {
        for src in [
            &b""[..],
            &b"hello"[..],
            &vec![b'x'; 70_000][..],
        ] {
            let xz = compress_store(src);
            let back = decompress(&xz).expect("自己写的 XZ 流要能解回来");
            assert_eq!(back, src);
        }
    }

    #[test]
    fn store_roundtrip_block34_fixture() {
        let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/sg_trna_block34.xz");
        let raw = decompress(&std::fs::read(&p).unwrap()).unwrap();
        let xz = compress_store(&raw);
        assert_eq!(decompress(&xz).unwrap(), raw);
    }

    #[test]
    fn crc_known_values() {
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc64(b"123456789"), 0x995D_C9BB_DF19_39FA);
    }
}
