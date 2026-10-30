//! 最小 JSON 解析器（零依赖）。
//!
//! 只为读 SnapGene 的压缩块服务（block 34 的二级结构 JSON），但实现是
//! 通用的：对象/数组/字符串/数字/布尔/null，允许任意空白。

use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(fields) => fields.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n),
            _ => None,
        }
    }

    pub fn as_usize(&self) -> Option<usize> {
        match self {
            Json::Num(n) if *n >= 0.0 && n.fract() == 0.0 => Some(*n as usize),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Json::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonError {
    pub offset: usize,
    pub msg: String,
}

impl fmt::Display for JsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "JSON 非法（偏移 {}）：{}", self.offset, self.msg)
    }
}

impl std::error::Error for JsonError {}

pub fn parse(bytes: &[u8]) -> Result<Json, JsonError> {
    let mut p = Parser { s: bytes, i: 0 };
    p.ws();
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(p.err("末尾有多余内容"));
    }
    Ok(v)
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> Parser<'a> {
    fn err(&self, msg: impl Into<String>) -> JsonError {
        JsonError {
            offset: self.i,
            msg: msg.into(),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }

    fn ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.i += 1;
        }
    }

    fn value(&mut self) -> Result<Json, JsonError> {
        match self.peek() {
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b't') => self.literal(b"true", Json::Bool(true)),
            Some(b'f') => self.literal(b"false", Json::Bool(false)),
            Some(b'n') => self.literal(b"null", Json::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => Err(self.err("需要值")),
        }
    }

    fn literal(&mut self, word: &[u8], v: Json) -> Result<Json, JsonError> {
        if self.s.get(self.i..self.i + word.len()) == Some(word) {
            self.i += word.len();
            Ok(v)
        } else {
            Err(self.err("字面量不完整"))
        }
    }

    fn object(&mut self) -> Result<Json, JsonError> {
        self.i += 1; // {
        let mut fields = Vec::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Ok(Json::Obj(fields));
        }
        loop {
            self.ws();
            if self.peek() != Some(b'"') {
                return Err(self.err("对象键必须是字符串"));
            }
            let k = self.string()?;
            self.ws();
            if self.peek() != Some(b':') {
                return Err(self.err("缺少 ':'"));
            }
            self.i += 1;
            self.ws();
            let v = self.value()?;
            fields.push((k, v));
            self.ws();
            match self.peek() {
                Some(b',') => {
                    self.i += 1;
                }
                Some(b'}') => {
                    self.i += 1;
                    return Ok(Json::Obj(fields));
                }
                _ => return Err(self.err("对象缺少 ',' 或 '}'")),
            }
        }
    }

    fn array(&mut self) -> Result<Json, JsonError> {
        self.i += 1; // [
        let mut items = Vec::new();
        self.ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Ok(Json::Arr(items));
        }
        loop {
            self.ws();
            items.push(self.value()?);
            self.ws();
            match self.peek() {
                Some(b',') => {
                    self.i += 1;
                }
                Some(b']') => {
                    self.i += 1;
                    return Ok(Json::Arr(items));
                }
                _ => return Err(self.err("数组缺少 ',' 或 ']'")),
            }
        }
    }

    fn string(&mut self) -> Result<String, JsonError> {
        self.i += 1; // "
        let mut out = String::new();
        loop {
            match self.peek() {
                None => return Err(self.err("字符串未闭合")),
                Some(b'"') => {
                    self.i += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.i += 1;
                    let esc = self.peek().ok_or_else(|| self.err("转义未完成"))?;
                    self.i += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let cp = self.hex4()?;
                            // 代理对
                            if (0xD800..0xDC00).contains(&cp) {
                                if self.s.get(self.i..self.i + 2) == Some(b"\\u") {
                                    self.i += 2;
                                    let lo = self.hex4()?;
                                    let c = 0x10000
                                        + ((cp - 0xD800) << 10)
                                        + (lo.wrapping_sub(0xDC00) & 0x3FF);
                                    out.push(char::from_u32(c).unwrap_or('\u{FFFD}'));
                                } else {
                                    out.push('\u{FFFD}');
                                }
                            } else {
                                out.push(char::from_u32(cp).unwrap_or('\u{FFFD}'));
                            }
                        }
                        other => return Err(self.err(format!("未知转义 \\{}", other as char))),
                    }
                }
                Some(b) if b < 0x20 => return Err(self.err("字符串里有裸控制字符")),
                Some(_) => {
                    // UTF-8：按字节起首长度取一段
                    let start = self.i;
                    let b = self.s[start];
                    let len = if b < 0x80 {
                        1
                    } else if b >> 5 == 0b110 {
                        2
                    } else if b >> 4 == 0b1110 {
                        3
                    } else if b >> 3 == 0b11110 {
                        4
                    } else {
                        return Err(self.err("UTF-8 首字节非法"));
                    };
                    let end = start + len;
                    let piece = self
                        .s
                        .get(start..end)
                        .ok_or_else(|| self.err("UTF-8 截断"))?;
                    let s = std::str::from_utf8(piece).map_err(|_| self.err("UTF-8 非法"))?;
                    out.push_str(s);
                    self.i = end;
                }
            }
        }
    }

    fn hex4(&mut self) -> Result<u32, JsonError> {
        let chunk = self
            .s
            .get(self.i..self.i + 4)
            .ok_or_else(|| self.err("\\u 不完整"))?;
        let s = std::str::from_utf8(chunk).map_err(|_| self.err("\\u 非 ASCII"))?;
        let v = u32::from_str_radix(s, 16).map_err(|_| self.err("\\u 不是十六进制"))?;
        self.i += 4;
        Ok(v)
    }

    fn number(&mut self) -> Result<Json, JsonError> {
        let start = self.i;
        if self.peek() == Some(b'-') {
            self.i += 1;
        }
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.i += 1;
        }
        if self.peek() == Some(b'.') {
            self.i += 1;
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.i += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.i += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.i += 1;
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.i += 1;
            }
        }
        let text = std::str::from_utf8(&self.s[start..self.i]).map_err(|_| self.err("数字非 ASCII"))?;
        text.parse::<f64>()
            .map(Json::Num)
            .map_err(|_| self.err("数字格式非法"))
    }
}
