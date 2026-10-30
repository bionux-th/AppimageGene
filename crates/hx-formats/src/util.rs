//! 写文件用的小工具：UUID v4、UTC 时间戳、XML 转义。
//!
//! 为了保持工作区“零外部依赖”，这些都不引第三方 crate：
//! - UUID v4 用 `/dev/urandom`（不可用时退化为时间+进程扰动）。
//! - 时间换算用 Howard Hinnant 的 `civil_from_days` 算法。

use std::io::Read;
use std::time::{SystemTime, UNIX_EPOCH};

/// 生成 UUID v4 字符串（`xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx`）。
pub fn uuid_v4() -> String {
    let mut b = [0u8; 16];
    if std::fs::File::open("/dev/urandom")
        .and_then(|mut f| f.read_exact(&mut b))
        .is_err()
    {
        // 退化路径：用系统时间 + 进程号做简单混淆（仅用于元数据，不作为安全随机源）。
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0);
        let pid = std::process::id() as u64;
        let mut x = nanos ^ (pid << 32) ^ 0x9E37_79B9_7F4A_7C15;
        for chunk in b.chunks_mut(8) {
            // xorshift64*
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            let v = x.wrapping_mul(0x2545_F491_4F6C_DD1D);
            let bytes = v.to_le_bytes();
            chunk.copy_from_slice(&bytes[..chunk.len()]);
        }
    }
    b[6] = (b[6] & 0x0f) | 0x40; // version 4
    b[8] = (b[8] & 0x3f) | 0x80; // variant
    let h = |r: &[u8]| r.iter().map(|b| format!("{b:02x}")).collect::<String>();
    format!(
        "{}-{}-{}-{}-{}",
        h(&b[0..4]),
        h(&b[4..6]),
        h(&b[6..8]),
        h(&b[8..10]),
        h(&b[10..16])
    )
}

/// 当前 UTC 时间：`(year, month, day, hour, minute, second)`。
pub fn utc_now() -> (i64, u32, u32, u32, u32, u32) {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    utc_from_unix(secs)
}

/// Unix 秒 → UTC 年月日时分秒（Hinnant `civil_from_days`）。
pub fn utc_from_unix(secs: i64) -> (i64, u32, u32, u32, u32, u32) {
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    let y = if m <= 2 { y + 1 } else { y };
    (
        y,
        m as u32,
        d as u32,
        (rem / 3600) as u32,
        ((rem % 3600) / 60) as u32,
        (rem % 60) as u32,
    )
}

/// SnapGene `<Created>/<LastModified>` 的日期文本：`YYYY.M.D`（无前导零）。
pub fn snapgene_date((y, m, d, _, _, _): (i64, u32, u32, u32, u32, u32)) -> String {
    format!("{y}.{m}.{d}")
}

/// SnapGene `<Created>/<LastModified>` 的 `UTC` 属性值：`H:M:S`（无前导零）。
pub fn snapgene_utc_attr((_, _, _, h, mi, s): (i64, u32, u32, u32, u32, u32)) -> String {
    format!("{h}:{mi}:{s}")
}

/// XML 文本/属性转义（`&`、`<`、`>`、`"`）。
pub fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            _ => out.push(c),
        }
    }
    out
}

/// 纯文本 → SnapGene 的“富文本”HTML 包装（`Description`/`Comments` 用）。
pub fn html_wrap(text: &str) -> String {
    format!(
        "&lt;html&gt;&lt;body&gt;{}&lt;/body&gt;&lt;/html&gt;",
        xml_escape(text)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uuid_shape() {
        let u = uuid_v4();
        assert_eq!(u.len(), 36);
        let parts: Vec<&str> = u.split('-').collect();
        assert_eq!(parts.iter().map(|p| p.len()).collect::<Vec<_>>(), vec![8, 4, 4, 4, 12]);
        assert!(u.chars().all(|c| c.is_ascii_hexdigit() || c == '-'));
    }

    #[test]
    fn unix_epoch_conversion() {
        assert_eq!(utc_from_unix(0), (1970, 1, 1, 0, 0, 0));
        assert_eq!(utc_from_unix(1_000_000_000), (2001, 9, 9, 1, 46, 40));
        // 2026-10-16 00:00:00 UTC
        assert_eq!(utc_from_unix(1_792_108_800), (2026, 10, 16, 0, 0, 0));
        // 2026-10-15 21:20:00 UTC
        assert_eq!(utc_from_unix(1_792_099_200), (2026, 10, 15, 21, 20, 0));
    }

    #[test]
    fn date_and_utc_text() {
        let t = utc_from_unix(1_792_108_800);
        assert_eq!(snapgene_date(t), "2026.10.16");
        assert_eq!(snapgene_utc_attr(t), "0:0:0");
    }

    #[test]
    fn escaping() {
        assert_eq!(xml_escape("a&b<c>\"d\""), "a&amp;b&lt;c&gt;&quot;d&quot;");
    }
}
