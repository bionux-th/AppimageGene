//! 端到端验证：用我们的匹配器 + 规范站点表**重建** block 3，与真 SnapGene
//! 写在该文件里的缓存**逐字节比对**。
//!
//! 这一条测试同时验证三件事（无需激活版 SnapGene）：
//!
//! 1. 酶匹配语义：IUPAC、允许重叠、环状跨原点命中；
//! 2. block 3 编码：站点表、u32/u16 双计数表、7 个**文件绝对偏移**、u16 LE 位置流；
//! 3. 分桶次序：`len≥9 / 8 / 7 / 6 / 5 / 4` 的 6 个边界与位置分组。
//!
//! 样例目录优先取 `HELIX_SAMPLES_DIR`；找不到时跳过（打印提示）。

use std::path::{Path, PathBuf};

const DEFAULT_SAMPLE_DIR: &str =
    "/home/mky/Software/snapgene_8.2.2_linux/data/opt/gslbiotech/snapgene/resources";

fn samples_root() -> Option<PathBuf> {
    let dir = std::env::var("HELIX_SAMPLES_DIR").unwrap_or_else(|_| DEFAULT_SAMPLE_DIR.to_string());
    let p = PathBuf::from(dir);
    p.is_dir().then_some(p)
}

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for e in rd.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect(&p, out);
        } else if matches!(
            p.extension().and_then(|s| s.to_str()),
            Some("dna" | "rna" | "prot")
        ) {
            out.push(p);
        }
    }
}

/// 首个不同字节的位置（用于诊断输出）。
fn first_diff(a: &[u8], b: &[u8]) -> String {
    let n = a.len().min(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            return format!("首个差异 @{i}（重建 {:#04x} vs 原文件 {:#04x}）", a[i], b[i]);
        }
    }
    format!("长度不同（重建 {} vs 原文件 {}）", a.len(), b.len())
}

#[test]
fn rebuilt_block3_matches_snapgene_cache() {
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录（设置 HELIX_SAMPLES_DIR 以启用）");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    files.sort();

    let mut checked = 0usize;
    let mut total_hits = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for p in &files {
        let bytes = std::fs::read(p).unwrap();
        let Ok(f) = hx_formats::DnaFile::parse(&bytes) else {
            continue;
        };
        let Some(b3) = f.block(3) else { continue };
        let Some(rec) = f.sequence() else { continue };
        let circular = rec.flags.has(hx_core::MoleculeFlags::CIRCULAR);
        let payload_offset = b3.offset as u64 + 5; // 跳过块头 [type:1][len:4]
        let rebuilt = hx_enzymes::block3::encode(&rec.sequence, circular, payload_offset);
        checked += 1;
        if let Ok(d) = hx_enzymes::block3::decode(&rebuilt) {
            total_hits += d.counts.iter().map(|c| *c as usize).sum::<usize>();
        }
        if rebuilt != b3.payload {
            let detail = first_diff(&rebuilt, &b3.payload);
            // 若两边都能解码，进一步指出第一个不一致的站点
            let hint = match (
                hx_enzymes::block3::decode(&rebuilt),
                hx_enzymes::block3::decode(&b3.payload),
            ) {
                (Ok(a), Ok(b)) => {
                    let sites: Vec<&str> = a.site_text.split(',').collect();
                    let mut msg = String::from("（结构相同但内容不同）");
                    if let Some(i) = (0..a.counts.len().min(b.counts.len()))
                        .find(|i| a.counts[*i] != b.counts[*i] || a.positions[*i] != b.positions[*i])
                    {
                        msg = format!(
                            "（站点 #{} {} 计数 {} vs {}；位置 重建{:?} 原{:?}）",
                            i,
                            sites.get(i).copied().unwrap_or("?"),
                            a.counts[i],
                            b.counts[i],
                            &a.positions[i][..a.positions[i].len().min(8)],
                            &b.positions[i][..b.positions[i].len().min(8)],
                        );
                    } else if a.offsets != b.offsets {
                        msg = format!("（偏移表 重建{:?} vs 原{:?}）", a.offsets, b.offsets);
                    }
                    msg
                }
                _ => "（结构不同）".to_string(),
            };
            failures.push(format!(
                "{} 环状={circular} {} {hint}",
                p.file_name().and_then(|s| s.to_str()).unwrap_or("?"),
                detail
            ));
        }
    }

    println!(
        "block 3 重建比对：{checked} 个文件，命中总计 {total_hits} 处，失败 {}",
        failures.len()
    );
    assert!(checked >= 190, "带 block 3 的样例只有 {checked} 个，偏少");
    if !failures.is_empty() {
        for f in failures.iter().take(20) {
            eprintln!("  ✗ {f}");
        }
        panic!("block 3 重建与原文件不一致：{} 个", failures.len());
    }
}
