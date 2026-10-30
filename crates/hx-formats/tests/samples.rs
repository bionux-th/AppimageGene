//! 黄金测试：对 SnapGene 8.2.2 自带 `sampleData` 全量做读→写往返（逐字节）。
//!
//! 样例目录优先取 `HELIX_SAMPLES_DIR`；否则用本机 SnapGene 安装包默认路径。
//! 找不到样例时跳过（打印提示，不失败），便于在无 SnapGene 的机器上运行。

use std::path::{Path, PathBuf};

const DEFAULT_SAMPLE_DIR: &str =
    "/home/mky/Software/snapgene_8.2.2_linux/data/opt/gslbiotech/snapgene/resources";

fn samples_root() -> Option<PathBuf> {
    let dir = std::env::var("HELIX_SAMPLES_DIR").unwrap_or_else(|_| DEFAULT_SAMPLE_DIR.to_string());
    let p = PathBuf::from(dir);
    if p.is_dir() {
        Some(p)
    } else {
        None
    }
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

fn find_by_name(root: &Path, name: &str) -> Option<PathBuf> {
    let mut all = Vec::new();
    collect(root, &mut all);
    all.into_iter()
        .find(|p| p.file_name().and_then(|s| s.to_str()) == Some(name))
}

#[test]
fn sample_files_roundtrip_byte_exact() {
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录（设置 HELIX_SAMPLES_DIR 以启用）");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    assert!(files.len() >= 200, "样例数量 {} 偏少，目录：{}", files.len(), root.display());

    let mut failed = Vec::new();
    let mut blocks_total = 0usize;
    for p in &files {
        let bytes = std::fs::read(p).unwrap();
        match hx_formats::DnaFile::parse(&bytes) {
            Ok(f) => {
                blocks_total += f.blocks.len();
                let out = f.to_bytes();
                if out != bytes {
                    let at = out
                        .iter()
                        .zip(bytes.iter())
                        .position(|(a, b)| a != b)
                        .unwrap_or(out.len().min(bytes.len()));
                    failed.push(format!(
                        "{}：首个差异 @{}（{} vs {} 字节）",
                        p.display(),
                        at,
                        out.len(),
                        bytes.len()
                    ));
                }
            }
            Err(e) => failed.push(format!("{}：解析失败 {e}", p.display())),
        }
    }
    assert!(
        failed.is_empty(),
        "往返失败 {}/{}：\n{}",
        failed.len(),
        files.len(),
        failed.join("\n")
    );
    eprintln!(
        "往返通过：{}/{} 文件，共 {} 个块",
        files.len(),
        files.len(),
        blocks_total
    );
}

/// block 1 快照解码验证（`2. mCherry BsaI.dna` 节点 0：711 bp mCherry）。
#[test]
fn block1_snapshot_mcherry_decodes() {
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录");
        return;
    };
    let Some(p) = find_by_name(&root, "2. mCherry BsaI.dna") else {
        eprintln!("跳过：样例文件缺失");
        return;
    };
    let f = hx_formats::DnaFile::parse(&std::fs::read(&p).unwrap()).unwrap();
    let node0 = f
        .blocks_of(11)
        .find(|b| b.payload.len() > 9 && b.payload[4] == 1)
        .expect("节点 0 快照（st=1）");
    let payload = &node0.payload;
    let cl = u32::from_be_bytes(payload[5..9].try_into().unwrap()) as usize;
    let block1_len = 4 + cl;

    let comp = hx_formats::block1::CompressedSeq::parse(&payload[5..5 + block1_len]).unwrap();
    assert_eq!(comp.uncompressed_len, 711);
    let seq = comp.to_sequence();
    assert_eq!(seq.len(), 711);
    // mCherry CDS 起始（大写，无小写区间）
    assert!(seq.starts_with(b"ATGGTGAGCAAGGGCGAGGAG"), "实际前缀：{:?}", &seq[..12]);
    // 保真：解析后再序列化必须与原字节一致
    assert_eq!(comp.to_payload(), &payload[5..5 + block1_len]);
}

/// 头部不变量：`.dna` = tos 1 / export 15 / import 19。
#[test]
fn sample_headers_are_canonical() {
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    for p in &files {
        let bytes = std::fs::read(p).unwrap();
        let f = hx_formats::DnaFile::parse(&bytes).unwrap();
        let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
        let expected: u16 = match ext {
            "dna" => 1,
            "prot" => 2,
            "rna" => 7,
            _ => continue,
        };
        assert_eq!(
            f.header.type_of_sequence,
            expected,
            "{} 的 type_of_sequence",
            p.display()
        );
        if ext == "dna" {
            assert_eq!(f.header.export_version, 15, "{}", p.display());
            assert_eq!(f.header.import_version, 19, "{}", p.display());
        }
    }
}
