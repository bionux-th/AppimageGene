//! XZ/LZMA2 解码器测试：固定样本（含 SnapGene 真实的 block 34 载荷）。

use hx_formats::xz;

fn fixture(name: &str) -> Vec<u8> {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    std::fs::read(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()))
}

#[test]
fn decodes_simple_text() {
    let got = xz::decompress(&fixture("hello.txt.xz")).unwrap();
    assert_eq!(got, fixture("hello.txt"));
}

#[test]
fn decodes_incompressible_uncompressed_chunks() {
    // preset 0 对这段伪随机数据会退化成未压缩块（控制字节 0x01/0x02 路径）。
    let got = xz::decompress(&fixture("random.bin.xz")).unwrap();
    assert_eq!(got, fixture("random.bin"));
}

#[test]
fn decodes_snapgene_block34() {
    // SnapGene 8.2.3 的 tRNA 样例：XZ + LZMA2 dict=8MiB + CRC64。
    let got = xz::decompress(&fixture("sg_trna_block34.xz")).unwrap();
    let want = fixture("sg_trna_block34.json");
    assert_eq!(got.len(), want.len());
    assert_eq!(got, want);
    assert!(got.starts_with(b"{\"bondedPairs\":[[0,67],"));
}

#[test]
fn rejects_non_xz() {
    assert!(matches!(
        xz::decompress(b"not an xz stream"),
        Err(xz::XzError::NotXz)
    ));
}
