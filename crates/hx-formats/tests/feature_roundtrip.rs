//! 特征编辑的语料级黄金测试：对样例库里**每一个** `block 10` 特征做
//! 「解析 → 原样 Update → 写回」，要求输出与输入**逐字节一致**。
//!
//! 这会把整份特征 XML 词汇表（富文本 /note、`noColor`、跨原点片段、多片段、
//! `prioritize`/`visible`/`swappingSegmentNumbering`…）都过一遍：任何解析-重写
//! 不对称都会被抓到。
//!
//! 样例目录优先取 `HELIX_SAMPLES_DIR`；找不到样例时跳过（不失败）。

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

#[test]
fn every_feature_update_is_byte_identical() {
    use hx_formats::feature_edit::{apply_feature_edits, parse_specs, FeatureEdit};
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录（设 HELIX_SAMPLES_DIR 指定）");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    let mut n_files = 0usize;
    let mut n_feats = 0usize;
    for p in files {
        let Ok(bytes) = std::fs::read(&p) else {
            continue;
        };
        let Ok(file) = hx_formats::DnaFile::parse(&bytes) else {
            continue;
        };
        let Some(b10) = file.block(10) else {
            continue;
        };
        let specs = parse_specs(&b10.payload);
        if specs.is_empty() {
            continue;
        }
        // 展示模型与可编辑模型必须一一对应（UI 用下标互相索引）
        let anns = hx_formats::annot::parse_features(&b10.payload);
        assert_eq!(
            anns.len(),
            specs.len(),
            "{}: parse_features 与 parse_specs 数量不一致",
            p.display()
        );
        for (a, s) in anns.iter().zip(&specs) {
            assert_eq!(a.name, s.name, "{}", p.display());
            assert_eq!(a.segments.len(), s.segments.len(), "{}", p.display());
        }
        n_files += 1;
        let edits: Vec<FeatureEdit> = specs.iter().cloned().map(FeatureEdit::Update).collect();
        let out = apply_feature_edits(&b10.payload, &edits);
        assert_eq!(
            String::from_utf8_lossy(&out),
            String::from_utf8_lossy(&b10.payload),
            "{}: 原样 Update 改动了字节",
            p.display()
        );
        n_feats += specs.len();
    }
    assert!(n_files > 0, "样例目录里没有带特征的文件？");
    eprintln!("特征往返：{n_files} 个文件、{n_feats} 个特征全部逐字节一致");
}

/// 「只改名」必须是一处**外科手术**：整份载荷只多 1 个字节，且被编辑特征之外的
/// 字节（前后缀）与特征自身的限定符元素都逐字节不动。
#[test]
fn renaming_a_feature_is_surgical() {
    use hx_formats::feature_edit::{apply_feature_edits, parse_specs, FeatureEdit};
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录（设 HELIX_SAMPLES_DIR 指定）");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    let mut n = 0usize;
    for p in files {
        let Ok(bytes) = std::fs::read(&p) else {
            continue;
        };
        let Ok(file) = hx_formats::DnaFile::parse(&bytes) else {
            continue;
        };
        let Some(b10) = file.block(10) else {
            continue;
        };
        let Ok(text) = std::str::from_utf8(&b10.payload) else {
            continue;
        };
        for sp in parse_specs(&b10.payload) {
            let idpat = format!("recentID=\"{}\"", sp.recent_id);
            let Some(ipos) = text.find(&idpat) else {
                continue;
            };
            let Some(fstart) = text[..ipos].rfind("<Feature") else {
                continue;
            };
            let Some(fend) = text[ipos..]
                .find("</Feature>")
                .map(|i| ipos + i + "</Feature>".len())
            else {
                continue; // 自闭合等罕见形态：跳过
            };
            let mut edited = sp.clone();
            edited.name.push('X');
            let out = apply_feature_edits(&b10.payload, &[FeatureEdit::Update(edited)]);
            let new_text = String::from_utf8_lossy(&out).to_string();
            assert_eq!(
                out.len(),
                b10.payload.len() + 1,
                "{}: 改名一个字节，输出长度却变了 {}",
                p.display(),
                out.len() as isize - b10.payload.len() as isize
            );
            assert!(
                new_text.starts_with(&text[..fstart]),
                "{}: 目标特征之前的内容被改动",
                p.display()
            );
            assert!(
                new_text.ends_with(&text[fend..]),
                "{}: 目标特征之后的内容被改动",
                p.display()
            );
            // 原文该特征里的限定符元素必须原样还在
            let orig_elem = &text[fstart..fend];
            let mut qs = Vec::new();
            let mut at = 0usize;
            while let Some(r) = orig_elem[at..].find("<Q name=\"") {
                let s = at + r;
                let e = orig_elem[s..]
                    .find("</Q>")
                    .map(|i| s + i + "</Q>".len())
                    .unwrap_or(orig_elem.len());
                qs.push(orig_elem[s..e].to_string());
                at = e;
            }
            for q in qs {
                assert!(
                    new_text.contains(&q),
                    "{}: 限定符被重写\n原: {q}",
                    p.display()
                );
            }
            n += 1;
        }
    }
    assert!(n > 0, "样例目录里没有带特征的文件？");
    eprintln!("外科手术式改名：{n} 个特征，全部只多 1 字节");
}

/// 逐特征做「隐藏 → 显示」往返：**必须回到原始字节**。这会把 `visible="0"` 属性的
/// 增删路径（原位更新还是追加）在全语料上过一遍。
#[test]
fn hide_then_show_round_trips_every_feature() {
    use hx_formats::feature_edit::{apply_feature_edits, parse_specs, FeatureEdit};
    let Some(root) = samples_root() else {
        eprintln!("跳过：未找到样例目录（设 HELIX_SAMPLES_DIR 指定）");
        return;
    };
    let mut files = Vec::new();
    collect(&root, &mut files);
    let mut n = 0usize;
    for p in files {
        let Ok(bytes) = std::fs::read(&p) else {
            continue;
        };
        let Ok(file) = hx_formats::DnaFile::parse(&bytes) else {
            continue;
        };
        let Some(b10) = file.block(10) else {
            continue;
        };
        for sp in parse_specs(&b10.payload) {
            let id = sp.recent_id;
            let was = sp.visible;
            let hide = [FeatureEdit::SetVisible {
                recent_id: id,
                visible: !was,
            }];
            let hidden = apply_feature_edits(&b10.payload, &hide);
            let back = apply_feature_edits(
                &hidden,
                &[FeatureEdit::SetVisible {
                    recent_id: id,
                    visible: was,
                }],
            );
            assert_eq!(
                back,
                b10.payload,
                "{}: 特征 recentID={id} 隐藏再显示后与原文不一致",
                p.display()
            );
            n += 1;
        }
    }
    assert!(n > 0, "样例目录里没有带特征的文件？");
    eprintln!("隐藏/显示往返：{n} 个特征全部回到原始字节");
}
