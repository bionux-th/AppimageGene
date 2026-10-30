//! 把 `vendor/viennarna` 里的 C 源码静态编进当前 crate。
//!
//! 目标：**不依赖任何运行时组件**（不需要 Python、不需要系统装的 ViennaRNA），
//! Linux / macOS / Windows 都一样，只要构建时有 C 编译器。
//!
//! 这里刻意不使用 ViennaRNA 自带的 `config.h`（它是 autoconf 在 Linux 上生成的，
//! 含 `HAVE_*` 平台假设）。不定义 `HAVE_CONFIG_H` 时，源码走的是自带的可移植
//! 回退分支；只额外定义版本号与 naview 布局开关。

use std::path::{Path, PathBuf};

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = manifest.join("../../vendor/viennarna");
    let src = vendor.join("src");

    let mut files = Vec::new();
    collect_c(&src, &mut files);
    files.sort();
    assert!(!files.is_empty(), "vendor/viennarna 里没有找到 C 源码");

    let mut build = cc::Build::new();
    build
        // 先放自己的可移植 config.h，避免上游 config.h 的平台假设。
        .include(manifest.join("shim"))
        .include(&src)
        .include(src.join("ViennaRNA"))
        .include(src.join("cephes"))
        .define("VRNA_VERSION", "\"2.7.2\"")
        .define("VRNA_WITH_NAVIEW_LAYOUT", "1")
        // MSVC 下 `<cmath>` 的 M_PI 需要这个宏；其它编译器无副作用。
        .define("_USE_MATH_DEFINES", None)
        .warnings(false)
        .flag_if_supported("-fno-strict-aliasing")
        .files(&files)
        .file(manifest.join("shim/vienna_shim.c"));

    // 只为了在 build.rs 里少写一次路径。
    let _ = Path::new(&vendor);

    build.compile("viennarna");

    println!("cargo:rerun-if-changed=shim/vienna_shim.c");
    println!("cargo:rerun-if-changed={}", vendor.display());
}

/// 递归收集 `.c`（ViennaRNA 的源码目录结构是固定的两层）。
fn collect_c(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_c(&path, out);
        } else if path.extension().is_some_and(|e| e == "c") {
            out.push(path);
        }
    }
}
