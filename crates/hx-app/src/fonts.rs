//! UI 字体：与真 SnapGene 对齐。
//!
//! SnapGene 是 Qt 程序，界面字体取自 fontconfig 的默认 sans —— 本机
//! `fc-match sans-serif` = **Noto Sans / Noto Sans CJK**（cap 高 ≈0.71 em，Regular 粗细）。
//! egui 自带的是 **Ubuntu-Light**：字形更矮更细，于是同样的字号在两边看起来
//! 「字号、粗细都偏小」（oracle 实测：同一行菜单，SnapGene 的字高 10 px、
//! 我们是 9 px，而字宽几乎相同 —— 差距主要来自字重与 x-height）。
//!
//! 这里把系统 sans 装成 proportional 首选、CJK 字体装成 fallback（顺带解决中文
//! 文件名显示方块），**monospace 保持 egui 自带的 Hack 不变**：DNA 序列视图的
//! 字距是前几轮按 oracle 对齐过的，换等宽字体会把它破坏掉。

use eframe::egui::{self, FontData, FontDefinitions, FontFamily};
use std::sync::Arc;

/// proportional 首选字体的键名。
const SANS_KEY: &str = "hb-sans";
/// CJK fallback 的键名。
const CJK_KEY: &str = "hb-cjk";
const BOLD_KEY: &str = "hb-sans-bold";

/// Use a real bold face from the installed regular face's family. Never
/// simulate weight by painting the same text twice with an offset.
pub fn document_title_font(ctx:&egui::Context,active:bool)->egui::FontId {
    let bold=FontFamily::Name(BOLD_KEY.into());
    let family=if active && ctx.data(|data|data.get_temp::<bool>(egui::Id::unique(BOLD_KEY)).unwrap_or(false)) {
        bold
    } else {FontFamily::Proportional};
    egui::FontId::new(13.0,family)
}

/// 系统 sans 候选路径（按顺序取第一个存在的）。
const SANS: &[&str] = &[
    "/usr/share/fonts/noto/NotoSans-Regular.ttf",
    "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
    "/usr/share/fonts/TTF/NotoSans-Regular.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    "/usr/share/fonts/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
];

/// CJK fallback 候选路径（中文文件名/中文文本不再显示成方块）。
const CJK: &[&str] = &[
    "/usr/share/fonts/noto-cjk/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/noto-cjk/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/wqy-microhei/wqy-microhei.ttc",
];

fn read_first(paths: &'static [&'static str]) -> Option<(&'static str, Vec<u8>)> {
    paths
        .iter()
        .find_map(|p| std::fs::read(p).ok().map(|b| (*p, b)))
}

/// 装载系统字体到 `ctx`；返回实际使用的 sans 路径（找不到则 `None`，界面仍用 egui 自带字体）。
pub fn install(ctx: &egui::Context) -> Option<&'static str> {
    #[cfg(windows)]
    let (sans_path,sans_bytes)={let root=std::env::var_os("WINDIR").map(std::path::PathBuf::from).unwrap_or_else(||"C:\\Windows".into());let bytes=std::fs::read(root.join("Fonts/segoeui.ttf")).ok()?;("Windows Segoe UI",bytes)};
    #[cfg(not(windows))]
    let (sans_path, sans_bytes) = read_first(SANS)?;
    #[cfg(windows)]
    let cjk={let root=std::env::var_os("WINDIR").map(std::path::PathBuf::from).unwrap_or_else(||"C:\\Windows".into());["msyh.ttc","msyh.ttf","simhei.ttf","simsun.ttc"].iter().find_map(|name|std::fs::read(root.join("Fonts").join(name)).ok().map(|b|("Windows CJK",b)))};
    #[cfg(not(windows))]
    let cjk = read_first(CJK);

    let mut fonts = FontDefinitions::default();
    fonts.font_data.insert(
        SANS_KEY.to_owned(),
        Arc::new(FontData::from_owned(sans_bytes)),
    );
    // Pair the chosen Regular face with its sibling Bold, not an unrelated
    // fallback family, so active/inactive document titles keep the same glyphs.
    #[cfg(windows)]
    let bold_path=std::env::var_os("WINDIR").map(std::path::PathBuf::from).unwrap_or_else(||"C:\\Windows".into()).join("Fonts/segoeuib.ttf").to_string_lossy().into_owned();
    #[cfg(not(windows))]
    let bold_path=sans_path.replace("-Regular.ttf","-Bold.ttf").replace("DejaVuSans.ttf","DejaVuSans-Bold.ttf");
    let has_bold=bold_path!=sans_path && std::fs::read(&bold_path).ok().map(|bytes|{
        fonts.font_data.insert(BOLD_KEY.to_owned(),Arc::new(FontData::from_owned(bytes)));
    }).is_some();
    if let Some((_, bytes)) = &cjk {
        fonts.font_data.insert(
            CJK_KEY.to_owned(),
            Arc::new(FontData::from_owned(bytes.clone())),
        );
    }

    // proportional：系统 sans → CJK → egui 自带（Ubuntu-Light / emoji）兜底
    if let Some(list) = fonts.families.get_mut(&FontFamily::Proportional) {
        let mut want = vec![SANS_KEY.to_owned()];
        if cjk.is_some() {
            want.push(CJK_KEY.to_owned());
        }
        want.append(list);
        *list = want;
    }
    // monospace：Hack 仍旧第一（DNA 序列字距不变），只在最后挂 CJK 兜底
    if cjk.is_some() {
        if let Some(list) = fonts.families.get_mut(&FontFamily::Monospace) {
            list.push(CJK_KEY.to_owned());
        }
    }

    if has_bold {
        let mut family=vec![BOLD_KEY.to_owned()];
        family.extend(fonts.families[&FontFamily::Proportional].iter().cloned());
        fonts.families.insert(FontFamily::Name(BOLD_KEY.into()),family);
    }

    ctx.data_mut(|data|data.insert_temp(egui::Id::unique(BOLD_KEY),has_bold));
    ctx.set_fonts(fonts);
    Some(sans_path)
}
