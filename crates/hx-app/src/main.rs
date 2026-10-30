//! HelixBench 桌面入口：`helixbench [file...]`（多个文件会各自打开一个标签页）。

#![forbid(unsafe_code)]
#![cfg_attr(target_os = "windows", windows_subsystem = "windows")]

fn main() -> eframe::Result<()> {
    let paths: Vec<std::path::PathBuf> = std::env::args_os()
        .skip(1)
        .map(std::path::PathBuf::from)
        .collect();
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1320.0, 840.0])
            .with_maximized(true)
            .with_title("AppimageGene")
            .with_app_id("appimagegene")
            .with_icon(hx_app::branding::app_icon()),
        ..Default::default()
    };
    eframe::run_native(
        "AppimageGene",
        options,
        Box::new(move |cc| Ok(Box::new(hx_app::HelixApp::new(cc, paths)))),
    )
}
