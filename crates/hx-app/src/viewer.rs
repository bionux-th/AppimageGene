//! 文本 / PDF 阅读视图。
//!
//! 这两类文件不是分子文档：它们有自己的页签内容，不进入 [`crate::doc::Document`]、
//! 不参与序列编辑，也不显示 Map/Sequence 等分子页签。
//!
//! * `.txt` / `.text`：**可直接编辑**，Ctrl+S 写回原文件；Ctrl+`+` / Ctrl+`-`
//!   调整字号。左右留出很大的页边距，没有顶部控件条。
//! * `.pdf`：**连续滚动**所有页面（不是一次一页），可见页在后台用系统现有的
//!   `pdfinfo` + `pdftoppm`（poppler-utils）渲染成 PPM，再解码成纹理；滚出视口的
//!   纹理会被回收。Ctrl+`+` / Ctrl+`-` 缩放。找不到工具时给出明确提示，并保留
//!   “用系统程序打开”的退路，不静默失败。

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, TryRecvError};

use eframe::egui;
use eframe::egui::Vec2;

/// 阅读器支持的扩展名（小写）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewerKind {
    Text,
    Pdf,
}

/// 扩展名 → 阅读器类型。大小写不敏感；其它扩展名返回 `None`。
pub fn kind_of_ext(ext: &str) -> Option<ViewerKind> {
    match ext.to_ascii_lowercase().as_str() {
        "txt" | "text" => Some(ViewerKind::Text),
        "pdf" => Some(ViewerKind::Pdf),
        _ => None,
    }
}

/// 路径 → 阅读器类型。
pub fn kind_of_path(path: &Path) -> Option<ViewerKind> {
    path.extension()
        .and_then(|e| e.to_str())
        .and_then(kind_of_ext)
}

/// 文本阅读上限：超过后只显示前 2 MiB，避免一个巨大日志把 UI 线程拖死。
const TEXT_LIMIT: usize = 2 * 1024 * 1024;

const TEXT_FONT_MIN: f32 = 9.0;
const TEXT_FONT_MAX: f32 = 40.0;
const TEXT_FONT_STEP: f32 = 1.5;
const TEXT_FONT_DEFAULT: f32 = 13.0;

/// 原始文本编码：保存时按原编码写回，避免“打开 UTF-16、保存变 UTF-8”。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextEncoding {
    Utf8,
    Utf8Bom,
    Utf16Le,
    Utf16Be,
}

pub struct TextView {
    pub path: PathBuf,
    pub title: String,
    pub text: String,
    pub bytes: usize,
    pub lines: usize,
    pub truncated: bool,
    pub lossy: bool,
    pub font_size: f32,
    /// 有未保存的编辑。
    pub dirty: bool,
    encoding: TextEncoding,
    /// 原文件用 CRLF；保存时还原。
    crlf: bool,
}

impl TextView {
    fn open(path: &Path) -> Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        let (mut text, lossy, encoding) = decode_text(&bytes);
        let truncated = bytes.len() > TEXT_LIMIT;
        if truncated {
            let mut cut = TEXT_LIMIT.min(text.len());
            while cut > 0 && !text.is_char_boundary(cut) {
                cut -= 1;
            }
            text.truncate(cut);
        }
        let crlf = text.contains("\r\n");
        if crlf {
            text = text.replace("\r\n", "\n");
        }
        let lines = line_count(&text);
        Ok(Self {
            path: path.to_path_buf(),
            title: file_title(path),
            text,
            bytes: bytes.len(),
            lines,
            truncated,
            lossy,
            font_size: TEXT_FONT_DEFAULT,
            dirty: false,
            encoding,
            crlf,
        })
    }

    /// 截断 / 非 UTF-8 的预览不能安全写回：编辑会把丢失的部分覆盖掉。
    pub fn editable(&self) -> bool {
        !self.truncated && !self.lossy
    }

    fn summary(&self) -> String {
        let mut s = format!(
            "Text · {} line(s) · {} · {}%",
            nfmt(self.lines),
            nfmt(self.bytes),
            (self.font_size / TEXT_FONT_DEFAULT * 100.0).round() as i32
        );
        if self.dirty {
            s.push_str(" · modified");
        }
        if !self.editable() {
            s.push_str(" · read-only preview");
        }
        s
    }

    fn zoom_in(&mut self) {
        self.font_size = (self.font_size + TEXT_FONT_STEP).min(TEXT_FONT_MAX);
    }

    fn zoom_out(&mut self) {
        self.font_size = (self.font_size - TEXT_FONT_STEP).max(TEXT_FONT_MIN);
    }

    fn reset_zoom(&mut self) {
        self.font_size = TEXT_FONT_DEFAULT;
    }

    fn draw(&mut self, ui: &mut egui::Ui) {
        // 大量左右页边距：像文字处理软件一样，正文不贴边。
        let avail = ui.available_width().max(160.0);
        let margin = (avail * 0.08).clamp(48.0, 160.0);
        let content_w = (avail - margin * 2.0).max(120.0);
        let mut changed = false;
        egui::ScrollArea::vertical()
            .id_salt("text_viewer_scroll")
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.add_space(18.0);
                if self.truncated || self.lossy {
                    ui.horizontal(|ui| {
                        ui.add_space(margin);
                        let note = if self.lossy {
                            "This file is not valid UTF-8; invalid bytes are shown as replacement characters, so editing is disabled."
                        } else {
                            "Preview limited to the first 2 MiB, so editing is disabled."
                        };
                        ui.label(
                            egui::RichText::new(crate::i18n::text(ui.ctx(),note))
                                .size(11.0)
                                .color(crate::theme::TEXT_DIM),
                        );
                    });
                    ui.add_space(6.0);
                }
                ui.horizontal(|ui| {
                    ui.add_space(margin);
                    ui.vertical(|ui| {
                        ui.set_width(content_w);
                        if self.editable() {
                            let resp = ui.add(
                                egui::TextEdit::multiline(&mut self.text)
                                    .font(egui::FontId::monospace(self.font_size))
                                    .text_color(crate::theme::TEXT)
                                    .desired_width(f32::INFINITY)
                                    .desired_rows(30)
                                    .margin(egui::Margin::ZERO)
                                    .frame(
                                        egui::Frame::new()
                                            .fill(egui::Color32::WHITE)
                                            .inner_margin(egui::Margin::symmetric(0, 2)),
                                    ),
                            );
                            changed = resp.changed();
                        } else {
                            ui.add(
                                egui::Label::new(
                                    egui::RichText::new(&self.text)
                                        .monospace()
                                        .size(self.font_size)
                                        .color(crate::theme::TEXT),
                                )
                                .selectable(true)
                                .wrap(),
                            );
                        }
                    });
                });
                ui.add_space(24.0);
            });
        if changed {
            self.dirty = true;
            self.lines = line_count(&self.text);
        }
    }

    /// 写回文件（`target` 为空 = 原路径）。返回写入字节数。
    fn save_to(&mut self, target: Option<&Path>) -> Result<usize, String> {
        if !self.editable() {
            return Err("this file is a read-only preview".to_string());
        }
        let path = target
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| self.path.clone());
        let body = if self.crlf {
            self.text.replace('\n', "\r\n")
        } else {
            self.text.clone()
        };
        let bytes = encode_text(&body, self.encoding);
        std::fs::write(&path, &bytes).map_err(|e| format!("write {}: {e}", path.display()))?;
        self.title = file_title(&path);
        self.path = path;
        self.bytes = bytes.len();
        self.lines = line_count(&self.text);
        self.dirty = false;
        Ok(bytes.len())
    }
}

fn line_count(text: &str) -> usize {
    if text.is_empty() {
        0
    } else {
        text.lines().count()
    }
}

// ---------------------------------------------------------------------------
// PDF
// ---------------------------------------------------------------------------

/// 连续滚动时同时保留的页面纹理上限（200 dpi 的 A4 约 15 MB/页）。
const PDF_TEXTURE_CACHE: usize = 8;
/// 同时渲染的页数上限：每页是一个 `pdftoppm` 进程，太多会抢内存。
const PDF_RENDER_PARALLEL: usize = 2;
/// 页面之间的间隙。
const PDF_PAGE_GAP: f32 = 14.0;
/// 页面两侧的最小留白。
const PDF_SIDE_MARGIN: f32 = 24.0;

struct PageTexture {
    handle: egui::TextureHandle,
    used: u64,
}

pub struct PdfView {
    pub path: PathBuf,
    pub title: String,
    pub pages: usize,
    /// 显示缩放；1.0 = 96 dpi 下的物理大小（纹理按 200 dpi 渲染）。
    pub zoom: f32,
    pub fit_width: bool,
    /// 滚动时最靠上的可见页（1-based，和状态栏一致）。
    pub page: usize,
    pub error: Option<String>,
    page_pt: Vec2,
    textures: HashMap<usize, PageTexture>,
    pending: Vec<RenderJob>,
    page_errors: HashMap<usize, String>,
    last_width: f32,
    frame: u64,
}

struct RenderJob {
    page: usize,
    rx: Receiver<Result<Ppm, String>>,
}

impl PdfView {
    fn open(path: &Path) -> Self {
        let mut view = Self {
            path: path.to_path_buf(),
            title: file_title(path),
            pages: 0,
            zoom: 1.0,
            fit_width: true,
            page: 1,
            error: None,
            page_pt: Vec2::new(595.0, 842.0),
            textures: HashMap::new(),
            pending: Vec::new(),
            page_errors: HashMap::new(),
            last_width: 900.0,
            frame: 0,
        };
        if crate::sysopen::which("pdfinfo").is_none() || crate::sysopen::which("pdftoppm").is_none() {
            view.error = Some(
                "PDF preview needs poppler-utils (`pdfinfo` and `pdftoppm`). Install it, or use Open externally."
                    .to_string(),
            );
            return view;
        }
        match pdf_info(path) {
            Ok(info) if info.pages > 0 => {
                view.pages = info.pages;
                view.page_pt = info.page_pt;
            }
            Ok(_) => view.error = Some("This PDF has no pages.".to_string()),
            Err(e) => view.error = Some(e),
        }
        view
    }

    fn summary(&self) -> String {
        if let Some(e) = &self.error {
            return format!("PDF · {e}");
        }
        format!(
            "PDF · page {} of {} · {}%",
            self.page,
            self.pages,
            (self.effective_zoom() * 100.0).round() as i32
        )
    }

    /// 页面在 96 dpi 下的自然尺寸（纹理是 200 dpi）。
    fn natural(&self) -> Vec2 {
        self.page_pt * (96.0 / 72.0)
    }

    fn zoom_in(&mut self) {
        let current = self.effective_zoom();
        self.fit_width = false;
        self.zoom = (current * 1.2).min(4.0);
    }

    fn zoom_out(&mut self) {
        let current = self.effective_zoom();
        self.fit_width = false;
        self.zoom = (current / 1.2).max(0.25);
    }

    fn fit(&mut self) {
        self.fit_width = true;
    }

    /// 当前实际缩放（fit-width 时由窗口宽度反推）。
    fn effective_zoom(&self) -> f32 {
        if self.fit_width {
            let n = self.natural();
            let page_w = (self.last_width - PDF_SIDE_MARGIN * 2.0).max(120.0);
            if n.x > 0.0 {
                (page_w / n.x).clamp(0.25, 4.0)
            } else {
                1.0
            }
        } else {
            self.zoom
        }
    }

    fn start_render(&mut self, page: usize) {
        let path = self.path.clone();
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(render_page(&path, page));
        });
        self.pending.push(RenderJob { page, rx });
    }

    fn poll_render(&mut self, ctx: &egui::Context) {
        if self.pending.is_empty() {
            return;
        }
        let mut done: Vec<(usize, usize, Result<Ppm, String>)> = Vec::new();
        for (i, job) in self.pending.iter().enumerate() {
            match job.rx.try_recv() {
                Ok(res) => done.push((i, job.page, res)),
                Err(TryRecvError::Empty) => {
                    ctx.request_repaint_after(std::time::Duration::from_millis(40));
                }
                Err(TryRecvError::Disconnected) => done.push((
                    i,
                    job.page,
                    Err("PDF renderer stopped unexpectedly.".to_string()),
                )),
            }
        }
        // 倒序删除，索引不会失效。
        for (i, page, res) in done.into_iter().rev() {
            self.pending.remove(i);
            match res {
                Ok(ppm) => {
                    let image =
                        egui::ColorImage::from_rgba_unmultiplied([ppm.width, ppm.height], &ppm.rgba);
                    let handle = ctx.load_texture(
                        format!("pdf-page-{}-{page}", self.path.display()),
                        image,
                        egui::TextureOptions::LINEAR,
                    );
                    self.page_errors.remove(&page);
                    self.textures.insert(
                        page,
                        PageTexture {
                            handle,
                            used: self.frame,
                        },
                    );
                }
                Err(e) => {
                    self.page_errors.insert(page, e);
                }
            }
        }
        if !self.pending.is_empty() {
            ctx.request_repaint_after(std::time::Duration::from_millis(40));
        }
    }

    /// 回收滚出视口的纹理，避免长文档吃光显存。
    fn evict(&mut self, visible: &[usize]) {
        if self.textures.len() <= PDF_TEXTURE_CACHE {
            return;
        }
        let mut stale: Vec<(u64, usize)> = self
            .textures
            .iter()
            .filter(|(p, _)| !visible.contains(p))
            .map(|(p, t)| (t.used, *p))
            .collect();
        stale.sort_unstable();
        for (_, page) in stale {
            if self.textures.len() <= PDF_TEXTURE_CACHE {
                break;
            }
            self.textures.remove(&page);
        }
    }

    fn draw(&mut self, ui: &mut egui::Ui) {
        self.frame = self.frame.wrapping_add(1);
        self.poll_render(ui.ctx());

        if let Some(e) = self.error.clone() {
            ui.vertical_centered(|ui| {
                ui.add_space(40.0);
                ui.label(egui::RichText::new(crate::i18n::message(ui.ctx(),&e)).size(13.0).color(crate::theme::SITE));
                ui.add_space(8.0);
                if ui.button(crate::i18n::text(ui.ctx(),"Open with system PDF viewer")).clicked() {
                    let _ = crate::sysopen::open_path(&self.path);
                }
                ui.add_space(6.0);
                ui.label(
                    egui::RichText::new(self.path.display().to_string())
                        .size(11.0)
                        .color(crate::theme::TEXT_DIM),
                );
            });
            return;
        }

        let avail = ui.available_width().max(160.0);
        if (self.last_width - avail).abs() > 0.5 {
            // 标题/状态栏在本帧绘制之前取的是上一帧的宽度；窗口尺寸变化后
            // 再排一帧，让 fit-width 的百分比立即跟上。
            self.last_width = avail;
            ui.ctx().request_repaint();
        } else {
            self.last_width = avail;
        }
        let natural = self.natural();
        let page_w = if self.fit_width {
            (avail - PDF_SIDE_MARGIN * 2.0).max(120.0)
        } else {
            (natural.x * self.zoom).max(80.0)
        };
        let page_h = page_w * natural.y / natural.x;
        let row_w = page_w.max(avail - 1.0);

        let mut visible: Vec<usize> = Vec::new();
        let mut top_page = 0usize;
        egui::ScrollArea::both()
            .id_salt("pdf_viewer_scroll")
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.add_space(10.0);
                for page in 1..=self.pages {
                    let (row, _) =
                        ui.allocate_exact_size(Vec2::new(row_w, page_h), egui::Sense::hover());
                    if ui.is_rect_visible(row) {
                        if top_page == 0 {
                            top_page = page;
                        }
                        visible.push(page);
                        let img = egui::Rect::from_center_size(row.center(), Vec2::new(page_w, page_h));
                        let p = ui.painter();
                        // 页面阴影 + 白纸。
                        p.rect_filled(img.expand(3.0), 2.0, egui::Color32::from_gray(196));
                        p.rect_filled(img, 0.0, egui::Color32::WHITE);
                        if let Some(tex) = self.textures.get_mut(&page) {
                            tex.used = self.frame;
                            p.image(
                                tex.handle.id(),
                                img,
                                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                                egui::Color32::WHITE,
                            );
                        } else if let Some(e) = self.page_errors.get(&page) {
                            p.text(
                                img.center(),
                                egui::Align2::CENTER_CENTER,
                                crate::ui_format!(ui.ctx(),"Page {}: {}",page,e),
                                egui::FontId::proportional(12.0),
                                crate::theme::SITE,
                            );
                        } else {
                            p.text(
                                img.center(),
                                egui::Align2::CENTER_CENTER,
                                crate::ui_format!(ui.ctx(),"Rendering page {}…",page),
                                egui::FontId::proportional(12.0),
                                crate::theme::TEXT_DIM,
                            );
                        }
                    }
                    ui.add_space(PDF_PAGE_GAP);
                }
                ui.add_space(10.0);
            });

        if top_page > 0 {
            self.page = top_page;
        }

        // 只为可见页启动渲染，并限制并发。
        for &page in &visible {
            if self.textures.contains_key(&page) || self.page_errors.contains_key(&page) {
                continue;
            }
            if self.pending.iter().any(|j| j.page == page) {
                continue;
            }
            if self.pending.len() >= PDF_RENDER_PARALLEL {
                break;
            }
            self.start_render(page);
        }
        self.evict(&visible);
        if !self.pending.is_empty() {
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_millis(40));
        }
    }
}

/// 一个已打开的文件。
pub enum ViewerDoc {
    Text(TextView),
    Pdf(PdfView),
}

impl ViewerDoc {
    /// 只接管 `.txt` / `.text` / `.pdf`；其它扩展名返回 `Ok(None)`，交给原文件流程。
    pub fn open(path: &Path) -> Result<Option<Self>, String> {
        match kind_of_path(path) {
            Some(ViewerKind::Text) => TextView::open(path).map(|v| Some(Self::Text(v))),
            Some(ViewerKind::Pdf) => Ok(Some(Self::Pdf(PdfView::open(path)))),
            None => Ok(None),
        }
    }

    pub fn path(&self) -> &Path {
        match self {
            Self::Text(v) => &v.path,
            Self::Pdf(v) => &v.path,
        }
    }

    pub fn title(&self) -> String {
        match self {
            Self::Text(v) => v.title.clone(),
            Self::Pdf(v) => v.title.clone(),
        }
    }

    pub fn summary(&self) -> String {
        match self {
            Self::Text(v) => v.summary(),
            Self::Pdf(v) => v.summary(),
        }
    }

    pub fn is_pdf(&self) -> bool {
        matches!(self, Self::Pdf(_))
    }

    pub fn is_dirty(&self) -> bool {
        matches!(self, Self::Text(v) if v.dirty)
    }

    /// 这个页签可以保存吗（只有可编辑文本可以；PDF 只读）。
    pub fn can_save(&self) -> bool {
        matches!(self, Self::Text(v) if v.editable())
    }

    pub fn zoom_in(&mut self) {
        match self {
            Self::Text(v) => v.zoom_in(),
            Self::Pdf(v) => v.zoom_in(),
        }
    }

    pub fn zoom_out(&mut self) {
        match self {
            Self::Text(v) => v.zoom_out(),
            Self::Pdf(v) => v.zoom_out(),
        }
    }

    pub fn fit(&mut self) {
        match self {
            Self::Text(v) => v.reset_zoom(),
            Self::Pdf(v) => v.fit(),
        }
    }

    /// 保存（`target` 为空 = 原路径）。返回写入字节数。
    pub fn save_to(&mut self, target: Option<&Path>) -> Result<usize, String> {
        match self {
            Self::Text(v) => v.save_to(target),
            Self::Pdf(_) => Err("PDF files are read-only".to_string()),
        }
    }

    pub fn draw(&mut self, ui: &mut egui::Ui) {
        match self {
            Self::Text(v) => v.draw(ui),
            Self::Pdf(v) => v.draw(ui),
        }
    }
}

/// 文件名的显示形式（不带目录）。
fn file_title(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| path.display().to_string())
}

/// 文本解码：UTF-8（可带 BOM）/ UTF-16 LE/BE BOM；其余按 UTF-8 lossy。
fn decode_text(bytes: &[u8]) -> (String, bool, TextEncoding) {
    if let Some(rest) = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]) {
        return (
            String::from_utf8_lossy(rest).into_owned(),
            false,
            TextEncoding::Utf8Bom,
        );
    }
    if let Some(rest) = bytes.strip_prefix(&[0xff, 0xfe]) {
        return (decode_utf16(rest, true), false, TextEncoding::Utf16Le);
    }
    if let Some(rest) = bytes.strip_prefix(&[0xfe, 0xff]) {
        return (decode_utf16(rest, false), false, TextEncoding::Utf16Be);
    }
    match std::str::from_utf8(bytes) {
        Ok(s) => (s.to_string(), false, TextEncoding::Utf8),
        Err(_) => (
            String::from_utf8_lossy(bytes).into_owned(),
            true,
            TextEncoding::Utf8,
        ),
    }
}

fn decode_utf16(bytes: &[u8], little_endian: bool) -> String {
    let units: Vec<u16> = bytes
        .chunks_exact(2)
        .map(|c| {
            if little_endian {
                u16::from_le_bytes([c[0], c[1]])
            } else {
                u16::from_be_bytes([c[0], c[1]])
            }
        })
        .collect();
    String::from_utf16_lossy(&units)
}

fn encode_text(text: &str, encoding: TextEncoding) -> Vec<u8> {
    match encoding {
        TextEncoding::Utf8 => text.as_bytes().to_vec(),
        TextEncoding::Utf8Bom => {
            let mut out = vec![0xef, 0xbb, 0xbf];
            out.extend_from_slice(text.as_bytes());
            out
        }
        TextEncoding::Utf16Le | TextEncoding::Utf16Be => {
            let little = encoding == TextEncoding::Utf16Le;
            let mut out = if little {
                vec![0xff, 0xfe]
            } else {
                vec![0xfe, 0xff]
            };
            for unit in text.encode_utf16() {
                out.extend_from_slice(&if little {
                    unit.to_le_bytes()
                } else {
                    unit.to_be_bytes()
                });
            }
            out
        }
    }
}

// ---------------------------------------------------------------------------
// poppler
// ---------------------------------------------------------------------------

struct PdfInfo {
    pages: usize,
    page_pt: Vec2,
}

/// `pdfinfo` 输出里的 `Pages:` 行。
fn parse_pdf_pages(stdout: &str) -> Option<usize> {
    stdout.lines().find_map(|line| {
        line.strip_prefix("Pages:")
            .and_then(|n| n.trim().parse::<usize>().ok())
    })
}

/// `pdfinfo` 输出里的 `Page size: 595.276 x 841.89 pts (A4)` 行（单位 pt）。
fn parse_pdf_page_size(stdout: &str) -> Option<Vec2> {
    stdout.lines().find_map(|line| {
        let rest = line.strip_prefix("Page size:")?.trim();
        let (w, h) = rest.split_once('x')?;
        let w = w.trim().parse::<f32>().ok()?;
        let h = h.trim().split_whitespace().next()?.parse::<f32>().ok()?;
        (w > 0.0 && h > 0.0).then(|| Vec2::new(w, h))
    })
}

fn pdf_info(path: &Path) -> Result<PdfInfo, String> {
    let out = crate::sysopen::helper_command(crate::sysopen::which("pdfinfo").ok_or("pdfinfo not found")?)
        .arg(path)
        .env("LANG", "C")
        .env("LC_ALL", "C")
        .output()
        .map_err(|e| format!("run pdfinfo: {e}"))?;
    if !out.status.success() {
        let msg = String::from_utf8_lossy(&out.stderr).trim().to_string();
        return Err(if msg.is_empty() {
            "pdfinfo could not read this PDF.".to_string()
        } else {
            msg
        });
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let pages = parse_pdf_pages(&text)
        .ok_or_else(|| "pdfinfo did not report a page count.".to_string())?;
    Ok(PdfInfo {
        pages,
        page_pt: parse_pdf_page_size(&text).unwrap_or(Vec2::new(595.0, 842.0)),
    })
}

/// 200 dpi 的 PPM（P6）页面。
struct Ppm {
    width: usize,
    height: usize,
    rgba: Vec<u8>,
}

fn render_page(path: &Path, page: usize) -> Result<Ppm, String> {
    // 这个 poppler 的 `pdftoppm` 把 `-` 当成文件名（生成 `-.ppm`），不写 stdout，
    // 所以渲染到唯一临时文件，再读回来删除。
    static NONCE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = NONCE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let stem = std::env::temp_dir().join(format!("hx-pdf-{}-{page}-{n}", std::process::id()));
    let rendered = stem.with_extension("ppm");
    let out = crate::sysopen::helper_command(crate::sysopen::which("pdftoppm").ok_or("pdftoppm not found")?)
        .args([
            "-r",
            "200",
            "-f",
            &page.to_string(),
            "-l",
            &page.to_string(),
            "-singlefile",
        ])
        .arg(path)
        .arg(&stem)
        .env("LANG", "C")
        .env("LC_ALL", "C")
        .output()
        .map_err(|e| format!("run pdftoppm: {e}"))?;
    if !out.status.success() {
        let msg = String::from_utf8_lossy(&out.stderr).trim().to_string();
        let _ = std::fs::remove_file(&rendered);
        return Err(if msg.is_empty() {
            format!("pdftoppm could not render page {page}.")
        } else {
            msg
        });
    }
    let bytes =
        std::fs::read(&rendered).map_err(|e| format!("read rendered page {page}: {e}"))?;
    let _ = std::fs::remove_file(&rendered);
    parse_ppm(&bytes)
}

/// 解析 P6 PPM（支持 `#` 注释）。返回 RGBA8，供 egui 纹理直接使用。
fn parse_ppm(data: &[u8]) -> Result<Ppm, String> {
    let mut at = 0usize;
    let mut token = |data: &[u8]| -> Option<String> {
        loop {
            while at < data.len() && data[at].is_ascii_whitespace() {
                at += 1;
            }
            if at < data.len() && data[at] == b'#' {
                while at < data.len() && data[at] != b'\n' {
                    at += 1;
                }
                continue;
            }
            break;
        }
        let start = at;
        while at < data.len() && !data[at].is_ascii_whitespace() && data[at] != b'#' {
            at += 1;
        }
        (at > start).then(|| String::from_utf8_lossy(&data[start..at]).to_string())
    };
    if token(data).as_deref() != Some("P6") {
        return Err("PDF renderer did not return a P6 PPM image.".to_string());
    }
    let width = token(data)
        .and_then(|t| t.parse::<usize>().ok())
        .ok_or("malformed PPM: width")?;
    let height = token(data)
        .and_then(|t| t.parse::<usize>().ok())
        .ok_or("malformed PPM: height")?;
    let max = token(data)
        .and_then(|t| t.parse::<usize>().ok())
        .ok_or("malformed PPM: max value")?;
    if width == 0 || height == 0 || max == 0 || max > 255 {
        return Err("unsupported PPM dimensions or sample depth.".to_string());
    }
    // 头部之后恰好一个空白字节，随后是像素。
    if at < data.len() && data[at].is_ascii_whitespace() {
        at += 1;
    }
    let need = width
        .checked_mul(height)
        .and_then(|p| p.checked_mul(3))
        .ok_or("PPM too large")?;
    let pixels = data
        .get(at..at + need)
        .ok_or_else(|| format!("truncated PPM: expected {need} pixel bytes"))?;
    let mut rgba = Vec::with_capacity(width * height * 4);
    for px in pixels.chunks_exact(3) {
        rgba.extend_from_slice(&[px[0], px[1], px[2], 255]);
    }
    Ok(Ppm {
        width,
        height,
        rgba,
    })
}

fn nfmt(n: usize) -> String {
    crate::app::nfmt(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extensions_are_case_insensitive_and_limited_to_reading_files() {
        assert_eq!(kind_of_ext("TXT"), Some(ViewerKind::Text));
        assert_eq!(kind_of_ext("pdf"), Some(ViewerKind::Pdf));
        assert_eq!(kind_of_ext("dna"), None);
        assert_eq!(kind_of_path(Path::new("/tmp/a.PDF")), Some(ViewerKind::Pdf));
    }

    #[test]
    fn text_decoding_handles_utf8_utf16_and_invalid_bytes() {
        assert_eq!(decode_text("héllo".as_bytes()).0, "héllo");
        let mut le = vec![0xff, 0xfe];
        for u in "hello".encode_utf16() {
            le.extend_from_slice(&u.to_le_bytes());
        }
        assert_eq!(decode_text(&le).0, "hello");
        assert_eq!(decode_text(&le).2, TextEncoding::Utf16Le);
        let mut be = vec![0xfe, 0xff];
        for u in "hello".encode_utf16() {
            be.extend_from_slice(&u.to_be_bytes());
        }
        assert_eq!(decode_text(&be).0, "hello");
        assert_eq!(decode_text(&be).2, TextEncoding::Utf16Be);
        let (s, lossy, _) = decode_text(b"ok\xffbad");
        assert!(lossy);
        assert!(s.starts_with("ok"));
    }

    #[test]
    fn text_round_trips_through_each_encoding() {
        for enc in [
            TextEncoding::Utf8,
            TextEncoding::Utf8Bom,
            TextEncoding::Utf16Le,
            TextEncoding::Utf16Be,
        ] {
            let bytes = encode_text("aé\n", enc);
            let (text, lossy, back) = decode_text(&bytes);
            assert_eq!(text, "aé\n", "{enc:?}");
            assert_eq!(back, enc);
            assert!(!lossy);
        }
    }

    #[test]
    fn ppm_parser_accepts_comments_and_rejects_truncated_pixels() {
        let ppm = b"P6\n# generated\n2 1\n255\n\xff\x00\x00\x00\xff\x00";
        let p = parse_ppm(ppm).unwrap();
        assert_eq!((p.width, p.height), (2, 1));
        assert_eq!(p.rgba, vec![255, 0, 0, 255, 0, 255, 0, 255]);
        assert!(parse_ppm(b"P6\n2 1\n255\n\xff").is_err());
        assert!(parse_ppm(b"P3\n1 1\n255\n0 0 0").is_err());
    }

    #[test]
    fn pdf_page_count_reads_the_english_pdfinfo_line() {
        assert_eq!(
            parse_pdf_pages("Title: x\nPages:          12\nEncrypted: no\n"),
            Some(12)
        );
        assert_eq!(parse_pdf_pages("Title: x\n"), None);
    }

    #[test]
    fn pdf_page_size_is_parsed_in_points() {
        let size = parse_pdf_page_size("Page size:       595.276 x 841.89 pts (A4)\n").unwrap();
        assert!((size.x - 595.276).abs() < 0.01);
        assert!((size.y - 841.89).abs() < 0.01);
        assert!(parse_pdf_page_size("Pages: 3\n").is_none());
    }

    #[test]
    fn pdf_zoom_is_relative_to_fit_width_and_clamped() {
        let mut view = PdfView {
            path: PathBuf::from("/tmp/paper.pdf"),
            title: "paper.pdf".into(),
            pages: 3,
            zoom: 1.0,
            fit_width: true,
            page: 1,
            error: None,
            page_pt: Vec2::new(595.0, 842.0),
            textures: HashMap::new(),
            pending: Vec::new(),
            page_errors: HashMap::new(),
            last_width: 900.0,
            frame: 0,
        };
        let fit = view.effective_zoom();
        view.zoom_in();
        assert!(!view.fit_width);
        assert!(view.zoom > fit);
        for _ in 0..100 {
            view.zoom_out();
        }
        assert!(view.zoom >= 0.25);
        view.fit();
        assert!(view.fit_width);
    }

    #[test]
    fn text_editor_saves_back_in_place_and_clears_dirty() {
        let path = std::env::temp_dir().join(format!(
            "hx-viewer-save-{}-{}.txt",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"one\ntwo\n").unwrap();
        let mut view = TextView::open(&path).unwrap();
        assert!(view.editable());
        assert!(!view.dirty);
        view.text.push_str("three\n");
        view.dirty = true;
        assert_eq!(view.save_to(None).unwrap(), 14);
        assert!(!view.dirty);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "one\ntwo\nthree\n");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn crlf_files_are_normalized_for_editing_and_restored_on_save() {
        let path = std::env::temp_dir().join(format!(
            "hx-viewer-crlf-{}-{}.txt",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"a\r\nb\r\n").unwrap();
        let mut view = TextView::open(&path).unwrap();
        assert_eq!(view.text, "a\nb\n");
        assert!(view.crlf);
        view.text.push_str("c\n");
        view.save_to(None).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"a\r\nb\r\nc\r\n");
        let _ = std::fs::remove_file(&path);
    }

    fn ui_with<R>(mut f: impl FnMut(&mut egui::Ui) -> R) -> R {
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx);
        let mut out = None;
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::Vec2::new(1000.0, 700.0),
                )),
                ..Default::default()
            },
            |ui| out = Some(f(ui)),
        );
        output.textures_delta.clear();
        out.unwrap()
    }

    #[test]
    fn text_viewer_draws_its_content_without_a_molecule_document() {
        let mut view = TextView {
            path: PathBuf::from("/tmp/notes.txt"),
            title: "notes.txt".into(),
            text: "first line\nsecond line\n".into(),
            bytes: 23,
            lines: 2,
            truncated: false,
            lossy: false,
            font_size: 13.0,
            dirty: false,
            encoding: TextEncoding::Utf8,
            crlf: false,
        };
        ui_with(|ui| view.draw(ui));
        assert_eq!(view.lines, 2);
        assert!(!view.dirty);
    }

    #[test]
    fn missing_pdf_renderer_is_reported_without_panicking() {
        let mut view = PdfView {
            path: PathBuf::from("/tmp/paper.pdf"),
            title: "paper.pdf".into(),
            pages: 0,
            zoom: 1.0,
            fit_width: true,
            page: 1,
            error: Some("PDF preview needs poppler-utils".into()),
            page_pt: Vec2::new(595.0, 842.0),
            textures: HashMap::new(),
            pending: Vec::new(),
            page_errors: HashMap::new(),
            last_width: 900.0,
            frame: 0,
        };
        ui_with(|ui| view.draw(ui));
    }
}
