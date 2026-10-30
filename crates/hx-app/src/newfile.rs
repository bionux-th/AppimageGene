//! 「New DNA File / New RNA File / New Protein File」对话框（SnapGene 8.x 的
//! *New File* 面板，见 `docs/new-file-and-features-2026-10-18.zh-CN.md`）。
//!
//! 三种模式（oracle：用户提供的 SnapGene 8.2.2 截图，中文界面「新建 DNA 文件」/
//! 「新蛋白质文件」）：
//!
//! ```text
//!  ┌ New DNA File ▾ ─────────────────────────────────┐
//!  │ ◉ Create a sequence                             │
//!  │   ┌───────────────────────────────────────────┐ │
//!  │   │ GAGCG...                                  │ │
//!  │   └───────────────────────────────────────────┘ │
//!  │   [ ] Reverse complement sequence  [x] Detect common features │
//!  │ ○ Create a sequence with multiple N's           │
//!  │ ○ Create a placeholder file with no sequence    │
//!  │ Strandedness   (•) Double-stranded ( ) Single-stranded │
//!  │ Topology       (•) Circular        ( ) Linear   │
//!  │ File Name      [Untitled.dna]                   │
//!  │                            [Cancel] [Create]    │
//!  └─────────────────────────────────────────────────┘
//! ```
//!
//! 打开对话框时立刻建立一个**预览标签页**（SnapGene 的 New File 面板会在主窗口里
//! 实时预览；页签名 `New Sequence`，图谱标签用文件名的词干），Cancel 关掉它并丢弃。
//! 预览文档是**未保存**文档：首次 Save/Save As 才落盘（与 SnapGene 一致）。

use eframe::egui::{self, RichText};
use hx_core::{MoleculeFlags, MoleculeKind, Strandedness, Topology};
use hx_formats::write::NewFile;

use crate::doc::Document;
use crate::theme;

/// 对话框里的三种「创建」模式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewMode {
    /// 创建序列（粘贴/键入）。
    Sequence,
    /// 创建含多个 N（DNA/RNA）或 X（蛋白）的序列。
    MultipleX,
    /// 创建一个不含序列的占位符文件。
    Placeholder,
}

/// New File 对话框的状态。存在即表示对话框打开。
pub struct NewFileState {
    pub kind: MoleculeKind,
    pub mode: NewMode,
    /// 序列文本（原样保留用户的换行/格式；写文件前才清洗）。
    pub seq: String,
    /// `Create a sequence with multiple N's` 的个数。
    pub count: usize,
    /// 文件名称（含扩展名）。
    pub name: String,
    /// 链型（DNA/RNA；蛋白无此项）。
    pub double_stranded: bool,
    /// 拓扑（DNA/RNA；蛋白无此项）。
    pub circular: bool,
    /// `Reverse complement sequence`（DNA/RNA）。
    pub reverse_complement: bool,
    /// `Detect common features`（DNA/RNA；本实现未实现检测引擎，见文档）。
    pub detect_common: bool,
    /// `Transfer features, translation numbering, and amino acid colors`（蛋白）。
    pub transfer: bool,
    /// 预览标签页下标（`None` = 还没建）。
    pub preview_tab: Option<usize>,
    /// 上次构建预览用的签名（只在内容变化时重建文档）。
    pub sig: String,
    /// 用户是否手动改过文件名（改过就不再随类型自动换扩展名）。
    pub name_touched: bool,
}

impl NewFileState {
    pub fn new(kind: MoleculeKind) -> Self {
        Self {
            kind,
            mode: NewMode::Sequence,
            seq: String::new(),
            count: 1000,
            name: default_name(kind),
            // SnapGene 实测默认：DNA 双链 + 环状；RNA 单链 + 线性；蛋白无链型/拓扑
            double_stranded: kind != MoleculeKind::Rna,
            circular: kind == MoleculeKind::Dna,
            reverse_complement: false,
            detect_common: true,
            transfer: false,
            preview_tab: None,
            sig: String::new(),
            name_touched: false,
        }
    }

    /// 换分子类型：文件名未改过时同步换扩展名。
    pub fn set_kind(&mut self, kind: MoleculeKind) {
        if kind == self.kind {
            return;
        }
        if !self.name_touched {
            self.name = default_name(kind);
        }
        if kind == MoleculeKind::Protein {
            // 蛋白无链型/拓扑；切回来时恢复成 SnapGene 的默认
            self.detect_common = false;
        } else {
            self.double_stranded = kind != MoleculeKind::Rna;
            self.circular = kind == MoleculeKind::Dna;
        }
        self.kind = kind;
    }

    /// 对话框标题（也是「序列类型」下拉的当前项）。
    pub fn title(&self) -> &'static str {
        match self.kind {
            MoleculeKind::Dna => "New DNA File",
            MoleculeKind::Rna => "New RNA File",
            MoleculeKind::Protein => "New Protein File",
        }
    }

    /// 预览标签页的文档名（图谱标签用；SnapGene 在 Create 之前显示不带扩展名的名字）。
    pub fn preview_stem(&self) -> String {
        file_stem(&self.name).unwrap_or_else(|| "Untitled".to_string())
    }

    /// 序列文本是否像 FASTA / GenBank（SnapGene 会解析其中的特征）。
    fn looks_like_record(&self) -> bool {
        let t = self.seq.trim_start();
        t.starts_with('>') || t.starts_with("LOCUS") || t.starts_with("ID ") || t.starts_with("ID   ")
    }

    /// 写文件用的序列字节。
    fn sequence_bytes(&self) -> Vec<u8> {
        match self.mode {
            NewMode::Placeholder => Vec::new(),
            NewMode::MultipleX => {
                let n = self.count.max(1);
                let c = if self.kind == MoleculeKind::Protein {
                    b'X'
                } else {
                    b'N'
                };
                vec![c; n]
            }
            NewMode::Sequence => {
                if self.kind == MoleculeKind::Protein {
                    clean_protein(&self.seq)
                } else {
                    clean_nucleotides(&self.seq)
                }
            }
        }
    }

    /// 生成文件的字节（+ 特征数）。
    pub fn build_bytes(&self) -> Vec<u8> {
        let mut seq = self.sequence_bytes();
        // 粘进来的是 FASTA / GenBank 文本：取它的序列与特征（SnapGene 实测行为）
        let mut features: Vec<hx_formats::write::Feature> = Vec::new();
        if self.mode == NewMode::Sequence && self.looks_like_record() {
            if let Ok(recs) = hx_formats::import::parse(self.seq.as_bytes()) {
                if let Some(rec) = recs.first() {
                    seq = rec.sequence.clone();
                    features = rec
                        .features
                        .iter()
                        .filter_map(|f| {
                            let (a, b) = *f.segments.first()?;
                            Some(hx_formats::write::Feature {
                                start: a,
                                end: b,
                                name: f.name.clone(),
                                kind: if f.kind.is_empty() {
                                    "misc_feature".to_string()
                                } else {
                                    f.kind.clone()
                                },
                                color: "#ff0000".to_string(),
                                directionality: if f.reverse { 2 } else { 1 },
                                extra_segments: f.segments.iter().skip(1).copied().collect(),
                            })
                        })
                        .collect();
                }
            }
        }
        if self.reverse_complement && self.kind != MoleculeKind::Protein {
            seq = revcomp(&seq);
        }
        let mut flags = MoleculeFlags::from_byte(0);
        if self.kind != MoleculeKind::Protein {
            flags = flags
                .with_topology(if self.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                })
                .with_strandedness(if self.double_stranded {
                    Strandedness::Double
                } else {
                    Strandedness::Single
                });
        }
        NewFile {
            kind: self.kind,
            flags,
            sequence: seq,
            // Default map titles follow the file name. Do not pin a generated
            // "Untitled" custom label before the user chooses the save path.
            name: None,
            description: None,
            created_by: Some("AppimageGene".to_string()),
            features,
            extra_blocks: Vec::new(),
        }
        .to_bytes()
    }

    /// 预览重建签名（内容变了才重建文档）。
    pub fn signature(&self) -> String {
        format!(
            "{:?}|{:?}|{}|{}|{}|{}|{}|{}|{}",
            self.kind,
            self.mode,
            self.count,
            self.double_stranded,
            self.circular,
            self.reverse_complement,
            self.preview_stem(),
            self.seq.len(),
            self.seq
        )
    }
}

/// 新建文档的默认文件名（SnapGene 实测：`Untitled.dna` / `Untitled.rna` / `Untitled.prot`）。
pub fn default_name(kind: MoleculeKind) -> String {
    match kind {
        MoleculeKind::Dna => "Untitled.dna".to_string(),
        MoleculeKind::Rna => "Untitled.rna".to_string(),
        MoleculeKind::Protein => "Untitled.prot".to_string(),
    }
}

fn file_stem(name: &str) -> Option<String> {
    let t = name.trim();
    if t.is_empty() {
        return None;
    }
    Some(match t.rfind('.') {
        Some(i) if i > 0 => t[..i].to_string(),
        _ => t.to_string(),
    })
}

/// 核酸序列清洗：只保留 ASCII 字母（大写），丢弃空白/数字/符号。
fn clean_nucleotides(s: &str) -> Vec<u8> {
    s.bytes()
        .filter(|b| b.is_ascii_alphabetic())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// 蛋白序列清洗：保留 ASCII 字母与终止符 `*`。
fn clean_protein(s: &str) -> Vec<u8> {
    s.bytes()
        .filter(|b| b.is_ascii_alphabetic() || *b == b'*')
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            other => other,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 对话框 UI
// ---------------------------------------------------------------------------

/// 画 New File 对话框；返回 `Some(true)` = Create，`Some(false)` = Cancel。
pub fn draw(
    st: &mut NewFileState,
    ctx: &egui::Context,
    area: Option<egui::Rect>,
) -> Option<bool> {
    let mut done: Option<bool> = None;
    let mut open = true;
    // 贴文档区左上角（SnapGene 实测），没有文档区矩形时退回固定位置
    let pos = area
        .map(|r| r.left_top() + egui::vec2(10.0, 8.0))
        .unwrap_or(egui::pos2(360.0, 60.0));

    egui::Window::new("new_file")
        .title_bar(false)
        .resizable(false)
        .collapsible(false)
        .open(&mut open)
        .fixed_size(egui::vec2(DLG_W, 0.0))
        // SnapGene 的这个面板不能拖动：每帧钉在文档区左上角
        .current_pos(pos)
        .frame(
            egui::Frame::window(&ctx.style_of(egui::Theme::Light))
                .fill(egui::Color32::WHITE)
                .stroke(egui::Stroke::new(0.7, theme::BORDER))
                .inner_margin(egui::Margin::symmetric(12, 10)),
        )
        .show(ctx, |ui| {
            ui.set_width(DLG_W);
            // ---- 头行：标题本身是「序列类型」下拉（SnapGene 8 的 New File 面板）----
            ui.horizontal(|ui| {
                let cur = st.title();
                // oracle 实测：SnapGene 的标题下拉框高 ≈23 px、文字 cap 11 px
                ui.spacing_mut().interact_size.y = 24.0;
                egui::ComboBox::from_id_salt("new_file_kind")
                    .selected_text(RichText::new(crate::i18n::text(ctx,cur)).size(theme::FONT_UI + 1.0).strong())
                    .width(200.0)
                    .show_ui(ui, |ui| {
                        for k in [
                            MoleculeKind::Dna,
                            MoleculeKind::Rna,
                            MoleculeKind::Protein,
                        ] {
                            let label = match k {
                                MoleculeKind::Dna => "New DNA File",
                                MoleculeKind::Rna => "New RNA File",
                                MoleculeKind::Protein => "New Protein File",
                            };
                            if ui.selectable_label(st.kind == k, crate::i18n::text(ctx,label)).clicked() {
                                st.set_kind(k);
                            }
                        }
                    });
            });
            ui.add_space(11.0);

            // ---- ① Create a sequence ----
            let card = |ui: &mut egui::Ui, sel: bool, label: &str, f: &mut dyn FnMut(&mut egui::Ui)| {
                let label=crate::i18n::text(ctx,label);
                let stroke = if sel {
                    egui::Stroke::new(1.2, theme::ACCENT)
                } else {
                    egui::Stroke::new(1.0, theme::BORDER)
                };
                let resp = egui::Frame::new()
                    .fill(egui::Color32::WHITE)
                    .stroke(stroke)
                    .corner_radius(3.0)
                    .inner_margin(egui::Margin::symmetric(8, 4))
                    .show(ui, |ui| {
                        ui.set_width(DLG_W - 24.0);
                        ui.horizontal(|ui| {
                            // 整个卡片可点（点击即选中）；radio 本身只做指示
                            let _ = ui.radio(sel, label);
                        });
                        if sel {
                            f(ui);
                        }
                    })
                    .response;
                resp
            };

            let sel = st.mode == NewMode::Sequence;
            let r = card(ui, sel, "Create a sequence", &mut |ui| {
                ui.add_space(3.0);
                // 固定高度的滚动框（SnapGene 的序列框不随内容变高；oracle 实测框高 ≈58 px）
                egui::ScrollArea::vertical()
                    .max_height(SEQ_BOX_H)
                    .id_salt("new_file_seq")
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut st.seq)
                                .desired_width(DLG_W - 44.0)
                                .desired_rows(4)
                                .hint_text(seq_hint(st.kind))
                                .font(egui::TextStyle::Monospace),
                        );
                    });
                ui.add_space(3.0);
                if st.kind == MoleculeKind::Protein {
                    ui.add_enabled(
                        false,
                        egui::Checkbox::new(
                            &mut st.transfer,
                            crate::i18n::text(ctx,"Transfer features, translation numbering, and amino acid colors"),
                        ),
                    )
                    .on_hover_text(
                        crate::i18n::text(ctx,"Paste a region or feature copied from an open sequence to enable transfer"),
                    );
                } else {
                    ui.horizontal(|ui| {
                        ui.add(egui::Checkbox::new(
                            &mut st.reverse_complement,
                            crate::i18n::text(ctx,"Reverse complement sequence"),
                        ));
                        ui.add_space(10.0);
                        ui.add_enabled(
                            false,
                            egui::Checkbox::new(&mut st.detect_common, crate::i18n::text(ctx,"Detect common features")),
                        )
                        .on_hover_text(
                            crate::i18n::text(ctx,"Common feature detection needs SnapGene's .ftrs engine (not implemented)"),
                        );
                    });
                }
            });
            if r.clicked() {
                st.mode = NewMode::Sequence;
            }
            ui.add_space(11.0);

            // ---- ② Create a sequence with multiple N's / X's ----
            let sel = st.mode == NewMode::MultipleX;
            let letter = if st.kind == MoleculeKind::Protein {
                "X"
            } else {
                "N"
            };
            let r = card(
                ui,
                sel,
                &crate::ui_format!(ctx,"Create a sequence with multiple {}'s",letter),
                &mut |ui| {
                    ui.add_space(2.0);
                    ui.horizontal(|ui| {
                        ui.label(crate::i18n::text(ctx,"Number of residues:"));
                        ui.add(egui::DragValue::new(&mut st.count).range(1..=1_000_000));
                    });
                },
            );
            if r.clicked() {
                st.mode = NewMode::MultipleX;
            }
            ui.add_space(11.0);

            // ---- ③ Create a placeholder file with no sequence ----
            let sel = st.mode == NewMode::Placeholder;
            let r = card(
                ui,
                sel,
                "Create a placeholder file with no sequence",
                &mut |_ui| {},
            );
            if r.clicked() {
                st.mode = NewMode::Placeholder;
            }
            ui.add_space(14.0);

            // ---- 链型 / 拓扑（DNA/RNA；蛋白没有）----
            if st.kind != MoleculeKind::Protein {
                ui.label(RichText::new(crate::i18n::text(ctx,"Strandedness")).size(10.5).color(theme::TEXT_DIM));
                egui::Frame::new()
                    .stroke(egui::Stroke::new(1.0, theme::BORDER))
                    .corner_radius(3.0)
                    .inner_margin(egui::Margin::symmetric(8, 3))
                    .show(ui, |ui| {
                        ui.set_width(DLG_W - 24.0);
                        ui.horizontal(|ui| {
                            ui.radio_value(
                                &mut st.double_stranded,
                                true,
                                 crate::i18n::text(ctx,"Double-stranded"),
                            );
                            ui.add_space(24.0);
                            ui.radio_value(
                                &mut st.double_stranded,
                                false,
                                 crate::i18n::text(ctx,"Single-stranded"),
                            );
                        });
                    });
                ui.add_space(10.0);
                ui.label(RichText::new(crate::i18n::text(ctx,"Topology")).size(10.5).color(theme::TEXT_DIM));
                egui::Frame::new()
                    .stroke(egui::Stroke::new(1.0, theme::BORDER))
                    .corner_radius(3.0)
                    .inner_margin(egui::Margin::symmetric(8, 3))
                    .show(ui, |ui| {
                        ui.set_width(DLG_W - 24.0);
                        ui.horizontal(|ui| {
                            ui.radio_value(&mut st.circular, true, crate::i18n::text(ctx,"Circular"));
                            ui.add_space(64.0);
                            ui.radio_value(&mut st.circular, false, crate::i18n::text(ctx,"Linear"));
                        });
                    });
                ui.add_space(12.0);
            }

            // ---- 文件名称 ----
            ui.label(RichText::new(crate::i18n::text(ctx,"File Name")).size(10.5).color(theme::TEXT_DIM));
            egui::Frame::new()
                .stroke(egui::Stroke::new(1.0, theme::BORDER))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(8, 3))
                .show(ui, |ui| {
                    ui.set_width(DLG_W - 24.0);
                    let te = ui.add(
                        egui::TextEdit::singleline(&mut st.name)
                            .desired_width(DLG_W - 44.0)
                            .hint_text(crate::i18n::text(ctx,"Untitled.dna")),
                    );
                    if te.changed() {
                        st.name_touched = true;
                    }
                });
            ui.add_space(16.0);

            // oracle 实测：SnapGene 的按钮 ≈24 px 高（比正文行矮），右下角 Cancel/Create
            ui.horizontal(|ui| {
                ui.spacing_mut().interact_size.y = 22.0;
                ui.spacing_mut().button_padding = egui::vec2(10.0, 2.0);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    // 序列模式下序列框为空时不可创建（SnapGene：空序列的 Create 是灰的；
                    // 要建空文件请用「Create a placeholder file with no sequence」）
                    let can_create = match st.mode {
                        NewMode::Sequence => !st.seq.trim().is_empty(),
                        NewMode::MultipleX | NewMode::Placeholder => true,
                    };
                    // SnapGene 的禁用 Create 是**浅灰底**（不是浅蓝底）
                    let fill = if can_create {
                        theme::ACCENT
                    } else {
                        egui::Color32::from_rgb(0xd6, 0xd6, 0xd6)
                    };
                    let create = ui.add_enabled(
                        can_create,
                        egui::Button::new(RichText::new(crate::i18n::text(ctx,"Create")).color(egui::Color32::WHITE))
                            .fill(fill)
                            .min_size(egui::vec2(80.0, 24.0)),
                    );
                    if create.clicked() {
                        done = Some(true);
                    }
                    if ui
                        .add(egui::Button::new(crate::i18n::text(ctx,"Cancel")).min_size(egui::vec2(80.0, 24.0)))
                        .clicked()
                    {
                        done = Some(false);
                    }
                });
            });
            if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
                done = Some(false);
            }
            if ui.input(|i| i.key_pressed(egui::Key::Enter) && i.modifiers.command) {
                done = Some(true);
            }
        });
    if !open {
        done = Some(false);
    }
    done
}

const DLG_W: f32 = 470.0;
/// 序列输入框的固定高度（oracle 实测 SnapGene 的框高 ≈79 px，约 4 行等宽文本）。
const SEQ_BOX_H: f32 = 76.0;

/// 序列框的灰色占位提示（SnapGene 的空框里有一行提示文字）。
fn seq_hint(kind: MoleculeKind) -> &'static str {
    match kind {
        MoleculeKind::Dna => "Enter plain or formatted DNA sequence…",
        MoleculeKind::Rna => "Enter plain or formatted RNA sequence…",
        MoleculeKind::Protein => "Enter plain or formatted protein sequence…",
    }
}

/// 预览文档：由当前对话框状态构建（未保存、无路径）。
pub fn build_doc(st: &NewFileState) -> Result<Document, String> {
    Document::from_imported(st.build_bytes(), st.preview_stem())
}
