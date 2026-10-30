//! 「Golden Gate Settings…」窗口，布局照真 SnapGene 8.2.3 实拍：
//!
//! ```text
//!  ┌ Golden Gate Settings ─────────────────────────────────────┐
//!  │ Include in the list of Golden Gate enzymes:               │
//!  │ ┌ Standard Enzymes for Golden Gate Assembly ────────────┐ │
//!  │ │ ☐ AarI (PaqCI isoschizomer)   ☑ BbsI   ☑ BfuAI …    │ │
//!  │ └───────────────────────────────────────────────────────┘ │
//!  │ ┌ Additional Type IIS enzymes ──────────────────────────┐ │
//!  │ │ ☐ Acc36I  ☐ AcuI  ☐ Alw26I …                         │ │
//!  │ └───────────────────────────────────────────────────────┘ │
//!  │ Default Golden Gate enzyme: [BsaI ▾]                     │
//!  │ ☐ Choose the Golden Gate enzyme sets when simulating …    │
//!  │                                      [Reset] [OK] [Cancel]│
//!  └───────────────────────────────────────────────────────────┘
//! ```
//!
//! 勾选的 Standard / Additional 酶清单已经驱动新 Golden Gate 对话框的
//! 酶下拉（`GgSettings::enzyme_names`）；「enzyme sets」选项仍未实现，因此
//! 在界面上明确禁用。Golden Gate 模拟使用独立的 Type IIS 双链酶切与末端
//! 相容性模型，不再使用同源重叠拼接；具体差异见
//! `docs/cloning-equivalence.md`。

use eframe::egui::{self, RichText};

use crate::theme;

/// SnapGene「Standard Enzymes for Golden Gate Assembly」的实测清单与默认勾选。
pub const STANDARD: &[(&str, bool)] = &[
    ("AarI (PaqCI isoschizomer)", false),
    ("BbsI", true),
    ("BfuAI", true),
    ("BpiI (BbsI isoschizomer)", true),
    ("BsaI", true),
    ("BsmBI", false),
    ("BspMI", false),
    ("BspQI (SapI isoschizomer)", true),
    ("BtgZI", true),
    ("Eco31I (BsaI isoschizomer)", false),
    ("Esp3I (BsmBI isoschizomer)", true),
    ("LguI (SapI isoschizomer)", true),
    ("PaqCI", true),
    ("SapI", true),
];

/// SnapGene「Additional Type IIS enzymes」的首屏清单（可滚动）。
pub const ADDITIONAL: &[&str] = &[
    "Acc36I", "AcuI", "Alw26I", "BbvI", "BcoDI", "BfuAI", "BpmI", "BpuAI", "BsaXI",
    "BseRI", "BseXI", "BsgI", "BsmFI", "BsmI", "BspMI", "BspQI", "BtgZI", "BtsI",
    "EarI", "EciI", "Eco31I", "Eco57I", "Esp3I", "FokI", "HgaI", "HphI", "LguI",
    "MboII", "MmeI", "NmeDI", "PaqCI", "PleI", "SapI", "SfaNI", "TaqI", "Tth111I",
];

pub const DEFAULTS: &[&str] = &["BsaI", "BbsI", "BsmBI", "SapI", "PaqCI"];

#[derive(Debug, Clone, PartialEq)]
pub struct GgSettings {
    pub standard: Vec<(String, bool)>,
    pub additional: Vec<(String, bool)>,
    pub default_enzyme: usize,
    pub choose_sets: bool,
}

impl Default for GgSettings {
    fn default() -> Self {
        Self {
            standard: STANDARD
                .iter()
                .map(|(n, on)| ((*n).to_string(), *on))
                .collect(),
            additional: ADDITIONAL
                .iter()
                .map(|n| ((*n).to_string(), false))
                .collect(),
            default_enzyme: 0,
            choose_sets: false,
        }
    }
}

impl GgSettings {
    pub fn default_enzyme_name(&self) -> &str {
        DEFAULTS
            .get(self.default_enzyme)
            .copied()
            .unwrap_or("BsaI")
    }

    /// 勾选的 Golden Gate 酶清单，按设置窗口的顺序排列。
    ///
    /// 「Standard Enzymes」里带同工酶说明的名字（`BpiI (BbsI isoschizomer)`）
    /// 取括号前的酶名；酶库里不存在的名字直接跳过。全部取消勾选时退回
    /// [`DEFAULTS`]，否则克隆对话框没有可用的酶。
    pub fn enzyme_names(&self) -> Vec<&'static str> {
        let mut out: Vec<&'static str> = Vec::new();
        for (display, on) in self.standard.iter().chain(self.additional.iter()) {
            if !*on {
                continue;
            }
            let base = display
                .split(" (")
                .next()
                .unwrap_or(display)
                .trim();
            let Some(enzyme) = hx_enzymes::by_name(base) else {
                continue;
            };
            if !out.contains(&enzyme.name) {
                out.push(enzyme.name);
            }
        }
        if out.is_empty() {
            out = DEFAULTS
                .iter()
                .filter_map(|name| hx_enzymes::by_name(name).map(|e| e.name))
                .collect();
        }
        out
    }

    /// 从持久化设置恢复默认酶（未知名字保持当前值）。
    pub fn set_default_enzyme_name(&mut self, name: &str) {
        if let Some(i) = DEFAULTS.iter().position(|n| *n == name) {
            self.default_enzyme = i;
        }
    }

    fn reset(&mut self) {
        let d = GgSettings::default();
        self.standard = d.standard;
        self.additional = d.additional;
        self.default_enzyme = d.default_enzyme;
        self.choose_sets = d.choose_sets;
    }
}

/// 设置窗口的关闭结果：OK 保存，Cancel / 标题栏 × 放弃本次修改。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GgClose {
    Ok,
    Cancel,
}

/// 返回 `Some(Ok)` / `Some(Cancel)` 表示窗口已关闭，`None` 表示仍打开。
pub fn draw(s: &mut GgSettings, ctx: &egui::Context) -> Option<GgClose> {
    let mut open = true;
    let mut close = None;
    egui::Window::new(crate::i18n::text(ctx,"Golden Gate Settings"))
        .id(egui::Id::unique("gg_settings"))
        .resizable(false)
        .collapsible(false)
        .open(&mut open)
        .fixed_size(egui::vec2(540.0, 600.0))
        .order(egui::Order::Foreground)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.set_width(520.0);
            ui.label(
                RichText::new(crate::i18n::text(ctx,"Include in the list of Golden Gate enzymes:"))
                    .size(11.5)
                    .strong(),
            );
            ui.add_space(4.0);

            group(ui, "Standard Enzymes for Golden Gate Assembly", 172.0, |ui| {
                grid(ui, &mut s.standard);
            });
            ui.add_space(6.0);
            group(ui, "Additional Type IIS enzymes", 150.0, |ui| {
                grid(ui, &mut s.additional);
            });
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ctx,"Default Golden Gate enzyme:"));
                egui::ComboBox::from_id_salt("gg_default")
                    .selected_text(DEFAULTS[s.default_enzyme.min(DEFAULTS.len() - 1)])
                    .width(150.0)
                    .show_ui(ui, |ui| {
                        for (i, n) in DEFAULTS.iter().enumerate() {
                            ui.selectable_value(&mut s.default_enzyme, i, *n);
                        }
                    });
            });
            ui.add_space(4.0);
            // SnapGene 的酶组（enzyme sets）尚未实现；保留持久化值，但明确
            // 禁用，避免出现勾选后不影响模拟的假控件。
            ui.add_enabled(
                false,
                egui::Checkbox::new(
                    &mut s.choose_sets,
                    "Choose the Golden Gate enzyme sets when simulating Golden Gate Assembly",
                ),
            )
            .on_disabled_hover_text(crate::i18n::text(ctx,"Golden Gate enzyme sets are not implemented yet."));
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                ui.add_space((ui.available_width() - 210.0).max(0.0));
                if ui
                    .add_sized(egui::vec2(66.0, 24.0), egui::Button::new(crate::i18n::text(ctx,"Reset")))
                    .clicked()
                {
                    s.reset();
                }
                if ui
                    .add_sized(egui::vec2(66.0, 24.0), egui::Button::new(crate::i18n::text(ctx,"OK")))
                    .clicked()
                {
                    close = Some(GgClose::Ok);
                }
                if ui
                    .add_sized(egui::vec2(66.0, 24.0), egui::Button::new(crate::i18n::text(ctx,"Cancel")))
                    .clicked()
                {
                    close = Some(GgClose::Cancel);
                }
            });
        });
    if close.is_none() && !open {
        close = Some(GgClose::Cancel);
    }
    close
}

fn group(ui: &mut egui::Ui, title: &str, height: f32, body: impl FnOnce(&mut egui::Ui)) {
    egui::Frame::NONE
        .fill(egui::Color32::from_rgb(0xf7, 0xf7, 0xf7))
        .stroke(egui::Stroke::new(0.8, theme::BORDER))
        .inner_margin(egui::Margin::symmetric(6, 5))
        .show(ui, |ui| {
            ui.set_width(500.0);
            ui.label(RichText::new(crate::i18n::text(ui.ctx(),title)).size(10.8).strong());
            ui.add_space(2.0);
            egui::ScrollArea::vertical()
                .id_salt(title)
                .max_height(height)
                .show(ui, body);
        });
}

fn grid(ui: &mut egui::Ui, items: &mut [(String, bool)]) {
    let cols = 4;
    egui::Grid::new(format!("gg_grid_{}", items.len()))
        .num_columns(cols)
        .spacing(egui::vec2(6.0, 2.0))
        .show(ui, |ui| {
            for (i, (name, on)) in items.iter_mut().enumerate() {
                ui.checkbox(on, RichText::new(name.as_str()).size(10.3));
                if (i + 1) % cols == 0 {
                    ui.end_row();
                }
            }
        });
}
