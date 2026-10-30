//! # hx-app — HelixBench egui/eframe 桌面壳
//!
//! V1 界面：
//! - **Map**：环状（外圈酶标签、环内刻度、特征扇环 + 方向箭头）与线性（轴、
//!   特征带、刻度、酶标签分层避让）两种图谱；
//! - **Sequence**：60 碱基/行、10 碱基/组的等宽序列视图，特征高亮 +
//!   酶位点下划线 + 点击定位；
//! - 侧栏：文件夹文件树（可折叠、路径持久化到 `~/.config/helixbench/settings`）
//!   与「其他打开的文件」；
//! - 打开方式：命令行参数、拖放、侧栏路径输入。
//!
//! 几何计算在 `hx-render`（与 egui 解耦，可复用于 SVG/PNG 导出）。

#![forbid(unsafe_code)]

pub mod app;
mod chrome;
pub mod clonedlg;
pub mod doc;
pub mod feat_card;
pub mod featuredlg;
pub mod fonts;
pub mod ggsettings;
pub mod icons;
mod info_panel;
pub mod newfile;
pub mod panels;
pub mod primer_tail;
pub mod primerdlg;
pub mod settings;
pub mod sidebar;
pub mod sysopen;
pub mod theme;
pub mod view_flat;
pub mod view_struct;
pub mod view_map;
pub mod view_seq;
pub mod viewer;

pub use app::{HelixApp, Tab};
mod history_view;
pub mod molecular_actions;
pub mod preferences;
pub mod branding;
mod rename;
mod enzyme_labels;
mod enzyme_manager;
mod native_view;
mod gel_view;
pub mod i18n;
