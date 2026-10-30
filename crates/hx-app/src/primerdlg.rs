//! 「Add Primer… / Edit Primer…」与「Design Primers…（PCR）」对话框。
//!
//! 「Add Primer」对应 SnapGene 的 `Primers ▸ Add Primer…`（Ctrl+R）。流程按
//! SnapGene 8.2.3 实测复刻：
//!   1. 序列视图里有选区时先问 `Make a primer from the selection?`
//!      （Top Strand / Bottom Strand / Cancel，见 `crate::app`），
//!      选区**锁定**为结合位点；
//!   2. 打开引物编辑器：`Primer` / `Description` / `5′ [序列] 3′`，
//!      三个插入面板（Codon / Enzyme site / Short feature sequence）+ 各自
//!      的 `Insert` 按钮、`5′ Phosphorylated` 勾选、`Reverse Complement` 链接、
//!      统计行 `Primer: N-mer · X% GC` 与
//!      `Binding site: (1/1) a .. b · N annealed bases · Tm = nn°C`，
//!      以及引物在模板上的预览（5′ 尾巴单独画在相邻一条线上，与序列视图一致）。
//!
//! 3 个插入面板对应 SnapGene 的密码子 / 酶切位点 / 短特征序列：在序列框的
//! 光标处插入（没有光标就在 5′ 端），并进入本对话框的 Undo/Redo 栈。
//!
//! 「Design Primers」对应 SnapGene 的 `Actions ▸ PCR…`（Viewer 里是付费锁）。
//! 我们实现的是它的核心：给定目标区间，按长度/Tm/GC/产物长度筛选，枚举
//! 正向引物（5′ 端在目标起点 ±窗口内，3′ 端朝内）与反向引物（5′ 端在目标
//! 终点 ±窗口内，3′ 端朝内），按 Tm 配对打分排序，一键把两条引物写进文件。

use eframe::egui::{self, RichText};
use hx_algo::primer::{self, Anneal, DesignOptions, PrimerPair};
use hx_formats::annot::{PrimerAnnotation, Range};
use hx_formats::primer_edit::{PrimerEdit, PrimerSiteSpec, PrimerSpec};

use crate::theme;

/// 对话框宽度（SnapGene 的 Add Primer 页更宽；我们按窗口取。
const DLG_W: f32 = 560.0;

/// 退火判定的最少连续配对数（SnapGene `minContinuousMatchLen` 默认 10）。
const MIN_MATCH: usize = 10;

// ---------------------------------------------------------------------------// ---------------------------------------------------------------------------
// Add / Edit Primer
// ---------------------------------------------------------------------------

/// 对话框模式。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mode {
    Add,
    Edit { name: String },
}

/// 对话框结果。
pub enum Outcome {
    /// 施加编辑（Add 或 Replace）。
    Apply(Box<PrimerEdit>),
    /// 关闭对话框。
    Cancel,
}

/// 从选区建引物时锁定的结合位点。
///
/// 锁定的意义：SnapGene 里「选中一段 → Add Primer」问「Top Strand / Bottom
/// Strand」后，**这段选区就一直作为结合位点**（即使引物 5′ 端加了一堆不配对的
/// 尾巴）。用启发式扫描会重新找到别的位点，所以这里显式记住。
#[derive(Debug, Clone, PartialEq)]
pub struct ForcedSite {
    pub range: Range,
    pub strand: u8,
    /// 退火碱基数。
    pub annealed: usize,
    /// 退火段在 `sequence` 的**前缀**（3′ 端多打了碱基的引物，少见但真文件里有）。
    pub tail3: bool,
}

/// 「Add/Edit Primer」对话框状态；存在即表示打开。
pub struct PrimerDlg {
    pub mode: Mode,
    pub name: String,
    pub description: String,
    /// 用户输入的序列（可含 5′ 尾巴、空格/数字）。
    pub sequence: String,
    /// 当前选中的结合位点下标。
    pub picked: usize,
    /// 5′ 磷酸化（`<Primer phosphorylated="1">`）。
    pub phosphorylated: bool,
    /// 文档名（标题写成 `Add Primer for pUC19`，与 SnapGene 一致）。
    pub doc_name: String,
    /// 选区锁定的位点（`Make a primer from the selection` 之后）。
    pub forced: Option<ForcedSite>,
    /// 序列框的光标（插入面板插在这里）。
    caret: usize,
    /// 序列编辑历史（Undo / Redo）；连锁定位点一起存，Reverse Complement
    /// 会同时改链和坐标，撤销时得一起还原。
    undo: Vec<(String, Option<ForcedSite>)>,
    redo: Vec<(String, Option<ForcedSite>)>,
    last_seq: String,
    // ---- 三个插入面板的状态 ----
    codon_aa: usize,
    codon_idx: usize,
    codon_rc: bool,
    enzyme: usize,
    enzyme_rc: bool,
    enzyme_add6: bool,
    tag: usize,
    tag_rc: bool,
    /// 校验失败的提示（画在按钮旁）。
    error: Option<String>,
    // ---- 打开时的原值 ----
    orig_phos: bool,
    orig_date_added: String,
    orig_sequence: String,
    orig_description: String,
    orig_name: String,
    orig_sites: Vec<PrimerSiteSpec>,
}

impl PrimerDlg {
    pub fn add(doc_name: &str) -> Self {
        Self {
            mode: Mode::Add,
            name: "Primer 1".into(),
            description: String::new(),
            sequence: String::new(),
            picked: 0,
            phosphorylated: false,
            doc_name: doc_name.to_string(),
            forced: None,
            caret: 0,
            undo: Vec::new(),
            redo: Vec::new(),
            last_seq: String::new(),
            codon_aa: 0,
            codon_idx: 0,
            codon_rc: false,
            enzyme: 0,
            enzyme_rc: false,
            enzyme_add6: false,
            tag: 0,
            tag_rc: false,
            error: None,
            orig_phos: false,
            orig_date_added: String::new(),
            orig_sequence: String::new(),
            orig_description: String::new(),
            orig_name: String::new(),
            orig_sites: Vec::new(),
        }
    }

    /// 从序列视图的选区新建一条引物（SnapGene 的
    /// `Make a primer from the selection? → Top/Bottom Strand`）。
    ///
    /// `strand`：0 = Top Strand（引物序列 = 选区原样），1 = Bottom Strand
    /// （引物序列 = 选区的反向互补）。结合位点锁定为这段选区。
    pub fn from_selection(
        doc_name: &str,
        template: &[u8],
        start: usize,
        end: usize,
        strand: u8,
        name: &str,
    ) -> Option<Self> {
        let end = end.min(template.len());
        if start >= end {
            return None;
        }
        let picked: Vec<u8> = template[start..end].to_vec();
        let seq: String = if strand == 1 {
            primer::revcomp(&picked)
        } else {
            picked
        }
        .iter()
        .map(|b| *b as char)
        .collect();
        let mut d = Self::add(doc_name);
        d.name = name.to_string();
        d.sequence = seq.clone();
        d.last_seq = seq.clone();
        d.forced = Some(ForcedSite {
            range: Range::new(start, end),
            strand,
            annealed: end - start,
            tail3: false,
        });
        Some(d)
    }

    pub fn edit(p: &PrimerAnnotation, template: &[u8], circular: bool) -> Self {
        let sites: Vec<PrimerSiteSpec> = p
            .sites
            .iter()
            .map(|s| PrimerSiteSpec {
                range: s.range,
                strand: s.strand,
                annealed: annealed_of(p, s, template),
                tm: s
                    .melting_temperature
                    .unwrap_or_else(|| primer::tm(annealed_of(p, s, template).as_bytes())),
            })
            .collect();
        // 3′ 端多打了碱基的引物（退火段是 `sequence` 前缀）：位点无法靠
        // 「3′ 后缀匹配」搜到，直接把文件里的位点锁住，编辑框才显示得出来。
        let forced = p.sites.first().and_then(|s| {
            let tail3 = !crate::primer_tail::annealed_is_suffix(&p.sequence, s, template);
            tail3.then(|| ForcedSite {
                range: s.range,
                strand: s.strand,
                annealed: s.annealed_bases.min(p.sequence.len()),
                tail3: true,
            })
        });
        let _ = circular;
        Self {
            mode: Mode::Edit { name: p.name.clone() },
            name: p.name.clone(),
            description: p.description.clone(),
            sequence: p.sequence.clone(),
            picked: 0,
            phosphorylated: p.phosphorylated,
            doc_name: String::new(),
            forced,
            caret: 0,
            undo: Vec::new(),
            redo: Vec::new(),
            last_seq: p.sequence.clone(),
            codon_aa: 0,
            codon_idx: 0,
            codon_rc: false,
            enzyme: 0,
            enzyme_rc: false,
            enzyme_add6: false,
            tag: 0,
            tag_rc: false,
            error: None,
            orig_phos: p.phosphorylated,
            orig_date_added: p.date_added.clone(),
            orig_sequence: p.sequence.clone(),
            orig_description: p.description.clone(),
            orig_name: p.name.clone(),
            orig_sites: sites,
        }
    }

    /// 窗口标题（SnapGene 写 `Add Primer for <文档名>`）。
    pub fn title(&self, doc_name: &str) -> String {
        let doc = if doc_name.is_empty() { self.doc_name.as_str() } else { doc_name };
        match self.mode {
            Mode::Add => {
                if doc.is_empty() {
                    "Add Primer".into()
                } else {
                    format!("Add Primer for {doc}")
                }
            }
            Mode::Edit { .. } => {
                if doc.is_empty() {
                    "Edit Primer".into()
                } else {
                    format!("Edit Primer for {doc}")
                }
            }
        }
    }

    /// 清洗后的序列（只留字母、转大写）——与 SnapGene 粘贴的行为一致。
    pub fn clean_sequence(&self) -> String {
        self.sequence
            .bytes()
            .filter(|b| b.is_ascii_alphabetic())
            .map(|b| b.to_ascii_uppercase() as char)
            .collect()
    }

    /// 当前序列在模板上的退火位点（极大匹配，按长度降序）。
    ///
    /// 有锁定位点时只返回它（引物的 5′ 尾巴不参与匹配）。
    fn bindings(&self, template: &[u8], circular: bool) -> Vec<Anneal> {
        let seq = self.clean_sequence();
        if let Some(f) = &self.forced {
            return vec![Anneal {
                start: f.range.start,
                end: f.range.end,
                strand: f.strand,
                len: f.annealed.min(seq.len()),
            }];
        }
        if seq.len() < MIN_MATCH {
            return Vec::new();
        }
        let mut v = primer::anneal_bindings(template, seq.as_bytes(), circular, MIN_MATCH);
        v.sort_by(|a, b| b.len.cmp(&a.len).then(a.start.cmp(&b.start)));
        v
    }

    /// 退火段（默认是引物 3′ 后缀；锁定的 3′ 尾巴位点取前缀）。
    fn annealed_for(&self, a: &Anneal) -> String {
        let seq = self.clean_sequence();
        let take = a.len.min(seq.len());
        if self.forced.as_ref().map(|f| f.tail3).unwrap_or(false) {
            seq[..take].to_string()
        } else {
            seq[seq.len() - take..].to_string()
        }
    }

    /// 当前生效的结合位点（预览与统计行用）。
    fn active_site(&self, template: &[u8], circular: bool) -> Option<Anneal> {
        let v = self.bindings(template, circular);
        v.get(self.picked.min(v.len().saturating_sub(1))).cloned()
    }

    /// 生成编辑；序列未改动时返回 `None`（调用方不做任何事）。
    fn to_edit(&self, template: &[u8], circular: bool) -> Result<Option<PrimerEdit>, String> {
        let seq = self.clean_sequence();
        if self.name.trim().is_empty() {
            return Err("Primer name is required".into());
        }
        if seq.is_empty() {
            return Err("Primer sequence is required".into());
        }
        if !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            return Err("Primer sequence must be A/C/G/T".into());
        }
        // 序列没改：保留文件里原有的结合信息（含错配、多结合位点），只改名字/描述
        if matches!(self.mode, Mode::Edit { .. })
            && seq == self.orig_sequence
            && self.forced.is_none()
        {
            if self.name == self.orig_name
                && self.description == self.orig_description
                && self.phosphorylated == self.orig_phos
            {
                return Ok(None);
            }
            let spec = PrimerSpec {
                name: self.name.trim().to_string(),
                sequence: seq,
                description: self.description.trim().to_string(),
                date_added: self.orig_date_added.clone(),
                phosphorylated: self.phosphorylated,
                sites: self.orig_sites.clone(),
            };
            return Ok(Some(PrimerEdit::Replace {
                name: self.orig_name.clone(),
                spec: Box::new(spec),
            }));
        }
        let sites: Vec<Anneal> = self.bindings(template, circular);
        let Some(site) = sites.get(self.picked.min(sites.len().saturating_sub(1))) else {
            return Err("No binding site found in this sequence".into());
        };
        let spec = PrimerSpec {
            name: self.name.trim().to_string(),
            sequence: seq,
            description: self.description.trim().to_string(),
            date_added: self.orig_date_added.clone(),
            phosphorylated: self.phosphorylated,
            sites: vec![PrimerSiteSpec {
                range: Range::new(site.start, site.end),
                strand: site.strand,
                annealed: self.annealed_for(site),
                tm: primer::tm(self.annealed_for(site).as_bytes()),
            }],
        };
        Ok(Some(match &self.mode {
            Mode::Add => PrimerEdit::Add(Box::new(spec)),
            Mode::Edit { .. } => PrimerEdit::Replace {
                name: self.orig_name.clone(),
                spec: Box::new(spec),
            },
        }))
    }

    // ---- 序列编辑 ----

    fn push_undo(&mut self) {
        self.undo
            .push((self.last_seq.clone(), self.forced.clone()));
        if self.undo.len() > 64 {
            self.undo.remove(0);
        }
        self.redo.clear();
    }

    /// 在光标处插入一段碱基（SnapGene 的 Insert 按钮）。
    fn insert_at_caret(&mut self, text: &str) {
        self.push_undo();
        let chars: Vec<char> = self.sequence.chars().collect();
        let at = self.caret.min(chars.len());
        let mut out: String = chars[..at].iter().collect();
        out.push_str(text);
        out.extend(chars[at..].iter());
        self.sequence = out;
        self.caret = at + text.chars().count();
        self.last_seq = self.sequence.clone();
    }

    fn undo_edit(&mut self) {
        if let Some((prev, site)) = self.undo.pop() {
            self.redo
                .push((self.sequence.clone(), self.forced.clone()));
            self.sequence = prev.clone();
            self.forced = site;
            self.last_seq = prev;
        }
    }

    fn redo_edit(&mut self) {
        if let Some((next, site)) = self.redo.pop() {
            self.undo
                .push((self.sequence.clone(), self.forced.clone()));
            self.sequence = next.clone();
            self.forced = site;
            self.last_seq = next;
        }
    }

    /// Reverse Complement：序列反向互补，锁定的结合位点跟着翻链并镜像坐标。
    fn reverse_complement(&mut self, template_len: usize) {
        self.push_undo();
        let seq = self.clean_sequence();
        let rc = primer::revcomp(seq.as_bytes());
        self.sequence = rc.iter().map(|b| *b as char).collect();
        self.last_seq = self.sequence.clone();
        self.caret = 0;
        if let Some(f) = self.forced.as_mut() {
            f.strand = 1 - f.strand;
            f.tail3 = !f.tail3;
            let n = template_len;
            if n > 0 {
                f.range = Range::new(n.saturating_sub(f.range.end), n.saturating_sub(f.range.start));
            }
        }
    }

    /// 引物 5′ 端外侧的 6 个模板碱基（「Add 6 upstream bases」）。
    fn six_upstream(&self, template: &[u8], circular: bool) -> Option<String> {
        let site = self.active_site(template, circular)?;
        let n = template.len();
        if n == 0 {
            return None;
        }
        let mut out: Vec<u8> = Vec::new();
        if site.strand == 0 {
            for k in 1..=6 {
                let idx = if site.start >= k {
                    site.start - k
                } else if circular {
                    (n + site.start - k) % n
                } else {
                    return None;
                };
                out.push(template[idx]);
            }
            out.reverse();
            Some(out.iter().map(|b| *b as char).collect())
        } else {
            for k in 0..6 {
                let idx = if site.end + k < n {
                    site.end + k
                } else if circular {
                    (site.end + k) % n
                } else {
                    return None;
                };
                out.push(template[idx]);
            }
            let rc = primer::revcomp(&out);
            Some(rc.iter().map(|b| *b as char).collect())
        }
    }
}

/// 从文件里的引物模型取退火段：长度 = 位点碱基数，端别由模板判定
/// （5′ 尾巴 = `sequence` 后缀退火；末尾多打的碱基 = 前缀退火）。
fn annealed_of(
    p: &PrimerAnnotation,
    s: &hx_formats::annot::PrimerSite,
    template: &[u8],
) -> String {
    let n = p.sequence.len();
    let take = s.annealed_bases.min(n);
    if take == 0 || n == 0 {
        return p.sequence.clone();
    }
    if crate::primer_tail::annealed_is_suffix(&p.sequence, s, template) {
        p.sequence[n - take..].to_string()
    } else {
        p.sequence[..take].to_string()
    }
}

/// 画对话框；返回 `Some(_)` 表示对话框应当关闭。
pub fn draw(
    dlg: &mut PrimerDlg,
    ctx: &egui::Context,
    template: &[u8],
    circular: bool,
) -> Option<Outcome> {
    let mut out: Option<Outcome> = None;
    let mut open = true;
    let seq = dlg.clean_sequence();
    let bindings = dlg.bindings(template, circular);
    let len = seq.len();
    let gc = primer::gc_percent(seq.as_bytes());
    let site = dlg.active_site(template, circular);
    let width = DLG_W;

    egui::Window::new("primer_dialog")
        .id(egui::Id::unique("primer_dialog"))
        .title_bar(false)
        .resizable(false)
        .collapsible(false)
        .open(&mut open)
        .fixed_size(egui::vec2(width, 0.0))
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -20.0))
        .frame(
            egui::Frame::window(&ctx.style_of(egui::Theme::Light))
                .fill(egui::Color32::WHITE)
                .stroke(egui::Stroke::new(0.8, theme::BORDER))
                .inner_margin(egui::Margin::symmetric(12, 10)),
        )
        .show(ctx, |ui| {
            ui.set_width(width);
            ui.horizontal(|ui| {
                let template=match dlg.mode{Mode::Add=>"Add Primer for {}",Mode::Edit{..}=>"Edit Primer for {}"};
                ui.label(RichText::new(crate::i18n::format_values(ctx,template,&[dlg.doc_name.clone()])).size(13.5).strong());
            });
            ui.separator();
            ui.add_space(2.0);

            // ---- 名称 / 描述 ----
            ui.horizontal(|ui| {
                ui.label(RichText::new(crate::i18n::text(ctx,"Primer")).size(11.5).color(theme::TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut dlg.name)
                        .desired_width(180.0)
                        .hint_text(crate::i18n::text(ctx,"Primer name")),
                );
                ui.add_space(10.0);
                ui.label(RichText::new(crate::i18n::text(ctx,"Description")).size(11.5).color(theme::TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut dlg.description)
                        .desired_width(width - 330.0)
                        .hint_text(crate::i18n::text(ctx,"Optional description")),
                );
            });

            // ---- 5′ [序列] 3′ ----
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.label(RichText::new("5′").size(11.5).color(theme::TEXT_DIM));
                let out = egui::TextEdit::singleline(&mut dlg.sequence)
                    .id(egui::Id::unique("primer_seq_edit"))
                    .desired_width(width - 70.0)
                    .font(egui::TextStyle::Monospace)
                    .hint_text(crate::i18n::text(ctx,"ACGT…  (5′ tail may not anneal)"))
                    .show(ui);
                if let Some(cr) = out.cursor_range {
                    dlg.caret = cr.primary.index.into();
                }
                if out.response.changed() {
                    dlg.push_undo();
                }
                dlg.last_seq = dlg.sequence.clone();
                ui.label(RichText::new("3′").size(11.5).color(theme::TEXT_DIM));
            });

            // ---- 三个插入面板 ----
            ui.add_space(6.0);
            let mut insert: Option<String> = None;
            ui.columns(3, |cols| {
                // Codon
                egui::Frame::group(cols[0].style())
                    .fill(theme::PANEL_ALT)
                    .show(&mut cols[0], |ui| {
                        ui.set_width(width / 3.0 - 26.0);
                        ui.label(RichText::new(crate::i18n::text(ctx,"Codon")).size(11.5).strong());
                        ui.add_space(2.0);
                        ui.horizontal(|ui| {
                            let aa_txt = format!(
                                "{} ({})",
                                primer::CODON_TABLE[dlg.codon_aa].aa3,
                                primer::CODON_TABLE[dlg.codon_aa].aa1
                            );
                            egui::ComboBox::from_id_salt("codon_aa")
                                .selected_text(aa_txt)
                                .width(96.0)
                                .show_ui(ui, |ui| {
                                    for (i, g) in primer::CODON_TABLE.iter().enumerate() {
                                        if ui
                                            .selectable_value(&mut dlg.codon_aa, i, g.aa3)
                                            .clicked()
                                        {
                                            dlg.codon_idx = 0;
                                        }
                                    }
                                });
                            let codons = primer::CODON_TABLE[dlg.codon_aa].codons;
                            let ci = dlg.codon_idx.min(codons.len() - 1);
                            egui::ComboBox::from_id_salt("codon_choice")
                                .selected_text(codons[ci])
                                .width(66.0)
                                .show_ui(ui, |ui| {
                                    for (i, c) in codons.iter().enumerate() {
                                        ui.selectable_value(&mut dlg.codon_idx, i, *c);
                                    }
                                });
                        });
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut dlg.codon_rc, RichText::new(crate::i18n::text(ctx,"Reverse")).size(11.0));
                            if ui.button(crate::i18n::text(ui.ctx(),"Insert")).clicked() {
                                let codons = primer::CODON_TABLE[dlg.codon_aa].codons;
                                let c = codons[dlg.codon_idx.min(codons.len() - 1)].as_bytes();
                                let text: String = if dlg.codon_rc {
                                    primer::revcomp(c).iter().map(|b| *b as char).collect()
                                } else {
                                    c.iter().map(|b| *b as char).collect()
                                };
                                insert = Some(text);
                            }
                        });
                    });

                // Enzyme site
                egui::Frame::group(cols[1].style())
                    .fill(theme::PANEL_ALT)
                    .show(&mut cols[1], |ui| {
                        ui.set_width(width / 3.0 - 26.0);
                        ui.label(RichText::new(crate::i18n::text(ctx,"Enzyme site")).size(11.5).strong());
                        ui.add_space(2.0);
                        let enz = hx_enzymes::ENZYMES
                            .get(dlg.enzyme.min(hx_enzymes::ENZYMES.len() - 1));
                        let enz_name = enz.map(|e| e.name).unwrap_or("—");
                        ui.horizontal(|ui| {
                            egui::ComboBox::from_id_salt("enz_choice")
                                .selected_text(enz_name)
                                .width(120.0)
                                .show_ui(ui, |ui| {
                                    egui::ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                                        for (i, e) in hx_enzymes::ENZYMES.iter().enumerate() {
                                            ui.selectable_value(&mut dlg.enzyme, i, e.name);
                                        }
                                    });
                                });
                        });
                        ui.horizontal(|ui| {
                            let site_txt = enz.map(|e| e.site).unwrap_or("—");
                            egui::ComboBox::from_id_salt("enz_site")
                                .selected_text(site_txt)
                                .width(110.0)
                                .show_ui(ui, |ui| {
                                    ui.label(site_txt);
                                });
                        });
                        ui.horizontal(|ui| {
                            if dlg.six_upstream(template, circular).is_some() {
                                ui.checkbox(
                                    &mut dlg.enzyme_add6,
                                    RichText::new(crate::i18n::text(ctx,"Add 6 upstream bases")).size(11.0),
                                );
                            } else {
                                ui.add_enabled(
                                    false,
                                    egui::Checkbox::new(&mut dlg.enzyme_add6, "Add 6 upstream bases"),
                                );
                            }
                        });
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut dlg.enzyme_rc, RichText::new(crate::i18n::text(ctx,"Reverse")).size(11.0));
                            if ui.button(crate::i18n::text(ui.ctx(),"Insert")).clicked() {
                                let site_seq = enz.map(|e| e.site).unwrap_or("");
                                let mut text = String::new();
                                if dlg.enzyme_add6 {
                                    if let Some(ctx6) = dlg.six_upstream(template, circular) {
                                        text.push_str(&ctx6);
                                    }
                                }
                                text.push_str(&site_seq.to_ascii_uppercase());
                                if dlg.enzyme_rc {
                                    let rc = primer::revcomp(text.as_bytes());
                                    text = rc.iter().map(|b| *b as char).collect();
                                }
                                insert = Some(text);
                            }
                        });
                    });

                // Short feature sequence
                egui::Frame::group(cols[2].style())
                    .fill(theme::PANEL_ALT)
                    .show(&mut cols[2], |ui| {
                        ui.set_width(width / 3.0 - 26.0);
                        ui.label(RichText::new(crate::i18n::text(ctx,"Short feature sequence")).size(11.5).strong());
                        ui.add_space(2.0);
                        let t = primer::PEPTIDES
                            .get(dlg.tag.min(primer::PEPTIDES.len() - 1));
                        let tag_name = t.map(|p| p.name).unwrap_or("—");
                        egui::ComboBox::from_id_salt("tag_choice")
                            .selected_text(tag_name)
                            .width(width / 3.0 - 40.0)
                            .show_ui(ui, |ui| {
                                let mut last_group = "";
                                for (i, p) in primer::PEPTIDES.iter().enumerate() {
                                    if p.group != last_group {
                                        if !last_group.is_empty() {
                                            ui.separator();
                                        }
                                        ui.label(
                                            RichText::new(p.group)
                                                .size(10.5)
                                                .color(theme::TEXT_DIM),
                                        );
                                        last_group = p.group;
                                    }
                                    ui.selectable_value(&mut dlg.tag, i, p.name);
                                }
                            });
                        let pep = t.map(|p| p.peptide).unwrap_or("");
                        ui.label(
                            RichText::new(pep).monospace().size(10.5).color(theme::TEXT_DIM),
                        );
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut dlg.tag_rc, RichText::new(crate::i18n::text(ctx,"Reverse")).size(11.0));
                            if ui.button(crate::i18n::text(ui.ctx(),"Insert")).clicked() {
                                let dna = primer::dna_for_peptide(pep);
                                let text: String = if dlg.tag_rc {
                                    primer::revcomp(dna.as_bytes())
                                        .iter()
                                        .map(|b| *b as char)
                                        .collect()
                                } else {
                                    dna
                                };
                                insert = Some(text);
                            }
                        });
                    });
            });
            if let Some(text) = insert {
                dlg.insert_at_caret(&text);
            }

            // ---- 5′ 磷酸化 / 反向互补 ----
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.checkbox(
                    &mut dlg.phosphorylated,
                    RichText::new(crate::i18n::text(ctx,"5′ Phosphorylated")).size(11.5),
                );
                ui.add_space(16.0);
                if ui
                    .link(RichText::new(crate::i18n::text(ctx,"Reverse Complement")).size(11.5))
                    .clicked()
                {
                    dlg.reverse_complement(template.len());
                }
            });

            // ---- 统计行 ----
            ui.add_space(3.0);
            ui.horizontal(|ui| {
                let stats = if len == 0 {
                    crate::ui_format!(ctx,"Primer:  {}-mer",0)
                } else {
                    crate::ui_format!(ctx,"Primer:  {}-mer  ·  {}% GC",len,format!("{gc:.0}"))
                };
                ui.label(RichText::new(stats).size(11.5));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(a) = &site {
                        let ann = dlg.annealed_for(a);
                        let tm = primer::tm(ann.as_bytes());
                        let (s, e) = if a.start < a.end {
                            (a.start + 1, a.end)
                        } else {
                            (a.start + 1, a.end)
                        };
                        let which = if bindings.len() > 1 {
                            format!("({}/{})", dlg.picked + 1, bindings.len())
                        } else {
                            "(1/1)".to_string()
                        };
                        ui.label(
                            RichText::new(crate::ui_format!(ctx,
                                "Binding site: {}  {} .. {}  ·  {} annealed bases  ·  Tm = {}°C",
                                which,s,e,ann.len(),format!("{tm:.0}")
                            ))
                            .size(11.5),
                        );
                    } else if len < MIN_MATCH {
                        ui.label(
                            RichText::new(crate::ui_format!(ctx,"Need ≥ {} bases to bind.",MIN_MATCH))
                                .size(11.5)
                                .color(theme::TEXT_DIM),
                        );
                    } else {
                        ui.label(
                            RichText::new(crate::i18n::text(ctx,"No binding site in this sequence."))
                                .size(11.5)
                                .color(theme::SITE_HOVER),
                        );
                    }
                });
            });

            // ---- 其它结合位点（多命中时可选）----
            if bindings.len() > 1 {
                ui.add_space(2.0);
                ui.horizontal_wrapped(|ui| {
                    for (i, a) in bindings.iter().take(6).enumerate() {
                        let label = format!(
                            "{}..{} {}",
                            a.start + 1,
                            a.end.max(a.start + a.len),
                            if a.strand == 0 { "forward" } else { "reverse" }
                        );
                        if ui
                            .selectable_label(dlg.picked == i, RichText::new(label).size(11.0))
                            .clicked()
                        {
                            dlg.picked = i;
                        }
                    }
                    if bindings.len() > 6 {
                        ui.label(
                        RichText::new(crate::ui_format!(ctx,"+{} more", bindings.len() - 6))
                                .size(11.0)
                                .color(theme::TEXT_DIM),
                        );
                    }
                });
            }

            // ---- 引物预览 ----
            ui.add_space(6.0);
            draw_preview(ui, &seq, dlg.phosphorylated, template, site.as_ref(), width);

            // ---- 底部按钮 ----
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(!dlg.undo.is_empty(), egui::Button::new(crate::i18n::text(ctx,"Undo")))
                    .clicked()
                {
                    dlg.undo_edit();
                }
                if ui
                    .add_enabled(!dlg.redo.is_empty(), egui::Button::new(crate::i18n::text(ctx,"Redo")))
                    .clicked()
                {
                    dlg.redo_edit();
                }
                if matches!(dlg.mode, Mode::Edit { .. })
                    && ui.button(crate::i18n::text(ui.ctx(),"Remove Primer")).clicked()
                {
                    out = Some(Outcome::Apply(Box::new(PrimerEdit::Remove {
                        name: dlg.orig_name.clone(),
                    })));
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let can_add = !bindings.is_empty() && !seq.is_empty();
                    let add_label = match dlg.mode {
                        Mode::Add => "Add Primer to Template",
                        Mode::Edit { .. } => "Save Changes",
                    };
                    if ui
                        .add_enabled(can_add, egui::Button::new(add_label))
                        .clicked()
                    {
                        match dlg.to_edit(template, circular) {
                            Ok(Some(e)) => out = Some(Outcome::Apply(Box::new(e))),
                            Ok(None) => out = Some(Outcome::Cancel),
                            Err(e) => dlg.error = Some(e),
                        }
                    }
                    if ui.button(crate::i18n::text(ui.ctx(),"Close")).clicked() {
                        out = Some(Outcome::Cancel);
                    }
                    if let Some(e) = &dlg.error {
                        ui.label(RichText::new(crate::i18n::message(ctx,e)).size(11.5).color(theme::SITE_HOVER));
                    }
                });
            });
        });
    if !open {
        return Some(Outcome::Cancel);
    }
    out
}

/// 引物预览：模板窗口 + 退火段 + 5′ 尾巴（不配对的那段单独画在相邻一条线上，
/// 与真 SnapGene 的序列视图一致）。
fn draw_preview(
    ui: &mut egui::Ui,
    seq: &str,
    phosphorylated: bool,
    template: &[u8],
    site: Option<&Anneal>,
    width: f32,
) {
    const N: usize = 60;
    let h = 108.0;
    let (rect, _resp) = ui.allocate_exact_size(egui::vec2(width, h), egui::Sense::hover());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 2.0, egui::Color32::WHITE);
    p.rect_stroke(
        rect.shrink(0.5),
        2.0,
        egui::Stroke::new(1.0, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    if template.is_empty() {
        p.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no template",
            egui::FontId::proportional(11.0),
            theme::TEXT_DIM,
        );
        return;
    }
    let n = template.len();
    // 窗口起点：以结合位点为中心
    let start = match site {
        Some(a) if a.start < a.end => a.start.saturating_sub((N - (a.end - a.start)) / 2),
        _ => 0,
    };
    let start = start.min(n.saturating_sub(1));
    let base_w = (width - 24.0) / N as f32;
    let font = egui::FontId::monospace((base_w * 1.35).clamp(9.0, 13.0));
    let x0 = rect.left() + 12.0;
    let y_top = rect.top() + 34.0;
    let y_bot = y_top + 18.0;

    let at = |i: usize| -> u8 { template[(start + i) % n] };
    let top: String = (0..N)
        .map(|i| if start + i < n { at(i) as char } else { ' ' })
        .collect();
    let bot: String = (0..N)
        .map(|i| {
            if start + i < n {
                primer::complement(at(i)) as char
            } else {
                ' '
            }
        })
        .collect();
    let gal_top = p.layout_no_wrap(top.clone(), font.clone(), theme::SEQ_LINE);
    p.galley(egui::Pos2::new(x0, y_top), gal_top, theme::SEQ_LINE);
    let gal_bot = p.layout_no_wrap(bot.clone(), font.clone(), theme::BASE_COMPLEMENT);
    p.galley(egui::Pos2::new(x0, y_bot), gal_bot, theme::BASE_COMPLEMENT);

    // 标尺（每 10 个碱基一个数字）
    for i in (0..N).step_by(10) {
        if start + i >= n {
            break;
        }
        p.text(
            egui::Pos2::new(x0 + i as f32 * base_w, y_top - 13.0),
            egui::Align2::LEFT_TOP,
            format!("{}", start + i + 1),
            egui::FontId::proportional(8.5),
            theme::TEXT_DIM,
        );
    }

    // 结合区高亮
    if let Some(a) = site {
        if a.start < a.end && a.end <= n && a.start >= start && a.end <= start + N {
            let x1 = x0 + (a.start - start) as f32 * base_w;
            let x2 = x0 + (a.end - start) as f32 * base_w;
            p.rect_filled(
                egui::Rect::from_min_max(
                    egui::Pos2::new(x1, y_top - 2.0),
                    egui::Pos2::new(x2, y_bot + 15.0),
                ),
                0.0,
                egui::Color32::from_rgb(0xdc, 0xe9, 0xf8),
            );
        }
    }

    // ---- 引物本体 ----
    let Some(a) = site else {
        p.text(
            egui::Pos2::new(x0, y_bot + 22.0),
            egui::Align2::LEFT_TOP,
            crate::i18n::text(ui.ctx(),"No binding site — the primer will not bind this sequence."),
            egui::FontId::proportional(10.5),
            theme::SITE_HOVER,
        );
        return;
    };
    let overhang = seq.len().saturating_sub(a.len.min(seq.len()));
    let primer_color = theme::PRIMER;
    let bar_h = 11.0;
    let step = bar_h * 0.45;
    // 退火段在模板上的格子
    let ax1 = x0 + a.start.saturating_sub(start) as f32 * base_w;
    let ax2 = x0 + a.end.saturating_sub(start) as f32 * base_w;
    let (ax1, ax2) = if a.start < a.end && a.end <= n && a.start >= start {
        (ax1, ax2.max(ax1 + base_w))
    } else {
        // 位点在窗口外：画在下方居中
        let w = a.len as f32 * base_w;
        let c = rect.center().x;
        (c - w * 0.5, c + w * 0.5)
    };
    let bar_y = y_bot + 20.0;
    // 5′ 尾巴（不配对）：正向画在左上、反向画在右下，与 SnapGene 一致
    let tail_w = overhang as f32 * base_w;
    let tail_src = &seq[..overhang.min(seq.len())];
    let tail: String = if a.strand == 0 {
        tail_src.to_string()
    } else {
        tail_src.chars().rev().collect()
    };
    let (tail_rect, body_rect) = if a.strand == 0 {
        (
            egui::Rect::from_min_max(
                egui::Pos2::new(ax1 - tail_w, bar_y - step),
                egui::Pos2::new(ax1, bar_y - step + bar_h),
            ),
            egui::Rect::from_min_max(
                egui::Pos2::new(ax1, bar_y),
                egui::Pos2::new(ax2, bar_y + bar_h),
            ),
        )
    } else {
        (
            egui::Rect::from_min_max(
                egui::Pos2::new(ax2, bar_y + step),
                egui::Pos2::new(ax2 + tail_w, bar_y + step + bar_h),
            ),
            egui::Rect::from_min_max(
                egui::Pos2::new(ax1, bar_y),
                egui::Pos2::new(ax2, bar_y + bar_h),
            ),
        )
    };
    for r in [body_rect, tail_rect] {
        if r.width() < 0.5 {
            continue;
        }
        p.rect_filled(r, 0.0, egui::Color32::WHITE);
        p.rect_stroke(
            r,
            0.0,
            egui::Stroke::new(1.2, primer_color),
            egui::StrokeKind::Inside,
        );
    }
    // 5′ 尾巴碱基
    if overhang > 0 {
        let gal = p.layout_no_wrap(tail.clone(), font.clone(), theme::TEXT);
        if gal.size().x <= tail_rect.width() + 2.0 {
            p.galley(
                egui::Pos2::new(tail_rect.left(), tail_rect.center().y - gal.size().y * 0.5),
                gal,
                theme::TEXT,
            );
        }
    }
    // 退火段碱基（与序列视图一致：画的是**显示链**的文字）
    if overhang < seq.len() {
        let ann = &seq[overhang..];
        let shown: String = if a.strand == 0 {
            ann.to_string()
        } else {
            ann.chars().rev().collect()
        };
        let gal = p.layout_no_wrap(shown, font.clone(), theme::TEXT);
        if gal.size().x <= body_rect.width() + 2.0 {
            p.galley(
                egui::Pos2::new(
                    body_rect.left() + (body_rect.width() - gal.size().x).max(0.0) * 0.5,
                    body_rect.center().y - gal.size().y * 0.5,
                ),
                gal,
                theme::TEXT,
            );
        }
    }
    // 3′ 箭头
    let (tip_x, dir) = if a.strand == 0 {
        (body_rect.right(), 1.0f32)
    } else {
        (body_rect.left(), -1.0f32)
    };
    let ty = body_rect.center().y;
    p.add(egui::Shape::line(
        vec![
            egui::Pos2::new(tip_x - dir * 7.0, ty - 6.0),
            egui::Pos2::new(tip_x, ty),
            egui::Pos2::new(tip_x - dir * 7.0, ty + 6.0),
        ],
        egui::Stroke::new(1.2, primer_color),
    ));
    // 5′ 磷酸化
    if phosphorylated {
        let c = if a.strand == 0 {
            egui::Pos2::new(body_rect.left() - 9.0, body_rect.center().y - step)
        } else {
            egui::Pos2::new(body_rect.right() + 9.0, body_rect.center().y + step)
        };
        p.circle_stroke(c, 6.0, egui::Stroke::new(1.2, primer_color));
        p.text(
            c,
            egui::Align2::CENTER_CENTER,
            "P",
            egui::FontId::proportional(9.0),
            primer_color,
        );
    }
    p.text(
        egui::Pos2::new(body_rect.center().x, body_rect.top() - 3.0),
        egui::Align2::CENTER_BOTTOM,
        format!("{} bp primer", seq.len()),
        egui::FontId::proportional(9.5),
        theme::TEXT_DIM,
    );
}

// ---------------------------------------------------------------------------
// Design Primers (PCR)
// ---------------------------------------------------------------------------

/// 设计对话框结果。
pub enum DesignOutcome {
    /// 把两条引物加入文件（正向、反向）。
    Add(Box<PrimerPair>),
    Cancel,
}

/// 「Design Primers」对话框状态。
pub struct DesignDlg {
    /// 目标区间（1-based 闭区间，界面用）。
    pub start: usize,
    pub end: usize,
    pub opts: DesignOptions,
    pub results: Vec<PrimerPair>,
    pub picked: usize,
    pub has_run: bool,
}

impl DesignDlg {
    pub fn new(len: usize, selection: Option<(usize, usize)>) -> Self {
        let (s, e) = match selection {
            Some((a, b)) if b > a => (a + 1, b),
            _ => (1, len),
        };
        Self {
            start: s,
            end: e.max(s),
            opts: DesignOptions {
                max_product: len.max(100),
                ..Default::default()
            },
            results: Vec::new(),
            picked: 0,
            has_run: false,
        }
    }

    pub fn run(&mut self, template: &[u8], circular: bool) {
        let len = template.len();
        let (s, e) = (self.start.min(len), self.end.min(len));
        let region = if s <= e {
            Range::new(s.saturating_sub(1), e)
        } else {
            Range::new(s.saturating_sub(1), e)
        };
        self.results = primer::design_pair(template, region, circular, &self.opts);
        self.picked = 0;
        self.has_run = true;
    }

    pub fn picked_pair(&self) -> Option<&PrimerPair> {
        self.results.get(self.picked.min(self.results.len().saturating_sub(1)))
    }
}

/// 画「Design Primers」；返回 `Some(_)` 表示关闭。
pub fn draw_design(
    dlg: &mut DesignDlg,
    ctx: &egui::Context,
    template: &[u8],
    circular: bool,
) -> Option<DesignOutcome> {
    let mut out: Option<DesignOutcome> = None;
    let mut open = true;
    let len = template.len();

    egui::Window::new("design_primers")
        .id(egui::Id::unique("design_primers"))
        .title_bar(false)
        .resizable(false)
        .collapsible(false)
        .open(&mut open)
        .fixed_size(egui::vec2(DLG_W + 120.0, 0.0))
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -30.0))
        .frame(
            egui::Frame::window(&ctx.style_of(egui::Theme::Light))
                .fill(egui::Color32::WHITE)
                .stroke(egui::Stroke::new(0.8, theme::BORDER))
                .inner_margin(egui::Margin::symmetric(12, 10)),
        )
        .show(ctx, |ui| {
            ui.set_width(DLG_W + 120.0);
            ui.horizontal(|ui| {
                ui.label(RichText::new(crate::i18n::text(ctx,"Design Primers")).size(13.5).strong());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                    RichText::new(crate::ui_format!(ctx,"{} bp template", crate::app::nfmt(len)))
                            .size(11.0)
                            .color(theme::TEXT_DIM),
                    );
                });
            });
            ui.separator();

            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ui.ctx(),"Target region (1-based):"));
                ui.add(
                    egui::DragValue::new(&mut dlg.start)
                        .range(1..=len.max(1))
                        .speed(1.0),
                );
                ui.label("..");
                ui.add(
                    egui::DragValue::new(&mut dlg.end)
                        .range(1..=len.max(1))
                        .speed(1.0),
                );
                ui.add_space(10.0);
                ui.label(crate::i18n::text(ui.ctx(),"Product:"));
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.min_product)
                        .range(1..=1_000_000)
                        .speed(5.0),
                );
                ui.label("..");
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.max_product)
                        .range(1..=1_000_000)
                        .speed(5.0),
                );
                ui.label(crate::i18n::text(ctx,"bp"));
            });
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ui.ctx(),"Primer length:"));
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.min_len)
                        .range(10..=40)
                        .speed(0.2),
                );
                ui.label("..");
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.max_len)
                        .range(10..=40)
                        .speed(0.2),
                );
                ui.label(crate::i18n::text(ctx,"nt"));
                ui.add_space(10.0);
                ui.label(crate::i18n::text(ui.ctx(),"Tm:"));
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.min_tm)
                        .range(0.0..=90.0)
                        .speed(0.5),
                );
                ui.label("..");
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.max_tm)
                        .range(0.0..=90.0)
                        .speed(0.5),
                );
                ui.label("°C");
            });
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                ui.label(crate::i18n::text(ui.ctx(),"GC content:"));
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.min_gc)
                        .range(0.0..=100.0)
                        .speed(0.5),
                );
                ui.label("..");
                ui.add(
                    egui::DragValue::new(&mut dlg.opts.max_gc)
                        .range(0.0..=100.0)
                        .speed(0.5),
                );
                ui.label("%");
                ui.add_space(10.0);
                if ui.button(crate::i18n::text(ui.ctx(),"Find Primers")).clicked() {
                    dlg.run(template, circular);
                }
                if !dlg.results.is_empty() {
                    ui.label(
                    RichText::new(crate::ui_format!(ctx,"{} candidate pairs", dlg.results.len()))
                            .size(11.5)
                            .color(theme::TEXT_DIM),
                    );
                }
            });

            ui.add_space(6.0);
            if dlg.has_run && dlg.results.is_empty() {
                ui.label(
                    RichText::new(crate::i18n::text(ctx,"No primer pair satisfies these constraints."))
                        .size(11.5)
                        .color(theme::SITE_HOVER),
                );
            }
            if let Some(p) = dlg.picked_pair() {
                let n = len;
                let f_end = if p.forward_start + p.forward.len() <= n {
                    p.forward_start + p.forward.len()
                } else {
                    p.forward_start + p.forward.len() - n
                };
                let r_start = p.reverse_site_start;
                let r_end = if r_start + p.reverse.len() <= n {
                    r_start + p.reverse.len()
                } else {
                    r_start + p.reverse.len() - n
                };
                let rows = [
                    (
                        "Forward",
                        p.forward.clone(),
                        format!("{}..{}", p.forward_start + 1, f_end),
                        "→",
                        p.forward_tm,
                        p.forward_gc,
                    ),
                    (
                        "Reverse",
                        p.reverse.clone(),
                        format!("{}..{}", r_start + 1, r_end),
                        "←",
                        p.reverse_tm,
                        p.reverse_gc,
                    ),
                ];
                for (name, seq, site, arrow, tm, gc) in rows {
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(name).strong().size(11.5));
                        ui.label(RichText::new(format!("{site} {arrow}")).monospace().size(11.0));
                        ui.label(
                            RichText::new(format!("{tm:.1} °C  ·  {gc:.0}% GC"))
                                .size(11.0)
                                .color(theme::TEXT_DIM),
                        );
                    });
                    ui.label(RichText::new(seq).monospace().size(11.5));
                }
                ui.add_space(2.0);
                ui.label(
                    RichText::new(crate::ui_format!(ctx,
                        "Product:  {}..{}  ({} bp)",
                        p.product_start + 1,
                        p.product_start + p.product_len(),
                        p.product_len()
                    ))
                    .size(11.5),
                );
            }

            if dlg.results.len() > 1 {
                ui.add_space(6.0);
                ui.label(RichText::new(crate::i18n::text(ctx,"Candidates")).strong().size(12.0));
                egui::ScrollArea::vertical()
                    .max_height(140.0)
                    .auto_shrink([false, true])
                    .show(ui, |ui| {
                        for (i, p) in dlg.results.iter().enumerate() {
                            let label = format!(
                                "{} bp   F {:.1} / R {:.1} °C   ΔTm {:.1}",
                                p.product_len(),
                                p.forward_tm,
                                p.reverse_tm,
                                p.tm_delta()
                            );
                            if ui
                                .selectable_label(dlg.picked == i, RichText::new(label).size(11.5))
                                .clicked()
                            {
                                dlg.picked = i;
                            }
                        }
                    });
            }

            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let can = dlg.picked_pair().is_some();
                    if ui
                        .add_enabled(can, egui::Button::new(crate::i18n::text(ctx,"Add Both Primers")))
                        .clicked()
                    {
                        if let Some(p) = dlg.picked_pair() {
                            out = Some(DesignOutcome::Add(Box::new(p.clone())));
                        }
                    }
                    if ui.button(crate::i18n::text(ui.ctx(),"Cancel")).clicked() {
                        out = Some(DesignOutcome::Cancel);
                    }
                });
            });
        });
    if !open {
        return Some(DesignOutcome::Cancel);
    }
    out
}

/// 从文件里的引物模型转成写侧模型（Duplicate / Edit 用）。
pub fn spec_from_annotation(p: &PrimerAnnotation, template: &[u8]) -> PrimerSpec {
    PrimerSpec {
        name: p.name.clone(),
        sequence: p.sequence.clone(),
        description: p.description.clone(),
        date_added: p.date_added.clone(),
        phosphorylated: p.phosphorylated,
        sites: p
            .sites
            .iter()
            .map(|s| PrimerSiteSpec {
                range: s.range,
                strand: s.strand,
                annealed: annealed_of(p, s, template),
                tm: s.melting_temperature.unwrap_or_else(|| {
                    primer::tm(annealed_of(p, s, template).as_bytes())
                }),
            })
            .collect(),
    }
}

/// 把一对设计好的引物转成写侧模型（结合位点按模板长度回绕）。
pub fn pair_specs(
    pair: &PrimerPair,
    template_len: usize,
    forward_name: &str,
    reverse_name: &str,
) -> (PrimerSpec, PrimerSpec) {
    let wrap = |start: usize, len: usize| -> Range {
        let end = start + len;
        if end <= template_len {
            Range::new(start, end)
        } else {
            Range::new(start, end - template_len)
        }
    };
    let f = PrimerSpec {
        name: forward_name.to_string(),
        sequence: pair.forward.clone(),
        description: String::new(),
        date_added: String::new(),
        phosphorylated: false,
        sites: vec![PrimerSiteSpec {
            range: wrap(pair.forward_start, pair.forward.len()),
            strand: 0,
            annealed: pair.forward.clone(),
            tm: pair.forward_tm,
        }],
    };
    let r = PrimerSpec {
        name: reverse_name.to_string(),
        sequence: pair.reverse.clone(),
        description: String::new(),
        date_added: String::new(),
        phosphorylated: false,
        sites: vec![PrimerSiteSpec {
            range: wrap(pair.reverse_site_start, pair.reverse.len()),
            strand: 1,
            annealed: pair.reverse.clone(),
            tm: pair.reverse_tm,
        }],
    };
    (f, r)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dlg_with(seq: &str) -> PrimerDlg {
        let mut d = PrimerDlg::add("test");
        d.sequence = seq.into();
        d.last_seq = d.sequence.clone();
        d
    }

    #[test]
    fn insert_at_caret_keeps_undo_history() {
        let mut d = dlg_with("ACGT");
        d.caret = 4;
        d.insert_at_caret("TT");
        assert_eq!(d.sequence, "ACGTTT");
        d.insert_at_caret("GG");
        assert_eq!(d.sequence, "ACGTTTGG");
        d.undo_edit();
        assert_eq!(d.sequence, "ACGTTT");
        d.undo_edit();
        assert_eq!(d.sequence, "ACGT");
        d.redo_edit();
        assert_eq!(d.sequence, "ACGTTT");
    }

    #[test]
    fn reverse_complement_flips_site_and_undo_restores_it() {
        let mut d = dlg_with("AAGC");
        d.forced = Some(ForcedSite {
            range: Range::new(0, 4),
            strand: 0,
            annealed: 4,
            tail3: false,
        });
        d.reverse_complement(100);
        assert_eq!(d.sequence, "GCTT");
        let f = d.forced.as_ref().unwrap();
        assert_eq!(f.strand, 1);
        assert_eq!((f.range.start, f.range.end), (96, 100));
        // 撤销要把序列**和**位点一起还原（曾经只还原序列，界面会出现
        // 正向序列配反向位点的自相矛盾状态）。
        d.undo_edit();
        assert_eq!(d.sequence, "AAGC");
        let f = d.forced.as_ref().unwrap();
        assert_eq!(f.strand, 0);
        assert_eq!((f.range.start, f.range.end), (0, 4));
    }
}
