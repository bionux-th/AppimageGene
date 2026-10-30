//! 轻量持久化设置（零依赖）：文件夹/最近文件、显示开关、工具栏/侧栏
//! 可见性、克隆默认参数与 Golden Gate 设置。
//!
//! 存放于 `$XDG_CONFIG_HOME/helixbench/settings`（默认 `~/.config/helixbench/settings`），
//! 纯文本 `key=value`，逐行解析。读取失败一律当作空设置——不能让配置问题挡住启动。

use std::path::{Path, PathBuf};

/// 同源重组类克隆的默认重叠参数（SnapGene `[actions]` 的
/// `gibsonAssembly` / `hifiAssembly`）。
#[derive(Debug, Clone, PartialEq)]
pub struct OverlapSettings {
    pub use_exact_length: bool,
    pub exact_length: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub overlap_tm: f32,
    pub amplify_vector_without_extensions: bool,
}

impl Default for OverlapSettings {
    fn default() -> Self {
        // SnapGene 8.2.3 本机 settings.ini 的出厂值（只读核对）。
        Self {
            use_exact_length: true,
            exact_length: 40,
            min_length: 15,
            max_length: 25,
            overlap_tm: 50.0,
            amplify_vector_without_extensions: false,
        }
    }
}

impl OverlapSettings {
    /// 设计 PCR 引物时加在片段两端的同源臂长度。
    ///
    /// `useExactLength=true` 时 SG 使用 `exactOverlapLength`；否则在
    /// `minOverlapLength..=maxOverlapLength` 之间按 `overlapTm` 优化。
    /// HX 没有 Tm 计算器，非精确模式取区间中值，保证可预测。
    pub fn primer_overlap(&self) -> usize {
        if self.use_exact_length {
            self.exact_length.clamp(8, 80)
        } else {
            ((self.min_length + self.max_length) / 2).clamp(8, 80)
        }
    }

    /// 判断两个片段是否真的重叠时使用的最小长度。
    pub fn minimum_overlap(&self) -> usize {
        self.min_length.max(1)
    }
}

/// SnapGene `[actions]` 里的克隆默认参数。
#[derive(Debug, Clone, PartialEq)]
pub struct CloningSettings {
    pub pcr_primer_tm: f32,
    pub overlap_extension_tm: f32,
    pub gibson: OverlapSettings,
    pub hifi: OverlapSettings,
    /// `inFusionCloning\overlapLength`
    pub in_fusion_overlap: usize,
    /// `inFusionCloning\multipleOverlapLength`
    pub in_fusion_multiple_overlap: usize,
    /// `goldenGateCloning\upstreamBases`
    pub golden_gate_upstream_bases: usize,
    /// Golden Gate Settings 里选择的默认酶（新对话框的初始值）。
    pub golden_gate_default_enzyme: String,
    /// `topoTaCloning\vector`（商品化载体的显示名）。
    pub topo_ta_vector: String,
}

impl Default for CloningSettings {
    fn default() -> Self {
        Self {
            pcr_primer_tm: 60.0,
            overlap_extension_tm: 60.0,
            gibson: OverlapSettings::default(),
            hifi: OverlapSettings::default(),
            in_fusion_overlap: 15,
            in_fusion_multiple_overlap: 20,
            golden_gate_upstream_bases: 3,
            golden_gate_default_enzyme: "BsaI".to_string(),
            topo_ta_vector: "pCR®2.1-TOPO®".to_string(),
        }
    }
}

/// 持久化设置。
///
/// SnapGene 把显示开关和工具栏/侧栏可见性保存在 `~/.snapgene/settings.ini`
/// 的 `[dnaFiles] show\regions`、`[toolbars] topToolbar\visible`、
/// `[General] collapseContentSidebar` 等键里（本机 8.2.3 只读核对）。
/// HX 用自己的纯文本格式保存同类状态，重启后恢复。
#[derive(Debug, Clone, Default)]
pub struct Settings {
    /// 侧栏文件树当前文件夹（SnapGene 的 Files 面板）。
    pub folder: Option<PathBuf>,
    /// 最近打开的文件（新的在前）。
    pub recent: Vec<PathBuf>,
    /// 已启用的显示开关（稳定名字，见 `HelixApp::display_preferences`）。
    /// `None` = 还没有保存过偏好，保留内置默认值。
    pub display: Option<Vec<String>>,
    /// 文件侧栏展开/折叠（SG `collapseContentSidebar`）。
    pub sidebar_open: Option<bool>,
    /// 顶部工具栏可见（SG `topToolbar\visible`）。
    pub show_toolbar: Option<bool>,
    /// 画布左侧工具条可见（SG `sideToolbar\visible`）。
    pub show_side_tools: Option<bool>,
    /// 克隆默认参数（SG `[actions]`）。
    pub cloning: CloningSettings,
    /// Golden Gate Settings 窗口状态（勾选的酶、默认酶）。
    pub gg: Option<crate::ggsettings::GgSettings>,
}

#[cfg(any(windows,test))]
fn windows_config_base(roaming:Option<std::ffi::OsString>,local:Option<std::ffi::OsString>)->Option<PathBuf>{
    roaming.filter(|s|!s.is_empty()).or_else(||local.filter(|s|!s.is_empty())).map(PathBuf::from)
}

impl Settings {
    /// 配置文件路径（`$XDG_CONFIG_HOME` 优先）。
    pub fn path() -> Option<PathBuf> {
        #[cfg(target_os="windows")]
        let base=windows_config_base(std::env::var_os("APPDATA"),std::env::var_os("LOCALAPPDATA"))?;
        #[cfg(not(target_os="windows"))]
        let base = std::env::var_os("XDG_CONFIG_HOME")
            .map(PathBuf::from)
            .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")))?;
        let current=base.join("appimagegene");
        let legacy=base.join("helixbench");
        // First upgrade: move the complete old config directory. Existing AG
        // configuration takes precedence; never overwrite it on every lookup.
        if !current.exists() && legacy.is_dir() {
            if std::fs::rename(&legacy,&current).is_err() && std::fs::create_dir_all(&current).is_ok(){
                for name in ["settings","preferences","workspace"]{
                    let from=legacy.join(name);let to=current.join(name);
                    if from.is_file()&&!to.exists(){let _=std::fs::copy(from,to);}
                }
            }
        }
        Some(current.join("settings"))
    }

    /// 读取（不存在/不可读 ⇒ 默认值）。
    pub fn load() -> Self {
        let Some(p) = Self::path() else {
            return Self::default();
        };
        let Ok(text) = std::fs::read_to_string(&p) else {
            return Self::default();
        };
        Self::decode(&text)
    }

    /// 解析设置文本（纯函数，测试不需要碰真实 HOME）。
    pub fn decode(text: &str) -> Self {
        let mut s = Settings::default();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let Some((k, v)) = line.split_once('=') else {
                continue;
            };
            let v = v.trim();
            match k.trim() {
                "folder" => {
                    let p = PathBuf::from(v);
                    if p.is_dir() {
                        s.folder = Some(p);
                    }
                }
                "recent" => {
                    let p = PathBuf::from(v);
                    if !p.as_os_str().is_empty() && !s.recent.contains(&p) && s.recent.len() < 12 {
                        s.recent.push(p);
                    }
                }
                "display" => {
                    let names: Vec<String> = v
                        .split(',')
                        .map(str::trim)
                        .filter(|n| !n.is_empty())
                        .map(str::to_string)
                        .collect();
                    s.display = Some(names);
                }
                "sidebar_open" => s.sidebar_open = parse_bool(v),
                "show_toolbar" => s.show_toolbar = parse_bool(v),
                "show_side_tools" => s.show_side_tools = parse_bool(v),
                "cloning.pcr_primer_tm" => set_f32(&mut s.cloning.pcr_primer_tm, v),
                "cloning.overlap_extension_tm" => {
                    set_f32(&mut s.cloning.overlap_extension_tm, v)
                }
                "cloning.gibson.use_exact" => {
                    set_bool(&mut s.cloning.gibson.use_exact_length, v);
                }
                "cloning.gibson.exact" => set_usize(&mut s.cloning.gibson.exact_length, v),
                "cloning.gibson.min" => set_usize(&mut s.cloning.gibson.min_length, v),
                "cloning.gibson.max" => set_usize(&mut s.cloning.gibson.max_length, v),
                "cloning.gibson.tm" => set_f32(&mut s.cloning.gibson.overlap_tm, v),
                "cloning.gibson.amplify_without_extensions" => {
                    set_bool(&mut s.cloning.gibson.amplify_vector_without_extensions, v)
                }
                "cloning.hifi.use_exact" => {
                    set_bool(&mut s.cloning.hifi.use_exact_length, v);
                }
                "cloning.hifi.exact" => set_usize(&mut s.cloning.hifi.exact_length, v),
                "cloning.hifi.min" => set_usize(&mut s.cloning.hifi.min_length, v),
                "cloning.hifi.max" => set_usize(&mut s.cloning.hifi.max_length, v),
                "cloning.hifi.tm" => set_f32(&mut s.cloning.hifi.overlap_tm, v),
                "cloning.hifi.amplify_without_extensions" => {
                    set_bool(&mut s.cloning.hifi.amplify_vector_without_extensions, v)
                }
                "cloning.infusion.overlap" => {
                    set_usize(&mut s.cloning.in_fusion_overlap, v)
                }
                "cloning.infusion.multiple" => {
                    set_usize(&mut s.cloning.in_fusion_multiple_overlap, v)
                }
                "cloning.goldengate.upstream_bases" => {
                    set_usize(&mut s.cloning.golden_gate_upstream_bases, v)
                }
                "cloning.goldengate.default_enzyme" => {
                    if !v.is_empty() {
                        s.cloning.golden_gate_default_enzyme = v.to_string();
                    }
                }
                "cloning.topo_ta.vector" => {
                    if !v.is_empty() {
                        s.cloning.topo_ta_vector = v.to_string();
                    }
                }
                "gg.standard" => {
                    s.gg
                        .get_or_insert_with(crate::ggsettings::GgSettings::default)
                        .standard = parse_checked_list(v);
                }
                "gg.additional" => {
                    s.gg
                        .get_or_insert_with(crate::ggsettings::GgSettings::default)
                        .additional = parse_checked_list(v);
                }
                "gg.default" => {
                    s.gg
                        .get_or_insert_with(crate::ggsettings::GgSettings::default)
                        .set_default_enzyme_name(v);
                }
                "gg.choose_sets" => {
                    let gg = s
                        .gg
                        .get_or_insert_with(crate::ggsettings::GgSettings::default);
                    set_bool(&mut gg.choose_sets, v);
                }
                _ => {}
            }
        }
        s
    }

    /// 序列化设置文本（纯函数）。
    pub fn encode(&self) -> String {
        let mut out = String::from("# AppimageGene settings\n");
        if let Some(f) = &self.folder {
            out.push_str(&format!("folder={}\n", f.display()));
        }
        for r in &self.recent {
            out.push_str(&format!("recent={}\n", r.display()));
        }
        if let Some(names) = &self.display {
            out.push_str(&format!("display={}\n", names.join(",")));
        }
        for (key, value) in [
            ("sidebar_open", self.sidebar_open),
            ("show_toolbar", self.show_toolbar),
            ("show_side_tools", self.show_side_tools),
        ] {
            if let Some(v) = value {
                out.push_str(&format!("{key}={}\n", if v { 1 } else { 0 }));
            }
        }
        let c = &self.cloning;
        out.push_str(&format!(
            "cloning.pcr_primer_tm={}\ncloning.overlap_extension_tm={}\n",
            c.pcr_primer_tm, c.overlap_extension_tm
        ));
        for (prefix, o) in [("gibson", &c.gibson), ("hifi", &c.hifi)] {
            out.push_str(&format!(
                "cloning.{prefix}.use_exact={}\ncloning.{prefix}.exact={}\ncloning.{prefix}.min={}\ncloning.{prefix}.max={}\ncloning.{prefix}.tm={}\ncloning.{prefix}.amplify_without_extensions={}\n",
                if o.use_exact_length { 1 } else { 0 },
                o.exact_length,
                o.min_length,
                o.max_length,
                o.overlap_tm,
                if o.amplify_vector_without_extensions { 1 } else { 0 }
            ));
        }
        out.push_str(&format!(
            "cloning.infusion.overlap={}\ncloning.infusion.multiple={}\ncloning.goldengate.upstream_bases={}\ncloning.goldengate.default_enzyme={}\ncloning.topo_ta.vector={}\n",
            c.in_fusion_overlap,
            c.in_fusion_multiple_overlap,
            c.golden_gate_upstream_bases,
            c.golden_gate_default_enzyme,
            c.topo_ta_vector
        ));
        if let Some(gg) = &self.gg {
            out.push_str(&format!(
                "gg.standard={}\ngg.additional={}\ngg.default={}\ngg.choose_sets={}\n",
                encode_checked_list(&gg.standard),
                encode_checked_list(&gg.additional),
                gg.default_enzyme_name(),
                if gg.choose_sets { 1 } else { 0 }
            ));
        }
        out
    }

    /// 写回（失败静默：设置丢失不该打断使用）。
    pub fn save(&self) {
        let Some(p) = Self::path() else { return };
        if let Some(dir) = p.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(&p, self.encode());
    }
}

fn parse_bool(v: &str) -> Option<bool> {
    match v {
        "1" | "true" | "True" | "yes" => Some(true),
        "0" | "false" | "False" | "no" => Some(false),
        _ => None,
    }
}

fn set_bool(slot: &mut bool, v: &str) {
    if let Some(b) = parse_bool(v) {
        *slot = b;
    }
}

fn set_usize(slot: &mut usize, v: &str) {
    if let Ok(n) = v.parse::<usize>() {
        *slot = n;
    }
}

fn set_f32(slot: &mut f32, v: &str) {
    if let Ok(n) = v.parse::<f32>() {
        *slot = n;
    }
}

/// `Name:1,Name:0` ⇄ `Vec<(String, bool)>`。
fn parse_checked_list(v: &str) -> Vec<(String, bool)> {
    v.split(',')
        .filter_map(|item| {
            let item = item.trim();
            if item.is_empty() {
                return None;
            }
            let (name, on) = item.rsplit_once(':')?;
            Some((name.trim().to_string(), parse_bool(on.trim()).unwrap_or(false)))
        })
        .collect()
}

fn encode_checked_list(items: &[(String, bool)]) -> String {
    items
        .iter()
        .map(|(name, on)| format!("{name}:{}", if *on { 1 } else { 0 }))
        .collect::<Vec<_>>()
        .join(",")
}

/// 目录项（文件树用）。
#[derive(Debug, Clone)]
pub struct Entry {
    pub path: PathBuf,
    pub name: String,
    pub is_dir: bool,
    /// 是否是我们能打开的序列文件（原生 + 可导入的文本格式，见 [`is_openable`]）。
    pub is_seq: bool,
    /// 是否能以只读阅读视图打开（`.txt` / `.pdf`）。
    pub is_view: bool,
    /// 小写扩展名（无扩展名为空串）。
    pub ext: String,
}

/// 列出目录（目录在前、各自按名称排序）；失败返回错误文本。
pub fn list_dir(dir: &Path, limit: usize) -> Result<Vec<Entry>, String> {
    let rd = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut out: Vec<Entry> = Vec::new();
    for ent in rd.flatten() {
        if out.len() >= limit {
            out.push(Entry {
                path: PathBuf::new(),
                name: format!("… (只显示前 {limit} 项)"),
                is_dir: false,
                is_seq: false,
                is_view: false,
                ext: String::new(),
            });
            break;
        }
        let path = ent.path();
        let name = ent.file_name().to_string_lossy().to_string();
        if name.starts_with('.') {
            continue;
        }
        let is_dir = path.is_dir();
        let ext = path
            .extension()
            .map(|e| e.to_string_lossy().to_ascii_lowercase())
            .unwrap_or_default();
        // 能被打开的序列文件：SnapGene 原生 + 可导入的文本格式
        let is_seq = crate::settings::is_openable(&ext);
        let is_view = crate::viewer::kind_of_ext(&ext).is_some();
        out.push(Entry {
            path,
            name,
            is_dir,
            is_seq,
            is_view,
            ext,
        });
    }
    out.sort_by(|a, b| {
        b.is_dir
            .cmp(&a.is_dir)
            .then_with(|| a.name.to_lowercase().cmp(&b.name.to_lowercase()))
    });
    Ok(out)
}

#[cfg(test)]
mod tests {
    #[test]fn windows_uses_roaming_appdata_then_local_without_home_or_cwd(){
        assert_eq!(super::windows_config_base(Some("C:/Users/Test/AppData/Roaming".into()),Some("C:/Local".into())),Some(std::path::PathBuf::from("C:/Users/Test/AppData/Roaming")));
        assert_eq!(super::windows_config_base(Some("".into()),Some("C:/Local".into())),Some(std::path::PathBuf::from("C:/Local")));
        assert_eq!(super::windows_config_base(None,None),None);
    }
    use super::*;

    #[test]
    fn list_dir_sorts_folders_first() {
        let dir = std::env::temp_dir().join(format!("hx_settings_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("zsub")).unwrap();
        std::fs::write(dir.join("b.dna"), b"x").unwrap();
        std::fs::write(dir.join("a.txt"), b"x").unwrap();
        let list = list_dir(&dir, 100).unwrap();
        assert_eq!(list[0].name, "zsub");
        assert!(list[0].is_dir);
        assert_eq!(list[1].name, "a.txt");
        assert_eq!(list[2].name, "b.dna");
        assert!(list[2].is_seq);
        assert!(!list[1].is_seq);
        assert!(list[1].is_view);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn settings_roundtrip() {
        // 不触碰真实 HOME：只验证序列化格式本身
        let mut gg = crate::ggsettings::GgSettings::default();
        gg.standard[1].1 = false;
        gg.choose_sets = true;
        gg.set_default_enzyme_name("SapI");
        let s = Settings {
            folder: Some(PathBuf::from("/tmp/x")),
            recent: vec![PathBuf::from("/tmp/a.dna"), PathBuf::from("/tmp/b.dna")],
            display: Some(vec!["features".into(), "primers".into()]),
            sidebar_open: Some(false),
            show_toolbar: Some(true),
            show_side_tools: Some(false),
            cloning: CloningSettings {
                gibson: OverlapSettings {
                    exact_length: 35,
                    min_length: 18,
                    ..OverlapSettings::default()
                },
                in_fusion_multiple_overlap: 22,
                golden_gate_upstream_bases: 5,
                golden_gate_default_enzyme: "SapI".into(),
                topo_ta_vector: "pCR4-TOPO".into(),
                ..CloningSettings::default()
            },
            gg: Some(gg.clone()),
        };
        let text = s.encode();
        assert!(text.contains("folder=/tmp/x"));
        assert!(text.contains("/tmp/a.dna"));
        assert!(text.contains("display=features,primers"));
        assert!(text.contains("cloning.gibson.exact=35"));
        assert!(text.contains("cloning.goldengate.upstream_bases=5"));
        let back = Settings::decode(&text);
        assert_eq!(back.display, s.display);
        assert_eq!(back.sidebar_open, Some(false));
        assert_eq!(back.show_toolbar, Some(true));
        assert_eq!(back.show_side_tools, Some(false));
        assert_eq!(back.recent, s.recent);
        assert_eq!(back.cloning, s.cloning);
        assert_eq!(back.gg, Some(gg));
    }

    #[test]
    fn absent_display_preferences_are_not_mistaken_for_all_off() {
        let s = Settings::decode("folder=/tmp\n");
        assert_eq!(s.display, None);
        assert_eq!(s.sidebar_open, None);
        let s = Settings::decode("display=\n");
        assert_eq!(s.display, Some(Vec::new()));
    }
}

/// 我们能打开的分子文件扩展名 → 分子类型。
///
/// - 原生：`.dna`（DNA）/`.rna`/`.prot`；原生但只有概览页的：`.dnaln`/`.praln`/`.gel`
///   （用 [`kind_of_ext`] 之外的单列判断，见 [`is_openable`]）。
/// - 导入：FASTA/GenBank 一族（内容嗅探决定真实类型，扩展名只用于侧栏图标）。
pub fn kind_of_ext(ext: &str) -> Option<hx_core::MoleculeKind> {
    match ext {
        "dna" | "gb" | "gbk" | "genbank" | "fasta" | "fa" | "fas" | "fna" | "ffn" | "faa"
        | "seq" | "xdna" | "snapgene" => Some(hx_core::MoleculeKind::Dna),
        "rna" | "frna" => Some(hx_core::MoleculeKind::Rna),
        "prot" | "prt" | "aa" | "gp" => Some(hx_core::MoleculeKind::Protein),
        _ => None,
    }
}

/// 序列**文本**格式（FASTA / GenBank / `.fa` / `.gbk` …）。
///
/// SnapGene 的 folder view 把这类文件与原生文档分开画：原生 `.dna`/`.rna`/`.prot`
/// 是**类型色圆环 + 双螺旋**，文本格式是**亮蓝圆环 + 折角页**（oracle `sg-sb.png`）。
pub fn is_text_seq(ext: &str) -> bool {
    kind_of_ext(ext).is_some() && !matches!(ext, "dna" | "rna" | "prot")
}

/// 能否打开（含只有概览页的原生格式）。
pub fn is_openable(ext: &str) -> bool {
    kind_of_ext(ext).is_some() || matches!(ext, "dnaln" | "praln" | "gel" | "ftrs")
}

/// 能否以只读阅读视图打开（`.txt` / `.pdf`）。
pub fn is_viewable(ext: &str) -> bool {
    crate::viewer::kind_of_ext(ext).is_some()
}
