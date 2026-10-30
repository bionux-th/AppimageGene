//! 系统集成：**系统文件管理器**与**系统原生文件对话框**。
//!
//! SnapGene 用 Qt 的 `QFileDialog`（原生对话框）与 `QDesktopServices`（在文件管理器中显示：
//! 英文版文案 `Show in Explorer` / `Show in Finder` / `Reveal in File Browser`）。
//! 我们保持"零外部 crate"的约定，改为调用系统上现成的工具：
//!
//! | 用途 | 使用的命令 | 说明 |
//! |---|---|---|
//! | 显示所在文件夹 | `dbus-send` → `org.freedesktop.FileManager1.ShowItems`，失败退 `xdg-open <dir>`、`gio open` | 优先**选中该文件**，退化为只打开文件夹 |
//! | 打开文件/文件夹（外部程序） | `xdg-open` / `gio open` | |
//! | 打开文件对话框 | `zenity --file-selection` / `kdialog --getopenfilename` | 探测不到就返回 `None`，调用方回退到内置路径输入小窗口 |
//! | 选择文件夹 | `zenity --file-selection --directory` / `kdialog --getexistingdirectory` | |
//! | 另存为对话框 | `zenity --file-selection --save` / `kdialog --getsavefilename` | |
//!
//! 全部通过 `std::process::Command` 运行，找不到工具时**安静地**失败（返回 `None`），
//! 界面仍可用内置路径输入完成同样操作。

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

/// 在 `PATH` 里找一个可执行文件（不引入 `which` crate）。
pub fn which(name: &str) -> Option<PathBuf> {
    #[cfg(windows)]
    {
        let exe=if name.ends_with(".exe"){name.to_owned()}else{format!("{name}.exe")};
        // Optional Poppler distribution lives beside AG, never in the current
        // working directory. Let its DLLs remain in its own bin directory.
        if matches!(name,"pdfinfo"|"pdftoppm"){
            if let Some(dir)=std::env::current_exe().ok().and_then(|p|p.parent().map(Path::to_path_buf)){
                for sub in ["poppler/bin","poppler/Library/bin","bin"]{let p=dir.join(sub).join(&exe);if p.is_file(){return Some(p);}}
            }
        }
        if let Some(path)=std::env::var_os("PATH"){for dir in std::env::split_paths(&path){let p=dir.join(&exe);if p.is_file(){return Some(p);}}}
    }
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        let cand = dir.join(name);
        if cand.is_file() {
            return Some(cand);
        }
    }
    None
}

/// Hide console windows for background helper processes on Windows.
pub fn helper_command(program:impl AsRef<std::ffi::OsStr>)->Command{
    #[allow(unused_mut)]
    let mut command=Command::new(program);
    #[cfg(windows)]{use std::os::windows::process::CommandExt;command.creation_flags(0x08000000);}
    command
}

/// 可用的原生对话框工具。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DialogTool {
    Zenity,
    Kdialog,
    /// 使用 `qarma`（zenity 的轻量替代，参数与 zenity 兼容）。
    Qarma,
}

impl DialogTool {
    fn exe(self) -> &'static str {
        match self {
            DialogTool::Zenity => "zenity",
            DialogTool::Kdialog => "kdialog",
            DialogTool::Qarma => "qarma",
        }
    }
}

/// 桌面环境是不是 KDE（`XDG_CURRENT_DESKTOP` 里含 KDE）。
///
/// KDE 上优先用 `kdialog`：选出来的对话框就是 Dolphin 那套（右侧 Places、搜索、缩略图），
/// 与 SnapGene 的 `QFileDialog` 观感一致；GTK 桌面上则优先 `zenity`。
fn kde_desktop() -> bool {
    std::env::var("XDG_CURRENT_DESKTOP")
        .map(|d| d.to_ascii_uppercase().contains("KDE"))
        .unwrap_or(false)
}

/// 探测系统上可用的原生对话框工具（结果缓存）。
pub fn dialog_tool() -> Option<DialogTool> {
    static TOOL: OnceLock<Option<DialogTool>> = OnceLock::new();
    *TOOL.get_or_init(|| {
        let order = if kde_desktop() {
            [DialogTool::Kdialog, DialogTool::Zenity, DialogTool::Qarma]
        } else {
            [DialogTool::Zenity, DialogTool::Kdialog, DialogTool::Qarma]
        };
        for t in order {
            if which(t.exe()).is_some() {
                return Some(t);
            }
        }
        None
    })
}

/// 是否存在可用的原生对话框（界面据此决定是否显示「使用系统对话框…」按钮）。
pub fn has_native_dialog() -> bool {
    #[cfg(windows)]{return true;}
    #[cfg(not(windows))]
    dialog_tool().is_some()
}

/// 原生对话框可用时的名字（状态栏提示用）。
pub fn dialog_tool_name() -> Option<&'static str> {
    #[cfg(windows)]{return Some("Windows native file dialog");}
    #[cfg(not(windows))]
    dialog_tool().map(|t| t.exe())
}

/// 运行命令并取其标准输出（首行、去掉尾随换行）；失败返回 `None`。
fn run_capture(cmd: &str, args: &[String]) -> Option<String> {
    let out = Command::new(cmd).args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
    if s.is_empty() {
        None
    } else {
        Some(s)
    }
}

/// 静默启动一个后台程序（不关心输出），返回是否成功**启动**。
fn spawn(cmd: &str, args: &[String]) -> bool {
    Command::new(cmd)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .is_ok()
}

/// 用 `xdg-open` / `gio open` 打开任意路径（文件或文件夹）。
pub fn open_path(path: &Path) -> Result<(), String> {
    #[cfg(windows)]{return Command::new("explorer.exe").arg(path).spawn().map(|_|()).map_err(|e|format!("Explorer: {e}"));}
    #[cfg(not(windows))]
    {
    let p = path.display().to_string();
    if which("xdg-open").is_some() && spawn("xdg-open", &[p.clone()]) {
        return Ok(());
    }
    if which("gio").is_some() && spawn("gio", &["open".into(), p.clone()]) {
        return Ok(());
    }
    Err(format!("no xdg-open/gio available to open {p}"))
    }
}

/// 在**系统文件管理器**里显示 `path`（文件则尽量选中它，目录则打开它）。
///
/// 对应 SnapGene 的 `Show in Finder` / `Show in Explorer` / `Reveal in File Browser`。
pub fn reveal(path: &Path) -> Result<(), String> {
    #[cfg(windows)]{
        let absolute=std::path::absolute(path).map_err(|e|e.to_string())?;
        let mut command=Command::new("explorer.exe");
        if absolute.is_dir(){command.arg(&absolute);}else{let mut argument=std::ffi::OsString::from("/select,");argument.push(&absolute);command.arg(argument);}
        return command.spawn().map(|_|()).map_err(|e|format!("Explorer: {e}"));
    }
    #[cfg(not(windows))]
    {
    let abs = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|d| d.join(path))
            .unwrap_or_else(|_| path.to_path_buf())
    };
    let dir = if abs.is_dir() {
        abs.clone()
    } else {
        abs.parent().map(|p| p.to_path_buf()).unwrap_or(abs.clone())
    };
    // 1) 有文件管理器 DBus 接口时优先：能在文件夹里**选中**该文件
    if which("dbus-send").is_some() {
        let uri = format!("file://{}", abs.display());
        let args: Vec<String> = vec![
            "--session".into(),
            "--dest=org.freedesktop.FileManager1".into(),
            "--type=method_call".into(),
            "/org/freedesktop/FileManager1".into(),
            "org.freedesktop.FileManager1.ShowItems".into(),
            format!("array:string:{uri}"),
            "string:".into(),
        ];
        if spawn("dbus-send", &args) {
            return Ok(());
        }
    }
    // 2) 退化为「打开所在文件夹」
    open_path(&dir)
    }
}

/// 弹出**系统原生**的文件打开对话框。
pub fn pick_file(title: &str, start: Option<&Path>) -> Option<PathBuf> {
    #[cfg(windows)]{return windows_dialog(title,start).pick_file();}
    #[cfg(not(windows))]
    {
    let tool = dialog_tool()?;
    let start_s = start
        .map(|p| {
            // zenity/kdialog 用目录路径 + 分隔符作为起始目录
            let mut s = p.display().to_string();
            if p.is_dir() && !s.ends_with('/') {
                s.push('/');
            }
            s
        })
        .unwrap_or_default();
    let out = match tool {
        DialogTool::Kdialog => {
            // kdialog 的顺序：先选项、再 `--getopenfilename`、最后起始目录
            // （实测 `kdialog --title T --getopenfilename <dir>` 才认这个位置参数）
            let mut args: Vec<String> = vec![
                "--title".into(),
                title.into(),
                "--getopenfilename".into(),
            ];
            if !start_s.is_empty() {
                args.push(start_s.clone());
            }
            run_capture(tool.exe(), &args)?
        }
        _ => {
            let mut args: Vec<String> = vec![
                "--file-selection".into(),
                format!("--title={title}"),
            ];
            if !start_s.is_empty() {
                args.push(format!("--filename={start_s}"));
            }
            run_capture(tool.exe(), &args)?
        }
    };
    Some(PathBuf::from(out))
    }
}

/// 弹出系统原生文件对话框并允许**多选**（SnapGene 的 `Open Files…`）。
///
/// `zenity` / `kdialog` 都用 `--multiple --separate-output` 让每行一个路径。
pub fn pick_files(title: &str, start: Option<&Path>) -> Vec<PathBuf> {
    #[cfg(windows)]{return windows_dialog(title,start).pick_files().unwrap_or_default();}
    #[cfg(not(windows))]
    {
    let Some(tool) = dialog_tool() else {
        return Vec::new();
    };
    let start_s = start
        .map(|p| {
            let mut s = p.display().to_string();
            if p.is_dir() && !s.ends_with('/') {
                s.push('/');
            }
            s
        })
        .unwrap_or_default();
    let out = match tool {
        DialogTool::Kdialog => {
            let mut args: Vec<String> = vec![
                "--title".into(),
                title.into(),
                "--multiple".into(),
                "--separate-output".into(),
                "--getopenfilename".into(),
            ];
            if !start_s.is_empty() {
                args.push(start_s.clone());
            }
            run_capture(tool.exe(), &args)
        }
        _ => {
            let mut args: Vec<String> = vec![
                "--file-selection".into(),
                "--multiple".into(),
                "--separate-output".into(),
                format!("--title={title}"),
            ];
            if !start_s.is_empty() {
                args.push(format!("--filename={start_s}"));
            }
            run_capture(tool.exe(), &args)
        }
    };
    out.map(|text| {
        text.lines()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .map(PathBuf::from)
            .collect()
    })
    .unwrap_or_default()
    }
}

/// 弹出**系统原生**的文件夹选择对话框。
pub fn pick_folder(title: &str, start: Option<&Path>) -> Option<PathBuf> {
    #[cfg(windows)]{return windows_dialog(title,start).pick_folder();}
    #[cfg(not(windows))]
    {
    let tool = dialog_tool()?;
    let start_s = start.map(|p| p.display().to_string()).unwrap_or_default();
    let out = match tool {
        DialogTool::Kdialog => {
            let mut args: Vec<String> = vec![
                "--title".into(),
                title.into(),
                "--getexistingdirectory".into(),
            ];
            if !start_s.is_empty() {
                args.push(start_s.clone());
            }
            run_capture(tool.exe(), &args)?
        }
        _ => {
            let mut args: Vec<String> = vec![
                "--file-selection".into(),
                "--directory".into(),
                format!("--title={title}"),
            ];
            if !start_s.is_empty() {
                args.push(format!("--filename={start_s}"));
            }
            run_capture(tool.exe(), &args)?
        }
    };
    Some(PathBuf::from(out))
    }
}

/// 弹出**系统原生**的另存为对话框；返回用户选择的路径。
pub fn pick_save_path(title: &str, default_path: &Path) -> Option<PathBuf> {
    #[cfg(windows)]{let mut dialog=windows_dialog(title,default_path.parent());if let Some(name)=default_path.file_name(){dialog=dialog.set_file_name(name.to_string_lossy());}return dialog.save_file();}
    #[cfg(not(windows))]
    {
    let tool = dialog_tool()?;
    let def = default_path.display().to_string();
    let out = match tool {
        DialogTool::Kdialog => run_capture(
            tool.exe(),
            &[
                "--title".into(),
                title.into(),
                "--getsavefilename".into(),
                def.clone(),
            ],
        )?,
        _ => run_capture(
            tool.exe(),
            &[
                "--file-selection".into(),
                "--save".into(),
                "--confirm-overwrite".into(),
                format!("--title={title}"),
                format!("--filename={def}"),
            ],
        )?,
    };
    Some(PathBuf::from(out))
    }
}

#[cfg(windows)]
fn windows_dialog(title:&str,start:Option<&Path>)->rfd::FileDialog{
    let mut dialog=rfd::FileDialog::new().set_title(title);
    if let Some(path)=start{if path.is_dir(){dialog=dialog.set_directory(path);}else if let Some(parent)=path.parent(){dialog=dialog.set_directory(parent);}}
    dialog
}

/// 系统文件对话框是否可用 + 用的是哪个（状态栏/提示文案）。
pub fn dialog_tool_label() -> Option<String> {
    #[cfg(windows)]{return Some("Windows native file dialog".into());}
    #[cfg(not(windows))]
    dialog_tool().map(|t| {
        if t == DialogTool::Kdialog {
            "KDE file dialog (kdialog)".to_string()
        } else {
            format!("{} file dialog", t.exe())
        }
    })
}
