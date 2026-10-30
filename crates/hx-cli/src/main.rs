//! HelixBench 调试驱动：样例往返验证、块结构 dump、集合统计。

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use hx_formats::block::{block_type_name, write_blocks};
use hx_formats::dna::DnaFile;
use hx_formats::write::{rewrite_after_edit, Edit, Feature, NewFile};
use hx_core::{MoleculeFlags, MoleculeKind};

const DEFAULT_SAMPLE_DIR: &str =
    "/home/mky/Software/snapgene_8.2.2_linux/data/opt/gslbiotech/snapgene/resources";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("roundtrip") => cmd_roundtrip(&args[1..]),
        Some("dump") => cmd_dump(&args[1..]),
        Some("stats") => cmd_stats(&args[1..]),
        Some("new") => cmd_new(&args[1..]),
        Some("edit") => cmd_edit(&args[1..]),
        Some("export") => cmd_export(&args[1..]),
        Some("-h") | Some("--help") | None => {
            print_usage();
            ExitCode::SUCCESS
        }
        Some(other) => {
            eprintln!("未知命令：{other}\n");
            print_usage();
            ExitCode::from(2)
        }
    }
}

fn print_usage() {
    println!(
        "AppimageGene 调试驱动\n\
         \n\
         用法：\n\
         \x20 hx-cli roundtrip [目录]   SnapGene 原生文件读→写逐字节往返验证\n\
         \x20 hx-cli dump <文件>        打印块结构摘要\n\
         \x20 hx-cli stats [目录]       样例集合统计\n\
         \x20 hx-cli new --seq <序列> --out <文件> [选项]    新建 .dna/.rna/.prot\n\
         \x20 hx-cli edit <输入> --out <文件> [编辑…]        编辑并写回 .dna\n\
         \x20 hx-cli export <输入> <输出.svg> [选项]        导出图谱/序列为 SVG\n\
         \n\
         new 选项：\n\
         \x20 --from <文件>           从既有文件取序列与 flags（与 --seq 二选一）\n\
         \x20 --kind dna|rna|protein  分子类型（默认 dna；--from 时继承源文件）\n\
         \x20 --name <文本>           图谱标签（CustomMapLabel）\n\
         \x20 --desc <文本>           描述\n\
         \x20 --circular | --linear   拓扑（默认线性；蛋白无拓扑）\n\
         \x20 --feature NAME:START-END[:TYPE[:COLOR]]   添加特征（可重复）\n\
         \n\
         edit 选项（坐标 1-based 闭区间；按命令行顺序施加）：\n\
         \x20 --delete A-B            删除 A..B\n\
         \x20 --insert A:SEQ          在位置 A 前插入 SEQ\n\
         \x20 --replace A-B:SEQ       用 SEQ 替换 A..B\n\
         \x20 --primer N:SEQ:A-B:S[:p]  加引物（A-B = 1-based 结合区间，S = 0/1 链，:p = 5′ 磷酸化）\n\
         \x20 --circular | --linear   改写拓扑\n\
         \x20 --no-history          不写历史树（block 7）\n\
         \x20 --resurrectable N     父节点 resurrectable 标记（实验用）\n\
         \n\
         export 选项：\n\
         \x20 --seq                   导出序列视图（默认导出图谱）\n\
         \x20 --enzyme NAME           指定酶（可重复；默认常用集）\n\
         \x20 --all6                  改用全部「6 切点、非 IIS」酶\n\
         \x20 --no-sites | --no-labels\n\
         \x20 --width N | --height N  画布尺寸（图谱默认 900×900 或 1200×460）\n\
         \n\
         样例目录默认：{DEFAULT_SAMPLE_DIR}\n\
         可用环境变量 HELIX_SAMPLES_DIR 覆盖。"
    );
}

fn sample_dir(arg: Option<&String>) -> PathBuf {
    if let Some(d) = arg {
        return PathBuf::from(d);
    }
    if let Ok(d) = std::env::var("HELIX_SAMPLES_DIR") {
        return PathBuf::from(d);
    }
    PathBuf::from(DEFAULT_SAMPLE_DIR)
}

fn is_sample(p: &Path) -> bool {
    matches!(
        p.extension().and_then(|s| s.to_str()),
        Some("dna" | "rna" | "prot" | "dnaln" | "praln" | "gel" | "ftrs")
    )
}

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = fs::read_dir(dir) else {
        return;
    };
    for e in rd.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect(&p, out);
        } else if is_sample(&p) {
            out.push(p);
        }
    }
}

fn cmd_roundtrip(args: &[String]) -> ExitCode {
    let dir = sample_dir(args.first());
    if !dir.is_dir() {
        eprintln!("目录不存在：{}", dir.display());
        return ExitCode::from(2);
    }
    let mut files = Vec::new();
    collect(&dir, &mut files);
    files.sort();

    let mut ok = 0usize;
    let mut failed: Vec<String> = Vec::new();
    let mut bytes_total: u64 = 0;
    let mut blocks_total: u64 = 0;

    for p in &files {
        let data = match fs::read(p) {
            Ok(d) => d,
            Err(e) => {
                failed.push(format!("{}：读取失败 {e}", p.display()));
                continue;
            }
        };
        bytes_total += data.len() as u64;
        match DnaFile::parse(&data) {
            Ok(f) => {
                blocks_total += f.blocks.len() as u64;
                let out = f.to_bytes();
                if out == data {
                    ok += 1;
                } else {
                    let at = out
                        .iter()
                        .zip(data.iter())
                        .position(|(a, b)| a != b)
                        .unwrap_or(out.len().min(data.len()));
                    failed.push(format!(
                        "{}：首个差异 @{}（{} vs {} 字节）",
                        p.display(),
                        at,
                        out.len(),
                        data.len()
                    ));
                }
            }
            Err(e) => failed.push(format!("{}：解析失败 {e}", p.display())),
        }
    }

    println!(
        "往返验证：{}/{} 通过；共 {} 个块，{:.1} MiB",
        ok,
        files.len(),
        blocks_total,
        bytes_total as f64 / 1048576.0
    );
    for f in &failed {
        println!("  失败 {f}");
    }
    if failed.is_empty() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

fn cmd_dump(args: &[String]) -> ExitCode {
    let Some(path) = args.first() else {
        eprintln!("用法：hx-cli dump <文件>");
        return ExitCode::from(2);
    };
    let data = match fs::read(path) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("读取失败：{e}");
            return ExitCode::from(2);
        }
    };
    let f = match DnaFile::parse(&data) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("解析失败：{e}");
            return ExitCode::from(1);
        }
    };
    println!("文件：{path}");
    println!(
        "  头部：type_of_sequence={} ({:?})  export_version={}  import_version={}",
        f.header.type_of_sequence,
        f.header.molecule_kind(),
        f.header.export_version,
        f.header.import_version
    );
    println!("  块数：{}（不含头部）", f.blocks.len());
    for b in &f.blocks {
        let mut line = format!(
            "  @{:<9} type {:>2}  {:<32} len {:>7}",
            b.offset,
            b.block_type,
            block_type_name(b.block_type),
            b.payload.len()
        );
        if matches!(b.block_type, 0 | 21 | 32) {
            let flags = b.payload.first().copied().unwrap_or(0);
            line.push_str(&format!(
                "  flags=0x{flags:02x}  序列 {} 字符",
                b.payload.len().saturating_sub(1)
            ));
        }
        println!("{line}");
        if !b.payload.is_empty() {
            let n = b.payload.len().min(24);
            let hex: String = b.payload[..n]
                .iter()
                .map(|c| format!("{c:02x}"))
                .collect::<Vec<_>>()
                .join(" ");
            let asc: String = b.payload[..n]
                .iter()
                .map(|&c| {
                    if (32..127).contains(&c) {
                        c as char
                    } else {
                        '.'
                    }
                })
                .collect();
            println!("        {hex}  |{asc}|");
        }
    }
    ExitCode::SUCCESS
}

fn cmd_stats(args: &[String]) -> ExitCode {
    let dir = sample_dir(args.first());
    if !dir.is_dir() {
        eprintln!("目录不存在：{}", dir.display());
        return ExitCode::from(2);
    }
    let mut files = Vec::new();
    collect(&dir, &mut files);
    files.sort();

    let mut by_ext: BTreeMap<String, usize> = BTreeMap::new();
    let mut files_with_type: BTreeMap<u8, usize> = BTreeMap::new();
    let mut blocks_of_type: BTreeMap<u8, usize> = BTreeMap::new();
    let mut flag_dist: BTreeMap<u8, usize> = BTreeMap::new();

    for p in &files {
        let ext = p
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        *by_ext.entry(ext).or_default() += 1;
        let data = match fs::read(p) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("{}：读取失败 {e}", p.display());
                continue;
            }
        };
        match DnaFile::parse(&data) {
            Ok(f) => {
                let types: BTreeSet<u8> = f.blocks.iter().map(|b| b.block_type).collect();
                for t in types {
                    *files_with_type.entry(t).or_default() += 1;
                }
                for b in &f.blocks {
                    *blocks_of_type.entry(b.block_type).or_default() += 1;
                }
                if let Some(rec) = f.sequence() {
                    *flag_dist.entry(rec.flags.to_byte()).or_default() += 1;
                }
            }
            Err(e) => eprintln!("{}：解析失败 {e}", p.display()),
        }
    }

    println!("样例：{} 个文件 {by_ext:?}", files.len());
    println!("块类型统计（type / 名称 / 出现文件数 / 总块数）：");
    for (t, n_files) in &files_with_type {
        let n_blocks = blocks_of_type.get(t).copied().unwrap_or(0);
        println!(
            "  {t:>3}  {:<36} 文件 {n_files:>3}   块 {n_blocks:>4}",
            block_type_name(*t)
        );
    }
    println!("序列 flags 分布（flag → 文件数）：");
    for (f, n) in &flag_dist {
        println!("  0x{f:02x}  {n:>3}");
    }
    ExitCode::SUCCESS
}

// ---------------------------------------------------------------------------
// new / edit
// ---------------------------------------------------------------------------

/// 去掉序列里的空白与数字（允许带行号/分组的粘贴输入）。
fn clean_seq(s: &str) -> Vec<u8> {
    s.bytes()
        .filter(|b| !b.is_ascii_whitespace() && !b.is_ascii_digit())
        .collect()
}

/// 解析 `NAME:START-END[:TYPE[:COLOR]]`。
fn parse_feature(s: &str) -> Result<Feature, String> {
    let mut it = s.splitn(4, ':');
    let name = it.next().filter(|x| !x.is_empty()).ok_or("缺少名称")?;
    let range = it.next().ok_or("缺少 START-END")?;
    let (a, b) = range.split_once('-').ok_or("START-END 需形如 10-20")?;
    let start: usize = a.trim().parse().map_err(|_| "START 不是整数")?;
    let end: usize = b.trim().parse().map_err(|_| "END 不是整数")?;
    if start == 0 || end < start {
        return Err("区间需满足 1 ≤ START ≤ END".into());
    }
    let mut f = Feature::new(name, start, end);
    if let Some(t) = it.next() {
        if !t.is_empty() {
            f.kind = t.to_string();
        }
    }
    if let Some(c) = it.next() {
        if !c.is_empty() {
            f.color = c.to_string();
        }
    }
    Ok(f)
}

fn cmd_new(args: &[String]) -> ExitCode {
    let mut seq_text: Option<String> = None;
    let mut from: Option<String> = None;
    let mut name: Option<String> = None;
    let mut desc: Option<String> = None;
    let mut out: Option<String> = None;
    let mut circular = false;
    let mut kind_flag: Option<MoleculeKind> = None;
    let mut features: Vec<Feature> = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        let take = |i: &mut usize| -> Option<String> {
            *i += 1;
            args.get(*i).cloned()
        };
        match args[i].as_str() {
            "--seq" => seq_text = take(&mut i),
            "--from" => from = take(&mut i),
            "--kind" => match take(&mut i).as_deref() {
                Some("dna") | Some("DNA") => kind_flag = Some(MoleculeKind::Dna),
                Some("rna") | Some("RNA") => kind_flag = Some(MoleculeKind::Rna),
                Some("protein") | Some("prot") | Some("aa") => {
                    kind_flag = Some(MoleculeKind::Protein)
                }
                other => {
                    eprintln!("--kind 只支持 dna | rna | protein（收到 {other:?}）");
                    return ExitCode::from(2);
                }
            },
            "--name" => name = take(&mut i),
            "--desc" => desc = take(&mut i),
            "--out" | "-o" => out = take(&mut i),
            "--circular" => circular = true,
            "--linear" => circular = false,
            "--feature" => match take(&mut i).as_deref().map(parse_feature) {
                Some(Ok(f)) => features.push(f),
                Some(Err(e)) => {
                    eprintln!("--feature 解析失败：{e}");
                    return ExitCode::from(2);
                }
                None => {
                    eprintln!("--feature 需要参数 NAME:START-END[:TYPE[:COLOR]]");
                    return ExitCode::from(2);
                }
            },
            other => {
                eprintln!("未知参数：{other}（hx-cli --help 查看用法）");
                return ExitCode::from(2);
            }
        }
        i += 1;
    }

    let Some(out) = out else {
        eprintln!("缺少 --out <文件>");
        return ExitCode::from(2);
    };

    // 序列来源：--seq 或 --from（二选一）；同时定分子类型与 flags
    let (seq, kind, flags) = match (&seq_text, &from) {
        (Some(_), Some(_)) => {
            eprintln!("--seq 与 --from 只能二选一");
            return ExitCode::from(2);
        }
        (Some(s), None) => {
            let seq = clean_seq(s);
            if seq.is_empty() {
                eprintln!("序列为空");
                return ExitCode::from(2);
            }
            let kind = kind_flag.unwrap_or(MoleculeKind::Dna);
            let mut f = MoleculeFlags::from_byte(0);
            if kind != MoleculeKind::Protein {
                f = f
                    .with_topology(if circular {
                        hx_core::Topology::Circular
                    } else {
                        hx_core::Topology::Linear
                    })
                    .with_strandedness(if kind == MoleculeKind::Rna {
                        hx_core::Strandedness::Single
                    } else {
                        hx_core::Strandedness::Double
                    });
            }
            (seq, kind, f)
        }
        (None, Some(p)) => {
            let data = match fs::read(p) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("读取失败 {p}：{e}");
                    return ExitCode::from(2);
                }
            };
            let src = match DnaFile::parse(&data) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("解析失败 {p}：{e}");
                    return ExitCode::from(1);
                }
            };
            let Some(rec) = src.sequence() else {
                eprintln!("{p} 没有序列块");
                return ExitCode::from(1);
            };
            // 未显式给 --kind 时继承源文件的类型
            let kind = kind_flag.or_else(|| src.header.molecule_kind()).unwrap_or(MoleculeKind::Dna);
            let mut f = rec.flags;
            if kind == MoleculeKind::Protein {
                f = MoleculeFlags::from_byte(0);
            }
            if circular {
                f = f.with_topology(hx_core::Topology::Circular);
            }
            (rec.sequence, kind, f)
        }
        (None, None) => {
            eprintln!("需要 --seq <序列> 或 --from <文件>");
            return ExitCode::from(2);
        }
    };

    // 特征区间越界检查（环状允许跨原点：START > END）
    let len = seq.len();
    for f in &features {
        if f.end > len || f.start > len {
            eprintln!(
                "特征 {} 的区间 {}-{} 超出序列长度 {len}",
                f.name, f.start, f.end
            );
            return ExitCode::from(2);
        }
    }

    let mut nf = NewFile::of(kind, seq).with_flags(flags);
    if let Some(n) = &name {
        nf = nf.with_name(n);
    }
    if let Some(d) = &desc {
        nf = nf.with_description(d);
    }
    nf = nf.with_features(features);

    if let Err(e) = fs::write(&out, nf.to_bytes()) {
        eprintln!("写入失败 {out}：{e}");
        return ExitCode::from(1);
    }
    println!(
        "已写出 {out}（{} {}，{}，{} 个特征）",
        nf.sequence.len(),
        match kind {
            MoleculeKind::Dna | MoleculeKind::Rna => "bp",
            MoleculeKind::Protein => "aa",
        },
        if kind == MoleculeKind::Protein {
            "无链性".to_string()
        } else if flags.has(MoleculeFlags::CIRCULAR) {
            "环状".to_string()
        } else {
            "线性".to_string()
        },
        nf.features.len()
    );
    ExitCode::SUCCESS
}

fn cmd_edit(args: &[String]) -> ExitCode {
    let mut input: Option<String> = None;
    let mut out: Option<String> = None;
    let mut edits: Vec<Edit> = Vec::new();
    let mut topology: Option<bool> = None; // Some(true) = 环状
    let mut no_history = false;
    let mut resurrectable: Option<u8> = None;
    let mut primer_specs: Vec<hx_formats::primer_edit::PrimerSpec> = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        let take = |i: &mut usize| -> Option<String> {
            *i += 1;
            args.get(*i).cloned()
        };
        match args[i].as_str() {
            "--out" | "-o" => out = take(&mut i),
            "--delete" => match take(&mut i).as_deref().map(parse_range) {
                Some(Ok((a, b))) => edits.push(Edit::Delete { start: a, end: b }),
                _ => {
                    eprintln!("--delete 需要 A-B（1-based 闭区间）");
                    return ExitCode::from(2);
                }
            },
            "--insert" => match take(&mut i).as_deref().map(parse_point_seq) {
                Some(Ok((a, s))) => edits.push(Edit::Insert {
                    at: a,
                    seq: clean_seq(&s),
                }),
                _ => {
                    eprintln!("--insert 需要 A:SEQ（A 为 1-based 插入点）");
                    return ExitCode::from(2);
                }
            },
            "--replace" => match take(&mut i).as_deref().and_then(|s| {
                let (range, seq) = s.split_once(':')?;
                let (a, b) = parse_range(range).ok()?;
                Some(Edit::Replace {
                    start: a,
                    end: b,
                    seq: clean_seq(seq),
                })
            }) {
                Some(e) => edits.push(e),
                None => {
                    eprintln!("--replace 需要 A-B:SEQ（1-based 闭区间）");
                    return ExitCode::from(2);
                }
            },
            "--primer" => match take(&mut i).as_deref().and_then(parse_primer_spec) {
                Some(spec) => primer_specs.push(spec),
                None => {
                    eprintln!("--primer 需要 NAME:SEQ:A-B:S[:p]（A-B 为 1-based 结合区间）");
                    return ExitCode::from(2);
                }
            },
            "--circular" => topology = Some(true),
            "--linear" => topology = Some(false),
            "--no-history" => no_history = true,
            "--resurrectable" => match take(&mut i).and_then(|v| v.parse().ok()) {
                Some(n) => resurrectable = Some(n),
                None => {
                    eprintln!("--resurrectable 需要 1 或 2");
                    return ExitCode::from(2);
                }
            },
            other if other.starts_with("--") => {
                eprintln!("未知参数：{other}");
                return ExitCode::from(2);
            }
            path => {
                if input.is_some() {
                    eprintln!("只接受一个输入文件");
                    return ExitCode::from(2);
                }
                input = Some(path.to_string());
            }
        }
        i += 1;
    }

    let (Some(input), Some(out)) = (input, out) else {
        eprintln!("用法：hx-cli edit <输入> --out <文件> [--delete A-B] [--insert A:SEQ] [--replace A-B:SEQ]");
        return ExitCode::from(2);
    };
    let data = match fs::read(&input) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("读取失败 {input}：{e}");
            return ExitCode::from(2);
        }
    };
    let src = match DnaFile::parse(&data) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("解析失败 {input}：{e}");
            return ExitCode::from(1);
        }
    };

    let seq_type = src
        .header
        .molecule_kind()
        .map(|k| k.sequence_block_type());
    let mut file = match rewrite_after_edit(&src, &edits) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("编辑失败：{e}");
            return ExitCode::from(1);
        }
    };
    if let (Some(circ), Some(t)) = (topology, seq_type) {
        for b in &mut file.blocks {
            if b.block_type == t && !b.payload.is_empty() {
                let f = MoleculeFlags::from_byte(b.payload[0]).with_topology(if circ {
                    hx_core::Topology::Circular
                } else {
                    hx_core::Topology::Linear
                });
                b.payload[0] = f.to_byte();
            }
        }
    }

    // 引物（block 5）：与 UI 走同一条局部手术路径。
    if !primer_specs.is_empty() {
        let edits: Vec<hx_formats::primer_edit::PrimerEdit> = primer_specs
            .into_iter()
            .map(|s| hx_formats::primer_edit::PrimerEdit::Add(Box::new(s)))
            .collect();
        let at = file.blocks.iter().position(|b| b.block_type == 5);
        let cur = at
            .map(|i| file.blocks[i].payload.clone())
            .unwrap_or_else(|| hx_formats::write::empty_primers_xml().as_bytes().to_vec());
        let new = hx_formats::primer_edit::apply_primer_edits(&cur, &edits);
        match at {
            Some(i) => file.blocks[i].payload = new,
            None => {
                let pos = file
                    .blocks
                    .iter()
                    .position(|b| b.block_type == 3)
                    .map(|i| i + 1)
                    .unwrap_or(1);
                file.blocks.insert(pos, hx_formats::block::Block::new(5, new));
            }
        }
    }

    // 历史树（block 7）：根 = 本次编辑后的状态，子 = 编辑前状态。
    // 实测块序：…, 6, 7, 17, …（9,0,[3,11…],8,10,5,6,7,17,13,14,28）
    let mut history_note = String::from("无历史树");
    if !no_history && !edits.is_empty() {
        let doc_name = Path::new(&out)
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| out.clone());
        let src_name = Path::new(&input)
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| input.clone());
        match hx_formats::write::insert_history_block(
            &mut file,
            &src,
            &edits,
            &doc_name,
            &src_name,
            resurrectable,
        ) {
            Ok(nodes) => history_note = format!("历史树 {nodes} 节点（新根→{src_name}）"),
            Err(e) => {
                eprintln!("历史树构造失败（已跳过）：{e}");
            }
        }
    }

    if let Err(e) = fs::write(&out, file.to_bytes()) {
        eprintln!("写入失败 {out}：{e}");
        return ExitCode::from(1);
    }
    let old_len = src.sequence().map(|r| r.len()).unwrap_or(0);
    let new_len = file.sequence().map(|r| r.len()).unwrap_or(0);
    let dropped: Vec<String> = hx_formats::write::DROPPED_ON_SEQUENCE_CHANGE
        .iter()
        .filter(|t| src.block(**t).is_some() && file.block(**t).is_none())
        .map(|t| t.to_string())
        .collect();
    println!(
        "已写出 {out}：{old_len} → {new_len} bp，{} 次编辑，丢弃派生块 [{}]，保留 {} 块，{history_note}",
        edits.len(),
        dropped.join(","),
        file.blocks.len()
    );
    ExitCode::SUCCESS
}

/// SVG 导出：`export <输入> <输出.svg> [--seq] [--enzyme NAME]… [--all6] …`
///
/// 几何/文本由 `hx_algo::export` + `hx_render::svg` 生成（与 UI 的 SVG 按钮同一条路径）。
fn cmd_export(args: &[String]) -> ExitCode {
    let mut i = 0;
    let mut input: Option<String> = None;
    let mut out: Option<String> = None;
    let mut seq_view = false;
    let mut no_sites = false;
    let mut no_labels = false;
    let mut all6 = false;
    let mut width = 0.0f32;
    let mut height = 0.0f32;
    let mut names: Vec<String> = Vec::new();
    let take = |i: &mut usize| -> Option<String> {
        *i += 1;
        args.get(*i).cloned()
    };
    while let Some(a) = args.get(i) {
        match a.as_str() {
            "--seq" => seq_view = true,
            "--no-sites" => no_sites = true,
            "--no-labels" => no_labels = true,
            "--all6" => all6 = true,
            "--width" => width = take(&mut i).and_then(|v| v.parse().ok()).unwrap_or(0.0),
            "--height" => height = take(&mut i).and_then(|v| v.parse().ok()).unwrap_or(0.0),
            "--enzyme" => {
                if let Some(n) = take(&mut i) {
                    names.push(n);
                }
            }
            other => {
                if input.is_none() {
                    input = Some(other.to_string());
                } else if out.is_none() {
                    out = Some(other.to_string());
                }
            }
        }
        i += 1;
    }
    let (Some(input), Some(out)) = (input, out) else {
        eprintln!("用法：hx-cli export <输入> <输出.svg> [--seq] [--enzyme NAME]… [--all6]");
        return ExitCode::from(2);
    };
    let bytes = match fs::read(&input) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("读取失败 {input}：{e}");
            return ExitCode::from(1);
        }
    };
    let file = match DnaFile::parse(&bytes) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("解析失败 {input}：{e}");
            return ExitCode::from(1);
        }
    };
    let session = match hx_algo::EditSession::new(file, bytes) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("会话建立失败：{e}");
            return ExitCode::from(1);
        }
    };

    // 酶集合：显式 `--enzyme` > `--all6` > 常用集
    let enzymes: Vec<usize> = if !names.is_empty() {
        let mut v: Vec<usize> = names
            .iter()
            .filter_map(|n| {
                hx_enzymes::ENZYMES
                    .iter()
                    .position(|e| e.name.eq_ignore_ascii_case(n))
            })
            .collect();
        v.sort_unstable();
        v.dedup();
        if v.is_empty() {
            eprintln!("--enzyme 未匹配到库中酶（示例：EcoRI）");
            return ExitCode::from(2);
        }
        v
    } else if all6 {
        hx_enzymes::six_cutter_indices()
    } else {
        hx_enzymes::common_indices()
    };

    let mut style = if seq_view {
        hx_algo::export::default_sequence_style()
    } else {
        hx_algo::export::default_map_style(hx_algo::export::auto_ring(&session))
    };
    if width > 0.0 {
        style.width = width;
    }
    if height > 0.0 {
        style.height = height; // 序列视图高度按行数自动计算，此项仅对图谱生效
    }
    if no_sites {
        style.show_sites = false;
    }
    if no_labels {
        style.show_labels = false;
    }
    let title = Path::new(&input)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("untitled")
        .to_string();
    let ex = if seq_view {
        hx_algo::export::sequence_svg(&session, &enzymes, &style, &title)
    } else {
        hx_algo::export::map_svg(&session, &enzymes, &style, &title)
    };
    if let Err(e) = fs::write(&out, ex.svg.as_bytes()) {
        eprintln!("写入失败 {out}：{e}");
        return ExitCode::from(1);
    }
    let view = if seq_view {
        "序列视图"
    } else if session.circular() {
        "环状图谱"
    } else {
        "线性图谱"
    };
    println!(
        "已导出 {out}：{} bp，{view}，{} 特征带，{} 位点，{} 酶，{} 字节",
        session.len(),
        ex.bands,
        ex.sites,
        enzymes.len(),
        ex.svg.len()
    );
    ExitCode::SUCCESS
}

/// 解析 `A-B` → 0-based 半开区间（A ≥ 1，B ≥ A）。
fn parse_range(s: &str) -> Result<(usize, usize), ()> {    let (a, b) = s.split_once('-').ok_or(())?;
    let a: usize = a.trim().parse().map_err(|_| ())?;
    let b: usize = b.trim().parse().map_err(|_| ())?;
    if a == 0 || b < a {
        return Err(());
    }
    Ok((a - 1, b))
}

/// 解析 `A:SEQ` → 0-based 插入点。
fn parse_point_seq(s: &str) -> Result<(usize, String), ()> {
    let (a, seq) = s.split_once(':').ok_or(())?;
    let a: usize = a.trim().parse().map_err(|_| ())?;
    if a == 0 {
        return Err(());
    }
    Ok((a - 1, seq.to_string()))
}

// 便于将来核对：以块序列重建文件字节。
#[allow(dead_code)]
fn rebuild(blocks: &[hx_formats::Block]) -> Vec<u8> {
    let mut out = Vec::new();
    write_blocks(&mut out, blocks);
    out
}

/// 解析 `--primer NAME:SEQ:A-B:S[:p]`。
///
/// A-B 是 1-based 闭区间的**退火**区间；`SEQ` 是完整引物（可带 5′ 尾巴），
/// 退火段取 `SEQ` 的 3′ 后缀、长度 = B-A+1。
fn parse_primer_spec(s: &str) -> Option<hx_formats::primer_edit::PrimerSpec> {
    let mut it = s.splitn(5, ':');
    let name = it.next()?.trim();
    let seq: String = it.next()?.bytes().filter(|b| b.is_ascii_alphabetic()).map(|b| b.to_ascii_uppercase() as char).collect();
    let range = it.next()?;
    // parse_range 返回 (0-based 起点, 半开终点)
    let (start, end) = parse_range(range).ok()?;
    let strand: u8 = it.next()?.trim().parse().ok()?;
    let extra = it.next().unwrap_or("");
    let n = end.saturating_sub(start);
    if name.is_empty() || n == 0 || seq.len() < n {
        return None;
    }
    let annealed = seq[seq.len() - n..].to_string();
    let tm = hx_algo::primer::tm(annealed.as_bytes());
    Some(hx_formats::primer_edit::PrimerSpec {
        name: name.to_string(),
        sequence: seq,
        description: String::new(),
        date_added: String::new(),
        phosphorylated: extra.eq_ignore_ascii_case("p"),
        sites: vec![hx_formats::primer_edit::PrimerSiteSpec {
            range: hx_formats::annot::Range::new(start, end),
            strand: strand.min(1),
            annealed,
            tm,
        }],
    })
}
