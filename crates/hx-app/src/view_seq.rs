//! 序列视图（对齐 SnapGene 的版式与分块）。
//!
//! 每个「行块」（默认 60 bp）自上而下：
//!
//! ```text
//!   5' 上链碱基（10 bp 一组，组间留缝）                 行末绝对坐标
//!   ─┬─┬─┬─ 骨架线（每碱基小刻度、每 10 碱基长刻度）
//!   3' 下链碱基（反向互补）
//!   每层特征：自己的氨基酸标尺 / 氨基酸 / 特征条（无向、正向、反向或双向）
//! ```
//!
//! 顶部固定「行内相对刻度」（10/20/…，与 SnapGene 一致；绝对坐标在每行右端）。
//! 底部固定两条：**全序列特征概览**（含可视范围框，点击/拖拽定位）与
//! **缩放滑块**（可视范围滑块 + 碱基刻度 + 缩放按钮）。
//!
//! 交互：单击 = 定位光标，拖拽 = 选段，滚轮 = 滚动，`Ctrl+滚轮` = 缩放。

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Vec2};

use hx_render::map::ruler_ticks_subdiv;
use hx_render::seq::SeqView;

use crate::app::{nfmt, nfmt_plain, HelixApp};
use crate::{feat_card, theme};

/// 每多少个氨基酸标一次号（SnapGene 惯例）。
const AA_LABEL_EVERY: usize = 5;
/// 特征条基准高度（再乘缩放）。
const BAR_H: f32 = 13.0;
/// 每个翻译行的基准高度（刻度 + aa + 名称条，再乘缩放）。
const AA_ROW_H: f32 = 30.0;
/// 独立 ORF 只有彩色氨基酸带，不占特征标尺和名称条的空间。
const ORF_ROW_H: f32 = 16.0;
/// 底部「概览 + 滑块」总高度。
const BOTTOM_H: f32 = 50.0;

/// 显示用互补碱基：`.rna` 文件里存的是 DNA 字母，SnapGene 的 RNA 视图
/// 把 `T` 写成 `U`，所以 `A` 的互补显示成 `U` 而不是 `T`。
fn complement_display(kind: hx_core::molecule::MoleculeKind, b: u8) -> char {
    let c = complement(b);
    if kind == hx_core::molecule::MoleculeKind::Rna && c == 'T' {
        'U'
    } else {
        c
    }
}

/// 互补碱基（IUPAC 尽量保持语义）。
fn complement(b: u8) -> char {
    match b.to_ascii_uppercase() {
        b'A' => 'T',
        b'T' | b'U' => 'A',
        b'C' => 'G',
        b'G' => 'C',
        b'R' => 'Y',
        b'Y' => 'R',
        b'S' => 'S',
        b'W' => 'W',
        b'K' => 'M',
        b'M' => 'K',
        b'B' => 'V',
        b'V' => 'B',
        b'D' => 'H',
        b'H' => 'D',
        b'N' => 'N',
        _ => 'N',
    }
}

/// 条/块底色对应的墨色（浅底黑字、深底白字）。
fn ink_on(fill: Color32) -> Color32 {
    if theme::luma(fill) > 0.55 {
        theme::TEXT
    } else {
        Color32::WHITE
    }
}

/// 由特征片段构造**编码序列**与「编码位 → 基因组位」映射。
///
/// - 正向特征：按基因组顺序取上链；
/// - 反向特征：按基因组**逆序**取反向互补（`map[0]` = 该特征的 3′ 端基因组位）。
fn coding_map(seq: &[u8], segs: &[(usize, usize)], forward: bool) -> (Vec<usize>, Vec<u8>) {
    let mut map = Vec::new();
    let mut code = Vec::new();
    if forward {
        for (a, b) in segs {
            for bp in *a..*b {
                map.push(bp);
                code.push(seq[bp]);
            }
        }
    } else {
        for (a, b) in segs.iter().rev() {
            for bp in (*a..*b).rev() {
                map.push(bp);
                code.push(complement(seq[bp]) as u8);
            }
        }
    }
    (map, code)
}

/// 可翻译特征（CDS）：编码序列 + 基因组映射 + 翻译结果。
#[derive(Debug)]
struct Cds {
    forward: bool,
    /// 编码序列第 `j` 个碱基对应的基因组 0-based 位置（跨原点已回绕）。
    map: Vec<usize>,
    /// 每个密码子的三字母缩写，包含终止标记。
    aa: Vec<&'static str>,
    /// 标注区内的密码子数；后续连接区氨基酸不继续标号。
    annotated_codons: usize,
}

impl Cds {
    fn codon_bases(&self, k: usize) -> Option<[usize; 3]> {
        self.map.get(3 * k..3 * k + 3)?.try_into().ok()
    }

    /// 跨行/跨原点/剪接密码子拆成实际可见片段，绝不以 min..max 虚构跨度。
    fn codons_in_row(&self, row_start: usize, row_end: usize) -> Vec<CodonFragment> {
        let mut out = Vec::new();
        for k in 0..self.aa.len() {
            let Some(bases) = self.codon_bases(k) else {
                break;
            };
            let mut visible: Vec<usize> = bases
                .into_iter()
                .filter(|bp| *bp >= row_start && *bp < row_end)
                .collect();
            visible.sort_unstable();
            for bp in visible {
                if let Some(last) = out.last_mut() {
                    let last: &mut CodonFragment = last;
                    if last.codon == k && last.end == bp {
                        last.end += 1;
                        continue;
                    }
                }
                let coding_slot = bases.iter().position(|b| *b == bp).unwrap();
                out.push(CodonFragment {
                    codon: k,
                    start: bp,
                    end: bp + 1,
                    visual_slot: if self.forward {
                        coding_slot
                    } else {
                        2 - coding_slot
                    },
                });
            }
        }
        out
    }
}

#[derive(Debug)]
struct CodonFragment {
    codon: usize,
    start: usize,
    end: usize,
    /// 片段起点在完整密码子中的可视列（0/1/2）。
    visual_slot: usize,
}

/// 完整的三碱基选择；跨原点保留回绕语义，不连续的剪接密码子没有连续选区。
/// Split an unwrapped strand interval without discarding its signed cut offset.
fn duplex_selection_parts(start:i64,end:i64,len:usize)->Vec<(usize,usize)>{
    if len==0||end<=start{return Vec::new();}let n=len as i64;
    if end-start>=n{return vec![(0,len)];}let a=start.rem_euclid(n)as usize;let b=end.rem_euclid(n)as usize;
    if a<b{vec![(a,b)]}else{let mut parts=vec![(a,len)];if b>0{parts.push((0,b));}parts}
}
fn codon_selection(mut bases: [usize; 3], len: usize) -> Option<(usize, usize)> {
    bases.sort_unstable();
    if bases[2] - bases[0] == 2 {
        Some((bases[0], bases[2] + 1))
    } else if len > 2 && bases[0] == 0 && bases[2] == len - 1 {
        if bases[1] == 1 {
            Some((bases[2], 2))
        } else if bases[1] == len - 2 {
            Some((bases[1], 1))
        } else {
            None
        }
    } else {
        None
    }
}

fn codon_parts(mut bases: [usize; 3]) -> Vec<(usize, usize)> {
    bases.sort_unstable();
    let mut parts = Vec::new();
    for bp in bases {
        if let Some((_, end)) = parts.last_mut() {
            if *end == bp {
                *end += 1;
                continue;
            }
        }
        parts.push((bp, bp + 1));
    }
    parts
}

struct FeatureSegment {
    start: usize,
    end: usize,
    left_arrow: bool,
    right_arrow: bool,
}

struct Feature {
    source: usize,
    name: String,
    /// 特征类型原文（工具卡芯片用）。
    kind: String,
    /// 元素上的 `subType`（`Site` 才有）。
    sub_type: String,
    /// `/note` 纯文本（工具卡的 `NOTE` 段）。
    note: String,
    /// 是否有反向箭头（工具卡标题行右端的方向箭头）。
    reverse: bool,
    segs: Vec<FeatureSegment>,
    color: Color32,
    cds: Option<Cds>,
    orf: bool,
}

/// 终止子在 ORF/特征本体之外，显示时追加，但不延长注释范围。
fn append_stop_codon(
    seq: &[u8],
    circular: bool,
    forward: bool,
    map: &mut Vec<usize>,
    code: &mut Vec<u8>,
) {
    if seq.is_empty()
        || code.len() < 3
        || code.len() % 3 != 0
        || hx_algo::translate::codon_aa(&code[code.len() - 3..]) == "***"
    {
        return;
    }
    let Some(&last) = map.last() else { return };
    let next: Option<Vec<usize>> = (1..=3)
        .map(|d| {
            if forward {
                let bp = last + d;
                if circular {
                    Some(bp % seq.len())
                } else {
                    (bp < seq.len()).then_some(bp)
                }
            } else if circular {
                Some((last + seq.len() - d % seq.len()) % seq.len())
            } else {
                last.checked_sub(d)
            }
        })
        .collect();
    if let Some(next) = next {
        let stop: Vec<u8> = next
            .iter()
            .map(|bp| {
                if forward {
                    seq[*bp]
                } else {
                    complement(seq[*bp]) as u8
                }
            })
            .collect();
        if hx_algo::translate::codon_aa(&stop) == "***" {
            map.extend(next);
            code.extend(stop);
        }
    }
}

fn orf_translation(seq: &[u8], o: &hx_algo::orf::Orf, circular: bool) -> Cds {
    let segs = hx_render::split_wrap(o.start, o.end, seq.len());
    let forward = o.strand > 0;
    let (mut map, mut code) = coding_map(seq, &segs, forward);
    let annotated_codons = code.len() / 3;
    append_stop_codon(seq, circular, forward, &mut map, &mut code);
    Cds {
        forward,
        map,
        aa: hx_algo::translate::translate(&code),
        annotated_codons,
    }
}

fn feature_translation(
    seq: &[u8],
    f: &hx_formats::annot::FeatureAnnotation,
    circular: bool,
    annotations: &[hx_formats::annot::FeatureAnnotation],
) -> Option<Cds> {
    let segs: Vec<_> = f
        .segments
        .iter()
        .filter(|s| s.translated)
        .flat_map(|s| hx_render::split_wrap(s.range.start, s.range.end, seq.len()))
        .collect();
    let forward = f.is_forward();
    let (mut map, mut code) = coding_map(seq, &segs, forward);
    let skip = (f.reading_frame.unsigned_abs() as usize)
        .saturating_sub(1)
        .min(code.len());
    map.drain(..skip);
    code.drain(..skip);
    if code.len() < 3 {
        return None;
    }
    let annotated_codons = code.len() / 3;
    // SnapGene 在同框的蛋白标签之间显示连接区翻译与虚线，下一特征重新编号。
    // 这里只使用文件明确给出的正向 maxRunOn，不从特征名称猜测融合关系。
    if forward
        && code.len() % 3 == 0
        && hx_algo::translate::codon_aa(&code[code.len() - 3..]) != "***"
    {
        let start = map.last()? + 1;
        if let Some(mut limit) = f.max_run_on.filter(|limit| *limit > start) {
            limit = limit.min(seq.len());
            for other in annotations.iter().filter(|other| other.is_forward()) {
                if let Some(next) = other
                    .segments
                    .iter()
                    .filter(|s| s.translated)
                    .map(|s| s.range.start)
                    .filter(|a| *a >= start && (*a - start) % 3 == 0)
                    .min()
                {
                    limit = limit.min(next);
                }
            }
            for (k, codon) in seq[start..limit.max(start)].chunks_exact(3).enumerate() {
                map.extend(start + 3 * k..start + 3 * k + 3);
                code.extend_from_slice(codon);
                if hx_algo::translate::codon_aa(codon) == "***" {
                    break;
                }
            }
        }
    }
    // 特征外紧随的终止子也显示为红色 *，但不延长特征条。
    append_stop_codon(seq, circular, forward, &mut map, &mut code);
    let mut aa = hx_algo::translate::translate(&code);
    if f.translate_first_codon_as_met && aa.first() != Some(&"***") {
        aa[0] = "Met";
    }
    Some(Cds {
        forward,
        map,
        aa,
        annotated_codons,
    })
}

struct FeatureSpan {
    start: usize,
    end: usize,
    left_arrow: bool,
    right_arrow: bool,
    continues_left: bool,
    continues_right: bool,
}

struct FeatureTrack {
    feature: usize,
    level: usize,
    spans: Vec<FeatureSpan>,
    codons: Vec<CodonFragment>,
    external_label: bool,
}

#[derive(Default)]
struct FeatureLevel {
    right: f32,
    translation: bool,
    external_label: bool,
    top: f32,
    height: f32,
}

/// 同坐标的重复特征保留独立轨道；相邻区间按坐标排序后共享一层。
fn feature_tracks(
    features: &[Feature],
    start: usize,
    end: usize,
    v: &SeqView,
    label_width: impl Fn(&str) -> f32,
) -> (Vec<FeatureTrack>, Vec<FeatureLevel>) {
    let mut tracks = Vec::new();
    for (fi, f) in features.iter().enumerate() {
        if f.orf {
            continue;
        }
        let spans: Vec<_> = f
            .segs
            .iter()
            .filter_map(|s| {
                let (a, b) = (s.start.max(start), s.end.min(end));
                (a < b).then_some(FeatureSpan {
                    start: a,
                    end: b,
                    left_arrow: s.left_arrow && a == s.start,
                    right_arrow: s.right_arrow && b == s.end,
                    continues_left: a > s.start,
                    continues_right: b < s.end,
                })
            })
            .collect();
        let codons = f
            .cds
            .as_ref()
            .map(|c| c.codons_in_row(start, end))
            .unwrap_or_default();
        if spans.is_empty() && codons.is_empty() {
            continue;
        }
        let label_w = label_width(&f.name);
        let external_label = !spans.is_empty()
            && spans
                .iter()
                .all(|s| v.x_of(s.end - 1) + v.base_w - v.x_of(s.start) < label_w + 10.0);
        tracks.push(FeatureTrack {
            feature: fi,
            level: 0,
            spans,
            codons,
            external_label,
        });
    }
    let bounds = |t: &FeatureTrack| {
        let a = t
            .spans
            .iter()
            .map(|s| s.start)
            .chain(t.codons.iter().map(|c| c.start))
            .min()
            .unwrap();
        let b = t
            .spans
            .iter()
            .map(|s| s.end)
            .chain(t.codons.iter().map(|c| c.end))
            .max()
            .unwrap();
        let (mut left, mut right) = (v.x_of(a), v.x_of(b - 1) + v.base_w);
        if t.external_label {
            let w = label_width(&features[t.feature].name);
            for s in &t.spans {
                let cx = (v.x_of(s.start) + v.x_of(s.end - 1) + v.base_w) * 0.5;
                left = left.min(cx - w * 0.5);
                right = right.max(cx + w * 0.5);
            }
        }
        (left, right)
    };
    tracks.sort_by(|a, b| {
        bounds(a)
            .0
            .total_cmp(&bounds(b).0)
            .then_with(|| a.feature.cmp(&b.feature))
    });
    let mut levels: Vec<FeatureLevel> = Vec::new();
    for t in &mut tracks {
        let (left, right) = bounds(t);
        let level = levels
            .iter()
            .position(|l| l.right <= left + 0.1)
            .unwrap_or_else(|| {
                levels.push(FeatureLevel::default());
                levels.len() - 1
            });
        t.level = level;
        levels[level].right = right;
        levels[level].translation |= !t.codons.is_empty();
        levels[level].external_label |= t.external_label;
    }
    (tracks, levels)
}

struct OrfTrack {
    feature: usize,
    level: usize,
    codons: Vec<CodonFragment>,
}

/// ORF 与特征分开排版，正链在前，反链在后；重叠的阅读框各占一行。
fn orf_tracks(features: &[Feature], start: usize, end: usize) -> (Vec<OrfTrack>, usize) {
    let mut tracks = Vec::new();
    let mut level_count = 0;
    for forward in [true, false] {
        let mut strand_tracks: Vec<_> = features
            .iter()
            .enumerate()
            .filter_map(|(fi, f)| {
                let c = f.cds.as_ref()?;
                if !f.orf || c.forward != forward {
                    return None;
                }
                let codons = c.codons_in_row(start, end);
                (!codons.is_empty()).then_some(OrfTrack {
                    feature: fi,
                    level: 0,
                    codons,
                })
            })
            .collect();
        strand_tracks.sort_by_key(|t| (t.codons.iter().map(|c| c.start).min().unwrap(), t.feature));
        let mut right_edges: Vec<usize> = Vec::new();
        for mut t in strand_tracks {
            let left = t.codons.iter().map(|c| c.start).min().unwrap();
            let right = t.codons.iter().map(|c| c.end).max().unwrap();
            let level = right_edges
                .iter()
                .position(|r| *r <= left)
                .unwrap_or_else(|| {
                    right_edges.push(0);
                    right_edges.len() - 1
                });
            right_edges[level] = right;
            t.level = level_count + level;
            tracks.push(t);
        }
        level_count += right_edges.len();
    }
    (tracks, level_count)
}

fn orf_fill(forward: bool, colors: bool) -> Color32 {
    if !colors {
        Color32::from_gray(224)
    } else if forward {
        Color32::from_rgb(255, 213, 128)
    } else {
        Color32::from_rgb(187, 246, 128)
    }
}

fn feature_polygon(r: Rect, left: bool, right: bool) -> Vec<Pos2> {
    let tip = (r.width() * if left && right { 0.3 } else { 0.45 }).min(r.height() * 0.75);
    let mut points = vec![
        Pos2::new(r.left() + if left { tip } else { 0.0 }, r.top()),
        Pos2::new(r.right() - if right { tip } else { 0.0 }, r.top()),
    ];
    if right {
        points.push(Pos2::new(r.right(), r.center().y));
    }
    points.push(Pos2::new(
        r.right() - if right { tip } else { 0.0 },
        r.bottom(),
    ));
    points.push(Pos2::new(
        r.left() + if left { tip } else { 0.0 },
        r.bottom(),
    ));
    if left {
        points.push(Pos2::new(r.left(), r.center().y));
    }
    points
}

fn codon_rect(c: &CodonFragment, v: &SeqView, left: f32, aa_y: f32, z: f32) -> (Rect, f32) {
    let x0 = left + v.x_of(c.start);
    let x1 = left + v.x_of(c.end - 1) + v.base_w;
    let cx = if c.end - c.start == 3 {
        (x0 + x1) * 0.5
    } else {
        x0 + (1.5 - c.visual_slot as f32) * v.base_w
    };
    (
        Rect::from_min_max(Pos2::new(x0, aa_y - 6.0 * z), Pos2::new(x1, aa_y + 6.0 * z)),
        cx,
    )
}

struct PrimerRow {
    name: String,
    segs: Vec<(usize, usize)>,
    forward: bool,
    /// 结合位点在模板上的 0-based 半开区间（5′ 端判定用）。
    site_start: usize,
    site_end: usize,
    /// 不配对的那一段，已按**显示方向**（屏幕左→右）排好；空 = 整条引物都退火。
    tail: String,
    /// 尾巴画在结合框左边（低坐标侧）还是右边。
    ///
    /// SnapGene 里尾巴跟着它在 `sequence` 里的位置走：5′ 尾巴（sequence 前缀）
    /// 正向画左、反向画右；在序列框末尾多打的碱基（sequence 后缀，反向引物的
    /// 3′ 端）则正向画右、反向画左。
    tail_left: bool,
    /// 5′ 磷酸化（`<Primer phosphorylated="1">`）。
    phosphorylated: bool,
}

/// 一个引物结合位点在当前序列行中的可见部分。
struct PrimerSpan {
    primer: usize,
    x0: f32,
    x1: f32,
    level: usize,
    /// 这一段的末端是引物 3′ 端（画箭头）。
    tip: bool,
    /// 这一段带引物 5′ 端（画 Ⓟ）。
    five_prime: bool,
    /// 这一段外侧就是不配对的那段尾巴（画尾巴框）。
    tail_here: bool,
}

/// 一个行块的排版（可变行高：注释轨道只在有内容时占位）。
struct Row {
    start: usize,
    end: usize,
    top: f32,
    h: f32,
    tracks: Vec<FeatureTrack>,
    levels: Vec<FeatureLevel>,
    features_h: f32,
    orfs: Vec<OrfTrack>,
    orfs_h: f32,
    /// 上链（正向）引物的车道数：画在上链上方（SnapGene 实测）。
    fwd_levels: usize,
    /// 下链（反向）引物的车道数：画在特征/ORF 轨道下方。
    rev_levels: usize,
    /// 上链引物band 的高度（`row.top` 是内容起点，band 在它上面）。
    fwd_h: f32,
    primers: Vec<PrimerSpan>,
}

/// 一条引物在某一序列行里的框（退火段 / 尾巴），绘制与命中测试共用同一套几何。
struct PrimerBoxes {
    y: f32,
    x0: f32,
    x1: f32,
    /// 退火段方框。
    arect: Rect,
    /// 不配对的那段（尾巴）的方框；整条引物都退火时为 `None`。
    trect: Option<Rect>,
}

#[allow(clippy::too_many_arguments)]
fn primer_boxes(
    span: &PrimerSpan,
    primer: &PrimerRow,
    row: &Row,
    top: f32,
    content_left: f32,
    base_w: f32,
    primer_h: f32,
    bar_gap: f32,
    strand_h: f32,
) -> PrimerBoxes {
    // 上链（正向）引物画在上链上方，下链（反向）引物画在特征/ORF 轨道下方。
    let step = primer_h * 0.45;
    let after_tracks = top + 2.0 * strand_h + row.features_h + row.orfs_h;
    let y = if primer.forward {
        top - row.fwd_h
            + step
            + (row.fwd_levels as f32 - 1.0 - span.level as f32) * (primer_h + bar_gap)
    } else {
        after_tracks + 2.0 + span.level as f32 * (primer_h + bar_gap)
    };
    let x0 = content_left + span.x0;
    let x1 = content_left + span.x1;
    let arect = Rect::from_min_max(Pos2::new(x0, y), Pos2::new(x1, y + primer_h));
    let trect = if span.tail_here && !primer.tail.is_empty() {
        let tw = primer.tail.len() as f32 * base_w;
        // 尾巴框 = 尾巴那几个碱基的网格（`tw` 宽）再朝退火框方向压进 `lap`。
        // 实测（`pET28b-SprB-CDT368.dna` Primer 1 的 13 nt 尾巴，SnapGene 8.2.3，
        // 100%）：框的左边缘正好落在尾巴第一格的格线上（613），右边缘只压进
        // 退火框（760）约 3 px —— 是 0.4 × step，不是一整个 step；框里的字母
        // 仍在模板网格上（见绘制处）。
        let lap = step * 0.4;
        let tx0 = if primer.tail_left { x0 - tw } else { x1 - lap };
        let ty = if primer.forward { y - step } else { y + step };
        Some(Rect::from_min_max(
            Pos2::new(tx0, ty),
            Pos2::new(tx0 + tw + lap, ty + primer_h),
        ))
    } else {
        None
    };
    PrimerBoxes { y, x0, x1, arect, trect }
}


/// 尾巴文字第一格的左边缘（字符中心 = `left + (i + 0.5) * base_w`）。
///
/// 尾巴碱基不在模板里，但 SnapGene 仍把它们排在**模板网格**上：左尾从退火框
/// 左侧向左数 `n` 格（与上方模板字母逐列对齐），右尾从退火框右边缘向右排。
fn tail_text_left(x0: f32, x1: f32, n: usize, base_w: f32, tail_left: bool) -> f32 {
    if tail_left {
        x0 - n as f32 * base_w
    } else {
        x1
    }
}

/// 退火框 + 尾巴框的**并集轮廓**（两框只在角上重叠，SnapGene 只描外轮廓，
/// 重叠区里的两条内边不画 —— 否则会有一条竖线压在末位碱基上）。
fn primer_union_outline(a: Rect, t: Rect, tail_left: bool, forward: bool) -> Vec<Pos2> {
    let (ax0, ay0, ax1, ay1) = (a.left(), a.top(), a.right(), a.bottom());
    let (tx0, ty0, tx1, ty1) = (t.left(), t.top(), t.right(), t.bottom());
    match (tail_left, forward) {
        (true, true) => vec![
            Pos2::new(tx0, ty0),
            Pos2::new(tx1, ty0),
            Pos2::new(tx1, ay0),
            Pos2::new(ax1, ay0),
            Pos2::new(ax1, ay1),
            Pos2::new(ax0, ay1),
            Pos2::new(ax0, ty1),
            Pos2::new(tx0, ty1),
        ],
        (false, true) => vec![
            Pos2::new(ax0, ay0),
            Pos2::new(tx0, ay0),
            Pos2::new(tx0, ty0),
            Pos2::new(tx1, ty0),
            Pos2::new(tx1, ty1),
            Pos2::new(ax1, ty1),
            Pos2::new(ax1, ay1),
            Pos2::new(ax0, ay1),
        ],
        (true, false) => vec![
            Pos2::new(tx0, ty0),
            Pos2::new(ax0, ty0),
            Pos2::new(ax0, ay0),
            Pos2::new(ax1, ay0),
            Pos2::new(ax1, ay1),
            Pos2::new(tx1, ay1),
            Pos2::new(tx1, ty1),
            Pos2::new(tx0, ty1),
        ],
        (false, false) => vec![
            Pos2::new(ax0, ay0),
            Pos2::new(ax1, ay0),
            Pos2::new(ax1, ty0),
            Pos2::new(tx1, ty0),
            Pos2::new(tx1, ty1),
            Pos2::new(tx0, ty1),
            Pos2::new(tx0, ay1),
            Pos2::new(ax0, ay1),
        ],
    }
}

pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) -> egui::Rect {
    // SnapGene 的序列视图没有顶部工具条：编辑走 Edit 菜单 / 右键菜单 / 快捷键。
    // 这里只有一条**选区信息条**（SnapGene 的 SelectionBar），画在序列区底部。
    let _ = &app;
    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    if rect.width() < 160.0 || rect.height() < 120.0 {
        return rect;
    }
    let Some(doc) = app.doc.as_ref() else {
        return rect;
    };
    let len = doc.len();
    if len == 0 {
        return rect;
    }
    let seq: Vec<u8> = doc.sequence().to_vec();
    let kind = doc.kind();
    let sequence_circular=doc.circular();

    let ends = doc.end_properties();
    let stored = doc.session.stored_sequence();
    let top_range = ends.top_range(stored.len());
    let left_bottom_tail = stored[..top_range.start].to_vec();
    let right_bottom_tail = stored[top_range.end..].to_vec();
    let bottom = ends.bottom_range(stored.len());
    let bottom_start = bottom.start.saturating_sub(top_range.start);
    let bottom_end = bottom.end.saturating_sub(top_range.start).min(len);

    // ---- 数据准备：特征 / 酶位点 / 可翻译 CDS ----
    let mut feats = Vec::new();
    if app.show_features {
        for (fi, f) in doc.features().iter().enumerate() {
            if app.hidden_features.contains(&fi) {
                continue;
            }
            let color = theme::hex(
                f.segments.first().map(|s| s.color.as_str()).unwrap_or(""),
                theme::FEATURE_FALLBACK,
            );
            let segs = f
                .segments
                .iter()
                .flat_map(|s| {
                    hx_render::split_wrap(s.range.start, s.range.end, len)
                        .into_iter()
                        .map(|(a, b)| FeatureSegment {
                            start: a,
                            end: b,
                            left_arrow: f.has_reverse_arrow() && a == s.range.start,
                            right_arrow: f.has_forward_arrow() && b == s.range.end,
                        })
                })
                .collect();
            feats.push(Feature {
                source: fi,
                name: f.name.clone(),
                kind: f.kind.clone(),
                sub_type: f.sub_type.clone(),
                note: f.note.clone(),
                reverse: f.has_reverse_arrow(),
                segs,
                color,
                cds: feature_translation(&seq, f, doc.circular(), doc.features()),
                orf: false,
            });
        }
    }
    // ORF 独立于特征显示：与 CDS 同框/重叠的 ORF 也保留。
    if app.show_translation {
        for (oi, o) in app.displayed_orfs(&seq, doc.circular())
            .iter()
            .enumerate()
        {
            let frame = if o.strand > 0 { o.start % 3 } else { (len.wrapping_sub(o.end)) % 3 };
            if app.translation_frames & (1 << (frame + if o.strand < 0 { 3 } else { 0 })) == 0 { continue; }
            let segs: Vec<_> = hx_render::split_wrap(o.start, o.end, len)
                .into_iter()
                .map(|(a, b)| FeatureSegment {
                    start: a,
                    end: b,
                    left_arrow: o.strand < 0 && a == o.start,
                    right_arrow: o.strand > 0 && b == o.end,
                })
                .collect();
            let cds = orf_translation(&seq, o, doc.circular());
            feats.push(Feature {
                source: doc.features().len() + oi,
                name: format!("ORF ({} aa)", cds.annotated_codons),
                kind: "ORF".into(),
                sub_type: String::new(),
                note: String::new(),
                reverse: o.strand < 0,
                segs,
                color: theme::FEATURE_FALLBACK,
                cds: Some(cds),
                orf: true,
            });
        }
    }
    let primers: Vec<PrimerRow> = if app.show_primers {
        let template = doc.sequence();
        let params=hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap_or_default());
        doc.primers()
            .iter()
            .flat_map(|p| {
                let primary=p.sites.iter().enumerate().max_by(|(_,a),(_,b)|a.melting_temperature.unwrap_or(f64::NEG_INFINITY).total_cmp(&b.melting_temperature.unwrap_or(f64::NEG_INFINITY))).map(|(i,_)|i).unwrap_or(0);
                p.sites.iter().enumerate().filter_map(move |(site_i,s)| {
                    let tm=s.melting_temperature.unwrap_or_else(||{
                        let n=s.annealed_bases.min(p.sequence.len());
                        hx_algo::primer::tm(p.sequence[p.sequence.len().saturating_sub(n)..].as_bytes())
                    });
                    if site_i!=primary && (!params.show_additional_five_prime_matches || tm<params.min_melting_temperature) {return None;}
                    let n = p.sequence.len();
                    let annealed = s.annealed_bases.min(n);
                    // 尾巴在 sequence 的哪一端由模板决定（见 `primer_tail`）。
                    let five_prime_tail =
                        crate::primer_tail::annealed_is_suffix(&p.sequence, s, template);
                    if site_i != primary && annealed < if five_prime_tail { params.min_continuous_match_len } else { params.minimum_five_prime_annealing } { return None; }
                    let tail_src = if five_prime_tail {
                        p.sequence.get(..n - annealed).unwrap_or("")
                    } else {
                        p.sequence.get(annealed..).unwrap_or("")
                    };
                    let forward = s.strand != 1;
                    let tail: String = if forward {
                        tail_src.to_string()
                    } else {
                        tail_src.chars().rev().collect()
                    };
                    Some(PrimerRow {
                        name: p.name.clone(),
                        segs: hx_render::split_wrap(s.range.start, s.range.end, len),
                        forward,
                        site_start: s.range.start,
                        site_end: s.range.end,
                        tail,
                        tail_left: forward == five_prime_tail,
                        phosphorylated: p.phosphorylated,
                    })
                })
            })
            .filter(|p| !p.segs.is_empty())
            .collect()
    } else {
        Vec::new()
    };
    let sites: Vec<(String, usize, bool)> = if app.show_sites {
        let mut v = Vec::new();
        for (k, &eidx) in app.sites.selected.iter().enumerate() {
            if app.hidden_enzymes.contains(&k) {
                continue;
            }
            let Some(e) = hx_enzymes::ENZYMES.get(eidx) else {
                continue;
            };
            let unique = app.sites.hits.get(k).map(|h| h.len() == 1).unwrap_or(false);
            for &pos in app.sites.hits.get(k).map(|h| h.as_slice()).unwrap_or(&[]) {
                v.push((e.name.to_string(), pos as usize, unique));
            }
        }
        // Stack identical recognition/cut geometry on one shared leader.
        let mut groups:std::collections::BTreeMap<(usize,String,Vec<(Option<i64>,Option<i64>)>),(String,bool)>=std::collections::BTreeMap::new();
        for(name,start,unique)in v{let index=hx_enzymes::catalog::CATALOG.iter().position(|e|e.name==name).unwrap();let g=crate::enzyme_labels::geometry(index,start,&seq,doc.circular());let entry=groups.entry((start,g.pattern,g.cuts)).or_insert((String::new(),true));if !entry.0.is_empty(){entry.0.push('\n');}entry.0.push_str(&name);entry.1&=unique;}
        groups.into_iter().map(|((start,_,_),(name,unique))|(name,start,unique)).collect()
    } else {
        Vec::new()
    };

    // ---- 版式参数 ----
    let z = app.seq_zoom.clamp(0.45, 2.6);
    let mut v = app.seq;
    // Keep the pitch close to SnapGene's sequence grid and reserve enough
    // room on the right for the absolute row-end coordinate.
    let base_w = 10.8 * z;
    let group_gap = 0.0;
    let strand_h = 30.0 * z;
    let bar_h = BAR_H * z;
    let bar_gap = 3.0 * z;
    // 引物框高度按真 SnapGene 实测：框高 ≈ 1.8 个碱基宽，名字单独占一条。
    let primer_h = 16.0 * z;
    let primer_label_h = 13.0 * z;
    let primer_pitch = primer_h + primer_label_h + 2.0 * z;
    let aa_h = AA_ROW_H * z;
    let block_gap = 18.0 * z;
    let pad_left = (70.0 * z).max(66.0) + left_bottom_tail.len() as f32 * base_w;
    let pad_right = 94.0 * z + right_bottom_tail.len() as f32 * base_w;
    // 留白不随缩放缩小；绘制、裁剪、自动换行和命中测试共用这一内容区域。
    let side_margin = 24.0;
    // 顶部给两组悬浮工具条留出完整安全区。视图页签已从停靠栏移入画布，
    // 因而这里补回原页签高度，并把标尺/第一行序列整体下移，避免右上控件遮住行末数字。
    let ruler_h = 62.0 + 22.0 * z;
    v.group_bp = 10;
    v.base_w = base_w;
    v.group_gap = group_gap;
    v.line_h = 1.0;
    let avail = (rect.width() - 2.0 * side_margin - pad_left - pad_right).max(80.0);
    let stride = 10.0 * base_w + group_gap;
    let groups = (((avail + group_gap) / stride).floor() as usize).clamp(2, 60);
    // 行宽：None = 换行到窗口（自适应，10 的整数倍）；Some(n) = 固定每行 n 碱基
    v.line_bp = match app.seq_line_bp {
        Some(n) => n.clamp(10, 2000),
        None => (groups * 10).max(10),
    };
    // 固定线宽可能超宽：内容宽度与可平移量
    let per_line_groups = v.line_bp.div_ceil(v.group_bp.max(1)) as f32;
    let content_w = per_line_groups * stride - group_gap;
    let pan_max = (content_w - avail).max(0.0);
    app.seq_pan_x = app.seq_pan_x.clamp(0.0, pan_max);
    v.pad_left = pad_left - app.seq_pan_x;
    v.pad_top = 0.0;
    let content_left = rect.left() + side_margin;
    let content_top = rect.top() + ruler_h;
    let content_right = rect.right() - side_margin;
    // 选区信息条（SnapGene 的 SelectionBar）：内容区底部一条
    let sel_bar_h = 18.0;
    let content_bottom = rect.bottom() - BOTTOM_H - sel_bar_h;
    let view_h = (content_bottom - content_top).max(40.0);

    let small_font = FontId::proportional((11.0 * z).max(8.5));
    // ---- 每层统一安放氨基酸与自己的特征条 ----
    let n_rows = v.row_count(len);
    let mut rows: Vec<Row> = Vec::with_capacity(n_rows);
    let mut enzyme_labels=Vec::with_capacity(n_rows);
    // 顶上留一点空隙：第一行的正向引物名字就画在内容起点上，不留缝会顶到标尺。
    let mut y = 6.0 * z;
    for row in 0..n_rows {
        let start = row * v.line_bp;
        let end = (start + v.line_bp).min(len);
        let (tracks, mut levels) = feature_tracks(&feats, start, end, &v, |name| {
            painter
                .layout_no_wrap(name.to_owned(), small_font.clone(), theme::TEXT)
                .size()
                .x
        });
        let mut features_h = 0.0;
        for level in &mut levels {
            level.top = features_h;
            level.height = bar_h
                + bar_gap
                + if level.translation { aa_h } else { 0.0 }
                + if level.external_label { 12.0 * z } else { 0.0 };
            features_h += level.height;
        }
        let (orfs, orf_levels) = orf_tracks(&feats, start, end);
        let orfs_h = if orf_levels == 0 {
            0.0
        } else {
            6.0 * z + orf_levels as f32 * ORF_ROW_H * z
        };
        let mut fwd_right: Vec<f32> = Vec::new();
        let mut rev_right: Vec<f32> = Vec::new();
        let mut primer_spans = Vec::new();
        for (pi, primer) in primers.iter().enumerate() {
            let right = if primer.forward {
                &mut fwd_right
            } else {
                &mut rev_right
            };
            for (si, &(s, e)) in primer.segs.iter().enumerate() {
                let a = s.max(start);
                let b = e.min(end);
                if a >= b {
                    continue;
                }
                let x0 = v.x_of(a);
                let x1 = v.x_of(b - 1) + base_w;
                let level = match right.iter().position(|r| *r + 5.0 * z <= x0) {
                    Some(level) => {
                        right[level] = x1;
                        level
                    }
                    None => {
                        right.push(x1);
                        right.len() - 1
                    }
                };
                primer_spans.push(PrimerSpan {
                    primer: pi,
                    x0,
                    x1,
                    level,
                    tip: if primer.forward {
                        si + 1 == primer.segs.len() && b == e
                    } else {
                        si == 0 && a == s
                    },
                    // 5′ 端：正向在位点左端、反向在位点右端（与真 SnapGene 一致）
                    five_prime: if primer.forward {
                        s == primer.site_start && a == s
                    } else {
                        e == primer.site_end && b == e
                    },
                    // 尾巴挂在结合框的哪一头由 `tail_left` 决定。
                    tail_here: if primer.tail_left {
                        s == primer.site_start && a == s
                    } else {
                        si + 1 == primer.segs.len() && b == e
                    },
                });
            }
        }
        let items=sites.iter().enumerate().filter(|(_,(_,p,_))|*p>=start&&*p<end).map(|(i,(name,p,_))|{
            let width=name.lines().map(|n|painter.layout_no_wrap(n.to_owned(),FontId::proportional((12.0*z).max(10.0)),theme::TEXT).size().x).fold(0.0,f32::max);
            (i,v.x_of(*p)+base_w*0.5,width)
        }).collect::<Vec<_>>();
        let(labels,lanes)=crate::enzyme_labels::layout(&items,v.x_of(start),v.x_of(end.saturating_sub(1))+base_w);
        let group_rows=*crate::enzyme_labels::lane_offsets(&labels,&sites).last().unwrap_or(&0);
        enzyme_labels.push(labels);
        let enzyme_h=if lanes==0{0.0}else{group_rows as f32*16.0*z+16.0*z};
        // 上链引物画在行首的上方（SnapGene：正向引物在上链之上），所以行内容起点
        // 从 band 下面算；行高里仍要含这段 band。
        // 引物框的 5′ 尾巴要往链外侧再错开半个框高，band 得带上这段。
        let primer_step = primer_h * 0.3;
        let band = |n: usize| {
            if n == 0 {
                0.0
            } else {
                (n - 1) as f32 * primer_pitch + primer_h + primer_label_h + primer_step
            }
        };
        let fwd_h = band(fwd_right.len());
        let rev_h = band(rev_right.len());
        let h = fwd_h
            + 2.0 * strand_h
            + features_h
            + orfs_h
            + rev_h
            + enzyme_h
            + block_gap;
        rows.push(Row {
            start,
            end,
            top: y + fwd_h + enzyme_h,
            h,
            tracks,
            levels,
            features_h,
            orfs,
            orfs_h,
            fwd_levels: fwd_right.len(),
            rev_levels: rev_right.len(),
            fwd_h,
            primers: primer_spans,
        });
        y += h;
    }
    let total_h = y;

    // ---- 滚动 / 缩放 ----
    let max_scroll = (total_h + ruler_h - view_h).max(0.0);
    if resp.hovered() {
        let (scroll, ctrl) = ui.input(|i| {
            (
                i.smooth_scroll_delta.y,
                i.modifiers.ctrl || i.modifiers.command,
            )
        });
        if ctrl {
            app.seq_zoom = (app.seq_zoom * (1.0 + scroll * 0.0015)).clamp(0.45, 2.6);
        } else if scroll.abs() > 0.1 {
            v.scroll_y = (v.scroll_y - scroll).clamp(0.0, max_scroll);
        }
    }
    v.scroll_y = v.scroll_y.clamp(0.0, max_scroll);

    // 行 → 屏幕 y（内容坐标 + 滚动）
    let row_y = |r: &Row| content_top + r.top - v.scroll_y;
    let row_at = |py: f32| -> Option<usize> {
        let cy = py - content_top + v.scroll_y;
        if cy < 0.0 || cy > total_h {
            return None;
        }
        // 行高不一，二分定位
        let mut lo = 0usize;
        let mut hi = rows.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if rows[mid].top + rows[mid].h <= cy {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        (lo < rows.len()).then_some(lo)
    };
    // x → bp（组间缝隙/行末也吸附到最近碱基）
    let bp_from_x = |row: usize, px: f32| -> Option<usize> {
        let line_start = row * v.line_bp;
        if line_start >= len {
            return None;
        }
        let line_end = (line_start + v.line_bp).min(len);
        let local_x = px - content_left - v.pad_left;
        let stride_g = v.group_bp as f32 * v.base_w + v.group_gap;
        let mut g = (local_x.max(0.0) / stride_g.max(1.0)).floor() as usize;
        let per_line = v.line_bp / v.group_bp.max(1);
        g = g.min(per_line.saturating_sub(1));
        let within = ((local_x - g as f32 * stride_g).max(0.0) / v.base_w).floor() as usize;
        let bp =
            (line_start + g * v.group_bp + within.min(v.group_bp - 1)).min(len.saturating_sub(1));
        (bp < line_end).then_some(bp)
    };
    let bp_at = |p: Pos2| -> Option<usize> {
        if p.x < content_left + pad_left - 6.0 || p.x > content_right {
            return None;
        }
        // 顶部刻度与底部「概览/滑块」区不参与选段
        if p.y < content_top || p.y > content_bottom {
            return None;
        }
        let row = row_at(p.y)?;
        bp_from_x(row, p.x)
    };

    // 命中区直接使用绘制时同一组密码子片段；反向、跨行、重复特征不会混淆。
    let mut residue_hits = Vec::new();
    for row in &rows {
        let after_primers = row_y(row) + 2.0 * strand_h;
        let feature_rows = row.tracks.iter().map(|t| {
            (
                t.feature,
                t.codons.as_slice(),
                after_primers + row.levels[t.level].top + 22.0 * z,
            )
        });
        let orf_rows = row.orfs.iter().map(|t| {
            (
                t.feature,
                t.codons.as_slice(),
                after_primers + row.features_h + 6.0 * z + (t.level as f32 + 0.5) * ORF_ROW_H * z,
            )
        });
        for (fi, codons, aa_y) in feature_rows.chain(orf_rows) {
            let f = &feats[fi];
            let Some(cds) = &f.cds else { continue };
            for fragment in codons {
                let (hit, _) = codon_rect(fragment, &v, content_left, aa_y, z);
                let hit = hit.intersect(Rect::from_min_max(
                    Pos2::new(content_left, content_top),
                    Pos2::new(content_right, content_bottom),
                ));
                residue_hits.push((
                    hit,
                    f.source,
                    fragment.codon,
                    cds.codon_bases(fragment.codon).unwrap(),
                    cds.aa[fragment.codon],
                ));
            }
        }
    }
    // 特征条与特征名的命中区（点条或点名字都选中该特征；用户 2026-10-18：
    // 「点击标签不会产生高亮」—— SnapGene 里点特征名等于选中该特征，图谱与序列视图都要有）。
    let mut feat_hits: Vec<(Rect, usize, bool)> = Vec::new();
    for row in &rows {
        let after_primers = row_y(row) + 2.0 * strand_h;
        for track in &row.tracks {
            let f = &feats[track.feature];
            if f.orf {
                continue;
            }
            let level = &row.levels[track.level];
            let bar_y = after_primers + level.top + if level.translation { aa_h } else { 0.0 };
            let name_size = painter
                .layout_no_wrap(f.name.clone(), small_font.clone(), theme::TEXT)
                .size();
            for span in &track.spans {
                let x0 = content_left + v.x_of(span.start);
                let x1 = content_left + v.x_of(span.end - 1) + base_w;
                let bar = Rect::from_min_max(Pos2::new(x0, bar_y), Pos2::new(x1, bar_y + bar_h));
                feat_hits.push((bar, f.source, false));
                if name_size.x + 10.0 * z > bar.width() && track.external_label {
                    // 条里放不下、写到条下方的名字：给它自己一块命中区
                    feat_hits.push((
                        Rect::from_min_size(
                            Pos2::new(bar.center().x - name_size.x * 0.5, bar.bottom() + 2.0 * z),
                            name_size,
                        ),
                        f.source,
                        true,
                    ));
                }
            }
        }
    }
    // 引物框（退火段 / 尾巴 / 引物名）的命中区：单击选中、双击进 `Edit Primer…`。
    let mut primer_hits: Vec<(Rect, usize)> = Vec::new();
    for row in &rows {
        let top = row_y(row);
        for span in &row.primers {
            let primer = &primers[span.primer];
            let g = primer_boxes(
                span, primer, row, top, content_left, base_w, primer_h, bar_gap, strand_h,
            );
            primer_hits.push((g.arect, span.primer));
            if let Some(t) = g.trect {
                primer_hits.push((t, span.primer));
            }
            let size = painter
                .layout_no_wrap(primer.name.clone(), small_font.clone(), theme::TEXT)
                .size();
            let cy = if primer.forward {
                g.y - 2.0 * z
            } else {
                g.y + primer_h + 1.5 * z
            };
            let ny = if primer.forward { cy - size.y } else { cy };
            primer_hits.push((
                Rect::from_min_size(
                    Pos2::new((g.x0 + g.x1) * 0.5 - size.x * 0.5, ny),
                    size,
                ),
                span.primer,
            ));
        }
    }
    // ---- 交互：光标 / 拖选 / 氨基酸三碱基选择 ----
    if resp.double_clicked() {
        // 双击引物框 / 引物名：选中并打开 `Edit Primer…`（SnapGene 实测 `t16`）。
        if let Some(p) = resp.interact_pointer_pos() {
            if let Some((_, pi)) = primer_hits.iter().find(|(r, _)| r.contains(p)) {
                app.selected_primer = Some(*pi);
                app.selected_feature = None;
                app.selected_residue = None;
                app.open_edit_primer_dialog();
            }
        }
    }
    if resp.clicked() {
        // 单击序列 = 把键盘从侧栏搜索框等文本控件手里拿回来，否则输入/快捷键会被
        // `egui_wants_keyboard_input()` 吞掉（「触发编辑不稳定」的根因）。
        ui.ctx().memory_mut(|m| m.stop_text_input());
        if let Some(p) = resp.interact_pointer_pos() {
            if let Some((_, pi)) = primer_hits.iter().find(|(r, _)| r.contains(p)) {
                // 单击引物 = 选中（状态栏照 SnapGene 的 `Selected: …` 文案）。
                if let Some(pr) = app.doc.as_ref().and_then(|d| d.primers().get(*pi)) {
                    let (a, b) = pr
                        .sites
                        .first()
                        .map(|s| (s.range.start + 1, s.range.end))
                        .unwrap_or((0, 0));
                    let ann = pr.sites.first().map(|s| s.annealed_bases).unwrap_or(0);
                    app.status = format!(
                        "Selected: {} ({} .. {} = {} annealed bases)",
                        pr.name, a, b, ann
                    );
                }
                app.selected_primer = Some(*pi);
                app.selected_feature = None;
                app.selected_residue = None;
                app.selection = None;
                app.cursor_bp = None;
                app.sel_anchor = None;
            } else if let Some((_, source, k, bases, aa)) =
                residue_hits.iter().find(|(r, ..)| r.contains(p))
            {
                app.selected_residue = Some((*source, *k));
                app.selection = codon_selection(*bases, len);
                app.cursor_bp = None;
                app.sel_anchor = None;
                app.status = format!(
                    "{} · {} aa {} · 3 bp",
                    feats
                        .iter()
                        .find(|f| f.source == *source)
                        .map(|f| f.name.as_str())
                        .unwrap_or(""),
                    aa,
                    k + 1
                );
            } else if let Some((_, source, _)) = feat_hits
                .iter()
                .filter(|(r, _, is_label)| *is_label && r.expand(1.0).contains(p))
                .next()
                .or_else(|| feat_hits.iter().find(|(r, _, _)| r.contains(p)))
            {
                // 点特征条 / 特征名 = 选中该特征（与图谱、单链序列视图一致）
                app.selected_feature = Some(*source);
                app.selected_residue = None;
                app.selection = None;
                app.cursor_bp = None;
                app.sel_anchor = None;
                app.selected_primer = None;
            } else if !p.y.is_nan() && p.y >= content_top && p.y <= content_bottom {
                app.cursor_bp = bp_at(p);
                app.selection = None;
                app.selected_residue = None;
                app.selected_feature = None;
                app.selected_primer = None;
                app.sel_anchor = None;
            }
        }
    }
    if resp.drag_started() {
        app.selected_residue = None;
        app.selected_feature = None;
        app.selected_primer = None;
        app.sel_anchor = resp.interact_pointer_pos().and_then(bp_at);
    }
    if resp.dragged() {
        if let (Some(a), Some(p)) = (app.sel_anchor, resp.interact_pointer_pos()) {
            if let Some(b) = bp_at(p) {
                let (s, e) = if a <= b { (a, b + 1) } else { (b, a + 1) };
                app.selection = Some((s, e));
                app.cursor_bp = Some(b);
            }
        }
    }
    if resp.drag_stopped() {
        app.sel_anchor = None;
    }
    // 悬停的特征（条或条名字 → 该特征）：oracle `dn7-chip.png` 实测 SnapGene 悬停时
    // 条描边变红 `(253,1,1)`、条内名字块变**白底红边红字**（选中才是蓝底白字）。
    let hover_src = resp.hover_pos().and_then(|p| {
        feat_hits
            .iter()
            .filter(|(r, _, is_label)| *is_label && r.expand(1.0).contains(p))
            .next()
            .or_else(|| feat_hits.iter().find(|(r, _, _)| r.contains(p)))
            .map(|(_, source, _)| *source)
    });
    if let Some(p) = resp.hover_pos() {
        app.hover_bp = bp_at(p);
        if residue_hits.iter().any(|(r, ..)| r.contains(p))
            || feat_hits.iter().any(|(r, ..)| r.contains(p))
            || primer_hits.iter().any(|(r, _)| r.contains(p))
        {
            ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
        }
    } else {
        app.hover_bp = None;
    }

    // SnapGene 的右键菜单（SequenceContextMenu）
    resp.context_menu(|ui| app.seq_context_menu(ui));

    let selected_bases = app.selected_residue.and_then(|(source, k)| {
        let cds = feats.iter().find(|f| f.source == source)?.cds.as_ref()?;
        let bases = cds.codon_bases(k)?;
        (codon_selection(bases, len) == app.selection && app.cursor_bp.is_none()).then_some(bases)
    });
    if selected_bases.is_none() {
        app.selected_residue = None;
    }

    // 内容区裁剪（不越界到顶部刻度/底部概览）
    let p = painter.with_clip_rect(Rect::from_min_max(
        Pos2::new(content_left, content_top),
        Pos2::new(content_right, content_bottom),
    ));
    let p = &p;
    let font = FontId::monospace((12.0 * z).max(7.0));
    let aa_font = FontId::proportional((10.5 * z).max(7.0));
    // ---- 逐行绘制 ----
    for row in &rows {
        let top = row_y(row);
        if top > content_bottom + 60.0 || top + row.h < content_top - 60.0 {
            continue;
        }
        let y_top_c = top + strand_h * 0.5;
        let y_line = top + strand_h;
        let y_bot_c = top + strand_h * 1.5;
        let after_dna = top + strand_h * 2.0;
        // 特征/ORF 轨道结束的位置（反向引物的 band 从这里往下排）
        let _after_tracks = after_dna + row.features_h + row.orfs_h;

        // 碱基（双链）
        // 选中特征的碱基先铺一层蓝底（与单链序列视图同一套：`#B5D5FF`）。
        if let Some(sfi) = app.selected_feature {
            for track in &row.tracks {
                let f = &feats[track.feature];
                if f.orf || f.source != sfi {
                    continue;
                }
                for span in &track.spans {
                    let x0 = content_left + v.x_of(span.start);
                    let x1 = content_left + v.x_of(span.end - 1) + base_w;
                    p.rect_filled(
                        Rect::from_min_max(
                            Pos2::new(x0, top + 3.0 * z),
                            Pos2::new(
                                x1,
                                top
                                    + if app.show_complement {
                                        2.0 * strand_h
                                    } else {
                                        strand_h
                                    }
                                    - 3.0 * z,
                            ),
                        ),
                        0.0,
                        theme::FEAT_SEL_BG,
                    );
                }
            }
        }
        let mut max_x = content_left + v.pad_left;
        for bp in row.start..row.end {
            let x = content_left + v.x_of(bp);
            let cx = x + base_w * 0.5;
            max_x = x + base_w;
            let b = seq[bp];
            let up = (kind.display_base(b) as char).to_ascii_uppercase();
            // 碱基保持中性墨色；Colors 控制注释图形的颜色。
            let top_col = if app.show_colors { sequence_color(b, app.sequence_color_mode) } else { theme::TEXT };
            p.text(
                Pos2::new(cx, y_top_c),
                Align2::CENTER_CENTER,
                up.to_string(),
                font.clone(),
                top_col,
            );
            if app.show_complement && bp >= bottom_start && bp < bottom_end {
                let bot_col = theme::BASE_COMPLEMENT;
                p.text(
                    Pos2::new(cx, y_bot_c),
                    Align2::CENTER_CENTER,
                    complement_display(kind, b).to_string(),
                    font.clone(),
                    bot_col,
                );
            }
        }

        // 底链独有的悬垂碱基绘在顶链坐标网格之外，不计入长度/光标坐标。
        if app.show_complement {
            if row.start == 0 {
                let x0 = content_left + v.x_of(0) - left_bottom_tail.len() as f32 * base_w;
                for (i, &base) in left_bottom_tail.iter().enumerate() {
                    p.text(Pos2::new(x0 + (i as f32 + 0.5) * base_w, y_bot_c),
                        Align2::CENTER_CENTER, complement_display(kind, base).to_string(), font.clone(), theme::BASE_COMPLEMENT);
                }
            }
            if row.end == len {
                for (i, &base) in right_bottom_tail.iter().enumerate() {
                    p.text(Pos2::new(max_x + (i as f32 + 0.5) * base_w, y_bot_c),
                        Align2::CENTER_CENTER, complement_display(kind, base).to_string(), font.clone(), theme::BASE_COMPLEMENT);
                }
            }
        }

        // 骨架线：每碱基小刻度、每 10 碱基长刻度（SnapGene 的链线）
        let x0 = content_left + v.x_of(row.start);
        p.line_segment(
            [Pos2::new(x0, y_line), Pos2::new(max_x, y_line)],
            Stroke::new(1.0, theme::SEQ_LINE),
        );
        let mut bp = row.start;
        while bp < row.end {
            let x = content_left + v.x_of(bp);
            let long = bp % 10 == 0;
            p.line_segment(
                [
                    Pos2::new(x, y_line - if long { 3.2 } else { 1.8 }),
                    Pos2::new(x, y_line + if long { 3.2 } else { 1.8 }),
                ],
                Stroke::new(if long { 0.9 } else { 0.6 }, theme::SEQ_LINE),
            );
            bp += 1;
        }

        // 行末绝对坐标
        p.text(
            Pos2::new(max_x + 8.0, y_line - 2.0),
            Align2::LEFT_CENTER,
            nfmt_plain(row.end),
            small_font.clone(),
            theme::RULER_TEXT,
        );

        // 行首 5'/3'
        p.text(
            Pos2::new(content_left + pad_left - 6.0 - app.seq_pan_x, y_top_c),
            Align2::RIGHT_CENTER,
            "5'",
            small_font.clone(),
            theme::TEXT_DIM,
        );
        if app.show_complement {
            p.text(
                Pos2::new(content_left + pad_left - 6.0 - app.seq_pan_x
                    - if row.start == 0 { left_bottom_tail.len() as f32 * base_w } else { 0.0 }, y_bot_c),
                Align2::RIGHT_CENTER,
                "3'",
                small_font.clone(),
                theme::TEXT_DIM,
            );
        }

        // ---- 引物结合位点 ----
        //
        // 版式按真 SnapGene 8.x 的序列视图实测（`tail4.dna` 的 T5F/T3F/T5R/T3R）：
        //   * 退火段是一个**透明底 + 引物紫（`#A020F0`）描边的方框**，框里是该段
        //     「显示链」的碱基，引物名画在框外侧（紫色）；
        //   * 不与模板配对的那一段（尾巴）**另画一个同高的框**：正向引物朝上错开、
        //     反向引物朝下错开 `step`，并且**横向压进退火框 `step` 宽**（两框共用
        //     那一小段边界，SnapGene 的台阶形状就是这么来的）；尾巴文字按显示方向排；
        //   * 3′ 端的斜旗是 45° 平行四边形，画在**整条引物 3′ 端那个框**的外侧 ——
        //     有 3′ 尾巴时是尾巴框，否则是退火框（实测 T3F / T3R）；
        //   * 5′ 磷酸化在 5′ 端外侧画一个圈起来的 `P`。
        for span in &row.primers {
            let primer = &primers[span.primer];
            let picked = app.selected_primer == Some(span.primer);
            let color = if app.show_colors {
                theme::PRIMER
            } else {
                theme::FEATURE_NO_COLOR
            };
            let g = primer_boxes(
                span, primer, row, top, content_left, base_w, primer_h, bar_gap, strand_h,
            );
            let PrimerBoxes { y, x0, x1, arect, trect } = g;
            let mid = y + primer_h * 0.5;
            let stroke = Stroke::new(1.2 * z, color);
            // 选中的引物框铺 SnapGene 的选中蓝底（oracle `t15` 实测 `#B5D5FF`）。
            let fill = if picked { theme::FEAT_SEL_BG } else { theme::CANVAS };
            let tip_font = font.clone();

            // ---- 退火段 + 尾巴（画成一个并集轮廓）----
            p.rect_filled(arect, 0.0, fill);
            if let Some(trect) = trect {
                p.rect_filled(trect, 0.0, fill);
            }
            let outline = match trect {
                Some(t) => primer_union_outline(arect, t, primer.tail_left, primer.forward),
                None => vec![
                    Pos2::new(x0, y),
                    Pos2::new(x1, y),
                    Pos2::new(x1, y + primer_h),
                    Pos2::new(x0, y + primer_h),
                ],
            };
            p.add(Shape::closed_line(outline, stroke));

            // ---- 3′ 端的斜旗（45° 平行四边形，贴在 3′ 端那个框的外侧）----
            if span.tip {
                let (ex, ey, up) = if primer.forward {
                    let end = match trect {
                        Some(r) if !primer.tail_left => r,
                        _ => arect,
                    };
                    (end.right(), end.top(), true)
                } else {
                    let end = match trect {
                        Some(r) if primer.tail_left => r,
                        _ => arect,
                    };
                    (end.left(), end.bottom(), false)
                };
                let w = 1.5 * base_w;
                let d = 1.3 * base_w;
                let hh = 1.2 * base_w;
                let flag = if up {
                    vec![
                        Pos2::new(ex - w, ey),
                        Pos2::new(ex, ey),
                        Pos2::new(ex - d, ey - hh),
                        Pos2::new(ex - w - d, ey - hh),
                    ]
                } else {
                    vec![
                        Pos2::new(ex, ey),
                        Pos2::new(ex + w, ey),
                        Pos2::new(ex + w + d, ey + hh),
                        Pos2::new(ex + d, ey + hh),
                    ]
                };
                p.add(Shape::convex_polygon(flag, theme::CANVAS, stroke));
            }

            // ---- 尾巴文字（碱基已按「显示方向」排好）----
            //
            // 字母落在**模板网格**上：SnapGene 实测（pET28b-Rep1.dna 610..640
            // 的 `CCC` 尾巴）尾巴三个字母与退火框左侧紧邻的三位模板字母逐像素
            // 同列 —— 尾巴只是被抬到另一行，横向仍是连续序列。不能按尾巴框居中
            // 排（框比网格宽 `step`，会把字母推开半格）。
            if let Some(trect) = trect {
                let chars: Vec<char> = primer.tail.chars().collect();
                let tx0 = tail_text_left(x0, x1, chars.len(), base_w, primer.tail_left);
                let ty_c = trect.center().y;
                // 尾巴超出可视区时按 SnapGene 的做法截断并补 `…`。
                let clip = p.clip_rect();
                let mut clipped = false;
                for (i, ch) in chars.iter().enumerate() {
                    let cx = tx0 + (i as f32 + 0.5) * base_w;
                    if cx < clip.left() + base_w * 0.5 || cx > clip.right() - base_w * 0.5 {
                        clipped = true;
                        continue;
                    }
                    p.text(
                        Pos2::new(cx, ty_c),
                        Align2::CENTER_CENTER,
                        ch.to_string(),
                        tip_font.clone(),
                        theme::TEXT,
                    );
                }
                if clipped {
                    let ex = if primer.tail_left {
                        clip.left() + 2.0 * z
                    } else {
                        clip.right() - 2.0 * z
                    };
                    p.text(
                        Pos2::new(ex, ty_c),
                        if primer.tail_left {
                            Align2::LEFT_CENTER
                        } else {
                            Align2::RIGHT_CENTER
                        },
                        "…",
                        tip_font.clone(),
                        theme::TEXT,
                    );
                }
            }
            // 框里的碱基：与序列行同一列的**显示链**字符
            for bp in row.start..row.end {
                let bx = content_left + v.x_of(bp);
                if bx + base_w <= x0 + 0.5 || bx >= x1 - 0.5 {
                    continue;
                }
                let ch = if primer.forward {
                    (kind.display_base(seq[bp]) as char).to_ascii_uppercase()
                } else {
                    complement_display(kind, seq[bp]).to_ascii_uppercase()
                };
                p.text(
                    Pos2::new(bx + base_w * 0.5, mid),
                    Align2::CENTER_CENTER,
                    ch.to_string(),
                    font.clone(),
                    theme::TEXT,
                );
            }

            // ---- 引物名（上链画在框上方、下链画在框下方，紫色）----
            let name_size = p
                .layout_no_wrap(primer.name.clone(), small_font.clone(), color)
                .size();
            if primer.forward {
                let c = Pos2::new((x0 + x1) * 0.5, y - 2.0 * z);
                p.text(c, Align2::CENTER_BOTTOM, primer.name.clone(), small_font.clone(), color);
                primer_hits.push((
                    Rect::from_min_size(
                        Pos2::new(c.x - name_size.x * 0.5, c.y - name_size.y),
                        name_size,
                    ),
                    span.primer,
                ));
            } else {
                let c = Pos2::new((x0 + x1) * 0.5, y + primer_h + 1.5 * z);
                p.text(c, Align2::CENTER_TOP, primer.name.clone(), small_font.clone(), color);
                primer_hits.push((
                    Rect::from_min_size(Pos2::new(c.x - name_size.x * 0.5, c.y), name_size),
                    span.primer,
                ));
            }

            // ---- 5′ 磷酸化标记 Ⓟ ----
            if span.five_prime && primer.phosphorylated {
                // 尾巴挡在 5′ 端外侧时把 Ⓟ 再往外让一格；尾巴在另一头就不用让。
                let tail_w = if span.tail_here {
                    primer.tail.len() as f32 * base_w
                } else {
                    0.0
                };
                let c = if primer.forward {
                    Pos2::new(x0 - tail_w - 9.0 * z, y + primer_h * 0.4)
                } else {
                    Pos2::new(x1 + tail_w + 9.0 * z, y + primer_h * 0.6)
                };
                p.circle_stroke(c, 6.0 * z, Stroke::new(1.1 * z, color));
                p.text(
                    c,
                    Align2::CENTER_CENTER,
                    "P",
                    FontId::monospace((9.0 * z).max(5.0)),
                    color,
                );
            }
        }

        // 反向引物的 band 在特征/ORF 轨道下方（SnapGene 实测），
        // 所以特征从这里直接开始。
        let after_primers = after_dna;
        // ---- 每个特征自己的标尺、氨基酸和特征条 ----
        for track in &row.tracks {
            let f = &feats[track.feature];
            let level = &row.levels[track.level];
            let track_top = after_primers + level.top;
            let bar_y = track_top + if level.translation { aa_h } else { 0.0 };
            let fill = if app.show_colors {
                f.color
            } else {
                theme::FEATURE_NO_COLOR
            };
            let edge = theme::mix(fill, Color32::BLACK, 0.6);
            // 选中特征：条只描边变蓝、填充保持本色（oracle：SnapGene 的选中效果）
            let picked = app.selected_feature == Some(f.source) && !f.orf;
            let hovered = !picked && !f.orf && hover_src == Some(f.source);
            let edge = if picked {
                theme::FEAT_SEL_EDGE
            } else if hovered {
                theme::SITE_HOVER
            } else {
                edge
            };
            let edge_w = if picked || hovered { 1.6 } else { 0.8 };
            for span in &track.spans {
                let x0 = content_left + v.x_of(span.start);
                let x1 = content_left + v.x_of(span.end - 1) + base_w;
                let bar = Rect::from_min_max(Pos2::new(x0, bar_y), Pos2::new(x1, bar_y + bar_h));
                p.add(Shape::convex_polygon(
                    feature_polygon(bar, span.left_arrow, span.right_arrow),
                    fill,
                    Stroke::new(edge_w, edge),
                ));
                let label = p.layout_no_wrap(f.name.clone(), small_font.clone(), ink_on(fill));
                if label.size().x + 10.0 * z <= bar.width() {
                    // 选中：名字画成蓝底白字圆角块（与单链视图、图谱里的选中标签一致）
                    if picked {
                        let chip = Rect::from_center_size(
                            bar.center(),
                            label.size() + Vec2::new(10.0 * z, 5.0 * z),
                        );
                        p.rect_filled(chip, 3.0 * z, theme::FEAT_SEL_EDGE);
                        p.galley(
                            bar.center() - label.size() * 0.5,
                            p.layout_no_wrap(
                                f.name.clone(),
                                small_font.clone(),
                                Color32::WHITE,
                            ),
                            Color32::WHITE,
                        );
                    } else if hovered {
                        // oracle `dn7-chip.png`：悬停 = 白底 + 红边 + 红字（不是选中蓝块）
                        let chip = Rect::from_center_size(
                            bar.center(),
                            label.size() + Vec2::new(10.0 * z, 5.0 * z),
                        );
                        p.rect_filled(chip, 3.0 * z, Color32::WHITE);
                        p.rect_stroke(
                            chip,
                            3.0 * z,
                            Stroke::new(1.0, theme::SITE_HOVER),
                            egui::StrokeKind::Inside,
                        );
                        p.galley(
                            bar.center() - label.size() * 0.5,
                            p.layout_no_wrap(
                                f.name.clone(),
                                small_font.clone(),
                                theme::SITE_HOVER,
                            ),
                            theme::SITE_HOVER,
                        );
                    } else {
                        p.galley(bar.center() - label.size() * 0.5, label, ink_on(fill));
                    }
                } else if track.external_label {
                    let label = p.layout_no_wrap(f.name.clone(), small_font.clone(), theme::TEXT);
                    p.line_segment(
                        [
                            bar.center_bottom(),
                            bar.center_bottom() + Vec2::new(0.0, 3.0 * z),
                        ],
                        Stroke::new(0.7, edge),
                    );
                    let lpos = Pos2::new(
                        bar.center().x - label.size().x * 0.5,
                        bar.bottom() + 2.0 * z,
                    );
                    if picked {
                        let chip = Rect::from_min_size(lpos, label.size())
                            .expand2(Vec2::new(5.0 * z, 2.5 * z));
                        p.rect_filled(chip, 3.0 * z, theme::FEAT_SEL_EDGE);
                        p.galley(
                            lpos,
                            p.layout_no_wrap(
                                f.name.clone(),
                                small_font.clone(),
                                Color32::WHITE,
                            ),
                            Color32::WHITE,
                        );
                    } else if hovered {
                        // 条外名字：红字 + 红圆角框（与蛋白序列视图的悬停标签同一套）
                        let chip = Rect::from_min_size(lpos, label.size())
                            .expand2(Vec2::new(5.0 * z, 2.5 * z));
                        p.rect_stroke(
                            chip,
                            3.0 * z,
                            Stroke::new(1.0, theme::SITE_HOVER),
                            egui::StrokeKind::Inside,
                        );
                        p.galley(
                            lpos,
                            p.layout_no_wrap(
                                f.name.clone(),
                                small_font.clone(),
                                theme::SITE_HOVER,
                            ),
                            theme::SITE_HOVER,
                        );
                    } else {
                        p.galley(lpos, label, theme::TEXT);
                    }
                }
                // 跨行延续用条外空心折线；真正的方向尖端只画在特征终点。
                for (yes, x, dir) in [
                    (span.continues_left, x0, -1.0),
                    (span.continues_right, x1, 1.0),
                ] {
                    if yes {
                        p.add(Shape::line(
                            vec![
                                Pos2::new(x + dir * 4.0 * z, bar.top() + 1.0),
                                Pos2::new(x + dir * 9.0 * z, bar.center().y),
                                Pos2::new(x + dir * 4.0 * z, bar.bottom() - 1.0),
                            ],
                            Stroke::new(0.8, theme::TEXT_DIM),
                        ));
                    }
                }
            }
            let Some(cds) = &f.cds else { continue };
            let run_on: Vec<_> = track
                .codons
                .iter()
                .filter(|c| c.codon >= cds.annotated_codons && cds.aa[c.codon] != "***")
                .collect();
            if let (Some(start), Some(end)) = (
                run_on.iter().map(|c| c.start).min(),
                run_on.iter().map(|c| c.end).max(),
            ) {
                let x0 = content_left + v.x_of(start);
                let x1 = content_left + v.x_of(end - 1) + base_w;
                let y = bar_y + bar_h * 0.5;
                let mut x = x0;
                while x < x1 - 5.0 * z {
                    p.line_segment(
                        [
                            Pos2::new(x, y),
                            Pos2::new((x + 3.0 * z).min(x1 - 5.0 * z), y),
                        ],
                        Stroke::new(0.8, theme::TEXT_DIM),
                    );
                    x += 6.0 * z;
                }
                p.add(Shape::convex_polygon(
                    vec![
                        Pos2::new(x1 - 5.0 * z, y - 2.0 * z),
                        Pos2::new(x1, y),
                        Pos2::new(x1 - 5.0 * z, y + 2.0 * z),
                    ],
                    theme::TEXT_DIM,
                    Stroke::NONE,
                ));
            }
            let y_ruler = track_top + 13.0 * z;
            let y_aa = track_top + 22.0 * z;
            // 标尺按实际片段绘制，剪接间隙不会被一条虚构的实线连接。
            for span in &track.spans {
                p.line_segment(
                    [
                        Pos2::new(content_left + v.x_of(span.start), y_ruler),
                        Pos2::new(content_left + v.x_of(span.end - 1) + base_w, y_ruler),
                    ],
                    Stroke::new(0.7, theme::TEXT_DIM),
                );
            }
            for fragment in &track.codons {
                let k = fragment.codon;
                let (cell, cx) = codon_rect(fragment, &v, content_left, y_aa, z);
                let cp = p.with_clip_rect(p.clip_rect().intersect(cell));
                let picked = app.selected_residue == Some((f.source, k));
                if picked {
                    cp.rect_filled(cell, 0.0, theme::SELECTION);
                }
                let aa = cds.aa[k];
                if aa == "***" {
                    let marker = Rect::from_center_size(Pos2::new(cx, y_aa), Vec2::splat(8.0 * z));
                    cp.rect_filled(marker, 0.0, Color32::RED);
                    cp.text(
                        marker.center(),
                        Align2::CENTER_CENTER,
                        "*",
                        small_font.clone(),
                        Color32::WHITE,
                    );
                } else {
                    let label = if app.aa_one_letter {
                        hx_algo::translate::one_letter(aa)
                    } else {
                        aa
                    };
                    cp.text(
                        Pos2::new(cx, y_aa),
                        Align2::CENTER_CENTER,
                        label,
                        aa_font.clone(),
                        if picked {
                            Color32::WHITE
                        } else if k == 0 && aa == "Met" {
                            Color32::RED
                        } else {
                            theme::TEXT
                        },
                    );
                    let idx = k + 1;
                    // 编号仍按编码方向递增，所以反链在屏幕上由左向右递减。
                    if k < cds.annotated_codons && (idx == 1 || idx % AA_LABEL_EVERY == 0) {
                        p.text(
                            Pos2::new(cx, y_ruler - 1.0 * z),
                            Align2::CENTER_BOTTOM,
                            nfmt(idx),
                            FontId::proportional((9.0 * z).max(7.0)),
                            theme::RULER_TEXT,
                        );
                    }
                }
                let tick_x = if cds.forward {
                    cell.left()
                } else {
                    cell.right()
                };
                p.line_segment(
                    [
                        Pos2::new(tick_x, y_ruler - 3.0 * z),
                        Pos2::new(tick_x, y_ruler),
                    ],
                    Stroke::new(0.6, theme::TEXT_DIM),
                );
            }
        }

        // ---- 独立 ORF：浅橙/浅绿氨基酸带，无特征名称、编号和标尺 ----
        for track in &row.orfs {
            let f = &feats[track.feature];
            let cds = f.cds.as_ref().unwrap();
            let aa_y = after_primers
                + row.features_h
                + 6.0 * z
                + (track.level as f32 + 0.5) * ORF_ROW_H * z;
            let fill = orf_fill(cds.forward, app.show_colors);
            for fragment in &track.codons {
                let k = fragment.codon;
                let (cell, cx) = codon_rect(fragment, &v, content_left, aa_y, z);
                let cp = p.with_clip_rect(p.clip_rect().intersect(cell));
                let picked = app.selected_residue == Some((f.source, k));
                cp.rect_filled(cell, 0.0, if picked { theme::SELECTION } else { fill });
                let aa = cds.aa[k];
                if aa == "***" {
                    let marker = Rect::from_center_size(Pos2::new(cx, aa_y), Vec2::splat(8.0 * z));
                    cp.rect_filled(marker, 0.0, Color32::RED);
                    cp.text(
                        marker.center(),
                        Align2::CENTER_CENTER,
                        "*",
                        small_font.clone(),
                        Color32::WHITE,
                    );
                } else {
                    cp.text(
                        Pos2::new(cx, aa_y),
                        Align2::CENTER_CENTER,
                        if app.aa_one_letter {
                            hx_algo::translate::one_letter(aa)
                        } else {
                            aa
                        },
                        aa_font.clone(),
                        if picked { Color32::WHITE } else { theme::TEXT },
                    );
                }
            }
            // 方向提示放在碱基区域外侧，不伪装成一个新的特征箭头。
            let (x, dir) = if cds.forward {
                (content_left + v.x_of(row.end - 1) + base_w + 16.0 * z, 1.0)
            } else {
                (content_left + v.x_of(row.start) - 16.0 * z, -1.0)
            };
            let color = theme::darken(fill, 0.36);
            p.arrow(
                Pos2::new(x - dir * 5.0 * z, aa_y),
                Vec2::new(dir * 10.0 * z, 0.0),
                Stroke::new(1.6 * z, color),
            );
        }
    }

    // ---- Enzyme label lanes above the sequence, never suppress by count. ----
    let mut hovered_enzyme:Option<usize>=None;
    let mut hovered_name:Option<String>=None;
    for(row_index,labels)in enzyme_labels.iter().enumerate(){
        let Some(row)=rows.get(row_index)else{continue};let top=row_y(row);
        let offsets=crate::enzyme_labels::lane_offsets(labels,&sites);
        for label in labels{
            let(name,_,unique)=&sites[label.site];
            let anchor=Pos2::new(content_left+label.left,top-row.fwd_h-12.0*z-offsets[label.lane+1]as f32*16.0*z);
            let hit=egui::Rect::from_min_size(anchor,egui::vec2(label.width,name.lines().count()as f32*16.0*z));
            for(line,enzyme)in name.lines().enumerate(){
                let target=Rect::from_min_size(hit.min+egui::vec2(0.0,line as f32*16.0*z),egui::vec2(hit.width(),16.0*z));
                let response=ui.interact(target,ui.scope_id().with(("enzyme_label",row_index,label.site,line)),Sense::hover());
                if response.hovered(){hovered_enzyme=Some(label.site);hovered_name=Some(enzyme.to_owned());}
                response.on_hover_ui(|ui|{ui.strong(format!("{} ({})",enzyme,sites[label.site].1+1));if let Some(e)=hx_enzymes::by_name(enzyme){ui.label(crate::ui_format!(ui.ctx(),"{} sites · {} bp",sites.iter().filter(|s|s.0.lines().any(|n|n==enzyme)).count(),e.site_len()));ui.label(crate::i18n::text(ui.ctx(),e.overhang.label()));ui.separator();ui.monospace(format!("5′  {}  3′",e.site));let complement=e.site.bytes().map(hx_algo::edit_session::complement).map(char::from).collect::<String>();ui.monospace(format!("3′  {complement}  5′"));ui.label(crate::ui_format!(ui.ctx(),"Cut: {} / {}",e.cut,e.site_len()as i16-e.reverse_cut));if let Some(i)=hx_enzymes::catalog::CATALOG.iter().position(|x|x.name==e.name){if hx_enzymes::details::DETAILS[i].variable_ends{ui.label(crate::i18n::text(ui.ctx(),"Sticky ends from different sites may not be compatible."));}}}ui.label(crate::i18n::text(ui.ctx(),if *unique{"Unique cutter"}else{"Multiple sites"}));});
            }
        }
    }
    if let Some(&hovered)=hovered_enzyme.as_ref(){
        let name=&sites[hovered].0;
        if let Some(index)=hx_enzymes::catalog::CATALOG.iter().position(|e|Some(e.name)==name.lines().next()){
            for(_,start,_)in std::iter::once(&sites[hovered]){
                let geometry=crate::enzyme_labels::geometry(index,*start,&seq,sequence_circular);
                for(k,b)in geometry.pattern.bytes().enumerate(){if b==b'N'{continue;}let at=(*start+k)%len;let ri=v.row_of(at);let Some(row)=rows.get(ri)else{continue};let y=row_y(row);let x=content_left+v.x_of(at);
                    for dy in [0.0,strand_h]{if y+dy>=content_top&&y+dy<content_bottom{p.rect_filled(Rect::from_min_size(Pos2::new(x,y+dy),egui::vec2(base_w,strand_h-5.0*z)),0.0,Color32::from_black_alpha(28));}}
                }
                for(top,bottom)in geometry.cuts{
                    let normalize=|cut:i64|if sequence_circular{Some(cut.rem_euclid(len as i64)as usize)}else if cut>=0&&cut<=len as i64{Some(cut as usize)}else{None};
                    let upper=top.and_then(normalize);let lower=bottom.and_then(normalize);
                    let at_point=|at:usize|{let ri=v.row_of(at.min(len-1));rows.get(ri).map(|row|(ri,content_left+v.x_of(at),row_y(row)))};
                    if let(Some(a),Some(b))=(upper,lower){if let(Some((ra,xa,ya)),Some((rb,xb,_)))=(at_point(a),at_point(b)){
                        if ra==rb&&ya>=content_top&&ya<content_bottom{
                            let middle=ya+strand_h;
                            p.add(egui::Shape::line(vec![Pos2::new(xa,ya+3.0*z),Pos2::new(xa,middle),Pos2::new(xb,middle),Pos2::new(xb,ya+strand_h*2.0-4.0*z)],Stroke::new(1.6,Color32::RED)));
                            p.arrow(Pos2::new(xa,ya-9.0*z),egui::vec2(0.0,10.0*z),Stroke::new(1.5,Color32::RED));
                            p.arrow(Pos2::new(xb,ya+strand_h*2.0+7.0*z),egui::vec2(0.0,-10.0*z),Stroke::new(1.5,Color32::RED));
                        }
                    }}else{
                        for at in [upper,lower].into_iter().flatten(){if let Some((_,x,y))=at_point(at){if y>=content_top&&y<content_bottom{p.arrow(Pos2::new(x,y-9.0*z),egui::vec2(0.0,10.0*z),Stroke::new(1.5,Color32::RED));}}}
                    }
                }
            }
        }
    }
    for(row_index,labels)in enzyme_labels.iter().enumerate(){
        let Some(row)=rows.get(row_index)else{continue};let top=row_y(row);
        if top<content_top-20.0||top-row.fwd_h-labels.len()as f32*14.0*z>content_bottom{continue;}
        let offsets=crate::enzyme_labels::lane_offsets(labels,&sites);
        for label in labels{
            let(name,_,unique)=&sites[label.site];let highlight=hovered_enzyme==Some(label.site);
            let anchor=Pos2::new(content_left+label.left,top-row.fwd_h-12.0*z-offsets[label.lane+1]as f32*16.0*z);
            let hit=egui::Rect::from_min_size(anchor,egui::vec2(label.width,name.lines().count()as f32*16.0*z));
            if hit.bottom()<content_top||hit.top()>content_bottom{continue;}
            let color=if highlight{Color32::RED}else if app.preferences.green_blunt&&hx_enzymes::by_name(name).is_some_and(|e|e.overhang==hx_enzymes::Overhang::Blunt){Color32::from_rgb(0,150,45)}else{theme::TEXT};
            let leader_x=if highlight{
                let start=sites[label.site].1;
                hx_enzymes::catalog::CATALOG.iter().position(|e|Some(e.name)==name.lines().next())
                    .and_then(|index|crate::enzyme_labels::geometry(index,start,&seq,sequence_circular).cuts.first().and_then(|c|c.0))
                    .and_then(|at|if sequence_circular{Some(at.rem_euclid(len as i64)as usize)}else{usize::try_from(at).ok().filter(|&a|a<=len)})
                    .filter(|&at|v.row_of(at.min(len-1))==row_index).map(|at|content_left+v.x_of(at)).unwrap_or(content_left+label.x)
            }else{content_left+label.x};
            let endpoint=Pos2::new(leader_x,top-3.0*z);
            let text_edge=if hit.right()<endpoint.x{hit.right()+2.0*z}else{hit.left()-2.0*z};
            p.add(egui::Shape::line(vec![Pos2::new(text_edge,hit.top()+7.0*z),Pos2::new(endpoint.x,hit.top()+13.0*z),endpoint],Stroke::new(if highlight{1.2}else{0.6},if highlight{Color32::RED}else{theme::TEXT_DIM})));
            if highlight{let line=name.lines().position(|n|Some(n)==hovered_name.as_deref()).unwrap_or(0);let target=Rect::from_min_size(hit.min+egui::vec2(0.0,line as f32*16.0*z),egui::vec2(hit.width(),16.0*z));p.rect_stroke(target.expand(1.0),2.0,Stroke::new(0.8,Color32::RED),egui::StrokeKind::Inside);p.arrow(endpoint-egui::vec2(0.0,7.0),egui::vec2(0.0,7.0),Stroke::new(1.2,Color32::RED));}
            let mut job=egui::text::LayoutJob::default();
            for(line,enzyme)in name.lines().enumerate(){if line>0{job.append("\n",0.0,egui::TextFormat::simple(FontId::proportional((12.0*z).max(10.0)),theme::TEXT));}let line_color=if highlight&&hovered_name.as_deref()==Some(enzyme){Color32::RED}else{theme::TEXT};let mut format=egui::TextFormat::simple(FontId::proportional((12.0*z).max(10.0)),line_color);format.italics=app.preferences.italic_enzymes;job.append(enzyme,0.0,format);}
            let galley=p.layout_job(job);
            if *unique&&app.preferences.bold_unique{p.galley(anchor+egui::vec2(0.4,0.0),galley.clone(),color);}p.galley(anchor,galley,color);
        }
    }

    // ---- 选段 / 光标 ----
    let digest=app.digest_selection.as_ref().filter(|d|Some(d.range)==app.selection&&selected_bases.is_none());
    if let Some(d)=digest.filter(|_|len>0) {
        for (start,end,lower) in [(d.top_start,d.top_end,false),(d.bottom_start,d.bottom_end,true)] {
            for(s,e)in duplex_selection_parts(start,end,len){
                for(row,x0,_)in v.rows_of_span(s,e,len){let Some(r)=rows.get(row)else{continue};let top=row_y(r);let x1=v.x_of(e.min(r.end)-1)+base_w;
                    let y=top+if lower{strand_h}else{0.0};
                    p.rect_filled(Rect::from_min_max(Pos2::new(content_left+x0-1.0,y+1.0),Pos2::new(content_left+x1+1.0,y+strand_h-1.0)),0.0,theme::with_alpha(theme::SELECTION,60));
                }
            }
            // One marker per actual cut, not per wrapped/line-split highlight.
            for gap in [start,end]{let gap=if sequence_circular{gap.rem_euclid(len as i64)as usize}else{gap.clamp(0,len as i64)as usize};let at=gap.min(len-1);let row=v.row_of(at);let Some(r)=rows.get(row)else{continue};let top=row_y(r);let x=content_left+v.x_of(at)+if gap==len{base_w}else{0.0};
                let y=if lower{top+strand_h*2.0+2.0*z}else{top-2.0*z};let direction=if lower{1.0}else{-1.0};
                p.add(egui::Shape::convex_polygon(vec![Pos2::new(x,y),Pos2::new(x-4.5*z,y+direction*7.0*z),Pos2::new(x+4.5*z,y+direction*7.0*z)],theme::ACCENT,Stroke::NONE));
            }
        }
    }
    let selection_spans = if digest.is_some(){Vec::new()}else if let Some(bases) = selected_bases {
        codon_parts(bases)
    } else {
        app.selection.into_iter().collect()
    };
    for (s, e) in selection_spans {
        for (row, x0, mut x1) in v.rows_of_span(s, e, len) {
            let Some(r) = rows.get(row) else { continue };
            if s < e {
                x1 = v.x_of(e.min(r.end) - 1) + base_w;
            }
            let top = row_y(r);
            let rect_sel = Rect::from_min_max(
                Pos2::new(content_left + x0 - 1.0, top + 1.0),
                Pos2::new(content_left + x1 + 1.0, top + strand_h * 2.0 - 1.0),
            );
            p.rect_filled(rect_sel, 0.0, theme::with_alpha(theme::SELECTION, 60));
            p.rect_stroke(
                rect_sel,
                0.0,
                Stroke::new(1.0, theme::with_alpha(theme::SELECTION, 160)),
                egui::StrokeKind::Inside,
            );
        }
    }
    if let Some(cbp) = app.cursor_bp {
        if len > 0 {
            // 光标 = **插入点**（在 `cbp` 那一位碱基之前）；`len` 表示序列末尾之后
            let at = cbp.min(len - 1);
            let row = v.row_of(at);
            if let Some(r) = rows.get(row) {
                let top = row_y(r);
                let x = content_left + v.x_of(at) + if cbp >= len { base_w } else { 0.0 };
                p.line_segment(
                    [
                        Pos2::new(x, top + 1.0),
                        Pos2::new(x, top + strand_h * 2.0 - 1.0),
                    ],
                    Stroke::new(1.6, Color32::from_rgb(0xe5, 0x2b, 0x2b)),
                );
            }
        }
    }

    // ---- 顶部固定刻度（行内相对坐标，绝对坐标见行末）----
    let ruler = Rect::from_min_size(
        Pos2::new(content_left, rect.top()),
        Vec2::new(content_right - content_left, ruler_h),
    );
    // 序列页整块底色是 `#F2F2F2`（oracle 实测），都用它，别在标尺带留白条。
    painter.rect_filled(ruler, 0.0, theme::SEQ_BG);
    let y_r = rect.top() + ruler_h - 6.0;
    painter.line_segment(
        [
            Pos2::new(content_left + pad_left, y_r),
            Pos2::new(content_right, y_r),
        ],
        Stroke::new(1.0, theme::ACCENT),
    );
    if app.show_ruler {
        // Dense blue minor ticks plus numbered 10 bp major ticks.
        for c in 0..=v.line_bp {
            let x = if c == v.line_bp {
                content_left + v.x_of(v.line_bp.saturating_sub(1)) + base_w
            } else {
                content_left + v.x_of(c)
            };
            if x > content_right + 1.0 {
                continue;
            }
            let major = c > 0 && c % v.group_bp == 0;
            let half = if major {
                6.0
            } else if c % 5 == 0 {
                4.0
            } else {
                2.5
            };
            painter.line_segment(
                [Pos2::new(x, y_r - half), Pos2::new(x, y_r + 1.0)],
                Stroke::new(if major { 1.0 } else { 0.65 }, theme::ACCENT),
            );
            if major {
                painter.text(
                    Pos2::new(x - 2.0, y_r - 2.0),
                    Align2::RIGHT_BOTTOM,
                    nfmt(c),
                    small_font.clone(),
                    theme::RULER_TEXT,
                );
            }
        }
    }

    // ---- 选区信息条（SnapGene 的 SelectionBar）----
    {
        let top = rect.bottom() - BOTTOM_H - sel_bar_h;
        let bar = Rect::from_min_max(
            Pos2::new(rect.left(), top),
            Pos2::new(rect.right(), top + sel_bar_h),
        );
        painter.rect_filled(bar, 0.0, theme::PANEL_ALT);
        painter.line_segment(
            [Pos2::new(bar.left(), top), Pos2::new(bar.right(), top)],
            Stroke::new(0.6, theme::BORDER),
        );
        let y = bar.center().y;
        let mut x = rect.left() + 8.0;
        let mut put = |txt: String, col: Color32| {
            let g = painter.layout_no_wrap(txt, FontId::proportional(11.0), col);
            let w = g.size().x;
            painter.galley(Pos2::new(x, y - g.size().y * 0.5), g, col);
            x += w + 14.0;
        };
        put(crate::i18n::text(ui.ctx(),"Selected:").to_string(), theme::TEXT_DIM);
        if let Some((s0, e0)) = app.selection {
            let bp = if selected_bases.is_some() {
                3
            } else {
                hx_render::span_len(s0, e0, len)
            };
            let mut n_feats = 0usize;
            let mut enzymes = 0usize;
            let selection_parts = hx_render::split_wrap(s0, e0, len);
            // 用已算好的特征片段（避免再借 `doc`）
            for f in &feats {
                if !f.orf
                    && f.segs.iter().any(|s| {
                        selection_parts
                            .iter()
                            .any(|&(a, b)| s.start < b && s.end > a)
                    })
                {
                    n_feats += 1;
                }
            }
            for (k, _) in app.sites.selected.iter().enumerate() {
                enzymes += app
                    .sites
                    .hits
                    .get(k)
                    .map(|h| {
                        h.iter()
                            .filter(|p| {
                                selection_parts
                                    .iter()
                                    .any(|&(a, b)| (**p as usize) >= a && (**p as usize) < b)
                            })
                            .count()
                    })
                    .unwrap_or(0);
            }
            let codons = bp / 3;
            put(crate::ui_format!(ui.ctx(),"{} bp", nfmt(bp)), theme::TEXT);
            put(crate::ui_format!(ui.ctx(),"{} codons", nfmt(codons)), theme::TEXT);
            put(crate::ui_format!(ui.ctx(),"{} features", nfmt(n_feats)), theme::TEXT);
            put(crate::ui_format!(ui.ctx(),"{} enzymes", nfmt(enzymes)), theme::TEXT);
            put(
                crate::ui_format!(ui.ctx(),"Insertion Point: {}", nfmt(s0 + 1)),
                theme::TEXT_DIM,
            );
        } else if selected_bases.is_some() {
            put(crate::ui_format!(ui.ctx(),"{} bp",3), theme::TEXT);
            put(crate::i18n::text(ui.ctx(),"1 codon (joined segments)").into(), theme::TEXT);
        } else if let Some(c) = app.cursor_bp {
            put(crate::ui_format!(ui.ctx(),"Insertion Point: {}", nfmt(c + 1)), theme::TEXT_DIM);
        }
        // 固定线宽超宽时的横向滚动条（SnapGene 的横向平移）
        if pan_max > 0.0 {
            let track = Rect::from_min_size(
                Pos2::new(bar.right() - 260.0, bar.center().y - 3.0),
                Vec2::new(250.0, 6.0),
            );
            painter.rect_filled(track, 3.0, theme::with_alpha(theme::BORDER, 140));
            let view_w = track.width() * (avail / content_w).min(1.0);
            let off = track.width() * (app.seq_pan_x / content_w);
            let thumb = Rect::from_min_size(
                Pos2::new(track.left() + off, track.top()),
                Vec2::new(view_w.max(20.0), 6.0),
            );
            painter.rect_filled(thumb, 3.0, theme::with_alpha(theme::TEXT_DIM, 200));
            if let Some(pos) = painter.ctx().pointer_interact_pos() {
                if track.contains(pos)
                    && painter
                        .ctx()
                        .input(|i| i.pointer.primary_down() || i.pointer.primary_clicked())
                {
                    let t = ((pos.x - track.left()) / track.width()).clamp(0.0, 1.0);
                    app.seq_pan_x = (t * content_w - view_w * 0.5).clamp(0.0, pan_max);
                }
            }
        }
    }

    // ---- 底部：全序列概览 + 缩放滑块 ----
    draw_bottom(
        &painter,
        rect,
        len,
        content_left,
        content_right,
        &mut v,
        &rows,
        view_h,
        max_scroll,
        &feats,
        app.show_colors,
        app.seq_zoom,
    );

    // 写回排版状态（自适应 line_bp + 滚动）
    app.seq = v;

    // ---- 悬停工具卡 ----
    // oracle `/tmp/opencode/sg/dn3c.png`、`dn4c.png`：DNA 序列视图的卡片与蛋白同一版式
    // （标题 + 类型芯片 + 区间 `867 .. 1455 · 589 bp` + 虚线 + NOTE），**多一个方向箭头**
    // （正向 `←`?? 实测 ori 是反向特征 ⇒ 箭头指左；蛋白没有方向箭头）。
    if let Some(p) = resp.hover_pos() {
        if let Some((_, source, _)) = feat_hits
            .iter()
            .filter(|(r, _, is_label)| *is_label && r.expand(1.0).contains(p))
            .next()
            .or_else(|| feat_hits.iter().find(|(r, _, _)| r.contains(p)))
        {
            if let Some(f) = feats.iter().find(|f| f.source == *source && !f.orf) {
                let segs: Vec<(usize, usize)> = f.segs.iter().map(|s| (s.start, s.end)).collect();
                let card = feat_card::Card {
                    name: &f.name,
                    kind: &f.kind,
                    sub_type: &f.sub_type,
                    note: &f.note,
                    segs: &segs,
                    prefix: "",
                    unit: "bp",
                    swatch: Some(f.color),
                    direction: f.reverse.then_some(true),
                };
                egui::Tooltip::for_enabled(&resp)
                    .at_pointer()
                    .gap(14.0)
                    .show(|ui| {
                        let frame = feat_card::frame();
                        frame.show(ui, |ui| {
                            ui.style_mut().override_font_id =
                                Some(egui::FontId::proportional(11.5));
                            feat_card::show(ui, &card);
                        });
                    });
            }
        }
    }

    rect
}

/// 底部两条：特征概览（点击/拖拽定位）+ 缩放滑块。
#[allow(clippy::too_many_arguments)]
fn draw_bottom(
    painter: &egui::Painter,
    rect: Rect,
    len: usize,
    content_left: f32,
    content_right: f32,
    v: &mut SeqView,
    rows: &[Row],
    view_h: f32,
    max_scroll: f32,
    feats: &[Feature],
    show_colors: bool,
    _seq_zoom: f32,
) {
    let top = rect.bottom() - BOTTOM_H;
    let bg = Rect::from_min_max(Pos2::new(rect.left(), top), rect.right_bottom());
    painter.rect_filled(bg, 0.0, theme::PANEL_ALT);
    painter.line_segment(
        [Pos2::new(rect.left(), top), Pos2::new(rect.right(), top)],
        Stroke::new(0.8, theme::BORDER),
    );
    // 轨道左侧留出信息文字区（SnapGene 的「碱基对 / 范围」标注区）
    let (ox0, ox1) = (content_left + 130.0, content_right - 24.0);
    let ow = (ox1 - ox0).max(20.0);
    let xb = |bp: usize| -> f32 { ox0 + (bp as f32 / len.max(1) as f32) * ow };

    // 可视范围（bp）
    let cy = v.scroll_y;
    let mut vis_start = len;
    let mut vis_end = len;
    for r in rows {
        let y = r.top - cy;
        if y + r.h > 0.0 && y < view_h {
            vis_start = vis_start.min(r.start);
            vis_end = r.end.max(vis_start + 1);
        }
    }
    // ---- 概览条：全部特征 ----
    let ov = Rect::from_min_size(Pos2::new(ox0, top + 4.0), Vec2::new(ow, 16.0));
    painter.rect_filled(ov, 2.0, theme::SEQ_BG);
    painter.rect_stroke(
        ov,
        2.0,
        Stroke::new(0.7, theme::BORDER),
        egui::StrokeKind::Inside,
    );
    {
        for f in feats.iter().filter(|f| !f.orf) {
            for seg in &f.segs {
                let (a, b) = (seg.start, seg.end);
                if a >= b {
                    continue;
                }
                let x0 = xb(a);
                let x1 = xb(b).max(x0 + 1.0);
                let r =
                    Rect::from_min_size(Pos2::new(x0, ov.top() + 3.0), Vec2::new(x1 - x0, 10.0));
                let fill = if show_colors {
                    f.color
                } else {
                    theme::FEATURE_NO_COLOR
                };
                painter.add(Shape::convex_polygon(
                    feature_polygon(r, seg.left_arrow, seg.right_arrow),
                    theme::mix(fill, Color32::WHITE, 0.25),
                    Stroke::new(0.6, theme::darken(fill, 0.3)),
                ));
            }
        }
    }
    // 可视范围框
    let view_rect = Rect::from_min_max(
        Pos2::new(xb(vis_start), ov.top() - 1.0),
        Pos2::new(xb(vis_end.min(len)), ov.bottom() + 1.0),
    );
    painter.rect_filled(view_rect, 2.0, theme::with_alpha(theme::SELECTION, 36));
    painter.rect_stroke(
        view_rect,
        2.0,
        Stroke::new(1.0, theme::with_alpha(theme::SELECTION, 170)),
        egui::StrokeKind::Inside,
    );

    // 概览条交互：点击/拖拽 → 让该 bp 居中
    if let Some(pos) = painter.ctx().pointer_interact_pos() {
        let hit = Rect::from_min_max(
            Pos2::new(ox0, ov.top() - 2.0),
            Pos2::new(ox1, ov.bottom() + 2.0),
        );
        if hit.contains(pos) {
            let bp = ((pos.x - ox0) / ow * len as f32).clamp(0.0, len as f32 - 1.0) as usize;
            if painter
                .ctx()
                .input(|i| i.pointer.primary_down() || i.pointer.primary_clicked())
            {
                scroll_to_bp(v, rows, bp, view_h, max_scroll);
            }
        }
    }

    // ---- 滑块条：碱基刻度 + 可视范围滑块 + 缩放按钮 ----
    let sl = Rect::from_min_size(Pos2::new(ox0, top + 26.0), Vec2::new(ow, 18.0));
    // 刻度（整齐步长，间隔 ≥ 90 px）
    let px_per_bp = ow / len.max(1) as f32;
    for (bp, major) in ruler_ticks_subdiv(len, px_per_bp, 90.0, 5) {
        let b = bp as usize;
        if b > len {
            break;
        }
        let x = xb(b);
        let h = if major { 6.0 } else { 3.0 };
        painter.line_segment(
            [Pos2::new(x, sl.top()), Pos2::new(x, sl.top() + h)],
            Stroke::new(0.8, theme::TICK_MAJOR),
        );
        if major && b > 0 {
            painter.text(
                Pos2::new(x, sl.top() + h + 1.0),
                Align2::CENTER_TOP,
                nfmt_plain(b),
                FontId::proportional(9.0),
                theme::RULER_TEXT,
            );
        }
    }
    painter.line_segment(
        [
            Pos2::new(ox0, sl.bottom() - 5.0),
            Pos2::new(ox1, sl.bottom() - 5.0),
        ],
        Stroke::new(1.0, theme::with_alpha(theme::BORDER, 200)),
    );
    // 滑块（可视范围）
    let thumb = Rect::from_min_max(
        Pos2::new(xb(vis_start), sl.top() + 5.0),
        Pos2::new(
            xb(vis_end.min(len)).max(xb(vis_start) + 12.0),
            sl.top() + 13.0,
        ),
    );
    painter.rect_filled(thumb, 3.0, theme::PANEL);
    painter.rect_stroke(
        thumb,
        3.0,
        Stroke::new(1.0, theme::with_alpha(theme::TEXT_DIM, 200)),
        egui::StrokeKind::Inside,
    );
    for i in 0..3 {
        let x = thumb.center().x - 3.0 + i as f32 * 3.0;
        painter.line_segment(
            [
                Pos2::new(x, thumb.top() + 2.0),
                Pos2::new(x, thumb.bottom() - 2.0),
            ],
            Stroke::new(0.7, theme::with_alpha(theme::TEXT_DIM, 160)),
        );
    }
    // 左侧信息（长度 + 当前可见范围 + 缩放）
    painter.text(
        Pos2::new(rect.left() + 6.0, top + 30.0),
        Align2::LEFT_TOP,
        format!(
            "{} bp  ·  {} … {}",
            nfmt(len),
            nfmt(vis_start + 1),
            nfmt(vis_end)
        ),
        FontId::proportional(10.0),
        theme::TEXT_DIM,
    );
    // 滑块交互：拖拽 → 平移
    if let Some(pos) = painter.ctx().pointer_interact_pos() {
        let hit = Rect::from_min_max(Pos2::new(ox0, sl.top()), Pos2::new(ox1, sl.bottom()));
        if hit.contains(pos)
            && painter
                .ctx()
                .input(|i| i.pointer.primary_down() || i.pointer.primary_clicked())
        {
            let bp = ((pos.x - ox0) / ow * len as f32).clamp(0.0, len as f32 - 1.0) as usize;
            scroll_to_bp(v, rows, bp, view_h, max_scroll);
        }
    }
}

/// 把视图滚动到让 `bp` 位于窗口中间（供底部概览/滑块定位用）。
fn scroll_to_bp(v: &mut SeqView, rows: &[Row], bp: usize, view_h: f32, max_scroll: f32) {
    let row = (bp / v.line_bp.max(1)).min(rows.len().saturating_sub(1));
    let Some(r) = rows.get(row) else { return };
    let target = r.top + r.h * 0.5 - view_h * 0.5;
    v.scroll_y = target.clamp(0.0, max_scroll);
}

#[cfg(test)]
mod tests {
    #[test]fn duplex_selection_preserves_staggered_wrapping_and_negative_offsets(){
        assert_eq!(super::duplex_selection_parts(212,308,219),vec![(212,219),(0,89)]);
        assert_eq!(super::duplex_selection_parts(216,312,219),vec![(216,219),(0,93)]);
        assert_eq!(super::duplex_selection_parts(-3,90,219),vec![(216,219),(0,90)]);
        assert_eq!(super::duplex_selection_parts(89,212,219),vec![(89,212)]);
        assert_eq!(super::duplex_selection_parts(93,216,219),vec![(93,216)]);
        assert_eq!(super::duplex_selection_parts(4,223,219),vec![(0,219)]);
        assert!(super::duplex_selection_parts(0,0,0).is_empty());
    }
    use super::*;

    #[test]
    fn forward_coding_map_is_genomic_order() {
        let seq = b"AAACCCGGGTTT";
        let (map, code) = coding_map(seq, &[(3, 6)], true);
        assert_eq!(map, vec![3, 4, 5]);
        assert_eq!(code, b"CCC");
    }

    #[test]
    fn reverse_coding_map_is_reversed_complement() {
        let seq = b"AAACCCGGGTTT";
        let (map, code) = coding_map(seq, &[(3, 9)], false);
        assert_eq!(map, vec![8, 7, 6, 5, 4, 3]);
        assert_eq!(code, b"CCCGGG");
    }

    #[test]
    fn cloned_product_extra_sites_respect_snapgene_hybridization_thresholds() {
        let Ok(bytes)=std::fs::read("/tmp/opencode/oracle/Cloned.dna") else {
            eprintln!("skip: saved SnapGene cloning fixture /tmp/opencode/oracle/Cloned.dna is not present");
            return;
        };
        let doc=crate::doc::Document::from_bytes(bytes).unwrap();
        assert_eq!(doc.len(),4581);
        let params=hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap());
        assert_eq!(params.min_melting_temperature,40.0);
        assert_eq!(doc.primers().len(),2);
        assert_eq!(doc.primers()[0].name,"Primer 1");
        assert_eq!(doc.primers()[1].name,"Primer 2");
        let actual=doc.primers().iter().flat_map(|p|p.sites.iter().map(|s|s.melting_temperature.unwrap())).collect::<Vec<_>>();
        assert_eq!(actual,vec![50.0,26.0,26.0,31.0,53.0]);
        assert_eq!(actual.iter().filter(|tm|**tm>=params.min_melting_temperature).count(),2);
    }

    #[test]
    fn reverse_cds_translates_cleanly() {
        // 反向链上的 CDS：基因组（上链）= ATGGCTTAA 的反向互补 = TTAAGCCAT
        let seq = b"TTAAGCCAT";
        let (map, code) = coding_map(seq, &[(0, 9)], false);
        assert_eq!(code, b"ATGGCTTAA");
        assert_eq!(map[0], 8);
        let aa = hx_algo::translate::translate_orf(&code, false);
        assert_eq!(aa, vec!["Met", "Ala"]);
    }

    fn annotation(attrs: &str, segments: &str) -> hx_formats::annot::FeatureAnnotation {
        hx_formats::annot::parse_features(
            format!("<Features><Feature {attrs}>{segments}</Feature></Features>").as_bytes(),
        )
        .remove(0)
    }

    #[test]
    fn feature_translation_uses_only_translated_segments_and_explicit_start_setting() {
        let f = annotation(
            "name=\"6xHis\" directionality=\"1\"",
            r#"<Segment range="1-3" translated="1"/><Segment range="4-6"/><Segment range="7-9" translated="1"/>"#,
        );
        let c = feature_translation(b"CATGCTCAT", &f, false, &[]).unwrap();
        assert_eq!(c.aa, vec!["His", "His"]);
        assert_eq!(c.map, vec![0, 1, 2, 6, 7, 8]);
        let f = annotation(
            "name=\"lacI\" directionality=\"1\" translateFirstCodonAsMet=\"1\" readingFrame=\"2\"",
            r#"<Segment range="1-7" translated="1"/>"#,
        );
        let c = feature_translation(b"AGTGAAA", &f, false, &[]).unwrap();
        assert_eq!(c.aa, vec!["Met", "Lys"]);
        assert_eq!(c.codon_bases(0), Some([1, 2, 3]));
    }

    #[test]
    fn stop_codons_are_retained_and_adjacent_stop_does_not_extend_feature_bar() {
        let f = annotation(
            "name=\"CDS\" directionality=\"1\"",
            r#"<Segment range="1-6" translated="1"/>"#,
        );
        let c = feature_translation(b"ATGGCTTAA", &f, false, &[]).unwrap();
        assert_eq!(c.aa, vec!["Met", "Ala", "***"]);
        assert_eq!(f.segments[0].range.end, 6);
        let reverse = annotation(
            "directionality=\"2\"",
            r#"<Segment range="4-9" translated="1"/>"#,
        );
        let c = feature_translation(b"TTAAGCCAT", &reverse, false, &[]).unwrap();
        assert_eq!(c.aa, vec!["Met", "Ala", "***"]);
        assert_eq!(c.codon_bases(2), Some([2, 1, 0]));
    }

    #[test]
    fn reverse_asp_selects_exactly_atc_on_the_top_strand() {
        let f = annotation(
            "name=\"KanR\" directionality=\"2\"",
            r#"<Segment range="1-6" translated="1"/>"#,
        );
        let c = feature_translation(b"ATCCAT", &f, false, &[]).unwrap();
        assert_eq!(c.aa, vec!["Met", "Asp"]);
        let bases = c.codon_bases(1).unwrap();
        assert_eq!(bases, [2, 1, 0]);
        assert_eq!(codon_selection(bases, 6), Some((0, 3)));
        assert_eq!(codon_parts(bases), vec![(0, 3)]);
    }

    #[test]
    fn fusion_linker_uses_run_on_limit_and_stops_before_next_in_frame_feature() {
        let xml = br#"<Features><Feature name="tag" directionality="1" maxRunOn="12"><Segment range="1-3" translated="1"/></Feature><Feature name="next CDS" directionality="1"><Segment range="7-9" translated="1"/></Feature></Features>"#;
        let features = hx_formats::annot::parse_features(xml);
        let c = feature_translation(b"CATAGCGCTTAA", &features[0], false, &features).unwrap();
        assert_eq!(c.aa, vec!["His", "Ser"]);
        assert_eq!(c.annotated_codons, 1);
        assert_eq!(c.map, (0..6).collect::<Vec<_>>());
        let alone =
            feature_translation(b"CATAGCGCTTAA", &features[0], false, &features[..1]).unwrap();
        assert_eq!(alone.aa, vec!["His", "Ser", "Ala", "***"]);
    }

    #[test]
    fn split_codons_keep_all_three_bases_without_selecting_a_gap() {
        let c = Cds {
            forward: true,
            map: vec![9, 10, 11],
            aa: vec!["Asp"],
            annotated_codons: 1,
        };
        assert_eq!(c.codons_in_row(0, 10)[0].end, 10);
        assert_eq!(c.codons_in_row(10, 20)[0].start, 10);
        assert_eq!(c.codons_in_row(10, 20)[0].visual_slot, 1);
        assert_eq!(
            codon_selection(c.codon_bases(0).unwrap(), 20),
            Some((9, 12))
        );
        assert_eq!(codon_selection([99, 0, 1], 100), Some((99, 2)));
        assert_eq!(codon_selection([98, 99, 0], 100), Some((98, 1)));
        assert_eq!(codon_parts([9, 30, 31]), vec![(9, 10), (30, 32)]);
        assert_eq!(codon_selection([9, 30, 31], 100), None);
        let wrap = Cds {
            forward: false,
            map: vec![0, 99, 98],
            aa: vec!["Asp"],
            annotated_codons: 1,
        };
        let fragments = wrap.codons_in_row(0, 100);
        assert_eq!(fragments.len(), 2);
        assert_eq!((fragments[0].start, fragments[0].end), (0, 1));
        assert_eq!((fragments[1].start, fragments[1].end), (98, 100));
    }

    fn feature(a: usize, b: usize, cds: bool) -> Feature {
        Feature {
            source: a,
            name: "F".into(),
            kind: "misc_feature".into(),
            sub_type: String::new(),
            note: String::new(),
            reverse: false,
            segs: vec![FeatureSegment {
                start: a,
                end: b,
                left_arrow: false,
                right_arrow: false,
            }],
            color: Color32::GREEN,
            cds: cds.then(|| Cds {
                forward: true,
                map: (a..b).collect(),
                aa: vec!["Ala"; (b - a) / 3],
                annotated_codons: (b - a) / 3,
            }),
            orf: false,
        }
    }

    #[test]
    fn adjacent_unsorted_features_share_a_track_and_duplicates_keep_their_own_translation() {
        let v = SeqView {
            line_bp: 160,
            ..SeqView::default()
        };
        // 与实际文件顺序相同：中间的RBS先于lac operator。
        let features = vec![
            feature(24, 43, false),
            feature(82, 105, false),
            feature(43, 68, false),
        ];
        let (tracks, levels) = feature_tracks(&features, 0, 160, &v, |_| 8.0);
        assert_eq!(levels.len(), 1);
        assert_eq!(
            tracks.iter().map(|t| t.feature).collect::<Vec<_>>(),
            vec![0, 2, 1]
        );
        let features = (0..6).map(|_| feature(20, 80, true)).collect::<Vec<_>>();
        let (tracks, levels) = feature_tracks(&features, 0, 160, &v, |_| 8.0);
        assert_eq!(levels.len(), 6); // 不再静默丢掉第5层以后的特征。
        assert_eq!(tracks.len(), 6);
        assert!(tracks.iter().all(|t| t.codons.len() == 20));
        assert!(levels.iter().all(|l| l.translation));
    }

    #[test]
    fn orfs_use_separate_compact_tracks_even_when_they_overlap_a_cds() {
        let mut features = vec![feature(0, 90, true)];
        for (a, b, forward) in [
            (0, 90, true),
            (3, 87, true),
            (90, 120, true),
            (0, 90, false),
        ] {
            let mut f = feature(a, b, true);
            f.orf = true;
            f.cds.as_mut().unwrap().forward = forward;
            features.push(f);
        }
        let v = SeqView {
            line_bp: 120,
            ..SeqView::default()
        };
        let (tracks, levels) = feature_tracks(&features, 0, 120, &v, |_| 8.0);
        assert_eq!((tracks.len(), levels.len()), (1, 1));
        let (orfs, count) = orf_tracks(&features, 0, 120);
        assert_eq!((orfs.len(), count), (4, 3));
        let lane = |fi| orfs.iter().find(|t| t.feature == fi).unwrap().level;
        assert_eq!(lane(1), lane(3)); // 相邻 ORF 共享一行。
        assert_ne!(lane(1), lane(2)); // 重叠框分行。
        assert_eq!(lane(4), 2); // 反链在正链之后。
        assert!(orf_tracks(&features, 120, 160).0.is_empty());
    }

    #[test]
    fn orf_stops_and_click_mapping_follow_the_coding_strand_including_origin_wrap() {
        for (seq, circular, start, end, strand, stop_bases) in [
            (b"ATGGCTTAA".as_slice(), false, 0, 6, 1, [6, 7, 8]),
            (b"TTAAGCCAT".as_slice(), false, 3, 9, -1, [2, 1, 0]),
            (b"TAAATGGCT".as_slice(), true, 3, 0, 1, [0, 1, 2]),
            (b"AGCCATTTA".as_slice(), true, 0, 6, -1, [8, 7, 6]),
        ] {
            let o = hx_algo::orf::Orf {
                start,
                end,
                strand,
                has_start: true,
            };
            let cds = orf_translation(seq, &o, circular);
            assert_eq!(cds.aa, vec!["Met", "Ala", "***"]);
            assert_eq!(cds.annotated_codons, 2);
            assert_eq!(cds.codon_bases(2), Some(stop_bases));
            assert!(codon_selection(cds.codon_bases(1).unwrap(), seq.len()).is_some());
        }
        let o = hx_algo::orf::Orf {
            start: 0,
            end: 6,
            strand: 1,
            has_start: true,
        };
        assert_eq!(orf_translation(b"ATGGCT", &o, false).aa, vec!["Met", "Ala"]);
        assert_eq!(
            orf_translation(b"ATGGCTGCT", &o, false).aa,
            vec!["Met", "Ala"]
        );
    }

    #[test]
    fn orf_colors_distinguish_strands_and_respect_the_color_toggle() {
        assert_ne!(orf_fill(true, true), orf_fill(false, true));
        assert_eq!(orf_fill(true, false), orf_fill(false, false));
        let gray = orf_fill(true, false);
        assert_eq!(gray.r(), gray.g());
        assert_eq!(gray.g(), gray.b());
    }

    #[test]
    fn feature_tips_stay_inside_the_actual_span() {
        let r = Rect::from_min_max(Pos2::new(10.0, 0.0), Pos2::new(30.0, 13.0));
        assert_eq!(feature_polygon(r, false, false).len(), 4);
        for (left, right) in [(true, false), (false, true), (true, true)] {
            assert!(feature_polygon(r, left, right)
                .iter()
                .all(|p| r.contains(*p)));
        }
    }

    fn empty_row(fwd_levels: usize, fwd_h: f32) -> Row {
        Row {
            start: 0,
            end: 60,
            top: 100.0,
            h: 200.0,
            tracks: Vec::new(),
            levels: Vec::new(),
            features_h: 0.0,
            orfs: Vec::new(),
            orfs_h: 0.0,
            fwd_levels,
            rev_levels: 1,
            fwd_h,
            primers: Vec::new(),
        }
    }

    fn primer_with(tail: &str, forward: bool, tail_left: bool) -> PrimerRow {
        PrimerRow {
            name: "P".into(),
            segs: vec![(10, 30)],
            forward,
            site_start: 10,
            site_end: 30,
            tail: tail.into(),
            tail_left,
            phosphorylated: false,
        }
    }

    const BW: f32 = 10.8;
    const PH: f32 = 16.0;

    /// SnapGene 的台阶：尾巴框与退火框同高、朝外侧错开 `step`、横向压进 `step`。
    #[test]
    fn tail_box_steps_out_and_overlaps_the_annealed_box() {
        let row = empty_row(1, 20.0);
        let span = PrimerSpan {
            primer: 0,
            x0: 0.0,
            x1: 20.0 * BW,
            level: 0,
            tip: true,
            five_prime: false,
            tail_here: true,
        };
        let step = PH * 0.45;
        let lap = step * 0.4;
        for forward in [true, false] {
            for tail_left in [true, false] {
                let p = primer_with("TTTTTTTTT", forward, tail_left);
                let g = primer_boxes(&span, &p, &row, 100.0, 50.0, BW, PH, 3.0, 30.0);
                let t = g.trect.expect("有尾巴就该有尾巴框");
                assert!((t.height() - PH).abs() < 0.01);
                // 框宽 = 尾巴碱基的网格宽 + 压进退火框的 `lap`（实测 ≈ 0.4 × step）
                assert!((t.width() - (9.0 * BW + lap)).abs() < 0.01);
                if forward {
                    assert!((g.arect.top() - t.top() - step).abs() < 0.01);
                } else {
                    assert!((t.top() - g.arect.top() - step).abs() < 0.01);
                }
                if tail_left {
                    // 左边缘 = 退火框左侧再往左 9 格
                    assert!((t.right() - (g.arect.left() + lap)).abs() < 0.01);
                    assert!((t.left() - (g.arect.left() - 9.0 * BW)).abs() < 0.01);
                } else {
                    assert!((t.left() - (g.arect.right() - lap)).abs() < 0.01);
                    assert!((t.right() - (g.arect.right() + 9.0 * BW)).abs() < 0.01);
                }
            }
        }
    }

    /// 尾巴字母与模板字母同列：左尾落在退火框左侧紧邻的 `n` 格，右尾从
    /// 退火框右边缘起排（SnapGene 的 `CCC` 尾巴实测与模板逐像素同列）。
    #[test]
    fn tail_text_sits_on_the_template_grid() {
        let x0 = 300.0;
        let x1 = 300.0 + 20.0 * BW;
        // 左尾：最后一格紧贴退火框左侧，与「上一行同一位置」的碱基同列
        let left = tail_text_left(x0, x1, 3, BW, true);
        assert!((left + 3.0 * BW - x0).abs() < 0.01);
        assert!((left + 2.5 * BW - (x0 - 0.5 * BW)).abs() < 0.01);
        // 右尾：第一格从退火框右边缘开始
        let right = tail_text_left(x0, x1, 3, BW, false);
        assert!((right - x1).abs() < 0.01);
        assert!((right + 0.5 * BW - (x1 + 0.5 * BW)).abs() < 0.01);
    }

    /// 两框的并集轮廓：8 个顶点、面积 = 两框之和 − 角上那块重叠，
    /// 且没有一条竖边落在另一个框的内部（否则会压住末位碱基）。
    #[test]
    fn union_outline_drops_the_internal_edges() {
        let a = Rect::from_min_max(Pos2::new(0.0, 100.0), Pos2::new(216.0, 116.0));
        let step = 7.2;
        let tw = 97.2;
        let area = |v: &[Pos2]| {
            let mut s = 0.0;
            for i in 0..v.len() {
                let p = v[i];
                let q = v[(i + 1) % v.len()];
                s += p.x * q.y - q.x * p.y;
            }
            (s * 0.5).abs()
        };
        for (tail_left, forward) in [(true, true), (false, true), (true, false), (false, false)] {
            let tx0 = if tail_left {
                a.left() + step - tw
            } else {
                a.right() - step
            };
            let ty0 = if forward {
                a.top() - step
            } else {
                a.top() + step
            };
            let t = Rect::from_min_max(Pos2::new(tx0, ty0), Pos2::new(tx0 + tw, ty0 + PH));
            let v = primer_union_outline(a, t, tail_left, forward);
            assert_eq!(v.len(), 8);
            let expect = a.area() + t.area() - step * (PH - step);
            assert!(
                (area(&v) - expect).abs() < 0.5,
                "并集面积 {} != {expect}（tail_left={tail_left} forward={forward}）",
                area(&v)
            );
            for w in v.windows(2) {
                let (p, q) = (w[0], w[1]);
                if (p.x - q.x).abs() >= 0.01 {
                    continue;
                }
                let (y0, y1) = (p.y.min(q.y), p.y.max(q.y));
                let inside_a = p.x > a.left() + 0.01
                    && p.x < a.right() - 0.01
                    && y0 > a.top() + 0.01
                    && y1 < a.bottom() - 0.01;
                let inside_t = p.x > t.left() + 0.01
                    && p.x < t.right() - 0.01
                    && y0 > t.top() + 0.01
                    && y1 < t.bottom() - 0.01;
                assert!(!inside_a && !inside_t, "内部竖边 {p:?}-{q:?}");
            }
        }
    }
}
fn sequence_color(base: u8, mode: u8) -> Color32 {
    let b = base.to_ascii_uppercase();
    match mode {
        1 => if matches!(b, b'G' | b'C') { Color32::from_rgb(25, 115, 190) } else { theme::TEXT },
        2 => if matches!(b, b'A' | b'G') { Color32::from_rgb(140, 55, 165) } else if matches!(b, b'C' | b'T' | b'U') { Color32::from_rgb(20, 135, 115) } else { theme::TEXT },
        3 => match b { b'A' => Color32::from_rgb(30, 140, 65), b'C' => Color32::from_rgb(40, 95, 200), b'G' => Color32::from_rgb(170, 110, 15), b'T' | b'U' => Color32::from_rgb(195, 45, 60), _ => theme::TEXT },
        _ => theme::TEXT,
    }
}
