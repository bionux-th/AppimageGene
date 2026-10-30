//! 图谱视图（对齐 SnapGene 的绘制习惯）：
//!
//! - **环状**：双线 DNA 环 → 刻度与 ORF 独立轨道 → 留白后分层特征弧带（极小特征画成
//!   小方块标记）→ 特征名沿弧旋转写在环内（放得下就写在带里）→ 径向刻度与旋转
//!   数字 → 环外酶位点标签（`Name (pos)`，左侧坐标在前；唯一者加粗，细灰引线）
//!   → 圆心名称 + `N bp`。
//! - **线性**：双线轴与刻度 → 独立 ORF 细箭头区 → 按重叠区间分行的特征带。
//!
//! 交互：环状图谱固定方向、沿 DNA 环拖选；线性拖拽 = 平移，滚轮 = 缩放，单击特征/ORF = 选区
//! （点空白取消）。环状图谱的原点只由明确的“居中特征/设置起点”操作改变。

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Vec2};

use hx_render::map::{ruler_ticks_subdiv, MapView};
use hx_render::circular_labels::{CircularLabel, CircularLabelLayout};
use hx_render::span_len;

use crate::app::{feat_group, nfmt, nfmt_plain, FeatGroup, HelixApp};
use crate::theme;

/// 特征带厚度下限（半径很小时也看得见）。
const BAND_MIN: f32 = 9.0;
/// 小于该弧长（像素）的特征画成「小标记 + 外置标签」。
const TINY_PX: f32 = 13.0;
/// 线性视图里「条已经看不见」的宽度阈值（像素）。低于它的特征不排条下标签，
/// 直接走轴上方：DNA/RNA 画小胶囊、蛋白画细刻度 + 纯文字（oracle 实测：
/// DNA 的 19 bp 启动子 2 px 是胶囊，蛋白的 2 aa Site 5 px 是细刻度，
/// 而 8.6 px 的 MCS / 7 px 的 3 aa 螺旋都是「条 + 条下标签」）。
const NO_BAR_PX: f32 = 6.0;
/// 箭头长度（像素，弧长）。
const TIP_PX: f32 = 9.0;
/// 外置标签层间距。
const LEVEL_DY: f32 = 16.0;
/// ORF 箭头颜色（正向 / 反向）——SnapGene 用细箭头把 ORF 与特征区分开。
const ORF_FWD: Color32 = Color32::from_rgb(0xd9, 0x8a, 0x2b);
const ORF_REV: Color32 = Color32::from_rgb(0x7f, 0xb0, 0x3a);
/// 环状骨架的线宽。oracle 实测（`sg-rna-ring.png`）：SnapGene 的黑色环线 **3 px**。
const RING_LINE_W: f32 = 2.6;
const LINEAR_ORF_TOP: f32 = 26.0;
const LINEAR_ORF_STEP: f32 = 14.0;

/// 与外置特征标签使用相同的字体和字号。
fn map_font() -> FontId {
    FontId::proportional(11.5)
}

/// 线性图谱「条下标签行」中的一项：写不进特征条里的名字排在条下方的行里，
/// 由一条细引线连回特征条（SnapGene 的线性图谱排布）。
struct BelowLabel {
    /// 特征所在车道（`interval_lanes` 的行号）。
    lane: usize,
    /// 车道内的标签行号（0 最靠上）。
    row: usize,
    /// 文本中心 x（避让后）。
    xc: f32,
    /// 特征条中心 x（引线起点）。
    cx: f32,
    /// 文本宽度。
    tw: f32,
    /// 特征下标。
    idx: usize,
    name: String,
}

#[derive(Clone, Copy)]
struct RingTracks {
    feature_radius: f32,
    feature_step: f32,
    band: f32,
    orf_radius: f32,
    orf_step: f32,
}

impl RingTracks {
    fn new(radius: f32, band: f32, features: usize, orfs: usize) -> Self {
        let orf_depth = if orfs == 0 {
            0.0
        } else {
            (orfs - 1) as f32 * 12.0 + 16.0
        };
        let feature_depth = 30.0 + orf_depth + band * 0.46;
        let depth = feature_depth + features.saturating_sub(1) as f32 * (band + 24.0) + band * 0.46;
        // 窄窗口或多层注释时一起缩放留白和带宽，不把内层压到圆心。
        let scale = ((radius - 28.0).max(1.0) / depth).min(1.0);
        Self {
            feature_radius: radius - feature_depth * scale,
            feature_step: (band + 24.0) * scale,
            band: band * scale,
            orf_radius: radius - 30.0 * scale,
            orf_step: 12.0 * scale,
        }
    }
}

/// 跟踪连续角位移，使跨原点、反向拖选和超过半圈的选择都保持用户实际拖动的路径。
#[derive(Clone, Copy)]
struct RingDrag {
    anchor: f64,
    previous: f64,
    sweep: f64,
}

impl RingDrag {
    fn new(bp: f64) -> Self {
        Self {
            anchor: bp.round(),
            previous: bp,
            sweep: 0.0,
        }
    }

    fn update(&mut self, bp: f64, len: usize) -> Option<(usize, usize)> {
        let n = len as f64;
        let mut delta = bp - self.previous;
        if delta > n * 0.5 {
            delta -= n;
        }
        if delta < -n * 0.5 {
            delta += n;
        }
        self.sweep += delta;
        self.previous = bp;
        if self.sweep.abs().round() >= n {
            return Some((0, len));
        }
        let end = (self.anchor + self.sweep).round() as i64;
        let anchor = self.anchor as i64;
        let (lo, hi) = (anchor.min(end), anchor.max(end));
        if lo == hi {
            return None;
        }
        let s = lo.rem_euclid(len as i64) as usize;
        let e = hi.rem_euclid(len as i64) as usize;
        Some((s, if e == 0 { len } else { e }))
    }
}

/// 覆盖一个多段特征的最短区间；环状时保留跨原点语义。
fn feature_range(spans: &[(usize, usize)], len: usize, circular: bool) -> Option<(usize, usize)> {
    let mut spans: Vec<_> = spans.iter().copied().filter(|(a, b)| a < b).collect();
    spans.sort_unstable();
    let start = spans.first()?.0;
    let end = spans.iter().map(|s| s.1).max()?;
    if !circular {
        return Some((start, end));
    }
    let mut best = (start, end);
    let mut gap = len - end + start;
    let mut right = spans[0].1;
    for &(a, b) in spans.iter().skip(1) {
        if a > right && a - right > gap {
            gap = a - right;
            best = (a, right);
        }
        right = right.max(b);
    }
    Some(best)
}

/// 扇环使用独立四边形组成网格，避免将非凸环带当作凸多边形填充产生蓝色楔形。
fn selection_mesh(c: Pos2, a0: f32, sweep: f32, inner: f32, outer: f32) -> egui::Mesh {
    let steps = ((sweep * outer / 4.0).ceil() as usize).clamp(2, 1000);
    let mut mesh = egui::Mesh::default();
    let color = theme::with_alpha(theme::SELECTION, 80);
    for k in 0..=steps {
        let a = a0 + sweep * k as f32 / steps as f32;
        mesh.colored_vertex(pt(c, a, inner), color);
        mesh.colored_vertex(pt(c, a, outer), color);
        if k > 0 {
            let i = (k * 2) as u32;
            mesh.add_triangle(i - 2, i - 1, i);
            mesh.add_triangle(i - 1, i + 1, i);
        }
    }
    mesh
}

/// 将完整注释当作一个实体排布；多段/跨原点的片段始终共用一行。
/// 相邻区间可以共行，但相同坐标的重复注释必须进入不同的行。
fn interval_lanes(spans: &[Vec<(usize, usize)>]) -> (Vec<usize>, usize) {
    let mut order: Vec<usize> = (0..spans.len()).collect();
    order.sort_by_key(|&i| (spans[i].iter().map(|s| s.0).min().unwrap_or(usize::MAX), i));
    let mut occupied: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut levels = vec![0; spans.len()];
    for i in order {
        if !spans[i].iter().any(|(a, b)| a < b) {
            continue;
        }
        let level = occupied
            .iter()
            .position(|row| {
                spans[i]
                    .iter()
                    .all(|&(a, b)| row.iter().all(|&(c, d)| a >= b || a >= d || c >= b))
            })
            .unwrap_or(occupied.len());
        if level == occupied.len() {
            occupied.push(Vec::new());
        }
        occupied[level].extend(spans[i].iter().copied().filter(|(a, b)| a < b));
        levels[i] = level;
    }
    (levels, occupied.len())
}

/// 正向在上、反向在下；每条链内部按实际区间避让，不把阅读框编号当成屏幕行号。
fn linear_orf_lanes(orfs: &[(usize, usize, i8, bool, u8)], len: usize) -> (Vec<usize>, usize) {
    let mut levels = vec![0; orfs.len()];
    let mut count = 0;
    for strand in [1, -1] {
        let indices: Vec<usize> = (0..orfs.len()).filter(|&i| orfs[i].2 == strand).collect();
        let spans: Vec<_> = indices
            .iter()
            .map(|&i| hx_render::split_wrap(orfs[i].0, orfs[i].1, len))
            .collect();
        let (rows, n) = interval_lanes(&spans);
        for (i, row) in indices.into_iter().zip(rows) {
            levels[i] = count + row;
        }
        count += n;
    }
    (levels, count)
}

/// 返回各片段及其是否包含真实的翻译终点，避免跨原点或画布裁剪处出现伪箭头。
fn linear_orf_segments(
    start: usize,
    end: usize,
    strand: i8,
    len: usize,
) -> Vec<(usize, usize, bool)> {
    hx_render::split_wrap(start, end, len)
        .into_iter()
        .map(|(a, b)| (a, b, if strand > 0 { b == end } else { a == start }))
        .collect()
}

/// 沿弧/沿轴画一条「细箭头」（ORF 用）：返回折线点与箭头三角。
fn orf_arrow(
    view: &MapView,
    len: usize,
    a: usize,
    b: usize,
    forward: bool,
    r: f32,
) -> (Vec<Pos2>, Vec<Pos2>) {
    let (s, e) = if b >= a { (a, b) } else { (a, b) };
    let span = if e > s {
        e - s
    } else {
        len.saturating_sub(s) + e
    };
    let steps = (span / 12).clamp(2, 120);
    let pt_at = |bp: f64, r: f32| -> Pos2 {
        let a = view.bp_angle(len, bp);
        Pos2::new(view.center[0] + a.sin() * r, view.center[1] - a.cos() * r)
    };
    let mut pts = Vec::new();
    for i in 0..=steps {
        let f = i as f64 / steps as f64;
        pts.push(pt_at(s as f64 + f * span as f64, r));
    }
    // 箭头三角（3′ 端）
    let tip_i = if forward { pts.len() - 1 } else { 0 };
    let tip = pts[tip_i];
    let prev = if forward {
        pts[pts.len().saturating_sub(2)]
    } else {
        pts[1.min(pts.len() - 1)]
    };
    let d = (tip - prev).normalized();
    let n = Vec2::new(-d.y, d.x);
    let head = vec![tip, tip - d * 7.0 + n * 3.4, tip - d * 7.0 - n * 3.4];
    (pts, head)
}

fn pt(c: Pos2, a: f32, r: f32) -> Pos2 {
    Pos2::new(c.x + a.sin() * r, c.y - a.cos() * r)
}

/// 环上标签的旋转角（保持可读：下半圈翻转 180°）。
fn readable_angle(a: f32) -> f32 {
    if a.cos() < 0.0 {
        a + std::f32::consts::PI
    } else {
        a
    }
}

/// 旋转文字（绕 `p` 处 `align` 指定的锚点）。
fn text_rot(
    painter: &egui::Painter,
    text: &str,
    p: Pos2,
    angle: f32,
    align: Align2,
    font: FontId,
    color: Color32,
    halo: bool,
) {
    let galley = painter.layout_no_wrap(text.to_string(), font, color);
    let size = galley.size();
    let a0 = align
        .pos_in_rect(&egui::Rect::from_min_size(Pos2::ZERO, size))
        .to_vec2();
    let pos = p - a0;
    if halo {
        for d in [
            Vec2::new(-1.1, 0.0),
            Vec2::new(1.1, 0.0),
            Vec2::new(0.0, -1.1),
            Vec2::new(0.0, 1.1),
        ] {
            painter.add(
                egui::epaint::TextShape::new(pos + d, galley.clone(), theme::HALO)
                    .with_override_text_color(theme::HALO)
                    .with_angle_and_anchor(angle, align),
            );
        }
    }
    painter
        .add(egui::epaint::TextShape::new(pos, galley, color).with_angle_and_anchor(angle, align));
}

/// 分片快照（避免绘制期间与 `&mut app` 借用冲突）。
struct FeatRow {
    idx: usize,
    name: String,
    /// SnapGene 的特征类型原文（`type="Site"` / `CDS` / `region`…）。
    kind: String,
    segs: Vec<(usize, usize)>,
    forward_arrow: bool,
    reverse_arrow: bool,
    color: Color32,
}

struct PrimerRow {
    /// `doc.primers()` 下标（点击选中 / 双击编辑用）。
    idx: usize,
    name: String,
    segs: Vec<(usize, usize)>,
    forward: bool,
}

struct Snap {
    len: usize,
    /// 分子拓扑（来自文件）。
    circular: bool,
    /// 是否按**环状布局**绘制（[`hx_render::MapLayout`] 解析结果；线性分子也可能画环）。
    ring: bool,
    name: String,
    /// 长度单位（DNA `bp` / RNA `bases` / 蛋白 `aa`）。
    unit: &'static str,
    /// DNA 文档（决定是否显示酶集标签）。
    is_dna: bool,
    /// 单链分子（RNA、单链 DNA、蛋白）。骨架画**一根**线；双链 DNA 画两根平行线
    /// （用户 2026-10-18 反馈的图一：「环不论单链双链都是双环，直线也是双直线」）。
    single_stranded: bool,
    /// 核苷酸分子（DNA/RNA）。非核苷酸（蛋白）的线性图谱更紧凑：
    /// 特征行间距小，标签优先写在条内（SnapGene 的蛋白图谱就是这样）。
    nucleotide: bool,
    feats: Vec<FeatRow>,
    primers: Vec<PrimerRow>,
    /// (酶名, 位点起点 0-based, 是否唯一, 是否高亮)
    sites: Vec<(String, usize, bool, bool)>,
    /// ORF（`show_translation` 打开时计算）：(起点, 终点, 链, 有起始密码子)
    orfs: Vec<(usize, usize, i8, bool, u8)>,
}

enum RingLabelKind { Feature(usize), Site(usize), Primer(usize) }

struct RingLabels {
    geometry: CircularLabelLayout,
    inputs: Vec<CircularLabel>,
    kinds: Vec<RingLabelKind>,
}

fn primer_label(primer: &PrimerRow, len: usize) -> String {
    let (a,b)=feature_range(&primer.segs,len,true).unwrap_or((0,0));
    format!("{} ({} .. {})",primer.name,nfmt(a+1),nfmt(b))
}

/// Every external kind is measured and packed together. Feature labels which
/// don't fit their actual arc (or collide with another interior label) migrate
/// outside, instead of being painted over neighbouring features.
fn ring_labels(snap:&Snap,p:&egui::Painter,bounds:Rect,rotation:f32,zoom:f32,
    lanes:&[usize],feature_count:usize,orf_count:usize)->RingLabels {
    let center=bounds.center();
    let b=[bounds.left(),bounds.top(),bounds.right(),bounds.bottom()];
    let mut external=vec![false;snap.feats.len()];
    let mut radius=bounds.width().min(bounds.height())*0.40;
    let mut result=None;
    for _ in 0..5 {
        let band=(radius*0.04).clamp(BAND_MIN,18.0);
        let tracks=RingTracks::new(radius,band,feature_count,orf_count);
        let title_size=p.layout_no_wrap(truncate(&snap.name,36),FontId::proportional(15.0),theme::TEXT).size();
        let mut occupied=vec![Rect::from_center_size(center,Vec2::new(title_size.x+16.0,48.0))];
        for (fi,f) in snap.feats.iter().enumerate() {
            if external[fi] {continue;}
            let Some(&(a,b))=f.segs.iter().max_by_key(|(a,b)|b-a) else{continue};
            let r=tracks.feature_radius-lanes[fi] as f32*tracks.feature_step;
            let angle=rotation+(a as f32+(b-a) as f32*0.5)/snap.len as f32*std::f32::consts::TAU;
            let size=p.layout_no_wrap(f.name.clone(),map_font(),theme::TEXT).size()+Vec2::new(8.0,4.0);
            let arc=(b-a) as f32/snap.len as f32*std::f32::consts::TAU*r;
            let text_r=if tracks.band>=12.0 && size.x<=arc {r}else{r-tracks.band-10.0};
            let at=pt(center,angle,text_r);
            let (s,c)=(angle.sin().abs(),angle.cos().abs());
            let box_=Rect::from_center_size(at,Vec2::new(size.x*c+size.y*s,size.x*s+size.y*c));
            if arc<TINY_PX || size.x>r*0.65 || occupied.iter().any(|q|q.intersects(box_)) {external[fi]=true;}
            else {occupied.push(box_);}
        }
        let mut inputs=Vec::new();let mut kinds=Vec::new();
        let internal_only=p.ctx().data(|d|d.get_temp::<crate::preferences::Preferences>(egui::Id::unique("enzyme_render_preferences"))).map(|prefs|prefs.map_internal_only[if snap.is_dna{0}else if snap.unit=="aa"{2}else{1}]).unwrap_or(false);
        for (fi,f) in snap.feats.iter().enumerate().filter(|(i,_)|external[*i]&&!internal_only) {
            let Some(&(a,b))=f.segs.iter().max_by_key(|(a,b)|b-a) else{continue};
            let size=p.layout_no_wrap(f.name.clone(),map_font(),theme::TEXT).size()+Vec2::new(10.0,2.0);
            inputs.push(CircularLabel{angle:rotation+(a as f32+(b-a) as f32*0.5)/snap.len as f32*std::f32::consts::TAU,width:size.x,height:size.y});
            kinds.push(RingLabelKind::Feature(fi));
        }
        for (i,(name,pos,_,_)) in snap.sites.iter().enumerate() {
            let size=p.layout_no_wrap(format!("({}) {}",nfmt(*pos),name),map_font(),theme::TEXT).size()+Vec2::new(4.0,0.0);
            inputs.push(CircularLabel{angle:rotation+*pos as f32/snap.len as f32*std::f32::consts::TAU,width:size.x,height:size.y});
            kinds.push(RingLabelKind::Site(i));
        }
        for (i,primer) in snap.primers.iter().enumerate() {
            let (a,b)=feature_range(&primer.segs,snap.len,true).unwrap_or((0,0));
            let span=span_len(a,b,snap.len);
            let size=p.layout_no_wrap(primer_label(primer,snap.len),map_font(),theme::PRIMER).size()+Vec2::new(4.0,0.0);
            inputs.push(CircularLabel{angle:rotation+(a as f32+span as f32*0.5)/snap.len as f32*std::f32::consts::TAU,width:size.x,height:size.y});
            kinds.push(RingLabelKind::Primer(i));
        }
        let geometry=CircularLabelLayout::fit(b,&inputs,zoom);
        let delta=(radius-geometry.radius).abs();radius=geometry.radius;
        result=Some(RingLabels{geometry,inputs,kinds});
        if delta<0.5 {break;}
    }
    result.unwrap()
}

fn snapshot(app: &HelixApp) -> Option<Snap> {
    let doc = app.doc.as_ref()?;
    let len = doc.len();
    let kind = doc.kind();
    let directional = kind.features_have_direction();
    let feats = doc
        .features()
        .iter()
        .enumerate()
        // 被隐藏的特征（Features 面板的「Toggle visibility」）不参与绘制；
        // 蛋白按 `Region` / `Site`+`Het` / `Bond` 三段开关过滤（`feature_visible`）。
        .filter(|(idx, f)| app.feature_visible(*idx, f))
        .map(|(idx, f)| {
            let mut segs = Vec::new();
            for s in &f.segments {
                for (a, b) in hx_render::split_wrap(s.range.start, s.range.end, len) {
                    segs.push((a, b));
                }
            }
            FeatRow {
                idx,
                name: f.name.clone(),
                kind: f.kind.clone(),
                segs,
                // 蛋白特征没有方向：画成 SnapGene 那样的方块
                forward_arrow: directional && f.has_forward_arrow(),
                reverse_arrow: directional && f.has_reverse_arrow(),
                color: theme::hex(
                    &f.segments
                        .first()
                        .map(|s| s.color.clone())
                        .unwrap_or_default(),
                    theme::FEATURE_FALLBACK,
                ),
            }
        })
        .collect();
    let primers = if app.show_primers {
        doc.primers()
            .iter()
            .enumerate()
            .flat_map(|(pi, p)| {
                let params=hx_formats::annot::parse_hybridization_params(doc.block_payload(5).unwrap_or_default());
                let primary=p.sites.iter().enumerate().max_by(|(_,a),(_,b)|a.melting_temperature.unwrap_or(f64::NEG_INFINITY).total_cmp(&b.melting_temperature.unwrap_or(f64::NEG_INFINITY))).map(|(i,_)|i).unwrap_or(0);
                p.sites.iter().enumerate().filter_map(move |(si,s)| {
                    let tm=s.melting_temperature.unwrap_or_else(||{
                        let n=s.annealed_bases.min(p.sequence.len());
                        hx_algo::primer::tm(p.sequence[p.sequence.len().saturating_sub(n)..].as_bytes())
                    });
                    if si!=primary && (!params.show_additional_five_prime_matches || tm<params.min_melting_temperature) {return None;}
                    let suffix = crate::primer_tail::annealed_is_suffix(&p.sequence, s, doc.sequence());
                    if si != primary && s.annealed_bases < if suffix { params.min_continuous_match_len } else { params.minimum_five_prime_annealing } { return None; }
                    Some(PrimerRow {
                    idx: pi,
                    name: p.name.clone(),
                    segs: hx_render::split_wrap(s.range.start, s.range.end, len),
                    forward: s.strand != 1,
                })})
            })
            .filter(|p| !p.segs.is_empty())
            .collect()
    } else {
        Vec::new()
    };
    // 位点：同一切点的多个酶合并成 SnapGene 的「A - B (pos)」写法；
    // 编号 = 识别起点(0-based) + 切点偏移（见 `hx_enzymes::display_pos`）。
    let circular = doc.circular();
    let mut raw: Vec<(String, usize, bool, bool)> = Vec::new();
    if app.show_sites && doc.is_dna() {
        for (k, &eidx) in app.sites.selected.iter().enumerate() {
            if app.hidden_enzymes.contains(&k) {
                continue;
            }
            let Some(e) = hx_enzymes::ENZYMES.get(eidx) else {
                continue;
            };
            let unique = app.sites.hits.get(k).map(|h| h.len() == 1).unwrap_or(false);
            for &pos in app.sites.hits.get(k).map(|v| v.as_slice()).unwrap_or(&[]) {
                raw.push((
                    e.name.to_string(),
                    hx_enzymes::display_pos(pos as usize, e.cut, len, circular),
                    unique,
                    app.highlighted_enzymes.contains(&k),
                ));
            }
        }
    }
    raw.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut sites: Vec<(String, usize, bool, bool)> = Vec::new();
    let mut i = 0usize;
    while i < raw.len() {
        let pos = raw[i].1;
        let mut names: Vec<&str> = Vec::new();
        let mut unique = true;
        let mut hl = false;
        let mut j = i;
        while j < raw.len() && raw[j].1 == pos {
            if !names.contains(&raw[j].0.as_str()) {
                names.push(&raw[j].0);
            }
            unique &= raw[j].2;
            hl |= raw[j].3;
            j += 1;
        }
        // SnapGene 标签最多并置 3 个名字；优先常用/偏好名，其次短名
        names.sort_by_key(|n| {
            let pref = hx_enzymes::PREFERRED
                .iter()
                .any(|p| p.eq_ignore_ascii_case(n));
            let com = hx_enzymes::COMMON.iter().any(|p| p.eq_ignore_ascii_case(n));
            (!pref as u8, !com as u8, n.len())
        });
        names.truncate(3);
        sites.push((names.join(" - "), pos, unique, hl));
        i = j;
    }
    // ORF（SnapGene 的 “Show translations and ORFs”）
    let orfs: Vec<(usize, usize, i8, bool, u8)> = if app.show_translation && doc.kind().is_nucleotide() {
        app.displayed_orfs(doc.sequence(), circular)
            .into_iter()
            .map(|o| {
                // 阅读框 0/1/2：正向按起点位置，反向按编码起点（= 该 ORF 的终点）
                let frame = if o.strand > 0 {
                    (o.start % 3) as u8
                } else {
                    ((len.wrapping_sub(o.end)) % 3) as u8
                };
                (o.start, o.end, o.strand, o.has_start, frame)
            })
            .filter(|o| app.translation_frames & (1 << (o.4 + if o.2 < 0 { 3 } else { 0 })) != 0)
            .collect()
    } else {
        Vec::new()
    };
    Some(Snap {
        len,
        circular,
        // 蛋白没有「环状布局」（oracle：SnapGene 的蛋白图谱只有直线，左侧也没有环/线按钮）；
        // **RNA 有**（用户 2026-10-18 指出的 图一：线性 RNA 也能切环状查看）。
        ring: kind != hx_core::MoleculeKind::Protein && app.map_layout.is_ring(circular, len),
        name: doc.map_name(),
        unit: doc.unit(),
        is_dna: doc.is_dna(),
        single_stranded: doc.is_single_stranded(),
        nucleotide: kind.is_nucleotide(),
        feats,
        primers,
        sites,
        orfs,
    })
}

pub fn draw(ui: &mut egui::Ui, app: &mut HelixApp) -> egui::Rect {
    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    if rect.width() < 60.0 || rect.height() < 60.0 {
        return rect;
    }
    let Some(snap) = snapshot(app) else {
        return rect;
    };
    if snap.len == 0 {
        return rect;
    }

    let rotation = app.rotation;
    let mut zoom = app.zoom;
    let mut pan_bp = app.pan_bp;
    let mut hover: Option<usize> = None;

    // ---- 视图参数 ----
    let half = rect.width().min(rect.height()) * 0.5;
    // Circular radius is resolved below from measured callout bounds.
    let r_out=(half-50.0).max(24.0);
    let mut band = if snap.ring {
        (r_out * 0.04).clamp(BAND_MIN, 18.0)
    } else {
        (rect.height() * 0.055).clamp(11.0, 22.0)
    };
    // 线性布局：细轴在上，特征带挂在轴**下方**（SnapGene 排布）；`AXIS_GAP` 留给轴下刻度数字。
    // 底部再留 26 px 给 SnapGene 的横向平移/缩放条。
    const AXIS_GAP: f32 = 16.0;
    const PAN_BAR_H: f32 = 26.0;
    // 克隆对话框里的线性质粒图不画底部的平移/缩放条（SnapGene 的对话框没有）。
    let pan_h = if app.cloning_canvas { 0.0 } else { PAN_BAR_H };
    let lin_bottom = if snap.ring {
        rect.bottom()
    } else {
        rect.bottom() - pan_h
    };
    let px_per_bp = zoom.clamp(0.002, 60.0);
    let x_origin = rect.left() + 40.0 - pan_bp as f32 * px_per_bp;
    // 蛋白图谱里 `type="Site"` / `"Het"` / `"Bond"` 的特征（SnapGene 的 Site 符号）
    // 不当作普通条：画细刻度 + 轴上方纯文字标签（oracle：THYX_MYCTU / KPYK1_ECOLI 的 Site 栈）。
    // 多片段的 Site 另有「外框 + 每片段一根刻度」（oracle：KPYK1 的
    // `domain interfaces [active]`）。
    let tick_mask: Vec<bool> = snap
        .feats
        .iter()
        .map(|f| {
            !snap.ring && !snap.nucleotide && !matches!(feat_group(&f.kind), FeatGroup::Region)
        })
        .collect();
    let lane_spans: Vec<Vec<(usize, usize)>> = snap
        .feats
        .iter()
        .enumerate()
        .map(|(i, f)| {
            if !tick_mask.get(i).copied().unwrap_or(false) {
                return f.segs.clone();
            }
            // 单片段 Site 只画一根刻度，不占车道；多片段的要一条车道画外框
            if f.segs.len() > 1 {
                let s = f.segs.iter().map(|s| s.0).min().unwrap_or(0);
                let e = f.segs.iter().map(|s| s.1).max().unwrap_or(s);
                vec![(s, e)]
            } else {
                Vec::new()
            }
        })
        .collect();
    let (feat_level, feature_count) = interval_lanes(&lane_spans);
    let level_count = feature_count.max(1);
    let (orf_levels, orf_count) = linear_orf_lanes(&snap.orfs, snap.len);
    // Keep callouts away from the floating view tabs/tools and footer. Fit from
    // actual external labels, including primers, rather than all feature count.
    let ring_labels=if snap.ring {Some(ring_labels(&snap,&painter,
        Rect::from_min_max(rect.min+Vec2::new(46.0,62.0),rect.max-Vec2::new(18.0,34.0)),
        rotation,zoom,&feat_level,feature_count,orf_count))}else{None};
    let r_out=ring_labels.as_ref().map(|l|l.geometry.radius).unwrap_or(r_out);
    if snap.ring {band=(r_out*0.04).clamp(BAND_MIN,18.0);}
    let tracks = RingTracks::new(r_out, band, feature_count, orf_count);
    let band = if snap.ring { tracks.band } else { band };
    let feature_offset = if orf_count == 0 {
        AXIS_GAP + band * 0.5
    } else {
        LINEAR_ORF_TOP + (orf_count - 1) as f32 * LINEAR_ORF_STEP + 12.0 + band * 0.5
    };
    // 条下标签行：SnapGene 把「写不进条里」的特征名排在条下方的行里（细引线连回条上），
    // 按文件顺序贪婪占位；两行都占满才退回轴上方分层（即 `Lab::Feat` / `Lab::LongFeat`）。
    const ROW_H: f32 = 13.0;
    const ROW_GAP: f32 = 5.0;
    const MAX_BELOW_ROWS: usize = 2;
    // 蛋白图谱在条下写特征起止数字（oracle：THYX_MYCTU）；核苷酸图谱不写，
    // 蛋白也可以用 `Show translation numbering` 关掉（oracle：SG 的 `(123)` 按钮）。
    let num_row_h: f32 = if snap.nucleotide || !app.show_trans_numbering {
        0.0
    } else {
        13.0
    };
    let label_font = FontId::proportional(11.5);
    let mut used_rows: Vec<Vec<Vec<(f32, f32)>>> =
        vec![vec![Vec::new(); MAX_BELOW_ROWS]; level_count];
    let mut rows_used: Vec<usize> = vec![0; level_count];
    let mut below: Vec<BelowLabel> = Vec::new();
    // 特征已在 `snapshot` 里按类别开关过滤（`app.feature_visible`），这里只看布局
    // （环状布局的标签走环外分层，不用条下行）。
    if !snap.ring {
        // 文件顺序：先来的特征先占行（与 oracle 的贪婪结果一致）
        for (fi, f) in snap.feats.iter().enumerate() {
            if tick_mask.get(fi).copied().unwrap_or(false) {
                continue;
            }
            let Some((s, e)) = f.segs.first().copied() else {
                continue;
            };
            if e <= s {
                continue;
            }
            let tw = painter
                .layout_no_wrap(f.name.clone(), label_font.clone(), theme::TEXT)
                .size()
                .x;
            let span_px = (e - s) as f32 * px_per_bp;
            // 条看不见（< NO_BAR_PX，只有 DNA 图谱才有「小胶囊」这一档）⇒ 轴上方的小胶囊；
            // 条里写得下 ⇒ 不用条下标签；名字长过画布 ⇒ 交给轴上方分层
            let no_bar_px = if snap.is_dna { NO_BAR_PX } else { 0.0 };
            if span_px < no_bar_px || tw + 8.0 <= span_px || tw + 12.0 > rect.width() - 10.0 {
                continue;
            }
            let lane = feat_level.get(fi).copied().unwrap_or(0).min(level_count - 1);
            let cx = x_origin + (s + e) as f32 * 0.5 * px_per_bp;
            'rows: for row in 0..MAX_BELOW_ROWS {
                for dx in [
                    0.0, -8.0, 8.0, -16.0, 16.0, -24.0, 24.0, -32.0, 32.0, -40.0, 40.0, -48.0, 48.0,
                ] {
                    let xc = (cx + dx)
                        .clamp(rect.left() + tw * 0.5 + 5.0, rect.right() - tw * 0.5 - 5.0);
                    let (l, r) = (xc - tw * 0.5 - 3.0, xc + tw * 0.5 + 3.0);
                    if used_rows[lane][row]
                        .iter()
                        .all(|&(a2, b2)| l >= b2 || r <= a2)
                    {
                        used_rows[lane][row].push((l, r));
                        rows_used[lane] = rows_used[lane].max(row + 1);
                        below.push(BelowLabel {
                            lane,
                            row,
                            xc,
                            cx,
                            tw,
                            idx: f.idx,
                            name: f.name.clone(),
                        });
                        break 'rows;
                    }
                }
            }
        }
    }
    // 车道偏移表：每车道 = 条带 + 条下（间隔 + 数字行 + 用到的标签行 + 4）
    let mut lane_dy = vec![0.0f32; level_count];
    for l in 1..level_count {
        lane_dy[l] = lane_dy[l - 1]
            + band
            + ROW_GAP
            + num_row_h
            + rows_used[l - 1] as f32 * ROW_H
            + 4.0;
    }
    let linear_depth = feature_offset
        + lane_dy[level_count - 1]
        + band * 0.42
        + ROW_GAP
        + num_row_h
        + rows_used[level_count - 1] as f32 * ROW_H
        + 4.0;
    let axis_y = (rect.top() + lin_bottom) * 0.5 - linear_depth * 0.5;
    let view = if snap.ring {
        let center = ring_labels.as_ref().map(|l|Pos2::from(l.geometry.center)).unwrap_or(rect.center());
        let mut v = MapView::circular([center.x, center.y], tracks.feature_radius);
        v.ring_width = band;
        v.rotation = rotation;
        v
    } else {
        let mut v = MapView::linear(x_origin, axis_y + feature_offset, px_per_bp);
        v.ring_width = band;
        v
    };
    let c = Pos2::new(view.center[0], view.center[1]);

    // 绘制与命中测试使用相同的轨道半径/行坐标。
    let view_at = |l: usize| -> MapView {
        let mut v = view;
        if snap.ring {
            v.radius -= l as f32 * tracks.feature_step;
        } else {
            v.axis_y += lane_dy.get(l).copied().unwrap_or(0.0);
        }
        v
    };
    let has_below = |idx: usize| below.iter().any(|b| b.idx == idx);
    // 车道条带中心 y / 条下标签行文字中心 y
    let lane_y = |lane: usize| -> f32 {
        axis_y + feature_offset + lane_dy.get(lane).copied().unwrap_or(0.0)
    };
    let below_y = |lane: usize, row: usize| -> f32 {
        lane_y(lane) + band * 0.42 + ROW_GAP + num_row_h + row as f32 * ROW_H + ROW_H * 0.5
    };
    let orf_radius = |i: usize| tracks.orf_radius - orf_levels[i] as f32 * tracks.orf_step;

    // ---- 交互 ----
    if resp.dragged() {
        let d = resp.drag_delta();
        // SnapGene 的环状图谱**固定不旋转**（改原点只在 Actions ▸ Set Origin…）；
        // 之前拖动会旋转质粒，标签层与特征径向层随角度重排，看起来像“特征莫名变化”。
        if !snap.ring {
            pan_bp = (pan_bp - d.x as f64 / view.px_per_bp.max(0.001) as f64).max(0.0);
        }
    }
    if resp.hovered() {
        let (scroll, shift) = ui.input(|i| (i.smooth_scroll_delta.y, i.modifiers.shift));
        if scroll.abs() > 0.1 {
            if snap.ring || !shift {
                zoom = (zoom * (1.0 + scroll * 0.0015)).clamp(0.05, 40.0);
            } else {
                pan_bp = (pan_bp - scroll as f64 / view.px_per_bp.max(0.001) as f64).max(0.0);
            }
        }
    }

    // 命中测试：必须靠近环带/轴带，且取包含该碱基、跨度最小的特征
    let hit_at = |p: Pos2| -> Option<usize> {
        let mut best: Option<(usize, usize)> = None;
        for (fi, f) in snap.feats.iter().enumerate() {
            let fv = view_at(feat_level[fi]);
            let (bp_f, track) = fv.hit(snap.len, [p.x, p.y]);
            if track.abs() > if snap.ring { 0.8 } else { 0.5 } {
                continue;
            }
            if !snap.ring && !(0.0..snap.len as f64).contains(&bp_f) {
                continue;
            }
            let bp = bp_f.floor().max(0.0) as usize;
            for &(a, b) in &f.segs {
                if a <= bp && bp < b {
                    let span = b - a;
                    if best.map(|(_, bs)| span < bs).unwrap_or(true) {
                        best = Some((f.idx, span));
                    }
                }
            }
        }
        best.map(|(i, _)| i)
    };
    let hit_orf = |pos: Pos2| -> Option<usize> {
        let (bp, _) = view.hit(snap.len, [pos.x, pos.y]);
        snap.orfs
            .iter()
            .enumerate()
            .filter(|(i, (s, e, ..))| {
                let near = if snap.ring {
                    (pos.distance(c) - orf_radius(*i)).abs() <= 5.0
                } else {
                    (pos.y - (axis_y + LINEAR_ORF_TOP + orf_levels[*i] as f32 * LINEAR_ORF_STEP))
                        .abs()
                        <= 5.0
                };
                near && hx_render::split_wrap(*s, *e, snap.len)
                    .iter()
                    .any(|&(a, b)| bp >= a as f64 && bp < b as f64)
            })
            .min_by_key(|(_, (s, e, ..))| span_len(*s, *e, snap.len))
            .map(|(i, _)| i)
    };
    // 引物：图谱上的细弧/细线（与绘制同一套半径），点中即选中、双击进 `Edit Primer…`。
    let hit_primer = |pos: Pos2| -> Option<usize> {
        let (bp, _) = view.hit(snap.len, [pos.x, pos.y]);
        snap.primers
            .iter()
            .filter(|pr| {
                let near = if snap.ring {
                    let radius = r_out + if pr.forward { 12.0 } else { 25.0 };
                    (pos.distance(c) - radius).abs() <= 6.0
                } else {
                    let y = if pr.forward { axis_y - 4.0 } else { axis_y + 4.0 };
                    (pos.y - y).abs() <= 6.0
                };
                near && pr
                    .segs
                    .iter()
                    .any(|&(a, b)| bp >= a as f64 && bp < b as f64)
            })
            .min_by_key(|pr| pr.segs.iter().map(|&(a, b)| b - a).sum::<usize>())
            .map(|pr| pr.idx)
    };
    let near_backbone = |pos: Pos2| (pos.distance(c) - r_out).abs() <= 12.0;
    let drag_id = resp.id.with("dna-ring-selection");
    // 标签命中区在 egui 临时数据里的键（标签矩形在绘制末尾写入，点击处理读上一帧的）。
    let label_hits_id = resp.id.with("map-label-hits");
    let primer_hits_id=resp.id.with("map-primer-label-hits");
    let primer_at=|pos:Pos2|hit_primer(pos).or_else(||ui.ctx().data(|d|
        d.get_temp::<Vec<(Rect,usize)>>(primer_hits_id).and_then(|v|v.into_iter().find(|(r,_)|r.contains(pos)).map(|(_,i)|i))));
    if snap.ring && resp.drag_started_by(egui::PointerButton::Primary) {
        if let Some(origin) = ui.input(|i| i.pointer.press_origin()) {
            if rect.contains(origin) && near_backbone(origin) {
                let (bp, _) = view.hit(snap.len, [origin.x, origin.y]);
                ui.ctx()
                    .data_mut(|d| d.insert_temp(drag_id, RingDrag::new(bp)));
                app.selected_feature = None;
                app.selected_residue = None;
                app.cursor_bp = None;
                app.sel_anchor = None;
            }
        }
    }
    if snap.ring && (resp.dragged_by(egui::PointerButton::Primary) || resp.drag_stopped()) {
        if let (Some(mut drag), Some(pos)) = (
            ui.ctx().data(|d| d.get_temp::<RingDrag>(drag_id)),
            resp.interact_pointer_pos(),
        ) {
            // 靠近圆心时角度没有操作意义，保留上一次的端点。
            if pos.distance(c) > r_out * 0.2 {
                let (bp, _) = view.hit(snap.len, [pos.x, pos.y]);
                app.selection = drag.update(bp, snap.len).map(|(s, e)| {
                    if !snap.circular && s > e {
                        (e, s)
                    } else {
                        (s, e)
                    }
                });
                ui.ctx().data_mut(|d| d.insert_temp(drag_id, drag));
                if let Some((s, e)) = app.selection {
                    app.status = format!(
                        "Selected {}..{} ({} bp)",
                        nfmt(s + 1),
                        nfmt(e),
                        nfmt(span_len(s, e, snap.len))
                    );
                }
            }
        }
    }
    if resp.drag_stopped() || !ui.input(|i| i.pointer.primary_down()) {
        ui.ctx().data_mut(|d| d.remove::<RingDrag>(drag_id));
    }
    if resp.double_clicked_by(egui::PointerButton::Primary) {
        if let Some(pos) = resp.interact_pointer_pos() {
            if let Some(pi) = primer_at(pos) {
                app.selected_primer = Some(pi);
                app.open_edit_primer_dialog();
            }
        }
    }
    if resp.clicked_by(egui::PointerButton::Primary) {
        if let Some(pos) = resp.interact_pointer_pos() {
            app.selected_residue = None;
            app.cursor_bp = None;
            app.sel_anchor = None;
            // 先看条/轴，再看**标签**（用户 2026-10-18：「点击标签不会产生高亮」——
            // SnapGene 里点特征名、点小特征的胶囊标签都等于选中该特征）。
            let label_hit = ui.ctx().data(|d| {
                d.get_temp::<Vec<(Rect, usize)>>(label_hits_id)
            });
            let feature = hit_at(pos).or_else(|| {
                label_hit.and_then(|v| {
                    v.iter()
                        .rev()
                        .find(|(r, _)| r.expand(1.5).contains(pos))
                        .map(|(_, i)| *i)
                })
            });
            app.selected_feature = feature;
            if let Some(f) = feature.and_then(|i| snap.feats.iter().find(|f| f.idx == i)) {
                // 蛋白图谱里 Site/Het/Bond 是**刻度轴**（`tick_mask`），选中它们时 SG 只把
                // 轴线与标签变蓝，**不铺选区蓝带**（用户 2026-10-19 裁定「改成与 SG 完全一致」；
                // 证据 `/tmp/opencode/sg/lb-02-click-site-label.png`：SG 状态栏只有
                // `Selected: … (p3 .. p414 / 39 aa)`，没有 Selection 段）。
                let fi = snap.feats.iter().position(|x| x.idx == f.idx);
                let is_tick = fi.and_then(|i| tick_mask.get(i).copied()).unwrap_or(false);
                app.selection = if is_tick {
                    None
                } else {
                    feature_range(&f.segs, snap.len, snap.circular)
                };
                app.status = app
                    .selection
                    .map(|(s, e)| {
                        format!(
                            "{} · {}..{} ({} bp)",
                            f.name,
                            nfmt(s + 1),
                            nfmt(e),
                            nfmt(span_len(s, e, snap.len))
                        )
                    })
                    .unwrap_or_else(|| f.name.clone());
            } else if let Some(pi) = primer_at(pos) {
                // 引物：与序列视图同一套选中文案（SnapGene 实测）。
                app.selected_feature = None;
                app.selection = None;
                if let Some(pr) = app.doc.as_ref().and_then(|d| d.primers().get(pi)) {
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
                app.selected_primer = Some(pi);
            } else if let Some(i) = hit_orf(pos) {
                let (s, e, strand, ..) = snap.orfs[i];
                app.selection = Some((s, e));
                app.status = format!(
                    "ORF ({}) · {}..{} ({} bp)",
                    if strand > 0 { "+" } else { "−" },
                    nfmt(s + 1),
                    nfmt(e),
                    nfmt(span_len(s, e, snap.len))
                );
            } else {
                app.selection = None;
                app.status = "Selection cleared".into();
            }
        }
    }
    if let Some(p) = resp.hover_pos() {
        let label_hit = ui
            .ctx()
            .data(|d| d.get_temp::<Vec<(Rect, usize)>>(label_hits_id));
        hover = hit_at(p).or_else(|| {
            label_hit.and_then(|v| {
                v.iter()
                    .rev()
                    .find(|(r, _)| r.expand(1.5).contains(p))
                    .map(|(_, i)| *i)
            })
        });
        let (bp_f, _) = view.hit(snap.len, [p.x, p.y]);
        let orf = hit_orf(p);
        app.hover_bp = (!snap.ring || near_backbone(p) || hover.is_some() || orf.is_some())
            .then(|| (bp_f.floor().max(0.0) as usize).min(snap.len - 1));
        if hover.is_some() || orf.is_some() || primer_at(p).is_some() {
            ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
        } else if snap.ring && near_backbone(p) {
            ui.ctx().set_cursor_icon(egui::CursorIcon::Crosshair);
        }
    } else {
        app.hover_bp = None;
    }

    let p = &painter;

    // 标签命中区（点标签 = 选中对应特征；用户在 2026-10-18 指摘「点击标签不会产生
    // 高亮」）。帧末存进 egui 临时数据，下一帧的点击处理读它 —— 标签位置只依赖布局，
    // 与选中态无关，所以用上一帧的矩形命中是稳定的。
    let mut label_hits: Vec<(Rect, usize)> = Vec::new();
    let mut primer_hits:Vec<(Rect,usize)>=Vec::new();

    // ---- 骨架 ----
    if snap.ring {
        // 骨架线宽 oracle 实测：SnapGene 环状骨架 **3 px**（我们原来是 1.35 px，用户
        // 2026-10-18 说「希望我们的黑色直线或者黑色圆线可以和 SnapGene 一样粗」）；
        // 双链 DNA 两根线（中心相距 4 px、各 3 px），单链（RNA / 单链 DNA）只有一根。
        let core = if snap.single_stranded {
            vec![0.0]
        } else {
            vec![-2.0, 2.0]
        };
        for dr in core {
            let radius = r_out + dr;
            if snap.circular {
                p.circle_stroke(c, radius, Stroke::new(RING_LINE_W, theme::BACKBONE));
            } else {
                // 线性分子按环状布局：正上方留缺口表示分子的两端（oracle 实测 17.3°）。
                let mut rv = view;
                rv.radius = radius;
                let (main, _) = rv.open_ring_arcs();
                p.add(Shape::line(
                    main.into_iter().map(|q| Pos2::new(q[0], q[1])).collect(),
                    Stroke::new(RING_LINE_W, theme::BACKBONE),
                ));
            }
        }
    } else {
        // 线性轴：两条 3 px 粗线（SnapGene 的线性轴比环状骨架粗），只画在序列长度范围内。
        // 蛋白图谱没有骨架轴（oracle：KPYK1_ECOLI 的图谱下方直接就是特征条）。
        let x0 = view.bp_x(0.0).max(rect.left() + 2.0);
        let x1 = view.bp_x(snap.len as f64).min(rect.right() - 2.0);
        // Protein backbone represents the complete molecule, not the union of
        // its annotations. Range features may cover only a small internal span.
        if !snap.nucleotide && x1>x0 {
            p.line_segment([Pos2::new(x0,axis_y),Pos2::new(x1,axis_y)],Stroke::new(2.4,theme::BACKBONE));
            for (at,text,align) in protein_end_labels(snap.len) {
                let x=view.bp_x(at as f64);
                if x>=rect.left()+2.0 && x<=rect.right()-2.0 {
                    p.line_segment([Pos2::new(x,axis_y-4.0),Pos2::new(x,axis_y+5.0)],Stroke::new(1.3,theme::BACKBONE));
                    p.text(Pos2::new(x,axis_y-9.0),align,text,FontId::proportional(12.0),theme::TEXT);
                }
            }
        }
        if snap.nucleotide && x1 > x0 {
            if snap.single_stranded {
                // 单链（RNA / 单链 DNA）没有双链骨架：只画一条标尺线（oracle：eco16s 的
                // 图谱上方是一根细线，刻度朝下、数字写在线上方）。
                p.line_segment(
                    [Pos2::new(x0, axis_y), Pos2::new(x1, axis_y)],
                    Stroke::new(2.0, theme::BACKBONE),
                );
            } else {
                for dy in [-2.0, 2.0] {
                    p.line_segment(
                        [Pos2::new(x0, axis_y + dy), Pos2::new(x1, axis_y + dy)],
                        Stroke::new(2.4, theme::BACKBONE),
                    );
                }
            }
        }
        // 两端竖向端标（起点/终点；SnapGene 的端标比刻度长）
        if snap.nucleotide {
            for x in [view.bp_x(0.0), view.bp_x(snap.len as f64)] {
                if (rect.left()..=rect.right()).contains(&x) {
                    p.line_segment(
                        [Pos2::new(x, axis_y - 4.0), Pos2::new(x, axis_y + 11.0)],
                        Stroke::new(2.4, theme::BACKBONE),
                    );
                }
            }
            // 环状分子按线性显示：两端画「⋯」，表示序列在环上继续（SnapGene 的写法）
            if snap.circular {
                for i in 0..3 {
                    let dx = 20.0 + i as f32 * 14.0;
                    for x in [view.bp_x(0.0) - dx, view.bp_x(snap.len as f64) + dx] {
                        p.rect_filled(
                            Rect::from_center_size(Pos2::new(x, axis_y), Vec2::new(5.0, 5.0)),
                            0.0,
                            theme::BACKBONE,
                        );
                    }
                }
            }
        }
    }

    // ---- 特征带 ----
    let selected = app.selected_feature;
    // 极小特征（外置标签）：(角度 或 x, 宽度, 特征下标, 引线起点)
    let mut tiny: Vec<(f32, f32, usize, Pos2)> = Vec::new();
    // 线性蛋白图谱：短到画不出条的特征（SnapGene 的 `Site`，1–2 aa）画成
    // 细刻度 + 轴上方**纯文字**标签（oracle：THYX_MYCTU 的 Site 栈）。
    let mut tick_feats: Vec<(f32, f32, usize)> = Vec::new();
    // 每个车道「条下数字行」已占用的横向区间（SnapGene 放不下就省略数字）
    let mut num_used: std::collections::HashMap<usize, Vec<Rect>> = std::collections::HashMap::new();
    for (fi, f) in snap.feats.iter().enumerate() {
        let feature_color = if app.show_colors {
            f.color
        } else {
            theme::FEATURE_NO_COLOR
        };
        let fill = theme::band_fill(feature_color);
        let edge = theme::darken(feature_color, 0.5);
        let is_sel = selected == Some(f.idx);
        let fv = view_at(feat_level.get(fi).copied().unwrap_or(0));
        let view = &fv;
        for (si, (a, b)) in f.segs.iter().enumerate() {
            let (a, b) = (*a, *b);
            if a >= b {
                continue;
            }
            let (mut a0_ang, mut a1_ang) = if snap.ring {
                let x0 = view.bp_angle(snap.len, a as f64);
                let mut x1 = view.bp_angle(snap.len, b as f64);
                while x1 <= x0 {
                    x1 += std::f32::consts::TAU;
                }
                (x0, x1)
            } else {
                (0.0, 0.0)
            };
            // 线性分子按环状布局：特征带在正上方缺口处**断开**（oracle 实测：SnapGene
            // 的绿带 rRNA 与骨架缺口同为 17.3°；我们原来整条带绕满一圈，看不出一端）。
            if snap.ring && !snap.circular {
                let (lo, hi) = view.open_ring_main_arc();
                if a1_ang <= lo || a0_ang >= hi {
                    continue;
                }
                a0_ang = a0_ang.max(lo);
                a1_ang = a1_ang.min(hi);
            }
            let arc_px = if snap.ring {
                (a1_ang - a0_ang) * view.radius
            } else {
                (b - a) as f32 * view.px_per_bp
            };
            let mid = (a + b) / 2;
            let mid_ang = if snap.ring {
                view.bp_angle(snap.len, mid as f64)
            } else {
                0.0
            };
            // 小特征**不改变几何**（SnapGene 按真实序列跨度按比例绘制，放大后能看到
            // 18 bp 与 25 bp 的长短差别）；这里只是登记一次外置标签。
            // 已经排进「条下标签行」的特征不再参与轴上方分层。
            let is_tick_feat = tick_mask.get(fi).copied().unwrap_or(false);
            if is_tick_feat && !tick_feats.iter().any(|t| t.2 == f.idx) {
                // 蛋白图谱的 Site：细刻度 + 轴上方纯文字标签
                tick_feats.push((view.bp_x(mid as f64), 92.0, f.idx));
            }
            let no_bar = if snap.ring {
                TINY_PX
            } else if snap.is_dna {
                NO_BAR_PX
            } else {
                // RNA / 蛋白图谱没有「小胶囊」这一档：小特征画小方块 + 条下标签
                0.0
            };
            if !is_tick_feat
                && arc_px < no_bar
                && !tiny.iter().any(|t| t.2 == f.idx)
                && !has_below(f.idx)
            {
                // DNA/RNA（以及环状视图）：带颜色的小胶囊标签
                let pos = if snap.ring {
                    pt(c, mid_ang, view.radius)
                } else {
                    Pos2::new(view.bp_x(mid as f64), view.axis_y)
                };
                // 宽度按真实文本量（`draw_chip` 的盒子 = 文本 + 左右各 5 px），
                // 否则窄名字会虚占位置、宽名字会与酶标签叠在一起
                let chip_w = p
                    .layout_no_wrap(f.name.clone(), map_font(), theme::TEXT)
                    .size()
                    .x
                    + 10.0;
                tiny.push((if snap.ring { mid_ang } else { pos.x }, chip_w, f.idx, pos));
            }

            if snap.ring {
                // SnapGene 的特征箭头是**单体五边形**：平口尾 → 等厚带体 →
                // **比带体宽的喇叭形箭头** → 尖端；整条轮廓连续。
                // 填充必须拆成凸块（带体 quads + 头部三角），轮廓用 `closed_line`
                // （允许非凸），否则 egui 会把非凸多边形三角化成楔形。
                let head_px = (arc_px * 0.45).clamp(2.0, 18.0);
                let head_count = usize::from(f.forward_arrow) + usize::from(f.reverse_arrow);
                let dth = (head_px / view.radius.max(8.0))
                    .min((a1_ang - a0_ang) * 0.7)
                    .min((a1_ang - a0_ang) / head_count.max(1) as f32);
                let c = Pos2::new(view.center[0], view.center[1]);
                let r_mid = view.radius;
                let half = view.ring_width * 0.46;
                let (r_in, r_out) = (r_mid - half, r_mid + half);
                // 头部基部比带体宽（喇叭）
                const FLARE: f32 = 1.55;
                let (h_in, h_out) = (r_mid - half * FLARE, r_mid + half * FLARE);
                let stroke = Stroke::new(
                    if is_sel { 2.0 } else { 0.9 },
                    if is_sel { theme::ACCENT } else { edge },
                );
                // 带体：自己切「相邻**略微重叠**的凸四边形」（重叠 ~0.6 px 消除抗锯齿接缝，
                // 否则一条条四边形之间会出现浅色细缝＝颜色不均）。
                let b0 = a0_ang + if f.reverse_arrow { dth } else { 0.0 };
                let b1 = a1_ang - if f.forward_arrow { dth } else { 0.0 };
                let body_px = (b1 - b0) * view.radius;
                let overlap = if body_px > 2.0 { 0.6 / body_px } else { 0.02 };
                if body_px > 0.5 {
                    let seg = ((body_px / 3.0) as usize).clamp(2, 240);
                    for i in 0..seg {
                        let t0 = ((i as f32 / seg as f32) - overlap).max(0.0);
                        let t1 = (((i + 1) as f32 / seg as f32) + overlap).min(1.0);
                        let (aa, ab) = (b0 + (b1 - b0) * t0, b0 + (b1 - b0) * t1);
                        p.add(Shape::convex_polygon(
                            vec![
                                pt(c, aa, r_out),
                                pt(c, ab, r_out),
                                pt(c, ab, r_in),
                                pt(c, aa, r_in),
                            ],
                            fill,
                            Stroke::NONE,
                        ));
                    }
                }
                // directionality=0 时两端均为平口；双向特征才在两端各画一个箭头。
                for (enabled, tip_ang, base_ang) in
                    [(f.reverse_arrow, a0_ang, b0), (f.forward_arrow, a1_ang, b1)]
                {
                    if enabled {
                        p.add(Shape::convex_polygon(
                            vec![
                                pt(c, tip_ang, r_mid),
                                pt(c, base_ang, h_out),
                                pt(c, base_ang, h_in),
                            ],
                            fill,
                            Stroke::NONE,
                        ));
                    }
                }
                // 整条轮廓（非凸 ⇒ 用折线）：外弧 → 喇叭外角 → 尖 → 喇叭内角 → 内弧
                let mut poly: Vec<Pos2> = Vec::new();
                let steps = ((body_px / 4.0) as usize).clamp(2, 160);
                // 顺着外弧、右端、内弧、左端走同一条轮廓，选择框也沿此描边。
                for i in 0..=steps {
                    let t = b0 + (b1 - b0) * i as f32 / steps as f32;
                    poly.push(pt(c, t, r_out));
                }
                if f.forward_arrow {
                    poly.push(pt(c, b1, h_out));
                    poly.push(pt(c, a1_ang, r_mid));
                    poly.push(pt(c, b1, h_in));
                }
                for i in (0..=steps).rev() {
                    let t = b0 + (b1 - b0) * i as f32 / steps as f32;
                    poly.push(pt(c, t, r_in));
                }
                if f.reverse_arrow {
                    poly.push(pt(c, b0, h_in));
                    poly.push(pt(c, a0_ang, r_mid));
                    poly.push(pt(c, b0, h_out));
                }
                if poly.len() >= 3 {
                    p.add(Shape::closed_line(poly, stroke));
                }
            } else if is_tick_feat {
                // 蛋白图谱的 Site：细刻度；多片段的先画一条圆角外框（SnapGene 的写法）
                let y = if f.segs.len() > 1 {
                    lane_y(feat_level.get(fi).copied().unwrap_or(0))
                } else {
                    lane_y(level_count - 1)
                };
                if f.segs.len() > 1 {
                    let s0 = f.segs.iter().map(|s| s.0).min().unwrap_or(0);
                    let e0 = f.segs.iter().map(|s| s.1).max().unwrap_or(s0);
                    p.add(Shape::rect_stroke(
                        Rect::from_min_max(
                            Pos2::new(view.bp_x(s0 as f64), y - band * 0.45),
                            Pos2::new(view.bp_x(e0 as f64), y + band * 0.45),
                        ),
                        5.0,
                        Stroke::new(1.0, theme::with_alpha(feature_color, 210)),
                        egui::StrokeKind::Middle,
                    ));
                }
                for (a, b) in &f.segs {
                    let x = view.bp_x((*a + *b) as f64 * 0.5);
                    let (up, dn) = if f.segs.len() > 1 {
                        (band * 0.45, band * 0.45)
                    } else {
                        (band * 1.05, band * 0.55)
                    };
                    p.line_segment(
                        [Pos2::new(x, y - up), Pos2::new(x, y + dn)],
                        Stroke::new(1.2, theme::with_alpha(theme::TEXT_DIM, 170)),
                    );
                }
            } else {
                // 箭头尖端留在真实区间内；无方向特征保留平直的两端。
                let bx0 = view.bp_x(a as f64);
                let bx1 = view.bp_x(b as f64);
                let points = linear_feature_polygon(
                    bx0,
                    bx1,
                    view.axis_y,
                    view.ring_width * 0.42,
                    f.forward_arrow,
                    f.reverse_arrow,
                );
                p.add(Shape::convex_polygon(
                    points,
                    fill,
                    Stroke::new(
                        if is_sel { 2.0 } else { 0.9 },
                        if is_sel { theme::ACCENT } else { edge },
                    ),
                ));
                // 蛋白图谱：条下写特征的起止数字（oracle：`1` / `238`；核苷酸图谱不写）。
                // SnapGene 的蛋白特征坐标从 1 数起（`range="1-250"` 显示 1/250），
                // 我们内部是 0 基半开区间 ⇒ 起点 +1。数字也要让位：同车道里放不下就省略。
                if !snap.nucleotide
                    && app.show_trans_numbering
                    && si == 0
                    && !tick_feats.iter().any(|t| t.2 == f.idx)
                    && feat_level.get(fi).is_some()
                {
                    let lane = feat_level.get(fi).copied().unwrap_or(0);
                    let y = view.axis_y + view.ring_width * 0.42 + ROW_GAP + num_row_h * 0.5;
                    let font = FontId::proportional(9.5);
                    let (ta, tb) = (nfmt_plain(a + 1), nfmt_plain(b));
                    let wa = p
                        .layout_no_wrap(ta.clone(), font.clone(), theme::TEXT)
                        .size()
                        .x;
                    let wb = p
                        .layout_no_wrap(tb.clone(), font.clone(), theme::TEXT)
                        .size()
                        .x;
                    let (ra, rb) = if wa + wb + 14.0 <= bx1 - bx0 {
                        (
                            Rect::from_min_max(
                                Pos2::new(bx0, y - 6.0),
                                Pos2::new(bx0 + wa, y + 6.0),
                            ),
                            Rect::from_min_max(
                                Pos2::new(bx1 - wb, y - 6.0),
                                Pos2::new(bx1, y + 6.0),
                            ),
                        )
                    } else {
                        let cx = (bx0 + bx1) * 0.5;
                        (
                            Rect::from_min_max(
                                Pos2::new(cx - 4.0 - wa, y - 6.0),
                                Pos2::new(cx - 4.0, y + 6.0),
                            ),
                            Rect::from_min_max(
                                Pos2::new(cx + 4.0, y - 6.0),
                                Pos2::new(cx + 4.0 + wb, y + 6.0),
                            ),
                        )
                    };
                    let used = num_used.entry(lane).or_default();
                    if !used.iter().any(|r: &Rect| r.intersects(ra))
                        && !used.iter().any(|r: &Rect| r.intersects(rb))
                    {
                        used.push(ra);
                        used.push(rb);
                        let (pa, pb) = if wa + wb + 14.0 <= bx1 - bx0 {
                            (
                                (Pos2::new(bx0, y), Align2::LEFT_CENTER),
                                (Pos2::new(bx1, y), Align2::RIGHT_CENTER),
                            )
                        } else {
                            let cx = (bx0 + bx1) * 0.5;
                            (
                                (Pos2::new(cx - 4.0, y), Align2::RIGHT_CENTER),
                                (Pos2::new(cx + 4.0, y), Align2::LEFT_CENTER),
                            )
                        };
                        p.text(pa.0, pa.1, ta, font.clone(), theme::TEXT);
                        p.text(pb.0, pb.1, tb, font, theme::TEXT);
                    }
                }
            }
        }
    }

    // ---- 引物结合位点 ----
    // 引物使用独立的细轨道，避免与特征箭头混在一起；其显示由同一个
    // `Show primers` 开关控制。颜色开关关闭时退回中性灰。
    for primer in &snap.primers {
        let primer_color = if app.show_colors {
            theme::PRIMER
        } else {
            theme::FEATURE_NO_COLOR
        };
        let radius = r_out + if primer.forward { 12.0 } else { 25.0 };
        // 线性图谱按真 SnapGene 的做法：小箭头直接骑在骨架线上（正向在上、
        // 反向在下），名字单独占骨架上方一条道 —— 原来和酶位点标签挤同一行，
        // 引物几乎看不见。
        let y = if snap.ring {
            axis_y - if primer.forward { 16.0 } else { 32.0 }
        } else if primer.forward {
            axis_y - 4.0
        } else {
            axis_y + 4.0
        };

        for (si, &(a, b)) in primer.segs.iter().enumerate() {
            let terminal = if primer.forward {
                si + 1 == primer.segs.len()
            } else {
                si == 0
            };
            let span_px = if snap.ring {
                (b - a) as f32 / snap.len as f32 * std::f32::consts::TAU * radius
            } else {
                (b - a) as f32 * view.px_per_bp
            };
            if snap.ring {
                let (points, head) = orf_arrow(&view, snap.len, a, b, primer.forward, radius);
                p.add(Shape::line(points, Stroke::new(2.0, primer_color)));
                if terminal {
                    p.add(Shape::convex_polygon(head, primer_color, Stroke::NONE));
                }
            } else {
                let x0 = view.bp_x(a as f64);
                let x1 = view.bp_x(b as f64);
                p.line_segment(
                    [Pos2::new(x0, y), Pos2::new(x1, y)],
                    Stroke::new(2.0, primer_color),
                );
                if terminal {
                    let tip = if primer.forward { x1 } else { x0 };
                    let direction = if primer.forward { 1.0 } else { -1.0 };
                    let size = span_px.clamp(1.0, 7.0);
                    p.add(Shape::convex_polygon(
                        vec![
                            Pos2::new(tip, y),
                            Pos2::new(tip - direction * size, y - 3.0),
                            Pos2::new(tip - direction * size, y + 3.0),
                        ],
                        primer_color,
                        Stroke::NONE,
                    ));
                }
            }
        }
    }

    // ---- 刻度与数字 ----
    // 蛋白图谱没有标尺轴（oracle：THYX_MYCTU / KPYK1_ECOLI 的图谱上方只有特征标签）
    if app.show_ruler && snap.nucleotide {
        let px_per_bp = if snap.ring {
            r_out * std::f32::consts::TAU / snap.len as f32
        } else {
            view.px_per_bp
        };
        // 环状图谱用小刻度细分；线性图谱只有主刻度（oracle：轴在 500/1000 之间是干净的）
        let (min_px, subdiv) = if snap.ring && !snap.is_dna {(300.0,1)}else if snap.ring { (150.0, 10) } else { (100.0, 1) };
        for (bp, major) in ruler_ticks_subdiv(snap.len, px_per_bp, min_px, subdiv) {
            let bp_i = (bp.round() as usize).max(0);
            if snap.ring {
                let a = view.bp_angle(snap.len, bp);
                // 线性分子按环状布局：缺口里不画刻度（那一小段是分子的两端之间）。
                if !snap.circular {
                    let (lo, hi) = view.open_ring_main_arc();
                    if a < lo || a > hi {
                        continue;
                    }
                }
                if snap.is_dna {p.line_segment(
                    [
                        pt(c, a, r_out + 2.0),
                        pt(c, a, r_out - if major { 10.0 } else { 5.0 }),
                    ],
                    Stroke::new(
                        if major { 1.4 } else { 1.0 },
                        if major {
                            theme::TICK_MAJOR
                        } else {
                            theme::TICK_MINOR
                        },
                    ),
                );}
                if major && bp > 0.5 {
                    text_rot(
                        p,
                        &nfmt_plain(bp_i),
                        pt(c, a, r_out - 16.0),
                        readable_angle(a),
                        Align2::CENTER_CENTER,
                        map_font(),
                        theme::TEXT,
                        false,
                    );
                }
            } else {
                let x = view.bp_x(bp);
                if x < rect.left() - 40.0 || x > rect.right() + 40.0 {
                    continue;
                }
                // 主刻度跨轴并向下伸出（SnapGene 线性图谱没有小刻度）；数字在轴下、右对齐到刻度。
                // RNA 的刻度从细标尺线**朝下**伸，数字写在线上方（oracle：eco16s）。
                if snap.is_dna {
                    p.line_segment(
                        [Pos2::new(x, axis_y - 4.0), Pos2::new(x, axis_y + 8.0)],
                        Stroke::new(1.2, theme::TICK_MAJOR),
                    );
                    if major && bp > 0.5 {
                        p.text(
                            Pos2::new(x - 3.0, axis_y + 9.0),
                            Align2::RIGHT_TOP,
                            nfmt_plain(bp_i),
                            FontId::proportional(10.0),
                            theme::RULER_TEXT,
                        );
                    }
                } else {
                    p.line_segment(
                        [Pos2::new(x, axis_y), Pos2::new(x, axis_y + 7.0)],
                        Stroke::new(1.2, theme::TICK_MAJOR),
                    );
                    if major && bp > 0.5 {
                        p.text(
                            Pos2::new(x - 3.0, axis_y - 3.0),
                            Align2::RIGHT_BOTTOM,
                            nfmt_plain(bp_i),
                            FontId::proportional(10.0),
                            theme::RULER_TEXT,
                        );
                    }
                }
            }
        }
    }

    // ---- ORF：环上位于 DNA 环与特征之间；线性位于标尺与特征之间 ----
    if !snap.orfs.is_empty() {
        for (oi, (s0, e0, strand, _has_start, _frame)) in snap.orfs.iter().enumerate() {
            let (a, b) = (*s0, *e0);
            let col = if !app.show_colors {
                theme::FEATURE_NO_COLOR
            } else if *strand > 0 {
                ORF_FWD
            } else {
                ORF_REV
            };
            if snap.ring {
                let r = orf_radius(oi);
                let (pts, head) = orf_arrow(&view, snap.len, a, b, *strand > 0, r);
                if pts.len() >= 2 {
                    p.add(Shape::line(pts, Stroke::new(1.6, col)));
                }
                p.add(Shape::convex_polygon(head, col, Stroke::NONE));
            } else {
                let y = axis_y + LINEAR_ORF_TOP + orf_levels[oi] as f32 * LINEAR_ORF_STEP;
                for (start, end, terminal) in linear_orf_segments(a, b, *strand, snap.len) {
                    let x0 = view.bp_x(start as f64);
                    let x1 = view.bp_x(end as f64);
                    if x1 <= x0 || x1 < rect.left() || x0 > rect.right() {
                        continue;
                    }
                    let tip = if *strand > 0 { x1 } else { x0 };
                    let dir = if *strand > 0 { 1.0 } else { -1.0 };
                    let draw_head = terminal && (rect.left()..=rect.right()).contains(&tip);
                    let head_len = 7.0_f32.min(x1 - x0);
                    let base = tip - dir * head_len;
                    let (sx, ex) = if draw_head && *strand > 0 {
                        (x0, base)
                    } else if draw_head {
                        (base, x1)
                    } else {
                        (x0, x1)
                    };
                    let shaft = [Pos2::new(sx, y), Pos2::new(ex, y)];
                    p.line_segment(shaft, Stroke::new(3.6, theme::with_alpha(col, 65)));
                    p.line_segment(shaft, Stroke::new(1.6, col));
                    if draw_head {
                        p.add(Shape::convex_polygon(
                            vec![
                                Pos2::new(tip, y),
                                Pos2::new(base, y - 3.2),
                                Pos2::new(base, y + 3.2),
                            ],
                            col,
                            Stroke::NONE,
                        ));
                    }
                }
            }
        }
    }

    // ---- 条下标签（线性：写不进条里的特征名 + 细引线）----
    for bl in &below {
        let ly = lane_y(bl.lane);
        let bar_bottom = ly + band * 0.42;
        let y = below_y(bl.lane, bl.row);
        let row_top = y - ROW_H * 0.5 + 1.0;
        let leader = Stroke::new(0.8, theme::with_alpha(theme::TEXT_DIM, 115));
        // 竖引线：条中心 → 标签行顶
        p.line_segment(
            [Pos2::new(bl.cx, bar_bottom), Pos2::new(bl.cx, row_top)],
            leader,
        );
        let (l, r) = (bl.xc - bl.tw * 0.5, bl.xc + bl.tw * 0.5);
        // 标签被挤到一侧时补一段横引线（SnapGene 的 L 形引线）
        if bl.xc - bl.cx > 6.0 {
            p.line_segment(
                [Pos2::new(bl.cx, row_top), Pos2::new(l - 1.0, row_top)],
                leader,
            );
        } else if bl.cx - bl.xc > 6.0 {
            p.line_segment(
                [Pos2::new(bl.cx, row_top), Pos2::new(r + 1.0, row_top)],
                leader,
            );
        }
        let col = if selected == Some(bl.idx) {
            theme::ACCENT
        } else {
            theme::TEXT
        };
        // 条下标签 = 特征名：点它也选中特征（oracle：SnapGene 里点特征名会高亮）。
        label_hits.push((
            Rect::from_min_max(
                Pos2::new(bl.xc - bl.tw * 0.5, y - ROW_H * 0.5),
                Pos2::new(bl.xc + bl.tw * 0.5, y + ROW_H * 0.5),
            ),
            bl.idx,
        ));
        p.text(
            Pos2::new(bl.xc, y),
            Align2::CENTER_CENTER,
            &bl.name,
            label_font.clone(),
            col,
        );
    }

    // ---- 线性外置标签；环状标签使用独立的屏幕矩形求解器 ----
    enum Lab {
        /// 小特征：特征色圆角盒。
        Feat {
            idx: usize,
            name: String,
            color: Color32,
        },
        /// 酶位点：名称加粗、坐标常规（SnapGene 写法）。
        Site {
            name: String,
            pos: usize,
            unique: bool,
            highlighted: bool,
            /// 线性图谱：标签框 `(左, 右, 标签在引线左侧)`——左右镜像排版，
            /// 「位置」永远靠近引线一侧（与 SnapGene 一致）。
            lin: Option<(f32, f32, bool)>,
        },
        /// 线性图谱的引物名（紫色纯文字，和特征/酶位点共用一套分层）。
        Primer { name: String, forward: bool },
        /// 线性视图里轴上放不下的特征名（纯文字，SnapGene 的写法）。
        LongFeat {
            idx: usize,
            name: String,
        },
        /// 线性图谱两端的 Start / End 标注。
        EndMark { start: bool },
    }
    // (键, 宽, 标签)
    let mut items: Vec<(f32, f32, Lab)> = Vec::new();
    for (key, _w, fidx) in &tick_feats {
        if let Some(f) = snap.feats.iter().find(|f| f.idx == *fidx) {
            let w = p
                .layout_no_wrap(f.name.clone(), FontId::proportional(11.5), theme::TEXT)
                .size()
                .x
                + 14.0;
            items.push((
                *key,
                w,
                Lab::LongFeat {
                    idx: f.idx,
                    name: truncate(&f.name, 40),
                },
            ));
        }
    }
    for (key, _w, fidx, _pos) in tiny.iter().filter(|_|!snap.ring) {
        if let Some(f) = snap.feats.iter().find(|f| f.idx == *fidx) {
            // 真实文本宽度（用于层避让）
            let w = p
                .layout_no_wrap(f.name.clone(), FontId::proportional(11.5), theme::TEXT)
                .size()
                .x
                + 14.0;
            items.push((
                *key,
                w,
                Lab::Feat {
                    idx: f.idx,
                    name: f.name.clone(),
                    color: if app.show_colors {
                        f.color
                    } else {
                        theme::FEATURE_NO_COLOR
                    },
                },
            ));
        }
    }
    for (name, pos, unique, hl) in snap.sites.iter().filter(|_|!snap.ring) {
        let sfont = FontId::proportional(11.5);
        let npos = nfmt(*pos);
        let w_name = p
            .layout_no_wrap(name.clone(), sfont.clone(), theme::SITE_UNIQUE)
            .size()
            .x;
        let w_pos = p
            .layout_no_wrap(format!("({npos})"), sfont, theme::TEXT)
            .size()
            .x;
        let tw = w_name + 4.0 + w_pos;
        {
            // 线性图谱：标签朝画布两侧展开（左右镜像），位置数字永远靠近引线一端；
            // 靠边放不下则换到另一侧。
            let tx = view.bp_x(*pos as f64);
            let (mut l, mut r) = if tx < rect.center().x {
                (tx - 6.0 - tw, tx - 6.0)
            } else {
                (tx + 6.0, tx + 6.0 + tw)
            };
            let mut side_left = tx < rect.center().x;
            if l < rect.left() + 2.0 {
                l = tx + 6.0;
                r = l + tw;
                side_left = false;
            } else if r > rect.right() - 2.0 {
                r = tx - 6.0;
                l = r - tw;
                side_left = true;
            }
            items.push((
                (l + r) * 0.5,
                tw + 12.0,
                Lab::Site {
                    name: name.clone(),
                    pos: *pos,
                    unique: *unique,
                    highlighted: *hl,
                    lin: Some((l, r, side_left)),
                },
            ));
        }
    }
    // Start / End 只在**真线性 DNA** 上写：环状分子按线性显示时两端是「⋯」，
    // RNA / 蛋白图谱写（oracle：lin1 环状不写、lin3 线性写、eco16s/THYX 都不写）。
    if !snap.ring && app.show_end_labels && snap.is_dna && !snap.circular {
        // 轴两端当作伪位点参与分层（SnapGene 写作「(0) Start」/「End (N)」）
        let w_start =
            p.layout_no_wrap(
                "(0) Start".to_string(),
                FontId::proportional(11.5),
                theme::TEXT,
            )
            .size()
            .x + 6.0;
        items.push((view.bp_x(0.0), w_start, Lab::EndMark { start: true }));
        let txt_end = format!("End ({})", nfmt(snap.len));
        let w_end = p
            .layout_no_wrap(txt_end, FontId::proportional(11.5), theme::TEXT)
            .size()
            .x
            + 6.0;
        items.push((
            view.bp_x(snap.len as f64),
            w_end,
            Lab::EndMark { start: false },
        ));
    }
    // 引物名也进同一套分层：原来单独画在酶位点那一行，几乎全被压住。
    if !snap.ring && app.show_primers {
        for primer in &snap.primers {
            let Some(&(a, b)) = primer.segs.first() else {
                continue;
            };
            // 与 SnapGene 一致：正向是「(起 .. 止) 名字」，反向是「名字 (起 .. 止)」。
            let label = if primer.segs.len() == 1 {
                if primer.forward {
                    format!("({} .. {}) {}", nfmt(a + 1), nfmt(b), primer.name)
                } else {
                    format!("{} ({} .. {})", primer.name, nfmt(a + 1), nfmt(b))
                }
            } else {
                primer.name.clone()
            };
            let tw = p
                .layout_no_wrap(label.clone(), FontId::proportional(11.5), theme::TEXT)
                .size()
                .x;
            // 两端引物的标签往视口里夹一下，避免被裁掉半截。
            let half = (tw + 10.0) * 0.5;
            let key = view
                .bp_x(((a + b) / 2) as f64)
                .clamp(rect.left() + half, rect.right() - half);
            items.push((
                key,
                tw + 10.0,
                Lab::Primer {
                    name: label,
                    forward: primer.forward,
                },
            ));
        }
    }
    // 线性视图下轴上写不进带里的特征名，也进分层（画在轴上方）；
    // 已排进「条下标签行」或已有小标记标签的不再参与。
    // RNA / 蛋白例外：放不进条、又占不到条下行的名字**不写**（oracle：eco16s 只保留
    // 排得下的二十来个名字；KPYK1 的小方块同样只保留十来个），只有蛋白的 `Site`
    // 类走轴上方纯文字。
    if !snap.ring && app.show_features {
        for (fi, f) in snap.feats.iter().enumerate() {
            if tiny.iter().any(|t| t.2 == f.idx)
                || tick_feats.iter().any(|t| t.2 == f.idx)
                || has_below(f.idx)
                || (!snap.is_dna && !tick_mask.get(fi).copied().unwrap_or(false))
            {
                continue;
            }
            let Some((s, e)) = f.segs.first() else {
                continue;
            };
            let (s, e) = (*s, *e);
            if e <= s {
                continue;
            }
            let tw = p
                .layout_no_wrap(f.name.clone(), FontId::proportional(11.5), theme::TEXT)
                .size()
                .x;
            let span_px = (e - s) as f32 * view.px_per_bp;
            if tw + 8.0 > span_px {
                let x = view.bp_x(((s + e) / 2) as f64);
                items.push((
                    x,
                    tw + 10.0,
                    Lab::LongFeat {
                        idx: f.idx,
                        name: truncate(&f.name, 40),
                    },
                ));
            }
        }
    }

    // 线性标签继续共用轴上方的严格行分配，不改变其布局规则。
    let mut order: Vec<usize> = (0..items.len()).collect();
    order.sort_by(|&i, &j| {
        items[i]
            .0
            .partial_cmp(&items[j].0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let pos_sorted: Vec<f32> = order
        .iter()
        .map(|&i| items[i].0)
        .collect();
    let width_sorted: Vec<f32> = order.iter().map(|&i| items[i].1).collect();
    // 放不下的标签明确隐藏，不强塞回已有行。
    let mut level_of = vec![0usize; items.len()];
    let mut hidden_labels = 0usize;
    if !snap.ring {
        let room = (((axis_y - rect.top() - 32.0) / LEVEL_DY).floor() as i32).clamp(1, 24) as usize;
        let lv = hx_render::assign_label_levels_strict(&pos_sorted, &width_sorted, room);
        for (k, &i) in order.iter().enumerate() {
            match lv[k] {
                Some(l) => level_of[i] = l,
                None => {
                    hidden_labels += 1;
                    level_of[i] = usize::MAX;
                }
            }
        }
    }
    let skip_label = |i: usize| level_of[i] == usize::MAX;

    for (i, (key, _w, lab)) in items.iter().enumerate() {
        if skip_label(i) {
            continue;
        }
        match lab {
            // 小特征：带特征色的圆角标签盒 + 细引线
            Lab::Feat { idx, name, color } => {
                let leader_from = {
                    let fi = snap.feats.iter().position(|f| f.idx == *idx).unwrap_or(0);
                    Pos2::new(*key, view_at(feat_level[fi]).axis_y - band * 0.42)
                };
                let picked = selected == Some(*idx);
                {
                    let y = axis_y - 6.0 - (level_of[i] as f32 + 0.5) * LEVEL_DY;
                    let r = draw_chip(
                        p,
                        name,
                        Pos2::new(*key, y),
                        Align2::CENTER_CENTER,
                        *color,
                        leader_from,
                        picked,
                    );
                    label_hits.push((r, *idx));
                }
            }
            // 引物名：紫色纯文字 + 细紫引线，指到骨架上的小箭头
            Lab::Primer { name, forward } => {
                let y = axis_y - 6.0 - (level_of[i] as f32 + 0.5) * LEVEL_DY;
                let from = Pos2::new(*key, if *forward { axis_y - 5.0 } else { axis_y + 6.0 });
                p.line_segment(
                    [Pos2::new(*key, y + 5.0), from],
                    Stroke::new(0.8, theme::with_alpha(theme::PRIMER, 130)),
                );
                p.text(
                    Pos2::new(*key, y),
                    Align2::CENTER_CENTER,
                    name,
                    FontId::proportional(11.5),
                    theme::PRIMER,
                );
            }
            // 酶位点：近黑文字（唯一者加粗、非唯一者灰）+ 细灰引线
            Lab::Site {
                name,
                pos,
                unique,
                highlighted,
                lin,
            } => {
                let col = if *highlighted {
                    theme::ACCENT
                } else if *unique {
                    theme::SITE_UNIQUE
                } else {
                    theme::SITE
                };
                let font = FontId::proportional(11.5);
                let npos = nfmt(*pos);
                {
                    let (l, r, side_left) = lin.unwrap_or((*key, *key, false));
                    let x = view.bp_x(*pos as f64);
                    let y = axis_y - 6.0 - (level_of[i] as f32 + 0.5) * LEVEL_DY;
                    // 酶切位点竖向小刻度（SnapGene 线性图谱：自轴上方向上）
                    p.line_segment(
                        [Pos2::new(x, axis_y - 4.0), Pos2::new(x, axis_y - 17.0)],
                        Stroke::new(1.0, theme::with_alpha(theme::BACKBONE, 170)),
                    );
                    // 引线：刻度顶端 → 标签近边（标签被挤开时是斜线，SnapGene 的扇形引线）
                    let near = if side_left {
                        Pos2::new(r, y + 2.0)
                    } else {
                        Pos2::new(l, y + 2.0)
                    };
                    p.line_segment(
                        [Pos2::new(x, axis_y - 17.0), near],
                        Stroke::new(0.8, theme::with_alpha(theme::TEXT_DIM, 120)),
                    );
                    // 位置数字靠近引线一侧（左右镜像）
                    let runs: [(String, bool); 2] = if side_left {
                        [(format!("({npos}) "), false), (name.clone(), *unique)]
                    } else {
                        [(name.clone(), *unique), (format!(" ({npos})"), false)]
                    };
                    draw_runs(
                        p,
                        &runs,
                        if side_left {
                            Pos2::new(r, y)
                        } else {
                            Pos2::new(l, y)
                        },
                        if side_left {
                            Align2::RIGHT_CENTER
                        } else {
                            Align2::LEFT_CENTER
                        },
                        &font,
                        col,
                    );
                }
            }
            // 线性视图里放不下的特征名（画在轴上方，分层避让；纯文字）
            Lab::LongFeat { idx, name } => {
                let y = axis_y - 6.0 - (level_of[i] as f32 + 0.5) * LEVEL_DY;
                // 标签必须留在视图内：蛋白图谱的 Site 可能在分子起点附近（如 `p3`），
                // 居中画会把前半截文字甩到内容区左边被裁掉（我们 e05 截图里
                // `domain interfaces [active]` 变成 `ain interfaces [active]`）。
                let key = if snap.ring {
                    *key
                } else {
                    let tw = p
                        .layout_no_wrap(name.clone(), FontId::proportional(11.5), theme::TEXT)
                        .size()
                        .x
                        + 10.0;
                    key.clamp(rect.left() + tw * 0.5 + 2.0, rect.right() - tw * 0.5 - 2.0)
                };
                let fi = snap.feats.iter().position(|f| f.idx == *idx).unwrap_or(0);
                let is_tick = tick_mask.get(fi).copied().unwrap_or(false);
                // Site 的引线接到刻度顶端（刻度和标签一起构成 SnapGene 的标记）
                let multi = snap
                    .feats
                    .get(fi)
                    .map(|f| f.segs.len() > 1)
                    .unwrap_or(false);
                let from_y = if is_tick && multi {
                    lane_y(feat_level[fi]) - band * 0.5
                } else if is_tick {
                    lane_y(level_count - 1) - band * 1.05
                } else {
                    view_at(feat_level[fi]).axis_y - band * 0.42
                };
                p.line_segment(
                    [Pos2::new(key, from_y), Pos2::new(key, y + 7.5)],
                    Stroke::new(0.8, theme::with_alpha(theme::TEXT_DIM, 120)),
                );
                p.text(
                    Pos2::new(key, y),
                    Align2::CENTER_CENTER,
                    name,
                    FontId::proportional(11.5),
                    theme::TEXT,
                );
                label_hits.push((
                    text_rect(
                        p,
                        name,
                        FontId::proportional(11.5),
                        Pos2::new(key, y),
                        Align2::CENTER_CENTER,
                    ),
                    *idx,
                ));
            }
            // 线性图谱两端：Start 左对齐到起点、End 右向延伸到终点（无引线）
            Lab::EndMark { start } => {
                let y = axis_y - 6.0 - (level_of[i] as f32 + 0.5) * LEVEL_DY;
                let font = FontId::proportional(11.5);
                if *start {
                    draw_runs(
                        p,
                        &[("(0) ".to_string(), false), ("Start".to_string(), true)],
                        Pos2::new(*key, y),
                        Align2::RIGHT_CENTER,
                        &font,
                        theme::TEXT,
                    );
                } else {
                    draw_runs(
                        p,
                        &[
                            ("End".to_string(), true),
                            (format!(" ({})", nfmt_plain(snap.len)), false),
                        ],
                        Pos2::new(*key, y),
                        Align2::LEFT_CENTER,
                        &font,
                        theme::TEXT,
                    );
                }
            }
        }
    }

    if let Some(labels)=&ring_labels {
        for (i,kind) in labels.kinds.iter().enumerate() {
            let Some(callout)=labels.geometry.labels[i] else{hidden_labels+=1;continue};
            let anchor=Pos2::from(callout.port);
            let align=if callout.right {Align2::LEFT_CENTER}else{Align2::RIGHT_CENTER};
            let (source_radius,color)=match kind {
                RingLabelKind::Feature(fi)=>(view_at(feat_level[*fi]).radius+band*0.46,theme::with_alpha(theme::TEXT_DIM,150)),
                RingLabelKind::Site(_)=>(r_out+4.0,theme::with_alpha(theme::TEXT_DIM,150)),
                RingLabelKind::Primer(pi)=>(r_out+if snap.primers[*pi].forward {12.0}else{25.0},
                    if app.show_colors {theme::PRIMER}else{theme::TEXT_DIM}),
            };
            let points=labels.geometry.leader(&labels.inputs[i],&callout,source_radius);
            p.add(Shape::line(points.into_iter().map(Pos2::from).collect(),Stroke::new(0.8,color)));
            let box_=Rect::from_min_max(Pos2::new(callout.rect[0],callout.rect[1]),Pos2::new(callout.rect[2],callout.rect[3]));
            match kind {
                RingLabelKind::Feature(fi)=>{
                    let f=&snap.feats[*fi];let picked=selected==Some(f.idx);
                    let color=if app.show_colors {f.color}else{theme::FEATURE_NO_COLOR};
                    p.rect_filled(box_,2.0,if picked {theme::ACCENT}else{theme::mix(color,Color32::WHITE,0.68)});
                    p.rect_stroke(box_,2.0,Stroke::new(0.7,if picked {theme::ACCENT}else{theme::darken(color,0.30)}),egui::StrokeKind::Inside);
                    p.text(box_.center(),Align2::CENTER_CENTER,&f.name,map_font(),if picked {Color32::WHITE}else{theme::TEXT});
                    label_hits.push((box_,f.idx));
                },
                RingLabelKind::Site(si)=>{
                    let (name,pos,unique,highlighted)=&snap.sites[*si];
                    let runs=if callout.right {[(name.clone(),*unique),(format!(" ({})",nfmt(*pos)),false)]}
                        else{[(format!("({}) ",nfmt(*pos)),false),(name.clone(),*unique)]};
                    draw_runs(p,&runs,anchor,align,&map_font(),if *highlighted {theme::ACCENT}else if *unique {theme::SITE_UNIQUE}else{theme::SITE});
                },
                RingLabelKind::Primer(pi)=>{
                    let primer=&snap.primers[*pi];
                    p.text(anchor,align,primer_label(primer,snap.len),map_font(),if app.show_colors {theme::PRIMER}else{theme::TEXT_DIM});
                    primer_hits.push((box_,primer.idx));
                },
            }
        }
    }

    // ---- 特征名（贴弧旋转写在环内 / 带里） ----
    let i_sel = selected;
    {
        let font = map_font();
        for (fi, f) in snap.feats.iter().enumerate() {
            // 名字跟着特征所在的那一层走
            let fv = view_at(feat_level.get(fi).copied().unwrap_or(0));
            let view = &fv;
            if (snap.ring && ring_labels.as_ref().is_some_and(|l|l.kinds.iter().any(|k|matches!(k,RingLabelKind::Feature(i) if *i==fi))))
                || (!snap.ring && tiny.iter().any(|t| t.2 == f.idx)) || tick_feats.iter().any(|t| t.2 == f.idx) {
                continue;
            }
            let segment=if snap.ring {f.segs.iter().max_by_key(|(a,b)|b-a)}else{f.segs.first()};
            let Some((s, e)) = segment else {
                continue;
            };
            let (s, e) = (*s, *e);
            if e <= s {
                continue;
            }
            let mid = (s + e) / 2;
            let tw = p
                .layout_no_wrap(f.name.clone(), font.clone(), theme::TEXT)
                .size()
                .x;
            let feature_color = if app.show_colors {
                f.color
            } else {
                theme::FEATURE_NO_COLOR
            };
            let col = theme::TEXT;
            // 带内名字：SnapGene 用中性灰（暗底转白），带色只负责填充
            let light_band = theme::luma(theme::band_fill(feature_color)) > 0.55;
            let in_band_col = if light_band {
                theme::TEXT
            } else {
                Color32::WHITE
            };
            if snap.ring {
                let ang = view.bp_angle(snap.len, mid as f64);
                let arc_px = (e - s) as f32 / snap.len as f32 * std::f32::consts::TAU * view.radius;
                let circ = std::f32::consts::TAU * view.radius;
                if tw + 8.0 <= arc_px && band >= 12.0 {
                    if Some(f.idx) == i_sel {
                        // 选中的特征：SnapGene 把名字画成**蓝底白字圆角框**（“rop”那样）
                        let center = pt(c, ang, view.radius);
                        let g = p.layout_no_wrap(f.name.clone(), font.clone(), Color32::WHITE);
                        let sz = g.size() + Vec2::new(12.0, 5.0);
                        let r = readable_angle(ang);
                        let (sa, ca) = (r.sin(), r.cos());
                        let rot = |v: Vec2| -> Pos2 {
                            Pos2::new(
                                center.x + v.x * ca - v.y * sa,
                                center.y + v.x * sa + v.y * ca,
                            )
                        };
                        let h = sz * 0.5;
                        p.add(Shape::convex_polygon(
                            vec![
                                rot(Vec2::new(-h.x, -h.y)),
                                rot(Vec2::new(h.x, -h.y)),
                                rot(Vec2::new(h.x, h.y)),
                                rot(Vec2::new(-h.x, h.y)),
                            ],
                            theme::ACCENT,
                            Stroke::new(0.8, theme::darken(theme::ACCENT, 0.25)),
                        ));
                        text_rot(
                            p,
                            &f.name,
                            center,
                            r,
                            Align2::CENTER_CENTER,
                            font.clone(),
                            Color32::WHITE,
                            false,
                        );
                    } else {
                        text_rot(
                            p,
                            &f.name,
                            pt(c, ang, view.radius),
                            readable_angle(ang),
                            Align2::CENTER_CENTER,
                            font.clone(),
                            in_band_col,
                            false,
                        );
                    }
                } else if tw <= circ * 0.3 {
                    // 沿弧写在环内
                    let anchor = pt(c, ang, view.radius - band - 10.0);
                    label_hits.push((
                        rot_text_rect(
                            p,
                            &f.name,
                            &font,
                            anchor,
                            readable_angle(ang),
                            Align2::CENTER_CENTER,
                        ),
                        f.idx,
                    ));
                    text_rot(
                        p,
                        &f.name,
                        anchor,
                        readable_angle(ang),
                        Align2::CENTER_CENTER,
                        font.clone(),
                        col,
                        true,
                    );
                } else {
                    // 太长：水平放在弧中点内侧
                    let t = truncate(&f.name, 34);
                    let inside = pt(c, ang, view.radius - band - 18.0);
                    let right = ang.sin() >= 0.0;
                    let anchor = inside;
                    p.line_segment(
                        [
                            pt(c, ang, view.radius - band - 3.0),
                            anchor
                                + if right {
                                    Vec2::new(6.0, 0.0)
                                } else {
                                    Vec2::new(-6.0, 0.0)
                                },
                        ],
                        Stroke::new(0.8, theme::with_alpha(theme::TEXT_DIM, 110)),
                    );
                    text_rot(
                        p,
                        &t,
                        anchor,
                        0.0,
                        if right {
                            Align2::LEFT_CENTER
                        } else {
                            Align2::RIGHT_CENTER
                        },
                        font.clone(),
                        col,
                        true,
                    );
                }
            } else {
                let x = view.bp_x(mid as f64);
                let span_px = (e - s) as f32 * view.px_per_bp;
                if tw + 8.0 <= span_px {
                    // 带内名字（SnapGene 的线性图谱：灰字、半粗）
                    for d in [Vec2::new(0.45, 0.0), Vec2::ZERO] {
                        p.text(
                            Pos2::new(x, view.axis_y) + d,
                            Align2::CENTER_CENTER,
                            &f.name,
                            font.clone(),
                            in_band_col,
                        );
                    }
                } else {
                    // 放不进带里：已由 `Lab::LongFeat` 在轴上分层绘制
                }
            }
        }
    }

    // 标签放不下时在图角提示（SnapGene: “N enzymes are hidden”）
    if hidden_labels > 0 {
        p.text(
            Pos2::new(rect.left() + 8.0, lin_bottom - 6.0),
            Align2::LEFT_BOTTOM,
            format!("{hidden_labels} label(s) hidden (crowded)"),
            FontId::proportional(10.5),
            theme::TEXT_DIM,
        );
    }

    // ---- 线性图谱的横向平移/缩放条（SnapGene 图 2 底部那条）----
    if !snap.ring && snap.len > 0 && !app.cloning_canvas {
        let len = snap.len;
        let bar = Rect::from_min_size(
            Pos2::new(rect.left() + 2.0, rect.bottom() - PAN_BAR_H),
            Vec2::new((rect.width() - 4.0).max(20.0), PAN_BAR_H),
        );
        p.rect_filled(bar, 0.0, theme::PANEL_ALT);
        p.line_segment(
            [bar.left_top(), bar.right_top()],
            Stroke::new(0.8, theme::BORDER),
        );
        let ox0 = bar.left() + 40.0;
        let ox1 = bar.right() - 40.0;
        let ow = (ox1 - ox0).max(20.0);
        let xb = |bp: f64| ox0 + (bp / len as f64) as f32 * ow;
        let px = view.px_per_bp.max(1e-6);
        let vis0 = (pan_bp as f32 - 40.0 / px).max(0.0);
        let vis1 = (pan_bp as f32 + (rect.width() - 40.0) / px).min(len as f32);
        // 刻度 + 数字
        for (bp, major) in hx_render::map::ruler_ticks_subdiv(len, ow / len as f32, 90.0, 5) {
            let b = bp as usize;
            if b > len {
                break;
            }
            let x = xb(bp);
            let h = if major { 6.0 } else { 3.0 };
            p.line_segment(
                [
                    Pos2::new(x, bar.top() + 3.0),
                    Pos2::new(x, bar.top() + 3.0 + h),
                ],
                Stroke::new(0.8, theme::TICK_MAJOR),
            );
            if major && b > 0 {
                p.text(
                    Pos2::new(x, bar.top() + 10.0),
                    Align2::CENTER_TOP,
                    nfmt_plain(b),
                    FontId::proportional(9.0),
                    theme::RULER_TEXT,
                );
            }
        }
        // 平移轨道 + 可视范围
        p.line_segment(
            [
                Pos2::new(ox0, bar.bottom() - 6.0),
                Pos2::new(ox1, bar.bottom() - 6.0),
            ],
            Stroke::new(1.0, theme::with_alpha(theme::BORDER, 200)),
        );
        let vrect = Rect::from_min_max(
            Pos2::new(xb(vis0 as f64), bar.top() + 16.0),
            Pos2::new(
                xb(vis1 as f64).max(xb(vis0 as f64) + 14.0),
                bar.bottom() - 2.0,
            ),
        );
        p.rect_filled(vrect, 3.0, theme::PANEL);
        p.rect_stroke(
            vrect,
            3.0,
            Stroke::new(1.0, theme::with_alpha(theme::TEXT_DIM, 200)),
            egui::StrokeKind::Inside,
        );
        for i in 0..3 {
            let x = vrect.center().x - 3.0 + i as f32 * 3.0;
            p.line_segment(
                [
                    Pos2::new(x, vrect.top() + 2.0),
                    Pos2::new(x, vrect.bottom() - 2.0),
                ],
                Stroke::new(0.7, theme::with_alpha(theme::TEXT_DIM, 160)),
            );
        }
        // 拖拽平移
        if let Some(pos) = p.ctx().pointer_interact_pos() {
            let hit = Rect::from_min_max(Pos2::new(ox0, bar.top()), Pos2::new(ox1, bar.bottom()));
            if hit.contains(pos) && p.ctx().input(|i| i.pointer.any_down()) {
                let f = ((pos.x - ox0) / ow).clamp(0.0, 1.0) as f64;
                let half = (rect.width() - 40.0) as f64 * 0.5 / px as f64;
                let target = f * len as f64 - half;
                let max_pan = (len as f64 - (rect.width() - 40.0) as f64 / px as f64).max(0.0);
                pan_bp = target.clamp(0.0, max_pan);
                app.pan_bp = pan_bp;
            }
        }
    }

    // ---- 选段 / 悬停 ----
    if let Some((s, e)) = app.selection {
        let digest=app.digest_selection.as_ref().filter(|d|d.range==(s,e));
        if let Some(d)=digest {
            // Paint each strand's selected span independently. The two half-bands
            // meet in a staircase, never a fabricated common radial/blunt edge.
            let color=theme::with_alpha(theme::SELECTION,80);
            if snap.ring {
                let mid=r_out+1.0;
                for (start,end,inner,outer) in [(d.top_start,d.top_end,mid,r_out+10.0),(d.bottom_start,d.bottom_end,r_out-8.0,mid)] {
                    let a=view.bp_angle(snap.len,start as f64);let sweep=(end-start) as f32/snap.len as f32*std::f32::consts::TAU;
                    if sweep>0.0{p.add(Shape::mesh(selection_mesh(c,a,sweep,inner,outer)));}
                }
                for (top,bottom) in [(d.top_start,d.bottom_start),(d.top_end,d.bottom_end)]{
                    let at=view.bp_angle(snap.len,top as f64);let ab=view.bp_angle(snap.len,bottom as f64);
                    p.line_segment([pt(c,at,r_out+10.0),pt(c,at,mid)],Stroke::new(1.2,theme::with_alpha(theme::SELECTION,170)));
                    let sweep=(bottom-top) as f32/snap.len as f32*std::f32::consts::TAU;let steps=((sweep.abs()*mid/2.0).ceil() as usize).max(1);
                    let points=(0..=steps).map(|i|pt(c,at+sweep*i as f32/steps as f32,mid)).collect();p.add(Shape::line(points,Stroke::new(1.2,theme::with_alpha(theme::SELECTION,170))));
                    p.line_segment([pt(c,ab,mid),pt(c,ab,r_out-8.0)],Stroke::new(1.2,theme::with_alpha(theme::SELECTION,170)));
                }
            }else{
                let mid=view.axis_y;
                for(start,end,y0,y1)in [(d.top_start,d.top_end,mid-12.0,mid),(d.bottom_start,d.bottom_end,mid,mid+12.0)]{
                    // Split a circular interval when its source is displayed linearly.
                    let n=snap.len as i64;let length=end-start;let start=start.rem_euclid(n)as usize;let end=if length>=n{snap.len}else{(start as i64+length).rem_euclid(n)as usize};let start=if length>=n{0}else{start};
                    for(a,b)in hx_render::split_wrap(start,end,snap.len){p.rect_filled(Rect::from_min_max(Pos2::new(view.bp_x(a as f64),y0),Pos2::new(view.bp_x(b as f64),y1)),0.0,color);}
                }
            }
        } else if snap.ring {
            let a0 = view.bp_angle(snap.len, s as f64);
            let sweep = span_len(s, e, snap.len) as f32 / snap.len as f32 * std::f32::consts::TAU;
            if sweep > 0.0 {
                p.add(Shape::mesh(selection_mesh(
                    c,
                    a0,
                    sweep,
                    r_out - 8.0,
                    r_out + 10.0,
                )));
                if span_len(s, e, snap.len) < snap.len {
                    let inner = app
                        .selected_feature
                        .and_then(|idx| snap.feats.iter().position(|f| f.idx == idx))
                        .map(|i| view_at(feat_level[i]).radius - band * 0.7)
                        .unwrap_or(view.radius - band);
                    for angle in [a0, a0 + sweep] {
                        p.line_segment(
                            [pt(c, angle, inner), pt(c, angle, r_out + 10.0)],
                            Stroke::new(1.4, theme::with_alpha(theme::SELECTION, 170)),
                        );
                    }
                }
            }
        } else {
            for (s, e) in hx_render::split_wrap(s, e, snap.len) {
                let x0 = view.bp_x(s as f64);
                let x1 = view.bp_x(e as f64);
                p.rect_filled(
                    egui::Rect::from_min_max(
                        Pos2::new(x0, view.axis_y - band * 0.5 - 9.0),
                        Pos2::new(
                            x1,
                            view_at(level_count.saturating_sub(1)).axis_y + band * 0.5 + 9.0,
                        ),
                    ),
                    1.0,
                    theme::with_alpha(theme::SELECTION, 80),
                );
            }
        }
    }
    if let (Some(hbp), true) = (app.hover_bp, snap.ring) {
        if hbp < snap.len {
            let a = view.bp_angle(snap.len, hbp as f64);
            p.line_segment(
                [pt(c, a, r_out + 3.0), pt(c, a, r_out - band - 14.0)],
                Stroke::new(0.9, theme::with_alpha(theme::SELECTION, 150)),
            );
        }
    }

    // ---- 圆心名称 + 长度 ----
    if snap.ring {
        let title=truncate(&snap.name,36);
        // SG map title is emphasized; use the same small-font weight technique
        // as enzyme names without replacing the user's native title.
        p.text(c+Vec2::new(0.35,-10.0),Align2::CENTER_CENTER,&title,FontId::proportional(14.0),theme::TEXT);
        p.text(
            c + Vec2::new(0.0, -10.0),
            Align2::CENTER_CENTER,
            title,
            FontId::proportional(14.0),
            theme::TEXT,
        );
        p.text(
            c + Vec2::new(0.0, 6.0),
            Align2::CENTER_CENTER,
            format!("{} {}", nfmt(snap.len), snap.unit),
            FontId::proportional(11.5),
            theme::TEXT,
        );
    } else {
        // 线性图谱：分子名 + 长度写在画布**底部中央**（SnapGene 的写法：粗体名一行，
        // 灰色「N bp」一行；oracle：lin3 的 T-Vector pMD20 / 2737 bp、eco16s / 1541 bases）。
        // 对话框里名字紧贴轴下方（SnapGene 的对话框），主窗口仍放在画布底部。
        let (cx, y) = if app.cloning_canvas {
            (rect.center().x, (axis_y + 118.0).min(rect.bottom() - 30.0))
        } else {
            (rect.center().x, rect.bottom() - 30.0)
        };
        p.text(
            Pos2::new(cx, y),
            Align2::CENTER_CENTER,
            truncate(&snap.name, 60),
            FontId::proportional(13.0),
            theme::TEXT,
        );
        p.text(
            Pos2::new(cx, y + 16.0),
            Align2::CENTER_CENTER,
            format!("{} {}", nfmt(snap.len), snap.unit),
            FontId::proportional(11.5),
            theme::TEXT_DIM,
        );
    }

    // ---- 悬停提示 ----
    if let (Some(i), Some(pos)) = (hover, resp.hover_pos()) {
        if let Some(f) = snap.feats.iter().find(|f| f.idx == i) {
            if let Some((s, e)) = f.segs.first() {
                let label = format!(
                    "{}  {}..{} ({} {})",
                    f.name,
                    nfmt(s + 1),
                    nfmt(*e),
                    nfmt(span_len(*s, *e, snap.len)),
                    snap.unit
                );
                let font = FontId::proportional(11.5);
                let size = p
                    .layout_no_wrap(label.clone(), font.clone(), theme::TEXT)
                    .size();
                let box_pos = pos + Vec2::new(14.0, 14.0);
                let r = egui::Rect::from_min_size(box_pos, size + Vec2::new(10.0, 6.0));
                p.rect_filled(r, 3.0, Color32::from_rgba_unmultiplied(255, 255, 255, 236));
                p.rect_stroke(
                    r,
                    3.0,
                    Stroke::new(1.0, theme::BORDER),
                    egui::StrokeKind::Inside,
                );
                p.text(
                    box_pos + Vec2::new(5.0, 3.0),
                    Align2::LEFT_TOP,
                    label,
                    font,
                    theme::TEXT,
                );
            }
        }
    }

    // ---- 右下角：当前酶显示集合（SnapGene 的 `Unique 6+ Cutters | Nonredundant`）----
    // 只有 DNA 文档有酶切位点；RNA/蛋白画这句会误导。
    if snap.is_dna {
        let label = crate::i18n::text(ui.ctx(),app.enzyme_preset.label());
        let anchor = if snap.ring {
            Pos2::new(rect.right() - 12.0, rect.bottom() - 12.0)
        } else {
            // 线性图谱底部有一条平移/缩放条，标签放在它上面
            Pos2::new(rect.right() - 12.0, rect.bottom() - PAN_BAR_H - 8.0)
        };
        p.text(
            anchor,
            Align2::RIGHT_BOTTOM,
            label,
            FontId::proportional(11.5),
            theme::TEXT_DIM,
        );
    }

    // ---- 写回状态 ----
    app.rotation = rotation;
    app.zoom = zoom;
    app.pan_bp = pan_bp;
    // 标签命中区留给下一帧的点击处理（见 `label_hits` 的说明）
    ui.ctx()
        .data_mut(|d| d.insert_temp(label_hits_id, label_hits));
    ui.ctx().data_mut(|d|d.insert_temp(primer_hits_id,primer_hits));

    rect
}



/// 名称截断（超长补 `…`）。
fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() > max {
        s.chars().take(max.saturating_sub(1)).collect::<String>() + "…"
    } else {
        s.to_string()
    }
}

/// 线性特征的连续轮廓，所有尖端均在其真实坐标内。
fn linear_feature_polygon(
    x0: f32,
    x1: f32,
    y: f32,
    half_height: f32,
    forward_arrow: bool,
    reverse_arrow: bool,
) -> Vec<Pos2> {
    let span = (x1 - x0).max(0.0);
    let heads = usize::from(forward_arrow) + usize::from(reverse_arrow);
    let tip = (span * 0.35)
        .clamp(1.0, TIP_PX)
        .min(span / heads.max(1) as f32);
    let body_start = x0 + if reverse_arrow { tip } else { 0.0 };
    let body_end = x1 - if forward_arrow { tip } else { 0.0 };
    let mut points = vec![
        Pos2::new(body_start, y - half_height),
        Pos2::new(body_end, y - half_height),
    ];
    if forward_arrow {
        points.push(Pos2::new(x1, y));
    }
    points.push(Pos2::new(body_end, y + half_height));
    points.push(Pos2::new(body_start, y + half_height));
    if reverse_arrow {
        points.push(Pos2::new(x0, y));
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rep1_dense_ring_keeps_all_callouts_and_primer_binding_labels_when_space_allows() {
        let Some(path)=std::env::var_os("HX_RING_LABEL_FIXTURE") else{
            eprintln!("Skipping optional local Rep1 fixture: set HX_RING_LABEL_FIXTURE");return;
        };
        let doc=crate::doc::Document::from_bytes(std::fs::read(path).unwrap()).unwrap();
        let ctx=egui::Context::default();
        let mut app=HelixApp::preview(&ctx,doc);app.show_primers=true;app.recompute_sites();
        let snap=snapshot(&app).unwrap();
        let (lanes,count)=interval_lanes(&snap.feats.iter().map(|f|f.segs.clone()).collect::<Vec<_>>());
        let (_,orfs)=linear_orf_lanes(&snap.orfs,snap.len);
        let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(Rect::from_min_size(Pos2::ZERO,Vec2::new(2000.0,1200.0))),..Default::default()},|ui|{
            for max in [Pos2::new(1982.0,1140.0),Pos2::new(1782.0,1040.0)] {
                let labels=ring_labels(&snap,ui.painter(),Rect::from_min_max(Pos2::new(346.0,168.0),max),0.0,1.0,&lanes,count,orfs);
                assert!(labels.geometry.labels.iter().all(Option::is_some));
                assert_eq!(labels.kinds.iter().filter(|k|matches!(k,RingLabelKind::Primer(_))).count(),snap.primers.len());
            }
        });output.textures_delta.clear();
    }

    fn callout_document()->crate::doc::Document {
        let bytes=hx_formats::write::NewFile::dna(b"A".repeat(9332)).to_bytes();
        let mut file=hx_formats::DnaFile::parse(&bytes).unwrap();
        file.blocks.iter_mut().find(|b|b.block_type==0).unwrap().payload[0]|=1;
        file.blocks.iter_mut().find(|b|b.block_type==10).unwrap().payload=br##"<Features><Feature name="6xHis" directionality="1" type="CDS"><Segment range="9240-9257" color="#a04070"/></Feature></Features>"##.to_vec();
        let primers=br#"<Primers><Primer name="Primer 1" sequence="AAAAAAAAAAAAAAAAAAAA"><BindingSite simplified="1" location="20-39" boundStrand="0" annealedBases="AAAAAAAAAAAAAAAAAAAA"/></Primer></Primers>"#.to_vec();
        if let Some(block)=file.blocks.iter_mut().find(|b|b.block_type==5) {block.payload=primers;}
        else {file.blocks.push(hx_formats::Block{block_type:5,offset:0,payload:primers});}
        crate::doc::Document::from_bytes(file.to_bytes()).unwrap()
    }

    #[test]
    fn short_primer_and_feature_callouts_are_drawn_and_clickable() {
        let ctx=egui::Context::default();
        let mut app=HelixApp::preview(&ctx,callout_document());app.show_translation=false;app.show_sites=false;
        let mut targets=Vec::new();
        let rect=Rect::from_min_size(Pos2::ZERO,Vec2::new(1800.0,1100.0));
        let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(rect),..Default::default()},|ui|{
            let snap=snapshot(&app).unwrap();
            assert_eq!(snap.primers.len(),1);
            let labels=ring_labels(&snap,ui.painter(),Rect::from_min_max(rect.min+Vec2::new(46.0,62.0),rect.max-Vec2::new(18.0,34.0)),0.0,1.0,&[0],1,0);
            for (i,k) in labels.kinds.iter().enumerate() {
                let c=labels.geometry.labels[i].unwrap();
                targets.push((matches!(k,RingLabelKind::Primer(_)),Pos2::new((c.rect[0]+c.rect[2])*0.5,(c.rect[1]+c.rect[3])*0.5)));
            }
            draw(ui,&mut app);
        });
        let text=output.shapes.iter().filter_map(|s|if let Shape::Text(t)=&s.shape {Some(t.galley.job.text.as_str())}else{None}).collect::<Vec<_>>().join("\n");
        assert!(text.contains("Primer 1 (21 .. 40)"));assert!(text.contains("6xHis"));output.textures_delta.clear();
        assert_eq!(targets.len(),2);
        for (primer,pos) in targets {
            for pressed in [true,false] {
                let mut output=ctx.run_ui(egui::RawInput{screen_rect:Some(rect),events:vec![egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton{pos,button:egui::PointerButton::Primary,pressed,modifiers:egui::Modifiers::NONE}],..Default::default()},|ui|{draw(ui,&mut app);});
                output.textures_delta.clear();
            }
            if primer {assert_eq!(app.selected_primer,Some(0));}else{assert_eq!(app.selected_feature,Some(0));}
        }
    }

    #[test]
    fn cloned_product_map_hides_low_tm_additional_primer_binding_labels() {
        let Ok(bytes)=std::fs::read("/tmp/opencode/oracle/Cloned.dna") else {
            eprintln!("skip: saved SnapGene cloning fixture /tmp/opencode/oracle/Cloned.dna is not present");
            return;
        };
        let doc=crate::doc::Document::from_bytes(bytes).unwrap();
        let ctx=egui::Context::default();let mut app=HelixApp::preview(&ctx,doc);app.show_primers=true;
        let snap=snapshot(&app).unwrap();
        assert_eq!(snap.primers.len(),2,"only the >=40°C primary binding sites should display");
        assert_eq!(snap.primers[0].name,"Primer 1");
        assert_eq!(snap.primers[1].name,"Primer 2");
        assert_eq!(snap.primers[0].segs,vec![(294,314)]);
        assert_eq!(snap.primers[1].segs,vec![(924,944)]);
    }

    #[test]
    fn all_external_kinds_use_real_measured_boxes_in_dense_layout() {
        let ctx=egui::Context::default();
        let app=HelixApp::preview(&ctx,callout_document());
        let mut snap=snapshot(&app).unwrap();
        snap.sites=(0..42).map(|i|(format!("Restriction site {i}"),8800+i*11,true,false)).collect();
        snap.primers=(0..20).map(|i|PrimerRow{idx:i,name:format!("Primer {}",i%2+1),segs:vec![(i*200,i*200+28)],forward:i%2==0}).collect();
        let mut output=ctx.run_ui(egui::RawInput::default(),|ui|{
            for rotation in [0.0,2.0,3.14] {
                for size in [Vec2::new(1800.0,1100.0),Vec2::new(1100.0,740.0)] {
                    let bounds=Rect::from_min_size(Pos2::ZERO,size);
                    let labels=ring_labels(&snap,ui.painter(),bounds,rotation,1.0,&[0],1,0);
                    assert_eq!(labels.kinds.len(),63);
                    if size.y>1000.0 {assert!(labels.geometry.labels.iter().all(Option::is_some));}
                    let boxes=labels.geometry.labels.iter().flatten().map(|l|Rect::from_min_max(Pos2::new(l.rect[0],l.rect[1]),Pos2::new(l.rect[2],l.rect[3]))).collect::<Vec<_>>();
                    for (i,a) in boxes.iter().enumerate() {
                        assert!(bounds.expand(0.01).contains_rect(*a));
                        for b in &boxes[..i] {assert!(!a.intersects(*b));}
                    }
                }
            }
        });output.textures_delta.clear();
    }

    #[test]
    fn ring_drag_follows_both_directions_across_origin() {
        let mut forward = RingDrag::new(950.0);
        assert_eq!(forward.update(20.0, 1000), Some((950, 20)));
        assert_eq!(forward.update(80.0, 1000), Some((950, 80)));
        let mut reverse = RingDrag::new(50.0);
        assert_eq!(reverse.update(980.0, 1000), Some((980, 50)));
        assert_eq!(reverse.update(900.0, 1000), Some((900, 50)));
        assert_eq!(reverse.update(0.0, 1000), Some((0, 50)));
    }

    #[test]
    fn ring_drag_can_select_more_than_half_a_circle_and_a_full_circle() {
        let mut drag = RingDrag::new(100.0);
        for bp in [300.0, 490.0, 700.0] {
            drag.update(bp, 1000);
        }
        assert_eq!(drag.update(900.0, 1000), Some((100, 900)));
        assert_eq!(drag.update(50.0, 1000), Some((100, 50)));
        assert_eq!(drag.update(100.0, 1000), Some((0, 1000)));
    }

    #[test]
    fn multi_segment_feature_selection_keeps_the_origin_crossing() {
        let spans = [(900, 1000), (0, 100)];
        assert_eq!(feature_range(&spans, 1000, true), Some((900, 100)));
        assert_eq!(feature_range(&spans, 1000, false), Some((0, 1000)));
        assert_eq!(
            feature_range(&[(10, 50), (40, 100), (80, 200)], 1000, true),
            Some((10, 200))
        );
        assert_eq!(feature_range(&[(0, 1000)], 1000, true), Some((0, 1000)));
    }

    #[test]
    fn ring_tracks_separate_orfs_features_and_nested_feature_rows() {
        for radius in [120.0, 410.0] {
            for rows in [2, 7] {
                let t = RingTracks::new(radius, 18.0, rows, 4);
                let deepest = t.feature_radius - (rows - 1) as f32 * t.feature_step - t.band * 0.46;
                assert!(deepest >= 27.9);
                assert!(t.feature_step > t.band);
                assert!(t.orf_radius - 3.0 * t.orf_step > t.feature_radius + t.band * 0.46);
            }
        }
    }

    #[test]
    fn large_selection_mesh_never_fills_the_center_of_the_map() {
        for sweep in [0.2, 4.8, std::f32::consts::TAU] {
            let mesh = selection_mesh(Pos2::ZERO, 5.8, sweep, 300.0, 318.0);
            assert!(mesh.is_valid());
            for triangle in mesh.indices.chunks_exact(3) {
                let centroid = triangle
                    .iter()
                    .map(|i| mesh.vertices[*i as usize].pos.to_vec2())
                    .fold(Vec2::ZERO, |sum, point| sum + point)
                    / 3.0;
                assert!(centroid.length() >= 299.9);
                assert!(centroid.length() <= 318.0);
            }
        }
    }

    #[test]
    fn duplicate_features_get_distinct_rows_but_adjacent_features_share_a_row() {
        // 模拟文件中的乱序注释：两份 KanR、紧挨着的 T7 promoter / lac operator。
        let spans = vec![
            vec![(560, 1375)],
            vec![(5002, 5027)],
            vec![(560, 1375)],
            vec![(4983, 5002)],
        ];
        let (rows, count) = interval_lanes(&spans);
        assert_eq!(count, 2);
        assert_ne!(rows[0], rows[2]);
        assert_eq!(rows[1], rows[3]);
    }

    #[test]
    fn feature_lanes_keep_wrapped_segments_together_without_a_row_cap() {
        let mut spans = vec![vec![(900, 1000), (0, 100)], vec![(100, 900)]];
        spans.extend((0..6).map(|_| vec![(950, 980)]));
        let (rows, count) = interval_lanes(&spans);
        assert_eq!(rows[0], rows[1]);
        assert_eq!(count, 7);
        assert_eq!(
            rows[2..]
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            6
        );
    }

    #[test]
    fn linear_orfs_pack_overlaps_and_reserve_separate_strand_groups() {
        let orfs = vec![
            (100, 800, 1, true, 0),
            (500, 900, 1, true, 1),
            (200, 700, -1, true, 0),
            (1000, 1300, -1, true, 2),
            (1400, 100, -1, true, 1),
        ];
        let (rows, count) = linear_orf_lanes(&orfs, 1500);
        assert_eq!(count, 3);
        assert_ne!(rows[0], rows[1]);
        assert_eq!(rows[2], rows[3]);
        assert_eq!(rows[2], rows[4]);
        assert!(rows[..2].iter().all(|r| *r < rows[2]));
    }

    #[test]
    fn wrapped_linear_orfs_have_a_head_only_at_the_true_translation_end() {
        assert_eq!(
            linear_orf_segments(900, 100, 1, 1000),
            vec![(900, 1000, false), (0, 100, true)]
        );
        assert_eq!(
            linear_orf_segments(900, 100, -1, 1000),
            vec![(900, 1000, true), (0, 100, false)]
        );
        assert_eq!(
            linear_orf_segments(100, 900, 1, 1000),
            vec![(100, 900, true)]
        );
        assert_eq!(
            linear_orf_segments(100, 900, -1, 1000),
            vec![(100, 900, true)]
        );
    }

    #[test]
    fn undirected_feature_has_flat_ends_at_its_coordinates() {
        let points = linear_feature_polygon(10.0, 30.0, 20.0, 4.0, false, false);
        assert_eq!(
            points,
            vec![
                Pos2::new(10.0, 16.0),
                Pos2::new(30.0, 16.0),
                Pos2::new(30.0, 24.0),
                Pos2::new(10.0, 24.0),
            ]
        );
    }

    #[test]
    fn directional_tips_stay_inside_the_feature_including_short_features() {
        for span in [0.2, 2.0, 20.0] {
            for (forward, reverse) in [(true, false), (false, true), (true, true)] {
                let points = linear_feature_polygon(10.0, 10.0 + span, 20.0, 4.0, forward, reverse);
                assert!(points.iter().all(|p| p.x >= 10.0 && p.x <= 10.0 + span));
                assert_eq!(points.contains(&Pos2::new(10.0 + span, 20.0)), forward);
                assert_eq!(points.contains(&Pos2::new(10.0, 20.0)), reverse);
            }
        }
    }
}

/// 外置标签（特征色浅底圆角盒 + 细引线）。文字颜色按底色明暗自选（SnapGene 的浅底色黑字）。
#[allow(clippy::too_many_arguments)]
fn draw_chip(
    p: &egui::Painter,
    text: &str,
    anchor: Pos2,
    align: Align2,
    color: Color32,
    leader_from: Pos2,
    picked: bool,
) -> egui::Rect {
    let font = map_font();
    let fill = theme::mix(color, Color32::WHITE, 0.68);
    let text_color = if theme::luma(fill) < 0.52 {
        Color32::WHITE
    } else {
        theme::TEXT
    };
    let size = p
        .layout_no_wrap(text.to_string(), font.clone(), text_color)
        .size();
    let pad = Vec2::new(5.0, 1.0);
    let box_size = size + pad * 2.0;
    let a0 = align
        .pos_in_rect(&egui::Rect::from_min_size(Pos2::ZERO, box_size))
        .to_vec2();
    let r = egui::Rect::from_min_size(anchor - a0, box_size);
    p.rect_filled(r, 3.0, fill);
    p.rect_stroke(
        r,
        3.0,
        if picked {
            Stroke::new(1.6, theme::ACCENT)
        } else {
            Stroke::new(0.9, theme::darken(color, 0.3))
        },
        egui::StrokeKind::Inside,
    );
    p.line_segment(
        [leader_from, anchor],
        Stroke::new(
            if picked { 1.4 } else { 0.8 },
            if picked {
                theme::ACCENT
            } else {
                theme::with_alpha(theme::TEXT_DIM, 120)
            },
        ),
    );
    p.text(r.min + pad, Align2::LEFT_TOP, text, font, text_color);
    r
}

/// 旋转文字块的屏幕矩形（命中测试用；取旋转后包围盒）。
fn rot_text_rect(
    p: &egui::Painter,
    text: &str,
    font: &FontId,
    anchor: Pos2,
    angle: f32,
    align: Align2,
) -> egui::Rect {
    let size = p.layout_no_wrap(text.to_string(), font.clone(), theme::TEXT).size();
    let a0 = align
        .pos_in_rect(&egui::Rect::from_min_size(Pos2::ZERO, size))
        .to_vec2();
    let center = anchor - a0 + size * 0.5;
    let (s, co) = (angle.sin().abs(), angle.cos().abs());
    let half = Vec2::new(
        size.x * co * 0.5 + size.y * s * 0.5,
        size.x * s * 0.5 + size.y * co * 0.5,
    );
    egui::Rect::from_center_size(center, half * 2.0)
}

/// 文字块的屏幕矩形（命中测试用；与 `p.text` 的对齐方式一致）。
fn text_rect(p: &egui::Painter, text: &str, font: FontId, anchor: Pos2, align: Align2) -> egui::Rect {
    let size = p.layout_no_wrap(text.to_string(), font, theme::TEXT).size();
    let a0 = align
        .pos_in_rect(&egui::Rect::from_min_size(Pos2::ZERO, size))
        .to_vec2();
    egui::Rect::from_min_size(anchor - a0, size)
}

/// 水平多段标签：`runs` 依次排列（各段可单独加粗），整块按 `align` 对齐到 `anchor`。
///
/// SnapGene 的酶位点标签即「名称加粗 + 坐标常规」两段。
fn protein_end_labels(len:usize)->Vec<(usize,String,Align2)>{
    if len==0{return Vec::new();}
    if len==1{return vec![(0,"1 aa".into(),Align2::LEFT_BOTTOM)];}
    vec![(0,"1".into(),Align2::LEFT_BOTTOM),(len,format!("{} aa",nfmt(len)),Align2::RIGHT_BOTTOM)]
}

#[cfg(test)]mod protein_backbone_tests{
    #[test]fn endpoints_use_full_sequence_length_not_feature_extent(){
        let labels=super::protein_end_labels(543);
        assert_eq!(labels[0].0,0);assert_eq!(labels[0].1,"1");
        assert_eq!(labels[1].0,543);assert_eq!(labels[1].1,"543 aa");
        assert!(super::protein_end_labels(0).is_empty());
        assert_eq!(super::protein_end_labels(1)[0].1,"1 aa");
    }
}

fn draw_runs(
    p: &egui::Painter,
    runs: &[(String, bool)],
    anchor: Pos2,
    align: Align2,
    font: &FontId,
    color: Color32,
) {
    let prefs=p.ctx().data(|d|d.get_temp::<crate::preferences::Preferences>(egui::Id::unique("enzyme_render_preferences"))).unwrap_or_default();
    let render_color=|text:&str| if prefs.green_blunt && text.split(" - ").all(|n|hx_enzymes::by_name(n).is_some_and(|e|e.overhang==hx_enzymes::Overhang::Blunt)){Color32::from_rgb(0,150,45)}else{color};
    let widths: Vec<f32> = runs
        .iter()
        .map(|(t, b)| {
            p.layout_no_wrap(t.clone(), font.clone(), color).size().x + if *b { 0.5 } else { 0.0 }
        })
        .collect();
    let total: f32 = widths.iter().sum();
    let mut x = match align.x() {
        egui::Align::Min => anchor.x,
        egui::Align::Max => anchor.x - total,
        egui::Align::Center => anchor.x - total * 0.5,
    };
    let one = egui::Align2::new(egui::Align::Min, align.y());
    for ((t, bold), w) in runs.iter().zip(widths) {
        if *bold && prefs.bold_unique {
            for d in [Vec2::new(0.5, 0.0), Vec2::new(0.0, 0.4)] {
                p.text(Pos2::new(x, anchor.y) + d, one, t, font.clone(), render_color(t));
            }
        }
        let mut format=egui::TextFormat::simple(font.clone(),render_color(t));format.italics=prefs.italic_enzymes && hx_enzymes::by_name(t.split(" - ").next().unwrap_or("")).is_some();
        let mut job=egui::text::LayoutJob::default();job.append(t,0.0,format);
        let galley=p.layout_job(job);let origin=one.anchor_size(Pos2::new(x,anchor.y),galley.size()).min;p.galley(origin,galley,render_color(t));
        x += w;
    }
}
