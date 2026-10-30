//! 图谱（环状 / 线性）几何：碱基 ↔ 屏幕坐标、特征扇环/矩形带、刻度、标签避让。

use crate::{norm_angle, split_wrap};

/// 分子拓扑。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    Circular,
    Linear,
}

/// 图谱**显示布局**（与分子拓扑分开）。
///
/// SnapGene 的 Map 默认对「环状分子 或 长度 > 1 kb 的分子」画成环：实测 717 bp 的
/// 环状分子仍画直线，3442 bp 环状与 5835 bp **线性**都画环（线性的环在正上方留一个
/// 缺口表示两端）。窗口标题/状态栏仍显示真实拓扑。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MapLayout {
    /// 自动：环状分子或长度 > [`MapLayout::AUTO_RING_MIN_BP`] 时画环。
    #[default]
    Auto,
    /// 强制画环。
    Circular,
    /// 强制画直线。
    Linear,
}

impl MapLayout {
    /// `Auto` 下画成环的最小长度（bp）。
    pub const AUTO_RING_MIN_BP: usize = 1000;

    /// 解析成「是否按环状布局绘制」。
    pub fn is_ring(self, circular: bool, len: usize) -> bool {
        match self {
            MapLayout::Circular => true,
            MapLayout::Linear => false,
            MapLayout::Auto => circular || len > Self::AUTO_RING_MIN_BP,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            MapLayout::Auto => "Auto (ring when > 1 kb)",
            MapLayout::Circular => "Ring (circular layout)",
            MapLayout::Linear => "Line (linear layout)",
        }
    }
}

/// 圆弧采样折线（角度自正上方起、顺时针；y 向下）。
fn arc_points(center: [f32; 2], radius: f32, a0: f32, a1: f32) -> Vec<[f32; 2]> {
    let n = (((a1 - a0).abs() * radius / 5.0) as usize).clamp(6, 480);
    (0..=n)
        .map(|k| {
            let a = a0 + (a1 - a0) * k as f32 / n as f32;
            [
                center[0] + a.sin() * radius,
                center[1] - a.cos() * radius,
            ]
        })
        .collect()
}

/// 图谱视图变换（屏幕坐标 y 向下）。
///
/// - 环状：以 `center` 为圆心、`radius` 为基线圆周；`track` 沿法线外推
///   （`track * ring_width` 像素，正 = 环外，负 = 环内）。
/// - 线性：以 `axis_y` 为轴，`origin_x` 为 bp 0 的横坐标，`px_per_bp` 为缩放；
///   `track` 沿 **-y 方向**（上）外推。
#[derive(Debug, Clone, Copy)]
pub struct MapView {
    pub topo: Topology,
    pub center: [f32; 2],
    pub radius: f32,
    pub origin_x: f32,
    pub axis_y: f32,
    pub px_per_bp: f32,
    /// 环带厚度（特征、位点刻度的径向/纵向尺寸）。
    pub ring_width: f32,
    /// 环状初始旋转（弧度；0 = bp0 在正上方，随后顺时针）。
    pub rotation: f32,
}

impl MapView {
    pub fn circular(center: [f32; 2], radius: f32) -> Self {
        Self {
            topo: Topology::Circular,
            center,
            radius,
            origin_x: 0.0,
            axis_y: 0.0,
            px_per_bp: 1.0,
            ring_width: 18.0,
            rotation: 0.0,
        }
    }

    pub fn linear(origin_x: f32, axis_y: f32, px_per_bp: f32) -> Self {
        Self {
            topo: Topology::Linear,
            center: [0.0, 0.0],
            radius: 0.0,
            origin_x,
            axis_y,
            px_per_bp,
            ring_width: 12.0,
            rotation: 0.0,
        }
    }

    /// 环状：bp → 角度（弧度，自正上方起顺时针）。
    pub fn bp_angle(&self, len: usize, bp: f64) -> f32 {
        if len == 0 {
            return self.rotation;
        }
        let frac = (bp / len as f64) as f32;
        norm_angle(frac * std::f32::consts::TAU + self.rotation)
    }

    /// 线性分子按环状布局时，bp 0（正上方）缺口的角度半宽（弧度）。
    ///
    /// oracle 实测（`/tmp/opencode/sg/sg-rna-ring.png`，ECO 16S rRNA 1541 bases 切环）：
    /// SnapGene 的缺口 = 261.2°..278.5° ⇒ **17.3°**（两端各 8.65°）。缺口是**固定角度**
    /// （与环半径无关），所以我们按角度而不是像素留缺口 —— 之前固定 30 px 在 r≈260 的环上
    /// 只有 6.6°，用户看不出两端（2026-10-18 反馈的图一）。
    pub const OPEN_RING_GAP: f32 = 0.302;

    /// 开环的「主弧」角度区间 `[lo, hi]`（弧度，可能超过 `TAU`；缺口以 `rotation` 为中心）。
    ///
    /// 线性分子按环状显示时，**特征带、骨架、刻度**都要被限制在这个区间里
    /// （oracle：SnapGene 的绿带 rRNA 在缺口处断开，实测缺口同为 17.3°）。
    pub fn open_ring_main_arc(&self) -> (f32, f32) {
        let g = Self::OPEN_RING_GAP * 0.5;
        (self.rotation + g, self.rotation + std::f32::consts::TAU - g)
    }

    /// 环状布局的「**开环**」几何：线性分子按环状画时用。
    ///
    /// 返回 `(主弧折线, [两端装饰弧; 2])`：主弧在 bp 0 的角度上留
    /// [`Self::OPEN_RING_GAP`] 的缺口（随 `rotation` 一起转），两端各有一段略靠内的
    /// 短平行弧（SnapGene 的端点样式）。
    pub fn open_ring_arcs(&self) -> (Vec<[f32; 2]>, [Vec<[f32; 2]>; 2]) {
        let r = self.radius.max(1.0);
        let (a0, a1) = self.open_ring_main_arc();
        let d = (20.0 / r).min(Self::OPEN_RING_GAP * 0.85);
        (
            arc_points(self.center, r, a0, a1),
            [
                arc_points(self.center, r - 4.5, a0, a0 + d),
                arc_points(self.center, r - 4.5, a1 - d, a1),
            ],
        )
    }

    /// 线性：bp → x。
    pub fn bp_x(&self, bp: f64) -> f32 {
        self.origin_x + (bp as f32) * self.px_per_bp
    }

    /// 线性：x → bp（可为小数，便于命中测试）。
    pub fn x_to_bp(&self, x: f32) -> f64 {
        if self.px_per_bp <= 0.0 {
            return 0.0;
        }
        ((x - self.origin_x) / self.px_per_bp) as f64
    }

    /// 碱基处的点：`track` 为轨道偏移（见结构体说明）。
    pub fn bp_point(&self, len: usize, bp: f64, track: f32) -> [f32; 2] {
        match self.topo {
            Topology::Circular => {
                let a = self.bp_angle(len, bp);
                let r = self.radius + track * self.ring_width;
                [
                    self.center[0] + a.sin() * r,
                    self.center[1] - a.cos() * r,
                ]
            }
            Topology::Linear => [self.bp_x(bp), self.axis_y - track * self.ring_width],
        }
    }

    /// 点 → 该点在圆周/轴上的最近碱基（用于命中测试）；返回 `(bp, 轨道偏移)`。
    pub fn hit(&self, len: usize, p: [f32; 2]) -> (f64, f32) {
        match self.topo {
            Topology::Circular => {
                let dx = p[0] - self.center[0];
                let dy = p[1] - self.center[1];
                let r = (dx * dx + dy * dy).sqrt();
                let track = if self.ring_width == 0.0 {
                    0.0
                } else {
                    (r - self.radius) / self.ring_width
                };
                let a = norm_angle(dy.atan2(dx) + std::f32::consts::FRAC_PI_2 - self.rotation);
                let frac = (a / std::f32::consts::TAU) as f64;
                (frac * len as f64, track)
            }
            Topology::Linear => {
                let original = self.x_to_bp(p[0]);
                // 线性命中允许超出 0..len 的“外部”，由调用方裁剪
                let track = if self.ring_width == 0.0 {
                    0.0
                } else {
                    (self.axis_y - p[1]) / self.ring_width
                };
                (original, track)
            }
        }
    }

    /// 环状：一段弧的采样点（含首末），`segments` 为细分段数。
    pub fn arc_points(&self, len: usize, start_bp: f64, end_bp: f64, track: f32, segments: usize) -> Vec<[f32; 2]> {
        let n = segments.max(1);
        (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64;
                let bp = start_bp + (end_bp - start_bp) * t;
                self.bp_point(len, bp, track)
            })
            .collect()
    }

    /// 特征带形状（环状 = 扇环轮廓；线性 = 矩形轮廓）。
    ///
    /// `track_in/track_out` 为内外轨道；`strand_reverse` 只影响线性箭头的
    /// 指向（由调用方用 [`Self::arrow_head`] 绘制）。
    pub fn feature_band(
        &self,
        len: usize,
        start: usize,
        end: usize,
        track_in: f32,
        track_out: f32,
    ) -> Vec<[f32; 2]> {
        let mut pts = Vec::new();
        for (s, e) in split_wrap(start, end, len) {
            if s >= e {
                continue;
            }
            match self.topo {
                Topology::Circular => {
                    // 同上：跨原点时必须把终点角度回绕到起点之后
                    let a0 = self.bp_angle(len, s as f64);
                    let mut a1 = self.bp_angle(len, e as f64);
                    while a1 <= a0 {
                        a1 += std::f32::consts::TAU;
                    }
                    let arc_px = (a1 - a0) * self.radius.max(1.0);
                    let seg = ((arc_px / 2.0) as usize).clamp(2, 512);
                    // 外弧正向 + 内弧反向
                    let outer: Vec<[f32; 2]> = (0..=seg)
                        .map(|i| {
                            let a = a0 + (a1 - a0) * (i as f32 / seg as f32);
                            let r = self.radius + track_out * self.ring_width;
                            [
                                self.center[0] + a.sin() * r,
                                self.center[1] - a.cos() * r,
                            ]
                        })
                        .collect();
                    let inner: Vec<[f32; 2]> = (0..=seg)
                        .rev()
                        .map(|i| {
                            let a = a0 + (a1 - a0) * (i as f32 / seg as f32);
                            let r = self.radius + track_in * self.ring_width;
                            [
                                self.center[0] + a.sin() * r,
                                self.center[1] - a.cos() * r,
                            ]
                        })
                        .collect();
                    pts.extend(outer);
                    pts.extend(inner);
                }
                Topology::Linear => {
                    let x0 = self.bp_x(s as f64);
                    let x1 = self.bp_x(e as f64);
                    let y_in = self.axis_y - track_in * self.ring_width;
                    let y_out = self.axis_y - track_out * self.ring_width;
                    pts.extend([
                        [x0, y_in],
                        [x1, y_in],
                        [x1, y_out],
                        [x0, y_out],
                    ]);
                }
            }
        }
        pts
    }

    /// 特征带的**凸四边形**细分（可直接交给 egui 填充；`feature_band` 的凹轮廓
    /// 不适合直接填充，这里按小段拆成四边形/三角形）。
    pub fn feature_quads(
        &self,
        len: usize,
        start: usize,
        end: usize,
        track_in: f32,
        track_out: f32,
    ) -> Vec<[[f32; 2]; 4]> {
        let mut out = Vec::new();
        for (s, e) in split_wrap(start, end, len) {
            if s >= e {
                continue;
            }
            match self.topo {
                Topology::Circular => {
                    // `bp_angle` 会把角度归一化到 [0, TAU)，跨原点（或视图旋转）时
                    // a1 < a0 —— 必须**回绕**成 a1 > a0，否则弧会朝反方向绕整圈。
                    let a0 = self.bp_angle(len, s as f64);
                    let mut a1 = self.bp_angle(len, e as f64);
                    while a1 <= a0 {
                        a1 += std::f32::consts::TAU;
                    }
                    let arc_px = (a1 - a0) * self.radius.max(1.0);
                    let seg = ((arc_px / 2.0) as usize).clamp(2, 512);
                    let r_in = self.radius + track_in * self.ring_width;
                    let r_out = self.radius + track_out * self.ring_width;
                    let pt = |a: f32, r: f32| {
                        [
                            self.center[0] + a.sin() * r,
                            self.center[1] - a.cos() * r,
                        ]
                    };
                    for i in 0..seg {
                        let aa = a0 + (a1 - a0) * (i as f32 / seg as f32);
                        let ab = a0 + (a1 - a0) * ((i + 1) as f32 / seg as f32);
                        out.push([pt(aa, r_out), pt(ab, r_out), pt(ab, r_in), pt(aa, r_in)]);
                    }
                }
                Topology::Linear => {
                    let x0 = self.bp_x(s as f64);
                    let x1 = self.bp_x(e as f64);
                    let y_in = self.axis_y - track_in * self.ring_width;
                    let y_out = self.axis_y - track_out * self.ring_width;
                    out.push([[x0, y_out], [x1, y_out], [x1, y_in], [x0, y_in]]);
                }
            }
        }
        out
    }

    /// 方向箭头（三角形）顶点：`bp` 处、指向 `forward`（5′→3′）。
    pub fn arrow_head(
        &self,
        len: usize,
        bp: f64,
        forward: bool,
        track: f32,
        size: f32,
    ) -> Vec<[f32; 2]> {
        match self.topo {
            Topology::Circular => {
                let a = self.bp_angle(len, bp);
                let r = self.radius + track * self.ring_width;
                let c = [
                    self.center[0] + a.sin() * r,
                    self.center[1] - a.cos() * r,
                ];
                // 切向（顺时针为正）
                let (tx, ty) = (a.cos(), a.sin());
                let (sx, sy) = (-ty, tx);
                let s = if forward { 1.0 } else { -1.0 };
                let half = size * 0.5;
                vec![
                    [c[0] + tx * s * size * 0.5, c[1] + ty * s * size * 0.5],
                    [c[0] - tx * s * size * 0.5 + sx * half, c[1] - ty * s * size * 0.5 + sy * half],
                    [c[0] - tx * s * size * 0.5 - sx * half, c[1] - ty * s * size * 0.5 - sy * half],
                ]
            }
            Topology::Linear => {
                let x = self.bp_x(bp);
                let y = self.axis_y - track * self.ring_width;
                let half = size * 0.5;
                let s = if forward { 1.0 } else { -1.0 };
                vec![
                    [x + s * half, y],
                    [x - s * half, y - half],
                    [x - s * half, y + half],
                ]
            }
        }
    }

    /// 中心 → 圆周的放射线端点（原点标记 / 引线用）。
    pub fn radial_line(&self, len: usize, bp: f64, track_in: f32, track_out: f32) -> [[f32; 2]; 2] {
        [
            self.bp_point(len, bp, track_in),
            self.bp_point(len, bp, track_out),
        ]
    }
}

/// 刻度间距：从 1/2/5×10^k 中选最小的、使 `px_per_bp * step ≥ min_px` 者。
pub fn nice_step(px_per_bp: f32, min_px: f32) -> f64 {
    if px_per_bp <= 0.0 {
        return 1.0;
    }
    let mut step = 1.0f64;
    let mut mult = 1.0f64;
    loop {
        let px = (step * mult) as f32 * px_per_bp;
        if px >= min_px {
            return step * mult;
        }
        mult = if (mult - 1.0).abs() < f64::EPSILON {
            2.0
        } else if (mult - 2.0).abs() < f64::EPSILON {
            5.0
        } else {
            step *= 10.0;
            1.0
        };
        if step > 1e12 {
            return step;
        }
    }
}

/// 直观刻度：返回 `(bp, 是否主刻度)`。
pub fn ruler_ticks(len: usize, px_per_bp: f32, min_px: f32) -> Vec<(f64, bool)> {
    let step = nice_step(px_per_bp, min_px);
    let mut out = Vec::new();
    let mut bp = 0.0f64;
    while bp <= len as f64 + 1e-9 {
        let major = (bp / (step * 5.0)).round() * step * 5.0;
        out.push((bp, (bp - major).abs() < 1e-6));
        bp += step;
    }
    out
}

/// 刻度（标注步长直接取「整齐值」1/2/5×10ⁿ，每步再细分 `subdiv` 个小刻度）。
///
/// 线性图谱用这个：SnapGene 的轴下数字落在 200/400/600 这类整齐值上，而不是
/// [`ruler_ticks`] 的「5 倍小刻度」（250/500）。
pub fn ruler_ticks_subdiv(len: usize, px_per_bp: f32, min_px: f32, subdiv: u32) -> Vec<(f64, bool)> {
    let major = nice_step(px_per_bp, min_px).max(1.0);
    let minor = major / subdiv.max(1) as f64;
    let mut out = Vec::new();
    let mut bp = 0.0f64;
    while bp <= len as f64 + 1e-9 {
        let k = (bp / major).round();
        out.push((bp, (bp - k * major).abs() < minor * 0.25));
        bp += minor;
    }
    out
}

/// 标签避让（严格版）：放不进 `max_levels` 层的条目返回 `None`。
///
/// SnapGene 在空间不足时会**隐藏**标签（并在图角提示 “N enzymes are hidden”），
/// 而不是把它们塞进已有层造成重叠。
pub fn assign_label_levels_strict(
    pos: &[f32],
    width: &[f32],
    max_levels: usize,
) -> Vec<Option<usize>> {
    assert_eq!(pos.len(), width.len());
    let mut level_right: Vec<f32> = Vec::new();
    let mut order: Vec<usize> = (0..pos.len()).collect();
    order.sort_by(|a, b| pos[*a].partial_cmp(&pos[*b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut lv: Vec<Option<usize>> = vec![None; pos.len()];
    for i in order {
        let left = pos[i] - width[i] / 2.0;
        let right = pos[i] + width[i] / 2.0;
        let mut placed = false;
        for (l, r) in level_right.iter_mut().enumerate() {
            if *r <= left {
                *r = right;
                lv[i] = Some(l);
                placed = true;
                break;
            }
        }
        if !placed && level_right.len() < max_levels {
            level_right.push(right);
            lv[i] = Some(level_right.len() - 1);
        }
    }
    lv
}

/// 标签避让：给定时序位置与宽度，贪心分配到若干“层”。
///
/// 返回每个条目的层号（0 起）。`pos` 需按升序给出（同层内从左到右不重叠）。
pub fn assign_label_levels(pos: &[f32], width: &[f32], max_levels: usize) -> Vec<usize> {
    assert_eq!(pos.len(), width.len());
    let mut level_right: Vec<f32> = Vec::new();
    let mut order: Vec<usize> = (0..pos.len()).collect();
    order.sort_by(|a, b| pos[*a].partial_cmp(&pos[*b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut lv = vec![0usize; pos.len()];
    for i in order {
        let left = pos[i] - width[i] / 2.0;
        let right = pos[i] + width[i] / 2.0;
        let mut placed = false;
        for (l, r) in level_right.iter_mut().enumerate() {
            if *r <= left {
                *r = right;
                lv[i] = l;
                placed = true;
                break;
            }
        }
        if !placed {
            if level_right.len() < max_levels {
                level_right.push(right);
                lv[i] = level_right.len() - 1;
            } else {
                // 层数用尽：放到“最先空出”的一层
                let l = level_right
                    .iter()
                    .enumerate()
                    .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|(l, _)| l)
                    .unwrap();
                level_right[l] = right;
                lv[i] = l;
            }
        }
    }
    lv
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn circular_bp0_is_top() {
        let v = MapView::circular([100.0, 100.0], 80.0);
        let p = v.bp_point(100, 0.0, 0.0);
        assert!((p[0] - 100.0).abs() < 1e-3, "{p:?}");
        assert!((p[1] - 20.0).abs() < 1e-3, "{p:?}"); // 上方
        // 顺时针：1/4 圈 → 右侧
        let q = v.bp_point(100, 25.0, 0.0);
        assert!((q[0] - 180.0).abs() < 1e-3, "{q:?}");
        assert!((q[1] - 100.0).abs() < 1e-3, "{q:?}");
    }

    #[test]
    fn circular_hit_round_trip() {
        let v = MapView::circular([0.0, 0.0], 100.0);
        for bp in [0.0, 10.0, 33.0, 66.5, 99.0] {
            let p = v.bp_point(100, bp, 0.0);
            let (got, track) = v.hit(100, p);
            assert!((got - bp).abs() < 0.05, "bp {bp} → {got}");
            assert!(track.abs() < 1e-3);
        }
    }

    #[test]
    fn wrap_band_has_two_lobes() {
        let v = MapView::circular([0.0, 0.0], 100.0);
        // 90→10 跨原点：两段（各 10 bp + 10 bp）
        let pts = v.feature_band(100, 90, 10, 0.2, 1.0);
        assert!(!pts.is_empty());
        // 全部点都在两半径之间
        for p in &pts {
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(r > 100.0 && r < 130.0, "{p:?} r={r}");
        }
        // 两段的采样点数量应远超单段
        let single = v.feature_band(100, 20, 30, 0.2, 1.0);
        assert!(pts.len() > single.len(), "{} vs {}", pts.len(), single.len());
    }

    #[test]
    fn linear_mapping_and_hit() {
        let v = MapView::linear(50.0, 300.0, 0.5);
        assert!((v.bp_x(100.0) - 100.0).abs() < 1e-4);
        assert!((v.x_to_bp(100.0) - 100.0).abs() < 1e-4);
        // track 向上
        let p = v.bp_point(1000, 200.0, 1.0);
        assert!((p[1] - (300.0 - 12.0)).abs() < 1e-4);
    }

    #[test]
    fn nice_step_picks_readable() {
        // 1 bp = 0.1 px，要求 ≥ 60 px ⇒ step ≥ 600 ⇒ 1000
        assert_eq!(nice_step(0.1, 60.0), 1000.0);
        // 1 bp = 1 px，要求 ≥ 60 px ⇒ 100
        assert_eq!(nice_step(1.0, 60.0), 100.0);
        assert_eq!(nice_step(1.0, 5.0), 5.0);
        assert_eq!(nice_step(0.4, 50.0), 200.0);
    }

    #[test]
    fn strict_levels_drop_overflow() {
        // 三个同宽标签重叠、只给 1 层 ⇒ 只放第一个，其余为 None
        let pos = [10.0f32, 12.0, 14.0];
        let w = [30.0f32, 30.0, 30.0];
        let lv = assign_label_levels_strict(&pos, &w, 1);
        assert_eq!(lv.iter().filter(|l| l.is_some()).count(), 1);
        let lv2 = assign_label_levels_strict(&pos, &w, 3);
        assert!(lv2.iter().all(|l| l.is_some()));
    }

    #[test]
    fn ruler_subdiv_labels_are_round() {
        // 1.28 px/bp、要求 ≥150 px ⇒ 标注步长 200（SnapGene 的 200/400/600）
        let t = ruler_ticks_subdiv(717, 1.28, 150.0, 10);
        let labels: Vec<f64> = t.iter().filter(|x| x.1).map(|x| x.0).collect();
        assert_eq!(labels, vec![0.0, 200.0, 400.0, 600.0]);
        assert!(t.len() >= 30, "{}", t.len());
    }

    #[test]
    fn ruler_ticks_cover_range() {
        let t = ruler_ticks(1000, 1.0, 60.0);
        assert_eq!(t.first().unwrap().0, 0.0);
        assert!(t.last().unwrap().0 <= 1000.0);
        assert!(t.len() >= 10);
        // 每 5 步一个主刻度
        assert!(t[0].1);
        assert!(t.iter().filter(|x| x.1).count() * 5 >= t.len() - 5);
    }

    #[test]
    fn label_levels_avoid_overlap() {
        let pos = [10.0, 20.0, 30.0, 100.0];
        let w = [30.0, 30.0, 30.0, 30.0];
        let lv = assign_label_levels(&pos, &w, 4);
        // 前三个相互重叠（区间 [-5,25] [5,35] [15,45]）⇒ 需三层
        assert_eq!(lv[0], 0);
        assert_eq!(lv[1], 1);
        assert_eq!(lv[2], 2);
        assert_eq!(lv[3], 0, "远处的应复用第 0 层");
    }

    #[test]
    fn label_levels_respect_order_independence() {
        // 输入未排序也应得到同样的层次
        let pos = [100.0, 10.0, 30.0];
        let w = [30.0, 30.0, 30.0];
        let lv = assign_label_levels(&pos, &w, 3);
        assert_eq!(lv[1], 0);
        assert_eq!(lv[2], 1);
        assert_eq!(lv[0], 0);
    }

    #[test]
    fn band_arc_wraps_without_full_circle() {
        // 环状、带一点旋转：起点角度落在终点之后的情况必须回绕成小弧（不能绕整圈）
        let mut v = MapView::circular([0.0, 0.0], 100.0);
        v.rotation = std::f32::consts::PI; // 让 0 bp 处角度落到 π
        let quads = v.feature_quads(1000, 990, 10, -0.4, 0.4);
        // 20 bp / 1000 bp 的弧在 r=100 上 ≈ 12.6 px ⇒ 每段 ≤2 px ⇒ 约 6–7 段
        assert!(quads.len() <= 12, "段数过多（绕整圈了）：{}", quads.len());
        let band = v.feature_band(1000, 990, 10, -0.4, 0.4);
        assert!(!band.is_empty());
        // 所有点都应靠近环（半径 ~100），没有落到圆心附近
        for p in &band {
            let d = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(d > 80.0, "轮廓点离环太远：{d}");
        }
    }

    #[test]
    fn arrow_head_direction() {
        let v = MapView::linear(0.0, 0.0, 1.0);
        let f = v.arrow_head(100, 50.0, true, 0.0, 10.0);
        let r = v.arrow_head(100, 50.0, false, 0.0, 10.0);
        // 正向箭头尖在右侧（x 更大）
        let tip_f = f.iter().map(|p| p[0]).fold(f32::MIN, f32::max);
        let tip_r = r.iter().map(|p| p[0]).fold(f32::MIN, f32::max);
        assert!(tip_f > 50.0 && tip_r <= 55.0);
        let min_f = f.iter().map(|p| p[0]).fold(f32::MAX, f32::min);
        assert!(min_f < 50.0);
        let _ = PI;
    }

    #[test]
    fn map_layout_auto_rule() {
        // SnapGene 实测（oracle 实拍）：600 bp 环状画环；717 bp 线性画直线；
        // 3442 bp 环状、5835 bp 线性都画环 ⇒ Auto 规则 = 环状 || 长度 > 1000
        assert!(MapLayout::Auto.is_ring(true, 600));
        assert!(MapLayout::Auto.is_ring(true, 3442));
        assert!(!MapLayout::Auto.is_ring(false, 717));
        assert!(MapLayout::Auto.is_ring(false, 5835));
        assert!(!MapLayout::Auto.is_ring(false, 1000));
        // 覆盖优先于自动规则
        assert!(MapLayout::Circular.is_ring(false, 100));
        assert!(!MapLayout::Linear.is_ring(true, 50_000));
    }

    #[test]
    fn open_ring_leaves_gap_at_bp0() {
        let v = MapView::circular([100.0, 100.0], 80.0);
        let (main, marks) = v.open_ring_arcs();
        assert!(main.len() > 20);
        assert_eq!(marks.len(), 2);
        // 主弧两端都不落在正上方（bp 0 的位置），角度差 ≈ 缺口角（oracle 实测 17.3°）
        let first = main.first().unwrap();
        let last = main.last().unwrap();
        let ang = |p: &[f32; 2]| (p[0] - 100.0).atan2(100.0 - p[1]);
        let norm = |a: f32| {
            if a < 0.0 {
                a + 2.0 * PI
            } else {
                a
            }
        };
        let a0 = ang(first);
        let a1 = ang(last);
        assert!(a0 > 0.15 && a1 < -0.1, "{a0} {a1}");
        let gap_deg = (a0 - a1) * 180.0 / PI;
        assert!((gap_deg - 17.3).abs() < 1.0, "缺口角 {gap_deg}°");
        // 装饰弧半径更小（在环内侧）
        let m0 = marks[0][0];
        let d0 = ((m0[0] - 100.0).powi(2) + (m0[1] - 100.0).powi(2)).sqrt();
        assert!(d0 < 78.0 && d0 > 70.0, "{d0}");
        // 缺口随旋转移动（旋转 π 后缺口在正下方；atan2 返回 (-π, π]，需回绕）
        let mut r = v;
        r.rotation = PI;
        let (m2, _) = r.open_ring_arcs();
        assert!(norm(ang(&m2[0])) > 3.0, "{:?}", m2[0]);
    }
}
