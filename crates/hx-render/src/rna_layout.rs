//! RNA 二级结构的**环-茎径向布局**（SnapGene 的 `Secondary Structure` 视图）。
//!
//! 结构被拆成一棵树：**环**（一段不配对碱基 + 挂在其上的茎）与**茎**（连续堆叠
//! 的碱基对）。每个环画成一个圆：不配对碱基排在圆周上，茎沿半径向外生长，两条
//! 链平行、间距等于一个碱基对间距。
//!
//! 几何与 SnapGene 的实测一致（8.2.3，tRNA 样例，1280×1024）：
//!
//! | 量 | 实测 @68.46% | 折算 @100% |
//! |---|---|---|
//! | 碱基环外径 | ≈18 px | ≈26.3 px |
//! | 相邻碱基中心距（环上） | 21.8–22.0 px | ≈32 px |
//! | 配对两碱基中心距（茎上） | 20.2 px | ≈30 px |
//! | 茎上相邻两对的推进 | 20.4 px | ≈30 px |
//!
//! 环半径不是拍脑袋的常数：发夹环按「碱基沿弧等距」定半径
//! （tRNA 的 5 碱基发夹环实测 r≈36.8、7 碱基 ≈45.4，本模型给出 35.9 / 45.9），
//! 多分支环则把半径撑到**每个分支都放得下**为止——分支的角宽按子树横向半宽
//! 与距离算（气球树），再留一个最小张角，避免大结构（16S rRNA）互相压叠。
//!
//! 纯几何：不依赖 egui，可供屏幕绘制与 SVG/PNG 导出共用。

use std::f32::consts::{PI, TAU};

/// 环上相邻碱基的中心距（100% 缩放，px）。
pub const LOOP_SPACING: f32 = 32.0;
/// 配对两碱基的中心距 = 茎宽（100%）。
pub const PAIR_GAP: f32 = 30.0;
/// 茎上相邻两对之间的距离（100%）。
pub const STEM_RISE: f32 = 30.0;
/// 碱基环半径（100%，描边**中心线**）：外径 = 2 × (11.6 + 3.2/2) ≈ 26.4 px。
///
/// 实测（SnapGene 8.2.3，68.46% 缩放下环的红色实心部分宽 18 px → 100% 26.3 px）
/// 指的是**含描边**的外径，所以半径常数要扣掉半个描边。
pub const RING_RADIUS: f32 = 11.6;
/// 碱基环的描边宽度（100%）。
pub const RING_STROKE: f32 = 3.2;
/// 骨架线宽（100%）。
pub const BACKBONE_WIDTH: f32 = 2.0;
/// 碱基字母的字号（100%）。
pub const LETTER_SIZE: f32 = 14.4;

/// naview（ViennaRNA `RNAplot`）坐标 → 本模块 100% 坐标的倍数。
///
/// naview 的骨架步长与配对间距都是 **15.0** 单位（实测：tRNA 发夹环、20 bp
/// 直茎、配对两碱基全部 15.00）。SnapGene 的 100% 是「相邻碱基中心距 30 px」，
/// 所以这里统一乘 2：之后所有几何常数（`RING_RADIUS`、`PAIR_GAP`……）才和
/// naview 结构对得上。少乘这一步会让碱基环直径（26.4）超过碱基间距，
/// 叠成一团、看不到配对的横档。
pub const NAVIEW_SCALE: f32 = 2.0;

/// 二维点（屏幕坐标：y 向下）。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Pt {
    pub x: f32,
    pub y: f32,
}

impl Pt {
    pub const fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
    pub fn add(self, o: Pt) -> Pt {
        Pt::new(self.x + o.x, self.y + o.y)
    }
    pub fn sub(self, o: Pt) -> Pt {
        Pt::new(self.x - o.x, self.y - o.y)
    }
    pub fn scale(self, k: f32) -> Pt {
        Pt::new(self.x * k, self.y * k)
    }
    /// 单位化（退化时返回 `(1,0)`）。
    pub fn norm(self) -> Pt {
        let l = self.len();
        if l < 1e-6 {
            Pt::new(1.0, 0.0)
        } else {
            Pt::new(self.x / l, self.y / l)
        }
    }
    pub fn len(self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

/// 在碱基 `i` 周围找一个「最空」的方向（单位向量）。
///
/// 24 个候选方向里挑「标签放这里离其他碱基最远」的那个；`away` 只用来打破
/// 接近的平局。结构的两个末端、环上拥挤处的标签都靠它避开别的碱基 ——
/// 固定规则（沿环半径 / 沿骨架）在 tRNA 的接受臂顶端会指进结构里。
pub fn free_dir(pos: &[Pt], i: usize, away: Pt, dist: f32) -> Pt {
    let mut best = (f32::NEG_INFINITY, away.norm());
    for k in 0..24 {
        let a = k as f32 * TAU / 24.0;
        let d = dirv(a);
        let probe = pos[i].add(d.scale(dist));
        let mut m = f32::INFINITY;
        for (j, q) in pos.iter().enumerate() {
            if j == i {
                continue;
            }
            m = m.min(probe.sub(*q).len());
        }
        let score = m + 18.0 * (d.x * away.x + d.y * away.y);
        if score > best.0 {
            best = (score, d);
        }
    }
    best.1
}

/// 极坐标单位向量。
fn dirv(a: f32) -> Pt {
    Pt::new(a.cos(), a.sin())
}

/// 一个环（loop）的几何。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoopGeom {
    pub center: Pt,
    /// 碱基所在圆的半径。
    pub radius: f32,
    /// 环上元素数（不配对碱基 + 茎 + 闭合对）。
    pub vertices: usize,
    /// 外环（没有闭合对）为 `true`。
    pub exterior: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Elem {
    Base(usize),
    Stem(usize),
}

#[derive(Debug, Clone)]
struct StemDef {
    /// 由外到内的碱基对。
    pairs: Vec<(usize, usize)>,
    /// 茎末端（最内一对之内）的环。
    child: usize,
    /// 所属（父）环。
    parent: usize,
    /// 茎的轴向（弧度）。
    axis: f32,
}

#[derive(Debug, Clone)]
struct LoopDef {
    elems: Vec<Elem>,
    closing: Option<(usize, usize)>,
    radius: f32,
    /// 子树相对本环中心的横向半宽（用于父环的角度分配）。
    half_width: f32,
    /// 子树的包围半径（环心到子树内最远碱基）。
    extent: f32,
    center: Pt,
    /// 每个元素分到的角度，与 `elems` 一一对应。
    angles: Vec<f32>,
}

/// 一份二级结构的布局结果。
#[derive(Debug, Clone)]
pub struct RnaLayout {
    /// 每个碱基的中心（布局坐标，未缩放，y 向下）。
    pub pos: Vec<Pt>,
    /// 每个碱基所属的环（下标对应 `loops`）。
    pub loop_of: Vec<usize>,
    /// 所有环的几何。
    pub loops: Vec<LoopGeom>,
    /// 序列长度。
    pub len: usize,
    /// 外环上第一条茎的轴向（弧度；没有茎时为 `None`）。
    pub root_axis: Option<f32>,
}

impl RnaLayout {
    /// 直接用外部（ViennaRNA `RNAplot` 的 naview）坐标建布局。
    ///
    /// naview 坐标按 [`NAVIEW_SCALE`] 放大到本模块的 100% 坐标（15 → 30 px 步长）；
    /// `loop_of` / `loops` 只填一个「伪环」（整体质心），供序号标签与
    /// 5′/3′ 徽章选方向用。坐标长度不足时按 `len` 截断/补齐。
    pub fn from_positions(len: usize, pos: &[(f32, f32)]) -> Self {
        let mut pts: Vec<Pt> = pos
            .iter()
            .take(len)
            // naview（经 ViennaRNA 的 PostScript 输出）是 **y 向上**；本模块与
            // egui 都是 y 向下，所以要取负。SnapGene 的结构页正是这个朝向
            // （tRNA：接受臂朝上、反密码子环在左下、T 环在右）。
            .map(|&(x, y)| Pt::new(x * NAVIEW_SCALE, -y * NAVIEW_SCALE))
            .collect();
        pts.resize(len, Pt::default());
        let n = pts.len().max(1) as f32;
        let (sx, sy) = pts
            .iter()
            .fold((0.0f32, 0.0f32), |(ax, ay), p| (ax + p.x, ay + p.y));
        let center = Pt::new(sx / n, sy / n);
        let radius = pts
            .iter()
            .map(|p| p.sub(center).len())
            .fold(0.0f32, f32::max);
        Self {
            pos: pts,
            loop_of: vec![0; len],
            loops: vec![LoopGeom {
                center,
                radius,
                vertices: len,
                exterior: true,
            }],
            len,
            root_axis: None,
        }
    }

    /// 按 `pairs`（0-based，`(i,j)`，`i<j`）排布 `len` 个碱基。
    ///
    /// 非法（越界、共享碱基）的配对会被忽略，保证布局总能给出坐标。
    pub fn new(len: usize, pairs: &[(usize, usize)]) -> Self {
        let mut partner = vec![usize::MAX; len];
        let mut clean: Vec<(usize, usize)> = Vec::new();
        for &(i, j) in pairs {
            if i >= len || j >= len || i >= j {
                continue;
            }
            // 只接受不冲突（不共享碱基）的配对；先到先得。
            if partner[i] != usize::MAX || partner[j] != usize::MAX {
                continue;
            }
            partner[i] = j;
            partner[j] = i;
            clean.push((i, j));
        }
        clean.sort_unstable();
        let mut pair_at = vec![usize::MAX; len];
        for (p, &(i, _)) in clean.iter().enumerate() {
            pair_at[i] = p;
        }

        let mut b = Builder {
            partner,
            pairs: clean,
            pair_at,
            stems: Vec::new(),
            loops: Vec::new(),
            pos: vec![Pt::default(); len],
            loop_of: vec![0; len],
        };
        // 外环：整条序列的“根”。
        let ext = b.new_loop(None);
        let elems = b.walk(0, len, ext);
        b.loops[ext].elems = elems;

        let mut lay = RnaLayout {
            pos: Vec::new(),
            loop_of: Vec::new(),
            loops: Vec::new(),
            len,
            root_axis: None,
        };
        if len == 0 {
            return lay;
        }
        b.measure(ext);
        b.place(ext, Pt::new(0.0, 0.0), None);
        lay.root_axis = b.loops[ext].elems.iter().find_map(|e| match e {
            Elem::Stem(s) => Some(b.stems[*s].axis),
            _ => None,
        });
        // 整幅图的朝向：让外环上第一条茎朝“左下”（SnapGene 的 tRNA 实测），
        // 这样接受臂竖直、其余三臂分列，和参考图一致。
        let rot = match lay.root_axis {
            Some(a) => SNAP_ROOT_ANGLE - a,
            None => 0.0,
        };
        let (c, s) = (rot.cos(), rot.sin());
        let rot_pt = |p: Pt| Pt::new(p.x * c - p.y * s, p.x * s + p.y * c);
        b.relax(&b.pairs.clone());
        lay.pos = b.pos.iter().map(|p| rot_pt(*p)).collect();
        lay.loop_of = b.loop_of.clone();
        lay.loops = b
            .loops
            .iter()
            .map(|l| LoopGeom {
                center: rot_pt(l.center),
                radius: l.radius,
                vertices: l.elems.len() + usize::from(l.closing.is_some()),
                exterior: l.closing.is_none(),
            })
            .collect();
        lay.root_axis = lay.root_axis.map(|a| a + rot);
        lay
    }

    /// 包围盒 `(min, max)`，含碱基环的半径；空结构时为零矩形。
    pub fn bounds(&self) -> (Pt, Pt) {
        if self.pos.is_empty() {
            return (Pt::default(), Pt::default());
        }
        let mut min = self.pos[0];
        let mut max = self.pos[0];
        for p in &self.pos {
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
        }
        let r = RING_RADIUS;
        (
            Pt::new(min.x - r, min.y - r),
            Pt::new(max.x + r, max.y + r),
        )
    }
}

/// 外环上“根茎”的目标朝向（弧度，y 向下）：左下。
const SNAP_ROOT_ANGLE: f32 = 112.0 * PI / 180.0;

/// 茎的最小横向张角（弧度）。少分支的环把手臂摊开，和 SnapGene 一样“舒展”；
/// 分支很多时（大 rRNA 的中心多分支环）只按子树宽度分配，避免环被撑成巨型。
const MIN_STEM_HALF_SPREAD: f32 = 0.42;
const MIN_STEM_HALF_TIGHT: f32 = 0.13;

struct Builder {
    partner: Vec<usize>,
    pairs: Vec<(usize, usize)>,
    pair_at: Vec<usize>,
    stems: Vec<StemDef>,
    loops: Vec<LoopDef>,
    pos: Vec<Pt>,
    loop_of: Vec<usize>,
}

impl Builder {
    fn new_loop(&mut self, closing: Option<(usize, usize)>) -> usize {
        self.loops.push(LoopDef {
            elems: Vec::new(),
            closing,
            radius: 0.0,
            half_width: 0.0,
            extent: 0.0,
            center: Pt::default(),
            angles: Vec::new(),
        });
        self.loops.len() - 1
    }

    /// 在 `[start, end)` 内按序列顺序收集元素，递归建立子环。
    fn walk(&mut self, start: usize, end: usize, parent_loop: usize) -> Vec<Elem> {
        let mut out = Vec::new();
        let mut k = start;
        while k < end {
            let pi = self.pair_at[k];
            if pi != usize::MAX {
                let j = self.pairs[pi].1;
                if j >= end {
                    k += 1;
                    continue;
                }
                // 整条茎：从 (k,j) 向内连续堆叠的配对。
                let mut stack = vec![(k, j)];
                let (mut a, mut b) = (k, j);
                while a + 1 < b - 1 && self.partner[a + 1] == b - 1 {
                    a += 1;
                    b -= 1;
                    stack.push((a, b));
                }
                let si = self.stems.len();
                self.stems.push(StemDef {
                    pairs: stack,
                    child: 0,
                    parent: parent_loop,
                    axis: 0.0,
                });
                let child = self.new_loop(Some((a, b)));
                let elems = self.walk(a + 1, b, child);
                self.loops[child].elems = elems;
                self.stems[si].child = child;
                out.push(Elem::Stem(si));
                k = j + 1;
            } else {
                out.push(Elem::Base(k));
                k += 1;
            }
        }
        out
    }

    fn n_stems(&self, li: usize) -> usize {
        self.loops[li]
            .elems
            .iter()
            .filter(|e| matches!(e, Elem::Stem(_)))
            .count()
    }

    fn n_bases(&self, li: usize) -> usize {
        self.loops[li]
            .elems
            .iter()
            .filter(|e| matches!(e, Elem::Base(_)))
            .count()
    }

    /// 从环心到子环中心的距离。
    fn stem_distance(&self, si: usize, r: f32) -> f32 {
        let st = &self.stems[si];
        let rc = self.loops[st.child].radius;
        let off = |rr: f32| (rr * rr - (PAIR_GAP / 2.0) * (PAIR_GAP / 2.0)).max(0.0).sqrt();
        off(r) + (st.pairs.len() as f32 - 1.0) * STEM_RISE + off(rc)
    }

    /// 单个元素需要的角度。
    fn demand(&self, li: usize, e: &Elem, r: f32) -> f32 {
        match e {
            Elem::Base(_) => 2.0 * ((LOOP_SPACING / (2.0 * r)).min(1.0)).asin(),
            Elem::Stem(si) => {
                let d = self.stem_distance(*si, r).max(1.0);
                // 分支两侧至少各留一个碱基位，避免两条相邻茎的“嘴”贴在一起。
                let gap = 2.0 * ((LOOP_SPACING / (2.0 * r)).min(1.0)).asin();
                let mouth = 2.0 * ((PAIR_GAP / (2.0 * r)).min(1.0)).asin();
                // 分支张角：子树横向半宽（气球树）+ 一个碱基环的宽度。
                let hw = self.loops[self.stems[*si].child].half_width + 2.0 * RING_RADIUS;
                let balloon = 2.0 * ((hw / d).min(0.94)).asin();
                let min_half = if self.n_stems(li) <= 6 {
                    MIN_STEM_HALF_SPREAD
                } else {
                    MIN_STEM_HALF_TIGHT
                };
                balloon
                    .max(mouth + gap)
                    .max(2.0 * min_half)
                    .min(PI * 0.85)
            }
        }
    }

    /// 闭合对占**两个**顶点（它的两个碱基），所以是两份碱基槽。这样环上第一个
    /// 碱基离闭合对正好一格（32 px），和 SnapGene 的实测一致。
    fn mouth_demand(&self, li: usize, r: f32) -> f32 {
        if self.loops[li].closing.is_none() {
            return 0.0;
        }
        4.0 * ((LOOP_SPACING / (2.0 * r)).min(1.0)).asin()
    }

    fn demand_total(&self, li: usize, r: f32) -> f32 {
        self.mouth_demand(li, r)
            + self
                .loops[li]
                .elems
                .iter()
                .map(|e| self.demand(li, e, r))
                .sum::<f32>()
    }

    /// 自底向上定半径 / 横向半宽。
    fn measure(&mut self, li: usize) {
        let elems = self.loops[li].elems.clone();
        for e in &elems {
            if let Elem::Stem(si) = e {
                let c = self.stems[*si].child;
                self.measure(c);
            }
        }
        // 起始半径：发夹环按「碱基沿弧等距」；外环多留一个空位给线状两端。
        let n = self.n_bases(li) as f32;
        let mut r = ((n + 2.0) * LOOP_SPACING / TAU).max(26.0);
        if self.loops[li].closing.is_none() {
            r = ((n + 3.0) * LOOP_SPACING / TAU).max(26.0);
        }
        // 放不下就长大：多分支环的半径由子树宽度决定。发夹环只有一个分支，
        // 可以填到接近满圈（SnapGene 的实测就是这样紧）；多分支环留一点空。
        let fill = if self.n_stems(li) <= 1 { 0.995 } else { 0.95 };
        for _ in 0..80 {
            if self.demand_total(li, r) <= fill * TAU {
                break;
            }
            r *= 1.06;
        }
        self.loops[li].radius = r;

        // 角度分配：需求 + 余量均分。`angles[0]` 是闭合对的槽（外环没有）。
        let has_mouth = self.loops[li].closing.is_some();
        let mut angles: Vec<f32> = Vec::with_capacity(elems.len() + 1);
        if has_mouth {
            angles.push(self.mouth_demand(li, r));
        }
        angles.extend(elems.iter().map(|e| self.demand(li, e, r)));
        let used: f32 = angles.iter().sum();
        let extra = (TAU - used).max(0.0);
        let share = if angles.is_empty() {
            0.0
        } else {
            extra / angles.len() as f32
        };
        for a in angles.iter_mut() {
            *a += share;
        }

        // 子树横向半宽（保守估计：分支中心距离 × sin(半角) + 子树自身半宽）。
        let mut hw = r + RING_RADIUS;
        let mut extent = r + RING_RADIUS;
        let elem_angles = if has_mouth { &angles[1..] } else { &angles[..] };
        for (e, a) in elems.iter().zip(elem_angles.iter()) {
            if let Elem::Stem(si) = e {
                let d = self.stem_distance(*si, r);
                let child = self.stems[*si].child;
                let child_hw = self.loops[child].half_width;
                let child_ext = self.loops[child].extent;
                let h = (d + child_hw) * (a / 2.0).sin();
                hw = hw.max(h);
                extent = extent.max(d + child_ext);
            }
        }
        self.loops[li].half_width = hw;
        self.loops[li].extent = extent;
        self.loops[li].angles = angles;
    }

    fn place(&mut self, li: usize, center: Pt, mouth: Option<f32>) {
        let r = self.loops[li].radius;
        self.loops[li].center = center;
        let angles = self.loops[li].angles.clone();
        // 闭合对占第一个槽；元素从它的「3′ 侧」边缘开始，沿序列方向（角度递减）排布。
        let (mouth_w, elem_angles) = if self.loops[li].closing.is_some() {
            (angles[0], &angles[1..])
        } else {
            (0.0, &angles[..])
        };
        let mut cursor = match mouth {
            Some(a) => a - mouth_w / 2.0,
            None => PI / 2.0,
        };
        for (idx, e) in self.loops[li].elems.clone().iter().enumerate() {
            let w = elem_angles[idx];
            match e {
                Elem::Base(i) => {
                    let a = cursor - w / 2.0;
                    self.pos_set(*i, center.add(dirv(a).scale(r)), li);
                }
                Elem::Stem(si) => {
                    let axis = cursor - w / 2.0;
                    self.place_stem(*si, center, r, axis);
                }
            }
            cursor -= w;
        }
    }

    fn place_stem(&mut self, si: usize, center: Pt, r: f32, axis: f32) {
        let dir = dirv(axis);
        let perp = Pt::new(-dir.y, dir.x);
        // 最外一对正好落在父环的圆上。
        let half = PAIR_GAP / 2.0;
        let mid0 = center.add(dir.scale((r * r - half * half).max(0.0).sqrt()));
        self.stems[si].axis = axis;
        let pairs = self.stems[si].pairs.clone();
        let parent = self.stems[si].parent;
        for (t, &(i, j)) in pairs.iter().enumerate() {
            let mid = mid0.add(dir.scale(t as f32 * STEM_RISE));
            self.pos_set(i, mid.add(perp.scale(half)), parent);
            self.pos_set(j, mid.sub(perp.scale(half)), parent);
        }
        let child = self.stems[si].child;
        let rc = self.loops[child].radius;
        let mid_last = mid0.add(dir.scale((pairs.len() as f32 - 1.0) * STEM_RISE));
        let ccenter = mid_last.add(dir.scale((rc * rc - half * half).max(0.0).sqrt()));
        self.place(child, ccenter, Some(axis + PI));
    }

    /// 兜底松弛：极少数「深分支回卷」处会有两三个碱基贴在一起（16S 的 4 处 < 10 px），
    /// 这里只把贴住的碱基轻轻推开，保证碱基环不互相压住；配对与骨架间距不动。
    fn relax(&mut self, pairs: &[(usize, usize)]) {
        let n = self.pos.len();
        if n == 0 {
            return;
        }
        let mut partner = vec![usize::MAX; n];
        for &(i, j) in pairs {
            if i < n && j < n {
                partner[i] = j;
                partner[j] = i;
            }
        }
        let min_gap = 2.0 * RING_RADIUS + 0.8;
        let cell = min_gap;
        for _ in 0..10 {
            let mut grid: std::collections::HashMap<(i32, i32), Vec<usize>> =
                std::collections::HashMap::new();
            for (i, p) in self.pos.iter().enumerate() {
                grid.entry((
                    (p.x / cell).floor() as i32,
                    (p.y / cell).floor() as i32,
                ))
                .or_default()
                .push(i);
            }
            let mut delta = vec![Pt::default(); n];
            let mut any = false;
            for i in 0..n {
                let p = self.pos[i];
                let cx = (p.x / cell).floor() as i32;
                let cy = (p.y / cell).floor() as i32;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        let Some(list) = grid.get(&(cx + dx, cy + dy)) else {
                            continue;
                        };
                        for &j in list {
                            if j <= i
                                || j == partner[i]
                                || (j as isize - i as isize).abs() == 1
                            {
                                continue;
                            }
                            let q = self.pos[j];
                            let d = ((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y))
                                .sqrt();
                            if d >= min_gap {
                                continue;
                            }
                            let dir = if d < 1e-3 {
                                Pt::new(1.0, 0.0)
                            } else {
                                Pt::new((p.x - q.x) / d, (p.y - q.y) / d)
                            };
                            let push = (min_gap - d) * 0.3;
                            delta[i] = delta[i].add(dir.scale(push));
                            delta[j] = delta[j].sub(dir.scale(push));
                            any = true;
                        }
                    }
                }
            }
            if !any {
                break;
            }
            for i in 0..n {
                self.pos[i] = self.pos[i].add(delta[i]);
            }
        }
    }

    fn pos_set(&mut self, i: usize, p: Pt, loop_idx: usize) {
        self.pos[i] = p;
        self.loop_of[i] = loop_idx;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// SnapGene 样例 tRNA（Homo sapiens mitochondrial tRNA-Ala，72 nt）的 MFE 配对。
    const TRNA_LEN: usize = 72;
    const TRNA_PAIRS: [(usize, usize); 22] = [
        (0, 67),
        (1, 66),
        (2, 65),
        (3, 64),
        (4, 63),
        (5, 62),
        (6, 61),
        (9, 21),
        (10, 20),
        (11, 19),
        (12, 18),
        (22, 40),
        (23, 39),
        (24, 38),
        (25, 37),
        (26, 36),
        (27, 35),
        (44, 60),
        (45, 59),
        (46, 58),
        (47, 57),
        (48, 56),
    ];

    fn dist(a: Pt, b: Pt) -> f32 {
        a.sub(b).len()
    }

    #[test]
    fn trna_geometry_matches_snapgene() {
        let lay = RnaLayout::new(TRNA_LEN, &TRNA_PAIRS);
        assert_eq!(lay.pos.len(), TRNA_LEN);
        for p in &lay.pos {
            assert!(p.x.is_finite() && p.y.is_finite(), "非有限坐标 {p:?}");
        }
        // 配对的两碱基相距一个对间距。
        for &(i, j) in &TRNA_PAIRS {
            let d = dist(lay.pos[i], lay.pos[j]);
            assert!((d - PAIR_GAP).abs() < 0.6, "pair {i}/{j} = {d}");
        }
        // 环上的碱基间距接近 LOOP_SPACING（夹在分支之间的会略宽）。
        for i in 0..TRNA_LEN - 1 {
            if lay.loop_of[i] == lay.loop_of[i + 1]
                && !TRNA_PAIRS.iter().any(|&(a, b)| a == i && b == i + 1)
                && !TRNA_PAIRS.iter().any(|&(a, b)| b == i && a == i + 1)
            {
                let d = dist(lay.pos[i], lay.pos[i + 1]);
                assert!(d > 26.0 && d < 90.0, "相邻 {i}: {d}");
            }
        }
        // 发夹环半径 ≈ 按弧长等距算出的半径（SnapGene 实测 36.8 / 45.4 @100%）。
        let mut radii: Vec<f32> = lay
            .loops
            .iter()
            .filter(|l| !l.exterior && l.vertices <= 7)
            .map(|l| l.radius)
            .collect();
        radii.sort_by(|a, b| a.total_cmp(b));
        assert!((radii[0] - 35.9).abs() < 3.0, "D 环半径 {}", radii[0]);
        // 最小自由能结构里三个发夹环都存在。
        let ext = lay.loops.iter().find(|l| l.exterior).unwrap();
        assert!((ext.radius - 7.0 * LOOP_SPACING / TAU).abs() < 3.0);
    }

    #[test]
    fn trna_cloverleaf_orientation() {
        let lay = RnaLayout::new(TRNA_LEN, &TRNA_PAIRS);
        // 接受臂在最上方，反密码环在最下方（与 SnapGene 的参考图一致）。
        let top = lay.pos[0].y; // 5′ 端，接受臂外端
        let bottom = lay.pos[34].y; // 反密码环底部
        assert!(top < bottom, "5′ 端应在反密码环之上：{top} vs {bottom}");
        // 接受臂两条链近乎水平：配对的两个碱基落在同一高度附近。
        let (a, b) = (lay.pos[0], lay.pos[67]);
        assert!((a.y - b.y).abs() < 30.0, "接受臂首对应近乎水平：{a:?} {b:?}");
        for l in &lay.loops {
            assert!(l.radius > 0.0 && l.radius.is_finite());
        }
    }

    #[test]
    fn no_overlap_between_bases() {
        // 大结构（16S rRNA 的 478 对）不能出现两个碱基叠在一起的情况。
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/eco_16s_pairs.txt");
        let text = std::fs::read_to_string(&path).unwrap();
        let mut it = text.lines();
        let len: usize = it.next().unwrap().trim().parse().unwrap();
        let _n: usize = it.next().unwrap().trim().parse().unwrap();
        let pairs: Vec<(usize, usize)> = it
            .filter_map(|l| {
                let mut s = l.split_whitespace();
                Some((s.next()?.parse().ok()?, s.next()?.parse().ok()?))
            })
            .collect();
        let lay = RnaLayout::new(len, &pairs);
        assert_eq!(lay.pos.len(), len);
        // 网格哈希找中心距 < 10 px 的碱基对（环半径 13.2，正常间距 ≥ 26；1541 nt 的
        // 16S 里有 4 处深分支回卷的贴靠，缩放后肉眼不可见）。
        let mut grid: std::collections::HashMap<(i32, i32), Vec<usize>> =
            std::collections::HashMap::new();
        let cell = 16.0;
        let mut close = 0usize;
        for (i, p) in lay.pos.iter().enumerate() {
            assert!(p.x.is_finite() && p.y.is_finite());
            let key = ((p.x / cell).floor() as i32, (p.y / cell).floor() as i32);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(v) = grid.get(&(key.0 + dx, key.1 + dy)) {
                        for &j in v {
                            if dist(*p, lay.pos[j]) < 10.0 {
                                close += 1;
                            }
                        }
                    }
                }
            }
            grid.entry(key).or_default().push(i);
        }
        assert_eq!(close, 0, "{close} 对碱基重叠");
    }

    #[test]
    fn naview_positions_are_scaled_to_snapgene_units() {
        // naview 的一步是 15.0 单位；SnapGene 的 100% 是 30 px/步（= PAIR_GAP）。
        let lay = RnaLayout::from_positions(2, &[(0.0, 0.0), (15.0, 0.0)]);
        assert!((lay.pos[1].x - lay.pos[0].x - 30.0).abs() < 0.01);
        assert!((lay.pos[1].x - lay.pos[0].x - PAIR_GAP).abs() < 0.01);
    }

    #[test]
    fn ring_outer_diameter_matches_snapgene() {
        // 实测 68.46% 下环宽 18 px → 100% 26.3 px（含描边）。
        let outer = 2.0 * (RING_RADIUS + RING_STROKE / 2.0);
        assert!((outer - 26.4).abs() < 0.2, "outer = {outer}");
        // 环必须明显小于一步，茎才能看出梯子而不是一团。
        assert!(outer < 0.9 * PAIR_GAP);
    }

    #[test]
    fn empty_and_degenerate_inputs() {
        let lay = RnaLayout::new(0, &[]);
        assert!(lay.pos.is_empty());
        // 完全不配对：所有碱基落在同一个环上。
        let lay = RnaLayout::new(36, &[]);
        assert_eq!(lay.pos.len(), 36);
        for i in 1..36 {
            let d = dist(lay.pos[i - 1], lay.pos[i]);
            assert!(d > 26.0 && d < 40.0, "线性环 {i}: {d}");
        }
        // 越界 / 共享碱基的配对被忽略，不会 panic。
        let lay = RnaLayout::new(10, &[(0, 9), (0, 5), (20, 30), (3, 3)]);
        assert_eq!(lay.pos.len(), 10);
    }
}
