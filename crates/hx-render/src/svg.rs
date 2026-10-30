//! # SVG 导出（无外部依赖）
//!
//! 把 [`crate::map`] / [`crate::seq`] 的几何翻译成 SVG 字符串：图谱（环状/线性）
//! 与序列视图各一个入口，供 CLI 与 UI 导出矢量图（发表/制图用）。
//!
//! - 输入是**几何无关**的 [`Segment`] / [`SiteMark`]（与 `.dna` 格式 crate 解耦）；
//! - 坐标系与屏幕一致（y 向下）；默认白底（[`Palette::default`]）；
//! - 数字保留 2 位小数；文本按 XML 转义；
//! - 环状扇环按 2 px 采样折线（与 UI 同一几何：`feature_quads`）。

use crate::map::{assign_label_levels, ruler_ticks, MapView, Topology};
use crate::{split_wrap, thousands};
use crate::circular_labels::{CircularLabel,CircularLabelLayout};

/// 一段特征片段（0-based 半开；环状 `start > end` 表示跨原点，`start == end` 为整圈）。
#[derive(Debug, Clone, Default)]
pub struct Segment {
    pub start: usize,
    pub end: usize,
    pub name: String,
    /// 特征类型（`CDS` / `promoter` …）；名称为空时用作标签。
    pub kind: String,
    /// `#rrggbb`（空 = 用调色板兜底色）。
    pub color: String,
    /// 1 正向 / -1 反向 / 0 无方向。
    pub strand: i8,
}

/// 酶切位点标记（`pos` 为 0-based 识别序列起点）。
#[derive(Debug, Clone, Default)]
pub struct SiteMark {
    pub pos: usize,
    pub name: String,
    pub site_len: usize,
    /// 该酶在本序列中是否唯一（唯一用强调色）。
    pub unique: bool,
}

/// 调色板。
#[derive(Debug, Clone)]
pub struct Palette {
    pub bg: String,
    pub backbone: String,
    pub backbone_thin: String,
    pub tick_major: String,
    pub tick_minor: String,
    pub text: String,
    pub text_dim: String,
    pub site: String,
    pub site_unique: String,
    pub feature_fallback: String,
    /// 碱基色（A/C/G/T/其它）。
    pub base_a: String,
    pub base_c: String,
    pub base_g: String,
    pub base_t: String,
    pub base_other: String,
}

impl Default for Palette {
    fn default() -> Self {
        Self {
            bg: "#ffffff".into(),
            backbone: "#3c4148".into(),
            backbone_thin: "#9aa0a8".into(),
            tick_major: "#4a5058".into(),
            tick_minor: "#a8adb5".into(),
            text: "#282c33".into(),
            text_dim: "#5a6068".into(),
            site: "#8f1fa0".into(),
            site_unique: "#c4145a".into(),
            feature_fallback: "#c0c0c0".into(),
            base_a: "#1f7a1f".into(),
            base_c: "#1b4fb8".into(),
            base_g: "#c87400".into(),
            base_t: "#c02828".into(),
            base_other: "#707070".into(),
        }
    }
}

impl Palette {
    /// 碱基颜色。
    pub fn base_color(&self, b: u8) -> &str {
        match b.to_ascii_uppercase() {
            b'A' => &self.base_a,
            b'C' => &self.base_c,
            b'G' => &self.base_g,
            b'T' | b'U' => &self.base_t,
            _ => &self.base_other,
        }
    }
}

/// 导出样式。
#[derive(Debug, Clone)]
pub struct Style {
    pub width: f32,
    pub height: f32,
    pub show_labels: bool,
    /// 图谱是否画酶切位点。
    pub show_sites: bool,
    /// 中心/顶部标题（一般传文件名/图谱标签）。`None` = 不画。
    pub title: Option<String>,
    /// 副标题（默认自动生成 `N bp · circular|linear`）；`None` = 用自动值。
    pub subtitle: Option<String>,
    /// 环状布局覆盖：`None` = 按调用方给的拓扑；`Some(true/false)` = 强制环/直线。
    ///
    /// [`map_svg`] 在 `None` 时按 [`crate::MapLayout::Auto`] 填充（与界面一致）。
    pub circular_layout: Option<bool>,
    /// 线性分子按环状布局绘制时，环在正上方留**缺口**（两端标记）。
    pub open_ring: bool,
    /// 双链分子（DNA）：骨架画两根平行线；单链（RNA / 单链 DNA / 蛋白）只画一根。
    pub double_strand: bool,
    pub palette: Palette,
}

impl Style {
    pub fn map(width: f32, height: f32) -> Self {
        Self {
            width,
            height,
            show_labels: true,
            show_sites: true,
            title: None,
            subtitle: None,
            circular_layout: None,
            open_ring: false,
            double_strand: true,
            palette: Palette::default(),
        }
    }

    pub fn sequence(width: f32) -> Self {
        Self {
            width,
            height: 0.0, // 序列视图高度按行数自动计算
            show_labels: true,
            show_sites: true,
            title: None,
            subtitle: None,
            circular_layout: None,
            open_ring: false,
            double_strand: true,
            palette: Palette::default(),
        }
    }
}

/// 自动副标题。
fn auto_subtitle(len: usize, circular: bool) -> String {
    let topo = if circular { "circular" } else { "linear" };
    format!("{} bp · {topo}", thousands(len))
}

/// 特征颜色（空/非法 ⇒ 调色板兜底）。
fn feature_fill(color: &str, fallback: &str) -> String {
    let t = color.trim().trim_start_matches('#');
    if t.len() >= 6 && t[..6].bytes().all(|b| b.is_ascii_hexdigit()) {
        format!("#{}", &t[..6].to_ascii_lowercase())
    } else {
        fallback.to_string()
    }
}

/// 混合两色（`t` = 第二色占比）；用于半透明填充的白底近似。
fn mix(a: &str, b: &str, t: f32) -> String {
    let pa = parse_hex(a);
    let pb = parse_hex(b);
    let t = t.clamp(0.0, 1.0);
    let c = |i: usize| (pa[i] as f32 * (1.0 - t) + pb[i] as f32 * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", c(0), c(1), c(2))
}

fn parse_hex(s: &str) -> [u8; 3] {
    let t = s.trim().trim_start_matches('#');
    if t.len() >= 6 && t[..6].bytes().all(|b| b.is_ascii_hexdigit()) {
        let v = u32::from_str_radix(&t[..6], 16).unwrap_or(0);
        [(v >> 16) as u8, (v >> 8) as u8, v as u8]
    } else {
        [0x80, 0x80, 0x80]
    }
}

/// `#rrggbb` 的感知亮度（0..1）。
fn luma_str(c: &str) -> f32 {
    let [r, g, b] = parse_hex(c);
    (0.299 * r as f32 + 0.587 * g as f32 + 0.114 * b as f32) / 255.0
}

/// 环上文字/数字的旋转角（保持可读：下半圈翻转 180°）。
fn readable_angle(a: f32) -> f32 {
    if a.cos() < 0.0 {
        a + std::f32::consts::PI
    } else {
        a
    }
}
struct Svg {
    w: f32,
    h: f32,
    body: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    Start,
    Middle,
    End,
}

impl Anchor {
    fn attr(self) -> &'static str {
        match self {
            Anchor::Start => "start",
            Anchor::Middle => "middle",
            Anchor::End => "end",
        }
    }
}

fn f(x: f32) -> String {
    let v = (x * 100.0).round() / 100.0;
    if v == v.trunc() && v.abs() < 1e9 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

impl Svg {
    fn new(w: f32, h: f32) -> Self {
        Self {
            w,
            h,
            body: String::with_capacity(4096),
        }
    }

    fn line(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, stroke: &str, width: f32) {
        self.body.push_str(&format!(
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
            f(x1),
            f(y1),
            f(x2),
            f(y2),
            stroke,
            f(width)
        ));
    }

    fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, fill: &str, stroke: Option<(&str, f32)>) {
        match stroke {
            Some((c, sw)) => self.body.push_str(&format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
                f(x), f(y), f(w), f(h), fill, c, f(sw)
            )),
            None => self.body.push_str(&format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\"/>",
                f(x),
                f(y),
                f(w),
                f(h),
                fill
            )),
        }
    }

    fn circle(&mut self, cx: f32, cy: f32, r: f32, fill: &str, stroke: Option<(&str, f32)>) {
        match stroke {
            Some((c, sw)) => self.body.push_str(&format!(
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{}\"/>",
                f(cx), f(cy), f(r), fill, c, f(sw)
            )),
            None => self.body.push_str(&format!(
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\"/>",
                f(cx),
                f(cy),
                f(r),
                fill
            )),
        }
    }

    fn poly(&mut self, pts: &[[f32; 2]], fill: &str, stroke: Option<(&str, f32)>) {
        if pts.len() < 2 {
            return;
        }
        let mut d = String::with_capacity(pts.len() * 14 + 16);
        d.push('M');
        for (i, p) in pts.iter().enumerate() {
            if i > 0 {
                d.push('L');
            }
            d.push_str(&f(p[0]));
            d.push(',');
            d.push_str(&f(p[1]));
        }
        if fill != "none" {
            d.push('Z');
        }
        match stroke {
            Some((c, sw)) => self.body.push_str(&format!(
                "<path d=\"{d}\" fill=\"{fill}\" stroke=\"{c}\" stroke-width=\"{}\" stroke-linejoin=\"round\"/>",
                f(sw)
            )),
            None if fill != "none" => self.body.push_str(&format!("<path d=\"{d}\" fill=\"{fill}\"/>")),
            None => self.body.push_str(&format!(
                "<path d=\"{d}\" fill=\"none\" stroke=\"{fill}\" stroke-width=\"1\"/>"
            )),
        }
    }

    fn text(
        &mut self,
        x: f32,
        y: f32,
        size: f32,
        anchor: Anchor,
        fill: &str,
        content: &str,
        bold: bool,
    ) {
        let weight = if bold { " font-weight=\"bold\"" } else { "" };
        self.body.push_str(&format!(
            "<text x=\"{}\" y=\"{}\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"{}\" text-anchor=\"{}\"{weight} fill=\"{}\">{}</text>",
            f(x),
            f(y),
            f(size),
            anchor.attr(),
            fill,
            esc(content)
        ));
    }

    /// 旋转文字（绕锚点 `(x, y)` 逆时针为正的 SVG 约定：`rot_deg` 用屏幕角度，即顺时针为正）。
    #[allow(clippy::too_many_arguments)]
    fn text_rot(
        &mut self,
        x: f32,
        y: f32,
        size: f32,
        anchor: Anchor,
        fill: &str,
        content: &str,
        rot_deg: f32,
        bold: bool,
    ) {
        let weight = if bold { " font-weight=\"bold\"" } else { "" };
        self.body.push_str(&format!(
            "<text x=\"{}\" y=\"{}\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"{}\" text-anchor=\"{}\"{weight} fill=\"{}\" transform=\"rotate({} {} {})\">{}</text>",
            f(x),
            f(y),
            f(size),
            anchor.attr(),
            fill,
            f(rot_deg),
            f(x),
            f(y),
            esc(content)
        ));
    }

    /// 多个字形共用一个 `<text>`（`xs` 与 `chars` 一一对应；等宽字体）。
    fn text_run(&mut self, xs: &[f32], y: f32, size: f32, fill: &str, chars: &str, monospace: bool) {
        if xs.is_empty() || chars.is_empty() {
            return;
        }
        let family = if monospace { "Courier New,Courier" } else { "Helvetica,Arial" };
        let mut xattr = String::with_capacity(xs.len() * 7);
        for (i, x) in xs.iter().enumerate() {
            if i > 0 {
                xattr.push(' ');
            }
            xattr.push_str(&f(*x));
        }
        self.body.push_str(&format!(
            "<text x=\"{xattr}\" y=\"{}\" font-family=\"{family},monospace\" font-size=\"{}\" fill=\"{}\">{}</text>",
            f(y),
            f(size),
            fill,
            esc(chars)
        ));
    }

    fn finish(self, bg: &str) -> String {
        let mut out = String::with_capacity(self.body.len() + 256);
        out.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">",
            f(self.w),
            f(self.h),
            f(self.w),
            f(self.h)
        ));
        out.push_str(&format!(
            "<rect x=\"0\" y=\"0\" width=\"{}\" height=\"{}\" fill=\"{}\"/>",
            f(self.w),
            f(self.h),
            bg
        ));
        out.push_str(&self.body);
        out.push_str("</svg>");
        out
    }
}

/// 环状图谱 SVG。
pub fn circular_map_svg(
    len: usize,
    segments: &[Segment],
    sites: &[SiteMark],
    style: &Style,
    title: Option<&str>,
) -> String {
    const RING: f32 = 22.0;
    let pal = &style.palette;
    let mut svg = Svg::new(style.width, style.height);
    if len == 0 {
        return svg.finish(&pal.bg);
    }
    // SVG has no font shaper. Use conservative glyph bounds, not the former
    // angle*radius layering; both backends share the same callout solver.
    let width=|s:&str|s.chars().map(|c|if c.is_ascii(){8.0}else{12.0}).sum::<f32>();
    let bounds=[18.0,28.0,style.width-18.0,style.height-28.0];
    let mut labels=Vec::new();let mut label_kinds=Vec::new();
    let mut external=vec![false;segments.len()];
    let mut radius=style.width.min(style.height)*0.36;
    let mut layout=CircularLabelLayout::fit(bounds,&labels,1.0);
    for _ in 0..5 {
        labels.clear();label_kinds.clear();
        if style.show_labels {
            for (i,seg) in segments.iter().enumerate() {
                let name=if seg.name.is_empty(){&seg.kind}else{&seg.name};
                let span=crate::span_len(seg.start,seg.end,len);
                if width(name)+8.0>std::f32::consts::TAU*(radius-RING)*span as f32/len as f32 {external[i]=true;}
                if external[i] && !name.is_empty() {
                    labels.push(CircularLabel{angle:(seg.start as f32+span as f32*0.5)/len as f32*std::f32::consts::TAU,width:width(name)+10.0,height:17.5});
                    label_kinds.push((true,i));
                }
            }
            if style.show_sites {
                for (i,site) in sites.iter().enumerate() {
                    labels.push(CircularLabel{angle:site.pos as f32/len as f32*std::f32::consts::TAU,
                        width:width(&format!("{} ({})",site.name,thousands(site.pos)))+4.0,height:15.5});
                    label_kinds.push((false,i));
                }
            }
        }
        layout=CircularLabelLayout::fit(bounds,&labels,1.0);
        let delta=(radius-layout.radius).abs();radius=layout.radius;
        if delta<0.5 {break;}
    }
    let center=layout.center;
    let mut view = MapView::circular(center, radius);
    view.ring_width = RING;
    let circ = std::f32::consts::TAU * radius;
    let pt = |a: f32, r: f32| -> [f32; 2] {
        [center[0] + a.sin() * r, center[1] - a.cos() * r]
    };

    // 骨架（线性分子画成环时留缺口 = 两端标记；单链只有一根线）
    let cores: &[f32] = if style.double_strand { &[-2.0, 2.0] } else { &[0.0] };
    if style.open_ring {
        for dr in cores {
            let mut rv = MapView::circular(center, radius + dr);
            rv.rotation = view.rotation;
            let (main, marks) = rv.open_ring_arcs();
            svg.poly(&main, "none", Some((&pal.backbone, 3.0)));
            for m in &marks {
                svg.poly(m, "none", Some((&pal.backbone, 3.0)));
            }
        }
        let inner = MapView::circular(center, (radius - RING).max(2.0));
        let (im, _) = inner.open_ring_arcs();
        svg.poly(&im, "none", Some((&pal.backbone_thin, 1.0)));
    } else {
        for dr in cores {
            svg.circle(center[0], center[1], radius + dr, "none", Some((&pal.backbone, 3.0)));
        }
        if style.double_strand {
            svg.circle(
                center[0],
                center[1],
                radius - RING,
                "none",
                Some((&pal.backbone_thin, 1.0)),
            );
        }
    }

    // 刻度 + 坐标（跨环 + 环内旋转数字，与界面一致）
    for (bp, major) in ruler_ticks(len, circ / len as f32, 46.0) {
        if bp >= len as f64 {
            continue;
        }
        let a = view.bp_angle(len, bp);
        let r0 = radius + 6.0;
        let r1 = radius - if major { 14.0 } else { 8.0 };
        let (c, w) = if major {
            (&pal.tick_major, 1.4)
        } else {
            (&pal.tick_minor, 1.0)
        };
        let p0 = pt(a, r0);
        let p1 = pt(a, r1);
        svg.line(p0[0], p0[1], p1[0], p1[1], c, w);
        if major {
            let tl = pt(a, radius - RING - 20.0);
            let label = thousands(bp as usize);
            let deg = readable_angle(a).to_degrees();
            svg.text_rot(
                tl[0],
                tl[1] + 3.5,
                10.5,
                if a.cos() < 0.0 { Anchor::End } else { Anchor::Start },
                &pal.text_dim,
                &label,
                deg,
                false,
            );
        }
    }

    // 特征（环**内侧**带 + 方向箭头 + 带内名字，与界面一致）
    for (seg_i,seg) in segments.iter().enumerate() {
        let fill = feature_fill(&seg.color, &pal.feature_fallback);
        for quad in view.feature_quads(len, seg.start, seg.end, -0.92, -0.08) {
            svg.poly(&quad, &fill, None);
        }
        let span = crate::span_len(seg.start, seg.end, len);
        let arc_px = circ * span as f32 / len as f32;
        if seg.strand != 0 {
            if arc_px > 18.0 {
                let bp = if seg.strand > 0 {
                    seg.end as f64
                } else {
                    seg.start as f64
                };
                let tri = view.arrow_head(len, bp, seg.strand > 0, -0.5, 11.0);
                svg.poly(&tri, &mix(&fill, "#ffffff", 0.2), Some((&mix(&fill, "#000000", 0.4), 1.0)));
            }
        }
        // 带内名字：灰字（暗底转白）、沿弧旋转、放得下才写
        let name = if seg.name.is_empty() {
            &seg.kind
        } else {
            &seg.name
        };
        if style.show_labels && !external[seg_i] && !name.is_empty() && arc_px > width(name) + 8.0 {
            let a = view.bp_angle(len, (seg.start + span / 2) as f64);
            let at = pt(a, radius - RING * 0.5);
            let light = luma_str(&fill) > 0.55;
            svg.text_rot(
                at[0],
                at[1] + 3.5,
                11.5,
                Anchor::Middle,
                if light { &pal.text_dim } else { "#ffffff" },
                name,
                readable_angle(a).to_degrees(),
                false,
            );
        }
    }

    // Every site remains represented; nearby sites are not silently discarded.
    if style.show_sites {
        for s in sites {
            let color = if s.unique {
                pal.site_unique.clone()
            } else {
                pal.site.clone()
            };
            let a = view.bp_angle(len, s.pos as f64);
            let p0 = pt(a, radius + RING * 0.04);
            let p1 = pt(a, radius + RING * 0.96);
            svg.line(p0[0], p0[1], p1[0], p1[1], &color, 1.6);
        }
    }
    for (i,&(feature,index)) in label_kinds.iter().enumerate() {
        let Some(callout)=layout.labels[i] else{continue};
        let (name,color,bold,source_radius)=if feature {
            let seg=&segments[index];
            (if seg.name.is_empty(){seg.kind.clone()}else{seg.name.clone()},pal.text.clone(),false,radius-RING*0.5)
        }else{
            let s=&sites[index];
            (if callout.right {format!("{} ({})",s.name,thousands(s.pos))}else{format!("({}) {}",thousands(s.pos),s.name)},
                if s.unique {pal.site_unique.clone()}else{pal.site.clone()},s.unique,radius+4.0)
        };
        svg.poly(&layout.leader(&labels[i],&callout,source_radius),"none",Some((&pal.text_dim,0.8)));
        if feature {
            let r=callout.rect;
            let fill=feature_fill(&segments[index].color,&pal.feature_fallback);
            svg.rect(r[0],r[1],r[2]-r[0],r[3]-r[1],&mix(&fill,"#ffffff",0.68),Some((&fill,0.7)));
            svg.text((r[0]+r[2])*0.5,callout.port[1]+4.0,11.5,Anchor::Middle,&color,&name,false);
        }else{
            svg.text(callout.port[0],callout.port[1]+4.0,11.5,if callout.right {Anchor::Start}else{Anchor::End},&color,&name,bold);
        }
    }
    let hidden=layout.labels.iter().filter(|l|l.is_none()).count();
    if hidden>0 {
        svg.text(12.0,style.height-10.0,11.0,Anchor::Start,&pal.text_dim,&format!("{hidden} label(s) hidden (crowded)"),false);
    }

    // 中心标注
    let title = title.unwrap_or("");
    if !title.is_empty() {
        svg.text(center[0], center[1] + 3.0, 17.0, Anchor::Middle, &pal.text, title, true);
    }
    let sub = style
        .subtitle
        .clone()
        .unwrap_or_else(|| auto_subtitle(len, !style.open_ring));
    svg.text(
        center[0],
        center[1] + 22.0,
        11.5,
        Anchor::Middle,
        &pal.text_dim,
        &sub,
        false,
    );
    svg.finish(&pal.bg)
}

/// 线性图谱 SVG（整条序列适配宽度；两端标 `1` 与 `len`）。
pub fn linear_map_svg(
    len: usize,
    segments: &[Segment],
    sites: &[SiteMark],
    style: &Style,
    title: Option<&str>,
) -> String {
    const LINE_RING: f32 = 26.0;
    let pal = &style.palette;
    let mut svg = Svg::new(style.width, style.height);
    if len == 0 {
        return svg.finish(&pal.bg);
    }
    let top = if title.is_some() { 40.0 } else { 16.0 };
    let axis_y = top + (style.height - top) * 0.40;
    let left = 34.0;
    let right = style.width - 20.0;
    let zoom = ((right - left) / len as f32).max(0.0005);
    let x_of = |bp: f64| -> f32 { left + bp as f32 * zoom };

    if let Some(t) = title {
        svg.text(16.0, 22.0, 14.0, Anchor::Start, &pal.text, t, true);
    }

    // 骨架轴（双链两根、单链一根）
    if style.double_strand {
        svg.line(left, axis_y - 2.0, right, axis_y - 2.0, &pal.backbone, 2.4);
        svg.line(left, axis_y + 2.0, right, axis_y + 2.0, &pal.backbone, 2.4);
    } else {
        svg.line(left, axis_y, right, axis_y, &pal.backbone, 2.0);
    }

    // 刻度 + 坐标（轴下）
    for (bp, major) in ruler_ticks(len, zoom, 56.0) {
        let x = x_of(bp);
        if x < left - 40.0 || x > right + 40.0 {
            continue;
        }
        let (c, w) = if major {
            (&pal.tick_major, 1.4)
        } else {
            (&pal.tick_minor, 1.0)
        };
        svg.line(x, axis_y, x, axis_y + if major { 8.0 } else { 4.0 }, c, w);
        if major {
            let label = thousands(if bp == 0.0 { 1 } else { bp as usize });
            svg.text(x, axis_y + 20.0, 10.5, Anchor::Middle, &pal.text_dim, &label, false);
        }
    }

    // 特征（轴上方带 + 方向箭头 + 名称）
    for seg in segments {
        let color = feature_fill(&seg.color, &pal.feature_fallback);
        for (rs, re) in split_wrap(seg.start, seg.end, len) {
            if rs >= re {
                continue;
            }
            let x0 = x_of(rs as f64);
            let x1 = x_of(re as f64);
            if x1 < left - 10.0 || x0 > right + 10.0 {
                continue;
            }
            let y1 = axis_y - LINE_RING * 0.15;
            let y0 = axis_y - LINE_RING * 1.05;
            svg.rect(
                x0,
                y0,
                (x1 - x0).max(1.0),
                y1 - y0,
                &color,
                Some((&color, 0.5)),
            );
            if seg.strand != 0 && (x1 - x0) > 14.0 {
                let size = ((x1 - x0) * 0.5).min(LINE_RING * 0.95);
                let (tip_x, base_x) = if seg.strand > 0 {
                    (x1, x1 - size)
                } else {
                    (x0, x0 + size)
                };
                let mid = (y0 + y1) / 2.0;
                svg.poly(
                    &[[tip_x, mid], [base_x, y0], [base_x, y1]],
                    &mix("#ffffff", &pal.bg, 0.0),
                    Some((&color, 1.0)),
                );
            }
            if (x1 - x0) > 46.0 {
                let label = if seg.name.is_empty() {
                    seg.kind.clone()
                } else {
                    seg.name.clone()
                };
                svg.text(x0 + 4.0, y0 - 5.0, 12.0, Anchor::Start, &color, &label, false);
            }
        }
    }

    // 酶位点（轴下）与标签
    let mut labels: Vec<(f32, f32, String, String)> = Vec::new();
    if style.show_sites {
        for s in sites {
            let x = x_of(s.pos as f64);
            if x < left - 30.0 || x > right + 30.0 {
                continue;
            }
            let color = if s.unique {
                pal.site_unique.clone()
            } else {
                pal.site.clone()
            };
            svg.line(x, axis_y + LINE_RING * 1.05, x, axis_y + LINE_RING * 1.6, &color, 1.6);
            if style.show_labels {
                let w = s.name.chars().count() as f32 * 6.2 + 6.0;
                labels.push((x, w, s.name.clone(), color));
            }
        }
    }
    if style.show_labels && !labels.is_empty() {
        let pos: Vec<f32> = labels.iter().map(|l| l.0).collect();
        let wid: Vec<f32> = labels.iter().map(|l| l.1).collect();
        let levels = assign_label_levels(&pos, &wid, 8);
        for (i, l) in labels.iter().enumerate() {
            let y = axis_y + LINE_RING * 1.6 + 12.0 + levels[i] as f32 * 13.0;
            svg.line(l.0, axis_y + LINE_RING * 1.6, l.0, y - 9.0, &l.3, 1.0);
            svg.text(l.0, y, 11.0, Anchor::Middle, &l.3, &l.2, false);
        }
    }

    // 端点标注
    svg.text(x_of(0.0), axis_y - LINE_RING * 1.9, 10.0, Anchor::Middle, &pal.text_dim, "1", false);
    svg.text(
        x_of(len as f64),
        axis_y - LINE_RING * 1.9,
        10.0,
        Anchor::Middle,
        &pal.text_dim,
        &thousands(len),
        false,
    );

    let sub = style
        .subtitle
        .clone()
        .unwrap_or_else(|| auto_subtitle(len, false));
    if title.is_some() {
        svg.text(16.0, 38.0, 10.5, Anchor::Start, &pal.text_dim, &sub, false);
    }
    svg.finish(&pal.bg)
}

/// 按拓扑分发的图谱 SVG。
pub fn map_svg(
    len: usize,
    segments: &[Segment],
    sites: &[SiteMark],
    style: &Style,
    topological: Topology,
    title: Option<&str>,
) -> String {
    let circular = topological == Topology::Circular;
    let ring = style
        .circular_layout
        .unwrap_or_else(|| crate::MapLayout::Auto.is_ring(circular, len));
    if !ring {
        return linear_map_svg(len, segments, sites, style, title);
    }
    // 线性分子按环状布局画 ⇒ 环上留缺口（副标题仍写真实拓扑）
    let mut st = style.clone();
    st.open_ring = !circular;
    circular_map_svg(len, segments, sites, &st, title)
}

/// 序列视图 SVG（60 bp/行、10 bp/组；高度按行数自动计算）。
///
/// `style.height` 忽略；返回的 SVG 高度 = 顶部留白 + 行数 × 行高 + 底部留白。
pub fn sequence_svg(
    seq: &[u8],
    segments: &[Segment],
    sites: &[SiteMark],
    style: &Style,
    title: Option<&str>,
) -> String {
    const LINE_BP: usize = 60;
    const GROUP_BP: usize = 10;
    const BASE_W: f32 = 11.0;
    const GROUP_GAP: f32 = 9.0;
    const LINE_H: f32 = 17.0;
    const PAD_LEFT: f32 = 84.0;
    const PAD_TOP: f32 = 32.0;
    const PAD_BOTTOM: f32 = 16.0;

    let pal = &style.palette;
    let len = seq.len();
    let rows = len.div_ceil(LINE_BP).max(1);
    let w = style.width;
    let h = PAD_TOP + rows as f32 * LINE_H + PAD_BOTTOM;
    let mut svg = Svg::new(w, h);

    let x_of = |bp: usize| -> f32 {
        let c = bp % LINE_BP;
        let g = c / GROUP_BP;
        PAD_LEFT + g as f32 * (GROUP_BP as f32 * BASE_W + GROUP_GAP) + (c % GROUP_BP) as f32 * BASE_W
    };
    let y_of = |row: usize| -> f32 { PAD_TOP + row as f32 * LINE_H };

    if let Some(t) = title {
        svg.text(12.0, 20.0, 13.0, Anchor::Start, &pal.text, t, true);
    }

    let last_row = rows.saturating_sub(1);
    for row in 0..rows {
        let start = row * LINE_BP;
        let end = (start + LINE_BP).min(len);
        let y = y_of(row);
        let base_y = y + LINE_H * 0.72;

        // 行号（该行起点，1-based）
        svg.text(
            PAD_LEFT - 8.0,
            base_y,
            10.5,
            Anchor::End,
            &pal.text_dim,
            &thousands(start + 1),
            false,
        );

        // 组间分隔
        let mut g = 1;
        while start + g * GROUP_BP < end {
            let x = x_of(start + g * GROUP_BP) - GROUP_GAP * 0.5;
            svg.line(x, y + 2.0, x, y + LINE_H - 2.0, &pal.tick_minor, 0.8);
            g += 1;
        }

        // 特征高亮（半透明近似：与白底混合 + 顶部实条）
        for seg in segments {
            for (rs, re) in split_wrap(seg.start, seg.end, len) {
                let a = rs.max(start);
                let b = re.min(end);
                if a >= b {
                    continue;
                }
                let color = feature_fill(&seg.color, &pal.feature_fallback);
                let x0 = x_of(a);
                let x1 = x_of(b - 1) + BASE_W;
                svg.rect(x0, y + 1.0, (x1 - x0).max(1.0), LINE_H - 2.0, &mix(&color, &pal.bg, 0.86), None);
                svg.rect(x0, y, (x1 - x0).max(1.0), 2.5, &color, None);
            }
        }

        // 碱基：按颜色分组，一组一个 `<text>`（x 列表逐个定位）
        let mut runs: [Vec<(f32, char)>; 5] = Default::default();
        for bp in start..end {
            let c = (seq[bp] as char).to_ascii_uppercase();
            let class = match c {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' | 'U' => 3,
                _ => 4,
            };
            runs[class].push((x_of(bp), c));
        }
        let colors = [
            pal.base_a.clone(),
            pal.base_c.clone(),
            pal.base_g.clone(),
            pal.base_t.clone(),
            pal.base_other.clone(),
        ];
        for (class, items) in runs.iter().enumerate() {
            if items.is_empty() {
                continue;
            }
            let xs: Vec<f32> = items.iter().map(|(x, _)| *x).collect();
            let chars: String = items.iter().map(|(_, c)| *c).collect();
            svg.text_run(&xs, base_y, 10.0, &colors[class], &chars, true);
        }

        // 酶位点下划线
        if style.show_sites {
            for s in sites {
                let a = s.pos.max(start);
                let b = (s.pos + s.site_len.max(1)).min(end);
                if a >= b {
                    continue;
                }
                let color = if s.unique {
                    pal.site_unique.clone()
                } else {
                    pal.site.clone()
                };
                let x0 = x_of(a);
                let x1 = x_of(b - 1) + BASE_W;
                svg.rect(x0, y + LINE_H - 3.0, (x1 - x0).max(1.0), 2.0, &color, None);
            }
        }

        // 行分隔（淡）
        if row != last_row {
            svg.line(PAD_LEFT - 60.0, y + LINE_H, w - 8.0, y + LINE_H, &mix(&pal.tick_minor, &pal.bg, 0.6), 0.6);
        }
    }

    svg.finish(&pal.bg)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(start: usize, end: usize, name: &str, strand: i8) -> Segment {
        Segment {
            start,
            end,
            name: name.to_string(),
            kind: "CDS".to_string(),
            color: "#22cc44".to_string(),
            strand,
        }
    }

    fn site(pos: usize, name: &str, unique: bool) -> SiteMark {
        SiteMark {
            pos,
            name: name.to_string(),
            site_len: 6,
            unique,
        }
    }

    /// 简易结构检查：所有 `<` 之后一定是标签（字母或 `/`），且以 `>` 收尾。
    fn well_formed(svg: &str) -> bool {
        let b = svg.as_bytes();
        let mut i = 0;
        let mut depth = 0i32;
        while i < b.len() {
            if b[i] == b'<' {
                if i + 1 >= b.len() || !(b[i + 1].is_ascii_alphabetic() || b[i + 1] == b'/') {
                    return false;
                }
                let close = b[i + 1] == b'/';
                let Some(rel) = svg[i..].find('>') else {
                    return false;
                };
                let tag = &svg[i..i + rel];
                if close {
                    depth -= 1;
                } else if !tag.ends_with('/') {
                    depth += 1;
                }
                i += rel + 1;
            } else {
                i += 1;
            }
        }
        depth == 0
    }

    #[test]
    fn circular_map_has_backbone_features_and_sites() {
        let segs = vec![seg(100, 400, "GFP", 1), seg(5000, 200, "wrap", -1)];
        let sites = vec![site(10, "EcoRI", true), site(2000, "BamHI", false)];
        let style = Style::map(900.0, 900.0);
        let svg = circular_map_svg(8000, &segs, &sites, &style, Some("pTEST"));
        assert!(svg.starts_with("<svg "), "{svg:.120}");
        assert!(svg.ends_with("</svg>"));
        assert!(well_formed(&svg), "结构不闭合");
        assert!(svg.contains(">pTEST<"), "标题缺失");
        assert!(svg.contains("8000 bp · circular"), "副标题缺失");
        assert!(svg.contains("EcoRI") && svg.contains("BamHI"), "酶标签缺失");
        assert!(svg.contains("#22cc44"), "特征色缺失");
        // 环内特征名（沿弧旋转，与界面一致）
        assert!(svg.contains(">GFP<"), "环状图应画特征名");
        assert!(svg.contains("transform=\"rotate("), "环上文字应带旋转");
        assert!(svg.matches("<path").count() >= 3, "扇环/箭头应有 path");
        assert!(svg.contains("#c4145a") && svg.contains("#8f1fa0"), "唯一/非唯一位点色");
    }

    #[test]
    fn dense_circular_svg_retains_nearby_sites_and_small_feature_names() {
        let sites=(0..45).map(|i|site(8900+i*8,&format!("Enzyme-{i}"),true)).collect::<Vec<_>>();
        let segments=(0..8).map(|i|seg(9200+i*10,9206+i*10,&format!("tag-{i}"),1)).collect::<Vec<_>>();
        let style=Style::map(1800.0,1100.0);
        let svg=circular_map_svg(10000,&segments,&sites,&style,Some("dense"));
        assert!(well_formed(&svg));assert!(!svg.contains("hidden (crowded)"));
        for i in 0..45 {assert!(svg.contains(&format!("Enzyme-{i} (")) || svg.contains(&format!(") Enzyme-{i}<")),"missing site {i}");}
        for i in 0..8 {assert!(svg.contains(&format!(">tag-{i}<")),"missing small feature {i}");}
        let mut hidden=style.clone();hidden.show_labels=false;
        let svg=circular_map_svg(10000,&segments,&sites,&hidden,None);
        assert!(!svg.contains("Enzyme-") && !svg.contains("tag-"));
    }

    #[test]
    fn linear_map_has_end_labels_and_ruler() {
        let segs = vec![seg(0, 600, "GFP", 1)];
        let style = Style::map(1000.0, 420.0);
        let svg = linear_map_svg(717, &segs, &[site(300, "NcoI", true)], &style, Some("GFP"));
        assert!(well_formed(&svg));
        assert!(svg.contains(">1<") && svg.contains(">717<"), "端点标注缺失");
        assert!(svg.contains(">GFP<"), "特征名缺失");
        assert!(svg.contains("NcoI"), "酶标签缺失");
        assert!(svg.contains("717 bp · linear"), "副标题缺失");
        assert!(svg.contains("<rect"), "特征带应为 rect");
    }

    #[test]
    fn sequence_view_rows_and_bases() {
        let seq: Vec<u8> = (0..130).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let segs = vec![seg(10, 40, "f", 1)];
        let style = Style::sequence(800.0);
        let svg = sequence_svg(&seq, &segs, &[site(61, "BamHI", true)], &style, Some("t"));
        assert!(well_formed(&svg));
        // 3 行（60/60/10）
        assert!(svg.contains("height=\"99\""), "高度按行数计算: {}", &svg[..160]);
        assert!(svg.contains(">1<") && svg.contains(">61<") && svg.contains(">121<"), "行号缺失");
        assert!(svg.contains(">BamHI<") == false, "序列视图不画酶名");
        assert!(svg.contains("#c4145a"), "位点下划线色");
        // 碱基 runs：A/C/G/T 各一段
        assert_eq!(svg.matches("font-family=\"Courier New,Courier,monospace\"").count(), 12, "3 行 × 4 色");
    }

    #[test]
    fn escaping_and_edge_cases() {
        let segs = vec![Segment {
            name: "A<B&C\"D".into(),
            kind: "x".into(),
            color: "not-a-color".into(),
            ..seg(0, 50, "", 0)
        }];
        let style = Style::map(400.0, 400.0);
        let svg = linear_map_svg(100, &segs, &[], &style, None);
        assert!(svg.contains("A&lt;B&amp;C&quot;D"), "文本未转义: {svg}");
        assert!(svg.contains(&Palette::default().feature_fallback), "非法色未兜底");
        assert!(well_formed(&svg));
        // 空序列不 panic
        let empty = circular_map_svg(0, &[], &[], &style, None);
        assert!(empty.ends_with("</svg>"));
        let empty_seq = sequence_svg(&[], &[], &[], &Style::sequence(400.0), None);
        assert!(empty_seq.ends_with("</svg>"));
    }

    #[test]
    fn thousands_and_mix() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        // SnapGene：4 位不分组，5 位才分组
        assert_eq!(thousands(1000), "1000");
        assert_eq!(thousands(9999), "9999");
        assert_eq!(thousands(25803), "25,803");
        assert_eq!(thousands(1234567), "1,234,567");
        assert_eq!(mix("#000000", "#ffffff", 0.5), "#808080");
        assert_eq!(mix("#ff0000", "#ffffff", 1.0), "#ffffff");
        assert_eq!(parse_hex("#00ff00"), [0, 255, 0]);
    }
}
