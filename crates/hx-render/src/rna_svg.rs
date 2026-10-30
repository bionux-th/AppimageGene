//! RNA 二级结构 → SVG（`Export as… ▸ Image…`）。
//!
//! 与屏幕上的画法一一对应：`#F2F2F2` 画布、空心彩色环 + 粗体字母、深色骨架、
//! 配对碱基之间的浅灰横档（氢键）、每 10 个碱基的序号标签、5′/3′ 深色小标签。
//! 坐标就是 `RnaLayout` 的 100% 坐标，只做「包围盒 → 画布」的平移。

use crate::rna_layout::{
    free_dir, Pt, RnaLayout, BACKBONE_WIDTH, LETTER_SIZE, RING_RADIUS, RING_STROKE,
};

/// 画布四周留白（给序号标签、5′/3′ 标签）。
const PAD: f32 = 46.0;

fn hex(c: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// 生成整张结构的 SVG。`colors` 与 `seq` 等长（每碱基一个 RGB）。
pub fn to_svg(
    layout: &RnaLayout,
    seq: &[u8],
    colors: &[[u8; 3]],
    pairs: &[(usize, usize)],
    show_index: bool,
) -> String {
    let len = layout.pos.len();
    let (bmin, bmax) = layout.bounds();
    let w = (bmax.x - bmin.x) + 2.0 * PAD;
    let h = (bmax.y - bmin.y) + 2.0 * PAD;
    let tf = |p: &Pt| {
        (
            p.x - bmin.x + PAD,
            p.y - bmin.y + PAD,
        )
    };
    let mut out = String::with_capacity(64 * len + 512);
    out.push_str(&format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w:.1}\" height=\"{h:.1}\" \
         viewBox=\"0 0 {w:.1} {h:.1}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#f2f2f2\"/>\n"
    ));

    // 碱基环
    for i in 0..len {
        let (x, y) = tf(&layout.pos[i]);
        let c = colors.get(i).copied().unwrap_or([0x25, 0x25, 0x25]);
        out.push_str(&format!(
            "<circle cx=\"{x:.2}\" cy=\"{y:.2}\" r=\"{RING_RADIUS}\" fill=\"none\" \
             stroke=\"{}\" stroke-width=\"{RING_STROKE}\"/>\n",
            hex(c)
        ));
    }
    // 字母
    let font = LETTER_SIZE.max(2.0);
    out.push_str(&format!(
        "<g font-family=\"Noto Sans, DejaVu Sans, sans-serif\" font-size=\"{font:.2}\" \
         font-weight=\"bold\" fill=\"#252525\" text-anchor=\"middle\" \
         dominant-baseline=\"central\">\n"
    ));
    for i in 0..len {
        let (x, y) = tf(&layout.pos[i]);
        // 与屏幕一致：一律大写，T 显示成 U。
        let ch = (match seq.get(i).copied().unwrap_or(b'?').to_ascii_uppercase() {
            b'T' => 'U',
            other => other as char,
        })
        .to_string();
        out.push_str(&format!(
            "<text x=\"{x:.2}\" y=\"{y:.2}\" dy=\"0.35\">{}</text>\n",
            esc(&ch)
        ));
    }
    out.push_str("</g>\n");

    // 骨架 + 配对横档：与屏幕一致，画在环**之上**、两端停在环内一点点；
    // 环是空心的，按外缘截断时描边会把缝隙吃光。
    let link = (RING_RADIUS - 1.0).max(RING_RADIUS * 0.5);
    out.push_str(&format!(
        "<g stroke=\"#252525\" stroke-width=\"{BACKBONE_WIDTH}\" \
         stroke-linecap=\"round\">\n"
    ));
    for i in 0..len.saturating_sub(1) {
        let (a, b) = (layout.pos[i], layout.pos[i + 1]);
        let d = b.sub(a);
        let l = d.len();
        if l <= 2.0 * link + 0.2 {
            continue;
        }
        let u = d.scale(1.0 / l);
        let (x1, y1) = tf(&a.add(u.scale(link)));
        let (x2, y2) = tf(&b.sub(u.scale(link)));
        out.push_str(&format!(
            "<line x1=\"{x1:.2}\" y1=\"{y1:.2}\" x2=\"{x2:.2}\" y2=\"{y2:.2}\"/>\n"
        ));
    }
    out.push_str("</g>\n");
    out.push_str(
        "<g stroke=\"#a8a8a8\" stroke-width=\"1.2\" stroke-linecap=\"round\">\n",
    );
    for &(i, j) in pairs {
        if i >= len || j >= len {
            continue;
        }
        let (a, b) = (layout.pos[i], layout.pos[j]);
        let d = b.sub(a);
        let l = d.len();
        if l <= 2.0 * link + 0.2 {
            continue;
        }
        let u = d.scale(1.0 / l);
        let (x1, y1) = tf(&a.add(u.scale(link)));
        let (x2, y2) = tf(&b.sub(u.scale(link)));
        out.push_str(&format!(
            "<line x1=\"{x1:.2}\" y1=\"{y1:.2}\" x2=\"{x2:.2}\" y2=\"{y2:.2}\"/>\n"
        ));
    }
    out.push_str("</g>\n");

    // 序号标签（与界面一致：挤在一起的跳过，避免把结构盖成一团）
    if show_index {
        out.push_str(
            "<g font-family=\"Noto Sans, DejaVu Sans, sans-serif\" font-size=\"11.5\" \
             fill=\"#252525\">\n",
        );
        let mut placed: Vec<Pt> = Vec::new();
        for i in 0..len {
            if i % 10 != 9 {
                continue;
            }
            let li = layout.loop_of.get(i).copied().unwrap_or(0);
            let c = layout.loops.get(li).map(|l| l.center).unwrap_or(layout.pos[i]);
            let d = free_dir(&layout.pos, i, layout.pos[i].sub(c).norm(), RING_RADIUS + 22.0);
            let tip = layout.pos[i].add(d.scale(RING_RADIUS + 5.0));
            let at = layout.pos[i].add(d.scale(RING_RADIUS + 18.0));
            if placed.iter().any(|q| q.sub(at).len() < 40.0) {
                continue;
            }
            placed.push(at);
            let (tx, ty) = tf(&tip);
            let (t2x, t2y) = tf(&tip.add(d.scale(9.0)));
            let (ax, ay) = tf(&at);
            out.push_str(&format!(
                "<line x1=\"{tx:.2}\" y1=\"{ty:.2}\" x2=\"{t2x:.2}\" y2=\"{t2y:.2}\" \
                 stroke=\"#8a8a8a\" stroke-width=\"0.8\"/>\n"
            ));
            let anchor = if d.x > 0.25 {
                "start"
            } else if d.x < -0.25 {
                "end"
            } else {
                "middle"
            };
            out.push_str(&format!(
                "<text x=\"{ax:.2}\" y=\"{ay:.2}\" text-anchor=\"{anchor}\" \
                 dominant-baseline=\"central\">{}</text>\n",
                i + 1
            ));
        }
        out.push_str("</g>\n");
    }

    // 5′ / 3′
    if len >= 2 {
        let centroid = layout
            .pos
            .iter()
            .fold(Pt::new(0.0, 0.0), |a, b| a.add(*b))
            .scale(1.0 / len.max(1) as f32);
        for (idx, _other, text) in [(0usize, 1usize, "5'"), (len - 1, len - 2, "3'")] {
            let away = layout.pos[idx].sub(centroid).norm();
            let d = free_dir(&layout.pos, idx, away, RING_RADIUS + 24.0);
            let at = layout.pos[idx].add(d.scale(RING_RADIUS + 24.0));
            let (x, y) = tf(&at);
            out.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"22\" height=\"16\" rx=\"3.5\" \
                 fill=\"#3a3a3a\"/>\n\
                 <text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"Noto Sans, sans-serif\" \
                 font-size=\"11.5\" font-weight=\"bold\" fill=\"#ffffff\" \
                 text-anchor=\"middle\" dominant-baseline=\"central\">{text}</text>\n",
                x - 11.0,
                y - 8.0
            ));
        }
    }
    out.push_str("</svg>\n");
    out
}
