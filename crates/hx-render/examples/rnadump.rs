//! 临时工具：把一份二级结构布局画成 SVG，用于与 SnapGene 参考图目视比对。
use hx_render::rna_layout::*;

fn main() {
    let path = std::env::args().nth(1);
    let (len, pairs): (usize, Vec<(usize, usize)>) = match path {
        Some(p) => {
            let text = std::fs::read_to_string(&p).unwrap();
            let mut it = text.lines();
            let len: usize = it.next().unwrap().trim().parse().unwrap();
            let _n: usize = it.next().unwrap().trim().parse().unwrap();
            let pairs = it
                .filter_map(|l| {
                    let mut s = l.split_whitespace();
                    Some((s.next()?.parse().ok()?, s.next()?.parse().ok()?))
                })
                .collect();
            (len, pairs)
        }
        None => (
            72,
            vec![
                (0, 67), (1, 66), (2, 65), (3, 64), (4, 63), (5, 62), (6, 61), (9, 21), (10, 20),
                (11, 19), (12, 18), (22, 40), (23, 39), (24, 38), (25, 37), (26, 36), (27, 35),
                (44, 60), (45, 59), (46, 58), (47, 57), (48, 56),
            ],
        ),
    };
    let lay = RnaLayout::new(len, &pairs);
    let (min, max) = lay.bounds();
    let w = max.x - min.x + 40.0;
    let h = max.y - min.y + 40.0;
    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="#F2F2F2"/>"##
    );
    let tx = |p: Pt| ((p.x - min.x) + 20.0, (p.y - min.y) + 20.0);
    for i in 0..len.saturating_sub(1) {
        let a = tx(lay.pos[i]);
        let b = tx(lay.pos[i + 1]);
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#252525" stroke-width="2"/>"##,
            a.0, a.1, b.0, b.1
        ));
    }
    let seq = std::env::args().nth(2).unwrap_or_default();
    let seq = seq.as_bytes();
    for i in 0..len {
        let (x, y) = tx(lay.pos[i]);
        let c = if i < seq.len() { char::from(seq[i]) } else { '?' };
        svg.push_str(&format!(
            r##"<circle cx="{x:.1}" cy="{y:.1}" r="13.2" fill="white" stroke="#F80000" stroke-width="3.2"/><text x="{x:.1}" y="{y:.1}" font-family="Noto Sans" font-size="15" font-weight="bold" text-anchor="middle" dominant-baseline="central" fill="#252525">{c}</text>"##
        ));
    }
    for i in [0usize, 1, 69, 70, 71] {
        let lo = lay.loop_of[i];
        let c = lay.loops[lo].center;
        eprintln!("base {i}: pos ({:.1},{:.1}) loop {lo} center ({:.1},{:.1}) r {:.1}", lay.pos[i].x, lay.pos[i].y, c.x, c.y, lay.loops[lo].radius);
    }
    svg.push_str("</svg>");
    std::fs::write("/tmp/opencode/rna/out.svg", svg).unwrap();
    println!("{w} {h}");
}
