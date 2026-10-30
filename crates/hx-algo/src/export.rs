//! # SVG 导出胶水（UI 无关）
//!
//! 把编辑会话的**工作状态**（序列/拓扑/特征）与酶位点命中翻译成
//! [`hx_render::svg`] 的输入并产出 SVG 字符串。UI（`hx-app` 的 SVG 按钮）
//! 与 CLI（`hx-cli export`）共用这一条路径，保证两边导出结果一致。
//!
//! 坐标：特征/位点均为 0-based；导出器内部按拓扑处理环状跨原点区间。

use hx_render::svg::{Segment, SiteMark, Style};
use hx_render::Topology;

use crate::EditSession;

/// 一次导出的产物与统计。
#[derive(Debug, Clone)]
pub struct Export {
    pub svg: String,
    /// 画出的特征片段数（多片段特征按段计）。
    pub bands: usize,
    /// 画出的酶切位点数。
    pub sites: usize,
}

/// 工作状态里的特征 → 渲染片段（`directionality` 0/1/2/3 → strand 0/1/-1/1）。
pub fn segments(session: &EditSession) -> Vec<Segment> {
    let mut out = Vec::new();
    for f in session.features() {
        let strand = if f.directionality == 0 {
            0
        } else if f.is_forward() {
            1
        } else {
            -1
        };
        for seg in &f.segments {
            out.push(Segment {
                start: seg.range.start,
                end: seg.range.end,
                name: f.name.clone(),
                kind: f.kind.clone(),
                color: seg.color.clone(),
                strand,
            });
        }
    }
    out
}

/// 酶下标 → 位点命中（0-based 起点；与 `hx-app` 的 [`SiteHits`] 语义一致）。
pub fn find_hits(seq: &[u8], circular: bool, enzymes: &[usize]) -> Vec<Vec<usize>> {
    enzymes
        .iter()
        .map(|&i| match hx_enzymes::ENZYMES.get(i) {
            Some(e) => hx_enzymes::find_site_positions(seq, e.site.as_bytes(), circular)
                .into_iter()
                .map(|p| p as usize)
                .collect(),
            None => Vec::new(),
        })
        .collect()
}

/// 命中 → 渲染标记（按坐标升序，供标签避让算法使用）。
pub fn site_marks(hits: &[Vec<usize>], enzymes: &[usize]) -> Vec<SiteMark> {
    let mut out: Vec<SiteMark> = Vec::new();
    for (k, &i) in enzymes.iter().enumerate() {
        let Some(e) = hx_enzymes::ENZYMES.get(i) else {
            continue;
        };
        let unique = hits.get(k).map(|v| v.len() == 1).unwrap_or(false);
        let site_len = e.site_len();
        for &pos in hits.get(k).map(|v| v.as_slice()).unwrap_or(&[]) {
            out.push(SiteMark {
                pos,
                name: e.name.to_string(),
                site_len,
                unique,
            });
        }
    }
    out.sort_by_key(|s| s.pos);
    out
}

/// 图谱 SVG（按会话拓扑自动选环状/线性）。
pub fn map_svg(session: &EditSession, enzymes: &[usize], style: &Style, title: &str) -> Export {
    let bands = segments(session);
    let hits = find_hits(session.sequence(), session.circular(), enzymes);
    let marks = site_marks(&hits, enzymes);
    let topo = if session.circular() {
        Topology::Circular
    } else {
        Topology::Linear
    };
    // 布局：调用方给了就用（界面用当前显示设置），否则按 Auto 规则
    let mut style = style.clone();
    if style.circular_layout.is_none() {
        style.circular_layout = Some(auto_ring(session));
    }
    style.open_ring = !session.circular() && style.circular_layout == Some(true);
    // 单链（RNA / 单链 DNA / 蛋白）只有一根骨架线，与界面一致。
    style.double_strand = session.kind() == hx_core::MoleculeKind::Dna
        && session.flags().strandedness() == hx_core::Strandedness::Double;
    let svg = hx_render::svg::map_svg(
        session.len(),
        &bands,
        &marks,
        &style,
        topo,
        Some(title),
    );
    Export {
        svg,
        bands: bands.len(),
        sites: marks.len(),
    }
}

/// 序列视图 SVG（行数由序列长度决定）。
pub fn sequence_svg(session: &EditSession, enzymes: &[usize], style: &Style, title: &str) -> Export {
    let bands = segments(session);
    let hits = find_hits(session.sequence(), session.circular(), enzymes);
    let marks = site_marks(&hits, enzymes);
    // `.rna` 里存的是 DNA 字母，导出图与屏幕一致显示 `U`。
    let display: Vec<u8> = session
        .sequence()
        .iter()
        .map(|&b| session.kind().display_base(b))
        .collect();
    let svg = hx_render::svg::sequence_svg(
        &display,
        &bands,
        &marks,
        style,
        Some(title),
    );
    Export {
        svg,
        bands: bands.len(),
        sites: marks.len(),
    }
}

/// `hx-cli` 与 UI 共用的默认宽度（图谱 900×900（环状）/ 1200×460（线性），序列 900）。
pub fn default_map_style(ring: bool) -> Style {
    let mut style = if ring {
        Style::map(900.0, 900.0)
    } else {
        Style::map(1200.0, 460.0)
    };
    // 与界面一致的布局决定（Auto 规则由调用方用 `auto_ring` 解析）
    style.circular_layout = Some(ring);
    style
}

/// `MapLayout::Auto` 下是否画成环（线性分子 >1 kb 也画环，与 SnapGene 一致）。
pub fn auto_ring(session: &EditSession) -> bool {
    hx_render::MapLayout::Auto.is_ring(session.circular(), session.len())
}

/// 序列视图默认样式。
pub fn default_sequence_style() -> Style {
    Style::sequence(900.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use hx_formats::write::NewFile;
    use hx_formats::DnaFile;

    /// 22 bp 骨架（`NewFile::dna` 默认**线性**，flags = `DOUBLE_STRANDED`）；
    /// 特征 3..8（0-based、1-based 闭 4-8）、`GGGCCC`（ApaI）在 11 处。
    fn build(circular: bool) -> EditSession {
        let seq = b"ACGTACGTACGGGGCCCCTTTT".to_vec();
        let bytes0 = NewFile::dna(seq).to_bytes();
        let mut f = DnaFile::parse(&bytes0).unwrap();
        let flags = hx_core::MoleculeFlags::from_byte(if circular {
            hx_core::MoleculeFlags::DOUBLE_STRANDED | hx_core::MoleculeFlags::CIRCULAR
        } else {
            hx_core::MoleculeFlags::DOUBLE_STRANDED
        })
        .to_byte();
        let feats = "<?xml version=\"1.0\"?><Features nextValidID=\"1\"><Feature recentID=\"0\" name=\"f1\" directionality=\"2\" type=\"CDS\"><Segment range=\"3-8\" color=\"#ff0000\" type=\"standard\" translated=\"0\"/></Feature></Features>";
        for b in f.blocks.iter_mut() {
            match b.block_type {
                0 => b.payload[0] = flags,
                10 => b.payload = feats.as_bytes().to_vec(),
                _ => {}
            }
        }
        let bytes = f.to_bytes();
        let parsed = DnaFile::parse(&bytes).unwrap();
        EditSession::new(parsed, bytes).unwrap()
    }

    /// 线性会话。
    fn sample() -> EditSession {
        let s = build(false);
        assert!(!s.circular());
        s
    }

    /// 环状会话。
    fn circular_session() -> EditSession {
        let s = build(true);
        assert!(s.circular());
        s
    }

    fn enz(names: &[&str]) -> Vec<usize> {
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
        v
    }

    #[test]
    fn segments_carry_strand_and_color() {
        let s = sample();
        let segs = segments(&s);
        assert_eq!(segs.len(), 1);
        assert_eq!((segs[0].start, segs[0].end), (2, 8));
        assert_eq!(segs[0].name, "f1");
        assert_eq!(segs[0].kind, "CDS");
        assert_eq!(segs[0].strand, -1, "directionality=2 → 反向");
        assert_eq!(segs[0].color, "#ff0000");
    }

    #[test]
    fn sticky_exports_use_top_coordinates_not_hidden_marker_bases() {
        let bytes = hx_formats::write::NewFile::dna(b"TGGGGAATTCCCCAA".to_vec()).to_bytes();
        let mut file = hx_formats::DnaFile::parse(&bytes).unwrap();
        file.blocks.iter_mut().find(|b| b.block_type == 8).unwrap().payload =
            b"<AdditionalSequenceProperties><UpstreamStickiness>-1</UpstreamStickiness><DownstreamStickiness>2</DownstreamStickiness></AdditionalSequenceProperties>".to_vec();
        file.blocks.iter_mut().find(|b| b.block_type == 10).unwrap().payload =
            b"<Features nextValidID=\"1\"><Feature recentID=\"0\" name=\"test\" type=\"misc_feature\"><Segment range=\"2-5\" color=\"#ff0000\"/></Feature></Features>".to_vec();
        let bytes = file.to_bytes();
        let s = EditSession::new(file, bytes).unwrap();
        assert_eq!(s.len(), 12);
        let bands = segments(&s);
        assert_eq!((bands[0].start, bands[0].end), (0, 4));
        let indices = enz(&["EcoRI"]);
        assert_eq!(find_hits(s.sequence(), false, &indices), vec![vec![3]]);
        let exported = sequence_svg(&s, &[], &default_sequence_style(), "sticky");
        let expected = hx_render::svg::sequence_svg(b"GGGGAATTCCCC", &bands, &[], &default_sequence_style(), Some("sticky"));
        assert_eq!(exported.svg, expected);
    }

    #[test]
    fn hits_and_marks_track_sites() {
        let s = sample();
        let e = enz(&["ApaI", "EcoRI"]);
        assert_eq!(e.len(), 2, "酶名需在库中");
        let hits = find_hits(s.sequence(), s.circular(), &e);
        let marks = site_marks(&hits, &e);
        // ApaI = GGGCCC 命中 1 处（唯一）；EcoRI 无命中
        assert_eq!(marks.len(), 1, "{marks:?}");
        assert_eq!(marks[0].pos, 11, "GGGCCC 起点（0-based）");
        assert!(marks[0].unique);
        assert_eq!(marks[0].site_len, 6);
        assert_eq!(marks[0].name, "ApaI");
        // 命中数校验（EcoRI 为空）
        assert!(hits[0].is_empty() || hits[1].is_empty());
    }

    #[test]
    fn map_export_is_circular_with_site_label() {
        let s = circular_session();
        let e = enz(&["ApaI"]);
        let ex = map_svg(&s, &e, &default_map_style(s.circular()), "t.dna");
        assert_eq!(ex.bands, 1);
        assert_eq!(ex.sites, 1);
        assert!(ex.svg.starts_with("<svg "), "{}", &ex.svg[..80.min(ex.svg.len())]);
        assert!(ex.svg.ends_with("</svg>"));
        assert!(ex.svg.contains(">t.dna<"), "标题缺失");
        assert!(ex.svg.contains("22 bp · circular"), "副标题缺失");
        assert!(ex.svg.contains("ApaI"), "酶标签缺失");
        assert!(ex.svg.contains("#ff0000"), "特征色缺失");
        assert!(!ex.svg.contains("<rect x=\"0\" y=\"0\" width=\"0\""), "尺寸异常");
    }

    #[test]
    fn sequence_export_has_rows_and_underline() {
        let s = sample();
        let e = enz(&["ApaI"]);
        let ex = sequence_svg(&s, &e, &default_sequence_style(), "t.dna");
        assert_eq!(ex.bands, 1);
        assert_eq!(ex.sites, 1);
        assert!(ex.svg.ends_with("</svg>"));
        // 22 bp → 1 行：4 个碱基色 run + 行号 + 标题
        assert_eq!(
            ex.svg
                .matches("font-family=\"Courier New,Courier,monospace\"")
                .count(),
            4,
            "碱基 run 数"
        );
        assert!(ex.svg.contains(">1<"), "行号缺失");
        assert!(ex.svg.contains(">t.dna<"), "标题缺失");
        assert!(ex.svg.contains("#c4145a"), "唯一酶位点下划线色");
    }

    #[test]
    fn linear_export_for_linear_session() {
        let s = sample();
        let ex = map_svg(&s, &[], &default_map_style(false), "lin");
        assert_eq!(ex.sites, 0);
        assert!(ex.svg.contains("22 bp · linear"), "线性副标题");
        assert!(ex.svg.contains(">1<") && ex.svg.contains(">22<"), "端点标注");
    }
}
