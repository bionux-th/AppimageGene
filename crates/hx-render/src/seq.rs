//! 序列视图排版：碱基 → 行列坐标、分组间隙、可见行裁剪、区间横向像素跨度。
//!
//! 默认排版对齐 SnapGene 的序列视图习惯：每行 60 碱基、每 10 碱基一组，
//! 组间留窄缝；行首标注该行起点坐标（1-based）。

use std::ops::Range;

/// 序列视图布局（屏幕坐标 y 向下；`scroll_y` 为像素滚动量）。
#[derive(Debug, Clone, Copy)]
pub struct SeqView {
    /// 每行碱基数（默认 60）。
    pub line_bp: usize,
    /// 每组的碱基数（默认 10）。
    pub group_bp: usize,
    /// 单个碱基格宽度（像素）。
    pub base_w: f32,
    /// 组间额外间隙（像素）。
    pub group_gap: f32,
    /// 行高（像素）。
    pub line_h: f32,
    /// 左侧留白（行号/坐标标注区）。
    pub pad_left: f32,
    /// 顶部留白。
    pub pad_top: f32,
    /// 垂直滚动（像素，向下为正）。
    pub scroll_y: f32,
}

impl Default for SeqView {
    fn default() -> Self {
        Self {
            line_bp: 60,
            group_bp: 10,
            base_w: 13.0,
            group_gap: 8.0,
            line_h: 20.0,
            pad_left: 76.0,
            pad_top: 28.0,
            scroll_y: 0.0,
        }
    }
}

impl SeqView {
    /// 第 `g` 组（0 起）的起点 x。
    #[inline]
    fn group_x(&self, g: usize) -> f32 {
        self.pad_left + g as f32 * (self.group_bp as f32 * self.base_w + self.group_gap)
    }

    /// 碱基在行内的列号。
    #[inline]
    pub fn col_of(&self, bp: usize) -> usize {
        bp % self.line_bp.max(1)
    }

    /// 碱基所在行号（屏幕行，未计滚动）。
    #[inline]
    pub fn row_of(&self, bp: usize) -> usize {
        bp / self.line_bp.max(1)
    }

    /// 碱基单元格左上角 x。
    pub fn x_of(&self, bp: usize) -> f32 {
        let c = self.col_of(bp);
        let g = c / self.group_bp.max(1);
        let within = c % self.group_bp.max(1);
        self.group_x(g) + within as f32 * self.base_w
    }

    /// 第 `row` 行的基线 y（含 `scroll_y`；即该行碱基文本的垂直中心）。
    pub fn y_of(&self, row: usize) -> f32 {
        self.pad_top + row as f32 * self.line_h - self.scroll_y
    }

    /// 碱基单元格左上角（含滚动）。
    pub fn xy_of(&self, bp: usize) -> [f32; 2] {
        [self.x_of(bp), self.y_of(self.row_of(bp))]
    }

    /// 命中测试：屏幕点 → 碱基序号（超出碱基区返回 `None`）。
    pub fn bp_at(&self, len: usize, p: [f32; 2]) -> Option<usize> {
        let row_f = (p[1] - self.pad_top + self.scroll_y) / self.line_h;
        if row_f < 0.0 {
            return None;
        }
        let row = row_f as usize;
        let row_start = row * self.line_bp;
        if row_start >= len {
            return None;
        }
        if p[0] < self.pad_left {
            return None;
        }
        // 反解组号：先按“无间隙”估算，再逐组修正
        let stride = self.group_bp.max(1) as f32 * self.base_w + self.group_gap;
        let mut g = (((p[0] - self.pad_left) / stride) as usize).min(self.line_bp / self.group_bp.max(1));
        // 修正到包含该 x 的组
        while g > 0 && self.group_x(g) > p[0] {
            g -= 1;
        }
        while self.group_x(g + 1) <= p[0] && (g + 1) * self.group_bp < self.line_bp {
            g += 1;
        }
        let local = p[0] - self.group_x(g);
        if local < 0.0 || local >= self.group_bp as f32 * self.base_w {
            return None; // 落在组间缝隙里
        }
        let within = (local / self.base_w) as usize;
        let bp = row_start + g * self.group_bp + within;
        (bp < len && bp < (row + 1) * self.line_bp).then_some(bp)
    }

    /// 可见行范围（含少量上下余量）。
    pub fn visible_rows(&self, len: usize, viewport_h: f32) -> Range<usize> {
        let total = self.row_count(len);
        if total == 0 {
            return 0..0;
        }
        let first = ((self.scroll_y - self.pad_top) / self.line_h).floor().max(0.0) as usize;
        let last = ((self.scroll_y + viewport_h - self.pad_top) / self.line_h).ceil().max(0.0) as usize + 1;
        first.min(total)..last.min(total)
    }

    /// 总行数。
    pub fn row_count(&self, len: usize) -> usize {
        len.div_ceil(self.line_bp.max(1))
    }

    /// 内容总高度（像素）。
    pub fn content_height(&self, len: usize) -> f32 {
        self.pad_top + self.row_count(len) as f32 * self.line_h
    }

    /// 区间 `[start, end)` 在第 `row` 行的横向跨度（像素）；该行无交集时返回空。
    pub fn spans_in_row(&self, start: usize, end: usize, row: usize) -> Vec<(f32, f32)> {
        let row_start = row * self.line_bp;
        let row_end = row_start + self.line_bp;
        let s = start.max(row_start);
        let e = end.min(row_end);
        if s >= e {
            return Vec::new();
        }
        // 行尾（e == row_end）取该行碱基区右缘，而不是下一行首 x
        let x1 = if e >= row_end {
            self.x_of(row_end - 1) + self.base_w
        } else {
            self.x_of(e)
        };
        vec![(self.x_of(s), x1)]
    }

    /// 把（可跨原点的）区间铺成 `(row, x0, x1)` 列表（用于特征/位点高亮条）。
    ///
    /// `start == end` 视为整圈（调用方对环状分子使用）。
    pub fn rows_of_span(&self, start: usize, end: usize, len: usize) -> Vec<(usize, f32, f32)> {
        let mut out = Vec::new();
        for (s, e) in crate::split_wrap(start, end, len) {
            if s >= e {
                continue;
            }
            for row in self.row_of(s)..=self.row_of(e.saturating_sub(1)) {
                out.extend(
                    self.spans_in_row(s, e, row)
                        .into_iter()
                        .map(|(x0, x1)| (row, x0, x1)),
                );
            }
        }
        out.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        });
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_and_col() {
        let v = SeqView::default();
        assert_eq!(v.row_of(0), 0);
        assert_eq!(v.col_of(0), 0);
        assert_eq!(v.row_of(59), 0);
        assert_eq!(v.row_of(60), 1);
        assert_eq!(v.col_of(60), 0);
        assert_eq!(v.col_of(65), 5);
    }

    #[test]
    fn groups_have_gaps() {
        let v = SeqView::default();
        // 组 0 起点 = pad_left；组 1 起点 = pad_left + 10*13 + 8
        assert!((v.x_of(0) - 76.0).abs() < 1e-4);
        assert!((v.x_of(10) - (76.0 + 138.0)).abs() < 1e-4);
        assert!((v.x_of(11) - (76.0 + 138.0 + 13.0)).abs() < 1e-4);
        // 组内相邻只差 base_w
        assert!((v.x_of(5) - v.x_of(4) - 13.0).abs() < 1e-4);
    }

    #[test]
    fn hit_test_round_trip() {
        let v = SeqView::default();
        for bp in [0usize, 7, 9, 10, 33, 59, 60, 61, 119, 259] {
            let p = v.xy_of(bp);
            let got = v.bp_at(1000, [p[0] + 1.0, p[1] + 1.0]);
            assert_eq!(got, Some(bp), "bp {bp} → {got:?} (p={p:?})");
        }
    }

    #[test]
    fn hit_test_gap_is_none() {
        let v = SeqView::default();
        // 组间缝隙（组 0 右缘之后、组 1 起点之前）
        let x = 76.0 + 10.0 * 13.0 + 2.0;
        let y = v.y_of(0) + 2.0;
        assert_eq!(v.bp_at(1000, [x, y]), None);
    }

    #[test]
    fn hit_test_outside() {
        let v = SeqView::default();
        assert_eq!(v.bp_at(100, [10.0, 30.0]), None); // 左留白
        assert_eq!(v.bp_at(100, [80.0, 5.0]), None); // 顶部之上
        assert_eq!(v.bp_at(10, [80.0, v.y_of(1) + 2.0]), None); // 行 1 超出长度
    }

    #[test]
    fn visible_rows_scrolling() {
        let v = SeqView::default();
        assert_eq!(v.row_count(600), 10);
        let r = v.visible_rows(600, 200.0);
        assert_eq!(r.start, 0);
        assert!(r.end >= 9);
        let mut v2 = v;
        v2.scroll_y = 4.0 * 20.0;
        let r2 = v2.visible_rows(600, 200.0);
        assert!(r2.start >= 2 && r2.start <= 4, "{r2:?}");
    }

    #[test]
    fn spans_single_row() {
        let v = SeqView::default();
        // 特征 5..15（行 0 内跨组）
        let sp = v.spans_in_row(5, 15, 0);
        assert_eq!(sp.len(), 1);
        assert!((sp[0].0 - v.x_of(5)).abs() < 1e-4);
        assert!((sp[0].1 - v.x_of(15)).abs() < 1e-4);
        // 行 1 无交集
        assert!(v.spans_in_row(5, 15, 1).is_empty());
    }

    #[test]
    fn spans_row_end_clamps_to_base_right_edge() {
        let v = SeqView::default();
        // 55..60：终点在行尾，右缘应为最后一格右缘（不含组间缝）
        let sp = v.spans_in_row(55, 60, 0);
        assert_eq!(sp.len(), 1);
        assert!((sp[0].1 - (v.x_of(59) + v.base_w)).abs() < 1e-4);
    }

    #[test]
    fn spans_multi_row() {
        let v = SeqView::default();
        let rows: Vec<usize> = v.rows_of_span(50, 70, 1000).iter().map(|r| r.0).collect();
        assert_eq!(rows, vec![0, 1]);
        // 行 0：50..60；行 1：60..70
        let s0 = v.rows_of_span(50, 70, 1000);
        assert!((s0[0].1 - v.x_of(50)).abs() < 1e-4);
        assert!((s0[0].2 - (v.x_of(59) + v.base_w)).abs() < 1e-4);
        assert!((s0[1].1 - v.x_of(60)).abs() < 1e-4);
        assert!((s0[1].2 - v.x_of(70)).abs() < 1e-4);
    }

    #[test]
    fn spans_wrapping_feature() {
        let v = SeqView::default();
        // 长度 100、特征 90..10（跨原点）→ 行 0 的 0..10 + 行 1 的 90..100
        let rows = v.rows_of_span(90, 10, 100);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].0, 0);
        assert_eq!(rows[1].0, 1);
    }

    #[test]
    fn content_height_and_scroll_clamp() {
        let v = SeqView::default();
        assert!(v.content_height(0) > 0.0);
        let h = v.content_height(600);
        assert!((h - (28.0 + 10.0 * 20.0)).abs() < 1e-4);
    }
}
