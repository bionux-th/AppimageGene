//! Screen-space callout layout for circular maps, shared by UI and SVG.
//!
//! Angles locate biological anchors, not text boxes. Horizontal text is packed
//! in circular order on each side of an elliptical rail. Isotonic regression
//! spreads crowded runs with minimum squared displacement, without swapping
//! labels. All label kinds must enter the SAME solve with measured box sizes.

#[derive(Debug, Clone, Copy)]
pub struct CircularLabel {
    /// Clockwise from twelve o'clock, in radians. Arbitrary rotations allowed.
    pub angle: f32,
    pub width: f32,
    pub height: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct Callout {
    /// Text/chip bounding box [left, top, right, bottom].
    pub rect: [f32; 4],
    /// Leader attaches to the near edge of the label, not through its text.
    pub port: [f32; 2],
    pub right: bool,
}

#[derive(Debug, Clone)]
pub struct CircularLabelLayout {
    pub center: [f32; 2],
    pub radius: f32,
    /// Same indices as the input. None means genuinely insufficient space.
    pub labels: Vec<Option<Callout>>,
    /// Ellipse enclosing leaders; labels are outside it, even at their corners.
    pub rail: [f32; 2],
}

const GAP: f32 = 1.0;

impl CircularLabelLayout {
    pub fn fit(bounds: [f32; 4], labels: &[CircularLabel], zoom: f32) -> Self {
        let [left, top, right, bottom] = bounds;
        let center = [(left + right) * 0.5, (top + bottom) * 0.5];
        let half_w = ((right - left) * 0.5).max(1.0);
        let half_h = ((bottom - top) * 0.5).max(1.0);
        if labels.is_empty() {
            let radius=((half_w.min(half_h)-4.0).max(1.0)*zoom.clamp(0.35,2.5)).min(half_w.min(half_h));
            return Self{center,radius,labels:Vec::new(),rail:[half_w,half_h]};
        }
        let mut sides = [Vec::new(), Vec::new()];
        let mut widths = [0.0_f32; 2];
        let mut demand = [0.0_f32; 2];
        for (i, label) in labels.iter().enumerate() {
            if !label.angle.is_finite() || !label.width.is_finite() || !label.height.is_finite()
                || label.width <= 0.0 || label.height <= 0.0 { continue; }
            let angle = crate::norm_angle(label.angle);
            let side = usize::from(angle <= std::f32::consts::PI);
            sides[side].push(i);
            widths[side] = widths[side].max(label.width);
            demand[side] += label.height + GAP;
        }
        // Reserve measured text width horizontally, and fan-out space according
        // to actual vertical demand (not a guessed number of radial levels).
        let max_width = widths[0].max(widths[1]);
        let density = (demand[0].max(demand[1]) / (2.0 * half_h)).min(1.0);
        let available_r = (half_w - max_width - 44.0).min(half_h - 44.0).max(12.0);
        // Leave a small radial clearance even for sparse maps; otherwise an
        // ellipse close to the poles could be inside the common leader stub.
        let desired=(half_h*(0.86-0.30*density)).min(available_r);
        let radius = (desired * zoom.clamp(0.35, 2.5)).clamp(12.0, available_r);
        let rail = [(radius + (half_h - radius).max(30.0) * 0.75)
            .min((half_w - max_width - 5.0).max(radius + 5.0)), half_h];
        let mut placed = vec![None; labels.len()];
        for (side, indices) in sides.iter_mut().enumerate() {
            // Top→bottom on each half circle. Angle, not cos(angle), breaks
            // near-pole ties deterministically and handles 0/2π correctly.
            indices.sort_by(|&a, &b| {
                let aa = crate::norm_angle(labels[a].angle);
                let ab = crate::norm_angle(labels[b].angle);
                let order = if side == 1 { aa.total_cmp(&ab) } else { ab.total_cmp(&aa) };
                order.then(a.cmp(&b))
            });
            // Preserve representative labels around the entire half-circle if
            // it is physically impossible to fit all of them at readable size.
            let height = indices.iter().map(|&i| labels[i].height).fold(0.0_f32, f32::max);
            let total_height=indices.iter().map(|&i|labels[i].height+GAP).sum::<f32>()-GAP;
            let capacity = if total_height<=2.0*half_h {indices.len()}else{
                (((2.0 * half_h + GAP) / (height + GAP)).floor() as usize).min(indices.len())};
            if capacity == 0 { continue; }
            let chosen: Vec<_> = (0..capacity).map(|k| indices[k * indices.len() / capacity]).collect();
            let mut offsets = Vec::with_capacity(chosen.len());
            let mut extent = 0.0;
            for &i in &chosen {
                offsets.push(extent + labels[i].height * 0.5);
                extent += labels[i].height + GAP;
            }
            extent -= GAP;
            let wanted: Vec<_> = chosen.iter().enumerate().map(|(k, &i)|
                -labels[i].angle.cos() * (half_h - labels[i].height * 0.5) - offsets[k]).collect();
            let values = isotonic(&wanted);
            let low = -half_h;
            let high = half_h - extent;
            for (k, &i) in chosen.iter().enumerate() {
                let y = values[k].clamp(low, high.max(low)) + offsets[k];
                // Near edge lies outside the ellipse over the entire box's
                // height. Thus a leader inside the ellipse cannot cut text.
                let near_y = (y.abs() - labels[i].height * 0.5).max(0.0);
                let x = rail[0] * (1.0 - (near_y / rail[1]).powi(2)).max(0.0).sqrt();
                let right_side = side == 1;
                let edge = center[0] + if right_side { x + 4.0 } else { -x - 4.0 };
                let (l, r) = if right_side { (edge, edge + labels[i].width) }
                    else { (edge - labels[i].width, edge) };
                let rect = [l, center[1] + y - labels[i].height * 0.5,
                    r, center[1] + y + labels[i].height * 0.5];
                if rect[0] < left || rect[2] > right || rect[1] < top - 0.01 || rect[3] > bottom + 0.01 { continue; }
                placed[i] = Some(Callout { rect, port: [edge, center[1] + y], right: right_side });
            }
        }
        Self { center, radius, labels: placed, rail }
    }

    /// A short radial stub followed by a fan-out leader ending before the text.
    /// Callers supply the real feature/primer radius rather than moving anchors.
    pub fn leader(&self, label: &CircularLabel, callout: &Callout, source_radius: f32) -> Vec<[f32; 2]> {
        let point = |r: f32| [self.center[0] + label.angle.sin() * r,
            self.center[1] - label.angle.cos() * r];
        let source = point(source_radius);
        // One common routing circle for features, enzymes AND both primer
        // tracks; different stub radii could reverse fan order near a cluster.
        let stub = point((self.radius + 34.0).max(source_radius + 5.0));
        let y = callout.port[1] - self.center[1];
        let rail_x = self.rail[0] * (1.0 - (y / self.rail[1]).powi(2)).max(0.0).sqrt();
        let edge = [self.center[0] + if callout.right { rail_x } else { -rail_x }, callout.port[1]];
        vec![source, stub, edge, [callout.port[0] + if callout.right { -2.0 } else { 2.0 }, callout.port[1]]]
    }
}

/// Pool-adjacent-violators: least-squares non-decreasing projection, O(n).
fn isotonic(values: &[f32]) -> Vec<f32> {
    let mut blocks: Vec<(f32, usize)> = Vec::new();
    for &v in values {
        blocks.push((v, 1));
        while blocks.len() >= 2 {
            let n = blocks.len();
            if blocks[n - 2].0 / blocks[n - 2].1 as f32 <= blocks[n - 1].0 / blocks[n - 1].1 as f32 { break; }
            let (sum, count) = blocks.pop().unwrap();
            blocks.last_mut().unwrap().0 += sum;
            blocks.last_mut().unwrap().1 += count;
        }
    }
    blocks.into_iter().flat_map(|(sum, count)| std::iter::repeat(sum / count as f32).take(count)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unannotated_rna_ring_uses_the_available_canvas_instead_of_callout_padding() {
        let layout=CircularLabelLayout::fit([46.0,62.0,1782.0,1040.0],&[],1.0);
        assert!(layout.radius>480.0);
        assert!(layout.labels.is_empty());
    }
    fn assert_clear(layout: &CircularLabelLayout, bounds: [f32; 4]) {
        let boxes: Vec<_> = layout.labels.iter().flatten().map(|p| p.rect).collect();
        for (i, a) in boxes.iter().enumerate() {
            assert!(a[0] >= bounds[0] && a[1] >= bounds[1] - 0.01 && a[2] <= bounds[2] && a[3] <= bounds[3] + 0.01, "{a:?}");
            for b in &boxes[..i] {
                assert!(a[2] <= b[0] || b[2] <= a[0] || a[3] <= b[1] || b[3] <= a[1], "{a:?} overlaps {b:?}");
            }
        }
    }
    #[test]
    fn dense_origin_cluster_and_primers_share_one_nonoverlapping_solve() {
        let mut labels: Vec<_> = (0..52).map(|i| CircularLabel { angle: 5.7 + i as f32 * 0.01, width: 50.0 + (i % 7) as f32 * 15.0, height: 14.0 }).collect();
        labels.extend((0..20).map(|i| CircularLabel { angle: 0.03 + i as f32 * 0.14, width: 160.0, height: 14.0 }));
        let bounds = [40.0, 60.0, 1760.0, 1040.0];
        let layout = CircularLabelLayout::fit(bounds, &labels, 1.0);
        assert_eq!(layout.labels.iter().flatten().count(), labels.len());
        assert_clear(&layout, bounds);
        assert!(layout.radius > 200.0 && layout.radius < 360.0);
        let mut left: Vec<_> = labels.iter().zip(&layout.labels).filter(|(_, p)| !p.unwrap().right).collect();
        left.sort_by(|a, b| b.0.angle.total_cmp(&a.0.angle));
        assert!(left.windows(2).all(|w| w[0].1.unwrap().port[1] < w[1].1.unwrap().port[1]));
    }
    #[test]
    fn rotations_zoom_resizes_and_top_bottom_seams_do_not_collide() {
        for rotation in [0.0, 0.001, 2.0, 3.14, 6.28] {
            let labels: Vec<_> = (0..65).map(|i| CircularLabel { angle: rotation + (i as f32 * 0.39), width: 55.0 + (i % 5) as f32 * 23.0, height: 14.0 }).collect();
            for (w, h) in [(1800.0, 1000.0), (1100.0, 740.0), (600.0, 400.0)] {
                let bounds = [0.0, 0.0, w, h];
                for zoom in [0.35, 1.0, 2.5] { assert_clear(&CircularLabelLayout::fit(bounds, &labels, zoom), bounds); }
            }
        }
    }

    #[test]
    fn mixed_text_and_chip_heights_remain_clear_and_deterministic() {
        let mut seed=73_u32;
        for _ in 0..120 {
            let mut labels=Vec::new();
            for _ in 0..85 {
                seed=seed.wrapping_mul(1664525).wrapping_add(1013904223);
                labels.push(CircularLabel{angle:seed as f32/u32::MAX as f32*8.0,
                    width:40.0+(seed%130) as f32,height:12.0+(seed%7) as f32});
            }
            let bounds=[15.0,30.0,1700.0,1000.0];
            let a=CircularLabelLayout::fit(bounds,&labels,1.0);
            let b=CircularLabelLayout::fit(bounds,&labels,1.0);
            assert_clear(&a,bounds);
            assert_eq!(a.labels.iter().map(|l|l.map(|p|p.rect)).collect::<Vec<_>>(),
                b.labels.iter().map(|l|l.map(|p|p.rect)).collect::<Vec<_>>());
        }
    }
    #[test]
    fn impossible_density_is_explicit_not_forced_into_overlapping_levels() {
        let labels = vec![CircularLabel { angle: 6.0, width: 110.0, height: 16.0 }; 400];
        let bounds = [0.0, 0.0, 600.0, 400.0];
        let layout = CircularLabelLayout::fit(bounds, &labels, 1.0);
        assert!(layout.labels.iter().any(Option::is_none));
        assert_clear(&layout, bounds);
    }

    #[test]
    fn fan_leaders_do_not_run_through_other_label_boxes() {
        let labels: Vec<_> = (0..60).map(|i| CircularLabel {
            angle: if i<40 {5.6+i as f32*0.015}else{(i-40) as f32*0.15},
            width: if i<40 {110.0}else{165.0}, height:14.0,
        }).collect();
        let layout=CircularLabelLayout::fit([0.0,0.0,1700.0,980.0],&labels,1.0);
        for (i,label) in labels.iter().enumerate() {
            let Some(callout)=layout.labels[i] else{continue};
            let line=layout.leader(label,&callout,layout.radius+4.0);
            for (j,other) in layout.labels.iter().enumerate() {
                let Some(other)=other else{continue};if i==j {continue;}
                for ends in line.windows(2) {
                    for k in 0..=100 {
                        let t=k as f32/100.0;
                        let p=[ends[0][0]*(1.0-t)+ends[1][0]*t,ends[0][1]*(1.0-t)+ends[1][1]*t];
                        assert!(!(p[0]>other.rect[0] && p[0]<other.rect[2] && p[1]>other.rect[1] && p[1]<other.rect[3]),"leader {i} crosses text {j}");
                    }
                }
            }
        }
    }

    #[test]
    fn sorted_fan_segments_do_not_cross_each_other() {
        for rotation in [0.0,0.5,2.0,3.14,5.0] {
        let labels:Vec<_>=(0..72).map(|i|CircularLabel{angle:rotation+if i<45 {4.0+i as f32*0.05}else{(i-45) as f32*0.12},width:150.0,height:16.0}).collect();
        let layout=CircularLabelLayout::fit([0.0,0.0,1800.0,1100.0],&labels,1.0);
        let cross=|a:[f32;2],b:[f32;2],c:[f32;2]|(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
        for i in 0..labels.len() {
            let Some(a)=layout.labels[i] else{continue};
            let a=layout.leader(&labels[i],&a,layout.radius+[-80.0,4.0,12.0,25.0][i%4]);
            for j in 0..i {
                let Some(b)=layout.labels[j] else{continue};let b=layout.leader(&labels[j],&b,layout.radius+[-80.0,4.0,12.0,25.0][j%4]);
                for aa in a.windows(2) {for bb in b.windows(2) {
                    assert!(!(cross(aa[0],aa[1],bb[0])*cross(aa[0],aa[1],bb[1]) < -0.01
                        && cross(bb[0],bb[1],aa[0])*cross(bb[0],bb[1],aa[1]) < -0.01),"leaders {i} and {j} cross");
                }}
            }
        }
        }
    }
}
