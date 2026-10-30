//! Main-window geometry measured from the supplied expanded/collapsed SG views.
//! Keep layout decisions testable rather than scattering offsets in painters.
pub const SIDEBAR_W:f32=230.0;
pub const MENU_H:f32=22.0;
pub const DOC_TABS_H:f32=34.0;
pub const TOOLBAR_H:f32=50.0;
pub const VIEW_TAB_H:f32=28.0;
pub const VIEW_TAB_FONT:f32=12.5;

/// Right inset shared by the sidebar `⋮` trigger and the Home icon above it, so
/// both centres sit in one column (user request 2026-10-29).
pub const NAV_RIGHT_INSET:f32=6.0;

pub fn navigation_width(expanded:bool,sidebar_width:f32)->f32 {
    if expanded {sidebar_width}else{88.0}
}

pub fn document_tab_width(available:f32,count:usize)->f32 {
    if count==0 {return 196.0;}
    ((available-68.0).max(1.0)/count as f32).clamp(62.0,196.0)
}

/// Horizontal separator leaves exactly one gap under the active document.
/// Navigation has its own bottom rule only when the project panel is open.
pub fn tab_separator_segments(left:f32,right:f32,active:Option<(f32,f32)>)->Vec<(f32,f32)> {
    match active {
        Some((a,b))=>[(left,a.max(left).min(right)),(b.max(left).min(right),right)]
            .into_iter().filter(|(a,b)|b>a).collect(),
        None=>vec![(left,right)],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn separator_spans_unused_row_and_inactive_tabs_but_not_active_tab() {
        assert_eq!(tab_separator_segments(0.0,2000.0,Some((230.0,426.0))),vec![(0.0,230.0),(426.0,2000.0)]);
        assert_eq!(tab_separator_segments(0.0,2000.0,Some((622.0,818.0))),vec![(0.0,622.0),(818.0,2000.0)]);
        assert_eq!(tab_separator_segments(0.0,2000.0,None),vec![(0.0,2000.0)]);
    }
    #[test]
    fn document_tabs_keep_standard_width_until_the_row_is_full() {
        assert_eq!(document_tab_width(1900.0,6),196.0);
        assert!(document_tab_width(1000.0,6)<196.0);
        assert_eq!(navigation_width(false,270.0),88.0);
        assert_eq!(navigation_width(true,270.0),270.0);
    }

    #[test]
    fn crowded_tabs_reserve_space_for_new_and_overflow_controls() {
        for count in [1,3,6,14,30] {
            let width=document_tab_width(1700.0,count);
            assert!(width>=62.0 && width<=230.0);
            if count<=14 {assert!(width*count as f32<=1700.0-68.0+0.1);}
        }
    }
}
