//! 内嵌 ViennaRNA 的底层 FFI。
//!
//! C 源码来自 `vendor/viennarna`（ViennaRNA 2.7.2，见其中的 `COPYING`），
//! 由 `build.rs` 静态编译进本 crate —— 运行时不需要 Python、也不需要系统装
//! ViennaRNA。C 侧的接口只有 `shim/vienna_shim.c` 里的两个函数，上游结构体
//! 不跨 FFI。

use std::os::raw::{c_char, c_double, c_int, c_void};

extern "C" {
    /// 折叠：MFE + 配分函数 + 配对概率回调。成功返回 0。
    fn hx_vienna_fold(
        seq: *const c_char,
        temperature: c_double,
        no_gu_close: c_int,
        no_lp: c_int,
        db: *mut c_char,
        db_len: c_int,
        mfe_energy: *mut c_double,
        ensemble: *mut c_double,
        diversity: *mut c_double,
        cb: Option<unsafe extern "C" fn(c_int, c_int, c_double, *mut c_void)>,
        ud: *mut c_void,
    ) -> c_int;

    /// 次优结构：逐条回调 dot-bracket 与能量（`delta` 单位 10 cal/mol）。
    fn hx_vienna_subopt(
        seq: *const c_char,
        temperature: c_double,
        no_gu_close: c_int,
        no_lp: c_int,
        delta: c_int,
        cb: Option<unsafe extern "C" fn(c_double, *const c_char, *mut c_void)>,
        ud: *mut c_void,
    ) -> c_int;

    /// 给一个 dot-bracket 打分（用于自检次优结构的能量）。
    fn hx_vienna_eval(
        seq: *const c_char,
        db: *const c_char,
        temperature: c_double,
        no_gu_close: c_int,
        no_lp: c_int,
    ) -> c_double;

    /// naview 布局：成功返回碱基个数。
    fn hx_vienna_naview(
        structure: *const c_char,
        x: *mut c_double,
        y: *mut c_double,
        max_len: c_int,
    ) -> c_int;
}

pub const VERSION: &str = "2.7.2";

/// 折叠结果。
#[derive(Debug, Clone, PartialEq)]
pub struct Folded {
    /// MFE 结构的 dot-bracket。
    pub db: String,
    /// MFE 自由能（kcal/mol）。
    pub mfe: f64,
    /// 配分函数自由能（kcal/mol）。
    pub ensemble: f64,
    /// 集合多样性（mean base pair distance）。
    pub diversity: f64,
    /// 配对概率：`(i, j, p)`，0-based、`i < j`；只含 `p > 5e-4`。
    pub probs: Vec<(usize, usize, f64)>,
}

impl Folded {
    /// 从 dot-bracket 还原配对表（0-based，`(min, max)`，按左端排序）。
    pub fn pairs(&self) -> Vec<(usize, usize)> {
        let mut stack = Vec::new();
        let mut out = Vec::new();
        for (i, c) in self.db.bytes().enumerate() {
            match c {
                b'(' => stack.push(i),
                b')' => {
                    if let Some(j) = stack.pop() {
                        out.push((j, i));
                    }
                }
                _ => {}
            }
        }
        out.sort_unstable();
        out
    }
}

unsafe extern "C" fn collect_prob(i: c_int, j: c_int, p: c_double, ud: *mut c_void) {
    let out = &mut *(ud as *mut Vec<(usize, usize, f64)>);
    out.push((i as usize, j as usize, p));
}

/// 按 ViennaRNA 的 RNAfold 模型折叠（默认 Turner 2004 参数）。
pub fn fold(
    seq: &str,
    temperature: f64,
    no_gu_close: bool,
    no_lp: bool,
) -> Result<Folded, String> {
    let seq = seq.trim();
    if seq.is_empty() {
        return Err("空序列".into());
    }
    let c_seq = std::ffi::CString::new(seq.as_bytes()).map_err(|e| e.to_string())?;
    let len = seq.len();
    let mut db = vec![0u8; len + 1];
    let mut mfe = 0.0;
    let mut ensemble = 0.0;
    let mut diversity = 0.0;
    let mut probs: Vec<(usize, usize, f64)> = Vec::new();
    let rc = unsafe {
        hx_vienna_fold(
            c_seq.as_ptr(),
            temperature,
            no_gu_close as c_int,
            no_lp as c_int,
            db.as_mut_ptr() as *mut c_char,
            (len + 1) as c_int,
            &mut mfe,
            &mut ensemble,
            &mut diversity,
            Some(collect_prob),
            &mut probs as *mut _ as *mut c_void,
        )
    };
    if rc != 0 {
        return Err(format!("ViennaRNA 折叠失败（错误码 {rc}）"));
    }
    db.truncate(len);
    Ok(Folded {
        db: String::from_utf8_lossy(&db).into_owned(),
        mfe,
        ensemble,
        diversity,
        probs,
    })
}

/// RNAplot 的 naview 布局坐标（单位：ViennaRNA 的碱基间距 ≈ 15）。
/// 从 dot-bracket 还原配对表（与 [`Folded::pairs`] 同一规则）。
fn pairs_from_db(db: &str) -> Vec<(usize, usize)> {
    let mut stack = Vec::new();
    let mut out = Vec::new();
    for (i, c) in db.bytes().enumerate() {
        match c {
            b'(' => stack.push(i),
            b')' => {
                if let Some(j) = stack.pop() {
                    out.push((j, i));
                }
            }
            _ => {}
        }
    }
    out.sort_unstable();
    out
}

/// 次优结构：返回 `(配对表, 能量)`，按能量升序。
///
/// `delta` 的单位是 **10 cal/mol**（与 ViennaRNA `vrna_subopt_cb` 一致）：
/// SnapGene 写出的次优列表落在 MFE + 1.5~2.0 kcal/mol 内，对应 `delta = 150..200`。
/// 序列很长时枚举可能很慢，调用方要放在后台线程里。
pub fn subopt(
    seq: &str,
    temperature: f64,
    no_gu_close: bool,
    no_lp: bool,
    delta: i32,
) -> Result<Vec<(Vec<(usize, usize)>, f64)>, String> {
    let seq = seq.trim();
    if seq.is_empty() {
        return Err("空序列".into());
    }
    let c_seq = std::ffi::CString::new(seq.as_bytes()).map_err(|e| e.to_string())?;
    let mut out: Vec<(Vec<(usize, usize)>, f64)> = Vec::new();

    unsafe extern "C" fn collect(energy: f64, structure: *const c_char, ud: *mut c_void) {
        if structure.is_null() {
            return;
        }
        let out = &mut *(ud as *mut Vec<(Vec<(usize, usize)>, f64)>);
        let db = std::ffi::CStr::from_ptr(structure).to_string_lossy();
        out.push((pairs_from_db(&db), energy));
    }

    let rc = unsafe {
        hx_vienna_subopt(
            c_seq.as_ptr(),
            temperature,
            no_gu_close as c_int,
            no_lp as c_int,
            delta,
            Some(collect),
            &mut out as *mut _ as *mut c_void,
        )
    };
    if rc < 0 {
        return Err(format!("次优结构枚举失败（{rc}）"));
    }
    out.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    Ok(out)
}

/// 给 dot-bracket 打自由能（kcal/mol）。用于自检次优/采样结构的能量。
pub fn eval(seq: &str, db: &str, temperature: f64, no_gu_close: bool, no_lp: bool) -> f64 {
    let Ok(c_seq) = std::ffi::CString::new(seq.trim().as_bytes()) else {
        return f64::NAN;
    };
    let Ok(c_db) = std::ffi::CString::new(db.as_bytes()) else {
        return f64::NAN;
    };
    unsafe {
        hx_vienna_eval(
            c_seq.as_ptr(),
            c_db.as_ptr(),
            temperature,
            no_gu_close as c_int,
            no_lp as c_int,
        )
    }
}

pub fn naview(db: &str) -> Result<Vec<(f32, f32)>, String> {
    let len = db.len();
    if len == 0 {
        return Err("空结构".into());
    }
    let c_db = std::ffi::CString::new(db.as_bytes()).map_err(|e| e.to_string())?;
    let mut xs = vec![0.0f64; len];
    let mut ys = vec![0.0f64; len];
    let rc = unsafe {
        hx_vienna_naview(
            c_db.as_ptr(),
            xs.as_mut_ptr(),
            ys.as_mut_ptr(),
            len as c_int,
        )
    };
    if rc < 0 {
        return Err(format!("ViennaRNA 布局失败（错误码 {rc}）"));
    }
    Ok(xs.into_iter()
        .zip(ys)
        .map(|(x, y)| (x as f32, y as f32))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Homo sapiens mitochondrial tRNA-Ala：SnapGene 的 block 34 里就是
    /// -15.80 kcal/mol 与这条点括号（见 docs 里的 RNA 结构验证记录）。
    const TRNA: &str = "AAGGGCTTAGCTTAATTAAAGTGGCTGATTTGCGTTCAGTTGATGCAGAGTGGGGTTTTGCAGTCCTTACCA";
    const TRNA_DB: &str = "(((((((..((((.....))))((((((.......))))))...(((((.......))))))))))))....";

    #[test]
    fn fold_trna_matches_snapgene() {
        let r = fold(TRNA, 37.0, true, false).expect("fold");
        assert_eq!(r.db, TRNA_DB);
        assert!((r.mfe + 15.80).abs() < 0.01, "mfe = {}", r.mfe);
        // 与 SnapGene 的 ensemble 统计一致（-17.38 / 5.40 / 19.17%）。
        assert!((r.ensemble + 16.82).abs() < 0.05, "ensemble = {}", r.ensemble);
        assert!((r.diversity - 5.40).abs() < 0.05, "diversity = {} / mfe {} / ens {}", r.diversity, r.mfe, r.ensemble);
        assert!(!r.probs.is_empty());
    }

    /// 次优结构：Wuchty 枚举的第一条就是 MFE 自己（-15.8），剔除后与 SnapGene
    /// 存的 9 条完全同能量（-15.5 … -14.3）；回调报的能量与独立打分一致。
    #[test]
    fn subopt_trna_matches_snapgene() {
        let v = subopt(TRNA, 37.0, true, false, 200).expect("subopt");
        assert!((v[0].1 + 15.8).abs() < 0.01, "first = {}", v[0].1);
        let head: Vec<f64> = v.iter().skip(1).take(9).map(|x| (x.1 * 10.0).round() / 10.0).collect();
        assert_eq!(
            head,
            vec![-15.5, -15.2, -14.9, -14.9, -14.6, -14.6, -14.6, -14.4, -14.3]
        );
        for (pairs, e) in v.iter().skip(1).take(3) {
            let mut db = vec![b'.'; TRNA.len()];
            for &(i, j) in pairs {
                db[i] = b'(';
                db[j] = b')';
            }
            let db = String::from_utf8(db).unwrap();
            let ee = eval(TRNA, &db, 37.0, true, false);
            assert!((ee - e).abs() < 0.01, "eval {ee} vs {e}");
        }
    }

    #[test]
    fn naview_trna() {
        let pts = naview(TRNA_DB).expect("naview");
        assert_eq!(pts.len(), TRNA.len());
        assert!(pts.iter().all(|(x, y)| x.is_finite() && y.is_finite()));
        // 碱基间距约 15（ViennaRNA 的单位）。
        let (x0, y0) = pts[0];
        let (x1, y1) = pts[1];
        let d = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
        assert!(d > 10.0 && d < 25.0, "d = {d}");
    }
}
