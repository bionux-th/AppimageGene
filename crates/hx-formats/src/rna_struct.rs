//! SnapGene block 34：RNA 二级结构（XZ 压缩 JSON）。
//!
//! 实测（8.2.3，两个官方 `.rna` 样例）的顶层键：
//!
//! ```json
//! {"bondedPairs":[[i,j],...],
//!  "parameters":{"noClosingGU":true,"noLP":false,"temperature":37},
//!  "probabilityMatrix":[[i,j,p],...],   // 稀疏上三角，p = 配对概率
//!  "revision":1,
//!  "stats":{"ensembleDiversity":…,"freeEnergyEnsemble":…,
//!           "freeEnergyOptimal":…,"freqOptimalEnsemble":…},
//!  "suboptimalStructures":[{"bondedPairs":[…],"freeEnergy":…},…]}
//! ```
//!
//! SnapGene 的 `Secondary Structure` 视图**完全由这份数据驱动**（能量、置信度
//! 着色、次优结构缩略图、底部统计），所以读它会得到与 SnapGene 一致的画面；
//! 文件没有这块时才需要本地计算。

use crate::error::FormatError;
use crate::json::{self, Json};

/// 一条次优结构。
#[derive(Debug, Clone, PartialEq)]
pub struct Suboptimal {
    pub pairs: Vec<(usize, usize)>,
    pub energy: f64,
}

/// 底部统计栏的五个数字（后四个来自配分函数）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StructStats {
    /// 系综自由能（`Free energy of ensemble`，kcal/mol）。
    pub free_energy_ensemble: f64,
    /// 最优（MFE）自由能（`Free energy of optimal`）。
    pub free_energy_optimal: f64,
    /// 最优结构在系综中的比例（0..1，界面显示百分比）。
    pub freq_optimal: f64,
    /// 系综多样性（Σ p(1-p)）。
    pub ensemble_diversity: f64,
}

/// 一个 `.rna` 文件里存的二级结构数据。
#[derive(Debug, Clone, PartialEq)]
pub struct RnaStructure {
    /// 最小自由能结构的配对（0-based，`(min,max)`）。
    pub pairs: Vec<(usize, usize)>,
    /// 稀疏配对概率矩阵：`(i, j, p)`，`i < j`。
    pub probs: Vec<(usize, usize, f64)>,
    /// 次优结构（不含 MFE 本身）。
    pub suboptimal: Vec<Suboptimal>,
    pub stats: Option<StructStats>,
    pub temperature: f64,
    pub no_closing_gu: bool,
    pub no_lonely_pairs: bool,
    /// 按碱基整理的配对概率（`probs` 的倒排索引，解析时构建）。
    by_base: Vec<Vec<(usize, f64)>>,
}

impl RnaStructure {
    /// 从 block 34 的载荷解析：先解 XZ，再读 JSON。
    pub fn parse_payload(payload: &[u8]) -> Result<Self, FormatError> {
        let raw = crate::xz::decompress(payload)
            .map_err(|e| FormatError::BadStructure(format!("解压 block 34：{e}")))?;
        Self::parse_json(&raw)
    }

    /// 从结构 JSON 文本解析。
    pub fn parse_json(raw: &[u8]) -> Result<Self, FormatError> {
        let root = json::parse(raw).map_err(|e| FormatError::BadStructure(e.to_string()))?;
        Self::from_json(&root)
    }

    /// 从已解析的 JSON 树取字段。缺字段时给出保守默认值（SnapGene 写出的块
    /// 这四项都在，但不要让一个缺项毁掉整个视图）。
    pub fn from_json(root: &Json) -> Result<Self, FormatError> {
        let pairs = parse_pairs(root.get("bondedPairs"))
            .ok_or_else(|| FormatError::BadStructure("缺少 bondedPairs".into()))?;
        let probs = parse_probs(root.get("probabilityMatrix")).unwrap_or_default();
        let suboptimal = root
            .get("suboptimalStructures")
            .and_then(Json::as_array)
            .map(|items| {
                items
                    .iter()
                    .filter_map(|it| {
                        Some(Suboptimal {
                            pairs: parse_pairs(it.get("bondedPairs"))?,
                            energy: it.get("freeEnergy").and_then(Json::as_f64)?,
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();
        let stats = root.get("stats").and_then(|s| {
            Some(StructStats {
                free_energy_ensemble: s.get("freeEnergyEnsemble")?.as_f64()?,
                free_energy_optimal: s.get("freeEnergyOptimal")?.as_f64()?,
                freq_optimal: s.get("freqOptimalEnsemble")?.as_f64()?,
                ensemble_diversity: s.get("ensembleDiversity")?.as_f64()?,
            })
        });
        let params = root.get("parameters");
        let temperature = params
            .and_then(|p| p.get("temperature"))
            .and_then(Json::as_f64)
            .unwrap_or(37.0);
        let no_closing_gu = params
            .and_then(|p| p.get("noClosingGU"))
            .and_then(Json::as_bool)
            .unwrap_or(true);
        let no_lonely_pairs = params
            .and_then(|p| p.get("noLP"))
            .and_then(Json::as_bool)
            .unwrap_or(false);
        let mut s = Self {
            pairs,
            probs,
            suboptimal,
            stats,
            temperature,
            no_closing_gu,
            no_lonely_pairs,
            by_base: Vec::new(),
        };
        s.reindex();
        Ok(s)
    }

    /// 由外部折叠器（界面层的内嵌 ViennaRNA）的结果构造。
    ///
    /// 与读 block 34 得到的结构完全同构：后面所有绘制/导出代码不用分叉。
    pub fn from_fold(
        pairs: Vec<(usize, usize)>,
        probs: Vec<(usize, usize, f64)>,
        suboptimal: Vec<Suboptimal>,
        stats: StructStats,
        temperature: f64,
        no_closing_gu: bool,
        no_lonely_pairs: bool,
    ) -> Self {
        let mut s = Self {
            pairs,
            probs,
            suboptimal,
            stats: Some(stats),
            temperature,
            no_closing_gu,
            no_lonely_pairs,
            by_base: Vec::new(),
        };
        s.reindex();
        s
    }

    fn reindex(&mut self) {
        let n = self
            .probs
            .iter()
            .flat_map(|(i, j, _)| [*i, *j])
            .max()
            .map(|m| m + 1)
            .unwrap_or(0)
            .max(self.pairs.iter().flat_map(|(i, j)| [*i, *j]).max().map(|m| m + 1).unwrap_or(0));
        self.by_base = vec![Vec::new(); n];
        for &(i, j, p) in &self.probs {
            self.by_base[i].push((j, p));
            self.by_base[j].push((i, p));
        }
        for v in self.by_base.iter_mut() {
            v.sort_by(|a, b| b.1.total_cmp(&a.1));
        }
    }

    /// 某碱基最可能的配对对象与概率。
    pub fn best_pair(&self, i: usize) -> Option<(usize, f64)> {
        self.by_base.get(i).and_then(|v| v.first().copied())
    }

    /// 某碱基的配对概率总和（系综里“这个碱基参与配对”的概率）。
    pub fn paired_prob(&self, i: usize) -> f64 {
        self.by_base
            .get(i)
            .map(|v| v.iter().map(|(_, p)| p).sum::<f64>())
            .unwrap_or(0.0)
            .min(1.0)
    }

    /// 某碱基“保持未配对”的概率。
    pub fn unpaired_prob(&self, i: usize) -> f64 {
        1.0 - self.paired_prob(i)
    }

    /// 某对碱基的配对概率（不在矩阵里时为 0）。
    pub fn pair_prob(&self, i: usize, j: usize) -> f64 {
        self.by_base
            .get(i)
            .and_then(|v| v.iter().find(|(k, _)| *k == j).map(|(_, p)| *p))
            .unwrap_or(0.0)
    }

    /// dot-bracket 表示（`Export as… ▸ Dot-Bracket…` 的文本内容）。
    pub fn dot_bracket(&self, pairs: &[(usize, usize)], len: usize) -> String {
        let mut out = vec![b'.'; len];
        for &(i, j) in pairs {
            if i < len && j < len && i < j {
                out[i] = b'(';
                out[j] = b')';
            }
        }
        String::from_utf8(out).unwrap_or_default()
    }

    /// 序列化成 SnapGene 写出的那种紧凑 JSON（键按字母序，与实测一致）。
    ///
    /// 用途：把本地折叠的结果写回 block 34，让 SnapGene 也能显示这份结构
    /// （SnapGene 的 `Secondary Structure` 页只显示块里存的结构，块不在就一直
    /// 转圈——见 2026-10-21 的实测）。
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(64 + self.pairs.len() * 8 + self.probs.len() * 26);
        s.push_str("{\"bondedPairs\":");
        push_pairs(&mut s, &self.pairs);
        s.push_str(",\"parameters\":{\"noClosingGU\":");
        s.push_str(if self.no_closing_gu { "true" } else { "false" });
        s.push_str(",\"noLP\":");
        s.push_str(if self.no_lonely_pairs { "true" } else { "false" });
        s.push_str(",\"temperature\":");
        s.push_str(&num(self.temperature));
        s.push_str("},\"probabilityMatrix\":[");
        for (k, (i, j, p)) in self.probs.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            s.push_str(&format!("[{i},{j},{}]", num(*p)));
        }
        s.push_str("],\"revision\":1");
        if let Some(st) = self.stats {
            s.push_str(",\"stats\":{\"ensembleDiversity\":");
            s.push_str(&num(st.ensemble_diversity));
            s.push_str(",\"freeEnergyEnsemble\":");
            s.push_str(&num(st.free_energy_ensemble));
            s.push_str(",\"freeEnergyOptimal\":");
            s.push_str(&num(st.free_energy_optimal));
            s.push_str(",\"freqOptimalEnsemble\":");
            s.push_str(&num(st.freq_optimal));
            s.push('}');
        }
        s.push_str(",\"suboptimalStructures\":[");
        for (k, sub) in self.suboptimal.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            s.push_str("{\"bondedPairs\":");
            push_pairs(&mut s, &sub.pairs);
            s.push_str(",\"freeEnergy\":");
            s.push_str(&num(sub.energy));
            s.push('}');
        }
        s.push_str("]}");
        s
    }

    /// block 34 的载荷（XZ 压缩的 JSON）。
    pub fn to_payload(&self) -> Vec<u8> {
        crate::xz::compress_store(self.to_json().as_bytes())
    }

    /// 界面上的结构列表：MFE 在前，其后是块里存的次优结构。
    pub fn structures(&self) -> Vec<(Vec<(usize, usize)>, f64)> {
        let mfe_energy = self.stats.map(|s| s.free_energy_optimal).unwrap_or(0.0);
        let mut v = vec![(self.pairs.clone(), mfe_energy)];
        v.extend(self.suboptimal.iter().map(|s| (s.pairs.clone(), s.energy)));
        v
    }
}

/// 紧凑写配对表：`[[i,j],[i,j],…]`。
fn push_pairs(s: &mut String, pairs: &[(usize, usize)]) {
    s.push('[');
    for (k, (i, j)) in pairs.iter().enumerate() {
        if k > 0 {
            s.push(',');
        }
        s.push_str(&format!("[{i},{j}]"));
    }
    s.push(']');
}

/// JSON 数字：最短往返写法（`37`、`-16.82`、`0.19169`），非有限值写 0。
fn num(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else {
        "0".to_string()
    }
}

fn parse_pairs(v: Option<&Json>) -> Option<Vec<(usize, usize)>> {
    let items = v?.as_array()?;
    let mut out = Vec::with_capacity(items.len());
    for it in items {
        let a = it.as_array()?;
        if a.len() < 2 {
            continue;
        }
        let i = a[0].as_usize()?;
        let j = a[1].as_usize()?;
        out.push((i.min(j), i.max(j)));
    }
    out.sort_unstable();
    out.dedup();
    Some(out)
}

fn parse_probs(v: Option<&Json>) -> Option<Vec<(usize, usize, f64)>> {
    let items = v?.as_array()?;
    let mut out = Vec::with_capacity(items.len());
    for it in items {
        let a = it.as_array()?;
        if a.len() < 3 {
            continue;
        }
        let i = a[0].as_usize()?;
        let j = a[1].as_usize()?;
        let p = a[2].as_f64()?;
        if p > 0.0 {
            out.push((i.min(j), i.max(j), p));
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> Vec<u8> {
        let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name);
        std::fs::read(&p).unwrap()
    }

    #[test]
    fn parses_snapgene_trna_block() {
        let s = RnaStructure::parse_payload(&fixture("sg_trna_block34.xz")).unwrap();
        assert_eq!(s.pairs.len(), 22);
        assert_eq!(s.pairs[0], (0, 67));
        assert_eq!(s.probs.len(), 321);
        assert_eq!(s.suboptimal.len(), 9);
        let st = s.stats.unwrap();
        assert!((st.free_energy_optimal + 15.8).abs() < 1e-9);
        assert!((st.free_energy_ensemble + 16.82).abs() < 1e-9);
        assert!((st.freq_optimal - 0.1917).abs() < 0.001);
        assert!((st.ensemble_diversity - 5.40).abs() < 0.001);
        assert_eq!(s.temperature, 37.0);
        assert!(s.no_closing_gu);
        assert!(!s.no_lonely_pairs);
        // 概率检索
        assert!((s.pair_prob(0, 67) - 0.885506394).abs() < 1e-6);
        assert_eq!(s.best_pair(0).unwrap().0, 67);
        assert!(s.paired_prob(0) > 0.88);
        assert!(s.unpaired_prob(71) > 0.99);
    }

    #[test]
    fn dot_bracket_shape() {
        let s = RnaStructure::parse_payload(&fixture("sg_trna_block34.xz")).unwrap();
        let db = s.dot_bracket(&s.pairs, 72);
        assert_eq!(db.len(), 72);
        assert_eq!(&db[0..1], "(");
        assert_eq!(&db[67..68], ")");
        assert_eq!(db.matches('(').count(), 22);
        assert_eq!(db.matches(')').count(), 22);
    }

    #[test]
    fn payload_roundtrip_keeps_every_field() {
        // 写 block 34（本地折叠结果）之后，SnapGene 与我们读回的必须是同一份数据。
        let s = RnaStructure::parse_payload(&fixture("sg_trna_block34.xz")).unwrap();
        let again = RnaStructure::parse_payload(&s.to_payload()).unwrap();
        assert_eq!(again.pairs, s.pairs);
        assert_eq!(again.probs, s.probs);
        assert_eq!(again.suboptimal, s.suboptimal);
        assert_eq!(again.stats, s.stats);
        assert_eq!(again.temperature, s.temperature);
        assert_eq!(again.no_closing_gu, s.no_closing_gu);
        assert_eq!(again.no_lonely_pairs, s.no_lonely_pairs);
        // JSON 里 `probabilityMatrix` 的顺序也要保持
        assert_eq!(again.to_json(), s.to_json());
    }

    #[test]
    fn structures_list_puts_mfe_first() {
        let s = RnaStructure::parse_payload(&fixture("sg_trna_block34.xz")).unwrap();
        let list = s.structures();
        assert_eq!(list.len(), 10);
        assert!((list[0].1 + 15.8).abs() < 1e-9);
        assert!((list[1].1 + 15.5).abs() < 1e-9);
    }
}
