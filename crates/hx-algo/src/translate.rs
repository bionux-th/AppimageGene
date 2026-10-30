//! 标准遗传密码翻译（NCBI 表 1），输出 SnapGene 序列视图用的**三字母**缩写。
//!
//! 只服务于显示：不处理可变起始密码子、稀有密码子或线粒体表；歧义/未知密码子
//! 返回 `Xaa`，终止密码子返回 `***`。编码序列（CDS）由调用方按特征链方向拼好。

/// 一个密码子的三字母缩写（`Xaa` = 含歧义碱基或长度不足；`***` = 终止）。
pub fn codon_aa(codon: &[u8]) -> &'static str {
    if codon.len() != 3 {
        return "Xaa";
    }
    let b = |i: usize| codon[i].to_ascii_uppercase();
    let (x, y, z) = (b(0), b(1), b(2));
    // 只接受 ACGT；其余（N/R/Y…）一律 Xaa，避免误导
    if ![x, y, z].iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
        return "Xaa";
    }
    match x {
        // T 起手：Phe/Leu/Ser/Tyr/Cys/Trp/终止
        b'T' => match (y, z) {
            (b'T', b'T') | (b'T', b'C') => "Phe",
            (b'T', _) => "Leu",
            (b'C', _) => "Ser",
            (b'A', b'T') | (b'A', b'C') => "Tyr",
            (b'G', b'T') | (b'G', b'C') => "Cys",
            (b'G', b'G') => "Trp",
            _ => "***", // TAA / TAG / TGA
        },
        // C 起手：Leu/Pro/His/Gln/Arg
        b'C' => match (y, z) {
            (b'T', _) => "Leu",
            (b'C', _) => "Pro",
            (b'A', b'T') | (b'A', b'C') => "His",
            (b'A', _) => "Gln",
            _ => "Arg", // CGN
        },
        // A 起手：Ile/Met/Thr/Asn/Lys/Ser/Arg
        b'A' => match (y, z) {
            (b'T', b'G') => "Met",
            (b'T', _) => "Ile",
            (b'C', _) => "Thr",
            (b'A', b'T') | (b'A', b'C') => "Asn",
            (b'A', _) => "Lys",
            (b'G', b'T') | (b'G', b'C') => "Ser",
            _ => "Arg", // AGA / AGG
        },
        // G 起手：Val/Ala/Asp/Glu/Gly
        b'G' => match (y, z) {
            (b'T', _) => "Val",
            (b'C', _) => "Ala",
            (b'A', b'T') | (b'A', b'C') => "Asp",
            (b'A', _) => "Glu",
            _ => "Gly", // GGN
        },
        _ => "Xaa",
    }
}

/// 单字母氨基酸（SnapGene 的 `Use 1-Letter Amino Acid Codes`）。
///
/// 不能简单取三字母首字母：Phe→F、Lys→K、Gln→Q、Arg→R、Asn→N 等都要显式映射。
pub fn one_letter(three: &'static str) -> &'static str {
    match three {
        "Phe" => "F",
        "Leu" => "L",
        "Ser" => "S",
        "Pro" => "P",
        "Tyr" => "Y",
        "His" => "H",
        "Gln" => "Q",
        "Arg" => "R",
        "Ile" => "I",
        "Met" => "M",
        "Thr" => "T",
        "Asn" => "N",
        "Lys" => "K",
        "Val" => "V",
        "Ala" => "A",
        "Asp" => "D",
        "Glu" => "E",
        "Gly" => "G",
        "Cys" => "C",
        "Trp" => "W",
        "Xaa" => "X",
        _ => "*",
    }
}

/// 单字母 → 三字母（`one_letter` 的反函数；未知字母返回 `Xaa`）。
///
/// SnapGene 的 `Copy 3-Letter Amino Acids` 用这个方向：蛋白文件里选的是一串
/// 单字母残基，复制时展开成三字母。
pub fn three_letter(one: u8) -> &'static str {
    match one.to_ascii_uppercase() {
        b'F' => "Phe",
        b'L' => "Leu",
        b'S' => "Ser",
        b'P' => "Pro",
        b'Y' => "Tyr",
        b'H' => "His",
        b'Q' => "Gln",
        b'R' => "Arg",
        b'I' => "Ile",
        b'M' => "Met",
        b'T' => "Thr",
        b'N' => "Asn",
        b'K' => "Lys",
        b'V' => "Val",
        b'A' => "Ala",
        b'D' => "Asp",
        b'E' => "Glu",
        b'G' => "Gly",
        b'C' => "Cys",
        b'W' => "Trp",
        b'*' => "***",
        _ => "Xaa",
    }
}

/// 按 3 碱基步长翻译；尾部不足一个密码子的碱基忽略。
pub fn translate(cds: &[u8]) -> Vec<&'static str> {
    cds.chunks_exact(3).map(codon_aa).collect()
}

/// 翻译并跳过终止密码子之前的序列（`include_stop` = 是否保留末尾的 `***`）。
pub fn translate_orf(cds: &[u8], include_stop: bool) -> Vec<&'static str> {
    let mut out = translate(cds);
    if !include_stop {
        if let Some(p) = out.iter().position(|a| *a == "***") {
            out.truncate(p);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codons_map_to_three_letter_codes() {
        assert_eq!(codon_aa(b"ATG"), "Met");
        assert_eq!(codon_aa(b"atg"), "Met");
        assert_eq!(codon_aa(b"TTT"), "Phe");
        assert_eq!(codon_aa(b"GCT"), "Ala");
        assert_eq!(codon_aa(b"TAA"), "***");
        assert_eq!(codon_aa(b"TAG"), "***");
        assert_eq!(codon_aa(b"TGA"), "***");
    }

    #[test]
    fn one_letter_matches_three_letter() {
        assert_eq!(one_letter(codon_aa(b"ATG")), "M");
        assert_eq!(one_letter(codon_aa(b"TTT")), "F");
        assert_eq!(one_letter(codon_aa(b"AAA")), "K");
        assert_eq!(one_letter(codon_aa(b"CAA")), "Q");
        assert_eq!(one_letter(codon_aa(b"TAA")), "*");
        assert_eq!(one_letter("Xaa"), "X");
        // 64 个密码子：1 字母与 3 字母都能对上参考表
        for (three, one) in [
            ("Phe", "F"), ("Leu", "L"), ("Ile", "I"), ("Met", "M"), ("Val", "V"),
            ("Ser", "S"), ("Pro", "P"), ("Thr", "T"), ("Ala", "A"), ("Tyr", "Y"),
            ("His", "H"), ("Gln", "Q"), ("Asn", "N"), ("Lys", "K"), ("Asp", "D"),
            ("Glu", "E"), ("Cys", "C"), ("Trp", "W"), ("Arg", "R"), ("Gly", "G"),
        ] {
            assert_eq!(one_letter(three), one, "{three}");
        }
    }

    #[test]
    fn all_sixty_four_codons_map_correctly() {
        // 参考表（NCBI 表 1，三字母）——整表比对，避免漏group（曾漏 CTN/AGC）
        const REF: &str = "\
AAA Lys AAC Asn AAG Lys AAT Asn \
ACA Thr ACC Thr ACG Thr ACT Thr \
AGA Arg AGC Ser AGG Arg AGT Ser \
ATA Ile ATC Ile ATG Met ATT Ile \
CAA Gln CAC His CAG Gln CAT His \
CCA Pro CCC Pro CCG Pro CCT Pro \
CGA Arg CGC Arg CGG Arg CGT Arg \
CTA Leu CTC Leu CTG Leu CTT Leu \
GAA Glu GAC Asp GAG Glu GAT Asp \
GCA Ala GCC Ala GCG Ala GCT Ala \
GGA Gly GGC Gly GGG Gly GGT Gly \
GTA Val GTC Val GTG Val GTT Val \
TAA *** TAC Tyr TAG *** TAT Tyr \
TCA Ser TCC Ser TCG Ser TCT Ser \
TGA *** TGC Cys TGG Trp TGT Cys \
TTA Leu TTC Phe TTG Leu TTT Phe";
        let toks: Vec<&str> = REF.split_whitespace().collect();
        assert_eq!(toks.len(), 128, "参考表应为 64 组");
        for pair in toks.chunks(2) {
            let (codon, want) = (pair[0], pair[1]);
            assert_eq!(codon_aa(codon.as_bytes()), want, "密码子 {codon}");
        }
    }

    #[test]
    fn ambiguous_and_short_codons_are_xaa() {
        assert_eq!(codon_aa(b"NNN"), "Xaa");
        assert_eq!(codon_aa(b"GCN"), "Xaa");
        assert_eq!(codon_aa(b"AT"), "Xaa");
        assert!(translate(b"AT").is_empty());
    }

    #[test]
    fn translates_frame_and_orf() {
        // KanR 起始（pET28b 特征 560-1375 的反向互补 5' 端）示例片段
        let aa = translate(b"ATGGCTAAA");
        assert_eq!(aa, vec!["Met", "Ala", "Lys"]);
        let orf = translate_orf(b"ATGGCCTAA GTG", false);
        assert_eq!(orf, vec!["Met", "Ala"]);
        let with_stop = translate_orf(b"ATGTAA", true);
        assert_eq!(with_stop, vec!["Met", "***"]);
    }
}
