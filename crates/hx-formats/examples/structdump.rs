//! 临时工具：打印某个 `.dna`/`.rna` 文件 block 34 里的二级结构（dot-bracket）。
use hx_formats::dna::DnaFile;
use hx_formats::rna_struct::RnaStructure;

fn main() {
    let path = std::env::args().nth(1).expect("usage: structdump <file>");
    let bytes = std::fs::read(&path).unwrap();
    let file = DnaFile::parse(&bytes).unwrap();
    let seq = file.sequence().map(|r| r.sequence.clone()).unwrap_or_default();
    let block = file.block(34).expect("no block 34");
    let st = RnaStructure::parse_payload(&block.payload).unwrap();
    let all = st.structures();
    println!("sequence len {} / {} structures", seq.len(), all.len());
    println!("{}", String::from_utf8_lossy(&seq));
    for (idx, (pairs, e)) in all.iter().enumerate() {
        println!("\n#{}  E={e:.2} kcal/mol  pairs={}", idx, pairs.len());
        println!("{}", st.dot_bracket(pairs, seq.len()));
    }
    println!("\nstats {:?}", st.stats);

    println!("MFE pairs: {:?}", all[0].0);
}
