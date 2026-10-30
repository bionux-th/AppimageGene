fn main() {
    use hx_enzymes::sets::{Parent,Filter};
    for parent in [Parent::All,Parent::Nonredundant,Parent::Neb,Parent::Thermo] {
        println!("{parent:?}: {}",hx_enzymes::sets::parent_indices(parent).len());
        for filter in [Filter::SixPlus,Filter::TypeIis,Filter::GoldenGate,Filter::Nicking] {
            println!("  {filter:?}: {}", hx_enzymes::sets::select(parent,filter,b"",false).len());
        }
    }
    if let Some(path) = std::env::args().nth(1) {
        let b = std::fs::read(path).unwrap(); let mut at = 0;
        while at + 5 <= b.len() {
            let kind = b[at]; let n = u32::from_be_bytes(b[at+1..at+5].try_into().unwrap()) as usize;
            if kind == 0 {
                let seq = &b[at+6..at+5+n]; let circular = b[at+5]&1 != 0;
                for f in [Filter::Unique, Filter::UniqueDual, Filter::UniqueSixPlus] {
                    println!("sequence {f:?}: {}",hx_enzymes::sets::select(Parent::Nonredundant,f,seq,circular).len());
                }
            }
            at += 5 + n;
        }
    }
}
