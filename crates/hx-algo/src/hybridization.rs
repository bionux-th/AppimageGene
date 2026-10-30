//! Ungapped primer/template binding search. A terminal mismatch is never
//! accepted; the optional single mismatch must have matched neighbours.
use hx_formats::annot::{HybridizationParams, PrimerSite, Range};

#[derive(Debug, Clone, PartialEq)]
pub struct Binding {
    pub site: PrimerSite,
    pub five_prime: bool,
    pub mismatch: Option<usize>,
}

pub fn search(template:&[u8],primer:&[u8],circular:bool,params:HybridizationParams)->Vec<Binding>{
    let n=template.len();let m=primer.len().min(n);
    if n==0||m==0{return Vec::new();}
    let mut out=Vec::new();
    for strand in 0..=1u8 {
        let oriented=if strand==0 {template.to_vec()}else{crate::reverse_complement(template)};
        for five_prime in [false,true] {
            if five_prime&&!params.show_additional_five_prime_matches{continue;}
            let min=if five_prime{params.minimum_five_prime_annealing}else{params.min_continuous_match_len}.max(1);
            if min>m{continue;}
            for anchor in 0..if circular{n}else{n+1} {
                let max=if circular{m}else if five_prime{m.min(n.saturating_sub(anchor))}else{m.min(anchor)};
                let mut matched=0;let mut mismatch=None;
                for k in 0..max {
                    let pi=if five_prime{k}else{primer.len()-1-k};
                    let ti=if five_prime{(anchor+k)%n}else{(anchor+n-1-k)%n};
                    if !crate::hybridization::matches(primer[pi],oriented[ti]) {
                        if !params.allow_mismatch||mismatch.is_some()||k==0||k+1>=max{break;}
                        let next_pi=if five_prime{pi+1}else{pi-1};
                        let next_ti=if five_prime{(ti+1)%n}else{(ti+n-1)%n};
                        if !crate::hybridization::matches(primer[next_pi],oriented[next_ti]){break;}
                        mismatch=Some(pi);
                    }
                    matched=k+1;
                }
                if matched<min{continue;}
                let part=if five_prime{&primer[..matched]}else{&primer[primer.len()-matched..]};
                // Exact-match NN Tm is only a screening estimate for an
                // ungapped mismatch, not a mismatch thermodynamic model.
                let tm=crate::primer::tm(part);
                if tm<params.min_melting_temperature{continue;}
                let start=if five_prime{anchor%n}else{(anchor+n-matched)%n};
                let end=if circular{(start+matched)%n}else{start+matched};
                let (start,end)=if strand==0{(start,end)}else if circular{((n-end)%n,(n-start)%n)}else{(n-end,n-start)};
                out.push(Binding{site:PrimerSite{range:Range::new(start,end),strand,annealed_bases:matched,melting_temperature:Some(tm)},five_prime,mismatch});
            }
        }
    }
    out.sort_by_key(|b|(b.site.strand,b.five_prime,std::cmp::Reverse(b.site.annealed_bases),b.site.range.start));
    let mut maximal:Vec<Binding>=Vec::new();
    for candidate in out{
        let covered=maximal.iter().any(|b|b.site.strand==candidate.site.strand&&b.five_prime==candidate.five_prime&&b.site.annealed_bases>=candidate.site.annealed_bases&&{
            let distance=(candidate.site.range.start+n-b.site.range.start)%n;
            distance+candidate.site.annealed_bases<=b.site.annealed_bases
        });
        if !covered{maximal.push(candidate);}
    }
    let mut out=maximal;
    out.sort_by_key(|b|(b.site.strand,b.site.range.start,b.site.range.end,b.five_prime));
    out.dedup_by(|a,b|a.site.range==b.site.range&&a.site.strand==b.site.strand);
    out
}

fn matches(primer:u8,template:u8)->bool{
    match(hx_enzymes::iupac_bits(primer),hx_enzymes::iupac_bits(template)){
        (Some(a),Some(b))=>a&b==b,_=>false,
    }
}

#[cfg(test)]mod tests{
    use super::*;
    fn params()->HybridizationParams{HybridizationParams{min_continuous_match_len:8,minimum_five_prime_annealing:8,min_melting_temperature:-100.0,show_additional_five_prime_matches:false,allow_mismatch:false}}
    #[test]fn linear_terminal_base_and_isolated_mismatch(){
        assert!(search(b"AACCGGTT",b"AACCGGTT",false,params()).iter().any(|b|b.site.range==Range::new(0,8)&&b.site.strand==0));
        assert!(search(b"AACAGGTT",b"AACCGGTT",false,params()).is_empty());
        let p=HybridizationParams{allow_mismatch:true,..params()};
        assert!(search(b"AACAGGTT",b"AACCGGTT",false,p).iter().any(|b|b.mismatch==Some(3)));
        assert!(!search(b"AACCGGTA",b"AACCGGTT",false,p).iter().any(|b|b.site.strand==0));
        assert!(search(b"AATTGGTT",b"AACCGGTT",false,p).is_empty());
    }
    #[test]fn prefix_reverse_and_circular(){
        let p=HybridizationParams{show_additional_five_prime_matches:true,..params()};
        assert!(search(b"AACCGGTTAAAA",b"AACCGGTTCCCC",false,p).iter().any(|b|b.five_prime&&b.site.strand==0));
        assert!(search(b"GGTTAAAAAACC",b"AACCGGTT",true,params()).iter().any(|b|b.site.range.start>b.site.range.end));
        assert!(search(b"AAGTCCGA",b"TCGGACTT",false,params()).iter().any(|b|b.site.strand==1));
    }
    #[test]fn iupac_primer_and_temperature_threshold(){
        assert!(search(b"AACCGGTT",b"AAYCGGTT",false,params()).iter().any(|b|b.site.strand==0));
        assert!(search(b"AACNGGTT",b"AACCGGTT",false,params()).is_empty());
        let p=HybridizationParams{min_melting_temperature:150.0,..params()};assert!(search(b"AACCGGTT",b"AACCGGTT",false,p).is_empty());
    }
}
