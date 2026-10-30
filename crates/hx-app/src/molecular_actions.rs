//! Independent, non-destructive molecular products from explicit user parameters.
use eframe::egui;
use hx_core::{MoleculeFlags, MoleculeKind};
use hx_formats::{Block, DnaFile};
use crate::doc::Document;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Action { Linearize, SingleStrand, Anneal, Rna, Protein }

impl Action {
    pub fn title(self) -> &'static str { match self {
        Self::Linearize => "Linearize", Self::SingleStrand => "Convert to Single-Stranded",
        Self::Anneal => "Anneal Oligos", Self::Rna => "Make RNA", Self::Protein => "Make Protein",
    } }
    pub fn allowed(self, doc: Option<&Document>) -> bool {
        if self == Self::Anneal { return true; }
        doc.is_some_and(|d| d.native_class().is_none() && d.len() > 0 && match self {
            Self::Linearize => d.kind() == MoleculeKind::Dna && d.circular(),
            Self::SingleStrand => d.kind().is_nucleotide() && d.session.flags().has(MoleculeFlags::DOUBLE_STRANDED),
            Self::Rna => d.kind() == MoleculeKind::Dna,
            Self::Protein => d.kind().is_nucleotide(), Self::Anneal => true,
        })
    }
}

pub struct Dialog {
    pub action: Action,
    source: Option<DnaFile>,
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub frame: usize,
    pub stop: bool,
    pub oligo1: String,
    pub oligo2: String,
    pub offset: i32,
    pub name: String,
    error: String,
}

impl Dialog {
    pub fn new(action: Action, doc: Option<&Document>, range: Option<(usize, usize)>) -> Result<Self, String> {
        if !action.allowed(doc) { return Err("This action is not applicable to the active molecule.".into()); }
        let source = if action == Action::Anneal { None } else {
            let d = doc.unwrap();
            if d.session.end_properties().upstream_stickiness != 0 || d.session.end_properties().downstream_stickiness != 0 {
                return Err("This action currently requires a source without sticky ends. Prepare a blunt source first.".into());
            }
            let (bytes, _) = d.session.save_bytes_with_history(&d.map_name(), &d.map_name(), None)?;
            Some(DnaFile::parse(&bytes).map_err(|e| e.to_string())?)
        };
        let len = doc.map(|d| d.len()).unwrap_or(0);
        let (start, end) = range.filter(|(s, e)| s < e).unwrap_or((0, len));
        Ok(Self { action, source, start: start + 1, end, reverse: false, frame: 1, stop: true,
            oligo1: String::new(), oligo2: String::new(), offset: 0,
            name: format!("{} product", action.title()), error: String::new() })
    }

    pub fn build(&self) -> Result<(Vec<u8>, String), String> {
        if self.name.trim().is_empty() { return Err("Please enter a product name.".into()); }
        let (mut file, operation, summary) = if self.action == Action::Anneal {
            (anneal(&self.oligo1, &self.oligo2, self.offset)?, "annealOligos", hx_formats::history::InputSummary::empty())
        } else {
            let src = self.source.as_ref().ok_or("Missing source")?;
            let rec = src.sequence().ok_or("Missing source sequence")?;
            let len = rec.len();
            let source = Document::from_bytes(src.to_bytes())?.cloning_source("source");
            let (seq, features, kind, flags, op) = match self.action {
                Action::Linearize => {
                    if self.start == 0 || self.start > len { return Err("First base must be within 1 .. sequence length.".into()); }
                    let part = source.slice(self.start - 1, self.start - 1, false);
                    (part.seq, part.features, rec.kind, rec.flags.with(MoleculeFlags::CIRCULAR, false), "linearize")
                }
                Action::SingleStrand => {
                    let part = source.slice(0, len, self.reverse);
                    (part.seq, part.features, rec.kind, rec.flags.with(MoleculeFlags::DOUBLE_STRANDED, false), "convertToSingleStranded")
                }
                Action::Rna | Action::Protein => {
                    if self.start == 0 || self.end == 0 || self.start > len || self.end > len || (self.start > self.end && !source.circular) {
                        return Err("Invalid region; wrapping regions require a circular source.".into());
                    }
                    let mut part = source.slice(self.start - 1, self.end, self.reverse);
                    if self.action == Action::Rna {
                        // Native RNA uses T; the RNA view displays U.
                        for b in &mut part.seq { if *b == b'U' { *b = b'T'; } if *b == b'u' { *b = b't'; } }
                        (part.seq, part.features, MoleculeKind::Rna, MoleculeFlags::default(), "makeRNA")
                    } else {
                        if !(1..=3).contains(&self.frame) { return Err("Reading frame must be 1, 2 or 3.".into()); }
                        let seq = part.seq.get(self.frame - 1..).ok_or("Region is shorter than the reading frame")?;
                        let aas = hx_algo::translate::translate(seq);
                        let mut protein = Vec::new();
                        for aa in aas {
                            let letter = hx_algo::translate::one_letter(aa).as_bytes()[0];
                            if letter == b'*' && self.stop { break; }
                            protein.push(letter);
                        }
                        if protein.is_empty() { return Err("No amino acids in the selected region/frame before the first stop.".into()); }
                        (protein, Vec::new(), MoleculeKind::Protein, MoleculeFlags::default(), "makeProtein")
                    }
                }
                Action::Anneal => unreachable!(),
            };
            let mut new = hx_formats::write::NewFile::of(kind, seq).with_flags(flags).with_name(&self.name);
            // Serialize the full supported feature model rather than reducing it.
            new.name = Some(self.name.clone());
            let mut file = DnaFile::parse(&new.to_bytes()).map_err(|e| e.to_string())?;
            let payload = hx_formats::write::feature_annotations_xml(&features, file.sequence().unwrap().len());
            if let Some(b) = file.blocks.iter_mut().find(|b| b.block_type == 10) { b.payload = payload; }
            else { file.blocks.push(Block::new(10, payload)); }
            (file, op, hx_formats::history::InputSummary::range("select", self.start.saturating_sub(1), self.end))
        };
        if let Some(src) = &self.source {
            let mut tree = hx_formats::history::edit_history(&file, src, &[], &self.name, "Source", None).map_err(|e| e.to_string())?;
            tree.operation = operation.into(); tree.input_summary = Some(summary);
            file.blocks.push(Block::new(7, hx_formats::history::to_block_payload(&tree)));
        } else {
            let make_source = |text: &str| -> Result<DnaFile, String> {
                let seq = text.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()).collect();
                DnaFile::parse(&hx_formats::write::NewFile::dna(seq).with_flags(MoleculeFlags::default()).to_bytes()).map_err(|e| e.to_string())
            };
            let first = make_source(&self.oligo1)?;
            let second = make_source(&self.oligo2)?;
            let mut tree = hx_formats::history::edit_history(&file, &first, &[], &self.name, "Oligo 1", None).map_err(|e| e.to_string())?;
            let other = hx_formats::history::edit_history(&file, &second, &[], &self.name, "Oligo 2", None).map_err(|e| e.to_string())?;
            let mut parent = other.children.into_iter().next().ok_or("Missing oligo history")?;
            parent.id = 1; tree.id = 2; tree.children.push(parent);
            tree.operation = operation.into(); tree.input_summary = Some(summary);
            file.blocks.push(Block::new(7, hx_formats::history::to_block_payload(&tree)));
        }
        if let Some(notes) = file.blocks.iter_mut().find(|b| b.block_type == 6) {
            notes.payload = hx_formats::notes::set_field(&notes.payload, "CustomMapLabel", &self.name);
        }
        let extension = match file.header.molecule_kind().unwrap() { MoleculeKind::Dna => "dna", MoleculeKind::Rna => "rna", MoleculeKind::Protein => "prot" };
        Ok((file.to_bytes(), format!("{}.{extension}", self.name.trim())))
    }

    pub fn show(&mut self, ctx: &egui::Context) -> Option<Result<(Vec<u8>, String), ()>> {
        let mut outcome = None;
        egui::Modal::new(egui::Id::unique("molecular_action")).show(ctx, |ui| {
            ui.set_width(440.0); ui.heading(crate::i18n::text(ctx,self.action.title())); ui.separator();
            match self.action {
                Action::Linearize => { ui.label(crate::i18n::text(ctx,"Open the circle immediately before this base (blunt cut):")); ui.add(egui::DragValue::new(&mut self.start).range(1..=usize::MAX)); }
                Action::SingleStrand => { ui.checkbox(&mut self.reverse, crate::i18n::text(ctx,"Use bottom strand (reverse complement, 5′ → 3′)")); ui.label(crate::i18n::text(ctx,"Creates a new single-stranded molecule; the source is unchanged.")); }
                Action::Anneal => {
                    ui.label(crate::i18n::text(ctx,"Oligo 1 (5′ → 3′)")); ui.text_edit_multiline(&mut self.oligo1);
                    ui.label(crate::i18n::text(ctx,"Oligo 2 (5′ → 3′)")); ui.text_edit_multiline(&mut self.oligo2);
                    ui.label(crate::i18n::text(ctx,"Offset of reverse-complemented Oligo 2 relative to Oligo 1:"));
                    ui.add(egui::DragValue::new(&mut self.offset));
                    ui.label(crate::i18n::text(ctx,"0 aligns left edges. Positive shifts it right. Only exact A/C/G/T duplexes are accepted; no mismatch or melting-temperature prediction."));
                }
                Action::Rna | Action::Protein => {
                    ui.horizontal(|ui| { ui.label(crate::i18n::text(ctx,"Region (1-based, inclusive):")); ui.add(egui::DragValue::new(&mut self.start).range(1..=usize::MAX)); ui.label(".."); ui.add(egui::DragValue::new(&mut self.end).range(1..=usize::MAX)); });
                    ui.checkbox(&mut self.reverse, crate::i18n::text(ctx,"Reverse-complement source region"));
                    if self.action == Action::Protein {
                        ui.horizontal(|ui| { ui.label(crate::i18n::text(ctx,"Reading frame:")); ui.add(egui::DragValue::new(&mut self.frame).range(1..=3)); });
                        ui.checkbox(&mut self.stop, crate::i18n::text(ctx,"Stop at first stop codon"));
                        ui.label(crate::i18n::text(ctx,"Standard genetic code only. Incomplete terminal codons are ignored. Nucleotide features are not copied into protein coordinates."));
                    } else { ui.label(crate::i18n::text(ctx,"The chosen strand is treated as the coding strand (T → U). No promoter recognition or splicing is predicted.")); }
                }
            }
            ui.horizontal(|ui| { ui.label(crate::i18n::text(ctx,"Product name:")); ui.text_edit_singleline(&mut self.name); });
            if !self.error.is_empty() { ui.colored_label(egui::Color32::RED, crate::i18n::text(ctx,&self.error)); }
            ui.horizontal(|ui| {
                if ui.button(crate::i18n::text(ctx,"Create")).clicked() { match self.build() { Ok(p) => outcome = Some(Ok(p)), Err(e) => self.error = e } }
                if ui.button(crate::i18n::text(ctx,"Cancel")).clicked() { outcome = Some(Err(())); }
            });
            if ui.input(|i| i.key_pressed(egui::Key::Escape)) { outcome = Some(Err(())); }
        });
        outcome
    }
}

fn anneal(a: &str, b: &str, offset: i32) -> Result<DnaFile, String> {
    let parse = |s: &str| -> Result<Vec<u8>, String> {
        let seq: Vec<_> = s.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()).collect();
        if seq.is_empty() || !seq.iter().all(|b| b"ACGT".contains(b)) { return Err("Enter non-empty oligos using A/C/G/T only.".into()); }
        Ok(seq)
    };
    let a = parse(a)?; let b = hx_algo::reverse_complement(&parse(b)?);
    let offset = offset as i64;
    let lo = 0.max(offset); let hi = (a.len() as i64).min(offset + b.len() as i64);
    if hi <= lo { return Err("The oligos do not overlap at this offset.".into()); }
    if (lo..hi).any(|x| a[x as usize] != b[(x - offset) as usize]) { return Err("The overlap is not complementary at this offset.".into()); }
    let start = 0.min(offset); let end = (a.len() as i64).max(offset + b.len() as i64);
    let seq = (start..end).map(|x| if x >= 0 && x < a.len() as i64 { a[x as usize] } else { b[(x - offset) as usize] }).collect();
    let mut file = DnaFile::parse(&hx_formats::write::NewFile::dna(seq).with_flags(MoleculeFlags::from_byte(2)).to_bytes()).map_err(|e| e.to_string())?;
    let props = format!("<AdditionalSequenceProperties><UpstreamStickiness>{}</UpstreamStickiness><DownstreamStickiness>{}</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>", offset, offset + b.len() as i64 - a.len() as i64);
    file.blocks.retain(|b| b.block_type != 8); file.blocks.push(Block::new(8, props.into_bytes()));
    Ok(file)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source() -> Document {
        Document::from_bytes(hx_formats::write::NewFile::dna(b"ATGAAATAG".to_vec()).with_flags(MoleculeFlags::from_byte(3)).to_bytes()).unwrap()
    }
    fn reload(d: &Dialog) -> Document {
        let (bytes, _) = d.build().unwrap();
        let doc = Document::from_bytes(bytes.clone()).unwrap();
        assert_eq!(doc.session.save_bytes().unwrap(), bytes);
        static NEXT:std::sync::atomic::AtomicUsize=std::sync::atomic::AtomicUsize::new(0);
        let dir = std::path::PathBuf::from("/tmp/opencode").join(format!("molecular-save-{}-{:?}-{}", std::process::id(), d.action,NEXT.fetch_add(1,std::sync::atomic::Ordering::Relaxed)));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("product.dna");
        assert!(!path.exists());
        let mut disk = Document::from_bytes(bytes).unwrap();
        disk.save_to(&path).unwrap();
        let reopened = Document::from_path(&path).unwrap();
        assert_eq!(reopened.sequence(), doc.sequence());
        assert_eq!(reopened.kind(), doc.kind());
        std::fs::remove_file(path).unwrap(); std::fs::remove_dir(dir).unwrap();
        doc
    }
    #[test]
    fn all_source_actions_build_and_reload_without_changing_source() {
        let source = source(); let before = source.session.bytes().to_vec();
        let mut d = Dialog::new(Action::Linearize, Some(&source), None).unwrap(); d.start = 4;
        let linear = reload(&d); assert_eq!(linear.sequence(), b"AAATAGATG"); assert!(!linear.circular());
        let mut d = Dialog::new(Action::SingleStrand, Some(&source), None).unwrap(); d.reverse = true;
        let ss = reload(&d); assert_eq!(ss.sequence(), b"CTATTTCAT"); assert!(!ss.session.flags().has(2));
        let d = Dialog::new(Action::Rna, Some(&source), Some((0, 6))).unwrap();
        let rna = reload(&d); assert_eq!(rna.kind(), MoleculeKind::Rna); assert_eq!(rna.sequence(), b"ATGAAA");
        let mut d = Dialog::new(Action::Protein, Some(&source), None).unwrap();
        let protein = reload(&d); assert_eq!(protein.kind(), MoleculeKind::Protein); assert_eq!(protein.sequence(), b"MK");
        d.stop = false; assert_eq!(reload(&d).sequence(), b"MK*");
        assert_eq!(source.session.bytes(), before);
        assert!(protein.session.file().block(7).is_some());
    }
    #[test]
    fn oligo_offsets_preserve_top_and_bottom_overhangs() {
        let mut dialog = Dialog::new(Action::Anneal, None, None).unwrap();
        dialog.oligo1 = "AAACGT".into(); dialog.oligo2 = "TTACGT".into(); dialog.offset = 2;
        let saved = reload(&dialog);
        let xml = hx_formats::history::decode_payload(&saved.session.file().block(7).unwrap().payload).unwrap();
        assert_eq!(hx_formats::history::parse_entries(&xml).unwrap().len(), 3);
        let blunt = anneal("ACGT", "ACGT", 0).unwrap(); assert_eq!(blunt.sequence().unwrap().sequence, b"ACGT");
        let a = anneal("AAACGT", "TTACGT", 2).unwrap();
        assert_eq!(a.sequence().unwrap().sequence, b"AAACGTAA");
        assert_eq!((a.end_properties().upstream_stickiness, a.end_properties().downstream_stickiness), (2, 2));
        let doc = Document::from_bytes(a.to_bytes()).unwrap(); assert_eq!(doc.sequence(), b"AAACGT");
        let b = anneal("ACGT", "ACGTAA", -2).unwrap();
        assert_eq!(b.end_properties().upstream_stickiness, -2);
        assert_eq!(Document::from_bytes(b.to_bytes()).unwrap().sequence(), b"ACGT");
        assert!(anneal("ACGT", "AAAA", 0).is_err());
        assert!(anneal("ACGT", "ACGT", 10).is_err());
        assert!(anneal("AC?GT", "ACGT", 0).is_err());
    }
    #[test]
    fn applicability_and_bad_regions_are_rejected() {
        let src = source();
        let p = reload(&Dialog::new(Action::Protein, Some(&src), None).unwrap());
        assert!(!Action::Rna.allowed(Some(&p))); assert!(!Action::SingleStrand.allowed(Some(&p)));
        assert!(Action::Anneal.allowed(None)); assert!(!Action::Linearize.allowed(None));
        let mut d = Dialog::new(Action::Rna, Some(&src), None).unwrap(); d.start = 0; assert!(d.build().is_err());
        let mut d = Dialog::new(Action::Linearize, Some(&src), None).unwrap(); d.start = 100; assert!(d.build().is_err());
    }
}
