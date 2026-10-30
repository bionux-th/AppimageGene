# SnapGene 8.2.3 cloning oracle — 2026-10-25

These FASTA records are extracted from products **generated and saved by the
installed SnapGene GUI**, not computed from HelixBench. All inputs are synthetic.
No license, user profile, credentials or private biological data are retained here.

Target: `/opt/gslbiotech/snapgene/snapgene`, SHA-256
`ba77b869c7f2648f2deb324eea371430a980da8f7212745a2f382fa291b9f02a`.
GUI: private HOME and X display, 2560×1440 desktop / 1800×1100 action window.
Temporary raw output and screenshots:
`/tmp/opencode/sg-cloning-oracle-aanuouze/` (not a durable dependency).

## Restriction cloning

- Vector: circular, `A × 100 + GAATTC + C × 30 + GGATCC + ACCG × 170`, 822 bp.
- Insert: linear, `T × 20 + GAATTC + ACCA × 40 + GGATCC + T × 20`, 212 bp.
- Both sources initially have blunt unmodified ends and no annotations.
- Actions → Restriction and Insertion Cloning → Insert Fragment.
- Vector: cut EcoRI and BamHI, default region **to replace** EcoRI(101)–BamHI(137),
  36 bp, leaving 786 bp backbone. Insert: EcoRI(21)–BamHI(187), 166 bp.
- No blunting, both orientations forward. SG displays Ready to clone / 952 bp.
- Clone, Save As `restriction-product.dna`.
- Sequence extracted from native block 0, without flags; circular 952 bp.
- Sequence SHA-256:
  `201837d81e72459092b5113c88b512abf6f65851d4b31434d64b3546d3b1bf84`.
- SG retains the vector's original position 1: the product starts with 100 A,
  **not** at the BamHI opening. Tests compare literal sequence, without rotation.
- Native block 8: stickiness 0/0; both modifications `FivePrimePhosphorylated`.
- Screenshot: `restriction-ready.png`; raw product: `restriction-product.dna`.

### Both 5′ ends filled in

Repeat the restriction case with the EcoRI and BamHI **fill in 5′ overhangs**
checkboxes selected for both sources. SG shows backbone 790 bp, insert 170 bp,
product 960 bp (and one warning, not yet interpreted).
The saved product `restriction-blunt-product.dna` is exactly the sticky-product
sequence with AATT inserted at boundary 105 and GATC at boundary 271 (boundaries
in the original 952 bp product, applied simultaneously). Its sequence SHA-256 is
`052079490ed6242c5f86ceceac0cdda27b296fdc19e23e16e44f8984b45bc928`.
This is a third saved SG product; it validates the existing fill-in boundary rule
and the corrected origin numbering, not every per-enzyme modification setting.

## Linear ligation

- Left: stored `T + G × 79 + AATT`, 84 bases; ends -1/+4, both Unmodified.
- Right: stored `AATT + C × 80 + AA`, 86 bases; ends +4/+2, both Unmodified.
- Actions → Linear Ligation; select left as Fragment 1, right as Fragment 2.
- SG rejects despite complementary AATT ends, displaying:
  `For ligation, at least one DNA end must be 5' phosphorylated at the junction between Fragments 1 and 2.`
- Change **only** left's DownstreamModification to `FivePrimePhosphorylated`
  in a new synthetic input `left-phos.dna`. Keep right unchanged.
- SG now displays Ready to ligate / 163 bp. Ligate, Save As `ligation-product.dna`.
- Native stored sequence 166 bases; top strand 163 bp; output is linear.
- Stored sequence SHA-256:
  `88295e8bfd0d133bfc4949152d1734aba35df12b8d4a57a83fdd5f1c068ff50b`.
- Native output ends -1/+2, both outer modifications Unmodified; SG additionally
  writes UpstreamEnzymeName=Start and DownstreamEnzymeName=End.
- Screenshots: `ligation-ready.png` (rejection), `ligation-phos-ready.png` (success).

## Coverage boundaries

This verifies three controlled products and one phosphorylation rejection.
It does not establish parity for other chemical modifications, partial digestion,
all default selections, IUPAC junction resolution, origin deletion/reversed-vector
numbering, annotated-product merging, intermediate files or complete native-file
metadata/history. No Ghidra inference is used for these observations.
