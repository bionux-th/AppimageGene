# Third-party notices

## ViennaRNA

The embedded ViennaRNA 2.7.2 source belongs to its authors and the Institute for
Theoretical Chemistry, University of Vienna. Its license is **not MIT**. Preserve
`vendor/viennarna/COPYING`, `AUTHORS`, source notices and `HELIXBENCH-PATCHES.md`.
The naview source has a separate embedded copyright notice. Review these terms
before redistribution, particularly distribution for a fee.

## Rust dependencies

Cargo.lock records exact package versions. egui/eframe, rfd and other dependencies
retain their own licenses; obtain notices from the corresponding crate packages
when preparing binary distributions. The root MIT license does not replace them.

## Enzyme and ladder data

Generated factual records in `crates/hx-enzymes/src/{catalog,data,details}.rs`
and `crates/hx-app/src/gel_markers.xml` were extracted from a locally installed
SnapGene 8.2.3 catalog. Source annotations retain provenance. The project makes
no claim of ownership or blanket MIT licensing over these third-party records.
Redistribution rights and any database rights require separate assessment by
the distributor; absence of URLs/marketing text is not a legal clearance.

## Optional external components

Poppler is not bundled with this source publication. Windows system fonts,
SnapGene executable/resources, activation data and commercial vector templates
are not included. Optional installed vector resources remain subject to their
owners' terms. SnapGene and other third-party product names are used to identify
compatibility and are not an endorsement or affiliation.

## Screenshots

The Image directory contains maintainer-supplied screenshots for README display,
not redistributed sequence/template files. Application icons are project assets;
maintainers should retain the provenance/permissions for any contributed images.
