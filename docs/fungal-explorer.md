# Fungal combination explorer

The explorer starts with 20 combinations spanning short ITS products, full ITS,
and longer rDNA spans. Choose the sequencing platform, a class or phylum, and
forward/reverse components; use the matrix to open a pair or compare up to five
product-length distributions. Disable the starter filter to search the full
candidate list. Documented assays and new combinations of published components
have separate evidence labels.

## Length and evidence

The project defaults are Illumina 500 bp, PacBio 2,000 bp, and Nanopore 3,500 bp.
The Illumina slider spans 452–500 bp; both longer-product ceilings are editable.
These settings express the intended experiment, not platform specifications.
Short products can meet more than one ceiling. Colors indicate the shortest
configured tier containing the median; eligibility counts inspect every
observed product, and a callout exposes ranges crossing tiers.

A supported product requires two unique, exact IUPAC-compatible binding sites
on the same raw template, facing inward without overlap. An additional
random-match guard excludes implausible highly degenerate sites. Product length
includes both primers and excludes adapters. Missing, ambiguous, or unsupported
sites remain separate from length-ineligible products. Exact-site support is
not a mismatch-tolerant amplification score or experimental PCR validation.

Counts always use the selected snapshot and lineage subset. The selected panel
is a bounded reference sample, not a prevalence sample or a universal fungal
coverage estimate. Amplicon lengths come from raw sequences, independent of
alignment gaps. MAFFT aligns the complete panel without truncating new records;
the full alignment does not establish homologous spacer columns across phyla.

## Reference acquisition and uncertainty

The initial expanded panel retains FN812768.2 and adds four targeted NCBI
nucleotide records for each of Ascomycota, Basidiomycota, Chytridiomycota,
Mucoromycota, and Glomeromycota. The query requires an 18S/small-subunit title,
an internal-transcribed-spacer title, and a 2,500–12,000 bp sequence range.
The retained baseline is shorter. The 21 templates span 739–10,646 bp and seven
resolved classes; unresolved taxonomy remains explicit. The exact queries,
returned accession-version sets, raw GenBank/taxonomy XML hashes, acquisition
time proxy, build time, and MAFFT command are frozen with the snapshot.

Only complete, uniquely annotated gene boundaries become landmarks. Fuzzy
endpoints and duplicated features are withheld, and spacer intervals require
their adjacent boundaries. The expanded panel has nine trustworthy individual
landmarks; every template lacks a fully verified five-region architecture.
Annotation notes accompany the sequence evidence. Primer-site product lengths
remain directly measurable even when a gene boundary is unresolved.

## Versions and recovery

The combination explorer also renders the selected snapshot's actual frozen
multiple-sequence alignment. A shared alignment-column axis shows the current
combination tracks on a chosen template and the selected templates with their
gaps and comparison products. Clicking a combination selects its evidence and
highlights its alignment span; the view can show the shortlist, up to five
compared combinations, or the selected pair. Lineage and version filters also
change the alignment rows and available reference templates.

Raw exact-site coordinates are mapped through the non-gap bases of each frozen
aligned row. Ungapped aligned sequence identity must match its frozen template
before a span can be shown. Missing, ambiguous, or inconsistent sites remain
explicitly unresolved. Product lengths and platform colors use raw bp, never
the number of alignment columns. The zoomable nucleotide window displays actual
bases and gaps with primer-site shading; a CSV exports both coordinate systems
and the snapshot identity. Cross-clade spacer alignment remains exploratory
and does not establish homology or PCR validation.

**Template bases (bp)** is the default axis for the upper combination chart.
Its chosen template's product, primer and gene spans use inclusive raw bases,
so product width corresponds to physical sequence length. **Gapped alignment
columns** restores the shared-column upper view. Cross-template comparisons
and nucleotide windows always retain their explicitly labeled alignment axis;
the alignment slider selects corresponding raw bases for the upper bp chart.

**Where the bp count comes from** shows the selected pair's inclusive raw
coordinate calculation, binding-site/interior lengths and alignment gap count.
For example, PZ962702.1 ITS1 + ITS2 is raw bases 1,775–2,109: 2,109 − 1,775 + 1
= 335 bp, including both primer sites. Its frozen alignment span is columns
7,995–9,567: 1,573 columns = 335 sequence bases + 1,238 inserted gap columns.
Gaps do not add bp to the PCR product or its platform length limit. The bp
chart also exposes each plotted pair's raw and alignment accounting.

**Show genes and regions** adds per-template rRNA gene and ITS spacer tracks
on the same alignment columns as the primer products. Slate identifies genes;
amber identifies spacers independently of sequencing-platform colors. The
18S → ITS1 → 5.8S → ITS2 → 28S guide shows which boundaries are available
on the selected combination-track template. Click an available region to zoom
the overview and move the nucleotide window to that region. The template
selector reports annotation counts and initially prefers an annotated
template with actual product support when one exists.
The overview repeats visible region names and their full alignment bounds
above the chart, marks regions only partly in the zoomed window, and gives a
sideways-scroll cue in narrow layouts.

Every boundary is projected through its own frozen template's non-gap map.
The nucleotide window adds labeled annotation rows and base tooltips; hover a
track for raw bp, alignment columns and its source note. Unavailable boundaries
remain disabled in the guide. V1's 18S/28S segments are labeled partial anchors
and its curated pilot boundaries use dashed borders. V2 displays only its nine
unambiguous accession-specific NCBI annotations, without assigning one
template's gene boundaries to another. **Gene and region coordinates (CSV)**
exports the selected lineages' frozen annotations, mapping status, both
coordinate systems, source notes and snapshot identity.

V1 captures the original catalog, FN812768.2 reference, landmarks, and eight
legacy geometries. Its frozen exact-site outputs were computed when the
baseline was captured; unresolved legacy placements do not supply platform
eligibility. V2 freezes the expanded panel. Later distinct content receives the
next ordinal, while an identical rebuild reuses its content identity.

Each snapshot includes sequences, alignment, taxonomy, catalog, landmarks,
pair results, per-template results, and source/code provenance. Artifacts have
SHA-256 checksums. Historical selection reads that snapshot's files together,
without joining today's catalog or classes into an older alignment.
Downloads include the snapshot identity; the combination CSV includes the
chosen ceilings. The version-change view compares accessions, represented
classes, and changed shared-pair length evidence against V1.

Release workers restore the previous complete panel and frozen history before
acquiring new records. Previously retained accession versions are preserved;
a conflicting sequence for the same identity stops the build. New distinct
records trigger a complete alignment rebuild. Unreachable or corrupt hosted
state stops an update instead of silently resetting a growing reference.

R2 releases use immutable identities. Publication verifies the anonymous public
manifest and each artifact's SHA-256 before any latest pointer is promoted.
A manual staging run defaults to leaving the current pointer active. The app
uses verified hosted history when configured and retains bundled snapshots as
an outage fallback. COI's broader growing-reference requirements are tracked
separately in [growing-reference-contract.md](growing-reference-contract.md).

## Rebuild and checks

```sh
python3 scripts/build_fungal_reference_panel.py
Rscript scripts/build_fungal_snapshots.R
Rscript tests/test_fungal_explorer.R
Rscript tests/test_fungal_alignment.R
Rscript tests/test_fungal_alignment_regions.R
Rscript tests/test_fungal_amplicon_bp_view.R
Rscript tests/test_fungal_versions.R
Rscript tests/test_fungal_workspace.R
python3 -m unittest discover -s tests -p test_fungal_reference_panel.py
Rscript scripts/build_release.R ITS_FUNGAL
Rscript scripts/validate_release.R ITS_FUNGAL
```

Use `--offline` for a cached acquisition rebuild and `--previous-dir` for a
restored growing panel. The pinned snapshot can be explored without reacquiring
NCBI data or installing MAFFT.
