# Focused investigation pages

The Primer Atlas now separates primer geometry, taxonomic reference support and reference records. Primary navigation is Start, Primer explorer, Taxonomic fit and Reference library. Location comparisons, temperatures and sources remain available under More, with their existing route identities and evidence readers preserved.

Start offers three marker cards: animal COI, fungal rDNA and the PR2 18S catalog. The explorer defaults to its map and selected products; filters, candidate tables, alignment internals, provenance and exports open on demand. Support logos and update diagnostics are available under Sources rather than occupying every investigation page. The COI fit view starts with five primer columns and a closed reading guide; selecting a heatmap cell opens its comparison inspector.

## Fungal reference support

Fungal investigation now has two explicit collections. **Broad ITS / UNITE** is the default for regional alignments, primer-site investigation and reference breadth. **Long rDNA / NCBI** retains the existing versioned 18S/28S maps and their gene/domain annotations. The two collections do not substitute fragments, coordinates or taxonomy for one another.

The first broad snapshot, `u1-afd349e6cdcd`, retains 15,402 references selected from all 102,137 paired standard/developer records in the pinned [UNITE general fungal FASTA release, 19 February 2025](https://doi.org/10.15156/BIO/3301229). This downloadable UNITE 10.0 release is distinct from the current live database. The source has 19 resolved phyla, 86 resolved classes and 293 resolved hierarchy-scoped order groups; unresolved and incertae groups remain explicit. Up to 100 distinct SH/sequence references are selected per order group, and every source-limited shortage is retained. Three conflicting standard/developer sequence relationships are withheld rather than reconciled by guesswork.

The developer source empirically retains additional flanks in 44,884 matched records. ITSx 1.1.3 with the fungal HMM profiles detects 15,286 of the 15,402 selected records. Complete, verified extracted intervals are available for 5,132 ITS1, 14,601 5.8S and 6,868 ITS2 records. SSU/LSU intervals are explicitly retained flanks, not claims of complete genes. All accepted extractions are checked against oriented raw source slices; missing emitted intervals, problematic notes and unknown boundaries remain separate. The archived parser and original tool execution are identified independently, including upstream batch-output quirks.

There are 1,022 independently built MAFFT regional alignments within resolved classes or orders. Each contains at most 32 diverse comparison exemplars; the interface shows at most 12 rows and a movable 60/80/120-column window. Display membership, complete-region membership and the full retained-panel denominator are reported separately. A missing region is not fabricated as an all-gap row. All 17,691 aligned comparison rows retain every extracted source base.

The broad collection freezes 1,725,024 exact-site searches for 112 primer components. Sparse absent rows are interpreted as no exact hit only after checking a completed exhaustive search proof against the sequence, catalog, membership and artifact hashes. The existing 884 pair identities use those indexed sites; all 15,402 references remain in their denominators. An exact-site result does not resolve whether a missing hit reflects a truncated record or a primer mismatch. For example, ITS1 + ITS4 has 1,633 exact inward products, 3,726 single-site records and 10,037 records with no exact sites; six further records have unresolved geometry/site evidence. This is reference evidence, not an amplification probability.

Source census, both source headers, raw published accessions, SH identities, original/oriented sequence hashes, taxonomy status, annotation evidence, tool identities, complete search proof and every regional alignment are frozen together. `data/pinned/UNITE_ITS/index.json` and `latest.json` retain immutable `uN` versions; historical selection appears inside the closed source details when additional versions exist. The source is attributed to Abarenkov and colleagues (2025), under CC BY 4.0. The original archive is identified by URL and checksum, and the full identity census is included in the snapshot.

The pipeline entry points are `build_unite_reference_panel.py`, `annotate_unite_regions.py`, `scan_unite_primer_sites.py`, `build_unite_regional_alignments.py` and `freeze_unite_snapshot.py`. Annotation requires the explicit verified ITSx/HMMER tool manifest; cached work can be reused only with matching inputs, tool/profile identities and captured command semantics. Frozen outputs supply the app without a runtime alignment-tool installation.

The following exact-product view describes the separate **Long rDNA / NCBI** collection:

Taxonomic fit has explicit Animal COI and Fungal rDNA modes. Fungal mode uses the selected immutable snapshot's raw sequences, taxonomy and exact-site product records. Phylum and class are always available with unresolved assignments retained; order and family are offered only when the snapshot has actual resolved labels at those ranks.

Each clickable taxon card shows exact product support as a count divided by **all retained reference templates in that hierarchy-scoped group**. Same-name taxa under different parent lineages remain distinct. Drilling down opens accession-level statuses and reasons. These states partition the denominator:

- exact product;
- one exact site / partial pair;
- ambiguous or otherwise unresolved sites;
- no exact sites found in the record;
- exact sites whose geometry or minimum-length rule withholds a product;
- missing, inconsistent or unsupported evidence.

Absence of an exact site in a retained record does not establish PCR failure or absence of the locus. This bounded fungal panel is not a representative survey of all fungi, and the view does not borrow COI mismatch penalties or amplification probabilities. The taxonomy-group count explicitly includes unresolved assignments.

Platform ceilings use encoded ungapped reference spans including both primers: 500 bp Illumina, 2,000 bp PacBio and 3,500 bp Nanopore. Selecting a ceiling never changes the exact-support numerator or its panel denominator. The planning count instead reports eligible encoded spans amongst exact-supported products. Internal ambiguous bases are disclosed; an encoded span cannot determine an assembly gap's unknown physical length.

## Context and historical versions

Explorer → fit and fit → explorer carry the exact version and pair through asynchronous Shiny control updates. An empty explorer selection produces a choose-combination notice and retains the current route, rather than opening unrelated evidence. User selections and filter resets cancel stale handoff acknowledgments.

The fungal reference library displays the full frozen panel from the fungal fit view. Metadata downloads retain each accession, source fields and snapshot ID. Taxon-count and reference-evidence CSVs use the displayed version, pair, rank and search scope. A historical V1 view uses its one retained record and unresolved taxonomy; the current V3 view uses its 21 records. Frozen snapshots and monthly acquisition/scoring logic are unchanged by this interface work.

## Validation

`tests/test_fungal_fit.R` reconciles all 18,564 V3 pair-template rows across 884 pairs and verifies status partitions, hierarchy ownership, unavailable evidence and historical behavior. Rendered module tests cover count bars, taxon selection, independent ceiling changes, empty searches, version changes and handoff release. Real application navigation tests cover both fit modes, historical pair/version round trips, actual versioned library CSV callbacks, preserved COI state and the original four mapped 18S defaults. Explorer tests cover delayed client acknowledgments and later direct selections.

UNITE regressions cover import eligibility, deterministic quotas, taxonomy ownership, partial and problematic regions, strand and source slices, exhaustive sparse-search proof, all frozen legacy-engine site/pair results, class/order alignment ownership, full FASTA/CSV exports, historical selections, interrupted snapshot registration and corrupt-state rejection. The first broad snapshot independently verifies 2,085 artifact hashes and all regional bases. Its 32.1 MiB payload brings the copied Connect application bundle to approximately 54.9 MiB, below the existing 100 MiB gate. NCBI V1–V3 sequence, scoring and annotation bytes remain unchanged.
