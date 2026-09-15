# Public Multi-Marker Primer Atlas

An open-source Shiny atlas for exploring where metabarcoding primers bind and
why particular lineage–primer combinations may be vulnerable to mismatch or
missing-reference bias. Animal COI is active; fungal ITS and eukaryotic 18S are
pilots. Vertebrate mitochondrial 12S, bacterial/archaeal 16S, animal
mitochondrial 16S, and 28S have stable registry entries for later releases.

The public architecture separates four layers:

- GitHub versions code, schemas, compact catalogs, manifests, and documentation.
- GitHub Pages is the project and citation landing page.
- Posit Connect Cloud runs the Shiny application from `main`.
- Cloudflare R2 stores immutable marker releases and lazy sequence evidence.

Compact pinned summaries keep the app usable during an R2 outage. Raw
sequences and large generated score tables are intentionally excluded from Git.
See [DEPLOYMENT.md](DEPLOYMENT.md) for the remaining account configuration.

For a concise statement of scope, completed work, operational status, and next
milestones, see [PROJECT_STATUS.md](PROJECT_STATUS.md). Human-readable changes
are recorded in [PATCHNOTES.md](PATCHNOTES.md); the synchronized machine-readable
record is [patchnotes.json](patchnotes.json).

## Marker-aware behavior

Markers and biological targets are separate many-to-many dimensions: fish can
route to 12S, mitochondrial 16S, COI, or complementary 18S evidence, while 28S
can route independently to fungi, phytoplankton, protists, or other eukaryotes.
The interface therefore starts with an icon-assisted organism-group choice and
then shows primary, complementary, active, pilot, or planned loci. Primer pairs
still belong to one marker. Applications are normalized to `barcoding`,
`bulk_community`, `edna`, and `diet`; environments, target taxa, and design
intent are separate many-to-many facets. Filters use OR within a facet and AND
between facets. Original use-case labels and claim wording remain preserved in
the source catalog.

The fungal ITS pilot imports an attributed snapshot of the
[UNITE primer resource](https://unite.ut.ee/primers.php) offline. Its map uses
the 18S–ITS1–5.8S–ITS2–28S landmark architecture on FN812768.2. It does not
scrape UNITE during user sessions, and its current single-reference placements
are explicitly distinguished from expanded fungal lineage evidence. All 121
sequence-valid imported primer records are searchable in the app. Each primer
uses the publication named by UNITE as its citation; UNITE is stored separately
as compilation provenance. Unpublished and unreported references stay visibly
labelled instead of being converted into a generic UNITE paper citation.

The 18S pilot imports PR2-primer 2.1.1 as a versioned snapshot: 321 individual
18S primers and 123 documented primer sets. Ninety-three sets have conventional
forward–reverse geometry that can be drawn on the 1,799 bp *Saccharomyces
cerevisiae* FU970071 coordinate map. All sets remain searchable, including those
without drawable coordinates, and every set retains its original reference and
DOI. The default map highlights TAReuk V4, E572F/E1009R V4, 1391F/EukBr V9,
and Uni18SF/Uni18SR as recognizable starting points rather than universal
recommendations.

The Region Comparison tab joins exact primer-binding sequences to accession-level
reference geography. It compares two countries or territories at each primer
position, retains exact sequences, and reports within-order summaries so a
taxonomic sampling difference is not mislabelled a regional PCR effect. Its
interactive globe shows the two activated subsets and distinguishes plotted
coordinates from the larger located/all denominator.

## COI reference-panel policy

General COI primer pairs are scored against full NC_001322.1-coordinate order
alignments, not the Gurten centroid panel. The initial production target is at
least 1,000 retained sequences for each listed arthropod order or composite
target group when NCBI can supply them. Every monthly release is append-only:
the target becomes `ceiling(previous retained × 1.01)`, and only unseen
accessions are added. Existing accessions are not replaced by a fresh sample.

BeePrime is a named exception: its Gurten et al. (2026) 99.5% author centroids
remain available as study-specific, publication-linked evidence. They are not
silently reused for other primers. Compact 97% alignments bundled in Git remain
UI/regression previews and are not the production scoring denominator. The
machine-readable contract is in
[`data/catalog/reference_panel_policy.csv`](data/catalog/reference_panel_policy.csv).

## Releases

`scripts/build_release.R` creates immutable Parquet artifacts and manifests;
`scripts/validate_release.R` enforces schema, IUPAC, citation, checksum,
coordinate, count, completeness, and accession gates. The monthly workflow
builds COI and fungal ITS independently, stages artifacts in R2, and replaces a
marker's `latest.json` only after its gates pass. Failed markers leave their
previous pointers active and open a GitHub issue.

The app deliberately keeps these concepts separate:

1. **Primer geography** — alignment-derived sites on a versioned, full-length
   COX1 reference.
2. **Amplicon + primers** — the full PCR product including both primer
   sequences.
3. **Amplicon − primers** — the primer-excluded, informative sequence used for
   downstream taxonomic inference.
4. **Reference coverage** — whether an order consensus contains sequence at
   both primer sites.
5. **Mismatch penalty** — the raw PrimerMiner position/type/adjacency score,
   never labelled as a PCR probability.
6. **Thermal compatibility** — primer Tm intervals and a configurable,
   transparent annealing-temperature starting window.

## Open the app

From this project folder, run:

```bash
Rscript scripts/run_app.R
```

The app opens in your web browser. If it does not open automatically, copy the
local URL printed in the terminal (normally `http://127.0.0.1:3838`) into a
browser.

Required app packages are `shiny`, `bslib`, `readr`, `dplyr`, `stringr`, `DT`,
and the project-local pinned `PrimerMiner`.
The expanded-reference rebuild additionally uses `jsonlite`, `readxl`, and
`xml2`; empirical-table extraction uses Python `python-docx`.

Run the custom-primer regression check with:

```bash
Rscript tests/test_custom_primer_gate.R
```

It confirms that every curated COI pair passes the placement gate, BeePrime is
rescued by homologous order-alignment evidence, lowercase alignment data remain
valid, and an unrelated 28S D2 pair is rejected with diagnostic feedback.

Run the expanded BeePrime reference checks with:

```bash
Rscript tests/test_beeprime_deep_reference.R
```

Run the target-claim, ZBJ-orientation, site-completeness, Acari/Collembola, and
Lepidoptera-group checks with:

```bash
Rscript tests/test_claimed_primer_reference.R
```

## BeePrime taxonomic drill-down

For **Order lens → Hymenoptera → BeePrime**, the app explicitly separates two
reference panels:

- The legacy order panel contains 50 randomly retrieved Hymenoptera records,
  clustered to 45 anonymous 97% consensuses. It remains visible as a
  reproducibility trace, but is not suitable for a bee-wide performance claim.
- The expanded panel uses Gurten et al. (2026) Supplement 5: 67,352 aligned COI
  centroids clustered by the authors at 99.5% from 316,254 NCBI sequences.
  Taxonomy from the authors' supplements and accession-level NCBI lineage
  records makes family, subfamily, and genus results traceable.

The expanded reference contains 590 mapped bee centroids representing 1,305
source records and all six listed bee families. Thirty of the 42 bee genera in
Supplement 7 have recoverable centroid labels; the other 12 are shown as
explicit mapping gaps instead of being omitted. Both BeePrime sites are
scorable in 541/590 centroids (91.7%), and the median pair penalty is 51.69.
The table flags small reference groups rather than treating them as stable
estimates. It can also show 54 Hymenoptera families and broad taxonomic
lineages (bee, ant, aculeate wasp, sawfly/woodwasp, and other Hymenoptera).
These are lineage groups—not ecological guild assignments.

Wet-lab results remain a separate evidence layer. The article reports 33/38
bee extracts amplified (86.8%); its Supplement 1 table contains 32/37 listed
target rows. Both sources identify five failures. The app retains this
one-row source discrepancy instead of silently reconciling it. The authors also
warn that their empirical validation is taxonomically and geographically
focused on Central European bees.

## Target claims, reference suitability, and ZBJ

The taxonomic drill-down uses a primer-appropriate evidence layer. ZBJ-Art and
its degenerate comparison retain their separate site-complete COX1 panel;
BeePrime retains the Gurten study panel. Other primer pairs do not inherit a
Gurten family/genus score: those views stay unavailable until the append-only
full-order artifacts and their taxonomy partitions have been published. A
published target claim is stored as a hypothesis in
`data/primer_target_claims.csv`; hybrid pairs do not inherit validation from
their component-primer papers.

Before any penalty is interpreted, the app audits whether both binding sites
are actually present. The 67,352-centroid Gurten alignment is suitable for the
later Spidprey/NoSpi/NoPlant sites (about 98% pair-scorable), but fewer than 1%
of its centroids span the early ZBJ/ANML forward sites. Those sparse values are
labelled **reference-site failure**, their penalties are withheld, and they are
never presented as primer coverage.

ZBJ therefore has a separate, reproducible site-complete reference. The build
samples complete mitochondrial records across NCBI search results, extracts
their annotated full COX1 coding features, aligns them to NC_001322.1, retrieves
accession-level NCBI taxonomy, and retains every raw denominator. The current
targeted panel contains:

- 198 Hymenoptera sequences from 137 species;
- 241 Lepidoptera sequences from 129 species;
- 228 Coleoptera sequences from 211 species;
- 244 Acari sequences from 134 species; and
- 61 Collembola sequences from 56 species.

Both ZBJ sites are scorable in 92.5–100% of these sequences. Original ZBJ shows
a high median mismatch penalty for Hymenoptera (653.86) but a low median for
Lepidoptera (14.65), resolving the earlier false 0.8% “coverage” result.
Coleoptera is no longer judged from the site-incomplete Gurten subset:
211/228 complete beetle COX1 records contain both scorable ZBJ sites, with an
original-ZBJ median penalty of 219.2. This establishes that the binding sites
are represented while retaining the substantial family- and species-level
variation visible in the drill-down.
Family and genus tables show which lineages drive each result; for example,
Apidae is strongly represented and has a high original-ZBJ median penalty.
These are in-silico mismatch-risk patterns, not claims that ZBJ cannot amplify
bulk samples.

The app also includes the more-degenerate ZBJ variant from Elbrecht et al.
(2019). The reverse sequence is stored in synthesized 5′→3′ orientation as
`WAYTARTCARTTWCCRAAHCCHCC`, not copied left-to-right from the paper's
left-pointing alignment arrow. Degeneracy lowers many in-silico mismatch
penalties, but the app explicitly retains Elbrecht et al.'s experimental result
that this variant had lower amplification efficiency.

Acari and Collembola are displayed as composite study groups rather than
mislabelled as single modern orders. Lepidoptera is split into butterflies,
a Macroheterocera macro-moth core, an operational (non-monophyletic)
microlepidopteran grade, and convention-sensitive boundary families. Exact
family assignments and unresolved records remain visible.

Family, genus, and subfamily summaries are searchable, filterable by column,
and sortable from their headers. Selecting a row opens the underlying
sequence-level evidence: organism label, representative accession, cluster
size, both primer-oriented binding sequences, mismatch counts, terminal-3′
mismatches, and pair penalty. Accessions link to their versioned NCBI records,
and the selected evidence can be downloaded as CSV or exact FASTA.

For the Gurten reference, “20 centroids · 69 records” means that 20 exact
99.5%-clustered representative sequences were scored and their author-supplied
cluster sizes sum to 69 source records. Supplement 5 provides the representative
accessions and cluster sizes but not every cluster-member identifier, so the app
does not pretend that all 69 underlying accessions can be reconstructed.

When a reference lacks sufficient sequence at one binding site, the app shows a
single explanation and `—` in unavailable penalty columns. These cells do not
mean zero penalty or primer failure.

## Try a custom primer pair

Open **Primer map → Add a custom primer pair** and enter both oligos in their
synthesized 5′→3′ orientation. Primers must contain 15–80 IUPAC DNA bases;
ambiguity codes are supported. An optional expected primer-inclusive length can
disambiguate repeated COI motifs. Each accepted pair is appended to the current
map, so several custom pairs can be compared in one session; **Clear custom
pairs** removes them together. The **Custom-pair statistics** selector keeps
the diagnostics accessible for every accepted session pair, not only the most
recent one. An optional expected target order adds a target-specific warning.
Every accepted custom pair is also added automatically to the multi-primer
selector and single-pair inspector in **Order lens**.

The app does not treat the least-bad local alignment as a binding site. The
first evidence route judges each selected site on the coordinate reference and
requires:

- at least 80% IUPAC-compatible identity;
- no more than one mismatch among the terminal three primer bases;
- at least one uninterrupted 7-bp compatible run; and
- a chance-match expectation of no more than 0.05 across NC_001322.1.

Because NC_001322.1 is a *Drosophila* sequence, a taxon-specific COI primer may
legitimately fail that first route. The app therefore checks the leading
coordinate candidates against the installed COI reference-group alignments. A site can
be rescued only when at least two templates and at least 10% of the available
templates in one order support it, and both primers must be supported together
at the selected coordinate pair. If the user supplies a primer-inclusive
product length, the alignment-rescue identity floor is relaxed to 70% while
retaining the paired-site, terminal-three, compatible-run, and multi-template
requirements. Without a product length, the 80% identity floor remains.

If neither evidence route supports both oligos, nothing is added. The feedback
shows the best accidental sites, incompatible-base counts, compatible identity,
terminal-three mismatches, longest compatible runs, chance-match expectations,
and the failed rules. This accepts the bee-specific BeePrime pair using
homologous COI evidence while continuing to reject an unrelated 28S D2 pair.

For accepted pairs, the app reports:

- forward and reverse alignment coordinates;
- amplicon + primers, amplicon − primers (informative), and Folmer-overlap
  lengths;
- reference and terminal-five mismatch counts;
- Tm intervals at 50 mM monovalent salt;
- leading paired and individual-primer placements; and
- explicit warnings when multiple sites are equally supported; and
- per-order counts of templates supporting both primer sites; and
- the complete on-demand PrimerMiner order matrix, example-template trace, and
  downloadable position-level components.

The order summary reports the strongest compatible orders, whether one order
appears preferred, and orders with ≤5% compatible templates. These are
alignment warnings rather than amplification claims. All plotted positions
remain on NC_001322.1; the order alignments provide behind-the-scenes
compatibility evidence and never replace the coordinate system.

Primer names sit beside their binding arrows on the map. Start/end COI
coordinates, both amplicon sizes, and clickable citations are shown in the same
row. A shared forward/reverse source is purple; different forward and reverse
sources are blue and red, respectively.

Custom inputs and their on-demand reference-group PrimerMiner scores are session-only
and are not silently added to `data/primers.csv`. Passing either gate is
evidence of a plausible COI location—not proof of PCR amplification or
taxonomic specificity. A permanent library addition should still be rebuilt
through the offline pipeline so its scores and provenance become reproducible
project data.

BeePrimeF/BeePrimeR are included in the curated library from Gurten et al.
(2026). On NC_001322.1 the pair maps to 132–161 / 361–383, giving 252 bp
including primers and 199 bp between primers. The source figure reports a
198-bp marker; the one-base difference is retained transparently rather than
overwriting the alignment-derived coordinate calculation.

The spider-prey library distinguishes the adapted Spidprey_F/Spidprey_R pair
(Melcher et al. 2024) from the original NoSpi2_F/Laurelin_R combination. The
latter keeps forward and reverse provenance separate: Lafage et al. (2020) and
Cuff et al. (2021), respectively.

## Rebuild everything

The checked-in derived tables make the app directly runnable. To reconstruct
them from public sequence data:

```bash
# macOS, once
brew install vsearch mafft

# project-local PrimerMiner 0.22 and its minimal NCBI-only dependencies
Rscript scripts/install_primerminer.R

# 20 orders plus Acari and Collembola: NCBI -> PrimerMiner/VSEARCH -> MAFFT
GB_SUBSET=50 Rscript scripts/build_22_group_alignments.R

# exhaustive, alignment-derived primer coordinates and Folmer overlap
Rscript scripts/build_primer_positions.R

# PrimerMiner order, template, primer, and per-position components
Rscript scripts/build_order_scores.R

# checksum-verified ODbL BeePrime supplements and centroid FASTA
Rscript scripts/download_gurten2026_reference.R

# accession-level NCBI family/subfamily/tribe lineage records
NCBI_EMAIL=you@example.org Rscript scripts/fetch_gurten2026_bee_taxonomy.R

# Gurten et al. 99.5% centroid alignment -> BeePrime lineage summaries
Rscript scripts/build_gurten2026_beeprime_scores.R

# audited target-claim summaries on the Gurten alignment
Rscript scripts/build_claimed_primer_scores.R

# ZBJ-specific site-complete COX1 features and family/genus summaries
NCBI_EMAIL=you@example.org Rscript scripts/download_targeted_complete_cox1.R
NCBI_EMAIL=you@example.org Rscript scripts/build_targeted_zbj_reference_scores.R

# Supplement 1 DOCX -> empirical species/genus validation tables
/path/to/python3 scripts/extract_gurten2026_empirical_table.py
```

Set `ORDERS=Plecoptera` for a one-order pilot or `FORCE=true` to rebuild
existing group files. The reference-group list and its rationale are explicit and
editable in `data/orders_22.csv`.

## Deep-research catalog additions

The 2016–2026 literature audit added ten sequence-resolved pairs that were
missing from the first atlas build:

- BF2 + BR1, BF1 + BR2, and BF3 + BR2;
- fwhF1 + fwhR1 and fwhF2 + EPTDr2n;
- HexCOIF4 + HexCOIR4;
- BR5 (B + ArR5) and F230R (LCO1490 + 230_R);
- FishF2 + FishR1; and
- MollCOI253.

Each pair has primary-source provenance, reported product and annealing
metadata where available, alignment-derived geometry, all 23 reference-group
scores, and the 67,352-centroid family/genus/accession drill-down. Fish and
mollusk claims are explicitly marked as off-target-only in the shared
arthropod centroid panel. A new 50-record Bivalvia build supplies the
coordinate-gate evidence needed for MollCOI253 without pretending that the
arthropod panel validates mollusk coverage.

## What is derived—not hardcoded

`data/primers.csv` contains primer oligos and reported metadata, but no plotting
coordinates. For every primer, `scripts/build_primer_positions.R` exhaustively
searches all 1,536 positions of the versioned *Drosophila yakuba* COX1
reference (NC_001322.1, feature 1474–3009), respecting IUPAC ambiguity.
Candidates are ranked by:

1. total incompatible bases;
2. incompatible bases in the terminal five positions;
3. a position/type-weighted localization score.

The script records the selected location, five leading candidates,
reference mismatches, ambiguity count, and PrimerMiner reference penalty.
If independently selected sites are reversed or inconsistent with the
published product length, the pair is re-evaluated jointly against
forward/reverse order and both primer-inclusive and primer-excluded length
conventions. The selected non-leading candidate is retained in the audit table.
LCO1490 and HCO2198 then define the displayed 658-bp Folmer interior between
their 3′ ends. Pair target length, primer-inclusive amplicon length, Folmer
overlap, and y-axis rank are calculated separately from these coordinates.
This reproduces distinctions such as Leray 313/365, fwh2 205/254, and
BF2–BR2 421/461 (targeted region/full amplicon).

The primer map uses a dynamic 100-bp coordinate window. Its default extent is
the Folmer region rounded to the surrounding 0–700 bp ticks; the axis expands
only when a selected curated or custom primer extends beyond that window. Thus
Folmer/HCO2198 selections use 0–800 bp, while a custom primer ending at 1,300
bp expands the view to 0–1,300 instead of permanently reserving the whole gene.
Rows can be ordered by COI binding-site position, primer-pair name, Folmer
overlap, or primer-inclusive amplicon length. Hover text reports normalized
method tags such as bulk, diet, freshwater, and eDNA metabarcoding together
with the target group. Each displayed row receives a dynamic number after
filtering and sorting; hovering the number identifies both the pair label and
its exact forward/reverse oligos. The same current number is prefixed to each
displayed pair in the clickable primer checklist; unselected or currently
filtered-out pairs remain unnumbered. A larger custom tooltip appears when
hovering a primer name, primer bar, arrow, or informative amplicon; clicking
that map element pins or unpins the panel, and Escape closes it. The
**Select all** and **Clear** actions sit below these view controls, before the
amplicon-size filter.

This QA step detected that the initial `fwhR2n` value was stored as a reverse
complement. It is now stored in the published 5′→3′ orientation,
`GTRATWGCHCCDGCTARWACWGG`.

## Current reference build

- 20 arthropod orders, Acari and Collembola composite groups, plus Bivalvia
- 1,993 NCBI COX1 records downloaded in the PrimerMiner group build
- 1,431 PrimerMiner 97% group consensus templates
- 23 reference-coordinate alignments
- 71,550 primer–template evaluations
- 1,658,529 per-position component records
- 35,775 paired-template evaluations
- 575 group × primer-pair summaries
- the separate ZBJ panel retains 972 complete COX1 features across five groups
- 67,352 Gurten et al. (2026) 99.5% COI centroids scored for all 25 curated pairs
- 67,352 Gurten et al. (2026) 99.5% COI centroids scored for BeePrime
- 6,797 mapped Hymenoptera centroids across 54 families
- 590 mapped bee centroids across 30 recoverable genus labels and six families;
  12 additional genera from the 42-genus supplement remain explicit mapping gaps

See `data/provenance/order_alignment_manifest.csv` for exact per-order counts,
queries, paths, status, and reference metadata, and
`data/provenance/software_versions.tsv` for pinned software provenance.

## Score interpretation

The order matrix presents:

- the median raw PrimerMiner pair penalty among templates for which both sites
  are scorable; and
- the percentage of order consensus templates with sequence present at both
  binding sites.

The in-app legend states these denominators explicitly. Green-to-coral colour
is a log-scaled relative encoding capped at the 95th percentile of the
currently displayed primers; hatched cells have no scorable pair. A
synchronized scrollbar remains pinned above the matrix while lower orders are
read. Clicking any scored cell selects that order–primer combination and moves
to its detailed PrimerMiner component trace. Target-specific panels then expose
filterable family/genus summaries and row-level organisms, accessions, and
binding sequences where those reference layers are available.

Gaps and `N` bases are unavailable evidence, not high mismatch penalties.
Forward and reverse penalties, p90 values, terminal mismatches, adjacent
mismatches, and every position component can be inspected and downloaded.
The colour scale is only a within-table visual encoding; there is no hidden
pass/fail cutoff.

PrimerMiner’s default penalty matrices are an in-silico prioritization model,
not a universal probability of amplification. The Stadhouders et al. evidence
supports strong position- and type-dependent effects in a controlled qPCR
setting, but empirical calibration is still required for each PCR chemistry,
primer concentration, template mixture, and annealing protocol.

## Project structure

- `app.R` — Shiny interface
- `R/functions.R` — IUPAC and thermal helpers
- `data/primers.csv` — primer library and source keys
- `data/orders_22.csv` — transparent order/composite-group selection
- `data/primerminer/` — raw, clustered, and reference-aligned order data
- `data/external/gurten2026/` — ODbL author supplements, retained unchanged
- `data/external/targeted_zbj/` — sampled mitochondrial CDS batches, complete
  COX1 features, accession taxonomy, alignments, and retrieval manifest
- `data/derived/` — geometry and complete score components
- `data/provenance/` — software and dataset manifests
- `data/citations.csv` — machine-readable evidence registry
- `scripts/` — reproducible installation, acquisition, alignment, scoring, and
  launch scripts
- `data/legacy_prototype/` and `scripts/legacy/` — preserved first-prototype
  inputs/outputs; the app never reads these files

## Contributing

Primer Atlas welcomes contributions to code, documentation, primer catalogs,
reference panels, tests, and scientific interpretation. Start with the
[contribution guide](CONTRIBUTING.md), use a structured issue form, and make
changes through a reviewed pull request. Early ideas and interpretation
questions belong in [GitHub Discussions](https://github.com/sven9r/primer-atlas/discussions).

Please read the [code of conduct](CODE_OF_CONDUCT.md). Primer and reference
contributions must retain primary sources, exact denominators, reference scope,
and unresolved ambiguity; an in-silico mismatch score is not an amplification
probability.

## License

Code is released under the MIT License. Third-party data retain their stated
licenses and attribution; see the app’s Evidence & citations tab.
