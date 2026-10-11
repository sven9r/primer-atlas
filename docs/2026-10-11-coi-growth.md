# COI acquisition and retention proof — 2026-10-11

The latest verified hosted stage, **38151290190-1-COI**, retains **13,150
arthropod accession-version records**, of which **12,267 pass marker-feature
and taxonomy verification**. Each of the 22 supported groups has at least
**105 verified records**. This exceeds the requested 10,000 verified COI
records overall and 100 per group. The scope is 20 supported orders plus the
Acari and Collembola composite display groups, not every arthropod order.

Two consecutive cloud builds passed, including recovery of the first hosted
cohort on a fresh worker, new acquisition, full alignment/scoring, immutable
publication and anonymous checksum verification. This exercises the monthly
workflow now; it does not claim that two calendar months have elapsed. The
local preview at `http://127.0.0.1:3838/` uses the second immutable hosted
cohort without changing the primary local panel or public latest pointer.

A second real acquisition added **139 records** to the frozen initial panel
of 12,721 retained / 11,838 verified records. Every group met
`ceiling(previous retained × 1.01)`. All 12,721 previous sequence/header hash
pairs were unchanged, and 327 cached NCBI source files passed checksum checks.
The primary local panel remains at 12,860 retained / 11,977 verified records.
The subsequent cloud versions are documented below.

## Frozen local proof group counts

Local counts come from `data/derived/coi_reference_proof/retention_proof.csv`
and `data/derived/coi_reference_panel/growth_qa.csv`.
Acari and Collembola are composite display groups; the other 20 are orders.

| Group | Initial retained | Final retained | Added | Final verified marker |
| --- | ---: | ---: | ---: | ---: |
| Amphipoda | 539 | 545 | 6 | 483 |
| Araneae | 1,022 | 1,033 | 11 | 1,020 |
| Blattodea | 133 | 135 | 2 | 131 |
| Coleoptera | 1,042 | 1,053 | 11 | 939 |
| Decapoda | 1,009 | 1,020 | 11 | 736 |
| Diptera | 1,047 | 1,058 | 11 | 1,035 |
| Ephemeroptera | 100 | 101 | 1 | 101 |
| Hemiptera | 1,033 | 1,044 | 11 | 1,028 |
| Hymenoptera | 1,299 | 1,312 | 13 | 1,260 |
| Isopoda | 100 | 101 | 1 | 101 |
| Lepidoptera | 1,299 | 1,312 | 13 | 1,146 |
| Megaloptera | 103 | 105 | 2 | 102 |
| Neuroptera | 159 | 161 | 2 | 158 |
| Odonata | 103 | 105 | 2 | 102 |
| Orthoptera | 1,010 | 1,021 | 11 | 957 |
| Plecoptera | 507 | 513 | 6 | 471 |
| Psocodea | 100 | 101 | 1 | 101 |
| Siphonaptera | 117 | 119 | 2 | 102 |
| Thysanoptera | 103 | 105 | 2 | 102 |
| Trichoptera | 107 | 109 | 2 | 102 |
| Acari | 1,299 | 1,312 | 13 | 1,307 |
| Collembola | 490 | 495 | 5 | 493 |
| **Total** | **12,721** | **12,860** | **139** | **11,977** |

## Provenance and limits

The preserved accession baseline combines 10,849 published arthropod
accessions with 1,943 local legacy arthropod accessions: their union is 12,224.
Published accession recovery cannot establish the bytes of historical raw
files that were never available publicly. Exact sequence/header preservation
is demonstrated between the frozen local initial panel and its second build.
The 50 local Bivalvia records remain a separate legacy archive; the 267
Bivalvia records in the older public release remain in that immutable release.
Neither contributes to the new arthropod totals.

The 883 inherited records that do not pass current marker-feature QC remain
explicitly distinguishable from the verified records in each cohort. Neither marker
verification nor panel inclusion establishes both primer binding sites or
PCR success. Unknown taxonomy ranks remain unresolved. Query result counts
precede QC, and acquisition is bounded rather than exhaustive or random.

The local `data/derived/coi_reference_proof/retention_proof.json`
records the panel and run-manifest hashes. The final proof panel was promoted
to `data/derived/coi_reference_panel`; its FASTA SHA-256 is
`85ef297a2ba07a7966aca42f5c204d923cb66530408858620492c2ee501bc010`.
The initial acquisition was provisional while the builder was refined; the
second build revalidated retained marker and taxonomy sources and captured
the loaded acquisition-code hash. Generated evidence lives with the release
artifacts rather than being committed as a large Git data bundle.

## Rollout at this validation cutoff

The monthly workflow is scheduled for the **first day of each month at
06:00 UTC**. Immutable full-state restoration, alignment input invalidation,
and coherent app reference loading have local fixture coverage. Full alignment
and scoring of the real panel finished successfully for all 22 groups. The
outputs contain 550 group/pair summary rows, 321,500 pair/template rows,
643,000 primer/template rows, and 14,904,740 position rows. The 321,500
pair/template rows reconcile to 12,860 references × 25 curated pairs. Position
traces are partitioned by pair for app loading.

The COI Primer map shows primer geometry on the fixed NC_001322.1 coordinate
reference. The growing full alignments supply custom placement support, score
denominators and per-template component traces; the map itself is not an
alignment matrix showing all retained records. BeePrime publication views keep
their separate study denominator.

Cloud staging initially failed before acquisition after frozen-seed recovery
passed. Expanded diagnostics identified an R-launched Python library-path
mismatch, which prevented the helper from importing the installed `requests`
package. The launcher now verifies and preserves the selected Python runtime.
[Corrected run 38144690584](https://github.com/sven9r/primer-atlas/actions/runs/38144690584)
completed successfully. It restored the frozen seed, added **145 distinct
accession versions**, and retained **13,005 records / 12,122 verified COI
markers**, with at least **103 verified markers per group**. Anonymous reads of
the actual raw FASTAs and ledgers confirm every seed sequence/header hash is
unchanged and each group meets `ceiling(previous retained × 1.01)`. All **988
hosted artifacts** passed anonymous SHA-256 read-back. The R-launched preflight
and NCBI helper both used Python 3.12.15 with the installed `requests` 2.34.2.
The [immutable manifest](https://pub-f0382349aff94e77a952e7d562ded123.r2.dev/releases/COI/38144690584-1-COI/manifest.json)
and [frozen growth QA](https://pub-f0382349aff94e77a952e7d562ded123.r2.dev/releases/COI/38144690584-1-COI/coi_state/data/derived/coi_reference_panel/growth_qa.csv)
are public evidence. The failed attempts published no release.

A second clean worker restored and extended this complete hosted release.
Its first attempt, run 38150641509, stopped at a syntax error in the
hosted-state restore command before acquisition or publication. That CLI path
was bypassed by the successful seed-based build. The command is now corrected,
covered by an actual subprocess recovery regression, and exercised against
the real hosted release in an independent clean directory. All 921 restored
state artifacts, including 841 source/archive files, match the hosted hashes;
the 13,005 retained records, 12,122 verified markers and exact raw
sequence/header hashes reconcile to the frozen manifest. The replacement
[run 38151290190](https://github.com/sven9r/primer-atlas/actions/runs/38151290190)
completed successfully at 16:29:38 UTC. It added **145 distinct accession
versions**, producing **13,150 retained / 12,267 verified COI records** with
minimum **105 verified per group**. All **1,087 hosted artifacts** passed
anonymous SHA-256 read-back. Its [immutable manifest](https://pub-f0382349aff94e77a952e7d562ded123.r2.dev/releases/COI/38151290190-1-COI/manifest.json)
has SHA-256 `c4641af857e536bde94fb88b436693af6bb47c4d9188a205622bede8efb943b2`;
the [complete frozen growth QA](https://pub-f0382349aff94e77a952e7d562ded123.r2.dev/releases/COI/38151290190-1-COI/coi_state/data/derived/coi_reference_panel/growth_qa.csv)
records all group counts and acquisition limits.

An independent anonymous first-to-second comparison passed at 16:30:58 UTC:
all **13,005 previous sequence/header hashes** and **431 required immutable
scientific source/archive hashes** were preserved exactly. Every group met
`ceiling(previous retained × 1.01)`. All group alignment input hashes changed,
and the second acquisition identifies the exact first frozen run-manifest
hash. The release's exact application reader also passed with 13,150
references and **631,200 ANML position rows**.
The temporary GitHub draft seed release was removed after these recovery
proofs passed; its local archive remains available. Future runs recover
immutable hosted state rather than depending on that temporary seed.

The application reader passed against the first actual hosted release: 13,005
references and the expected 624,240 ANML position rows. A clean 171-file,
16.3 MB Connect bundle also loaded this immutable candidate without any local
raw panel or full alignments; its real Ephemeroptera sequence CSV and full
aligned FASTA downloads matched the frozen cohort. Each session can check a
configured hosted pointer manually and every ten minutes. Complete verified
cohorts replace the active reference together; invalid candidates preserve it.
Custom COI scoring uses the unchanged pinned PrimerMiner evaluator and tables
included with the application, without the installed build-only package.

The actual 171-file Connect bundle passed a same-session first-to-second
hosted refresh with no local raw panel or full alignment files. Changing only
the copied bundle's pointer and triggering the real refresh input updated the
version, all group counts, score denominators and download cache identity.
The complete 22-row, 15-field QA exports reconciled to both frozen cohorts;
Ephemeroptera sequence CSV and exact 1,536-column aligned FASTA downloads
increased from 103 to 105 records and retained all previous sequence/header
hashes. BeePrime publication evidence stayed unchanged.

The live preview's **Sources and downloads** view now shows the second
hosted cohort's whole-group previous, added, retained and marker-verified
counts. Its actual **Complete growth QA CSV** download reconciles to the
frozen hosted QA: **13,150 retained, 12,267 verified, minimum 105 per group,
145 additions, 22 rows and all 15 fields**. Missing values are serialized as
blank cells in the app export; all values otherwise match. The QA export
keeps the complete release scope even when a group is selected. The counts
and downloads follow the same verified session reference after an update.

The public
`latest.json` pointer is unchanged; public app activation
requires [PR #22](https://github.com/sven9r/primer-atlas/pull/22) review and merge.
Family coverage inventories, aquatic reference selection, accession geography,
and ecological comparisons remain open under the
[growing-reference contract](growing-reference-contract.md).
