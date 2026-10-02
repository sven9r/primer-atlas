# Geographic range and endemism evidence for the Hawaiʻi–Madagascar insect COI study

Status: foundation only, 2026-09-23. No taxon-level range assertions have been reviewed or joined to COI sequences. The empty evidence table is intentional.

## Study regions and comparison unit

[`regions.csv`](../data/catalog/range/regions.csv) fixes the geographic vocabulary. `HI_ARCH` is the Hawaiian Archipelago, with eight named main-island children and `HI_NWHI` for the northwestern island group. `MG_COUNTRY` is Madagascar at country resolution, including its offshore islands. A Madagascar country record must not be described as a record from the main island. A Hawaiian record with only “Hawaiian Islands” maps to `HI_ARCH`, never to a named island. `HI_NWHI` retains the specific atoll or island in the source note when reported.

Archipelago endemicity and single-island endemicity are different assertions. “Endemic to Hawaiʻi” means restricted to `HI_ARCH`; it does not imply restriction to Hawaiʻi Island (`HI_HAWAII`) or any other child. The planned primary comparison is `HI_ARCH` versus `MG_COUNTRY`; island-level Hawaiian summaries are secondary and require island-specific evidence. These regions are analytical definitions, not new controls in the current app.

## Evidence files and interpretation

[`taxon_region_evidence.csv`](../data/catalog/range/taxon_region_evidence.csv) holds one source-backed assertion per accepted taxon, region, and source record. Its column contract and controlled values are in [`taxon_region_evidence.schema.json`](../data/catalog/range/taxon_region_evidence.schema.json). `accepted_taxon_id` must be a namespace-prefixed ID from `taxon_authority` at `taxon_authority_version`; preserve the range source's original name separately in `source_taxon_name`. Do not force an unresolved name or COI accession into an accepted taxon ID. Source IDs, versions, access dates, scope, and limitations are in [`sources.csv`](../data/catalog/range/sources.csv). Every populated row needs the exact source-record URL, a short verbatim assertion, and a review state. Multiple sources or disagreements remain separate rows.

`establishment_status` is `native`, `introduced`, or `uncertain`. `endemicity_assertion` is `endemic`, `non_endemic`, `unresolved`, or `not_applicable`; use `not_applicable` for introduced taxa because regional endemism is a native-range claim. `confidence` describes strength of the *stated assertion* (`low`, `medium`, `high`), not sequence quality. `review_state` is `candidate`, `reviewed`, or `contested`; a candidate is not publishable as a range conclusion. Species-level endemicity is preferred. Genus and family rows provide context only and may not be copied onto member species.

## Decision rules

1. Normalize the sequence taxon to an accepted name and identifier with an explicit taxonomy source and version. Ambiguous genus-only, synonym-only, or conflicting matches stay outside species-level endemism summaries.
2. Transcribe a source's positive regional occurrence at its actual resolution. An accession collection location, specimen, type locality, or checklist entry establishes a reported record, not native status by itself.
3. Assign `native` or `introduced` only from a source that explicitly supports that status or from a documented taxonomic review. Bishop Museum's `endemic` and `indigenous` map to native; `adventive` and `purposely introduced` map to introduced. Unknown, dubious, and quarantine-only entries remain uncertain pending review.
4. Assign `endemic` only when a source explicitly asserts restriction to the *same defined region* and the current accepted taxon concept has been checked. A source's Hawaiian-archipelago endemic label does not justify an island-endemic label. A Madagascar country list, AntWeb regional list, or Afromoths distribution list is not an endemicity assertion.
5. Assign `non_endemic` only with explicit native occurrence outside the defined region or a source's explicit wider-native-range statement. Never infer endemicity from a lack of records elsewhere. Absence from the other study region is also not evidence of absence globally.
6. Conflicting origin, distribution, or taxonomy sources become `contested` rows with notes; the downstream summary remains unresolved until a reviewer resolves the conflict. Preserve source versions and retrieval dates so later updates can be compared rather than silently overwritten.
7. For comparison, count distinct accepted species with reviewed assertions and report unresolved and unmatched species separately. Do not treat sequence counts, COI centroids, collection effort, or primer binding as species-range denominators. Report taxonomic and geographic coverage alongside any eventual percentage.

## Source selection and limits

- [Nishida (ed.), *Hawaiian terrestrial arthropod checklist*, 4th ed. (2002)](https://hbs.bishopmuseum.org/publications/pub2002.html) is the published regional authority. The [Bishop Museum checklist query](https://hbs.bishopmuseum.org/checklist/query.asp) labels its data revision as 9 April 2002 and exposes mode-of-origin and island fields. The [Bishop Museum explanation of checklist terms](https://data.bishopmuseum.org/HBS/checklist/abtarthrocklist.html) defines endemic as limited to the Hawaiian Islands, indigenous as naturally present with a wider range, and adventive as immigrant. The edition is historic; reconcile names and later revisions taxon by taxon.
- [AntWeb's Madagascar species page](https://www.antweb.org/taxonomicPage.do?countryName=Madagascar&rank=species&statusSet=valid+extant), labelled version 8.114 on access, is a specialist Formicidae distribution and name lead. Its region list combines specimen records, type localities, and curator additions; the page does not by itself resolve native status or endemicity.
- [Afromoths](https://www.afromoths.net/), labelled last updated 16 August 2026 on access, is a specialist Afrotropical moth name and distribution source. Its continuously updated range statements need taxon-level inspection and another explicit basis before native or endemic labels are assigned.

All three web resources and the checklist citation were checked on 2026-09-23. Source-page version labels are recorded as displayed, not presented as immutable downloadable snapshots. The direct Bishop report PDF and individual taxon records were not audited in this pass; source-specific rows remain empty until they are.

## Relationship to the current app and sequence audit

The Shiny Region Comparison currently filters accession-level `country_or_territory` from `data/derived/reference_geography.csv` and compares exact primer-binding evidence with located/all denominators and within-order summaries. Its `USA` subset is not Hawaiʻi-specific, and the view does not use this new range catalog. The COI pinned geography release reports only 0.0074 geography-resolved fraction. No app claims or filters change in this pass.

The current `data/derived/reference_geography.csv` has 67,352 rows, 500 with a country/territory, no USA rows whose raw or parsed locality mentions Hawaiʻi, and one Madagascar country row (checked 2026-09-23). Those counts describe this file, not the full universe of possible COI accessions.

The sequence audit must enumerate eligible insect COI records, check accepted-name and rank resolution, identify which have source geography at the study-region resolution, match those taxa to reviewed range assertions, and report unmatched/uncertain/conflicting counts. Only then can a Hawaiʻi–Madagascar species comparison or endemicity breakdown be computed. It must retain the distinction between accession locations, species ranges, and primer-site evidence.
