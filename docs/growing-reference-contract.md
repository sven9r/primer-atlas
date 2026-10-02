# Growing COI references, families, and ecological comparisons

**Requirement clarified and corrected by the user:** 2026-10-02.
**Implementation status:** required next scientific milestone; incomplete in
the published 2026.10.02 app. This document defines the work and its acceptance
criteria; it does not claim that a new sequence pipeline has been deployed.

## Core purpose

Every monthly release must extend the retained reference sequences, update
their alignments, and evaluate primer fit for each supported order and at
least its families. Ecological questions then refine that evidence by family,
region, and finer taxa when available. Family evaluation belongs to the core
pipeline and must work across every supported order.

The retained bee alignment is a foundation for **BeePrime only**. Its evolving
bee reference can acquire further eligible bee accessions. General arthropod
primers require their own broad arthropod reference across supported orders
and families. Aquatic arthropods require a separate selectable ecological
scope. The retained bee panel must not be used as the shared reference for
general or aquatic arthropod primers.

## Reference routing required for a useful decision tool

| User question / primer scope | Reference used | Required refinement |
| --- | --- | --- |
| BeePrime / bee communities | Retained bee reference plus additional eligible bee sequences, preserving its original study panel. | Bee families, finer taxa and source-supported regions. |
| General arthropod primers / all arthropods | Broad arthropod reference across every supported arthropod order. | Order, family, region and application; disclose lineage gaps and exclusions. |
| Aquatic arthropod communities | Explicit aquatic arthropod scope using eligible aquatic taxa in the broad reference or a documented dedicated panel. | Freshwater/marine, aquatic family or finer taxon, region and sampling context. |

The tool must combine the ecological target, documented primer use, and
matching scoring reference. Changing a target must change eligible primers
and the reference subset used to calculate its summaries. A primer name or
an environment tag alone cannot establish reference suitability.

Monthly acquisition grows the appropriate references while preserving their
identities and source versions. Independent broad arthropod acquisition can
include Hymenoptera; that does not authorize substituting the retained
bee-focused study panel as a broad reference. Off-target results, if exposed,
have a separate evidence role and denominator.

Aquatic membership needs explicit taxonomic/ecological evidence. Whole orders
such as Diptera or Coleoptera contain mixed habitats, so an order-wide
`environment_group = both` label is only a discovery hint. Preserve the
family, finer taxon or life-stage resolution that supports aquatic membership
and retain uncertain assignments. Aquatic arthropods are distinct from fish,
mollusks, and the broader zooplankton grouping already in the app.

## Required outputs for each monthly release

| Layer | Required result |
| --- | --- |
| Sequence acquisition | Retained previous accessions plus new eligible accessions, with sequence checksums, source, first release, and acquisition/QC status. |
| Alignment | All retained eligible sequences aligned to the versioned coordinate reference; new input must invalidate a stale alignment. |
| Taxonomy | Accession-to-accepted-taxon hierarchy including actual order and family, plus unresolved/conflicting assignments and taxonomy source/version. |
| General primer fit | General arthropod COI pairs scored on the growing broad arthropod panel, with order and family summaries and exact sequence drill-down. BeePrime uses its separate bee reference. |
| Coverage | For every supported order, report families represented, poorly represented, unresolved, and lacking eligible references. Preserve zero-reference families from an explicit taxonomic inventory. |
| Geography | Location metadata for the retained accession population, with country, locality and coordinates where reported, explicit missingness, and release provenance. |
| Ecological comparisons | Route to the correct reference for bee, general arthropod or aquatic arthropod questions, then recalculate summaries from its accession-level records by order, family, region and finer supported taxa. |
| Reproducible state | An immutable release containing the retained sequence/state files and their checksums, so the next runner can reconstruct the complete previous panel. |

An order-wide total is insufficient to establish family coverage. Additional
sequences should be able to fill missing families and regional gaps rather
than only increase the size of already common lineages. The existing minimum
1% per-order growth policy is a volume floor, not evidence that every family
or ecological region has been sampled adequately. The working planning
assumption is to fill family and regional coverage gaps first, then increase
overall counts. An optional user question about that acquisition priority is
pending; this assumption does not change the current 1% pipeline settings.

“Every family” means every family in the explicitly versioned inventory for
the supported order. A family with no eligible sequence must appear as a
coverage gap; it cannot receive a fabricated fit score. Acari and Collembola
are composite display groups: retain their constituent orders and families
so the hierarchy remains biologically meaningful.

## Bee foundation and continuing growth

1. Retain the original Gurten et al. alignment, clustering information,
   source version and study-specific BeePrime scores for reproducibility.
2. Use its accession/taxonomy information to seed or reconcile the evolving
   BeePrime bee reference, recording how study sequences enter that panel.
3. Add further eligible bee accessions each month, including new taxa and
   regions. Identify bees from taxonomic evidence and retain unresolved
   assignments explicitly.
4. Show both the evolving BeePrime bee panel and its original study panel,
   with a clear source and denominator for each. Retained bee-reference
   scoring is exclusive to BeePrime; general arthropod primers use the broad
   arthropod reference.
5. Deduplicate accession identities across sources. Publication cluster sizes
   describe historical represented records; they are not counts of newly
   downloaded sequences and must not inflate monthly growth.

The BeePrime detail, lineage, and regional views share a **Reference
alignment** choice. The publication option reproduces the 590 mapped bee
centroids. The expanded option retains those exact accessions and adds each
new sequence only after its NCBI taxonomy lineage places it within Anthophila.
The release publishes a checksummed aligned FASTA, accession-level score
table, geography table, added-sequence FASTA, and a taxonomy audit. The UI
reports the reference version and total added bees, changes score denominators
and downloads with the selection, and withholds the expanded option's results
when its release artifacts are absent or fail validation. This selection is
scoped to BeePrime; it does not reroute any general arthropod primer.

Short COI/barcode records can be useful for ecological coverage. Define QC
and inclusion separately from primer-site scoreability: a sequence that lacks
a binding site can contribute to reference coverage while its pair penalty
remains unavailable. Do not make the current 1,200–1,800 bp query a claim that
shorter bee or other arthropod records are unsuitable for all purposes.

## Regional work required for ecological questions

- Use the selected scoring reference's own accession geography and scores. The current
  geography fetch is based on the fixed Gurten centroid list; it must be
  extended to new references across the supported orders.
- Provide family selection within each order and regional comparison within
  the same order/family. Report families present in only one subset and
  unequal taxonomic composition before interpreting regional summaries.
- Recompute penalties after filtering the actual sequences. Adding local
  counts beside an unchanged global family score is not a regional score.
- Show total retained, located, region-matched, family-assigned, and
  both-site-scorable denominators. Preserve missing location and taxonomy.
- Support the documented geographic hierarchy (including Hawaiʻi islands
  and the archipelago) when the source resolution permits it; do not infer
  finer locations from a country label.
- Keep specimen locality, reviewed species-range assertions, and empirical
  PCR evidence as separate source-backed fields. Regional sequence fit is
  conditional on the reference sample, not an amplification probability or
  an endemicity conclusion.

The 805-record Crambidae audit remains a useful bounded pilot. Completion of
the regional system requires reusable comparisons across orders and families;
the pilot alone does not complete that requirement.

## Code and public-artifact audit on 2026-10-02

- `.github/workflows/monthly-release.yml` restores the accession ledger from
  mutable R2 state, but does not restore its retained raw FASTA files.
  `scripts/update_order_reference_panels.R` uses the ledger to calculate
  prior counts while its retained sequences come from local files. A fresh
  runner therefore cannot safely reconstruct the previous sequence panel.
- `scripts/build_22_group_alignments.R` may reuse an existing alignment
  without checking whether its source sequences changed.
- `scripts/build_order_scores.R` creates order summaries and full-order
  sequence scores, but does not attach family taxonomy.
- `scripts/build_expanded_beeprime_reference.R` now builds a separate
  BeePrime-only panel by retaining the 590 publication accessions and adding
  current Hymenoptera records whose NCBI lineage is under Anthophila. The app
  offers publication and expanded alignments from the release manifest; this
  does not address the missing taxonomy or family partitions for general
  arthropod panels.
- `scripts/build_release.R` publishes `families`, `subfamilies`, `genera`,
  `taxonomy`, and `geography` from the fixed Gurten-derived tables. These are
  not partitions of the growing full-order panel.
- `app.R` routes BeePrime exact evidence to its study panel, which is the
  correct reference separation. General arthropod family views still need
  partitions of their own broad growing panels. Regional joins use fixed
  study-panel geography, and some geographic family views retain global
  scores while adding local counts. Aquatic arthropods have no explicit
  dedicated choice; the existing environment facet does not complete that
  ecological reference-routing requirement.
- The public [COI manifest](https://pub-f0382349aff94e77a952e7d562ded123.r2.dev/releases/COI/latest.json)
  names release `36864394345-COI`. Its full-order BeePrime score artifact has
  11,116 rows with accessions but no `family`, `genus`, or `taxid` columns.
  Its growth QA records zero previous retained accessions for all 22 groups;
  this release cannot prove cross-month persistence. No immutable reference
  state is described by that manifest.

## Completion gates

1. Run two consecutive builds from clean runner directories. The second
   restores and retains every prior accepted accession and its sequence;
   state/download errors must not silently reset the panel.
2. Add a new bee accession to the BeePrime reference and broad-reference
   records from at least two orders. Verify the respective primer/family
   summaries update, the original bee study panel is reproducible, and the
   bee panel is never used as the scoring reference for general primers. Verify
   the app's publication/expanded selector changes the BeePrime evidence and
   denominators without changing other pairs.
3. Require family coverage output for every supported order, with explicit
   absent/unresolved/poorly represented families. Family totals reconcile to
   order totals and unknown taxonomy is accounted for.
4. Preserve total and scorable denominators for every pair/order/family;
   absent sites have an explicit reason and unavailable penalty.
5. Use a regional fixture containing different families and missing metadata.
   Verify family-specific comparisons are computed from the selected
   sequences and are not copied from global scores.
6. Test explicit BeePrime, all-arthropod and aquatic-arthropod routing. Aquatic
   comparisons must restrict the eligible accession population by documented
   ecology, support freshwater/marine refinement, and expose unknown habitat.
   Test mixed-habitat orders to prevent treating all their families as aquatic.
7. Publish only after immutable state, full-panel taxonomy, family summaries,
   matching geography, exact evidence, and app filters pass together. A
   successful order-score build alone does not complete this milestone.
