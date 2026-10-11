# Primer Atlas project status

**Status date:** 2026-10-11

**Current phase:** Public hosting retains the verified 2026.10.02 app and
October 1 marker releases. The October 11 local COI acquisition now retains
12,860 records / 11,977 verified markers, with at least 101 verified records
in each supported group. A second real acquisition added 139 records and
preserved all 12,721 prior sequence/header hash pairs. Full alignment and
real-panel scoring finished successfully for all 22 groups. A cloud worker
restored the frozen seed and published a verified immutable staging release
with 13,005 retained / 12,122 verified records; all 988 artifacts passed public
checksum checks. A second hosted-state restoration proof and public activation
remain pending. The
broader scientific milestone still requires family coverage
inventories, continued BeePrime-specific bee acquisition, a separate aquatic
arthropod option, and family-level regional comparisons.
The measured local growth and rollout limits are in
[`docs/2026-10-11-coi-growth.md`](docs/2026-10-11-coi-growth.md); the last published
app evidence is in [`docs/2026-10-02-release.md`](docs/2026-10-02-release.md).
The clarified scope, implementation gaps, and completion gates are in
[`docs/growing-reference-contract.md`](docs/growing-reference-contract.md).

## Project goal

Primer Atlas is an open, reproducible Shiny application for investigating where
metabarcoding primers bind and why a lineage–primer combination may be
vulnerable to mismatch or missing-reference bias. It is intended to connect
primer sequences, published claims, coordinate geometry, reference-panel
composition, mismatch evidence, and exact source records without turning any
single evidence layer into an unsupported amplification-probability claim.

Its core sequence reference must grow with monthly releases. Every supported
order needs general primer scoring and family-level evaluation. Ecological
questions refine those same accession-level scores by taxon and region. The
bee study alignment is retained for BeePrime only, with continued bee
acquisition and preserved study provenance. General arthropod primers use a
broad arthropod reference; aquatic arthropods require their own ecological
selection and matching reference subset.

The project separates four operational responsibilities:

- GitHub versions code, schemas, compact catalogs, manifests, tests, and
  documentation.
- GitHub Pages provides the public project and citation landing page.
- Posit Connect Cloud runs the Shiny application.
- Cloudflare R2 stores immutable marker releases and lazily loaded sequence
  evidence.

## Scientific and product principles

- Keep markers and biological target groups as separate many-to-many concepts.
- Preserve the original primer sequences, citations, claims, and ambiguities.
- Distinguish binding-site representation, mismatch penalty, reference
  coverage, empirical amplification, and geographic sampling.
- Show exact denominators and retain accession- or sequence-level drill-down.
- Treat monthly retained-sequence growth, family evaluation for every order,
  and taxonomically controlled regional comparisons as core requirements.
- Use marker-appropriate reference panels; do not silently transfer evidence
  from one primer, lineage, or study panel to another.
- Make generated releases immutable, validated, attributable, and recoverable.

## What has been achieved

### Application and data model

- Built the multi-marker Shiny atlas with marker-aware target routing, primer
  maps, searchable evidence tables, custom-primer evaluation, and
  sequence-level drill-down.
- Normalized primer applications, environments, target taxa, sources, claims,
  oligos, pairs, and reference-panel policy into auditable catalogs.
- Kept organism group, marker, primer pair, region, and evidence type as
  distinct filters rather than collapsing them into a single recommendation.

### Marker coverage

- **COI — active implementation:** curated pair geometry, order-level scoring,
  custom-primer gating, BeePrime study-specific evidence, targeted complete
  COX1 panels for ZBJ, claimed-primer evaluation, and reference geography.
- **Fungal rDNA — bounded versioned explorer:** 20 starter combinations,
  full candidate list, pairing matrix, five-pair comparison, platform ceilings,
  per-template lengths, and class/phylum filters. The selected combinations are
  linked to a shared frozen alignment map and nucleotide inspection window.
  V1 preserves the original
  reference; V2 contains 21 accession-versioned templates. Snapshot history
  freezes catalogs, sequences, alignments, taxonomy, and exact-site evidence.
  The October 11 overhaul is implemented and tested; immutable hosted staging
  and app-loader read-back succeeded. Public deployment awaits PR #22 review
  and merge; the active public pointer was retained.
- **18S — pilot:** PR2-primer 2.1.1 snapshot with 321 individual primers and 123
  documented sets; 93 sets have drawable conventional geometry.
- **Planned registry coverage:** vertebrate mitochondrial 12S,
  bacterial/archaeal 16S, animal mitochondrial 16S, and 28S LSU have explicit
  marker identities and target-routing documentation but are not production
  releases.

### Reference evidence and safeguards

- Established full-reference and study-specific panel policies so Gurten bee
  centroids are not silently reused as general COI evidence.
- Acquired and locally verified a growing general COI panel across 20 orders
  plus the Acari and Collembola composite groups. The goal is 10,000 verified
  COI records overall and at least 100 per group, with retained legacy records
  reported separately from marker-verified records. The October 11 proof
  verifies actual additions, previous sequence/header retention, and 327
  checksummed NCBI source files; it does not establish primer-site coverage.
- Added accession-level NCBI taxonomy with explicit unresolved ranks and a
  coherent reference consumer that validates full alignments, score-table
  membership and denominators. Pair-specific position files load lazily;
  corrupt artifacts cannot silently borrow a different release's evidence.
  Consumer tests pass; the real panel's full alignment and scoring are
  complete, and hosted activation is still pending.
- Added a local BeePrime reference-alignment choice for the unchanged
  publication panel or a Bee-only expanded panel. The monthly builder verifies
  bee lineage with NCBI taxonomy, preserves the exact 590 publication
  accessions, and publishes score, sequence, alignment, and geography
  artifacts. This change is tested locally but not yet published.
- Added explicit site-completeness checks that withhold mismatch penalties when
  a reference does not span both primer sites.
- Preserved the reported BeePrime wet-lab denominator discrepancy and labelled
  empirical results separately from in-silico evidence.
- Added located/all geographic denominators and within-order comparisons so
  reference sampling is not mislabelled a regional PCR effect.

### Reproducibility and QA

- Added reproducible import, download, alignment, geometry, scoring, release,
  validation, and R2 publication scripts.
- Added regression coverage for catalog integrity, marker runtime, release
  layers, reference-panel policy, map domains, custom-primer gates, geographic
  comparisons, order interactions, and sequence drill-down.
- Implemented marker-independent matrix jobs, hard release gates, immutable
  release directories, atomic `latest.json` promotion, QA-artifact upload, and
  automatic failure issues.
- Implemented checksummed immutable COI raw/state restoration, exact ledger
  validation, and input-hash invalidation of changed alignments. Clean-worker
  fixtures verify successive restoration, retention, additions, and rejection
  of incomplete or corrupt state. These checks complement the real local
  acquisition proof; two successive hosted clean-runner releases remain open.

### Collaboration foundation

- Published an MIT software license and `CITATION.cff` citation record.
- Added structured primer, reference-panel, and bug-report forms.
- Documented local setup, evidence standards, tests, pull-request review, and
  community conduct for outside contributors.
- Assigned review ownership for scientific catalogs, provenance, scripts,
  workflows, and the pinned R environment.
- Enabled GitHub Discussions and added the public Pages homepage, repository
  description, scientific topics, and workflow labels.
- Protected `main` with pull requests, code-owner review, one approval, resolved
  conversations, and the required `test` check; force-push and deletion are
  disabled.

## Current operational status

- Public app: https://01a0cdbb-4448-22b4-500c-b5329e3b1904.share.connect.posit.cloud/.
- [Monthly run 36864394345](https://github.com/sven9r/primer-atlas/actions/runs/36864394345)
  completed successfully on 2026-10-01 for both COI and ITS_FUNGAL, including
  hard gates and atomic R2 promotion. The public pointers currently name
  `36864394345-COI` and `36864394345-ITS_FUNGAL`.
- On 2026-10-02 both public manifests were fetched, and a published artifact
  from each marker passed its SHA-256 check. This supersedes the older
  COI-blocked status; the old incident issues still need an evidence-linked
  closure review.
- The 2026.10.02 update fixes marker-switch selection, adds organism and dietary
  study tags, groups sources, and downloads the current mapped-pair metadata.
  PR #20 merged at `648dd55`; both PR and merged-main CI passed all 16 tests.
  Posit publication, public switch/filter checks, and the 126-row/46-column
  metadata download were verified. Details are recorded in the release note.
- The published October 2 runtime excludes build-only PrimerMiner and readxl;
  its custom COI reference scoring is unavailable on Connect Cloud. The working
  branch now bundles the unchanged pinned PrimerMiner 0.22 evaluator, default
  tables and GPL-3 license in an isolated runtime. Numerical parity and both
  custom-scoring paths pass without the installed PrimerMiner or seqinr package.
- Research foundations are versioned independently of public evidence. The
  regional COI audit and 48-accession ITS scoring outputs are not loaded into
  the app or promoted as broad performance results.
- The first-day monthly schedule remains 06:00 UTC. The October 11 COI
  staging attempts initially failed before acquisition after frozen-seed
  recovery passed. Diagnostics identified the Python library selected when
  launched from R; the corrected
  [run 38144690584](https://github.com/sven9r/primer-atlas/actions/runs/38144690584)
  completed successfully: 13,005 retained / 12,122 verified COI records,
  minimum 103 verified per group, 145 additions with exact seed sequence/header
  retention, and all 988 public artifact checksums verified. Those earlier
  failures published no release. A second clean-worker hosted-state restoration
  proof remains pending.
  All 22
  real local full alignments and scoring finished successfully: 550
  group/pair summaries, 321,500 pair/template rows, 643,000 primer/template
  rows and 14,904,740 position rows.
  Public pointers remain unchanged, and PR #22 review/merge is required before
  app activation. Hosted-to-hosted recovery and public app activation remain
  unverified at this cutoff.

## What still needs to be done

### P0 — Growing references and family evaluation across all orders

- [x] Implement immutable retained sequence/state restoration with exact
      sequence/header checks and clean-worker regression coverage.
- [x] Demonstrate real local growth with at least 10,000 verified COI records
      overall, 100 per group, and at least 1% retained growth in all 22 groups.
- [x] Invalidate stale alignments when retained inputs change; test the input
      hash and reference-coordinate contract.
- [x] Attach accession-level order/family taxonomy to the growing panel,
      retaining unresolved ranks and source checksums.
- [x] Complete full reference-coordinate alignments for all 22 real groups.
- [x] Complete real-panel scoring for all 22 groups and 25 curated pairs.
- [ ] Complete two successive hosted clean-runner builds that restore every
      prior sequence and its metadata.
- [ ] Publish coverage gaps and primer summaries for every supported order
      and its families, including zero-reference families from a versioned
      inventory, unresolved assignments and exact denominators.
- [ ] Publish and verify the new BeePrime publication/expanded alignment
      choice. The local build added 30 NCBI-classified bees to the 590-accession
      publication baseline; monthly clean-runner retention and the public app
      still need verification. Preserve the publication panel and never route
      the new bee records to general primers.
- [x] Implement and test coherent general COI reference routing with explicit
      full-panel versus legacy-preview identity and pair-lazy evidence.
- [ ] Verify that routing with the real scored release in the public app,
      including the retained versus marker-verified counts for every group.
- [ ] Add a separate aquatic arthropod option with freshwater/marine
      refinement, matching primer uses, and an ecology-supported reference
      subset at family/finer-taxon resolution. Mixed-habitat orders cannot
      establish aquatic membership for every sequence.
- [ ] Publish and verify family evaluation for every curated pair/order in
      the app. General-panel taxonomy joins and family consumers are tested
      locally; real scored artifacts, coverage inventories and public routing
      still need verification.
- [ ] Extend accession geography to the growing panels and recalculate
      ecological comparisons within order/family and source-resolved regions.
- [ ] Pass the complete gates in `docs/growing-reference-contract.md` before
      claiming that monthly growth or ecological family evaluation is done.

### P1 — Release documentation and operational reliability

- [x] Configure R2 and pass both marker builds, hard gates, and promotions.
- [x] Verify public COI and ITS_FUNGAL pointers and one artifact checksum each.
- [x] Confirm the newly published application and metadata download.
- [ ] Review and close historic release-failure issues with successful evidence.

- [ ] Add a low-cost scheduled or dispatchable smoke test for configuration and
      public manifest reachability.
- [ ] Decide whether repeated matrix failures should update existing marker
      issues instead of opening duplicates each month.
- [ ] Test the pinned-summary fallback during an intentional R2 outage.
- [ ] Document release ownership, credential rotation, rollback, and recovery
      from a partially uploaded staging directory.
- [ ] Verify hosted release recovery retains accession identities, sequence
      bytes and header bytes across two real clean-runner releases. Local
      fixtures and the two-acquisition proof pass; ledger persistence alone
      does not satisfy this operational gate.

### P2 — Finish active pilots

- [x] Merge and verify the ITS primer-first navigator in a clean browser session.
- [ ] Add 18S to the release matrix only after its artifact contract and hard
      release gates are defined.
- [ ] Expand fungal ITS beyond the single pilot reference while continuing to
      label reference scope and missing lineage evidence explicitly.
- [x] Show per-marker readiness and limitations on the app landing page.
- [ ] Complete source-backed organism/application tags for further catalog entries.
- [ ] Reconcile Crambidae taxa and coordinates, then score both binding sites
      with panel and scorable denominators.
- [ ] Populate reviewed taxon-region range assertions; zero reviewed rows exist.
- [ ] Independently annotate fungal ITS boundaries and terminal binding sites
      before expanding performance comparisons.

### P3 — Expand scientific coverage without weakening evidence standards

- [ ] Expand taxonomic depth beyond the required family level where supported
      (subfamilies, genera and species), using the same growing panels.
- [ ] Prioritize the next marker release among 12S, prokaryotic 16S, animal
      mitochondrial 16S, and 28S based on reference availability and user need.
- [ ] Define marker-specific coordinate references, catalog schemas, reference
      policies, and acceptance gates before enabling each new marker.
- [ ] Add empirical evidence only as a separate, source-linked layer with exact
      experimental denominators and scope limitations.

### P4 — Public deployment and project maturity

- [ ] Verify the Posit Connect deployment against promoted R2 manifests.
- [ ] Confirm the GitHub Pages launch URL and citation instructions.
- [ ] Add versioned release notes for every promoted marker release.
- [x] Establish contribution, review, and data-provenance expectations for new
      primers and reference panels.
- [ ] Invite the first outside collaborators and appoint a second trusted
      maintainer before enforcing branch rules for administrators.

## Definition of the next milestone

The next scientific milestone is monthly retained-reference growth with
general primer scores for every supported order, family scores and coverage
gaps within every order, continued BeePrime-specific bee acquisition, a
separate aquatic arthropod option, and regional comparisons computed within
order/family from the appropriate reference. It requires
annotated binding sites, reviewed taxa, explicit reference-panel and scorable
denominators, traceable sources, and two consecutive clean-runner releases.
Regional locality and a documented dietary use must remain distinct from
species range, endemicity, and PCR performance.
