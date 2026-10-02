# Primer Atlas project status

**Status date:** 2026-10-02

**Current phase:** Public Shiny hosting and both monthly marker releases are operational.
The 2026.10.02 interface/catalog update is published and verified at commit
`648dd55`. The next required scientific milestone is a growing sequence
reference with order and family scores across every supported order, including
continued BeePrime-specific bee acquisition, a separate aquatic arthropod
option, and family-level regional comparisons.
Current release evidence and the remaining backlog are in
[`docs/2026-10-02-release.md`](docs/2026-10-02-release.md).
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
- **Fungal ITS — pilot:** attributed offline UNITE primer snapshot, 121
  sequence-valid primer records, eight documented pair combinations, and an
  18S–ITS1–5.8S–ITS2–28S coordinate map.
- **18S — pilot:** PR2-primer 2.1.1 snapshot with 321 individual primers and 123
  documented sets; 93 sets have drawable conventional geometry.
- **Planned registry coverage:** vertebrate mitochondrial 12S,
  bacterial/archaeal 16S, animal mitochondrial 16S, and 28S LSU have explicit
  marker identities and target-routing documentation but are not production
  releases.

### Reference evidence and safeguards

- Established full-reference and study-specific panel policies so Gurten bee
  centroids are not silently reused as general COI evidence.
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
- Cloud runtime excludes build-only PrimerMiner and readxl. A new custom pair
  can be mapped, but its interactive COI reference scoring needs PrimerMiner
  and remains unavailable on Connect Cloud.
- Research foundations are versioned independently of public evidence. The
  regional COI audit and 48-accession ITS scoring outputs are not loaded into
  the app or promoted as broad performance results.

## What still needs to be done

### P0 — Growing references and family evaluation across all orders

- [ ] Persist and restore complete retained sequence files with their ledger
      in immutable releases. Verify two successive builds on clean runners;
      the current workflow restores only the ledger.
- [ ] Rebuild or extend alignments when retained sequence input changes.
- [ ] Attach accession-level order/family taxonomy to the growing panels and
      publish coverage gaps plus primer summaries for every supported order
      and its families, retaining unresolved assignments and denominators.
- [ ] Publish and verify the new BeePrime publication/expanded alignment
      choice. The local build added 30 NCBI-classified bees to the 590-accession
      publication baseline; monthly clean-runner retention and the public app
      still need verification. Preserve the publication panel and never route
      the new bee records to general primers.
- [ ] Route general arthropod primers to their broad arthropod reference,
      covering all supported orders and families, with explicit panel identity.
- [ ] Add a separate aquatic arthropod option with freshwater/marine
      refinement, matching primer uses, and an ecology-supported reference
      subset at family/finer-taxon resolution. Mixed-habitat orders cannot
      establish aquatic membership for every sequence.
- [ ] Expose family evaluation for every curated pair/order in the app. Current
      full-order score artifacts have no family taxonomy, while current
      lineage family views are restricted to BeePrime study evidence.
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
- [ ] Verify release recovery retains both accession identities and sequences;
      ledger persistence alone does not demonstrate sequence-panel growth.

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
