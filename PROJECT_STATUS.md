# Primer Atlas project status

**Status date:** 2026-09-15

**Current phase:** Scientific and interface foundation implemented; first
production release automation and public Shiny hosting are being brought online.

## Project goal

Primer Atlas is an open, reproducible Shiny application for investigating where
metabarcoding primers bind and why a lineage–primer combination may be
vulnerable to mismatch or missing-reference bias. It is intended to connect
primer sequences, published claims, coordinate geometry, reference-panel
composition, mismatch evidence, and exact source records without turning any
single evidence layer into an unsupported amplification-probability claim.

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

- The repository-side application, build scripts, validation gates, tests, and
  primer-first ITS navigator exist on `main`. ITS_FUNGAL has a public release;
  COI remains blocked by a release gate.
- The first scheduled monthly workflow ran on 2026-09-01. COI exposed a missing
  locked `readxl` dependency. ITS_FUNGAL built and validated but could not
  publish because the R2 repository configuration was absent.
- The Cloudflare R2 bucket, four credential names, and public data URL are
  configured. In live run 35005218773, ITS_FUNGAL built, passed its hard gates,
  and was promoted publicly. The subsequent run 35009155171 built COI and
  passed its initial validation, but the COI hard gate failed on a
  Lepidoptera study-group label assertion. Its public pointer was not promoted.
- The public Posit Connect Cloud deployment was initiated on 2026-09-23. Its
  first build stopped because the generated manifest omitted source coordinates
  for the vendored `PrimerMiner` package; Connect Cloud still refused the
  package download. The live app does not need `PrimerMiner` to start or browse
  its precomputed evidence, so the package is now excluded from the cloud
  runtime and installed explicitly only by the marker-build workflow.
- A subsequent Connect Cloud publish reached `readxl` and failed because its
  `cellranger` dependency was missing from the app manifest. `readxl` is used
  by monthly data-build scripts, so the cloud manifest now excludes it while
  the release lockfile retains it. After merge `42165ea`, the public app loaded
  without authentication at
  https://01a0cdbb-4448-22b4-500c-b5329e3b1904.share.connect.posit.cloud/;
  the COI map populated and the URL returned the `Primer Atlas` page title.
- All 15 scripts in the GitHub Actions regression loop pass in the current
  working tree, including release-layer, marker-runtime, reference-policy,
  geography, map, and sequence-drill-down checks.
- The ITS navigator and release recovery were separated into scoped commits,
  verified on `codex/release-recovery`, and merged through pull request 3 as
  `764f83d`. Post-merge GitHub Test atlas run 34912852240 passed the clean
  environment restore, full regression suite, and bundle-size gate on `main`.
- Therefore, **the production monthly release path is not yet operational**.
  Repository readiness must not be reported as a completed R2 deployment.

## What still needs to be done

### P0 — Complete the first production release

- [x] Diagnose both failures from workflow run 33502385055.
- [x] Record `readxl` in `renv.lock`.
- [x] Add a fail-fast release-configuration preflight.
- [x] Configure the four R2 GitHub Actions secrets.
- [x] Configure the `ATLAS_DATA_BASE_URL` repository variable.
- [x] Separate the ITS navigator and release recovery into scoped commits and
      push `codex/release-recovery`.
- [x] Pass the full GitHub Test atlas workflow on the recovery branch.
- [x] Integrate the verified recovery branch into `main` through pull request 3.
- [ ] Merge the `xml2` marker-build dependency repair, then manually rerun the
      monthly workflow and require the COI job to pass (ITS_FUNGAL is already
      promoted).
- [ ] Verify the COI and ITS_FUNGAL public `latest.json` pointers, checksums,
      immutable artifact URLs, and application loading.
- [ ] Close failure issues 1 and 2 with links to the successful evidence.

### P1 — Stabilize release operations

- [ ] Add a low-cost scheduled or dispatchable smoke test for configuration and
      public manifest reachability.
- [ ] Decide whether repeated matrix failures should update existing marker
      issues instead of opening duplicates each month.
- [ ] Test the pinned-summary fallback during an intentional R2 outage.
- [ ] Document release ownership, credential rotation, rollback, and recovery
      from a partially uploaded staging directory.
- [ ] Confirm that the monthly append-only accession ledger persists correctly
      across two consecutive successful COI releases.

### P2 — Finish active pilots

- [x] Merge and verify the ITS primer-first navigator in a clean browser session.
- [ ] Add 18S to the release matrix only after its artifact contract and hard
      release gates are defined.
- [ ] Expand fungal ITS beyond the single pilot reference while continuing to
      label reference scope and missing lineage evidence explicitly.
- [ ] Publish a clear per-marker readiness table in the public documentation.

### P3 — Expand scientific coverage without weakening evidence standards

- [ ] Complete append-only full-order COI production panels and taxonomy
      partitions for general primer pairs.
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

The recovery milestone is complete only when a pushed revision produces green
COI and ITS_FUNGAL jobs, both public manifests and their checksums resolve, the
deployed application loads those releases, and the two incident issues contain
links to that evidence. Until then, the release is **prepared**, not
**completed**.
