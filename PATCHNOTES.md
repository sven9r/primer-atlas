# Primer Atlas patchnotes

This is the human-readable change record. Entries are newest first and use the
same `id` as their counterparts in `patchnotes.json`. A patchnote distinguishes
repository changes from external deployment actions so that prepared code is
not mistaken for a completed release.

## 2026-09-23 — Connect Cloud readxl dependency repair

**ID:** `2026-09-23-connect-readxl-runtime-repair`

**Status:** Ready to republish from Posit Connect Cloud; public deployment is
not claimed until the hosted URL responds successfully.

- The next Connect Cloud build reached `readxl` but could not install it because
  its `cellranger` dependency was absent from the app manifest.
- `readxl` is used only by monthly data-build scripts. The Connect manifest
  generator now removes it from the Shiny runtime while the release workflow
  continues to use the project lockfile.
- The bundle regression check now requires the app's runtime packages and
  rejects build-only `readxl` and `PrimerMiner` dependencies.

## 2026-09-23 — Connect Cloud runtime dependency separation

**ID:** `2026-09-23-connect-runtime-dependency-separation`

**Status:** Ready to republish from Posit Connect Cloud; public deployment is
not claimed until the hosted URL responds successfully.

- The source-coordinate repair did not resolve Connect Cloud's refusal to
  install `PrimerMiner` from GitHub.
- `PrimerMiner` is not needed to start or browse the Atlas: its published
  evidence is precomputed. It is now excluded from the cloud runtime manifest,
  along with its build-only source folder.
- The monthly build workflow installs the vendored package explicitly before
  operations that require custom-primer scoring or new data generation.
- The regenerated Connect manifest has 88 standard dependencies, 111 files,
  and a 11.5 MB bundle; it no longer contains `PrimerMiner`.

## 2026-09-23 — Connect Cloud deployment manifest repair

**ID:** `2026-09-23-connect-manifest-primerminer`

**Status:** Ready to republish from Posit Connect Cloud; public deployment is
not claimed until the hosted URL responds successfully.

- Diagnosed the Connect Cloud dependency failure: its generated `manifest.json`
  named the vendored `PrimerMiner` package but omitted the GitHub repository
  URL, commit, and `vendor/PrimerMiner` subdirectory required to download it.
- Updated the manifest generator to preserve those source coordinates, then
  regenerated `manifest.json` (165 deployable files; 13 MB).
- Confirmed the app starts locally after the manifest repair.

## 2026-09-15 — COI retry dependency repair

**ID:** `2026-09-15-coi-xml2-retry`

**Status:** ITS_FUNGAL is promoted and public. COI completed its reference
build, then stopped before release validation because its final taxonomy import
requires `xml2`, which was not installed on the clean GitHub runner.

- [Run 35005218773](https://github.com/sven9r/primer-atlas/actions/runs/35005218773)
  proved the R2 configuration and publication path: ITS_FUNGAL built, passed
  its hard gates, and atomically promoted its public `latest.json` pointer.
- COI completed the expanded reference, alignment, geography, and scoring
  stages, then failed only at `scripts/fetch_gurten2026_bee_taxonomy.R` with
  `there is no package called 'xml2'`.
- `scripts/install_primerminer.R` now installs `xml2` explicitly with the
  other marker-build dependencies. The affected taxonomy script was executed
  locally against its cached inputs successfully.
- The workflow now also removes CR/LF characters before its early COI R2 state
  read, so the append-only ledger can be recovered correctly on the retry.

COI has not been promoted; its previous public pointer remains unchanged. The
next controlled retry must pass COI validation and public-artifact checks before
the first production release is declared complete.

## 2026-09-15 — First live R2 release attempt

**ID:** `2026-09-15-r2-credential-normalization`

**Status:** Live configuration reached the R2 publisher; promotion is still
pending a retry with line-ending-safe credential handling.

- Created the public `primer-atlas` Cloudflare R2 bucket and configured its
  public base URL as the GitHub `ATLAS_DATA_BASE_URL` variable.
- Configured the four required GitHub Actions secret names and dispatched
  [monthly release run 35003502393](https://github.com/sven9r/primer-atlas/actions/runs/35003502393).
- Both matrix jobs passed the release-configuration preflight. ITS_FUNGAL also
  built and passed its hard release gates, proving that the repository and R2
  configuration are connected.
- R2 publication failed before any pointer promotion because a browser-copied
  S3 access key included a line ending, making AWS Signature V4 generate an
  invalid HTTP Authorization header.
- `scripts/publish_r2.sh` now strips CR/LF characters from the R2 account,
  bucket, access-key, and secret-key inputs before invoking the S3 client. This
  makes GitHub-secret copy/paste line endings harmless on the next controlled
  retry.

No marker release is claimed as promoted until the retried jobs pass and both
public `latest.json` pointers resolve.

## 2026-09-14 — Collaboration foundation

**ID:** `2026-09-14-collaboration-foundation`

**Status:** Integrated and active; collaborator invitations await GitHub
usernames.

- Expanded `CONTRIBUTING.md` with local setup, contribution workflow, complete
  test commands, scientific evidence requirements, and review expectations.
- Added a code of conduct, default review ownership, and a pull-request evidence
  checklist.
- Added structured bug and reference-panel proposal forms alongside the existing
  primer-proposal form.
- Added an issue-form router directing early scientific ideas to GitHub
  Discussions.
- Added contributor entry points to the README and corrected stale project-status
  language now that the ITS navigator is on `main`.

The files were merged through
[pull request 4](https://github.com/sven9r/primer-atlas/pull/4). GitHub
Discussions, a public description, the Pages homepage, eight scientific topics,
and nine workflow labels are active. `main` now requires an up-to-date green
`test` check, a pull request, one approval, code-owner review, and resolved
conversations; force-push and deletion are disabled. Administrator enforcement
remains off until a second trusted maintainer exists, avoiding an owner lockout.

The test workflow now runs branch work through the pull-request event and direct
`main` updates through the push event. This removes the duplicate full-suite
runs previously triggered by every branch push plus its pull request.

## 2026-09-14 — Recovery branch separated, pushed, and verified

**ID:** `2026-09-14-recovery-branch-verified`

**Status:** Integrated and verified on `main`; external R2 configuration and a
successful marker release remain outstanding.

The previous working-tree changes were separated into two reviewable commits on
`codex/release-recovery`:

- `c77b648` — primer-first ITS navigator and its regression coverage.
- `6901b36` — monthly-release recovery, dependency lock, deployment guidance,
  project status, and patchnotes.

The branch was pushed to GitHub. The resulting
[Test atlas run 34912008403](https://github.com/sven9r/primer-atlas/actions/runs/34912008403)
restored the project environment on a clean Ubuntu runner and passed the full
regression suite and deployable-bundle-size gate in 4 minutes 58 seconds. This
provides remote confirmation that the updated `renv.lock`, including `readxl`
1.5.0, is usable by CI.

[Pull request 3](https://github.com/sven9r/primer-atlas/pull/3) was then merged
into `main` as `764f83d`. The post-merge
[Test atlas run 34912852240](https://github.com/sven9r/primer-atlas/actions/runs/34912852240)
again restored the clean environment and passed the full regression and bundle
checks in 4 minutes 45 seconds.

No monthly marker release was dispatched because GitHub still has none of the
four required R2 secrets or the `ATLAS_DATA_BASE_URL` variable. Running it in
that state would intentionally fail the new preflight and open more failure
issues without testing publication.

## 2026-09-01 — First monthly-release incident and recovery preparation

**ID:** `2026-09-01-monthly-release-recovery`

**Status:** Partially resolved; repository fixes are prepared, external R2
configuration and a successful rerun remain outstanding.

The first scheduled **Monthly marker release** run
([run 33502385055](https://github.com/sven9r/primer-atlas/actions/runs/33502385055))
failed in both matrix jobs for independent reasons:

- **COI:** the build reached `scripts/fetch_gurten2026_bee_taxonomy.R` and
  stopped because `readxl` was used by the build scripts but absent from
  `renv.lock`. A globally installed local copy had masked the incomplete
  project environment.
- **ITS_FUNGAL:** the marker artifacts built and passed the release gates, but
  R2 publication stopped because `R2_ACCOUNT_ID` was empty. Repository
  inspection confirmed that none of the four required Actions secrets or the
  `ATLAS_DATA_BASE_URL` repository variable had been configured.

No marker was promoted to R2. The workflow retained its safety property: it did
not update a `latest.json` pointer after either failure. QA artifacts were
uploaded, and the workflow opened
[issue 1](https://github.com/sven9r/primer-atlas/issues/1) for ITS_FUNGAL and
[issue 2](https://github.com/sven9r/primer-atlas/issues/2) for COI.

### Repository changes prepared

- Recorded `readxl` 1.5.0 in `renv.lock` so clean CI environments install the
  dependency used by the expanded COI reference build.
- Added an early workflow preflight that reports missing release-configuration
  names before installing tools or building marker references.
- Expanded `DEPLOYMENT.md` with configuration checks, a controlled manual rerun,
  and post-release acceptance criteria.
- Added this human-readable patchnote, the synchronized `patchnotes.json`, and
  `PROJECT_STATUS.md`.

### Still required outside the repository

1. Create or select the Cloudflare R2 bucket and its public read URL.
2. Add GitHub Actions secrets `R2_ACCOUNT_ID`, `R2_BUCKET`,
   `R2_ACCESS_KEY_ID`, and `R2_SECRET_ACCESS_KEY`.
3. Add repository variable `ATLAS_DATA_BASE_URL`.
4. Commit and push the recovery changes, manually run the monthly workflow, and
   verify both public `latest.json` pointers and application loading.
5. Close issues 1 and 2 only after that end-to-end verification succeeds.

### Verification completed locally

- Confirmed that both expanded-reference scripts declare `readxl`.
- Confirmed that `renv.lock` now contains `readxl` 1.5.0.
- Parsed `patchnotes.json` and the monthly workflow successfully.
- Ran the complete regression loop from `.github/workflows/test.yml`; all 15
  test scripts passed in the current working tree.
- Ran `git diff --check` successfully.
- Confirmed that the human and machine patchnotes describe the same incident
  and remaining external actions.

This patchnote does not claim that R2 has been configured, that the workflow has
been rerun, or that a production marker release has completed.
