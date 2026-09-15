# Primer Atlas patchnotes

This is the human-readable change record. Entries are newest first and use the
same `id` as their counterparts in `patchnotes.json`. A patchnote distinguishes
repository changes from external deployment actions so that prepared code is
not mistaken for a completed release.

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
