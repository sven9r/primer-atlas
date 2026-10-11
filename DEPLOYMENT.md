# Deployment

## Posit Connect Cloud

Connect this public GitHub repository in Connect Cloud and deploy `app.R` from
`main`. Set `ATLAS_MANIFEST_BASE_URL` to the public R2 URL ending in `/releases`.
Runtime writes are caches only and may disappear when a worker stops.
The app checks this hosted pointer when a session opens, every ten minutes,
and through **Check for COI updates**. It verifies a complete candidate before
swapping that session's reference; a failed check preserves the current cohort.
Without the environment variable, it can only inspect the bundled pointer.
The Connect bundle includes the unchanged pinned PrimerMiner evaluator and
tables, so custom COI scoring does not require the installed build-only package.

## Cloudflare R2

Create a bucket and public custom-domain or `r2.dev` read URL. Add GitHub secrets
`R2_ACCOUNT_ID`, `R2_BUCKET`, `R2_ACCESS_KEY_ID`, and `R2_SECRET_ACCESS_KEY`, plus
repository variable `ATLAS_DATA_BASE_URL`. The monthly workflow uploads a staged
immutable release and replaces the small marker `latest.json` only after gates
pass.

Confirm that the five names exist without printing their values:

```bash
gh secret list --app actions
gh variable list
```

After configuration, start a controlled release and follow it to completion:

```bash
gh workflow run monthly-release.yml -f marker=all -f promote=false
gh run watch
```

Do not close a release-failure issue until both marker jobs pass, their
`latest.json` pointers resolve from the public base URL, and the application can
load the promoted manifests. The workflow has an early configuration gate so a
missing secret or variable is reported before the reference builds consume
runner time.

The schedule runs both markers on the first day of each month at 06:00 UTC.
Manual dispatch defaults to ITS_FUNGAL; select `COI` or `all` explicitly when
testing the growing COI panel. Each marker's publisher is serialized. COI
promotion also checks the currently hosted latest panel for exact retained
sequence/header hashes and scientific source archives before replacing its
pointer. A failed or stale current-pointer read prevents promotion.

The legacy public COI release lacks recoverable full raw state. Its first
activation must therefore start from an explicitly verified full-state staging
release, after the app implementation is reviewed and merged. Use
`coi_previous_release_id` to select that immutable release and `promote=true`
to build, validate, host and activate its next append-only version. Ordinary
scheduled runs then restore the complete state selected by `latest.json`.
See the dated [COI validation record](docs/2026-10-11-coi-growth.md) for actual
staging release identities and proof limits. The temporary frozen seed is only
for controlled recovery tests; it is not the recurring storage backend.

## GitHub Pages

Set Pages source to GitHub Actions. The documentation workflow publishes `docs/`.
Update the Launch Atlas URL in `docs/index.html` after the first Connect deploy.
