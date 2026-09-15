# Deployment

## Posit Connect Cloud

Connect this public GitHub repository in Connect Cloud and deploy `app.R` from
`main`. Set `ATLAS_MANIFEST_BASE_URL` to the public R2 URL ending in `/releases`.
Runtime writes are caches only and may disappear when a worker stops.

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
gh workflow run monthly-release.yml
gh run watch
```

Do not close a release-failure issue until both marker jobs pass, their
`latest.json` pointers resolve from the public base URL, and the application can
load the promoted manifests. The workflow has an early configuration gate so a
missing secret or variable is reported before the reference builds consume
runner time.

## GitHub Pages

Set Pages source to GitHub Actions. The documentation workflow publishes `docs/`.
Update the Launch Atlas URL in `docs/index.html` after the first Connect deploy.
