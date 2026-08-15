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

## GitHub Pages

Set Pages source to GitHub Actions. The documentation workflow publishes `docs/`.
Update the Launch Atlas URL in `docs/index.html` after the first Connect deploy.
