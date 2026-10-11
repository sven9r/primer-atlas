#!/usr/bin/env bash
set -euo pipefail

: "${R2_ACCOUNT_ID:?R2_ACCOUNT_ID is required}"
: "${R2_BUCKET:?R2_BUCKET is required}"
: "${AWS_ACCESS_KEY_ID:?AWS_ACCESS_KEY_ID is required}"
: "${AWS_SECRET_ACCESS_KEY:?AWS_SECRET_ACCESS_KEY is required}"
: "${MARKER_ID:?MARKER_ID is required}"
: "${RELEASE_ID:?RELEASE_ID is required}"
PROMOTE_RELEASE="${PROMOTE_RELEASE:-false}"
if [ "$PROMOTE_RELEASE" != true ] && [ "$PROMOTE_RELEASE" != false ]; then
  echo "PROMOTE_RELEASE must be true or false" >&2
  exit 1
fi

# Browser copy controls can put a trailing line ending in a GitHub Actions
# secret. AWS Signature V4 treats that line ending as part of the credential,
# which produces an invalid Authorization header. R2 identifiers and S3 keys
# cannot contain line endings, so remove them before constructing the client.
R2_ACCOUNT_ID="$(printf '%s' "$R2_ACCOUNT_ID" | tr -d '\r\n')"
R2_BUCKET="$(printf '%s' "$R2_BUCKET" | tr -d '\r\n')"
AWS_ACCESS_KEY_ID="$(printf '%s' "$AWS_ACCESS_KEY_ID" | tr -d '\r\n')"
AWS_SECRET_ACCESS_KEY="$(printf '%s' "$AWS_SECRET_ACCESS_KEY" | tr -d '\r\n')"

endpoint="https://${R2_ACCOUNT_ID}.r2.cloudflarestorage.com"
source_dir="data/releases/${MARKER_ID}/${RELEASE_ID}"
staging="s3://${R2_BUCKET}/staging/${MARKER_ID}/${RELEASE_ID}"
release="s3://${R2_BUCKET}/releases/${MARKER_ID}/${RELEASE_ID}"

if [ ! -f "$source_dir/manifest.json" ]; then
  echo "Validated release manifest is missing" >&2
  exit 1
fi
existing=$(aws s3api list-objects-v2 --bucket "$R2_BUCKET" \
  --prefix "releases/${MARKER_ID}/${RELEASE_ID}/" --max-keys 1 \
  --endpoint-url "$endpoint" --query KeyCount --output text)
if [ "$existing" != 0 ]; then
  echo "Refusing to overwrite an existing immutable R2 release: ${MARKER_ID}/${RELEASE_ID}" >&2
  exit 1
fi

aws s3 sync "$source_dir" "$staging" --endpoint-url "$endpoint" --checksum-algorithm SHA256
aws s3 sync "$source_dir" "$release" --endpoint-url "$endpoint" --checksum-algorithm SHA256 \
  --cache-control 'public,max-age=31536000,immutable'

# Verify the exact published bytes anonymously before updating any pointer.
# A successful authenticated upload does not establish app-readable hosting.
python3 - "$source_dir/manifest.json" <<'PY'
import hashlib, json, subprocess, sys
from urllib.parse import urlsplit
manifest_bytes = open(sys.argv[1], 'rb').read()
manifest = json.loads(manifest_bytes)
first_url = next(iter(manifest['artifacts'].values()))['url']
public_base = first_url.split('/releases/', 1)[0]
manifest_url = public_base + '/releases/' + manifest['marker_id'] + '/' + manifest['release_version'] + '/manifest.json'
published_manifest = subprocess.check_output(['curl', '-fsS', '--retry', '3', '--retry-delay', '2', '--max-time', '120', manifest_url])
if hashlib.sha256(published_manifest).digest() != hashlib.sha256(manifest_bytes).digest():
    raise SystemExit('Public immutable manifest checksum differs')
for name, artifact in manifest['artifacts'].items():
    url = artifact['url']
    parsed = urlsplit(url)
    if parsed.scheme != 'https' or parsed.username or parsed.password:
        raise SystemExit('Artifact needs a public HTTPS URL')
    payload = subprocess.check_output(['curl', '-fsS', '--retry', '3', '--retry-delay', '2', '--max-time', '120', url])
    if hashlib.sha256(payload).hexdigest() != artifact['sha256']:
        raise SystemExit('Public artifact checksum differs: ' + name)
print('Public artifact read-back SHA-256 checks passed for', len(manifest['artifacts']), 'artifacts')
PY

if [ "$PROMOTE_RELEASE" != true ]; then
  echo "Immutable release uploaded and verified; latest.json retained (PROMOTE_RELEASE=false)."
  exit 0
fi
aws s3 cp "$source_dir/manifest.json" "s3://${R2_BUCKET}/releases/${MARKER_ID}/latest.json" \
  --endpoint-url "$endpoint" --content-type application/json --cache-control no-cache

# COI raw sequences and provenance are already immutable release artifacts.
# The single latest pointer selects their complete state atomically.
