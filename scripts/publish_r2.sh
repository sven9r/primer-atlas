#!/usr/bin/env bash
set -euo pipefail

: "${R2_ACCOUNT_ID:?R2_ACCOUNT_ID is required}"
: "${R2_BUCKET:?R2_BUCKET is required}"
: "${AWS_ACCESS_KEY_ID:?AWS_ACCESS_KEY_ID is required}"
: "${AWS_SECRET_ACCESS_KEY:?AWS_SECRET_ACCESS_KEY is required}"
: "${MARKER_ID:?MARKER_ID is required}"
: "${RELEASE_ID:?RELEASE_ID is required}"

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

aws s3 sync "$source_dir" "$staging" --endpoint-url "$endpoint" --checksum-algorithm SHA256
aws s3 sync "$source_dir" "$release" --endpoint-url "$endpoint" --checksum-algorithm SHA256
aws s3 cp "$source_dir/manifest.json" "s3://${R2_BUCKET}/releases/${MARKER_ID}/latest.json" \
  --endpoint-url "$endpoint" --content-type application/json --cache-control no-cache

if [ "$MARKER_ID" = "COI" ] && [ -f data/provenance/coi_reference_accession_ledger.csv ]; then
  aws s3 cp data/provenance/coi_reference_accession_ledger.csv "s3://${R2_BUCKET}/state/COI/coi_reference_accession_ledger.csv" \
    --endpoint-url "$endpoint" --content-type text/csv
  aws s3 cp data/provenance/coi_reference_growth_qa.csv "s3://${R2_BUCKET}/state/COI/coi_reference_growth_qa.csv" \
    --endpoint-url "$endpoint" --content-type text/csv
fi
