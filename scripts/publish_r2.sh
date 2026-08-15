#!/usr/bin/env bash
set -euo pipefail

: "${R2_ACCOUNT_ID:?R2_ACCOUNT_ID is required}"
: "${R2_BUCKET:?R2_BUCKET is required}"
: "${AWS_ACCESS_KEY_ID:?AWS_ACCESS_KEY_ID is required}"
: "${AWS_SECRET_ACCESS_KEY:?AWS_SECRET_ACCESS_KEY is required}"
: "${MARKER_ID:?MARKER_ID is required}"
: "${RELEASE_ID:?RELEASE_ID is required}"

endpoint="https://${R2_ACCOUNT_ID}.r2.cloudflarestorage.com"
source_dir="data/releases/${MARKER_ID}/${RELEASE_ID}"
staging="s3://${R2_BUCKET}/staging/${MARKER_ID}/${RELEASE_ID}"
release="s3://${R2_BUCKET}/releases/${MARKER_ID}/${RELEASE_ID}"

aws s3 sync "$source_dir" "$staging" --endpoint-url "$endpoint" --checksum-algorithm SHA256
aws s3 sync "$source_dir" "$release" --endpoint-url "$endpoint" --checksum-algorithm SHA256
aws s3 cp "$source_dir/manifest.json" "s3://${R2_BUCKET}/releases/${MARKER_ID}/latest.json" \
  --endpoint-url "$endpoint" --content-type application/json --cache-control no-cache
