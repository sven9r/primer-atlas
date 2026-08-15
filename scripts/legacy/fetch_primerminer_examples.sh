#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target_dir="$project_root/data/reference_alignments"
work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT

git clone --depth 1 --branch v0.22 \
  https://github.com/VascoElbrecht/PrimerMiner.git \
  "$work_dir/PrimerMiner"

mkdir -p "$target_dir"
source_dir="$work_dir/PrimerMiner/Sample_Data/1 COI alignments (unprocessed)"

for source_path in "$source_dir"/*_sub.fasta; do
  cp "$source_path" "$target_dir/$(basename "$source_path")"
done

cp "$work_dir/PrimerMiner/PrimerMiner/inst/Position_v1.csv" \
  "$project_root/data/Position_v1.csv"
cp "$work_dir/PrimerMiner/PrimerMiner/inst/Type_v1.csv" \
  "$project_root/data/Type_v1.csv"

printf 'Fetched PrimerMiner v0.22 examples into %s\n' "$target_dir"
