#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
for tool in python3 mafft Rscript; do
  command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 1; }
done

python3 scripts/build_fungal_its_pilot.py
mafft --quiet --addfragments data/derived/its_fungal_pilot/panel.fasta --keeplength \
  data/reference/its_fungal_reference_FN812768.2.fasta \
  > data/derived/its_fungal_pilot/reference_aligned.fasta
Rscript scripts/score_fungal_its_pilot.R

python3 - <<'PY'
import hashlib
import json
import pathlib
import subprocess
out = pathlib.Path('data/derived/its_fungal_pilot')
items = ('panel.fasta', 'panel_manifest.csv', 'reference_aligned.fasta',
         'site_scores.csv', 'pair_scores.csv', 'site_summary.csv',
         'pair_summary.csv', 'clade_pair_summary.csv')
manifest = {
    'source_doi': '10.15156/BIO/3301229',
    'mafft_version': subprocess.run(['mafft', '--version'], capture_output=True,
                                    text=True, check=True).stderr.strip(),
    'primerminer_version': subprocess.run(
        ['Rscript', '-e', 'cat(as.character(packageVersion("PrimerMiner")))'],
        capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1],
    'artifacts_sha256': {name: hashlib.sha256((out / name).read_bytes()).hexdigest()
                         for name in items},
}
(out / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('Local pilot outputs:', out)
PY
