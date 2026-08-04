#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
GLOBAL_ENV="${GLOBAL_ENV:-$ROOT/config/global.env}"
SPIKE_ENV="${SPIKE_ENV:-}"

fail=0
check_file() {
  if [[ ! -s "$1" ]]; then
    echo "[ERROR] Missing or empty: $1" >&2
    fail=1
  else
    echo "[OK] $1"
  fi
}

required=(
  README.md
  config/global.env.example
  config/dataset.env.example
  spikes/spikein.env.example
  spikes/spike_panel.tsv
  taxonomy/run_profiling.sh
  taxonomy/run_profiling_array.sbatch
  workflows/metaphlan4/profile.sbatch
  workflows/metaphlan4/postprocess_local.sh
  workflows/kraken2_bracken/classify_bracken.sbatch
  workflows/kraken2_bracken/postprocess_local.sh
  spikes/scripts/spikein/spikein_prepare_pools.sh
  spikes/scripts/spikein/spikein_run_independent.sh
  spikes/scripts/spikein/spikein_run_community.sh
  spikes/scripts/spikein/stable_seed.py
  datasets/yachida/build_manifest.py
  datasets/yachida/build_pilot_design.py
  datasets/yachida/stream_sample.py
  scripts/select_samples_deterministically.py
  scripts/assign_processing_batches.py
  spikes/scripts/spikein/organize_profile_tables_by_spike.py
  spikes/scripts/spikein/merge_profile_tables.py
)
for path in "${required[@]}"; do
  check_file "$ROOT/$path"
done

check_file "$GLOBAL_ENV"
if [[ -n "$SPIKE_ENV" ]]; then
  check_file "$SPIKE_ENV"
fi

if ((fail)); then
  echo "[FAIL] Required files are incomplete." >&2
  exit 1
fi

for file in \
  "$ROOT"/taxonomy/*.sh "$ROOT"/taxonomy/*.sbatch \
  "$ROOT"/workflows/metaphlan4/*.sh "$ROOT"/workflows/metaphlan4/*.sbatch \
  "$ROOT"/workflows/kraken2_bracken/*.sh "$ROOT"/workflows/kraken2_bracken/*.sbatch \
  "$ROOT"/spikes/scripts/spikein/*.sh "$ROOT"/spikes/scripts/spikein/*.sbatch
do
  bash -n "$file"
done
echo "[OK] Shell syntax"

for file in \
  "$ROOT"/spikes/scripts/spikein/*.py \
  "$ROOT"/workflows/*/*.py "$ROOT"/workflows/utils/*.py \
  "$ROOT"/datasets/yachida/*.py \
  "$ROOT"/scripts/select_samples_deterministically.py \
  "$ROOT"/scripts/assign_processing_batches.py
do
  python3 -m py_compile "$file"
done
echo "[OK] Python syntax"

expected_seed=873630871
actual_seed="$(python3 "$ROOT/spikes/scripts/spikein/stable_seed.py" \
  --base 13 --namespace spike-independent-v1 DRR127476 Fnuc 0.0001)"
[[ "$actual_seed" == "$expected_seed" ]] || {
  echo "[ERROR] stable-seed-v1 contract changed: expected=$expected_seed observed=$actual_seed" >&2
  exit 1
}
echo "[OK] Stable seed contract"

# shellcheck source=/dev/null
source "$GLOBAL_ENV"
: "${SIF:?Set SIF in $GLOBAL_ENV}"
: "${K2_DB:?Set K2_DB in $GLOBAL_ENV}"
: "${MPA_DB:?Set MPA_DB in $GLOBAL_ENV}"

[[ -s "$SIF" ]] || { echo "[ERROR] Container image missing/empty: $SIF" >&2; exit 1; }
[[ -d "$K2_DB" ]] || { echo "[ERROR] Kraken2 database missing: $K2_DB" >&2; exit 1; }
[[ -d "$MPA_DB" ]] || { echo "[ERROR] MetaPhlAn database missing: $MPA_DB" >&2; exit 1; }
[[ -s "$K2_DB/hash.k2d" && -s "$K2_DB/opts.k2d" ]] || {
  echo "[ERROR] Kraken2 database lacks hash.k2d or opts.k2d: $K2_DB" >&2
  exit 1
}
find "$K2_DB" -maxdepth 1 -type f -name '*100*mers.kmer_distrib' -print -quit |
  grep -q . || {
    echo "[ERROR] No Bracken 100-nt distribution found under: $K2_DB" >&2
    exit 1
  }

if command -v apptainer >/dev/null 2>&1; then
  CTR=apptainer
elif command -v singularity >/dev/null 2>&1; then
  CTR=singularity
else
  echo "[ERROR] Apptainer/Singularity is unavailable." >&2
  exit 1
fi

"$CTR" exec "$SIF" micromamba run -n taxonomic_tools bash -c '
  set -euo pipefail
  command -v art_illumina
  command -v seqtk
  command -v fastp
  command -v kraken2
  command -v bracken
  command -v metaphlan
' >/dev/null
echo "[OK] Required container executables"

if [[ -n "$SPIKE_ENV" ]]; then
  # shellcheck source=/dev/null
  source "$SPIKE_ENV"
  : "${WORK:?Set WORK in $SPIKE_ENV}"
  : "${SAMPLES_TSV:?Set SAMPLES_TSV in $SPIKE_ENV}"
  : "${SPIKE_PANEL:?Set SPIKE_PANEL in $SPIKE_ENV}"
  check_file "$SAMPLES_TSV"
  check_file "$SPIKE_PANEL"
  ((fail == 0)) || exit 1
  python3 - "$SAMPLES_TSV" "$SPIKE_PANEL" <<'PY'
import csv
import sys
from pathlib import Path

samples_path = Path(sys.argv[1])
panel_path = Path(sys.argv[2])
expected = ["Fnuc", "Pmic", "Pana", "Psto", "Dpne", "Csym", "Hhat", "Bfrag", "Porp", "Pint"]

with panel_path.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
required = {"label", "taxon_name", "assembly", "fasta", "weight"}
missing = required.difference(rows[0] if rows else {})
if missing:
    raise SystemExit(f"[ERROR] Spike panel lacks columns: {sorted(missing)}")
labels = [r["label"] for r in rows]
if labels != expected:
    raise SystemExit(f"[ERROR] Spike panel labels/order differ from publication design: {labels}")
for row in rows:
    if not row["assembly"].strip():
        raise SystemExit(f"[ERROR] Missing assembly for {row['label']}")
    if not Path(row["fasta"]).is_file():
        raise SystemExit(f"[ERROR] Missing FASTA for {row['label']}: {row['fasta']}")

with samples_path.open(newline="") as fh:
    samples = csv.DictReader(fh, delimiter="\t")
    missing = {"sample_id", "fastq1", "fastq2"}.difference(samples.fieldnames or [])
    if missing:
        raise SystemExit(f"[ERROR] Sample manifest lacks columns: {sorted(missing)}")
    n = 0
    for row in samples:
        n += 1
        for col in ("fastq1", "fastq2"):
            if not Path(row[col]).is_file():
                raise SystemExit(f"[ERROR] Missing {col} for {row['sample_id']}: {row[col]}")
if n == 0:
    raise SystemExit("[ERROR] Sample manifest has no data rows")
print(f"[OK] Publication spike panel ({len(rows)} taxa) and sample manifest ({n} rows)")
PY
fi

echo "[PASS] Upstream workflow preflight passed."
