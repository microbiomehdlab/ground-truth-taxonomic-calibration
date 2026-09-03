#!/usr/bin/env bash
# Read-only gate for a frozen FengQ or ZellerG production configuration.
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --env cohort.env --manifest production_manifest.tsv" >&2
}

env_file=""
manifest=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) env_file="$2"; shift 2 ;;
    --manifest) manifest="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -s "$env_file" && -s "$manifest" ]] || { usage; exit 2; }
env_file="$(cd "$(dirname "$env_file")" && pwd -P)/$(basename "$env_file")"
manifest="$(cd "$(dirname "$manifest")" && pwd -P)/$(basename "$manifest")"

# Site configuration is trusted shell input and is intentionally kept outside Git.
source "$env_file"

required_variables=(
  METASHOTGUNPREP_ROOT METASHOTGUNPREP_COMMIT UPSTREAM_SIF HOST_INDEX
  PREPROCESS_THREADS K2_THREADS MPA_THREADS PROFILE_CONCURRENCY
  KRAKEN2_DB METAPHLAN_DB METAPHLAN_INDEX READ_LEN BRACKEN_THRESHOLD
  CRC_SCRATCH_ROOT CRC_STATE_DIR PERSISTENT_QC_ROOT PERSISTENT_RESULTS_ROOT
  INDEPENDENT_MANIFEST SPIKE_ENV YACHIDA_POOLS_DIR
)
for variable in "${required_variables[@]}"; do
  [[ -n "${!variable:-}" ]] || { echo "[ERROR] Missing configuration: $variable" >&2; exit 1; }
done

[[ "$METASHOTGUNPREP_COMMIT" == a717a105d56934e205c21ef59a316b8613b6d1c1 ]] || {
  echo "[ERROR] Unexpected MetaShotgunPrep commit: $METASHOTGUNPREP_COMMIT" >&2; exit 1;
}
[[ "$(git -C "$METASHOTGUNPREP_ROOT" rev-parse HEAD)" == "$METASHOTGUNPREP_COMMIT" ]] || {
  echo "[ERROR] MetaShotgunPrep checkout does not match configured commit" >&2; exit 1;
}
[[ -z "$(git -C "$METASHOTGUNPREP_ROOT" status --porcelain --untracked-files=no)" ]] || {
  echo "[ERROR] MetaShotgunPrep checkout has tracked modifications" >&2; exit 1;
}
[[ "$PREPROCESS_THREADS" == 8 && "$K2_THREADS" == 8 && "$MPA_THREADS" == 8 ]] || {
  echo "[ERROR] Frozen production thread settings must all equal 8" >&2; exit 1;
}
[[ "$PROFILE_CONCURRENCY" == 1 && "$READ_LEN" == 100 && "$BRACKEN_THRESHOLD" == 10 ]] || {
  echo "[ERROR] Frozen profiling parameters do not match the production contract" >&2; exit 1;
}

for file in "$UPSTREAM_SIF" "$INDEPENDENT_MANIFEST" "$SPIKE_ENV" \
  "$YACHIDA_POOLS_DIR/pool_pair_counts.tsv" "$YACHIDA_POOLS_DIR/pool_pair_counts.tsv.sha256"; do
  [[ -s "$file" ]] || { echo "[ERROR] Missing or empty required file: $file" >&2; exit 1; }
done
for file in hash.k2d opts.k2d taxo.k2d database100mers.kmer_distrib; do
  [[ -s "$KRAKEN2_DB/$file" ]] || { echo "[ERROR] Missing Kraken2/Bracken asset: $file" >&2; exit 1; }
done
[[ -s "$METAPHLAN_DB/$METAPHLAN_INDEX.pkl" ]] || {
  echo "[ERROR] Missing MetaPhlAn database pickle" >&2; exit 1;
}
for suffix in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
  [[ -s "${HOST_INDEX}${suffix}" || -s "${HOST_INDEX}${suffix}l" ]] || {
    echo "[ERROR] Missing Bowtie2 index component: ${HOST_INDEX}${suffix}[l]" >&2; exit 1;
  }
done
(cd "$YACHIDA_POOLS_DIR" && sha256sum -c pool_pair_counts.tsv.sha256 >/dev/null)

python3 - "$manifest" "$INDEPENDENT_MANIFEST" <<'PY'
import csv, pathlib, sys

production, independent = map(pathlib.Path, sys.argv[1:])
required = {
    "sample_id", "Study", "Target_Condition", "run_count", "run_accessions",
    "fastq1_urls", "fastq2_urls", "fastq1_md5s", "fastq2_md5s",
    "fastq1_bytes", "fastq2_bytes",
}

def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise SystemExit(f"[ERROR] {path}: missing columns: {','.join(missing)}")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"[ERROR] Empty manifest: {path}")
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit(f"[ERROR] Duplicate sample_id in {path}")
    return rows

prod = read(production)
ind = read(independent)
prod_ids = {row["sample_id"] for row in prod}
if len(ind) != 30 or not {row["sample_id"] for row in ind}.issubset(prod_ids):
    raise SystemExit("[ERROR] Independent manifest must be a 30-sample subset of production")
counts = {}
for row in ind:
    counts[row["Target_Condition"]] = counts.get(row["Target_Condition"], 0) + 1
if counts != {"Control": 10, "Adenoma": 10, "CRC": 10}:
    raise SystemExit(f"[ERROR] Independent condition counts are not 10/10/10: {counts}")
print(f"[PASS] Production manifest: {len(prod)} unique samples")
print("[PASS] Independent subset: 30 samples (10 Control, 10 Adenoma, 10 CRC)")
PY

echo "[PASS] CRC cohort production preflight completed"
