#!/usr/bin/env bash
# Profile one already-preprocessed Yachida sample without deleting any input.
set -euo pipefail
IFS=$'\n\t'
umask 002

: "${YACHIDA_ENV:?set YACHIDA_ENV to the local ignored configuration}"
# shellcheck source=/dev/null
source "$YACHIDA_ENV"

: "${SAMPLE_ID:?set one frozen-pilot sample ID}"
: "${SAMPLE_WORK:?set the retained preprocessing smoke-test directory}"
: "${UPSTREAM_SIF:?set in YACHIDA_ENV}"
: "${KRAKEN2_DB:?set in YACHIDA_ENV}"
: "${METAPHLAN_DB:?set in YACHIDA_ENV}"

K2_THREADS="${K2_THREADS:-16}"
MPA_THREADS="${MPA_THREADS:-8}"
READ_LEN="${READ_LEN:-100}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"
METAPHLAN_INDEX="${METAPHLAN_INDEX:-mpa_vJan25_CHOCOPhlAnSGB_202503}"
PROFILE_ROOT="${YACHIDA_BASELINE_SMOKE_ROOT:?set the persistent baseline-smoke output root}"

[[ "$SAMPLE_ID" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "[ERROR] Unsafe sample ID: $SAMPLE_ID" >&2; exit 1; }
[[ -s "$UPSTREAM_SIF" ]] || { echo "[ERROR] Missing upstream image: $UPSTREAM_SIF" >&2; exit 1; }
command -v apptainer >/dev/null 2>&1 || { echo "[ERROR] apptainer is unavailable" >&2; exit 1; }

handoff="$SAMPLE_WORK/metashotgunprep_outputs.env"
[[ -s "$handoff" ]] || { echo "[ERROR] Missing preprocessing handoff: $handoff" >&2; exit 1; }
# shellcheck source=/dev/null
source "$handoff"
for file in "$CLEAN_R1" "$CLEAN_R2"; do
  [[ -s "$file" ]] || { echo "[ERROR] Missing cleaned mate: $file" >&2; exit 1; }
  gzip -t -- "$file"
done

for file in hash.k2d opts.k2d taxo.k2d; do
  [[ -s "$KRAKEN2_DB/$file" ]] || { echo "[ERROR] Missing Kraken2 database file: $KRAKEN2_DB/$file" >&2; exit 1; }
done
bracken_kmers="$KRAKEN2_DB/database${READ_LEN}mers.kmer_distrib"
[[ -s "$bracken_kmers" ]] || {
  echo "[ERROR] Bracken distribution for read length $READ_LEN is missing: $bracken_kmers" >&2
  exit 1
}
[[ -s "$METAPHLAN_DB/${METAPHLAN_INDEX}.pkl" ]] || {
  echo "[ERROR] Missing MetaPhlAn index metadata: $METAPHLAN_DB/${METAPHLAN_INDEX}.pkl" >&2
  exit 1
}
for suffix in 1 2 3 4 rev.1 rev.2; do
  [[ -s "$METAPHLAN_DB/${METAPHLAN_INDEX}.${suffix}.bt2" ||
     -s "$METAPHLAN_DB/${METAPHLAN_INDEX}.${suffix}.bt2l" ]] || {
    echo "[ERROR] Missing MetaPhlAn Bowtie2 component: ${METAPHLAN_INDEX}.${suffix}.bt2[l]" >&2
    exit 1
  }
done

outdir="$PROFILE_ROOT/$SAMPLE_ID"
mkdir -p "$outdir"
kraken_report="$outdir/${SAMPLE_ID}.kraken2.report"
bracken_species="$outdir/${SAMPLE_ID}.bracken.S.tsv"
metaphlan_profile="$outdir/${SAMPLE_ID}.metaphlan.tsv"
versions="$outdir/tool_versions.txt"
parameters="$outdir/profiling_parameters.tsv"
receipt="$outdir/retained_outputs.tsv"
rm -f -- "$outdir/SUCCESS"
mapout="$SAMPLE_WORK/${SAMPLE_ID}.baseline_smoke.mapout.bz2"
container_home="$SAMPLE_WORK/baseline_profile_home"
mkdir -p "$container_home/.cache" "$container_home/.config" "$container_home/.local/share"

appt=(apptainer exec --cleanenv --home "$container_home" "$UPSTREAM_SIF")

{
  printf 'tool\tversion\n'
  printf 'kraken2\t%s\n' "$("${appt[@]}" kraken2 --version 2>&1 | head -n 1)"
  printf 'bracken\t%s\n' "$("${appt[@]}" bracken -v 2>&1 | head -n 1)"
  printf 'metaphlan\t%s\n' "$("${appt[@]}" metaphlan --version 2>&1 | head -n 1)"
} > "$versions"

cat > "$parameters" <<EOF
field	value
sample_id	$SAMPLE_ID
kraken2_db	$KRAKEN2_DB
metaphlan_db	$METAPHLAN_DB
metaphlan_index	$METAPHLAN_INDEX
read_length	$READ_LEN
bracken_threshold	$BRACKEN_THRESHOLD
kraken2_threads	$K2_THREADS
metaphlan_threads	$MPA_THREADS
upstream_sif	$UPSTREAM_SIF
EOF

echo "[STEP 1/2] Kraken2 and Bracken: $SAMPLE_ID"
"${appt[@]}" kraken2 \
  --db "$KRAKEN2_DB" \
  --threads "$K2_THREADS" \
  --paired "$CLEAN_R1" "$CLEAN_R2" \
  --report "$kraken_report" \
  --output /dev/null
"${appt[@]}" bracken \
  -d "$KRAKEN2_DB" \
  -i "$kraken_report" \
  -o "$bracken_species" \
  -r "$READ_LEN" \
  -l S \
  -t "$BRACKEN_THRESHOLD"
[[ -s "$kraken_report" ]] || { echo "[ERROR] Empty Kraken2 report" >&2; exit 1; }
awk 'NR > 1 {found=1} END {exit !found}' "$bracken_species" || {
  echo "[ERROR] Bracken species output contains no data rows" >&2
  exit 1
}

echo "[STEP 2/2] MetaPhlAn 4: $SAMPLE_ID"
"${appt[@]}" metaphlan \
  "$CLEAN_R1,$CLEAN_R2" \
  --input_type fastq \
  --db_dir "$METAPHLAN_DB" \
  --index "$METAPHLAN_INDEX" \
  --nproc "$MPA_THREADS" \
  --ignore_eukaryotes \
  --ignore_archaea \
  --mapout "$mapout" \
  -o "$metaphlan_profile"
[[ -s "$metaphlan_profile" ]] || { echo "[ERROR] Empty MetaPhlAn profile" >&2; exit 1; }
grep -qE '(^|[|])s__[^|[:space:]]+' "$metaphlan_profile" || {
  echo "[ERROR] MetaPhlAn profile contains no species-level rows" >&2
  exit 1
}
rm -f -- "$mapout"

python3 - "$outdir" "$receipt" <<'PY'
import hashlib
import pathlib
import sys

root = pathlib.Path(sys.argv[1]).resolve()
receipt = pathlib.Path(sys.argv[2]).resolve()
files = sorted(path for path in root.iterdir() if path.is_file() and path != receipt and path.stat().st_size)
with receipt.open("w", encoding="utf-8") as handle:
    handle.write("path\tsha256\tbytes\n")
    for path in files:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            while block := source.read(4 * 1024 * 1024):
                digest.update(block)
        handle.write(f"{path}\t{digest.hexdigest()}\t{path.stat().st_size}\n")
print(f"[OK] Checksummed {len(files)} retained baseline outputs")
PY

touch "$outdir/SUCCESS"
echo "[PASS] Baseline profiler smoke test completed: $SAMPLE_ID"
echo "[KEEP] Cleaned reads retained: $SAMPLE_WORK"
echo "[INFO] Profiles: $outdir"
