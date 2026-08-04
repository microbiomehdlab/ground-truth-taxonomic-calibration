#!/usr/bin/env bash
# Portable adapter from stream_sample.py's R1/R2 contract to MetaShotgunPrep.
set -euo pipefail
IFS=$'\n\t'
umask 002

if [[ -n "${YACHIDA_ENV:-}" ]]; then
  [[ -s "$YACHIDA_ENV" ]] || { echo "[ERROR] Missing YACHIDA_ENV: $YACHIDA_ENV" >&2; exit 1; }
  # shellcheck source=/dev/null
  source "$YACHIDA_ENV"
fi

: "${SAMPLE_ID:?provided by stream_sample.py}"
: "${RAW_R1:?provided by stream_sample.py}"
: "${RAW_R2:?provided by stream_sample.py}"
: "${SAMPLE_WORK:?provided by stream_sample.py}"
: "${STUDY:?provided by stream_sample.py}"
: "${METASHOTGUNPREP_ROOT:?set in YACHIDA_ENV}"
: "${METASHOTGUNPREP_COMMIT:?set in YACHIDA_ENV}"
: "${HOST_INDEX:?set Bowtie2 index prefix in YACHIDA_ENV}"
: "${PERSISTENT_QC_ROOT:?set in YACHIDA_ENV}"

METASHOTGUNPREP_PYTHON="${METASHOTGUNPREP_PYTHON:-python3}"
PREPROCESS_THREADS="${PREPROCESS_THREADS:-16}"
pipeline="$METASHOTGUNPREP_ROOT/pipeline_preprocess.py"
[[ -s "$pipeline" ]] || { echo "[ERROR] Missing MetaShotgunPrep entry point: $pipeline" >&2; exit 1; }
[[ -s "$RAW_R1" && -s "$RAW_R2" ]] || { echo "[ERROR] Raw mates are missing or empty" >&2; exit 1; }
[[ "$SAMPLE_ID" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "[ERROR] Unsafe sample ID: $SAMPLE_ID" >&2; exit 1; }
[[ "$STUDY" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "[ERROR] Unsafe study ID: $STUDY" >&2; exit 1; }

"$METASHOTGUNPREP_PYTHON" - "$SAMPLE_WORK" "$PERSISTENT_QC_ROOT" <<'PY'
import pathlib, sys
scratch, persistent = (pathlib.Path(value).resolve() for value in sys.argv[1:])
if persistent == scratch or scratch in persistent.parents:
    raise SystemExit("[ERROR] PERSISTENT_QC_ROOT cannot be inside disposable SAMPLE_WORK")
PY

for suffix in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
  if [[ ! -s "${HOST_INDEX}${suffix}" && ! -s "${HOST_INDEX}${suffix}l" ]]; then
    echo "[ERROR] Bowtie2 index component is missing: ${HOST_INDEX}${suffix}[l]" >&2
    exit 1
  fi
done

if ! command -v git >/dev/null 2>&1; then
  echo "[ERROR] git is required to verify the MetaShotgunPrep revision" >&2
  exit 1
fi
observed_commit="$(git -C "$METASHOTGUNPREP_ROOT" rev-parse HEAD)"
if [[ "$observed_commit" != "$METASHOTGUNPREP_COMMIT" ]]; then
  echo "[ERROR] MetaShotgunPrep revision mismatch" >&2
  echo "        expected: $METASHOTGUNPREP_COMMIT" >&2
  echo "        observed: $observed_commit" >&2
  exit 1
fi

# MetaShotgunPrep expects <batch>/01.RawData/<sample>/<mates> and derives names
# from this layout. Symlinks avoid duplicating the already verified downloads.
layout_root="$SAMPLE_WORK/metashotgunprep_input"
input_dir="$layout_root/01.RawData/$SAMPLE_ID"
output_root="$SAMPLE_WORK/metashotgunprep_output"
mkdir -p "$input_dir" "$output_root"
ln -sfn "$RAW_R1" "$input_dir/${SAMPLE_ID}_1.fastq.gz"
ln -sfn "$RAW_R2" "$input_dir/${SAMPLE_ID}_2.fastq.gz"

"$METASHOTGUNPREP_PYTHON" "$pipeline" \
  --sequences_dir "$input_dir" \
  --genome_index "$HOST_INDEX" \
  --output_root "$output_root" \
  --mode all \
  --threads "$PREPROCESS_THREADS"

sample_output="$output_root/metashotgunprep_input/$SAMPLE_ID"
clean_r1="$sample_output/preprocessed/${SAMPLE_ID}_1.fastq.gz"
clean_r2="$sample_output/preprocessed/${SAMPLE_ID}_2.fastq.gz"
for file in "$clean_r1" "$clean_r2"; do
  [[ -s "$file" ]] || { echo "[ERROR] Cleaned mate missing/empty: $file" >&2; exit 1; }
  gzip -t -- "$file"
done
grep -q 'Pipeline completed successfully' "$sample_output/pipeline.log" || {
  echo "[ERROR] MetaShotgunPrep completion marker absent: $sample_output/pipeline.log" >&2
  exit 1
}

# Preserve compact QC and provenance outside disposable scratch. Do not copy
# preprocessed FASTQs; downstream stages consume them in place.
qc_out="$PERSISTENT_QC_ROOT/$STUDY/$SAMPLE_ID"
mkdir -p "$qc_out"
cp -f "$sample_output/pipeline.log" "$qc_out/pipeline.log"
for directory in qc_before qc_after; do
  if [[ -d "$sample_output/$directory" ]]; then
    mkdir -p "$qc_out/$directory"
    cp -a "$sample_output/$directory/." "$qc_out/$directory/"
  fi
done
cat > "$qc_out/metashotgunprep_provenance.tsv" <<EOF
field	value
metashotgunprep_commit	$observed_commit
host_index_prefix	$HOST_INDEX
mode	all
threads	$PREPROCESS_THREADS
minimum_length	60
EOF

# Shell-safe handoff consumed by the site runner.
handoff="$SAMPLE_WORK/metashotgunprep_outputs.env"
{
  printf 'CLEAN_R1=%q\n' "$clean_r1"
  printf 'CLEAN_R2=%q\n' "$clean_r2"
  printf 'METASHOTGUNPREP_SAMPLE_OUTPUT=%q\n' "$sample_output"
  printf 'PERSISTENT_SAMPLE_QC=%q\n' "$qc_out"
  printf 'METASHOTGUNPREP_COMMIT_OBSERVED=%q\n' "$observed_commit"
} > "$handoff"

echo "[OK] MetaShotgunPrep completed: $SAMPLE_ID"
echo "[OK] Cleaned mates: $clean_r1 $clean_r2"
echo "[OK] Handoff: $handoff"
