#!/usr/bin/env bash
# Profile a deterministic prefix of one finalized pure target pool.
set -euo pipefail
IFS=$'\n\t'
umask 002

: "${PROJECT:?set PROJECT}"
: "${YACHIDA_ENV:?set YACHIDA_ENV}"
: "${POOL_LABEL:?set POOL_LABEL}"
source "$YACHIDA_ENV"

: "${YACHIDA_POOLS_DIR:?set in YACHIDA_ENV}"
: "${PURE_POOL_AUDIT_ROOT:?set PURE_POOL_AUDIT_ROOT}"
PURE_POOL_AUDIT_PAIRS="${PURE_POOL_AUDIT_PAIRS:-1000000}"
[[ "$POOL_LABEL" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "[ERROR] Unsafe pool label" >&2; exit 1; }
[[ "$PURE_POOL_AUDIT_PAIRS" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] PURE_POOL_AUDIT_PAIRS must be positive" >&2; exit 1; }

index="$YACHIDA_POOLS_DIR/pool_pair_counts.tsv"
seal="$YACHIDA_POOLS_DIR/pool_pair_counts.tsv.sha256"
(cd "$YACHIDA_POOLS_DIR" && sha256sum -c "$(basename "$seal")")
available="$(awk -F'\t' -v label="$POOL_LABEL" '$1 == label {print $4}' "$index")"
[[ "$available" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Unknown pool: $POOL_LABEL" >&2; exit 1; }
(( PURE_POOL_AUDIT_PAIRS <= available )) || { echo "[ERROR] Requested more pairs than pool contains" >&2; exit 1; }

root="$PURE_POOL_AUDIT_ROOT/$POOL_LABEL"
sample_work="$root/work"
profiles="$root/profiles"
mkdir -p "$sample_work" "$profiles"
r1="$sample_work/${POOL_LABEL}.pure_pool_1.fastq.gz"
r2="$sample_work/${POOL_LABEL}.pure_pool_2.fastq.gz"
python3 "$PROJECT/scripts/systematic_sample_paired_fastq.py" \
  --r1 "$YACHIDA_POOLS_DIR/${POOL_LABEL}.pool_1.fq" \
  --r2 "$YACHIDA_POOLS_DIR/${POOL_LABEL}.pool_2.fq" \
  --output-r1 "$r1.tmp" --output-r2 "$r2.tmp" \
  --total-pairs "$available" --sample-pairs "$PURE_POOL_AUDIT_PAIRS"
mv -f "$r1.tmp" "$r1"
mv -f "$r2.tmp" "$r2"
python3 "$PROJECT/scripts/validate_paired_fastq.py" \
  --r1 "$r1" --r2 "$r2" --minimum-pairs "$PURE_POOL_AUDIT_PAIRS" \
  --output "$root/fastq_integrity.tsv"

cat > "$root/audit_design.tsv" <<EOF
field\tvalue
pool_label\t$POOL_LABEL
selection_rule\tdeterministic_systematic_midpoints_across_full_pool
selected_pairs\t$PURE_POOL_AUDIT_PAIRS
available_pairs\t$available
pool_index_sha256\t$(sha256sum "$index" | awk '{print $1}')
source_r1\t$YACHIDA_POOLS_DIR/${POOL_LABEL}.pool_1.fq
source_r2\t$YACHIDA_POOLS_DIR/${POOL_LABEL}.pool_2.fq
EOF

export SAMPLE_ID="pure_pool_${POOL_LABEL}"
export SAMPLE_WORK="$sample_work"
export PROFILE_R1="$r1"
export PROFILE_R2="$r2"
export YACHIDA_BASELINE_SMOKE_ROOT="$profiles"
export ALLOW_NO_METAPHLAN_SPECIES=true
bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"

metaphlan_profile="$profiles/$SAMPLE_ID/${SAMPLE_ID}.metaphlan.tsv"
if grep -qE '(^|[|])s__[^|[:space:]]+' "$metaphlan_profile"; then
  metaphlan_species_status="species_rows_present"
else
  metaphlan_species_status="no_species_rows"
fi
printf 'pool_label\tmetaphlan_species_status\n%s\t%s\n' \
  "$POOL_LABEL" "$metaphlan_species_status" > "$root/profile_status.tsv"

sha256sum "$root/audit_design.tsv" "$root/fastq_integrity.tsv" "$root/profile_status.tsv" > "$root/audit_inputs.sha256"
touch "$root/SUCCESS"
echo "[PASS] Pure-pool audit completed: $POOL_LABEL"
