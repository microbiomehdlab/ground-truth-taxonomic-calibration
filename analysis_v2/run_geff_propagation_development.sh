#!/usr/bin/env bash
# Controlled DEVELOPMENT_ONLY profiler-scale propagation.
# Agents never execute this on the cluster; André runs it personally.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
set -euo pipefail

: "${TARGET_GENOME_SIZES:?Set TARGET_GENOME_SIZES to the newly generated and validated target-genome-size table}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the verified analysis-v2.1 image}"
test -s "$ANALYSIS_SIF" || { echo "FAIL missing image: $ANALYSIS_SIF"; exit 1; }

analysis_python() {
  apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" \
    "$ANALYSIS_SIF" python3 "$@"
}

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
CANONICAL="work/analysis_v2_three_cohort_input_dev_20260913_180614/yachida/canonical_input.tsv"
ALIASES="examples/spike_taxon_aliases.csv"
AUDIT_ROOT="work/yachida_geff_audit_20260917T144215Z"
GEFF_PRIMARY="$AUDIT_ROOT/geff_primary_cov95/effective_genome_size.tsv"
GEFF_PRIMARY_SUCCESS="$AUDIT_ROOT/geff_primary_cov95/SUCCESS"
GEFF_SENSITIVITY="$AUDIT_ROOT/geff_sensitivity_cov90/effective_genome_size.tsv"
GEFF_SENSITIVITY_SUCCESS="$AUDIT_ROOT/geff_sensitivity_cov90/SUCCESS"
CANONICAL_SHA256="251d0ed2df12d49cad12907d66e1273e29c9416de5a2223ca613756b0c4434cd"
DB_SHA256="e7d23a73a7959b4f41af0bbe403f4b5bbb7c1879528d376d146ee5294515df9a"
RUN_ROOT="work/geff_propagation_dev_${STAMP}"

# ---------- preflight: every check is a hard gate -------------------------
for f in "$CANONICAL" "$GEFF_PRIMARY" "$GEFF_PRIMARY_SUCCESS" \
         "$GEFF_SENSITIVITY" "$GEFF_SENSITIVITY_SUCCESS" \
         "$TARGET_GENOME_SIZES" "$ALIASES"; do
  test -s "$f" || { echo "FAIL missing or empty: $f"; exit 1; }
done
test ! -e "$RUN_ROOT" || { echo "FAIL run root exists"; exit 1; }

actual="$(sha256sum "$CANONICAL" | cut -d' ' -f1)"
[ "$actual" = "$CANONICAL_SHA256" ] || { echo "FAIL canonical checksum: $actual"; exit 1; }
grep -Rqs "$DB_SHA256" "$AUDIT_ROOT" || { echo "FAIL database provenance"; exit 1; }
awk -F '\t' 'NR>1{n++} END{exit !(n==201)}' "$GEFF_PRIMARY" \
  || { echo "FAIL 95% G_eff is not 201 rows"; exit 1; }

# Target table: exactly ten expected labels, non-empty hashes, positive lengths.
awk -F '\t' '
  NR==1 { for (i=1;i<=NF;i++) c[$i]=i
          if (!(("target_label" in c) && ("fasta_sha256" in c) && ("genome_size_bp" in c)))
            { print "FAIL target table schema"; exit 1 }
          next }
  { n++; seen[$c["target_label"]]=1
    if (length($c["fasta_sha256"]) != 64) { print "FAIL empty sha256"; exit 1 }
    if (($c["genome_size_bp"]+0) <= 0)    { print "FAIL non-positive length"; exit 1 } }
  END { split("Bfrag Csym Dpne Fnuc Hhat Pmic Pana Psto Porp Pint", want, " ")
        for (i in want) if (!(want[i] in seen)) { print "FAIL missing " want[i]; exit 1 }
        if (n != 10) { print "FAIL target count " n; exit 1 } }' "$TARGET_GENOME_SIZES" \
  || exit 1

mkdir -p "$RUN_ROOT"; echo DEVELOPMENT_ONLY > "$RUN_ROOT/DEVELOPMENT_ONLY.txt"
export EFFECTIVE_GENOME_SIZES="$GEFF_PRIMARY"

# ---------- stage 0: deterministic complete native-abundance input --------
analysis_python analysis_v2/tests/test_biomarker_abundance_input.py
analysis_python analysis_v2/scripts/build_biomarker_abundance_input.py \
  --canonical "$CANONICAL" --aliases "$ALIASES" \
  --outdir "$RUN_ROOT/native_abundance"
test -s "$RUN_ROOT/native_abundance/SUCCESS" \
  || { echo "FAIL native abundance build"; exit 1; }
ABUNDANCE_LONG="$RUN_ROOT/native_abundance/biomarker_abundance_long.tsv"
test -s "$ABUNDANCE_LONG" || { echo "FAIL empty abundance table"; exit 1; }
head -1 "$ABUNDANCE_LONG" | grep -q $'profiler\tsource_profile\tfeature\tabundance_fraction' \
  || { echo "FAIL ABUNDANCE_LONG schema"; exit 1; }

# ---------- stage 1: paired endpoints, primary + sensitivity --------------
bash -c '
  set -euo pipefail
  ROOT="$PWD"; source analysis_v2/lib/metaphlan_reference.sh
  derive_endpoints_with_references "'"$CANONICAL"'" "'"$RUN_ROOT"'/endpoints"
'
PRIMARY="$RUN_ROOT/endpoints"
SENS="$RUN_ROOT/endpoints_read_reference_sensitivity"
for d in "$PRIMARY" "$SENS"; do
  test -s "$d/SUCCESS" || { echo "FAIL no SUCCESS in $d"; exit 1; }
done

# ---------- stage 2: MANDATORY Bracken regression -------------------------
bracken_cols() {
  awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
               $c["profiler"]=="kraken2_bracken" {
                 print $c["cohort"], $c["analysis_population"],
                       $c["sample_id"], $c["target_label"], $c["spike_fraction_target"],
                       $c["read_proportional_reference"], $c["recovered_spike_signal"],
                       $c["response_ratio"] }' "$1" | sort
}
if ! diff <(bracken_cols "$PRIMARY/paired_endpoints.tsv") \
          <(bracken_cols "$SENS/paired_endpoints.tsv") \
          > "$RUN_ROOT/bracken_regression.diff"; then
  echo "FAIL Bracken differs between reference arms; see $RUN_ROOT/bracken_regression.diff"; exit 1
fi
echo "[PASS] Bracken identical between profiler-scale primary and read-reference sensitivity"

# ---------- stage 3: perturbation-response input --------------------------
analysis_python analysis_v2/scripts/build_perturbation_response_input.py \
  --profile-manifest "$CANONICAL" \
  --endpoints "$PRIMARY/paired_endpoints.tsv" \
  --abundance "$ABUNDANCE_LONG" \
  --outdir "$RUN_ROOT/response_input" \
  --metaphlan-reference genome_equivalent \
  --target-genome-sizes "$TARGET_GENOME_SIZES" \
  --effective-genome-sizes "$GEFF_PRIMARY"
test -s "$RUN_ROOT/response_input/SUCCESS" || { echo "FAIL response input"; exit 1; }
RESPONSES="$RUN_ROOT/response_input/paired_feature_responses.parquet"

# ---------- stage 4: downstream on the primary scale ----------------------
analysis_python analysis_v2/scripts/summarize_target_recovery.py \
  --responses "$RESPONSES" --outdir "$RUN_ROOT/target_recovery" \
  --reference-scale profiler_scale
test -s "$RUN_ROOT/target_recovery/SUCCESS" || { echo "FAIL target recovery"; exit 1; }

analysis_python analysis_v2/scripts/analyze_perturbation_response.py \
  --responses "$RESPONSES" --outdir "$RUN_ROOT/response_analysis" \
  --reference-scale profiler_scale
test -s "$RUN_ROOT/response_analysis/SUCCESS" || { echo "FAIL response analysis"; exit 1; }

# ---------- stage 5: MetaPhlAn read-reference sensitivity -----------------
analysis_python analysis_v2/scripts/summarize_target_recovery.py \
  --responses "$RESPONSES" \
  --outdir "$RUN_ROOT/target_recovery_read_sensitivity" \
  --reference-scale read_proportional
test -s "$RUN_ROOT/target_recovery_read_sensitivity/SUCCESS" \
  || { echo "FAIL sensitivity"; exit 1; }

# ---------- reporting ------------------------------------------------------
for d in "$PRIMARY" "$SENS" "$RUN_ROOT/target_recovery"; do
  f="$d/paired_endpoints.tsv"; [ -f "$f" ] || f="$d/target_recovery_observations.tsv"
  echo "-- $f"
  awk -F '\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="reference_type")r=i;if($i=="profiler")p=i};next}
               {n++;k[$p"/"$r]++} END{printf "  rows %d\n",n;
                 for(x in k) printf "  %-42s %d\n",x,k[x]}' "$f"
done

{
  printf 'field\tvalue\n'
  printf 'stamp\t%s\n' "$STAMP"
  printf 'canonical_input\t%s\n' "$(readlink -f "$CANONICAL")"
  printf 'canonical_sha256\t%s\n' "$actual"
  printf 'geff_audit_root\t%s\n' "$(readlink -f "$AUDIT_ROOT")"
  printf 'geff_primary\t%s\n' "$GEFF_PRIMARY"
  printf 'geff_sensitivity\t%s\n' "$GEFF_SENSITIVITY"
  printf 'target_genome_sizes\t%s\n' "$(readlink -f "$TARGET_GENOME_SIZES")"
  printf 'abundance_long\t%s\n' "$(readlink -f "$ABUNDANCE_LONG")"
  printf 'abundance_sha256\t%s\n' "$(sha256sum "$ABUNDANCE_LONG" | cut -d' ' -f1)"
  printf 'bracken_comparator\t%s\n' "$(readlink -f "$SENS/paired_endpoints.tsv")"
  printf 'metaphlan_db_sha256\t%s\n' "$DB_SHA256"
  printf 'git_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'git_dirty_files\t%s\n' "$(git status --porcelain | wc -l)"
  printf 'host\t%s\n' "$(hostname)"
  printf 'user\t%s\n' "$USER"
} > "$RUN_ROOT/run_provenance.tsv"
find "$RUN_ROOT" -type f \( -name '*.tsv' -o -name '*.parquet' \) -print0 \
  | sort -z | xargs -0 sha256sum > "$RUN_ROOT/run_checksums.sha256"
echo "[DONE] $RUN_ROOT"
