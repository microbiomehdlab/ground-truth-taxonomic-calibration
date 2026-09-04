#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${SENSITIVITY_ROOT:?Set SENSITIVITY_ROOT to the sealed clean-arm run root}"
: "${ORIGINAL_ROOT:?Set ORIGINAL_ROOT to the sealed strict original-arm run root}"
: "${BASELINE_ROOT:=$ORIGINAL_ROOT}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty analysis-v2 run directory}"
MANIFEST="${MANIFEST:-work/yachida_67x3/metadata/independent_10_per_condition.tsv}"
ARMS="${ARMS:-datasets/yachida/assembly_sensitivity_arms.tsv}"
SPIKE_PANEL="${SPIKE_PANEL:-spikes/spike_panel.tsv}"
ALIASES="${ALIASES:-examples/spike_taxon_aliases.csv}"

if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2
  exit 1
fi
for path in "$MANIFEST" "$ARMS" "$SPIKE_PANEL" "$ALIASES" "$ANALYSIS_SIF"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
mkdir -p "$OUTDIR"/{canonical,profiler_semantics,endpoints,models,provenance}

{
  printf 'field\tvalue\n'
  printf 'status\tIN_PROGRESS\n'
  printf 'created_at\t%s\n' "$(date -Iseconds)"
  printf 'repository_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'manifest\t%s\n' "$(realpath "$MANIFEST")"
  printf 'arms\t%s\n' "$(realpath "$ARMS")"
  printf 'spike_panel\t%s\n' "$(realpath "$SPIKE_PANEL")"
  printf 'aliases\t%s\n' "$(realpath "$ALIASES")"
  printf 'sensitivity_root\t%s\n' "$(realpath "$SENSITIVITY_ROOT")"
  printf 'original_root\t%s\n' "$(realpath "$ORIGINAL_ROOT")"
  printf 'baseline_root\t%s\n' "$(realpath "$BASELINE_ROOT")"
  printf 'analysis_sif\t%s\n' "$(realpath "$ANALYSIS_SIF")"
} > "$OUTDIR/provenance/run_manifest.tsv"

python3 analysis_v2/scripts/build_assembly_sensitivity_input.py \
  --manifest "$MANIFEST" --arms "$ARMS" --spike-panel "$SPIKE_PANEL" \
  --aliases "$ALIASES" --sensitivity-root "$SENSITIVITY_ROOT" \
  --original-root "$ORIGINAL_ROOT" --baseline-root "$BASELINE_ROOT" \
  --outdir "$OUTDIR/canonical"

audit_args=(--outdir "$OUTDIR/profiler_semantics")
while IFS= read -r profile; do
  case "$profile" in
    *.bracken.S.tsv) audit_args+=(--bracken "$profile") ;;
    *.metaphlan.tsv) audit_args+=(--metaphlan "$profile") ;;
    *) echo "[ERROR] Unexpected native profile: $profile" >&2; exit 1 ;;
  esac
done < <(awk -F '\t' 'NR > 1 {print $20}' "$OUTDIR/canonical/canonical_input.tsv" | sort -u)
python3 analysis_v2/scripts/audit_profiler_semantics.py "${audit_args[@]}"

python3 analysis_v2/scripts/derive_paired_endpoints.py \
  --input "$OUTDIR/canonical/canonical_input.tsv" --outdir "$OUTDIR/endpoints"

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_assembly_sensitivity_sample_level.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/fit_assembly_sensitivity_sample_level.R \
  --input "$OUTDIR/endpoints/paired_endpoints.tsv" \
  --outdir "$OUTDIR/models/assembly_sensitivity_primary" \
  --cohort yachida --population independent

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_assembly_sensitivity_model.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/fit_assembly_sensitivity.R \
  --input "$OUTDIR/endpoints/paired_endpoints.tsv" \
  --outdir "$OUTDIR/models/assembly_sensitivity_gam_secondary" \
  --cohort yachida --population independent

sha256sum "$MANIFEST" "$ARMS" "$SPIKE_PANEL" "$ALIASES" "$ANALYSIS_SIF" \
  "$OUTDIR/canonical/canonical_input.tsv" \
  "$OUTDIR/profiler_semantics/profile_semantics_audit.tsv" \
  "$OUTDIR/endpoints/paired_endpoints.tsv" \
  "$OUTDIR/models/assembly_sensitivity_primary/primary_assembly_effects.tsv" \
  "$OUTDIR/models/assembly_sensitivity_primary/sample_paired_slope_differences.tsv" \
  "$OUTDIR/models/assembly_sensitivity_gam_secondary/assembly_profiler_response_slopes.tsv" \
  "$OUTDIR/models/assembly_sensitivity_gam_secondary/assembly_slope_contrasts.tsv" \
  > "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
sed -i 's/^status\tIN_PROGRESS$/status\tPASS/' "$OUTDIR/provenance/run_manifest.tsv"
printf 'analysis\tyachida_assembly_choice_sensitivity\nstatus\tPASS\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed assembly-choice sensitivity analysis: $OUTDIR"
