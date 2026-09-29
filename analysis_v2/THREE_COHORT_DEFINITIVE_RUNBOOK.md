# Shared three-cohort downstream workflow

Use `run_cohort_definitive_analysis.sbatch` for **each** cohort. It calls the
same `analysis_v2/run_cohort_definitive_analysis.sh`, which reads the verified
`production_seal_v2` manifest pair for Yachida, Feng, or Zeller. The code and
output layout are shared; sample counts and the Yachida-only sealed
assembly-sensitivity experiment are explicit cohort settings. Analyses remain
cohort-stratified. Do not pool samples merely because the runner is common.

The old cohort-specific launchers remain for historical runs, but new
definitive runs should use this shared entry point. A complete cross-cohort
synthesis is a separate step after all three definitive packages pass.

From the project root, set the common inputs:

```bash
export PROJECT="$PWD"
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif
export ASSEMBLY_SENSITIVITY_ROOT="$PWD/work/yachida_assembly_sensitivity_20260901"
unset COHORT_ENV CRC_ENV SPIKE_PANEL CANONICAL_INPUT CANONICAL_VALIDATION_SUCCESS
```

For a preflight, set `COHORT` to `yachida`, `feng`, or `zeller`, choose a fresh
`RUN_ROOT`, and submit the **same** job script:

```bash
export COHORT=yachida
export PREFLIGHT_ONLY=1
export RUN_ROOT="$PWD/work/analysis_v2_${COHORT}_preflight_$(date -u +%Y%m%dT%H%M%SZ)"
JOB="$(sbatch --parsable --export=ALL run_cohort_definitive_analysis.sbatch)"
echo "$JOB $RUN_ROOT"
```

Repeat for Feng and Zeller by changing only `COHORT` and selecting a new
`RUN_ROOT`. Each preflight verifies the unified seal and identity policy,
builds the canonical input, checks native production receipts and profiles,
and checks readiness. It stops before fitting models. A passing job has both
`$RUN_ROOT/canonical/validation/SUCCESS` and `$RUN_ROOT/readiness/SUCCESS`.

For a definitive run, use the same command with `PREFLIGHT_ONLY=0` and a new
`RUN_ROOT`. Do this only after the corresponding preflight passes. The runner
does not repair upstream data. Yachida additionally reruns the sealed
assembly-sensitivity analysis. A complete package has top-level `SUCCESS`
and `provenance/definitive_run.sha256`.

After all three definitive packages exist, set `YACHIDA_RUN`, `FENG_RUN`,
`ZELLER_RUN`, `ANALYSIS_SIF`, and a fresh `OUTDIR`, then run
`analysis_v2/run_three_cohort_definitive_synthesis.sh`. It checks all three
packages before producing the cross-cohort synthesis.
