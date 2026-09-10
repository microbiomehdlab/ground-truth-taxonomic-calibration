# Definitive Yachida v2 runbook

This workflow is blocked until the upstream audit reports 201/201 verified
samples and `run_yachida_production_audit.sbatch` creates a valid production
seal. Development outputs must never be passed to this runner.
The production audit independently enforces 1 baseline and 7 community
profiles for every sample, plus 60 independent profiles for each member of the
frozen nested 30-sample subset, and records these counts in the sealed dataset
completion table.

The runner deterministically builds and validates the canonical input from the
sealed native profiles. The table covers all 201 samples for community spikes
and the frozen nested 30-sample subset for independent spikes. Community target
fractions use the exact integer read allocation recorded by the frozen
community-generation algorithm; they are not assumed to be equal fractions.

Required environment variables are `PROJECT`, `YACHIDA_ENV`, `ANALYSIS_SIF`,
`ASSEMBLY_SENSITIVITY_ROOT`, and a new `RUN_ROOT`. `CANONICAL_INPUT` is
optional: when omitted, the runner creates
`$RUN_ROOT/canonical/canonical_input.tsv` from `PERSISTENT_RESULTS_ROOT`, the
frozen manifests, spike panel, and taxon aliases. An explicitly supplied
canonical table must already exist and have its adjacent validation marker.
The sensitivity seal defaults to
`$ASSEMBLY_SENSITIVITY_ROOT/experiment_seal/SUCCESS`.

Before computation, submit with `PREFLIGHT_ONLY=1`. A successful preflight
creates the canonical table and readiness report but runs no models. Then
choose a fresh `RUN_ROOT`, unset `PREFLIGHT_ONLY`, and submit the same Slurm
entry point. Alternatively, explicitly reuse the validated preflight canonical
table by setting both `CANONICAL_INPUT` and `CANONICAL_VALIDATION_SUCCESS`
while still using a fresh definitive `RUN_ROOT`. The completed package
contains native-profile semantics, paired endpoints, detection and continuous
models for independent and community populations, pooled-primary and
phenotype-stratified artificial-biomarker analyses, native disease-biomarker
analysis, integrated reports, calibration linkage, provenance hashes, and a
top-level `SUCCESS` marker.

The clean-versus-original Pana/Pint experiment is required as sealed evidence
and is rerun inside the definitive package as a separate sensitivity analysis.
It does not replace either the original independent panel or the community
experiment.

## Exact Lobo launch sequence

After all upstream sample jobs finish, submit the production audit from the
project root. A valid audit has `SUCCESS` and `production_seal.sha256` and no
`AUDIT_IN_PROGRESS`. Any failed, cancelled, or interrupted audit invalidates
the publication markers before checking the cohort and must be rerun.

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.strict-production.env"
AUDIT_JOB="$(sbatch --parsable --export=ALL run_yachida_production_audit.sbatch)"
```

After that job completes successfully, verify the seal locally:

```bash
source "$YACHIDA_ENV"
SEAL="$YACHIDA_STATE_DIR/production_seal"
test -s "$SEAL/SUCCESS"
test -s "$SEAL/production_seal.sha256"
test ! -e "$SEAL/AUDIT_IN_PROGRESS"
(cd "$SEAL" && sha256sum -c production_seal.sha256)
```

First run the downstream preflight, which builds and validates the canonical
table but fits no models:

```bash
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v1.sif
export ASSEMBLY_SENSITIVITY_ROOT="$PWD/work/yachida_assembly_sensitivity_20260901"
export RUN_ROOT="$PWD/work/analysis_v2_yachida_preflight_$(date +%Y%m%d_%H%M%S)"
export PREFLIGHT_ONLY=1
PREFLIGHT_JOB="$(sbatch --parsable --export=ALL run_yachida_definitive_analysis.sbatch)"
```

When the preflight has completed with both `$RUN_ROOT/canonical/validation/SUCCESS`
and `$RUN_ROOT/readiness/SUCCESS`, preserve its path and submit the definitive
run to a fresh directory:

```bash
PREFLIGHT_ROOT="$RUN_ROOT"
export CANONICAL_INPUT="$PREFLIGHT_ROOT/canonical/canonical_input.tsv"
export CANONICAL_VALIDATION_SUCCESS="$PREFLIGHT_ROOT/canonical/validation/SUCCESS"
export RUN_ROOT="$PWD/work/analysis_v2_yachida_definitive_$(date +%Y%m%d_%H%M%S)"
unset PREFLIGHT_ONLY
DEFINITIVE_JOB="$(sbatch --parsable --export=ALL run_yachida_definitive_analysis.sbatch)"
```

Do not submit the definitive job unless the production audit, assembly-
sensitivity experiment, and preflight all pass. The definitive top-level
`SUCCESS` and `provenance/definitive_run.sha256` are the final completion gate.
