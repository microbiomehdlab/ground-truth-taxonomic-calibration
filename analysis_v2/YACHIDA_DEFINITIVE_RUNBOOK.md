# Definitive Yachida v2 runbook

This workflow is blocked until the upstream audit reports 201/201 verified
samples and `run_yachida_production_audit.sbatch` creates a valid production
seal. Development outputs must never be passed to this runner.

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
