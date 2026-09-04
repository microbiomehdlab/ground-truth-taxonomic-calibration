# Paired downstream analysis (v2)

This directory is the isolated home of the revised downstream analysis. It is
deliberately separate from the historical `original_unpaired` workflow at the
repository root and under `scripts/`. Do not modify or overwrite that workflow
or its completed run directories while developing v2.

## Why v2 exists

The experiment repeatedly measures the same biological samples after controlled
paired-read perturbations. The final analysis must therefore use the biological
sample as the unit of replication and explicitly account for repeated fractions,
targets, and profilers. The controlled truth is implanted read identity, count,
and final-library fraction—not cellular abundance.

The v2 primary quantitative endpoint is the paired, baseline-adjusted
dose-response of profiler-native species abundance. Threshold-based recovery
classes and the historical unpaired biomarker analysis are secondary historical
comparators, not the primary inference.

## Isolation contract

- New code, specifications, tests, and documentation belong under `analysis_v2/`.
- New generated results must use a new run root such as
  `RUNS_paired_dose_response_v2_<timestamp>`.
- Never point a v2 command at an existing `RUNS_publication_original_unpaired_*`
  directory for output.
- Legacy code remains reproducible from Git; its authoritative entry points are
  `run_publication_original_unpaired_q010.sh` and
  `run_original_unpaired_q010_cluster.sbatch`.
- Input staging may be shared only through immutable, checksummed profiler
  outputs and manifests. No v2 script may edit upstream evidence.

## Current status

`STATISTICAL_ANALYSIS_PLAN.md` is a draft specification. No definitive v2 model
should be fitted and no manuscript figure should be replaced until every item
marked **TO FREEZE** is resolved and the final upstream input manifests are
sealed. Development with synthetic fixtures is allowed before then.

Planned implementation order:

1. define and validate a canonical long-format input contract;
2. audit native profiler denominators and all transformations;
3. compute baseline-adjusted outcomes without model fitting;
4. implement detection dose-response models;
5. implement continuous response models and nonlinear diagnostics;
6. implement biomarker-propagation endpoints and multiplicity control;
7. add cohort-specific estimation and cross-cohort synthesis;
8. integrate Pana/Pint assembly-choice sensitivity;
9. generate tables, figures, diagnostics, and a machine-readable run manifest.

Every implementation step requires fixture-based tests before use on final
cohort results.

Step 1 is specified in `INPUT_CONTRACT.md` and enforced by
`scripts/validate_canonical_input.py`. The command-level part of step 2 is implemented in
`scripts/audit_profiler_semantics.py` and documented in
`PROFILER_SEMANTICS.md`. Run it on native cohort outputs before freezing zero
and transformation rules.

Step 3 is specified in `ENDPOINTS.md` and implemented by
`scripts/derive_paired_endpoints.py`. It distinguishes total-community dilution
from the target-specific implanted fraction and produces checksummed derived
evidence without fitting models.

Step 4 is prespecified in `DETECTION_MODEL.md` and implemented by
`scripts/fit_detection_dose_response.R`. Its synthetic execution test must run
inside the frozen analysis image, where the pinned `mgcv` dependency is present.

Step 5 is prespecified in `CONTINUOUS_MODEL.md` and implemented by
`scripts/fit_continuous_dose_response.R`. It estimates profiler-specific linear
response slopes and a categorical-dose nonlinearity diagnostic.

The Yachida Pana/Pint integration path is implemented by
`scripts/build_assembly_sensitivity_input.py`. It extracts exact frozen aliases
from native profiles, joins clean and optional original assembly arms to the
same strict baseline profiles, and validates the resulting canonical table.
`ASSEMBLY_SENSITIVITY_MODEL.md` and
`scripts/fit_assembly_sensitivity_sample_level.R` implement the primary paired
biological-sample-level response-slope analysis. The random-slope model in
`scripts/fit_assembly_sensitivity.R` is a secondary trajectory diagnostic.
The fail-closed end-to-end entry point and Lobo instructions are
`run_assembly_sensitivity.sh` and `ASSEMBLY_SENSITIVITY_RUNBOOK.md`.

Step 6's model-independent evaluation contract is documented in
`BIOMARKER_PROPAGATION.md` and implemented by
`scripts/evaluate_biomarker_propagation.py`. `PAIRED_BIOMARKER_MODEL.md`,
`scripts/build_biomarker_abundance_input.py`, and
`scripts/fit_paired_biomarker_models.R` implement the paired perturbation model.
The fail-closed entry point is `run_paired_biomarker_propagation.sh`. A separate
native disease-contrast model remains deliberately unfrozen.

`METHODS_DECISION_LOG.md` is the dated, tracked record of decisions and
development-only findings needed when drafting the manuscript. Update it when
an estimand, model, multiplicity family, or interpretation rule changes.
