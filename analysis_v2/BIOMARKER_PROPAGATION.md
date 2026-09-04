# Biomarker-propagation evaluation contract

**Status:** paired perturbation model, feature policy, and evaluation metrics
implemented and fixture-tested before definitive fitting.

The study's central question is how profiler output and controlled read
perturbations propagate into biomarker conclusions. The evaluation script is
therefore separated from the program that will fit differential-abundance
models. It consumes a standardized table of already fitted feature-level
results and cannot choose covariates, filters, contrasts, or favourable calls.
Target identity is resolved through the frozen spike panel and profiler alias
table rather than by fuzzy name matching.

## Input contract

One row represents one feature in one cohort, population, target, assembly arm,
profiler, dose, and declared contrast. Required columns are:

`cohort, study, analysis_population, target_label, assembly_arm, profiler,`
`spike_fraction_target, contrast, feature, effect, p_value, q_value, include,`
`exclusion_reason`.

Rows with `include=0` require a reason and are retained in the exclusion
ledger. Included q-values must already have been adjusted within a frozen
multiplicity family by the upstream DA step. The evaluator never recomputes or
changes them.

## Metrics

For each positive dose, enriched calls have `effect > 0` and `q_value <=` the
threshold. The intended target is matched through the frozen profiler-specific
alias table. The evaluator reports:

- binary target recall;
- number of enriched calls and off-target enriched calls;
- precision and F1, with zero-call cases explicitly defined as zero;
- target effect and target q-value, including a valid non-significant target;
- change in target effect from the corresponding zero-dose result;
- Jaccard stability of the enriched biomarker set against zero dose.

The primary threshold is `q <= 0.05`; `q <= 0.10` is a labelled sensitivity
analysis retained for continuity with the historical paper. Thresholds are
evaluated together and are never selected after inspecting results. These
metrics do not themselves validate the DA model. The frozen paired model and
filtering policy are specified in `PAIRED_BIOMARKER_MODEL.md`. Native baseline
disease contrasts and cross-cohort synthesis remain separate gates.
