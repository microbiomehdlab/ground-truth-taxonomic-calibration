# Downstream v2 methods decision log

This tracked log records decisions that affect manuscript methods or
interpretation. Generated run directories preserve the corresponding inputs,
diagnostics, provenance, and checksums.

## 2026-09-04 — native profiler semantics

- Kraken2/Bracken `fraction_total_reads` and MetaPhlAn relative abundance are
  retained as profiler-native outputs; MetaPhlAn percentages are converted to
  fractions only by division by 100.
- Outputs are not claimed to measure identical cellular abundance.
- Baseline-adjusted within-profiler response is primary.
- The corrected audit passed all 720 native clean-assembly Yachida profiles.
  Definitive runs must repeat it on their exact inputs.

## 2026-09-04 — incomplete-data development run

- Ten completed strict samples were used only for integration testing while the
  prespecified 30-sample subset remained incomplete.
- The run is labelled `DEVELOPMENT_ONLY`; its estimates and p-values are not
  manuscript evidence and cannot be used for model selection.
- Its gates passed: 560 canonical rows, 500 unique native profiles, 480 paired
  positive-dose endpoints, and a converged model.

## 2026-09-04 — assembly-choice estimands

- Cleaner Pana and Pint assemblies are additive sensitivity arms; original
  results remain preserved.
- Assembly effects are estimated separately for each target and profiler.
  Equal-weight pooled effects are secondary.
- Claims concern assembly choice/quality, not a causal contamination effect,
  because strain, completeness, contiguity, and database representation differ.
- Four target-by-profiler clean-minus-original tests form the primary BH family.
  Target-specific profiler differences and pooled results are secondary.

## 2026-09-04 — repeated-measures uncertainty

- A sample random intercept alone was rejected for final assembly inference
  because it does not represent heterogeneous dose-response trajectories.
- The revised model includes sample-target-profiler random intercepts, shared
  random dose slopes, and paired clean-arm random slope deviations.
- This decision precedes the sealed 30-sample fit. Earlier random-intercept-only
  output remains diagnostic history and must not be reported.
- Literal zero p-values from numerical underflow are prohibited; extreme
  evidence is additionally represented as `-log10(p)`.

## 2026-09-04 — biological-sample-level primary assembly inference

- Development fits showed that even the random-slope GAM could give implausibly
  narrow model-based uncertainty for clean-minus-original slopes with only ten
  completed biological samples. Those development p-values are not evidence.
- The primary assembly estimator is now a two-stage paired analysis: estimate a
  six-positive-dose slope within every sample, target, profiler, and arm; then
  compare clean and original slopes within each biological sample.
- Primary uncertainty is a deterministic biological-sample bootstrap and the
  primary null test is a two-sided sample-level sign-flip test. The four
  target-by-profiler tests form one BH family.
- The GAM is retained as a secondary trajectory diagnostic. This decision was
  made before fitting the sealed 30-sample dataset and must not be changed based
  on which method produces more favorable significance.

## 2026-09-04 — paired perturbation biomarker model

- Biomarker propagation is estimated from within-sample spiked-minus-baseline
  log2 abundance changes, separately within phenotype backgrounds. Biological
  samples are the replicates; duplicated profiles are never independent.
- A fixed `1e-8` fraction pseudocount replaces outcome-dependent minimum-value
  pseudocounts. The off-target feature universe is fixed across all doses and
  requires 10% nonzero prevalence; the intended target is always retained.
- BH correction is performed across species within each target, arm, profiler,
  background, and dose. Positive calls at q <= 0.05 are primary; q <= 0.10 is a
  prespecified sensitivity analysis.
- This controlled perturbation contrast is not called CRC-versus-control.
  Baseline disease contrasts require a separate, still-unfrozen model.

## 2026-09-04 — native disease-biomarker model

- Actual disease contrasts are fitted separately by cohort, target, assembly
  arm, profiler, and dose; cohorts are not pooled for primary inference.
- CRC versus Control is primary and Adenoma versus Control is secondary. The
  primary adjustment set is age and sex; adding BMI on complete cases is a
  prespecified sensitivity analysis.
- Native abundance fractions are transformed as `log2(x + 1e-8)`. HC3 robust
  standard errors are used. The cross-dose species universe requires 10%
  prevalence, with the intended target always retained.
- BH correction is across species within the exact cohort, population, target,
  assembly, profiler, dose, contrast, and model context. Baseline disease calls
  are observed calls, so baseline-to-dose Jaccard stability is meaningful.
- Development data may exercise this model but cannot support manuscript
  estimates. Definitive claims require sealed complete cohort inputs.
