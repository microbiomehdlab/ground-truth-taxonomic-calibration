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
