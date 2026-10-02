# DA2 paired inference amendment — 2026-10-02

Agreed after inspecting the engineering pilot; not preregistration and not a
claim that definitive results or null validation are complete.

## Scope and rationale

DA2 primary inference is a two-sided one-sample t test of biological-person
differences (spiked minus original). DA1 and DA3 retain MaAsLin2 primary.
The paired MaAsLin2 random-intercept analysis remains a separately reported
sensitivity, not a fallback selected by significance or optimizer outcome.

Evidence: all 12 DA2 engineering contexts failed the strict mixed-model gate;
the diagnostic refit job 3097924 retained 45 features and 225 model records.
Optimizers did not uniformly resolve warnings and degrees-of-freedom instability.
Higher mixed-model degrees of freedom are not inherently proof of incorrect
pairing; this amendment chooses an explicit person-level difference estimand.

## Frozen statistical contract

- Use unmodified native profiler abundance fractions, not recovery corrections.
- Transform each observation as log2(1+a/1e-8), identical to the previous
  shifted log2(a+1e-8) contract. Pair by biological ID, never by row position.
- Require exactly one original (group 0) and one spiked (group 1) per person;
  forbid missing observations and silent pair removal. Require at least two pairs.
- Effect is mean transformed difference; SE is sample SD divided by sqrt(n);
  degrees of freedom n-1. Report two-sided p and 95% t confidence interval.
- Preserve the pre-existing frozen family and apply external BH to its entirety.
  Positive discovery requires an estimable positive effect and q<=0.05.
- Degenerate differences: raw p and CI are NA, retain the effect/SE and status,
  and use 1 only for BH bookkeeping. Never drop a species or invent a tiny p.
  Numerical SE tolerance is 100 times machine epsilon times max(1,absolute
  transformed response), also respecting the t-test constant-input guard.
  Separate no-change from nonzero constant differences. This conservative
  treatment is a limitation, not evidence that a constant shift was absent.
- No new independent-person analysis substitutes for paired DA2. Genuine
  distinct-person discovery remains DA3.

## Validation and output

The separate parallel pilot reuses the sealed 84-context plan but runs only
the 12 DA2 contexts. Each context retains native inputs, person differences,
results, method, session information, code/input hashes and engineering status.
Prior model outputs remain unchanged. Runtime source/image hash checks and
resumption gates remain enforced. This is not a production launch.

Before production: review all 12 outputs including non-estimable counts,
validate paired-null behavior with sufficient replicates and explicitly justified
real-data null construction, evaluate pseudocount sensitivity, and freeze the
final method/aggregation contract. Synthetic fixture success alone cannot prove
false-positive control in these data. Normality of differences, small n=10,
zero-heavy abundances and compositional changes remain relevant limitations.

## Cluster pilot

Run `tests/test_paired_difference_context.R` inside the pinned analysis image.
Use `run_paired_difference_pilot.sbatch` with the existing `DA_PILOT_PLAN`,
fresh `DA2_PAIRED_RESULTS`, and DA2 indices derived from `contexts.tsv`.
The worker refuses non-DA2 or unpaired contexts. Review results before proceeding.
