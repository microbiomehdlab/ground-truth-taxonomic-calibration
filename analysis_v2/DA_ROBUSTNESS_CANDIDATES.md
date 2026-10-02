# Concurrent DA robustness candidates — 2026-10-02

This development can proceed while the exact DA3 5+5 batch runs. It does not
adopt either candidate or authorize production. Uniform-spike pilot detection
is a technical benchmark, not a desired outcome or clinical validation.

## Larger DA3 randomized-label null

360 shards: 12 cohort/profiler/background strata times n=10,15,20 times 10
blocks of 100 allocations. Reuse all 1,000 frozen full-pool allocations per
stratum/sample size. Compare existing parametric full-family BH with 9,999
fixed-budget Monte Carlo permutations, full-family BH and separately labelled
family-maximum adjusted p. Statistic/standardization matches the exact n=5
candidate. Use one label assignment across all species to preserve dependence.
Seed namespace excludes profiler so matched biological allocations use the same
random assignments. Finite samples of permutations use (exceedances+1)/(B+1),
not zero p values, and conservative floating-point ties. No optional stopping.

BH p resolution is 1e-4 at this budget: this cannot resolve an isolated feature
at q=.05 over 3,471 species. Conservative BH discovery is NOT evidence of
insufficient biological signal. The family-maximum comparison avoids that
specific BH resolution requirement, but addresses a different error criterion.
Its validity here is under a global exchangeable-label null, not arbitrary
heteroscedastic mean-null populations or observational clinical confounding.
Individual adjusted p-values near .05 have Monte Carlo uncertainty; no definitive
feature decisions or production conclusions from this candidate budget.

This batch is a null-calibration comparison, not a complete larger-n power
study. Larger-n spiked scenarios must follow the approved biological matrix.
Only the first allocation's feature results are stored per shard; all 100 draw
discovery counts and their deterministic seeds are retained. Report summary
rates separately by cohort, profiler, background and n, with complete denominators.

## DA2 direction-based sensitivity, NOT a replacement for mean inference

12 saved paired pilot contexts: exact binomial sign test for
Pr(D>0 | D!=0)=.5. Save positive/negative/tied counts, informative n and every
feature p/q. Numerically tied differences are explicitly counted, never silently
treated as independent extra samples. Globally unchanged features remain in the
family with raw p NA and bookkeeping p=1. Independence is biological-person
independence. No symmetry assumption is required for this sign-null, but it is
NOT the null of zero mean difference and not a universal skewness repair.

16 synthetic comparisons: n=10,42,47,67; 1,000 replicates; Gaussian, sparse
symmetric, exponential centered at its median (valid sign null), and exponential
centered at its mean (different null). The last scenario is labelled NOT_A_SIGN_NULL
and rejection must not be called sign-test error. Exact sign discreteness also
limits full-family discovery with only 10 pairs. Never narrow the hypothesis
family after observing results to rescue significance.

## Execution and next decision

One submission: 388 shards, bounded concurrency and after-any report. Source is
snapshotted so active jobs survive subsequent repository development/pulls.
All outputs are separate and hash verified. These tests compare inference
assumptions and resolution limits, not optimize spike detection. User must approve
estimands/error criteria before a new primary protocol is frozen. DA1 remains a
separate clinical/covariate validation task; neither candidate certifies it.

Implementation reference for Monte Carlo plus-one and tie conventions:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html
