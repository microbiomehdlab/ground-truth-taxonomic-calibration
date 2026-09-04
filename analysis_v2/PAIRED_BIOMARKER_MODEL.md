# Paired perturbation biomarker model

**Status:** frozen and fixture-tested before definitive execution.

This model asks whether adding controlled reads produces a reproducible
species-level differential-abundance call relative to the matched unmodified
library. It is run separately by cohort, phenotype background, target,
assembly arm, profiler, and dose. It is a controlled perturbation-recovery
analysis. It is not itself a CRC-versus-control contrast and must not be
described as one.

Every spiked profile is paired to the baseline profile from the same biological
sample. Profiler-native species values are converted only to fractions:
Bracken `fraction_total_reads` is unchanged and MetaPhlAn percentages are
divided by 100. For each species and sample, the response is:

`log2(spiked fraction + 1e-8) - log2(baseline fraction + 1e-8)`.

The fixed pseudocount is declared in fraction units and never learned from the
observed minimum. The feature effect is the mean paired log2 change. Its
two-sided one-sample t-test uses biological samples—not profiles—as the
replicates. BH correction is applied across all included species separately
within each cohort, study, population, phenotype background, target, assembly,
profiler, and dose.

The species universe is fixed across baseline and all doses within an analysis
family. Off-target species require nonzero prevalence in at least 10% of all
family libraries. The intended target is always retained, even when absent, so
failure to detect it cannot disappear through filtering. A target with no
variation is represented conservatively rather than silently omitted. Every
excluded feature is retained with a controlled exclusion reason.

Primary biomarker-call threshold is `q <= 0.05` with positive effect. `q <=
0.10` is a prespecified sensitivity threshold. The downstream evaluator reports
target recall, off-target burden, precision, F1, target-effect change, and call
set stability. It does not refit models or alter q-values.

The exact achieved fractions vary slightly among samples. Profiles are grouped
by their within-sample prespecified dose rank, while the output records the
median achieved fraction for that context. This avoids treating harmless
denominator differences as separate experiments.

A complementary baseline disease-contrast model is still required before
making claims about native CRC biomarkers. That model must be cohort-specific,
covariate-adjusted where justified, and frozen independently of the controlled
perturbation results.
