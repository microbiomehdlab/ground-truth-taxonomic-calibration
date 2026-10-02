# Robustness review — 2026-10-02

Batch: `work/da_robustness_20261002T163305Z`, array 3098046,
collector 3098047. All 388 downloaded task bundles and REPORT checksums
verified locally. No computational failures. This is a methods review,
not authorization to launch definitive production.

## DA3

Across the 36 larger-sample strata (three cohorts, two profilers, adenoma
and CRC backgrounds, 10/15/20 people per group), Monte Carlo maximum-statistic
family-wide testing produced global-null any-discovery rates 0.033–0.059.
These estimate family-wise error under the tested exchangeable global null,
NOT general false-discovery rates in partially non-null clinical data.
The reassuring global-null behavior does not by itself establish strong
family-wise control under arbitrary mixtures of null and non-null features.

The parametric full-family BH comparator produced ranges 0.043–0.125 at
10/group, 0.029–0.066 at 15/group, and 0.025–0.047 at 20/group.
Under the complete null any discovery is false, so this diagnostic also
equals FDR in that specific setting. It is not a blanket conclusion about
MaAsLin2 or covariate-adjusted DA1.

Monte Carlo BH rates 0–0.022 are resolution-limited with 9,999 permutations:
the minimum attainable p is 0.0001. Do not interpret low rejection rates
as adequate power or superior calibration. Exact enumeration at 5/group
has its own discrete-resolution constraint.

Decision: permutation inference is the leading DA3 candidate; retain
parametric results as a method comparison. Freeze the multiplicity scope
(target panel versus all eligible species), null assumptions, permutation
budget and reporting before production. Maximum-statistic FWER and BH FDR
are different criteria and must not be given interchangeable labels.

## DA2

The direction candidate tests equal positive/negative sign probability among
nonzero within-person changes. It does not test zero mean change, and it
discards magnitude for its inferential statistic. Appropriate synthetic-null
any-discovery rates were 0–0.032. The skewed zero-mean scenario is explicitly
NOT a sign-null test; its large-sample rejection rate of 1 is not an error
rate for the sign null.

All six ten-pair individual contexts had zero full-family direction discoveries.
Discreteness plus multiplicity limits power; this is not evidence of failed
measurement. Community direction discovery totals include non-target species
and cannot be treated as target recovery counts or automatic false positives.

Decision: retain paired measurement/change magnitudes and direction consistency
as distinct outputs. Mean-change tests remain comparisons with their documented
skewness limitations. Do not silently substitute the sign estimand for mean
recovery or biological association.

## Next gate

The clinical target exporter has been delivered (commit 91c130c) and submitted
as 3098435, output `work/clinical_target_reference_20261002T172936Z`.
Its completed output has not yet been downloaded/reviewed locally.
Review its checksums and complete control/adenoma/CRC distributions, then
freeze biologically motivated arms without selecting doses to force failure
or success. No new upstream experiment or final production run is authorized
by this review. DA1 clinical covariate models require their own validation.
