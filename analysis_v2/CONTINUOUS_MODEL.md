# Continuous recovered-signal dose-response model

**Status:** prespecified and implemented against synthetic fixtures; final-data
fit pending completion of the upstream and canonical-input gates.

## Primary estimand

The response is the baseline-adjusted recovered spike signal
`r = a - (1 - F)o`, where `a` and `o` are profiler abundances converted to
fractions by unit scaling only, `F` is total implanted fraction, and `f` is the
target-specific implanted fraction. The primary estimand is the
profiler-specific linear slope of `r` against `f`; one is read-proportional.

For numerical stability, the model uses percentage points (`dose_pp = 100 f`).
Its ideal slope is 0.01 abundance-fraction units per percentage point. Reported
`response_slope` divides the fitted slope by 0.01, returning the scale where one
is ideal.

## Primary model

Models are separate by cohort and population, with the main analysis restricted
to the `original`/`not_applicable` arm:

`r ~ profiler * dose_pp + condition + s(sample_id, bs="re") + s(target_label, bs="re")`.

The Gaussian model uses `mgcv::gam(..., method="ML")`. Sample and target random
intercepts account for repeated observations and target heterogeneity. The
primary between-profiler contrast is the profiler-by-dose coefficient.
Profiler-specific slopes are tested against one, with BH adjustment across the
two tests within each fitted cohort/population/arm family.

The intercept is estimated rather than forced to zero, allowing systematic
offset to be diagnosed. It is not baseline detection: only positive-dose rows
enter this model.

## Nonlinearity and diagnostics

A prespecified categorical-dose model replaces linear dose while retaining the
same adjustment and random effects. Its likelihood-ratio comparison diagnoses
nonlinearity but cannot replace the primary linear slope after results are seen.

The implementation records convergence, Hessian status, residuals by dose and
profiler, coefficients, fitted rows, model comparison, session information,
checksums, and `SUCCESS`. Negative recovered signals remain. No pseudocount,
truncation, detected-only filtering, or response-ratio filtering is applied.
Failure of convergence or Hessian checks stops the run; there is no silent
fallback.
