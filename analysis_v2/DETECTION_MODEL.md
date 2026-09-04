# Detection dose-response model

**Status:** prespecified and implemented against synthetic fixtures; final-data
fit pending completion of the upstream and canonical-input gates.

## Endpoint

The model consumes the validated canonical table, including its zero-dose rows,
and uses `detected_native_nonzero`, the mechanical indicator that the native
target abundance is greater than zero. This threshold is provisional
until the empirical native-output audit confirms that zero has the same
operational meaning throughout each frozen profiler workflow. Any change after
that audit must be documented before comparative results are inspected.

## Primary model

Models are fitted separately by cohort and analysis population, with the main
analysis restricted to the `original`/`not_applicable` assembly arm. Implanted
target fraction is a categorical dose, including zero. The binomial logit model
is:

`detected ~ profiler * dose + condition + s(sample_id, bs="re") + s(target_label, bs="re")`.

This is fitted by `mgcv::gam(..., method="ML")`, allowing comparison of the
full and reduced fixed-effect structures. The sample random intercept
accounts for repeated targets, profilers, and doses within a biological sample;
the target random intercept accounts for pooled target heterogeneity. Cohorts
are not treated as exchangeable rows in one pooled primary model.

The primary profiler comparison is the likelihood-ratio test of the complete
profiler-by-dose interaction against the otherwise identical model without that
interaction. Dose-specific profiler contrasts are secondary and receive
Benjamini-Hochberg adjustment within each fitted cohort/population/arm family.
Predicted detection curves are population-level predictions with both random
effects excluded and condition averaged with observed-condition weights.

## Diagnostics and failure policy

Rows with `include=0` are excluded explicitly; the canonical exclusion ledger
remains authoritative. The implementation requires both profilers, at least two samples and targets,
at least two doses including zero, both outcome classes, finite coefficients,
successful full convergence, and a positive-definite outer Hessian. It writes no `SUCCESS`
marker unless these checks and all expected outputs pass. There is no automatic
fixed-effect fallback. Separation, singular information, or convergence failure
is retained as a failed diagnostic and requires a prespecified remedial or
descriptive analysis—not silent model substitution.

The assembly-choice experiment is fitted separately with the same structure;
it is not mixed into the primary original-assembly model. Continuous-dose and
nonlinear descriptions remain secondary sensitivity analyses to be added only
after the categorical primary model is validated.
