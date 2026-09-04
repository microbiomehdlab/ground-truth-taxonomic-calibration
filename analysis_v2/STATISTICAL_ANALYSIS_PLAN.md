# Statistical analysis plan — paired dose-response v2

**Status:** DRAFT; created before the final three-cohort comparative analysis.

This document separates decisions already made from choices that must be frozen
before examining definitive v2 comparisons. Updating a **TO FREEZE** item after
viewing its corresponding final result requires a dated deviation record.

## 1. Scientific estimand

**DECIDED.** Synthetic paired reads create controlled species-specific
perturbations to the final sequencing library. They do not establish cellular
abundance, biomass, extraction efficiency, or an identical biological abundance
estimand across profilers.

For biological sample `i`, target `t`, profiler `p`, and implanted fraction `f`:

- `o_itp`: native profiler abundance in the unmodified baseline;
- `a_itpf`: native profiler abundance after implantation;
- `e_RP = (1-f)o + f`: read-proportional reference abundance;
- `r = a - (1-f)o`: recovered spike signal;
- `R = r/f`: read-perturbation response ratio, defined only for `f > 0`.

**DECIDED.** Profiler-native output is the primary practical representation.
Denominator-harmonized quantities are sensitivity analyses and will not be
described as making Kraken2/Bracken and MetaPhlAn biologically identical.

## 2. Analysis populations

**DECIDED.** The community analysis uses every outcome-independent eligible
sample that passes the frozen technical pipeline. The independent-target analysis
uses the frozen nested 10/10/10 subset in each cohort. Transfer or processing
failures are retried and are not biological exclusions.

**TO VERIFY.** Produce a final participant/sample flow table for each cohort and
a zero-unexplained-missingness ledger before model fitting.

## 3. Endpoint hierarchy

### Primary quantitative endpoint

**DECIDED.** Estimate the dose-response of `r` against implanted fraction `f`,
with biological sample treated as a repeated-measures cluster. Report the slope,
95% confidence interval, and deviation from the read-proportional slope of one.

### Co-primary detection endpoint

**TO FREEZE.** Define detection from the unmodified native species output and
model detection probability across implanted fractions. Freeze whether the
primary contrast is profiler difference in the whole curve, a prespecified weak
dose, or the fraction giving a prespecified detection probability.

### Secondary endpoints

- response ratio `R` distributions;
- absolute percentage-point deviation from `e_RP`;
- monotonicity and evidence of nonlinear response or saturation;
- detected-only response, explicitly labelled conditional;
- descriptive agreement classes under prespecified thresholds;
- off-target abundance changes and de novo detections;
- biomarker-call propagation, effect-size changes, and call stability.

## 4. Repeated-measures model

**DECIDED.** Samples, not profile rows, are the independent biological units.
Inference must account for repeated fractions and targets within samples.

**TO FREEZE.** Select the exact model family and formula after simulation and
diagnostic work on synthetic fixtures, not by choosing whichever gives the most
favourable profiler comparison. Candidate continuous models include a linear
mixed model on a justified transformed scale and a flexible spline model with
sample-clustered uncertainty. Candidate detection models include binomial mixed
models or marginal models with sample-clustered standard errors.

Prespecify:

- fixed effects and interactions for profiler, fraction, target, condition, and
  cohort;
- random-intercept and any random-slope structure;
- nonlinear terms and the test used to retain them;
- covariates justified independently of outcomes;
- convergence criteria and fallback models;
- influence, residual, calibration, and sensitivity diagnostics.

Report cohort-specific estimates before any pooled or hierarchical summary.

## 5. Denominators, zeros, and transformations

**TO VERIFY.** Audit frozen commands and representative raw outputs to document:

- the denominator of Bracken species fractions and treatment of unclassified or
  undistributed reads;
- MetaPhlAn's unclassified output and whether reported rows sum to 100%;
- every downstream renormalization and filtering operation.

**TO FREEZE.** Define detection thresholds, pseudocount policy, transformation,
handling of baseline zero, and unconditional versus detected-only summaries.
Primary unconditional quantitative summaries retain non-detections as zero.
Ratios with zero or nearly zero references will not be used without explicit
stability rules.

## 6. Biomarker-discovery propagation

**DECIDED.** The purpose is to quantify how profiler choice and controlled
sequencing evidence propagate into biomarker conclusions—not merely to rank
profilers by abundance recovery.

**TO FREEZE.** Define primary disease contrasts, model formula, covariates,
feature filtering, effect-size scale, and multiplicity families before the final
run. Report target recall, off-target burden, precision, F1, effect sizes,
confidence intervals, and biomarker-set stability. Include `q <= 0.05` and
`q <= 0.10` only in a prespecified hierarchy; do not select the threshold after
comparing results.

Feature filters, artefact lists, and calibration rules learned from outcomes must
be restricted to discovery data and frozen before validation. Grouped resampling
must keep all profiles from one biological sample in the same fold.

## 7. Cross-cohort inference

**TO FREEZE.** Treat Yachida, Feng, and Zeller as separate cohorts. Estimate each
cohort first, quantify heterogeneity, and then use a prespecified hierarchical or
meta-analytic synthesis. A pooled analysis that ignores cohort is not primary.

## 8. Assembly-choice sensitivity

**DECIDED.** Compare the original and cleaner Pana/Pint arms within the same 30
Yachida samples and fractions using matched seed identities. Analyse detection,
dose-response, off-target assignments, and biomarker propagation. Describe this
as assembly-choice/quality sensitivity, not a causal contamination experiment.

## 9. Multiplicity and uncertainty

**TO FREEZE.** Define families across targets, profilers, fractions, conditions,
cohorts, and endpoints. Distinguish confirmatory from exploratory tests. Report
effect sizes and uncertainty even when multiplicity-adjusted significance is not
reached. Bootstrap or resampling procedures must operate at sample level.

## 10. Required outputs and reproducibility

The definitive workflow must write:

- canonical input and exclusion manifests with checksums;
- analysis specification and code commit;
- container/database/software identities;
- endpoint tables with estimates, intervals, raw and adjusted p-values;
- model diagnostics and convergence records;
- cohort-specific and synthesized results;
- denominator and zero-handling audit;
- assembly-sensitivity results;
- figure-source tables;
- a machine-readable run manifest and checksum seal;
- a dated protocol-deviation ledger, including an explicit `none` entry when
  applicable.

The v2 run must fail rather than silently drop samples, taxa, models, or files.
