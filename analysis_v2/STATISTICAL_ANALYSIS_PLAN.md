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

For biological sample `i`, target `t`, profiler `p`, total implanted fraction
`F`, and target-specific implanted fraction `f`:

- `o_itp`: native profiler abundance in the unmodified baseline;
- `a_itpf`: native profiler abundance after implantation;
- `e_RP = (1-F)o + f`: read-proportional reference abundance;
- `r = a - (1-F)o`: recovered spike signal;
- `R = r/f`: read-perturbation response ratio, defined only for `f > 0`.

For individual spikes, `F = f`; for community spikes, `F` is the complete
community fraction and `f` is the member-specific fraction.

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

**IMPLEMENTED.** `INPUT_CONTRACT.md` freezes the canonical row structure and
`scripts/validate_canonical_input.py` enforces pairing, uniqueness, provenance,
unit conversion, and explicit exclusions. Building the final cohort tables
remains pending upstream seals.

## 3. Endpoint hierarchy

### Primary quantitative endpoint

**DECIDED.** Estimate the dose-response of `r` against implanted fraction `f`,
with biological sample treated as a repeated-measures cluster. Report the slope,
95% confidence interval, and deviation from the read-proportional slope of one.

**IMPLEMENTED.** `scripts/derive_paired_endpoints.py` deterministically derives
the reference, `r`, `R`, and reference errors from validated canonical input.
It applies no fitted model, pseudocount, truncation, or detected-only filtering.

**PRESPECIFIED.** `CONTINUOUS_MODEL.md` defines the cohort-specific Gaussian
GAM for `r`, using profiler-by-linear-dose, condition, and sample/target random
intercepts. Primary estimands are profiler-specific slopes relative to one and
their difference. Categorical dose is a frozen nonlinearity diagnostic.

### Co-primary detection endpoint

**PRESPECIFIED, PENDING EMPIRICAL SEMANTICS GATE.** Detection is provisionally
native abundance greater than zero. Cohort-specific binomial GAMs use
categorical dose, profiler-by-dose effects, condition, and random intercepts for
biological sample and target. The primary contrast is the whole
profiler-by-dose interaction; dose-specific profiler contrasts are secondary
and BH-adjusted across positive doses. `DETECTION_MODEL.md` records the complete
model, diagnostic requirements, and failure policy. The zero definition is
frozen only after representative final native outputs pass the semantics audit.

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

**FROZEN FOR SYNTHETIC VALIDATION.** Detection and continuous models are
specified in `DETECTION_MODEL.md` and `CONTINUOUS_MODEL.md`. Final fitting still
requires upstream, semantic, canonical-input, and container execution gates;
model choice must not change in response to favourable comparative results.

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

**PARTLY VERIFIED.** The frozen-command audit is recorded in
`PROFILER_SEMANTICS.md`. A fixture-tested audit program preserves native fields
and measures their totals. It must still be run on representative final outputs
from every cohort to document:

- the denominator of Bracken species fractions and treatment of unclassified or
  undistributed reads;
- MetaPhlAn's unclassified output and whether reported rows sum to 100%;
- every downstream renormalization and filtering operation.

On 2026-09-04, the audit passed all 720 native profiles in the sealed Yachida
clean-assembly sensitivity experiment. This resolves the previously observed
MetaPhlAn hierarchy/rounding parsing issues in the audit itself. It does not
replace the required exact-input audit in each definitive cohort run.

**DECIDED.** MetaPhlAn percentages will not be multiplied by FASTQ totals and
interpreted as assigned reads. Missing Bracken species mass and the complement
of MetaPhlAn species rows will not automatically be called unclassified.
Species-closed renormalization, if retained, is sensitivity-only. Native units
and baseline-adjusted within-profiler response are primary.

**TO FREEZE.** Define detection thresholds, pseudocount policy, transformation,
handling of baseline zero, and unconditional versus detected-only summaries.
Primary unconditional quantitative summaries retain non-detections as zero.
Ratios with zero or nearly zero references will not be used without explicit
stability rules.

## 6. Biomarker-discovery propagation

**DECIDED.** The purpose is to quantify how profiler choice and controlled
sequencing evidence propagate into biomarker conclusions—not merely to rank
profilers by abundance recovery.

**FROZEN FOR CONTROLLED PERTURBATION.** `PAIRED_BIOMARKER_MODEL.md` defines the
within-sample spiked-versus-baseline contrast, fixed log2 pseudocount,
cross-dose feature universe, prevalence filter, target exception,
biological-sample test, and per-context BH family. This is not a disease-group
contrast. Report target recall, off-target burden, precision, F1, effect sizes,
confidence intervals, and biomarker-set stability. Include `q <= 0.05` and
`q <= 0.10` only in a prespecified hierarchy; do not select the threshold after
comparing results.

Feature filters, artefact lists, and calibration rules learned from outcomes must
be restricted to discovery data and frozen before validation. Grouped resampling
must keep all profiles from one biological sample in the same fold.

**IMPLEMENTED FOR CONTROLLED PERTURBATION.**
`scripts/build_biomarker_abundance_input.py` expands exact native profiles,
`scripts/fit_paired_biomarker_models.R` fits the paired feature models, and
`BIOMARKER_PROPAGATION.md` freezes the standardized result interface. The
fail-closed entry point is `run_paired_biomarker_propagation.sh`.
`scripts/evaluate_biomarker_propagation.py` deterministically computes target
recall, off-target burden, precision, F1, target-effect change, and biomarker-set
Jaccard stability at both q-value thresholds. It does not fit DA models or
recompute q-values.

**PRESPECIFIED AND IMPLEMENTED FOR DISEASE CONTRASTS.**
`DISEASE_BIOMARKER_MODEL.md` defines cohort-specific native-abundance models
with CRC versus Control primary and Adenoma versus Control secondary. The
primary adjustment set is age and sex; BMI is a complete-case sensitivity.
Models use a fixed `1e-8` fraction pseudocount, a cross-dose 10% prevalence
filter with target exception, HC3 robust uncertainty, and BH correction across
species within each exact analysis context. The fail-closed entry point is
`run_disease_biomarker_propagation.sh`. Definitive fitting still requires the
sealed final cohort inputs. Disease propagation uses the dedicated
`scripts/evaluate_disease_biomarker_propagation.py`: retained, lost, and gained
baseline disease biomarkers, Jaccard stability, direction flips, and effect
changes are the endpoints. Target significance is a spurious-association
diagnostic, not recall, because implantation is performed across phenotypes.

## 7. Cross-cohort inference

**PRESPECIFIED AND IMPLEMENTED.** Treat Yachida, Feng, and Zeller as separate
cohorts. Cohort-specific estimates are primary. Shared estimands are summarized
with REML random effects and modified Hartung–Knapp uncertainty; report tau²,
Q, I², prediction intervals, direction agreement, and leave-one-cohort-out
results. Only features estimable in every prespecified cohort are synthesized,
with incomplete coverage retained in a ledger. BH correction is across features
within each exact synthesis context. `CROSS_COHORT_SYNTHESIS.md` defines the
contract. Definitive execution remains pending all three cohort seals.

## 8. Assembly-choice sensitivity

**DECIDED.** Compare the original and cleaner Pana/Pint arms within the same 30
Yachida samples and fractions using matched seed identities. Analyse detection,
dose-response, off-target assignments, and biomarker propagation. Describe this
as assembly-choice/quality sensitivity, not a causal contamination experiment.

**IMPLEMENTED FOR QUANTITATIVE RESPONSE.**
`scripts/build_assembly_sensitivity_input.py` builds validated paired original
and clean canonical rows. The primary
`scripts/fit_assembly_sensitivity_sample_level.R` first estimates each
sample-target-profiler-arm six-dose slope and performs inference on paired
clean-minus-original sample differences. Bootstrap resampling and sign flips
therefore operate at biological-sample level. The random-slope GAM in
`scripts/fit_assembly_sensitivity.R` is a secondary trajectory diagnostic. The
fail-closed runner executes both. Detection, off-target, and
biomarker-propagation assembly comparisons remain pending.

## 9. Multiplicity and uncertainty

**PARTLY FROZEN.** For assembly sensitivity, the four target-by-profiler paired
slope differences are one primary BH family; target-pooled and profiler
difference-in-differences are secondary. For the main detection, continuous,
and biomarker analyses, define families across targets, profilers, fractions, conditions,
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

**IMPLEMENTED IN PART.** `REPORTING.md` and
`run_disease_biomarker_report.sh` implement the first modular publication
report. Every graphic has a preserved TSV source, status inheritance is
fail-closed, and the package includes diagnostics, captions, provenance, and a
checksum seal. Detection, continuous-response, assembly-sensitivity, and final
cross-cohort reporting modules remain to be integrated after their definitive
inputs are available.
