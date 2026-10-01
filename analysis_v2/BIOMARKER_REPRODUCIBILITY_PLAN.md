# Controlled biomarker discovery and reproducibility

Date: 30 September 2026.
Status: scientific direction AGREED; detailed new protocol PROPOSED pending the
design review below. This is a dated revision informed by prior development
results, not an untouched prospective or preregistered study.

## Central question and scope

Under what conditions can a known microbial addition be reproducibly discovered
across taxonomic profilers and cohorts, particularly in adenoma backgrounds?

The main paper will emphasize native observations, controlled recovery,
artificial-marker replication, and controlled discovery experiments. Learned
feature removal and abundance calibration are removed from the required main
analysis. Their historical code and negative results remain available. The
profiler-specific expected-abundance correction remains necessary: choosing an
appropriate measurement reference is distinct from correcting profiler output.

The implanted truth comprises read identity, counts, and achieved added fraction.
It does not establish native biological associations, cellular abundance, or
the reasons CRC-associated taxa may differ in adenomas. The same addition can
have different effects against different biological backgrounds. Report the
conditions of recovery rather than a universal abundance threshold for current
tools. Conclusions apply to the tested workflow/database configurations.

Full downstream submissions are on hold. Existing upstream data are sufficient
to develop this design; no new upstream profiling is proposed.

## Evidence checkpoint

All 511 samples are upstream-sealed under the common v2 interface. Shared
preflights passed for Yachida (3097753, 201 samples, 36,360 canonical rows),
Feng (3097754, 154 samples, 28,840 rows), and Zeller (3097755, 156 samples,
29,160 rows). These user-reported results validate input preparation, not the
unexecuted model stages. The prior full-cohort runner is an implementation
starting point, not the complete protocol defined here.

## Analysis matrix

| ID | Scientific question and population | Endpoint and method | Uncertainty / validation | Intended display | Current implementation |
|---|---|---|---|---|---|
| A1 | How do native candidate associations differ? Full unspiked cohorts, CRC/Control and Adenoma/Control | Existing age/sex-adjusted log-abundance model, HC3 intervals and baseline-frozen feature universe. CRC is the existing primary disease contrast; adenoma is central contextual evidence in this revised narrative, without silently changing its testing hierarchy. | Effect intervals, counts, BMI sensitivity, common/evaluable taxon mapping; report both contrasts, including null findings | Candidate effect and prevalence panel by cohort/profiler | Models exist; definitive baseline rebuild/review needed |
| A2 | Can known additions be detected and quantified? Separate independent and community populations | Native nonzero detection; unconditional paired response and error against declared profiler-specific reference. Existing GAMs plus descriptive taxon/dose summaries. | Sample dependence respected; real-model diagnostics; read-reference sensitivity; assess whether pooled models support the intended heterogeneity claims | Detection and corrected recovery curves | Implemented; final reference coverage, execution and consumer migration incomplete |
| A3 | Do artificial discoveries agree across cohorts and profilers? All ten prespecified targets at matched design/background/dose rank | Existing paired log2-change tests; positive effect and primary BH q <= 0.05, q <= 0.10 sensitivity. Count target discoveries and pairwise/joint recovery patterns. | Retain all ten intended targets in eligibility ledger; compare matched sample sizes; report effects/p-values/q-values and feature-family sizes | Dose versus joint recovery; taxon agreement atlas | Paired fits exist; explicit replication evaluator is new |
| A4 | How do dose and study size affect discovery in adenoma backgrounds? Community-spiked and baseline profiles from distinct biological samples | Randomly allocated pseudo cases/controls; one state per sample per experiment; same allocations across profilers/doses. New between-sample DA module and repeated discovery/replication summaries. | Frozen seeds and sample-size grid; null-label baseline arm; full feature/BH fitting; report conditional resampling frequency separately from population uncertainty | Signal-by-sample-size reproducibility map | New design and implementation |
| A5 | Which measured factors accompany reproducibility failure, and how does this inform native observations? | Prespecified summaries of baseline, prevalence, paired response variability, effect, uncertainty and multiplicity. Native association comparison is contextual. | No causal attribution from correlations; no many-predictor model over ten taxa; avoid outcome-derived classification of the same observations | Mechanistic summary and qualified native-context comparison | Existing descriptive pieces; new integration and interpretation rules |

Non-implanted responses and community reconstruction support A2–A5 where useful.
Native disease perturbation/ideal-counterfactual distortion remains an optional
mechanistic module, mandatory only if the paper makes formal claims about
phenotype-dependent distortion under perturbation. Cross-cohort meta-analysis
does not replace explicit discovery/replication assessment.

## Definitions for the new replication layer

Discovery is a positive effect with q <= 0.05 in a prespecified eligible test.
The q <= 0.10 result is reported as sensitivity, never chosen after comparison.

For target t, cohort c, profiler p, background b and dose d, retain a binary
discovery indicator along with its effect, interval, p, q, sample count and
number of features in the BH family. Pairwise cohort results are both, first
only, second only, neither. Profiler results use the same four-way structure.
Three-cohort results preserve the eight exact discovery patterns plus the
summary count 0/1/2/3.

Report two different denominators:

- Joint recovery: targets discovered in both cohorts divided by all evaluable
  intended targets. Show numerator and denominator; normally ten targets.
- Directional replication: targets discovered in both divided by targets
  discovered in the designated discovery cohort and evaluable in validation.
  No discovery denominator means NA, not zero or perfect replication.

An observed absent target is a recovery failure if its profile is complete;
a missing/corrupt profile is an input failure, never a negative discovery.
Use exact frozen taxon identities across workflows. Do not assert that their
whole native feature universes or q-values are identical. Preserve p-values and
BH-family sizes to help distinguish measurement from multiplicity differences.

PROPOSED primary cross-cohort display: all three pairwise comparisons, with
all-three recovery as a stricter secondary summary. Do not treat the three
pairs as independent observations. The ten targets are a selected panel, not
a random sample of all microbial species. Avoid generalizing uncertainty over
this panel to all taxa.

## A4: proposed controlled discovery protocol

### Feasible sample populations

The sealed adenoma populations are Feng 47, Zeller 42, and Yachida 67 before any
model-specific exclusion. Validate these counts and full profile availability
from the canonical manifests. Each independent-spike background has only ten
samples per cohort, allowing at most five per synthetic group without reuse.

Use the community experiment as the main adenoma discovery experiment. Its
ten taxa are added together; this must not be described as a single-marker
intervention. Treat the individual-spike pseudo-group analysis as a small-sample
sensitivity, if retained. Control and CRC backgrounds are contextual checks.

PROPOSED common grid: 5, 10, 15 and 20 biological samples per group, subject to
input and covariate feasibility checks. Cohort-specific larger panels are
secondary and must not be used to claim superior cohort recoverability.

### Allocation and observations

For each cohort/background/replicate, draw 2n distinct biological samples without
replacement and randomly allocate n to pseudo cases and n to pseudo controls.
Use post-spike profiles only for pseudo cases and baseline profiles only for
pseudo controls. Never put a sample's spiked and baseline versions in opposite
groups within the same experiment. Reuse the allocation across both profilers
and all doses; nested sample-size allocations are recommended for paired
comparisons. Draw allocations independently across cohorts.

Use cohort plus sample ID as the physical sample key. Retain a manifest of
every allocation, dose, source profile, seed and model context. Reusing people
across repeated experiments is allowed but creates dependent summaries.

### Statistical model and feature family

PROPOSED primary model: feature-wise `log2(abundance_fraction + 1e-8) ~
pseudo_case`, with HC3 uncertainty, followed by BH over the complete eligible
native feature family. Random allocation makes this a controlled artificial
group contrast, not the existing paired test and not a native disease contrast.
Review small-n behavior, constant features and failed fits before freezing it.
Keep a rank/nonparametric sensitivity optional until its scientific purpose
is agreed; do not add methods merely to select favorable results.

PROPOSED feature universe: the cohort/profiler baseline-frozen species set
(10% baseline prevalence plus ten prespecified targets), shared across
replicates, doses and model arms. This describes conditional performance in
these observed cohorts. It differs from prospective training-only feature
selection and must be labelled accordingly. Alternative within-replicate
filtering would answer a different pipeline-performance question.

Quantify how the whole native feature universe affects target q-values. A
predeclared target-panel-only multiplicity analysis could be a secondary
diagnostic, never substituted for genome-wide discovery after results are seen.

### Null and ideal-reference comparators

- Null: use baseline profiles in both pseudo groups with the same allocation.
  All label discoveries are false discoveries for this randomized-label
  experiment. Report any-discovery frequency and counts; do not equate the
  fraction of features significant with FDR.
- PROPOSED ideal-reference arm: replace pseudo-case observations by expected
  profiler-scale post-spike profiles derived from the same baseline and dose.
  Keep allocation, model and BH family identical. This tests discovery under
  the stipulated ideal reference, not biological truth. Validate the complete
  feature-table construction, including compositional changes and zeroes,
  before using it. Its observed-versus-ideal difference must not be called an
  exact causal partition of all real-world non-replication.

No rules are learned for the basic resampling experiment, so three-cohort
comparison is not automatically a training/validation split. Any learned
prediction of the reproducibility boundary must train on other cohorts and be
evaluated in the held-out cohort with fixed thresholds.

### Repetition and uncertainty

PROPOSED 1,000 allocations per feasible cell after a small implementation/runtime
pilot. Use fixed seeds and a predeclared numerical precision/runtime rule, not
stopping based on favorable significance. Conditional discovery frequency is
the fraction of allocations recovering a target. Paired cross-cohort joint
frequency uses independently generated cohort allocations.

More allocations improve Monte Carlo precision, not the number of independent
biological samples. Label Monte Carlo error separately from biological
uncertainty. Freeze an outer biological-sample resampling scheme if population
intervals are claimed; prohibit the same original sample appearing in both
groups under that scheme. If a defensible scheme is not implemented, report
conditional frequencies explicitly and cohort-specific variation, without
presenting binomial intervals over allocations as population confidence limits.

### Boundary interpretation

The primary result is a frequency surface over tested dose and n, with all
targets visible. Do not impose a universal threshold or assume monotonicity.
If an 80% or 90% reproducibility contour is desired, choose and label it as an
operational summary before inspecting the new outcomes. Report unbracketed
boundaries as below/above the tested range or not established. No extrapolation
beyond available samples or dose support.

Do not convert a native disease log-fold change directly to an added read
fraction. Any overlay needs a baseline-aware measurement-scale mapping and a
validated interpretation. Until then native associations and the experimental
surface are presented alongside each other, without assigning each native
marker a biological truth category. Sequencing depth is observed background
variation here; no read-depth intervention or universal required depth is inferred.

## Implementation and acceptance plan

| Work item | Type | Acceptance gate |
|---|---|---|
| Reconcile scientific protocol with frozen policy and old model documents | Design/documentation | Detailed new choices signed off and dated; old results stay development |
| Resolve endpoint-directory conflict and stale quantitative consumers | Development | Real caller/helper fixture; selected-scale results verified; old reference explicitly sensitivity |
| Extend preflight to target/effective-genome inputs | Development + cluster verification | All three final sample panels covered; exact input hashes |
| Audit corrected references across cohorts | Execution | Feng/Zeller replication of scale checks; full sample/target coverage |
| Add artificial discovery/replication evaluator | Development | Hand-checked both/one/neither, zero denominator, missing-versus-zero and unequal feature families |
| Add synthetic allocation manifest generator | Development | Distinct samples per experiment, reproducible nested allocations, identical profiler/dose pairing |
| Add complete pseudo-group model and null/ideal comparators | Development | Known synthetic positive/null cases, full BH refits, explicit degenerate-fit handling |
| Add conditional-frequency aggregation and uncertainty | Design + development | No pseudoreplicated population intervals; denominators and Monte Carlo precision preserved |
| Execute agreed modules on final three-cohort inputs | Execution | Diagnostics, sample panels, source-data provenance and stage seals pass |
| Build common plots and final manuscript evidence | Reporting | One style, full species order, units, n, uncertainty and TSV source per panel |

The six execution/documentation gaps in `DOWNSTREAM_SCOPE_AUDIT_20260929.md`
remain open. Modules removed from the main story no longer need integration as
mandatory run stages. Do not implement every historical roadmap item by default.

## Decisions to resolve at the next design review

Recommended defaults, not yet implemented/frozen:

1. Main pseudo-group population: adenoma community experiment; individual spikes
   as small-n sensitivity, other backgrounds as context.
2. Shared n grid: 5/10/15/20 per group, with eligibility review.
3. Main replication reporting: all pairwise cohort patterns; all-three secondary.
4. Model: unadjusted randomized pseudo-group log-abundance/HC3; shared baseline
   feature universe, rather than changing filters in every replicate.
5. Replicates: 1,000 after a computational pilot; uncertainty labels and optional
   outer-resampling design resolved before final interpretation.
6. Main output: full frequency surfaces; no mandatory arbitrary cutoff.
7. Comparator: baseline null required; ideal-reference arm pending construction
   validation and feasible compute requirements.

No new statistical analyses, cluster submissions, or manuscript edits were
performed in preparing this plan.
