# Calibration-aware false-discovery and robustness analysis

**Status:** prospective extension for the revised analysis. The Feng and Zeller
legacy outputs may be used for development. No Yachida filtering result is
confirmatory unless a timestamped rule lock predates first inspection of the
relevant Yachida outputs; otherwise it is exploratory. This analysis was not
prespecified in the original preprint.

## Scientific distinction

In the paired artificial contrast, the implanted taxon is the known positive
and every significantly enriched non-target is a ground-truth false discovery.
In a native CRC contrast, biological truth is unknown. A call that disappears
after implantation is therefore **perturbation-sensitive**, not necessarily
false: it may cross the q-value boundary through effect attenuation, changed
uncertainty, multiplicity, zeros, or profiler cross-talk.

## A. Leakage-safe artificial-call filter

For each profiler and taxon, learn an artefact-susceptibility score from
non-target artificial contexts in training cohorts only. A context in which the
taxon is the intended target is excluded from score estimation. Correlated dose
levels are nested within implanted target; implanted targets, not individual
doses, are the independent recurrence units.

Use inner leave-one-target-out validation to construct a Pareto curve of
held-out off-target discovery rate against held-out target recall. The rule is
selected in training data subject to target recall of at least 95%. Apply the
frozen rule to the held-out cohort without protecting its true targets. Refit
the complete abundance model and BH family after prefiltering; deleting rows
from an already significant result table is not validation.

Primary validation is leave-one-cohort-out. The planned decisive test learns
the rule with Feng and Zeller, records the lock, and applies it once to Yachida
only if the relevant filtering outcomes have not already been inspected. If
that condition is not met, Yachida validation is labelled exploratory. Report absolute feature removal, eligible
non-target feature tests, off-target discovery rate, precision, F1, target
recall, target-effect bias, and uncertainty. Report all threshold sensitivity,
including an unfiltered comparator.

`scripts/screen_cross_cohort_artifact_filter.py` provides a deliberately
conservative, development-only post-fit falsification screen. It learns
recurrence across implanted targets (not correlated doses), transfers the score
between Feng and Zeller, and never protects the held-out implanted target. Its
automatic choice table enforces at least 95% relative target-recall retention.
It cannot support a manuscript filtering claim because the abundance model and
complete BH family have not yet been refitted after prefiltering.
It writes `heldout_context_metrics.tsv`, `heldout_filter_summary.tsv`,
`recall_preserving_choices.tsv`, and `training_artifact_scores.tsv`.
`scripts/plot_cross_cohort_artifact_filter_screen.R` adds a separately
checksummed Pareto display without changing the screen's development-only
status.

### Development result: a static blacklist is not supported

The blind Feng-to-Zeller and Zeller-to-Feng screen rejected a general static
taxon blacklist. For both profilers in the community experiment and for
Kraken2 + Bracken in the independent experiment, the only rule satisfying the
95% relative-recall constraint was no filtering. The one sizeable exception
was Zeller-to-Feng MetaPhlAn 4 in the independent experiment: a recurrence
threshold of 0.10 removed 34.7% of held-out off-target calls while retaining
all target calls. The reverse transfer removed only 0.6%. This asymmetric,
post-fit exception is insufficient to freeze a rule. The negative result is
retained because it shows why the earlier target-protected exclusion appeared
more successful than an honest blind application.

## B. Native disease-effect distortion

For original abundance `o`, total implanted fraction `F`, target implanted
fraction `f`, and observed spiked abundance `a`, define the ideal
read-proportional reference

- non-target: `e = (1 - F)o`;
- intended target: `e = (1 - F)o + f`.

With the frozen `1e-8` fraction pseudocount, compute

`u = log2(a + 1e-8) - log2(e + 1e-8)`.

For every feature and perturbation context, fit the same frozen disease model
to `u`, with CRC versus Control primary and age/sex adjustment. The condition
coefficient tests phenotype-dependent departure from the read-proportional
reference: a paired difference-in-differences estimand.
Use HC3 uncertainty and BH correction over the frozen feature universe.
`IDEAL_COUNTERFACTUAL_DID.md`, `scripts/fit_ideal_counterfactual_did.R`, and
`run_ideal_counterfactual_did.sh` implement this model and its fail-closed
execution contract.

The disease transition ledger records observed baseline/dose effects, q values,
signs, and call transitions. The separate ideal-reference result table records
the residual effect with HC3 uncertainty and context-specific BH correction.
A future joined table may classify calls as:

1. stable within a prespecified equivalence margin;
2. q-threshold-only loss with effect change inside that margin;
3. attenuated or amplified;
4. direction reversed;
5. significant condition-dependent profiler distortion;
6. newly gained.

Until an equivalence margin and joined classifier are frozen, a nonsignificant
residual test means only “no detected departure”; it is not evidence of no
distortion.

The implanted target is reported separately. Target addition can mechanically
compress its pre-existing disease contrast, so the primary native-signature
robustness endpoint is target-excluded bystander distortion.

## C. Does filtering improve disease inference?

Apply the artificial rule learned in other cohorts before disease-model fitting
in the held-out cohort. Label removed calls `calibration-sensitive`, not false
positive. A useful filter should reduce label-permutation discoveries and
increase directional cross-cohort replication while preserving strong effects.
Perturbation stability and the difference-in-differences residual are secondary
validation endpoints and cannot be used to train the same rule.

### Implemented abundance-level calibration experiment

`scripts/calibrate_abundance_cross_cohort.py` learns only positive, recurrent
departures from the ideal dilution reference in training cohorts and subtracts
their predicted contribution in a held-out cohort. Independent spikes retain
target-specific coefficients. Community spikes are calibrated as one joint
mixture because their component fractions are not independently identifiable.
Direct target features are protected, negative or weakly supported corrections
are set to zero, and every altered abundance is bounded below by zero.

The script emits uncorrected and corrected disease-model inputs with identical
held-out contexts. The existing frozen disease model is then completely refit,
including HC3 uncertainty and each BH family. `evaluate_calibration_restoration.py`
compares both fits with the unchanged unspiked baseline and reports abundance
error restoration, baseline biomarkers rescued or harmed, induced calls
removed or created, and disease-effect restoration. The complete fail-closed
entry point is `run_cross_cohort_abundance_calibration.sh`; its figures are
generated by `plot_calibration_restoration.R`.

This remains an experimental calibration estimator. Positive held-out
restoration is evidence that the learned measurement response transfers; it
does not prove that every altered taxon was biologically absent or that the
same correction applies outside the tested spike design.

### Truth-anchored test of the fragility hypothesis

The statement “biomarkers that disappear are enriched for false positives” is
tested rather than assumed. Use the same permuted phenotype assignment for a
biological sample's baseline and every perturbed profile, then refit every
disease model and its complete BH family. Discoveries under these null-label
fits are known false discoveries. Raw label shuffling is used only when its
exchangeability assumptions hold; otherwise use a covariate-preserving
Freedman-Lane or wild-bootstrap null.

Contrast three prespecified groups:

1. null-label discoveries (known false discoveries);
2. directionally replicated baseline biomarkers in an external cohort
   (higher-confidence, not guaranteed true);
3. cohort-specific baseline biomarkers (unknown truth).

The primary hypothesis is that null discoveries have lower target-excluded
perturbation persistence and larger ideal-counterfactual distortion than the
externally replicated group. Estimate the complete stability-score ROC/PR
curve in training cohorts and evaluate it once in a held-out cohort. If that
separation does not replicate, disappearance cannot be used as a false-positive
filter.

For interpretation, report four evidence tiers rather than a binary truth
label: replicated and stable; replicated but perturbation-sensitive;
cohort-specific but stable; and cohort-specific plus distorted. These are
confidence/robustness tiers, not declarations of biological truth.

### Implemented perturbation-reliability score

`scripts/build_perturbation_reliability_scores.py` implements a transparent
descriptive score for each target-excluded baseline CRC biomarker. It averages
three equally weighted components: call persistence across eligible stress
tests, disease-effect direction stability, and a bounded symmetric
effect-fidelity measure. When a separately learned artificial off-target score
is supplied, off-target resistance is a fourth equally weighted component.
The score is not a posterior probability that the biomarker is true.

External directional replication is deliberately excluded from score
calculation and added as an independent evidence label. The resulting evidence
tiers distinguish replicated/stable, replicated/sensitive, cohort-specific/
stable, cohort-specific/sensitive, and directionally discordant biomarkers.
The default 0.80 stability boundary is a development visualization aid, not a
validated clinical or feature-removal threshold. The fail-closed runner
`run_perturbation_reliability_analysis.sh` currently refuses definitive status
until thresholds and held-out validation are frozen.
Supply the full disease-model result table with `--disease-results` to
distinguish taxa tested but nonsignificant in the other cohort from taxa absent
from its frozen feature universe. Without it, the sparse transition ledger
cannot support a complete external-replication classification.

## Main displays

1. Held-out specificity-versus-target-recall Pareto curve, with the frozen rule.
2. Taxon-by-profiler artefact-susceptibility atlas with cross-cohort transfer.
3. Native call-transition composition versus dose: stable, threshold-only,
   distorted, reversed and gained.
4. Null discoveries versus externally replicated biomarkers across the
   perturbation-stability score.
5. Disease-effect replication before and after the externally learned filter,
   shown only if the held-out filter succeeds.

For engineering and manuscript-layout development,
`run_disease_biomarker_robustness_figures.sh` implements the target-excluded
retention trajectories, induced-call rate, implanted-species stress-test atlas,
and descriptive transition effect-change panel. These displays inherit the
source run's status; only the ideal-counterfactual model supplies a formal test
of phenotype-dependent measurement distortion.

The previous post-hoc claim that artefact exclusion “removes false positives”
must be retained only as historical proof of concept unless this held-out,
pre-fit and fully refitted validation succeeds.
