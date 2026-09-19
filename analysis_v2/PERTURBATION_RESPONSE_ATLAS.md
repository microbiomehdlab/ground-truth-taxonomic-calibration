# Target-agnostic perturbation-response atlas

This analysis treats controlled spike-ins as interventions used to characterize
the measurement system, rather than only as a calibration standard for named
taxa. For each profiler separately, it estimates a target-to-reported-feature
response operator from independent perturbations and asks whether those learned
responses generalize to mixtures, unseen targets, and unseen cohorts.

The exact composition counterfactual for feature *j* is
`(1 - F) * baseline_j + implanted_j`. The response signal is observed abundance
minus the dilution-retained baseline. Non-reported taxa are structural zeroes;
the implementation restores these sparse cells in every denominator.

Reliability certificates exclude every perturbation in which the evaluated
feature was directly implanted. They combine three equally weighted components:
resistance to unexpected appearance, resistance to dropout, and quantitative
stability when both expected and observed abundance are positive. Component
rates use Jeffreys Beta(0.5, 0.5) shrinkage. The score is explicitly a
measurement-reliability index, not a probability that a taxon is biologically
real and not an empirical false-positive rate, because no sham reprocessing arm
is available.

Validation is deliberately target-agnostic: leave-one-cohort-out tests transport,
leave-one-target-and-cohort-out tests generalization beyond a priori calibration
taxa, and panel-size saturation tests whether a small sentinel panel preserves
the full-panel feature-risk ranking. Community superposition uses ordinal dose
matching, not floating-point equality, and is interpreted as an ensemble system
response because independent and community read draws are not byte-identical.

Legacy Feng and Zeller fractions are nominal; Yachida fractions are derived from
integer read allocation. Therefore current mixed-source output must remain
`DEVELOPMENT_ONLY` until the analysis and display decisions are frozen and rerun
on the definitive cohort inputs.
