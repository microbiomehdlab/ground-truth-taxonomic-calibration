# Reference assumptions and scope decision — 5 October 2026

Update: the half/double conversion stress is executed and verified; see
REFERENCE_SCALING_REVIEW_20261005.md. The low-dose contrast persists, but Bfrag,
Csym and Hhat expected detectability is scale-sensitive. Construction sensitivity
for this declared global factor is no longer pending. Next recommended scope:
all 100 saved allocations with the same grid and uniform arm.

Reviewed expected_profile() and preparation in reference_response_pilot.py,
existing formula fixtures, verified pilot outputs and verified NE audit. This is
a conceptual/code review, not independent validation of true cell abundances.

## What is defensible

The comparator uses the same saved people/group assignments, existing baseline
profiles, frozen families, transform and full-family correction. Integer inserted
counts and achieved fractions are checked. Expected native profiles are built
before family selection. Observed primary fits are reused. Analytical zeros and
perfect fixed fits are not conflated. Thus it is a controlled, baseline-anchored
expected-response comparison, useful to locate conditional gaps in this design.

## Important assumptions that remain

1. Bracken uses S+N, with S the estimated baseline species count total, not R+N.
   It assumes baseline estimated background counts remain unchanged and every
   inserted pair contributes to the intended target. It is native-scale response,
   not the all-input recovery estimand. Unclassified/non-species observations do
   not enter S. All-input recovery is separate evidence and must not be substituted
   into this DA comparison without changing the abundance estimand.
2. MetaPhlAn genome-equivalent conversion uses baseline G_eff and target genome
   size. It approximates a relative cell/genome-equivalent response, not marker
   breadth, copy number, alignment, genome representation or perfect measurement.
   Large low-dose gaps demonstrate deviation from this stipulated comparator,
   not an independently measured error in true cell abundance.
3. Both references retain baseline profiling bias. Between-tool references are
   not one truth scale. A reference-only result is not a causal fraction of absent
   natural clinical biomarkers, nor a correction to apply to disease abundance.
4. The log2(1+a/1e-8) transform magnifies changes near zero. MetaPhlAn zero versus
   positive reference can therefore create large coefficient differences; these
   are not disease fold changes. Retain native recovery/zero evidence alongside
   effects and existing pseudocount sensitivities. No pseudocount selected by
   desired significance.
5. The reference modifies background abundances through dilution and can alter
   the distribution of non-target p-values. Full-family BH is correct for the
   declared procedure, but reference/observed q differences do not isolate target
   attenuation alone. Inspect target coefficients, SE and raw p as well as q.
6. Zero within-group variance after printed-output quantization can prevent a
   conventional test even with positive separation. The NE audit establishes the
   observed input pattern, not the cause of quantization. Never replace NE p=1
   bookkeeping with a valid nonsignificant test or infer absence of signal.
7. Three doses and 20 allocations are an explicitly post-hoc pilot, not clinical
   population power or dose thresholds. No independent disease/control comparison
   is introduced by resampling adenoma backgrounds into artificial groups.

## Scope decision

Implementation update: REFERENCE_SCALING_SENSITIVITY.md prespecifies half/double
global MetaPhlAn conversion stress scenarios with fixed people/doses/models. These
are broad assumption probes, not empirically estimated uncertainty bounds. The
original pilot is preserved. Scientific review of these outputs precedes the
proposed 100-allocation expansion below.

The completed pilot is suitable for descriptive reporting with these caveats.
Do not launch a giant reference expansion simply because the contrast is large.
Next validate construction-level sensitivity using saved values/designs, focusing
on MetaPhlAn G_eff/reference scaling and the low-dose zero boundary; record all
parameter choices before examining new significance results. Such sensitivity
is an assumption check, not optimization of detection.

After that review, the efficient first expansion is the same declared grid and
uniform arm U with all 100 existing allocations (3600 contexts, 360 ten-context
batches). Adding heterogeneous U/P25/P50/P75 arms or other n/doses is a distinct
scope extension, not automatically authorized by completion of the pilot.
Current gated submission scripts remain 720-context pilot scripts; changing only
an environment variable will not expand them. No expanded jobs were submitted.

## Compact figure draft

Script: plot_reference_response_compact.R. Local output outside Git:
`../reference_response_compact_20261005_v2/`. Four panels:

- A: observed/reference positive significance versus tested dose, n=20/group,
  pooled across ten targets and three cohorts (600 comparisons/point).
- B/C: target-specific observed/reference frequency at the lowest dose for each
  profiler (60 comparisons/target), preserving heterogeneity.
- D: separate verified NE states at the lowest dose: Bracken perfect fits versus
  MetaPhlAn all-zero input (32/600 versus 112/600).

The source tables, caption, input hashes and private style snapshot are retained.
The pooled conditional frequencies are descriptive, not independent patient
rates. Complete cohort-specific results remain in the 30-page review booklet.
This is a draft figure, not silently integrated into the manuscript or submitted.
