# DA2: within-person dose response

DA2 describes how a known addition changes the same person's measured profile.
It does not test disease association or select a method to force detection or
nondetection. This implementation does not modify any DA1/DA3 jobs or sources.

## Complete grid

Three cohorts × two profilers × three backgrounds (Control/Adenoma/CRC).
Community: all cohort members, seven total doses 0.0001/0.0005/0.001/0.005/0.01/
0.05/0.1 and all ten targets. Independent: the frozen ten people per condition,
ten separate targets and six doses (same list without 0.1). Every selected
person has the original and the exact spiked profile. There is no silent
intersection or dropping of samples, missing doses or targets.

The grid contains 1,206 paired contexts, 82,340 person–target responses and
2,340 target summary cells. Eighteen tasks (cohort/profiler/background) reuse
baselines within each task. No profiling or new spike generation is required.

## Primary descriptive outputs

`target_responses.tsv.gz` retains individual baseline/spiked native fractions,
signed native and transformed changes, achieved dose, reference-adjusted
recovery and the legacy recovery ratio. `dose_response_summary.tsv` reports
medians, IQRs, extrema, positive/negative/zero changes and reported prevalence.
No recovery ratio is clipped or discarded for being negative or above one.

Bracken primary recovery uses estimated counts divided by **all input pairs**,
including unclassified inputs in the denominator. MetaPhlAn primary recovery
uses the previously sealed genome-equivalent reference; its native output is
unchanged. Legacy native/read-reference ratios remain sensitivities. These
ratios have different profiler-specific reference definitions: they are not
one common absolute abundance scale. Native changes remain available separately.
Nominal target dose is community total/10 for the frozen equal-weight panel,
versus the full independent dose; achieved read fractions are also retained.

The six multipage PDFs use a common colour scheme, one page per target,
separate profiler/background panels, cohort lines and median/IQR intervals.
Community and independent experiments are separate files because participant
sets differ. IQRs describe person variation, not confidence intervals; these
plots do not establish a statistically matched community/independent contrast.
Complete underlying plot values remain in the summary table. Final manuscript
layout and visual review follow the cluster report.

`contexts.tsv` and `pairing_manifest.tsv.gz` identify the exact people and
hashed original/spiked profiles for every context. Input-bundle and task
checksum inventories retain the link to the unmodified source evidence.

## Supplementary paired inference (existing method, no new selection)

Use `paired_difference_context.R` unchanged: log2(1+native fraction/1e-8),
spiked minus original within biological person, two-sided one-sample t test,
mean/SE and nominal 95% t CI, n−1 df. Preserve constant/no-change statuses,
raw p/CI NA and p=1 only for BH bookkeeping. Apply BH across each **entire
frozen baseline family**, not just the spike targets. The independent subset
uses the same pre-existing condition-level family, not a newly selected one.
`full_family_results.tsv.gz` retains every feature. `target_results.tsv`
extracts only the actually implanted target(s) without recalculating q-values.
Non-target changes are not automatically false positives (composition matters).

Small n=10, zeros, pseudocount dependence, non-normal differences and existing
null diagnostics remain limitations. Paired significance is supplementary to
magnitude/consistency; no additional sign-flip campaign or model selection is
triggered. Old mixed-model pilot failures remain preserved, not silently
replaced or relabelled. These technical perturbations do not prove that an
adenoma disease effect exists or isolate profiling's causal contribution.

## Execution and verification

The launcher snapshots committed source and the image, submits preparation,
18 parallel tasks and an automatic collector. Each task verifies native
profile hashes and exact target endpoint/native values, retains profile
provenance and full-family results, and publishes only after success. Existing
verified tasks can be resumed; corrupt results are refused. Failures retain
attempt directories. The collector refuses incomplete contexts/families.
The input recovery bundles and DA inventory remain untouched.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA_INVENTORY_ROOT="$PWD/work/three_cohort_da_inventory_20261001T151231Z"
export BRACKEN_RECOVERY_ROOT="$PWD/work/bracken_all_input_recovery_20260930T153115Z"
export METAPHLAN_RECOVERY_ROOT="$PWD/work/three_cohort_metaphlan_recovery_20261001T140855Z"
export DA2_ROOT="$PWD/work/da2_dose_response_$(date -u +%Y%m%dT%H%M%SZ)"
export DA2_CONCURRENCY=18
bash analysis_v2/submit_da2_dose_response.sh
```

At completion, inspect `$DA2_ROOT/REPORT/status.json`, then run
`(cd "$DA2_ROOT/REPORT" && sha256sum -c --quiet SHA256SUMS)`.
PASS means computationally complete, not that every spike should be detected.
Review magnitudes and limitations alongside DA1 and DA3 before manuscript claims.
