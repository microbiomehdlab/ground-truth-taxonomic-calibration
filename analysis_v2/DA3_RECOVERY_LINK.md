# Exact-person recovery and significance

This is the next analysis in ADENOMA_SIGNAL_EXPLANATION_20261003.md. It connects
known additions and recovered signal to statistical significance in the exact
selected people. It reuses all completed DA3 experiments; it fits no new models,
changes no significance criteria and requires no upstream processing.

Inputs are the immutable 100-allocation DA3 plan, all 4,800 completed result
batches, three-cohort Bracken all-input recovery endpoints, and the three-cohort
MetaPhlAn genome-equivalent recovery endpoints. Summary-only plot reports are
insufficient inputs because they do not identify individual recovery.

The worker checks checksum coverage, result-to-plan identity, target aliases,
physical profile paths, dose, biological identity, group membership and frozen
family size. It reconstructs each target's transformed mean difference from
the endpoint native fractions and requires agreement with the saved effect.
It retains standard parametric p/full-family BH q and separately named
permutation FWER p. Parametric SE is reconstructed from within-group transformed
variances; it is blank when the saved parametric p is non-estimable. No CI or
new p-value is manufactured from permutation outputs.

Outputs in report/:

- person_target_response.tsv.gz: unique community sample/dose/target endpoints,
  both native measurements and explicitly named quantitative recovery scale.
- context_target_explanation.tsv.gz: 1,200,000 target-context rows. Includes
  effect, p/q, selected-person baseline variation, exposure counts and continuous
  recovery summaries for the exact exposed cases and assigned doses.
- conditional_significance_frequencies.tsv: 12,000 cells, each with 100
  allocation results, standard and permutation positive-significance frequencies,
  and frequency of positive raw p<=.05 without BH significance.
- status.json, input_hashes.tsv, code_hashes.tsv, SUCCESS and SHA256SUMS.

Recovery summaries exclude unexposed cases. Zero-dose null arms have blank
recovery summaries, not a ratio of zero. Negative and zero recoveries are
preserved. Family corrections remain those saved from the full eligible
species family. No family is reduced to the ten targets.

The known allocation is recovered via context ID in the checksummed source
plan; the output does not expand every repeated person into millions of
duplicated endpoint rows. All ten targets remain present. Failures stop the
job without a completed SUCCESS; a missing endpoint is never replaced with a
background-pool median.

Interpretation is descriptive: correspondence between recovery and significance
does not isolate a causal fraction attributable to profiling. Native DA fractions
and recovery references are separate quantities. In particular, the Bracken
all-input denominator is not substituted into its native DA model. Conditional
allocation frequencies are not independent biological replications or universal
detection thresholds.

## Cluster command for André

Run only once the large exploratory DA3 collection has completed. The launcher
checks its report and refuses incomplete coverage. If DA3_RUN_ROOT is unset,
it accepts exactly one work/da3_exploratory_100_* directory; multiple candidates
require an explicit path and are never resolved by choosing the newest.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA3_BATCH_PLAN="$PWD/work/preproduction_checks_20261002T192903Z/draft_da3_plan"
export BRACKEN_RECOVERY_ROOT="$PWD/work/bracken_all_input_recovery_20260930T153115Z"
export METAPHLAN_RECOVERY_ROOT="$PWD/work/three_cohort_metaphlan_recovery_20261001T140855Z"
export RECOVERY_LINK_ROOT="$PWD/work/da3_recovery_link_$(date -u +%Y%m%dT%H%M%SZ)"
bash analysis_v2/submit_da3_recovery_link.sh
```

The launcher snapshots committed source, records the image digest and job ID,
and requests one CPU, 16 GB and a two-hour limit. That is a resource limit, not
a measured runtime prediction. This job is a streamed evidence join; it does
not repeat the expensive permutation work. Existing running snapshots are
unchanged. Read the final report's status.json and the job's stderr before
scientific interpretation.

## Verification

python3 -B analysis_v2/tests/test_da3_recovery_link.py exercises both reference
adapters, exact variable doses, partial exposure, empty null recovery,
zero/negative recovery, real checksum verification and all-cohort fixture input.
It also rejects missing batches/endpoints, corrupted checksums, changed native
effects, duplicate targets/people, dose identity errors and overlapping output.
The real cluster endpoint join remains pending until André runs this job.
