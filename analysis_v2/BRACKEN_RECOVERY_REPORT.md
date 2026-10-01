# Three-cohort Bracken recovery comparison report

Consume the successful recovery outputs from job3097807, at
work/bracken_all_input_recovery_20260930T153115Z. Build descriptive summaries and
figures to evaluate the effect of changing the denominator. Primary DA tables
and upstream profiles are untouched. No recovery-quality conclusion is assumed.

## Summaries

Each cell is cohort x population x clinical condition x target x nominal dose.
Individual and community interventions are separate. Nominal doses identify
designed levels; achieved total and target fraction ranges are also reported.
The summary never pools species or mixes dose levels due to integer rounding.

All observations enter medians and Q1–Q3, including count/printed-fraction zeros
and negative recovered signals. Quantiles use linear interpolation (equivalent
to R quantile type7). No sample-interquartile-range bar is labelled a confidence
interval. Repeated doses/targets from the same people do not create independent
biological replicates; no inferential p-values are computed by this report.

For both methods, absolute relative error is abs(recovery_ratio-1). It uses
that method's own preserved target reference; exact versus printed fractions
may contribute small differences. Paired changes are calculated per sample
before aggregation, not as differences of independently summarized medians.
Negative all-input-minus-old error means lower error under the all-input
reference, not evidence of statistically significant improvement.

## Outputs

- tables/recovery_summary.tsv: median, IQR, min/max recovery, error IQR, n,
  zero/negative counts, and achieved-dose ranges for both methods.
- tables/paired_method_differences.tsv: paired ratio/error changes, IQR and
  counts with lower/equal/higher error.
- figures/: eighteen comparison figures, each PDF and PNG: recovery ratio,
  absolute relative error, paired error change for three cohorts x two populations.
- tables/input_hashes.tsv, session_info.txt, source_commit.txt and code/container
  hashes: provenance. SHA256SUMS covers report artifacts except the final status
  marker; SUCCESS is written only after tables and plots complete.

Figures use the same taxon/condition order and method colors. Each has ten species
columns and three condition rows. Sample medians and IQR remain visible on a
signed pseudo-log y-axis (linear close to zero, compressed in the tails), labelled
with actual values. Dose x-axis is logarithmic; community dose label explicitly
means total addition. Ratio1/error0 reference lines define ideal recovery.
Figures are inspection outputs; publication panel selection should follow review
of the real results. Species-specific achieved fractions are in the source TSVs.

## Input validation

Verify source checksums, PASS_ENDPOINT_CONSTRUCTION status, cohort identity,
production sample/endpoint counts, unique target rows, finite numeric values and
valid achieved fractions. Duplicate samples within a summary cell fail. Missing
data are not silently dropped. Output must be fresh and outside the recovery
input tree. A scaled-fixture switch exists only for synthetic tests, not launcher use.

## Cluster run

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
python3 -B analysis_v2/tests/test_bracken_recovery_report.py
export PROJECT="$PWD"
export RECOVERY_ROOT="$PWD/work/bracken_all_input_recovery_20260930T153115Z"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export REPORT_ROOT="$PWD/work/bracken_recovery_report_$(date -u +%Y%m%dT%H%M%SZ)"
REPORT_JOB="$(sbatch --parsable --export=ALL analysis_v2/run_bracken_recovery_report.sbatch)"
printf 'Job: %s\nReport: %s\n' "$REPORT_JOB" "$REPORT_ROOT"
```

Submit only after tests report OK. The launcher checks R dependencies before
building tables and renders in the pinned analysis image (ggplot2/scales/jsonlite).
The script hashes that image for provenance. Check accounting and stdout/stderr,
then the report checksum file. Synthetic fixtures verify rendering locally;
only real cluster output can establish the observed denominator effects.
