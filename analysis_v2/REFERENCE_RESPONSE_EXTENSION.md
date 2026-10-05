# Reference-response extension: implementation and execution

This is a new, explicitly post-hoc study extension. It does not replace the
completed observed DA3 fits or authorize a final production analysis. No new
FASTQs or taxonomic profiles are generated. Interpretation must not equate a
reference-only discovery with a causal explanation of natural adenoma disease.

## Three deliverables

1. A baseline-anchored reference-response DA3 pilot: all ten targets, three
   cohorts, both profilers, Adenoma backgrounds, uniform exposure, total
   community doses 0.0001/0.001/0.01, n=10/20 per group, the first 20 saved
   allocation IDs in lexical order. There are 720 contexts and 7,200 target
   comparisons. The total community doses correspond approximately to
   0.001%/0.01%/0.1% per target; exact integer counts are used in construction.
2. Conditional cross-cohort consistency and exact-allocation profiler agreement
   from all 120,000 existing observed contexts; no new models.
3. A matched native-response assembly comparison for Pana/Pint: 720 paired
   tool-level observations in 30 people. No new DA models or legacy recovery
   ratios. Assembly-specific genome-corrected recovery remains separate work.

## Reference contract

All native background features are diluted before the original frozen DA family
is selected. No new prevalence filtering, normalization, target-only BH family,
or outcome-based selection is introduced. Comparisons retain effect, standard
error, p, full-family q, estimability and numerical sensitivity.

Bracken: let S be the sum of baseline estimated species counts, N the total
inserted pairs and Nt the target-specific inserted pairs. The baseline-native
anchored expectation for feature j is `b_j*S/(S+N) + Nt/(S+N)` for an implanted
target, and `b_j*S/(S+N)` for every other feature. This is deliberately NOT the
all-input recovery denominator R+N. Baseline printed fractions are retained:
we do not replace baseline native values with reconstructed counts or silently
renormalize their rounding deficit. The assumption is unchanged background
estimated counts and complete assignment of the inserted signal. It is a
stipulated reference, not the internals of a perfect Bracken run.

MetaPhlAn: F=N/(R+N), f_t=Nt/(R+N), q_t=f_t*G_eff/G_t and
D=(1-F)+sum(q_t). Expected abundance is `b_j*(1-F)/D + q_t/D` for targets,
and `b_j*(1-F)/D` otherwise. G_eff comes from the verified sealed-baseline
cov95 bundle. This approximation does not model marker breadth or establish
true cell abundance; baseline measurement bias and genome-size assumptions
remain. The two workflows do not share an independently measured truth scale.

Canonical inserted counts must reconcile exactly to the recorded community
design N_total. Canonical total fractions must reproduce the recorded design
f_hat (tolerance 1e-12); historical f_hat is rounded to eight decimals and must
agree with N/(R+N) within 5.01e-9. Canonical community target fractions retain
their exact integer-count definition (tolerance 1e-12, not the rounded-total
tolerance). Reference calculations always use exact N/(R+N), and both exact and
recorded totals plus their difference are audited. Genome-size and source
profile checksums, ten-target coverage, baseline identities and observed-plan
identity are required; no fallback, interpolation or silently dropped sample.

Reference fits use the existing MaAsLin2 1.18.0 group-only LM wrapper, NONE
normalization, `log2(1+a/1e-8)`, and full-family BH with non-estimable p=1.
The observed primary results are reused, not refitted. Group-only HC3 sandwich
standard errors with a residual-df t approximation are computed on both arms.
This is an inferential sensitivity, NOT a proof of calibration under arbitrary
skew and NOT randomization inference. The reconstructed observed ordinary
coefficients/SE/p/q must agree numerically with saved MaAsLin2 target results.

Categories: BOTH_POSITIVE, REFERENCE_ONLY, OBSERVED_ONLY, NEITHER_POSITIVE,
OBSERVED_NON_ESTIMABLE, REFERENCE_NON_ESTIMABLE, BOTH_NON_ESTIMABLE. The first
four apply only if both primary fits are estimable. Separate HC3 categories
are retained. No category is labeled a causal fraction of biomarker failure.

## Cluster: one gated engineering submission

First transfer/publish the new code to the cluster repository. These changes
are not automatically committed or pushed by local implementation.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
export PROJECT="$PWD"
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif
export DA3_BATCH_PLAN="$PWD/work/preproduction_checks_20261002T192903Z/draft_da3_plan"
export DA_INVENTORY_ROOT="$PWD/work/three_cohort_da_inventory_20261001T151231Z"
export GEFF_ROOT="$PWD/work/three_cohort_sealed_geff_20261001T135712Z"
export DA3_OBSERVED_REPORT="$PWD/work/da3_direct_maaslin_20261003T170501Z/REPORT_20261003T170503Z"
export REFERENCE_ROOT="$PWD/work/reference_response_pilot_$(date -u +%Y%m%dT%H%M%SZ)"
export REFERENCE_CONCURRENCY=30
bash analysis_v2/submit_reference_extension.sh
```

Preparation tests the REAL pinned MaAsLin2 backend in the analysis image,
checks the input contracts and snapshots the calculation code and prepared
native/reference values. Only then may 72 ten-context array tasks run, with up
to 30 concurrently. Collection depends on every task succeeding. Later edits
to the working repository do not alter already prepared worker code. Do not
edit source during the preparation stage itself. Preparation failure leaves
dependent jobs pending: inspect reference_prepare.JOB.err rather than removing
the dependency. jobs are recorded outside the immutable plan at
`$REFERENCE_ROOT.jobs.tsv` so they can be recovered after reconnecting.

```bash
cat "$REFERENCE_ROOT.jobs.tsv"
jobs="$(awk -F '\t' '{sub(/;.*/,"",$2); print $2}' "$REFERENCE_ROOT.jobs.tsv" | paste -sd, -)"
sacct -X -j "$jobs" --format=JobID,State,Elapsed,ExitCode
cat "$REFERENCE_ROOT/REPORT/status.json"
(cd "$REFERENCE_ROOT/REPORT" && sha256sum -c --quiet SHA256SUMS)
```

Workers use atomic publication, locks and image/plan identity checks. Failed
scratch attempts remain for diagnosis. Successful scratch models/matrices are
removed only after compressed full-family statistics, target statistics and
checksums have been published. Completed tasks can be retried with the same
array indices; valid batches resume without fitting again. The convenience
submission refuses an existing run/ledger to prevent duplicate submission.

## Existing-output cohort summary

Can run locally against the downloaded report, or on the cluster:

```bash
python3 -B analysis_v2/scripts/summarize_da3_consistency.py \
  --report "$DA3_OBSERVED_REPORT" \
  --out "$PWD/work/da3_conditional_consistency_$(date -u +%Y%m%dT%H%M%SZ)"
```

Tables include cohort frequencies, exact-allocation profiler agreement and
cross-cohort summaries. Joint cohort frequencies are products of empirical
frequencies under independently drawn allocations from each cohort's saved
allocation distribution. Allocation numbers across cohorts are NOT paired
people. They are not prospective clinical replication rates or population
power. A zero discovery denominator gives NA, not perfect/zero conditional
replication. All three pairwise cohort comparisons are retained, not selected
after observing favorable results. Directional positive evidence uses the
existing positive-effect/full-family q<=0.05 criterion and unequal cohort
feature families remain visible.

## Assembly native comparison

Use the existing builder to select both original and replacement profiles:

```bash
export ASSEMBLY_COMPARISON_ROOT="$PWD/work/assembly_native_comparison_$(date -u +%Y%m%dT%H%M%SZ)"
python3 -B analysis_v2/scripts/build_assembly_sensitivity_input.py \
  --manifest "$PWD/work/yachida_strict_final_20260823/state/production_seal/independent_10_per_condition.tsv" \
  --arms datasets/yachida/assembly_sensitivity_arms.tsv \
  --spike-panel spikes/spike_panel.tsv --aliases examples/spike_taxon_aliases.csv \
  --sensitivity-root "$PWD/work/yachida_assembly_sensitivity_20260901" \
  --original-root "$PWD/work/yachida_strict_final_20260823" \
  --baseline-root "$PWD/work/yachida_strict_final_20260823" \
  --outdir "$ASSEMBLY_COMPARISON_ROOT/canonical"
python3 -B analysis_v2/scripts/compare_assembly_native_response.py \
  --canonical "$ASSEMBLY_COMPARISON_ROOT/canonical/canonical_input.tsv" \
  --aliases examples/spike_taxon_aliases.csv \
  --out "$ASSEMBLY_COMPARISON_ROOT/REPORT"
```

This checks native values against profiles, identical baseline, exact achieved
dose, mate counts, original pair count, inserted count and seed. Missing or
unmatched data fail, not silently shrink the panel. It produces per-person
differences and condition-specific descriptive summaries, not genome-corrected
recovery or a causal contamination claim. The output explicitly remains pending
source-seal/receipt verification and scientific review before manuscript use.

## Validation and remaining work

Cluster preparation 3108871 (5 October) passed the actual pinned MaAsLin2
fixture but exposed an incorrect new 1e-10 total-fraction tolerance against
historical eight-decimal f_hat. The reconciliation contract above corrects that
interface bug without changing data or reference formulas. Unit tests reject
out-of-rounding design differences, canonical/design disagreements, rounded
community target fractions and incorrect integer counts. The complete pipeline
fixture now uses eight-decimal totals rather than unrealistically exact totals.
Failed preparation leaves the dependent array/collector unstarted; cancel those
pending jobs and retry in a fresh run directory after pulling the correction.

Local: Python formula/design tests, complete prepare/run/resume/collect fixture
with a mocked R process, R full-family/backend-interface fixture with an LM test
double, assembly pairing gates and shell syntax. The actual pinned backend is
NOT available locally; the cluster preparation job tests it and fails closed.
The mock tests are not represented as real MaAsLin2 execution.

### Collector-only recovery (2026-10-05)

Run `reference_response_pilot_20261005T202107Z` completed preparation (3108874)
and all 72 analysis tasks (3108875), but collector 3108876 rejected its own
designated `REPORT` child as overlapping output. Collection now permits only
that designated child inside the run, still refuses existing outputs and other
nested destinations, and records the collector hash in the sealed report.
The collection sbatch uses the current `PROJECT` checkout; the immutable
calculation snapshot and completed batches remain unchanged. After pulling the
fix, export `PROJECT` and the existing `REFERENCE_ROOT`, then submit only
`analysis_v2/collect_reference_response.sbatch`. Do not resubmit preparation or
the analysis array. Regression checks cover the actual nested report path,
unsafe destinations, checksum preservation and refusal to overwrite a report.

Next, without selecting settings by favorable significance:

- Review all taxa's matched primary/HC3 effects and categories, including every
  non-estimable case; validate the complete-family reference construction.
- Run 100 saved allocations and U/P25/P50/P75 once construction/operational
  gates are sound. The CLI already supports these plans; this is a scientific
  scope expansion and is not automatically launched. Use --batch-size 10.
- Add guide-consistent figure scripts for reference/observed effects and
  frequencies after auditing outputs. Do not update manuscript conclusions yet.
- Verify assembly receipt/source coverage and prepare assembly-specific target
  genome sizes for corrected quantitative recovery (especially MetaPhlAn).
- Inventory retained Kraken per-read classifications and MetaPhlAn alignment
  files; define a limited read-pool/marker diagnostic before any mixed-library
  reruns. Do not reprofile the entire experiment.
- Design robust assignment-based inference only with an explicit sharp/weak
  null, valid exchangeability assumptions and adequate p-value resolution.
  HC3 alone is not the end of this work; avoid an open-ended DA-method search.
- Defer new strains and exact patient-specific CRC-informed doses until the
  existing analyses identify a specific unresolved question.
