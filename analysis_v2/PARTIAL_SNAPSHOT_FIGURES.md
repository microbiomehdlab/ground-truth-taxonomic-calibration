# Figures from the finished-sample three-cohort snapshot

This development-only handoff fills candidate Figures 3–6 using the existing
finished-sample snapshot. It does **not** imply completion of Feng or Zeller
production, definitive cohort inference, or final manuscript readiness.

Snapshot already established by André: 507 sample-wide MetaPhlAn baselines
(Feng 151, Yachida 201, Zeller 155), with corrected three-cohort paired
endpoints under `three_cohort_recoverability_20260921T154031Z/endpoints`.
The response runner reselects and counts these baselines, rebuilds native
abundance and feature responses from the same canonical input, requires the
MetaPhlAn genome-equivalent reference, and explicitly passes
`--reference-scale profiler_scale` to both downstream analyses. It requires
three-cohort holdout validation. The disease runner does **not** refit the
932-MB disease model: it reevaluates its existing partial-snapshot calls with
the corrected implanted-target exclusion and makes exploratory robustness and
directional-replication plots. Every output is `DEVELOPMENT_ONLY` and has
checksums. A missing preview remains a placeholder in the book.

Only André may execute cluster commands. Run these in a suitable **compute
allocation**, not on a login node. They can take hours and each `OUTDIR` must
be new. If a stage fails, retain its directory and log; do not claim success
or point the book at an incomplete stage.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
FIGURE_ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
SOURCE_RUN="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010"
ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
CANONICAL_INPUT="$PWD/work/analysis_v2_three_cohort_input_dev_20260913_180614/combined/canonical_input.tsv"
ENDPOINTS="$FIGURE_ROOT/endpoints/paired_endpoints.tsv"
TARGET_GENOME_SIZES="$PWD/work/target_genome_sizes_20260919T224340Z/target_genome_sizes.tsv"
EFFECTIVE_GENOME_SIZES="$PWD/work/three_cohort_geff_audit_20260921T153427Z/geff_primary_cov95/effective_genome_size.tsv"
EXPECTED_BASELINES=507
export FIGURE_ROOT SOURCE_RUN ANALYSIS_SIF CANONICAL_INPUT ENDPOINTS
export TARGET_GENOME_SIZES EFFECTIVE_GENOME_SIZES EXPECTED_BASELINES

OUTDIR="$FIGURE_ROOT/partial_snapshot_response" \
  bash analysis_v2/run_partial_snapshot_response_figures.sh
test -s "$FIGURE_ROOT/partial_snapshot_response/SUCCESS"

OUTDIR="$FIGURE_ROOT/partial_snapshot_biomarker" \
  bash analysis_v2/run_partial_snapshot_biomarker_figures.sh
test -s "$FIGURE_ROOT/partial_snapshot_biomarker/SUCCESS"
```

The biomarker runner is independent of the long response analyzer and can run
in a separate compute allocation while that analyzer continues. It uses the
existing disease-model calls and does not switch the response reference scale.
For a talk-sized, feature-level stress-test plot, after its evaluation stage
passes, run:

```bash
python3 analysis_v2/scripts/plot_biomarker_stress_test.py \
  --ledger "$FIGURE_ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv" \
  --panel spikes/spike_panel.tsv \
  --cohort yachida --dose-percent 0.1 --top 10 \
  --outdir "$FIGURE_ROOT/biomarker_stress_test_yachida_0p1"
```

The SVG selects the top ten baseline CRC calls separately per profiler,
ranked **only** by unspiked BH q. Each independent target is a column; the
direct target is excluded from its own challenge. Green means same-direction
significance retained, gold means significance lost without reversal, and
magenta means effect direction reversed. This is technical robustness, not a
false-positive verdict or biological validation. The script requires a
complete selected target grid and refuses to overwrite an existing output.

Then assemble a **new** provisional book using
`analysis_v2/PROVISIONAL_FIGURE_BOOK.md`. The inventory points to four response
previews (Figure 3A/B and 4A/B), two target-excluded biomarker previews
(Figure 5A/B), a baseline directional-replication matrix (Figure 6A), and an
exploratory reliability-versus-replication panel (Figure 6B). Figure 5B is
**not** a substitute for the continuous model or DiD; both remain pending.

Scientific limitations for review: incomplete sample accrual can bias
prevalence and case/control composition; the disease model is an existing
partial-snapshot fit rather than a definitive rerun; cross-cohort validation
uses available and evaluable features, never treats a missing feature as a
negative result; and the response atlas depicts profiler behavior, not
ecological interactions. Inspect cohort/sample counts and source TSVs before
interpreting patterns.
