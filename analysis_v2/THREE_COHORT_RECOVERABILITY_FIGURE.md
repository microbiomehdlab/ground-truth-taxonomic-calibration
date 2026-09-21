# Three-cohort biomarker recoverability figure

This is a **development-only** successor to the two-cohort three-panel figure.
It is not a cosmetic extension of the old `joint_correlate_01_02.R` output:
the old recovery drivers use a read-proportional reference for MetaPhlAn.
The builder refuses that reference and requires corrected profiler-scale paired
endpoints for **all three** cohorts.

Inputs:

- stratified paired-model `biomarker_propagation_metrics.tsv`, from a sealed
  three-cohort model run; independent, original-assembly contexts only;
- one three-cohort `paired_endpoints.tsv` derived with MetaPhlAn
  `--metaphlan-reference genome_equivalent` and the frozen measured genome-size
  tables; Bracken remains read-proportional;
- frozen `spikes/spike_panel.tsv`.

The builder demands all ten targets, all three conditions, both profilers,
three cohorts, and all six independent-dose levels in **both** input tables.
Missing/duplicate contexts, wrong reference types, nonfinite measurements, or
invalid observed-detection flags fail closed. Do not concatenate endpoints
from different reference policies.

Definitions:

- **A:** first nominal dose with target enrichment at BH q ≤ 0.05 and positive
  fitted effect; `NR` means none of the six tested doses. The stratified
  biomarker evaluator supplies this binary call.
- **B:** at nominal 0.01%, biomarker strength is `-log10(target q)` for each
  target/cohort/condition/profiler. The four drivers are log10(median baseline
  fraction + 1e-6), observed target-detection fraction, median absolute
  `(recovered−implanted)/implanted` error, and IQR of the signed
  recovered/implanted ratio. Both recovery drivers use the selected
  profiler-scale signal from the paired endpoints.
- **C:** descriptive Spearman correlation between each driver and biomarker
  strength over ten targets, separately for every cohort/condition/profiler.
  Undefined correlations are missing, never set to zero. These 10-target
  correlations are exploratory, not validation or causal effects.

The three panels and combined PDF/PNG are produced by
`plot_three_cohort_recoverability.R`. The source builder writes TSVs, SHA-256
provenance, `DEVELOPMENT_ONLY.txt`, and `SUCCESS`. Local fixture tests exercise
all panels and reject a read-reference MetaPhlAn row.

Cluster execution is **André-only**. Agents must not access the cluster. Once
the corrected three-cohort endpoint table exists and its own seal is verified,
run in a compute allocation with fresh output directories:

```bash
python3 analysis_v2/scripts/build_three_cohort_recoverability_figure.py \
  --biomarkers "$STRATIFIED_BIOMARKER_METRICS" \
  --endpoints "$CORRECTED_THREE_COHORT_ENDPOINTS" \
  --panel spikes/spike_panel.tsv \
  --outdir "$FIGURE_ROOT/source"

test -s "$FIGURE_ROOT/source/SUCCESS"

apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_three_cohort_recoverability.R \
    --source-dir "$FIGURE_ROOT/source" --outdir "$FIGURE_ROOT/figures"
```

Do not run this against the older mixed-provenance map-reduce endpoint table
unless its MetaPhlAn rows actually carry `reference_type=genome_equivalent`.
The builder's reference gate will reject the old read reference. Do not label
the resulting development figures as definitive Nature Communications figures.
