# Integrated *F. nucleatum* adenoma measurement figure

This one-page SVG brings the existing matched-species and spike-to-biomarker audits together. It displays, for each cohort and profiler: unspiked detection in Control/Adenoma/CRC; adenoma positive-only median abundance and all-sample IQR; adjusted adenoma and CRC effects with 95% CIs and BH q; **the first significant paired spiked-vs-own-unspiked target call in the independent single-species dose series**, with its BH q in each Control, Adenoma and CRC background; two equal-community-spike disease-call tests, at 0.001% and 0.1% per member; paired rescue and quantitative recovery at both community doses; and gained/lost CRC calls among related *Fusobacterium* labels and other bystanders at both doses. Independent doses begin at 0.01%; there is no independent 0.001% dose. Optional paired related-species response is omitted from this call-focused figure.

The figure is an **audit of measured explanations**, not a verdict. The 0.1% dose is per species, so the ten-species mixture is approximately **1% total**; it is a strong perturbation and not an endogenous-equivalent detection limit. A changed call cannot be attributed to *F. nucleatum* alone. A successful spike does not prove native abundance is measured without bias. Nonsignificant adenoma q is not biological absence. Direct CRC–adenoma contrast/power, read-level QC and marker support, and orthogonal native validation remain separate work; they are named on the figure rather than imputed.

From the repository root, use the existing successful run directories (do not pass a timestamped parent):

```bash
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
MATCHED="$ROOT/matched_crc_adenoma_species_20260922T113325Z"
AUDIT="$ROOT/fnuc_spike_biomarker_audit_20260922T144908Z"
PAIRED="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/artificial_stratified/evaluation/biomarker_propagation_metrics.tsv"
HIGH="$ROOT/fnuc_spike_biomarker_audit_0p1_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_spike_biomarker_audit.py \
  --ledger "$ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --member-dose-percent 0.1 --outdir "$HIGH"
OUTDIR="$ROOT/fnuc_integrated_methods_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_integrated_methods.py \
  --matched-dir "$MATCHED" --audit-dir "$AUDIT" \
  --high-audit-dir "$HIGH" --paired-metrics "$PAIRED" --outdir "$OUTDIR"
test -s "$OUTDIR/SUCCESS" && echo "$OUTDIR/fnuc_integrated_methods.svg"
```

The paired row requires the phenotype-stratified independent metrics, not the pooled analysis. It reports the first tested significant independent dose and its target BH q from each condition's matched spiked-versus-own-unspiked model. This matches the heatmap's experiment but comes from the complete six-dose source table, not the heatmap itself. It is distinct from the community disease-call rows below it, and is a first significant **tested dose**, not a continuous detection limit. The script fails if any of the 18 independent contexts lacks a complete six-dose grid. This plot uses only already-computed TSVs and does not refit disease models. The 0.1% community audit validates its dose in `dose_metadata.tsv`; never substitute a 0.01% independent-spike output for this equal-community-spike disease-model comparison. Bystander call changes in the ten-member mixture are not attributable solely to *F. nucleatum* and are not proven false-positive biological biomarkers.
