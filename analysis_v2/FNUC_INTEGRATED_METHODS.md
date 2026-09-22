# Integrated *F. nucleatum* adenoma measurement figure

This one-page SVG brings the existing matched-species and spike-to-biomarker audits together. It displays, for each cohort and profiler: unspiked detection in Control/Adenoma/CRC; adenoma positive-only median abundance and all-sample IQR; adjusted adenoma and CRC effects with 95% CIs and BH q; **two equal-community-spike call tests, at 0.001% and 0.1% per member**; paired rescue and quantitative recovery at both community doses; and 0.01% independent-spike recovery/detection. An optional paired related-*Fusobacterium* response appears only if the low-dose source audit included `--paired-features`.

The figure is an **audit of measured explanations**, not a verdict. The 0.1% dose is per species, so the ten-species mixture is approximately **1% total**; it is a strong perturbation and not an endogenous-equivalent detection limit. A changed call cannot be attributed to *F. nucleatum* alone. A successful spike does not prove native abundance is measured without bias. Nonsignificant adenoma q is not biological absence. Direct CRC–adenoma contrast/power, read-level QC and marker support, and orthogonal native validation remain separate work; they are named on the figure rather than imputed.

From the repository root, use the existing successful run directories (do not pass a timestamped parent):

```bash
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
MATCHED="$ROOT/matched_crc_adenoma_species_20260922T113325Z"
AUDIT="$ROOT/fnuc_spike_biomarker_audit_20260922T144908Z"
HIGH="$ROOT/fnuc_spike_biomarker_audit_0p1_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_spike_biomarker_audit.py \
  --ledger "$ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --member-dose-percent 0.1 --outdir "$HIGH"
OUTDIR="$ROOT/fnuc_integrated_methods_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_integrated_methods.py \
  --matched-dir "$MATCHED" --audit-dir "$AUDIT" \
  --high-audit-dir "$HIGH" --outdir "$OUTDIR"
test -s "$OUTDIR/SUCCESS" && echo "$OUTDIR/fnuc_integrated_methods.svg"
```

If you rerun the low-dose spike-to-biomarker audit with `--paired-features`, point `AUDIT` to that new directory to populate the related-species row. This plot uses only already-computed TSVs and does not refit disease models. The 0.1% audit validates its dose in `dose_metadata.tsv`; never substitute a 0.01% independent-spike output for this equal-community-spike disease-model comparison.
