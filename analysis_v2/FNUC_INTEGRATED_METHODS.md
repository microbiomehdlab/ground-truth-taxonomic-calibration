# Integrated *F. nucleatum* adenoma measurement figure

This one-page SVG brings the existing matched-species and spike-to-biomarker audits together. It displays, for each cohort and profiler: unspiked detection in Control/Adenoma/CRC; adenoma positive-only median abundance and all-sample IQR; adjusted adenoma and CRC effects with 95% CIs and BH q; 0.001%-per-member community-spike rescue and quantitative recovery; 0.01% independent-spike recovery and detection; and the equal-community-spike changes in adenoma effect/q, CRC call, and other bystander calls. An optional paired related-*Fusobacterium* response appears only if the source audit included `--paired-features`.

The figure is an **audit of measured explanations**, not a verdict. The equal spike is a ten-species mixture, and a changed call cannot be attributed to *F. nucleatum* alone. A successful spike does not prove native abundance is measured without bias. Nonsignificant adenoma q is not biological absence. Direct CRC–adenoma contrast/power, read-level QC and marker support, and orthogonal native validation remain separate work; they are named on the figure rather than imputed.

From the repository root, use the existing successful run directories (do not pass a timestamped parent):

```bash
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
MATCHED="$ROOT/matched_crc_adenoma_species_20260922T113325Z"
AUDIT="$ROOT/fnuc_spike_biomarker_audit_20260922T144908Z"
OUTDIR="$ROOT/fnuc_integrated_methods_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_integrated_methods.py \
  --matched-dir "$MATCHED" --audit-dir "$AUDIT" --outdir "$OUTDIR"
test -s "$OUTDIR/SUCCESS" && echo "$OUTDIR/fnuc_integrated_methods.svg"
```

If you rerun the spike-to-biomarker audit with `--paired-features`, point `AUDIT` to that new directory to populate the last related-species row. This plot uses only already-computed TSVs and does not refit disease models.
