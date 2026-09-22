# Provisional cross-cohort MetaPhlAn CRC candidates

This talk preview uses the **original, unspiked community abundance fits** from
the sealed partial three-cohort disease-model run. It is separate from the
independent-spike stress-test heatmap and from the genome-equivalent correction
to expected spike abundance. It does not pool cohorts, refit models, or claim
that a candidate is biologically validated.

On lobo, after pulling the commit containing this script:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
CALLS="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/disease/models/primary_disease_da_results.tsv"
OUTDIR="$PWD/work/three_cohort_metaphlan_crc_candidates_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_three_cohort_metaphlan_crc_candidates.py \
  --calls "$CALLS" --top 12 --outdir "$OUTDIR"
cat "$OUTDIR/candidate_summary.tsv"
echo "$OUTDIR/metaphlan_crc_candidates.svg"
```

The script streams the large TSV, selects included `primary_age_sex`
`CRC_vs_Control` fits for `community`, `original` assembly, dose zero and
`metaphlan4`, then collapses repeated identical baseline fits. It fails on
contradictory duplicates. A feature enters if BH q <= 0.05 in at least one
cohort. All three cohort results, including nonsignificant and not-evaluable
cells, are preserved in `metaphlan_crc_candidates_all.tsv`. Rows are ranked by
the number of cohorts significant in the **same effect direction**, then by
three-cohort direction agreement, evaluability and best q. The SVG shows the
first 12; `candidate_summary.tsv` reports how many candidates have consistent
significant calls in two or three cohorts. Missing fits are not negative
results. Exact names must be checked before relating this preview to the
poster; the partial snapshot may have different counts.

This plot answers which baseline MetaPhlAn calls recur across cohorts. The
independent-spike Feng heatmap answers a different question: technical
stability under a known perturbation. Do not call the former spike-validated
or the latter cross-cohort replicated.
