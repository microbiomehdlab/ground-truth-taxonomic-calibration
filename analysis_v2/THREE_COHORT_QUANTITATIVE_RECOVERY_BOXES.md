# Three-cohort quantitative recovery with visible variability

The screenshot-era plot had cohort medians but its interquartile bars were hard
to discern. `scripts/plot_three_cohort_quantitative_recovery.R` is the tracked
development replacement using the **corrected** three-cohort paired endpoints.
It plots the same three community target doses (0.010%, 0.050%, 0.100%), ten
targets, three cohorts, three conditions, and two profilers. Each colored,
condition-shaped point is the sample median of `log2(observed /
expected_abundance_profiler_scale)`. A narrow translucent box behind it spans
Q1–Q3 across biological samples; the vertical outline makes a narrow IQR
visible. The plot does not imply that the box is a confidence interval.

Only original-assembly community rows enter. MetaPhlAn must have
`reference_type=genome_equivalent`, and Bracken must have
`reference_type=read_proportional`. A zero observed abundance has undefined
log2 recovery, so it is excluded from the positive-only quartiles and counted
explicitly in `quantitative_recovery_box_summary.tsv` (`n_zero`). No
pseudocount is introduced. The script requires all 540
cohort-condition-profiler-target-dose contexts, writes its input/output hashes,
and never changes the paired endpoint table.

André-only cluster command, using the already-derived corrected endpoints:

```bash
FIGURE_ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_three_cohort_quantitative_recovery.R \
    --endpoints "$FIGURE_ROOT/endpoints/paired_endpoints.tsv" \
    --outdir "$FIGURE_ROOT/quantitative_recovery_boxes"
```

This is a development figure. Review sample counts, zero counts, and plotting
limits before using it in a manuscript.
