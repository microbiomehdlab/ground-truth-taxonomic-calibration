# Preserved historical analysis boundary

The historical downstream implementation is preserved at Git commit
`15c48d552d959f6d4b5047305389a8e284e9df24`. Revised analysis development began
after this boundary. Git can reproduce the exact source even if later maintenance
is required.

Historical entry points:

- `run_original_unpaired_q010_cluster.sbatch`
- `run_publication_original_unpaired_q010.sh`
- `cluster_preflight_original_unpaired.sh`
- `scripts/00_build_spike_design.R`
- `scripts/01_compute_spike_metrics.R`
- `scripts/02_run_spike_biomarker_benchmark.R`
- `scripts/99_make_current_manuscript_figures.R`
- the manuscript plotting scripts under `scripts/`
- reusable historical helpers under `R/`

Historical generated outputs remain in their existing
`RUNS_publication_original_unpaired_q010_*` directories and are not copied into
Git. Preserve their existing input/output checksum manifests. Do not use those
directories as v2 output roots.

To inspect the exact historical source without changing the working tree:

```bash
git show 15c48d552d959f6d4b5047305389a8e284e9df24:run_publication_original_unpaired_q010.sh
git show 15c48d552d959f6d4b5047305389a8e284e9df24:scripts/02_run_spike_biomarker_benchmark.R
```

This boundary preserves the previous analysis for comparison and audit. It does
not endorse the historical unpaired model as the final inferential analysis.
