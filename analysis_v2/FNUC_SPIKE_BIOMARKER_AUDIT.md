# Does an equal F. nucleatum community spike alter disease-biomarker calls?

`plot_fnuc_spike_biomarker_audit.py` connects the 0.001%-per-member community spike to the existing paired disease-model transition ledger. It tests the question: **if the same known *F. nucleatum* addition is made across clinical groups, is the post-spike CRC-vs-Control or Adenoma-vs-Control *F. nucleatum* call recovered, lost, or changed?** The added ten-species mixture totals about 0.01%. It does not create a new biological disease association; changes diagnose the response of the profiling-and-modeling workflow to a controlled addition.

Outputs:

- `fnuc_spike_biomarker_audit.svg`: baseline→spiked target effect, BH q and call status for both contrasts in all three cohorts and both profilers; gained/lost calls among related *Fusobacterium* species and other bystanders; paired target rescue and recovery separately in Control, Adenoma and CRC; optional related-species response summary.
- `fnuc_call_transition_detail.tsv`: source ledger rows for target and all changed/retained other features, with a `category` (`target`, `related_Fusobacterium`, `other_bystander`, `other_implanted_target`).
- `fnuc_call_transition_summary.tsv`: one row per cohort/profiler/contrast with target baseline/spiked effect and q and counts of gained/lost/retained related and other bystander calls.
- `fnuc_target_condition_response.tsv`: paired 0.001% target detection rescue and quantitative recovery by clinical group.
- `fnuc_related_species_response.tsv`: optional condition-specific response of other *Fusobacterium* labels from the large paired-feature parquet. The response is a profiler-scale abundance change, not a fraction of reads misassigned.

The transition ledger is sparse: an uncalled bystander is omitted if it is uncalled both before and after spiking. The summary therefore counts **changes in the call set**, not every tested feature. Newly significant bystander calls are **spike-induced call changes**, not proven false-positive biological biomarkers. Any change in BH q depends on the frozen feature universe and the complete disease-model pipeline, not on profiler measurement alone. The mixture adds nine other taxa; related-species changes cannot automatically be attributed specifically to *F. nucleatum*. More specific attribution needs independent single-species paired responses and/or read-level evidence.

Run from the repository root after setting `ROOT` and `SOURCE`:

```bash
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
SOURCE="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010"
OUTDIR="$ROOT/fnuc_spike_biomarker_audit_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_spike_biomarker_audit.py \
  --ledger "$ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --outdir "$OUTDIR"
```

This standard-library command already gives the target-call and related-call audit. For condition-specific related-species response, use the analysis image, which contains DuckDB:

```bash
SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
OUTDIR="$ROOT/fnuc_spike_biomarker_audit_with_bystanders_$(date -u +%Y%m%dT%H%M%SZ)"
apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" "$SIF" \
  python3 analysis_v2/scripts/plot_fnuc_spike_biomarker_audit.py \
  --ledger "$ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --paired-features "$ROOT/partial_snapshot_response/response_input/paired_feature_responses.parquet" \
  --outdir "$OUTDIR"
```

The 27-million-row paired-feature parquet may take time to scan. Both runs write to new directories and leave source data untouched. Review `SUCCESS` and the TSVs before using the SVG in a talk or paper. A `dose_called=1` result after equal spiking means the *association test* is significant under the perturbed profiles; it is **not** evidence that the added DNA is a true natural CRC biomarker.
