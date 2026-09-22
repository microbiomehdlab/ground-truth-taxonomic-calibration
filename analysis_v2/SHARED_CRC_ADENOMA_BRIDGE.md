# Shared CRC candidates: do they appear in adenomas?

This provisional two-stage analysis starts with the **same-direction,
BH-significant CRC-vs-Control candidates in all three cohorts** separately
for Kraken2 + Bracken and MetaPhlAn 4. It deliberately does **not** select on
the Adenoma-vs-Control q-value: a weak or uncertain adenoma effect is precisely
what we want to inspect.

For each profiler it writes:

- `*_crc_anchored_adenoma.svg`: original unspiked CRC and adenoma effect,
  95% confidence interval, BH q, and model n in Feng, Yachida, and Zeller.
  Adenoma cells distinguish a same-direction call, same-direction non-call,
  opposite direction, and unavailable fit.
- `*_adenoma_measurement.svg`: full-cohort unspiked Control/Adenoma/CRC
  detectability and positive-only median [Q1,Q3] abundance. For species
  individually implanted in the ten-taxon panel, each cohort also shows
  0.01% direct-spike detection and recovery median [Q1,Q3] separately in
  Control/Adenoma/CRC. Non-panel species say **not directly spiked**.

The model, full-baseline sample, summary, direct-spike sample and summary TSVs
are written alongside the SVGs. The direct-spike TSV contains all 0.01%,
0.05%, and 0.10% doses. Source checksums and a `DEVELOPMENT_ONLY` marker are
included. This analysis does not refit disease models or rerun the expensive
feature-response analysis.

Run on lobo in a compute allocation, after pulling the relevant commit:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
SOURCE="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010"
OUTDIR="$ROOT/shared_crc_adenoma_bridge_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_shared_crc_adenoma_bridge.py \
  --calls "$SOURCE/models/disease/models/primary_disease_da_results.tsv" \
  --manifest "$ROOT/partial_snapshot_response/native_abundance/biomarker_profile_manifest.tsv" \
  --abundance "$ROOT/partial_snapshot_response/native_abundance/biomarker_abundance_long.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --panel spikes/spike_panel.tsv \
  --aliases examples/spike_taxon_aliases.csv \
  --outdir "$OUTDIR"
cat "$OUTDIR/candidate_summary.tsv"
echo "$OUTDIR/kraken2_bracken_crc_anchored_adenoma.svg"
echo "$OUTDIR/metaphlan4_crc_anchored_adenoma.svg"
echo "$OUTDIR/kraken2_bracken_adenoma_measurement.svg"
echo "$OUTDIR/metaphlan4_adenoma_measurement.svg"
```

Scientific guardrails: baseline prevalence and positive-only abundance answer
different questions; a missing taxon from a profile is treated as zero, and
zeros remain in detection counts. The full-cohort baseline sample set is
required to match between profilers, but model n may differ because of
covariate/evaluability filters. Direct recovery is measured only in the
available paired-spike subset and is on each profiler's own reference scale
(read-proportional Kraken, genome-equivalent MetaPhlAn). Non-panel taxa have
no direct recovery estimate. A non-significant adenoma q-value is not proof
of biological absence, and a low-dose recovery deficit is only evidence of a
plausible measurement limitation—not a causal explanation of disease biology.
