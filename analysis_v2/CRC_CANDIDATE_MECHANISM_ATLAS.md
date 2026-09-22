# Focused CRC candidate mechanism atlas (development only)

This is the focused follow-up to the broad replication–reliability plot. It
compares the two specific discordant calls: Yachida *Dialister pneumosintes*
(reported as *Allisonella pneumosintes* by Kraken) and Zeller *Parvimonas
micra*. Each SVG shows both profilers and three distinct sources of evidence:

1. The original unspiked, original-assembly community CRC-vs-Control model:
   signed effect, 95% CI, BH q, and model n.
2. Every available full-cohort unspiked community baseline: one dot per sample,
   zeros displayed separately, detection count and prevalence, and the median
   [Q1, Q3] abundance **among positive samples**. MetaPhlAn and Kraken values
   are native profiler-relative abundance estimates, not interchangeable
   absolute cell counts.
3. The smaller same-sample independent single-species spike subset: per-sample
   recovery ratios at 0.01%, 0.05%, and 0.10% in Control and CRC, plus median
   and IQR. The 0.01% readout explicitly gives post-spike detection and the
   fraction recovering less than half the expected signal.

Run on lobo after pulling the relevant commit. The inputs were already made
for the three-cohort partial snapshot; do not rerun the 27-million-row response
analysis. Use a **new** output directory:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
SOURCE="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010"
OUTDIR="$ROOT/crc_candidate_mechanism_atlas_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_crc_candidate_mechanisms.py \
  --calls "$SOURCE/models/disease/models/primary_disease_da_results.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --manifest "$ROOT/partial_snapshot_response/native_abundance/biomarker_profile_manifest.tsv" \
  --abundance "$ROOT/partial_snapshot_response/native_abundance/biomarker_abundance_long.tsv" \
  --aliases examples/spike_taxon_aliases.csv --outdir "$OUTDIR"
echo "$OUTDIR/yachida_Dpne_mechanism_atlas.svg"
echo "$OUTDIR/zeller_Pmic_mechanism_atlas.svg"
```

The script uses only the Python standard library. The two SVGs open in VS
Code. Auditable source tables include `full_baseline_samples.tsv`,
`full_baseline_and_low_dose_summary.tsv`, `paired_spike_samples.tsv`,
`paired_spike_summary.tsv`, and `baseline_crc_model.tsv`; input checksums are
recorded. It refuses conflicting duplicate baseline profiles and requires the
same full-baseline sample sets for both profilers within each condition.

Interpretation guardrails: the full-cohort baseline sample count and paired
spike sample count are deliberately separate. Native abundance is on a
profiler-specific scale. Low-dose recovery tests analytical behavior in
controlled samples, not the causal origin of an unspiked disease association.
The positive-only abundance IQR excludes zeros, which are displayed and
counted separately. This is exploratory, provisional evidence, not a final
manuscript claim.
