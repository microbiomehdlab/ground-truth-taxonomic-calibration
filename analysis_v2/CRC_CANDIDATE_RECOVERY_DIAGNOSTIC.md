# Why do two CRC calls differ between profilers or cohorts?

This development-only diagnostic tests whether *measurement sensitivity is a
plausible contributor* to two specific discrepancies. It does not establish
the cause of an unspiked disease association.

- Zeller *Parvimonas micra*: compare Kraken2 + Bracken with MetaPhlAn 4.
- Yachida *Dialister pneumosintes*: Kraken reports this implanted species under
  the alias *Allisonella pneumosintes*. The canonical taxon and original
  reported feature are both retained; frozen feature identifiers are not
  silently rewritten.

The top-10 CRC impact plots now display the Kraken alias as
`Dialister pneumosintes [Allisonella in Kraken]`; source TSVs retain the
original `feature` as the join key and add `display_feature`.

Run on lobo after pulling the commit:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
CALLS="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/disease/models/primary_disease_da_results.tsv"
ENDPOINTS="$PWD/work/three_cohort_recoverability_20260921T154031Z/endpoints/paired_endpoints.tsv"
OUTDIR="$PWD/work/crc_candidate_recovery_diagnostic_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/diagnose_crc_candidate_recovery.py \
  --calls "$CALLS" --endpoints "$ENDPOINTS" \
  --aliases examples/spike_taxon_aliases.csv --outdir "$OUTDIR"
echo "$OUTDIR/zeller_Pmic_diagnostic.svg"
echo "$OUTDIR/yachida_Dpne_diagnostic.svg"
```

The upper cards show each profiler's original, unspiked *community* CRC model
effect, 95% interval, BH q, and model sample count. The lower grid uses
*independent single-species spikes* at 0.01%, 0.05%, and 0.1% on the paired
samples available to **both** profilers, separately for Control and CRC.
It reports baseline-positive fraction, post-spike-positive fraction, and the
median [sample Q1, Q3] of recovered spike signal divided by implanted signal
on the selected profiler scale (read-proportional for Kraken, genome-equivalent
for MetaPhlAn). One is ideal; zeros and negative recovered signals remain in
the calculation, with no pseudocount. All source summaries and paired-sample
values are TSVs; unmatched profiler samples are counted, not silently pooled.

Poor or variable recovery at the relevant low dose would support measurement
sensitivity as a contributor. Good recovery would point toward model
uncertainty, cohort composition, or biology instead. Neither outcome proves
why a baseline CRC association differs; the community model and independent
perturbation are related but distinct estimands. These figures are not final
manuscript evidence.
