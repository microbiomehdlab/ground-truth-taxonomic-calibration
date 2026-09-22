# Three-cohort CRC candidates and spike impact, both profilers

This development-only figure answers two separate questions in each cohort:
was a species a CRC-versus-Control call in the **original unspiked community
abundance**, and what happened to that call after a **physical ten-species
community spike**? It does not treat a changed call as a biological false
positive. It does not use independent single-species spikes or refit models.

The input disease calls are the same frozen partial snapshot as the poster-style
directional-replication matrix. The transition ledger is the corrected
community-mix analysis that excludes all ten implanted taxa from bystander
robustness. A baseline candidate that is itself an implanted target is shown
as *direct target excluded*, not as retained or lost.

On lobo after pulling this script:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
CALLS="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/disease/models/primary_disease_da_results.tsv"
LEDGER="$ROOT/partial_snapshot_biomarker/evaluation/disease_biomarker_transition_ledger.tsv"
OUTDIR="$ROOT/three_cohort_crc_spike_impact_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_three_cohort_crc_spike_impact.py \
  --calls "$CALLS" --ledger "$LEDGER" --dose-percent 0.1 --top 10 \
  --outdir "$OUTDIR"
cat "$OUTDIR/candidate_summary.tsv"
echo "$OUTDIR/kraken2_bracken_top_candidates.svg"
echo "$OUTDIR/metaphlan4_top_candidates.svg"
echo "$OUTDIR/kraken2_bracken_shared_candidates.svg"
echo "$OUTDIR/metaphlan4_shared_candidates.svg"
```

Each profiler has its own SVG so names remain readable. The full selected
feature-by-cohort records, including original effect, BH q, baseline call and
spike fate, are in `top_candidate_spike_fates.tsv`. Top ten is **per profiler**,
ranked by the number of same-direction significant cohorts, then all-cohort
effect direction, evaluability and strongest q; this is a descriptive display
rule, not a biomarker-validation score. Missing fits are not negatives.
Kraken's `Allisonella pneumosintes` feature is displayed as
`Dialister pneumosintes [Allisonella in Kraken]` using the validated spike
alias mapping; the original feature remains unchanged in the `feature` column
and the display label is recorded separately.

The `*_shared_candidates.svg` files are the slide-sized focused view: only
features with a significant effect in the same direction in **all three**
cohorts. For each cohort they separate the original signed effect/BH q from
the post-spike fate/BH q. The companion
`shared_candidate_spike_fates.tsv` retains the exact numbers. Post-spike q is
hidden for direct implanted targets, because their abundance was intentionally
changed; those cells cannot assess bystander robustness. The focused plots do
not imply that different profilers report identical feature identities.

The community spike is **one physical mixture per cohort and dose**, not ten
independent challenges. The selected median member dose is 0.1%; the SVG
obtains and labels the actual total mixture fraction separately for each
cohort and profiler. Totals are required to agree *within* a cohort/profiler
challenge, but need not be the same across cohorts.
Dark green means a baseline q<=0.05 call remains significant in the same
direction, gold means significance is lost, magenta means a significant call
reverses direction, blue means a previously nonsignificant call is gained, and
grey means direct implanted target excluded. A missing fate for a significant
baseline call fails the plot. The output is `DEVELOPMENT_ONLY` and should be
visually reviewed before a talk.
