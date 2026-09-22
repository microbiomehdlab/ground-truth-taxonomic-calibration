# Provisional CRC replication–reliability bridge

This analysis asks whether **local spike-measured measurement reliability**
tracks the external replication of an **unspiked CRC biomarker call**. It uses
the already-written `reliability_certificates.tsv` and existing disease-model
calls; it does **not** repeat the 27-million-row response calculation.

Run on lobo in a suitable compute allocation, using a new output directory:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
CALLS="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/disease/models/primary_disease_da_results.tsv"
CERTS="$ROOT/partial_snapshot_response/response_analysis/reliability_certificates.tsv"
OUTDIR="$ROOT/crc_replication_reliability_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/bridge_crc_replication_reliability.py \
  --calls "$CALLS" --certificates "$CERTS" \
  --aliases examples/spike_taxon_aliases.csv --outdir "$OUTDIR"
```

Open `replication_by_reliability_tertile.png`. The underlying
`candidate_destination_bridge.tsv` retains every source–destination pair,
including missing destination fits and missing spike scores. The figure shows
replication fractions in **within-profiler, within-destination-cohort**
reliability tertiles, with raw numerators/denominators. The score comes from
the certificate's *evaluation* field, measured in the destination cohort;
its training score is kept separately in the table. Each source discovery is
an original-assembly, unspiked community CRC-vs-Control call with BH q ≤ 0.05
from the primary age/sex model. Replication requires a destination fit with
the same effect direction and BH q ≤ 0.05. Opposite-direction results and
non-evaluable fits are distinct.

This is **descriptive**, not a prospective prediction or a causal explanation.
One source feature can contribute twice (one for each destination), and one
destination fit can be counted from multiple source discoveries. No p-value
is attached to the tertile display. Before claiming an association, inspect
the raw table for baseline prevalence, destination sample size, effect size,
cohort imbalance, and the fraction without a local spike score. Directly
implanted taxa are not automatically excluded: each certificate is built from
*other* implanted targets, so its score describes bystander measurement
reliability, not that taxon's direct recovery.

The pre-existing optional `biomarker_replication.tsv` in the response analyzer
uses a different source/holdout orientation and should **not** be interpreted
as this bridge. The analyzer now accepts the spike-panel CSV alias schema, but
its absence does not prevent this standalone analysis. Do not rerun the large
response calculation solely to obtain this plot.
