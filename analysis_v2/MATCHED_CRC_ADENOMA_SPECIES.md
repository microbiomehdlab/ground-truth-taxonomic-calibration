# Matched CRC/adenoma species across profilers

`plot_matched_crc_adenoma_species.py` adds a cross-profiler comparison for four directly spiked species: *Fusobacterium nucleatum*, *Parvimonas micra*, *Dialister pneumosintes*, and *Peptostreptococcus stomatis*. Unlike the shared-candidate plots, inclusion does **not** require a significant CRC association in a given profiler. Thus the MetaPhlAn *F. nucleatum* and Kraken *P. micra* and *D. pneumosintes* rows remain visible even if they are not three-cohort CRC calls. The alias file explicitly maps Kraken `Allisonella pneumosintes` to the implanted *D. pneumosintes*; the source feature name is preserved in the output.

The SVG displays separately, for each cohort and profiler, the age/sex-adjusted CRC-vs-Control and Adenoma-vs-Control effects, 95% CIs and BH q-values; full-community baseline detection and positive-only median abundance by Control/Adenoma/CRC; and direct 0.01% independent-spike recovery by condition (median, IQR, fraction with recovery below 0.5, and post-spike detection). Box colors indicate association-call status only: green = significant same-direction CRC and adenoma calls; amber = significant CRC call but adenoma not called despite same-direction effect; pink = opposite-direction adenoma effect or reversed call; blue = adenoma called without a CRC call; grey = neither called or not evaluable. Colors do not diagnose the reason for a call difference. Auditable TSVs also retain all sample-level baseline values and paired spike responses at 0.01%, 0.05%, and 0.1%.

Run from the repository root after setting `ROOT` and `SOURCE` as in `SHARED_CRC_ADENOMA_BRIDGE.md`:

```bash
OUTDIR="$ROOT/matched_crc_adenoma_species_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_matched_crc_adenoma_species.py \
  --calls "$SOURCE/models/disease/models/primary_disease_da_results.tsv" \
  --manifest "$ROOT/partial_snapshot_response/native_abundance/biomarker_profile_manifest.tsv" \
  --abundance "$ROOT/partial_snapshot_response/native_abundance/biomarker_abundance_long.tsv" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --panel spikes/spike_panel.tsv \
  --aliases examples/spike_taxon_aliases.csv \
  --outdir "$OUTDIR"
```

Open `matched_crc_adenoma_species.svg` only after `SUCCESS` exists. This is descriptive evidence, not a formal CRC-vs-Adenoma test. The two model contrasts share a Control reference; their difference and confidence interval cannot be calculated from the exported marginal CIs because the coefficient covariance is not exported. A valid direct contrast requires refitting or retaining the model covariance. Spike recovery can reveal a measurement problem but cannot prove whether an unspiked association is biological.
