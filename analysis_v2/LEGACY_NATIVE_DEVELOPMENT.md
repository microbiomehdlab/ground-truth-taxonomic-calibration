# Historical native profiles for development

Historical Feng and Zeller native profiler outputs may be used to develop and
exercise the v2 workflow while strict-production cohorts are incomplete. They
are **not manuscript evidence** and must never be silently mixed with final
strict-production results.

`scripts/build_legacy_native_development_input.py` reads native Bracken species
tables and MetaPhlAn tables from `results_old_v2`, maps ERR run accessions to
the frozen cohort manifests, keeps the frozen dose grid, and writes a validated
canonical input plus an exclusion ledger and checksums. By default it retains
only baseline–spike profiles available for both profilers, so direct profiler
comparisons use the same biological samples and perturbations.
Both the historical `kraken_bracken` directory name and the current
`kraken2_bracken` name are recognized and mapped to the canonical
`kraken2_bracken` profiler identifier.

Some frozen biological samples combine multiple ERR runs. Historical outputs
profiled those runs separately, so they cannot be reconstructed as the final
concatenated sample. For development only, the adapter deterministically keeps
the lexicographically first ERR accession per biological sample and records
the remaining run directories as `non_primary_run_for_multirun_sample`. This
prevents pseudoreplication; it is not a substitute for the definitive combined
sample.

The historical run did not preserve exact implanted pair-allocation records.
Consequently, its directory-name fractions are treated as nominal doses and
the generated output is stamped `DEVELOPMENT_ONLY`. Community (`CRCpanel`)
profiles are provisionally interpreted as an equal ten-member mixture. Neither
assumption is permitted in the definitive analysis, which must use sealed
strict-production design tables and exact target fractions.

Example on Lobo:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration

DEV_OUT="$PWD/work/analysis_v2_legacy_native_dev_$(date +%Y%m%d_%H%M%S)"

python3 analysis_v2/scripts/build_legacy_native_development_input.py \
  --results-root /mnt/nfs/microbiomehd/tax_benchmarking_for_biomarker/results_old_v2 \
  --feng-manifest datasets/fengq/manifests/production_manifest.tsv \
  --zeller-manifest datasets/zellerg/manifests/production_manifest.tsv \
  --spike-panel spikes/spike_panel.tsv \
  --aliases examples/spike_taxon_aliases.csv \
  --profiler-coverage paired \
  --outdir "$DEV_OUT"

test -s "$DEV_OUT/DEVELOPMENT_ONLY.txt" &&
test -s "$DEV_OUT/validation/SUCCESS" &&
test -s "$DEV_OUT/development_input.sha256" &&
echo "[PASS] Historical development input is ready"
```

Use `--profiler-coverage available` only for profiler-specific engineering
checks. It is not appropriate for direct profiler comparisons because sample
availability may differ. All reports derived from this adapter must inherit
the development-only status. When strict outputs are complete, replace this
adapter output with the canonical input produced by the definitive cohort
builder; downstream model and reporting interfaces remain unchanged.
