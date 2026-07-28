# Spike-in workflow

This directory implements the two designs used in the paper:

- independent insertion of each of ten CRC-associated taxa at final read-pair
  fractions `0.0001,0.0005,0.001,0.005,0.01,0.05`;
- equal-weight insertion of all ten taxa at total final fractions
  `0.0001,0.0005,0.001,0.005,0.01,0.05,0.10`.

For a background containing `R` synchronized read pairs and target final
fraction `f`, the workflow inserts:

```text
N = round(f / (1 - f) * R)
```

so the achieved fraction is `N / (R + N)`. Community `N` is divided among the
ten taxa using the weights in `spike_panel.tsv`; the publication design uses
equal weights.

## Configuration

```bash
cp spikes/spikein.env.example spikes/spikein.env
```

Edit the local file. The independent and community fraction lists are separate
to prevent accidental use of the 10% community-only fraction in the
independent design.

## Reference genomes

`spike_panel.tsv` records the canonical taxon, assembly accession, local FASTA
path, and community weight. Download missing assemblies:

```bash
bash spikes/scripts/spikein/download_refs_ncbi.sh \
  --panel spikes/spike_panel.tsv \
  --outdir references/genomes
```

## Pools

```bash
bash spikes/scripts/spikein/spikein_prepare_pools.sh \
  --env spikes/spikein.env
```

This submits ART simulation followed by fastp preprocessing for every taxon.
Use a new `POOLS_DIR` whenever pool-generation parameters change. Existing
paired pools are reused only after consistency checks.

## Independent design

```bash
bash spikes/scripts/spikein/spikein_run_independent.sh \
  --env spikes/spikein.env
```

One Slurm array is submitted per target taxon. Every task handles one
background sample and all independent fractions.

## Community design

```bash
bash spikes/scripts/spikein/spikein_run_community.sh \
  --env spikes/spikein.env \
  --community-label CRCpanel
```

Every task handles one background sample and all total-community fractions.
When rerunning only selected fractions, seed indices remain anchored to the
full fraction list.

## Outputs and provenance

Each array task records:

- background read-pair counts;
- requested and inserted spike-pair counts;
- nominal and achieved fractions;
- deterministic seeds;
- pool and output paths; and
- a manifest fragment used to build the combined profiler manifest.

Temporary files are deleted after successful or failed tasks unless
`KEEP_TMP=1`.

## Profile-table organization

After profiling, create the primary open-world matrices:

```bash
python3 spikes/scripts/spikein/organize_profile_tables_by_spike.py \
  --inroot results \
  --outroot profile_results \
  --kind open_world
```

Then merge across labels and fractions if needed:

```bash
python3 spikes/scripts/spikein/merge_profile_tables.py \
  --root profile_results \
  --outdir merged_profiles \
  --kind open_world
```

Open-world values are never renormalized. Closed-world tables remain available
as an optional sensitivity output through `--kind closed_world` or
`--kind both`.

