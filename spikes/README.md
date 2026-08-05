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
  --outdir references/genomes \
  --image /path/to/taxonomic-tools.sif
```

Use a new, empty `references/genomes` directory for a from-scratch audit. The
downloader requests the exact versioned assembly accessions, retains NCBI data
reports and the pinned CLI version under `references/genomes/provenance/`, and
writes both byte-level and header/line-wrapping-independent sequence checksums
to `references/genomes/reference_genome_checksums.tsv`.

Compare freshly downloaded FASTAs against an archived study-genome directory:

```bash
python3 scripts/audit_reference_genomes.py \
  --panel spikes/spike_panel.tsv \
  --fasta-root /path/to/archived/study/genomes \
  --output work/reference_audit/original_reference_genome_checksums.tsv \
  --expected references/genomes/reference_genome_checksums.tsv \
  --comparison-output work/reference_audit/fresh_vs_original.tsv
```

The comparison fails on any contig-sequence, contig-count, or base-count
difference. Byte-level SHA differences alone are reported but accepted because
FASTA headers, wrapping, or contig order may change without changing biological
sequence content.

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
Rerunning selected fractions or reordered sample manifests preserves the exact
seed for every sample/fraction/member combination.

## Deterministic seed contract

`stable_seed.py` implements the frozen `stable-seed-v1` SHA-256 derivation.
ART pool simulation uses the label, assembly, and simulation parameters.
Independent subsampling uses sample, target, and fraction; community
subsampling additionally uses the community and member label. Scheduler task
IDs, manifest row positions, fraction positions, batch sizes, and filesystem
paths are excluded. Exact seeds are recorded in the design/provenance tables.

Changing `SEED_BASE` intentionally creates a different realization. Changing
pool parameters requires a new `POOLS_DIR`; do not reuse pools created before
the seeded ART workflow when claiming byte-for-byte reproduction.

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
