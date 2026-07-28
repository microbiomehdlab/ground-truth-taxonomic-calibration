# Ground-truth taxonomic calibration workflow

This repository contains the complete computational workflow for the
accompanying CRC-biomarker benchmarking study: public-read acquisition,
controlled spike construction, taxonomic profiling, quantitative recovery,
MaAsLin2 biomarker analysis, artefact filtering, calibratability analysis, and
manuscript Figures 2–7. It supports the two profilers used in the paper:

- MetaPhlAn 4
- Kraken2 with Bracken

The primary quantitative analyses use open-world species abundances so that
observed abundances and designed spike fractions retain the same all-read
denominator.

## Workflow

```text
public paired-end FASTQs
        |
        v
reference genomes -> simulated read pools
        |                    |
        +---- independent and 10-member community spike-ins
                             |
                             v
                   MetaPhlAn 4 / Kraken2+Bracken
                             |
                             v
                sample-level open-world profiles
                             |
                             v
             fraction-level matrices in profile_results/
                             |
                             v
                  downstream statistical analysis
```

The repository does not redistribute raw reads, profiler databases, or
generated results. Users obtain or generate those external inputs as described
below.

## Repository layout

```text
config/                       host and dataset configuration templates
datasets/                     public-study metadata and download helpers
installation/                 upstream profiling-container definition
spikes/
  spike_panel.tsv             ten canonical spike targets
  spikein.env.example         spike-generation configuration template
  scripts/spikein/            pool, spike-in, organization, and merge scripts
taxonomy/                     Slurm profiling submission drivers
workflows/
  metaphlan4/                 MetaPhlAn profiling and post-processing
  kraken2_bracken/            Kraken2/Bracken profiling and post-processing
R/                            reusable downstream R functions
scripts/                      analysis/plot entry points and repository checks
containers/                   downstream MaAsLin2 analysis environment
examples/                     small manifest and taxon-alias examples
tests/                        dependency-free workflow smoke tests
*.sh and *.sbatch             top-level staging, validation, and run entry points
```

Generated FASTQs, databases, containers, results, logs, and local `.env` files
are intentionally excluded from Git.

### Main entry points

| Command | Purpose |
|---|---|
| `scripts/preflight.sh` | Validate the configured upstream profiling environment |
| `taxonomy/run_profiling.sh` | Submit taxonomic profiling for a sample manifest |
| `stage_required_inputs.sh` | Stage completed profiler outputs for analysis |
| `preflight_all.sh` | Validate the integrated upstream and downstream workflow |
| `run_original_unpaired_q010_cluster.sbatch` | Submit the complete statistical analysis and figure regeneration |
| `run_publication_original_unpaired_q010.sh` | Run the same downstream workflow inside an allocated job |
| `transfer_analysis_to_cluster.sh` | Transfer the reviewed analysis code and staged inputs |

`R/` and `scripts/` intentionally have different roles: files under `R/`
define shared functions, while files under `scripts/` are executable analysis
stages that source those functions. Keeping this conventional separation makes
dependencies explicit and permits individual stages to be tested or rerun.

## Requirements

- Linux and Slurm
- Apptainer or Singularity
- Python 3
- `rsync`, `awk`, `gzip`, and standard Unix command-line tools
- an image containing ART, fastp, MetaPhlAn 4, Kraken2, and Bracken
- a MetaPhlAn 4 database
- a Kraken2/Bracken database built for the configured read length
- R 4.3.3 and MaAsLin2 1.18.0 through the included analysis-container
  definition

The exact database release and image digest used for a publication run should
be recorded with the release metadata. Database contents are too large to
store in Git.

## 1. Configure the host

```bash
cp config/global.env.example config/global.env
cp spikes/spikein.env.example spikes/spikein.env
```

Edit both copied files. They are ignored by Git because they contain
machine-specific paths.

The study spike panel contains assembly accessions and the ten canonical target
labels. Its `fasta` column must point to locally available reference genomes:

```text
Bfrag, Csym, Dpne, Fnuc, Hhat, Pmic, Pana, Psto, Porp, Pint
```

## 2. Obtain public reads

See [datasets/README.md](datasets/README.md). Verify the applicable study data
access conditions before downloading or redistributing metadata.

Create a tab-separated sample manifest with:

```text
sample_id    fastq1    fastq2
```

Paths may be absolute or valid from the directory where jobs run.

## 3. Run preflight

```bash
GLOBAL_ENV="$PWD/config/global.env" \
SPIKE_ENV="$PWD/spikes/spikein.env" \
bash scripts/preflight.sh
```

Do not submit jobs unless this ends with:

```text
[PASS] Upstream workflow preflight passed.
```

Run the dependency-free utility smoke tests:

```bash
bash tests/test_upstream_tools.sh
bash scripts/check_repository.sh
```

## 4. Generate deterministic spike pools

Download the exact reference assemblies declared in the panel:

```bash
bash spikes/scripts/spikein/download_refs_ncbi.sh \
  --panel spikes/spike_panel.tsv \
  --outdir references/genomes
```

Record checksums after download:

```bash
find references/genomes -type f -name '*.fa' -print0 |
  sort -z |
  xargs -0 sha256sum > reference_genomes.sha256
```

Then generate the pools:

```bash
bash spikes/scripts/spikein/spikein_prepare_pools.sh \
  --env spikes/spikein.env
```

ART generates paired reads from each reference genome. fastp then creates the
final pools. Pool parameters and submitted job IDs are recorded under
`POOLS_DIR`.

## 5. Generate spike-in datasets

Independent spikes:

```bash
bash spikes/scripts/spikein/spikein_run_independent.sh \
  --env spikes/spikein.env
```

Ten-member community spikes:

```bash
bash spikes/scripts/spikein/spikein_run_community.sh \
  --env spikes/spikein.env \
  --community-label CRCpanel
```

The configured random seed and full fraction list determine all subsampling.
Subset reruns retain the seed associated with each original fraction.

The publication design uses six independent final fractions
(`0.01%`–`5%`) and seven total-community fractions (`0.01%`–`10%`).
These are configured separately as `INDEPENDENT_FRACTIONS` and
`COMMUNITY_FRACTIONS`; the effective per-taxon community fraction is one tenth
of the configured total fraction.

## 6. Profile samples

Prepare `config/dataset.env` from the example, or submit the array driver with a
three-column sample manifest:

```bash
N=$(awk 'NR>1 && NF && $1 !~ /^#/ {n++} END {print n+0}' samples.tsv)

SAMPLES_TSV="$PWD/samples.tsv" \
REPO_ROOT="$PWD" \
sbatch --array="1-${N}%5" taxonomy/run_profiling_array.sbatch
```

Each array task submits the enabled MetaPhlAn and Kraken2/Bracken jobs, waits
for them, and post-processes only jobs Slurm reports as completed. Any failed
profiling or post-processing step produces a non-zero exit status.

## 7. Organize profiler outputs

The primary downstream input is one open-world abundance matrix per spike
target, profiler, and fraction:

```bash
python3 spikes/scripts/spikein/organize_profile_tables_by_spike.py \
  --inroot results \
  --outroot profile_results \
  --kind open_world
```

Merge sample/fraction matrices when a combined profiler table is required:

```bash
python3 spikes/scripts/spikein/merge_profile_tables.py \
  --root profile_results \
  --outdir merged_profiles \
  --kind open_world
```

The output contract is:

```text
profile_results/
  Bfrag/
    metaphlan4/spike_f0p0001.open_world.csv
    kraken2_bracken/spike_f0p0001.open_world.csv
  ...
  CRCpanel/
```

## 8. Stage the downstream analysis inputs

The statistical workflow consumes:

- the independent and community sample manifests;
- `spike_panel.tsv`;
- the complete profiler-specific taxon-alias table;
- unspiked MetaPhlAn and Kraken2/Bracken abundance tables;
- study/condition metadata; and
- the complete `profile_results/` directory above.

The upstream workflow generates the spike manifests and `profile_results`.
The two original profiler tables and combined study metadata are generated from
the corresponding unspiked FengQ/ZellerG runs. The complete alias mapping is
maintained with the analysis inputs; the publication-design template is
provided at `examples/spike_taxon_aliases.csv`.

Stage the seven tables and `profile_results/` into the repository root:

```bash
SOURCE_DATA_ROOT=/path/to/completed/upstream_outputs \
DESTINATION="$PWD" \
bash stage_required_inputs.sh
```

The exact requirements and the upstream-to-analysis output contract are
documented in [REQUIRED_INPUTS.md](REQUIRED_INPUTS.md).

## 9. Build or verify the analysis container

The downstream benchmark is pinned to R 4.3.3 and MaAsLin2 1.18.0:

```bash
SIF=/path/to/crc_spike_original_unpaired_q010_r433_maaslin2_1180_v1.sif \
BUILD_TMPDIR=/tmp \
bash build_crc_spike_maaslin2_container.sh
```

The helper builds under local temporary storage, installs the final image at
`SIF`, and verifies the R environment.

To transfer a staged analysis to a cluster:

```bash
LOCAL_PROJECT="$PWD" \
DATA_ROOT="$PWD" \
REMOTE=user@cluster \
REMOTE_PROJECT=/path/to/project \
bash transfer_analysis_to_cluster.sh
```

## 10. Validate the complete repository

```bash
GLOBAL_ENV="$PWD/config/global.env" \
SPIKE_ENV="$PWD/spikes/spikein.env" \
ANALYSIS_SIF=/path/to/crc_spike_original_unpaired_q010_r433_maaslin2_1180_v1.sif \
bash preflight_all.sh
```

This checks the upstream container and databases, spike panel and FASTQs,
analysis inputs, complete taxon-alias mapping, R helpers, R packages, R
version, and MaAsLin2 version.

## 11. Run all statistical analyses and regenerate Figures 2–7

Submit one clean sequential Slurm job:

```bash
RUN_ROOT="RUNS_publication_original_unpaired_q010_$(date +%Y%m%d_%H%M%S)"

PROJECT="$PWD" \
SIF=/path/to/crc_spike_original_unpaired_q010_r433_maaslin2_1180_v1.sif \
RUN_ROOT="$RUN_ROOT" \
sbatch run_original_unpaired_q010_cluster.sbatch
```

The authoritative runner:

1. validates the complete 10-target × 2-profiler alias mapping;
2. builds the spike design;
3. computes independent/community recovery metrics;
4. runs the original unpaired MaAsLin2 benchmark at `q ≤ 0.10`;
5. computes biomarker-recoverability drivers;
6. regenerates manuscript Figures 2–5;
7. regenerates the artefact-exclusion and calibratability figures; and
8. validates outputs and records checksums/provenance.

Outputs are written under:

```text
RUN_ROOT/
  design_auto/
  spike_metrics/
  maaslin_spike/
  species_driver_and_thresholds/
  manuscript_figures/
  INPUT_AND_CODE_SHA256.txt
  EXPECTED_OUTPUTS.txt
  PROVENANCE.txt
  logs/pipeline.log
```

## Reproducibility records

For every released run, retain:

- the Git commit or release tag;
- the container SHA-256 digest;
- MetaPhlAn database identifier and checksum;
- Kraken2/Bracken database source, build command, date, and checksum inventory;
- reference-genome assembly accessions and FASTA checksums;
- input-manifest checksums;
- Slurm logs and job IDs;
- generated parameter/provenance files; and
- checksums for the final abundance matrices.

## License and citation

Code is released under the [MIT License](LICENSE). Add the final repository URL
and paper DOI to [CITATION.cff](CITATION.cff) when they are assigned.
