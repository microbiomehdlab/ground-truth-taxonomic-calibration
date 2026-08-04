# Profiling and spike-generation container

The container definition pins the software used by the upstream workflow:

| Software | Version |
|---|---:|
| Kraken2 | 2.1.6 |
| Bracken | 3.1 |
| MetaPhlAn | 4.2.2 |
| Bowtie2 | 2.5.4 |
| ART | 2016.06.05 |
| seqtk | 1.4 |
| fastp | 0.23.4 |
| FastQC | 0.12.1 |

Reference databases are kept outside the image.

## Build directly with Apptainer

The source-controlled `installation/Apptainer.def` does not require a Docker
daemon. Build under local temporary storage, install into the site-approved
image directory, and verify the installed image:

```bash
SIF=/approved/image/directory/taxonomic-tools_1.0.0.sif \
bash build_taxonomic_tools_container.sh
```

The script uses `apptainer build --fakeroot` and never executes the temporary
`/tmp` image. It records adjacent `.sha256` and `.inspect.txt` provenance
files. If a site explicitly provides an unprivileged default/remote builder
instead of fakeroot, set `BUILD_WITHOUT_FAKEROOT=true`; the script does not
automatically retry recipe failures under a different build mode.

Slurm is not required. Both images can be built in the current shell:

```bash
UPSTREAM_SIF=/approved/image/directory/taxonomic-tools_1.0.0.sif \
ANALYSIS_SIF=/approved/image/directory/crc-spike-maaslin2_1.0.0.sif \
bash build_reproducibility_containers.sh
```

On clusters where long builds should not run on login nodes, submit the same
build through Slurm:

```bash
PROJECT="$PWD" \
UPSTREAM_SIF=/approved/image/directory/taxonomic-tools_1.0.0.sif \
ANALYSIS_SIF=/approved/image/directory/crc-spike-maaslin2_1.0.0.sif \
sbatch --export=ALL build_reproducibility_containers.sbatch
```

This writes a build manifest under `work/container_build/` containing image,
definition, complete recipe-input, Git-commit, and Apptainer-version
provenance. Use new versioned image names instead of overwriting prior images.

## Alternative OCI/Docker build

Docker:

```bash
docker build \
  --tag ground-truth-taxonomic-calibration:1.0.0 \
  installation/
```

Convert a locally built Docker image to Apptainer:

```bash
apptainer build taxonomic-tools_1.0.0.sif \
  docker-daemon://ground-truth-taxonomic-calibration:1.0.0
```

If Docker is unavailable on the cluster, build and publish the OCI image on a
trusted build host, then use:

```bash
apptainer pull taxonomic-tools_1.0.0.sif \
  docker://REGISTRY/ground-truth-taxonomic-calibration:1.0.0
```

Record the immutable OCI digest and:

```bash
sha256sum taxonomic-tools_1.0.0.sif
```

## Verify

```bash
apptainer exec taxonomic-tools_1.0.0.sif \
  micromamba run -n taxonomic_tools bash -c '
    set -euo pipefail
    art_illumina --version 2>&1 | head -n 2
    seqtk 2>&1 | head -n 2
    fastp --version
    kraken2 --version
    bracken -v
    metaphlan --version
    bowtie2 --version
  '
```

## Kraken2/Bracken database

The paper used UHGG v2.0.2:

```text
https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/
```

The database must include the Bracken distribution matching the fixed
100-nucleotide read-length model. Confirm:

```bash
test -s "$K2_DB/hash.k2d"
test -s "$K2_DB/opts.k2d"
find "$K2_DB" -maxdepth 1 -name '*100*mers.kmer_distrib' -print
```

Record the source URL, download date, total size, and a sorted checksum
inventory of all database files.

## MetaPhlAn database

Install the exact release used for the publication into a persistent directory:

```bash
apptainer exec --bind "$MPA_DB:$MPA_DB" taxonomic-tools_1.0.0.sif \
  micromamba run -n taxonomic_tools metaphlan \
    --install \
    --db_dir "$MPA_DB"
```

Record the database identifier printed in MetaPhlAn outputs and checksums for
the Bowtie2 index and taxonomy files. Do not rely only on “latest,” because the
default database can change over time.

## Configuration

Copy `config/global.env.example` to `config/global.env` and set the image and
database paths. Run `scripts/preflight.sh` before submitting jobs.
