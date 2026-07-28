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

## Build

Docker:

```bash
docker build \
  --tag ground-truth-taxonomic-calibration:1.0.0 \
  installation/
```

Apptainer:

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

