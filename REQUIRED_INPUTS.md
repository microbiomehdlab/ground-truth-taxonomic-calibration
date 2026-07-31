# Required inputs for a self-contained project folder

The cluster project must contain the following **seven top-level input files**:

```text
samples_spiked_all_independent.tsv
samples_spiked_all_community.tsv
spike_panel.tsv
spike_taxon_aliases.csv
metadata_w_study.tsv
metaphlan4_merged_unspecified.csv
kraken2_bracken_merged_unspecified.csv
```

It must also contain this complete generated-data directory:

```text
profile_results/
```

The expected profiler-result layout is:

```text
profile_results/
├── Bfrag/
├── Csym/
├── Dpne/
├── Fnuc/
├── Hhat/
├── Pana/
├── Pint/
├── Pmic/
├── Porp/
├── Psto/
└── CRCpanel/
```

Each spike-label directory contains profiler subdirectories such as
`kraken2_bracken/` and `metaphlan4/`, holding the spike-fraction abundance
tables.

## Upstream-to-analysis output contract

The independent-spike manifest must contain exactly these six final fractions:

```text
0.0001  0.0005  0.001  0.005  0.01  0.05
```

The community-spike manifest must contain exactly these seven total fractions:

```text
0.0001  0.0005  0.001  0.005  0.01  0.05  0.10
```

For every spike label, profiler, and expected fraction, the upstream workflow
must produce:

```text
profile_results/<label>/<tool>/spike_<fraction>.open_world.csv
```

Required spike labels:

```text
Bfrag Csym Dpne Fnuc Hhat Pmic Pana Psto Porp Pint CRCpanel
```

Required profiler directories:

```text
metaphlan4
kraken2_bracken
```

Each abundance CSV must contain a unique `sample` column, numeric profiler
features, and zero for absent features. Profiler-native labels must remain
unchanged and the tables must not be renormalized after feature selection.
The unspiked profiler tables must cover every biological background sample,
with sample IDs mapping unambiguously to those in the spike manifests.

The upstream run should retain its spike-design logs, including inserted read
pairs, achieved fractions, and random seeds. The downstream preflight validates
the complete 10-target × 2-profiler alias table against the unspiked abundance
table headers.

## What these inputs represent

| Path | Purpose |
|---|---|
| `samples_spiked_all_independent.tsv` | Independent spike-in sample manifest |
| `samples_spiked_all_community.tsv` | Community spike-in sample manifest |
| `spike_panel.tsv` | Canonical spike design and target membership |
| `spike_taxon_aliases.csv` | Complete 10-target × 2-profiler harmonization table (20 rows, including identity mappings) |
| `metadata_w_study.tsv` | Sample study and diagnostic-background metadata |
| `metaphlan4_merged_unspecified.csv` | Unspiked MetaPhlAn 4 abundance table |
| `kraken2_bracken_merged_unspecified.csv` | Unspiked Kraken2/Bracken abundance table |
| `profile_results/` | Profiled abundance tables for all spike-in experiments |

The repository already includes the shared `R/` helper code; it is not an
external input.

No prior `RUNS*` output directory is required for a clean run. The following
are generated outputs and should not be copied as inputs:

```text
design_auto/
spike_metrics/
maaslin_spike/
species_driver_and_thresholds/
manuscript_figures/
RUNS*/
```

## Create the self-contained local folder

From the repository root, point `SOURCE_DATA_ROOT` to the completed upstream
output directory:

```bash
SOURCE_DATA_ROOT=".." bash stage_required_inputs.sh
```

The script copies only the required inputs into the current folder, preserves
the existing reviewed code, and then prints checksums and sizes.

After staging, make the wrappers executable (file modes may be lost during a
manual copy) and run the structural check with the analysis container:

```bash
chmod +x ./*.sh

SIF=/path/to/crc_spike_original_unpaired_q010_r433_maaslin2_1180_mash23_v2.sif \
PROJECT="$PWD" \
bash cluster_preflight_original_unpaired.sh
```

Once this passes, the folder no longer depends on the old project directory and
can be transferred with:

```bash
LOCAL_PROJECT="$PWD" \
DATA_ROOT="$PWD" \
REMOTE=user@cluster \
REMOTE_PROJECT=/path/to/project \
bash transfer_analysis_to_cluster.sh
```
