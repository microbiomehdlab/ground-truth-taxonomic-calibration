# Reproducing the taxonomic-calibration study

This is the authoritative navigation document for a clean reproduction. Large
raw reads, reference databases, Apptainer images, scheduler state, and generated
results are intentionally not stored in Git. Their releases, paths, checksums,
and connections to outputs must be captured in the run evidence.

## Scientific scope

The experiment implants known paired reads into already sequenced metagenomes.
Its controlled ground truth is therefore read identity, implanted-pair count,
and final-library read fraction—not cellular abundance, biomass, extraction
efficiency, or an identical abundance estimand across profilers.

Profiler-native MetaPhlAn and Kraken2/Bracken species outputs are the primary
practical measurements. Denominator-harmonized representations are sensitivity
analyses. The primary quantitative endpoint for the final analysis is a paired
baseline-adjusted dose response, with biological sample treated as the unit of
replication.

## Reproduction map

| Stage | Authoritative entry point | Required evidence |
|---|---|---|
| Build software | `build_reproducibility_containers.sh` | image hashes, definitions, build manifest, tool versions |
| Configure host | `config/*.env.example`, `datasets/crc.env.example`, `datasets/yachida/yachida.env.example` | ignored site-specific env plus recorded non-secret values |
| Obtain target genomes | `spikes/scripts/spikein/download_refs_ncbi.sh` | versioned accessions, NCBI reports, FASTA checksums |
| Build and freeze pools | `spikes/README.md` | settings, seeds, pair counts, checksums, preproduction audit |
| Freeze cohorts | `datasets/CRC_MANIFESTS.md`, `datasets/yachida/README.md` | eligibility, run accessions, URLs, sizes, MD5s, selection provenance |
| Strict preprocessing and profiling | `datasets/submit_crc_batch.sh`, `datasets/yachida/submit_batch.sh` | pinned MetaShotgunPrep commit, host/database/image identity, receipts and output checksums |
| Audit batches | `run_yachida_batch_audit.sbatch` and cohort completeness checks | verified sample markers, sealed batches, failure/retry ledger |
| Pure-pool assignment | `run_yachida_pure_pool_audit.sbatch` | systematic full-pool sample, native profiles, checksummed summary |
| Target/reference validation | `run_phase3_target_validation.sh` | assembly-integrity and NCBI quality tables, UHGG/MetaPhlAn table, provenance, checksum seal |
| Replacement-assembly sensitivity | `scripts/rank_candidate_assemblies.py` | frozen NCBI snapshots, eligibility ledger, selected accessions, checksums |
| Stage downstream inputs | `stage_required_inputs.sh`, `REQUIRED_INPUTS.md` | self-contained input contract and checksums |
| Statistical analysis | `run_original_unpaired_q010_cluster.sbatch` | frozen analysis specification, run manifest, model diagnostics and outputs |
| Rebuild figures | `rerun_current_figures_only.sh`, `rerun_supplementary_figures_only.sh` | source tables, vector figures, checksums |

The statistical entry-point names retain `original_unpaired` for compatibility
with the historical analysis. They must not be interpreted as approval of that
model for the final manuscript: the final paired/repeated-measures analysis
specification must be frozen before the definitive comparison is run.

## 1. Clone and validate the source

```bash
git clone git@github.com:microbiomehdlab/ground-truth-taxonomic-calibration.git
cd ground-truth-taxonomic-calibration
git rev-parse HEAD
bash scripts/check_repository.sh
bash tests/test_upstream_tools.sh
```

Record the full commit, not only its short form. Do not run publication jobs
from a checkout with unexplained tracked modifications.

## 2. Build or identify the pinned images

Follow [`installation/README.md`](installation/README.md). For a final run,
record both Apptainer image paths and SHA-256 hashes. A pre-existing image is
acceptable only when its hash and verification output are retained.

## 3. Configure without committing cluster paths

Copy the relevant examples, fill in absolute database, image, scratch, state,
result, and QC paths, and keep these `.env` files ignored. Before submission,
print and archive the non-secret configuration together with its SHA-256.

Never infer the software revision from a mutable directory name. Configuration
must include the expected MetaShotgunPrep commit, and the runner must verify the
actual checkout against it.

## 4. Freeze inputs before outcome inspection

- Use `spikes/spike_panel.tsv` for the versioned target assemblies.
- Use the committed Feng and Zeller production manifests.
- Use the audited Yachida 67/67/67 design and nested 10/10/10 subset.
- Never replace a biological sample because of profiler or recovery results.
- Treat transfer/MD5 failures as retryable operational failures, not exclusions.

## 5. Require gates, not merely successful Slurm states

A Slurm `COMPLETED` state is necessary but insufficient. A final sample must
have verified input downloads, synchronized mates, construction receipts,
both native profiler outputs, retained-output checksums, and the applicable
verified marker. Final datasets require a completeness report and sealed batch
markers. Use `test -f` for markers created with `touch`.

## 6. Preserve a publication evidence package

For every final claim, retain:

1. scientific decision and rationale;
2. exact repository and MetaShotgunPrep commits;
3. image and database identities/checksums;
4. frozen input and output manifests;
5. automated validation output;
6. failure/retry and protocol-deviation records;
7. manuscript-ready method wording and stated limitations.

Private evidence may remain under ignored `work/`, but it must be backed up and
must not exist only in a terminal scrollback. Public code, example configuration,
schemas, and non-sensitive methods belong in Git.

## Current qualification

The repository contains reproducible upstream production and audit components,
but a fully frozen paper reproduction is not yet claimed until all three final
cohorts are complete, all batches are audited, the remaining sensitivity checks
are resolved, and the paired statistical-analysis specification and downstream
figures are finalized.
