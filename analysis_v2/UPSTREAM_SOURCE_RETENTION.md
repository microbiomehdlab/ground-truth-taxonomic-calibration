# Upstream source files for the public release

**Status:** provisional source-retention decision, 29 September 2026. No files are deleted by this plan. The verified `upstream_evidence_v1` package has its own exact 56-file `MANIFEST.tsv`; that manifest governs data, tables, figures, and provenance to deposit. This document governs which source files should remain understandable in the public Git tree.

## Keep for upstream reproduction

| Role | Source paths |
|---|---|
| Entry and explanation | `README.md`, `REPRODUCING.md`, `analysis_v2/UPSTREAM_METHODS_FOR_MANUSCRIPT.md`, `analysis_v2/UNIFIED_UPSTREAM_SEAL.md`, `analysis_v2/UPSTREAM_EVIDENCE_PACKAGE.md`, `analysis_v2/TAXON_IDENTITY_FREEZE.md` |
| Frozen target identity | `spikes/spike_panel.tsv`, `examples/spike_taxon_aliases.csv`, `analysis_v2/feature_equivalence_aliases.tsv`, `analysis_v2/taxon_identity_freeze.sha256` |
| Cohort selection and manifests | `datasets/yachida/` production selection and audit code, `datasets/fengq/manifests/`, `datasets/zellerg/manifests/`, `datasets/build_crc_production_manifest.py`, `scripts/select_samples_deterministically.py`, and their tests |
| Preprocessing and production | `datasets/yachida/run_metashotgunprep.sh`, `datasets/yachida/extract_strict_unmapped_pairs.py`, `datasets/yachida/run_full_streaming_sample.sh`, `datasets/yachida/stream_sample.py`, `datasets/yachida/submit_batch.sh`, `datasets/submit_crc_batch.sh`, `run_crc_batch_sample.sbatch`, `workflows/kraken2_bracken/classify_bracken.sbatch`, `workflows/metaphlan4/profile.sbatch`, and their direct dependencies |
| Spike construction | `spikes/scripts/spikein/` generation, seed, allocation, pool, and read-sampling code plus `spikes/README.md`; retain validation tests and frozen reference metadata |
| Audit and package | `analysis_v2/scripts/seal_cohort_upstream.py`, `analysis_v2/run_cohort_upstream_audit.sbatch`, `analysis_v2/scripts/build_upstream_evidence_package.py`, `analysis_v2/scripts/plot_upstream_evidence.R`, `analysis_v2/run_upstream_evidence_package.sh`, `analysis_v2/templates/upstream_evidence/`, and their tests |
| Rebuild environment | relevant `containers/` recipes and validation code, installation instructions, public configuration examples, and database/reference retrieval instructions |
| Reported sensitivity | Yachida assembly-sensitivity configuration, sample runner, audit code, and tests used for the 360-profile comparison |

The definitive release inventory must expand these groups to individual tracked files and verify the dependency graph after downstream runs and figures are final. Files named `legacy`, `smoke`, `pilot`, or `postprocess_local` are review candidates, not automatically removable: some establish selection or validation provenance. Keep the historical preprint code clearly labelled if it is cited as a comparator. Presentation-only, abandoned development, and superseded duplicate scripts can be removed from the release branch after their outputs are confirmed absent from the final paper. Git history retains removed files.

## Share from the upstream evidence package

Deposit the package's release-safe `manuscript_supplement/`, `zenodo/`, `source_seal_projections/`, and `provenance/` content with its top-level `README.md`, `MANIFEST.tsv`, `SHA256SUMS`, and `SUCCESS`, subject to the final manual figure/source-table review. Its `MANIFEST.tsv` is the authoritative per-file release list and checksum record. Do not put generated package data in Git. The default package excludes raw or cleaned reads, databases, container images, scratch, scheduler logs, cluster-absolute paths, and individual age/sex/BMI values. Record external asset versions and checksums so they can be obtained independently.

Before tagging, run the public path/secret scan, validate a clean clone using the release instructions and archived inputs, and have the final retained-file list reviewed against every final manuscript panel and table.
