# Public release and reproducibility checklist

**Status:** active, not yet release-complete
**Last updated:** 23 September 2026

This document separates work that can be documented now from information that
must be filled only after the definitive three-cohort run. It is the concise
release-facing checklist. Detailed scientific decisions remain in
`analysis_v2/METHODS_DECISION_LOG.md`; live operational detail remains in
`NEXT_ACTIONS.md` and `WORK_HANDOFF.md`.

## Current reproducibility claim

The repository contains source-controlled upstream production, validation,
paired-analysis, and reporting workflows. It does not yet constitute the final
frozen reproduction package for the revised paper. Current multi-cohort
figures are development evidence unless a definitive seal explicitly says
otherwise.

The controlled truth is implanted paired-read identity, pair count, and final-
library read fraction. It is not cellular abundance, biomass, extraction
efficiency, or biological presence/absence.

## Documentation completed now

- [x] Identify `analysis_v2` as the current paired publication workflow.
- [x] Retain the original unpaired workflow as a labeled historical comparator.
- [x] Provide installation, required-input, cohort-manifest, spike-pool, and
  definitive-run documentation.
- [x] Record the distinction between development and definitive evidence.
- [x] Record the 23 September Feng checkpoint without presenting it as final.
- [x] Keep machine-specific paths and large generated assets outside Git.
- [x] Maintain a license and `CITATION.cff`.
- [x] Separate the GUTBIOME pilot from the paper's spike benchmark.

## Required before definitive analysis

- [x] Finish Feng and verify 154/154 `.verified` markers and retained-output
  receipts.
- [x] Repair and regression-test the CRC seal's zero-byte sentinel handling,
  retained-output hashing, persistent-root constraint, and stale-seal
  invalidation locally.
- [x] Harden the CRC seal to enforce baseline, independent, and community
  profile topology separately and validate retained completion tables.
- [x] Transfer the topology-hardened seal to the cluster and reseal Feng and
  Zeller without rerunning profiling.
- [x] Complete and checksum-verify the first Feng and Zeller seals against all
  retained receipt contents.
- [x] Revalidate the existing Yachida production and assembly-sensitivity seals.
- [ ] Freeze and record the final taxon identity and alias policy.
- [ ] Record exact profiler image, analysis image, database, spike FASTA,
  manifest, and source-commit identities.

## Required definitive executions

- [ ] Run the Yachida definitive paired analysis from its sealed input.
- [ ] Run Feng and Zeller independently through the CRC definitive driver.
- [ ] Confirm native-profile semantics and canonical-input validation in every
  cohort package.
- [ ] Run both primary profiler-scale and prespecified sensitivity analyses.
- [ ] Run the continuous, detection, artificial-biomarker, disease-biomarker,
  response-atlas, community-reconstruction, and robustness stages required by
  the frozen statistical plan.
- [ ] Require top-level `SUCCESS` and checksum provenance for every cohort.
- [ ] Run the three-cohort definitive synthesis using only the three sealed
  cohort packages.
- [ ] Run the final publication report from the definitive synthesis.

## Manuscript and source-data package

- [ ] Regenerate every main and supplementary figure from the same definitive
  run family.
- [ ] Export a machine-readable source table for every plotted panel.
- [ ] Produce the final participant, sample, and profile flow table with an
  explicit reason for every exclusion or missing profile.
- [ ] Export model formulas, covariates, contrasts, multiple-testing families,
  effect estimates, confidence intervals, raw P values, and adjusted q values.
- [ ] Check every manuscript number against a definitive source table.
- [ ] Replace provisional presentation figures and text with definitive outputs.
- [ ] Update the abstract, Results, Discussion, and limitations only after the
  final estimates are sealed.
- [ ] Preserve careful terminology: non-implanted response is not automatically
  a false positive, perturbation-sensitive is not proven biologically false,
  and successful spike recovery does not prove native biological absence.

## Public archive and sharing metadata

- [ ] Choose the public code release commit and create a signed/versioned tag.
- [ ] Record Apptainer image SHA-256 hashes and decide whether images or build
  recipes plus registries will be archived.
- [ ] Record database releases and checksum inventories without redistributing
  assets whose licenses prohibit redistribution.
- [ ] Create a non-sensitive configuration manifest replacing local paths with
  documented placeholders.
- [ ] Archive frozen manifests, spike-panel metadata, reference checksums,
  analysis policies, logs, validation reports, source-data tables, and final
  checksums.
- [ ] Write the final Data Availability and Code Availability statements,
  including public read accessions and controlled-access limitations, if any.
- [ ] Deposit the release package in the selected archive and add its DOI.
- [ ] Update `CITATION.cff`, README citation text, repository URL, preprint DOI,
  article DOI when available, authors, and release version.
- [ ] Test a clean clone using only public instructions and archived inputs.
- [ ] Have a second person reproduce at least the validation and final-figure
  stages from the release candidate.

## Repository cleanup before tagging

- [ ] Make the root README's numbered downstream sections follow the definitive
  paired route; move detailed legacy execution into a clearly labeled section.
- [ ] Consolidate duplicated status prose and retain dates only for operational
  snapshots.
- [ ] Remove or ignore private handoffs, scheduler logs, local paths, temporary
  figures, and presentation-only products from the public release.
- [ ] Decide whether the separate GUTBIOME pilot belongs on another branch or
  in a separately labeled directory; do not mix it into paper inputs.
- [ ] Review every modified and untracked file before committing.
- [ ] Run repository checks, unit tests, containerized R tests, link checks, and
  a secret/path scan on the release candidate.
- [ ] Confirm that no definitive output directory inherits a
  `DEVELOPMENT_ONLY` status and that no development product is labeled final.

## Final release gate

The revised study is release-ready only when all three upstream seals, all
three definitive cohort packages, the three-cohort synthesis, publication
report, source-data package, and archive manifest are present and their
checksums verify from a clean checkout. Slurm completion alone is never a
release gate.
