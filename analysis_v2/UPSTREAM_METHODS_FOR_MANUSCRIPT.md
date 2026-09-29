# Upstream Methods source for the revised manuscript

**Status:** upstream methods and evidence are verified; downstream analysis is pending
**Evidence source:** `upstream_evidence_v1`, job `3097746`, source commit `9df9ae7dce670ec83fbbfbb709c10a0bb51c4947`
**Purpose:** one authoritative source for drafting the main Methods, Supplementary Methods, sample-flow text, and reproducibility statements

This document translates the frozen upstream contracts into manuscript-facing language. It does not replace the machine-readable seals or evidence package. Exact software, image, database, reference, manifest, and checksum identities must be taken from the evidence package rather than copied from terminal history.

## Study design and controlled ground truth

We evaluated taxonomic profilers using in silico addition of simulated paired-end reads to previously sequenced human faecal metagenomes. This design preserves the biological and technical complexity of each background metagenome while fixing the identity of the added organism, the number of added read pairs, and the resulting read fraction. The controlled truth is therefore sequence identity, implanted-pair count, and final-library read fraction. It is not cellular abundance, biomass, extraction efficiency, or proof of biological presence or absence.

Ten colorectal-cancer-associated taxa were represented by versioned reference assemblies recorded in `spikes/spike_panel.tsv`. Reference downloads, sequence identities, pool settings, random seeds, and checksums were retained as provenance. Taxon naming and profiler-specific equivalences were resolved using the frozen alias table rather than by post hoc string matching.

## Cohorts and frozen sample sets

Three public colorectal-cancer metagenomic cohorts were included: Yachida, Feng, and Zeller. Eligibility and condition labels were frozen before evaluation of profiling or recovery outcomes. The final upstream population comprised 511 samples:

| Cohort | Control | Adenoma | CRC | Total |
|---|---:|---:|---:|---:|
| Yachida | 67 | 67 | 67 | 201 |
| Feng | 61 | 47 | 46 | 154 |
| Zeller | 61 | 42 | 53 | 156 |
| Total | 189 | 156 | 166 | 511 |

Within each cohort, a deterministic nested subset of 30 samples was selected for the independent spike experiment, balanced as 10 Control, 10 Adenoma, and 10 CRC. Selection provenance was retained in the sealed manifests. Samples were not replaced or excluded on the basis of profiler output, spike recovery, or downstream statistical significance.

## Read preprocessing and taxonomic profiling

Paired-end reads were processed through the strict production workflow using pinned software environments and reference assets. Adapter and quality processing used paired-read detection, four-base sliding-window end trimming at mean Phred 20, a minimum retained length of 60 bases, poly-G trimming, and low-complexity filtering at threshold 30. Human reads were removed by alignment with Bowtie2 in `--very-sensitive` mode. Only synchronized pairs for which both primary mates were unmapped were retained; secondary and supplementary alignments were excluded. Input provenance included run/accession identity, paired-file metadata, checksums or receipt-verified integrity records, and sample-level processing evidence.

Taxonomic profiles were generated independently from the paired compressed FASTQ files. MetaPhlAn 4 used its pinned marker database and reported bacterial profiles while excluding eukaryotes and archaea. Kraken2 classified paired reads against its pinned database, after which Bracken estimated species-level abundance using the read-length and abundance-threshold settings recorded for each production run in the evidence package. Profiler-native species abundance outputs were retained as the primary practical measurements; the workflow did not force the two profilers onto an identical biological abundance estimand.

The exact tool versions, Apptainer image digests, database releases, reference checksums, parameters, and native abundance fields used for the definitive upstream run are indexed in the verified evidence package. These machine-readable records are authoritative if abbreviated manuscript text differs from a software or database identifier.

## In silico spike construction

For a background library containing `R` synchronized read pairs and a target final fraction `f`, the number of inserted pairs was

```text
N = round[f / (1 - f) x R],
```

which gives an achieved final fraction of `N / (R + N)`. Nominal and achieved fractions, inserted-pair counts, and deterministic seeds were recorded for each generated library.

Two perturbation designs were used:

1. **Independent spikes:** each taxon was added separately at final read-pair fractions 0.01%, 0.05%, 0.1%, 0.5%, 1%, and 5% in the nested 30-sample subset of each cohort.
2. **Community spikes:** all ten taxa were added together with equal target weights at total final fractions 0.01%, 0.05%, 0.1%, 0.5%, 1%, 5%, and 10% in every cohort sample.

Randomness followed the frozen `stable-seed-v1` SHA-256 derivation. Seeds depended on the scientific identifiers relevant to each operation, including sample, target, fraction, and community member where applicable; scheduler task number, manifest position, batch size, and filesystem path were excluded. Requested fractions consumed deterministic paired-read streams, with the selected first mate defining the synchronized second mate.

## Upstream profile topology

Every retained sample was required to have exactly one unspiked baseline profile and seven community-spike profiles. Each member of the nested independent subset was additionally required to have exactly 60 independent profiles (10 taxa x 6 fractions); samples outside that subset were required to have none. The verified topology was:

| Cohort | Baseline profiles | Community profiles | Independent profiles |
|---|---:|---:|---:|
| Yachida | 201 | 1,407 | 1,800 |
| Feng | 154 | 1,078 | 1,800 |
| Zeller | 156 | 1,092 | 1,800 |
| Total | 511 | 3,577 | 5,400 |

Counts were checked separately by design so that an extra profile in one design could not compensate for a missing profile in another.

## Unified validation and sealing

One configurable auditor applied the same scientific and integrity contract to all three cohorts and emitted the common `upstream_seal_v2` interface. Cohort differences—such as source condition-column names, fixed versus manifest-provided study labels, batch metadata, historical independent-selection headers, and input-provenance storage—were handled as explicit adapters rather than separate scientific workflows.

For every sample, the auditor cross-checked manifest identity, condition, independent-subset membership, expected and observed profile topology, completion metadata, success markers, input provenance, and retained-output receipts. Every retained file had to exist, be nonempty where appropriate, match its recorded byte count and SHA-256, lie outside disposable scratch, and fall under an explicitly permitted results or QC root. The independent subset had to be nested exactly within its production cohort and retain the frozen 10/10/10 condition balance.

The unified audits completed successfully for Yachida (job `3097679`; 201 samples), Feng (`3097680`; 154 samples), and Zeller (`3097681`; 156 samples). All returned exit code `0:0`, empty error logs, checksum-verifying manifests, and `PASS`. The canonical output schema is shared across cohorts and contains aggregate covariate missingness rather than individual covariate values by default.

## Assembly-choice sensitivity

A prespecified Yachida sensitivity experiment compared the original and clean replacement assemblies while holding the selected samples, dose grid, and deterministic seeds fixed. It used the frozen 30-sample independent subset, two assembly arms, and six fractions per arm, yielding 360 expected and 360 observed profiles. The matched-seed audit contains 126 applicable audited combinations; every recorded original/clean seed pair agreed and every row passed. This audit file post-dates the historical checksum manifest and is therefore reported explicitly as `AUDITED_UNSEALED`, with its own package-recorded digest, rather than being represented as covered by that manifest.

This sensitivity analysis evaluates dependence on assembly choice. It does not by itself isolate contamination, strain divergence, reference representation, or database coverage as causal mechanisms.

## Evidence package and reproducibility

The three common cohort seals and the assembly-sensitivity seal were assembled into `upstream_evidence_v1`. Cluster job `3097746` completed successfully and produced a top-level `PASS` package containing 56 files and 57 independently verified checksum entries. The package records source commit `9df9ae7dce670ec83fbbfbb709c10a0bb51c4947`, the analysis-container digest, software and database identities, reference and manifest checksums, cohort/sample/profile flow, covariate missingness summaries, audit jobs, build parameters, release-safe seal projections, supplementary source tables, and figure assets.

Raw reads, host-cleaned reads, databases, container images, scratch files, unrestricted scheduler logs, and individual covariate values are not redistributed in the default package. Large external assets are identified by stable version information and checksums. Public deposition should include the release-safe package, source commit or tagged archive, code and build recipes, frozen manifests, non-sensitive provenance, and final checksum manifest, subject to source-data and database redistribution terms.

## Manuscript reporting rules

When drafting the paper:

- report the exact cohort and profile counts above from the verified seals;
- describe the spike as an artificial biomarker with known read-level dose, not a biological abundance standard;
- keep nominal and achieved spike fractions distinct;
- state that MetaPhlAn and Kraken2/Bracken provide profiler-native measurements with different database and abundance semantics;
- distinguish analytical non-detection from biological absence;
- do not call every non-implanted taxon response a false positive without an explicit null definition;
- do not infer extraction or cellular-abundance performance from an in silico read-addition experiment;
- treat assembly-choice results as sensitivity evidence rather than a single-mechanism diagnosis; and
- take exact software/database identifiers and checksums from the evidence package at submission time.

## Material still pending outside upstream Methods

Upstream computation, validation, sealing, and evidence assembly are complete. The following are not upstream gaps: definitive downstream modelling, multiple-testing and uncertainty outputs, three-cohort synthesis, final main and supplementary figures, claim-to-source-data verification, clean-clone release testing, manual visual review of packaged figure assets, and Zenodo deposition. Results and Discussion text must be updated only from those definitive downstream products.
