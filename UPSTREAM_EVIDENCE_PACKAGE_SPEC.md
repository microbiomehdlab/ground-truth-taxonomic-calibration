# Upstream evidence package specification

**Status:** design frozen for implementation review
**Scope:** sealed Yachida, Feng, Zeller, and Yachida assembly-sensitivity
evidence only

**Prerequisite:** all three cohorts must first pass the common contract in
`UNIFIED_UPSTREAM_SEAL_SPEC.md`. This package consumes the resulting three
`production_seal_v2` directories. Native cohort-specific seals remain
provenance evidence, not separate scientific input interfaces.
**Out of scope:** taxonomic re-profiling, downstream model fitting, and any
change to a scientific estimand

## 1. Purpose

Build one immutable, checksummed package that proves which upstream inputs were
complete before definitive downstream analysis began. The package must support
three uses without reconstructing evidence from terminal history:

1. manuscript supplementary figures and tables;
2. a machine-readable Zenodo release;
3. downstream readiness and provenance checks.

The package is an evidence index, not a duplicate of raw FASTQs, reference
databases, containers, or every native profile. It may copy compact sealed
tables and must record checksums and stable identifiers for large external
assets.

## 2. Authoritative inputs

The builder must receive every path explicitly. It must not recursively search
`work/`, select the newest directory, infer a cohort from a filename, or use a
hard-coded Lobo path.

Required inputs:

- Yachida `production_seal_v2` directory;
- Feng `production_seal_v2` directory;
- Zeller `production_seal_v2` directory;
- Yachida assembly-sensitivity `experiment_seal` directory;
- the three frozen production manifests;
- the three frozen independent-subset manifests;
- `spikes/spike_panel.tsv`;
- the frozen taxon-alias/equivalence file used by the definitive workflow;
- a non-sensitive provenance metadata TSV described below;
- repository root or explicit source commit;
- a new output directory.

Every production v2 seal must contain the identical canonical member set,
`SUCCESS` schema, and `production_seal.sha256` contract defined in
`UNIFIED_UPSTREAM_SEAL_SPEC.md`. Assembly sensitivity must contain `SUCCESS` and
`experiment_inputs_and_summary.sha256`.

## 3. Required validation

The builder must fail before writing final `SUCCESS` if any condition is false.

### Seal integrity

- Verify every checksum listed by each source seal.
- Reject absolute or parent-traversal member names in checksum manifests.
- Require each `production_seal_v2` to hold exactly the seven documented
  regular files and no subdirectory, and require its `production_seal.sha256`
  to cover exactly the other six: an omitted or additional entry fails.
- Require the assembly seal's `SUCCESS` to be covered by
  `experiment_inputs_and_summary.sha256`. Its one historical post-seal sidecar,
  `matched_seed_audit.tsv`, must be present with the exact documented schema,
  126 unique audited combinations using the documented taxa and fractions,
  equal original and clean seeds, and `PASS` throughout. Reject every other unlisted member.
  Record this sidecar as `AUDITED_UNSEALED`, never as checksum-covered.
- Verify every packaged projection or exact copy independently, from inside
  its own directory, before granting `SUCCESS`.
- Require the cohort identity and expected sample count in each `SUCCESS` file.
- Require no `AUDIT_IN_PROGRESS` beside any production seal.
- Reject duplicate sample identifiers within or across cohorts.
- Require exact cohort sizes: Yachida 201, Feng 154, Zeller 156.
- Require exact independent-subset size 30 per cohort and condition balance
  10 Control, 10 Adenoma, and 10 CRC.
- Require assembly sensitivity to report 30 samples, 2 arms, 6 fractions per
  arm, and 360 observed of 360 expected profiles.

### Sample and profile topology

- Yachida: 201 baseline, 1,407 community, and 1,800 independent profiles.
- Feng: 154 baseline, 1,078 community, and 1,800 independent profiles.
- Zeller: 156 baseline, 1,092 community, and 1,800 independent profiles.
- Preserve per-sample expected and observed profile counts from sealed ledgers.
- Require every sample status to be `PASS`.
- Require manifest condition labels to use only Control, Adenoma, and CRC.
- Compare manifest and seal sample sets exactly, not only their counts.
- Require the canonical manifest columns `sample_id`, `condition`, `study`, and
  `independent_subset` to be the **first four, in that order**, in both
  canonical manifests. Native column-name differences must already have been
  resolved and audited by the unified seal workflow.
- Cross-check `sample_id`, `condition`, `study` and `independent_subset` for
  every sample between `sample_flow.tsv` and the canonical production manifest.
- Require the independent manifest's sample set to equal exactly the production
  rows flagged `independent_subset=1`.
- For every sample require expected baseline 1, expected community 7, expected
  independent 60 for subset members and 0 otherwise, observed equal to expected
  for each design, `expected_profiles` and `observed_profiles` equal to the
  component sums, every boolean audit field 1, `status` `PASS`, and
  `failure_reasons`, `receipt_error`, `completion_error` and `provenance_error`
  all blank. Aggregate totals alone are insufficient: an offset-preserving pair
  of per-sample errors must still fail.
- Require `input_provenance_mode` to be exactly
  `sealed_manifest_and_qc_receipt` for Yachida and `state_file` for Feng and
  Zeller.
- Validate `source_seal_inventory.tsv` against its exact schema, requiring a
  valid SHA-256, a positive integer byte count and `VERIFIED` status on every
  row.

### Provenance and privacy

- Record a full 40-character Git commit.
- Reject blank required software, database, reference, or container identities.
- Reject obvious credentials and private keys.
- Do not copy raw FASTQs, host-cleaned FASTQs, databases, container images,
  scratch files, or unrestricted scheduler logs.
- Do not publish cluster-absolute paths in release-facing tables. Store a
  logical asset identifier, basename, release, checksum, and availability note.
- Preserve pseudonymous sample IDs. Individual covariates may only be copied to
  the Zenodo package after a separate redistribution review. Default outputs
  contain aggregate missingness and counts, not individual age, sex, or BMI.

## 4. Non-sensitive provenance metadata contract

Input file: `upstream_provenance_metadata.tsv`

Required columns:

```text
category	asset_id	name	version_or_release	sha256	availability	notes
```

Required categories:

- `source_code`;
- `container_upstream`;
- `container_analysis`;
- `database_kraken2`;
- `database_metaphlan`;
- `host_reference`;
- `spike_reference` for each implanted target;
- `upstream_parameter` for frozen read length, Bracken threshold, profiler
  threads, spike fractions, pool coverage, and relevant seeds.

If an asset is a directory, its `sha256` must identify a committed checksum
inventory rather than an ad hoc hash of directory metadata.

Additional requirements:

- `asset_id` values are unique, and an undocumented category is refused.
- `source_code`, `container_upstream`, `container_analysis`,
  `database_kraken2`, `database_metaphlan`, `host_reference` and
  `spike_reference` rows require a nonblank 64-character hexadecimal SHA-256.
- **Checksum policy for `upstream_parameter`:** those rows record frozen
  settings rather than files, so their `sha256` must be the documented
  sentinel `not_a_file`. A blank cell is refused, making the policy explicit
  rather than silent.
- A `source_code` row must identify the exact 40-character commit supplied to
  the builder.
- One `spike_reference` row per implanted target in the validated spike panel.
- The frozen parameter names `read_length`, `bracken_threshold`,
  `profiler_threads`, `spike_fractions`, `pool_coverage` and `seeds` must all
  be present.
- The `container_analysis` digest must equal the analysis-container SHA-256
  the runner passes as `--analysis-image-sha256`, which is also recorded in
  `provenance/build_parameters.tsv`.

The spike panel and taxon-alias inputs are parsed and validated, not copied
blindly. The frozen alias file is comma-delimited, so it is parsed by its real
delimiter and rewritten as a genuine tab-delimited
`zenodo/taxon_aliases.tsv`. The panel must carry ten unique implanted targets
with nonblank taxon and assembly identities and positive numeric weights.

Templates for the two manually authored inputs live in
`analysis_v2/templates/upstream_evidence/`. They contain schemas and
placeholders, never invented hashes.

## 5. Package layout

```text
upstream_evidence_<UTC timestamp>/
├── README.md
├── MANIFEST.tsv
├── SHA256SUMS
├── SUCCESS
├── manuscript_supplement/
│   ├── upstream_completeness.pdf
│   ├── upstream_completeness.svg
│   ├── upstream_completeness.png
│   ├── supplementary_sample_flow.tsv
│   ├── supplementary_profile_completeness.tsv
│   └── figure_source_data.tsv
├── zenodo/
│   ├── README.md
│   ├── cohort_sample_flow.tsv
│   ├── cohort_condition_counts.tsv
│   ├── profile_completeness.tsv
│   ├── independent_subset_balance.tsv
│   ├── covariate_missingness.tsv
│   ├── assembly_sensitivity_completeness.tsv
│   ├── spike_panel.tsv
│   ├── taxon_aliases.tsv
│   ├── upstream_provenance_metadata.tsv
│   ├── audit_ledger.tsv
│   └── source_seal_inventory.tsv
├── source_seal_projections/          # default; see section 5a
│   ├── yachida/
│   ├── feng/
│   ├── zeller/
│   └── yachida_assembly_sensitivity/
└── provenance/
    ├── build_parameters.tsv
    ├── input_checksums.tsv
    ├── repository_commit.txt
    └── validation_report.tsv
```

Only compact source-seal files are copied. Any source path recorded internally
during the build must be converted to a logical identifier in release-facing
outputs.

## 5a. Packaged seal copies (corrected 28 September 2026)

A packaged directory must never carry a checksum manifest it cannot satisfy.
The default privacy-safe output is therefore an explicitly named **projection**,
not something presented as a source seal:

```text
source_seal_projections/<component>/
├── README.md                              states this is NOT the original seal
├── projection.sha256                      covers exactly the files present
├── original_source_seal_inventory.tsv     every original member, SHA-256, bytes
└── <release-safe projected or copied members>
```

- The original `production_seal.sha256` / `experiment_inputs_and_summary.sha256`
  is **not** copied into a projection, because it indexes members the
  projection does not contain.
- Original hashes stay available through
  `zenodo/source_seal_inventory.tsv` and `provenance/input_checksums.tsv`, and
  now also through each projection's own inventory.
- `projection.sha256` must verify from inside the projection directory, and
  must cover every file there.

When `--include-individual-covariates` is explicitly approved with a recorded
`--redistribution-review`, a complete **byte-identical** seal is copied to
`source_seals/<component>/` instead. It must contain every original member, and
its original checksum manifest must verify unchanged from that directory;
nothing else is added beside it. The two semantics are never mixed in one
package.

## 6. Required tables

### `cohort_sample_flow.tsv`

One row per cohort and condition plus totals:

```text
cohort	condition	manifest_samples	sealed_samples	independent_samples	status
```

### `cohort_condition_counts.tsv`

```text
cohort	condition	samples	fraction_of_cohort
```

### `profile_completeness.tsv`

```text
cohort	design	expected_profiles	observed_profiles	completion_fraction	status
```

Design values are `baseline`, `independent`, and `community`.

### `independent_subset_balance.tsv`

```text
cohort	condition	expected_samples	observed_samples	status
```

### `covariate_missingness.tsv`

```text
cohort	field	missing	present	distinct_nonmissing
```

Aggregate only. Do not include individual covariate values by default.

### `assembly_sensitivity_completeness.tsv`

```text
experiment	samples	assembly_arms	fractions_per_arm	expected_profiles	observed_profiles	status
```

### `audit_ledger.tsv`

Manually supplied operational evidence, not reconstructed from unrestricted
logs:

```text
cohort	audit_job_id	audit_date_utc	state	exit_code	samples	seal_sha256	notes
```

The unified `upstream_seal_v2` audits this package consumes are Yachida
`3097679`, Feng `3097680` and Zeller `3097681`, each `COMPLETED` with exit code
`0:0` and an empty error log. Those identifiers are frozen in the builder and
overridable only through explicit `--<cohort>-audit-job` arguments, so an
arbitrary job identifier is never silently accepted as released provenance.
The ledger must carry exactly one row per cohort, with no unknown or duplicate
cohort, `state=COMPLETED`, `exit_code=0:0`, `samples` equal to the frozen
cohort size, a numeric job identifier, a valid explicit UTC date, and a
`seal_sha256` equal to the actual SHA-256 of that cohort's
`production_seal.sha256`. The registry's `source_audit_job` comes only from
this validated ledger.

The earlier topology-hardened CRC audits were Feng `3097587` and Zeller
`3097588` on 24 September 2026; they are the native-seal provenance recorded in
each v2 seal's inventory, not the unified audits this package validates.
Historical Yachida predates a native Slurm audit identifier, so its inventory
legitimately records blank `source_audit_job` values. Within one inventory,
these values must be consistently blank or one numeric native job; they are
not compared with the later unified audit recorded in `audit_ledger.tsv`.

### `source_seal_inventory.tsv`

```text
component	logical_file	sha256	bytes	source_seal_checksum	status
```

## 7. Supplementary figure

Create one four-panel figure titled **Upstream cohort and profile
completeness**.

- Panel A: samples by cohort and condition.
- Panel B: expected and observed profiles by cohort and design.
- Panel C: independent-subset samples by cohort and condition with the expected
  value of 10 marked.
- Panel D: covariate missingness counts by cohort and field.

Requirements:

- derive every plotted value from exported `figure_source_data.tsv`;
- display exact counts, not only percentages;
- use color-blind-safe colors and remain legible in grayscale;
- avoid implying that completion validates biological accuracy;
- write PDF and SVG as publication masters and PNG as a review copy;
- use a fixed plotting theme and explicit dimensions;
- do not embed cluster paths or timestamps in the plot;
- include a draft caption in the package README.

## 8. Implementation interfaces

Preferred implementation:

- Python builder/validator for seal parsing, tables, copying, checksums, and
  package assembly;
- R/ggplot2 plotter for the supplementary figure;
- shell entry point that runs both inside the pinned analysis container;
- dependency-free Python fixture tests for the builder and fail-closed gates;
- a small R fixture or structural output test for plotting.

Suggested files:

```text
analysis_v2/scripts/build_upstream_evidence_package.py
analysis_v2/scripts/plot_upstream_evidence.R
analysis_v2/run_upstream_evidence_package.sh
analysis_v2/tests/test_upstream_evidence_package.py
analysis_v2/UPSTREAM_EVIDENCE_PACKAGE.md
```

The builder CLI must use named arguments for every authoritative input. The
shell runner must require `ANALYSIS_SIF`, use a new output directory, record
the image checksum, and never mutate source seals.

## 9. Atomicity and sealing

- Build into `<outdir>.incomplete` or an internal staging directory.
- Refuse an existing nonempty final output directory.
- Write no final `SUCCESS` until every validation and plot has passed.
- Generate `MANIFEST.tsv` with relative path, bytes, SHA-256, media type, and
  release role for every archived file except `MANIFEST.tsv` and
  `SHA256SUMS`; the final `SUCCESS` declaration is evidence and must be
  included.
- Generate `SHA256SUMS` from relative paths, covering `MANIFEST.tsv`,
  `SUCCESS`, and every other archived file except `SHA256SUMS` itself.
- Recheck `SHA256SUMS` before atomically promoting the package.
- The builder itself, not only the shell wrapper, must snapshot every
  authoritative source seal and revalidate it immediately before granting
  `SUCCESS`, detecting changed, removed and newly added files and an
  `AUDIT_IN_PROGRESS` appearing during the build. The wrapper compares
  complete before/after inventories for the same reason: checking the original
  list alone would miss an added file.
- `SUCCESS` must record package version, source commit, three cohort sizes,
  assembly-sensitivity profile count, file count, and status `PASS`.

The downstream definitive runners may reference the package's `SUCCESS` and
`SHA256SUMS` as provenance, but this package must not replace their existing
cohort-readiness gates.

## 10. Required tests

At minimum:

1. valid three-cohort fixture succeeds with zero-byte sentinel markers;
2. missing sample fails;
3. wrong condition balance fails;
4. mismatched sample identity with correct count fails;
5. wrong profile topology fails;
6. corrupted source-seal member fails checksum validation;
7. missing or stale seal fails;
8. duplicate sample across cohorts fails;
9. incomplete assembly sensitivity fails;
10. blank provenance identity fails;
11. absolute paths do not appear in release-facing outputs;
12. attempted overwrite fails;
13. source seals remain byte-identical after success and failure;
14. tampering with a completed package makes `sha256sum -c` fail;
15. figure-source data exactly match the plotted aggregate table.

## 11. Definition of done

Implementation is complete only after:

- fixture tests pass locally;
- shell and Python syntax checks pass;
- Codex review finds no weakening of source seals or privacy rules;
- the code is committed and pulled on Lobo;
- containerized tests pass on Lobo;
- the real package builds in a new directory;
- every package checksum verifies;
- aggregate counts match the sealed evidence exactly;
- the figure and source data are manually reviewed;
- the package is committed only as code and small schemas, while generated
  evidence remains outside Git for later Zenodo deposition.
