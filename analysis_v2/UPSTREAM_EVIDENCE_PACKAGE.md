# Upstream evidence package

**Status:** implemented and fixture-tested; not yet built from real cohort data
**Package version:** `upstream_evidence_v1`
**Consumes:** `upstream_seal_v2`

One immutable, checksummed package proving which upstream inputs were complete
before definitive downstream analysis began. It serves manuscript supplementary
material, a machine-readable Zenodo release, and downstream provenance checks
without reconstructing evidence from terminal history.

It is an evidence **index**, not a copy of the data: raw reads, host-cleaned
reads, reference databases, container images, scratch files and scheduler logs
are never copied. Large external assets are recorded as logical identifiers
with checksums.

Implements `UPSTREAM_EVIDENCE_PACKAGE_SPEC.md`. Where the two disagree, the
specification wins.

## Prerequisite: satisfied

The package consumes three `production_seal_v2` directories, not the three
different native seal formats. All three cohorts passed the common contract in
`UNIFIED_UPSTREAM_SEAL_SPEC.md` using the same public auditor:

| Cohort | Slurm job | Samples | Independent | State |
|---|---|---:|---:|---|
| Yachida | 3097679 | 201 | 30 | COMPLETED 0:0 |
| Feng | 3097680 | 154 | 30 | COMPLETED 0:0 |
| Zeller | 3097681 | 156 | 30 | COMPLETED 0:0 |

Every audit verified each `production_seal.sha256` member, produced an empty
error log, reported `PASS`, and used the identical seven-member layout,
`sample_flow.tsv` schema, `covariate_audit.tsv` schema and canonical manifest
core (`sample_id`, `condition`, `study`, `independent_subset`). A separate
cross-cohort verification confirmed all three share the unified interface.

## Entry points

| Purpose | Path |
|---|---|
| Builder and validator | `analysis_v2/scripts/build_upstream_evidence_package.py` |
| Supplementary figure | `analysis_v2/scripts/plot_upstream_evidence.R` |
| Container entry point | `analysis_v2/run_upstream_evidence_package.sh` |
| Fixture suite | `analysis_v2/tests/test_upstream_evidence_package.py` |

Every authoritative input is an explicit named argument. Nothing is discovered
by searching `work/`, taking the newest directory, or inferring a cohort from a
filename.

## What is validated before anything is sealed

**Seal integrity.** Every checksum listed by every source seal is verified.
Absolute and parent-traversal member names are rejected before a file is
opened. Each `SUCCESS` must declare the `upstream_seal_v2` contract, the right
cohort, `PASS`, and the frozen counts. No `AUDIT_IN_PROGRESS` may sit beside a
seal. The assembly-sensitivity seal must carry `SUCCESS` and
`experiment_inputs_and_summary.sha256`.

**Seal membership.** Each `production_seal_v2` must hold exactly the seven
documented regular files, no subdirectory, and a `production_seal.sha256`
covering exactly the other six — an omitted *or* additional entry fails. The
assembly seal's `SUCCESS` must be checksummed. The historical post-seal
`matched_seed_audit.tsv` sidecar is required and audited independently for its
30 x 2 x 6 design, matched seeds and `PASS` statuses; it is explicitly recorded
as `AUDITED_UNSEALED`. Every other unlisted member is rejected.

**Identity and topology.** Cohort sizes exactly 201 / 154 / 156; 30 independent
samples each, balanced 10 Control / 10 Adenoma / 10 CRC; sample **sets**
compared exactly against the frozen manifests; no duplicate identifier within
or across cohorts; only `Control`, `Adenoma`, `CRC`; the canonical four columns
**first and in order** in both manifests; the common `sample_flow.tsv` and
`covariate_audit.tsv` schemas unchanged; and the `source_seal_inventory.tsv`
schema with a valid SHA-256, positive byte count and `VERIFIED` status on every
row.

Per **sample**, not merely per cohort: `sample_id`, `condition`, `study` and
`independent_subset` cross-checked between the ledger and the canonical
manifest; expected baseline 1, community 7, independent 60 for subset members
and 0 otherwise; observed equal to expected for each design;
`expected_profiles` and `observed_profiles` equal to the component sums; every
boolean audit field 1; `status` `PASS`; `failure_reasons`, `receipt_error`,
`completion_error` and `provenance_error` all blank; and
`input_provenance_mode` exactly `sealed_manifest_and_qc_receipt` for Yachida
and `state_file` for Feng and Zeller. An **offset-preserving** pair of
per-sample errors that leaves the cohort total correct still fails. The
independent manifest must equal exactly the rows flagged
`independent_subset=1`.

Assembly sensitivity must report 30 samples, 2 arms, 6 fractions per arm and
360 of 360 profiles, from the documented `SUCCESS` schema (`experiment`,
`samples`, `assembly_arms`, `fractions_per_arm`, `expected_profiles`,
`observed_profiles`, `status`), which matches the real seal.

**Audit ledger.** Exactly one row per cohort, no unknown or duplicate cohort,
`state=COMPLETED`, `exit_code=0:0`, the frozen cohort size, a numeric job
identifier equal to the expected unified audit job (Yachida 3097679, Feng
3097680, Zeller 3097681, overridable with `--<cohort>-audit-job`), a valid
explicit UTC stamp, and a `seal_sha256` equal to that cohort's actual
`production_seal.sha256`. The registry's `source_audit_job` comes only from
this validated ledger.

The `source_audit_job` inside each seal's `source_seal_inventory.tsv` describes
the earlier native/source seal, not this unified-v2 audit. Historical Yachida
therefore legitimately has consistently blank values, while Feng and Zeller
record one numeric native audit job. The builder validates that distinction
but does not require a native job to equal the unified ledger job.

**Provenance metadata.** A full 40-character commit is required. Categories
must be documented ones, `asset_id` values unique, and file-backed assets need
a real 64-character SHA-256. `upstream_parameter` rows record settings rather
than files, so the documented policy is that their `sha256` is the sentinel
`not_a_file`; a blank cell is refused. A `source_code` row must name the
supplied commit, there must be one `spike_reference` per implanted target, and
the frozen parameter names `read_length`, `bracken_threshold`,
`profiler_threads`, `spike_fractions`, `pool_coverage` and `seeds` must all be
present. The `container_analysis` digest must equal the analysis-image SHA-256
the runner passes as `--analysis-image-sha256`, which is recorded in
`provenance/build_parameters.tsv`. Rows are screened for credentials and
private keys.

**Panel and aliases.** Both inputs are parsed and validated, not copied
blindly. The frozen alias file is comma-delimited, so it is read by its real
delimiter and rewritten as a genuine tab-delimited
`zenodo/taxon_aliases.tsv`. The spike panel must carry ten unique targets with
nonblank taxon and assembly identities and positive numeric weights.

**Privacy.** No product may contain a cluster-absolute path; the whole staged
tree is screened, with URLs excluded so legitimate asset locations survive.

**Mid-build mutation.** Before first parsing any seal, the builder snapshots
every authoritative seal and revalidates it immediately before granting
`SUCCESS`, detecting changed, removed and newly added files and an
`AUDIT_IN_PROGRESS` appearing mid-run.
The shell runner independently compares complete before/after inventories, so
a newly added file is caught there too.

The final `SUCCESS` declaration is protected by both `MANIFEST.tsv` and
`SHA256SUMS`; only the checksum/index files are self-exempt where necessary.

## Privacy: individual covariates

Sample identifiers are the pseudonymous accessions already in the sealed
manifests. Individual age, sex and BMI are **not** released by default:

- `zenodo/covariate_missingness.tsv` is aggregate counts only;
- the manifests projected into `source_seal_projections/` keep only
  `sample_id, condition, study, independent_subset` and are named
  `production_manifest.canonical.tsv`, because the sealed originals carry
  `source_age`, `source_sex` and `source_bmi`.

Releasing individual covariates requires **both**
`--include-individual-covariates` and `--redistribution-review <reference>`,
and the reference is recorded in `provenance/build_parameters.tsv`. Neither
flag has a default.

## Layout

```text
upstream_evidence_<UTC timestamp>/
├── README.md, MANIFEST.tsv, SHA256SUMS, SUCCESS
├── manuscript_supplement/         figure (PDF/SVG/PNG), figure_source_data.tsv,
│                                  two supplementary tables
├── zenodo/                        release tables, incl. cohort_registry.tsv
├── source_seal_projections/       release-safe projection per component
└── provenance/                    build parameters, input checksums, commit,
                                   validation report
```

### Packaged seal copies

A packaged directory never carries a checksum manifest it cannot satisfy.
By default each component gets an explicitly named **projection**, not
something presented as a source seal:

```text
source_seal_projections/<component>/
├── README.md                            states this is NOT the original seal
├── projection.sha256                    covers exactly the files present
├── original_source_seal_inventory.tsv   every original member, SHA-256, bytes
└── <release-safe projected or copied members>
```

The original `production_seal.sha256` and
`experiment_inputs_and_summary.sha256` are deliberately **not** copied into a
projection, because they index members a projection does not contain and would
fail to verify there. Original hashes remain in each projection's inventory,
in `zenodo/source_seal_inventory.tsv` and in
`provenance/input_checksums.tsv`. The builder verifies every projection from
inside its own directory before granting `SUCCESS`.

With `--include-individual-covariates` and a recorded
`--redistribution-review`, a complete **byte-identical** seal is copied to
`source_seals/<component>/` instead, containing every original member with its
original manifest verifying unchanged and nothing else added beside it. The
two semantics are never mixed in one package.

`zenodo/cohort_registry.tsv` is the machine-readable one-row-per-cohort index:
`cohort, study, samples, independent_samples, baseline_profiles,
community_profiles, independent_profiles, seal_contract, seal_status,
input_provenance_mode, source_audit_job`. It carries no cluster path.

## Atomicity and sealing

The build refuses a nonempty final output directory and a leftover
`<outdir>.incomplete`. Everything is staged in `<outdir>.incomplete`; any
failure removes that staging directory, so a partial package is never left to
be mistaken for a finished one. `MANIFEST.tsv` (relative path, bytes, SHA-256,
media type, release role) and `SHA256SUMS` are written only after every
validation and the figure have passed, `SHA256SUMS` is rechecked entry by
entry, and only then is the staging directory atomically renamed into place.
`SUCCESS` records the package version, source commit, three cohort sizes, the
assembly-sensitivity profile count, the file count and `PASS`.

Verify a built package with `sha256sum -c SHA256SUMS` from its root.

## Supplementary figure

One four-panel figure, **Upstream cohort and profile completeness**: (A)
samples by cohort and condition; (B) expected and observed profiles by cohort
and design; (C) independent-subset samples with the expected value of 10
marked; (D) aggregate covariate missingness. Every plotted value is read from
`manuscript_supplement/figure_source_data.tsv`, so the figure and the released
tables cannot disagree. Colours are Okabe-Ito, exact counts are printed on
every bar, axes use integer breaks, and no cluster path or timestamp is
embedded. PDF and SVG are the publication masters; PNG is a review copy. A
draft caption is in the package README and states that completeness of
analytical inputs does not establish biological accuracy.

The plotter composes panels with `patchwork` or `gridExtra` when present and
otherwise with base `grid`, so the figure does not depend on an optional
package.

## Running it

```bash
export ANALYSIS_SIF=/path/to/pinned/analysis.sif
export OUTDIR="$PWD/work/upstream_evidence_$(date -u +%Y%m%dT%H%M%SZ)"
export YACHIDA_SEAL="$PWD/work/yachida_strict_final_20260823/state/production_seal_v2"
export FENG_SEAL="$PWD/work/feng_strict_final_20260911/state/production_seal_v2"
export ZELLER_SEAL="$PWD/work/zeller_strict_final_20260911/state/production_seal_v2"
export ASSEMBLY_SENSITIVITY_SEAL="$PWD/work/yachida_assembly_sensitivity_20260901/experiment_seal"
export YACHIDA_MANIFEST=...            # and the five other frozen manifests
export PROVENANCE_METADATA=...         # upstream_provenance_metadata.tsv
export AUDIT_LEDGER=...                # manually supplied operational evidence
bash analysis_v2/run_upstream_evidence_package.sh
```

The runner requires `ANALYSIS_SIF`, computes the image checksum and passes it
to the builder for cross-checking against `container_analysis`, and compares
complete before/after seal inventories, failing if any file changed, vanished
or appeared.

Templates for the two manually authored inputs are in
`analysis_v2/templates/upstream_evidence/`: schemas and placeholders, never
invented hashes.

## Local verification

```bash
python3 analysis_v2/tests/test_upstream_evidence_package.py
python3 -m py_compile analysis_v2/scripts/build_upstream_evidence_package.py
bash -n analysis_v2/run_upstream_evidence_package.sh
```

The suite drives the real builder against synthetic seals shaped exactly like
the real ones and covers all 15 specification cases plus the additional gates.
It uses production-scale cohort sizes, so the frozen counts are genuinely
exercised.

## Definition of done

Not complete until the fixture tests pass on the cluster, Codex finds no
weakening of source seals or privacy rules, the real package builds in a new
directory, every checksum verifies, the aggregate counts match the sealed
evidence exactly, and the figure and its source data are manually reviewed.
Generated real-data packages stay outside Git for later Zenodo deposition;
only code and small schemas are committed.

**As of this commit no real-data package has been built.**
