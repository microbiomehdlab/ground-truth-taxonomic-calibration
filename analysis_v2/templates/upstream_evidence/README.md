# Upstream evidence package — manual input templates

Two inputs to `analysis_v2/scripts/build_upstream_evidence_package.py` cannot
be derived from the seals and must be authored by hand. These are **schema
templates with placeholders**, not usable inputs: every `<...>` must be
replaced with a real value. No hash here is invented or real.

Copy them outside Git, fill them in, and pass the copies as
`--provenance-metadata` and `--audit-ledger`.

## `upstream_provenance_metadata.tsv`

Columns: `category asset_id name version_or_release sha256 availability notes`.

Rules the builder enforces:

- every category must be one of `source_code`, `container_upstream`,
  `container_analysis`, `database_kraken2`, `database_metaphlan`,
  `host_reference`, `spike_reference`, `upstream_parameter`; an unknown
  category is refused rather than passed through;
- `asset_id` values are unique;
- `asset_id`, `name`, `version_or_release` and `availability` are nonblank;
- **file-backed** categories need a real 64-character hexadecimal `sha256`.
  For a directory asset, hash a committed checksum inventory, never ad hoc
  directory metadata;
- **`upstream_parameter`** rows record frozen settings rather than files, so
  their `sha256` must be the documented sentinel `not_a_file`. A blank cell is
  refused, so the policy is explicit rather than silent;
- required parameter names: `read_length`, `bracken_threshold`,
  `profiler_threads`, `spike_fractions`, `pool_coverage`, `seeds`;
- one `spike_reference` row per implanted target in the validated spike panel,
  keyed by the panel label;
- a `source_code` row must contain the exact 40-character commit the builder
  was given;
- the `container_analysis` digest must equal the runner's
  `--analysis-image-sha256`.

## `audit_ledger.tsv`

Columns: `cohort audit_job_id audit_date_utc state exit_code samples
seal_sha256 notes`.

Exactly one row for each of `yachida`, `feng` and `zeller`. The builder
requires `state=COMPLETED`, `exit_code=0:0`, the frozen cohort size, a numeric
job identifier matching the expected unified audit job (Yachida 3097679, Feng
3097680, Zeller 3097681, overridable with `--<cohort>-audit-job`), a valid UTC
`YYYY-MM-DD` or `YYYY-MM-DDThh:mm:ssZ` stamp, and a `seal_sha256` equal to the
actual SHA-256 of that cohort's `production_seal.sha256`:

```bash
sha256sum "$YACHIDA_SEAL/production_seal.sha256"
```

The registry's `source_audit_job` comes only from this validated ledger.
