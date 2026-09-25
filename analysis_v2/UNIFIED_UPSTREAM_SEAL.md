# Unified upstream cohort seal

**Status:** implemented, not yet executed against real cohort outputs
**Contract:** `upstream_seal_v2`
**Scope:** Yachida, Feng, Zeller

One configurable auditor validates all three paper cohorts against the same
scientific and integrity contract and emits one canonical seal layout. It is a
validation and interface migration: it never recomputes a profile, never
rewrites a retained output, and never modifies a native seal.

Implements `UNIFIED_UPSTREAM_SEAL_SPEC.md` (frozen 24 September 2026). Where
this document and that specification disagree, the specification wins.

## Entry points

| Purpose | Path |
|---|---|
| Auditor | `analysis_v2/scripts/seal_cohort_upstream.py` |
| Cluster job | `analysis_v2/run_cohort_upstream_audit.sbatch` |
| Fixture suite | `analysis_v2/tests/test_unified_upstream_seal.py` |

The native auditors `datasets/yachida/audit_production.py` and
`analysis_v2/scripts/seal_crc_cohort_upstream.py` are unchanged and remain
available for provenance. They produce the **native source seals** this
workflow consumes; they are not replaced and must not be rerun to satisfy it.

## What is validated

For every manifest sample:

- exactly one baseline profile, exactly seven community profiles, and exactly
  sixty independent profiles when the sample is in the frozen independent
  subset (otherwise exactly zero);
- a sample-level success sentinel, and a nonempty verified marker,
  retained-output receipt, input-provenance table and sample-completion table;
- exact agreement between manifest identity, condition, subset membership,
  expected topology, observed topology and completion metadata;
- every retained output exists, is nonempty, matches its recorded byte count
  and SHA-256, lies outside disposable scratch, and lies under that sample's
  permitted results or QC roots.

Each topology component is checked separately, never only the total. A
community profile miscounted as a baseline profile leaves the total unchanged
and is still rejected.

For every cohort: unique production identifiers, an exactly nested independent
subset of 30 samples balanced 10 Control / 10 Adenoma / 10 CRC, explicit
expected sample and condition counts, exact aggregate topology, and a verified
native source seal whose manifest copies are byte-identical to the
authoritative manifests.

Frozen expectations, unchanged by this workflow:

| Cohort | Samples | Conditions | Baseline | Community | Independent |
|---|---:|---|---:|---:|---:|
| Yachida | 201 | `Control=67,Adenoma=67,CRC=67` | 201 | 1,407 | 1,800 |
| Feng | 154 | `Control=61,Adenoma=47,CRC=46` | 154 | 1,078 | 1,800 |
| Zeller | 156 | `Control=61,Adenoma=42,CRC=53` | 156 | 1,092 | 1,800 |

`--expected-conditions` is required. When `--expected-samples` equals a frozen
cohort size the supplied counts are cross-checked against the table above and a
disagreement is refused, so a production run cannot be sealed against drifted
expectations. Scaled fixtures, which never use the production sizes, are
unaffected.

## Configuration, not cohort-specific science

Nothing is discovered by search. Every manifest, seal and root is named
explicitly. Cohort differences are boundary settings only:

| Setting | Yachida | Feng / Zeller |
|---|---|---|
| Condition column | `Target_Condition` | `condition` |
| Study | fixed `YachidaS_2019` | `study` column |
| Batch column | `batch_id` | none (blank in output) |
| Completion table has `study` | no | yes |
| Covariates audited | `age`, `sex`, `bmi`, `batch_id` | `age`, `sex`, `bmi` |
| Production-only columns | batch provenance | none |
| Independent header schema | `yachida_historical_duplicate_selection` | `unique` |
| Native SUCCESS schema | `yachida_dataset` | `cohort_profiles` |
| Seal manifest member | `pilot_batched.tsv` | `production_manifest.tsv` |
| Seal subset member | `independent_10_per_condition.tsv` | `production_manifest.independent.tsv` |
| Native audit job | not recorded | Feng `3097587`, Zeller `3097588` |

Defaults live in the adapter table and every one can be overridden on the
command line. The auditor refuses unknown cohorts, ambiguous study
configuration (both or neither of a study column and a fixed study value),
missing required columns, unexpected condition counts, duplicate manifest
columns, duplicate receipt paths, relative paths, an output directory that
overlaps the native seal or scratch, and an output directory holding
unexpected members.

## Canonical seal

```text
production_seal_v2/
├── sample_flow.tsv
├── covariate_audit.tsv
├── production_manifest.tsv
├── production_manifest.independent.tsv
├── source_seal_inventory.tsv
├── SUCCESS
└── production_seal.sha256
```

`sample_flow.tsv` uses one schema for every cohort; `batch_id` is always
present and blank where no upstream batch exists. `covariate_audit.tsv` holds
aggregate missingness only, read through the explicit original-to-canonical
column mapping so a field renamed to `source_<name>` is still audited against
its real values. It always covers `condition`, `study` and
`independent_subset`, plus the cohort's configured covariates: `age`, `sex` and
`bmi` for Feng and Zeller, and `age`, `sex`, `bmi` and `batch_id` for Yachida,
all four of which the sealed 29-column production manifest carries.
Unavailable fields are never invented; use
`--covariate-column` when a manifest genuinely gains one. Each audited field is
named for a real column of `production_manifest.tsv`, so a covariate whose name
collides with a canonical one is reported as `source_<name>` from its own
values instead of being shadowed by the normalized column. The canonical manifest copies normalize
`sample_id`, `condition`, `study` and `independent_subset` and preserve every
original column under a collision-safe `source_<name>` alias; the independent
manifest is an exact row-consistent subset. `source_seal_inventory.tsv` lists
the native seal's logical members, hashes, sizes, verification status and the
recorded native audit job, with no release-facing absolute paths.

The independent manifest does **not** have to repeat the production manifest's
header exactly. It may carry additional deterministic-selection provenance
columns, which are preserved collision-safely in the canonical copy. Required
columns are also split: the sealed 26-column Yachida independent manifest
genuinely carries no `batch_id` or other batch field, so production-only batch
provenance (`batch_id`, `batch_hash`, `batch_position`, `batch_size`,
`batch_seed`, `processing_order`) is required in the production manifest only.
Everything else needed for identity and scientific metadata must be present in
both, and every column the two share must hold exactly the same value row by
row. A missing required column, an altered shared value, an unexpected
duplicate header, or a sample that is not nested in the production manifest all
fail closed.

### The historical Yachida independent header

The sealed Yachida independent manifest holds `selection_rank`,
`selection_hash` and `selection_seed` **twice**: the first occurrence is the
pilot value it inherited, the second was appended by the independent-subset
selection step. That file is checksummed and is never rewritten.

The Yachida adapter's `yachida_historical_duplicate_selection` schema accepts
**exactly two formats**:

| Format | Header | Handling |
|---|---|---|
| A, historical sealed | the bare triplet exactly twice each | read positionally and translated by occurrence: first → `pilot_selection_*`, second → `independent_selection_*` |
| B, corrected selector | exactly one `pilot_selection_*` triplet and one `independent_selection_*` triplet, no bare `selection_*` | already unambiguous, used as written |

Everything else is refused: a **single bare triplet** (ambiguous provenance),
three or more occurrences, an incomplete duplicated triplet, a partial `pilot_*`
or `independent_*` triplet, a mixture of bare and prefixed selection fields,
missing pilot or missing independent provenance, and any duplicated field
outside the exact historical triplet.

This is a narrowly specified translation of one known checksummed format into
an unambiguous canonical representation — **not** tolerance of duplicate
headers. Production manifests always reject duplicates, and Feng and Zeller
independent manifests always reject duplicates. In both accepted formats the
pilot values must equal the production manifest's own selection values, the
independent values are preserved in the canonical copy, and the canonical
header is unique. Byte-identity checks against the native seal and the
authoritative manifests always use the original unmodified bytes.

`scripts/select_samples_deterministically.py` no longer creates this ambiguity.
Without a name collision its output is unchanged; with one it fails unless the
caller supplies `--existing-selection-prefix` and `--new-selection-prefix`. The
Yachida reproduction command in `datasets/yachida/README.md` now passes
`pilot` and `independent`, so a re-derived manifest is format B. The auditor
accepts format A and format B, and nothing in between.

`SUCCESS`:

```text
seal_contract	upstream_seal_v2
cohort	<cohort>
samples	<n>
independent_samples	30
baseline_profiles	<n>
independent_profiles	1800
community_profiles	<7*n>
status	PASS
```

`production_seal.sha256` covers every other member.

## Ordering and atomicity

1. All paths are resolved and overlap with the native seal is refused.
2. Any previous `SUCCESS` and `production_seal.sha256` are deleted
   **immediately** on entering the output directory, before its contents are
   inspected at all — so even an invocation that is about to be refused for an
   unexpected member cannot leave an older seal looking current. Recognized
   `<member>.tmp` leftovers are then removed, any other unexpected member is
   refused, and only then is `AUDIT_IN_PROGRESS` written.
3. The native source seal is verified: no `AUDIT_IN_PROGRESS`, a nonempty
   `SUCCESS`, a complete checksum manifest that verifies member by member, and
   manifest copies byte-identical to the authoritative manifests. The native
   `SUCCESS` is then read semantically through an explicit per-cohort schema
   adapter — `yachida_dataset` (`dataset`, `samples`, `independent_subset`,
   `profiles`, `conditions`, `status`) or `cohort_profiles` (`cohort`,
   `samples`, `independent_samples`, `baseline_profiles`,
   `independent_profiles`, `community_profiles`, `status`) — and its identity,
   counts, profile totals and `PASS` status must agree with the selected cohort
   and the explicit expectations. Every field listed for the schema is
   **required**: a missing topology field fails even though the remaining file
   checksums correctly. The Yachida `profiles` summary must carry baseline,
   independent and community, and its `conditions` summary must carry Control,
   Adenoma and CRC. Nothing is inferred from the file's shape.
4. Cohort-level and per-sample validation runs. The ledgers are written even
   when the cohort fails, so the failure is auditable.
5. The native seal is verified **again**: still no `AUDIT_IN_PROGRESS`, the
   same checksum manifest and members as the inventory taken in step 3, and the
   same byte-identical manifest copies. A seal that changed while the audit was
   running cannot be sealed over.
6. Only then is `production_seal.sha256` written, then `SUCCESS`, then
   `AUDIT_IN_PROGRESS` removed. The checksum lands before the authority marker
   it covers.

**Atomicity is per member, not per directory.** Each file is replaced in one
step by writing a temporary sibling and renaming it, so no reader sees a
partial member; the directory as a whole is *not* swapped atomically. Directory
safety comes from the ordering above: stale authority is dropped before any
validation, the checksum manifest precedes the `SUCCESS` it covers, and an
interrupted run leaves `AUDIT_IN_PROGRESS` behind with no authority files. A
crash can leave a recognized `<member>.tmp` sibling; the next run removes those
specific temporaries and refuses any other unexpected file, so the directory is
recoverable without ever trusting a partial result.

Path safety is settled before anything is created or removed: all roots are
resolved and overlap between `--outdir` and the native seal — in either
direction, including equality — is refused before `mkdir`, before any stale
authority is deleted, and before any write.

The native seal is opened read-only throughout. The fixture suite proves it is
byte-identical after successful audits, failed audits and refused invocations.

## Running it

Run on a compute node: the audit rehashes every retained output.

```bash
export PROJECT="$PWD"
export COHORT=feng
export COHORT_ENV="$PWD/config/feng.strict-production.env"
export PRODUCTION_MANIFEST="$PWD/datasets/fengq/manifests/production_manifest.tsv"
export INDEPENDENT_MANIFEST="$PWD/datasets/fengq/manifests/production_manifest.independent.tsv"
export SOURCE_SEAL="<verified native seal directory>"
export EXPECTED_SAMPLES=154
export EXPECTED_CONDITIONS="Control=61,Adenoma=47,CRC=46"
export OUTDIR="$PWD/work/feng_production_seal_v2_$(date -u +%Y%m%dT%H%M%SZ)"
sbatch --export=ALL analysis_v2/run_cohort_upstream_audit.sbatch
```

Zeller uses `COHORT=zeller`, the `datasets/zellerg` manifests,
`EXPECTED_SAMPLES=156` and `EXPECTED_CONDITIONS="Control=61,Adenoma=42,CRC=53"`.
Yachida uses `COHORT=yachida`,
`COHORT_ENV=config/yachida.strict-production.env`, the
`work/yachida_67x3/metadata` manifests, `EXPECTED_SAMPLES=201` and
`EXPECTED_CONDITIONS="Control=67,Adenoma=67,CRC=67"`. `OUTDIR` must be a new
directory and must not overlap the native seal in either direction.

Verify afterwards:

```bash
test -s "$OUTDIR/SUCCESS"
test ! -e "$OUTDIR/AUDIT_IN_PROGRESS"
(cd "$OUTDIR" && sha256sum -c production_seal.sha256)
```

## Local verification

```bash
python3 analysis_v2/tests/test_unified_upstream_seal.py
python3 -m py_compile analysis_v2/scripts/seal_cohort_upstream.py
bash -n analysis_v2/run_cohort_upstream_audit.sbatch
```

The suite drives the real auditor against synthetic cohort trees — as a
subprocess, except for the mid-audit mutation case, which wraps the
revalidation hook in process so that removing the call makes the test fail. It
covers all 18 cases in section 7 of the specification plus the Codex-review
corrections (path-safety ordering, covariate mapping, independent-manifest
selection provenance, interrupted temporaries, native-SUCCESS semantics and
source-seal revalidation). It uses scaled-down cohorts, so it proves the
contract's logic, not the production counts.

## Acceptance gate

Not complete until the fixture suite passes on the cluster, Codex reviews the
implementation, the auditor runs against the real Yachida, Feng and Zeller
outputs, all three v2 seals verify, counts and sample sets agree with the
native evidence, and no source seal or profile changed. **Only then** may
`UPSTREAM_EVIDENCE_PACKAGE_SPEC.md` be implemented.

As of this commit none of the real-cohort steps has been performed. Standard
library only, with postponed annotation evaluation for the cluster's older
Python.
