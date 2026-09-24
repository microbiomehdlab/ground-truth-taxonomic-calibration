# Unified upstream cohort-seal specification

**Status:** frozen implementation contract, 24 September 2026  
**Scope:** Yachida, Feng, and Zeller production outputs  
**Purpose:** one public validation and sealing workflow for all paper cohorts

## 1. Decision

The public workflow must use one configurable cohort auditor that validates
the same scientific and integrity contract for all three cohorts and emits one
canonical seal schema.

This is a validation/interface migration only. It must not recompute profiles,
rewrite retained outputs, or replace the already verified native seals. Native
seals remain immutable provenance inputs. The unified seal becomes the input
contract for the evidence package and definitive downstream workflows.

## 2. Shared scientific contract

For every manifest sample, require:

- exactly one baseline profile;
- exactly seven community profiles;
- exactly sixty independent profiles when the sample belongs to the frozen
  independent subset, otherwise zero;
- a sample-level success sentinel;
- a nonempty verified marker, retained-output receipt, input-provenance table,
  and sample-completion table;
- exact agreement between manifest identity, condition, independent-subset
  membership, expected topology, observed topology, and completion metadata;
- every retained output exists, is nonempty, has the recorded byte count and
  SHA-256 digest, is outside disposable scratch, and is under that sample's
  explicitly permitted results or QC roots.

For every cohort, require unique production-manifest sample identifiers, an
exact nested independent subset of 30 samples balanced 10 Control, 10 Adenoma,
and 10 CRC, explicit expected cohort and condition counts, exact aggregate
profile topology, and verified authoritative manifests or a verified native
source seal containing their byte-identical copies.

Frozen expectations:

| Cohort | Samples | Conditions | Baseline | Community | Independent |
|---|---:|---|---:|---:|---:|
| Yachida | 201 | 67/67/67 | 201 | 1,407 | 1,800 |
| Feng | 154 | frozen manifest counts | 154 | 1,078 | 1,800 |
| Zeller | 156 | frozen manifest counts | 156 | 1,092 | 1,800 |

## 3. Configuration, not cohort-specific scientific logic

The common command must receive explicit named arguments. No recursive search
for manifests, seals, or cohort roots is allowed. Configuration must cover the
cohort identifier; production and independent manifests; native source seal;
state, results, QC, scratch, and output roots; expected sample and condition
counts; condition-column name; study-column name or fixed study value; optional
batch-column name; and whether source completion tables contain a study field.

Yachida uses `Target_Condition`, a fixed study value of `YachidaS_2019`, and an
optional `batch_id`. Feng and Zeller use `condition` and `study`. These are
adapter settings at the boundary, not separate scientific implementations.

The command must refuse unknown cohorts, ambiguous study configuration,
missing required columns, unexpected conditions, duplicate fields, duplicate
receipt paths, and paths outside permitted roots.

## 4. Canonical unified seal

Every cohort must emit exactly this top-level layout:

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

`sample_flow.tsv` must use one schema for all cohorts and include at least:

```text
sample_id study condition independent_subset
batch_id
expected_baseline_profiles observed_baseline_profiles
expected_independent_profiles observed_independent_profiles
expected_community_profiles observed_community_profiles
expected_profiles observed_profiles retained_output_files
verified_marker retained_output_receipt input_provenance sample_success
baseline_profile_count independent_profile_count community_profile_count
completion_table manifest_independent_flag status failure_reasons
receipt_error completion_error
```

`batch_id` is always present and is blank for cohorts without upstream batches,
so the sample-flow schema remains identical. `covariate_audit.tsv` contains
aggregate missingness only.
Canonical manifest copies normalize `condition`, `study`, and
`independent_subset` while preserving original columns with collision-safe
names. The independent manifest must be an exact row-consistent subset.

`source_seal_inventory.tsv` records logical source-seal members, hashes, sizes,
and verification status without release-facing absolute paths.

`SUCCESS` uses the same schema for every cohort:

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

`production_seal.sha256` covers every other top-level member. `SUCCESS` and
the checksum are created only after all validation succeeds.

## 5. Source-seal handling

- Verify every checksum in the native source seal before trusting it.
- Require its success sentinel and reject `AUDIT_IN_PROGRESS`.
- Never modify, rename, or repair a native seal.
- Snapshot hashes and prove byte identity after success and failure tests.
- Record Feng audit job `3097587` and Zeller job `3097588`; do not infer a
  Yachida job identifier.

Native format differences are allowed only in the adapter/parser layer. All
cohorts must pass the same retained-output, topology, identity, and canonical
schema checks.

## 6. Public entry points

Create:

```text
analysis_v2/scripts/seal_cohort_upstream.py
analysis_v2/run_cohort_upstream_audit.sbatch
analysis_v2/tests/test_unified_upstream_seal.py
analysis_v2/UNIFIED_UPSTREAM_SEAL.md
```

The public runbook must use these common entry points for all cohorts. Existing
native audit scripts remain available for provenance. They may become
compatibility wrappers only if their existing CLI behavior and tests remain
intact. The common auditor uses only the Python standard library and supports
the cluster's older Python through postponed annotation evaluation.

## 7. Required fixture tests

The suite must execute the real common auditor and prove:

1. successful Yachida configuration;
2. successful Feng configuration;
3. successful Zeller configuration;
4. identical seal member, sample-flow, covariate, and `SUCCESS` schemas across
   modes;
5. normalization of `Target_Condition` versus `condition`;
6. fixed study value versus manifest study column;
7. exact condition counts and 10/10/10 subset balance;
8. missing baseline failure;
9. community-to-baseline offset failure despite unchanged total;
10. missing or extra independent profile failure;
11. sample-completion mismatch failure;
12. missing input provenance failure;
13. receipt path, size, digest, duplicate, and scratch failures;
14. corrupt, stale, or incomplete native source-seal failure;
15. exact manifest/seal identity mismatch failure;
16. stale unified authority invalidated before a failed rerun;
17. native source seals remain byte-identical after success and failure;
18. emitted checksum manifest verifies completely.

## 8. Acceptance gate

Implementation is not complete until the fixture suite passes locally and on
the cluster, Codex independently reviews it, the common auditor runs against
the real Yachida, Feng, and Zeller outputs, all three v2 seals verify, counts
and sample sets agree with native evidence, and no source seal or profile was
changed. Only then may `UPSTREAM_EVIDENCE_PACKAGE_SPEC.md` be implemented.
