# Claude implementation brief: unified upstream cohort seal

Implement the frozen contract in `UNIFIED_UPSTREAM_SEAL_SPEC.md`.

## Scope

Implement code, documentation, and fixture tests only. Do not run against real
cluster outputs, recompute profiles, modify native seals, change scientific
expectations, implement the evidence package, or touch the separate GUTBIOME
pilot. Do not commit or push.

Read completely before editing:

1. `UNIFIED_UPSTREAM_SEAL_SPEC.md`;
2. `datasets/yachida/audit_production.py`;
3. `analysis_v2/scripts/seal_crc_cohort_upstream.py`;
4. `run_yachida_production_audit.sbatch`;
5. `run_crc_cohort_upstream_audit.sbatch`;
6. `analysis_v2/YACHIDA_DEFINITIVE_RUNBOOK.md`;
7. `analysis_v2/CRC_COHORT_DEFINITIVE_RUNBOOK.md`;
8. `analysis_v2/METHODS_DECISION_LOG.md`.

## Deliverables

Create:

```text
analysis_v2/scripts/seal_cohort_upstream.py
analysis_v2/run_cohort_upstream_audit.sbatch
analysis_v2/tests/test_unified_upstream_seal.py
analysis_v2/UNIFIED_UPSTREAM_SEAL.md
```

Update only minimal navigation, runbook, and decision documentation. Preserve
both native auditors unless a fully backward-compatible wrapper is justified
and tested.

## Non-negotiable behavior

- One implementation validates all modes through explicit configuration.
- Verify native source seals first and preserve them byte-for-byte.
- Rehash every retained output and enforce permitted per-sample roots.
- Enforce each topology component, never only the total count.
- Compare exact identities, conditions, subset membership, and completion data.
- Normalize the output and `SUCCESS` schemas exactly as specified.
- Invalidate stale output authority before work and fail closed.
- Use atomic output creation and reject ambiguous or unsafe paths.
- Use only the Python standard library and support the cluster's older Python.
- Isolate Yachida format translation in an adapter; do not weaken checks.

## Verification

Implement all 18 fixture cases and run:

```bash
python3 analysis_v2/tests/test_unified_upstream_seal.py
python3 -m py_compile analysis_v2/scripts/seal_cohort_upstream.py
bash -n analysis_v2/run_cohort_upstream_audit.sbatch
git diff --check
```

End with changed files, exact test output, design choices, and assumptions that
require scientific judgment. Do not claim real-cohort validation. Codex will
review the patch and prepare cluster commands.
