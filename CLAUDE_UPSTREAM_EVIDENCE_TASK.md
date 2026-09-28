# Claude implementation brief: upstream evidence package

> **PREREQUISITE SATISFIED (28 September 2026).** The unified three-cohort
> seal in `UNIFIED_UPSTREAM_SEAL_SPEC.md` is implemented, independently
> reviewed, and has run successfully against all three real cohorts with the
> same public auditor and schema: Yachida job **3097679** (201 samples, 30
> independent), Feng job **3097680** (154, 30) and Zeller job **3097681**
> (156, 30), each COMPLETED 0:0, each verifying every `production_seal.sha256`
> member, each with an empty error log and `PASS`, and each using the identical
> seven-member layout, `sample_flow.tsv` schema, `covariate_audit.tsv` schema
> and canonical manifest core. A separate cross-cohort verification confirmed
> all three share the unified upstream interface. This task therefore consumes
> the three canonical `production_seal_v2` directories and never parses the
> native seal formats.

Implement the frozen design in `UPSTREAM_EVIDENCE_PACKAGE_SPEC.md`.

## Scope

Write code, documentation, and fixture tests only. Do not run on the real
cluster, modify source seals, alter upstream outputs, fit downstream models, or
change the scientific design. Do not include the separate GUTBIOME pilot.

Read first:

1. `UPSTREAM_EVIDENCE_PACKAGE_SPEC.md`;
2. `REPRODUCING.md`;
3. `PUBLIC_RELEASE_CHECKLIST.md`;
4. `analysis_v2/CRC_COHORT_DEFINITIVE_RUNBOOK.md`;
5. `analysis_v2/YACHIDA_DEFINITIVE_RUNBOOK.md`;
6. `datasets/yachida/audit_production.py`;
7. `analysis_v2/scripts/seal_crc_cohort_upstream.py`;
8. `analysis_v2/METHODS_DECISION_LOG.md`.

## Deliverables

Create:

```text
analysis_v2/scripts/build_upstream_evidence_package.py
analysis_v2/scripts/plot_upstream_evidence.R
analysis_v2/run_upstream_evidence_package.sh
analysis_v2/tests/test_upstream_evidence_package.py
analysis_v2/UPSTREAM_EVIDENCE_PACKAGE.md
```

Update only the minimum navigation and decision documents needed to expose the
new workflow and record implementation decisions. Do not rewrite unrelated
documentation.

## Non-negotiable behavior

- All authoritative paths are explicit named CLI arguments.
- Verify the source seals before reading their scientific tables.
- Compare exact sample identities, not only counts.
- Treat `SUCCESS` as an existence sentinel where the producing workflow uses
  `touch`; require substantive tables to be nonempty.
- Do not expose absolute cluster paths in release-facing products.
- Do not publish individual covariate values by default.
- Never copy raw reads, databases, containers, scratch, or unrestricted logs.
- Use an atomic staging directory and refuse overwrite.
- Produce `MANIFEST.tsv`, `SHA256SUMS`, and final `SUCCESS` only after all
  validations and plots pass.
- Preserve source seals byte-for-byte on success and failure.
- Keep generated real-data packages outside Git.
- Make the plot consume exported figure source data rather than recomputing
  aggregates independently.

## Testing expectations

Implement all 15 tests listed in the specification where practical in one
dependency-free Python fixture suite. The test must exercise the real builder,
not a duplicate implementation. Plotting may be skipped only when R/ggplot2 is
unavailable, with an explicit skip; structural source-data validation must
still run.

Run and report:

```bash
python3 analysis_v2/tests/test_upstream_evidence_package.py
python3 -m py_compile analysis_v2/scripts/build_upstream_evidence_package.py
bash -n analysis_v2/run_upstream_evidence_package.sh
git diff --check
```

Do not claim the real package is built or publication-ready. End with a concise
implementation summary, test output, changed-file list, and any questions that
require scientific judgment. Codex will independently review the patch and
prepare cluster execution commands.
