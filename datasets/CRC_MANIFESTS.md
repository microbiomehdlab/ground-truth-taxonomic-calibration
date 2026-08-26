# Frozen Feng and Zeller production manifests

The version-controlled files under `datasets/fengq/manifests/` and
`datasets/zellerg/manifests/` are the complete public input specifications for
the final cohort runs. They contain public study metadata, ENA accessions and
download metadata, deterministic independent-subset membership, provenance,
and SHA-256 checksums. They do not contain raw reads, credentials, host-specific
configuration, or private filesystem paths.

For each cohort:

- `eligible_cohort.tsv` freezes metadata eligibility;
- `eligible_cohort.excluded.tsv` retains every exclusion and reason;
- `independent_10_per_condition.tsv` freezes the nested independent subset;
- `selection_audit/` records deterministic-selection and balance checks;
- `PRJ*.read_run.tsv` is the official ENA Portal API snapshot retrieved on
  2026-08-26;
- `production_manifest.tsv` has one row per biological sample;
- `production_manifest.runs.tsv` has one row per paired sequencing run;
- `production_manifest.independent.tsv` contains exactly the selected 10
  samples per condition and is consumed by the strict production runner;
- `production_manifest.provenance.tsv` records all inputs and frozen rules.

Rebuild every derived manifest from the committed metadata and ENA snapshots:

```bash
bash datasets/build_crc_public_manifests.sh
```

The rebuild is deterministic. Afterward, `git diff --exit-code` must report no
changes. The ENA snapshots are committed so that future changes in the external
API cannot silently alter the frozen study inputs.

Machine-specific files remain private and ignored:

- `config/*.env`;
- raw FASTQs and disposable scratch;
- scheduler logs and state;
- absolute cluster paths;
- private validation outputs under `work/`.
