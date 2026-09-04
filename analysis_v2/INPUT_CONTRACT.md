# Canonical analytical input contract

**Status:** frozen structural contract; final detection rule and cohort files
remain pending.

The canonical v2 table has one row per biological sample, target taxon,
profiler, assembly arm, and implanted dose. It includes a zero-dose row for the
matched unmodified baseline. Consequently, fractions and targets are repeated
measurements—not independent samples.

## Required columns

| Column | Meaning |
|---|---|
| `schema_version` | Exactly `paired-dose-response-v2.1`. |
| `cohort` | Stable cohort key: `yachida`, `feng`, or `zeller`. |
| `study` | Published study identifier retained from the frozen manifest. |
| `sample_id` | Stable biological-sample identifier. |
| `condition` | Frozen phenotype label. |
| `analysis_population` | `community` or `independent`. |
| `target_label` | Stable spike-panel label, such as `Fnuc`. |
| `target_taxon` | Canonical target species name. |
| `assembly_arm` | `original`, `clean`, or `not_applicable`. |
| `profiler` | `kraken2_bracken` or `metaphlan4`. |
| `profile_id` | Unique upstream profile identifier. |
| `baseline_profile_id` | Profile identifier for the matched zero-dose library. |
| `spike_fraction_total` | Fraction of the final library occupied by the complete spike. |
| `spike_fraction_target` | Fraction of the final library implanted for this target. |
| `implanted_read_pairs_target` | Exact implanted paired-read count for this target. |
| `native_abundance` | Unmodified native profiler value for this row. |
| `native_unit` | `fraction_total_reads` or `relative_abundance_pct`. |
| `abundance_fraction` | Native value converted only by unit scaling to 0–1. |
| `detected_native_nonzero` | Mechanical indicator `abundance_fraction > 0`; not yet the frozen primary detection endpoint. |
| `source_profile` | Immutable native profile path or publication-relative identifier. |
| `source_design` | Immutable spike-design path or publication-relative identifier. |
| `include` | `1` for analysis or `0` for a prespecified exclusion. |
| `exclusion_reason` | Empty when included; otherwise a controlled reason. |

## Invariants

- The analytical identity (cohort, sample, population, target, assembly arm,
  profiler, and profile) is unique. A physical community profile may
  legitimately appear once for each implanted target.
- Every included positive-dose row has exactly one included matched zero-dose
  row with the same cohort, sample, population, target, profiler, and arm.
- At zero dose, both fractions and implanted pairs are zero, `profile_id` equals
  `baseline_profile_id`, and the design source may be `BASELINE`.
- Positive-dose rows have positive target fraction and implanted pairs.
- `spike_fraction_target <= spike_fraction_total < 1`.
- Bracken native values use `fraction_total_reads` and are copied unchanged to
  `abundance_fraction`. MetaPhlAn values use `relative_abundance_pct` and are
  divided by 100. This is unit conversion, not compositional closure.
- Missing measurements are not represented as zero. They require `include=0`
  and a non-empty reason, and the source ledger must retain them.
- No row is silently discarded. The validator writes a complete status report
  and fails if an invariant is violated.

Derived endpoints (`r`, `R`, and the read-proportional reference) do not belong
in the canonical evidence table. A later deterministic step will compute them
from this validated table and retain its checksum.
