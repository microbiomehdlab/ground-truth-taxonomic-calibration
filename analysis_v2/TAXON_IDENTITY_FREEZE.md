# Taxon identity freeze for definitive downstream analysis

**Status:** source identity fixed; real-profile semantic audit remains a cohort-run gate.
**Lock file:** `analysis_v2/taxon_identity_freeze.sha256`.

From the repository root, verify the exact bytes with:

```bash
sha256sum -c analysis_v2/taxon_identity_freeze.sha256
```

The lock covers `spikes/spike_panel.tsv` (the 10 implanted labels, species names, assemblies, and weights), `examples/spike_taxon_aliases.csv` (one exact target feature per implanted species and profiler), and `analysis_v2/feature_equivalence_aliases.tsv` (four additional source-to-canonical feature mappings for the response atlas). The first two files were inputs to the verified upstream evidence package built from source commit `9df9ae7dce670ec83fbbfbb709c10a0bb51c4947`. The third file is a downstream feature-policy input and is locked here separately.

A local structural check found 10 unique panel labels and species, 20 distinct `(canonical, profiler)` target mappings, and four distinct `(profiler, source_feature)` feature-equivalence mappings. Every equivalence destination is one of the 10 panel species. This confirms completeness and unambiguous keys in these files; it does not establish that every declared feature is observed in every native profile.

The target aliases and feature-equivalence mappings serve different questions. For example, the direct B. fragilis target uses the exact `Bacteroides fragilis` row in each native profiler, while the response atlas can additionally group Bracken's `Bacteroides fragilis_A` under the canonical B. fragilis species group. Do not silently sum the split label into the direct target measurement. Both profilers use the exact species-level `Fusobacterium nucleatum` feature for the implanted `F. nucleatum subsp. nucleatum` assembly. Keep `Fusobacterium nucleatum_E` distinct; this freeze has no rule collapsing it into Fnuc. No fuzzy or prefix matching is authorized.

Before any definitive cohort model is accepted, verify this lock, build the canonical input, and run `analysis_v2/scripts/audit_profiler_semantics.py` against the included native files. Archive the audit table, input hashes, and `SUCCESS` with the cohort run. Inspect target detection and non-detection examples from every cohort and profiler. If an expected target feature is absent from a profile, the exact-match reader records zero; absence must be interpreted using the native profile and detection policy, not repaired by changing aliases after viewing results. A proposed alias change requires a new identity-lock revision, documented scientific justification, and regeneration of affected downstream outputs.

The frozen profiler settings are described in `analysis_v2/PROFILER_SEMANTICS.md`. All three actual strict-production environment files specify Bracken read length 100 and abundance threshold 10. A user-run, read-only scan on Lobo found those values in every retained `profiling_parameters.tsv` encountered: 3,415 for Yachida, 3,045 for Feng, and 3,054 for Zeller. That broad scan included 7, 13, and 6 more parameter files than the sealed production profile counts, so it is evidence for the settings, not a replacement for the exact profile-topology seals. The generic profiling script has a fallback of 150, but it was not the recorded setting in these files. Exact image/database identities remain in the upstream evidence package.
