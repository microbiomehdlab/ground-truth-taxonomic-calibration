# MaAsLin2 implementation checkpoint

The user agreed to a fixed-pseudocount primary transformation on 2026-10-01,
with sensitivity checks and baseline-null validation, not a claim of unique
correctness. DA_PROTOCOL_AMENDMENT_20261001.md is the current design amendment.
The original DA plan/statistical protocol now carry supersession notices.

Implemented locally: lib/maaslin_contract.R and tests/test_maaslin_contract.R.
Pure contracts enforce native fractions, finite data, complete frozen family,
exact sample intersection, deterministic metadata order, and one original plus
one spiked observation per person in paired contexts. The shifted fixed log
transform preserves group coefficients relative to log2(a+1e-8), with zero
input transformed exactly to zero. Missing metadata is an error, not silently
removed. Matrix preparation does not select the family or parse raw profiles.

full_family_bh requires one accounted-for result for every frozen member;
declared non-estimable raw p is NA and its bookkeeping p is 1. It refuses
missing/duplicate feature results, nonfinite valid p and unjustified numeric p
on a non-estimable result. This contract is a proposed wrapper-multiplicity
rule, not validated handling of every possible MaAsLin2 fit failure.

Local R fixtures pass for transform identity, LM group-coefficient invariance,
metadata order, pairing defects, invalid inputs and hand-calculated full-family
BH. MaAsLin2 is absent from local R; no package-backend or mixed-model claims
are made. The repository image specification pins Maaslin2 1.18.0; verify the
actual container at the next user-run gate.

Next implementation packages, in order:

Development context wrapper now exists in lib/maaslin_context.R. Local mock
fixtures passed; actual MaAsLin2 backend remains UNVERIFIED. Run
test_maaslin_context.R with REQUIRE_MAASLIN_BACKEND=1 inside the pinned image
to exercise actual unpaired and paired synthetic comparisons. Saved-model
schema and warnings may reveal pinned-package compatibility issues; no study
model or definitive SUCCESS is authorized by local mock tests. The wrapper
emits DEVELOPMENT_NOT_DEFINITIVE on successful contexts.

1. Inspect the existing native-species input builder for completeness, finite
   values, profile provenance and metadata; freeze eligibility, prevalence
   families and dose/sample feasibility. Preserve all ten targets in ledgers.
2. Implement MaAsLin2 context wrapper with explicit models/references, package
   version checks and native-versus-wrapper q fields. Check exact 1.18.0 model
   output, internal filtering, covariance/df and warning behavior in container.
3. Add direct LM/mixed-model fixtures and conservative constant/near-zero
   variance handling. Establish singular/convergence diagnostic policy before
   real outcomes. No automatic unpaired fallback.
4. Add deterministic nested allocations, 252-subset enumeration and null arm.
5. Authorized schema/runtime pilot, then full baseline-null review gates;
   definitive spike fits remain on hold.

No cluster model submission, package installation, manuscript edit, commit or
push was performed by this checkpoint. Existing historic models/policies were
not overwritten. Recovery results are saved separately in
METAPHLAN_RECOVERY_INTERPRETATION_20261001.md and RECOVERY_DECISIONS_20261001.md.
