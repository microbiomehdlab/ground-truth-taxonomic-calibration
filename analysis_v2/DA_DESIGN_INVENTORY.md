# DA2/DA3 design inventory checkpoint

Local implementation: scripts/build_da_design_inventory.py and da_allocations.py.
No study models or cluster runs. Real-data inventory review is still pending.

Shared cluster launcher: run_three_cohort_da_inventory.sbatch. Export PROJECT,
ANALYSIS_SIF, fresh DA_INVENTORY_ROOT and PREFLIGHT_STAMP=20260929T155304Z.
The job fixes 1,000 full-pool allocations; it only writes their memberships,
not 1,000 model fits. Baseline families include clinical controls for DA2, while
DA3 allocations are only adenoma and CRC. Each ten-sample subset has 252 splits.
All cohort configurations are loaded in separate subshells to prevent cross-
cohort environment leakage. Inputs must be visible through /mnt/nfs binding.
On success verify root SHA256SUMS and inspect eligibility and family sizes;
the root marker explicitly says definitive model readiness is not established.

Inputs: one definitive canonical input, matching current cohort state directory,
expected sample count and fresh output directory. Verify checksummed v2 seal,
PASS flow, sample membership and retained receipt hashes for every native profile.
No recursive scan of obsolete/quarantined profiles. All canonical inputs must
be included; unexpected exclusions require explicit policy review.

Outputs: physical profile inventory with nominal doses/hashes; baseline-frozen
DA2/DA3 clinical-background species families (10% baseline prevalence union ten
exact target aliases); eligibility for n=5/10/15/20; allocations; source hashes;
SUCCESS and SHA256SUMS. Family construction reads original raw species abundances,
never post-spike response. DA1 families remain PENDING_COVARIATE_REVIEW because
age/sex complete-case eligibility is not yet defined by this tool.

For DA3 adenoma and CRC pools, assignments are independent by stable SHA256 seed
namespace (cohort/background/pool/version), nested across n, and reusable across
all profilers/doses/null arms. Full-pool default 1,000 independent random draws;
each ten-sample independent subset enumerates 252 labelled five/five allocations.
This ledger defines biological assignments, not extra biological replication.
Python version and runtime must be recorded when freezing definitive outputs.

The inventory currently validates the retained total-dose grid per target series;
it is not the final sample-specific achieved-target dose mapping or proof of
individual/community matched-dose support. These ledgers and model-context
cost/selection manifests remain next implementation tasks. No family ID/config
freeze or definitive-run readiness follows automatically from inventory SUCCESS.

Local tests: test_da_allocations.py (3 tests), test_da_design_inventory.py (1
integrated synthetic inventory/corruption test). Allocation tests cover sorted-
input invariance, nested/disjoint groups, feasible n, independent namespaces,
252 unique assignments and duplicate biological IDs. Synthetic inventory covers
all retained dose series, receipt checks, families, count expectations and source
corruption refusal before output. It does not validate real cohort eligibility.
