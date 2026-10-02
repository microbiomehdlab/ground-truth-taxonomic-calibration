# DA3 100-allocation exploratory launch

User-approved next step after the completed canary, global-null and partial-null
reviews. DA1 remains unresolved and is NOT run by this submission.

Reuse the checksummed `draft_da3_plan` inside
`work/preproduction_checks_20261002T192903Z`. Do not regenerate the plan or
reuse earlier canary outputs under the changed resume contract.

`submit_da3_exploratory.sh` verifies the complete plan before any submission,
snapshots committed code, fingerprints the analysis image and submits six
800-task arrays. Each batch holds 25 contexts; each array allows five simultaneous
tasks, giving a combined maximum of 30. Arrays have offsets 0 through 4000.
Collector runs after any termination and reports all missing/corrupt batches.

The exploratory run has 120,000 contexts. Based on the compact canary,
approximately 2.6 hours at concurrency 30 was estimated, before scheduler,
image checksum and storage overhead. Image hashing per worker adds read traffic;
this is a reproducibility check, not statistical computation. Do not promise
an exact completion time. Each task requests one CPU and 8 GB; 30 concurrent
tasks may reserve 240 GB across the allocated nodes.

Existing native profiles are read-only. Workers check code/image/source identity,
lock outputs and atomically publish finished batches. Do not resubmit the entire
run into a fresh root merely because a few batches fail. Review the failure
ledger, then retry those batches with the same snapshot, image, plan and results
directory. Never merge outputs from different policies or snapshots.

All results are exploratory finite-cohort simulation evidence: FWER versus BH
remain separate, nominal read doses are not cell abundances, artificial groups
are not real disease/control groups, and repeated allocations do not increase
biological sample size. Collector completeness does not automatically authorize
publication claims or select the DA1 primary statistical method.
