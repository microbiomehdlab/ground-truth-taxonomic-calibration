# Real-data MetaPhlAn recovery checkpoint

Report: work/metaphlan_recovery_report_20261001T142003Z, job 3097832,
source ca549e9, completed 1m08s with exit 0. Local copy inspected in
../bracken_report/metaphlan_recovery_report_20261001T142003Z.
All local bundle checksums passed. Summary: 41,170 target endpoints,
1,170 context cells, non-fixture. Both arms use unchanged native abundances.

Primary recovery assessment: genome-equivalent reference using exact implanted
assembly lengths and baseline-derived sample-wide effective genome sizes.
Sensitivity: the historical read-fraction reference. This is an expectation
change, not learned calibration or modified profiler output.

## Findings

Across the 1,170 cell medians, median absolute relative error is 0.108 under
the primary reference versus 0.481 under the sensitivity reference. Cell median
ratios lie within 0.8–1.2 in 936 cells versus 154. These are descriptive summaries
across cells, NOT pooled sample estimates, biological CIs or inferential tests.
Paired sample-level error is lower in 30,935 endpoints, equal in 4,631 and higher
in 5,604; equal-error endpoints are native nondetections in this report.

At the weakest community total dose 0.0001 (0.01%, approximately 0.001% per
equal-weight target), 3,818 of 5,110 sample-target endpoints report zero native
target abundance. The corresponding primary median cell error is 1.0.
This is native non-detection, not direct identification of implanted read origins.
At community total doses >=0.005, all 90 cell median ratios per dose fall within
0.8–1.2; individual samples still vary. Independent-spike doses >=0.001 likewise
have all 90 cell median ratios within this descriptive band.

Residual differences remain: median across primary cell median ratios is
approximately 1.11 for Yachida individual spikes, 1.024 for Feng individual
spikes and 1.097 for Zeller individual spikes. Selected plots show stable
over-recovery in some contexts. Do not claim perfect calibration or uniform
performance merely because the corrected reference reduces apparent error.

Suggested interpretation: “Matching the expected signal to the profiler's
abundance scale substantially reduces apparent quantitative bias, while
low-abundance detection remains a major limitation.” Effective-size estimation
and coverage limits remain assumptions; no exact causal decomposition follows.

## What is frozen and what remains

Bracken primary: all-input estimated count recovery; native-fraction sensitivity.
MetaPhlAn primary: genome-equivalent expected-reference recovery;
read-fraction sensitivity. Both retain zeros and negative responses.
DA input is native observed abundance, not corrected recovery ratios.
No upstream rerun is required. No DA models were executed by these reports.
DA protocol implementation, null validation and small end-to-end pilots remain.
