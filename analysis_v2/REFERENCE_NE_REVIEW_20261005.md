# Verified non-estimability review — 5 October 2026

Evidence: cluster job 3108950 completed in 3 seconds, exit 0:0, after the audit
schema fix a88881a. Downloaded bundle:
`../reference_ne_audit_20261005T215250Z/`. SHA256SUMS verified locally. No model
refits or profiling runs. Original sealed runs and review figure bundle unchanged.

All 439 target/context rows pass arithmetic checks against saved model status.
The audit includes 11820 person-target observations; people recur across draws.

## MetaPhlAn: analytical zeros, exclusively at the lowest tested dose

All 291 constant-input cases contain exactly one unique native value: zero.
Every artificial control and every exposed person reports zero for that target.
All occur at nominal total community fraction 0.0001 (approximately 0.001% per
target). By n/group: 179 cases at n=10, 112 at n=20. Cohorts: Yachida 114, Feng
63, Zeller 114. There are no such cases at the two higher pilot doses.

| Target | All-zero target/context cases |
|---|---:|
| Pint | 109 |
| Pana | 68 |
| Porp | 63 |
| Fnuc | 35 |
| Psto | 14 |
| Pmic | 1 |
| Dpne | 1 |
| Bfrag, Csym, Hhat | 0 |

This verifies analytical non-reporting despite the known addition, not just a
valid nonsignificant coefficient. Zero native abundance does not show absence of
organisms, absence of alignments or which internal marker/threshold rule caused
the zero. This audit alone is not a marker-level mechanistic explanation.

## Bracken: complete reported separation, not absent signal

All 148 perfect-fixed-fit cases have exactly two unique native values. Every
artificial control reports zero; every exposed person reports the same strictly
positive target abundance within that context. All group differences are positive.
This creates zero within-group residual variance and is prescreened as a perfect
fixed fit by the existing inferential contract. Do not assign p=0 or count it as
ordinary nonsignificance or analytical non-detection.

| Target | Perfect-fixed-fit target/context cases |
|---|---:|
| Porp | 94 |
| Pint | 20 |
| Dpne | 19 |
| Fnuc | 15 |

118 cases occur at total fraction 0.0001; 30 at 0.001; none at 0.01. Yachida
64, Feng 42, Zeller 42; n=10 has 108, n=20 has 40. Identical printed abundances
could reflect output precision, but this audit does not establish that cause.

## Figure/caption and manuscript implications

- The existing separate NE pages are retained, but their caption must now say
  “Bracken perfect fixed fits: zero controls/identical positive exposed values”
  versus “MetaPhlAn all-zero native inputs.” Never call both missing signals.
- In standard-error/effect summaries, m=0 means no jointly estimable comparisons,
  not zero uncertainty. Lines must not bridge missing effect estimates across an
  unobserved dose; inspect final panels with zero matched-estimable cells.
- Every significance frequency retains all allocations in its denominator and
  must be shown alongside NE composition. For Bracken, a lower significance rate
  can coexist with complete reported separation excluded by the model guard.
- Cautious Results wording: “At the lowest tested addition, 291 MetaPhlAn
  target/allocation comparisons had zero reported abundance in every exposed and
  unexposed sample. In contrast, 148 Bracken comparisons were non-estimable
  because all unexposed samples reported zero and all exposed samples reported
  an identical positive value.” Denominators, doses and repeated-draw scope must
  accompany this statement. These counts are not clinical disease-association
  failures or independent patient counts.

Next: reference-construction sensitivity and broader saved-allocation/exposure
review, then a compact paper figure and private manuscript integration. The
zero/nonzero input question is now answered; no repeat NE audit is needed.
Marker/read-origin inspection is optional targeted mechanistic follow-up, not
an automatic new profiling campaign or prerequisite for describing these zeros.
