# Scaling sensitivity: verified results, 5 October 2026

Downloaded source: ../reference_scaling_20261005T222723Z_REPORT/. All report
checksums verified locally. 10800 target/scenario rows, 3600 unique target/context
pairs, exactly three scenarios each. Independently checked: every observed field
is identical across scenarios; every primary-scenario field matches the original
downloaded pilot report. No observed fits were replaced.

At the lowest nominal dose (0.001% per target), n=20/group, pooled across three
cohorts × ten targets × 20 saved allocations (600 comparisons per scenario):

| Conversion multiplier | Reference positive, primary | Reference positive, HC3 | Observed positive, primary |
|---|---:|---:|---:|
| 0.5 | 475/600 | 470/600 | 1/600 |
| 1.0 | 495/600 | 489/600 | 1/600 |
| 2.0 | 512/600 | 510/600 | 1/600 |

The large low-dose observed/reference contrast persists under this prespecified
half/double global conversion stress. Each cohort retains the same qualitative
contrast: half-scaling reference counts are Yachida 148/200, Feng 172/200 and
Zeller 155/200, versus observed 0/200, 1/200 and 0/200 respectively.
At n=10 the reference counts are 427/600, 441/600, 452/600; observed 0/600.

At n=20 and the lowest dose, Dpne, Fnuc, Pana, Pint, Pmic, Porp and Psto are
reference-positive in 60/60 target allocations under all three scenarios and
HC3. Other targets are meaningfully sensitive:

| Target | Reference positive, 0.5 | 1.0 | 2.0 |
|---|---:|---:|---:|
| Bfrag | 8/60 | 8/60 | 12/60 |
| Csym | 31/60 | 42/60 | 47/60 |
| Hhat | 16/60 | 25/60 | 33/60 |

Do not claim that every species' expected detectability is invariant. At the
higher doses the reference and observed rates are closer and changing the scale
can move the reference rate above or below the observed rate. At n=20, 0.01%
per-target, reference counts are 535/558/565 versus observed 505/600; at 0.1%
they are 580/590/592 versus observed 590/600.

This tests one uniform scaling assumption only. Multipliers are broad stress
scenarios, not empirical confidence bounds. It does not establish a biological
truth scale, marker-level failure mechanism, statistical population power or
causal explanation of natural adenoma non-association. Repeated allocations and
targets are not independent biological replicates. Verified all-zero observed
inputs remain unchanged under every reference multiplier.

## Next scope

The sensitivity supports extending the same pilot grid and uniform arm to all
100 saved allocations, to examine conditional sample-selection stability. This
is the next recommended implementation, not an automatic production authorization.
Keep all targets/cohorts, doses and primary settings fixed; do not select favorable
taxa or multipliers. Heterogeneous exposure and other sample sizes remain later
explicit extensions. Manuscript integration should retain the species-sensitive
results alongside the robust low-dose contrast.
