# Deterministic paired endpoints

**Status:** implemented and fixture-tested; final cohort derivation pending.

For a matched sample, target, profiler, and assembly arm, define:

- `o`: baseline abundance, expressed as a fraction after unit scaling only;
- `a`: abundance in the perturbed profile, on the same scale;
- `F`: total fraction of the final library occupied by all implanted reads;
- `f`: target-specific implanted fraction of the final library.

The read-proportional reference is

`e_RP = (1 - F)o + f`.

The baseline-adjusted recovered spike signal and response ratio are

`r = a - (1 - F)o`

and

`R = r / f`.

For an individual spike, `F = f`. For a community spike, `F` is the complete
community fraction while `f` is the fraction contributed by the specific
target. This distinction prevents over-diluting or over-crediting individual
community members.

`R = 1` represents agreement with the read-proportional reference. Values below
or above one represent lower or higher response on that profiler's native
abundance scale. Negative `r` and `R` values are retained because they are
possible after baseline adjustment and must not be silently truncated.

No pseudocount, log transform, compositional closure, or detected-only filter is
applied at this stage. These are deterministic derived data, not fitted results.
The derivation writes the validated-input report, endpoint table, explicit
exclusion ledger, checksums, and a success marker.

## Profiler-specific reference scale (17 September 2026)

Bracken keeps the read-proportional reference `e_RP = (1 - F)o + f` exactly.

MetaPhlAn reports a genome-equivalent-like composition, so the implanted read
fraction is rescaled and the composition renormalised. For implanted target `t`
and non-implanted feature `j`, with audited sample-wide `G_eff,i`:

```
q_it   = f_it * G_eff,i / G_t
Q_i    = sum_k(f_ik * G_eff,i / G_k)          over every implanted member k
e_MP,it = ((1 - F_i) o_it + q_it) / ((1 - F_i) + Q_i)
e_MP,ij = ((1 - F_i) o_ij)        / ((1 - F_i) + Q_i)
```

`derive_paired_endpoints.py` implements the implanted-target case additively:
the read-proportional fields keep their names and meaning, and the profiler
scale is carried by `reference_type` plus the `*_profiler_scale` fields.
`scripts/build_perturbation_response_input.py` implements the non-implanted case
`e_MP,ij`, computing `Q_i` once per perturbed profile over every implanted
member. Each row also carries `signed_bounded_error_profiler_scale` and
`absolute_bounded_error_profiler_scale`, computed from the profiler-scale
expected value and validated within [-1, 1] and [0, 1].

Each row also carries `implanted_signal_profiler_scale`, so that
`expected_abundance_profiler_scale = retained_baseline_profiler_scale +
implanted_signal_profiler_scale` holds exactly and recovery ratios divide two
quantities on the same scale. Consumers select the scale with a required
`--reference-scale {profiler_scale|read_proportional}`;
`expected_abundance_fraction` remains the read-proportional reference for every
profiler and is aliased as `read_proportional_reference`.

Runners select the scale through `analysis_v2/lib/metaphlan_reference.sh`.
There is no default: MetaPhlAn rows without both `TARGET_GENOME_SIZES` and
`EFFECTIVE_GENOME_SIZES` fail closed. The primary genome-equivalent result and
the labelled read-proportional sensitivity are written to separate directories.
