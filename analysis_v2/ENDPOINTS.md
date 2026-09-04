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
