# Rolling CRC production

`submit_crc_rolling.sh` separates network staging from compute so cohort runs can
use many compute slots without opening the same number of ENA connections.

- `--lanes 42` bounds active sample lanes and compute concurrency.
- `--download-concurrent 8` bounds downloads within a generation.
- Every lane stages, verifies, computes, verifies retained outputs, cleans its
  disposable inputs, and only then admits the next sample assigned to that lane.
- Download generations cannot overlap. Compute generations may overlap, but the
  per-lane `aftercorr` dependency ensures at most one compute task per lane.
- A failed stage blocks its matching compute/lane; it does not invalidate other
  lanes. Missing samples are recovered from persistent `.verified` markers.

Example:

```bash
export PROJECT="$PWD"
export CRC_ENV="$PWD/config/zeller.strict-production.env"
bash datasets/submit_crc_rolling.sh \
  --manifest datasets/zellerg/manifests/production_manifest.tsv \
  --lanes 42 \
  --download-concurrent 8 \
  --exclude compute-1,compute-10 \
  --delete-verified-inputs
```

The submitter writes all Slurm job IDs and generation offsets to the reported
`jobs.tsv`. Production completion must still be established from the persistent
per-sample receipts and verified markers, followed by the cohort production
audit; Slurm completion alone is not a production seal.
