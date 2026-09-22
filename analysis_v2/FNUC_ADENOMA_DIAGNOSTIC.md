# F. nucleatum adenoma diagnostic evidence map

This figure is a compact follow-up to the matched-species atlas, not a claim that one cause of nonsignificance has been identified. It displays, for each cohort and profiler: adjusted Adenoma-vs-Control and CRC-vs-Control effects; full-cohort unspiked adenoma detection and abundance spread; paired zero-to-positive rescue and quantitative recovery of a 0.001% per-member **community-mixture** spike; quantitative recovery of a separate 0.01% **single-species** spike; and the strongest *Fusobacterium* bystander response slope for an independent *F. nucleatum* perturbation. Each metric retains its own denominator and scale. The bystander slope is not a percentage of reads misclassified.

The community mixture totals about 0.01% across ten species. Its lowest per-member dose is 0.001%, an order of magnitude below the lowest independent dose shown in the previous matched-species plot. These are different perturbation designs; do not connect them as one dose-response curve. A strong spike response weakens that measured failure mode at that dose, but cannot establish biological absence at lower native abundance. Exact-looking recovery ratios should be checked against unrounded source endpoints and profile provenance before being presented as precision.

Run from the repository root on lobo after setting `ROOT`:

```bash
MATCHED="$ROOT/matched_crc_adenoma_species_20260922T113325Z"
OUTDIR="$ROOT/fnuc_adenoma_diagnostic_$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/plot_fnuc_adenoma_diagnostic.py \
  --matched-dir "$MATCHED" \
  --endpoints "$ROOT/endpoints/paired_endpoints.tsv" \
  --operator "$ROOT/partial_snapshot_response/response_analysis/response_operator.tsv" \
  --outdir "$OUTDIR"
```

Use the **newer four-species matched output directory** if available; the three-species output also contains *F. nucleatum* and is sufficient for this focused plot. Open `fnuc_adenoma_diagnostic.svg` after `SUCCESS` exists. The companion TSV contains every displayed number.

The most defensible presentation claim is narrow: at the tested 0.001% member dose, low-dose *F. nucleatum* rescue differs sharply by profiler; a failure shared across clinical backgrounds cannot alone explain an adenoma-specific non-call. Review the source design and QC before attributing any difference to biology.
