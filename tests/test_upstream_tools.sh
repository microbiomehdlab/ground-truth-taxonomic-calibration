#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP_ROOT="$(mktemp -d)"
trap 'rm -rf "$TMP_ROOT"' EXIT

calc="$ROOT/spikes/scripts/spikein/calc_spike_reads.py"
out="$(python3 "$calc" --R 999900 --f 0.0001)"
grep -q $'^N\t100$' <<<"$out"
grep -q $'^f_hat\t0.00010000$' <<<"$out"

seed="$ROOT/spikes/scripts/spikein/stable_seed.py"
[[ "$(python3 "$seed" --base 13 --namespace spike-independent-v1 DRR127476 Fnuc 0.0001)" == 873630871 ]]
! grep -Eq 'SEED_BASE[[:space:]]*\+[[:space:]]*task|task[[:space:]]*\*[[:space:]]*1000' \
  "$ROOT/spikes/scripts/spikein/spike_one_taxon_array.sbatch" \
  "$ROOT/spikes/scripts/spikein/spike_community_array.sbatch"
grep -q -- '-rs "$ART_SEED"' "$ROOT/spikes/scripts/spikein/make_spike_pool.sbatch"

mkdir -p \
  "$TMP_ROOT/results/S1_Bfrag_f0p0001/metaphlan4" \
  "$TMP_ROOT/results/S2_Bfrag_f0p0001/metaphlan4"

printf 'sample,Bacteroides fragilis\nS1_Bfrag_f0p0001,0.01\n' \
  > "$TMP_ROOT/results/S1_Bfrag_f0p0001/metaphlan4/S1_Bfrag_f0p0001.open_world.csv"
printf 'sample,Bacteroides fragilis\nS2_Bfrag_f0p0001,0.02\n' \
  > "$TMP_ROOT/results/S2_Bfrag_f0p0001/metaphlan4/S2_Bfrag_f0p0001.open_world.csv"

python3 "$ROOT/spikes/scripts/spikein/organize_profile_tables_by_spike.py" \
  --inroot "$TMP_ROOT/results" \
  --outroot "$TMP_ROOT/profile_results" \
  --kind open_world

organized="$TMP_ROOT/profile_results/Bfrag/metaphlan4/spike_f0p0001.open_world.csv"
[[ -s "$organized" ]]
[[ "$(wc -l < "$organized")" -eq 3 ]]

python3 "$ROOT/spikes/scripts/spikein/merge_profile_tables.py" \
  --root "$TMP_ROOT/profile_results" \
  --outdir "$TMP_ROOT/merged" \
  --kind open_world

merged="$TMP_ROOT/merged/metaphlan4.open_world.merged.csv"
[[ -s "$merged" ]]
[[ "$(wc -l < "$merged")" -eq 3 ]]

echo "[PASS] Upstream utility smoke tests passed."
