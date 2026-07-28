#!/usr/bin/env bash
set -euo pipefail

LOCAL_PROJECT="${LOCAL_PROJECT:-$PWD}"
DATA_ROOT="${DATA_ROOT:-$LOCAL_PROJECT}"
: "${REMOTE:?Set REMOTE to user@cluster}"
: "${REMOTE_PROJECT:?Set REMOTE_PROJECT to the destination directory}"

LOCAL_PROJECT="$(cd "$LOCAL_PROJECT" && pwd -P)"
DATA_ROOT="$(cd "$DATA_ROOT" && pwd -P)"
input_files=(
  samples_spiked_all_independent.tsv samples_spiked_all_community.tsv
  spike_panel.tsv spike_taxon_aliases.csv metadata_w_study.tsv
  metaphlan4_merged_unspecified.csv kraken2_bracken_merged_unspecified.csv
)
code_files=(
  README.md REQUIRED_INPUTS.md stage_required_inputs.sh
  transfer_analysis_to_cluster.sh
  run_publication_original_unpaired_q010.sh
  run_original_unpaired_q010_cluster.sbatch
  cluster_preflight_original_unpaired.sh rscript_in_container.sh
  build_crc_spike_maaslin2_container.sh
)
for file in "${input_files[@]}"; do
  [[ -s "$DATA_ROOT/$file" ]] || {
    echo "[ERROR] Missing or empty input file: $DATA_ROOT/$file" >&2
    exit 1
  }
done
for file in "${code_files[@]}"; do
  [[ -s "$LOCAL_PROJECT/$file" ]] || {
    echo "[ERROR] Missing or empty code file: $LOCAL_PROJECT/$file" >&2
    exit 1
  }
done
[[ -d "$LOCAL_PROJECT/scripts" ]] || {
  echo "[ERROR] Missing code directory: $LOCAL_PROJECT/scripts" >&2
  exit 1
}
PROFILE_ROOT="$LOCAL_PROJECT/profile_results"
if [[ ! -d "$PROFILE_ROOT" ]]; then
  PROFILE_ROOT="$DATA_ROOT/profile_results"
fi
[[ -d "$PROFILE_ROOT" ]] || {
  echo "[ERROR] Missing profile_results in LOCAL_PROJECT and DATA_ROOT." >&2
  exit 1
}
if ! find "$PROFILE_ROOT" -type f -print -quit | grep -q .; then
  echo "[ERROR] profile_results is empty: $PROFILE_ROOT" >&2
  exit 1
fi
R_ROOT="$LOCAL_PROJECT/R"
if [[ ! -d "$R_ROOT" ]]; then
  R_ROOT="$(cd "$LOCAL_PROJECT/.." && pwd -P)/R"
fi
[[ -d "$R_ROOT" ]] || {
  echo "[ERROR] Missing R helper directory in LOCAL_PROJECT and its parent." >&2
  exit 1
}
required_r_helpers=(
  common_utils.R io_utils.R maaslin_utils.R spike_design_utils.R
  spike_metrics_utils.R spike_plot_utils.R
  spike_plot_utils_independent_patch.R classifier_utils.R
)
for file in "${required_r_helpers[@]}"; do
  [[ -s "$R_ROOT/$file" ]] || {
    echo "[ERROR] Missing or empty R helper: $R_ROOT/$file" >&2
    exit 1
  }
done
if [[ -d "$LOCAL_PROJECT/containers" ]]; then
  include_containers=true
else
  include_containers=false
fi

inventory="$LOCAL_PROJECT/TRANSFER_SHA256.txt"
{
  for file in "${input_files[@]}"; do
    sha256sum "$DATA_ROOT/$file"
  done
  for file in "${code_files[@]}"; do
    sha256sum "$LOCAL_PROJECT/$file"
  done
  find "$LOCAL_PROJECT/scripts" "$R_ROOT" "$PROFILE_ROOT" \
    -type f -print0 | sort -z | xargs -0 sha256sum
  if [[ "$include_containers" == "true" ]]; then
    find "$LOCAL_PROJECT/containers" -type f -print0 | sort -z | xargs -0 sha256sum
  fi
} > "$inventory"

ssh "$REMOTE" mkdir -p "$REMOTE_PROJECT"
rsync -avP --checksum \
  "${input_files[@]/#/$DATA_ROOT/}" \
  "$REMOTE:$REMOTE_PROJECT/"
rsync -avP --checksum \
  "${code_files[@]/#/$LOCAL_PROJECT/}" \
  "$inventory" \
  "$REMOTE:$REMOTE_PROJECT/"
rsync -avP --checksum "$LOCAL_PROJECT/scripts/" "$REMOTE:$REMOTE_PROJECT/scripts/"
rsync -avP --checksum "$R_ROOT/" "$REMOTE:$REMOTE_PROJECT/R/"
rsync -avP --checksum "$PROFILE_ROOT/" "$REMOTE:$REMOTE_PROJECT/profile_results/"
if [[ "$include_containers" == "true" ]]; then
  rsync -avP --checksum "$LOCAL_PROJECT/containers/" "$REMOTE:$REMOTE_PROJECT/containers/"
fi

echo "[PASS] Transfer completed: $REMOTE:$REMOTE_PROJECT"
echo "[INFO] Code root:       $LOCAL_PROJECT"
echo "[INFO] Input data root: $DATA_ROOT"
echo "[INFO] Profile root:    $PROFILE_ROOT"
echo "[INFO] R helper root:   $R_ROOT"
echo "[INFO] Local inventory: $inventory"
