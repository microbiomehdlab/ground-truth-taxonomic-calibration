#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${DEV_INPUT_ROOT:?}"; : "${ANALYSIS_SIF:?}"; : "${OUTDIR:?}"
test ! -e "$OUTDIR" || { echo "[ERROR] OUTDIR exists: $OUTDIR" >&2; exit 1; }
read -r paired_n disease_n < <(python3 analysis_v2/scripts/shard_biomarker_model_input.py --profile-manifest "$DEV_INPUT_ROOT/canonical_input.tsv" --counts-only)
mkdir -p "$OUTDIR/slurm"; export PROJECT="$ROOT" DEV_INPUT_ROOT ANALYSIS_SIF OUTDIR
S=${SBATCH_BIN:-sbatch}; worker=analysis_v2/run_legacy_biomarker_mapreduce_stage.sbatch
one(){ local stage=$1 dep=$2 cpu=$3 mem=$4 time=$5; shift 5; local x=(--parsable --job-name="mr_$stage" --cpus-per-task="$cpu" --mem="$mem" --time="$time" --output="$OUTDIR/slurm/$stage.%A_%a.out" --error="$OUTDIR/slurm/$stage.%A_%a.err" --export="ALL,PIPELINE_STAGE=$stage"); test -z "$dep" || x+=(--dependency="afterok:$dep"); "$S" "${x[@]}" "$@" "$worker"; }
setup=$(one setup '' 2 32G 08:00:00)
paired=$(one paired_map "$setup" 2 24G 1-00:00:00 --array="1-${paired_n}%${MAP_CONCURRENCY:-12}")
disease=$(one disease_map "$setup" 2 32G 1-00:00:00 --array="1-${disease_n}%${DISEASE_CONCURRENCY:-4}")
paired_reduce=$(one paired_reduce "$paired" 2 24G 06:00:00)
disease_reduce=$(one disease_reduce "$disease" 2 24G 06:00:00)
art_report=$(one artificial_report "$paired_reduce" 2 24G 12:00:00)
dis_report=$(one disease_report "$disease_reduce" 2 24G 12:00:00)
linkage=$(one linkage "$paired_reduce" 2 24G 12:00:00)
seal=$(one seal "$art_report:$dis_report:$linkage" 1 4G 01:00:00)
printf '[SUBMITTED] map-reduce run\nOUTDIR=%s\npaired_shards=%s\ndisease_shards=%s\nsetup=%s\npaired_array=%s\ndisease_array=%s\npaired_reduce=%s\ndisease_reduce=%s\nartificial_report=%s\ndisease_report=%s\nlinkage=%s\nseal=%s\n' "$OUTDIR" "$paired_n" "$disease_n" "$setup" "$paired" "$disease" "$paired_reduce" "$disease_reduce" "$art_report" "$dis_report" "$linkage" "$seal"
