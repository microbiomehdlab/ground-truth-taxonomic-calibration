#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
: "${SIF:?Set SIF to the analysis container image containing Mash}"

PROJECT="$(cd "$PROJECT" && pwd -P)"
binds=("$PROJECT:$PROJECT")

# Mash receives the sketch and query FASTAs as ordinary path arguments. Bind
# their parent directories explicitly so this works even when the cluster does
# not mount all of /mnt inside containers by default.
for argument in "$@"; do
  if [[ "$argument" = /* && -e "$argument" ]]; then
    parent="$(dirname "$argument")"
    candidate="$parent:$parent"
    already=false
    for existing in "${binds[@]}"; do
      [[ "$existing" == "$candidate" ]] && already=true
    done
    [[ "$already" == "true" ]] || binds+=("$candidate")
  fi
done

bind_csv="$(IFS=,; echo "${binds[*]}")"
exec apptainer exec --cleanenv \
  --bind "$bind_csv" \
  "$SIF" \
  mash "$@"
