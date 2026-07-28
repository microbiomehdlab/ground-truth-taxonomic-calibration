#!/bin/bash
# Download ZellerG_2014 FASTQ files from EBI, organized by Target Condition
#
# Usage:
#   ./download_zellerg2014.sh <output_dir>
#
# Arguments:
#   <output_dir>     Root directory where condition subdirectories will be created
#
# Expected in the same directory as this script:
#   - Public_study__ZellerG_2014.tsv   (metadata with 'Name' and 'Target Condition' columns)
#   - zellerg_wgets.sh                 (wget script with '# Sample: <name>' blocks)
#
# Output structure:
#   <output_dir>/
#       zellerg/
#          ├── Control/
#          │   ├── CCIS00146684ST-4-0/
#          │   │   ├── ERR478960_1.fastq.gz
#          │   │   └── ...
#          ├── colorectal_carcinoma/
#          │   ├── CCIS02379307ST-4-0/
#          │   └── ...
#          └── Adenoma/
#              └── ...

set -e

if [ -z "$1" ]; then
  echo "Usage: $0 <output_dir>"
  exit 1
fi

OUTPUT_DIR="$1/zellerg"

# Locate metadata and wget script relative to this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA_TSV="$SCRIPT_DIR/Public_study__ZellerG_2014.tsv"
WGET_SCRIPT="$SCRIPT_DIR/zellerg_wgets.sh"


# Reads TSV header to find column indices for 'Name' and 'Target Condition'
declare -A SAMPLE_CONDITION
while IFS=$'\t' read -r name condition; do
  # Sanitize condition for use as directory name (spaces -> underscores)
  safe_condition="${condition// /_}"
  SAMPLE_CONDITION["$name"]="$safe_condition"
done < <(awk -F'\t' '
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i == "Name") name_col = i
      if ($i == "Target Condition") cond_col = i
    }
    next
  }
  { print $name_col "\t" $cond_col }
' "$METADATA_TSV")




current_sample=""
declare -a current_urls=()

flush_sample() {
  if [ -z "$current_sample" ]; then
    return
  fi

  condition="${SAMPLE_CONDITION[$current_sample]}"

  sample_dir="$OUTPUT_DIR/$condition/$current_sample"
  mkdir -p "$sample_dir"

  echo "[$condition] $current_sample (${#current_urls[@]} files)"

  for url in "${current_urls[@]}"; do
    filename=$(basename "$url")
    if [ -f "$sample_dir/$filename" ]; then
      echo "Skipping $filename (already exists)"
    else
      echo "Downloading $filename"
      wget -nc -q -P "$sample_dir" "$url"
    fi
  done

  current_sample=""
  current_urls=()
}

while IFS= read -r line; do
  if [[ "$line" =~ ^#\ Sample:\ (.+)$ ]]; then
    flush_sample
    current_sample="${BASH_REMATCH[1]}"
  elif [[ "$line" =~ ^wget ]]; then
    # Extract URL (last field)
    url=$(echo "$line" | awk '{print $NF}')
    current_urls+=("$url")
  fi
done < "$WGET_SCRIPT"

flush_sample

echo "----------------------------------------"
echo "All downloads complete. Files are in '$OUTPUT_DIR/'"