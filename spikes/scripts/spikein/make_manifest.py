#!/usr/bin/env python3

import csv
import os
import re
import sys

def extract_fastq_id(fastq_path):
    """
    Extract sample accession from FASTQ filename.
    Examples:
      ERR688505_1.fastq.gz -> ERR688505
      ERR688505_2.fastq.gz -> ERR688505
      ERR688505_R1.fastq.gz -> ERR688505
    """
    fname = os.path.basename(fastq_path)

    patterns = [
        r"(.+?)[._-][Rr]?1\.f(ast)?q\.gz$",
        r"(.+?)[._-][Rr]?2\.f(ast)?q\.gz$",
        r"(.+?)\.f(ast)?q\.gz$",
    ]

    for pattern in patterns:
        m = re.match(pattern, fname)
        if m:
            return m.group(1)

    raise ValueError(f"Could not extract sample_id from FASTQ filename: {fname}")


def load_metadata(metadata_file):
    """
    Load metadata into a dictionary keyed by Name.
    """
    meta = {}
    with open(metadata_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row["Name"]
            meta[key] = row
    return meta


def main(seq_file, metadata_file, output_file):
    metadata = load_metadata(metadata_file)

    with open(seq_file, newline="") as fin, open(output_file, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = ["sample_id", "fastq1", "fastq2", "original_id", "Target_Condition"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            original_id = row["sample_id"]
            fastq1 = row["fastq1"]
            fastq2 = row["fastq2"]

            new_sample_id = extract_fastq_id(fastq1)

            if original_id not in metadata:
                print(f"Warning: {original_id} not found in metadata", file=sys.stderr)
                target_condition = ""
            else:
                target_condition = metadata[original_id].get("Target Condition", "")

            writer.writerow({
                "sample_id": new_sample_id,
                "fastq1": fastq1,
                "fastq2": fastq2,
                "original_id": original_id,
                "Target_Condition": target_condition,
            })


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} sequencing.tsv metadata.tsv output.tsv")
        sys.exit(1)

    seq_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_file = sys.argv[3]

    main(seq_file, metadata_file, output_file)
