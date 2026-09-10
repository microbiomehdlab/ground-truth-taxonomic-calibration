#!/usr/bin/env python3
import csv
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "MANUSCRIPT_RESULTS_CHECKLIST.tsv"
with path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
required = {"claim_id", "question", "estimand", "primary_source", "candidate_output", "required_seal", "status"}
assert rows and set(rows[0]) == required
assert len({row["claim_id"] for row in rows}) == len(rows)
assert all(row["status"] in {"IMPLEMENTED", "PENDING"} for row in rows)
assert all(all(row[field].strip() for field in required) for row in rows)
print(f"[PASS] manuscript-results checklist: {len(rows)} claims")
