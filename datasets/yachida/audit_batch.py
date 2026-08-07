#!/usr/bin/env python3
"""Verify every persistent sample receipt and seal a completed batch."""
from __future__ import annotations
import argparse, csv, hashlib, pathlib
from stream_sample import verify_receipt

def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=pathlib.Path)
    p.add_argument("--scratch-root", required=True, type=pathlib.Path)
    p.add_argument("--state-dir", required=True, type=pathlib.Path)
    args=p.parse_args()
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle, delimiter="\t"))
    if not rows: raise SystemExit("[ERROR] Empty batch manifest")
    batch_ids={row.get("batch_id", "") for row in rows}
    if len(batch_ids)!=1 or not next(iter(batch_ids)):
        raise SystemExit(f"[ERROR] Expected one nonempty batch_id; found {batch_ids}")
    batch_id=next(iter(batch_ids)); sample_state=args.state_dir.resolve()/"samples"
    summary=[]
    for row in rows:
        sample=row["sample_id"]; receipt=sample_state/f"{sample}.retained_outputs.tsv"
        marker=sample_state/f"{sample}.verified"
        if not marker.is_file(): raise SystemExit(f"[ERROR] Missing verified marker: {sample}")
        outputs=verify_receipt(receipt, args.scratch_root.resolve()/sample)
        summary.append({"sample_id":sample,"condition":row["Target_Condition"],"outputs":str(len(outputs)),"status":"PASS"})
        print(f"[OK] {sample}: {len(outputs)} retained outputs")
    out=args.state_dir.resolve()/"batches"/batch_id; out.mkdir(parents=True, exist_ok=True)
    table=out/"batch_completion.tsv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer=csv.DictWriter(handle, fieldnames=["sample_id","condition","outputs","status"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    (out/"batch_completion.tsv.sha256").write_text(f"{digest}  batch_completion.tsv\n", encoding="utf-8")
    (out/"SUCCESS").write_text(f"batch_id\t{batch_id}\nsamples\t{len(summary)}\n", encoding="utf-8")
    print(f"[PASS] Batch sealed: {batch_id} ({len(summary)} samples)")

if __name__ == "__main__": main()
