#!/usr/bin/env python3
import csv, subprocess, tempfile
from pathlib import Path
repo = Path(__file__).resolve().parents[2]; script = repo / "analysis_v2/scripts/build_development_sample_metadata.py"
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp); canonical = root/"canonical.tsv"; feng=root/"feng.tsv"; zeller=root/"zeller.tsv"; output=root/"metadata.tsv"
    canonical.write_text("cohort\tsample_id\tinclude\nfeng\tF1\t1\nzeller\tZ1\t1\nfeng\tF2\t0\n",encoding="utf-8"); header="sample_id\tcondition\tstudy\tage\tsex\tbmi\n"
    feng.write_text(header+"F1\tControl\tFeng\t60\tFemale\t25\nF2\tCRC\tFeng\t65\tMale\t26\n",encoding="utf-8"); zeller.write_text(header+"Z1\tCRC\tZeller\t70\tMale\t24\n",encoding="utf-8")
    subprocess.run(["python3",str(script),"--canonical",str(canonical),"--feng-manifest",str(feng),"--zeller-manifest",str(zeller),"--output",str(output)],check=True)
    with output.open(newline="",encoding="utf-8") as handle: rows=list(csv.DictReader(handle,delimiter="\t"))
    assert [row["sample_id"] for row in rows]==["F1","Z1"]; assert output.with_suffix(".tsv.sha256").is_file()
print("[PASS] development sample-metadata fixture")
