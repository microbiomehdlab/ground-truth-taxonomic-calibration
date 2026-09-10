#!/usr/bin/env python3
"""Verify scientific-table equivalence between sequential and map-reduce runs."""
import argparse, csv, hashlib, math
from pathlib import Path

TABLES = [
 "models/artificial_pooled/models/paired_da_results.tsv",
 "models/artificial_pooled/models/sample_feature_log2_changes.tsv",
 "models/artificial_pooled/evaluation/biomarker_propagation_metrics.tsv",
 "models/artificial_stratified/models/paired_da_results.tsv",
 "models/artificial_stratified/models/sample_feature_log2_changes.tsv",
 "models/artificial_stratified/evaluation/biomarker_propagation_metrics.tsv",
 "models/disease/models/primary_disease_da_results.tsv",
 "models/disease/models/sensitivity_bmi_disease_da_results.tsv",
 "models/disease/models/disease_da_exclusions.tsv",
 "models/disease/evaluation/disease_biomarker_propagation_metrics.tsv",
]

def normalized(path, digits):
    with path.open(newline='', encoding='utf-8') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if reader.fieldnames is None: raise ValueError(f'empty table: {path}')
        fields=tuple(reader.fieldnames); rows=[]
        for row in reader:
            values=[]
            for field in fields:
                value=row[field]
                try: number=float(value)
                except ValueError: values.append(value); continue
                values.append(format(number,f'.{digits}g') if math.isfinite(number) else value)
            rows.append(tuple(values))
    rows.sort(); digest=hashlib.sha256()
    digest.update(('\t'.join(fields)+'\n').encode())
    for row in rows: digest.update(('\t'.join(row)+'\n').encode())
    return fields,rows,digest.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True); p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--significant-digits',type=int,default=12); a=p.parse_args()
    if a.significant_digits < 8 or a.significant_digits > 17: p.error('significant digits must be 8..17')
    records=[]; failed=False
    for relative in TABLES:
        left=a.reference/relative; right=a.candidate/relative
        if not left.is_file() or not right.is_file():
            records.append((relative,'NA','NA','NA','MISSING')); failed=True; continue
        lf,lr,lh=normalized(left,a.significant_digits); rf,rr,rh=normalized(right,a.significant_digits)
        status='PASS' if lf==rf and lr==rr else 'FAIL'; failed |= status!='PASS'
        records.append((relative,len(lr),len(rr),f'{lh}:{rh}',status))
    a.outdir.mkdir(parents=True,exist_ok=False)
    report=a.outdir/'equivalence_audit.tsv'
    with report.open('w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['table','reference_rows','candidate_rows','normalized_sha256_pair','status']); w.writerows(records)
    (a.outdir/'audit_settings.tsv').write_text(f'setting\tvalue\nsignificant_digits\t{a.significant_digits}\nrow_order\tignored\nmissing_tables\tfail\n')
    if failed: raise SystemExit(f'[ERROR] Development runs are not equivalent; see {report}')
    (a.outdir/'SUCCESS').write_text(f'tables\t{len(TABLES)}\nstatus\tPASS\n')
    print(f'[PASS] {len(TABLES)} scientific tables are equivalent: {report}')
if __name__=='__main__': main()
