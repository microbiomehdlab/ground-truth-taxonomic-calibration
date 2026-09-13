#!/usr/bin/env python3
import csv, subprocess, tempfile
from pathlib import Path
from test_canonical_input import HEADER, row

repo=Path(__file__).resolve().parents[2]
script=repo/'analysis_v2/scripts/combine_canonical_development_inputs.py'
validator=repo/'analysis_v2/scripts/validate_canonical_input.py'
with tempfile.TemporaryDirectory() as name:
    root=Path(name); inputs=[]
    for cohort,sample in [('feng','F1'),('yachida','Y1')]:
        d=root/cohort; d.mkdir(); table=d/'canonical_input.tsv'; rows=[]
        for profiler,base,spike in [('kraken2_bracken','0.01','0.011'),('metaphlan4','1','1.2')]:
            rows.extend([row(profiler,f'{sample}-spike-{profiler}','0',base,'0'),row(profiler,f'{sample}-spike-{profiler}','0.001',spike,'1000')])
        for values in rows:
            values[1]=cohort; values[3]=sample
        with table.open('w',newline='') as h:
            w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(HEADER); w.writerows(rows)
        subprocess.run(['python3',str(validator),'--input',str(table),'--outdir',str(d/'validation')],check=True)
        inputs.append(table)
    out=root/'combined'
    command=['python3',str(script),'--input',str(inputs[0]),'--input',str(inputs[1]),'--outdir',str(out)]
    subprocess.run(command,check=True)
    with (out/'canonical_input.tsv').open(newline='') as h: combined=list(csv.DictReader(h,delimiter='\t'))
    assert len(combined)==8 and {r['cohort'] for r in combined}=={'feng','yachida'}
    assert (out/'validation/SUCCESS').is_file() and (out/'DEVELOPMENT_ONLY.txt').is_file()
print('[PASS] combined canonical development input fixture')
