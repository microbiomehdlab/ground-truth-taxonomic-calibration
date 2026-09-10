#!/usr/bin/env python3
import pathlib,subprocess,tempfile
root=pathlib.Path(__file__).resolve().parents[2]; script=root/'analysis_v2/scripts/compare_biomarker_development_runs.py'
tables=[]
for line in script.read_text().splitlines():
    if line.strip().startswith('"models/'): tables.append(line.strip().strip('",'))
with tempfile.TemporaryDirectory() as t:
    t=pathlib.Path(t); left=t/'left'; right=t/'right'
    for rel in tables:
        for base,rows in [(left,['a\tb','x\t1.0000000000001','y\t2']), (right,['a\tb','y\t2.0','x\t1.0000000000002'])]:
            path=base/rel; path.parent.mkdir(parents=True,exist_ok=True); path.write_text('\n'.join(rows)+'\n')
    subprocess.run(['python3',str(script),'--reference',str(left),'--candidate',str(right),'--outdir',str(t/'pass')],check=True)
    assert (t/'pass/SUCCESS').is_file()
    target=right/tables[0]; target.write_text('a\tb\nx\t9\ny\t2\n')
    failed=subprocess.run(['python3',str(script),'--reference',str(left),'--candidate',str(right),'--outdir',str(t/'fail')])
    assert failed.returncode!=0 and not (t/'fail/SUCCESS').exists()
print('[PASS] biomarker development-run equivalence fixture')
