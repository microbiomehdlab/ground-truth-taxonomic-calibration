import copy, csv, gzip, json, math, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from link_da3_recovery import (build, normalize, link_context, validate_points,
                               manifest, load_endpoints, digest, write_table, dose)
from plan_da3_canary import finish

def fixture(root):
    plan=root/'plan';plan.mkdir();(plan/'batches').mkdir()
    b=root/'bracken';mp=root/'metaphlan'
    catalog=dict(profiles={},families={},targets={t:{'Target':'Fnuc'} for t in ('kraken2_bracken','metaphlan4')})
    for cohort in ('yachida','feng','zeller'):
        br=[];mr=[]
        for tool in catalog['targets']:
            catalog['families']['|'.join((cohort,tool,'Adenoma'))]=['Target']
            for i in range(6):
                sid='s%d'%i;baseline=.001*(i+1)
                for d in ('0','0.001','0.01'):
                    path=str(root/'profiles'/cohort/tool/sid/d)
                    catalog['profiles']['|'.join((cohort,tool,sid,d))]=dict(source_profile=path,sha256='fixture')
                    if d=='0':continue
                    signal=float(d)/10
                    recovered=signal*(0 if i==1 else -0.5 if i==0 and d=='0.01' else 1)
                    observed=baseline*(1-float(d))+recovered
                    row=dict(cohort=cohort,sample_id=sid,condition='Adenoma',analysis_population='community',
                             target_label='Fnuc',profile_id=sid+'_CRCpanel_f'+d.replace('.','p'),source_profile=path)
                    basepath=str(root/'profiles'/cohort/tool/sid/'0')
                    if tool=='kraken2_bracken':
                        row.update(target_alias='Target',source_baseline=basepath,nominal_total_dose=d,
                            legacy_baseline_native=baseline,legacy_observed_native=observed,
                            total_fraction_exact=d,implanted_fraction_target_exact=signal,
                            recovered_signal_all_input=recovered,recovery_ratio_all_input=recovered/signal,
                            signed_error_all_input=recovered-signal);br.append(row)
                    else:
                        row.update(profiler=tool,reference_type='genome_equivalent',source_baseline_profile=basepath,
                            baseline_abundance_fraction=baseline,observed_abundance_fraction=observed,
                            spike_fraction_total=d,spike_fraction_target=signal,implanted_signal_profiler_scale=signal,
                            recovered_spike_signal_profiler_scale=recovered,response_ratio_profiler_scale=recovered/signal,
                            response_residual_profiler_scale=recovered-signal);mr.append(row)
        folder=b/cohort;folder.mkdir(parents=True)
        write_table(folder/'bracken_recovery_comparison.tsv',br)
        (folder/'summary.json').write_text(json.dumps(dict(cohort=cohort,status='PASS_ENDPOINT_CONSTRUCTION')))
        finish(folder,'fixture')
        folder=mp/cohort/'endpoints';folder.mkdir(parents=True)
        write_table(folder/'paired_endpoints.tsv',mr)
    finish(mp,'PASS_ENDPOINTS')
    contexts=[]
    for tool in catalog['targets']:
        for arm,ds in [('N',['0']*3),('P25',['0.001','0','0']),
                       ('V',['0.01','0.001','0.01']),('PV',['0.01','0.001','0'])]:
            obs=[dict(sample_id='s%d'%i,group=int(i<3),dose=ds[i] if i<3 else '0',
                key='|'.join(('feng',tool,'s%d'%i,ds[i] if i<3 else '0'))) for i in range(6)]
            contexts.append(dict(context_id='c%d'%len(contexts),cohort='feng',profiler=tool,background='Adenoma',
                n=3,arm=arm,anchor='',allocation_id='a',requested_exposure='0.5',observations=obs))
    (plan/'catalog.json').write_text(json.dumps(catalog))
    (plan/'batches/a.json').write_text(json.dumps(contexts))
    write_table(plan/'tasks.tsv',[dict(index=0,file='batches/a.json',contexts=len(contexts))])
    finish(plan,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
    points=load_endpoints(b,mp,{})
    targets=[]
    for c in contexts:
        values=[]
        for o in c['observations']:
            key=('feng',c['profiler'],o['sample_id'],dose(o['dose'] or '0'),'Fnuc')
            if o['dose']!='0':v=points[key]['spiked_native']
            else:v=points[('feng',c['profiler'],o['sample_id'],dose('0.001'),'Fnuc')]['baseline_native']
            values.append(math.log2(1+v/1e-8))
        targets.append(dict(context_id=c['context_id'],target_label='Fnuc',feature='Target',
            beta=sum(values[:3])/3-sum(values[3:])/3,parametric_p=.02,parametric_bh_q=.2,max_statistic_fwer_p=.3))
    folder=root/'results/batch_00000';folder.mkdir(parents=True)
    write_table(folder/'summary.tsv',[dict(context_id=c['context_id'],family_n=1) for c in contexts])
    (folder/'resume_identity.json').write_text(json.dumps(dict(plan_sha256=digest(plan/'SHA256SUMS'))))
    with gzip.open(folder/'targets.tsv.gz','wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(targets[0]),delimiter='\t');w.writeheader();w.writerows(targets)
    finish(folder,'PASS_BATCH_COMPUTATION_NOT_PRODUCTION_AUTHORIZATION')
    return plan,root/'results',b,mp,contexts,targets

class LinkTests(unittest.TestCase):
    def test_end_to_end_partial_variable_null_and_denominators(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p,r,b,m,_,_=fixture(root)
            build(p,r,b,m,root/'out',fixture=True)
            manifest(root/'out',['SUCCESS','context_target_explanation.tsv.gz'],{})
            with gzip.open(root/'out/context_target_explanation.tsv.gz','rt') as f:
                rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),8)
            for row in rows:
                if row['arm']=='N':
                    self.assertEqual(row['exposed_cases'],'0')
                    self.assertEqual(row['exposed_recovery_ratio_median'],'')
                elif row['arm']=='P25':
                    self.assertEqual(row['exposed_cases'],'1')
                    self.assertEqual(float(row['exposed_recovery_ratio_median']),1)
                elif row['arm']=='PV':
                    self.assertEqual(row['exposed_cases'],'2')
                    self.assertEqual(float(row['exposed_recovery_ratio_median']),-.25)
                    self.assertEqual(row['exposed_negative_recovery'],'1')
                    self.assertEqual(row['exposed_zero_recovery'],'1')
            self.assertEqual(json.loads((root/'out/status.json').read_text())['context_target_rows'],8)
            with self.assertRaises(ValueError):build(p,r,b,m,p/'overlap',fixture=True)

    def test_missing_batch_and_tampered_endpoint_fail(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p,r,b,m,_,_=fixture(root)
            with self.assertRaises(FileNotFoundError):build(p,root/'missing',b,m,root/'out',fixture=True)
            self.assertFalse((root/'out').exists())
            f=b/'feng/bracken_recovery_comparison.tsv';f.write_text(f.read_text()+'\n')
            with self.assertRaises(ValueError):build(p,r,b,m,root/'out',fixture=True)

    def test_mismatch_and_no_pool_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p,r,b,m,contexts,targets=fixture(root)
            catalog=json.loads((p/'catalog.json').read_text())
            from decimal import Decimal
            catalog['_profiles_by_identity']={tuple(k.split('|')[:3])+(Decimal(k.split('|')[3]),):v for k,v in catalog['profiles'].items()}
            points=load_endpoints(b,m,{});baseline=validate_points(points,catalog)
            c=contexts[3];f=[targets[3]]
            broken=copy.deepcopy(c);broken['observations'][0]['dose']='0.001'
            with self.assertRaises(ValueError):link_context(broken,f,1,catalog,points,baseline)
            incomplete=dict(points);del incomplete[('feng','kraken2_bracken','s0',dose('.01'),'Fnuc')]
            with self.assertRaises(KeyError):link_context(c,f,1,catalog,incomplete,baseline)
            wrong=copy.deepcopy(f);wrong[0]['beta']+=1
            with self.assertRaises(ValueError):link_context(c,wrong,1,catalog,points,baseline)
            with self.assertRaises(ValueError):link_context(c,f+f,1,catalog,points,baseline)
            wrong=copy.deepcopy(c);wrong['observations'][-1]=wrong['observations'][0]
            with self.assertRaises(ValueError):link_context(wrong,f,1,catalog,points,baseline)
            row=next(v for v in points.values() if v['profiler']=='metaphlan4')
            self.assertEqual(row['reference_type'],'genome_equivalent')
            bad=copy.deepcopy(points)
            next(v for v in bad.values() if v['profiler']=='kraken2_bracken')['source_target_alias']='Unrelated'
            with self.assertRaises(ValueError):validate_points(bad,catalog)

    def test_wrong_metaphlan_scale_and_duplicate_endpoint(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p,r,b,m,_,_=fixture(root)
            from audit_bracken_denominators import table
            file=m/'feng/endpoints/paired_endpoints.tsv'
            rows=table(file);rows[0]['reference_type']='read_proportional'
            write_table(file,rows);finish(m,'PASS_ENDPOINTS')
            with self.assertRaises(ValueError):load_endpoints(b,m,{})
            rows[0]['reference_type']='genome_equivalent'
            write_table(file,rows+rows[:1]);finish(m,'PASS_ENDPOINTS')
            with self.assertRaises(ValueError):load_endpoints(b,m,{})

if __name__=='__main__':unittest.main()
