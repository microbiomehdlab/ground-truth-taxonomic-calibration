#!/usr/bin/env python3
"""Read-only, receipt-selected denominator audit; never certifies spike recovery."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def table(path):
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        require(len(fields) == len(set(fields)), 'Duplicate columns: ' + str(path))
        rows = list(reader)
        require(all(None not in r and all(v is not None for v in r.values()) for r in rows),
                'Malformed TSV: ' + str(path))
        return rows


def verify_file(path, receipt):
    require(path in receipt, 'Absent from receipt: ' + str(path))
    row = receipt[path]
    require(path.is_file() and path.stat().st_size > 0, 'Missing/empty: ' + str(path))
    require(path.stat().st_size == int(row['bytes']), 'Size mismatch: ' + str(path))
    require(digest(path) == row['sha256'], 'Checksum mismatch: ' + str(path))


def parse_counts(report, bracken):
    roots = {}
    with report.open() as handle:
        for line in handle:
            f = line.rstrip('\n').split('\t')
            require(len(f) == 6, 'Expected six-column Kraken report: ' + str(report))
            if f[4].strip() in ('0', '1'):
                taxid = f[4].strip()
                require(taxid not in roots, 'Duplicate Kraken root')
                require(f[3].strip() == ('U' if taxid == '0' else 'R'), 'Wrong root rank')
                roots[taxid] = int(f[1])
    require(set(roots) == {'0', '1'}, 'Missing Kraken root/unclassified row')
    u, c = roots['0'], roots['1']
    require(u >= 0 and c >= 0 and u+c > 0, 'Invalid Kraken counts')
    rows = table(bracken)
    require(bool(rows), 'Empty Bracken table')
    seen = set()
    estimated, fraction_sum = 0, 0.0
    for row in rows:
        require(row['taxonomy_lvl'].strip() == 'S', 'Non-species Bracken row')
        taxid = row['taxonomy_id']
        require(taxid not in seen, 'Duplicate Bracken taxid')
        seen.add(taxid)
        n, f = int(row['new_est_reads']), float(row['fraction_total_reads'])
        require(n >= 0 and math.isfinite(f) and 0 <= f <= 1, 'Invalid Bracken number')
        estimated += n
        fraction_sum += f
    require(estimated <= c, 'Estimated species counts exceed classified observations')
    total = u+c
    return dict(total_observations=total, classified=c, unclassified=u,
                estimated_species_counts=estimated, classified_not_in_species=c-estimated,
                unclassified_fraction=u/total, represented_fraction=estimated/total,
                native_fraction_sum=fraction_sum, species_rows=len(rows))


def write_table(path, rows):
    require(bool(rows), 'No rows for ' + str(path))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def audit(cohort, state, results, out):
    state, results, out = (p.resolve() for p in (state, results, out))
    for source in (state, results):
        require(out != source and source not in out.parents and out not in source.parents,
                'Output must not overlap state/results')
    require(not out.exists(), 'Output already exists; use a fresh directory')
    seal = state / 'production_seal_v2'
    require(not (seal / 'AUDIT_IN_PROGRESS').exists(), 'Source audit in progress')
    checked = set()
    for line in (seal / 'production_seal.sha256').read_text().splitlines():
        checksum, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        require(Path(name).name == name and name not in checked, 'Unsafe/duplicate seal member')
        require(digest(seal / name) == checksum, 'Seal checksum mismatch: ' + name)
        checked.add(name)
    require({'SUCCESS', 'sample_flow.tsv', 'production_manifest.tsv'} <= checked,
            'Required seal members not checksummed')
    success = dict(line.split('\t', 1) for line in (seal / 'SUCCESS').read_text().splitlines())
    require(success.get('status') == 'PASS' and success.get('cohort') == cohort, 'Wrong/failed seal')
    manifest, flow = table(seal/'production_manifest.tsv'), table(seal/'sample_flow.tsv')
    ids = [r['sample_id'] for r in manifest]
    require(len(ids) == len(set(ids)) == int(success['samples']), 'Manifest sample count/duplicates')
    require(len(flow) == len(ids) and {r['sample_id'] for r in flow} == set(ids), 'Flow mismatch')
    flow = {r['sample_id']: r for r in flow}
    output, sources, designs = [], [], []
    for sample in manifest:
        sid, study = sample['sample_id'], sample['study']
        require(Path(sid).name == sid and Path(study).name == study, 'Unsafe sample/study')
        require(flow[sid]['status'] == 'PASS', 'Sample flow not PASS: ' + sid)
        receipt_path = state/'samples'/(sid+'.retained_outputs.tsv')
        receipt = {}
        sources.append(dict(path=str(receipt_path), sha256=digest(receipt_path), role='receipt'))
        for row in table(receipt_path):
            path = Path(row['path'])
            require(path.is_absolute(), 'Receipt path is not absolute')
            path = path.resolve()
            require(path not in receipt, 'Duplicate receipt path')
            receipt[path] = row
        root = results/study/sid/'profiles'
        counts = Counter()
        sample_rows = []
        for path in sorted(receipt):
            if not path.name.endswith('.bracken.S.tsv'):
                continue
            require(root in path.parents, 'Profile outside expected sample root: ' + str(path))
            category = path.relative_to(root).parts[0]
            require(category in ('baseline', 'community', 'independent'), 'Unknown profile category')
            counts[category] += 1
            report = path.with_name(path.name[:-len('.bracken.S.tsv')]+'.kraken2.report')
            for source in (path, report):
                verify_file(source, receipt)
                sources.append(dict(path=str(source), sha256=receipt[source]['sha256'], role='profile'))
            row = dict(cohort=cohort, sample_id=sid, condition=sample['condition'],
                       category=category, profile=str(path), **parse_counts(report, path))
            sample_rows.append(row)
        for category in ('baseline', 'community', 'independent'):
            require(counts[category] == int(flow[sid]['expected_'+category+'_profiles']),
                    'Profile count mismatch: ' + sid + ' ' + category)
        baseline = [r for r in sample_rows if r['category'] == 'baseline']
        require(len(baseline) == 1, 'Expected exactly one baseline')
        for row in sample_rows:
            row['delta_unclassified_fraction'] = row['unclassified_fraction']-baseline[0]['unclassified_fraction']
            row['delta_represented_fraction'] = row['represented_fraction']-baseline[0]['represented_fraction']
            row['delta_total_observations'] = row['total_observations']-baseline[0]['total_observations']
            row['spike_pair_reconciliation'] = 'NOT_YET_VERIFIED'
        output.extend(sample_rows)
        # Inventory design evidence without guessing its read/pair or dose semantics.
        for path in sorted(receipt):
            if path.suffix == '.tsv' and 'design' in str(path.relative_to(results)) if results in path.parents else False:
                verify_file(path, receipt)
                with path.open() as handle:
                    header = handle.readline().rstrip('\n')
                designs.append(dict(cohort=cohort, sample_id=sid, path=str(path),
                                    sha256=receipt[path]['sha256'], columns=header.replace('\t', '|')))
    require(output, 'No profiles selected')
    # No output directory is created until all input validations pass.
    out.mkdir(parents=True, exist_ok=False)
    write_table(out/'profile_denominators.tsv', output)
    write_table(out/'input_hashes.tsv', sources)
    if designs:
        write_table(out/'design_evidence_inventory.tsv', designs)
    summary = dict(cohort=cohort, samples=len(ids), profiles=len(output),
                   categories=dict(Counter(r['category'] for r in output)),
                   status='PASS_DENOMINATOR_AUDIT_ONLY',
                   spike_pair_reconciliation='NOT_YET_VERIFIED',
                   recovery_validated=False, script_sha256=digest(Path(__file__)))
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.is_file()))
    print(json.dumps(summary))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True, choices=['yachida', 'feng', 'zeller'])
    for flag in ('state-dir', 'results-root', 'outdir'):
        p.add_argument('--'+flag, required=True, type=Path)
    a = p.parse_args()
    try:
        audit(a.cohort, a.state_dir, a.results_root, a.outdir)
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, '[FAIL] '+str(error)+'\n')
