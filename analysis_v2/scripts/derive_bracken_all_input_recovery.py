#!/usr/bin/env python3
"""Add all-input Bracken recovery beside legacy endpoints, without modifying inputs.

Community target counts are explicitly reconstructed from the frozen ordered
panel and integer allocator, then cross-checked against canonical target counts.
They are not independently measured per-target read-origin counts.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter
from decimal import Decimal
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table, parse_counts
from build_crc_cohort_canonical_input import allocate


def verify_checksums(root, filename, required):
    seen = set()
    for line in (root/filename).read_text().splitlines():
        checksum, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        require(Path(name).name == name and name not in seen, 'Unsafe/duplicate checksum member')
        require(digest(root/name) == checksum, 'Checksum mismatch: '+str(root/name))
        seen.add(name)
    require(set(required) <= seen, 'Required audit members not checksummed')


def close(a, b, tolerance=1e-12):
    return math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a)-float(b)) <= tolerance


def endpoints(baseline_count, observed_count, original_pairs, total_pairs,
              target_pairs, baseline_native, observed_native, legacy_total, legacy_target):
    require(original_pairs > 0 and total_pairs > original_pairs, 'Invalid pair totals')
    require(0 < target_pairs <= total_pairs-original_pairs, 'Invalid target pairs')
    b, o = baseline_count/original_pairs, observed_count/total_pairs
    retained = baseline_count/total_pairs
    signal = target_pairs/total_pairs
    expected = retained+signal
    recovered = o-retained
    old_retained = (1-legacy_total)*baseline_native
    old_expected = old_retained+legacy_target
    require(legacy_target > 0, 'Invalid legacy target fraction')
    return dict(baseline_all_input_abundance=b, observed_all_input_abundance=o,
                baseline_retained_all_input=retained, expected_all_input_abundance=expected,
                implanted_fraction_target_exact=signal,
                recovered_signal_all_input=recovered, recovery_ratio_all_input=recovered/signal,
                signed_error_all_input=o-expected, absolute_error_all_input=abs(o-expected),
                absolute_relative_error_all_input=abs(o-expected)/signal,
                legacy_baseline_native=baseline_native, legacy_observed_native=observed_native,
                legacy_expected_read_proportional=old_expected,
                legacy_recovered_signal=observed_native-old_retained,
                legacy_recovery_ratio=(observed_native-old_retained)/legacy_target,
                legacy_signed_error=observed_native-old_expected,
                legacy_absolute_error=abs(observed_native-old_expected),
                detected_native_nonzero=int(observed_native > 0),
                estimated_count_nonzero=int(observed_count > 0))


def run(a):
    audit, canonical, project, out = [p.resolve() for p in (a.audit, a.canonical, a.project, a.outdir)]
    require(not out.exists(), 'Output exists; use a fresh directory')
    require(a.community_allocation == 'frozen-panel-reconstruction', 'Explicit allocation mode required')
    verify_checksums(audit, 'SHA256SUMS', ['summary.json', 'profile_denominators.tsv',
                                         'input_hashes.tsv', 'design_evidence_inventory.tsv'])
    summary = json.loads((audit/'summary.json').read_text())
    require(summary['status'] == 'PASS_DENOMINATOR_AUDIT_ONLY' and summary['cohort'] == a.cohort,
            'Wrong/failed denominator audit')
    frozen = {}
    for line in (project/'analysis_v2/taxon_identity_freeze.sha256').read_text().splitlines():
        checksum, name = line.split(maxsplit=1)
        path = (project/name).resolve()
        require(project in path.parents and digest(path) == checksum, 'Identity freeze mismatch')
        frozen[path] = checksum
    panel_path, alias_path = project/'spikes/spike_panel.tsv', project/'examples/spike_taxon_aliases.csv'
    require(panel_path in frozen and alias_path in frozen, 'Panel/aliases not frozen')
    panel = table(panel_path)
    require(len(panel) == 10 and len({r['label'] for r in panel}) == 10, 'Expected ten unique panel targets')
    require(all(Decimal(r['weight']) > 0 and Decimal(r['weight']).is_finite() for r in panel), 'Invalid weights')
    targets = {r['label']: r for r in panel}
    aliases = {}
    with alias_path.open(newline='') as handle:
        for row in csv.DictReader(handle):
            if row['tool'] == 'kraken2_bracken':
                require(row['canonical'] not in aliases, 'Duplicate alias')
                aliases[row['canonical']] = row['alias']
    require(all(r['taxon_name'] in aliases for r in panel), 'Missing alias')

    profiles = table(audit/'profile_denominators.tsv')
    paths = [Path(r['profile']).resolve() for r in profiles]
    require(len(paths) == len(set(paths)) == summary['profiles'], 'Audit profile count/duplicates')
    sources = dict(frozen)
    sources[canonical] = digest(canonical)
    sources[Path(__file__).resolve()] = digest(Path(__file__))
    allocator = project/'analysis_v2/scripts/build_crc_cohort_canonical_input.py'
    sources[allocator] = digest(allocator)
    for row in table(audit/'input_hashes.tsv'):
        path = Path(row['path']).resolve()
        if row['role'] == 'profile':
            require(path not in sources or sources[path] == row['sha256'], 'Conflicting hashes')
            sources[path] = row['sha256']
    design_index = {}
    for row in table(audit/'design_evidence_inventory.tsv'):
        path = Path(row['path']).resolve()
        require(digest(path) == row['sha256'], 'Design hash mismatch')
        sources[path] = row['sha256']
        for design in table(path):
            category = 'community' if 'community' in design else 'independent'
            label = design['community'] if category == 'community' else design['label']
            key = (design['sample_id'], category, label, Decimal(design['fraction']))
            require(key not in design_index, 'Duplicate design')
            design_index[key] = (path, design)
    for source in (audit/'SHA256SUMS', audit/'summary.json', audit/'profile_denominators.tsv',
                   audit/'design_evidence_inventory.tsv', audit/'input_hashes.tsv'):
        sources[source] = digest(source)
    # Reject overlap with any consumed input tree before creating output.
    for source in sources:
        require(out != source and out not in source.parents and source.parent not in out.parents,
                'Output overlaps input directory: '+str(source.parent))

    canonical_rows = [r for r in table(canonical) if r['profiler'] == 'kraken2_bracken']
    canonical_index = {}
    for row in canonical_rows:
        require(row['include'] == '1' and row['cohort'] == a.cohort and row['assembly_arm'] == 'original',
                'Unexpected excluded/cohort/assembly canonical row')
        key = (Path(row['source_profile']).resolve(), row['analysis_population'], row['target_label'])
        require(key not in canonical_index, 'Duplicate canonical key')
        canonical_index[key] = row
    require({key[0] for key in canonical_index} == set(paths), 'Canonical/audit profile coverage mismatch')
    baselines = {}
    feature_cache = {}
    for row, path in zip(profiles, paths):
        report = path.with_name(path.name[:-len('.bracken.S.tsv')]+'.kraken2.report')
        for p in (path, report):
            require(p in sources and digest(p) == sources[p], 'Profile hash mismatch/absent')
        counts = parse_counts(report, path)
        for name, value in counts.items():
            require(close(row[name], value), 'Audit counts changed: '+name)
        features = {}
        for feature in table(path):
            name = feature['name'].strip()
            require(name not in features, 'Duplicate feature name')
            features[name] = (int(feature['new_est_reads']), float(feature['fraction_total_reads']))
        feature_cache[path] = features
        if row['category'] == 'baseline':
            require(row['sample_id'] not in baselines, 'Duplicate baseline')
            baselines[row['sample_id']] = (path, row)
    require(len(baselines) == summary['samples'], 'Baseline sample count mismatch')

    output, reconciliation, used_design, used_canonical = [], [], set(), set()
    for row, path in zip(profiles, paths):
        if row['category'] == 'baseline':
            continue
        sid, category = row['sample_id'], row['category']
        stem = path.name[:-len('.bracken.S.tsv')]
        require(stem.startswith(sid+'_'), 'Profile ID mismatch')
        label, fraction = stem[len(sid)+1:].rsplit('_f', 1)
        key = (sid, category, label, Decimal(fraction.replace('p', '.')))
        require(key in design_index and key not in used_design, 'Missing/duplicate design match')
        used_design.add(key)
        design_path, design = design_index[key]
        original = int(design['R'])
        added = int(design['N_total'] if category == 'community' else design['N_inserted'])
        total = int(row['total_observations'])
        base_path, base = baselines[sid]
        require(int(base['total_observations']) == original and total == original+added,
                'Pair count reconciliation failed')
        require(int(design['R1']) == int(design['R2']) == original and added > 0, 'Mate/count mismatch')
        require(close(design['f_hat'], added/total, 5.01e-9), 'Achieved dose outside eight-decimal rounding')
        allocations = allocate(added, panel) if category == 'community' else {label: added}
        require(sum(allocations.values()) == added and all(n > 0 for n in allocations.values()),
                'Nonpositive/incorrect allocation')
        reconciliation.append(dict(cohort=a.cohort, sample_id=sid, profile=str(path),
                                   original_pairs=original, inserted_pairs=added, total_pairs=total,
                                   status='PASS', allocation_basis='frozen_panel_integer_reconstruction'
                                   if category == 'community' else 'recorded_N_inserted'))
        for target, inserted in allocations.items():
            require(target in targets, 'Unknown target')
            ckey, bkey = (path, category, target), (base_path, category, target)
            require(ckey in canonical_index and bkey in canonical_index, 'Missing canonical target')
            c, b = canonical_index[ckey], canonical_index[bkey]
            used_canonical.update((ckey, bkey))
            for cr in (c, b):
                require(cr['sample_id'] == sid and cr['condition'] == row['condition'] and
                        cr['target_taxon'] == targets[target]['taxon_name'], 'Canonical identity mismatch')
            require(Path(c['source_design']).resolve() == design_path, 'Canonical design mismatch')
            require(int(c['implanted_read_pairs_target']) == inserted, 'Canonical target count mismatch')
            require(close(c['spike_fraction_total'], design['f_hat']), 'Canonical total dose mismatch')
            # Independent canonical fractions preserve eight-decimal f_hat; community
            # target fractions come from exact integer allocation / total pairs.
            expected_canonical = float(design['f_hat']) if category == 'independent' else inserted/total
            require(close(c['spike_fraction_target'], expected_canonical), 'Canonical target fraction mismatch')
            alias = aliases[targets[target]['taxon_name']]
            bn, bf = feature_cache[base_path].get(alias, (0, 0.0))
            on, of = feature_cache[path].get(alias, (0, 0.0))
            require(close(b['abundance_fraction'], bf) and close(c['abundance_fraction'], of),
                    'Canonical native abundance mismatch')
            require(float(b['spike_fraction_total']) == float(b['spike_fraction_target']) == 0,
                    'Nonzero canonical baseline dose')
            values = endpoints(bn, on, original, total, inserted, bf, of,
                               float(c['spike_fraction_total']), float(c['spike_fraction_target']))
            output.append(dict(cohort=a.cohort, sample_id=sid, condition=row['condition'],
                               analysis_population=category, target_label=target, target_alias=alias,
                               profile_id=stem, source_profile=str(path), source_baseline=str(base_path),
                               source_design=str(design_path), nominal_total_dose=str(key[-1]),
                               total_fraction_exact=added/total, target_pairs=inserted,
                               baseline_estimated_pairs=bn, observed_estimated_pairs=on,
                               baseline_total_pairs=original, observed_total_pairs=total,
                               allocation_basis=reconciliation[-1]['allocation_basis'], **values))
    require(used_design == set(design_index), 'Unmatched design rows')
    require(used_canonical == set(canonical_index), 'Unconsumed canonical target rows')
    require(output, 'No recovery endpoints')
    for source, checksum in sources.items():
        require(digest(source) == checksum, 'Input changed during run: '+str(source))
    out.mkdir(parents=True, exist_ok=False)
    write_table(out/'bracken_recovery_comparison.tsv', output)
    write_table(out/'pair_reconciliation.tsv', reconciliation)
    write_table(out/'input_hashes.tsv', [dict(path=str(p), sha256=h) for p, h in sorted(sources.items())])
    result = dict(cohort=a.cohort, samples=len(baselines), profiles=len(profiles),
                  spiked_profiles=len(reconciliation), target_endpoints=len(output),
                  populations=dict(Counter(r['analysis_population'] for r in output)),
                  status='PASS_ENDPOINT_CONSTRUCTION',
                  interpretation='Validation of construction, not proof of accurate recovery',
                  community_allocation='frozen_panel_integer_reconstruction_crosschecked_with_canonical',
                  native_profiles_modified=False, primary_DA_inputs_modified=False)
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.is_file()))
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True, choices=['yachida', 'feng', 'zeller'])
    for name in ('audit', 'canonical', 'project', 'outdir'):
        p.add_argument('--'+name, required=True, type=Path)
    p.add_argument('--community-allocation', required=True, choices=['frozen-panel-reconstruction'])
    args = p.parse_args()
    try:
        run(args)
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, '[FAIL] '+str(error)+'\n')
