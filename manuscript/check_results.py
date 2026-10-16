#!/usr/bin/env python3
"""Check saved-result arithmetic, provenance and the compiled review manuscript."""
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / 'thermochemistry/benchmarks/liu2025'))
from report import rows


def profile(root, name, species):
    filename = {'I': 'element_I', 'PbI2': 'PbI2_g', 'BiI3': 'BiI3_g'}[species]
    data = rows(root / 'profiles' / name / (filename + '_native.dat'))
    return np.array([r['x'] for r in data]), np.array([
        r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])


def amounts(folder, filename='element_I'):
    data = rows(folder / (filename + '_native.dat'))
    return np.array([r['x'] for r in data]), np.array([r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])


def l1(a, b):
    return 100 * np.abs(a - b).sum() / np.abs(b).sum()


def check_production(full, follow):
    """19 production cases: gates, stoichiometry of the saved profiles."""
    prod = REPO / 'run/024-liu2025/production_576x72_limitedLinear'
    analysis = json.loads((prod / 'analysis.json').read_text())
    ledger = json.loads((prod / 'ledger.json').read_text())
    assert (analysis['PASS'], analysis['FAIL'], analysis['NOTRUN']) == (15, 0, 0)
    assert analysis['cells'] == 41472 and analysis['scheme'] == 'Gauss limitedLinear 1'
    assert ledger['n_pass'] == 19 and set(ledger['status'].values()) == {'PASS'}
    names = [r['case'] for r in analysis['case_metrics']] + analysis['silica_alignment_cases']
    assert sorted(names) == sorted(ledger['status']) and len(names) == 19
    for name in names:
        _, iodine = profile(prod, name, 'I')
        _, lead = profile(prod, name, 'PbI2')
        _, bismuth = profile(prod, name, 'BiI3')
        assert np.max(np.abs(iodine - 2 * lead - 3 * bismuth)) / iodine.max() < 1e-12, name
    return len(names)


def check_scheme_study(full, follow):
    """Recompute every successive and cross-scheme profile difference of the study from its profiles."""
    study = REPO / 'run/024-liu2025/scheme_mesh_study'
    result = json.loads((study / 'results.json').read_text())
    bases = {'pure': 'PbI2_SS_Q1_baseline', 'steel': 'LBE-I_SS_II_Q2_shift29k_sat0.01'}
    def folder(key, level, scheme):
        if scheme == 'up' and level < 4:
            base = bases[key]
            return {1: follow / 'profiles' / (base + '__grid256x32'), 2: full / 'profiles' / base,
                    3: follow / 'profiles' / (base + '__grid576x72')}[level]
        return study / 'profiles' / f'{key}_L{level}_{scheme}'
    checked = 0
    for family, record in result.items():
        if family == 'cross_scheme_I_profile_L1_pct':
            continue
        key, scheme = family.split('.')
        levels = record['levels']
        for a, b in zip(levels[:-1], levels[1:]):
            value = l1(amounts(folder(key, a, scheme))[1], amounts(folder(key, b, scheme))[1])
            assert np.isclose(value, record['successive'][f'L{a}->L{b}']['I_profile_L1_pct'], rtol=1e-10), (family, a)
            checked += 1
    for name, expected in result['cross_scheme_I_profile_L1_pct'].items():
        key, level, pair = name.split('.')
        scheme = pair.split('_vs_')[0]
        level = int(level[1:])
        value = l1(amounts(folder(key, level, scheme))[1], amounts(folder(key, level, 'up'))[1])
        assert np.isclose(value, expected, rtol=1e-10), name
        checked += 1
    for key, base in bases.items():   # the L3 limitedLinear runs are the production runs
        for f in ['element_I', 'PbI2_g', 'BiI3_g']:
            assert (folder(key, 3, 'll') / (f + '_native.dat')).read_bytes() == (
                REPO / 'run/024-liu2025/production_576x72_limitedLinear/profiles' / base / (f + '_native.dat')).read_bytes()
    table = (HERE / 'generated/numerical-table.tex').read_text()
    for family in ['pure.ll', 'steel.ll', 'steel.up', 'steel.vl', 'pure.vl']:
        for step in result[family]['successive'].values():
            assert f"{step['I_profile_L1_pct']:.2f}" in table, family
    return checked


def check_upwind_l4(follow, generated):
    """Upwind 93,312-cell package runs: recompute the 41,472 -> 93,312 changes and their table cells."""
    root = REPO / 'run/024-liu2025/scheme_mesh_study/upwind_L4'
    summary = json.loads((root / 'summary.json').read_text())
    assert summary['solver_sha256'] == generated['solver_sha256'], 'upwind L4 solver differs'
    assert len(summary['cases']) == 4
    table = (HERE / 'generated/numerical-table.tex').read_text()
    for case in summary['cases']:
        base = case['base_case']
        assert case['status'] == 'PASS' and case['scheme'] == 'Gauss upwind' and case['cells'] == 93312
        assert set(case['element_gates'].values()) == {'PASS'} and case['max_respeciation_failures'] == 0
        for species in ['element_I', 'PbI2_g', 'BiI3_g']:
            x, fine = amounts(root / 'profiles' / base, species)
            xc, coarse = amounts(follow / 'profiles' / (base + '__grid576x72'), species)
            assert np.allclose(x, xc)
            record = case['L3_to_L4'][species]
            if fine.sum() > 0:
                assert np.isclose(l1(coarse, fine), record['profile_L1_pct'], rtol=1e-9), (base, species)
            else:
                assert record is None
        assert f"{case['L3_to_L4']['element_I']['profile_L1_pct']:.2f}" in table, base
    repeat = summary['cases'][3]['repeat_of_study_steel_L4_up']
    assert all(v['profile_L1_pct'] == 0 for v in repeat.values()), 'package L4 repeat differs from steel_L4_up'
    return len(summary['cases'])


def check_iteration(full):
    """Eight versus three correctors over 60 s: recompute the profile differences."""
    root = REPO / 'run/024-liu2025/iteration_check'
    summary = json.loads((root / 'summary.json').read_text())
    assert len(summary['cases']) == 4
    assert max(c['eight_minus_three']['element_I']['profile_L1_pct'] for c in summary['cases']) == \
        summary['max_iodine_profile_L1_pct']
    for case in summary['cases']:
        assert case['status'] == 'PASS' and case['species_correctors'] == 8
        for species in ['element_I', 'PbI2_g', 'BiI3_g']:
            x, a = amounts(root / 'profiles' / case['case'], species)
            xb, b = amounts(full / 'profiles' / case['base_case'], species)
            assert np.allclose(x, xb)
            record = case['eight_minus_three'][species]
            if b.sum() > 0:
                assert np.isclose(l1(a, b), record['profile_L1_pct'], rtol=1e-9, atol=1e-18), case['case']
            else:
                assert record is None
    return len(summary['cases'])


def check_lobresco(generated):
    """Independent recomputation of the deposition-onset table from the repository copies."""
    bench = REPO / 'thermochemistry/benchmarks/lobresco2026'
    run = REPO / 'run/023-lobresco2026'
    provenance = json.loads((bench / 'digitized/provenance.json').read_text())
    rows_out = {r['case']: r for r in generated['lobresco2026']['rows']}
    def first_crossing(x, y, wall_x, wall_t):
        threshold = 0.01 * y.max()
        k = next(i for i, v in enumerate(y) if v >= threshold)
        xo = x[k - 1] + (threshold - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1])
        return np.interp(xo, wall_x, wall_t)
    squares = []
    inlet_squares = []
    recomputed = 0
    for level in ['prod', 'coarse']:
        for case_dir in sorted((run / 'cases' / level).iterdir()):
            info = json.loads((case_dir / 'case.json').read_text())
            name = info['case']
            assert info['solver_sha256'].startswith('16fd804e') and info['carrier_audit']['status'] == 'PASS'
            ledger = json.loads((case_dir / 'ledger_summary.json').read_text())
            assert ledger['max_abs_closure_over_initial'] < 1e-15
            wall = np.loadtxt(run / 'inputs' / f"wallT_{info['table1']['flow_mL_min']}mLmin.csv", delimiter=',', skiprows=1)
            profile_file = case_dir / 'final_profiles/PbI2_g_gamma.dat'
            data = rows(profile_file)
            x = np.array([r['x'] for r in data]); deposit = np.array([r['wall'] for r in data])
            keep = x > info['boat']['x_m'][1]
            onset = first_crossing(x[keep], deposit[keep], wall[:, 0], wall[:, 1])
            if level == 'prod':
                inlet = first_crossing(x, deposit, wall[:, 0], wall[:, 1])
                assert abs(inlet - rows_out[name]['ours_from_inlet_K']) < 1e-9, name
                inlet_squares.append((inlet - provenance['table2_verbatim'][name]['experiment_K']) ** 2)
            table2 = provenance['table2_verbatim'][name]
            key = 'ours_minus_exp_K' if level == 'prod' else 'ours_coarse_minus_exp_K'
            assert abs(onset - table2['experiment_K'] - rows_out[name][key]) < 1e-9, (level, name)
            recomputed += 1
            if level == 'prod':
                squares.append((onset - table2['experiment_K']) ** 2)
    assert len(squares) == 5
    rms = (sum(squares) / 5) ** 0.5
    assert abs(rms - generated['lobresco2026']['summary']['rms_ours_K']) < 1e-9
    assert abs((sum(inlet_squares) / 5) ** 0.5 - generated['lobresco2026']['summary']['rms_ours_from_inlet_K']) < 1e-9
    table = (HERE / 'generated/lobresco-table.tex').read_text()
    assert f'& {rms:.0f} &' in table
    for check in provenance['tdep_checks']:
        assert abs(check['difference_K']) < 1.1, check['case']
    for path in list(run.rglob('*')) + list(bench.rglob('*')):
        if path.is_file() and path.suffix in {'.json', '.md', '.py', '.sh', '.csv'}:
            assert 'claude' not in path.read_text(errors='ignore').lower(), path
    assert recomputed == 10
    return recomputed


def check_graetz(generated):
    """Graetz-Robin verification: recompute observed orders and GCI from the tracked per-run fits and compare
    them with the solver-test summary table (table.md) and the manuscript tables."""
    root = HERE / 'data/graetz'
    provenance = json.loads((root / 'provenance.json').read_text())
    assert provenance['runs_made_with_this_solver'], 'Graetz runs not made with the recorded solver'
    assert provenance['solver_sha256_at_collection'] == generated['solver_sha256'], 'Graetz solver differs'
    criteria = (root / 'criteria.log').read_text().splitlines()
    assert criteria and all(line.startswith('PASS') for line in criteria), 'Graetz criterion not passed'
    reference = {r['Bi']: r for r in json.loads((root / 'reference.json').read_text())}
    for r in reference.values():
        values = list(r['methods'].values())
        assert max(abs(v / r['m'] - 1) for v in values) <= 1.01 * r['bound'] and r['bound'] < 1e-10
    fits = {(f['scheme'], f['outlet'], f['Bi'], f['law'], f['nr']): f
            for f in json.loads((root / 'fits.json').read_text()) if f['scheme'] in ('linear', 'SuperBee')}
    summary = {}
    for line in (root / 'table.md').read_text().splitlines():
        cells = [c.strip() for c in line.strip('|').split('|')]
        if len(cells) == 10 and '/' in cells[4]:
            summary[(cells[0], cells[1], float(cells[2]), cells[3], cells[4])] = (float(cells[5]), cells[8])
    main_table = (HERE / 'generated/graetz-table.tex').read_text()
    checked = 0
    for (scheme, outlet, bi, law, levels), (p_table, gci_table) in summary.items():
        mu = [fits[(scheme, outlet, bi, law, int(n))]['mu'] for n in levels.split('/')]
        e21, e32 = mu[1] - mu[2], mu[0] - mu[1]
        p = abs(np.log(abs(e32 / e21))) / np.log(2.0)
        gci = 125 * abs((mu[2] - mu[1]) / mu[2]) / (2.0 ** p - 1)
        assert abs(p - p_table) < 6e-5, (scheme, outlet, bi, law, levels)
        assert f'{gci:.3g}' == f"{float(gci_table.split()[0]):.3g}", (scheme, outlet, bi, law, levels, gci)
        assert e32 / e21 > 0, 'non-monotone Graetz triplet'
        if scheme == 'linear' and outlet == 'zeroGradient' and levels == '64/128/256':
            assert f'& {p:.3f} & {gci:.3g} \\\\' in main_table, (bi, law)
        checked += 1
    for (scheme, outlet, bi, law, nr), f in fits.items():
        assert abs(f['mu'] / f['muH'] - 1) < 1e-7, (scheme, outlet, bi, law, nr)
    assert checked == 20, checked
    return checked


def macro_values():
    text = (HERE / 'generated/numbers.tex').read_text()
    return dict(re.findall(r'\\newcommand\{\\(\w+)\}\{(.*)\}', text))


def check_damkohler(generated):
    """Damkoehler sweep (case 016): the record against case 016's own expectations, the csv against the
    json record, and the manuscript numbers recomputed from the per-column records."""
    root = HERE / 'data/damkohler'
    provenance = json.loads((root / 'provenance.json').read_text())
    assert provenance['solver_sha256'] == generated['solver_sha256'], 'Damkoehler solver differs'
    record = json.loads((root / 'sweep_summary.json').read_text())
    runs, refined = record['sweep']['runs'], record['refined_100x48']['runs']
    assert len(runs) == 13 and len(refined) == 5
    expect = {c['id'].split(' (')[0]: c for c in json.loads(
        (REPO / 'run/016-hks-channel-halfcell/expect.json').read_text())['checks']}
    a1 = next(r for r in runs if r['alpha'] == 1.0)
    n0 = expect['released [mol] = n0']['value']
    assert abs(a1['deposited_fraction'] * n0 / expect['wall inventory [mol]']['value'] - 1) < 1e-9
    assert abs(a1['T_onset'] - expect['1 % onset of the PbI2 deposit']['value']) < 0.01
    assert abs(a1['T_peak'] - expect['peak bin temperature']['value']) < 1e-6
    assert max(r['max_relative_closure'] for r in runs + refined) < 1e-14
    table = {float(r['alpha']): r for r in csv.DictReader((root / 'damkohler_sweep.csv').open())}
    for r in runs:
        t = table[float(f"{r['alpha']:.5g}")]
        for key, column in [('deposited_fraction', 'deposited_fraction'), ('centroid_mm', 'centroid_mm'),
                            ('Da_at_561K', 'Da_561K'), ('width_10_90_mm', 'width_10_90_mm')]:
            assert f'{r[key]:.5g}' == t[column], (r['alpha'], key)
        assert f"{r['t1.5']['Jhat'][84]:.5g}" == t['Jhat_84p5'], r['alpha']
    plateau = np.pi ** 2 / 4
    kinetic = max(abs(r['t1.5']['Jhat'][i] / r['t1.5']['Da_loc'][i] - 1)
                  for r in runs for i in (84, 94) if r['t1.5']['Da_loc'][i] <= 0.01)
    diffusion = max(abs(r['t1.5']['Jhat'][i] / plateau - 1)
                    for r in runs for i in (84, 94) if r['t1.5']['Da_loc'][i] >= 250)
    series = max(abs(r['t1.5']['Jhat'][i] * (1 / r['t1.5']['Da_loc'][i] + 1 / plateau) - 1)
                 for r in runs for i in (84, 94))
    by_alpha = {r['alpha']: r for r in runs}
    ref = np.array(a1['wall_profile'])
    prof = lambda a: 100 * np.abs(np.array(by_alpha[a]['wall_profile']) - ref).sum() / ref.sum()
    m = macro_values()
    assert m['DaKinDev'] == f'{100 * kinetic:.1f}' and m['DaDiffDev'] == f'{100 * diffusion:.1f}'
    assert m['DaSeriesDev'] == f'{100 * series:.0f}'
    assert m['DaProfOne'] == f'{prof(0.1):.1f}' and m['DaProfTwo'] == f'{prof(0.01):.0f}'
    assert m['DaCaptureLow'] == f"{100 * runs[0]['deposited_fraction']:.1f}"
    assert m['DaCaptureHigh'] == f"{100 * a1['deposited_fraction']:.1f}"
    assert m['DaFluxDeficit'] == f"{100 * (1 - by_alpha[0.01]['t1.5']['Jhat'][84] / plateau):.0f}"
    onsets = [r['T_onset'] for r in runs]
    assert max(onsets) - min(onsets) < 0.2
    refine = max(max(abs(q[k] / by_alpha[q['alpha']][k] - 1) for k in ('deposited_fraction', 'centroid_mm', 'width_10_90_mm'))
                 for q in refined)
    assert refine < 2e-3, refine
    return len(runs) + len(refined)


def check_production_sensitivities(follow):
    """Production-scheme D and accommodation variants: L1 recomputed from the tracked profiles against their
    production bases, the upwind columns against the follow-up record, and the manuscript table."""
    run = REPO / 'run/024-liu2025'
    root = run / 'production_sensitivities'
    rows_ = list(csv.DictReader((root / 'comparison.csv').open()))
    ledger = json.loads((root / 'ledger.json').read_text())
    up = {(r['case'], r['species']): r for r in follow['comparisons']}
    table = (HERE / 'generated/production-sensitivity-table.tex').read_text()
    for r in rows_:
        assert ledger[r['case']]['ended'] and ledger[r['case']]['end_time'] == 60.0, r['case']
        for species in ['I', 'PbI2', 'BiI3']:
            x, a = profile(root, r['case'], species)
            xb, b = profile(run / 'production_576x72_limitedLinear', r['base'], species)
            assert np.allclose(x, xb)
            if b.sum():
                assert np.isclose(l1(a, b), float(r[f'{species}_L1_pct']), rtol=1e-12), (r['case'], species)
        assert np.isclose(100 * up[(r['case'], 'I')]['profile_L1_relative_to_base'], float(r['upwind_I_L1_pct']),
                          rtol=1e-12), r['case']
        assert f"& {float(r['I_L1_pct']):.2f} & {float(r['upwind_I_L1_pct']):.2f} &" in table, r['case']
    m = macro_values()
    l1s = {r['variant']: [float(q['I_L1_pct']) for q in rows_ if q['variant'] == r['variant']] for r in rows_}
    assert m['PSDMin'] == f"{min(l1s['D5e-05'] + l1s['D0.0002']):.2f}"
    assert m['PSDMax'] == f"{max(l1s['D5e-05'] + l1s['D0.0002']):.2f}"
    assert m['PSAlphaTwoMax'] == f"{max(l1s['alpha0.01']):.2f}"
    return len(rows_)


def check_fields(generated):
    """PbI2 fields of the pure PbI2 column (data/fields): each record against the tracked wall profile of the same
    run, and the manuscript numbers recomputed from the npz without the helper of prepare_results."""
    root = HERE / 'data/fields'
    provenance = json.loads((root / 'provenance.json').read_text())
    m = macro_values()
    cells = 0
    for name, prefix in [('production', 'Field'), ('upwind_L4', 'FieldUp')]:
        record = provenance[name]
        assert record['solver_sha256'] == generated['solver_sha256'], (name, 'field solver differs')
        path = root / f'pbi2_field_{name}.npz'
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['npz_sha256'], name
        d = np.load(path)
        tracked = REPO / record['tracked_profile'] / 'element_I_native.dat'
        wall_I = float(re.search(r'# inventory gas \S+ wall (\S+)', tracked.read_text()).group(1))
        dx = (d['x_cm'][1] - d['x_cm'][0]) / 100
        assert abs(2 * d['wall_mol_per_m'].sum() * dx / wall_I - 1) < 1e-12, (name, 'deposit differs from profile')
        # the same 1 cm bins as the tracked profile; a wall face straddling a bin edge is split by overlap length,
        # as in axialProfiles.H
        data = rows(tracked)
        lo, hi = d['x_cm'] / 100 - dx / 2, d['x_cm'] / 100 + dx / 2
        overlap = lambda r: np.clip(np.minimum(hi, r['x_hi']) - np.maximum(lo, r['x_lo']), 0, None)
        binned = np.array([2 * (d['wall_mol_per_m'] * overlap(r)).sum() for r in data])
        reference = np.array([r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])
        assert np.abs(binned - reference).sum() / reference.sum() < 1e-12, (name, 'binned deposit differs')
        S = lambda key: d[key] / d['mole_eq_wall']
        zone = np.flatnonzero(d['wall_mol_per_m'] >= 0.01 * d['wall_mol_per_m'].max())
        z = slice(zone[0], zone[-1] + 1)
        assert m[prefix + 'WallSMax'] == f"{S('mole_wall_row')[z].max():.3f}", name
        assert m[prefix + 'CupSMax'] == f"{S('mole_cup')[z].max():.2f}", name
        assert m[prefix + 'AxisSMax'] == f"{S('mole_axis')[z].max():.2f}", name
        if name == 'production':
            assert m['FieldCupSMin'] == f"{S('mole_cup')[z].min():.2f}"
            assert m['FieldZoneLow'] == f"{d['x_cm'][zone[0]]:.1f}" and m['FieldZoneHigh'] == f"{d['x_cm'][zone[-1]]:.1f}"
        assert np.flatnonzero(S('mole_cup') >= 1)[0] == np.flatnonzero(d['wall_mol_per_m'] > 0)[0], (name, 'onset')
        assert np.all(S('mole_wall_row')[z] >= 1) and np.all(S('mole_axis')[z] >= S('mole_cup')[z]), name
        cells += int(d['mole_window'].size)
    return cells


def main():
    full = REPO / 'run/024-liu2025/full_analysis'
    follow = REPO / 'run/024-liu2025/followup_results'
    original = json.loads((full / 'analysis.json').read_text())
    result = json.loads((follow / 'analysis.json').read_text())
    assert (original['PASS'], original['FAIL'], original['NOTRUN']) == (15, 0, 1)
    assert result['completed_cases'] == 64 and result['execution_status'] == 'COMPLETE'
    assert not result['missing_or_failed']
    by_name = {r['case']: r for r in result['case_metrics']}
    nonzero = 0
    zero = 0
    for record in result['comparisons']:
        name = record['case']
        x, current = profile(follow, name, record['species'])
        xb, baseline = profile(full, by_name[name]['base_case'], record['species'])
        assert np.allclose(x, xb), name
        if baseline.sum():
            value = np.abs(current - baseline).sum() / baseline.sum()
            assert np.isclose(value, record['profile_L1_relative_to_base'], rtol=2e-12, atol=1e-15), name
            nonzero += 1
        else:
            assert record['profile_L1_relative_to_base'] is None, name
            zero += 1
    assert (nonzero, zero) == (180, 12)
    stoichiometric = 0
    for root, metrics in [(full, original['case_metrics']), (follow, result['case_metrics'])]:
        for case in metrics:
            _, iodine = profile(root, case['case'], 'I')
            _, lead = profile(root, case['case'], 'PbI2')
            _, bismuth = profile(root, case['case'], 'BiI3')
            defect = np.max(np.abs(iodine - 2 * lead - 3 * bismuth)) / max(iodine.max(), 1e-300)
            assert defect < 1e-12, (case['case'], defect)
            stoichiometric += 1
    hashes = 0
    for root, filename in [(follow, 'artifact_sha256.json'), (HERE, 'generated/artifact_sha256.json')]:
        for name, expected in json.loads((root / filename).read_text()).items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
            hashes += 1
    generated = json.loads((HERE / 'generated/results.json').read_text())
    for name, expected in generated['input_sha256'].items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected, name
        hashes += 1
    production = check_production(full, follow)
    scheme_checks = check_scheme_study(full, follow)
    iteration = check_iteration(full)
    upwind_l4 = check_upwind_l4(follow, generated)
    lobresco = check_lobresco(generated)
    graetz = check_graetz(generated)
    damkohler = check_damkohler(generated)
    psens = check_production_sensitivities(result)
    fields = check_fields(generated)
    source = '\n'.join(p.read_text() for p in [HERE / 'main.tex', HERE / 'supplementary.tex',
                                               *sorted((HERE / 'sections').glob('*.tex'))])
    from assess_resolution import evaluate
    resolution = json.loads((HERE / 'generated/resolution-audit.json').read_text())
    assert resolution == evaluate(), 'Saved-profile resolution reassessment is stale'
    for name, expected in resolution['input_sha256'].items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected, name
        hashes += 1
    assert not re.search(r'\\(?:todo|decide)\b', source)
    for latex_log in [HERE / 'build/main.log', HERE / 'build/supplementary.log']:
        log = latex_log.read_text()
        assert not re.search(r'undefined|Overfull|LaTeX Error|^!', log, re.MULTILINE), ('LaTeX issue', latex_log)
    for pof_log in [HERE / 'pof/build/lesto-pof.log', HERE / 'pof/build/supplementary.log']:
        pof_text = pof_log.read_text()
        assert not re.search(r'undefined|^Overfull|LaTeX Error|^!', pof_text, re.MULTILINE | re.IGNORECASE), pof_log
    # IJHMT working draft (arxiv.sty, line numbers) and its supplement: clean.
    for ijhmt_log in [HERE / 'ijhmt/build/lesto-ijhmt.log', HERE / 'ijhmt/build/supplementary.log']:
        assert not re.search(r'undefined|^Overfull|LaTeX Error|^!', ijhmt_log.read_text(),
                             re.MULTILINE | re.IGNORECASE), ijhmt_log
    # IJHMT Elsevier wrapper (build.sh --elsarticle), checked when built: one
    # 2.6 pt overfull box while page 1 is shipped out is a quirk of the class
    # (present without our content, invisible in the PDF); nothing else.
    elsarticle_log = HERE / 'ijhmt/build/lesto-ijhmt-elsarticle.log'
    elsarticle_overfull = []
    if elsarticle_log.exists():
        text = elsarticle_log.read_text()
        elsarticle_overfull = re.findall(r'^Overfull .*$', text, re.MULTILINE)
        assert len(elsarticle_overfull) <= 1 and all(
            re.match(r'Overfull \\hbox \([0-3]\.\d+pt too wide\) has occurred while \\output is active', o)
            for o in elsarticle_overfull), ('IJHMT elsarticle overfull', elsarticle_overfull)
        assert not re.search(r'undefined|LaTeX Error|^!', text, re.MULTILINE | re.IGNORECASE), elsarticle_log
    for blg in [HERE / 'build/main.blg', HERE / 'build/supplementary.blg']:
        bibliography = blg.read_text()
        assert 'Warning' not in bibliography and 'error' not in bibliography.lower(), ('Bibliography issue', blg)
    pdf = HERE / 'lesto-iodine-methods.pdf'
    pdf_info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)', pdf_info, re.MULTILINE).group(1))
    text = subprocess.check_output(['pdftotext', str(pdf), '-'], text=True)
    assert '[?]' not in text and '[TODO' not in text and '[DECIDE' not in text
    supplement = HERE / 'lesto-iodine-supplementary.pdf'
    supplement_pages = int(re.search(r'^Pages:\s+(\d+)', subprocess.check_output(['pdfinfo', str(supplement)], text=True),
                                     re.MULTILINE).group(1))
    supplement_text = subprocess.check_output(['pdftotext', str(supplement), '-'], text=True)
    assert '[?]' not in supplement_text and '??' not in supplement_text and 'TODO' not in supplement_text
    assert 'TODO' not in text and '??' not in text
    checks = {
        'status': 'PASS', 'baseline_cases': 15, 'followup_cases': 64,
        'independently_recomputed_profile_comparisons': nonzero,
        'undefined_zero_reference_comparisons_confirmed': zero,
        'iodine_stoichiometry_profiles_checked': stoichiometric,
        'production_cases_checked': production, 'scheme_mesh_differences_recomputed': scheme_checks,
        'iteration_check_cases_recomputed': iteration, 'upwind_L4_cases_recomputed': upwind_l4,
        'lobresco_onsets_recomputed': lobresco,
        'graetz_triplets_recomputed': graetz, 'damkohler_runs_checked': damkohler,
        'production_scheme_sensitivities_recomputed': psens, 'field_cells_checked': fields, 'supplement_pdf_pages': supplement_pages,
        'supplement_pdf_sha256': hashlib.sha256(supplement.read_bytes()).hexdigest(),
        'resolution_aggregation_comparisons': len(resolution['aggregation']),
        'thermodynamic_resolution_pairs': len(resolution['thermodynamic_resolution']),
        'input_and_artifact_hashes_checked': hashes, 'pdf_pages': pages,
        'undefined_latex_references': 0, 'overfull_latex_boxes': 0,
        'bibliography_warnings': 0, 'draft_placeholders': 0,
        'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'pof_overfull_latex_boxes': 0, 'pof_undefined_latex_references': 0,
        'ijhmt_overfull_latex_boxes': 0, 'ijhmt_undefined_latex_references': 0,
        'ijhmt_elsarticle_class_page1_overfull_boxes': len(elsarticle_overfull),
        'scope': 'Saved-results and manuscript consistency checks; not new CFD runs or experimental validation.'}
    (HERE / 'review-checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
