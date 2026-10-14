#!/usr/bin/env python3
"""Check saved-result arithmetic, provenance and the compiled review manuscript."""
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
    prod = REPO / 'run/034-liu2025/production_576x72_limitedLinear'
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
    study = REPO / 'run/034-liu2025/scheme_mesh_study'
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
                REPO / 'run/034-liu2025/production_576x72_limitedLinear/profiles' / base / (f + '_native.dat')).read_bytes()
    table = (HERE / 'generated/numerical-table.tex').read_text()
    for family in ['pure.ll', 'steel.ll', 'steel.up', 'steel.vl', 'pure.vl']:
        for step in result[family]['successive'].values():
            assert f"{step['I_profile_L1_pct']:.2f}" in table, family
    return checked


def check_iteration(full):
    """Eight versus three correctors over 60 s: recompute the profile differences."""
    root = REPO / 'run/034-liu2025/iteration_check'
    summary = json.loads((root / 'summary.json').read_text())
    assert len(summary['cases']) == 3
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
    run = REPO / 'run/035-lobresco2026'
    provenance = json.loads((bench / 'digitized/provenance.json').read_text())
    rows_out = {r['case']: r for r in generated['lobresco2026']['rows']}
    def first_crossing(x, y, wall_x, wall_t):
        threshold = 0.01 * y.max()
        k = next(i for i, v in enumerate(y) if v >= threshold)
        xo = x[k - 1] + (threshold - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1])
        return np.interp(xo, wall_x, wall_t)
    squares = []
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
            table2 = provenance['table2_verbatim'][name]
            key = 'ours_minus_exp_K' if level == 'prod' else 'ours_coarse_minus_exp_K'
            assert abs(onset - table2['experiment_K'] - rows_out[name][key]) < 1e-9, (level, name)
            recomputed += 1
            if level == 'prod':
                squares.append((onset - table2['experiment_K']) ** 2)
    assert len(squares) == 5
    rms = (sum(squares) / 5) ** 0.5
    assert abs(rms - generated['lobresco2026']['summary']['rms_ours_K']) < 1e-9
    table = (HERE / 'generated/lobresco-table.tex').read_text()
    assert f'& {rms:.0f} &' in table
    for check in provenance['tdep_checks']:
        assert abs(check['difference_K']) < 1.1, check['case']
    for path in list(run.rglob('*')) + list(bench.rglob('*')):
        if path.is_file() and path.suffix in {'.json', '.md', '.py', '.sh', '.csv'}:
            assert 'claude' not in path.read_text(errors='ignore').lower(), path
    assert recomputed == 10
    return recomputed


def main():
    full = REPO / 'run/034-liu2025/full_analysis'
    follow = REPO / 'run/034-liu2025/followup_results'
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
    lobresco = check_lobresco(generated)
    source = '\n'.join(p.read_text() for p in [HERE / 'main.tex', *sorted((HERE / 'sections').glob('*.tex'))])
    assert not re.search(r'\\(?:todo|decide)\b', source)
    log = (HERE / 'build/main.log').read_text()
    assert not re.search(r'undefined|Overfull|LaTeX Error|^!', log, re.MULTILINE), 'LaTeX issue'
    for pof_log in [HERE / 'pof/build/lesto-pof.log']:
        pof_text = pof_log.read_text()
        assert not re.search(r'undefined|^Overfull|LaTeX Error|^!', pof_text, re.MULTILINE | re.IGNORECASE), pof_log
    bibliography = (HERE / 'build/main.blg').read_text()
    assert 'Warning' not in bibliography and 'error' not in bibliography.lower(), 'Bibliography issue'
    pdf = HERE / 'lesto-iodine-methods.pdf'
    pdf_info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)', pdf_info, re.MULTILINE).group(1))
    text = subprocess.check_output(['pdftotext', str(pdf), '-'], text=True)
    assert '[?]' not in text and '[TODO' not in text and '[DECIDE' not in text
    checks = {
        'status': 'PASS', 'baseline_cases': 15, 'followup_cases': 64,
        'independently_recomputed_profile_comparisons': nonzero,
        'undefined_zero_reference_comparisons_confirmed': zero,
        'iodine_stoichiometry_profiles_checked': stoichiometric,
        'production_cases_checked': production, 'scheme_mesh_differences_recomputed': scheme_checks,
        'iteration_check_cases_recomputed': iteration, 'lobresco_onsets_recomputed': lobresco,
        'input_and_artifact_hashes_checked': hashes, 'pdf_pages': pages,
        'undefined_latex_references': 0, 'overfull_latex_boxes': 0,
        'bibliography_warnings': 0, 'draft_placeholders': 0,
        'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'pof_overfull_latex_boxes': 0, 'pof_undefined_latex_references': 0,
        'scope': 'Saved-results and manuscript consistency checks; not new CFD runs or experimental validation.'}
    (HERE / 'review-checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
