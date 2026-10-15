#!/usr/bin/env python3
"""Add one finished full-run corrector case of the lesto-submission-run package to summary.json.

    python3 collect.py <package output directory> <case>      # e.g. ~/lesto-submission-followup

Copies the final (60 s) native wall profiles to profiles/<case>/ and writes the case record: the
eight-minus-three profile change against ../full_analysis/profiles/<base case> (relative L1 of the
1 cm wall amounts, total change, centroid shift), the element-I gate and the outer-iteration counts.
The three earlier records of summary.json are reproduced by record() from the same package output.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FULL = HERE.parent / 'full_analysis/profiles'
SPECIES = ['element_I', 'PbI2_g', 'BiI3_g']


def amounts(path):
    rows = [l.split() for l in path.read_text().splitlines() if l.strip() and not l.startswith('#')]
    head = next(l for l in path.read_text().splitlines() if l.startswith('# columns')).split()[2:]
    data = np.array([[float(v) for v in r] for r in rows])
    col = {k: i for i, k in enumerate(head)}
    x = data[:, col['x']]
    return x, data[:, col['wall']] * (data[:, col['x_hi']] - data[:, col['x_lo']])


def change(folder, base, species):
    x, a = amounts(folder / f'{species}_native.dat')
    xb, b = amounts(base / f'{species}_native.dat')
    assert np.allclose(x, xb)
    if b.sum() <= 0:
        return None
    return {'profile_L1_pct': 100 * np.abs(a - b).sum() / np.abs(b).sum(),
            'total_change_pct': 100 * (a.sum() - b.sum()) / b.sum(),
            'centroid_shift_cm': 100 * ((x * a).sum() / a.sum() - (x * b).sum() / b.sum())}


def record(output, name):
    case = output / 'cases' / name
    bench = json.loads((case / 'benchmark.json').read_text())
    audit = json.loads((case / 'log-audit.json').read_text())
    report = json.loads((case / 'numerical-report.json').read_text())
    state = json.loads((output / 'run_state.json').read_text())['cases'][name]
    sha = json.loads((case / 'solver-process.json').read_text()).get('solver_sha256')
    base = bench['followup']['base_case']
    assert bench['followup']['changes']['correctors'] == 8 and audit['normal_end']
    profiles = case / 'postProcessing/axialProfiles/60'
    target = HERE / 'profiles' / name
    target.mkdir(parents=True, exist_ok=True)
    for species in SPECIES:
        shutil.copyfile(profiles / f'{species}_native.dat', target / f'{species}_native.dat')
    return {
        'case': name, 'base_case': base, 'status': state['status'], 'finished_utc': state['finished_utc'],
        'species_correctors': 8, 'reference_species_correctors': 3, 'cells': 18432,
        'scheme': 'Gauss upwind', 'dt_s': bench['dt_s'], 'end_time_s': int(bench['end_time_s']),
        'reference_profiles': f'run/024-liu2025/full_analysis/profiles/{base}',
        'element_I_gate': report['element_gates']['I'],
        'BiI3_share_of_iodide_I': report['BiI3_fraction_of_modelled_iodide_wall_iodine'],
        'steps_without_outer_convergence': {
            pair: {'unconverged': v['unconverged_steps'], 'steps': v['steps'],
                   'max_initial_residual': v['max_initial_residual'], 'tolerance': v['tolerance']}
            for pair, v in audit['pairs'].items()},
        'eight_minus_three': {s: change(target, FULL / base, s) for s in SPECIES},
        **({'solver_sha256': sha} if sha else {}),
    }


def main():
    output, name = Path(sys.argv[1]).expanduser(), sys.argv[2]
    summary_path = HERE / 'summary.json'
    summary = json.loads(summary_path.read_text())
    new = record(output, name)
    new.pop('solver_sha256', None)
    summary['cases'] = [c for c in summary['cases'] if c['case'] != name] + [new]
    summary.pop('not_run', None)
    summary_path.write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(new['eight_minus_three'], indent=1))


if __name__ == '__main__':
    main()
