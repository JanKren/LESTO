#!/usr/bin/env python3
"""Collect the four upwind 864x108 (L4) runs of the lesto-submission-run package.

    python3 collect.py <package output directory>      # e.g. ~/lesto-submission-followup

For each base case of the iteration and mesh checks it copies the final (60 s) native wall profiles of
`<case>__grid864x108` to profiles/<case>/ and writes summary.json: the 41,472 -> 93,312-cell change of
the iodine, PbI2 and BiI3 wall profiles (relative L1 of the 1 cm wall amounts, normalised by the finer
level as in ../results.json), the total and centroid change, against the upwind 576x72 follow-up case
in ../../followup_results/profiles, plus the gates of the run.  The shifted-BiI case repeats the
upwind L4 run of the study (../profiles/steel_L4_up) independently; that difference is recorded too.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FOLLOW = HERE.parents[1] / 'followup_results/profiles'
STUDY_L4_UP = HERE.parent / 'profiles/steel_L4_up'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
SPECIES = ['element_I', 'PbI2_g', 'BiI3_g']


def amounts(path):
    lines = path.read_text().splitlines()
    head = next(l for l in lines if l.startswith('# columns')).split()[2:]
    data = np.array([[float(v) for v in l.split()] for l in lines if l.strip() and not l.startswith('#')])
    col = {k: i for i, k in enumerate(head)}
    return data[:, col['x']], data[:, col['wall']] * (data[:, col['x_hi']] - data[:, col['x_lo']])


def change(fine_dir, coarse_dir, species):
    """Relative change of the coarse profile, normalised by the fine (reference) level."""
    x, fine = amounts(fine_dir / f'{species}_native.dat')
    xc, coarse = amounts(coarse_dir / f'{species}_native.dat')
    assert np.allclose(x, xc)
    if fine.sum() <= 0:
        return None
    return {'profile_L1_pct': 100 * np.abs(fine - coarse).sum() / fine.sum(),
            'total_change_pct': 100 * (fine.sum() - coarse.sum()) / fine.sum(),
            'centroid_shift_cm': 100 * ((x * fine).sum() / fine.sum() - (x * coarse).sum() / coarse.sum())}


def main():
    output = Path(sys.argv[1]).expanduser()
    state = json.loads((output / 'run_state.json').read_text())['cases']
    solver = hashlib.sha256((output / 'bin/rhoFixedFlowFoam').read_bytes()).hexdigest()
    cases = []
    for base in BASES:
        name = f'{base}__grid864x108'
        case = output / 'cases' / name
        bench = json.loads((case / 'benchmark.json').read_text())
        report = json.loads((case / 'numerical-report.json').read_text())
        audit = json.loads((case / 'log-audit.json').read_text())
        schemes = (case / 'system/fvSchemes').read_text()
        assert bench['followup']['changes']['study'] == 'grid864x108' and audit['normal_end']
        assert 'Gauss upwind' in schemes and 'limitedLinear' not in schemes
        target = HERE / 'profiles' / base
        target.mkdir(parents=True, exist_ok=True)
        for species in SPECIES:
            shutil.copyfile(case / f'postProcessing/axialProfiles/60/{species}_native.dat',
                            target / f'{species}_native.dat')
        record = {
            'case': name, 'base_case': base, 'status': state[name]['status'],
            'finished_utc': state[name].get('finished_utc'), 'cells': 93312, 'scheme': 'Gauss upwind',
            'dt_s': bench['dt_s'], 'end_time_s': int(bench['end_time_s']),
            'carrier_audit': bench['carrier_audit']['status'],
            'element_gates': {k: v['status'] for k, v in report['element_gates'].items()},
            'max_relative_closure': max(v['max_relative_closure'] for v in report['element_gates'].values()),
            'max_respeciation_failures': report['max_respeciation_failures'],
            'BiI3_share_of_iodide_I': report['BiI3_fraction_of_modelled_iodide_wall_iodine'],
            'coarse_reference': f'run/024-liu2025/followup_results/profiles/{base}__grid576x72',
            'L3_to_L4': {s: change(target, FOLLOW / f'{base}__grid576x72', s) for s in SPECIES},
        }
        if base == BASES[3]:
            record['repeat_of_study_steel_L4_up'] = {s: change(STUDY_L4_UP, target, s) for s in SPECIES}
        cases.append(record)
    summary = {
        'what': 'Upwind 864x108 (93,312-cell) runs of the four iteration/mesh-check bases, made by the '
                'lesto-submission-run package; compared with the upwind 576x72 (41,472-cell) follow-up cases.',
        'solver_sha256': solver,
        'runner': 'lesto-submission-run package (run_study.py), output directory ~/lesto-submission-followup '
                  'on the workstation',
        'scope': 'Mesh sensitivity of the first-order upwind discretisation only; the production scheme is '
                 'limitedLinear, refined to L4 for the shifted-BiI case only (../profiles/steel_L4_ll).',
        'cases': cases,
    }
    (HERE / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    for c in cases:
        print(c['base_case'], c['status'], f"{c['L3_to_L4']['element_I']['profile_L1_pct']:.3f}")


if __name__ == '__main__':
    main()
