#!/usr/bin/env python3
"""Reassess saved profiles; no new CFD, fitted inputs or continuum error bounds."""
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
FULL = REPO / 'run/034-liu2025/full_analysis/profiles'
FOLLOW = REPO / 'run/034-liu2025/followup_results/profiles'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
LABELS = ['Pure PbI2', 'Steel Q2, 0.1\\%', 'Steel Q2, 1\\%', 'Steel Q2, shifted BiI']


def read_profile(path):
    lines = path.read_text().splitlines()
    columns = next(s[10:].split() for s in lines if s.startswith('# columns '))
    data = np.loadtxt(path)
    c = {k: data[:, i] for i, k in enumerate(columns)}
    amount = c['wall'] * (c['x_hi'] - c['x_lo'])
    if not np.isfinite(data).all() or np.any(amount < 0) or amount.sum() <= 0:
        raise ValueError(f'Invalid positive wall profile: {path}')
    if not np.allclose(np.diff(c['x']), .01, rtol=0, atol=1e-12):
        raise ValueError(f'Expected native 1 cm bins: {path}')
    return c['x'], amount


def aggregate(values, width, offset):
    # Include the edge groups: no cropping, smoothing, or mass loss.
    groups = np.floor_divide(np.arange(len(values)) - offset, width)
    return np.bincount(groups - groups.min(), weights=values)


def relative_l1(current, reference):
    return float(100 * np.abs(current - reference).sum() / reference.sum())


def shape_metrics(current, reference, x):
    p, q = current / current.sum(), reference / reference.sum()
    return {'total_variation_pct': float(50 * np.abs(p-q).sum()),
            'binned_W1_cm': float(100*np.sum(np.abs(np.cumsum(p-q)[:-1])*np.diff(x)))}


def evaluate():
    hashes, profiles, bins = {}, {}, []
    for base in BASES:
        profiles[base] = {}
        for level, root, suffix in [('coarse', FOLLOW, '__grid256x32'),
                                    ('medium', FULL, ''),
                                    ('fine', FOLLOW, '__grid576x72'),
                                    ('half_step', FOLLOW, '__dt1ms')]:
            folder = root / (base + suffix)
            values = {}
            for species in ['element_I', 'BiI3_g']:
                path = folder / (species + '_native.dat')
                hashes[path.relative_to(REPO).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
                if species == 'BiI3_g' and base == BASES[0]:
                    continue
                x, a = read_profile(path)
                values[species] = a
            values['x'] = x
            profiles[base][level] = values
        medium = profiles[base]['medium']
        for level in ['coarse', 'fine', 'half_step']:
            assert np.allclose(profiles[base][level]['x'], medium['x'], atol=1e-12, rtol=0)
        for width in [1, 2, 4]:
            for offset in range(width):
                a = {k: aggregate(v['element_I'], width, offset) for k,v in profiles[base].items()}
                for k in a:
                    assert np.isclose(a[k].sum(), profiles[base][k]['element_I'].sum(), rtol=2e-15, atol=0)
                fm = relative_l1(a['fine'], a['medium'])
                cm = relative_l1(a['coarse'], a['medium'])
                bins.append({'base_case': base, 'width_cm': width, 'offset_native_bins': offset,
                             'coarse_medium_L1_pct_medium_denominator': cm,
                             'fine_medium_L1_pct': fm,
                             'successive_difference_ratio': fm/cm if cm else None})
    shapes = []
    for base in BASES:
        p = profiles[base]
        shapes.append({'base_case':base, **shape_metrics(p['fine']['element_I'],p['medium']['element_I'],p['medium']['x'])})
    thermo = []
    for level in ['coarse', 'medium', 'fine', 'half_step']:
        shares = []
        for base in BASES[2:]:
            a = profiles[base][level]
            shares.append(float(300*a['BiI3_g'].sum()/a['element_I'].sum()))
        thermo.append({'resolution':level, 'baseline_BiI3_iodine_pct':shares[0],
                       'shifted_BiI3_iodine_pct':shares[1], 'separation_percentage_points':shares[0]-shares[1]})
    hashes[Path(__file__).relative_to(REPO).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return {'scope':'Saved-profile reassessment only; bin aggregation is not mesh refinement and tested-grid ranges are not continuum uncertainty bounds.',
            'input_sha256':hashes, 'aggregation':bins, 'fine_medium_shapes':shapes,
            'thermodynamic_resolution':thermo,
            'minimum_tested_thermodynamic_separation_pp':min(r['separation_percentage_points'] for r in thermo)}


def main():
    report = evaluate()
    out = HERE / 'generated/resolution-audit.json'
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    rows = []
    for base,label in zip(BASES,LABELS):
        row = [label]
        for width in [1,2,4]:
            values = [r['fine_medium_L1_pct'] for r in report['aggregation'] if r['base_case']==base and r['width_cm']==width]
            row.append(f'{min(values):.2f}' if width==1 else f'{min(values):.2f}--{max(values):.2f}')
        shape = next(r for r in report['fine_medium_shapes'] if r['base_case']==base)
        row.append(f"{shape['binned_W1_cm']:.4f}")
        rows.append(' & '.join(row)+r' \\')
    table = '\\begin{tabular}{lrrrr}\n\\toprule\nCase & 1 cm & 2 cm & 4 cm & $W_1$ (cm) \\\\\n\\midrule\n'
    (HERE/'generated/resolution-audit-table.tex').write_text(table+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n')
    print(json.dumps({'thermodynamic_resolution':report['thermodynamic_resolution'],
                      'fine_medium_shapes':report['fine_medium_shapes']},indent=2))


if __name__ == '__main__':
    main()
