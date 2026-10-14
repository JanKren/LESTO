#!/usr/bin/env python3
"""Compare a finished Lobresco case with the digitised experiment (Fig. 13/14).

Usage: analyse.py <case> [<case> ...]  -> <case>/comparison.json, comparison.png
Profiles: postProcessing/axialProfiles/<time>/PbI2_g_gamma.dat (wedge, mol/m) x 360/5.
T_dep: paper rule (first bin from the inlet whose value exceeds 1 % of the peak,
linear interpolation between the bracketing 1 cm bins), wall T from the digitised curve.
"""
import csv, json, sys, re
from pathlib import Path
import numpy as np
G = Path(__file__).resolve().parent
DIGITIZED = G.parents[1] / 'thermochemistry/benchmarks/lobresco2026/digitized'
SCALE = 360 / 5


def read_profile(path):
    meta, rows = {}, []
    for line in open(path):
        if line.startswith('#'):
            m = re.match(r'# Tdep (\w+) status (\w+)(?: T (\S+) x (\S+))?', line)
            if m:
                meta['Tdep_' + m.group(1)] = dict(status=m.group(2), T=float(m.group(3)) if m.group(3) else None,
                                                   x=float(m.group(4)) if m.group(4) else None)
            if line.startswith('# columns'):
                meta['columns'] = line.split()[2:]
            if line.startswith('# inventory'):
                p = line.split(); meta['inventory'] = dict(gas=float(p[3]), wall=float(p[5]), sample=float(p[7]))
            continue
        if line.strip():
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    return meta, {c: a[:, i] for i, c in enumerate(meta['columns'])}


def onset(x, y, Tw, fraction=0.01, from_peak=False):
    pk = y.max()
    if not pk > 0:
        return None
    thr = fraction * pk
    if from_peak:
        k = int(np.argmax(y))
        while k > 0 and y[k - 1] >= thr:
            k -= 1
    else:
        k = int(np.flatnonzero(y >= thr)[0])
    if k == 0:
        return None
    xo = x[k - 1] + (thr - y[k - 1]) / (y[k] - y[k - 1]) * (x[k] - x[k - 1])
    return dict(x_m=float(xo), T_K=float(np.interp(xo, *Tw)), x_peak_m=float(x[np.argmax(y)]),
                peak_mol_per_m=float(pk), T_peak_K=float(np.interp(x[np.argmax(y)], *Tw)))


def main():
    out_all = {}
    for case in map(Path, sys.argv[1:]):
        info = json.loads((case / 'case.json').read_text())
        name = info['case']
        flow = info['table1']['flow_mL_min']
        wt = list(csv.DictReader(open(G / 'inputs' / f'wallT_{flow}mLmin.csv')))
        Tw = (np.array([float(r['x_m']) for r in wt]), np.array([float(r['T_K']) for r in wt]))
        times = sorted((p for p in (case / 'postProcessing/axialProfiles').iterdir()), key=lambda p: float(p.name))
        last = times[-1]
        meta, prof = read_profile(last / 'PbI2_g_gamma.dat')
        x = prof['x']
        obs = prof['observable'] * SCALE
        wall = prof['wall'] * SCALE
        ex = list(csv.DictReader(open(DIGITIZED / f'{name}_experiment.csv')))
        xe = np.array([float(r['x_bin_m']) for r in ex]); ye = np.array([float(r['n_mol_per_m']) for r in ex])
        n0 = info['table1']['moles']
        # experiment on the same 1 cm centres (missing scan points outside the plotted range = 0)
        ye_full = np.array([ye[np.argmin(abs(xe - xi))] if np.min(abs(xe - xi)) < 1e-6 else 0.0 for xi in x])
        res = dict(case=name, profile_time=float(last.name), solver_onsets=meta,
                   model_wallPlusGas_fromInlet=onset(x, obs, Tw), model_wallPlusGas_fromPeak=onset(x, obs, Tw, from_peak=True),
                   model_wall_fromInlet=onset(x, wall, Tw), model_wall_fromPeak=onset(x, wall, Tw, from_peak=True),
                   experiment=onset(xe, ye, Tw), table2_K=dict(T106_6=814, T106_60=723, T106_300=667, T5_60=847, T540_60=700)[name],
                   model_integral_over_n0=float(obs.sum() * 0.01 / n0),
                   L1_obs_minus_exp_over_n0=float(np.abs(obs - ye_full).sum() * 0.01 / n0),
                   inventory_full_column_mol={k: v * SCALE for k, v in meta['inventory'].items()})
        if res['model_wallPlusGas_fromInlet'] and res['experiment']:
            res['dT_dep_fromInlet_K'] = res['model_wallPlusGas_fromInlet']['T_K'] - res['table2_K']
            res['dx_peak_m'] = res['model_wallPlusGas_fromInlet']['x_peak_m'] - res['experiment']['x_peak_m']
        bal = sorted((case / 'postProcessing/phaseChangeBalance').glob('*/PbI2_g.dat'))
        if bal:
            lines = [l for l in open(bal[-1]) if not l.startswith('#') and l.strip()]
            hdr = [l for l in open(bal[-1]) if l.startswith('#')][-1].lstrip('#').split()
            res['ledger_last'] = dict(zip(hdr, lines[-1].split()))
        (case / 'comparison.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
        out_all[str(case)] = res
        try:
            import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(xe, ye, 'bo-.', ms=4, label='experiment (digitised)')
            ax.step(x, obs, where='mid', color='k', label='HKS, wall+gas, 1 cm bins')
            ax.step(x, wall, where='mid', color='g', ls=':', label='HKS, wall only')
            ax.set_xlabel('x [m]'); ax.set_ylabel('PbI2 [mol/m]'); ax.set_title(f'{name} at t = {last.name} s')
            ax2 = ax.twinx(); ax2.plot(*Tw, 'r--', lw=.8); ax2.set_ylabel('T_wall [K]')
            ax.legend(fontsize=8); fig.savefig(case / 'comparison.png', dpi=130, bbox_inches='tight'); plt.close(fig)
        except Exception as e:  # plotting is optional
            print('plot skipped', e)
        print(name, json.dumps({k: res.get(k) for k in ['dT_dep_fromInlet_K', 'dx_peak_m', 'model_integral_over_n0', 'L1_obs_minus_exp_over_n0']}))
    return out_all


if __name__ == '__main__':
    main()
