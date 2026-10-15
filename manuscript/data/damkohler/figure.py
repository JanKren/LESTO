#!/usr/bin/env python3
"""Figure of the Damkoehler sweep (manuscript Fig. damkohler_sweep) from sweep_summary.json.

    python3 manuscript/data/damkohler/figure.py   # writes manuscript/figures/damkohler_sweep.{pdf,png}
"""
import json, os, numpy as np, matplotlib
HERE = os.path.dirname(os.path.abspath(__file__))
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = json.load(open(os.path.join(HERE, 'sweep_summary.json')))['sweep']; runs = d['runs']; T = np.array(d['T_wall_cell']); x = np.array(d['x_mm'])
y = np.array(d['y_mm']); h = 2.4
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.labelsize': 8, 'legend.fontsize': 8,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.grid': True, 'grid.alpha': .3, 'lines.linewidth': 1.2})
fig, ax = plt.subplots(2, 2, figsize=(5.4, 6.4), constrained_layout=True)
bands = [(1.4e3, 2.9e3, r'$\alpha$=1'), (76, 150, r'$\alpha$=0.1'), (7, 15, r'$\alpha$=0.01')]
def shade(a, label=False):
    for lo, hi, t in bands:
        a.axvspan(lo, hi, color='0.85', zorder=0, lw=0)
        if label: a.text(np.sqrt(lo*hi), 0.02, t, rotation=90, ha='center', va='bottom', fontsize=8, color='0.35', transform=a.get_xaxis_transform())
# (a) normalised wall flux vs local Da at three stations
A = ax[0, 0]; shade(A)
Da = np.logspace(-3.3, 3.6, 200); Sh = np.pi**2/4
A.plot(Da, Da, 'k--', lw=.9, label=r'kinetic limit $\hat J$=Da')
A.axhline(Sh, color='k', ls=':', lw=.9, label=r'diffusion limit $\pi^2/4$')
A.plot(Da, 1/(1/Da+1/Sh), color='0.5', lw=.8, label='series resistances')
for i, mk, col in [(72, 's', '#9ecae1'), (84, 'o', '#3182bd'), (94, '^', '#08306b')]:
    dl = [r['t1.5']['Da_loc'][i] for r in runs]; jh = [r['t1.5']['Jhat'][i] for r in runs]
    A.plot(dl, jh, mk, ms=4, mfc=col, mec='k', mew=.4, ls='none', label=f'x = {x[i]:.1f} mm, {T[i]:.0f} K')
A.set(xscale='log', yscale='log', xlabel='wall Damköhler number Da (local $T_w$)',
      ylabel=r'$\hat J = J\,h/[D(\bar c-c_{eq})]$', ylim=(5e-4, 20), xlim=(5e-4, 5e3))
A.legend(loc='upper left', bbox_to_anchor=(0, -0.22), ncol=1, fontsize=8, frameon=False); A.set_title('(a) normalised wall flux', loc='left', fontsize=8)
# (b) wall-normal profiles
B = ax[0, 1]
sel = [(1e-5, '#c6dbef'), (1e-3, '#6baed6'), (1e-2, '#2171b5'), (1.0, '#08306b')]
for a, col in sel:
    r = min(runs, key=lambda r: abs(np.log(r['alpha']/a)))
    for i, ls in [(72, '--'), (84, '-')]:
        p = np.array(r['t1.5']['profiles'][str(i)]); yy = y[:len(p)]
        half = yy <= h
        B.plot(yy[half]/h, p[half], ls=ls, color=col, label=f"Da={r['t1.5']['Da_loc'][84]:.2g}" if ls == '-' else None)
B.plot([], [], 'k-', label=f'x = {x[84]:.1f} mm'); B.plot([], [], 'k--', label=f'x = {x[72]:.1f} mm')
B.set(xlabel='distance from wall $y/h$', ylabel=r'$c/\bar c$', xlim=(0, 1), ylim=(0, 1.65))
B.legend(fontsize=8, ncol=2, loc='upper left', bbox_to_anchor=(0, -0.22), frameon=False); B.set_title('(b) profiles at t = 1.5 s', loc='left', fontsize=8)
# (c) deposited fraction vs Da(561 K)
C = ax[1, 0]; shade(C)
Dr = np.array([r['Da_at_561K'] for r in runs]); f = np.array([r['deposited_fraction'] for r in runs])
# plug-flow kinetic limit: capture exponent Lambda = int_{x_on}^{L} G_k rho /(G h) dx, rho = p M_He/(R T)
Ru, M, MHe, p, G, hm = 8.314462618, .46100894, 4.0026e-3, 101325., .01743, 2.4e-3
on = x >= runs[-1]['x_onset_mm']
lam1 = np.sum(np.sqrt(Ru*T[on]/(2*np.pi*M))*(p*MHe/(Ru*T[on]))/(G*hm))*1e-3   # per unit 2a/(2-a)
Dk = np.logspace(-3.3, 1, 100); C.plot(Dk, 1-np.exp(-lam1*Dk/runs[-1]['Da_at_561K']*2), 'k--', lw=.9, label='plug-flow kinetic limit')
C.plot(Dr, f, 'o-', ms=4, color='#08519c', mec='k', mew=.4, label='deposited at 3 s')
C.set(xscale='log', yscale='log', xlabel='Da at onset temperature (561 K)', ylabel='deposited fraction of release', ylim=(1e-3, 2), xlim=(5e-4, 5e3))
C.legend(fontsize=8, loc='center right'); C.set_title('(c) captured fraction', loc='left', fontsize=8)
# (d) positions
Dd = ax[1, 1]; shade(Dd, True)
Dd.plot(Dr, [r['x_onset_mm'] for r in runs], 'v-', ms=4, color='0.4', label='1 % onset')
Dd.plot(Dr, [r['x_peak_mm'] for r in runs], 's-', ms=4, color='#6baed6', mec='k', mew=.4, label='peak bin')
Dd.plot(Dr, [r['centroid_mm'] for r in runs], 'o-', ms=4, color='#08306b', mec='k', mew=.4, label='centroid')
Dd.set(xscale='log', xlabel='Da at onset temperature (561 K)', ylabel='axial position (mm)', ylim=(52, 102), xlim=(5e-4, 5e3))
Dd.legend(fontsize=8, loc='upper right'); Dd.set_title('(d) deposit position', loc='left', fontsize=8)
for e in ['pdf', 'png']:
    fig.savefig(os.path.join(HERE, f'../../figures/damkohler_sweep.{e}'), dpi=300,
                metadata={'CreationDate': None} if e == 'pdf' else None)
print('capture exponent per unit Da/2:', lam1)
