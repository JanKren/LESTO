#!/usr/bin/env python3
"""Wall-temperature inputs per flow rate from the digitised Lobresco (2026) curves."""
import csv, json, hashlib
from pathlib import Path
import numpy as np
G = Path(__file__).resolve().parent
D = G.parents[1] / 'thermochemistry/benchmarks/lobresco2026/digitized'
def load(c):
    r = list(csv.DictReader(open(D / f'{c}_wall_temperature.csv')))
    return np.array([float(q['x_m']) for q in r]), np.array([float(q['T_K']) for q in r])
grid = np.round(np.arange(0.0, 0.690001, 0.001), 4)
prof = {}
x6, t6 = load('T106_6'); x3, t3 = load('T106_300')
prof['106'] = 0.5 * (np.interp(grid, x6, t6) + np.interp(grid, x3, t3))
for flow, c in [('5', 'T5_60'), ('540', 'T540_60')]:
    x, t = load(c); prof[flow] = np.interp(grid, x, t)
meta = {}
for k, T in prof.items():
    p = G / 'inputs' / f'wallT_{k}mLmin.csv'
    with open(p, 'w') as f:
        f.write('x_m,T_K\n'); f.writelines(f'{x:.3f},{t:.2f}\n' for x, t in zip(grid, T))
    meta[k] = dict(file=str(p.relative_to(G)), sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                   T_max_K=float(T.max()), x_T_max_m=float(grid[np.argmax(T)]),
                   T_at_0_K=float(T[0]), T_at_0p69_K=float(T[-1]))
meta['106']['source'] = 'mean of the T106_6 and T106_300 digitisations (Fig. 14 i, ii); T106_60 (Fig. 13) agrees within 3.3 K over 0.2-0.6 m but its curve is hidden by the legend at 0.105-0.135 m'
meta['5']['source'] = 'T5_60 digitisation (Fig. 14 iii)'
meta['540']['source'] = 'T540_60 digitisation (Fig. 14 iv)'
meta['note'] = 'Values below x = 0.001 m (first digitised point) are the value at 0.001 m, interpolated onto a 1 mm grid; held constant by the carrier tool outside [0, 0.69] m.'
(G / 'inputs' / 'wallT_provenance.json').write_text(json.dumps(meta, indent=1))
print(json.dumps(meta, indent=1))
