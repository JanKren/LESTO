#!/usr/bin/env python3
"""Production-scheme (576 x 72, Gauss limitedLinear 1) sensitivities against their own production
baselines, with the follow-up metrics of report_followup.py, next to the upwind 18,432-cell values.

    python3 analyse.py      # reads the tracked profiles, writes comparison.csv and analysis.json

Variant profiles: profiles/<case>/; base profiles: ../production_576x72_limitedLinear/profiles/<base>/;
upwind values: ../followup_results/analysis.json and ../full_analysis/case_metrics.csv.
"""
import csv, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
BASES = RUN/'production_576x72_limitedLinear/profiles'
FILES = {'I': 'element_I_native.dat', 'PbI2': 'PbI2_g_native.dat', 'BiI3': 'BiI3_g_native.dat'}

UP = json.loads((RUN/'followup_results/analysis.json').read_text())
up_cmp = {(r['case'], r['species']): r for r in UP['comparisons']}
up_met = {r['case']: r for r in UP['case_metrics']}
up_base = {r['case']: float(r['BiI3_iodide_I_pct'] or 0) for r in csv.DictReader((RUN/'full_analysis/case_metrics.csv').open())}
ledger = json.loads((HERE/'ledger.json').read_text())


def profile(case_dir, sp):
    """Bin centres [m] and wall amount per bin [mol] of the 60 s native profile."""
    a = np.loadtxt(case_dir/FILES[sp], comments='#', ndmin=2)
    return a[:, 2], a[:, 4]*(a[:, 1] - a[:, 0])


def l1(a, b): return float(np.abs(a - b).sum()/np.abs(b).sum())
def centroid(x, a): return float((x*a).sum()/a.sum())


def share(case_dir):
    """BiI3 share 3 Bi/(2 Pb + 3 Bi) of wall iodide iodine [%]."""
    pb = profile(case_dir, 'PbI2')[1].sum(); bi = profile(case_dir, 'BiI3')[1].sum()
    return 100*3*bi/(2*pb + 3*bi) if 2*pb + 3*bi > 0 else None


rows = []
for case in sorted((HERE/'profiles').iterdir()):
    base, var = case.name.split('__'); b = BASES/base
    row = dict(case=case.name, base=base, variant=var, last_time=ledger[case.name]['end_time'],
               exec_h=ledger[case.name]['execution_h'], ended=ledger[case.name]['ended'])
    for sp in FILES:
        x, a = profile(case, sp); xb, ab = profile(b, sp); assert np.allclose(x, xb)
        if ab.sum() == 0:
            continue
        row[f'{sp}_L1_pct'] = 100*l1(a, ab)
        row[f'{sp}_total_change_pct'] = 100*(a.sum() - ab.sum())/ab.sum()
        row[f'{sp}_centroid_shift_cm'] = 100*(centroid(x, a) - centroid(xb, ab))
    row['BiI3_share_pct'] = share(case); row['BiI3_share_base_pct'] = share(b)
    row['BiI3_share_change_pp'] = row['BiI3_share_pct'] - row['BiI3_share_base_pct']
    u = up_cmp[(case.name, 'I')]
    row['upwind_I_L1_pct'] = 100*u['profile_L1_relative_to_base']
    row['upwind_I_total_change_pct'] = 100*u['deposit_relative_change']
    row['upwind_I_centroid_shift_cm'] = u['centroid_shift_cm']
    row['upwind_BiI3_share_pct'] = up_met[case.name]['BiI3_iodide_I_pct']
    row['upwind_BiI3_share_change_pp'] = row['upwind_BiI3_share_pct'] - up_base[base]
    rows.append(row)

(HERE/'analysis.json').write_text(json.dumps(rows, indent=1) + '\n')
keys = list(dict.fromkeys(k for r in rows for k in r))
with (HERE/'comparison.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, keys, lineterminator='\n'); w.writeheader(); w.writerows(rows)
fmt = lambda v, p='.2f': '' if v is None else format(v, p)
print(f"{'case':46s} | {'L1 LL':>6} {'L1 up':>6} | {'tot LL':>7} {'tot up':>7} | {'cen LL':>7} {'cen up':>7} | {'dBiI3 LL':>8} {'up':>8}")
for r in rows:
    print(f"{r['case']:46s} | {fmt(r.get('I_L1_pct')):>6} {fmt(r['upwind_I_L1_pct']):>6} | "
          f"{fmt(r.get('I_total_change_pct'), '+.3f'):>7} {fmt(r['upwind_I_total_change_pct'], '+.3f'):>7} | "
          f"{fmt(r.get('I_centroid_shift_cm'), '+.4f'):>7} {fmt(r['upwind_I_centroid_shift_cm'], '+.4f'):>7} | "
          f"{fmt(r['BiI3_share_change_pp'], '+.3f'):>8} {fmt(r['upwind_BiI3_share_change_pp'], '+.3f'):>8}")
