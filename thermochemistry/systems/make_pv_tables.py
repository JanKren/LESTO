#!/usr/bin/env python3
"""
make_pv_tables.py -- phase-resolved vapour-pressure tables of the M10
Pb-Bi-I-He dataset for the Hertz-Knudsen-Schrage model of rhoFixedFlowFoam,
in the format of thermochemistry/systems/data/pv_PbI2_phases.csv
(vapourPressureTable.H): '#' comment lines, the '#' line right before the
first data row names the columns (comma separated), column 0 is T [K], each
further column is log10(p/Pa) of the gas species over one pure condensed
phase over the whole range (crystal extrapolated above its melting point,
liquid below it); the solver takes p_eq = the minimum over the columns, so
the kink at the melting point (a node) is exact.  Nodes every 5 K from 270
to 1250 K plus the crossing of the columns (rounded to 1e-6 K unless it
falls on a 5 K node).

Tables (gas over condensate, reaction nu*cond = gas):
  pv_BiI3_phases.csv   BiI3(g) over BiI3(cr), BiI3(l)       T_fus 680.85 K (TKV)
  pv_Bi_phases.csv     Bi(g)   over Bi(cr), Bi(l)            T_m of Barin, 544.51 K
  pv_Bi2_phases.csv    Bi2(g)  over 2 Bi(cr), 2 Bi(l)        same melting point
  pv_BiI_phases.csv    BiI(g)  over BiI(cr)                  one phase
  pv_Pb_phases.csv     Pb(g)   over Pb(cr), Pb(l)            (extra; NASA/Gurvich)
and, as a regression test only, pv_PbI2_phases.csv, which must equal the
repository's file byte for byte (--regress FILE).

The data are those of the PbBiIHe GEMS3K system (make_gems3k_input.pbbii_system:
NASA-9 records of thermo.inp for Pb-I-He, data/bi_i_m10-thermofun.json for
Bi-I), evaluated with the same functions as its DCH grid; p0 = 1 bar,
R = 8.314462618 J/(mol K).

Usage
  make_pv_tables.py [--nasa thermo.inp] [--tf bi_i_m10-thermofun.json]
                    [--out DIR] [--check] [--regress pv_PbI2_phases.csv
                    [--regress-dir DIR]]
                    [--reference-dir DIR --reference-step K]

  --check      per table: the crossing of the columns and the largest error
               of the table interpolated as the solver does (log10 p linear
               in 1/T per phase, then the minimum) against the direct
               evaluation on a 0.1 K grid, 270-1250 K (limit 2e-3 decades)
  --reference-dir  also write pv_<gas>_reference.csv: T and log10(p_eq/Pa)
               (minimum over the phases) every --reference-step K (0.25)

Licences: the Bi tables derive from HERACLES-TDB (= Barin 1995; provenance
to be confirmed with PSI, internal until then) and, for BiI3(l), TKV
literature values; the Pb and PbI2 tables from NASA thermo.inp (Apache-2.0).
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_gems3k_input import R, P0, species_props, pbbii_system  # noqa: E402

T_MIN, T_MAX, STEP = 270.0, 1250.0, 5.0

HER = ('HERACLES-TDB v0.2 (= Barin 1995) via data/bi_i_m10-thermofun.json; '
       'internal until PSI confirms the provenance')
TABLES = [
    # file, gas, [(column, DC)], nu, crossing bracket, source lines
    ('pv_BiI3_phases.csv', 'BiI3(g)', [('BiI3(cr)', 'BiI3(cr)'), ('BiI3(l)', 'BiI3(l)')], 1.0,
     (600.0, 750.0),
     ['BiI3(g), BiI3(cr): ' + HER + ';',
      'BiI3(cr) cut at T_fus = 680.85 K (TKV); BiI3(l) constructed: BiI3(cr) at T_fus + '
      'dH_fus 31.798 kJ/mol (TKV) + Cp(l) 150.624 J/(mol K) (Barin)']),
    ('pv_Bi_phases.csv', 'Bi(g)', [('Bi(cr)', 'Bi(cr)'), ('Bi(l)', 'Bi(l)')], 1.0, (500.0, 600.0),
     ['Bi(g), Bi(cr), Bi(l): ' + HER + ';',
      'the vapour also holds Bi2(g) (pv_Bi2_phases.csv); Bi4(g) has no open record']),
    ('pv_Bi2_phases.csv', 'Bi2(g)', [('Bi(cr)', 'Bi(cr)'), ('Bi(l)', 'Bi(l)')], 2.0,
     (500.0, 600.0),
     ['Bi2(g) over 2 Bi(cr) and 2 Bi(l) (2 Bi(cond) = Bi2(g)): ' + HER + ';',
      'the monomer is in pv_Bi_phases.csv; Bi4(g) has no open record']),
    ('pv_BiI_phases.csv', 'BiI(g)', [('BiI(cr)', 'BiI(cr)')], 1.0, None,
     ['BiI(g), BiI(cr): ' + HER + ';',
      'one phase (no open BiI(l)); BiI(cr) record 273-600 K, extrapolated outside']),
    ('pv_Pb_phases.csv', 'Pb(g)', [('Pb(cr)', 'Pb(cr)'), ('Pb(l)', 'Pb(l)')], 1.0, (550.0, 650.0),
     ['Gurvich et al. (1991) through the NASA-9 records Pb, Pb(cr), Pb(L) of NASA Glenn '
      'thermo.inp;', 'Apache-2.0 (NASA CEA); Pb2(g) is not in thermo.inp']),
]


class Data:
    def __init__(self, nasa, tf):
        sysd = pbbii_system(nasa, tf)
        self.sp = {s['name']: s for ph in sysd['phases'] for s in ph['species']}

    def G(self, name, T):
        return species_props(self.sp[name], T)[3]


def pv_phase(data, gas, dc, nu, T):
    """log10(p/Pa) of `gas` over the pure phase `dc` (nu dc = gas) at T [K]."""
    return -(data.G(gas, T) - nu * data.G(dc, T)) / (R * T * math.log(10.0)) + math.log10(P0)


def crossing(data, gas, cols, nu, lo, hi):
    """T [K] where the two columns cross (bisection to 1e-12 K)."""
    def f(T):
        return pv_phase(data, gas, cols[0][1], nu, T) - pv_phase(data, gas, cols[1][1], nu, T)
    flo = f(lo)
    if flo * f(hi) > 0:
        raise ValueError(f'the phase curves of {gas} do not cross in [{lo}, {hi}] K')
    while hi - lo > 1e-12:
        mid = 0.5 * (lo + hi)
        if f(mid) * flo > 0:
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


def nodes(t_cross):
    n = int(round((T_MAX - T_MIN) / STEP))
    grid = [T_MIN + STEP * k for k in range(n + 1)]
    if t_cross is not None:
        t_node = round(t_cross, 6)
        if all(abs(t - t_node) > 1e-9 for t in grid):
            grid.append(t_node)
    return sorted(grid)


def write_table(data, path, gas, cols, nu, t_cross, source_lines, title=None):
    grid = nodes(t_cross)
    with open(path, 'w') as f:
        f.write(title or ('# %s equilibrium vapour pressure over its condensed phases, '
                          'log10(p/Pa) of %s\n' % (gas.replace('(g)', ''), gas)))
        for ln in source_lines:
            f.write('# %s\n' % ln)
        f.write(f'# p0 = 1 bar, R = {R} J/(mol K)\n')
        if len(cols) > 1:
            f.write('# each column is one phase over the whole range (crystal extrapolated '
                    'above T_m, liquid below);\n')
            f.write('# p_eq = the minimum over the phases.  Melting point (crossing of the '
                    f'curves) T_m = {t_cross:.6f} K, a node.\n')
        else:
            f.write('# one column: p_eq over the only condensed phase of the data\n')
        f.write('# written by make_pv_tables.py (M10 design, scratch)\n')
        f.write('# T_K,' + ','.join(c for c, _ in cols) + '\n')
        for T in grid:
            vals = [pv_phase(data, gas, dc, nu, T) for _, dc in cols]
            f.write(f'{T:.6f},' + ','.join(f'{v:.12f}' for v in vals) + '\n')
    return grid


def interpolate(grid, columns, T):
    """log10 p_eq as vapourPressureTable.H: per phase linear in 1/T (end
    intervals extrapolated), then the minimum over the phases."""
    n = len(grid)
    i = 0
    if T >= grid[-1]:
        i = n - 2
    elif T > grid[0]:
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if grid[mid] > T:
                hi = mid
            else:
                lo = mid
        i = lo
    w = (1.0 / T - 1.0 / grid[i]) / (1.0 / grid[i + 1] - 1.0 / grid[i])
    return min((1 - w) * c[i] + w * c[i + 1] for c in columns)


def check(data, gas, cols, nu, grid):
    columns = [[float(f'{pv_phase(data, gas, dc, nu, T):.12f}') for T in grid] for _, dc in cols]
    worst, where, k = 0.0, None, 0
    while True:
        T = T_MIN + 0.1 * k
        if T > T_MAX + 1e-9:
            break
        exact = min(pv_phase(data, gas, dc, nu, T) for _, dc in cols)
        err = abs(interpolate(grid, columns, T) - exact)
        if err > worst:
            worst, where = err, T
        k += 1
    return worst, where


def regress_pbi2(data, ref_path, out_dir):
    """pv_PbI2_phases.csv as the repository writes it; compare byte for byte."""
    cols = [('PbI2(cr)', 'PbI2(cr)'), ('PbI2(l)', 'PbI2(l)')]
    t_m = crossing(data, 'PbI2(g)', cols, 1.0, 600.0, 750.0)
    path = os.path.join(out_dir, 'pv_PbI2_phases.regress.csv')
    title = ('# PbI2 equilibrium vapour pressure over its condensed phases, log10(p/Pa) of '
             'PbI2(g)\n')
    grid = nodes(t_m)
    with open(path, 'w') as f:
        f.write(title)
        f.write('# Gurvich et al. (1991) through the NASA-9 records PbI2, PbI2(cr), PbI2(L) of '
                'NASA Glenn thermo.inp\n')
        f.write(f'# (file nasa_thermo.inp); p0 = 1 bar, R = {R} J/(mol K); Apache-2.0 '
                '(NASA CEA)\n')
        f.write('# each column is one phase over the whole range (crystal extrapolated '
                'above T_m, liquid below);\n')
        f.write('# p_eq = the minimum over the phases.  Melting point (crossing of the '
                f'curves) T_m = {t_m:.6f} K, a node.\n')
        f.write('# written by thermochemistry/systems/make_pv_table.py\n')
        f.write('# T_K,PbI2(cr),PbI2(l)\n')
        for T in grid:
            vals = [pv_phase(data, 'PbI2(g)', dc, 1.0, T) for _, dc in cols]
            f.write(f'{T:.6f},' + ','.join(f'{v:.12f}' for v in vals) + '\n')
    same = open(path, 'rb').read() == open(ref_path, 'rb').read()
    return same, path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--nasa', required=True)
    ap.add_argument('--tf', required=True)
    ap.add_argument('--out', required=True, help='output directory (external for restricted data)')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--regress')
    ap.add_argument('--regress-dir', help='where the regression copy goes (default --out)')
    ap.add_argument('--reference-dir')
    ap.add_argument('--reference-step', type=float, default=0.25)
    a = ap.parse_args()
    if not math.isfinite(a.reference_step) or a.reference_step <= 0:
        ap.error("--reference-step must be finite and positive")
    os.makedirs(a.out, exist_ok=True)
    data = Data(a.nasa, a.tf)
    ok = True
    for fn, gas, cols, nu, bracket, src in TABLES:
        t_c = crossing(data, gas, cols, nu, *bracket) if bracket else None
        path = os.path.join(a.out, fn)
        grid = write_table(data, path, gas, cols, nu, t_c, src)
        print('Wrote %s: %d nodes%s' % (path, len(grid),
                                        (', T_m = %.6f K' % t_c) if t_c else ''))
        if a.check:
            worst, where = check(data, gas, cols, nu, grid)
            good = worst < 2e-3
            ok &= good
            print('   check: max |table - direct| = %.3e decades at %.1f K (0.1 K grid, '
                  '270-1250 K) %s' % (worst, where, 'ok' if good else 'FAILED'))
        if a.reference_dir:
            os.makedirs(a.reference_dir, exist_ok=True)
            rp = os.path.join(a.reference_dir, fn.replace('_phases.csv', '_reference.csv'))
            with open(rp, 'w') as f:
                f.write('# log10(p_eq/Pa) of %s, minimum over %s (M10 dataset)\n'
                        % (gas, ', '.join(c for c, _ in cols)))
                f.write('# written by make_pv_tables.py --reference-dir\n')
                f.write('# T_K,log10p\n')
                k = 0
                while True:
                    T = T_MIN + a.reference_step * k
                    if T > T_MAX + 1e-9:
                        break
                    v = min(pv_phase(data, gas, dc, nu, T) for _, dc in cols)
                    f.write(f'{T!r},{v!r}\n')
                    k += 1
    if a.regress:
        rdir = a.regress_dir or a.out
        os.makedirs(rdir, exist_ok=True)
        same, path = regress_pbi2(data, a.regress, rdir)
        ok &= same
        print('regression: %s %s the repository file %s' % (path, 'EQUALS' if same else
                                                              'DIFFERS FROM', a.regress))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
