#!/usr/bin/env python3
"""
make_pv_table.py -- vapour-pressure tables of PbI2 for rhoFixedFlowFoam.

PURPOSE

Writes the two p_v tables that the Hertz-Knudsen-Schrage model of
rhoFixedFlowFoam reads (constant/thermochemistryProperties, pairs { PbI2_g {
vapourPressure { file ...; } } }; applications/rhoFixedFlowFoam/
vapourPressureTable.H, plan section 2, E9):

  data/pv_PbI2_phases.csv   phase-resolved Gurvich (1991) data through the
                            NASA-9 records of NASA Glenn's thermo.inp: one
                            column per condensed phase, log10(p/Pa) of
                            PbI2(g) over PbI2(cr) and over PbI2(l), every
                            5 K from 270 to 1250 K plus the node T_m where
                            the two columns cross (683 K, the melting point
                            of the data).  The solver takes p_eq = the
                            minimum over the listed phases (the stable
                            condensate), so the kink at T_m is exact.
  data/pv_PbI2_TFlows.csv   the GEMS_VAPOR_PRESSURE column [bar] of the
                            T-Flows reference implementation (Types.f90 of
                            Tests/Laminar/Scalar_Transport_Gpu/LESTO_HKS,
                            branch advanced_scalar_model, pinned to the
                            revision TFLOWS_REVISION below, which the csv
                            header records) at its 95 nodes 283.15 + 10 k
                            K, copied digit for digit, for the code-to-code
                            comparison (interpolation linearT, outOfRange
                            hold, units bar in the solver).

Thermodynamics.  For a phase phi of the condensate, the reaction
PbI2(phi) -> PbI2(g) has the standard Gibbs energy dG(T) = G_g(T) - G_phi(T)
(NASA-9 polynomials, standard pressure p0 = 1 bar), and the equilibrium
partial pressure of the monomer over the pure condensate is

  p_phi(T) = p0 exp(-dG/(R T)),   log10(p_phi/Pa) = -dG/(R T ln 10) + 5.

Each column is the phase's own curve over the whole range: the crystal
polynomial (fitted 298.15-683 K) is extrapolated above T_m (metastable
crystal), the liquid one (683-6000 K, constant Cp) below T_m (supercooled
liquid), and the gas polynomial (300-1000-6000 K) below 300 K.  The stable
phase has the lower p_v, so min(p_cr, p_l) is the equilibrium over the
stable condensate, and the crossing of the two columns is the melting point
of the data.  The dimer Pb2I4(g) is not in thermo.inp and not included.

The script reads a LOCAL thermo.inp (it never downloads anything; the
default is the copy in ~/opt/gems-data, fetch_nasa_thermo.sh gets the
pinned public file) and the T-Flows source through a local read-only
'git show' of a fixed commit (not of the branch, whose tip moves with every
fetch: a changed Types.f90 must not change the regenerated table).

Usage

  make_pv_table.py [--nasa thermo.inp] [--tflows Types.f90 | --tflows-git
                   <repository>] [--out data] [--check]
                   [--reference FILE --reference-step K] [--sources]

  --nasa         NASA Glenn thermo.inp (default
                 ~/opt/gems-data/research-2026-09-23/format/nasa_thermo.inp)
  --tflows       a copy of Types.f90; otherwise it is read with
                 git -C <repository> show <spec> (default repository
                 ~/src/T-Flows, spec <TFLOWS_REVISION>:Tests/Laminar/
                 Scalar_Transport_Gpu/LESTO_HKS/User_Mod/Types.f90)
  --out          directory of the two csv files (default: data next to this
                 script)
  --check        print the melting point (the crossing of the NASA-9
                 curves) and the largest error of the 5 K table,
                 interpolated as the solver does (log10 p linear in 1/T per
                 phase, then the minimum), against NASA-9 on a 0.1 K grid
  --reference    also write FILE: T and log10(p_eq/Pa) of NASA-9 (the minimum
                 over the phases) every --reference-step K (default 0.25)
                 from 270 to 1250 K, for the solver's table tests
                 (tests/testVapourPressureTable.C)
  --sources      only check that both sources can be read (thermo.inp, and
                 the pinned T-Flows revision or the --tflows file); prints
                 what is missing and exits with 1 then, else 0 (the tests
                 report NOTRUN on a machine without them)

As a module: pv_phase(records, phase, T) -> log10(p/Pa), melting_point(
records), and read_tflows(text) -> (T list, p list [bar], tokens).

File format (vapourPressureTable.H): '#' lines are comments; the '#' line
right before the first data row names the columns, separated by commas;
column 0 is T [K].

Licences: thermo.inp and the data derived from it are Apache-2.0 (NASA
CEA); the T-Flows values are copied from its source for comparison only.
"""

import argparse
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_gems3k_input import R, P0, read_nasa9, nasa_props   # noqa: E402

GAS = 'PbI2'
PHASES = (('PbI2(cr)', 'PbI2(cr)'), ('PbI2(l)', 'PbI2(L)'))   # column, record
T_MIN, T_MAX, STEP = 270.0, 1250.0, 5.0

DEFAULT_NASA = os.path.expanduser(
    '~/opt/gems-data/research-2026-09-23/format/nasa_thermo.inp')
DEFAULT_TFLOWS_REPO = os.path.expanduser('~/src/T-Flows')
# the tip of branch advanced_scalar_model when the table was made (last
# change of Types.f90: 2026-01-25)
TFLOWS_REVISION = '35c1f6406e8eef050b894633bcff9db63e21f06a'
DEFAULT_TFLOWS_SPEC = (TFLOWS_REVISION + ':Tests/Laminar/'
                       'Scalar_Transport_Gpu/LESTO_HKS/User_Mod/Types.f90')


def read_records(path):
    """The NASA-9 records of the gas and of the condensed phases."""
    return read_nasa9(path, [GAS] + [record for _, record in PHASES])


def pv_phase(records, record, T):
    """log10(p/Pa) of PbI2(g) over the pure phase `record` at T [K]."""
    g_gas = nasa_props(records[GAS], T)[3]
    g_phase = nasa_props(records[record], T)[3]
    return -(g_gas - g_phase)/(R*T*math.log(10.0)) + math.log10(P0)


def pv_stable(records, T):
    """log10(p_eq/Pa) over the stable condensate: the minimum over phases."""
    return min(pv_phase(records, record, T) for _, record in PHASES)


def melting_point(records, lo=600.0, hi=750.0):
    """T [K] where the crystal and liquid curves cross (bisection to 1e-12 K)."""
    cr, liq = PHASES[0][1], PHASES[1][1]

    def f(T):
        return pv_phase(records, cr, T) - pv_phase(records, liq, T)

    flo = f(lo)
    if flo*f(hi) > 0:
        raise ValueError(f'the phase curves do not cross in [{lo}, {hi}] K')
    while hi - lo > 1e-12:
        mid = 0.5*(lo + hi)
        if f(mid)*flo > 0:
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5*(lo + hi)


def nodes(t_melt):
    """Every 5 K from 270 to 1250 K, plus the melting node (rounded to 1e-6 K
    unless it falls on a 5 K node)."""
    n = int(round((T_MAX - T_MIN)/STEP))
    grid = [T_MIN + STEP*k for k in range(n + 1)]
    t_node = round(t_melt, 6)
    if all(abs(t - t_node) > 1e-9 for t in grid):
        grid.append(t_node)
    return sorted(grid)


def write_phases(records, path, source):
    t_melt = melting_point(records)
    grid = nodes(t_melt)
    with open(path, 'w') as f:
        f.write('# PbI2 equilibrium vapour pressure over its condensed '
                'phases, log10(p/Pa) of PbI2(g)\n')
        f.write('# Gurvich et al. (1991) through the NASA-9 records '
                'PbI2, PbI2(cr), PbI2(L) of NASA Glenn thermo.inp\n')
        f.write(f'# ({source}); p0 = 1 bar, R = {R} J/(mol K); '
                'Apache-2.0 (NASA CEA)\n')
        f.write('# each column is one phase over the whole range '
                '(crystal extrapolated above T_m, liquid below);\n')
        f.write('# p_eq = the minimum over the phases.  Melting point '
                f'(crossing of the curves) T_m = {t_melt:.6f} K, a node.\n')
        f.write('# written by thermochemistry/systems/make_pv_table.py\n')
        f.write('# T_K,' + ','.join(name for name, _ in PHASES) + '\n')
        for T in grid:
            values = [pv_phase(records, record, T) for _, record in PHASES]
            f.write(f'{T:.6f},' + ','.join(f'{v:.12f}' for v in values)
                    + '\n')
    return t_melt, grid


def read_tflows(text):
    """GEMS_TEMPERATURE [K] and GEMS_VAPOR_PRESSURE [bar] of Types.f90: the
    values as floats and as the source tokens."""
    def array(name):
        start = text.index(name)
        body = text[text.index('(/', start) + 2:text.index('/)', start)]
        tokens = [t.strip() for t in body.replace('&', ' ').replace('\n', ' ')
                  .split(',') if t.strip()]
        return tokens

    t_tokens = array('GEMS_TEMPERATURE')
    p_tokens = array('GEMS_VAPOR_PRESSURE')
    if len(t_tokens) != len(p_tokens):
        raise ValueError(f'{len(t_tokens)} temperatures, {len(p_tokens)} '
                         'pressures')
    return ([float(t) for t in t_tokens], [float(p) for p in p_tokens],
            (t_tokens, p_tokens))


def missing_sources(args):
    """The sources that cannot be read, as text (empty: both present)."""
    missing = []
    if not os.path.isfile(args.nasa):
        missing.append(f'thermo.inp {args.nasa}')
    if args.tflows:
        if not os.path.isfile(args.tflows):
            missing.append(f'Types.f90 {args.tflows}')
    else:
        found = subprocess.run(['git', '-C', args.tflows_git, 'cat-file',
                                '-e', DEFAULT_TFLOWS_SPEC],
                               capture_output=True).returncode == 0
        if not found:
            missing.append(f'{DEFAULT_TFLOWS_SPEC} in {args.tflows_git}')
    return missing


def tflows_source(args):
    if args.tflows:
        with open(args.tflows) as f:
            return f.read(), args.tflows
    text = subprocess.run(['git', '-C', args.tflows_git, 'show',
                           DEFAULT_TFLOWS_SPEC], check=True,
                          capture_output=True, text=True).stdout
    return text, f'git -C {args.tflows_git} show {DEFAULT_TFLOWS_SPEC}'


def write_tflows(text, path, source):
    temps, pressures, (t_tokens, p_tokens) = read_tflows(text)
    for k, T in enumerate(temps):
        if abs(T - (283.15 + 10*k)) > 1e-9:
            raise ValueError(f'node {k}: {T} K is not 283.15 + 10 k')
    with open(path, 'w') as f:
        f.write('# PbI2 vapour pressure of the T-Flows reference '
                'implementation, GEMS_VAPOR_PRESSURE [bar]\n')
        f.write(f'# at its {len(temps)} nodes 283.15 + 10 k K, copied '
                'digit for digit from\n')
        f.write(f'# {source}\n')
        f.write('# (T-Flows interpolates p linearly in T and holds the end '
                'values: interpolation linearT; outOfRange hold; units '
                'bar)\n')
        f.write('# written by thermochemistry/systems/make_pv_table.py\n')
        f.write('# T_K,GEMS_VAPOR_PRESSURE\n')
        for t, p in zip(t_tokens, p_tokens):
            f.write(f'{t},{p}\n')
    return temps, pressures


def interpolate_table(grid, columns, T):
    """log10 p_eq of the table as vapourPressureTable.H evaluates it: per
    phase log10 p linear in 1/T (end intervals extrapolated), then the
    minimum over the phases."""
    n = len(grid)
    i = 0
    if T >= grid[-1]:
        i = n - 2
    elif T > grid[0]:
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi)//2
            if grid[mid] > T:
                hi = mid
            else:
                lo = mid
        i = lo
    w = (1.0/T - 1.0/grid[i])/(1.0/grid[i + 1] - 1.0/grid[i])
    return min((1 - w)*c[i] + w*c[i + 1] for c in columns)


def check(records, grid, t_melt):
    columns = [[pv_phase(records, record, T) for T in grid]
               for _, record in PHASES]
    # the values as written (12 decimals)
    columns = [[float(f'{v:.12f}') for v in c] for c in columns]
    worst, where = 0.0, None
    k = 0
    while True:
        T = T_MIN + 0.1*k
        if T > T_MAX + 1e-9:
            break
        err = abs(interpolate_table(grid, columns, T) - pv_stable(records, T))
        if err > worst:
            worst, where = err, T
        k += 1
    print(f'melting point (crossing of the NASA-9 curves): {t_melt:.6f} K')
    print(f'5 K table (+ T_m node), log10 p linear in 1/T per phase, min: '
          f'max |error| vs NASA-9 = {worst:.3e} decades at {where:.1f} K '
          f'(0.1 K grid, {T_MIN:g}-{T_MAX:g} K)')
    return worst


def write_reference(records, path, step):
    with open(path, 'w') as f:
        f.write('# NASA-9 (Gurvich 1991) log10(p_eq/Pa) of PbI2 over the '
                'stable condensate (minimum over PbI2(cr), PbI2(l))\n')
        f.write('# written by thermochemistry/systems/make_pv_table.py '
                '--reference\n')
        f.write('# T_K,log10p\n')
        k = 0
        while True:
            T = T_MIN + step*k
            if T > T_MAX + 1e-9:
                break
            f.write(f'{T!r},{pv_stable(records, T)!r}\n')
            k += 1


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--nasa', default=DEFAULT_NASA)
    ap.add_argument('--tflows')
    ap.add_argument('--tflows-git', default=DEFAULT_TFLOWS_REPO)
    ap.add_argument('--out', default=os.path.join(here, 'data'))
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--reference')
    ap.add_argument('--reference-step', type=float, default=0.25)
    ap.add_argument('--no-write', action='store_true',
                    help='do not (re)write the csv files')
    ap.add_argument('--sources', action='store_true')
    args = ap.parse_args()

    if args.sources:
        missing = missing_sources(args)
        print('sources: ' + ('; '.join(missing) + ' not found' if missing
                             else f'{args.nasa} and {DEFAULT_TFLOWS_SPEC} '
                             f'of {args.tflows_git}'))
        sys.exit(1 if missing else 0)

    records = read_records(args.nasa)
    t_melt = melting_point(records)
    grid = nodes(t_melt)

    if not args.no_write:
        phases = os.path.join(args.out, 'pv_PbI2_phases.csv')
        write_phases(records, phases, f'file {os.path.basename(args.nasa)}')
        print(f'Wrote {phases}: {len(grid)} nodes, T_m = {t_melt:.6f} K')
        text, source = tflows_source(args)
        tflows = os.path.join(args.out, 'pv_PbI2_TFlows.csv')
        temps, _ = write_tflows(text, tflows, source)
        print(f'Wrote {tflows}: {len(temps)} nodes {temps[0]}-{temps[-1]} K')

    if args.reference:
        write_reference(records, args.reference, args.reference_step)
        print(f'Wrote {args.reference}')

    if args.check:
        worst = check(records, grid, t_melt)
        ok = worst < 2e-3 and abs(t_melt - 683.0) <= 0.1
        print('check ' + ('passed' if ok else 'FAILED')
              + ' (limits: 2e-3 decades; T_m 683.0 +- 0.1 K)')
        sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
