#!/usr/bin/env python3
"""
utilities/axialProfiles.py

PURPOSE

Reader and plotter of the axial profiles of rhoFixedFlowFoam
(postProcessing/axialProfiles/<time>/<gas>_<set>.dat, written by
axialProfiles.H at the write times and at the profile times of
constant/thermochemistryProperties).  For every file it prints the time,
the pair, the bin set, T_dep of every search mode (as written by the solver
and recomputed here from the columns, an independent implementation of
plan E10), and the integrals sum_b n'(b) dx_b M of the gas, the wall and
the sample next to the inventories of the file header.  For a file of the
observable wallPlusGas it checks that the observable is the wall deposit
plus the gas outside the sample cells (wall + gas - sampleGas; plan E10 as
decided with the sign-off of M3, section 19, item 6: the gas of the cells
of every sample is excluded) and also prints the onset of the wall column
alone ('wall column: T_dep ...').

Significance (axialProfiles.H): an onset search whose peak bin holds at
most 'minimumAmount' (significance x the pair inventory, header line
'# significance') reports 'insignificant' instead of an onset; the
recomputation applies the same rule.

Usage

  axialProfiles.py <case> [--time T ...] [--set NAME] [--pair GAS]
                   [--plot FILE] [--column observable|gas|wall|sample]

  <case>     the case directory (or its postProcessing/axialProfiles)
  --time     only these time directories (names as written)
  --set      only this bin set (e.g. native, gamma)
  --pair     only this gas species (e.g. PbI2_g)
  --plot     plot the chosen column [mol/m] of every selected file against
             x into FILE (png, pdf, ...; needs matplotlib), with T_dep
             marked; one figure, one curve per file
  --column   the column to plot (default observable)

Only the Python standard library is needed; matplotlib only for --plot.

As a module

  import axialProfiles as ap
  p = ap.read('case/postProcessing/axialProfiles/7/PbI2_g_gamma.dat')
  p.column('wall'), p.time, p.inventory['gas'], p.tdep['fromPeak']['T']
  ap.integrals(p)        # {'gas': mol, 'wall': mol, 'sample': mol}
  ap.onset(p, fromPeak=True)   # recomputed onset: dict with T, x, peak, ...

File format (axialProfiles.H): header lines '# key values', e.g.

  # time 7.0000000000000009 directory 7 reason profile time 7
  # molarMass 0.46100894
  # inventory gas <mol> wall <mol> sample <mol>
  # binned gas ... / # outside gas ...
  # observable wallPlusGas
  # fraction 0.01
  # significance 1e-12 minimumAmount <mol>
  # Tdep fromInlet status found T <K> x <m> peak <mol/m> xPeak <m>
    threshold <mol/m> [searchFrom <m>]
    (status noPeak | notBracketed | insignificant: no T and x)
  # columns x_lo x_hi x gas wall sample observable Twall wallArea sampleGas

then one row per bin (files of the M3 build have no sampleGas column and
no significance line; the reader accepts them).
"""

import argparse
import math
import os
import sys


class Profile:
    """One profile file: header values and the columns of the bins."""

    def __init__(self, path):
        self.path = path
        self.pair = None
        self.condensate = None
        self.set = None
        self.time = None
        self.directory = None
        self.reason = ''
        self.molarMass = None
        self.nBins = None
        self.inventory = {}
        self.binned = {}
        self.outside = {}
        self.observable = None
        self.fraction = None
        self.significance = None
        self.minimumAmount = None
        self.tdep = {}
        self.columns = []
        self.rows = []

    def column(self, name):
        i = self.columns.index(name)
        return [row[i] for row in self.rows]


def _number(text):
    return None if text == 'none' else float(text)


def _pairs(words):
    """'key value key value ...' -> dict (values as floats where possible)"""
    out = {}
    for key, value in zip(words[0::2], words[1::2]):
        try:
            out[key] = _number(value)
        except ValueError:
            out[key] = value
    return out


def read(path):
    """Reads one profile file."""
    p = Profile(path)
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if not line.startswith('#'):
                p.rows.append([float(x) for x in line.split()])
                continue
            words = line[1:].split()
            if not words:
                continue
            key = words[0]
            if key == 'axial' and 'pair' in words:
                i = words.index('pair')
                p.pair = words[i + 1]
                p.condensate = words[i + 3].rstrip(',')
                p.set = words[words.index('set') + 1]
            elif key == 'time':
                p.time = float(words[1])
                p.directory = words[3]
                p.reason = ' '.join(words[5:])
            elif key == 'molarMass':
                p.molarMass = float(words[1])
            elif key == 'nBins':
                p.nBins = int(words[1])
            elif key in ('inventory', 'binned', 'outside'):
                getattr(p, key).update(_pairs(words[1:]))
            elif key == 'observable':
                p.observable = words[1]
            elif key == 'fraction':
                p.fraction = float(words[1])
            elif key == 'significance':
                p.significance = float(words[1])
                p.minimumAmount = float(words[3])
            elif key == 'Tdep':
                p.tdep[words[1]] = _pairs(words[2:])
            elif key == 'columns':
                p.columns = words[1:]
    if p.nBins is not None and len(p.rows) != p.nBins:
        raise ValueError(f'{path}: {len(p.rows)} rows, nBins {p.nBins}')
    return p


def integrals(p):
    """sum_b n'(b) (x_hi - x_lo) of the gas, wall and sample columns [mol]"""
    lo, hi = p.column('x_lo'), p.column('x_hi')
    out = {}
    for name in ('gas', 'wall', 'sample'):
        values = p.column(name)
        out[name] = math.fsum(v*(b - a) for v, a, b in zip(values, lo, hi))
    return out


def onset(p, fromPeak, searchFrom=None, fraction=None, values=None):
    """
    Onset of the deposit recomputed from the columns (plan E10): the
    crossing of fraction*peak of the observable (or the given values),
    interpolated linearly between the bin centres.  fromInlet searches the
    bins whose centre is at or after searchFrom (default: all); fromPeak
    walks upstream from the first maximum of all bins.  A peak bin holding
    at most the file's minimumAmount [mol] is insignificant.  Returns a
    dict with 'status' (found | noPeak | notBracketed | insignificant) and,
    when found, 'x' and 'T'.
    """
    x = p.column('x')
    widths = [b - a for a, b in zip(p.column('x_lo'), p.column('x_hi'))]
    T = p.column('Twall')
    area = p.column('wallArea')
    O = values if values is not None else p.column('observable')
    f = p.fraction if fraction is None else fraction
    n = len(O)

    first = 0
    if not fromPeak and searchFrom is not None:
        while first < n and x[first] < searchFrom:
            first += 1

    peakBin = None
    for b in range(first, n):
        if peakBin is None or O[b] > O[peakBin]:
            peakBin = b
    if peakBin is None or not O[peakBin] > 0:
        return {'status': 'noPeak'}
    peak = O[peakBin]
    threshold = f*peak
    if p.minimumAmount is not None \
            and not peak*widths[peakBin] > p.minimumAmount:
        return {'status': 'insignificant', 'peak': peak,
                'xPeak': x[peakBin], 'threshold': threshold}

    if fromPeak:
        b = peakBin
        while b > first and O[b - 1] >= threshold:
            b -= 1
    else:
        b = first
        while O[b] < threshold:
            b += 1
    result = {'peak': peak, 'xPeak': x[peakBin], 'threshold': threshold}
    if b == first:
        result['status'] = 'notBracketed'
        return result

    w = (threshold - O[b - 1])/(O[b] - O[b - 1])
    result['status'] = 'found'
    result['x'] = x[b - 1] + w*(x[b] - x[b - 1])
    result['T'] = (T[b - 1] + w*(T[b] - T[b - 1])
                   if area[b - 1] > 0 and area[b] > 0 else None)
    return result


def files(case, times=None, binSet=None, pair=None):
    """The profile files of a case, sorted by time, set and pair."""
    root = case
    candidate = os.path.join(case, 'postProcessing', 'axialProfiles')
    if os.path.isdir(candidate):
        root = candidate
    found = []
    for name in os.listdir(root):
        directory = os.path.join(root, name)
        if not os.path.isdir(directory):
            continue
        try:
            t = float(name)
        except ValueError:
            continue
        if times and name not in times:
            continue
        for f in sorted(os.listdir(directory)):
            if not f.endswith('.dat') or '_' not in f:
                continue
            stem = f[:-4]
            if binSet is not None and not stem.endswith('_' + binSet):
                continue
            if pair is not None and not stem.startswith(pair + '_'):
                continue
            found.append((t, f, os.path.join(directory, f)))
    return [path for _, _, path in sorted(found)]


def describe(p):
    """Printed summary of one file."""
    lines = [f'{p.pair} -> {p.condensate}, set {p.set}, time {p.time!r} '
             f'(directory {p.directory}; {p.reason})']
    total = integrals(p)
    parts = []
    for name in ('gas', 'wall', 'sample'):
        inv = p.inventory.get(name, 0.0)
        out = p.outside.get(name, 0.0)
        diff = total[name] + out - inv
        rel = 0.0 if diff == 0 else diff/abs(inv) if inv else float('inf')
        parts.append(f'{name} {total[name]:.10g} mol (+ outside {out:.3g}; '
                     f'inventory {inv:.10g}; relative {rel:.2e})')
    lines.append('  integrals: ' + '; '.join(parts))
    for mode in ('fromInlet', 'fromPeak'):
        if mode not in p.tdep:
            continue
        written = p.tdep[mode]
        again = onset(p, mode == 'fromPeak', written.get('searchFrom'))
        if written.get('status') == 'found':
            T = written.get('T')
            text = (f'T_dep {mode}: {T:.6f} K at x = {written["x"]:.6g} m'
                    if T is not None else
                    f'T_dep {mode}: (no wall face) at x = {written["x"]:.6g} m')
        else:
            text = f'T_dep {mode}: none ({written.get("status")})'
        if again['status'] == 'found' and again.get('T') is not None \
                and written.get('T') is not None:
            text += (f'; recomputed {again["T"]:.6f} K '
                     f'(difference {again["T"] - written["T"]:.2e} K)')
        else:
            text += f'; recomputed: {again["status"]}'
        lines.append('  ' + text)
    if p.observable == 'wallPlusGas' and 'sampleGas' in p.columns:
        worst = 0.0
        scale = max([abs(v) for v in p.column('observable')] + [1e-300])
        for o, w, g, sg in zip(p.column('observable'), p.column('wall'),
                               p.column('gas'), p.column('sampleGas')):
            worst = max(worst, abs(o - (w + g - sg))/scale)
        lines.append(f'  observable = wall + gas outside the sample cells '
                     f'(wall + gas - sampleGas) to {worst:.1e} of its '
                     f'maximum; gas in the sample cells '
                     f'{math.fsum(v*(b - a) for v, a, b in zip(p.column("sampleGas"), p.column("x_lo"), p.column("x_hi"))):.10g} mol')
    if p.observable == 'wallPlusGas':
        for mode in ('fromInlet', 'fromPeak'):
            if mode not in p.tdep:
                continue
            wall = onset(p, mode == 'fromPeak',
                         p.tdep[mode].get('searchFrom'),
                         values=p.column('wall'))
            if wall['status'] == 'found':
                T = wall.get('T')
                text = (f'{T:.6f} K' if T is not None else '(no wall face)')
                text += (f' at x = {wall["x"]:.6g} m (peak'
                         f' {wall["peak"]:.6g} mol/m at x ='
                         f' {wall["xPeak"]:.6g} m)')
            else:
                text = f'none ({wall["status"]})'
            lines.append(f'  wall column: T_dep {mode}: {text}')
    return '\n'.join(lines)


def plot(profiles, column, output):
    """One figure of the column against x for every profile."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit('--plot needs matplotlib')
    fig, ax = plt.subplots(figsize=(7, 4))
    for p in profiles:
        label = f'{p.pair} {p.set} t = {p.time:.6g} s'
        ax.step(p.column('x'), p.column(column), where='mid', label=label)
        for mode, style in (('fromInlet', ':'), ('fromPeak', '--')):
            onsetValues = p.tdep.get(mode, {})
            if onsetValues.get('status') == 'found':
                ax.axvline(onsetValues['x'], linestyle=style, linewidth=0.8,
                           color='grey')
    ax.set_xlabel('x [m]')
    ax.set_ylabel(f'{column} [mol/m]')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f'Wrote {output}')


def main():
    parser = argparse.ArgumentParser(
        description='Read the axial profiles of rhoFixedFlowFoam and print '
                    'T_dep and the integrals.')
    parser.add_argument('case')
    parser.add_argument('--time', nargs='*')
    parser.add_argument('--set')
    parser.add_argument('--pair')
    parser.add_argument('--plot')
    parser.add_argument('--column', default='observable',
                        choices=['observable', 'gas', 'wall', 'sample'])
    args = parser.parse_args()

    paths = files(args.case, args.time, args.set, args.pair)
    if not paths:
        sys.exit(f'No profile files in {args.case}')
    profiles = [read(path) for path in paths]
    for p in profiles:
        print(describe(p))
    if args.plot:
        plot(profiles, args.column, args.plot)


if __name__ == '__main__':
    main()
