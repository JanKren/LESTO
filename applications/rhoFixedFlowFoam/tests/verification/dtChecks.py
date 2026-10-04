#!/usr/bin/env python3
"""
tests/verification/dtChecks.py

PURPOSE

The measures of the time-step studies of milestone M7 (acceptance 3 of
doc/phase-change-plan.md: "Deposit profile L1 change < 1 % between dt = 2
and 1 ms"), run by tests/verification/longStudies ('tests/Alltest M7long')
on three runs at dt = 4, 2 and 1 ms of one case: copies of run/012
(temperature model, release sample) and run/013 (HKS, inventory sample
with removal) on the 007 carrier (criteria M7.3b-d), and run/014 on the
paper-mode carrier (M7.3pb; longChecks.py imports this module).  Every
subcommand prints its lines and exits with 0 (criterion met) or 1.  The
profile files (postProcessing/axialProfiles) are read with
utilities/axialProfiles.py; the runs themselves (closure, solverDefect)
are checked by longChecks.py runs.

Measures, for the wall-deposit line density n'_w(b) [mol/m] (the column
'wall') of a bin set at one profile time:

  inventory  I = sum_b n'_w(b) dx_b [mol] (the wall inventory)
  L1(a, b)   sum_b |n'_w,a(b) - n'_w,b(b)| dx_b / I_b, the L1 difference of
             the runs a and b normalised by the deposit inventory of the
             finer run b (dimensionless; 0.01 = 1 %)
  peak       the first maximum of n'_w: its height and its bin centre,
             and its position from the parabola through the peak bin and
             its neighbours (a continuous measure of a shift)
  centroid   the first moment sum x n'_w dx/I of the deposit
  T_dep      the onset at 'fraction' of the peak, both search modes
             (fromInlet, fromPeak; plan E10), of the wall column, recomputed
             with axialProfiles.onset (the solver writes T_dep of its
             observable, wallPlusGas by default, which is reported too)
  order      p = log2(D(4,2)/D(2,1)) from the three time steps, D the L1
             difference or |f(dt) - f(dt/2)| of a scalar f; the
             Richardson estimate of the error left at 1 ms is
             D(2,1)/(2^p - 1)

  study <name> <time> <limit> <case dt=4> <case dt=2> <case dt=1>
        [<set> ...]
      the table of the measures above for the bin sets (default native
      gamma) at the profile directory <time> ('latest': the last profile
      time common to the three runs); 0 when L1(2 ms, 1 ms) < limit for
      every set
  history <set> <case dt=4> <case dt=2> <case dt=1>
      L1(4,2), L1(2,1) and the wall inventory at every profile time common
      to the three runs (always 0)
  plot <file> <time> <set> <case dt=4> <case dt=2> <case dt=1> [<title>]
      the three wall line densities and their differences at <time>
      ('latest' as in study; needs matplotlib; always 0)

  study prints, besides the tables, one line 'computed <name> <set>: yes'
  per bin set whose three L1 differences are finite, and one line
  'summary <name> at <time> s: ...' with the measures between 2 and 1 ms
  (L1, observed order, Richardson estimate, peak height and bin shift,
  T_dep shift of both modes) for the criterion lines of longStudies.
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'utilities'))

import axialProfiles as ap  # noqa: E402

M_PBI2 = 0.46100894


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------

def profile(case, time, binSet, pair='PbI2_g'):
    path = os.path.join(case, 'postProcessing', 'axialProfiles', str(time),
                        f'{pair}_{binSet}.dat')
    if not os.path.isfile(path):
        raise SystemExit(f'No profile file {path}')
    return ap.read(path)


def widths(p):
    return [b - a for a, b in zip(p.column('x_lo'), p.column('x_hi'))]


def wallInventory(p):
    return math.fsum(v*w for v, w in zip(p.column('wall'), widths(p)))


def l1(pa, pb):
    """L1 difference of the wall columns, relative to the inventory of b"""
    if pa.column('x_lo') != pb.column('x_lo') \
            or pa.column('x_hi') != pb.column('x_hi'):
        raise SystemExit(f'Different bins: {pa.path} {pb.path}')
    inventory = wallInventory(pb)
    diff = math.fsum(abs(a - b)*w for a, b, w in
                     zip(pa.column('wall'), pb.column('wall'), widths(pb)))
    return diff/inventory if inventory > 0 else math.inf


def peak(p):
    w = p.column('wall')
    b = max(range(len(w)), key=lambda i: (w[i], -i))
    return w[b], p.column('x')[b], b


def peakParabolic(p):
    """the peak position from the parabola through the peak bin and its
    neighbours (bins of equal width), a continuous measure of its shift"""
    w, x = p.column('wall'), p.column('x')
    _, xb, b = peak(p)
    if b == 0 or b == len(w) - 1:
        return xb
    curvature = w[b - 1] - 2*w[b] + w[b + 1]
    if not curvature < 0:
        return xb
    return xb + 0.5*(x[b + 1] - x[b])*(w[b - 1] - w[b + 1])/curvature


def centroid(p):
    """the first moment of the deposit, sum x n'_w dx/sum n'_w dx [m]"""
    wd = [v*d for v, d in zip(p.column('wall'), widths(p))]
    total = math.fsum(wd)
    return (math.fsum(x*v for x, v in zip(p.column('x'), wd))/total
            if total > 0 else math.nan)


def wallOnsets(p):
    """T_dep and x_dep of the wall column, both modes (recomputed)"""
    out = {}
    for mode in ('fromInlet', 'fromPeak'):
        searchFrom = p.tdep.get(mode, {}).get('searchFrom')
        out[mode] = ap.onset(p, mode == 'fromPeak', searchFrom,
                             values=p.column('wall'))
    return out


def order(d42, d21):
    """observed order and the Richardson estimate of the error at dt/2"""
    if not (d42 > 0 and d21 > 0):
        return math.nan, math.nan
    p = math.log2(d42/d21)
    if p <= 0:
        return p, math.nan
    return p, d21/(2**p - 1)


def fmtT(o):
    if o.get('status') != 'found':
        return f'none ({o.get("status")})'
    T = o.get('T')
    return (f'{T:.3f} K at {o["x"]*1e3:.3f} mm' if T is not None
            else f'(no wall face) at {o["x"]*1e3:.3f} mm')


def studyTable(name, time, cases, binSet):
    """the measures of one bin set; returns (lines, L1(2,1))"""
    ps = [profile(c, time, binSet) for c in cases]
    labels = ['4 ms', '2 ms', '1 ms']
    lines = [f'{name}, bin set {binSet} ({ps[0].nBins} bins), t = {time} s:']
    inv = [wallInventory(p) for p in ps]
    pair = [max(sum(p.inventory[k] for k in ('gas', 'wall', 'sample')),
                1e-300) for p in ps]
    lines.append('  wall inventory [mol]: ' + ', '.join(
        f'{l} {v:.8g}' for l, v in zip(labels, inv))
        + '; of the pair inventory ' + ', '.join(
        f'{v/s:.6f}' for v, s in zip(inv, pair)))
    d42, d21 = l1(ps[0], ps[1]), l1(ps[1], ps[2])
    p, rich = order(d42, d21)
    lines.append(f'  L1(4,2) {d42:.4e}  L1(2,1) {d21:.4e}  observed order '
                 f'{p:.3f}  Richardson error at 1 ms {rich:.3e}')
    pk = [peak(q) for q in ps]
    lines.append('  peak [mol/m] at x [mm]: ' + ', '.join(
        f'{l} {h:.8g} at {x*1e3:.3f}' for l, (h, x, _) in zip(labels, pk)))
    hp, _ = order(abs(pk[0][0] - pk[1][0]), abs(pk[1][0] - pk[2][0]))
    h42 = (pk[1][0] - pk[0][0])/pk[1][0] if pk[1][0] else math.nan
    h21 = (pk[2][0] - pk[1][0])/pk[2][0] if pk[2][0] else math.nan
    lines.append(f'  peak height change: 4->2 {h42:+.4e}, 2->1 {h21:+.4e}'
                 f' (relative); order {hp:.3f}; peak bin shift 4->2 '
                 f'{pk[1][2] - pk[0][2]:+d}, 2->1 {pk[2][2] - pk[1][2]:+d}')
    for label, f in (('parabolic peak position', peakParabolic),
                     ('centroid', centroid)):
        v = [f(q) for q in ps]
        vp, _ = order(abs(v[0] - v[1]), abs(v[1] - v[2]))
        lines.append(f'  {label} [mm]: ' + ', '.join(
            f'{l} {u*1e3:.4f}' for l, u in zip(labels, v))
            + f'; change 4->2 {(v[1] - v[0])*1e3:+.3e} mm, 2->1 '
            f'{(v[2] - v[1])*1e3:+.3e} mm, order {vp:.3f}')
    onsets = [wallOnsets(q) for q in ps]
    for mode in ('fromInlet', 'fromPeak'):
        Ts = [o[mode].get('T') if o[mode].get('status') == 'found' else None
              for o in onsets]
        text = ', '.join(f'{l} {fmtT(o[mode])}' for l, o in zip(labels, onsets))
        if all(T is not None for T in Ts):
            tp, _ = order(abs(Ts[0] - Ts[1]), abs(Ts[1] - Ts[2]))
            text += (f'; change 4->2 {Ts[1] - Ts[0]:+.3e} K, 2->1 '
                     f'{Ts[2] - Ts[1]:+.3e} K, order {tp:.3f}')
        lines.append(f'  T_dep {mode} (wall column): {text}')
    for mode in ('fromInlet', 'fromPeak'):
        written = [q.tdep.get(mode, {}) for q in ps]
        text = ', '.join(
            f'{l} ' + (f'{w["T"]:.3f} K' if w.get('status') == 'found'
                       and w.get('T') is not None else
                       f'none ({w.get("status")})')
            for l, w in zip(labels, written))
        lines.append(f'  T_dep {mode} (solver, observable {ps[0].observable}):'
                     f' {text}')
    def shift(mode):
        if all(o[mode].get('status') == 'found' and o[mode].get('T')
               is not None for o in onsets[1:]):
            return f'{onsets[2][mode]["T"] - onsets[1][mode]["T"]:+.2e} K'
        return 'n/a'
    summary = (f'{binSet}: L1(2,1) {d21:.2e}, order {p:.2f}, Richardson '
               f'{rich:.1e}; peak 2->1 {h21:+.1e} relative, '
               f'{pk[2][2] - pk[1][2]:+d} bins; T_dep 2->1 fromInlet '
               f'{shift("fromInlet")}, fromPeak {shift("fromPeak")}')
    return lines, d21, summary, math.isfinite(d42) and math.isfinite(d21)


def cmdStudy(name, time, limit, c4, c2, c1, *sets):
    limit = float(limit)
    sets = sets or ('native', 'gamma')
    if time == 'latest':
        times = commonTimes((c4, c2, c1))
        if not times:
            done(False, f'{name}: no profile time common to the three runs')
        time = times[-1]
    ok = True
    worst = []
    summaries = []
    for s in sets:
        lines, d21, summary, finite = studyTable(name, time, (c4, c2, c1), s)
        print('\n'.join(lines))
        ok = ok and d21 < limit
        worst.append(f'{s} {d21:.3e}')
        summaries.append(summary)
        print(f'computed {name} {s}: {"yes" if finite else "no"}')
    print(f'summary {name} at {time} s: ' + '; '.join(summaries))
    done(ok, f'{name}: L1 change of the wall deposit between dt = 2 and 1 ms '
             f'at t = {time} s: {", ".join(worst)} (limit {limit:g})')


def commonTimes(cases):
    sets = []
    for c in cases:
        base = os.path.join(c, 'postProcessing', 'axialProfiles')
        sets.append({d for d in os.listdir(base)
                     if re.fullmatch(r'[0-9.eE+-]+', d)})
    return sorted(set.intersection(*sets), key=float)


def cmdHistory(binSet, c4, c2, c1):
    print(f'# set {binSet}: time [s], wall inventory of dt=1 ms [mol], '
          'L1(4,2), L1(2,1), observed order')
    for t in commonTimes((c4, c2, c1)):
        ps = [profile(c, t, binSet) for c in (c4, c2, c1)]
        inv = wallInventory(ps[2])
        if not inv > 0:
            print(f'{t} {inv:.6e} - - -')
            continue
        d42, d21 = l1(ps[0], ps[1]), l1(ps[1], ps[2])
        p, _ = order(d42, d21)
        print(f'{t} {inv:.6e} {d42:.4e} {d21:.4e} {p:.3f}')
    sys.exit(0)


def cmdPlot(output, time, binSet, c4, c2, c1, title=''):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('plot: no matplotlib')
        sys.exit(0)
    if time == 'latest':
        times = commonTimes((c4, c2, c1))
        if not times:
            print('plot: no common profile time')
            sys.exit(0)
        time = times[-1]
    ps = [profile(c, time, binSet) for c in (c4, c2, c1)]
    x = [v*1e3 for v in ps[0].column('x')]
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for p, label, style in zip(ps, ('dt = 4 ms', 'dt = 2 ms', 'dt = 1 ms'),
                               ('-', '--', ':')):
        ax.step(x, p.column('wall'), where='mid', linestyle=style,
                label=label)
    ax.set_ylabel("wall deposit n'_w [mol/m]")
    ax.legend()
    ax.set_title(title or f'{binSet} bins, t = {time} s')
    for a, b, label in ((ps[0], ps[1], '4 ms - 2 ms'),
                        (ps[1], ps[2], '2 ms - 1 ms')):
        bx.step(x, [u - v for u, v in zip(a.column('wall'), b.column('wall'))],
                where='mid', label=label)
    bx.set_xlabel('x [mm]')
    bx.set_ylabel('difference [mol/m]')
    bx.legend()
    for axis in (ax, bx):
        axis.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f'Wrote {output}')
    sys.exit(0)


COMMANDS = {
    'study': cmdStudy,
    'history': cmdHistory,
    'plot': cmdPlot,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    COMMANDS[sys.argv[1]](*sys.argv[2:])
