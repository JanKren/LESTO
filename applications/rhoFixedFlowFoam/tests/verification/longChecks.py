#!/usr/bin/env python3
"""
tests/verification/longChecks.py

PURPOSE

Numerical checks of the long runs of milestone M7
(tests/verification/longStudies, 'tests/Alltest M7long'): the runs of both
families (the time-step study of acceptance 3 on the 007 carrier, run/012
and run/013; on the paper-mode carrier the h vs h/2 first-layer study of
acceptance 2, the axial x2 study and the time-step study of acceptance 3
on run/014).  The profiles
(postProcessing/axialProfiles) are read with utilities/axialProfiles.py,
the balance files and the fields with tests/phaseChange/checks.py, the
bin drop of the onset with tests/paperMode/checks.py, and the time-step
measures are those of tests/verification/dtChecks.py.

Measures of the wall-deposit line density n'_w(b) [mol/m] (the column
'wall' of a bin set) at one profile time:

  T_dep, x_dep  the onset fromInlet at 1 % of the peak of the wall
                column, interpolated between the bin centres (plan E10,
                recomputed with axialProfiles.onset); the solver writes the
                same for its observable (run/014: 'wall')
  dT_bin        the drop of the wall temperature across the native bin of
                the onset (tests/paperMode/checks.py binDrop)
  peak mDep_    the largest deposit per wall area [kg/m2] on the WALL faces
                (the field mDep_PbI2_s of the written time, all processors)
  peak          the largest n'_w of a bin set, its bin centre and the
                parabola through it and its neighbours
  L1(a, b)      sum_b |n'_w,a - n'_w,b| dx_b relative to the deposit of b

Subcommands; each prints its lines, the last one 'OK: ...' or 'FAILED:
...', and exits with 0 (met), 1 (not met) or 2 (cannot be computed: a
file is missing).  A crash (an exception) prints 'CRASHED: <exception>'
last and exits with 3, never with 1 as an uncaught exception of python
would, which a reported criterion took for a computed 'not met' (plan
section 35, item N; longStudies requires the status 0 or 1 and a verdict
line).

  runs <closure tol> <defect tol> <case>[=listed] ...
      per run, over all its segments (log.<n>, one balance file
      postProcessing/phaseChangeBalance/<t0>/PbI2_g.dat per start time):
      ended ('End' in the last log), the steps, max |closure|/reference
      in every step of every balance file, |solverDefect|/reference at the
      end and its largest value, the inventories at the end relative to
      the released amount (gas, wall, transport per patch), the time from
      which the gas stays below 1 % of it, and the solver's end-of-run
      summary.  Met when every run ended, the closure of every run not
      marked '=listed' (a revalidated run: listed, not gated, as M4.17)
      is at most <closure tol> and |solverDefect|/reference at the end is
      at most <defect tol> for every such run
  vd <value> <case of A | -> <case of depositionVelocity> ...
      the start-up line 'Temperature model: v_d = ...' of the run with A
      prints <value> (with '(from A = ... on this layer)'), and every
      flux-form run prints the same v_d (given); also the layer lines.
      '-' instead of the case of A: the flux-form runs only
  pair <case h> <case h/2> <time> <max bins> <peak tol> <label> [<la> <lb>]
      acceptance 2 between the base mesh (h) and the h/2 first-layer mesh
      at <time>: the onset shift |x_dep(h/2) - x_dep(h)| in native bins
      (met: at most <max bins>) and in K against dT_bin, the peak mDep_
      ratio (met: within <peak tol>), and, reported, the peaks, L1 of the
      native and 1 cm wall profiles, the deposits and the gas left; <la>
      and <lb> name the two runs in the lines (default h and h/2)
  axial <case base> <case x2> <time> <label> [<max bins> <peak tol> <L1>]
      the axial x2 mesh against the base mesh at <time>: T_dep on the
      native bins (the 424 axial cells of the base mesh, 2 cells of the x2
      mesh each) and on the cells of each mesh (set 'cells' of the x2
      run), the peak height and position of the 1 cm and native profiles,
      L1 of the 1 cm and native wall profiles, the peak mDep_.  With the
      three limits, the verdict of the proposed criterion (onset shift at
      most <max bins> native bins, 1 cm peak height within <peak tol>, L1
      of the 1 cm profiles at most <L1>): exit 0/1; without, exit 0
  dtHistory <case dt=4> <case dt=2> <case dt=1> <t_ev>
      the intermediate-time error of a time-step study of the temperature
      model (plan section 28, item 7: reported, not gated; t_ev the end of
      the release, 6 s on the 007 carrier and 15 s on the paper-mode
      carrier): L1(4,2) and L1(2,1) of the wall
      column (native bins) at every profile time, the largest L1(2,1) over
      the times with a deposit above 1 % of the released amount and its
      observed order, the values at t_ev and t_ev + 1, and L1 of the gas
      column at t_ev + 1 (Fig. 9's gas); exit 0 when computed
  bitwise <case A> <case B> <times> <sets>
      the profile files of two runs at the listed times (comma-separated;
      'all': every common time) and bin sets are identical byte for byte;
      prints the first differing time otherwise (exit 1)
"""

import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', '..', 'utilities'))
sys.path.insert(2, os.path.join(HERE, '..', 'phaseChange'))

import axialProfiles as ap   # noqa: E402
import checks as pc          # noqa: E402  (tests/phaseChange/checks.py)
import dtChecks as dc        # noqa: E402

_spec = importlib.util.spec_from_file_location(
    'paperModeChecks', os.path.join(HERE, '..', 'paperMode', 'checks.py'))
pm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pm)

M_PBI2 = 0.46100894


def finish(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def cannot(text):
    print('FAILED: ' + text)
    sys.exit(2)


# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------

def segmentLogs(case):
    logs = [f for f in os.listdir(case) if re.fullmatch(r'log\.\d+', f)]
    return [os.path.join(case, f)
            for f in sorted(logs, key=lambda f: int(f.split('.')[1]))]


def balanceFiles(case):
    base = os.path.join(case, 'postProcessing', 'phaseChangeBalance')
    if not os.path.isdir(base):
        return []
    return [os.path.join(base, t0, 'PbI2_g.dat')
            for t0 in sorted(os.listdir(base), key=float)
            if os.path.isfile(os.path.join(base, t0, 'PbI2_g.dat'))]


def runFacts(case):
    files = balanceFiles(case)
    if not files:
        raise FileNotFoundError(f'no balance file in {case}')
    rows, names = [], None
    closure = 0.0
    for path in files:
        names, part = pc.readBalance(path)
        t0 = float(part[0][0])
        rows = [r for r in rows if float(r[0]) < t0] + part
        for r in part:
            closure = max(closure, abs(pc.column(names, r, 'closure'))
                          / max(pc.reference(names, r), 1e-300))
    col = {n: i for i, n in enumerate(names)}
    last = rows[-1]
    defects = [abs(r[col['solverDefect']])/max(pc.reference(names, r), 1e-300)
               for r in rows]
    supplied = last[col['initial']] + last[col['released']]
    below = None
    for r in rows:
        if r[col['gas']] >= 0.01*supplied:
            below = None
        elif below is None:
            below = float(r[0])
    logs = segmentLogs(case)
    ended, summary, vdLine = False, '', ''
    for log in logs:
        with open(log, errors='replace') as f:
            for line in f:
                if line.startswith('Phase change summary'):
                    summary = line.strip()
                elif line.startswith('Temperature model: v_d'):
                    vdLine = line.strip()
    if logs:
        with open(logs[-1], errors='replace') as f:
            ended = any(line.strip() == 'End' for line in f)
    facts = {'steps': len(rows), 'end': float(last[0]), 'closureMax': closure,
             'defectEnd': defects[-1], 'defectMax': max(defects),
             'supplied': supplied, 'gasBelow1pct': below, 'ended': ended,
             'summary': summary, 'segments': len(logs), 'vdLine': vdLine}
    for name in ('gas', 'wall', 'sample', 'removed'):
        facts[name] = last[col[name]]/supplied if supplied else math.nan
    facts['transport'] = {n[10:]: last[col[n]]/supplied for n in names
                          if n.startswith('transport_')}
    return facts


M8_DEFECT = 1e-8


def cmdRuns(closureTol, defectTol, *cases):
    closureTol, defectTol = float(closureTol), float(defectTol)
    ok, worstC, worstD, listed, missing = True, (0.0, '-'), (0.0, '-'), [], []
    for item in cases:
        case, _, mark = item.partition('=')
        name = os.path.basename(case.rstrip('/'))
        try:
            f = runFacts(case)
        except (OSError, SystemExit, IndexError, ValueError) as e:
            print(f'{name}: {e}')
            missing.append(name)
            continue
        gated = mark != 'listed'
        good = f['ended'] and (not gated or (f['closureMax'] <= closureTol
                                             and f['defectEnd'] <= defectTol))
        ok = ok and good
        if gated:
            worstC = max(worstC, (f['closureMax'], name))
            worstD = max(worstD, (f['defectEnd'], name))
        else:
            listed.append(f'{name} closure {f["closureMax"]:.1e}, solverDefect '
                          f'{f["defectEnd"]:.1e}')
        transport = ' '.join(f'{k} {v:.3g}' for k, v in f['transport'].items())
        print(f'{name}: {"ended" if f["ended"] else "NOT ENDED"} at t = '
              f'{f["end"]:g} s, {f["steps"]} steps in {f["segments"]} '
              f'segment(s){"" if gated else " (revalidated: listed)"}; max '
              f'|closure|/ref {f["closureMax"]:.2e}, |solverDefect|/ref at '
              f'the end {f["defectEnd"]:.2e} (max {f["defectMax"]:.2e}); of '
              f'the released {f["supplied"]:.6g} kg: gas {f["gas"]:.3e}, wall '
              f'{f["wall"]:.10f}, transport {transport}; gas below 1 % from '
              f't = {f["gasBelow1pct"]} s')
        if f['summary']:
            print(f'  {f["summary"]}')
    if missing:
        cannot(f'no balance file or log for: {", ".join(missing)}')
    if worstC[1] != '-':
        text = (f'{len(cases)} runs ended; largest |closure|/reference '
                f'{worstC[0]:.1e} ({worstC[1]}) in every step (limit '
                f'{closureTol:g}), |solverDefect|/reference at the end '
                f'{worstD[0]:.1e} ({worstD[1]}; limit {defectTol:g})')
    else:
        text = (f'{len(cases)} runs ended; none was made with the current '
                f'sources, so the closure and the solverDefect are listed, '
                f'not gated')
    if listed:
        text += f'; revalidated, not gated: {"; ".join(listed)}'
    finish(ok, text if ok else 'not all runs ended within the limits: '
           + text)


def cmdVd(value, caseA, *fluxCases):
    target = float(value)
    parts, ok = [], True
    runs = [(c, False) for c in fluxCases]
    if caseA != '-':
        runs.insert(0, (caseA, True))
    for case, fromA in runs:
        logs = segmentLogs(case)
        line = ''
        layer = ''
        if logs:
            with open(logs[0], errors='replace') as f:
                for text in f:
                    if text.startswith('Temperature model: v_d') and not line:
                        line = text.strip()
                    elif text.startswith('Interface (layer') and not layer:
                        layer = text.strip()
        m = re.match(r'Temperature model: v_d = (\S+) m/s( \(from A = (\S+) '
                     r'1/s on this layer\))?', line)
        if not m:
            cannot(f'no line "Temperature model: v_d" in the first log of '
                   f'{case}')
        vd, a = float(m.group(1)), m.group(3)
        good = vd == target and (a is not None) == fromA
        ok = ok and good
        lay = re.search(r'A_I/V_I = (\S+) 1/m', layer)
        parts.append(f'{os.path.basename(case.rstrip("/"))} v_d = {m.group(1)}'
                     + (f' from A = {a}' if a else ' given')
                     + (f', A_I/V_I = {lay.group(1)} 1/m' if lay else ''))
    finish(ok, f'the flux form uses depositionVelocity {value} m/s, the v_d '
               f'that A gives on the base mesh as the solver prints it: '
               + '; '.join(parts))


# ---------------------------------------------------------------------------
# Profiles and fields
# ---------------------------------------------------------------------------

def profile(case, time, binSet):
    path = os.path.join(case, 'postProcessing', 'axialProfiles', str(time),
                        f'PbI2_g_{binSet}.dat')
    if not os.path.isfile(path):
        cannot(f'no profile file {path}')
    return ap.read(path)


def wallOnset(p):
    return ap.onset(p, False, values=p.column('wall'), fraction=p.fraction)


def binWidthAt(p, x):
    lo, hi = p.column('x_lo'), p.column('x_hi')
    for a, b in zip(lo, hi):
        if a <= x < b:
            return b - a
    return hi[-1] - lo[-1]


def peakMdep(case, time, patch='WALL', field='mDep_PbI2_s'):
    """the largest value of the field on the patch at <time> [kg/m2], over
    the processor directories (or the case for a serial run)"""
    dirs = sorted(d for d in os.listdir(case)
                  if re.fullmatch(r'processor\d+', d))
    paths = [os.path.join(case, d, str(time), field) for d in dirs] \
        or [os.path.join(case, str(time), field)]
    best = -math.inf
    for path in paths:
        if not os.path.isfile(path):
            cannot(f'no field {path}')
        _, patches = pc.readScalarField(path)
        values = patches.get(patch, [])
        if values:
            best = max(best, max(values))
    return best, len(paths)


def onsetText(o):
    if o.get('status') != 'found':
        return f'none ({o.get("status")})'
    return f'{o["T"]:.4f} K at x = {o["x"]*1e3:.3f} mm'


def compareOnsets(pa, pb):
    """onsets of the wall columns of two profiles (the same bins): the
    shift in m, in bins (the bin of a) and in K, and dT_bin of a"""
    oa, ob = wallOnset(pa), wallOnset(pb)
    if oa.get('status') != 'found' or ob.get('status') != 'found':
        return oa, ob, None
    dx = ob['x'] - oa['x']
    dT = ob['T'] - oa['T']
    width = binWidthAt(pa, oa['x'])
    drop, _ = pm.binDrop(pa, oa['x'])
    return oa, ob, {'dx': dx, 'bins': dx/width, 'dT': dT, 'dTbin': drop,
                    'width': width}


def peakText(p):
    h, x, _ = dc.peak(p)
    return (f'{h:.6e} mol/m at {x*1e3:.3f} mm (parabola '
            f'{dc.peakParabolic(p)*1e3:.3f} mm)')


def cmdPair(caseH, caseH2, time, maxBins, peakTol, label, la='h', lb='h/2'):
    maxBins, peakTol = float(maxBins), float(peakTol)
    nh, nh2 = profile(caseH, time, 'native'), profile(caseH2, time, 'native')
    gh, gh2 = profile(caseH, time, 'gamma'), profile(caseH2, time, 'gamma')
    oh, oh2, s = compareOnsets(nh, nh2)
    if s is None:
        cannot(f'{label}: no onset (h: {onsetText(oh)}; h/2: {onsetText(oh2)})')
    mh, nfh = peakMdep(caseH, time)
    mh2, nfh2 = peakMdep(caseH2, time)
    ratio = mh2/mh - 1 if mh > 0 else math.inf
    l1n, l1g = dc.l1(nh2, nh), dc.l1(gh2, gh)
    print(f'{label}, t = {time} s:')
    print(f'  T_dep fromInlet (wall column, native bins): {la} '
          f'{onsetText(oh)}, '
          f'{lb} {onsetText(oh2)}; shift {s["dx"]*1e3:+.4f} mm = '
          f'{s["bins"]:+.4f} native bins ({s["width"]*1e3:.4f} mm), '
          f'{s["dT"]:+.4f} K (dT_bin {s["dTbin"]:.3f} K)')
    print(f'  T_dep fromInlet (solver, observable {nh.observable}): {la} '
          f'{nh.tdep.get("fromInlet", {}).get("T")} K, {lb} '
          f'{nh2.tdep.get("fromInlet", {}).get("T")} K')
    print(f'  peak mDep_PbI2_s (largest over the WALL faces): {la} {mh:.6e} '
          f'kg/m2, {lb} {mh2:.6e} kg/m2 (from {nfh} and {nfh2} field files): '
          f'{ratio:+.4%}')
    print(f'  peak, native bins: {la} {peakText(nh)}; {lb} {peakText(nh2)}')
    print(f'  peak, 1 cm bins: {la} {peakText(gh)}; {lb} {peakText(gh2)}')
    print(f'  L1 of the wall profiles relative to the deposit of {la}: native '
          f'{l1n:.4e}, 1 cm {l1g:.4e}')
    print(f'  deposit [mol]: {la} {dc.wallInventory(nh):.8e}, {lb} '
          f'{dc.wallInventory(nh2):.8e}; gas left {la} '
          f'{nh.inventory["gas"]:.3e}, '
          f'{lb} {nh2.inventory["gas"]:.3e}')
    met = abs(s['bins']) <= maxBins and abs(ratio) <= peakTol
    finish(met, f'{label}: T_dep shift {s["bins"]:+.3f} native bins '
           f'({s["dx"]*1e3:+.3f} mm, {s["dT"]:+.3f} K; dT_bin '
           f'{s["dTbin"]:.2f} K; limit {maxBins:g} bin), peak mDep_ '
           f'{ratio:+.2%} (limit {peakTol*100:g} %); T_dep {la} '
           f'{oh["T"]:.3f} K, {lb} '
           f'{oh2["T"]:.3f} K; peak mDep_ {la} {mh:.4e}, {lb} {mh2:.4e} kg/m2; '
           f'L1 native {l1n:.2e}, 1 cm {l1g:.2e}')


def cmdAxial(caseB, caseX, time, label, maxBins=None, peakTol=None, l1Tol=None):
    nb, nx = profile(caseB, time, 'native'), profile(caseX, time, 'native')
    gb, gx = profile(caseB, time, 'gamma'), profile(caseX, time, 'gamma')
    cx = profile(caseX, time, 'cells')
    ob, ox, s = compareOnsets(nb, nx)
    ocx = wallOnset(cx)
    if s is None or ocx.get('status') != 'found':
        cannot(f'{label}: no onset')
    hb, xb, _ = dc.peak(gb)
    hx, xx, _ = dc.peak(gx)
    peakG = hx/hb - 1
    shiftG = dc.peakParabolic(gx) - dc.peakParabolic(gb)
    hnb, _, _ = dc.peak(nb)
    hnx, _, _ = dc.peak(nx)
    shiftN = dc.peakParabolic(nx) - dc.peakParabolic(nb)
    l1g, l1n = dc.l1(gx, gb), dc.l1(nx, nb)
    mb, _ = peakMdep(caseB, time)
    mx, _ = peakMdep(caseX, time)
    print(f'{label}, t = {time} s:')
    print(f'  T_dep fromInlet, native bins (the 424 axial cells of the base): '
          f'base {onsetText(ob)}, x2 {onsetText(ox)}; shift '
          f'{s["dx"]*1e3:+.4f} mm = {s["bins"]:+.4f} bins, {s["dT"]:+.4f} K '
          f'(dT_bin {s["dTbin"]:.3f} K)')
    print(f'  T_dep fromInlet on the cells of each mesh: base (424) '
          f'{onsetText(ob)}, x2 (848) {onsetText(ocx)}: '
          f'{ocx["T"] - ob["T"]:+.4f} K, {(ocx["x"] - ob["x"])*1e3:+.4f} mm')
    print(f'  peak, 1 cm bins: base {peakText(gb)}; x2 {peakText(gx)}: '
          f'height {peakG:+.4%}, parabolic position {shiftG*1e3:+.4f} mm')
    print(f'  peak, native bins: base {peakText(nb)}; x2 {peakText(nx)}: '
          f'height {hnx/hnb - 1:+.4%}, parabolic position '
          f'{shiftN*1e3:+.4f} mm')
    print(f'  L1 of the wall profiles relative to the base deposit: 1 cm '
          f'{l1g:.4e}, native {l1n:.4e}')
    print(f'  peak mDep_PbI2_s: base {mb:.6e}, x2 {mx:.6e} kg/m2: '
          f'{mx/mb - 1:+.4%}')
    print(f'  deposit [mol]: base {dc.wallInventory(nb):.8e}, x2 '
          f'{dc.wallInventory(nx):.8e}; gas left base '
          f'{nb.inventory["gas"]:.3e}, x2 {nx.inventory["gas"]:.3e}')
    text = (f'{label}: T_dep shift {s["bins"]:+.3f} native bins '
            f'({s["dT"]:+.3f} K; on the own cells {ocx["T"] - ob["T"]:+.3f} '
            f'K), 1 cm peak height {peakG:+.2%} and position '
            f'{shiftG*1e3:+.3f} mm, L1 of the 1 cm profiles {l1g:.2e} '
            f'(native {l1n:.2e}), peak mDep_ {mx/mb - 1:+.2%}')
    if maxBins is None:
        finish(True, text)
    met = (abs(s['bins']) <= float(maxBins) and abs(peakG) <= float(peakTol)
           and l1g <= float(l1Tol))
    finish(met, text + f'; proposed criterion (shift <= {float(maxBins):g} '
           f'bin, 1 cm peak within {float(peakTol)*100:g} %, L1 <= '
           f'{float(l1Tol):g}): {"met" if met else "not met"}')


def cmdDtHistory(c4, c2, c1, tev):
    tev = float(tev)
    times = dc.commonTimes((c4, c2, c1))
    if not times:
        cannot('no profile time common to the three runs')
    rows = []
    print('# time [s], wall inventory (dt 1 ms) [mol], L1(4,2), L1(2,1), '
          'observed order (native bins)')
    released = None
    for t in times:
        ps = [profile(c, t, 'native') for c in (c4, c2, c1)]
        inv = dc.wallInventory(ps[2])
        pair = sum(ps[2].inventory[k] for k in ('gas', 'wall', 'sample'))
        released = pair
        if not inv > 0:
            print(f'{t} {inv:.6e} - - -')
            continue
        d42, d21 = dc.l1(ps[0], ps[1]), dc.l1(ps[1], ps[2])
        p, _ = dc.order(d42, d21)
        rows.append((float(t), t, inv, pair, d42, d21, p))
        print(f'{t} {inv:.6e} {d42:.4e} {d21:.4e} {p:.3f}')
    n0 = max(r[3] for r in rows) if rows else released
    big = [r for r in rows if r[2] > 0.01*n0]
    if not big:
        cannot('no profile time with a deposit above 1 % of the released '
               'amount')
    worst = max(big, key=lambda r: r[5])
    at = {}
    for target in (tev, tev + 1):
        r = [q for q in rows if abs(q[0] - target) < 1e-9]
        at[target] = (f'{r[0][5]:.3e} (order {r[0][6]:.2f})' if r else 'n/a')
    gasText = 'n/a'
    tg = [t for t in times if abs(float(t) - (tev + 1)) < 1e-9]
    if tg:
        gs = [profile(c, tg[0], 'native') for c in (c4, c2, c1)]
        diff = [math.fsum(abs(a - b)*w for a, b, w in
                          zip(x.column('gas'), y.column('gas'), dc.widths(y)))
                for x, y in ((gs[0], gs[1]), (gs[1], gs[2]))]
        norm = math.fsum(abs(v)*w for v, w in
                         zip(gs[2].column('gas'), dc.widths(gs[2])))
        g42, g21 = diff[0]/norm, diff[1]/norm
        gp, _ = dc.order(g42, g21)
        gasText = (f'{g21:.3e} (L1(4,2) {g42:.3e}, order {gp:.2f}; relative to '
                   f'the gas of the 1 ms run)')
    finish(True, f'intermediate-time error of the wall deposit (native bins, '
           f'reported): largest L1(2,1) {worst[5]:.3e} at {worst[1]} s (order '
           f'{worst[6]:.2f}; deposit {worst[2]/n0:.1%} of the released '
           f'amount), at t_ev = {tev:g} s {at[tev]}, at t_ev + 1 '
           f'{at[tev + 1]}; '
           f'gas at t_ev + 1: L1(2,1) {gasText}')


def cmdBitwise(caseA, caseB, times, sets):
    sets = sets.split(',')
    if times == 'all':
        times = dc.commonTimes((caseA, caseB))
    else:
        times = times.split(',')
    n = 0
    for t in times:
        for s in sets:
            a = os.path.join(caseA, 'postProcessing', 'axialProfiles', t,
                             f'PbI2_g_{s}.dat')
            b = os.path.join(caseB, 'postProcessing', 'axialProfiles', t,
                             f'PbI2_g_{s}.dat')
            if not (os.path.isfile(a) and os.path.isfile(b)):
                cannot(f'missing {a} or {b}')
            if open(a, 'rb').read() != open(b, 'rb').read():
                finish(False, f'{a} and {b} differ ({n} identical before)')
            n += 1
    finish(True, f'{n} profile files identical byte for byte ({len(times)} '
           f'times, sets {",".join(sets)})')


COMMANDS = {
    'runs': cmdRuns,
    'vd': cmdVd,
    'pair': cmdPair,
    'axial': cmdAxial,
    'dtHistory': cmdDtHistory,
    'bitwise': cmdBitwise,
}

CRASHED = 3


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    try:
        COMMANDS[sys.argv[1]](*sys.argv[2:])
    except SystemExit:
        raise
    except BaseException as e:      # noqa: B902 (every crash, reported)
        import traceback
        traceback.print_exc()
        sys.stdout.flush()
        print(f'CRASHED: {type(e).__name__}: {e}')
        sys.exit(CRASHED)
