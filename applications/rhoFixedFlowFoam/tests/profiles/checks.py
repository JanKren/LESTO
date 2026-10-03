#!/usr/bin/env python3
"""
tests/profiles/checks.py

PURPOSE

Numerical checks of the M3 tests (tests/profiles/Allrun): the axial
profiles, T_dep, mDep_ and exch_ of rhoFixedFlowFoam.  Every subcommand
prints one summary line and exits with 0 (criterion met) or 1.  The
profile files are read with utilities/axialProfiles.py; the balance files
and the ascii/binary fields with tests/phaseChange/checks.py.

  integrals <case> <tol> [<set> ...]
      every profile file of the case (of the given bin sets, default all):
      sum_b n'(b) dx_b M of the gas, the wall and the sample, plus the
      amounts outside the bins, equals the inventory of the balance file
      at the exact time of the profile (the row with the same time) to tol
      relative (an inventory of 0 must give exactly 0); prints the
      largest relative difference and the number of files
  inventoryHeader <case> <tol>
      the inventories written in the profile headers equal the balance
      file at the same time to tol (they are the ledger's own sums)
  outside <case> <set> <lastWidth>
      the bin set does not cover the mesh: in every file the amount
      outside the bins is positive for the gas, and the last bin has the
      given width
  compareProfiles <case A> <case B> <tol> [<Ttol> [<xtol>]]
      the same profile files in both cases (same times, sets, pairs); the
      columns gas, wall, sample and observable of B differ from A by at
      most tol times the largest |value| of the column in A; of every
      onset search the status and the peak position xPeak are equal, T_dep
      differs by at most Ttol K (default 1e-9), x_dep by at most xtol m
      (default 1e-12), the peak and the threshold by at most tol times the
      largest |observable| in A; prints the maxima
  overlap <case> <native set> <cut set> <tol>
      bins that cut cells: in every file of the cut set, every bin's gas,
      wall, sample and observable equal the exact overlap of the native
      profile of the same time directory (whose bins are the cells along
      the axis), sum_c n'_c |c ^ b|/dx_b, to tol times the largest |value|
      of the native column; Twall equals the wall-area-weighted mean
      sum_c A'_c |c ^ b| T_c/sum_c A'_c |c ^ b| (A' = wallArea [m2/m]) to
      tol relative, and wallArea equals sum_c A'_c |c ^ b|/dx_b to tol
      times the largest native wallArea
  times <case> <dt> <write times> <profile times> [<start>]
      comma-separated lists.  postProcessing/axialProfiles holds exactly
      one directory per write time and per profile time (a profile time is
      written at the step end nearest to it, ties to the earlier one, and
      several on one step give one directory); the exact time in every
      file header lies within dt/2 of its requested times; the time
      directories of the case (and of every processor*) are exactly the
      write times (and 0 and <start>); the reason of a write time's
      directory starts with 'write time'; with <start> (a restart) only the
      entries after it are written; entries after the last step (the end
      of the balance file) are not expected
  tdep <case> <set> <mode> <expected T [K]> <tol [K]>
      the T_dep of the search mode in every file of the set equals the
      expected value to tol, and the recomputation from the columns
      (utilities/axialProfiles.py) agrees with the file to tol
  rampField <case> <spike 0|1>
      writes 0/Y_PbI2_g of the synthetic ramp from 0/rho: rho Y = r(x_c)
      with r = 1.2e-3 (x - 0.0305) kg/m3 on cells 30 to 80 of the 100 cells
      of 1 mm, zero elsewhere, and with spike 1 r = 0.3 r(x_80) in cell 10
  rampExpected <mode> <observable> <set> <spike 0|1> [<searchFrom>]
      prints the analytic T_dep [K] of the ramp: T = 1000 - 6000 x at the
      crossing of 0.01 of the peak of the analytic line densities, for the
      native bins (1 mm, the cells) and the coarse bins (5 mm: the means of
      the analytic cell values), interpolated between bin centres; the
      spike (T >= 900 K, f = 0) is gas only, so it is not seen by the
      observable wall
  mDepFactor <case> <factor>
      at every written time: mDep_PbI2_s on WALL equals factor times
      mw_PbI2_s exactly (ascii, 17 digits; factor 2 is a power of two)
  areaIntegral <case> <function object dat> <column> <tol> [<control min>]
      the column (0 = the first value column) of a surfaceFieldValue file
      (areaIntegrate on WALL, one row per step) equals the wall inventory
      of the balance file in every step to tol relative; with <control
      min>, the check is inverted: the relative difference must reach at
      least <control min> (e.g. mw_ with areaPerVolume layer)
  exchange <case> <dt> <tol>
      closed box, writeExchangeField yes: at every written time n >= 1,
      exch_PbI2_g = rho (Y^n - Y^n-1)/dt in every cell to tol relative to
      the largest |exch|, and sum_c exch V dt = -(wall^n - wall^n-1) of the
      balance file to tol relative
  uniformPatch <field file> <patch> [noValue]
      the patch holds 'value uniform <x>' with at most 6 significant
      digits (a rounded single-face value in a binary file); with noValue:
      the patch writes no value at all (fixedGradient of v2412: the value
      is recomputed from the cell when read)
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'utilities'))
sys.path.insert(0, os.path.join(HERE, '..', 'phaseChange'))

import axialProfiles as ap  # noqa: E402
import checks as pc         # noqa: E402

M_PBI2 = 0.46100894


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def balanceRows(case):
    """{time string: row} of the first balance file of a case"""
    base = os.path.join(case, 'postProcessing', 'phaseChangeBalance')
    t0 = sorted(os.listdir(base), key=float)[0]
    names, rows = pc.readBalance(os.path.join(base, t0, 'PbI2_g.dat'))
    return names, {row[0]: row for row in rows}


def rowAt(rows, time):
    """the balance row whose exact time equals the profile's time"""
    for key, row in rows.items():
        if float(key) == time:
            return row
    raise SystemExit(f'No balance row at t = {time!r}')


def profileFiles(case, sets=None):
    paths = ap.files(case)
    if sets:
        paths = [p for p in paths
                 if any(p.endswith('_' + s + '.dat') for s in sets)]
    if not paths:
        raise SystemExit(f'No profile files in {case}')
    return [ap.read(p) for p in paths]


def relative(value, reference):
    difference = value - reference
    if difference == 0:
        return 0.0
    return abs(difference)/abs(reference) if reference else math.inf


def cmdIntegrals(case, tol, *sets):
    tol = float(tol)
    names, rows = balanceRows(case)
    worst = {'gas': 0.0, 'wall': 0.0, 'sample': 0.0}
    profiles = profileFiles(case, sets)
    for p in profiles:
        row = rowAt(rows, p.time)
        total = ap.integrals(p)
        for q in worst:
            inventory = pc.column(names, row, q)/p.molarMass
            value = total[q] + p.outside.get(q, 0.0)
            worst[q] = max(worst[q], relative(value, inventory))
    ok = all(v <= tol for v in worst.values())
    names = ', '.join(sorted({p.set for p in profiles}))
    done(ok, f'{len(profiles)} profile files ({names}):'
             f' max |sum n\' dx M + outside - inventory|/inventory gas'
             f' {worst["gas"]:.2e}, wall {worst["wall"]:.2e}, sample'
             f' {worst["sample"]:.2e} (limit {tol:g})')


def cmdInventoryHeader(case, tol):
    tol = float(tol)
    names, rows = balanceRows(case)
    worst = 0.0
    profiles = profileFiles(case)
    for p in profiles:
        row = rowAt(rows, p.time)
        for q in ('gas', 'wall', 'sample'):
            worst = max(worst, relative(p.inventory[q],
                                        pc.column(names, row, q)/p.molarMass))
    done(worst <= tol, f'{len(profiles)} headers: inventories vs the balance'
                       f' file, max relative {worst:.2e} (limit {tol:g})')


def cmdOutside(case, binSet, lastWidth):
    lastWidth = float(lastWidth)
    profiles = profileFiles(case, [binSet])
    ok = True
    for p in profiles:
        lo, hi = p.column('x_lo'), p.column('x_hi')
        ok = ok and p.outside['gas'] > 0 \
            and abs(hi[-1] - lo[-1] - lastWidth) <= 1e-9*lastWidth
    p = profiles[-1]
    lo, hi = p.column('x_lo'), p.column('x_hi')
    done(ok, f'set {binSet}: {len(profiles)} files, bins [{lo[0]:g},'
             f' {hi[-1]:g}] m, last bin {hi[-1] - lo[-1]:g} m, outside'
             f' gas {p.outside["gas"]:.4g} mol, wall {p.outside["wall"]:.4g}'
             f' mol at t = {p.time:g}')


def cmdCompareProfiles(caseA, caseB, tol, Ttol='1e-9', xtol='1e-12'):
    tol, Ttol, xtol = float(tol), float(Ttol), float(xtol)
    A = {(os.path.basename(os.path.dirname(p)), os.path.basename(p)): p
         for p in ap.files(caseA)}
    B = {(os.path.basename(os.path.dirname(p)), os.path.basename(p)): p
         for p in ap.files(caseB)}
    if not A or set(A) != set(B):
        done(False, f'different profile files: {sorted(set(A) ^ set(B))}')
    worst = {c: (0.0, None) for c in ('gas', 'wall', 'sample', 'observable')}
    onsetWorst = {'T': 0.0, 'x': 0.0, 'peak': 0.0, 'threshold': 0.0}
    problems = []
    nOnsets = 0

    def ratio(d, scale):
        return 0.0 if d == 0 else (d/scale if scale else math.inf)

    for key in sorted(A):
        a, b = ap.read(A[key]), ap.read(B[key])
        name = '/'.join(key)
        for c in worst:
            va, vb = a.column(c), b.column(c)
            scale = max(abs(x) for x in va)
            d = max(abs(x - y) for x, y in zip(va, vb))
            r = ratio(d, scale)
            if r > worst[c][0]:
                worst[c] = (r, name)
        scale = max(abs(x) for x in a.column('observable'))
        for mode in sorted(set(a.tdep) | set(b.tdep)):
            ta, tb = a.tdep.get(mode, {}), b.tdep.get(mode, {})
            nOnsets += 1
            if ta.get('status') != tb.get('status') \
                    or ta.get('xPeak') != tb.get('xPeak'):
                problems.append(f'{name} {mode}: status/xPeak'
                                f' {ta.get("status")}/{ta.get("xPeak")!r} vs'
                                f' {tb.get("status")}/{tb.get("xPeak")!r}')
                continue
            for q in ('T', 'x', 'peak', 'threshold'):
                va, vb = ta.get(q), tb.get(q)
                if va is None or vb is None:
                    if va != vb:
                        problems.append(f'{name} {mode}: {q} {va} vs {vb}')
                    continue
                d = abs(va - vb)
                onsetWorst[q] = max(onsetWorst[q],
                                    d if q in ('T', 'x') else ratio(d, scale))
    ok = (all(v[0] <= tol for v in worst.values()) and not problems
          and onsetWorst['T'] <= Ttol and onsetWorst['x'] <= xtol
          and onsetWorst['peak'] <= tol and onsetWorst['threshold'] <= tol)
    text = ', '.join(f'{c} {v[0]:.2e}' for c, v in worst.items())
    done(ok, f'{len(A)} profile files: max |B - A|/max|A| per column {text}'
             f' (limit {tol:g}); {nOnsets} onset searches (status and xPeak'
             f' equal): max |T_dep B - T_dep A| {onsetWorst["T"]:.2e} K'
             f' (limit {Ttol:g} K), |x_dep B - x_dep A|'
             f' {onsetWorst["x"]:.2e} m (limit {xtol:g} m), peak'
             f' {onsetWorst["peak"]:.2e} and threshold'
             f' {onsetWorst["threshold"]:.2e} of max|observable|'
             + (f'; problems: {problems}' if problems else ''))


def cmdOverlap(case, nativeSet, cutSet, tol):
    tol = float(tol)
    natives = {(p.directory, p.pair): p
               for p in profileFiles(case, [nativeSet])}
    cuts = profileFiles(case, [cutSet])
    columns = ('gas', 'wall', 'sample', 'observable')
    worst = {c: 0.0 for c in columns + ('Twall', 'wallArea')}
    nCut, nBins = 0, 0
    for q in cuts:
        p = natives.get((q.directory, q.pair))
        if p is None:
            done(False, f'no {nativeSet} file for {q.path}')
        lo, hi = p.column('x_lo'), p.column('x_hi')
        area, T = p.column('wallArea'), p.column('Twall')
        native = {c: p.column(c) for c in columns}
        scale = {c: max(abs(v) for v in native[c]) for c in columns}
        areaScale = max(area)
        cut = {c: q.column(c) for c in columns + ('Twall', 'wallArea')}
        qlo, qhi = q.column('x_lo'), q.column('x_hi')
        for b, (a0, a1) in enumerate(zip(qlo, qhi)):
            dx = a1 - a0
            pieces = [(i, min(hi[i], a1) - max(lo[i], a0))
                      for i in range(len(lo))]
            pieces = [(i, length) for i, length in pieces if length > 0]
            if any(length < (hi[i] - lo[i])*(1 - 1e-12)
                   for i, length in pieces):
                nCut += 1
            nBins += 1
            for c in columns:
                expected = math.fsum(native[c][i]*length
                                     for i, length in pieces)/dx
                d = abs(cut[c][b] - expected)
                r = 0.0 if d == 0 else (d/scale[c] if scale[c] else math.inf)
                worst[c] = max(worst[c], r)
            wallArea = math.fsum(area[i]*length for i, length in pieces)
            expectedT = (math.fsum(area[i]*length*T[i]
                                   for i, length in pieces)/wallArea
                         if wallArea > 0 else 0.0)
            worst['Twall'] = max(worst['Twall'],
                                 relative(cut['Twall'][b], expectedT)
                                 if expectedT else abs(cut['Twall'][b]))
            d = abs(cut['wallArea'][b] - wallArea/dx)
            worst['wallArea'] = max(worst['wallArea'],
                                    0.0 if d == 0 else d/areaScale)
    ok = bool(cuts) and nCut > 0 and all(v <= tol for v in worst.values())
    text = ', '.join(f'{c} {v:.2e}' for c, v in worst.items())
    done(ok, f'{len(cuts)} {cutSet} files, {nBins} bins ({nCut} cut a'
             f' {nativeSet} bin): max difference from the exact overlap'
             f' {text} (limit {tol:g})')


def timeDirectories(directory):
    out = set()
    if not os.path.isdir(directory):
        return out
    for name in os.listdir(directory):
        try:
            float(name)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(directory, name)):
            out.add(name)
    return out


def cmdTimes(case, dt, writeTimes, profileTimes, start=None):
    dt = float(dt)
    writes = [w for w in writeTimes.split(',') if w]
    requested = sorted(float(t) for t in profileTimes.split(',') if t)
    startTime = float(start) if start is not None else None

    profiles = profileFiles(case)
    directories = timeDirectories(os.path.join(case, 'postProcessing',
                                               'axialProfiles'))
    byDirectory = {}
    for p in profiles:
        byDirectory.setdefault(p.directory, []).append(p)

    problems = []
    names, rows = balanceRows(case)
    tLast = max(float(k) for k in rows)
    notReached = [t for t in requested if t - 0.5*dt > tLast]
    # every requested time (after the start, up to the end) is covered by
    # one directory whose exact time is the nearest step end (within dt/2,
    # ties to the earlier one)
    for t in requested:
        if startTime is not None and t - 0.5*dt <= startTime:
            continue
        if t in notReached:
            continue
        hits = [d for d, ps in byDirectory.items()
                if -0.5*dt < ps[0].time - t <= 0.5*dt + 1e-12*dt]
        if len(hits) != 1:
            problems.append(f'profile time {t:g}: {len(hits)} directories')
        else:
            reason = byDirectory[hits[0]][0].reason
            if f'profile time {t:g}' not in reason:
                problems.append(f'profile time {t:g}: reason "{reason}"')
    for w in writes:
        if w not in directories:
            problems.append(f'write time {w}: no profile directory')
        elif not all(p.reason.startswith('write time')
                     for p in byDirectory.get(w, [])):
            problems.append(f'write time {w}: reason'
                            f' "{byDirectory[w][0].reason}"')
    # no other directory
    expected = set(writes)
    for d, ps in byDirectory.items():
        reason = ps[0].reason
        if 'profile time' not in reason and d not in expected:
            problems.append(f'unexpected profile directory {d}')
    # the case's (and processors') time directories are the write times only
    allowed = set(writes) | {'0'} | ({start} if start is not None else set())
    caseDirs = timeDirectories(case)
    extra = caseDirs - allowed
    procExtra = set()
    for name in os.listdir(case):
        if name.startswith('processor'):
            procExtra |= timeDirectories(os.path.join(case, name)) - allowed
    if extra or procExtra:
        problems.append(f'time directories besides the write times:'
                        f' {sorted(extra | procExtra)}')
    if startTime is not None:
        early = [d for d, ps in byDirectory.items()
                 if ps[0].time <= startTime + 0.5*dt]
        if early:
            problems.append(f'written at or before the restart: {early}')
    done(not problems,
         f'profile directories {sorted(directories, key=float)}; case time'
         f' directories {sorted(caseDirs, key=float)}'
         + (f'; not reached (after the end): {notReached}'
            if notReached else '')
         + (f'; problems: {problems}' if problems else ''))


def cmdTdep(case, binSet, mode, expected, tol):
    expected, tol = float(expected), float(tol)
    profiles = profileFiles(case, [binSet])
    worst, worstAgain, text = 0.0, 0.0, ''
    ok = True
    for p in profiles:
        written = p.tdep.get(mode, {})
        if written.get('status') != 'found' or written.get('T') is None:
            ok = False
            text = f'status {written.get("status")}'
            continue
        again = ap.onset(p, mode == 'fromPeak', written.get('searchFrom'))
        worst = max(worst, abs(written['T'] - expected))
        worstAgain = max(worstAgain,
                         abs(again.get('T', math.inf) - written['T']))
        text = (f'T_dep {mode} {written["T"]!r} K at x {written["x"]!r} m')
    ok = ok and worst <= tol and worstAgain <= tol
    done(ok, f'set {binSet}: {text}, analytic {expected!r} K: difference'
             f' {worst:.2e} K, recomputed from the columns {worstAgain:.2e} K'
             f' (limit {tol:g} K)')


# the synthetic ramp: 100 cells of 1 mm, ramp on cells 30..80, spike in 10
RAMP_CELLS = 100
RAMP_DX = 1e-3
RAMP_FIRST, RAMP_PEAK, RAMP_SPIKE = 30, 80, 10
RAMP_SLOPE = 1.2e-3


def rampValues(spike):
    """rho Y of the cells [kg/m3], analytic"""
    x = [(i + 0.5)*RAMP_DX for i in range(RAMP_CELLS)]
    x0 = x[RAMP_FIRST]
    r = [RAMP_SLOPE*(xi - x0) if RAMP_FIRST <= i <= RAMP_PEAK else 0.0
         for i, xi in enumerate(x)]
    if spike:
        r[RAMP_SPIKE] = 0.3*r[RAMP_PEAK]
    return x, r


def cmdRampField(case, spike):
    rho = pc.readInternal(os.path.join(case, '0', 'rho'), RAMP_CELLS)
    if len(rho) != RAMP_CELLS:
        done(False, f'0/rho has {len(rho)} cells, not {RAMP_CELLS}')
    _, r = rampValues(spike == '1')
    Y = [ri/rhoi for ri, rhoi in zip(r, rho)]
    path = os.path.join(case, '0', 'Y_PbI2_g')
    text = open(path).read()
    values = '\n'.join(repr(y) for y in Y)
    new, n = re.subn(r'^internalField .*?;$',
                     f'internalField nonuniform List<scalar> {len(Y)}\n(\n'
                     f'{values}\n);', text, count=1, flags=re.M)
    if n != 1:
        done(False, f'no internalField in {path}')
    open(path, 'w').write(new)
    done(True, f'{path}: ramp on cells {RAMP_FIRST}-{RAMP_PEAK}'
               + (f', spike in cell {RAMP_SPIKE}' if spike == '1' else ''))


def rampOnset(values, centres, fromPeak, searchFrom=None, fraction=0.01):
    n = len(values)
    first = 0
    if not fromPeak and searchFrom is not None:
        while first < n and centres[first] < searchFrom:
            first += 1
    peakBin = max(range(first, n), key=lambda b: (values[b], -b))
    threshold = fraction*values[peakBin]
    if fromPeak:
        b = peakBin
        while b > first and values[b - 1] >= threshold:
            b -= 1
    else:
        b = first
        while values[b] < threshold:
            b += 1
    w = (threshold - values[b - 1])/(values[b] - values[b - 1])
    return centres[b - 1] + w*(centres[b] - centres[b - 1])


def cmdRampExpected(mode, observable, binSet, spike, searchFrom=None):
    x, r = rampValues(spike == '1')
    if observable == 'wall' and spike == '1':
        r[RAMP_SPIKE] = 0.0            # f = 0 above Tdep: gas only
    if binSet == 'coarse':
        k = 5
        centres = [(b + 0.5)*k*RAMP_DX for b in range(RAMP_CELLS//k)]
        values = [sum(r[b*k:(b + 1)*k])/k for b in range(RAMP_CELLS//k)]
    else:
        centres, values = x, r
    sf = float(searchFrom) if searchFrom is not None else None
    xDep = rampOnset(values, centres, mode == 'fromPeak', sf)
    print(repr(1000.0 - 6000.0*xDep))


def fieldPatch(path, patch):
    """values of a patch of an ascii field (17 digits)"""
    return pc.readPatch(path, patch)


def cmdMDepFactor(case, factor):
    factor = float(factor)
    times = pc.timeDirs(case)
    worst, n = 0.0, 0
    ok = bool(times)
    for t in times:
        mw = fieldPatch(os.path.join(case, t, 'mw_PbI2_s'), 'WALL')
        md = fieldPatch(os.path.join(case, t, 'mDep_PbI2_s'), 'WALL')
        if len(mw) != len(md) or not mw:
            ok = False
            continue
        for a, b in zip(mw, md):
            n += 1
            if b != factor*a:
                ok = False
                worst = max(worst, abs(b - factor*a)/max(abs(b), 1e-300))
        ok = ok and max(mw) > 0
    done(ok, f'{len(times)} written times, {n} WALL faces: mDep_PbI2_s ='
             f' {factor:g} x mw_PbI2_s exactly (max relative difference'
             f' {worst:.2e})')


def readFunctionObject(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            rows.append((parts[0], [float(v) for v in parts[1:]]))
    return rows


def cmdAreaIntegral(case, dat, col, tol, control=None):
    col, tol = int(col), float(tol)
    names, rows = balanceRows(case)
    fo = readFunctionObject(dat)
    worst, n = 0.0, 0
    for t, values in fo:
        row = rowAt(rows, float(t)) if any(float(k) == float(t)
                                            for k in rows) else None
        if row is None:
            # function objects print the rounded time name; match by value
            candidates = [k for k in rows if abs(float(k) - float(t))
                          <= 1e-9*max(1.0, abs(float(t)))]
            if len(candidates) != 1:
                done(False, f'no balance row for the function-object time {t}')
            row = rows[candidates[0]]
        wall = pc.column(names, row, 'wall')
        worst = max(worst, relative(values[col], wall))
        n += 1
    if control is None:
        done(n > 0 and worst <= tol,
             f'{n} steps: max |areaIntegrate - wall inventory|/wall'
             f' {worst:.2e} (limit {tol:g})')
    else:
        done(n > 0 and worst >= float(control),
             f'{n} steps: control, max |areaIntegrate - wall inventory|/wall'
             f' {worst:.3g} (at least {float(control):g}: not the inventory)')


def cmdExchange(case, dt, tol):
    dt, tol = float(dt), float(tol)
    names, rows = balanceRows(case)
    times = ['0'] + pc.timeDirs(case)
    nCells = 10
    rho = pc.readInternal(os.path.join(case, '0', 'rho'), nCells)
    V = 1e-3*57.46e-6*1e-3
    worstCell, worstSum, n = 0.0, 0.0, 0
    ok = len(times) > 1
    previousWall = 0.0
    for a, b in zip(times[:-1], times[1:]):
        Y0 = pc.readInternal(os.path.join(case, a, 'Y_PbI2_g'), nCells)
        Y1 = pc.readInternal(os.path.join(case, b, 'Y_PbI2_g'), nCells)
        exch = pc.readInternal(os.path.join(case, b, 'exch_PbI2_g'), nCells)
        expected = [r*(y1 - y0)/dt for r, y0, y1 in zip(rho, Y0, Y1)]
        scale = max(abs(e) for e in exch)
        worstCell = max(worstCell, max(abs(e - x) for e, x in
                                       zip(exch, expected))/scale)
        # the balance row of this written time (names are rounded)
        key = min(rows, key=lambda k: abs(float(k) - float(b)))
        wall = pc.column(names, rows[key], 'wall')
        total = math.fsum(e*V*dt for e in exch)
        worstSum = max(worstSum, relative(-total, wall - previousWall))
        previousWall = wall
        n += 1
    ok = ok and worstCell <= tol and worstSum <= tol
    done(ok, f'{n} steps: max |exch - rho dY/dt|/max|exch| {worstCell:.2e},'
             f' max |sum exch V dt + d(wall)|/d(wall) {worstSum:.2e}'
             f' (limit {tol:g})')


def cmdUniformPatch(path, patch, mode=None):
    data = open(path, 'rb').read()
    m = re.search(rb'\n\s*' + patch.encode() + rb'\s*\{(.*?)\}', data, re.S)
    if not m:
        done(False, f'no patch {patch} in {path}')
    if mode == 'noValue':
        t = re.search(rb'type\s+(\w+)\s*;', m.group(1))
        has = re.search(rb'\bvalue\s', m.group(1)) is not None
        done(not has, f'{os.path.basename(path)} on {patch}'
                      f' ({t.group(1).decode() if t else "?"}): '
                      + ('a value entry' if has else 'no value written'))
    v = re.search(rb'value\s+uniform\s+([-+0-9.eE]+)\s*;', m.group(1))
    if not v:
        done(False, f'{path}: {patch} holds no "value uniform"')
    text = v.group(1).decode()
    digits = len(re.sub(r'[eE].*$', '', text).replace('-', '')
                 .replace('.', '').lstrip('0'))
    done(digits <= 6, f'{os.path.basename(path)} on {patch}: "value uniform'
                      f' {text}" ({digits} significant digits)')


commands = {
    'integrals': cmdIntegrals,
    'inventoryHeader': cmdInventoryHeader,
    'outside': cmdOutside,
    'compareProfiles': cmdCompareProfiles,
    'overlap': cmdOverlap,
    'times': cmdTimes,
    'tdep': cmdTdep,
    'rampField': cmdRampField,
    'rampExpected': cmdRampExpected,
    'mDepFactor': cmdMDepFactor,
    'areaIntegral': cmdAreaIntegral,
    'exchange': cmdExchange,
    'uniformPatch': cmdUniformPatch,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in commands:
        print(__doc__)
        sys.exit(2)
    commands[sys.argv[1]](*sys.argv[2:])
