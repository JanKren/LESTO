#!/usr/bin/env python3
"""
tests/hks/checks.py

PURPOSE

Numerical checks of the M5 tests (tests/hks/Allrun): the Hertz-Knudsen-
Schrage law, inventory samples and their removal.  Every subcommand prints
one summary line and exits with 0 (criterion met) or 1.  The model values
are computed here independently of the solver: the p_eq table is read and
interpolated as vapourPressureTable.H specifies (log10 p linear in 1/T per
phase, the minimum over the phases), G from the formula of plan E2.

Box cases (tests/phaseChange/closedBox at 700 K: 10 cells of height h =
57.46 um, one WALL face each, A/V = 1/h, at rest): with
lambda = (1/h) G R T/M and Y_eq = p_eq/beta, beta = rho R T/M, a step of
an IMPLICIT element gives Y1 = (Y0 + lambda dt Y_eq)/(1 + lambda dt).

  model     <case> [halfCell]
      prints lambda dt and Y_eq of the box (with the half cell in series:
      G' = G/(1 + G beta d/(rho D)), d = h/2)
  approach  <case> <Y0> <tol> [halfCell] [deposit=<tol>]
      every written time n (1, 2, ...): every cell of Y_PbI2_g equals
      Y_eq + (Y0 - Y_eq)/(1 + lambda dt)^n (deposition from a bare wall in
      the first step, from a wall with a deposit afterwards) to tol, and the
      deposit m'' = (Y0 - Y^n) rho h to 1e3 tol (or the tolerance given)
  unchanged <case> <Y0> <tol>
      bare wall below Y_eq: at every written time mw_PbI2_s, Y_PbI2_s and
      the realised exchange exch_PbI2_g are exactly 0 (deposit only; the
      case writes exch_), and Y_PbI2_g equals Y0 to tol (relative).  Y is
      not bit for bit Y0: on a uniform field the normalisation factor of
      OpenFOAM's residual collapses to round-off, so the linear solver
      iterates once per step and changes Y in its last bits (the maximum
      deviation is printed)
  evaporate <case> <time> <deposit [kg/m2]> <tol>
      one step from Y = 0 with a wall deposit m'' (a restart): IMPLICIT when
      the deposit suffices, Y1 = lambda dt Y_eq/(1 + lambda dt), m''_1 =
      m'' - Y1 rho h to tol; otherwise CAPPED, Y1 = m''/(rho h) to tol and
      the deposit, Y_PbI2_s and mw_PbI2_s exactly 0; prints the regime
  closureStep <dat> <tol>
      the closure of every row of a balance file relative to the
      inventory at most tol, and the closure of the last row, relative to
      the reference max(|initial + released|, |held|, inflow)

Channel and pipe:

  removal <dat> <removeTime> <dt> <tStart> <nSteps>
      on the exact clock t_0 = tStart, t_k+1 = t_k + dt: the removal belongs
      to the first step whose start t_k satisfies t_k >= removeTime - dt/2;
      the row of its end (time t_k + dt) is the first with removed > 0,
      removed equals the sample column of the row t_k (the reservoir at the
      start of that step) to 1e-15, and the sample column is 0 exactly from
      that row on
  decreasing <dat> <column> [<until time>]
      the column decreases strictly from row to row (up to the given time)
  logCounts <log>
      prints the numbers of steps with significant regime changes in the
      deciding corrector and in the final predictor, with projections
      (any, significant, a total above the tolerance), guard passes and
      complementarity violations (any, significant, a total above the
      tolerance); exit 0 always (the Allrun decides)
  totalLeft <log> <dat> <guardTolerance>
      the projected mass and the missed mass of the complementarity
      violations of every step (exchange lines) against guardTolerance
      times the inventory at the start of the step, recomputed from the
      balance file (a fresh run of one pair; the inflow scale of the
      solver's tolerance is left out, which only lowers the limit); exit 0
      when no step's complementarity mass exceeds it and the log flags
      none (review of M5, round 3: the guard accepted a sum of
      sub-threshold misses above it)
  layoutSample <uniform/phaseChangeLayout> <sample> <key> <value>
      the samples record of the layout holds the entry (a quoted string)
  sameRegimes [steps=<n>] <log> <log> ...
      the regimes of every step (the counts of the exchange line, wall and
      sample elements; of the first n steps only with steps=<n>) are the
      same in all logs; prints the number of steps, how many of them are
      equal, and the sample regimes of the first steps
  inventories <datA> <datB> <reference [kg]> <tol>
      gas, wall and sample of the rows of equal time differ by at most tol
      x reference (a fixed reference: the sample column ends at 0)
  complementarity <case> <limit> box <x0 y0 z0 x1 y1 z1>
  complementarity <case> <limit> cylinder <x0 y0 z0 x1 y1 z1 r>
      the bounded HKS law of the inventory sample, recomputed from the
      written fields independently of the solver's own count (review of
      M5, round 2).  The run writes every step and exch_PbI2_g; 0/V, 0/Cx,
      0/Cy and 0/Cz come from postProcess -func writeCellVolumes and
      writeCellCentres.  The sample cells are those whose centre boxToCell
      (inclusive) or cylinderToCell (strict, its arithmetic) selects.  At
      every written time t_k and in every sample cell, the realised flow
      exch_PbI2_g V must equal the admissible flow min(max(q, 0), upper)
      with q = a_s V G (p_eq(T) - rho R T Y/M) at the written Y (G of the
      sample's kineticScale and Ce) and upper = Cs V/dt, Cs of the previous
      written time (the fill n0 M/V_Z before the first step of a fresh
      start).  An element misses |admissible - realised| dt; exit 0 when no
      element misses more than limit x n0 M.  Also prints the steps with
      misses, the largest miss, and the cells that took their bound
      (exhausted, or no flow) in supersaturated or undersaturated gas
  wallLaw <case> <limit> [<patch>]
      the bounded HKS law of the wall elements (the faces of the interface
      patch, default WALL; faceCells, areaMultiplier 1, no half cell),
      recomputed from the written fields (review of M5, round 4: the guard
      re-solve re-assembled the SuperBee matrix and let a bare wall element
      cycle between IMPLICIT and NONE).  The run writes every step.  At
      every written time t_k and for every element e (face f, owner cell
      c): the realised flow q_e = (mw_e(t_k-1) - mw_e(t_k)) A_f/dt from the
      deposits must equal the admissible flow min(a - g Y_c, upper) with a =
      A_f G p_eq(T_e), g = A_f G rho_c R T_c/M, G of HKSCoeffs at T_e (the
      patch value of T with interfaceTemperature wall, T_c otherwise) and
      upper = mw_e(t_k-1) A_f/dt (a bare wall only deposits).  The face
      areas come from constant/polyMesh (points, faces, owner, boundary;
      ascii or binary).  An element misses |admissible - q| dt; exit 0 when
      no element misses more than limit x the initial inventory of the
      balance file (at least that of the samples, n0 M).  Also prints the
      steps with misses and the largest miss

Balance-file columns are addressed by the names of its '#' header line.
"""

import math
import os
import re
import struct
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'phaseChange'))
from checks import readBalance, column, reference, readInternal, \
    readScalarField, timeDirs, foamFormat, parseList   # noqa: E402

R = 8.314462618
H_BOX = 57.46e-6


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def dictValue(path, key, default=None, cast=float):
    text = open(path).read()
    m = re.search(r'(?<![\w.])' + re.escape(key) + r'\s+([^;\s]+)\s*;', text)
    if not m:
        if default is None:
            raise SystemExit(f'No {key} in {path}')
        return default
    return cast(m.group(1).strip('"'))


def readTable(path, phases):
    """Nodes and log10 p of the phases of a csv table (vapourPressureTable.H
    format; values in log10 Pa)."""
    header, rows = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if not rows:
                    header = [c.strip() for c in line[1:].split(',')]
                continue
            rows.append([float(x) for x in line.split(',')])
    cols = [header.index(p) for p in phases]
    return [r[0] for r in rows], [[r[c] for r in rows] for c in cols]


def pEq(T, nodes, columns):
    """p_eq [Pa]: per phase log10 p linear in 1/T, the minimum."""
    n = len(nodes)
    for k, t in enumerate(nodes):
        if t == T:
            return min(10.0**c[k] for c in columns)
    i = 0
    if T >= nodes[-1]:
        i = n - 2
    elif T > nodes[0]:
        while nodes[i + 1] <= T:
            i += 1
    w = (1/T - 1/nodes[i])/(1/nodes[i + 1] - 1/nodes[i])
    return min(10.0**((1 - w)*c[i] + w*c[i + 1]) for c in columns)


def boxModel(case, halfCell=False):
    """lambda dt, Y_eq, rho, lambda of the box (uniform carrier)."""
    thermo = open(os.path.join(case, 'constant',
                               'thermochemistryProperties')).read()
    table = re.search(r'file\s+"<constant>/([^"]+)"', thermo).group(1)
    phases = re.findall(r'"([^"]+)"',
                        re.search(r'phases\s*\((.*?)\)\s*;', thermo).group(1))
    nodes, columns = readTable(os.path.join(case, 'constant', table), phases)
    species = os.path.join(case, 'constant', 'speciesTransportProperties')
    M = dictValue(species, 'molarMass')
    D = dictValue(species, 'D')
    dt = dictValue(os.path.join(case, 'system', 'controlDict'), 'deltaT')
    T = readInternal(os.path.join(case, '0', 'T'), 1)[0]
    rho = readInternal(os.path.join(case, '0', 'rho'), 1)[0]

    def coefficient(key, default):
        m = re.search(r'HKSCoeffs\s*\{[^}]*?\b' + key + r'\s+([^;\s]+);',
                      thermo)
        return float(m.group(1)) if m else default

    sigma = coefficient('accommodation', 1.0)
    Ce = coefficient('Ce', 1.0)
    scale = coefficient('kineticScale', 1.0)
    G = scale*Ce*2*sigma/(2 - sigma)*math.sqrt(M/(2*math.pi*R*T))
    beta = rho*R*T/M
    if halfCell:
        d = H_BOX/2
        G = G/(1 + G*beta*d/(rho*D))
    lam = G*beta/rho/H_BOX
    Yeq = pEq(T, nodes, columns)/beta
    return lam*dt, Yeq, rho, lam


def cmdModel(case, *options):
    lamDt, Yeq, rho, lam = boxModel(case, 'halfCell' in options)
    done(True, f'lambda = {lam:.6e} 1/s, lambda dt = {lamDt:.6f}, '
               f'Y_eq = {Yeq!r}, rho = {rho!r}')


def cmdApproach(case, Y0, tol, *options):
    Y0, tol = float(Y0), float(tol)
    depositTol = 1e3*tol
    for o in options:
        if o.startswith('deposit='):
            depositTol = float(o.split('=', 1)[1])
    lamDt, Yeq, rho, _ = boxModel(case, 'halfCell' in options)
    times = timeDirs(case)
    worst, worstDeposit = 0.0, 0.0
    for n, t in enumerate(times, 1):
        Y = readInternal(os.path.join(case, t, 'Y_PbI2_g'), 10)
        want = Yeq + (Y0 - Yeq)/(1 + lamDt)**n
        worst = max(worst, max(abs(y - want)/abs(want) for y in Y))
        mw = readScalarField(os.path.join(case, t, 'mw_PbI2_s'))[1]['WALL']
        deposit = (Y0 - want)*rho*H_BOX
        worstDeposit = max(worstDeposit,
                           max(abs(m - deposit)/abs(deposit) for m in mw))
    ok = times and worst <= tol and worstDeposit <= depositTol
    done(ok, f'Y^n = Y_eq + (Y0 - Y_eq)/(1 + lambda dt)^n, lambda dt = '
             f'{lamDt:.6g}, Y_eq = {Yeq:.6g}, Y0 = {Y0:g}: max relative '
             f'error {worst:.2e} over {len(times)} steps (limit {tol:g}); '
             f'deposit (Y0 - Y^n) rho h to {worstDeposit:.1e} (limit '
             f'{depositTol:g})')


def cmdUnchanged(case, Y0, tol):
    Y0, tol = float(Y0), float(tol)
    times = timeDirs(case)
    exact = bool(times)
    worst, nSame = 0.0, 0
    for t in times:
        Y = readInternal(os.path.join(case, t, 'Y_PbI2_g'), 10)
        mw = readScalarField(os.path.join(case, t, 'mw_PbI2_s'))[1]['WALL']
        Ys = readInternal(os.path.join(case, t, 'Y_PbI2_s'), 10)
        exch = readInternal(os.path.join(case, t, 'exch_PbI2_g'), 10)
        exact = exact and all(m == 0 for m in mw) \
            and all(y == 0 for y in Ys) and all(q == 0 for q in exch)
        worst = max(worst, max(abs(y - Y0)/Y0 for y in Y))
        nSame += sum(y == Y0 for y in Y)
    _, Yeq, _, _ = boxModel(case)
    ok = exact and worst <= tol
    done(ok, f'bare wall with Y0 = {Y0:g} < Y_eq = {Yeq:.6g}: exch_PbI2_g, '
             f'mw_PbI2_s and Y_PbI2_s exactly 0 at all {len(times)} written '
             f'times: {exact}; Y_PbI2_g = Y0 to {worst:.1e} (limit {tol:g}; '
             f'{nSame} of {10*len(times)} values bit for bit, the rest '
             f'changed by the linear solver)')


def cmdEvaporate(case, t, deposit, tol):
    deposit, tol = float(deposit), float(tol)
    lamDt, Yeq, rho, _ = boxModel(case)
    Y = readInternal(os.path.join(case, t, 'Y_PbI2_g'), 10)
    mw = readScalarField(os.path.join(case, t, 'mw_PbI2_s'))[1]['WALL']
    Ys = readInternal(os.path.join(case, t, 'Y_PbI2_s'), 10)
    Yi = lamDt*Yeq/(1 + lamDt)
    if Yi*rho*H_BOX <= deposit:
        regime = 'IMPLICIT'
        wantY, wantW = Yi, deposit - Yi*rho*H_BOX
        errY = max(abs(y - wantY)/wantY for y in Y)
        errW = max(abs(m - wantW)/wantW for m in mw)
        ok = errY <= tol and errW <= 1e3*tol
        extra = (f'Y1 = lambda dt Y_eq/(1 + lambda dt) = {wantY:.10g} to '
                 f'{errY:.1e}, deposit m\'\' - Y1 rho h to {errW:.1e}')
    else:
        regime = 'CAPPED'
        wantY = deposit/(rho*H_BOX)
        errY = max(abs(y - wantY)/wantY for y in Y)
        exact = all(m == 0 for m in mw) and all(y == 0 for y in Ys)
        ok = errY <= tol and exact
        extra = (f'Y1 = m\'\'/(rho h) = {wantY:.10g} to {errY:.1e}; deposit '
                 f'mw_PbI2_s and Y_PbI2_s exactly 0: {exact}')
    done(ok, f'{regime} (implicit evaporation would take {Yi*rho*H_BOX:.4g} '
             f'kg/m2 of the deposit {deposit:g} kg/m2): {extra} '
             f'(limit {tol:g})')


def cmdClosureStep(dat, tol):
    names, rows = readBalance(dat)
    tol = float(tol)
    worst = 0.0
    for r in rows:
        held = abs(column(names, r, 'gas') + column(names, r, 'wall')
                   + column(names, r, 'sample'))
        worst = max(worst, abs(column(names, r, 'closure'))/max(held, 1e-300))
    last = abs(column(names, rows[-1], 'closure'))/reference(names, rows[-1])
    done(worst <= tol, f'closure relative to the inventory: max {worst:.2e} '
                       f'over {len(rows)} rows, last {last:.2e} (limit '
                       f'{tol:g})')


def cmdRemoval(dat, removeTime, dt, tStart, nSteps):
    removeTime, dt, t = float(removeTime), float(dt), float(tStart)
    names, rows = readBalance(dat)
    clock = []
    for _ in range(int(nSteps) + 1):
        clock.append(t)
        t = t + dt
    k = next(i for i, tk in enumerate(clock) if tk >= removeTime - 0.5*dt)
    byTime = {float(r[0]): r for r in rows}
    end = clock[k] + dt
    if end not in byTime:
        done(False, f'no balance row at the exact end {end!r} of the removal '
                    f'step')
    first = next(i for i, r in enumerate(rows)
                 if column(names, r, 'removed') > 0)
    removed = column(names, rows[first], 'removed')
    before = byTime.get(clock[k])
    sampleBefore = column(names, before, 'sample') if before else None
    ok = float(rows[first][0]) == end and before is not None \
        and abs(removed - sampleBefore) <= 1e-15*abs(sampleBefore) \
        and all(column(names, r, 'sample') == 0 for r in rows[first:]) \
        and all(column(names, r, 'removed') == removed for r in rows[first:])
    done(ok, f'removal in the step starting at t_{k} = {clock[k]!r} (the '
             f'first t_k >= removeTime - dt/2 = {removeTime - 0.5*dt!r}): '
             f'first removed > 0 at {rows[first][0]} (expected {end!r}), '
             f'removed {removed!r} kg = the reservoir at the step start '
             f'{sampleBefore!r} kg; sample 0 exactly and removed constant '
             f'in the {len(rows) - first} rows from there')


def cmdDecreasing(dat, name, until=None):
    names, rows = readBalance(dat)
    if until is not None:
        rows = [r for r in rows if float(r[0]) <= float(until) + 1e-12]
    values = [column(names, r, name) for r in rows]
    ok = len(values) > 1 and all(b < a for a, b in zip(values, values[1:]))
    done(ok, f'{name} decreases strictly over {len(values)} rows, from '
             f'{values[0]:.10g} to {values[-1]:.10g} kg')


# the counts of an exchange line; a total above the guard's tolerance is
# named after the mass, '(<kg> kg > tolerance <kg> kg, <n> significant)'
EXCHANGE = re.compile(
    r'regime changes (\d+), \S+ kg, (\d+) significant\)(?:; final '
    r'predictor: regime changes (\d+), \S+ kg, (\d+) significant)?; '
    r'guard passes (\d+) \((\d+) elements\); projected (\d+) \((\S+) kg'
    r'( > tolerance \S+ kg)?, (\d+) significant\); complementarity (\d+) '
    r'\((\S+) kg( > tolerance \S+ kg)?, (\d+) significant\)')


def exchangeLines(log):
    """The exchange lines of a log: dicts of their counts and masses."""
    lines = []
    for line in open(log):
        m = EXCHANGE.search(line)
        if not m:
            continue
        g = m.groups()
        n = [int(g[i]) if g[i] is not None else 0
             for i in (0, 1, 2, 3, 4, 5, 6, 9, 10, 13)]
        lines.append(dict(
            decidingSignificant=n[1], finalPredictorSignificant=n[3],
            guard=n[4], projected=n[6], projectedMass=float(g[7]),
            projectedAbove=g[8] is not None, projectedSignificant=n[7],
            complementarity=n[8], complementarityMass=float(g[11]),
            complementarityAbove=g[12] is not None,
            complementaritySignificant=n[9]))
    return lines


def cmdLogCounts(log):
    counts = dict(steps=0, decidingSignificant=0,
                  finalPredictorSignificant=0, projected=0,
                  projectedSignificant=0, projectedAbove=0, guard=0,
                  complementarity=0, complementarityAbove=0,
                  complementaritySignificant=0)
    for x in exchangeLines(log):
        counts['steps'] += 1
        for k in counts:
            if k != 'steps':
                counts[k] += x[k] > 0
    print(' '.join(f'{k}={v}' for k, v in counts.items()))
    sys.exit(0)


def cmdTotalLeft(log, dat, guardTolerance):
    """The masses the guard left in every step (projected, and the missed
    mass of the complementarity violations) against guardTolerance times
    the inventory at the start of the step, recomputed from the balance
    file (gas + wall + sample of the row before, the initial inventory
    before the first row, plus the release of the step; the convective
    inflow of massTolerance() is not in the file and is left out, which
    only makes the limit smaller).  A fresh run of one pair."""
    guardTolerance = float(guardTolerance)
    lines = exchangeLines(log)
    names, rows = readBalance(dat)
    held = column(names, rows[0], 'initial')
    released = 0.0
    worst = {'projected': (0.0, ''), 'complementarity': (0.0, '')}
    above = {'projected': 0, 'complementarity': 0}
    flagged = {'projected': 0, 'complementarity': 0}
    for x, row in zip(lines, rows):
        releasedNow = column(names, row, 'released')
        limit = guardTolerance*(abs(held) + releasedNow - released)
        for k in worst:
            ratio = x[k + 'Mass']/limit if limit > 0 else 0.0
            if ratio > worst[k][0]:
                worst[k] = (ratio, row[0])
            above[k] += ratio > 1
            flagged[k] += x[k + 'Above']
        held = sum(column(names, row, k) for k in ('gas', 'wall', 'sample'))
        released = releasedNow
    ok = len(lines) == len(rows) > 0 and above['complementarity'] == 0 \
        and flagged['complementarity'] == 0
    done(ok, f'{len(lines)} exchange lines, {len(rows)} balance rows: steps '
             f'whose complementarity mass exceeds guardTolerance x inventory '
             f'{above["complementarity"]} (largest ratio '
             f'{worst["complementarity"][0]:.3g} at t = '
             f'{worst["complementarity"][1] or "-"}), flagged in the log '
             f'{flagged["complementarity"]}; projected mass: '
             f'{above["projected"]} steps above (largest ratio '
             f'{worst["projected"][0]:.3g}), flagged {flagged["projected"]}')


def cmdLayoutSample(path, sample, key, value):
    text = open(path).read()
    m = re.search(r'\bsamples\s*\{(.*?)\n\}', text, re.S)
    block = m.group(1) if m else ''
    s = re.search(r'\b' + re.escape(sample) + r'\s*\{([^}]*)\}', block)
    found = re.search(r'\b' + re.escape(key) + r'\s+"([^"]*)"\s*;',
                      s.group(1)) if s else None
    got = found.group(1) if found else None
    done(got == value, f'layout record of sample {sample}: {key} '
                       f'"{got}" (expected "{value}")')


def regimesOfLog(log):
    """The regime counts of the exchange line of every step."""
    pat = re.compile(r'exchange: implicit (\d+), capped (\d+), none (\d+) '
                     r'elements(?: \(sample: implicit (\d+), capped (\d+), '
                     r'none (\d+)\))?')
    return [m.groups() for m in map(pat.search, open(log)) if m]


def cmdSameRegimes(*args):
    steps = None
    if args and args[0].startswith('steps='):
        steps, args = int(args[0][6:]), args[1:]
    runs = [regimesOfLog(log) for log in args]
    compared = [r[:steps] for r in runs]
    same = all(r == compared[0] for r in compared) and len(compared[0]) > 0 \
        and (steps is None or len(compared[0]) == steps)
    equal = sum(all(r[k] == runs[0][k] for r in runs)
                for k in range(min(len(r) for r in runs)))
    first = '; '.join(f'{g[3]}/{g[4]}/{g[5]}' for g in runs[0][:4])
    which = 'every step' if steps is None else f'the first {steps} steps'
    done(same, f'regimes of {which} equal in {len(args)} runs: {same} '
               f'({[len(r) for r in runs]} steps, {equal} of them equal; '
               f'sample implicit/capped/none in the first steps: {first})')


def cmdInventories(datA, datB, ref, tol):
    ref, tol = float(ref), float(tol)
    namesA, rowsA = readBalance(datA)
    namesB, rowsB = readBalance(datB)
    byTime = {r[0]: r for r in rowsB}
    worst = {k: 0.0 for k in ('gas', 'wall', 'sample')}
    matched = 0
    for rA in rowsA:
        rB = byTime.get(rA[0])
        if rB is None:
            continue
        matched += 1
        for k in worst:
            worst[k] = max(worst[k], abs(column(namesA, rA, k)
                                         - column(namesB, rB, k))/ref)
    ok = matched == len(rowsA) == len(rowsB) and max(worst.values()) <= tol
    done(ok, f'{matched} steps, largest differences relative to '
             f'{ref:.6g} kg: ' + ', '.join(f'{k} {v:.2e}'
                                            for k, v in worst.items())
             + f' (limit {tol:g})')


def stripComments(text):
    """An OpenFOAM dictionary without its /* */ and // comments."""
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    return re.sub(r'//[^\n]*', ' ', text)


def internalField(path, nCells):
    """Internal values of a scalar field, a 'uniform' entry expanded."""
    values = readScalarField(path)[0]
    return values*nCells if len(values) == 1 and nCells > 1 else values


def selectCells(case, shape, coords):
    """Cells whose centre boxToCell (inclusive) or cylinderToCell (strict,
    with its arithmetic) selects; centres from 0/Cx, 0/Cy, 0/Cz."""
    C = [readScalarField(os.path.join(case, '0', 'C' + x))[0] for x in 'xyz']
    v = [float(x) for x in coords]
    cells = []
    if shape == 'box':
        lo, hi = v[0:3], v[3:6]
        for c in range(len(C[0])):
            if all(lo[i] <= C[i][c] <= hi[i] for i in range(3)):
                cells.append(c)
    elif shape == 'cylinder':
        p1, p2, r = v[0:3], v[3:6], v[6]
        axis = [p2[i] - p1[i] for i in range(3)]
        magAxis2 = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]
        orad2 = r*r
        for c in range(len(C[0])):
            d = [C[i][c] - p1[i] for i in range(3)]
            magD = d[0]*axis[0] + d[1]*axis[1] + d[2]*axis[2]
            if 0 < magD < magAxis2:
                d2 = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) - magD*magD/magAxis2
                if d2 < orad2:
                    cells.append(c)
    else:
        raise SystemExit(f'Unknown selection {shape}')
    return cells


def cmdComplementarity(case, limit, shape, *coords):
    limit = float(limit)
    thermo = stripComments(open(os.path.join(
        case, 'constant', 'thermochemistryProperties')).read())
    table = re.search(r'file\s+"<constant>/([^"]+)"', thermo).group(1)
    phases = re.findall(r'"([^"]+)"', re.search(r'phases\s*\((.*?)\)\s*;',
                                                thermo).group(1))
    nodes, columns = readTable(os.path.join(case, 'constant', table), phases)
    hks = re.search(r'HKSCoeffs\s*\{([^}]*)\}', thermo)
    m = re.search(r'\baccommodation\s+([^;\s]+);', hks.group(1)) \
        if hks else None
    sigma = float(m.group(1)) if m else 1.0
    # the dictionary of the (first) inventory sample: from the brace that
    # opens it to the matching one
    samples = thermo[re.search(r'\bsamples\s*\{', thermo).end():]
    mode = re.search(r'\bmode\s+inventory\s*;', samples).start()
    depth, begin = 0, mode
    while begin > 0:
        begin -= 1
        if samples[begin] == '}':
            depth += 1
        elif samples[begin] == '{':
            if depth == 0:
                break
            depth -= 1
    depth, end = 0, begin + 1
    while depth >= 0:
        depth += {'{': 1, '}': -1}.get(samples[end], 0)
        end += 1
    body = samples[begin + 1:end - 1]

    def key(name, default=None):
        k = re.search(r'\b' + name + r'\s+([^;\s]+)\s*;', body)
        if not k:
            if default is None:
                raise SystemExit(f'No {name} in the inventory sample')
            return default
        return float(k.group(1))

    n0, aS = key('amount'), key('areaPerVolume')
    scaleS, CeS = key('kineticScale', 1.0), key('Ce', 1.0)
    M = dictValue(os.path.join(case, 'constant', 'speciesTransportProperties'),
                  'molarMass')
    dt = dictValue(os.path.join(case, 'system', 'controlDict'), 'deltaT')
    G0 = scaleS*CeS*2*sigma/(2 - sigma)

    cells = selectCells(case, shape, coords)
    V = readScalarField(os.path.join(case, '0', 'V'))[0]
    T = internalField(os.path.join(case, '0', 'T'), len(V))
    rho = internalField(os.path.join(case, '0', 'rho'), len(V))
    VZ = sum(V[c] for c in cells)
    inventory = n0*M
    coeff = {}
    for c in cells:
        A = aS*V[c]
        G = G0*math.sqrt(M/(2*math.pi*R*T[c]))
        beta = rho[c]*R*T[c]/M
        coeff[c] = (A*G*pEq(T[c], nodes, columns), A*G*beta)

    times = timeDirs(case)
    Cs = {c: inventory/VZ for c in cells}
    nMiss, stepsMiss, worst, total, perStep = 0, 0, 0.0, 0.0, []
    exhaustedSuper, idleUnder = 0, 0
    for t in times:
        Y = internalField(os.path.join(case, t, 'Y_PbI2_g'), len(V))
        exch = internalField(os.path.join(case, t, 'exch_PbI2_g'), len(V))
        missStep = 0
        for c in cells:
            a, g = coeff[c]
            q = a - g*Y[c]
            upper = Cs[c]*V[c]/dt
            flow = exch[c]*V[c]
            qc = min(max(q, 0.0), upper)
            miss = abs(qc - flow)*dt
            total += miss
            worst = max(worst, miss/inventory)
            if miss > limit*inventory:
                missStep += 1
                if upper > 0 and abs(flow - upper) <= 1e-12*upper and q < 0:
                    exhaustedSuper += 1
                if upper > 0 and flow == 0 and q > 0:
                    idleUnder += 1
        nMiss += missStep
        stepsMiss += missStep > 0
        if missStep:
            perStep.append(f'{t}: {missStep}')
        CsNew = internalField(os.path.join(case, t, 'Cs_PbI2_s'), len(V))
        Cs = {c: CsNew[c] for c in cells}
    done(nMiss == 0 and len(times) > 0,
         f'bounded law of the {len(cells)} sample cells (V_Z {VZ:.10g} m3) '
         f'over {len(times)} written steps: {nMiss} elements in {stepsMiss} '
         f'steps miss more than {limit:g} of n0 M (exhausted in '
         f'supersaturated gas {exhaustedSuper}, no flow in undersaturated gas '
         f'{idleUnder}; per step {", ".join(perStep[:5]) or "none"}); '
         f'largest miss {worst:.2e} of n0 M, all misses '
         f'{total/inventory:.2e}')


def meshBody(path):
    """(bytes, start of the data after the FoamFile header, binary?,
    label size) of a polyMesh file."""
    data = open(path, 'rb').read()
    header = data.index(b'}', data.index(b'FoamFile')) + 1
    label = 8 if re.search(rb'label=64', data[:header]) else 4
    return data, header, foamFormat(data) == 'binary', label


def readLabels(path):
    """A labelList of constant/polyMesh (owner)."""
    data, header, binary, label = meshBody(path)
    start = re.search(rb'\n\s*\d+\s*\(', data[header:]).start() + header
    return parseList(data, start, binary, label, 'q' if label == 8 else 'i')[0]


def readPoints(path):
    """The points of constant/polyMesh as (x, y, z) tuples."""
    data, header, binary, _ = meshBody(path)
    m = re.compile(rb'\n\s*(\d+)\s*\(').search(data, header)
    n, start = int(m.group(1)), m.end()
    if binary:
        flat = struct.unpack('<%dd' % (3*n), data[start:start + 24*n])
        return [flat[3*i:3*i + 3] for i in range(n)]
    return [tuple(float(x) for x in t.split()) for t in
            re.findall(rb'\(([^()]*)\)', data[start:])[:n]]


def readFaces(path):
    """The faces of constant/polyMesh (faceCompactList, or faceList in
    ascii) as lists of point labels."""
    data, header, binary, label = meshBody(path)
    code = 'q' if label == 8 else 'i'
    start = re.search(rb'\n\s*\d+\s*\(', data[header:]).start() + header
    if binary:
        offsets, end = parseList(data, start, True, label, code)
        labels = parseList(data, end, True, label, code)[0]
        return [labels[offsets[i]:offsets[i + 1]]
                for i in range(len(offsets) - 1)]
    body = data[data.index(b'(', start) + 1:]
    return [[int(x) for x in f.split()]
            for f in re.findall(rb'\d+\s*\(([^()]*)\)', body)]


def faceArea(points, face):
    """|S_f| of a face (sum of the triangles about its first point)."""
    p0 = points[face[0]]
    S = [0.0, 0.0, 0.0]
    for a, b in zip(face[1:-1], face[2:]):
        u = [points[a][i] - p0[i] for i in range(3)]
        v = [points[b][i] - p0[i] for i in range(3)]
        S[0] += 0.5*(u[1]*v[2] - u[2]*v[1])
        S[1] += 0.5*(u[2]*v[0] - u[0]*v[2])
        S[2] += 0.5*(u[0]*v[1] - u[1]*v[0])
    return math.sqrt(S[0]*S[0] + S[1]*S[1] + S[2]*S[2])


def patchFaces(case, patch):
    """(owner cells, face areas [m2]) of the faces of a patch."""
    poly = os.path.join(case, 'constant', 'polyMesh')
    boundary = stripComments(open(os.path.join(poly, 'boundary')).read())
    m = re.search(r'\b' + re.escape(patch) + r'\s*\{([^}]*)\}', boundary)
    if not m:
        raise SystemExit(f'No patch {patch} in {poly}/boundary')
    start = int(re.search(r'startFace\s+(\d+)', m.group(1)).group(1))
    n = int(re.search(r'nFaces\s+(\d+)', m.group(1)).group(1))
    owner = readLabels(os.path.join(poly, 'owner'))
    faces = readFaces(os.path.join(poly, 'faces'))
    points = readPoints(os.path.join(poly, 'points'))
    return ([owner[f] for f in range(start, start + n)],
            [faceArea(points, faces[f]) for f in range(start, start + n)])


def patchValues(path, patch, n):
    """The values of a patch of a written field, a 'uniform' one expanded;
    zeros when the file does not exist (a fresh start's 0/)."""
    if not os.path.exists(path):
        return [0.0]*n
    values = readScalarField(path)[1][patch]
    return values*n if len(values) == 1 and n > 1 else values


def cmdWallLaw(case, limit, patch='WALL'):
    limit = float(limit)
    thermo = stripComments(open(os.path.join(
        case, 'constant', 'thermochemistryProperties')).read())
    table = re.search(r'file\s+"<constant>/([^"]+)"', thermo).group(1)
    phases = re.findall(r'"([^"]+)"', re.search(r'phases\s*\((.*?)\)\s*;',
                                                thermo).group(1))
    nodes, columns = readTable(os.path.join(case, 'constant', table), phases)
    hks = re.search(r'HKSCoeffs\s*\{([^}]*)\}', thermo)

    def coefficient(name):
        k = re.search(r'\b' + name + r'\s+([^;\s]+);', hks.group(1)) \
            if hks else None
        return float(k.group(1)) if k else 1.0

    sigma = coefficient('accommodation')
    G0 = coefficient('kineticScale')*coefficient('Ce')*2*sigma/(2 - sigma)
    atWall = re.search(r'\binterfaceTemperature\s+wall\s*;', thermo) \
        is not None
    M = dictValue(os.path.join(case, 'constant', 'speciesTransportProperties'),
                  'molarMass')
    dt = dictValue(os.path.join(case, 'system', 'controlDict'), 'deltaT')

    cells, areas = patchFaces(case, patch)
    n = len(cells)
    Tfield = readScalarField(os.path.join(case, '0', 'T'))
    nCells = max(cells) + 1
    Tc = Tfield[0]*nCells if len(Tfield[0]) == 1 else Tfield[0]
    Tpatch = Tfield[1][patch]
    Tpatch = Tpatch*n if len(Tpatch) == 1 else Tpatch
    rho = readScalarField(os.path.join(case, '0', 'rho'))[0]
    rho = rho*nCells if len(rho) == 1 else rho
    coeff = []
    for k, (c, A) in enumerate(zip(cells, areas)):
        Te = Tpatch[k] if atWall else Tc[c]
        G = G0*math.sqrt(M/(2*math.pi*R*Te))
        coeff.append((A*G*pEq(Te, nodes, columns), A*G*rho[c]*R*Tc[c]/M))

    names, rows = readBalance(os.path.join(
        case, 'postProcessing', 'phaseChangeBalance', '0', 'PbI2_g.dat'))
    inventory = abs(column(names, rows[0], 'initial'))
    samples = re.search(r'\bsamples\s*\{', thermo)
    if samples:
        amounts = re.findall(r'\bamount\s+([^;\s]+)\s*;',
                             thermo[samples.end():])
        inventory = max(inventory, sum(float(a) for a in amounts)*M)

    times = timeDirs(case)
    previous = patchValues(os.path.join(case, '0', 'mw_PbI2_s'), patch, n)
    nMiss, stepsMiss, worst, total, perStep = 0, 0, 0.0, 0.0, []
    nImplicit = nNone = nExhausted = 0
    for t in times:
        Y = readScalarField(os.path.join(case, t, 'Y_PbI2_g'))[0]
        Y = Y*nCells if len(Y) == 1 else Y
        mw = patchValues(os.path.join(case, t, 'mw_PbI2_s'), patch, n)
        missStep = 0
        for k, (c, A) in enumerate(zip(cells, areas)):
            a, g = coeff[k]
            upper = previous[k]*A/dt
            flow = (previous[k] - mw[k])*A/dt
            admissible = min(a - g*Y[c], upper)
            miss = abs(admissible - flow)*dt
            total += miss
            worst = max(worst, miss/inventory)
            if miss > limit*inventory:
                missStep += 1
            if a - g*Y[c] < upper:
                nImplicit += 1
            elif upper > 0:
                nExhausted += 1
            else:
                nNone += 1
        nMiss += missStep
        stepsMiss += missStep > 0
        if missStep:
            perStep.append(f'{t}: {missStep}')
        previous = mw
    done(nMiss == 0 and len(times) > 0,
         f'bounded law of the {n} elements of {patch} (interfaceTemperature '
         f'{"wall" if atWall else "cell"}) over {len(times)} written steps: '
         f'{nMiss} element-steps in {stepsMiss} steps miss more than '
         f'{limit:g} of the inventory {inventory:.6g} kg (per step '
         f'{", ".join(perStep[:5]) or "none"}); largest miss {worst:.2e}, '
         f'all misses {total/inventory:.2e}; element-steps deposit/'
         f'evaporate {nImplicit}, exhausted {nExhausted}, bare below Y_eq '
         f'{nNone}')


COMMANDS = {
    'model': cmdModel,
    'approach': cmdApproach,
    'unchanged': cmdUnchanged,
    'evaporate': cmdEvaporate,
    'closureStep': cmdClosureStep,
    'removal': cmdRemoval,
    'decreasing': cmdDecreasing,
    'logCounts': cmdLogCounts,
    'totalLeft': cmdTotalLeft,
    'layoutSample': cmdLayoutSample,
    'complementarity': cmdComplementarity,
    'sameRegimes': cmdSameRegimes,
    'inventories': cmdInventories,
    'wallLaw': cmdWallLaw,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    COMMANDS[sys.argv[1]](*sys.argv[2:])
