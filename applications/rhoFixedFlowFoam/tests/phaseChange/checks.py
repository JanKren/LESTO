#!/usr/bin/env python3
"""
tests/phaseChange/checks.py

PURPOSE

Numerical checks of the phase-change tests (tests/phaseChange/Allrun).  Every
subcommand prints one summary line and exits with 0 (criterion met) or 1.

  closure   <dat> <tol>
      max over all steps of |closure|/reference of a balance file
      postProcessing/phaseChangeBalance/<t0>/<gas>.dat, with the solver's
      reference max(|initial + released|, |held|, inflow) (phaseChangeLedger.H;
      it is initial + released unless the pair is fed through a patch)
  released  <dat> <n0 [mol]> <M [kg/mol]> <tol>
      |released(last step)/(n0 M) - 1|
  releasedPerStep <dat> <n0> <M> <t start> <dt> <window start> <duration> <tol>
      the release of every step (difference of the cumulative column)
      equals n0 M |[t_k, t_k + dt] ^ window|/length to tol, relative to
      the release of a full step, with the exact step times t_k of the
      solver's clock (t_0 = t start, t_k+1 = t_k + dt in double precision)
      and the window length as the solver represents it, length =
      (start + duration) - start (phaseChangeKinetics.H,
      releaseWindowLength), and the time column holds exactly these t_k+1
  defect    <dat> <tol>
      |solverDefect|/reference at the last step (reference as for closure)
  nonzero   <dat> <column>
      the column is non-zero at the last step (e.g. solverDefect, restart)
  compare   <dat A> <dat B> <tol> <column> [<column> ...]
      relative difference of the given columns at the last common time
  compareSteps <dat A> <dat B> <tol> <column> [<column> ...]
      the same at every step: both files hold the same times, and the
      maximum over them of the relative difference of each column (with
      the time of the maximum) is at most tol
  closedForm <case> <Y0> <lambdaDt> <tol>
      every cell of Y_PbI2_g at the n-th written time equals
      Y0/(1 + lambdaDt)^n (ascii fields, 17 digits)
  keepsDeposit <case> <element>
      review item B2: the deposit of the element (face of WALL, mw_PbI2_s)
      is positive after the first step and never decreases, and at the last
      time Y_PbI2_g of its cell is negative (a negative predictor)
  residualNotUsed <log> <dat>
      parallel: gSum(fvMatrix::residual()) printed in the log differs from
      the ledger's own per-cell defect by more than 1e3 times the largest
      |closure| (i.e. the ledger would not close with it), while the
      closure itself stays below 1e-13
  compareFields <case A> <case B> <time> <field> <tol>
      max relative difference of an ascii internal field between two cases
  closedFormMisses <case> <Y0> <lambdaDt> <min>
      the opposite of closedForm: the max relative error is at least min
      (a control run that must NOT follow the closed form)
  depositsNonDecreasing <case>
      every element deposit (mw_PbI2_s on WALL, ascii) never decreases from
      one written time to the next, and at least one is positive
  derivedCondensate <case> <M [kg/mol]> <h [m]>
      closed box (one wall face per cell of height h): at every written
      time Y_PbI2_s = mw_PbI2_s/(rho h) to 1e-14 (the derived condensate of
      the element deposits of the same step) and c_PbI2_s = rho Y_PbI2_s/M
      to 1e-15 (c_ written from the same state: plan section 12, NB5)
  sameRows <dat A> <dat B>
      every data line of B is byte-identical to the line of A with the same
      time (e.g. a restart continues the ledger of the continuous run)
  logClosure <log> <tol>
      max |relative closure| over the Balance lines of a log (the relative
      value the solver prints)
  sameAsColumn <dat> <column> <function object .dat> <tol> [<factor>]
      the last column of a function-object file (one row per step, e.g.
      surfaceFieldValue areaIntegrate of mw_PbI2_s on WALL), times factor
      (default 1), equals the balance-file column in every step, to tol
      relative to the largest |value|, and it changes between write times
      (i.e. it was not sampled from a field refreshed at write times only)
  stripLayout <file> <key> [<key> ...]
      removes the top-level entries <key> (sub-dictionaries) from a
      uniform/phaseChangeLayout file, as written by an earlier build
  processorPatches <case> <time> <field>
      decomposed case (uncollated, ascii or binary files): on every
      processor patch procBoundaryAtoB of every rank A the written values
      of the scalar field are bit for bit the internal values of rank B in
      the cells of its patch procBoundaryBtoA (the faces of the two patches
      are in the same order), i.e. the processor patches were evaluated
      before the field was written; at least one such value is non-zero

Balance-file columns are addressed by the names of its '#' header line.
A balance file with a value that is not finite (nan, inf) is refused by
every subcommand that reads one (readBalance: SystemExit naming the column,
the time and the file, exit status 1, no verdict line), and logClosure
refuses a printed closure that is not finite: their largest values are a
max(), which passes over a NaN that is not the first value (max(0, nan) is
0), so such a ledger passed (plan section 35, item X).
"""

import math
import os
import re
import struct
import sys


def readBalance(path):
    """Header names and rows (lists of floats; time kept as a string)."""
    names = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('# time'):
                names = line[1:].split()
            elif line.startswith('#') or not line.strip():
                continue
            else:
                parts = line.split()
                rows.append([parts[0]] + [float(x) for x in parts[1:]])
    if names is None or not rows:
        raise SystemExit(f'No header or no data in {path}')
    for r in rows:
        for name, value in zip(names[1:], r[1:]):
            if not math.isfinite(value):
                raise SystemExit(f'Not finite: {name} {value} at t = {r[0]} '
                                 f'in {path}')
    return names, rows


def column(names, row, name):
    return row[names.index(name)]


def reference(names, row):
    """max(|initial + released|, |held|, inflow): phaseChangeLedger.H"""
    supplied = abs(column(names, row, 'initial') + column(names, row, 'released'))
    held = abs(column(names, row, 'gas') + column(names, row, 'wall')
               + column(names, row, 'sample'))
    inflow = sum(max(-row[i], 0.0) for i, n in enumerate(names)
                 if n.startswith('transport_'))
    return max(supplied, held, inflow)


def readInternal(path, nCells=None):
    """Internal field of an ascii OpenFOAM scalar field."""
    text = open(path).read()
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)',
                  text)
    if m:
        return [float(x) for x in m.group(2).split()]
    m = re.search(r'internalField\s+uniform\s+(\S+);', text)
    if m and nCells:
        return [float(m.group(1))]*nCells
    raise SystemExit(f'Cannot read the internal field of {path}')


def readPatch(path, patch):
    """Values of an ascii patch field (calculated, 'value' entry)."""
    text = open(path).read()
    m = re.search(patch + r'\s*\{[^}]*?value\s+(nonuniform\s+List<scalar>\s*\d+\s*'
                  r'\(([^)]*)\)|uniform\s+(\S+));', text)
    if not m:
        raise SystemExit(f'No patch {patch} with a value in {path}')
    if m.group(2) is not None:
        return [float(x) for x in m.group(2).split()]
    return [float(m.group(3))]


def foamFormat(data):
    """'ascii' or 'binary' from the FoamFile header of a file's bytes."""
    m = re.search(rb'format\s+(\w+)\s*;', data)
    return m.group(1).decode() if m else 'ascii'


def parseList(data, pos, binary, itemSize, code):
    """A list 'N (...)' at pos (text or raw block); returns (values, end)."""
    m = re.compile(rb'\s*(\d+)\s*(\(\s*\))?').match(data, pos)
    if m and int(m.group(1)) == 0:
        return [], m.end()            # '0;' (binary) or '0()' (ascii)
    m = re.compile(rb'\s*(\d+)\s*\(').match(data, pos)
    if not m:
        raise SystemExit('No list at byte %d' % pos)
    n, start = int(m.group(1)), m.end()
    if binary:
        end = start + n*itemSize
        values = list(struct.unpack('<%d%s' % (n, code), data[start:end]))
        if data[end:end + 1] != b')':
            raise SystemExit('Binary list of %d items not closed' % n)
        return values, end + 1
    end = data.index(b')', start)
    cast = float if code == 'd' else int
    return [cast(x) for x in data[start:end].split()], end + 1


def parseEntry(data, pos, binary, nValues=None):
    """'uniform <v>;' or 'nonuniform List<scalar> N (...)' at pos."""
    m = re.compile(rb'\s*uniform\s+([^;\s]+)\s*;').match(data, pos)
    if m:
        return [float(m.group(1))]*(nValues if nValues is not None else 1)
    m = re.compile(rb'\s*nonuniform\s+List<scalar>').match(data, pos)
    if not m:
        raise SystemExit('No field entry at byte %d' % pos)
    return parseList(data, m.end(), binary, 8, 'd')[0]


def readScalarField(path):
    """(internal values, {patch: values}) of an ascii or binary scalar field
    (patches with a 'value' entry only)."""
    data = open(path, 'rb').read()
    binary = foamFormat(data) == 'binary'
    m = re.search(rb'\binternalField\b', data)
    internal = parseEntry(data, m.end(), binary)
    boundary = data.index(b'boundaryField', m.end())
    patches = {}
    for pm in re.finditer(rb'\n\s*([A-Za-z]\w*)\s*\n\s*\{\s*\n\s*type\s+\w+\s*;',
                          data[boundary:]):
        name = pm.group(1).decode()
        start = boundary + pm.end()
        close = data.index(b'}', start) if not binary else None
        vm = re.compile(rb'(?:[^{}v]|v(?!alue\s))*value\s+', re.S).match(data, start)
        if vm and (close is None or vm.end() < close):
            patches[name] = parseEntry(data, vm.end(), binary)
    return internal, patches


def readPolyMesh(procDir):
    """Processor patches {name: (startFace, nFaces, neighbProcNo)} and the
    owner list of a processor mesh."""
    text = open(os.path.join(procDir, 'constant', 'polyMesh', 'boundary')).read()
    patches = {}
    for m in re.finditer(r'(\w+)\s*\{([^}]*)\}', text):
        body = m.group(2)
        if re.search(r'type\s+processor\s*;', body):
            get = lambda k: int(re.search(k + r'\s+(\d+)\s*;', body).group(1))
            patches[m.group(1)] = (get('startFace'), get('nFaces'),
                                   get('neighbProcNo'))
    data = open(os.path.join(procDir, 'constant', 'polyMesh', 'owner'), 'rb').read()
    binary = foamFormat(data) == 'binary'
    header = data.index(b'}', data.index(b'FoamFile')) + 1
    label = 8 if re.search(rb'label=64', data[:header]) else 4
    owner = parseList(data, re.search(rb'\n\s*\d+\s*\(', data[header:]).start()
                      + header, binary, label, 'q' if label == 8 else 'i')[0]
    return patches, owner


def timeDirs(case):
    out = []
    for d in os.listdir(case):
        try:
            if float(d) > 0 and os.path.isdir(os.path.join(case, d)):
                out.append(d)
        except ValueError:
            pass
    return sorted(out, key=float)


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def cmdClosure(dat, tol):
    names, rows = readBalance(dat)
    worst = max(abs(column(names, r, 'closure'))/max(reference(names, r), 1e-300)
                for r in rows)
    done(worst <= float(tol),
         f'max |closure|/reference = {worst:.3e} over {len(rows)} '
         f'steps (limit {tol})')


def cmdReleased(dat, n0, M, tol):
    names, rows = readBalance(dat)
    released = column(names, rows[-1], 'released')
    err = abs(released/(float(n0)*float(M)) - 1)
    done(err <= float(tol),
         f'released {released!r} kg vs n0 M = {float(n0)*float(M)!r} kg: '
         f'relative {err:.3e} (limit {tol})')


def cmdReleasedPerStep(dat, n0, M, tStart, dt, start, duration, tol):
    names, rows = readBalance(dat)
    n0M, dt, start, duration = float(n0)*float(M), float(dt), float(start), \
        float(duration)
    length = (start + duration) - start
    full = n0M*dt/length
    t0, previous = float(tStart), 0.0
    worst, clockOk, nPartial = 0.0, True, 0
    for r in rows:
        t1 = t0 + dt
        clockOk = clockOk and float(r[0]) == t1
        overlap = max(min(t1, start + duration) - max(t0, start), 0.0)
        if 0 < overlap < dt:
            nPartial += 1
        released = column(names, r, 'released')
        worst = max(worst, abs((released - previous) - n0M*overlap/length)/full)
        previous, t0 = released, t1
    done(worst <= float(tol) and clockOk,
         f'{len(rows)} steps ({nPartial} partly in the window): max |released '
         f'per step - n0 M overlap/t_ev|/(n0 M dt/t_ev) = {worst:.3e} (limit {tol});'
         f' time column = exact step times: {clockOk}')


def cmdDefect(dat, tol):
    names, rows = readBalance(dat)
    r = rows[-1]
    value = abs(column(names, r, 'solverDefect'))/max(reference(names, r), 1e-300)
    done(value <= float(tol),
         f'|solverDefect|/reference = {value:.3e} at t = {r[0]} '
         f'(limit {tol})')


def cmdNonzero(dat, name):
    names, rows = readBalance(dat)
    value = column(names, rows[-1], name)
    done(value != 0, f'{name} = {value!r} kg at t = {rows[-1][0]}')


def cmdCompare(datA, datB, tol, *cols):
    namesA, rowsA = readBalance(datA)
    namesB, rowsB = readBalance(datB)
    common = sorted(set(r[0] for r in rowsA) & set(r[0] for r in rowsB), key=float)
    if not common:
        done(False, 'no common time')
    t = common[-1]
    rA = next(r for r in rowsA if r[0] == t)
    rB = next(r for r in rowsB if r[0] == t)
    worst = 0
    parts = []
    for name in cols:
        a, b = column(namesA, rA, name), column(namesB, rB, name)
        d = abs(a - b)/max(abs(a), abs(b), 1e-300)
        worst = max(worst, d)
        parts.append(f'{name} {d:.2e}')
    done(worst <= float(tol),
         f'at t = {t}: relative differences ' + ', '.join(parts) + f' (limit {tol})')


def cmdCompareSteps(datA, datB, tol, *cols):
    namesA, rowsA = readBalance(datA)
    namesB, rowsB = readBalance(datB)
    rowsOfB = {r[0]: r for r in rowsB}
    sameTimes = [r[0] for r in rowsA] == [r[0] for r in rowsB]
    worst = {name: (0.0, None) for name in cols}
    for rA in rowsA:
        rB = rowsOfB.get(rA[0])
        if rB is None:
            continue
        for name in cols:
            a, b = column(namesA, rA, name), column(namesB, rB, name)
            d = abs(a - b)/max(abs(a), abs(b), 1e-300)
            if worst[name][1] is None or d > worst[name][0]:
                worst[name] = (d, rA[0])
    overall = max(w[0] for w in worst.values())
    done(sameTimes and overall <= float(tol),
         f'{len(rowsA)} steps (same times: {sameTimes}): max relative '
         f'differences ' + ', '.join(f'{n} {worst[n][0]:.2e} at t = {worst[n][1]}'
                                     for n in cols) + f' (limit {tol})')


def cmdClosedForm(case, Y0, lambdaDt, tol):
    Y0, lambdaDt = float(Y0), float(lambdaDt)
    worst = 0
    times = timeDirs(case)
    for n, t in enumerate(times, 1):
        exact = Y0/(1 + lambdaDt)**n
        Y = readInternal(os.path.join(case, t, 'Y_PbI2_g'), 10)
        worst = max(worst, max(abs(y/exact - 1) for y in Y))
    done(len(times) > 0 and worst <= float(tol),
         f'Y^n = Y0/(1 + lambda dt)^n, lambda dt = {lambdaDt:.6f}: max relative '
         f'error {worst:.3e} over {len(times)} steps and 10 cells (limit {tol})')


def cmdKeepsDeposit(case, element):
    e = int(element)
    times = timeDirs(case)
    deposits = [readPatch(os.path.join(case, t, 'mw_PbI2_s'), 'WALL')[e]
                for t in times]
    lastY = readInternal(os.path.join(case, times[-1], 'Y_PbI2_g'), 10)[e]
    increasing = all(b >= a for a, b in zip(deposits, deposits[1:]))
    ok = deposits[0] > 0 and increasing and lastY < 0
    done(ok, f'element {e}: deposit {deposits[0]!r} kg/m2 after the first step, '
             f'{deposits[-1]!r} at t = {times[-1]} (never decreasing: '
             f'{increasing}); Y of its cell at t = {times[-1]}: {lastY!r}')


def cmdResidualNotUsed(log, dat):
    own, of = [], []
    for line in open(log):
        m = re.search(r'solver defect (\S+) kg/s; gSum\(fvMatrix::residual\(\)\) '
                      r'(\S+) kg/s', line)
        if m:
            own.append(float(m.group(1)))
            of.append(float(m.group(2)))
    names, rows = readBalance(dat)
    closures = [abs(column(names, r, 'closure')) for r in rows]
    refs = [max(reference(names, r), 1e-300) for r in rows]
    worstClosure = max(c/r for c, r in zip(closures, refs))
    # closure the ledger would show with gSum(residual()) instead of its own
    # defect: accumulated (own - OF)*dt, dt from the time column
    times = [float(r[0]) for r in rows]
    dts = [times[0]] + [b - a for a, b in zip(times, times[1:])]
    shift, worstAlt = 0, 0
    for i in range(min(len(own), len(rows))):
        shift += (own[i] - of[i])*dts[i]
        worstAlt = max(worstAlt, abs(closures[i] + shift)/refs[i])
    ok = len(own) > 0 and worstClosure <= 1e-13 and worstAlt > 1e3*worstClosure
    done(ok, f'{len(own)} steps: own defect up to {max(map(abs, own)):.3e} kg/s, '
             f'gSum(fvMatrix::residual()) up to {max(map(abs, of)):.3e} kg/s; '
             f'closure {worstClosure:.3e}, with gSum(residual()) it would be '
             f'{worstAlt:.3e}')


def cmdCompareFields(caseA, caseB, time, field, tol):
    a = readInternal(os.path.join(caseA, time, field))
    b = readInternal(os.path.join(caseB, time, field))
    scale = max(max(map(abs, a)), 1e-300)
    worst = max(abs(x - y) for x, y in zip(a, b))/scale
    done(len(a) == len(b) and worst <= float(tol),
         f'{field} at {time}: max |A - B|/max|A| = {worst:.3e} (limit {tol})')


def cmdClosedFormMisses(case, Y0, lambdaDt, minimum):
    Y0, lambdaDt = float(Y0), float(lambdaDt)
    worst = 0
    times = timeDirs(case)
    for n, t in enumerate(times, 1):
        exact = Y0/(1 + lambdaDt)**n
        Y = readInternal(os.path.join(case, t, 'Y_PbI2_g'), 10)
        worst = max(worst, max(abs(y/exact - 1) for y in Y))
    done(len(times) > 0 and worst >= float(minimum),
         f'control: max relative deviation from Y0/(1 + {lambdaDt:.6f})^n is '
         f'{worst:.3e} (must be at least {minimum})')


def cmdDepositsNonDecreasing(case):
    times = timeDirs(case)
    deposits = [readPatch(os.path.join(case, t, 'mw_PbI2_s'), 'WALL')
                for t in times]
    nElements = len(deposits[0]) if deposits else 0
    decreasing = sum(1 for a, b in zip(deposits, deposits[1:])
                     for x, y in zip(a, b) if y < x)
    positive = sum(1 for x in deposits[-1] if x > 0) if deposits else 0
    done(len(times) > 1 and decreasing == 0 and positive > 0,
         f'{nElements} element deposits over {len(times)} times: '
         f'{decreasing} decreases, {positive} positive at t = {times[-1]}')


def cmdDerivedCondensate(case, M, h):
    M, h = float(M), float(h)
    rho = readInternal(os.path.join(case, '0', 'rho'), 10)
    worstY, worstC = 0.0, 0.0
    times = timeDirs(case)
    for t in times:
        Ys = readInternal(os.path.join(case, t, 'Y_PbI2_s'), 10)
        cs = readInternal(os.path.join(case, t, 'c_PbI2_s'), 10)
        mw = readPatch(os.path.join(case, t, 'mw_PbI2_s'), 'WALL')
        for i in range(10):
            expectedY = mw[i]/(rho[i]*h)
            worstY = max(worstY, abs(Ys[i] - expectedY)/max(abs(expectedY), 1e-300))
            expectedC = rho[i]*Ys[i]/M
            worstC = max(worstC, abs(cs[i] - expectedC)/max(abs(expectedC), 1e-300))
    done(len(times) > 0 and worstY <= 1e-14 and worstC <= 1e-15,
         f'{len(times)} write times: Y_PbI2_s vs mw_PbI2_s/(rho h) {worstY:.2e} '
         f'(limit 1e-14), c_PbI2_s vs rho Y_PbI2_s/M {worstC:.2e} (limit 1e-15)')


def dataLines(path):
    lines = {}
    for line in open(path):
        if line.strip() and not line.startswith('#'):
            lines[line.split()[0]] = line
    return lines


def cmdSameRows(datA, datB):
    a, b = dataLines(datA), dataLines(datB)
    missing = [t for t in b if t not in a]
    different = [t for t in b if t in a and a[t] != b[t]]
    done(len(b) > 0 and not missing and not different,
         f'{len(b)} rows of the restart: {len(b) - len(missing) - len(different)}'
         f' byte-identical to the continuous run, {len(different)} different, '
         f'{len(missing)} without a row at the same time')


def cmdLogClosure(log, tol):
    values = [float(m.group(1)) for m in
              (re.search(r'^Balance .* closure \S+ \(relative (\S+)\)', line)
               for line in open(log)) if m]
    finite = all(math.isfinite(v) for v in values)
    worst = (max(map(abs, values)) if values and finite
             else float('nan') if values else float('inf'))
    done(len(values) > 0 and finite and worst <= float(tol),
         f'printed relative closure up to {worst:.3e} in {len(values)} '
         f'Balance lines (limit {tol})'
         + ('' if finite else '; a closure that is not finite'))


def cmdSameAsColumn(dat, name, foDat, tol, factor='1'):
    names, rows = readBalance(dat)
    ledger = [column(names, r, name) for r in rows]
    fo = [float(factor)*float(line.split()[-1]) for line in open(foDat)
          if line.strip() and not line.startswith('#')]
    scale = max(max(map(abs, ledger)), 1e-300)
    n = min(len(ledger), len(fo))
    worst = max(abs(a - b) for a, b in zip(ledger[:n], fo[:n]))/scale
    nChanged = sum(1 for a, b in zip(fo, fo[1:]) if a != b)
    scaled = '' if float(factor) == 1 else f' ({factor} x the function object)'
    done(len(fo) == len(ledger) and worst <= float(tol)
         and nChanged >= len(fo)//2,
         f'{len(fo)} function-object rows{scaled} vs {len(ledger)} ledger '
         f'rows: max |difference|/max|{name}| = {worst:.3e} (limit {tol}); '
         f'the value changes in {nChanged} of {len(fo) - 1} steps')


def cmdStripLayout(path, *keys):
    text = open(path, 'rb').read().decode('latin-1')
    for key in keys:
        m = re.search(r'\n' + key + r'\s*\{', text)
        if not m:
            done(False, f'no entry {key} in {path}')
        depth, i = 1, m.end()
        while depth:
            depth += {'{': 1, '}': -1}.get(text[i], 0)
            i += 1
        text = text[:m.start() + 1] + text[i:].lstrip('\n')
    open(path, 'wb').write(text.encode('latin-1'))
    done(True, f'removed {" ".join(keys)} from {path}')


def cmdProcessorPatches(case, time, field):
    procs = sorted((d for d in os.listdir(case) if re.fullmatch(r'processor\d+', d)),
                   key=lambda d: int(d[9:]))
    meshes = {int(d[9:]): readPolyMesh(os.path.join(case, d)) for d in procs}
    fields = {int(d[9:]): readScalarField(os.path.join(case, d, time, field))
              for d in procs}
    nPatches, nFaces, nNonZero, nWrong = 0, 0, 0, 0
    for a in meshes:
        for name, (start, n, b) in meshes[a][0].items():
            if n == 0:
                continue
            values = fields[a][1][name]
            if len(values) == 1 and n > 1:
                values = values*n
            other = 'procBoundary%dto%d' % (b, a)
            startB, nB, _ = meshes[b][0][other]
            ownerB = meshes[b][1]
            internalB = fields[b][0]
            if len(internalB) == 1:
                internalB = internalB*(max(ownerB) + 1)
            expected = [internalB[ownerB[startB + i]] for i in range(nB)]
            nPatches += 1
            nFaces += n
            nNonZero += sum(1 for v in values if v != 0)
            nWrong += sum(1 for v, w in zip(values, expected) if v != w) \
                + abs(len(values) - len(expected))
    done(nPatches > 0 and nWrong == 0 and nNonZero > 0,
         f'{field} at {time}: {nPatches} processor patches, {nFaces} faces, '
         f'{nWrong} values differ from the neighbour cells, {nNonZero} non-zero')


commands = {
    'closure': cmdClosure,
    'released': cmdReleased,
    'releasedPerStep': cmdReleasedPerStep,
    'defect': cmdDefect,
    'nonzero': cmdNonzero,
    'compare': cmdCompare,
    'compareSteps': cmdCompareSteps,
    'closedForm': cmdClosedForm,
    'keepsDeposit': cmdKeepsDeposit,
    'residualNotUsed': cmdResidualNotUsed,
    'compareFields': cmdCompareFields,
    'closedFormMisses': cmdClosedFormMisses,
    'depositsNonDecreasing': cmdDepositsNonDecreasing,
    'derivedCondensate': cmdDerivedCondensate,
    'sameRows': cmdSameRows,
    'logClosure': cmdLogClosure,
    'sameAsColumn': cmdSameAsColumn,
    'stripLayout': cmdStripLayout,
    'processorPatches': cmdProcessorPatches,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in commands:
        print(__doc__)
        sys.exit(2)
    commands[sys.argv[1]](*sys.argv[2:])
