#!/usr/bin/env python3
"""
tests/gems/checks.py

PURPOSE

Numerical checks of the M9 tests (tests/gems/Allrun): the GEMS3K backend of
the HKS law (equilibrium GEMS, modes frozen and local).  Every subcommand
prints one summary line and exits with 0 (criterion met) or 1.  The values
are recomputed here from the written fields, independently of the solver:
the p_eq table is read and interpolated as vapourPressureTable.H specifies
(tests/hks/checks.py), the guards of the local mode are evaluated from the
written T, p, rho and Y.

The runs write pEq_PbI2_g [Pa] and gemsSource_PbI2_g (GEMSCoeffs
writeEquilibriumField yes): in the cells of the elements the p_eq of the
law and its source (0 GEMS3K, 1 below minTemperature, 2 below
minMoleFraction, 3 outside the GEMS grid or table, 4 rejected, 5 failed;
-1 elsewhere, or before the first update), on the WALL faces the value of
the element of the face.  On the channel every exchange cell holds one
element (one WALL face, or one cell of the boat).

  frozenElements <case> <time> <Tlo> <Thi> <limit> <minTemperature>
      mode frozen: every element (every cell with a source >= 0, and every
      WALL face) at <time>: |log10 p_eq - log10 p_eq,table(T_e)| <= limit
      for T_e in [Tlo, Thi] (T_e = the cell T, interfaceTemperature cell);
      the source 0 at and above minTemperature, and below it 1 with p_eq
      the table's to 1e-12 (the fallback); prints the largest deviation in
      the range and overall, and the counts
  sources <case> <log> <minTemperature> <minMoleFraction>
      mode local, a run that writes every step: for every written time t_k
      and every element, the source written at t_k against the guards
      evaluated at the composition of the start of that step, Y of t_k-1
      (0 at the first step): 1 if T_e < minTemperature, else 2 if x =
      c/(c + p/(R T)) < minMoleFraction (c = rho Y/M, no gas: c <= 0), else
      0; the counts of every source equal those of the step's GEMS line of
      the log, and no step has a failure or a rejection
  gemsLines <log>
      mode local: every step has a GEMS line, none with a failure or an
      engine re-created, and the summary exists; prints the totals (calls,
      warm, cold, failures, guarded per kind, rejected, the largest
      deviation) and the cost line
  engines <case> <nRanks> <log>
      mode local on nRanks ranks: the log names nRanks engines; every
      processorN/gemsLog holds gemsEquilibrium.log, whose first line is
      'rank N of nRanks', and GEMS3K's ipmlog.txt; the calls of the ranks
      (their update lines) add up to the calls of the log's summary
  profiles <caseA> <caseB> <tol>
      the wall column (the deposit) of every profile file of the native bin
      set at every time both cases wrote: sum |dn'| dx <= tol sum |n'_A| dx
      and max |dn'| <= tol max |n'_A|; prints both and the times
  warmCost <log> <limit [us]> <load text>
      the mean time of a warm call of the summary <= limit; prints it with
      the number of calls and the load text
  fault <bridge log> <solver log> <rank log> <nFatal> <nZeroPgas>
        <nNanOmega> <nSuppressed>
      the fault-injecting bridge (faultBridge.c) in mode local: nFatal
      injected GEMSB_ERR_FATAL, nZeroPgas converged calls whose p_g was
      reported 0 and nNanOmega whose log10 Omega was reported NaN (each
      found in the bridge log; the latter two converged, status 0); after
      each fatal fault the solver destroyed the engine, created a new one,
      suppressed nSuppressed condensed species on it (bounds 0 0), and the
      first call of every element afterwards was cold (state not valid);
      before the first fault every element called again was warm; the
      solver log counts every injected fault as a failure (nFatal +
      nZeroPgas + nNanOmega: a call whose result was not used, plan section
      35, item F) and nFatal re-creations, in the steps and in the summary;
      the rank log names each: ERR_FATAL, 'p_g 0 Pa' and 'max log10 Omega
      nan' of a converged and balanced call
  frozenFault <bridge log> <solver log> <rank log> <nFatal> <nZeroPgas>
        <nNanOmega> <nSuppressed> <nNodes>
      the same in mode frozen (plan section 35, items E and Y): after each
      of the nFatal injected GEMSB_ERR_FATAL the engine is destroyed and
      created again with nSuppressed condensed species suppressed, and
      every later call goes to the new engine; all calls are cold; nZeroPgas
      converged node calls reported with p_g = 0 and nNanOmega with log10
      Omega = NaN fail their node: the rank log names each, a p_g of 0 with
      the node's real max log10 Omega (finite, not the start value -1e+300
      of the build before item Y), a NaN Omega as nan; the start-up line of
      the pair counts nNodes nodes, nNodes - nFailed from GEMS3K, nFailed
      failed (nFatal + nZeroPgas + nNanOmega) and nFatal re-creations, and
      the rank log names the re-creations
  updates <log> <interval>
      updateInterval: the GEMS line of step n (1, 2, ...) says 'update'
      exactly when (n - 1) is a multiple of interval, 'no update this step'
      otherwise
  restartSchedule <continuous log> <restart log> <first step> <interval>
        <restored: yes|no>
      a restart in mode local (plan section 35, item C): the steps of the
      restart, from the time index <first step> on, update exactly when
      (n - 1) is a multiple of interval and, without a restored lagged
      state (restored no), at the first step too; with it, every step line
      before the first update equals the continuous run's line of the same
      step; the first update of the restart is cold (no warm call)
  rankLogs <case> <nRanks> <log>
      every gemsEquilibrium.log below the case is the log of exactly one
      rank (one line 'rank N of nRanks', its first), the nRanks ranks have
      one each in different directories, each with GEMS3K's ipmlog.txt
      (plan section 35, item D), and the main log's 'logs in <directory>'
      names the directory of every rank (processor<N> its number; relative
      to the case, or absolute; plan section 35, item Y)
  iterations <log> <case> <nRanks> <value>
      64-bit counts (plan section 35, item A), the fault bridge reporting
      <value> IPM iterations per call: every update line with calls and
      the summary print <value> iterations per call, and every update line
      of every rank log the total calls x <value> exactly
  speciation <caseNone> <caseLagged> <logLagged>
      speciation lagged as a separate run (reported, plan section 6 M9
      acceptance 4): chi of the summary lies in [1e-12, 1] and below 1
      (the speciation acts); prints chi, the largest relative difference of
      the deposit and the gas inventories from the run without speciation,
      the L1 difference of the deposit profiles, and T_dep of both
"""

import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'phaseChange'))
import checks as pc   # noqa: E402  (tests/phaseChange/checks.py)

_spec = importlib.util.spec_from_file_location(
    'hksChecks', os.path.join(HERE, '..', 'hks', 'checks.py'))
hks = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hks)

R = 8.314462618
NAMES = ['GEMS', 'below minTemperature', 'below minMoleFraction',
         'outside the GEMS grid', 'rejected', 'failed']


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def expand(values, n):
    return values*n if len(values) == 1 and n > 1 else values


def caseTable(case):
    """Nodes and log10 p columns of the pair's table (thermochemistry-
    Properties: the file below <constant> and the phases)."""
    thermo = hks.stripComments(open(os.path.join(
        case, 'constant', 'thermochemistryProperties')).read())
    table = re.search(r'file\s+"<constant>/([^"]+)"', thermo).group(1)
    phases = re.findall(r'"([^"]+)"', re.search(
        r'phases\s*\((.*?)\)\s*;', thermo).group(1))
    return hks.readTable(os.path.join(case, 'constant', table), phases)


def field(case, time, name):
    return pc.readScalarField(os.path.join(case, time, name))


def gemsLines(log):
    """The GEMS lines of the steps, in order: (update, counts dict)."""
    out = []
    pattern = re.compile(
        r'GEMS local: (update \(calls (\d+): warm (\d+), retried (\d+), '
        r'cold (\d+); \S+ iterations per call; failures (\d+), engine '
        r're-created (\d+)\)|no update this step \(updateInterval \d+\)); '
        r'elements[^:]*: GEMS (\d+), table (\d+) \(below minTemperature '
        r'(\d+), below minMoleFraction (\d+), outside the GEMS grid (\d+), '
        r'rejected (\d+), failed (\d+)\); max \|dlog10 p_eq/p_eq,table\| '
        r'([^;\s]+)')
    for line in open(log):
        m = pattern.search(line)
        if not m:
            continue
        update = m.group(2) is not None
        g = lambda i: int(m.group(i)) if m.group(i) is not None else 0
        out.append((update, {
            'calls': g(2), 'warm': g(3), 'retried': g(4), 'cold': g(5),
            'failures': g(6), 'recreated': g(7),
            'sources': [g(8), g(10), g(11), g(12), g(13), g(14)],
            'deviation': float(m.group(15)), 'line': line.strip()}))
    return out


def summary(log):
    text = open(log).read()
    m = re.search(r'GEMS3K summary \(mode local\): (\d+) updates, (\d+) calls '
                  r'\(warm (\d+), retried (\d+), cold (\d+)\), (\d+) failures, '
                  r'engine re-created (\d+) times; (\S+) iterations per call; '
                  r'(?:(\S+) us per warm call|no warm call)', text)
    if not m:
        return None
    return {'updates': int(m.group(1)), 'calls': int(m.group(2)),
            'warm': int(m.group(3)), 'retried': int(m.group(4)),
            'cold': int(m.group(5)), 'failures': int(m.group(6)),
            'recreated': int(m.group(7)), 'iterations': float(m.group(8)),
            'usWarm': float(m.group(9)) if m.group(9) else None,
            'line': m.group(0)}


#------------------------------------------------------------------------------

def cmdFrozenElements(case, time, Tlo, Thi, limit, minT):
    Tlo, Thi, limit, minT = float(Tlo), float(Thi), float(limit), float(minT)
    nodes, columns = caseTable(case)
    T = pc.readScalarField(os.path.join(case, '0', 'T'))[0]
    pEqI, pEqP = field(case, time, 'pEq_PbI2_g')
    srcI, srcP = field(case, time, 'gemsSource_PbI2_g')
    n = len(srcI)
    T = expand(T, n)
    pEqI = expand(pEqI, n)

    worstIn, worstAll, nIn, nGems, nTable, bad = 0.0, 0.0, 0, 0, 0, []
    TinLo, TinHi = 1e300, -1e300
    for c in range(n):
        if srcI[c] < 0:
            continue
        table = hks.pEq(T[c], nodes, columns)
        d = abs(math.log10(pEqI[c]) - math.log10(table))
        if T[c] < minT:
            nTable += 1
            if srcI[c] != 1 or abs(pEqI[c]/table - 1) > 1e-12:
                bad.append(f'cell {c} T {T[c]:.6g}: source {srcI[c]:g}, '
                           f'p_eq/table - 1 = {pEqI[c]/table - 1:.3g}')
            continue
        nGems += 1
        if srcI[c] != 0:
            bad.append(f'cell {c} T {T[c]:.6g}: source {srcI[c]:g}')
        worstAll = max(worstAll, d)
        if Tlo <= T[c] <= Thi:
            nIn += 1
            worstIn = max(worstIn, d)
            TinLo, TinHi = min(TinLo, T[c]), max(TinHi, T[c])

    # the WALL faces: the value of the element of the face (owner cell)
    cells, _ = hks.patchFaces(case, 'WALL')
    wall = expand(pEqP['WALL'], len(cells))
    nFace = sum(1 for k, c in enumerate(cells) if wall[k] != pEqI[c])

    done(not bad and nFace == 0 and nIn > 0 and worstIn <= limit,
         f'mode frozen, {nGems + nTable} elements at {time}: {nIn} with T_e '
         f'in [{Tlo:g}, {Thi:g}] K (from {TinLo:.6g} to {TinHi:.6g} K), '
         f'max |dlog10 p_eq - table| {worstIn:.3e} (limit {limit:g}); over '
         f'all {nGems} GEMS elements {worstAll:.3e}; {nTable} below '
         f'{minT:g} K take the table (source 1, p_eq = table to 1e-12); '
         f'{len(bad)} wrong sources or fallbacks{": " + bad[0] if bad else ""}'
         f'; {nFace} WALL faces differ from their cell')


def cmdSources(case, log, minT, minX):
    minT, minX = float(minT), float(minX)
    times = pc.timeDirs(case)
    lines = gemsLines(log)
    T = pc.readScalarField(os.path.join(case, '0', 'T'))[0]
    p = pc.readScalarField(os.path.join(case, '0', 'p'))[0]
    rho = pc.readScalarField(os.path.join(case, '0', 'rho'))[0]
    M = hks.dictValue(os.path.join(case, 'constant',
                                   'speciesTransportProperties'), 'molarMass')
    previous = '0'
    nMismatch, nSteps, nCountDiff, examples = 0, 0, 0, []
    totals = [0]*6
    for k, t in enumerate(times):
        src = field(case, t, 'gemsSource_PbI2_g')[0]
        n = len(src)
        Tn, pn, rn = expand(T, n), expand(p, n), expand(rho, n)
        Y = expand(field(case, previous, 'Y_PbI2_g')[0], n)
        counts = [0]*6
        for c in range(n):
            if src[c] < 0:
                continue
            cGas = rn[c]*Y[c]/M
            cCarrier = pn[c]/(R*Tn[c])
            if Tn[c] < minT:
                expected = 1
            elif not (cGas > 0 and cGas/(cGas + cCarrier) >= minX):
                expected = 2
            else:
                expected = 0
            counts[int(src[c])] += 1
            if src[c] != expected:
                nMismatch += 1
                if len(examples) < 3:
                    examples.append(f'{t} cell {c}: {src[c]:g} (expected '
                                    f'{expected})')
        logged = lines[k][1]['sources'] if k < len(lines) else None
        if logged != counts:
            nCountDiff += 1
            if len(examples) < 3:
                examples.append(f'{t}: fields {counts}, log {logged}')
        totals = [a + b for a, b in zip(totals, counts)]
        nSteps += 1
        previous = t
    failures = sum(1 for _, c in lines if c['failures'] or c['sources'][4]
                   or c['sources'][5])
    done(nSteps > 0 and nMismatch == 0 and nCountDiff == 0 and failures == 0,
         f'sources of the elements in {nSteps} written steps against the '
         f'guards recomputed from the fields (T < {minT:g} K, x < {minX:g} '
         f'at the composition of the start of the step): {nMismatch} '
         f'elements differ, {nCountDiff} steps whose counts differ from the '
         f'GEMS line of the log; element-steps '
         + ', '.join(f'{NAMES[i]} {totals[i]}' for i in range(6))
         + f'; steps with a failure or rejection {failures}'
         + (f' ({"; ".join(examples)})' if examples else ''))


def cmdGemsLines(log):
    lines = gemsLines(log)
    s = summary(log)
    nSteps = len(re.findall(r'^Time = ', open(log).read(), re.M))
    bad = [c['line'] for _, c in lines if c['failures'] or c['recreated']
           or c['sources'][4] or c['sources'][5]]
    totals = [sum(c['sources'][i] for _, c in lines) for i in range(6)]
    deviation = max((c['deviation'] for _, c in lines), default=0.0)
    done(s is not None and len(lines) == nSteps and nSteps > 0 and not bad,
         f'{len(lines)} GEMS lines in {nSteps} steps, {len(bad)} with a '
         f'failure, a rejection or a re-created engine; element-updates '
         + ', '.join(f'{NAMES[i]} {totals[i]}' for i in range(6))
         + f'; max |dlog10 p_eq/p_eq,table| {deviation:.3e}; '
         + (s['line'] if s else 'no summary'))


def cmdEngines(case, nRanks, log):
    nRanks = int(nRanks)
    text = open(log).read()
    engines = re.search(r'; (\d+) engines? \(one per rank\)', text)
    s = summary(log)
    calls, problems = 0, []
    for r in range(nRanks):
        d = os.path.join(case, f'processor{r}', 'gemsLog')
        lf = os.path.join(d, 'gemsEquilibrium.log')
        if not os.path.isfile(lf) or not os.path.isfile(
                os.path.join(d, 'ipmlog.txt')):
            problems.append(f'no log or ipmlog.txt in {d}')
            continue
        lines = open(lf).read().splitlines()
        if not lines or not lines[0].startswith(f'rank {r} of {nRanks},'):
            problems.append(f'{lf} does not start with rank {r} of {nRanks}')
        for line in lines:
            m = re.match(r'time \S+ \(index \d+\): calls (\d+) ', line)
            if m:
                calls += int(m.group(1))
    extra = os.path.isdir(os.path.join(case, f'processor{nRanks}'))
    ok = (engines is not None and int(engines.group(1)) == nRanks
          and s is not None and calls == s['calls'] and not problems
          and not extra)
    done(ok, f'{engines.group(1) if engines else "no"} engines in the log, '
         f'{nRanks} rank logs with GEMS3K\'s ipmlog.txt, the calls of the '
         f'rank logs {calls} = the summary\'s '
         f'{s["calls"] if s else "(none)"}'
         + (f'; {"; ".join(problems)}' if problems else ''))


def profileRows(case, t):
    rows = []
    for line in open(os.path.join(case, 'postProcessing', 'axialProfiles', t,
                                  'PbI2_g_native.dat')):
        if not line.startswith('#'):
            rows.append([float(x) for x in line.split()])
    return rows


def profileTdep(case, t):
    for line in open(os.path.join(case, 'postProcessing', 'axialProfiles', t,
                                  'PbI2_g_native.dat')):
        m = re.match(r'# Tdep fromInlet status (\w+)(?: T (\S+))?', line)
        if m:
            return m.group(1) + (f' {float(m.group(2)):.4f} K'
                                 if m.group(2) else '')
    return 'none'


def cmdProfiles(caseA, caseB, tol):
    tol = float(tol)
    base = 'postProcessing/axialProfiles'
    times = sorted(set(os.listdir(os.path.join(caseA, base)))
                   & set(os.listdir(os.path.join(caseB, base))), key=float)
    worstL1, worstMax, per = 0.0, 0.0, []
    for t in times:
        a, b = profileRows(caseA, t), profileRows(caseB, t)
        dx = [r[1] - r[0] for r in a]
        l1 = sum(abs(r[4])*w for r, w in zip(a, dx))
        if l1 == 0:
            continue
        dl1 = sum(abs(r[4] - s[4])*w for r, s, w in zip(a, b, dx))/l1
        dmax = max(abs(r[4] - s[4]) for r, s in zip(a, b)) \
            / max(abs(r[4]) for r in a)
        worstL1, worstMax = max(worstL1, dl1), max(worstMax, dmax)
        per.append(f'{t}: {dl1:.2e}/{dmax:.2e}')
    done(len(per) > 0 and worstL1 <= tol and worstMax <= tol,
         f'deposit profiles (native bins) of {len(per)} times: largest L1 '
         f'difference {worstL1:.3e} and largest bin difference {worstMax:.3e} '
         f'of the peak (limit {tol:g}); per time L1/max {", ".join(per)}')


def cmdWarmCost(log, limit, load):
    limit = float(limit)
    s = summary(log)
    if s is None or s['usWarm'] is None:
        done(False, f'no warm call in the summary of {log}')
    done(s['usWarm'] <= limit,
         f'{s["usWarm"]:.4g} us per warm call (limit {limit:g} us), '
         f'{s["warm"]} warm calls of {s["calls"]}, {s["iterations"]:.3g} '
         f'iterations per call; {load}')


def injected(lines, kind):
    """the indices of the lines '<kind> injected' of a bridge log"""
    return [i for i, l in enumerate(lines) if l[:-1] == kind.split()
            and l[-1:] == ['injected']]


def cmdFault(bridgeLog, solverLog, rankLog, nFatal, nZero, nNan, nSuppressed):
    nFaults, nZero, nNan, nSuppressed = int(nFatal), int(nZero), int(nNan), \
        int(nSuppressed)
    lines = [l.split() for l in open(bridgeLog).read().splitlines()]
    faults = injected(lines, 'fatal')
    zeros, nans = injected(lines, 'zero p_g'), injected(lines, 'nan omega')
    problems = []
    for k in zeros + nans:
        # the injected call converged: its 'status 0' line before the note
        if lines[k - 1][:2] != ['status', '0'] and \
                lines[k - 2][:2] != ['status', '0']:
            problems.append(f'the call before line {k + 1} of the bridge log '
                            'did not converge')
    # before the first fatal fault every element called again is warm, but
    # the next call of an element whose result was not used (an injected
    # p_g or Omega), which must be cold: its warm state was invalidated
    seen, invalid, last = set(), set(), None
    rewarmedBefore, notCold, recalled = 0, 0, 0
    first = faults[0] if faults else len(lines)
    for i, l in enumerate(lines[:first]):
        if l[0] == 'equilibrate':
            state, valid = l[5], int(l[7])
            if state in invalid:
                recalled += 1
                notCold += valid != 0
                invalid.discard(state)
            elif state in seen and valid != 1:
                rewarmedBefore += 1
            seen.add(state)
            last = state
        elif i in zeros or i in nans:
            invalid.add(last)
    if rewarmedBefore:
        problems.append(f'{rewarmedBefore} repeated calls before the first '
                        'fault were not warm')
    if notCold or (zeros + nans and recalled == 0):
        problems.append(f'of the elements whose result was not used, '
                        f'{recalled} called again before the first fatal '
                        f'fault, {notCold} of them warm')
    for n, f in enumerate(faults):
        engine = lines[f - 1][3]
        after = lines[f + 1:faults[n + 1] if n + 1 < len(faults) else None]
        steps = [l[0] for l in after[:2 + nSuppressed]]
        if steps[:2] != ['destroy', 'create'] or after[0][1] != engine:
            problems.append(f'fault {n + 1}: not destroy {engine} then create '
                            f'({" ".join(steps[:2])})')
            continue
        new = after[1][1]
        bounds = [l for l in after[2:2 + nSuppressed] if l[0] == 'bounds'
                  and l[1] == new and l[3:] == ['0', '0']]
        if len(bounds) != nSuppressed:
            problems.append(f'fault {n + 1}: {len(bounds)} species suppressed '
                            f'on the new engine, not {nSuppressed}')
        cold, warm, called = 0, 0, set()
        for l in after:
            if l[0] != 'equilibrate':
                continue
            if l[3] != new:
                problems.append(f'fault {n + 1}: a call on engine {l[3]}')
                break
            if l[5] not in called:
                called.add(l[5])
                if int(l[7]) == 0:
                    cold += 1
                else:
                    warm += 1
        if warm:
            problems.append(f'fault {n + 1}: {warm} elements called warm after '
                            'the re-creation')
    gl = gemsLines(solverLog)
    s = summary(solverLog)
    stepsWith = sum(1 for _, c in gl if c['recreated'])
    stepFailures = sum(c['failures'] for _, c in gl)
    rank = open(rankLog).read().splitlines()
    named = {
        'fatal': sum(1 for l in rank if ': ERR_FATAL, ' in l),
        'zero': sum(1 for l in rank if 'converged and balanced' in l
                    and ' p_g 0 Pa' in l),
        'nan': sum(1 for l in rank if 'converged and balanced' in l
                   and 'max log10 Omega nan' in l)}
    expected = nFaults + nZero + nNan
    ok = (len(faults) == nFaults and len(zeros) == nZero
          and len(nans) == nNan and not problems and s is not None
          and s['recreated'] == nFaults and s['failures'] == expected
          and stepFailures == expected and stepsWith == nFaults
          and named == {'fatal': nFaults, 'zero': nZero, 'nan': nNan})
    done(ok, f'injected: {len(faults)} GEMSB_ERR_FATAL (expected {nFaults}), '
         f'{len(zeros)} p_g = 0 and {len(nans)} log10 Omega = NaN of a '
         f'converged call (expected {nZero}, {nNan}); after each fatal one, '
         f'expected: the engine destroyed and created again, {nSuppressed} '
         f'condensed species suppressed again on it, the first call of every '
         f'element cold; found: '
         + ('; '.join(problems[:4]) if problems else 'so')
         + f'; {recalled} elements whose result was not used called again '
         f'before the first fatal fault, {notCold} of them warm'
         + f'; the log: {stepsWith} steps with a re-created engine, '
         f'failures {stepFailures} in the steps and '
         f'{s["failures"] if s else "?"} in the summary (expected '
         f'{expected}), {s["recreated"] if s else "?"} re-creations; the rank '
         f'log names {named["fatal"]} ERR_FATAL, {named["zero"]} p_g 0 and '
         f'{named["nan"]} NaN Omega')


def cmdFrozenFault(bridgeLog, solverLog, rankLog, nFatal, nZero, nNan,
                   nSuppressed, nNodes):
    nFatal, nZero, nNan, nSuppressed, nNodes = int(nFatal), int(nZero), \
        int(nNan), int(nSuppressed), int(nNodes)
    nFailed = nFatal + nZero + nNan
    lines = [l.split() for l in open(bridgeLog).read().splitlines()]
    faults = injected(lines, 'fatal')
    zeros, nans = injected(lines, 'zero p_g'), injected(lines, 'nan omega')
    problems = []
    for k in zeros + nans:
        # the injected call converged: its 'status 0' line before the note
        if lines[k - 1][:2] != ['status', '0'] and \
                lines[k - 2][:2] != ['status', '0']:
            problems.append(f'the call before line {k + 1} of the bridge log '
                            'did not converge')
    for n, f in enumerate(faults):
        engine = lines[f - 1][3]
        after = lines[f + 1:]
        steps = [l[0] for l in after[:2]]
        if steps != ['destroy', 'create'] or after[0][1] != engine:
            problems.append(f'fault {n + 1}: not destroy {engine} then create '
                            f'({" ".join(steps)})')
            continue
        new = after[1][1]
        bounds = [l for l in after[2:2 + nSuppressed] if l[0] == 'bounds'
                  and l[1] == new and l[3:] == ['0', '0']]
        if len(bounds) != nSuppressed:
            problems.append(f'fault {n + 1}: {len(bounds)} species suppressed '
                            f'on the new engine, not {nSuppressed}')
        later = [l for l in after if l[0] == 'equilibrate']
        if any(l[3] != new for l in later[:1]):
            problems.append(f'fault {n + 1}: the next call not on engine {new}')
    calls = [l for l in lines if l[0] == 'equilibrate']
    warm = sum(1 for l in calls if l[5] != '(nil)' or l[7] != '-1')
    if warm:
        problems.append(f'{warm} calls with a warm state')
    text = open(solverLog).read()
    m = re.search(r'GEMS frozen: p_eq at the (\d+) nodes of the table in [^:]+: '
                  r'(\d+) from GEMS3K \([^)]*\), (\d+) failed, (\d+) rejected '
                  r'\(taken from the table\), engine re-created (\d+) times?;',
                  text)
    found = (tuple(int(m.group(k)) for k in range(1, 6)) if m else None)
    rank = open(rankLog).read()
    named = rank.count('engine re-created after GEMSB_ERR_FATAL')
    # the unusable nodes: p_g 0 with the node's real max log10 Omega, and a
    # NaN Omega
    omegas = re.findall(r'converged, but p_g 0 Pa, max log10 Omega (\S+)\)',
                        rank)
    real = [w for w in omegas
            if math.isfinite(float(w)) and abs(float(w)) < 1e3]
    nanNamed = len(re.findall(r'converged, but p_g \S+ Pa, max log10 Omega '
                              r'nan\)', rank))
    ok = (len(faults) == nFatal and len(zeros) == nZero
          and len(nans) == nNan and not problems
          and found == (nNodes, nNodes - nFailed, nFailed, 0, nFatal)
          and named == nFatal and len(omegas) == nZero
          and len(real) == nZero and nanNamed == nNan)
    done(ok, f'mode frozen: {len(faults)} GEMSB_ERR_FATAL, {len(zeros)} p_g '
         f'= 0 and {len(nans)} log10 Omega = NaN injected (expected '
         f'{nFatal}, {nZero}, {nNan}) among {len(calls)} cold calls; after '
         f'each fatal one, expected: the '
         f'engine destroyed and created again with {nSuppressed} condensed '
         f'species suppressed, the next calls on it; found: '
         + ('; '.join(problems[:4]) if problems else 'so')
         + '; the start-up line: '
         + (f'{found[0]} nodes, {found[1]} from GEMS3K, {found[2]} failed, '
            f'{found[3]} rejected, engine re-created {found[4]}' if found
            else 'no frozen line with re-creations')
         + f'; the rank log names {named} re-creations, the p_g 0 nodes with '
         f'max log10 Omega {" ".join(omegas) or "-"} (expected {nZero}, '
         f'finite) and {nanNamed} NaN Omega')


def cmdUpdates(log, interval):
    interval = int(interval)
    lines = gemsLines(log)
    wrong = [n + 1 for n, (u, _) in enumerate(lines)
             if u != (n % interval == 0)]
    done(len(lines) > interval and not wrong,
         f'{len(lines)} steps, updates at steps '
         f'{[n + 1 for n, (u, _) in enumerate(lines) if u][:6]} ..., '
         f'{len(wrong)} steps wrong ({wrong[:5]})')


def cmdRestartSchedule(contLog, restLog, first, interval, restored):
    first, interval, restored = int(first), int(interval), restored == 'yes'
    cont, rest = gemsLines(contLog), gemsLines(restLog)
    wrong, differ, cold = [], [], None
    updates = []
    for k, (u, c) in enumerate(rest):
        n = first + k
        expected = (n - 1) % interval == 0 or (k == 0 and not restored)
        if u != expected:
            wrong.append(n)
        if u:
            updates.append(n)
            if cold is None:
                cold = c['warm'] == 0 and c['calls'] == c['cold'] + c['retried']
        elif restored and not updates:
            if n - 1 >= len(cont) or cont[n - 1][1]['line'] != c['line']:
                differ.append(n)
    before = len([1 for k, (u, _) in enumerate(rest)
                  if not updates or first + k < updates[0]])
    ok = (len(rest) > 0 and not wrong and not differ and cold is not None
          and cold)
    done(ok, f'{len(rest)} steps from time index {first} (updateInterval '
         f'{interval}, lagged state {"restored" if restored else "missing"}):'
         f' updates at {updates[:6]}, {len(wrong)} steps off the schedule '
         f'{wrong[:5]}; '
         + (f'{before} steps before the first update, their lines equal to the'
            f' continuous run\'s but {len(differ)} {differ[:3]}; '
            if restored else '')
         + f'the first update {"cold" if cold else "not cold"}')


def cmdRankLogs(case, nRanks, log):
    nRanks = int(nRanks)
    found, problems = {}, []
    for root, _, files in os.walk(case):
        if 'gemsEquilibrium.log' not in files:
            continue
        lines = open(os.path.join(root, 'gemsEquilibrium.log')).read() \
            .splitlines()
        ranks = [l for l in lines if re.match(r'rank \d+ of \d+,', l)]
        m = re.match(r'rank (\d+) of (\d+),', lines[0]) if lines else None
        if len(ranks) != 1 or not m or int(m.group(2)) != nRanks:
            problems.append(f'{os.path.relpath(root, case)}: {len(ranks)} '
                            f'rank lines')
            continue
        if not os.path.isfile(os.path.join(root, 'ipmlog.txt')):
            problems.append(f'{os.path.relpath(root, case)}: no ipmlog.txt')
        found.setdefault(int(m.group(1)), []).append(
            os.path.relpath(root, case))
    where = re.search(r'one per rank\), logs in (\S+)', open(log).read())
    # the line names the directory of every rank: processor<N> its number,
    # relative to the case (processor<N>/...) or absolute
    if not where:
        problems.append('no line with logs in')
    else:
        for r in sorted(found):
            d = where.group(1).replace('<N>', str(r))
            d = os.path.normpath(d if os.path.isabs(d)
                                 else os.path.join(case, d))
            if os.path.relpath(d, case) != os.path.normpath(found[r][0]):
                problems.append(f'rank {r} writes in {found[r][0]}, the line '
                                f'names {os.path.relpath(d, case)}')
    ok = (sorted(found) == list(range(nRanks))
          and all(len(v) == 1 for v in found.values()) and not problems)
    done(ok, f'{sum(len(v) for v in found.values())} rank logs below the case '
         f'for {nRanks} ranks: '
         + ', '.join(f'{r} {found[r][0]}' for r in sorted(found)[:2])
         + (', ...' if len(found) > 2 else '')
         + f'; the log says: logs in {where.group(1) if where else "?"}'
         + (f'; {"; ".join(problems[:3])}' if problems else ''))


def cmdIterations(log, case, nRanks, value):
    nRanks, value = int(nRanks), int(value)
    lines = [c for u, c in gemsLines(log) if u and c['calls'] > 0]
    text = open(log).read()
    perCall = re.findall(r'update \(calls [1-9]\d*: [^;]*; (\S+) iterations per '
                         r'call', text)
    s = summary(log)
    bad = [x for x in perCall if float(x) != value]
    totals, wrong = 0, []
    for r in range(nRanks):
        lf = os.path.join(case, f'processor{r}', 'gemsLog',
                          'gemsEquilibrium.log')
        for line in open(lf):
            m = re.match(r'time \S+ \(index \d+\): calls (\d+) .*, iterations '
                         r'(\d+)', line)
            if m:
                totals += 1
                if int(m.group(2)) != int(m.group(1))*value:
                    wrong.append(f'rank {r}: {m.group(1)} calls, iterations '
                                 f'{m.group(2)}')
    largest = max((int(m) for m in re.findall(
        r'iterations (\d+)', ''.join(open(os.path.join(
            case, f'processor{r}', 'gemsLog', 'gemsEquilibrium.log')).read()
            for r in range(nRanks)))), default=0)
    ok = (len(lines) > 0 and len(perCall) == len(lines) and not bad
          and s is not None and s['iterations'] == value and totals > 0
          and not wrong)
    done(ok, f'{value} IPM iterations reported per call: {len(perCall)} update '
         f'lines with calls print {set(perCall) if perCall else "none"} '
         f'iterations per call, the summary {s["iterations"] if s else "?"}; '
         f'{totals} update lines of {nRanks} rank logs, {len(wrong)} totals '
         f'other than calls x {value} (the largest {largest}, '
         f'{largest/2**31:.1f} times 2^31)'
         + (f': {"; ".join(wrong[:2])}' if wrong else ''))


def cmdSpeciation(caseNone, caseLagged, log):
    text = open(log).read()
    m = re.search(r'GEMS3K summary PbI2_g: .*; chi (\S+)-(\S+) \(clipped '
                  r'(\d+)\)', text)
    if not m:
        done(False, f'no chi in the summary of {log}')
    lo, hi = float(m.group(1)), float(m.group(2))
    namesA, rowsA = pc.readBalance(os.path.join(
        caseNone, 'postProcessing', 'phaseChangeBalance', '0', 'PbI2_g.dat'))
    namesB, rowsB = pc.readBalance(os.path.join(
        caseLagged, 'postProcessing', 'phaseChangeBalance', '0',
        'PbI2_g.dat'))
    a, b = rowsA[-1], rowsB[-1]
    ref = pc.reference(namesA, a)
    dWall = (pc.column(namesB, b, 'wall') - pc.column(namesA, a, 'wall'))/ref
    dGas = (pc.column(namesB, b, 'gas') - pc.column(namesA, a, 'gas'))/ref
    t = sorted(os.listdir(os.path.join(caseNone, 'postProcessing',
                                       'axialProfiles')), key=float)[-1]
    ra, rb = profileRows(caseNone, t), profileRows(caseLagged, t)
    dx = [r[1] - r[0] for r in ra]
    l1 = sum(abs(r[4])*w for r, w in zip(ra, dx))
    dl1 = sum(abs(r[4] - s[4])*w for r, s, w in zip(ra, rb, dx))/max(l1, 1e-300)
    done(1e-12 <= lo <= hi <= 1 and lo < 1,
         f'speciation lagged: chi {lo:.4g}-{hi:.6g} over the GEMS elements '
         f'(clipped {m.group(3)}); at {a[0]} s, against speciation none: '
         f'wall {dWall:+.3e} and gas {dGas:+.3e} of the inventory, deposit '
         f'profile L1 {dl1:.3e}; T_dep (fromInlet, wall) {profileTdep(caseNone, t)}'
         f' without, {profileTdep(caseLagged, t)} with speciation')


COMMANDS = {
    'frozenElements': cmdFrozenElements,
    'sources': cmdSources,
    'gemsLines': cmdGemsLines,
    'engines': cmdEngines,
    'profiles': cmdProfiles,
    'warmCost': cmdWarmCost,
    'fault': cmdFault,
    'frozenFault': cmdFrozenFault,
    'updates': cmdUpdates,
    'restartSchedule': cmdRestartSchedule,
    'rankLogs': cmdRankLogs,
    'iterations': cmdIterations,
    'speciation': cmdSpeciation,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print(__doc__)
        sys.exit(2)
    COMMANDS[sys.argv[1]](*sys.argv[2:])
