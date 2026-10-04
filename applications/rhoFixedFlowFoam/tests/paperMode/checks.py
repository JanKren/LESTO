#!/usr/bin/env python3
"""
tests/paperMode/checks.py

PURPOSE

Numerical checks of milestone M4 (tests/paperMode/Allrun): the paper-mode
carrier of utilities/setFrozenCarrier, the interaction layer 'distance' and
the code-to-code case run/014 of the temperature-based model of Lobresco et
al. (2026).  Every subcommand prints one summary line and exits with 0
(criterion met) or 1.

  continuity <log> <tol> [<min lines>]
      every 'Continuity (...)' line of a setFrozenCarrier log reports a
      max |sum_f phi| per cell of at most tol times the inflow (and at least
      <min lines> such lines, default 1)
  continuityAbove <log> <min>
      the opposite (negative evidence): a 'Continuity' line reports more
      than min times the inflow
  carrier <log> <key=value> ...
      the carrier line of the log: U_b, the discrete inflow (relative to Q
      sector), the facet-area ratio, the discrete bulk velocity; each
      key=value is checked to 1e-9 relative (keys Ub, inflowRatio,
      facetRatio, bulk)
  layerCounts <log> <rule> <cells> <A/V> <tol>
      the start-up line 'Interface (layer <rule>)': the number of elements
      (cells) and A_I/V_I to tol relative
  mdep <balance .dat> <function object .dat> <tol>
      areaIntegrate of mDep_PbI2_s on the interface patch (surfaceFieldValue,
      every step) equals the wall inventory of the balance file in every
      step to tol relative to the largest wall inventory, and the wall is
      non-zero at the end
  facesAgree <field A> <field B> <patch> <tol>
      the values of a scalar field on a patch agree between two cases (e.g.
      mDep_ of a serial and a reconstructed parallel run) to tol relative
      to the largest value, and a symmetric channel deposits the same on
      both walls (the first and the second half of the patch faces in the
      blockMesh order, to tol)
  plateau <case> <time> <t_ev> <n0> <x lo> <x hi> <tol> <U_b,d> <U_b> <paper>
      Fig. 9 identity: on the native gas profile at <time> every bin whose
      centre lies in [x lo, x hi] has n' t_ev = n0/U_b,d to tol (U_b,d the
      discrete bulk velocity of the carrier: its inflow over its inlet
      area); prints the comparison with n0/U_b (the analytic U_b) and with
      the paper's value
  solidAgree <case A> <time A> <case B> <time B> <tol>
      the solid (wall) profiles of the native set agree: sum_b |n'_A - n'_B|
      dx_b <= tol sum_b |n'_B| dx_b (the L1 norm relative to the deposit of
      B); also prints the largest difference relative to the peak of B, the
      deposits, the gas left and T_dep of both
  tdep <case> <time> <t end> <expected> <tol> <label>
      T_dep fromInlet of the wall column of the native set (fraction of the
      file) at <time>, and, when <time> is before the end of the experiment
      <t end>, the interval in which T_dep at <t end> lies: the deposit is
      irreversible and only the gas inventory G left at <time> can still
      deposit, so every bin ends between n'_w and n'_w + G/(dx); the first
      crossing of fraction*peak moves upstream at most with every bin at its
      upper bound and the threshold of the lower peak, downstream at most
      with no gain and the threshold of the upper peak.  The gate is the
      one restated by the researcher (2026-10-02, plan section 26): both
      ends of the interval within max(tol, dT_bin) of <expected>, dT_bin
      the drop of the wall temperature across the native bin that holds
      the onset (binDrop); prints the raw difference, dT_bin, the verdict,
      and the plain +-tol comparison for information
  tdepDecided <case> <time> <t end> <expected> <tol> <label>
      as 'tdep' for a run that has not reached <t end>: exit 0 when the
      interval lies within the gate of <expected> (the criterion is
      decided: met), 2 when it lies entirely outside (decided: not met), 1
      when it straddles a limit or is not yet bounded (undecided)
  positions <case> <gas time> <solid time> <paper values ...>
      the positions of Figs. 9 and 10 of the paper against the run (not a
      gate; exit 0 when they could be computed): at <gas time> the gas
      plateau (median of the native bins in 0.30-0.42 m) and where the gas
      falls to 50 % and 10 % of it downstream (linear between the bin
      centres); at <solid time> the onset of the wall column (T_dep
      fromInlet, x and T) and its peak (height and position).  The paper's
      values are key=value: gas50, gas10, solidOnset, solidPeakX,
      solidPeak [m, mol/m]; the differences are printed in mm
  exactLayer <case> <patches> <distance> <log>
      the cells whose centre lies within <distance> of the faces of the
      patches (comma-separated), computed here from constant/polyMesh and
      the cell centres of 0/C (postProcess writeCellCentres) with the exact
      point-polygon distance (the triangles of every face about its
      centre), against the start-up line of the solver log: the same
      number of cells
  perturbWedge <case> <mode>
      moves the points of the wedge mesh of tests/paperMode/wedge (R 2.4
      mm, L 20 mm) in constant/polyMesh/points (ascii or binary, rewritten
      in its format).  interior: the inner points radially, r' = r (1 +
      0.15 sin(pi r/R) sin(2 pi x/L)), and axially, x' = x + 1.5e-4 sin(pi
      x/L) sin(pi r/R) (the axis, the wall, the inlet and the outlet stay;
      a radial move keeps a point on its wedge plane), so the inner faces
      are no longer parallel or normal to the axis while the wall stays
      parallel to it; wall: every point radially by the factor 1 - 0.1 (1 -
      cos(2 pi x/L)), so the wall faces are tilted against the axis
  faceCentreRule <case> <inflow>
      the continuity of the face-centre flux u(C_f) (S_f . x) of the
      parabolic profile (U_b from the carrier log) on the mesh of the case:
      prints max |sum_f phi| per cell relative to <inflow>; exit 0 when it
      exceeds 1e-10 (the mesh is not extruded along the axis: negative
      evidence that the exact integral is needed)
  wallFlux <log> <max> [<min>]
      the 'Wall flux' lines of a setFrozenCarrier log: the sum over the
      wall faces of |phi| relative to the inflow is at most <max> (and at
      least <min>, for the negative evidence)
  mdepPatches <field> <zero patches> <positive patches> [<field B>]
      mDep_ of a written time: every face of the zero patches holds exactly
      0, every face of the positive patches more than 0; with <field B>,
      the values of every listed patch equal those of field B bit for bit
      (the same layer in two boundary orders)
  samePatches <field A> <field B> <patches> [<tol>]
      the values of the listed patches (comma-separated) of two written
      fields are equal bit for bit, or with <tol> to tol times the largest
      value of field B
  patchValue <field> <patch> <index>
      prints the value of one face of a patch ('0' for an exact zero)
  boundWidth <case> <time> <t end> <width>
      the interval of 'tdep' is at most <width> K wide (0 when <time> is the
      end): exit 0, else 1 (the driver extends the run)
  closureRuns <tol> <case> ...
      the closure of every step in every balance file of each run (one
      file per start time of a run restarted in segments): max
      |closure|/reference per run, every one at most tol (M4.17: the drift
      gate of the long runs; review of M4, round 2)
  sameBalance <file A> <file B> <excluded columns>
      two balance files (postProcessing/phaseChangeBalance) hold the same
      header lines and rows, and every column but the excluded ones
      (comma-separated names, e.g. solverDefect,clamped,closure: the
      booking that a change of the ledger alters on purpose) equal as
      printed, i.e. bit for bit; prints in how many rows each excluded
      column differs (the revalidation of the long runs, plan section 27)
  mdepCount <log> <remote elements> <faces with two or more ranks>
      the start-up line of layer distance: '<n> elements with their nearest
      wall face on another rank, <m> interface faces with elements of two
      or more other ranks' with the given numbers (review of M4, round 3:
      the count of round 2 took a face with the elements of one other rank
      twice when two ranks map elements to each other's faces)
  explicitRemainder <log> <min lines> [nonzero]
      the debug report of balance { reportMatrixResidual yes; } on the
      transport line: every 'explicit remainder <E> kg/s in <n> cells (of
      an explicit source sum |b - s| <S> kg/s)' has E = 0 and n = 0 while
      S > 0 (the face-flux correction of the transport-only matrix is
      present and enters the defect by its face values: matrixDefect.H),
      in at least <min lines> lines; with 'nonzero', the positive control:
      some line has n > 0 (a correction in source form, e.g. linearUpwind)
  transportIdentity <log> <min lines> <limit> [violated]
      the transport identity of the same report (matrixDefect.H, plan
      section 33): every 'transport identity <D> kg/s (|D| <r> of its
      terms, <T> kg/s)' has r <= <limit>, in at least <min lines> lines:
      the cell defects in flux form telescope exactly, so D is 0 to the
      rounding of exact sums (measured at most 1.1e-36 of the terms); with
      'violated', the negative control: some line has r > <limit> (a cell
      defect with the face-flux correction in source form leaves its
      double-precision residue, measured 5e-21 to 8e-20 of the terms)
  sourceFormMutant <matrixDefect.H>
      edits a COPY of matrixDefect.H (the negative control of M4.18) so that
      cellDefectT builds the cell defect with the face-flux correction in
      source form, as round 2 did: r_c starts from b_c - s_c (the
      correction's source contribution c_c included) instead of E_c, and
      the face values C_f of internal, coupled and non-coupled faces are
      left out, while the reported explicit remainder E is unchanged; it
      stops (exit 1) unless every edited statement is found exactly once
      (twice for the two boundary corrections)
  moleFraction <case> <time> <log> [<tol>]
      the passive-carrier monitor: the 'max mole fraction' of the step of
      <time> in the log equals the maximum over the cells of c/(c + p/(R
      T)), c = rho Y/M, recomputed from the fields written at <time> (R =
      8.314462618, M of PbI2 = 0.46100894 kg/mol) to tol relative (default
      1e-5, the printed digits); prints the old form (Y/M)/(Y/M +
      1/M_carrier) for comparison (review of M4, round 3)
  inflowExact <phi> <patch> <rho> <flow [mL/min]> <sector [deg]> <tol>
      the inflow through a patch of a written phi, summed exactly (rational
      arithmetic), against rho Q sector/360 with Q and the sector the exact
      decimals given: |inflow/target - 1| <= tol
  restartExact <log> <tol>
      the restart line of a solver log ('Restart of pair ...: the inventory
      re-read differs from the ledger by ... (relative r; <sources>;
      <verdict>)') has the verdict 'exact' and |r| <= tol, and the run
      ended
  outExcerpt ...  (see tests/regress.sh for the excerpts)
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'utilities'))
import axialProfiles as ap   # noqa: E402


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def native(case, time):
    path = os.path.join(case, 'postProcessing', 'axialProfiles', time,
                        'PbI2_g_native.dat')
    if not os.path.exists(path):
        done(False, f'no profile file {path}')
    return ap.read(path)


# -----------------------------------------------------------------------------
# carrier
# -----------------------------------------------------------------------------

def cmdContinuity(log, tol, minLines='1'):
    tol = float(tol)
    found = []
    for line in open(log):
        m = re.search(r'Continuity \((.*?)\): max over the cells of \|sum_f phi\| '
                      r'= (\S+) kg/s = (\S+) of the inflow', line)
        if m:
            found.append((m.group(1), float(m.group(3))))
    ok = len(found) >= int(minLines) and all(r <= tol for _, r in found)
    done(ok, f'{log}: ' + '; '.join(f'{what}: {r:.3g} of the inflow'
                                     for what, r in found)
         + f' (limit {tol:g})')


def cmdContinuityAbove(log, minimum):
    minimum = float(minimum)
    worst = 0.0
    for line in open(log):
        m = re.search(r'Continuity \((.*?)\): .*?= (\S+) of the inflow', line)
        if m:
            worst = max(worst, float(m.group(2)))
    done(worst > minimum, f'{log}: largest continuity error {worst:.3g} of the '
         f'inflow (must exceed {minimum:g}: the check sees the rounding)')


def carrierValues(log):
    text = open(log).read()
    out = {}
    m = re.search(r'U_b = Q/\(pi R\^2\) = (\S+) m/s', text)
    out['UbAnalytic'] = float(m.group(1)) if m else None
    m = re.search(r'rescaled to the facets by \S+: U_b = (\S+) m/s', text)
    out['Ub'] = float(m.group(1)) if m else out['UbAnalytic']
    m = re.search(r'discrete inflow through \S+ = (\S+) m3/s = (\S+) of Q sector',
                  text)
    if m:
        out['inflow'], out['inflowRatio'] = float(m.group(1)), float(m.group(2))
    m = re.search(r'inlet area (\S+) m2, facet-area ratio \(inlet area/\(pi R\^2 '
                  r'sector\)\) (\S+); discrete bulk velocity (\S+) m/s', text)
    if m:
        out['inletArea'] = float(m.group(1))
        out['facetRatio'] = float(m.group(2))
        out['bulk'] = float(m.group(3))
    return out


def cmdCarrier(log, *checks):
    values = carrierValues(log)
    ok = True
    parts = []
    for check in checks:
        key, expected = check.split('=')
        expected = float(expected)
        value = values.get(key)
        good = value is not None and abs(value - expected) <= 1e-9*abs(expected)
        ok = ok and good
        parts.append(f'{key} {value!r} (expected {expected!r})')
    done(ok, f'{log}: ' + ', '.join(parts))


def cmdLayerCounts(log, rule, cells, perVolume, tol):
    cells, perVolume, tol = int(cells), float(perVolume), float(tol)
    for line in open(log):
        m = re.search(r'Interface \(layer ' + rule + r'\): .*?(\d+) elements in '
                      r'(\d+) cells; A_I = (\S+) m2, V_I = (\S+) m3, A_I/V_I = '
                      r'(\S+) 1/m', line)
        if m:
            n, nc = int(m.group(1)), int(m.group(2))
            av = float(m.group(5))
            ok = nc == cells and abs(av/perVolume - 1) <= tol
            done(ok, f'layer {rule}: {n} elements in {nc} cells (expected '
                 f'{cells}), A_I/V_I = {av} 1/m (expected {perVolume} to '
                 f'{tol:g})')
    done(False, f'no line "Interface (layer {rule})" in {log}')


# -----------------------------------------------------------------------------
# mDep_ and the nearest-face map
# -----------------------------------------------------------------------------

def readBalance(path):
    names, rows = None, []
    for line in open(path):
        if line.startswith('# time'):
            names = line[1:].split()
        elif line.startswith('#') or not line.strip():
            continue
        else:
            parts = line.split()
            rows.append([float(x) for x in parts])
    return names, rows


def readFunctionObject(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[-1])))
    return rows


def cmdMdep(dat, foDat, tol):
    tol = float(tol)
    names, rows = readBalance(dat)
    iw = names.index('wall')
    wall = {round(r[0], 9): r[iw] for r in rows}
    fo = readFunctionObject(foDat)
    scale = max(abs(v) for v in wall.values())
    worst, n = 0.0, 0
    for t, value in fo:
        key = round(t, 9)
        if key in wall:
            worst = max(worst, abs(value - wall[key]))
            n += 1
    ok = n > 0 and scale > 0 and worst <= tol*scale
    done(ok, f'areaIntegrate of mDep_ on the wall vs the wall inventory over '
         f'{n} steps: max difference {worst:.3e} kg = {worst/max(scale, 1e-300):.2e}'
         f' of the largest wall inventory {scale:.6e} kg (limit {tol:g})')


def cmdCompareInventories(datA, datB, tol, *cols):
    """the columns of two balance files at every common step, relative to
    the ledger's reference of B (max(|initial + released|, |held|, inflow))"""
    pc = phaseChangeChecks()
    tol = float(tol)
    namesA, rowsA = pc.readBalance(datA)
    namesB, rowsB = pc.readBalance(datB)
    byTime = {r[0]: r for r in rowsA}
    worst = {c: (0.0, None) for c in cols}
    n = 0
    for row in rowsB:
        other = byTime.get(row[0])
        if other is None:
            continue
        n += 1
        ref = pc.reference(namesB, row)
        for c in cols:
            d = abs(pc.column(namesA, other, c) - pc.column(namesB, row, c))/ref
            if d > worst[c][0]:
                worst[c] = (d, row[0])
    ok = n == len(rowsB) == len(rowsA) and all(w <= tol for w, _ in
                                               worst.values())
    done(ok, f'{n} common steps of {len(rowsA)} and {len(rowsB)}: max |A - B|/'
         f'reference ' + ', '.join(f'{c} {w:.2e} (t = {t})'
                                   for c, (w, t) in worst.items())
         + f' (limit {tol:g})')


def cmdFacesAgree(fieldA, fieldB, patch, tol):
    sys.path.insert(0, os.path.join(HERE, '..', 'phaseChange'))
    import checks as pc
    tol = float(tol)
    _, pa = pc.readScalarField(fieldA)
    _, pb = pc.readScalarField(fieldB)
    a, b = pa[patch], pb[patch]
    scale = max(max(abs(v) for v in b), 1e-300)
    worst = max(abs(x - y) for x, y in zip(a, b))
    half = len(b)//2
    sym = max(abs(x - y) for x, y in zip(b[:half], b[half:]))
    ok = len(a) == len(b) and worst <= tol*scale and sym <= tol*scale
    done(ok, f'{patch} of {fieldA} vs {fieldB}: {len(a)} faces, max difference '
         f'{worst/scale:.2e} of the largest value {scale:.6e}; the two walls of '
         f'the symmetric channel differ by {sym/scale:.2e} (limit {tol:g})')


# -----------------------------------------------------------------------------
# Fig. 9 identity and agreement of the solid profiles
# -----------------------------------------------------------------------------

def cmdPlateau(case, time, tev, n0, xlo, xhi, tol, UbD, Ub, paper):
    tev, n0, xlo, xhi, tol = map(float, (tev, n0, xlo, xhi, tol))
    UbD, Ub, paper = float(UbD), float(Ub), float(paper)
    p = native(case, time)
    x, gas = p.column('x'), p.column('gas')
    target = n0/UbD
    ratios = [g*tev/target for xc, g in zip(x, gas) if xlo <= xc <= xhi]
    if not ratios:
        done(False, f'no bin in [{xlo}, {xhi}] m')
    worst = max(abs(r - 1) for r in ratios)
    mean = math.fsum(ratios)/len(ratios)*target
    done(worst <= tol,
         f'{case} at {time} s (t_ev {tev:g} s): n\' t_ev on {len(ratios)} '
         f'native bins in [{xlo}, {xhi}] m = {mean:.6e} mol s/m (from '
         f'{min(ratios)*target:.6e} to {max(ratios)*target:.6e}); n0/U_b,d = '
         f'{target:.6e} (U_b,d = {UbD:.8g} m/s, the discrete bulk velocity of '
         f'the carrier): max |ratio - 1| = {worst:.2e} (limit {tol:g}); '
         f'n0/U_b = {n0/Ub:.6e} with the analytic U_b = {Ub:.8g} m/s; the '
         f'paper: {paper:g} (plateau/paper - 1 = {mean/paper - 1:+.3f})')


def tdepOf(p, values, fraction):
    r = ap.onset(p, False, values=values, fraction=fraction)
    return r


def cmdSolidAgree(caseA, timeA, caseB, timeB, tol):
    tol = float(tol)
    a, b = native(caseA, timeA), native(caseB, timeB)
    wa, wb = a.column('wall'), b.column('wall')
    dx = [hi - lo for lo, hi in zip(b.column('x_lo'), b.column('x_hi'))]
    l1 = math.fsum(abs(x - y)*d for x, y, d in zip(wa, wb, dx))
    norm = math.fsum(abs(y)*d for y, d in zip(wb, dx))
    peak = max(wb)
    worst = max(abs(x - y) for x, y in zip(wa, wb))
    ta = tdepOf(a, wa, a.fraction)
    tb = tdepOf(b, wb, b.fraction)
    text = (f'solid profiles {caseA}@{timeA} vs {caseB}@{timeB}: L1 difference '
            f'{l1/norm:.3e} of the deposit of B (limit {tol:g}); largest bin '
            f'difference {worst/peak:.3e} of the peak {peak:.6e} mol/m; '
            f'deposits {a.inventory.get("wall"):.6e} and '
            f'{b.inventory.get("wall"):.6e} mol, gas left '
            f'{a.inventory.get("gas"):.3e} and {b.inventory.get("gas"):.3e} mol; '
            f'T_dep fromInlet {ta.get("T")} and {tb.get("T")} K')
    done(l1 <= tol*norm, text)


def releaseEnd(case):
    """the end of the last release window of the case's samples [s]
    (startTime + duration of every sample with a duration), 0 without"""
    path = os.path.join(case, 'constant', 'thermochemistryProperties')
    text = re.sub(r'/\*.*?\*/', ' ', open(path).read(), flags=re.S)
    text = re.sub(r'//[^\n]*', ' ', text)
    end = 0.0
    for m in re.finditer(r'samples\s*\{(.*)', text, re.S):
        body = m.group(1)
        starts = [float(x) for x in re.findall(r'startTime\s+([^;\s]+)\s*;', body)]
        durations = [float(x) for x in re.findall(r'duration\s+([^;\s]+)\s*;', body)]
        for a, b in zip(starts, durations):
            end = max(end, a + b)
    return end


def tdepInterval(p, final):
    """T_dep fromInlet of the wall column at the file's time, and the
    interval of T_dep at a later time (irreversible deposit, gas left G)."""
    w = p.column('wall')
    f = p.fraction
    point = tdepOf(p, w, f)
    if final:
        return point, point.get('T'), point.get('T')
    G = p.inventory.get('gas', 0.0)
    dx = [hi - lo for lo, hi in zip(p.column('x_lo'), p.column('x_hi'))]
    upper = [v + G/d for v, d in zip(w, dx)]
    peakLow, peakHigh = max(w), max(upper)
    if not peakLow > 0 or f*peakHigh > peakLow:
        # no deposit yet, or the remaining gas could make a new peak whose
        # threshold no bin reaches now: the onset is not bounded
        return point, None, None
    early = tdepOf(p, upper, f*peakLow/peakHigh)
    late = tdepOf(p, w, f*peakHigh/peakLow)
    Ts = [r.get('T') for r in (early, late)]
    if any(T is None for T in Ts):
        return point, None, None
    return point, min(Ts), max(Ts)


def binDrop(p, x):
    """The drop of the wall temperature [K] across the bin that holds x,
    with the wall temperature at its edges interpolated linearly between
    the bin centres (the Twall column; extrapolated with the end
    intervals), and the bin's index."""
    lo, hi, xc, T = (p.column('x_lo'), p.column('x_hi'), p.column('x'),
                     p.column('Twall'))
    n = len(xc)
    b = next((i for i in range(n) if lo[i] <= x < hi[i]), n - 1)

    def wallT(z):
        j = 1
        while j < n - 1 and xc[j] < z:
            j += 1
        return T[j - 1] + (z - xc[j - 1])*(T[j] - T[j - 1])/(xc[j] - xc[j - 1])
    return abs(wallT(lo[b]) - wallT(hi[b])), b


def gateText(lo, hi, expected, tol, dT):
    """The verdict of the restated gate of acceptance 5 (max(tol, dT_bin))
    and the plain +-tol comparison, for both ends of an interval."""
    gate = max(tol, dT)
    ok = abs(lo - expected) <= gate and abs(hi - expected) <= gate
    plain = abs(lo - expected) <= tol and abs(hi - expected) <= tol
    diff = (f'{lo - expected:+.2f} K' if lo == hi
            else f'{lo - expected:+.2f} to {hi - expected:+.2f} K')
    return ok, (f'difference {diff}; dT_bin {dT:.2f} K (the drop of the wall '
                f'temperature across the native bin of the onset); gate '
                f'max({tol:g}, dT_bin) = {gate:.2f} K: '
                f'{"met" if ok else "not met"}; plain +-{tol:g} K (for '
                f'information): {"within" if plain else "outside"}')


def cmdTdep(case, time, tEnd, expected, tol, label):
    expected, tol = float(expected), float(tol)
    p = native(case, time)
    final = abs(float(time) - float(tEnd)) <= 1e-9*max(1.0, float(tEnd))
    point, lo, hi = tdepInterval(p, final)
    if point.get('status') != 'found' or lo is None:
        done(False, f'{label}: T_dep fromInlet at {time} s: {point}')
    dT, _ = binDrop(p, point['x'])
    ok, verdict = gateText(lo, hi, expected, tol, dT)
    if final:
        what = f'T_dep {point["T"]:.4f} K at x = {point["x"]:.6f} m at the end ' \
               f'of the experiment, {time} s'
    else:
        what = (f'T_dep {point["T"]:.4f} K at x = {point["x"]:.6f} m at '
                f'{time} s with {p.inventory.get("gas"):.3e} mol of gas left '
                f'(of {p.inventory.get("gas") + p.inventory.get("wall"):.4e}); '
                f'at {tEnd} s it lies in [{lo:.4f}, {hi:.4f}] K (the deposit '
                f'is irreversible)')
    done(ok, f'{label}: {what}; the paper (Table 2): {expected:g} K, '
         f'{verdict}')


def cmdTdepDecided(case, time, tEnd, expected, tol, label):
    expected, tol = float(expected), float(tol)
    if float(time) < releaseEnd(case):
        print(f'{label}: at {time} s the release (until {releaseEnd(case):g} s) '
              f'has not ended, so T_dep at {tEnd} s is not bounded yet')
        sys.exit(1)
    p = native(case, time)
    final = abs(float(time) - float(tEnd)) <= 1e-9*max(1.0, float(tEnd))
    point, lo, hi = tdepInterval(p, final)
    gas = p.inventory.get('gas')
    if point.get('status') != 'found' or lo is None:
        print(f'{label}: at {time} s T_dep is not bounded yet ({gas:.3e} mol '
              f'of gas left; {point.get("status")})')
        sys.exit(1)
    dT, _ = binDrop(p, point['x'])
    gate = max(tol, dT)
    inside, verdict = gateText(lo, hi, expected, tol, dT)
    outside = hi < expected - gate or lo > expected + gate
    state = 'met' if inside else 'not met' if outside else 'undecided'
    print(f'{label}: at {time} s with {gas:.3e} mol of gas left, T_dep at '
          f'{tEnd} s lies in [{lo:.4f}, {hi:.4f}] K (the deposit is '
          f'irreversible; now {point["T"]:.4f} K at x = {point["x"]:.6f} m); '
          f'the paper (Table 2): {expected:g} K, {verdict}; decided: {state}')
    sys.exit(0 if inside else 2 if outside else 1)


def cmdBoundWidth(case, time, tEnd, width):
    if float(time) < releaseEnd(case):
        print(f'{case} at {time} s: the release has not ended')
        sys.exit(1)
    p = native(case, time)
    final = abs(float(time) - float(tEnd)) <= 1e-9*max(1.0, float(tEnd))
    point, lo, hi = tdepInterval(p, final)
    w = None if lo is None else hi - lo
    print(f'{case} at {time} s: T_dep interval width {w} K '
          f'(gas left {p.inventory.get("gas"):.3e} mol)')
    sys.exit(0 if w is not None and w <= float(width) else 1)


# -----------------------------------------------------------------------------
# carrier files and provenance
# -----------------------------------------------------------------------------

def phaseChangeChecks():
    sys.path.insert(0, os.path.join(HERE, '..', 'phaseChange'))
    import checks as pc
    return pc


def cmdSamePhi(caseA, caseB, nProcs):
    pc = phaseChangeChecks()
    worst, nValues = 0.0, 0
    for i in range(int(nProcs)):
        a = pc.readScalarField(os.path.join(caseA, f'processor{i}', '0', 'phi'))
        b = pc.readScalarField(os.path.join(caseB, f'processor{i}', '0', 'phi'))
        pairs = [(a[0], b[0])] + [(a[1][k], b[1][k]) for k in a[1] if k in b[1]]
        if set(a[1]) != set(b[1]):
            done(False, f'processor{i}: patches differ: {sorted(a[1])} and '
                 f'{sorted(b[1])}')
        for x, y in pairs:
            if len(x) != len(y):
                done(False, f'processor{i}: lists of {len(x)} and {len(y)} values')
            for u, v in zip(x, y):
                worst = max(worst, abs(u - v))
                nValues += 1
    done(worst == 0.0, f'phi written on {nProcs} ranks equals the decomposed '
         f'serial phi bit for bit: {nValues} values (internal and patch '
         f'faces), largest difference {worst:g} kg/s')


def cmdHeldWall(case, x0, T0):
    pc = phaseChangeChecks()
    T0 = float(T0)
    _, patches = pc.readScalarField(os.path.join(case, '0', 'T'))
    wall = patches['WALL']
    held = sum(1 for v in wall if v == T0)
    ok = min(wall) == T0 and held > 0
    done(ok, f'outOfRange hold: the wall T beyond x = {x0} m is the last row '
         f'read, {T0:g} K, on {held} faces (the minimum of the wall {min(wall):g} '
         f'K; with all 70 rows it would be 305 K)')


def cmdSameWallColumns(caseA, caseB):
    files = ap.files(caseA)
    if not files:
        done(False, f'no profile files in {caseA}')
    worst, n = 0.0, 0
    for path in files:
        other = os.path.join(caseB, os.path.relpath(path, caseA))
        a, b = ap.read(path), ap.read(other)
        for name in ('Twall', 'wallArea'):
            for u, v in zip(a.column(name), b.column(name)):
                worst = max(worst, abs(u - v))
                n += 1
    done(worst == 0.0, f'Twall and wallArea of {len(files)} profile files '
         f'equal bit for bit to those of {caseB} ({n} values): the wall '
         f'temperature of the interface faces for both layer rules')


def cmdConductivityTable(dictPath):
    """Types.f90 on stdin: k_he at t_values_k against the conductivity list
    of a setFrozenCarrierDict."""
    source = sys.stdin.read()

    def array(name):
        m = re.search(name + r'\s*&[^\n]*\n(.*?)/\)', source, re.S)
        body = re.sub(r'![^\n]*', '', m.group(1)).replace('(/', '')
        return [float(x) for x in
                re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', body)]

    k, T = array('k_he'), array('t_values_k')
    text = open(dictPath).read()
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    m = re.search(r'conductivity\s*\((.*?)\)\s*;', text, re.S)
    pairs = re.findall(r'\(\s*([^\s()]+)\s+([^\s()]+)\s*\)', m.group(1))
    table = [(float(a), float(b)) for a, b in pairs]
    ok = len(table) == len(k) == len(T) == 51 \
        and all(x == (t, v) for x, t, v in zip(table, T, k))
    done(ok, f'the conductivity table of {dictPath} ({len(table)} rows) is '
         f'k_he at t_values_k of User_Mod/Types.f90 ({len(k)} values)')


def cmdRoundProcessorFace(case, digits):
    """negative evidence of the decomposed continuity check: the first value
    of the first processor patch of processor0/0/phi (an ascii file) is
    rewritten with <digits> significant digits, as decomposePar writes a
    'uniform' processor patch; prints the relative change"""
    path = os.path.join(case, 'processor0', '0', 'phi')
    text = open(path).read()
    m = re.search(r'(procBoundary\w+\s*\{[^}]*?value\s+nonuniform\s+'
                  r'List<scalar>\s*\d+\s*\(\s*)(\S+)', text)
    if not m:
        done(False, f'no processor patch with values in {path} (ascii?)')
    old = float(m.group(2))
    new = float('%.*g' % (int(digits), old))
    text = text[:m.start(2)] + repr(new) + text[m.end(2):]
    open(path, 'w').write(text)
    done(new != old, f'{path}: the first value of the first processor patch '
         f'{old!r} -> {new!r} ({digits} digits, relative change '
         f'{(new - old)/old:.2e})')


def cmdRestartExact(log, tol):
    tol = float(tol)
    text = open(log, errors='replace').read()
    m = re.search(r'Restart of pair (\S+): the inventory re-read differs from '
                  r'the ledger by (\S+) kg \(relative (\S+); (.*?); '
                  r'(exact|rounded|inventory exact, continuation rounded)\), '
                  r'booked as restart', text)
    if not m:
        done(False, f'no restart line in {log}')
    relative = float(m.group(3))
    ended = re.search(r'^End\s*$', text, re.M) is not None
    done(m.group(5) == 'exact' and abs(relative) <= tol and ended,
         f'{os.path.basename(os.path.dirname(os.path.abspath(log)))}: '
         f'{m.group(1)} re-read with {m.group(2)} kg (relative {relative:.3g}, '
         f'limit {tol:g}; {m.group(4)}; {m.group(5)}), '
         f'{"ended" if ended else "did not end"}')


# -----------------------------------------------------------------------------
# Figs. 9 and 10 (not a gate), the exact layer, perturbed wedges
# -----------------------------------------------------------------------------

def cmdPositions(case, gasTime, solidTime, *paper):
    ref = dict(item.split('=') for item in paper)
    ref = {k: float(v) for k, v in ref.items()}
    g = native(case, gasTime)
    x, gas = g.column('x'), g.column('gas')
    inner = sorted(v for xc, v in zip(x, gas) if 0.30 <= xc <= 0.42)
    plateau = inner[len(inner)//2]
    start = max(i for i, xc in enumerate(x) if xc <= 0.42)

    def falls(frac):
        level = frac*plateau
        for i in range(start + 1, len(x)):
            if gas[i] < level <= gas[i - 1]:
                return x[i - 1] + (level - gas[i - 1])*(x[i] - x[i - 1]) \
                    /(gas[i] - gas[i - 1])
        return None
    gas50, gas10 = falls(0.5), falls(0.1)
    s = native(case, solidTime)
    w = s.column('wall')
    onset = tdepOf(s, w, s.fraction)
    peakBin = max(range(len(w)), key=lambda i: w[i])
    xs = s.column('x')

    def mm(value, key):
        if value is None or key not in ref:
            return ''
        return f' (paper {ref[key]:g}: {1e3*(value - ref[key]):+.1f} mm)'
    text = (f'{case}: gas at {gasTime} s: plateau {plateau:.4e} mol/m, 50 % '
            f'at {gas50:.4f} m{mm(gas50, "gas50")}, 10 % at {gas10:.4f} m'
            f'{mm(gas10, "gas10")}; solid at {solidTime} s: onset (1 %) at '
            f'{onset.get("x", float("nan")):.4f} m'
            f'{mm(onset.get("x"), "solidOnset")}, T_dep '
            f'{onset.get("T", float("nan")):.2f} K; peak '
            f'{w[peakBin]:.3e} mol/m at {xs[peakBin]:.4f} m'
            f'{mm(xs[peakBin], "solidPeakX")}'
            + (f' (paper {ref["solidPeak"]:.3g} mol/m)'
               if 'solidPeak' in ref else ''))
    done(gas50 is not None and gas10 is not None
         and onset.get('status') == 'found', text)


def hksChecks():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'hkschecks', os.path.join(HERE, '..', 'hks', 'checks.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def patchRange(case, patch):
    """(startFace, nFaces) of a patch of constant/polyMesh/boundary"""
    text = open(os.path.join(case, 'constant', 'polyMesh', 'boundary')).read()
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    m = re.search(r'\b' + re.escape(patch) + r'\s*\{([^}]*)\}', text)
    if not m:
        done(False, f'no patch {patch} in {case}/constant/polyMesh/boundary')
    return (int(re.search(r'startFace\s+(\d+)', m.group(1)).group(1)),
            int(re.search(r'nFaces\s+(\d+)', m.group(1)).group(1)))


def triangleDistance(p, a, b, c):
    """The distance of point p from the triangle (a, b, c) (Ericson)."""
    sub = lambda u, v: (u[0] - v[0], u[1] - v[1], u[2] - v[2])
    dot = lambda u, v: u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    ab, ac, ap = sub(b, a), sub(c, a), sub(p, a)
    d1, d2 = dot(ab, ap), dot(ac, ap)
    if d1 <= 0 and d2 <= 0:
        q = a
    else:
        bp = sub(p, b)
        d3, d4 = dot(ab, bp), dot(ac, bp)
        cp = sub(p, c)
        d5, d6 = dot(ab, cp), dot(ac, cp)
        vc, vb, va = d1*d4 - d3*d2, d5*d2 - d1*d6, d3*d6 - d5*d4
        if d3 >= 0 and d4 <= d3:
            q = b
        elif vc <= 0 and d1 >= 0 and d3 <= 0:
            t = d1/(d1 - d3)
            q = tuple(a[i] + t*ab[i] for i in range(3))
        elif d6 >= 0 and d5 <= d6:
            q = c
        elif vb <= 0 and d2 >= 0 and d6 <= 0:
            t = d2/(d2 - d6)
            q = tuple(a[i] + t*ac[i] for i in range(3))
        elif va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
            t = (d4 - d3)/((d4 - d3) + (d5 - d6))
            q = tuple(b[i] + t*(c[i] - b[i]) for i in range(3))
        else:
            den = 1/(va + vb + vc)
            v, w = vb*den, vc*den
            q = tuple(a[i] + v*ab[i] + w*ac[i] for i in range(3))
    return math.sqrt(sum((p[i] - q[i])**2 for i in range(3)))


def faceCentre(points, face):
    """The area-weighted centre of a polygon (triangles about the mean)."""
    n = len(face)
    if n == 3:
        return tuple(sum(points[v][i] for v in face)/3 for i in range(3))
    m = tuple(sum(points[v][i] for v in face)/n for i in range(3))
    total, acc = 0.0, [0.0, 0.0, 0.0]
    for j in range(n):
        a, b = points[face[j]], points[face[(j + 1) % n]]
        u = [a[i] - m[i] for i in range(3)]
        v = [b[i] - m[i] for i in range(3)]
        c = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2],
             u[0]*v[1] - u[1]*v[0])
        area = math.sqrt(sum(x*x for x in c))
        for i in range(3):
            acc[i] += area*(a[i] + b[i] + m[i])/3
        total += area
    return tuple(acc[i]/total for i in range(3))


def cmdExactLayer(case, patches, distance, log):
    hk = hksChecks()
    pc = phaseChangeChecks()
    distance = float(distance)
    poly = os.path.join(case, 'constant', 'polyMesh')
    points = hk.readPoints(os.path.join(poly, 'points'))
    faces = hk.readFaces(os.path.join(poly, 'faces'))
    tris = []
    for name in patches.split(','):
        start, n = patchRange(case, name)
        for f in range(start, start + n):
            face = faces[f]
            if len(face) == 3:
                tris.append(tuple(points[v] for v in face))
            else:
                ctr = faceCentre(points, face)
                for j in range(len(face)):
                    tris.append((points[face[j]],
                                 points[face[(j + 1) % len(face)]], ctr))
    C = [pc.readScalarField(os.path.join(case, '0', 'C' + x))[0]
         for x in 'xyz']
    # bins of the triangles (a uniform grid of the search radius)
    h = distance
    grid = {}
    for k, (a, b, c) in enumerate(tris):
        lo = [min(a[i], b[i], c[i]) - h for i in range(3)]
        hi = [max(a[i], b[i], c[i]) + h for i in range(3)]
        rng = [range(int(math.floor(lo[i]/h)), int(math.floor(hi[i]/h)) + 1)
               for i in range(3)]
        for i in rng[0]:
            for j in rng[1]:
                for l in rng[2]:
                    grid.setdefault((i, j, l), []).append(k)
    selected, dmax = 0, 0.0
    for p in zip(*C):
        key = tuple(int(math.floor(p[i]/h)) for i in range(3))
        best = min((triangleDistance(p, *tris[k]) for k in grid.get(key, [])),
                   default=math.inf)
        if best <= distance:
            selected += 1
            dmax = max(dmax, best)
    text = open(log, errors='replace').read()
    m = re.search(r'Interface \(layer distance\): .*?(\d+) elements in (\d+) '
                  r'cells', text)
    solver = int(m.group(2)) if m else None
    done(solver == selected,
         f'{case}: {selected} cells within {distance:g} m of {patches} '
         f'(the exact point-polygon distance from constant/polyMesh, '
         f'{len(tris)} triangles; largest selected {dmax:.6g} m); the solver '
         f'({os.path.basename(log)}): {solver}')


def cmdPerturbWedge(case, mode):
    path = os.path.join(case, 'constant', 'polyMesh', 'points')
    data = open(path, 'rb').read()
    header = data.index(b'}', data.index(b'FoamFile')) + 1
    binary = phaseChangeChecks().foamFormat(data) == 'binary'
    hk = hksChecks()
    pts = hk.readPoints(path)
    R, L = 2.4e-3, 0.02
    out, moved = [], 0
    for x, y, z in pts:
        r = math.hypot(y, z)
        nx, scale = x, 1.0
        if mode == 'interior':
            inner = r < R*(1 - 1e-9) and 1e-12 < x < L - 1e-12 and r > 0
            if inner:
                scale = 1 + 0.15*math.sin(math.pi*r/R)*math.sin(2*math.pi*x/L)
                nx = x + 1.5e-4*math.sin(math.pi*x/L)*math.sin(math.pi*r/R)
        elif mode == 'wall':
            scale = 1 - 0.1*(1 - math.cos(2*math.pi*x/L))
        else:
            done(False, f'unknown mode {mode}')
        if scale != 1.0 or nx != x:
            moved += 1
        out.append((nx, y*scale, z*scale))
    m = re.compile(rb'\n\s*(\d+)\s*\(').search(data, header)
    n, start = int(m.group(1)), m.end()
    if binary:
        import struct
        body = struct.pack('<%dd' % (3*n), *[v for q in out for v in q])
        data = data[:start] + body + data[start + 24*n:]
    else:
        end = start
        for _ in range(n):
            end = data.index(b')', end) + 1
        body = b'\n' + b'\n'.join(('(%r %r %r)' % q).encode() for q in out)
        data = data[:start] + body + data[end:]
    open(path, 'wb').write(data)
    done(moved > 0, f'{path}: {moved} of {n} points moved ({mode}; '
         f'{"binary" if binary else "ascii"})')


def cmdFaceCentreRule(case, inflow):
    hk = hksChecks()
    inflow = float(inflow)
    log = open(os.path.join(case, 'log.setFrozenCarrier')).read()
    Ub = float(re.search(r'rescaled to the facets by \S+: U_b = (\S+) m/s',
                         log).group(1)) \
        if 'rescaled to the facets' in log else \
        float(re.search(r'U_b = Q/\(pi R\^2\) = (\S+) m/s', log).group(1))
    poly = os.path.join(case, 'constant', 'polyMesh')
    points = hk.readPoints(os.path.join(poly, 'points'))
    faces = hk.readFaces(os.path.join(poly, 'faces'))
    owner = hk.readLabels(os.path.join(poly, 'owner'))
    neighbour = hk.readLabels(os.path.join(poly, 'neighbour'))
    R = 2.4e-3
    total = [0.0]*(max(owner) + 1)
    for f, face in enumerate(faces):
        ctr = faceCentre(points, face)
        p0 = points[face[0]]
        Sx = 0.0
        for a, b in zip(face[1:-1], face[2:]):
            u = [points[a][i] - p0[i] for i in range(3)]
            v = [points[b][i] - p0[i] for i in range(3)]
            Sx += 0.5*(u[1]*v[2] - u[2]*v[1])
        phi = 2*Ub*(1 - (ctr[1]**2 + ctr[2]**2)/R**2)*Sx
        total[owner[f]] += phi
        if f < len(neighbour):
            total[neighbour[f]] -= phi
    worst = max(abs(t) for t in total)
    done(worst > 1e-10*inflow,
         f'{case}: the face-centre flux u(C_f) S_f,x misses continuity by '
         f'{worst/inflow:.2e} of the inflow (max |sum_f phi| per cell; the '
         f'exact integral of setFrozenCarrier: round-off)')


def cmdWallFlux(log, maximum, minimum='-1'):
    maximum, minimum = float(maximum), float(minimum)
    found = []
    for line in open(log):
        m = re.search(r'Wall flux \((.*?)\): sum over the faces of the walls of '
                      r'\|phi\| = (\S+) kg/s = (\S+) of the inflow', line)
        if m:
            found.append((m.group(1), float(m.group(3))))
    ok = found and all(minimum < r <= maximum for _, r in found)
    done(bool(ok), f'{log}: ' + '; '.join(f'{w}: {r:.3g} of the inflow'
                                          for w, r in found)
         + f' (within ({minimum:g}, {maximum:g}])')


def cmdMdepPatches(field, zero, positive, other=None):
    pc = phaseChangeChecks()
    _, patches = pc.readScalarField(field)
    zeros = [z for z in zero.split(',') if z]
    positives = [q for q in positive.split(',') if q]
    ok, parts = True, []
    for name in zeros:
        values = patches[name]
        bad = sum(1 for v in values if v != 0)
        ok = ok and bad == 0
        parts.append(f'{name}: {bad} faces not 0')
    for name in positives:
        values = patches[name]
        bad = sum(1 for v in values if not v > 0)
        ok = ok and bad == 0 and len(values) > 0
        parts.append(f'{name}: {bad} faces not > 0 (of {len(values)})')
    if other:
        _, b = pc.readScalarField(other)
        for name in zeros + positives:
            same = patches[name] == b[name]
            ok = ok and same
            parts.append(f'{name} {"equal" if same else "DIFFERENT"} in '
                         f'{other}')
    done(ok, f'{field}: ' + '; '.join(parts))


def patchList(field, patch):
    pc = phaseChangeChecks()
    _, patches = pc.readScalarField(field)
    return patches[patch]


def cmdSamePatches(fieldA, fieldB, patches, tol=None):
    parts, ok = [], True
    names = patches.split(',')
    lists = {n: (patchList(fieldA, n), patchList(fieldB, n)) for n in names}
    scale = max(max(abs(v) for v in b) for _, b in lists.values())
    for name in names:
        a, b = lists[name]
        if len(a) == 1 and len(b) > 1:
            a = a*len(b)
        if len(b) == 1 and len(a) > 1:
            b = b*len(a)
        if tol is None:
            same = a == b
            parts.append(f'{name} {"equal" if same else "DIFFERENT"}')
        else:
            worst = max(abs(x - y) for x, y in zip(a, b))
            same = len(a) == len(b) and worst <= float(tol)*scale
            parts.append(f'{name} to {worst/max(scale, 1e-300):.1e}')
        ok = ok and same
    limit = 'bit for bit' if tol is None else f'limit {float(tol):g} of the ' \
        f'largest value {scale:.6e}'
    done(ok, f'{fieldA} vs {fieldB}: ' + ', '.join(parts) + f' ({limit})')


def cmdPatchValue(field, patch, index):
    values = patchList(field, patch)
    value = values[0] if len(values) == 1 else values[int(index)]
    print('0' if value == 0 else repr(value))


def cmdClosureRuns(tol, *cases):
    """the closure of every step of every balance file of the runs (a run
    restarted in segments writes one file per start time): max
    |closure|/reference per run, all at most tol"""
    pc = phaseChangeChecks()
    tol = float(tol)
    parts, ok = [], bool(cases)
    for case in cases:
        base = os.path.join(case, 'postProcessing', 'phaseChangeBalance')
        worst, n, where = 0.0, 0, None
        for t0 in sorted(os.listdir(base), key=float):
            names, rows = pc.readBalance(os.path.join(base, t0, 'PbI2_g.dat'))
            for r in rows:
                rel = abs(pc.column(names, r, 'closure'))/max(
                    pc.reference(names, r), 1e-300)
                n += 1
                if rel >= worst:
                    worst, where = rel, r[0]
        ok = ok and n > 0 and worst <= tol
        parts.append(f'{os.path.basename(case.rstrip("/"))} {worst:.1e}'
                     f' ({n} steps, largest at t = {float(where):.6g})')
    done(ok, 'max |closure|/reference per run: ' + '; '.join(parts)
         + f' (limit {tol:g})')


def cmdSameBalance(fileA, fileB, excluded):
    """two balance files: the same header lines and rows, every column but
    the excluded ones identical as printed (17 significant digits, so bit
    for bit)"""
    skip = [c for c in excluded.split(',') if c]

    def read(path):
        header, names, rows = [], None, []
        for line in open(path):
            if line.startswith('#'):
                header.append(line)
                if line.startswith('# time'):
                    names = line[2:].split()
            elif line.strip():
                rows.append(line.split())
        return header, names, rows

    ha, na, ra = read(fileA)
    hb, nb, rb = read(fileB)
    if na is None or ha != hb or len(ra) != len(rb) \
            or any(c not in na for c in skip):
        done(False, f'{fileA} vs {fileB}: different header lines or numbers '
             f'of rows ({len(ra)} and {len(rb)}), or an excluded column '
             f'missing ({", ".join(skip)})')
    kept = [i for i, n in enumerate(na) if n not in skip]
    differing = sum(1 for a, b in zip(ra, rb)
                    if any(a[i] != b[i] for i in kept))
    changed = {c: sum(1 for a, b in zip(ra, rb)
                      if a[na.index(c)] != b[na.index(c)]) for c in skip}
    done(differing == 0,
         f'{len(ra)} rows: the columns {" ".join(na[i] for i in kept)} '
         f'{"identical" if differing == 0 else f"differ in {differing} rows"}'
         f'; left out (rows in which they differ): '
         + ', '.join(f'{c} {n}' for c, n in changed.items()))


# -----------------------------------------------------------------------------
# The nearest-face map with two-way maps, the explicit remainder of the
# transport matrix, the mole-fraction monitor, the exact inflow (review of
# M4, round 3; plan section 29)
# -----------------------------------------------------------------------------

def cmdMdepCount(log, remote, multiple):
    text = open(log, errors='replace').read()
    m = re.search(r'(\d+) elements with their nearest wall face on another '
                  r'rank, (\d+) interface faces with elements of two or more '
                  r'other ranks', text)
    ok = m is not None and m.group(1) == remote and m.group(2) == multiple
    done(ok, f'{log}: ' + (f"'{m.group(0)}'" if m else 'no count line')
         + f' (expected {remote} and {multiple})')


def cmdExplicitRemainder(log, minLines, mode='zero'):
    pattern = re.compile(r'explicit remainder (\S+) kg/s in (\d+) cells \(of '
                         r'an explicit source sum \|b - s\| (\S+) kg/s\)')
    found = [(float(a), int(n), float(b)) for a, n, b in
             pattern.findall(open(log, errors='replace').read())]
    if not found:
        done(False, f'{log}: no explicit remainder in the transport lines')
    worstE = max(abs(e) for e, _, _ in found)
    cells = max(n for _, n, _ in found)
    sources = [b for _, _, b in found]
    summary = (f'{log}: {len(found)} transport lines: explicit remainder up '
               f'to {worstE:.3e} kg/s in up to {cells} cells, explicit '
               f'source sum |b - s| {min(sources):.3e} to {max(sources):.3e} '
               f'kg/s')
    if mode == 'nonzero':
        done(cells > 0, summary + ' (positive control: a remainder in some '
             'cells)')
    ok = len(found) >= int(minLines) and worstE == 0 and cells == 0 \
        and min(sources) > 0
    done(ok, summary + f' (required: 0 in 0 cells in at least {minLines} '
         'lines, with an explicit source > 0)')


def cmdTransportIdentity(log, minLines, limit, mode='hold'):
    pattern = re.compile(r'transport identity (\S+) kg/s \(\|D\| (\S+) of '
                         r'its terms, (\S+) kg/s\)')
    found = [(float(d), float(r), float(t)) for d, r, t in
             pattern.findall(open(log, errors='replace').read())]
    if not found:
        done(False, f'{log}: no transport identity in the transport lines')
    limit = float(limit)
    above = [f for f in found if f[1] > limit]
    worstD = max(abs(d) for d, _, _ in found)
    rels = [r for _, r, _ in found]
    terms = [t for _, _, t in found]
    summary = (f'{log}: {len(found)} transport lines: |D| up to {worstD:.2e} '
               f'kg/s, {min(rels):.2e} to {max(rels):.2e} of its terms '
               f'({min(terms):.3g} to {max(terms):.3g} kg/s); {len(above)} '
               f'lines above {limit:g}')
    if mode == 'violated':
        done(len(above) > 0, summary + ' (negative control: the identity '
             'violated in some line)')
    done(len(found) >= int(minLines) and not above,
         summary + f' (required: none above {limit:g}, in at least '
         f'{minLines} lines)')


def cmdSourceFormMutant(path):
    text = open(path).read()
    try:
        a = text.index('inline void cellDefectT (')
        b = text.index('inline void cellDefect (', a)
    except ValueError:
        done(False, f'{path}: no cellDefectT followed by cellDefect')
    body = text[a:b]
    edits = [
        ('      r[celli] = E;\n',
         '      r[celli] = cellRemainder(b[celli], ddtSource[celli], 0);\n', 1),
        ('    if (correction) {\n', '    if (false && correction) {\n', 1),
        ('if (Cp) {', 'if (false && Cp) {', 2),
    ]
    for old, new, count in edits:
        if body.count(old) != count:
            done(False, f'{path}: the statement {old.strip()!r} occurs '
                 f'{body.count(old)} times in cellDefectT, not {count}: '
                 'update sourceFormMutant to the source')
        body = body.replace(old, new)
    open(path, 'w').write(text[:a] + body + text[b:])
    done(True, f'{path}: cellDefectT edited into the source form of the '
         'face-flux correction (r_c from b_c - s_c, no face values C_f)')


def cmdMoleFraction(case, time, log, tol='1e-5'):
    pc = phaseChangeChecks()
    R, M = 8.314462618, 0.46100894
    Y, _ = pc.readScalarField(os.path.join(case, time, 'Y_PbI2_g'))
    n = len(Y)

    def field(name):
        path = os.path.join(case, time, name)
        if not os.path.exists(path):
            path = os.path.join(case, '0', name)
        values, _ = pc.readScalarField(path)
        return values*n if len(values) == 1 else values

    rho, p, T = field('rho'), field('p'), field('T')
    W = float(re.search(r'molWeight\s+([^;\s]+)\s*;', open(os.path.join(
        case, 'constant', 'thermophysicalProperties')).read()).group(1))/1000
    new = max([0.0] + [(r*y/M)/(r*y/M + pp/(R*t))
                       for r, y, pp, t in zip(rho, Y, p, T)])
    old = max([0.0] + [(y/M)/(y/M + 1/W) for y in Y])
    printed = None
    for block in open(log, errors='replace').read().split('\nTime = ')[1:]:
        if float(block.split('\n', 1)[0].strip()) == float(time):
            m = re.search(r'max mole fraction (\S+)', block)
            printed = float(m.group(1)) if m else None
    ok = printed is not None and abs(printed - new) <= float(tol)*new
    done(ok, f'{case} at {time}: max mole fraction printed {printed}, '
         f'recomputed rho Y/M against p/(R T) {new:.6g} (relative '
         f'{abs((printed or 0) - new)/max(new, 1e-300):.1e}, limit {tol}); '
         f'the old form (Y/M)/(Y/M + 1/M_carrier) gives {old:.6g}, '
         f'{new/max(old, 1e-300):.3g} times less')


def cmdInflowExact(phi, patch, rho, flow, sector, tol):
    from fractions import Fraction
    pc = phaseChangeChecks()
    _, patches = pc.readScalarField(phi)
    inflow = -sum(Fraction(v) for v in patches[patch])
    target = Fraction(rho)*Fraction(flow)/60000000*Fraction(sector)/360
    deviation = abs(float(inflow/target - 1))
    done(deviation <= float(tol),
         f'{phi}: the inflow through {patch} is rho Q {sector}/360 to '
         f'{deviation:.1e} (exact sums; limit {float(tol):g})')


COMMANDS = {
    'restartExact': cmdRestartExact,
    'compareInventories': cmdCompareInventories,
    'roundProcessorFace': cmdRoundProcessorFace,
    'samePhi': cmdSamePhi,
    'heldWall': cmdHeldWall,
    'sameWallColumns': cmdSameWallColumns,
    'conductivityTable': cmdConductivityTable,
    'continuity': cmdContinuity,
    'continuityAbove': cmdContinuityAbove,
    'carrier': cmdCarrier,
    'layerCounts': cmdLayerCounts,
    'mdep': cmdMdep,
    'facesAgree': cmdFacesAgree,
    'plateau': cmdPlateau,
    'solidAgree': cmdSolidAgree,
    'tdep': cmdTdep,
    'tdepDecided': cmdTdepDecided,
    'boundWidth': cmdBoundWidth,
    'positions': cmdPositions,
    'exactLayer': cmdExactLayer,
    'perturbWedge': cmdPerturbWedge,
    'faceCentreRule': cmdFaceCentreRule,
    'wallFlux': cmdWallFlux,
    'mdepPatches': cmdMdepPatches,
    'samePatches': cmdSamePatches,
    'patchValue': cmdPatchValue,
    'closureRuns': cmdClosureRuns,
    'sameBalance': cmdSameBalance,
    'mdepCount': cmdMdepCount,
    'explicitRemainder': cmdExplicitRemainder,
    'transportIdentity': cmdTransportIdentity,
    'sourceFormMutant': cmdSourceFormMutant,
    'moleFraction': cmdMoleFraction,
    'inflowExact': cmdInflowExact,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print(__doc__)
        sys.exit(2)
    COMMANDS[sys.argv[1]](*sys.argv[2:])
