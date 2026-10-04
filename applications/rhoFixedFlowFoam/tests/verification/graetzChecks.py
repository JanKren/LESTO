#!/usr/bin/env python3
"""
tests/verification/graetzChecks.py

PURPOSE

Numerical checks of the Graetz-Robin verification (tests/verification/
graetzRobin; doc/phase-change-plan.md, section 6, M7 acceptance 1).  Every
subcommand except fit and table prints one summary line (plus detail lines)
and exits with 0 (criterion met) or 1.

  reference <out.json> <Pe> <L/R> <tol> <Bi> [<Bi> ...]
      the reference decay rate m = mu R of every Bi by three independent
      methods (graetz.py: power series in mpmath, DOP853 shooting, Newton-
      refined Chebyshev collocation N = 32 and 48), its error bound (the
      largest relative difference of the methods) and the fit window;
      met when every bound <= tol
  mesh <tol> <mesh case> [<mesh case> ...]
      every level: checkMesh 'Mesh OK', non-orthogonality 0, both wedge
      patches at theta/2, nCells = Nx Nr; consecutive levels: Nr and Nx
      double (refinement ratio 2 in r and x); the carrier
      (log.graetzCarrier): max |sum_f phi| per cell <= tol x the inflow, and
      the inlet mean velocity equals U_b to 1e-12
  setup <case> [<case> ...]
      the log of every run: D_PbI2_g uniform and equal to the correlation
      at (T, p) to 1e-15; v_d = k_w; every element below Tdep (f = 1);
      halfCell: the applied rate sum_e A_e v'_e/V_I equals
      v' A_I/V_I with v' = 1/(1/k_w + d/D), d = R - y_last, to 1e-12 with
      y_last the centroid of the wall row on the unit radius, times R
      (graetz.half_cell; the agreement with OpenFOAM's centroid of the wall
      row, 0/C of the mesh case, is printed too)
  fit <ref.json> <case> [<case> ...]
      the decay rate of every run (written to <case>/graetzFit.json and
      <case>/graetzLocalRate.dat; always exit 0):
        mu_fit   least-squares slope of ln n'_g (the gas line density, the
                 gas column of postProcessing/axialProfiles at the last
                 write, one bin per axial cell) over the bins whose centres
                 lie in the window of the reference record
        variants the same over the window with its ends moved by -1.5 R,
                 0 or +1.5 R each (the 8 combinations but the window
                 itself: one end moved, or both; the outlet end at most
                 L - 1 R)
        mixing   the mixing-cup concentration sum U V Y/sum U V per axial
                 slab (Y_PbI2_g of the last write, U and the cell volumes of
                 the mesh case), same window
        wall     the wall uptake per unit length (the increment of the wall
                 column over the last step), same window
        mu_h     the discrete predictor (graetz.discrete_mode) of the run's
                 mesh, wall law and scheme
        steady   max over the window bins of |n'_g(last)/n'_g(previous) - 1|
  steady <tol> <case> [<case> ...]
      every run: steady (above) <= tol; the outer loop of the last step
      converged (initial residual < its tolerance); no guard pass and no
      significant projection or complementarity violation in the run
      (end-of-run summary); |relative closure| <= 1e-13 in every Balance
      line
  discrete <tol> <case> [<case> ...]
      every run: |mu_fit/mu_h - 1| <= tol (the solver solves the documented
      discretisation, wall law included)
  windows <tol> <case> [<case> ...]
      every run: the window variants (one end moved, and both), the
      mixing-cup and the wall-uptake fits within tol of mu_fit; the
      largest change of each kind is printed
  outlet <tol> <zeroGradient case> <fixedValue case> [<zg> <fv> ...]
      every pair (same mesh, Bi and law): |mu_fit(fv)/mu_fit(zg) - 1| <= tol

  discrete, windows and outlet fail on a deviation that is not finite (a
  fit of fewer than 3 bins or of a non-positive profile is NaN): their
  largest deviation is a max(), and max(0, nan) is 0, so a NaN fit passed
  (plan section 35, item O); steady fails on a closure that is not finite
  (item X).  The runs with one are named ('not finite: <run>').
  gci <ref.json> <gci|extrapolation|order> <value> <coarse> <medium> <fine>
      Celik et al. (2008) on mu_fit of three runs with refinement ratio 2:
        gci            GCI_fine <= value (e.g. 0.005)
        extrapolation  |mu_ext - mu_ref| <= GCI_fine mu_fine + bound mu_ref
                       (value unused)
        order          |p - value| <= 0.25 (value = the formal order) and
                       monotone convergence
        all            gci (value), extrapolation and order (the formal
                       order of the run's wall law: 2 halfCell, 1 none)
  table <ref.json> <out.md> <case> [<case> ...]
      a markdown table of every run and the GCI of every complete triplet
      (by scheme, outlet, Bi and law), for the readme and the plan

------------------------------------------------------------------------------
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graetz  # noqa: E402

GAS = 'PbI2_g'


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

def params(case):
    out = {}
    for line in open(os.path.join(case, 'graetzParameters')):
        k, v = line.split(None, 1)
        v = v.strip()
        try:
            out[k] = int(v)
        except ValueError:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def times(case):
    base = os.path.join(case, 'postProcessing', 'axialProfiles')
    ts = [os.path.basename(d) for d in glob.glob(os.path.join(base, '*'))]
    return sorted(ts, key=float)


def profile(case, t):
    path = os.path.join(case, 'postProcessing', 'axialProfiles', t,
                        GAS + '_native.dat')
    return np.loadtxt(path)


def readField(path):
    """Internal field of an ascii OpenFOAM vol field: scalars as (n,),
    vectors as (n, 3); 'uniform' as a single value."""
    text = open(path).read()
    m = re.search(r'internalField\s+uniform\s+([^;]+);', text)
    if m:
        v = m.group(1).strip()
        if v.startswith('('):
            return np.array([float(x) for x in v.strip('()').split()])
        return float(v)
    m = re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(',
                  text)
    if not m:
        raise ValueError("no internalField in %s" % path)
    kind, n = m.group(1), int(m.group(2))
    body = text[m.end():]
    if kind == 'scalar':
        vals = np.array(body.split(None, n)[:n], dtype=float)
        return vals
    if kind == 'vector':
        nums = re.findall(r'\(([^()]*)\)', body[:body.find('\n)')+2])[:n]
        return np.array([[float(x) for x in s.split()] for s in nums])
    raise ValueError(kind)


def meshCase(case):
    poly = os.path.join(case, 'constant', 'polyMesh')
    if os.path.islink(poly):
        return os.path.dirname(os.path.dirname(os.path.realpath(poly)))
    return case


def reference(refFile, Bi):
    for r in json.load(open(refFile)):
        if abs(r['Bi'] - Bi) <= 1e-12*abs(Bi):
            return r
    raise SystemExit("no reference record for Bi = %r in %s" % (Bi, refFile))


def logText(case):
    return open(os.path.join(case, 'log')).read()


def label(case):
    """.../run/[SuperBee/]<outlet>/bi<Bi>/<law>/nr<Nr> -> the part after run/;
    .../transient/... and .../mesh/... keep their first component"""
    parts = os.path.normpath(case).split(os.sep)
    for k in range(len(parts) - 1, -1, -1):
        if parts[k] == 'run':
            return '/'.join(parts[k + 1:])
        if parts[k] in ('transient', 'mesh'):
            return '/'.join(parts[k:])
    return '/'.join(parts[-4:])


# ----------------------------------------------------------------------------
# Fits
# ----------------------------------------------------------------------------

def slope(x, y):
    xm = x.mean()
    return np.sum((x - xm)*(y - y.mean()))/np.sum((x - xm)**2)


def decay(x, values, xs, xe):
    sel = (x >= xs) & (x <= xe)
    if sel.sum() < 3 or np.any(values[sel] <= 0):
        return float('nan'), int(sel.sum())
    return -slope(x[sel], np.log(values[sel])), int(sel.sum())


def mixingCup(case, p, t):
    mc = meshCase(case)
    try:
        C = readField(os.path.join(mc, '0', 'C'))
        V = readField(os.path.join(mc, '0', 'V'))
        U = readField(os.path.join(mc, '0', 'U'))
        Y = readField(os.path.join(case, t, 'Y_' + GAS))
    except (OSError, ValueError):
        return None, None
    j = np.floor(C[:, 0]/p['dx']).astype(int)
    w = U[:, 0]*V
    num = np.bincount(j, weights=w*Y, minlength=p['nx'])
    den = np.bincount(j, weights=w, minlength=p['nx'])
    xc = (np.arange(p['nx']) + 0.5)*p['dx']
    return xc, num/den


def fitCase(ref, case):
    p = params(case)
    R = p['R']
    ts = times(case)
    last, prev = ts[-1], ts[-2]
    d = profile(case, last)
    dp = profile(case, prev)
    x, gas, wall = d[:, 2], d[:, 3], d[:, 4]
    xs, xe = ref['window']
    LR = p['L_over_R']
    mu, n = decay(x, gas, xs*R, xe*R)
    variants = {}
    for da in (-1.5, 0.0, 1.5):
        for db in (-1.5, 0.0, 1.5):
            if da == 0 and db == 0:
                continue
            name = 'start%+.1fR,end%+.1fR' % (da, db)
            variants[name] = decay(x, gas, (xs + da)*R,
                                   min(xe + db, LR - 1)*R)[0]
    muWall = decay(x, wall - dp[:, 4], xs*R, xe*R)[0]
    xc, cup = mixingCup(case, p, last)
    muCup = decay(xc, cup, xs*R, xe*R)[0] if cup is not None else None
    sel = (x >= xs*R) & (x <= xe*R)
    steady = float(np.max(np.abs(gas[sel]/dp[sel, 3] - 1)))
    steadyAll = float(np.max(np.abs(gas/dp[:, 3] - 1)))
    history = []
    for t0, t1 in zip(ts[:-1], ts[1:]):
        g0, g1 = profile(case, t0)[:, 3], profile(case, t1)[:, 3]
        history.append(float(np.max(np.abs(g1[sel]/g0[sel] - 1))))
    muH = graetz.discrete_mode(p['Pe'], p['Bi'], p['nr'],
                               p['aspect']/p['nr'], p['law'],
                               p['scheme'])/R
    local = -np.diff(np.log(gas))/np.diff(x)
    with open(os.path.join(case, 'graetzLocalRate.dat'), 'w') as f:
        f.write('# x/R  local decay rate -dln(n_g)/dx [1/m]  relative to '
                'mu_h (predictor)  n_g/n_g(first bin)\n')
        for k in range(len(local)):
            xm = 0.5*(x[k] + x[k + 1])
            f.write('%.10g %.17g %.6e %.6e\n' % (xm/R, local[k],
                                                  local[k]/muH - 1,
                                                  gas[k]/gas[0]))
    out = dict(case=case, law=p['law'], Bi=p['Bi'], nr=p['nr'], nx=p['nx'],
               scheme=p['scheme'], outlet=p['outlet'], dt=p['dt'], Co=p['Co'],
               time=last, previous=prev, window=[xs, xe], nBins=n,
               mu=mu, m=mu*R, variants=variants, muWall=muWall, muCup=muCup,
               muH=muH, steady=steady, steadyAll=steadyAll,
               history=history, times=ts,
               muRef=ref['m']/R, bound=ref['bound'])
    json.dump(out, open(os.path.join(case, 'graetzFit.json'), 'w'), indent=1)
    return out


def fitOf(case):
    return json.load(open(os.path.join(case, 'graetzFit.json')))


# ----------------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------------

def cmdReference(args):
    out, Pe, LR, tol = args[0], float(args[1]), float(args[2]), float(args[3])
    recs = [graetz.reference_record(Pe, float(b), LR) for b in args[4:]]
    json.dump(recs, open(out, 'w'), indent=1)
    ok = True
    for r in recs:
        ok &= r['bound'] <= tol
        print("  Bi = %g: m = mu R = %s (series), %s; error bound %.2e; fit "
              "window [%.3f, %.3f] R" % (r['Bi'], r['m_text'], ", ".join(
                  "%s %.17g" % kv for kv in r['methods'].items()
                  if kv[0] != 'series'), r['bound'], *r['window']))
    print("reference: %d Biot numbers, largest error bound %.2e (limit %g)"
          % (len(recs), max(r['bound'] for r in recs), tol))
    return ok


def cmdMesh(args):
    tol = float(args[0])
    cases = args[1:]
    ok = True
    prev = None
    for c in cases:
        p = params(c)
        cm = open(os.path.join(c, 'log.checkMesh')).read()
        nCells = int(re.search(r'cells:\s+(\d+)', cm).group(1))
        nonOrth = float(re.search(r'Mesh non-orthogonality Max: (\S+)',
                                  cm).group(1))
        wedges = [float(a) for a in re.findall(r'Wedge \w+ with angle (\S+)'
                                               r' degrees', cm)]
        meshOk = 'Mesh OK.' in cm
        lc = open(os.path.join(c, 'log.graetzCarrier')).read()
        net = float(re.search(r'max \|sum_f phi\| per cell = (\S+)',
                              lc).group(1))
        inflow = -float(re.search(r'patch INLET: sum phi = (\S+)',
                                  lc).group(1))
        umean = float(re.search(r'patch INLET: .*mean velocity (\S+)',
                                lc).group(1))
        good = (meshOk and nonOrth == 0 and nCells == p['nx']*p['nr']
                and len(wedges) == 2
                and all(abs(w - p['theta']/2) < 1e-6 for w in wedges)
                and net <= tol*inflow
                and abs(umean/p['Ub'] - 1) <= 1e-12)
        if prev is not None:
            good &= p['nr'] == 2*prev['nr'] and p['nx'] == 2*prev['nx']
        ok &= good
        print("  %s: %d x %d = %d cells, %s, non-orthogonality %g, wedges %s"
              " deg; carrier max|sum phi| %.2e kg/s = %.2e of the inflow, "
              "inlet mean velocity/U_b - 1 = %.1e%s"
              % (label(c), p['nx'], p['nr'], nCells,
                 'Mesh OK' if meshOk else 'checkMesh FAILED', nonOrth,
                 wedges, net, net/inflow, umean/p['Ub'] - 1,
                 '' if good else '  <-- FAIL'))
        prev = p
    print("meshes: %d levels, ratio 2 in r and x, carrier continuity <= %g "
          "of the inflow" % (len(cases), tol))
    return ok


_wallRowCentre = {}


def wallRowCentre(case):
    """y of the cell centres of the wall row as OpenFOAM computes them (the
    largest y of 0/C of the run's mesh case), cached per mesh"""
    mc = meshCase(case)
    if mc not in _wallRowCentre:
        C = readField(os.path.join(mc, '0', 'C'))
        _wallRowCentre[mc] = float(np.max(C[:, 1]))
    return _wallRowCentre[mc]


def cmdSetup(args):
    ok = True
    worst = dict(exact=0.0, openfoam=0.0)
    for c in args:
        p = params(c)
        log = logText(c)
        dmin, dmax = (float(v) for v in re.search(
            r'D_%s internal min/max \[m2/s\]: (\S+) (\S+)' % GAS,
            log).groups())
        vd = float(re.search(r'Temperature model: v_d = (\S+) m/s',
                             log).group(1))
        below = re.search(r'(\d+) of (\d+) elements below Tdep', log)
        good = (dmin == dmax and abs(dmin/p['D'] - 1) <= 1e-15
                and vd == p['kw'] and below.group(1) == below.group(2))
        extra = ''
        if p['law'] == 'halfCell':
            # the half-cell distance d = R - y_last of the wall row: from the
            # centroid (2/3)(r2^3 - r1^3)/(r2^2 - r1^2) on the unit radius
            # times R (graetz.half_cell, the gate) and from OpenFOAM's
            # centroid of the wall row (both agree to a few 1e-14, as the
            # centroid in exact arithmetic; the formula in metres loses
            # digits in r2^3 - r1^3)
            applied = float(re.search(r"applied rate sum_e A_e v'_e/V_I = "
                                      r"(\S+) 1/s", log).group(1))
            aiv = float(re.search(r'A_I/V_I = (\S+) 1/m', log).group(1))
            rels = {}
            for kind, d in (('exact', graetz.half_cell(p['nr'])*p['R']),
                            ('openfoam', p['R'] - wallRowCentre(c))):
                vprime = 1.0/(1.0/p['kw'] + d/p['D'])
                rels[kind] = applied/(vprime*aiv) - 1
                worst[kind] = max(worst[kind], abs(rels[kind]))
            good &= abs(rels['exact']) <= 1e-12
            extra = (", halfCell applied rate/(v' A_I/V_I) - 1 = %.1e (d "
                     "from the centroid of the unit radius), %.1e (from "
                     "OpenFOAM's)"
                     % (rels['exact'], rels['openfoam']))
        ok &= good
        print("  %s: D = %r (min = max: %s), v_d = %r m/s, %s/%s elements "
              "below Tdep%s%s" % (label(c), dmin, dmin == dmax, vd,
                                    below.group(1), below.group(2), extra,
                                    '' if good else '  <-- FAIL'))
    print("set-up: %d runs, uniform D = the PbI2-He correlation at (T, p), "
          "f = 1, the half-cell law of the mesh: the applied rate equals "
          "v'(d) A_I/V_I, d = R - y_last of the wall row, to %.1e with the "
          "centroid of the unit radius (limit 1e-12) and %.1e with "
          "OpenFOAM's"
          % (len(args), worst['exact'], worst['openfoam']))
    return ok


def cmdFit(args):
    refFile = args[0]
    ok = True
    for c in args[1:]:
        try:
            p = params(c)
            f = fitCase(reference(refFile, p['Bi']), c)
        except Exception as e:      # report and go on with the other runs
            ok = False
            print("  %s: no fit (%s: %s)" % (label(c), type(e).__name__, e))
            continue
        print("  %s: mu_fit = %.15g 1/m (%d bins), mu_h = %.15g, "
              "fit/mu_h - 1 = %+.2e, steady %.1e"
              % (label(c), f['mu'], f['nBins'], f['muH'],
                 f['mu']/f['muH'] - 1, f['steady']))
    return ok


def cmdSteady(args):
    tol = float(args[0])
    ok, bad = True, []
    for c in args[1:]:
        f = fitOf(c)
        log = logText(c)
        lastEx = re.findall(r'exchange: .*outer iterations (\d+), (\w+) '
                            r'\(initial residual (\S+), tolerance (\S+);',
                            log)[-1]
        summary = re.search(r'Phase change summary Y_%s: (\d+) steps without'
                            r' outer convergence, (\d+) guard passes, (\d+) '
                            r'projected elements \((\d+) significant\), (\d+)'
                            r' complementarity violations \((\d+) '
                            r'significant\)' % GAS, log)
        closures = [abs(float(v)) for v in re.findall(
            r'^Balance %s .* closure \S+ \(relative (\S+)\)' % GAS, log,
            re.M)]
        # a closure that is not finite fails the run: max() passes over a
        # NaN that is not the first value (plan section 35, item X)
        finite = all(math.isfinite(x) for x in closures)
        if not finite:
            bad.append(label(c))
        worstClosure = (max(closures) if closures and finite
                        else float('nan'))
        good = (f['steady'] <= tol and lastEx[1] == 'converged'
                and summary is not None
                and all(int(summary.group(k)) == 0 for k in (2, 4, 6))
                and len(closures) > 0 and finite
                and worstClosure <= 1e-13)
        ok &= good
        print("  %s: profile change %.1e (window), %.1e (all bins) between "
              "%s and %s (history %s); last step %s (initial residual %s); "
              "summary: %s steps without convergence, %s guard passes, %s "
              "projected (%s significant), %s complementarity (%s "
              "significant); max |closure| %.1e%s"
              % (label(c), f['steady'], f['steadyAll'], f['previous'],
                 f['time'], ' '.join('%.1e' % h for h in f['history']),
                 lastEx[1], lastEx[2],
                 *(summary.group(k) for k in range(1, 7)), worstClosure,
                 '' if good else '  <-- FAIL'))
    print("steady state: %d runs, profile change <= %g between the last two "
          "writes%s" % (len(args) - 1, tol, notFinite(bad)))
    return ok


def notFinite(bad):
    """the text naming the runs with a deviation that is not finite"""
    return ('; not finite: %s' % ', '.join(bad[:4])
            + (' and %d more' % (len(bad) - 4) if len(bad) > 4 else '')
            if bad else '')


def cmdDiscrete(args):
    tol = float(args[0])
    worst, bad = 0, []
    for c in args[1:]:
        f = fitOf(c)
        dev = f['mu']/f['muH'] - 1
        if math.isfinite(dev):
            worst = max(worst, abs(dev))
        else:
            bad.append(label(c))
        print("  %s: mu_fit %.15g, mu_h %.15g, relative %+.2e"
              % (label(c), f['mu'], f['muH'], dev))
    print("discrete: %d runs, max |mu_fit/mu_h - 1| = %.2e (limit %g)%s"
          % (len(args) - 1, worst, tol, notFinite(bad)))
    return not bad and worst <= tol


def windowKind(name):
    """'one' when a variant moves one end of the window, 'both' when it
    moves both"""
    return 'both' if '+0.0R' not in name else 'one'


def cmdWindows(args):
    tol = float(args[0])
    worst = dict(one=0.0, both=0.0, mixing=0.0, wall=0.0)
    bad = []
    for c in args[1:]:
        f = fitOf(c)
        devs = {k: v/f['mu'] - 1 for k, v in f['variants'].items()}
        kinds = {k: windowKind(k) for k in devs}
        devs['wall'] = f['muWall']/f['mu'] - 1
        kinds['wall'] = 'wall'
        if f['muCup'] is not None:
            devs['mixing'] = f['muCup']/f['mu'] - 1
            kinds['mixing'] = 'mixing'
        for k, v in devs.items():
            if math.isfinite(v):
                worst[kinds[k]] = max(worst[kinds[k]], abs(v))
            else:
                bad.append('%s %s' % (label(c), k))
        print("  %s: %s" % (label(c), ", ".join(
            "%s %+.1e" % kv for kv in devs.items())))
    print("windows: %d runs, largest relative change of mu: one window end "
          "moved by -+1.5 R %.2e, both ends %.2e, mixing cup %.2e, wall "
          "uptake %.2e (limit %g)%s"
          % (len(args) - 1, worst['one'], worst['both'], worst['mixing'],
             worst['wall'], tol, notFinite(bad)))
    return not bad and max(worst.values()) <= tol


def cmdOutlet(args):
    tol = float(args[0])
    pairs = args[1:]
    worst, bad = 0, []
    for a, b in zip(pairs[0::2], pairs[1::2]):
        fa, fb = fitOf(a), fitOf(b)
        dev = fb['mu']/fa['mu'] - 1
        if math.isfinite(dev):
            worst = max(worst, abs(dev))
        else:
            bad.append('%s vs %s' % (label(a), label(b)))
        print("  %s vs %s: %+.2e" % (label(a), label(b), dev))
    print("outlet: %d pairs, max |mu(fixedValue 0)/mu(zeroGradient) - 1| = "
          "%.2e (limit %g)%s" % (len(pairs)//2, worst, tol, notFinite(bad)))
    return not bad and worst <= tol


def cmdGci(args):
    refFile, what, value = args[0], args[1], float(args[2])
    fs = [fitOf(c) for c in args[3:6]]
    ref = reference(refFile, fs[0]['Bi'])
    g = graetz.gci(fs[0]['mu'], fs[1]['mu'], fs[2]['mu'])
    label = "%s Bi %g, %s, Nr %d/%d/%d" % (fs[0]['law'], fs[0]['Bi'],
                                           fs[0]['scheme'], fs[0]['nr'],
                                           fs[1]['nr'], fs[2]['nr'])
    muRef = fs[0]['muRef']
    errExt = g['f_ext'] - muRef
    allowed = g['gci']*fs[2]['mu'] + ref['bound']*muRef
    print("  mu = %.12g, %.12g, %.12g 1/m; p = %.4f; mu_ext = %.12g; "
          "GCI_fine = %.3g %%; mu_ref = %.12g (bound %.1e); mu_ext - mu_ref "
          "= %.3e (%.2e relative), allowed %.3e; errors of the fine run "
          "%.2e"
          % (fs[0]['mu'], fs[1]['mu'], fs[2]['mu'], g['p'], g['f_ext'],
             100*g['gci'], muRef, ref['bound'], errExt, errExt/muRef,
             allowed, fs[2]['mu']/muRef - 1))
    formal = 2.0 if fs[0]['law'] == 'halfCell' else 1.0
    if what == 'gci':
        print("%s: GCI_fine = %.3g %% (limit %g %%)"
              % (label, 100*g['gci'], 100*value))
        return g['gci'] <= value
    if what == 'all':
        ok = (g['gci'] <= value and abs(errExt) <= allowed
              and abs(g['p'] - formal) <= 0.25 and g['monotone'])
        print("%s: GCI_fine = %.3g %% (limit %g %%), p = %.4f (formal %g), "
              "|mu_ext - mu_ref| = %.2e <= %.2e: %s"
              % (label, 100*g['gci'], 100*value, g['p'], formal, abs(errExt),
                 allowed, ok))
        return ok
    if what == 'extrapolation':
        print("%s: |mu_ext - mu_ref| = %.2e <= GCI_fine mu_fine + bound "
              "mu_ref = %.2e" % (label, abs(errExt), allowed))
        return abs(errExt) <= allowed
    if what == 'order':
        print("%s: observed order p = %.4f (formal %g, accepted within 0.25;"
              " monotone %s)" % (label, g['p'], value, g['monotone']))
        return abs(g['p'] - value) <= 0.25 and g['monotone']
    raise SystemExit("unknown gci check " + what)


def cmdTable(args):
    refFile, out = args[0], args[1]
    fs = [fitOf(c) for c in args[2:]]
    lines = ['| scheme | outlet | Bi | law | Nr | Nx | mu_fit [1/m] | '
             'mu_fit/mu_ref - 1 | mu_fit/mu_h - 1 | window spread | '
             'mixing cup | steady | Co |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for f in fs:
        spread = max(abs(v/f['mu'] - 1) for v in f['variants'].values())
        cup = ('%+.1e' % (f['muCup']/f['mu'] - 1)) if f['muCup'] else '-'
        lines.append('| %s | %s | %g | %s | %d | %d | %.12g | %+.3e | %+.1e |'
                     ' %.1e | %s | %.0e | %.1e |'
                     % (f['scheme'], f['outlet'], f['Bi'], f['law'], f['nr'],
                        f['nx'], f['mu'], f['mu']/f['muRef'] - 1,
                        f['mu']/f['muH'] - 1, spread, cup, f['steady'],
                        f['Co']))
    lines += ['', '| scheme | outlet | Bi | law | Nr | p | mu_ext [1/m] | '
              'mu_ext/mu_ref - 1 | GCI_fine | allowed/mu_ref |',
              '|---|---|---|---|---|---|---|---|---|---|']
    groups = {}
    for f in fs:
        groups.setdefault((f['scheme'], f['outlet'], f['Bi'], f['law']),
                          {})[f['nr']] = f
    for key in sorted(groups):
        g = groups[key]
        nrs = sorted(g)
        for k in range(len(nrs) - 2):
            a, b, c = (g[n] for n in nrs[k:k + 3])
            if not (b['nr'] == 2*a['nr'] and c['nr'] == 2*b['nr']):
                continue
            r = graetz.gci(a['mu'], b['mu'], c['mu'])
            ref = reference(refFile, a['Bi'])
            muRef = a['muRef']
            allowed = (r['gci']*c['mu'] + ref['bound']*muRef)/muRef
            lines.append('| %s | %s | %g | %s | %d/%d/%d | %.4f | %.12g | '
                         '%+.2e | %.3g %% | %.2e |'
                         % (key[0], key[1], key[2], key[3], a['nr'], b['nr'],
                            c['nr'], r['p'], r['f_ext'],
                            r['f_ext']/muRef - 1, 100*r['gci'], allowed))
    open(out, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return True


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    cmds = dict(reference=cmdReference, mesh=cmdMesh, setup=cmdSetup,
                fit=cmdFit, steady=cmdSteady, discrete=cmdDiscrete,
                windows=cmdWindows, outlet=cmdOutlet, gci=cmdGci,
                table=cmdTable)
    if argv[1] not in cmds:
        print(__doc__)
        return 2
    return 0 if cmds[argv[1]](argv[2:]) else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
