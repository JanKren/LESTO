#!/usr/bin/env python3
"""
tests/verification/graetzCase.py

PURPOSE

Fills the case template graetzWedge/ of the M7 Graetz-Robin verification
(tests/verification/graetzRobin) for one mesh or one run.  The physical
parameters are fixed for the study (doc/phase-change-plan.md, section 6, M7
acceptance 1; plan section 33 gives the reasons):

  R  = 2.4e-3 m          pipe radius; the WALL faces lie at y = R
  T  = 500 K, p = 101325 Pa, helium (perfect gas): rho = p M/(R_gas T)
  D  = D_PbI2He(T, p)    uniform: PbI2HeDiffusivity.H at constant T and p
                         (createVariableDiffusivity.H), 7.6963204441542405e-05
                         m2/s
  Pe = U_b R/D = 5       U_b = Pe D/R
  Bi = k_w R/D           k_w = Bi D/R = depositionVelocity, with k = 1 1/K
                         and Tdep = T + 100 K: f = 1 - exp(-100) = 1 exactly
  L  = 24 R              wedge angle 2 degrees; uniform cells, dx = 2 dr,
                         Nx = 12 Nr

Usage
  graetzCase.py <template> <case> --nr N --bi BI --law none|halfCell
                [--outlet zeroGradient|fixedValue] [--scheme linear|SuperBee]
                [--dt 1e6] [--steps 4] [--write-interval N]
                [--mesh <mesh case>]

  --steps           end time = steps x dt
  --write-interval  steps between field writes; default steps (the fields
                    are written at the end only and the profiles after
                    every step, as profile times); with any other value the
                    profiles are written at the write times only

  --mesh  link constant/polyMesh and the frozen carrier 0/{T,p,U,rho,phi}
          of an existing mesh case (made by this script, blockMesh and
          graetzCarrier) instead of keeping the template's own

The derived parameters are written to <case>/graetzParameters (key value),
which graetzChecks.py reads.
------------------------------------------------------------------------------
"""
import argparse
import math
import os
import shutil

R = 2.4e-3
T = 500.0
P = 101325.0
PE = 5.0
THETA_DEG = 2.0
L_OVER_R = 24
ASPECT = 2           # dx/dr
L_TEXT = "0.0576"    # 24 R: the same text in blockMeshDict and the bins
M_CARRIER = 4.0026e-3
R_GAS = 8.314462618
CARRIER = ('T', 'p', 'U', 'rho', 'phi')


def diffusivity(T, p, pressureScale=101325.0):
    """PbI2HeDiffusivity.H, operation by operation (the same double as the
    C++ header for these inputs)."""
    sigma = 4.2584682449728355
    epsilonByK = 71.11785702515257
    reducedT = T/epsilonByK
    omega = (0.5313 - 1.2946*math.exp(-0.9922*reducedT) + 1.5591/reducedT
             - 0.1918/(reducedT*reducedT))
    return (1.86e-7*math.pow(T, 1.5)*math.sqrt(1.0/4.0026 + 1.0/461.0)
            / ((p/pressureScale)*sigma*sigma*omega))


def derived(nr, bi):
    D = diffusivity(T, P)
    L = float(L_TEXT)
    nx = L_OVER_R*nr//ASPECT
    return dict(R=R, T=T, p=P, D=D, Pe=PE, Bi=bi, Ub=PE*D/R, kw=bi*D/R,
                rho=P*M_CARRIER/(R_GAS*T), L=L, L_over_R=L_OVER_R, nr=nr,
                nx=nx, dr=R/nr, dx=L/nx, aspect=ASPECT, theta=THETA_DEG,
                t=math.tan(math.radians(THETA_DEG/2)), umax=2*PE*D/R)


def solverEntry():
    return ("solver GAMG; smoother GaussSeidel; directSolveCoarsest yes;\n"
            "    tolerance 1e-13; relTol 0; maxIter 2000;")


def fill(template, case, values):
    if os.path.exists(case):
        shutil.rmtree(case)
    shutil.copytree(template, case)
    for root, _, files in os.walk(case):
        for f in files:
            path = os.path.join(root, f)
            text = open(path).read()
            for k, v in values.items():
                text = text.replace('@' + k + '@', v)
            left = [w for w in text.split() if w.count('@') >= 2]
            if left:
                raise SystemExit("unfilled placeholder in %s: %s"
                                 % (path, left[:3]))
            open(path, 'w').write(text)


def linkMesh(case, mesh):
    mesh = os.path.abspath(mesh)
    poly = os.path.join(case, 'constant', 'polyMesh')
    if os.path.exists(poly):
        shutil.rmtree(poly)
    os.symlink(os.path.join(mesh, 'constant', 'polyMesh'), poly)
    for f in CARRIER:
        dst = os.path.join(case, '0', f)
        if os.path.lexists(dst):
            os.remove(dst)
        src = os.path.join(mesh, '0', f)
        if not os.path.exists(src):
            raise SystemExit("no carrier field %s" % src)
        os.symlink(src, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('template')
    ap.add_argument('case')
    ap.add_argument('--nr', type=int, required=True)
    ap.add_argument('--bi', type=float, required=True)
    ap.add_argument('--law', choices=('none', 'halfCell'), required=True)
    ap.add_argument('--outlet', choices=('zeroGradient', 'fixedValue'),
                    default='zeroGradient')
    ap.add_argument('--scheme', choices=('linear', 'SuperBee'),
                    default='linear')
    ap.add_argument('--dt', type=float, default=1e6)
    ap.add_argument('--steps', type=int, default=4)
    ap.add_argument('--write-interval', type=int, default=None)
    ap.add_argument('--mesh', default=None)
    a = ap.parse_args()

    d = derived(a.nr, a.bi)
    co = d['umax']*a.dt/d['dx']
    summary = ("Pe = %r, Bi = %r, D = %r m2/s, U_b = %r m/s, k_w = %r m/s,"
               " Nr = %d, Nx = %d, law %s, outlet %s, scheme %s, dt = %r s"
               " (Co_max = %.3g)" % (d['Pe'], d['Bi'], d['D'], d['Ub'],
                                     d['kw'], d['nr'], d['nx'], a.law,
                                     a.outlet, a.scheme, a.dt, co))
    outlet = ('type zeroGradient;' if a.outlet == 'zeroGradient'
              else 'type fixedValue; value uniform 0;')
    nCorr = '1' if a.scheme == 'linear' else '20'
    interval = a.write_interval or a.steps
    ptimes = (' '.join(repr(k*a.dt) for k in range(1, a.steps + 1))
              if interval == a.steps else '')
    values = {
        'L': L_TEXT, 'R': repr(d['R']), 'RT': repr(d['R']*d['t']),
        'MRT': repr(-d['R']*d['t']), 'NX': str(d['nx']), 'NR': str(d['nr']),
        'UB': repr(d['Ub']), 'T': repr(d['T']), 'P': repr(d['p']),
        'ENDTIME': repr(a.dt*a.steps), 'DT': repr(a.dt),
        'WRITEINTERVAL': str(interval), 'PTIMES': ptimes,
        'OUTERTOL': '1e-12', 'NCORR': nCorr, 'SCHEME': a.scheme,
        'SOLVER': solverEntry(),
        'SOLVERNOTE': ('GAMG to 1e-13 (relTol 0, at most 2000 iterations: '
                       'the first solve on 786432 cells stops there at '
                       '1.4e-13); outer loop nCorr %s, tolerance 1e-12'
                       % nCorr),
        'LAW': a.law, 'KW': repr(d['kw']), 'TDEP': repr(d['T'] + 100.0),
        'OUTLETBC': outlet, 'SUMMARY': summary,
    }
    fill(a.template, a.case, values)
    if a.mesh:
        linkMesh(a.case, a.mesh)
    with open(os.path.join(a.case, 'graetzParameters'), 'w') as f:
        for k in ('R', 'T', 'p', 'D', 'Pe', 'Bi', 'Ub', 'kw', 'rho', 'L',
                  'L_over_R', 'nr', 'nx', 'dr', 'dx', 'aspect', 'theta', 't',
                  'umax'):
            f.write('%s %r\n' % (k, d[k]))
        f.write('law %s\noutlet %s\nscheme %s\ndt %r\nsteps %d\n'
                'writeInterval %d\nCo %r\n'
                % (a.law, a.outlet, a.scheme, a.dt, a.steps, interval, co))
    open(os.path.join(a.case, 'case.foam'), 'w').close()
    print(summary)


if __name__ == '__main__':
    main()
