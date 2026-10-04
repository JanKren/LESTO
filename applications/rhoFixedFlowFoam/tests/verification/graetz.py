#!/usr/bin/env python3
"""
tests/verification/graetz.py

PURPOSE

Reference of the M7 Graetz-Robin verification (doc/phase-change-plan.md,
section 6, M7 acceptance 1): the extended Graetz problem (fully developed
laminar pipe flow, axial diffusion included) with a first-order wall
uptake, its discrete counterpart on the wedge meshes of graetzRobin, and
the grid convergence index.

Continuous problem (derived from the solver's equation and wall law).
For a paired gas species rhoFixedFlowFoam solves (plan E1)

    d(rho Y)/dt + div(phi Y) - div(rho D grad Y) = -(1/V_c) sum_e g_e Y_c,

and the temperature law gives every wall element (face) e the uptake
g_e = A_e rho v with v = v_d f(T_e), f = 1 - exp(-k (Tdep - T_e)) (plan E2,
readme "Interface elements and the temperature law").  With k (Tdep - T) =
100, f = 1 exactly in double precision, and v = k_w = depositionVelocity:
the mass flux into the wall is rho k_w Y per unit wall area.  The two wall
laws are two discretisations of ONE continuous condition,

    -rho D dY/dr (R) = rho k_w Y(R):

  none      the flux k_w Y_c with the gas of the wall cell (first order:
            Y_c - Y(R) = O(d) with d the half-cell distance);
  halfCell  the flux Y_c/(1/k_w + d/D): the face value of the Robin
            condition with a linear profile over the half cell (second
            order).

With constant rho, T, p (so D uniform) and u = 2 U_b (1 - r^2/R^2) the steady
equation u dY/dx = D (d2Y/dx2 + (1/r) d/dr (r dY/dr)) has the separated
solutions Y = phi(r) exp(-mu x):

    phi'' + phi'/r + (mu u(r)/D + mu^2) phi = 0,
    phi'(0) = 0,   D phi'(R) = -k_w phi(R),

i.e. with eta = r/R, m = mu R, Pe = U_b R/D, Bi = k_w R/D

    phi'' + phi'/eta + (2 m Pe (1 - eta^2) + m^2) phi = 0,
    phi'(0) = 0,   phi'(1) + Bi phi(1) = 0.

The smallest positive m (the slowest downstream-decaying mode; the m^2
term is the axial diffusion) is the decay rate of every cross-sectional
measure (gas line density, mixing-cup concentration) far from the inlet
and the outlet.

Methods (independent; the error bound is their largest relative
difference):

  series     phi = sum a_k eta^(2k), a_0 = 1, a_(k+1) = -(A a_k +
             B a_(k-1))/(2k+2)^2, A = 2 m Pe + m^2, B = -2 m Pe (an entire
             function), in mpmath at 50 digits; the first root of
             F(m) = phi'(1) + Bi phi(1) after a scan from m = 0 (F(0) = Bi
             > 0), refined by mpmath (illinois, 1e-45)
  shooting   scipy DOP853 (rtol 1e-13, atol 1e-16) from eta0 = 1e-3 (start
             values from three series terms, truncation < 1e-18); brentq
  chebyshev  collocation in s = eta^2 (4 s psi'' + 4 psi' + (A + B s) psi =
             0, regular at s = 0; 2 psi'(1) + Bi psi(1) = 0) on N + 1
             Chebyshev points, companion linearisation of the quadratic
             eigenproblem (scipy.linalg.eig), the smallest positive real
             eigenvalue refined by Newton on the bordered collocation system

Discrete predictor (not the reference: it designs the meshes and checks
that the solver solves the documented discretisation).  The wedge of
graetzWedge has its vertices at (x, r, -+r tan(theta/2)), so its face
areas (2 t r dx radial, t (r2^2 - r1^2) axial), volumes t (r2^2 - r1^2) dx
and centroids y = (2/3)(r2^3 - r1^3)/(r2^2 - r1^2) are those of an
axisymmetric finite-volume scheme in the radius y (t = tan(theta/2)
cancels).  On uniform radial (Nr) and axial (dx) spacing, with phi the
exact flux of u(y) through each axial face, Gauss linear corrected
diffusion (orthogonal mesh), Gauss linear or SuperBee convection and the
wall law, the developed steady solution is Y_ij = phi_i lambda^j; lambda
makes the tridiagonal matrix

    M(lambda) = diag(Q_i g(lambda) - c_i h(lambda)) + dx K

singular, with Q_i the volumetric flux of ring i, c_i = S_i/dx,
h = lambda - 2 + 1/lambda, K the radial operator plus the wall uptake
2 R k_eff in the last ring (k_eff = k_w for none, 1/(1/k_w + d/D) for
halfCell, d = R - y_last), and g = (lambda - 1/lambda)/2 for Gauss linear
or (3 - 4/lambda + 1/lambda^2)/2 for SuperBee in the developed region
(r = exp(mu dx) in [1, 2]: SuperBee's psi = r is the linear-upwind face
value Y_C + (Y_C - Y_U)/2).  mu_h = -ln(lambda)/dx of the largest root
lambda < 1, found as the first sign change of det M(exp(-mu dx)) in mu
(continuant recurrence) and bisected to machine precision.

GCI (Celik et al. 2008, J. Fluids Eng. 130, 078001), constant ratio r:
p = |ln|e32/e21||/ln r, f_ext = (r^p f1 - f2)/(r^p - 1),
e_a = |f1 - f2|/|f1|, GCI_fine = 1.25 e_a/(r^p - 1).

Usage
  graetz.py reference <Pe> <Bi> [--no-axial]
  graetz.py record    <Pe> <L/R> <Bi> [<Bi> ...]       (JSON on stdout)
  graetz.py discrete  <Pe> <Bi> <Nr> <dx/R> <none|halfCell> [linear|SuperBee]
  graetz.py table     <Pe> <Bi> <Nr0> <levels> <dx/dr> [linear|SuperBee]
  graetz.py gci       <f_coarse> <f_medium> <f_fine> [<ratio>]

------------------------------------------------------------------------------
"""
import json
import sys

import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import mpmath as mp

mp.mp.dps = 50


# ----------------------------------------------------------------------------
# Method 1: power series (mpmath; a double-precision copy for the scan)
# ----------------------------------------------------------------------------

def series_F(m, Pe, Bi, kmax=4000):
    """F(m) = phi'(1) + Bi phi(1) of the even series with phi(0) = 1."""
    m, Pe, Bi = mp.mpf(m), mp.mpf(Pe), mp.mpf(Bi)
    A = 2*m*Pe + m*m
    B = -2*m*Pe
    tol = mp.mpf(10)**(-mp.mp.dps + 5)
    a_prev, a = mp.mpf(0), mp.mpf(1)
    phi = dphi = mp.mpf(0)
    small = 0
    for k in range(kmax):
        phi += a
        dphi += 2*k*a
        a_prev, a = a, -(A*a + B*a_prev)/((2*k + 2)**2)
        if abs(a) < tol*(1 + abs(phi)) and abs(a_prev) < tol*(1 + abs(phi)):
            small += 1
            if small > 3:
                break
        else:
            small = 0
    return dphi + Bi*phi


def series_F_float(m, Pe, Bi, kmax=4000):
    A = 2*m*Pe + m*m
    B = -2*m*Pe
    a_prev, a, phi, dphi = 0.0, 1.0, 0.0, 0.0
    for k in range(kmax):
        phi += a
        dphi += 2*k*a
        a_prev, a = a, -(A*a + B*a_prev)/((2*k + 2)**2)
        if k > 4 and abs(a) < 1e-30*(1 + abs(phi)) \
                and abs(a_prev) < 1e-30*(1 + abs(phi)):
            break
    return dphi + Bi*phi


def series_root(Pe, Bi, mmax=None, n_scan=4000):
    """The smallest positive root of F (an mpmath number)."""
    mmax = mmax or max(20.0, 40.0/Pe)
    grid = np.linspace(1e-9, mmax, n_scan)
    vals = [series_F_float(g, Pe, Bi) for g in grid]
    for i in range(len(grid) - 1):
        if vals[i]*vals[i + 1] < 0:
            f = lambda m: series_F(m, Pe, Bi)
            return mp.findroot(f, (mp.mpf(grid[i]), mp.mpf(grid[i + 1])),
                               solver='illinois', tol=mp.mpf(10)**-45)
    raise RuntimeError("no root found below m = %g" % mmax)


# ----------------------------------------------------------------------------
# Method 2: shooting (DOP853)
# ----------------------------------------------------------------------------

def shoot_F(m, Pe, Bi, eta0=1e-3, rtol=1e-13, atol=1e-16):
    A = 2*m*Pe + m*m
    B = -2*m*Pe
    a0 = 1.0
    a1 = -A*a0/4.0
    a2 = -(A*a1 + B*a0)/16.0
    phi0 = a0 + a1*eta0**2 + a2*eta0**4
    dphi0 = 2*a1*eta0 + 4*a2*eta0**3
    s = solve_ivp(lambda eta, y: [y[1], -y[1]/eta - (A + B*eta*eta)*y[0]],
                  [eta0, 1.0], [phi0, dphi0], method='DOP853',
                  rtol=rtol, atol=atol)
    return s.y[1, -1] + Bi*s.y[0, -1]


def shooting_root(Pe, Bi, bracket):
    return brentq(shoot_F, bracket[0], bracket[1], args=(Pe, Bi),
                  xtol=1e-15, rtol=4*np.finfo(float).eps, maxiter=200)


# ----------------------------------------------------------------------------
# Method 3: Chebyshev collocation in s = eta^2
# ----------------------------------------------------------------------------

def cheb(N):
    """Trefethen's Chebyshev differentiation matrix on [-1,1]."""
    x = np.cos(np.pi*np.arange(N + 1)/N)
    c = np.hstack([2, np.ones(N - 1), 2])*(-1)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    Dm = np.outer(c, 1/c)/(dX + np.eye(N + 1))
    Dm = Dm - np.diag(Dm.sum(axis=1))
    return Dm, x


def _collocation(Pe, Bi, N, axial):
    Dm, x = cheb(N)
    s = (x + 1)/2              # s[0] = 1 (wall), s[N] = 0 (axis)
    D1 = 2*Dm
    D2 = D1 @ D1
    n = N + 1
    L0 = 4*np.diag(s) @ D2 + 4*D1
    L1 = np.diag(2*Pe*(1 - s))
    L2 = np.eye(n) if axial else np.zeros((n, n))
    L0[0, :] = 2*D1[0, :]      # wall row: 2 psi'(1) + Bi psi(1) = 0
    L0[0, 0] += Bi
    L1[0, :] = 0
    L2[0, :] = 0
    return L0, L1, L2


def chebyshev_modes(Pe, Bi, N=48, axial=True):
    """All real eigenvalues m of the collocated problem, sorted."""
    L0, L1, L2 = _collocation(Pe, Bi, N, axial)
    n = N + 1
    if axial:
        Z = np.zeros((n, n))
        I = np.eye(n)
        w = la.eig(np.block([[Z, I], [-L0, -L1]]),
                   np.block([[I, Z], [Z, L2]]), right=False)
    else:
        w = la.eig(-L0, L1, right=False)
    w = w[np.isfinite(w)]
    real = w[np.abs(w.imag) <= 1e-9*np.maximum(1, np.abs(w.real))].real
    return np.sort(real)


def chebyshev_refine(Pe, Bi, m0, N=48, axial=True, iters=30):
    """Newton on the collocated problem bordered with psi(axis) = 1."""
    L0, L1, L2 = _collocation(Pe, Bi, N, axial)
    n = N + 1
    m = float(m0)
    psi = np.linalg.solve(L0 + m*L1 + m*m*L2 + 1e-14*np.eye(n), np.ones(n))
    psi = psi/psi[N]
    for _ in range(iters):
        Mm = L0 + m*L1 + m*m*L2
        J = np.zeros((n + 1, n + 1))
        J[:n, :n] = Mm
        J[:n, n] = (L1 + 2*m*L2) @ psi
        J[n, N] = 1.0
        delta = np.linalg.solve(J, -np.concatenate([Mm @ psi,
                                                    [psi[N] - 1.0]]))
        psi = psi + delta[:n]
        m_new = m + delta[n]
        done = abs(m_new - m) <= 4*np.finfo(float).eps*abs(m_new)
        m = m_new
        if done:
            break
    return m


# ----------------------------------------------------------------------------
# Reference record and fit window
# ----------------------------------------------------------------------------

def reference(Pe, Bi, axial=True, N_cheb=(32, 48)):
    out = {}
    if axial:
        m_series = series_root(Pe, Bi)
        out['series'] = m_series
        m0 = float(m_series)
        out['shooting'] = shooting_root(Pe, Bi, (m0*(1 - 1e-6),
                                                 m0*(1 + 1e-6)))
    cheb_vals = {}
    for N in N_cheb:
        modes = chebyshev_modes(Pe, Bi, N, axial)
        pos = modes[modes > 1e-12]
        neg = modes[modes < -1e-12]
        cheb_vals[N] = (chebyshev_refine(Pe, Bi, pos[0], N, axial),
                        pos[1] if len(pos) > 1 else np.nan,
                        neg[-1] if len(neg) else np.nan)
    out['chebyshev'] = cheb_vals
    return out


def fit_window(Pe, Bi, L_over_R, contamination=1e-10, dynamic=100.0):
    """The developed fit window [xs, xe] in units of R: the second
    downstream mode m2 changes the local decay rate by at most
    `contamination` (unit amplitude ratio) after xs; the first upstream
    mode m(-1), through which the outlet acts, by at most `contamination`
    before xe; and the signal drops by at most `dynamic` over the window
    (the absolute error floor of the linear solver is relative to the
    inlet value)."""
    modes = chebyshev_modes(Pe, Bi, 48)
    pos = modes[modes > 1e-12]
    neg = modes[modes < -1e-12]
    m1, m2, mneg = pos[0], pos[1], neg[-1]
    xs = np.log((m2 - m1)/m1/contamination)/(m2 - m1)
    xe = min(L_over_R - np.log(1/contamination)/(abs(mneg) + m1),
             xs + np.log(dynamic)/m1)
    return float(xs), float(xe), dict(m1=float(m1), m2=float(m2),
                                      mneg=float(mneg))


def reference_record(Pe, Bi, L_over_R):
    """The reference m by the three methods, its error bound and the fit
    window of one (Pe, Bi)."""
    ref = reference(Pe, Bi, True)
    m = float(ref['series'])
    methods = {'series': m, 'shooting': float(ref['shooting'])}
    for N, v in ref['chebyshev'].items():
        methods['chebyshev%d' % N] = float(v[0])
    bound = max(abs(v - m)/m for v in methods.values())
    xs, xe, spec = fit_window(Pe, Bi, L_over_R)
    return dict(Pe=Pe, Bi=Bi, m=m, m_text=mp.nstr(ref['series'], 25),
                methods=methods, bound=bound, window=[xs, xe],
                spectrum=spec, L_over_R=L_over_R)


# ----------------------------------------------------------------------------
# Discrete predictor: the developed mode of the wedge finite-volume scheme
# ----------------------------------------------------------------------------

def wedge_rings(Nr):
    """Uniform radial rings of the unit radius: faces, x-face areas (per
    t), centroids."""
    r = np.linspace(0.0, 1.0, Nr + 1)
    S = r[1:]**2 - r[:-1]**2
    y = (2.0/3.0)*(r[1:]**3 - r[:-1]**3)/S
    return r, S, y


def half_cell(Nr):
    """d/R of the wall face of a uniform radial mesh."""
    return 1.0 - wedge_rings(Nr)[2][-1]


def discrete_operator(Pe, Bi, Nr, dx, law):
    """Per ring (R = 1, D = 1, U_b = Pe; rho cancels): Q (volumetric flux
    through an x-face, per t), c = S/dx, and the tridiagonal radial
    operator K (diagonal kd, off-diagonal ko) with the wall uptake."""
    r, S, y = wedge_rings(Nr)
    Q = 4*Pe*((r[1:]**2/2 - r[1:]**4/4) - (r[:-1]**2/2 - r[:-1]**4/4))
    kd = np.zeros(Nr)
    ko = np.zeros(Nr - 1)
    for i in range(Nr - 1):
        c = 2*r[i + 1]/(y[i + 1] - y[i])
        kd[i] += c
        kd[i + 1] += c
        ko[i] = -c
    d = 1.0 - y[-1]
    if law == 'none':
        keff = Bi
    elif law == 'halfCell':
        keff = 1.0/(1.0/Bi + d)
    else:
        raise ValueError(law)
    kd[-1] += 2*keff
    return Q, S/dx, kd, ko


def _convection(lam, scheme):
    if scheme == 'linear':
        return (lam - 1/lam)/2
    if scheme == 'SuperBee':           # developed branch: linear upwind
        return (3 - 4/lam + 1/(lam*lam))/2
    raise ValueError(scheme)


def _det_sign(Q, cax, kd, ko, dx, lam, scheme):
    a = Q*_convection(lam, scheme) - cax*(lam - 2 + 1/lam) + dx*kd
    off = dx*ko
    sign = 1.0
    rprev = None
    for i in range(len(a)):
        ri = a[i] if i == 0 else a[i] - off[i - 1]*off[i - 1]/rprev
        if ri == 0.0:
            ri = 1e-300
        sign *= np.sign(ri)
        rprev = ri
    return sign


def discrete_mode(Pe, Bi, Nr, dx, law, scheme='linear', mu_max=20.0,
                  n_scan=400):
    """mu_h R of the slowest downstream mode (dx = axial spacing/R)."""
    Q, cax, kd, ko = discrete_operator(Pe, Bi, Nr, dx, law)
    f = lambda mu: _det_sign(Q, cax, kd, ko, dx, np.exp(-mu*dx), scheme)
    grid = np.linspace(1e-9, mu_max, n_scan)
    s0 = f(grid[0])
    for k in range(1, len(grid)):
        if f(grid[k]) != s0:
            lo, hi = grid[k - 1], grid[k]
            for _ in range(200):
                mid = 0.5*(lo + hi)
                if mid <= lo or mid >= hi:
                    break
                if f(mid) == s0:
                    lo = mid
                else:
                    hi = mid
            return 0.5*(lo + hi)
    raise RuntimeError("no sign change below mu = %g" % mu_max)


def discrete_mode_eig(Pe, Bi, Nr, dx, law):
    """Gauss linear only, by a dense companion eigensolve (a check of
    discrete_mode for small Nr)."""
    Q, cax, kd, ko = discrete_operator(Pe, Bi, Nr, dx, law)
    n = Nr
    K = np.diag(kd) + np.diag(ko, 1) + np.diag(ko, -1)
    A2 = np.diag(Q/2 - cax)
    A1 = dx*K + np.diag(2*cax)
    A0 = np.diag(-Q/2 - cax)
    Z = np.zeros((n, n))
    I = np.eye(n)
    w = la.eig(np.block([[Z, I], [-A0, -A1]]), np.block([[I, Z], [Z, A2]]),
               right=False)
    w = w[np.isfinite(w)]
    lam = w[(np.abs(w.imag) < 1e-10) & (w.real > 0) & (w.real < 1)].real
    return -np.log(np.max(lam))/dx


# ----------------------------------------------------------------------------
# GCI (Celik et al. 2008), constant refinement ratio
# ----------------------------------------------------------------------------

def gci(f3, f2, f1, r=2.0):
    """f3 coarse, f2 medium, f1 fine."""
    e21 = f2 - f1
    e32 = f3 - f2
    p = abs(np.log(abs(e32/e21)))/np.log(r)
    f_ext = (r**p*f1 - f2)/(r**p - 1)
    e_a = abs((f1 - f2)/f1)
    return dict(p=p, f_ext=f_ext, e_a=e_a,
                e_ext=abs((f_ext - f1)/f_ext),
                gci=1.25*e_a/(r**p - 1),
                monotone=bool(e32/e21 > 0))


# ----------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------

def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    cmd = argv[1]
    if cmd == 'reference':
        Pe, Bi = float(argv[2]), float(argv[3])
        axial = '--no-axial' not in argv
        ref = reference(Pe, Bi, axial, N_cheb=(32, 48, 64))
        print("Pe = %r, Bi = %r, axial diffusion %s" % (Pe, Bi, axial))
        if axial:
            print("  series    m = mu R = %s" % mp.nstr(ref['series'], 30))
            print("  shooting  m = mu R = %.17g" % ref['shooting'])
        for N, (m1, m2, mneg) in ref['chebyshev'].items():
            print("  chebyshev N = %3d: m1 = %.17g, m2 = %.10g, upstream "
                  "m(-1) = %.10g" % (N, m1, m2, mneg))
        if axial:
            mr = float(ref['series'])
            dev = [abs(ref['shooting'] - mr)/mr]
            dev += [abs(v[0] - mr)/mr for v in ref['chebyshev'].values()]
            print("  error bound (largest relative difference of the "
                  "methods) = %.3e" % max(dev))
        return 0
    if cmd == 'record':
        Pe, LR = float(argv[2]), float(argv[3])
        recs = [reference_record(Pe, float(b), LR) for b in argv[4:]]
        print(json.dumps(recs, indent=1))
        return 0
    if cmd == 'discrete':
        scheme = argv[7] if len(argv) > 7 else 'linear'
        print("%.17g" % discrete_mode(float(argv[2]), float(argv[3]),
                                      int(argv[4]), float(argv[5]), argv[6],
                                      scheme))
        return 0
    if cmd == 'table':
        Pe, Bi, Nr0, levels, aspect = (float(argv[2]), float(argv[3]),
                                       int(argv[4]), int(argv[5]),
                                       float(argv[6]))
        scheme = argv[7] if len(argv) > 7 else 'linear'
        mref = float(series_root(Pe, Bi))
        print("Pe %g, Bi %g, dx/dr %g, %s: reference m = %.15g"
              % (Pe, Bi, aspect, scheme, mref))
        for law in ('none', 'halfCell'):
            Ns = [Nr0*2**l for l in range(levels)]
            vals = [discrete_mode(Pe, Bi, N, aspect/N, law, scheme)
                    for N in Ns]
            for N, v in zip(Ns, vals):
                print("  %-8s Nr = %4d: m_h = %.15g, error %+.3e"
                      % (law, N, v, (v - mref)/mref))
            for l in range(levels - 2):
                g = gci(vals[l], vals[l + 1], vals[l + 2])
                print("    Nr %d/%d/%d: p = %.4f, GCI_fine = %.4f %%, "
                      "f_ext error %+.2e" % (Ns[l], Ns[l + 1], Ns[l + 2],
                                             g['p'], 100*g['gci'],
                                             (g['f_ext'] - mref)/mref))
        return 0
    if cmd == 'gci':
        r = float(argv[5]) if len(argv) > 5 else 2.0
        print(json.dumps(gci(float(argv[2]), float(argv[3]),
                             float(argv[4]), r)))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
