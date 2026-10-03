/*------------------------------------------------------------------------------

phaseChangeKinetics.H

PURPOSE

Local gas/condensate exchange kinetics of rhoFixedFlowFoam, written as small
standalone C++ functions in the style of PbI2HeDiffusivity.H.  They use no
OpenFOAM types, so tests/testPhaseChange.C checks them with plain g++.

Sign convention: a positive exchange is net EVAPORATION into the gas.

Interface elements.  The exchange happens at interface elements e (one wall
face of the interaction layer, or later one sample cell).  Element e lies in
gas cell c(e) and has the area A_e [m2].  Its mass flow into the gas is
affine in the gas mass ratio Y of its cell,

  q_e(Y) = a_e - g_e*Y   [kg/s],
  a_e    = A_e*alpha_e >= 0,   g_e = A_e*gamma_e >= 0,

where alpha_e [kg/(m2 s)] and gamma_e [kg/(m2 s)] are the coefficients of
the unconstrained law (plan section 2, E2):

  temperature model (Lobresco et al. 2026, Eq. 5, flux form, irreversible)
    alpha_e = 0,   gamma_e = rho_c v_d f(T_e),
    f(T)    = [1 - exp(-k (Tdep - T))]^+          (depositionActivation)

  Hertz-Knudsen-Schrage (milestone M5), SI
    alpha_e = G_e p_eq(T_e),   gamma_e = G_e beta_c,   beta_c = rho_c R T_c/M,
    G_e     = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T_e))
                                                   (hksConductance, [s/m]),
  so that q_e = A_e G_e (p_eq - p_i) with the partial pressure
  p_i = beta_c Y_c of the carrier's actual density; with the half cell in
  series, G'_e = G_e/(1 + G_e beta_c d_e/(rhoD)_e) (halfCellConductance).

The realised flow of an element must respect its bounds for the step,

  lower_e <= q_e <= upper_e   [kg/s]

  wall element     lower = -infinity (deposition is always possible),
                   upper = deposit/deltaT (0 for a bare wall; always 0 for
                   the irreversible temperature model);
  sample element   one per cell of an inventory sample (HKS):
                   lower = 0 (no deposition on the sample),
                   upper = reservoir/deltaT.

Regimes.  Every element is in one of three regimes; its contribution to the
gas equation of its cell is

  IMPLICIT  q = a - g*Y     Su += a,     Sp += g   (per cell volume)
  CAPPED    q = upper       Su += upper, Sp += 0   (the deposit is exhausted)
  NONE      q = 0           nothing

so Su >= 0 and Sp >= 0 always hold (fvm::Sp(Sp,Y) plus Su; never
fvm::SuSp): the exchange only adds V*Sp >= 0 to the diagonal and V*Su >= 0
to the source, so it never weakens the diagonal of the transport matrix nor
adds a negative source (Patankar positivity of the source term).  Whether
the whole matrix is an M-matrix is decided by the transport terms: upwind
convection without a non-orthogonal correction has non-positive
off-diagonal coefficients; a limited scheme such as SuperBee puts weight on
the downwind cell wherever its (lagged) limiter is positive, so with
SuperBee the transport matrix itself is not an M-matrix.

The regimes of the elements of one cell are chosen by the LOCAL IMPLICIT
PREDICTOR (cellPredictor, E4 of the plan) from the diagonal A_T and the
off-diagonal/source part H_T of the cell's transport-only gas equation
(fvMatrix::A() and fvMatrix::H()).  An explicit test at the old iterate would
misjudge the regime when the exchange is stiff (lambda*deltaT >> 1); the
local implicit test does not.  An element never switches directly between
its two bounds (CAPPED and NONE of a sample element) from one corrector to
the next: it passes IMPLICIT (acrossBounds), since the lagged neighbours of
a block of sample cells would otherwise make the regimes alternate.  After
the solve, the guard (interfaceExchange) enforces the bounded law at the
solved gas values, in both directions (admissibleFlow).

------------------------------------------------------------------------------*/

#ifndef phaseChangeKinetics_H
#define phaseChangeKinetics_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace LESTO {

  /* CODATA 2018 molar gas constant [J/(mol K)] */
  constexpr double universalGasConstant = 8.314462618;

  /*--------------------------------------------------------------------------
  Regime of one interface element.  UNKNOWN is used only before the first
  corrector of a fresh start, so that the first predictor counts as a change.
  --------------------------------------------------------------------------*/
  enum exchangeRegime { UNKNOWN = -1, NONE = 0, IMPLICIT = 1, CAPPED = 2 };

  /*--------------------------------------------------------------------------
  Activation of the temperature model, f = [1 - exp(-k (Tdep - T))]^+ [-].
  T and Tdep in K, k in 1/K.
  --------------------------------------------------------------------------*/
  inline double depositionActivation (const double T,
                                      const double k,
                                      const double Tdep) {
    return T < Tdep ? 1.0 - std::exp(-k*(Tdep - T)) : 0.0;
  }

  /*--------------------------------------------------------------------------
  Series combination of a kinetic transfer velocity v [m/s] with the
  diffusion resistance of the half cell between the cell centre and the wall
  face, resistance = rho_c d_e/(rhoD)_e [s/m]:

    v' = 1/(1/v + resistance)       (wallResistance halfCell)

  This is algebraically the Robin (mixed) boundary condition of the
  face-value law, evaluated at the cell centre.  v = 0 gives 0.
  --------------------------------------------------------------------------*/
  inline double halfCellVelocity (const double v,
                                  const double resistance) {
    return v > 0.0 ? 1.0/(1.0/v + resistance) : 0.0;
  }

  /*--------------------------------------------------------------------------
  Hertz-Knudsen-Schrage conductance G [s/m]: the net evaporation mass flux
  [kg/(m2 s)] per pressure difference p_eq - p_i [Pa],

    G = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T)),

  M in kg/mol, T in K, sigma the accommodation (condensation) coefficient
  in (0, 1].  kineticScale = 1 is SI; kineticScale = 1e-5/sqrt(1000) =
  3.16228e-7 reproduces a code that inserts p in bar and M in g/mol into the
  same expression (the T-Flows reference implementation).  At 700 K and
  sigma = 1 for PbI2: G R T/M = 2 sqrt(R T/(2 pi M)) = 89.65 m/s.
  --------------------------------------------------------------------------*/
  inline double hksConductance (const double M,
                                const double T,
                                const double accommodation,
                                const double Ce,
                                const double kineticScale) {
    const double pi = 3.14159265358979323846;
    return kineticScale*Ce*2.0*accommodation/(2.0 - accommodation)
         * std::sqrt(M/(2.0*pi*universalGasConstant*T));
  }

  /*--------------------------------------------------------------------------
  HKS conductance G [s/m] in series with the diffusion of the half cell
  between the cell centre and the wall face (wallResistance halfCell):

    G' = G/(1 + G beta resistance),   resistance = d_e/(rhoD)_e [m2 s/kg],

  from G (p_eq - beta Y_w) = (rhoD/d) (Y_w - Y_c): the Robin (mixed)
  condition of the face-value law, evaluated at the cell centre.  G = 0
  gives 0.
  --------------------------------------------------------------------------*/
  inline double halfCellConductance (const double G,
                                     const double beta,
                                     const double resistance) {
    return G > 0.0 ? G/(1.0 + G*beta*resistance) : 0.0;
  }

  /*--------------------------------------------------------------------------
  Duration [s] of the overlap of the time step [t0, t1] with the release
  window [start, start + duration].  Summed over consecutive steps (t1 of one
  step is t0 of the next) it telescopes to the window length for any step
  size, so a finite-duration release injects exactly its amount.
  --------------------------------------------------------------------------*/
  inline double releaseOverlap (const double t0,
                                const double t1,
                                const double start,
                                const double duration) {
    const double end   = start + duration;
    const double begin = t0 > start ? t0 : start;
    const double stop  = t1 < end ? t1 : end;
    return stop > begin ? stop - begin : 0.0;
  }

  /*--------------------------------------------------------------------------
  Length [s] of the release window as releaseOverlap() represents it:
  (start + duration) - start, with the window end rounded as there.  The
  overlaps of consecutive steps telescope to this length, not to the
  nominal duration: the rounded end differs from start + duration by up to
  half an ulp of the end time (5.7e-14 s near t = 1000 s, i.e. 3e-12 of a
  window of 0.0198 s).  A release rate n0 M overlap/(length dt V_Z)
  therefore injects n0 M to round-off at any start time.
  --------------------------------------------------------------------------*/
  inline double releaseWindowLength (const double start,
                                     const double duration) {
    return (start + duration) - start;
  }

  /*--------------------------------------------------------------------------
  Regime of an element whose flow q [kg/s] violates its bounds:

    q > upper : CAPPED when there is something to exhaust (upper > 0),
                NONE otherwise.  An element with upper = 0 (bare wall, or
                any wall element of the irreversible temperature model)
                never becomes CAPPED, so a negative undershoot of the gas
                can never delete an existing deposit (review item B2);
    q < lower : NONE (e.g. deposition onto the sample).

  Returns IMPLICIT when q is within the bounds.  The same rule is used by
  the predictor and by the guard after the solve.
  --------------------------------------------------------------------------*/
  inline exchangeRegime boundedRegime (const double q,
                                       const double lower,
                                       const double upper) {
    if (q > upper) {
      return upper > 0.0 ? CAPPED : NONE;
    }
    if (q < lower) {
      return NONE;
    }
    return IMPLICIT;
  }

  /*--------------------------------------------------------------------------
  Admissible flow [kg/s] of an element whose unconstrained (implicit) flow
  is q: the projection min(max(q, lower), upper) onto its bounds, i.e. the
  flow the bounded law prescribes at the gas value that gives q.

  The regime of an element is COMPLEMENTARY (consistent with the gas value
  Y of its cell) when its realised flow equals the admissible flow of
  q = a - g Y:

    IMPLICIT  lower <= q <= upper
    CAPPED    q >= upper                    (the condensate is exhausted)
    NONE      q <= lower = 0                (a sample: no deposition), or
              q >= upper = 0                (a bare wall: deposition only;
                                             the temperature model), or
              lower = upper = 0             (an empty sample)

  e.g. a CAPPED sample element in supersaturated gas (q < 0) or a NONE
  sample element in undersaturated gas (0 < q) is not complementary.
  --------------------------------------------------------------------------*/
  inline double admissibleFlow (const double q,
                                const double lower,
                                const double upper) {
    return q < lower ? lower : (q > upper ? upper : q);
  }

  /*--------------------------------------------------------------------------
  Cell value of the local implicit predictor for given regimes (E4),

    Y_loc = (rhs + sum_IMPLICIT a + sum_CAPPED upper)/(diag + sum_IMPLICIT g),

  with diag, rhs, a, g and upper as in cellPredictor() below.
  --------------------------------------------------------------------------*/
  inline double cellYloc (const double  diag,
                          const double  rhs,
                          const int     n,
                          const double* a,
                          const double* g,
                          const double* upper,
                          const int*    regime) {

    double numerator   = rhs;
    double denominator = diag;

    for (int i = 0; i < n; ++i) {
      if (regime[i] == IMPLICIT) {
        numerator   += a[i];
        denominator += g[i];
      } else if (regime[i] == CAPPED) {
        numerator   += upper[i];
      }
    }

    return numerator/denominator;
  }

  /*--------------------------------------------------------------------------
  True when a switch from the regime previous to next would jump between
  the two bounds of an element, CAPPED (q = upper > 0) and NONE at the
  lower bound (q = lower = 0: a sample element), without passing IMPLICIT.
  A wall element has no finite lower bound: its NONE is the upper bound 0
  of a bare wall, on the same side as CAPPED (a wall that has a deposit
  again, e.g. set through mw_ at a restart, turns CAPPED directly).  With
  upper = 0 there is only one bound (a bare wall, an exhausted or empty
  sample), and NONE is the only admissible regime left.
  --------------------------------------------------------------------------*/
  inline bool acrossBounds (const int    previous,
                            const int    next,
                            const double lower,
                            const double upper) {
    return lower == 0.0 && upper > 0.0
        && ((previous == CAPPED && next == NONE)
         || (previous == NONE && next == CAPPED));
  }

  /*--------------------------------------------------------------------------
  Local implicit predictor of one exchange cell (plan section 2, E4).

    diag      V*A_T         [kg/s]  transport-only diagonal (fvMatrix::A()*V)
    rhs       V*(H_T+Srel)  [kg/s]  transport-only H (fvMatrix::H()*V) plus
                                    the release source of the cell
    n         number of elements of the cell
    a, g      element coefficients a_e, g_e [kg/s]
    lower,
    upper     element bounds [kg/s]
    regime    output: the regime of every element
    previous  the regimes of the previous corrector (at the first corrector
              of a step: of the previous step), or nullptr (a fresh start)

  All elements start IMPLICIT (inactive ones, a = g = 0, are NONE).  The
  cell value Y_loc (cellYloc) is evaluated and every IMPLICIT element whose
  flow at Y_loc violates a bound is switched with boundedRegime().  This is
  repeated until nothing changes.  Elements only ever leave IMPLICIT, so at
  most n+1 passes are needed.

  With the neighbour cells fixed (H_T), this is the exact solution of the
  cell.  But H_T lags the neighbours by one corrector, and a block of
  exchange cells that all switch at once (the cells of an inventory sample)
  invalidates it: a small reservoir, CAPPED in every sample cell, dumps
  into the gas of all of them, so at the next corrector every cell sees
  supersaturated neighbours and turns NONE; with NONE everywhere the gas
  drains, and every cell turns CAPPED again.  The regimes then alternate
  between the two bounds at every corrector and never settle (review of
  M5, round 2).  Therefore an element never switches DIRECTLY between its
  two bounds: when the previous regime is CAPPED and the predictor gives
  NONE at the lower bound, or the reverse, and IMPLICIT is admissible
  (upper > 0), the element is set IMPLICIT instead (acrossBounds).  The
  coupled solve then resolves the exchange of the whole block implicitly,
  and the next corrector starts from neighbour values that belong to it.
  A single cell never needs the direct switch, because its local solution
  does not change once the neighbours are right; a CAPPED element that
  exhausted its condensate (upper = 0 in the next step) turns NONE as
  before.  Returns Y_loc of the final regimes.
  --------------------------------------------------------------------------*/
  inline double cellPredictor (const double  diag,
                               const double  rhs,
                               const int     n,
                               const double* a,
                               const double* g,
                               const double* lower,
                               const double* upper,
                               int*          regime,
                               const int*    previous = nullptr) {

    for (int i = 0; i < n; ++i) {
      regime[i] = (a[i] == 0.0 && g[i] == 0.0) ? NONE : IMPLICIT;
    }

    double Yloc = 0.0;

    for (int pass = 0; pass <= n; ++pass) {

      Yloc = cellYloc(diag, rhs, n, a, g, upper, regime);

      bool changed = false;
      for (int i = 0; i < n; ++i) {
        if (regime[i] == IMPLICIT) {
          const exchangeRegime r =
            boundedRegime(a[i] - g[i]*Yloc, lower[i], upper[i]);
          if (r != IMPLICIT) {
            regime[i] = r;
            changed = true;
          }
        }
      }

      if (!changed) {
        break;
      }
    }

    if (previous) {
      bool damped = false;
      for (int i = 0; i < n; ++i) {
        if (acrossBounds(previous[i], regime[i], lower[i], upper[i])) {
          regime[i] = IMPLICIT;
          damped = true;
        }
      }
      if (damped) {
        Yloc = cellYloc(diag, rhs, n, a, g, upper, regime);
      }
    }

    return Yloc;
  }

  /*--------------------------------------------------------------------------
  Guard (E5): which CAPPED or NONE elements go back to IMPLICIT.  misses
  are the missed masses [kg] of the CAPPED and NONE elements (all ranks)
  whose realised flow differs from their admissible flow, budget the mass
  the guard accepts, guardTolerance times the inventory [kg].  Returns the
  threshold tau: the elements whose miss is at least tau are switched back,
  the largest first, so that the misses left add up to at most budget.
  When all misses add up to at most budget nothing is switched (+infinity),
  and every miss above budget is always switched.

  The guard thus acts on the total, as it is triggered by the total (review
  of M5, round 3): before, it switched back only elements whose own miss
  exceeded budget, so several smaller misses (two bare wall elements just
  above Y_eq, 1.43 times budget together) were accepted without a guard
  pass.  The smallest misses, those of elements at their bound to
  round-off, stay last: switched back, they would only return with a bound
  violation of the same size.  misses is taken by value and sorted.
  --------------------------------------------------------------------------*/
  inline double switchBackThreshold (std::vector<double> misses,
                                     const double        budget) {
    std::sort(misses.begin(), misses.end());
    double left = 0.0;
    for (const double miss : misses) {
      if (left + miss > budget) {
        return miss;
      }
      left += miss;
    }
    return HUGE_VAL;
  }

  /*--------------------------------------------------------------------------
  Realised flow [kg/s] of one element in the given regime at the cell value
  Y of the solved gas equation.
  --------------------------------------------------------------------------*/
  inline double elementFlow (const int    regime,
                             const double a,
                             const double g,
                             const double upper,
                             const double Y) {
    if (regime == IMPLICIT) {
      return a - g*Y;
    }
    if (regime == CAPPED) {
      return upper;
    }
    return 0.0;
  }
}

#endif
