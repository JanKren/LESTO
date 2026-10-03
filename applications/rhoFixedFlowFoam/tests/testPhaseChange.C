/*------------------------------------------------------------------------------

tests/testPhaseChange.C

PURPOSE

Standalone test of phaseChangeKinetics.H and compensatedSum.H (no OpenFOAM
libraries needed): the temperature-model activation, the half-cell series
velocity, the exact release overlap, the bound rule of the regimes (review
item B2), the local implicit predictor of one exchange cell (plan section
2, E4), including the CAPPED branch that the temperature model never
reaches, the order independence of the compensated inventory sums
(binary restart bound after decomposePar or renumberMesh), and the total
of a release over many steps on the solver's clock: the per-step releases
n0 M overlap/length, accumulated with compensatedAdd as the ledger does,
sum to n0 M to a few eps for long windows (6100 and 136364 steps) and for
windows far from t = 0 (start 60 s and 1000 s), where a plain running
total or the nominal duration miss by up to 1.6e-12 (review M2, round 2).

Milestone M5 (criterion M5.0): the Hertz-Knudsen-Schrage conductance
(G R T/M = 89.65 m/s for PbI2 at 700 K; kineticScale 1e-5/sqrt(1000)
reproduces the T-Flows expression in bar and g/mol), the half-cell series
conductance against the Robin condition solved for the wall value, and the
predictor of HKS wall and sample elements: a bare wall deposits only, a
wall with a deposit evaporates up to it (CAPPED), a sample element
evaporates only (lower = 0: NONE above Y_eq), up to its reservoir
(CAPPED), and is NONE when empty.  Review of M5, round 2: the admissible
flow min(max(q, lower), upper) of the complementarity check and the
guard; the predictor with the previous regimes never jumps between the
two bounds of a sample element (CAPPED <-> NONE passes IMPLICIT; an
exhausted element still turns NONE); and a row of 10 stiff sample cells
with a small reservoir, run as the corrector loop runs (predictor with
the neighbours of the previous corrector, exact coupled solve): without
the previous regimes the regimes alternate between all CAPPED and all
NONE at every corrector and the final state is not complementary; with
them they settle after a few correctors, complementary.  Review of M5,
round 3: the guard switches back CAPPED and NONE elements, the largest
misses first, until the misses left add up to at most its tolerance
(switchBackThreshold; before, two sub-threshold misses of 1.43 times the
tolerance together were accepted).

Build and run from the solver directory:

  g++ -std=c++17 -Wall -Wextra tests/testPhaseChange.C -o /tmp/testPhaseChange
  /tmp/testPhaseChange

Criteria M2.0 of tests/Alltest M2 and M5.0 of tests/Alltest M5.

------------------------------------------------------------------------------*/

#include "../phaseChangeKinetics.H"
#include "../compensatedSum.H"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

namespace {

  int nFailed = 0;
  int nChecks = 0;

  void check(const bool ok, const char* what) {
    ++nChecks;
    if (!ok) {
      ++nFailed;
      std::printf("FAILED: %s\n", what);
    }
  }

  bool close(const double a, const double b, const double tol) {
    return std::abs(a - b) <= tol*std::abs(b);
  }
}

int main() {

  using namespace LESTO;

  /*--------------------------------------------------------------------------
  Activation f = 1 - exp(-k (Tdep - T)), zero at and above Tdep.
  --------------------------------------------------------------------------*/
  check(close(depositionActivation(500, 8e-3, 680), 1 - std::exp(-1.44), 1e-15),
        "activation at 500 K");
  check(depositionActivation(680, 8e-3, 680) == 0, "activation at Tdep");
  check(depositionActivation(900, 8e-3, 680) == 0, "activation above Tdep");

  /* closed box of the plan: lambda dt = A f dt = 0.671504 */
  check(close(220*0.004*depositionActivation(500, 8e-3, 680), 0.671504, 1e-6),
        "box lambda dt = 0.671504");

  /*--------------------------------------------------------------------------
  Half cell in series: 1/(1/v + r); zero for v = 0.
  --------------------------------------------------------------------------*/
  check(close(halfCellVelocity(0.01, 50), 1/(100.0 + 50.0), 1e-15),
        "half-cell series velocity");
  check(halfCellVelocity(0, 50) == 0, "half-cell with v = 0");

  /*--------------------------------------------------------------------------
  Release overlap: a window of 40.5 steps starting inside a step sums to
  its duration.
  --------------------------------------------------------------------------*/
  {
    const double dt = 1e-3, start = 3e-4, duration = 0.0405;
    double sum = 0;
    for (int n = 0; n < 100; ++n) {
      sum += releaseOverlap(n*dt, (n + 1)*dt, start, duration);
    }
    check(close(sum, duration, 1e-13), "release overlap sums to duration");
    check(releaseOverlap(0.0, dt, start, duration) == dt - start,
          "partial first step");
    check(releaseOverlap(0.05, 0.051, start, duration) == 0,
          "step after the window");
  }

  /*--------------------------------------------------------------------------
  Bound rule (B2): above upper -> CAPPED only when upper > 0.
  --------------------------------------------------------------------------*/
  check(boundedRegime(1.0, -1e300, 0.0) == NONE,
        "q > upper = 0: NONE, a deposit is never deleted");
  check(boundedRegime(1.0, -1e300, 0.5) == CAPPED, "q > upper > 0: CAPPED");
  check(boundedRegime(-1.0, 0.0, 0.5) == NONE, "q < lower: NONE");
  check(boundedRegime(0.2, 0.0, 0.5) == IMPLICIT, "within bounds: IMPLICIT");

  /*--------------------------------------------------------------------------
  Predictor, temperature model (a = 0, upper = 0): positive H deposits
  implicitly, Y_loc = H/(A + g).
  --------------------------------------------------------------------------*/
  {
    const double a[1] = {0}, g[1] = {2}, lower[1] = {-1e300}, upper[1] = {0};
    int regime[1];
    const double Y = cellPredictor(3, 5, 1, a, g, lower, upper, regime);
    check(regime[0] == IMPLICIT, "temperature, H > 0: IMPLICIT");
    check(close(Y, 1.0, 1e-15), "temperature, H > 0: Y_loc = H/(A+g)");
  }

  /*--------------------------------------------------------------------------
  B2: an element with a deposit but upper = 0 (temperature model) and a
  negative predictor (H < 0) is NONE, not CAPPED.
  --------------------------------------------------------------------------*/
  {
    const double a[1] = {0}, g[1] = {2}, lower[1] = {-1e300}, upper[1] = {0};
    int regime[1];
    const double Y = cellPredictor(3, -5, 1, a, g, lower, upper, regime);
    check(regime[0] == NONE, "temperature, H < 0: NONE (deposit kept)");
    check(close(Y, -5.0/3.0, 1e-15), "temperature, H < 0: Y_loc = H/A");
  }

  /*--------------------------------------------------------------------------
  HKS-like element with a small deposit (upper > 0) and strong evaporation:
  CAPPED, and Y_loc = (H + upper)/A.
  --------------------------------------------------------------------------*/
  {
    const double a[1] = {100}, g[1] = {1}, lower[1] = {-1e300};
    const double upper[1] = {0.5};
    int regime[1];
    const double Y = cellPredictor(2, 1, 1, a, g, lower, upper, regime);
    check(regime[0] == CAPPED, "evaporation above the deposit: CAPPED");
    check(close(Y, 1.5/2.0, 1e-15), "CAPPED: Y_loc = (H + upper)/A");
    check(elementFlow(regime[0], a[0], g[0], upper[0], Y) == upper[0],
          "CAPPED flow = upper");
  }

  /*--------------------------------------------------------------------------
  Bare wall (upper = 0) that would evaporate: NONE.
  --------------------------------------------------------------------------*/
  {
    const double a[1] = {100}, g[1] = {1}, lower[1] = {-1e300};
    const double upper[1] = {0};
    int regime[1];
    cellPredictor(2, 1, 1, a, g, lower, upper, regime);
    check(regime[0] == NONE, "bare wall, evaporation: NONE");
  }

  /*--------------------------------------------------------------------------
  Sample (lower = 0) with deposition predicted: NONE.
  --------------------------------------------------------------------------*/
  {
    const double a[1] = {1}, g[1] = {1}, lower[1] = {0}, upper[1] = {10};
    int regime[1];
    cellPredictor(1, 9, 1, a, g, lower, upper, regime);
    check(regime[0] == NONE, "sample, deposition: NONE");
  }

  /*--------------------------------------------------------------------------
  Two elements in one cell: element 0 exhausts its deposit (CAPPED),
  element 1 deposits (IMPLICIT); inactive element 2 is NONE.  The returned
  Y_loc is consistent with the final regimes.
  --------------------------------------------------------------------------*/
  {
    const double a[3] = {50, 0, 0}, g[3] = {1, 3, 0};
    const double lower[3] = {-1e300, -1e300, -1e300};
    const double upper[3] = {0.25, 0, 0};
    int regime[3];
    const double Y = cellPredictor(2, 4, 3, a, g, lower, upper, regime);
    check(regime[0] == CAPPED && regime[1] == IMPLICIT && regime[2] == NONE,
          "two elements: CAPPED + IMPLICIT, inactive NONE");
    check(close(Y, (4 + 0.25)/(2 + 3.0), 1e-15),
          "two elements: Y_loc = (H + upper_0)/(A + g_1)");
    const double q1 = elementFlow(regime[1], a[1], g[1], upper[1], Y);
    check(q1 <= upper[1] && q1 >= lower[1], "two elements: IMPLICIT in bounds");
  }

  /*--------------------------------------------------------------------------
  HKS conductance (M5).  G R T/M = 2 sqrt(R T/(2 pi M)) for sigma = 1: 89.65
  m/s for PbI2 at 700 K (plan E2).  The T-Flows expression is the molar flux
  k (p_eq - p) with k = 2 sigma/(2 - sigma)/sqrt(2 pi R M[g/mol] T) and p in
  bar; kineticScale 1e-5/sqrt(1000) gives the same molar flux G p/M.
  --------------------------------------------------------------------------*/
  {
    const double M = 0.46100894, T = 700, R = universalGasConstant;
    const double pi = 3.14159265358979323846;
    const double G = hksConductance(M, T, 1.0, 1.0, 1.0);
    check(close(G*R*T/M, 89.65, 1e-4), "HKS: G R T/M = 89.65 m/s at 700 K");
    check(close(G, 2*std::sqrt(M/(2*pi*R*T)), 1e-15),
          "HKS: G = 2 sqrt(M/(2 pi R T)) for sigma = 1");
    check(close(hksConductance(M, T, 0.5, 2.0, 3.0),
                3*2*(2*0.5/1.5)*std::sqrt(M/(2*pi*R*T)), 1e-15),
          "HKS: kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T))");

    const double scale = 1e-5/std::sqrt(1000.0);
    const double p = 3.51254e-4;                          /* bar */
    const double kTFlows = 2.0/std::sqrt(2*pi*R*(M*1000)*T);
    check(close(hksConductance(M, T, 1.0, 1.0, scale)*(p*1e5)/M,
                kTFlows*p, 1e-14),
          "HKS: kineticScale 1e-5/sqrt(1000) = the T-Flows expression");
  }

  /*--------------------------------------------------------------------------
  Half cell in series with HKS: G' (p_eq - beta Y_c) equals the flux of the
  Robin condition G (p_eq - beta Y_w) = (rhoD/d)(Y_w - Y_c) solved for Y_w.
  --------------------------------------------------------------------------*/
  {
    const double G = 3e-3, beta = 880, pEq = 35.9, Yc = 0.01;
    const double rhoD = 6e-6, d = 2.873e-5;
    const double Yw = (G*pEq + rhoD/d*Yc)/(G*beta + rhoD/d);
    const double robin = rhoD/d*(Yw - Yc);
    const double Gs = halfCellConductance(G, beta, d/rhoD);
    check(close(Gs*(pEq - beta*Yc), robin, 1e-13),
          "HKS half cell: G' = G/(1 + G beta d/rhoD) is the Robin flux");
    check(halfCellConductance(0, beta, d/rhoD) == 0, "half cell with G = 0");
  }

  /*--------------------------------------------------------------------------
  HKS elements in the predictor.  Law q = a - g Y with a = g Y_eq, Y_eq =
  0.04, g = 1000 (stiff: g >> A = 1).
  --------------------------------------------------------------------------*/
  {
    const double g[1] = {1000}, a[1] = {40};
    const double wall[1] = {-1e300}, sampleLower[1] = {0};
    int regime[1];

    /* bare wall (upper 0), gas above Y_eq: deposits (IMPLICIT) */
    const double none[1] = {0};
    double Y = cellPredictor(1, 0.1, 1, a, g, wall, none, regime);
    check(regime[0] == IMPLICIT && close(Y, (0.1 + 40)/1001.0, 1e-15),
          "HKS bare wall, Y > Y_eq: deposition, IMPLICIT");

    /* bare wall, gas below Y_eq: NONE (deposit only) */
    Y = cellPredictor(1, 0.01, 1, a, g, wall, none, regime);
    check(regime[0] == NONE && close(Y, 0.01, 1e-15),
          "HKS bare wall, Y < Y_eq: NONE (no evaporation)");

    /* wall with an ample deposit, gas below Y_eq: evaporation, IMPLICIT */
    const double ample[1] = {1e3};
    Y = cellPredictor(1, 0.0, 1, a, g, wall, ample, regime);
    check(regime[0] == IMPLICIT && close(Y, 40/1001.0, 1e-15),
          "HKS wall with an ample deposit: evaporation, IMPLICIT");

    /* wall with a small deposit: CAPPED, Y_loc = upper/A */
    const double small[1] = {0.01};
    Y = cellPredictor(1, 0.0, 1, a, g, wall, small, regime);
    check(regime[0] == CAPPED && close(Y, 0.01, 1e-15),
          "HKS wall with a small deposit: CAPPED");

    /* sample element: evaporation up to the reservoir */
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, ample, regime);
    check(regime[0] == IMPLICIT, "sample, Y < Y_eq: evaporation, IMPLICIT");
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, small, regime);
    check(regime[0] == CAPPED && close(Y, 0.01, 1e-15),
          "sample with a small reservoir: CAPPED");

    /* sample element above Y_eq: no deposition onto the sample (NONE) */
    Y = cellPredictor(1, 0.1, 1, a, g, sampleLower, ample, regime);
    check(regime[0] == NONE && close(Y, 0.1, 1e-15),
          "sample, Y > Y_eq: NONE (no deposition onto the sample)");

    /* empty sample (upper 0): NONE */
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, none, regime);
    check(regime[0] == NONE, "empty sample: NONE");

    /*------------------------------------------------------------------------
    The previous regimes (review of M5, round 2): an element never jumps
    between its two bounds, CAPPED and NONE; it passes IMPLICIT.
    ------------------------------------------------------------------------*/
    const int wasNone[1] = {NONE}, wasCapped[1] = {CAPPED};
    const int wasImplicit[1] = {IMPLICIT}, wasUnknown[1] = {UNKNOWN};

    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, small, regime, wasNone);
    check(regime[0] == IMPLICIT && close(Y, 40/1001.0, 1e-15),
          "sample, NONE before, CAPPED now: IMPLICIT (Y_loc of IMPLICIT)");
    Y = cellPredictor(1, 0.1, 1, a, g, sampleLower, small, regime, wasCapped);
    check(regime[0] == IMPLICIT && close(Y, (0.1 + 40)/1001.0, 1e-15),
          "sample, CAPPED before, NONE now: IMPLICIT");
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, small, regime,
                      wasImplicit);
    check(regime[0] == CAPPED, "sample, IMPLICIT before: CAPPED as before");
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, small, regime, wasUnknown);
    check(regime[0] == CAPPED, "sample, fresh start (UNKNOWN): CAPPED");
    Y = cellPredictor(1, 0.1, 1, a, g, sampleLower, ample, regime, wasNone);
    check(regime[0] == NONE, "sample, NONE before and now: NONE");
    Y = cellPredictor(1, 0.0, 1, a, g, sampleLower, none, regime, wasCapped);
    check(regime[0] == NONE,
          "exhausted sample (CAPPED before, upper 0 now): NONE, not damped");
    Y = cellPredictor(1, 0.1, 1, a, g, wall, none, regime, wasNone);
    check(regime[0] == IMPLICIT, "bare wall, NONE before: deposits as before");
    Y = cellPredictor(1, 0.01, 1, a, g, wall, none, regime, wasImplicit);
    check(regime[0] == NONE, "bare wall, IMPLICIT before: NONE as before");

    Y = cellPredictor(1, 0.0, 1, a, g, wall, small, regime, wasNone);
    check(regime[0] == CAPPED,
          "wall, NONE before (bare), a small deposit now: CAPPED, not damped "
          "(both regimes at the upper bound)");

    check(acrossBounds(CAPPED, NONE, 0.0, 0.5)
          && acrossBounds(NONE, CAPPED, 0.0, 0.5)
          && !acrossBounds(CAPPED, NONE, 0.0, 0.0)
          && !acrossBounds(NONE, CAPPED, -1e300, 0.5)
          && !acrossBounds(IMPLICIT, CAPPED, 0.0, 0.5)
          && !acrossBounds(UNKNOWN, CAPPED, 0.0, 0.5)
          && !acrossBounds(NONE, IMPLICIT, 0.0, 0.5),
          "acrossBounds: only CAPPED <-> NONE at lower = 0 with upper > 0");

    /* admissible flow and complementarity */
    check(admissibleFlow(-1.0, 0.0, 0.5) == 0.0
          && admissibleFlow(2.0, 0.0, 0.5) == 0.5
          && admissibleFlow(0.2, 0.0, 0.5) == 0.2
          && admissibleFlow(2.0, -1e300, 0.0) == 0.0
          && admissibleFlow(-2.0, -1e300, 0.0) == -2.0
          && admissibleFlow(3.0, 0.0, 0.0) == 0.0,
          "admissibleFlow = min(max(q, lower), upper)");

    /* the guard's switch back of CAPPED and NONE elements (review of M5,
       round 3): the largest misses first, until the rest is within the
       budget; the misses of the step of M5.2b that the round-2 guard
       accepted (two of 0.7 and 0.73 budgets, 1.43 together) */
    const double budget = 2.2313e-18;
    const std::vector<double> m5b ({0.7*budget, 0.73*budget});
    const double tau = switchBackThreshold(m5b, budget);
    check(tau == 0.73*budget,
          "switch back: two sub-threshold misses above the budget together, "
          "the larger is switched back");
    check(std::isinf(switchBackThreshold({0.4*budget, 0.5*budget}, budget))
          && std::isinf(switchBackThreshold({}, budget)),
          "switch back: misses within the budget, nothing switched");
    const std::vector<double> mixed ({3*budget, 0.2*budget, 0.9*budget,
                                      1e-12*budget, 0.5*budget});
    const double tauMixed = switchBackThreshold(mixed, budget);
    double left = 0;
    int nBack = 0;
    for (const double miss : mixed) {
      if (miss >= tauMixed) {
        ++nBack;
      } else {
        left += miss;
      }
    }
    check(tauMixed == 0.9*budget && nBack == 2 && left <= budget,
          "switch back: the misses above the budget and the next largest, "
          "the rest (0.7 budget, incl. the round-off one) left");
    check(switchBackThreshold({0.2*budget, 5*budget}, budget) == 5*budget,
          "switch back: a miss above the budget alone is switched, as in "
          "round 2");
  }

  /*--------------------------------------------------------------------------
  A block of sample cells (review of M5, round 2).  A row of n cells, each
  with one stiff sample element (g = 1000, Y_eq = 0.04) and a small
  reservoir, upper = 0.2: more than the storage of one cell at equilibrium
  (s Y_eq = 0.04), less than what a cell next to empty neighbours would
  evaporate (about (s + 2D) Y_eq = 0.44).  Storage s = 1 (Y^n = 0),
  diffusion D = 5 to each neighbour and to empty gas (Y = 0) beyond both
  ends.  Every corrector runs the predictor of every cell with the
  neighbours of the previous corrector (H = D (Y_left + Y_right), the
  transport-only matrix of the latest iterate), then solves the coupled
  tridiagonal system exactly, as the solver's corrector loop does.

  Without the previous regimes (the round-2 predictor) the regimes alternate:
  all CAPPED (neighbours empty), then all NONE (neighbours supersaturated),
  and so on; neither state is complementary, and the last corrector's
  parity decides.  With them the block settles within a few correctors in
  a complementary state (IMPLICIT inside, CAPPED at the ends).
  --------------------------------------------------------------------------*/
  {
    const int n = 10, nCorr = 20;
    const double s = 1, D = 5, gs = 1000, Yeq = 0.04, as = gs*Yeq;
    const double lowerS = 0, upperS = 0.2;

    /* returns the number of correctors with regime changes (the last one
       counted), the final regimes and whether the final state is
       complementary */
    auto run = [&](const bool damped, std::vector<int>& regimes,
                   int& lastChange, bool& complementary) {
      std::vector<double> Y(n, 0.0);
      std::vector<int> previous(n, UNKNOWN);
      regimes.assign(n, UNKNOWN);
      lastChange = -1;
      for (int corr = 0; corr < nCorr; ++corr) {
        std::vector<double> Su(n), Sp(n);
        bool changed = false;
        for (int i = 0; i < n; ++i) {
          const double left = i > 0 ? Y[i-1] : 0.0;
          const double right = i < n - 1 ? Y[i+1] : 0.0;
          int r;
          cellPredictor(s + 2*D, D*(left + right), 1, &as, &gs, &lowerS,
                        &upperS, &r, damped ? &previous[i] : nullptr);
          changed = changed || r != regimes[i];
          regimes[i] = r;
          Su[i] = r == IMPLICIT ? as : r == CAPPED ? upperS : 0.0;
          Sp[i] = r == IMPLICIT ? gs : 0.0;
        }
        previous = regimes;
        if (changed) {
          lastChange = corr;
        }
        /* (s + 2D + Sp) Y_i - D (Y_i-1 + Y_i+1) = Su_i, Thomas algorithm */
        std::vector<double> c(n), d(n);
        for (int i = 0; i < n; ++i) {
          const double b = s + 2*D + Sp[i] - (i > 0 ? -D*c[i-1] : 0.0);
          c[i] = -D/b;
          d[i] = (Su[i] + (i > 0 ? D*d[i-1] : 0.0))/b;
        }
        Y[n-1] = d[n-1];
        for (int i = n - 2; i >= 0; --i) {
          Y[i] = d[i] - c[i]*Y[i+1];
        }
      }
      complementary = true;
      for (int i = 0; i < n; ++i) {
        const double q = as - gs*Y[i];
        const double flow = elementFlow(regimes[i], as, gs, upperS, Y[i]);
        complementary = complementary
          && std::abs(admissibleFlow(q, lowerS, upperS) - flow) <= 1e-12;
      }
    };

    std::vector<int> regimes;
    int lastChange;
    bool complementary;

    run(false, regimes, lastChange, complementary);
    check(lastChange == nCorr - 1 && !complementary,
          "sample block, round-2 predictor: regimes alternate at every "
          "corrector, final state not complementary");

    run(true, regimes, lastChange, complementary);
    const int nCapped = int(std::count(regimes.begin(), regimes.end(), CAPPED));
    check(lastChange >= 0 && lastChange <= 5 && complementary
          && regimes.front() == CAPPED && regimes.back() == CAPPED
          && regimes[n/2] == IMPLICIT && nCapped < n,
          "sample block, predictor with the previous regimes: settled in at "
          "most 6 correctors, complementary, CAPPED at the ends only");
    std::printf("sample block of %d cells: settled after corrector %d, "
                "%d CAPPED\n", n, lastChange, nCapped);
  }

  /*--------------------------------------------------------------------------
  Compensated sums.  (1) Terms of very different size: 1 followed by 1e5
  terms of 1e-16.  A plain sum in this order loses every small term (1 +
  1e-16 rounds to 1); in the reverse order it keeps them.  The compensated
  sum gives 1 + 1e-11 in both orders.
  --------------------------------------------------------------------------*/
  {
    const int n = 100000;
    compensatedSum forward, reverse;
    double plainForward = 1.0, plainReverse = 0.0;
    forward.add(1.0);
    for (int i = 0; i < n; ++i) {
      forward.add(1e-16);
      plainForward += 1e-16;
    }
    for (int i = 0; i < n; ++i) {
      reverse.add(1e-16);
      plainReverse += 1e-16;
    }
    reverse.add(1.0);
    plainReverse += 1.0;

    check(plainForward == 1.0 && !close(plainReverse, 1.0, 1e-12),
          "plain sums depend on the order (1 + 1e5 x 1e-16)");
    check(close(forward.value(), 1.0 + 1e-11, 2.3e-16)
       && close(reverse.value(), 1.0 + 1e-11, 2.3e-16),
          "compensated sum of 1 + 1e5 x 1e-16 = 1 + 1e-11 in both orders");
  }

  /*--------------------------------------------------------------------------
  (2) 4e5 terms like rho*Y*V of a mesh (positive, spread over 4 decades,
  a few tiny negative ones), summed forward, backward, and as 8 blocks
  whose (sum, compensation) pairs are added in block order, as the
  reduction over 8 ranks does (interfaceExchange.C, globalSum).  The
  compensated results agree to 2 ulp; the plain sums are printed.
  --------------------------------------------------------------------------*/
  {
    const int n = 400000;
    std::mt19937_64 generator(20260924);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
      x[i] = unit(generator)*std::pow(10.0, -4.0*unit(generator))*1e-12;
      if (i % 1000 == 0) {
        x[i] = -1e-40*unit(generator);
      }
    }

    compensatedSum forward, backward;
    double plainForward = 0, plainBackward = 0;
    for (int i = 0; i < n; ++i) {
      forward.add(x[i]);
      plainForward += x[i];
    }
    for (int i = n - 1; i >= 0; --i) {
      backward.add(x[i]);
      plainBackward += x[i];
    }

    const int nBlocks = 8;
    compensatedSum total;
    double plainBlocks = 0;
    for (int b = 0; b < nBlocks; ++b) {
      compensatedSum block;
      double plainBlock = 0;
      for (int i = b*n/nBlocks; i < (b + 1)*n/nBlocks; ++i) {
        block.add(x[i]);
        plainBlock += x[i];
      }
      total.add(block.sum());
      total.add(block.compensation());
      plainBlocks += plainBlock;
    }

    const double sum = forward.value();
    const double ulp2 = 2*std::nextafter(sum, 2*sum) - 2*sum;
    const double worst = std::max (
      std::abs(backward.value() - sum), std::abs(total.value() - sum)
    );
    const double worstPlain = std::max (
      std::abs(plainBackward - plainForward),
      std::abs(plainBlocks - plainForward)
    );
    std::printf("4e5 terms: compensated orders differ by %.3e, plain orders "
                "by %.3e (relative)\n", worst/sum, worstPlain/sum);
    check(worst <= ulp2,
          "compensated sums of 4e5 terms: forward, backward and 8 blocks "
          "agree to 2 ulp");
  }

  /*--------------------------------------------------------------------------
  compensatedAdd on two doubles is the step of compensatedSum.
  --------------------------------------------------------------------------*/
  {
    compensatedSum reference;
    double sum = 0, compensation = 0;
    const double terms[5] = {1.0, 1e-16, -3e-17, 1e-16, 2.5e-17};
    for (const double x : terms) {
      reference.add(x);
      compensatedAdd(sum, compensation, x);
    }
    check(sum == reference.sum() && compensation == reference.compensation(),
          "compensatedAdd is the step of compensatedSum");
  }

  /*--------------------------------------------------------------------------
  Release totals on the solver's clock (interfaceExchange::beginStep and
  balance): t_0 = start time, t_k+1 = t_k + dt in double precision, the
  release of a step n0 M overlap_k/length with length =
  releaseWindowLength(), accumulated with compensatedAdd (the ledger).  The
  total is n0 M to 4 eps; the overlaps sum to length to 2 ulp.  Controls:
  a plain running total, and the nominal duration instead of length.
  --------------------------------------------------------------------------*/
  {
    struct window { double t0, dt, start, duration, plainMin, nominalMin; };
    const double n0M = 4.84e-6*0.46100894;
    const double eps = 2.220446049250313e-16;
    /* plainMin/nominalMin: the error the control must at least show
       (0: no requirement) */
    const window windows[4] = {
      {0.0,    1e-3,    3e-4,      6.0,    1e-14, 0.0},    /* run/012 */
      {0.0,    2.2e-4,  3e-4,      15.0,   1e-13, 0.0},    /* M8 scale */
      {60.0,   2.2e-4,  60.0001,   0.0198, 0.0,   1e-13},
      {1000.0, 2.2e-4,  1000.0001, 0.0198, 0.0,   1e-12}
    };

    for (const window& w : windows) {
      const double length = releaseWindowLength(w.start, w.duration);
      const double end = w.start + w.duration;

      double released = 0, releasedC = 0;      /* ledger: compensated */
      double overlaps = 0, overlapsC = 0;
      double plain = 0, nominal = 0, nominalC = 0;
      long nSteps = 0;

      for (double t = w.t0; t < end + 2*w.dt; ++nSteps) {
        const double t1 = t + w.dt;
        const double overlap = releaseOverlap(t, t1, w.start, w.duration);
        compensatedAdd(released, releasedC, n0M*overlap/length);
        compensatedAdd(overlaps, overlapsC, overlap);
        compensatedAdd(nominal, nominalC, n0M*overlap/w.duration);
        plain += n0M*overlap/length;
        t = t1;
      }

      const double total = released + releasedC;
      const double error = std::abs(total/n0M - 1);
      const double errorPlain = std::abs(plain/n0M - 1);
      const double errorNominal = std::abs((nominal + nominalC)/n0M - 1);
      const double ulp = std::nextafter(length, 2*length) - length;

      std::printf("release window %.10g s + %.10g s, dt %.3g s, %ld steps: "
                  "compensated %.2e, plain %.2e, nominal duration %.2e "
                  "(relative to n0 M)\n", w.start, w.duration, w.dt, nSteps,
                  error, errorPlain, errorNominal);

      check(error <= 4*eps, "release total = n0 M to 4 eps (compensated, "
                            "window length as represented)");
      check(std::abs(overlaps + overlapsC - length) <= 2*ulp,
            "the overlaps of all steps sum to releaseWindowLength");
      check(errorPlain >= w.plainMin,
            "control: a plain running total misses n0 M");
      check(errorNominal >= w.nominalMin,
            "control: the nominal duration misses n0 M");
    }
  }

  std::printf("%d of %d phase-change kinetics checks passed.\n",
              nChecks - nFailed, nChecks);
  return nFailed == 0 ? 0 : 1;
}
