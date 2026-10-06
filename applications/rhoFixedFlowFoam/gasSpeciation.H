/*------------------------------------------------------------------------------

gasSpeciation.H  (M10b gas-only equilibrium kernel)

PURPOSE

Homogeneous ideal-gas speciation of the trace elements Pb, Bi, I in an inert
carrier, at fixed T and p, for given element concentrations b_k [mol/m3]:
the re-speciation step of plan section 36 (option R1, the 'kernel').  It
solves the same minimisation as a GEMS3K call with every condensed phase
suppressed, for an ideal gas, in double precision, without the absolute
cut-offs of an interior-point method.

Model.  Every gas species s has p_s = exp(lnKf_s + sum_k nu_sk pi_k) [bar],
pi_k = ln p of the monatomic gas of element k, lnKf_s its formation
constant from the monatomic gases at the cell temperature (T is frozen:
precomputed once per cell).  Concentrations c_s = p_s 1e5/(R T) [mol/m3]
(the carrier adds the rest of p; the dilute gas does not change p).

Method.  No species contains both Pb and Bi (none in the data), so for a
given iodine potential pi_I the Pb and Bi balances decouple, each a
monotone 1-D problem in its own potential (Newton on the
logarithm, globally convergent from above with a bracket).  The iodine
balance with the metal balances enforced is strictly increasing in pi_I
(its derivative is a Schur complement of the positive definite Gram matrix
sum_s c_s nu_s nu_s^T), solved by a safeguarded Newton with that exact
derivative.  The API accepts a warm start; the solver adapter starts cold on every call
so no unsaved numerical state influences a restart.
A short extended-precision refinement removes cancellation in trace
mixtures. Exponentials use double libm with an extended argument correction.
Extremely dilute bulks bypass the double predictor, and extreme exponents
use extended libm so pressure underflow cannot degrade concentration accuracy.
The element balances hold to round-off by construction (the amounts are
the exponentials of the converged potentials; the residual is printed).

------------------------------------------------------------------------------*/

#ifndef gasSpeciation_H
#define gasSpeciation_H

#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace LESTO {

struct gasSpecies {
  int nPb, nBi, nI;      /* formula */
};

class gasSpeciation {

public:

  static constexpr int PB = 0, BI = 1, I = 2;

  explicit gasSpeciation(const std::vector<gasSpecies>& sp) : sp_(sp) {
    for (const auto& g : sp_) {
      if (g.nPb < 0 || g.nBi < 0 || g.nI < 0 || (g.nPb && g.nBi)
          || !(g.nPb || g.nBi || g.nI))
        throw std::invalid_argument("gasSpeciation: require nonnegative Pb/Bi/I formulas, without mixed metals");
    }
  }

  int size() const { return int(sp_.size()); }

  /*--------------------------------------------------------------------------
  lnKf: per species at the cell T; scale = 1e5/(R T) [mol/m3 per bar];
  b: element concentrations [mol/m3] (b_k <= 0: element absent);
  pi: in = warm start (NaN: cold), out = converged potentials;
  c: out, species concentrations [mol/m3].  Returns the number of outer
  iterations (negative: no convergence).
  --------------------------------------------------------------------------*/
  int solve (
    const double* lnKf,
    const double scale,
    const std::array<double, 3>& b,
    std::array<double, 3>& pi,
    double* c
  ) const {

    if (!(scale>0) || !std::isfinite(scale)) return -1;
    for (double amount:b) if (!std::isfinite(amount)) return -1;
    for (int j=0;j<size();++j) if (!std::isfinite(lnKf[j])) return -1;
    bool dilute = false;
    for (double amount : b)
      dilute = dilute || (amount > 0 && static_cast<long double>(amount) < 1e-280L*scale);
    std::array<double,3> fastPi = pi;
    const int iterations = dilute ? 0 : solveImpl<double>(lnKf, scale, b, fastPi, c);
    if (iterations < 0) return iterations;
    std::array<long double,3> precisePi = {fastPi[0], fastPi[1], fastPi[2]};
    const int refined = solveImpl<long double>(lnKf, scale, b, precisePi, c);
    pi = {double(precisePi[0]), double(precisePi[1]), double(precisePi[2])};
    return refined < 0 ? refined : (dilute ? refined : iterations);
  }

private:
  template<class Real> static Real exponential(Real x) {
    if (std::numeric_limits<Real>::digits > std::numeric_limits<double>::digits
        && (x < -700 || x > 700)) return std::exp(x);
    const double base = static_cast<double>(x);
    return static_cast<Real>(std::exp(base))*(1 + (x-static_cast<Real>(base)));
  }
  template<class Real>
  int solveImpl(const double* lnKf, double scale,
                const std::array<double,3>& b, std::array<Real,3>& pi,
                double* c) const {
    const Real NaN = std::numeric_limits<Real>::quiet_NaN();
    const bool hasPb = b[PB] > 0, hasBi = b[BI] > 0, hasI = b[I] > 0;

    /* inner solve for a metal E at given piI: returns pi_E and the
       derivative d pi_E/d pi_I (implicit function theorem) */
    auto metal = [&](int E, Real piI, Real start, Real& dpi) {
      Real x = std::isfinite(start) ? start : std::log(static_cast<Real>(b[E])/scale);
      Real lo = -1e300, hi = 1e300;
      for (int it = 0; it < 200; ++it) {
        Real s1 = 0, s2 = 0, sI = 0;
        for (int j = 0; j < size(); ++j) {
          const int nu = nuOf(j, E);
          if (!nu) continue;
          const int nI = sp_[j].nI;
          const Real e =
            exponential(Real(lnKf[j]) + nu*x + (nI && hasI ? nI*piI : 0.0))
            *(nI && !hasI ? 0.0 : 1.0);
          s1 += nu*e; s2 += nu*nu*e; sI += nu*nI*e;
        }
        const Real f = std::log(scale*s1) - std::log(static_cast<Real>(b[E]));
        if (f > 0) hi = std::min(hi, x); else lo = std::max(lo, x);
        Real dx = -f*s1/s2;
        if (dx > 30) dx = 30;
        if (dx < -30) dx = -30;
        Real xn = x + dx;
        if (!(xn > lo && xn < hi) && lo > -1e299 && hi < 1e299) {
          xn = 0.5*(lo + hi);
        }
        dpi = -sI/s2;
        if (std::fabs(f) < 2e-16L) break;      /* balance to round-off */
        const bool done = xn == x;
        x = xn;
        if (done) break;
      }
      return x;
    };

    auto gasI = [&](Real piI, Real piPb, Real piBi, Real dPb, Real dBi,
                    Real& deriv) {
      Real s = 0, d = 0;
      for (int j = 0; j < size(); ++j) {
        const int nI = sp_[j].nI;
        if (!nI) continue;
        const int nPb = sp_[j].nPb, nBi = sp_[j].nBi;
        if ((nPb && !hasPb) || (nBi && !hasBi)) continue;
        const Real e = exponential(Real(lnKf[j]) + nI*piI + (nPb ? nPb*piPb : 0.0)
                                  + (nBi ? nBi*piBi : 0.0));
        s += nI*e;
        d += nI*e*(nI + nPb*dPb + nBi*dBi);
      }
      deriv = d;
      return s;
    };

    Real piI = hasI ? (std::isfinite(pi[I]) ? pi[I] : std::log(static_cast<Real>(b[I])/scale)) : NaN;
    Real piPb = hasPb ? pi[PB] : NaN, piBi = hasBi ? pi[BI] : NaN;
    int outer = 0;

    if (!hasI) {
      Real d;
      if (hasPb) piPb = metal(PB, NaN, piPb, d);
      if (hasBi) piBi = metal(BI, NaN, piBi, d);
    } else {
      Real lo = -1e300, hi = 1e300;
      for (outer = 1; outer <= 200; ++outer) {
        Real dPb = 0, dBi = 0, deriv = 0;
        if (hasPb) piPb = metal(PB, piI, piPb, dPb);
        if (hasBi) piBi = metal(BI, piI, piBi, dBi);
        const Real s = gasI(piI, piPb, piBi, dPb, dBi, deriv);
        const Real f = std::log(scale*s) - std::log(static_cast<Real>(b[I]));
        if (f > 0) hi = std::min(hi, piI); else lo = std::max(lo, piI);
        Real dx = -f*s/deriv;
        if (dx > 30) dx = 30;
        if (dx < -30) dx = -30;
        Real xn = piI + dx;
        if (!(xn > lo && xn < hi) && lo > -1e299 && hi < 1e299) {
          xn = 0.5*(lo + hi);
        }
        const bool done = std::fabs(f) < 4e-16L || xn == piI;
        if (!done) piI = xn;
        if (done) {
          if (hasPb) piPb = metal(PB, piI, piPb, dPb);
          if (hasBi) piBi = metal(BI, piI, piBi, dBi);
          break;
        }
      }
      if (outer > 200) outer = -200;
    }

    pi = {Real(piPb), Real(piBi), Real(piI)};
    for (int j = 0; j < size(); ++j) {
      const gasSpecies& g = sp_[j];
      if ((g.nPb && !hasPb) || (g.nBi && !hasBi) || (g.nI && !hasI)) {
        c[j] = 0;
        continue;
      }
      c[j] = scale*exponential(Real(lnKf[j]) + (g.nPb ? g.nPb*piPb : 0.0)
                            + (g.nBi ? g.nBi*piBi : 0.0)
                            + (g.nI ? g.nI*piI : 0.0));
    }
    return outer;
  }


  int nuOf(int j, int E) const {
    return E == PB ? sp_[j].nPb : E == BI ? sp_[j].nBi : sp_[j].nI;
  }

  std::vector<gasSpecies> sp_;
};

}

#endif
