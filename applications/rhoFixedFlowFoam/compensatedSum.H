/*------------------------------------------------------------------------------

compensatedSum.H

PURPOSE

Compensated (Kahan-Babuska-Neumaier) summation for the inventories and the
cumulative entries of the phase-change ledger of rhoFixedFlowFoam.  It uses
no OpenFOAM types, so tests/testPhaseChange.C checks it with plain g++.

Why.  A plain sum of n terms carries a rounding error that depends on the
order of the terms, up to about n*eps of the sum (eps = 2.2e-16).  The
gas inventory of a pair, sum_c rho_c Y_c V_c, has one term per cell, and a
restart after decomposePar, reconstructPar or renumberMesh sums the same
terms in another order.  On the 228960-cell pipe the forward and the
reverse order already differ by 1.7e-14 of the inventory, above the 1e-14
bound with which a binary restart is compared with the ledger (plan
section 12, non-blocking item 1).

A compensated sum carries the rounding error of every addition in a second
variable.  Its error is about eps of the sum (for terms of one sign, and in
general eps times the sum of the magnitudes, which is the same here: the
negative Y of a paired species is solver noise), whatever the number and
the order of the terms.  Two compensated sums of the same terms in any two
orders therefore agree to a few eps.

Usage

  LESTO::compensatedSum s;
  for (...) { s.add(x); }
  s.value()          // the sum, rounded once

  In parallel, every rank sums its own terms; the pairs (sum(),
  compensation()) of all ranks are then gathered and added, in rank order,
  into one more compensatedSum (interfaceExchange.C, globalSum), so that
  the reduction over the ranks is compensated too and every rank holds the
  same value.

  compensatedAdd(sum, compensation, x) is the same step on two doubles kept
  elsewhere, e.g. the cumulative entries of the phase-change ledger
  (phaseChangeLedger.H): a running total of one increment per time step
  would otherwise carry a rounding error that grows with the number of
  steps (6e-14 of the release after the 6100 steps of run/012, 5e-13 after
  136364 steps of 0.22 ms), while the compensated total stays at a few eps
  whatever the number of steps.

Reference: A. Neumaier, Rundungsfehleranalyse einiger Verfahren zur
Summation endlicher Summen, ZAMM 54 (1974) 39-51.  The compiler must not
reassociate floating-point operations (no -ffast-math), which holds for the
OpenFOAM wmake rules.

------------------------------------------------------------------------------*/

#ifndef compensatedSum_H
#define compensatedSum_H

#include <cmath>

namespace LESTO {

  /*--------------------------------------------------------------------------
  One Neumaier step: sum += x, with the rounding error of the addition,
  recovered exactly (two-sum of the larger and the smaller magnitude),
  accumulated in compensation.  The compensated value is
  sum + compensation.
  --------------------------------------------------------------------------*/
  inline void compensatedAdd (double&      sum,
                              double&      compensation,
                              const double x) {
    const double t = sum + x;
    if (std::fabs(sum) >= std::fabs(x)) {
      compensation += (sum - t) + x;
    } else {
      compensation += (x - t) + sum;
    }
    sum = t;
  }

  class compensatedSum {

    double sum_;             /* running sum */
    double compensation_;    /* accumulated rounding errors of sum_ */

  public:

    compensatedSum() : sum_(0.0), compensation_(0.0) {}

    /* adds x (compensatedAdd) */
    void add(const double x) {
      compensatedAdd(sum_, compensation_, x);
    }

    /* the running sum and its compensation (for a reduction over ranks) */
    double sum() const { return sum_; }
    double compensation() const { return compensation_; }

    /* the compensated sum */
    double value() const { return sum_ + compensation_; }
  };
}

#endif
