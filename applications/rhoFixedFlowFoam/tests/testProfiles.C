/*------------------------------------------------------------------------------

tests/testProfiles.C

PURPOSE

Standalone test of profileBinning.H (no OpenFOAM libraries needed): the bin
edges from a number of bins or a width, the exact overlap weights of an
interval with the bins (inside, across several bins, partly or entirely
outside, of zero extent, and for random intervals and bins, where the
binned amounts plus the amount outside must reproduce the total), and the
onset search of T_dep (plan section 2, E10): a linear ramp gives the
analytic crossing in both search modes, an upstream spike separates
fromInlet from fromPeak, searchFrom skips it, and the edge cases (no peak,
the threshold reached from the first bin, fraction 1, a bin exactly at the
threshold) are reported as specified, and a column of noise whose peak bin
holds at most the minimum amount is 'insignificant' (M5, carry-over (b)).

Build and run from the solver directory:

  g++ -std=c++17 -Wall -Wextra tests/testProfiles.C -o /tmp/testProfiles
  /tmp/testProfiles

Criterion M3.0 of tests/Alltest M3 (and M5.7a of tests/Alltest M5).

------------------------------------------------------------------------------*/

#include "../profileBinning.H"
#include <cstdio>
#include <cmath>
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

  /* weights of one interval as a dense vector, and the outside weight */
  std::vector<double> weights (const double lo,
                               const double hi,
                               const std::vector<double>& edges,
                               double& outside) {
    const int n = int(edges.size()) - 1;
    std::vector<double> w(n, 0.0);
    outside = LESTO::intervalOverlap (
      lo, hi, edges.data(), n,
      [&](const int b, const double x) { w[b] += x; }
    );
    return w;
  }

  double sum(const std::vector<double>& v) {
    double s = 0;
    for (const double x : v) {
      s += x;
    }
    return s;
  }
}

int main() {

  using namespace LESTO;

  /*--------------------------------------------------------------------------
  Bin edges
  --------------------------------------------------------------------------*/
  {
    const std::vector<double> e = uniformBinEdges(0, 0.69, 424);
    bool ascending = true;
    for (std::size_t k = 1; k < e.size(); ++k) {
      ascending = ascending && e[k] > e[k - 1];
    }
    check(e.size() == 425 && e.front() == 0 && e.back() == 0.69 && ascending,
          "424 uniform bins: 425 ascending edges from 0 to exactly 0.69");
    check(close(e[1], 0.69/424, 1e-15), "uniform bin width 0.69/424");

    const std::vector<double> g = widthBinEdges(0, 0.69, 0.01);
    check(g.size() == 70 && g.back() == 0.69
          && close(g[69] - g[68], 0.01, 1e-12),
          "width 0.01 over [0, 0.69]: 69 bins (0.69/0.01 = 68.99999999999999)");

    const std::vector<double> h = widthBinEdges(0, 0.695, 0.01);
    check(h.size() == 71 && h.back() == 0.695
          && close(h[70] - h[69], 0.005, 1e-12),
          "width 0.01 over [0, 0.695]: 70 bins, the last one 5 mm");

    const std::vector<double> one = widthBinEdges(0.1, 0.2, 1.0);
    check(one.size() == 2 && one[0] == 0.1 && one[1] == 0.2,
          "a width larger than the range: one bin");

    /* a width of 1e6 range and more (ratio within the 1e-6 tolerance of
       0) gave no bin at all before (review of M3, round 1) */
    bool single = true;
    for (const double width : {1e5, 2e5, 1e6, 1e12}) {
      const std::vector<double> w = widthBinEdges(0, 0.1, width);
      single = single && w.size() == 2 && w[0] == 0 && w[1] == 0.1;
    }
    check(single, "a width of 1e5 to 1e12 over [0, 0.1]: one bin [0, 0.1]");

    const std::vector<double> exact = widthBinEdges(0, 0.1, 0.1);
    check(exact.size() == 2 && exact[1] == 0.1,
          "a width equal to the range: one bin");
  }

  /*--------------------------------------------------------------------------
  Overlap weights of single intervals
  --------------------------------------------------------------------------*/
  {
    const std::vector<double> e = {0.0, 1.0, 2.0, 3.0, 4.0};
    double out = -1;

    std::vector<double> w = weights(1.25, 1.75, e, out);
    check(w[1] == 1.0 && sum(w) == 1.0 && out == 0.0,
          "interval inside bin 1: weight 1");

    w = weights(0.5, 3.0, e, out);
    check(close(w[0], 0.2, 1e-15) && close(w[1], 0.4, 1e-15)
          && close(w[2], 0.4, 1e-15) && w[3] == 0.0 && out == 0.0,
          "interval [0.5, 3] over bins 0-2: weights 0.2 0.4 0.4");

    w = weights(-1.0, 1.0, e, out);
    check(close(w[0], 0.5, 1e-15) && close(out, 0.5, 1e-15),
          "interval [-1, 1]: half in bin 0, half outside");

    w = weights(3.5, 6.0, e, out);
    check(close(w[3], 0.2, 1e-15) && close(out, 0.8, 1e-15),
          "interval [3.5, 6]: 0.2 in bin 3, 0.8 outside");

    w = weights(-2.0, 6.0, e, out);
    check(close(sum(w), 0.5, 1e-15) && close(out, 0.5, 1e-15)
          && close(w[2], 0.125, 1e-15),
          "interval [-2, 6] around all bins: 0.125 per bin, 0.5 outside");

    w = weights(5.0, 6.0, e, out);
    check(sum(w) == 0.0 && out == 1.0, "interval entirely above: outside 1");

    w = weights(2.0, 2.0, e, out);
    check(w[2] == 1.0 && out == 0.0, "zero extent at the edge 2: bin 2");

    w = weights(4.0, 4.0, e, out);
    check(w[3] == 1.0 && out == 0.0, "zero extent at the last edge: last bin");

    w = weights(-0.5, -0.5, e, out);
    check(sum(w) == 0.0 && out == 1.0, "zero extent below the bins: outside");

    w = weights(1.0, 2.0, e, out);
    check(w[1] == 1.0 && w[0] == 0.0 && w[2] == 0.0,
          "interval equal to bin 1: no weight in its neighbours");
  }

  /*--------------------------------------------------------------------------
  Random intervals and bins: weights = overlap/length, and conservation of
  amounts binned plus outside
  --------------------------------------------------------------------------*/
  {
    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<double> u(0.0, 1.0);

    std::vector<double> edges(1, -0.1 + 0.05*u(rng));
    for (int k = 0; k < 57; ++k) {
      edges.push_back(edges.back() + 1e-3 + 0.05*u(rng));
    }
    const int n = int(edges.size()) - 1;
    const double span = edges.back() - edges.front();

    double maxWeightError = 0, maxSumError = 0;
    double total = 0, binned = 0, outside = 0;
    std::vector<double> bins(n, 0.0);

    for (int i = 0; i < 20000; ++i) {
      const double lo = edges.front() - 0.2*span + 1.4*span*u(rng);
      const double hi = lo + (i % 10 == 0 ? 0.0 : 0.3*span*u(rng));
      double out = 0;
      const std::vector<double> w = weights(lo, hi, edges, out);

      maxSumError = std::max(maxSumError, std::abs(sum(w) + out - 1.0));
      if (hi > lo) {
        for (int b = 0; b < n; ++b) {
          const double a = std::max(lo, edges[b]);
          const double c = std::min(hi, edges[b + 1]);
          const double exact = c > a ? (c - a)/(hi - lo) : 0.0;
          maxWeightError = std::max(maxWeightError, std::abs(w[b] - exact));
        }
      }

      const double amount = 1e-6*(u(rng) - 0.1);
      total += amount;
      outside += out*amount;
      for (int b = 0; b < n; ++b) {
        bins[b] += w[b]*amount;
      }
    }
    for (const double x : bins) {
      binned += x;
    }

    std::printf("random overlaps: max |sum w + outside - 1| %.3g, max weight"
                " error %.3g, (binned + outside - total)/total %.3g\n",
                maxSumError, maxWeightError,
                (binned + outside - total)/total);
    check(maxSumError <= 1e-15,
          "random intervals: weights plus outside sum to 1 to 1e-15");
    check(maxWeightError <= 1e-13,
          "random intervals: every weight is the overlap/length to 1e-13");
    check(std::abs(binned + outside - total) <= 1e-13*std::abs(total),
          "random amounts: binned plus outside reproduce the total to 1e-13");
  }

  /*--------------------------------------------------------------------------
  Onset: a linear ramp from x0 (a bin centre) to its peak, zero before
  --------------------------------------------------------------------------*/
  {
    const int n = 100;
    const double dx = 1e-3;
    std::vector<double> x(n), O(n, 0.0);
    const int k0 = 30, kp = 80;
    const double x0 = (k0 + 0.5)*dx, slope = 7.0;
    for (int b = 0; b < n; ++b) {
      x[b] = (b + 0.5)*dx;
      if (b >= k0 && b <= kp) {
        O[b] = slope*(x[b] - x0);
      }
    }
    const double f = 0.01;
    const double exact = x0 + f*(x[kp] - x0);

    for (const bool fromPeak : {false, true}) {
      const onsetCrossing c = onsetSearch(O.data(), n, 0, fromPeak, f);
      const double xc = x[c.lower] + c.weight*(x[c.lower + 1] - x[c.lower]);
      std::printf("ramp %s: x_dep %.17g, analytic %.17g\n",
                  fromPeak ? "fromPeak" : "fromInlet", xc, exact);
      check(c.status == ONSET_FOUND && c.peakBin == kp && c.lower == k0
            && std::abs(xc - exact) <= 1e-15,
            fromPeak ? "ramp, fromPeak: the analytic crossing to 1e-15 m"
                     : "ramp, fromInlet: the analytic crossing to 1e-15 m");
    }

    /* an upstream spike of 0.3 of the peak at bin 10 */
    std::vector<double> S(O);
    S[10] = 0.3*O[kp];
    const onsetCrossing inlet = onsetSearch(S.data(), n, 0, false, f);
    const onsetCrossing peak = onsetSearch(S.data(), n, 0, true, f);
    const double xInlet =
      x[inlet.lower] + inlet.weight*(x[inlet.lower + 1] - x[inlet.lower]);
    const double xPeak =
      x[peak.lower] + peak.weight*(x[peak.lower + 1] - x[peak.lower]);
    const double exactSpike = x[9] + f*O[kp]/S[10]*dx;
    check(inlet.status == ONSET_FOUND && inlet.lower == 9
          && std::abs(xInlet - exactSpike) <= 1e-15,
          "spike, fromInlet: the crossing below the spike");
    check(peak.status == ONSET_FOUND && std::abs(xPeak - exact) <= 1e-15,
          "spike, fromPeak: the crossing of the ramp");

    /* searchFrom past the spike */
    const onsetCrossing from = onsetSearch(S.data(), n, 20, false, f);
    const double xFrom =
      x[from.lower] + from.weight*(x[from.lower + 1] - x[from.lower]);
    check(from.status == ONSET_FOUND && std::abs(xFrom - exact) <= 1e-15,
          "spike, fromInlet from bin 20: the crossing of the ramp");

    /* the spike higher than the ramp: fromPeak finds the spike */
    std::vector<double> H(O);
    H[10] = 2.0*O[kp];
    const onsetCrossing high = onsetSearch(H.data(), n, 0, true, f);
    check(high.status == ONSET_FOUND && high.peakBin == 10 && high.lower == 9,
          "a spike above the ramp: fromPeak walks up from the spike");

    /* the range from the ramp's peak on: nothing below the threshold
       upstream of the peak within the range */
    const onsetCrossing late = onsetSearch(O.data(), n, kp, false, f);
    check(late.status == ONSET_NOT_BRACKETED,
          "search range starting at the peak: not bracketed");

    const onsetCrossing latePeak = onsetSearch(O.data(), n, kp, true, f);
    check(latePeak.status == ONSET_NOT_BRACKETED,
          "fromPeak, the peak at the first bin of the range: not bracketed");

    /* fraction 1: the crossing at the peak itself (weight 1) */
    const onsetCrossing full = onsetSearch(O.data(), n, 0, false, 1.0);
    check(full.status == ONSET_FOUND && full.lower == kp - 1
          && full.weight == 1.0, "fraction 1: the onset at the peak bin");

    /* a bin exactly at the threshold counts as reached */
    std::vector<double> T(n, 0.0);
    T[50] = 1.0;
    T[49] = 0.01;
    T[48] = 0.0;
    const onsetCrossing tie = onsetSearch(T.data(), n, 0, true, 0.01);
    check(tie.status == ONSET_FOUND && tie.lower == 48 && tie.weight == 1.0,
          "a bin at exactly the threshold is reached (fromPeak)");
    const onsetCrossing tieInlet = onsetSearch(T.data(), n, 0, false, 0.01);
    check(tieInlet.status == ONSET_FOUND && tieInlet.lower == 48,
          "a bin at exactly the threshold is reached (fromInlet)");

    /* no positive peak */
    std::vector<double> Z(n, 0.0);
    Z[3] = -1.0;
    check(onsetSearch(Z.data(), n, 0, false, f).status == ONSET_NO_PEAK
          && onsetSearch(Z.data(), n, 0, true, f).status == ONSET_NO_PEAK,
          "no positive value: no peak in both modes");
    check(onsetSearch(O.data(), n, n, false, f).status == ONSET_NO_PEAK,
          "an empty search range: no peak");

    /* significance (M5 carry-over (b)): a noise-level column whose peak
       bin holds at most the minimum amount is insignificant; the same
       column above it, and any column without widths, is found */
    std::vector<double> noise(n), widths(n, 0.002);
    for (int b = 0; b < n; ++b) {
      noise[b] = 1e-30*(1 + ((b*7919) % 13));
    }
    const onsetCrossing quiet =
      onsetSearch(noise.data(), n, 0, false, f, widths.data(), 1e-28);
    check(quiet.status == ONSET_INSIGNIFICANT && quiet.peakBin >= 0
       && quiet.peak == 13e-30 && quiet.threshold == f*13e-30,
          "noise column below the minimum amount: insignificant, with its "
          "peak and threshold");
    check(onsetSearch(noise.data(), n, 0, true, f, widths.data(), 1e-28)
            .status == ONSET_INSIGNIFICANT,
          "noise column: insignificant fromPeak as well");
    check(onsetSearch(noise.data(), n, 0, false, f, widths.data(), 2e-32)
            .status != ONSET_INSIGNIFICANT,
          "the same column above the minimum amount: searched");
    check(onsetSearch(O.data(), n, 0, false, f, widths.data(), 1e-28).status
            == ONSET_FOUND
       && onsetSearch(noise.data(), n, 0, false, f).status
            != ONSET_INSIGNIFICANT,
          "a significant ramp is found; without widths no significance test");
    check(onsetSearch(Z.data(), n, 0, false, f, widths.data(), 1e-28).status
            == ONSET_NO_PEAK,
          "no positive value: no peak before the significance test");
  }

  std::printf("%d of %d checks passed.\n", nChecks - nFailed, nChecks);
  return nFailed == 0 ? 0 : 1;
}
