/*------------------------------------------------------------------------------

profileBinning.H

PURPOSE

Standalone arithmetic of the axial profiles of rhoFixedFlowFoam (plan
section 2, E10; axialProfiles.H), written as small C++ functions in the
style of phaseChangeKinetics.H.  They use no OpenFOAM types, so
tests/testProfiles.C checks them with plain g++.

Bins.  A bin set is a list of ascending edges e_0 < e_1 < ... < e_n along
the profile axis; bin b is [e_b, e_b+1] with the width dx_b = e_b+1 - e_b
and the centre x_b = (e_b + e_b+1)/2.  uniformBinEdges() builds them from
min, max and a number of bins, widthBinEdges() from min, max and a bin
width (the last bin is then shorter if the width does not divide the
range; a width at or above the range gives one bin).

Exact overlap binning.  An item (a cell, or the wall face of an interface
element) covers the axial interval [lo, hi] of its points.  Its amount is
distributed over the bins in proportion to the overlap,

  w_b = |[lo, hi] ^ [e_b, e_b+1]| / (hi - lo),

and the part of [lo, hi] outside [e_0, e_n] is returned as the weight
outside the bins, so that the weights plus the outside weight sum to 1
(to a few ulp: every piece is divided by the sum of the pieces, not by
hi - lo).  The sum over the bins of the binned amounts plus the amount
outside is therefore the total amount, whatever the mesh and the bins; a
profile in mol/m integrates to the inventory.  An item of zero axial extent
(e.g. a face normal to the axis) goes entirely into the bin that contains
its position (the last bin for e_n itself), or outside.

Onset of the deposit (T_dep).  Given the line density O_b of the observable
per bin, the onset is where O first reaches the threshold f*peak:

  fromInlet  search range: the bins from 'first' on (the bins whose centre
             lies at or after searchFrom); peak = the maximum of O over the
             range; the onset is the first bin b of the range with
             O_b >= f*peak;
  fromPeak   the peak bin p (the first maximum of the range, normally the
             whole set); walking upstream from p, the onset is the last bin
             b <= p whose upstream neighbour b-1 lies below f*peak, i.e.
             the upstream edge of the peak's own deposit.

Both bracket the crossing between the bins b-1 (O < f*peak) and b
(O >= f*peak); the position and any other quantity (the wall temperature)
are interpolated linearly between the two bin centres with the weight

  w = (f*peak - O_b-1)/(O_b - O_b-1),   0 < w <= 1.

The result is 'not bracketed' when the crossing bin is the first bin of the
range (the observable is at or above the threshold from the first bin on,
so the onset lies upstream of the range), 'no peak' when the peak is not
positive, and 'insignificant' when the amount in the peak bin, O_p dx_p
[mol], is at most a given minimum amount (axialProfiles.H: the significance
times the inventory of the pair): a column of solver noise (e.g. a deposit
of 1e-30 mol/m ahead of the gas front) has a peak and a crossing, but no
onset (review of M3, carry-over (b) of M5).

------------------------------------------------------------------------------*/

#ifndef profileBinning_H
#define profileBinning_H

#include <cmath>
#include <vector>

namespace LESTO {

  /*--------------------------------------------------------------------------
  Edges of nBins equal bins in [min, max] (nBins > 0, max > min); the last
  edge is max itself.
  --------------------------------------------------------------------------*/
  inline std::vector<double> uniformBinEdges (const double min,
                                              const double max,
                                              const int    nBins) {
    std::vector<double> edges(nBins + 1);
    for (int b = 0; b < nBins; ++b) {
      edges[b] = min + (max - min)*b/nBins;
    }
    edges[nBins] = max;
    return edges;
  }

  /*--------------------------------------------------------------------------
  Edges of bins of the given width from min (width > 0, max > min, and
  (max - min)/width within the range of int): n bins of the width when
  (max - min)/width is an integer n to 1e-6, otherwise the whole bins plus
  a shorter last bin that ends at max.  Always at least one bin: a width at
  or above the range gives the single bin [min, max].
  --------------------------------------------------------------------------*/
  inline std::vector<double> widthBinEdges (const double min,
                                            const double max,
                                            const double width) {
    const double ratio = (max - min)/width;
    const double whole = std::floor(ratio + 0.5);
    const int rounded = std::abs(ratio - whole) <= 1e-6*(whole > 1 ? whole : 1)
      ? int(whole)
      : int(std::ceil(ratio));
    const int nBins = rounded > 1 ? rounded : 1;
    std::vector<double> edges(nBins + 1);
    for (int b = 0; b < nBins; ++b) {
      edges[b] = min + b*width;
    }
    edges[nBins] = max;
    return edges;
  }

  /*--------------------------------------------------------------------------
  First bin whose upper edge lies above x (bins 0 ... nBins-1 of the edges
  e_0 ... e_nBins); nBins if x >= e_nBins.
  --------------------------------------------------------------------------*/
  inline int firstBinAbove (const double  x,
                            const double* edges,
                            const int     nBins) {
    int lo = 0, hi = nBins;
    while (lo < hi) {
      const int mid = (lo + hi)/2;
      if (edges[mid + 1] > x) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  /*--------------------------------------------------------------------------
  Exact overlap of the interval [lo, hi] (lo <= hi) with the bins.  Calls
  add(b, w) for every bin b with a non-zero weight w and returns the weight
  outside [e_0, e_nBins]; the weights and the outside weight sum to 1 to a
  few ulp.
  --------------------------------------------------------------------------*/
  template<class AddFunction>
  inline double intervalOverlap (const double  lo,
                                 const double  hi,
                                 const double* edges,
                                 const int     nBins,
                                 AddFunction   add) {

    const double first = edges[0];
    const double last  = edges[nBins];

    /* zero extent: the bin that holds the point, or outside */
    if (!(hi > lo)) {
      if (lo < first || lo > last) {
        return 1.0;
      }
      const int b = firstBinAbove(lo, edges, nBins);
      add(b < nBins ? b : nBins - 1, 1.0);
      return 0.0;
    }

    const double below = lo < first ? (hi < first ? hi : first) - lo : 0.0;
    const double above = hi > last ? hi - (lo > last ? lo : last) : 0.0;

    /* pieces inside, and their sum with the outside parts */
    double total = below + above;
    const int b0 = firstBinAbove(lo, edges, nBins);
    int b1 = b0;
    for (int b = b0; b < nBins && edges[b] < hi; ++b) {
      const double a = lo > edges[b] ? lo : edges[b];
      const double c = hi < edges[b + 1] ? hi : edges[b + 1];
      if (c > a) {
        total += c - a;
      }
      b1 = b + 1;
    }

    for (int b = b0; b < b1; ++b) {
      const double a = lo > edges[b] ? lo : edges[b];
      const double c = hi < edges[b + 1] ? hi : edges[b + 1];
      if (c > a) {
        add(b, (c - a)/total);
      }
    }

    return (below + above)/total;
  }

  /*--------------------------------------------------------------------------
  Onset of the deposit (see PURPOSE).  status: FOUND, NO_PEAK (the peak of
  the range is not positive, or the range is empty) or NOT_BRACKETED (the
  crossing bin is the first bin of the range).
  --------------------------------------------------------------------------*/
  enum onsetStatus {
    ONSET_FOUND = 0,
    ONSET_NO_PEAK = 1,
    ONSET_NOT_BRACKETED = 2,
    ONSET_INSIGNIFICANT = 3
  };

  struct onsetCrossing {
    int    status;
    int    peakBin;    /* first maximum of the range, -1 if it is empty */
    double peak;       /* O at peakBin */
    double threshold;  /* fraction*peak */
    int    lower;      /* the crossing lies between the centres of the bins */
    double weight;     /* lower and lower + 1, at this fraction (0, 1] */
  };

  /*--------------------------------------------------------------------------
  widths (optional, nBins values [m]) and minimumAmount [mol]: the peak is
  significant only when O at the peak bin times its width exceeds
  minimumAmount; without widths every positive peak is.
  --------------------------------------------------------------------------*/
  inline onsetCrossing onsetSearch (const double* O,
                                    const int     nBins,
                                    const int     first,
                                    const bool    fromPeak,
                                    const double  fraction,
                                    const double* widths = nullptr,
                                    const double  minimumAmount = 0.0) {

    onsetCrossing c = {ONSET_NO_PEAK, -1, 0.0, 0.0, -1, 0.0};

    for (int b = first; b < nBins; ++b) {
      if (c.peakBin < 0 || O[b] > c.peak) {
        c.peakBin = b;
        c.peak = O[b];
      }
    }
    if (c.peakBin < 0 || !(c.peak > 0.0)) {
      return c;
    }
    c.threshold = fraction*c.peak;

    if (widths && !(c.peak*widths[c.peakBin] > minimumAmount)) {
      c.status = ONSET_INSIGNIFICANT;
      return c;
    }

    int b = first;
    if (fromPeak) {
      b = c.peakBin;
      while (b > first && O[b - 1] >= c.threshold) {
        --b;
      }
    } else {
      while (O[b] < c.threshold) {
        ++b;
      }
    }

    if (b == first) {
      c.status = ONSET_NOT_BRACKETED;
      return c;
    }

    c.status = ONSET_FOUND;
    c.lower = b - 1;
    c.weight = (c.threshold - O[b - 1])/(O[b] - O[b - 1]);
    return c;
  }
}

#endif
