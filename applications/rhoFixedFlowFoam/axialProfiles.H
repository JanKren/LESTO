/*------------------------------------------------------------------------------

axialProfiles.H

PURPOSE

Axial profiles of the gas/condensate pairs of rhoFixedFlowFoam (plan
section 2, E10; milestone M3): line densities in mol/m along an axis, with
exact overlap binning, the wall temperature per bin and the onset of the
deposit, T_dep.  They are written to postProcessing/axialProfiles only.

The profiles are the output that compares with the gamma scans of the
experiments (Lobresco et al. 2026, Figs. 9-14): the amount of PbI2 per unit
length of the tube, on the native axial cells and on the 1 cm bins of the
scans.

Dictionary (constant/thermochemistryProperties; the sub-dictionary is
optional, and without it nothing is computed, written or printed)

  profiles
  {
      axis     (1 0 0);          // direction of the profile axis
      origin   (0 0 0);          // x = (point - origin) & axis [m]
      bins                       // one or more bin sets, each written to
      {                          // its own file
          native { min 0; max 0.69; nBins 424; }  // = the axial cells
          gamma  { min 0; max 0.69; width 0.01; } // gamma-scan bins
      }                          // min/max default to the extent of the
                                 // mesh along the axis
      times      (7 12 16);      // [s], in addition to the write times
      onset
      {
          fraction    0.01;                  // of the peak (default 0.01)
          search      (fromInlet fromPeak);  // one or both (default both)
          searchFrom  0;                     // [m], fromInlet only
                                             // (default: the first bin)
          significance 1e-12;                // relative to the pair
      }                                      // inventory (default: balance
                                             // guardTolerance)
      observable wallPlusGas;    // wall | wallPlusGas (default)
  }

Line densities of bin b of width dx_b [mol/m], with the overlap weights w
of profileBinning.H (cells: from their points; interface elements of layer
faceCells: from the points of their wall face; of layer distance: from the
points of their cell, like the gas of the cell, plan section 25):

  gas      n'_g(b) = sum_c w_cb rho_c Y_c V_c/(M dx_b)
  wall     n'_w(b) = sum_e w_eb m''_e A_e/(M dx_b)    (A_e the effective
                                                        element area)
  sample   n'_s(b) = sum_c w_cb Cs_c V_c/(M dx_b)   (the reservoir of the
           inventory samples; zero for release samples, which inject
           directly into the gas)
  sampleGas  the gas line density of the cells of every sample (release
           and inventory), part of n'_g

so that sum_b n'(b) dx_b M plus the amount outside the bins is the gas,
wall and sample inventory of the pair.  The wall temperature of a bin is
the area-weighted temperature of the faces of the interface patches,
T_w(b) = sum_f w_fb |S_f| T_f/sum_f w_fb |S_f| (frozen: computed once), and
the wall area per length sum_f w_fb |S_f|/dx_b is written with it (layer
faceCells: these are the faces of its elements; layer distance: all faces
of the interface patches, so that T_dep is the same wall temperature for
both rules).

Onset (profileBinning.H, onsetSearch): the observable is n'_w (wall) or
wallPlusGas = n'_w + n'_g - n'_sampleGas: the wall deposit plus the gas
OUTSIDE the cells of every sample (the gas left in the tube at shutdown
condenses locally; the source cloud in the sample cells is not a deposit).
The sample reservoir n'_s is always excluded.  (Decided by the researcher
with the sign-off of M3, plan section 19, item 6: before, wallPlusGas
counted all gas, and while gas was flowing its peak and onsets lay on the
gas cloud of the sample, e.g. 949.5 K on the sample box of the test channel
at t = 0.03 s; at shutdown, with the sample removed, both definitions are
the same.)  T_dep and x_dep are interpolated linearly between the centres
of the two bins that bracket fraction*peak; fromInlet searches from
searchFrom downstream, fromPeak walks upstream from the peak.
utilities/axialProfiles.py prints the onset of the wall column next to that
of a wallPlusGas file.

Significance: a column whose peak bin holds at most significance times the
inventory of the pair (gas + wall + sample, mol) has no onset: its status
is 'insignificant' (a column of solver noise has a peak and a crossing but
no deposit).  significance defaults to the guardTolerance of the balance
dictionary (1e-12), the scale below which the exchange counts a mass as
round-off.

When the profiles are written.  At every write of Time, and at every entry
t_i of 'times'.
  - Writes of Time: the write times of controlDict and the writes that
    function objects start between two steps (Time::writeNow(),
    Time::writeAndEnd(), e.g. runTimeControl with satisfiedAction end or
    nWriteStep).  interfaceExchange hooks the profiles into Time's write
    (an AUTO_WRITE object whose writeObject() writes them), so every
    state directory has its profiles, with the reason 'write time'.
  - 'times': at the end of the first step whose end time t1 (the exact
    clock of the ledger) satisfies t1 >= t_i - deltaT/2, i.e. at the step
    end nearest to t_i.
Entries that fall on one step, or on the step of a write, give one set of
files (reason e.g. 'write time; profile time 7'; a function object that
writes after the files of a profile time rewrites them with that reason).
The entries reached so far are recorded with the state (uniform/
phaseChangeLayout, profiles { written (...); }, exact), so a restart skips
exactly the entries the run that wrote the restart state has written,
whatever the deltaT of either run (review of M3, carry-over (d) of M5:
before, the rule t_i - deltaT/2 <= start used the new run's deltaT, so a
restart with a smaller deltaT wrote an entry twice and one with a larger
deltaT lost one), and entries before the start that it did not write
(e.g. added to 'times' since); every other entry is written by the new run
at its first step end t1 >= t_i - deltaT/2.  On a fresh start (and on a
restart from a layout without the record, of an earlier build) an entry
whose nearest step end is the start state or earlier (t_i - deltaT/2 <=
start) is skipped: it lies before the first step (also an entry up to
deltaT/2 after the start): the initial state has profiles only if Time
writes it.  At a profile time that is not a write
time only postProcessing/axialProfiles/<time>/ is written: no field, no
uniform/ state, no time directory of the case (review item B4).

Files (master only): postProcessing/axialProfiles/<time>/<gas>_<set>.dat,
header lines '# key values' (the pair, the bin set, the exact time, the
inventories, the binned and outside totals, the onset of every search
mode), then one row per bin:

  x_lo x_hi x gas wall sample observable Twall wallArea sampleGas
  [m]  [m]  [m] [mol/m] [mol/m] [mol/m] [mol/m] [K] [m2/m] [mol/m]

with 17 significant digits, so that sum_b (x_hi - x_lo) n'(b) M reproduces
the inventories to round-off.  Twall and wallArea are 0 in a bin without
interface faces.  utilities/axialProfiles.py reads them.

Log: one line per pair and bin set at every profile time, with the totals
in the bins, the amounts outside them, the relative difference of binned
plus outside from the inventory (round-off), and T_dep of every search
mode.

Parallel.  The weights are local to each rank.  The binned amounts are
summed over the ranks with one Pstream::listCombineReduce per pair and bin
set, of a list whose size depends only on the dictionary (never on the
processor patches: review item B1), and every rank calls it at every
profile time, also a rank without interface elements.  The frozen wall
temperature and area per bin are reduced once at construction.

------------------------------------------------------------------------------*/

#ifndef axialProfiles_H
#define axialProfiles_H

#include "fvMesh.H"
#include "volFields.H"
#include "PtrList.H"
#include "FixedList.H"

namespace LESTO {

class axialProfiles {

  /*--------------------------------------------------------------------------
  One bin set: its edges, the overlap weights of the local cells and of the
  local interface elements (compressed rows: the entries of item i are
  start[i] ... start[i+1]-1), their weights outside the bins, and the
  frozen wall area and temperature per bin (all ranks).
  --------------------------------------------------------------------------*/
  class binSet {

  public:

    Foam::word        name;
    Foam::scalarField edges;           /* nBins + 1, ascending [m] */

    Foam::labelList   cellStart;
    Foam::labelList   cellBin;
    Foam::scalarField cellWeight;
    Foam::scalarField cellOutside;     /* per local cell */

    Foam::labelList   elementStart;
    Foam::labelList   elementBin;
    Foam::scalarField elementWeight;
    Foam::scalarField elementOutside;  /* per local element */

    Foam::scalarField wallArea;        /* sum_f w_fb |S_f| of the wall
                                          faces [m2] */
    Foam::scalarField wallT;           /* area-weighted T_f [K]; 0 without
                                          wall area */
    Foam::label       firstSearchBin;  /* first bin at or after searchFrom */

    Foam::label nBins() const { return edges.size() - 1; }
  };

  const Foam::fvMesh&    mesh_;
  Foam::vector           axis_;
  Foam::point            origin_;
  Foam::PtrList<binSet>  sets_;
  Foam::boolList         sampleCell_;  /* cells of every sample (excluded
                                          from wallPlusGas) */
  Foam::scalarList       times_;       /* ascending, distinct */
  Foam::boolList         done_;        /* entry written or skipped */
  Foam::DynamicList<Foam::scalar> written_;  /* entries written, by this
                                          run and the runs before it */
  Foam::scalar           fraction_;
  Foam::scalar           significance_;
  bool                   fromInlet_;
  bool                   fromPeak_;
  bool                   haveSearchFrom_;
  Foam::scalar           searchFrom_;
  bool                   wallPlusGas_;

  /* the exact overlap weights of intervals [lo_i, hi_i] with the bins */
  static void binIntervals (
    const Foam::scalarField& edges,
    const Foam::scalarField& lo,
    const Foam::scalarField& hi,
    Foam::labelList&         start,
    Foam::labelList&         bin,
    Foam::scalarField&       weight,
    Foam::scalarField&       outside
  );

  /* "T_dep <mode>" of one search, for the log and the file;
     minimumAmount [mol] is the significance of the peak bin */
  Foam::string onsetText (
    const binSet&            set,
    const Foam::scalarField& observable,
    const bool               fromPeak,
    const bool               forFile,
    const Foam::scalar       minimumAmount
  ) const;

public:

  /*--------------------------------------------------------------------------
  Reads the profiles dictionary, builds the bins and the overlap weights of
  the cells and of the interface elements (patch and face of the wall face
  of every wall element, or patch -1 and the element's cell: layer
  distance), the frozen wall temperature and area per bin (T on the wall
  faces wallPatch/wallFace), and skips the profile times already written
  or before the start (tStart, the clock of the ledger; written: the
  record of the entries written by the runs before a restart, nullptr on a
  fresh start or without a record).  sampleCells marks the cells of every
  sample (excluded from wallPlusGas); significance is the default of
  onset { significance; }.  Prints one start-up line.
  --------------------------------------------------------------------------*/
  axialProfiles (
    const Foam::fvMesh&         mesh,
    const Foam::dictionary&     dict,
    const Foam::labelList&      elementPatch,
    const Foam::labelList&      elementFace,
    const Foam::labelList&      elementCell,
    const Foam::labelList&      wallPatch,
    const Foam::labelList&      wallFace,
    const Foam::volScalarField& T,
    const Foam::boolList&       sampleCells,
    const Foam::scalar          significance,
    const Foam::scalar          tStart,
    const Foam::scalar          deltaT,
    const Foam::scalarList*     written,
    const bool                  restarted
  );

  axialProfiles(const axialProfiles&) = delete;
  void operator=(const axialProfiles&) = delete;

  /*--------------------------------------------------------------------------
  The 'times' entries reached by the state at stepEnd, the end of a step of
  deltaT (t_i - deltaT/2 <= stepEnd), as the text of the reason (e.g.
  "profile time 7; profile time 7.001"; empty: none).  Consumes them.  The
  same on every rank.
  --------------------------------------------------------------------------*/
  Foam::string reached (
    const Foam::scalar stepEnd,
    const Foam::scalar deltaT
  );

  /* the entries written so far (the restart record), ascending */
  Foam::scalarList written() const;

  /*--------------------------------------------------------------------------
  Bins, reduces, logs and writes the profiles of one pair, for every bin
  set.  Collective: every rank calls it.

    gasMass     per local cell [kg]      rho Y V of the gas
    wallMass    per local element [kg]   m''_e A_e
    sampleMass  per local cell [kg]      the sample reservoir (empty: 0)
    inventory   gas, wall, sample [kg], all ranks (the ledger's sums)
  --------------------------------------------------------------------------*/
  void write (
    const Foam::word&        timeName,
    const Foam::scalar       time,
    const Foam::string&      reason,
    const Foam::word&        gasName,
    const Foam::word&        condensateName,
    const Foam::scalar       molarMass,
    const Foam::scalarField& gasMass,
    const Foam::scalarField& wallMass,
    const Foam::scalarField& sampleMass,
    const Foam::FixedList<Foam::scalar, 3>& inventory
  ) const;
};

}

#endif
