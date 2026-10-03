/*------------------------------------------------------------------------------

phaseChangeLedger.H

PURPOSE

Cumulative mass ledger of the gas/condensate pairs of rhoFixedFlowFoam (plan
section 2, E8).  It starts with one global entry,

  TIME           the exact time [s] at the end of the last step booked: the
                 clock of the release schedule (see below)

and holds for every pair, in kg:

  INITIAL        inventory at the fresh start (gas + wall + sample)
  RELEASED       mass injected by release samples
  REMOVED        condensate withdrawn with an inventory sample at its
                 removeTime (the reservoir sum Cs V, booked at the start of
                 the removal step; milestone M5)
  CLAMPED        round-off residue of the condensate: removed when a
                 CAPPED element is set to exactly zero, and the rounding
                 of every other condensate update (the stored condensate
                 changes by the rounded update, not by exactly dt q; plan
                 section 27, item 1); never a clip of Y
  SOLVER_DEFECT  sum of R*deltaT, R the per-cell defect of the solved gas
                 equation in flux form: the transport-only matrix (its face
                 fluxes, its face-flux correction included, and the storage
                 measured as the gas inventory is), plus the realised
                 exchange and the release (matrixDefect.H)
  RESTART        sum over restarts of (inventory re-read from the fields -
                 inventory held by the ledger): field rounding of an ascii
                 restart (plan section 8, risk 4)
  HELD           inventory at the end of the last step
  PROJECTED      diagnostic only: mass moved by the conservative same-cell
                 projection after the guard (not part of the closure)
  EVAPORATED_ABOVE_WARNING
                 diagnostic only: mass that entered the gas in cells whose
                 mole fraction exceeded moleFractionWarning
  REMOVAL_DONE   the number of inventory samples of the pair removed so
                 far (a count; which ones is recorded per sample in
                 uniform/phaseChangeLayout, and a restart checks that both
                 agree; milestone M5)

followed by one TRANSPORT entry per non-coupled, non-empty boundary patch
(outward positive), so that e.g. back-diffusion through the inlet is booked
separately from the outlet.  The closure of a pair is

  closure = held - supplied,
  supplied = INITIAL + RELEASED - sum_p TRANSPORT_p - REMOVED - CLAMPED
             - SOLVER_DEFECT + RESTART.

Accumulation.  The cumulative entries (RELEASED, TRANSPORT, SOLVER_DEFECT,
CLAMPED, REMOVED, RESTART and the diagnostics) grow by one increment per
time step.  A plain running total rounds every addition against the
growing total, so its error grows with the number of steps: 6.2e-14 of the
release after the 6100 steps of a 6 s window at dt = 1 ms, 5.2e-13 and a
closure of 1.3e-12 after 136364 steps of 0.22 ms (review M2, round 2).
Every entry therefore keeps its Neumaier compensation (compensatedAdd,
compensatedSum.H) next to it; add() and addTransport() update both, the
accessors return value + compensation, and supplied() adds all terms
compensated.  The totals then stay at a few eps whatever the number of
steps.  The closure is then the round-off of the inventory sums and of
one step's booking, and does not drift: SOLVER_DEFECT is the defect of the
transport-only matrix in flux form plus the realised exchange
(matrixDefect.H), so neither the stiff exchange coefficients nor the
double-precision imbalances of the assembled transport operator (its
column sums, the storage products, the source form of its face-flux
correction) enter the closure, and CLAMPED takes
the rounding of the stored condensate (review of M4, round 2: before, the
closure drifted linearly where the gas or a deposit stayed in place for
many steps, to 7.6e-13 of the inventory in 6000 steps).  TIME is not
compensated: it repeats the arithmetic of Time::operator++ (below).

Relative values (closure, solverDefect, restart) are printed relative to

  reference = max(|INITIAL + RELEASED|, |HELD|, inflow),
  inflow    = sum_p max(-TRANSPORT_p, 0),

so that a pair fed only through an inlet (no initial gas, no sample) is not
divided by zero.

Clock.  The release schedule needs the exact start and end of every step,
[t, t + deltaT].  The time-directory names are rounded to timePrecision
(e.g. 10.0007 for 10.00066 at timePrecision 6), and after a restart
OpenFOAM continues from the rounded name, so neither the names nor the
time value of a restarted run can be used.  TIME is set to the start time
on a fresh start, advanced by deltaT in every step (the same arithmetic as
Time::operator++ of a continuous run) and restored from the ledger on a
restart, so a continuous run and any binary restart use the same step
times to the last bit.  The balance files print it in the time column.

Patch lists in parallel (review item B1).  The transport entries are sized
by the non-coupled patches only.  In OpenFOAM every rank holds the same
non-processor patches first, in the same order, followed by its own
processor patches, whose number differs between ranks.  A list sized by
mesh.boundary().size() would therefore have a different length on every
rank, and a reduction over it fails (MPI_ERR_TRUNCATE).

Input and output.  The values are a list <time>/uniform/
phaseChangeBalance that OpenFOAM writes with the fields, always in BINARY
(binaryIOList.H), whatever the case's writeFormat, so that a restart
continues the balance with full precision.  It is registered AUTO_WRITE:
every write of Time stores it, also the writes that function objects start
between two steps (Time::writeNow(), Time::writeAndEnd(), e.g. runTimeControl
with satisfiedAction end); between two steps the ledger holds the balance at
the end of the last step, consistent with the fields.  It holds the running
values of all entries (TIME first, then the pairs) followed by their
compensations in the same order, so a binary restart continues the
compensated sums bit for bit.  dictionary::set() is deliberately not used
for the values: a primitiveEntry stores a scalar through a text round trip
at the default output precision.

Global objects.  In parallel every rank holds the same (reduced) values and
the same layout (below), so both are GLOBAL objects (binaryGlobalIOList,
IOdictionary): the master reads them and broadcasts them.  Uncollated,
every rank writes its identical copy into processor<i>/<time>/uniform/; with
-fileHandler collated the master writes one plain file into
processors<N>/<time>/uniform/, which reconstructPar and decomposePar copy
and which a serial run or a run on any number of ranks reads (review of M5,
round 3: as per-rank objects they were decomposedBlockData containers of N
blocks, and every restart after reconstructPar, or on another number of
ranks, stopped with a FatalIOError of the file handler).

Layout.  The ledger is addressed by position: pair i, entry k, transport
patch p.  Next to it, <time>/uniform/phaseChangeLayout (a global
dictionary, written with the fields as well) names what the positions mean:

  pairs             (PbI2_g);          the gas species of the pairs
  condensates       (PbI2_s);          their condensates, in the same order
  transportPatches  (INLET OUTLET WALL);  the transport entries

A restart compares them with the current pairs and patches and stops with
a FatalError that lists both when they differ (e.g. the pairs reordered in
thermochemistryProperties, or patches renamed or reordered with the same
count), instead of booking the entries to the wrong pair or patch.

The layout is also the record of what the restart state means beyond the
ledger: interfaceExchange adds the interface (the element areas of its
per-element deposits), the samples (the schedule of a release sample; the
pair, amount, cells, volume, removeTime and removal of an inventory
sample) and the profile times written so far, compares them on a restart
and stops when a change would create, destroy or misattribute mass
(interfaceExchange.H, "Restart record"), and the rounded 'uniform'
entries of the fields written with the state, merged over the ranks
(uniformValues).  Values that must compare exactly are stored as strings of
exactDecimal(): a scalar entry of a dictionary passes through text at the
output precision.

Balance files.  The master writes one line per step and pair to
postProcessing/phaseChangeBalance/<start time>/<gas species>.dat, with 17
significant digits; the time column is the clock TIME.

------------------------------------------------------------------------------*/

#ifndef phaseChangeLedger_H
#define phaseChangeLedger_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "PtrList.H"
#include "binaryIOList.H"

namespace LESTO {

class phaseChangeLedger {

public:

  /* global entries, before the pairs */
  enum globalEntry {
    TIME,
    nGlobalEntries
  };

  /* entries of every pair */
  enum entry {
    INITIAL,
    RELEASED,
    REMOVED,
    CLAMPED,
    SOLVER_DEFECT,
    RESTART,
    HELD,
    PROJECTED,
    EVAPORATED_ABOVE_WARNING,
    REMOVAL_DONE,
    nBaseEntries
  };

private:

  const Foam::fvMesh&               mesh_;
  Foam::wordList                    pairNames_;
  Foam::labelList                   patchIds_;  /* non-coupled, non-empty */
  Foam::label                       nEntries_;  /* entries, without their
                                                   compensations */
  binaryGlobalIOList<Foam::scalar>  values_;    /* entries, then
                                                   compensations */
  Foam::IOdictionary                layout_;    /* uniform/
                                                   phaseChangeLayout */
  bool                              restarted_;
  Foam::PtrList<Foam::OFstream>     files_;     /* master only */

  Foam::label index(const Foam::label pairi, const Foam::label k) const {
    return nGlobalEntries + (nBaseEntries + patchIds_.size())*pairi + k;
  }

  /* entry i: its running value plus its compensation */
  Foam::scalar value(const Foam::label i) const {
    return values_[i] + values_[nEntries_ + i];
  }

  /* entry i += x, compensated */
  void addTo(const Foam::label i, const Foam::scalar x);

  /* entry i := x, compensation 0 */
  void setTo(const Foam::label i, const Foam::scalar x) {
    values_[i] = x;
    values_[nEntries_ + i] = 0;
  }

public:

  /*--------------------------------------------------------------------------
  Reads <start time>/uniform/phaseChangeBalance if it exists (restarted()
  is then true) and checks its layout (uniform/phaseChangeLayout) against
  the pairs and patches of this run; otherwise all entries and
  compensations start at zero.  pairNames are the gas species of the
  pairs, condensateNames their condensates.
  --------------------------------------------------------------------------*/
  phaseChangeLedger (
    const Foam::fvMesh&   mesh,
    const Foam::wordList& pairNames,
    const Foam::wordList& condensateNames
  );

  bool restarted() const { return restarted_; }

  /* exact time at the end of the last step booked [s] (the clock) */
  Foam::scalar& time() { return values_[TIME]; }
  Foam::scalar time() const { return values_[TIME]; }

  /*--------------------------------------------------------------------------
  The shortest decimal representation of x (15 to 17 significant digits)
  that reads back exactly: 0.004 for 0.004, 0.10000000000000007 for the sum
  of one hundred steps of 0.001.  Equal strings mean equal values.
  --------------------------------------------------------------------------*/
  static Foam::word exactDecimal(const Foam::scalar x);

  /* uniform/phaseChangeLayout: read on a restart, written with the ledger;
     other classes add their own records to it (see Layout above) */
  Foam::dictionary& layout() { return layout_; }

  /* ids of the non-coupled, non-empty patches: the transport entries */
  const Foam::labelList& patchIds() const { return patchIds_; }

  /* entry e of pair pairi (value plus compensation) */
  Foam::scalar operator()(const Foam::label pairi, const entry e) const {
    return value(index(pairi, e));
  }

  /* entry e := x (e.g. INITIAL, HELD) */
  void set(const Foam::label pairi, const entry e, const Foam::scalar x) {
    setTo(index(pairi, e), x);
  }

  /* entry e += x, compensated (the cumulative entries) */
  void add(const Foam::label pairi, const entry e, const Foam::scalar x) {
    addTo(index(pairi, e), x);
  }

  /* transport through the k-th entry of patchIds() [kg, outward] */
  Foam::scalar transport(const Foam::label pairi, const Foam::label k) const {
    return value(index(pairi, nBaseEntries + k));
  }

  /* transport through the k-th entry of patchIds() += x, compensated */
  void addTransport (
    const Foam::label  pairi,
    const Foam::label  k,
    const Foam::scalar x
  ) {
    addTo(index(pairi, nBaseEntries + k), x);
  }

  /* sum of the transport entries [kg], compensated */
  Foam::scalar totalTransport(const Foam::label pairi) const;

  /* INITIAL + RELEASED - transport - REMOVED - CLAMPED - SOLVER_DEFECT
     + RESTART [kg], every value and compensation added compensated */
  Foam::scalar supplied(const Foam::label pairi) const;

  /* sum over the patches of max(-TRANSPORT, 0) [kg]: cumulative inflow */
  Foam::scalar inflow(const Foam::label pairi) const;

  /* max(|INITIAL + RELEASED|, |HELD|, inflow) [kg]: the reference of
     relative values */
  Foam::scalar reference(const Foam::label pairi) const;

  /*--------------------------------------------------------------------------
  Master only: opens postProcessing/phaseChangeBalance/<start>/<pair>.dat
  and writes the column header.
  --------------------------------------------------------------------------*/
  void openFiles();

  /* one line of the balance file of pair pairi (no-op on other ranks) */
  void writeLine (
    const Foam::label  pairi,
    const Foam::scalar gas,
    const Foam::scalar wall,
    const Foam::scalar sample,
    const Foam::scalar closure
  );
};

}

#endif
