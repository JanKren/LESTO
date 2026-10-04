/*------------------------------------------------------------------------------

gemsEquilibrium.H

PURPOSE

Optional GEMS3K backend of the Hertz-Knudsen-Schrage law of rhoFixedFlowFoam
(milestone M9 of doc/phase-change-plan.md: plan E9, sections 4, 8 risks 7
and 13, 9 items 7-8, 12 item 7 and 32).  The equilibrium vapour pressure
p_eq of a gas/condensate pair comes from a Gibbs-energy minimisation of the
gas with GEMS3K instead of the pair's vapourPressure table.  The table
stays required: it validates every GEMS3K result and replaces it wherever
a guard applies.

Build.  The GEMS3K bridge (thermochemistry/gemsbridge, a plain-C ABI in
libgemsbridge.so) is an OPTIONAL build dependency:

  - it is compiled in only when LESTO_GEMSBRIDGE names a directory holding
    gemsbridge.h and libgemsbridge.so at build time (see Make/options);
  - Make/gemsConfig.sh records that decision in the generated header
    gemsConfig.H, which defines LESTO_HAVE_GEMS or not.  gemsEquilibrium.C
    includes it, so a toggle of LESTO_GEMSBRIDGE recompiles it and relinks
    the executable even with a plain wmake;
  - this header never includes gemsbridge.h.  Only gemsEquilibrium.C does,
    and only under #ifdef LESTO_HAVE_GEMS; the engine and every call of
    the bridge live there.

available(), bridgeDirectory() and buildState() report how the executable
was built (rhoFixedFlowFoam -help).  requireBridge() stops a run that
selects equilibrium GEMS in a build without the bridge with a FatalIOError
and the rebuild hint; in a bridge build it returns.  The class itself is
compiled in both builds but constructed only in a bridge build.

How p_eq is obtained (plan E9).  Every condensed species of the GEMS3K
system is suppressed (metastability bounds 0 and 0), so a call returns the
homogeneous gas equilibrium of the bulk composition, and the dual
thermodynamic activity of a suppressed pure condensate is its log10
saturation ratio log10 Omega (gemsbridge.h).  For a pair whose gas species
is g and whose condensates are phi,

  p_eq = p_g 10^(-max_phi log10 Omega_phi),

the pressure of g over the most stable condensate.  For an ideal gas with
pure condensates (the PbI2-He system) p_eq depends on the temperature only
and equals the minimum over the phases of the table: GEMS3K reproduces
the table to its G0 interpolation (|dlog10 p| about 1e-5 to 6e-5 at
500-1250 K in cold calls).

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

gemsEquilibrium.H, continued (the description is split into several
comment blocks: wmkdepend of OpenFOAM v2412 stops scanning a file at a
comment longer than 16 KB; plan section 18, item 13).

Dictionary (constant/thermochemistryProperties; model HKS)

  HKSCoeffs  { equilibrium GEMS; speciation none; }   // or lagged

  GEMSCoeffs
  {
      system            "<constant>/gems/PbI2He-dat.lst";  // GEMS3K input set
      mode              local;      // frozen | local
      carrier           "He(g)";    // the carrier: a gas species of the system
      referenceComposition { He 1; Pb 1e-3; I 2.000002e-3; } // frozen [mol]
      referencePressure 101325;     // frozen [Pa]
      minMoleFraction   1e-6;       // local: below it the table
      minTemperature    500;        // [K] both modes: below it the table
      maxLog10Deviation 0.01;       // reject a result off the table by more
      updateInterval    1;          // local: steps between updates
      allowNonPhysicalCarrier no;   // local: see below
      logDirectory      "gemsLog";  // per rank: processorN/gemsLog
      logLevel          4;          // GEMS3K logging: 4 errors, 6 off
      writeEquilibriumField no;     // pEq_<gas>, gemsSource_<gas> (chi_<gas>)
  }

  pairs { PbI2_g { condensed PbI2_s; vapourPressure { ... }
    gems
    {
        gas          "PbI2(g)";              // the gas species of the pair
        condensates  ("PbI2(cr)" "PbI2(l)"); // its condensates (p_eq: the
                                             // most stable of them)
        elements     { Pb 1; I 2; }          // its formula: checked against
                                             // the system, also for every
                                             // condensate (one formula)
        excess       { I 1e-6; }             // optional: relative excess of
                                             // an element in the bulk
    }
  } }

  system, mode and carrier are required, and referenceComposition in mode
  frozen; every other key has the default shown (the values of plan
  section 4), and the start-up line prints the values used.  The system
  file names its DCH, IPM and DBR files relative to its own directory;
  <constant> and the like are expanded.  A key of the other mode
  (referenceComposition and referencePressure in local; minMoleFraction,
  updateInterval and allowNonPhysicalCarrier in frozen) is named in a
  warning and not read: it is neither validated nor can it stop the run.
  A key that is neither, e.g. a misspelled minMolFraction, stops the run
  with the list of the valid keys, in GEMSCoeffs and in a gems block
  (plan section 35, items G and H).

The excess.  At the exact formula PbI2 (b_I = 2 b_Pb) the interior-point
method of GEMS3K occasionally stalls (thermochemistry/README.md, pitfalls;
the scan of plan section 32: 1 failed warm call of 2364 at x >= 1e-6 and
T >= 500 K, none with a relative excess of 1e-6 iodine, and more failures
outside the guards).  excess adds b_i excess_i to the bulk of every
formula unit: b = c (nu_i (1 + excess_i)).  It changes p_eq of an ideal
gas not at all and the speciation by about the excess; the frozen
referenceComposition of the plan carries the same 1e-6 (I = 2.000002e-3).

mode frozen (plan section 12, item 7).  At start-up every rank creates an
engine; the master evaluates p_eq at the nodes of the pair's table that
lie in [max(minTemperature, T_min), T_max] of the GEMS3K grid, cold, with
the referenceComposition at referencePressure.  A node whose call fails, or
whose p_eq lies further than maxLog10Deviation from the table, takes the
table's value (counted).  A call that returns GEMSB_ERR_FATAL re-creates
the engine as mode local does (destroyed, created again from the system
file, every condensed species suppressed again; the calls of mode frozen
are cold, so no warm state survives it), and the start-up line counts the
re-creations (plan section 35, item E).  The values form a table of one
column, the 'GEMS table', interpolated as the pair's table (logInverseT
or linearT) and broadcast to every rank; the engines are then destroyed.
An element whose T_e lies within the GEMS table takes p_eq from it, every
other element from the pair's table (below minTemperature, outside the
grid): the same guards and fallback as local.  frozen is a consistency
check of the data (plan section 9, item 7) and costs nothing per step.

mode local.  One engine per rank for the whole run.  At every step whose
time index n (Time::timeIndex(), restored from uniform/time on a restart)
has (n - 1) mod updateInterval = 0, and at the first step of any run that
has no lagged state (below), beginStep() evaluates, before the
correctors, every element of every pair, a sample element for its own
pair only.  The schedule is that of the absolute time index: a run that
starts at index 0 updates at its first step and then every updateInterval
steps; a run that starts from a time directory whose uniform/time holds an
index k without a lagged state updates at its first step (k + 1) and then
at the next n of the schedule, so its first interval may be shorter (plan
section 35, item Y).  An element is evaluated at its temperature T_e,
the pressure p of its cell and the composition of its cell at the start
of the step, the carrier c_car = p/(R T) [mol/m3] (the frozen p and T, as
the mole-fraction monitor) times the carrier's formula plus c_k = rho
Y_k/M_k of every pair k times its formula (with its excess).  One call
serves every pair of the element.  A pair takes the GEMS3K p_eq unless a
guard applies, in this order, and then the table's (the 'source' of the
element):

  BELOW_TEMPERATURE    T_e < minTemperature (rc = 3 at 300-490 K, trace
                       Pb at 300 K does not converge; plan risk 7)
  OUTSIDE_GRID         T_e or p outside the T/P grid of the system
  BELOW_MOLE_FRACTION  x_k = c_k/(c_k + c_car) < minMoleFraction (silent
                       errors of several decades at x <= 1e-8)
  FAILED               the call returns another status than GEMSB_OK or
                       GEMSB_OK_RETRIED, or its element balance is off by
                       more than 1e-6 relative (pa_PSTALL may report a
                       stall as a success; thermochemistry/README.md), or
                       it gives the pair p_g <= 0 or a max log10 Omega
                       that is not finite; each such call counts as a
                       failure and is written to the log of the rank
  REJECTED             |log10 p_eq - log10 p_eq,table| > maxLog10Deviation

Warm starts: every element keeps its own GEMS3K state (GEMSB_WARM_CELL);
the state of a failed or rejected element is invalidated, so its next call
is cold.  GEMSB_ERR_FATAL re-creates the engine: destroyed, created again
from the system file, every condensed species suppressed again, and the
warm states of all elements of the rank invalidated (plan section 12,
item 7).  A step without an update keeps p_eq (and chi) of the last: the
lagged state.  It is written with the per-rank restart state
(uniform/phaseChange/phaseChangeEquilibrium, interfaceExchange.H: p_eq,
chi and the source of every element, and the counts of the update for its
step line, lastUpdate()), and a restart onto the same elements restores it
(restoreLastUpdate()), so up to its next update a restart repeats the
continuous run bit for bit.  The warm states are not written: that update
is cold and follows the continuous run to the tolerance of GEMS3K only
(plan section 8, risk 13).  A run without a lagged state (a fresh start; a
restart on another decomposition or after mode frozen or equilibrium
table) updates at its first step too and then keeps the schedule of the
time index (plan section 35, item C).

speciation lagged (a separate modelling option, off by default; plan E9).
GEMS3K also gives the speciation of the gas: chi_e = n_g/c_k, the fraction
of the pair's formula units present as the gas species g (PbI2(g); the
rest is PbI, I, Pb, I2), at the composition of the update, clipped to
[1e-12, 1].  The law counts only the molecules of g in the partial
pressure, q_e = A_e G'_e (p_eq - chi_e beta_c Y_c): g_e is multiplied by
chi_e, and with wallResistance halfCell the half cell enters with chi_e
beta_c (G'_e = G_e/(1 + G_e chi_e beta_c d_e/(rhoD)_e)).  An element
whose p_eq comes from the table has chi_e = 1.  local mode only.

allowNonPhysicalCarrier.  local evaluates the composition of the cell from
the carrier's p and T, which needs a carrier that is a gas at its p and T:
the start-up refuses (FatalIOError) a carrier whose density differs from
p W/(R T) by more than 1 % anywhere (paper mode, rho = 1), unless
allowNonPhysicalCarrier yes (then a warning).

Logs.  Every rank with an engine writes gemsEquilibrium.log into its own
log directory, logDirectory relative to its case directory (the case in a
serial run, processorN in parallel): the engine, every update with its
counts and its timing on the rank, every failure (element, T, p, x,
status, the message of GEMS3K, at most 20 per update), every re-creation;
GEMS3K's own ipmlog.txt goes there too.  In a parallel run a logDirectory
that does not lie below processorN (an absolute path, <case>/...,
$FOAM_CASE/..., or one that leaves processorN with ..) gets a subdirectory
processor<N> per rank: two ranks never share gemsEquilibrium.log or the
spdlog sink of ipmlog.txt (plan section 35, item D).  The decision is
collective, one rule for every rank (../processor0/gemsLog lies inside
processor0 for rank 0 only, and every rank takes its subdirectory), so
the start-up line names the directory of every rank ('logs in
<directory>/processor<N>'; plan section 35, item Y).  The main log has the
engine line at start-up, in local mode one GEMS line per pair in the step
report (deterministic: counts, iterations and deviations, no timing) and
the summary at the end of the run (with the mean time of a warm call, the
plan's cost criterion).  The counts are 64-bit integers, reduced in 64
bits: a long run exceeds a 32-bit label (9e8 calls of the 15264 wall
elements in 60 s at 1 ms, and about eight IPM iterations per call; plan
section 35, item A).

------------------------------------------------------------------------------*/

#ifndef gemsEquilibrium_H
#define gemsEquilibrium_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fileName.H"
#include "string.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "FixedList.H"
#include "vapourPressureTable.H"

#include <cstdint>
#include <fstream>

namespace LESTO {

class gemsEquilibrium {

public:

  enum class modeType { frozen, local };

  /* where the p_eq of an element of a pair comes from (gemsSource_<gas>) */
  enum sourceType {
    FROM_GEMS           = 0,
    BELOW_TEMPERATURE   = 1,
    BELOW_MOLE_FRACTION = 2,
    OUTSIDE_GRID        = 3,
    REJECTED            = 4,
    FAILED              = 5,
    nSourceTypes        = 6
  };

  /* the lower clip of the speciation factor chi */
  static constexpr double chiMin = 1e-12;

private:

  /* RAII owner of one engine of the bridge (gemsEquilibrium.C) */
  class engine;

  /*--------------------------------------------------------------------------
  One pair: the names of its gems block, and, once the engine exists, the
  indices of its species in the system and its bulk per formula unit.
  --------------------------------------------------------------------------*/
  struct pairSpec {
    Foam::word             gasName;      /* solver species, e.g. PbI2_g */
    Foam::scalar           molarMass;    /* [kg/mol], solver */
    Foam::word             gas;          /* GEMS3K species, e.g. PbI2(g) */
    Foam::wordList         condensates;  /* GEMS3K species */
    Foam::wordList         elementNames; /* formula: names and */
    Foam::scalarList       elementCoeffs;/*   stoichiometric coefficients */
    Foam::wordList         excessNames;  /* excess: names and */
    Foam::scalarList       excessValues; /*   relative excess */
    Foam::label            jGas;
    Foam::labelList        jCondensates;
    Foam::scalarList       bulk;         /* per system element: nu (1 + e) */
  };

  /* the counts of one update (or of the run): 64-bit integers, reduced in
     64 bits (endUpdate()); a 32-bit label overflows in a long run (the
     IPM iterations of 9e8 calls; header, Logs) */
  struct counts {
    std::int64_t calls;          /* calls of the engine */
    std::int64_t warm;           /* GEMSB_OK from a valid warm state */
    std::int64_t cold;           /* GEMSB_OK without a valid state */
    std::int64_t retried;        /* GEMSB_OK_RETRIED (warm failed, cold) */
    std::int64_t failed;         /* failed calls: a status other than
                                    OK and OK_RETRIED, an element balance
                                    off by more than 1e-6, or p_g <= 0 or
                                    a max log10 Omega that is not finite
                                    for a pair; a rejected result (not
                                    used either) is counted per pair, in
                                    REJECTED, not here */
    std::int64_t recreated;      /* engines re-created */
    std::int64_t iterations;     /* IPM iterations, all calls */
    std::int64_t updates;
    Foam::scalar warmSeconds;
    Foam::scalar coldSeconds;    /* cold and retried calls */
  };

  /* the sources of the elements of one pair in an update (or the run) */
  struct pairCounts {
    Foam::FixedList<std::int64_t, nSourceTypes> sources;
    Foam::scalar maxDeviation;   /* |dlog10 p_eq - table| of FROM_GEMS */
    Foam::scalar chiMin;
    Foam::scalar chiMax;
    std::int64_t chiClipped;
  };

  const Foam::fvMesh&    mesh_;
  const Foam::dictionary dict_;          /* GEMSCoeffs */

  /* GEMSCoeffs */
  Foam::fileName  systemFile_;           /* as written */
  Foam::fileName  system_;               /* expanded */
  modeType        mode_;
  bool            speciation_;
  Foam::word      carrier_;
  Foam::wordList  referenceNames_;
  Foam::scalarList referenceAmounts_;
  Foam::scalar    referencePressure_;
  Foam::scalar    minMoleFraction_;
  Foam::scalar    minTemperature_;
  Foam::scalar    maxLog10Deviation_;
  Foam::label     updateInterval_;
  bool            allowNonPhysicalCarrier_;
  Foam::fileName  logDirectory_;
  Foam::label     logLevel_;
  bool            writeField_;

  Foam::List<pairSpec> pairs_;

  /* the system, read from the engine at start() */
  Foam::label      nElements_;           /* of the system (ICs) */
  Foam::label      nSpecies_;            /* of the system (DCs) */
  Foam::label      nPhases_;
  Foam::label      stateSize_;
  Foam::wordList   elementNames_;
  Foam::wordList   speciesNames_;
  Foam::boolList   isGas_;
  Foam::scalarList stoichiometry_;       /* nSpecies x nElements */
  Foam::labelList  suppressed_;          /* the condensed species */
  Foam::scalar     Tmin_, Tmax_, Pmin_, Pmax_;
  Foam::label      jCarrier_;
  Foam::scalarList carrierFormula_;      /* per system element */

  /* the engine of this rank (local; frozen: during start() only) */
  Foam::autoPtr<engine> engine_;
  Foam::label           nEngines_;       /* created at start, all ranks */

  /* local: the warm state of every element of the rank */
  Foam::scalarList states_;
  Foam::List<int>  valid_;

  /* scratch of a call */
  Foam::scalarList b_, pGas_, log10a_, amounts_;

  /* frozen: the GEMS table of every pair, and its start-up report */
  Foam::PtrList<vapourPressureTable> frozen_;
  Foam::List<Foam::string>           frozenReport_;

  /* local: update control (the schedule follows the time index) */
  bool        updated_;                  /* this step */
  bool        lagged_;                   /* the elements hold the p_eq of
                                            an update (of this run, or
                                            restored) */

  /* local: counts of the last update (reduced) and of the run */
  counts                 update_;
  counts                 run_;
  Foam::List<pairCounts> pairUpdate_;
  Foam::List<pairCounts> pairRun_;
  Foam::label            failuresLogged_; /* this update, this rank */

  /* the log of this rank */
  Foam::fileName logFile_;
  std::ofstream  log_;

  static counts noCounts();
  static pairCounts noPairCounts();

  /*--------------------------------------------------------------------------
  The keys a dictionary may hold: stops (FatalIOError against d, naming
  what) at any other key, with the list of the valid ones (header,
  Dictionary; plan section 35, item H)
  --------------------------------------------------------------------------*/
  static void checkKeys (
    const Foam::dictionary& d,
    const Foam::wordList&   valid,
    const Foam::string&     what
  );

  /*--------------------------------------------------------------------------
  The log directory of this rank: logDirectory expanded, relative to the
  rank's case directory (Time::path(): the case, or processorN), cleaned;
  in a parallel run with a subdirectory processor<N> when it does not lie
  below processorN on some rank (header, Logs; collective).  shared tells
  whether that subdirectory was added (the same on every rank).
  --------------------------------------------------------------------------*/
  Foam::fileName logDirectory(bool& shared) const;

  /* the name and the coefficient of every entry of a sub-dictionary */
  static void readAmounts (
    const Foam::dictionary& d,
    const Foam::word&       key,
    Foam::wordList&         names,
    Foam::scalarList&       values,
    const bool              required
  );

  /* index of a system element or species; FatalIOError naming what if
     absent (dict is the dictionary of the name) */
  Foam::label elementIndex (
    const Foam::word&       name,
    const Foam::dictionary& dict,
    const Foam::string&     what
  ) const;
  Foam::label speciesIndex (
    const Foam::word&       name,
    const Foam::dictionary& dict,
    const Foam::string&     what
  ) const;

  /* creates the engine of this rank and suppresses every condensed
     species; false with the message of GEMS3K on failure */
  bool createEngine(Foam::string& message);

  /* the system description from the engine, the carrier and the pairs
     resolved and checked (FatalIOError otherwise) */
  void readSystem();

  /* frozen: the GEMS table of every pair from its table (master), and the
     broadcast of its values */
  void buildFrozenTables (
    const Foam::UList<const vapourPressureTable*>& tables
  );

  /* the element balance of the last call: true when every element with
     an amount above 1e-12 of the bulk is matched to 1e-6 relative */
  bool balanced() const;

  /* writes a line to the log of this rank (gemsEquilibrium.log) */
  void logLine(const std::string& line);

public:

  /* true when the executable was built with the GEMS3K bridge */
  static bool available();

  /* bridge directory used at build time; empty without the bridge */
  static Foam::string bridgeDirectory();

  /* one-line summary of the build state */
  static Foam::string buildState();

  /*--------------------------------------------------------------------------
  Without the bridge: stops the run with a FatalIOError reported against
  dict, the dictionary that selected equilibrium GEMS, with the rebuild
  hint (M9.1).  With the bridge: returns.
  --------------------------------------------------------------------------*/
  static void requireBridge(const Foam::dictionary& dict);

  /*--------------------------------------------------------------------------
  Reads GEMSCoeffs (dict); speciation: HKSCoeffs speciation lagged.  No
  engine yet (start()).  Refuses unknown keys and invalid values of the
  keys of its mode; names the keys of the other mode, which it does not
  read, in a warning (header, Dictionary).
  --------------------------------------------------------------------------*/
  gemsEquilibrium (
    const Foam::fvMesh&     mesh,
    const Foam::dictionary& dict,
    const bool              speciation
  );

  ~gemsEquilibrium();

  gemsEquilibrium(const gemsEquilibrium&) = delete;
  void operator=(const gemsEquilibrium&) = delete;

  /* adds a pair from its gems sub-dictionary, in the order of the pairs
     (unknown keys refused) */
  void addPair (
    const Foam::word&       gasName,
    const Foam::dictionary& gemsDict,
    const Foam::scalar      molarMass
  );

  /*--------------------------------------------------------------------------
  Creates the engine of every rank, reads and checks the system against the
  carrier and the pairs, and prints the engine lines.  frozen: builds the
  GEMS table of every pair from its vapourPressure table (tables, in the
  order of the pairs; the master evaluates, every rank receives it) and
  destroys the engines.  local: keeps the engine and sizes the warm states
  for nElements elements of this rank.  Collective.
  --------------------------------------------------------------------------*/
  void start (
    const Foam::label                              nElements,
    const Foam::UList<const vapourPressureTable*>& tables
  );

  modeType mode() const { return mode_; }
  bool local() const { return mode_ == modeType::local; }
  bool speciation() const { return speciation_; }
  bool allowNonPhysicalCarrier() const { return allowNonPhysicalCarrier_; }
  bool writeEquilibriumField() const { return writeField_; }
  const Foam::dictionary& dict() const { return dict_; }

  /*--------------------------------------------------------------------------
  frozen: p_eq [Pa] of pair k at T [K]: from its GEMS table inside the
  table's range and at or above minTemperature; otherwise tableValue, the
  pair's table at T.  source returns the sourceType.
  --------------------------------------------------------------------------*/
  Foam::scalar frozenPEq (
    const Foam::label  k,
    const Foam::scalar T,
    const Foam::scalar tableValue,
    Foam::label&       source
  ) const;

  /* local: whether the step of time index n updates by the schedule,
     (n - 1) mod updateInterval = 0, n the absolute time index (header,
     mode local) */
  bool updateDueAt(const Foam::label timeIndex) const;

  /* local: whether the step of time index n (Time::timeIndex() of the
     step) updates: by the schedule, or because the elements hold no
     lagged state yet; sets updated() */
  bool updateDue(const Foam::label timeIndex);

  /* local: whether this step updated */
  bool updated() const { return updated_; }

  /* local: whether the elements hold the p_eq of an update (lagged) */
  bool lagged() const { return lagged_; }

  Foam::label updateInterval() const { return updateInterval_; }

  /*--------------------------------------------------------------------------
  local: the counts of the last update per pair (its step line; the same on
  every rank, reduced), nPairState values per pair in the order of the
  pairs: the elements of each source, the largest deviation, the range of
  chi and the number clipped.  lastUpdate() writes them for the per-rank
  restart state (interfaceExchange: uniform/phaseChange/
  phaseChangeEquilibrium, with p_eq, chi and the source of every element);
  restoreLastUpdate() takes them back on a restart that restores the
  elements' lagged state, so that its steps without an update print the
  lines of the continuous run.  Values are exact: counts as doubles (below
  2^53).
  --------------------------------------------------------------------------*/
  static constexpr Foam::label nPairState = nSourceTypes + 4;
  Foam::scalarList lastUpdate() const;
  void restoreLastUpdate(const Foam::UList<Foam::scalar>& values);

  /* local: an update is beginUpdate(), evaluate() of every element of the
     rank, endUpdate() (collective: the counts are reduced) */
  void beginUpdate();

  /*--------------------------------------------------------------------------
  local: evaluates element e (index among the elements of the rank) at its
  temperature T [K], the pressure p [Pa] of its cell and the composition of
  its cell, the carrier concentration cCarrier [mol/m3] and per pair k the
  concentration c[k] = rho Y_k/M_k [mol/m3] of its formula units; active[k]
  tells whether the element belongs to pair k.  For every active pair it
  returns pEq[k] [Pa], chi[k] and source[k]: the GEMS3K values, or
  pEqTable[k] (the pair's table at T) and chi 1 where a guard applies.
  The other entries are left unchanged.
  --------------------------------------------------------------------------*/
  void evaluate (
    const Foam::label          e,
    const Foam::scalar         T,
    const Foam::scalar         p,
    const Foam::scalar         cCarrier,
    const Foam::UList<Foam::scalar>& c,
    const Foam::UList<bool>&   active,
    const Foam::UList<Foam::scalar>& pEqTable,
    Foam::UList<Foam::scalar>& pEq,
    Foam::UList<Foam::scalar>& chi,
    Foam::UList<Foam::label>&  source
  );

  void endUpdate();

  /* local: the GEMS part of the step report of pair k (deterministic) */
  Foam::string stepReport(const Foam::label k) const;

  /* end-of-run summary (local: the counts and the cost of the calls) */
  void summary() const;
};

}

#endif
