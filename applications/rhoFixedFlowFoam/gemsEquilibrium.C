/*------------------------------------------------------------------------------

gemsEquilibrium.C

PURPOSE

Implementation of the optional GEMS3K backend of the HKS law (see
gemsEquilibrium.H for the model, the dictionary and the guards).  Every
call of the bridge is inside #ifdef LESTO_HAVE_GEMS, which the generated
gemsConfig.H defines only for a build with LESTO_GEMSBRIDGE set.  Without
it the dictionary is still read and checked here, so that the code of both
builds is compiled, but the engine, start() and evaluate() are stubs that
are never reached: requireBridge() stops a run that selects equilibrium
GEMS before the class is constructed.

gemsConfig.H must be included here, before any test of LESTO_HAVE_GEMS: the
include is what makes wmake recompile this file when the bridge is toggled.

The bridge header is included with angle brackets on purpose.  wmkdepend
follows only quoted includes and ignores #ifdef, so a quoted include would
make it report a missing gemsbridge.h in every build without the bridge.
wmkdepend therefore does not record gemsbridge.h itself.  Instead, a bridge
build's gemsConfig.H carries the checksum of gemsbridge.h, so an edit of the
bridge header rewrites gemsConfig.H and recompiles this file.

Parallel hygiene.  The engines never communicate: there is no MPI inside
the bridge, and one engine serves one rank.  Everything collective is
called on every rank: the outcome of the creation in start() is reduced,
so that a failure on one rank stops all of them; the frozen tables are
evaluated on the master and broadcast; endUpdate() reduces the counts.
evaluate() is local to its rank, including the re-creation of the engine
after GEMSB_ERR_FATAL.

------------------------------------------------------------------------------*/

#include "gemsConfig.H"
#include "gemsEquilibrium.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "StringStream.H"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>

#ifdef LESTO_HAVE_GEMS
#include <gemsbridge.h>
#endif

using namespace Foam;


/*------------------------------------------------------------------------------
The engine: RAII owner of one gemsb_engine of the bridge.  It is created
from the system file and destroyed with its owner; it cannot be copied.
Without the bridge it is an empty class that is never constructed.
------------------------------------------------------------------------------*/

class LESTO::gemsEquilibrium::engine {

#ifdef LESTO_HAVE_GEMS
  gemsb_engine* handle_;
  std::string   message_;

public:

  explicit engine(const fileName& system)
  :
    handle_(nullptr),
    message_()
  {
    char message[1024] = "";
    handle_ = gemsb_create_from_lst(system.c_str(), message, sizeof(message));
    if (!handle_) {
      message_ = message[0] ? message : "no message";
    }
  }

  ~engine() {
    if (handle_) {
      gemsb_destroy(handle_);
    }
  }

  engine(const engine&) = delete;
  void operator=(const engine&) = delete;

  bool valid() const { return handle_ != nullptr; }
  const std::string& message() const { return message_; }
  gemsb_engine* handle() const { return handle_; }
#endif
};


/*------------------------------------------------------------------------------
Local helpers
------------------------------------------------------------------------------*/

namespace {

#ifdef LESTO_HAVE_GEMS
  /* the name of a return code of gemsb_equilibrate() */
  const char* statusName(const int rc) {
    switch (rc) {
      case GEMSB_OK: return "OK";
      case GEMSB_OK_RETRIED: return "OK_RETRIED";
      case GEMSB_BAD_QUALITY: return "BAD_QUALITY";
      case GEMSB_ERR_NOCONV: return "ERR_NOCONV";
      case GEMSB_ERR_FATAL: return "ERR_FATAL";
      case GEMSB_ERR_INPUT: return "ERR_INPUT";
      default: return "unknown";
    }
  }
#endif

  /* text of a value at the default precision of Info */
  template<class Type>
  std::string text(const Type& value) {
    OStringStream os;
    os << value;
    return os.str();
  }

  /* "name1 value1 name2 value2 ..." */
  std::string amounts(const wordList& names, const scalarList& values) {
    std::string s;
    forAll(names, i) {
      s += (i ? " " : "") + names[i] + " " + text(values[i]);
    }
    return s;
  }

  /*--------------------------------------------------------------------------
  A 64-bit count as text.  Not through Foam::Ostream: OpenFOAM v2412 writes
  an int64_t as label(val) (int64IO.C), which an Int32 build truncates to
  32 bits (v2606 writes it in full); std::to_string is exact in both.
  --------------------------------------------------------------------------*/
  std::string num(const std::int64_t n) {
    return std::to_string(n);
  }

  /* "(n1 n2 ...)" of 64-bit counts */
  template<unsigned N>
  std::string nums(const FixedList<std::int64_t, N>& n) {
    std::string s("(");
    for (unsigned i = 0; i < N; ++i) {
      s += (i ? " " : "") + num(n[i]);
    }
    return s + ")";
  }

  /* "name1 name2 ..." */
  std::string joined(const wordList& names) {
    std::string s;
    forAll(names, i) {
      s += (i ? " " : "") + names[i];
    }
    return s;
  }

#ifdef LESTO_HAVE_GEMS
  /*--------------------------------------------------------------------------
  The largest log10 Omega of the condensates j of a pair (log10a: the dual
  activities of the call); false, with lgOmega the first value that is not
  finite, when one of them is not.  A value is compared only once it is
  known to be finite: Foam::max is (a > b) ? a : b, which passes over a NaN
  (max(-VGREAT, NaN) is NaN, then max(NaN, x) is x) and, as an ordered
  comparison, may raise FE_INVALID, which OpenFOAM traps (FOAM_SIGFPE): the
  build before plan section 35 stopped with a floating-point exception at a
  NaN log10 Omega of a converged call.
  --------------------------------------------------------------------------*/
  bool maxLog10Omega (
    const labelList&  condensates,
    const scalarList& log10a,
    scalar&           lgOmega
  ) {
    lgOmega = -VGREAT;
    for (const label j : condensates) {
      if (!std::isfinite(log10a[j])) {
        lgOmega = log10a[j];
        return false;
      }
      lgOmega = max(lgOmega, log10a[j]);
    }
    return true;
  }

  /* p_g [Pa] of a pair is usable: finite (checked first, see above) and
     positive */
  bool positive(const scalar pg) {
    return std::isfinite(pg) && pg > 0;
  }
#endif
}


/*------------------------------------------------------------------------------
Build state
------------------------------------------------------------------------------*/

bool LESTO::gemsEquilibrium::available() {
#ifdef LESTO_HAVE_GEMS
  return true;
#else
  return false;
#endif
}


Foam::string LESTO::gemsEquilibrium::bridgeDirectory() {
  return Foam::string(LESTO_GEMSBRIDGE_DIR);
}


Foam::string LESTO::gemsEquilibrium::buildState() {
  if (available()) {
    return "GEMS3K bridge: compiled in from " + bridgeDirectory();
  }
  return "GEMS3K bridge: not compiled in (set LESTO_GEMSBRIDGE and rerun "
         "Allwmake to enable it)";
}


void LESTO::gemsEquilibrium::requireBridge(const dictionary& dict) {

  if (available()) {
    return;
  }

  FatalIOErrorInFunction(dict)
    << "equilibrium GEMS needs rhoFixedFlowFoam built with the GEMS3K "
    << "bridge, but this executable was built without it." << nl
    << "Rebuild with" << nl
    << "    export LESTO_GEMSBRIDGE=<directory with gemsbridge.h and "
    << "libgemsbridge.so>" << nl
    << "    ./Allwmake" << nl
    << "or select equilibrium table." << exit(FatalIOError);
}


/*------------------------------------------------------------------------------
Counts
------------------------------------------------------------------------------*/

LESTO::gemsEquilibrium::counts LESTO::gemsEquilibrium::noCounts() {
  counts c;
  c.calls = 0;
  c.warm = 0;
  c.cold = 0;
  c.retried = 0;
  c.failed = 0;
  c.recreated = 0;
  c.iterations = 0;
  c.warmSeconds = 0;
  c.coldSeconds = 0;
  c.updates = 0;
  return c;
}


LESTO::gemsEquilibrium::pairCounts LESTO::gemsEquilibrium::noPairCounts() {
  pairCounts c;
  c.sources = std::int64_t(0);
  c.maxDeviation = 0;
  c.chiMin = GREAT;
  c.chiMax = -GREAT;
  c.chiClipped = 0;
  return c;
}


/*------------------------------------------------------------------------------
Construction: GEMSCoeffs (the engine follows in start())
------------------------------------------------------------------------------*/

void LESTO::gemsEquilibrium::readAmounts (
  const dictionary& d,
  const word&       key,
  wordList&         names,
  scalarList&       values,
  const bool        required
) {
  names.clear();
  values.clear();

  if (!d.found(key, keyType::LITERAL)) {
    if (required) {
      FatalIOErrorInFunction(d) << "Missing " << key << " { <element> "
        << "<amount>; ... } (elements of the GEMS3K system)."
        << exit(FatalIOError);
    }
    return;
  }

  const dictionary& sub = d.subDict(key);
  for (const entry& e : sub) {
    const scalar value = sub.get<scalar>(e.keyword());
    if (!std::isfinite(value)) {
      FatalIOErrorInFunction(sub) << key << ": " << e.keyword() << " "
        << value << " is not finite." << exit(FatalIOError);
    }
    names.append(e.keyword());
    values.append(value);
  }

  if (required && names.empty()) {
    FatalIOErrorInFunction(sub) << key << " lists no element."
      << exit(FatalIOError);
  }
}


void LESTO::gemsEquilibrium::checkKeys (
  const dictionary& d,
  const wordList&   valid,
  const string&     what
) {
  wordList unknown;
  for (const entry& e : d) {
    if (!valid.found(e.keyword())) {
      unknown.append(e.keyword());
    }
  }
  if (unknown.size()) {
    FatalIOErrorInFunction(d) << what.c_str() << ": unknown "
      << (unknown.size() == 1 ? "key " : "keys ") << joined(unknown).c_str()
      << ".  The valid keys are " << joined(valid).c_str() << "."
      << exit(FatalIOError);
  }
}


/*------------------------------------------------------------------------------
The keys of GEMSCoeffs.  Only the keys of the mode are read: a key of the
other mode is named in a warning and never validated, so it cannot stop
the run (plan section 35, item G).  Every key must be one of these (item
H): a misspelled key would otherwise leave its default in place silently.
------------------------------------------------------------------------------*/

LESTO::gemsEquilibrium::gemsEquilibrium (
  const fvMesh&     mesh,
  const dictionary& dict,
  const bool        speciation
)
:
  mesh_(mesh),
  dict_(dict),
  systemFile_(),
  system_(),
  mode_(modeType::local),
  speciation_(speciation),
  carrier_(),
  referenceNames_(),
  referenceAmounts_(),
  referencePressure_(101325),
  minMoleFraction_(1e-6),
  maxMoleFraction_(3.2e-4),
  minTemperature_(500),
  maxLog10Deviation_(0.01),
  updateInterval_(1),
  allowNonPhysicalCarrier_(false),
  logDirectory_("gemsLog"),
  logLevel_(4),
  writeField_(false),
  pairs_(),
  nElements_(0),
  nSpecies_(0),
  nPhases_(0),
  stateSize_(0),
  elementNames_(),
  speciesNames_(),
  isGas_(),
  stoichiometry_(),
  suppressed_(),
  Tmin_(0),
  Tmax_(0),
  Pmin_(0),
  Pmax_(0),
  jCarrier_(-1),
  carrierFormula_(),
  engine_(),
  nEngines_(0),
  states_(),
  valid_(),
  b_(),
  pGas_(),
  log10a_(),
  amounts_(),
  frozen_(),
  frozenReport_(),
  updated_(false),
  lagged_(false),
  update_(noCounts()),
  run_(noCounts()),
  pairUpdate_(),
  pairRun_(),
  failuresLogged_(0),
  logFile_(),
  log_()
{
  checkKeys (
    dict,
    wordList ({
      "system", "mode", "carrier", "referenceComposition",
      "referencePressure", "minMoleFraction", "minTemperature",
      "maxLog10Deviation", "updateInterval", "allowNonPhysicalCarrier",
      "logDirectory", "logLevel", "writeEquilibriumField", "maxMoleFraction"
    }),
    "GEMSCoeffs"
  );

  const word mode(dict.get<word>("mode"));
  if (mode == "frozen") {
    mode_ = modeType::frozen;
  } else if (mode != "local") {
    FatalIOErrorInFunction(dict) << "GEMSCoeffs: unknown mode " << mode
      << ".  Choose frozen (start-up: p_eq(T) on the nodes of the table) "
      << "or local (per step and element, warm started)."
      << exit(FatalIOError);
  }

  /* the keys of both modes */
  systemFile_ = dict.get<fileName>("system");
  system_ = systemFile_;
  system_.expand();
  carrier_ = dict.get<word>("carrier");
  minTemperature_ = dict.getOrDefault<scalar>("minTemperature", 500);
  maxLog10Deviation_ = dict.getOrDefault<scalar>("maxLog10Deviation", 0.01);
  logDirectory_ = dict.getOrDefault<fileName>("logDirectory", "gemsLog");
  logLevel_ = dict.getOrDefault<label>("logLevel", 4);
  writeField_ = dict.getOrDefault<bool>("writeEquilibriumField", false);

  if (!(minTemperature_ >= 0) || !std::isfinite(minTemperature_)
    || !(maxLog10Deviation_ > 0) || !std::isfinite(maxLog10Deviation_)
    || logLevel_ < 0 || logLevel_ > 6) {
    FatalIOErrorInFunction(dict) << "GEMSCoeffs: require minTemperature "
      << ">= 0, maxLog10Deviation > 0 (both finite) and logLevel 0 to 6, "
      << "not " << minTemperature_ << ", " << maxLog10Deviation_ << ", "
      << logLevel_ << "." << exit(FatalIOError);
  }

  maxMoleFraction_=dict.getOrDefault<scalar>("maxMoleFraction",3.2e-4);
  if(!std::isfinite(maxMoleFraction_) || maxMoleFraction_<=0 || maxMoleFraction_>1)
    FatalIOErrorInFunction(dict)<<"Require finite maxMoleFraction in (0,1]"<<exit(FatalIOError);

  /* the keys of the mode */
  if (mode_ == modeType::frozen) {
    referencePressure_ =
      dict.getOrDefault<scalar>("referencePressure", 101325);
    if (!(referencePressure_ > 0) || !std::isfinite(referencePressure_)) {
      FatalIOErrorInFunction(dict) << "GEMSCoeffs (mode frozen): require "
        << "referencePressure > 0 (finite), not " << referencePressure_
        << "." << exit(FatalIOError);
    }
    readAmounts(dict, "referenceComposition", referenceNames_,
                referenceAmounts_, true);
    scalar sum = 0;
    for (const scalar amount : referenceAmounts_) {
      if (amount < 0) {
        FatalIOErrorInFunction(dict) << "GEMSCoeffs: referenceComposition "
          << "holds a negative amount." << exit(FatalIOError);
      }
      sum += amount;
    }
    if (!(sum > 0)) {
      FatalIOErrorInFunction(dict) << "GEMSCoeffs: referenceComposition "
        << "holds no positive amount." << exit(FatalIOError);
    }
  } else {
    minMoleFraction_ = dict.getOrDefault<scalar>("minMoleFraction", 1e-6);
    updateInterval_ = dict.getOrDefault<label>("updateInterval", 1);
    allowNonPhysicalCarrier_ =
      dict.getOrDefault<bool>("allowNonPhysicalCarrier", false);
    if (!(minMoleFraction_ >= 0) || !std::isfinite(minMoleFraction_)
      || updateInterval_ < 1) {
      FatalIOErrorInFunction(dict) << "GEMSCoeffs (mode local): require "
        << "minMoleFraction >= 0 (finite) and updateInterval >= 1, not "
        << minMoleFraction_ << " and " << updateInterval_ << "."
        << exit(FatalIOError);
    }
  }

  /*--------------------------------------------------------------------------
  speciation lagged takes the speciation of the cell's gas from the local
  update; frozen has no cell composition.
  --------------------------------------------------------------------------*/
  if (speciation_ && mode_ == modeType::frozen) {
    FatalIOErrorInFunction(dict) << "speciation lagged needs GEMSCoeffs "
      << "mode local: it takes the gas speciation of every element from "
      << "the GEMS3K state of its last update." << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  Keys of the other mode are not read; they may stay in the file (to switch
  the mode back), but they are named.
  --------------------------------------------------------------------------*/
  wordList ignored;
  if (mode_ == modeType::local) {
    for (const word key : {"referenceComposition", "referencePressure"}) {
      if (dict.found(key, keyType::LITERAL)) {
        ignored.append(key);
      }
    }
  } else {
    for (const word key : {"minMoleFraction", "updateInterval",
                           "allowNonPhysicalCarrier"}) {
      if (dict.found(key, keyType::LITERAL)) {
        ignored.append(key);
      }
    }
  }
  if (ignored.size()) {
    WarningInFunction << "GEMSCoeffs: mode " << mode << " does not read "
      << ignored << " (keys of the other mode)." << nl << endl;
  }
}


LESTO::gemsEquilibrium::~gemsEquilibrium() {}


/*------------------------------------------------------------------------------
A pair: its gems block (resolved against the system in start())
------------------------------------------------------------------------------*/

void LESTO::gemsEquilibrium::addPair (
  const word&       gasName,
  const dictionary& gemsDict,
  const scalar      molarMass
) {
  checkKeys (
    gemsDict,
    wordList({"gas", "condensates", "elements", "excess"}),
    "Pair " + gasName + ", gems"
  );

  pairSpec p;
  p.gasName = gasName;
  p.molarMass = molarMass;
  p.gas = gemsDict.get<word>("gas");
  p.condensates = gemsDict.get<wordList>("condensates");
  readAmounts(gemsDict, "elements", p.elementNames, p.elementCoeffs, true);
  readAmounts(gemsDict, "excess", p.excessNames, p.excessValues, false);
  p.jGas = -1;

  if (p.condensates.empty()) {
    FatalIOErrorInFunction(gemsDict) << "Pair " << gasName << ": gems "
      << "condensates lists no species." << exit(FatalIOError);
  }
  for (const scalar nu : p.elementCoeffs) {
    if (!(nu > 0)) {
      FatalIOErrorInFunction(gemsDict) << "Pair " << gasName << ": gems "
        << "elements must be positive coefficients of the formula."
        << exit(FatalIOError);
    }
  }
  forAll(p.excessNames, i) {
    if (!(p.excessValues[i] >= 0 && p.excessValues[i] < 1)) {
      FatalIOErrorInFunction(gemsDict) << "Pair " << gasName << ": gems "
        << "excess " << p.excessNames[i] << " " << p.excessValues[i]
        << ": a relative excess must lie in [0, 1)." << exit(FatalIOError);
    }
    if (!p.elementNames.found(p.excessNames[i])) {
      FatalIOErrorInFunction(gemsDict) << "Pair " << gasName << ": gems "
        << "excess names " << p.excessNames[i] << ", which is not an "
        << "element of its formula " << p.elementNames << "."
        << exit(FatalIOError);
    }
  }

  pairs_.append(p);
}


/*------------------------------------------------------------------------------
Indices of the system
------------------------------------------------------------------------------*/

Foam::label LESTO::gemsEquilibrium::elementIndex (
  const word&       name,
  const dictionary& dict,
  const string&     what
) const {
  const label i = elementNames_.find(name);
  if (i < 0) {
    FatalIOErrorInFunction(dict) << what.c_str() << ": " << name
      << " is not an element of the GEMS3K system " << systemFile_
      << " (elements " << elementNames_ << ")." << exit(FatalIOError);
  }
  return i;
}


Foam::label LESTO::gemsEquilibrium::speciesIndex (
  const word&       name,
  const dictionary& dict,
  const string&     what
) const {
  const label j = speciesNames_.find(name);
  if (j < 0) {
    FatalIOErrorInFunction(dict) << what.c_str() << ": " << name
      << " is not a species of the GEMS3K system " << systemFile_
      << " (species " << speciesNames_ << ")." << exit(FatalIOError);
  }
  return j;
}


/*------------------------------------------------------------------------------
The engine of this rank, with every condensed species suppressed (bounds 0
and 0): the call then returns the homogeneous gas equilibrium.  Used at
start-up and to re-create the engine after GEMSB_ERR_FATAL (the bounds are
re-applied, plan section 12, item 7).
------------------------------------------------------------------------------*/

bool LESTO::gemsEquilibrium::createEngine(string& message) {
#ifdef LESTO_HAVE_GEMS
  engine_.reset(nullptr);
  engine_.reset(new engine(system_));
  if (!engine_->valid()) {
    message = engine_->message();
    engine_.reset(nullptr);
    return false;
  }
  gemsb_engine* h = engine_->handle();
  if(gemsb_element_index(h,"Bi")>=0)gemsb_set_warm_iteration_limit(h,200);
  suppressed_.clear();
  for (int j = 0; j < gemsb_num_species(h); ++j) {
    if (!gemsb_species_is_gas(h, j)) {
      gemsb_set_species_bounds(h, j, 0.0, 0.0);
      suppressed_.append(j);
    }
  }
  return true;
#else
  message = "this executable was built without the GEMS3K bridge";
  return false;
#endif
}


/*------------------------------------------------------------------------------
The system description, and the carrier, the reference composition and the
pairs resolved against it.  Every rank reads its own engine; the checks are
the same on every rank, so a FatalIOError stops all of them.
------------------------------------------------------------------------------*/

void LESTO::gemsEquilibrium::readSystem() {
#ifdef LESTO_HAVE_GEMS
  gemsb_engine* h = engine_->handle();

  nElements_ = gemsb_num_elements(h);
  nSpecies_ = gemsb_num_species(h);
  nPhases_ = gemsb_num_phases(h);
  stateSize_ = gemsb_state_size(h);

  elementNames_.resize(nElements_);
  for (label i = 0; i < nElements_; ++i) {
    elementNames_[i] = word(gemsb_element_name(h, int(i)));
  }
  dualPairs_=elementNames_.found("Bi");
  speciesNames_.resize(nSpecies_);
  isGas_.resize(nSpecies_);
  stoichiometry_.resize(nSpecies_*nElements_);
  for (label j = 0; j < nSpecies_; ++j) {
    speciesNames_[j] = word(gemsb_species_name(h, int(j)));
    isGas_[j] = gemsb_species_is_gas(h, int(j)) != 0;
    for (label i = 0; i < nElements_; ++i) {
      stoichiometry_[j*nElements_ + i] = gemsb_stoich(h, int(j), int(i));
    }
  }
  gemsb_TP_range(h, &Tmin_, &Tmax_, &Pmin_, &Pmax_);

  b_.resize(nElements_, 0);
  pGas_.resize(nSpecies_, 0);
  log10a_.resize(nSpecies_, 0);
  amounts_.resize(nSpecies_, 0);
#endif

  /*--------------------------------------------------------------------------
  The carrier: a gas species; its formula is its row of the stoichiometry.
  --------------------------------------------------------------------------*/
  jCarrier_ = speciesIndex(carrier_, dict_, "GEMSCoeffs carrier");
  if (!isGas_[jCarrier_]) {
    FatalIOErrorInFunction(dict_) << "GEMSCoeffs carrier " << carrier_
      << " is not a gas species of the system " << systemFile_ << "."
      << exit(FatalIOError);
  }
  carrierFormula_.resize(nElements_);
  for (label i = 0; i < nElements_; ++i) {
    carrierFormula_[i] = stoichiometry_[jCarrier_*nElements_ + i];
  }

  forAll(referenceNames_, n) {
    elementIndex(referenceNames_[n], dict_, "GEMSCoeffs referenceComposition");
  }

  /*--------------------------------------------------------------------------
  The pairs: the gas species (a gas), its condensates (condensed), one
  formula for all of them (the mass balance of a pair assumes it), equal to
  the elements of the gems block, and the molar mass of the gas species
  equal to the solver's molarMass (to 1e-6: the atomic weights of the
  system and of speciesTransportProperties).
  --------------------------------------------------------------------------*/
  forAll(pairs_, k) {
    pairSpec& p = pairs_[k];
    const string where("Pair " + p.gasName + ", gems");

    p.jGas = speciesIndex(p.gas, dict_, where + " gas");
    if (!isGas_[p.jGas]) {
      FatalIOErrorInFunction(dict_) << where.c_str() << " gas " << p.gas
        << " is not a gas species of the system " << systemFile_ << "."
        << exit(FatalIOError);
    }

    scalarList formula(nElements_, 0.0);
    forAll(p.elementNames, n) {
      formula[elementIndex(p.elementNames[n], dict_, where + " elements")] =
        p.elementCoeffs[n];
    }

    auto sameFormula = [&](const label j) {
      for (label i = 0; i < nElements_; ++i) {
        if (mag(stoichiometry_[j*nElements_ + i] - formula[i])
          > 1e-12*max(mag(formula[i]), scalar(1))) {
          return false;
        }
      }
      return true;
    };

    auto row = [&](const label j) {
      std::string s;
      for (label i = 0; i < nElements_; ++i) {
        if (stoichiometry_[j*nElements_ + i] != 0) {
          s += (s.empty() ? "" : " ") + elementNames_[i] + " "
             + text(stoichiometry_[j*nElements_ + i]);
        }
      }
      return s;
    };

    if (!sameFormula(p.jGas)) {
      FatalIOErrorInFunction(dict_) << where.c_str() << " elements { "
        << amounts(p.elementNames, p.elementCoeffs).c_str() << " } differ "
        << "from the formula of " << p.gas << " in the system, "
        << row(p.jGas).c_str() << "." << exit(FatalIOError);
    }

    p.jCondensates.clear();
    for (const word& name : p.condensates) {
      const label j = speciesIndex(name, dict_, where + " condensates");
      if (isGas_[j]) {
        FatalIOErrorInFunction(dict_) << where.c_str() << " condensates: "
          << name << " is a gas species." << exit(FatalIOError);
      }
      if (!sameFormula(j)) {
        FatalIOErrorInFunction(dict_) << where.c_str() << " condensates: "
          << name << " (" << row(j).c_str() << ") has another formula than "
          << p.gas << " (" << row(p.jGas).c_str() << "): the mass balance "
          << "of a pair assumes one formula." << exit(FatalIOError);
      }
      p.jCondensates.append(j);
    }

#ifdef LESTO_HAVE_GEMS
    const scalar M = gemsb_species_molar_mass(engine_->handle(), int(p.jGas));
    if (!(mag(M/p.molarMass - 1) <= 1e-6)) {
      FatalIOErrorInFunction(dict_) << where.c_str() << " gas " << p.gas
        << ": the molar mass of the system, " << M << " kg/mol, differs "
        << "from molarMass " << p.molarMass << " kg/mol of "
        << p.gasName << " (speciesTransportProperties) by more than 1e-6."
        << exit(FatalIOError);
    }
#endif

    p.bulk = formula;
    forAll(p.excessNames, n) {
      const label i = elementIndex(p.excessNames[n], dict_, where + " excess");
      p.bulk[i] = formula[i]*(1.0 + p.excessValues[n]);
    }
  }
}


/*------------------------------------------------------------------------------
The log directory of this rank (header, Logs; plan section 35, item D).
Time::path() is the rank's case directory: the case in a serial run,
<case>/processorN in a parallel one.  A logDirectory outside it (absolute,
<case>/... and $FOAM_CASE/... after the expansion, or ../ out of
processorN) would be shared by every rank.  The decision is collective:
when the directory lies outside its processorN on any rank, every rank
takes its subdirectory processor<N>, so that one rule holds for all ranks
and the start-up line can name it (with ../processor0/gemsLog rank 0 alone
would stay inside; plan section 35, item Y).
------------------------------------------------------------------------------*/

Foam::fileName LESTO::gemsEquilibrium::logDirectory(bool& shared) const {

  fileName own(mesh_.time().path());
  own.clean();

  fileName dir(logDirectory_);
  dir.expand();
  if (!dir.isAbsolute()) {
    dir = own/dir;
  }
  dir.clean();

  const bool outside = dir != own && !dir.starts_with(own + '/');
  shared = Pstream::parRun() && returnReduceOr(outside);
  if (shared) {
    dir = dir/("processor" + Foam::name(Pstream::myProcNo()));
  }
  return dir;
}


/*------------------------------------------------------------------------------
Start: the engines, the checks, the frozen tables or the warm states
------------------------------------------------------------------------------*/

void LESTO::gemsEquilibrium::start (
  const label                              nElements,
  const UList<const vapourPressureTable*>& tables
) {
#ifdef LESTO_HAVE_GEMS

  /*--------------------------------------------------------------------------
  The log directory of this rank (logDirectory()), for GEMS3K's ipmlog.txt
  and for gemsEquilibrium.log.  The directory and the level of GEMS3K's
  loggers are global settings of the process: they must be set before the
  first engine is created (GEMS3K opens ipmlog.txt then).
  --------------------------------------------------------------------------*/
  bool shared = false;
  const fileName logDir(logDirectory(shared));
  mkDir(logDir);
  gemsb_set_log_directory(logDir.c_str());
  gemsb_set_log_level(int(logLevel_));

  logFile_ = logDir/"gemsEquilibrium.log";
  log_.open(logFile_.c_str(), std::ios::out | std::ios::app);

  string message;
  const bool created = createEngine(message);
  nEngines_ = returnReduce(label(created ? 1 : 0), sumOp<label>());

  if (nEngines_ < Pstream::nProcs()) {
    if (!created) {
      Pout<< "GEMS3K engine not created from " << system_ << ": "
        << message.c_str() << nl;
      logLine("engine not created from " + system_ + ": " + message);
    }
    FatalIOErrorInFunction(dict_) << "The GEMS3K engine could not be "
      << "created from " << systemFile_ << " (" << system_ << ") on "
      << Pstream::nProcs() - nEngines_ << " of " << Pstream::nProcs()
      << " ranks" << (created ? "" : (": " + message).c_str()) << "."
      << exit(FatalIOError);
  }

  readSystem();

  /* the line of the system, in the main log and in the log of every rank */
  OStringStream os;
  os << "system " << systemFile_.c_str() << " (" << nElements_
    << " elements " << joined(elementNames_).c_str() << ", " << nSpecies_
    << " species, " << nPhases_
    << " phases; T " << Tmin_ << "-" << Tmax_ << " K, P " << Pmin_ << "-"
    << Pmax_ << " Pa; " << stateSize_ << " values of warm state per "
    << "element); carrier " << carrier_ << "; suppressed condensed species";
  for (const label j : suppressed_) {
    os << ' ' << speciesNames_[j];
  }
  const string systemLine(os.str());

  logLine
  (
    "rank " + text(Pstream::myProcNo()) + " of " + text(Pstream::nProcs())
  + ", start time " + mesh_.time().timeName() + ", mode "
  + (local() ? "local" : "frozen") + ": " + systemLine
  );

  Info<< "GEMS3K backend (equilibrium GEMS, mode "
    << (local() ? "local" : "frozen") << "): " << systemLine.c_str() << nl
    << "  guards: minTemperature " << minTemperature_ << " K";
  if (local()) {
    Info<< ", minMoleFraction " << minMoleFraction_;
    if(dualPairs_)Info<<", maxMoleFraction "<<maxMoleFraction_<<", dual p_eq, warm iteration cap 200";
  }
  Info<< ", maxLog10Deviation " << maxLog10Deviation_;
  if (local()) {
    Info<< "; updateInterval " << updateInterval_
      << "; allowNonPhysicalCarrier "
      << (allowNonPhysicalCarrier_ ? "yes" : "no") << "; speciation "
      << (speciation_ ? "lagged" : "none");
  } else {
    Info<< "; reference composition "
      << amounts(referenceNames_, referenceAmounts_).c_str() << " mol at "
      << referencePressure_ << " Pa";
  }
  /*--------------------------------------------------------------------------
  Where the logs are, as this rank resolved them, with processor<N> for its
  number in a parallel run: a directory with a subdirectory per rank
  (shared, the same decision on every rank) as the full path ending in
  processor<N>, checked first since rank 0 may lie inside its own
  processor0 then; otherwise below the rank's directory as before
  ('processor<N>/gemsLog').
  --------------------------------------------------------------------------*/
  std::string where(logDir);
  {
    fileName own(mesh_.time().path());
    own.clean();
    const std::string rank("processor" + Foam::name(Pstream::myProcNo()));
    if (shared) {
      where = where.substr(0, where.size() - rank.size()) + "processor<N>";
    } else if (where.compare(0, own.size() + 1, own + '/') == 0) {
      where = (Pstream::parRun() ? "processor<N>/" : "")
            + where.substr(own.size() + 1);
    }
  }
  Info<< "; " << nEngines_ << (nEngines_ == 1 ? " engine" : " engines")
    << " (one per rank), logs in " << where.c_str() << nl;

  forAll(pairs_, k) {
    const pairSpec& p = pairs_[k];
    Info<< "Pair " << p.gasName << ", GEMS: gas " << p.gas
      << ", condensates " << joined(p.condensates).c_str()
      << " (p_eq over the most stable),"
      << " formula " << amounts(p.elementNames, p.elementCoeffs).c_str();
    if (p.excessNames.size()) {
      Info<< ", excess in the bulk "
        << amounts(p.excessNames, p.excessValues).c_str();
    }
    Info<< nl;
  }

  pairUpdate_.resize(pairs_.size(), noPairCounts());
  pairRun_.resize(pairs_.size(), noPairCounts());

  if (mode_ == modeType::frozen) {
    if (referencePressure_ < Pmin_ || referencePressure_ > Pmax_) {
      FatalIOErrorInFunction(dict_) << "GEMSCoeffs referencePressure "
        << referencePressure_ << " Pa lies outside the pressure grid ["
        << Pmin_ << ", " << Pmax_ << "] Pa of the system."
        << exit(FatalIOError);
    }
    buildFrozenTables(tables);
    engine_.reset(nullptr);
    logLine("frozen: GEMS tables built, engine destroyed");
    return;
  }

  /* local: the warm states of the elements of this rank, all cold; no
     lagged state until an update, or until interfaceExchange restores one
     (restoreLastUpdate()) */
  states_.resize(nElements*stateSize_, 0.0);
  valid_.resize(nElements, 0);
  lagged_ = false;

#else
  FatalIOErrorInFunction(dict_) << "This executable was built without the "
    << "GEMS3K bridge." << exit(FatalIOError);
#endif
}


/*------------------------------------------------------------------------------
frozen: the GEMS table of every pair.  The master evaluates p_eq at the
nodes of the pair's table within [max(minTemperature, T_min), T_max] (cold
calls, the reference composition, referencePressure); a failed or rejected
node takes the table's value.  GEMSB_ERR_FATAL re-creates the engine, as
in mode local (the bridge's contract: the engine must be re-created; the
condensed species are suppressed again by createEngine(), and the calls
are cold, so there is no warm state to invalidate; plan section 35, item
E).  Every rank receives the values.
------------------------------------------------------------------------------*/

void LESTO::gemsEquilibrium::buildFrozenTables (
  const UList<const vapourPressureTable*>& tables
) {
#ifdef LESTO_HAVE_GEMS
  frozen_.resize(pairs_.size());
  frozenReport_.resize(pairs_.size());

  scalarList bReference(nElements_, 0.0);
  forAll(referenceNames_, n) {
    bReference[elementIndex(referenceNames_[n], dict_, "referenceComposition")]
      += referenceAmounts_[n];
  }

  forAll(pairs_, k) {

    const pairSpec& p = pairs_[k];
    const vapourPressureTable& table = *tables[k];

    std::vector<double> T;
    for (const double t : table.nodes()) {
      if (t >= minTemperature_ && t >= Tmin_ && t <= Tmax_) {
        T.push_back(t);
      }
    }
    const label n = label(T.size());

    if (n < 2) {
      FatalIOErrorInFunction(dict_) << "Pair " << p.gasName << ": the table "
        << table.file() << " has " << n << " nodes in [max(minTemperature "
        << minTemperature_ << ", " << Tmin_ << "), " << Tmax_ << "] K, the "
        << "range of the GEMS table of mode frozen; it needs two."
        << exit(FatalIOError);
    }

    /* log10 p_eq [Pa] at the nodes, and the source of each */
    scalarList lgp(n, 0.0);
    labelList source(n, label(FROM_GEMS));
    label iterations = 0;
    label recreated = 0;

    if (Pstream::master()) {
      gemsb_engine* h = engine_->handle();
      for (label i = 0; i < n; ++i) {
        const scalar lgTable = std::log10(table(T[i]));
        b_ = bReference;
        const int rc = gemsb_equilibrate(h, T[i], referencePressure_,
                                         b_.cdata(), GEMSB_COLD, nullptr,
                                         nullptr);
        iterations += gemsb_last_iterations(h);
        bool ok = rc == GEMSB_OK || rc == GEMSB_OK_RETRIED;
        scalar lgOmega = -VGREAT;
        std::string why(statusName(rc));
        if (ok) {
          gemsb_species_amounts(h, amounts_.data());
          gemsb_gas_partial_pressures(h, pGas_.data());
          gemsb_log10_activities(h, log10a_.data());
          /* log10 Omega first, as evaluate() does: the line of a node with
             p_g <= 0 names its real max log10 Omega, not the start value
             -VGREAT (plan section 35, item Y) */
          const bool finite =
            maxLog10Omega(p.jCondensates, log10a_, lgOmega);
          if (!balanced()) {
            ok = false;
            why = "element balance off by more than 1e-6";
          } else if (!positive(pGas_[p.jGas]) || !finite) {
            ok = false;
            why = "converged, but p_g " + text(pGas_[p.jGas])
                + " Pa, max log10 Omega " + text(lgOmega);
          }
        }
        if (!ok) {
          source[i] = FAILED;
          lgp[i] = lgTable;
          logLine("frozen " + p.gasName + ": the call at T = " + text(T[i])
                + " K failed (" + why + "): " + gemsb_last_error(h));
          if (rc == GEMSB_ERR_FATAL) {
            string message;
            if (!createEngine(message)) {
              FatalErrorInFunction << "The GEMS3K engine could not be "
                << "re-created after GEMSB_ERR_FATAL from " << system_
                << ": " << message.c_str() << nl << abort(FatalError);
            }
            h = engine_->handle();
            ++recreated;
            logLine("frozen " + p.gasName + ": engine re-created after "
                  + "GEMSB_ERR_FATAL; condensed species suppressed again");
          }
          continue;
        }
        scalar lg = std::log10(pGas_[p.jGas]) - lgOmega;
        if(dualPairs_) {
          std::vector<int> condensed;for(label j:p.jCondensates)condensed.push_back(int(j));
          double pressure=0;int best=-1;
          if(gemsb_pair_peq(h,int(p.jGas),condensed.data(),int(condensed.size()),1,&pressure,&best)!=GEMSB_OK || !positive(pressure)) {
            source[i]=FAILED;lgp[i]=lgTable;continue;
          }
          lg=std::log10(pressure);
        }
        if (!(mag(lg - lgTable) <= maxLog10Deviation_)) {
          source[i] = REJECTED;
          lgp[i] = lgTable;
          logLine("frozen " + p.gasName + ": T = " + text(T[i]) + " K "
                + "rejected, log10 p_eq " + text(lg) + " (table "
                + text(lgTable) + ")");
        } else {
          lgp[i] = lg;
        }
      }
    }

    Pstream::broadcast(lgp);
    Pstream::broadcast(source);
    Pstream::broadcast(iterations);
    Pstream::broadcast(recreated);

    try {
      frozen_.set (
        k,
        new vapourPressureTable (
          T,
          std::vector<std::vector<double>>(1, std::vector<double>(
            lgp.cbegin(), lgp.cend())),
          std::vector<std::string>(1, "p_eq"),
          vapourPressureTable::LOG10_PA,
          table.interpolationMode(),
          vapourPressureTable::EXTRAPOLATE,
          "the GEMS table of pair " + p.gasName
        )
      );
    } catch (const std::exception& error) {
      FatalIOErrorInFunction(dict_) << "Pair " << p.gasName << ": "
        << error.what() << exit(FatalIOError);
    }

    /* report: nodes, sources, the largest deviation of the GEMS nodes */
    label nGems = 0, nFailed = 0, nRejected = 0;
    scalar maxDeviation = 0;
    for (label i = 0; i < n; ++i) {
      if (source[i] == FROM_GEMS) {
        ++nGems;
        maxDeviation =
          max(maxDeviation, mag(lgp[i] - std::log10(table(T[i]))));
      } else if (source[i] == FAILED) {
        ++nFailed;
      } else {
        ++nRejected;
      }
    }

    OStringStream os;
    os << "Pair " << p.gasName << ", GEMS frozen: p_eq at the " << n
      << " nodes of the table in " << T.front() << "-" << T.back() << " K: "
      << nGems << " from GEMS3K (" << scalar(iterations)/n
      << " iterations per cold call), " << nFailed << " failed, "
      << nRejected << " rejected (taken from the table), engine re-created "
      << recreated << (recreated == 1 ? " time" : " times")
      << "; max |dlog10 p_eq - table| " << maxDeviation << " at the GEMS "
      << "nodes; an element "
      << "below " << T.front() << " K or above " << T.back() << " K takes "
      << "the table";
    frozenReport_[k] = os.str();
    Info<< frozenReport_[k].c_str() << nl;
  }
#endif
}


/*------------------------------------------------------------------------------
Element balance of the last call.  The sums of the amounts are compared only
once they are known to be finite (a NaN amount is unbalanced, and never
enters an ordered comparison: maxLog10Omega()).
------------------------------------------------------------------------------*/

bool LESTO::gemsEquilibrium::balanced() const {
  scalar total = 0;
  for (const scalar b : b_) {
    total += max(b, scalar(0));
  }
  for (label i = 0; i < nElements_; ++i) {
    if (!(b_[i] > 1e-12*total)) {
      continue;
    }
    scalar s = 0;
    for (label j = 0; j < nSpecies_; ++j) {
      s += stoichiometry_[j*nElements_ + i]*amounts_[j];
    }
    if (!std::isfinite(s) || !(mag(s - b_[i]) <= 1e-6*b_[i])) {
      return false;
    }
  }
  return true;
}


void LESTO::gemsEquilibrium::logLine(const std::string& line) {
  if (log_.is_open()) {
    log_ << line << '\n';
    log_.flush();
  }
}


/*------------------------------------------------------------------------------
frozen: p_eq of an element
------------------------------------------------------------------------------*/

Foam::scalar LESTO::gemsEquilibrium::frozenPEq (
  const label  k,
  const scalar T,
  const scalar tableValue,
  label&       source
) const {
  if (!(T >= minTemperature_)) {
    source = BELOW_TEMPERATURE;
    return tableValue;
  }
  const vapourPressureTable& gems = frozen_[k];
  if (!(T >= gems.Tmin() && T <= gems.Tmax())) {
    source = OUTSIDE_GRID;
    return tableValue;
  }
  source = FROM_GEMS;
  return gems(T);
}


/*------------------------------------------------------------------------------
local: the update of a step
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
The schedule follows the time index n of the step (Time::timeIndex(), which
a restart takes from uniform/time): updates at n = 1, 1 + N, 1 + 2N, ...
(N = updateInterval), so a restart keeps the steps of the continuous run
(plan section 35, item C; before, every run updated at its own first step
and counted from there).  A run whose elements hold no lagged state (a
fresh start; a restart that could not restore it) also updates at its
first step.
------------------------------------------------------------------------------*/

bool LESTO::gemsEquilibrium::updateDueAt(const label timeIndex) const {
  const label n = timeIndex - 1;
  return n < 0 || n % updateInterval_ == 0;
}


bool LESTO::gemsEquilibrium::updateDue(const label timeIndex) {
  updated_ = !lagged_ || updateDueAt(timeIndex);
  return updated_;
}


/*------------------------------------------------------------------------------
The counts of the last update per pair, for the per-rank restart state
(interfaceExchange: uniform/phaseChange/phaseChangeEquilibrium), and back
------------------------------------------------------------------------------*/

Foam::scalarList LESTO::gemsEquilibrium::lastUpdate() const {
  scalarList values(pairUpdate_.size()*nPairState);
  forAll(pairUpdate_, k) {
    const pairCounts& u = pairUpdate_[k];
    scalar* v = values.data() + k*nPairState;
    for (label s = 0; s < nSourceTypes; ++s) {
      v[s] = scalar(u.sources[s]);
    }
    v[nSourceTypes] = u.maxDeviation;
    v[nSourceTypes + 1] = u.chiMin;
    v[nSourceTypes + 2] = u.chiMax;
    v[nSourceTypes + 3] = scalar(u.chiClipped);
  }
  return values;
}


void LESTO::gemsEquilibrium::restoreLastUpdate(const UList<scalar>& values) {
  forAll(pairUpdate_, k) {
    pairCounts& u = pairUpdate_[k];
    const scalar* v = values.cdata() + k*nPairState;
    for (label s = 0; s < nSourceTypes; ++s) {
      u.sources[s] = std::int64_t(v[s]);
    }
    u.maxDeviation = v[nSourceTypes];
    u.chiMin = v[nSourceTypes + 1];
    u.chiMax = v[nSourceTypes + 2];
    u.chiClipped = std::int64_t(v[nSourceTypes + 3]);
  }
  update_ = noCounts();
  updated_ = false;
  lagged_ = true;
}


void LESTO::gemsEquilibrium::beginUpdate() {
  update_ = noCounts();
  update_.updates = 1;
  for (pairCounts& c : pairUpdate_) {
    c = noPairCounts();
  }
  failuresLogged_ = 0;
}


void LESTO::gemsEquilibrium::evaluate (
  const label          e,
  const scalar         T,
  const scalar         p,
  const scalar         cCarrier,
  const UList<scalar>& c,
  const UList<bool>&   active,
  const UList<scalar>& pEqTable,
  UList<scalar>&       pEq,
  UList<scalar>&       chi,
  UList<label>&        source
) {
#ifdef LESTO_HAVE_GEMS

  /* the table for pair k, with the reason */
  auto useTable = [&](const label k, const label reason) {
    pEq[k] = pEqTable[k];
    chi[k] = 1;
    source[k] = reason;
    ++pairUpdate_[k].sources[reason];
  };

  /*--------------------------------------------------------------------------
  Guards before the call: the temperature, the grid of the system, the
  mole fraction of every pair (a pair without gas, c <= 0, is below it).
  One call serves every pair of the element that passes them.
  --------------------------------------------------------------------------*/
  label guard = FROM_GEMS;
  if (!(T >= minTemperature_)) {
    guard = BELOW_TEMPERATURE;
  } else if (!(T >= Tmin_ && T <= Tmax_ && p >= Pmin_ && p <= Pmax_)) {
    guard = OUTSIDE_GRID;
  }

  bool call = false;
  forAll(pairs_, k) {
    if (!active[k]) {
      continue;
    }
    label reason = guard;
    if (reason == FROM_GEMS
      && !(c[k] > 0 && c[k]/(c[k] + cCarrier) >= minMoleFraction_)) {
      reason = BELOW_MOLE_FRACTION;
    }
    if(reason==FROM_GEMS && dualPairs_ && c[k]/(c[k]+cCarrier)>maxMoleFraction_)reason=REJECTED;
    if (reason == FROM_GEMS) {
      source[k] = FROM_GEMS;
      call = true;
    } else {
      useTable(k, reason);
    }
  }
  if (!call) {
    return;
  }

  /*--------------------------------------------------------------------------
  The bulk composition of the cell [mol/m3]: the carrier and the formula
  units of every pair present (with their excess).
  --------------------------------------------------------------------------*/
  for (label i = 0; i < nElements_; ++i) {
    b_[i] = cCarrier*carrierFormula_[i];
  }
  forAll(pairs_, k) {
    if (c[k] > 0) {
      const scalarList& bulk = pairs_[k].bulk;
      for (label i = 0; i < nElements_; ++i) {
        b_[i] += c[k]*bulk[i];
      }
    }
  }

  gemsb_engine* h = engine_->handle();
  int& valid = valid_[e];
  const bool warmStart = valid != 0;
  double* state = states_.data() + e*stateSize_;

  const auto t0 = std::chrono::steady_clock::now();
  const int rc = gemsb_equilibrate(h, T, p, b_.cdata(), GEMSB_WARM_CELL,
                                   state, &valid);
  const bool converged = rc == GEMSB_OK || rc == GEMSB_OK_RETRIED;
  if (converged) {
    gemsb_species_amounts(h, amounts_.data());
    gemsb_gas_partial_pressures(h, pGas_.data());
    gemsb_log10_activities(h, log10a_.data());
  }
  const scalar seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - t0).count();

  ++update_.calls;
  update_.iterations += gemsb_last_iterations(h);
  if (rc == GEMSB_OK && warmStart) {
    ++update_.warm;
    update_.warmSeconds += seconds;
  } else if (rc == GEMSB_OK) {
    ++update_.cold;
    update_.coldSeconds += seconds;
  } else if (rc == GEMSB_OK_RETRIED) {
    ++update_.retried;
    update_.coldSeconds += seconds;
  }

  /*--------------------------------------------------------------------------
  A failure: a status other than OK and OK_RETRIED, or an element balance
  off by more than 1e-6.  The table for every pair that called, a cold
  start next time, and after GEMSB_ERR_FATAL a new engine (bounds
  re-applied, every warm state invalidated).
  --------------------------------------------------------------------------*/
  if (!converged || !balanced()) {
    ++update_.failed;
    valid = 0;
    forAll(pairs_, k) {
      if (active[k] && source[k] == FROM_GEMS) {
        useTable(k, FAILED);
      }
    }
    if (failuresLogged_ < 20) {
      ++failuresLogged_;
      std::string fractions;
      forAll(pairs_, k) {
        fractions += " x(" + pairs_[k].gasName + ") "
          + text(c[k]/(c[k] + cCarrier));
      }
      logLine("time " + mesh_.time().timeName() + ": element " + text(e)
            + ", T " + text(T) + " K, p " + text(p) + " Pa," + fractions
            + ": " + (converged ? "element balance off by more than 1e-6"
                                : statusName(rc))
            + ", " + text(gemsb_last_iterations(h)) + " iterations: "
            + gemsb_last_error(h));
    }
    if (rc == GEMSB_ERR_FATAL) {
      string message;
      if (!createEngine(message)) {
        FatalErrorInFunction << "Rank " << Pstream::myProcNo() << ": the "
          << "GEMS3K engine could not be re-created after GEMSB_ERR_FATAL "
          << "from " << system_ << ": " << message.c_str() << nl
          << abort(FatalError);
      }
      valid_ = 0;
      ++update_.recreated;
      logLine("time " + mesh_.time().timeName() + ": engine re-created "
            + "after GEMSB_ERR_FATAL; condensed species suppressed again, "
            + text(valid_.size()) + " warm states invalidated");
    }
    return;
  }

  /*--------------------------------------------------------------------------
  p_eq = p_g 10^(-max log10 Omega) of every pair that called, validated
  against the table; the speciation factor chi = n_g/c.  A converged and
  balanced call that gives a pair p_g <= 0 or a max log10 Omega that is not
  finite is a failed call as well: the table, a cold start next time, the
  call counted and written to the log of the rank (plan section 35, item
  F).
  --------------------------------------------------------------------------*/
  std::string unusable;
  forAll(pairs_, k) {
    if (!active[k] || source[k] != FROM_GEMS) {
      continue;
    }
    const pairSpec& ps = pairs_[k];
    scalar lgOmega = -VGREAT;
    const bool finite = maxLog10Omega(ps.jCondensates, log10a_, lgOmega);
    const scalar pg = pGas_[ps.jGas];
    if (!positive(pg) || !finite) {
      useTable(k, FAILED);
      valid = 0;
      unusable += " " + ps.gasName + " p_g " + text(pg)
                + " Pa, max log10 Omega " + text(lgOmega) + ";";
      continue;
    }
    scalar lg = std::log10(pg) - lgOmega;
    if(dualPairs_) {
      std::vector<int> condensed;for(label j:ps.jCondensates)condensed.push_back(int(j));
      double pressure=0;int best=-1;
      if(gemsb_pair_peq(h,int(ps.jGas),condensed.data(),int(condensed.size()),1,&pressure,&best)!=GEMSB_OK || !positive(pressure)) {
        useTable(k,FAILED);valid=0;unusable+=" "+ps.gasName+" unusable dual pressure;";continue;
      }
      lg=std::log10(pressure);
      if(speciation_ && pg/p<1e-7){useTable(k,REJECTED);valid=0;continue;}
    }
    const scalar deviation = mag(lg - std::log10(pEqTable[k]));
    if (!(deviation <= maxLog10Deviation_)) {
      useTable(k, REJECTED);
      valid = 0;
      continue;
    }

    pairCounts& counted = pairUpdate_[k];
    pEq[k] = std::pow(10.0, lg);
    chi[k] = 1;
    if (speciation_) {
      const scalar x = amounts_[ps.jGas]/c[k];
      chi[k] = min(max(x, scalar(chiMin)), scalar(1));
      counted.chiClipped += chi[k] != x;
      counted.chiMin = min(counted.chiMin, chi[k]);
      counted.chiMax = max(counted.chiMax, chi[k]);
    }
    ++counted.sources[FROM_GEMS];
    counted.maxDeviation = max(counted.maxDeviation, deviation);
  }

  if (!unusable.empty()) {
    ++update_.failed;
    if (failuresLogged_ < 20) {
      ++failuresLogged_;
      std::string fractions;
      forAll(pairs_, k) {
        fractions += " x(" + pairs_[k].gasName + ") "
          + text(c[k]/(c[k] + cCarrier));
      }
      logLine("time " + mesh_.time().timeName() + ": element " + text(e)
            + ", T " + text(T) + " K, p " + text(p) + " Pa," + fractions
            + ": converged and balanced (" + statusName(rc) + ", "
            + text(gemsb_last_iterations(h)) + " iterations), but"
            + unusable + " the table");
    }
  }
#endif
}


void LESTO::gemsEquilibrium::endUpdate() {

  /* the log of this rank: its own counts */
  {
    OStringStream os;
    os << "time " << mesh_.time().timeName() << " (index "
      << mesh_.time().timeIndex() << "): calls "
      << num(update_.calls).c_str() << " (warm " << num(update_.warm).c_str()
      << ", retried " << num(update_.retried).c_str() << ", cold "
      << num(update_.cold).c_str() << "), failures "
      << num(update_.failed).c_str() << ", re-created "
      << num(update_.recreated).c_str() << ", iterations "
      << num(update_.iterations).c_str();
    if (update_.warm) {
      os << ", " << 1e6*update_.warmSeconds/scalar(update_.warm)
        << " us per warm call";
    }
    forAll(pairs_, k) {
      os << "; " << pairs_[k].gasName << " sources "
        << nums(pairUpdate_[k].sources).c_str();
    }
    logLine(os.str());
  }

  /*--------------------------------------------------------------------------
  The counts over all ranks: every integer count in one reduction of 64-bit
  integers (a label would overflow in a long run: plan section 35, item
  A), the seconds in one of doubles, the extremes of the pairs in one of
  their maxima (chiMin negated).
  --------------------------------------------------------------------------*/
  const label nPairs = pairUpdate_.size();
  List<std::int64_t> n(7 + nPairs*(nSourceTypes + 1));
  {
    label i = 0;
    n[i++] = update_.calls;
    n[i++] = update_.warm;
    n[i++] = update_.cold;
    n[i++] = update_.retried;
    n[i++] = update_.failed;
    n[i++] = update_.recreated;
    n[i++] = update_.iterations;
    for (const pairCounts& u : pairUpdate_) {
      for (const std::int64_t m : u.sources) {
        n[i++] = m;
      }
      n[i++] = u.chiClipped;
    }
  }
  reduce(n.data(), int(n.size()), sumOp<std::int64_t>(), UPstream::msgType(),
         UPstream::worldComm);

  FixedList<scalar, 2> seconds({update_.warmSeconds, update_.coldSeconds});
  reduce(seconds, sumOp<scalar>());

  scalarList extremes(3*nPairs);
  forAll(pairUpdate_, k) {
    extremes[3*k] = pairUpdate_[k].maxDeviation;
    extremes[3*k + 1] = pairUpdate_[k].chiMax;
    extremes[3*k + 2] = -pairUpdate_[k].chiMin;
  }
  reduce(extremes.data(), int(extremes.size()), maxOp<scalar>(),
         UPstream::msgType(), UPstream::worldComm);

  {
    label i = 0;
    update_.calls = n[i++];
    update_.warm = n[i++];
    update_.cold = n[i++];
    update_.retried = n[i++];
    update_.failed = n[i++];
    update_.recreated = n[i++];
    update_.iterations = n[i++];
    for (pairCounts& u : pairUpdate_) {
      for (std::int64_t& m : u.sources) {
        m = n[i++];
      }
      u.chiClipped = n[i++];
    }
  }
  update_.warmSeconds = seconds[0];
  update_.coldSeconds = seconds[1];
  forAll(pairUpdate_, k) {
    pairUpdate_[k].maxDeviation = extremes[3*k];
    pairUpdate_[k].chiMax = extremes[3*k + 1];
    pairUpdate_[k].chiMin = -extremes[3*k + 2];
  }

  run_.calls += update_.calls;
  run_.warm += update_.warm;
  run_.cold += update_.cold;
  run_.retried += update_.retried;
  run_.failed += update_.failed;
  run_.recreated += update_.recreated;
  run_.iterations += update_.iterations;
  run_.warmSeconds += update_.warmSeconds;
  run_.coldSeconds += update_.coldSeconds;
  run_.updates += 1;

  forAll(pairUpdate_, k) {
    const pairCounts& u = pairUpdate_[k];
    pairCounts& r = pairRun_[k];
    forAll(r.sources, s) {
      r.sources[s] += u.sources[s];
    }
    r.maxDeviation = max(r.maxDeviation, u.maxDeviation);
    r.chiMin = min(r.chiMin, u.chiMin);
    r.chiMax = max(r.chiMax, u.chiMax);
    r.chiClipped += u.chiClipped;
  }

  /* the elements hold the p_eq of this update (the lagged state) */
  lagged_ = true;
}


/*------------------------------------------------------------------------------
Reports
------------------------------------------------------------------------------*/

Foam::string LESTO::gemsEquilibrium::stepReport(const label k) const {

  const pairCounts& u = pairUpdate_[k];
  OStringStream os;

  if (updated_) {
    os << "GEMS local: update (calls " << num(update_.calls).c_str()
      << ": warm " << num(update_.warm).c_str() << ", retried "
      << num(update_.retried).c_str() << ", cold "
      << num(update_.cold).c_str() << "; "
      << scalar(update_.iterations)
        /scalar(std::max(update_.calls, std::int64_t(1)))
      << " iterations per call; failures " << num(update_.failed).c_str()
      << ", engine re-created " << num(update_.recreated).c_str()
      << "); elements: ";
  } else {
    os << "GEMS local: no update this step (updateInterval "
      << updateInterval_ << "); elements of the last update: ";
  }

  std::int64_t nTable = 0;
  for (label s = 1; s < nSourceTypes; ++s) {
    nTable += u.sources[s];
  }
  os << "GEMS " << num(u.sources[FROM_GEMS]).c_str() << ", table "
    << num(nTable).c_str() << " (below minTemperature "
    << num(u.sources[BELOW_TEMPERATURE]).c_str() << ", below "
    << "minMoleFraction " << num(u.sources[BELOW_MOLE_FRACTION]).c_str()
    << ", outside the GEMS grid " << num(u.sources[OUTSIDE_GRID]).c_str()
    << ", rejected " << num(u.sources[REJECTED]).c_str() << ", failed "
    << num(u.sources[FAILED]).c_str() << "); max |dlog10 p_eq/p_eq,table| "
    << u.maxDeviation;
  if (speciation_) {
    os << "; speciation lagged: chi ";
    if (u.sources[FROM_GEMS] > 0) {
      os << u.chiMin << "-" << u.chiMax << " (clipped "
        << num(u.chiClipped).c_str() << ")";
    } else {
      os << "1 (no GEMS element)";
    }
  }
  return os.str();
}


void LESTO::gemsEquilibrium::summary() const {

  if (!local()) {
    return;
  }

  Info<< "GEMS3K summary (mode local): " << num(run_.updates).c_str()
    << " updates, " << num(run_.calls).c_str() << " calls (warm "
    << num(run_.warm).c_str() << ", retried " << num(run_.retried).c_str()
    << ", cold " << num(run_.cold).c_str() << "), "
    << num(run_.failed).c_str() << " failures, engine re-created "
    << num(run_.recreated).c_str() << " times; "
    << scalar(run_.iterations)/scalar(std::max(run_.calls, std::int64_t(1)))
    << " iterations per call; ";
  if (run_.warm > 0) {
    Info<< 1e6*run_.warmSeconds/scalar(run_.warm) << " us per warm call";
  } else {
    Info<< "no warm call";
  }
  if (run_.cold + run_.retried > 0) {
    Info<< ", " << 1e6*run_.coldSeconds/scalar(run_.cold + run_.retried)
      << " us per cold or retried call";
  }
  Info<< " (all ranks)" << nl;

  forAll(pairs_, k) {
    const pairCounts& r = pairRun_[k];
    std::int64_t nTable = 0;
    for (label s = 1; s < nSourceTypes; ++s) {
      nTable += r.sources[s];
    }
    Info<< "GEMS3K summary " << pairs_[k].gasName << ": element-updates "
      << "GEMS " << num(r.sources[FROM_GEMS]).c_str() << ", table "
      << num(nTable).c_str() << " (below minTemperature "
      << num(r.sources[BELOW_TEMPERATURE]).c_str() << ", below "
      << "minMoleFraction " << num(r.sources[BELOW_MOLE_FRACTION]).c_str()
      << ", outside the GEMS grid " << num(r.sources[OUTSIDE_GRID]).c_str()
      << ", rejected " << num(r.sources[REJECTED]).c_str() << ", failed "
      << num(r.sources[FAILED]).c_str() << "); max |dlog10 p_eq/p_eq,table| "
      << r.maxDeviation;
    if (speciation_ && r.sources[FROM_GEMS] > 0) {
      Info<< "; chi " << r.chiMin << "-" << r.chiMax << " (clipped "
        << num(r.chiClipped).c_str() << ")";
    }
    Info<< nl;
  }
}

/* M10c homogeneous gas backend; no bridge symbols in the default build. */
#include "gemsGasSpeciation.H"
#include <fstream>
#include <sstream>
#include <iomanip>

struct LESTO::gemsGasSpeciation::implementation {
  fileName system;
  word mode;
  string fingerprint;
  scalar minimum=1e-12, maximum=3.2e-4;
  label warmLimit=200;
#ifdef LESTO_HAVE_GEMS
  gemsb_engine* h=nullptr;
  std::vector<int> gases, atoms;
  int carrier=-1, nElements=0, nSpecies=0, stateSize=0;
  std::vector<double> states, bulk, potentials;
  std::vector<int> valid;
  ~implementation() { if(h)gemsb_destroy(h); }
  void create() {
    char message[2048]={0};
    h=gemsb_create_from_lst(system.c_str(),message,sizeof(message));
    if(!h)FatalErrorInFunction<<"Cannot create gas GEMS engine: "<<message<<exit(FatalError);
    gemsb_suppress_condensed(h,1);gemsb_set_warm_iteration_limit(h,int(warmLimit));
  }
#endif
};

LESTO::gemsGasSpeciation::gemsGasSpeciation(const fvMesh& mesh,const dictionary& cfg,
  const wordList& names,const std::vector<gasSpecies>& formulas,const scalarList& masses)
: impl_(new implementation) {
  gemsEquilibrium::requireBridge(cfg);
  auto& d=*impl_;
  const wordList known({"system","mode","carrier","gases","minMoleFraction","maxMoleFraction","warmIterationLimit","logDirectory","logLevel",
    "referenceComposition","referencePressure","minTemperature","updateInterval","allowNonPhysicalCarrier","writeEquilibriumField","maxLog10Deviation"});
  for(const entry& e:cfg)if(!known.found(e.keyword()))
    FatalIOErrorInFunction(cfg)<<"Unknown gas GEMS key "<<e.keyword()<<exit(FatalIOError);
  d.mode=cfg.getOrDefault<word>("mode","frozen");
  d.minimum=cfg.getOrDefault<scalar>("minMoleFraction",1e-12);
  d.maximum=cfg.getOrDefault<scalar>("maxMoleFraction",3.2e-4);
  d.warmLimit=cfg.getOrDefault<label>("warmIterationLimit",200);
  if((d.mode!="frozen" && d.mode!="local") || !std::isfinite(d.minimum)
    || !std::isfinite(d.maximum) || d.minimum<0 || d.maximum<=d.minimum
    || d.maximum>1 || d.warmLimit<1 || d.warmLimit>32000)
    FatalIOErrorInFunction(cfg)<<"Invalid gas GEMS mode, mole-fraction guards or warm iteration limit"<<exit(FatalIOError);
  d.system=cfg.get<fileName>("system");
  d.system.replaceAll("<constant>",mesh.time().globalPath()/mesh.time().constant());d.system.expand();
  // Include every document named by the export list in the restart signature.
  // Both JSON and key-value exports use quoted relative document filenames.
  std::ifstream list(d.system.c_str());if(!list)FatalIOErrorInFunction(cfg)<<"Cannot read "<<d.system<<exit(FatalIOError);
  std::ostringstream data;data<<list.rdbuf();const std::string lst=data.str();
  std::istringstream tokens(lst);std::string token;
  while(tokens>>std::quoted(token)) {
    if(token.empty() || token[0]=='-')continue;
    const fileName path=fileName(token).isAbsolute()?fileName(token):d.system.path()/fileName(token);
    std::ifstream input(path.c_str());if(!input)FatalIOErrorInFunction(cfg)<<"Cannot read GEMS document "<<path<<exit(FatalIOError);
    data<<'\0'<<input.rdbuf();
  }
  const word carrier=cfg.getOrDefault<word>("carrier","He(g)");
  data<<'\0'<<d.mode<<'\0'<<std::setprecision(17)<<d.minimum<<' '<<d.maximum<<' '<<d.warmLimit<<' '<<carrier;
#ifdef LESTO_HAVE_GEMS
  fileName log=mesh.time().path()/cfg.getOrDefault<fileName>("logDirectory","gasGemsLog");
  mkDir(log);gemsb_set_log_directory(log.c_str());gemsb_set_log_level(cfg.getOrDefault<label>("logLevel",4));
  d.create();d.nElements=gemsb_num_elements(d.h);d.nSpecies=gemsb_num_species(d.h);
  d.carrier=gemsb_species_index(d.h,carrier.c_str());
  if(d.carrier<0 || !gemsb_species_is_gas(d.h,d.carrier))
    FatalIOErrorInFunction(cfg)<<"GEMS carrier must be an inert gas species"<<exit(FatalIOError);
  d.atoms={gemsb_element_index(d.h,"Pb"),gemsb_element_index(d.h,"Bi"),gemsb_element_index(d.h,"I")};
  for(int a:d.atoms)if(a>=0 && gemsb_stoich(d.h,d.carrier,a)!=0)
    FatalIOErrorInFunction(cfg)<<"Carrier contains a reactive element"<<exit(FatalIOError);
  forAll(names,i) {
    word gas;
    if(cfg.isDict("gases"))gas=cfg.subDict("gases").get<word>(names[i]);
    else { gas=names[i];if(gas.size()>2 && gas.substr(gas.size()-2)=="_g")gas=gas.substr(0,gas.size()-2)+"(g)"; }
    const int j=gemsb_species_index(d.h,gas.c_str());
    if(j<0 || !gemsb_species_is_gas(d.h,j) || std::find(d.gases.begin(),d.gases.end(),j)!=d.gases.end())
      FatalIOErrorInFunction(cfg)<<"Invalid or duplicate mapped GEMS gas "<<gas<<exit(FatalIOError);
    const auto& f=formulas[i];const int nu[]={f.nPb,f.nBi,f.nI};
    for(int e=0;e<d.nElements;++e) {
      int expected=0;for(int a=0;a<3;++a)if(d.atoms[a]==e)expected=nu[a];
      if(gemsb_stoich(d.h,j,e)!=expected)FatalIOErrorInFunction(cfg)<<"GEMS formula differs for "<<gas<<exit(FatalIOError);
    }
    for(int a=0;a<3;++a)if(nu[a] && d.atoms[a]<0)
      FatalIOErrorInFunction(cfg)<<"Missing GEMS element for "<<gas<<exit(FatalIOError);
    if(std::fabs(gemsb_species_molar_mass(d.h,j)/masses[i]-1)>1e-6)
      FatalIOErrorInFunction(cfg)<<"GEMS molar mass differs for "<<gas<<exit(FatalIOError);
    d.gases.push_back(j);data<<'\0'<<gas;
  }
  // An omitted gas with reactive atoms would make the two engines solve
  // different chemical systems. Refuse it instead of losing its inventory.
  for(int j=0;j<d.nSpecies;++j)if(gemsb_species_is_gas(d.h,j) && j!=d.carrier
      && std::find(d.gases.begin(),d.gases.end(),j)==d.gases.end())
    FatalIOErrorInFunction(cfg)<<"Unmapped GEMS gas "<<gemsb_species_name(d.h,j)<<exit(FatalIOError);
  d.stateSize=gemsb_state_size(d.h);d.bulk.resize(d.nElements);d.potentials.resize(d.nElements);
  if(local()){d.states.resize(mesh.nCells()*d.stateSize);d.valid.resize(mesh.nCells(),0);}
#else
  (void)names;(void)formulas;(void)masses;
#endif
  std::uint64_t hash=14695981039346656037ULL;for(unsigned char c:data.str()){hash^=c;hash*=1099511628211ULL;}
  d.fingerprint=std::to_string(hash);
  Info<<"Gas GEMS: mode "<<d.mode<<"; system "<<d.system<<"; minMoleFraction "<<d.minimum
      <<"; maxMoleFraction "<<d.maximum<<"; warmIterationLimit "<<d.warmLimit<<nl;
}
LESTO::gemsGasSpeciation::~gemsGasSpeciation()=default;
bool LESTO::gemsGasSpeciation::local()const{return impl_->mode=="local";}
string LESTO::gemsGasSpeciation::signature()const{return impl_->fingerprint;}
bool LESTO::gemsGasSpeciation::constants(double T,double P,std::vector<double>& lnK)const {
#ifdef LESTO_HAVE_GEMS
  auto& d=*impl_;std::vector<double> g(d.nSpecies);
  double tmin,tmax,pmin,pmax;gemsb_TP_range(d.h,&tmin,&tmax,&pmin,&pmax);
  if(T<tmin || T>tmax || P<pmin || P>pmax)return false;
  if(gemsb_species_g0(d.h,T,P,g.data())!=GEMSB_OK)return false;
  lnK.resize(d.gases.size());
  for(std::size_t i=0;i<d.gases.size();++i){double k=-g[d.gases[i]];
    for(int a:d.atoms)if(a>=0){const char* element=gemsb_element_name(d.h,a);
      const std::string atom=std::string(element)+"(g)";
      const int j=gemsb_species_index(d.h,atom.c_str());
      const double nu=gemsb_stoich(d.h,d.gases[i],a);
      if(nu && j<0)return false;
      if(nu)k+=nu*g[j];}
    if(!std::isfinite(k))return false;
    lnK[i]=k;}
  return true;
#else
  (void)T;(void)P;(void)lnK;return false;
#endif
}
bool LESTO::gemsGasSpeciation::warm(label cell,double T,double P,
  const std::array<double,3>& b,double trace,std::array<double,3>& pi,
  int& iterations,bool& guarded) {
  iterations=0;guarded=false;
#ifdef LESTO_HAVE_GEMS
  auto& d=*impl_;const double carrier=P/(8.314462618*T);
  const double x=trace/(carrier+trace);
  double tmin,tmax,pmin,pmax;gemsb_TP_range(d.h,&tmin,&tmax,&pmin,&pmax);
  if(x<d.minimum || x>d.maximum || T<tmin || T>tmax || P<pmin || P>pmax){guarded=true;return false;}
  std::fill(d.bulk.begin(),d.bulk.end(),0);
  for(int e=0;e<d.nElements;++e)d.bulk[e]=carrier*gemsb_stoich(d.h,d.carrier,e);
  for(int a=0;a<3;++a)if(d.atoms[a]>=0)d.bulk[d.atoms[a]]+=b[a];
  int& valid=d.valid[cell];
  const int rc=gemsb_equilibrate(d.h,T,P,d.bulk.data(),GEMSB_WARM_CELL,
                                d.states.data()+cell*d.stateSize,&valid);
  iterations=gemsb_last_iterations(d.h);
  if((rc!=GEMSB_OK && rc!=GEMSB_OK_RETRIED) || gemsb_balance_error(d.h,d.bulk.data(),1e-12)>1e-6){
    valid=0;if(rc==GEMSB_ERR_FATAL){gemsb_destroy(d.h);d.h=nullptr;d.create();std::fill(d.valid.begin(),d.valid.end(),0);}return false;}
  gemsb_element_potentials(d.h,d.potentials.data());
  std::vector<double> g(d.nSpecies);if(gemsb_species_g0(d.h,T,P,g.data())!=GEMSB_OK){valid=0;return false;}
  const char* atoms[]={"Pb(g)","Bi(g)","I(g)"};
  for(int a=0;a<3;++a)if(b[a]>0 && d.atoms[a]>=0){const int j=gemsb_species_index(d.h,atoms[a]);
    if(j<0 || !std::isfinite(d.potentials[d.atoms[a]]-g[j])){valid=0;return false;}
    pi[a]=d.potentials[d.atoms[a]]-g[j];}
  return true;
#else
  (void)cell;(void)T;(void)P;(void)b;(void)trace;(void)pi;return false;
#endif
}
