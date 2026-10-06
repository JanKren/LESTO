/*------------------------------------------------------------------------------

interfaceExchange.C

PURPOSE

Implementation of the gas/condensate exchange of rhoFixedFlowFoam; see
interfaceExchange.H for the model, the step algorithm and the state, and
phaseChangeKinetics.H for the local laws and the predictor.

Parallel hygiene: every reduction is called on all ranks, including ranks
without interface elements or sample cells; per-step counters are summed
locally first and reduced once; files are written by the master only.
Demand-driven mesh quantities (deltaCoeffs, Cf, magSf) are obtained on every
rank before a loop over elements or samples, never inside it: their first
construction may be collective (deltaCoeffs builds the lduAddressing and
polyMesh::globalData), and a rank without elements would skip it while the
others wait for it (review round 4: wallResistance halfCell deadlocked on
the ranks without WALL faces).

------------------------------------------------------------------------------*/

#include "interfaceExchange.H"
#include "gasRespeciation.H"
#include "gemsEquilibrium.H"
#include "matrixDefect.H"
#include "fvmDdt.H"
#include "compensatedSum.H"
#include "cellBitSet.H"
#include "emptyFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "processorFvPatch.H"
#include "fileOperation.H"
#include "decomposedBlockData.H"
#include "SpanStream.H"
#include "StringStream.H"
#include "ISstream.H"
#include "FixedList.H"
#include "Time.H"
#include "globalIndex.H"
#include "PstreamBuffers.H"
#include "boundBox.H"
#include "treeBoundBox.H"
#include "primitivePatch.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <tuple>

using namespace Foam;


/*------------------------------------------------------------------------------
Local helpers
------------------------------------------------------------------------------*/

namespace {

  /*----------------------------------------------------------------------------
  The directory of the per-rank element state below a time: a subdirectory
  of uniform/, not uniform/ itself (readState(); reconstructPar and
  decomposePar copy or link it with uniform/)
  ----------------------------------------------------------------------------*/
  const fileName stateLocal("uniform/phaseChange");

  /*----------------------------------------------------------------------------
  Compensated sum over all ranks: the (sum, compensation) pairs of the ranks
  are gathered and added in rank order, so every rank gets the same value,
  accurate to about eps whatever the decomposition (compensatedSum.H).
  ----------------------------------------------------------------------------*/
  scalar globalSum(const LESTO::compensatedSum& local) {

    if (!Pstream::parRun()) {
      return local.value();
    }

    List<FixedList<scalar, 2>> parts(Pstream::nProcs());
    parts[Pstream::myProcNo()][0] = local.sum();
    parts[Pstream::myProcNo()][1] = local.compensation();
    Pstream::allGatherList(parts);

    LESTO::compensatedSum total;
    for (const FixedList<scalar, 2>& part : parts) {
      total.add(part[0]);
      total.add(part[1]);
    }
    return total.value();
  }

  /*----------------------------------------------------------------------------
  64-bit FNV-1a hash of raw bytes (fingerprint of the interface elements).
  ----------------------------------------------------------------------------*/
  class fnv1a {

    std::uint64_t hash_;

  public:

    fnv1a() : hash_(14695981039346656037ULL) {}

    void add(const void* data, const std::size_t n) {
      const unsigned char* bytes = static_cast<const unsigned char*>(data);
      for (std::size_t i = 0; i < n; ++i) {
        hash_ ^= bytes[i];
        hash_ *= 1099511628211ULL;
      }
    }

    std::uint64_t value() const { return hash_; }
  };

  /*----------------------------------------------------------------------------
  True when the values are one repeated non-zero value.  OpenFOAM writes such
  a field or patch as 'uniform <value>', in text at writePrecision, also in
  a binary file (Field::writeEntry), so values read back like this may be
  rounded.  An exact 0 is written as "0" and reads back exactly.
  ----------------------------------------------------------------------------*/
  bool repeatedNonZero(const scalarField& values) {
    return values.uniform() && values.front() != 0;
  }

  /*----------------------------------------------------------------------------
  True when the values will be written as a ROUNDED 'uniform' entry: one
  repeated non-zero value whose text at the given precision does not read
  back to it.  A value that itself came from a text entry of at most that
  many digits (an inlet value 1e-3, a pressure 101325, a carrier written at
  that precision) reads back exactly.  Only the writer can decide this: the
  value re-read from such an entry always reproduces its own text.
  ----------------------------------------------------------------------------*/
  bool roundedUniform(const scalarField& values, const label precision) {

    if (!repeatedNonZero(values)) {
      return false;
    }

    OStringStream os;
    os.precision(precision);
    os << values.front();
    IStringStream is(os.str());
    scalar readBack = 0;
    is >> readBack;
    return readBack != values.front();
  }

  /*----------------------------------------------------------------------------
  True when the stored values of a patch are read back by a restart, so
  that a rounded 'uniform' entry of them reaches the next steps.  Not for a
  volume patch field whose value is recomputed when it is read:

    zeroGradient, fixedGradient, mixed   evaluated from the cells and the
                                         patch parameters (gradient;
                                         refValue, refGradient,
                                         valueFraction);
    symmetry, symmetryPlane, wedge       evaluated from the cells;
    empty                                no values;
    coupled patches (cyclic, processor)  evaluated from the neighbour cells.

  inletOutlet (a mixed condition) reads its value back and is checked.
  The parameters themselves (gradient, refValue, inletValue, ...) are
  written as 'uniform' text as well; they are NOT checked, so a restart
  verdict covers the inventory, the internal fields and the patch values
  only (review of the final M2 pass, carried over to M3).  A face field
  (phi) reads the values of all its patches back; empty ones have none.
  ----------------------------------------------------------------------------*/
  bool valueReadBack(const fvPatchScalarField& patchField) {
    static const wordHashSet recomputed ({
      "zeroGradient", "fixedGradient", "mixed", "symmetry", "symmetryPlane",
      "wedge", "empty"
    });
    return !patchField.coupled() && !recomputed.found(patchField.type());
  }

  bool valueReadBack(const fvsPatchScalarField& patchField) {
    return patchField.type() != "empty";
  }

  /*----------------------------------------------------------------------------
  The fewest significant decimal digits (1 to 17) whose text reads back to
  x exactly: 1 for 0.6, 10 for 0.6012345679, 16 or 17 for most doubles.  A
  value written as text with p significant digits needs at most p (a
  decimal of at most 15 digits survives the round trip through a double),
  so the most any value of a set needs is a lower bound of the digits it
  was written with.
  ----------------------------------------------------------------------------*/
  label significantDigits(const scalar x) {
    char text[40];
    for (int digits = 1; digits < 17; ++digits) {
      std::snprintf(text, sizeof(text), "%.*e", digits - 1, x);
      if (std::strtod(text, nullptr) == x) {
        return digits;
      }
    }
    return 17;
  }

  /*----------------------------------------------------------------------------
  Lexicographic order of two points (x, then y, then z): the tie rule of the
  nearest wall face of layer distance (readDistanceLayer()).
  ----------------------------------------------------------------------------*/
  bool lessPoint(const point& a, const point& b) {
    return std::lexicographical_compare(a.cbegin(), a.cend(),
                                        b.cbegin(), b.cend());
  }

  /*----------------------------------------------------------------------------
  The push of the nearest-face map of layer distance (nearestFaceSums() and
  the start-up count of readDistanceLayer()): the values of the compact
  slots of every rank are added onto the faces they stand for, on the rank
  that holds the face.  The sum of a face is

    (its own slot) + (the slot of rank 0) + (rank 1) + ... ,

  its own contribution first, then those of the other ranks in ascending
  rank order, every contribution exactly once: an order fixed by the
  decomposition alone, so identical runs give identical sums.  The values
  of the remote slots travel in PstreamBuffers (each rank sends to every
  rank that holds faces of its slots; a rank's message to proci lists the
  slots of constructMap[proci], which the receiver adds to its faces
  subMap[sender] in the same order).

  Why not a reverse mapDistributeBase::distribute with plusEqOp (review of
  M4, round 3; plan section 29, item 1).  Its nonBlocking branch adds the
  received slots in the order in which the receives complete, a sum that
  can differ in its last bit between identical runs.  Its scheduled
  branch, used in round 2 for a fixed order, takes the schedule of
  mapDistributeBase::schedule(), which lists DIRECTED pairs: when two ranks
  map elements to each other's faces, both (a, b) and (b, a) are in the
  schedule of each, and the branch exchanges both ways for each entry, so
  every remote contribution was added twice (v2412 and v2606:
  mapDistributeBase.C, schedule(); mapDistributeBaseTemplates.C, the
  scheduled branch).  mDep_ then exceeded the wall inventory, by +50 % on
  a channel split so that two ranks hold each other's second rows, and the
  start-up count of faces with elements of two or more other ranks counted
  every face with elements of one such rank.  Collective.
  ----------------------------------------------------------------------------*/
  template<class Type>
  List<Type> pushToFaces (
    const mapDistributeBase& map,
    const label              nFaces,
    const UList<Type>&       slots
  ) {

    if (map.constructHasFlip() || map.subHasFlip()) {
      FatalErrorInFunction << "The nearest-face map has flips; it is built "
        << "from global face indices and should have none." << nl
        << abort(FatalError);
    }

    const label comm = map.comm();
    const label myRank = UPstream::myProcNo(comm);
    const labelListList& constructMap = map.constructMap();
    const labelListList& subMap = map.subMap();

    /* its own slots first */
    List<Type> sums(nFaces, Zero);
    {
      const labelList& faces = subMap[myRank];
      const labelList& own = constructMap[myRank];
      forAll(faces, i) {
        sums[faces[i]] += slots[own[i]];
      }
    }

    if (!UPstream::parRun()) {
      return sums;
    }

    PstreamBuffers buffers(UPstream::commsTypes::nonBlocking,
                           UPstream::msgType(), comm);
    for (const int proci : UPstream::allProcs(comm)) {
      const labelList& remote = constructMap[proci];
      if (proci != myRank && remote.size()) {
        UOPstream os(proci, buffers);
        os << List<Type>(UIndirectList<Type>(slots, remote));
      }
    }
    buffers.finishedSends();

    /* the other ranks in ascending order */
    for (const int proci : UPstream::allProcs(comm)) {
      const labelList& faces = subMap[proci];
      if (proci != myRank && faces.size()) {
        UIPstream is(proci, buffers);
        List<Type> values;
        is >> values;
        if (values.size() != faces.size()) {
          FatalErrorInFunction << "Rank " << proci << " sent "
            << values.size() << " values for " << faces.size()
            << " faces of the nearest-face map." << nl << abort(FatalError);
        }
        forAll(faces, i) {
          sums[faces[i]] += values[i];
        }
      }
    }
    return sums;
  }
}


/*------------------------------------------------------------------------------
pairData
------------------------------------------------------------------------------*/

LESTO::interfaceExchange::pairData::pairData (
  const fvMesh& mesh,
  const word&   gasFieldName,
  const label   nElements
)
:
  gas(-1),
  condensed(-1),
  molarMass(0),
  formula(),
  accommodation(1), Ce(1), kineticScale(1),
  solverName(gasFieldName),
  finalSolver(false),
  table(),
  a(nElements, 0.0),
  g(nElements, 0.0),
  lower(nElements, -VGREAT),
  upper(nElements, 0.0),
  flow(nElements, 0.0),
  regime(nElements, label(UNKNOWN)),
  deposit(nElements, 0.0),
  reservoir(),
  inventory(false),
  Su (
    IOobject("phaseChangeSu_" + gasFieldName,
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false),
    mesh,
    dimensionedScalar(dimDensity/dimTime, Zero)
  ),
  Sp (
    IOobject("phaseChangeSp_" + gasFieldName,
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false),
    mesh,
    dimensionedScalar(dimDensity/dimTime, Zero)
  ),
  Srel (
    IOobject("phaseChangeSrel_" + gasFieldName,
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false),
    mesh,
    dimensionedScalar(dimDensity/dimTime, Zero)
  ),
  depositField(),
  wallDepositField(),
  exchangeField(),
  reservoirField(),
  releasedStep(0),
  inflowStep(0),
  nChanged(0),
  nChangedSignificant(0),
  transport(),
  defect(0),
  residualOF(0),
  explicitRemainder(0),
  nExplicitRemainder(0),
  explicitSource(0),
  transportIdentity(0),
  identityTerms(0),
  identityRelative(0),
  clampedStep(0),
  projectedStep(0),
  evaporatedAboveWarningStep(0),
  changedMass(0),
  nChangedOuter(0),
  nChangedSignificantOuter(0),
  changedMassOuter(0),
  nGuardPasses(0),
  nGuardSwitched(0),
  nUnconverged(0),
  nGuardPassesRun(0),
  nProjectedRun(0),
  nProjectedSignificantRun(0),
  nComplementarityRun(0),
  nComplementaritySignificantRun(0),
  nFinalChangesSignificantRun(0),
  nProjectedAboveRun(0),
  nComplementarityAboveRun(0)
{}


void LESTO::interfaceExchange::pairData::resizeElements (
  const label nElements,
  const label nWall,
  const label nSample
) {
  a.resize(nElements, 0.0);
  g.resize(nElements, 0.0);
  lower.resize(nElements, -VGREAT);
  upper.resize(nElements, 0.0);
  flow.resize(nElements, 0.0);
  regime.resize(nElements, label(UNKNOWN));
  deposit.resize(nWall, 0.0);
  reservoir.resize(nSample, 0.0);
}


/*------------------------------------------------------------------------------
Construction
------------------------------------------------------------------------------*/

LESTO::interfaceExchange::interfaceExchange (
  const fvMesh&                  mesh,
  const fluidThermo&             thermo,
  const volScalarField&          rho,
  const surfaceScalarField&      phi,
  PtrList<volScalarField>&       species,
  const wordList&                speciesNames,
  const List<word>&              state,
  const List<scalar>&            molarMass,
  const PtrList<volScalarField>& rhoD,
  const label                    nCorr,
  const scalar                   tolerance
)
:
  mesh_(mesh),
  thermo_(thermo),
  rho_(rho),
  phi_(phi),
  species_(species),
  speciesNames_(speciesNames),
  rhoD_(rhoD),
  nCorr_(nCorr),
  tolerance_(tolerance),
  dict_ (
    IOobject("thermochemistryProperties",
             mesh.time().constant(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE)
  ),
  model_(modelType::mock),
  pairOfSpecies_(speciesNames.size(), -1),
  pairs_(),
  accommodation_(1),
  Ce_(1),
  kineticScale_(1),
  gems_(),
  interfacePatches_(),
  layer_("faceCells"),
  layerDistance_(0),
  areaPerVolume_("cell"),
  areaMultiplier_(1),
  interfaceTemperature_("cell"),
  wallResistance_("none"),
  layerArea_(0),
  layerVolume_(0),
  nLayerCells_(0),
  nElementsTotal_(0),
  elementCell_(),
  elementPatch_(),
  elementFace_(),
  elementArea_(),
  elementWallArea_(),
  elementT_(),
  elementSample_(),
  exchangeCells_(),
  cellStart_(1, 0),
  nWall_(0),
  nWallCells_(0),
  nSampleElementsTotal_(0),
  interfaceFacePatch_(),
  interfaceFaceFace_(),
  elementFaceSlot_(),
  faceMap_(),
  depositionVelocity_(0),
  k_(0),
  Tdep_(0),
  samples_(),
  sampleCell_(),
  nGuard_(3),
  guardTolerance_(1e-12),
  moleFractionWarning_(0.05),
  writeFile_(true),
  reportMatrixResidual_(false),
  writeExchangeField_(false),
  profiles_(),
  profileWriter_(),
  profilesWritten_(false),
  profilesClock_(0),
  profilesAtWrite_(false),
  profilesTimes_(),
  profilesPending_(),
  carrierConcentration_(),
  stepStart_(0),
  stepEnd_(0),
  stepDeltaT_(mesh.time().deltaTValue()),
  ledger_(),
  fingerprint_(),
  elementsFingerprint_(),
  samplesFingerprint_(),
  regimesIO_(),
  depositsIO_(),
  reservoirsIO_(),
  equilibriumIO_()
{
  /*--------------------------------------------------------------------------
  model mock, or no dictionary: nothing else is read or printed, so the
  solver output is unchanged.
  --------------------------------------------------------------------------*/
  const word model = dict_.getOrDefault<word>("model", "mock");

  if (model == "mock") {
    return;
  }

  if (model != "temperature" && model != "HKS") {
    FatalIOErrorInFunction(dict_)
      << "Unknown model " << model << ". Choose mock, temperature or HKS."
      << exit(FatalIOError);
  }

  model_ = model == "HKS" ? modelType::HKS : modelType::temperature;

  Info<< "Phase change: model " << model
    << " (constant/thermochemistryProperties)" << nl;

  /*--------------------------------------------------------------------------
  The wall elements (readInterface), the pairs with their p_eq tables, the
  coefficients of the model, and the samples, whose inventory samples add
  one sample element per cell after the wall elements; the element lists of
  the pairs are sized once all elements are known.
  --------------------------------------------------------------------------*/
  readInterface();
  if (model_ == modelType::HKS) {
    readHKSModel();
  }
  readPairs(state, molarMass);
  if (model_ == modelType::temperature) {
    readTemperatureModel();
  }
  readSamples();

  /*--------------------------------------------------------------------------
  Inputs of the other model are not read.  They may stay in the file (to
  switch the model back), but they are named, so that a model keyword left
  at the wrong value does not go unnoticed (review of M5).
  --------------------------------------------------------------------------*/
  {
    OStringStream ignored;
    const word other =
      model_ == modelType::HKS ? "temperatureCoeffs" : "HKSCoeffs";
    if (dict_.found(other, keyType::LITERAL)) {
      ignored << ' ' << other;
    }
    if (model_ == modelType::temperature) {
      if (dict_.found("GEMSCoeffs", keyType::LITERAL)) {
        ignored << " GEMSCoeffs";
      }
      for (const entry& e : dict_.subDict("pairs")) {
        for (const word key : {"vapourPressure", "gems"}) {
          if (e.isDict() && e.dict().found(key, keyType::LITERAL)) {
            ignored << " pairs/" << e.keyword() << '/' << key;
          }
        }
      }
    }
    if (!ignored.str().empty()) {
      WarningInFunction
        << "Model " << model << " does not read" << ignored.str().c_str()
        << " (inputs of the other model) in constant/"
        << "thermochemistryProperties." << nl << endl;
    }

    /* model HKS with equilibrium table: the inputs of the GEMS3K backend
       (M9) are not read either */
    if (model_ == modelType::HKS && !gems_) {
      OStringStream gemsInputs;
      if (dict_.found("GEMSCoeffs", keyType::LITERAL)
        && !(dict_.isDict("respeciation") && dict_.subDict("respeciation").getOrDefault<word>("engine","none")=="GEMS")) {
        gemsInputs << " GEMSCoeffs";
      }
      for (const entry& e : dict_.subDict("pairs")) {
        if (e.isDict() && e.dict().found("gems", keyType::LITERAL)) {
          gemsInputs << " pairs/" << e.keyword() << "/gems";
        }
      }
      if (!gemsInputs.str().empty()) {
        WarningInFunction
          << "equilibrium table does not read" << gemsInputs.str().c_str()
          << " (inputs of equilibrium GEMS) in constant/"
          << "thermochemistryProperties." << nl << endl;
      }
    }
  }

  forAll(pairs_, pairi) {
    pairs_[pairi].resizeElements (
      elementCell_.size(), nWall_, elementCell_.size() - nWall_
    );
  }
  for (const sampleData& s : samples_) {
    if (s.inventory) {
      pairs_[s.pair].inventory = true;
    }
  }

  /*--------------------------------------------------------------------------
  Balance controls
  --------------------------------------------------------------------------*/
  const dictionary& balanceDict = dict_.optionalSubDict("balance");
  nGuard_ = balanceDict.getOrDefault<label>("nGuard", 3);
  guardTolerance_ = balanceDict.getOrDefault<scalar>("guardTolerance", 1e-12);
  moleFractionWarning_ =
    balanceDict.getOrDefault<scalar>("moleFractionWarning", 0.05);
  writeFile_ = balanceDict.getOrDefault<bool>("writeFile", true);
  reportMatrixResidual_ =
    balanceDict.getOrDefault<bool>("reportMatrixResidual", false);
  writeExchangeField_ =
    balanceDict.getOrDefault<bool>("writeExchangeField", false);

  if (nGuard_ < 0 || !(guardTolerance_ >= 0) || !(moleFractionWarning_ > 0)) {
    FatalIOErrorInFunction(balanceDict)
      << "Require nGuard >= 0, guardTolerance >= 0 and "
      << "moleFractionWarning > 0." << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  equilibrium GEMS (M9): the engine of every rank, checked against the
  system; frozen builds the GEMS table of every pair from its table, local
  sizes the warm states of the elements of the rank (gemsEquilibrium.H).
  Before setCoefficients(), which takes p_eq from it.  Collective.
  --------------------------------------------------------------------------*/
  if (gems_) {
    List<const vapourPressureTable*> tables(pairs_.size());
    forAll(pairs_, pairi) {
      tables[pairi] = pairs_[pairi].table.get();
    }
    gems_->start(elementCell_.size(), tables);
  }

  /*--------------------------------------------------------------------------
  Start-up checks of the paired species, element coefficients (T is frozen,
  so they are evaluated once), and the transport flux of the gas equation
  for the ledger.  The processor patches of the paired gas are evaluated
  from the neighbour cells: the file may hold them as 'uniform <value>' at
  writePrecision (a single-face patch), and they enter the first matrix of
  the run (header: 'Processor patches of Y_<gas>').  Serial runs have none.
  --------------------------------------------------------------------------*/
  forAll(pairs_, pairi) {
    pairData& pair = pairs_[pairi];
    checkPairedSpecies(pair);
    setCoefficients(pair);
    mesh_.setFluxRequired(pair.solverName);
    species_[pair.gas].boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  }

  /*--------------------------------------------------------------------------
  Carrier: its molar concentration p/(R T) for the mole-fraction monitor,
  and the consistency of the frozen T with the perfect-gas law (plan
  section 2, E11).  thermo.W() is in kg/kmol.  The concentration comes
  from p and T, not from rho/W: in paper mode rho = 1 kg/m3 is no density
  of helium, and rho/W would count about 250 mol/m3 of carrier instead of
  10 to 40 (review of M4, round 3; plan section 29, item 4).  For a
  carrier that obeys the perfect-gas law both agree; for the frozen
  carrier of run/007 within its deviation printed here (0.3 %).
  --------------------------------------------------------------------------*/
  {
    const tmp<volScalarField> tW(thermo_.W());
    const scalarField W(tW().primitiveField()/1000.0);

    const scalarField& T = thermo_.T().primitiveField();
    const scalarField& p = thermo_.p().primitiveField();
    const scalarField& rhoI = rho_.primitiveField();

    carrierConcentration_.resize(T.size());
    scalar maxDeviation = 0;
    forAll(T, celli) {
      carrierConcentration_[celli] = p[celli]/(universalGasConstant*T[celli]);
      const scalar Teos = p[celli]*W[celli]
                        /(universalGasConstant*rhoI[celli]);
      maxDeviation = max(maxDeviation, mag(T[celli] - Teos)/T[celli]);
    }
    reduce(maxDeviation, maxOp<scalar>());

    Info<< "Carrier: max |T - p W/(R rho)|/T = " << maxDeviation
      << " (frozen T vs perfect gas)" << nl;

    /*------------------------------------------------------------------------
    equilibrium GEMS, mode local, evaluates GEMS3K at the composition of
    every cell from the carrier's p/(R T): it needs a carrier that is a gas
    at its p and T.  |T - p W/(R rho)|/T is |rho - p W/(R T)|/rho; above 1 %
    (paper mode, rho = 1) the run stops unless allowNonPhysicalCarrier
    (plan section 4).
    ------------------------------------------------------------------------*/
    if (gems_ && gems_->local() && maxDeviation > 0.01) {
      if (!gems_->allowNonPhysicalCarrier()) {
        FatalIOErrorInFunction(gems_->dict())
          << "GEMSCoeffs mode local evaluates GEMS3K at the composition of "
          << "every cell, from the carrier's p/(R T), but the frozen carrier "
          << "is not a gas at its p and T: max |rho - p W/(R T)|/rho = "
          << maxDeviation << " > 0.01 (paper mode has rho = 1)." << nl
          << "Use mode frozen or equilibrium table, or set "
          << "allowNonPhysicalCarrier yes to accept it." << exit(FatalIOError);
      }
      WarningInFunction
        << "GEMSCoeffs allowNonPhysicalCarrier yes: mode local runs on a "
        << "carrier that is not a gas at its p and T (max |rho - p W/(R T)|"
        << "/rho = " << maxDeviation << " > 0.01); the composition of every "
        << "cell is taken from p/(R T)." << nl << endl;
    }
  }

  /*--------------------------------------------------------------------------
  equilibrium GEMS with writeEquilibriumField: the output fields pEq_<gas>
  [Pa], gemsSource_<gas> and, with speciation lagged, chi_<gas> (refreshed
  in balance(), never read; the sixth comment block of interfaceExchange.H).
  --------------------------------------------------------------------------*/
  if (gems_ && gems_->writeEquilibriumField()) {
    auto outputField = [this](const word& name, const dimensionSet& dims) {
      return autoPtr<volScalarField> (
        new volScalarField (
          IOobject(name,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar(dims, Zero),
          calculatedFvPatchScalarField::typeName
        )
      );
    };
    forAll(pairs_, pairi) {
      pairData& pair = pairs_[pairi];
      const word& gasName = speciesNames_[pair.gas];
      pair.equilibriumField = outputField("pEq_" + gasName, dimPressure);
      pair.sourceField = outputField("gemsSource_" + gasName, dimless);
      if (gems_->speciation()) {
        pair.chiField = outputField("chi_" + gasName, dimless);
      }
      updateEquilibriumFields(pair);
    }
  }

  readState();
  recordFormulas();

  if (writeFile_) {
    ledger_->openFiles();
  }

  /*--------------------------------------------------------------------------
  Axial profiles (optional; axialProfiles.H).  A restart skips exactly the
  profile times recorded as written by the runs before it (uniform/
  phaseChangeLayout, profiles { written (...); }) and those at or before
  its start; a fresh start (or a layout without the record) skips the
  entries whose nearest step end is not after the start.  The cells of
  every sample are excluded from the observable wallPlusGas; the
  significance of an onset defaults to guardTolerance.  profileWriter hooks
  them into every write of Time (interfaceExchange.H).
  --------------------------------------------------------------------------*/
  if (dict_.found("profiles")) {

    autoPtr<scalarList> written;
    const dictionary* record = ledger_->restarted()
      ? ledger_->layout().findDict("profiles", keyType::LITERAL)
      : nullptr;
    if (record) {
      written.reset(new scalarList());
      for (const string& t : record->get<List<string>>("written")) {
        written().append(readScalar(t));
      }
    }

    /* the wall elements (layer distance: binned by their cells), and the
       faces of the wall temperature per bin: the faces of the elements
       (faceCells), or every face of the interface patches (distance) */
    const labelList wallCells(SubList<label>(elementCell_, nWall_));
    profiles_.reset (
      new axialProfiles (
        mesh_,
        dict_.subDict("profiles"),
        elementPatch_,
        elementFace_,
        wallCells,
        distanceLayer() ? interfaceFacePatch_ : elementPatch_,
        distanceLayer() ? interfaceFaceFace_ : elementFace_,
        thermo_.T(),
        sampleCell_,
        guardTolerance_,
        ledger_->time(),
        mesh_.time().deltaTValue(),
        written.get(),
        ledger_->restarted()
      )
    );
    recordProfileTimes();
    profileWriter_.reset (
      new profileWriter (
        IOobject("phaseChangeProfiles",
                 mesh_.time().timeName(),
                 mesh_.time(),
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE),
        *this
      )
    );
  }
}


/*------------------------------------------------------------------------------
Pairs: one sub-dictionary per paired gas species, naming its condensate
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readPairs (
  const List<word>&   state,
  const List<scalar>& molarMass
) {

  const IOdictionary properties(IOobject(
    "speciesTransportProperties", mesh_.time().constant(), mesh_,
    IOobject::MUST_READ, IOobject::NO_WRITE, false));
  const dictionary& pairsDict = dict_.subDict("pairs");
  std::map<word, scalar> elements;
  label formulaPairs = 0;

  if (pairsDict.empty()) {
    FatalIOErrorInFunction(pairsDict) << "No pair is defined."
      << exit(FatalIOError);
  }

  labelList condensateOwner(speciesNames_.size(), -1);
  pairs_.resize(pairsDict.size());
  label pairi = 0;

  for (const entry& e : pairsDict) {

    const word& gasName = e.keyword();

    if (!e.isDict()) {
      FatalIOErrorInFunction(pairsDict) << "Pair " << gasName
        << " must be a dictionary." << exit(FatalIOError);
    }

    const dictionary& pd = e.dict();
    const word condensedName = pd.get<word>("condensed");
    const label gi = speciesNames_.find(gasName);
    const label ci = speciesNames_.find(condensedName);

    if (gi < 0 || state[gi] != "gas") {
      FatalIOErrorInFunction(pd) << "Pair " << gasName
        << ": not a gas species of speciesTransportProperties."
        << exit(FatalIOError);
    }
    if (ci < 0 || state[ci] != "solid") {
      FatalIOErrorInFunction(pd) << "Pair " << gasName << ": condensed "
        << condensedName << " is not a solid species of "
        << "speciesTransportProperties." << exit(FatalIOError);
    }
    if (condensateOwner[ci] >= 0) {
      FatalIOErrorInFunction(pd) << "Condensate " << condensedName
        << " belongs to two pairs." << exit(FatalIOError);
    }
    if (!(molarMass[gi] > 0)) {
      FatalIOErrorInFunction(pd) << "Pair " << gasName << " -> "
        << condensedName << ": the gas species " << gasName
        << " needs molarMass [kg/mol] in speciesTransportProperties."
        << exit(FatalIOError);
    }
    if (molarMass[ci] > 0 && mag(molarMass[ci]/molarMass[gi] - 1) > 1e-12) {
      FatalIOErrorInFunction(pd) << "Pair " << gasName << " -> "
        << condensedName << ": the condensate molarMass " << molarMass[ci]
        << " differs from the gas molarMass " << molarMass[gi]
        << ".  The mass balance of a pair assumes the same formula."
        << exit(FatalIOError);
    }

    const word gasFieldName = species_[gi].name();
    pairs_.set(pairi, new pairData(mesh_, gasFieldName, elementCell_.size()));
    pairData& pair = pairs_[pairi];
    pair.gas = gi;
    pair.condensed = ci;
    pair.molarMass = molarMass[gi];
    pair.formula = readFormula(properties.subDict(gasName), pair.molarMass);
    const auto condensedFormula = readFormula(
      properties.subDict(condensedName),
      molarMass[ci] > 0 ? molarMass[ci] : pair.molarMass
    );
    if (pair.formula != condensedFormula) {
      FatalIOErrorInFunction(pd) << "Pair " << gasName
        << " requires identical gas and condensate formulas."
        << exit(FatalIOError);
    }
    if (!pair.formula.empty()) ++formulaPairs;
    for (const auto& atom : pair.formula) elements.emplace(atom);
    pair.accommodation = accommodation_;
    pair.Ce = Ce_;
    pair.kineticScale = kineticScale_;
    if (pd.found("HKS")) {
      pair.pairHKS = true;
      if (model_ != modelType::HKS) {
        FatalIOErrorInFunction(pd) << "Per-pair HKS needs model HKS."
          << exit(FatalIOError);
      }
      const dictionary& hks = pd.subDict("HKS");
      pair.accommodation = hks.getOrDefault<scalar>("accommodation", accommodation_);
      pair.Ce = hks.getOrDefault<scalar>("Ce", Ce_);
      pair.kineticScale = hks.getOrDefault<scalar>("kineticScale", kineticScale_);
      if (!std::isfinite(pair.accommodation) || pair.accommodation <= 0
        || pair.accommodation > 1 || !std::isfinite(pair.Ce) || pair.Ce < 0
        || !std::isfinite(pair.kineticScale) || pair.kineticScale < 0) {
        FatalIOErrorInFunction(hks) << "Require finite 0 < accommodation <= 1 "
          << "and Ce, kineticScale >= 0." << exit(FatalIOError);
      }
    }

    pairOfSpecies_[gi] = pairi;
    condensateOwner[ci] = pairi;

    if (model_ == modelType::HKS) {
      readVapourPressure(pair, pd);
    }

    /* equilibrium GEMS: the species of the pair in the GEMS3K system */
    if (gems_) {
      if (!pd.isDict("gems")) {
        FatalIOErrorInFunction(pd) << "Pair " << gasName << ": equilibrium "
          << "GEMS needs gems { gas \"<GEMS3K species>\"; condensates "
          << "(...); elements { <element> <coefficient>; ... } } (optional: "
          << "excess { <element> <relative excess>; })." << exit(FatalIOError);
      }
      gems_->addPair(gasName, pd.subDict("gems"), pair.molarMass);
    }
    ++pairi;
  }

  if (formulaPairs && formulaPairs != pairs_.size()) {
    FatalIOErrorInFunction(pairsDict) << "Element accounting needs formulas "
      << "for every pair." << exit(FatalIOError);
  }
  elementNames_.resize(elements.size());
  label ai = 0;
  for (const auto& atom : elements) elementNames_[ai++] = atom.first;

  /*--------------------------------------------------------------------------
  A solid species without a pair would never change: there is no model for
  it outside the mock path.
  --------------------------------------------------------------------------*/
  forAll(state, si) {
    bool channelSolid = false;
    if(dict_.isDict("wallChannels")) for(const entry& channel:dict_.subDict("wallChannels"))
      if(channel.isDict() && channel.dict().getOrDefault<word>("condensed",word::null)==speciesNames_[si])channelSolid=true;
    if (state[si] == "solid" && condensateOwner[si] < 0 && !channelSolid) {
      FatalIOErrorInFunction(pairsDict) << "Solid species "
        << speciesNames_[si] << " is not the condensate of any pair; "
        << "outside model mock every solid species must belong to a pair."
        << exit(FatalIOError);
    }
  }
}


/*------------------------------------------------------------------------------
Interface elements of the interaction layer
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readInterface() {

  const dictionary& d = dict_.subDict("interface");

  /*--------------------------------------------------------------------------
  The interaction layer: faceCells (the cells of the wall faces, one element
  per face) or distance (the cells whose centre lies within 'distance' of
  the wall, one element per cell; fifth comment block of
  interfaceExchange.H).  The area per volume of a distance layer is the
  layer mean by definition: its default is 'layer', and 'cell' is refused.
  --------------------------------------------------------------------------*/
  layer_ = d.getOrDefault<word>("layer", "faceCells");
  if (layer_ != "faceCells" && layer_ != "distance") {
    FatalIOErrorInFunction(d) << "Unknown layer " << layer_
      << ". Choose faceCells or distance." << exit(FatalIOError);
  }
  if (distanceLayer()) {
    layerDistance_ = d.get<scalar>("distance");
    if (!(layerDistance_ > 0) || !std::isfinite(layerDistance_)) {
      FatalIOErrorInFunction(d) << "layer distance needs a positive, finite "
        << "distance [m], not " << layerDistance_ << "." << exit(FatalIOError);
    }
  }

  areaPerVolume_ = d.getOrDefault<word> (
    "areaPerVolume", distanceLayer() ? "layer" : "cell"
  );
  areaMultiplier_ = d.getOrDefault<scalar>("areaMultiplier", 1);
  interfaceTemperature_ = d.getOrDefault<word>("interfaceTemperature", "cell");
  wallResistance_ = d.getOrDefault<word>("wallResistance", "none");

  if (areaPerVolume_ != "cell" && areaPerVolume_ != "layer") {
    FatalIOErrorInFunction(d) << "areaPerVolume must be cell or layer."
      << exit(FatalIOError);
  }
  if (distanceLayer() && areaPerVolume_ != "layer") {
    FatalIOErrorInFunction(d) << "layer distance needs areaPerVolume layer "
      << "(the default), not " << areaPerVolume_ << "." << nl
      << "Its element areas are A_e = V_c (A_I/V_I) areaMultiplier: every "
      << "cell of the layer has the mean area per volume A_I/V_I of the "
      << "layer, and a cell of the second layer has no wall face of its own "
      << "whose area could be used (areaPerVolume cell)." << exit(FatalIOError);
  }
  if (!(areaMultiplier_ > 0)) {
    FatalIOErrorInFunction(d) << "areaMultiplier must be positive."
      << exit(FatalIOError);
  }
  if (interfaceTemperature_ != "cell" && interfaceTemperature_ != "wall") {
    FatalIOErrorInFunction(d) << "interfaceTemperature must be cell or wall."
      << exit(FatalIOError);
  }
  if (wallResistance_ != "none" && wallResistance_ != "halfCell") {
    FatalIOErrorInFunction(d) << "wallResistance must be none or halfCell."
      << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  wallResistance halfCell is the series law of the kinetics and the
  diffusion through the half cell, which is an OpenFOAM mixed (Robin)
  condition per unit PHYSICAL wall area |S_f|.  The element law is

    g_e = A_e rho v',   v' = 1/(1/v + R),   v = v_d f(T_e),
    R   = rho d_e/(rhoD)_e   [s/m], the resistance of the half cell.

  With a scaled element area A_e = m |S_f| (areaMultiplier m, or the
  rescaling of areaPerVolume layer) this is

    g_e = m |S_f| rho/(1/v + R) = |S_f| rho/(1/(m v) + R/m):

  the scale multiplies the kinetic velocity AND divides the half-cell
  resistance by m, i.e. it makes the diffusion through the half cell m
  times more conductive.  The half cell belongs to the mesh: the diffusive
  conductance between the cell centre and a face of area |S_f| is
  |S_f| rhoD/d_e, whatever the kinetic area.  A Robin condition with the
  kinetic area m |S_f| on that face gives

    g   = |S_f| rho/(1/(m v) + R),

  which differs from the scaled element law unless R << 1/(m v) (the
  kinetic limit, where both are m |S_f| rho v).  Scaling A_e by m is
  therefore not Robin-equivalent, and the combination is refused until a
  scaling is decided (M7; review, final M2 pass).

  layer distance has no half cell at all: a cell of its second layer has no
  wall face, and the area V_c A_I/V_I of an element is no mesh face.
  --------------------------------------------------------------------------*/
  if (distanceLayer() && wallResistance_ == "halfCell") {
    FatalIOErrorInFunction(d)
      << "wallResistance halfCell is not available with layer distance." << nl
      << "The half cell is the diffusion between a wall face and the centre "
      << "of its cell (a Robin condition per wall face).  The elements of "
      << "layer distance are cells, most of them without a wall face (the "
      << "second layer), and their areas V_c A_I/V_I are no mesh faces." << nl
      << "Use wallResistance none, or layer faceCells for halfCell."
      << exit(FatalIOError);
  }
  if (wallResistance_ == "halfCell"
    && (areaMultiplier_ != 1 || areaPerVolume_ != "cell")) {
    FatalIOErrorInFunction(d)
      << "wallResistance halfCell needs areaMultiplier 1 and areaPerVolume "
      << "cell, not areaMultiplier " << areaMultiplier_ << " and "
      << "areaPerVolume " << areaPerVolume_ << "." << nl
      << "The element law is g_e = A_e rho v' with v' = 1/(1/v + R), "
      << "R = rho d/rhoD the resistance of the half cell.  With a scaled "
      << "element area A_e = m |S_f| it is |S_f| rho/(1/(m v) + R/m): the "
      << "scale divides the half-cell resistance by m as well, whereas the "
      << "half cell of the mesh face is the same for any kinetic area (the "
      << "Robin condition with the kinetic area m |S_f| is "
      << "|S_f| rho/(1/(m v) + R)).  The two agree only in the kinetic "
      << "limit R << 1/(m v)." << nl
      << "Use wallResistance none with a scaled area." << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  Patches: by name; they must be ordinary (non-coupled, non-empty) patches.
  --------------------------------------------------------------------------*/
  const wordList names(d.get<wordList>("patches"));
  if (names.empty()) {
    FatalIOErrorInFunction(d) << "No interface patch listed."
      << exit(FatalIOError);
  }

  interfacePatches_.resize(names.size(), -1);
  forAll(names, i) {
    const label patchi = mesh_.boundaryMesh().findPatchID(names[i]);
    if (patchi < 0) {
      FatalIOErrorInFunction(d) << "Interface patch " << names[i]
        << " not found." << exit(FatalIOError);
    }
    if (mesh_.boundary()[patchi].coupled()
      || isA<emptyFvPatch>(mesh_.boundary()[patchi])) {
      FatalIOErrorInFunction(d) << "Interface patch " << names[i]
        << " is coupled or empty." << exit(FatalIOError);
    }
    if (interfacePatches_.found(patchi)) {
      FatalIOErrorInFunction(d) << "Interface patch " << names[i]
        << " is listed twice." << exit(FatalIOError);
    }
    interfacePatches_[i] = patchi;
  }

  /*--------------------------------------------------------------------------
  layer faceCells: one element per face of the interface patches, sorted
  by cell so that the elements of one cell are contiguous.  layer distance:
  one element per cell of the layer (readDistanceLayer()), whose areas
  follow from A_I/V_I below.
  --------------------------------------------------------------------------*/
  scalarField physicalArea;
  scalarField interfaceFaceArea;
  string distanceReport;

  if (distanceLayer()) {
    distanceReport = readDistanceLayer(d, interfaceFaceArea);
  } else {

    DynamicList<label>  cells, patches, faces;
    DynamicList<scalar> areas;
    const surfaceScalarField::Boundary& magSfb =
      mesh_.magSf().boundaryField();

    for (const label patchi : interfacePatches_) {
      const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
      const scalarField& magSf = magSfb[patchi];
      forAll(faceCells, facei) {
        cells.append(faceCells[facei]);
        patches.append(patchi);
        faces.append(facei);
        areas.append(magSf[facei]);
      }
    }

    const labelList order(sortedOrder(cells));
    const label nElements = order.size();

    elementCell_.resize(nElements);
    elementPatch_.resize(nElements);
    elementFace_.resize(nElements);
    physicalArea.resize(nElements);

    forAll(order, e) {
      elementCell_[e]  = cells[order[e]];
      elementPatch_[e] = patches[order[e]];
      elementFace_[e]  = faces[order[e]];
      physicalArea[e]  = areas[order[e]];
    }

    DynamicList<label> exchange, start;
    start.append(0);
    forAll(elementCell_, e) {
      if (e == 0 || elementCell_[e] != elementCell_[e-1]) {
        if (e > 0) {
          start.append(e);
        }
        exchange.append(elementCell_[e]);
      }
    }
    if (nElements > 0) {
      start.append(nElements);
    }
    exchangeCells_.transfer(exchange);
    cellStart_.transfer(start);
    nWall_ = nElements;
    nWallCells_ = exchangeCells_.size();
  }

  const label nElements = elementCell_.size();
  elementArea_.resize(nElements);

  /*--------------------------------------------------------------------------
  Layer: A_I (physical wall area: the faces of the interface patches), V_I
  (the volume of the layer cells) and A_I/V_I.

  With areaPerVolume layer (and so with layer distance) every element area
  is scaled by A_I/V_I, so A_I and V_I are compensated sums over all ranks,
  independent of the decomposition to about eps: plain sums over the 30528
  layer cells of the pipe differ between decompositions by up to some 1e-14
  (relative), which moved the wall inventory of a restart onto another
  decomposition by as much (M4: 2.0e-14 from 16 ranks to serial, beyond the
  round-off bound 1e-14 of an exact restart: FatalError; 7.3e-15 from
  serial to 4 ranks; with compensated sums 0).  With areaPerVolume cell the
  element areas are the face areas and A_I, V_I enter only v_d = A V_I/A_I;
  there the plain sums are kept, so that the runs of M2, M3 and M5 (e.g. the
  out_excerpt of run/012, whose round-off digits follow the last bit of
  v_d) are unchanged.
  --------------------------------------------------------------------------*/
  const scalarField& V = mesh_.V();
  if (areaPerVolume_ == "layer") {
    compensatedSum area;
    for (const scalar a : distanceLayer() ? interfaceFaceArea : physicalArea) {
      area.add(a);
    }
    compensatedSum volume;
    for (const label celli : exchangeCells_) {
      volume.add(V[celli]);
    }
    layerArea_ = globalSum(area);
    layerVolume_ = globalSum(volume);
  } else {
    layerArea_ = gSum(physicalArea);
    layerVolume_ = 0;
    for (const label celli : exchangeCells_) {
      layerVolume_ += V[celli];
    }
    reduce(layerVolume_, sumOp<scalar>());
  }
  nLayerCells_ = returnReduce(exchangeCells_.size(), sumOp<label>());
  nElementsTotal_ = returnReduce(nElements, sumOp<label>());

  if (nElementsTotal_ == 0) {
    FatalIOErrorInFunction(d) << "The interface patches have no face"
      << (distanceLayer() ? ", or no cell centre lies within the distance" : "")
      << "." << exit(FatalIOError);
  }

  const scalar layerAreaPerVolume = layerArea_/layerVolume_;

  /*--------------------------------------------------------------------------
  layer distance: every element represents the wall area V_c A_I/V_I (they
  add up to A_I); its area per volume is the layer mean.
  --------------------------------------------------------------------------*/
  if (distanceLayer()) {
    physicalArea.resize(nElements);
    forAll(elementCell_, e) {
      physicalArea[e] = V[elementCell_[e]]*layerAreaPerVolume;
    }
  }

  /*--------------------------------------------------------------------------
  Element areas: |S_f|*areaMultiplier, and for areaPerVolume layer scaled by
  (A_I/V_I)/(A_w,c/V_c), so that every cell of the layer has the mean area
  per volume A_I/V_I (layer distance: V_c (A_I/V_I) areaMultiplier, already
  the layer mean).
  --------------------------------------------------------------------------*/
  elementWallArea_ = physicalArea;
  elementArea_ = areaMultiplier_*physicalArea;

  if (areaPerVolume_ == "layer" && !distanceLayer()) {
    forAll(exchangeCells_, k) {
      const label celli = exchangeCells_[k];
      scalar cellArea = 0;
      for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
        cellArea += physicalArea[e];
      }
      const scalar scale = layerAreaPerVolume*V[celli]/cellArea;
      for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
        elementArea_[e] *= scale;
      }
    }
  }

  /*--------------------------------------------------------------------------
  Element temperature (frozen): the cell or the wall face (layer distance:
  the nearest wall face, possibly on another rank; collective).
  --------------------------------------------------------------------------*/
  const volScalarField& T = thermo_.T();
  elementT_.resize(nElements);

  if (distanceLayer() && interfaceTemperature_ == "wall") {
    scalarField faceT(interfaceFacePatch_.size());
    forAll(faceT, i) {
      faceT[i] =
        T.boundaryField()[interfaceFacePatch_[i]][interfaceFaceFace_[i]];
    }
    elementT_ = nearestFaceValues(faceT);
  } else {
    forAll(elementT_, e) {
      elementT_[e] = interfaceTemperature_ == "wall"
        ? T.boundaryField()[elementPatch_[e]][elementFace_[e]]
        : T[elementCell_[e]];
    }
  }

  Info<< "Interface (layer " << layer_ << "): patches (";
  forAll(names, i) {
    Info<< (i ? " " : "") << names[i];
  }
  Info<< "), ";
  if (distanceLayer()) {
    Info<< "distance " << layerDistance_ << " m: ";
  }
  Info<< nElementsTotal_ << " elements in " << nLayerCells_
    << " cells; A_I = "
    << layerArea_ << " m2, V_I = " << layerVolume_ << " m3, A_I/V_I = "
    << layerAreaPerVolume << " 1/m" << nl
    << "  areaPerVolume " << areaPerVolume_ << ", areaMultiplier "
    << areaMultiplier_ << ", interfaceTemperature " << interfaceTemperature_
    << ", wallResistance " << wallResistance_ << nl;
  if (distanceLayer()) {
    Info<< distanceReport.c_str();
  }
}


/*------------------------------------------------------------------------------
layer distance: the cells within the distance of the interface patches,
their nearest wall faces and the face map (fifth comment block of
interfaceExchange.H)
------------------------------------------------------------------------------*/

Foam::string LESTO::interfaceExchange::readDistanceLayer (
  const dictionary& d,
  scalarField&      interfaceFaceArea
) {

  const fvBoundaryMesh& boundary = mesh_.boundary();
  const polyBoundaryMesh& patches = mesh_.boundaryMesh();
  const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
  const pointField& meshPoints = mesh_.points();
  const vectorField& C = mesh_.C().primitiveField();

  /*--------------------------------------------------------------------------
  The interface faces of this rank, in the order of interfacePatches_: their
  local index is the slot of the face in the face map, and globalIndex
  numbers them over all ranks (collective).
  --------------------------------------------------------------------------*/
  label nFaces = 0;
  for (const label patchi : interfacePatches_) {
    nFaces += boundary[patchi].size();
  }
  interfaceFacePatch_.resize(nFaces);
  interfaceFaceFace_.resize(nFaces);
  interfaceFaceArea.resize(nFaces);
  List<pointField> facePoints(nFaces);
  List<boundBox> faceBoxes(nFaces);
  {
    label i = 0;
    for (const label patchi : interfacePatches_) {
      forAll(boundary[patchi], facei) {
        interfaceFacePatch_[i] = patchi;
        interfaceFaceFace_[i] = facei;
        interfaceFaceArea[i] = magSfb[patchi][facei];
        facePoints[i] = patches[patchi][facei].points(meshPoints);
        faceBoxes[i] = boundBox(facePoints[i], false);
        ++i;
      }
    }
  }
  const globalIndex faceNumbering(nFaces);

  /*--------------------------------------------------------------------------
  The wall distance of a cell is the exact distance of its centre from the
  nearest face of the interface patches, over the faces of ALL patches and
  ALL ranks (face::nearestPoint: the triangles of the face about its
  centre, exact for a planar face), and its nearest wall face f(e) is that
  face.  This is the definition of plan section 2; it replaces the
  patchDataWave<wallPointData<label>> wave of the first M4 round, whose
  near-wall correction (cellDistFuncs) is wrong in three ways (review of
  M4, round 1; plan section 26):

    - with the default combined wall patch it stores the label
      nWalls + minFacei of the combined patch where patchDataWave expects
      a mesh face, so the nearest faces of the cells next to the wall are
      unrelated faces (OpenFOAM v2412 and v2606);
    - per patch (useCombinedWallPatch false) every interface patch
      overwrites the distance of a cell with faces on several patches (the
      last patch wins), while the label of the first patch is kept, so a
      corner cell of two interface patches is dropped from the layer, or
      mapped to the farther face, depending on the order of the patches;
    - either way it sees only the wall faces of the cell's own rank: a
      cell that touches the wall by a point or an edge whose wall faces lie
      on another rank keeps the wave distance (to the nearest face centre),
      which is larger, so the layer depended on the decomposition (a tet
      mesh lost a cell on 4 ranks, and its serial state restarted on the 4
      ranks was refused).

  The search is local to every rank once every rank holds the interface
  faces that can lie within the search radius (2 x distance, which also
  gives the margin of the rule below) of one of its cell centres: the
  bounding boxes of the cell centres of the ranks, grown by the radius, are
  gathered, and every rank sends the faces whose bounding box overlaps them
  (PstreamBuffers; their global index and their points, in their own
  order).  An octree over these faces (indexedOctree of
  treeDataPrimitivePatch) gives the faces within the radius of a cell
  centre, and the distances are computed here.  Ties (a cell centre
  equidistant from several faces, e.g. on the diagonal of a corner of two
  patches, or nearest to an edge or a point shared by several faces;
  distances within 1e-10 relative of the smallest) are broken by the
  lexicographically smallest face centre, computed from the face's points
  in their own order: a boundary face keeps its points and their order in
  every decomposition, so the layer, the nearest faces and A_I/V_I depend
  on the decomposition only through the last bits of the cell centres of
  cells at processor faces (a cell exactly at the distance; the margins
  below show how far the rule is from that).  Collective.
  --------------------------------------------------------------------------*/
  const scalar searchRadius = 2*layerDistance_;
  const label myProc = Pstream::myProcNo();

  boundBox reach(C, false);
  if (reach.good()) {
    reach.grow(searchRadius);
  }
  List<boundBox> reaches(Pstream::nProcs());
  reaches[myProc] = reach;
  Pstream::allGatherList(reaches);

  /* the search set: the faces of this rank within reach, then those of the
     other ranks (their global indices and points) */
  DynamicList<label> setGlobal;
  DynamicList<pointField> setPoints;
  forAll(facePoints, i) {
    if (faceBoxes[i].overlaps(reach)) {
      setGlobal.append(faceNumbering.toGlobal(i));
      setPoints.append(facePoints[i]);
    }
  }

  if (Pstream::parRun()) {
    PstreamBuffers buffers(UPstream::commsTypes::nonBlocking);
    for (const int proci : Pstream::allProcs()) {
      if (proci == myProc) {
        continue;
      }
      DynamicList<label> ids;
      DynamicList<pointField> points;
      forAll(facePoints, i) {
        if (faceBoxes[i].overlaps(reaches[proci])) {
          ids.append(faceNumbering.toGlobal(i));
          points.append(facePoints[i]);
        }
      }
      UOPstream os(proci, buffers);
      os << labelList(ids) << List<pointField>(points);
    }
    buffers.finishedSends();
    for (const int proci : Pstream::allProcs()) {
      if (proci == myProc) {
        continue;
      }
      UIPstream is(proci, buffers);
      labelList ids;
      List<pointField> points;
      is >> ids >> points;
      setGlobal.append(ids);
      setPoints.append(points);
    }
  }

  /* the faces of the search set with their own points (a face list of
     local point labels), and their centres for the tie rule */
  label nSetPoints = 0;
  for (const pointField& points : setPoints) {
    nSetPoints += points.size();
  }
  pointField allPoints(nSetPoints);
  faceList allFaces(setPoints.size());
  pointField allCentres(setPoints.size());
  {
    label p = 0;
    forAll(setPoints, k) {
      face& f = allFaces[k];
      f.resize(setPoints[k].size());
      forAll(f, j) {
        allPoints[p] = setPoints[k][j];
        f[j] = p++;
      }
      allCentres[k] = f.centre(allPoints);
    }
  }

  /*--------------------------------------------------------------------------
  Every cell of this rank: its distance and nearest face, when a face lies
  within the search radius.  The layer is every cell with distance <=
  layerDistance_, in cell order (one element and one exchange cell per
  cell); the margin of the rule is the range of the distances of the
  selected cells and the smallest distance outside (within the radius).
  --------------------------------------------------------------------------*/
  DynamicList<label> cells, nearestFaces;
  scalar dMin = GREAT, dMax = -GREAT, dOutside = GREAT;

  if (allFaces.size()) {

    const primitivePatch searchPatch(SubList<face>(allFaces), allPoints);
    treeBoundBox domain(allPoints);
    domain.inflate(0.01);
    const indexedOctree<treeDataPrimitivePatch<primitivePatch>> tree (
      treeDataPrimitivePatch<primitivePatch> (
        false,
        searchPatch,
        indexedOctree<treeDataPrimitivePatch<primitivePatch>>::perturbTol()
      ),
      domain,
      8,
      10,
      3.0
    );

    const scalar radiusSqr = sqr(searchRadius);
    DynamicList<scalar> distances;

    forAll(C, celli) {

      const labelList hits(tree.findSphere(C[celli], radiusSqr));
      if (hits.empty()) {
        continue;
      }

      distances.resize(hits.size());
      scalar dCell = GREAT;
      forAll(hits, h) {
        distances[h] =
          allFaces[hits[h]].nearestPoint(C[celli], allPoints).distance();
        dCell = min(dCell, distances[h]);
      }

      /* the nearest face; ties by the smallest face centre */
      label nearest = -1;
      const scalar tied = dCell*(1 + 1e-10);
      forAll(hits, h) {
        if (distances[h] <= tied && (
              nearest < 0
           || lessPoint(allCentres[hits[h]], allCentres[nearest]))) {
          nearest = hits[h];
        }
      }

      if (dCell <= layerDistance_) {
        cells.append(celli);
        nearestFaces.append(setGlobal[nearest]);
        dMin = min(dMin, dCell);
        dMax = max(dMax, dCell);
      } else {
        dOutside = min(dOutside, dCell);
      }
    }
  }
  reduce(dMin, minOp<scalar>());
  reduce(dMax, maxOp<scalar>());
  reduce(dOutside, minOp<scalar>());

  const label n = cells.size();
  elementCell_.transfer(cells);
  elementPatch_.resize(n, -1);
  elementFace_.resize(n, -1);
  exchangeCells_ = elementCell_;
  cellStart_ = identity(n + 1);
  nWall_ = n;
  nWallCells_ = n;

  /*--------------------------------------------------------------------------
  The face map: the global face indices of the nearest faces become slots of
  the compact numbering (local faces first, then the remote faces this rank
  needs).  Collective.
  --------------------------------------------------------------------------*/
  labelList slots(nearestFaces);
  List<Map<label>> compactMap;
  faceMap_.reset(new mapDistributeBase(faceNumbering, slots, compactMap));
  elementFaceSlot_.transfer(slots);

  /* elements whose nearest face lies on another rank (for the log) */
  label nRemote = 0;
  for (const label slot : elementFaceSlot_) {
    if (slot >= nFaces) {
      ++nRemote;
    }
  }
  reduce(nRemote, sumOp<label>());

  /*--------------------------------------------------------------------------
  Interface faces that receive elements of two or more other ranks (for the
  log): their sums in nearestFaceSums() add the contributions of several
  ranks, in ascending rank order (pushToFaces()).  Every rank flags the
  remote slots it uses; the flags are summed per face by the same push, so
  a face counts the other ranks with elements on it, each once (review of
  M4, round 3: the scheduled reverse distribute of round 2 added a flag
  twice between two ranks that map elements to each other's faces, and the
  line counted faces with the elements of ONE other rank, e.g. 200 instead
  of 0 on a channel split by rows).  Collective.
  --------------------------------------------------------------------------*/
  label nMultiple = 0;
  {
    const mapDistributeBase& map = faceMap_();
    List<label> flags(map.constructSize(), Zero);
    for (const label slot : elementFaceSlot_) {
      if (slot >= nFaces) {
        flags[slot] = 1;
      }
    }
    for (const label nRanks : pushToFaces(map, nFaces, flags)) {
      if (nRanks >= 2) {
        ++nMultiple;
      }
    }
    reduce(nMultiple, sumOp<label>());
  }

  OStringStream report;
  report << "  layer distance (the exact distance of the cell centres from "
    << "the faces of the interface patches, all ranks): selected cells at "
    << dMin << " to " << dMax << " m, the nearest cell outside at ";
  if (dOutside < GREAT) {
    report << dOutside << " m";
  } else {
    report << "more than " << searchRadius << " m";
  }
  report << "; " << nRemote << " elements with their nearest wall face on "
    << "another rank, " << nMultiple << " interface faces with elements of "
    << "two or more other ranks" << nl;

  /*--------------------------------------------------------------------------
  Comparison with the centre-distance rule of T-Flows (plan section 12,
  item 10): R - r <= distance, r the distance of the cell centre from the
  axis, R the given radius.  Printed only; the layer is the set above.
  --------------------------------------------------------------------------*/
  if (d.found("centreDistance")) {
    const dictionary& cd = d.subDict("centreDistance");
    vector axis(cd.getOrDefault<vector>("axis", vector(1, 0, 0)));
    const point origin(cd.getOrDefault<point>("origin", Zero));
    const scalar radius = cd.get<scalar>("radius");
    if (!(mag(axis) > 0) || !(radius > 0)) {
      FatalIOErrorInFunction(cd) << "centreDistance needs an axis of "
        << "non-zero length and a positive radius." << exit(FatalIOError);
    }
    axis /= mag(axis);

    boolList inLayer(mesh_.nCells(), false);
    for (const label celli : elementCell_) {
      inLayer[celli] = true;
    }
    label nCentre = 0, nOnlyCentre = 0, nOnlyLayer = 0;
    scalar rMax = -GREAT, rOutside = GREAT;
    forAll(C, celli) {
      const vector offset = C[celli] - origin;
      const vector arm = offset - (offset & axis)*axis;
      const scalar wallDistance = radius - mag(arm);
      const bool centre = wallDistance <= layerDistance_;
      if (centre) {
        ++nCentre;
        rMax = max(rMax, wallDistance);
      } else {
        rOutside = min(rOutside, wallDistance);
      }
      nOnlyCentre += centre && !inLayer[celli];
      nOnlyLayer += !centre && inLayer[celli];
    }
    reduce(nCentre, sumOp<label>());
    reduce(nOnlyCentre, sumOp<label>());
    reduce(nOnlyLayer, sumOp<label>());
    reduce(rMax, maxOp<scalar>());
    reduce(rOutside, minOp<scalar>());

    report << "  centre-distance rule R - r <= " << layerDistance_
      << " m (R = " << radius << " m about the axis " << axis
      << " through " << origin << "; the T-Flows test, a comparison "
      << "only): " << nCentre << " cells (largest R - r selected " << rMax
      << " m, smallest outside " << rOutside << " m); " << nOnlyCentre
      << " cells only in this rule, " << nOnlyLayer << " only in the "
      << "layer" << nl;
  }

  return report.str();
}


/*------------------------------------------------------------------------------
layer distance: the values of the nearest faces at the elements (pull), and
the sums of element values per face (push: pushToFaces(), its own elements
first, then the other ranks in ascending rank order)
------------------------------------------------------------------------------*/

Foam::tmp<Foam::scalarField> LESTO::interfaceExchange::nearestFaceValues (
  const scalarField& faceValues
) const {

  /* local faces first, then the remote faces, in the compact numbering */
  List<scalar> slots(faceValues);
  faceMap_->distribute(slots);

  tmp<scalarField> tValues(new scalarField(nWall_));
  scalarField& values = tValues.ref();
  forAll(values, e) {
    values[e] = slots[elementFaceSlot_[e]];
  }
  return tValues;
}


Foam::tmp<Foam::scalarField> LESTO::interfaceExchange::nearestFaceSums (
  const scalarField& elementValues
) const {

  const mapDistributeBase& map = faceMap_();

  /* the sum of the elements of this rank per compact slot */
  List<scalar> slots(map.constructSize(), Zero);
  forAll(elementFaceSlot_, e) {
    slots[elementFaceSlot_[e]] += elementValues[e];
  }

  /*--------------------------------------------------------------------------
  The slots of remote faces go to the rank of the face and are added to its
  own slot there: its own elements first, then those of the other ranks in
  ascending rank order, every contribution once (pushToFaces()).  Not a
  reverse mapDistributeBase::distribute: reverseDistribute() assigns (a
  face with elements on several ranks would keep one contribution), the
  nonBlocking branch with plusEqOp adds them in the order the receives
  complete (mDep_ could differ in its last bit between identical runs;
  review of M4, round 2), and the scheduled branch adds every contribution
  twice between two ranks that map elements to each other's faces (review
  of M4, round 3: mDep_ +50 % on a channel split so; plan section 29).
  --------------------------------------------------------------------------*/
  return tmp<scalarField> (
    new scalarField(pushToFaces(map, interfaceFacePatch_.size(), slots))
  );
}


/*------------------------------------------------------------------------------
Temperature model: v_d (given, or from A [1/s] on this layer), k, Tdep
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readTemperatureModel() {

  const dictionary& c = dict_.subDict("temperatureCoeffs");
  const bool haveVelocity = c.found("depositionVelocity");
  const bool haveRate = c.found("A");

  if (haveVelocity == haveRate) {
    FatalIOErrorInFunction(c) << "Give exactly one of depositionVelocity "
      << "[m/s] or A [1/s]." << exit(FatalIOError);
  }

  const scalar layerAreaPerVolume = layerArea_/layerVolume_;

  if (haveVelocity) {
    depositionVelocity_ = c.get<scalar>("depositionVelocity");
  } else {
    depositionVelocity_ = c.get<scalar>("A")/layerAreaPerVolume;
  }
  k_ = c.get<scalar>("k");
  Tdep_ = c.get<scalar>("Tdep");

  if (!(depositionVelocity_ >= 0) || !(k_ > 0) || !(Tdep_ > 0)) {
    FatalIOErrorInFunction(c) << "Require depositionVelocity (or A) >= 0, "
      << "k > 0 and Tdep > 0." << exit(FatalIOError);
  }

  label nActive = 0;
  forAll(elementT_, e) {
    if (elementT_[e] < Tdep_) {
      ++nActive;
    }
  }
  reduce(nActive, sumOp<label>());

  /*--------------------------------------------------------------------------
  areaMultiplier scales the element areas A_e of every model (plan E2), so
  the rate applied to the gas of the layer is areaMultiplier*v_d*A_I/V_I
  (times f(T)), not the input A when areaMultiplier differs from 1.  With
  wallResistance halfCell this is the kinetic rate only; the rate applied
  with the half-cell resistance is printed per pair (setCoefficients()).
  --------------------------------------------------------------------------*/
  const bool halfCell = wallResistance_ == "halfCell";

  Info<< "Temperature model: v_d = " << depositionVelocity_ << " m/s";
  if (haveRate) {
    Info<< " (from A = " << c.get<scalar>("A") << " 1/s on this layer)";
  }
  Info<< ", k = " << k_ << " 1/K, Tdep = " << Tdep_ << " K; "
    << (halfCell ? "kinetic" : "applied") << " rate "
    << "areaMultiplier*v_d*A_I/V_I = "
    << areaMultiplier_*depositionVelocity_*layerAreaPerVolume
    << " 1/s (layer mean, times f(T)"
    << (halfCell ? "; the half-cell resistance reduces it: pair line" : "")
    << "); " << nActive << " of " << nElementsTotal_
    << " elements below Tdep" << nl;
}


/*------------------------------------------------------------------------------
Model HKS: HKSCoeffs (the p_eq tables are read with the pairs)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readHKSModel() {

  const dictionary c(dict_.subOrEmptyDict("HKSCoeffs"));

  accommodation_ = c.getOrDefault<scalar>("accommodation", 1);
  Ce_ = c.getOrDefault<scalar>("Ce", 1);
  kineticScale_ = c.getOrDefault<scalar>("kineticScale", 1);
  const word equilibrium(c.getOrDefault<word>("equilibrium", "table"));
  const word speciation(c.getOrDefault<word>("speciation", "none"));

  if (!(accommodation_ > 0 && accommodation_ <= 1) || !(Ce_ >= 0)
    || !(kineticScale_ >= 0) || !std::isfinite(Ce_)
    || !std::isfinite(kineticScale_)) {
    FatalIOErrorInFunction(dict_) << "HKSCoeffs: require 0 < accommodation "
      << "<= 1, Ce >= 0 and kineticScale >= 0 (finite), not accommodation "
      << accommodation_ << ", Ce " << Ce_ << ", kineticScale "
      << kineticScale_ << "." << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  The GEMS3K backend (M9; gemsEquilibrium.H).  Without the bridge
  requireBridge() stops with the rebuild hint.  speciation lagged takes the
  gas speciation of the local updates, so it needs equilibrium GEMS (and
  mode local, checked by gemsEquilibrium).
  --------------------------------------------------------------------------*/
  if (equilibrium == "GEMS") {
    gemsEquilibrium::requireBridge(dict_.subDict("HKSCoeffs"));
  }
  if (equilibrium != "table" && equilibrium != "GEMS") {
    FatalIOErrorInFunction(dict_) << "HKSCoeffs: unknown equilibrium "
      << equilibrium << ".  Choose table or GEMS." << exit(FatalIOError);
  }
  if (speciation != "none" && speciation != "lagged") {
    FatalIOErrorInFunction(dict_) << "HKSCoeffs: unknown speciation "
      << speciation << ".  Choose none or lagged (equilibrium GEMS, "
      << "GEMSCoeffs mode local)." << exit(FatalIOError);
  }
  if (speciation == "lagged" && equilibrium != "GEMS") {
    FatalIOErrorInFunction(dict_) << "HKSCoeffs: speciation lagged takes the "
      << "gas speciation from GEMS3K: it needs equilibrium GEMS with "
      << "GEMSCoeffs mode local." << exit(FatalIOError);
  }

  if (equilibrium == "GEMS") {
    if (!dict_.isDict("GEMSCoeffs")) {
      FatalIOErrorInFunction(dict_) << "equilibrium GEMS needs GEMSCoeffs "
        << "{ system \"<...>-dat.lst\"; mode frozen | local; carrier "
        << "\"<GEMS3K gas species>\"; ... } (plan section 4)."
        << exit(FatalIOError);
    }
    gems_.reset (
      new gemsEquilibrium (
        mesh_, dict_.subDict("GEMSCoeffs"), speciation == "lagged"
      )
    );
  }

  Info<< "HKS model: accommodation sigma = " << accommodation_ << ", Ce = "
    << Ce_ << ", kineticScale = " << kineticScale_
    << (kineticScale_ == 1 ? " (SI)" : "") << "; G = kineticScale Ce "
    << "2 sigma/(2 - sigma) sqrt(M/(2 pi R T_e)), p_eq from ";
  if (gems_) {
    Info<< "GEMS3K (equilibrium GEMS, mode "
      << (gems_->local() ? "local" : "frozen") << "), with the table of "
      << "every pair as its check and fallback; speciation " << speciation;
  } else {
    Info<< "the table of every pair";
  }
  Info<< nl;
}


/*------------------------------------------------------------------------------
The p_eq table of a pair (model HKS; vapourPressureTable.H)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readVapourPressure (
  pairData&         pair,
  const dictionary& pd
) {

  const word& gasName = speciesNames_[pair.gas];

  if (!pd.found("vapourPressure")) {
    FatalIOErrorInFunction(pd) << "Pair " << gasName << ": model HKS needs "
      << "vapourPressure { file; phases; units; } (optional: interpolation "
      << "logInverseT | linearT, outOfRange extrapolate | hold | fatal)."
      << exit(FatalIOError);
  }
  const dictionary& vp = pd.subDict("vapourPressure");

  fileName file(vp.get<fileName>("file"));
  file.expand();

  const wordList phaseWords(vp.get<wordList>("phases"));
  if (phaseWords.empty()) {
    FatalIOErrorInFunction(vp) << "Pair " << gasName << ": vapourPressure "
      << "lists no phase." << exit(FatalIOError);
  }
  std::vector<std::string> phases;
  for (const word& w : phaseWords) {
    phases.push_back(w);
  }

  try {
    pair.table.reset (
      new vapourPressureTable (
        vapourPressureTable::read (
          file,
          phases,
          vapourPressureTable::units(vp.get<word>("units")),
          vapourPressureTable::interpolation (
            vp.getOrDefault<word>("interpolation", "logInverseT")
          ),
          vapourPressureTable::outOfRange (
            vp.getOrDefault<word>("outOfRange", "extrapolate")
          )
        )
      )
    );
  } catch (const std::exception& error) {
    FatalIOErrorInFunction(vp) << "Pair " << gasName << ": " << error.what()
      << exit(FatalIOError);
  }

  const vapourPressureTable& table = pair.table();

  Info<< "Pair " << gasName << ": p_eq from "
    << vp.get<fileName>("file") << ", phases "
    << phaseWords << " (the minimum), " << label(table.nNodes())
    << " nodes " << table.Tmin() << "-" << table.Tmax() << " K; units "
    << vp.get<word>("units") << ", interpolation "
    << vp.getOrDefault<word>("interpolation", "logInverseT")
    << ", outOfRange " << vp.getOrDefault<word>("outOfRange", "extrapolate");
  for (const vapourPressureTable::crossing& x : table.crossings()) {
    Info<< "; stable phase " << table.phases()[x.below].c_str() << " -> "
      << table.phases()[x.above].c_str() << " at " << x.T << " K";
  }
  Info<< nl;
}


/*------------------------------------------------------------------------------
Samples: keyed to a pair; cells selected in memory or from a cellZone
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readSamples() {

  sampleCell_.resize(mesh_.nCells(), false);

  if (!dict_.found("samples")) {
    return;
  }

  dictionary expanded;
  for (const entry& entry : dict_.subDict("samples")) {
    if (!entry.isDict() || !entry.dict().found("amounts")) {
      expanded.add(entry);
      continue;
    }
    const dictionary& sd = entry.dict();
    if (sd.found("pair") || sd.found("amount") || sd.subDict("amounts").empty()) {
      FatalIOErrorInFunction(sd) << "Give nonempty amounts OR pair/amount."
        << exit(FatalIOError);
    }
    for (const Foam::entry& amount : sd.subDict("amounts")) {
      const word name(entry.keyword() + "__" + amount.keyword());
      if (dict_.subDict("samples").found(name) || expanded.found(name)) {
        FatalIOErrorInFunction(sd) << "Expanded sample name collision: " << name
          << exit(FatalIOError);
      }
      dictionary single(sd);
      single.remove("amounts");
      single.set("pair", amount.keyword());
      single.set("amount", sd.subDict("amounts").get<scalar>(amount.keyword()));
      expanded.add(name, single);
    }
  }
  dict_.set("samples", expanded);
  const dictionary& all = dict_.subDict("samples");
  const scalarField& V = mesh_.V();
  const vectorField& C = mesh_.C().primitiveField();
  const scalarField& T = thermo_.T().primitiveField();

  boolList layerCell(mesh_.nCells(), false);
  for (const label celli : exchangeCells_) {
    layerCell[celli] = true;
  }
  boolList& sampleCell = sampleCell_;
  List<boolList> pairSampleCell(pairs_.size());
  for (boolList& cells : pairSampleCell) cells.resize(mesh_.nCells(), false);

  DynamicList<sampleData> samples;

  /* the sample elements of the inventory samples, appended after the wall
     elements */
  DynamicList<label>  elementCells(elementCell_);
  DynamicList<scalar> elementAreas(elementArea_);
  DynamicList<scalar> elementTemperatures(elementT_);
  DynamicList<label>  elementSamples;
  DynamicList<label>  exchange(exchangeCells_);
  DynamicList<label>  start(cellStart_);

  for (const entry& e : all) {

    if (!e.isDict()) {
      FatalIOErrorInFunction(all) << "Sample " << e.keyword()
        << " must be a dictionary." << exit(FatalIOError);
    }

    const dictionary& sd = e.dict();
    sampleData s;
    s.name = e.keyword();

    const word pairName = sd.get<word>("pair");
    const label si = speciesNames_.find(pairName);
    if (si < 0 || pairOfSpecies_[si] < 0) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": pair "
        << pairName << " is not a paired gas species." << exit(FatalIOError);
    }
    s.pair = pairOfSpecies_[si];

    /*------------------------------------------------------------------------
    Cells: an existing cellZone, or a topoSet-like selection evaluated in
    memory (never topoSet on a symlinked mesh).
    ------------------------------------------------------------------------*/
    const bool haveSelection = sd.found("selection");
    const bool haveZone = sd.found("cellZone");
    if (haveSelection == haveZone) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": give exactly "
        << "one of selection or cellZone." << exit(FatalIOError);
    }

    if (haveZone) {
      const word zoneName = sd.get<word>("cellZone");
      const label zonei = mesh_.cellZones().findZoneID(zoneName);
      if (returnReduceAnd(zonei < 0)) {
        FatalIOErrorInFunction(sd) << "Sample " << s.name << ": cellZone "
          << zoneName << " not found." << exit(FatalIOError);
      }
      if (zonei >= 0) {
        s.cells = mesh_.cellZones()[zonei];
      }
    } else {
      s.cells = cellBitSet::select(mesh_, sd.subDict("selection")).sortedToc();
    }

    const word mode = sd.get<word>("mode");
    if (mode != "release" && mode != "inventory") {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": unknown mode "
        << mode << ". Choose release or inventory." << exit(FatalIOError);
    }
    s.inventory = mode == "inventory";

    if (s.inventory && model_ != modelType::HKS) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": mode "
        << "inventory needs model HKS (its reservoir evaporates by the HKS "
        << "law); with model " << dict_.get<word>("model") << " use mode "
        << "release." << exit(FatalIOError);
    }

    /*------------------------------------------------------------------------
    Keys of the other mode would be ignored silently (e.g. a removeTime of
    a release sample, which never removes anything): refused.
    ------------------------------------------------------------------------*/
    const wordList inventoryKeys ({
      "removeTime", "areaPerVolume", "kineticScale", "Ce"
    });
    const wordList releaseKeys({"startTime", "duration"});
    for (const word& key : s.inventory ? releaseKeys : inventoryKeys) {
      if (sd.found(key)) {
        FatalIOErrorInFunction(sd) << "Sample " << s.name << ": " << key
          << " belongs to mode " << (s.inventory ? "release" : "inventory")
          << ", not to mode " << mode << "; the sample would not use it."
          << nl << "Remove the entry." << exit(FatalIOError);
      }
    }

    s.amount = sd.get<scalar>("amount");
    if (!std::isfinite(s.amount)) {
      FatalIOErrorInFunction(sd) << "Sample amount must be finite."
        << exit(FatalIOError);
    }
    s.startTime = 0;
    s.duration = 0;
    s.length = 0;
    s.firstElement = elementCells.size();
    s.areaPerVolume = 0;
    s.kineticScale = 1;
    s.Ce = 1;
    s.removable = false;
    s.removeTime = 0;
    s.removed = false;
    s.removalClock = 0;
    s.volumeTolerance = 0;
    s.centreTolerance = 0;
    s.meshPoints = "binary";
    s.volumeLimit = 0;
    s.centreLimit = 0;

    if (s.inventory) {
      s.areaPerVolume = sd.get<scalar>("areaPerVolume");
      s.kineticScale = sd.getOrDefault<scalar>("kineticScale",
        pairs_[s.pair].pairHKS ? pairs_[s.pair].kineticScale : 1);
      s.Ce = sd.getOrDefault<scalar>("Ce",
        pairs_[s.pair].pairHKS ? pairs_[s.pair].Ce : 1);
      s.removable = sd.readIfPresent("removeTime", s.removeTime);

      if (!(s.amount >= 0) || !(s.areaPerVolume >= 0)
        || !(s.kineticScale >= 0) || !(s.Ce >= 0)
        || !std::isfinite(s.areaPerVolume) || !std::isfinite(s.Ce)
        || !std::isfinite(s.kineticScale)
        || (s.removable && !std::isfinite(s.removeTime))) {
        FatalIOErrorInFunction(sd) << "Sample " << s.name
          << ": require amount >= 0, areaPerVolume >= 0, kineticScale >= 0 "
          << "and Ce >= 0 (finite), and a finite removeTime."
          << exit(FatalIOError);
      }
    } else {
      s.startTime = sd.getOrDefault<scalar>("startTime", 0);
      s.duration = sd.get<scalar>("duration");
      s.length = releaseWindowLength(s.startTime, s.duration);

      if (!(s.amount >= 0) || !(s.duration > 0) || !(s.length > 0)) {
        FatalIOErrorInFunction(sd) << "Sample " << s.name
          << ": require amount >= 0 and duration > 0 (resolvable at the "
          << "startTime " << s.startTime << ")." << exit(FatalIOError);
      }
    }

    /*------------------------------------------------------------------------
    V_Z of the release rate and the fill: a plain sum reduced over the
    ranks.  For the restart record, V_Z and the volume-weighted centre as
    compensated sums over all ranks, which do not depend on the
    decomposition beyond the round-off of the cell volumes and centres
    (review of M5: the plain V_Z differs in its last bits after
    decomposePar, and the record compared it exactly).
    ------------------------------------------------------------------------*/
    label nInLayer = 0, nInOther = 0;
    s.volume = 0;
    compensatedSum volume;
    FixedList<compensatedSum, 3> moment;
    for (const label celli : s.cells) {
      if (layerCell[celli]) {
        ++nInLayer;
      }
      if (pairSampleCell[s.pair][celli]) {
        ++nInOther;
      }
      sampleCell[celli] = true;
      pairSampleCell[s.pair][celli] = true;
      s.volume += V[celli];
      volume.add(V[celli]);
      for (direction d = 0; d < vector::nComponents; ++d) {
        moment[d].add(V[celli]*C[celli][d]);
      }
    }
    reduce(nInLayer, sumOp<label>());
    reduce(nInOther, sumOp<label>());
    reduce(s.volume, sumOp<scalar>());
    s.recordVolume = globalSum(volume);
    for (direction d = 0; d < vector::nComponents; ++d) {
      s.centre[d] = globalSum(moment[d])/max(s.recordVolume, VSMALL);
    }
    s.nCells = returnReduce(s.cells.size(), sumOp<label>());

    if (s.nCells == 0) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name
        << " selects no cell." << exit(FatalIOError);
    }
    if (nInLayer > 0) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": " << nInLayer
        << " of its " << s.nCells << " cells lie in the interaction layer "
        << (
             distanceLayer()
           ? "(cells within the layer distance of the interface patches)"
           : "(cells of the interface patches)"
           )
        << ".  Keep samples and the layer disjoint." << exit(FatalIOError);
    }
    if (nInOther > 0) {
      FatalIOErrorInFunction(sd) << "Sample " << s.name << ": " << nInOther
        << " cells belong to another sample of the same pair." << exit(FatalIOError);
    }

    /*------------------------------------------------------------------------
    Inventory: one sample element per cell, of the interfacial area
    A_e = areaPerVolume V_c, at the cell temperature; each is its own
    exchange cell (the samples are disjoint from the layer and from each
    other).
    ------------------------------------------------------------------------*/
    if (s.inventory) {
      s.elements.resize(s.cells.size());
      forAll(s.elements, i) s.elements[i] = s.firstElement + i;
      for (const label celli : s.cells) {
        elementCells.append(celli);
        elementAreas.append(s.areaPerVolume*V[celli]);
        elementTemperatures.append(T[celli]);
        elementSamples.append(samples.size());
        exchange.append(celli);
        start.append(elementCells.size());
      }
    }

    Info<< "Sample " << s.name << " (pair " << pairName << "): " << s.nCells
      << " cells, V_Z = " << s.volume << " m3; ";
    if (s.inventory) {
      Info<< "inventory of " << s.amount << " mol (Cs = n0 M/V_Z = "
        << s.amount*pairs_[s.pair].molarMass/s.volume << " kg/m3 on a fresh "
        << "start), areaPerVolume " << s.areaPerVolume << " 1/m, "
        << "kineticScale " << s.kineticScale << ", Ce " << s.Ce << ", ";
      if (s.removable) {
        Info<< "removed at the step starting at t >= removeTime - deltaT/2, "
          << "removeTime " << s.removeTime << " s";
      } else {
        Info<< "never removed (no removeTime)";
      }
      Info<< nl;
    } else {
      Info<< "release of " << s.amount << " mol from " << s.startTime
        << " s over " << s.duration << " s" << nl;
    }

    samples.append(s);
  }

  std::map<label, DynamicList<label>> groups;
  DynamicList<label> unique;
  bool overlap = false;
  for (label e = nWall_; e < elementCells.size(); ++e) {
    const label celli = elementCells[e];
    auto found = groups.find(celli);
    if (found == groups.end()) unique.append(celli);
    else overlap = true;
    groups[celli].append(e);
  }
  if (overlap) {
    const labelList oldCells(elementCells);
    const scalarField oldAreas(elementAreas), oldT(elementTemperatures);
    const labelList oldSamples(elementSamples);
    labelList newIndex(elementCells.size(), -1);
    elementCells.resize(nWall_);
    elementAreas.resize(nWall_);
    elementTemperatures.resize(nWall_);
    elementSamples.clear();
    exchange.resize(nWallCells_);
    start.resize(nWallCells_ + 1);
    for (const label celli : unique) {
      exchange.append(celli);
      for (const label old : groups[celli]) {
        newIndex[old] = elementCells.size();
        elementCells.append(oldCells[old]);
        elementAreas.append(oldAreas[old]);
        elementTemperatures.append(oldT[old]);
        elementSamples.append(oldSamples[old - nWall_]);
      }
      start.append(elementCells.size());
    }
    for (sampleData& sample : samples) {
      for (label& e : sample.elements) e = newIndex[e];
    }
  }
  samples_.transfer(samples);

  /* a start list without wall elements is (0); it gains one entry per
     sample element */
  if (start.empty()) {
    start.append(0);
  }
  elementCell_.transfer(elementCells);
  elementArea_.transfer(elementAreas);
  elementT_.transfer(elementTemperatures);
  elementSample_.transfer(elementSamples);
  exchangeCells_.transfer(exchange);
  cellStart_.transfer(start);
  nSampleElementsTotal_ =
    returnReduce(elementCell_.size() - nWall_, sumOp<label>());
}


/*------------------------------------------------------------------------------
Start-up checks of a paired gas species (plan section 3)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::checkPairedSpecies(pairData& pair) {

  const volScalarField& Y = species_[pair.gas];
  const word& name = pair.solverName;

  /*--------------------------------------------------------------------------
  Implicit Euler: the predictor and the ledger assume it.  A conservative
  convection scheme: the ledger books the transport as face fluxes.
  --------------------------------------------------------------------------*/
  checkDdtScheme(pair, true);
  checkConvectionScheme(pair, true);

  /*--------------------------------------------------------------------------
  Zero species flux through the interface patches: the exchange enters the
  gas only through the source.
  --------------------------------------------------------------------------*/
  for (const label patchi : interfacePatches_) {
    const word& type = Y.boundaryField()[patchi].type();
    if (type != zeroGradientFvPatchScalarField::typeName) {
      FatalErrorInFunction
        << "Paired species " << name << " needs zeroGradient on the "
        << "interface patch " << mesh_.boundary()[patchi].name()
        << ", not " << type << "." << nl
        << "The exchange enters the gas through the source term only."
        << exit(FatalError);
    }
  }

  /*--------------------------------------------------------------------------
  The final solve (with <field>Final if it exists, otherwise the last
  <field> solve) must be converged to its absolute tolerance: relTol 0.
  --------------------------------------------------------------------------*/
  checkFinalSolver(pair, true);
  const word finalName(name + "Final");

  /*--------------------------------------------------------------------------
  Paired species are never relaxed.
  --------------------------------------------------------------------------*/
  for (const word& n : wordList({name, finalName})) {
    if (mesh_.solution().relaxEquation(n)) {
      WarningInFunction << "Equation relaxation factor "
        << mesh_.solution().equationRelaxationFactor(n) << " of " << n
        << " is ignored: paired species are never relaxed." << nl << endl;
    }
  }

  /*--------------------------------------------------------------------------
  The initial residual of an outer corrector is about the final residual of
  the previous solve, which stops at the absolute tolerance of <field>.  An
  outer tolerance below it is rarely reached: the loop then runs all
  nCorr + 1 correctors in (almost) every step.
  --------------------------------------------------------------------------*/
  const scalar solverTolerance =
    mesh_.solverDict(name).getOrDefault<scalar>("tolerance", 1e-6);

  if (tolerance_ < solverTolerance) {
    WarningInFunction << "The outer tolerance " << tolerance_
      << " (speciesTransport) is below the tolerance " << solverTolerance
      << " of the solver " << name << ": the outer loop of the paired "
      << "species will rarely converge and then runs all nCorr + 1 = "
      << nCorr_ + 1 << " correctors per step." << nl << endl;
  }

  /* the solvers entry used by the final solve, named if it is a pattern */
  const word finalUsed(pair.finalSolver ? finalName : name);
  const entry* finalEntry = solverEntry(finalUsed);

  Info<< "Pair " << speciesNames_[pair.gas] << " -> "
    << speciesNames_[pair.condensed] << ": M = " << pair.molarMass
    << " kg/mol; at most " << nCorr_ + 1 << " outer correctors to the "
    << "initial residual " << tolerance_ << "; final solve with "
    << finalUsed << (pair.finalSolver ? "" : " (relTol 0)");
  if (finalEntry && finalEntry->keyword() != finalUsed) {
    Info<< " (fvSolution entry " << finalEntry->keyword() << ")";
  }
  Info<< nl;
}


/*------------------------------------------------------------------------------
The solvers entry of fvSolution that OpenFOAM uses for a field name, as
solution::solverDict(): the literal keyword, or else the last regular
expression that matches it.  nullptr if there is none.
------------------------------------------------------------------------------*/

const Foam::entry* LESTO::interfaceExchange::solverEntry (
  const word& fieldName
) const {
  return mesh_.solution().solversDict().csearch(fieldName).ptr();
}


/*------------------------------------------------------------------------------
ddt scheme of a paired species (start-up and every step)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::checkDdtScheme (
  const pairData& pair,
  const bool      startup
) const {

  const word& name = pair.solverName;
  const word ddtName("ddt(" + rho_.name() + ',' + name + ')');
  ITstream& ddtStream = mesh_.ddtScheme(ddtName);
  ddtStream.rewind();
  const word ddtScheme(ddtStream);
  ddtStream.rewind();

  if (ddtScheme != "Euler") {
    FatalErrorInFunction
      << "Paired species " << name << " needs the Euler ddtScheme for "
      << ddtName << ", not " << ddtScheme;
    if (!startup) {
      FatalError << " (fvSchemes was modified during the run)";
    }
    FatalError << "." << nl
      << "The exchange predictor and the mass ledger assume implicit Euler."
      << exit(FatalError);
  }
}


/*------------------------------------------------------------------------------
Convection scheme of a paired species (start-up and every step)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::checkConvectionScheme (
  const pairData& pair,
  const bool      startup
) const {

  /*--------------------------------------------------------------------------
  The ledger books the transport of the gas as the face fluxes of the solved
  matrix (fvMatrix::flux() on the non-coupled patches; the internal faces
  telescope).  The Gauss convection scheme is such a sum of face fluxes for
  every interpolation scheme.  'bounded Gauss ...' (boundedConvectionScheme)
  subtracts fvm::Sp(fvc::div(phi),Y): a cell source, not a face flux, which
  the ledger does not book; with a carrier whose phi is not exactly
  divergence-free it creates or destroys the species (the run/012 pipe
  closed to 8.7e-12 instead of 3e-16 in 5 steps).  Every other convection
  scheme is refused as well: only Gauss is known to be a sum of face
  fluxes.  The scheme name is the one of solvePairedGasSpecies.H.
  --------------------------------------------------------------------------*/
  const word& name = pair.solverName;
  const word divName("div(" + phi_.name() + ',' + name + ')');
  ITstream& divStream = mesh_.divScheme(divName);
  divStream.rewind();
  const word convectionScheme(divStream);
  const std::string schemeText(divStream.toString());
  divStream.rewind();

  if (convectionScheme != "Gauss") {
    FatalIOErrorInFunction(mesh_.divSchemes())
      << "Paired species " << name << " needs a conservative convection "
      << "scheme 'Gauss <interpolation>' for " << divName << ", not '"
      << schemeText.c_str() << "'";
    if (!startup) {
      FatalIOError << " (fvSchemes was modified during the run)";
    }
    FatalIOError << "." << nl;
    if (convectionScheme == "bounded") {
      FatalIOError << "'bounded' subtracts fvm::Sp(fvc::div(phi), " << name
        << "), a cell source that is not a face flux: the mass ledger does "
        << "not book it, and where phi is not exactly divergence-free it "
        << "creates or destroys the species.";
    } else {
      FatalIOError << "Only the Gauss convection scheme is known to be a sum "
        << "of face fluxes, which the mass ledger books as transport.";
    }
    FatalIOError << nl << "Use Gauss with any interpolation scheme, e.g. "
      << "'Gauss SuperBee'." << exit(FatalIOError);
  }
}


/*------------------------------------------------------------------------------
Final solver dictionary of a pair (start-up and every step)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::checkFinalSolver (
  pairData&  pair,
  const bool startup
) {

  /*--------------------------------------------------------------------------
  A separate Final solve exists when <field>Final selects another solvers
  entry than <field>: a literal Y_PbI2_gFinal, or a pattern such as
  "Y_.*Final" defined after "Y_.*".  A single pattern "Y_.*" matches both
  names; it gives no separate Final solve (dictionary::found(), a pattern
  match, took it for one: one redundant solve and a 'final predictor' per
  step, review round 4).
  --------------------------------------------------------------------------*/
  const word& name = pair.solverName;
  const word finalName(name + "Final");
  const entry* finalEntry = solverEntry(finalName);
  const bool finalSolver = finalEntry && finalEntry != solverEntry(name);

  if (!startup && finalSolver != pair.finalSolver) {
    Info<< name << ": fvSolution was modified, the final solve now uses "
      << (finalSolver ? finalName : name) << nl;
  }
  pair.finalSolver = finalSolver;

  const dictionary& finalDict =
    mesh_.solverDict(finalSolver ? finalName : name);
  const scalar relTol = finalDict.getOrDefault<scalar>("relTol", 0);

  if (relTol != 0) {
    FatalIOErrorInFunction(finalDict)
      << "The final solve of the paired species " << name << " uses "
      << (finalSolver ? finalName : name) << " with relTol " << relTol
      << "; it needs relTol 0";
    if (!startup) {
      FatalIOError << " (fvSolution was modified during the run)";
    }
    FatalIOError << "." << nl
      << "Either add a solver " << finalName << " with relTol 0 (loose "
      << name << " correctors plus a tight final solve), or set relTol 0 "
      << "in " << name << ".";
    if (!finalSolver && finalEntry && finalEntry->keyword().isPattern()) {
      FatalIOError << nl << finalName << " selects the solvers entry "
        << finalEntry->keyword() << ", as " << name << " does: of several "
        << "matching patterns OpenFOAM uses the last one, so a pattern for "
        << "the Final solve must follow " << finalEntry->keyword() << ".";
    }
    FatalIOError << exit(FatalIOError);
  }
}


/*------------------------------------------------------------------------------
Element coefficients a_e, g_e (T is frozen: once)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::setCoefficients(pairData& pair) const {

  const scalarField& rhoI = rho_.primitiveField();
  const scalarField& Tc = thermo_.T().primitiveField();
  const scalarField& pc = thermo_.p().primitiveField();
  const scalarField& V = mesh_.V();
  const volScalarField& rhoDi = rhoD_[pair.gas];
  const bool HKS = model_ == modelType::HKS;
  const scalar M = pair.molarMass;
  const scalar R = universalGasConstant;

  /*--------------------------------------------------------------------------
  halfCell: the half-cell distance d_e = 1/deltaCoeffs of the wall face.
  The delta coefficients are built on first use, collectively (lduAddr,
  globalData), so they are obtained here on every rank, also on a rank
  without interface elements (see Parallel hygiene above).
  --------------------------------------------------------------------------*/
  const surfaceScalarField* deltaCoeffs =
    wallResistance_ == "halfCell" ? &mesh_.deltaCoeffs() : nullptr;

  /* the half cell of wall element e: d_e/(rhoD)_e [m2 s/kg], and whether
     rhoD > 0 on its face */
  auto halfCell = [&](const label e, bool& conducting) -> scalar {
    const label patchi = elementPatch_[e];
    const label facei = elementFace_[e];
    const scalar distance = 1.0/deltaCoeffs->boundaryField()[patchi][facei];
    const scalar rhoDface = rhoDi.boundaryField()[patchi][facei];
    conducting = rhoDface > 0;
    return conducting ? distance/rhoDface : 0;
  };

  if (!HKS) {

    /* sum_e A_e v and sum_e A_e v' [m3/s], for the halfCell report */
    scalar kineticFlow = 0, appliedFlow = 0;

    for (label e = 0; e < nWall_; ++e) {
      const label celli = elementCell_[e];
      scalar v =
        depositionVelocity_*depositionActivation(elementT_[e], k_, Tdep_);
      kineticFlow += elementArea_[e]*v;

      if (deltaCoeffs) {
        /* no diffusion to the wall (rhoD = 0, e.g. D 0): an infinite
           resistance of the half cell, v' = 0, without dividing by zero */
        bool conducting = false;
        const scalar resistance = halfCell(e, conducting);
        v = conducting ? halfCellVelocity(v, rhoI[celli]*resistance) : 0;
      }

      appliedFlow += elementArea_[e]*v;
      pair.a[e] = 0;
      pair.g[e] = elementArea_[e]*rhoI[celli]*v;
    }

    /*------------------------------------------------------------------------
    halfCell: the rate applied to the gas of the layer, sum_e A_e v'_e/V_I
    with f(T) and the half-cell resistance, next to the kinetic rate
    sum_e A_e v_d f_e/V_I of the same elements (the model line above gives
    the kinetic rate per unit f).  Collective: every rank calls
    setCoefficients().
    ------------------------------------------------------------------------*/
    if (deltaCoeffs) {
      reduce(kineticFlow, sumOp<scalar>());
      reduce(appliedFlow, sumOp<scalar>());
      Info<< "Pair " << speciesNames_[pair.gas] << ", wallResistance "
        << "halfCell: applied rate sum_e A_e v'_e/V_I = "
        << appliedFlow/layerVolume_ << " 1/s (layer mean with f(T) and the "
        << "half-cell resistance; kinetic sum_e A_e v_d f_e/V_I = "
        << kineticFlow/layerVolume_ << " 1/s)" << nl;
    }
    return;
  }

  /*--------------------------------------------------------------------------
  HKS (plan E2).  Wall element e in cell c:

    G_e    = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T_e)),
    beta_c = rho_c R T_c/M,  a_e = A_e G_e p_eq(T_e),  g_e = A_e G_e beta_c;

  halfCell: G'_e = G_e/(1 + G_e beta_c d_e/(rhoD)_e), 0 without diffusion
  to the wall.  Sample element (one per cell of an inventory sample of this
  pair): the same with A_e = areaPerVolume V_c, T_e = T_c and the sample's
  own kineticScale and Ce; the elements of the samples of other pairs stay
  inactive (a = g = 0).  p_eq outside the table follows outOfRange; fatal
  is checked for all ranks first, so that every rank stops.
  --------------------------------------------------------------------------*/
  const vapourPressureTable& table = pair.table();
  const word& gasName = speciesNames_[pair.gas];

  label nOutside = 0;
  scalar Tlo = GREAT, Thi = -GREAT;
  forAll(elementCell_, e) {
    const bool own = e < nWall_ || samples_[elementSample_[e - nWall_]].pair
      == pairOfSpecies_[pair.gas];
    if (own && !table.inRange(elementT_[e])) {
      ++nOutside;
      Tlo = min(Tlo, elementT_[e]);
      Thi = max(Thi, elementT_[e]);
    }
  }
  reduce(nOutside, sumOp<label>());
  reduce(Tlo, minOp<scalar>());
  reduce(Thi, maxOp<scalar>());

  if (nOutside > 0) {
    string range;
    {
      OStringStream os;
      os << nOutside << " element temperatures from " << Tlo << " to " << Thi
        << " K lie outside the table [" << table.Tmin() << ", "
        << table.Tmax() << "] K";
      range = os.str();
    }
    try {
      table(Tlo < table.Tmin() ? Tlo : Thi);
    } catch (const std::exception& error) {
      FatalIOErrorInFunction(dict_.subDict("pairs").subDict(gasName))
        << "Pair " << gasName << ": " << range.c_str() << " (outOfRange "
        << "fatal)." << nl << error.what() << exit(FatalIOError);
    }
  }

  /* report: ranges over the wall elements, the layer rate and the samples */
  scalar TeMin = GREAT, TeMax = -GREAT, pMin = GREAT, pMax = -GREAT;
  scalar YeqMin = GREAT, YeqMax = -GREAT, velocityMin = GREAT;
  scalar velocityMax = -GREAT;
  label nBoiling = 0;
  scalar layerRate = 0, layerMass = 0;
  scalar lambdaMin = GREAT, lambdaMax = -GREAT;
  List<FixedList<scalar, 4>> sampleRate(samples_.size());
  for (FixedList<scalar, 4>& r : sampleRate) {
    r[0] = 0;       /* sum g_e [kg/s] */
    r[1] = 0;       /* sum rho V [kg] */
    r[2] = GREAT;   /* min Y_eq */
    r[3] = -GREAT;  /* max Y_eq */
  }

  /*--------------------------------------------------------------------------
  The per-element data of the law (interfaceExchange.H, sixth comment
  block): the kinetic conductance, the half cell, p_eq of the table and of
  the law (equilibrium GEMS, mode frozen: the GEMS table of the pair; mode
  local: the table until the first update), chi = 1.
  --------------------------------------------------------------------------*/
  const label nElements = elementCell_.size();
  const bool frozenGems = gems_ && !gems_->local();
  pair.conductance.resize(nElements, 0.0);
  pair.resistance.resize(nElements, 0.0);
  pair.pEqTable.resize(nElements, 0.0);
  pair.pEq.resize(nElements, 0.0);
  pair.chi.resize(nElements, 1.0);
  pair.source.resize(nElements, label(-1));
  pair.chi = 1.0;
  pair.source = label(-1);

  forAll(elementCell_, e) {

    const label celli = elementCell_[e];
    const bool wall = e < nWall_;
    const label samplei = wall ? -1 : elementSample_[e - nWall_];

    if (!wall && samples_[samplei].pair != pairOfSpecies_[pair.gas]) {
      pair.a[e] = 0;
      pair.g[e] = 0;
      continue;
    }

    const scalar Te = elementT_[e];
    const scalar beta = rhoI[celli]*R*Tc[celli]/M;
    pair.pEqTable[e] = table(Te);
    pair.pEq[e] = frozenGems
      ? gems_->frozenPEq(pairOfSpecies_[pair.gas], Te, pair.pEqTable[e],
                         pair.source[e])
      : pair.pEqTable[e];
    const scalar pEq = pair.pEq[e];

    pair.conductance[e] = wall
      ? hksConductance(M, Te, pair.accommodation, pair.Ce, pair.kineticScale)
      : hksConductance(M, Te, pair.accommodation, samples_[samplei].Ce,
                       samples_[samplei].kineticScale);

    if (wall && deltaCoeffs) {
      bool conducting = false;
      const scalar resistance = halfCell(e, conducting);
      pair.resistance[e] = conducting ? resistance : -1;
    }

    const scalar G = hksCoefficients(pair, e, deltaCoeffs != nullptr);

    const scalar Yeq = pEq/beta;
    if (wall) {
      TeMin = min(TeMin, Te);
      TeMax = max(TeMax, Te);
      pMin = min(pMin, pEq);
      pMax = max(pMax, pEq);
      YeqMin = min(YeqMin, Yeq);
      YeqMax = max(YeqMax, Yeq);
      velocityMin = min(velocityMin, G*beta/rhoI[celli]);
      velocityMax = max(velocityMax, G*beta/rhoI[celli]);
      layerRate += pair.g[e];
      if (pEq > pc[celli]) {
        ++nBoiling;
      }
    } else {
      sampleRate[samplei][0] += pair.g[e];
      sampleRate[samplei][1] += rhoI[celli]*V[celli];
      sampleRate[samplei][2] = min(sampleRate[samplei][2], Yeq);
      sampleRate[samplei][3] = max(sampleRate[samplei][3], Yeq);
      if (pEq > pc[celli]) {
        ++nBoiling;
      }
    }
  }

  /* lambda of the layer cells: sum_{e in c} g_e/(rho_c V_c) */
  for (label k = 0; k < nWallCells_; ++k) {
    const label celli = exchangeCells_[k];
    scalar g = 0;
    for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
      g += pair.g[e];
    }
    const scalar lambda = g/(rhoI[celli]*V[celli]);
    lambdaMin = min(lambdaMin, lambda);
    lambdaMax = max(lambdaMax, lambda);
    layerMass += rhoI[celli]*V[celli];
  }

  reduce(TeMin, minOp<scalar>());
  reduce(TeMax, maxOp<scalar>());
  reduce(pMin, minOp<scalar>());
  reduce(pMax, maxOp<scalar>());
  reduce(YeqMin, minOp<scalar>());
  reduce(YeqMax, maxOp<scalar>());
  reduce(velocityMin, minOp<scalar>());
  reduce(velocityMax, maxOp<scalar>());
  reduce(lambdaMin, minOp<scalar>());
  reduce(lambdaMax, maxOp<scalar>());
  reduce(layerRate, sumOp<scalar>());
  reduce(layerMass, sumOp<scalar>());
  reduce(nBoiling, sumOp<label>());

  Info<< "Pair " << gasName << ", HKS on the interface: T_e " << TeMin << "-"
    << TeMax << " K, p_eq " << pMin << "-" << pMax << " Pa, Y_eq = "
    << "p_eq/beta_c " << YeqMin << "-" << YeqMax << "; G_e beta_c/rho_c "
    << velocityMin << "-" << velocityMax << " m/s"
    << (deltaCoeffs ? " (with the half cell)" : "") << "; lambda = "
    << "sum_e g_e/(rho_c V_c) " << lambdaMin << "-" << lambdaMax
    << " 1/s, layer mean " << layerRate/max(layerMass, VSMALL) << " 1/s; "
    << nOutside << " elements outside the table ("
    << (nOutside > 0 ? "extrapolated or held, outOfRange" : "none") << "); "
    << nBoiling << " elements with p_eq > p (the condensate would boil; "
    << "not capped)" << nl;

  forAll(samples_, samplei) {
    const sampleData& s = samples_[samplei];
    if (!s.inventory || s.pair != pairOfSpecies_[pair.gas]) {
      continue;
    }
    FixedList<scalar, 4>& r = sampleRate[samplei];
    reduce(r[0], sumOp<scalar>());
    reduce(r[1], sumOp<scalar>());
    reduce(r[2], minOp<scalar>());
    reduce(r[3], maxOp<scalar>());
    Info<< "Sample " << s.name << ", HKS: lambda_s = sum g_e/(rho V) "
      << r[0]/max(r[1], VSMALL) << " 1/s (mean over its cells), Y_eq "
      << r[2] << "-" << r[3] << nl;
  }

  /*--------------------------------------------------------------------------
  equilibrium GEMS, mode frozen: the elements of each source of p_eq, and
  the largest deviation of the elements of the GEMS table from the pair's
  table (plan section 6, M9 acceptance 2).
  --------------------------------------------------------------------------*/
  if (frozenGems) {
    FixedList<label, gemsEquilibrium::nSourceTypes> nSource(label(0));
    scalar maxDeviation = 0, TgemsMin = GREAT, TgemsMax = -GREAT;
    forAll(elementCell_, e) {
      if (e >= nWall_
        && samples_[elementSample_[e - nWall_]].pair
        != pairOfSpecies_[pair.gas]) {
        continue;
      }
      ++nSource[pair.source[e]];
      if (pair.source[e] == gemsEquilibrium::FROM_GEMS) {
        maxDeviation = max (
          maxDeviation,
          mag(std::log10(pair.pEq[e]) - std::log10(pair.pEqTable[e]))
        );
        TgemsMin = min(TgemsMin, elementT_[e]);
        TgemsMax = max(TgemsMax, elementT_[e]);
      }
    }
    for (label& n : nSource) {
      reduce(n, sumOp<label>());
    }
    reduce(maxDeviation, maxOp<scalar>());
    reduce(TgemsMin, minOp<scalar>());
    reduce(TgemsMax, maxOp<scalar>());

    Info<< "Pair " << gasName << ", GEMS frozen on the elements: "
      << nSource[gemsEquilibrium::FROM_GEMS] << " from the GEMS table";
    if (nSource[gemsEquilibrium::FROM_GEMS] > 0) {
      Info<< " (T_e " << TgemsMin << "-" << TgemsMax << " K; max |dlog10 "
        << "p_eq/p_eq,table| " << maxDeviation << ")";
    }
    Info<< ", from the table " << nSource[gemsEquilibrium::BELOW_TEMPERATURE]
      << " below minTemperature and "
      << nSource[gemsEquilibrium::OUTSIDE_GRID] << " outside the GEMS "
      << "table" << nl;
  }
}


/*------------------------------------------------------------------------------
HKS law of an element from its p_eq and chi (interfaceExchange.H, sixth
comment block).  With chi = 1 the arithmetic is that of M5: (A G') p_eq and
(A G') beta, G' = G/(1 + G beta d/(rhoD)) with the half cell (0 without
diffusion to the wall), so equilibrium table gives the same doubles.
------------------------------------------------------------------------------*/

Foam::scalar LESTO::interfaceExchange::hksCoefficients (
  pairData&   pair,
  const label e,
  const bool  halfCell
) const {

  const label celli = elementCell_[e];
  const scalar beta = rho_.primitiveField()[celli]*universalGasConstant
                    *thermo_.T().primitiveField()[celli]/pair.molarMass;
  const scalar chiBeta = pair.chi[e]*beta;

  scalar G = pair.conductance[e];
  if (halfCell && e < nWall_) {
    G = pair.resistance[e] >= 0
      ? halfCellConductance(G, chiBeta, pair.resistance[e])
      : 0;
  }

  pair.a[e] = elementArea_[e]*G*pair.pEq[e];
  pair.g[e] = elementArea_[e]*G*chiBeta;
  return G;
}


/*------------------------------------------------------------------------------
equilibrium GEMS, mode local: the update of p_eq and chi of every element
(beginStep(), at the start of the steps of the schedule of the time index:
gemsEquilibrium::updateDue()).  The composition is that of the cell at the
start of the step: the carrier p/(R T) and c_k = rho Y_k/M_k of every pair.
A sample element belongs to its own pair only; every pair's gas enters the
bulk.  The rank's elements are evaluated locally, the counts reduced at the
end (collective).
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::updateEquilibrium() {

  gemsEquilibrium& gems = gems_();
  const label nPairs = pairs_.size();
  const scalarField& rhoI = rho_.primitiveField();
  const scalarField& pc = thermo_.p().primitiveField();
  const bool halfCell = wallResistance_ == "halfCell";

  scalarList c(nPairs), pEqTable(nPairs), pEq(nPairs), chi(nPairs);
  labelList source(nPairs);
  boolList active(nPairs);

  gems.beginUpdate();

  forAll(elementCell_, e) {

    const label celli = elementCell_[e];
    const label own =
      e < nWall_ ? -1 : samples_[elementSample_[e - nWall_]].pair;

    forAll(pairs_, k) {
      const pairData& pair = pairs_[k];
      active[k] = own < 0 || own == k;
      c[k] = rhoI[celli]*species_[pair.gas][celli]/pair.molarMass;
      pEqTable[k] = pair.pEqTable[e];
      pEq[k] = pair.pEq[e];
      chi[k] = pair.chi[e];
      source[k] = pair.source[e];
    }

    gems.evaluate (
      e, elementT_[e], pc[celli], carrierConcentration_[celli], c, active,
      pEqTable, pEq, chi, source
    );

    forAll(pairs_, k) {
      if (active[k]) {
        pairData& pair = pairs_[k];
        pair.pEq[e] = pEq[k];
        pair.chi[e] = chi[k];
        pair.source[e] = source[k];
        hksCoefficients(pair, e, halfCell);
      }
    }
  }

  gems.endUpdate();
}


/*------------------------------------------------------------------------------
pEq_<gas>, gemsSource_<gas> and chi_<gas> (writeEquilibriumField): the
element of each face of the interface patches (layer faceCells), the mean
(the largest source) of the elements of an exchange cell, 0 (source -1)
elsewhere; processor patches from the neighbour cells
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::updateEquilibriumFields(pairData& pair) const {

  const label pairi = pairOfSpecies_[pair.gas];
  volScalarField& pEqField = pair.equilibriumField();
  volScalarField& sourceField = pair.sourceField();
  volScalarField* chiField = pair.chiField.get();

  pEqField == dimensionedScalar(pEqField.dimensions(), Zero);
  sourceField == dimensionedScalar(dimless, -1);
  if (chiField) {
    *chiField == dimensionedScalar(dimless, Zero);
  }

  forAll(exchangeCells_, k) {
    const label celli = exchangeCells_[k];
    scalar sumP = 0, sumChi = 0;
    label n = 0, source = -1;
    for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
      if (e >= nWall_ && samples_[elementSample_[e - nWall_]].pair != pairi) {
        continue;
      }
      sumP += pair.pEq[e];
      sumChi += pair.chi[e];
      source = max(source, pair.source[e]);
      ++n;
    }
    if (n > 0) {
      pEqField[celli] = sumP/n;
      sourceField[celli] = source;
      if (chiField) {
        (*chiField)[celli] = sumChi/n;
      }
    }
  }

  for (label e = 0; e < nWall_; ++e) {
    const label patchi = elementPatch_[e];
    if (patchi < 0) {
      continue;
    }
    const label facei = elementFace_[e];
    pEqField.boundaryFieldRef()[patchi][facei] = pair.pEq[e];
    sourceField.boundaryFieldRef()[patchi][facei] = pair.source[e];
    if (chiField) {
      chiField->boundaryFieldRef()[patchi][facei] = pair.chi[e];
    }
  }

  pEqField.boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  sourceField.boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  if (chiField) {
    chiField->boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  }
}


/*------------------------------------------------------------------------------
Format of a field file at its instance.  An uncollated file states it in its
header ('format').  With the collated file handler every field file is a
decomposedBlockData container whose header always says 'format binary' (the
container); the format of the data inside is 'data.format', or, in a
container without that entry, the 'format' of the header of its first
block.  Every rank must call this function: the master-reading file
handlers (masterUncollated, collated) open the stream collectively.
------------------------------------------------------------------------------*/

bool LESTO::interfaceExchange::writtenBinary (
  const IOobject& io,
  const word&     typeName
) {

  const fileName path(io.localFilePath(typeName, false));
  if (returnReduceOr(path.empty())) {
    return false;
  }

  autoPtr<ISstream> isPtr(fileHandler().NewIFstream(path));
  IOobject header(io);
  dictionary headerDict;
  if (!isPtr || !isPtr->good() || !header.readHeader(headerDict, *isPtr)) {
    return false;
  }

  if (!decomposedBlockData::isCollatedType(header)) {
    return isPtr->format() == IOstreamOption::BINARY;
  }

  word dataFormat;
  if (headerDict.readIfPresent("data.format", dataFormat)) {
    return IOstreamOption::formatEnum(dataFormat) == IOstreamOption::BINARY;
  }

  List<char> block;
  decomposedBlockData::readBlockEntry(*isPtr, block);
  ISpanStream blockStream(block);
  IOobject blockHeader(io);
  return blockHeader.readHeader(blockStream)
      && blockStream.format() == IOstreamOption::BINARY;
}


/*------------------------------------------------------------------------------
Whether this run can read a per-rank list of uniform/phaseChange/.
Uncollated, every rank reads its own processor<i>/<time>/uniform/
phaseChange/<name>, a plain list.  With the collated file handler the
per-rank lists of a parallel run are decomposedBlockData containers, one
block per rank: one container per time (processors<N>/<time>/...) with one
IO rank, one per group of ranks (processors<N>_<lo>-<hi>/<time>/..., the
ranks lo to hi) with several IO ranks (-ioRanks, FOAM_IORANKS).  Each rank
reads its own block of the container of its group.  reconstructPar and
decomposePar copy uniform/ verbatim: such a container reaches a serial
case, or processors<M>/ of another number of ranks, where the file handler
stops the run ('could not detect processor number', or a missing block;
review of M5, round 3).  A container is therefore read only in a parallel
run of as many ranks as its directory names, and only by the ranks it was
written for: as many ranks as it has blocks, all of them with one IO rank,
the ranks of its group (and of this run's IO group) with several (review
of M5, round 4: the size of the whole run was compared, so every
restart under -ioRanks read the lists as foreign).  The fingerprint at the
start of the list then decides whether it belongs to these elements
(another decomposition of the same number of ranks, a renumbered mesh).
A file that some ranks lack is FOREIGN as well: the ranks without it could
not fit.  A file whose header cannot be read is left to the reader, which
reports it.  Collective: every rank calls it.
------------------------------------------------------------------------------*/

LESTO::interfaceExchange::stateFile LESTO::interfaceExchange::perRankFile (
  const IOobject& io,
  const word&     typeName,
  string&         why
) {

  const fileName path(io.localFilePath(typeName, false));
  if (!returnReduceOr(!path.empty())) {
    return stateFile::ABSENT;
  }
  if (!returnReduceAnd(!path.empty())) {
    why = "missing on some ranks";
    return stateFile::FOREIGN;
  }

  /* the number of blocks of a collated container, -1 for a plain list */
  label nBlocks = -1;
  {
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(path));
    IOobject header(io);
    if (
        isPtr && isPtr->good() && header.readHeader(*isPtr)
     && decomposedBlockData::isCollatedType(header)
    ) {
      nBlocks = decomposedBlockData::getNumBlocks(*isPtr);
    }
  }

  bool readable = nBlocks < 0;
  if (!readable && UPstream::parRun()) {
    fileName root, processorsDir, local;
    fileOperation::procRangeType group;
    label nProcs = -1;
    fileOperation::splitProcessorPath (
      path, root, processorsDir, local, group, nProcs
    );
    readable = nProcs == UPstream::nProcs() && (
      group.empty()
    ? nBlocks == UPstream::nProcs()
    : nBlocks == group.size() && group.contains(UPstream::myProcNo())
   && nBlocks == UPstream::nProcs(fileHandler().comm())
    );
  }
  if (returnReduceAnd(readable)) {
    return stateFile::READABLE;
  }

  reduce(nBlocks, maxOp<label>());
  OStringStream os;
  os << "a container of the collated file handler written on " << nBlocks
     << " ranks, which ";
  if (UPstream::parRun()) {
    os << "this run on " << UPstream::nProcs() << " ranks";
  } else {
    os << "a serial run";
  }
  os << " cannot read (reconstructPar and decomposePar copy uniform/ "
     << "verbatim)";
  why = os.str();
  return stateFile::FOREIGN;
}


/*------------------------------------------------------------------------------
State at the start time: the ledger, the element deposits and regimes, and
the restart booking
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::readState() {

  const Time& runTime = mesh_.time();
  const word startName = runTime.timeName();
  const label nPairs = pairs_.size();
  const label nElements = elementCell_.size();
  const label nSample = nElements - nWall_;
  const bool haveInventory = nSampleElementsTotal_ > 0;

  wordList gasNames(nPairs), condensateNames(nPairs);
  forAll(pairs_, pairi) {
    gasNames[pairi] = speciesNames_[pairs_[pairi].gas];
    condensateNames[pairi] = speciesNames_[pairs_[pairi].condensed];
  }
  ledger_.reset(new phaseChangeLedger(mesh_, gasNames, condensateNames));
  phaseChangeLedger& ledger = *ledger_;

  forAll(pairs_, pairi) {
    pairs_[pairi].transport.resize(ledger.patchIds().size(), 0.0);
  }

  /*--------------------------------------------------------------------------
  The record of the rounded 'uniform' values of the files at the start time
  (recordUniformValues(), in uniform/phaseChangeLayout), read before it is
  replaced by the record of this run.
  --------------------------------------------------------------------------*/
  bool haveUniformRecord = false;
  wordList checkedRecord, roundedRecord;
  if (ledger.restarted()) {
    const dictionary* record =
      ledger.layout().findDict("uniformValues", keyType::LITERAL);
    if (record) {
      haveUniformRecord = true;
      checkedRecord = record->get<wordList>("checked");
      roundedRecord = record->get<wordList>("rounded");
    }
  }

  /* on a restart the interface must be the one the deposits were kept on */
  recordInterface();

  /*--------------------------------------------------------------------------
  Per-rank element state in uniform/phaseChange/ (binaryIOList.H): the
  regimes (review item B5) of all elements, the deposits of the wall
  elements and, with inventory samples, the reservoirs of the sample
  elements, each behind the fingerprint of those elements of this rank.
  Not in uniform/ itself: the master-reading file handlers
  (masterUncollated, and collated, which derives from it) resolve every
  object whose local directory is exactly "uniform" to the master's file
  and broadcast that path (masterUncollatedFileOperation::filePath(): they
  take uniform/ for global data), so every rank read rank 0's list (or the
  container of rank 0's IO group), the fingerprints refused it, and every
  restart on the same decomposition re-read the rounded mw_ and Cs_ (review
  of M5, round 4).  A list is read only when this run can read it
  (perRankFile()): a collated container of another number of ranks is
  present, but written for other elements, and its state comes from mw_
  and Cs_ like that of another decomposition.
  --------------------------------------------------------------------------*/
  fingerprint_ = elementsFingerprint(0);
  elementsFingerprint_ = elementsFingerprint(1);
  samplesFingerprint_ = elementsFingerprint(2);
  const label nHeader = fingerprint_.size();
  const label nState = nHeader + nPairs*nWall_;
  const label nRegimeState = nHeader + nPairs*nElements;
  const label nReservoirState = nHeader + nPairs*nSample;

  auto stateObject = [&](const word& name, const word& typeName,
                         stateFile& file, string& why) {
    IOobject io(name, startName, stateLocal, runTime, IOobject::NO_READ,
                IOobject::AUTO_WRITE);
    file = perRankFile(io, typeName, why);
    if (file == stateFile::READABLE) {
      io.readOpt(IOobject::READ_IF_PRESENT);
    }
    return io;
  };

  stateFile regimesFile = stateFile::ABSENT;
  stateFile depositsFile = stateFile::ABSENT;
  stateFile reservoirsFile = stateFile::ABSENT;
  string regimesWhy, depositsWhy, reservoirsWhy;

  regimesIO_.reset (
    new binaryIOList<label> (
      stateObject("phaseChangeRegimes", IOList<label>::typeName,
                  regimesFile, regimesWhy)
    )
  );
  depositsIO_.reset (
    new binaryIOList<scalar> (
      stateObject("phaseChangeDeposits", IOList<scalar>::typeName,
                  depositsFile, depositsWhy)
    )
  );
  if (haveInventory) {
    reservoirsIO_.reset (
      new binaryIOList<scalar> (
        stateObject("phaseChangeReservoirs", IOList<scalar>::typeName,
                    reservoirsFile, reservoirsWhy)
      )
    );
  }

  /* equilibrium GEMS, mode local: the lagged state of the last update
     (restoreEquilibrium(); gemsEquilibrium.H, mode local) */
  stateFile equilibriumFile = stateFile::ABSENT;
  string equilibriumWhy;
  if (gems_ && gems_->local()) {
    equilibriumIO_.reset (
      new binaryIOList<scalar> (
        stateObject("phaseChangeEquilibrium", IOList<scalar>::typeName,
                    equilibriumFile, equilibriumWhy)
      )
    );
  }
  const labelList& regimes = *regimesIO_;
  const scalarList& deposits = *depositsIO_;
  const bool haveDeposits =
    depositsFile == stateFile::FOREIGN || returnReduceOr(!deposits.empty());
  const bool haveReservoirs = haveInventory && (
    reservoirsFile == stateFile::FOREIGN
 || returnReduceOr(!reservoirsIO_().empty())
  );

  /*--------------------------------------------------------------------------
  Deposit fields mw_<condensate> (output): read at the start time if
  present, otherwise created with zero deposits.
  --------------------------------------------------------------------------*/
  label nPresent = 0;
  wordList fieldNames(nPairs);

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];
    fieldNames[pairi] = "mw_" + speciesNames_[pair.condensed];

    IOobject io(fieldNames[pairi],
                startName,
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE);

    const bool here = io.typeHeaderOk<volScalarField>(true, false);
    const bool present = returnReduceOr(here);

    if (present != returnReduceAnd(here)) {
      FatalErrorInFunction << fieldNames[pairi] << " exists at time "
        << startName << " in some processor directories only."
        << exit(FatalError);
    }

    if (present) {
      ++nPresent;
      pair.depositField.reset(new volScalarField(io, mesh_));
      if (pair.depositField().dimensions() != dimMass/dimArea) {
        FatalErrorInFunction << fieldNames[pairi]
          << " must have the dimensions of kg/m2." << exit(FatalError);
      }
    } else {
      io.readOpt(IOobject::NO_READ);
      pair.depositField.reset (
        new volScalarField (
          io,
          mesh_,
          dimensionedScalar(dimMass/dimArea, Zero),
          calculatedFvPatchScalarField::typeName
        )
      );
    }
  }

  /*--------------------------------------------------------------------------
  Output fields that are never read (NO_READ: a file of the same name at
  the start time is ignored and overwritten): mDep_<condensate>, the
  deposit per wall area, and with writeExchangeField exch_<gas>, the
  realised exchange of the last step.
  --------------------------------------------------------------------------*/
  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];

    pair.wallDepositField.reset (
      new volScalarField (
        IOobject("mDep_" + speciesNames_[pair.condensed],
                 startName,
                 mesh_,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE),
        mesh_,
        dimensionedScalar(dimMass/dimArea, Zero),
        calculatedFvPatchScalarField::typeName
      )
    );

    if (writeExchangeField_) {
      pair.exchangeField.reset (
        new volScalarField (
          IOobject("exch_" + speciesNames_[pair.gas],
                   startName,
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar(dimDensity/dimTime, Zero),
          calculatedFvPatchScalarField::typeName
        )
      );
    }
  }

  /*--------------------------------------------------------------------------
  Reservoir fields Cs_<condensate> [kg/m3] of the pairs with inventory
  samples (output): read at the start time if present (the fallback of the
  reservoirs, as mw_ for the deposits), otherwise created with zero.
  --------------------------------------------------------------------------*/
  label nReservoirPairs = 0, nReservoirFields = 0;
  wordList reservoirNames;

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];
    if (!pair.inventory) {
      continue;
    }
    ++nReservoirPairs;
    const word name("Cs_" + speciesNames_[pair.condensed]);
    reservoirNames.append(name);

    IOobject io(name, startName, mesh_, IOobject::MUST_READ,
                IOobject::AUTO_WRITE);
    const bool here = io.typeHeaderOk<volScalarField>(true, false);
    const bool present = returnReduceOr(here);

    if (present != returnReduceAnd(here)) {
      FatalErrorInFunction << name << " exists at time " << startName
        << " in some processor directories only." << exit(FatalError);
    }

    if (present) {
      ++nReservoirFields;
      pair.reservoirField.reset(new volScalarField(io, mesh_));
      if (pair.reservoirField().dimensions() != dimDensity) {
        FatalErrorInFunction << name << " must have the dimensions of "
          << "kg/m3." << exit(FatalError);
      }
    } else {
      io.readOpt(IOobject::NO_READ);
      pair.reservoirField.reset (
        new volScalarField (
          io,
          mesh_,
          dimensionedScalar(dimDensity, Zero),
          calculatedFvPatchScalarField::typeName
        )
      );
    }
  }

  /*--------------------------------------------------------------------------
  Fresh start = no ledger and no condensate state.  The state without the
  ledger is refused: a restart directory that lost uniform/ must never
  re-fill or silently re-initialise the balance.
  --------------------------------------------------------------------------*/
  if (
      !ledger.restarted()
   && (nPresent > 0 || haveDeposits || nReservoirFields > 0 || haveReservoirs)
  ) {
    FatalErrorInFunction
      << "Time " << startName << " holds phase-change state (";
    bool first = true;
    auto separator = [&first]() {
      const char* text = first ? "" : ", ";
      first = false;
      return text;
    };
    if (nPresent > 0) {
      FatalError<< separator() << "the deposit field(s) " << fieldNames;
    }
    if (haveDeposits) {
      FatalError<< separator() << "uniform/phaseChange/phaseChangeDeposits";
    }
    if (nReservoirFields > 0) {
      FatalError<< separator() << "the reservoir field(s) " << reservoirNames;
    }
    if (haveReservoirs) {
      FatalError<< separator() << "uniform/phaseChange/phaseChangeReservoirs";
    }
    FatalError<< ") but not the ledger uniform/phaseChangeBalance." << nl
      << "This is neither a fresh start (no ledger, no deposits) nor a "
      << "complete restart.  Copy the uniform/ directory of the run that "
      << "wrote the state, or remove the state for a fresh start."
      << exit(FatalError);
  }

  /*--------------------------------------------------------------------------
  Clock of the release schedule: the start time on a fresh start; on a
  restart the exact time kept in the ledger.  OpenFOAM restarts from the
  time-directory name, which is rounded to timePrecision; the ledger's
  clock is not.
  --------------------------------------------------------------------------*/
  const scalar dt = runTime.deltaTValue();

  if (!ledger.restarted()) {
    ledger.time() = runTime.value();
  } else {
    const scalar offset = runTime.value() - ledger.time();
    if (mag(offset) > 1e-6*dt) {
      WarningInFunction
        << "The start time " << startName << " differs from the time "
        << phaseChangeLedger::exactDecimal(ledger.time()) << " s of the "
        << "phase-change ledger by " << offset << " s (" << offset/dt
        << " deltaT)." << nl
        << "The time-directory names are rounded (raise timePrecision), or "
        << "uniform/ belongs to another time.  The release schedule "
        << "continues on the clock of the ledger." << nl << endl;
    }
  }

  /*--------------------------------------------------------------------------
  The samples against the record of the run that wrote the state (and
  their record for this run), on the clock just set and before the restart
  booking below: a sample dropped or changed at a restart is named as such,
  instead of failing the round-off bound of the inventory re-read.
  --------------------------------------------------------------------------*/
  recordSamples();

  /*--------------------------------------------------------------------------
  Element state of a restart.  The deposits come from uniform/phaseChange/
  phaseChangeDeposits when it was written for the elements of every rank
  (exact); otherwise, e.g. after decomposePar, reconstructPar or
  renumberMesh, from the interface-patch values of mw_ (a warning).  The
  regimes decide whether the first corrector of the next step counts
  regime changes, so a restart continues with them to repeat the
  continuous run bit for bit; they are restored only onto the same
  elements.
  --------------------------------------------------------------------------*/
  bool depositsFromList = false;

  /* the value of mw_ that holds the deposit of wall element e: its face on
     the interface patch (layer faceCells), or its cell (layer distance) */
  auto depositFieldValue = [this](const volScalarField& mw, const label e) {
    return distanceLayer()
      ? mw[elementCell_[e]]
      : mw.boundaryField()[elementPatch_[e]][elementFace_[e]];
  };

  /* true when the values of mw_ that hold the deposits of this rank may be
     a rounded 'uniform' entry (one repeated non-zero value): the interface
     patches, or (layer distance) the internal field */
  auto depositFieldRepeated = [this](const volScalarField& mw) {
    if (distanceLayer()) {
      return repeatedNonZero(mw.primitiveField());
    }
    bool repeated = false;
    for (const label patchi : interfacePatches_) {
      repeated = repeated || repeatedNonZero(mw.boundaryField()[patchi]);
    }
    return repeated;
  };

  if (ledger.restarted()) {

    bool fits = deposits.size() == nState;
    for (label i = 0; fits && i < nHeader; ++i) {
      fits = deposits[i] == scalar(fingerprint_[i]);
    }
    depositsFromList = returnReduceAnd(fits);

    if (depositsFromList) {

      forAll(pairs_, pairi) {
        pairData& pair = pairs_[pairi];
        forAll(pair.deposit, e) {
          pair.deposit[e] = deposits[nHeader + pairi*nWall_ + e];
        }
      }

    } else if (nPresent == nPairs) {

      WarningInFunction
        << "uniform/phaseChange/phaseChangeDeposits at time " << startName
        << " is "
        << (
             depositsFile == stateFile::FOREIGN ? depositsWhy
           : haveDeposits ? string("written for other interface elements "
               "(another decomposition or a renumbered mesh?)")
           : string("missing")
           ).c_str()
        << ".  The deposits are re-read from the "
        << (distanceLayer() ? "layer cells" : "interface patches") << " of "
        << fieldNames << "." << nl << endl;

      forAll(pairs_, pairi) {
        pairData& pair = pairs_[pairi];
        forAll(pair.deposit, e) {
          pair.deposit[e] = depositFieldValue(pair.depositField(), e);
        }
      }

    } else {

      FatalErrorInFunction
        << "Time " << startName << " holds the ledger uniform/"
        << "phaseChangeBalance, but neither uniform/phaseChange/"
        << "phaseChangeDeposits of these interface elements nor all of the "
        << "deposit fields "
        << fieldNames << "." << nl
        << "A restart needs the element deposits written with the ledger."
        << exit(FatalError);
    }

    bool valid = regimes.size() == nRegimeState;
    for (label i = 0; valid && i < nHeader; ++i) {
      valid = regimes[i] == elementsFingerprint_[i];
    }
    for (label i = nHeader; valid && i < regimes.size(); ++i) {
      valid = regimes[i] >= UNKNOWN && regimes[i] <= CAPPED;
    }
    if (returnReduceAnd(valid)) {
      forAll(pairs_, pairi) {
        forAll(pairs_[pairi].regime, e) {
          pairs_[pairi].regime[e] = regimes[nHeader + pairi*nElements + e];
        }
      }
    } else {
      WarningInFunction
        << "uniform/phaseChange/phaseChangeRegimes at time " << startName
        << " is "
        << (
             regimesFile == stateFile::FOREIGN ? regimesWhy
           : string("missing or belongs to other interface elements "
               "(another decomposition, a renumbered mesh, or inventory "
               "samples added or dropped?)")
           ).c_str()
        << ".  The first corrector counts every active element as a "
        << "regime change." << nl << endl;
    }

    if (equilibriumIO_) {
      restoreEquilibrium(equilibriumFile, equilibriumWhy, startName);
    }
  }

  /*--------------------------------------------------------------------------
  Reservoirs of the inventory samples.  Fresh start: Cs = n0 M/V_Z in every
  cell of the sample (plan E7), part of the initial inventory.  Restart:
  from uniform/phaseChange/phaseChangeReservoirs when it was written for
  the sample elements of every rank (exact); otherwise from the fields
  Cs_<condensate> (a warning, e.g. after decomposePar, reconstructPar or
  renumberMesh).
  When every inventory sample has been removed, the reservoirs are 0 and
  need neither: e.g. a removed sample that was dropped at a restart (its
  record is kept, recordSamples()) and is listed again later, after runs
  without inventory samples wrote neither the list nor Cs_.
  --------------------------------------------------------------------------*/
  bool reservoirsFromList = true;
  bool reservoirsRemoved = false;

  if (haveInventory && !ledger.restarted()) {

    for (const sampleData& s : samples_) {
      if (!s.inventory) {
        continue;
      }
      pairData& pair = pairs_[s.pair];
      const scalar Cs = s.amount*pair.molarMass/s.volume;
      forAll(s.cells, i) {
        pair.reservoir[s.elements[i] - nWall_] = Cs;
      }
      Info<< "Fresh start of sample " << s.name << ": reservoir Cs = n0 M/V_Z"
        << " = " << Cs << " kg/m3 of " << speciesNames_[pair.condensed]
        << " in its " << s.nCells << " cells (" << s.amount << " mol)" << nl;
    }

  } else if (haveInventory) {

    const scalarList& reservoirs = *reservoirsIO_;
    bool fits = reservoirs.size() == nReservoirState;
    for (label i = 0; fits && i < nHeader; ++i) {
      fits = reservoirs[i] == scalar(samplesFingerprint_[i]);
    }
    reservoirsFromList = returnReduceAnd(fits);

    if (reservoirsFromList) {

      forAll(pairs_, pairi) {
        pairData& pair = pairs_[pairi];
        forAll(pair.reservoir, k) {
          pair.reservoir[k] = reservoirs[nHeader + pairi*nSample + k];
        }
      }

    } else if (nReservoirFields == nReservoirPairs) {

      WarningInFunction
        << "uniform/phaseChange/phaseChangeReservoirs at time " << startName
        << " is "
        << (
             reservoirsFile == stateFile::FOREIGN ? reservoirsWhy
           : haveReservoirs ? string("written for other sample elements "
               "(another decomposition, a renumbered mesh, or inventory "
               "samples reordered or dropped?)")
           : string("missing")
           ).c_str()
        << ".  The reservoirs are re-read from " << reservoirNames << "."
        << nl << endl;

      forAll(pairs_, pairi) {
        pairData& pair = pairs_[pairi];
        if (!pair.inventory) {
          continue;
        }
        for (label k = 0; k < nSample; ++k) {
          if (samples_[elementSample_[k]].pair == pairi) {
            pair.reservoir[k] =
              pair.reservoirField()[elementCell_[nWall_ + k]];
          }
        }
      }

    } else if (
      std::all_of(samples_.begin(), samples_.end(),
                  [](const sampleData& s) { return !s.inventory || s.removed; })
    ) {

      reservoirsRemoved = true;
      forAll(pairs_, pairi) {
        pairs_[pairi].reservoir = 0;
      }
      Info<< "Every inventory sample was removed by the run that wrote the "
        << "state at time " << startName << ": the reservoirs are 0 (no "
        << "uniform/phaseChange/phaseChangeReservoirs or " << reservoirNames
        << " needed)" << nl;

    } else {

      FatalErrorInFunction
        << "Time " << startName << " holds the ledger uniform/"
        << "phaseChangeBalance, but neither uniform/phaseChange/"
        << "phaseChangeReservoirs of these sample elements nor all of the "
        << "reservoir fields "
        << reservoirNames << "." << nl
        << "A restart with inventory samples needs the reservoirs written "
        << "with the ledger." << exit(FatalError);
    }
  }

  /*--------------------------------------------------------------------------
  Inventory: the initial value on a fresh start; on a restart, the
  difference between the inventory re-read and the one held by the ledger
  is booked as 'restart'.
  --------------------------------------------------------------------------*/
  const scalar roundingBound =
    std::pow(10.0, 1.0 - scalar(IOstream::defaultPrecision()));

  /*--------------------------------------------------------------------------
  How the state of a restart was stored.  Every call below is collective
  (writtenBinary() for the master-reading file handlers, the reductions;
  the number of reductions never depends on the data of a rank).

  The inventory re-read is sum rho Y V plus the deposits.  It is exact when
  Y_<gas> is written in binary and its internal field is no rounded
  'uniform' entry on any rank, the same for rho when it is re-read, and the
  deposits come from uniform/phaseChange/phaseChangeDeposits (or from an
  mw_ in binary without a rounded interface patch); the difference from the
  ledger must then be round-off.

  The next steps also use values that the inventory does not contain (the
  review of the final M2 pass: the verdict ignored them): Y_<gas> on the
  non-coupled patches (fixedValue, inletOutlet: the inflow and the
  diffusion through the patch; a patch whose value is recomputed when it is
  read, e.g. zeroGradient, fixedGradient, mixed, symmetryPlane or cyclic,
  is not concerned: valueReadBack()) and, with writeFrozenFields yes, where
  the carrier
  is read from the restart directory, the patches of rho (rhoD on the
  boundary), T and p (element temperatures of interfaceTemperature wall,
  the PbI2He diffusivity) and phi, including its processor patches, which
  are not evaluated.  A rounded one does not change the inventory, so the
  round-off bound still holds, but the restart does not repeat the
  continuous run bit for bit: the verdict is 'inventory exact, continuation
  rounded', and the values are listed.  With writeFrozenFields no the
  carrier is the one the run started from, the same in both runs.

  What the verdict covers: the inventory, the internal fields and the
  patch VALUES that are read back.  Other boundary parameters that
  OpenFOAM writes as 'uniform' text as well (inletValue of inletOutlet,
  refValue/refGradient/valueFraction of mixed, gradient of fixedGradient)
  are not checked; the restart line says so.

  Which 'uniform' entries are rounded is known only to the run that wrote
  them (roundedUniform); it records them, merged over its ranks
  (recordUniformValues).  The record is used when it belongs to these
  files: present and the element deposits of uniform/ fit the elements (the
  same decomposition and mesh; not after decomposePar, reconstructPar or
  renumberMesh, which rewrite the fields).  Otherwise, and for a field the
  record does not cover, every repeated non-zero value counts as rounded
  ('possibly rounded').
  --------------------------------------------------------------------------*/
  const bool useRecord =
    returnReduceAnd(haveUniformRecord) && depositsFromList;

  /* the rounded values of a field on the ranks: its internal field, and the
     names of the patches ('processor' for all processor patches, only
     checked when processorPatches) */
  auto roundedValues = [&](const auto& field, const bool processorPatches,
                           bool& internal, string& patches) {

    const word& name = field.name();
    const bool fromRecord = useRecord && checkedRecord.found(name);
    auto rounded = [&](const word& item, const scalarField& values) {
      return fromRecord ? roundedRecord.found(item) : repeatedNonZero(values);
    };

    internal = returnReduceOr(rounded(name, field.primitiveField()));
    patches.clear();

    const auto& bf = field.boundaryField();
    const label nNonProcessor = mesh_.boundaryMesh().nNonProcessor();

    for (label patchi = 0; patchi < nNonProcessor; ++patchi) {
      const word& patchName = bf[patchi].patch().name();
      const bool here =
        valueReadBack(bf[patchi])
     && rounded(word(name + ':' + patchName), bf[patchi]);
      if (returnReduceOr(here)) {
        patches += (patches.empty() ? "" : " ") + patchName;
      }
    }

    if (processorPatches) {
      bool here = false;
      if (fromRecord) {
        here = roundedRecord.found(word(name + ":processor"));
      } else {
        for (label patchi = nNonProcessor; patchi < bf.size(); ++patchi) {
          here = here || repeatedNonZero(bf[patchi]);
        }
      }
      if (returnReduceOr(here)) {
        patches += (patches.empty() ? "" : " ") + string("processor");
      }
    }

    return returnReduceAnd(fromRecord);
  };

  /* "<name> (written in ascii, a 'uniform' internal field, 'uniform' on
     <patches>)" for what applies; empty when nothing is rounded */
  auto roundedItem = [](const word& name, const bool binary,
                        const bool internal, const string& patches) {
    string what;
    if (!binary) {
      what = "written in ascii";
    }
    if (internal) {
      what += string(what.empty() ? "" : ", ") + "a 'uniform' internal field";
    }
    if (!patches.empty()) {
      what += string(what.empty() ? "" : ", ") + "'uniform' on " + patches;
    }
    if (what.empty()) {
      return what;
    }
    return string(name + " (" + what + ")");
  };

  auto append = [](string& list, const string& item) {
    if (!item.empty()) {
      list += (list.empty() ? "" : ", ") + item;
    }
  };

  const bool carrierReread =
    ledger.restarted() && returnReduceOr(rho_.instance() == startName);
  bool carrierBinary = true, carrierRepeated = false;
  bool carrierRecorded = true;
  string carrierRounded;

  if (carrierReread) {

    bool internal = false;
    string patches;

    carrierBinary = returnReduceAnd (
      writtenBinary(IOobject(rho_.name(), startName, mesh_))
    );
    carrierRecorded = roundedValues(rho_, false, carrierRepeated, patches);
    append(carrierRounded, roundedItem(rho_.name(), true, false, patches));

    for (const volScalarField* fieldPtr : {&thermo_.T(), &thermo_.p()}) {
      const volScalarField& field = *fieldPtr;
      const bool binary = returnReduceAnd (
        writtenBinary(IOobject(field.name(), startName, mesh_))
      );
      carrierRecorded =
        roundedValues(field, false, internal, patches) && carrierRecorded;
      append(carrierRounded, roundedItem(field.name(), binary, internal,
                                         patches));
    }

    const bool phiBinary = returnReduceAnd (
      writtenBinary (
        IOobject(phi_.name(), startName, mesh_),
        surfaceScalarField::typeName
      )
    );
    carrierRecorded =
      roundedValues(phi_, true, internal, patches) && carrierRecorded;
    append(carrierRounded, roundedItem(phi_.name(), phiBinary, internal,
                                       patches));
  }

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];
    const word& gasName = speciesNames_[pair.gas];
    const volScalarField& Y = species_[pair.gas];
    const volScalarField& mw = pair.depositField();
    const scalar held =
      gasInventory(pair) + wallInventory(pair) + sampleInventory(pair);

    if (!ledger.restarted()) {

      ledger.set(pairi, phaseChangeLedger::INITIAL, held);
      ledger.set(pairi, phaseChangeLedger::HELD, held);

      Info<< "Fresh start of pair " << gasName << ": initial inventory "
        << held << " kg (" << held/pair.molarMass << " mol)" << nl;

      const volScalarField& Ys = species_[pair.condensed];
      if (gMax(mag(Ys.primitiveField())) > 0) {
        WarningInFunction << "The initial values of " << Ys.name()
          << " are ignored: outside model mock it is derived output, "
          << "(deposit*area + reservoir*V)/(rho*V), recomputed every step."
          << nl << endl;
      }

      continue;
    }

    const scalar adjustment = held - ledger(pairi, phaseChangeLedger::HELD);
    ledger.add(pairi, phaseChangeLedger::RESTART, adjustment);
    ledger.set(pairi, phaseChangeLedger::HELD, held);

    const scalar reference = max(ledger.reference(pairi), VSMALL);
    const scalar relative = adjustment/reference;

    /* the gas: its internal field (inventory) and its non-coupled patches
       (the next steps) */
    const bool gasBinary =
      returnReduceAnd(writtenBinary(IOobject(Y.name(), startName, mesh_)));
    bool gasRepeated = false;
    string gasPatches;
    const bool recorded =
      roundedValues(Y, false, gasRepeated, gasPatches) && carrierRecorded;

    string continuation(roundedItem(Y.name(), true, false, gasPatches));
    append(continuation, carrierRounded);

    bool depositsBinary = true, depositsRepeated = false;
    if (!depositsFromList) {
      depositsBinary =
        returnReduceAnd(writtenBinary(IOobject(mw.name(), startName, mesh_)));
      depositsRepeated = returnReduceOr(depositFieldRepeated(mw));
    }

    /* the reservoirs: exact from the list or when all are removed (0),
       otherwise as mw_ */
    bool reservoirsBinary = true, reservoirsRepeated = false;
    if (pair.inventory && !reservoirsFromList && !reservoirsRemoved) {
      const volScalarField& Cs = pair.reservoirField();
      reservoirsBinary =
        returnReduceAnd(writtenBinary(IOobject(Cs.name(), startName, mesh_)));
      reservoirsRepeated = returnReduceOr(repeatedNonZero(Cs.primitiveField()));
    }

    const bool exact =
      gasBinary && !gasRepeated && carrierBinary && !carrierRepeated
   && depositsBinary && !depositsRepeated && reservoirsBinary
   && !reservoirsRepeated;

    auto describe = [](const word& name, const bool binary,
                       const bool repeated) -> string {
      return name + " written in " + (binary ? "binary" : "ascii")
        + (repeated ? ", a 'uniform' value on some rank (text at "
                      "writePrecision)" : "");
    };

    /* the sources of the inventory, e.g. "gas: Y_PbI2_g written in binary;
       carrier: rho written in binary; deposits: uniform/phaseChange/
       phaseChangeDeposits" (no carrier item when it is not re-read) */
    string sources("gas: " + describe(Y.name(), gasBinary, gasRepeated));
    if (carrierReread) {
      sources += "; carrier: "
        + describe(rho_.name(), carrierBinary, carrierRepeated);
    }
    sources += "; deposits: " + (
      depositsFromList
        ? string("uniform/phaseChange/phaseChangeDeposits")
        : describe(mw.name(), depositsBinary, depositsRepeated)
    );
    if (pair.inventory) {
      sources += "; reservoirs: " + (
        reservoirsFromList ? string("uniform/phaseChange/phaseChangeReservoirs")
      : reservoirsRemoved ? string("0, every inventory sample removed")
      : describe(pair.reservoirField().name(), reservoirsBinary,
                 reservoirsRepeated)
      );
    }

    const char* verdict = !exact ? "rounded"
      : continuation.empty() ? "exact"
      : "inventory exact, continuation rounded";

    Info<< "Restart of pair " << gasName << ": the inventory re-read "
      << "differs from the ledger by " << adjustment << " kg (relative "
      << relative << "; " << sources.c_str() << "; ";
    if (!continuation.empty()) {
      Info<< (recorded ? "" : "possibly ") << "rounded for the next steps: "
        << continuation.c_str() << "; ";
    }
    Info<< verdict << "), booked as restart; checked: the inventory, the "
      << "internal fields and the patch values re-read, not other boundary "
      << "parameters (e.g. inletValue)" << nl;

    if (exact && mag(relative) > 1e-14) {
      FatalErrorInFunction
        << "Restart of pair " << gasName << ": the inventory at time "
        << startName << " is stored exactly (" << sources.c_str()
        << "), but it differs from the ledger by " << relative
        << " (relative), more than the round-off bound 1e-14." << nl
        << "The files of time " << startName << " and its uniform/ do not "
        << "belong to the same run and time." << exit(FatalError);
    }

    if (!exact && mag(relative) > roundingBound) {
      WarningInFunction
        << "Restart of pair " << gasName << ": the rounded state at time "
        << startName << " (" << sources.c_str() << ") changes the "
        << "inventory by " << relative << " (relative), more than "
        << "10^(1 - writePrecision) = " << roundingBound << ".  The "
        << "difference is booked as restart; use writeFormat binary (and "
        << "writeFrozenFields no) for exact restarts." << nl << endl;
    }

    /*------------------------------------------------------------------------
    mw_ is output: with the deposits of uniform/ it is only compared (a
    warning if it differs by more than its rounding) and rewritten from
    them.
    ------------------------------------------------------------------------*/
    if (depositsFromList && nPresent == nPairs) {

      scalar maxDifference = 0, maxDeposit = 0;
      forAll(pair.deposit, e) {
        const scalar value = depositFieldValue(mw, e);
        maxDifference = max(maxDifference, mag(value - pair.deposit[e]));
        maxDeposit = max(maxDeposit, mag(pair.deposit[e]));
      }
      const bool repeated = depositFieldRepeated(mw);
      reduce(maxDifference, maxOp<scalar>());
      reduce(maxDeposit, maxOp<scalar>());

      const bool fieldExact =
        returnReduceAnd(writtenBinary(IOobject(mw.name(), startName, mesh_)))
     && !returnReduceOr(repeated);
      const scalar tolerance = fieldExact ? 0 : roundingBound*maxDeposit;

      if (maxDifference > tolerance) {
        WarningInFunction
          << mw.name() << " at time " << startName << " differs from the "
          << "element deposits of uniform/phaseChange/phaseChangeDeposits by "
          << "up to " << maxDifference << " kg/m2 (largest deposit "
          << maxDeposit << " kg/m2)." << nl << "The deposits of uniform/ are "
          << "used; " << mw.name() << " is output only and is rewritten from "
          << "them."
          << nl << endl;
      }
    }

    /*------------------------------------------------------------------------
    Cs_ likewise (plan section 20, item 5; review of M5, round 2: a Cs_
    edited by hand was overwritten silently): with the reservoirs of
    uniform/ it is compared in the sample cells of the pair, with the same
    tolerance (0 for a binary field without a repeated value, otherwise its
    rounding), and rewritten from them.
    ------------------------------------------------------------------------*/
    if (pair.inventory && reservoirsFromList
     && nReservoirFields == nReservoirPairs) {

      const volScalarField& Cs = pair.reservoirField();
      scalar maxDifference = 0, maxReservoir = 0;
      for (label k = 0; k < nSample; ++k) {
        if (samples_[elementSample_[k]].pair == pairi) {
          const scalar value = Cs[elementCell_[nWall_ + k]];
          maxDifference = max(maxDifference, mag(value - pair.reservoir[k]));
          maxReservoir = max(maxReservoir, mag(pair.reservoir[k]));
        }
      }
      reduce(maxDifference, maxOp<scalar>());
      reduce(maxReservoir, maxOp<scalar>());

      const bool fieldExact =
        returnReduceAnd(writtenBinary(IOobject(Cs.name(), startName, mesh_)))
     && !returnReduceOr(repeatedNonZero(Cs.primitiveField()));
      const scalar tolerance = fieldExact ? 0 : roundingBound*maxReservoir;

      if (maxDifference > tolerance) {
        WarningInFunction
          << Cs.name() << " at time " << startName << " differs from the "
          << "reservoirs of uniform/phaseChange/phaseChangeReservoirs by up to "
          << maxDifference << " kg/m3 (largest reservoir " << maxReservoir
          << " kg/m3)." << nl << "The reservoirs of uniform/ are used; "
          << Cs.name() << " is output only and is rewritten from them."
          << nl << endl;
      }
    }
  }

  /* output and uniform/ lists consistent with the state from the start */
  forAll(pairs_, pairi) {
    updateDepositField(pairs_[pairi]);
    updateReservoirField(pairs_[pairi]);
  }
  storeElementState();
}


/*------------------------------------------------------------------------------
Restart record of the interface (uniform/phaseChangeLayout): what the areas
A_e of the element deposits depend on
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::recordInterface() {

  phaseChangeLedger& ledger = *ledger_;
  dictionary& layout = ledger.layout();
  const word startName = mesh_.time().timeName();

  wordList patchNames(interfacePatches_.size());
  forAll(interfacePatches_, i) {
    patchNames[i] = mesh_.boundary()[interfacePatches_[i]].name();
  }
  const string multiplier(phaseChangeLedger::exactDecimal(areaMultiplier_));

  /* layer distance only: the distance, exact (not recorded for faceCells,
     whose record is unchanged) */
  const string distance (
    distanceLayer() ? phaseChangeLedger::exactDecimal(layerDistance_) : ""
  );

  /*--------------------------------------------------------------------------
  The deposits of a restart (uniform/phaseChange/phaseChangeDeposits, or
  mw_) are kept per effective element area A_e = |S_f| areaMultiplier
  (rescaled for areaPerVolume layer) on the faces of the interface patches.
  Another interface would re-scale them by the ratio of the areas, or drop
  the deposits of a patch, and the ledger would book the change as field
  rounding (ascii) or refuse the files as foreign (binary; review round
  4).  A restart must therefore keep the interface of the run that wrote
  the state; the record names what has changed.
  --------------------------------------------------------------------------*/
  if (ledger.restarted()) {

    const dictionary* old = layout.findDict("interface", keyType::LITERAL);

    if (!old) {
      WarningInFunction
        << "uniform/phaseChangeLayout at time " << startName << " does "
        << "not record the interface (written by an earlier build of the "
        << "solver): a change of the interface patches, areaPerVolume or "
        << "areaMultiplier since that run cannot be detected." << nl << endl;
    } else {

      OStringStream changes;
      const word oldLayer(old->get<word>("layer"));
      const wordList oldPatches(old->get<wordList>("patches"));
      const word oldAreaPerVolume(old->get<word>("areaPerVolume"));
      const string oldMultiplier(old->get<string>("areaMultiplier"));

      if (oldLayer != layer_) {
        changes << "    layer " << layer_ << " (written with " << oldLayer
          << ")" << nl;
      } else if (distanceLayer()) {
        const string oldDistance (
          old->getOrDefault<string>("distance", "unknown")
        );
        if (oldDistance != distance) {
          changes << "    distance " << distance.c_str() << " (written with "
            << oldDistance.c_str() << ")" << nl;
        }
      }
      if (oldPatches != patchNames) {
        changes << "    patches " << patchNames << " (written with "
          << oldPatches << ")" << nl;
      }
      if (oldAreaPerVolume != areaPerVolume_) {
        changes << "    areaPerVolume " << areaPerVolume_ << " (written with "
          << oldAreaPerVolume << ")" << nl;
      }
      if (oldMultiplier != multiplier) {
        changes << "    areaMultiplier " << multiplier.c_str()
          << " (written with " << oldMultiplier.c_str() << ")" << nl;
      }

      if (!changes.str().empty()) {
        FatalIOErrorInFunction(dict_.subDict("interface"))
          << "Restart from time " << startName << ": the interface differs "
          << "from the one of the run that wrote the state "
          << "(uniform/phaseChangeLayout):" << nl << changes.str().c_str()
          << "The element deposits of the restart are kept per effective "
          << "element area A_e = |S_f| areaMultiplier (rescaled for "
          << "areaPerVolume layer) on the faces of the interface patches; "
          << "with another interface the wall inventory would change by "
          << "the ratio of the areas, or lose the deposits of a patch."
          << nl;
        if (distanceLayer() || oldLayer != "faceCells") {
          FatalIOError << "With layer distance they are kept per cell of "
            << "the layer (A_e = V_c (A_I/V_I) areaMultiplier): another "
            << "layer rule or distance selects other cells and another "
            << "A_I/V_I." << nl;
        }
        FatalIOError
          << "Restart with the interface of that run, or start fresh."
          << exit(FatalIOError);
      }
    }
  }

  dictionary record;
  record.add("layer", layer_);
  if (distanceLayer()) {
    record.add("distance", distance);
  }
  record.add("patches", patchNames);
  record.add("areaPerVolume", areaPerVolume_);
  record.add("areaMultiplier", multiplier);
  layout.set("interface", record);
}


/*------------------------------------------------------------------------------
Samples: the release windows on a fresh start, the schedule and the
inventory samples on a restart, and their restart record (uniform/
phaseChangeLayout)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::recordSamples() {

  phaseChangeLedger& ledger = *ledger_;
  dictionary& layout = ledger.layout();
  const word startName = mesh_.time().timeName();
  const scalar tStart = ledger.time();
  const auto exact = &phaseChangeLedger::exactDecimal;

  /* the record of a sample: its mode and pair, the release or the
     inventory (scalars exact) and, for inventory, the removal */
  auto record = [&](const sampleData& s) -> dictionary {
    dictionary d;
    d.add("pair", string(speciesNames_[pairs_[s.pair].gas]));
    d.add("mode", string(s.inventory ? "inventory" : "release"));
    d.add("amount", string(exact(s.amount)));
    if (s.inventory) {
      d.add("cells", string(Foam::name(s.nCells)));
      d.add("volume", string(exact(s.recordVolume)));
      d.add("centre", string("(" + exact(s.centre.x()) + " "
                             + exact(s.centre.y()) + " "
                             + exact(s.centre.z()) + ")"));
      d.add("volumeTolerance", string(exact(s.volumeTolerance)));
      d.add("centreTolerance", string(exact(s.centreTolerance)));
      d.add("meshPoints", s.meshPoints);
      d.add("removeTime", string(s.removable ? exact(s.removeTime)
                                             : word("never")));
      d.add("removed", string(s.removed ? "yes" : "no"));
      if (s.removed) {
        d.add("removalClock", string(exact(s.removalClock)));
      }
    } else {
      d.add("startTime", string(exact(s.startTime)));
      d.add("duration", string(exact(s.duration)));
    }
    return d;
  };

  /* the window [start, start + duration] of a sample (record or current) */
  auto window = [](const scalar start, const scalar duration) -> string {
    OStringStream os;
    os << "[" << start << ", " << start + duration << "] s";
    return os.str();
  };

  /* the mode of a record; a record of an earlier build is a release */
  auto modeOf = [](const dictionary& d) -> string {
    return d.getOrDefault<string>("mode", "release");
  };

  const dictionary* old = ledger.restarted()
    ? layout.findDict("samples", keyType::LITERAL)
    : nullptr;

  /*--------------------------------------------------------------------------
  The error bounds of V_Z and the centre of every inventory sample that the
  rounding of this run's mesh points implies (sampleGeometryTolerance()),
  recorded with the sample, so that a restart can add those of the run
  that wrote the state to its own (review of M5, round 3: with the bound
  of the reading run's points alone, a restart on the binary mesh after
  reconstructPar of an ascii-decomposed run was refused, and so was an
  ascii-decomposed sample far from the origin).  The format of the points
  (constant/polyMesh/points of the case or of every processor) is read
  collectively, once.
  --------------------------------------------------------------------------*/
  bool anyInventory = false;
  for (const sampleData& s : samples_) {
    anyInventory = anyInventory || s.inventory;
  }
  const bool exactPoints = !anyInventory || returnReduceAnd (
    writtenBinary (
      IOobject (
        "points",
        mesh_.pointsInstance(),
        polyMesh::meshSubDir,
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
      ),
      pointIOField::typeName
    )
  );
  for (sampleData& s : samples_) {
    if (s.inventory) {
      sampleGeometryTolerance(s, exactPoints);
    }
  }

  if (!ledger.restarted()) {

    /*------------------------------------------------------------------------
    A fresh start at t_start releases only the part of a window after
    t_start.  The part before it would never be released (plan E7: the
    release sums to n0 M), so a window that began before the start stops
    the run.
    ------------------------------------------------------------------------*/
    for (const sampleData& s : samples_) {
      if (s.inventory) {

        /* filled and removed at the start of the first step, before any
           evaporation (removeSamples(): t0 >= removeTime - deltaT/2) */
        const scalar dt = mesh_.time().deltaTValue();
        if (s.removable && tStart >= s.removeTime - 0.5*dt) {
          FatalIOErrorInFunction(dict_.subDict("samples").subDict(s.name))
            << "Sample " << s.name << " (mode inventory): the run starts "
            << "fresh at t = " << tStart << " s, at or after removeTime - "
            << "deltaT/2 = " << s.removeTime - 0.5*dt << " s.  Its "
            << "reservoir (" << s.amount << " mol) would be filled and "
            << "removed at the start of the first step, without any "
            << "evaporation." << nl << "Set removeTime later than "
            << tStart + 0.5*dt << " s, or remove the entry (never "
            << "removed)." << exit(FatalIOError);
        }
        continue;
      }
      const scalar missed =
        releaseOverlap(s.startTime, tStart, s.startTime, s.duration);
      if (missed > 0) {
        FatalIOErrorInFunction(dict_.subDict("samples").subDict(s.name))
          << "Sample " << s.name << ": the run starts fresh at t = " << tStart
          << " s, after the start " << s.startTime << " s of its release "
          << "window " << window(s.startTime, s.duration).c_str() << "." << nl
          << missed/s.length*100 << " % of its amount ("
          << s.amount*missed/s.length << " mol) would never be released."
          << nl << "Set startTime to " << tStart << " or later, or start "
          << "the run at t <= " << s.startTime << "." << exit(FatalIOError);
      }
    }

  } else if (!old) {

    /*------------------------------------------------------------------------
    A layout of an earlier build: its samples were release samples (mode
    inventory did not exist), so an inventory sample now is new.
    ------------------------------------------------------------------------*/
    for (const sampleData& s : samples_) {
      if (s.inventory) {
        FatalIOErrorInFunction(dict_.subDict("samples").subDict(s.name))
          << "Sample " << s.name << " (mode inventory) is new at the restart "
          << "from time " << startName << ": uniform/phaseChangeLayout "
          << "(an earlier build) records no inventory sample." << nl
          << "The reservoir of an inventory sample is filled on a fresh "
          << "start only; start fresh to add it." << exit(FatalIOError);
      }
    }

    WarningInFunction
      << "uniform/phaseChangeLayout at time " << startName << " does not "
      << "record the release samples (written by an earlier build of the "
      << "solver): a change of their schedule since that run cannot be "
      << "detected." << nl << endl;

  } else {

    /*------------------------------------------------------------------------
    The removals recorded per sample must add up to the count of the ledger
    (REMOVAL_DONE of every pair): both are written together, so a
    difference means a uniform/ assembled from different runs.
    ------------------------------------------------------------------------*/
    forAll(pairs_, pairi) {
      label nRemoved = 0;
      for (const entry& e : *old) {
        if (e.isDict() && modeOf(e.dict()) == "inventory"
          && e.dict().get<string>("pair") == speciesNames_[pairs_[pairi].gas]
          && e.dict().get<string>("removed") == "yes") {
          ++nRemoved;
        }
      }
      const scalar counted = ledger(pairi, phaseChangeLedger::REMOVAL_DONE);
      if (scalar(nRemoved) != counted) {
        FatalErrorInFunction
          << "Restart from time " << startName << ": the ledger counts "
          << counted << " sample removals of pair "
          << speciesNames_[pairs_[pairi].gas] << ", uniform/"
          << "phaseChangeLayout records " << nRemoved << "." << nl
          << "The files of uniform/ do not belong to the same run and time."
          << exit(FatalError);
      }
    }

    /*------------------------------------------------------------------------
    Restart at the ledger's clock tStart.  The ledger has booked the release
    up to tStart under the schedule recorded by the run that wrote the
    state, so a release sample whose window has begun must keep its pair,
    amount, startTime and duration; once its window has ended it may be
    removed.  A release sample added now must not have begun (as on a fresh
    start).  Release samples whose windows lie after tStart may be added,
    changed or removed (review round 4: such changes were accepted
    silently, and the total release was no longer n0 M).

    An inventory sample holds a reservoir in its cells: it must exist in
    the record (its fill belongs to a fresh start), keep its pair, amount,
    cells and volume, and may be dropped only after its removal; its
    removeTime and kinetics may change.  No sample changes its mode.
    ------------------------------------------------------------------------*/

    const boundBox& bb = mesh_.bounds();
    const scalar extent = mag(bb.span()) + max(mag(bb.min()), mag(bb.max()));

    for (sampleData& s : samples_) {

      const dictionary& sd = dict_.subDict("samples").subDict(s.name);
      const dictionary now(record(s));
      const dictionary* was = old->findDict(s.name, keyType::LITERAL);

      if (!was) {
        if (s.inventory) {
          FatalIOErrorInFunction(sd)
            << "Sample " << s.name << " (mode inventory) is new at the "
            << "restart from time " << startName << " (t = " << tStart
            << " s; not in uniform/phaseChangeLayout)." << nl
            << "The reservoir of an inventory sample is filled on a fresh "
            << "start only; start fresh to add it." << exit(FatalIOError);
        }
        const scalar begun =
          releaseOverlap(s.startTime, tStart, s.startTime, s.duration);
        if (begun > 0) {
          FatalIOErrorInFunction(sd)
            << "Sample " << s.name << " is new at the restart from time "
            << startName << " (t = " << tStart << " s; not in uniform/"
            << "phaseChangeLayout), but its release window "
            << window(s.startTime, s.duration).c_str() << " began before."
            << nl << begun/s.length*100 << " % of its amount ("
            << s.amount*begun/s.length << " mol) would never be released."
            << nl << "Set startTime to " << tStart << " or later."
            << exit(FatalIOError);
        }
        continue;
      }

      const string wasMode(modeOf(*was));
      if (wasMode != now.get<string>("mode")) {
        FatalIOErrorInFunction(sd)
          << "Sample " << s.name << " changes its mode at the restart from "
          << "time " << startName << ": mode " << now.get<string>("mode")
          << " (written with " << wasMode << ")." << nl
          << "Give the new sample another name." << exit(FatalIOError);
      }

      if (s.inventory) {

        /*--------------------------------------------------------------------
        pair, amount and the number of cells are compared exactly.  V_Z and
        the centre are sums over the cells of all ranks, whose last bits
        depend on the decomposition (decomposePar, reconstructPar,
        renumberMesh), 1e-12 (relative, and of the mesh extent), and on the
        rounding of the mesh points of both runs: the bounds recorded by
        the run that wrote the state (0 in records of an earlier build) plus
        those of this run (sampleGeometryTolerance()), but at most the
        limits of this run's sample (what a selection changed by one cell
        can cause; review of M5, round 4: bounds from the digits of round
        coordinates accepted a sample moved by 20 mm).  The centre names a
        selection moved to other cells of the same number and volume
        (records of an earlier build have no centre).
        --------------------------------------------------------------------*/
        OStringStream changes;
        for (const word key : {"pair", "amount", "cells"}) {
          const string before(was->get<string>(key));
          const string after(now.get<string>(key));
          if (before != after) {
            changes << "    " << key << " " << after.c_str()
              << " (written with " << before.c_str() << ")" << nl;
          }
        }

        const string wasPoints =
          was->getOrDefault<string>("meshPoints", "not recorded");
        bool limited = false;
        auto tolerance = [&](const word& key, const scalar here,
                             const scalar summation, const scalar limit,
                             const char* why, const char* unit) {
          const scalar written =
            readScalar(was->getOrDefault<string>(key, "0"));
          const scalar rounding = min(written + here, limit);
          OStringStream note;
          note << "allowed " << summation + rounding << unit
            << ": 1e-12 of the summation order, " << written << unit
            << " for the mesh points of the run that wrote the state ("
            << wasPoints.c_str() << "), " << here << unit << " for those "
            << "of this run (" << s.meshPoints.c_str() << ")";
          if (written + here > limit) {
            note << ", limited to " << limit << unit << " (" << why << ")";
          }
          /* the allowed difference, the note, and the difference that the
             rounding alone would allow */
          return std::make_tuple (
            summation + rounding, note.str(), summation + written + here
          );
        };

        const scalar wasVolume = readScalar(was->get<string>("volume"));
        const auto volumeBound = tolerance (
          "volumeTolerance", s.volumeTolerance, 1e-12*mag(wasVolume),
          s.volumeLimit, "half the smallest cell volume of the sample", " m3"
        );
        const scalar volumeDifference = mag(s.recordVolume - wasVolume);
        if (volumeDifference > std::get<0>(volumeBound)) {
          changes << "    volume " << now.get<string>("volume").c_str()
            << " m3 (written with " << was->get<string>("volume").c_str()
            << "; difference " << s.recordVolume - wasVolume << " m3, "
            << std::get<1>(volumeBound).c_str() << ")" << nl;
          limited = limited || volumeDifference <= std::get<2>(volumeBound);
        }
        if (was->found("centre", keyType::LITERAL)) {
          IStringStream is(was->get<string>("centre"));
          vector wasCentre(Zero);
          is >> wasCentre;
          const auto centreBound = tolerance (
            "centreTolerance", s.centreTolerance, 1e-12*extent,
            s.centreLimit, "a quarter of the centre shift of a selection "
            "with one cell replaced", " m"
          );
          const scalar distance = mag(s.centre - wasCentre);
          if (distance > std::get<0>(centreBound)) {
            changes << "    centre " << now.get<string>("centre").c_str()
              << " m (written with " << was->get<string>("centre").c_str()
              << "; distance " << distance << " m, "
              << std::get<1>(centreBound).c_str() << "; mesh extent "
              << extent << " m)" << nl;
            limited = limited || distance <= std::get<2>(centreBound);
          }
        }
        if (!changes.str().empty()) {
          FatalIOErrorInFunction(sd)
            << "Inventory sample " << s.name << " was changed at the "
            << "restart from time " << startName << ":" << nl
            << changes.str().c_str()
            << "Its reservoir (uniform/phaseChange/phaseChangeReservoirs) "
            << "belongs to the cells and the amount of the run that wrote the "
            << "state."
            << nl << "Restore the sample of uniform/phaseChangeLayout, or "
            << "start fresh.  A selection whose boundary passes through "
            << "cell centres (e.g. a cylinder radius on a ring of centres) "
            << "selects other cells when the mesh points are rewritten "
            << "(ascii, another precision); keep the boundary between the "
            << "centres."
            << (
                 limited
               ? "  The difference lies within the rounding of the ascii mesh "
                 "points at the digits above, which cannot be told from a "
                 "changed selection; write the mesh points in binary, or with "
                 "more digits (writePrecision)."
               : ""
               )
            << exit(FatalIOError);
        }

        s.removed = was->get<string>("removed") == "yes";
        if (s.removed) {
          s.removalClock = readScalar(was->get<string>("removalClock"));
        }
        continue;
      }

      OStringStream changes;
      for (const word key : {"pair", "amount", "startTime", "duration"}) {
        const string before(was->get<string>(key));
        const string after(now.get<string>(key));
        if (before != after) {
          changes << "    " << key << " " << after.c_str()
            << " (written with " << before.c_str() << ")" << nl;
        }
      }
      if (changes.str().empty()) {
        continue;
      }

      const scalar begun =
        releaseOverlap(s.startTime, tStart, s.startTime, s.duration);
      const scalar wasStart = readScalar(was->get<string>("startTime"));
      const scalar wasDuration = readScalar(was->get<string>("duration"));
      const scalar wasBegun =
        releaseOverlap(wasStart, tStart, wasStart, wasDuration);

      if (begun > 0 || wasBegun > 0) {
        FatalIOErrorInFunction(sd)
          << "Sample " << s.name << " was changed at the restart from time "
          << startName << " (t = " << tStart << " s), after the start of "
          << "its release window (recorded "
          << window(wasStart, wasDuration).c_str() << ", now "
          << window(s.startTime, s.duration).c_str() << "):" << nl
          << changes.str().c_str()
          << "The ledger has booked its release up to the restart; the "
          << "release of a sample must not change once its window has "
          << "begun, or its total is no longer the amount n0." << nl
          << "Restore the sample of uniform/phaseChangeLayout, or give a "
          << "new release another sample name and a window after t = "
          << tStart << " s." << exit(FatalIOError);
      }
    }

    for (const entry& e : *old) {

      if (!e.isDict()) {
        continue;
      }
      bool present = false;
      for (const sampleData& s : samples_) {
        present = present || s.name == e.keyword();
      }
      if (present) {
        continue;
      }

      if (modeOf(e.dict()) == "inventory") {
        if (e.dict().get<string>("removed") != "yes") {
          FatalIOErrorInFunction(dict_)
            << "Inventory sample " << e.keyword() << " of the run that "
            << "wrote the state at time " << startName << " (uniform/"
            << "phaseChangeLayout) is missing, but it has not been removed "
            << "(removeTime " << e.dict().get<string>("removeTime").c_str()
            << "): its reservoir would vanish from the balance." << nl
            << "Keep the sample until its removal (removeTime)."
            << exit(FatalIOError);
        }
        continue;
      }

      const scalar start = readScalar(e.dict().get<string>("startTime"));
      const scalar duration = readScalar(e.dict().get<string>("duration"));
      const scalar amount = readScalar(e.dict().get<string>("amount"));
      const scalar length = releaseWindowLength(start, duration);
      const scalar released = releaseOverlap(start, tStart, start, duration);
      const scalar remaining = releaseOverlap(tStart, VGREAT, start, duration);

      if (released > 0 && remaining > 0) {
        FatalIOErrorInFunction(dict_)
          << "Sample " << e.keyword() << " of the run that wrote the state "
          << "at time " << startName << " (uniform/phaseChangeLayout) is "
          << "missing, but its release window "
          << window(start, duration).c_str() << " has not ended at the "
          << "restart (t = " << tStart << " s)." << nl
          << remaining/length*100 << " % of its amount ("
          << amount*remaining/length << " mol) would never be released."
          << nl << "Keep the sample until its window has ended."
          << exit(FatalIOError);
      }
    }
  }

  dictionary samples;
  for (const sampleData& s : samples_) {
    samples.add(s.name, record(s));
    if (s.inventory && s.removed) {
      Info<< "Sample " << s.name << " was removed at t = " << s.removalClock
        << " s by the run that wrote the state; it stays empty" << nl;
    }
  }

  /*--------------------------------------------------------------------------
  An inventory sample dropped after its removal keeps its record for the
  rest of the run: the ledger's REMOVAL_DONE still counts it, and every
  later restart checks the count against the records (review of M5: the
  record was lost, and the next restart stopped).  Re-adding it under the
  same name finds it removed.
  --------------------------------------------------------------------------*/
  if (old) {
    for (const entry& e : *old) {
      if (e.isDict() && !samples.found(e.keyword(), keyType::LITERAL)
        && modeOf(e.dict()) == "inventory"
        && e.dict().get<string>("removed") == "yes") {
        samples.add(e.keyword(), e.dict());
        Info<< "Sample " << e.keyword() << " (mode inventory, removed at t = "
          << e.dict().get<string>("removalClock").c_str() << " s) is no "
          << "longer in the dictionary; its record is kept in uniform/"
          << "phaseChangeLayout" << nl;
      }
    }
  }

  layout.set("samples", samples);
}


/*------------------------------------------------------------------------------
Error bounds of the recorded V_Z and centre of an inventory sample from the
rounding of the mesh points (interfaceExchange.H)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::sampleGeometryTolerance (
  sampleData& s,
  const bool  exactPoints
) const {

  const pointField& points = mesh_.points();
  const faceList& faces = mesh_.faces();
  const cellList& cells = mesh_.cells();
  const vectorField& faceAreas = mesh_.faceAreas();
  const scalarField& V = mesh_.V();

  /* the face areas of the sample cells (a bound of its surface), the
     largest coordinate and distance from the centre of their points, the
     digits their coordinates need (ascii points only), and the thinnest
     cell (its volume over its largest face) and the smallest volume */
  scalar surface = 0, coordinate = 0, radius = 0;
  scalar thinnest = GREAT, smallest = GREAT;
  label digits = 0;
  for (const label celli : s.cells) {
    scalar largestFace = 0;
    for (const label facei : cells[celli]) {
      surface += mag(faceAreas[facei]);
      largestFace = max(largestFace, mag(faceAreas[facei]));
      for (const label pointi : faces[facei]) {
        const point& p = points[pointi];
        radius = max(radius, mag(p - s.centre));
        for (direction d = 0; d < vector::nComponents; ++d) {
          coordinate = max(coordinate, mag(p[d]));
          if (!exactPoints && p[d] != 0) {
            digits = max(digits, significantDigits(p[d]));
          }
        }
      }
    }
    thinnest = min(thinnest, V[celli]/max(largestFace, VSMALL));
    smallest = min(smallest, V[celli]);
  }
  reduce(surface, sumOp<scalar>());
  reduce(coordinate, maxOp<scalar>());
  reduce(radius, maxOp<scalar>());
  reduce(digits, maxOp<label>());
  reduce(thinnest, minOp<scalar>());
  reduce(smallest, minOp<scalar>());

  s.centreLimit =
    s.nCells > 0 ? 0.25*thinnest*smallest/max(s.recordVolume, VSMALL) : 0;
  s.volumeLimit = s.nCells > 0 ? 0.5*smallest : 0;

  if (exactPoints || coordinate == 0) {
    s.volumeTolerance = 0;
    s.centreTolerance = 0;
    s.meshPoints = "binary";
    return;
  }

  /* the decimal exponent e of the largest coordinate, 10^e <= c < 10^(e+1)
     (log10 may round at a power of ten) */
  label exponent = label(std::floor(std::log10(coordinate)));
  if (std::pow(10.0, exponent + 1) <= coordinate) {
    ++exponent;
  } else if (std::pow(10.0, exponent) > coordinate) {
    --exponent;
  }
  const label fewest =
    min(label(6), label(max(IOstream::defaultPrecision(), 1u)));
  digits = max(digits, fewest);
  const scalar rounding = 0.5*std::pow(10.0, exponent + 1 - digits);
  const scalar displacement = 2*std::sqrt(3.0)*rounding;

  s.volumeTolerance = displacement*surface;
  s.centreTolerance =
    displacement*surface*radius/max(s.recordVolume, VSMALL);
  s.meshPoints = "ascii, " + Foam::name(digits) + " significant digits";
}


/*------------------------------------------------------------------------------
Inventories [kg]: compensated sums, independent of the summation order
------------------------------------------------------------------------------*/

Foam::scalar LESTO::interfaceExchange::gasInventory (
  const pairData& pair
) const {
  const scalarField& rhoI = rho_.primitiveField();
  const scalarField& Y = species_[pair.gas].primitiveField();
  const scalarField& V = mesh_.V();

  compensatedSum sum;
  forAll(Y, celli) {
    sum.add(rhoI[celli]*Y[celli]*V[celli]);
  }
  return globalSum(sum);
}


Foam::scalar LESTO::interfaceExchange::wallInventory (
  const pairData& pair
) const {
  compensatedSum sum;
  forAll(pair.deposit, e) {
    sum.add(pair.deposit[e]*elementArea_[e]);
  }
  return globalSum(sum);
}


Foam::scalar LESTO::interfaceExchange::reservoirMass (
  const pairData& pair,
  const label     e
) const {
  return pair.reservoir[e - nWall_]*mesh_.V()[elementCell_[e]];
}


Foam::scalar LESTO::interfaceExchange::sampleInventory (
  const pairData& pair
) const {
  compensatedSum sum;
  for (label e = nWall_; e < elementCell_.size(); ++e) {
    sum.add(reservoirMass(pair, e));
  }
  return globalSum(sum);
}


/*------------------------------------------------------------------------------
Round-off threshold of the mass of a pair [kg]
------------------------------------------------------------------------------*/

Foam::scalar LESTO::interfaceExchange::massTolerance(const label pairi) const {
  const pairData& pair = pairs_[pairi];
  const scalar inventory = mag((*ledger_)(pairi, phaseChangeLedger::HELD))
                         + pair.releasedStep + pair.inflowStep;
  return guardTolerance_*max(inventory, VSMALL);
}


/*------------------------------------------------------------------------------
Fingerprints of the elements of this rank: their number and a 64-bit hash
(two labels of 31 bits).  A wall element contributes its cell, patch, patch
face and the exact face centre; a sample element its cell, the name of its
sample and the exact cell centre.  A restart onto another decomposition or a
renumbered mesh changes them, even when the numbers of elements are the same,
and so does a changed set of inventory samples.  The name, not the position
in the samples dictionary, identifies the sample: release samples added in
front of an inventory sample leave its elements unchanged (review of M5).
------------------------------------------------------------------------------*/

Foam::labelList LESTO::interfaceExchange::elementsFingerprint (
  const int which
) const {

  fnv1a hash;
  const surfaceVectorField::Boundary& Cfb = mesh_.Cf().boundaryField();
  const volVectorField& C = mesh_.C();

  const label first = which == 2 ? nWall_ : 0;
  const label last = which == 0 ? nWall_ : elementCell_.size();

  /*--------------------------------------------------------------------------
  layer distance: the rule and the distance first (plan section 16, item 1:
  a restart with another rule or distance does not fit), then per element
  its cell and the exact cell centre.  layer faceCells hashes exactly as
  before (its restart format is unchanged).
  --------------------------------------------------------------------------*/
  if (distanceLayer() && first < nWall_) {
    const char rule[] = "distance";
    hash.add(rule, sizeof(rule));
    hash.add(&layerDistance_, sizeof(layerDistance_));
  }

  for (label e = first; e < last; ++e) {
    if (e < nWall_ && distanceLayer()) {
      const vector& Cc = C[elementCell_[e]];
      const label ids[1] = {elementCell_[e]};
      const scalar centre[3] = {Cc.x(), Cc.y(), Cc.z()};
      hash.add(ids, sizeof(ids));
      hash.add(centre, sizeof(centre));
    } else if (e < nWall_) {
      const label patchi = elementPatch_[e];
      const label facei = elementFace_[e];
      const vector& Cf = Cfb[patchi][facei];
      const label ids[3] = {elementCell_[e], patchi, facei};
      const scalar centre[3] = {Cf.x(), Cf.y(), Cf.z()};
      hash.add(ids, sizeof(ids));
      hash.add(centre, sizeof(centre));
    } else {
      const vector& Cc = C[elementCell_[e]];
      const word& sample = samples_[elementSample_[e - nWall_]].name;
      const label ids[2] = {elementCell_[e], label(sample.size())};
      const scalar centre[3] = {Cc.x(), Cc.y(), Cc.z()};
      hash.add(ids, sizeof(ids));
      hash.add(sample.data(), sample.size());
      hash.add(centre, sizeof(centre));
    }
  }

  const std::uint64_t h = hash.value();
  return labelList ({
    last - first,
    label(h & 0x7fffffffULL),
    label((h >> 31) & 0x7fffffffULL)
  });
}


/*------------------------------------------------------------------------------
Solver dictionary of the next solve
------------------------------------------------------------------------------*/

const Foam::dictionary& LESTO::interfaceExchange::solverDict (
  const label pairi,
  const bool  final
) const {
  const pairData& pair = pairs_[pairi];
  return mesh_.solverDict (
    final && pair.finalSolver ? word(pair.solverName + "Final")
                              : pair.solverName
  );
}


/*------------------------------------------------------------------------------
Start of a physical time step
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::beginStep() {

  if (mock()) {
    return;
  }

  const Time& runTime = mesh_.time();
  const scalar dt = runTime.deltaTValue();

  /*--------------------------------------------------------------------------
  The exact step [t0, t1] on the clock of the ledger (phaseChangeLedger.H):
  t1 = t0 + deltaT, the arithmetic of Time::operator++, restored exactly on
  a restart.  The time-directory names are rounded to timePrecision (and a
  restarted OpenFOAM continues from the rounded name), so they would give
  a wrong overlap of the step with a release window.
  --------------------------------------------------------------------------*/
  stepStart_ = ledger_->time();
  stepEnd_ = stepStart_ + dt;
  stepDeltaT_ = dt;
  const scalar t0 = stepStart_;
  const scalar t1 = stepEnd_;

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];

    /* fvSolution and fvSchemes may have been modified (runTimeModifiable) */
    checkFinalSolver(pair, false);
    checkDdtScheme(pair, false);
    checkConvectionScheme(pair, false);

    pair.Srel.field() = 0;
    pair.releasedStep = 0;
    pair.defect = 0;
    pair.residualOF = 0;
    pair.explicitRemainder = 0;
    pair.nExplicitRemainder = 0;
    pair.explicitSource = 0;
    pair.transportIdentity = 0;
    pair.identityTerms = 0;
    pair.identityRelative = 0;
    pair.transport = 0;
    pair.clampedStep = 0;
    pair.projectedStep = 0;
    pair.evaporatedAboveWarningStep = 0;
    pair.nGuardPasses = 0;
    pair.nGuardSwitched = 0;
  }

  /*--------------------------------------------------------------------------
  Removal of the inventory samples due at the start of this step (plan E7),
  before the bounds: a removed sample is empty (NONE) from this step on.
  --------------------------------------------------------------------------*/
  removeSamples(t0, dt);

  /*--------------------------------------------------------------------------
  Bounds of the step (plan E3).  Wall elements: lower = -infinity; the
  temperature model is irreversible (upper = 0, whatever the deposit), HKS
  may evaporate the deposit of the element, upper = m''_e A_e/dt (0 for a
  bare wall: deposit only).  Sample elements: no deposition (lower = 0),
  evaporation up to the reservoir, upper = Cs V/dt.
  --------------------------------------------------------------------------*/
  const scalarField& V = mesh_.V();

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];

    for (label e = 0; e < nWall_; ++e) {
      pair.lower[e] = -VGREAT;
      pair.upper[e] = model_ == modelType::HKS
        ? pair.deposit[e]*elementArea_[e]/dt
        : 0;
    }
    for (label e = nWall_; e < elementCell_.size(); ++e) {
      pair.lower[e] = 0;
      pair.upper[e] = reservoirMass(pair, e)/dt;
    }
  }

  /*--------------------------------------------------------------------------
  equilibrium GEMS, mode local (M9): p_eq (and chi) of every element from
  GEMS3K at the composition of the start of the step, at the steps of the
  schedule of the absolute time index n ((n - 1) mod updateInterval = 0;
  also at the first step of a run without a lagged state;
  gemsEquilibrium::updateDue()), and a_e, g_e with them; the step
  then uses one law throughout.  Collective (the counts).
  --------------------------------------------------------------------------*/
  if (gems_ && gems_->local() && gems_->updateDue(runTime.timeIndex())) {
    updateEquilibrium();
  }

  /*--------------------------------------------------------------------------
  Release: Srel = n0 M |[t0,t1] ^ [start, start+duration]|/(length dt V_Z)
  in every cell of the sample.  length is the window as the overlaps
  represent it (releaseWindowLength), so that the overlaps of all steps sum
  to it and the release to n0 M, also for a window far from t = 0.
  --------------------------------------------------------------------------*/
  for (const sampleData& s : samples_) {
    if (s.inventory) {
      continue;
    }
    pairData& pair = pairs_[s.pair];
    const scalar overlap = releaseOverlap(t0, t1, s.startTime, s.duration);
    const scalar rate =
      s.amount*pair.molarMass*overlap/(s.length*dt*s.volume);

    for (const label celli : s.cells) {
      pair.Srel[celli] = rate;
      pair.releasedStep += rate*V[celli]*dt;
    }
  }

  /*--------------------------------------------------------------------------
  Convective inflow of the step through the non-coupled patches, from the
  boundary values at its start: only a scale for massTolerance().
  --------------------------------------------------------------------------*/
  forAll(pairs_, pairi) {
    pairData& pair = pairs_[pairi];
    const volScalarField& Y = species_[pair.gas];

    scalar inflow = 0;
    for (const label patchi : ledger_->patchIds()) {
      const scalarField& phip = phi_.boundaryField()[patchi];
      const scalarField& Yp = Y.boundaryField()[patchi];
      forAll(phip, facei) {
        inflow += max(-phip[facei], scalar(0))*mag(Yp[facei]);
      }
    }
    pair.inflowStep = returnReduce(inflow, sumOp<scalar>())*dt;

    reduce(pair.releasedStep, sumOp<scalar>());
    pair.Su.field() = pair.Srel.field();
    pair.Sp.field() = 0;
  }
}


/*------------------------------------------------------------------------------
Removal of inventory samples (plan E7): once, at the start of the first step
whose start t0 on the ledger's clock satisfies t0 >= removeTime - dt/2
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::removeSamples (
  const scalar t0,
  const scalar dt
) {

  phaseChangeLedger& ledger = *ledger_;

  for (sampleData& s : samples_) {

    if (!s.inventory || !s.removable || s.removed
      || !(t0 >= s.removeTime - 0.5*dt)) {
      continue;
    }

    /* the reservoir of the sample, compensated over all ranks; Cs := 0 */
    pairData& pair = pairs_[s.pair];
    compensatedSum sum;
    forAll(s.cells, i) {
      const label e = s.elements[i];
      sum.add(reservoirMass(pair, e));
      pair.reservoir[e - nWall_] = 0;
    }
    const scalar removed = globalSum(sum);

    ledger.add(s.pair, phaseChangeLedger::REMOVED, removed);
    ledger.add(s.pair, phaseChangeLedger::REMOVAL_DONE, 1);
    s.removed = true;
    s.removalClock = t0;

    /* the restart record: this sample is removed (uniform/
       phaseChangeLayout, written with the ledger) */
    dictionary& record = ledger.layout().subDict("samples").subDict(s.name);
    record.set("removed", string("yes"));
    record.set("removalClock", string(phaseChangeLedger::exactDecimal(t0)));

    Info<< "Sample " << s.name << " removed at the start of the step at t = "
      << t0 << " s (removeTime " << s.removeTime << " s): its reservoir of "
      << removed << " kg (" << removed/pair.molarMass << " mol) of "
      << speciesNames_[pair.condensed] << " is booked as removed" << nl;
  }
}


/*------------------------------------------------------------------------------
Su, Sp of exchange cell k from the regimes of its elements
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::cellCoefficients (
  pairData&   pair,
  const label k
) const {

  const label celli = exchangeCells_[k];
  scalar su = 0, sp = 0;

  for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
    if (pair.regime[e] == IMPLICIT) {
      su += pair.a[e];
      sp += pair.g[e];
    } else if (pair.regime[e] == CAPPED) {
      su += pair.upper[e];
    }
  }

  const scalar V = mesh_.V()[celli];
  pair.Su[celli] = pair.Srel[celli] + su/V;
  pair.Sp[celli] = sp/V;
}


/*------------------------------------------------------------------------------
Local implicit predictor (E4) at an outer corrector
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::linearise (
  const label           pairi,
  const fvScalarMatrix& transportEqn
) {

  pairData& pair = pairs_[pairi];

  const tmp<volScalarField> tA(transportEqn.A());
  const tmp<volScalarField> tH(transportEqn.H());
  const scalarField& A = tA().primitiveField();
  const scalarField& H = tH().primitiveField();
  const scalarField& V = mesh_.V();
  const scalar dt = mesh_.time().deltaTValue();
  const scalar significant = massTolerance(pairi);

  label nChanged = 0, nSignificant = 0;
  scalar changedFlow = 0;
  List<int> cellRegime, cellPrevious;

  forAll(exchangeCells_, k) {

    const label celli = exchangeCells_[k];
    const label first = cellStart_[k];
    const int   n = int(cellStart_[k+1] - first);

    /*------------------------------------------------------------------------
    The regimes of the previous corrector (of the previous step at the
    first corrector): an element does not jump between its two bounds
    (cellPredictor, acrossBounds).  UNKNOWN at a fresh start.
    ------------------------------------------------------------------------*/
    cellRegime.resize(n);
    cellPrevious.resize(n);
    for (int i = 0; i < n; ++i) {
      cellPrevious[i] = int(pair.regime[first + i]);
    }

    const scalar Yloc = cellPredictor (
      A[celli]*V[celli],
      (H[celli] + pair.Srel[celli])*V[celli],
      n,
      pair.a.cdata() + first,
      pair.g.cdata() + first,
      pair.lower.cdata() + first,
      pair.upper.cdata() + first,
      cellRegime.data(),
      cellPrevious.cdata()
    );

    /*------------------------------------------------------------------------
    A change of regime is counted, and the change of the element's flow at
    the predictor value Y_loc is summed for the log (e.g. sign changes of Y
    at round-off level ahead of a gas front carry no mass).  A change whose
    own mass exceeds the round-off threshold is significant.
    ------------------------------------------------------------------------*/
    for (int i = 0; i < n; ++i) {
      const label e = first + i;
      if (pair.regime[e] != cellRegime[i]) {
        const scalar change = mag (
          elementFlow(cellRegime[i], pair.a[e], pair.g[e], pair.upper[e], Yloc)
        - elementFlow(pair.regime[e], pair.a[e], pair.g[e], pair.upper[e], Yloc)
        );
        ++nChanged;
        changedFlow += change;
        if (change*dt > significant) {
          ++nSignificant;
        }
        pair.regime[e] = cellRegime[i];
      }
    }

    cellCoefficients(pair, k);
  }

  pair.nChanged = returnReduce(nChanged, sumOp<label>());
  pair.nChangedSignificant = returnReduce(nSignificant, sumOp<label>());
  pair.changedMass = returnReduce(changedFlow, sumOp<scalar>())*dt;
}


/*------------------------------------------------------------------------------
Decision of the outer loop: residual, and regime changes that carry no mass
do not keep the loop going
------------------------------------------------------------------------------*/

bool LESTO::interfaceExchange::outerConverged (
  const label  pairi,
  const scalar initialResidual
) {
  pairData& pair = pairs_[pairi];
  pair.nChangedOuter = pair.nChanged;
  pair.nChangedSignificantOuter = pair.nChangedSignificant;
  pair.changedMassOuter = pair.changedMass;
  return initialResidual < tolerance_
      && pair.changedMass <= massTolerance(pairi);
}


/*------------------------------------------------------------------------------
Transport per non-coupled patch and solver defect of the final matrix
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::recordSolve (
  const label           pairi,
  const fvScalarMatrix& eqn,
  const fvScalarMatrix& transportEqn,
  const volScalarField& Y
) {

  pairData& pair = pairs_[pairi];
  const labelList& patchIds = ledger_->patchIds();
  const scalarField& V = mesh_.V();
  const scalarField& YI = Y.primitiveField();

  /*--------------------------------------------------------------------------
  R' per cell: the defect of the transport-only matrix, plus the realised
  flow of every element of the cell at the solved Y, computed exactly as
  realise() computes the flow that updates the condensate, plus the release
  Srel V.  In exact arithmetic this is the defect b - A Y of the assembled
  matrix; in double precision it avoids the rounding of the stiff assembled
  diagonal and source (about eps lambda deltaT of the gas of an exchange
  cell per step), which differs from that of q_e = a_e - g_e Y and made the
  closure drift linearly near equilibrium (matrixDefect.H).

  The transport-only defect is the flux form E - (m - m0)/deltaT - sum_f
  (F_f + C_f) (matrixDefect.H): s the source of the ddt term (taken from
  fvm::ddt itself, the same doubles as in the solved matrix), C the
  face-flux correction of the matrix (the non-orthogonal correction of the
  laplacian: the paired species is fluxRequired) entered by its face
  values, E = b - fl(s + c) what is left of the source besides the storage
  and the correction's contribution c (0 for the solver's schemes), m0 and
  m the masses of the cells as gasInventory() computes them.  Its face
  fluxes and its storage telescope exactly, so the double-precision
  imbalances of the discrete operator (the column sums of the assembled
  ddt + div - laplacian are not exactly 0, the Euler scheme rounds the new
  and the old storage differently, and the correction's source form
  V div(C) telescopes only to round-off) are booked in solverDefect instead
  of drifting the closure (review of M4, rounds 2 and 3; plan section 27,
  item 1, and section 29, item 2).  Every cell is accumulated in long
  double and the cells in a compensated sum, so the rounding of R', which
  repeats in every step of a steady state, stays far below the closure.
  Collective (the fluxes of processor faces, the divergence of C).
  --------------------------------------------------------------------------*/
  const scalarField& rhoI = rho_.primitiveField();
  const scalarField& Y0 = Y.oldTime().primitiveField();
  scalarField mass(YI.size()), oldMass(YI.size());
  forAll(YI, celli) {
    mass[celli] = rhoI[celli]*YI[celli]*V[celli];
    oldMass[celli] = rhoI[celli]*Y0[celli]*V[celli];
  }
  const tmp<fvScalarMatrix> tDdt(fvm::ddt(rho_, Y));
  List<long double> defect;
  scalarField remainder;
  cellDefect (
    transportEqn, tDdt().source(), oldMass, mass,
    mesh_.time().deltaTValue(), Y, defect, remainder
  );

  forAll(exchangeCells_, k) {
    const label celli = exchangeCells_[k];
    for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
      defect[celli] += static_cast<long double> (
        elementFlow (
          pair.regime[e], pair.a[e], pair.g[e], pair.upper[e], YI[celli]
        )
      );
    }
  }

  compensatedSum total;
  forAll(defect, celli) {
    defect[celli] += static_cast<long double>(pair.Srel[celli]*V[celli]);
    total.add(static_cast<scalar>(defect[celli]));
  }
  pair.defect = globalSum(total);

  const tmp<surfaceScalarField> tFlux(eqn.flux());
  forAll(patchIds, k) {
    pair.transport[k] = sum(tFlux().boundaryField()[patchIds[k]]);
  }
  Pstream::listCombineReduce(pair.transport, plusEqOp<scalar>());

  /* the explicit remainder E of the transport-only matrix (unbooked: the
     closure changes by dt sum_c E_c; 0 for the solver's schemes;
     matrixDefect.H) and the explicit source it is left of, sum_c |b - s|
     (the face-flux correction, in source form), and the transport identity
     D of the cell defects, summed exactly (0 in the flux form, the
     double-precision residue of the correction's source form otherwise;
     matrixDefect.H), for the debug report.  Collective. */
  if (reportMatrixResidual_) {
    pair.residualOF = gSum(eqn.residual());
    const scalarField& b = transportEqn.source();
    const scalarField& s = tDdt().source();
    compensatedSum explicitSum, explicitSource;
    label nExplicit = 0;
    forAll(remainder, celli) {
      explicitSource.add(mag(b[celli] - s[celli]));
      if (remainder[celli] != 0) {
        explicitSum.add(remainder[celli]);
        ++nExplicit;
      }
    }
    pair.explicitRemainder = globalSum(explicitSum);
    pair.nExplicitRemainder = returnReduce(nExplicit, sumOp<label>());
    pair.explicitSource = globalSum(explicitSource);

    long double D = 0, terms = 0;
    transportIdentity (
      transportEqn, tDdt().source(), oldMass, mass,
      mesh_.time().deltaTValue(), Y, D, terms
    );
    pair.transportIdentity = static_cast<scalar>(D);
    pair.identityTerms = static_cast<scalar>(terms);
    pair.identityRelative =
      terms > 0 ? static_cast<scalar>(std::fabs(D)/terms) : 0;
  }
}


/*------------------------------------------------------------------------------
Guard (E5): the bounded law at the solved gas values.  An element violates
it when its realised flow differs from the admissible flow of its implicit
flow q = a - g Y (admissibleFlow): an IMPLICIT element outside its bounds,
a CAPPED element with q < upper (it would not exhaust its condensate, e.g.
a sample cell whose gas is supersaturated), a NONE element with a bound it
should take or leave (a sample element in undersaturated gas, a bare wall
in supersaturated gas).  The mass of the violations is |admissible -
realised| dt; for IMPLICIT elements this is the bound violation of the
monotone guard of the plan, so its sum is unchanged wherever the other
elements are complementary.
------------------------------------------------------------------------------*/

bool LESTO::interfaceExchange::guard (
  const label           pairi,
  const volScalarField& Y
) {

  pairData& pair = pairs_[pairi];
  const scalar dt = mesh_.time().deltaTValue();
  const scalar tolerance = massTolerance(pairi);

  /* the violation of every element [kg]: of an IMPLICIT one its bound
     violation, of a CAPPED or NONE one its missed mass */
  scalarList missed(elementCell_.size(), Zero);
  scalar violation = 0;
  forAll(elementCell_, e) {
    const scalar Yc = Y[elementCell_[e]];
    const scalar flow = elementFlow (
      pair.regime[e], pair.a[e], pair.g[e], pair.upper[e], Yc
    );
    const scalar q = pair.a[e] - pair.g[e]*Yc;
    missed[e] = mag(admissibleFlow(q, pair.lower[e], pair.upper[e]) - flow)*dt;
    violation += missed[e];
  }
  reduce(violation, sumOp<scalar>());

  if (violation <= tolerance || pair.nGuardPasses >= nGuard_) {
    return false;
  }

  /*--------------------------------------------------------------------------
  Switch the violating elements (review of M5, round 2: the monotone guard
  never switched back, so a sample left CAPPED in supersaturated gas, or
  NONE in undersaturated gas, by an unsettled corrector loop was accepted):

    IMPLICIT outside its bounds   -> CAPPED or NONE (boundedRegime), as the
                                     predictor would at this Y;
    CAPPED or NONE, violating     -> IMPLICIT, never directly to the other
                                     bound; the largest misses first, until
                                     the misses left add up to at most the
                                     tolerance (switchBackThreshold()).

  The guard is triggered by the total violation, and it acts on the total
  (review of M5, round 3): the misses left after a pass are at most
  guardTolerance x inventory in sum, so a pass is never triggered without
  a switch.  Round 2 switched back only elements whose own miss exceeded
  the tolerance; several smaller ones were accepted without a guard pass
  (the channel of M5.2b: two bare wall elements just above Y_eq, 1.43
  times the tolerance together, in the step ending at 0.167 s).  Where the
  misses left by the old rule add up to at most the tolerance, the same
  elements are switched as before.  The smallest misses go last: an
  element whose implicit flow is at its bound to round-off (a bare wall at
  Y_eq, a noise-level Y ahead of a gas front) would only return with a
  bound violation of the same size, left to the projection.

  The second rule is the adjacency of the predictor (acrossBounds): with
  lambda deltaT of 1e3 to 1e4 the gas value of a CAPPED cell and that of
  the same cell NONE differ by far more than the width of the bounds, so a
  direct switch to the regime admissible at the solved Y would alternate
  between the two bounds and never reach the implicit solution between
  them.  A single element that leaves CAPPED (or NONE) for IMPLICIT stays
  IMPLICIT or moves on to the other bound in the next pass, which then is
  complementary.

  Both rules assume that the passes solve the same linear system with only
  Su and Sp changed: switching one violating IMPLICIT element of a fixed
  matrix A to its bound gives the implicit flow q_N = q* (1 + g (A^-1)_ee),
  on the same side of the bound as q*.  A guard re-solve therefore solves
  the transport-only matrix of the recorded solve again
  (solvePairedGasSpecies.H).  Re-assembled, the SuperBee limiter was
  re-evaluated at the solved Y between the passes, and near a bound a single
  element cycled: a bare wall element at Y_eq on a moving gas front went
  IMPLICIT, NONE, IMPLICIT until nGuard was used up, and the step was
  accepted with a miss 1e5 times the tolerance (review of M5, round 4).
  --------------------------------------------------------------------------*/
  DynamicList<scalar> misses;
  forAll(elementCell_, e) {
    if (pair.regime[e] != IMPLICIT && missed[e] > 0) {
      misses.append(missed[e]);
    }
  }
  List<scalarList> perRank(UPstream::nProcs());
  perRank[UPstream::myProcNo()] = misses;
  Pstream::gatherList(perRank);
  scalar threshold = 0;
  if (UPstream::master()) {
    std::vector<double> all;
    for (const scalarList& rankMisses : perRank) {
      all.insert(all.end(), rankMisses.begin(), rankMisses.end());
    }
    threshold = switchBackThreshold(all, tolerance);
  }
  Pstream::broadcast(threshold);

  label nSwitched = 0;

  forAll(exchangeCells_, k) {
    bool changed = false;
    for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
      const scalar Yc = Y[elementCell_[e]];
      const scalar q = pair.a[e] - pair.g[e]*Yc;
      exchangeRegime r = exchangeRegime(pair.regime[e]);
      if (r == IMPLICIT) {
        r = boundedRegime(q, pair.lower[e], pair.upper[e]);
      } else if (missed[e] > 0 && missed[e] >= threshold) {
        r = IMPLICIT;
      }
      if (r != pair.regime[e]) {
        pair.regime[e] = r;
        changed = true;
        ++nSwitched;
      }
    }
    if (changed) {
      cellCoefficients(pair, k);
    }
  }

  reduce(nSwitched, sumOp<label>());

  if (nSwitched == 0) {
    return false;
  }

  ++pair.nGuardPasses;
  pair.nGuardSwitched += nSwitched;
  return true;
}


/*------------------------------------------------------------------------------
Realised exchange of the step (E5 projection, E6 condensate, E11 monitors)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::realise (
  const label     pairi,
  volScalarField& Y,
  const label     nOuter,
  const bool      converged,
  const scalar    initialResidual
) {

  pairData& pair = pairs_[pairi];
  const word& fieldName = pair.solverName;
  const scalar dt = mesh_.time().deltaTValue();
  const scalar significant = massTolerance(pairi);
  const scalarField& V = mesh_.V();
  const scalarField& rhoI = rho_.primitiveField();
  scalarField& YI = Y.primitiveFieldRef();

  /*--------------------------------------------------------------------------
  Realised flows of the last solved matrix at the final Y.
  --------------------------------------------------------------------------*/
  forAll(elementCell_, e) {
    pair.flow[e] = elementFlow (
      pair.regime[e], pair.a[e], pair.g[e], pair.upper[e],
      YI[elementCell_[e]]
    );
  }

  /*--------------------------------------------------------------------------
  Conservative same-cell projection of what the guard left (round-off, or
  guard passes exhausted): the inadmissible part of the flow is taken back
  out of the gas cell, so gas + condensate stay exact.  Y is never clipped.
  --------------------------------------------------------------------------*/
  label nProjected = 0, nProjectedSignificant = 0;
  scalar projected = 0;

  forAll(elementCell_, e) {
    if (pair.regime[e] == IMPLICIT) {
      const scalar q = pair.flow[e];
      const scalar qa = min(max(q, pair.lower[e]), pair.upper[e]);
      if (qa != q) {
        const label celli = elementCell_[e];
        YI[celli] += (qa - q)*dt/(rhoI[celli]*V[celli]);
        projected += mag(qa - q)*dt;
        if (mag(qa - q)*dt > significant) {
          ++nProjectedSignificant;
        }
        pair.flow[e] = qa;
        ++nProjected;
      }
    }
  }

  reduce(nProjected, sumOp<label>());
  reduce(nProjectedSignificant, sumOp<label>());
  reduce(projected, sumOp<scalar>());

  if (nProjected > 0) {
    Y.correctBoundaryConditions();
  }

  pair.projectedStep = projected;
  pair.nProjectedRun += nProjected;
  pair.nProjectedSignificantRun += nProjectedSignificant;
  pair.nProjectedAboveRun += projected > significant;
  pair.nGuardPassesRun += pair.nGuardPasses;

  /*--------------------------------------------------------------------------
  Condensate (plan E6): the wall deposit m''^{n+1} = m''^n - dt q/A, the
  reservoir of a sample element Cs^{n+1} = Cs^n - dt q/V.  A CAPPED element
  is set to exactly zero; its residue must be round-off (review item B2) and
  is booked as clamped.  The regimes are counted for the wall elements and
  for the sample elements of this pair's inventory samples.

  The stored condensate of every other element changes by the rounded
  update, not by exactly dt q: the residue (the mass held before, minus
  the mass held after, minus dt q, both masses computed as the inventories
  compute them, condensate*measure) is booked as clamped too (review of
  M4, round 2; plan section 27, item 1).  It is round-off, but systematic
  where a deposit grows by the same increment for many steps, a few ulp of
  the deposit apart: the increment is rounded the same way at every step,
  and the closure of the M8-scale channel drifted to -1.2e-13 of the
  inventory in 68000 steps of a steady release (M2.2h).  A residue of
  exactly 0 (no flow) adds nothing.
  --------------------------------------------------------------------------*/
  scalar clamped = 0;
  FixedList<label, 3> nWallRegime(0), nSampleRegime(0);   /* NONE,
                                                             IMPLICIT, CAPPED */

  forAll(elementCell_, e) {

    const bool wall = e < nWall_;
    if (!wall && samples_[elementSample_[e - nWall_]].pair != pairi) {
      continue;
    }

    /* the condensate of the element [kg] and its measure (A or V) */
    const scalar measure = wall ? elementArea_[e] : V[elementCell_[e]];
    scalar& condensate = wall ? pair.deposit[e] : pair.reservoir[e - nWall_];
    FixedList<label, 3>& count = wall ? nWallRegime : nSampleRegime;

    if (pair.regime[e] == CAPPED) {
      const scalar held = condensate*measure;
      const scalar residue = held - dt*pair.flow[e];
      if (mag(residue) > 1e-12*held) {
        FatalErrorInFunction
          << fieldName << ": the CAPPED " << (wall ? "wall" : "sample")
          << " element " << e << " (cell " << elementCell_[e] << ") leaves "
          << "the residue " << residue << " kg of its condensate " << held
          << " kg, more than round-off." << abort(FatalError);
      }
      clamped += residue;
      condensate = 0;
    } else {
      const scalar held = condensate*measure;
      condensate -= dt*pair.flow[e]/measure;
      clamped += (held - condensate*measure) - dt*pair.flow[e];
    }
    ++count[pair.regime[e] == CAPPED ? 2 : pair.regime[e] == IMPLICIT ? 1 : 0];
  }

  reduce(clamped, sumOp<scalar>());
  for (label i = 0; i < 3; ++i) {
    reduce(nWallRegime[i], sumOp<label>());
    reduce(nSampleRegime[i], sumOp<label>());
  }
  const label nNone = nWallRegime[0];
  const label nImplicit = nWallRegime[1];
  const label nCapped = nWallRegime[2];
  pair.clampedStep = clamped;

  /*--------------------------------------------------------------------------
  Complementarity (review, non-blocking 3): elements that are not IMPLICIT
  and whose realised flow differs from the admissible flow of their
  implicit flow at the final Y, qc = min(max(a - g Y, lower), upper)
  (admissibleFlow), with the missed mass |qc - flow| dt: e.g. a NONE bare
  wall with Y > Y_eq (missed deposition), a CAPPED element whose implicit
  flow no longer exhausts it, also a CAPPED sample element in
  supersaturated gas (qc = 0) and a NONE sample element in undersaturated
  gas (qc > 0; review of M5, round 2: before, only elements whose implicit
  flow lay within the bounds were counted, which missed both).  Counted and
  reported, not changed; the guard has switched them unless they are
  round-off or its passes were exhausted.
  --------------------------------------------------------------------------*/
  label nComplementarity = 0, nComplementaritySignificant = 0;
  scalar complementarity = 0;

  forAll(elementCell_, e) {
    if (pair.regime[e] != IMPLICIT && (pair.a[e] != 0 || pair.g[e] != 0)) {
      const scalar qi = pair.a[e] - pair.g[e]*YI[elementCell_[e]];
      const scalar missed =
        mag(admissibleFlow(qi, pair.lower[e], pair.upper[e]) - pair.flow[e]);
      if (missed != 0) {
        ++nComplementarity;
        complementarity += missed*dt;
        if (missed*dt > significant) {
          ++nComplementaritySignificant;
        }
      }
    }
  }

  reduce(nComplementarity, sumOp<label>());
  reduce(nComplementaritySignificant, sumOp<label>());
  reduce(complementarity, sumOp<scalar>());
  pair.nComplementarityRun += nComplementarity;
  pair.nComplementaritySignificantRun += nComplementaritySignificant;
  pair.nComplementarityAboveRun += complementarity > significant;

  /*--------------------------------------------------------------------------
  Monitors (E11): min Y, negative mass, carrier-validity mole fraction
  x = c/(c + p/(R T)), c = rho Y/M the molar concentration of the species
  and p/(R T) that of the carrier (of the frozen p and T; the same as
  (Y/M)/(Y/M + 1/M_carrier) for a perfect-gas carrier, but in paper mode,
  rho = 1, that form measured the species against an imaginary helium of
  1 kg/m3 and was 6 to 25 times too small; review of M4, round 3), and the
  mass that entered the gas in cells with x > moleFractionWarning.  The
  maximum starts at 0: where Y <= 0 everywhere (solver noise ahead of a gas
  front) there is no gas to load the carrier, and a negative "mole
  fraction" means nothing.
  --------------------------------------------------------------------------*/
  const scalar M = pair.molarMass;
  scalar negativeMass = 0;
  scalar maxMoleFraction = 0;
  label nAbove = 0;
  boolList above(mesh_.nCells(), false);

  forAll(YI, celli) {
    negativeMass += min(rhoI[celli]*YI[celli], scalar(0))*V[celli];
    const scalar c = rhoI[celli]*YI[celli]/M;
    const scalar x = c/(c + carrierConcentration_[celli]);
    maxMoleFraction = max(maxMoleFraction, x);
    if (x > moleFractionWarning_) {
      above[celli] = true;
      ++nAbove;
    }
  }

  scalar entered = 0;
  if (returnReduceOr(nAbove > 0)) {
    forAll(elementCell_, e) {
      if (above[elementCell_[e]] && pair.flow[e] > 0) {
        entered += pair.flow[e]*dt;
      }
    }
    forAll(YI, celli) {
      if (above[celli]) {
        entered += pair.Srel[celli]*V[celli]*dt;
      }
    }
  }

  reduce(negativeMass, sumOp<scalar>());
  reduce(maxMoleFraction, maxOp<scalar>());
  reduce(nAbove, sumOp<label>());
  reduce(entered, sumOp<scalar>());
  pair.evaporatedAboveWarningStep = entered;

  if (!converged) {
    ++pair.nUnconverged;
  }
  if (pair.finalSolver) {
    pair.nFinalChangesSignificantRun += pair.nChangedSignificant;
  }

  /*--------------------------------------------------------------------------
  Log.  For paired species the outer-convergence status is part of the
  exchange line (no WarningInFunction per step); end() counts the steps.
  The regime changes are those of the corrector that decided the outer
  loop and, if a Final solve followed, those of its own predictor.
  --------------------------------------------------------------------------*/
  Info<< fieldName << " exchange: implicit " << nImplicit << ", capped "
    << nCapped << ", none " << nNone << " elements";
  if (pair.inventory) {
    Info<< " (sample: implicit " << nSampleRegime[1] << ", capped "
      << nSampleRegime[2] << ", none " << nSampleRegime[0] << ")";
  }
  Info<< "; outer iterations "
    << nOuter << (converged ? ", converged" : ", NOT converged")
    << " (initial residual " << initialResidual << ", tolerance "
    << tolerance_ << "; regime changes " << pair.nChangedOuter << ", "
    << pair.changedMassOuter << " kg, " << pair.nChangedSignificantOuter
    << " significant)";
  if (pair.finalSolver) {
    Info<< "; final predictor: regime changes " << pair.nChanged << ", "
      << pair.changedMass << " kg, " << pair.nChangedSignificant
      << " significant";
  }
  /* a total above the tolerance is named (the guard left more than it
     accepts: its passes were exhausted, or nGuard 0) */
  auto aboveTolerance = [significant](const scalar total) {
    OStringStream os;
    if (total > significant) {
      os << " > tolerance " << significant << " kg";
    }
    return os.str();
  };
  Info<< "; guard passes "
    << pair.nGuardPasses << " (" << pair.nGuardSwitched << " elements); "
    << "projected " << nProjected << " (" << projected << " kg"
    << aboveTolerance(projected).c_str() << ", " << nProjectedSignificant
    << " significant); complementarity " << nComplementarity << " ("
    << complementarity << " kg" << aboveTolerance(complementarity).c_str()
    << ", " << nComplementaritySignificant << " significant)" << nl;

  /* equilibrium GEMS, mode local: the update of this step (M9) */
  if (gems_ && gems_->local()) {
    Info<< fieldName << ' ' << gems_->stepReport(pairi).c_str() << nl;
  }

  Info<< fieldName << " transport [kg/s]:";
  forAll(ledger_->patchIds(), k) {
    Info<< ' ' << mesh_.boundary()[ledger_->patchIds()[k]].name() << ' '
      << pair.transport[k];
  }
  Info<< "; solver defect " << pair.defect << " kg/s";
  if (reportMatrixResidual_) {
    Info<< "; gSum(fvMatrix::residual()) " << pair.residualOF
      << " kg/s (not used); explicit remainder " << pair.explicitRemainder
      << " kg/s in " << pair.nExplicitRemainder << " cells (of an explicit "
      << "source sum |b - s| " << pair.explicitSource << " kg/s); "
      << "transport identity " << pair.transportIdentity << " kg/s (|D| "
      << pair.identityRelative << " of its terms, " << pair.identityTerms
      << " kg/s)";
  }
  Info<< nl;

  Info<< fieldName << " monitors: min Y " << gMin(YI) << ", negative mass "
    << negativeMass << " kg, max mole fraction " << maxMoleFraction
    << " (" << nAbove << " cells above " << moleFractionWarning_ << ", "
    << entered << " kg entered there)" << nl;
}


/*------------------------------------------------------------------------------
Ledger, closure and derived condensate fields (E8)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::balance() {

  if (mock()) {
    return;
  }

  const Time& runTime = mesh_.time();
  const scalar dt = runTime.deltaTValue();
  phaseChangeLedger& ledger = *ledger_;
  const scalarField& V = mesh_.V();
  const scalarField& rhoI = rho_.primitiveField();

  /* the ledger now holds the balance at the end of this step */
  ledger.time() = stepEnd_;

  forAll(pairs_, pairi) {

    pairData& pair = pairs_[pairi];
    const scalar M = pair.molarMass;

    /* cumulative entries: compensated sums (phaseChangeLedger.H) */
    ledger.add(pairi, phaseChangeLedger::RELEASED, pair.releasedStep);
    forAll(pair.transport, k) {
      ledger.addTransport(pairi, k, pair.transport[k]*dt);
    }
    ledger.add(pairi, phaseChangeLedger::SOLVER_DEFECT, pair.defect*dt);
    ledger.add(pairi, phaseChangeLedger::CLAMPED, pair.clampedStep);
    ledger.add(pairi, phaseChangeLedger::PROJECTED, pair.projectedStep);
    ledger.add (
      pairi,
      phaseChangeLedger::EVAPORATED_ABOVE_WARNING,
      pair.evaporatedAboveWarningStep
    );

    const scalar gas = gasInventory(pair);
    const scalar wall = wallInventory(pair);
    const scalar sample = sampleInventory(pair);
    pair.elementGas = gas;
    pair.elementWall = wall;
    pair.elementSample = sample;
    const scalar held = gas + wall + sample;
    ledger.set(pairi, phaseChangeLedger::HELD, held);

    const scalar closure = held - ledger.supplied(pairi);
    const scalar reference = max(ledger.reference(pairi), VSMALL);
    const scalar solverDefect = ledger(pairi, phaseChangeLedger::SOLVER_DEFECT);

    /*------------------------------------------------------------------------
    Derived Y_<condensate> = (sum_e m''_e A_e + Cs V)/(rho V), before
    runTime.write(), so that Y_ and c_ of the condensate are consistent
    (plan E6; Cs the reservoir of an inventory sample cell).
    ------------------------------------------------------------------------*/
    volScalarField& Ys = species_[pair.condensed];
    scalarField& YsI = Ys.primitiveFieldRef();
    YsI = 0;
    for (label e = 0; e < nWall_; ++e) {
      YsI[elementCell_[e]] += pair.deposit[e]*elementArea_[e];
    }
    for (label e = nWall_; e < elementCell_.size(); ++e) {
      YsI[elementCell_[e]] += reservoirMass(pair, e);
    }
    for (const label celli : exchangeCells_) {
      YsI[celli] /= rhoI[celli]*V[celli];
    }
    Ys.correctBoundaryConditions();

    Info<< "Balance " << speciesNames_[pair.gas] << " [mol]: gas " << gas/M
      << " wall " << wall/M << " sample " << sample/M << " released "
      << ledger(pairi, phaseChangeLedger::RELEASED)/M << " transport "
      << ledger.totalTransport(pairi)/M << " removed "
      << ledger(pairi, phaseChangeLedger::REMOVED)/M << " clamped "
      << ledger(pairi, phaseChangeLedger::CLAMPED)/M << " restart "
      << ledger(pairi, phaseChangeLedger::RESTART)/M << " solverDefect "
      << solverDefect/M << " (relative " << solverDefect/reference
      << ") closure " << closure/M << " (relative " << closure/reference
      << ")" << nl;

    ledger.writeLine(pairi, gas, wall, sample, closure);

    /*------------------------------------------------------------------------
    Deposit fields mw_<condensate> and mDep_<condensate> (output),
    refreshed every step so that function objects sampling them between
    write times see the current deposit; written by runTime.write() with
    the fields (review item B4).
    ------------------------------------------------------------------------*/
    updateDepositField(pair);
    updateReservoirField(pair);

    /*------------------------------------------------------------------------
    Realised exchange of the step per cell volume, exch_<gas> [kg/(m3 s)]
    (output, optional): sum over the elements of a cell of the realised
    flow q_e of realise(), divided by V_c; positive = evaporation.
    ------------------------------------------------------------------------*/
    if (pair.exchangeField) {
      volScalarField& exch = pair.exchangeField();
      scalarField& exchI = exch.primitiveFieldRef();
      exchI = 0;
      forAll(elementCell_, e) {
        exchI[elementCell_[e]] += pair.flow[e];
      }
      for (const label celli : exchangeCells_) {
        exchI[celli] /= V[celli];
      }
      exch.correctBoundaryConditions();
    }

    /* equilibrium GEMS with writeEquilibriumField: p_eq of the step (M9) */
    if (pair.equilibriumField) {
      updateEquilibriumFields(pair);
    }
  }

  /*--------------------------------------------------------------------------
  The profile times reached by this step (written after runTime.write(), or
  by it at a write time) and their record in the layout, so that the state
  written by runTime.write() records them (restart: carry-over (d)).
  --------------------------------------------------------------------------*/
  if (respeciation_) respeciation_->book();
  writeElementBalance();

  if (profiles_) {
    const string reached = profiles_->reached(stepEnd_, stepDeltaT_);
    profilesPending_ += string(profilesPending_.empty() || reached.empty()
                               ? "" : "; ") + reached;
    recordProfileTimes();
  }

  /*--------------------------------------------------------------------------
  The ledger (updated above), the element regimes and the deposits now hold
  the end of this step.  They are registered AUTO_WRITE (binaryIOList.H),
  so every write of Time stores them with the fields: runTime.write() of
  this step, and a write started by a function object before the next.
  --------------------------------------------------------------------------*/
  storeElementState();
}


/*------------------------------------------------------------------------------
mw_<condensate> from the element deposits: the deposits on the interface
patches, the area mean of the elements of a cell in the internal field, and
the processor patches from the neighbour cells
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::updateDepositField(pairData& pair) const {
  writeWallStore(pair.depositField(), pair.wallDepositField(), pair.deposit);
}

void LESTO::interfaceExchange::writeWallStore(volScalarField& mw, volScalarField& mDep, const scalarField& deposit) const {
  if(deposit.size()==nWall_ && elementArea_.size()>nWall_) {
    // A shared wall store has no inventory on a legacy pair's boat samples.
    scalarField full(elementArea_.size(),0);
    forAll(deposit,e)full[e]=deposit[e];
    writeWallStore(mw,mDep,full);
    return;
  }
  scalarField& mwI = mw.primitiveFieldRef();
  scalarField& mDepI = mDep.primitiveFieldRef();
  mwI = 0;
  mDepI = 0;
  for (label k = 0; k < nWallCells_; ++k) {
    scalar area = 0, wallArea = 0, mass = 0;
    for (label e = cellStart_[k]; e < cellStart_[k+1]; ++e) {
      area += elementArea_[e];
      wallArea += elementWallArea_[e];
      mass += deposit[e]*elementArea_[e];
    }
    mwI[exchangeCells_[k]] = mass/area;
    mDepI[exchangeCells_[k]] = mass/wallArea;
  }

  /*--------------------------------------------------------------------------
  layer distance: one element per cell, and the cell value of mw_ is the
  restart record of its deposit (the fallback after a re-decomposition), so
  it holds m''_e itself, not (m''_e A_e)/A_e, which differs from m''_e by an
  ulp for many values (review of M4, round 1: every same-decomposition
  restart warned that mw_ differed from uniform/phaseChange/
  phaseChangeDeposits, and a re-decomposed restart re-read rounded
  deposits).
  --------------------------------------------------------------------------*/
  if (distanceLayer()) {
    for (label k = 0; k < nWallCells_; ++k) {
      mwI[exchangeCells_[k]] = deposit[cellStart_[k]];
    }
  }

  /*--------------------------------------------------------------------------
  Interface patches: mw_ holds the element deposit m''_e, mDep_ the sum
  over the elements of a face of m''_e A_e/|S_f|, evaluated as
  m''_e (A_e/|S_f|) so that a power-of-two areaMultiplier scales it
  exactly.  (layer faceCells has one element per face.)

  layer distance: the elements are cells.  mw_ is 0 on the interface
  patches (its deposits are the cell values), and mDep_ on face f is the
  mass sum_{e -> f} m''_e A_e of the elements whose nearest wall face is f,
  on any rank (nearestFaceSums(), collective), over |S_f|.
  --------------------------------------------------------------------------*/
  for (const label patchi : interfacePatches_) {
    mDep.boundaryFieldRef()[patchi] = 0;
  }

  if (distanceLayer()) {

    scalarField elementMass(nWall_);
    forAll(elementMass, e) {
      elementMass[e] = deposit[e]*elementArea_[e];
    }
    const tmp<scalarField> tFaceMass(nearestFaceSums(elementMass));
    const scalarField& faceMass = tFaceMass();
    const surfaceScalarField::Boundary& magSfb =
      mesh_.magSf().boundaryField();

    for (const label patchi : interfacePatches_) {
      mw.boundaryFieldRef()[patchi] = 0;
    }
    forAll(faceMass, i) {
      const label patchi = interfaceFacePatch_[i];
      const label facei = interfaceFaceFace_[i];
      mDep.boundaryFieldRef()[patchi][facei] =
        faceMass[i]/magSfb[patchi][facei];
    }

  } else {

    for (label e = 0; e < nWall_; ++e) {
      const label patchi = elementPatch_[e];
      const label facei = elementFace_[e];
      mw.boundaryFieldRef()[patchi][facei] = deposit[e];
      mDep.boundaryFieldRef()[patchi][facei] +=
        deposit[e]*(elementArea_[e]/elementWallArea_[e]);
    }
  }

  /* calculated patches keep the values just set; processor patches take
     the neighbour's cells */
  mw.correctBoundaryConditions();
  mDep.correctBoundaryConditions();
}


/*------------------------------------------------------------------------------
Cs_<condensate> from the reservoirs of the inventory samples of a pair: the
reservoir in the sample cells, 0 elsewhere; processor patches from the
neighbour cells
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::updateReservoirField(pairData& pair) const {

  if (!pair.reservoirField) {
    return;
  }

  volScalarField& Cs = pair.reservoirField();
  scalarField& CsI = Cs.primitiveFieldRef();
  CsI = 0;
  const label pairi = pairOfSpecies_[pair.gas];
  for (label e = nWall_; e < elementCell_.size(); ++e) {
    if (samples_[elementSample_[e - nWall_]].pair == pairi) {
      CsI[elementCell_[e]] = pair.reservoir[e - nWall_];
    }
  }
  Cs.correctBoundaryConditions();
}


/*------------------------------------------------------------------------------
equilibrium GEMS, mode local: the lagged state of a restart (plan section
35, item C).  uniform/phaseChange/phaseChangeEquilibrium holds the
fingerprint of all elements (as exact scalars), then per pair the counts
of the last update (gemsEquilibrium::lastUpdate(), nPairState values),
then per pair and element p_eq [Pa], chi and the source of the element (-1
for a sample element of another pair, which keeps a = g = 0).  Onto the
same elements every element's law is recomputed from it as the update did
(hksCoefficients(): the same doubles), so the steps up to the next update
repeat the continuous run bit for bit.  With speciation none chi is 1
whatever the file holds.
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::restoreEquilibrium (
  const stateFile file,
  const string&   why,
  const word&     startName
) {

  gemsEquilibrium& gems = gems_();
  const scalarList& state = *equilibriumIO_;
  const label nPairs = pairs_.size();
  const label nElements = elementCell_.size();
  const label nHeader = elementsFingerprint_.size();
  const label nCounts = nPairs*gemsEquilibrium::nPairState;
  const label nState = nHeader + nCounts + 3*nPairs*nElements;

  bool fits = state.size() == nState;
  for (label i = 0; fits && i < nHeader; ++i) {
    fits = state[i] == scalar(elementsFingerprint_[i]);
  }
  for (label i = nHeader + nCounts + 2; fits && i < nState; i += 3) {
    fits = state[i] == std::floor(state[i]) && state[i] >= -1
        && state[i] < scalar(gemsEquilibrium::nSourceTypes);
  }

  const label next = mesh_.time().timeIndex() + 1;
  if (returnReduceAnd(fits)) {

    gems.restoreLastUpdate(SubList<scalar>(state, nCounts, nHeader));
    const bool halfCell = wallResistance_ == "halfCell";

    forAll(pairs_, k) {
      pairData& pair = pairs_[k];
      forAll(elementCell_, e) {
        const label own =
          e < nWall_ ? -1 : samples_[elementSample_[e - nWall_]].pair;
        if (own >= 0 && own != k) {
          continue;
        }
        const scalar* v =
          state.cdata() + nHeader + nCounts + 3*(k*nElements + e);
        pair.pEq[e] = v[0];
        pair.chi[e] = gems.speciation() ? v[1] : 1.0;
        pair.source[e] = label(v[2]);
        hksCoefficients(pair, e, halfCell);
      }
      if (pair.equilibriumField) {
        updateEquilibriumFields(pair);
      }
    }

    label update = next;
    while (!gems.updateDueAt(update)) {
      ++update;
    }
    Info<< "GEMS3K (mode local): the lagged state of the last update (p_eq"
      << (gems.speciation() ? ", chi" : "") << " and the source of every "
      << "element) restored from uniform/phaseChange/phaseChangeEquilibrium;"
      << " the next update at the step of time index " << update << nl;

  } else if (!gems.updateDueAt(next)) {

    WarningInFunction
      << "uniform/phaseChange/phaseChangeEquilibrium at time " << startName
      << " is "
      << (
           file == stateFile::FOREIGN ? why
         : file == stateFile::ABSENT ? string("missing")
         : string("empty or written for other interface elements (another "
             "decomposition or a renumbered mesh?)")
         ).c_str()
      << ": the first step updates p_eq of every element (equilibrium "
      << "GEMS, mode local), although updateInterval "
      << gems.updateInterval() << " schedules no update at the step of "
      << "time index " << next << "; the schedule follows the time index "
      << "from then on." << nl << endl;
  }
}


/*------------------------------------------------------------------------------
Per-rank element state of uniform/: the fingerprint, then the regimes
(phaseChangeRegimes) or the deposits (phaseChangeDeposits, the fingerprint
as exact scalars) of every element of every pair
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::storeElementState() {

  recordUniformValues();

  const label nElements = elementCell_.size();
  const label nSample = nElements - nWall_;
  const label nHeader = fingerprint_.size();
  labelList& regimes = *regimesIO_;
  scalarList& deposits = *depositsIO_;

  regimes.resize(nHeader + pairs_.size()*nElements);
  deposits.resize(nHeader + pairs_.size()*nWall_);

  forAll(fingerprint_, i) {
    regimes[i] = elementsFingerprint_[i];
    deposits[i] = scalar(fingerprint_[i]);
  }
  forAll(pairs_, pairi) {
    const pairData& pair = pairs_[pairi];
    forAll(pair.regime, e) {
      regimes[nHeader + pairi*nElements + e] = pair.regime[e];
    }
    forAll(pair.deposit, e) {
      deposits[nHeader + pairi*nWall_ + e] = pair.deposit[e];
    }
  }

  if (reservoirsIO_) {
    scalarList& reservoirs = *reservoirsIO_;
    reservoirs.resize(nHeader + pairs_.size()*nSample);
    forAll(samplesFingerprint_, i) {
      reservoirs[i] = scalar(samplesFingerprint_[i]);
    }
    forAll(pairs_, pairi) {
      const pairData& pair = pairs_[pairi];
      forAll(pair.reservoir, k) {
        reservoirs[nHeader + pairi*nSample + k] = pair.reservoir[k];
      }
    }
  }

  /* equilibrium GEMS, mode local: the lagged state (restoreEquilibrium());
     empty until the first update */
  if (equilibriumIO_) {
    scalarList& state = *equilibriumIO_;
    if (gems_->lagged()) {
      const label nAll = elementsFingerprint_.size();
      const scalarList counts(gems_->lastUpdate());
      const label nCounts = counts.size();
      state.resize(nAll + nCounts + 3*pairs_.size()*nElements);
      forAll(elementsFingerprint_, i) {
        state[i] = scalar(elementsFingerprint_[i]);
      }
      forAll(counts, i) {
        state[nAll + i] = counts[i];
      }
      forAll(pairs_, pairi) {
        const pairData& pair = pairs_[pairi];
        forAll(elementCell_, e) {
          scalar* v =
            state.data() + nAll + nCounts + 3*(pairi*nElements + e);
          v[0] = pair.pEq[e];
          v[1] = pair.chi[e];
          v[2] = scalar(pair.source[e]);
        }
      }
    } else {
      state.clear();
    }
  }
}


/*------------------------------------------------------------------------------
Record of the rounded 'uniform' entries of the fields that a restart re-reads
(uniform/phaseChangeLayout, merged over the ranks)
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::recordUniformValues() {

  /*--------------------------------------------------------------------------
  OpenFOAM writes a field or patch whose values are all equal as 'uniform
  <value>' in text at writePrecision, also in a binary file.  Whether such
  an entry loses digits is decided here, on the values about to be written
  (roundedUniform); a restart cannot tell, since the value it reads always
  reproduces its own text.  'checked' names the fields, 'rounded' the
  internal fields (<field>) and patches (<field>:<patch>, and
  <field>:processor for the processor patches of phi) whose entry is
  rounded on some rank.  Checked are the paired gases and, when the
  carrier is written (writeFrozenFields yes), rho, T, p and phi.
  Patches whose value is recomputed when read (valueReadBack(): e.g.
  zeroGradient, fixedGradient, mixed, symmetryPlane, cyclic) and the
  processor patches of the volume fields (evaluated at start-up) are not
  recorded.  Called with the element state before every write
  (storeElementState()).
  --------------------------------------------------------------------------*/
  const label precision = IOstream::defaultPrecision();
  const label nNonProcessor = mesh_.boundaryMesh().nNonProcessor();
  DynamicList<word> checked, rounded;

  auto check = [&](const auto& field, const bool processorPatches) {

    const word& name = field.name();
    checked.append(name);
    if (roundedUniform(field.primitiveField(), precision)) {
      rounded.append(name);
    }

    const auto& bf = field.boundaryField();
    bool processor = false;
    forAll(bf, patchi) {
      if (patchi >= nNonProcessor) {
        processor = processor
          || (processorPatches && roundedUniform(bf[patchi], precision));
      } else if (
          valueReadBack(bf[patchi])
       && roundedUniform(bf[patchi], precision)
      ) {
        rounded.append(word(name + ':' + bf[patchi].patch().name()));
      }
    }
    if (processor) {
      rounded.append(word(name + ":processor"));
    }
  };

  forAll(pairs_, pairi) {
    check(species_[pairs_[pairi].gas], false);
  }
  if (rho_.writeOpt() == IOobject::AUTO_WRITE) {
    check(rho_, false);
    check(thermo_.T(), false);
    check(thermo_.p(), false);
    check(phi_, true);
  }

  /*--------------------------------------------------------------------------
  The layout is one global file (phaseChangeLedger.H: with the collated
  file handler the master writes it for all ranks), so the record holds the
  rounded entries of all ranks, merged in the order of the ranks.  A
  restart uses it as it used the per-rank records before: every item is
  reduced with 'or' over the ranks (readState(), roundedValues), which the
  merged list gives on every rank.  'checked' is the same on all ranks.
  --------------------------------------------------------------------------*/
  List<wordList> perRank(UPstream::nProcs());
  perRank[UPstream::myProcNo()] = rounded;
  Pstream::allGatherList(perRank);
  DynamicList<word> merged;
  for (const wordList& items : perRank) {
    for (const word& item : items) {
      if (!merged.found(item)) {
        merged.append(item);
      }
    }
  }

  dictionary record;
  record.add("checked", wordList(checked));
  record.add("rounded", wordList(merged));
  ledger_->layout().set("uniformValues", record);
}


/*------------------------------------------------------------------------------
Axial profiles at the writes of Time (profileWriter) and at the profile times
(writeProfiles(), after runTime.write())
------------------------------------------------------------------------------*/

bool LESTO::interfaceExchange::profileWriter::writeObject (
  IOstreamOption,
  const bool
) const {

  /*--------------------------------------------------------------------------
  Time calls this for every write, on every rank, with writeTime() true:
  the write times of controlDict (runTime.write() of the solver) and the
  writes of function objects between two steps (Time::writeNow(),
  Time::writeAndEnd()).  An explicit write() of the object by anything
  else (e.g. a writeObjects function object between two write times) is
  ignored.  No file of the time directory is written.
  --------------------------------------------------------------------------*/
  if (time().writeTime()) {
    exchange_.writeProfilesOfState(true);
  }
  return true;
}


void LESTO::interfaceExchange::writeProfiles() {
  writeProfilesOfState(false);
}


/*------------------------------------------------------------------------------
Record of the profile times written so far (uniform/phaseChangeLayout,
profiles { written (...); }, exact decimals), refreshed by balance() before
runTime.write(), so every write of Time stores the entries its state covers
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::recordProfileTimes() {
  List<string> written;
  for (const scalar t : profiles_->written()) {
    written.append(string(phaseChangeLedger::exactDecimal(t)));
  }
  dictionary record;
  record.add("written", written);
  ledger_->layout().set("profiles", record);
}


void LESTO::interfaceExchange::recordFormulas() {
  dictionary record;
  for (const pairData& pair : pairs_) {
    if (pair.formula.empty()) continue;
    dictionary formula;
    for (const auto& atom : pair.formula) formula.add(atom.first, atom.second);
    record.add(speciesNames_[pair.gas], formula);
  }
  const dictionary* old = ledger_->layout().findDict("formulas", keyType::LITERAL);
  if (old) {
    bool same = old->size() == record.size();
    for (const entry& e : record) {
      const dictionary* previous = old->findDict(e.keyword(), keyType::LITERAL);
      if (!previous || previous->size() != e.dict().size()) {
        same = false;
        continue;
      }
      for (const entry& atom : e.dict()) {
        if (!previous->found(atom.keyword())
          || previous->get<scalar>(atom.keyword())
             != e.dict().get<scalar>(atom.keyword())) same = false;
      }
    }
    if (!same) {
      FatalErrorInFunction << "Species formulas changed across restart."
        << exit(FatalError);
    }
  }
  if (record.empty()) return;
  ledger_->layout().set("formulas", record);
  if (dict_.isDict("respeciation") && dict_.subDict("respeciation").getOrDefault<word>("engine","none") != "none") return;
  if (!writeFile_ || !Pstream::master()) return;
  const fileName dir(mesh_.time().globalPath()/"postProcessing"
    /"phaseChangeElements"/mesh_.time().timeName());
  mkDir(dir);
  elementFiles_.resize(elementNames_.size());
  forAll(elementNames_, i) {
    elementFiles_.set(i, new OFstream(dir/(elementNames_[i] + ".dat")));
    OFstream& os = elementFiles_[i];
    os.precision(17);
    os << "# element " << elementNames_[i] << "; amounts [mol of atoms]" << nl
      << "# Derived from pair ledgers with stoichiometric coefficients nu/M; "
      << "no independent accumulated state." << nl
      << "# columns time initial released transport removed clamped solverDefect "
      << "restart gas wall sample held supplied closure pairClosure reference" << nl;
  }
}

void LESTO::interfaceExchange::writeElementBalance() {
  if (respeciation_) { respeciation_->writeElements(); return; }
  if (!writeFile_ || !Pstream::master()) return;
  const phaseChangeLedger& ledger = *ledger_;
  forAll(elementFiles_, ei) {
    FixedList<compensatedSum, 15> totals;
    forAll(pairs_, pi) {
      const pairData& pair = pairs_[pi];
      const auto atom = pair.formula.find(elementNames_[ei]);
      if (atom == pair.formula.end()) continue;
      const scalar factor = atom->second/pair.molarMass;
      const phaseChangeLedger::entry entries[] = {
        phaseChangeLedger::INITIAL, phaseChangeLedger::RELEASED,
        phaseChangeLedger::REMOVED, phaseChangeLedger::CLAMPED,
        phaseChangeLedger::SOLVER_DEFECT, phaseChangeLedger::RESTART
      };
      const int columns[] = {0, 1, 3, 4, 5, 6};
      for (int k = 0; k < 6; ++k) totals[columns[k]].add(factor*ledger(pi, entries[k]));
      totals[2].add(factor*ledger.totalTransport(pi));
      // Inventories require reductions: ledger HELD is already global. Pair
      // gas/wall/sample are booked by balance(), avoiding collectives here.
      totals[7].add(factor*pair.elementGas);
      totals[8].add(factor*pair.elementWall);
      totals[9].add(factor*pair.elementSample);
      const scalar held = ledger(pi, phaseChangeLedger::HELD);
      const scalar supplied = ledger.supplied(pi);
      totals[10].add(factor*held);
      totals[11].add(factor*supplied);
      totals[13].add(factor*(held - supplied));
      totals[14].add(factor*ledger.reference(pi));
    }
    OFstream& os = elementFiles_[ei];
    os << ledger.time();
    for (int k = 0; k < 12; ++k) os << ' ' << totals[k].value();
    os << ' ' << totals[10].value() - totals[11].value()
      << ' ' << totals[13].value() << ' ' << totals[14].value() << nl;
    os.flush();
  }
}

void LESTO::interfaceExchange::writeProfilesOfState(const bool timeWrite) {

  if (mock() || !profiles_) {
    return;
  }

  /*--------------------------------------------------------------------------
  The state is that of the ledger's clock: the end of the last step booked
  by balance() (or the start, before the first step); the same on all
  ranks.  Its profiles are written once, and once more when Time writes
  it after the profiles of its profile times, so that the reason (and the
  file) is that of a write time of controlDict at this step.
  --------------------------------------------------------------------------*/
  const scalar clock = ledger_->time();
  const bool sameState = profilesWritten_ && clock == profilesClock_;

  if (sameState && (profilesAtWrite_ || !timeWrite)) {
    return;
  }

  /* the profile times reached by the step of this state (balance()) */
  string times = profilesPending_;
  profilesPending_.clear();
  if (sameState) {
    times = profilesTimes_
      + (profilesTimes_.empty() || times.empty() ? "" : "; ") + times;
  }
  if (!timeWrite && times.empty()) {
    return;
  }

  const string reason = timeWrite
    ? string("write time") + (times.empty() ? "" : "; ") + times
    : times;

  profilesWritten_ = true;
  profilesClock_ = clock;
  profilesAtWrite_ = timeWrite;
  profilesTimes_ = times;

  const scalarField& V = mesh_.V();
  const scalarField& rhoI = rho_.primitiveField();

  forAll(pairs_, pairi) {

    const pairData& pair = pairs_[pairi];
    const scalarField& Y = species_[pair.gas].primitiveField();

    /* the masses of the inventories, cell by cell and element by element,
       with the same products as gasInventory() and wallInventory() */
    scalarField gasMass(Y.size());
    forAll(gasMass, celli) {
      gasMass[celli] = rhoI[celli]*Y[celli]*V[celli];
    }
    scalarField wallMass(pair.deposit.size());
    forAll(wallMass, e) {
      wallMass[e] = pair.deposit[e]*elementArea_[e];
    }

    /* the reservoirs of the inventory samples, per cell */
    scalarField sampleMass;
    if (pair.inventory) {
      sampleMass.resize(Y.size(), 0.0);
      for (label e = nWall_; e < elementCell_.size(); ++e) {
        sampleMass[elementCell_[e]] += reservoirMass(pair, e);
      }
    }
    FixedList<scalar, 3> inventory;
    inventory[0] = gasInventory(pair);
    inventory[1] = wallInventory(pair);
    inventory[2] = sampleInventory(pair);

    // Gas belonging to another pair's sample remains part of this pair's
    // wallPlusGas observable. Elemental profiles use the union of samples.
    boolList ownSamples;
    if (pairs_.size() > 1) {
      ownSamples.resize(mesh_.nCells(), false);
      for (const sampleData& sample : samples_) {
        if (sample.pair == pairi) {
          for (const label celli : sample.cells) ownSamples[celli] = true;
        }
      }
    }
    profiles_->write (
      mesh_.time().timeName(),
      clock,
      reason,
      speciesNames_[pair.gas],
      speciesNames_[pair.condensed],
      pair.molarMass,
      gasMass,
      wallMass,
      sampleMass,
      inventory,
      !pair.formula.empty(),
      pairs_.size() > 1 ? &ownSamples : nullptr
    );
  }
  if(respeciation_)respeciation_->writeWallProfiles(mesh_.time().timeName(),clock,reason);
  for (const word& element : elementNames_) {
    scalarField gas(V.size(), 0.0), wall(nWall_, 0.0), sample(V.size(), 0.0);
    FixedList<scalar, 3> inventory(Zero);
    for (const pairData& pair : pairs_) {
      const auto atom = pair.formula.find(element);
      if (atom == pair.formula.end()) continue;
      const scalar factor = atom->second/pair.molarMass;
      const scalarField& Y = species_[pair.gas].primitiveField();
      forAll(gas, c) gas[c] += factor*rhoI[c]*Y[c]*V[c];
      forAll(wall, e) wall[e] += factor*pair.deposit[e]*elementArea_[e];
      for (label e = nWall_; e < elementCell_.size(); ++e) {
        sample[elementCell_[e]] += factor*reservoirMass(pair, e);
      }
      inventory[0] += factor*gasInventory(pair);
      inventory[1] += factor*wallInventory(pair);
      inventory[2] += factor*sampleInventory(pair);
    }
    if (respeciation_) { respeciation_->addGasProfiles(element, gas, inventory[0]);
      respeciation_->addWallProfiles(element,wall,inventory[1]); }
    profiles_->write(mesh_.time().timeName(), clock, reason,
      word("element_" + element), element, 1.0, gas, wall, sample, inventory, true);
  }
}


/*------------------------------------------------------------------------------
End-of-run summary
------------------------------------------------------------------------------*/

void LESTO::interfaceExchange::end() const {

  if (mock()) {
    return;
  }

  forAll(pairs_, pairi) {
    const pairData& pair = pairs_[pairi];
    /* steps whose total was above the tolerance, named only if any */
    auto steps = [](const label n) {
      OStringStream os;
      if (n > 0) {
        os << "; " << n << " step" << (n > 1 ? "s" : "")
           << " above the tolerance in total";
      }
      return os.str();
    };
    Info<< "Phase change summary " << pair.solverName << ": "
      << pair.nUnconverged << " steps without outer convergence, "
      << pair.nGuardPassesRun << " guard passes, " << pair.nProjectedRun
      << " projected elements (" << pair.nProjectedSignificantRun
      << " significant" << steps(pair.nProjectedAboveRun).c_str() << "), "
      << pair.nComplementarityRun << " complementarity violations ("
      << pair.nComplementaritySignificantRun << " significant"
      << steps(pair.nComplementarityAboveRun).c_str() << ")";
    if (pair.finalSolver) {
      Info<< ", " << pair.nFinalChangesSignificantRun
        << " significant regime changes in final predictors";
    }
    Info<< nl;
  }

  /* equilibrium GEMS, mode local: the counts and the cost of the run */
  if (gems_) {
    gems_->summary();
  }
}

Foam::scalar LESTO::interfaceExchange::stepExchange(const Foam::label pi) const {
  compensatedSum sum;
  const pairData& pair = pairs_[pi];
  for (const scalar q : pair.flow) sum.add(q);
  return globalSum(sum);
}

void LESTO::interfaceExchange::writeChannelProfile(const word& timeName,scalar clock,
  const string& reason,const word& condensate,scalar mass,const scalarField& wall,scalar held) {
  if(!profiles_)return;
  scalarField gas(mesh_.nCells(),0),sample;
  FixedList<scalar,3> inventory({0,held,0});
  profiles_->write(timeName,clock,reason,word("condensate_"+condensate),condensate,mass,gas,wall,sample,inventory,true);
}

Foam::scalarField LESTO::interfaceExchange::wallDiffusionResistance(label si) const {
  scalarField resistance(nWall_,0);
  if(wallResistance_!="halfCell")return resistance;
  const surfaceScalarField& delta=mesh_.deltaCoeffs();
  forAll(resistance,e) {
    const label p=elementPatch_[e],f=elementFace_[e];
    const scalar diffusion=rhoD_[si].boundaryField()[p][f];
    resistance[e]=diffusion>0?1/(delta.boundaryField()[p][f]*diffusion):-1;
  }
  return resistance;
}
