/*------------------------------------------------------------------------------

rhoFixedFlowFoam -- fixed carrier flow, transient species transport,
OpenCFD OpenFOAM v2412 and v2606.

Derived from the initialization and program structure of rhoSimpleFoam:
Copyright (C) 2011-2017 OpenFOAM Foundation.
Distributed under GNU GPL version 3 or (at your option) any later version.
https://www.gnu.org/licenses/gpl-3.0.html

PURPOSE

rhoFixedFlowFoam advances gaseous and solid species on an already converged
carrier flow.  The carrier fields U, p, T/e, rho and phi are read from disk and
remain fixed; no momentum, pressure, continuity or energy equation is solved.

Solid species may also be listed.  Their Y_<speciesName> fields are read,
registered and written.  What happens to them is selected by the model of
constant/thermochemistryProperties (interfaceExchange.H):

  - model mock (default, also without the file): the coupled mock of the
    main branch (tag 0.1, case run/011-coupled-species-sources).  One
    thermochemistry routine, a mock-up of a future GEMS call, receives the
    local temperature, pressure and all species values at the beginning of
    every physical time step and returns a source for every species: the
    precipitation PbI2_g -> PbI2_s, applied only in cells adjacent to the
    WALL boundary patch.  The gas equation treats the sink implicitly; a
    solid species has no convection or diffusion and receives the realised
    gas sink.  The integrated sources and the balance of every species are
    printed, and the source fields source_<species> are written.  This is
    the behaviour of the pristine build of tag 0.1, unchanged to the last
    bit;
  - model temperature or HKS: gas/condensate PAIRS exchange mass through
    the interface elements of the interaction layer (the faces of the
    listed wall patches): the temperature law of Lobresco et al. (2026,
    Eq. 5, irreversible deposition), or the Hertz-Knudsen-Schrage law with
    the equilibrium vapour pressure from a table (deposition and
    evaporation).  The exchange is an implicit, linearised source of the
    paired gas equation (solvePairedGasSpecies.H); the condensate is a
    deposit per wall element, and Y_<condensate> is derived output.
    Samples release a given amount of the gas over a finite time window,
    or (HKS) hold a condensate reservoir that evaporates and is removed at
    a given time.  A cumulative ledger closes the mass balance of every
    pair to round-off in every step (see PHASE CHANGE below).

For every gaseous species i, the solver reads the mass-fraction field
Y_<speciesName> and solves

  ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = S_i,

where rhoD_i = rho*D_i.  The molecular diffusivity D_i is selected separately
for each species in constant/speciesTransportProperties.  The currently
supported models are:

  - constant : a constant molecular diffusivity D_i in m2/s;
  - PbI2He   : the temperature- and pressure-dependent PbI2-in-He correlation.

The species themselves are listed at run time.  One species is therefore just
a list containing one entry; adding another passive gaseous species does not
require another solver.

If a species dictionary also contains molarMass in kg/mol, the solver writes
the derived molar concentration c_<speciesName> = rho*Y_i/M_i in mol/m3.

FROZEN CARRIER INPUT AND OUTPUT

The switch writeFrozenFields in the speciesTransport dictionary of
system/controlDict selects how the frozen carrier is read and written:

  - yes (default): the carrier fields are read from the start time and
    written at every write time together with the species, exactly as
    before the switch existed;
  - no: the carrier is read from the latest time directory, at or before the
    start time, that holds rho (normally 0), and it is never written.  The
    write directories then hold Y_*, c_*, uniform/ and, in model mock,
    source_*; a run can be restarted from a time directory with only Y_*
    fields.  T and p are registered before fluidThermo is constructed, so
    that fluidThermo adopts them instead of reading its own copies from the
    start time.  The log states the instance the carrier was read from.

In parallel the processor patches of T, p, rho and U are evaluated once from
the neighbour cells after reading (decomposePar writes a uniform processor
patch with about 10 digits only), so both sides of a processor face use the
same carrier values and the species transport is conservative across ranks.
This applies to every model, including mock: a decomposed run whose carrier
has such patches differs from a build without the evaluation in the last
digits (serial runs are unchanged).

PHASE CHANGE

For a paired gas species the transport equation receives the source
Su - Sp*Y with Su >= 0 and Sp >= 0 (fvm::Sp plus Su, never fvm::SuSp).  At
every outer corrector a local implicit predictor chooses the regime of every
interface element (IMPLICIT, CAPPED or NONE) from the transport-only matrix;
the paired equation is never relaxed, its final solve is converged to relTol
0, and a guard re-solves when the realised exchange violates a bound.  The
ledger books the release, the transport through every non-coupled patch and
the per-cell solver defect (of the transport-only matrix plus the realised
exchange: matrixDefect.H), so its closure is round-off for any linear-
solver tolerance.  The sign convention is: positive = evaporation into the
gas.  See interfaceExchange.H, phaseChangeKinetics.H, phaseChangeLedger.H
and readme.md.  The equilibrium vapour pressure of the HKS law comes from a
table, or from GEMS3K (equilibrium GEMS, a build with the bridge; mode
frozen at start-up, or mode local per step and element: gemsEquilibrium.H).

Outputs of the exchange: the deposit per effective element area mw_ and
per wall area mDep_ of every condensate, the reservoir Cs_ of the inventory
samples, optionally the realised exchange exch_ of every paired gas, the
ledger in postProcessing/phaseChangeBalance, and, with a profiles
dictionary, the axial profiles in mol/m and T_dep in postProcessing/
axialProfiles at the write times and at the profile times (axialProfiles.H).

PROGRAM FLOW

  - create case/time database and mesh
  - select the frozen-carrier instance and write policy (writeFrozenFields)
  - read the frozen thermodynamic and carrier-flow fields
  - read the species list and their properties
  - create Y_i for every species and rhoD_i for gaseous species
  - create the phase-change object (nothing at all for model mock)
  - model mock: identify cells adjacent to the WALL boundary patch (other
    models use the interface patches of thermochemistryProperties instead)

  In model mock a solid Y_i field is advanced only by
  ddt(rho,Y_i) = source.  No rhoD_i is created for a solid species.

  - while physical time advances {
      phase change: the exact step on the clock of the ledger, the removal
        of inventory samples, bounds, equilibrium GEMS in mode local: p_eq
        of every element from GEMS3K (every updateInterval steps), release
        rates and a check of the final solver of the step
      model mock: evaluate thermochemistry once in every cell and suppress
        the precipitation sources outside WALL-adjacent cells
      solve all gaseous species: paired ones with the exchange, the others
        with the thermochemistry source terms (zero outside model mock)
      model mock {
        print the integrated thermochemistry source of every species
        solve all solid species using the updated source terms
        report the balance of every species and of their total
      }
      otherwise: mass ledger and closure, derived Y_<condensate>, the
        deposit and reservoir fields (mw_, mDep_, Cs_), the optional
        exchange field (exch_) and the per-element state of uniform/
      write requested fields (the frozen carrier only if writeFrozenFields;
        the source fields in model mock only; the phase-change ledger and
        element state of uniform/ with them)
      write times: molar concentrations c_<species>
      phase change: axial profiles and T_dep at the write times and the
        profile times (postProcessing/axialProfiles only)
    }
    phase change: end-of-run summary

BUILD

Allwmake builds the solver with wmake.  The GEMS3K bridge is an optional
build dependency (LESTO_GEMSBRIDGE, see Make/options and gemsEquilibrium.H);
rhoFixedFlowFoam -help reports whether it is compiled in.  So is the
xGEMS start-up call of the main branch, an optional build dependency as
well (LESTO_XGEMS, an xGEMS installation prefix; Make/options and
Make/xgemsConfig.sh; see the start of main() below).

OPENFOAM STYLE

The program follows the usual OpenFOAM solver structure.  Typical OF objects
used in this program include:
  - runTime manages physical time and output;
  - mesh is the finite-volume mesh and object database;
  - IOobject describes how fields are read and written;
  - IOdictionary reads OpenFOAM dictionaries from disk;
  - volScalarField, volVectorField and surfaceScalarField store cell/face data;
  - PtrList stores a run-time-sized list of OpenFOAM objects;
  - fluidThermo provides the thermodynamic state;
  - Info/WarningInFunction/FatalErrorInFunction are OpenFOAM's standard
    reporting streams.

The included .H files are code fragments inserted directly at the #include
locations, following common OpenFOAM solver practice.

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "processorFvPatch.H"
#include "PbI2HeDiffusivity.H"
#include "ChapmanEnskogDiffusivity.H"
#include "chemicalFormula.H"
#include "evaluateThermochemistry.H"
#include "gemsEquilibrium.H"
#include "interfaceExchange.H"
#include "gasRespeciation.H"

/*------------------------------------------------------------------------------
Optional xGEMS start-up call of the main branch (71bf192, see the start of
main()).  It is compiled only in a build with LESTO_XGEMS, an xGEMS
installation prefix: Make/xgemsConfig.sh then defines LESTO_HAVE_XGEMS in
the generated header xgemsConfig.H, and Make/options adds the include and
library flags of the prefix.  Without LESTO_XGEMS neither the xGEMS nor
the Eigen headers are needed.
------------------------------------------------------------------------------*/
#include "xgemsConfig.H"

#ifdef LESTO_HAVE_XGEMS
#include <xGEMS/ChemicalEngine.hpp>   /* xGEMS engine */
#include <eigen3/Eigen/Dense>         /* library for vectors and matrices */
#endif

int main(int argc, char *argv[]) {

  argList::addNote (
    "Transient species transport on a frozen rhoSimpleFoam solution."
  );

  /*--------------------------------------------------------------------------
  Build information only; notes are shown by -help and never in a run log.
  --------------------------------------------------------------------------*/
  argList::addNote(LESTO::gemsEquilibrium::buildState());
#ifdef LESTO_HAVE_XGEMS
  argList::addNote("xGEMS start-up call: compiled in from " LESTO_XGEMS_DIR);
#else
  argList::addNote (
    "xGEMS start-up call: not compiled in (set LESTO_XGEMS and rerun"
    " Allwmake to enable it)"
  );
#endif

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

#ifdef LESTO_HAVE_XGEMS
  /*--------------------------------------------------------------------------
  The first call of xGEMS (main branch, 71bf192), kept as written there and
  compiled only with LESTO_XGEMS (see above).  It equilibrates a cement demo
  system of the xGEMS sources, at the path of that author's machine, 1 K
  above the system's own temperature, and prints the status.  It tests the
  link only: nothing below uses its result.  Elsewhere, point the path at
  the demo system of the local xGEMS sources.
  --------------------------------------------------------------------------*/

  /* Call xGEMS */
  {
    xGEMS::ChemicalEngine gems;
    Info<< "xGEMS ChemicalEngine constructed successfully" << nl;

    gems.initialize(
      "/home/niceno/Development/GEMS-Related/xgems/demos/resources/"
      "CemGEMS-keyvalue/CemHyds-dat.lst"
    );

    const Eigen::VectorXd bulk = gems.elementAmounts();

    const int status = gems.equilibrate(
      gems.temperature() + 1.0,
      gems.pressure(),
      bulk
    );

    Info<< "xGEMS status: " << status
        << ", converged: " << gems.converged()
        << ", iterations: " << gems.numIterations() << nl;
  }
#endif

  /*--------------------------------------------------------------------------
  Numerical controls shared by all gaseous species are kept in the
  speciesTransport dictionary of controlDict.  It is read here, before the
  carrier, because writeFrozenFields decides where the carrier is read from.
  --------------------------------------------------------------------------*/
  const dictionary& transportControls =
    runTime.controlDict().subDict("speciesTransport");

  /*--------------------------------------------------------------------------
  Frozen-carrier input and output (writeFrozenFields, default yes).

  frozenWrite is the write policy of every field that never changes during
  the run: the carrier (U, p, T, e, rho, phi), mu, kappa and the diffusion
  coefficients.  frozenInstance is the time directory the carrier is read
  from:

    - yes: the start time, and the fields are written at every write time,
      exactly as before this switch existed;
    - no:  the latest time directory, searching back from the start time,
      that holds rho.  Nothing frozen is written, so a later restart finds
      the carrier in the same directory again.

  fluidThermo would read T and p from the start time.  With "no" they are
  therefore read here from frozenInstance, NO_WRITE, and handed to the
  object registry with store().  fluidThermo finds them there and adopts
  them instead of reading its own copies.
  --------------------------------------------------------------------------*/
  const bool writeFrozenFields =
    transportControls.getOrDefault<bool>("writeFrozenFields", true);

  const IOobject::writeOption frozenWrite =
    writeFrozenFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;

  const word frozenInstance = writeFrozenFields
    ? runTime.timeName()
    : runTime.findInstance(mesh.dbDir(), "rho");

  if (!writeFrozenFields) {

    for (const word& frozenName : wordList({"p", "T"})) {
      regIOobject::store (
        new volScalarField (
          IOobject(frozenName,
                   frozenInstance,
                   mesh,
                   IOobject::MUST_READ,
                   IOobject::NO_WRITE),
          mesh
        )
      );
    }

    Info<< "Frozen carrier read from " << frozenInstance
      << " and not written (writeFrozenFields no)" << nl;

    /*------------------------------------------------------------------------
    findInstance returns the LATEST directory at or before the start time
    that holds rho.  If an earlier directory holds rho as well, the carrier
    found is most likely a copy written by an earlier run with
    writeFrozenFields yes (rounded to writePrecision in ascii), and the
    original carrier is not used.
    ------------------------------------------------------------------------*/
    for (const instant& t : runTime.times()) {
      if (t.name() == runTime.constant() || t.name() == frozenInstance) {
        if (t.name() == frozenInstance) {
          break;
        }
        continue;
      }
      IOobject earlier("rho", t.name(), mesh);
      if (returnReduceOr(earlier.typeHeaderOk<volScalarField>(false, false))) {
        WarningInFunction
          << "The frozen carrier is read from " << frozenInstance
          << ", but the earlier time directory " << t.name()
          << " also holds rho." << nl
          << "    " << frozenInstance << " may be a copy written by an "
          << "earlier run with writeFrozenFields yes (rounded to "
          << "writePrecision in ascii); the carrier in " << t.name()
          << " is not used." << nl << endl;
        break;
      }
    }
  }

  /*--------------------------------------------------------------------------
  fluidThermo is OpenFOAM's thermodynamic model object.  The New(mesh) factory
  reads thermophysicalProperties and constructs the model selected there.
  autoPtr is OpenFOAM's owning smart pointer.

  Please note that "he" is OpenFOAM's enthalpy (h) or energy (e) variable.
  --------------------------------------------------------------------------*/
  Info << "Reading frozen thermophysical state" << nl;
  autoPtr     <fluidThermo> pThermo(fluidThermo::New(mesh));
  fluidThermo & thermo = pThermo();

  thermo.validate(args.executable(), "e");
  thermo.he().writeOpt(frozenWrite);

  /*---------------------------------------------------------------------------
  Read the frozen carrier fields.  IOobject specifies the field name, the time
  directory, the mesh database, and the read/write policy.  MUST_READ means the
  saved field must exist; frozenWrite (AUTO_WRITE or NO_WRITE, see above)
  decides whether runTime.write() writes it later.

  The saved, pressure-corrected mass flux phi is required explicitly.
  Reconstructing phi from U would no longer reproduce the converged carrier
  solution.
  ---------------------------------------------------------------------------*/
  volScalarField rho (
    IOobject("rho",
              frozenInstance,
              mesh,
              IOobject::MUST_READ,
              frozenWrite),
    mesh
  );

  volVectorField U (
    IOobject("U",
              frozenInstance,
              mesh,
              IOobject::MUST_READ,
              frozenWrite),
    mesh
  );

  surfaceScalarField phi (
    IOobject("phi",
              frozenInstance,
              mesh,
              IOobject::MUST_READ,
              frozenWrite),
    mesh
  );

  if (rho.dimensions() != dimDensity || phi.dimensions() != dimMass/dimTime) {
    FatalErrorInFunction
      << "Require rho in kg/m3 and saved phi in kg/s."
      << exit(FatalError);
  }

  if (gMin(rho.primitiveField()) <= 0) {
    FatalErrorInFunction << "Density must be positive."
      << exit(FatalError);
  }

  /*--------------------------------------------------------------------------
  Processor patches of the frozen carrier.  The value of a processor patch
  must be the neighbour's cell value, but decomposePar writes a patch whose
  values are all equal as "uniform <value>" with about 10 significant
  digits, also in a binary file, while the neighbour's cell values are
  written exactly.  The two sides of a processor face would then see
  different T, p and rho (5e-11 relative in the plug-flow channel of
  tests/phaseChange), hence a different rhoD in the diffusion term, and
  the diffusive flux leaving one rank would not be the flux entering the
  other.  Evaluating the processor patches once from the neighbour cells
  makes both sides exact.  Only processor patches are touched, so a serial
  run is unchanged; a decomposed run changes in the last digits for every
  model, mock included (plan section 15.7); phi is a face field and is not
  affected.
  --------------------------------------------------------------------------*/
  thermo.T().boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  thermo.p().boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  rho.boundaryFieldRef().evaluateCoupled<processorFvPatch>();
  U.boundaryFieldRef().evaluateCoupled<processorFvPatch>();

  rho.oldTime();

  /*----------------------------------------------------------------------
  Evaluate viscosity and thermal conductivity once from the frozen thermo
  state.  They are volume-scalar fields written with the frozen carrier
  (frozenWrite), but are never updated because thermo.correct() is never
  called.
  ----------------------------------------------------------------------*/
  volScalarField mu (
    IOobject("mu",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              frozenWrite),
    thermo.mu()
  );

  volScalarField kappa (
    IOobject("kappa",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              frozenWrite),
    thermo.kappa()
  );

  /*--------------------------------------------------------------------------
  speciesTransportProperties contains the species list and one sub-dictionary
  for each species.  Gaseous species select a diffusivity model and transport
  properties.  Solid species appear in the same list but do not select or use
  a diffusivity model.

  Their source terms are supplied later by evaluateThermochemistry(), which is
  deliberately kept separate from the transport solver.
  --------------------------------------------------------------------------*/
  IOdictionary speciesProperties (
    IOobject("speciesTransportProperties",
              runTime.constant(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE)
  );

  const wordList speciesNames(speciesProperties.lookup("species"));

  if (speciesNames.size() == 0) {
    FatalErrorInFunction << "The species list must contain at least one species."
      << exit(FatalError);
  }

  /*-----------------------------------------------------------------------
  Numerical controls shared by all gaseous species are kept in controlDict
  (transportControls, read above).  The linear-solver and under-relaxation
  settings remain field-specific in fvSolution, while convection schemes
  remain field-specific in fvSchemes.
  -----------------------------------------------------------------------*/
  const label nCorr = transportControls.getOrDefault<label>("nCorr", 0);
  const scalar tolerance =
    transportControls.getOrDefault<scalar>("tolerance", 1);

  if (nCorr < 0 || !std::isfinite(tolerance) || tolerance < 0) {
    FatalErrorInFunction << "Require nCorr >= 0 and finite tolerance >= 0."
      << exit(FatalError);
  }

  /*-------------------------------------------------------------------------
  PtrList is OpenFOAM's owning list of objects.  The number of fields is not
  known when the solver is compiled; it is determined by the species list in
  speciesTransportProperties when the case starts.

  species[i]   stores Y_<speciesName>.
  rhoD[i]      stores rho*D_i, the coefficient used in the diffusion term.
  molecularD[i] is used only by models which create a spatial D_i field, such
                as PbI2He; unused PtrList entries remain null.
  molarMass[i] is negative when no molar-concentration output is requested.
  -------------------------------------------------------------------------*/

  /*
  For a solid species only species[i] is populated.  Its rhoD and molecularD
  entries remain null because the solid equation has no diffusion term.

  speciesSource[i] stores the thermochemistry source for species i in
  kg/(m3 s).  The source fields are evaluated at the beginning of each
  physical time step.  The PbI2 sources are adjusted after the gas solve to
  pass the realized gas sink to the solid equation.  This is model mock;
  with a phase-change model the source fields stay zero and are not written
  (see below the creation of the phase-change object).
  */
  PtrList <volScalarField> species        (speciesNames.size());
  PtrList <volScalarField> rhoD           (speciesNames.size());
  PtrList <volScalarField> molecularD     (speciesNames.size());
  PtrList <volScalarField> speciesSource  (speciesNames.size());
  List <scalar>            molarMass      (speciesNames.size(), -1.0);
  List <word>              diffusionModel (speciesNames.size());
  List <word>              state          (speciesNames.size());

  forAll(speciesNames, speciesi) {

    const word &       speciesName = speciesNames[speciesi];
    const word         fieldName("Y_" + speciesName);
    const word         rhoDName("rhoD_" + speciesName);
    const dictionary & speciesDict = speciesProperties.subDict(speciesName);
    state[speciesi]                = speciesDict.get<word>("state");

    /*----------------------
    Check for allowed states
    ----------------------*/
    if (state[speciesi] != "gas" && state[speciesi] != "solid") {
      FatalErrorInFunction
        << "Unknown state " << state[speciesi]
        << " for species " << speciesName
        << ". Choose gas or solid."
        << exit(FatalError);
    }

    /*---------------------------------------------------------------------
    A solid species has no diffusivity.  Its source is not configured here;
    evaluateThermochemistry() determines the source from the complete local
    thermochemical state during the physical-time loop.
    ---------------------------------------------------------------------*/
    /*
    Read and register the species field for both allowed states.  AUTO_WRITE
    ensures that an unchanged solid field is still written at output times.
    */
    species.set (
      speciesi,
      new volScalarField (
        IOobject(fieldName,
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE),
        mesh
      )
    );

    if (species[speciesi].dimensions() != dimless) {
      FatalErrorInFunction << fieldName << " must be dimensionless."
        << exit(FatalError);
    }

    /*--------------------------------------
    Create a source field for every species.
    --------------------------------------*/
    speciesSource.set (
      speciesi,
      new volScalarField (
        IOobject("source_" + speciesName,
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("zeroSource", dimDensity/dimTime, 0)
      )
    );

    if (speciesDict.found("molarMass")) {
      molarMass[speciesi] = speciesDict.get<scalar>("molarMass");
      if (!std::isfinite(molarMass[speciesi]) || molarMass[speciesi] <= 0) {
        FatalErrorInFunction << "Species " << speciesName
          << " requires molarMass > 0 when specified. Got "
          << molarMass[speciesi] << exit(FatalError);
      }
    }
    LESTO::readFormula(speciesDict, molarMass[speciesi]);

    if (state[speciesi] == "gas") {

      diffusionModel[speciesi] = speciesDict.getOrDefault<word> (
        "diffusivityModel", "constant"
      );

      if (diffusionModel[speciesi] != "constant"
        && diffusionModel[speciesi] != "PbI2He"
        && diffusionModel[speciesi] != "ChapmanEnskog") {
        FatalErrorInFunction << "Unknown diffusivityModel "
          << diffusionModel[speciesi] << " for species " << speciesName
          << ". Choose constant, PbI2He or ChapmanEnskog." << exit(FatalError);
      }

      if (diffusionModel[speciesi] == "constant") {

        const scalar D = speciesDict.get<scalar>("D");
        if (!std::isfinite(D) || D < 0) {
          FatalErrorInFunction << "Species " << speciesName
            << " requires finite D >= 0. Got " << D << exit(FatalError);
        }

        /*--------------------------------------------------------------------
        dimensionedScalar carries both the numerical value and the OpenFOAM
        dimensions.  D is the molecular diffusivity [m2/s]; multiplying by rho
        creates the equation coefficient rhoD [kg/(m s)].
        --------------------------------------------------------------------*/
        const dimensionedScalar constantD (
          "D", dimViscosity, D
        );

        rhoD.set (
          speciesi,
          new volScalarField (
            IOobject(rhoDName,
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     frozenWrite),
            rho*constantD
          )
        );

        Info<< "Species " << speciesName << ": field " << fieldName
          << ", diffusivityModel = constant, D = " << D << " m2/s" << nl;

      } else if (diffusionModel[speciesi] == "ChapmanEnskog") {

        #include "createChapmanEnskogDiffusivity.H"

      } else {

        if (speciesName != "PbI2_g") {
          FatalErrorInFunction << "The PbI2He model is specific to PbI2_g, not "
            << speciesName << "." << exit(FatalError);
        }

        if (speciesDict.found("D")) {
          FatalErrorInFunction << "Species " << speciesName
            << " uses diffusivityModel PbI2He: remove constant D."
            << exit(FatalError);
        }

        /*---------------------------------------------------------------------
        Insert the setup code which evaluates the same PbI2-He correlation used
        by the original single-species case 009.  It creates D_PbI2_g(T,p) and
        rhoD_PbI2_g = rho*D_PbI2_g once from the frozen carrier fields.
        ---------------------------------------------------------------------*/
        #include "createVariableDiffusivity.H"
      }

    }

    /*----------------------------------------------------------------------
    Store the starting species field as the old-time field used by the Euler
    transient term.  This is required by both the gas transport equation and
    the solid accumulation equation.  For PbI2_g this remains the same
    operation as in the original single-species case 009.
    ----------------------------------------------------------------------*/
    species[speciesi].oldTime();


  }

  /*-------------------------------------------------------------------------
  Gas/condensate exchange (constant/thermochemistryProperties).  For model
  mock, or without the file, the object is inert and prints nothing.  It is
  created before the WALL mask below, which only the mock path uses.
  -------------------------------------------------------------------------*/
  #include "createPhaseChange.H"
  LESTO::gasRespeciation respeciation(mesh, thermo, rho, species,
    speciesNames, state, molarMass, phaseChange);

  /*-------------------------------------------------------------------------
  The source fields source_<species> belong to the thermochemistry of model
  mock, which writes them (AUTO_WRITE above).  With a phase-change model
  evaluateThermochemistry() is not called, the fields stay zero, and they
  are not written: the outputs of the phase-change models are unchanged by
  them.
  -------------------------------------------------------------------------*/
  if (!phaseChange.mock()) {
    forAll(speciesSource, speciesi) {
      speciesSource[speciesi].writeOpt(IOobject::NO_WRITE);
    }
  }

  /*-------------------------------------------------------------------------
  model mock: build a mask for the first layer of cells adjacent to the
  boundary region named WALL.  OpenFOAM already stores the owner cell of
  every boundary face, so no geometrical search is needed: faceCells() gives
  the adjacent internal cells directly.

  The mesh and carrier flow are fixed, so this connectivity is inspected only
  once at startup.  The thermochemistry routine remains unaware of the mesh;
  when solid species exist, the CFD solver later suppresses the complete
  thermochemistry source vector in cells for which wallAdjacentCell is
  false.

  With a phase-change model the mask is neither built nor printed: the
  condensate lives on the interface patches named in constant/
  thermochemistryProperties (interface { patches (...); }), which need not
  be called WALL.
  -------------------------------------------------------------------------*/
  bool haveSolidSpecies = false;
  forAll(state, speciesi) {
    if (state[speciesi] == "solid") {
      haveSolidSpecies = true;
      break;
    }
  }

  /*-----------------------
  Find wall-adjencent cells
  -----------------------*/
  boolList wallAdjacentCell(mesh.nCells(), false);

  if (phaseChange.mock() && haveSolidSpecies) {

    const label wallPatchi = mesh.boundaryMesh().findPatchID("WALL");

    if (wallPatchi < 0) {
      FatalErrorInFunction
        << "Solid species are configured, but boundary patch WALL was not found."
        << exit(FatalError);
    }

    const labelUList& wallCells =
      mesh.boundary()[wallPatchi].faceCells();

    forAll(wallCells, facei) {
      wallAdjacentCell[wallCells[facei]] = true;
    }

    label nWallAdjacentCells = 0;
    forAll(wallAdjacentCell, celli) {
      if (wallAdjacentCell[celli]) {
        ++nWallAdjacentCells;
      }
    }
    reduce(nWallAdjacentCells, sumOp<label>());

    Info<< "Cells adjacent to WALL: " << nWallAdjacentCells << nl;
  }

  Info<< "Frozen fields: U, p, e, T, rho, phi, mu, kappa" << nl
    << "Configured species: " << speciesNames << nl
    << "No SIMPLE/PIMPLE, pressure, momentum or energy solve." << nl
    << "Species controls are read at startup; restart to change them." << nl
    << "Starting physical time loop" << endl;

  while (runTime.loop()) {

    /* Lists for checking the species mass conservation */
    List <scalar> speciesBoundaryInRate (species.size(), 0);
    List <scalar> speciesBoundaryOutRate(species.size(), 0);
    List <scalar> realizedSourceRate    (species.size(), 0);

    Info<< "Time = " << runTime.timeName() << nl << endl;

    /*-------------------------------------------------------------------------
    Phase change: release rates, element bounds and counters of the step
    (nothing for model mock).
    -------------------------------------------------------------------------*/
    phaseChange.beginStep();
    respeciation.beginStep();

    if (phaseChange.mock()) {

      /*----------------------------------------------------------------------
      Evaluate the complete thermochemistry state once in every cell at the
      beginning of the physical time step.  The routine receives T, p, all
      species names and the current species values, and returns one
      volumetric source for every species.  After the gas solve, the PbI2_s
      source is set to the realized PbI2_g sink, making their exchange
      conservative.

      The present routine is deliberately a simple mock-up.  Its interface
      is the part intended to survive when the implementation is eventually
      replaced by a GEMS call.  It belongs to model mock: with a phase-change
      model it is not called, and every source field stays zero.
      ----------------------------------------------------------------------*/
      const volScalarField& T = thermo.T();
      const volScalarField& p = thermo.p();

      scalarField localSpecies(speciesNames.size(), 0.0);
      scalarField localSources(speciesNames.size(), 0.0);

      forAll(T, celli) {

        forAll(speciesNames, speciesi) {
          localSpecies[speciesi] = species[speciesi][celli];
        }

        evaluateThermochemistry (
          T[celli],
          p[celli],
          speciesNames,
          localSpecies,
          localSources
        );

        forAll(speciesNames, speciesi) {

          scalar sourceValue = localSources[speciesi];

          /*------------------------------------------------------------------
          At this development step every non-zero thermochemistry source
          belongs to the PbI2 precipitation reaction.  The complete source
          vector is therefore suppressed outside the first cell layer
          adjacent to WALL, so the PbI2_g sink and PbI2_s source always act
          in exactly the same cells.
          ------------------------------------------------------------------*/
          if (haveSolidSpecies && !wallAdjacentCell[celli]) {
            sourceValue = 0;
          }

          speciesSource[speciesi][celli] = sourceValue;
        }
      }
    }

    /*-------------------------------------------------------------------------
    The carrier flow remains frozen.  Advance all gaseous species.  A gas
    paired with a condensate (not in model mock) is solved with the exchange
    source (solvePairedGasSpecies.H).  Every other gas species is solved by
    solveGasSpecies.H with the thermochemistry source included on the
    right-hand side of the transport equation.  Outside model mock that
    source is zero, which leaves the solution of such a species, e.g. a
    passive tracer, unchanged to the last bit.
    -------------------------------------------------------------------------*/
    forAll(species, speciesi) {
      if (state[speciesi] == "gas") {
        if (phaseChange.pairOfGas(speciesi) >= 0) {
          #include "solvePairedGasSpecies.H"
        } else if (respeciation.wallGas(speciesi)) {
          #include "solveWallChannelGasSpecies.H"
        } else if (respeciation.gasOnly(speciesi)) {
          #include "solveLedgeredGasSpecies.H"
        } else {
          #include "solveGasSpecies.H"
        }
      }
    }

    respeciation.applyWallProducts();
    respeciation.apply();

    /*-------------------------------------------------------------------------
    model mock: the integrated sources, the solid solve and the balance of
    every species.  Otherwise the condensate is updated inside the paired gas
    solve and the ledger closes the balance of every pair here, before
    runTime.write(), so that the derived Y_<condensate> is written
    consistently.  The per-species balance of model mock is not reported
    with a phase-change model: the ledger books the balance of every pair.
    -------------------------------------------------------------------------*/
    if (phaseChange.mock()) {

      /*----------------------------------------------------------------------
      Print the volume-integrated thermochemistry source for every species.
      Equal and opposite PbI2_g/PbI2_s values provide a direct conservation
      check for the new gas-solid coupling.
      ----------------------------------------------------------------------*/
      scalar netThermochemistrySource = 0;

      forAll(speciesSource, speciesi) {
        const scalar integratedSource =
          fvc::domainIntegrate(speciesSource[speciesi]).value();

        Info<< "Integrated source " << speciesNames[speciesi] << " = "
          << integratedSource << " kg/s" << nl;

        netThermochemistrySource += integratedSource;
      }

      Info<< "Net thermochemistry source = " << netThermochemistrySource
        << " kg/s" << nl;

      /*----------------------------------------------------------------------
      Solid species have no convection or diffusion.  They are advanced only
      by the adjusted thermochemistry source retained in cells adjacent to
      WALL; the source has been set to zero in all other cells.
      ----------------------------------------------------------------------*/
      forAll(species, speciesi) {
        if (state[speciesi] == "solid") {
          #include "solveSolidSpecies.H"
        }
      }

      #include "checkSpeciesConservation.H"

    } else {
      phaseChange.balance();
    }

    runTime.write();

    if (runTime.writeTime()) {

      /*----------------------------------------------------------------------
      Molar concentration is a derived output only.  It is written for species
      that define molarMass in speciesTransportProperties.
      ----------------------------------------------------------------------*/
      forAll(species, speciesi) {
        if (molarMass[speciesi] > 0) {
          const dimensionedScalar Mi (
            "Mi", dimensionSet(1, 0, 0, 0, -1, 0, 0), molarMass[speciesi]
          );

          const word concentrationName("c_" + speciesNames[speciesi]);

          volScalarField concentration (
            IOobject(concentrationName,
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE),
            rho*species[speciesi]/Mi
          );

          concentration.write();
          Info<< "Wrote " << concentration.name() << " [mol/m3]" << nl;
        }
      }
    }

    /*-------------------------------------------------------------------------
    Phase change: axial profiles and T_dep at the profile times of
    constant/thermochemistryProperties reached by this step, written to
    postProcessing/axialProfiles only (nothing for model mock or without a
    profiles dictionary).  At every write of Time, runTime.write() above or
    a write of a function object before the next step, Time itself writes
    them (interfaceExchange.H, profileWriter).
    -------------------------------------------------------------------------*/
    phaseChange.writeProfiles();

    runTime.printExecutionTime(Info);
  }

  phaseChange.end();

  Info<< "End" << endl;
  return 0;
}
