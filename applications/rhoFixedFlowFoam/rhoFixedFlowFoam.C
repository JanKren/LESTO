/*------------------------------------------------------------------------------

rhoFixedFlowFoam -- fixed carrier flow, transient species transport, OF v2606.

Derived from the initialization and program structure of rhoSimpleFoam:
Copyright (C) 2011-2017 OpenFOAM Foundation.
Distributed under GNU GPL version 3 or (at your option) any later version.
https://www.gnu.org/licenses/gpl-3.0.html

PURPOSE

rhoFixedFlowFoam advances gaseous and solid species on an already converged
carrier flow.  The carrier fields U, p, T/e, rho and phi are read from disk and
remain fixed; no momentum, pressure, continuity or energy equation is solved.

Solid species may also be listed.  Their Y_<speciesName> fields are read,
registered and written.  At this development step a solid species is advanced
only by a volumetric source; it has no convection or diffusion.  The source is
applied only in cells adjacent to the WALL boundary patch.

The source terms are supplied by one thermochemistry routine which receives the
local temperature, pressure and all species values.  The present implementation
is only a mock-up of a future GEMS call and is evaluated at every physical time
step.

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

PROGRAM FLOW

  - create case/time database and mesh
  - read the frozen thermodynamic and carrier-flow fields
  - read the species list and their properties
  - create Y_i for every species and rhoD_i for gaseous species
  - identify cells adjacent to the WALL boundary patch

  At the present solid-species development step, a solid Y_i field is advanced
  only by ddt(rho,Y_i) = source.  No rhoD_i is created for a solid species.

  - while physical time advances {
      evaluate thermochemistry once in every cell
      suppress precipitation sources outside WALL-adjacent cells
      solve all gaseous species using the thermochemistry source terms
      solve all solid species using the updated source terms
      write requested fields
    }

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
#include "PbI2HeDiffusivity.H"
#include "evaluateThermochemistry.H"

#include <xGEMS/ChemicalEngine.hpp>   /* xGEMS engine */
#include <eigen3/Eigen/Dense>         /* library for vectors and matrices */

int main(int argc, char *argv[]) {

  argList::addNote (
    "Transient species transport on a frozen rhoSimpleFoam solution."
  );

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

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
  thermo.he().writeOpt() = IOobject::AUTO_WRITE;

  /*---------------------------------------------------------------------------
  Read the frozen carrier fields.  IOobject specifies the field name, the time
  directory, the mesh database, and the read/write policy.  MUST_READ means the
  saved field must exist; AUTO_WRITE lets runTime.write() write it later.

  The saved, pressure-corrected mass flux phi is required explicitly.
  Reconstructing phi from U would no longer reproduce the converged carrier
  solution.
  ---------------------------------------------------------------------------*/
  volScalarField rho (
    IOobject("rho",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
    mesh
  );

  volVectorField U (
    IOobject("U",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
    mesh
  );

  surfaceScalarField phi (
    IOobject("phi",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
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

  rho.oldTime();

  /*----------------------------------------------------------------------
  Evaluate viscosity and thermal conductivity once from the frozen thermo
  state.  They are volume-scalar fields and are marked AUTO_WRITE, but are
  never updated because thermo.correct() is never called.
  ----------------------------------------------------------------------*/
  volScalarField mu (
    IOobject("mu",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
    thermo.mu()
  );

  volScalarField kappa (
    IOobject("kappa",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
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
  Numerical controls shared by all gaseous species are kept in controlDict.
  The linear-solver and under-relaxation settings remain field-specific in
  fvSolution, while convection schemes remain field-specific in fvSchemes.
  -----------------------------------------------------------------------*/
  const dictionary& transportControls =
    runTime.controlDict().subDict("speciesTransport");

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
  pass the realized gas sink to the solid equation.
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

    if (state[speciesi] == "gas") {

      diffusionModel[speciesi] = speciesDict.getOrDefault<word> (
        "diffusivityModel", "constant"
      );

      if (diffusionModel[speciesi] != "constant"
        && diffusionModel[speciesi] != "PbI2He") {
        FatalErrorInFunction << "Unknown diffusivityModel "
          << diffusionModel[speciesi] << " for species " << speciesName
          << ". Choose constant or PbI2He." << exit(FatalError);
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
                     IOobject::AUTO_WRITE),
            rho*constantD
          )
        );

        Info<< "Species " << speciesName << ": field " << fieldName
          << ", diffusivityModel = constant, D = " << D << " m2/s" << nl;

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

    if (speciesDict.found("molarMass")) {
      molarMass[speciesi] = speciesDict.get<scalar>("molarMass");
      if (!std::isfinite(molarMass[speciesi]) || molarMass[speciesi] <= 0) {
        FatalErrorInFunction << "Species " << speciesName
          << " requires molarMass > 0 when specified. Got "
          << molarMass[speciesi] << exit(FatalError);
      }
    }
  }

  /*-------------------------------------------------------------------------
  Build a mask for the first layer of cells adjacent to the boundary region
  named WALL.  OpenFOAM already stores the owner cell of every boundary face,
  so no geometrical search is needed: faceCells() gives the adjacent internal
  cells directly.

  The mesh and carrier flow are fixed, so this connectivity is inspected only
  once at startup.  The thermochemistry routine remains unaware of the mesh;
  the CFD solver later suppresses solid-source terms in cells for which
  wallAdjacentCell is false.
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

  if (haveSolidSpecies) {

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

    /*----------------------------------------------------------------------
    Evaluate the complete thermochemistry state once in every cell at the
    beginning of the physical time step.  The routine receives T, p, all
    species names and the current species values, and returns one volumetric
    source for every species.  After the gas solve, the PbI2_s source is set
    to the realized PbI2_g sink, making their exchange conservative.

    The present routine is deliberately a simple mock-up.  Its interface is
    the part intended to survive when the implementation is eventually
    replaced by a GEMS call.
    -------------------------------------------------------------------------*/
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

        /*--------------------------------------------------------------------
        At this development step every non-zero thermochemistry source belongs
        to the PbI2 precipitation reaction.  The complete source vector is
        therefore suppressed outside the first cell layer adjacent to WALL, so
        the PbI2_g sink and PbI2_s source always act in exactly the same cells.
        --------------------------------------------------------------------*/
        if (haveSolidSpecies && !wallAdjacentCell[celli]) {
          sourceValue = 0;
        }

        speciesSource[speciesi][celli] = sourceValue;
      }
    }

    /*------------------------------------------------------------------------
    Advance gaseous species with the thermochemistry source included on the
    right-hand side of the transport equation.
    ------------------------------------------------------------------------*/
    forAll(species, speciesi) {
      if (state[speciesi] == "gas") {
        #include "solveGasSpecies.H"
      }
    }

    /*------------------------------------------------------------------------
    Print the volume-integrated thermochemistry source for every species.
    Equal and opposite PbI2_g/PbI2_s values provide a direct conservation check
    for the new gas-solid coupling.
    ------------------------------------------------------------------------*/
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

    /*------------------------------------------------------------------------
    Solid species have no convection or diffusion.  They are advanced only by
    the adjusted thermochemistry source retained in cells adjacent to WALL; the
    source has been set to zero in all other cells.
    ------------------------------------------------------------------------*/
    forAll(species, speciesi) {
      if (state[speciesi] == "solid") {
        #include "solveSolidSpecies.H"
      }
    }

    #include "checkSpeciesConservation.H"

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

    runTime.printExecutionTime(Info);
  }

  Info<< "End" << endl;
  return 0;
}
