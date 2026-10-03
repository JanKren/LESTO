/*------------------------------------------------------------------------------

makeCarrier.C

PURPOSE

Test helper of tests/phaseChange (not part of the solver).  Writes a frozen
carrier for the blockMesh test cases of the phase-change tests into the
start time directory:

  T    linear in x from Tin (x <= 0) to Tout (x >= L), clamped outside;
  p    uniform;
  rho  p M/(R T), the perfect gas of constant/thermophysicalProperties;
  U    G/rho in x (plug flow of the axial mass flux G [kg/(m2 s)]);
  phi  G*(S_f . x), exactly divergence-free on a mesh extruded along x.

system/carrierDict

  Tin     1000;         // [K]
  Tout    400;          // [K]
  L       0.1;          // [m]
  p       101325;       // [Pa]
  molarMass 4.0026e-3;  // [kg/mol]
  massFlux 0.01743;     // G [kg/(m2 s)]; 0 for a closed box
  Twall    450;         // optional [K]: T (and rho) on the patch WALL,
                        // e.g. to tell interfaceTemperature wall from cell

T, p and U must exist in the start time as templates (they give the patch
types, normally calculated); rho and phi are created.  The maximum of
|sum_f phi| per cell is printed as a continuity check.

------------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[]) {

  argList::addNote("Frozen test carrier for tests/phaseChange");

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  IOdictionary carrierDict (
    IOobject("carrierDict",
             runTime.system(),
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE)
  );

  const scalar Tin = carrierDict.get<scalar>("Tin");
  const scalar Tout = carrierDict.get<scalar>("Tout");
  const scalar L = carrierDict.get<scalar>("L");
  const scalar p0 = carrierDict.get<scalar>("p");
  const scalar M = carrierDict.get<scalar>("molarMass");
  const scalar G = carrierDict.get<scalar>("massFlux");
  const scalar R = 8.314462618;
  const bool haveTwall = carrierDict.found("Twall");
  const scalar Twall = carrierDict.getOrDefault<scalar>("Twall", 0);
  const label wallPatch = mesh.boundaryMesh().findPatchID("WALL");

  auto temperature = [&](const scalar x) -> scalar {
    return Tin + (Tout - Tin)*min(max(x/L, scalar(0)), scalar(1));
  };

  volScalarField T (
    IOobject("T", runTime.timeName(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
  );
  volScalarField p (
    IOobject("p", runTime.timeName(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
  );
  volVectorField U (
    IOobject("U", runTime.timeName(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
  );
  volScalarField rho (
    IOobject("rho", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimDensity, 1)
  );
  surfaceScalarField phi (
    IOobject("phi", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimMass/dimTime, 0)
  );

  const volVectorField& C = mesh.C();

  forAll(T, celli) {
    T[celli] = temperature(C[celli].x());
    p[celli] = p0;
    rho[celli] = p0*M/(R*T[celli]);
    U[celli] = vector(G/rho[celli], 0, 0);
  }

  forAll(T.boundaryField(), patchi) {
    const vectorField& Cf = mesh.C().boundaryField()[patchi];
    forAll(Cf, facei) {
      const scalar Tf =
        haveTwall && patchi == wallPatch ? Twall : temperature(Cf[facei].x());
      T.boundaryFieldRef()[patchi][facei] = Tf;
      p.boundaryFieldRef()[patchi][facei] = p0;
      rho.boundaryFieldRef()[patchi][facei] = p0*M/(R*Tf);
      U.boundaryFieldRef()[patchi][facei] = vector(G*R*Tf/(p0*M), 0, 0);
    }
  }

  phi.setOriented();
  const surfaceVectorField& Sf = mesh.Sf();
  forAll(phi, facei) {
    phi[facei] = G*Sf[facei].x();
  }
  forAll(phi.boundaryField(), patchi) {
    forAll(phi.boundaryField()[patchi], facei) {
      phi.boundaryFieldRef()[patchi][facei] =
        G*Sf.boundaryField()[patchi][facei].x();
    }
  }

  const scalarField netFlux(fvc::div(phi)().primitiveField()*mesh.V());
  Info<< "max |sum_f phi| per cell = " << gMax(mag(netFlux)) << " kg/s" << nl;

  T.write();
  p.write();
  U.write();
  rho.write();
  phi.write();

  Info<< "End" << endl;
  return 0;
}
