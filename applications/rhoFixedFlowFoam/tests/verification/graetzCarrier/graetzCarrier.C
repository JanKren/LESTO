/*------------------------------------------------------------------------------

graetzCarrier.C

PURPOSE

Test helper of the M7 Graetz-Robin verification (tests/verification; not
part of the solver).  Writes the frozen carrier of a blockMesh wedge of a
pipe (axis along x, radius along y, the wedge planes z = +-y tan(theta/2))
into the start time directory:

  T    uniform;
  p    uniform;
  rho  p M/(R_gas T), the perfect gas of constant/thermophysicalProperties
       (uniform);
  U    (u(y), 0, 0) at the cell centres and on the patches, with the fully
       developed profile u(y) = 2 U_b (1 - y^2/R^2) of the radial coordinate
       y (the distance from the axis in the wedge's centre plane z = 0);
  phi  rho times the EXACT flux of U through every face:
         int_f u(y) dS_x = sum over the fan triangles (p0, pk, pk+1) of the
         face of their signed area projected on the y-z plane times the
         mean of u at their three edge midpoints.
       That quadrature is exact for the quadratic u, and since U is
       divergence-free the flux through a face depends only on its boundary
       loop, so the fluxes of every closed cell sum to zero up to
       round-off, whatever the triangulation.  The fluxes through the wedge
       planes and the wall (S_x = 0 exactly in exact geometry) are set to 0.

Why u(y) and not u(sqrt(y^2 + z^2)).  With the vertices at
(x, r, +-r tan(theta/2)) (as in OpenFOAM's wedge tutorials, e.g.
SandiaD_LTS), the wedge cells have the face areas 2 t r dx (radial faces),
t (r2^2 - r1^2) (axial faces), the volumes t (r2^2 - r1^2) dx and the
centroids y = (2/3)(r2^3 - r1^3)/(r2^2 - r1^2), with t = tan(theta/2):
exactly the metrics of an axisymmetric finite-volume scheme in the radius
y, times the constant t.  With u(y) the discrete problem is therefore a
consistent discretisation of the axisymmetric problem with the wall at
y = R for any wedge angle (no O(theta^2) geometry error).

Why not utilities/setFrozenCarrier (milestone M4).  It evaluates the
parabola at the distance from the axis, sqrt(y^2 + z^2) = y/cos(theta/2)
on the wedge planes: on these wedges that is the axisymmetric problem with
the wall at R cos(theta/2), a fixed bias of the decay rate of +6.9e-5 (Bi
1) and +2.3e-4 (Bi 50) at 2 degrees, larger than the GCI of
wallResistance halfCell (about 1e-6), and its wall points lie beyond the
radius R, which it refuses unless allowBackflow (plan section 33).

system/carrierDict

  R         2.4e-3;       // [m] wall radius: the WALL faces lie at y = R
  Ub        0.16;         // [m/s] bulk (mean) velocity of the section
  T         500;          // [K]
  p         101325;       // [Pa]
  molarMass 4.0026e-3;    // [kg/mol] of the carrier

T, p and U must exist in the start time as templates (they give the patch
types, normally calculated); rho and phi are created.  The helper prints
max |sum_f phi| per cell (absolute and relative to the largest face flux of
the cell), the inflow and outflow through every patch and the bulk
velocity they imply (inflow/(rho A_inlet)).

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wedgeFvPatch.H"

int main(int argc, char *argv[]) {

  argList::addNote("Frozen parabolic pipe-flow carrier of the Graetz wedge");

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

  const scalar Rw = carrierDict.get<scalar>("R");
  const scalar Ub = carrierDict.get<scalar>("Ub");
  const scalar T0 = carrierDict.get<scalar>("T");
  const scalar p0 = carrierDict.get<scalar>("p");
  const scalar M = carrierDict.get<scalar>("molarMass");
  const scalar Rgas = 8.314462618;
  const scalar rho0 = p0*M/(Rgas*T0);
  const label wallPatch = mesh.boundaryMesh().findPatchID("WALL");

  if (wallPatch < 0) {
    FatalErrorInFunction << "No patch WALL" << exit(FatalError);
  }

  /* the parabolic profile of the radius y */
  auto u = [&](const scalar y) -> scalar {
    return 2*Ub*(1 - (y*y)/(Rw*Rw));
  };

  /* exact int_f u dS_x over a face given by its points (fan from p0) */
  const pointField& points = mesh.points();
  auto faceFlux = [&](const face& f) -> scalar {
    scalar flux = 0;
    const point& a = points[f[0]];
    for (label k = 1; k + 1 < f.size(); ++k) {
      const point& b = points[f[k]];
      const point& c = points[f[k+1]];
      const scalar Sx =
        0.5*((b.y() - a.y())*(c.z() - a.z()) - (b.z() - a.z())*(c.y() - a.y()));
      const scalar mean =
        (u(0.5*(a.y() + b.y())) + u(0.5*(b.y() + c.y()))
       + u(0.5*(c.y() + a.y())))/3;
      flux += Sx*mean;
    }
    return flux;
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
    dimensionedScalar(dimDensity, rho0)
  );
  surfaceScalarField phi (
    IOobject("phi", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimMass/dimTime, 0)
  );

  const volVectorField& C = mesh.C();

  forAll(T, celli) {
    T[celli] = T0;
    p[celli] = p0;
    rho[celli] = rho0;
    U[celli] = vector(u(C[celli].y()), 0, 0);
  }

  forAll(T.boundaryField(), patchi) {
    const vectorField& Cf = mesh.C().boundaryField()[patchi];
    const bool wedge = isA<wedgeFvPatch>(mesh.boundary()[patchi]);
    if (wedge) {
      continue;   /* evaluated below from the cells */
    }
    forAll(Cf, facei) {
      T.boundaryFieldRef()[patchi][facei] = T0;
      p.boundaryFieldRef()[patchi][facei] = p0;
      rho.boundaryFieldRef()[patchi][facei] = rho0;
      U.boundaryFieldRef()[patchi][facei] =
        patchi == wallPatch ? vector::zero : vector(u(Cf[facei].y()), 0, 0);
    }
  }
  T.correctBoundaryConditions();
  p.correctBoundaryConditions();
  rho.correctBoundaryConditions();
  U.correctBoundaryConditions();

  phi.setOriented();
  const faceList& faces = mesh.faces();
  forAll(phi, facei) {
    phi[facei] = rho0*faceFlux(faces[facei]);
  }
  forAll(phi.boundaryField(), patchi) {
    const fvPatch& fvp = mesh.boundary()[patchi];
    const bool zero = patchi == wallPatch || isA<wedgeFvPatch>(fvp);
    const label start = fvp.patch().start();
    forAll(phi.boundaryField()[patchi], facei) {
      phi.boundaryFieldRef()[patchi][facei] =
        zero ? 0 : rho0*faceFlux(faces[start + facei]);
    }
  }

  /* continuity per cell, absolute and relative to the cell's largest flux */
  scalarField net(mesh.nCells(), 0);
  scalarField largest(mesh.nCells(), 0);
  const labelUList& own = mesh.owner();
  const labelUList& nei = mesh.neighbour();
  forAll(nei, facei) {
    net[own[facei]] += phi[facei];
    net[nei[facei]] -= phi[facei];
    largest[own[facei]] = max(largest[own[facei]], mag(phi[facei]));
    largest[nei[facei]] = max(largest[nei[facei]], mag(phi[facei]));
  }
  forAll(phi.boundaryField(), patchi) {
    const labelUList& fc = mesh.boundary()[patchi].faceCells();
    forAll(fc, facei) {
      const scalar f = phi.boundaryField()[patchi][facei];
      net[fc[facei]] += f;
      largest[fc[facei]] = max(largest[fc[facei]], mag(f));
    }
  }
  scalar maxRel = 0;
  forAll(net, celli) {
    if (largest[celli] > 0) {
      maxRel = max(maxRel, mag(net[celli])/largest[celli]);
    }
  }
  reduce(maxRel, maxOp<scalar>());

  Info().precision(17);
  Info<< "rho = " << rho0 << " kg/m3, T = " << T0 << " K, p = " << p0
    << " Pa, U_b = " << Ub << " m/s, R = " << Rw << " m" << nl
    << "max |sum_f phi| per cell = " << gMax(mag(net)) << " kg/s"
    << ", relative to the largest face flux of the cell "
    << maxRel << nl;

  forAll(phi.boundaryField(), patchi) {
    const scalar flow = gSum(phi.boundaryField()[patchi]);
    const scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
    Info<< "patch " << mesh.boundary()[patchi].name() << ": sum phi = "
      << flow << " kg/s, area " << area << " m2";
    if (area > 0 && mag(flow) > 0) {
      Info<< ", mean velocity " << mag(flow)/(rho0*area) << " m/s";
    }
    Info<< nl;
  }

  T.write();
  p.write();
  U.write();
  rho.write();
  phi.write();

  Info<< "End" << endl;
  return 0;
}
