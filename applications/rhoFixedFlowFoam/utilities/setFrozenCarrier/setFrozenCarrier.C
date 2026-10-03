/*------------------------------------------------------------------------------

setFrozenCarrier.C

PURPOSE

Writes the frozen carrier of rhoFixedFlowFoam in 'paper mode' (plan sections
1, 5 and 11, milestone M4): a fully developed, non-expanding laminar pipe
flow of constant density, as T-Flows prescribes it for the code-to-code
cases of Lobresco et al. (2026), on any mesh of an axisymmetric pipe: the
gmsh pipe of run/001 or an OpenFOAM wedge.  It is a utility, not solver
logic: the solver only reads the fields it writes (writeFrozenFields,
findInstance of rho).

  rho  the constant 'density' (1 kg/m3 in paper mode: T-Flows sets
       Flow % density = 1 in Beginning_Of_Iteration);
  U    2 U_b (1 - r^2/R^2) along the axis, r the distance from the axis, at
       the cell centres and the face centres (0 on the wall patches, and 0
       where r > R);
  phi  rho times the EXACT flux of that velocity through every face: the
       face is split into the triangles of a fan from its first point, and
       on every triangle the three-edge-midpoint rule, exact for the
       quadratic profile, integrates u (axis . n) dA.  The triangles of the
       faces of a cell close its surface, and the profile has no divergence
       (it does not vary along the axis), so sum_f phi of every cell is 0 to
       round-off on ANY mesh, not only on a mesh extruded along the axis
       (tests/paperMode: a wedge whose inner points are moved radially and
       axially, where a face-centre flux misses continuity by 1e-5);
  p    the constant 'pressure' (1e5 Pa in paper mode: with pressureUnit bar
       the PbI2-He diffusivity is evaluated at p_ref = 1 bar, as in
       T-Flows);
  T    uniform, the wall table T(x) on the whole cross-section ('profile'),
       or the steady temperature of the frozen flow with the wall table as
       the wall condition ('solved', what T-Flows does; below).

U_b = Q/(pi R^2), with Q the volumetric flow of the whole pipe, given in
m3/s (flowRate; T-Flows FLOW_RATE 1.67e-6) or in mL/min (flowRateMLPerMin;
the paper's standard flow rates, used as the actual flow in paper mode).
On a wedge the inlet carries the fraction theta/(2 pi) of Q (theta from the
two wedge patches, or sectorAngle).  From the wedge patches theta is pi -
acos(n1 . n2) of the normals of the two wedge planes, exact only to their
rounding: with rescaleToFacets the inflow of the 5-degree wedge of the
tests is Q theta/(2 pi) to 2.5e-13, with sectorAngle 5 to 5e-17 (both
summed exactly from the written phi; review of M4, round 3, criterion
M4.1c).  Give sectorAngle where the inflow must be exact to round-off.
The faceted cross-section of the mesh
is smaller than the circle (by 0.5 % for the 36 sides of run/001): with
rescaleToFacets the amplitude is scaled so that the discrete inflow is
exactly Q theta/(2 pi); without it (T-Flows: u_max = 2 FLOW_RATE/(pi R^2))
the discrete inflow is the exact flux of the analytic profile through the
faceted inlet.  The log gives U_b, Q, the discrete inflow, the facet-area
ratio (inlet area over the area of the circular sector) and the continuity
check, max over the cells of |sum_f phi| relative to the inflow, and the
net flux through all patches.

The walls.  The exact flux of the axial profile through a wall face is 0
only when the face is parallel to the axis, as on a mesh extruded along it
(the gmsh pipe of run/001, blockMesh wedges).  On a wall face tilted
against the axis (a tetrahedral mesh) it is not, and the wall would be
permeable: the paired species of the solver, whose wall condition must be
zeroGradient, would be carried through it.  The log therefore gives the
flux through the 'walls' patches (sum and largest |phi| relative to the
inflow; review of M4, round 1), and the utility stops when the sum exceeds
1e-12 of the inflow, unless permeableWalls is yes (then a warning).  On the
gmsh pipe the wall faces carry +-1e-25 kg/s (8e-16 of the inflow in total):
the exact flux through faces that gmsh wrote with the last bits of their
coordinates differing along the axis.  They are kept, as the exact
integral of every face (the continuity of the cells next to the wall
depends on them at that level), and the solver books them as the
transport through WALL.

The input against the mesh (review of M4, round 2).  The utility stops
without writing when the analytic profile does not enter the domain
through 'inlet' (its flux there is not positive: an axis pointing out of
the inlet, or 'inlet' naming the outflow patch), and when the radius is
smaller than the largest distance of a wall point from the axis by more
than 1e-9 relative (the parabola would be negative near the wall, with
backflow through those faces while U is clipped to 0 there), unless
allowBackflow is yes.  Before, both wrote a reversed or partly reversed
carrier, and rescaleToFacets hid them (a negative amplitude; an inflow of
exactly Q again).

How T-Flows applies the wall temperature (origin/advanced_scalar_model,
Tests/Laminar/Scalar_Transport_Gpu/LESTO_HKS; the same in Lesto_Thesis and
in multiple_scalars_gpu_francesco:Tests/Laminar/LESTO-project): control
has HEAT_TRANSFER yes and BOUNDARY_CONDITION internal_wall, TYPE wall,
VARIABLES x u v w t q_01 q_02, FILE wall_x_profile_complete.dat, i.e. a
Dirichlet wall temperature interpolated linearly in x at every wall face
(Read_Controls_Mod/Boundary_Conditions.f90); Main_Pro_Frozen.f90 solves the
energy equation (Process % Compute_Energy) in every outer iteration of
every time step, with the frozen velocity, density 1, HEAT_CAPACITY 5193
J/(kg K) and the NIST helium conductivity of User_Mod/Types.f90 (k_he at 20
K steps, linear); the inlet is a Dirichlet 900 K (LESTO_HKS control), the
outlet an outflow (zero gradient).  T is therefore SOLVED with a wall
condition, not imposed on the cross-section.

Beyond the last row of the table T-Flows has no wall temperature: its
count line says 66 rows (x <= 0.65 m) of the 70 in the file, a wall face
whose centre lies in no interval of the table is never assigned
(Read_Controls_Mod/Boundary_Conditions.f90:538-607 of the GPU solver), so
it keeps the boundary value 0 of its allocation (Var_Mod/
Create_Variable.f90:38), which the copy of all boundary values (:670-683)
makes the wall temperature, and the energy equation treats the whole
region as a Dirichlet wall (Process_Mod/Insert_Energy_Bc.f90,
Form_Energy_Matrix.f90): in LESTO_HKS the wall beyond x = 0.65 m is a 0 K
wall.  The conductivity lookup of User_Mod/Beginning_Of_Iteration.fpp:33-36
clamps only above (dataIndex <= 50), so a cell cooled below 250 K
extrapolates k_He below the table and one below 230 K reads it out of
bounds.  Neither 'rows count' with 'outOfRange hold' (the wall beyond 0.65
m held at T(0.65 m) = 324 K) nor 'rows all' (the measured rows to 0.68 m,
run/014) reproduces that; it is not reproduced here (review of M4, round
3; plan section 29, item 5).  From its initial 500 K the
field reaches its steady state within a few thermal diffusion times R^2/
alpha = 0.1 s (alpha = k/(rho cp) with rho = 1), long before the gas
reaches the deposition zone (a few seconds), so mode 'solved' computes
that steady state once:

  div(rho cp phi_v T) - laplacian(k(T), T) = 0   (phi_v the volume flux),

Picard iterations on k(T) and on the convection scheme, fixedValue on the
wall patches (the table at the face centre) and the inlet, zeroGradient on
every other patch (constraint patches keep their type).  The written T has
calculated patches holding the solved boundary values.

system/setFrozenCarrierDict (-dict for another file)

  axis          (1 0 0);     // direction of the pipe axis
  origin        (0 0 0);     // a point on the axis
  radius        2.4e-3;      // R [m] (T-Flows R_pipe)
  flowRate      1.67e-6;     // Q [m3/s] of the whole pipe, or
  // flowRateMLPerMin 106;   // Q [mL/min]
  rescaleToFacets no;        // yes: the discrete inflow is exactly
                             // Q theta/(2 pi)
  // sectorAngle 5;          // [deg], default: from two wedge patches
                             // (to the rounding of their normals), else 360
  inlet         INLET;       // the inflow patch
  walls         (WALL);      // U = 0, and the wall table of T
  permeableWalls no;         // yes: accept a flux through the walls
                             // above round-off (a tilted wall face)
  allowBackflow no;          // yes: accept a radius below the largest
                             // distance of a wall point from the axis
  density       1;           // rho [kg/m3]
  pressure      1e5;         // p [Pa]
  temperature
  {
      mode      solved;      // uniform | profile | solved
      value     700;         // uniform [K]
      table                  // profile and solved: the wall table T(x)
      {
          file     "<constant>/wall_x_profile_complete.dat";
          format   TFlows;   // TFlows (a count line, then the rows) | columns
          rows     all;      // TFlows: all rows, or 'count' (the rows
                             // T-Flows reads; it leaves the wall beyond
                             // them at 0 K, see above)
          xColumn  0;        // 0-based columns (TFlows default 0 and 4,
          TColumn  4;        // columns default 0 and 1)
          xScale   1;        // x in the file times xScale = x [m]
          outOfRange fatal;  // fatal | hold (beyond 1e-9 of the range)
      }
      // solved:
      inletTemperature 900;  // [K]
      heatCapacity     5193; // cp [J/(kg K)]
      conductivity     ((250 0.13754) (270 0.14503) ...);
                             // (T [K] k [W/(m K)]), linear, extrapolated
                             // linearly beyond the ends
      convection       "Gauss linearUpwind grad(T)"; // default
      laplacian        "Gauss linear corrected";     // default
      solver { solver PBiCGStab; preconditioner DILU;
               tolerance 1e-14; relTol 0; }
      tolerance        1e-9; // [K], Picard: max |T - T_previous|
      maxIter          200;
  }

Usage

  setFrozenCarrier [-dict <file>]        write rho, U, phi, p and T (with
                                         -parallel: into every processor
                                         directory) into the time where the
                                         solver reads the carrier, the
                                         instance of rho that findInstance
                                         finds from the start time
                                         (writeFrozenFields no), or the
                                         start time when there is no rho
  setFrozenCarrier -check [-parallel]    only the continuity check (and the
                                         wall flux) of the phi that
                                         findInstance finds from the start
                                         time (0 for a carrier that the
                                         solver does not write), e.g. of a
                                         case decomposed by decomposePar,
                                         which writes every processor patch
                                         and, in ascii, every value at its
                                         precision (plan section 15, item
                                         7), and a patch with a single face
                                         on a rank as a 'uniform' value of
                                         10 digits (review of M4, round 1):
                                         run it after every decomposition,
                                         or write the carrier with -parallel

The fields are written with 17 significant digits for their text entries
(a 'uniform' entry, and everything with writeFormat ascii), so that a
uniform value reads back exactly.  In parallel, the two sides of a
processor face evaluate the same triangles (the reversed face starts at
the same point), each triangle in a canonical order of its vertices and
midpoints, and the triangles of a face in a canonical order, so that the
two fluxes are exactly opposite and the continuity of every cell holds to
round-off on the decomposed mesh as well.

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wedgePolyPatch.H"
#include "processorFvPatch.H"
#include "unitConversion.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "convectionScheme.H"
#include "laplacianScheme.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "SpanStream.H"

#include <algorithm>
#include <array>
#include <vector>

using namespace Foam;


/*------------------------------------------------------------------------------
Local helpers
------------------------------------------------------------------------------*/

namespace {

  /*----------------------------------------------------------------------------
  The analytic profile u(p) = 2 U_b (1 - r^2/R^2), r the distance of p from
  the axis; a quadratic polynomial of the coordinates.
  ----------------------------------------------------------------------------*/
  struct parabolicProfile {

    vector axis;     /* unit vector */
    point  origin;
    scalar radius;
    scalar Ub;

    scalar radiusOf(const point& p) const {
      const vector d = p - origin;
      return mag(d - (d & axis)*axis);
    }

    scalar operator()(const point& p) const {
      const vector d = p - origin;
      const vector perp = d - (d & axis)*axis;
      return 2*Ub*(1 - magSqr(perp)/sqr(radius));
    }
  };

  /* lexicographic order of two points */
  bool lessPoint(const point& a, const point& b) {
    return std::lexicographical_compare(a.cbegin(), a.cend(),
                                        b.cbegin(), b.cend());
  }

  /*----------------------------------------------------------------------------
  The flux of u axis through the oriented triangle (a, b, c):
  (axis . S_abc) (u(m_ab) + u(m_bc) + u(m_ca))/3, S_abc = (b - a) x (c - a)/2,
  exact for the quadratic profile.  Evaluated in a canonical order: the
  vertices sorted lexicographically (the area vector of the sorted triangle
  times the sign of the permutation) and the three midpoint values sorted,
  so that the same triangle seen with the opposite orientation (the
  neighbour's copy of a processor face) gives exactly the opposite value.
  ----------------------------------------------------------------------------*/
  scalar triangleFlux (
    const point& a,
    const point& b,
    const point& c,
    const parabolicProfile& u
  ) {
    std::array<const point*, 3> v{{&a, &b, &c}};
    int parity = 1;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2 - i; ++j) {
        if (lessPoint(*v[j + 1], *v[j])) {
          std::swap(v[j], v[j + 1]);
          parity = -parity;
        }
      }
    }
    const point& p0 = *v[0];
    const point& p1 = *v[1];
    const point& p2 = *v[2];

    const vector S = 0.5*((p1 - p0) ^ (p2 - p0));
    std::array<scalar, 3> m{{
      u(0.5*(p0 + p1)), u(0.5*(p1 + p2)), u(0.5*(p0 + p2))
    }};
    std::sort(m.begin(), m.end());

    return parity*(u.axis & S)*((m[0] + m[1] + m[2])/3.0);
  }

  /*----------------------------------------------------------------------------
  The exact flux of u axis through a face: the triangles of the fan from
  its first point (the same triangles for a face and its reverse, which
  starts at the same point), added in a canonical order (by their sorted
  vertices), so that a face and its reverse give exactly opposite values.
  ----------------------------------------------------------------------------*/
  scalar faceFlux (
    const face&             f,
    const pointField&       points,
    const parabolicProfile& u
  ) {
    const label n = f.size();
    if (n < 3) {
      return 0;
    }

    /* per triangle: its canonical key (the vertices sorted) and flux */
    struct item {
      std::array<point, 3> key;
      scalar flux;
    };
    std::vector<item> items;
    items.reserve(n - 2);

    const point& p0 = points[f[0]];
    for (label i = 1; i < n - 1; ++i) {
      const point& p1 = points[f[i]];
      const point& p2 = points[f[i + 1]];
      item t;
      t.key = {{p0, p1, p2}};
      std::sort(t.key.begin(), t.key.end(), lessPoint);
      t.flux = triangleFlux(p0, p1, p2, u);
      items.push_back(t);
    }

    std::sort(items.begin(), items.end(), [](const item& x, const item& y) {
      for (int k = 0; k < 3; ++k) {
        if (lessPoint(x.key[k], y.key[k])) {
          return true;
        }
        if (lessPoint(y.key[k], x.key[k])) {
          return false;
        }
      }
      return false;
    });

    scalar sum = 0;
    for (const item& t : items) {
      sum += t.flux;
    }
    return sum;
  }

  /*----------------------------------------------------------------------------
  A piecewise-linear table y(x), x ascending.  Within the range (and up to
  1e-9 of it beyond the ends: the round-off of mesh coordinates) it
  interpolates linearly, as T-Flows does for a profile file: for x in
  [x_m, x_m+1], w = (x_m+1 - x)/(x_m+1 - x_m), y = w y_m + (1 - w) y_m+1.
  Beyond: extrapolated with the end intervals, held, or a FatalError.
  ----------------------------------------------------------------------------*/
  class linearTable {

    std::vector<double> x_, y_;
    word   name_;
    word   outOfRange_;     /* fatal | hold | extrapolate */

  public:

    linearTable() {}

    linearTable (
      const std::vector<double>& x,
      const std::vector<double>& y,
      const word& name,
      const word& outOfRange
    )
    :
      x_(x), y_(y), name_(name), outOfRange_(outOfRange)
    {
      if (x_.size() < 2 || x_.size() != y_.size()) {
        FatalErrorInFunction << name_ << ": a table needs at least two "
          << "rows of (x y), not " << label(x_.size()) << "."
          << exit(FatalError);
      }
      for (std::size_t i = 1; i < x_.size(); ++i) {
        if (!(x_[i] > x_[i - 1])) {
          FatalErrorInFunction << name_ << ": the abscissae must increase "
            << "strictly; row " << label(i + 1) << " has " << x_[i]
            << " after " << x_[i - 1] << "." << exit(FatalError);
        }
      }
    }

    double xMin() const { return x_.front(); }
    double xMax() const { return x_.back(); }
    label  size() const { return label(x_.size()); }

    double operator()(double x) const {
      const double tol = 1e-9*(x_.back() - x_.front());
      if (x < x_.front() - tol || x > x_.back() + tol) {
        if (outOfRange_ == "fatal") {
          FatalErrorInFunction << name_ << ": " << x << " lies outside the "
            << "table [" << x_.front() << ", " << x_.back() << "] "
            << "(outOfRange fatal; hold keeps the end values)."
            << exit(FatalError);
        }
        if (outOfRange_ == "hold") {
          x = min(max(x, x_.front()), x_.back());
        }
      } else {
        x = min(max(x, x_.front()), x_.back());
      }

      /* the interval [x_m, x_m+1] that holds x (the end intervals beyond) */
      std::size_t m = std::upper_bound(x_.begin(), x_.end(), x) - x_.begin();
      m = m == 0 ? 0 : m - 1;
      m = std::min(m, x_.size() - 2);
      const double w = (x_[m + 1] - x)/(x_[m + 1] - x_[m]);
      return w*y_[m] + (1 - w)*y_[m + 1];
    }
  };

  /*----------------------------------------------------------------------------
  The wall table of T(x) from a file.  Lines whose first character (after
  blanks) is '#', '!' or '%' are comments (T-Flows File_Mod Read_Line).
  format TFlows: the first remaining line is the number of rows n (T-Flows
  reads exactly n rows after it; rows count does the same, rows all reads
  every remaining row and reports a mismatch); format columns: every
  remaining line is a row.  x and T are the columns xColumn and TColumn
  (0-based), x times xScale.
  ----------------------------------------------------------------------------*/
  linearTable readWallTable(const dictionary& d) {

    fileName file(d.get<fileName>("file"));
    file.expand();
    const word format(d.getOrDefault<word>("format", "TFlows"));
    const bool tflows = format == "TFlows";
    if (!tflows && format != "columns") {
      FatalIOErrorInFunction(d) << "Unknown format " << format
        << ".  Choose TFlows or columns." << exit(FatalIOError);
    }
    const word rows(d.getOrDefault<word>("rows", "all"));
    if (rows != "all" && rows != "count") {
      FatalIOErrorInFunction(d) << "Unknown rows " << rows
        << ".  Choose all or count." << exit(FatalIOError);
    }
    const label xColumn = d.getOrDefault<label>("xColumn", 0);
    const label TColumn = d.getOrDefault<label>("TColumn", tflows ? 4 : 1);
    const scalar xScale = d.getOrDefault<scalar>("xScale", 1);
    const word outOfRange(d.getOrDefault<word>("outOfRange", "fatal"));
    if (outOfRange != "fatal" && outOfRange != "hold") {
      FatalIOErrorInFunction(d) << "Unknown outOfRange " << outOfRange
        << ".  Choose fatal or hold." << exit(FatalIOError);
    }

    IFstream is(file);
    if (!is.good()) {
      FatalIOErrorInFunction(d) << "Cannot read the wall table " << file
        << exit(FatalIOError);
    }

    std::vector<double> x, T;
    label count = -1, lineNo = 0, rowsInFile = 0;
    string line;
    while (is.getLine(line)) {
      ++lineNo;
      const std::size_t first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line[first] == '#'
        || line[first] == '!' || line[first] == '%') {
        continue;
      }
      ISpanStream ls(line);
      std::vector<double> values;
      while (true) {
        token t(ls);
        if (!t.good() || !t.isNumber()) {
          break;
        }
        values.push_back(t.number());
      }
      if (tflows && count < 0) {
        if (values.size() != 1) {
          FatalIOErrorInFunction(d) << file << ", line " << lineNo << ": "
            << "format TFlows expects the number of rows here." << nl
            << exit(FatalIOError);
        }
        count = label(values[0]);
        continue;
      }
      ++rowsInFile;
      if (rows == "count" && label(x.size()) == count) {
        continue;
      }
      if (label(values.size()) <= max(xColumn, TColumn)) {
        FatalIOErrorInFunction(d) << file << ", line " << lineNo << ": "
          << "fewer than " << max(xColumn, TColumn) + 1 << " numbers."
          << exit(FatalIOError);
      }
      x.push_back(xScale*values[xColumn]);
      T.push_back(values[TColumn]);
    }

    Info<< "Wall table " << d.get<fileName>("file") << " (format " << format;
    if (tflows) {
      Info<< ", count line " << count << ", " << rowsInFile << " rows in the "
        << "file, rows " << rows << ": " << label(x.size()) << " read";
      if (count != rowsInFile) {
        Info<< "; the count differs from the rows in the file: T-Flows reads "
          << count << " rows";
      }
    }
    Info<< "): " << label(x.size()) << " rows, x "
      << (x.empty() ? 0 : x.front()) << " to " << (x.empty() ? 0 : x.back())
      << " m, outOfRange " << outOfRange << nl;

    return linearTable(x, T, file, outOfRange);
  }

  /*----------------------------------------------------------------------------
  The flux through the wall patches: the sum over their faces of |phi| and
  the largest |phi| [kg/s].  A wall face parallel to the axis carries none
  (header: 'The walls').  Collective.
  ----------------------------------------------------------------------------*/
  FixedList<scalar, 2> wallFlux (
    const surfaceScalarField& phi,
    const labelList&          walls
  ) {
    FixedList<scalar, 2> flux({0, 0});
    for (const label patchi : walls) {
      for (const scalar value : phi.boundaryField()[patchi]) {
        flux[0] += mag(value);
        flux[1] = max(flux[1], mag(value));
      }
    }
    reduce(flux[0], sumOp<scalar>());
    reduce(flux[1], maxOp<scalar>());
    return flux;
  }

  /*----------------------------------------------------------------------------
  The continuity of a flux field: max over the cells of |sum_f phi| (all
  faces of the cell, processor and wedge faces included; empty patches
  carry no values), the inflow through the inlet patch and the net flux
  through all non-processor patches; and the flux through the wall patches
  (wallFlux()), relative to the inflow, which it returns.  Collective.
  ----------------------------------------------------------------------------*/
  scalar checkContinuity (
    const fvMesh&              mesh,
    const surfaceScalarField&  phi,
    const label                inletPatch,
    const labelList&           walls,
    const string&              what
  ) {
    scalarField sum(mesh.nCells(), Zero);
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(neighbour, facei) {
      sum[owner[facei]] += phi[facei];
      sum[neighbour[facei]] -= phi[facei];
    }

    scalar inflow = 0, net = 0;
    forAll(phi.boundaryField(), patchi) {
      const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
      const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
      forAll(pphi, facei) {
        sum[faceCells[facei]] += pphi[facei];
      }
      if (!isA<processorFvPatch>(mesh.boundary()[patchi])) {
        net += Foam::sum(pphi);
      }
      if (patchi == inletPatch) {
        inflow -= Foam::sum(pphi);
      }
    }
    reduce(inflow, sumOp<scalar>());
    reduce(net, sumOp<scalar>());

    const scalar maxError = gMax(mag(sum));
    const scalar reference = max(mag(inflow), VSMALL);

    const FixedList<scalar, 2> wall(wallFlux(phi, walls));

    Info.stream().precision(6);
    Info<< "Continuity (" << what.c_str() << "): max over the cells of "
      << "|sum_f phi| = " << maxError << " kg/s = " << maxError/reference
      << " of the inflow " << inflow << " kg/s; net flux through the "
      << "patches " << net << " kg/s (" << net/reference << " of the inflow)"
      << nl
      << "Wall flux (" << what.c_str() << "): sum over the faces of the "
      << "walls of |phi| = " << wall[0] << " kg/s = " << wall[0]/reference
      << " of the inflow; largest |phi| " << wall[1] << " kg/s" << nl;

    return wall[0]/reference;
  }
}


/*------------------------------------------------------------------------------
Main
------------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {

  argList::addNote (
    "Frozen paper-mode carrier of rhoFixedFlowFoam: constant rho, the exact\n"
    "flux of a parabolic pipe profile, uniform p, and T uniform, from a\n"
    "wall table, or solved with the wall table (system/setFrozenCarrierDict)"
  );
  argList::addBoolOption (
    "check",
    "only check the continuity of the phi of the start time"
  );
  argList::addOption (
    "dict",
    "file",
    "the dictionary (default system/setFrozenCarrierDict)"
  );

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  const fileName dictPath (
    args.getOrDefault<fileName>("dict", runTime.system()/"setFrozenCarrierDict")
  );
  IOdictionary dict (
    IOobject(dictPath.name(),
             dictPath.path(),
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE)
  );

  const word inletName(dict.get<word>("inlet"));
  const label inletPatch = mesh.boundaryMesh().findPatchID(inletName);
  if (inletPatch < 0) {
    FatalIOErrorInFunction(dict) << "No patch " << inletName << "."
      << exit(FatalIOError);
  }

  labelHashSet walls;
  for (const word& name : dict.get<wordList>("walls")) {
    const label patchi = mesh.boundaryMesh().findPatchID(name);
    if (patchi < 0) {
      FatalIOErrorInFunction(dict) << "No wall patch " << name << "."
        << exit(FatalIOError);
    }
    walls.insert(patchi);
  }
  const labelList wallList(walls.sortedToc());
  const bool permeableWalls = dict.getOrDefault<bool>("permeableWalls", false);

  /* the flux through the walls, relative to the inflow, that the utility
     accepts as round-off (header: 'The walls') */
  const scalar wallFluxTolerance = 1e-12;

  /*--------------------------------------------------------------------------
  -check: the continuity of the phi that the solver reads, the instance of
  phi that findInstance finds from the start time (not necessarily the
  start time itself: a case that has run with writeFrozenFields no holds
  the carrier only in 0; review of M4, round 1), e.g. decomposed.
  --------------------------------------------------------------------------*/
  if (args.found("check")) {
    const word instance(runTime.findInstance(mesh.dbDir(), "phi"));
    surfaceScalarField phi (
      IOobject("phi", instance, mesh, IOobject::MUST_READ,
               IOobject::NO_WRITE),
      mesh
    );
    const scalar wall = checkContinuity (
      mesh, phi, inletPatch, wallList, "phi read at time " + instance
      + (instance == runTime.timeName() ? string("")
         : " (findInstance from " + runTime.timeName() + ")")
      + (Pstream::parRun()
         ? " on " + Foam::name(Pstream::nProcs()) + " ranks" : " in serial")
    );
    if (wall > wallFluxTolerance) {
      WarningInFunction
        << "The walls carry " << wall << " of the inflow (more than "
        << wallFluxTolerance << "): the wall faces are not parallel to the "
        << "axis, and the wall is permeable." << nl << endl;
    }
    Info<< "End" << endl;
    return 0;
  }

  /*--------------------------------------------------------------------------
  The pipe and the flow rate
  --------------------------------------------------------------------------*/
  parabolicProfile u;
  u.axis = dict.getOrDefault<vector>("axis", vector(1, 0, 0));
  u.origin = dict.getOrDefault<point>("origin", Zero);
  u.radius = dict.get<scalar>("radius");
  if (!(mag(u.axis) > 0) || !(u.radius > 0)) {
    FatalIOErrorInFunction(dict) << "Require an axis of non-zero length and "
      << "radius > 0." << exit(FatalIOError);
  }
  u.axis /= mag(u.axis);

  const bool haveQ = dict.found("flowRate");
  const bool haveQml = dict.found("flowRateMLPerMin");
  if (haveQ == haveQml) {
    FatalIOErrorInFunction(dict) << "Give exactly one of flowRate [m3/s] or "
      << "flowRateMLPerMin [mL/min]." << exit(FatalIOError);
  }
  const scalar Q = haveQ
    ? dict.get<scalar>("flowRate")
    : dict.get<scalar>("flowRateMLPerMin")*1e-6/60.0;
  if (!(Q > 0)) {
    FatalIOErrorInFunction(dict) << "The flow rate must be positive."
      << exit(FatalIOError);
  }

  const scalar density = dict.getOrDefault<scalar>("density", 1);
  const scalar pressure = dict.getOrDefault<scalar>("pressure", 1e5);
  const bool rescale = dict.getOrDefault<bool>("rescaleToFacets", false);
  if (!(density > 0) || !(pressure > 0)) {
    FatalIOErrorInFunction(dict) << "Require density > 0 and pressure > 0."
      << exit(FatalIOError);
  }

  /*--------------------------------------------------------------------------
  Sector of the mesh: sectorAngle [deg], else from two wedge patches (the
  angle between the wedge planes, pi - the angle between their normals),
  else the whole pipe.
  --------------------------------------------------------------------------*/
  scalar sectorAngle = constant::mathematical::twoPi;
  word sectorFrom("the whole pipe");
  if (dict.found("sectorAngle")) {
    sectorAngle = degToRad(dict.get<scalar>("sectorAngle"));
    sectorFrom = "sectorAngle";
  } else {
    DynamicList<vector> normals;
    for (const polyPatch& pp : mesh.boundaryMesh()) {
      if (isA<wedgePolyPatch>(pp)) {
        normals.append(refCast<const wedgePolyPatch>(pp).n());
      }
    }
    if (normals.size() == 2) {
      const scalar c = min(max(normals[0] & normals[1], scalar(-1)), scalar(1));
      sectorAngle = constant::mathematical::pi - std::acos(c);
      sectorFrom = "the two wedge patches";
    } else if (normals.size() != 0) {
      FatalIOErrorInFunction(dict) << normals.size() << " wedge patches: "
        << "give sectorAngle [deg]." << exit(FatalIOError);
    }
  }
  const scalar sector = sectorAngle/constant::mathematical::twoPi;
  const scalar Qsector = Q*sector;

  /*--------------------------------------------------------------------------
  The time the carrier is written into: where the solver reads it with
  writeFrozenFields no, the instance of rho that findInstance finds from
  the start time (e.g. 0 of a case that has run), or the start time when
  there is no rho yet (review of M4, round 1: the carrier went into the
  latest time of a case with startFrom latestTime).
  --------------------------------------------------------------------------*/
  word instance (
    runTime.findInstance(mesh.dbDir(), "rho", IOobject::READ_IF_PRESENT,
                         word::null, false)
  );
  const bool newCarrier = instance.empty();
  const word startName(runTime.timeName());
  if (newCarrier) {
    instance = startName;
  }

  /*--------------------------------------------------------------------------
  phi: the exact flux of the profile through every face (above), for the
  analytic U_b = Q/(pi R^2), then (rescaleToFacets) again with U_b scaled
  so that the discrete inflow is Q sector.
  --------------------------------------------------------------------------*/
  const scalar circleArea = constant::mathematical::pi*sqr(u.radius);
  u.Ub = Q/circleArea;

  const pointField& points = mesh.points();
  const faceList& faces = mesh.faces();

  surfaceScalarField phi (
    IOobject("phi", instance, mesh, IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimMass/dimTime, Zero)
  );
  phi.setOriented();

  auto setFluxes = [&]() {
    forAll(phi, facei) {
      phi[facei] = density*faceFlux(faces[facei], points, u);
    }
    forAll(phi.boundaryField(), patchi) {
      fvsPatchScalarField& pphi = phi.boundaryFieldRef()[patchi];
      const label start = mesh.boundaryMesh()[patchi].start();
      forAll(pphi, i) {
        pphi[i] = density*faceFlux(faces[start + i], points, u);
      }
    }
  };

  /* the discrete inflow [m3/s] and the inlet area */
  auto inflowOf = [&]() {
    return -gSum(phi.boundaryField()[inletPatch])/density;
  };

  /*--------------------------------------------------------------------------
  The geometry of the input against the mesh (review of M4, round 2): the
  radius must reach the wall, and the profile must enter through 'inlet'.
  Either mistake gave a reversed carrier without a message (an axis
  pointing out of the inlet, or 'inlet' naming the outflow patch: a
  negative inflow, which rescaleToFacets turned into a negative amplitude;
  a radius below the wall: the parabola negative near the wall, backflow
  in phi while U is clipped to 0 there), and the rescale hid both.
  --------------------------------------------------------------------------*/

  /* the largest distance of a point of the wall patches from the axis */
  scalar wallRadius = 0;
  for (const label patchi : walls) {
    for (const label pointi : mesh.boundaryMesh()[patchi].meshPoints()) {
      wallRadius = max(wallRadius, u.radiusOf(points[pointi]));
    }
  }
  reduce(wallRadius, maxOp<scalar>());

  /* beyond round-off of the point coordinates: a radius below the wall */
  const scalar radiusTolerance = 1e-9;
  const bool allowBackflow = dict.getOrDefault<bool>("allowBackflow", false);
  if (wallRadius > u.radius*(1 + radiusTolerance) && !allowBackflow) {
    FatalIOErrorInFunction(dict)
      << "The radius " << u.radius << " m is smaller than the largest "
      << "distance of a point of the wall patches from the axis, "
      << wallRadius << " m (relative " << wallRadius/u.radius - 1
      << ", more than " << radiusTolerance << "): the parabola "
      << "2 U_b (1 - r^2/R^2) would be negative near the wall, with backflow "
      << "through those faces." << nl << "Give the radius of the mesh, or "
      << "set allowBackflow yes to accept it." << nl << exit(FatalIOError);
  }

  setFluxes();
  const scalar inflowAnalytic = inflowOf();
  if (!(inflowAnalytic > 0)) {
    FatalIOErrorInFunction(dict)
      << "The flux of the profile through the patch " << inletName << " is "
      << inflowAnalytic << " m3/s, not an inflow: the axis " << u.axis
      << " must point downstream, into the domain through " << inletName
      << ", and 'inlet' must name the inflow patch (not the outlet)." << nl
      << exit(FatalIOError);
  }
  scalar scale = 1;
  if (rescale) {
    scale = Qsector/inflowAnalytic;
    u.Ub *= scale;
    setFluxes();
  }
  const scalar inflow = inflowOf();
  const scalar inletArea =
    gSum(mag(u.axis & mesh.Sf().boundaryField()[inletPatch]));

  Info.stream().precision(10);
  Info<< nl << "Paper-mode carrier (setFrozenCarrier): R = " << u.radius
    << " m (largest distance of a wall point from the axis " << wallRadius
    << " m), axis " << u.axis << " through " << u.origin << nl
    << "  Q = " << Q << " m3/s ("
    << Q*60e6 << " mL/min) of the whole pipe; sector " << radToDeg(sectorAngle)
    << " deg (" << sectorFrom.c_str() << "): Q sector = " << Qsector
    << " m3/s" << nl
    << "  U_b = Q/(pi R^2) = " << Q/circleArea << " m/s";
  if (rescale) {
    Info<< ", rescaled to the facets by " << scale << ": U_b = " << u.Ub
      << " m/s";
  } else {
    Info<< " (not rescaled to the facets)";
  }
  Info<< ", centreline velocity 2 U_b = " << 2*u.Ub << " m/s" << nl
    << "  discrete inflow through " << inletName << " = " << inflow
    << " m3/s = " << inflow/Qsector << " of Q sector (analytic profile: "
    << inflowAnalytic/Qsector << "); inlet area " << inletArea << " m2, "
    << "facet-area ratio (inlet area/(pi R^2 sector)) "
    << inletArea/(circleArea*sector) << "; discrete bulk velocity "
    << inflow/inletArea << " m/s" << nl
    << "  rho = " << density << " kg/m3, p = " << pressure << " Pa" << nl;

  /*--------------------------------------------------------------------------
  The walls must carry no flux beyond round-off (header: 'The walls'), or
  the carrier is not written (unless permeableWalls).
  --------------------------------------------------------------------------*/
  {
    const FixedList<scalar, 2> wall(wallFlux(phi, wallList));
    const scalar relative = wall[0]/max(mag(inflow*density), VSMALL);
    if (relative > wallFluxTolerance) {
      if (!permeableWalls) {
        FatalIOErrorInFunction(dict)
          << "The exact flux of the axial profile through the walls is "
          << wall[0] << " kg/s in total (" << relative << " of the inflow, "
          << "more than " << wallFluxTolerance << "; largest face "
          << wall[1] << " kg/s): the wall faces are not parallel to the "
          << "axis (e.g. a tetrahedral mesh), and the wall would be "
          << "permeable.  The solver assumes a zero flux through the wall "
          << "of a paired species." << nl << "Use a mesh extruded along the "
          << "axis, or set permeableWalls yes to accept it." << nl
          << exit(FatalIOError);
      }
      WarningInFunction
        << "The walls carry " << relative << " of the inflow (largest face "
        << wall[1] << " kg/s): the wall faces are not parallel to the axis, "
        << "and the wall is permeable (permeableWalls yes)." << nl << endl;
    }
  }

  /*--------------------------------------------------------------------------
  rho, p, U
  --------------------------------------------------------------------------*/
  volScalarField rho (
    IOobject("rho", instance, mesh, IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimDensity, density),
    calculatedFvPatchScalarField::typeName
  );
  volScalarField p (
    IOobject("p", instance, mesh, IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimPressure, pressure),
    calculatedFvPatchScalarField::typeName
  );
  volVectorField U (
    IOobject("U", instance, mesh, IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero),
    calculatedFvPatchVectorField::typeName
  );

  const vectorField& C = mesh.C().primitiveField();
  forAll(U, celli) {
    U[celli] = max(u(C[celli]), scalar(0))*u.axis;
  }
  forAll(U.boundaryField(), patchi) {
    fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];
    const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
    forAll(pU, i) {
      pU[i] = walls.found(patchi)
        ? vector(Zero)
        : vector(max(u(Cf[i]), scalar(0))*u.axis);
    }
  }

  /*--------------------------------------------------------------------------
  T: uniform, the wall table on the whole cross-section (profile), or the
  steady energy equation with the wall table (solved; header)
  --------------------------------------------------------------------------*/
  const dictionary& Tdict = dict.subDict("temperature");
  const word mode(Tdict.get<word>("mode"));
  if (mode != "uniform" && mode != "profile" && mode != "solved") {
    FatalIOErrorInFunction(Tdict) << "Unknown temperature mode " << mode
      << ".  Choose uniform, profile or solved." << exit(FatalIOError);
  }

  volScalarField T (
    IOobject("T", instance, mesh, IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimTemperature,
                      mode == "uniform" ? Tdict.get<scalar>("value") : 300),
    calculatedFvPatchScalarField::typeName
  );

  linearTable wallT;
  if (mode != "uniform") {
    wallT = readWallTable(Tdict.subDict("table"));
  }

  /* the axial coordinate of a point */
  auto axial = [&u](const point& x) { return (x - u.origin) & u.axis; };

  if (mode == "profile") {
    forAll(T, celli) {
      T[celli] = wallT(axial(C[celli]));
    }
    forAll(T.boundaryField(), patchi) {
      fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
      const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
      forAll(pT, i) {
        pT[i] = wallT(axial(Cf[i]));
      }
    }
  }

  if (mode == "solved") {

    const scalar Tin = Tdict.get<scalar>("inletTemperature");
    const scalar cp = Tdict.get<scalar>("heatCapacity");
    const List<Tuple2<scalar, scalar>> kTable (
      Tdict.get<List<Tuple2<scalar, scalar>>>("conductivity")
    );
    std::vector<double> kT, kk;
    for (const Tuple2<scalar, scalar>& row : kTable) {
      kT.push_back(row.first());
      kk.push_back(row.second());
    }
    const linearTable conductivity(kT, kk, "conductivity", "extrapolate");
    const scalar tolerance = Tdict.getOrDefault<scalar>("tolerance", 1e-9);
    const label maxIter = Tdict.getOrDefault<label>("maxIter", 200);
    dictionary schemes (
      IStringStream (
        "convection "
        + Tdict.getOrDefault<string>("convection",
                                     "Gauss linearUpwind grad(T)")
        + "; laplacian " + Tdict.getOrDefault<string>("laplacian",
                                                      "Gauss linear corrected")
        + ";"
      )()
    );
    const dictionary solverDict (
      Tdict.found("solver")
      ? Tdict.subDict("solver")
      : dictionary(IStringStream("solver PBiCGStab; preconditioner DILU; "
                                 "tolerance 1e-14; relTol 0;")())
    );

    /* the scheme entries, read from the start at every iteration */
    auto schemeStream = [&schemes](const word& name) -> ITstream& {
      ITstream& is = schemes.lookup(name);
      is.rewind();
      return is;
    };

    /* fixedValue on the walls and the inlet, zeroGradient elsewhere;
       constraint patches keep their own types */
    wordList types(mesh.boundary().size());
    forAll(types, patchi) {
      const polyPatch& pp = mesh.boundaryMesh()[patchi];
      types[patchi] = polyPatch::constraintType(pp.type()) ? pp.type()
        : (walls.found(patchi) || patchi == inletPatch)
        ? fixedValueFvPatchScalarField::typeName
        : zeroGradientFvPatchScalarField::typeName;
    }

    volScalarField Ts (
      IOobject("Tsolved", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar(dimTemperature, Tin),
      types
    );
    forAll(Ts.boundaryField(), patchi) {
      fvPatchScalarField& pT = Ts.boundaryFieldRef()[patchi];
      const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
      if (walls.found(patchi)) {
        forAll(pT, i) {
          pT[i] = wallT(axial(Cf[i]));
        }
      } else if (patchi == inletPatch) {
        pT = Tin;
      }
    }
    /* the initial guess: the wall table along the axis */
    forAll(Ts, celli) {
      Ts[celli] = wallT(axial(C[celli]));
    }
    Ts.correctBoundaryConditions();

    volScalarField kappa (
      IOobject("kappa", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
      calculatedFvPatchScalarField::typeName
    );
    const surfaceScalarField phiCp(phi*dimensionedScalar(dimEnergy/dimMass
                                   /dimTemperature, cp));

    label iter = 0;
    scalar change = GREAT;
    for (; iter < maxIter && change > tolerance; ++iter) {

      forAll(kappa, celli) {
        kappa[celli] = conductivity(Ts[celli]);
      }
      forAll(kappa.boundaryField(), patchi) {
        fvPatchScalarField& pk = kappa.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pT = Ts.boundaryField()[patchi];
        forAll(pk, i) {
          pk[i] = conductivity(pT[i]);
        }
      }
      kappa.boundaryFieldRef().evaluateCoupled<processorFvPatch>();

      const scalarField previous(Ts.primitiveField());

      tmp<fv::convectionScheme<scalar>> convection (
        fv::convectionScheme<scalar>::New(mesh, phiCp,
                                          schemeStream("convection"))
      );
      tmp<fv::laplacianScheme<scalar, scalar>> laplacian (
        fv::laplacianScheme<scalar, scalar>::New(mesh,
                                                 schemeStream("laplacian"))
      );
      fvScalarMatrix TEqn (
        convection.ref().fvmDiv(phiCp, Ts)
      - laplacian.ref().fvmLaplacian(kappa, Ts)
      );
      TEqn.solve(solverDict);

      change = gMax(mag(Ts.primitiveField() - previous));
    }

    Info<< "Temperature solved (steady, rho cp = " << density*cp << " J/(m3 K),"
      << " k(T) from " << label(kTable.size()) << " rows, inlet "
      << Tin << " K, the wall table on " << walls.size() << " patch(es)): "
      << iter << " Picard iterations, last change " << change << " K"
      << (change > tolerance ? " (NOT converged)" : "") << nl;

    if (change > tolerance) {
      FatalErrorInFunction << "The temperature did not converge in "
        << maxIter << " iterations (change " << change << " K > tolerance "
        << tolerance << " K)." << exit(FatalError);
    }

    T.primitiveFieldRef() = Ts.primitiveField();
    forAll(T.boundaryField(), patchi) {
      T.boundaryFieldRef()[patchi] = Ts.boundaryField()[patchi];
    }
  }

  /*--------------------------------------------------------------------------
  The wall temperature written against the table, and the range of T
  --------------------------------------------------------------------------*/
  if (mode != "uniform") {
    scalar maxDeviation = 0;
    for (const label patchi : walls) {
      const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
      forAll(Cf, i) {
        maxDeviation = max(maxDeviation,
          mag(T.boundaryField()[patchi][i] - wallT(axial(Cf[i]))));
      }
    }
    reduce(maxDeviation, maxOp<scalar>());
    Info<< "Temperature " << mode << ": T " << gMin(T.primitiveField())
      << " to " << gMax(T.primitiveField()) << " K; the wall faces hold the "
      << "table at their centres (max deviation " << maxDeviation << " K)"
      << nl;
  } else {
    Info<< "Temperature uniform " << Tdict.get<scalar>("value") << " K" << nl;
  }

  /*--------------------------------------------------------------------------
  Write with 17 significant digits for every text entry, so that a uniform
  value (and an ascii file) reads back exactly
  --------------------------------------------------------------------------*/
  IOstream::defaultPrecision(max(17u, IOstream::defaultPrecision()));

  /* regIOobject::write() moves an object whose instance is an older time
     to the current time: the Time is set to the instance first */
  if (instance != runTime.timeName()) {
    runTime.setTime(instant(instance), runTime.timeIndex());
  }

  rho.write();
  p.write();
  U.write();
  T.write();
  phi.write();

  checkContinuity(mesh, phi, inletPatch, wallList, "phi as written"
    + string(Pstream::parRun()
             ? " on " + Foam::name(Pstream::nProcs()) + " ranks" : ""));

  const string reason (
    newCarrier ? string("the start time: no rho yet")
    : "the instance of rho found from the start time " + startName
  );
  Info<< "Wrote rho, p, U, T and phi into " << instance << " ("
    << reason.c_str() << ")" << nl << "End" << endl;

  return 0;
}
