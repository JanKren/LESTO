# rhoFixedFlowFoam

## Purpose

`rhoFixedFlowFoam` advances gaseous species and accumulates solid species on an
already converged `rhoSimpleFoam` carrier flow.

The carrier fields `U`, `p`, `e/T`, `rho` and `phi` are read from the saved
steady solution and remain frozen throughout the transient species calculation.
No momentum, pressure, continuity or energy equation is solved, and there is no
SIMPLE/PIMPLE loop.

Gaseous and solid species are defined at run time in
`constant/speciesTransportProperties`.  The solver therefore does not need to
be rebuilt when another passive gaseous species is added.  A one-species case
is simply a species list containing one entry.

Gas species can be PAIRED with a condensate (a solid species) through
`constant/thermochemistryProperties`; the pair then exchanges mass with the
wall and with sample cells, with an exactly closed mass balance (see Phase
change below).  Without that file, or with `model mock;`, the solver runs
the coupled mock of the main branch (case 011, below) and behaves exactly
like the pristine build of its tag 0.1, to the last bit of its logs and
fields.

## Transported variables and scalar equation

For every gaseous species `i`, the transported field is

    Y_<speciesName>

and represents a dimensionless mass fraction:

    Y_i = kg species i / kg mixture.

Each species is advanced in physical time with

    ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = S_i

where

    rhoD_i = rho * D_i.

`D_i` is the molecular diffusivity in m2/s.  The coefficient entering the
finite-volume diffusion term is therefore `rhoD_i`, with units kg/(m s).

Solid species have no convection or diffusion and satisfy `ddt(rho,Y_i) = S_i`.
In case 011 a mock precipitation source transfers PbI2 from gas to solid in
cells adjacent to `WALL`.  It depends on temperature and local `Y_PbI2_g`;
the solid source is adjusted after the gas solve to match the realized sink.
After each time step, the solver reports the inventory, boundary fluxes,
sources and balance error for each configured species and their total.

This is `model mock`, the default (case `run/011-coupled-species-sources`;
`evaluateThermochemistry.C`, `solveGasSpecies.H`, `solveSolidSpecies.H`,
`collectSpeciesBalance.H`, `checkSpeciesConservation.H`).  It also writes
the source fields `source_<species>` [kg/(m3 s)] at every write time.  With
a phase-change model (`temperature` or `HKS`, see Phase change below) the
gas/condensate pairs exchange mass through the interface instead:
`evaluateThermochemistry()` is not called, every source field stays zero and
is not written, a paired gas is solved by `solvePairedGasSpecies.H`, a
condensate is a deposit per interface element, and the mass ledger of every
pair replaces the per-species report.  A gas species that is not paired,
e.g. a passive tracer, is still solved by `solveGasSpecies.H`; its source is
then zero, so its solution is the same to the last bit as without the
source term.  That was shown once against the build of the branch before
the rebase (plan section 31); tests M0.1c, M0.1d and M2.9 check the printed
lines and the equivalence with tag 0.1.

## Frozen carrier fields

`U`, `p`, `e/T`, `rho` and `phi` are read from the converged carrier-flow
solution.  The saved, pressure-corrected mass flux `phi` is required directly;
it is not reconstructed from `U`, because doing so could change the converged
carrier solution.

The solver also evaluates `mu` and `kappa` once from the frozen thermodynamic
state.  They are not updated because `thermo.correct()` is never called.

### Writing the frozen carrier: `writeFrozenFields`

The switch `writeFrozenFields` in the `speciesTransport` dictionary of
`system/controlDict` (see below) selects where the frozen carrier is read
from and whether it is written again.  All fields that never change during
the run share this policy: `U`, `p`, `T`, `e`, `rho`, `phi`, `mu`, `kappa`,
`D_<species>` and `rhoD_<species>`.

- `writeFrozenFields yes` (default) keeps the original behaviour.  The
  carrier is read from the start time and all these fields are written at
  every write time.  For case 010 a write directory then holds 18 files,
  about 39.6 MB (15 files and 36 MB before the source fields of the coupled
  mock of tag 0.1, which `model mock` writes).
- `writeFrozenFields no` reads the carrier from the latest time directory,
  searching back from the start time, that holds `rho`
  (`runTime.findInstance(mesh.dbDir(), "rho")`), normally `0`.  None of these
  fields is ever written.  A write directory then holds only `Y_*`, `c_*`,
  in `model mock` the source fields `source_*`, and `uniform/`: 11.1 MB for
  case 010 (7.5 MB for the four binary `Y_*`/`c_*` scalar fields on 228,960
  cells, 3.7 MB for the source fields; `source_tracer`, zero everywhere, is
  written as a uniform field).  The log states where the carrier came from:

      Frozen carrier read from 0 and not written (writeFrozenFields no)

  `fluidThermo` would read `T` and `p` from the start time.  With `no` the
  solver therefore reads them itself from the frozen instance, `NO_WRITE`,
  and registers them before `fluidThermo` is constructed; `fluidThermo`
  adopts registered fields instead of reading its own copies.
- The switch changes only the I/O: the `Y_*`, `c_*` and `source_*` fields
  written with `no` are byte-identical to those written with `yes` (case
  010).

**Restart.**  With `writeFrozenFields no` a restart needs nothing but the
`Y_*` fields in the restart directory, e.g. `startFrom latestTime` with

    0/        carrier (U, p, T, rho, phi) and the initial Y_*
    0.024/    Y_PbI2_g  Y_PbI2_s  Y_tracer

The carrier is found in `0` again.

- **Parallel.**  The carrier instance is looked up in every
  `processor*/` directory, so each must hold `0/` as well as the restart
  time.  `decomposePar -latestTime` alone produces processor directories
  without `0/`, and the run then stops with `Cannot find file "rho" ...`.
  Decompose both times instead (`decomposePar -time 0,<t>`), or keep the
  existing `processor*/0`.
- **Bit-for-bit.**  The restarted fields reproduce the continuous run bit
  for bit at the final time, in serial and in parallel (`tests/Alltest
  M1`), when the `Y_*` fields are written losslessly: `writeFormat binary`
  (as in case 010, the tested case) or ascii with `writePrecision` of at
  least 17.  Ascii at the default precision (6) rounds the restart state,
  so the restarted run is slightly different (plan section 8, risk 4; booked
  as the `restart` ledger term from M2 on).
- **`uniform/time`.**  Without `uniform/` in the restart directory the time
  index counts from the restart; the solver does not use it.  The time
  value is always taken from the directory name (`0.036` read as the
  nearest double), while a continuous run accumulates `value + deltaT`, so
  the `value` written in `uniform/time` can differ from the continuous
  run's in the last bits: at 0.048, after a restart from 0.036, it is
  `0.048000000000000001`, and in the continuous run
  `0.0480000000000000149`.  The results do not depend on it (M1.3d).

With `writeFrozenFields yes` a restart directory must hold the carrier as
well, as before; every write directory does.  The carrier is then re-read
from the restart directory, where the solver wrote it: in ascii rounded to
`writePrecision`, and a field or patch whose values are all equal (a
uniform `rho`, a uniform inlet value) as `uniform <value>` in text at
`writePrecision` even in a binary file.  Such a restart does not continue
the run bit for bit; with a phase-change model the rounding of `rho` is
booked as the `restart` ledger term (see "State, output and restart").

**Recommendation.**  Keep the default in existing cases (their output stays
unchanged); use `writeFrozenFields no` in new cases.

### Paper-mode carrier: `utilities/setFrozenCarrier` (M4)

The code-to-code cases of Lobresco et al. (2026) use the carrier of the
T-Flows reference code, not a `rhoSimpleFoam` solution: a fully developed,
non-expanding laminar pipe flow of density 1 ("paper mode", plan section
11).  `setFrozenCarrier` writes it where the solver reads it with
`writeFrozenFields no` (the instance of `rho` that findInstance finds from
the start time, e.g. 0 of a case that has run; the start time when there
is no `rho` yet):

- `rho` = `density` (1 kg/m3, T-Flows `Flow % density = 1`), `p` =
  `pressure` (1e5 Pa: with `pressureUnit bar` the PbI2-He diffusivity is
  evaluated at T-Flows' p_ref = 1 bar), both uniform;
- `U = 2 U_b (1 - r^2/R^2)` along `axis`, r the distance from the axis,
  `U_b = Q/(pi R^2)` with Q from `flowRate` [m3/s] (T-Flows `FLOW_RATE
  1.67e-6`) or `flowRateMLPerMin` (the paper's flow rates, used as the
  actual flow of the non-expanding parabola);
- `phi` = rho times the **exact** flux of that velocity through every face
  (a fan of triangles per face, the three-edge-midpoint rule on each, exact
  for the quadratic profile, in a canonical order): `max |sum_f phi|` per
  cell is round-off on any mesh (2e-18 of the inflow on the pipe of run/001,
  6e-17 on a 5-degree wedge, also on a wedge whose inner points are moved
  radially and axially, where a face-centre flux misses continuity by
  1e-5; M4.1h), also on a decomposed case (the two sides of a processor
  face are exactly opposite);
- `T`: `uniform`, `profile` (the wall table T(x) on the whole
  cross-section) or `solved` (the default of run/014): the steady
  temperature of the frozen flow, `div(rho cp phi_v T) = laplacian(k(T),
  T)`, with the wall table on the wall patches (Dirichlet at the face
  centres), `inletTemperature` at the inlet and zero gradient elsewhere,
  Picard on k(T).  This is what T-Flows does: it solves its energy
  equation in every iteration (`Main_Pro_Frozen.f90`) with the table as a
  wall condition, density 1, cp 5193 J/(kg K) and the NIST helium
  conductivity of `User_Mod/Types.f90`; its transient from 500 K reaches
  this steady state within 0.1-0.2 s.

The wall table (`format TFlows`: a count line, then rows of x u v w T q_01
q_02; or `columns`) is interpolated linearly in x; `rows count` reads only
as many rows as the count line says, as T-Flows does (66 of the 70 rows of
`wall_x_profile_complete.dat`, x <= 0.65 m), `rows all` reads them all and
names the mismatch; `outOfRange fatal | hold`.  Beyond the rows it reads,
T-Flows assigns the wall no temperature: those faces keep the boundary
value 0 of their allocation, a 0 K Dirichlet wall of its energy equation
(and its `k_He` lookup leaves the table for the cells it cools), which
neither `hold` (324 K) nor `rows all` reproduces (plan section 29, item
5).  On a wedge the sector angle comes from the two wedge patches, exact
to the rounding of their normals (2.5e-13 for the wedge of M4.1c), or
from `sectorAngle`, exact, and the inlet carries the fraction theta/(2 pi)
of Q.  The faceted cross-section is smaller than the circle (0.5 % for the
36 sides of run/001): `rescaleToFacets yes` scales the amplitude so that
the discrete inflow is Q (exactly, to round-off, with `sectorAngle`); `no`
(T-Flows, run/014) keeps U_b = Q/(pi R^2).  The discrete flows of the two
codes still differ: T-Flows takes its face fluxes from the parabola at the
vertex-mean centres of its cells and faces, which on the mesh of run/001
carries 1.0038 Q (a discrete bulk velocity of 0.09311 m/s against the
0.09276 m/s of the exact integral here), so at the same `FLOW_RATE`
run/014's carrier is 0.39 % slower than T-Flows' (plan section 27, item
9).  A carrier with T-Flows' total flow on this mesh is `flowRate
1.67641e-6` with `rescaleToFacets yes`.

The log prints U_b, Q, the discrete inflow, the facet-area ratio, the
discrete bulk velocity (inflow over inlet area) and the continuity check:

    U_b = Q/(pi R^2) = 0.09228776214 m/s (not rescaled to the facets), ...
    discrete inflow through INLET = 1.669948533e-06 m3/s = 0.9999691814 of
      Q sector (analytic profile: 0.9999691814); inlet area 1.800384306e-05
      m2, facet-area ratio (inlet area/(pi R^2 sector)) 0.99493077; discrete
      bulk velocity 0.09275511497 m/s
    Continuity (phi as written): max over the cells of |sum_f phi| =
      3.67127e-24 kg/s = 2.19843e-18 of the inflow 1.66995e-06 kg/s; ...

**The walls.**  The exact flux of the axial profile through a wall face
vanishes only when the face is parallel to the axis (a mesh extruded along
it, as the pipe of run/001 and blockMesh wedges).  Every continuity line is
followed by a `Wall flux` line (the sum over the faces of the `walls` of
|phi| relative to the inflow, and the largest face), and the utility stops
without writing when the sum exceeds 1e-12 of the inflow (a tetrahedral
mesh: 9.5e-4), because the solver assumes a zero flux through the wall of
a paired species; `permeableWalls yes` accepts it with a warning (M4.1i).
On the pipe of run/001 the wall faces carry +-1e-25 kg/s (8.4e-16 of the
inflow in total): the exact flux through faces that gmsh wrote tilted by
1.5e-15; they are kept, as the exact integral, and the solver books them as
transport through WALL.

**The input against the mesh.**  The utility stops without writing when
the analytic profile does not enter the domain through `inlet` (its flux
there is not positive: an axis pointing out of the inlet, or `inlet`
naming the outflow patch), and when `radius` is smaller than the largest
distance of a wall point from the axis by more than 1e-9 relative (the
parabola would be negative near the wall, with backflow through those
faces while U is clipped to 0 there), unless `allowBackflow yes`.  Before
the M4 review, round 2, both wrote a reversed or partly reversed carrier
without a message, and `rescaleToFacets yes` hid them (a negative
amplitude; an inflow of exactly Q again) (M4.1k).

`setFrozenCarrier -check [-parallel]` only checks the continuity (and the
wall flux) of the `phi` the solver reads, the instance that findInstance
finds from the start time (0 for a carrier the solver does not write; the
line names it), e.g. of a case decomposed by `decomposePar`.  decomposePar
writes every patch whose values are all equal as a `uniform` value of
about 10 digits, in a binary file too: processor patches, and an ordinary
patch with a single face on a rank (plan section 15, item 7, and section
26: an exact flux can lose its round-off continuity there; M4.1e and
M4.1d show the check reporting such values).  Run the check after every
decomposition (e.g. scotch on Merlin7), or write the carrier on the
decomposed case: `-parallel` without `-check` writes it on the decomposed
mesh, exactly (M4.1b, M4.1d).  The dictionary
(`system/setFrozenCarrierDict`) is documented in the header of
`setFrozenCarrier.C`; run/014 holds the T-Flows parity version.  `Allwmake`
builds the utility with the solver.

## Numerical controls

The corrector controls shared by all gaseous species are read from the
`speciesTransport` dictionary in `system/controlDict`, for example:

    speciesTransport {
      tolerance          1.0e-7;
      nCorr              47;
      writeFrozenFields  no;      /* optional, default yes */
    }

These controls are read once at start-up; restart to change them.

## Molar concentration output

If a species dictionary contains

    molarMass <value>;

with the molar mass in kg/mol, the solver also writes

    c_<speciesName> = rho * Y_i / M_i

in mol/m3 at every scheduled output time.

For gaseous PbI2,

    molarMass 0.46100894;

corresponds to

    MPbI2 = 207.2 + 2*126.90447 = 461.00894 g/mol.

The concentration field is derived output only.  It is not transported and no
additional equation is solved.  No `0/c_<speciesName>` input file is required.
The field is reconstructed from local `rho` and the current mass fraction at
write time, including boundary values, and is recalculated after a restart.
Existing time directories are not retroactively updated.

For PbI2 the resulting field `c_PbI2_g` has dimensions
`[0 -3 0 0 1 0 0]`; its numerical values are in mol/m3, not kmol/m3.  Select
`c_PbI2_g` in ParaView when molar concentration is desired for post-processing.

## PbI2-He variable diffusivity

For

    diffusivityModel PbI2He;

the solver evaluates a spatial molecular-diffusivity field `D_PbI2_g` in m2/s
from the frozen local `T` and absolute `p`, including boundary values, and then
forms

    rhoD_PbI2_g = rho * D_PbI2_g

in kg/(m s).

Because the carrier state is frozen, `D_PbI2_g` and `rhoD_PbI2_g` are evaluated
once when the solver starts and remain fixed during the transient species run.
Both are written at normal output times when `writeFrozenFields` is `yes`
(never with `no`) and are recalculated when the solver is restarted.

`PbI2HeDiffusivity.H` reproduces the active fit used in Francesco's Python
`diffusivity.py` / `interpolation.py` calculation:

    reducedT = T / 71.11785702515257

    omega = 0.5313
            - 1.2946*exp(-0.9922*reducedT)
            + 1.5591/reducedT
            - 0.1918/reducedT^2

    D = 1.86e-7 * T^1.5 * sqrt(1/4.0026 + 1/461)
        / ((p/pressureScale) * 4.2584682449728355^2 * omega)

The original correlation coefficient `1.86e-3` gives diffusivity in cm2/s;
the factor `1e-4` converting cm2/s to m2/s is already included above, hence
`1.86e-7`.

The collision parameters are derived from `alpha(PbI2)=17.23 Angstrom^3`,
`N=16` and the tabulated helium properties.  The physical interpretation of
`N` still needs clarification; its numerical value is preserved exactly from
the Python implementation.  The spline evaluation in the Python utilities is
not used here because the script overwrites it with the explicit fit above.

The diffusivity calculation preserves the Python script's molecular mass
`461 g/mol` so that the OpenFOAM implementation reproduces the same curve.  The
more precise `0.46100894 kg/mol` is used separately for the derived molar-
concentration output.

## Phase change: `constant/thermochemistryProperties`

Milestones M2 to M5 of `doc/phase-change-plan.md`: the exchange core,
the temperature model, release samples and the mass ledger (M2); the
axial profiles, T_dep, `mDep_` and `exch_` (M3); the paper-mode carrier
(`utilities/setFrozenCarrier`, above), the `distance` interaction layer
and the code-to-code case run/014 of the temperature model (M4); the
Hertz-Knudsen-Schrage model with vapour-pressure tables, inventory samples
and their removal (M5); the GEMS3K backend of its vapour pressure,
`equilibrium GEMS` in mode frozen or local, with the optional speciation
(M9; a build with the GEMS3K bridge, see Building).  The verification of
M7 (the Graetz-Robin problem, the refined meshes and the time-step
studies) is described under Tests.

### Models

| `model` | meaning |
|---|---|
| `mock` (default; also without the file) | the coupled mock of the main branch (tag 0.1, case 011): the thermochemistry at the beginning of every step, the PbI2_g -> PbI2_s precipitation in the cells next to `WALL`, the gas solve with the implicit sink, the solid solve with the realised sink, the integrated sources and the balance of every species, the fields `source_<species>`; unchanged to the last bit.  The phase-change object reads, creates and prints nothing.  The mask of the cells next to the patch `WALL` (and the line `Cells adjacent to WALL`) belongs to this model only |
| `temperature` | Lobresco et al. (2026) Eq. 5 as a wall flux, irreversible deposition |
| `HKS` | Hertz-Knudsen-Schrage (M5): `mdot'' = G (p_eq(T) - p_i)` with p_eq from a table, deposition and evaporation; inventory samples |

**Sign convention:** a positive exchange is net evaporation into the gas.

### Dictionary

    model temperature;

    pairs {                          // one entry per paired GAS species
      PbI2_g { condensed PbI2_s; }   // a listed solid species
    }

    interface {
      patches              (WALL);   // ordinary patches, any name (WALL is type patch)
      layer                faceCells;// faceCells: one element per wall face;
                                     //   distance: one per cell whose centre
                                     //   lies within 'distance' (M4, below)
      // distance          1e-4;     // [m], required with layer distance
      // centreDistance {            // layer distance, optional: also print
      //   axis (1 0 0); origin (0 0 0); radius 2.4e-3;  // the T-Flows test
      // }                           //   R - r <= distance (a comparison only)
      areaPerVolume        cell;     // cell: A_e = |S_f|; layer: uniform A_I/V_I
                                     //   (layer distance: layer, the default)
      areaMultiplier       1;
      interfaceTemperature cell;     // cell (paper, T-Flows) | wall (face T;
                                     //   distance: of the nearest wall face)
      wallResistance       none;     // none | halfCell (series half-cell diffusion;
                                     //   needs areaMultiplier 1, areaPerVolume cell,
                                     //   layer faceCells); the default none keeps
                                     //   paper-mode parity, the physical runs of
                                     //   M8 set halfCell (plan sections 28 and 33)
    }

    temperatureCoeffs {
      A     220;     // [1/s], or depositionVelocity v_d [m/s]; exactly one
      k     8e-3;    // [1/K]
      Tdep  680;     // [K]
    }

    samples {                        // optional
      boat {
        pair       PbI2_g;
        selection  { boat { action use; source cylinder;
                     point1 (0.04 0 0); point2 (0.05 0 0); radius 5e-4; } }
                                     // (run/012; an inventory sample needs
                                     // its boundary between the cell
                                     // centres: run/013 has 4.8e-4)
        // cellZone sample;          // alternative: an existing cellZone
        mode       release;          // release | inventory (model HKS)
        amount     3.16e-5;          // n0 [mol]
        startTime  0;                // [s]
        duration   6;                // t_ev [s]
        // removeTime, areaPerVolume, kineticScale, Ce: mode inventory;
        // a release sample stops with "belongs to mode inventory"
      }
    }

With `model HKS` (M5; section "The HKS law" below):

    model HKS;

    pairs {
      PbI2_g {
        condensed PbI2_s;
        vapourPressure {             // required for HKS
          file          "<constant>/pv_PbI2_phases.csv"; // '#' header names the columns
          phases        ("PbI2(cr)" "PbI2(l)");         // p_eq = their minimum
          units         log10Pa;       // log10Pa | Pa | bar (required)
          interpolation logInverseT;   // logInverseT (default) | linearT (T-Flows)
          outOfRange    extrapolate;   // extrapolate (default) | hold (T-Flows) | fatal
        }
      }
    }

    HKSCoeffs {                      // optional, defaults shown
      accommodation 1;               // sigma, 0 < sigma <= 1
      Ce            1;               // multiplier
      kineticScale  1;               // 1 = SI; 3.16228e-7 = the T-Flows expression
      equilibrium   table;           // table | GEMS (M9; a build with the bridge)
      speciation    none;            // none | lagged (equilibrium GEMS, local)
    }

    samples {
      boat {
        pair          PbI2_g;
        selection     { ... };       // as above
        mode          inventory;     // a reservoir of the condensate
        amount        4.84e-6;       // n0 [mol], filled on a fresh start
        areaPerVolume 4000;          // a_s [1/m] (2/R_sample), required
        kineticScale  1;             // optional, of this sample
        Ce            1;             // optional, of this sample
        removeTime    6;             // [s], optional: default never
        // startTime, duration: mode release; an inventory sample stops
        // with "belongs to mode release"
      }
    }

    balance {                        // optional, defaults shown
      nGuard               3;
      guardTolerance       1e-12;    // relative to the pair inventory
      moleFractionWarning  0.05;     // passive-carrier validity monitor
      writeFile            yes;      // postProcessing/phaseChangeBalance
      reportMatrixResidual no;       // debug: also print gSum(fvMatrix::
                                     // residual()), the explicit
                                     // remainder of the transport matrix
                                     // and the transport identity of the
                                     // cell defects (summed exactly)
      writeExchangeField   no;       // exch_<gas>: realised exchange [kg/(m3 s)]
    }

    profiles {                       // optional: axial profiles and T_dep
      axis     (1 0 0);              // default (1 0 0)
      origin   (0 0 0);              // default (0 0 0)
      bins {                         // one or more sets, one file each
        native { min 0; max 0.69; nBins 424; }  // = the axial cells
        gamma  { min 0; max 0.69; width 0.01; } // the 1 cm gamma-scan bins
      }                              // min/max default to the mesh extent
      times      (7 12 16);          // [s], in addition to the write times
      onset {
        fraction    0.01;                  // of the peak (default)
        search      (fromInlet fromPeak);  // one or both (default both)
        searchFrom  0;                     // [m], fromInlet only (default:
                                           //   the first bin)
        significance 1e-12;                // of the pair inventory (default:
      }                                    //   balance guardTolerance)
      observable wallPlusGas;        // wall | wallPlusGas (default): the
    }                                //   deposit plus the gas outside the
                                     //   cells of every sample

A paired gas needs `molarMass` in `speciesTransportProperties`.  If the
condensate has `molarMass` too, it must be the same (the balance of a pair
is kept in kg).  Outside `model mock` every solid species must be the
condensate of a pair, and the mesh needs no patch named `WALL`: the
interface patches may have any name (the mock mask of the `WALL` cells is
not built; before the final M2 pass a mesh without `WALL` stopped with
"boundary patch WALL was not found", M2.1o).

### Interface elements and the temperature law

`layer faceCells`: every face f of the interface patches is an element e in
its cell c, with the area `A_e = |S_f| * areaMultiplier` (for
`areaPerVolume layer` rescaled to the layer mean `A_I/V_I`).  The start-up
report prints the layer: elements, cells, `A_I`, `V_I` and `A_I/V_I`; for
the pipe mesh `A_I/V_I = 17681.9 1/m`.

The element law is affine in the gas mass ratio of its cell,

    q_e = a_e - g_e Y_c   [kg/s],   a_e >= 0,  g_e >= 0,

and for the temperature model

    a_e = 0,  g_e = A_e rho_c v_d f(T_e),  f = [1 - exp(-k (Tdep - T_e))]^+,

with `v_d = A V_I/A_I` when `A` is given (`v_d = 0.0124421 m/s` on the
pipe).  With `wallResistance halfCell` the kinetic velocity is put in series
with the diffusion of the half cell, `v' = 1/(1/(v_d f) + rho_c d_e/(rhoD)_e)`,
`d_e = 1/deltaCoeffs` of the face; `rhoD = 0` on the face (e.g. `D 0`)
gives `v' = 0` (M2.1n).  It works on ranks without interface faces (M2.4h;
before, the delta coefficients were built inside the element loop, which is
collective, and such a decomposition hung at start-up).  The series law
is an OpenFOAM mixed/Robin condition per physical wall area `|S_f|`, so
`halfCell` needs `areaMultiplier 1` and `areaPerVolume cell` and stops with
a FatalIOError otherwise (M2.6l).  The element law is `g_e = A_e rho v'`
with `v' = 1/(1/v + R)`, `v = v_d f(T_e)` and `R = rho d_e/(rhoD)_e` the
resistance of the half cell.  With a scaled area `A_e = m |S_f|` it is

    g_e = m |S_f| rho/(1/v + R) = |S_f| rho/(1/(m v) + R/m):

the scale divides the half-cell resistance by m as well (the diffusion
through the half cell becomes m times more conductive).  The half cell
belongs to the mesh face, whose diffusive conductance `|S_f| rhoD/d_e` does
not depend on the kinetic area; the Robin condition with the kinetic area
`m |S_f|` is `|S_f| rho/(1/(m v) + R)`.  The two agree only in the kinetic
limit `R << 1/(m v)`, so scaling `A_e` is not Robin-equivalent.

`halfCell` is the second-order discretisation of the Robin wall law and
`none` the first-order one: on the Graetz-Robin problem of M7 (Pe 5, Bi 1
and 50) halfCell converges at order 2.0 with errors of the developed decay
rate of 8e-7 to 2e-4 on the meshes of the study (256 to 32 radial cells)
and at most 6e-5 at the resolution of the production first layer, against
0.3 % (temperature model) to 1.7 % (HKS) for `none` (M7.1, plan sections
28 and 33).  It is
therefore the wall law of the physical runs (M8: their cases set
`wallResistance halfCell`), while the code default stays `none`, the law
of the paper and of T-Flows, which the code-to-code cases (run/014, M4long,
M6) reproduce; every case of `run/` up to 015 keeps `none`.

`areaMultiplier` scales `A_e`, and hence the rate, of every model (plan E2):
the temperature model applies `areaMultiplier*v_d*A_I/V_I` (times f), i.e.
`areaMultiplier*A` when `A` is given, not `A`.  The start-up line prints it,

    Temperature model: v_d = 0.0124421 m/s (from A = 220 1/s on this layer),
      k = 0.008 1/K, Tdep = 680 K; applied rate areaMultiplier*v_d*A_I/V_I =
      220 1/s (layer mean, times f(T)); 8280 of 15264 elements below Tdep

With `wallResistance halfCell` this is the kinetic rate only: the line says
`kinetic rate ... (layer mean, times f(T); the half-cell resistance reduces
it: pair line)`, and the rate applied with the half cell is printed per
pair (M2.1d),

    Pair PbI2_g, wallResistance halfCell: applied rate sum_e A_e v'_e/V_I =
      ... 1/s (layer mean with f(T) and the half-cell resistance; kinetic
      sum_e A_e v_d f_e/V_I = ... 1/s)

Keep `areaMultiplier 1` for the temperature model: the factor 200 of the
T-Flows parity layer (plan section 11) belongs to its HKS source.

### The interaction layer `distance` (M4)

    interface {
      patches        (WALL);
      layer          distance;
      distance       1e-4;               // [m], the paper's 0.1 mm (T-Flows DELTA)
      centreDistance {                   // optional, printed only
        axis (1 0 0); origin (0 0 0); radius 2.4e-3;
      }
    }

**Elements.**  One element per cell whose centre lies within `distance`
of the interface patches: the exact distance of the cell centre from the
nearest face of any interface patch on any rank (`face::nearestPoint`),
and f(e) that face, identified by its global index.  Every rank receives
the interface faces within reach of its cells (twice the distance around
their bounding box; PstreamBuffers) and searches them with an octree; a
tie (a centre nearest to an edge or a point shared by several faces) goes
to the face with the smallest centre (x, then y, then z), so the layer and
its nearest faces depend neither on the order of the interface patches nor
on the decomposition (M4.13, M4.14).  The first M4 round used
patchDataWave, whose near-wall correction let the last of several patches
overwrite the distance of a corner cell and saw only the wall faces of the
cell's own rank (plan section 26).  Every element has the mean area per
volume of the layer (paper Eq. 7: A_I/V_I applies uniformly to the
interaction region):

    A_e = V_c (A_I/V_I) areaMultiplier,  A_I = sum |S_f| (interface faces),
                                         V_I = sum V_c (layer cells), all ranks

so `areaPerVolume` is `layer` (the default here; `cell` stops: a cell of
the second layer has no wall face of its own), and the temperature model
applies `areaMultiplier A f(T)` in every layer cell.  `wallResistance
halfCell` stops as well (no half cell of a wall face belongs to a cell of
the second layer).  `interfaceTemperature wall` uses the temperature of
the nearest wall face.  Both models work with the layer (HKS: M4.10).

The start-up lines give both counts of plan section 12, item 10: the
layer of the exact distance with the distances of its cells and of the
nearest cell outside within twice the distance (the margin of the rule),
and with `centreDistance` the T-Flows test `R - r <= distance` about the
given axis (a comparison only).  On the pipe of run/001 (run/014):

    Interface (layer distance): patches (WALL), distance 0.0001 m: 30528
      elements in 30528 cells; A_I = 0.0103918 m2, V_I = 1.21773e-06 m3,
      A_I/V_I = 8533.68 1/m
      layer distance (the exact distance of the cell centres from the
      faces of the interface patches, all ranks): selected cells at
      2.85045e-05 to 8.85789e-05 m, the nearest cell outside at
      0.000154656 m; 0 elements with their nearest wall face on another
      rank, 0 interface faces with elements of two or more other ranks
      centre-distance rule R - r <= 0.0001 m (...): 30528 cells (largest
      R - r selected 9.77116e-05 m, smallest outside 0.000163789 m); 0
      cells only in this rule, 0 only in the layer

two cell layers, against one (15264 cells, `A_I/V_I = 17681.9 1/m`) for
`layer faceCells`.  The paper's text ("all control volumes lying entirely
within" 0.1 mm) describes the one layer, T-Flows' `R - r <= DELTA` the two
(`manuscript/notes-for-authors.md`, item C); run/014 switches between them.

**Why not patchDataWave.**  OpenFOAM's wall distance with the nearest face
(`patchDataWave<wallPointData<label>>`) is wrong for this purpose in three
ways (v2412 and v2606; plan sections 25 and 26): with the combined wall
patch (the default) `correctBoundaryCells` stores a label of the combined
patch for the cells next to the wall (on the test channel the first row
of one wall got the faces of the other wall); per patch, every interface
patch overwrites the distance of a cell with faces on several patches while
the face of the first is kept (a corner cell dropped or mapped to the
farther face, depending on the order of the patches in the boundary
file); and either way only the wall faces of the cell's own rank are seen
(a cell touching the wall by an edge of a face on another rank kept the
larger wave distance: the layer depended on the decomposition).

**mDep_ through the nearest face** (plan section 12, B3).  A
`mapDistributeBase` built from the global indices of the nearest faces
defines which faces every rank needs; the mass `m''_e A_e` of every
element is pushed to its face, so that on every interface face

    mDep_ = sum_{e -> f} m''_e A_e/|S_f|,

and `areaIntegrate` of `mDep_` over the interface patches is the wall
inventory to round-off in serial and in parallel, also for the elements
whose nearest face lies on another rank (the start-up line counts them:
15264 of the pipe's elements when its first cell layer is on a rank of
its own, M4.7c; 200 on the channel and 161 on the pipe of M4.7h, whose
maps run both ways between two ranks).  The push sends the sums of the
remote slots in `PstreamBuffers` to the rank of each face, which adds them
to its own: its own elements first, then the other ranks in ascending
order, every contribution once, so a face with the elements of two or
more other ranks (counted at start-up by the same push) gets the same sum
in every run of a decomposition (M4.7g).  Neither branch of a reverse
`mapDistributeBase::distribute` with `plusEqOp` does both: with the
default nonBlocking communication it adds the contributions in the order
in which the receives complete, and `mDep_` could differ in its last bit
between identical runs (review of M4, round 2); its scheduled branch,
used in round 2, follows `mapDistributeBase::schedule()`, which lists
directed pairs, and exchanges both ways for each, so between two ranks
that map elements to each other's faces every contribution was added
twice (review of M4, round 3: `mDep_` +50 % of the wall inventory on the
channel of M4.7h, +7.8e-4 on the pipe split by hierarchical (1 1 4),
+0.92 % on the tetrahedral pipe of M4.14b, and the start-up count gave 200,
13 and 14 faces instead of 0; plan section 29, item 1).  `mDep_` is output
only: the ledger, the inventories, the deposits, the profiles and T_dep
were never affected (with the build before the fix and the build after
it, the fields, the per-rank state and the profiles of those runs are
identical but for `mDep_` itself; their balance files differ only in the
columns of the ledger's booking, `solverDefect` and `closure`, which the
flux form of the defect changed in the same round, and in the monitor
`evaporatedAboveWarning`).  The pull of the face values (the temperature
of the nearest face for `interfaceTemperature wall`) is tested on 4 ranks
with model HKS (M4.10, second line).  `mw_` is a cell field here: the
deposit m''_e of every layer cell itself (0 on the patches), equal to the
deposits of `uniform/phaseChange/` bit for bit (M4.7d, M4.7f).

**State and restart.**  The deposits and regimes are kept in the per-rank
lists of `uniform/phaseChange/`, like those of `layer faceCells`, behind a
fingerprint that hashes the rule and the distance first, then the cell and
the exact centre of every element (the faceCells fingerprint and format
are unchanged).  The layout records `layer` and `distance` (exact
decimal): a restart with the other rule or another distance stops with a
FatalIOError that names the change (M4.8).  A binary restart repeats the
continuous run bit for bit, in serial and on 4 ranks (M4.7d); a restart
onto another decomposition re-reads the deposits from the cells of `mw_`
(exact in binary), with a warning, and finds the ledger's inventory
exactly (M4.7e): with `areaPerVolume layer` the element areas follow
A_I/V_I, and A_I and V_I are compensated sums, independent of the
decomposition (plain sums over the 30528 cells differed by up to 2e-14
between decompositions, beyond the round-off bound of an exact restart).

**Profiles** (decision of plan section 25).  The deposit of a layer cell is
binned by the axial overlap of its own cell, like the gas of that cell,
not by the extent of its nearest wall face: bins + outside = inventory
exactly and locally on every rank, and the line density is T-Flows' sum
over the cells of an axial slice.  On a mesh extruded along the axis (the
pipe) a layer cell and its nearest face have the same axial extent, so
the choice does not matter there.  The `Twall` and `wallArea` columns come
from all faces of the interface patches, as for `layer faceCells` (equal
bit for bit, M4.9), so T_dep is the wall temperature at the onset for both
rules.

### The HKS law and the vapour-pressure tables

`model HKS` (M5; plan E2, SI) applies the Hertz-Knudsen-Schrage law to
every interface element, in both directions:

    q_e = A_e G_e (p_eq(T_e) - beta_c Y_c),   beta_c = rho_c R T_c/M,
    G_e = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T_e))  [s/m],

i.e. `a_e = A_e G_e p_eq(T_e)` and `g_e = A_e G_e beta_c`: the partial
pressure `p_i = beta_c Y_c` uses the carrier's actual density (with rho = 1
in paper mode), and the equilibrium mass ratio of an element is
`Y_eq = p_eq/beta_c`.  For PbI2 at 700 K `G R T/M = 89.65 m/s`; on the
closed test box (57.46 um cells) `lambda = sum_e g_e/(rho V) = 1.56e6 1/s`,
`lambda dt = 6241` at 4 ms (M5.1a).  `kineticScale 1e-5/sqrt(1000) =
3.16228e-7` reproduces the T-Flows expression (p in bar, M in g/mol;
M5.0).  With `wallResistance halfCell` (faceCells, `areaMultiplier 1`,
`areaPerVolume cell`) the conductance is put in series with the half cell,
`G'_e = G_e/(1 + G_e beta_c d_e/(rhoD)_e)`, the Robin condition of the
face-value law (M5.0, M5.1e).

Bounds of a wall element in a step (plan E3): a **bare wall** (no deposit)
only deposits (`upper = 0`: an element that would evaporate is NONE, never
CAPPED; M5.1d); an element with a deposit also evaporates, at most its
deposit, `upper = m''_e A_e/dt`, and becomes CAPPED when the law would take
more (the deposit is then set to exactly 0 and the round-off residue is
booked as `clamped`; M5.1b/c).  The predictor, the settled-regime exit, the
guard, the projection and the complementarity check are those of the
temperature model.

The start-up prints, per pair, the table (file, phases, nodes, range,
options, and the temperature where the stable phase changes) and the HKS
ranges on the interface:

    Pair PbI2_g: p_eq from "<constant>/pv_PbI2_phases.csv", phases
      2(PbI2(cr) PbI2(l)) (the minimum), 198 nodes 270-1250 K; units log10Pa,
      interpolation logInverseT, outOfRange extrapolate; stable phase
      PbI2(cr) -> PbI2(l) at 683 K
    Pair PbI2_g, HKS on the interface: T_e 300.825-999.175 K, p_eq ... Pa,
      Y_eq = p_eq/beta_c ...; G_e beta_c/rho_c 58.77-107.1 m/s; lambda =
      sum_e g_e/(rho_c V_c) 1.04e+06-1.89e+06 1/s, layer mean 1.42e+06 1/s;
      0 elements outside the table (none); 0 elements with p_eq > p (the
      condensate would boil; not capped)

**Tables** (`vapourPressureTable.H`, standalone, tested by
`tests/testVapourPressureTable.C`): a csv file whose `#` line right before
the first data row names the columns, column 0 is T [K]; `phases` selects
columns, each the vapour pressure over one condensed phase over the whole
range, and `p_eq` is their minimum (the stable phase), so the kink of a
phase transition is exact.  `logInverseT`: log10 p linear in 1/T
(Clausius-Clapeyron), `linearT`: p linear in T (T-Flows).  The node values
are returned exactly in either mode and unit.  Outside the table
`extrapolate` continues the end interval, `hold` keeps the end value, and
`fatal` stops the start-up with the range of the element temperatures
outside it (M5.5).

`thermochemistry/systems/make_pv_table.py` writes the two tables of the
repository from local sources (no download): `data/pv_PbI2_phases.csv`,
Gurvich (1991) through the NASA-9 records of `thermo.inp` (PbI2(g) over
PbI2(cr) and over PbI2(l), log10 Pa, every 5 K from 270 to 1250 K plus the
crossing of the two phases, T_m = 683.000002 K), and
`data/pv_PbI2_TFlows.csv`, the `GEMS_VAPOR_PRESSURE` column [bar] of the
T-Flows reference implementation at its 95 nodes 283.15 + 10 k K, copied
digit for digit (use `interpolation linearT; outOfRange hold; units bar;`).
`make_pv_table.py --check` compares the 5 K table, interpolated as the
solver does, with NASA-9 on a 0.1 K grid: 4.4e-5 decades at most (M5.4b).
The T-Flows source is read from a fixed commit (`TFLOWS_REVISION`, recorded
in the header of the csv file), not from the moving branch; `--sources`
checks that both sources are on the machine (M5.4b is NOTRUN otherwise,
while M5.4a still tests the repository's tables).

### p_eq from GEMS3K: `equilibrium GEMS` (M9)

With a solver built with the GEMS3K bridge (see Building), `HKSCoeffs {
equilibrium GEMS; }` takes the equilibrium vapour pressure of the HKS law
from a Gibbs-energy minimisation of the gas (`gemsEquilibrium.{H,C}`; plan
E9 and section 32).  The law keeps its form; the pair's table stays
required, as the check and the fallback of every GEMS3K result.  A build
without the bridge stops at `equilibrium GEMS` with a FatalIOError and the
rebuild hint (M9.1).

    HKSCoeffs  { equilibrium GEMS; speciation none; }   // or lagged (local)

    GEMSCoeffs {
      system            "<constant>/gems/PbI2He-dat.lst";  // required
      mode              local;        // frozen | local (required)
      carrier           "He(g)";      // the carrier's gas species (required)
      referenceComposition { He 1; Pb 1e-3; I 2.000002e-3; }  // frozen [mol]
      referencePressure 101325;       // frozen [Pa]
      minMoleFraction   1e-6;         // local: below it the table
      minTemperature    500;          // [K]: below it the table
      maxLog10Deviation 0.01;         // further from the table: rejected
      updateInterval    1;            // local: steps between updates
      allowNonPhysicalCarrier no;     // local refuses a carrier that is no gas
      logDirectory      "gemsLog";    // per rank: processorN/gemsLog
      logLevel          4;            // GEMS3K logging: errors only
      writeEquilibriumField no;       // pEq_<gas>, gemsSource_<gas>, chi_<gas>
    }

    pairs { PbI2_g { condensed PbI2_s; vapourPressure { ... }
      gems {
        gas          "PbI2(g)";              // the gas species in the system
        condensates  ("PbI2(cr)" "PbI2(l)"); // p_eq over the most stable
        elements     { Pb 1; I 2; }          // the formula (checked)
        excess       { I 1e-6; }             // optional, see below
      } } }

The values shown are the defaults (the plan's), but for the required keys.
`constant/gems` is a link to `thermochemistry/systems/PbI2He` in run/015;
the system file names its DCH, IPM and DBR files relative to itself.  A
key of the other mode (`referenceComposition` and `referencePressure` in
mode local; `minMoleFraction`, `updateInterval` and
`allowNonPhysicalCarrier` in mode frozen) is named in a start-up warning
and not read, so an invalid value of it cannot stop the run.  Any other
key, e.g. a misspelled `minMolFraction`, stops the run with the list of
the valid keys, in `GEMSCoeffs` and in a pair's `gems` block (plan section
35, items G and H; M9.7, M9.11).

**What a call computes.**  Every condensed species of the system is
suppressed (metastability bounds 0 and 0), so GEMS3K returns the gas at the
bulk composition, and the activity of a suppressed pure condensate is its
saturation ratio Omega: `p_eq = p_g 10^(-max log10 Omega)` over the pair's
condensates.  For the ideal PbI2-He gas this is the table's minimum over
the phases, independent of the composition, to the G0 interpolation of
GEMS3K: within 1e-5 to 6e-5 decades in cold calls at 500-1250 K.

**mode frozen** (a consistency check, plan section 9 item 7): at start-up
the master evaluates p_eq at the nodes of the pair's table within
[max(minTemperature, 270 K), 1250 K] (the grid of the system), cold, at the
reference composition and pressure; a failed node, or one off the table by
more than `maxLog10Deviation`, keeps the table's value.  This GEMS table is
broadcast to every rank and interpolated as the pair's table; an element
below `minTemperature` or outside it takes the pair's table.  A call that
returns `GEMSB_ERR_FATAL` re-creates the engine as mode local does (M9.5;
plan section 35, item E).  Nothing changes per step, and a restart is
exact (M9.2b).  The start-up says, e.g.,

    Pair PbI2_g, GEMS frozen: p_eq at the 152 nodes of the table in
      500-1250 K: 152 from GEMS3K (...), 0 failed, 0 rejected (taken from
      the table), engine re-created 0 times; max |dlog10 p_eq - table|
      5.68e-05 at the GEMS nodes; ...
    Pair PbI2_g, GEMS frozen on the elements: 206 from the GEMS table (T_e
      505-997 K; max |dlog10 p_eq/p_eq,table| 5.53e-05), from the table 34
      below minTemperature and 0 outside the GEMS table

**mode local**: one engine per rank for the whole run.  At the steps whose
time index n has (n - 1) divisible by `updateInterval` (the absolute index,
restored on a restart: a run that starts at index 0 updates at its first
step and then every `updateInterval` steps; one that starts from a time
whose `uniform/time` holds an index k, without a lagged state, at its
first step k + 1 and then at the schedule's next, so its first interval
may be shorter; plan section 35, item Y), before the correctors, every
element (wall and sample; a
sample element for its own pair) is evaluated at its T_e, the pressure of
its cell and the composition of its cell at the start of the step: the
carrier `p/(R T)` [mol/m3] (as the mole-fraction monitor) plus
`c = rho Y/M` of every pair, times the formulas.  The step then uses one
law throughout.  A pair takes the GEMS3K p_eq unless, in this order,
`T_e < minTemperature`, T_e or p lies outside
the grid of the system, `x = c/(c + p/(R T)) < minMoleFraction`, the call
failed (a status other than OK and OK_RETRIED, an element balance off by
more than 1e-6, or a converged call that gives the pair a partial pressure
p_g <= 0 or a log10 Omega that is not finite: each counted as a failure
and written to the rank's log, plan section 35, item F), or the result
lies further than `maxLog10Deviation` from the table (rejected): then the
table.  Every element keeps its own GEMS3K state for a warm start; a
failure or a rejection makes its next call cold, and `GEMSB_ERR_FATAL`
re-creates the engine (destroyed, created again, every condensed species
suppressed again, every warm state invalidated; M9.5).
A step without an update keeps p_eq (and chi) of the last one: this lagged
state, with the counts of the update's step line, is written with the
per-rank state (`uniform/phaseChange/phaseChangeEquilibrium`), so that a
restart onto the same elements continues the schedule of the time index
and repeats the continuous run bit for bit up to its next update (M9.6;
plan section 35, item C).  The warm states are not written: that update
is cold and follows the continuous run to the tolerance of GEMS3K (5.6e-8
of the inventory in M9.9; plan section 8, risk 13).  Without the lagged
state (another decomposition, or a run in mode frozen or with equilibrium
table before) the first step updates too, with a warning when the
schedule had none there.  Every step has a line per pair after the
exchange line (deterministic, so excerpts compare):

    Y_PbI2_g GEMS local: update (calls 110: warm 108, retried 0, cold 2;
      2.4 iterations per call; failures 0, engine re-created 0); elements:
      GEMS 110, table 130 (below minTemperature 34, below minMoleFraction 96,
      outside the GEMS grid 0, rejected 0, failed 0); max |dlog10
      p_eq/p_eq,table| 8.8e-05

and the end of the run a summary with the mean time of a warm call (its
counts are 64-bit integers: a production run exceeds a 32-bit label; plan
section 35, item A).  Every rank writes `gemsLog/gemsEquilibrium.log` below
its case directory (the case, or `processorN`): the engine, every update
with its counts and timing on the rank, every failure (element, T, p, x,
status, the message of GEMS3K), every re-creation; GEMS3K's own
`ipmlog.txt` lies next to it.  In a parallel run a `logDirectory` outside
`processorN` (an absolute path, `<case>/...`, `$FOAM_CASE/...`, or `..`
out of it) gets a subdirectory `processor<N>` per rank, so no two ranks
share a log (M9.12; plan section 35, item D).  The decision is the same
on every rank (`../processor0/gemsLog` lies inside `processor0` for rank
0 alone, and every rank takes its subdirectory), so the start-up line
names the directory of every rank, e.g. `logs in
/path/gemsLog/processor<N>` (plan section 35, item Y).

Mode local evaluates the cell's gas from the carrier's p and T, so it
refuses (FatalIOError) a carrier that is not a gas at its p and T, `max
|rho - p W/(R T)|/rho` (the carrier line's `max |T - p W/(R rho)|/T`) above
1 %, i.e. paper mode, unless `allowNonPhysicalCarrier yes` (M9.3e).  Mode
frozen needs no composition and runs on any carrier.

**The key `excess`.**  At the exact formula (b_I = 2 b_Pb) the
interior-point method of GEMS3K occasionally stalls (thermochemistry
README, pitfalls): 1 of 2364 warm calls inside the guards in the scan of
plan section 32, item 3, none with an iodine excess of 1e-6, which the
prototypes and the reference composition use as well.  `excess { I 1e-6; }`
adds `nu_I 1e-6 c` iodine to the bulk; it changes p_eq of an ideal gas not
at all and the speciation by about 1e-6.  The default is none (the exact
formula); run/015 and the tests set it.

**speciation lagged** (a separate modelling change, off by default; mode
local only; M9.4): GEMS3K also gives the speciation of the gas, chi = n_g/c,
the fraction of the pair's formula units present as the gas species (the
rest is PbI, I, Pb, I2).  The law then counts only those molecules in the
partial pressure, `q = A G (p_eq - chi beta Y)`: `g_e` is multiplied by chi
(and with `wallResistance halfCell` the half cell enters with chi beta), so
the equilibrium mass ratio becomes `p_eq/(chi beta)`.  chi comes from the
last update ('lagged') and is clipped to [1e-12, 1]; an element on the table
has chi = 1.  In the ramp channel chi is 0.46 to 1 (the dilute gas at the
hot end dissociates).

**Cost.**  A warm call took 42-45 us in a scan of the bridge on the idle
workstation, and 60-70 us in the solver with most other cores busy (M9.3b);
most of it is GEMS3K's re-interpolation of G0(T) for the new temperature of
every element (consecutive elements lie at different T; the bridge clamps
Ttol to 1e-3 K against the Ttol trap of thermochemistry/README.md).  Mode
local adds the calls of the elements with gas above `minMoleFraction` and T
above `minTemperature`: in run/015 about 700 calls per step during the
first steps, 0.05 s of a 2.3 s step.

**A known limit: the gas that is not dilute.**  At high PbI2 mole fractions
GEMS3K does not converge: on the channel from 1300 to 600 K (the boat at
the hot end, 10 steps) 37 calls failed with ERR_NOCONV, at x = 0.05 to 0.73
and 1170-1240 K, after 6 to 7452 IPM iterations (mean 2043), and even the
converged warm calls there took about 4 ms (1000-2000 iterations) instead
of about 50 us: 0.52 s per step against 0.006 s with the table.  Every
failed call falls back to the table and is counted and logged, so the
results stay correct; but such a gas lies outside the passive-carrier
(dilute) assumption of the solver anyway, which the mole-fraction monitor
flags (cells above `moleFractionWarning`, 0.05).  No guard skips GEMS3K
there (plan section 35, item M: a `maxMoleFraction` guard is left for the
qualification of the guards with Pb-Bi-I in M10).

**Output** (`writeEquilibriumField yes`): `pEq_<gas>` [Pa], the p_eq of the
law, `gemsSource_<gas>`, where it came from (0 GEMS3K, 1 below
minTemperature, 2 below minMoleFraction, 3 outside the GEMS grid or table,
4 rejected, 5 failed; -1 elsewhere, or the table before the first update),
and with speciation lagged `chi_<gas>`: the value of the element of each
interface face, the mean (the largest source) of the elements of a cell,
refreshed every step and written with the fields, never read.

### Inventory samples and their removal

A sample with `mode inventory` (model HKS only) holds a **reservoir** of the
condensate in its cells, `Cs` [kg/m3]:

- on a fresh start `Cs = n0 M/V_Z` in every cell of the sample (part of the
  ledger's `initial`);
- every cell is a sample element of the interfacial area `A_e =
  areaPerVolume V_c` at the cell temperature, with the sample's own
  `kineticScale` and `Ce`: `Su = a_s G_s p_eq(T_s)`, `Sp = a_s G_s beta_s`
  per volume.  It only evaporates (`lower = 0`: NONE when the gas is above
  `Y_eq`, no deposition onto the sample), at most its reservoir (`upper =
  Cs V/dt`, CAPPED sets `Cs` to exactly 0); `Cs = 0` gives NONE;
- `Cs^{n+1} = Cs^n - dt q_e/V_c`;
- **removal** at the start of the first step whose start `t^n`, on the
  ledger's exact clock, satisfies `t^n >= removeTime - dt/2`: the
  reservoir `sum Cs V` (compensated over the ranks) is booked as `removed`,
  `Cs := 0`, once.  The log says

      Sample boat removed at the start of the step at t = 0.25 s (removeTime
        0.25 s): its reservoir of 2.08194e-06 kg (4.51606e-06 mol) of PbI2_s
        is booked as removed

  The ledger counts the removals of a pair (`REMOVAL_DONE`), the layout
  records the sample as removed with the clock of the removal, and a
  restart after it never removes again (M5.2a/c/e).  Without `removeTime`
  the sample is never removed.  A fresh start at or after `removeTime -
  dt/2` stops: the reservoir would be filled and removed in the first step
  without evaporating (M5.5).

The exchange line adds the regimes of the sample elements, e.g. `implicit
9360, capped 0, none 5904 elements (sample: implicit 132, capped 0, none
12)`; the balance line and file hold `sample = sum Cs V`; `Y_<condensate>`
(derived) includes `Cs/rho` in the sample cells; `Cs_<condensate>` [kg/m3]
is written with the fields (output), and the profiles bin the reservoir as
the sample line density.

**Round-off of the stiff exchange.**  With `lambda dt` of order 1e3 to 6e3
the diagonal and the source of an exchange cell are that many times its
storage term, and their rounding when the matrix is assembled is `eps
lambda dt` of the gas in that cell per step, while the condensate changes
by `dt q_e` with `q_e = a_e - g_e Y`, a different rounding of the same
stiff quantity.  The ledger therefore books as `solverDefect` the defect of
the **transport-only** matrix plus the realised flows `q_e` (exactly those
that update the condensate) plus the release (`matrixDefect.H`): equal to
the defect of the assembled matrix in exact arithmetic, but its closure
stays at the round-off of the storage and transport terms (and since the
M4 review, round 2, the transport-only defect is formed in flux form, so
that their round-off does not accumulate either: "Mass ledger" below).
With the
defect of the assembled matrix the closure of the closed box at 700 K (all
gas in stiff cells near equilibrium) grew by 6.7e-14 of the inventory per
step, 2.7e-12 after 50 steps; now it is 2e-16 after 3 and 4e-16 after 50
steps (M5.1a, M5.1f), and about 2e-16 in the channel and 4e-16 on the pipe
(M5.2a, M5.3a).

Near equilibrium the exchange is still not exact: at a fixed point the
linear solver no longer iterates (its residual is round-off), but `q_e =
fl(a_e - g_e Y)` is not 0, so the condensate changes by about `eps lambda
dt Y_eq` of the cell's gas per step while the gas does not.  This is a real
imbalance of double precision, and the ledger books it as `solverDefect`
(in the 50-step box the deposit drifts by 4.3e-11, `solverDefect` 1.5e-11
of the reference; M5.1f).  At M8 scale (lambda dt about 350, 1.4e5 steps)
it may matter wherever gas sits near `Y_eq` for many steps; see plan
section 20, item 9.

### Solving a paired species

`solvePairedGasSpecies.H` replaces `solveGasSpecies.H` for a paired gas
(the latter is unchanged for all other species):

- The exchange is a linearised source, `fvm::Sp(Sp, Y)` plus `Su`, with
  `Su >= 0` and `Sp >= 0` (never `fvm::SuSp`): it adds `V Sp >= 0` to the
  diagonal and `V Su >= 0` to the source, so it preserves the diagonal
  contribution of the transport matrix and never adds a negative source.
  Whether the whole matrix is an M-matrix depends on the transport terms:
  `Gauss upwind` without a non-orthogonal correction has non-positive
  off-diagonal coefficients; `Gauss SuperBee` weights the downwind cell
  wherever its (lagged) limiter is positive, so with SuperBee the transport
  matrix itself is not an M-matrix.
- At every outer corrector a local implicit predictor chooses the regime
  of every element from the transport-only matrix (`fvMatrix::A()` and
  `H()`): IMPLICIT (the law), CAPPED (the deposit is exhausted: `q = upper`)
  or NONE.  The bounds of a step are `lower <= q <= upper`; for the
  irreversible temperature model `upper = 0`, so an element whose gas turns
  negative (a negative predictor) becomes NONE and its deposit is kept; an
  element becomes CAPPED only when `upper > 0` (review item B2).  The
  predictor starts from the regimes of the previous corrector (of the
  previous step at the first corrector): an element never switches
  directly between CAPPED and NONE, it passes IMPLICIT (see "A block of
  sample cells" below).
- A paired equation is never relaxed; a relaxation factor in `fvSolution`
  is ignored with a start-up warning.
- The outer loop is converged when the initial residual is below
  `tolerance` and the regimes are settled (see below).  Then one more solve
  with `Y_<gas>Final` follows if that solver exists; it must have `relTol
  0`, and it runs its own predictor (the "final predictor").  Without it
  the last `Y_<gas>` solve is the final one, and `Y_<gas>` must have
  `relTol 0` (FatalIOError otherwise).  Both are checked at start-up and
  again at every step, because `fvSolution` may be edited during a run
  (`runTimeModifiable`).  A separate Final solve exists when `Y_<gas>Final`
  selects its own `solvers` entry, as OpenFOAM selects it: the literal
  keyword, else the last matching pattern.  A single `"Y_.*"` serves both
  names, so there is no separate Final solve and `"Y_.*"` needs `relTol 0`;
  a pattern `"Y_.*Final"` must follow `"Y_.*"` in `fvSolution` to be used
  (M2.6j).  The start-up line names a pattern entry, e.g. `final solve with
  Y_PbI2_gFinal (fvSolution entry "Y_.*Final")`.
- After the final solve a guard checks the bounded law at the solved Y:
  the realised flow of every element against its admissible flow
  `min(max(a - g Y, lower), upper)`.  If the violations exceed
  `guardTolerance` times the inventory in total, violating IMPLICIT
  elements are switched to CAPPED or NONE, and violating CAPPED or NONE
  elements (e.g. a sample cell CAPPED in supersaturated gas, or NONE in
  undersaturated gas, a bare wall NONE in supersaturated gas) back to
  IMPLICIT, the largest misses first, until the misses left add up to at
  most that tolerance; the equation is then re-solved with the final
  solver, at most `nGuard` times.  The guard thus acts on the total that
  triggers it (review of M5, round 3: it switched back only elements whose
  own miss exceeded the tolerance, and accepted two smaller ones of 1.43
  times the tolerance together without a pass, M5.2b).
  A guard re-solve solves the **same** transport matrix as the solve it
  corrects (the transport-only copy of the recorded solve), with only the
  exchange coefficients changed; it never re-assembles `ddt + div -
  laplacian`.  Switching a violating element of a fixed linear system to
  its bound leaves its implicit flow on the same side of the bound, so the
  switches are monotone.  A re-assembled matrix re-evaluates the SuperBee
  limiter (and the non-orthogonal correction) at the solved Y, so a switch
  decided at the solution of one matrix was applied to another: next to a
  colder wall a bare wall element at Y_eq cycled IMPLICIT, NONE, IMPLICIT
  until `nGuard` was used up, and the step was accepted with a miss up to
  1e6 times the tolerance, or with projections, depending on `nCorr`
  (review of M5, round 4; M5.2p, and the 2 projected steps of M5.2b).
  More transport convergence comes from `nCorr` or the outer tolerance,
  never from the guard.
  What is left of the bound violations (round-off) is projected in the
  same cell, conservatively.  Y is never clipped.
- The status of the outer loop is part of the exchange line; paired species
  print no per-step warning, and the end of the run counts the steps
  without outer convergence.

Recommended: loose correctors plus a tight final solve, e.g.

    Y_PbI2_g      { solver GAMG; smoother GaussSeidel; tolerance 1e-8;  relTol 0.1; }
    Y_PbI2_gFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0;   }

and for model HKS a final tolerance of 1e-14 (run/013, run/015; plan
section 28, item 8; presumably its stiff wall cells dominate the
normalisation of the residual, so that the same normalised tolerance
leaves a larger absolute defect): with 1e-12 the booked solverDefect of
run/013 grew by about -2e-18
kg per step once the boat was empty, to 2e-8 to 6e-8 of the inventory
after 16 s (the residual of the final solve, of one sign), and 1e-14 cuts
it a hundred times for about 20 % more time per step (plan section 33).

An outer tolerance below the `tolerance` of the `Y_<gas>` solver is rarely
reached: the initial residual of a corrector is about the final residual of
the previous solve.  The start-up warns about it; the loop then runs all
`nCorr + 1` correctors in every step (the M2 pipe test does this on
purpose, so that its serial and parallel runs take the same path).
Serial and parallel runs agree to about 1e-9 only when the lagged terms
(the SuperBee limiter, the non-orthogonal correction) are converged, e.g.
with enough correctors as in the pipe test (21 per step).  A fixed number
of correctors alone is not enough: the final solve converges only the
linear system of the last iterate, whose lagged terms differ between
decompositions.  On the test channel with 4 correctors per step, serial and
4 ranks differ by 7e-8 (simple 4 1 1; 3.4e-7 in the wall inventory, M2.4f)
and 2e-6 (simple 1 4 1) in the gas inventory; with `Gauss upwind` (a linear
matrix) by 1e-16 (4 1 1) and 6e-12 (1 4 1, the tolerance of the final
solve).  With run/012's own outer loop (`nCorr 10`, `tolerance 1e-5`),
which converges, the pipe on 4 ranks agrees with serial to 4e-11 (gas) and
1.4e-10 (wall) (M2.4e).  The closure is at round-off (about 2e-15) in all
of them.  Plan section 16, item 7 records this condition of the 1e-9
agreement for sign-off.

**Settled regimes (plan section 14, item 1, signed off by the researcher
on 2026-09-25).**  The plan required zero regime changes at the last
corrector.  In cells without gas, Y is solver noise of order 1e-40 whose
sign flips elements between IMPLICIT and NONE at every corrector (run/012:
about a thousand elements carrying 1e-37 kg); with a strict count the
outer loop never converges.  The regimes are therefore settled when the
flow changed by the regime changes, times dt, is at most `guardTolerance`
times the inventory, the threshold at which the guard accepts a violation.

**Significant events.**  Regime changes, projections and complementarity
violations are counted in total and as *significant*: an element whose
own mass (changed, projected or missed flow times dt) exceeds
`guardTolerance` times the inventory of the step (held + released + the
convective inflow through the patches).  Noise-level elements cannot hide
a real event, so "0 significant" is a meaningful check.  The projected
and the missed mass of a step are also compared in total: above the
tolerance, the exchange line names it (`projected 2 (3.1e-18 kg >
tolerance 2.2e-18 kg, 0 significant)`) and the end-of-run summary counts
the steps (`...; 1 step above the tolerance in total`); with passes left
the guard never leaves that much (its passes solve one fixed linear
system, above).  The exchange line reports the regime
changes of the corrector that decided the outer loop and, separately,
those of the final predictor.

**A block of sample cells (review of M5, round 2).**  The predictor sees
the neighbour cells as they were after the previous corrector.  In the
cells of an inventory sample with a small reservoir (e.g. run/013 with 1e-7
mol, `Cs` = 5.7 kg/m3, or the test channel with 1.7e-8 mol), every cell
switched at once at every corrector: all CAPPED (the neighbours are
empty, so each cell would lose more than its reservoir), which dumps the
reservoir into supersaturated gas, then all NONE (the neighbours are now
supersaturated), which drains the gas, and so on.  The parity of the last
corrector decided the step: run/013 at `nCorr 10` never evaporated the
boat, at `nCorr 9` or 11 it dumped it in the first step, and neither state
obeyed the law.  The complementarity check did not see it, since it
counted only elements whose law lay within the bounds.  Now

- the predictor never jumps between the two bounds of a sample element
  (CAPPED <-> NONE at `lower = 0` passes IMPLICIT, whose coupled solve
  resolves the whole block; `acrossBounds` in `phaseChangeKinetics.H`; a
  wall element, whose NONE is the upper bound 0 of a bare wall, is never
  damped), so the loop settles within
  a few correctors whatever `nCorr` is: run/013 at 1e-7 mol converges in
  every step with 6 to 8 correctors, identically for `nCorr` 9, 10 and 11,
  and the channel's regimes are the same in every step for `nCorr` 2, 3,
  10 and 11 (M5.2k, M5.3d);
- the guard switches in both directions (above), so even one corrector
  (`nCorr 0`) ends in a state that obeys the law (M5.2l), and missed
  deposition on a bare wall is repaired too (M5.2b: significant in 56 of
  350 steps with the round-2 guard, now in none);
- `complementarity` counts every element that is not IMPLICIT and whose
  realised flow differs from its admissible flow, e.g. a CAPPED sample
  cell in supersaturated gas: with `nCorr 0` and `nGuard 0` the channel
  reports the 34 such cells of its first step, as a check recomputed from
  the written fields finds them (M5.2l).

Only the adjacent switch back to IMPLICIT is safe at the solved Y: with
`lambda dt` of 1e3 to 1e4 the gas of a cell CAPPED and of the same cell
NONE differ by far more than the width of its bounds, so switching to the
bound that is admissible at the solved Y would alternate between the two
and never reach the implicit solution between them.

### Start-up checks

The run stops with a clear message for: a ddt scheme other than `Euler`
for a paired species (also when `fvSchemes` is changed so during the run);
a convection scheme of `div(phi,Y_<gas>)` other than `Gauss
<interpolation>`, e.g. `bounded Gauss upwind`, which subtracts
`fvm::Sp(div(phi), Y)`, a cell source the ledger does not book (the run/012
pipe then closed to 8.7e-12 instead of 3e-16 in 5 steps; also when
`fvSchemes` is changed so during the run; M2.6k); a boundary condition
other than `zeroGradient` of the paired Y on an interface patch; a final
solve with `relTol > 0` (also when `fvSolution` is changed so during the
run); `wallResistance halfCell` with `areaMultiplier` other than 1 or
`areaPerVolume layer` (M2.6l); a sample that overlaps the interaction layer
or another sample, or selects no cell; a pair without `molarMass`; a fresh
start after the release window of a sample has begun (that part of n0
would never be released; the message gives the missed fraction);
`equilibrium GEMS` in a build without the bridge (the rebuild hint; M9.1),
and with it a pair without its `gems` block, a formula, species, molar
mass or carrier that does not fit the GEMS3K system, a condensed species
as the pair's gas or a condensate of another formula, `speciation lagged`
without `equilibrium GEMS` in mode local, a missing `GEMSCoeffs`, an
unknown mode, an unknown key of `GEMSCoeffs` or of a `gems` block,
`referencePressure` outside the pressure grid (mode frozen), a system file
without an engine (M9.7), and mode local on a carrier that is not a gas at
its p and T unless `allowNonPhysicalCarrier yes` (M9.3e); `layer distance` with
`wallResistance halfCell` or `areaPerVolume cell`, without `distance`, or
with a distance that selects no cell (M4.8); a release sample with a key
of mode inventory (`removeTime`,
`areaPerVolume`, `kineticScale`, `Ce`: "belongs to mode inventory, not to
mode release"; before M2's final pass they were ignored silently; M2.6m),
an inventory sample with a key of mode release (`startTime`, `duration`),
and `mode inventory` with `model temperature` (M5.5, M2.6f); for `model
HKS`: a pair without `vapourPressure`, an unknown `units`,
`interpolation` or `outOfRange`, a phase that is not a column of the
table, a malformed table (non-increasing T, a missing or non-numeric
value, p <= 0 with `logInverseT`), `outOfRange fatal` with element
temperatures outside the table, and `HKSCoeffs` outside `0 < accommodation
<= 1`, `Ce >= 0`, `kineticScale >= 0` (M5.5); an inventory sample whose
removal would fall into the first step of a fresh start (M5.5).  It warns
about the inputs of the other model, which are not read (`HKSCoeffs` and a
pair's `vapourPressure` under `model temperature`, `temperatureCoeffs`
under `model HKS`; M5.9; with equilibrium table the inputs of
equilibrium GEMS, and in `GEMSCoeffs` the keys of the other mode; M9.11),
about a relaxation factor of a paired species
(ignored) and an outer tolerance below the solver tolerance, and prints the
carrier consistency `max |T - p W/(R rho)|/T` (0.30 % for the 007 carrier).
`constant/thermochemistryProperties` is read at start-up only.

### Mass ledger and balance line

For every pair the ledger books, in kg: `initial` (gas, wall and the
reservoirs filled on a fresh start), `released`, `transport` through every
non-coupled patch (outward; the INLET back-diffusion is a separate entry),
`removed` (the reservoirs of inventory samples at their removal), `clamped`
(the round-off residue of the condensate: of a CAPPED element set to 0,
and the rounding of every other condensate update, which changes the
stored condensate by the rounded update rather than by exactly `dt q`),
`solverDefect` and `restart`, and keeps the inventory held at the last
step.  The closure is

    closure = gas + wall + sample
            - (initial + released - transport - removed - clamped
               - solverDefect + restart)

`solverDefect` is the sum over the steps of `R*dt`, R the defect of the
solved equation formed per cell: the defect of the transport-only matrix in
flux form, plus the realised element flows (the values that update the
condensate) plus the release (`matrixDefect.H`; in exact arithmetic the
defect b - A Y of the assembled matrix, which rounds the stiff exchange
differently, see "Round-off of the stiff exchange" below;
`fvMatrix::residual()` is wrong in parallel and not used).  The flux form
computes every face flux once and enters it into the two cells of the face
with opposite signs (a processor face with the mean of the fluxes of its
two sides, which are computed separately), and measures the storage as the
gas inventory does, `((rho Y) V - (rho Y^n) V)/dt`; every cell is summed in
long double.  The face flux includes the explicit face-flux correction of
the matrix, the non-orthogonal correction of the laplacian (a paired
species is `fluxRequired`, so OpenFOAM keeps its face values), and of the
matrix source only the explicit remainder `E = b - fl(s + c)` is kept, `s`
the storage source of the ddt term and `c = -V div(C)` the source form of
the correction: 0 in every cell for the solver's schemes (Gauss with
upwind or a limited scheme such as SuperBee).  The face fluxes then
telescope and the storage telescopes with the inventory, so the
double-precision imbalances of the discrete operator (the column sums of
the assembled `ddt + div - laplacian`, whose diagonal is the rounded
negSumDiag of the convection and laplacian coefficients plus the rounded
sums of the three matrices; the storage products of the Euler scheme,
`((rho/dt) V) Y` and `((rho/dt) Y^n) V`; the source form `V div(C)` of the
correction, which telescopes only to round-off) are booked in
`solverDefect` with the stiff exchange's, not in the closure.  Until the M4
review, round 3, the source `b - s` stayed in the defect, so the closure
changed by `dt sum_c (b - s)_c`, up to 1.7e-24 kg per step on the
non-orthogonal pipe of run/014, booked nowhere (plan section 29, item 2;
M4.18).  An interpolation scheme with an explicit correction (linearUpwind,
LUST) adds it to the source only, without face values in the matrix: it
stays in `E`, and its rounding reaches the closure (a few eps of the
correction per step).  `balance { reportMatrixResidual yes; }` prints `E`
in the transport line (`explicit remainder <E> kg/s in <n> cells (of an
explicit source sum |b - s| <S> kg/s)`), and the transport identity of the
cell defects themselves (`transport identity <D> kg/s (|D| <r> of its
terms, <T> kg/s)`): D = sum_c r_c + dM/dt + sum_p F_p - sum_c E_c over all
ranks, with the cell defects of the ledger's own code accumulated and
summed exactly (`matrixDefect.H`).  In the flux form every face flux
cancels between its two cells, so D is 0 to the rounding of the exact
sums (at most 1.1e-36 of its terms on the drift pipe of M4.15), while a
defect with the correction in source form leaves the double-precision
residue of `V div(C)` (5e-21 to 8e-20 of its terms there): M4.18 gates
the identity, with a build of that form as its negative control (the
report of `E` alone does not see how the defect is formed; plan section
33).
The closure is therefore at round-off for any linear-solver tolerance and
any run length; `|solverDefect|` measures the error of an incompletely
converged solve (plus those imbalances, a few eps of the inventory per
step) and needs the tight final solve.  The round-off of the closure is
that of the inventory sums and of one step's booking: in the closed box
2e-16, on the pipe 2e-15 of the inventory after 25 steps (a few 1e-14
before the inventories became compensated sums, see below).

The cumulative entries (released, transport, solverDefect, clamped,
removed, restart and the diagnostics) grow by one increment per step.  Each
keeps its Neumaier compensation in the ledger (`compensatedSum.H`,
`phaseChangeLedger.H`), so their totals stay at a few eps however many
steps a run takes.  With plain running totals the release of the 6 s
window of run/012 (6100 steps of 1 ms) missed n0 M by 6.2e-14, and a 15 s
window at 0.22 ms (136364 steps, M8 scale) by 5.2e-13 with a closure of
1.3e-12; compensated, the test channel gives 2e-16 and 2.3e-16 (M2.2g), and
4e-16 and 1.1e-15 (M2.2h, the largest over its 136364 steps).  Until the M4 review, round 2, M2.2h closed only
to 1.6e-13, and the T5_60 runs of M4long (6000 steps of 10 ms) to 7.6e-13:
not a random walk but a linear drift wherever gas or a deposit stays in
place (in M2.2h while the gas stood bitwise steady during the release),
from the imbalances above and the rounding of the stored deposits (one
increment per step, rounded alike in every step), which were booked
nowhere (plan section 27, item 1).  M2.2h (1e-14), M4.15 (a pipe at large
steps, 2e-15) and M4.17 (the long runs, 1e-14) gate it now.

Every step prints, per paired species,

    Y_PbI2_g exchange: implicit ..., capped ..., none ... elements; outer
      iterations ..., converged (initial residual ..., tolerance ...; regime
      changes N, M kg, S significant); final predictor: regime changes N,
      M kg, S significant; guard passes ...; projected N (M kg, S
      significant); complementarity N (M kg, S significant)
    Y_PbI2_g transport [kg/s]: INLET ... OUTLET ... WALL ...; solver defect ...
    Y_PbI2_g monitors: min Y ..., negative mass ... kg, max mole fraction ...
    Balance PbI2_g [mol]: gas ... wall ... sample ... released ... transport
      ... removed ... clamped ... restart ... solverDefect ... (relative ...)
      closure ... (relative ...)

(one line each).  Relative values are relative to max(|initial +
released|, |held|, cumulative inflow), so a pair fed only through an inlet
has a finite reference too.  `complementarity` counts elements that are
not IMPLICIT and whose realised flow differs from the admissible flow
`min(max(a - g Y, lower), upper)` at the final Y (e.g. a missed deposition,
a CAPPED sample cell in supersaturated gas), with the missed mass; a
projected or missed mass above the guard's tolerance in total is printed
as `M kg > tolerance T kg`.  The
monitors give the passive-carrier validity, the mole fraction `x = c/(c +
p/(R T))` with `c = rho Y/M` the molar concentration of the species and
`p/(R T)` that of the carrier at its frozen p and T (its maximum, 0 where
Y <= 0 everywhere), and the mass that entered the gas where `x >
moleFractionWarning`.  For a perfect-gas carrier this is `(Y/M)/(Y/M +
1/M_carrier)`, the form used until the M4 review, round 3; in paper mode
(rho = 1, not the density of helium) that form measured the species
against an imaginary carrier of 250 mol/m3 and was 6 to 25 times too
small (plan section 29, item 4; M4.19).  Measured with the carrier's own
concentration, the gas of run/014 is not dilute: at the end of the
release of Figs. 8-10 (15 s) x exceeds 0.05 in 148530 of the 228960 cells
and 0.1 in 54783, at most 0.12 (the plateau of Fig. 9, 1.25 mol/m3,
against 10 to 16 mol/m3 of helium at 1 bar), so the passive-carrier
assumption of the paper mode, and of T-Flows, holds to about 10 % there.
The end of the run sums the steps without outer convergence, the guard
passes, and the
projections, complementarity violations and final-predictor regime changes
with their significant ones (and the steps above the tolerance in total,
if any).

The master writes the ledger every step to
`postProcessing/phaseChangeBalance/<start time>/<gas>.dat` (kg, 17 digits).
Its time column is the exact time at the end of the step (the clock of the
ledger, below), not the rounded time name.

The gas and wall inventories are compensated (Neumaier) sums
(`compensatedSum.H`), per rank and over the ranks.  They do not depend on
the order of the cells, which matters for restarts onto another
decomposition (below).

### State, output and restart

The restart state of the exchange is kept in `<time>/uniform/` (the
global objects) and `<time>/uniform/phaseChange/` (the per-rank lists):

- `phaseChangeBalance` (global): the ledger (identical on every rank);
- `phaseChangeLayout` (global): a dictionary naming the ledger's pairs, their
  condensates and its transport patches.  A restart with other pairs, or
  with the patches renamed or reordered, stops with a FatalError that
  lists both layouts (the entries are addressed by position).  It also
  records the interface (`layer`, `patches`, `areaPerVolume`,
  `areaMultiplier`) and the release samples (`pair`, `amount`, `startTime`,
  `duration`), with the scalars as exact strings (see "Restart record"
  below);
- `phaseChange/phaseChangeDeposits` (per rank): the deposit m''_e [kg/m2]
  of every interface element of every pair, **the restart state of the
  condensate**;
- `phaseChange/phaseChangeRegimes` (per rank): the regime of every element
  (review item B5), so that a binary restart repeats the continuous run bit
  for bit;
- `phaseChange/phaseChangeReservoirs` (per rank, only with inventory
  samples): the
  reservoir `Cs` [kg/m3] of every sample element of every pair, **the
  restart state of the samples** (M5), behind a fingerprint of the sample
  elements (their cells, samples and cell centres).  The regimes cover the
  wall and the sample elements (their fingerprint is that of the wall
  elements when there is no inventory sample, so the files of earlier runs
  still fit); the deposits the wall elements only.

The ledger, the deposits and the regimes are always written in BINARY, as
raw blocks that read back bit for bit whatever `writeFormat` and
`writePrecision` are (`binaryIOList.H`).  All four are registered with
OpenFOAM and written **whenever it writes the fields**: at the write times
of `controlDict`, and also when a function object writes between two steps
(`Time::writeNow()`, `Time::writeAndEnd()`, e.g. `runTimeControl` with
`satisfiedAction end` or its `nWriteStep`; M2.5p).  The axial profiles
are hooked into the same writes (M3.4d).  Nothing is written unless
OpenFOAM writes (review item B4).  The per-rank lists start with a
fingerprint of the rank's elements (their number and a hash of their cells,
faces and face centres); they are used only for the same elements.

The ledger and the layout are **global** objects, the same on every rank:
the master reads them and broadcasts them.  This matters for
`-fileHandler collated`.  There a per-rank object of a parallel run is one
`decomposedBlockData` container per time,
`processors<N>/<time>/uniform/<name>`, with one block per rank, while a
global object is a plain file written by the master (like OpenFOAM's own
`uniform/time`).  `reconstructPar` and `decomposePar` copy `uniform/`
verbatim, so a container reaches a serial case or a case of another number
of ranks, whose file handler cannot read it (review of M5, round 3: with
the ledger as a per-rank list every such restart stopped with 'could not
detect processor number' or a missing block).  The per-rank lists are
therefore read only when the file handler can read them for this run (a
plain list, or a container of as many blocks as ranks); a container of
another number of ranks counts as written for other elements, with a
warning that names it, and the deposits and reservoirs come from `mw_` and
`Cs_` (M5.2n: the collated channel on 2 ranks, reconstructed and continued
in serial, and re-decomposed onto 4 ranks; both exact).  With several IO
ranks (`-ioRanks`, `FOAM_IORANKS`) there is one container per group of
ranks, `processors<N>_<lo>-<hi>/`, and each rank reads the block of its
own group.  Uncollated, every rank writes all five objects into its own
`processor<i>/<time>/`.

The per-rank lists lie in the subdirectory `uniform/phaseChange/`, never in
`uniform/` itself: the file handlers that read on the master
(`masterUncollated`, and `collated`, which derives from it) take every
object whose local directory is exactly `uniform` for global data and
resolve it to the master's file on all ranks
(`masterUncollatedFileOperation::filePath`, v2412 and v2606).  In
`uniform/` every rank read rank 0's list (or the container of rank 0's IO
group), the fingerprints refused it, and every restart on the same ranks
and decomposition lost the exact state: the deposits and reservoirs came
from the rounded `mw_` and `Cs_`, and 14 (masterUncollated) and 9 (collated,
`-ioRanks '(0 2)'`) files at the next write time differed from the
continuous run (review of M5, round 4).  Both restarts are now exact and
cmp-identical (M5.2q, with a rank that holds a single WALL face).
`decomposePar` and `reconstructPar` copy (or link) the subdirectory with
`uniform/`.

Output fields:

- `mw_<condensate>` [kg/m2]: the element deposits on the interface patches
  (the internal field is the area mean of a cell's elements, for display).
  It is the deposit per **effective element area** `A_e = |S_f| *
  areaMultiplier` (rescaled for `areaPerVolume layer`), not per wall area:
  with `areaPerVolume cell`, `surfaceFieldValue` `areaIntegrate` of `mw_`
  on WALL is the wall inventory divided by `areaMultiplier` (M2.1m), and
  the wall inventory itself only with `areaMultiplier 1` (M2.2j).  The
  deposit per wall area, `sum_e m''_e A_e/|S_f|`, is `mDep_` (M3).  `mw_`
  is refreshed every step, so function objects may sample it between write
  times; its processor patches hold the neighbour cells' values (M2.4g).
  It is read back only when `uniform/phaseChange/phaseChangeDeposits` does
  not fit the elements (below).
- `mDep_<condensate>` [kg/m2]: the deposit **per wall area**.  On the
  interface patches face f holds `sum_{e -> f} m''_e A_e/|S_f|` (one
  element per face for `layer faceCells`), evaluated as `m''_e
  (A_e/|S_f|)`, so that `areaMultiplier 2` gives exactly `2 mw_` (M3.5a);
  `surfaceFieldValue` `areaIntegrate` of `mDep_` on the interface patches
  is the wall inventory for any `areaMultiplier` and `areaPerVolume`
  (M3.5a/b).  The internal field is the deposit per wall area of a cell,
  `(sum_e m''_e A_e)/(sum_e |S_f|)`.  Refreshed every step, processor
  patches evaluated (M3.5c), written with the fields, never read.
- `exch_<gas>` [kg/(m3 s)], only with `balance { writeExchangeField yes;
  }`: the realised exchange of the last step per cell volume, `sum_{e in
  c} q_e/V_c`, positive = evaporation into the gas; the release of the
  samples is not included (M3.6).
- `Cs_<condensate>` [kg/m3] (only with inventory samples): the reservoirs in
  the sample cells, 0 elsewhere; refreshed every step, processor patches
  evaluated, written with the fields; read back only when
  `uniform/phaseChange/phaseChangeReservoirs` does not fit the sample
  elements
  (another decomposition, a renumbered mesh, inventory samples reordered
  or dropped), with a warning, as `mw_` for the deposits (M5.2f, M5.2g).
  The fingerprint of a sample element holds the sample's name, so release
  samples added in front of an inventory sample keep the list (M5.2i).
  When every inventory sample has been removed, neither is needed (the
  reservoirs are 0; M5.2h).
- `Y_<condensate>`: derived output, `(sum_e m''_e A_e + Cs V)/(rho V)`,
  recomputed every step before `runTime.write()`, so that
  `c_<condensate>` is consistent.  Its initial values are ignored (with a
  warning if non-zero).
- the axial profiles, `postProcessing/axialProfiles/` (next section).

Why the deposits are not kept in `mw_` itself: OpenFOAM writes a field or
patch whose values are all equal as `uniform <value>`, **in text at
writePrecision, also in a binary file**.  A rank with a single WALL face
(common with scotch) thus lost the deposit's digits, and the binary restart
was refused with "do not belong to the same run" (review, round 3).  The
same happens to a single-face processor patch of `Y_<gas>`; those values
enter the first matrix of a restart (the SuperBee limiter), so the solver
evaluates the processor patches of every paired gas from the neighbour
cells at start-up.  Without that, the restart of the 3-rank test channel
(M2.5o) was not cmp-identical, even with exact deposits.

- **Release.**  A step releases `n0 M overlap/length`, spread uniformly
  over the sample cells, with `length = (start + duration) - start`, the
  window as the overlaps represent it (its end rounded).  The overlaps of
  all steps sum to this length, so the release is n0 M to round-off also
  far from t = 0; with the nominal duration a window of 0.0198 s at
  t = 1000 s released 1.6e-12 too much (M2.2i).
- **Clock.**  The release of a step is the overlap of the exact step
  `[t, t + deltaT]` with the release window.  The ledger keeps this time
  (`TIME`, set to the start time on a fresh start and advanced by `deltaT`
  in every step), because the time-directory names are rounded to
  `timePrecision` between write times (10.0002, 10.0004, 10.0007 for
  `deltaT = 0.22 ms`) and a restarted OpenFOAM continues from the rounded
  name.  A binary restart continues the clock and repeats the continuous
  run bit for bit (M2.5a/b), also from a rounded name (M2.5g; a warning
  gives the offset).  OpenFOAM itself finds a start directory only if
  `timePrecision` is large enough to name it (10.0108 needs 6); with less
  it silently starts from another directory.

A **fresh start** has neither the ledger nor any condensate state
(`uniform/phaseChange/phaseChangeDeposits`, `mw_`,
`uniform/phaseChange/phaseChangeReservoirs`, `Cs_`) at the start time.
State without the ledger stops the run (e.g. a `Cs_PbI2_s` in `0/`, M5.5),
and so does the ledger with neither the deposits of these elements nor all
`mw_` fields, or with inventory samples but neither the reservoirs of these
sample elements nor all `Cs_` fields.  On a **restart** the deposits come
from `uniform/phaseChange/phaseChangeDeposits`; after `decomposePar`,
`reconstructPar` or `renumberMesh` (other elements) they are re-read from
the interface patches of `mw_`, with a warning.  An `mw_` that differs from the deposits
of `uniform/` is reported and ignored (it is output; M2.5n).  The inventory
re-read minus the inventory held by the ledger is booked as `restart` and
printed, e.g.

    Restart of pair PbI2_g: the inventory re-read differs from the ledger by
      0 kg (relative 0; gas: Y_PbI2_g written in binary; deposits:
      uniform/phaseChange/phaseChangeDeposits; exact), booked as restart

When the inventory re-read is stored exactly (the deposits of `uniform/`,
or an `mw_` in binary without a repeated non-zero interface-patch value;
`Y_<gas>` in binary, and no rank whose `Y_<gas>` internal field is a rounded
`uniform` entry) the term must be below 1e-14 (FatalError otherwise: the
files do not belong together).  With `writeFrozenFields yes` the carrier is
re-read from the restart directory, so `rho` is part of the gas inventory
re-read and of this rule; the restart line then names it, e.g. for a
uniform carrier computed with 17 digits and written at 6

    Restart of pair PbI2_g: ... (relative -3.24463e-08; gas: Y_PbI2_g written
      in binary; carrier: rho written in binary, a 'uniform' value on some
      rank (text at writePrecision); deposits:
      uniform/phaseChange/phaseChangeDeposits; rounded for the next steps: rho ('uniform' on INLET OUTLET WALL TOP);
      rounded), booked as restart

(before, such a restart was refused as "not the same run"; M2.5s).  Otherwise
(`rounded`, e.g. `writeFormat ascii`) the rounding is booked (about 1e-7
relative at the default precision 6), and a warning is printed above
`10^(1 - writePrecision)`.

The verdict also covers the values the next steps re-read beyond the
inventory (review, final M2 pass; before, the verdict ignored them and
said `exact` for restarts that did not repeat the run): `Y_<gas>` on the
non-coupled patches (a `fixedValue` or `inletOutlet` value; a
`zeroGradient` patch is evaluated from the cells), and with
`writeFrozenFields yes` the patches of `rho` and `T`, `p` and `phi`
(format, internal field, patches, the processor patches of `phi`).  A
rounded one leaves the inventory exact, so the bound 1e-14 still applies,
but the continuation is not the continuous run:

    Restart of pair PbI2_g: the inventory re-read differs from the ledger by
      0 kg (relative 0; gas: Y_PbI2_g written in binary; deposits:
      uniform/phaseChange/phaseChangeDeposits; rounded for the next steps:
      Y_PbI2_g
      ('uniform' on INLET); inventory exact, continuation rounded), booked
      as restart

(an inlet value 1.234567891e-3 written as `uniform 0.00123457`; the files
at the next write time differ from the continuous run, M2.5x; the same for
a carrier written at 12 digits and re-read at 6, M2.5y).

What the verdict covers (M3, plan section 18): the inventory, the
internal fields and the patch **values** that a restart reads back.  A
patch whose value is recomputed when it is read does not count:
`zeroGradient`, `fixedGradient` (v2412 does not even write its value),
`mixed` (evaluated from `refValue`, `refGradient` and `valueFraction`),
`symmetry`, `symmetryPlane`, `wedge`, `empty` and every coupled patch
(`cyclic`, processor).  `inletOutlet` reads its value back and is checked.
Boundary **parameters** that OpenFOAM writes as `uniform` text as well
(`inletValue` of `inletOutlet`; `refValue`, `refGradient`,
`valueFraction` of `mixed`; `gradient` of `fixedGradient`) are not
checked; the restart line ends with

    ..., booked as restart; checked: the inventory, the internal fields and
      the patch values re-read, not other boundary parameters (e.g.
      inletValue)

Before, a single-face `fixedGradient` or `mixed` patch (a rounded
`uniform` value in memory or in the file) gave `inventory exact,
continuation rounded` although the continuation was identical (M3.8).  Only the run
that writes a `uniform` entry can tell whether it lost digits (the value
re-read always reproduces its own text), so the solver records, per rank
and before every write, the fields it checked and their rounded entries in
`uniform/phaseChangeLayout` (`uniformValues`).  A value that itself came
from text of at most `writePrecision` digits, such as an inlet value
`1e-3`, `p 101325`, or a carrier written at that precision, is written
exactly and the restart stays `exact` (M2.5o, M2.5s).  Without a record
that belongs to the files (an earlier build, or after `decomposePar`,
`reconstructPar` or `renumberMesh`, which rewrite the fields) every
repeated non-zero value counts: `possibly rounded for the next steps`
(M2.5x).  Use `writeFormat binary`.  The format is read
from the header of the field files; with `-fileHandler collated` every file
is a `decomposedBlockData` container whose header always says `format
binary`, so the format of its data is taken from `data.format` (or from the
header of its first block).  Binary and ascii restarts work with the
collated handler (M2.5l/m).

In parallel the solver evaluates the processor patches of the frozen
carrier (`T`, `p`, `rho`, `U`) once from the neighbour cells after reading
it.  `decomposePar` writes a processor patch whose values are all equal as
`uniform <value>` with about 10 significant digits, also in a binary file,
while the neighbour's cell values are exact; the two sides of a processor
face then saw a different `rhoD`, and the diffusive flux leaving one rank
was not the flux entering the other (the 4-rank test channel closed only
to 2.3e-13, now 2e-15; M2.4d).  Serial runs have no processor patches and
are unchanged.  The evaluation is done for **every model, including
`mock`**: a decomposed mock run whose carrier has such patches (e.g. a plug
flow) differs from a build without it in the last digits (plan section
15.7).  `phi` is a face field and cannot be evaluated; a restart with
`writeFrozenFields yes` re-reads it from the restart directory, where a
single-face processor patch is rounded, so bit-identical parallel restarts
need `writeFrozenFields no` (the carrier is then read from `0` in both
runs).

**Restart record.**  The deposits of a restart are per effective element
area on the faces of the interface patches, and the ledger has booked the
release up to the restart under the schedule of that run.  A restart
therefore stops with a FatalIOError that names the change

- when the interface differs from the recorded one (`patches`,
  `areaPerVolume`, `areaMultiplier`): the wall inventory would change by
  the ratio of the areas; in ascii that change used to be booked silently
  as rounding (M2.5t);
- when a sample whose window has begun is changed (`pair`, `amount`,
  `startTime`, `duration`), when a sample is removed before its window has
  ended, and when a new sample's window began before the restart (M2.5u);
- for inventory samples (M5.5): when one is new (its reservoir is filled on
  a fresh start only), changed (`pair`, `amount` or its number of cells,
  compared exactly; its volume or its volume-weighted centre, compensated
  sums over the ranks compared to 1e-12, since their last bits depend on
  the decomposition, plus the error bounds that the rounding of the mesh
  points implies in the run that wrote the state and in this one: the
  reservoir belongs to those cells, and the centre names a selection
  moved to other cells of the same number and volume), dropped before its
  removal, or when any sample changes its mode.
  `removeTime` and the kinetics (`areaPerVolume`, `kineticScale`, `Ce`) may
  change.  The removals recorded per sample must add up to the ledger's
  count (`REMOVAL_DONE`), otherwise the files of `uniform/` belong to
  different runs.  A sample dropped after its removal keeps its record for
  the rest of the run, so later restarts still find the count (M5.2h);
  listed again under its name, it stays empty.  A restart after
  `decomposePar`, `reconstructPar` or `renumberMesh` passes, in binary
  (M5.2f, M5.2g) and in ascii (M5.2m, M5.2o).  With `writeFormat ascii`
  these utilities store the mesh points as text (`decomposePar` with at
  least 10 significant digits), which moves V_Z and the centre of the same
  cells: by 7e-10 (relative) for the boat of run/013 at x = 0.04 m, by
  5.8e-9 at 10 digits for the same boat moved to x = 0.6 m, since the
  error grows with |x| over the cell size.  Every run therefore records
  with the sample the error bounds of its own points (`volumeTolerance`,
  `centreTolerance`, `meshPoints` in the layout): 0 for binary points; for
  ascii points, with p significant digits, every point is off by at most
  sqrt(3) 0.5 10^(e + 1 - p), e the exponent of the largest coordinate,
  and V_Z by at most that times the surface of the sample (the face areas
  of its cells; the centre times the largest distance of a point from it,
  over V_Z; both doubled).  p is the most significant digits any
  coordinate of the sample's cells needs, but at least min(6,
  writePrecision), the fewest digits OpenFOAM's utilities write points
  with (`foamFormatConvert`, `renumberMesh` at writePrecision, default 6;
  `blockMesh`, `gmshToFoam`, `decomposePar` at least 10): the digits a
  value needs are no bound on their own, since a round coordinate
  (`0.011`) needs 2 (review of M5, round 4: on an ascii mesh of round
  coordinates the bounds exceeded the channel, and a restart with the
  sample moved by 20 mm was accepted, its reservoir re-read as 0 from
  `Cs_` in the new cells).  A restart allows the bounds of both runs, so
  the state of an ascii-decomposed run reconstructed onto the binary mesh
  passes, and so do writePrecision 10 and 12 for the boat at x = 0.6 m
  (M5.2o; the round-3 build used the reading run's points only, with
  10^(1 - writePrecision), and refused all three).  It never allows more
  than what a selection changed by one cell causes: for the centre a
  quarter of h_min V_min/V_Z (h_min the thinnest cell of the sample, its
  volume over its largest face, V_min the smallest cell volume; one cell
  replaced by another of the same volume moves the centre by about a cell
  distance times V_min/V_Z), for the volume half of V_min.  A moved
  selection, and one with a single cell replaced by its neighbour (5e-6 m,
  within the rounding bound of 6-digit points), are thus refused on an
  ascii mesh of round coordinates too (M5.5), while the unchanged restart
  passes; when the difference lies within the rounding bound but above the
  limit, the refusal says that the points are too coarse to tell and
  suggests binary points or more digits.  For binary points the
  comparison stays at 1e-12.  The
  selection itself must not pass through cell centres: a cylinder radius
  on a ring of centres selects other cells once the points are rounded
  (the boat of run/013 had `radius 5e-4`, with a ring of 18 centres
  2e-12 outside it, and selected 146 cells on an ascii-decomposed mesh; it
  has `4.8e-4` now, the same 144 cells; run/013/readme.md).  A restart then
  stops with the changed number of cells.

These checks run before the restart booking, so a dropped or changed
sample is named as such, not reported as a failed round-off bound.
Release samples whose windows lie after the restart may be added, changed
or removed, and a release sample may be removed once its window has ended
(M2.5v).  A layout written by an earlier build, without these records,
gives a warning (M2.5w).

In parallel the restart directories must hold `uniform/` of every
processor (with `-fileHandler collated`, of `processors<N>/`).  A restart
onto another decomposition (`decomposePar`, `reconstructPar`) or a
renumbered mesh passes the binary bound, because the inventories are
compensated sums (plain sums differed by 1.5e-14 on the pipe and stopped
such a restart), re-reads the deposits from `mw_` and starts with unknown
regimes (warnings), even when a rank happens to have as many elements as
before; with the collated handler as well (M5.2n, "State, output and
restart" above).

### Axial profiles and T_dep

With a `profiles` dictionary (M3; `axialProfiles.H`, `profileBinning.H`)
the solver writes, at every write of OpenFOAM and at every entry of
`times`, the amount of every pair per unit length along the axis, `x =
(point - origin) & axis`, in mol/m, for every bin set:

    n'_g(b) = sum_c w_cb rho_c Y_c V_c/(M dx_b)    gas
    n'_w(b) = sum_e w_eb m''_e A_e/(M dx_b)        wall deposit
    n'_s(b) = sum_c w_cb Cs_c V_c/(M dx_b)         sample reservoir (inventory
                                                   samples; 0 for release)
    sampleGas(b)                                   the gas of the cells of every
                                                   sample (part of n'_g)

The weights are the **exact overlap** of the axial extent of a cell (from
its points), or of the wall face of an interface element, with the bin,
divided by the extent; the part outside all bins is reported separately.
So `sum_b n'(b) dx_b M` plus the amount outside is the gas, wall and
sample inventory of the ledger to round-off (1e-16 in the tests, gate
1e-12: M3.1a-c), for the native bins (424 on the pipe, the axial cells),
the 1 cm gamma-scan bins, and bins that do not cover the mesh or whose
width does not divide the range (a shorter last bin; a width at or above
the range gives one bin, and more than 1e6 bins stop the run).  Bins that
cut cells hold the exact overlap of the cell-aligned profile (M3.1d).  The
wall temperature of a bin is the area-weighted wall-face temperature of
the interface patches (frozen).

**T_dep.**  The observable is the wall deposit (`wall`) or the deposit
plus the gas left in the tube outside the cells of every sample
(`wallPlusGas`, the default: that gas condenses locally at shutdown;
`wall + gas - sampleGas` per bin; decided with the sign-off of M3, plan
section 19, item 6, and implemented in M5); the sample reservoir is always
excluded.  The onset is where the observable first reaches `fraction`
(0.01) of its peak:

- `fromInlet`: the first bin, from `searchFrom` on (the bins whose centre
  lies at or after it; default all), that reaches the fraction of the
  peak of that range;
- `fromPeak`: from the peak, upstream to the last bin that still reaches
  it: the upstream edge of the main deposit (e.g. T5_60 with an upstream
  spike, where the two differ).

x_dep and T_dep are interpolated linearly between the centres of the two
bins that bracket the threshold.  An onset at the first bin of the range
(the observable is above the threshold from the start) or without a
positive peak is reported as `none`, never extrapolated.  On a synthetic
linear ramp deposit both modes give the analytic crossing (814 K) to
1e-9 K, and with an upstream spike fromInlet finds the spike (942.8 K)
and fromPeak the ramp (M3.3a-c).

`wallPlusGas` excludes the gas **in** the sample cells, but not the gas
cloud around them.  While that cloud is still in the tube its peak and
onsets lie on it: in the test channel at t = 0.03 `wallPlusGas` gives
949.9 K upstream of the sample box (949.5 K when all gas counted; M5.6),
because the cells of the same axial slices around the sample hold most of
the cloud.  At shutdown, with the sample removed and the gas left, both
definitions give the deposition onset.  At intermediate profile times use
`observable wall`, or read the onset of the wall column that
`utilities/axialProfiles.py` prints next to a `wallPlusGas` file (`wall
column: T_dep ...`).

**Significance.**  A search whose peak bin holds at most `significance`
times the inventory of the pair (gas + wall + sample, in mol; default the
`guardTolerance` of `balance`, 1e-12) reports `none (insignificant: the
peak bin holds ... mol, at most ...)`, status `insignificant` in the file:
a column of solver noise has a peak and a crossing but no deposit.  In
run/012 at 4 ms the wall column holds 1.5e-37 mol (noise ahead of the gas
front): `insignificant`; with `significance 0` it is searched (M5.7).  The
wall column of the test channel at t = 0.03 (7e-22 mol of 3.5e-6) is
insignificant as well; the review of M3 had read 684.9 K from it.

**When.**  At every write of OpenFOAM: the write times of `controlDict`
and the writes that function objects start between two steps
(`Time::writeNow()`, `Time::writeAndEnd()`, e.g. `runTimeControl` with
`satisfiedAction end` or `nWriteStep`).  The profiles are hooked into
OpenFOAM's write like the ledger (an object registered `AUTO_WRITE` whose
write produces them), so every state directory has its profiles, also
the one where a stop-on-condition run ends (M3.4d; review of M3, round
1).  An entry t_i of `times` is written at the end of the step whose end,
on the ledger's exact clock, is nearest to it (the first step end t1 >=
t_i - deltaT/2); entries on one step or on the step of a write give one
set of files (a function object writing after the files of a profile
time rewrites them with the reason `write time; profile time ...`,
byte-identical to a write time of `controlDict` at that step).  At a
profile time that is not a write time **only**
`postProcessing/axialProfiles/<time>/` is written: no field, no
`uniform/`, no time directory (review B4; M3.4a/c).  The entries written so
far are recorded with the state (`uniform/phaseChangeLayout`, `profiles {
written ("0.02" ...); }`, exact).  A restart skips exactly those, and the
entries at or before its start that were not written (e.g. added since),
and writes the others, byte-identical to the continuous run (M3.4b),
whatever the deltaT of either run (M5.8; before, the skip rule used the new
run's deltaT: a restart with a smaller deltaT wrote an entry twice, one
with a larger deltaT lost one).  On a fresh start (or from a layout without
the record) the entries whose nearest step end is not after the start are
skipped, also an entry up to deltaT/2 after the start (the initial state
has no profiles unless OpenFOAM writes it).

**Files**: `postProcessing/axialProfiles/<time>/<gas>_<set>.dat`, e.g.

    # axial profiles of the pair PbI2_g -> PbI2_s, bin set gamma (...)
    # time 0.10000000000000007 directory 0.1 reason write time; profile time 0.1
    # molarMass 0.46100893999999998
    # inventory gas 4.8392776707465876e-06 wall 7.2108506060177571e-10 sample 0
    # binned gas ... / # outside gas ...
    # observable wallPlusGas
    # fraction 0.01
    # significance 1e-12 minimumAmount 4.84e-18
    # Tdep fromInlet status found T 908.03846016952514 x 0.015326923305079192 peak ... xPeak ... threshold ...
    # Tdep fromPeak status found T ...
    # columns x_lo x_hi x gas wall sample observable Twall wallArea sampleGas
    # units m m m mol/m mol/m mol/m mol/m K m2/m mol/m
    0 0.01 0.0050000000000000001 7.2150123097554222e-10 0 0 ...

with 17 significant digits, and two log lines per pair and set,

    Profiles PbI2_g, set gamma (write time; profile time 0.1): in the bins
      gas 4.83928e-06 wall 7.21085e-10 sample 0 mol, outside 0 0 0 mol;
      relative difference of bins + outside from the inventory 0
      -1.55519e-16 0
      T_dep of wallPlusGas at 0.01 of the peak: fromInlet 908.038 K at
      x = 0.0153269 m (peak 0.000230558 mol/m at x = 0.035 m); fromPeak ...

**Reader**: `utilities/axialProfiles.py <case> [--time ...] [--set ...]
[--pair ...] [--plot file.png] [--column observable|gas|wall|sample]`
prints, for every file, T_dep of every mode (as written, and recomputed
from the columns with the same significance), for a `wallPlusGas` file the
check `observable = wall + gas - sampleGas` and the onset of the wall
column alone, and the integrals next to the inventories; `--plot` needs
matplotlib.  As a module: `read()`, `integrals()`, `onset()`,
`files()`.

**Parallel.**  The binned amounts are summed over the ranks with one
`Pstream::listCombineReduce` per pair and bin set, of a list sized by the
dictionary only (review B1), on every rank.  4 ranks equal serial to
6e-14 of the peak (gas) and 4e-14 (wall) on the pipe, `simple (4 1 1)`
and `hierarchical (2 2 1)`, and so do 3 ranks of an uneven manual split
along x (M3.2c, carry-over of the M3 sign-off), **when the outer loop is converged** to the
floor of the linear solver (correctors 1e-13, relTol 0; outer tolerance
1e-12; M3.2), with T_dep to about 1e-11 K and the onset's position, peak
and status equal.  The profile is a stricter measure than the inventory:
with
the forced loose correctors of M2.4a-c (inventories equal to 2e-12) the
gas profile differed by 7e-9 of its peak (plan section 18, item 9).

## Current limitations

- Gaseous species are transported; solid species accumulate without transport.
- Species interact only within a gas/condensate pair: in `model mock` the
  only interaction is `PbI2_g` -> `PbI2_s`; with a phase-change model the
  pairs of `thermochemistryProperties` exchange mass.  Otherwise the species
  equations are uncoupled from each other.
- The precipitation source of `model mock` is a mock model, restricted to
  cells adjacent to `WALL`.
- Solid species have no convection or diffusion equation.
- GEMS3K enters only through the vapour pressure of the HKS law and,
  optionally, the speciation factor of its partial pressure (`equilibrium
  GEMS`, M9); the transported species stay one gas per pair, its formula
  that of the condensate.  The guards keep GEMS3K away from trace
  compositions and from below 500 K, where it fails or errs silently (plan
  section 8, risk 7); in a gas that is not dilute (PbI2 mole fractions of
  0.05 and more at 1170-1240 K) it does not converge and its calls are
  about 100 times dearer, which no guard avoids (the failed calls take the
  table; plan section 35, item M).  Its warm states are not written, so
  the first update of a restart in mode local is cold and the restart
  follows the continuous run from then on to the tolerance of GEMS3K only.  The
  optional xGEMS start-up call of the main branch (`LESTO_XGEMS`, see
  Building) only equilibrates a demo system.  Several pairs are read (one
  GEMS3K call serves all pairs of an element), but the balance of a pair
  requires the same formula for gas and condensate. M10a adds independent
  PbI2/BiI3 pairs and element accounting; re-speciation and shared/reactive
  condensates remain future M10 work.
- `layer distance` has no `wallResistance halfCell` (a cell of the second
  layer has no half cell of a wall face).
- The carrier of `setFrozenCarrier` is the T-Flows parity carrier of the
  code-to-code cases (density 1, a non-expanding parabola,
  `manuscript/notes-for-authors.md`, item A), not a physical one (M8).
- The closure of a stiff HKS exchange is limited by `eps lambda dt` of the
  gas in the exchanging cells per step (see "Inventory samples" above).
- The carrier flow is frozen and is not modified by species transport
  (check the mole-fraction monitor).
- `PbI2He` is currently the only variable-diffusivity model and is restricted
  to `PbI2_g`.
- Other gaseous species currently use constant molecular diffusivity unless a
  new diffusivity model is implemented.
- Species transport currently has no `fvOptions` source-term path.
- Carrier flow must already be converged.

## Verification of the PbI2-He correlation

Numerical reference values are recorded in `tests/testDiffusivity.C`.  This is
an independent small C++ test and does not require OpenFOAM libraries.

From the solver directory:

    g++ -std=c++17 -Wall -Wextra tests/testDiffusivity.C -o /tmp/testDiffusivity
    /tmp/testDiffusivity

The test checks reference diffusivity values from the same explicit Python
correlation and verifies the inverse pressure scaling, including the bar/atm
conversion factor.  The numerical reference curve corresponds to a correlation
pressure value of 1; the current CFD cases use the atm convention described
above.

It also runs as criterion M0.0 of `tests/Alltest M0`.

## Building

From the solver directory, with OpenFOAM loaded:

    ./Allwmake

builds `$FOAM_USER_APPBIN/rhoFixedFlowFoam` and the utility
`$FOAM_USER_APPBIN/setFrozenCarrier` (`utilities/setFrozenCarrier`, M4)
with `wmake`.  `Allwmake` accepts OpenCFD OpenFOAM v2412 and v2606 (`WM_PROJECT_VERSION` 2412, v2412,
2606 or v2606) and stops for any other version.  The build must produce no
compiler warning.

Besides v2412, the solver has been compiled with 0 warnings against the
v2506 development tree and against the official v2606 source release
(OpenFOAM-v2606.tgz, api 2606), built locally with gcc 7.5 (not a Merlin7
installation: Merlin7 has no OpenFOAM module).  Against v2606 it was also
linked, and a short `writeFrozenFields no` run and its restart were run.
`tests/Alltest M0` reports M0.4 as `NOTRUN` unless it runs under a v2606
installation.

`Allwmake` keeps three stamps of the last successful build in
`Make/$WM_OPTIONS`: `gems.state` (the GEMS3K bridge state, below),
`xgems.state` (the xGEMS state, below) and `openfoam.state`
(`$WM_PROJECT_DIR`).  When one of them differs from the current one, or
objects exist without it, `Allwmake` runs `wclean` first, e.g.
`OpenFOAM installation changed (<old> -> <new>): wclean` (of the solver
and of the utility).  v2412 and v2606 share `WM_OPTIONS`, so without this check objects and dependency files of
one installation would be reused by the other.  After switching OpenFOAM
installations, build with `Allwmake`, not a plain `wmake`.

Keep every comment block of a header below 16 KB.  wmkdepend of v2412
reads a comment as one token into a 16 KB buffer; at a longer one it
prints `wmkdepend: buffer full while scanning '<file>'` and stops
scanning that file, so the `#include` lines after the comment are missing
from the dependencies and a changed header no longer recompiles its users.
The description of `interfaceExchange.H` is therefore split into several
comment blocks (each below 16 KB), and the tests count that message as a
warning (BUILD).

### Optional GEMS3K bridge

GEMS3K is reached through the plain-C bridge `thermochemistry/gemsbridge`
(`gemsbridge.h`, `libgemsbridge.so`).  It is an optional build dependency and
has no default path:

    export LESTO_GEMSBRIDGE=<LESTO>/thermochemistry/gemsbridge
    ./Allwmake

- Unset, empty, or a directory without both files: the solver is built
  without GEMS3K.  `Allwmake` warns when the variable is set but incomplete.
- A relative path is resolved by `Allwmake` against the directory it is
  called from (and printed).  A plain `wmake` runs in the solver directory,
  so it refuses a relative path with a clear message instead of resolving
  it against the wrong directory.
- `GEMS_INC` and `GEMS_LIBS` are always assigned by `Make/options`, so
  variables of the same name in the environment never reach the compile or
  link line of a build without the bridge.
- `rhoFixedFlowFoam -help` states the build: `GEMS3K bridge: compiled in
  from <dir>` or `GEMS3K bridge: not compiled in`.  Run logs are unchanged.
- `ldd $(which rhoFixedFlowFoam) | grep gemsbridge` lists the library only in
  a bridge build.  The bridge directory is also the run path (`-rpath`).
- Toggling is safe in both ways of building:
  - `Allwmake` keeps the stamp `Make/$WM_OPTIONS/gems.state` and runs
    `wclean` before building when the state has changed (or is unknown
    although objects exist): `GEMS3K bridge state changed (none -> <dir>):
    wclean`.
  - A plain `wmake` is safe too.  `Make/options` calls `Make/gemsConfig.sh`,
    which writes the generated header `Make/$WM_OPTIONS/gemsConfig.H`
    (`LESTO_HAVE_GEMS` defined or not) only when the state changes.
    `gemsEquilibrium.C` includes it, so wmake recompiles it after a toggle
    and relinks the executable with or without `-lgemsbridge`.
- An edit of the bridge header is picked up as well.  In a bridge build
  `gemsConfig.H` records the checksum of `gemsbridge.h`, so an edited
  header rewrites `gemsConfig.H` and a plain `wmake` or `Allwmake`
  recompiles `gemsEquilibrium.C`.  wmake's own dependency tracking cannot
  see `gemsbridge.h`, which is included with angle brackets (see
  `gemsEquilibrium.C`).
- The bridge serves `equilibrium GEMS` (M9; section "p_eq from GEMS3K"
  above): `gemsEquilibrium.C` calls it only under `LESTO_HAVE_GEMS`.  The
  default build stays without it: `requireBridge()` stops a run that
  selects `equilibrium GEMS` with a FatalIOError and the rebuild hint, and
  the M9 tests build their own copy with the bridge in the work area
  (`FOAM_USER_APPBIN` pointed there), so the solver in `$FOAM_USER_APPBIN`
  and the regression gate are untouched.
- `thermochemistry/gemsbridge/libgemsbridge.so` is git-ignored; it is
  built with `thermochemistry/gemsbridge/build-gems3k-static.sh` and its
  `Makefile` (the conda toolchain of `~/opt/gems-env`;
  `thermochemistry/README.md`).  It depends on libc and libm only, so the
  solver built with gcc 7.5 links it.  The library of this workstation is
  not portable as it stands: it is linked against the glibc 2.34 of the
  conda-forge sysroot (symbol versions up to `GLIBC_2.34`), so it needs
  glibc >= 2.34 where it runs (this workstation has 2.38), and it carries
  an RPATH into the local conda environment (`/home/jan/opt/gems-env/lib`;
  harmless where that directory does not exist, since it needs only libc
  and libm, but tied to this machine).  On another machine, e.g. Merlin7,
  rebuild it there with `build-gems3k-static.sh` and `make` instead of
  copying it (plan section 35, item L).

### Optional xGEMS start-up call (main branch)

The main branch calls xGEMS once at start-up (71bf192):
`rhoFixedFlowFoam.C` constructs an `xGEMS::ChemicalEngine`, initialises it
with a cement demo system of the xGEMS sources
(`CemGEMS-keyvalue/CemHyds-dat.lst`, at a path of that author's machine),
equilibrates it 1 K above its own temperature and prints

    xGEMS ChemicalEngine constructed successfully
    xGEMS status: <status>, converged: <0|1>, iterations: <n>

It tests the link: nothing uses the result.  On the main branch the xGEMS
include and library flags were unconditional
(`$(HOME)/Development/GEMS-Related/install`, 077ca12), so the solver built
only where that installation exists.  Here the call and its flags are an
optional build dependency without a default path, handled like the bridge:

    export LESTO_XGEMS=$HOME/Development/GEMS-Related/install   # an xGEMS prefix
    ./Allwmake

- The prefix must hold `include/xGEMS/ChemicalEngine.hpp` and
  `lib/libxGEMS.so`, and `<eigen3/Eigen/Dense>` must be found (in the
  prefix or in a system directory).  Unset, empty or incomplete: the solver
  is built without the call, and needs neither xGEMS nor Eigen.  `Allwmake`
  warns when the variable is set but incomplete.
- `Make/options` then adds the flags of the main branch with the prefix:
  `-I$LESTO_XGEMS/include` and `-L$LESTO_XGEMS/lib -lxGEMS
  -Wl,-rpath,$LESTO_XGEMS/lib`.  `Make/xgemsConfig.sh` defines
  `LESTO_HAVE_XGEMS` in the generated header `Make/$WM_OPTIONS/xgemsConfig.H`,
  which `rhoFixedFlowFoam.C` includes; the call, unchanged with its demo
  path, is compiled only then.
- Toggling is safe as for the bridge.  `Allwmake` keeps the stamp
  `Make/$WM_OPTIONS/xgems.state` and runs `wclean` when it changes
  (`xGEMS state changed (none -> <prefix>): wclean`); a plain `wmake`
  recompiles `rhoFixedFlowFoam.C` through `xgemsConfig.H`, which also
  records the checksum of `ChemicalEngine.hpp`.  A relative path is
  resolved by `Allwmake` and refused by a plain `wmake`.
- `rhoFixedFlowFoam -help` states the build: `xGEMS start-up call:
  compiled in from <prefix>` or `not compiled in`; `ldd` lists `libxGEMS.so`
  only in an xGEMS build.  Without xGEMS the run logs are unchanged.
- Point the demo path of the call at the demo of the local xGEMS sources.
  On the development workstation, with the conda build of xGEMS 2.1.2
  (`LESTO_XGEMS=~/opt/gems-env`, built with gcc 16) and the path pointed at
  `~/src/xGEMS/demos/resources/` in a scratch copy, the v2412 build (gcc
  7.5) compiled, linked and ran the call (`xGEMS status: 2, converged: 1,
  iterations: 157`; the rest of a log of case 008 unchanged).  The Eigen
  headers of that prefix are not system headers there and give 129
  `-Wold-style-cast` warnings, and GEMS3K writes `ipmlog.txt` into the case
  directory.
- Tests: M0.3j-l, with a dummy xGEMS installation.

### `applications/xgemsFoamProbe` (main branch)

A separate small application of the main branch (64c1c51):
`xgemsFoamProbe <system-dat.lst>` equilibrates a GEMS3K system with xGEMS
and prints the status.  It is not built by `Allwmake` of
`rhoFixedFlowFoam`.  Its `Make/options` takes xGEMS from
`$(HOME)/Development/GEMS-Related/install` unconditionally (edit it for
another installation).  Build it with OpenFOAM loaded:

    cd applications/xgemsFoamProbe
    wmake

which installs `$FOAM_USER_APPBIN/xgemsFoamProbe`.

## Tests

The acceptance tests of each milestone of `doc/phase-change-plan.md` run
with

    set +u; source ~/OpenFOAM/OpenFOAM-v2412/etc/bashrc
    tests/Alltest M0          # or M1, M2, M3, M4, M5, M7, M9, several: tests/Alltest M0 M1, or all
    tests/Alltest M4long      # the code-to-code runs of M4 (hours; not in all)
    tests/Alltest M7long      # the long verification runs of M7 (hours; not in all)

`Alltest` first builds the solver from scratch (`wclean`, then
`./Allwmake`), so every source is compiled and the 0-warnings check is real.
It then runs the milestone's criteria and, for every milestone, the M0
regression gate.  It prints one `PASS`, `FAIL` or `NOTRUN` line per
criterion.  `-fresh` also re-runs the cached pristine reference cases.

The exit status is 0 only when **every expected criterion was recorded**
and none failed.  A test script that stops early (a failed pristine build,
mesh conversion, `decomposePar`, ...) records nothing for the criteria it
did not reach.  `Alltest` therefore reports the script itself as a `FAIL`,
and at the end it records a `FAIL` for every expected criterion id that is
missing.  Criterion `SELF` tests this mechanism on every run: it runs the
gate with a pristine revision that does not exist and requires the abort and
the five missing gate criteria to be reported.

Everything runs in the git-ignored work area `tests/work/`; nothing is
written below `run/`:

- `mesh/`: `gmshToFoam` of `run/001-first-pipe-flow-simpleFoam/pipe.msh`,
  converted once (the stamp is the checksum of `pipe.msh`).  The mesh is
  shared by all OpenFOAM versions.
- `<version>-<WM_OPTIONS>/` (e.g. `v2412-linux64GccDPInt32Opt/`): everything
  else, so that a run under another OpenFOAM (e.g. v2606 for M0.4) never
  reuses builds, compiled coded boundary conditions or cached logs of
  this one:
  - `pristine/`: the solver of the tag `0.1` of the main branch, pinned
    (`PRISTINE_REF`, default `0.1`, a tag, commit or branch; plan section
    31), extracted with `git archive` and built with `wmake` into
    `pristine/bin`.  Rebuilt when the revision's commit or the OpenFOAM
    installation differs.  The main branch's later 71bf192 links xGEMS
    unconditionally from its author's installation and cannot be built
    here; without xGEMS it equals 0.1 (`git diff --stat 0.1 71bf192`:
    `Make/options`, the xGEMS start-up block of `rhoFixedFlowFoam.C`,
    `applications/xgemsFoamProbe`).
  - `dynamicCode/`: the coded `T`/`U` boundary conditions of the 007
    carrier, compiled once by a seed run and linked into every test case.
    The seed is repeated when the inputs of 008-011 change.
  - `regression/`, `frozenCarrier/`, `build/`: the test cases and builds.
    The cached pristine runs of the gate are keyed on the pristine build,
    the mesh, the checksum of every input taken from `run/008`-`run/011`
    (including the linked carrier of `run/007/261`) and the set-up
    functions, so a changed input re-runs them.
  - `results`, `log.Allwmake`: the result lines and the solver build log.
- test cases: copies of `run/` cases (`run/008`-`run/011` for the gate) with
  links to the carrier of `run/007-variable-properties-high-temperature/261`
  and to the mesh.

| script | purpose |
|---|---|
| `tests/Alltest` | driver: `M0`, `M1`, `M2`, `M3`, `M4`, `M4long`, `M5`, `M7`, `M7long`, `M9`, `all` (`M0` to `M5`, `M7` and `M9`, without `M4long` and `M7long`) |
| `tests/functions` | shared shell functions (paths, mesh, cases, builds, results, the load sampler of a timing run that never outlives its script) |
| `tests/regress.sh` | filtered log comparison: complete logs, or `-excerpt` against a case's `out_excerpt`; `-noMockReport` leaves out the per-step report of the coupled mock (the excerpts of 008-010 were recorded before tag 0.1 added it), `-ignore <regex>` the lines of a pattern |
| `tests/cmpFields.sh` | byte-for-byte comparison of the written time directories (including `processor*/`) |
| `tests/regression/Allrun` | the M0 regression gate (008-011 against the pristine build of the tag 0.1 and their `out_excerpt`s) |
| `tests/build/Allrun` | build plumbing: versions, GEMS3K bridge and xGEMS opt-in (a dummy xGEMS), stamps, `gemsConfig.H`, `xgemsConfig.H`, `ldd` |
| `tests/frozenCarrier/Allrun` | M1: `writeFrozenFields`, restarts in serial and on 4 ranks |
| `tests/phaseChange/Allrun` | M2: closed boxes, channels, pipe (serial and 4 ranks), restarts (also re-decomposed, single WALL face, written by function objects, rounded boundary values), start-up and per-step checks, `model mock`, run/012, an unpaired gas in model temperature |
| `tests/phaseChange/checks.py` | the numerical checks of M2 on the written fields and balance files (also the ledger of every later milestone's runs); a balance file or a printed closure with a value that is not finite is refused (plan section 35, item X) |
| `tests/phaseChange/{closedBox,channel}/` | blockMesh templates of the M2 cases |
| `tests/phaseChange/makeCarrier/` | test helper: linear-T plug-flow carrier for the blockMesh cases (built into the work area) |
| `tests/testPhaseChange.C` | standalone g++ test of `phaseChangeKinetics.H` and `compensatedSum.H` (M2.0) |
| `tests/profiles/Allrun` | M3: axial profiles (channel, pipe serial and 4 ranks, synthetic ramp deposits), bins that cut cells, profile times, restarts and writes of function objects, `mDep_`, `exch_`, run/012 with profiles, the restart verdict with recomputed patch values |
| `tests/profiles/checks.py` | the numerical checks of M3 on the profile files, fields and balance files |
| `tests/testProfiles.C` | standalone g++ test of `profileBinning.H`: bin edges, exact overlap weights, onset search (M3.0) |
| `utilities/axialProfiles.py` | reader and plotter of `postProcessing/axialProfiles` (prints T_dep and the integrals; used by the M3 checks) |
| `tests/hks/Allrun` | M5: HKS boxes at 700 K (bare wall, CAPPED, ample condensate, below Y_eq, halfCell), the channel with an inventory sample and its removal (nCorr 3 and 0, binary restarts across the removal, re-decomposed, 4 ranks simple and uneven, the exact clock), run/013 (25 steps and its `out_excerpt`), the tables, refusals, and the carry-overs of the M3 sign-off |
| `tests/hks/checks.py` | the numerical checks of M5 (independent model values: the table, G, lambda, Y_eq) |
| `tests/testVapourPressureTable.C` | standalone g++ test of `vapourPressureTable.H` with synthetic tables and the generated Gurvich and T-Flows tables (M5.4a; with the NASA-9 reference M5.4b) |
| `tests/paperMode/Allrun` | M4: `setFrozenCarrier` on the pipe and a wedge (serial, decomposed, written in parallel; the ascii negative evidence; the T-Flows table; inputs that do not fit the mesh), the layer counts of run/014, layer distance on the channel (serial, 4 ranks split by rows and with three ranks per face, restarts, profiles, HKS in serial and on 4 ranks, refusals) and on the pipe of run/014 (three 4-rank splits; the closure without drift at large steps), run/014's `out_excerpt`, provenance, the v2606 compile |
| `tests/paperMode/longRuns` | M4long: the runs of acceptance 3-5 and the T-Flows parameter set as edits of run/014 (reused when their key is unchanged: the sources, the inputs of run/014, the mesh checksum, the set-up functions, the run's spec and the OpenFOAM installation; run again when it changed, or, only with `M4LONG_REVALIDATE=yes`, revalidated against the binaries they were made with by a 5-step probe from their last write (fields, profiles and the balance file but its booking columns and the diagnostic `evaporatedAboveWarning` identical) and then named so in every criterion line; resumed when stopped; built from the current sources before their binaries are kept; `M4LONG_RUNS`, `M4LONG_JOBS`, `M4LONG_RANKS`), criteria M4.3-M4.5b, M4.16, M4.17 and M4.20, NOTRUN with the cost for runs that did not end; a recorded run of another key that is not selected is kept (NOTRUN, naming the way to revalidate it or to run it again; before the M7 integration an evaluation discarded it) |
| `tests/paperMode/checks.py` | the numerical checks of M4 (continuity lines, carrier values, layer counts, mDep_, inventories, faces, plateau, solid profiles, T_dep with its bound, provenance, the closure of segmented runs, balance files without their booking columns, the explicit remainder and the transport identity of the debug report; `sourceFormMutant` edits a copy of `matrixDefect.H` into the negative control of M4.18) |
| `tests/paperMode/wedge/` | blockMesh template of the 5-degree wedge (40 x 10 cells) and its `setFrozenCarrierDict` (M4.1c, M4.1g, M4.1h, M4.1i) |
| `tests/paperMode/corners/` | an L-shaped channel of 75 cells with three interface patches meeting at a concave and a re-entrant corner, a carrier at rest (makeCarrier) (M4.13, M4.14; the tetrahedral pipe of M4.14 uses its dictionaries) |
| `utilities/setFrozenCarrier/` | the paper-mode carrier (M4; section "Paper-mode carrier" above) |
| `thermochemistry/systems/make_pv_table.py` | writes `data/pv_PbI2_phases.csv` and `data/pv_PbI2_TFlows.csv` from local sources (the T-Flows source at a pinned commit); `--check`, `--reference`, `--sources` (M5.4b) |
| `tests/gems/Allrun` | M9: a copy of the solver built with the GEMS3K bridge into the work area (the default solver stays bridge-free); the channel of M5 (1000 -> 400 K, the boat) with `equilibrium GEMS` in mode local (serial, 4 ranks, speciation lagged, every step written, `updateInterval 3` and its restarts, a fault-injecting bridge, 1e9 iterations per call on 4 ranks, a restart, log directories outside `processorN`) against the table run, the 1150 -> 450 K channel in mode frozen (serial, restarted, 4 ranks, with the fault-injecting bridge), paper mode, start-up refusals and warnings, run/015 (also with a stand-in for another bridge); the closure and the solver defect of every run, and stand-in ledgers with values that are not finite |
| `tests/gems/checks.py` | the numerical checks of M9: p_eq of every element against the table interpolated in Python, the guards recomputed from the fields, the step lines and the summary of the log, the rank logs and their directories, the deposit profiles, the log of the fault-injecting bridge (both modes), the schedule of a restart, the 64-bit totals |
| `tests/gems/faultBridge.c` | a GEMS3K bridge that forwards every call to the real one (dlopen) and returns `GEMSB_ERR_FATAL` at chosen calls, reports p_g = 0 or log10 Omega = NaN after chosen converged calls, a fixed number of IPM iterations, or every p_g scaled by a factor (a stand-in for another bridge), logging the engines (by their creation number), the bounds, the warm states and the statuses (M9.5, M9.8, M9.13; loaded through `LD_LIBRARY_PATH` over the RUNPATH of the bridge build) |
| `tests/build/dummyBridge.c` | the dummy GEMS3K bridge of the build tests: every function of `gemsbridge.h` without GEMS3K (M0.3; since M9 the solver calls the whole API) |
| `tests/verification/graetzRobin` | M7 (acceptance 1): the Graetz-Robin verification: blockMesh wedges of 32 to 256 radial cells with exactly axisymmetric metrics, a parabolic carrier with exact face fluxes, wallResistance none and halfCell at Bi 1 and 50, Gauss linear (both outlet conditions) and SuperBee, a time-accurate control; GCI, observed orders, the discrete mode of the documented scheme (runs reused while the solver, the helper, the template, the scripts and the shell functions that set the runs up are unchanged; `check` evaluates only) |
| `tests/verification/graetz.py` | the reference of M7.1 (the extended Graetz-Robin eigenvalue by a power series in mpmath, DOP853 shooting and Chebyshev collocation), the fit window from the spectrum, the discrete predictor of the developed mode of the wedge scheme, the GCI of Celik et al. (2008) |
| `tests/verification/graetzChecks.py` | the numerical checks of M7.1 (reference, meshes and carrier, set-up from the log, fit, steady state, discrete consistency, fit window and measure, outlet, GCI, the summary table); a fit or a closure that is not finite fails, naming the run |
| `tests/verification/graetzCase.py`, `graetzWedge/` | the case template of M7.1 (`@...@` placeholders) and its filler with the fixed parameters of the study |
| `tests/verification/graetzCarrier/` | test helper: the parabolic carrier of the wedge in its radius y with exact face fluxes (built into the work area; not `setFrozenCarrier`, whose polar radius would bias mu by up to 2.3e-4 on these wedges) |
| `tests/verification/pipeMeshes` | M7long (M7.2a-b): the h/2 first-layer mesh (`refineWallLayer`) and the axial x2 mesh (`refineMesh` along the axis) of the test mesh, with checkMesh and their facts |
| `tests/verification/meshChecks.py` | its checks and the table of the mesh facts |
| `tests/verification/pipeMeshFacts/` | test helper (built into the work area): the interaction layers (faceCells; the cells within 1e-4 m of the wall), the first-layer thickness, the half-cell distance, the axial cells |
| `tests/verification/longStudies` | M7long: the long runs of M7, managed as those of M4long (keys, kept binaries, reuse, resume, `M7LONG_REVALIDATE=yes`, `M7LONG_RUNS`, `M7LONG_JOBS`, `M7LONG_RANKS`; `longStudies keys` prints the keys): run/012 and run/013 on the 007 carrier at dt 4, 2 and 1 ms to 16 s, and 11 edits of run/014 to 30 s (h vs h/2, axial x2, dt 4, 2 and 1 ms); criteria M7.3a-d, M7.2c-f, M7.5a-c, M7.3pa-pc, and its own mechanics M7.6a-b (the verdicts of the criteria, a resume) |
| `tests/verification/longChecks.py` | its checks: the runs (ended, closure in every step, solverDefect), the v_d of the flux form against the solver's print-out, the h vs h/2 measures (onset shift in native bins, peak `mDep_`, profile L1), the axial x2 measures, the intermediate-time error of a time-step study, a byte-for-byte comparison of profiles; a crash reported as `CRASHED` with exit status 3 |
| `tests/verification/dtChecks.py` | the measures of the time-step studies (L1 of the wall deposit between two dt relative to the deposit, peak, centroid, T_dep of both search modes, observed order and Richardson estimate; history over the profile times; a plot) |

Criteria:

| id | criterion |
|---|---|
| BUILD | `wclean` + `./Allwmake` of the solver: every source of `Make/files` compiled, 0 warnings |
| SELF | a test script that stops early fails `Alltest` (self-test of the completeness check) |
| CACHE | the cache keys of the pristine runs and the dynamicCode seed change with a copied dictionary and with the linked carrier of `run/008`-`run/011` |
| LIST | every criterion id recorded (PASS or NOTRUN) is in the expected list of the invocation, with a stand-in id as its control (plan section 35, item U) |
| M0.0 | `testDiffusivity.C` (g++, 0 warnings) passes |
| M0.1a | 008, 009, 010, 011 (10 steps): complete filtered logs identical to the pristine build of the tag 0.1 |
| M0.1b | 011 (`run/011-coupled-species-sources`, the coupled mock): the first block of `out_excerpt` reproduced (519 filtered lines, 10 steps; before tag 0.1 this was 010's, 402 lines) |
| M0.1c | 008, 009: `out_excerpt` (recorded before tag 0.1), without the per-step report of the coupled mock, differs from both builds only in the banner line renamed in e2ddbfa |
| M0.1d | 010: `out_excerpt` (recorded with the mock before tag 0.1), without the report of the coupled mock and without the lines of `Y_PbI2_g` and `Y_PbI2_s`, which the coupled mock changes, reproduced by both builds: 150 lines, the start-up block and every `Y_tracer` line |
| M0.2 | 010, writeInterval 3: all 18 fields (the 15 of before tag 0.1 and the 3 source fields of the coupled mock) and `uniform/` of every write time cmp-identical to the pristine build (60 files) |
| M0.3a-e | `Allwmake` without and with `LESTO_GEMSBRIDGE` (wclean on the toggle, `ldd`, `-help`), plain `wmake` toggle through `gemsConfig.H`, no-op rebuild, accepted versions; 0 warnings |
| M0.3f | an edit of `gemsbridge.h` recompiles `gemsEquilibrium.C`, and only it, with a plain `wmake` |
| M0.3g | `Allwmake` runs `wclean` when the OpenFOAM installation differs from its stamp, or objects (of the solver or of the utility) exist without it, and only then |
| M0.3h | `GEMS_INC`/`GEMS_LIBS` in the environment never reach a build without the bridge (`Make/options` assigns both) |
| M0.3i | a relative `LESTO_GEMSBRIDGE`: `Allwmake` resolves it against the directory it is called from, a plain `wmake` refuses it with a clear message |
| M0.3j | `Allwmake` with `LESTO_XGEMS` (a dummy xGEMS installation): `wclean` on the toggle, 0 warnings, `ldd` lists its `libxGEMS.so`, `-help` 'compiled in', the start-up call in the executable |
| M0.3k | a plain `wmake` after unsetting `LESTO_XGEMS` recompiles `rhoFixedFlowFoam.C`, and only it, through `xgemsConfig.H`: no `libxGEMS`, `-help` 'not compiled in', no start-up call, 0 warnings |
| M0.3l | a relative `LESTO_XGEMS`: `Allwmake` resolves it, a plain `wmake` refuses it with a clear message |
| M0.4 | compiles on v2606: `NOTRUN` unless the tests run under v2606 |
| M1.2a | `writeFrozenFields no`: write directories hold only `Y_*`, `c_*`, the source fields `source_*` of the coupled mock (case 010 is `model mock`) and `uniform/` |
| M1.2b | their size: at most 14 MB and a third of `yes` (see below) |
| M1.2c | `no` and `yes` write cmp-identical `Y_*`, `c_*` and `source_*` |
| M1.3a-c | restart from 0.024, from a directory holding only `Y_*`, cmp-identical to the continuous run at the final time: serial, 4 ranks `simple (4 1 1)`, 4 ranks `simple (2 2 1)` |
| M1.3d | the same, serial, from 0.036, where the restart's time value (from the directory name) differs from the continuous run's accumulated one |
| M1.4 | every `no` log prints the frozen instance; the `yes` log does not |
| M1.5 | `no` with a carrier copy in a later time directory (written by a `yes` run): the log warns that the earlier directory holds `rho` too; ordinary `no` runs do not |
| M2.0 | `testPhaseChange.C` (g++, 0 warnings): activation, half cell, release overlap, the B2 bound rule, the E4 predictor with CAPPED, compensated sums independent of the order, release totals of long and late windows on the solver's clock (compensated, represented length) to 4 eps |
| M2.1a | closed box (10 cells of 57.46 um, T = 500 K, A = 220, dt = 4 ms): Y^n = Y0/(1 + lambda dt)^n, lambda dt = 0.671504, in every cell to 1e-14 |
| M2.1b | closed box: closure <= 1e-15 in every step |
| M2.1c | closed box with a negative pulse (review B2): the deposit made in the first step is kept when its cell turns negative; no element CAPPED; the monitor prints no negative max mole fraction |
| M2.1d | closed box, `wallResistance halfCell`: the closed form with the series velocity to 1e-14; the start-up lines give the kinetic rate and the applied rate `v'/h` to 1e-12 |
| M2.1e | closed box, wall 450 K and cells 500 K: `interfaceTemperature wall` follows the closed form with f(450 K), `cell` the one with f(500 K), to 1e-14 |
| M2.1f | closed box, one corrector and a negative pulse next to the positive cell: the guard switches the elements whose gas turns negative to NONE and re-solves; no CAPPED, 0 projections, closure <= 1e-15 |
| M2.1g | inclined box (cell heights 42 to 78 um), `areaPerVolume layer`: every cell follows the closed form with lambda = A f to 1e-14; the control with `cell` does not |
| M2.1h | closed box, `nGuard 0`: significant projections and a significant complementarity violation (missed deposition) in the first step, closure <= 1e-15, no deposit ever decreases |
| M2.1i | closed box, `depositionVelocity` = A h instead of A: the closed form to 1e-14 |
| M2.1j | closed box with a condensate `molarMass`: at every write time `Y_PbI2_s` = m'' A/(rho V) and `c_PbI2_s` = rho `Y_PbI2_s`/M |
| M2.1k | closed box with a relaxation factor 0.5: the warning that it is ignored, the closed form to 1e-14 |
| M2.1l | closed box with `areaMultiplier 2`: the applied rate 440 1/s printed, the closed form with lambda = 2 A f to 1e-14 |
| M2.1m | the same box: 2 x `surfaceFieldValue areaIntegrate` of `mw_PbI2_s` on WALL = the ledger's wall inventory to 1e-14 in every step (`mw_` is per effective element area) |
| M2.1n | closed box, `wallResistance halfCell` with `D 0`: no floating-point exception, `v' = 0`, Y stays Y0 exactly, closure <= 1e-15 (review round 4) |
| M2.1o | closed box whose interface patch is named `wall` (no `WALL`): runs (before: "boundary patch WALL was not found"), no `Cells adjacent to WALL` line, the closed form to 1e-14, closure <= 1e-15 (review, final M2 pass) |
| M2.2a | channel 1000 to 400 K, release over 40.5 steps from inside a step: released = n0 M to 1e-14 |
| M2.2b | channel, loose solves: closure <= 1e-13 in every step, solverDefect reported (non-zero) |
| M2.2c | channel, `Y_PbI2_gFinal` 1e-12: solverDefect <= 1e-11, closure <= 1e-13 |
| M2.2d | channel from t = 10 s, `deltaT` 0.22 ms, `timePrecision` 6: the release of every step = n0 M overlap/length on the exact step times to 1e-12; total n0 M to 1e-14 |
| M2.2e | channel, outer tolerance 0.05: every step converged, deciding corrector and final predictor reported, closure <= 1e-13, 0 significant events |
| M2.2f | channel fed through its inlet only: printed relative closure <= 1e-13, 0 significant events |
| M2.2g | channel, release window 6 s at 1 ms (6100 steps, as run/012): released = n0 M to 1e-14, closure <= 1e-13 in every step |
| M2.2h | channel at M8 scale, window 15 s at 0.22 ms, 30 s (136364 steps, run in the background): released = n0 M to 1e-14, closure <= 1e-14 in every step (until the M4 review, round 2, the M8 gate 1e-12: the closure drifted linearly to 1.6e-13 while the gas stood steady; plan section 27, item 1) |
| M2.2i | channel from t = 1000 s, window 1000.0001 s over 0.0198 s: per-step release = n0 M overlap/length to 1e-12, total n0 M to 1e-14 |
| M2.2j | `surfaceFieldValue` integral of `mw_PbI2_s` over WALL every step = the ledger's wall inventory to 1e-12 (`mw_` refreshed every step) |
| M2.3a | pipe (run/012, 007 carrier): prints `A_I/V_I = 17681.9` and `v_d = 0.0124421` |
| M2.3b | pipe, 25 steps: closure <= 1e-13; min Y and negative mass in every step; 0 significant projections and complementarity violations |
| M2.4a-c | pipe on 4 ranks, simple (4 1 1), simple (2 2 1), hierarchical (2 2 1) (ranks with 3, 2, 2, 3 processor patches, review B1): closure <= 1e-13, gas and wall inventories = serial to 1e-9 in every step, `gSum(fvMatrix::residual())` printed and not usable |
| M2.4d | channel on 4 ranks, simple (4 1 1), uncollated and collated, with `uniform` 10-digit processor-patch values of the decomposed carrier: closure <= 1e-13 in every step |
| M2.4e | pipe with run/012's own outer loop (nCorr 10, tolerance 1e-5, converging), 4 ranks simple (4 1 1) and hierarchical (2 2 1): gas and wall inventories = serial to 1e-9 in every step, closure <= 1e-13 |
| M2.4f | channel with its own outer loop (nCorr 3, tolerance 1e-7, exits unconverged), 4 ranks vs serial: the difference is printed and below 1e-5 (1e-9 needs a converged loop) |
| M2.4g | 4-rank channel: the processor patches of `mw_PbI2_s` hold the neighbour rank's cell values, bit for bit |
| M2.4h | channel, `wallResistance halfCell`, 4 ranks simple (1 4 1), two ranks without WALL faces: ends within 300 s (before: a deadlock at start-up) and closes to 1e-13 (review round 4, blocking) |
| M2.5a/b | channel, binary restart from 0.1, serial and 4 ranks: every file at 0.2 cmp-identical (review B5) |
| M2.5c | channel, ascii restart without `writePrecision`: the restart term is non-zero and printed, closure <= 1e-13 |
| M2.5d/e | FatalError: state without the ledger; an exact (binary) restart whose deposits (`uniform/phaseChange/phaseChangeDeposits`) belong to another time (term above 1e-14) |
| M2.5f | B5 negative control: the converging channel restarted with `phaseChangeRegimes` is cmp-identical, without it it differs |
| M2.5g | restart from a rounded time name (10.0108 for t = 10.01078): clock warning, every ledger row identical to the continuous run |
| M2.5h | ascii restart of the box with foreign deposits: counted as rounded, the warning above 10^(1 - writePrecision), booked as restart, closure <= 1e-15 |
| M2.5i | FatalError: the ledger with neither `uniform/phaseChange/phaseChangeDeposits` nor `mw_` |
| M2.5j | pipe (Y 0.1) restarted from serial onto 4 ranks: deposits re-read from `mw_` (exact), restart term <= 1e-15 (compensated sums), regimes and deposits of the serial elements not used |
| M2.5k | 4-rank channel re-decomposed with reversed slabs (same element counts): regimes and deposits not used (fingerprint), deposits from `mw_` |
| M2.5l | `-fileHandler collated`, binary, 4 ranks: restart from 0.1 reported binary with term 0, `processors4/0.2` cmp-identical; also with `data.format` removed from the container headers |
| M2.5m | `-fileHandler collated`, ascii: restart reported ascii (not held against the binary bound), term booked, closure <= 1e-13; also without `data.format` |
| M2.5n | restart with an `mw_` of another time (uniform/ intact): warning, and cmp-identical at 0.2 (`mw_` is output) |
| M2.5o | 3 ranks, one with a single WALL face and one joined by a single processor face: the binary files hold `uniform` 6-digit values, yet the restart is exact and cmp-identical at 0.3 (review round 3, blocking) |
| M2.5p | directories written by `runTimeControl` between two steps (`writeNow`, `writeAndEnd`) hold the ledger and the element state; restarts from both are cmp-identical at 0.2 (review round 3) |
| M2.5q | restart without `uniform/phaseChange/phaseChangeDeposits`: warning, deposits from `mw_` (exact), cmp-identical at 0.2 |
| M2.5r | FatalError: a ledger layout (`uniform/phaseChangeLayout`) with the transport patches in another order |
| M2.5s | binary restarts with `writeFrozenFields yes`: the box with a uniform `rho` (`uniform 0.0975561`) is `rounded`, names the carrier and its rounded patches, books the rounding and closes to 1e-15 (before: FatalError); the channel (its uniform patch values came from 6-digit text) is `exact`, term 0, cmp-identical at 0.2 (review round 4) |
| M2.5t | FatalIOError naming the change: restart with `areaMultiplier 2` (binary) and with `areaPerVolume layer` (ascii; before: booked silently as rounding) |
| M2.5u | FatalIOError: at a restart inside the release window the amount changed, or the sample removed; a new sample whose window began before the restart |
| M2.5v | accepted: the sample removed after its window (cmp-identical at 0.2 but the layout); a new sample from 0.15 s (released = (n0 + 1e-6 mol) M to 1e-14, closure <= 1e-13) |
| M2.5w | a layout of an earlier build (no interface, samples or uniform values): two warnings, cmp-identical at 0.2 |
| M2.5x | an inlet value 1.234567891e-3 (written `uniform 0.00123457`): the restart reports `rounded for the next steps: Y_PbI2_g ('uniform' on INLET); inventory exact, continuation rounded`, term 0, and the files at 0.1 differ, as stated (before: `exact`); an inlet 1e-3: `exact`, and `possibly rounded ...` without the record `uniformValues` (review, final M2 pass) |
| M2.5y | `writeFrozenFields yes`, carrier written at 12 digits, the run at 6: the restart names `rho ('uniform' on INLET)`, term 0, and the files at 0.1 differ, as stated (before: `exact`; review, final M2 pass) |
| M2.6a-e | FatalError: `backward` ddt; fixedValue WALL; final solve with relTol > 0 (two cases); sample overlapping the layer; pair without molarMass |
| M2.6f | the models and layers of later milestones run or stop with their own messages: `model HKS` without a table stops ("vapourPressure"), `mode inventory` with `model temperature` stops; the closed box with `layer distance; distance 1e-4` runs with its 10 cells, and `layer distance` without `distance` stops with "Entry 'distance' not found" (plan section 25, item 8) |
| M2.6g | FatalError: fresh start after the release window has begun |
| M2.6h | FatalError at a later step: `fvSolution` changed to relTol 0.1, and `fvSchemes` changed to the ddt scheme `backward`, during the run |
| M2.6i | start-up warning for an outer tolerance below the solver tolerance (pipe), none for run/012 |
| M2.6j | solvers given by patterns: `"Y_.*"` alone gives no separate Final solve (no final predictor; before: one redundant solve per step), `"Y_.*"` then `"Y_.*Final"` gives one; both follow the closed form to 1e-14; `"Y_.*Final"` before `"Y_.*"` is shadowed: FatalIOError with that hint |
| M2.6k | FatalIOError: the convection scheme `bounded Gauss upwind`, at start-up and when `fvSchemes` is changed so during the run (review, final M2 pass) |
| M2.6l | FatalIOError: `wallResistance halfCell` with `areaMultiplier 2`, and with `areaPerVolume layer` (review, final M2 pass) |
| M2.6m | FatalIOError: `removeTime`, `areaPerVolume`, `kineticScale`, `Ce` in a release sample (review, final M2 pass) |
| M2.7 | `model mock` present: 008-011 logs and 010's 60 written files identical to the pristine build (tag 0.1) |
| M2.8 | run/012 reproduces its `out_excerpt` (5 steps), without the line `Cells adjacent to WALL: 15264`, which is printed for `model mock` only since the final M2 pass (run/012 is not edited by the tests; the line may be removed from its `out_excerpt`) |
| M2.9 | an unpaired gas in model temperature (run/012 with the tracer of run/010, 5 steps): the tracer, solved by `solveGasSpecies.H` with a zero source, has the solver lines and the field `Y_tracer` of the pristine build (tag 0.1, model mock there); the balance file of `PbI2_g` equals run/012's without the tracer; no report of the coupled mock, no source field written |
| M3.0 | `testProfiles.C` (g++, 0 warnings): bin edges (widths of 1e5 to 1e12 m over a range of 0.1 m give one bin); overlap weights inside, across, outside, of zero extent, and for 20000 random intervals (weights = overlap/length to 1e-13, binned + outside = total to 1e-13); onset search on a ramp, with a spike, with searchFrom and in the edge cases |
| M3.1a | channel: sum n' dx M + outside = the gas, wall and sample inventories of the ledger at the same exact time to 1e-12, native (100 bins) and 1 cm bins, every profile file; the header inventories equal the balance file exactly |
| M3.1b | pipe, serial and 4 ranks: the same for the native (424) and 1 cm bins |
| M3.1c | channel, bins [0.015, 0.0875] m of width 0.01 (outside the mesh range, last bin 2.5 mm): binned + outside = inventories to 1e-12 |
| M3.1d | bins that cut cells ([0.0003, 0.1] m, width 2.5 mm) on the channel and the ramps: every bin's gas, wall, sample and observable = the exact overlap of the cell-aligned native profile to 1e-13 of the column's maximum, Twall and wallArea likewise (review of M3, round 1) |
| M3.2a/b | pipe on 4 ranks, `simple (4 1 1)` and `hierarchical (2 2 1)`, with the outer loop converged (correctors 1e-13, relTol 0, outer tolerance 1e-12): every profile column within 1e-12 of its peak; of every onset search the status and xPeak equal, T_dep within 1e-9 K, x_dep within 1e-12 m, peak and threshold within 1e-12 of the observable's maximum, of serial |
| M3.2c | the same on 3 ranks of an uneven manual split along x (x < 0.1 m, 0.1-0.25 m, >= 0.25 m; carry-over (c) of the M3 sign-off, implemented in M5; `tests/Alltest M5` runs it as well, `tests/profiles/Allrun M3.2c`) |
| M3.3a | synthetic linear-ramp deposit, `wallPlusGas`: T_dep fromInlet = fromPeak = the analytic crossing 814 K to 1e-9 K (native bins) and the analytic value of the 5 mm bins |
| M3.3b | the ramp with an upstream spike: fromInlet finds the spike (942.8 K), fromPeak the ramp (814 K), both analytic to 1e-9 K; `searchFrom 0.02` makes fromInlet find the ramp |
| M3.3c | the ramp with the spike, observable `wall`: the spike is gas only, both modes find the ramp to 1e-9 K |
| M3.4a | channel: profiles at the write times and at every `times` entry (two entries on one step give one directory; an entry after the end is not reached); time directories of the case = the write times only (review B4) |
| M3.4b | channel restarted from 0.05: only the later entries are written, byte-identical to the continuous run; the start-up names the 3 entries written by the run that wrote the restart state (its record) |
| M3.4c | 4-rank pipe: the profile time 0.012 written, no time directory 0.012 in the case or any processor directory |
| M3.4d | channel written only by `runTimeControl` between two steps (`writeNow` at 0.056, `writeAndEnd` at 0.057, profile time 0.057), serial, 4 ranks, and 4 ranks collated: profiles at 0.056 (`write time`) and 0.057 (`write time; profile time 0.057`), byte-identical to `controlDict` writes at the same steps (serial) and between collated and uncollated; integrals to 1e-12 (review of M3, round 1) |
| M3.5a | box, `areaMultiplier 2`: `mDep_` = 2 `mw_` exactly on WALL; areaIntegrate of `mDep_` = the wall inventory to 1e-14 every step |
| M3.5b | inclined box, `areaPerVolume layer`, non-uniform gas: areaIntegrate of `mDep_` = the wall inventory to 1e-14; that of `mw_` misses it by more than 1e-3 |
| M3.5c | box, `areaMultiplier 1`: `mDep_` = `mw_` exactly; the processor patches of `mDep_` hold the neighbour cells (4-rank pipe) |
| M3.6 | box, `writeExchangeField yes`: `exch_PbI2_g` = rho (Y^n - Y^n-1)/dt to 1e-12 and sum exch V dt = -(change of the wall inventory) to 1e-12; no `exch_` without the switch |
| M3.7 | run/012 with profiles at 0.002 and 0.004: the log without the profile lines reproduces its `out_excerpt` as is; integrals to 1e-12; no time directory; `utilities/axialProfiles.py` prints T_dep |
| M3.8 | restart of a box whose single-face INLET (`fixedGradient`) and OUTLET (`mixed`) values are rounded `uniform` values recomputed when read: not recorded as rounded, `exact` with and without the writer's record, the verdict's coverage stated, continuation cmp-identical |
| M4.1a | the paper-mode carrier of run/014 on the pipe (serial): `max |sum_f phi|` per cell <= 1e-15 of the inflow (2.2e-18); U_b = Q/(pi R^2) to 1e-9; the discrete inflow, the facet-area ratio and the discrete bulk velocity printed; T solved (Picard converged), the wall faces holding the table at their centres exactly |
| M4.1b | the same carrier decomposed by `decomposePar`, simple (4 1 1) and hierarchical (2 2 1): `setFrozenCarrier -check -parallel` <= 1e-15 (plan section 15, item 7); written by `setFrozenCarrier -parallel` on the decomposed mesh: <= 1e-15, its `phi` equal bit for bit to the decomposed serial one on every rank |
| M4.1c | a wedge (blockMesh, 5 degrees, 40 x 10 cells, the meshes of M7): <= 1e-15; the sector from the wedge patches; `rescaleToFacets yes`: the discrete inflow Q 5/360, summed exactly from the written `phi` to 2.5e-13 with the sector from the wedge-plane normals and to 5e-17 (limit 1e-15) with `sectorAngle 5` (review of M4, round 3); the facet-area ratio sin(theta)/theta; a `columns` table on the whole cross-section (`mode profile`) |
| M4.1d | the wedge decomposed (simple (4 1 1)): `-check -parallel` <= 1e-15; decomposed with a single INLET face on rank 1: the `uniform` value decomposePar writes for that patch is reported above 1e-15, and the carrier written with `-parallel` on that decomposition is exact |
| M4.1e | negative evidence: the pipe decomposed with `writeFormat ascii` stays exact (both axial faces of a cell carry the same flux and round alike); one processor value rewritten as `decomposePar` writes a `uniform` patch is reported above 1e-15 |
| M4.1f | the T-Flows wall table (count line 66, 70 rows): `rows all` reads 70 and names the mismatch; `rows count` reads 66 (x <= 0.65 m), stops at the faces beyond with `outOfRange fatal` and holds T(0.65 m) there with `hold` |
| M4.1g | the solver on the wedge carrier (temperature model, layer distance, a release near the axis), 20 steps: closure <= 1e-13 every step, released n0 M to 1e-14 |
| M4.1h | a wedge not extruded along the axis (inner points moved radially and axially, `checks.py perturbWedge interior`; checkMesh OK): <= 1e-15 in serial and decomposed, the wall flux at round-off, and the face-centre flux of the profile misses continuity by more than 1e-10 (negative evidence) |
| M4.1i | a wedge whose wall faces are tilted (`perturbWedge wall`): the utility stops and writes nothing; with `permeableWalls yes` it writes with a warning, continuity <= 1e-15, the wall flux reported above 1e-12 |
| M4.1j | the wedge after its run (startFrom latestTime, 0.02): `-check` reads the `phi` of 0 found by findInstance (<= 1e-15); the utility writes the carrier into 0, unchanged, and nothing into 0.02 |
| M4.1k | inputs that do not fit the mesh (the wedge, `rescaleToFacets yes`): an axis pointing out of the inlet, `inlet` naming the outlet, and a radius below the wall each stop the utility with its message and without writing a field; with `allowBackflow yes` the last writes |
| M4.2 | run/014 on the pipe (acceptance 2): layer distance 30528 cells, `A_I/V_I = 8533.7` (to 1e-5), the T-Flows centre-distance rule the same cells; layer faceCells 15264 cells, 17681.9 |
| M4.3 | `M4long`: Fig. 9 identity n' t_ev = n0/U_b,d (U_b,d the discrete bulk velocity) on every native bin of the plateau 0.25-0.45 m at t_ev + 1 to 1 %, t_ev 6 and 15 s (one line each), printed against the paper's 3.38e-4 mol s/m |
| M4.4 | `M4long`: the solid profiles at 2 t_ev of t_ev 6 and 15 s agree (a property of the model; Fig. 9 has no 6 s curve): L1 difference of the native wall columns <= 1 % of the deposit |
| M4.5 | `M4long`: Table 2, temperature-based column, with the gate restated by the researcher (2026-10-02): \|T_dep - T_paper\| <= max(2 K, dT_bin), dT_bin the drop of the wall temperature across the native bin of the onset; T_dep fromInlet at 1 % of the peak of the wall column, for every case and both layer rules (10 lines), at the end of the experiment or bounded there through the irreversible deposit (a run that has not ended is decided from that bound at its latest profiles, else NOTRUN); the raw difference, dT_bin and the plain +-2 K comparison are printed.  All 10 lines meet the gate; T5_60 (+2.25 / +2.22 K, dT_bin 2.77 K) is outside the plain +-2 K (plan section 26, item 5) |
| M4.5b | `M4long`: T5_60 at dt 10 and 5 ms: T_dep within 0.5 K (the larger dt of T5_60) |
| M4.6 | the sources (solver and `setFrozenCarrier`) compile and link with 0 warnings against OpenFOAM v2606 (a copy without the build directories of any platform, wcleaned and built in the work area, every source of both Make/files compiled; `LESTO_OPENFOAM_V2606`; NOTRUN without it) |
| M4.7a | channel, layer distance 3.5e-4 m (400 cells, `A_I/V_I = 2500 1/m`), tight solves, 200 steps: closure <= 1e-13 every step; `areaIntegrate` of `mDep_PbI2_s` on WALL = the wall inventory to 1e-12 every step; both walls of the symmetric channel deposit the same to 1e-12 |
| M4.7b | the channel on 4 ranks split by rows (the first row of each wall on its own rank: 200 elements with their nearest face on another rank): closure <= 1e-13, inventories = serial to 1e-9 every step, `mDep_` of every wall face = serial to 1e-9 of its largest value, `areaIntegrate` = the wall inventory to 1e-12 |
| M4.7c | the pipe of run/014 (SuperBee, background Y = 0.1, 10 steps of 4 ms, every corrector solved to 1e-15, 3 correctors per step) in serial and on 4 ranks: simple (4 1 1), hierarchical (2 2 1), and an uneven manual split with the first cell layer on rank 0 (15264 elements with their nearest face on another rank): closure <= 1e-13, inventories = serial to 1e-9 every step, `areaIntegrate` of `mDep_` = the wall inventory to 1e-12 every step |
| M4.7d | binary restarts in distance mode from 0.1 (channel, serial and 4 ranks): every file at 0.2 cmp-identical to the continuous run, no warning that `mw_` differs from the deposits |
| M4.7e | restarts of the pipe of M4.7c at 0.04 onto another decomposition (the serial state on 4 ranks, hierarchical (2 2 1); the uneven split reconstructed and continued in serial): the deposits re-read from `mw_` give the ledger's inventory to 1e-15, `exact` (A_I and V_I are compensated sums; with plain sums 7.3e-15 and 2.8e-15, and 2.0e-14 from 16 ranks to serial, a FatalError) |
| M4.7f | the serial pipe of M4.7c restarted from 0.02 (deposits that do not all survive (m''_e A_e)/A_e): every file at 0.04 cmp-identical, exact (relative 0), no `mw_` warning |
| M4.7g | the channel with layer distance 5.5e-4 m (three rows per wall) on 4 ranks, the three rows of the lower wall on ranks 0, 1 and 2 (each of its 100 faces adds the deposits of two other ranks; printed at start-up): two runs on that decomposition write every file, the balance, profile and `mDep_` files bit for bit alike (the push of the nearest-face map adds the ranks in a fixed order; reviews of M4, rounds 2 and 3); closure <= 1e-13, inventories = serial to 1e-9, `mDep_` = serial to 1e-9, `areaIntegrate` = the wall inventory to 1e-12 |
| M4.7h | maps that run both ways between two ranks (review of M4, round 3), 2 lines: the channel of M4.7a on 4 ranks with rows 0 and 22 on rank 0 and rows 1 and 23 on rank 1, and the pipe of M4.7c on hierarchical (1 1 4): the start-up counts (200 and 161 elements with their nearest face on another rank, 0 faces with elements of two or more other ranks), closure <= 1e-13, inventories = serial to 1e-9, `areaIntegrate` of `mDep_` = the wall inventory to 1e-12 every step (the build of round 2: +50 % and +7.8e-4, 200 and 13 faces counted); the channel's `mDep_` = serial to 1e-9 |
| M4.8 | refusals (7): `halfCell` and `areaPerVolume cell` with layer distance; layer distance without `distance`; a distance that selects no cell; restarts with another distance and with the other rule (both ways) |
| M4.9 | profiles in distance mode (channel): bins + outside = inventory to 1e-12, every file and set; T_dep found; the `Twall` and `wallArea` columns equal to those of layer faceCells bit for bit |
| M4.10 | model HKS with layer distance (the inventory channel of M5, `interfaceTemperature wall` from the nearest faces), 2 lines: serial: closure <= 1e-13, 0 significant projections and complementarity violations, the wall deposits, one removal; on 4 ranks split by rows (200 elements pull the temperature of their nearest face from another rank): closure <= 1e-13, gas, wall and sample = serial to 1e-9 every step, `areaIntegrate` of `mDep_` = the wall inventory to 1e-12, the same counts, one removal |
| M4.11 | run/014 reproduces its `out_excerpt` (5 steps) |
| M4.12 | provenance (NOTRUN without the T-Flows repository at 35c1f6406e): the wall table of run/014 is the T-Flows file byte for byte, and the conductivity table of its `setFrozenCarrierDict` is `k_he` of `User_Mod/Types.f90` |
| M4.13 | layer distance with several interface patches (`tests/paperMode/corners`, an L-shaped channel; the patches in both orders of the boundary file), 2 lines: distance 7.5e-5 m, 5 cells (the corner cell of WALLB and WALLC included), `A_I/V_I = 20000 1/m`, `mDep_` 0 on WALLA and WALLC and > 0 on WALLB, the same in both orders; 1.25e-4 m, 15 cells (the cell at the re-entrant corner included), 6666.67 1/m, `mDep_` the same in both orders, 0 on the corner face of WALLC |
| M4.14 | the layer independent of the decomposition, 2 lines: the corners on 3 ranks, one block each (the rank of the cell at the re-entrant corner has no interface face): 15 cells, inventories and `mDep_` = serial to 1e-12, restarts serial -> 3 ranks and 3 ranks -> serial exact; a tetrahedral pipe (gmsh; NOTRUN without it): the same layer in serial, on 4 ranks and recomputed from the mesh (`checks.py exactLayer`), the serial state restarted on 4 ranks exact, `areaIntegrate` of `mDep_` = the wall inventory to 1e-12 in serial and on the 4 ranks (maps both ways; the start-up count exactly 21 elements with their nearest face on another rank and 0 faces with elements of two or more other ranks, gated since the M7 integration; round 2: +0.92 %, 14 faces) |
| M4.15 | the closure does not drift: the pipe of run/014 with its own solver settings, background Y = 0.1, 20 steps of 0.1 s on 4 ranks: closure <= 2e-15 of the reference in every step (4.2e-16; the binary before the M4 review, round 2: 1.16e-14, growing linearly) |
| M4.18 | the face-flux correction of the transport matrix in flux form (review of M4, round 3): the drift pipe of M4.15 with layer distance and with layer faceCells (4 ranks, `reportMatrixResidual yes`) prints an explicit remainder of 0 kg/s in 0 cells in every step, of an explicit source sum \|b - s\| of 1.4e-6 to 1.8e-6 kg/s (the non-orthogonal correction), and the transport identity of the cell defects, D = sum_c r_c + dM/dt + sum_p F_p - sum_c E_c summed exactly, at most 1e-30 of its terms in every step (measured at most 1.1e-36; plan section 33); faceCells closure <= 2e-15; positive control: the channel with `Gauss linearUpwind` leaves a remainder (its correction is in the source only) and D stays 0; negative control: a copy of the sources whose cell defect keeps the correction in source form while the report of E is unchanged (`checks.py sourceFormMutant`, built into the work area) passes the check of E but violates the identity in every step (1e-21 to 8e-20 of its terms) |
| M4.19 | the mole-fraction monitor in paper mode (review of M4, round 3): the max mole fraction printed for the last step of the serial pipe of M4.7c = max c/(c + p/(R T)), c = rho Y/M, recomputed from its fields, to 1e-5 (0.0675; the old form (Y/M)/(Y/M + 1/M_He) printed 0.00287) |
| M4.16 | `M4long`, not a gate: the positions of Figs. 9-10 (the gas front at 16 s, the solid onset and peak at 30 s) against run/014 (fig9-tev15) and against the T-Flows parameter set (fig9-tev15-TFlows: Tdep 638.15 K, A 150, k 0.02, distance 2.4e-4), printed in mm |
| M4.17 | `M4long`: the closure of every long run made with the current sources, at most 1e-14 of the reference in every step of every balance file (a run restarted in segments has several); a run revalidated from older binaries is listed with the closure of their booking, not gated; NOTRUN when no run was made with the current sources |
| M4.20 | `M4long`: the key of a long run follows the mesh (the checksum of `tests/work/mesh`) and the shell functions that set a run up, as the cache keys of the gate do: a stand-in mesh checksum and a redefined `setupLong` each give another key (review of M4, round 3: neither was in the key) |
| M5.0 | `testPhaseChange.C` with the HKS checks: G R T/M = 89.65 m/s at 700 K, kineticScale 1e-5/sqrt(1000) = the T-Flows expression, the HKS half cell = the Robin flux, the predictor of HKS wall and sample elements, the admissible flow, the predictor with the previous regimes (no jump between CAPPED and NONE), a row of 10 stiff sample cells run as the corrector loop: alternating without the previous regimes, settled and complementary with them; the guard's switch back, largest misses first until the rest is within the tolerance |
| M5.1a | box at 700 K, bare wall, Y0 = 0.1: Y1 = (Y0 + lambda dt Y_eq)/(1 + lambda dt), lambda dt = 6241, and the next steps Y^n = Y_eq + (Y0 - Y_eq)/(1 + lambda dt)^n, to 1e-14 (1.7e-16); closure <= 1e-13 (2e-16) |
| M5.1b | box, a deposit of 5e-8 kg/m2 and Y = 0 (through `mw_`, restart): all elements CAPPED, Y1 = m''/(rho h) to 1e-14, deposit, `mw_` and `Y_PbI2_s` exactly 0, the closure of the step 0 |
| M5.1c | box, a deposit of 1e-5 kg/m2 and Y = 0: IMPLICIT evaporation, Y1 = lambda dt Y_eq/(1 + lambda dt) to 1e-14; closure <= 1e-13 |
| M5.1d | box, bare wall below Y_eq: NONE; `exch_`, deposit and `Y_PbI2_s` exactly 0; Y = Y0 to 1e-14 (5.2e-16: not bit for bit, the linear solver iterates once per step on the uniform field); closure <= 1e-13 |
| M5.1e | box, `wallResistance halfCell`: the closed form with G' (lambda dt = 233) to 1e-14 |
| M5.1f | box, bare wall, 50 steps near equilibrium: closure <= 1e-13 in every step (4e-16; 2.7e-12 with the defect of the assembled matrix), the closed form of Y to 1e-14, the deposit to 1e-10 (its fixed-point drift is booked as solver defect) |
| M5.2a | channel with an inventory sample removed at 0.25 s, nCorr 3: closure <= 1e-13 every step, 0 projections, 0 significant complementarity violations and none above `guardTolerance` x inventory in total (`checks.py totalLeft`), the removal once in the step starting at t^n = 0.25 (exact clock) of exactly the reservoir then, sample 0 afterwards, recorded in the layout, sample decreasing before |
| M5.2b | the same with nCorr 0, acceptance 2 as written: guard passes > 0 and 0 projections; also 0 significant complementarity violations, no step whose complementarity violations exceed `guardTolerance` x inventory in total (recomputed from the balance file, and not flagged by the solver), closure <= 1e-13 (round 2: 0 projections, but significant missed deposition in 56 steps that the guard did not switch back; round 3: two misses of 1.43 times the tolerance together accepted without a guard pass; round 4: 2 steps with projections below the tolerance, left by guard re-solves of a re-assembled SuperBee matrix, which passed only under the proposed rewording '0 significant projections', now withdrawn) |
| M5.2c | binary restarts from 0.2, 0.25 and 0.3 (before, at and after the removal): every file at 0.35 cmp-identical, one removal (none after it), restart term 0 |
| M5.2d | 4 ranks, `simple (4 1 1)` and an uneven manual split (288/56/1096/960 cells, a rank without WALL faces, the sample on three ranks): closure <= 1e-13 every step, one removal |
| M5.2e | the channel from t = 10 s with dt = 0.22 ms: the removal in the step whose exact start is the first t_k >= removeTime - dt/2 |
| M5.2f | the serial state re-decomposed onto 4 ranks: the reservoirs re-read from `Cs_PbI2_s` (binary, exact), inventory equal to the ledger to 1e-15 |
| M5.2g | the serial pipe state (run/013 at 0.1) decomposed onto 4 ranks (hierarchical 1 2 2; the boat on all ranks) and reconstructed back to serial: both restarts accepted, exact, closure <= 1e-13 |
| M5.2h | the removed sample dropped at a restart, a second restart, a third with the sample re-added: all accepted and exact, the record kept, the re-added sample empty |
| M5.2i | a release sample added in front of the inventory sample at a restart: the reservoirs and regimes of `uniform/` still fit (no warning), fields cmp-identical to the continuous run |
| M5.2j | the channel of acceptance 2 with the recommended solves (loose correctors, tight Final solve), 500 steps, nCorr 3 and 0: closure <= 1e-13, 0 significant projections, guard passes > 0 with nCorr 0.  Not acceptance 2 as written (M5.2b): 1 and 4 steps end with projections below the tolerance, in steps without guard passes, which the guard leaves to the projection by design (E5.4) |
| M5.2k | a small reservoir (1.7e-8 mol, `Cs` = 0.98 kg/m3), nCorr 2, 3, 10, 11: the same regimes in every step, inventories equal to 1e-3 of n0 M between 2 and 3 and between 10 and 11 (round 2: 0.65), the bounded law in every sample cell recomputed from the fields, 0 significant projections and complementarity violations, none above the tolerance in total, closure <= 1e-13 |
| M5.2l | the same with nCorr 0: with `nGuard 0` the first step reports as many significant complementarity violations as the fields show (34 CAPPED sample cells in supersaturated gas; round 2: 0), and the exchange line and the summary name the total above the tolerance; with `nGuard 3` the guard alone leaves 0 (none above the tolerance in total), the regimes of nCorr 3 in the first step, the inventories within 2e-2 of n0 M of nCorr 3 (7e-3: one corrector per step of the lagged SuperBee transport; until round 4 the guard re-assembled the matrix, which acted as extra correctors) |
| M5.2m | the pipe state of M5.2g decomposed with `writeFormat ascii`, 4 ranks: the restart is accepted and booked as rounded, the boat keeps its 144 cells, closure <= 1e-13 (round 2: refused as 'changed') |
| M5.2n | `-fileHandler collated`: the inventory channel on 2 ranks to 0.2 (ledger and layout plain files, the per-rank lists containers), reconstructed and continued in serial, and re-decomposed onto 4 ranks (collated), both to 0.3 across the removal: the containers named as unreadable for this run (3 warnings each), deposits and reservoirs from `mw_` and `Cs_`, exact, one removal, closure <= 1e-13 (round 3: both stopped with a FatalIOError of the file handler) |
| M5.2o | the geometry record within the point rounding of both runs: the ascii-decomposed state of M5.2m reconstructed and continued on the binary mesh; the boat moved to x = 0.6 m, decomposed in ascii at writePrecision 10 and 12: all accepted with 144 cells, the layout records the points (round 3: all three refused as 'changed') |
| M5.2p | a sample next to a colder wall (WALL at 700 K, `interfaceTemperature wall`, 4.25e-7 mol in the second cell row), tight solves with nCorr 2, 3, 10, 11 and loose ones with nCorr 0, 1, 3, 10: the bounded law of every wall element (`checks.py wallLaw`) and every sample cell recomputed from the fields (0 misses above 1e-12 of the inventory), 0 significant projections and complementarity violations, none above the tolerance in total, closure <= 1e-13, guard passes in some runs (round 4: the guard re-solve re-assembled the SuperBee matrix; a bare wall element at Y_eq cycled between IMPLICIT and NONE, misses up to 1e6 times the tolerance or projections above it, in 6 of the 8 runs) |
| M5.2q | the inventory channel on 4 ranks (rank 1 with a single WALL face) with `-fileHandler masterUncollated` and with `collated -ioRanks '(0 2)'`, restarted from 0.05 on the same ranks: the per-rank lists of `uniform/phaseChange/` are read (no warning, exact), every file at 0.1 cmp-identical to the continuous run (round 4: in `uniform/` the lists resolved to the master's file on every rank; 14 and 9 files differed) |
| M5.3a | run/013 (background Y = 0.1, 25 steps of 4 ms): 0 significant regime changes in the deciding corrector and in the final predictor, 0 significant projections, none above the tolerance in total, closure <= 1e-13, solverDefect <= 1e-11, the sample decreasing strictly, the wall depositing |
| M5.3b | run/013 reproduces its `out_excerpt` (5 steps; regenerated in the M7 integration with the final tolerance 1e-14 of HKS cases, plan section 33) |
| M5.3c | run/013 at its own settings, 25 steps: 0 significant regime changes (deciding corrector, final predictor), 0 significant projections, closure <= 1e-13, solverDefect <= 1e-11, the sample decreasing (the wall stays empty) |
| M5.3d | run/013 with a small reservoir (1e-7 mol), nCorr 9 and 10, 12 steps: every step converged, identical balance files, the bounded law of the boat recomputed from the fields, 0 significant projections and complementarity violations, closure <= 1e-13, solverDefect <= 1e-11 (round 2: nCorr 10 never evaporated the boat, nCorr 9 dumped it) |
| M5.4a | the repository's tables through `vapourPressureTable.H` (`testVapourPressureTable.C`): nodes exact, T_m from the phase crossing 683.0 +- 0.1 K, the T-Flows table at its nodes to 1e-12 (exact) |
| M5.4b | the tables from their sources (thermo.inp, the pinned T-Flows revision; NOTRUN without them): `make_pv_table.py` regenerates the csv files byte for byte, its check passes, the 5 K table vs NASA-9 < 2e-3 decades (4.4e-5) |
| M5.5 | refusals (15): a release key in an inventory sample, unknown units, a phase not in the table, outOfRange fatal, equilibrium GEMS without the bridge, restarts with a new, changed, dropped, mode-changed or moved (same number and volume of cells) inventory sample, a fresh start with `Cs_PbI2_s`, a fresh start at `removeTime`; on an ascii mesh of round coordinates (writePrecision 6) the moved selection and one with a single cell replaced by its neighbour are refused and the unchanged restart is accepted (round 4: the bounds followed the digits the points need and accepted the moved sample) |
| M5.6 | carry-over (a): observable = wall + gas - sampleGas in every file, sampleGas only in the bins of the sample and equal to the gas of its cells |
| M5.7 | carry-over (b): a noise-level wall column (run/012 at 4 ms) is `insignificant` in both searches (the reader agrees), and searched with `significance 0` |
| M5.8 | carry-over (d): restarts with deltaT 0.2 ms and 2 ms of a run written at 0.05 with dt 1 ms: every profile time written exactly once over the two runs |
| M5.9 | the inputs of the other model (`HKSCoeffs` and `vapourPressure` under `model temperature`, `temperatureCoeffs` under `model HKS`) named in a start-up warning; the runs end |
| M7.1a | reference: the decay rate m = mu R of the extended Graetz-Robin problem at Pe 5, Bi 1 and 50 by a power series (mpmath), DOP853 shooting and Chebyshev collocation agrees to 1e-10 (the error bound of the reference) |
| M7.1b | the wedge meshes (Nr 32 to 256, dx = 2 dr): checkMesh OK, non-orthogonality 0, wedges of 1 degree, ratio 2 in r and x; the carrier: max \|sum_f phi\| per cell <= 1e-15 of the inflow, the inlet mean velocity U_b to 1e-12 |
| M7.1c | the set-up of every Gauss linear run from its log: D uniform and equal to the correlation to 1e-15, v_d = k_w, f(T) = 1 on every element; halfCell: the applied rate = v'(d) A_I/V_I to 1e-12, d = R - y_last with y_last the centroid formula (2/3)(r2^3 - r1^3)/(r2^2 - r1^2) of the wall row evaluated in double precision on the unit radius and times R (`graetz.half_cell`; agreement 4.0e-14; the centroid in exact arithmetic would give 3.4e-14), the agreement with OpenFOAM's centroid of the wall row printed too (3.7e-14) |
| M7.1d | steady state of every Gauss linear run: the gas line density in the fit window changes by <= 1e-12 over the last step, the last outer loop converged, no guard pass and no significant projection or complementarity violation, \|closure\| <= 1e-13 in every step |
| M7.1e | discrete consistency: the fitted mu equals mu_h, the decay rate of the developed discrete mode of the documented scheme on the run's mesh (`graetz.py`), to 1e-7, for both wall laws |
| M7.1f | the fit window with one end or both ends moved by -+1.5 R, the mixing-cup concentration and the wall uptake give mu within 1e-7 (the largest change of each kind printed) |
| M7.1g | the outlet fixedValue 0 instead of zeroGradient changes mu by <= 1e-8 |
| M7.1h | acceptance 1: GCI_fine <= 0.5 % on Nr 64/128/256, halfCell and none, Bi 1 and 50 (4 lines) |
| M7.1i | acceptance 1: \|mu_ext - mu_ref\| <= GCI_fine mu_fine + bound mu_ref (4 lines) |
| M7.1j | acceptance 1: the observed order, monotone, within 0.25 of 2 (halfCell) and 1 (none) (4 lines; the orders reported) |
| M7.1k | the time-accurate control (Nr 32, Bi 50, halfCell, Co 1, its final state written) reaches mu_h of its mesh to 1e-9, its profile steady to 1e-12 (the pseudo-transient run printed against mu_h too) |
| M7.1l | SuperBee (Nr 64/128/256): set-up, steady state and mu_fit = mu_h of SuperBee's developed branch (the linear-upwind face value) to 1e-7 |
| M7.1m | SuperBee: M7.1h-j together (4 lines) |
| M7.1n | the checks of M7.1d-g fail on a fit or a closure that is not finite: `graetzChecks.py` discrete, windows, outlet and steady on stand-in runs with mu, a window variant, the wall uptake, the mixing cup or a closure of the log NaN, each refused with 'not finite: <run>'; a checker that crashes at a NaN instead (the control) is told apart (plan section 35, items O and X; no run) |
| M7.1o | the reuse keys of the Graetz runs follow their set-up: a changed `meshCase`, `runCase` or `transientCase` changes the key of a mesh, a run and the time-accurate run (plan section 35, item Q; no run) |
| M7.2a | M7long: the h/2 first-layer mesh (`refineWallLayer '(WALL)' 0.5` of the test mesh): checkMesh OK, one more cell layer, the first layer half of the base's to 1e-9 in every element, A_I/V_I of layer faceCells 1.9 to 2.1 times the base's, one more layer within 1e-4 m of the wall with the same volume |
| M7.2b | M7long: the axial x2 mesh (`refineMesh`, `directions (tan1)`, every cell): checkMesh OK, twice the cells and the axial cells, half the axial extent, the first layer, A_I/V_I and the distance layer's volume unchanged |
| M7.3a | M7long: the six runs of the time-step study on the 007 carrier (run/012 and run/013 at dt 4, 2 and 1 ms to 16 s) end, the closure at most 1e-14 of the reference in every step (runs made with the current sources; a revalidated run is listed) and \|solverDefect\| at most 1e-8 at the end |
| M7.3b | M7long, acceptance 3, temperature model on the 007 carrier: the L1 change of the wall deposit between dt 2 and 1 ms at 16 s, relative to the deposit, below 1 %, native and 1 cm bins |
| M7.3c | the same for HKS |
| M7.3d | M7long, reported: the observed order and the Richardson estimate, the peak, the centroid and T_dep (both search modes) between 2 and 1 ms of both models, and the intermediate-time error of the temperature model while gas is present |
| M7.2c | M7long: the six runs of the h vs h/2 study (106 mL/min, layer faceCells, to 30 s) end, closure and solverDefect as M7.3a; the flux form uses `depositionVelocity 0.0124421` m/s, the v_d that A gives on the base mesh as the solver prints it, and on the base mesh the flux and the volumetric form give the same deposit (onset within 0.01 native bins, peak `mDep_` within 1e-5) |
| M7.2d | M7long, acceptance 2, flux form, wallResistance none: between the base and the h/2 first-layer mesh at 30 s the onset (T_dep fromInlet, 1 % of the peak of the wall column, native bins) moves by at most 1 native bin and the largest `mDep_` on WALL by at most 5 % |
| M7.2e | M7long, reported: the same for the volumetric form (A 220 1/s: v_d = A V_I/A_I halves on h/2), expected to miss them; the line says whether it does |
| M7.2f | M7long, acceptance 2 with wallResistance halfCell (the wall law of physical runs): the limits of M7.2d |
| M7.5a | M7long: the two runs on the axial x2 mesh end, closure and solverDefect as M7.3a, the flux form's v_d |
| M7.5b | M7long, the axial x2 study, flux form, none, against the base mesh at 30 s: T_dep shifts by at most 1 native bin, the 1 cm peak height by at most 5 %, the L1 change of the 1 cm wall profiles is at most 1 % (the criterion proposed by the second M7 study, a gate pending sign-off; `M7LONG_AXIAL_GATE=no` reports it) |
| M7.5c | the same with wallResistance halfCell |
| M7.3pa | M7long: the three runs of the time-step study on the paper-mode carrier (run/014 as distributed at dt 4, 2 and 1 ms to 30 s) end, closure and solverDefect as M7.3a |
| M7.3pb | M7long, acceptance 3 on the paper-mode carrier: L1(2 ms, 1 ms) of the wall deposit at 30 s below 1 %, native and 1 cm bins |
| M7.3pc | M7long, reported: the observed order and the Richardson estimate at 30 s, and the intermediate-time error while gas is present (the largest L1(2 ms, 1 ms) over the profile times with a deposit above 1 % of the release, the values at t_ev and t_ev + 1, the gas at t_ev + 1) |
| M7.6a | M7long: the verdicts of the criteria with a stand-in checker (met, not met, cannot be computed, a crash with python's exit 1, a crash reported with exit 3) in every mode (gate, expectFail, report, axial gated and reported): a crash or 'cannot be computed' is a FAIL; `longChecks.py` reports a crash as `CRASHED` with exit 3 (plan section 35, item N; no run) |
| M7.6b | M7long: a stopped run of the 007 carrier resumes from its latest time: run/012 set up by `setupM7`, stopped at 0.008 s, resumed by `runM7` to 0.016 s (2 steps from 0.008 s, an exact restart; plan section 35, item P) |
| M9.0 | the bridge build: a copy of the sources built with `LESTO_GEMSBRIDGE` into the work area, every source compiled, 0 warnings, `ldd` lists the bridge, `-help` 'compiled in'; the default solver stays bridge-free (2 lines) |
| M9.1 | acceptance 1: without the bridge `equilibrium GEMS` (also with `speciation lagged`) stops with a FatalIOError and the rebuild hint; the bridge build runs the same case (closure, solverDefect) |
| M9.2a | acceptance 2: mode frozen on the 1150 -> 450 K channel, every element with T_e in 500-1100 K within 1e-3 decades of the 5 K table (the written `pEq_PbI2_g` against the table in Python), the elements below 500 K on the table exactly, 152 GEMS nodes without failure or rejection; closure, solverDefect |
| M9.2b | mode frozen is deterministic: a restart cmp-identical to the continuous run, p_eq of every element on 4 ranks equal to serial bit for bit |
| M9.3a | acceptance 3, the 1000 -> 400 K ramp channel in mode local: the source of every element in 30 steps equal to the guards recomputed from the fields, the counts of the log equal to those of the fields, 0 failures; 350 steps with 0 failures and the guarded elements counted, closure <= 1e-13, solverDefect <= 1e-11 (2 lines) |
| M9.3b | acceptance 3: <= 70 us per warm call, with the load average; NOTRUN on a busy machine (load >= 4) when slower; and the load sampler of that run ends with its script, terminated or killed (the loop without a check of its parent that it replaced as the control, bounded to 30 s; the test stops what it started when interrupted) (2 lines) |
| M9.3c | acceptance 3: 4 ranks give 4 engines and 4 logs (`processorN/gemsLog`), the calls of the rank logs add up to the summary's; closure; inventories within 1e-6 of serial |
| M9.3d | acceptance 3: the deposit profile of mode local within 0.3 % of the table run (L1 and the largest bin of the peak; speciation none) |
| M9.3e | acceptance 3: paper mode (rho = 1) refused in mode local unless `allowNonPhysicalCarrier yes`; mode frozen runs on it |
| M9.4 | acceptance 4: `speciation lagged` as a separate run, reported (chi, the deposit and gas against speciation none, T_dep); 0 failures, closure |
| M9.5 | the fault-injecting bridge: in mode local the engine re-created after `GEMSB_ERR_FATAL` (destroyed, created, the condensed species suppressed again, the first call of every element cold), and a converged call with p_g = 0 or log10 Omega = NaN counted as a failure and logged; in mode frozen the engine re-created at a node, and a node with p_g = 0 (logged with its real max log10 Omega) or a NaN Omega failed (2 lines) |
| M9.6 | `updateInterval 3`: updates at the steps 1, 4, 7, 10 only; a restart from 0.004 restores the lagged state and keeps the schedule (the files of 0.005 and 0.006 identical to the continuous run's); without the lagged state the first step updates, then the schedule (3 lines) |
| M9.7 | start-up refusals of the bridge build (16 lines: no gems block, a wrong formula, an unknown species, a gas among the condensates, a carrier that is no gas, speciation lagged with mode frozen and with equilibrium table, no GEMSCoeffs, no engine; the molar mass of the system, a condensed species as the pair's gas, a condensate of another formula, referencePressure outside the grid, an unknown mode, an unknown key of GEMSCoeffs and of a gems block) |
| M9.8 | run/015: its `out_excerpt` reproduced exactly with the bridge that wrote it (the md5 of its `libgemsbridge.so`); with another bridge its start-up and step 1, which no GEMS3K call reaches, exactly and the rest NOTRUN (every round-off token of steps 2-5 moves with the last bits of GEMS3K's p_eq; plan section 35, item V); 25 steps at its own settings with 0 failures, closure <= 1e-13, solverDefect <= 1e-11; the verdicts on a stand-in for another bridge (p_g scaled by 1 + 4.4e-16): NOTRUN, and FAIL with step 1 of the excerpt changed (3 lines) |
| M9.9 | a restart in mode local: exact at its start, closure, the inventories within 1e-6 of the continuous run (the first update after a restart is cold) |
| M9.10 | `wallResistance halfCell` in mode local (the law recomputed with the stored half cell at every update): the deposit within 0.3 % of the table run with halfCell, 0 failures, closure; with speciation lagged as well (the half cell with chi beta): 0 failures, closure |
| M9.11 | start-up warnings: the keys of the other mode named and not read (invalid values of them do not stop the run); equilibrium table naming `GEMSCoeffs` and a `gems` block (2 lines) |
| M9.12 | the log directory of every rank on 4 ranks with `logDirectory` `"<case>/..."`, `"$FOAM_CASE/..."`, `"../..."` and `"../processor0/..."`: every rank its own `gemsEquilibrium.log` and `ipmlog.txt`, and the start-up line names the directory of every rank |
| M9.13 | 64-bit counts: the fault bridge reports 1e9 IPM iterations per call on 4 ranks; every step line and the summary print 1e+09 per call, the rank logs the exact totals (above 2^31) |
| M9.14 | the ledger refuses a value that is not finite: stand-in balance files with the closure NaN in a later row, the solverDefect NaN or the gas inf refused, the finite one met; `logClosure` refuses a NaN closure in a later Balance line (plan section 35, item X; no run, also without a bridge) |

The pipe of M2.3/M2.4 runs 25 steps of 4 ms with a background
`Y_PbI2_g = 0.1`, so that the wall deposits from the first step, and with
a fixed 21 correctors per step: `nCorr 20` and an outer `tolerance 1e-9`
below the corrector tolerance 1e-8, so the loop never converges (the
start-up warns, M2.6i) and serial and parallel runs take the same path.
With a loop that exits unconverged (e.g. the channel's `nCorr 3` with an
unreachable tolerance) serial and parallel runs differ by 1e-7 to 1e-6
(M2.4f): the SuperBee limiter and the non-orthogonal correction are lagged,
and the final solve converges only the linear system of the last iterate.
With run/012's converging loop the pipe agrees to about 1e-10 in every
step, simple (4 1 1) and hierarchical (2 2 1) (M2.4e).  The
converged exit with its final predictor is tested on the channel (M2.2e,
M2.5f), by run/012 (M2.8) and by M2.4e.

The bridge of M0.3b/c is the repository's `thermochemistry/gemsbridge` when
its `libgemsbridge.so` exists (it is git-ignored, so a fresh clone has
none), otherwise a throw-away dummy bridge built in the work area; M0.3f
always uses a dummy.  `LESTO_GEMSBRIDGE_TEST=<dir>` selects another bridge.

In M1.3, `uniform/time` is compared in `name`, `deltaT` and `deltaT0`, not in
`index` and `value` (see Restart above).

The regression filter (`regress.sh`) removes only what identifies a run
(Exec, Date, Time, Host, PID and Case header lines, ExecutionTime lines), the
dynamicCode compilation chatter, and replaces the case directory by
`<case>`.  It keeps the thermophysical selection block and the `Using
dynamicCode` lines, which the prototype filter dropped; hence 402 compared
lines of the 010 excerpt instead of 393 (M0.1b before tag 0.1).  Against
the pristine tag 0.1 (plan section 31) the excerpt of 011 is compared as it
is (519 lines, M0.1b).  The excerpts of 008-010 were recorded before tag 0.1
added the per-step report of the coupled mock: `-noMockReport` leaves out
its lines (`Integrated source ...`, `Net thermochemistry source ...` and the
`Species balance` blocks with `Tracked species total`) in both files (M0.1c,
M0.1d), and `-ignore 'Y_PbI2_[gs]'` the lines of the PbI2 species of 010,
whose results the coupled mock changed (M0.1d: 150 lines compared, 252 left
out, all 93 `Y_tracer` lines among the compared).

The plan estimated "at most 6 MB" for a `writeFrozenFields no` write
directory of 010.  Case 010 writes four fields there (`Y_PbI2_g`,
`Y_tracer`, `Y_PbI2_s`, `c_PbI2_g`); four binary scalar fields on 228,960
cells need at least 7.33 MB, and 7.5 MB were written (36.0 MB with `yes`),
so M1.2b checked at most 8 MB and at most a quarter of the `yes`
directory.  Since tag 0.1 case 010 (`model mock`) also writes the three
source fields of the coupled mock, 3.7 MB more in both (`source_tracer` is
zero everywhere and written as a uniform field): 11.1 MB with `no`, 39.6 MB
with `yes`.  M1.2b checks at most 14 MB and at most a third (plan section
31).

**Verification (M7; plan section 33).**  `tests/Alltest M7` (in `all`)
runs `tests/verification/graetzRobin`, the Graetz-Robin verification of
acceptance 1: the extended Graetz problem (fully developed laminar pipe
flow, axial diffusion included) with the temperature law as a constant
first-order wall uptake (`k (Tdep - T) = 100`, so f = 1 exactly, and
`depositionVelocity` = k_w), R = 2.4 mm, T = 500 K, p = 1 atm, D =
7.6963e-5 m2/s, Pe = U_b R/D = 5, Bi = k_w R/D = 1 and 50, L = 24 R.  The
reference decay rate m = mu R comes from three methods (`graetz.py`); the
meshes are blockMesh wedges of 2 degrees with the vertices at (x, r, -+r
tan 1 deg), so that the cell metrics are exactly axisymmetric, Nr 32 to
256 with dx = 2 dr, and the carrier (`graetzCarrier`) has the exact face
fluxes of u = 2 U_b (1 - y^2/R^2).  mu is fitted to the gas line density
of `postProcessing/axialProfiles` over a window fixed by the spectrum.

    tests/Alltest M7                          # or directly:
    tests/verification/graetzRobin            # all runs, then the criteria
    tests/verification/graetzRobin check      # the criteria of existing runs

Measured on v2412 (25 PASS; every fitted mu equals the first M7 study's
bit for bit), zeroGradient outlet, mu_ref = 107.593854287 (Bi 1) and
273.173651541 1/m (Bi 50):

| scheme | Bi | law | Nr | p | mu_ext/mu_ref - 1 | GCI_fine |
|---|---|---|---|---|---|---|
| linear | 1 | halfCell | 64/128/256 | 2.0010 | +1.1e-9 | 0.00019 % |
| linear | 1 | none | 64/128/256 | 1.0077 | +8.8e-6 | 0.155 % |
| linear | 50 | halfCell | 64/128/256 | 1.9919 | +9.5e-9 | 0.00010 % |
| linear | 50 | none | 64/128/256 | 1.0045 | +1.1e-5 | 0.341 % |
| linear | 50 | none | 32/64/128 | 1.0081 | +3.9e-5 | 0.678 % (above 0.5 %, reported) |
| SuperBee | 1 | halfCell | 64/128/256 | 2.0024 | +4.7e-9 | 0.00042 % |
| SuperBee | 1 | none | 64/128/256 | 1.0143 | +1.6e-5 | 0.154 % |
| SuperBee | 50 | halfCell | 64/128/256 | 2.0104 | +6.4e-8 | 0.0013 % |
| SuperBee | 50 | none | 64/128/256 | 1.0231 | +5.8e-5 | 0.336 % |

halfCell converges at order 2 and `none` at order 1, as the discrete
predictor of the documented scheme says; every run's mu equals that
predictor's mu_h to 7.1e-9 (M7.1e; the GAMG tolerance at the small values
of the Bi 50 window), so the GCI measures discretisation error only.  On
32/64/128 `none` at Bi 50 misses 0.5 %: that triplet is too coarse for the
first-order law (the table of every triplet is in
`verification/graetz/table.md`; `GRAETZ_MAIN="32 64 128"` makes M7.1h
fail for that line).  The time-accurate control (Co 1, its final state
written) reaches mu_h of its mesh to 7e-11 (M7.1k); the fit window moved
by -+1.5 R at one end or at both, the mixing cup and the wall uptake give
mu within 1.2e-8 (M7.1f).  Settings: `GRAETZ_LEVELS`, `GRAETZ_MAIN`,
`GRAETZ_OUTLET_LEVELS`, `GRAETZ_SUPERBEE`, `GRAETZ_JOBS` (4),
`GRAETZ_TRANSIENT_END` (1.5 s).  Needs python3 with numpy, scipy and
mpmath.

`tests/Alltest M7long` (not in `all`) runs `tests/verification/pipeMeshes`
(the h/2 first-layer mesh and the axial x2 mesh of the test mesh) and
`tests/verification/longStudies`, the long runs, managed as those of
M4long (keyed, kept with their binaries, reused, resumed; `M7LONG_RUNS`
selects runs, `ended` only evaluates; `M7LONG_REVALIDATE=yes` revalidates
a run of other sources by a 5-step probe; a run of another key that is
not selected is kept; `M7LONG_JOBS` x `M7LONG_RANKS`, default 2 x 8;
`M7LONG_CPUS` binds them, below):

    tests/Alltest M7long
    M7LONG_RUNS=ended tests/verification/longStudies   # evaluate only
    tests/verification/longStudies keys                # the key of every run

| study | runs | changed from the case |
|---|---|---|
| dt, 007 carrier (acceptance 3) | temperature-dt1ms/2ms/4ms, HKS-dt1ms/2ms/4ms | run/012 and run/013 (final tolerance 1e-14) at dt 1, 2 and 4 ms to 16 s, writes every second, a profiles dictionary (native and 1 cm bins, times every 0.5 s) |
| h vs h/2 (acceptance 2) | h-flux-none, h2-flux-none, h-A-none, h2-A-none, h-flux-halfCell, h2-flux-halfCell | run/014 at 106 mL/min, layer faceCells; the flux form (`depositionVelocity 0.0124421` m/s, the v_d that A gives on the base mesh) or the volumetric form (A 220 1/s); none or halfCell; the base or the h/2 mesh; to 30 s |
| axial x2 | x2-flux-none, x2-flux-halfCell | as the flux-form runs, on the axial x2 mesh (a bin set `cells` of its 848 axial cells besides the 424 native ones) |
| dt, paper-mode carrier (acceptance 3) | dt4ms, dt2ms, dt1ms | run/014 as distributed (100.2 mL/min, layer distance, A) at dt 4, 2 and 1 ms to 30 s |

Measured on v2412 (plan section 33; the runs on the paper-mode carrier are
the second M7 study's, revalidated by their probes):

- **Acceptance 3 on the 007 carrier** (M7.3a-d), at 16 s, the L1 change
  of the wall deposit (limit 1 % between 2 and 1 ms):

  | model | bins | L1(4 ms, 2 ms) | L1(2 ms, 1 ms) | order | Richardson error at 1 ms | peak height 2 -> 1 ms | T_dep 2 -> 1 ms |
  |---|---|---|---|---|---|---|---|
  | temperature (run/012) | native | 5.1e-7 | 2.7e-7 | 0.90 | 3.2e-7 | -2.4e-7 | +3.5e-7 K |
  | temperature (run/012) | 1 cm | 4.4e-7 | 2.2e-7 | 0.99 | 2.2e-7 | -4.9e-8 | +1.5e-6 K |
  | HKS (run/013) | native | 8.9e-4 | 4.5e-4 | 0.98 | 4.6e-4 | -2.2e-5 | +8.5e-4 K |
  | HKS (run/013) | 1 cm | 8.7e-4 | 4.4e-4 | 0.98 | 4.5e-4 | -6.0e-5 | +1.7e-3 K |

  The peak stays in its bin.  The gas of the temperature model has gone
  by 16 s (2e-14 to 4e-14 of the release left), but its final deposit
  still depends on dt at first order, weakly: L1(2 ms, 1 ms) 2.7e-7 (order
  0.90), as on the paper-mode carrier below, 1600 times below HKS (plan
  section 35, item Z; an earlier version said that it does not depend on
  dt); HKS still holds 3.6 % of its inventory in the gas at 16 s, and its
  deposit carries the first-order error of the transient.  While gas is present the error of
  the temperature model is first order (reported): 7.4e-3 at 4 s (3.2 %
  of the release deposited), 1.3e-5 at 6 s (the end of the release),
  6.2e-7 at 7 s; the gas at 7 s 5.4e-4.  solverDefect at the end: 5e-13 to
  4e-12 of the release (temperature), 1.4e-10 to 5.4e-10 of the inventory
  (HKS with the final tolerance 1e-14; 2.2e-8 to 6.0e-8 with 1e-12).

- **Acceptance 2** (M7.2c-f), h vs h/2 at 30 s, the release deposited:

  | form | T_dep h -> h/2 (native bins) | onset shift | largest `mDep_` on WALL | L1 of the wall profiles (native, 1 cm) | verdict |
  |---|---|---|---|---|---|
  | flux, none (M7.2d) | 679.754 -> 679.758 K | -0.002 bins | -0.024 % | 1.1e-3, 1.1e-3 | met (limits 1 bin, 5 %) |
  | flux, halfCell (M7.2f) | 679.754 -> 679.759 K | -0.002 bins | +0.006 % | 1.2e-3, 1.2e-3 | met |
  | volumetric A, none (M7.2e) | 679.754 -> 679.676 K | +0.032 bins | -32.0 % | 0.47, 0.47 | misses the 5 % (as expected) |

  In the flux form the deposition velocity per wall area is the same on
  both meshes; what remains is the cell temperature in f(T_c)
  (`interfaceTemperature cell`): the first cell lies 0.097 K (h) and
  0.048 K (h/2) above the wall at the onset.  With A, v_d = A V_I/A_I
  halves on h/2: the peak falls by 32 % and moves 7 mm downstream, while
  the onset moves by 0.03 bins only.
- **Axial x2** (M7.5a-c; criterion: onset within 1 native bin, 1 cm peak
  within 5 %, L1 of the 1 cm profiles below 1 %): onset +0.51 K (-0.21
  bins), 1 cm peak +0.0007 %, L1 9.8e-4, alike for none and halfCell.  The
  onset moves because axial2 resolves where the cell temperature crosses
  Tdep inside a base cell; the profiles and the peak do not move.
- **Acceptance 3 on the paper-mode carrier** (M7.3pa-pc), run/014 as
  distributed, at 30 s: L1(2 ms, 1 ms) 1.29e-7 (native bins) and 7.14e-8
  (1 cm), observed order 0.94 and 1.08.  Once the gas has gone, the final
  deposit of the temperature model still depends on dt at first order,
  but weakly: L1(2 ms, 1 ms) about 1e-7 here and 2.7e-7 on the 007
  carrier at 16 s (order 0.90), against 4.5e-4 for HKS on the same
  carrier (plan section 35, item R; an earlier version said that it does
  not depend on dt).  While gas is present the error is first order
  (reported, plan section 28, item 7): 7.1e-3 at 5 s (the largest with
  more than 1 % deposited), 1.2e-4 at 20 s; the gas at 16 s 3.7e-4.
- Every run: closure at most 6.6e-16 of the reference in every step,
  solverDefect at most 2.6e-9 at the end (M8's gate is 1e-8).

On this workstation (2 sockets of 10 cores, 40 logical CPUs with
hyperthreading) bind the runs to logical CPUs of distinct physical cores,
e.g. `M7LONG_CPUS=0-15 tests/Alltest M7long` (the driver re-executes
itself under `taskset -c 0-15`, and the runs inherit the set): on
2026-10-03 two unbound 8-rank runs (`mpirun --bind-to none`, as every
test) next to four serial Graetz runs took about 1 s per step, with ranks
on both hyperthreads of a core; bound, 0.15 to 0.3 s.

A complete `tests/Alltest M0 M1` takes 6 to 10 minutes on a loaded
40-core node (4 serial runs at a time, at most 8 MPI ranks); M2 adds about
10 minutes (its M8-scale channel, about 4 minutes, runs in the background).
`tests/Alltest -fresh all` (M0 to M5 with M4) took 70 minutes on the
40-core workstation at a load of 60-70 (other users' jobs), of which M4
about 10 minutes; after the M4 review, round 1, `tests/Alltest -fresh all
M4long` (the long runs reused) took 50 minutes at a load of 30-40; after
round 2, `tests/Alltest -fresh all` took 48 minutes at a load of 30-43
(three 4-rank long runs alongside), M4 about 8 of them, and `tests/Alltest
M4long` with every run reused 3 minutes; after round 3, `tests/Alltest
-fresh all` took 30 minutes at a load of 15-16 (one 4-rank long run
alongside), M4 about 8 of them; after the rebase onto the main branch
(plan section 31), 26 minutes at a load of 4-7 (231 passed, 1 not run:
M0.4); after M9 (plan section 32), 47 minutes at a load of 17-19 (258
passed); after the M7 integration (plan section 33), `tests/Alltest -fresh
all M7long` with every long run and the Graetz runs reused took 29 minutes
at a load of 1-17, bound to the 20 physical cores (`taskset -c 0-19`): 299
passed (283 of `all`, 16 of M7long), 1 not run (M0.4).  After the cleanup
of plan section 35, `tests/Alltest -fresh all` (with `LESTO_OPENFOAM_V2606`,
bound with `taskset -c 0-19`) took 41 minutes at a load of 2 to 11 (an
8-rank long run of M7long alongside for 18 of them): 300 passed, 1 not run
(M0.4); of it M7 (the Graetz runs made again, since the solver changed)
10 minutes and M9 2 minutes (M9.3b: 42.9 us per warm call).  Then
`M4LONG_REVALIDATE=yes M4LONG_RUNS=ended M7LONG_REVALIDATE=yes
M7LONG_RUNS=ended tests/Alltest M9 M4long M7long` revalidated the 29
recorded long runs that the change of the sources had left under other
keys and evaluated all of them in 10 minutes (91 passed), after the
cheapest run of each family, T5_60-distance and temperature-dt4ms, had
been made again with the final sources (plan section 35, item S: about
20 and 17 minutes on 8 bound cores; both reproduce their records bit for
bit).  After the follow-up of that cleanup (plan section 35, items U to
Z), on the workstation loaded by other users' jobs (a 1-minute load of 18
to 38), `tests/Alltest -fresh all` (with `LESTO_OPENFOAM_V2606`) took 49
minutes: 303 passed, 1 not run (M0.4); of it M9 about 3 minutes (M9.3b:
64.2 us per warm call at a load of 18.8) and M7 2, its Graetz runs made
just before with the final sources (20 minutes on 4 cores, since the
solver changed).  T5_60-distance and temperature-dt4ms, made again with
the final sources, took 38 and 27 minutes on 8 bound cores, one after the
other, and reproduce their records bit for bit; then
`M4LONG_REVALIDATE=yes M4LONG_RUNS=ended M7LONG_REVALIDATE=yes
M7LONG_RUNS=ended tests/Alltest M4long M7long` took 9 minutes (49 passed,
the other 29 runs revalidated).
`tests/Alltest M4long` takes hours when its runs have to be made (plan
section 25, item 7); with `M4LONG_REVALIDATE=yes` a run made with older
sources is revalidated in a few minutes instead, for a change known to act
alike in every step (plan section 26, item 11, and section 27).  So does
`tests/Alltest M7long`: its runs on the 007 carrier took 1 h 37 min on 2 x
8 bound ranks (the 1, 2 and 4 ms runs of both models side by side: 47, 30
and 20 minutes, 0.18 to 0.29 s per step; 26 core-hours), those on the
paper-mode carrier 4 h 22 min (the second M7 study); their revalidation
took 2 minutes.

## M10a: independent species pairs and element accounting

M10a adds simultaneous PbI2/BiI3 transport with different formulas and
kinetics, prescribed multi-species release, and separate inventory
reservoirs in a shared boat region. Each pair remains independent; gas
re-speciation, reactive deposition and shared condensates belong to later
M10 milestones. `run/030-pbi2-bii3-channel` is the first qualification case.

Add formulas in `constant/speciesTransportProperties`:

```foam
species (PbI2_g PbI2_s BiI3_g BiI3_s);
PbI2_g { state gas; diffusivityModel PbI2He;
         molarMass 0.46100894; formula { Pb 1; I 2; } }
PbI2_s { state solid; molarMass 0.46100894; formula { Pb 1; I 2; } }
BiI3_g { state gas; diffusivityModel constant; D 1e-4;
         molarMass 0.58969381; formula { Bi 1; I 3; } }
BiI3_s { state solid; molarMass 0.58969381; formula { Bi 1; I 3; } }
```

Formulas are optional for legacy cases. When used, every pair needs a
matching gas and condensate formula; each species declaring a formula
must declare `molarMass`. Positive integer atom counts and the calculated
mass are checked, with relative tolerance 1e-6. The explicit atomic-mass
registry currently covers Pb, Bi, I, He, Ar, N, O, H and K. Formula records
are written to `uniform/phaseChangeLayout` and changes are refused on
restart. Element accounting covers the configured gas/condensate pairs.

An optional per-pair `HKS` dictionary overrides `HKSCoeffs`:

```foam
pairs {
  BiI3_g {
    condensed BiI3_s;
    HKS { accommodation 0.8; Ce 1; kineticScale 1; }
    vapourPressure {
      file "$LESTO_M10_DATA/pv_BiI3_phases.csv";
      phases ("BiI3(cr)" "BiI3(l)");
      units log10Pa; interpolation logInverseT; outOfRange fatal;
    }
  }
  // PbI2_g has its own condensed/vapourPressure entries as before.
}
```

The accommodation must be finite and in (0, 1]; Ce and kineticScale must
be finite and nonnegative. Inventory samples can override their Ce and
kineticScale. Without per-pair `HKS`, legacy inventory defaults remain 1,
while wall kinetics keep using `HKSCoeffs`.

A shared source schedule uses species amounts in mol:

```foam
samples {
  boat {
    amounts { PbI2_g 4.84e-8; BiI3_g 6.17e-8; }
    selection {
      boat { action use; source box; box (.01 .002 -1) (.02 .0028 1); }
    }
    mode release; startTime 0; duration 2;
  }
}
```

Use `mode inventory; areaPerVolume 4000;` instead of the release schedule
for HKS reservoirs; an optional `removeTime` removes all listed reservoirs
on the same clock. Legacy `pair`/`amount` still works. Do not mix the two
syntaxes. Internally a multi-species source has names
`<sample>__<gasSpecies>`, which appear in restart records. These names must
not collide with explicit sample names. Different pairs may share cells;
samples of the same pair must be disjoint, and all samples remain disjoint
from the wall interaction layer. Shared cells have one exchange-cell
entry containing each pair's own reservoir element. Pair profiles exclude
only the gas in that pair's own sample cells from `wallPlusGas`.

`diffusivityModel ChapmanEnskog` is available for dilute binary transport.
Replace `D` with explicitly sourced Lennard-Jones parameters:

```foam
BiI3_g {
  state gas; molarMass 0.58969381; formula { Bi 1; I 3; }
  diffusivityModel ChapmanEnskog;
  ChapmanEnskogCoeffs {
    sigma 5; epsilonOverK 300;             // EXAMPLE estimates, not validated BiI3 data
    carrierSigma 2.576; carrierEpsilonOverK 10.22;
    carrierMolarMass 0.004002602;
  }
}
```

Sigma is in angstrom, epsilon/k in K, molar masses in kg/mol. The frozen
T and absolute p are in K and Pa. Parameters and states must be finite and
positive. Lorentz–Berthelot mixing and the Neufeld diffusion collision
integral produce D in m²/s once at startup; `pressureUnit` does not change
this model's SI inputs. The existing PbI2He correlation is unchanged. The
example parameters illustrate syntax and must be replaced by documented
estimates for research runs. The implementation references
[Neufeld et al. (1972)](https://doi.org/10.1063/1.1678363) and the SI
expression in [Langenberg et al. (2020)](https://doi.org/10.5194/acp-20-3669-2020).

With formulas and balance `writeFile yes`, element ledgers are written to
`postProcessing/phaseChangeElements/<startTime>/<element>.dat`, in mol of
atoms. Every entry is a compensated stoichiometric `nu/M` sum of the
existing pair ledgers; no second accumulated restart state is introduced.
The files include gas, wall, sample, supplied, closure, the independently
summed pairClosure and reference. Transport is positive outward, as in
the pair ledger. Restart files begin a new output segment.

With profiles enabled, `element_<element>_<binSet>.dat` uses the existing
conservative binning, with gas/wall/sample and observable in mol of
atoms/m; `molarMass 1` is the unit conversion for these already elemental
amounts. Formula-enabled pair and element profiles also record `Tpeak`,
the maximum bin of the selected observable and its wall temperature.
Insignificant peaks or bins without wall temperature report none. `Tpeak`
is separate from the existing 1% onset `Tdep`; peak locations depend on
bin size. Legacy cases without formulas keep their original output format.

Run `tests/Alltest M10a` (or `all`). The code gates use synthetic Bi tables.
Set `LESTO_M10_NASA`, `LESTO_M10_THERMOFUN` and `LESTO_M10_DATA` for the
external-data gates; otherwise those gates report NOTRUN. Generator usage
and redistribution boundaries are documented in
`../../thermochemistry/systems/README.md`.

Final M10a verification on 2026-10-06: `tests/Alltest all` on v2412,
with external data, reports 310 passed, 0 failed, 1 not run. All seven
M10a gates pass. M4.6 builds every source against local v2606 with zero
warnings; M0.4 is not run because this invocation uses v2412. Original
regression logs/fields, run/012–015 excerpts, M7 and M9 gates pass.
The separate Liu benchmark preparation suite passes all 6 tests.
The 31 completed long studies also pass five-step old/new binary probes:
35 additional checks pass, none fail. Their fresh whole-trajectory closure
gate M4.17 remains NOTRUN; the historical trajectories are revalidated,
not rerun in full.

## M10b: homogeneous gas re-speciation

`respeciation` is optional. With no entry, or `engine none`, the existing
M0–M10a paths keep their arithmetic and outputs. `engine kernel` solves the
ideal-gas equilibrium of Pb, Bi and I after all gas transport equations,
once per scheduled physical step. Condensed phases are suppressed in this
calculation; the bounded HKS exchange continues to handle deposition.

```foam
respeciation
{
    engine kernel;
    species (Pb_g PbI_g PbI2_g Bi_g Bi2_g BiI_g BiI3_g I_g I2_g);
    formationConstants "<constant>/log10Kf_PbBiI.csv";
    minTemperature 350;
    updateInterval 1;
}
```

Every participating gas needs its integer `formula` and checked
`molarMass`. Include every paired gas and the monatomic gas of every
represented element. The kernel supports gases containing Pb or Bi with
iodine, plus iodine species; mixed Pb-Bi molecules are rejected. Formation
constants use monatomic gases and the 1 bar standard state. The generator
is `thermochemistry/systems/make_kf_tables.py`. Standard-state constants
are cached for the frozen cell temperatures. Below `minTemperature`, gas
values remain those of the transport solve. Exactly absent elements stay
absent; an invalid negative element total leaves that cell unchanged and
is counted. Failed equilibrium or a balance residue above 1e-14 stops the
run. A slightly negative species can recover through equilibrium when its
element totals remain valid; no species is clipped.

Gas-only participants use the same unrelaxed, conservative Euler
transport equation as the paired path with zero exchange. Their final
solve needs `relTol 0`; optional `<field>Final` entries are supported.
The flux-form matrix defect is booked independently of the chemistry.
They can enter through initial fields or boundaries; sample release still
uses the existing condensate-pair sample interface. Combining re-speciation
with `HKSCoeffs speciation lagged` is refused. `engine GEMS` is implemented by the M10c backend described below.

The existing pair ledger retains its binary layout and balance-file
columns. Its supplied amount includes the reaction amount from a separate
species account when the new engine is active. New binary global lists
are written with every field write:

- `uniform/phaseChangeSpecies`: initial gas, release, boundary transport,
  exchange, solver defect, restart, REACTION, COPRODUCT, held gas and
  cumulative absolute reaction variation, with compensated sums. COPRODUCT
  is zero until reactive channels are added. Species references include
  reaction variation so a small net reaction does not erase a large
  cumulative transfer history.
- `uniform/phaseChangeElements`: compensated CONVERSION amounts in mol of
  Pb, Bi and I atoms. Each is the measured element residue of the stored
  before/after gas masses. Element supplied amounts use this rounding
  transfer, while the element reaction sum is also printed separately.
- `uniform/phaseChangeCondensates`: per paired condensate, initial wall and
  sample amounts, current wall and sample amounts, exchange, removal,
  clamping and restart. Entries follow the participating pair order.

The layout records the engine, ordered species, formulas, molar masses, temperature
threshold, update interval, ledger version and formation-table content
hash. A restart with changed settings or missing new accounts stops.
The kernel starts cold on every call, so no unsaved potential state changes
restart results. Same-decomposition binary restarts are byte-identical in
serial and on four ranks. Gas-only species contribute to the element
profiles, including the iodine observable. Each update reports calls,
iterations, invalid cells, the largest residue and mole-fraction/negative
mass monitors for every participating gas.

`tests/Alltest M10b` adds independent Decimal references, two initial
species splits in closed boxes at 450/1000 K, disabled-mode comparisons,
serial/MPI closure and restart, input refusals, written-field equilibrium
checks, and optional external thermodynamic checks. Written-field and
idempotence differences are normalized by the largest element amount;
reference species comparisons use relative errors for species above
1e-10 of that amount. The external exact GEMS oracle is polished to a
2e-15 residual target before comparing all 160740 stored gas-map inputs;
the original 1e-13 oracle stopping tolerance was insufficient for some
minor iodine species. No external thermodynamic records are bundled.

The reproducible preliminary case is
[`run/031-pbbi-gas-speciation`](../../run/031-pbbi-gas-speciation/readme.md).
Its common diffusivity and carrier ramp are provisional. In the 1173 K
prescribed-source study, the 1 ms BiI3 deposit differs from M10a by 0.1332%
in 1 cm bins, and the 2 ms to 1 ms change is 0.0575%; both are below the
1% gates. PbI2 differences are smaller. Experimental validation still
needs the physical thermal/carrier inputs and source constraints.

On 2026-10-06, `tests/Alltest all` completed with 319 PASS, 0 FAIL and
1 NOTRUN. All nine M10b gates passed with the external inputs enabled,
including the complete 160740-point oracle comparison and the 4/2/1 ms
channel study. Representative cold kernel calls took 16.2 us. Both the
v2412 build and the separate v2606 Debug compile had zero warnings.
M0.4 was NOTRUN because this invocation used v2412; M4.6 tested the
v2606 compilation. The hours-long M4long/M7long studies are separate
invocations.


## M10c: GEMS standard states and guarded local gas calls

The optional bridge now supports `respeciation { engine GEMS; }`. Configure
its backend in `respeciation/GEMS`; `GEMSCoeffs` can supply shared defaults:

```foam
respeciation
{
    engine GEMS;
    species (Pb_g PbI_g PbI2_g Bi_g Bi2_g BiI_g BiI3_g I_g I2_g);
    minTemperature 350;
    updateInterval 1;
    GEMS
    {
        mode frozen; // or local
        system "<constant>/gems/PbBiIHe-dat.lst";
        carrier "He(g)";
        minMoleFraction 1e-12;
        maxMoleFraction 3.2e-4;
        warmIterationLimit 200;
        logDirectory gasGemsLog;
        logLevel 4;
    }
}
```

`frozen` caches GEMS standard-state formation constants at each distinct
frozen cell T/P state, per rank, and uses the exact gas kernel every step.
`local` also makes per-cell warm IPM calls where the temperature/pressure
grid and total reactive-gas mole-fraction guards allow them. Valid duals
start a kernel refinement against the physical, unfloored element amounts.
A final cold refinement makes the fields independent of unsaved warm state,
including nearly exhausted iodine. Rejected or failed calls use the same
kernel and cached G0. Fatal calls recreate the engine and invalidate all
warm states of the rank. Counts, guards, failures and the combined IPM plus
refinement cost are reported every update.

Gas names map automatically from `Pb_g` to `Pb(g)`; an optional `gases`
dictionary maps custom solver names. Every reactive gas of the GEMS system
must participate. Formula, molar mass, carrier, guard and control checks
stop invalid configurations. Unknown keys are refused. Below the temperature
threshold, transport values are retained. `formationConstants` is optional
for this engine: provide an open-data table to allow a kernel fallback at
frozen T/P states outside the GEMS grid; without it those states stop startup.

The existing M10b species/condensate/element account layouts are retained.
Restart metadata also fingerprints the backend mode, guards, cap, mapping,
and every GEMS system document. Changed backend settings/data are refused.
Frozen-mode binary restarts are exact in serial and on four ranks; local
mode is qualified to the required 1e-12 inventory tolerance and has exact
fields after the canonical refinement.

For HKS pairs using `equilibrium GEMS` in a system containing Bi, equilibrium
pressure comes from dual fugacities. `GEMSCoeffs/maxMoleFraction` (default
3.2e-4) joins the guards, warm calls are capped at 200 iterations before
cold retry, and `speciation lagged` requires the pair's gas mole fraction
at least 1e-7. Upper-bound and minor-gas refusals use the existing REJECTED
source code and table fallback. Original Pb-I HKS behavior is retained.
Re-speciation and `speciation lagged` remain mutually exclusive.

`tests/Alltest M10c` includes the original ABI comparison, open Pb-I and
external Pb-Bi-I serial/MPI channels, restarts, invalid inputs, injected
fatal/NaN-dual faults, the HKS guards, a 936-point reduced map and the
prescribed-source deposit comparison. The map records raw IPM errors by
species mole-fraction decade. A general 70-digit Decimal Newton solve
refines the external oracle to the exact fixed-volume inputs: its double
KKT residue alone is insufficient for some nearly exhausted minor iodine
species. The gate compares both cold and warm kernel refinements with
that independent reference, with a 1e-9 relative tolerance for species
above 1e-10 of the largest element amount.

The reproducible case is
[`run/032-pbbi-gems-speciation`](../../run/032-pbbi-gems-speciation/readme.md).
Its thermal ramp, diffusivity and prescribed source remain provisional.
External inputs are selected using the M10 environment variables and
`LESTO_M10_SYSTEM`. `LESTO_M10C_SCIENCE_RECHECK` explicitly rechecks an
already completed three-mode channel study; the default test invocation
runs the study afresh. The criterion identifies such a recheck.


## M10d: shared wall condensates and Q2 sources

With `model HKS` and an enabled re-speciation engine, opt into `wallChannels`
for metal channels sharing a physical deposit. Declare the gas/solid fields,
molar masses and formulas in `speciesTransportProperties`, and select every
reactant, gas product and source gas in `respeciation/species`.

```foam
wallChannels
{
    Bi {reactant Bi_g; condensed Bi_s;
        vapourPressure {file "<constant>/pv_Bi_phases.csv";
            phases ("Bi(cr)" "Bi(l)"); units log10Pa;}}
    Bi2 {reactant Bi2_g; condensed Bi_s; stoichiometry 2;
        vapourPressure {file "<constant>/pv_Bi2_phases.csv";
            phases ("Bi(cr)" "Bi(l)"); units log10Pa;}}
    BiI {reactant BiI_g; condensed Bi_s;
        stoichiometry {BiI_g -3; Bi_s 2; BiI3_g 1;}
        equilibrium {file "<constant>/logK_BiI_disp.csv"; column logK;}
        reversible no;}
}
gasSources
{
    boat {amounts {Pb_g 2e-12; Bi_g 3e-12; I_g 4e-12;} // mol of gas
        box ((.01 0 -1) (.02 .0048 1)); startTime 0; duration .006;}
}
```

The source example is synthetic; run/033 provides the provisional LBE Q2
source. Source amounts are distributed uniformly over the selected volume,
with exact overlap of the release window and each time step. Paired gases
use the existing sample release mechanism instead.

Plain channels are reversible by default. A scalar stoichiometry gives the
number of condensate formula units per gas molecule; the dimer requires
two. A dictionary describes an atom-balanced reaction, with one consumed
gas, one positive condensate and optional positive gas products. Reactive
channels are irreversible. Their equilibrium table contains dimensionless
log10(K), with gas standard pressure 1 bar. The lagged product pressure gives
`p*_BiI = 1 bar * ((p_BiI3 / 1 bar) / K)^(1/3)`. Gas products are added only
after all gas transport solves, then the homogeneous re-speciation runs.
The pressure is lagged at every time step, regardless of the homogeneous
update interval.

Each channel's optional `HKS` block accepts `accommodation`, `Ce` and
`kineticScale`, independently defaulting to one. The interface geometry,
wall temperature and `wallResistance halfCell` are shared with the existing
pairs. Channel tables use log/inverse-T interpolation by default and refuse
temperatures outside the table; `interpolation` and `outOfRange` can be set
in the table block. Reaction tables use `column`, without `units` or `phases`.

Channels run in the species transport order. After each gas solve, the
shared deposit changes before the next gas obtains its evaporation bound.
Duplicate reactant/condensate channels and stores also owned by a legacy
pair are refused. The same gas may deposit into different stores, allowing
an optional plain BiI_g/BiI_s channel alongside disproportionation. The
predictor and at most eight guard re-solves bound each wall event; CAPPED
empties the store exactly, with rounding residues recorded separately.

`mw_<solid>`, `mDep_<solid>` and derived `Y_<solid>` describe each physical
store. Its profile is `condensate_<solid>`; element profiles count it once.
`phaseChangeChannelCondensates` balances initial mass, exchange, clamping,
restart adjustment and held mass. Its reference includes cumulative absolute
exchange, so a nearly empty store is compared against the mass processed.
Species references also include shared-store exchange and gas products.

The binary `uniform/phaseChangeWallChannels` holds exact per-rank stores;
`uniform/phaseChangeChannelCondensates` holds global compensated accounts.
Same-decomposition restarts reproduce fields and accounts byte for byte.
After re-decomposition, `mw_` reconstructs the stores and any inventory
rounding adjustment is booked. Fingerprints refuse changes to channel
definitions, source schedules or table contents across restart.

[`run/033-lbe-metal-vapour`](../../run/033-lbe-metal-vapour/readme.md)
contains the reproducible channel. `tests/Alltest M10d` includes analytical
boxes, optional BiI(cr), half-cell resistance, source/ledger/profile checks,
serial/four-rank and changed-decomposition restarts, and input refusals.
With `LESTO_M10D_NASA`, `LESTO_M10D_THERMOFUN` and `LESTO_M10D_REFERENCE`,
it also compares both BiI variants and two time steps with fractional plug
flow. External thermodynamic records remain outside the repository.

[`run/034-liu2025`](../../run/034-liu2025/readme.md) prepares M10e's physical
helium wedge carriers and Q1/Q2 column matrix from digitised paper curves.
Fifteen 60 s coarse cases and five full-size smoke cases pass numerical
conservation checks. Saved peak/profile comparisons and selected time-step
checks remain preliminary: a temperature profile is missing, Figure 7C
contains a curve/label discrepancy, and transport, source and thermodynamic
inputs still need qualification. The full-size 60 s matrix is prepared but
unrun. Reactive wall equilibrium handles zero/subnormal product pressures
in log space, with an FPE-enabled regression test.
