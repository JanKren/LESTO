# Case 014 — code-to-code case of the temperature-based model (paper mode)

Case of milestone M4 of
`applications/rhoFixedFlowFoam/doc/phase-change-plan.md`: the
temperature-based deposition model of Lobresco et al. (2026), Nuclear
Engineering and Design 459, 115201 (Section 2.3.1, Eq. 5), in the paper
mode of the T-Flows reference code: a frozen, non-expanding parabolic
helium flow of density 1, the temperature of that flow with the measured
wall table, and the release of the sample over the injection interval.
Every value of the T-Flows parity checklist (plan section 11) is pinned
below with its source; plan sections 25 and 26 record the decisions.

The case as distributed has the configuration of Figs. 8-10 as the paper
states it (Section 3.1 and the Fig. 8 caption: A = 220 1/s, k = 8e-3 1/K,
Tdep = 680 K, 14.6 mg = 3.16e-5 mol, 100 mL/min, t_ev = 15 s), run to
2 t_ev = 30 s (Fig. 8 shows 16 s, Fig. 9 t_ev + 1 = 16 s and 2 t_ev = 30 s,
Fig. 10 30 s).  **It does not reproduce the positions of those figures**:
its deposit begins 2.6 cm upstream of the published one (section "Figs.
8-10" below; plan section 26).  The other runs of acceptance 3-5 are edits
of it (table at the end); `tests/Alltest M4long` makes and runs them.

## Sources

- **Paper**: Lobresco et al. (2026), `papers/Lobresco2026.pdf`; cited as
  section, equation, figure or table.
- **T-Flows**: `/home/jan/src/T-Flows` (read-only), cited as
  `branch:file:line`, with
  - `HKS` = `origin/advanced_scalar_model` (35c1f6406e, 2026-01-25),
    `Tests/Laminar/Scalar_Transport_Gpu/LESTO_HKS/` — the paper-mode
    carrier, the layer and the sample of the HKS runs;
  - `Thesis` = the same branch, `Tests/Laminar/Scalar_Transport_Gpu/
    Lesto_Thesis/` — a temperature-model case with other parameters;
  - `Project` = `origin/multiple_scalars_gpu_francesco` (fc5065643c),
    `Tests/Laminar/LESTO-project/` — likewise;
  - `Transport` = `origin/scalar_transport_francesco` (2f574f9c40),
    `Tests/Laminar/Scalar_Transport_Gpu/LESTO/` — a scalar-transport test,
    not a temperature-model case: no `User_Mod` (so no deposition source),
    the scalar enters through the inlet as a flux (`q_01 1.0`,
    `Transport:control:88-91`), 12 time steps (`:13`); its control file,
    geometry (`LESTO_pipe.geo`) and wall table remain sources for the
    advection keyword, the tolerances and the domain;
  - `Gpu` = the GPU solver of these branches, `Sources/Process/Gpu/` (the
    LESTO cases are GPU cases: `Main_Pro_Frozen.f90` uses `Gpu_Mod`).

  No branch holds the temperature-model case of the paper's figures with
  the paper's parameters (`manuscript/notes-for-authors.md`, items A, C,
  D): the carrier, the layer and the sample follow `HKS`, the model and its
  parameters follow the paper.  The items of the checklist that belong to
  the HKS model only (the condensation coefficient sigma, `kineticScale`,
  the vapour-pressure table, the sample's 2/R_s x 5 area per volume) do
  not apply to this case; they are pinned with the HKS code-to-code case
  of milestone M6.

## Pins (plan section 11)

**Carrier** (`system/setFrozenCarrierDict`, written into `0` by
`setFrozenCarrier`, `applications/rhoFixedFlowFoam/utilities/`)

| Item | run/014 | Source |
|---|---|---|
| density | `density 1` kg/m3, uniform and constant | `HKS:User_Mod/Beginning_Of_Iteration.fpp:41` (`Flow % density(c) = 1.0`) |
| velocity | `U = 2 U_b (1 - r^2/R^2) e_x`, uniform along x (non-expanding), `U_b = Q/(pi R^2)`, R = 2.4 mm | `HKS:User_Mod/Initialize_Variables.f90:23-36` (`u_max = 2 FLOW_RATE/(pi R_PIPE^2)`), `HKS:User_Mod/Types.f90:5` (`R_pipe = 2.40e-3`) |
| flow rate, Figs. 8-10 | `flowRate 1.67e-6` m3/s (100.2 mL/min): U_b = 0.0922878 m/s.  The temperature-model codes prescribe instead `U_MAX = 0.176` m/s, used as the centreline velocity (`Thesis:User_Mod/Initialize_Variables.f90:25`, `u = U_MAX (1 - r^2/R^2)`; `Types.f90:6` calls it a bulk velocity): U_b = 0.088 m/s, Q = 95.5 mL/min, whose plateau n/(U_b t_ev) = 2.40e-5 mol/m at t_ev 15 s Fig. 9 (2.26e-5) rules out | `HKS:User_Mod/Types.f90:7` (`FLOW_RATE = 1.67e-6`); paper Section 3.1 (100 mL/min); plan section 27, item 14 |
| flow rate, Table 2 | `flowRateMLPerMin` 106, 5 or 540, used as a volumetric rate of the non-expanding parabola | paper Table 1 |
| facet-area rescale | off (`rescaleToFacets no`): the discrete inflow is 0.99997 Q, the faceted inlet area 0.99493 pi R^2, the discrete bulk velocity 0.0927551 m/s (all printed).  T-Flows' own discrete flow is 1.0038 Q on this mesh (its fluxes evaluate the parabola at the vertex-mean centres of its cells and faces; U_b,d = 0.0931139 m/s), so at the same FLOW_RATE this carrier is 0.39 % slower than T-Flows'; Fig. 9's plateau is that of T-Flows' carrier (row "gas plateau" below).  A carrier with T-Flows' total flow on this mesh: `flowRate 1.67641e-6` with `rescaleToFacets yes` (for the profile comparisons of M6; not used here) | T-Flows uses the analytic R in `u_max` (`Initialize_Variables.f90:24`), `v_flux = sx 0.5 (u(c1) + u(c2))` with u at the cell centres (`:41-46`) and at the boundary cells of the inflow faces (`:51-59`), centres that are node averages (`Sources/Shared/Grid_Mod/Calculate/Cell_Centers.f90`, `Face_Centers.f90`); plan section 27, item 9 |
| face fluxes | the exact integral of the analytic velocity over every face (`max |sum_f phi|` per cell 2e-18 of the inflow); the wall faces carry +-1e-25 kg/s (8e-16 of the inflow in total), the exact flux through faces that gmsh tilted by 1.5e-15 (printed as the wall flux) | T-Flows interpolates cell velocities to the faces (`Initialize_Variables.f90:41-46`); plan sections 25 and 26 |
| pressure | `pressure 1e5` Pa, uniform: the diffusivity at 1 bar | `HKS:User_Mod/Beginning_Of_Iteration.fpp:22` (`p_ref = 1.0` bar) |
| temperature | `mode solved`: the steady temperature of the frozen flow, with the wall table on WALL (Dirichlet at the face centres), 900 K at the inlet, zero gradient at the outlet, `rho cp = 5193` J/(m3 K), k(T) of helium | T-Flows solves the energy equation in every iteration (`HKS:Main_Pro_Frozen.f90:190`) with the table as the wall condition (`HKS:control:80-83`), inlet `t 9.0e2` (`HKS:control:85-88`), `HEAT_CAPACITY 5.193e3` (`HKS:control:25`), density 1 (above), k = `k_He` of `HKS:User_Mod/Types.f90` interpolated linearly (`HKS:User_Mod/Beginning_Of_Iteration.fpp:33-36,46-47`) |
| temperature, convection | `Gauss linearUpwind grad(T)`, not the paper's SuperBee: the steady Picard iteration does not converge with SuperBee (200 iterations, a change of 29 K left); upwind, linearUpwind and central differ by at most 0.023 K (mean 0.003 K) in the layer cells at x = 0.38-0.52 m | paper Section 2.2 ("Convective terms in the temperature and scalar-transport equations were discretized using the Superbee TVD scheme"); the T-Flows GPU solver runs its energy equation with upwind in `HKS`, `Thesis` and `Transport` and with central differencing in `Project` (row "convection" below); plan sections 26 and 27 |
| wall table | `constant/wall_x_profile_complete.dat`, a byte copy of `HKS:wall_x_profile_complete.dat` (criterion M4.12); column 5 (T), x in m although the header says cm; all 70 rows (the count line says 66, which T-Flows reads: x <= 0.65 m).  T-Flows assigns no temperature to the wall faces beyond the last row it reads: they keep the boundary value 0 of their allocation, so in `HKS` the wall beyond x = 0.65 m is a 0 K Dirichlet wall of its energy equation, and the conductivity of the cells it cools is looked up below the `k_He` table (extrapolated below 250 K, out of bounds below 230 K).  Neither `rows all` (the measured 315 to 305 K at 0.66-0.69 m) nor `rows count` with `outOfRange hold` (324 K) reproduces that; this case keeps `rows all`.  The difference matters beyond 0.65 m only, where T540_60 deposits 1.9 % (distance) and 7.1 % (faceCells) of its inventory by 8 s; its T_dep, at 0.48 m, is not affected | `HKS:control:83`; `Gpu:Read_Controls_Mod/Boundary_Conditions.f90:538-607` (a face in no interval of the table is not assigned) and `:670-683` (every boundary value copied to the field), `Gpu:Var_Mod/Create_Variable.f90:38` (`phi % b = 0.0`), `Gpu:Process_Mod/Insert_Energy_Bc.f90:49,68` and `Form_Energy_Matrix.fpp:189-199` (the whole WALL region Dirichlet), `HKS:User_Mod/Beginning_Of_Iteration.fpp:33-36` (`dataIndex` clamped above only); plan sections 25, 27 and 29 |
| wall table per case | the same table for every flow rate: the only measured table available.  The paper measured a dedicated wall profile for each operating condition (Section 2.4) and plots it in Figs. 13-14; for T5_60 it differs from this table (Fig. 14 iii: about 316 K at x = 0, 1050 K near 0.13-0.15 m, 793 K at 0.40 m; the table: 952 K, 1210 K at 0.04-0.05 m, 807 K at 0.40 m) | paper Section 2.4, Figs. 13-14; `manuscript/notes-for-authors.md`, item F; plan section 26 |

T-Flows integrates its energy equation in time from 500 K
(`HKS:control:72-74`); with density 1 its radial relaxation time
`R^2 rho cp/k` is 0.1-0.2 s, short compared with the 4-5 s the gas needs
to reach the deposition zone, so the steady temperature of `mode solved`
is the temperature the deposit sees in both codes.  The log line `Carrier:
max |T - p W/(R rho)|/T = 0.96` is expected: the paper mode's density 1 is
not the perfect-gas density of helium (0.04-0.17 kg/m3), and the solver
never uses the equation of state here.  The mole-fraction monitor measures
the gas against the carrier's p/(R T) (plan section 29, item 4): the gas
of this case is not dilute.  At the end of the release (15 s) x exceeds
0.05 in 148530 of the 228960 cells and 0.1 in 54783, at most 0.12 (the
plateau, 1.25 mol/m3, against 10 to 16 mol/m3 of helium at 1 bar), so the
log reports cells above `moleFractionWarning` and the mass that entered
the gas there; the model, like T-Flows, treats the gas as passive, an
approximation to about 10 % in this configuration.

**Species and transport** (`constant/speciesTransportProperties`,
`system/fvSchemes`, `system/fvSolution`, `system/controlDict`, `0/`)

| Item | run/014 | Source |
|---|---|---|
| species | `PbI2_g` (gas, M = 0.46100894 kg/mol) and its condensate `PbI2_s`; c = rho Y/M | paper Section 2.2; `HKS:User_Mod/Types.f90:17` (461.01 g/mol) |
| diffusivity | `diffusivityModel PbI2He`: the published Omega fit, `pressureUnit bar` with p = 1e5 Pa (p_ref = 1 bar, where Eq. 2 of the paper says atm); D = 3.17e-5 to 3.48e-4 m2/s | paper Eq. 2 and Fig. 3; `HKS:User_Mod/Beginning_Of_Iteration.fpp:22-27,61-72` |
| gas constant | R = 8.314 J/(mol K) in T-Flows; the temperature model uses no R (Eq. 5 has no partial pressure), so it matters only for HKS (M6) and for the printed carrier monitor (OpenFOAM's R) | `HKS:User_Mod/Types.f90:18`; plan section 12, item 9 |
| convection | `Gauss SuperBee`, as the paper states; diffusion central (`Gauss linear corrected`).  The T-Flows code does not use SuperBee: its GPU solver never reads `ADVECTION_SCHEME_FOR_SCALARS` (the keyword of `Thesis:control:45`, `Project:control:47` and `Transport:control:46` is ignored) and blends central, linear-upwind and upwind with `BLENDING_COEFFICIENTS_FOR_SCALARS` (and `..._FOR_ENERGY`), whose default depends on the branch: 0 0 1, pure upwind, on `advanced_scalar_model` and `scalar_transport_francesco` (`HKS`, `Thesis`, `Transport`), but 1 0 0, pure central differencing, on `multiple_scalars_gpu_francesco` (`Project`), for the scalars and the energy alike; no control sets it but `HKS`, 0.8 0.1 0.1 for the scalars (ramped in over the first 240 steps).  The plateau and T_dep do not depend on it; the shape of the profiles does (the numerical diffusion of upwind, U dx/2, about 1.4e-4 m2/s on the axis, is comparable to D, about 1e-4 m2/s at the onset): an upwind, a blended (0.8 0.1 0.1) and a central variant are planned for the profile comparisons of M6 | paper Section 2.2; `Gpu:Read_Controls_Mod/Numerical_Schemes.f90:101-105` ("Line reserved for advection scheme"), `Sources/Shared/Control_Mod/Numerics/Blending_Coefficients_For_Scalars.f90:19` and `Blending_Coefficients_For_Energy.f90:19` (`data def / 0.0, 0.0, 1.0 /` on advanced_scalar_model and scalar_transport_francesco, `data def / 1.0, 0.0, 0.0 /` on multiple_scalars_gpu_francesco), `Gpu:Field_Mod/Core/Add_Advection_Term.fpp:41-46` (advanced_scalar_model: central, linear-upwind and upwind weights) and `:95-100` (multiple_scalars_gpu_francesco: the first coefficient weights the central face value), `HKS:control:45`; plan sections 26 and 27 |
| time step | `deltaT 0.001` s (Courant number 0.11 at 100 mL/min, the case as distributed; 0.12 at the 106 mL/min of Table 2); implicit Euler | paper Section 2.2; `HKS:control:12` |
| inlet species BC | zero total flux (convective + diffusive) through INLET: `codedMixed`, refValue 0, valueFraction Pe/(1 + Pe) of every inflow face (Danckwerts) | `HKS:control:85-88` (`inflow`, `q_01 0.0`: a zero flux of the scalar) |
| outlet, wall | `zeroGradient` (the exchange enters the gas as a source) | `HKS:control:90-93` (`outflow`) |
| outer iterations | `nCorr 10`, `tolerance 1e-5`: the loop converges in every step but the first of the runs at dt 10, 5 and 4 ms (11 correctors, the limit, as in run/012), with loose GAMG correctors and a final PBiCGStab solve to 1e-12.  The correctors per step of the M4long runs (plan section 29, item 6): Fig. 9 at t_ev 6 s 2 in 89 % of the steps, at t_ev 15 s 1, 2 or 3 (37, 38 and 24 %); T106 2 in 67-78 %, 3 in 9-27 %; T540 (0.22 ms) 1 in 60 %, 2 in 35 %; T5_60 (10 ms) 3 or 4 (45-56 % and 28-41 %), at 5 ms 2 in 80 %; the T-Flows set (4 ms) 1 in 34 %, 3 in 30 %, 5 or more in 24 % | T-Flows: `MIN_SIMPLE_ITERATIONS 3`, `MAX_SIMPLE_ITERATIONS 6` (`HKS:control:14-15`; a second `MIN_SIMPLE_ITERATIONS 15` at `:66`), with `TOLERANCE_FOR_SIMPLE_ALGORITHM 1e-8` and a scalar solver tolerance of 1e-17 in `HKS` (`:58,65`) and 1e-3 for both in `Transport` (`:61,68`) and in the paper (Section 2.2).  The temperature-model source is implicit in T-Flows too (`Thesis:User_Mod/Source.f90:141-158`), so a converged step is the common target |
| domain | the pipe of run/001, x = 0-0.69 m, 424 axial cells (1.627 mm), a closed inlet at x = 0 and the sample at 0.04-0.05 m (the geometry of `HKS`) | `HKS`; the two temperature-model cases and the `Transport` test start the domain at the sample position (`Thesis:lesto_pipe.geo:16-17`, `Project:Stainless_Steel.geo:16-17`, `Transport:LESTO_pipe.geo:16-17`: `L1 = 0.04` "sample position", `L2 = 0.69`), with tables from x = 0.04 m and the release into the cells of the inflow faces (`Thesis:User_Mod/Source.fpp:61-73`).  Their mesh has 424 axial cells over 0.04-0.69 m, dx = 1.533 mm (`Thesis:lesto_pipe.geo:14-17,40,482`: `N_STREAM = L*400 + 1` = 425 nodes with L = 1.06; `Source.f90:28`, `inlet_volume = 2.7600231e-8` m3 = the faceted section 1.800384e-5 m2 times 0.65/424), so their native bins quantise T_dep and dT_bin differently.  In run/014 the gas can diffuse upstream of the sample, which the T-Flows temperature-model cases do not allow: at 5 mL/min (T5_60) 1.1e-6 mol/m of gas is left at x = 0.001 m at 60 s.  The rear flank of the gas at 16 s in Fig. 9 (whose x axis starts at 0.04 m) lies 3-7 mm upstream of run/014's (10/50/90 % at 0.099/0.129/0.160 m against 0.102/0.135/0.167 m; science review), consistent with the release positions: T-Flows at x = 0.04 m, run/014 over the boat, centre 0.045 m.  Plan sections 26 and 27 |

The diffusivity of the solver uses M_He = 4.0026 g/mol where T-Flows uses
4.0 (`HKS:User_Mod/Beginning_Of_Iteration.fpp:24`): D is 0.03 % lower;
the other constants agree to their printed digits.

**Layer and model** (`constant/thermochemistryProperties`)

| Item | run/014 | Source |
|---|---|---|
| model | `model temperature`: S = -A [1 - exp(-k (Tdep - T))]^+ phi in every layer cell (Eq. 5; zero above Tdep), irreversible | paper Eq. 5, Section 2.3.1 |
| A, k | `A 220` 1/s, `k 8e-3` 1/K | paper Section 3.1 (calibrated against Fig. 10, "adopted in all subsequent simulations") and the Fig. 8 caption; `Thesis:User_Mod/Types.f90:7-8` has 150 and 0.02 (the sensitivity variant below); `Project` lists the same 150 and 0.02 as the "Parameter used in the reference case" (`Project:reference.txt:4-5`), while its `User_Mod/Types.f90:6-7` holds 1.35 and 2.7 (both of commit 49806abdd8, "update with different reaction rate models") |
| Tdep | 680 K for Figs. 8-10, as the paper states; for Table 2 the case's T_dep^calc of Table 1 (775, 712, 666, 788, 671 K).  The published figures place the onset where this wall table holds about 640 K (section "Figs. 8-10"): both T-Flows temperature-model codes use 638.15 K | paper Section 3.1 (680 K from the semi-empirical formula), Section 2.3.1 ("evaluated using the semi-empirical formula given in Liu et al. (2025)"), Table 1; `Thesis:User_Mod/Source.fpp:25`, `Project:User_Mod/Source.f90:28` (`pb_t_dep = 638.15`) |
| layer | `layer distance`, `distance 1e-4`: the cells whose centre lies within 0.1 mm of WALL (the exact distance from the nearest wall face), 2 layers, 30528 cells, A_I/V_I = 8533.7 1/m; switch to `layer faceCells` (the cells of the wall faces, 1 layer, 15264 cells, 17681.9 1/m) | `HKS:User_Mod/Types.f90:14` (`DELTA = 1.00e-4`), `HKS:User_Mod/Initialize_Variables.f90:114-120` (`R_PIPE - r <= DELTA`, printed as the centre-distance rule: the same 30528 cells); paper Section 2.3.2 ("all control volumes lying entirely within" 0.1 mm: one layer).  The two temperature-model branches use r > 0.9 R (`Thesis:User_Mod/Source.f90:114`, `Project:User_Mod/Source.f90:111`): R - r < 0.24 mm, 4 layers, 61056 cells on this mesh (`layer distance; distance 2.4e-4` selects the same cells) |
| A_I/V_I | the faceted wall area over the layer volume (8533.7 1/m); T-Flows divides the analytic area 2 pi R L (8544.5 1/m) | `HKS:User_Mod/Initialize_Variables.f90:122` |
| rate per layer cell | `areaMultiplier 1`, `areaPerVolume` the default of the rule: `layer` with layer distance (by definition: v_d = A V_I/A_I = 0.0257802 m/s, the rate A f(T) in every layer cell), `cell` with layer faceCells (A_e = |S_f|: the rate per cell differs from A f(T) by |S_f|/V_c/(A_I/V_I) - 1, about 1e-6 on this mesh; the faceCells runs of M4long use it) | plan section 15, item 3; paper Eq. 5 is a volumetric rate |
| interface temperature | `interfaceTemperature cell` (the cell temperature) | `Thesis:User_Mod/Source.f90:115` (`flow_t_n(c)`) |

**Sample and amounts**

| Item | run/014 | Source |
|---|---|---|
| sample cells | `boat`: the cylinder x = 0.04-0.05 m, radius 4.8e-4 m: 144 cells (the rings of cell centres at r <= 0.458 mm; the next ring, at 0.5 mm, lies 5e-16 outside 0.5 mm) | `HKS:User_Mod/Types.f90:11-13`, `HKS:User_Mod/Initialize_Variables.f90:85-87`; plan section 22, item 5 |
| release | `mode release`: n0 released uniformly over t_ev from t = 0 into the sample cells | paper Section 2.3.1 (a constant release rate over the injection interval "at the inlet"); the temperature-model branches inject into the cells of the inflow faces at x = 0.04 m over half the run (`Thesis:User_Mod/Source.fpp:51-73`) |
| source position, Table 2 | x = 0.04-0.05 m for every case.  The arrows of Figs. 13-14 mark "the initial position of the PbI2 source term used in the simulations": about 0.13 m (T106) and 0.15 m (T5, T540) | paper Figs. 13-14; plan section 26 |
| n0, t_ev | 3.16e-5 mol (14.6 mg) and t_ev = 15 s for Figs. 8-10 (Fig. 8 caption; Fig. 9 shows t_ev = 15, 30, 45, 60, 75 and 90 s); Table 1's moles and t_ev = 6 s for Table 2.  Both T-Flows temperature-model codes release n = 3.17e-5 mol (14.6 mg/461 g/mol), not the printed 3.16e-5 | paper Section 3.1, Figs. 8-9, Section 3.3 ("an injection duration of 6 s was selected for all temperature-based simulations"), Table 1; `Thesis:User_Mod/Source.f90:27`, `Project:User_Mod/Source.f90:26` (`n = 3.17e-5`) |

## Outputs and T_dep

`profiles` (milestone M3): the native bins (the 424 axial cells) and the
1 cm gamma bins, at the write times (every second) and at `times (16 30)`:
t_ev + 1 (the gas of Fig. 9, top) and 2 t_ev (the solid, bottom).  T_dep
is the wall temperature of the first position from the inlet where the
condensed phase exceeds 1 % of its peak, interpolated linearly between
bins (paper Section 3.3): `onset { fraction 0.01; search (fromInlet ...);
}`, `observable wall` (the condensed phase, as Section 3.3 defines T_dep;
right at intermediate times, while `wallPlusGas`, the default, belongs to
the comparisons at shutdown).  `mDep_PbI2_s` maps the deposit of every
layer cell to its nearest wall face (plan sections 25 and 26); the profiles
bin a layer cell's deposit by its own cell.

## Figs. 8-10: what this case reproduces and what it does not

Measured on the native profiles of the M4long run `fig9-tev15` (this case
as distributed, t_ev = 15 s) against the paper, digitised by the science
review of M4 (criterion M4.16, not a gate):

| Quantity | paper | run/014 | difference |
|---|---|---|---|
| gas plateau at 16 s (Fig. 9) | 2.26e-5 mol/m | 2.2712e-5 mol/m | +0.5 % (acceptance 3: n0/U_b of the carrier) |
| gas at 50 % of the plateau, 16 s (Fig. 9) | 0.507 m | 0.4987 m | -8.3 mm |
| gas at 10 % of the plateau, 16 s (Fig. 9) | 0.518 m | 0.5212 m | +3.2 mm |
| solid onset (1 %), 30 s (Fig. 10; Fig. 9 starts at 0.500 m) | 0.502 m | 0.4761 m (T_dep 679.82 K) | -25.9 mm |
| solid peak, 30 s (Fig. 10) | 1.14e-3 mol/m at 0.513 m | 9.11e-4 mol/m at 0.4955 m | -17.5 mm |

The plateau (the identity n' t_ev = n0/U_b of acceptance 3) agrees, and
its +0.5 % is the slower carrier of this case: T-Flows' own carrier at
FLOW_RATE 1.67e-6 carries 1.0038 Q on this mesh, n0/U_b,d = 3.3937e-4 mol
s/m, a plateau of 2.2625e-5 mol/m at t_ev 15 s, within 0.4 % of the
paper's 3.38e-4 mol s/m and within -0.1 to +0.6 % of the plateaus of Fig.
9 as digitised (n' t_ev = 3.373e-4 to 3.397e-4 mol s/m for t_ev 15 to 90
s; row "facet-area rescale"; plan section 27, item 9, and section 29,
item 6).  The deposition zone
does not agree.  In this wall table 680 K lies at x = 0.476 m,
where run/014's onset is (T_dep = 679.82 K, the input Tdep to a bin), and
641 K at 0.50 m, where the published onset is.  Both T-Flows
temperature-model codes use Tdep = 638.15 K and the layer r > 0.9 R;
`Thesis` uses A 150, k 0.02 (`Project`: A 1.35, k 2.7 in its code, A 150,
k 0.02 in the reference.txt of its reference case).  The run of the
`Thesis` parameter set (the sensitivity variant
`fig9-tev15-TFlows` of M4long: this case with `Tdep 638.15`, `A 150`, `k
0.02` and `distance 2.4e-4`, 61056 layer cells, at dt 4 ms):

| Quantity | paper | T-Flows set | difference |
|---|---|---|---|
| gas at 50 % of the plateau, 16 s | 0.507 m | 0.5134 m | +6.4 mm |
| gas at 10 % of the plateau, 16 s | 0.518 m | 0.5272 m | +9.2 mm |
| solid onset (1 %), 30 s | 0.502 m | 0.5011 m (T_dep 639.14 K) | -0.9 mm |
| solid peak, 30 s | 1.14e-3 mol/m at 0.513 m | 1.62e-3 mol/m at 0.5102 m | -2.8 mm |

It puts the deposit where Figs. 9-10 show it (onset and peak within two
bins), with a gas front 6-9 mm downstream of Fig. 9's and a higher peak.
Whether Figs. 8-10 were made with Tdep of about 638 K (and the thesis
layer), or with another wall table, is a question for the authors (plan
section 26); until it is answered, this case reproduces the configuration
the paper states, not its figures.  Fig. 9's solid spike
(3.5e-3 mol/m at 0.501-0.505 m) and Fig. 10's peak (1.14e-3 mol/m at 0.513
m) show the same parameters at the same time with different shapes
(`manuscript/notes-for-authors.md`, item D).

## How to run

Never inside `run/`.  The test helper `makeCase 014 <dir>` of
`applications/rhoFixedFlowFoam/tests/functions` makes a working copy (the
small inputs copied, the mesh converted once by gmshToFoam into
`tests/work/mesh` and linked, the table linked); by hand:

    W=<work directory>
    mkdir -p $W/0 $W/constant $W/system
    cp run/014-paper-temperature-model/system/* $W/system/
    cp run/014-paper-temperature-model/constant/*Properties \
       run/014-paper-temperature-model/constant/wall_x_profile_complete.dat \
       $W/constant/
    cp run/014-paper-temperature-model/0/* $W/0/
    gmshToFoam -case $W run/001-first-pipe-flow-simpleFoam/pipe.msh
    cd $W
    setFrozenCarrier                  # rho, U, p, T, phi into 0
    decomposePar                      # 16 slabs along x
    mpirun -np 16 setFrozenCarrier -check -parallel   # continuity of phi
    mpirun -np 16 rhoFixedFlowFoam -parallel > log

`decomposePar` writes a patch whose values are all equal as a `uniform`
value of about 10 digits, in a binary file too: a processor patch, and an
ordinary patch with a single face on a rank (which a general decomposition
such as scotch can produce).  The check after `decomposePar` shows such a
loss of continuity; writing the carrier on the decomposed case instead
(`setFrozenCarrier -parallel`, after removing the decomposed carrier
fields) is exact (criteria M4.1b and M4.1d).  `setFrozenCarrier -check`
reads the carrier where the solver reads it (the `phi` that findInstance
finds from the start time: 0 of a case that has run), and
`setFrozenCarrier` writes into that instance as well.

The run to 30 s (30000 steps) takes about 0.26 s per step on 16 ranks of
an idle workstation (0.3-0.4 s on a loaded one), i.e. two to three hours.
`utilities/axialProfiles.py <case> --set gamma` prints T_dep and the
integrals.

## The runs of acceptance 3-5

`tests/Alltest M4long` (`tests/paperMode/longRuns`) runs these edits of
the case (`M4LONG_JOBS` runs of `M4LONG_RANKS` ranks at a time; serial runs
are the most efficient on a loaded machine), reuses a run whose inputs did
not change, runs again a run made with older sources (or, only with
`M4LONG_REVALIDATE=yes`, revalidates it: the binaries it was made with and
the current ones must write the same files from its state; the criterion
lines then name it; plan sections 26 and 27), resumes a stopped one, and
evaluates criteria M4.3-M4.5, M4.16 and M4.17, the closure of the runs
made with the current sources (`M4LONG_RUNS` selects the runs; plan
sections 25 to 27):

| Run | Flow | n0 [mol] | Tdep [K] | t_ev [s] | End | dt |
|---|---|---|---|---|---|---|
| Fig. 9, t_ev 6 s | 1.67e-6 m3/s | 3.16e-5 | 680 | 6 | 12 s (profiles at 7 and 12 s) | 1 ms |
| Fig. 9, t_ev 15 s (this case) | 1.67e-6 m3/s | 3.16e-5 | 680 | 15 | 30 s (profiles at 16 and 30 s) | 1 ms |
| T106_6 | 106 mL/min | 4.84e-6 | 775 | 6 | 6 s | 1 ms |
| T106_60 | 106 mL/min | 4.43e-6 | 712 | 6 | 10 s, then +2 s until T_dep at 60 s is bounded to 0.05 K | 1 ms |
| T106_300 | 106 mL/min | 2.97e-6 | 666 | 6 | 10 s, then +2 s until T_dep at 300 s is bounded to 0.05 K | 1 ms |
| T5_60 | 5 mL/min | 3.15e-6 | 788 | 6 | 60 s | 10 ms (5 ms for the dt check) |
| T540_60 | 540 mL/min | 3.93e-6 | 671 | 6 | 4 s, then +1 s until bounded | 0.22 ms |
| T-Flows set (sensitivity, not a gate) | 1.67e-6 m3/s | 3.16e-5 | 638.15, with A 150 1/s, k 0.02 1/K, distance 2.4e-4 m | 15 | 30 s (profiles at 16 and 30 s) | 4 ms (Courant 0.45; a 4 ms run of this case reproduces the 1 ms gas front to 4 digits) |

Each Table 2 run is made for both layer rules.  The deposit of the
temperature model is irreversible, so after the gas inventory G has
fallen far enough the deposit at the end of the experiment is bounded by
the deposit now and the deposit now plus G in any one bin; the runs stop
when that bound fixes T_dep to 0.05 K (plan section 25).  The t_ev = 6 s
run is no curve of the paper (Fig. 9 has none): with the 15 s run it
checks that the final deposit does not depend on t_ev (acceptance 4, a
property of the model).

Acceptance 5 is gated as the researcher restated it (2026-10-02):
|T_dep - T_paper| <= max(2 K, dT_bin), dT_bin the drop of the wall
temperature across the native bin that holds the onset (2.4-2.8 K here).
All ten lines meet it; T5_60 (+2.25 K with layer distance, +2.22 K with
faceCells; dT_bin 2.77 K) is outside the plain +-2 K, printed for
information.  Its offset has two parts.  +1.25 K is the bin placement
relative to the input Tdep = 788 K: the last bin without deposit has its
wall at 789.37 K, the first with deposit already 22 % of the peak, so the
interpolated 1 % crossing lies 1.25 K above Tdep, on a wall falling by 19
K/cm upstream of x = 0.41 m.  The other +1.0 K is the paper's own value
lying 1 K below its input (787 against 788 K), a cause this case cannot
establish.  The set-up of this case also differs from the paper's T5_60
run in ways that can move the onset by about that much: the wall table
(the paper measured one for 5 mL/min, Fig. 14 iii, about 15 K/cm at 788
K), the source position (about 0.15 m in the paper, 0.04-0.05 m here), the
axial resolution (1.627 mm here, 1.533 mm in the T-Flows temperature-model
mesh) and the state compared (the paper's curve holds about the whole
inventory, a complete deposit; run/014's T5_60 at 60 s holds 33 % on the
wall and 67 % in the gas, spread back to the closed inlet).  The evidence
does not single out one of them (plan sections 26 and 27).

Whether T5_60 is evaluated at 60 s or at completion changes the shape of
its profile, not its T_dep: no layer cell upstream of x = 0.41009 m (the
bins before the onset) is below Tdep (the coldest is at 789.38 K), so no
deposit can ever form there, and from 30 to 60 s the deposit grows
self-similarly (its peak per deposited fraction constant to 1 %, the
first bin with deposit at 22 % of the peak throughout), so at completion
the 1 % crossing stays where it is, to about 0.01 K (plan section 27,
item 10).

## Tests

`out_excerpt` holds the log of the first 5 steps (`endTime 0.005`) of this
case with the carrier of `setFrozenCarrier`, written on OpenFOAM v2412 by
the build of the M4 review, round 3 (it differs from that of round 2 only
in the digits of `solver defect` and `solverDefect`, the face-flux
correction now booked in flux form, and in the max mole fraction of the
monitor, now from the carrier's p/(R T), 25 times the old value: rho = 1
is no density of helium; plan section 29; round 2 differed from round 1 in
the digits of the ledger's booking and the start-up count of interface
faces with elements of two or more other ranks); `tests/Alltest M4`
(criterion M4.11) checks that the current build reproduces it.
Criterion M4.2 checks the layer counts of both rules (30528 cells, A_I/V_I
= 8533.7 1/m; 15264
cells, 17681.9 1/m) and that the T-Flows centre-distance rule selects the
same cells; M4.7c runs the case (background Y = 0.1, 10 steps of 4 ms) in
serial and on three 4-rank splits, M4.7e restarts its states onto another
decomposition and M4.7f on the same one; M4.15 runs it with 20 steps of
0.1 s on 4 ranks and requires the closure not to drift (at most 2e-15 of
the reference); M4.12 checks the provenance of the wall table and of the
conductivity table against the T-Flows revision.
