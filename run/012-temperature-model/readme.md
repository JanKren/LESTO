# Case 012 — Temperature-model deposition with a sample release

Development case of milestone M2 of
`applications/rhoFixedFlowFoam/doc/phase-change-plan.md`.  It replaces the
mock thermochemistry of case 010 by the gas/condensate exchange of
`rhoFixedFlowFoam` (`constant/thermochemistryProperties`, `model
temperature`).

## Carrier and mesh

The frozen helium carrier of case 007 (time 261, linked from `0/`) and the
pipe mesh of case 001 (linked `constant/polyMesh`).  This carrier is
28 mL/min, and its `T`/`U` were overwritten in 3ce3791, so it is a
development carrier only (plan section 8, risk 10); the log prints its
consistency, `max |T - p W/(R rho)|/T = 0.0030285`.

`writeFrozenFields no`: the carrier is read from `0` and never written; the
write directories hold only `Y_*`, `c_*`, `mw_PbI2_s` and `uniform/`.

## Species and model

- `PbI2_g`: gas, PbI2-He diffusivity, `molarMass 0.46100894` (required for
  a paired gas).  WALL must be `zeroGradient`: the exchange enters the gas
  as a source, not through the boundary.
- `PbI2_s`: its condensate.  With `model temperature` `Y_PbI2_s` is
  derived output, `(deposit*area)/(rho*V)`; the state is the deposit of
  every wall face, kept exactly in `uniform/phaseChangeDeposits` and shown
  as `mw_PbI2_s` [kg/m2] on the WALL patch (per effective element area,
  here the face area, since `areaMultiplier` is 1).

Temperature model (Lobresco et al. 2026, Eq. 5) as a wall flux of every
face of `WALL` (one interaction layer, `layer faceCells`, `areaMultiplier
1`: the multiplier would scale the applied rate too, and the start-up line
prints it as `applied rate areaMultiplier*v_d*A_I/V_I = 220 1/s`):

    mdot'' = -rho v_d [1 - exp(-k (Tdep - T))]^+ Y     (deposition, irreversible)

with `A = 220 1/s`, i.e. `v_d = A V_I/A_I = 0.0124421 m/s` on this layer
(`A_I/V_I = 17681.9 1/m`), `k = 8e-3 1/K`, `Tdep = 680 K`.

Sample `boat`: the 144 cells of the cylinder x = 0.04-0.05 m, r <= 0.5 mm
(selected in memory; no `topoSet` on the linked mesh) release
`n0 = 3.16e-5 mol` of `PbI2_g` uniformly over `t_ev = 6 s` from `t = 0`.

## Numerics

- `dt = 1 ms` (plan section 13; Courant number about 0.12 on this carrier),
  implicit Euler (required), `Gauss SuperBee` for `PbI2_g`.
- Outer loop: `nCorr 10`, `tolerance 1e-5`; a step is converged when the
  initial residual is below the tolerance and the regime changes of the
  exchange elements carry no significant mass.  Loose correctors
  (`Y_PbI2_g`: GAMG 1e-8, relTol 0.1) and one tight final solve
  (`Y_PbI2_gFinal`: 1e-12, relTol 0); no relaxation.
- `writeFormat binary`: restarts continue the mass ledger and the gas
  exactly (the deposits are exact in any format).  In parallel,
  `writeFrozenFields no` keeps the carrier's `phi` out of the restart
  directories, so a binary restart repeats the continuous run bit for bit.

## Output

- Every step, one exchange line, one transport line (per patch), one
  monitor line (min Y, negative mass, carrier-validity mole fraction) and
  the balance line

      Balance PbI2_g [mol]: gas ... wall ... released ... transport ...
      solverDefect ... (relative ...) closure ... (relative ...)

  whose closure stays at round-off: the cumulative ledger entries are
  compensated sums, so it does not grow with the number of steps (the
  6 s release window of this case at dt = 1 ms, 6100 steps, closes to
  1.8e-14 and releases n0 M to 2e-16 on the test channel, criterion M2.2g
  of `tests/Alltest M2`);
- `postProcessing/phaseChangeBalance/0/PbI2_g.dat`: the ledger, every
  step, in kg with 17 digits;
- at write times (and whenever a function object writes, e.g.
  `runTimeControl`): `Y_PbI2_g`, `Y_PbI2_s` (derived), `c_PbI2_g`,
  `c_PbI2_s`, `mw_PbI2_s`, `mDep_PbI2_s` (the deposit per wall area,
  milestone M3; equal to `mw_PbI2_s` here, since `areaMultiplier` is 1),
  and in `uniform/` the ledger `phaseChangeBalance`, its
  `phaseChangeLayout`, and the element state `phaseChangeDeposits` and
  `phaseChangeRegimes` (always binary, exact).

## Axial profiles (milestone M3)

The case dictionaries have no `profiles` block, so that the log keeps
reproducing `out_excerpt` (criterion M2.8).  For the axial profiles of
the paper (Figs. 9-14), add to `constant/thermochemistryProperties`

    profiles {
      bins {
        native { min 0; max 0.69; nBins 424; }   // the axial cells
        gamma  { min 0; max 0.69; width 0.01; }  // the 1 cm gamma-scan bins
      }
      times (6 7 12 16);          // [s], in addition to the write times
    }

(defaults: axis `(1 0 0)`, onset `fraction 0.01` with both search modes
`fromInlet` and `fromPeak`, observable `wallPlusGas`).  The solver then
writes `postProcessing/axialProfiles/<time>/PbI2_g_{native,gamma}.dat`
[mol/m] and two log lines per set with T_dep; the profiles integrate to
the inventories of the ledger to round-off.  Read them with

    ../../applications/rhoFixedFlowFoam/utilities/axialProfiles.py . \
      --set gamma --plot profiles.png

Criterion M3.7 of `tests/Alltest M3` runs a copy of this case with such a
dictionary for the five steps of `out_excerpt`; its log without the
profile lines reproduces the excerpt.

The full run (16 s, 16000 steps of about 2.5 s serial) is meant for a
parallel machine.  `out_excerpt` holds the first 5 steps written on
OpenFOAM v2412 by the build of the M4 review, round 3, which books the
solver defect of the transport-only matrix in flux form, its face-flux
correction included, plus the realised exchange, and computes the mole
fraction of the monitor from the carrier's p/(R T).  The earlier builds
wrote the same log with other round-off digits of the solver defect,
`clamped`, `solverDefect` and the closure (M5, M4 rounds 2 and 3), a max
mole fraction 1.6e-4 higher in relative terms (the form (Y/M)/(Y/M +
1/M_carrier); the frozen carrier of run/007 deviates from the perfect-gas
law by up to 0.3 %; plan section 29, item 4), and the M2 build the line
'Cells adjacent to WALL'; `tests/Alltest M2` (criterion M2.8) checks that
the current build reproduces it.
