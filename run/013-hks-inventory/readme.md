# Case 013 — Hertz-Knudsen-Schrage exchange with an evaporating sample

Development case of milestone M5 of
`applications/rhoFixedFlowFoam/doc/phase-change-plan.md`.  It replaces the
temperature model of case 012 by the Hertz-Knudsen-Schrage (HKS) law with
the phase-resolved vapour pressure of PbI2, and the release of case 012 by
an **inventory** sample: solid PbI2 in the boat evaporates by the same law
and is pulled out at `removeTime`.

## Carrier and mesh

The frozen helium carrier of case 007 (time 261, linked from `0/`) and the
pipe mesh of case 001 (linked `constant/polyMesh`), as in case 012.  This
carrier is 28 mL/min, and its `T`/`U` were overwritten in 3ce3791, so it
is a development carrier only (plan section 8, risk 10); the log prints its
consistency, `max |T - p W/(R rho)|/T = 0.0030285`.

`writeFrozenFields no`: the carrier is read from `0` and never written; the
write directories hold only `Y_*`, `c_*`, `mw_PbI2_s`, `mDep_PbI2_s`,
`Cs_PbI2_s` and `uniform/`.

## Species and model

- `PbI2_g`: gas, PbI2-He diffusivity, `molarMass 0.46100894` (required for
  a paired gas); `zeroGradient` on WALL (the exchange enters the gas as a
  source).
- `PbI2_s`: its condensate.  `Y_PbI2_s` is derived output, `(deposit*area
  + reservoir*V)/(rho*V)`.  The state is the deposit of every wall face
  (`uniform/phaseChangeDeposits`, shown as `mw_PbI2_s` [kg/m2]) and the
  reservoir of every boat cell (`uniform/phaseChangeReservoirs`, shown as
  `Cs_PbI2_s` [kg/m3]).

HKS law (plan E2, SI; `constant/thermochemistryProperties`, `HKSCoeffs`:
`accommodation 1`, `Ce 1`, `kineticScale 1`) on every face of `WALL` (one
interaction layer, `layer faceCells`, `A_I/V_I = 17681.9 1/m`):

    mdot'' = G (p_eq(T) - rho R T Y/M),   G = 2 sigma/(2 - sigma) sqrt(M/(2 pi R T))

A bare wall face only deposits; a face with a deposit also evaporates, at
most its deposit per step.  At 700 K, `G R T/M = 89.65 m/s`; the start-up
line prints `lambda = sum_e g_e/(rho_c V_c)` of the wall cells, 1.0e6 to
1.9e6 1/s, i.e. `lambda dt` = 1000 to 1900 at 1 ms: the exchange is stiff.
It is solved implicitly (`fvm::Sp` plus `Su`), and a local implicit
predictor chooses the regime of every face at every corrector.

p_eq: `constant/pv_PbI2_phases.csv`, a link to
`thermochemistry/systems/data/pv_PbI2_phases.csv` (written by
`thermochemistry/systems/make_pv_table.py` from the NASA-9 records of
Gurvich (1991)): log10 p over PbI2(cr) and over PbI2(l) every 5 K from 270
to 1250 K plus the melting point 683 K, `p_eq` = the minimum of the two
phases (the stable one), log10 p interpolated linearly in 1/T.  The
start-up line reports the phase change at 683 K.  Above 683 K the stable
liquid has a lower `p_v` than the crystal extrapolated through the melting
point (2 to 7 times lower in the evaporation zone; `thermochemistry/
README.md`).

Sample `boat` (`mode inventory`): the 144 cells of the cylinder x =
0.04-0.05 m, radius 0.48 mm (selected in memory; no `topoSet` on the linked
mesh) hold `n0 = 4.84e-6 mol` (2.23 mg, the T-Flows `MASS_SAMPLE`; T106_6)
as a reservoir `Cs = n0 M/V_Z = 274.9 kg/m3`, filled on a fresh start.
Every cell evaporates by the HKS law with the interfacial area per volume
`areaPerVolume 4000 1/m` (2/R_sample; the T-Flows factor 5 is not applied)
and never takes up gas.  At `removeTime 6 s` (the first step that starts at
t >= 6 s - dt/2 on the ledger's exact clock) the boat is pulled out: what
is left of the reservoir is booked as `removed`, once.

The radius lies between two rings of cell centres.  In 0.04 < x < 0.05 m
the centres lie on rings at r = 0.1, 0.265, 0.3, 0.436 and 0.458 mm (144
cells) and, next, at 0.5 mm (18 cells; computed 5.000000000000005e-4 to
...7e-4 m).  A radius of exactly 0.5 mm put the selection on that ring:
cylinderToCell's strict test excluded it by 2e-12 relative, so a mesh
rewritten at lower precision (decomposePar or renumberMesh with
`writeFormat ascii`) selected 146 cells and a restart stopped as 'changed'.
0.48 mm selects the same 144 cells with a margin of 4 % to both rings.  An
inclusive centre test r <= 0.5 mm on the exact centres, as the T-Flows
parity checklist words it, would include that ring (162 cells); the
code-to-code case of M6 must pin which of the two it compares (plan
section 11).

## Numerics

As case 012: `dt = 1 ms`, implicit Euler (required), `Gauss SuperBee`,
outer loop `nCorr 10`, `tolerance 1e-5` (a step is converged when the
initial residual is below the tolerance and the regime changes carry no
significant mass), loose correctors (`Y_PbI2_g`: GAMG 1e-8, relTol 0.1)
and one tight final solve (`Y_PbI2_gFinal`: 1e-14, relTol 0), no
relaxation, `writeFormat binary` (exact restarts, also across the removal).
The final tolerance 1e-14 is that of HKS cases (plan section 28, item 8;
1e-12 until the M7 integration, plan section 33): with 1e-12 the booked
solver defect grew by about -2e-18 kg per step once the boat was empty,
to 2e-8 to 6e-8 of the inventory after 16 s (the residual of the final
solve, one sign); 1e-14 cuts it a hundred times, for about 20 % more time
per step.

## Output

- Every step: the exchange line (regimes of the wall faces and, in
  brackets, of the boat cells; outer loop; regime changes of the deciding
  corrector and of the final predictor; guard; projections;
  complementarity), the transport line, the monitor line (min Y, negative
  mass, carrier-validity mole fraction: the gas near the boat reaches
  x = 0.1 and more, where the passive-carrier assumption fails; plan
  section 8, risk 5) and the balance line

      Balance PbI2_g [mol]: gas ... wall ... sample ... released 0
      transport ... removed ... clamped ... restart ... solverDefect ...
      closure ...

  whose closure stays at round-off;
- `postProcessing/phaseChangeBalance/0/PbI2_g.dat`: the ledger, every step;
- at write times: `Y_PbI2_g`, `Y_PbI2_s`, `c_*`, `mw_PbI2_s`,
  `mDep_PbI2_s`, `Cs_PbI2_s`, and in `uniform/` the ledger, its layout (with
  the record of the boat and its removal) and the element state.

For axial profiles and T_dep add a `profiles` dictionary as described in
`run/012-temperature-model/readme.md` (the observable `wallPlusGas` counts
the gas outside the boat cells).

## Tests

`out_excerpt` holds the first 5 steps written on OpenFOAM v2412 by the
build of the M7 integration (plan section 33), with the final tolerance
1e-14.  Against the excerpt of the M4 review, round 3 (final tolerance
1e-12), 61 of its 1656 tokens differ (the tokens that identify a run,
dates and timings, aside), all at or below the residual of the old final
solve: the final solve's residual and iterations (10 tokens), the INLET
and OUTLET transport of the transport line (10; 1.6e-18 and 1e-60 kg/s,
by up to 1e-3 relative), the solver defect of that line (5), the
cumulative `transport` of the balance line (5), `clamped` (5; the
rounding of the stored condensate), `solverDefect` and its relative value
(10; 2.6e-16 instead of 3.0e-14 of the inventory after 5 steps), the
closure and its relative value in 3 steps (6; by one ulp of the
inventory), and the minimum Y of the monitor and min/max lines (8) and the
negative mass (2) of the overshoot (-6e-7 to -5e-9, by up to 2e-4
relative); the inventories, the max mole fraction and every other token
are unchanged.  All 61 come from the tolerance: the build of the cleanup
of plan section 35 run with the old final tolerance 1e-12 reproduces the
old excerpt (`regress.sh -excerpt`: IDENTICAL, 111 lines, 5 steps; plan
section 35, item R).  `tests/Alltest M5` (criterion M5.3b) checks that
the current build reproduces it, and criterion M5.3a runs this case for 25
steps of 4 ms with a background gas `Y = 0.1` (so that the wall deposits
from the first step): closure below 1e-13, solver defect below 1e-11, no
significant regime change in the deciding corrector or the final
predictor, no significant projection, and the boat's reservoir decreasing
in every step.  Criterion M5.2g restarts that state on 4 ranks
(hierarchical 1 2 2) and back in serial.  The full run (16 s, 16000 steps
of about 1.5-2.5 s serial) is meant for a parallel machine.
