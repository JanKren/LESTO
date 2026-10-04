# Case 015 — HKS exchange with p_eq from GEMS3K (mode local)

Development case of milestone M9 of
`applications/rhoFixedFlowFoam/doc/phase-change-plan.md` (section 32).  It
is case 013 (the Hertz-Knudsen-Schrage law with an evaporating inventory
sample on the 007 carrier) with the equilibrium vapour pressure of the law
taken from a Gibbs-energy minimisation with GEMS3K, per step and element,
instead of the table: `HKSCoeffs { equilibrium GEMS; }`, `GEMSCoeffs { mode
local; }`.  The table stays: it checks every GEMS3K result and replaces it
where a guard applies.

## Build

The case needs `rhoFixedFlowFoam` built with the GEMS3K bridge:

    export LESTO_GEMSBRIDGE=<LESTO>/thermochemistry/gemsbridge
    cd applications/rhoFixedFlowFoam && ./Allwmake

(`libgemsbridge.so` is git-ignored; `thermochemistry/README.md` builds it.)
A build without the bridge stops at start-up with a FatalIOError and this
hint.  The tests build such a copy in their work area (`tests/Alltest M9`).

## Carrier, mesh and species

As case 013: the frozen helium carrier of case 007 (time 261, linked from
`0/`; 28 mL/min, a development carrier, plan section 8, risk 10), the pipe
mesh of case 001 (linked `constant/polyMesh`), `PbI2_g` with its condensate
`PbI2_s`, the HKS law in SI on every face of `WALL` (`layer faceCells`), the
boat (the 144 cells of the cylinder x = 0.04-0.05 m, radius 0.48 mm) with
`n0 = 4.84e-6 mol` removed at 6 s, `writeFrozenFields no`.  The carrier is a
gas at its p and T to 0.3 % (`max |T - p W/(R rho)|/T = 0.0030285`), as
mode local requires (paper mode would be refused unless
`allowNonPhysicalCarrier yes`).

## GEMS3K

`constant/gems` is a link to `thermochemistry/systems/PbI2He`, the PbI2-He
system written by `thermochemistry/systems/make_gems3k_input.py` (Gurvich
(1991) via the NASA-9 records: an ideal gas He, I, I2, Pb, PbI, PbI2 and the
pure condensates PbI2(cr), PbI2(l), Pb(cr), Pb(l); T 270-1250 K, P 0.8-1.2
bar).  The pair's `gems` block names the gas species `PbI2(g)`, its
condensates `PbI2(cr)` and `PbI2(l)` (p_eq over the most stable of them)
and its formula `Pb 1 I 2`, which the start-up checks against the system;
`excess { I 1e-6; }` makes the bulk composition slightly iodine rich, since
the interior-point method of GEMS3K occasionally stalls at the exact PbI2
stoichiometry (plan section 32, item 3).

Every condensed species is suppressed, so a call returns the gas at the
composition of the cell, and `p_eq = p_PbI2 10^(-max log10 Omega)`: for this
ideal gas the table's value to the G0 interpolation of GEMS3K (1e-5 to 1e-4
decades).  Mode local evaluates every element at the start of every step
(`updateInterval 1`) at its temperature, the cell's pressure and the
cell's composition (the carrier p/(R T) and the gas rho Y/M), warm started
from its own state of the previous update.  The guards take the table
instead below 500 K (`minTemperature`), below the mole fraction 1e-6
(`minMoleFraction`), outside the grid of the system, for a failed call, and
for a result more than 0.01 decades off the table (rejected).  With the
default `speciation none` the law is that of case 013; `speciation lagged`
would also take the fraction of the gas present as PbI2(g) molecules from
GEMS3K (a separate modelling change).

## Numerics

As case 013 (`dt = 1 ms`, Euler, `Gauss SuperBee`, outer loop `nCorr 10`,
`tolerance 1e-5`, loose GAMG correctors, `writeFormat binary`), but the
final solve `Y_PbI2_gFinal` at **1e-14**, the final tolerance of HKS cases
decided in plan section 28, item 8 (case 013 too since the M7 integration,
plan section 33): the solver defect of the first 25 steps stays at 2e-16
of the inventory (3e-14 with 1e-12), for about 20 % more time per step.

## Output

As case 013, and every step a GEMS line after the exchange line:

    Y_PbI2_g GEMS local: update (calls 576: warm 0, retried 0, cold 576;
      ... iterations per call; failures 0, engine re-created 0); elements:
      GEMS 576, table 14832 (below minTemperature 4356, below
      minMoleFraction 10476, outside the GEMS grid 0, rejected 0, failed 0);
      max |dlog10 p_eq/p_eq,table| ...

At the start no element has gas: the update of step 1, at the composition
of t = 0, makes no call (the first GEMS line of the excerpt: every element
on the table, 4356 wall faces below 500 K and the rest without gas).  From
step 2 on the boat's 144 cells and the wall cells around it are evaluated
(576 to 720 calls per step in steps 2 to 5), the other 14700 elements are
guarded: 4356 wall faces below 500 K (the cold end), the rest without gas.
The end of the run prints the summary with the mean time of a warm call
(60-90 us here on the loaded workstation).  `gemsLog/` in the case (`processorN/gemsLog/` in parallel)
holds the log of every rank (`gemsEquilibrium.log`: the engine, every
update, every failure) and GEMS3K's `ipmlog.txt`; it is git-ignored.
`writeEquilibriumField yes` adds the fields `pEq_PbI2_g` (the p_eq of every
element) and `gemsSource_PbI2_g` (where it came from).

## Tests

`out_excerpt` holds the first 5 steps, written on OpenFOAM v2412 by the
build of plan section 32 with the bridge `thermochemistry/gemsbridge` of
that workstation (`libgemsbridge.so` with the md5 fb67103e...); the lines
after its `...` (the end-of-run summary) are not compared, since the
summary holds the timing of the GEMS3K calls.  `tests/Alltest M9`
(criterion M9.8) checks that the current build reproduces it: with that
bridge (its md5 is `EXCERPT_BRIDGE_MD5` in `tests/gems/Allrun`) the whole
excerpt exactly.  With another bridge (another compiler, GEMS3K version or
libm; on Merlin7 it is rebuilt) only the start-up and step 1 are compared
exactly, since no GEMS3K call reaches them (no element has gas at t = 0),
and the rest is NOTRUN: GEMS3K's p_eq then differ in their last bits, and
from step 2 on every round-off token of the solver moves with them (the
GAMG final residuals, the solver defect, clamped, solverDefect, min Y),
not only GEMS3K's own numbers, the IPM iterations per call and the
largest `|dlog10 p_eq/p_eq,table|`; masking those (the first cleanup, plan
section 35, item J) did not suffice (item V).  An excerpt written again
with another bridge needs that bridge's md5 in `EXCERPT_BRIDGE_MD5`.  M9.8
also runs the case for 25 steps at its own settings: a GEMS line in every
step, 0 failures, closure below 1e-13, solverDefect below 1e-11.  The full
run (16 s, about 2.3 s per step serial, of which 0.05 s for GEMS3K in the
first steps) is meant for a parallel machine.
