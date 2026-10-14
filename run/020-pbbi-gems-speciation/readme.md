# Case 020 — Gas re-speciation with GEMS3K

This case is based on case 019 (and no later case builds on it: GEMS3K is
an optional build).

Case 019 re-speciates the gas with an exact kernel and formation constants
from a table.  Here GEMS3K supplies the standard states of all nine gases.
Mode `frozen` caches the formation constants at the frozen cell
temperatures and pressures and then uses the exact kernel every step; mode
`local` calls GEMS3K per cell and step (warm started), with physical guards
and the kernel as refinement and fallback.  For this ideal gas both must
reproduce case 019.

## Differences from case 019

In `constant/thermochemistryProperties`

    respeciation {engine kernel; ...}

became

    respeciation {engine GEMS; ...
      GEMS {mode frozen; system "$LESTO_M10_SYSTEM"; carrier "He(g)";
            minMoleFraction 1e-12; maxMoleFraction 3.2e-4;
            warmIterationLimit 200;} ...}

The GEMS3K system must hold the carrier and exactly the nine gases;
formulas and molar masses are checked at start-up.

## Modified files

    constant/thermochemistryProperties
    benchmark.json, verification.json, qualification.json

## Run

    make -C ../../thermochemistry/gemsbridge check
    export LESTO_GEMSBRIDGE="$PWD/../../thermochemistry/gemsbridge"
    ../../applications/rhoFixedFlowFoam/Allwmake   # the bridge build
    export LESTO_M10_SYSTEM=/external/PbBiIHe/PbBiIHe-dat.lst
    export LESTO_M10_DATA=/external/m10-tables     # as case 019
    ./Allrun            # mode frozen, 3 s (a few minutes), expect.json
    ./Allrun local      # mode local (about 15 min), against frozen
    ./Allclean

Point `LESTO_SOLVER` to the bridge build.  Without it, or without the
external data, `./Allrun` prints NOTRUN.  `tests/Alltest M10c` qualifies the
ABI, the reduced map, serial/MPI closure, restarts, guards and faults.

## What to check

From `verification.json` (the three-mode study):

- 1 cm deposit L1, frozen against the kernel of case 019: 1.0332e-10
  (PbI2) and 1.0637e-7 (BiI3); local against frozen: 0.
- Element closure 6.7e-16.
- Local calls, both refinements included: 132.2 us on average over
  5,814,420 calls on this workstation.

`./Allrun` checks the closures and, when case 019 has been run, the L1
against its deposits; `./Allrun local` checks the L1 against frozen.
`qualification.json` keeps the reduced-map IPM errors, guard counts and
reference diagnostics by dilution and mole-fraction decade: they show why
local results need the exact kernel refinement.
