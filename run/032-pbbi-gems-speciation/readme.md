# Preliminary M10c GEMS gas re-speciation channel

This case uses the prescribed PbI2/BiI3 source and provisional carrier of
run/031, with GEMS3K standard-state data for all nine gas species. Its
`frozen` mode caches formation constants at the frozen cell temperatures
and pressures, then uses the exact gas kernel every step. Select `local`
in the `respeciation/GEMS` block to exercise per-cell warm GEMS3K calls,
physical guards and kernel refinement/fallback.

Build the bridge and opt into it when building the solver:

```sh
make -C ../../thermochemistry/gemsbridge check
export LESTO_GEMSBRIDGE="$PWD/../../thermochemistry/gemsbridge"
../../applications/rhoFixedFlowFoam/Allwmake
export LESTO_M10_SYSTEM=/external/PbBiIHe/PbBiIHe-dat.lst
export LESTO_M10_DATA=/external/m10-tables
./Allrun
```

`LESTO_M10_DATA` holds the external BiI3 vapour-pressure table and the
formation table generated for run/031. The GEMS system must contain the
carrier plus exactly the participating reactive gases; formulas and molar
masses are checked. No licensed Bi records are bundled. The default solver
build remains bridge-free and gives a rebuild hint if this engine is selected.

The carrier ramp, common diffusivity and release history remain provisional.
This is numerical qualification of the implementation. Experimental inputs
and later reactive-condensate stages are still required for validation.
`tests/Alltest M10c` qualifies the ABI, reduced map, serial/MPI closure,
restarts, input guards, faults, HKS dual pressures and deposit comparison
with run/031. The `native` profiles use 1 mm bins; `cm` profiles use 1 cm bins.

The preliminary three-mode study is recorded in `verification.json`.
Frozen-versus-open-kernel deposit L1 differences in 1 cm bins are
1.0332e-10 (PbI2) and 1.0637e-7 (BiI3); local-versus-frozen differences
are zero. The maximum element closure is 6.7e-16. Local calls including
both refinements average 132.2 us over 5,814,420 calls on this workstation.
These are numerical results under the prescribed provisional inputs.

`qualification.json` preserves the reduced-map raw IPM errors, guard counts
and reference diagnostics by dilution and individual gas mole-fraction
decade. These raw errors explain why local results require the exact kernel
refinement. The combined M0/M9/M10c suite passed 78 checks with zero failures
and two documented skips; all nine M10c gates passed. Its scientific gate
rechecked the completed three-mode study; other M10c gates used a fresh
solver build. Build and skip details are recorded in `verification.json`.
