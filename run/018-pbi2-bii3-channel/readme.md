# Case 018 — Two independent pairs: PbI2 and BiI3

This case is based on case 014.

The channel, carrier, numerics and the PbI2 pair are those of case 014.  A
second gas/condensate pair, BiI3_g/BiI3_s, is released from the same boat.
The two pairs do not interact, so the PbI2 results must be those of case
014 byte for byte, and each pair deposits where its own vapour pressure is
met: PbI2 near 561 K, BiI3 near 450 K.  The element ledger now has three
elements, Pb, Bi and I.

The source ratio follows the inferred PbI2/BiI3 deposits of Liu et al.
(2025), LBE-I_SS_II, scaled to an iodine mole fraction of 3.9e-6 in this
carrier (`benchmark.json`).  This is a qualification case, not the paper's
geometry or a prediction: KI, metal vapour and source speciation are left
out, and measured T(x), flow conditions and the release history would be
needed.

## Differences from case 014

A second pair in `constant/thermochemistryProperties`

    BiI3_g {condensed BiI3_s; HKS {accommodation 1; Ce 1; kineticScale 1;}
            vapourPressure {file "$LESTO_M10_DATA/pv_BiI3_phases.csv";
            phases ("BiI3(cr)" "BiI3(l)"); units log10Pa; outOfRange fatal;}}

and a second amount from the boat

    amounts {PbI2_g 2.7992536908413373e-11;
             BiI3_g 3.5684583649025576e-11;}      // mol over 2 s

What follows from it: the species BiI3_g (D 1e-4 m2/s, `formula {Bi 1;
I 3;}`) and BiI3_s, their initial fields, and their entries in
`system/fvSchemes` and `system/fvSolution`.

## Modified files

    constant/speciesTransportProperties, constant/thermochemistryProperties
    0/Y_BiI3_g, 0/Y_BiI3_s           (new)
    system/fvSchemes, system/fvSolution
    benchmark.json, verification.json (new: inputs and the recorded run)

## Run

    export LESTO_M10_DATA=/external/m10-tables   # holds pv_BiI3_phases.csv
    ./Allrun            # blockMesh, makeCarrier and the 3 s (about 10 s),
                        # then the numbers of expect.json
    ./Allclean

`pv_BiI3_phases.csv` is made by `thermochemistry/systems/make_pv_tables.py`
from restricted Bi thermodynamic records; neither is redistributed.
Without `LESTO_M10_DATA`, `./Allrun` prints NOTRUN.

## What to check

At 3 s (`postProcessing/axialProfiles/3/*_native.dat`, `verification.json`):

| pair | 1 % onset | plug-flow saturation crossing | difference |
|---|---|---|---|
| PbI2 | 561.164 K | 567.057 K | -5.89 K |
| BiI3 | 449.979 K | 450.383 K | -0.40 K |

- Element closures in every step: Pb 4.05e-16, Bi 4.56e-16, I 4.37e-16;
  the released amounts equal the supplied ones to 1e-14.
- The PbI2 ledger, profiles and the element-Pb ledger equal those of case
  014 byte for byte (`./Allrun` compares them when case 014 has been run).
- The onsets differ from the peak temperatures; `Tpeak` is in the profile
  headers.

Measured again on 2026-10-14 with the M10e solver and the tables of
`m10-design-2026-10-04` (pv_BiI3_phases.csv md5 0429d3fc...): every number
above is reproduced exactly.
