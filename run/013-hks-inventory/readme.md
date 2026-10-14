# Case 013 — Evaporating sample and Hertz-Knudsen-Schrage exchange

This case is based on case 012.

The deposit of case 012 follows an empirical temperature law and the PbI2
is injected at a prescribed rate.  Here both ends are physical: solid PbI2
in the boat evaporates, and the gas condenses on the wall, by the same
Hertz-Knudsen-Schrage (HKS) law with the equilibrium vapour pressure of
PbI2.  The carrier, mesh, species and numerics are those of case 012.

## Differences from case 012

The exchange model in `constant/thermochemistryProperties`

    model temperature;
    temperatureCoeffs { A 220; k 8e-3; Tdep 680; }

was replaced by

    model HKS;
    HKSCoeffs { accommodation 1; Ce 1; kineticScale 1;
                equilibrium table; speciation none; }

with, for the pair `PbI2_g`, the vapour-pressure table

    vapourPressure { file "<constant>/pv_PbI2_phases.csv";
                     phases ("PbI2(cr)" "PbI2(l)");
                     units log10Pa; interpolation logInverseT; }

p_eq is the lower of the crystal and liquid curves (the stable phase),
from the NASA-9 records of Gurvich (1991); the phases cross at 683 K.

The sample was changed from a release

    mode release; amount 3.16e-5; duration 6;

to an inventory that evaporates by the HKS law and is pulled out at 6 s

    mode inventory; amount 4.84e-6;   // mol, 2.23 mg (T106_6)
    areaPerVolume 4000;               // 1/m, 2/R_sample
    kineticScale 1; Ce 1; removeTime 6;

and its radius from 5e-4 to 4.8e-4 m (the same 144 cells, but away from
the next ring of cell centres, so that a rewritten mesh selects them too).

The final tolerance of the PbI2 solve is 1e-14 instead of 1e-12
(`Y_PbI2_gFinal` in `system/fvSolution`): with 1e-12 the solver defect of
the ledger grew to 2e-8 of the inventory over 16 s (plan section 28, item 8).

## Modified files

    constant/thermochemistryProperties
    constant/pv_PbI2_phases.csv   (new, a link into thermochemistry/)
    system/fvSolution
    0/Y_PbI2_s                    (comment only)

## Run

    ./Allrun            # 5 steps into work/excerpt, compared with
                        # out_excerpt (about 30 s)
    ./Allrun full 16    # the 16 s case on 16 ranks (about 1.5-2.5 s
                        # per step serial)
    ./Allclean

The mesh is converted once from run/001's `pipe.msh`; the carrier is case
007 at time 261 (links in `0/`).

## What to check

At start-up:

    Pair PbI2_g: ... 198 nodes 270-1250 K ... stable phase
      PbI2(cr) -> PbI2(l) at 683 K
    Sample boat (pair PbI2_g): 144 cells ... inventory of 4.84e-06 mol
      (Cs = n0 M/V_Z = 274.864 kg/m3 on a fresh start) ...
    ... lambda = sum_e g_e/(rho_c V_c) 1.03917e+06-1.89388e+06 1/s

lambda dt is 1000 to 1900: the exchange is stiff and is solved implicitly.

After 5 steps (`out_excerpt`, reproduced exactly, 111 lines):

    Balance PbI2_g [mol]: gas 1.04761e-07 wall 0 sample 4.73524e-06 ...
      closure ... (relative -1.9e-16)

The wall stays empty (the walls near the boat are far above the saturation
temperature), the sample decreases in every step, and the closure stays at
round-off.  In the full run the boat is removed at 6 s; what is left is
booked once as `removed`.

Details: notes.md (the HKS law and its regimes, the ring-of-centres
argument, the tolerance study, the excerpt's history, and the tests M5.2g
and M5.3a-d that use this case).
