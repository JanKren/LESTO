# Case 012 — Temperature-model deposition with a sample release

This case is based on case 011.

Case 011 moves PbI2 from the gas to the solid with a mock rate that only
tests the coupling.  Here the exchange is a wall law: every face of `WALL`
takes up the gas next to it with the temperature-based model of Lobresco
et al. (2026), Eq. 5, and the PbI2 comes from a boat that releases a given
amount over a given time.  Every kilogram is booked in a mass ledger whose
closure stays at round-off.  Carrier, mesh and diffusivity are those of
case 011.

## Differences from case 011

The mock thermochemistry is replaced by `constant/thermochemistryProperties`
(new):

    model temperature;
    pairs { PbI2_g { condensed PbI2_s; } }
    interface { patches (WALL); layer faceCells; ... }
    temperatureCoeffs { A 220; k 8e-3; Tdep 680; }      // 1/s, 1/K, K
    samples { boat { pair PbI2_g; selection { ... cylinder x 0.04-0.05 m,
              radius 5e-4 ... } mode release; amount 3.16e-5; duration 6; } }
    balance { nGuard 3; guardTolerance 1e-12; moleFractionWarning 0.05; }

The wall flux is `mdot'' = -rho v_d [1 - exp(-k (Tdep - T))]^+ Y`
(irreversible), with `v_d = A V_I/A_I`.

What follows from it:
- `PbI2_s` gets a molar mass and becomes derived output (the deposit of
  every wall face is the state); the passive `tracer` is dropped.
- `system/controlDict`: dt 1 ms (Courant 0.12), outer loop `nCorr 10`,
  `tolerance 1e-5`, `writeFrozenFields no` (the carrier is read, never
  written).
- `system/fvSolution`: loose GAMG correctors and one tight final solve
  (`Y_PbI2_gFinal`, 1e-12, relTol 0).

## Modified files

    constant/thermochemistryProperties   (new)
    constant/speciesTransportProperties
    0/Y_PbI2_g, 0/Y_PbI2_s               (0/Y_tracer removed)
    system/controlDict, system/fvSchemes, system/fvSolution

## Run

    ./Allrun                 # 5 steps into work/excerpt, compared with
                             # out_excerpt (about 30 s)
    ./Allrun steps 100       # the first 100 steps
    ./Allrun full 16         # the 16 s case on 16 ranks (about 2.2-2.5 s
                             # per step serial)
    ./Allclean

The mesh is converted once from run/001's `pipe.msh` into `work/pipeMesh`
(or set `LESTO_PIPE_MESH`); the carrier is case 007 at time 261 (links in
`0/`).

## What to check

At start-up:

    Interface (layer faceCells): patches (WALL), 15264 elements ...
      A_I/V_I = 17681.9 1/m
    Temperature model: v_d = 0.0124421 m/s (from A = 220 1/s on this
      layer) ... applied rate ... = 220 1/s ...
    Sample boat (pair PbI2_g): 144 cells ... release of 3.16e-05 mol
      from 0 s over 6 s
    Carrier: max |T - p W/(R rho)|/T = 0.0030285

The last line is the consistency of the 007 carrier (28 mL/min, a
development carrier; plan section 8, risk 10).

Every step ends with the balance line, whose closure stays at round-off:

    Balance PbI2_g [mol]: gas 2.63333e-08 wall 5.72124e-37 ... released
      2.63333e-08 ... closure -3.58857e-24 (relative -1.36275e-16)

(step 5, `out_excerpt`; the default `./Allrun` reproduces its 123 lines
exactly).  Over the 6 s release window (6100 steps) the closure stays at
1.8e-14 and the released amount equals n0 M to 2e-16 (criterion M2.2g).

Axial profiles and T_dep are not part of this case (so that the log keeps
reproducing `out_excerpt`); they enter the ladder in case 014.  How to add
them, the outputs and the excerpt's history: notes.md.
