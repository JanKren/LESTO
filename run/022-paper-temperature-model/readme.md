# Case 022 — Code-to-code benchmark of the temperature model (paper mode)

This case is based on case 012 (not 021: it is the first of the three
benchmark studies at the top of the ladder, and it needs only case 012's
temperature model).

Lobresco et al. (2026), Nucl. Eng. Des. 459, 115201, computed their
temperature-based deposition model (Section 2.3.1, Eq. 5) with the T-Flows
code, in a "paper mode": a frozen, non-expanding parabolic helium flow of
density 1, the temperature of that flow with the measured wall table, and
the release of the sample over the injection interval.  This case puts
case 012 into exactly that mode, with every value of the T-Flows parity
checklist pinned (plan section 11; the pins and their sources: notes.md),
and compares with Figs. 8-10 and Table 2 of the paper.

## Differences from case 012

- The carrier is written by `setFrozenCarrier` (`system/setFrozenCarrierDict`,
  new): density 1, the exact flux of the parabola at 1.67e-6 m3/s
  (100 mL/min), uniform p, and T solved in that flow with the T-Flows wall
  table `constant/wall_x_profile_complete.dat` (new); no 007 links in `0/`.
- `interface { layer distance; distance 1e-4; }`: the T-Flows layer of
  the cell centres within 0.1 mm of the wall (DELTA), 30528 cells,
  `A_I/V_I = 8533.7 1/m`.
- The sample releases 3.16e-5 mol over t_ev = 15 s; end time 30 s (2 t_ev).
- The inlet is closed to the species, as in T-Flows (zero total flux,
  Danckwerts: `codedMixed` in `0/Y_PbI2_g`); `profiles` (424 native and
  1 cm bins, `observable wall`); `system/decomposeParDict` (16 slabs).

## Modified files

    constant/thermochemistryProperties, constant/speciesTransportProperties
    constant/wall_x_profile_complete.dat, system/setFrozenCarrierDict (new)
    0/Y_PbI2_g, 0/Y_PbI2_s           (the carrier fields are written)
    system/controlDict, system/fvSchemes, system/fvSolution,
    system/decomposeParDict

## Run

    ./Allrun            # 5 steps into work/excerpt, compared with
                        # out_excerpt (about 1 min)
    ./Allrun full 16    # the 30 s case on 16 ranks (2-3 h); the carrier is
                        # written on the decomposed case
    ./Allrun long       # the runs of acceptance 3-5: tests/Alltest M4long
    ./Allclean

`utilities/axialProfiles.py <case> --set gamma` prints T_dep and the
integrals of a finished run.

## What to check

At start-up: `Interface (layer distance) ... 30528 elements ... A_I/V_I =
8533.68 1/m` and the T-Flows centre-distance rule selecting the same cells.
The default `./Allrun` reproduces `out_excerpt` (126 lines, 5 steps).

The full run against the paper (M4long run `fig9-tev15`, criterion M4.16):

| quantity | paper | case 022 | difference |
|---|---|---|---|
| gas plateau at 16 s (Fig. 9) | 2.26e-5 mol/m | 2.2712e-5 mol/m | +0.5 % |
| gas at 50 % of the plateau, 16 s | 0.507 m | 0.4987 m | -8.3 mm |
| solid onset (1 %), 30 s (Fig. 10) | 0.502 m | 0.4761 m (T_dep 679.82 K) | -25.9 mm |
| solid peak, 30 s | 1.14e-3 mol/m at 0.513 m | 9.11e-4 mol/m at 0.4955 m | -17.5 mm |

The plateau (n0/U_b) agrees.  The deposit starts where the wall table
crosses the input Tdep = 680 K, 2.6 cm upstream of the published onset;
the T-Flows parameter set (Tdep 638.15 K, A 150, k 0.02, distance 2.4e-4)
puts it at 0.5011 m (-0.9 mm).  Which set produced Figs. 8-10 is a question
for the authors (plan section 26).

Table 2 (M4long, both layer rules, five experiments): all ten lines meet
|T_dep - T_paper| <= max(2 K, dT_bin).

Details: notes.md (sources and pins, outputs and T_dep, the decomposed
carrier, the runs of acceptance 3-5, T5_60, the tests M4.2, M4.7c, M4.11,
M4.12 and M4.15 that use this case).
