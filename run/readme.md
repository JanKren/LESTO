# The case ladder

Each case adds one thing to the case it is based on; its readme starts
with "This case is based on case N", shows the changed dictionary lines,
and ends with the numbers to check.  Cases 001-011 (Bojan Ničeno) build
the pipe, the helium carrier and the species transport; cases 012-021 add
the gas/condensate exchange step by step; cases 022-024 are the three
benchmark studies.  New cases go on top (025 and later).

| No. | case | based on | the one change | default check, cost |
|---|---|---|---|---|
| 001 | first-pipe-flow-simpleFoam | - | laminar pipe, gmsh mesh, `simpleFoam` | in place |
| 002 | parabolic-inlet | 001 | parabolic inlet profile (coded) | in place |
| 003 | parabolic-initialization | 002 | parabolic initial field (`setExprFields`) | in place |
| 004 | energy-with-rhoSimpleFoam | 003 | energy equation, `rhoSimpleFoam` | in place |
| 005 | varying-wall-temperature | 004 | linear wall temperature | in place |
| 006 | variable-density | 005 | perfect-gas density | in place |
| 007 | variable-properties-high-temperature | 006 | Sutherland transport at high temperature; the carrier at 261 | in place |
| 008 | one-species-transport-rhoFixedFlowFoam | 007 | `rhoFixedFlowFoam`: gaseous species on the frozen carrier | in place |
| 009 | variable-diffusivity | 008 | PbI2-He diffusivity | in place |
| 010 | handles-solid-species | 009 | solid species, mock thermochemistry | in place |
| 011 | coupled-species-sources | 010 | coupled gas-to-solid mock source | in place |
| 012 | temperature-model | 011 | wall law: temperature model (Eq. 5), release sample, mass ledger | excerpt, 30 s |
| 013 | hks-inventory | 012 | HKS law with the PbI2 p_v table; evaporating inventory boat | excerpt, 30 s |
| 014 | hks-channel | 013 | carrier: the synthetic 1000-350 K plane channel; profiles, T_dep | 3 s run, 10 s |
| 015 | hks-gems-local | 013 | p_eq from GEMS3K, mode local (optional build) | excerpt, 30 s |
| 016 | hks-channel-halfcell | 014 | `wallResistance halfCell` | 3 s run, 10 s |
| 017 | physical-carrier | 016 | the physical carrier of Lobresco T106_6 (coarse wedge) | 6 s run, 8 min |
| 018 | pbi2-bii3-channel | 014 | a second pair, BiI3 | 3 s run, 10 s |
| 019 | pbbi-gas-speciation | 018 | homogeneous Pb-Bi-I gas re-speciation (kernel) | 3 s run, 4 min |
| 020 | pbbi-gems-speciation | 019 | re-speciation by GEMS3K (optional build) | 3 s run, 5 min |
| 021 | lbe-metal-vapour | 019 | metal vapour, shared condensates, reactive wall channels | 1.2 s run, 40 s |
| 022 | paper-temperature-model | 012 | benchmark: paper mode, code-to-code with T-Flows (Figs. 8-10, Table 2) | excerpt, 1 min |
| 023 | lobresco2026 | 017 | benchmark: the five Lobresco experiments, two meshes | coarse T106_6, 8 min |
| 024 | liu2025 | 021 (+017) | benchmark: the Liu (2025) matrix and its studies | smoke, minutes |

## Running the cases

Cases 001-011 run in place, as their readmes say: for 001-003
`gmsh -3 pipe.geo; gmshToFoam pipe.msh; checkMesh`, then `simpleFoam > out`
(003: `setExprFields -time 0` first); for 004-007
`rm -rf 0; cp -r 0.orig 0; setExprFields; rhoSimpleFoam > out`; 008-011
read the mesh of 001 and the carrier of 007 at time 261 through links and
run `rhoFixedFlowFoam > out` (the test gate `tests/Alltest M0` runs copies
of them).

Cases 012-024 have `Allrun` and `Allclean` and share `RunFunctions`:

    ./Allrun            # the default check (table above), into work/
    ./Allrun <mode>     # the longer modes in the comment of each Allrun
    ./Allclean          # removes work/ (024: only work/smoke, work/coarse)

Everything is written into the case's git-ignored `work/`; no Allrun
writes a tracked file, and no Allclean removes one.  An Allrun ends with
one line `PASS`, `FAIL` or `NOTRUN` (exit status 0, 1, 0); it compares with
`out_excerpt` (012, 013, 015, 022) or with the numbers of `expect.json`
(`lestoCheck.py`).  `./Allrun-ladder 012-024` runs every default check and
prints one line per case.

Load OpenFOAM v2412 (or v2606) and build `applications/rhoFixedFlowFoam`
(`./Allwmake`) first.  Optional environment:

| variable | used by | what |
|---|---|---|
| `LESTO_SOLVER` | all | the rhoFixedFlowFoam to run (default `$FOAM_USER_APPBIN/rhoFixedFlowFoam`); for 015 and 020 a build with the GEMS3K bridge |
| `LESTO_CPUS` | all | a taskset core list, e.g. `4-11` |
| `LESTO_PIPE_MESH` | 012, 013, 015, 022 | an existing gmshToFoam conversion of 001's `pipe.msh` (default: made once in `work/pipeMesh`) |
| `LESTO_M10_DATA` | 018-020 | external `pv_BiI3_phases.csv`, `log10Kf_PbBiI.csv` |
| `LESTO_M10_SYSTEM` | 020 | external GEMS3K system `PbBiIHe-dat.lst` |
| `LESTO_M10D_DATA` | 021 | external tables of `make_wall_tables.py` |
| `LESTO_M10E_NASA`, `LESTO_M10E_THERMOFUN`, `LESTO_M10E_FLOW_T`, `LESTO_M10E_FLOW_P` | 024 | external thermodynamic inputs, flow reference conditions |

Without its external input a case prints NOTRUN and the sweep continues.
The restricted Bi thermodynamic records and the tables made from them are
not redistributed.

## Renumbering of 2026-10 (plan section 40)

| old | new |
|---|---|
| 014-paper-temperature-model | 022-paper-temperature-model |
| 030-pbi2-bii3-channel | 018-pbi2-bii3-channel |
| 031-pbbi-gas-speciation | 019-pbbi-gas-speciation |
| 032-pbbi-gems-speciation | 020-pbbi-gems-speciation |
| 033-lbe-metal-vapour | 021-lbe-metal-vapour |
| 035-lobresco2026 | 023-lobresco2026 |
| 034-liu2025 | 024-liu2025 |
| (new) | 014-hks-channel, 016-hks-channel-halfcell, 017-physical-carrier |

012, 013 and 015 kept their numbers.  The number 014 was reused: the plan's
sections 0-39 and the run records (JSON files with absolute paths) keep
the old numbers of their date.  The tests look the cases up by name
(`makeCase paper-temperature-model`).
