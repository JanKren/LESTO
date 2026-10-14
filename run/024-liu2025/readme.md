# Case 024 — Liu et al. (2025): Pb-Bi-I deposition on physical carriers

This case is based on case 021, with the carrier technique of case 017.

Liu et al. (2025) measured the iodine deposits of PbI2 and of iodine
released from LBE in thermochromatography columns of steel and silica.
This study computes their three configurations with the Pb-Bi-I chemistry
of case 021 (gas re-speciation, PbI2 and BiI3 pairs, metal wall channels,
the Q2 metal-vapour source) on physical helium carriers of the 1.25 m
columns (`rhoSimpleFoam`, the published wall temperatures, flux projection
and audit, as case 017), with the half-cell wall law of case 016.  It is
provisional numerical work: the release history, the flow reference
conditions and several coordinates are assumptions (notes.md).

## Differences from case 021

`thermochemistry/benchmarks/liu2025/cases.py` writes every case from case
021's `speciesTransportProperties`, `fvSchemes`, `fvSolution` and
`thermochemistryProperties`, and changes:

- the carrier: a 5-degree wedge of the column (steel radius 2.4 mm, silica
  2.5 mm; 384 x 48 = 18,432 cells), the measured wall temperature, a
  mass-flow inlet; sources and flow both use the 1/72 wedge fraction;
- `interface { ... interfaceTemperature wall; wallResistance halfCell; }`;
- the sources: Q1 (the deposited amount over three hours) and Q2 (LBE
  metal vapour at 0.1, 1 and 10 % saturation, two BiI datasets);
- dt 2 ms, 60 s.

## Run

    export LESTO_M10E_NASA=/path/to/nasa_thermo.inp
    export LESTO_M10E_THERMOFUN=/path/to/bi_i_m10-thermofun.json
    export LESTO_M10E_FLOW_T=298.15      # assumed until the authors confirm
    export LESTO_M10E_FLOW_P=101325
    ./Allrun            # smoke: 128 x 16 cells, 1 s at 20 ms (minutes)
    ./Allrun coarse     # 128 x 16, 60 s
    ./Allrun full       # 384 x 48, 60 s at 2 ms (about 6 h per case serial;
                        # LESTO_M10E_JOBS cases at a time)
    ./Allrun followup status --out <study work>   # the follow-up study
    ./Allclean          # work/smoke and work/coarse only

`work/` holds the raw output of every study (33 GB: `prepared/`, the 15
full-size runs, the only raw CFD evidence of M10e; `followup-20261007/`;
the attempts and probes).  `./Allclean --all` removes it only after a
confirmation, and never while a `run_state.json` reports a live RUNNING
process.  Licensed thermodynamic inputs and generated tables stay outside
git.  The follow-up package of 2026-10-12 (`~/lesto-submission-followup`)
names this case by its former number (034) in its `study.json`: resume it
with `--prepared <LESTO>/run/024-liu2025/work/prepared` and `--carrier-tool
<LESTO>/run/024-liu2025/work/followup-20261007/bin/liuCarrierTool`.

## What to check

- Every case: element closure <= 1e-12 and solverDefect <= 1e-8 of the
  element reference in every recorded step (recorded at most 8.63e-16),
  zero equilibrium failures (`full_size_completion.json`).
- Production, 576 x 72 = 41,472 cells, `Gauss limitedLinear 1`
  (`production_576x72_limitedLinear/`): 19 PASS; pure PbI2 binned peak at
  36 cm / 366.5 C against the published 36 cm / 368 +- 18 C; steel PbI2
  288.1 C binned.
- Scheme and mesh study (`scheme_mesh_study/`): the iodine profile changes
  by 1.98 % between 41,472 and 93,312 cells (limitedLinear, steel); no
  formal GCI (inconsistent apparent orders).
- Iteration check (`iteration_check/`): eight against three correctors,
  iodine profile L1 at most 0.021 %.
- `python3 manuscript/check_results.py` recomputes these numbers from the
  saved profiles: PASS.

## Where the results are

| folder or file | what |
|---|---|
| `full_analysis/` | the 15-case matrix on 18,432 cells, upwind: report, figures, CSVs, profiles |
| `followup_results/`, `followup_launch.json` | the 64-case follow-up (refinement and provisional sensitivities) |
| `followup_interim_20261008/` | the interim report of the follow-up |
| `production_576x72_limitedLinear/` | the 19 production runs (manuscript) |
| `scheme_mesh_study/`, `iteration_check/` | the numerical studies of plan section 39 |
| `physical_input_audit/` | evidence and gaps of the physical inputs |
| `profiles/`, `preliminary_*` | the preliminary (2026-10-06) profiles |
| `verification.json`, `full_size_*.json` | the run records (absolute paths of their date) |

Details: notes.md (the run modes, numerical checks and physical
assumptions, what remains, and the 2026-10-06/07 execution history).
