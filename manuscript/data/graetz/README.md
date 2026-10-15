# Graetz-Robin verification of the wall exchange (M7.1): summary records

Small text records of the verification reported in the manuscript
subsection "Verification of the wall-exchange discretisation" (Table III of
the arXiv build). The runs, their meshes and fields (1.2 GB) stay in the
git-ignored test work area
`applications/rhoFixedFlowFoam/tests/work/v2412-linux64GccDPInt32Opt/verification/graetz`.

## What was run

`applications/rhoFixedFlowFoam/tests/verification/graetzRobin` (criteria
M7.1a-o of `applications/rhoFixedFlowFoam/doc/phase-change-plan.md`, M7
acceptance 1; description in `applications/rhoFixedFlowFoam/readme.md`,
"Verification (M7)"). Extended Graetz problem: fully developed laminar
pipe flow with axial diffusion, constant first-order wall uptake
(temperature law with f = 1, `depositionVelocity` = k_w), R = 2.4 mm,
T = 500 K, p = 1 atm, D = 7.6963e-5 m2/s, Pe = U_b R/D = 5,
Bi = k_w R/D = 1 and 50, L = 24 R. Wedges of 2 degrees with exactly
axisymmetric metrics, Nr = 32, 64, 128, 256 radial cells, dx = 2 dr,
exact parabolic face fluxes; Gauss linear (zeroGradient and fixedValue
outlets) and SuperBee advection; wall laws `wallResistance halfCell` and
`none`; steady states by implicit Euler at dt = 1e6 s; one time-accurate
control run. The decay rate mu is fitted to the gas line density in a
window fixed by the spectrum. Reference: the smallest positive eigenvalue
by power series (mpmath), DOP853 shooting and Chebyshev collocation
(`graetz.py`), relative error bound 1.69e-12 (Bi 1) and 4.5e-14 (Bi 50).

## Identity of the runs

- Runs made on 11 October 2026 with OpenFOAM v2412
  (`_8dbc61e11c-20241220`).
- Solver `rhoFixedFlowFoam`, SHA-256
  `16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58`,
  the executable of the production, follow-up and Lobresco runs. Evidence:
  every run key (first field: md5 of the solver and helper executables,
  the case template and `graetz.py`/`graetzCase.py`) equals the key
  recomputed from that executable and the current scripts
  (`provenance.json`, `runs_made_with_this_solver: true`).
- Solver sources last changed in commit `fe214c6`; test scripts
  (`graetzRobin`, `graetz*.py`, `graetzWedge/`, `graetzCarrier/`)
  unchanged since commit `6311edf`; records collected at
  repository HEAD `3cea399` (14 October 2026).
- On 14 October 2026 `graetzRobin check` (evaluation of the existing runs,
  no new runs) gave 27 PASS (`criteria.log`), and reproduced `table.md`,
  `log.fit` and `reference.json` byte for byte.

## Files

| file | content |
|---|---|
| `reference.json` | reference decay rates m = mu R, the three methods, error bound, fit window |
| `fits.json` | per run: set-up, fitted mu, mu_h of the discrete mode, window variants, mixing-cup and wall-uptake fits, steady-state measure |
| `table.md` | the summary table of `graetzChecks.py` (every run; GCI per triplet) |
| `criteria.log` | the PASS/FAIL lines of `graetzRobin check` |
| `provenance.json` | solver hash, OpenFOAM build, run date, key check, HEAD |
| `collect.py` | writes the files above from the work area |

`manuscript/prepare_results.py` builds the manuscript table from
`fits.json` and `reference.json` (GCI of Celik et al. 2008, ratio 2,
F_s = 1.25), and `manuscript/check_results.py` recomputes it and checks
it against `table.md` and the file hashes.

## Regenerate

With OpenFOAM v2412 sourced and the solver built
(`applications/rhoFixedFlowFoam/Allwmake`), python3 with numpy, scipy and
mpmath:

    cd applications/rhoFixedFlowFoam
    tests/verification/graetzRobin          # runs (about 20 min on 4 cores), then the criteria
    tests/verification/graetzRobin check > graetz-check.log   # criteria of existing runs only
    cd ../..
    python3 manuscript/data/graetz/collect.py \
        applications/rhoFixedFlowFoam/tests/work/v2412-linux64GccDPInt32Opt/verification/graetz \
        applications/rhoFixedFlowFoam/graetz-check.log

## Scope

These records verify the discretisation of a linear (Robin) wall exchange
law and its coupling to advection and diffusion at Bi = 1 and 50. They do
not verify the nonlinear HKS law with tabulated equilibrium pressures, the
gas equilibrium, or the spatial convergence of the full iodine
calculations.
