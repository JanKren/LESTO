# Damköhler sweep of the half-cell HKS channel: summary records

Small text records of the sweep reported in the manuscript subsection
"Transport regime" (Fig. `damkohler_sweep`).  The runs and their fields
(about 110 MB) stay in the git-ignored
`run/016-hks-channel-halfcell/work/damkohler`.

## What was run

Case `run/016-hks-channel-halfcell` unchanged except the accommodation
coefficient α of the HKS law: plane channel, L = 0.1 m, full height
4.8 mm, 100 × 24 cells (1 mm × 0.2 mm), both walls condensing, frozen
plug-flow helium carrier (G = 0.01743 kg m⁻² s⁻¹, wall and gas temperature
ramp 1000 → 350 K), PbI2 released in a box at mid-height (x = 10–20 mm)
over 0–2 s, run to 3 s, dt = 1 ms, `wallResistance halfCell`, upwind
advection (the case's own scheme), constant D = 1e-4 m²/s.

- α = 10^0, 10^-0.5, …, 10^-6 (13 runs).  Both `accommodation` entries
  (the PbI2_g pair and `HKSCoeffs`) are changed; `Ce` and `kineticScale`
  stay 1.
- Wall-normal refinement: 100 × 48 cells at α = 1, 0.1, …, 1e-4 (5 runs).
  Deposited fraction, centroid, 10–90 % width and J^ change by at most
  0.11 %; the peak bin does not move.
- α = 1 reproduces `run/016-hks-channel-halfcell/expect.json`: wall
  inventory 2.7849908018699395e-11 mol, onset 561.1593 K, peak bin
  541.75 K.  Pair and element closures stay at 3.9e-16 (relative) in
  every run.
- Solver and build: `provenance.json` (the executable of the production,
  follow-up, Graetz and Lobresco runs).

## Damköhler number

    Da = C_e s (2α/(2 − α)) sqrt(R T_w / (2π M)) h / D

with C_e = s = 1, M = 0.46100894 kg/mol, D = 1e-4 m²/s and h = 2.4 mm,
the channel half-height (the wall-to-centreline diffusion length, the
analogue of the tube radius in the manuscript equation `damkohler`).
The factor before h/D is the solver's kinetic conductance G_e β_c/ρ_c.
Panel (a) uses the local wall temperature of each axial column; the
table, panels (c) and (d) use 561.16 K, the 1 % onset of the α = 1 run
(the onset is 561.1–561.2 K for every α).  At the production values
(400–1000 K) the channel gives Da = 1.6e3–2.6e3 (α = 1), 86–135 (α = 0.1)
and 8.2–12.9 (α = 0.01), the tube values of the manuscript.

## Quantities

- deposited, outlet and gas fractions: ledger at 3 s / released amount;
- normalised wall flux J^ = J h / [D (c̄ − c_eq)] at t = 1.5 s (during the
  release) per axial column: J from the exchange field, c̄ the
  cross-section mean, c_eq from the solver's own half-cell law
  (c_w = c_cell − J dy/(2D), c_eq = c_w − J/G_k); the table interpolation
  agrees with this law-implied c_eq to 0.26 %.  Kinetic limit J^ = Da;
  developed plug-flow diffusion limit between absorbing plates
  J^ = π²/4 = 2.467; series resistances 1/J^ = 1/Da + 4/π²;
- η = (c_w − c_eq)/(c̄ − c_eq): 1 for kinetic control, → 0 for diffusion
  control;
- 1 % onset (from x = 20 mm), peak bin, centroid, σ and 10–90 % width of
  the wall profile at 3 s; `wall_profile_L1_vs_alpha1` =
  Σ|w − w_α=1| / Σ w_α=1.

Near the onset (72.5 mm, c_eq/c̄ = 0.1–0.2, developing profile) the plateau
of J^ is higher (2.8–2.9), as expected for an entry region.  For Da < 0.05
the peak bin is the outlet bin, because the deposit runs off the end of
the channel; the centroid is the smooth measure.

## Files

| file | content |
|---|---|
| `damkohler_sweep.csv` | the table of the 13 runs |
| `sweep_summary.json` | per run (sweep and 100 × 48 refinement): the scalars above, the wall profile at 3 s, per-column J, c̄, c_w, c_eq, J^, η and local Da at 1.5 s, and the wall-normal profiles at 72.5, 84.5 and 94.5 mm |
| `provenance.json` | solver hash, OpenFOAM build, HEAD, the α = 1 check |
| `collect.py`, `foamio.py` | write the csv and json from the runs |
| `figure.py` | writes `manuscript/figures/damkohler_sweep.{pdf,png}` from `sweep_summary.json` |

`manuscript/prepare_results.py` reads `sweep_summary.json` for the numbers
of the text, and `manuscript/check_results.py` recomputes them from the
csv.

## Regenerate

With OpenFOAM v2412 sourced and the solver built (or `LESTO_SOLVER` set):

    run/016-hks-channel-halfcell/Allrun-damkohler 24    # about 2 min
    run/016-hks-channel-halfcell/Allrun-damkohler 48    # about 1 min
    python3 manuscript/data/damkohler/collect.py \
        run/016-hks-channel-halfcell/work/damkohler/ny24 \
        run/016-hks-channel-halfcell/work/damkohler/ny48
    python3 manuscript/data/damkohler/figure.py

On 15 October 2026 this sequence reproduced every metric of the first
(scratch) sweep of 14 October exactly; only the execution times differ.

## Scope

A plane channel with a plug-flow carrier, a linear temperature ramp and a
constant diffusivity, not the tube cases.  It shows that the discretised
wall law reproduces both limits of the exchange and the transition
between them, and where the production accommodation coefficients sit on
that curve.  It does not vary the interfacial conductance on the Liu et
al. or Lobresco et al. cases.
