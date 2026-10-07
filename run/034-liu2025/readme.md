# M10e: Liu (2025) physical carriers and provisional validation

**Latest status (2026-10-07):** the available full-size matrix is complete:
15 PASS, zero FAIL, one missing-profile case NOTRUN. Mesh/time-step and
experimental qualification remain open. See the
[full analysis](full_analysis/report.html) and
[completion evidence](full_size_completion.json).

This case prepares physical laminar helium carriers with `rhoSimpleFoam`,
then runs Q1 deposition and Q2 metal-vapour sensitivities using the M10d
kernel and shared wall condensates. The full mesh is 384 × 48 × 1 (18,432
cells), a 5° wedge of the 1.25 m column. Steel radius is 2.4 mm; silica
radius is 2.5 mm. Sources and mass flow both use the 1/72 wedge fraction.

This is **preliminary numerical preparation, not completed M10e experimental
validation**. Three published temperature curves and iodine histograms are
digitised in [the benchmark package](../../thermochemistry/benchmarks/liu2025/README.md).
The paper does not plot `BiI3_SiO2_800` temperature data. Figure 7C's curve
and printed peak temperatures disagree by approximately 63 K (PbI2) and
108 K (BiI3); the trace is retained without fitting it to the labels.

## Run

Load OpenFOAM v2412 and build the current `rhoFixedFlowFoam`. Set external
thermodynamic paths and **explicit** flow reference conditions:

```bash
export LESTO_M10E_NASA=/path/to/nasa_thermo.inp
export LESTO_M10E_THERMOFUN=/path/to/bi_i_m10-thermofun.json
export LESTO_M10E_FLOW_T=298.15
export LESTO_M10E_FLOW_P=101325
# These reference conditions are an assumption until confirmed by the authors.
bash run/034-liu2025/Allrun smoke
```

Each mode uses a separate fresh `work/<mode>` directory. Existing case
directories are refused, preserving previous results. `smoke` runs 1 s at
20 ms on 128 × 16 cells; `coarse` runs 60 s on the same coarse mesh;
`prepare` prepares full-size carriers and all available 60 s / 2 ms cases;
`full` also runs those cases. Use `LESTO_M10E_JOBS` for concurrent independent
serial cases. Full transport runs can take hours per case. Licensed
thermodynamic inputs and generated thermodynamic tables stay outside git.

Sixteen cases are listed: four Q1 cases, plus two LBE cases × two BiI
variants × three metal saturation fractions (0.001, 0.01, 0.1). Fifteen can
currently run; the missing BiI3 profile remains `NOTRUN`. The +29 kJ/mol
BiI variant is a sensitivity surrogate, not a second qualified dataset.

## Numerical checks and physical assumptions

The carrier tool initialises T and the fixed wall-face values from the
published profile. It retains a mass-flow inlet, ideal-gas helium,
Sutherland transport and constant Cp from case 007. The pressure solve
has a small finite-precision flux residue at these low Mach numbers.
A zero-based finite-volume flux-potential solve projects that residue away:
it refuses original cell/global defects above 1e-6 of inlet mass flow and
flux/velocity corrections above 1e-5. The original defects and correction
sizes are recorded. The subsequent carrier audit requires inlet flow,
global balance, cell balance and wall/axis leakage within 1e-12 of inflow.
Wedge boundary fluxes cancel and are checked as a net flux.

The transient report independently requires **element closure ≤1e-12**
and **solverDefect ≤1e-8** of the element reference at every recorded step,
and refuses incomplete/nonfinite outputs. Peak and iodine comparisons
are provisional observations; experimental uncertainty and the proposed
20/40 K science targets are reported separately. Wall deposits are binned
in 1 cm bins centred on the gamma scan coordinates. The BiI3 split uses
iodine moles, `3 n(BiI3)/(2 n(PbI2)+3 n(BiI3))`. Pb/Bi metal deposits do
not enter that iodide split. KI's published 6.6/7% share remains excluded
and the iodine profile comparison retains that exclusion explicitly.

For runs with at least three profile writes, the report compares integrated
net deposition rates and their spatial L1 difference in the last two windows.
Changes below 1% are labelled preliminary stationarity. This diagnostic does
not replace mesh/time-step refinement or experimental input qualification.
The peak table also reports the shift from the published constant-H/S
Eq. (3) threshold; that historical threshold is not a phase-resolved
equilibrium onset for the evolving Q2 gas.

For Q1, the prescribed rate is the three-hour sample/deposit amount divided
by 10800 s. A shorter CFD window releases only its corresponding amount;
it never releases the whole sample in 60 s. LBE Q1 conditions the source
on the measured deposited split and cannot validate a prediction of that
split. Q2 adds ideal-activity LBE metal vapour and re-speciates an atomic
Pb/Bi/I source, still conditioning iodine release on the deposited amount.
Both retain the source until the simulation ends; no clean-helium flush
is appended. The experiment stopped the flow and cooled before scanning.
Cooling/shutdown and its effect on the deposits are not yet modelled.

Unconfirmed assumptions include plot-to-mesh origin, constant temperature
extension over the unplotted ends, source placement at 9–10.5 cm in the
interior radial cells, omitted boat obstruction, constant gas diffusivity
1e-4 m²/s, HKS accommodation unity, and the release history. The flow
reference of 298.15 K / 101325 Pa approximately reproduces the paper's
Eq. (3), but that inference is not experimental confirmation.

## What remains for completion

Obtain the numerical wall-temperature profiles and coordinate origins,
especially `BiI3_SiO2_800` and the Figure 7C clarification; raw 1 cm gamma
scans; flow-controller reference conditions; boat placement/geometry and
release histories. Jörg Neuhausen and Lu Liu are the paper's listed contacts.
No request has been sent. Qualify diffusivity and thermodynamic provenance,
check quasi-steady spatial profiles, perform mesh/time-step refinement,
and run the complete full-size matrix with both qualified BiI datasets.
`verification.json` records completed preliminary work and explicit gaps.

## Preliminary results (2026-10-06)

The supplied paper is the only experimental input available. Three full-size
and three coarse physical carriers pass the 1e-12 mass-flow audits after the
bounded flux projection. Fifteen coarse cases ran for 60 s at 20 ms; all pass
the numerical conservation gates, with maximum relative element closure
6.11e-16 and solver defect 1.39e-13. Homogeneous re-speciation failures are
zero. Fourteen cases meet the selected 1% stationarity diagnostic; the silica
Q1 case does not, because its extremely small Pb metal deposit remains
sensitive to numerical rate differences.

Five selected full-size cases pass 1 s smoke checks. Two selected coarse
cases also pass at 10 ms for 60 s. Their iodide peak temperatures are unchanged,
but the deposit-profile L1 changes are 0.0243% for pure PbI2 and 3.17% / 2.00%
for PbI2 / BiI3 in the steel Q2 baseline case at 1% metal saturation. The
full-size 60 s / 2 ms matrix was subsequently launched locally; its completion
and conservation results are recorded separately from these coarse results.

Available Q1 peak errors are -18.43 K for pure PbI2; -17.37 K / +11.56 K
for steel PbI2 / BiI3; and -17.41 K / +21.73 K for silica PbI2 / BiI3.
These fall inside the proposed 20/40 K science targets, but pure PbI2 narrowly
exceeds its printed experimental uncertainty. Iodine-profile L1 differences
span 47.6–185.1 percentage points over the digitised scan range across the
coarse matrix. Peak agreement therefore does not establish profile validation.
The Q2 iodide split varies strongly with the BiI sensitivity and saturation;
neither variant has been fitted to the experimental profiles.

The full silica smoke exposed a reactive-pressure underflow: dividing a
positive subnormal product pressure by 1 bar rounded it to zero before
`log()`, causing an FPE. The solver now performs that conversion in log space,
and an independent closed-form test runs with floating-point traps enabled.
The final v2412 `M0 M10a M10b M10c M10d` regression records **43 PASS,
0 FAIL, 13 NOTRUN**. The omitted external studies and separate v2606 check
are listed in the evidence. Eleven Python input/failure-path tests pass.

Saved outputs:

- [Peak temperatures, targets and explicit NOTRUN rows](preliminary_peaks.csv)
- [Iodine profile comparison plot](preliminary_profiles.png)
- [Per-case iodine comparisons](profiles/)
- [Numerical results, carrier audits, assumptions and source hashes](verification.json)

The full 3-hour experimental inventories, cooling stage, qualified transport
properties, missing temperature profile and second thermodynamic dataset
remain outside this preliminary qualification.

## Full-size matrix execution (2026-10-06)

All fifteen available cases were launched concurrently at 18,432 cells,
60 s end time and 2 ms time steps. The local runner continues independently
of the chat and writes `work/prepared/run_state.json` every 30 s, with each
case's status, process id and latest simulated time. On completion it runs
the numerical report and creates `work/prepared/cases/validation.json`,
`peaks.csv`, per-case iodine comparisons and `preliminary_profiles.png`.
`RUNNING` is not a completed conservation or scientific validation result.
The missing `BiI3_SiO2_800` case remains explicitly NOTRUN.

The first silica Q1 step exposed a second precision issue: a bulk Pb
concentration of 3.0847e-310 mol/m3 lost relative accuracy when the equilibrium
kernel evaluated its smaller pressure with double libm. The kernel now
bypasses the double predictor for extremely dilute bulks and uses extended
exponentials for extreme arguments. It keeps every species in the equations
and retains the existing 1e-14 element-residue gate. Twenty-four independent
monatomic-limit tests cover the smallest representable double with FPE traps,
and the original 140-point equilibrium reference remains unchanged.
The corrected silica 10-step smoke passes with maximum element closure
5.32e-16. All cases restarted from fresh initial fields with a private copy
of that corrected binary; the first attempt is retained under
`work/attempt-1-pressure-underflow`.

The v2412 regression session recorded 22 PASS, zero FAIL and nine external
studies NOTRUN before receiving SIGTERM in the last scientific half-step run.
That remaining M10d.3 check was completed separately: both thermodynamic
variants and both time steps retain all five onsets within the 10 K gate and
all ledgers within 1e-14. Thus 23 numerical/build regression checks completed
successfully. The session did not reach its SELF/CACHE/LIST framework checks;
this is not claimed as a completed aggregate `Alltest` invocation. The
default and optional solver builds have zero warnings; no new v2606 check
was performed.

[The launch record](full_size_launch.json) captures the configuration,
source/binary hashes and regression evidence. It is a timestamped launch
snapshot; consult the local `run_state.json` for live completion status.

Six Q2 cases in that second attempt stopped at 1.8–2.6 s with equilibrium
nonconvergence. Near a strongly bound iodide, subtraction in the iodine
Schur derivative could make a positive derivative negative; Newton steps
could also cycle between opposite sides of a root. The kernel now evaluates
that derivative as nonnegative pairwise weighted variances, bisects finite
brackets when steps fail to contract, and treats the double predictor as an
accelerator with a cold extended-precision fallback. The element gate remains
1e-14 and no species or positive bulk is removed.

The revised default v2412 build has zero warnings. The unchanged 140-point
Decimal reference and 24 tiny-bulk checks pass. A new 54-case synthetic
strong-binding regression passes with FPE traps; its first case fails with
the previous kernel. A local 600,000-input stress check using the six failed
cells' thermodynamic constants has zero convergence or element-balance
failures (maximum relative residue 5.03e-16); an additional 30,000 synthetic
affinity cases converge. Closed-box checks, serial/four-rank agreement,
4,800 written-field equilibrium points, exact binary restarts and analytical
shared-wall checks pass. Two full-size silica cases pass ten steps with
zero re-speciation failures and maximum relative element closure 6.39e-16.
These are targeted checks, not a new complete `Alltest` invocation.

The second attempt is preserved under
`work/attempt-2-equilibrium-nonconvergence`. All fifteen cases were launched
fresh for attempt three with one corrected private binary and identical
thermodynamic tables, meshes, sources, end time and time step. Field/profile
writes now occur every 1 s (500 steps), with the latest two field checkpoints
retained; the case manifests record these execution overrides. The complete
60 s numerical and experimental comparisons remain pending.

## Completed full-size matrix and analysis (2026-10-07)

Attempt three finished at 00:28 UTC (02:28 Zurich). All fifteen available
cases reached 60 s, with 18,432 cells and 2 ms steps. The maximum accounted
element-ledger closure is 8.63e-16 relative to its reference; independently,
the raw material residue before numerical corrections is at most 5.98e-13
relative to supplied atoms. Gas-equilibrium solver failures are zero.

All five available Q1 printed peaks lie within the paper's uncertainties:
350.73 °C for pure PbI2; 288.09 / 166.41 °C for steel PbI2 / BiI3;
287.59 / 152.77 °C for silica PbI2 / BiI3. Q1 LBE sources prescribe an
observed-deposit split, so this is conditional evidence. Silica's Q1 peaks
remain 5 / 4 cm downstream of the observed scan maxima, and its normalised
profile overlap is only about 7%. Figure 7C's +63 / +108 K curve/label
discrepancy remains unresolved; no coordinates or temperatures were fitted.

The low-loading steel Q2 baseline has the lowest raw iodine-profile error
among the tested scenarios (36.80 percentage points L1), but predicts a
56.47% BiI3 share of iodide iodine versus the paper's 65.63%. At 1% / 10%
metal saturation, the baseline gives roughly 94.5% BiI3 in both columns;
the +29 kJ/mol BiI surrogate gives roughly 2%. Source and thermodynamic
assumptions therefore remain consequential and unqualified.

The largest coarse-to-full iodide profile change is 67.95% L1. Mesh,
time step and solver revision changed together, so this comparison does
not establish independent convergence. First-order upwind advection has
estimated numerical diffusion comparable to or greater than the imposed
1e-4 m²/s physical diffusivity even on the full mesh. Isolated time-step
and fixed-step mesh comparisons are the next numerical priority.

Most paired-gas steps miss the strict 1e-16 outer criterion, with median
residuals near 1–2e-15 and no significant regime changes. Twelve short
matched restarts across four cases compare nCorr 3 / 8 and outer tolerances
1e-16 / 1e-14. The largest species-scaled gas-field change is 3.60e-8;
new iodide/iodine wall-profile increments differ by at most 9.38e-14 L1.
Three five-step startup replays reproduce transient invalid-element skips
and bound iodine in skipped cells below 4e-225 of supply. These probes
support a round-off interpretation without changing gates or adding floors;
they do not replace full-run refinement.

Relevant deposition-rate profiles are steady: comparing 20–30 and 50–60 s
windows changes iodide/iodine rates by less than 0.001%. Cumulative profiles
retain startup differences, and the experiment's shutdown/cooling stage is
not represented. Fifteen Python benchmark/analysis tests pass.

The report includes all case metrics, source-conditioned splits, independent
material audits, transient diagnostics, refinement limits, input discrepancies
and prioritised next work. Its figures are embedded for standalone viewing;
[a six-page figure PDF](full_analysis/analysis_figures.pdf), CSV tables,
machine-readable results and final numerical profiles are saved beside it.
Licensed thermodynamic inputs and raw CFD fields/logs remain outside git.

Regenerate the analysis from retained local CFD outputs:

```bash
MPLCONFIGDIR=/tmp/lesto-matplotlib python3 thermochemistry/benchmarks/liu2025/analyse_matrix.py \
  --full-root run/034-liu2025/work/prepared/cases \
  --coarse-root /tmp/lesto-m10e-final-coarse \
  --half-root /tmp/lesto-m10e-final-half-step \
  --corrector-summary run/034-liu2025/work/corrector-probe/comparison.json \
  --positivity-summary run/034-liu2025/work/positivity-probe/comparison.json \
  --out run/034-liu2025/full_analysis
```

`corrector_probe.py` and `positivity_probe.py` reproduce the diagnostics
using the unchanged private solver in an OpenFOAM v2412 environment; each
requires a fresh `--out` directory. Both preserve the original sources.
