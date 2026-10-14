# Author review: 14 October 2026 (update of 9 October)

## Scientific reassessment on 12 October

Saved profiles now support 28 mass-preserving bin-aggregation comparisons
and four matched-resolution thermodynamic comparisons. A new verification
table shows that 4 cm grouping does not robustly remove the shifted-BiI
mesh sensitivity (0.52--10.02%, depending on bin origin). Binned transport
distances are small, but do not bound continuum error. The BiI perturbation
changes iodide partition by at least 92.555 percentage points at every
tested matched resolution. These are new analyses of existing CFD results.

The current Windows workspace and its accessible Ubuntu environment lack
the original OpenFOAM installation, private executable, prepared case
directories and external Bi thermodynamic input. No new CFD checks were launched in this workspace during that assessment.
The pulled 14 October results now include the production, scheme/mesh
and three full-duration corrector checks described below. See `scientific-gap-register.md` and
`next-study.json` for the outstanding evidence and a specified next study.

## Editorial revision on 12 October

The revised manuscript states three research questions, clarifies the
transport-source units and the limiting wall conductance, and distinguishes
the three profile comparison metrics. A new discussion connects numerical
resolution, input identifiability and the evidence needed for predictive
validation. The conclusions now state the principal findings directly.
No new simulations, experimental inputs or fitted parameters were introduced;
the scientific submission issues below remain open.

The working draft is framed as a conservative CFD method plus
numerical-resolution and input-sensitivity assessment, now with an
independent single-species comparison. It reports completed computational
cases, not independent experiments, and makes no claim of independent
predictive validation. The missing BiI3_SiO2_800 temperature profile is
stated explicitly; that case was not run.

## What changed on 14 October 2026

All new evidence is in the repository; plan section 39 summarises it and
awaits sign-off.

- **Production runs** (`run/024-liu2025/production_576x72_limitedLinear`):
  the 15 baselines and the 4 silica temperature-coordinate cases re-run on
  41,472 cells with `Gauss limitedLinear 1`; all 19 PASS. The baseline
  tables, Figs. 4 and 7 and the abstract numbers now come from them. The 64
  follow-up sensitivities stay 18,432-cell upwind and are reported as
  differences from their own upwind baselines (stated in Secs. 4.2, 5.3, 7.3).
  Changes that matter: the pure-PbI2 binned peak is 36 cm / 366.5 C
  (published 36 cm / 368 C) but its 3 cm smoothed maximum and centroid
  temperature are 350.5 / 348.7 C; steel PbI2 smoothed / centroid
  temperatures are 269.2 / 274.5 C (lower edge of 288 +- 18 C) while the
  binned value is 288.1 C; the silica overlap ranking of the two coordinate
  shifts reverses (Q1 77.80 vs 73.26 %, Q2 77.88 vs 75.54 %; upwind had the
  opposite order), so neither shift is preferred. The BiI3 share, totals and
  ledger are unchanged within the stated tolerances (largest raw residue
  1.41e-12).
- **Scheme and mesh study** (`run/024-liu2025/scheme_mesh_study`):
  Table III and Fig. 3 replaced. Steel shifted-BiI iodine profile, successive
  levels 8,192 / 18,432 / 41,472 / 93,312 cells: limitedLinear 15.35 / 8.95 /
  1.98 %, upwind 10.47 / 10.58 / 2.13 %. Apparent orders 1.3 and 3.7 are
  inconsistent: no formal GCI; the 1.98 % L3-L4 change is quoted, with
  indicative 2.5 % (apparent order) and 3.6 % (second order) estimates of the
  41,472-cell error. Upwind vs limitedLinear at equal mesh 20.8 / 14.1 / 11.3 /
  7.7 %; van Leer vs limitedLinear at most 2.8 % (steel) and 1.8 % (pure).
- **Iteration check** (`run/024-liu2025/iteration_check`): eight vs three
  correctors over the full 60 s, three cases. Iodine profile change 1.1e-13 %
  and 4.9e-13 % (pure PbI2, steel 0.1 %) and 0.021 % (steel 1 %). The
  "full-run iteration sensitivity unqualified" caveat is replaced.
  The fourth case (shifted BiI) was still queued and is not used.
- **Lobresco et al. (2026)** (`run/023-lobresco2026`,
  `thermochemistry/benchmarks/lobresco2026`): new Sec. 6.3, Table V,
  Fig. 5, one sentence each in abstract and conclusions. Onset temperature
  minus measurement +26 / +19 / +16 / -11 / -14 K, RMS 18 K, against 35 K for
  both published models; coarse mesh within 0.3 K; peaks 1-4 cm upstream;
  L1/n0 0.4-0.9. T5_60 uses the main deposit (from the inlet the rule hits a
  back-diffused upstream deposit at 548 K). Not represented: flow stop and
  withdrawal, the non-dilute source region (mole fraction 0.33-0.36), the
  boat. The T-Flows test-case wall table differs from the published curves;
  it is described, not copied.
- **Figures** redrawn at printed size (5.4 in wide or less, 8-9 pt text) and
  included without scaling, so the PoF layout has no text below 8 pt.
- **PoF title** changed to the transport framing (decision below).

## Numerical evidence: current state

Time step 2 to 1 ms: 0.0043-0.3228 %. Correctors 3 to 8 over 60 s: at most
0.021 %. Space: see above; the upwind numerical diffusion is the dominant
spatial error of the original campaign, and the production profiles differ
from the 18,432-cell upwind profiles by 9.1-22.3 %. Integral deposition is
much more stable than local shape. Mesh independence of local profile shapes
is not claimed; the L3-L4 check exists for one case only.

Not done: the upwind 93,312-cell runs of the other follow-up cases in
`~/lesto-submission-followup` were still running and are not used; a
manufactured-solution test; a T-Flows run with identical inputs.

## Title (PoF only)

`pof/lesto-pof.tex` now reads "Competition between radial diffusion and
interfacial exchange in thermal-gradient deposition columns: conservative
modeling of iodine transport". Basis: the wall Damkoehler number is
1.4e3-2.9e3 at alpha = 1; alpha = 0.1 changes profiles by at most 0.41 %,
alpha = 0.01 by at most 4.21 %, halving or doubling D by 3.83-14.20 %; and
the untuned HKS wall law with the near-wall diffusion resistance gives the
best onset RMS on the Lobresco cases. The arXiv title is unchanged. Revert
if the authors prefer the descriptive title.

## Experimental and thermodynamic limits

Everything currently available experimentally came from the supplied Liu
paper. The required independent inputs are not present in the workspace:
original T(x)/scan coordinates, controller reference conditions, source release
history and a qualified BiI record. Do not describe these as measurements
already obtained or silently supply substitutes.

Q1 species masses come from observed deposits, and Q2 retains that iodine
loading while assuming ideal activities and unmeasured metal saturation.
The +29 kJ/mol BiI gas shift is a sensitivity surrogate. The constructed BiI3
liquid record and Bi thermodynamic provenance need independent qualification.
KI, water/oxygen chemistry, surface-material chemistry, gas reaction kinetics
and boat obstruction are outside the represented physics.

The 60 s source-driven cases retain the three-hour mean release rate.
The four 60–80 s restarts hold a nonuniform temperature profile fixed after
stopping flow and source. They do not reproduce three-hour evaporation or
the experimental cooling protocol. These omissions are explicit in the draft.

The Lobresco comparison uses digitised wall curves, not the thermocouple
tables; those tables, the boat geometry and position, whether the figure
origin is the column inlet, and the C_e behind their Table 2
vapour-pressure column (the paper calls C_e = 5 illustrative, not
calibrated, and does not say which value Table 2 used) remain open
questions for the authors of that study.

## Before circulation as a journal submission

- Confirm author list/order, affiliation, funding, acknowledgments and any
  contribution/conflict statements with the coauthors. The existing two-author
  list is retained for this working draft; no messages have been sent.
- Review the scientific scope and whether further numerical verification is
  sufficient for the intended methods contribution.
- Confirm thermodynamic record provenance and the redistribution terms for any
  future CFD input package. External coefficient files are not included here.
- A Physics of Fluids version exists (`pof/`, checklist in
  `pof/SUBMISSION-CHECKLIST.md`); nothing has been submitted or published.

The prior `notes-for-authors.md` and `literature-review.md` are historical
planning notes, not evidence that their proposed experiments or code comparisons
have been performed. This file and the manuscript describe the current state.
