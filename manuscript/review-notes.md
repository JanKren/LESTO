# Author review: 9 October 2026

The completed working draft is framed as a conservative CFD method plus
numerical-resolution and input-sensitivity assessment. It reports 79 completed
computational cases, not 79 independent experiments, and makes no claim of
independent predictive validation. The missing BiI3_SiO2_800 temperature profile
is stated explicitly; that case was not run.

## What changed from the planning draft

- The title, abstract, model description, verification/results tables and
  conclusions now describe completed work.
- All reported cases use `respeciation { engine kernel; }`. GEMS3K is a software
  capability, rather than the per-cell engine used for these results.
- The executed scheme is implicit Euler and first-order upwind on a 5-degree
  wedge. The main mesh has 18,432 cells; separate 8,192/41,472-cell and 1 ms
  studies are reported. No SuperBee, 228,960-cell mesh, unrun manufactured
  solution or completed T-Flows comparison is claimed.
- Actual case manifests use 4.8 mm steel and **5.0 mm silica** diameters. Earlier
  presentation wording about a common 4.8 mm diameter should be corrected when
  that deck is next revised. Boat obstruction is still omitted.
- The original silica mismatch remains visible. Its two temperature-coordinate
  shifts are input-derived hypotheses, not a confirmed correction or fitted
  validation result.
- Funding statements and grant identifiers from the earlier draft were removed
  until the authors confirm that they apply to this work.

## Numerical evidence still needed

The isolated 2/1 ms comparisons change iodine profiles by 0.0043–0.3228%.
The 18,432/41,472-cell comparisons change them by 2.06–10.58%. Integral
deposition is much more stable than local shape. Three-grid estimates are
conditional or refused for oscillation/nondecreasing differences, so the
manuscript does not claim mesh independence or confirmed GCI uncertainty.

Further spatial evidence should focus on the shifted-BiI steel case with
10.58% fine/medium change and the pure-PbI2/low-loading branches whose binned
maxima switch. A further refinement or a separately verified bounded
higher-order scheme could address this. Full-run corrector sensitivity remains
needed if iteration independence is to be claimed; the saved 40 ms restarts
only diagnose the tested late-time checkpoints. No new long CFD study was
launched while preparing this draft.

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

## Before circulation as a journal submission

- Confirm author list/order, affiliation, funding, acknowledgments and any
  contribution/conflict statements with the coauthors. The existing two-author
  list is retained for this working draft; no messages have been sent.
- Review the scientific scope and whether further numerical verification is
  sufficient for the intended methods contribution.
- Confirm thermodynamic record provenance and the redistribution terms for any
  future CFD input package. External coefficient files are not included here.
- Select a journal/template after that review. No journal is selected and
  nothing has been submitted or published.

The prior `notes-for-authors.md` and `literature-review.md` are historical
planning notes, not evidence that their proposed experiments or code comparisons
have been performed. This file and the manuscript describe the current state.
