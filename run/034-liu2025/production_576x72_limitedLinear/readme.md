# Production runs: 41,472 cells, limitedLinear (plan sections 38-39)

The 15 completed baseline cases of `../full_analysis` and the four silica
temperature-coordinate cases of the follow-up study, re-run on the 576 x 72
(41,472-cell) mesh with `Gauss limitedLinear 1` for the advection of all nine
gases.  Everything else is byte-identical to the 18,432-cell upwind inputs
(sources, carriers of the follow-up `__grid576x72` cases, coefficients, 2 ms,
60 s).  Solver: the M10e binary, sha256 `16fd804e...` (as `../followup_results`).
All 19 runs PASS their numerical gates.

- `analysis.json`: the manuscript metrics (`q1_peaks` with 1 cm binned
  maximum, 3 cm smoothed maximum and centroid wall temperature;
  `case_metrics` of the 15 baselines, same keys as `../full_analysis`).
- `ledger.json`: per-case closure, raw material residue, solver defect,
  outer-iteration counts and run time, next to the 18,432-cell values.
- `profiles/<case>/`: final element-I, PbI2 and BiI3 wall profiles and the
  iodine comparison with the digitised scan.

The 64 follow-up sensitivities stay on 18,432 cells with upwind and are
compared with their own 18,432-cell upwind baselines in `../full_analysis`.
