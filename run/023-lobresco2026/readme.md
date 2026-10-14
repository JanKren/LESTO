# Case 023 — Head-to-head on Lobresco et al. (2026): five PbI2 experiments

This case is based on case 017.

Case 017 computes one experiment of Lobresco et al., Nucl. Eng. Des. 459
(2026) 115201, on one mesh.  Here are all five PbI2 thermosublimatography
experiments of their Table 1, T106_6, T106_60, T106_300, T5_60 and T540_60
(flow in standard mL/min, duration in s), on two meshes, compared with the
measured deposits and with the paper's own two models.  No parameter is
adjusted to these data.  The comparison is Sec. 6.3, Table V and Fig. 5 of
the manuscript (`manuscript/prepare_results.py` reads this folder).

## Differences from case 017

The study, not the model: five experiments (flow 5, 106 or 540 mL/min;
duration 6-300 s; 2 ms steps, 0.5 ms for T540_60; the boat at the
digitised arrow, 0.13 m at 106 and 0.15 m at 5 and 540 mL/min), each on
the production wedge (318 x 72 = 22,896 cells, the spacing of the
41,472-cell Liu mesh) and on the coarse wedge of case 017 (212 x 48).
Wall law, carrier technique, numerics and observable are those of case
017; the set-up table (what is given and what is assumed) is in notes.md.

## Files

    make_inputs.py, inputs/   wall temperatures from the digitised curves
    make_cases.py             carriers (built and audited) and cases
    analyse.py                the paper's 1 %-of-peak rule, comparison.json
    launch_full.sh            how the recorded runs were made (in place,
                              one core per case)
    cases/<level>/<case>/     the dictionaries, case.json, comparison.json,
                              ledger_summary.json and final_profiles/ of
                              the recorded runs (meshes, fields, logs not kept)

The digitised measurements, wall curves and both published model curves
are in `thermochemistry/benchmarks/lobresco2026` (CC BY 4.0).

## Run

    ./Allrun                       # coarse T106_6 = case 017 (about 8 min)
    ./Allrun prod T106_6 T540_60   # selected cases on the production mesh
    ./Allrun coarse                # all five (T106_300: 150,000 steps;
                                   # T540_60: 120,000 steps; hours per case)
    ./Allclean

Every case runs into `work/<level>/cases/<case>` and is checked against its
tracked record: T_dep, wall inventory, closure and the final profiles.

## What to check

Deposition-onset temperature (the paper's rule: 1 % of the peak of the
deposit, from the inlet, on the measured wall curve) minus the measurement
of Table 2, production mesh:

| case | T106_6 | T106_60 | T106_300 | T5_60 | T540_60 | RMS |
|---|---|---|---|---|---|---|
| this work | +26 K | +19 K | +16 K | -11 K (main deposit) | -14 K | 18 K |
| paper, both models | | | | | | 35 K |

- The coarse mesh agrees within 0.3 K (T106_6: 839.086 K against 839.112
  K, case 017).
- At 5 mL/min vapour diffuses upstream and deposits near the cold inlet;
  applied from the inlet the rule finds that small deposit (548 K), so the
  table uses the main deposit downstream of the boat.  `comparison.json`
  keeps the value from the inlet.
- Computed peaks lie 1-4 cm upstream of the measured ones.
- Largest |closure| over all steps and runs: 3.3e-16 of the inventory.

Solver of the recorded runs: the M10e binary, sha256
`16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58`.

Details: notes.md (the set-up table, what is not represented, the files).
