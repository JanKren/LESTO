# Upwind 864x108 runs of the four check bases (plan section 39, item 9)

The `lesto-submission-run` package (output `~/lesto-submission-followup`,
not in git) ran the four bases of the iteration and mesh checks on level L4
(864x108 = 93,312 cells) with the species advection of the original campaign,
`Gauss upwind`, at 2 ms to 60 s, solver sha256 `16fd804e...` (the M10e binary).
All four PASS: element gates (closure at most 1.9e-15), no re-speciation
failure, carrier continuity audit PASS.

- `profiles/<base>/{element_I,PbI2_g,BiI3_g}_native.dat`: final (60 s) profiles.
- `summary.json`: per base, the 41,472 -> 93,312-cell change of the iodine,
  PbI2 and BiI3 wall profiles (relative L1 of the 1 cm wall amounts,
  normalised by the finer level, as in `../results.json`), the total and
  centroid change, against the upwind 576x72 follow-up case
  (`../../followup_results/profiles/<base>__grid576x72`).
- `collect.py`: `python3 collect.py ~/lesto-submission-followup` rebuilds both.

| base | iodine L1, L3 -> L4 | total change | centroid shift |
|---|---|---|---|
| `PbI2_SS_Q1_baseline` | 3.68 % | 0.008 % | -0.39 mm |
| `LBE-I_SS_II_Q2_baseline_sat0.001` | 7.30 % | 0.015 % | -0.42 mm |
| `LBE-I_SS_II_Q2_baseline_sat0.01` | 6.68 % | 0.015 % | -0.36 mm |
| `LBE-I_SS_II_Q2_shift29k_sat0.01` | 2.13 % | 0.018 % | -0.21 mm |

The shifted-BiI run repeats the study's own upwind L4 run
(`../profiles/steel_L4_up`) and reproduces it bitwise.  The other three bases
change more between L3 and L4 than the shifted-BiI case, so that case, chosen
for the study because it had the largest upwind change between L2 and L3, is
not the most mesh-sensitive case at the finest step.  These are upwind
results only; the production scheme (limitedLinear) is refined to L4 for the
shifted-BiI case only.
