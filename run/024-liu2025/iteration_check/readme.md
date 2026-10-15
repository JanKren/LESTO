# Full-run corrector check (plan section 39)

Eight instead of three species correctors over the whole 60 s, with every
other input identical (18,432 cells, `Gauss upwind`, 2 ms, M10e solver sha256
`16fd804e...`).  All four cases completed (PASS): `PbI2_SS_Q1_baseline`,
`LBE-I_SS_II_Q2_baseline_sat0.001`, `LBE-I_SS_II_Q2_baseline_sat0.01` and
`LBE-I_SS_II_Q2_shift29k_sat0.01` (the shifted-BiI branch, added 2026-10-15).

- `summary.json`: per case, the iodine/PbI2/BiI3 profile change (relative L1
  of the 1 cm wall profile against `../full_analysis/profiles/<case>`), total
  and centroid change, the element-I gate and the outer-iteration counts.
- `profiles/<case>__correctors8_full60s/`: the final (60 s) profiles.
- `collect.py`: adds a finished case from the package output to `summary.json`
  (`python3 collect.py ~/lesto-submission-followup <case>`); it reproduces the
  three earlier records from the same output.

The runs were made with the `lesto-submission-run` package on the workstation
(output `~/lesto-submission-followup`, not in git).  The two first cases agree
with their three-corrector references to round-off (iodine profile L1 about
1e-13 %); the 1 % steel case differs by 0.021 %, and the shifted-BiI case by 1.1e-7 %
(BiI3 profile 1.5e-6 %).  With eight correctors the
largest outer initial residual is still 2.4e-15 to 4.4e-15 and nearly every
step stops at the corrector limit, as with three: the requested 1e-16 is
below the attainable residual floor.
