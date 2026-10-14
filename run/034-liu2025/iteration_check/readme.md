# Full-run corrector check (plan section 39)

Eight instead of three species correctors over the whole 60 s, with every
other input identical (18,432 cells, `Gauss upwind`, 2 ms, M10e solver sha256
`16fd804e...`).  Three cases completed (PASS): `PbI2_SS_Q1_baseline`,
`LBE-I_SS_II_Q2_baseline_sat0.001`, `LBE-I_SS_II_Q2_baseline_sat0.01`.
The fourth (`..._shift29k_sat0.01`) was not finished and is not used.

- `summary.json`: per case, the iodine/PbI2/BiI3 profile change (relative L1
  of the 1 cm wall profile against `../full_analysis/profiles/<case>`), total
  and centroid change, the element-I gate and the outer-iteration counts.
- `profiles/<case>__correctors8_full60s/`: the final (60 s) profiles.

The runs were made with the `lesto-submission-run` package on the workstation
(output `~/lesto-submission-followup`, not in git).  The two first cases agree
with their three-corrector references to round-off (iodine profile L1 about
1e-13 %); the 1 % steel case differs by 0.021 %.  With eight correctors the
largest outer initial residual is still 2.4e-15 to 4.4e-15 and nearly every
step stops at the corrector limit, as with three: the requested 1e-16 is
below the attainable residual floor.
