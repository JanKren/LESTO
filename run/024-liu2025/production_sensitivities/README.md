# Production-scheme sensitivities: 576 x 72, limitedLinear (plan section 39)

The diffusivity and accommodation sensitivities of the follow-up study,
repeated with the production scheme for two of the 15 baselines: the pure
PbI2 steel column (`PbI2_SS_Q1_baseline`) and the steel LBE column with Q2
at 1 % metal saturation (`LBE-I_SS_II_Q2_baseline_sat0.01`).  Each variant is
compared with its own production baseline in
`../production_576x72_limitedLinear`, and next to the upwind 18,432-cell
value of the same variant in `../followup_results`.

## Set-up

`make_cases.py` writes the eight cases from the two production base cases.
They change only what `thermochemistry/benchmarks/liu2025/followup_study.py`
changed for the upwind follow-up cases (confirmed by diffing those cases
against their bases):

- `__D5e-05`, `__D0.0002`: the diffusivity of all nine gases halved or
  doubled (all nine `D 0.0001;` entries);
- `__alpha0.1`, `__alpha0.01`: all three `accommodation 1;` entries (the
  PbI2_g and BiI3_g pairs and `HKSCoeffs`).

The table paths are repointed to the baseline tables; the tables are
byte-identical to those of the production runs.  576 x 72 = 41,472 cells,
`Gauss limitedLinear 1`, 2 ms, 60 s, serial.  Solver: the M10e binary,
sha256 `16fd804e...` (as `../production_576x72_limitedLinear`).  Run
14-15 October 2026; 17.4-22.2 h (pure PbI2) and 24.7-28.0 h (steel) each.
All eight reached 60 s and ended normally; the end-of-run ledgers are in
`ledger.json` (pair closures at most 3.7e-16, element residues of the
re-speciation at most 5.1e-16; the counts of steps without outer
convergence equal those of the bases).

## Results (`comparison.csv`, `analysis.json`)

Iodine wall profile at 60 s against the variant's own base: relative L1,
change of the total deposit, centroid shift and change of the BiI3 share
3 Bi/(2 Pb + 3 Bi) of wall iodide iodine (base 94.58 % limitedLinear,
94.56 % upwind).  "LL" is this study, "upwind" the 18,432-cell follow-up.

| case | L1 LL | L1 upwind | total LL | total upwind | centroid LL [cm] | centroid upwind [cm] | BiI3 share LL [pp] | upwind [pp] |
|---|---|---|---|---|---|---|---|---|
| PbI2 D/2 | 16.65 % | 14.20 % | -0.016 % | -0.016 % | +0.230 | +0.230 | 0 | 0 |
| PbI2 2D | 8.98 % | 7.25 % | +0.007 % | +0.007 % | -0.085 | -0.088 | 0 | 0 |
| PbI2 alpha 0.1 | 0.49 % | 0.41 % | -0.000 % | -0.000 % | +0.005 | +0.005 | 0 | 0 |
| PbI2 alpha 0.01 | 5.14 % | 4.21 % | -0.005 % | -0.005 % | +0.056 | +0.057 | 0 | 0 |
| steel D/2 | 1.07 % | 4.61 % | -0.020 % | -0.024 % | +0.003 | +0.032 | -0.028 | -0.067 |
| steel 2D | 12.53 % | 6.68 % | +0.017 % | +0.023 % | +0.069 | +0.036 | +0.002 | +0.012 |
| steel alpha 0.1 | 0.15 % | 0.36 % | -0.001 % | -0.001 % | +0.001 | +0.002 | -0.000 | -0.000 |
| steel alpha 0.01 | 1.90 % | 3.67 % | -0.008 % | -0.009 % | +0.012 | +0.024 | -0.000 | -0.000 |

Steel, per species (L1 LL / upwind): D/2 PbI2 9.26 / 8.22 %, BiI3 0.60 /
4.41 %; 2D PbI2 4.33 / 2.65 %, BiI3 13.00 / 6.91 %.

- The ordering holds with the production scheme: accommodation changes
  move the iodine profile by 0.15-0.49 % (alpha 0.1) and 1.90-5.14 %
  (alpha 0.01), the diffusivity by 1.07-16.65 %.  Total deposits change by
  at most 0.02 % and the BiI3 share by at most 0.03 percentage points.
- The pure PbI2 response is close to the upwind one and nearly symmetric
  in D.  In the steel case it is not: halving D hardly moves the
  BiI3-dominated profile (BiI3 L1 0.60 %), doubling it gives 12.53 %.  The
  logs confirm the halved value was read.  The D response is therefore
  case-dependent and not linear; an intermediate D would be needed for a
  per-case statement about its size.

## Files

| file | content |
|---|---|
| `make_cases.py` | writes the eight cases from the two production bases |
| `collect.py` | copies the 60 s profiles and the end-of-run ledger from the runs |
| `profiles/<case>/` | element-I, PbI2 and BiI3 wall profiles at 60 s |
| `ledger.json` | end time, run time, final pair balances, outer-convergence counts, element residue |
| `analyse.py` | writes `comparison.csv` and `analysis.json` from the tracked profiles |

`manuscript/prepare_results.py` builds the manuscript table from
`comparison.csv`, and `manuscript/check_results.py` recomputes it from the
profiles.
