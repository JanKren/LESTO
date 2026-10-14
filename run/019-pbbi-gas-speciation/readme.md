# Case 019 — Homogeneous Pb-Bi-I gas re-speciation

This case is based on case 018.

In case 018 PbI2 and BiI3 travel as fixed molecules.  At the temperatures
of an LBE furnace the gas re-speciates: PbI2 and BiI3 partly dissociate
into PbI, Pb, BiI, Bi, Bi2, I and I2, and recombine as the gas cools.  Here
nine gases are kept in homogeneous equilibrium in every cell and step by an
exact kernel with tabulated formation constants; PbI2 and BiI3 still
condense by their HKS pairs.  The source is that of case 018.

## Differences from case 018

In `constant/thermochemistryProperties`

    respeciation {engine kernel;
      species (Pb_g PbI_g PbI2_g Bi_g Bi2_g BiI_g BiI3_g I_g I2_g);
      formationConstants "$LESTO_M10_DATA/log10Kf_PbBiI.csv";
      minTemperature 350; updateInterval 1;}

What follows from it:
- seven more gases (D 1e-4 m2/s, each with its `formula`), their fields,
  schemes and solvers (`tolerance 1e-18`);
- the inlet at 1173 K (`Tin` in `system/carrierDict`), the LBE furnace, so
  that there is something to re-speciate; the comparison baseline is case
  018's set-up at 1173 K;
- 1 cm `cm` bins besides the native 1 mm bins, for the deposit comparison.

## Modified files

    constant/thermochemistryProperties, constant/speciesTransportProperties
    0/Y_Pb_g, 0/Y_PbI_g, 0/Y_Bi_g, 0/Y_Bi2_g, 0/Y_BiI_g, 0/Y_I_g, 0/Y_I2_g
    system/carrierDict, system/controlDict, system/fvSchemes,
    system/fvSolution

## Run

    export LESTO_M10_DATA=/external/m10-tables
    python3 ../../thermochemistry/systems/make_kf_tables.py \
      --nasa /external/thermo.inp --tf /external/bi_i_m10-thermofun.json \
      --out "$LESTO_M10_DATA/log10Kf_PbBiI.csv"
    ./Allrun            # the 3 s at 1 ms (about 3-4 min), then expect.json
    ./Allrun dt         # the same at 2 and 4 ms, deposit L1 against 1 ms
    ./Allclean

`LESTO_M10_DATA` must also hold `pv_BiI3_phases.csv` (case 018).  The
default solver build is enough (no GEMS3K bridge).  `tests/Alltest M10b`
repeats the whole 4/2/1 ms study and the conservation checks.

## What to check

At 3 s, the 1 mm profiles (`verification.json`, `run031_native_profiles`):

| pair | 1 % onset | peak bin |
|---|---|---|
| PbI2 | 559.778 K | 543.405 K |
| BiI3 | 444.563 K | 436.415 K |

- Element closure in every step at round-off (6.3e-16 recorded; 4.3-5.0e-16
  measured on 2026-10-14 with the M10e solver).
- Against the independent pairs of case 018 at the same 1173 K inlet, the
  1 cm deposits differ by 0.0000844 % (PbI2) and 0.1332 % (BiI3) at 1 ms;
  2 ms against 1 ms gives 0.0331 % and 0.0575 % (`./Allrun dt`).

These are numerical results under provisional inputs (synthetic ramp,
common diffusivity, prescribed release), not a prediction of the source
split.  Metal condensation and reactive wall channels come in case 021.
