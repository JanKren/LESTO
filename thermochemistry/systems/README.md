# M10a generators and external Bi inputs

`make_gems3k_input.py --system PbBiIHe --tf /external/bi_i_m10-thermofun.json`
adds Bi-I ThermoFun records to the existing NASA Pb-I-He system. Its default
`PbI2He` outputs retain the original bytes. The code checks record formulas,
Cp integrals and transition consistency. `--lbe-solution` prepares an ideal
Pb–Bi liquid variant for later work; M10a's transport case uses pure paired
condensates and does not activate a liquid-solution model.

Generate tables into an external directory:

```sh
python3 make_pv_tables.py --nasa /external/thermo.inp \
  --tf /external/bi_i_m10-thermofun.json --out /external/m10-tables \
  --check --regress data/pv_PbI2_phases.csv
```

The generator writes BiI3, Bi, Bi2, BiI and Pb tables, with each condensed
phase evaluated separately and crossings included as nodes. The solver
interpolates each column in inverse temperature and chooses the minimum.
Bi2's table describes `2 Bi(cond) = Bi2(g)`; its shared-condensate coupling
is **not supported by M10a**, which requires identical formulas in each
one-to-one pair. Generating a table does not qualify the corresponding
transport channel. Bi-I source records and derived tables remain external
until their provenance and redistribution rights are confirmed. This
repository contains the generator code, not those records or tables.

Acceptance checks accept `LESTO_M10_NASA`, `LESTO_M10_THERMOFUN` and
`LESTO_M10_DATA` as explicit external input paths; unavailable data are
reported as NOTRUN. The code tests use labelled synthetic tables.

## M10b gas formation constants

```sh
python3 make_kf_tables.py --nasa /external/thermo.inp \
  --tf /external/bi_i_m10-thermofun.json \
  --out /external/m10-tables/log10Kf_PbBiI.csv
```

Each column is `log10 Kf` for a named solver gas, formed from monatomic
Pb, Bi and I gases at the 1 bar standard state. The generator evaluates
exactly the same NASA/ThermoFun functions as the GEMS3K system generator.
The kernel interpolates in inverse temperature; monatomic constants are
zero. The table covers 270–1250 K, with 5 K nodes by default. Without
`--tf`, the generator writes only the open Pb-I gases. `--synthetic`
produces explicitly labelled invented constants for code tests.

M10b checks additionally accept `LESTO_M10_DESIGN` for external stored
reference points and gas-map inputs, and `LESTO_M10_ORACLE_LIB` for the
extended design bridge used only as a test oracle. Production GEMS engines
are a separate M10c milestone; ordinary solver builds remain bridge-free.
