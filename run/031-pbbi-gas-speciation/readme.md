# Preliminary M10b gas re-speciation channel

This case adds homogeneous Pb-Bi-I gas equilibrium to the prescribed
PbI2/BiI3 source of run/030. It uses nine gases (including Bi2), two plain
condensation pairs, a synthetic 1173–350 K plane-channel ramp and a common
provisional diffusivity of 1e-4 m2/s. The source lasts 2 s; the calculation
ends at 3 s. The `cm` profiles use 1 cm bins for the deposit comparison; `native`
profiles use 1 mm bins for onset and peak reporting.

It qualifies the transport/speciation implementation. It does not predict
the experimental source split: the release history, the actual carrier
and wall temperatures, and iodine activity in LBE are still missing.
Metal condensation and reactive BiI wall channels belong to M10d.

Load OpenFOAM, build `applications/rhoFixedFlowFoam/Allwmake`, and supply
external Bi data. Generate the formation table alongside the vapour table:

```sh
export LESTO_M10_DATA=/external/m10-tables
python3 ../../thermochemistry/systems/make_kf_tables.py \
  --nasa /external/thermo.inp --tf /external/bi_i_m10-thermofun.json \
  --out "$LESTO_M10_DATA/log10Kf_PbBiI.csv"
./Allrun
```

`LESTO_M10_DATA` must also contain `pv_BiI3_phases.csv`, made with
`make_pv_tables.py`. These derived Bi data remain outside the repository.
`Allrun` prepares and runs `work/case`, leaving the tracked inputs intact.
Results include species and element accounts, condensate profiles and a
re-speciation report. The default solver build requires no GEMS bridge.

The preliminary study in `verification.json` compares this configuration
with an independent-pair baseline at the same 1173 K inlet. At 1 ms the
1 cm deposit profiles differ by 0.0000844% for PbI2 and 0.1332% for BiI3.
The 2 ms versus 1 ms differences are 0.0331% and 0.0575%, respectively.
These are numerical comparisons under the specified provisional inputs.
Run `tests/Alltest M10b` with the external-data environment variables to
repeat the full 4/2/1 ms study and the conservation checks.

The 1 mm profiles report preliminary 1% onsets of 559.78 K (PbI2) and
444.56 K (BiI3), with bin peaks at 543.405 K and 436.415 K. The maximum
element closure in this run is 6.3e-16 relative. These values depend on
the prescribed inputs; they are not an experimental validation.
