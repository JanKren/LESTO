# Preliminary M10d metal-vapour channel

The Q2 source releases monatomic Pb, Bi and I into helium. Its iodine amount
uses the deposited PbI2/BiI3 amounts of Liu's LBE-I_SS_II case; the added
metal vapour is 1% of the ideal-LBE saturation pressure at 1173 K. The
release rate and plane-channel carrier are prescribed numerical inputs.
This case qualifies the new wall channels; measured T(x), carrier flow and
release information remain necessary for experimental validation in M10e.

Bi_g and Bi2_g share one Bi_s deposit. The dimer deposits two Bi atoms.
The irreversible channel `3 BiI_g -> 2 Bi_s + BiI3_g` shares that store,
and introduces its gas product after every gas has completed transport.
Pb_g has its own Pb_s deposit. PbI2 and BiI3 retain their existing pairs.
One field stores the metal across its solid/liquid transition; the stable
phase is selected by the minimum of the two equilibrium pressure curves.

Generate the tables outside the repository, then run:

```sh
python3 ../../thermochemistry/systems/make_wall_tables.py \
  --nasa /external/nasa_thermo.inp --tf /external/bi_i_m10-thermofun.json \
  --out /external/m10d-baseline
export LESTO_M10D_DATA=/external/m10d-baseline
../../applications/rhoFixedFlowFoam/Allwmake
./Allrun
```

`--bii-shift 29000` generates a consistent HSC-like sensitivity surrogate:
the gas formation, pure BiI pressure and disproportionation tables all use
the same shifted Gibbs energy. It is not a licensed HSC record. No external
thermodynamic records are bundled. `benchmark.json` describes the inputs.

The uniform 200 x 2 mesh and zero molecular diffusion approximate the
fractional plug-flow limit, with HKS accommodation and kinetic coefficients
of one. The measurement at 1.2 s is during the two-second source window.
A subsequent clean-helium flush evaporates hot deposits and changes their
onsets, so it is a separate physical comparison. Onset is the hottest bin
above 0.1% of each condensate's peak, with metal phases classified at the
melting temperature. Bin temperatures come from the written profiles.

The automated qualification also uses a 0.5 ms step, two BiI variants and
an independent constrained-equilibrium/fractional-condensation algorithm.
It shares the record evaluator with table generation and restricts the
reference to the same nine gases and eight condensed phases. Run it with:

```sh
export LESTO_M10D_NASA=/external/nasa_thermo.inp
export LESTO_M10D_THERMOFUN=/external/bi_i_m10-thermofun.json
export LESTO_M10D_REFERENCE=/external/m10-design
../../applications/rhoFixedFlowFoam/tests/Alltest M10d
```

The reference directory contains `design/py/pbbii.py` and its dependencies.
`verification.json` and `out_excerpt` record the completed qualification.
Shared-store profiles are `condensate_Bi_s` and `condensate_Pb_s`; element
profiles include each store once. Each condensate has a separate balance
in `phaseChangeChannelCondensates`, in addition to species/element balances.

All four M10d gates pass. The combined M0/M10a–d invocation records 43 PASS,
0 FAIL and 13 NOTRUN: M0.4 requires a v2606 suite invocation, and earlier
external M10a–c studies were not selected for this regression run. A separate
final v2606 Debug compilation and startup pass with zero warnings. The
default and GEMS v2412 builds also have zero warnings.

Both 1 ms and 0.5 ms onset comparisons differ from plug flow by at most
3.89 K. The half-step changes no baseline onset and shifts the sensitivity
case's Pb onset by one 4.115 K bin. Element closure is below 4.42e-16 across
all four runs. Deposit-profile L1 changes range from 0.29% to 6.45%; these
are reported, not gated. Quantitative experimental profiles need further
time-step refinement as well as the physical inputs for M10e.
