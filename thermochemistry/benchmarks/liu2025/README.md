# Preliminary M10 benchmarks: Liu et al. (2025)

Reproducible experimental targets, a published-equation baseline, and four
axisymmetric mesh/input staging directories for M10e. These are preliminary
benchmarks: no CFD validation result is claimed. The current solver still needs
the M10 Bi-species extensions, and the measured wall-temperature profiles and
physical carriers are not supplied by this package.

Source: Lu Liu et al., *Chemical speciation of iodine evaporated from liquid
lead–bismuth eutectic investigated by thermosublimatography*, J. Radioanal.
Nucl. Chem. **334**, 5201–5216 (2025),
[doi:10.1007/s10967-025-10246-4](https://doi.org/10.1007/s10967-025-10246-4).
Tabulated values are transcribed and adapted under the paper's CC BY 4.0 licence;
retain this attribution. No HERACLES/HSC input files, unpublished Bi data, or PDF
are copied into this package. The four Table 1 sublimation constants are the
published Eq. (3) baseline, not a production thermodynamic dataset.

## Primary benchmarks

| Experiment | Column inner diameter | He flow, reported | Prescribed source for the preliminary mean-release scenario | Measured wall-deposit peak |
|---|---|---|---|---|
| PbI2_SS | 4.8 mm, SS316L | 100 mL/min | 14.6 mg PbI2 | 368 ± 18 °C |
| BiI3_SiO2_800 | 5.0 mm, silica | 100 mL/min | 19.2 mg BiI3 | 206 ± 13 °C |
| LBE-I_SS_II | 4.8 mm, SS316L | 45 mL/min | 0.111 mg PbI2 + 0.181 mg BiI3 | PbI2 288 ± 18 °C; BiI3 154 ± 14 °C |
| LBE-I_SiO2_I | 5.0 mm, silica | 45 mL/min | 0.243 mg PbI2 + 0.085 mg BiI3 | PbI2 306 ± 24 °C; BiI3 127 ± 39 °C |

All experiments lasted 10800 s; both column types were 1.25 m long. The boat is
approximately 3 mm wide and 15 mm long (p. 5204); its precise axial position,
depth, exposed area, and mapping to an axisymmetric source are still required.
The blank `wall_temperature.csv` files must be filled with measured `x_m,T_K`
data in the mesh coordinate system. No linear temperature ramp is substituted.

The seven other Table 2 experiments and all 17 peaks in Tables 3–5 are also
recorded, including KI as a diagnostic outside the Pb–Bi–I model. The typo
`LEB-I_SS_II` in Table 5 is normalised to Table 2's `LBE-I_SS_II`.

The initial preparation passed six standalone checks and all four meshes
passed OpenFOAM v2412 `checkMesh` (18432 cells, maximum non-orthogonality 0°,
maximum skewness 0.331). `prepared/verification.json` records that evidence
against the reference/generator hashes. This is input and mesh verification;
all six primary CFD peak comparisons remain `NOTRUN`.

## Source and observable conventions

- The pure-compound source masses are the initial samples (Table 2). For this
  initial preparation, release them at their mean rate over 3 h (Q1). The later
  physical benchmark uses evaporating inventories (Q3), qualified against the
  unknown evaporation history. The elapsed experimental time is not proof of a
  constant release rate.
- For LBE, the source masses are the **inferred deposited species masses** in
  Table 5, prescribed as Q1 sources to study deposition temperatures. They are
  neither the full LBE sample mass nor an independently measured source split.
  Such a run cannot validate a prediction of the PbI2/BiI3 split. That requires
  Q2, metal vapour, re-speciation, and M10d/e.
- Table 5 assigns 6.6% and 7.0% of the total deposited iodine to KI in the two
  primary LBE runs. Keep these exclusions in the manifest; do not renormalise
  the remaining PbI2/BiI3 yields to 100%. The tables round masses, so derived
  iodine inventories agree only to the rounding of those masses.
- All amounts and rates carry explicit **full-column** and **wedge** values.
  The 5° wedge receives 1/72 of the column's total flow and source. Using the
  full sample amount in a wedge would overstate its loading by 72 times.
- The measured `Tdep` is the temperature at the **maximum deposit**, p. 5205.
  Use `T_peak` from the wall deposit, not the solver's 1% onset. The experiment
  used 1 cm scan steps; uncertainties correspond to temperature variation
  within ±1 cm of the peak. Gamma activity measures iodine and must be
  compared with iodine-weighted profiles, not total compound mass.

## Run the preparation

From the repository root, no OpenFOAM or third-party Python packages needed:

```bash
python3 thermochemistry/benchmarks/liu2025/prepare.py --check
python3 thermochemistry/benchmarks/liu2025/prepare.py \
    --out thermochemistry/benchmarks/liu2025/prepared
python3 -m unittest discover -s thermochemistry/benchmarks/liu2025 -p 'test_*.py'
```

`prepared/` contains the Eq. (3) reference calculation, a comparison table with
six `NOTRUN` entries, provenance hashes, and four `benchmark.json` manifests,
empty temperature templates, and mesh dictionaries. The 384 × 48 × 1 mesh has
18432 cells and a 5° wedge. `controlDict`, `fvSchemes`, and `fvSolution` are for
meshing/checking only; these directories
are **not runnable rhoFixedFlowFoam cases**. Mesh/time-step refinement and the
physical helium carrier remain required. For a generated staging directory:

```bash
blockMesh -case thermochemistry/benchmarks/liu2025/prepared/PbI2_SS
checkMesh -constant -case thermochemistry/benchmarks/liu2025/prepared/PbI2_SS
```

Generation rewrites the named output files, including blank temperature
templates. Preserve filled profiles elsewhere and regenerate into a fresh
directory. A custom `--out` can keep working files under `tests/work/`.

### Flow-reference uncertainty

The Liu paper names standard flow but does not specify a numerical reference
temperature/pressure in its methods. The package therefore leaves physical
helium mass flow unset. Provide both values explicitly after confirmation, or
to study a clearly labelled assumption:

```bash
python3 thermochemistry/benchmarks/liu2025/prepare.py --out /tmp/liu-stp-assumption \
    --flow-reference-temperature 273.15 --flow-reference-pressure 101325
```

`eq3_baseline.csv` reports two conventions. An inferred molar volume of
24.465 L/mol (approximately 25 °C, 1 atm) reproduces all 17 rounded published
Eq. (3) temperatures to 0.55 K; this is an inference, not calibration evidence.
The second column uses 0 °C, 1 atm. Eq. (3) uses constant sublimation H/S and
dilute partial pressure relative to its 1 bar thermodynamic standard state.
These temperatures represent an equilibrium threshold, not a computed CFD
peak, and the constant-H/S baseline omits the stable-liquid phase treatment.
The reproduction tolerance is a transcription/convention check, not a
criterion for agreement with the experiment.

## Compare future CFD results

Copy `prepared/predictions-template.csv` to a separate results file and add
rows with columns `experiment,species,peak_temperature_K` (species `PbI2` or
`BiI3`, without `_g`/`_s`). Supply the six primary peaks; partial results remain
`NOTRUN`. Evaluate wall-deposit peaks consistently with the 1 cm experiment,
and record the binning/interpolation and sampling time with the CFD run.

```bash
python3 thermochemistry/benchmarks/liu2025/prepare.py --out /tmp/liu-comparison \
    --predictions my-cfd-peaks.csv
```

`comparison.csv` reports differences against both the experimental uncertainty
and the proposed M10e science targets (±20 K PbI2, ±40 K BiI3). These science
targets are reported, not code gates. Missing entries do not count as success;
duplicate/unknown entries and non-finite/non-positive temperatures are errors.
This peak comparison does not check element conservation, profile agreement,
or predictive speciation. Those remain separate M10 acceptance criteria.

## Inputs still needed

Request the per-run numerical wall-temperature profiles and coordinate origins,
raw iodine gamma scans, evaporated amounts/release histories, exact boat
placement and geometry, and flow-controller reference conditions from Jörg
Neuhausen and Lu Liu. For production Bi cases, use the external licensed data
only after provenance confirmation; the repository's PbI2 phase-resolved table
is already available. The proposed M10e case block (run/030–034) remains reserved
for qualified solver cases; this preliminary package does not allocate those
case numbers.
