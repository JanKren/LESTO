# Case 023 — notes

The readme of this case as it stood before the renumbering of the case
ladder (plan section 40) (case 035 until then), with the paths updated.  `readme.md` is
the short version.  Original title: Head-to-head on Lobresco et al. (2026): five PbI2 cases (plan section 39)

Our single-species HKS model run, without any parameter adjusted to these
data, on the five PbI2 thermosublimatography experiments of Lobresco et al.,
Nucl. Eng. Des. 459 (2026) 115201: T106_6, T106_60, T106_300, T5_60 and
T540_60 (flow in standard mL/min, duration in s).  The digitised measurements,
wall temperatures and the paper's own model curves are in
`thermochemistry/benchmarks/lobresco2026` (CC BY 4.0, attribution there).
The comparison is reported in Sec. 6.3 of the manuscript
(`manuscript/prepare_results.py` builds its table and figure from this folder).

## Set-up (what is given and what is assumed)

| input | value | status |
|---|---|---|
| flow, duration, sample amount | Table 1 of the paper | given |
| flow reference | 273.15 K, 101,325 Pa (paper, standard flow) | given |
| tube | 4.8 mm inner diameter, 0-0.69 m from the figure origin, 5-degree wedge | given / reconstructed |
| wall temperature | `inputs/wallT_<flow>mLmin.csv`, digitised curves; 106 mL/min is the mean of the T106_6 and T106_300 panels (`make_inputs.py`) | reconstructed |
| sample | inventory of the Table 1 moles x 5/360, a 1 cm long, 0.48 mm radius cylinder on the axis at the digitised arrow (0.13 m at 106, 0.15 m at 5 and 540 mL/min), evaporating by the same HKS law, removed at the end time | position reconstructed, geometry assumed |
| carrier | expanding laminar helium (`rhoSimpleFoam`, run/007 properties), audited with `liuCarrierTool` to 1e-12 of the inflow | ours |
| wall law | `model HKS`, accommodation 1, Ce 1, `wallResistance halfCell` | ours, untuned |
| PbI2 vapour pressure | `thermochemistry/systems/data/pv_PbI2_phases.csv` (phase-resolved) | ours, not the paper's |
| PbI2-He diffusivity | `diffusivityModel PbI2He` (the paper's correlation, as run/013) | as paper |
| inlet species condition | zero total flux (Danckwerts): clean helium enters, vapour may diffuse upstream | ours |
| advection, time step | `Gauss limitedLinear 1`; 2 ms, 0.5 ms for T540_60 | ours |
| meshes | `prod` 318 x 72 = 22,896 cells (the spacing of the 41,472-cell Liu mesh), `coarse` 212 x 48 = 10,176 cells | ours |

Not represented: the boat obstruction, the non-dilute source region (the
computed PbI2 mole fraction next to the sample reaches 0.33-0.36), and the
end of the experiment (flow stop, sample withdrawal, column pulled out of the
furnace).  The deposit is taken at the end time; `observable wallPlusGas`
also counts the vapour still in the tube where it is.

Solver: `rhoFixedFlowFoam`, sha256
`16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58` (the M10e
binary of run/024, unchanged).

## Files

- `make_inputs.py`, `inputs/`: wall-temperature inputs from the digitised curves.
- `make_cases.py`: builds and audits the carriers and writes the cases
  (`--levels prod coarse`); expects `bin/rhoFixedFlowFoam` and
  `bin/liuCarrierTool` (build from `applications/rhoFixedFlowFoam` and
  `thermochemistry/benchmarks/liu2025/carrierTool`, identical source) and
  OpenFOAM v2412.  `LESTO_CPUS` sets the taskset core list.
- `launch_full.sh <level> <first core>`: one serial run per case and core.
- `analyse.py <case> ...`: the paper's 1 %-of-peak rule on the wall-plus-gas
  observable, peak shift, L1 against the experiment, `comparison.json`.
- `cases/<level>/<case>/`: the dictionaries of each case (`system/`,
  `constant/` without mesh, wall table or vapour-pressure copy, `0/Y_PbI2_*`),
  `case.json` (inputs, carrier audit, hashes; paths made relative),
  `comparison.json` (as written by `analyse.py`), `ledger_summary.json`
  (largest closure and solver defect over every step, mole-fraction and
  negative-mass monitors) and the final-time profiles
  `final_profiles/PbI2_g_{gamma,native}.dat` (copied from
  `postProcessing/axialProfiles/<end time>/`, whose name is in
  `final_profiles/time`; numbered directories are git-ignored).
  Meshes, carrier fields, time directories and logs are not kept.

## Results (production mesh; manuscript Table V)

Deposition-onset temperature of the deposit minus the measurement (Table 2
of the paper), main deposit downstream of the sample:
+26, +19, +16, -11 and -14 K (T106_6, T106_60, T106_300, T5_60, T540_60),
RMS 18 K; the paper's temperature-based and vapour-pressure-based models give
RMS 35 K each.  The coarse mesh agrees within 0.3 K.  At 5 mL/min vapour
diffuses upstream and deposits near the cold inlet; applied from the inlet
the rule finds that small deposit (548 K), so it is applied to the main
deposit.  `comparison.json` uses the wall-plus-gas observable from the inlet
and therefore reports that upstream value for T5_60.  Computed peaks lie 1-4 cm
upstream of the measured ones.  Largest |closure| over all steps and runs:
3.3e-16 of the initial inventory.
