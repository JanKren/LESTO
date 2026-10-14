# Case 017 — Physical carrier of a measured experiment (Lobresco T106_6)

This case is based on case 016.

The channels of cases 014-016 have an imposed plug flow and an imposed
temperature.  Here the carrier is that of a measured experiment: helium
flows at 106 mL/min (standard, 273.15 K and 101325 Pa) through the 4.8 mm
column of Lobresco et al. (2026), experiment T106_6, whose wall carries
the digitised temperature profile.  The flow is computed with
`rhoSimpleFoam` (expanding, laminar), projected to an exactly
divergence-free flux and audited, then frozen for the PbI2 transport.  The
wall law is that of case 016 (HKS with the half-cell resistance), and the
deposit at the end of the 6 s experiment is compared with the measurement.
It is one of the five experiments of case 023, on its coarse mesh.

## Differences from case 016

The carrier (made by `make_cases.py` of case 023; nothing tracked here):

    5-degree wedge of the 0.69 m column, radius 2.4 mm, 212 x 48 cells
    helium of case 007 (Sutherland, Cp 5193, perfect gas at 101325 Pa)
    INLET flowRateInletVelocity, massFlowRate 4.381729182205084e-09 kg/s
      (106 mL/min x 5/360); WALL T from inputs/wallT_106mLmin.csv
    rhoSimpleFoam to 3000 iterations; liuCarrierTool -project, -audit

What follows from it:
- `interfaceTemperature wall`: wall and wall-cell temperatures now differ.
- The inventory boat of case 013 returns, on the axis at the digitised
  arrow x = 0.13 m (0.125-0.135 m, radius 0.48 mm, `areaPerVolume 4000`),
  with the wedge's share of Table 1, 4.84e-6 x 5/360 = 6.7222e-8 mol,
  removed at 6 s.
- `diffusivityModel PbI2He` (as case 013) instead of the constant D.
- A zero-total-flux (Danckwerts) inlet for the gas: clean helium enters,
  vapour may diffuse upstream (`0/Y_PbI2_g`, codedMixed).
- `Gauss limitedLinear 1`, dt 2 ms, end time 6 s, outer loop `nCorr 10`,
  `tolerance 1e-5`, the final tolerance 1e-14 of case 013.
- Profiles in 1 cm bins centred on the gamma-scan points, `observable
  wallPlusGas` (the comparison at shutdown).

## Modified files

    0/Y_PbI2_g, 0/Y_PbI2_s                  (wedge patches, Danckwerts)
    constant/speciesTransportProperties     (that of case 013)
    constant/thermophysicalProperties       (that of case 007)
    constant/thermochemistryProperties
    constant/turbulenceProperties           (new, laminar)
    system/fvSolution                       (that of case 013)
    system/blockMeshDict, system/controlDict, system/fvSchemes
    physicalCarrier.py   (new: runs the tooling of case 023 into work/)
    (removed: 0/T, 0/U, 0/p templates, system/carrierDict,
     system/decomposeParDict)

These dictionaries are exactly those that `make_cases.py` of case 023
writes for its coarse T106_6; `./Allrun` checks that.

## Run

    ./Allrun            # carrier (about 5 min), 3000 steps (about 3 min
                        # serial), analyse.py and expect.json: about 8 min
    ./Allrun prod       # the same on 318 x 72 cells (expect-prod.json)
    ./Allclean

`LESTO_CPUS` pins the carrier steps (default: all cores).  Repository data
only.

## What to check

- `work/coarse/carriers/He106/carrier-audit.json`: `status PASS`; the
  largest cell flux defect is 1.2e-17 of the inflow.
- At 6 s (`work/coarse/case/comparison.json`):

| quantity | value |
|---|---|
| T_dep, 1 % of the peak of wall + gas from the inlet (solver, 1 cm bins) | 839.086 K at x = 0.35018 m |
| T_dep on the digitised wall curve, minus the measurement (814 K) | 839.734 K, +25.73 K |
| computed minus measured peak | -4 cm (0.37 m against 0.41 m) |
| wall / gas / sample, wedge | 6.59509e-8 / 1.27133e-9 / 0 mol |
| largest closure over all steps | 2.1e-16 of the inventory |

- The final profiles equal case 023's tracked
  `cases/coarse/T106_6/final_profiles/` (the check allows 1e-9; with the
  M10e solver, sha256 16fd804e..., they are identical byte for byte, and
  so is `comparison.json`).
- `./Allrun prod`: T_dep 839.112 K, +25.75 K; production and coarse agree
  within 0.03 K.

Not represented: the boat obstruction, the end of the experiment (flow
stop, withdrawal) and the non-dilute gas next to the boat (mole fraction
up to 0.36; the monitor reports it).  Case 023 has the five experiments
and the comparison with both published models.
