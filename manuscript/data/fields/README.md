# PbI2 gas field of the pure PbI2 steel column: tracked records

The records behind Fig. `deposition_field` (manuscript Sec. 4.2, "Transport
regime").  The case directories and their fields stay outside git.

## Runs

Both are `PbI2_SS_Q1_baseline`, 2 ms to 60 s, solver sha256 `16fd804e...`.

| record | mesh | species advection | tracked wall profile of the same run |
|---|---|---|---|
| `production` (figure, `Field*` numbers) | 576 x 72 = 41,472 | `Gauss limitedLinear 1` | `run/024-liu2025/production_576x72_limitedLinear/profiles/PbI2_SS_Q1_baseline` |
| `upwind_L4` (cross-check, `FieldUp*`) | 864 x 108 = 93,312 | `Gauss upwind` | `run/024-liu2025/scheme_mesh_study/upwind_L4/profiles/PbI2_SS_Q1_baseline` |

The case directories at collection are recorded in `provenance.json`.  The
production case was run in a session scratch directory; its 60 s
`axialProfiles` file is byte-identical to the tracked production profile.
It now lives, git-ignored, in
`archive/liu2025-production-576x72/cases/` (see `archive/README.md`)
(fields re-verified against `provenance.json` on 2026-10-16).  The upwind
case is in `~/lesto-submission-followup/cases/`, outside the repository.

## Contents

- `pbi2_field_<record>.npz`
  - `x_cm`, `r_mm`: cell centres.
  - `window`, `mole_window`: the PbI2(g) mole fraction x = c R T / p for
    30-50 cm, at full radial resolution.
  - `wall_T_K`, `mole_axis`, `mole_wall_row`, `mole_cup`, `mole_eq_wall`:
    whole-column series:
    - the temperature;
    - the mole fraction on the axis and in the wall row;
    - the mole fraction flow-weighted over the section (weight u_x r);
    - the mole fraction at equilibrium at the wall temperature: the
      case's own `pv_PbI2_phases.csv`, minimum over the phases, divided
      by p.
  - `wall_mol_per_m`: the PbI2 wall inventory per unit length of the
    5 degree wedge, `mDep_PbI2_s / M * 2 R tan(2.5 deg)`.  This is the
    chord width of the blockMesh wedge, not R theta.
- `provenance.json`: per record, the case, solver sha256, mesh, scheme and
  the sha256 of every file read.  The wall inventory, doubled for iodine,
  equals the tracked element-I profile of the same run to round-off.
- `collect.py`: rebuilds a record, for example

      python3 manuscript/data/fields/collect.py upwind_L4 \
          ~/lesto-submission-followup/cases/PbI2_SS_Q1_baseline__grid864x108

`check_results.py` (`check_fields`) does three things for each record:
- checks the npz against its provenance;
- bins the wall inventory into the 1 cm bins of the tracked profile, with a
  relative L1 below 1e-12;
- recomputes the `Field*` and `FieldUp*` numbers of
  `generated/numbers.tex`.

Both runs show the same picture:
- the deposit begins where the flow-weighted mean first saturates;
- the wall cell stays at S <= 1.009;
- the core is supersaturated, at most 1.64/2.03 (production) and
  1.77/2.26 (upwind) for the mean/axis.
