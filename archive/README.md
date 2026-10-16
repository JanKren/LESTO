# archive/: raw simulation cases kept outside git

The case directories (meshes, time directories, logs, solver binaries) behind
tracked results whose own folders keep only small records.  Everything here
except this README is git-ignored (`/archive/*` in `.gitignore`).  About
20 GB.

Nothing here is needed to rebuild the manuscript: `manuscript/build.sh` and
`manuscript/check_results.py` read only tracked records.  These directories
are for re-plotting fields, extending a run, or auditing a number back to the
solver output.

## Where it came from

All of it was produced between 2026-10-09 and 2026-10-15 by the agents of the
review round (plan sections 38 and 39 of
`applications/rhoFixedFlowFoam/doc/phase-change-plan.md`), each in its own
task directory of a temporary session scratchpad:

    ~/.claude/tmp/claude-1000/-home-jan-LESTO/2b1263bd-4432-4d2e-922b-8c06c5c08659/scratchpad/sonnet-tasks/<task>/

They were moved here on 2026-10-16.  Logs, `results.json` files and `RESULTS.md`
files inside still name the old paths; the table below maps them.  The rules
the agents worked under are in `agent-rules-2026-10-09.md`; each task's brief
is the `TASK.md` (or `DESIGN.md`) in its directory.

Common to all Liu (2025) runs: OpenFOAM v2412, `rhoFixedFlowFoam` M10e,
executable sha256 `16fd804e593f39ed...` (a copy is in each `bin/`), 5 degree
wedge, R = 2.4 mm, L = 1.25 m, frozen expanding helium carrier from
`liuCarrierTool`, 2 ms steps to 60 s, serial runs pinned with `taskset`.

## Contents

| directory | old task | what it is | tracked record it backs |
|---|---|---|---|
| `liu2025-production-576x72/` | E | **Production runs**: the 15 baselines and the 4 silica alignment cases at 576 x 72 = 41,472 cells, `Gauss limitedLinear 1` (19 cases in `cases/`) | `run/024-liu2025/production_576x72_limitedLinear/`, `manuscript/data/fields/` (`PbI2_SS_Q1_baseline` is the field figure) |
| `liu2025-scheme-mesh-L1-L3/` | A | Scheme and mesh study: pure PbI2 and steel, levels L1-L3 under limitedLinear (`*_ll`); `steel_L4_*` are the prepared, unrun inputs that task F copied; `carriers/grid864x108/` the audited L4 carrier | `run/024-liu2025/scheme_mesh_study/` (`levels1-3_*`) |
| `liu2025-mesh-L4-and-vanLeer/` | F | Fourth mesh level 864 x 108 = 93,312 cells for the steel case under upwind and limitedLinear, and `Gauss vanLeer` on L1-L3 for both cases | `run/024-liu2025/scheme_mesh_study/` (`results.json`, `vanLeer_vs_limitedLinear.txt`) |
| `liu2025-production-sensitivities/` | I (part A) | Diffusivity (D/2, 2D) and accommodation (0.1, 0.01) variants of two production baselines, 576 x 72, limitedLinear; `ref/` and `tables/` are the input copies the cases point to | `run/024-liu2025/production_sensitivities/` |
| `liu2025-gems-engine-check/` | H | Same steel case (`LBE-I_SS_II_Q2_baseline_sat0.01`) with GEMS called per cell (`local8`), with frozen GEMS standard states (`frozen8`) and with the tabulated kernel (`kernel8`), 8 species correctors.  local vs frozen: bitwise identical wall profile; frozen vs kernel: relative L1 8.6e-5 | `run/024-liu2025/gems_engine_check/` |
| `lobresco2026-head-to-head/` | G | Five PbI2 experiments of Lobresco et al. (2026), T106_6/60/300, T5_60, T540_60, on the production wedge (`cases/prod`, 318 x 72) and the coarse wedge (`cases/coarse`, 212 x 48); `digitized/` the curves taken from the paper | `run/023-lobresco2026/` (inputs and `final_profiles` only) |
| `damkohler-first-pass/` | I (part B) | Notes, CSV and figures of the first Damköhler sweep.  Its cases were dropped: the final sweep is in `run/016-hks-channel-halfcell/work/damkohler/` | superseded by `manuscript/data/damkohler/` |
| `test-suite-run-2026-10-12/` | D | Logs of `tests/Alltest all` on HEAD f7c11cc (v2412, plus the v2606 pass) | plan section 38, item 4 |
| `case-ladder-t008/` | (none) | Timing run of case 008 | `presentations/case-ladder/runs/runs.json` names its old path |

Integrity: the five fields of `liu2025-production-576x72/cases/PbI2_SS_Q1_baseline`
read by the field figure still have the sha256 recorded in
`manuscript/data/fields/provenance.json` (checked after the move).

## Related data that lives elsewhere

- `~/lesto-submission-followup/` (4.6 GB): the follow-up package.  It holds the
  upwind L4 field cross-check (`cases/PbI2_SS_Q1_baseline__grid864x108`) and the
  8-corrector runs of `run/024-liu2025/iteration_check`.  Tracked
  `collect.py` scripts and READMEs name this path, so it stays where it is.
- `run/024-liu2025/work/` (git-ignored): prepared inputs and tables, and the
  18,432-cell upwind follow-up study (`followup-20261007/`, 26 GB).
- `run/016-hks-channel-halfcell/work/damkohler/`: the final Damköhler sweep.
- `applications/rhoFixedFlowFoam/tests/work/`: test and verification output,
  including the Graetz runs of `manuscript/data/graetz/`.
- `~/opt/gems-data/m10-design-2026-10-04/`: the GEMS maps and synthesis of the
  M10 design (plan section 36).

## What was not kept

The rest of the scratchpad was deleted on 2026-10-16:
- repository clones and source snapshots, whose code is now in git
  (`impl/`, `ladder/clone*`, `m7verify/clone`, `v9r1/`, ...);
- the M4/M5/M7/M9 prototype harnesses (`m7study*`, `m7int`), superseded by
  `tests/verification/` and the signed-off plan sections 20-35;
- `mirror/`, a scratch overlay for re-running the case ladder;
- `m10design/`, a duplicate of `~/opt/gems-data/m10-design-2026-10-04`;
- smoke tests and aborted runs, page renders, and one-off edit scripts.
