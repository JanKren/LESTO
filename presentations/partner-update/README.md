# LESTO Task 3.4 partner update — PSI Beamer

"CFD–GEMS3K modelling of iodine transport and deposition — LESTO Task 3.4
progress", by Jan Kren with Bojan Ničeno, dated 15 October 2026. For the
LESTO project partners (reactor safety, radiochemistry; not OpenFOAM
developers): a 30-minute talk of 18 slides, followed by a backup section
(divider and two slides) that is not part of the 30 minutes.

Open `lesto-partner-update.pdf`. Edit `main.tex`; every slide has a Beamer
speaker note, extracted into `speaker-notes.txt` by `extract_notes.py`,
which also sums the minute budget given at the start of each note.

## Minute budget (main part, 30 min)

| page | slide | min |
|---|---|---|
| 1 | Title | 0.5 |
| 2 | From column experiments to the reactor hall (Task 3.4 goal, where the tool fits) | 2 |
| 3 | Two sets of published column experiments | 2 |
| 4 | Model, part 1: transport and wall exchange | 2 |
| 5 | Model, part 2: chemistry with GEMS3K | 2 |
| 6 | Built step by step: 24 test cases | 1 |
| 7 | Verified: wall law and conservation | 2 |
| 8 | Deposition is limited by diffusion (Damköhler sweep) | 1.5 |
| 9 | Inside the column: the wall stays saturated (PbI2 field) | 1.5 |
| 10 | What matters for the deposit (sensitivities) | 2 |
| 11 | Lobresco et al.: onset within 18 K | 2 |
| 12 | Liu et al.: the steel columns | 2 |
| 13 | Liu et al.: the silica column | 1.5 |
| 14 | Thermodynamic data decide the iodide | 1.5 |
| 15 | D38: the coupled tool and what it offers | 2 |
| 16 | What the model does not yet cover | 1.5 |
| 17 | Where partners could help | 2 |
| 18 | Next steps, thank you | 1 |
| 19–21 | Backup: divider; spatial and iteration accuracy; omitted effects (buoyancy, Soret) | — |

## Sources of the numbers

Every number on a slide comes from a tracked repository file:
`manuscript/generated/numbers.tex` and `results.json` (prepared by
`manuscript/prepare_results.py` and checked by `check_results.py`), the
generated manuscript tables, the manuscript sections, the run records of
`run/023-lobresco2026` and `run/024-liu2025`, and
`deliverables/D38/main.tex`. The partner descriptions on page 2 are the
one-line task descriptions agreed for this deck; nothing else is said
about partners' work, and no individuals are named.

## Figures

`make_figures.py` draws every figure into `figures/` at its size on the
slide (3.05 x 1.72 in for a 0.55 column, 5.6 in across for the Liu
profiles), so the 7–8 pt figure text is not scaled down. It reads only
tracked records and writes `figure_sources.json` with the sha256 of each
input:

| figure | slide | input |
|---|---|---|
| `temperature_profiles` | 3 | `thermochemistry/benchmarks/liu2025/digitized/*_temperature.csv` |
| `model_sketch` | 4 | `thermochemistry/benchmarks/liu2025/digitized/PbI2_SS_temperature.csv` |
| `graetz` | 7 | `manuscript/data/graetz/fits.json`, `reference.json` |
| `damkohler` | 8 | `manuscript/data/damkohler/damkohler_sweep.csv` |
| `field` | 9 | `manuscript/data/fields/pbi2_field_production.npz` (numbers on the slide: `manuscript/generated/numbers.tex`, `Field*`) |
| `sensitivity` | 10 | `run/024-liu2025/production_sensitivities/comparison.csv` |
| `lobresco_onset` | 11 | `manuscript/generated/results.json` (lobresco2026) |
| `liu_profiles` | 12 | `run/024-liu2025/production_576x72_limitedLinear/profiles/*/iodine-comparison.csv` |
| `silica` | 13 | the same, silica Q1 cases, and `results.json` (alignment) |
| `thermodynamics` | 14 | `results.json` (chemistry), production `analysis.json` |

Published data: Lobresco et al. (2026) digitised measurements and model
curves are CC BY 4.0 (attribution on the slide); Liu et al. (2025) scans
and temperatures are adapted with citation, as in the manuscript. No
HERACLES or Barin table is shown.

## Build

    bash build.sh              # Tectonic -> lesto-partner-update.pdf (20 pages)
    bash build.sh --refresh    # redraw the figures first
    python3 extract_notes.py   # speaker-notes.txt

The PSI theme (`beamer*PSI.sty`, `psi-assets/`, `template-README.md`) is
copied from `presentations/case-ladder`, so this folder is self-contained.
