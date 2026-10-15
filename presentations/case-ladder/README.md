# The LESTO case ladder — PSI Beamer

A 32-page deck by Jan Kren, dated 15 October 2026: one slide per case of
`run/` (001–024), with a title page, an overview of the ladder (the
"based on" links), a slide on how every case is organised, four section
dividers and a closing slide. Cases 001–011 are Bojan Ničeno's; 012–024
are Jan Kren's. The slides and the section dividers say so.

Open `lesto-case-ladder.pdf`. Edit `main.tex` for the slides; every
content frame has a Beamer speaker note, extracted into `speaker-notes.txt`
by `package_deck.py`. The PSI theme (`beamer*PSI.sty`, `psi-assets/`,
`template-README.md`) is stored locally in this folder. It was originally
copied from the older project-review deck, which is no longer retained.

## Every case slide

Left: the case it is based on and its author, the one change, the key
dictionary lines, and the numbers to check. Right: a figure of the case's
result, its source, and how to run the case. The numbers come from the
case's `readme.md`, `notes.md`, dictionaries, `out`/`out_excerpt` logs,
`expect*.json`, `verification.json`, `qualification.json`, the tracked
records of 023 and 024, or from the runs made for this deck.

## Figures and runs

`make_figures.py` draws every figure into `figures/` (PDF for the deck, PNG
for viewing) and writes `figure_sources.json`, the files each figure reads
with their sha256, and `figures/values.tex`, the few numbers the slides
take from those runs. The figures read tracked repository files or
`runs/`, the small outputs of the runs made for this deck.

The runs were made on 15 October 2026 in a scratch copy of `run/` (cases
copied, `applications/` and `thermochemistry/` linked), never inside
`run/`, with OpenFOAM v2412 and the M10e solver (`rhoFixedFlowFoam`,
sha256 `16fd804e...`), pinned to cores 4, 5 and 11. `runs/runs.json` keeps
each command, its cores, start time, wall time and exit status.
`extract_runs.py <scratch>/run` copied the profiles, ledgers and the
foamToVTK output needed into `runs/`:

| case | run | figure |
|---|---|---|
| 001, 007 | `foamToVTK` of run/007/261 on the gmshToFoam mesh of run/001 | inlet mesh, carrier temperature |
| 008, 009, 011 | to t = 1 s (the cases' own 16 s is hours), 3 ranks | axial profiles |
| 010 | to t = 1 s, kept in `runs/010` but not shown: today's solver carries the coupled mock source of 011 in the same routine, so the run is not case 010; its figure uses the tracked `out_excerpt` | — |
| 012, 013 | `./Allrun steps 100` and the default excerpt check (PASS) | ledgers |
| 014, 016, 018, 019, 021 | `./Allrun` (PASS) | wall profiles |
| 022 | `./Allrun` (excerpt PASS) | the wall table (tracked) |

015 and 020 need the GEMS3K bridge build, and 017, 023 and 024 take
minutes to hours per case: their figures use the tracked records
(`out_excerpt`, `qualification.json`, `023/cases/*/final_profiles`,
`024/production_576x72_limitedLinear`). 001–006 use the tracked `out`
logs, 010 its `out_excerpt`. 018, 019 and 021 need the restricted Bi tables
(`LESTO_M10_DATA`, `LESTO_M10D_DATA`), which are not in the repository.

## Build

    bash build.sh              # Tectonic -> lesto-case-ladder.pdf (32 pages)
    bash build.sh --refresh    # make_figures.py first
    python3 package_deck.py    # speaker-notes.txt, manifest, source zip

`build.sh` uses `/home/jan/.local/opt/tectonic-0.17.0/tectonic` (or
`TECTONIC_BIN`), `XDG_CACHE_HOME=/tmp/lesto-beamer-cache`, and checks the
page count and reports overfull boxes. Build products stay in the ignored
`build/`.
