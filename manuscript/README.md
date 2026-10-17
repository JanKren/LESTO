# Conservative CFD modelling of iodine deposition

Working manuscript for author review, updated 14 October 2026. The draft
reports the executed model and completed results for the three published
experimental configurations of Liu et al.:

- 15 baseline and 64 follow-up cases (18,432 cells, upwind);
- 19 production cases (the 15 baselines and the 4 silica temperature-coordinate
  cases on 41,472 cells with `Gauss limitedLinear 1`);
- a scheme and mesh study of two cases (upwind, limitedLinear and van Leer on
  8,192-41,472 cells, one case also on 93,312 cells; 12 runs beyond the
  existing upwind levels and the production runs);
- a full-run iteration check (3 cases, eight instead of three correctors);

and, separately, 10 runs of the five single-species PbI2 experiments of
Lobresco et al. (2026) on two meshes, compared with their measurements and
their two published models.  A pending-sign-off summary is in section 39 of
the phase-change plan.

Read [the compiled manuscript](lesto-iodine-methods.pdf). The scope is a
methods and input-sensitivity study. Experimental comparisons are conditional;
spatial convergence, release inputs and experimental coordinate registration
remain insufficient for independent predictive validation. See
[author review notes](review-notes.md) for the remaining decisions and evidence.
The [source and saved-results package](lesto-iodine-review.zip) is portable.

The manuscript contains eight figures, nine tables and a property/provenance
appendix. Six tables and all figures are generated from the saved completed
reports by `prepare_results.py`; a seventh generated table comes from
`assess_resolution.py`; `generated/results.json` records the values
and hashes of their inputs. No CFD cases are rerun or altered by this script.

## Build

The 12 October scientific reassessment adds a mass-preserving bin-width
audit and a matched-resolution thermodynamic comparison, using the original
79 CFD cases. These analyses concern the original upwind campaign; the
14 October production, scheme/mesh and iteration results are also retained.
`assess_resolution.py` reproduces those calculations and
`test_resolution.py` checks them against analytical examples. Build now
requires Python and NumPy even without figure regeneration. The precise
remaining dependencies and next numerical study are recorded in
`scientific-gap-register.md` and `next-study.json`; their current evidence
notes distinguish completed checks from the original study specification.
No CFD run was launched while merging these results.

From the repository root:

```bash
bash manuscript/build.sh
```

This compiles the existing TeX, generated tables and vector PDF figures with
Tectonic. Set `TECTONIC_BIN` if it is outside PATH. The first build may need
network access to download TeX packages. Intermediate files go into `build/`;
the review PDF is `lesto-iodine-methods.pdf`. On another TeX installation,
`cd manuscript` and `latexmk -pdf main.tex` also builds the article.

To regenerate figures, tables and numerical provenance before compiling:

```bash
bash manuscript/build.sh --figures
```

Regeneration requires Python, NumPy and Matplotlib, plus the accompanying
`run/024-liu2025/{full_analysis,followup_results,production_576x72_limitedLinear,scheme_mesh_study,iteration_check}`
reports and profiles, `run/023-lobresco2026` (case set-ups, final profiles),
and the `thermochemistry/benchmarks/{liu2025,lobresco2026}` scripts and
digitised CSV inputs. Compilation alone
does not require those inputs or the external thermodynamic coefficients.

After compilation (both builds), run `python3 manuscript/check_results.py` to
check profile arithmetic, conservation stoichiometry, the scheme/mesh,
iteration and Lobresco tables (recomputed from the saved profiles), hashes and
TeX/PDF consistency. It
requires `pdfinfo` and `pdftotext` (Poppler) and writes `review-checks.json`.
Run `python3 manuscript/package_review.py` to rebuild the review ZIP.

## Files

- `main.tex`, `sections/`: manuscript source; no submission template selected.
- `figures/`: eight scientific figures, as vector PDFs and PNG previews. They
  are drawn at their printed size (5.4 in wide or less, text 8-9 pt) and
  included without scaling.
- `generated/`: seven result tables, `numbers.tex` (macros used in the text),
  numerical values and SHA-256 manifests.
- `ijhmt/`: **target journal** (since 16 October 2026): International
  Journal of Heat and Mass Transfer (Elsevier `elsarticle`) version sharing
  `sections/`, `generated/` and `figures/`.  `bash manuscript/ijhmt/build.sh`
  builds the working draft (arxiv.sty with line numbers); `--elsarticle`
  the Elsevier copy for submission.  Highlights, cover letter and
  `SUBMISSION-CHECKLIST.md` are in the folder.
  The journal abstract and data statement are `sections/*-journal.tex`.
  IJHMT is the only target journal (17 October 2026); the earlier Physics
  of Fluids version was removed and stays in the git history (commit
  `bbceaf7`).
- `prepare_results.py`: results-to-manuscript pipeline.
- `check_results.py`: independent profile arithmetic, input hashes and TeX checks.
- `package_review.py`: portable review archive and its SHA-256 manifest.
- `build.sh`: local or portable Tectonic compilation.
- `review-notes.md`: current scope and unresolved author/scientific issues.
- `notes-for-authors.md`, `literature-review.md`: historical planning material;
  earlier hypotheses there are not current results.  The head-to-head on
  Lobresco et al. is now reported (Sec. 6.3); it compares with published
  measurements and published model results, not with T-Flows runs.
- `arxiv.sty`: existing MIT-licensed preprint style.

The source-and-results ZIP preserves repository-relative paths and includes
the saved profiles needed for figure regeneration. It omits source-paper PDFs,
raw CFD time directories and external thermodynamic coefficient files. The
package is for review and result reproduction; it is not a complete CFD input
distribution or a journal submission.
