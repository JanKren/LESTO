# Conservative CFD modelling of iodine deposition

Working manuscript for author review, updated 9 October 2026. The draft now
reports the executed model and completed results: 15 baseline cases and
64 follow-up cases, representing three published experimental configurations.

Read [the compiled manuscript](lesto-iodine-methods.pdf). The scope is a
methods and input-sensitivity study. Experimental comparisons are conditional;
spatial convergence, release inputs and experimental coordinate registration
remain insufficient for independent predictive validation. See
[author review notes](review-notes.md) for the remaining decisions and evidence.
The [source and saved-results package](lesto-iodine-review.zip) is portable.

The manuscript contains seven figures, seven tables and a property/provenance
appendix. Five tables and all figures are generated from the saved completed
reports by `prepare_results.py`; `generated/results.json` records the values
and hashes of their inputs. No CFD cases are rerun or altered by this script.

## Build

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
`run/034-liu2025/{full_analysis,followup_results}` reports and profiles and
`thermochemistry/benchmarks/liu2025` scripts/CSV inputs. Compilation alone
does not require those inputs or the external thermodynamic coefficients.

After compilation, run `python3 manuscript/check_results.py` to check profile
arithmetic, conservation stoichiometry, hashes and TeX/PDF consistency. It
requires `pdfinfo` and `pdftotext` (Poppler) and writes `review-checks.json`.
Run `python3 manuscript/package_review.py` to rebuild the review ZIP.

## Files

- `main.tex`, `sections/`: manuscript source; no submission template selected.
- `figures/`: seven scientific figures, as vector PDFs and PNG previews.
- `generated/`: five result tables, numerical values and SHA-256 manifests.
- `prepare_results.py`: results-to-manuscript pipeline.
- `check_results.py`: independent profile arithmetic, input hashes and TeX checks.
- `package_review.py`: portable review archive and its SHA-256 manifest.
- `build.sh`: local or portable Tectonic compilation.
- `review-notes.md`: current scope and unresolved author/scientific issues.
- `notes-for-authors.md`, `literature-review.md`: historical planning material;
  unexecuted T-Flows comparisons and earlier hypotheses are not current results.
- `arxiv.sty`: existing MIT-licensed preprint style.

The source-and-results ZIP preserves repository-relative paths and includes
the saved profiles needed for figure regeneration. It omits source-paper PDFs,
raw CFD time directories and external thermodynamic coefficient files. The
package is for review and result reproduction; it is not a complete CFD input
distribution or a journal submission.
