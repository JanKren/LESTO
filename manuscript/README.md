# Manuscript: CFD–GEMS coupling for volatile fission-product deposition

Follow-up to Lobresco et al., *Nucl. Eng. Des.* 459 (2026) 115201
(`../papers/Lobresco2026.pdf`), built on the OpenFOAM solver in
`../applications/rhoFixedFlowFoam`.

## Status

Skeleton and planning draft.

- `\todo{}` marks missing content.
- `\decide{}` marks open scientific or authorship decisions. They render in
  orange so they are hard to miss.

## Layout

    main.tex                 arxiv.sty preprint; title, authors, abstract, \input list
    sections/01-introduction.tex  ... 08-conclusions.tex, A-properties.tex
    references.bib           every entry DOI-checked against Crossref / doi.org
    literature-review.md     annotated review for the authors (not journal prose)
    figures/                 figure files (\graphicspath)
    arxiv.sty                https://github.com/kourgeorge/arxiv-style (MIT)

## Building

No LaTeX is installed on this workstation. Build on Overleaf, or with
`latexmk -pdf main.tex` on a machine that has TeX Live. Required packages:
natbib, cleveref, siunitx, mhchem, authblk, doi.

## Target journal

Not decided yet; the aim is above *Nucl. Eng. Des.* The arxiv template keeps
switching cheap: change the documentclass and the bibliography style.
