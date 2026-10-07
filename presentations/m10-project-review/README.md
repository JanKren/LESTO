# Pb–Bi–I project review — PSI Beamer

A first draft of a **20-slide / 20-page** presentation for Jan Kren, dated
7 October 2026. It uses the supplied `/home/jan/Desktop/PSI-Beamer.zip` theme,
with the original theme files, logos and background assets preserved.

Open `lesto-m10-review.pdf` for the presentation. Edit `main.tex` for slide
content; every slide includes a Beamer speaker note. `speaker-notes.txt`
provides those notes separately. `lesto-m10-review-source.zip` contains the
editable deck, figures and PSI theme for transfer to Overleaf or another
LaTeX installation.

## Narrative

1. Title
2. Project status and scientific question
3. Experimental setup
4. Available temperature profiles
5. Frozen-carrier modelling workflow
6. Transport and conservative phase change
7. Gas species and shared reactive wall inventories
8. Q1 and Q2 source assumptions
9. Verification and launch checks
10. Completed 15-case baseline matrix
11. Q1 peak-temperature comparisons
12. Steel deposition profiles
13. Silica input discrepancy
14. BiI thermodynamic and source sensitivity
15. Limits of the earlier mixed refinement comparison
16. Independent spatial and temporal refinement
17. The 64-case follow-up matrix and dated status
18. Exposure duration and flow-stop scope
19. Decisions and remaining work
20. Takeaways and primary references

## Build and update

Run from this directory:

```bash
bash build.sh
```

The script uses `latexmk`, `pdflatex` or Tectonic and checks for exactly
20 PDF pages. Automatic section-divider slides and note pages are disabled.
The PDF uses the supplied theme's Helvetica-class fallback (Liberation Sans
under XeTeX); it does not require Microsoft's Aptos font. All PDF fonts are
embedded. Build logs and intermediate LaTeX files stay in ignored `build/`.

Inside the LESTO checkout, refresh the figures and frozen execution counts
explicitly, then rebuild:

```bash
bash build.sh --refresh
```

This requires NumPy and Matplotlib and reads the saved baseline analysis,
published-curve digitisation, and local follow-up state. It does not change
the solver, input tables or running simulations. The standalone source ZIP
already includes all plotted figures and `status_snapshot.tex`, so it can be
compiled without access to CFD output or the external thermodynamic database.
After editing or refreshing, regenerate the transfer bundle with:

```bash
python3 package_deck.py
```

## Evidence and scope

- Numerical result slides use the completed 15-case analysis in
  `run/034-liu2025/full_analysis/`.
- Published scans retain their digitised coordinates and heights. Curves are
  redrawn as vector graphics; no experimental alignment or parameter fitting
  is performed. Liu et al. (2025) is credited on the relevant slides.
- The flow-stop simulations and three-hour analytical projections are
  explicitly distinguished from experimental cooling and native three-hour
  CFD runs.
- The silica CFD uses a common 4.8 mm empty-tube carrier; the paper's silica
  tube is 5.0 mm. Boat obstruction and surface-specific kinetics remain
  provisional.
- `status_snapshot.json` freezes the execution counts and their timestamp.
  These counts are not live in the PDF. Unfinished follow-up scientific
  results are not presented.
- `figure_sources.json` fingerprints numerical and digitised inputs.
  `deck_manifest.json` fingerprints the prepared PDF and source bundle.
- The theme and PSI graphics retain their ownership and usage terms in
  `template-README.md`. No external evaluated thermodynamic coefficients,
  raw CFD fields, solver binaries or experimental PDFs are included.

Before an external talk, update the completed follow-up analysis, check the
intended audience and presentation title, and review the scientific claims
with the project collaborators.
