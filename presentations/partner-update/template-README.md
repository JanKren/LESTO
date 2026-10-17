# PSI Beamer Theme

A LaTeX **Beamer** theme reproducing the corporate design of
**PSI – Center for Scientific Computing, Theory and Data**, derived from the
official PowerPoint template (`VO-9722-867_0_SCD Presentation Content Slides.potx`).

![preview](preview.png)

It reproduces the template's title slides (blue / dark‑blue / red particle
backgrounds + a white variant), the coloured section dividers, the content
layouts, the PSI logo lockups and the page‑number / institute / date footer,
using the real PSI logos and background images extracted from the `.potx`.

---

## Quick start

```bash
pdflatex main
pdflatex main      # run twice: section pages & the logo overlay use references
```

or simply `latexmk -pdf main`. The theme also compiles with **XeLaTeX** and
**LuaLaTeX** (see *Fonts* below).

A minimal document:

```latex
\documentclass[aspectratio=169]{beamer}
\usetheme{PSI}

\title{Your talk title}
\subtitle{An optional subtitle}
\author{Your Name}
\institute{Villigen}          % shown as the "location" on the title slide
\date{\today}

\psititlestyle{blue}          % blue | darkblue | red | white

\begin{document}

% Title slide (plain frame, no header/footer)
{\setbeamertemplate{footline}{}\setbeamertemplate{headline}{}
\begin{frame}[plain,noframenumbering]\titlepage\end{frame}}

\section{Introduction}         % automatically shows a coloured section page
\begin{frame}{A content slide}
  Your content here.
\end{frame}

\end{document}
```

---

## Files

| File | Purpose |
|------|---------|
| `beamerthemePSI.sty` | main theme (`\usetheme{PSI}`); table & figure helpers |
| `beamercolorthemePSI.sty` | the PSI colour palette and element colours |
| `beamerfontthemePSI.sty` | fonts (Fira Sans lookalike / Aptos) |
| `beamerouterthemePSI.sty` | header logo, footer, frame title |
| `beamerinnerthemePSI.sty` | title page, section pages, blocks, itemize, TOC |
| `main.tex` | a demonstration deck showing every layout |
| `psi-assets/` | PSI logos and particle backgrounds (see below) |

Keep the `.sty` files and the `psi-assets/` folder next to your `.tex` file
(or install them into your local `texmf` tree).

---

## The colour palette

| Name | Hex | Role |
|------|-----|------|
| `PSIblue`    | `#0014E6` | primary accent (default) |
| `PSIgreen`   | `#00F0A0` | accent 2 |
| `PSIred`     | `#DC005A` | accent 3 |
| `PSIpurple`  | `#6E14DC` | accent 4 |
| `PSIyellow`  | `#F0F500` | accent 5 |
| `PSImagenta` | `#F050FA` | accent 6 |
| `PSIdarkblue`| `#000073` | dark‑blue title background |
| `PSIgray`    | `#4B4B4B` | body / footer text |
| `PSIlightgray`| `#B9B9B9` | light grey |

Re‑skin the **whole deck** to a different accent with one line in your preamble:

```latex
\colorlet{PSIaccent}{PSIred}   % bullets, blocks, table tints, TOC numbers, ...
```

---

## Customisation

**Title‑slide style** — `\psititlestyle{blue|darkblue|red|white}` before the
title frame. `blue`/`darkblue`/`red` use the full‑bleed particle backgrounds;
`white` gives a white slide with a coloured bottom band (uses `PSIaccent`).

**Section pages** — a coloured section page is shown automatically at every
`\section{...}`. Pick its colour just before the `\section`:

```latex
\psisectionstyle{green}   % blue | darkblue | red | purple | green | yellow | magenta
\section{Results}
```

**Footer** — the institute line and date are macros:

```latex
\renewcommand{\PSIfootername}{PSI Center for Scientific Computing, Theory and Data}
\renewcommand{\PSIfootdate}{\today}   % default is Swiss-style DD.MM.YYYY
```

**Bold in‑slide headings** (the template's "Heading (formatted bold)"):

```latex
\psiheading{Incompressible Navier--Stokes}
Body text follows on the next line.
```

**Tables** — call `\psitablecolors{<accent>}` right before `\begin{tabular}`,
start the header row with `\rowcolor{<accent>}`, and wrap header cells in
`\psihead{...}`:

```latex
\psitablecolors{PSIblue}
\begin{tabular}{@{}lll@{}}
  \rowcolor{PSIblue}\psihead{Method} & \psihead{Order} & \psihead{Cost} \\
  Finite volume & 2 & low \\
  Spectral      & $\infty$ & medium \\
\end{tabular}
```

Swap `PSIblue` for `PSIpurple` or `PSIred` for the other table styles.

**Figures** — captioned figures and a photo placeholder:

```latex
\psifigure[0.8\linewidth]{myplot.pdf}{A caption}   % by width
\psifigureh{2cm}{myplot.pdf}{A caption}            % by fixed height (for grids)
\psiplaceholder{5cm}{3cm}{drop figure here}        % grey stand-in box
```

**Layout metrics** — the 5.7 % margins, logo size, footer baseline, etc. are
`\newcommand`s at the top of `beamerouterthemePSI.sty`; edit them to nudge the
whole grid.

---

## Fonts

The PSI template uses **Aptos** (a Microsoft Office font that is not freely
redistributable for LaTeX). This theme therefore ships configured for
**Fira Sans** — a free, humanist grotesque that closely matches Aptos and
compiles under plain **pdfLaTeX**. If Fira Sans is not installed it falls back
to Helvetica (Nimbus Sans), the template's own listed Arial‑class fallback.

To use the **real Aptos**, install it on your system and compile with
XeLaTeX or LuaLaTeX — the font theme picks Aptos up automatically. To force a
different font, edit the one relevant line in `beamerfontthemePSI.sty`.

---

## Assets (`psi-assets/`)

| File | Content |
|------|---------|
| `psi-full-white.png` / `psi-full-black.png` | full lockup: sphere + **PSI** + tagline |
| `psi-mark-white.png` / `psi-mark-black.png` | compact mark: sphere + **PSI** |
| `psi-bg-blue.png` / `psi-bg-darkblue.png` / `psi-bg-red.png` | particle title backgrounds |

These are the original PSI assets extracted from the `.potx`. They remain PSI
property — use them only for genuine PSI material, in line with PSI's corporate
design guidelines.
