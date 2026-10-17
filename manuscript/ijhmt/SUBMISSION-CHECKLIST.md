# International Journal of Heat and Mass Transfer: submission checklist

Created 16 October 2026, when IJHMT replaced Physics of Fluids as the
target journal (reason: the paper's content is laminar internal mass
transfer with a nonlinear wall law, the readership that will use it, and
a "larger fish than NED").  Nothing has been submitted.  `../pof/` is kept
as a fallback; both builds share `../sections/`, including the abstract
(`abstract-journal.tex`) and the data statement
(`data-availability-journal.tex`).

**Sources.**  Elsevier's official guide for authors,
<https://www.sciencedirect.com/journal/international-journal-of-heat-and-mass-transfer/publish/guide-for-authors>,
refused automated access (HTTP 403) on 16 October 2026; only search
snippets of it were read.  Items marked *(unverified)* come from
third-party summaries and **must be checked on the official page before
submitting**.

Status: **done**, **placeholder** (red `[AUTHOR TO CONFIRM: ...]` in the
PDF), **check** (verify on the official guide), **decision** (authors).

## Build

| Command | Output |
|---|---|
| `bash manuscript/ijhmt/build.sh` | **working draft** `lesto-ijhmt.pdf`: `lesto-ijhmt.tex`, the arxiv.sty layout of `../main.tex` with line numbers (26 pages, 0 overfull boxes), and `lesto-ijhmt-supplementary.pdf` (6 pages) |
| `bash manuscript/ijhmt/build.sh --elsarticle` | `lesto-ijhmt-elsarticle.pdf`: Elsevier `elsarticle` review copy (12 pt, 1.5 spacing, line numbers, 1 in margins; 41 pages), for the submission |
| `bash manuscript/ijhmt/build.sh --print` | `lesto-ijhmt-print.pdf`: Elsevier `final,5p,twocolumn`, to estimate the printed length (23 pages); not for upload |
| `python3 manuscript/check_results.py` | also checks the IJHMT logs |

Both wrappers read the same IJHMT-specific parts: `nomenclature.tex`,
`declarations.tex` (Elsevier declarations, data statement, supplement
note), the shared journal abstract and the shared sections; only the
title block, layout and bibliography style differ (`unsrtnat` in the
draft, `elsarticle-num` for Elsevier).  Every build checks the abstract
(<= 250 words) and the highlights (3-5, <= 85 characters).  The `elsarticle`
builds keep one 2.6 pt overfull box while page 1 is shipped out: it comes
from the class (it stays with the e-mail, footnote, line numbers, spacing
and font size each removed) and is invisible in the PDF; the checker
allows exactly that one there and none in the draft.  The draft patches
`lineno` so that lines around amsmath displays are numbered.

## Blockers (as for PoF)

| # | Item |
|---|---|
| B1 | Author declarations: CRediT roles, competing interest, funding and acknowledgements, corresponding e-mail (all placeholders). |
| B2 | Data availability: the branch is public on GitHub since 16 October; an archived release with a DOI (e.g. Zenodo) and the HERACLES redistribution licence are still open. |
| B3 | Ničeno's approval of the manuscript, the CRediT statement and the change of journal. |
| B4 | **Length (check).** A third-party summary gives a limit of 40 "equivalent pages" (1.5 spacing, 12 pt, each figure or table one page) *(unverified)*.  Estimate: about 25 text pages (966 numbered lines at about 38 per page), about 5 pages of front matter, nomenclature and references, and 18 floats (9 figures, 9 tables): **about 48**.  If the limit holds, move about 8 floats to the supplement; candidates are decided by the authors. |

## Elsevier / IJHMT requirements

| Requirement | Status | Notes |
|---|---|---|
| Class | **done** | Submission: `lesto-ijhmt-elsarticle.tex`, `elsarticle` `preprint,12pt`, numbered references `elsarticle-num` (`numbers,sort&compress`).  Working draft until then: `lesto-ijhmt.tex` (arxiv.sty, line numbers). |
| Title | **done** | The PoF title (American "modeling"); the body uses British spelling (decision below). |
| Abstract <= 250 words, no references | **done** | 250 words, shared with PoF.  A third-party summary suggests "about 200 words" *(unverified)*; the official snippet says 250. |
| Highlights, 3-5 bullets, <= 85 characters each | **done** (check) | `highlights.txt`, uploaded as a separate file; 76-85 characters. |
| Keywords | **done** (check the maximum, usually 1-7) | 7 keywords in the `keyword` environment. |
| Nomenclature | **done** | Customary in IJHMT; list after the keywords.  Only $R$ carries two meanings (see below). |
| Graphical abstract | **decision** | Encouraged, optional *(unverified)*.  Candidates: the Damköhler sweep panel (a) or the geometry schematic. |
| CRediT authorship statement | **placeholder** | 14 CRediT roles; not assigned here. |
| Declaration of competing interest | **placeholder** | Elsevier's standard wording inside the placeholder. |
| Declaration of generative AI use | **placeholder / decision** | Elsevier requires a statement at the end of the manuscript if such tools were used, under the heading "Declaration of generative AI and AI-assisted technologies in the writing process"; template wording is in the placeholder.  Remove the section if none were used. |
| Funding, acknowledgements | **placeholder** | Grant numbers to be confirmed (EU LESTO). |
| Data availability | **placeholder (B2)** | Shared statement. |
| Supplementary material | **done** | `lesto-ijhmt-supplementary.pdf` (the shared supplement), uploaded separately. |
| Cover letter | **done** (draft) | `cover-letter.md`; IJHMT accepts one (PoF does not). |
| ORCID iDs | **placeholder** | Entered in the submission system. |
| Suggested reviewers | **decision** | None proposed here. |
| Review model | note | Single-anonymised *(unverified)*: names stay on the manuscript. |
| Open access | **decision** | Hybrid; Elsevier quotes the article charge at submission.  Check whether PSI or the EU grant (Horizon open-access obligation) covers it; an arXiv copy gives green open access either way. |

## Shared-text issues that IJHMT readers will notice (decision)

- **Symbols**: resolved on 16 October 2026.  The Graetz decay rate is now
  $\kappa$ ($\mu$ stays the viscosity) and the Graetz Péclet number
  $\mathrm{Pe}_R=U_bR/D$ (Pe stays diameter-based, Re Sc), in
  `sections/05-verification.tex`, the supplement and the generated full
  Graetz table.  $R$ remains both the gas constant and the tube radius;
  the text says which is meant where both appear.
- **Spelling**: British in the body, American in the title.  Elsevier
  copy-editing accepts either if consistent; decide.
- **Framing**: the introduction already presents the column as a laminar
  heated-tube mass-transfer problem.  For IJHMT the Damköhler sweep and the
  Graetz-Robin verification are the heat/mass-transfer core; reviewers
  will expect them early and prominently.

## Changes made to shared files on 16 October 2026

- `sections/abstract-journal.tex`, `sections/data-availability-journal.tex`:
  moved verbatim out of `pof/lesto-pof.tex` (PoF output unchanged: 28
  pages, 0 overfull boxes).
- `sections/07-multispecies.tex`: Table VIII (production-scheme
  sensitivities) is now `table*`, like the other wide tables; it was 46 pt
  too wide in a two-column layout.  The arXiv and PoF builds are unchanged
  in layout (0 overfull boxes).
- `check_results.py`: checks the IJHMT logs.
