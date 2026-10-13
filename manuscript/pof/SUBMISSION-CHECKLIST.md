# Physics of Fluids submission checklist

Prepared 13 October 2026 for `manuscript/pof/lesto-pof.tex`. Nothing has
been submitted. Status: **done**, **placeholder** (red `[AUTHOR TO
CONFIRM: ...]` in the PDF), **blocked** (needs results or external input
first), **decision** (author's call).

Sources for AIP requirements: AIP Publishing author instructions,
<https://publishing.aip.org/resources/researchers/author-instructions/>,
and AIP's statement on AI language models (10 February 2023),
<https://publishing.aip.org/about/news/on-the-use-of-ai-language-models-in-scholarly-communications-at-aip-publishing/>,
both read on 13 October 2026. Re-check them before submitting.

## 0. Blockers (do not submit until these are cleared)

| # | Item | Why it blocks |
|---|------|---------------|
| B1 | Pending simulation campaigns (Section 3) | Numbers in the abstract, Sections 4–8, all generated tables and five figures will change. |
| B2 | Data are not yet public | `git branch -r --contains HEAD` is empty: the branch holding `run/034-liu2025/*` results and `manuscript/` is not on any remote. The Data Availability statement says "openly available"; it is only true after a push (author action) and ideally an archived release with a DOI (e.g. Zenodo). |
| B3 | Author declarations | CRediT roles, conflict-of-interest statement, acknowledgments/funding, corresponding e-mail: all placeholders. Required by AIP. |
| B4 | Figure text below 8 pt | In the PoF layout every legend renders at 5.3–6.6 pt (Section 5). AIP: "Legends or labels within figures should be a minimum of 8-point type size". |
| B5 | Co-author approval | Ničeno must see and approve the manuscript and the CRediT statement (AIP: corresponding author ensures all authors agree). `review-notes.md`: no messages have been sent. |

## 1. AIP / PoF requirements

| Requirement | Status | Notes |
|---|---|---|
| REVTeX 4-2 `aip,pof,reprint` | **done** | `pof/lesto-pof.tex`, `bash manuscript/pof/build.sh`. `--preprint` gives a single-column double-spaced copy. |
| Layout | **done** (note) | REVTeX's own `aip4-2.rtx` sets PoF to **single column at reprint stage** ("POF journal to be one column at Reprint stage"); the log says "pof journal style Single column". The two-column assumption does not apply to PoF. Shared sections still use `figure*`/`table*` and two-line equations, and a two-column AIP test build (`jap`) compiles with 0 overfull boxes, so a later switch is safe. |
| `\author`/`\affiliation` | **done** | Affiliation kept as in `main.tex`: "Paul Scherrer Institute, 5232 Villigen PSI, Switzerland". **decision**: PSI's current preferred affiliation string (e.g. centre/laboratory name). |
| Corresponding author e-mail | **placeholder** | `\email{...}` reads verbatim, so the placeholder is plain bracketed text, not red. |
| ORCIDs | **placeholder** | Entered in the submission system, not in the TeX. Not invented here. |
| Title | **done** (kept) | Current title kept; alternatives in Section 6. |
| Abstract ≤ 250 words, one paragraph, no citations/equations/footnotes | **done** | 242 words (math counted as one word; 245 counting units separately). No citations or displayed equations. PoF abstract is separate from the arXiv abstract; it adds the transport framing and keeps every number and caveat. |
| Keywords | **decision** | AIP's instructions state no keyword requirement; the TeX has none. If the submission form asks: thermal-gradient tube flow; vapor deposition; Hertz–Knudsen–Schrage; laminar mass transfer; finite-volume method; OpenFOAM; numerical verification. |
| Acknowledgments | **placeholder** | `acknowledgments` environment. The earlier draft's funding statements and grant numbers were removed pending confirmation; none invented. |
| Author Declarations: Conflict of Interest | **placeholder** | AIP standard text "The authors have no conflicts to disclose." inside the placeholder; confirm or replace with AIP's disclosure template. |
| Author Declarations: Author Contributions (CRediT) | **placeholder** | Roles must be chosen from the 14 CRediT roles (Conceptualization, Data curation, Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing – original draft, Writing – review & editing). Not assigned here. |
| Ethics approval | **done** (n/a) | No human or animal research. |
| Data Availability | **placeholder / blocked (B2)** | Based on the existing `main.tex` availability text, worded in AIP's "openly available in [repository]" form, plus the third-party/licence facts: digitised curves from Liu et al. under CC BY 4.0, external thermodynamic coefficient files not redistributed. The repository URL `https://github.com/JanKren/LESTO` answers HTTP 200 anonymously (public), but the relevant branch is not pushed. AIP's third-party template ("Restrictions apply to the availability of these data, which were used under license...") may be preferred for the coefficient files once their licence is known (requests-draft item 2). |
| Reference style | **done** | `aipnum4-2` (REVTeX default for `aip`), numbered by citation order, superscript in reprint. BibTeX's only warning, `jnrlst (dependency: not reversed)`, is REVTeX's internal control entry. |
| Cover letter | **done** (note) | AIP: "Cover letters are not accepted with submissions for Physics of Fluids." `cover-letter.md` is drafted as requested but must **not** be uploaded to PoF. |
| Alt text for every figure | **placeholder** | AIP: "Providing alt text is a requirement for publication." Seven figures need alt text, entered at submission. Not written. |
| Figures: format, fonts | **done** | Vector PDF, one embedded subset font (DejaVuSans, CID TrueType) per figure; no unembedded fonts. |
| Figures: size and minimum 8 pt | **blocked (B4)** | Section 5. |
| AI-use disclosure | **decision** | Section 7. |
| American spelling | **decision** | Text uses British spelling (modelling, vapour, digitised). AIP copy-editing normally standardises; decide whether to convert first. |

## 2. What was changed in shared files (affects the arXiv build too)

- `sections/*.tex`: all seven figures are now `figure*`; the three wide
  tables (experiments, study matrix, numerical) are `table*`. In the
  one-column arXiv layout both behave as before.
- Six equations are now set on two lines (`gathered`/`aligned`):
  eq:equilibrium, eq:halfcell, eq:balance, the helium property equation,
  the Ω(T*) line of the PbI2–He correlation (with `\nonumber`, so the
  numbering is unchanged), and the new eq:damkohler.
- Framing (no new results): one Introduction paragraph on the column as
  a laminar heated-tube mass-transfer problem; one Sec. 4.1 paragraph
  with Re, R²/D and a wall Damköhler number (Section 4); two sentences
  in the Conclusions.
- Citations added from entries already in `references.bib`:
  Hertz 1882, Knudsen 1915, Schrage 1953 (HKS); Shah & London 1978
  (Graetz-type duct mass transfer); Eichler 1975, Zvára 1985, Zvára 2008
  (thermochromatography); Walker et al. 1979, Batchelor & Shen 1985
  (thermophoretic deposition, listed as excluded); Rosner et al. 2000,
  Rosner & Arias-Zugasti 2007 (Soret transport of heavy vapours, listed
  as excluded); Roache 1994, Celik et al. 2008 (grid convergence). The
  reference list grows from 5 to 18.
- `references.bib`: capitals protected in two titles (`{H}eavy`,
  `{GEM-Selektor}`), Shah & London initials spaced ("R. K.", "A. L.").
- `main.tex`: `table*` given arxiv.sty's caption spacing. Abstract and
  everything else unchanged.

## 3. Pending simulation results that must update numbers before submission (B1)

Do not anticipate them; the current text describes the executed
first-order-upwind campaign only.

1. **Production limitedLinear runs on the 41,472-cell mesh.** If these
   become the reported results, the scheme statements (Sec. 4.2, 5.3:
   "first-order upwind", D_num ≈ |u_x|Δx/2) and every reported number
   change: abstract, Tables III–VI (and II if cases are added), Figs. 3–7, Conclusions, and
   `generated/*`. Cite the limiter (Sweby 1984 is in the bib; the
   OpenFOAM limitedLinear definition should be cited, e.g. Jasak's 1996
   thesis, also in the bib).
2. **Fourth mesh level and van Leer study.** Changes the three-grid
   statements (Sec. 5.2: 17 oscillatory / 4 nondecreasing / 1 roundoff /
   12 conditional), the 2.06–10.58 % fine/medium range, and possibly the
   "no mesh independence" scope. Add van Leer (1974) for the scheme; it
   is not in `references.bib`.
3. **8-versus-3 corrector check over a full run.** Sec. 4.2 and 5.1
   currently state that full-run iteration independence is
   "unqualified"; update only if the check supports it.
4. **Head-to-head on Lobresco et al.'s five PbI2 cases.** The paper
   currently says it makes no T-Flows comparison and that cross-code
   comparisons are "planned" (Sec. 5.1). If added, this becomes new
   results and a new section; the T-Flows case inputs requested in
   `requests-draft.md` item 3 may be needed.

After each update: `bash manuscript/build.sh --figures` (regenerates
tables and figures), `bash manuscript/pof/build.sh`, and
`python3 manuscript/check_results.py`. Recount the abstract (≤ 250 words).

## 4. Dimensionless numbers (added, Sec. 4.1)

Computed in the session scratchpad (`pof/dimensionless.py`, not part of
the repository) from inputs stated in the paper only: d = 4.8/5.0 mm,
ṁ = 2.727e-7 and 1.227e-7 kg/s, Sutherland A_s = 1.4683e-6,
T_s = 83.40 K, D_i = 1e-4 m²/s, α ∈ {1, 0.1, 0.01}, M(PbI2), M(BiI3),
T = 400 K (lowest printed peak label) to 1173.15 K (T_h).

| Quantity | Value in text | Computed |
|---|---|---|
| Re = 4ṁ/(πdμ), 100 mL/min steel | 1.5–3.0 | 1.54–2.98 |
| Re, 45 mL/min steel / silica | 0.7–1.3 | 0.69–1.34 / 0.67–1.29 |
| R²/D_i | 58 / 63 ms | 57.6 / 62.5 ms |
| Da, α = 1 | 1.4e3–2.9e3 | 1.44e3–2.90e3 |
| Da, α = 0.1 | about 76–150 | 75.7–153 |
| Da, α = 0.01 | about 7–15 | 7.23–14.6 |

**Decision (Péclet and Schmidt numbers).** These need the carrier density,
hence the column's absolute pressure, which the paper does not state.
The case files use `internalField uniform 101325` in `0/p`. If you add
"The column operates at 101,325 Pa" to Sec. 4.1, this sentence becomes
supportable: "At 101,325 Pa the baseline diffusivity gives Sc = 2.0–11.3
and a diameter-based Péclet number Pe = Re Sc of 5.9–17.4 at 100 mL/min
and 2.6–7.8 at 45 mL/min (400–1173 K); axial diffusion is therefore not
negligible next to advection." It was not added.

## 5. Figures (B4)

All seven PDFs are 6.62 in wide, drawn by `prepare_results.py` at
`figsize` width 6.5 in, 10 pt body text and 8 pt legends. The PoF
single-column text width is 5.48 in (396 pt), so LaTeX scales them:

| Figure | Width in PoF | Scale | 10 pt text | 8 pt legend |
|---|---|---|---|---|
| temperature_profiles | 0.88 × 5.48 in | 0.73 | 7.3 pt | 5.8 pt (legend default 10 pt → 7.3) |
| experiment_schematic | 0.94 × | 0.78 | 7.8 pt | – |
| refinement | 0.96 × | 0.79 | 7.9 pt | 6.4 pt |
| source_thermodynamics | 0.96 × | 0.79 | 7.9 pt | 6.4 pt |
| silica_alignment | 1.00 × | 0.83 | 8.3 pt | 6.6 pt |
| baseline_profiles | 1.00 × | 0.83 | 8.3 pt | 6.6 pt |
| physical_sensitivity | 0.80 × | 0.66 | 6.6 pt | – |

At AIP's one-column print width (3.37 in) the text would be about
5 pt, so none of these can be printed one-column as is.

Suggested changes to `prepare_results.py` (not edited here):

- Draw at the final width: `figsize` width 5.4 in (PoF text width) or
  3.37 in for the simple single-panel plots (temperature_profiles,
  physical_sensitivity, source_thermodynamics), and include them at
  `width=\linewidth` without scaling.
- `font.size` 9 and every `fontsize=8` legend at ≥ 8 pt *at final size*.
- Axis labels: "Temperature (°C)" instead of "degree C"; subscripts in
  legends and titles ("PbI$_2$", "BiI$_3$"); "5° wedge" in the schematic.
- baseline_profiles (6.92 in tall) and silica_alignment (6.59 in) fit the
  8.25 in maximum depth but will be the largest items on their pages.
- Keep `pdf.fonttype 42` (fonts are already embedded).

## 6. Title

Current title kept: "Conservative CFD modelling of iodine deposition in
thermal-gradient columns: numerical verification and input sensitivity".
It is accurate but reads as nuclear chemistry plus CFD; "CFD" in a PoF
title is common but generic. Alternatives (decision):

1. "Vapor deposition in laminar thermal-gradient tube flow: a
   conservative multispecies model, numerical verification and input
   sensitivity"
2. "Competition between radial diffusion and interfacial exchange in
   thermal-gradient deposition columns: conservative modeling of
   iodine transport"

## 7. AI-assisted writing (decision)

AIP's 10 February 2023 statement says that AI language models "should
not be listed as an author", and that their use "should be disclosed to
editors and reviewers — particularly if it's used to generate significant
amounts of text in the manuscript. Authors should provide this
information in the appropriate section of their manuscript and to the
editor with their submission." The current author instructions page,
read the same day, does not repeat a separate AI policy. The manuscript
contains no statement either way. Whether and where to disclose is the
authors' decision; if disclosed, the Acknowledgments or a methods note
are the usual places.

## 8. Open external inputs (`manuscript/requests-draft.md`; none sent)

1. Liu et al. group: original T(x) and scan coordinates (silica shift),
   controller reference conditions, source release history, boat
   geometry and position; enables the missing BiI3_SiO2_800 case.
2. HERACLES maintainers: Bi–I provenance and redistribution licence
   (needed for Data Availability and a public CFD input package), BiI3(l)
   record, BiI(g) value behind the +29 kJ/mol surrogate, Pb(g) record.
3. T-Flows authors: parameter sets, scheme and per-case wall tables;
   needed only if the head-to-head comparison (Section 3, item 4) is
   added.

## 9. References

- 18 cited entries; all articles have journal, volume, pages or article
  number, year and DOI. Books and chapters have publisher and year.
- **Incomplete:** `NASAGridConvergence` (NPARC Alliance web page) has
  `year = {n.d.}`; it prints "(n.d.), accessed 9 October 2026". AIP
  accepts web references but a year or "last modified" date is
  preferable; not guessed.
- `aipnum4-2` prints full journal names as stored ("Journal of
  Radioanalytical and Nuclear Chemistry"); AIP copy-editing abbreviates.
  DOIs are not printed by this style (expected).
- Uncited entries elsewhere in the bib have unspaced initials
  ("Martens, J.A.", "Piro, M.H.A."); fix if they become cited.
- **Suggested missing references** (not added; verify before citing):
  L. Graetz, Ann. Phys. (1883/1885), the original Graetz problem;
  R. Aris, Proc. R. Soc. A 235 (1956), Taylor–Aris dispersion, since
  Pe is O(1–10) (Taylor 1953 is in the bib, uncited);
  B. van Leer, J. Comput. Phys. 14 (1974), if the van Leer study is
  reported; a mixed-convection reference if buoyancy is discussed (below).

## 10. What a PoF reviewer is likely to criticise

1. **No spatial convergence and first-order upwind.** For a PoF methods
   paper, 2.06–10.58 % fine/medium profile changes, oscillatory
   three-grid sequences and numerical diffusion comparable to the
   physical diffusivity will draw the strongest objection. The pending
   limitedLinear, fourth-level and van Leer results matter most here.
2. **Constant diffusivity D = 1e-4 m²/s for all nine gases and all
   temperatures.** Real gas diffusivities scale roughly as T^1.75/p over
   400–1173 K (a factor of about 6), and Sc would be 2–11 at 1 atm.
   Only one variant uses a T- and p-dependent PbI2–He value.
3. **Uniform cross-sectional temperature and frozen carrier.** No
   conjugate heat transfer, no thermal entrance, no radial temperature
   profile; the carrier sees the wall temperature everywhere. The estimated
   wall-to-bulk temperature lag of developed flow is <= 2.6 K (Sec. 4.1),
   which supports the approximation away from the furnace entrance.
4. **Buoyancy: bounded (Sec. 4.1, `manuscript/estimates/`).** Gr <= 2.2,
   Gr/Re^2 0.15-0.78 in the 400-750 K deposition range and about 1.2 in the
   steepest cold silica section; compositional Gr/Re^2 <= 0.32. Stated as weak
   but not excluded in a horizontal column. The column orientation is still
   unknown: ask the experimenters (add to `manuscript/requests-draft.md`).
5. **Soret diffusion: bounded (Sec. 4.1).** With alpha_T = 1 as a generous
   bound, Soret drift <= 2.2 % of the carrier velocity and radial Soret flux
   <= 0.7 % of a Fickian wall flux.
6. **Axial diffusion and inlet species conditions.** Pe is O(1–10)
   (at 1 atm), so upstream diffusion from the source and the species
   inlet condition matter; these are not discussed.
7. **Fluid-dynamics novelty.** Most of the paper is chemistry, provenance
   and input reconstruction. Even with the framing, an editor may judge
   the fit marginal. Alternatives are a heat/mass-transfer or nuclear
   engineering journal.
8. **Validation is weak by design.** 60 s simulated against a 3 h
   experiment, digitised data, a silica coordinate inconsistency and
   deposit-conditioned source masses. The text is candid about all of
   these, which helps, but the comparison sections will be read as
   inconclusive.
9. **79 cases versus three experiments.** The paper already says these
   are computational cases, not replicates; keep that wording.
10. **Length.** About 21 pages in PoF reprint layout, with seven tables
    and seven figures; a reviewer may ask to move the property appendix
    and the case-matrix details to supplementary material.

## 11. Build record (13 October 2026)

| Build | Command | Pages | Overfull | Undefined refs | BibTeX |
|---|---|---|---|---|---|
| PoF reprint | `bash manuscript/pof/build.sh` | 21 | 0 | 0 | 1 benign (REVTeX `jnrlst`) |
| PoF preprint | `bash manuscript/pof/build.sh --preprint` | 33 | 0 | 0 | same |
| Two-column test (`jap`, not kept) | – | 16 | 0 | 0 | – |
| arXiv | `bash manuscript/build.sh` | 19 (was 17) | 0 | 0 | 0 |
| Consistency | `python3 manuscript/check_results.py` | PASS | | | rewrote `review-checks.json` |
