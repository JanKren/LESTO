# Physics of Fluids submission checklist

Prepared 13 October 2026 for `manuscript/pof/lesto-pof.tex`; updated
14 October 2026 after integrating the completed simulation campaigns
(production runs, scheme/mesh study, iteration check, Lobresco
comparison). Nothing has been submitted. Status: **done**, **placeholder** (red `[AUTHOR TO
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
| B1 | ~~Pending simulation campaigns (Section 3)~~ **resolved 14 Oct** | All four campaigns are integrated (abstract, Secs. 1 and 4–8, Tables III–VII, Figs. 3–7). Remaining: the researcher's sign-off of plan section 39. |
| B2 | Data are not yet public | `git branch -r --contains HEAD` is empty: the branch holding `run/024-liu2025/*` results and `manuscript/` is not on any remote. The Data Availability statement says "openly available"; it is only true after a push (author action) and ideally an archived release with a DOI (e.g. Zenodo). |
| B3 | Author declarations | CRediT roles, conflict-of-interest statement, acknowledgments/funding, corresponding e-mail: all placeholders. Required by AIP. |
| B4 | ~~Figure text below 8 pt~~ **resolved 14 Oct** | Figures are drawn at printed size (≤ 5.4 in) and included unscaled; text 8–9 pt (Section 5). Only mathtext sub/superscripts are smaller. |
| B5 | Co-author approval | Ničeno must see and approve the manuscript and the CRediT statement (AIP: corresponding author ensures all authors agree). `review-notes.md`: no messages have been sent. |
| B6 | ~~Running studies (review-A, review-B)~~ **resolved 15 Oct** | The production-scheme sensitivities (Table VIII) and the Damkoehler sweep (Fig. 3) are in the abstracts and Secs. 1, 4, 6.3, 7.3, 8 and 9; every `% TODO(review-A/B)` comment is removed. Title kept. Remaining: the researcher's sign-off of plan section 39, item 8. |

## 1. AIP / PoF requirements

| Requirement | Status | Notes |
|---|---|---|
| REVTeX 4-2 `aip,pof,reprint` | **done** | `pof/lesto-pof.tex`, `bash manuscript/pof/build.sh`. `--preprint` gives a single-column double-spaced copy. |
| Layout | **done** (note) | REVTeX's own `aip4-2.rtx` sets PoF to **single column at reprint stage** ("POF journal to be one column at Reprint stage"); the log says "pof journal style Single column". The two-column assumption does not apply to PoF. Shared sections still use `figure*`/`table*` and two-line equations, and a two-column AIP test build (`jap`) compiles with 0 overfull boxes, so a later switch is safe. |
| `\author`/`\affiliation` | **done** | Affiliation kept as in `main.tex`: "Paul Scherrer Institute, 5232 Villigen PSI, Switzerland". **decision**: PSI's current preferred affiliation string (e.g. centre/laboratory name). |
| Corresponding author e-mail | **placeholder** | `\email{...}` reads verbatim, so the placeholder is plain bracketed text, not red. |
| ORCIDs | **placeholder** | Entered in the submission system, not in the TeX. Not invented here. |
| Title | **done** (changed 14 Oct); **blocked** on review-B | Alternative 2 of Section 6 adopted; see Section 6 for the basis. arXiv title unchanged. A `% TODO(review-B)` above both titles: do not change them before the Damkoehler sweep is in. |
| Abstract ≤ 250 words, one paragraph, no citations/equations/footnotes | **done** | 250 words after the last package runs (15 Oct; inline math counted as one word; `%` comment lines excluded; note that `\%` is text, not a comment). Adds the Damkoehler sweep and the production-scheme sensitivities; the silica-coordinate sentence was moved out to make room. Recount after any edit. |
| Keywords | **decision** | AIP's instructions state no keyword requirement; the TeX has none. If the submission form asks: thermal-gradient tube flow; vapor deposition; Hertz–Knudsen–Schrage; laminar mass transfer; finite-volume method; OpenFOAM; numerical verification. |
| Acknowledgments | **placeholder** | `acknowledgments` environment. The earlier draft's funding statements and grant numbers were removed pending confirmation; none invented. |
| Author Declarations: Conflict of Interest | **placeholder** | AIP standard text "The authors have no conflicts to disclose." inside the placeholder; confirm or replace with AIP's disclosure template. |
| Author Declarations: Author Contributions (CRediT) | **placeholder** | Roles must be chosen from the 14 CRediT roles (Conceptualization, Data curation, Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing – original draft, Writing – review & editing). Not assigned here. |
| Ethics approval | **done** (n/a) | No human or animal research. |
| Data Availability | **placeholder / blocked (B2)** | Now also names the production, scheme/mesh and iteration results and the Lobresco et al. digitised data (CC BY 4.0) and comparison runs. Based on the existing `main.tex` availability text, worded in AIP's "openly available in [repository]" form, plus the third-party/licence facts: digitised curves from Liu et al. under CC BY 4.0, external thermodynamic coefficient files not redistributed. The repository URL `https://github.com/JanKren/LESTO` answers HTTP 200 anonymously (public), but the relevant branch is not pushed. AIP's third-party template ("Restrictions apply to the availability of these data, which were used under license...") may be preferred for the coefficient files once their licence is known (requests-draft item 2). |
| Reference style | **done** | `aipnum4-2` (REVTeX default for `aip`), numbered by citation order, superscript in reprint. BibTeX's only warning, `jnrlst (dependency: not reversed)`, is REVTeX's internal control entry. |
| Cover letter | **done** (note) | AIP: "Cover letters are not accepted with submissions for Physics of Fluids." `cover-letter.md` is drafted as requested but must **not** be uploaded to PoF. |
| Alt text for every figure | **placeholder** | AIP: "Providing alt text is a requirement for publication." Eight figures need alt text, entered at submission. Not written. |
| Figures: format, fonts | **done** | Vector PDF, one embedded subset font (DejaVuSans, CID TrueType) per figure; no unembedded fonts. |
| Figures: size and minimum 8 pt | **done** (14 Oct) | Section 5. |
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

- 14 October 2026: sections 01, 04–08 updated for the campaigns (new
  Sec. 6.3 on Lobresco et al.); figures included at natural size;
  `tab:peaks` is now `table*`; `main.tex` abstract, date and data
  statement updated; `references.bib` gains van Leer (1974).

## 3. Simulation campaigns (B1): integrated 14 October 2026

Evidence and numbers: plan section 39 (pending the researcher's sign-off)
and `manuscript/review-notes.md`.

1. **Production limitedLinear runs on 41,472 cells**: done. 19 cases
   (`run/024-liu2025/production_576x72_limitedLinear`). Sec. 4.2 states the
   scheme with Jasak 1996 and Sweby 1984; the follow-up sensitivities are
   stated as 18,432-cell upwind differences from their own baselines.
   Tables III, IV, VI, VII and Figs. 3, 4, 6, 7 regenerated.
2. **Fourth mesh level and van Leer**: done
   (`run/024-liu2025/scheme_mesh_study`). Table III is now the scheme/mesh
   table; Sec. 5.2 quotes the 1.98 % L3–L4 change and no formal GCI
   (inconsistent apparent orders). van Leer (1974) added to the bib.
3. **8-versus-3 correctors over a full run**: done for three cases
   (`run/024-liu2025/iteration_check`); Secs. 4.2 and 5.1 updated. The
   fourth (shifted-BiI) case did not finish and is not used.
4. **Head-to-head on Lobresco et al.**: done as a comparison with their
   published measurements and published model results (new Sec. 6.3,
   Table V, Fig. 5; `run/023-lobresco2026`,
   `thermochemistry/benchmarks/lobresco2026`). No T-Flows inputs were
   needed or copied; the T-Flows wall table is only described.

After each further change: `bash manuscript/build.sh --figures`,
`bash manuscript/pof/build.sh`, `python3 manuscript/check_results.py`
(it now also recomputes the scheme/mesh, iteration and Lobresco tables
and checks the PoF log), and recount the abstract (≤ 250 words).

## 4. Dimensionless numbers (added, Sec. 4.1)

Computed by `manuscript/estimates/dimensionless.py` (output in
`dimensionless.txt`) from inputs stated in the paper only: d = 4.8/5.0 mm,
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

## 5. Figures (B4): resolved 14 October 2026

The PoF reprint text width is 396 pt = 5.48 in (measured with
`\typeout{\the\textwidth}` for `aip,pof,reprint`; preprint 468 pt,
arXiv 6.5 in). `prepare_results.py` now draws every figure at its printed
size, at most 5.4 in wide (exact `figsize`, no `bbox_inches='tight'`), and
the sections include them with `\includegraphics{...}` at natural size,
so the font sizes are the printed sizes in all three layouts:

- body text and axis labels 9 pt; tick labels, legends and the
  annotations of the schematic 8–9 pt;
- `Temperature (°C)`, `PbI$_2$`/`BiI$_3$` subscripts, `5° wedge`;
- `pdf.fonttype 42`, DejaVu Sans embedded (checked with `pdffonts`).

| Figure | Size (in) |
|---|---|
| temperature_profiles | 4.2 × 2.7 |
| experiment_schematic | 5.4 × 1.9 |
| refinement (new content: scheme/mesh study) | 5.4 × 4.9 |
| baseline_profiles | 5.4 × 6.6 |
| lobresco_comparison (new) | 5.4 × 6.9 |
| source_thermodynamics | 5.4 × 2.9 |
| silica_alignment | 5.4 × 5.3 |
| physical_sensitivity | 4.6 × 2.3 |

Remaining detail: mathtext subscripts (the "2" of PbI$_2$) and the
log-axis exponents render at about 70 % of the surrounding size
(about 6 pt). If the production editor counts them, raise `font.size`
to 11 for those labels.

## 6. Title

**Changed 14 October 2026** to alternative 2 below, because the completed
results support the transport-limited reading: the wall Damköhler number
is 1.4e3–2.9e3 at α = 1; α = 0.1 changes iodine profiles by at most
0.41 % and α = 0.01 by at most 4.21 %, against 3.83–14.20 % for halving
or doubling D; and the untuned HKS wall law with near-wall diffusion
resistance gives a smaller onset RMS (18 K) on the Lobresco cases than
both published models (35 K). The arXiv build keeps the descriptive title.
Previous PoF title: "Conservative CFD modelling of iodine deposition in
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
3. T-Flows / Lobresco et al. authors: the thermocouple wall tables per
   flow rate (the test-case table differs from the published curves), boat
   geometry and position, whether the figure origin is the inlet, and the
   C_e behind their Table 2 vapour-pressure column. The comparison is done
   with digitised curves; these would sharpen it.

## 9. References

- 23 cited entries (Vijaya Kumar 2021 and Di Ronco 2021 added after the internal review; Jasak 1996, Sweby 1984 and van Leer 1974 added for
  the schemes); all articles have journal, volume, pages or article
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
  a mixed-convection reference if buoyancy is discussed (below).

## 10. What a PoF reviewer is likely to criticise

1. **Spatial convergence.** Partly answered: production runs use
   limitedLinear on 41,472 cells; the most sensitive case changes by
   1.98 % from 41,472 to 93,312 cells, and van Leer agrees within
   2.8 %. Still open: no consistent observed order (no formal GCI), the
   fourth level exists for one case only, and the 64 sensitivities are
   upwind differences on 18,432 cells.
2. **Constant diffusivity D = 1e-4 m²/s for all nine gases and all
   temperatures.** Real gas diffusivities scale roughly as T^1.75/p over
   400–1173 K (a factor of about 6), and Sc would be 2–11 at 1 atm.
   Only one variant uses a T- and p-dependent PbI2–He value.
3. **Uniform cross-sectional temperature and frozen carrier.** No
   conjugate heat transfer, no thermal entrance, no radial temperature
   profile; the carrier sees the wall temperature everywhere. The estimated
   wall-to-bulk temperature lag of developed flow is <= 2.6 K (Sec. 4.1),
   which supports the approximation away from the furnace entrance.
4. **Buoyancy: estimated, not excluded (Sec. 4.2, `manuscript/estimates/`).**
   Gr <= 2.2, Gr/Re^2 0.15-0.78 in the 400-750 K deposition range and about
   1.2 in the steepest cold silica section; compositional Gr/Re^2 <= 0.32.
   Since the internal review the text no longer calls this weak: Richardson
   numbers of order 0.1-1 do not show buoyancy to be negligible, it may be
   significant locally in a horizontal column, and it is listed as an
   unquantified limitation (orientation unknown; ask the experimenters).
5. **Soret diffusion: bounded (Sec. 4.1).** With alpha_T = 1 as a generous
   bound, Soret drift <= 2.2 % of the carrier velocity and radial Soret flux
   <= 0.7 % of a Fickian wall flux.
6. **Axial diffusion and inlet species conditions.** Now stated (Sec.
   4.1): inlet Y = 0 fixed, outlet p = 101,325 Pa fixed (assumed) with
   zero-gradient species, wall zero gradient plus face fluxes; Sc 2.0-11.3
   and Pe 5.9-17.4 (100 mL/min) / 2.6-7.8 (45 mL/min) at that pressure; the
   diffusive loss through the inlet is in the ledger. Not quantified: its
   size.
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
10. **Length.** 27 pages in PoF reprint layout after review part 2 (25
    before; the Damkoehler figure takes a float page of its own), eight
    tables and eight figures; property appendix, case inventory, ledger
    diagnostics, bin-width audit, Lobresco set-up details, the full
    Graetz table and the sensitivity bar chart are in the 6-page
    supplementary material (`supplementary.tex`).

## 11. Build record (15 October 2026, after review part 2)

| Build | Command | Pages | Overfull | Undefined refs | BibTeX |
|---|---|---|---|---|---|
| PoF reprint | `bash manuscript/pof/build.sh` | 27 | 0 | 0 | 1 benign (REVTeX `jnrlst`) |
| PoF supplement | (built by the same command) `pof/lesto-pof-supplementary.pdf` | 6 | 0 | 0 | 0 |
| PoF preprint | `bash manuscript/pof/build.sh --preprint` | 42 | 0 | 0 | same |
| arXiv | `bash manuscript/build.sh --figures` | 23 | 0 | 0 | 0 |
| arXiv supplement | (same command) `lesto-iodine-supplementary.pdf` | 6 | 0 | 0 | 0 |
| Consistency | `python3 manuscript/check_results.py` | PASS | | | rewrote `review-checks.json` |

The checker now also recomputes the 20 Graetz-Robin GCI triplets from
`manuscript/data/graetz/fits.json`, checks them against the solver test's
`table.md`, the criteria log (27 PASS) and the solver hash, recomputes the
Lobresco from-inlet onsets and RMS, and checks the supplement logs and PDF.
Since review part 2 it also checks the Damkoehler sweep record against case
016's `expect.json` and recomputes the manuscript numbers of the sweep, and
recomputes the eight production-scheme sensitivities from their tracked
profiles.

## 12. Internal review revision (14 October 2026)

See `manuscript/review-notes.md`, top section. Supplementary material:
PoF accepts it as a separate file; upload `pof/lesto-pof-supplementary.pdf`
(arXiv-style layout; references are self-contained). References added:
Vijaya Kumar et al. 2021 (doi:10.1016/j.applthermaleng.2021.116546) and
Di Ronco et al. 2021 (doi:10.3389/fenrg.2021.692627, `pages` added to the
entry); 23 entries cited in the main text. Open: `% TODO(review-A)` and
`% TODO(review-B)` comments (running simulation sets) must be resolved
before submission; `grep -rn "TODO(review" manuscript/` lists them.
