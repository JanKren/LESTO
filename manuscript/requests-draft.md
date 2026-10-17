# Draft requests for external inputs (not sent)

Status: drafts only, 2026-10-09 (section 4 added 2026-10-15). Nothing has been sent. Greeting names are
placeholders. Sources: `review-notes.md`, `notes-for-authors.md`,
`sections/02-experiments.tex`, `06-validation.tex`, `07-multispecies.tex`,
plan sections 26, 27, 36.1, 36.7, `thermochemistry/README.md`.

## Which manuscript statements each answer would unlock

| Request | Statement now unsupported | What an answer allows |
|---|---|---|
| 1. Liu et al. group: original T(x), scan coordinates | Sec. 2 "origin relative to the physical boat and furnace has not been independently confirmed"; Sec. 7 silica shifts of 4.096/4.653 cm are "input-derived hypotheses" (Q1 overlap 7.20 % to 75.29/81.00 %; low-loading Q2 79.16/83.33 %) | Confirm or reject a shift; replace digitised curves (about 4 K, 1 mm) by measured ones; run the missing BiI3_SiO2_800 case |
| 1. Controller reference conditions, flow reference | "supplied paper does not resolve ... flow-reference conditions" | Fix the carrier-gas state and temperature coordinate used for the 45 and 100 mL/min cases |
| 1. Source release history, boat geometry and position | Sec. 2/3: 60 s runs keep the three-hour mean rate (1/180 of inventory); boat obstruction omitted; source interval assumed | Time-resolved source instead of a constant mean; justify or change the source position; boat obstruction model |
| 2. HERACLES maintainers: Bi-I provenance and licence | Review notes: "confirm thermodynamic record provenance and redistribution terms"; Bi-I data cannot enter the repository (plan 36.1.2) | Publish a reproducible CFD input package; cite the data correctly |
| 2. BiI(g) and BiI3(l) records | The +29 kJ/mol BiI(g) shift is a "sensitivity surrogate" (BiI3 iodine share on refined steel: 94.57 % to 2.00 %, i.e. PbI2 share about 5 % to 98 %); BiI3(l) is "constructed" | Replace the surrogate branch by a qualified datum, or keep both and state which one holds |
| 2. Pb(g) heat-capacity defect | Plan 36.1: HERACLES Pb(g) is 3.9 kJ/mol off at 1000 K; NASA record used instead | State that the correction is justified, or that the defect is fixed in a newer release |
| 3. T-Flows authors (D, E, F) | The paper makes no T-Flows comparison and claims none. These answers matter only if a comparison is added later (see relevance marks) | Optional benchmark against Lobresco et al. (2026) |

Not requested (outside this task, still open in `notes-for-authors.md`):
items A, B, C (non-expanding carrier velocity, HKS unit conventions and
multipliers, interaction-layer rule).

---

## 1. To the Liu et al. (2025) experimental group, PSI

**Subject:** Request for original temperature profiles and run conditions, Liu et al. (2025) PbI2 and LBE-I columns

Dear [Name],

I am reconstructing your PbI2_SS, LBE-I_SS_II and LBE-I_SiO2_I experiments
for a CFD study of fission-product transport in a Pb-Bi-I system. So far I
have only your published figures, digitised to about 1 mm and 4 K. Could you
share the following, in any machine-readable form (CSV or Excel)?

1. **Axial T(x) profiles** (thermocouple readings, 1 cm or finer) and the
   **scan coordinate definition** (origin relative to boat and furnace) for
   each experiment, including **BiI3_SiO2_800**. Its profile was not in the
   paper, so I have not run that case. Why: in the silica LBE scan the
   labelled peak temperatures (306 and 127 C) occur 4.1-4.65 cm away from
   the labelled positions on the digitised curve. Shifting the curve by that
   amount moves the profile overlap from 7 % to 75-83 %. I cannot tell
   whether this is a registration issue in my input, and I do not want to
   present a fitted shift as a result.
2. **Controller reference conditions**: set points and where the reference
   temperature and gas flow were defined (standard conditions, location).
   Why: the carrier-gas state and the flow reference are currently
   unresolved.
3. **Source release history** (mass or mole rate against time, or any
   post-test mass of the boat/sample). Why: I use the three-hour mean rate
   inferred from the deposits, and cannot test the release shape.
4. **Boat dimensions and axial position** relative to the scan coordinate.
   Why: boat obstruction is omitted and the source position is assumed.
5. **Column orientation** (horizontal, vertical or inclined) and the
   absolute gas pressure in the column. Why: at the low flow rates the
   estimated Richardson number reaches about 1 where the wall temperature
   falls most steeply, so buoyancy-driven secondary flow matters only if the
   column is horizontal; the pressure fixes the Peclet number.

I will cite the data as you prefer and send you the manuscript draft before
anything is circulated.

Kind regards,
Jan Kren

---

## 2. To the HERACLES thermodynamic-database maintainers, PSI

**Subject:** Provenance and redistribution terms of the Bi-I records in HERACLES

Dear [Name],

I use the Bi-I records of HERACLES-TDB (v0.2, ThermoHub) for a CFD
speciation study of Pb-Bi-I vapour transport. Four questions:

1. **Provenance and licence.** The values for Bi(g), Bi2(g), BiI(g), BiI3(g),
   Bi(cr, l), BiI(cr) and BiI3(cr) agree with Barin (1995). Are they taken
   from Barin (1995) or an independent source, and may derived tables
   (vapour-pressure tables, GEMS3K input files) be redistributed with a
   journal paper and an open case package? The repository is GPL-3.0, but
   Barin's compilation is copyrighted, so until I know the terms I keep all
   derived files internal.
2. **BiI3(l).** Is there a liquid record? I found none and built one from
   the crystal at T_fus = 680.85 K, a fusion enthalpy of 31.798 kJ/mol and
   Cp(l) = 150.624 J/(mol K). It gives a normal boiling point of about
   800 K against 815 +- 1 K in the literature, which matters only above
   about 9 kPa.
3. **BiI(g).** Barin and TKV give 74.6 kJ/mol. Which value does the HSC
   MainDB hold, and what is its primary source? A shift of +29 kJ/mol moves
   the predicted PbI2 share of deposited iodide iodine from about 5 % to
   98 % on my steel case, so this single value decides the result.
4. **Pb(g).** The HERACLES record carries a heat capacity that resembles
   condensed Pb and is 3.9 kJ/mol off at 1000 K. I replaced it by the NASA
   record. Is this a known defect, and is a corrected record available?

Format: a short reply, or the reference list of the records, is enough.

Kind regards,
Jan Kren

---

## 3. To the T-Flows authors (Lobresco et al., Nucl. Eng. Des. 459:115201)

**Subject:** Questions on the temperature-based deposition model and cases in your 2026 NED paper

Dear [Name],

I have been comparing my OpenFOAM implementation with your NED 2026 paper
and the T-Flows branch `advanced_scalar_model` (LESTO_HKS, commit 35c1f6406e).
I cannot tell whether that branch produced the published figures, so these
are questions rather than claims.

**Relevance note:** the current OpenFOAM manuscript is framed as a method and
sensitivity study and makes no T-Flows comparison. Items 1 and 3 are only
needed if such a comparison is added. Item 2 and the wall tables are the most
useful, since they also give measured temperature input.

1. **Temperature-based model (still relevant for a future comparison).**
   Figs 9 and 10 share A = 220 1/s, k = 8e-3 1/K and T_dep = 680 K, but
   differ: Fig. 9 has a spike at x = 0.502 m, Fig. 10 rises to a peak at
   0.512 m. Which parameter set and mesh produced each? Which values are
   final (the thesis branch has A = 150, k = 0.02)?
2. **Numerics (relevant for any comparison).** Which advection scheme
   was used for the scalars? The GPU branch appears to apply a blending of
   central, linear-upwind and upwind whose default is pure upwind, not
   SuperBee.
3. **Case set-up (relevant for any comparison).** Per-case thermocouple wall
   tables (5, 106 and 540 mL/min and the Liu reference case; the repository
   has only `wall_x_profile_complete.dat`), source positions of Figs 13-14,
   whether the published runs' domain started at x = 0.04 m (as the
   temperature-model cases in the repository, 424 axial cells over
   0.04-0.69 m) or at x = 0, and the simulation time at which the Table 2
   curves were taken.

Format: text or CSV tables are sufficient. Without them I will say
explicitly that no comparison was done.

Kind regards,
Jan Kren

---

## 4. To Jörg Neuhausen, PSI (supersedes section 1)

Status: draft, 2026-10-15, not sent. Combines request 1 with a request for
further validation cases. Section 1 stays as the itemised record.

**Subject:** Iodine thermosublimatography (Liu et al. 2025): input data and possible further validation cases for the LESTO CFD model

Dear Jörg,

In LESTO Task 3.4 I have built a CFD model (OpenFOAM) of vapour transport
and wall deposition in thermosublimatography columns, coupled to GEMS3K for
the gas-phase speciation. I have reconstructed three of the experiments of
Liu et al. (2025): PbI2_SS, LBE-I_SS_II and LBE-I_SiO2_I. In the steel
columns the model reproduces the main deposition zones; the silica case
still disagrees (see A.1). A manuscript for the International
Journal of Heat and Mass Transfer is nearly ready. All inputs so far come from the published figures, so the
comparison is limited by what I could digitise. I would like to ask whether
you could help with two things.

**A. Inputs for the three experiments already modelled**

1. The measured axial temperature profiles (1 cm thermocouple readings) and
   the definition of the scan coordinate (origin relative to the boat and the
   furnaces). In the silica LBE scan, the labelled peak temperatures sit
   4.1-4.65 cm away from the labelled positions on my digitised curve. I
   cannot tell whether this is a registration issue on my side, and I do not
   want to present a fitted shift as a result.
2. Any information on the release history: a time-resolved signal, or the
   boat mass before and after the run. At present I use the mean rate over
   the three hours.
3. The reference conditions of the flow controller (standard temperature and
   pressure), the absolute pressure in the column, and the column
   orientation (horizontal or vertical). The last one decides whether
   buoyancy-driven secondary flow can matter in the steep cold sections.
4. The temperature profile of the BiI3_SiO2_800 experiment, which is not in
   the paper. With it I can add the fourth case.
5. Repeat measurements of any of these experiments, if they exist. They
   would give an experimental scatter to compare the model against.

**B. Further validation cases**

Would any of the following data be available, published or unpublished?

1. Runs in the same set-up at different carrier flow rates or evaporation
   times. The model predicts that deposition in these columns is limited by
   radial diffusion rather than by the wall reaction, so a flow-rate series
   would be a direct test of that prediction.
2. The raw temperature and deposition profiles of the Karlsson et al.
   experiments with iodine from LBE in fused silica (2020) and 316L steel
   (2021) columns, as a second, independent dataset.
3. Single-compound runs, for example pure BiI3 in steel and silica, with a
   known released amount. They test the wall model without the speciation
   uncertainty.

Machine-readable files (CSV or Excel) would be ideal, but scans of
lab-book tables would also help. I will of course cite the data as you
prefer. If unpublished data go into the paper, I would be glad to discuss
co-authorship, and I will send you the manuscript before it is submitted.

One related question: the Bi-I records I use from HERACLES agree with Barin
(1995). Do you know whether derived tables from them may be published with
the paper, or whom I should ask?

I would be happy to show you the current results in a short meeting if that
is easier.

Kind regards,
Jan
