# Cover letter for IJHMT (draft, unsigned)

> Every number below is taken from the shared journal abstract
> (`../sections/abstract-journal.tex`); update both together.  Statements
> marked **[AUTHOR TO CONFIRM]** must be checked before use.

---

Dear Editor,

We submit the manuscript "Competition between radial diffusion and
interfacial exchange in thermal-gradient deposition columns: conservative
modelling of iodine transport" by Jan Kren and Bojan Ničeno (Paul Scherrer
Institute) for consideration as a research paper in the *International
Journal of Heat and Mass Transfer*.

**The mass-transfer problem.** In a thermal-gradient column a slow laminar
helium flow (Re ≈ 0.7–3) carries trace vapours along a tube whose wall
temperature falls by several hundred kelvin. Where they deposit is set by
the competition between axial advection, radial diffusion to the wall and
a temperature-dependent interfacial exchange law: an extended Graetz
problem with a nonlinear, reversible wall condition. The application is
fission-product (iodine) transport in lead-cooled reactor safety, but the
questions are those of laminar internal mass transfer.

**What the paper contributes.**

- A conservative finite-volume treatment of wall exchange: a
  Hertz–Knudsen–Schrage law in series with the near-wall diffusive
  resistance, coupled to element-conserving gas equilibrium of nine
  Pb–Bi–I gases and shared wall inventories. On a Graetz–Robin problem with
  a reference solution, this face flux converges at second order, whereas
  the same law applied as a cell-centre sink converges at first order.
- A wall Damköhler sweep that reproduces the kinetic and diffusion limits
  and the series-resistance transition, and shows that radial diffusion,
  not interfacial kinetics, limits wall exchange in the studied cases.
- A separation of numerical error (time step, correctors, advection
  scheme, four mesh levels) from input sensitivity (diffusivity,
  accommodation, thermodynamic data, source and coordinate assumptions),
  with element balances closed to 1.41 × 10⁻¹² of the supplied material.
- An untuned comparison with five published deposition experiments: an RMS
  onset-temperature deviation of 18 K, against 35 K for two published
  models, reported with the reasons it does not yet isolate the wall
  treatment.

**Scope.** The work is a verified method and sensitivity study. It does
not claim independent predictive validation, reports the remaining mesh
sensitivity of local profiles, and states which experimental inputs are
missing. Solver, inputs, results and provenance records are openly
available **[AUTHOR TO CONFIRM: archived release / DOI]**.

This manuscript has not been published and is not under consideration for
publication elsewhere **[AUTHOR TO CONFIRM; state the arXiv preprint if
posted]**. All authors have approved the submitted version **[AUTHOR TO
CONFIRM]**. The authors declare no competing interests **[AUTHOR TO
CONFIRM]**.

Suggested reviewers (optional): **[AUTHOR TO CONFIRM: names, affiliations,
e-mails; none are proposed here]**

Sincerely,

[unsigned draft]
Corresponding author: Jan Kren, Paul Scherrer Institute, 5232 Villigen PSI,
Switzerland, **[AUTHOR TO CONFIRM: e-mail]**
