# Cover letter (draft, unsigned)

> **Read first.** AIP's author instructions state: "Cover letters are not
> accepted with submissions for Physics of Fluids." (publishing.aip.org,
> author instructions, checked 13 October 2026.) Do not upload this as a
> cover letter to PoF. It is kept as (a) a source for any short "comments
> to the editor" field the submission system offers, and (b) a ready
> letter if another journal is chosen. All statements marked
> AUTHOR TO CONFIRM must be checked before any use.

---

Dear Editor,

We submit the manuscript "Conservative CFD modelling of iodine deposition
in thermal-gradient columns: numerical verification and input sensitivity"
by Jan Kren and Bojan Ničeno (Paul Scherrer Institute) for consideration
as a regular article in *Physics of Fluids*.

**What the paper does.** Thermal-gradient (thermochromatography)
columns separate volatile species by carrying them in a slow laminar
helium flow along a tube whose wall temperature falls by several hundred
kelvin. Where they deposit is set by the competition between axial
advection, radial diffusion to the wall and a temperature-dependent
interfacial exchange law. We present a finite-volume model in OpenFOAM
that transports nine Pb–Bi–I gases in a frozen laminar carrier, couples
them through element-conserving gas equilibrium, and exchanges material
with shared condensed wall inventories through a Hertz–Knudsen–Schrage
law in series with the near-wall diffusive resistance. In 79 completed
computational cases built on three published experiments, we separate
numerical resolution (time step, three mesh levels, correctors) from
input sensitivities (diffusivity, accommodation, flow reference, source
loading, thermodynamic data and experimental coordinate registration).

**Why it fits Physics of Fluids.** The paper treats a reacting,
low-Reynolds-number (Re ≈ 0.7–3) internal flow with a nonlinear wall
mass-transfer condition, a Graetz-type problem in which the wall
concentration is set by the local wall temperature and deposit. Its
emphasis is on conservative discretisation of transport coupled to wall
exchange, and on an explicit separation of discretisation error from
input uncertainty. Integral deposition proves robust under refinement,
while local profile shapes are not yet spatially converged; the paper
reports this openly rather than claiming validation.

**Scope and limitations.** The work is a methods and input-sensitivity
study. It does not claim mesh independence or independent predictive
validation, and it states which experimental inputs are missing.

This manuscript has not been published and is not under consideration
for publication elsewhere. **[AUTHOR TO CONFIRM]** All authors have
approved the submitted version. **[AUTHOR TO CONFIRM]** The authors
declare no competing interests. **[AUTHOR TO CONFIRM]**

Suggested reviewers (optional): **[AUTHOR TO CONFIRM: names, affiliations,
e-mails; no names are proposed here]**

Sincerely,

[unsigned draft]
Corresponding author: Jan Kren, Paul Scherrer Institute, 5232 Villigen PSI,
Switzerland, **[AUTHOR TO CONFIRM: e-mail]**
