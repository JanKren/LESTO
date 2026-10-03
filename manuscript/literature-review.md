# Literature review: transport and wall deposition of volatile fission products under strong thermal gradients

*Working document for the authors, not journal prose. Scope: the planned follow-up to Lobresco et al. [@Lobresco2026assessment], built on `rhoFixedFlowFoam`. Every citation is a key in `references.bib`, which holds the verified corpus: 297 entries, each checked against Crossref, doi.org, DataCite or a publisher/repository record. Numbers marked **(our estimate)** are back-of-the-envelope values computed for this review, not literature values. Their assumptions are listed in Section 10.*

---

## 0. Key messages

1. **Every published interpretation of thermochromatography (TC) and thermosublimatography (TSG) columns is 0-D or 1-D.** These are Zvara-type Monte Carlo models [@Zvara1985simulation; @Dietzel2024extended], analytical kinetic models [@Zhuikov2021kinetic], and the saturation formula of Liu et al. [@Liu2025chemical]. The only CFD of such a column is Lobresco 2026. It used a prescribed incompressible velocity profile, a molar-concentration formulation, and a single mesh. Liu et al. state explicitly that GEM must be coupled with CFD to describe transport along the gradient [@Liu2025chemical].
2. **Nuclear lumped codes fail in exactly this regime.** ASTEC/SOPHAEROS and VICTORIA miss wall condensation in developing laminar flow and entrance regions, and they miss heterogeneous chemistry [@Cousin2013modelling; @Kissane2006interpretation; @Gouello2021interaction; @Shepherd1995modeling]. A resolved CFD model of a well-characterized thermal-gradient tube has clear value beyond LFRs.
3. **In consistent SI units, PbI2 deposition in the column should be transport-limited.** The Hertz–Knudsen–Schrage (HKS) kinetic velocity is about 90 m/s, against a laminar mass-transfer coefficient of about 0.1 m/s (our estimate). Interface kinetics therefore matter only if σ is below about 10⁻³ [@Bouvet2024multiscale; @Aoki1998vapor]. Under this regime, C_e cannot legitimately change a converged deposition profile. A face-flux wall condition is the right, mesh-independent implementation [@VijayaKumar2021implementation].
4. **The onset temperature is an equilibrium quantity.** The Liu formula is the mass balance of a saturated, *expanding* carrier. With JANAF data it reproduces the measured T_dep of the long runs to within about 15 K (our estimate). The remaining discrepancies occur at 6 s and at 5 mL/min, and those are the cases where transport and non-dilute loading matter.
5. **Of the four "missing physics" candidates:**
   - Soret diffusion is negligible for the wall flux in this column, because the gas and wall temperatures differ by only 1–5 K.
   - Non-dilute loading matters near the source for T106_6 and T5_60, and during the release pulse of every run.
   - Homogeneous nucleation is plausible at 540 mL/min: the bulk saturation ratio is about 7 there, against about 1.2 at 106 mL/min.
   - The collision-integral fit error (about 9% in D) is a verification issue, not a leading uncertainty.
6. **GEMS3K is designed for per-cell coupling** [@Kulik2013gemselektor]. Nuclear CFD–GEM precedents exist, but only for liquid LBE [@Marino2020multiphysics] and molten salt [@Scuro2024coupledb; @Scuro2024coupleda]. No surrogate or tabulation study treats high-temperature gas–solid equilibria with condensed phases that appear and disappear.

---

## 1. Safety context: volatile fission and spallation products in lead and LBE systems

**State of the art.** The GIF LFR literature argues that lead is safe because it is chemically inert, operates at low pressure, boils at 1745 °C, and retains fission products [@Alemberti2014overview; @Cinotti2011lead; @Smith2016lead; @Smith2023lead; @Alemberti2014lead]. Lead is preferred to LBE mainly because it produces far less ²¹⁰Po [@Cinotti2011lead].

The main projects are:
- ALFRED and FALCON [@Tarantino2021overview];
- MYRRHA, an LBE-cooled ADS [@AitAbderrahim2012myrrha];
- SEALER, whose design claims no relocation beyond 1 km even after core melt [@Wallenius2018design];
- Russian LBE experience [@Gromov1997use; @Pankratov2004analysis].

LESTO continues PASCAL on ALFRED and MYRRHA licensing [@Gianfelici2025research].

The GIF safety documents rely on coolant retention of I, Cs and Po, and they state that it still needs to be verified by R&D [@Alemberti2020lead; @Alemberti2021safety]. The quoted volatilized fractions (about 10⁻⁶ for Cs-137 and 4×10⁻⁶ for I-131 at 700 °C) trace back to EU project deliverables, not to peer-reviewed models.

The source-side data come mostly from PSI and SCK CEN [@Neuhausen2020radionuclide; @Aerts2020behavior]:
- **Henry's-law constants** for Po [@Ohno2006equilibrium; @GonzalezPrieto2014equilibrium; @Aerts2021thermochemistry], Te [@Ivan2026evaporation], Hg/Tl [@Neuhausen2005study] and I [@Neuhausen2006investigations]. Iodine release from LBE becomes significant only above about 800 K and is controlled by desorption and evaporation.
- **Speciation.** Po is carried as PbPo and atomic Po, and as H₂Po when steam is present [@Mertens2019po; @Buongiorno2003speciation]. Te is carried as PbTe(g) [@Zobnin2026speciation]. Iodine is carried as BiI at trace level [@Karlsson2020thermochromatographic; @Karlsson2021thermochromatographic]. At x_I = 10⁻⁵–10⁻³ it deposits as PbI₂ and BiI₃, and BiI disproportionates [@Liu2025chemical]. An independent platform found PbI₂(g) dominant in dry Ar [@Liu2025exploring].
- **Integral in-facility data.** Most ¹²⁹I in the MEGAPIE target ended up at the LBE/steel interfaces [@HammerRotzler2015radiochemical].

System-level source-term tools treat partitioning as lumped equilibrium:
- equilibrium cover-gas models for MEGAPIE [@Neuhausen2011vapour] and ALFRED (more than 99.995% retained in normal operation; I, Cs and noble gases significant after fuel failure) [@Ivan2025assessment];
- SAS4A-FATE [@Hua2020development];
- TRAIL, a flow network [@Miyahara2023development], with Thermo-Calc speciation [@Miyahara2019analytical];
- Henry's-law constants plus aerosol-size distributions [@Huang2026source];
- a data-gap review toward a mechanistic LFR source term [@Trivedi2017retention].

Pure-lead retention data are sparse. Iodine volatilized from molten lead and formed PbI₂ [@Fitzhugh2024inventory]. Lead evaporation is now being measured for the mechanistic source term [@BaniMelhem2026experimental; @BaniMelhem2024integral]. BREST work addresses iodine compounds in lead [@Dubenkov2021potential] and lead-vapor transport, which is negligible at operating conditions [@Beznosov2001mass]. Cs/I/Pb/Bi condensed-phase thermodynamics is being built with CALPHAD [@VanHattem2023ternary; @VanHattem2026chemistry], and Cs release potential has been quantified [@Kou2027quantifying].

A second line of work concerns **fog**. Metal fog in the cover-gas boundary layer can raise fission-product (FP) vaporization by up to an order of magnitude [@Epstein2020enhancement]. Cs/I transport from lead to cooled steel was 24–109 times the pure-Pb baseline, and PbI₂ fluxes agreed with the fog model to within a factor of 2.2 [@BaniMelhem2026experimentalb] (preprint). Nanometer-scale PbTe-doped LBE aerosols [@Wang2025deposition] and LBE aerosol in Ar cover gas [@Li2026experimental] have been measured.

**Established:** coolants retain FPs strongly; speciation depends on concentration, atmosphere and wall material; Henry's-law correlations are consistent at high temperature, but below about 350–600 °C surface enrichment and redox effects distort them [@Ivan2026evaporation; @Huang2026source].
**Contested or open:**
- what fraction of the evaporated inventory reaches the gas space as vapor and what fraction as fog;
- Kou 2027 and the PNE lead-evaporation paper were verified only from metadata, so their content is unknown;
- whether the retention fractions quoted by GIF survive a mechanistic analysis.

**Implication for us.** The literature supplies the *source* boundary condition (Henry constants, speciation) and *lumped* models. It supplies no resolved link from source to wall, which is the LESTO Task 3.4 gap. The pure-lead versus LBE distinction matters: Bi iodides are absent in ALFRED-type systems, so PbI₂ is the right first species for lead, and BiI₃/BiI are the LBE extension.

---

## 2. Thermochromatography and thermosublimatography: technique and interpretation

**State of the art.** Temperature-gradient tubes go back to Merinis [@Merinis1969etude] and Eichler & Domanov [@Eichler1975verfluchtigung]. The field has two regimes.

*Trace (sub-monolayer) regime.* T_dep is set by the adsorption enthalpy. Eichler & Zvára derived second-law and quasi-third-law routes from the time and flow dependence of T_dep [@Eichler1982evaluation]. The forward model is Zvára's Monte Carlo [@Zvara1985simulation], extended to decay chains and chemisorption [@Dietzel2024extended]. Diffusion approximations exist for vacuum TC [@Zvara2010fundamentals; @Zvara2014vacuum], as does an analytical kinetic peak-shape model [@Zhuikov2021kinetic], and reviews cover the method [@Zvara2008inorganic; @Gaggeler2014gasphase]. Empirical correlations between ΔH_ads and ΔH_subl close the inverse problem [@Eichler1975verfluchtigung; @Eichler2005empirical; @Eichler2014thermochemical; @Zhuikov2019theoretical]. Surface state matters: examples are Po on hydroxylated quartz [@Hermainski2025reactivity], Po in various gases [@Maugeri2014thermochromatography], and Po/Bi on Au and 316L [@Maugeri2016adsorption; @Maugeri2018adsorption].

*Weighable-amount regime (TSG).* This is the regime of our PbI₂ runs. Eichler et al. introduced TSG [@Eichler1992volatilization] and placed the adsorption-to-desublimation transition at about one monolayer [@Eichler1993complex]. Liu et al. set the deposition onset at the temperature where p_v(T_dep) equals m R T₀/(v₀ t M), that is, the evaporated mass diluted in the total carrier volume at standard conditions. They report T_dep to be independent of the column material (silica versus 316L) [@Liu2025chemical]. Liu C. et al. found almost the same peak temperatures on silica and 316L in Ar [@Liu2025exploring]. Serov et al. used macro and carrier-free experiments to cross-check speciation [@Serov2011gas].

Analogues from outside the field:
- 3-D finite-element models of thermal-gradient GC columns (in which the gradient barely mattered) [@Singh2022method];
- stochastic models of GC columns with a thermal gradient [@Avila2021comparison];
- Monte Carlo TC for forensics, noting the lack of validation data [@Garrison2011monte].

**How the works relate.** The Monte Carlo line and the Liu formula answer different questions: kinetics of trace adsorption versus mass balance of macro desublimation. Lobresco 2026 took T_dep for its Eq. 5 from the Liu formula, so its temperature-based model inherits the formula's 0-D assumptions.

**Contested.** The influence of flow and duration on T_dep is disputed. Helas et al. found no effect of flow rate or amount [@Helas1978investigation]. Eichler & Zvára predict a logarithmic dependence [@Eichler1982evaluation]. The PSI data show T_dep falling with duration and with flow rate [@Lobresco2026assessment].

**Implication (our estimate).** The Liu formula is exactly the equilibrium front of a saturated carrier *expanding* with temperature, with T₀ = 273 K. Anchored to JANAF p_sat(700 K) = 27.5 Pa and ΔH_subl = 172 kJ/mol [@Chase1998nist], it gives 798/731/682/808/688 K for T106_6/60/300, T5_60 and T540_60. The measured values are 814/723/667/847/700 K [@Lobresco2026assessment]. The long runs agree within about 15 K. The 6 s and 5 mL/min runs deviate by 16–39 K, and a transport model should explain those two. The formula ignores radial resistance, finite release time, axial back-diffusion (Péclet number about 0.4 at 5 mL/min) and non-dilute loading; the CFD should quantify each of these.

---

## 3. Nuclear analogues: fission-product vapor deposition in thermal-gradient systems

**LWR primary circuit.** Phébus FPT0–FPT3 remain the reference [@Clement2003lwr; @Haste2013transport; @Girault2006iodine; @Girault2013insights]. Deposition concentrates where wall and gas temperatures drop steeply, and developing flows and diameter reductions enhance it. The processes that act together are vapor chemistry, wall condensation, homogeneous nucleation, condensation on aerosols, thermophoresis and revaporization [@Haste2013transport; @Girault2013insights].

The codes treat the circuit as a chain of control volumes [@Cantrel2014astec; @Cousin2008new; @Bowsher1987fission]. They reproduce the overall I, Cs and Mo retention, but they:
- disagree with the data where laminar flow is not developed [@Cousin2013modelling];
- underestimate condensation of supersaturated vapor on structures when entrance effects are ignored [@Kissane2006interpretation];
- needed CFD with particle tracking to explain under-predicted deposition in the steam generator [@Jones2003validation].

Falcon thermal-gradient tubes showed that condensation temperatures shift with speciation and that ideal-solution assumptions under-predict volatility [@Shepherd1995modeling].

**Is local equilibrium valid?** Several studies say no. Evidence of kinetic limitation at about 1000 K and below, in steam/H₂ atmospheres:
- the early kinetic transport model [@Cantrel2003reaction];
- I–O–H kinetics that keep gaseous iodine present at 428 K [@Xerri2012ab];
- HI predicted by equilibrium but I₂ measured [@Gouello2013analysis];
- the explicit statement that speciation "cannot be calculated assuming equilibrium" [@Gregoire2017study];
- JAEA kinetic datasets for Cs–I–B–Mo–O–H [@Miyahara2019chemical; @Miwa2020development].

**Thermal-gradient-tube (TGT) experiments** are the closest precedent for our column:
- VTT separate-effect tests on deposit chemistry at 400–650 °C [@Kalilainen2014chemical; @Gouello2018scoping], and a flowing TGT at about 1030→450 K in which the reaction rate fell with increasing flow rate and SOPHAEROS lacked the heterogeneous chemistry [@Gouello2021interaction];
- JAEA TeRRa [@Miyahara2020experimental; @Miwa2020boron; @Rizaal2021revaporization];
- CsI condensation and thermophoretic deposition in gradient pipes [@Maruyama1999vapor].

**Reactive walls.** CsOH chemisorbs on stainless steel through Si impurities [@DiLemma2016surface], and Mo in 316 steel adds further phases [@DiLemma2017experimental]. The rate constant depends on T, H₂/H₂O, concentration and Si content [@Nishioka2019experimental]. Nakajima et al. combined penetration-theory mass transfer with surface reaction [@Nakajima2020study]. CsI deposits react with Cr₂O₃ and release iodine [@Miyahara2020experimental; @Rizaal2021revaporization]. Metal iodides on Fe/Cr oxides can release I₂ [@Hu2021dft], and HI attacks the oxide on 304 steel [@Bowsher1985high].

**Liquid-metal cover gas.**
- *Sodium.* Deposition in shield-plug annuli came mainly from mist rather than vapor [@Himeno1979sodium]. One nucleation-related parameter fitted a series of tests [@Ford1993sodium]. Cs enrichment follows the vapor-pressure ratio, with extra enrichment at cold spots [@Minges1994experiments]. In He, a heavy vapor makes the gas density nearly temperature-independent and diffusion becomes more important [@Hotchkiss1977sodium].
- *Modern CFD.* Huang & He coupled condensation, aerosol and radiation [@Huang2019numerical]. Ohira et al. added RANS mist equations to OpenFOAM [@Ohira2022numerical]. Patel et al. simulated aerosol transport [@Patel2023mechanistic; @Patel2025mechanistic]. Anoop & Mangarjuna Rao built a nucleation-aware Eulerian mixture model [@Anoop2026development]. All of these are sodium in Ar, mostly in turbulent natural convection. Po adhesion filters for LBE are a related topic [@Obara2008polonium].

**Established:** thermal gradients and developing flow control where material deposits, and lumped codes are weakest there. **Contested:** whether local equilibrium can be assumed, which depends on system, carrier and residence time.
**Implication:**
1. Our resolved model is the missing tool for TGT-type validation.
2. Tabulated single-species p_v(T) is defensible, but multi-species GEMS coupling needs a quench-temperature or Damköhler check.
3. The SS316L column may not be inert for PbI₂. Liu reports T_dep independent of wall material in the desublimation regime, but the reactive-wall literature says this should be checked for trace amounts and for long durations.

---

## 4. Deposition physics: interface kinetics, near-wall transport and nucleation

**Hertz–Knudsen–Schrage.** The flux law goes back to Hertz, Knudsen and Schrage [@Hertz1882verdunstung; @Knudsen1915maximale; @Schrage1953theoretical]. Schrage's monograph already treats gas–solid interfaces and states when interfacial equilibrium can be assumed. Measured coefficients scatter over three orders of magnitude [@Persad2016expressions] and depend on the surface [@Marek2001analysis]. Schrage overpredicts fluxes by about 15% at σ ≈ 1 [@Vaartstra2022revisiting]. Kinetic boundary conditions are reviewed by Frezzotti [@Frezzotti2011boundary]; vacuum limits are discussed by Safarian & Engh [@Safarian2012vacuum].

With a non-condensable gas present, net phase change is suppressed and transport through the gas controls [@Aoki1998vapor], and the apparent coefficients fall [@Ohashi2020evaporation].

For halide solids:
- σ ≈ 0.3 for CsI [@Rothberg1959free];
- dislocations raise α of NaCl about twofold [@Lester1968studies];
- BCF step theory gives the mechanism [@Burton1951growth].

These are a physical basis for a C_e-like factor, but only if the process is kinetics-limited. Bouvet et al. find the kinetic-to-diffusive crossover near α ≈ 10⁻⁴ for snow [@Bouvet2024multiscale].

**Near-wall boundary condition.**
- In CVD, the Motz–Wise-type correction to kinetic wall fluxes errs by up to 45% at high reactant mass fractions, and a fitted correction exists [@Dorsmann2007general].
- Species wall conditions for precipitating mixtures have been formulated [@Johnsen2017wall], as has a wall function built on Maxwell–Stefan diffusion with relaxation to equilibrium [@Johnsen2015wallfunction].
- Fick's law belongs on a mass basis when the mixture molar mass varies [@Liao2007generalized].
- Diffusion with a condensing surface drives flow [@Greenwell1981numerical].
- Face fluxes beat volumetric first-cell sinks in accuracy, cost and grid independence [@VijayaKumar2021implementation].
- Enforcing saturation in interface cells is a CFD precedent for an equilibrium wall limit [@Pan2016saturated].

**Boundary-layer deposition theory.**
- Rosner's frozen boundary layer combines Fick and Soret diffusion with variable properties [@Rosner1979chemically]. Kohl et al. validated the onset from the equilibrium dew point [@Kohl1979theoretical].
- Frozen and local-equilibrium limits bound the rate [@Gokoglu1988significance].
- Thermophoresis augments deposition of heavy vapors and clusters on cooled walls [@Gokoglu1986thermophoretically].
- Fog formation within the boundary layer, collected by thermophoresis, makes deposition *low* and dependent on wall temperature far below the dew point [@Castillo1989theory]. Pre-existing particles scavenge vapor [@Castillo1988nonequilibrium].
- Fog criteria [@Rosner1968fog; @Epstein1970enhancement]; nucleation in cooled laminar tubes [@Pesthy1983theory]; maximum saturation on the axis when Le > 1 [@Barrett2000aerosol].
- Heterogeneous and classical nucleation theory [@Fletcher1958size; @Kashchiev2000nucleation]. PbI₂ deposit morphology depends on the substrate [@KoffmanFrischknecht2018tuning].
- Combustion analogues: both direct vapor deposition and condensation within the boundary layer matter [@Pyykonen2003modelling]; KCl condensation initiates deposits [@Zhou2007dynamic].

**Thermophoretic particle deposition** (if fog forms):
- Talbot force law [@Talbot1980thermophoresis];
- tube theory [@Walker1979thermophoretic; @Housiadas2005thermophoretic], including fine-aerosol deposition at high temperature with large gas-to-wall temperature differences [@MunozBueno2005deposition] (abstract not verified);
- MCVD efficiency E ≈ 0.8(1 − T_e/T_r) [@Walker1980thermophoretic];
- flux proportional to heat flux [@Batchelor1985thermophoretic];
- failure of the heat-transfer analogy in nuclear codes [@Fernandes1996modeling];
- mesh sensitivity in FLUENT [@Gutti2009thermophoretic];
- sectional nucleation in compressible PISO [@Frederix2017application];
- textbook [@Friedlander2000smoke];
- continuity between thermal diffusion of large molecules and thermophoresis [@GarciaYbarra1989thermophoretic; @Mason1962motion].

### Table B. Deposition and condensation models

| Model | Formulation | Parameters | Reversible? | Where used |
|---|---|---|---|---|
| Monte Carlo adsorption–desorption | random displacements in idealized laminar flow plus adsorption residence times | ΔH_ads (fitted), ΔS_ads (modeled) | yes (adsorption equilibrium) | trace TC [@Zvara1985simulation; @Dietzel2024extended; @Karlsson2021polonium; @Liu2025exploring] |
| Kinetic peak-shape model | adsorption–desorption kinetic equation, analytic solution | ΔH_ads, duration, gradient, release history | yes | trace TC [@Zhuikov2021kinetic] |
| TSG saturation formula (0-D) | p_v(T_dep) = m R T₀/(v₀ t M) | ΔH_subl, ΔS_subl | onset only | [@Liu2025chemical]; supplies T_dep for Lobresco Eq. 5 |
| Temperature-based relaxation | S = −A[1 − e^{−k(T_dep−T)}]φ in a 0.1 mm layer | A, k, T_dep (calibrated) | no | [@Lobresco2026assessment] |
| HKS kinetic source | j = C_e (2σ/(2−σ))(p_v − p)/√(2πMRT), volumetric via A_I/V_I | σ, C_e, p_v(T) | yes | [@Lobresco2026assessment; @Persad2016expressions] |
| Equilibrium (Dirichlet) wall | c_w = c_sat(T_w); transport-limited flux | p_v(T) | yes | analytical limit [@Tandon2004extended; @Shah1978circular]; CFD analogue [@Pan2016saturated] |
| Mixed/Robin wall with sticking coefficient | kinetic flux × non-dilute correction(w_s) | γ, near-wall w | optional | CVD [@Dorsmann2007general; @Johnsen2017wall]; FP deposition in OpenFOAM [@DiRonco2021eulerian] |
| Face-flux condensation with NCG | gas-side diffusion-layer flux at the wall face | none (resolved) | condensation | containmentFOAM [@VijayaKumar2021implementation; @Liao2007generalized] |
| Sub-grid wall function | 1-D Maxwell–Stefan plus relaxation to equilibrium | relaxation time | yes | fouling [@Johnsen2015wallfunction] |
| Frozen or equilibrium boundary layer | Fick + Soret, variable properties, equilibrium dew point | α_T, D | onset from equilibrium | salt deposition [@Rosner1979chemically; @Kohl1979theoretical; @Gokoglu1988significance] |
| Fog within the boundary layer + thermophoresis | vapor in equilibrium with aerosol; critical supersaturation or nucleation rate | S_crit, thermophoretic coefficient | — | [@Castillo1989theory; @Epstein2020enhancement; @Ford1993sodium] |
| Lumped vapor condensation | Sherwood correlations per control volume + equilibrium speciation | Sh, database | revaporization | SOPHAEROS/VICTORIA [@Cantrel2014astec; @Shepherd1995modeling] |
| Reactive-wall chemisorption | first-order surface reaction k(T, H₂/H₂O, c, Si) + mass transfer | k | partly | CsOH on steel [@Nishioka2019experimental; @Nakajima2020study] |
| Sectional aerosol in CFD | nucleation, condensation (and evaporation) in size sections | nucleation model | yes | [@Frederix2017application; @Anoop2026development] |

**Implication.** Estimate the Damköhler number before tuning anything. If deposition is transport-limited, as in Section 10, then (i) HKS with σ = 1 and the Dirichlet limit give the same answer, (ii) C_e has no physical meaning, and (iii) the remaining freedom lies in p_v(T), the kinematics and gas-phase nucleation.

---

## 5. Transport properties

**Binary diffusivity.** No measured D exists for PbI₂, BiI₃ or CsI in any carrier. The standard route is first-order Chapman–Enskog theory with an LJ(12-6) potential [@Hirschfelder1954molecular; @Chapman1970mathematical]. Collision integrals come from Neufeld (better than 0.1%) [@Neufeld1972empirical] or Kim & Monroe (arbitrary precision) [@Kim2014high].

- **Lobresco 2026** took the PbI₂ parameters from Brandt's polarizability correlation [@Brandt1956calculation], which Brandt himself calls unreliable for dipolar, non-spherical molecules. He parameters come from Svehla [@Svehla1962estimated], combined with the Lorentz–Berthelot rule.
- **Better routes:**
  - polarizability-based "improved LJ" potentials [@Cambi1991generalized; @Laricchiuta2007classical; @Pirani2026collision], automated from ab initio polarizabilities [@BellasChatzigeorgis2022transport];
  - a direct unlike-pair diameter [@Sharipov2025estimating];
  - quantum-chemistry parameters for heavy clusters [@Sharipov2014theoretical];
  - alternative combining rules [@Kong1973combining] and orientation averaging [@Mo2022determination];
  - arbitrary-order D and α_T from any potential with Peng [@Zhai2023peng].
- LJ can fail badly for weakly bound pairs (up to 40%) [@Jasper2014first].
- **Empirical checks:** Fuller [@Fuller1966new] (±30% for inorganics [@Tang2014compilation]), and compilations [@Marrero1972gaseous; @Poling2001properties].
- **Chemically close benchmarks:** I₂–He [@Gardner1992binarya], SnBr₂/SnI₂–Ar with T exponents 1.69–1.75 [@Gardner1992binaryb], and ZnSe (Zn + ½Se₂) in He at 1140–1280 K [@Schonherr1996gaseous].

**Helium** is known to better than 1% [@Kestin1984equilibrium; @Bich2007ab; @Arp1998thermophysical; @Petersen1970properties]. Kestin also gives He–Xe α_T, a verification case for a Soret implementation. Mixture viscosity follows Wilke [@Wilke1950viscosity].

**Formulation.** Binary Fick diffusion is exact in mass fraction with the mass-averaged velocity. The molar-concentration form is not equivalent in a non-isothermal gas [@Curtiss1999multicomponent; @Bird2002transport]. The multi-species extension needs Maxwell–Stefan or corrected mixture-averaged fluxes [@Krishna1997maxwell; @Giovangigli1999multicomponent]. The carrier flow is a low-Mach, variable-density problem [@Majda1985derivation].

**Established:** the Chapman–Enskog framework and the He properties. **Open:** the PbI₂–He potential, uncertain by roughly ±10–30%, which is far larger than the collision-integral issue.
**Implication:**
- Use Neufeld in physical mode and keep the published fit only in the T-Flows code-to-code mode.
- Bracket D using Brandt, Cambi/Laricchiuta and Fuller, validated on I₂–He and SnI₂–Ar.
- Propagate ±15% in D. Section 10 argues that this barely moves T_dep.

---

## 6. PbI₂ thermochemistry and data provenance

- Konings et al. measured the PbI₂ vapor pressure by Knudsen effusion. The sublimation enthalpies evaluated across datasets span 165–183 kJ/mol, and their own value is 173.1 ± 1.6 kJ/mol [@Konings1996infrared].
- JANAF gives ΔsubH°₂₉₈ ≈ 172.2 kJ/mol and T_m = 683 K; p_sat is 0.28 Pa at 600 K, 27.5 Pa at 700 K and about 1 bar at 1100 K [@Chase1998nist].
- Saturated vapor contains Pb₂I₄ dimers [@Hilpert1985molecular], and PbI₂(g) dissociation becomes significant above 1000–1300 K [@Rybak2002equilibrium].
- Lobresco 2026 compared GEMS/HSC p_v only with Binnewies & Milke [@Binnewies2002thermochemical]. The independent assessed sources are IVTAN/Gurvich [@Gurvich1991thermodynamic; @Belov1999ivtanthermo; @Belov2018ivtanthermo] and Barin [@Barin1995thermochemical].
- ThermoFun makes the GEMS dataset explicit and swappable [@Miron2023thermofun].
- Coolant properties, including Pb/Bi vapor pressures, which are among the least certain properties [@Sobolev2007thermophysical], come from the handbook [@OECD2015handbook], the Pb/Bi/Bi₂ vapor model [@Morita2006thermophysical] and lbh15 [@Panico2023lbh15].
- For the Cs and Bi extension: CsI–PbI₂–BiI₃, in which PbI₂ and BiI₃ form a continuous solid solution [@VanHattem2023ternary; @VanHattem2026chemistry].

**Implication.** T_dep shifts by RT²/ΔH ≈ 17–31 K per e-fold of p_v between 600 and 800 K (our estimate). The spread between data sources is therefore the leading uncertainty in the predicted onset. It outweighs D. By our model-draft digitization of Lobresco's Fig. 6, GEMS/HSC lies 1.3–10 times above Binnewies–Milke, which is worth 15–25 K in T_dep.

Two further points:
- The measured T_dep range (667–847 K) straddles T_m. Above 683 K the stable condensate is liquid, so the metastable gas–solid branch used by Lobresco 2026 overestimates p_v there. Four of the five runs deposit above T_m.
- Dissociation matters only at the source.

---

## 7. Coupling CFD with Gibbs-energy minimization

**Architecture.** Operator splitting is the standard: transport the component or element totals, then call the equilibrium solver in every control volume whose state has changed [@Kulik2013gemselektor; @Steefel2015reactive; @Leal2017overview].
- Yeh & Tripathi recommend sequential iteration on total component concentrations [@Yeh1989critical].
- Splitting error is well characterized:
  - a continuous-flux boundary gives an inherent mass-balance error, requiring kΔt < 0.1 for less than 5% error, and alternating the split order helps [@Valocchi1992accuracy];
  - symmetric sequential iteration has no splitting error, and Strang SNI is second order [@Carrayrou2004operatorsplitting; @Strang1968construction];
  - under stiffness, the stiff operator should go last [@Sportisse2000analysis].
- Local-equilibrium validity criteria exist [@Valocchi1985validity].
- GEM handles phases appearing and disappearing without a posteriori stability tests, and partial equilibrium mixes equilibrium with kinetics [@Leal2017overview].

**How the works relate.** GEMS3K was built for this job and has been coupled to OpenGeoSys [@Kosakowski2014opengeosysgem], CSMP++ (kinetics as metastability constraints, density feedback) [@Yapparova2017reactive] and COMSOL [@Azad2016comsolgems]. Other precedents are Reaktoro–FEniCS [@Damiani2020framework] and OpenFOAM–PHREEQC [@Soulaine2021porousmedia4foam; @Pavuluri2022reactive]; multi-code benchmarks exist [@Poonoosamy2021benchmarking]. In nuclear engineering:
- Fluent + HSC for LBE coolant chemistry, the SCK CEN route [@Marino2020multiphysics];
- OpenFOAM + Thermochimica for molten salts [@Scuro2024coupledb; @Scuro2024coupleda; @Piro2013thermochemistry];
- lumped GEMS/Thermochimica couplings to MELCOR [@Nichenko2021modelling; @Nichenko2022msr; @Kalilainen2020evaporation; @Gelbard2023application].

Kalilainen et al. show that mixture non-ideality lowers release compared with pure-compound p_v, which matters for LBE. Related work comes from CVD, where equilibrium–transport coupling showed thermal diffusion changing deposit composition [@Rouch1996thermodynamic], and from boilers (CFD + equilibrium + particle dynamics) [@Leppanen2014numerical].

### Table A. CFD/transport–thermodynamics coupling approaches

| Work | Transport code | Chemistry solver / database | Coupling scheme | Application | Acceleration |
|---|---|---|---|---|---|
| [@Marino2020multiphysics] | ANSYS Fluent (FV) | HSC Chemistry GEM + HSC DB | sequential splitting, iterative or non-iterative; calls only in cells away from equilibrium | LBE coolant chemistry, corrosion | selective calls, in-memory coupler |
| [@Scuro2024coupledb] | OpenFOAM (Euler–Euler) | Thermochimica + JRC-MSD | two-way, local equilibrium, transport-limited kinetics | MSFR FP retention and release | CPU cost reported (not read) |
| [@Scuro2024coupleda] | OpenFOAM | Thermochimica + JRC-MSD | local equilibrium with phase change | fluoride volatility (105 vs 180 min measured) | none reported |
| [@Nichenko2021modelling; @Nichenko2022msr] | MELCOR (lumped) | GEMS + HERACLES | equilibrium per control volume | VERDON-1 release; MSR aerosol | — |
| [@Kalilainen2020evaporation] | MELCOR | GEMS | GEMS p_v into pool-evaporation model | MSR fuel salt | — |
| [@Gelbard2023application] | MELCOR | Thermochimica or user tables | p_v and solubility per volume | MSR source term | user-specified tables |
| [@Cousin2008new; @Cantrel2014astec] | ASTEC/SOPHAEROS (1-D) | internal equilibrium | per control volume | LWR primary circuit | — |
| [@Kosakowski2014opengeosysgem] | OpenGeoSys (FE) | GEMS3K | per-node calls after transport | repository materials | threads + domain decomposition |
| [@Yapparova2017reactive] | CSMP++ (FE–FV) | GEMS3K | mass-conservative splitting; kinetics as metastability constraints | dolomitization | — |
| [@Azad2016comsolgems] | COMSOL | GEMS | SNIA | cement | — |
| [@Damiani2020framework] | FEniCS | Reaktoro | SNIA | benchmarks | — |
| [@Soulaine2021porousmedia4foam; @Pavuluri2022reactive] | OpenFOAM | PHREEQC | sequential | porous media | — |
| [@Leppanen2014numerical] | CFD | equilibrium chemistry + particle dynamics | — | recovery-boiler fume | — |
| [@Rouch1996thermodynamic] | CVD transport | local equilibrium (gas + surface) | coupled | Si₁₋ₓGeₓ deposition | — |
| [@Leal2020accelerating] | reactive transport | Reaktoro | SNIA | geochemistry | on-demand ML, 1–2 orders |
| [@Prasianakis2020neural; @Laloy2022speeding] | pore/Darcy; HPx | NN from speciation; DNN/kNN | replace solver | geochemistry, cement | NN beats LUT; 3–33× |
| [@DeLucia2021dectree; @DeLucia2021poet] | POET | equation-based geochemistry | master/worker | reactive-transport benchmarks | XGBoost with mass-balance check; DHT cache |
| **This work (planned)** | rhoFixedFlowFoam (OpenFOAM, variable density) | GEMS3K (ThermoFun; currently HSC-derived MainDB) | tabulated p_v(T) → per-interfacial-cell calls with warm start → surrogate | PbI₂/He column → Pb–Bi–I | TODO: tabulation / ISAT / NN |

**Established:** splitting theory, the GEMS3K interface, and nuclear precedents for liquids.
**Open:** gas–solid deposition, compressible mass-fraction to element-total conversion, and the splitting error at a wall sink whose driving p_v changes by orders of magnitude over a few cells.
**Implication:**
- Couple GEMS only where phase change can occur (interfacial and sample cells).
- Put the stiff phase-change step last.
- Report Strang versus SNI mass-balance errors against exact gas-plus-solid conservation.

---

## 8. Acceleration: tabulation, surrogates, reduced-order models

Chemistry costs 10–10⁴ times the transport step [@Leal2020accelerating]. The options, from least to most aggressive:
1. **Warm-started GEM.** One or two iterations per call [@Leal2014efficient]; GEMS3K offers a smart-initial-approximation mode (documented in its source, not in the corpus paper).
2. **Tables and caches.**
   - look-up tables [@Huang2018new];
   - ISAT at about 1000× [@Pope1997computationally], also on constrained-equilibrium manifolds (~500×) [@Tang2002implementation] and for algebraic surface chemistry [@Mazumder2005adaptation];
   - ISAT in wall-reaction CFD: 5–15× overall [@Bracconi2017in];
   - DHT caching and scattered work packages [@DeLucia2021poet];
   - OpenFOAM chemistry load balancing [@Tekgul2021dlbfoam].
3. **ML surrogates.**
   - early demonstrations [@Jatnieks2016datadriven; @Guerillot2020geochemical];
   - PSI NN surrogates [@Prasianakis2020neural] and the 2025 benchmark (1–4 orders) [@Prasianakis2025geochemistry];
   - local emulators needed for 7 components [@Laloy2022speeding];
   - mass-balance-checked fallback [@DeLucia2021dectree];
   - KANs for solid solutions [@Boledi2026kolmogorovarnold];
   - NN best among GEM approximations for gas-phase atmospheres [@Himes2023toward];
   - training-set design [@Bracconi2020training].

   *Cautions:* multi-step rollout drift [@Silva2024rapid]; DNN emulators bias Bayesian calibration while Gaussian processes stay reliable [@Laloy2019emulation].
4. **ROMs of the fields.**
   - finite-volume POD-Galerkin from PoliMi/CIRTEN and SISSA [@Lorenzi2016podgalerkin; @Stabile2018finite; @Stabile2026ithacafv];
   - buoyant liquid-metal flow, about 10⁵× [@Star2021podgalerkin];
   - multiphysics MSR [@German2022genrom];
   - bias correction [@Riva2024multiphysics];
   - POD-NN [@Hesthaven2018nonintrusive], a CVD reactor ROM [@Gkinis2019building] and DeepONet [@Lu2021learning].

**Implication.** For PbI₂ alone, exact 1-D tabulation of p_v(T) costs nothing and is exact, so a surrogate adds nothing. The value of surrogates appears only for Pb–Bi–I, where the state space is (T, p, b_Pb, b_Bi, b_I, …) with appearing phases. The planned staged comparison (table → ISAT/on-demand learning → NN/KAN → direct warm-started GEMS3K) is novel. It must report element conservation at the wall and compare the surrogate error with the discretization error (GCI).

---

## 9. Numerical methods and verification

**State of the art.**
- *Codes.* OpenFOAM is the FOAM finite-volume C++ library [@Weller1998tensorial]. Its unstructured discretization, error sources and bounded convection go back to Jasak [@Jasak1996error], and source-term handling is documented in textbook form [@Moukalled2016finite]. T-Flows is a collocated, unstructured finite-volume code [@Niceno2005unstructured]; recent T-Flows work verified phase change against analytical solutions and showed that the mesh type (hexahedral versus polyhedral) biases near-interface gradients [@Kren2026sharp]. Both codes share the same basis, so a code-to-code comparison can isolate modeling differences, provided the discretization is aligned first. A tailored OpenFOAM solver with a structured verification and validation hierarchy already exists in nuclear safety [@Kelm2021tailored].
- *Convection.* Lobresco 2026 used Superbee, the most compressive admissible TVD limiter [@Sweby1984high]. Bounded schemes trade the resolution of sharp fronts against accuracy on smooth solutions [@Waterson2007design], so the limiter may sharpen the deposition onset.
- *Stiff sources.* The relaxation sink S = −A(φ − Φ_eq) and the HKS wall term are stiff (the HKS relaxation rate in the wall layer is about 10⁶ s⁻¹, our model draft). Implicit linearization of the sink keeps the diagonal positive [@Patankar1980numerical]. Modified Patankar schemes are unconditionally positive and conserve gas plus solid exactly [@Burchard2003high].
- *Variable density.* Pressure–velocity coupling uses PISO/PIMPLE [@Issa1986solution]. The boundedness of transported scalars is not automatic in segregated solvers [@Herrmann2006flux].
- *Verification methodology.* The method of manufactured solutions [@Roache2002code], with an OpenFOAM toolchain [@Ramoa2022semi]. The GCI and Richardson extrapolation [@Roache1994perspective; @Celik2008procedure], with a least-squares variant for scattered convergence [@Eca2014procedure]. The V&V hierarchy of Oberkampf & Roy [@Oberkampf2010verification].
- *Analytical benchmarks.*
  - Graetz tube, the Dirichlet or transport-limited limit [@Shah1978circular];
  - radially variable diffusivity [@Tandon2004extended];
  - Graetz problem with a finite-rate surface reaction [@Gupta2001heat];
  - transient pulses with first-order wall loss [@Sankarasubramanian1973unsteady] and Taylor dispersion [@Taylor1953dispersion];
  - analytical verification of Eulerian fission-product deposition with mixed boundary conditions in OpenFOAM [@DiRonco2021eulerian]. Its follow-up showed that the near-wall gradient, not the turbulence model, sets the mesh requirement [@DiRonco2022multiphysics].

**What Lobresco 2026 did:** one mesh (228,960 cells), Δt = 1 ms, relative tolerance 10⁻³, sequential under-relaxed scalar solves, and no grid-convergence or analytical checks.

**Established:** the methodology. **Open:** it has never been applied to a thermochromatography geometry with an imposed axial wall-temperature ramp and a saturation-type wall condition.

**Implication.** The verification hierarchy should be:
1. conservation to round-off (closed box; open tube with gas + deposit + outflow);
2. Dirichlet and Robin Graetz cases, and a transient pulse case;
3. MMS for the variable-density species equation with the nonlinear source and the Soret term;
4. code-to-code comparison with T-Flows in "paper mode";
5. GCI on T_dep, peak position and width, and deposited fraction;
6. sensitivity to limiter and time step;
7. a splitting-error study (Section 7).

---

## 10. Physics we may be missing: literature and order-of-magnitude estimates

*Basis (our estimate).*
- He at 1 atm and 700 K: μ = 3.6×10⁻⁵ Pa·s, ρ = 0.070 kg/m³, ν = 5.2×10⁻⁴ m²/s, α = 7.7×10⁻⁴ m²/s, with Pr = 2/3 from the Eucken relation.
- D(PbI₂–He) = 1.28×10⁻⁴ m²/s (Chapman–Enskog with Neufeld Ω), giving **Sc ≈ 4.0 and Le ≈ 6.0**.
- d = 4.8 mm; k_m = 3.66D/d ≈ 0.10 m/s.
- Physical bulk velocity at 700 K: 0.012, 0.25 and 1.27 m/s at 5, 106 and 540 mL/min STP. The corresponding Re are 0.1, 2.3 and 12, and the mass Péclet numbers Ud/D are 0.4, 9.4 and 48.
- Wall gradient |dT_w/dx| ≈ 1300 K/m in the deposition zone.
- d ln p_v/dx = (ΔH/RT²)|dT_w/dx| ≈ 55 m⁻¹, i.e. a saturation decay length of 18 mm.

### 10.1 Thermal diffusion (Soret) of PbI₂ in He

**Literature.** Heavy species undergo significant Ludwig–Soret transport. In the heavy limit α_T·D/ν → 0.47–0.54 [@Rosner2000heavy]. Factors for heavy dilute molecules are larger in He than in H₂ [@Holstein1988thermal]. Soret changed CVD growth rates by 7–20% [@Jenkinson1984thermal] and depleted WF₆ at the wafer [@Kleijn1991transport]. It restructures the near-wall layer [@GarciaYbarra1997mass], and it reportedly shifts dew points at cooled surfaces [@Rosner2007soret] (the "tens of kelvin" figure is from a search snippet and is unverified). The Chapman–Cowling approximate model is the most accurate simple closure [@Zirwes2025assessment; @Chapman1970mathematical], and flame work shows that it can be evaluated cheaply [@Ern1998thermal]. In these applications the gas–wall temperature differences are large (hundreds of kelvin in CVD reactors).

**Estimate.** α_T ≈ 0.5·Sc ≈ 2. The gas core leads the wall temperature by only 3UGR²/(8α) ≈ 0.04, 0.9 and 4.6 K at the three flows. The radial wall gradient is about 25, 500 and 2600 K/m, the same order as the axial gradient (not much larger, contrary to our initial expectation). The radial Soret drift D·α_T·|∇T|/T is therefore 0.2% (106 mL/min) to 1% (540 mL/min) of k_m. The implied dew-point shift is about (RT²/ΔH)·α_T·ΔT/T, below 0.5 K. Axially, the Soret flux is about 7% of the Fickian flux at the front.

**Verdict.** Soret diffusion is not quantitatively important for PbI₂ wall deposition in this column. Include the term (it is cheap), verify it on He–Xe [@Kestin1984equilibrium], and report it as a negligible correction. It becomes relevant only for aerosol particles, via thermophoresis, and near steep wall-temperature steps (furnace edges), where the gas–wall ΔT may reach tens of kelvin. **Check the measured wall profiles for such steps.**

### 10.2 Non-dilute loading

**Literature.** Molar-basis diffusion-layer models err when the molar mass varies [@Liao2007generalized]. Stefan flow and diffusion-driven convection appear in physical vapor transport [@Greenwell1981numerical]. The kinetic wall-flux correction depends on the near-wall mass fraction [@Dorsmann2007general]. A heavy vapor in He changes both the mixture density and the transport regime [@Hotchkiss1977sodium]. Deposition fluxes in Pb–Cs–I are strongly non-additive [@BaniMelhem2026experimentalb].

**Estimate.** The mean mole fraction at the deposition front equals p_v(T_dep)/p. It is about 1.0% for T106_6 and 1.4% for T5_60, which means mass fractions of 0.54 and 0.62 and mixture densities 2.2 and 2.6 times that of He. The other runs stay below 0.1% (w < 0.1). During the release pulse, our model draft estimates peak x ≈ 3% (106 mL/min) and 12% (5 mL/min). These give w ≈ 0.78 and 0.94, density ratios of 4 and 15, and Stefan corrections of 3% and 14%. Near the 1073 K source at 5 mL/min, Gr/Re² ≈ 200 for Δρ/ρ ≈ 1.6. The dense pulse may therefore stratify or drive secondary flow if the column is horizontal (orientation to confirm).

**Verdict.** Non-dilute effects are important for T106_6 and T5_60 and for every run's source region, and negligible downstream for T106_300 and T540_60. This supports the mass-fraction formulation with mixture density. It also coincides with the two runs where the 0-D formula fails (Section 2).

### 10.3 Nucleation, fog and surface roughness (what C_e may be hiding)

**Literature.** C_e was introduced to stand for roughness, deposit morphology and nucleation [@Lobresco2026assessment]. Defects roughly double the evaporation coefficients of halides [@Lester1968studies], but that helps only if the process is kinetics-limited. Mist, not vapor, explained sodium deposits [@Himeno1979sodium]. A single nucleation parameter fitted cover-gas tests [@Ford1993sodium]. Fog can either lower or raise deposition [@Castillo1989theory; @Epstein2020enhancement]. With Le > 1, the maximum supersaturation sits on the axis [@Barrett2000aerosol]. MCVD shows efficient thermophoretic collection of particles formed in the gas [@Walker1980thermophoretic].

**Estimate.** In fully developed flow the bulk lags the saturation curve. S_bulk ≈ 1/(1 − N), with N = (ΔH/RT²)|dT_w/dx|·Ud/(4k_m). This gives N ≈ 0.008, 0.17 and 0.86, so **S_bulk ≈ 1.01, 1.2 and ≈ 7** at 5, 106 and 540 mL/min, and higher on the axis. Nucleation is therefore implausible at 5 mL/min, marginal at 106 mL/min and likely at 540 mL/min. At 540 mL/min Lobresco's vapor-pressure model overestimated the broadening most and let material escape the column. Any particles formed would deposit thermophoretically, or leave the column, depending on size. The HKS/Damköhler ratio is about 900, so C_e cannot represent faster wall kinetics in a converged model. If C_e improves agreement, it is compensating for something else: kinematics, units, the wall layer, or nucleation (see `notes-for-authors.md` §A–C).

**Verdict.**
- Replace C_e by a nucleation *indicator*: map S and the classical-theory rate [@Kashchiev2000nucleation] as a post-processing check.
- Add a sectional or moment fog model [@Frederix2017application] only if S exceeds about 3–5 at 540 mL/min.
- The critical supersaturation of PbI₂ is unknown (no surface-tension data in the corpus).

### 10.4 Accuracy of the Ω(1,1) fit and of D

**Literature.** Neufeld is accurate to better than 0.1% [@Neufeld1972empirical]. The dominant uncertainty is the pair potential [@Brandt1956calculation; @Jasper2014first; @Cambi1991generalized].

**Estimate.** The four-term fit in `PbI2HeDiffusivity.H` lies 2.0% below Neufeld at 356 K, 5.2% at 500 K, 7.6% at 700 K and 8.7% at 1000–1200 K. D is therefore up to 9.5% too high. In the transport-limited regime T_dep is an equilibrium quantity and does not depend on D. Only the relaxation length Ud/(4k_m), about 3 mm at 106 mL/min and 16 mm at 540 mL/min, scales as 1/D. A 9% error moves profiles by 0.3–1.4 mm, well below the 1 cm bins. It also changes N, and hence the supersaturation estimate, by 9%.

**Verdict.** Quantitatively minor, but it must be corrected and documented as a verification item. A potential uncertainty of ±15–30% is a larger but still second-order issue. p_v(T) and the kinematics dominate.

| Candidate | Controlling group | 5 mL/min | 106 mL/min | 540 mL/min | Important? |
|---|---|---|---|---|---|
| Soret (radial) | D α_T ∇T/T ÷ k_m | ~0 | 0.2% | 1% | no |
| Soret (axial) | α_T ∇lnT ÷ ∇ln c | 7% | 7% | 7% | minor |
| Non-dilute (front) | w at T_dep | 0.62 (T5_60) | 0.54 / 0.10 / 0.01 | 0.02 | yes for T106_6 and T5_60 |
| Nucleation | S_bulk | 1.01 | 1.2 | ~7 | likely at 540 |
| Ω fit | ΔD/D | 5–9% | 5–9% | 5–9% | minor (verification) |
| p_v data | RT²/ΔH per e-fold | 17–31 K | same | same | **yes** |

Two further omissions are less exotic than the four above but may matter more: the liquid branch above T_m = 683 K [@Chase1998nist], and wall chemistry on SS316L (Section 3).

---

## 11. Research gaps and how our paper addresses them

| Aspect | Lobresco et al. 2026 | This paper (planned) | Supporting literature |
|---|---|---|---|
| Carrier flow | prescribed fully developed laminar profile; species passive, no feedback on flow | compressible, variable-density He; physical expansion 2.6–3.9× in the deposition zone | [@Majda1985derivation; @Issa1986solution] |
| Species variable | molar concentration φ, flux −D∇φ | mass fraction, flux −ρD∇Y | [@Curtiss1999multicomponent; @Bird2002transport; @Liao2007generalized] |
| Deposition | Eq. 5 relaxation (A, k, T_dep); HKS with σ = 1 and illustrative C_e = 5 | Eq. 5 and HKS in flux form; Dirichlet limit; Damköhler analysis; C_e retired or reinterpreted | [@Persad2016expressions; @Bouvet2024multiscale; @Aoki1998vapor] |
| Wall treatment | volumetric source in a 0.1 mm layer, uniform A_I/V_I | face flux with exact gas–solid bookkeeping; mesh independence shown | [@VijayaKumar2021implementation; @Gutti2009thermophoretic] |
| Thermodynamics | offline GEMS/HSC gas–solid p_v(T) | tabulated → per-cell GEMS3K; liquid branch; dataset sensitivity | [@Kulik2013gemselektor; @Miron2023thermofun; @Chase1998nist; @Konings1996infrared] |
| Transport properties | Brandt + LJ + fitted Ω | Neufeld Ω; D bracket; optional Soret | [@Neufeld1972empirical; @Cambi1991generalized; @Rosner2000heavy] |
| Numerics / V&V | one mesh, Δt = 1 ms, tolerance 10⁻³; "assessment, not validation" | analytical, MMS and code-to-code verification; GCI; splitting-error study | [@Roache2002code; @Celik2008procedure; @Carrayrou2004operatorsplitting] |
| Validation | Liu reference + five runs, T_dep and profiles | same data, with uncertainty in T_w, D and p_v; profile metrics | [@Liu2025chemical; @Oberkampf2010verification] |
| Chemistry scope | PbI₂ only | PbI₂ → Pb–Bi–I (BiI₃, BiI) | [@Liu2025chemical; @VanHattem2023ternary] |
| Missing physics | Soret, non-dilute, nucleation not treated | quantified (Section 10) | Section 10 |

**Gaps and responses:**

- **G1. No verified, variable-density CFD of a TC/TSG column exists.** This paper provides the first, with a full V&V hierarchy [@Zvara1985simulation; @Liu2025chemical].
- **G2. The assumptions of the 0-D T_dep formula have never been tested.** We quantify where it fails: short durations, low flow with Pe < 1, and non-dilute pulses. We also explain why T_dep falls with duration and flow, which may reconcile the Helas and Eichler observations for the desublimation regime [@Helas1978investigation; @Eichler1982evaluation].
- **G3. HKS parameters are unconstrained and the wall layer is mesh-bound.** We use Damköhler scaling, a face flux and the Dirichlet limit. Test whether C_e > 1 survives consistent SI units and physical kinematics. This is an open question with the T-Flows authors, `notes-for-authors.md` §A–B; do not state it as a result until confirmed.
- **G4. No CFD–GEM coupling exists for gas–solid deposition.** We implement it, analyze the splitting error at the wall sink, and apply it to multiple species [@Marino2020multiphysics; @Scuro2024coupledb; @Nichenko2022msr].
- **G5. Surrogates are untested for high-temperature gas–solid equilibria.** Deferred to a follow-up or to a section with a staged benchmark; decision pending (Section 8).
- **G6. Missing physics is unquantified.** Section 10 quantifies it. Soret is negligible, non-dilute loading is case-dependent, and nucleation is likely at 540 mL/min.
- **G7. p_v provenance has never been propagated into T_dep.** Report T_dep bands for HSC, JANAF, Binnewies–Milke and IVTAN [@Binnewies2002thermochemical; @Gurvich1991thermodynamic].
- **G8. Nothing links column-scale physics to reactor-scale source terms.** Provide wall-flux and nucleation closures for SAS4A-FATE, TRAIL, SOPHAEROS and ROMs, a LESTO/CIRTEN interface [@Hua2020development; @Miyahara2023development; @Lorenzi2016podgalerkin].

---

## 12. Must-cite works per section of the future paper

- **Introduction:**
  - LFR/LBE safety: [@Alemberti2014overview; @Alemberti2021safety; @Gianfelici2025research; @Neuhausen2020radionuclide; @OECD2015handbook]
  - Source data: [@Neuhausen2006investigations; @Ivan2025assessment; @Hua2020development]
  - Speciation and experiments: [@Karlsson2020thermochromatographic; @Karlsson2021thermochromatographic; @Liu2025chemical; @Liu2025exploring]
  - Models: [@Zvara1985simulation; @Eichler1993complex; @Cousin2013modelling; @Kissane2006interpretation; @Epstein2020enhancement]
  - Coupling and acceleration: [@Marino2020multiphysics; @Scuro2024coupledb; @Nichenko2022msr; @Kulik2013gemselektor; @Prasianakis2020neural]
  - Precursor: [@Lobresco2026assessment]
- **Methods:**
  - Flow and formulation: [@Weller1998tensorial; @Majda1985derivation; @Curtiss1999multicomponent; @Bird2002transport]
  - Properties: [@Neufeld1972empirical; @Brandt1956calculation; @Svehla1962estimated; @Kestin1984equilibrium; @Wilke1950viscosity]
  - Soret (if included): [@Rosner2000heavy; @Zirwes2025assessment]
  - HKS: [@Schrage1953theoretical; @Persad2016expressions; @Vaartstra2022revisiting]
  - Wall treatment: [@VijayaKumar2021implementation; @Dorsmann2007general]
  - Thermodynamics: [@Kulik2013gemselektor; @Miron2023thermofun; @Chase1998nist; @Konings1996infrared]
  - Splitting: [@Valocchi1992accuracy; @Carrayrou2004operatorsplitting; @Sportisse2000analysis]
  - Discretization: [@Patankar1980numerical; @Sweby1984high]
- **Verification:**
  - Methodology: [@Roache2002code; @Ramoa2022semi; @Roache1994perspective; @Celik2008procedure; @Eca2014procedure]
  - Analytical cases: [@Shah1978circular; @Tandon2004extended; @Sankarasubramanian1973unsteady; @Taylor1953dispersion; @DiRonco2021eulerian]
  - Code-to-code: [@Niceno2005unstructured; @Lobresco2026assessment; @Kren2026sharp]
  - Numerics to justify: [@Sweby1984high; @Patankar1980numerical; @Burchard2003high; @Issa1986solution]
- **Validation:**
  - Data: [@Liu2025chemical; @Lobresco2026assessment]
  - Framework: [@Oberkampf2010verification]
  - Comparisons: [@Eichler1982evaluation; @Helas1978investigation; @Liu2025exploring; @Gouello2021interaction; @Maruyama1999vapor; @Barrett2000aerosol]
  - p_v sources: [@Binnewies2002thermochemical; @Konings1996infrared]

---

## 13. Corpus weaknesses (to fix before submission)

- **Content not verified beyond metadata:**
  - Kou 2027 (Cs release; possibly GEMS-based) and Bani-Melhem 2026 PNE (lead evaporation; fog claim unverified);
  - Li 2026 (LBE aerosol), Maugeri 2014, Muñoz-Bueno 2005, Rosner & Epstein 1968, Rosner & Arias-Zugasti 2007 (the size of the Soret shift), Sharipov 2025, Waterson 2007, Kim 2014;
  - Clément 2003 (its INIS abstract is mis-attached);
  - Scuro 2024b's CPU-cost section.
  
  Read these before citing specifics.
- **Preprints and gray literature:**
  - Bani Melhem 2026 (SSRN, the only PbI₂ fog evidence), Johnsen 2017, Silva 2024, Kren 2026 (self-citation);
  - GIF reports (Alemberti 2014/2020/2021), Bowsher 1985, Trivedi 2017, the ITHACA-FV software entry.
- **Standard references missing from the corpus:**
  - Wagner et al. 2012 (the second GEMS reference, cited by Lobresco), Deng et al. 2021 (Po transport in LFR) and Lu et al. 2025 (LBE impurity CFD), all three cited by Lobresco;
  - Karlsson 2021a (silver capture), the HSC Chemistry reference, Graetz/Lévêque originals and Motz–Wise 1960.
- **Data missing from the corpus:** no surface tension or nucleation data for PbI₂, no measured α_T or D for any metal iodide in He, and no mixed-convection study of dense vapor pulses in horizontal tubes.
- **Metadata to check:** "Dorsmann" versus "Dorsman"; "Sugimoto Jun/Jum"; print versus online years (Hua 2019/2020, Poonoosamy 2018/2021, Liao 2006/2007, Garrison 2011/2012, Kulik 2012/2013).
