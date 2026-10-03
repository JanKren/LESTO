# rhoFixedFlowFoam phase change: final implementation plan

Temperature model, Hertz-Knudsen-Schrage (HKS) with p_v tables, and optional
GEMS3K. The exchange is implicit and exactly conservative, and 008-010 still
reproduce their logs.

Scope: `applications/rhoFixedFlowFoam` in /home/jan/LESTO. It must build on
OpenFOAM v2412 locally and on v2606, with gcc 7.5 and C++17.

Prototype (grafted, built, run): `scratchpad/design/final/src/rhoFixedFlowFoam`
- Test cases: `scratchpad/design/final/cases/`
- Diff against the incremental prototype: `scratchpad/design/final/grafts-vs-incremental-v2.patch`

All paths are relative to
`/home/jan/.claude/tmp/claude-1000/-home-jan-LESTO/983db48b-a50c-42b7-84fe-2c298dbe71b5/`.

---

## 0. How this plan was assembled

| design | numerics | OpenFOAM | science | total |
|---|---|---|---|---|
| incremental | 7.5 | 8.5 | 6.5 | **22.5** |
| physics | 8.5 | 7.0 | 7.0 | **22.5** |
| architecture | 5.5 | 5.5 | 5.5 | 16.5 |

Incremental and physics tie. **The base is incremental.** It is the only design
built into the real solver, and its mock path is verified bit-identical. It also
decides the file and hook structure, which the milestones implement. Physics
supplies most of the numerics:
- the per-cell solver defect;
- the unrelaxed, tight final solve;
- the monotone guard;
- per-face deposit storage;
- the halfCell (Robin) option;
- overlap binning;
- the start-up checks;
- the Graetz verification.

Architecture supplies:
- the phase-resolved p_v table;
- the GEMS `frozen` mode and the RAII engine;
- the bridge-state stamp;
- T_dep `fromPeak`.

Science supplies:
- the milestone order;
- the parity checklist;
- the physical-carrier validation;
- p_i computed from the actual rho;
- both interaction-layer rules.

Five grafts were prototyped in this round, on top of the incremental v2 code,
and run on v2412:
- the per-cell defect in the ledger;
- no relaxation for paired species;
- a regime-aware exit;
- a guard re-solve;
- a binary ledger with restart booking, and a tight `<Y>Final` solve.

See section 10.

---

## 1. Decisions at a glance

| topic | decision | from |
|---|---|---|
| mock path | `constant/thermochemistryProperties` absent, or `model mock`: the existing code path, no new Info lines, `evaluateThermochemistry()` and `solveSolidSpecies.H` unchanged | incremental |
| exchange form | affine per interface element: `fvm::Sp(Sp,Y)` plus `Su`, with Su >= 0 and Sp >= 0; **never fvm::SuSp** | all |
| regime choice | local implicit predictor `Y_loc = (H+Su*)/(A+Sp*)` from `YEqn.A()`/`H()` at every corrector. Converged = residual < tol **and** 0 regime changes | incremental + numerics must-fix |
| bounds | monotone guard (IMPLICIT to CAPPED/NONE, then an unrelaxed re-solve). Last resort: same-cell conservative projection, booked. **Y is never clipped** | physics + incremental |
| relaxation | paired species are never relaxed (a configured factor is ignored, with a start-up warning) | numerics/OF must-fix |
| final solve | the final matrix is solved with relTol 0: `<Y>Final` if present (one extra tight solve), otherwise `<Y>` must have relTol 0 (FatalIOError otherwise) | physics |
| balance | per-cell defect R summed into the ledger, so the closure is round-off for ANY tolerance. Per-patch transport from `YEqn.flux()`. Never `fvMatrix::residual()` | physics (verified) |
| ledger I/O | `scalarIOList uniform/phaseChangeBalance`, NO_WRITE, written explicitly in BINARY. Restart books field-rounding differences (amended in section 16, item 4: registered AUTO_WRITE, BINARY forced, so every write of Time stores it) | incremental + new finding |
| condensate state | wall deposit per interface element m''_e [kg/m2] (I/O via the `mw_<solid>` WALL patch values; amended in section 16, item 1: kept per rank in `uniform/phaseChangeDeposits`, `mw_` is output), sample reservoir `Cs_<solid>` [kg/m3]. `Y_<solid>` becomes derived output | physics/architecture |
| partial pressure | p_i = rho_c Y_c R T_c / M_i with the carrier's actual rho (rho = 1 in paper mode) | science must-fix |
| interaction layer | `faceCells` (default, 1 layer, A_I/V_I = 17681.9 1/m) or `distance` (centre distance <= delta, T-Flows, 2 layers, 8533.7 1/m). v_d = A V_I/A_I of the selected layer | science must-fix |
| wall law options | `interfaceTemperature cell|wall` (default cell); `wallResistance none|halfCell` (default none, halfCell opt-in, decided after M7; halfCell needs areaMultiplier 1 and areaPerVolume cell, section 17, item 12) | physics |
| p_eq | phase-resolved table, p_eq = min over phases, log10 p linear in 1/T (default). T-Flows table at its own nodes (bar, linear in T, ends held) for code-to-code. GEMS optional (`frozen` or `local`) | architecture + incremental |
| samples | `samples { <name> { pair PbI2_g; ... } }`, keyed to a pair. In-memory `cellBitSet::select` or an existing cellZone; never topoSet on symlinked meshes | incremental + OF must-fix |
| frozen carrier | `writeFrozenFields` default **yes** (unchanged output). `no`: read from `findInstance(mesh.dbDir(),"rho")`, T/p pre-registered NO_WRITE | incremental |
| outputs | `mDep_<solid>` [kg/m2] on the wall patches; axial profiles [mol/m] with exact overlap binning (native 424 bins plus 1 cm gamma-scan bins); profile `times`; T_dep interpolated, `fromInlet`/`fromPeak`/`searchFrom` | physics + incremental + architecture |
| GEMS build | opt-in: the `LESTO_GEMSBRIDGE` env var and the header present. `#ifdef LESTO_HAVE_GEMS`; stub FatalIOError; Allwmake stamp and wclean on toggle | incremental + architecture |
| code shape | one concrete class `LESTO::interfaceExchange` (no RTS) plus a thin include-hook in the species loop: a paired gas is solved by `solvePairedGasSpecies.H`, every other species by `solveGasSpecies.H`, which stays byte-identical to main (amended in section 17, item 8). Standalone headers for the kinetics and tables (g++-testable) | new (OF judge: "scope sprawl") |

---

## 2. Discrete equations

### Notation
All quantities are SI.
- Cell c: V_c [m3]; frozen rho_c [kg/m3]; Y_c = gas mass ratio (kg species per
  kg carrier; it can exceed 1); M = molar mass [kg/mol]; R = 8.314462618 J/(mol K).
- Time step: dt.
- Sign convention: S > 0 is net **evaporation** into the gas [kg/(m3 s)], and a
  wall flux mdot''_e > 0 is evaporation [kg/(m2 s)].

**Interface elements e** (built once). Each element has a cell c(e), an area
A_e [m2], a temperature T_e, a half-cell distance d_e and an output face f(e).
- `layer faceCells`: e = each face of the listed wall patches, and
  A_e = |S_f| * areaMultiplier.
  - `areaPerVolume layer`: A_e is scaled by (A_I/V_I)/(A_w,c/V_c), which is
    about 1 on this mesh.
- `layer distance`: e = each cell whose centre is within delta of the wall
  patches (patchDataWave with wallPointData<label>, which also gives the
  nearest wall face f(e) for output; amended in section 26, item 1: the
  exact distance from the nearest face of all interface patches on all
  ranks, an octree search, because patchDataWave's near-wall correction
  depends on the order of the patches and on the decomposition). A_e = V_c (A_I/V_I) * areaMultiplier,
  with A_I = sum of the wall-face areas and V_I = the layer volume. halfCell is
  not allowed in this mode.

**Sample cells s** in Z_j of sample j (keyed to a pair): a_s [1/m], with
reservoir Cs_s [kg/m3].

### E1 Gas equation (paired gas species; implicit Euler; corrector k)

    rho_c V_c (Y_c - Y_c^n)/dt + sum_f F_f(Y) + V_c Sp_c Y_c = V_c (Su_c + Srel_c)

    YEqn  = fvm::ddt(rho,Y) + fvm::div(phi,Y,"div(phi,Y_i)") - fvm::laplacian(rhoD_i,Y);
    A_T = YEqn.A(); H_T = YEqn.H();          // transport only: for the predictor
    YEqn += fvm::Sp(Sp, Y);                  // diag += V*Sp   (verified v2412 fvmSup.C)
    YEqn -= (SuExch + Srel);                 // source += V*Su (verified v2412 fvMatrix.C)
    // no relax() for paired species
    YEqn.solve(mesh.solverDict(final ? "Y_iFinal" : "Y_i"));

Su and Sp are `volScalarField::Internal` [kg/(m3 s)]. Both are >= 0 on every
branch, so the exchange adds V*Sp >= 0 to the diagonal and V*Su >= 0 to the
source: it preserves the diagonal contribution of the transport matrix and
makes the source term Patankar-positive and L-stable. (Corrected in section
17, item 14: the matrix itself is an M-matrix only if the transport terms make
it one, e.g. upwind without a non-orthogonal correction; SuperBee weights the
downwind cell where its limiter is positive, so with SuperBee it is not.)
SuSp is not needed and not used: it would make evaporating cells explicit at
lambda*dt of about 6e3.

### E2 Unconstrained element law, affine in Y_c

mdot''_e(Y) = alpha_e - gamma_e Y_c, where beta_c = rho_c R T_c / M, so that
p_i = beta_c Y_c.

- **Temperature (paper Eq. 5, flux form, irreversible):**
  - alpha_e = 0 and gamma_e = rho_c v_d f(T_e), with f = max(0, 1 - exp(-k (Tdep - T_e))).
  - v_d is given directly, or v_d = A V_I/A_I of the selected layer. faceCells:
    A = 220 1/s gives 1.24421e-2 m/s. distance: 2.578e-2 m/s.
- **HKS (paper Eqs. 6-7, SI):**
  - G_e = kineticScale * Ce * 2 sigma/(2 - sigma) * sqrt(M/(2 pi R T_e)) [s/m];
    alpha_e = G_e p_eq(T_e) and gamma_e = G_e beta_c.
  - At 700 K: G R T/M = 89.65 m/s, and lambda = sum A_e gamma_e/(rho V) = 1.585e6 1/s
    in the WALL layer, so lambda*dt = 6.3e3 at 4 ms.
  - Y_eq = p_eq/beta_c is about 0.041 at 700 K and about 20 at 960 K.
  - kineticScale = 1e-5/sqrt(1000) = 3.16228e-7 reproduces the T-Flows
    expression (bar, g/mol).
- **halfCell option** (faceCells only; kinetics in series with half-cell diffusion):
  - HKS: G'_e = G_e / (1 + G_e beta_c d_e/(rhoD)_e).
  - Temperature: v'_e = 1/(1/(v_d f) + rho_c d_e/(rhoD)_e).
  - Here d_e = 1/deltaCoeffs of the face and (rhoD)_e = the face value of
    rhoD_i. This is algebraically an OpenFOAM mixed BC (verified to 1.3e-12 by
    the physics prototype). With `none` the law is evaluated at the cell centre
    (paper and T-Flows form).
- **Sample (inventory, HKS):**
  - Su*_s = a_s G_s p_eq(T_s) and Sp*_s = a_s G_s beta_s.
  - G_s uses the sample's own kineticScale and Ce.
  - Evaporation only: no deposition on the sample.

### E3 Bounds per element or sample cell (this step)

- **Wall element:**
  - upper_e = m''_e^n A_e/dt [kg/s], i.e. evaporation cannot exceed the deposit.
    A bare element has upper = 0 (deposit only).
  - lower = -inf.
  - Temperature model: upper = 0 (irreversible).
- **Sample cell:**
  - lower = 0 (no deposition) and upper = Cs^n V/dt.
  - Cs = 0 gives NONE.

### E4 Regime choice (every corrector, before the solve): local implicit predictor

For each exchange cell c, with transport-only A_T and H_T + Srel_c:
1. Start with all elements of c IMPLICIT.
2. Y_loc = (H_T + Srel + sum_K A_e alpha_e/V + sum_C upper_e/V) / (A_T + sum_K A_e gamma_e/V)
3. Compute each IMPLICIT element's flux at Y_loc.
   - Above upper: CAPPED if m'' > 0, otherwise NONE.
   - Below lower: NONE.
4. Repeat until nothing changes: at most n_e + 1 passes, which is 1 on this mesh.

The regimes contribute:
- IMPLICIT: Su += A_e alpha_e/V and Sp += A_e gamma_e/V;
- CAPPED: Su += upper_e/V and Sp += 0;
- NONE: nothing.

A bare wall deposits exactly when H/A > Y_eq; the sample evaporates exactly
when H/A < Y_eq. H/A is parallel-consistent, because `fvScalarMatrix::H()`
includes the processor coupling.

`nChanged` = the number of elements whose regime differs from the previous
corrector (at the start of a step: from the previous step), reduced over
ranks. **The outer loop is converged when initialResidual < tolerance AND
nChanged == 0.** (Amended in section 14, item 1: the changes must carry at
most guardTolerance x inventory; pending sign-off.)

### E5 Final solve and guard (after the loop)

1. Final solve.
   - If `<Y>Final` exists, one extra solve with it (relTol 0).
   - Otherwise the last solve is final, and `<Y>` must have relTol 0 (start-up FatalIOError).
2. Realised flux of each element from THAT matrix's coefficients and the final Y.
3. Guard.
   - If the total bound violation exceeds guardTolerance (1e-12) * inventory,
     switch only the violating IMPLICIT elements to CAPPED or NONE (monotone:
     at most #elements passes), then re-solve unrelaxed with the final dict.
     At most nGuard = 3 passes.
4. Anything left (round-off) is projected in the same cell:
   Y_c += (S'_c - S_c) dt/rho_c. Gas plus condensate stay exact.
   - There is **no clamp of Y** (negative Y is monitored, not clipped).
   - The projected mass is booked as a diagnostic.

### E6 Condensate update (realised exchange)

- m''_e^{n+1} = m''_e^n - dt * mdot''_e.
  - CAPPED: m'' := 0 exactly, and the round-off residue is booked as `clamped`.
- Cs_s^{n+1} = Cs_s^n - dt * S_s.
- Derived for visualisation: Y_s,c = (sum_{e in c} m''_e A_e + Cs_c V_c)/(rho_c V_c).
- `mDep_<s>` on face f = sum_{e -> f} m''_e A_e/|S_f|. This is identical to
  m''_e for faceCells with areaMultiplier 1 and areaPerVolume cell only
  (corrected in section 16, item 6: m''_e is per effective element area).

### E7 Sources

- **Release**, in the sample cells of a `release` sample:
  Srel = n0 M |[t^n, t^{n+1}] intersected with [t0, t0+t_ev]| / (t_ev dt V_Z).
  The sum over steps is n0 M exactly for any dt.
- **Inventory**, on a fresh start: Cs = n0 M/V_Z.
- **Removal**, at the first step with t^n >= t_r - dt/2: `removed += sum Cs V`,
  then Cs = 0. This happens exactly once; the flag is kept in the ledger.

### E8 Global balance per pair [kg], every step

- R = sum_c (b - A Y)_c is formed **per cell**, then summed (`matrixDefect.H`):
  - non-coupled patch: + bc - ic*Y;
  - coupled patch: + bc*pnf - ic*Y.
  (Section 21, item 1: the transport-only matrix plus the realised element
  flows; section 27, item 1: in flux form, every face flux once with
  opposite signs in its two cells, a processor face with the mean of its
  two sides, and the storage measured as G measures it, so that the
  double-precision imbalances of the operator are booked in solverDefect.
  Section 29, item 2: the face-flux correction of the matrix enters by its
  face values as well, and of the source only the explicit remainder E =
  b - fl(s + c) is kept, 0 for the solver's schemes.)
- F_p = sum over the faces of `YEqn.flux()` on non-coupled patch p (outward
  positive), then reduced. This needs `mesh.setFluxRequired(Y.name())`.
- Ledger entries:
  - initial, released;
  - transport_p (one per non-coupled patch, so INLET back-diffusion is separate);
  - removed, clamped (the round-off residue of the condensate: of a CAPPED
    element set to 0, and since section 27, item 1, the rounding of every
    condensate update);
  - solverDefect (sum R dt);
  - restart (sum of held_read - held_ledger at each restart);
  - heldAtLastStep.

Closure:

    closure = G + W + P - [initial + released - sum_p transport_p - removed - clamped - solverDefect + restart]

where G = sum rho Y V, W = sum m''_e A_e and P = sum Cs V. The closure is at
round-off for any linear tolerance (verified) and, since the M4 review,
round 2, for any run length (section 27, item 1: before, it drifted
linearly where the gas or a deposit stayed in place, to 7.6e-13 of the
inventory in 6000 steps). The **physical** defect
|solverDefect|/(initial + released) is gated separately; it needs the tight
final solve.

### E9 p_eq

- **Table:** per phase phi, log10 p_phi is linear in 1/T (Clausius-Clapeyron),
  with the end intervals extrapolated. p_eq = min over phi, so the 683 K kink is
  exact.
- **T-Flows emulation:** `linearT` (p linear in T), ends held, units bar, at
  the T-Flows nodes 283.15 + 10k K.
- **GEMS:** condensed DCs are suppressed with bounds (0,0), and
  p_eq = p_i^GEMS * 10^(-max_phi lgOmega_phi).
  - `frozen`: one call per element at start-up with a reference composition,
    validated against the table (exact for an ideal gas with pure condensates).
  - `local`: per step and element, warm-started (GEMSB_WARM_CELL), with guards
    x >= 1e-6, T >= 500 K and rc in {OK, OK_RETRIED}; otherwise the table.
  - Optional speciation factor chi = p_i^GEMS/p_i, lagged and clipped to (0,1],
    applied to gamma. It is a separately reported modelling change, off by default.

### E10 Profiles and T_dep

- Overlap weights w_cb = |[xmin,xmax]_c intersected with bin_b| / (xmax - xmin)_c,
  from the cell points (and the face points for elements).
- Line densities:
  - n'_g(b) = sum w rho Y V/(M dx_b);
  - n'_w(b) = sum w m''_e A_e/(M dx_b);
  - n'_s likewise.
- T_w(b) = the area-weighted wall-face T.
- T_dep:
  - Find the first bin with n'_w >= f * peak (f = 0.01), searching from
    `searchFrom` (`fromInlet`) or walking upstream from the peak (`fromPeak`).
  - Interpolate x_dep and T_dep linearly between the bin centres b-1 and b.
- Observable `wall` or `wallPlusGas` (gas left at shutdown condenses locally).
  The sample is always excluded: its reservoir n'_s, and (decided with the
  sign-off of M3, section 19, item 6; implemented in M5, section 20) the
  gas in the cells of every sample, so `wallPlusGas` = the wall deposit +
  the gas outside the sample cells.
- A search whose peak bin holds at most `significance` x the pair inventory
  (default guardTolerance) reports `none (insignificant)` (section 20).

### E11 Monitors, every step

- min Y and the negative mass sum min(rho Y, 0) V.
- Max carrier-validity mole fraction x = (Y/M)/(Y/M + 1/M_car), with M_car from
  `thermo.W()` / 1000 (kg/kmol to kg/mol). Mass evaporated while x > xWarn (0.05).
  (Amended in section 29, item 4: x = c/(c + p/(R T)), c = rho Y/M, the
  carrier's molar concentration from its frozen p and T; the same for a
  perfect-gas carrier, 6 to 25 times the old value in paper mode.)
- At start-up: max |T - p M_car/(R rho)|/T (frozen-T consistency) and the number
  of cells with p_eq > p.

### E12 Time step

The exchange imposes no stability limit. Accuracy sets dt:
- Co <= 0.5 with SuperBee: dt = 1 ms at 106 mL/min and 0.22 ms at 540 mL/min (corrected in review: Co 0.5).
- The 007 carrier is 28 mL/min (8.41e-8 kg/s) with Co_max 0.466 at 4 ms, so it
  is a development carrier only.
- dt study at 4/2/1 ms in M7.

---

## 3. Algorithm per time step

    startup:
      frozen carrier (writeFrozenFields; findInstance; T,p pre-registered when "no")
      species setup (unchanged)
      #include createPhaseChange.H:
        interfaceExchange phaseChange(mesh, thermo, rho, phi, species, names, state, molarMass)
          mock/absent -> nothing else, no output
          else: pairs; elements (faceCells | distance); samples (pair-keyed, disjoint from the layer)
                start-up checks (Euler ddt, Gauss convection, zero-flux wall BC of Y, relTol 0 on the
                  final solve, relaxation entry -> warning, states/molarMass, rho vs EOS; GEMS local
                  needs a physical carrier)
                setFluxRequired(Y); per-element G, p_eq table (T frozen), gamma, alpha
                mw_<s>, Cs_<s> READ_IF_PRESENT; inventory fill on a fresh start
                ledger READ_IF_PRESENT; fresh -> initial; restart -> book (held_read - held_ledger)
                GEMS (optional): engine per rank, frozen validation or local warm-state arrays
                master: postProcessing/phaseChangeBalance/<t0>/<pair>.dat
      WALL mask: model mock only (section 17, item 9; before: always, and a mesh without
        a patch WALL stopped every model that has a condensate)
    while runTime.loop():
      Info Time
      phaseChange.beginStep(): release rates (exact overlap), removal, bounds (m''^n, Cs^n),
                               GEMS local update, reset per-step counters, checks of the final
                               solver, the Euler ddt and the Gauss convection scheme
      for each gas species i (existing loop):
        pairi = phaseChange.pairOfGas(i)                // -1: solveGasSpecies.H, byte-identical
                                                        // >= 0: solvePairedGasSpecies.H (17.8)
        for corr = 0..:
          YEqn = ddt + div - laplacian
          if pairi >= 0:
            if no guard pass: phaseChange.linearise(pairi, YEqn.A(), YEqn.H(), Y)  // E4 predictor, nChanged
            YEqn += fvm::Sp(Sp,Y); YEqn -= Su
            solve (dict Y_i, or Y_iFinal on the final or guard solve); no relax
            F_p = flux() per patch; R = matrixDefect(YEqn, Y)
          else: relax; solve (unchanged)
          loop control: converged <=> res < tol && nChanged == 0; exhausted <=> corr == nCorr
                        after converged/exhausted: one Final solve if Y_iFinal exists
                        then guard (E5): violation > tol*inventory -> switch (monotone), re-solve, <= nGuard
        if pairi >= 0: phaseChange.realise(pairi, Y)  // E5 projection, E6 condensate, monitors, Info
        Info Y min/max (unchanged)
      if phaseChange.mock(): existing evaluateThermochemistry + WALL mask + solveSolidSpecies.H (verbatim)
      else: phaseChange.balance()   // E8 ledger, closure, one Info line per pair, balance.dat
      runTime.write()               // Y_* ; frozen fields only if writeFrozenFields yes
      if writeTime: c_<species> (unchanged)
      if !mock and (writeTime or t in profiles.times):
        mw_/Cs_/derived Y_s/mDep_ ; ledger.writeObject(BINARY) ; axial profiles and T_dep
        (amended in section 16, item 4: the ledger and the element state are
         registered AUTO_WRITE and refreshed in balance(); runTime.write() stores them;
         amended in section 18: mDep_ and exch_ are refreshed in balance() too, and
         phaseChange.writeProfiles() after the c_ output writes only
         postProcessing/axialProfiles, at the write times and the profile times)
      printExecutionTime

---

## 4. Dictionaries

### `constant/thermochemistryProperties` (new; absent = `model mock`)

    model            HKS;            // mock (default) | temperature | HKS

    pairs                            // one entry per GAS species that exchanges with a condensate
    {
        PbI2_g
        {
            condensed      PbI2_s;   // listed species with 'state solid'; its Y_ becomes derived output
            // molarMass taken from speciesTransportProperties (required for pairs)
            vapourPressure           // HKS (and fallback of GEMS)
            {
                file          "<constant>/pv_PbI2_phases.csv"; // '#' header names the columns; column 0 = T [K]
                phases        ("PbI2(cr)" "PbI2(l)");         // p_eq = min over the listed columns
                units         log10Pa;                        // log10Pa | Pa | bar
                interpolation logInverseT;                    // logInverseT | linearT (T-Flows)
                outOfRange    extrapolate;                    // extrapolate | hold (T-Flows) | fatal
            }
            gems                     // GEMS only
            {
                gas          "PbI2(g)";
                condensates  ("PbI2(cr)" "PbI2(l)");
                elements     { Pb 1; I 2; }
            }
        }
    }

    interface
    {
        patches              (WALL);   // by name (WALL is type 'patch'; wallDist is not used)
        layer                faceCells;// faceCells (1 layer = the paper's 'entirely within 0.1 mm') | distance
        distance             1e-4;     // [m], layer distance only (T-Flows: centre distance <= 0.1 mm)
        areaPerVolume        cell;     // cell (A_e = |S_f|) | layer (uniform A_I/V_I; paper Eq. 7)
        areaMultiplier       1;        // T-Flows emulation: 200
        interfaceTemperature cell;     // cell (paper/T-Flows) | wall (face T)
        wallResistance       none;     // none (paper/T-Flows) | halfCell (Robin-equivalent; faceCells only)
    }

    temperatureCoeffs
    {
        depositionVelocity 1.2442e-2;  // v_d [m/s]; give exactly one of depositionVelocity / A
        // A               220;        // [1/s]; v_d = A V_I/A_I of the selected layer (printed)
        k                  8e-3;       // [1/K]
        Tdep               680;        // [K]
    }

    HKSCoeffs
    {
        accommodation 1;               // sigma, 0 < sigma <= 1
        Ce            1;               // effective multiplier
        kineticScale  1;               // 1 = SI; 3.16228e-7 = T-Flows bar/(g/mol) expression
        equilibrium   table;           // table | GEMS (needs a build with LESTO_GEMSBRIDGE)
        speciation    none;            // none | lagged (GEMS only; separate modelling change)
    }

    GEMSCoeffs
    {
        system            "<constant>/gems/PbI2He-dat.lst";
        mode              frozen;      // frozen: start-up, validated vs the table | local: per step, warm start
        carrier           "He(g)";
        referenceComposition { He 1; Pb 1e-3; I 2.000002e-3; }  // frozen mode [mol]
        minMoleFraction   1e-6;        // local: below it the table (GEMS off by up to 5 decades at 1e-8)
        minTemperature    500;         // [K] local: below it the table (rc=3 at 300-490 K)
        maxLog10Deviation 0.01;        // reject and use the table if |dlog10| exceeds this
        updateInterval    1;           // local: steps between updates
        allowNonPhysicalCarrier no;    // local mode refuses rho != pM/(RT) by more than 1 % (paper mode)
        logDirectory      "gemsLog";   // processorN/gemsLog in parallel
        logLevel          4;
    }

    samples
    {
        boat
        {
            pair           PbI2_g;
            selection      { boat { action use; source cylinder;
                             point1 (0.04 0 0); point2 (0.05 0 0); radius 5e-4; } }  // in memory: 144 cells
            // cellZone    sample;     // alternative: an existing zone (never run topoSet on a symlinked mesh)
            mode           inventory;  // release (temperature model) | inventory (HKS)
            amount         4.84e-6;    // [mol] (T106_6: 2.23 mg)
            startTime      0;          // release [s]
            duration       6;          // release [s] (t_ev)
            areaPerVolume  4000;       // inventory [1/m] = 2/R_sample (T-Flows: x5)
            kineticScale   1;          // inventory
            Ce             1;          // inventory
            removeTime     6;          // [s]; default never
        }
    }

    balance
    {
        nGuard          3;
        guardTolerance  1e-12;         // relative to the pair inventory
        moleFractionWarning 0.05;      // passive-carrier validity monitor
        writeFile       yes;           // postProcessing/phaseChangeBalance/<t0>/<pair>.dat
    }

    profiles
    {
        axis     (1 0 0);
        origin   (0 0 0);
        bins
        {
            native { min 0; max 0.69; nBins 424; }   // = axial cells (1.62736 mm)
            gamma  { min 0; max 0.69; width 0.01; }  // gamma-scan comparison bins (exact overlap)
        }
        times      (7 12 16);                        // [s] in addition to write times
        onset      { fraction 0.01; search fromInlet; searchFrom 0; }  // fromInlet | fromPeak
                                                     // significance (section 20)
        observable wall;                             // wall | wallPlusGas (without
                                                     // the gas of the sample cells)
        // as implemented (section 18): search takes one mode or a list,
        // default (fromInlet fromPeak); observable default wallPlusGas
        // (section 13); min/max of a bin set default to the mesh extent
    }

### `system/controlDict`

    speciesTransport
    {
        tolerance          1e-7;  // unchanged meaning (outer initial residual)
        nCorr              3;     // unchanged; paired species also need 0 regime changes
        writeFrozenFields  no;    // NEW, default yes (unchanged output); no = read from findInstance(rho)
    }

### `system/fvSolution` / `fvSchemes` (paired species)

    solvers
    {
      Y_PbI2_g      { solver GAMG; smoother GaussSeidel; tolerance 1e-8;  relTol 0.1; directSolveCoarsest yes; }
      Y_PbI2_gFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0;   directSolveCoarsest yes; }
    }
    // no relaxationFactors entry for paired species (ignored with a warning)
    ddtSchemes       { default Euler; }                  // enforced for paired species
    divSchemes       { div(phi,Y_PbI2_g) Gauss SuperBee; }
    laplacianSchemes { default Gauss linear corrected; } // or 'limited corrected 0.5' if min Y matters

Per-case `0/Y_PbI2_g` settings:
- INLET: fixedValue 0 is the documented default. The INLET transport is booked
  separately.
- OUTLET: zeroGradient.
- WALL: zeroGradient (enforced).

---

## 5. File changes

| file | change |
|---|---|
| `rhoFixedFlowFoam.C` | Header rewritten (PURPOSE/PROGRAM FLOW). `writeFrozenFields`: `frozenWrite` (AUTO/NO_WRITE) and `frozenInstance` (`timeName` or `runTime.findInstance(mesh.dbDir(),"rho")`); T and p pre-registered (MUST_READ, NO_WRITE, `store()`) before `fluidThermo::New` when frozen. `#include "createPhaseChange.H"` before the WALL mask, which is built for model mock only (section 17, item 9); the species loop includes `solvePairedGasSpecies.H` for a paired gas, `solveGasSpecies.H` otherwise. `phaseChange.beginStep()` after the Time line. Existing mock block wrapped in `if (phaseChange.mock())` (re-indent only); solid solve only in mock. `phaseChange.balance()` before `runTime.write()` (ledger and closure, derived `Y_<s>`, `mw_`/`mDep_`/`exch_`, and the per-element state of `uniform/`, which `runTime.write()` stores; section 16.4); `phaseChange.writeProfiles()` after the `c_` output (M3: axial profiles and T_dep at the write times and the profile times, postProcessing only; section 18); `phaseChange.end()` after the time loop (end-of-run summary). EXE unchanged. |
| `solveGasSpecies.H` | Unchanged, byte-identical to main (section 17, item 8): the unpaired path (verified: filtered 008/009/010 logs identical to the pristine build). |
| `solvePairedGasSpecies.H` (new) | The paired-species branch, included from the species loop instead of `solveGasSpecies.H`: predictor, `fvm::Sp`/`Su`, no relax, Final dict, flux per patch, `LESTO::matrixDefect`, regime-aware convergence, guard re-solves, `phaseChange.realise()`. |
| `solveSolidSpecies.H` | Unchanged (mock only; header comment notes this). |
| `evaluateThermochemistry.{H,C}` | Unchanged (the mock model). |
| `createVariableDiffusivity.H` | `AUTO_WRITE` of D_/rhoD_ becomes `frozenWrite`. |
| `phaseChangeKinetics.H` (new) | Standalone, namespace LESTO, no OpenFOAM types. Contents: gas constant; `enum exchangeRegime {NONE, IMPLICIT, CAPPED}`; `depositionActivation`; `hksConductance`; `halfCellConductance`; `cellPredictor` (E4, per-cell element loop); `admissibleExchange`. |
| `vapourPressureTable.H` (new) | Standalone: '#' header CSV, named phase columns, `min` over phases, logInverseT/linearT, log10Pa/Pa/bar, extrapolate/hold/fatal. |
| `matrixDefect.H` (new) | Per-cell residual sum with coupled terms counted once. **Prototyped.** |
| `interfaceExchange.{H,C}` (new) | Concrete class `LESTO::interfaceExchange`: dictionary parsing; elements (faceCells / distance via patchDataWave<wallPointData<label>>; section 26, item 1: an exact octree search over the interface faces of all ranks); samples; per-pair Su/Sp Internal fields, regimes, alpha/gamma, upper bounds, m''_e, Cs; `beginStep`, `linearise`, `guard`, `realise`, `balance`, `write`; I/O of `mw_<s>` (WALL patch values; internal = area mean, display only), `Cs_<s>`, derived `Y_<s>`, optional `exch_<gas>`. State lives in plain lists; the fields are I/O vehicles only. |
| `phaseChangeLedger.{H,C}` (new) | `scalarIOList` NO_WRITE, `writeObject(IOstreamOption(BINARY), true)` at write times (section 16, item 4: a `binaryIOList`, AUTO_WRITE; item 5: `uniform/phaseChangeLayout`); layout per pair: base entries plus one per non-coupled patch; restart booking; Info balance line; master-only `OFstream` under `runTime.globalPath()`. |
| `axialProfiles.{H,C}` (new) | Overlap weights (cells, elements), several bin sets, `Pstream::listCombineReduce`, T_dep (fromInlet/fromPeak/searchFrom, interpolation), `times` list, integral check. The arithmetic (bin edges, exact overlap weights, onset search) is in the standalone `profileBinning.H` (g++-testable, `tests/testProfiles.C`); the reader/plotter is `utilities/axialProfiles.py` (section 18). |
| `gemsEquilibrium.{H,C}` (new) | All bridge code under `#ifdef LESTO_HAVE_GEMS`; otherwise a FatalIOError with a rebuild hint. RAII engine per rank, per-processor log dir, suppression bounds, frozen/local modes, guards, statistics reduced over ranks, recreate on GEMSB_ERR_FATAL, fallback counts. |
| `createPhaseChange.H` (new) | Constructs `phaseChange`; runs the start-up checks; prints the layer report (cells, faces, A_I, V_I, A_I/V_I, v_d, sample cells/volume). In mock mode prints nothing. |
| `Make/files`, `Make/options` | Add the .C files. Opt-in conditional `ifneq (,$(wildcard $(LESTO_GEMSBRIDGE)/gemsbridge.h))` adds `-DLESTO_HAVE_GEMS -I...` and `-L... -lgemsbridge -Wl,-rpath,...`. No default path. |
| `Allwmake` | Accept 2412/v2412/2606/v2606. Bridge-state stamp in `Make/$WM_OPTIONS/gems.state`: runs `wclean` on toggle; echoes the state. |
| `tests/` | `testPhaseChange.C` and `testVapourPressureTable.C` (g++ only); `regress.sh` (filtered full log vs pristine build and vs `out_excerpt`); `cmpFields.sh` (binary cmp of written fields); `phaseChange/{closedBox,channel,pipe,restartAscii,parallel,graetzWedge,layerMesh}` with an `Allrun` that asserts the milestone thresholds; `graetz.py`, `gci.py`, `summarise.py`. |
| `utilities/setFrozenCarrier/` (new app) | Paper-mode carrier: rho = 1, parabolic u from Q/A with an optional facet-area rescale, phi exactly divergence-free, T from a wall table. Not solver logic. |
| `thermochemistry/systems/make_pv_table.py` + `data/pv_PbI2_phases.csv`, `data/pv_PbI2_TFlows.csv` | Direct NASA-9 (Gurvich): per-phase cr/l columns every 5 K (270-1250 K) plus the T_m = 683 K node; T-Flows GEMS_VAPOR_PRESSURE at its 95 nodes (283.15 + 10k K, bar). |
| `run/012-...` to `run/017-...` (new) | 012 temperature model (007 carrier, dev); 013 HKS inventory (dev); 014/015 paper-mode temperature/HKS (code-to-code); 016+ physical carriers (M8). Each has a readme and an `out_excerpt` baseline. 008-010 untouched; 011 is the coupled-mock case of the main branch (section 31). |
| `readme.md`, manuscript `04-numerics.tex` | Documentation: models, sign convention, balance line, outputs, restart, GEMS build, traps. Manuscript: replace "fvm::SuSp for the sign-changing HKS term" with "fvm::Sp + Su (affine law, Su, Sp >= 0)"; add the defect, guard and final-solve text. |

---

## 6. Milestones

Every milestone ends with `tests/Alltest <milestone>`, which re-runs the M0 gate.

### M0: Regression gate and build (0.5-1 d)
Deliverables: `regress.sh`, `cmpFields.sh`, pristine reference build, Allwmake
(versions, stamp), opt-in `Make/options`, test scaffolding.

Acceptance:
1. 008/009/010 (10 steps): filtered logs identical to the pristine build (only
   dynamicCode case-path lines may differ). 010 `out_excerpt` identical (393
   lines). 008/009 differ from both builds only in the banner line renamed in
   e2ddbfa.
2. 010 with writeInterval 3: all 15 written files and `uniform/` are
   cmp-identical to the pristine build.
3. `ldd` shows libgemsbridge only with `LESTO_GEMSBRIDGE` set; toggling triggers
   wclean; 0 warnings.
4. Compiles on v2606 (Merlin7).

### M1: Frozen-carrier I/O (0.5-1 d)
Acceptance:
1. M0 gate (default `writeFrozenFields yes`).
2. 010 with `no`: write directories hold only `Y_*`, `c_*`, `uniform/` (at most
   6 MB instead of about 36 MB).
3. A restart from a directory holding only `Y_*` is cmp-identical at the final
   time to the continuous run, in serial and on 4 ranks.
4. The log prints the frozen instance.

### M2: Exchange core, temperature model, release, ledger (3-4 d)
Deliverables:
- `interfaceExchange` (faceCells), per-element m'', Cs;
- temperature model (v_d or A); release samples;
- ledger (per-cell defect, per-patch transport, binary, restart booking);
- guard and projection;
- start-up checks and monitors.

Acceptance:
1. Closed box (10 cells of 57.46 um, WALL on one side, phi = 0, T = 500 K,
   A = 220, dt = 4 ms): Y^n = Y0/(1+lambda dt)^n to <= 1e-14 relative
   (lambda dt = 0.671504); closure <= 1e-15.
2. Channel, 1000 to 400 K, release with t_ev/dt not an integer:
   released = n0 M to 1e-14.
   - Closure <= 1e-13 relative every step with a **loose** `Y` dict
     (tolerance 1e-8, relTol 0.1), with solverDefect reported.
   - With `Y...Final` at 1e-12: solverDefect <= 1e-11 relative.
3. Pipe: prints `A_I/V_I = 17681.9` and `v_d = 0.0124421`. 25 steps: closure
   <= 1e-13; min Y and negative mass printed.
4. 4 ranks: closure <= 1e-13; inventories match serial to <= 1e-9 (Final dict
   1e-12; with a converged outer loop, section 16, item 7). The test also
   prints `gSum(residual())` and asserts that it is NOT used.
5. Binary restart: cmp-identical. Ascii restart without writePrecision: closure
   <= 1e-13, and the `restart` term is non-zero and printed.
6. FatalError cases:
   - `backward` ddt;
   - fixedValue WALL for Y;
   - final dict relTol > 0;
   - sample overlapping the layer;
   - pair with a missing molarMass.
7. M0 gate.

### M3: Outputs and diagnostics (1.5 d)
Deliverables: axial profiles (native plus 1 cm), `times`, T_dep
(fromInlet/fromPeak/searchFrom), `mDep_`, derived `Y_<s>`, `exch_`.

Acceptance:
1. sum n' dx M equals the gas, wall and sample inventories to 1e-12 for both bin
   sets.
2. 4 ranks equal serial to 1e-12 (with a converged outer loop, section 16,
   item 7).
3. Synthetic linear-ramp deposit: T_dep matches the analytic crossing to 1e-9 K
   (both search modes).
4. Profiles are written at every `times` entry.
5. M0 gate.

### M4: Paper-mode carrier, distance layer, temperature code-to-code (2-3 d)
Deliverables: `setFrozenCarrier`; `layer distance` (patchDataWave, nearest-face
mDep mapping; section 26, item 1: the exact distance instead of
patchDataWave); run/014.

Acceptance:
1. Carrier: max |sum phi| per cell <= 1e-15 of the inflow (measured 4e-18).
2. Distance layer on the pipe: 30528 cells, `A_I/V_I = 8533.7`.
3. Fig. 9 identity n' t_ev = n0/U_b of the carrier used, within 1 % on the
   plateau, compared with their 3.38e-4 mol s/m (n0 = 3.16e-5 mol; U_b is
   pinned by the section 11 checklist).
4. Solid profiles at 2 t_ev for t_ev = 6 and 15 s agree within 1 %.
5. Table 2 temperature-based column: T_dep within +-2 K of the paper for both
   layer rules (the rule that matches is recorded).  (Restated by the
   researcher on 2026-10-02: |T_dep - T_paper| <= max(2 K, dT_bin), dT_bin
   the drop of the wall temperature across the native bin of the onset;
   section 26, item 5.)
6. M0 gate.

### M5: HKS with tables, inventory, removal (3-4 d)
Deliverables: HKS law, E4 predictor, CAPPED, guard, inventory and removal,
kineticScale/Ce/sigma, `make_pv_table.py`, run/013.

Acceptance:
1. Box, 700 K:
   - bare wall: Y1 = (Y0 + lambda dt Y_eq)/(1 + lambda dt) with lambda dt = 6241,
     to 1e-14;
   - CAPPED: Y_s = 0 exactly, closure 0;
   - ample condensate: Y1 = lambda dt Y_eq/(1 + lambda dt).
2. Channel with inventory and removal:
   - nCorr 3: closure <= 1e-13 and 0 projections.
   - nCorr 0: guard passes > 0 and 0 projections (verified: 5 passes, 0 projections).
3. Pipe, SI, SuperBee, 25 steps:
   - 0 regime changes in the last corrector and 0 projections (proposed
     rewording in section 14, item 1: 0 significant regime changes in the
     deciding corrector and in the final predictor, 0 significant
     projections);
   - closure <= 1e-13 and solverDefect <= 1e-11 (Final 1e-12);
   - sample monotonically decreasing.
4. Tables:
   - nodes exact;
   - melting point from the phase crossing 683.0 +- 0.1 K;
   - 5 K table vs NASA-9 < 2e-3 decades;
   - T-Flows table reproduced at its nodes to 1e-12.
5. M0 gate.

### M6: HKS code-to-code with T-Flows (2-3 d + one reference run)
Deliverables: an emulation case built from the parity checklist (section 11);
the reference T-Flows LESTO_HKS run by the researcher with the same pinned set,
or digitised Figs 11-12.  (Section 29, items 5 and 6: T-Flows zeroes the
sample at the start of step t_ev/dt, which `removeTime` t_ev - dt
reproduces; and its wall beyond x = 0.65 m is a 0 K Dirichlet wall, which
no option of setFrozenCarrier reproduces.)

Acceptance:
1. Conservation to round-off, with a total of MASS_SAMPLE/MM = 4.837e-6 mol.
   The stale 3.1669595e-5 constant printed by Source.fpp is not the inventory.
2. Gas and solid axial profiles vs the T-Flows run: L1 difference <= 5 % of the
   peak, peak position within 1 axial cell.
3. Fig. 12 front positions within 1 cell.
4. The same case in SI (`kineticScale 1`, one layer, multiplier 1) is reported
   as the "units effect" figure.

### M7: Verification suite (3 d)
Deliverables: closed box; Graetz-Robin (`graetz.py`, extended with axial
diffusion) on 3 wedge meshes; GCI; halfCell vs none; h vs h/2 first-layer mesh;
dt 4/2/1 ms; axial x2.

Acceptance:
1. Fitted decay rate mu:
   - GCI_fine <= 0.5 %;
   - the Richardson-extrapolated mu lies within GCI plus the reference's own
     error bound;
   - halfCell reported with its observed order.
2. Flux form: T_dep shift <= 1 bin and peak m'' within 5 % between h and h/2
   (the volumetric A form is reported, and expected to fail).
3. Deposit profile L1 change < 1 % between dt = 2 and 1 ms.

### M8: Physical carriers and validation (4-5 d + compute; needs measured T_w(x))
Deliverables:
- rhoSimpleFoam carriers for 5/106/540 mL/min (mass-flow inlet from STP, NIST He,
  measured T_w(x));
- the five Table 1 cases (n0, removal at the experiment end, both models);
- an automatic validation table.

Acceptance (code gates):
1. Carrier mass flow within 0.1 % of target and continuity <= 5e-7 per cell.
2. All five runs complete with closure <= 1e-12 and solverDefect <= 1e-8
   relative.
3. The table is generated: onset and peak position, peak height, integral/n0,
   T_dep, on the 1 cm n0-normalised profiles.

Science targets (reported, not a code gate): peak and onset within +-1 cm, and
T_dep within +-15 K.

### M9: Optional GEMS3K backend (2-3 d)
Acceptance:
1. Build without the bridge: `equilibrium GEMS` stops with a FatalIOError and a
   rebuild hint.
2. `frozen`: |dlog10 p| <= 1e-3 per element vs the 5 K table over 500-1100 K.
3. `local`:
   - 1000 to 400 K ramp: 0 failures; guarded cells counted.
   - <= 70 us per warm call.
   - 4 ranks give 4 engines and 4 logs.
   - Deposit profile vs the table <= 0.3 % (speciation none).
   - Paper mode refused unless `allowNonPhysicalCarrier`.
4. `speciation lagged` reported as a separate run.

### M10 (follow-on, not in the core scope): Pb-Bi-I
Several pairs and samples, per-element balance, BiI3 and Bi tables, optional
operator-split re-speciation, and qualification of gemsGuard at trace
composition. 5+ d.

Documentation (`readme.md`, `out_excerpt` of each new case, and the manuscript
§4 correction) is part of each milestone's deliverable.

---

## 7. Must-fix items and where they are handled

| must-fix (judge) | handled by |
|---|---|
| never `fvMatrix::residual()` in parallel (numerics, OF) | `matrixDefect.H`. Verified: 4 ranks give gSum(residual()) = 5e-5 to 1.1e-4 kg/s vs own 1e-10 kg/s |
| final corrector unrelaxed, relTol 0, flux and exchange from that matrix (numerics, OF) | no relax for paired species; Final dict or relTol 0 check; realise uses the last matrix (E5) |
| implicit branch choice; no acceptance with pending switches (numerics) | E4 predictor plus `nChanged == 0` exit plus guard |
| bounds conservative; no unbalanced clip (numerics, OF) | monotone guard re-solve, then same-cell projection (booked); Y never clipped; only the CAPPED residue is booked as `clamped` |
| Euler, zero-flux wall BC checks (numerics, OF) | start-up FatalErrors (M2 test 6) |
| full-precision ledger for any writeFormat/writePrecision (numerics, OF) | binary `writeObject` plus **restart booking of field rounding** (new; verified) |
| carrier and dt (numerics) | 1 ms default; M7 dt study; M8 regenerated carriers |
| outer loop: zero changes plus residual; monitor min Y (numerics) | E4/E11 |
| profiles on the real 0.69 m mesh, overlap, interpolated T_dep (numerics, science) | E10, M3 |
| kg balance only for mass ratio 1 (numerics) | same-formula pairs only; per-element balance in M10 |
| findInstance default (numerics, OF) | M1 |
| manuscript §4 SuSp text (numerics) | documentation item (M5) |
| no new Info lines in mock; full-log plus cmp gate; v2606 compile (OF) | M0 (verified on the grafted build) |
| no topoSet on the symlinked mesh (OF, science) | in-memory `selection`; cellZone only if it already exists |
| GEMS opt-in, stub, stamp, per-rank engine, fallback (OF) | M0/M9 |
| parallel hygiene (OF) | counters reset before a single reduce; master-only files under `globalPath()`; `listCombineReduce`; collective calls unconditional |
| samples keyed to a pair (OF) | `samples { <name> { pair ...; } }` |
| physical-mode validation milestone (science) | M8 |
| quantitative code-to-code plus parity checklist (science) | M4, M6, section 11 |
| p_i = rho Y R T/M with the actual rho (science) | E2 |
| both layer definitions (science) | `layer faceCells|distance` (M4) |
| T_dep definition, fromPeak, 1 cm bins, observable (science) | E10, M3 |
| inlet BC documented; inlet loss separate (science) | per-patch transport in the ledger (verified INLET line) |
| passive-carrier validity monitor (science) | E11 |
| what GEMS contributes (science) | `frozen` = consistency check; `speciation lagged` separate; `local` refuses paper mode |
| milestone order (science) | regression, I/O, temperature, outputs, temperature code-to-code, HKS, HKS code-to-code, verification, validation, GEMS. Verification comes before validation, so validation results are trusted |
| GCI-based verification (science) | M7 |

---

## 8. Risks

1. **Physical defect vs cost.**
   - Loose solves everywhere (1e-8, relTol 0.1) give solverDefect 6.8e-5 of the
     inventory in 10 pipe steps. Tight solves everywhere cost 7x.
   - Mitigation (verified): loose correctors plus a tight `Final` solve give
     4.4e-13 at 2.4 s/step instead of 6.6 s/step (serial, load ~30).
2. **Negative Y** from SuperBee plus the lagged non-orthogonal correction next
   to the sample.
   - Measured min Y: -6.6e-14 (tight everywhere), -2.4e-7 (loose plus Final),
     -3.5e-6 (loose).
   - Monitored, never clipped. If it grows: `limited corrected 0.5`, more
     correctors, or upwind as a bound check.
3. **Stiff active set.** An explicit pre-solve test dumps the whole sample (two
   designs verified this independently). The E4 predictor plus the guard were
   verified: 0 regime changes in the last corrector on the pipe; with nCorr 0
   the guard replaced all projections.
4. **Ascii restart rounding.** Fields at writePrecision 6 change the inventory
   by 1.75e-7 relative on restart. It is booked as `restart`, and binary is
   recommended.
5. **Passive-carrier validity.** Y is about 14 at the sample on the pipe (SI
   HKS). The frozen-carrier assumption fails there; the manuscript decision needs
   the E11 numbers.
6. **T-Flows parity is uncertain.**
   - Units and multipliers are inferred from their code.
   - Their source is Picard-lagged (explicit in b, iterated). Our implicit form
     equals their converged fixed point only if their outer loop converges;
     lambda dt is about 0.05 in their units.
   - Their inventory print uses a stale constant.
   - They use FLOW_RATE = 1.67e-6 m3/s.
7. **GEMS at trace composition**: silent errors of up to 4 decades at x <= 1e-8;
   rc=3 below about 490 K. Guards, `frozen` mode and table fallback mitigate it;
   Pb-Bi-I needs qualification.
8. **Axial resolution.** The HKS decay length is about 2 mm vs 1.63 mm cells, so
   the M7 refinement study is required.
9. **v2606 untested locally.** APIs to check at every milestone:
   `listCombineReduce`, `cellBitSet::select`, `fvMatrix::A/H/flux`,
   `setFluxRequired`, `findInstance`, `solverDict`, `writeObject(IOstreamOption,bool)`,
   `patchDataWave<wallPointData<label>>` and `basicThermo::lookupOrConstruct`.
10. **007 carrier.** It is 28 mL/min, and its U/T were overwritten in 3ce3791
    (T off pM/(R rho) by up to 0.91 K). Development only.
11. **WarningInFunction spam** when the outer tolerance is not met. For paired
    species this is replaced by a status in the exchange line plus an
    end-of-run count; the unpaired text is unchanged.
12. **Corner cells.** Per-element storage is exact on this O-grid (1 face per
    cell). At corners of other meshes, interfaceTemperature `cell` distributes
    uniformly.
13. The **GEMS warm state** is not persisted, so restarts in GEMS local mode are
    reproducible only to the GEMS tolerance.
14. **Cost at 1 ms.** About 2.4 s/step serial, i.e. 6e4 steps for 60 s, which is
    about 2-3 h on 20 ranks. T106_300 needs about 5x that.

---

## 9. Open decisions (recommendation first)

1. **Interaction layer**: `faceCells` (1 layer; the paper's 'entirely within
   0.1 mm') or `distance` (T-Flows, 2 layers; halves A_I/V_I). **Rec:**
   faceCells for physics; run both in code-to-code until the authors pin the
   rule.
2. **Area per volume**: `cell` (flux form) or `layer` (paper Eq. 7). They are
   identical to about 1e-6 on this mesh. **Rec:** cell; layer for emulation.
3. **Wall resistance**: `none` (default, paper) or `halfCell` (Robin,
   second-order in the stiff limit). **Rec:** decide after M7; expected halfCell
   for physical mode.
4. **Interface temperature**: cell (default) or wall. The difference is < 0.03 K
   here. **Rec:** cell; wall as a sensitivity run.
5. **HKS units**: SI (default) or kineticScale 3.16228e-7 plus x200/x5
   multipliers. **Rec:** SI for physics; emulation only for the code-to-code
   section, after the authors confirm.
6. **p_v data**: Gurvich phase-resolved (stable liquid above 683 K, 2-7x lower
   p_v in the evaporation zone) or the paper's solid-extrapolated curve and
   T-Flows table. **Rec:** Gurvich for physics; the T-Flows table for
   code-to-code.
7. **Role of GEMS for PbI2-He**: `frozen` consistency check (equals the table by
   construction), `speciation lagged` as a separate modelling change, `local`
   for Pb-Bi-I. **Rec:** table for production; frozen as a check; speciation
   reported separately; local only with guards.
8. **GEMS guards**: minMoleFraction 1e-6 and minTemperature 500 K, or invest in
   GEMS numerics (He-free subsystem, wider P grid) first. **Rec:** keep the
   guards; revisit with Pb-Bi-I.
9. **Condensate state**: per element m'' (plan) or per-cell Y_s (the incremental
   prototype). **Rec:** per element (a mesh-independent restart state, and it
   allows halfCell and wall T); Y_s stays as derived output.
10. **Sample models**: release (temperature) and inventory (HKS), as in the
    paper; deposition on the sample forbidden. **Rec:** keep; allow release with
    HKS for isolating deposition physics.
11. **Final solve**: one relTol-0 dict (simple) or loose `Y` plus tight `YFinal`
    (2.8x cheaper, verified). **Rec:** loose plus Final.
12. **Y >= 0 policy**: monitor only, or conservative redistribution. **Rec:**
    monitor only; fix the causes.
13. **Frozen output**: default yes (unchanged). **Rec:** yes by default; `no` in
    all new cases.
14. **Time step**: 1 ms (paper, Co <= 0.5 at 106 mL/min) or 4 ms (development
    on 007). **Rec:** 1 ms for results, 0.22 ms at 540 mL/min (corrected in review: Co 0.5), plus the dt
    study.
15. **Inlet species BC**: fixedValue 0, zero flux (T-Flows q = 0), or
    Danckwerts. **Rec:** fixedValue 0 for physics and a T-Flows match for
    code-to-code. The INLET loss is reported separately; it matters for T5_60.
16. **gamma-scan observable**: wall deposit, or wall plus gas at shutdown.
    **Rec:** wall plus gas (the gas condenses locally at shutdown); the sample
    is excluded.
17. **T_dep search**: fromInlet (paper) or fromPeak (main deposit). **Rec:**
    report both; fromPeak for T5_60.
18. **p_v > p above about 1110 K** (PbI2 boils): cap or warn. **Rec:** warn,
    no cap in paper mode; decide with the two-way-coupling question.
19. **Temperature-model v_d in physical kinematics**: explicit input or hidden
    T/273 scaling. **Rec:** explicit input and a documented recalibration.
20. **Carrier for physics runs**: regenerated per flow rate (M8) or 007.
    **Rec:** regenerate; 007 is 28 mL/min with overwritten T.

---

## 10. Prototype evidence

### New in this round
The incremental v2 code was grafted in `scratchpad/design/final/`, built on
v2412 with gcc 7.5 (0 warnings), and run.

- **Grafts**
  - `matrixDefect.H`: per-cell residual.
  - Paired species unrelaxed.
  - Exit rule: `res < tol && nRegimeChanges == 0`.
  - Monotone guard re-solve.
  - `<Y>Final` tight final solve.
  - Ledger: NO_WRITE, explicit BINARY write, and LEDGER_DEFECT, LEDGER_HELD and
    LEDGER_RESTART entries.
  - Per-patch transport line; negative-mass monitor.
  - Diff: `final/grafts-vs-incremental-v2.patch`.
- **Mock regression**: 008, 009 and 010 (10 steps) give filtered logs identical
  to the pristine build (256, 398 and 419 lines; only the dynamicCode case-path
  `link:` lines differ). 010 reproduces its `out_excerpt` (393 lines). 008/009
  differ from both builds only in the banner line (e2ddbfa).
- **Pipe HKS, loose solves** (GAMG 1e-8, relTol 0.1, relaxation entry 0.995
  present and ignored):
  - solverDefect -3.3e-10 mol after 10 steps (-6.8e-5 relative);
  - closure 1.1e-15 to 1.1e-14 relative every step;
  - 0.89 s/step.
- **Pipe HKS, tight everywhere** (1e-12, relTol 0): solverDefect 4.9e-14,
  closure 3.8e-16, min Y -6.6e-14, 6.6 s/step.
- **Pipe HKS, loose correctors plus `Y_PbI2_gFinal` 1e-12**: solverDefect
  4.4e-13, closure 2.3e-15, min Y -2.4e-7, 2.4 s/step (5 solves per step; the
  final solve takes 74 GAMG iterations).
- **4 ranks** (simple 2x2x1, loose solves):
  - gSum(fvMatrix::residual()) = 5.1e-5 to 1.1e-4 kg/s;
  - own per-cell defect -3.9e-8 to 3.3e-11 kg/s;
  - closure 0 and -1.9e-16 relative.

  This confirms the v2412 double counting of coupled terms. 0.28 s/step.
- **Ascii restart without writePrecision** (channel, HKS inventory with removal
  at 0.25 s, restart at 0.1 s):
  - binary ledger alone: the closure jumps to -3.9e-13 kg (-1.75e-7 relative),
    because the Y fields are rounded to 6 digits;
  - with restart booking: 'inventory re-read ... differs by -3.89925e-13 kg,
    booked as restart', and the closure stays at -9.2e-21 kg (continuous run:
    2.0e-21).
- **Guard** (channel, nCorr 0, 500 steps): 5 guard passes (2-4 cells each),
  0 projections, closure 2.5e-15. The incremental prototype projected Y in the
  same 5 steps (up to 1.9e-12 kg).
- **T-Flows LESTO_HKS source read**:
  - concentration in mol/m3;
  - p = c R T/1e5 (bar);
  - k_HKS = 2 sigma/(2 - sigma)/sqrt(2 pi R M[g/mol] T);
  - source explicit in b (Picard-lagged) with gas_limit/solid_limit clamps;
  - layer = centre distance <= 0.1 mm;
  - a_v x200 (wall) and 2/R_sample x5 (sample);
  - table at 283.15 + 10k K, linear in T, ends held;
  - MASS_SAMPLE 2.23e-3 g, but the conservation print uses the stale
    3.1669595e-5 mol;
  - FLOW_RATE 1.67e-6 m3/s.
- **v2412 APIs present**: `patchDataWave` and `wallPointData`.

### Inherited from the designs
- **Incremental**:
  - mock bit-identity (15 fields cmp);
  - 21/21 g++ kinetics and table checks;
  - closed-box closed forms (1e-15 to 1e-16);
  - 250-step pipe HKS with 0 regime changes;
  - `dictionary::set` 6-digit trap (1.47e-6);
  - GEMS gas-only p/Omega within 1.2e-4 for x >= 1e-6, 54.7 us warm;
  - analytic paper-mode phi divergence 4e-18;
  - mesh facts (A_I/V_I 17681.9, 144 sample cells, 30528 distance-layer cells).
- **Physics**:
  - Robin equivalence of halfCell (1.3e-12);
  - Graetz reference and wedge series;
  - adjacent-branch cycling found with explicit branches;
  - relTol 0.05 gives a 1e-7 per-step defect;
  - stored T equals the initial field (0.30 %);
  - chi = 0.36 at 1100 K.
- **Architecture**:
  - Make/options conditionals;
  - stale objects after a toggle (needs a stamp);
  - fluidThermo adopts pre-registered T/p;
  - 482/1890 silent GEMS errors at trace composition;
  - phase-resolved 25 K table error 9.7e-4;
  - T_m from the crossing 683.01 K.

---

## 11. T-Flows parity checklist (paper mode, code-to-code)

Pin each item before M4 and M6. Record it in each case readme.

**Carrier**
- rho = 1.
- U = 2 U_b (1 - r^2/R^2), with U_b = Q/(pi R^2) and Q = 1.67e-6 m3/s (T-Flows
  FLOW_RATE), or Q_std = 106 mL/min (paper). Record which.
- Facet-area rescale on or off.
- T from `wall_x_profile_complete.dat`.
  (Section 29, item 5: T-Flows reads the 66 rows of the count line, x <=
  0.65 m, and leaves the wall faces beyond at 0 K, a Dirichlet wall of its
  energy equation; `rows count` with `hold` gives 324 K there.)

**Species and transport**
- Species: concentration c = rho Y/M. With rho = 1, p_i = c R T.
- D: the published Omega fit, with p_ref = 1 bar in T-Flows vs p/atm in
  OpenFOAM (`pressureUnit`).
- SuperBee; dt = 1 ms; the inlet species BC as in T-Flows; the outer-iteration
  count of T-Flows (their Picard source equals implicit only at convergence).

**Layer and multipliers**
- Layer: `distance 1e-4` (2 layers, A_I/V_I = 8533.7), `areaPerVolume layer`,
  `areaMultiplier 200`.
- Sample: cylinder x in 0.04-0.05 m, r <= 0.5 mm (centre test, 144 cells),
  areaPerVolume 2/R_s = 4000, x5 (via `Ce 5` or areaMultiplier).
  (Section 22, item 5: on the pipe mesh the centres in 0.04-0.05 m lie on
  rings at r = 0.1, 0.265, 0.3, 0.436, 0.458 and 0.5 mm; the ring at 0.5
  mm (18 cells) is computed 5e-16 outside a radius of 5e-4.  "144 cells"
  is the strict test, 162 the inclusive one on the exact centres.  Pin
  which one T-Flows uses, and give a radius between two rings: 4.8e-4 for
  144 cells (run/013), 5.5e-4 for 162.)
- Removal (HKS): T-Flows zeroes the sample at the start of step t_ev/dt
  (section 12, item 9; `HKS:User_Mod/Beginning_Of_Iteration.fpp:78-94`, the
  step whose end Curr_Dt dt equals T_EVAPORATION), i.e. it evaporates over
  [0, t_ev - dt]; the solver removes an inventory sample at the step that
  starts at t >= removeTime - dt/2, so `removeTime` t_ev - dt reproduces it
  (section 29, item 6).

**HKS and vapour pressure**
- HKS: sigma = 1, kineticScale 3.16228e-7, M = 461.01 g/mol, R = 8.314.
- p_v: the T-Flows GEMS_VAPOR_PRESSURE column (bar) at its nodes, `linearT`,
  `hold`.

**Amounts and model parameters**
- n0: 3.16e-5 mol for Figs 8-11 (temperature model), 4.837e-6 mol for
  MASS_SAMPLE (HKS); Table 1 masses for the experiments.
- Temperature model: A = 220 1/s, k = 8e-3 1/K, Tdep = 680 K, and t_ev per
  figure.
  (Section 25, items 4 and 9: Tdep = 680 K holds for Figs 8-10 as the paper
  states it (section 26, item 6: their published positions correspond to
  about 640 K in this wall table, and both T-Flows temperature-model codes
  use 638.15 K: a question for the authors); the Table 2 runs use the
  case's T_dep^calc of Table 1 (775, 712, 666, 788, 671 K),
  with t_ev = 6 s and Table 1's n0.  The `areaMultiplier 200` above belongs
  to the HKS source of T-Flows; the temperature model uses 1 (section 15,
  item 3).  `distance 1e-4` is the T-Flows layer (`R - r <= DELTA`); the
  paper's text describes one layer (`faceCells`); run/014 runs both.)

---

## 12. Amendments from the adversarial review (binding for implementation)

### Blocking fixes

**B1. Per-patch lists in parallel.**
- Size every per-patch list (transport, ledger entries, reports) by the number of
  **non-coupled** patches only; processor patches always come after them.
- Alternatively, reduce each non-coupled patch separately with `gSum`.
- `Pstream::listCombineReduce` on a list sized by `mesh.boundary().size()`
  crashes with MPI_ERR_TRUNCATE when ranks have different numbers of processor
  patches.
- M2(4) and M3(2) must include `simple (4 1 1)` and `scotch` decompositions
  (if scotch is unavailable locally, use `hierarchical` with an uneven split).

**B2. Temperature model: never delete a deposit on a negative undershoot.**
- An element becomes CAPPED only when `upper_e > 0`. Otherwise an element whose
  flux would be positive becomes NONE, and m'' stays unchanged. The same rule
  applies in the guard.
- Assert that a CAPPED residue is round-off: `|m''^{n+1}| <= 1e-12 * m''^n`,
  otherwise FatalError.
- Box test: an existing deposit with a set-up that forces `H_T < 0` (negative
  predictor) must keep the deposit.

**B3. The `distance` layer needs per-element restart storage.**
- Store the state per element cell, e.g. a cell field `W_<s>` [kg].
- Keep face storage for `faceCells`.
- `mDep_` for distance mode: patchDataWave seeded with globalIndex face labels,
  plus a parallel map.
- `layer distance` itself belongs to M4, which is not in this round. Design the
  storage so it can be added without changing the faceCells restart format.

**B4. Profile times write only postProcessing files.**
- At profile `times` that are not write times, write only
  `postProcessing/axialProfiles`.
- Write `mw_/Cs_/Y_s/mDep_` and the ledger only when `runTime.writeTime()`.

**B5. Persist the regimes.**
- Write the per-element regimes as a `labelIOList` under `uniform/`, together
  with the ledger, and restore them on restart.
- A binary restart must then be cmp-identical (M2(5)).

### Adopted non-blocking points

1. **Bound the `restart` ledger term.** In binary format it must be ≤ 1e-14
   relative, otherwise FatalError. In ascii, warn above 10^(1 - writePrecision).
2. **Fresh start** means no `Cs_`/`mw_` state fields. If the state fields exist
   but the ledger is missing, stop with FatalError; never re-fill the sample
   silently.
3. **Complementarity check after the final solve and guard.** Count and report
   the mass of: NONE bare-wall elements with final Y > Y_eq (missed deposition),
   and CAPPED elements whose implicit flux is now within bounds.
4. **No clipping of Y < 0** anywhere. Remove the prototype clamp in
   `realisePhaseChange.H`. "clamped" means only the CAPPED round-off residue.
5. **Step order.** Derive `Y_<solid>` in `balance()` before `runTime.write()`,
   so that `c_<solid>` is consistent.
6. **Time step at 540 mL/min:** about 0.22 ms for Co ≤ 0.5.
7. **GEMS.** Re-apply the species bounds and invalidate the warm states after
   an engine is recreated. `frozen` mode gets the same minTemperature/grid
   guards and table fallback as `local` (M9, not this round).
8. **Bridge toggle.** A generated `gemsConfig.H`, rewritten only when the state
   changes, so that a plain `wmake` picks up a toggle of `LESTO_GEMSBRIDGE`.
9. **T-Flows parity additions (M4/M6).**
   - a_v_wall uses the analytic 2πR L (8544.5 vs 8533.7 1/m);
   - R_GAS = 8.314;
   - the sample is zeroed at the start of step t_ev/dt;
   - x5 applies to the sample only.
10. **Layer-rule fragility.** Print both counts at start-up: centre distance
    R − r, and the patchDataWave distance.

---

## 13. This round (decided 2026-09-24)

**Scope:** M0, M1, M2, M3, M5. M4, M6, M7, M8 and M9 are the next round.
- `layer distance` (M4) and the GEMS backend (M9) are not implemented now. Their
  dictionary keywords may be parsed, but they stop with a clear
  "not yet implemented" message.
- The optional-bridge build plumbing (M0) and the stub are in scope.

**Defaults:** all section 9 recommendations are accepted:
- faceCells, cell area-per-volume, wallResistance none, interfaceTemperature cell;
- SI HKS units;
- Gurvich phase-resolved p_v;
- per-element deposit storage;
- loose correctors plus a tight Final solve;
- negative Y monitored, never clipped;
- writeFrozenFields yes by default, no in new cases;
- dt = 1 ms for results;
- fixedValue-0 inlet;
- observable wallPlusGas;
- both T_dep search modes.

**Git:** branch `phase-change`. **No commits**; the researcher commits.

**Working rules for implementers:**
- Edit only:
  - `applications/rhoFixedFlowFoam/**`
  - `thermochemistry/systems/` (make_pv_table.py and data)
  - new cases `run/012-*` and later
  - `.gitignore` entries
- Do not modify `run/001-010`, `manuscript/` or `papers/`.
- Do not restore `run/007/261/{T,U}`. It would break the digit-for-digit
  regression of 009/010, and it is the researcher's decision.
- Run cases only in a gitignored work directory
  (`applications/rhoFixedFlowFoam/tests/work/`), never inside `run/`.
  - Generate the mesh there once (gmshToFoam from `run/001.../pipe.msh`) and
    symlink it into test cases.
- Build the solver with `wmake` (EXE stays `$(FOAM_USER_APPBIN)/rhoFixedFlowFoam`).
  - Build the pristine regression reference from `git show main:` into
    `tests/work/pristine/`.
- Follow the existing code style: rich PURPOSE header comments, 2-space indent,
  the existing brace style and Info formats; 0 compiler warnings.
- Every milestone ends with `tests/Alltest <milestone>` passing, plus the
  M0 regression gate.

---

## 14. Amendments from the M2 review, round 1 (signed off by the researcher, 2026-09-25)

These points came out of the M2 verification and review. Each one is
implemented and tested (`tests/Alltest M2`), and each changes or sharpens a
binding statement above, so it needs the researcher's sign-off.

**M2 status: accepted (2026-09-25).** The researcher signed off sections 14-17,
including item 1 (settled regimes instead of `nChanged == 0`) and item 9
(the sample reservoir `Cs_<s>` moves to M5). Final pass: `tests/Alltest -fresh all`
gives 114 PASS, 0 FAIL, 0 warnings; M0 gate identical to main; v2606 compile
clean. Carried over to M3: restart-verdict wording for boundary parameters
(item 3 of the final pass), halfCell refusal explanation, stale section 5 row.

1. **Settled regimes (E4, section 1, section 7 "no acceptance with pending
   switches").**
   - The outer loop is converged when initialResidual < tolerance and the
     regime changes found by the predictor of that corrector carry at most
     guardTolerance x inventory: sum over the changed elements of
     |q_new(Y_loc) - q_old(Y_loc)| dt. This is the threshold at which the
     guard accepts a bound violation. It replaces `nChanged == 0`.
   - Reason: in cells without gas, Y is solver noise of about 1e-40, and its
     sign flips elements between IMPLICIT and NONE at every corrector.
     With the strict count, run/012 never converges: 1145 to 3177 changes
     per step carry about 1e-37 kg.
   - A **significant** event is an element whose own mass exceeds
     guardTolerance x inventory. The exchange line counts significant
     regime changes, projections and complementarity violations, both for
     the deciding corrector and, separately, for the predictor of the Final
     solve (the Final solve runs its own predictor). The end-of-run summary
     sums them.
   - Proposed rewording of M5 acceptance 3: "0 significant regime changes in
     the deciding corrector and in the final predictor, and 0 significant
     projections".
2. **Release clock (E7).** The overlap `|[t^n, t^{n+1}] ^ window|` uses the
   exact step times: t^{n+1} = t^n + dt, starting at the start time. The
   clock is kept in the binary ledger (a global `TIME` entry) and continues
   exactly after a restart. The time-directory names are rounded to
   timePrecision between write times (10.0002, 10.0004, 10.0007 for
   dt = 0.22 ms), so the old name-based overlap released 0.91 and 1.36
   times the correct amount in alternate steps. A restart whose start-time
   name differs from the clock by more than 1e-6 dt prints a warning. The
   balance files print the clock in their time column. M5's removal time
   test (`t^n >= t_r - dt/2`) must use the same clock.
3. **Fresh start after a release window has begun (E7).** On a fresh start
   (no ledger), a sample whose window begins before the start time stops the
   run with a FatalIOError that states the missed fraction of n0. On a
   restart the ledger has already booked that part.
4. **Binary restart bound (section 12, non-blocking item 1).** The bound
   stays at 1e-14 relative. The gas and wall inventories are now compensated
   (Neumaier) sums per rank, and the per-rank (sum, compensation) pairs are
   gathered and added in rank order (`compensatedSum.H`). They are accurate
   to about eps whatever the summation order, so a legitimate restart after
   decomposePar, reconstructPar or renumberMesh passes. Measured on the pipe
   with Y = 0.1 everywhere: the serial state at 0.008 s restarted on 4 ranks
   (hierarchical 2 2 1) differed by 1.54e-14 with plain sums, and the old
   build refused it with the FatalError; with compensated sums it differs by
   0 (M2.5j). As a side effect the closure of the 25-step pipe dropped from
   2.6e-14 to 1.8e-15 of the inventory.
5. **Regimes file (B5).** `uniform/phaseChangeRegimes` starts with a
   fingerprint per rank: the number of elements and a 62-bit hash of their
   cells, patch faces and face centres. Regimes are restored only onto the
   same elements. After a re-decomposition or renumberMesh with an equal
   element count they are discarded, with a warning.
6. **Reference of relative values (E8).** max(|initial + released|, |held|,
   cumulative inflow), with inflow = sum over the patches of
   max(-transport_p, 0). A pair fed only through an inlet no longer prints
   x/VSMALL.
7. **Final solver checked every step.** Whether `<Y>Final` exists and that
   its relTol is 0 are checked at every step as well as at start-up, since
   fvSolution may be modified during a run (runTimeModifiable). The
   exchange model itself (`constant/thermochemistryProperties`) is read at
   start-up only.
8. **Start-up warning:** an outer tolerance below the absolute tolerance of
   the `<Y>` solver is rarely reached, and the loop then runs all
   nCorr + 1 correctors.
9. **Sample reservoir `Cs_<s>`** (listed in the M2 deliverables) comes in M5
   with inventory samples. A release sample has no reservoir; it injects
   directly into the gas.

---

## 15. Amendments from the M2 review, round 2 (signed off by the researcher, 2026-09-25)

Implemented and tested in `tests/Alltest M2`; each sharpens a statement
above.

1. **Compensated cumulative ledger entries (E8).** The cumulative entries
   (released, transport per patch, solverDefect, clamped, removed,
   restart and the diagnostics) each keep a Neumaier compensation, stored
   in the binary ledger after the entries, so a binary restart continues
   them bit for bit. `supplied()` adds all terms compensated. With plain
   running totals the release of a 6 s window at dt = 1 ms missed n0 M by
   6.2e-14, and a 15 s window at 0.22 ms (136364 steps, M8 scale) by
   5.2e-13 with a closure of 1.3e-12, above the M8 gate. Compensated: 2e-16
   and 1.8e-14 (M2.2g), 4e-16 and 1.6e-13 (M2.2h). The ledger file of the
   round-1 build (without compensations) is refused with a clear message.
   The clock (`TIME`) stays a plain sum: it repeats Time::operator++.
2. **Release rate (E7).** Srel = n0 M overlap/(length dt V_Z) with
   length = (start + duration) - start, the window as the overlaps
   represent it, instead of the nominal duration. The overlaps telescope
   to this length, so the release is n0 M to round-off at any start time;
   with the nominal duration a 0.0198 s window at t = 1000 s released
   1.6e-12 too much (M2.2i).
3. **areaMultiplier applies to the temperature model (E2).** A_e =
   |S_f| areaMultiplier for every model, so the applied rate is
   areaMultiplier A (layer mean), and the start-up line prints it
   (M2.1l). Use areaMultiplier 1 with the temperature model; the factor
   200 of section 11 belongs to the T-Flows HKS source (to be confirmed
   for run/014, M4).
4. **Ddt scheme checked every step (section 3 start-up checks).** Like
   the final solver (14.7), the Euler ddt scheme is checked again at
   every step, because fvSchemes may be modified during a run
   (runTimeModifiable). Before, a change to `backward` went unnoticed and
   the closure left round-off (1.9e-2).
5. **Restart format with the collated file handler (section 12,
   non-blocking item 1).** The binary/ascii decision reads `data.format`
   of a decomposedBlockData container (its own `format` is always
   binary), or the header of its first block. An ascii restart under
   `-fileHandler collated` was refused with the binary bound (M2.5l/m).
6. **Deposit field every step (B4).** `mw_<s>` is refreshed every step
   (function objects may sample it) and written at write times only
   (M2.2j).
7. **Processor patches of the frozen carrier.** decomposePar writes a
   uniform processor patch as `uniform <value>` with about 10 digits,
   also in binary files. The solver evaluates the processor patches of T,
   p, rho and U once after reading the carrier, so rhoD is the same on
   both sides of a processor face. Before, the 4-rank channel closed only
   to 2.3e-13 (M2.4d). Serial runs are unchanged; decomposed runs of
   every model, mock included, change in the last digits (section 16,
   item 8). The same rounding would apply to `phi` if its processor-patch
   values were uniform; it is a face field and cannot be evaluated, so
   M4's continuity check (`setFrozenCarrier`) must be run on the
   decomposed case too.


---

## 16. Amendments from the M2 review, round 3 (signed off by the researcher, 2026-09-25)

Implemented and tested in `tests/Alltest M2`; each changes or sharpens a
statement above.

1. **Restart state of the condensate (section 1 "condensate state", E6,
   sections 3 and 5; B3, B5, NB1, NB2).**
   - The element deposits m''_e are kept per rank in
     `uniform/phaseChangeDeposits`: a binary list (raw block, exact) behind
     the element fingerprint of 14.5, written with the ledger.
     `mw_<solid>` becomes output only.
   - Reason: OpenFOAM writes a field or patch whose values are all equal
     as `uniform <value>` in text at writePrecision, also in a binary file
     (`Field::writeEntry`).  A rank with a single WALL face, common with
     scotch, stored its deposit with 6 digits, and its binary restart was
     refused with "do not belong to the same run" (relative 8.6e-11 on the
     3-rank test channel, M2.5o).
   - `mw_` is read back only when the deposits file does not fit the
     elements (decomposePar, reconstructPar, renumberMesh), with a warning
     (M2.5j/k/q).  An `mw_` that differs from `uniform/` is reported and
     ignored (M2.5n).
   - Fresh start (NB2) = no ledger and no condensate state (the deposits
     file or `mw_`).  The state without the ledger stops the run, and so
     does the ledger with neither the deposits of these elements nor all
     `mw_` fields (M2.5d/i).
   - B3 (`layer distance`, M4) adds its elements to the same per-rank
     lists; its display field is a cell field.  M5's sample reservoir
     `Cs_<s>` should be kept the same way (a per-rank binary list), since a
     cell field has the same `uniform` rounding.

2. **Exact or rounded restart state (NB1).**
   - The binary bound 1e-14 applies when everything re-read is stored
     exactly: the deposits of `uniform/` (or an `mw_` in binary), `Y_<gas>`
     in binary, and no rank whose `Y_<gas>` internal field or `mw_`
     interface patch is one repeated non-zero value (such a value may be a
     rounded `uniform` text entry, e.g. a rank of a single cell).
   - Otherwise the state is "rounded": the difference is booked, with a
     warning above 10^(1 - writePrecision), as for ascii (M2.5c/h/m).
   - The restart line names the sources and the verdict, e.g. "(relative
     0; gas: Y_PbI2_g written in binary; deposits:
     uniform/phaseChangeDeposits; exact)".

3. **Processor patches of a paired gas** are evaluated from the neighbour
   cells at start-up (as the carrier in 15.7).
   - A single-face processor patch of `Y_<gas>` is written as 6-digit
     text.  Its values enter the first matrix of a restart (the SuperBee
     limiter).  Without the evaluation, the 3-rank channel restarted from
     0.2 with exact deposits was not cmp-identical at 0.3 (gas 1.4e-16 and
     wall 7.1e-16 relative, 15 files differ); with it all 27 files are
     identical (M2.5o).
   - The processor patches of `mw_` are evaluated after every refresh
     (they stayed 0; output only, M2.4g).

4. **Writes started by function objects (B4, B5).**
   - Before: the ledger and the regimes were NO_WRITE and written by
     `phaseChange.write()` when `runTime.writeTime()`.  A function object
     that writes between two steps (`Time::writeNow()`,
     `Time::writeAndEnd()`, e.g. `runTimeControl` with `satisfiedAction
     end` or its `nWriteStep`) stored the fields without them, and the
     restart from that directory was refused.  Stopping a long run on a
     condition and continuing it later is the M8 use case.
   - Now the ledger, its layout, the deposits and the regimes are
     registered AUTO_WRITE, with BINARY forced in `writeObject`
     (`binaryIOList.H`), and refreshed every step in `balance()` before
     `runTime.write()`.  Every write of Time stores them with the fields;
     between two steps they hold the end of the last step, consistent with
     the fields.  B4 still holds: nothing is written unless Time writes.
   - `phaseChange.write()` is removed (section 3).  M3 adds its own hook
     for the profile times.
   - M2.5p: writeNow at 0.056 and writeAndEnd at 0.057; the restarts from
     both are cmp-identical at 0.2 to the continuous run.  The derived
     `c_<species>` fields are still written by the solver at the write
     times of controlDict only (unchanged from main); they are not state.

5. **Ledger layout.**  `uniform/phaseChangeLayout` names the pairs, their
   condensates and the transport patches.  A restart with another layout
   (pairs reordered, patches renamed or reordered with the same count)
   stops with a FatalError that lists both, instead of booking entries to
   the wrong pair or patch (M2.5r).  This matters from M10 (several pairs).

6. **`mw_` is per effective element area (E6 corrected).**
   - `mw_` holds m''_e per A_e = |S_f| areaMultiplier (times the scale of
     `areaPerVolume layer`), not per wall area.  `surfaceFieldValue
     areaIntegrate` of `mw_` on the wall is the wall inventory divided by
     areaMultiplier (M2.1m; areaMultiplier 2: exactly one half).
   - E6's "mDep_ ... identical to m''_e for faceCells" holds for
     areaMultiplier 1 and areaPerVolume cell only.  M3 computes
     `mDep_ = sum_e m''_e A_e/|S_f|` and must test it with areaMultiplier
     different from 1 (a trap otherwise: off by the factor 200 in the
     T-Flows parity layer of section 11).

7. **Serial/parallel agreement needs a converged outer loop (M2(4),
   M3(2)).**
   - The final solve converges only the linear system of the last
     iterate; its lagged terms (the SuperBee limiter, the non-orthogonal
     correction) depend on the decomposition.  Agreement to 1e-9 (M2(4))
     or 1e-12 (M3(2)) therefore holds only when the outer loop converges.
   - M2.4a-c reach 1e-9 with a forced 21 correctors per step (outer
     tolerance 1e-9 below the corrector tolerance): a test condition.
   - M2.4e: the pipe with run/012's own outer loop (nCorr 10, tolerance
     1e-5; it converges in every step but the first, which stops at an
     initial residual of 1.2e-5 in serial and on 4 ranks alike) on 4
     ranks, simple (4 1 1): gas 3.8e-11 and wall 1.4e-10 relative to
     serial.
   - M2.4f: the test channel's own loop (nCorr 3, tolerance 1e-7, never
     reached; it exits after 4 correctors at an initial residual of about
     1e-4): serial and 4 ranks differ by 7.0e-8 (gas) and 3.4e-7 (wall),
     asserted below 1e-5.  simple (1 4 1) gave 2e-6 (review).  The closure
     stays at round-off in all of them.
   - Proposed rewording of M2(4): "inventories match serial to 1e-9 when
     the outer loop converges (run/012 settings, Final dict 1e-12)"; the
     same condition for M3(2).

8. **15.7 applies to model mock too.**  The carrier's processor patches are
   evaluated before, and independently of, the model selection, so a
   decomposed mock run whose carrier has uniform processor patches (e.g. a
   plug flow) changes in the last digits compared with main.  It is a
   correctness improvement; the M0 gate (serial) is unaffected.  `phi`
   cannot be evaluated: with `writeFrozenFields yes` a parallel restart
   re-reads `phi` from the restart directory, where single-face processor
   patches are rounded, so bit-identical parallel restarts need
   `writeFrozenFields no` (the recommendation for new cases).


---

## 17. Amendments from the M2 review, round 4 (signed off by the researcher, 2026-09-25)

Implemented and tested in `tests/Alltest M2`; each changes or sharpens a
statement above.  Section 14's status note still applies: M2 is complete
only after the researcher signs off (or reverts) 14.1 and 14.9, and the
amendments of sections 15 to 17.

1. **`wallResistance halfCell` in parallel (blocking; section 1 "wall law
   options").**
   - `setCoefficients()` built the delta coefficients inside the loop over
     the interface elements.  Their first construction is collective in
     v2412 (`deltaCoeffs` -> `lduAddr` -> `polyMesh::globalData`), so a rank
     without WALL faces never called it and the run hung at start-up
     (channel on `simple (1 4 1)`, ranks 1 and 2 without WALL faces).
   - The delta coefficients (and `Cf`, `magSf`) are now obtained on every
     rank before the loop; the rule is stated in the parallel-hygiene
     header of `interfaceExchange.C`.  M2.4h: that decomposition ends and
     closes to 1e-13, run with a 300 s timeout so that a hang fails.
   - `rhoD = 0` on a wall face (`D 0`) gave a floating-point exception;
     it now means an infinite half-cell resistance, `v' = 0` (M2.1n).

2. **Restart with `writeFrozenFields yes` (NB1, M2(5), section 16 items 2
   and 8).**
   - With `yes` the carrier is re-read from the restart directory, where
     the solver wrote it; a uniform `rho` is a `uniform <value>` text entry
     at writePrecision also in a binary file.  The gas inventory
     sum rho Y V then changed by the rounding of `rho` (3.2e-8 on the box),
     and the exact-restart bound refused a legitimate binary restart as
     "not the same run".
   - When `rho` is read from the restart time, it is now part of the
     exactness decision, with the rule of `Y_<gas>` (binary, no rank with a
     repeated non-zero value), and the restart line names it
     (`carrier: rho written in binary, a 'uniform' value ...; rounded`).
   - M2(5) "binary restart: cmp-identical" holds for `writeFrozenFields no`
     (the recommendation for new cases, section 13), or for a carrier
     without rounded `uniform` entries (M2.5s: the channel with `yes` is
     exact and cmp-identical; the uniform-rho box is booked as rounded and
     closes to 1e-15).  Section 16.8 said this for parallel `phi` only.

3. **Restart record of the interface (E6, section 16 items 1 and 6).**
   - The deposits of a restart are m''_e per effective element area A_e
     (`areaMultiplier`, `areaPerVolume`) on the faces of the interface
     patches.  A restart with another `areaMultiplier` or `areaPerVolume`
     re-scaled the wall inventory by the ratio of the areas: in ascii the
     change was booked silently as `restart` rounding, in binary it was
     refused with the misleading "not the same run".
   - `uniform/phaseChangeLayout` now records the interface (`layer`,
     `patches`, `areaPerVolume`, `areaMultiplier`), and a restart with
     another one stops with a FatalIOError that names the change (M2.5t).
     The scalar is stored as a string of its exact decimal representation
     (a scalar dictionary entry passes through text at the output
     precision).
   - Chosen over storing the deposit mass per element [kg]: the `mw_`
     fallback after a re-decomposition is per area as well, and a change
     of the interface model in mid-run should be explicit, not silently
     mass-preserving.

4. **Restart record of the release samples (E7, section 14 items 2-3).**
   - On a restart, samples were not checked: changing the amount, window or
     pair of a sample whose window had begun, or adding one whose window
     had begun, was accepted, and the total release was no longer n0 M.
   - The layout now records every sample (`pair`, `amount`, `startTime`,
     `duration`, exact).  A restart stops (FatalIOError naming the change)
     when a sample whose window has begun is changed, when a sample is
     removed before its window has ended, and when a new sample's window
     began before the restart (M2.5u).  Samples whose windows lie after the
     restart may be added, changed or removed; a sample may be removed once
     its window has ended (M2.5v).
   - A layout of an earlier build without these records gives a warning
     (M2.5w).

5. **Final solver with patterns (section 1 "final solve").**
   - `dictionary::found()` matches patterns, so a single `"Y_.*"` entry was
     taken for a `Y_PbI2_gFinal` solver: one redundant solve and a "final
     predictor" per step, and a misleading start-up line.
   - A separate Final solve exists now when `<Y>Final` selects another
     `solvers` entry than `<Y>`, as `solution::solverDict()` selects it (the
     literal keyword, else the last matching pattern).  `"Y_.*Final"` after
     `"Y_.*"` is recognised; before it, it is shadowed by OpenFOAM's own
     rule, and the relTol error says so.  The start-up line names a pattern
     entry (M2.6j).

6. **Tests sharpened.**
   - M2.4a-c and M2.4e compare the serial and parallel inventories at every
     step (before: the last step only), and M2.4e runs the converged loop
     on `hierarchical (2 2 1)` as well as `simple (4 1 1)`.
   - The labels of `checks.py closure/defect` name the reference actually
     used, max(|initial + released|, |held|, inflow) (14.6).

7. **Observations for later milestones (no change now).**
   - Serial/parallel agreement with run/012's outer tolerance 1e-5 is
     about 1e-10 on 4 ranks, but reached 1.3e-9 at intermediate steps on an
     uneven 6-rank `hierarchical (3 2 1)` split (verifier, round 4).  The
     agreement is set by the outer tolerance (16.7); M3(2) (1e-12) needs a
     tighter outer loop, and M8 should state the tolerance it uses.
   - With a Final tolerance of 1e-12 the physical solverDefect accumulates
     with the number of steps (2.4e-10 relative after 6100 steps, 1.2e-9
     after 136364).  It is booked, so the closure stays at round-off, and it
     is below M8's 1e-8; M8 should report it next to the closure.

Items 8 to 14 come from the review of the final M2 pass.

8. **Code shape (sections 1, 3 and 5).**  The paired-species branch lives
   in its own fragment `solvePairedGasSpecies.H`, which the species loop of
   `rhoFixedFlowFoam.C` includes for a paired gas
   (`phaseChange.pairOfGas(i) >= 0`) instead of `solveGasSpecies.H`.
   `solveGasSpecies.H` is not touched: it is byte-identical to main
   (`git diff main -- solveGasSpecies.H` is empty), so no unpaired species
   can change.  The plan said "a paired branch in `solveGasSpecies.H`";
   sections 1, 3 and 5 now say what was built.

9. **WALL mask for model mock only (section 3).**
   - The mask of the cells next to the patch `WALL` (and the line `Cells
     adjacent to WALL`) was built for every model with a solid species,
     although only the mock path uses it.  A mesh whose interface patch has
     another name (`interface { patches (wall); }`) stopped with "boundary
     patch WALL was not found".
   - `createPhaseChange.H` is now included before the mask, which is built
     only for `phaseChange.mock()`.  The mock path is unchanged (M0 gate and
     M2.7 identical to the pristine build).  M2.1o: the box with the patch
     `wall` runs, follows the closed form and closes.
   - The log of the other models loses the line `Cells adjacent to WALL`.
     `run/012/out_excerpt` (line 45, `Cells adjacent to WALL: 15264`) is
     outside the files of this pass and was not edited: M2.8 compares the
     excerpt without that line, and the researcher may remove it.

10. **Conservative convection scheme (section 3 start-up checks).**
    - The ledger books the transport as face fluxes.  `bounded Gauss ...`
      subtracts `fvm::Sp(fvc::div(phi), Y)`, a cell source that is not a
      face flux: the run/012 pipe with `bounded Gauss SuperBee` closed to
      8.7e-12 in 5 steps (3.4e-16 with `Gauss SuperBee`).
    - The convection scheme of `div(phi,Y_<gas>)` must be `Gauss
      <interpolation>` (any interpolation scheme is a face interpolation);
      anything else stops with a FatalIOError, at start-up and at every step
      (fvSchemes may change during the run, like the ddt check of 15.4).
      M2.6k.

11. **Restart verdict with the values of the next steps (NB1, sections 16.2
    and 17.2).**
    - The verdict `exact` ignored values that the restart re-reads but the
      inventory does not contain: `Y_<gas>` on the non-coupled patches
      (fixedValue, inletOutlet) and, with `writeFrozenFields yes`, the
      patches of `rho`, and `T`, `p` and `phi` (including the processor
      patches of `phi`).  OpenFOAM writes a patch of equal values as
      `uniform <value>` at writePrecision.  An inlet value 1.234567891e-3
      (written `uniform 0.00123457`), or a carrier generated at 12 digits,
      gave `exact` although the files of the continuation differed.
    - They are part of the verdict now.  A rounded one leaves the inventory
      exact, so the bound 1e-14 still applies, but the verdict is
      `inventory exact, continuation rounded` and the line lists them
      (`rounded for the next steps: Y_PbI2_g ('uniform' on INLET)`).
    - Only the writer can tell whether a `uniform` entry lost digits (a value
      re-read always reproduces its own text).  The solver therefore
      records, per rank and before every write, the fields checked and
      their rounded entries in `uniform/phaseChangeLayout` (`uniformValues`;
      a value is rounded when its text at writePrecision does not read back
      to it).  A restart uses the record when it belongs to the files (on
      every rank, and the element deposits of `uniform/` fit the elements);
      otherwise (an earlier build, decomposePar, reconstructPar,
      renumberMesh) every repeated non-zero value counts, as in 16.2
      (`possibly rounded`).  An inlet value `1e-3` or `p 101325` is written
      exactly and stays `exact` (M2.5o, M2.5s).
    - M2.5x (gas patch), M2.5y (carrier patch); M2.5s now lists the rounded
      `rho` patches of the box; M2.5w strips `uniformValues` too.
    - Qualified in section 18, item 10: the verdict covers the patch
      values that are read back, not the boundary parameters, and patches
      whose value is recomputed when read are exempt.

12. **`wallResistance halfCell` with a scaled area (section 1 "wall law
    options", E2).**
    - The series law is a Robin condition per physical wall area
      `|S_f|`.  The element law is `g_e = A_e rho v'` with
      `v' = 1/(1/v + R)`, `v = v_d f(T_e)` and `R = rho d_e/(rhoD)_e`
      the resistance of the half cell.  With `A_e = m |S_f|`
      (areaMultiplier m, or the rescaling of `areaPerVolume layer`) it is
      `m |S_f| rho/(1/v + R) = |S_f| rho/(1/(m v) + R/m)`: the scale
      divides the half-cell resistance by m as well.  The half cell
      belongs to the mesh face (conductance `|S_f| rhoD/d_e` for any
      kinetic area), and the Robin condition with the kinetic area
      `m |S_f|` is `|S_f| rho/(1/(m v) + R)`.  The two agree only in the
      kinetic limit `R << 1/(m v)`, so scaling `A_e` is not
      Robin-equivalent.  (Corrected in section 18: an earlier version of
      this item gave the scaled law as `|S_f| rho/(1/(m v) + m R)`.)
    - The combination stops with a FatalIOError (M2.6l) until a scaling is
      decided (M7).
    - With halfCell the model line calls `areaMultiplier*v_d*A_I/V_I` the
      kinetic rate, and a line per pair prints the applied rate
      `sum_e A_e v'_e/V_I` with f(T) and the half cell (M2.1d: `v'/h` to
      1e-12).

13. **Keys of mode inventory in a release sample (sections 4 and 13).**
    `removeTime`, `areaPerVolume`, `kineticScale` and `Ce` belong to mode
    inventory (M5).  A release sample ignored them silently (a `removeTime`
    that removes nothing); it now stops with "belongs to mode inventory,
    which is not yet implemented" (M2.6m).

14. **M-matrix wording (E1).**  `Su, Sp >= 0` means that the exchange adds
    `V Sp >= 0` to the diagonal and `V Su >= 0` to the source: it preserves
    the diagonal contribution of the transport matrix and keeps the source
    Patankar-positive.  It does not make the matrix an M-matrix: that is
    decided by the transport terms, and SuperBee (downwind weights where
    its limiter is positive) is not one.  E1, `phaseChangeKinetics.H` and
    the readme are corrected.


---

## 18. M3 implementation decisions (signed off by the researcher, 2026-09-28)

M3 (outputs and diagnostics) is implemented on top of the accepted M2 code
and tested by `tests/Alltest M3` (`tests/profiles/Allrun`, criteria
M3.0 to M3.8).  The points below were decided during the implementation;
each sharpens or deviates from a statement above and needs the
researcher's sign-off.

1. **Profile times (E10; section 4 `times`).**
   - The profiles are written at every write of Time, also the writes of
     function objects between two steps (amended in section 19, item 1),
     and at the entries of `times`.
   - An entry t_i is written at the end of the first step whose end t1,
     on the exact clock of the ledger, satisfies t1 >= t_i - deltaT/2: the
     step end nearest to t_i, ties to the earlier one.  Entries on one
     step, or on a write time, give one set of files; the log line and
     the file header list them (`reason write time; profile time 0.1`).
   - A restart skips the entries with t_i - deltaT/2 <= its start (the run
     that wrote the restart state wrote them) and writes the others
     byte-identically to the continuous run (M3.4b).  Entries whose
     nearest step end is not after a fresh start (also an entry up to
     deltaT/2 after it; section 19, item 4) and entries after the end are
     never written.
   - Review item B4: at a profile time that is not a write time only
     `postProcessing/axialProfiles/<time>/` is written; no field, no
     `uniform/`, no time directory of the case or of a processor
     directory (M3.4a, M3.4c).  The directory is named after the step's
     time name; the header holds the exact clock.

2. **Onset search (E10; section 4 `onset`).**
   - `search` takes one mode or a list; the default is both,
     `(fromInlet fromPeak)` (section 13: "both T_dep search modes"); both
     are logged and written.
   - `searchFrom` [m along the axis] applies to fromInlet only.  Its range
     is the bins whose centre lies at or after searchFrom, and its peak is
     the maximum of that range, so an upstream spike before searchFrom does
     not set the threshold.  fromPeak takes the first maximum of all bins
     and walks upstream while the upstream neighbour reaches the
     threshold.  A bin exactly at the threshold counts as reached;
     fraction must lie in (0, 1].
   - An onset at the first bin of the range (the observable is at or above
     the threshold from the first bin on, so the onset lies upstream of
     the range) is reported as `none (notBracketed)`, a range without a
     positive peak as `none (noPeak)`; neither is extrapolated.  A
     bracketing bin without interface faces gives x_dep but no T_dep.
   - `wallPlusGas` counts all gas, also the source cloud around a sample;
     its T_dep is the deposition onset only once the gas has left the tube
     (section 19, item 6).

3. **Bin sets (E10).**  `{ min; max; nBins | width; }`, with min and max
   defaulting to the extent of the mesh along the axis.  A `width` that
   does not divide the range (to 1e-6) gives a shorter last bin, written
   with its own edges; a width at or above the range gives one bin, and a
   set of more than 1e6 bins stops the run (section 19, item 3).  Cells
   and wall faces are distributed by the exact overlap of their axial
   extent (from their points) with the bins; an
   item of zero axial extent goes into the bin that holds its position.
   The part outside all bins is summed and reported ('outside'), so
   binned + outside = inventory holds for any bin set, also one that does
   not cover the mesh (M3.1c).  The arithmetic is in the standalone
   `profileBinning.H` (g++ test `tests/testProfiles.C`, M3.0).  On the
   pipe the native set `{ min 0; max 0.69; nBins 424; }` is the axial
   cells to about 6e-10: gmsh wrote the axial node coordinates within
   9e-13 m of k 0.69/424 (23668 distinct values instead of 425), so a
   cell puts up to that fraction of its amount into the neighbouring bin,
   and the start-up line reports a round-off volume outside [0, 0.69].

4. **Files and log.**  `postProcessing/axialProfiles/<time>/<gas>_<set>.dat`
   with '# key values' header lines (pair and set, exact time, directory
   and reason, molarMass, axis and origin, nBins, the inventories and the
   binned and outside totals [mol], observable, fraction, one `Tdep` line
   per mode with status, T, x, peak, xPeak, threshold [, searchFrom]),
   then one row per bin, `x_lo x_hi x gas wall sample observable Twall
   wallArea` [m m m mol/m mol/m mol/m mol/m K m2/m], with 17 significant
   digits.  Two log lines per pair, set and profile time (the totals, the
   relative difference of binned + outside from the inventory, T_dep of
   each mode), and one start-up line; nothing without a `profiles`
   dictionary.  `utilities/axialProfiles.py` reads the files, prints T_dep
   (as written and recomputed from the columns) and the integrals, and
   plots with matplotlib on request.

5. **Sample line density.**  n'_s is written and checked, but it is zero
   until the inventory samples of M5: a release sample has no reservoir
   (section 14, item 9).  It uses the cell binning of the gas, which the
   gas checks exercise.

6. **Wall temperature per bin.**  The area-weighted boundary value of T on
   the faces of the interface elements, with the overlap weights
   (computed once, since T is frozen).  With `interfaceTemperature cell`
   the exchange uses the cell T, while T_dep reports the wall-face T (the
   paper's T_dep is a wall temperature).

7. **mDep_<condensate>** is written for every model other than mock (one
   more file per write directory).  Patch values:
   sum_{e -> f} m''_e (A_e/|S_f|), with the area ratio as one factor, so
   that areaMultiplier 2 gives exactly 2 mw_ (M3.5a).  Internal field: the
   deposit per wall area of the cell, (sum_e m''_e A_e)/(sum_e |S_f|).
   Refreshed every step in balance(), processor patches evaluated (M3.5c),
   never read.  With `areaPerVolume layer` the rescaled element areas sum
   to A_I, so areaIntegrate of mw_ on the wall differs from the wall
   inventory only where the deposit per effective area is non-uniform;
   M3.5b therefore uses a non-uniform gas, where mw_ misses the inventory
   by more than 1e-3 and mDep_ matches it to 1e-14.

8. **exch_<gas>** [kg/(m3 s)] is optional: `balance { writeExchangeField
   yes; }` (default no).  It is the realised exchange of the last step per
   cell volume, sum_{e in c} q_e/V_c, positive = evaporation; the release
   of the samples is not included (M3.6).

9. **Serial/parallel measure and settings of M3(2) (section 16, item 7).**
   - Measure: for every profile file and column (gas, wall, sample,
     observable), the largest difference over the bins relative to the
     largest value of the serial column (its peak); of every onset search
     the status and xPeak equal, T_dep to 1e-9 K, x_dep to 1e-12 m, the
     peak and the threshold to 1e-12 of the observable's maximum
     (tightened in section 19, item 2; before: T_dep to 1e-6 K only).
   - The per-bin profile is a stricter measure than the inventory.  On the
     pipe (background Y = 0.1, the run/012 release, 4 ms steps):
     - the forced 21 loose correctors of M2.4a-c (inventories equal to
       2e-12): gas profile 7.3e-9 (simple 4 1 1) and 3.4e-9
       (hierarchical 2 2 1) of its peak, wall 3e-13 to 7e-13;
     - run/012's own loop (nCorr 10, tolerance 1e-5, converged in every
       step): gas 1.2e-7 and 4.9e-7, wall 6e-11 and 2.5e-10, T_dep 5e-7
       and 1.4e-6 K;
     - converged to the floor of the linear solver: correctors to 1e-13
       with relTol 0 (GAMG stalls near 1.25e-14 on this mesh, so 1e-14 is
       unreachable), outer tolerance 1e-12, nCorr 40, Y_PbI2_gFinal
       1e-13.  24 to 29 correctors per step, serial and on 4 ranks, and
       the profiles agree to 6e-14 (gas) and 4e-14 (wall) of the peak,
       T_dep to 1e-11 K.  Cost: 35 s per step serial, 10 s on 4 ranks.
   - M3.2 uses the last setting.  The differences come from the lagged
     SuperBee limiter in the steep release region, which loose correctors
     do not converge.  A production run with run/012's settings is
     therefore reproducible between decompositions to about 1e-7 of the
     peak of the gas profile, not to round-off; M8 should state its outer
     tolerance (section 17, item 7).

10. **Restart verdict (carry-over (a) of the final M2 pass; section 17,
    item 11).**
    - The verdict covers the inventory, the internal fields and the patch
      values that are read back.  Boundary parameters that OpenFOAM
      writes as 'uniform' text as well (inletValue of inletOutlet;
      refValue, refGradient and valueFraction of mixed; gradient of
      fixedGradient) are not checked, and the restart line says so:
      `..., booked as restart; checked: the inventory, the internal
      fields and the patch values re-read, not other boundary parameters
      (e.g. inletValue)`.
    - A patch whose value is recomputed when read is exempt, in the
      writer's record and in the restart's own check (`valueReadBack()`):
      zeroGradient, fixedGradient (v2412 does not even write its value),
      mixed (evaluated from refValue, refGradient and valueFraction),
      symmetry, symmetryPlane, wedge, empty, and every coupled patch
      (cyclic; processor patches as before).  inletOutlet reads its value
      back and stays checked.  Before, a single-face fixedGradient or
      mixed patch was recorded as rounded and gave `inventory exact,
      continuation rounded` although the continuation was identical
      (M3.8).

11. **halfCell refusal (carry-over (b); section 17, item 12, corrected
    there).**  The element law with a scaled area is
    |S_f| rho/(1/(m v) + R/m), not |S_f| rho/(1/(m v) + m R): the scale
    divides the half-cell resistance by m.  The refusal stays; the code
    comment, the FatalIOError text, interfaceExchange.H and the readme
    give the correct law and why it is not Robin-equivalent.

12. **Section 5 row of `rhoFixedFlowFoam.C` (carry-over (c)).**  The hooks
    are `phaseChange.balance()` before `runTime.write()`,
    `phaseChange.writeProfiles()` after the `c_` output, and
    `phaseChange.end()` after the time loop (section 5 corrected).

13. **wmkdepend trap (build).**  wmkdepend of v2412 reads a comment as one
    token into a 16 KB buffer; at a longer comment it prints
    `wmkdepend: buffer full while scanning` and stops scanning the file,
    so the #include lines after it are missing from the dependencies and
    a changed header would not recompile its users.  The header comment
    of `interfaceExchange.H` passed 16 KB with M3 and is now split into
    two blocks; `tests/functions` counts `wmkdepend: buffer full` as a
    warning (criterion BUILD).  v2606 uses 64 KB.

14. **run/012** keeps its dictionaries: a `profiles` dictionary adds log
    lines, and M2.8 compares its `out_excerpt` as is.  M3.7 runs a copy
    with profiles for the excerpt steps and compares its log without the
    profile lines with the excerpt; the case readme gives the dictionary
    to add for the 16 s run.

15. **Integral check in the log** is a diagnostic (relative difference of
    binned + outside from the inventory, 1e-16 in all tests), not a
    FatalError.

M3 status: implemented and tested (`tests/Alltest M3`, and `tests/Alltest
-fresh all` for the M0 gate and M1 to M3); pending the researcher's
sign-off of this section.

---

## 19. Amendments from the M3 review, round 1 (signed off by the researcher, 2026-09-28)

The M3 verification passed.  The review found one major defect
(reported by two reviewers), one minor point and four nits.  All are
addressed and tested (`tests/Alltest M3`: new criteria M3.1d and M3.4d,
M3.0 and M3.2 extended).

1. **Profiles at every write of Time (major).**
   - Defect: the profiles were written at the end of a step when
     `runTime.writeTime()` was true.  A write that a function object
     starts between two steps (`Time::writeNow()`, `Time::writeAndEnd()`,
     e.g. `runTimeControl` with `satisfiedAction end` or `nWriteStep`)
     happens in `Time::run()` before the next step, and `Time::operator++`
     then resets the flag, so the solver never saw it: the state
     directory was written, its profiles were not, and a restart from it
     could not produce them either.  This is exactly the end of a
     stop-on-condition run (section 16, item 4; M8).
   - Fix: a small `regIOobject` (`interfaceExchange::profileWriter`,
     object `phaseChangeProfiles`) is registered on Time with AUTO_WRITE,
     like the ledger.  Its `writeObject()` writes the profiles of the
     state just written (reason `write time`) instead of a file of the
     time directory, and ignores a call when Time is not writing.  Every
     write of Time therefore has its profiles, on every rank (the
     reduction stays collective).  `writeProfiles()` after
     `runTime.write()` handles only the `times` entries reached by the
     step.
   - The profiles of one state (the ledger's clock) are written once.  If
     the per-step profiles of a `times` entry were written and Time then
     writes the same state, they are written again with the reason
     `write time; profile time ...`, as a write time of controlDict at
     that step would give (the log shows both).
   - At a write time of controlDict the profile log lines now come from
     `runTime.write()`, i.e. before the `Wrote c_...` lines.
   - Test M3.4d: the channel without regular write times, written by
     `runTimeControl` at 0.056 (`writeNow`) and 0.057 (`writeAndEnd`, a
     profile time too), serial, on 4 ranks and on 4 ranks with the
     collated file handler: profile directories 0.02, 0.056, 0.057 only,
     reasons `write time` and `write time; profile time 0.057`, time
     directories 0, 0.056, 0.057 only; the serial files are
     byte-identical to those of a controlDict write at every step, the
     collated ones to the uncollated ones; integrals to 1e-12.  The
     binary of the M3 round-1 sources fails it (no 0.056; 0.057 without
     'write time').
   - Not changed (the code path of main, M0 gate): the `c_<species>` fields
     are still written only at the write times of controlDict, not at a
     write of a function object between two steps.  They are derived
     output and no restart reads them.

2. **M3(2) gate (minor).**  `checks.py compareProfiles` compared only T_dep,
   to 1e-6 K.  It now also requires equal status and xPeak of every onset
   search, x_dep to 1e-12 m, the peak and the threshold to 1e-12 of the
   observable's maximum, and T_dep to 1e-9 K (about 1e-12 relative).
   Measured: T_dep 1.7e-12 K (simple) and 1.1e-11 K (hierarchical), x_dep
   1.7e-15 and 1.0e-14 m, peak 9.8e-16 of the maximum.

3. **Bin count (nit).**  `widthBinEdges()` gave 0 bins for a width of 1e6
   times the range or more (undefined behaviour afterwards).  It now
   gives at least one bin (a width at or above the range is the single
   bin [min, max]), and the constructor stops with a FatalIOError above
   1e6 bins (nBins or range/width), before the edges are built (a width
   in the wrong unit would overflow the bin count).  M3.0 tests widths of
   1e5 to 1e12 m over [0, 0.1] m (one bin each; 0 bins before).

4. **Profile times near a fresh start (nit).**  An entry up to deltaT/2
   after a fresh start has the initial state as its nearest step end and
   is skipped, like the entries before the start.  The start-up line said
   'at or before the start time skipped'; it now says `(N skipped: their
   nearest step end is not after the start time)`, and axialProfiles.H,
   the readme and item 18.1 say so.  Writing the initial state for such
   an entry was not added (the initial fields are known).

5. **Bins that cut cells (nit, coverage).**  New criterion M3.1d: the
   bin set cut { min 0.0003; max 0.1; width 0.0025; } of the channel (every
   profile time) and of the ramps (observables wallPlusGas and wall):
   every bin's gas, wall, sample and observable equal the exact overlap
   of the cell-aligned native profile, sum_c n'_c |c ^ b|/dx_b, to 1e-13
   of the column's maximum; Twall (the wall-area-weighted mean) to 1e-13
   relative and wallArea to 1e-13 of its maximum.  This checks the
   solver-side axial extents of the cells and of the wall faces, which
   M3.1a-c (binned + outside = inventory) cannot.

6. **wallPlusGas at intermediate times (nit).**  E10's 'the sample is
   always excluded' means the sample reservoir n'_s; the gas of the
   source cloud is counted.  wallPlusGas is the shutdown observable (gas
   left in the tube condenses locally), so its T_dep is the deposition
   onset only once the gas has left the tube; at t = 0.03 in the channel
   it lies on the sample box (949.5 K) while the wall column gives 684.9 K
   (model Tdep 680 K).  Decision: the observable is unchanged and the
   limitation is documented (axialProfiles.H, readme, item 18.2);
   `utilities/axialProfiles.py` prints the onset of the wall column next
   to a wallPlusGas file.  Option for the researcher: exclude the gas of
   the sample cells from wallPlusGas (a change of E10).

M3 status after round 1: `tests/Alltest -fresh all` (2026-09-25) gives
134 PASS, 0 FAIL, 1 NOTRUN (M0.4, the v2606 compile on Merlin7); 0
warnings; M0 gate identical to main.  The sources also compile and link
with 0 warnings against a local minimal v2606 build (Debug), where the
hook writes the profiles at the write times of controlDict as well
(runTimeControl is not part of that build); pending the researcher's
sign-off of sections 18 and 19.

**M3 status: accepted (2026-09-28).** The researcher signed off sections 18-19 and
decided §19.6: `wallPlusGas` shall EXCLUDE the gas in the cells of every sample
(wall deposit + gas outside the sample cells), so its T_dep is meaningful at any
time; at shutdown, with the sample removed, it equals the old definition.  This
E10 change is implemented as an M5 carry-over.

---

## 20. M5 implementation decisions (signed off by the researcher, 2026-10-02)

(Items 3, 5, 6, 7, 9, 10, 11, 12 and 19 are corrected or extended by the
M5 review, section 21; items 4, 5, 15 and 17 by its round 2, section 22;
items 5, 7 and 11 by its round 3, section 23; items 5, 7 and 11 by its
round 4, section 24: the per-rank lists moved to uniform/phaseChange/, the
geometry bounds are floored and limited, and acceptance 2 holds as written,
so the proposed rewording of item 11 is withdrawn for it.)

M5 (HKS with tables, inventory samples and their removal) is implemented on
top of the accepted M2 and M3 code and tested by `tests/Alltest M5`
(`tests/hks/Allrun`, criteria M5.0 to M5.8; M3.2c in `tests/profiles`).
The points below were decided during the implementation; each sharpens or
deviates from a statement above and needs the researcher's sign-off.

1. **HKS law (E2).**
   - Wall element: `G_e = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi
     R T_e))`, `a_e = A_e G_e p_eq(T_e)`, `g_e = A_e G_e beta_c` with
     `beta_c = rho_c R T_c/M` (the cell temperature: p_i is the partial
     pressure of the cell; T_e, cell or wall face, enters G and p_eq).
     halfCell: `G'_e = G_e/(1 + G_e beta_c d_e/(rhoD)_e)`, 0 for rhoD = 0
     (the Robin flux, M5.0; closed form M5.1e).  The halfCell restrictions
     of section 17, item 12 apply unchanged.
   - `HKSCoeffs` is optional with the SI defaults (accommodation 1, Ce 1,
     kineticScale 1, equilibrium table, speciation none).  `equilibrium
     GEMS` calls `gemsEquilibrium::requireBridge` (rebuild hint without the
     bridge, "not implemented yet" with it: M9); `speciation lagged` stops
     with "not yet implemented" (M9).
   - `vapourPressure` is required per pair for HKS, and so is its `units`
     (a wrong unit is a factor 1e5 or a misread logarithm, so there is no
     default); `interpolation` defaults to logInverseT, `outOfRange` to
     extrapolate.  With `outOfRange fatal` the element temperatures outside
     the table are counted on all ranks first, so every rank stops with the
     same message (the range and the number of elements).
   - Start-up lines per pair: the table (file as written, phases, nodes,
     range, options, the temperatures where the stable phase changes) and
     the HKS ranges on the interface (T_e, p_eq, Y_eq, G beta/rho, lambda =
     sum_e g_e/(rho_c V_c) of the layer cells, the elements outside the
     table and those with p_eq > p, reported and not capped: section 9,
     item 18); per inventory sample its mean lambda and the range of Y_eq.
     The HKS coefficient line comes before the pairs' table lines.

2. **Table (E9; `vapourPressureTable.H`).**
   - The column names are the `#` line right before the first data row, so
     a file may carry provenance comments above it.
   - A node returns exactly the value of the file (the value is converted
     to Pa once when read) in either interpolation and unit.
   - linearT with extrapolate continues the end line but not below 0;
     logInverseT refuses p <= 0 at a node.  A malformed file (non-increasing
     T, a missing or non-numeric value, an unknown column) stops the
     start-up with the file and line.
   - `crossings()` returns where the stable phase changes (bisection on the
     interpolated curves); for the Gurvich table 683.000002 K.

3. **Tables of the repository (`make_pv_table.py`).**
   - Only local sources: thermo.inp of `~/opt/gems-data` (the pinned public
     file of `fetch_nasa_thermo.sh` works as well) and the T-Flows source
     through a read-only `git show`; the NASA-9 reader of
     `make_gems3k_input.py` is reused.
   - `pv_PbI2_phases.csv`: log10 p [Pa] of PbI2(g) over PbI2(cr) and over
     PbI2(l), each over the whole range (crystal above T_m, liquid below:
     metastable extrapolations of the NASA fits; the gas fit, 300-6000 K,
     is extrapolated to 270 K), every 5 K from 270 to 1250 K plus the node
     T_m where the two NASA-9 curves cross, 683.000002 K (the plan's 683 K;
     the fitted records meet 2e-6 K above the nominal melting point), 12
     decimals.  Interpolated as the solver does it the table is within
     4.4e-5 decades of NASA-9 on a 0.1 K grid (M5.4; limit 2e-3).
   - `pv_PbI2_TFlows.csv`: the source tokens of GEMS_VAPOR_PRESSURE [bar]
     at its 95 nodes (283.15 + 10 k K, checked), unchanged (e.g. `1.06416e1`),
     so the nodes are reproduced exactly (M5.4; limit 1e-12).
   - (Section 21, item 7.)  The T-Flows source is read at the pinned commit
     35c1f6406e8e (the tip of advanced_scalar_model when the table was
     made; the csv header records it), not at the branch.  M5.4 is split:
     M5.4a tests the repository's tables through vapourPressureTable.H on
     any machine; M5.4b (regeneration byte for byte, NASA-9 comparison) is
     NOTRUN where thermo.inp or that commit is missing.

4. **Elements of inventory samples (E3, E6).**
   - One sample element per cell of an inventory sample, appended after the
     wall elements, each its own exchange cell (samples are disjoint from
     the layer and from each other), of the area `A_e = areaPerVolume V_c`
     at the cell temperature, with the sample's own kineticScale and Ce.
     The sample elements of a pair's samples are inactive (a = g = 0) for
     the other pairs.
   - Bounds: `lower = 0`, `upper = Cs V/dt`; CAPPED sets Cs to exactly 0
     with the residue booked as clamped and the round-off assertion of
     review item B2, like a wall element.  The regime machinery (predictor,
     guard, projection, complementarity, significance) is shared.
   - Only with model HKS: `mode inventory` with model temperature stops.
     The keys of the other mode are refused both ways (an inventory sample
     stops at startTime or duration, as a release sample at removeTime,
     areaPerVolume, kineticScale or Ce).

5. **Restart state of the reservoirs (section 16, item 1).**
   - `uniform/phaseChangeReservoirs`: a per-rank binary list, the
     fingerprint of the sample elements (number, hash of cell, sample name
     and cell centre; the name since section 21, item 8), then Cs of every
     sample element of every pair.  When every inventory sample has been
     removed, a restart needs neither the list nor `Cs_` (the reservoirs
     are 0; section 21, item 2).  It and
     the output field `Cs_<condensate>` exist only when inventory samples
     exist.  `Cs_` is read back only when the list does not fit (decomposed,
     reconstructed, renumbered), with a warning, and enters the restart
     verdict like `mw_` (M5.2f: re-decomposed onto 4 ranks, exact, relative
     0).
   - The regimes list now covers the wall and the sample elements behind a
     fingerprint of all elements, which equals the wall fingerprint when
     there is no inventory sample: files of M2/M3 runs still fit.  The
     deposits list keeps the wall fingerprint and format.
   - Fresh start = no ledger and none of deposits, `mw_`, reservoirs,
     `Cs_`; a `Cs_` field without the ledger stops the run (M5.5).

6. **Removal (E7, section 14.2).**
   - At the start of the step, in `beginStep()` before the bounds, when the
     step's start t^n on the ledger's exact clock satisfies t^n >=
     removeTime - dt/2: `removed += sum Cs V` (compensated over the ranks),
     Cs := 0.  The ledger entry `REMOVAL_DONE` of section 2 ("1 once a
     removal has happened") is a **count** of the removed samples of the
     pair (a pair may have several inventory samples); which sample was
     removed, and at which clock, is recorded per sample in
     `uniform/phaseChangeLayout`, and a restart checks that the count and
     the records agree.  Without removeTime a sample is never removed.
   - Tested at 0.25 s (M5.2a, exact reservoir, sample 0 afterwards), across
     restarts before, at and after it (M5.2c, cmp-identical, one removal),
     on 4 ranks (M5.2d) and at a rounded time name (M5.2e: dt = 0.22 ms
     from t = 10 s, the removal in the step starting at 10 + 5 dt).

7. **Restart record of inventory samples (section 17, item 4).**
   - Recorded: pair, mode, amount, cells, volume, centre (section 21,
     items 3 and 5), removeTime, removed, removalClock.  A restart stops
     when an inventory sample is new (its fill belongs to a fresh start),
     changed (pair, amount, cells exactly; volume and centre to 1e-12:
     section 21, item 3), dropped before its removal, or when any sample
     changes its mode (M5.5); removeTime and the kinetics may change.  A
     sample dropped after its removal keeps its record (section 21, item
     2).  A layout of an
     earlier build (no records) with an inventory sample now: stops (the
     sample is new).
   - `recordSamples()` now runs before the restart booking (after the
     clock), for all samples: a dropped inventory sample was otherwise
     reported as "the files ... do not belong to the same run" by the
     round-off bound.  For release samples only the order changes.

8. **Output.**
   - Exchange line: `(sample: implicit N, capped N, none N)` after the wall
     counts, only for a pair with inventory samples (other logs unchanged;
     run/012's `out_excerpt` is reproduced, M2.8).  The wall counts are the
     wall elements only.
   - Balance line and file: `sample` = sum Cs V; derived `Y_<condensate>`
     includes Cs/rho in the sample cells; `exch_` includes the sample
     flows; the profiles bin the reservoir as n'_s.

9. **Round-off of the stiff exchange (E8).**  (Corrected by section 21,
   item 1: the statement "no remedy" below was wrong.)  With lambda dt of
   1e3 to 6e3 the diagonal and the source of an exchange cell are that
   many times its storage term.  Booked as the defect of the ASSEMBLED
   matrix, their rounding (about eps lambda dt of the gas of the cell per
   step) differed from that of `q_e = a_e - g_e Y`, which updates the
   condensate, and near equilibrium the closure grew linearly with the
   number of steps (the closed box at 700 K: 3.4e-13 after 3 steps, +6.7e-14
   per step, 2.7e-12 after 50).  The ledger now books R' = the defect of
   the transport-only matrix plus the realised flows q_e plus the release
   (equal to R in exact arithmetic): the box closes to 2e-16 after 3 and
   4e-16 after 50 steps (M5.1a, M5.1f, gated at 1e-13).  What remains is a
   real imbalance of double precision, now booked in solverDefect: at a
   fixed point the linear solver no longer iterates, but fl(a - g Y) is not
   0 (section 21, item 1).

10. **Acceptance 1 (box).**  There is no input for an initial wall deposit
    (a fresh start has no condensate state, NB2), so the CAPPED and the
    ample-condensate steps restart the base box (Y = 0, bare wall, closure
    exactly 0) from a state whose deposit is set through `mw_PbI2_s`
    without `uniform/phaseChangeDeposits` (the documented fallback; the
    difference is booked as restart).  CAPPED: all elements CAPPED, Y1 =
    m''/(rho h) to 2.8e-16, deposit, `mw_` and `Y_PbI2_s` exactly 0, the
    closure of the step exactly 0 (M5.1b).  Ample: Y1 = lambda dt Y_eq/(1 +
    lambda dt) to 1.7e-16 (M5.1c).  Bare wall: lambda dt = 6240.88 at 4 ms
    (M5.1a); below Y_eq a bare wall stays NONE with no exchange and no
    deposit (exch_, mw_ and Y_PbI2_s exactly 0), and Y = Y0 to round-off
    (5.2e-16; not bit for bit, as claimed before: section 21, item 4)
    (M5.1d, review B2).

11. **Acceptance 2 (channel).**  "0 projections" holds literally with the
    solves of the prototype evidence (PBiCGStab 1e-15, relTol 0, no Final
    dictionary; nCorr 3 with tolerance 1e-12, i.e. 4 correctors, and nCorr
    0 with 4 guard passes), which the test uses.  With the loose correctors
    and the tight Final solve of the M2 channel, 1 step (nCorr 3) and 4
    steps (nCorr 0) of 500 had projections of at most 6.4e-19 kg, all below
    guardTolerance x inventory (2.2e-18 kg), where the guard by design
    leaves a violation to the projection (E5.4): "0 significant
    projections" holds for both settings.  Proposed rewording, as for
    acceptance 3: "0 significant projections".  One step of the nCorr 3
    run had 2 significant regime changes (1.9e-16 kg) in its fourth and
    last corrector (the loop is not converged by design with tolerance
    1e-12); they carry no projection and the closure stays 1.5e-15.
    Criterion M5.2j (section 21, item 10) now automates the rewording with
    the loose correctors and the tight Final solve, 500 steps, nCorr 3 and
    nCorr 0.

12. **Acceptance 3 (pipe).**  run/013 with its own solver settings at 1 ms
    does not bring gas to the wall within 25 steps (the wall inventory
    stays 0, every wall element NONE), so the test runs 25 steps of 4 ms
    with a background Y = 0.1, as M2.3: 9360 wall elements deposit, 0
    significant regime changes in the deciding corrector and in the final
    predictor, 0 projections, closure 3.6e-16, solverDefect 4.4e-12, the
    sample decreasing in every step (M5.3a).  Criterion M5.3c (section 21,
    item 10) runs run/013 at its own settings for 25 steps as well: the
    criteria hold literally, and the wall stays empty.

13. **Carry-over (a): wallPlusGas without the gas of the sample cells
    (section 19, item 6).**  Implemented as decided: the observable is
    `wall + gas - sampleGas` per bin, where sampleGas is the binned gas of
    the cells of every sample (release and inventory), a new last column
    of the profile files; in a bin without sample cells the observable is
    exactly `wall + gas` as before.  The reader checks the identity.
    Observation: in the test channel at t = 0.03 the onset of wallPlusGas
    is 949.9 K (949.5 K when all gas counted): the gas cloud lies mostly
    in the cells of the same axial slices around the sample box, which are
    not sample cells, so wallPlusGas is still not a deposition onset while
    the cloud is in the tube.  Options: keep (at shutdown, with the sample
    removed and the gas left, both definitions coincide); exclude the
    whole axial extent of the samples; or use `observable wall` at
    intermediate times (with the significance of item 14 the wall column
    of that example, 7e-22 mol, is insignificant; the review's 684.9 K was
    read from it).

14. **Carry-over (b): significance of T_dep.**  A search is `insignificant`
    when its peak bin holds at most `significance` x the pair inventory
    (gas + wall + sample, mol): the amount `peak x dx_peak`, compared after
    the no-peak test.  `onset { significance ...; }`, default the
    guardTolerance of `balance`.  The file writes `status insignificant`
    with peak, xPeak and threshold, and a header line `# significance <s>
    minimumAmount <mol>` that the reader applies (M3.0, M5.7).

15. **Carry-over (c): uneven split of the M3(2) gate.**  M3.2c: 3 ranks of a
    manual decomposition of the pipe by x (< 0.1 m, 0.1-0.25 m, >= 0.25 m;
    cell counts and processor patches per rank are printed), with the
    converged loop of M3.2a/b, to the same tolerances.

16. **Carry-over (d): profile times of a restart.**  The entries reached
    are recorded with the state (`uniform/phaseChangeLayout`, `profiles {
    written (...); }`, exact decimals, cumulative over restarts), updated in
    `balance()` before `runTime.write()`, where the entries reached by the
    step are now consumed (before: when their profiles were written; the
    files are the same).  A restart skips exactly the recorded entries and
    those at or before its start that were not written; every other entry
    is written at the first step end t1 >= t_i - deltaT/2 of the new run.
    Fresh starts and layouts without the record keep the rule of section
    18, item 1.  The start-up line says `(N skipped: written by the run
    that wrote the restart state; M skipped: at or before the start time,
    not written by that run)` (M3.4b updated; M5.8: restarts with deltaT
    0.2 ms and 2 ms of a run with 1 ms write every entry exactly once).

17. **run/013-hks-inventory.**  The 007 carrier (development only), the
    Gurvich table as a link `constant/pv_PbI2_phases.csv` to
    `thermochemistry/systems/data` (`makeCase` of the tests now links a
    case's `constant/*.csv`), an inventory boat of 4.84e-6 mol (T106_6,
    2.23 mg) in the 144 cells at x = 0.04-0.05 m, areaPerVolume 4000,
    removeTime 6 s, run/012's numerics, `writeFrozenFields no`; readme and
    a 5-step `out_excerpt` (M5.3b).

18. **v2606.**  The sources compile and link with 0 warnings against the
    local minimal v2606 build (Debug, dummy MPI), and that binary
    reproduces the box closed form (1.7e-16) and the inventory channel
    (closure 1.3e-15, the removal at 0.25 s).  M0.4 (the compile on
    Merlin7) remains NOTRUN.

19. **Tests changed by M5.**  M2.6f (model HKS and mode inventory are no
    longer "not yet implemented": model HKS without a table and mode
    inventory with model temperature stop with their own messages; layer
    distance unchanged), M2.6m (the message of the inventory keys in a
    release sample), M3.4b (the start-up text of item 16), M3.2c (new),
    M3.0 and M2.0 (extended: significance; HKS kinetics).

M5 status: implemented and tested (`tests/Alltest M5`, and `tests/Alltest
all` for the M0 gate and M1 to M3); pending the researcher's sign-off of
this section.

---

## 21. Amendments from the M5 review, round 1 (signed off by the researcher, 2026-10-02)

The M5 verification failed on one criterion (M5.1d recorded a false PASS),
and the review found two blocking and three major defects, two minor
points and three nits.  Each finding was assessed first; all were real and
are fixed and tested (`tests/Alltest M5`: new criteria M5.1f, M5.2g-j,
M5.3c, M5.9, M5.4 split into M5.4a/b, two more M5.5 refusals).  The
binary of the M5 round-1 sources fails M5.1a (closure), M5.1f, M5.2g,
M5.2h, M5.2i, M5.3b, two M5.5 refusals and M5.9 (item 12).

1. **Closure drift of the stiff exchange near equilibrium (major).**
   - Defect: `recordSolve()` booked R = the defect of the ASSEMBLED
     matrix, whose diagonal and source carry `g` and `a` (lambda dt times
     the storage term), while `realise()` changes the condensate by
     `dt q_e`, `q_e = a_e - g_e Y`: two roundings of the same stiff
     quantity.  At a steady state both repeat identically every step, so
     the closure drifted linearly, not as a random walk.  Section 20,
     item 9 called this unavoidable; it is not.
   - Fix (`matrixDefect.H`, `solvePairedGasSpecies.H`, `recordSolve()`):
     a copy of the transport-only matrix is kept before `fvm::Sp`/`Su` are
     added (only for a solve that may be recorded: the Final solve and the
     guard re-solves, and every corrector when there is no Final
     dictionary), and the ledger books per cell R' = defect of the
     transport-only matrix + sum over the cell's elements of
     `elementFlow(regime, a, g, upper, Y)` (exactly the value `realise()`
     computes) + `Srel V`, summed over the cells and ranks.  R' = R in
     exact arithmetic.
   - Only the booking changes: the fields of the 50-step box are
     cmp-identical to those of the round-1 binary (257 files), and the
     logs of run/012 and run/013 differ only in the digits of `solver
     defect`, `solverDefect` and `closure` (checked field by field).  The
     `out_excerpt` of run/012 (M2.8, M3.7) and run/013 (M5.3b) are
     therefore rewritten from the new build.
   - Result: the closed box at 700 K closes to 2.0e-16 after 3 steps
     (was 3.4e-13) and 4.0e-16 after 50 steps (was 2.7e-12, +6.7e-14 per
     step); the channel of M5.2a to 1.9e-16 (was 1.5e-15); the pipe of
     M5.3a is unchanged at 3.6e-16.  M5.1a and M5.1c are now gated at
     1e-13; new M5.1f: the box over 50 steps, closure <= 1e-13 in every
     step.
   - What remains, now visible in solverDefect: at a fixed point the
     linear solver no longer iterates (its residual is round-off), but
     `fl(a - g Y)` is not 0; the gas does not change while the condensate
     changes by up to eps lambda dt Y_eq of the cell's gas per step.  This
     is a real imbalance of double precision (no double Y makes a - g Y
     exactly 0), and R' books it as solverDefect: in the 50-step box the
     deposit drifts by 4.3e-11 of itself and solverDefect reaches 1.5e-11
     of the reference (M5.1f checks the deposit to 1e-10).  At M8 scale
     (lambda dt about 350 at 0.22 ms, 1.4e5 steps) the systematic part
     can reach about 5e-9 of the gas that sits near Y_eq in exchange cells
     for the whole run; in a transient the sign varies from step to step.
     Decision needed before M8 (not implemented): accept it (solverDefect
     is the honest record of it), or define the realised exchange of an
     exchange cell from its gas balance (the transport-only defect), which
     removes the drift but changes E5/E6 (the element split of a cell's
     flow, and an unconverged solve would then be booked as exchange).

2. **A removed sample dropped at a restart blocked every later restart
   (blocking).**
   - Defect: dropping an inventory sample after its removal is allowed
     (section 20, item 7), but `recordSamples()` rebuilt the layout's
     `samples` from the current samples only.  The record of the removed
     sample was lost while the ledger's `REMOVAL_DONE` still counted it,
     so the next restart stopped with "the ledger counts 1 sample
     removals ..., uniform/phaseChangeLayout records 0" (reproduced).
   - Fix: the records of removed inventory samples that are no longer in
     the dictionary are carried into the new layout for the rest of the
     run (log: "Sample boat (mode inventory, removed at t = ...) is no
     longer in the dictionary; its record is kept").  Listed again under
     its name, the sample is found removed and stays empty.  Since runs
     without inventory samples write neither `phaseChangeReservoirs` nor
     `Cs_`, a restart whose inventory samples are all removed now needs
     neither (the reservoirs are 0; the restart line says `reservoirs: 0,
     every inventory sample removed`).
   - The regimes warning of such a restart now names inventory samples
     added or dropped as a possible cause, besides another decomposition
     or a renumbered mesh.
   - Test M5.2h: the channel state at 0.3 (boat removed at 0.25) without
     the sample to 0.31, a second restart to 0.32, a third with the sample
     re-added to 0.33: 3 runs end, the record kept twice, 3 exact
     restarts, the re-added sample empty, closure 1.2e-17.

3. **The restart record compared V_Z exactly (blocking).**
   - Defect: V_Z is a plain sum reduced over the ranks, so its last bits
     depend on the decomposition; the record compared its exact decimal,
     and every restart after decomposePar, reconstructPar, renumberMesh or
     with another number of ranks was refused as "was changed"
     (reproduced on run/013: 8.117762652854522e-09 vs ...516), before the
     `Cs_` fallback of section 20, item 5 could apply.  M5.2f missed it:
     the channel's cells have equal volumes and the boat lay on one rank.
   - Fix: the record holds V_Z and the volume-weighted centre (item 5) as
     compensated sums over the ranks (compensatedSum.H: independent of
     the decomposition up to the round-off of the cell geometry), compared
     with a relative tolerance of 1e-12 (the centre: of the mesh extent);
     the number of cells, the pair and the amount stay exact.  A changed
     selection differs by at least one cell, many orders above 1e-12.  The
     V_Z of the release rate and the fill is unchanged (the plain sum), so
     no result changes.
   - Test M5.2g: run/013's state at 0.1 decomposed onto 4 ranks
     (hierarchical 1 2 2, the boat of cells of different volumes on all 4
     ranks), 2 steps, reconstructPar at 0.108 and 1 serial step: both
     restarts accepted, reservoirs from `Cs_PbI2_s`, exact; the recorded
     V_Z is 8.117762652854519e-09 in the serial, the parallel and the
     reconstructed run alike; closures 5.4e-16 and 3.6e-16.

4. **Vacuous exit-status checks in tests/hks/Allrun (major), and the
   false claim of M5.1d.**
   - Defect: M5.1d and M5.1e gated on `eval '[ $? -eq 0 ] ...'`; inside
     `check`, `$?` is the status of its own `shift`, so the checkers'
     results were ignored.  M5.1c computed its closure without gating it.
     Run properly, M5.1d FAILS: `Y_PbI2_g` is not Y0 bit for bit (up to
     5.2e-16 relative), as the verifier found; section 20, item 10, the
     readme and the Allrun header claimed it was.
   - Why Y is not bit for bit: on a uniform field the normalisation factor
     of OpenFOAM's residual collapses to round-off (normalised initial
     residual 1.5e-6 for a field that satisfies its equation), so
     DILUPBiCGStab iterates once per step and changes Y in its last bits.
     The exchange itself is exactly zero.
   - Fix: every checker's status is captured right after the call (`s1`,
     `s2`, ...), as the other criteria did.  M5.1d is restated as what
     holds and matters for review item B2: every element NONE, and the
     realised exchange `exch_PbI2_g` (now written in that case), the
     deposit and `Y_PbI2_s` exactly 0; Y = Y0 to 1e-14 (5.2e-16; 8 of 30
     values bit for bit); closure <= 1e-13.  M5.1c gates its closure at
     1e-13.
   - A copy of the Allrun with impossible limits (closure < 0 for M5.1a
     and M5.1c, 1e-17 for M5.1d, M5.1e and M5.1f) records exactly these 5
     criteria as FAIL (item 12).

5. **A moved selection was not named (minor).**  The record compared only
   the number of cells and V_Z; a selection moved to other cells of the
   same number and volume failed later as "the files ... do not belong to
   the same run".  The record now holds the volume-weighted centre
   (compensated, 1e-12 of the mesh extent), and the refusal lists
   `centre (...) m (written with (...))`.  New M5.5 refusal: the channel's
   boat box moved from x 10-20 mm to 30-40 mm (40 cells, equal volume).
   A cell-set hash was not used: cell centres and volumes differ in their
   last bits between decompositions, so an exact hash would refuse
   legitimate restarts, as V_Z did.

6. **Silently accepted inputs (nit).**
   - A fresh start at or after `removeTime - deltaT/2` of an inventory
     sample filled the reservoir and removed it in the first step without
     evaporation.  It now stops with a FatalIOError, like a release window
     that began before a fresh start (new M5.5 refusal: removeTime 0.4 ms
     with deltaT 1 ms).
   - The inputs of the other model (`HKSCoeffs` and a pair's
     `vapourPressure` under `model temperature`, `temperatureCoeffs` under
     `model HKS`) are named in one start-up warning.  A warning, not a
     refusal: keeping the coefficients of both models in one file to
     switch between them is the OpenFOAM convention for `<model>Coeffs`,
     while the sample keys of the other mode (refused) are never valid in
     a sample.  New M5.9.

7. **make_pv_table.py read a moving branch (minor).**  The T-Flows source
   is read at the pinned commit 35c1f6406e8e (the csv header records it;
   only that header line of `pv_PbI2_TFlows.csv` changed).
   `make_pv_table.py --sources` reports missing sources.  M5.4a tests the
   repository's csv files through vapourPressureTable.H on any machine
   (testVapourPressureTable with `-` for the NASA-9 reference: nodes
   exact, T_m, the T-Flows nodes); M5.4b regenerates both files from the
   sources byte for byte and compares with NASA-9 (4.4e-5 decades), and is
   NOTRUN when a source is missing.

8. **Fingerprint by the sample's index (nit).**  A sample element's
   fingerprint hashed its sample's position among all samples, so a
   release sample added in front of an inventory sample (allowed at a
   restart when its window lies after it) discarded `uniform/
   phaseChangeReservoirs` and `phaseChangeRegimes`.  It now hashes the
   sample's name (the wall fingerprint, and so the files of M2/M3 runs,
   are unchanged).  The reservoirs warning names inventory samples
   reordered or dropped as a possible cause.  Test M5.2i: the restart from
   0.2 with a release sample (window from 0.31 s) listed before the boat:
   no warning, reservoirs from the list, exact, and the fields at 0.25 and
   0.3 cmp-identical to the continuous run.  (A window from 0.3 s would
   not do: the exact clock of the step ending at 0.3 is
   0.30000000000000021, so that step releases.)

9. **The 117-column header line of rhoFixedFlowFoam.C (nit)** is
   re-wrapped.

10. **Points of the verification that were not failures.**
    - Acceptance 2 with the recommended solves: new M5.2j runs the
      inventory channel with the channel's loose correctors (GAMG 1e-8,
      relTol 0.1), the tight Final solve and outer tolerance 1e-7 for 500
      steps, with nCorr 3 and nCorr 0: closure 3.8e-16 in both, 0
      significant projections (1 and 4 steps with projections below
      guardTolerance x inventory), one removal, and 2 guard passes with
      nCorr 0, as the verifier found.  This automates the rewording "0
      significant projections" of section 20, item 11, which is still
      pending the researcher's sign-off.  Not gated, and printed: the final
      predictor has significant regime changes in 1 step (nCorr 3) and 58
      steps (nCorr 0: the single corrector is not converged, the Final
      solve's own predictor corrects the regimes), and nCorr 0 has 1 step
      with a significant complementarity violation; acceptance 2 does not
      ask for them (the final-predictor criterion of section 14.1 belongs
      to acceptance 3).
    - Acceptance 3 at run/013's own settings: new M5.3c runs run/013 for
      25 steps and checks the criteria of acceptance 3 literally (0
      regime changes of any size, 0 projections, closure 1.9e-16,
      solverDefect 1.2e-15, the sample decreasing; the wall stays empty:
      section 20, item 12).
    - The manuscript (defect, guard and final-solve text of sections 5
      and 7) is not edited by implementers (section 13); M0.4 (v2606 on
      Merlin7) remains NOTRUN.

11. **Tests changed.**  M5.1a (closure gate 1e-12 -> 1e-13), M5.1c
    (closure gated), M5.1d (restated, item 4), M5.3b (the regenerated
    excerpt), M5.4 (split into M5.4a/b), M5.5 (12 refusals);
    new M5.1f, M5.2g, M5.2h, M5.2i, M5.2j, M5.3c, M5.9; M2.8 and M3.7 use
    the regenerated run/012 excerpt.  `checks.py`: `unchanged` takes a
    tolerance and requires `exch_` = 0; `approach` takes `deposit=<tol>`.

12. **Negative evidence** (each new or corrected gate can fail).
    - `tests/hks/Allrun` with the binary of the M5 round-1 sources
      (`FOAM_USER_APPBIN` pointed at a copy of it): 9 of 38 criterion
      lines FAIL, exactly those this round changes in the solver: M5.1a
      (closure 3.4e-13 > 1e-13), M5.1f (2.7e-12), M5.2g ("Inventory sample
      boat was changed ...: volume 8.117762652854522e-09 (written with
      8.117762652854516e-09)" on 4 ranks), M5.2h (the second restart stops
      with "the ledger counts 1 sample removals ..."), M5.2i (2 warnings
      about uniform/ at 0.2, reservoirs from `Cs_`), M5.3b (the excerpt's
      solver-defect digits), the moved-selection refusal (it did stop,
      but only because V_Z differed in its last bit, 8.000000000000005e-09
      vs 8e-09: the defect of item 3, not a named change), the
      removeTime refusal (the run did not stop) and M5.9 (no warning).
      M5.1c, M5.1d, M5.1e, M5.2j and M5.3c pass with that binary: their
      solver behaviour is unchanged.
    - A copy of `tests/hks/Allrun` with impossible limits (closure < 0 for
      M5.1a and M5.1c, 1e-17 for Y in M5.1d, for the closed form in
      M5.1e and for the closure in M5.1f) with the new binary: exactly
      these 5 criteria FAIL, the other 33 lines PASS; with the round-1
      script, M5.1d and M5.1e passed whatever their checkers returned.

---

## 22. Amendments from the M5 review, round 2 (signed off by the researcher, 2026-10-02)

(Items 1, 4 and 9 are corrected or extended by the M5 review, round 3,
section 23; item 9 by round 4, section 24, item 1: the projections of M5.2b
came from guard re-solves of a re-assembled matrix, not from a coupled pair
that need not land on its bound, and they are gone.)

The M5 verification of round 2 passed.  The review found one blocking and
one major defect, both in the regime choice of inventory samples, two
minor points and three nits.  Each was assessed first; all were real and
are fixed and tested (`tests/Alltest M5`: new criteria M5.2k, M5.2l, M5.2m
and M5.3d, M5.0 extended, M5.1b sharpened; M3.2c now runs in `Alltest M5`
as well).

1. **The regimes of a block of sample cells never settled (blocking).**
   - Defect: the E4 predictor solves every exchange cell with its
     neighbours lagged by one corrector (H of the transport-only matrix).
     In the cells of an inventory sample with a small reservoir (run/013
     at 1e-7 mol, Cs = 5.68 kg/m3; the test channel at 1.7e-8 mol, Cs =
     0.98 kg/m3), every cell switched at once: all CAPPED (empty
     neighbours: each cell would lose more than its reservoir), which
     dumps the reservoir into supersaturated gas; then all NONE (the
     neighbours are supersaturated), which drains the gas; and so on at
     every corrector.  The guard only switched IMPLICIT elements, never
     back, so the parity of the last corrector decided the step: run/013
     at nCorr 10 never evaporated the boat (and converged in no step), at
     nCorr 9 or 11 it dumped it in the first step; the channel lost 70 %
     of n0 into supersaturated gas in one step at nCorr 2 and kept half of
     it for 30 steps at nCorr 3.  Neither state obeys the bounded law, and
     the step was accepted with the whole inventory as pending switches
     (section 7, "no acceptance with pending switches").
   - Fix, two rules (`phaseChangeKinetics.H`, `interfaceExchange.C`):
     - the predictor starts from the regimes of the previous corrector (of
       the previous step at the first corrector) and never switches an
       element directly between its two bounds: CAPPED <-> NONE at the
       lower bound (lower = 0, a sample element) with upper > 0 becomes
       IMPLICIT (`acrossBounds`).  The coupled solve of the IMPLICIT block
       resolves the exchange of all its cells at once, and the next
       corrector starts from neighbours that belong to it.  A single cell
       never needs the direct switch (its local solution does not change
       once the neighbours are right); an element that exhausted its
       condensate (upper = 0 in the next step) turns NONE as before, the
       fresh start (UNKNOWN) runs the predictor as before, and a wall
       element is never damped: its NONE is the upper bound 0 of a bare
       wall, on the same side as CAPPED (a first version damped a bare
       wall that had a deposit again, the CAPPED box of M5.1b restarted
       with a deposit through `mw_`; the closure of that step left 0);
     - the guard checks the bounded law at the solved Y for every element
       (the realised flow against the admissible flow min(max(a - g Y,
       lower), upper), `admissibleFlow`) and switches in both directions:
       a violating IMPLICIT element to CAPPED or NONE as before, a CAPPED
       or NONE element whose own missed mass is significant (above
       guardTolerance x inventory) back to IMPLICIT.  An element whose
       implicit flow sits at its bound to round-off (a bare wall at Y_eq,
       a noise-level Y ahead of a gas front) is not switched back: it
       would return with a bound violation of the same size.  The guard
       acts, as before, when the total violation exceeds guardTolerance x
       inventory; for IMPLICIT elements the violation mass is the old one
       bit for bit, so a step whose other elements are complementary
       guards exactly as before.
   - Why not the reviewer's "switch to the regime admissible at the
     solved Y" (CAPPED directly to NONE and back): with lambda dt of 1e3 to
     1e4 it alternates in a single stiff cell.  With the transport-only
     diagonal d, the IMPLICIT solution Y* with flow q* in (0, upper), and
     the CAPPED solution Y_C, d (Y_C - Y*) = upper - q*, so q(Y_C) = q* -
     lambda dt (upper - q*): negative unless q* lies within upper/(lambda
     dt) of the cap.  The admissible regime at Y_C is then NONE, whose
     solution has q > upper, and so on; the implicit solution between the
     bounds is never reached.  Only the adjacent switch (to IMPLICIT) is
     safe, and one pass of it reaches the exact solution of such a cell.
   - Standalone evidence (M5.0, `testPhaseChange.C`): a row of 10 stiff
     sample cells (g = 1000, Y_eq = 0.04, upper 0.2, diffusion 5 to each
     neighbour and to empty gas at both ends), run as the corrector loop
     runs (predictor with the neighbours of the previous corrector, exact
     tridiagonal solve): without the previous regimes the regimes
     alternate at every one of 20 correctors and the final state is not
     complementary; with them they settle after corrector 2, CAPPED at
     the two ends and IMPLICIT inside, complementary.
   - Results (new criteria, item 10):
     - the channel at 1.7e-8 mol, tight solves, nCorr 2, 3, 10 and 11: the
       same regimes in every one of 30 steps (sample implicit/capped/none
       28/12/0, 16/10/14, 12/6/22, 10/2/28, ...), inventories equal to
       1.3e-4 of n0 M between nCorr 2 and 3 and 2.7e-6 between 10 and 11
       (the SuperBee correctors are not converged at the outer tolerance
       1e-12; the round-2 build: 0.65 of n0 M and 0.5), the bounded law of
       every sample cell recomputed from the fields to 7e-17 of n0 M in
       total, closure 4e-16 to 8e-16 (M5.2k);
     - the same with one corrector per step: with nGuard 0 the first step
       reports 34 significant complementarity violations, the 34 CAPPED
       sample cells in supersaturated gas that the fields show (the
       round-2 check reported 0); with nGuard 3 the guard alone repairs
       them (2 passes in the first step, 9 steps with passes) and the
       regimes are those of nCorr 3 in every step (M5.2l);
     - run/013 at 1e-7 mol, its own settings: every step converges with 6
       to 8 correctors, identically for nCorr 9, 10 and 11 (balance files
       identical); the boat evaporates over 7 steps (sample 2.8e-8, 1.9e-8,
       1.4e-8, ... 8.0e-10 kg, then 0), the law holds in its 144 cells to
       3e-18 of n0 M, closure 4.3e-16, solverDefect 5.2e-13 (M5.3d; the
       round-2 build: nCorr 10 converged in no step and never evaporated
       the boat, nCorr 9 dumped it in the first step);
     - every other M5 case is unchanged bit for bit: rerun with the
       binary of the round-2 sources, the written files and the balance
       files of the boxes (M5.1a-f), the channel with nCorr 3 (M5.2a) and
       with the recommended solves (M5.2j, nCorr 3), the pipe (M5.3a) and
       run/013 (M5.3b/c) are identical, and run/013 reproduces its
       out_excerpt.  Only the two channels with nCorr 0 change (M5.2b,
       M5.2j with nCorr 0), where the guard now switches back the missed
       deposition that the round-2 guard left (item 9).

2. **The complementarity check was blind to wrong sample regimes
   (major).**
   - Defect: it counted non-IMPLICIT elements whose implicit flow lay
     within [lower, upper], so a CAPPED sample element in supersaturated
     gas (implicit flow < 0) and a NONE sample element in undersaturated
     gas (implicit flow > upper) were not counted: the log said
     "complementarity 0" exactly when the whole inventory was misplaced.
   - Fix: a non-IMPLICIT element is a violation when its realised flow
     differs from the admissible flow of its implicit flow, with the
     missed mass |admissible - realised| dt.  Where the old check counted
     an element, the value is the same; a NONE bare wall below Y_eq
     (upper = 0) and an empty sample still give 0.  The guard uses
     the same criterion, so a step accepted by the guard below its
     tolerance has no significant violation.
   - Independent check (`tests/hks/checks.py complementarity`): the law
     recomputed from the written fields (Y, exch_, Cs_ of the previous
     step, T, rho, the cell volumes and centres of `postProcess
     writeCellVolumes/writeCellCentres`, the table read and interpolated
     in Python), for every sample cell selected as boxToCell or
     cylinderToCell select it.

3. **M5.1b gated on "closure <= 1e-15" (nit).**  Acceptance 1 says
   "closure 0"; the gate is now 0 (the measured value).

4. **An ascii restart after decomposePar or renumberMesh was refused
   (minor).**  With `writeFormat ascii` these utilities store the mesh
   points at writePrecision, which moves V_Z and the centre of the same
   cells (7e-10 relative for the boat at 6 digits), and the restart
   record compared them to 1e-12: "Inventory sample boat was changed".
   The tolerance now follows the format of the mesh points
   (`constant/polyMesh/points` of the case or of every processor, read
   with `writtenBinary`, collective): 1e-12 when they are binary,
   10^(1 - writePrecision) (the bound of a rounded restart, section 12,
   NB1) when they are ascii.  The pair, the amount and the number of
   cells stay exact.  The refusal names the tolerance and the format and
   warns about selections whose boundary passes through cell centres.
   New M5.2m.

5. **run/013's boat lay on a ring of cell centres (minor).**  The radius
   5e-4 was 2e-12 (relative) inside a ring of 18 centres at
   5.000000000000005e-4 m; a mesh rewritten at lower precision (ascii
   decomposePar) selected 146 cells, so a parallel fresh start used
   another sample than the serial run and a restart stopped ("cells 146
   (written with 144)").  run/013 now has `radius 4.8e-4`: the same 144
   cells, 4 % from both rings (0.458 and 0.5 mm).  Its readme lists the
   rings; section 11 now says that "r <= 0.5 mm, 144 cells" is the strict
   test and must be pinned for M6 (162 cells with the inclusive test).
   Only the input changes: the fields and logs of run/013 are the same
   (M5.3b reproduces its `out_excerpt`).

6. **A Cs_ that differs from uniform/phaseChangeReservoirs was ignored
   silently (nit).**  It is now compared in the sample cells of the pair,
   as `mw_` is with the deposits (tolerance 0 for a binary field without
   a repeated value, otherwise its rounding), reported with a warning and
   rewritten from the list (section 20, item 5: "enters the restart
   verdict like mw_").

7. **Style (nits).**  The header block "State, restart and output" of
   `interfaceExchange.H` names the three per-rank lists and their
   fingerprints; `pairData::removedStep` (never read) is removed; the
   85-column line of `phaseChangeLedger.H` is re-wrapped; the second
   comment block of `interfaceExchange.H` (13.7 KB, 2.3 KB below the
   wmkdepend limit of section 18, item 13) is split in two (10.5 and 3.9
   KB).  The block of step algorithm and the readme describe the new
   predictor, guard and complementarity rules.

8. **Carry-over M3.2c in `Alltest M5` (non-blocking point of the
   verification).**  Section 20, item 15 lists the uneven 3-rank split as
   part of M5, but only `Alltest M3` ran it.  The pipe criteria of
   `tests/profiles/Allrun` moved into `tests/profiles/pipeSection`
   (sourced), and `tests/profiles/Allrun M3.2c` runs only the serial pipe
   and the 3-rank split; `Alltest M5` runs it unless M3 is run too.

9. **Acceptance 2 with nCorr 0 (M5.2b; proposed rewording, pending
   sign-off).**  The literal "0 projections" of section 6 held in round 2
   because the monotone guard never switched back: its log showed
   significant complementarity violations (missed deposition of bare wall
   elements, NONE in supersaturated gas) in 56 of 350 steps, which no gate
   checked.  The two-sided guard repairs them (57 steps with guard
   passes, 0 significant complementarity violations); in 1 step the two
   repaired elements end 1.0e-18 kg past their bound after the re-solve
   (0.46 of guardTolerance x inventory; the re-solve re-evaluates the
   SuperBee limiter, so a coupled pair need not land exactly on its
   bound), which E5 leaves to the booked projection.  M5.2b therefore
   gates "guard passes > 0, 0 significant projections and 0 significant
   complementarity violations", the rewording proposed for acceptances 2
   and 3 in section 20, item 11 and section 14, item 1; M5.2a (nCorr 3)
   still has 0 projections and now also gates 0 significant
   complementarity violations.

10. **Tests.**  New: M5.2k (item 1, the small reservoir in the channel),
    M5.2l (items 1 and 2, the guard alone and the complementarity count
    against the fields), M5.2m (item 4, with the radius of item 5), M5.3d
    (item 1, run/013 at 1e-7 mol); M3.2c in `Alltest M5` (item 8).
    Changed: M5.0 (the admissible flow, the predictor with the previous
    regimes, the row of sample cells), M5.1b (closure 0), M5.2a
    (complementarity), M5.2b (item 9).  `tests/hks/checks.py` has
    `complementarity` (the law recomputed from the fields), `sameRegimes`
    and `inventories`.

11. **Negative evidence.**  The cases of the new criteria rerun with the
    binary of the round-2 sources, checked as the criteria check them:
    - M5.2k: the regimes differ between nCorr 2, 3, 10 and 11 (step 1:
      sample implicit/capped/none 4/36/0 with nCorr 2); the inventories
      of nCorr 2 and 3 differ by 0.65 of n0 M, those of 10 and 11 by
      0.61; the fields show 28, 682, 24 and 678 sample elements missing
      the law, while every log reports 0 significant complementarity
      violations;
    - M5.2l: the first step logs 'complementarity 0 (0 kg, 0
      significant)' where the fields show 34 CAPPED sample cells in
      supersaturated gas; with nGuard 3 the same 34 remain (the guard
      never switched back);
    - M5.3d: nCorr 9 dumps the boat in the first step (sample 0 from then
      on, 144 cells in supersaturated gas), nCorr 10 converges in no step
      and never evaporates (sample 4.61e-8 kg in all 12 steps, 1728
      element-steps missing the law); the balance files differ;
    - M5.2m: 'Inventory sample boat was changed at the restart from time
      0.1: volume 8.117762658585557e-09 m3 (written with
      8.117762652854519e-09; tolerance 1e-12 relative)';
    - item 5, with the new binary: a fresh start of run/013 on the
      ascii-decomposed mesh selects 146 cells (V_Z 8.23051e-09 m3) with
      radius 5e-4 against 144 in serial, and 144 in both with 4.8e-4.

The other points of the verification (the decision of section 21, item 1
on the stiff round-off imbalance, the rewording of acceptance 2 in
section 20, item 11, the 4 ms steps of M5.3a, and M0.4 on Merlin7) are
unchanged and still pending.

M5 status after round 2: `tests/Alltest M5` gives 53 PASS, 0 FAIL (7 min;
M0 gate identical to main, M3.2c included); `tests/Alltest all` gives 177
PASS, 0 FAIL, 1 NOTRUN (M0.4, the v2606 compile on Merlin7); 0 warnings.
The sources compile and link with 0 warnings against the local minimal
v2606 build (Debug), and that binary reproduces the small-reservoir
channel of M5.2k (nCorr 2) bit for bit (balance file identical).  Pending
the researcher's sign-off of this section, with sections 20 and 21.

---

## 23. Amendments from the M5 review, round 3 (signed off by the researcher, 2026-10-02)

(Items 1, 2 and 4 are corrected or extended by the M5 review, round 4,
section 24: the per-rank lists lie in uniform/phaseChange/ and a container
is read by the ranks of its IO group; the geometry bounds are floored at
min(6, writePrecision) digits and limited by a one-cell change of the
selection; acceptance 2 holds as written, so the decision of item 4 is
moot.)

The M5 verification of round 3 passed (53 PASS), with acceptance 2 at
nCorr 0 passing only under the pending rewording of section 22, item 9.
The review found two major defects (restarts with the collated file handler
onto another number of ranks; the geometry tolerance of the inventory
record), two minor points (the guard acting on single misses only; an
unsigned rewording gated as PASS) and a nit.  Each was assessed first; all
were real and are fixed and tested (`tests/Alltest M5`: new criteria
M5.2n and M5.2o, M5.2b with a second result line, M5.0, M5.2a, M5.2k,
M5.2l and M5.3a extended).

1. **Collated restarts onto another decomposition (major).**
   - Defect: the ledger (`uniform/phaseChangeBalance`), its layout and the
     three per-rank lists were per-rank objects.  With `-fileHandler
     collated` a per-rank object is one `decomposedBlockData` container per
     time (`processors<N>/<t>/uniform/<name>`, one block per rank), and
     `reconstructPar` and `decomposePar` copy `uniform/` verbatim.  A
     serial restart after `reconstructPar -fileHandler collated` stopped
     with 'could not detect processor number from objectPath', and a run
     on another number of ranks after `decomposePar -fileHandler collated`
     with 'incorrect first token ... error' (a missing block), both before
     the documented `mw_`/`Cs_` fallback (section 16, item 1; section 20,
     item 5) could apply.  Reproduced (the reviewer's expE and expH, and
     M5.2n with the binary of the round-3 sources, item 8).
   - Fix (`binaryIOList.H`, `phaseChangeLedger.H`, `interfaceExchange.C`):
     - the ledger, which is the same on every rank, is a global object
       (`binaryGlobalIOList`, a `GlobalIOList` written in BINARY), and so is
       the layout (`IOdictionary` instead of `localIOdictionary`): the
       master reads and broadcasts them; uncollated every rank still writes
       its identical copy, collated the master writes one plain file, which
       any number of ranks reads, like OpenFOAM's `uniform/time`;
     - the only per-rank content of the layout, the record of rounded
       'uniform' entries (`uniformValues`), is merged over the ranks in
       rank order; every use of it was already reduced with 'or' over the
       ranks, so the verdicts are those of the per-rank records;
     - the per-rank lists (deposits, regimes, reservoirs) are read only
       when the file handler can read them for this run
       (`perRankFile()`: a plain list, or a container of as many blocks as
       ranks; the fingerprint then decides as before).  A container of
       another number of ranks, or a file that some ranks lack, counts as
       written for other elements: the warning names it ('... is a
       container of the collated file handler written on 2 ranks, which a
       serial run cannot read (reconstructPar and decomposePar copy
       uniform/ verbatim)'), and the deposits and reservoirs come from
       `mw_` and `Cs_`, the regimes start unknown.
   - Results (M5.2n): the inventory channel on 2 ranks, collated, to 0.2
     (in `processors2/0.2/uniform/` the ledger and the layout are plain
     files, the three lists containers); reconstructed and continued in
     serial, and re-decomposed with the collated handler onto 4 ranks, both
     to 0.3 across the removal: 3 warnings each, 'relative 0; ...
     reservoirs: Cs_PbI2_s written in binary; exact', one removal, closure
     3.8e-16 in every step.  A collated restart on the same ranks stays
     cmp-identical to the continuous run (M2.5l, 4 ranks; also checked by
     hand: every file of `processors4/0.1` of a restart from 0.05 but
     `uniform/time`, whose time value differs in its last bits).
     Uncollated and serial files are unchanged but for the new keys of the
     sample record (item 2): the headers of the ledger (`class
     scalarList`) and the layout (`class dictionary`) are the same.
   - Not covered: a collated state written by an earlier build of this
     branch (the ledger as a container) is not readable by the new build;
     no such state exists outside the test work area.

2. **The geometry tolerance of the inventory record (major).**
   - Defect: the tolerance of the recorded V_Z and centre of an inventory
     sample (section 22, item 4) followed the format of the mesh points of
     the reading run only (1e-12 for binary points, 10^(1 -
     writePrecision) for ascii).  (a) A run on ascii points (decomposePar
     with writeFormat ascii) records a V_Z that the restart on the
     original binary mesh, after reconstructPar, compared to 1e-12: the
     boat of run/013 was refused ('volume 8.117762652854519e-09 m3 (written
     with 8.117762658585557e-09; relative tolerance 1e-12; mesh points in
     binary)', 7.1e-10).  (b) The V_Z error of rounded points is not a
     relative rounding: it grows with |x| over the cell size, and
     decomposePar writes points with at least 10 digits whatever
     writePrecision says; the same boat moved to x = 0.60-0.61 m was
     refused at writePrecision 10 (5.8e-9 > 1e-9) and 12 (4.2e-11 >
     1e-11).  Reproduced (the reviewer's expI and expG; M5.2o with the
     round-3 binary, item 8).
   - Fix (`sampleGeometryTolerance()`, `recordSamples()`): every run
     records with the sample the error bounds of V_Z and the centre that
     the rounding of its own mesh points implies (`volumeTolerance`,
     `centreTolerance` [m3, m] as exact decimals, and `meshPoints`,
     'binary' or 'ascii, <p> significant digits'): 0 for binary points; for
     ascii points, p is the most significant digits any coordinate of the
     sample's cells needs (`significantDigits()`; a value written with p
     digits never needs more, so p is a lower bound of the digits written,
     whatever utility and precision wrote them), every point is off by at
     most eta = sqrt(3) 0.5 10^(e + 1 - p), e the decimal exponent of the
     largest coordinate, and to first order V_Z changes by at most eta
     times the surface of the sample and the centre by at most eta times
     that surface times the largest distance of a point from the centre,
     over V_Z; the surface is bounded by the sum of the face areas of the
     sample cells, and both bounds are doubled (the face centres of the
     decomposition, second order).  A restart compares with 1e-12 (the
     summation order) plus the bounds of BOTH runs, the writer's from the
     record (0 in records of an earlier build).  The refusal prints the
     difference, the allowed value and both point formats.
   - Results (M5.2o): (a) the state of M5.2m (4 ranks, ascii points)
     reconstructed and continued for one step on the binary mesh: accepted,
     144 cells; the ascii run recorded 'ascii, 10 significant digits',
     volumeTolerance 3.42e-15 m3 (4.2e-7 of V_Z) against a difference of
     5.7e-18 m3 (7.1e-10).  (b) The boat at x = 0.60-0.61 m (144 cells),
     one serial step, decomposed with writeFormat ascii: at writePrecision
     10 the bound is 3.42e-14 m3 (4.2e-6) against 4.7e-17 m3 (5.8e-9), at
     12 it is 3.42e-16 m3 (4.2e-8) against 3.4e-19 m3 (4.2e-11); both
     accepted with 144 cells ('ascii, 10' and 'ascii, 12 significant
     digits'); at writePrecision 15 (by hand) 3.4e-19 m3 against 3.5e-23
     m3.  The centre bounds are 2.1e-8 m (10 digits) and 2.1e-10 m (12
     digits).
   - What the check still detects: pair, amount and the number of cells
     are exact; with binary points the comparison stays at 1e-12 (the
     moved selection of M5.5 is refused as before); a selection changed by
     one cell of the boat moves its centre by about a cell size over the
     144 cells (about 3e-6 m), 100 times the bound of 10-digit points.  The
     bounds are 700 to 1000 times the differences observed (a worst case
     without cancellation); a looser bound only weakens the named
     diagnostic, while the inventory re-read from `Cs_` still enters the
     restart booking and its bound.

3. **The guard acted on single misses only (minor).**
   - Defect: the guard was triggered when the TOTAL violation (IMPLICIT
     bound violations plus the missed mass of CAPPED and NONE elements)
     exceeded guardTolerance x inventory, but switched a CAPPED or NONE
     element back only when its OWN miss exceeded it.  Several smaller
     misses were accepted without a pass: in M5.2b (nCorr 0) the step
     ending at 0.167 s ended with 'guard passes 0 (0 elements); ...
     complementarity 2 (3.18453e-18 kg, 0 significant)', two bare wall
     elements just above Y_eq, 1.43 times the tolerance of 2.23e-18 kg
     together (recomputed from the balance file by `checks.py totalLeft`).
     No mass is lost (the ledger closes); the law is violated above the
     tolerance and the log said '0 significant'.
   - Fix: the guard acts on the total (`switchBackThreshold()`,
     `phaseChangeKinetics.H`): the CAPPED and NONE elements with a miss go
     back to IMPLICIT with the largest misses first, until the misses left
     add up to at most the tolerance (every miss above the tolerance is
     always switched, as before; where the misses left by the old rule add
     up to at most the tolerance, the same elements are switched).  A pass
     is therefore never triggered without a switch, and with passes left
     the step ends within the tolerance.  The smallest misses, those of
     elements at their bound to round-off, go last.
   - Reporting: a projected or missed mass of a step above the tolerance is
     named in the exchange line ('complementarity 34 (6.66158e-09 kg >
     tolerance 7.83715e-21 kg, 34 significant)') and the end-of-run summary
     counts such steps ('46 complementarity violations (46 significant; 6
     steps above the tolerance in total)', only when there are any, so the
     logs and `out_excerpt`s of runs without them are unchanged).
   - Results: M5.2b now switches the larger of the two elements back in
     the step ending at 0.167 s (58 steps with guard passes, was 57); the
     re-solve leaves it 9.2e-21 kg past its bound, projected (below the
     tolerance); complementarity in 5 steps (was 6), the largest total 0.28
     of the tolerance.  The balance file of M5.2b is identical up to that
     step; every other M5 case is unchanged bit for bit (item 8).
   - Tests: M5.0 (`switchBackThreshold`: the two misses of M5.2b, misses
     within the budget, a mix with a round-off miss, a single miss above
     the budget); M5.2b gates 'no step whose complementarity mass exceeds
     guardTolerance x inventory' (`totalLeft`, from the balance file, and
     no flag in the log), as do M5.2a, M5.2k, M5.2l (nGuard 3) and M5.3a;
     M5.2l checks that with nGuard 0 the exchange line and the summary name
     the total above the tolerance.

4. **Acceptance 2 as written (minor).**  Section 6 asks for 'guard passes
   > 0 and 0 projections' at nCorr 0; M5.2b gated the rewording '0
   significant projections' (section 20, item 11; section 22, item 9),
   which awaits the researcher's sign-off, so `Alltest M5` showed PASS on
   an unsigned criterion.  M5.2b now records two lines: the reworded gate
   (PASS/FAIL), and acceptance 2 literally, PASS when it holds and
   otherwise NOTRUN with the projections: 'projections in 2 steps
   (1.0147e-18 9.23326e-21 kg), each below guardTolerance x inventory; it
   passes only as reworded ..., which awaits the researcher's sign-off'.
   `Alltest M5` therefore lists it among the NOTRUN criteria until the
   decision.  The second projection is the element that the guard of item
   3 switches back.  Decision needed: sign off the rewording (then the
   literal line is dropped), or keep the literal criterion (then M5.2b
   fails: the projections are the round-off that E5.4 leaves to the
   projection by design).  With nCorr 3 (M5.2a) the literal criterion
   holds (0 projections).

5. **Style (nit).**  The 99-column comment line of `phaseChangeKinetics.H`
   is re-wrapped; the only lines of the sources above 80 columns are the
   two lines of 81 columns of `rhoFixedFlowFoam.C` inherited from main.

6. **Points of the verification.**  The decision on the stiff round-off
   imbalance (section 21, item 1) is still pending and not gated; M0.4
   (v2606 on Merlin7) remains NOTRUN (no remote access); the local compile
   against the minimal v2606 build was not repeated this round (the
   command was refused in this session).

7. **Tests.**  New: M5.2n (item 1), M5.2o (item 2).  Changed: M5.0,
   M5.2a, M5.2b (two lines; item 4), M5.2k, M5.2l, M5.3a (item 3).
   `checks.py`: `logCounts` parses the new flag and counts it
   (`projectedAbove`, `complementarityAbove`); new `totalLeft`.
   `tests/hks/Allrun`: `decomposeWith` passes decomposePar options.

8. **Negative evidence.**  `tests/hks/Allrun` with the binary of the
   round-3 sources (`FOAM_USER_APPBIN` pointed at a copy of it): 4
   criterion lines FAIL, exactly those this round changes in the solver:
   - M5.2b: 'steps whose complementarity mass exceeds guardTolerance x
     inventory 1 (largest ratio 1.43 at t = 0.16700000000000012), flagged
     in the log 0';
   - M5.2l: with nGuard 0 the exchange line names no total above the
     tolerance and the summary counts none;
   - M5.2n: the 2-rank collated run writes the ledger and the layout as
     containers, and the serial restart after reconstructPar stops with
     'could not detect processor number from objectPath:
     ".../0.2/uniform/phaseChangeBalance"'; the restart on 4 ranks after
     `decomposePar -fileHandler collated` stops with 'incorrect first
     token, expected <int> or '(', found on line 28: error' (reproduced by
     hand: the Allrun chained both branches, so the second was not reached;
     they are now independent, and every case of M5.2n is removed first);
   - M5.2o: all three restarts refused, '(a) ... volume
     8.117762652854519e-09 m3 (written with 8.117762658585557e-09; relative
     tolerance 1e-12; mesh points in binary)', '(b) writePrecision 10: ...
     relative tolerance 1e-09 ...', '(b) writePrecision 12: ... relative
     tolerance 1e-11 ...'.
   The other 40 lines PASS; the literal line of M5.2b is NOTRUN with 1
   projection step.  Of the 46 balance files both binaries wrote (runs the
   round-3 binary could finish), 45 are identical and only that of M5.2b
   differs, from the step ending at 0.167 s (item 3).

M5 status after round 3: `tests/Alltest -fresh M5` gives 55 PASS, 0 FAIL,
1 NOTRUN (the literal line of M5.2b, item 4; 8 min; the M0 gate identical
to main, M3.2c included; the tracked and untracked files of the repository
are unchanged by the run); `tests/Alltest all` gives 179 PASS, 0 FAIL, 2
NOTRUN (M0.4, the v2606 compile on Merlin7, and the literal line of M5.2b;
19 min); 0 warnings.  Pending the researcher's sign-off of this section,
with sections 20, 21 and 22, and the decision of item 4.

---

## 24. Amendments from the M5 review, round 4 (signed off by the researcher, 2026-10-02)

The M5 verification of round 4 failed on one line only: acceptance 2 at
nCorr 0 as written in section 6 ('guard passes > 0 and 0 projections';
M5.2b had projections in 2 steps, 1.0e-18 and 9.2e-21 kg), which passed only
under the rewording proposed in section 20, item 11 and section 22, item 9.
The review found one blocking defect (the cause of those projections), two
major defects (the per-rank restart lists under the file handlers that read
on the master; the geometry bounds of the inventory record on a mesh of
round coordinates) and two nits.  Each was assessed first; all were real
and are fixed and tested (`tests/Alltest M5`: new criteria M5.2p and M5.2q,
three more M5.5 lines, M5.2b a single line with the literal acceptance 2).

1. **The guard re-solved another matrix (blocking).**
   - Defect: every guard re-solve re-assembled `ddt + div(phi, Y, SuperBee)
     - laplacian` at the latest Y, so the lagged SuperBee limiter (and the
     explicit non-orthogonal correction) changed between the passes.  A
     switch decided at the solution of one matrix was applied to another.
     The monotonicity the guard relies on holds for a FIXED linear matrix
     A: switching one violating IMPLICIT element e to its bound gives the
     implicit flow q_N = q* (1 + g (A^-1)_ee), on the same side of the
     bound as q*; with a changing matrix a single element cycled.  Next to
     a colder wall (the reviewer's case: the small-reservoir channel with
     the boat in the second cell row above the WALL at 700 K,
     `interfaceTemperature wall`, 4.25e-7 mol) a bare wall element at Y_eq
     on a moving gas front went IMPLICIT -> NONE -> IMPLICIT -> NONE until
     nGuard was used up, and the step was accepted with a bounded-law miss
     far above guardTolerance x inventory (flagged in the log), depending
     on nCorr; mass stayed conserved (closure <= 3e-16).  With `Gauss
     upwind` (a linear matrix) the same runs had no guard pass at all.  The
     same mechanism left the two projections of M5.2b, which section 22,
     item 9 attributed to 'the re-solve re-evaluates the SuperBee limiter'
     and proposed to accept by rewording the criterion instead of removing
     the cause.
   - Fix (`solvePairedGasSpecies.H`): the transport-only matrix of the
     solve that `recordSolve()` booked is kept (`recordedTransport`), and
     a guard re-solve starts from a copy of it instead of re-assembling;
     only Su and Sp change.  The ledger identity is unaffected: the defect
     and `flux()` are those of the matrix actually solved.  Transport
     convergence comes from nCorr and the outer tolerance, never from the
     guard.  The guard comment (`interfaceExchange.C`), the header of
     `solvePairedGasSpecies.H` and the readme say so.
   - Results (new M5.2p: the reviewer's case, tight solves with nCorr 2,
     3, 10 and 11 for 30 steps, the channel's loose correctors and tight
     Final solve with nCorr 0, 1, 3 and 10 for 100 steps): in all 8 runs 0
     projections and 0 complementarity violations of any size, 16 steps
     with guard passes in total, closure at most 4.1e-16; the bounded law
     recomputed from the written fields (new `checks.py wallLaw` for the
     200 wall elements, `complementarity` for the 10 sample cells) misses
     by at most 9e-18 of the inventory.  The build before the fix fails
     6 of the 8 runs (item 6).  In the reviewer's tight nCorr 3 run the
     fields show the two misses (2.99e-7 of the inventory, 5.9e-14 kg, at
     0.014 and 0.017 s; `wallLaw`), which the fix removes.  The balance
     files of every other M5 case are identical to those of the round-3
     sources (48 of 51); only the three runs whose guard re-solved change:
     M5.2b, M5.2j with nCorr 0 and M5.2l with nGuard 3.
   - Consequences:
     - acceptance 2 holds as written (M5.2b, nCorr 0: guard passes in 61
       steps, 0 projections; complementarity violations in 6 steps, the
       largest total 0.64 of the tolerance, none significant; closure
       3.8e-16), so the proposed rewording of acceptance 2 (section 20,
       item 11; section 22, item 9; section 23, item 4) is withdrawn and
       M5.2b gates the literal criterion as a single line (no NOTRUN).  The
       rewording of acceptance 3 (section 14, item 1) was signed off and is
       unchanged;
     - M5.2j, the same channel with the recommended loose correctors (an
       additional check, section 21, item 10), still has projections in 1
       step (nCorr 3, 1.4e-20 kg) and 4 steps (nCorr 0, at most 7.8e-19
       kg), all below guardTolerance x inventory and all in steps WITHOUT
       a guard pass: E5.4 leaves a violation below the tolerance to the
       projection by design.  M5.2j keeps its gate '0 significant
       projections'; it is not acceptance 2 as written, and its header
       now says so.  (Making it literal would need the guard to act below
       guardTolerance x inventory, at the level of the round-off of the
       stiff cells, section 21, item 1: a decision for the researcher, not
       a defect.);
     - M5.2l re-baselined: with one corrector per step the guard alone
       still leaves 0 misses (fields and log), but the regimes of nCorr 3
       held in every step only because the re-assembling guard acted as
       extra correctors (9 steps with passes).  Now the regimes are those
       of nCorr 3 in the first step, the one the guard repairs (the 34
       CAPPED sample cells of nGuard 0), and in 23 of 30 steps; the sample
       exhausts its last 2 cells one step later (steps 14 and 15) and 2
       wall elements start depositing one step earlier (steps 21-25).  The
       inventories agree with nCorr 3 to 7.2e-3 of n0 M (the round-3 build:
       1.5e-3; nCorr 2 vs 3: 1.3e-4), the transport error of one corrector
       of the lagged SuperBee scheme per step.  M5.2l gates the first
       step's regimes and the inventories to 2e-2 of n0 M (the round-2
       failure it guards against differed by 0.65).

2. **The per-rank lists under the file handlers that read on the master
   (major).**
   - Defect: in v2412 (and v2606) `masterUncollatedFileOperation::filePath`
     resolves every object whose local directory is exactly "uniform" to
     the MASTER's file and broadcasts that path, because it takes uniform/
     for global data; the collated handler derives from it.  With
     `-fileHandler masterUncollated`, and with `collated` and several IO
     ranks (`-ioRanks`, `FOAM_IORANKS`), every rank read rank 0's
     `uniform/phaseChangeDeposits`, `phaseChangeRegimes` and
     `phaseChangeReservoirs` (or the container of rank 0's IO group).  The
     fingerprints refused them, so no state was misattributed, but every
     restart on the same ranks and decomposition discarded the exact
     deposits, reservoirs and regimes (review item B5) and re-read the
     rounded `mw_` and `Cs_`: not cmp-identical (14 files differed at 0.1
     with masterUncollated, 9 with `-ioRanks '(0 2)'`), with a warning that
     gave a false cause.  `perRankFile()` also compared the blocks of a
     container with the size of the whole run, not of its IO group.
   - Fix (`interfaceExchange.C`): the three per-rank lists lie in
     `<time>/uniform/phaseChange/` (the local directory
     "uniform/phaseChange", which the exact comparison of `filePath` does
     not match; OpenFOAM's own `uniform/functionObjects/` is the precedent);
     `decomposePar` and `reconstructPar` copy (or link) the subdirectory
     with `uniform/`, so the fallbacks of sections 16, 20 and 23 are
     unchanged.  `perRankFile()` splits the path of the container
     (`fileOperation::splitProcessorPath`): `processors<N>/` is read by a
     run on N ranks with one block per rank, `processors<N>_<lo>-<hi>/` by
     the ranks lo..hi of a run on N ranks whose IO group has as many ranks
     as the container has blocks.  The ledger and the layout stay global
     objects in `uniform/`.  The warnings, the FatalErrors and the
     documentation (`binaryIOList.H`, `interfaceExchange.H`, readme) name
     the new paths.
   - Results (new M5.2q): the inventory channel on 4 ranks, a manual split
     in which rank 1 holds a single WALL face (40, 1, 79 and 80 WALL faces),
     to 0.1 with `-fileHandler masterUncollated` and with `-fileHandler
     collated -ioRanks '(0 2)'` (containers `processors4_0-1`,
     `processors4_2-3`), each restarted from 0.05 with the same handler:
     'deposits: uniform/phaseChange/phaseChangeDeposits; reservoirs:
     uniform/phaseChange/phaseChangeReservoirs; exact', no warning, and
     every file at 0.1 cmp-identical to the continuous run (48 and 24
     files; `uniform/time` excluded as before).  By hand also: the
     uncollated and the collated single-IO-rank restarts of the same case
     (exact, 0 files differ), and a state written with one IO rank
     restarted with `-ioRanks '(0 2)'` (exact).  M5.2n (collated onto
     another number of ranks) and every other restart criterion of M2, M3
     and M5 pass with the new location.
   - Not covered: a state written by an earlier build of this branch has
     the lists in `uniform/`; the new build finds them missing, warns, and
     re-reads `mw_` and `Cs_` (exact when they are binary without a
     repeated value).  No such state exists outside the test work area.
     A restart with another IO grouping than the writing run is refused by
     OpenFOAM itself (it looks for `processors4/` of a state written as
     `processors4_0-1/`; checked by hand), so it needs no rule here; the
     other direction (one IO rank written, `-ioRanks` read) reads the
     single container and is exact (checked by hand).

3. **The geometry bounds of the inventory record on a mesh of round
   coordinates (major).**
   - Defect: `sampleGeometryTolerance()` (section 23, item 2) took p = the
     most significant digits any coordinate of the sample's cells needs as
     the precision the points were written with.  On an ascii mesh with
     round coordinates (blockMesh with round spacing; the channel after
     `foamFormatConvert -constant` at writePrecision 6) p is 2 ('0.011'),
     and the bounds exceeded the sample and the domain (volumeTolerance
     1.94e-7 m3 = 24 V_Z, centreTolerance 0.122 m on a channel of 0.1 m).
     A restart with the boat moved by 20 mm was accepted; its reservoir was
     re-read as 0 from `Cs_` in the new cells and 98 % of the inventory
     booked as 'restart', with only a warning: the misattribution the
     sample record exists to refuse (section 17, item 4; section 20, item
     7).
   - Fix (`sampleGeometryTolerance()`, `recordSamples()`): p is at least
     min(6, writePrecision), the fewest digits OpenFOAM's utilities write
     points with (`foamFormatConvert`, `renumberMesh` write at
     writePrecision, default 6; `blockMesh`, `gmshToFoam`, `decomposePar`
     at least 10, `IOstream::minPrecision(10)`), and the allowed
     differences never exceed what a selection changed by one cell causes:
     for the centre a quarter of h_min V_min/V_Z (h_min the thinnest sample
     cell, its volume over its largest face; one cell replaced by another
     of the same volume moves the centre by about a cell distance times
     V_min/V_Z), for the volume half of V_min.  A restart allows 1e-12 plus
     the rounding bounds of both runs, at most these limits; a refusal
     whose difference lies within the rounding bounds but above the limit
     says that the points are too coarse to tell and suggests binary
     points or more digits (the safe side: refused, never misattributed).
   - Results (three new M5.5 lines): the inventory channel with its mesh
     and `0/` converted by `foamFormatConvert -constant` to ascii at
     writePrecision 6 records 'ascii, 6 significant digits' (was 2),
     centreTolerance 1.22e-5 m (was 0.122 m), volumeTolerance 1.94e-11 m3
     (was 1.94e-7 m3); the limits of its boat are 1.25e-6 m (a quarter of
     2e-4 m x (2e-10/8e-9)) and 1e-10 m3.  The restart with the boat moved
     by 20 mm is refused ('distance 0.02 m, allowed 1.25e-06 m ...
     limited to 1.25e-06 m'), and so is a selection with the cell (10, 10)
     replaced by the one below it: its centre moves by 5e-6 m, within the
     rounding bounds of both runs (2.4e-5 m), and the refusal says that the
     difference lies within the rounding of the ascii points.  The
     unchanged restart is accepted (40 cells, 'rounded').  The binary
     refusal of M5.5 is unchanged.
   - What the check still cannot do: on ascii points too coarse for the
     limits (few digits far from the origin) an unchanged restart is
     refused; the message names the remedy.  A changed selection whose
     centre moves by less than the limit (cells of different volumes
     replaced so that their moments cancel) is not detected by the
     geometry; its reservoir re-read from `Cs_` still enters the restart
     booking and its bound.  M5.2m and M5.2o (ascii-decomposed pipe
     states, 10 and 12 digits) pass unchanged.

4. **Nits.**
   - The readme called the settled-regimes rule 'pending sign-off'; section
     14 was signed off on 2026-09-25, and the heading says so now.
   - `tests/hks/log.decomposePar`, a stray artifact of a failed invocation
     of a helper with an empty case path (it held 'decomposePar: command
     not found'), is deleted.  `requireCaseDir` (`tests/functions`)
     refuses a case that is not an absolute path of an existing directory;
     `runSerial`, `runParallel`, `decomposeCase` and `decomposeWith` call
     it, so no helper can run or write in the source tree (checked: an
     empty and a relative path are refused with status 1).

5. **Tests.**  New: M5.2p (item 1), M5.2q (item 2), three M5.5 lines (item
   3: the moved and the one-cell-swapped selection refused, the unchanged
   restart accepted, on the ascii mesh of round coordinates).  Changed:
   M5.2b (one line, acceptance 2 as written), M5.2l (item 1: the first
   step's regimes and the inventories instead of the regimes of every
   step), M5.2j (its header: not acceptance 2 as written); every
   test that names a per-rank list uses `uniform/phaseChange/` (M2, M3.8,
   M5).  `tests/hks/checks.py` has `wallLaw`, the bounded law of the wall
   elements recomputed from the written fields (the face areas from the
   polyMesh files, `interfaceTemperature wall` or cell).

6. **Negative evidence.**  `tests/hks/Allrun` with the binary of the
   round-3 sources (`FOAM_USER_APPBIN` pointed at a copy of it): 37 PASS,
   12 FAIL.  Six lines fail by construction, because they name or remove
   the lists at their new place (M5.1b, M5.1c, M5.2c, M5.2f, M5.2i, M5.2n);
   the other six fail on the defects:
   - M5.2b: 'projected=2 ... guard=58' (1.0e-18 and 9.2e-21 kg);
   - M5.2p: 6 of the 8 runs: tight nCorr 3 'complementarity 2 ... above 2'
     (the fields: 2 wall misses, 2.99e-7 of the inventory), tight nCorr 10
     and loose nCorr 1 and 3 a projection above the tolerance, loose nCorr
     0 and 10 a complementarity miss above it (the fields: 1 wall miss
     each, at 0.017 s); only tight nCorr 2 and 11 pass;
   - M5.2q: both restarts 'deposits: mw_PbI2_s written in binary, a
     'uniform' value on some rank ...; reservoirs: Cs_PbI2_s ...; rounded'
     (relative -1.14e-10), and `Y_PbI2_s` at 0.1 differs from the
     continuous run;
   - M5.5, ascii mesh: the moved and the swapped selection are accepted
     ('the run did not stop'; the moved one books -0.983 of the inventory
     as restart, 'sample 0 ... restart -4.75788e-06' mol), and the record
     says 'ascii, 2 significant digits'.
   M5.2l passes with both binaries (its re-baselined gates hold for the
   old guard too).

7. **v2606.**  The sources compile and link with 0 warnings (6 of 6
   sources) against the local minimal v2606 Debug build, whose
   `masterUncollatedFileOperation::filePath` has the same exact
   comparison; that binary reproduces the balance file of the tight nCorr
   3 run of M5.2p bit for bit.  M0.4 on Merlin7 remains NOTRUN (no remote
   access).

M5 status after round 4: `tests/Alltest -fresh M5` gives 60 PASS, 0 FAIL,
0 NOTRUN (8 min; the M0 gate identical to main, M3.2c included; acceptance
2 now as written); `tests/Alltest all` gives 184 PASS, 0 FAIL, 1 NOTRUN
(M0.4, the v2606 compile on Merlin7; 18 min); 0 warnings.  The tracked and
untracked files of the repository (content and git status) are unchanged
by both runs.  Still pending: the researcher's sign-off of sections 20 to
24, the decision on the stiff round-off imbalance (section 21, item 1; not
gated), and M0.4 on Merlin7.

---

## 25. M4 implementation decisions (signed off by the researcher, 2026-10-03)

(Items 1 to 6 and the M4long status are corrected or extended by the M4
review, round 1, section 26: the layer of item 2 is the exact distance over
the interface faces of all patches and all ranks, patchDataWave and its
consistency check are gone, and mw_ holds the deposits themselves; the
utility checks and writes the carrier where the solver reads it and
reports the flux through the walls; run/014 as distributed is the t_ev =
15 s configuration of Figs. 8-10, which it does not reproduce; items 4a,
4b, 4e, 4h and 6 are corrected; acceptance 5 is gated as the researcher
restated it.  The M4 review, round 2, section 27, corrects items 1, 4d and
5 further: the discrete flow of T-Flows' own carrier, the wall table beyond
0.65 m, and three properties of the T-Flows temperature-model codes.)

M4 (the paper-mode carrier, `layer distance`, the code-to-code case of the
temperature model) is implemented on top of the M2, M3 and M5 code and
tested by `tests/Alltest M4` (`tests/paperMode/Allrun`, criteria M4.1a to
M4.12 without M4.3-M4.5b) and `tests/Alltest M4long`
(`tests/paperMode/longRuns`, criteria M4.3, M4.4, M4.5 and M4.5b; not part
of `all`).  The points below were decided during the implementation; each
sharpens or deviates from a statement above and needs the researcher's
sign-off.  Sources: the paper (`papers/Lobresco2026.pdf`, section,
equation, figure, table) and the T-Flows repository (read-only), cited as
in `run/014-paper-temperature-model/readme.md`: `HKS` = origin/
advanced_scalar_model (35c1f6406e) `Tests/Laminar/Scalar_Transport_Gpu/
LESTO_HKS`, `Thesis` = the same branch, `.../Lesto_Thesis`, `Project` =
origin/multiple_scalars_gpu_francesco `Tests/Laminar/LESTO-project`,
`Transport` = origin/scalar_transport_francesco `.../LESTO`.

1. **`utilities/setFrozenCarrier` (deliverable 1).**
   - A separate application (`utilities/setFrozenCarrier/`, built by
     `Allwmake`, cleaned with the solver by the `wclean` of `Allwmake`'s
     state changes and of the BUILD step; its `Make/<platform>` is
     git-ignored).  No solver logic: the solver reads its fields with
     `writeFrozenFields no` (`findInstance` of `rho`), like any carrier.
   - rho = `density` (1: `HKS:User_Mod/Beginning_Of_Iteration.fpp:41`),
     p = `pressure` (1e5 Pa: the diffusivity at p_ref = 1 bar with
     `pressureUnit bar`, `.fpp:22`), U = 2 U_b (1 - r^2/R^2) along the axis,
     uniform in x (`HKS:User_Mod/Initialize_Variables.f90:23-36`).
   - phi is the exact flux of the analytic velocity: every face is split
     into the fan of triangles from its first point, and on every triangle
     the three-edge-midpoint rule (exact for quadratics) integrates u (axis
     . n) dA.  The triangles of a cell close its surface and the profile is
     divergence-free, so `max |sum_f phi|` per cell is round-off on any
     mesh: 2.2e-18 of the inflow on the pipe, 6.0e-17 on the 5-degree
     wedge (M4.1a/c).  Every triangle is evaluated in a canonical order
     (sorted vertices with the permutation sign, sorted midpoint values),
     and the triangles of a face in a canonical order, so the two sides of
     a processor face are exactly opposite: written with `-parallel` on a
     decomposed mesh the continuity holds to the same 2.2e-18, and the
     fluxes equal those of `decomposePar` bit for bit (M4.1b).  T-Flows
     interpolates the cell velocities to the faces
     (`Initialize_Variables.f90:41-46`), which is not exactly
     divergence-free on a general mesh.
   - Q: `flowRate` [m3/s] or `flowRateMLPerMin` (the paper's rates, used as
     the volumetric flow of the non-expanding parabola), of the whole pipe;
     a wedge carries the fraction theta/(2 pi), theta from its two wedge
     patches or `sectorAngle`.  `rescaleToFacets no` is the default and
     T-Flows' choice (u_max from the analytic R): the discrete inflow is
     then 0.99997 Q on the pipe, its faceted inlet 0.99493 pi R^2, so the
     discrete bulk velocity U_b,d = inflow/area is 1.0051 U_b.  With `yes`
     the inflow is exactly Q (M4.1c on the wedge).  T-Flows' own discrete
     flow differs: its fluxes evaluate the parabola at the vertex-mean
     centres of its cells and faces, 1.0038 Q on this mesh (U_b,d =
     0.09311 m/s against 0.09276 m/s here), so at the same FLOW_RATE
     run/014's carrier is 0.39 % slower than T-Flows' (corrected in
     section 27, item 9).
   - Temperature.  T-Flows solves its energy equation in every outer
     iteration of every step (`HKS:Main_Pro_Frozen.f90:190`) with the wall
     table as a Dirichlet wall condition (`HKS:control:80-83`), 900 K at
     the inlet (`:85-88`), an outflow outlet, density 1, cp 5193
     (`:25`) and the NIST helium conductivity `k_he` of
     `User_Mod/Types.f90`, interpolated linearly in 20 K steps
     (`Beginning_Of_Iteration.fpp:33-36,46-47`); its initial field is 500 K
     (`:72-74`).  T is therefore solved with a wall condition, not imposed
     on the cross-section.  With density 1 the radial relaxation time
     R^2 rho cp/k is 0.1-0.2 s, short against the 4-5 s the gas needs to
     reach the deposition zone, so `mode solved` computes the steady state
     of that equation once: Picard on k(T), the table at the wall-face
     centres, `inletTemperature` 900 K, zero gradient elsewhere; on the
     pipe 26 iterations to 9e-10 K, 305.2 to 1210.0 K, the wall faces
     exactly on the table (M4.1a).  The default convection scheme is
     `Gauss linearUpwind grad(T)`: with `Gauss limitedLinear 1` the Picard
     iteration fell into a 2-cycle (the limiter switches between
     iterates).  `mode profile` (the table on the whole cross-section) and
     `uniform` remain for other uses (the wedge test, M7).
   - The wall table `wall_x_profile_complete.dat` (format TFlows): x in
     column 0 in m although the header says cm, T in column 4.  Its count
     line says 66, the file has 70 rows; T-Flows reads the count (x <= 0.65
     m).  `rows all` (run/014) reads the 70 rows and prints the mismatch;
     `rows count` reads 66 and stops at the faces beyond 0.65 m with
     `outOfRange fatal`, or holds T(0.65 m) there with `hold` (M4.1f).  The
     difference concerns only x > 0.65 m (T < 325 K), where T540_60
     deposits 1.9 % (distance) and 7.1 % (faceCells) of its inventory by 8
     s; T_dep, at 0.48 m, is not affected (corrected in section 27, item
     13: 'far downstream of every deposit'; section 29, item 5: T-Flows
     itself holds no temperature there but 0 K, its zero-initialised
     boundary value, which neither `rows count` with `hold` nor `rows all`
     reproduces).  The conductivity table of run/014 is `k_he` at
     `t_values` of `HKS:User_Mod/Types.f90` (M4.12).
   - `setFrozenCarrier -check [-parallel]` evaluates the continuity of the
     `phi` of the start time without writing (section 15, item 7: on the
     decomposed case).  On the pipe, a decomposition in ascii stays exact
     (the two axial faces of a cell carry the same flux and round alike);
     a processor value written as a 6-digit `uniform` entry is reported
     (M4.1e, the negative evidence).

2. **`layer distance` (deliverable 2; section 2 Notation, section 12 B3
   and items 9-10, section 16 item 1).**
   - Elements: one per cell whose wall distance d_c <= `distance`, d_c from
     `patchDataWave<wallPointData<label>>` over the interface patches
     (`correctWalls`), the label being the global index (`globalIndex` over
     the interface faces of all ranks) of the nearest wall face f(e).
     `A_I` = sum |S_f| of the interface faces and `V_I` = sum V_c of the
     layer cells, both global; `A_e = V_c (A_I/V_I) areaMultiplier`.
     `areaPerVolume` defaults to `layer` and `cell` stops; `wallResistance
     halfCell` stops with "is not available with layer distance" (M4.8).
     `interfaceTemperature wall` takes the temperature of f(e) through the
     map below.
   - A_I and V_I are compensated sums (`globalSum` of `compensatedSum.H`)
     whenever the element areas follow them (`areaPerVolume layer`, hence
     every distance layer).  With plain sums (`gSum`, as before) the 30528
     cells of the pipe gave A_I/V_I with up to some 1e-14 of decomposition
     dependence, and a restart onto another decomposition re-read a wall
     inventory that far from the ledger: 7.3e-15 from serial to 4 ranks,
     2.8e-15 from the uneven split to serial, and 2.0e-14 from 16 ranks to
     serial, which the exact-restart bound 1e-14 refused with a FatalError
     (found when the Fig. 9 run was moved between rank counts).  With the
     compensated sums both test restarts re-read the ledger's inventory
     exactly (relative 0; M4.7e, limit 1e-15).  With `areaPerVolume cell`
     the plain sums are kept: A_I and V_I enter only v_d there, and the last
     bit of v_d sets round-off digits of the logs (run/012's `out_excerpt`
     changed with compensated sums and is reproduced with the plain ones,
     M2.8).
   - (Superseded by section 26, item 1: the per-patch correction chosen
     below lets the last of several interface patches overwrite the
     distance of a corner cell, and every correction sees only the wall
     faces of the cell's own rank; the layer is now the exact distance.)
     OpenFOAM defect.  With `cellDistFuncs::useCombinedWallPatch` (true by
     default in v2412 and v2606) `patchDataWave::correctBoundaryCells`
     stores `nWalls + minFacei`, a label of the combined patch, instead of
     the face's data for the near-wall cells: on the pipe 15263 of the
     15264 wall cells got an unrelated face (the distances are right), and
     on the test channel the first row of one wall got the faces of the
     other wall (its mDep_ landed on the wrong wall).  The layer sets the flag to false for its wave and
     restores it, and checks every element: a nearest face farther from
     its cell than the wall distance allows stops the run.  The two walls
     of the symmetric channel now deposit the same to 1e-12 (M4.7a).
   - Both counts at start-up (section 12, item 10): the patchDataWave layer
     with the distance range of its cells and the nearest cell outside, and
     with the optional `centreDistance { axis; origin; radius; }` the
     T-Flows test R - r <= distance (`HKS:User_Mod/Initialize_Variables.
     f90:114-120`) with the cells only in one rule.  On the pipe: 30528
     cells, A_I = 0.0103918 m2, V_I = 1.21773e-6 m3, A_I/V_I = 8533.68
     1/m; selected wall distances 28.5-88.6 um, the next cell at 154.7 um;
     R - r of the selected cells up to 97.7 um, the next at 163.8 um; the
     same cells by both rules (M4.2).  T-Flows divides the analytic wall
     area 2 pi R L by the layer volume (`:122`), 8544.5 1/m (8544.52 on
     this mesh, 1.0013 times ours; corrected in section 26, item 12): irrelevant for the temperature model with A given (the rate
     per layer cell is A f(T) in both codes), a factor to remember for the
     HKS code-to-code case (M6).
   - A distance on a row of cell centres selects by round-off (3e-4 on the
     test channel, whose second rows lie at exactly 0.3 mm, selected 295
     cells); the tests use 3.5e-4 (400 cells, two rows per wall).  The
     start-up margins show how far the rule is from such a row.
   - Nearest-face map (B3): a `mapDistributeBase` from the global face
     indices of the elements (the local faces first, then the remote ones
     each rank needs).  Pull: face values to the elements
     (`interfaceTemperature wall`).  Push: the masses m''_e A_e of the
     elements summed onto their faces (reverse distribute, `plusEqOp`), so
     on every interface face `mDep_ = sum_{e -> f} m''_e A_e/|S_f|`, and
     `areaIntegrate(mDep_)` is the wall inventory to round-off in serial
     and in parallel, also for elements whose face lies on another rank
     (M4.7a-c: 200 such elements on the channel split by rows, 15264 on the
     pipe with its first cell layer on rank 0; section 29, item 1: in
     round 2 only for maps that run one way, its scheduled reverse
     distribute adding every contribution twice between two ranks that
     map elements to each other's faces).  The internal field of
     `mDep_` is the deposit per wall area that the cell represents,
     m''_e A_e/(V_c A_I/V_I).  `mw_` is a cell field (m''_e of every layer
     cell, 0 on the patches): the fallback of a restart onto another
     decomposition (section 26, item 2: assigned, not computed as m''_e
     A_e/A_e, so that it equals the deposits bit for bit).
   - State (section 16, item 1; section 24): the deposits and regimes of
     the elements are kept in the per-rank lists of `<time>/uniform/
     phaseChange/` as for faceCells; the fingerprint hashes the rule and
     the distance first, then the cell and the exact centre of every
     element.  The faceCells fingerprint and file format are unchanged (the
     M2, M3 and M5 restart tests pass with them).  The layout records
     `layer` and `distance` (exact decimal); a restart with another
     distance or the other rule stops with a FatalIOError naming both
     values (M4.8).  Binary restarts are cmp-identical in serial and on 4
     ranks (M4.7d); a restart onto another decomposition re-reads the
     deposits from `mw_` (exact in binary; used by the long runs, item 7).
   - Profiles and T_dep (decision asked for by the focus): the deposit of a
     layer cell is binned by the axial overlap of **its own cell**, like
     the gas of that cell, not by the extent of its nearest wall face.
     Bins + outside = inventory exactly and locally (1e-16 on the channel,
     M4.9), and the line density is T-Flows' sum over an axial slice.  On a
     mesh extruded along the axis (the pipe) a layer cell and its nearest
     face share their axial extent, so the choice does not change the
     pipe's profiles.  `Twall` and `wallArea` per bin come from the
     interface faces, as for faceCells (equal bit for bit, M4.9), so
     T_dep is the wall temperature at the onset for both rules.
   - Both models: HKS with layer distance runs on the inventory channel of
     M5 (`interfaceTemperature wall` from the nearest faces), closure 1e-13,
     no significant projection or complementarity violation (M4.10).

3. **run/014 (deliverable 3).**  `run/014-paper-temperature-model/`: the
   pipe mesh of run/001 (linked), the carrier of `setFrozenCarrier`
   (`system/setFrozenCarrierDict`, item 1), `model temperature`, the boat
   release, `layer distance` (switch to `faceCells` in one line), the
   T-Flows inlet (zero total flux, section 9 item 15: `codedMixed` with
   valueFraction Pe/(1 + Pe), Danckwerts), SuperBee, dt 1 ms, the outer
   loop of run/012 (nCorr 10, tolerance 1e-5: converged in every step, 2
   correctors in most), loose GAMG correctors and a PBiCGStab final solve
   (0.26 instead of 0.40 s per step on 16 ranks of the idle machine with a
   GAMG final solve), `writeFormat binary`, the profiles of Fig. 9
   (`times (7 12)`, native and 1 cm bins, observable `wall`, fraction
   0.01, fromInlet).  (Section 26, item 6: as distributed now the t_ev =
   15 s configuration of Figs. 8-10, to 30 s, `times (16 30)`, which does
   not reproduce their positions.)  Its readme pins every item of section 11 with its
   source; `out_excerpt` is the 5-step log (M4.11).

4. **The pinning notes, verified.**
   a. **Tdep per case.**  Section 2.3.1: "The deposition temperature Tdep
      is not treated as a free fitting parameter.  Its value is evaluated
      using the semi-empirical formula given in Liu et al. (2025)"; Section
      3.1: 680 K for the reference case (14.6 mg, 100 mL/min); Table 1 lists
      the formula's T_dep^calc per case (775, 712, 666, 788, 671 K), and
      Table 2's temperature-based column (775, 713, 666, 787, 670 K) lies
      within 1 K of it (the interpolated 1 % onset lies within a native
      bin of Tdep, on either side; corrected in section 26, item 5).  So the model's Tdep was set per case: section 11's
      `Tdep = 680 K` applies to Figs. 8-10 only, and the Table 2 runs use
      Table 1's value (section 11 amended).  Acceptance 5 thus tests that
      the simulated onset reproduces the input Tdep within 2 K, through the
      transport, the 1 % threshold and the binning, not an independent
      prediction.
   b. **Injection.**  Section 3.3: "an injection duration of 6 s was
      selected for all temperature-based simulations in this comparative
      analysis"; n0 from Table 1 (moles: 4.84, 4.43, 2.97, 3.15, 3.93e-6
      mol).  Figs. 8-10: n0 = 3.16e-5 mol (Section 3.1) and t_ev = 15 s
      (Fig. 8 at 16 s; Fig. 9 shows t_ev = 15, 30, 45, 60, 75 and 90 s, the
      gas at t_ev + 1 and the solid at 2 t_ev; Fig. 10 at 30 s); the t_ev =
      6 s run of M4.4 is no curve of the paper (corrected in section 26,
      item 7).
   c. **Q per case.**  Figs. 8-10: 1.67e-6 m3/s (T-Flows `FLOW_RATE`,
      `HKS:User_Mod/Types.f90:7`; 100.2 mL/min, the paper's 100 mL/min of
      Section 3.1).  Table 2: 106, 5 and 540 mL/min of Table 1, used as the
      volumetric flow of the non-expanding parabola, as T-Flows uses
      FLOW_RATE (paper mode; no T/273 expansion, `manuscript/
      notes-for-authors.md`, item A).
   d. **Parameters.**  The paper's set A = 220 1/s, k = 8e-3 1/K (Section
      3.1, calibrated against Fig. 10 and "adopted in all subsequent
      simulations"; Fig. 8 caption), with `areaMultiplier 1` (section 15,
      item 3).  The T-Flows temperature-model branches carry other values
      and are not the paper's: `Thesis` A = 150, k = 0.02, Tdep = 638.15 K,
      U_MAX 0.176 m/s, the layer r > 0.9 R and the release into the inlet
      cells over half the run (`Thesis:User_Mod/Types.f90:5-8`,
      `Source.f90:27-29,55-86,114`); `Project` A = 1.35, k = 2.7 with the
      same layer test (`Project:User_Mod/Types.f90:5-8`, `Source.f90:111`).
      Both use U_MAX as the centreline velocity (Q = 95.5 mL/min), release
      n = 3.17e-5 mol and run on 424 axial cells over 0.04-0.69 m (section
      27, item 14).
   e. **The wall table** is a byte copy of `HKS:wall_x_profile_complete.dat`
      (35c1f6406e; M4.12 compares it with `git show`).  It is the only
      measured table available (notes-for-authors.md, F: the tables for 5,
      106 and 540 mL/min are requested), and run/014 uses it for every
      flow rate.  The paper measured a dedicated wall profile for each
      operating condition (Section 2.4) and plots them in Figs. 13-14, so
      the Table 2 runs very likely used them; a different gradient at the
      onset changes the offset of the binned onset from Tdep by up to
      about a bin (1-2 K), the width of the window of acceptance 5
      (corrected in section 26, item 5).
   f. **R = 8.314** (`HKS:User_Mod/Types.f90:18`, section 12 item 9): Eq. 5
      contains no partial pressure, so the temperature model and run/014 do
      not use R; it matters for HKS (p_i = c R T, M6).  The carrier monitor
      of the solver (`max |T - p W/(R rho)|/T = 0.96`) uses OpenFOAM's R and
      is meaningless in paper mode (rho = 1 is not the perfect-gas density).
   g. **The diffusivity** (Eq. 2, Fig. 3) is the T-Flows fit
      (`Beginning_Of_Iteration.fpp:22-27,61-72`) with p_ref = 1 bar
      (`pressureUnit bar`; the paper's Eq. 2 says atm).  The solver uses
      M_He = 4.0026 where T-Flows uses 4.0: D is 0.03 % lower.
   h. **The release** (Section 2.3.1: "a prescribed source term acting over
      a finite time interval ... at the inlet") is the boat cylinder x =
      0.04-0.05 m, radius 4.8e-4 m (144 cells, section 22 item 5), released
      uniformly over t_ev; the `Thesis` code injects into the inlet cells.
      The plateau of Fig. 9 (0.25-0.45 m) does not depend on it; T_dep of
      the transient T106_6 and of the incomplete T5_60 does, and the paper's
      Table 2 runs release at about 0.13 m (T106) and 0.15 m (T5, T540),
      the arrows of Figs. 13-14 (corrected in section 26, item 5).
   i. **Observable `wall`.**  Section 3.3 defines T_dep for the simulated
      profiles on the condensed phase; run/014 and the long runs use
      `observable wall`, not the default `wallPlusGas` (section 19, item 6),
      which adds the gas of the tube as condensing at shutdown, a property
      of the experiment.

5. **Acceptance 3-5: definitions and interpretation.**
   - Acceptance 3 (Fig. 9 identity): on the plateau, upstream of the
     deposition (680 K at x = 0.476 m) and downstream of the release, the
     gas of a release at rate n0/t_ev moves at the bulk velocity of the
     discrete flow, so n' t_ev = n0/U_b,d with U_b,d = inflow/inlet area of
     the carrier used ("of the carrier used").  M4.3 checks every native
     bin with its centre in [0.25, 0.45] m at t_ev + 1 (7 and 16 s), to 1
     %, and prints n0/U_b with the analytic U_b and the paper's 3.38e-4 mol
     s/m: n0/U_b,d = 3.4068e-4 (+0.8 %), n0/U_b = 3.4241e-4 (+1.3 %).  The
     paper's value is that of T-Flows' own carrier at the same FLOW_RATE:
     its face fluxes carry 1.0038 Q on this mesh, n0/U_b,d = 3.3937e-4 mol
     s/m, within 0.4 % of 3.38e-4 (corrected in section 27, item 9; the
     inference 'U_b = 0.0935 m/s, 101.5 mL/min' is withdrawn).
   - Acceptance 4 (the norm): the native wall columns at 2 t_ev of the two
     runs, sum_b |n'_6(b) - n'_15(b)| dx_b <= 0.01 sum_b |n'_15(b)| dx_b
     (L1, relative to the deposit of the t_ev = 15 s run); M4.4 also prints
     the largest bin difference relative to the peak, the deposits and the
     gas left at 2 t_ev.
   - Acceptance 5: T_dep fromInlet at 1 % of the peak of the wall column of
     the native set, at the end of the experiment duration (6, 60, 300 s;
     Table 1).  T106_6 ends at 6 s (the end of the experiment, with most of
     the gas still in the tube).  The 60 s and 300 s runs stop early:
     the deposit of the temperature model is irreversible and only the gas
     inventory G left at a time t can still deposit, so every bin's final
     deposit lies between n'_w(b) and n'_w(b) + G/dx_b.  The earliest
     possible onset takes every bin at its upper bound against the
     threshold of the lowest peak, the latest every bin at its present
     value against the threshold of the highest peak; a run is extended
     (from 10 s by 2 s at 106 mL/min, from 4 s by 1 s at 540 mL/min, the
     bound only after the release has ended) until that interval is
     narrower than 0.05 K (the T106 runs stopped at 14 s with widths below
     2e-4 K, the T540 runs at 8 s with no gas left), and M4.5 requires both
     ends within +-2 K.  A run that has not ended is evaluated from the
     bound at its latest profiles (`checks.py tdepDecided`): met when the
     interval lies within +-2 K, not met when it lies entirely outside,
     NOTRUN otherwise.  T5_60 runs to 60 s.  Both layer rules are run and
     reported per line (10 lines); the reading of "for both layer rules"
     is that each rule must match, and a rule that does not is a FAIL
     line, which records which rule matches.  (Restated by the researcher
     on 2026-10-02: |T_dep - T_paper| <= max(2 K, dT_bin), dT_bin the drop
     of the wall temperature across the native bin of the onset; section
     26, item 5.)
   - dt: 1 ms (Courant 0.12 at 106 mL/min in paper mode), 0.22 ms at 540
     mL/min (section 12, item 6; Courant 0.13), 10 ms at 5 mL/min (Courant
     0.06), with the dt-insensitivity of T_dep demonstrated against 5 ms
     (M4.5b, 0.5 K).
   - M4.3 and M4.4 read the profiles of fixed times (7, 12, 16 and 30 s),
     final once written, so they are evaluated as soon as those exist.

6. **M4.7c: serial and parallel agreement with SuperBee.**  With run/014's
   own loop (nCorr 10, tolerance 1e-5) the pipe (background Y = 0.1, 10
   steps of 4 ms) agrees with serial to 1.2e-11 (gas) and 6.9e-12 (wall) on
   simple (4 1 1) but to 9.8e-9 and 1.3e-9 on hierarchical (2 2 1), above
   the 1e-9 of M2(4): the lagged SuperBee limiter (section 16, item 7).
   The loop cannot be converged to the floor with SuperBee: with GAMG
   correctors to 1e-13 the initial residual of the first step stalls at
   2.3e-6 after 16 correctors (scratch measurement).  M4.7c therefore
   solves every corrector to the round-off floor (PBiCGStab 1e-15, relTol
   0) with a fixed 3 correctors per step (nCorr 2, an outer tolerance never
   reached), so that the iterates depend on the decomposition only through
   round-off (qualified in section 26, item 12: the SuperBee limiter
   amplifies that round-off to 1e-5 relative per cell, while the
   inventories agree to 1e-10), and runs simple (4 1 1), hierarchical (2 2 1) and the uneven
   manual split with the first cell layer on rank 0 against serial at
   1e-9: gas 1.7e-12, 4.3e-11 and 5.8e-12, wall 1.7e-12, 5.3e-12 and
   9.4e-13 (the remaining difference is the floor of the linear solves
   and the lagged non-orthogonal correction, not round-off of the
   exchange); closure 1.7e-16 to 3.3e-16; areaIntegrate of mDep_ = the
   wall inventory to 1.4e-15 to 2.4e-15, also with the 15264 elements of
   the second layer mapped across ranks.  The production loop's 1e-8 on a
   hierarchical split is a property of the outer tolerance, as in section
   17, item 7, and far below what T_dep or the profiles resolve.

7. **Cost of the code-to-code runs and their status.**  A pipe step costs
   0.08-0.15 s on 16 ranks of the idle 40-core workstation, but 2 s and
   more when the machine is oversubscribed (load 50-80 during this round,
   from other users' jobs): the ranks of OpenMPI 1.10 busy-wait, and each of
   the many global reductions of a step waits for the slowest descheduled
   rank.  Measured at load 50 on the 7 s state of the Fig. 9 run: 16 ranks
   2 s, 8 ranks 0.66 s, 4 ranks 0.9 s per step; at load 75 four 4-rank
   runs took 1.8-2.9 s per step each.  Per core, serial runs are the most
   efficient on such a machine (OpenFOAM refuses `-parallel` on one rank,
   so `longRuns` runs them serially in the case).  `longRuns` therefore
   runs `M4LONG_JOBS` runs of `M4LONG_RANKS` ranks at a time (default 2 x
   8 for a quiet machine; this round 7, then 13, serial runs), resumes a
   stopped run from its latest second, and moves it to another number of
   ranks (reconstructPar/decomposePar, an exact restart in binary: item 2;
   the rank count is not part of the run's key).  The fix of item 2
   changed the solver after the runs had started, so their key changed and
   all were restarted from 0 (the first 7 s of the Fig. 9 run, made on 16
   ranks, were discarded).  The runs need at least 154,000 steps (Fig. 9
   12,000 and 30,000; T106 6,000 x 2 and at least 10,000 x 4; T5 6,000 x
   2 and 12,000; T540 at least 18,182 x 2): about 11 h at the 4 steps/s of
   four 4-rank runs, more than the ~8 h of the focus, so the T106 cases of
   both rules and the Fig. 9 run with t_ev = 6 s went first, the others
   after the full test suite had freed the cores (fig9-tev15 on 8 ranks,
   the T5 and T540 runs serially).  When the load fell (to 15-30), runs
   were moved onto 4 ranks (the last onto 8) right after a write (an exact
   restart, relative 0, on another decomposition).  In the end the 13 runs
   took 207,000 steps (the T106 runs to 14 s, the T540 runs to 8 s) and
   about 6 h of wall time.  The first checks of the T540 runs, at 4
   and 5 s, fell inside the 6 s release and were skipped (checks.py asks
   for the end of the release before it bounds T_dep); their t0 of 4 s
   should have been 6 s, a choice that costs checks, not steps.

8. **Tests changed or added.**
   - `tests/Alltest`: `M4` (in `all`) and `M4long` (not in `all`), their
     usage text and expected criteria (M4.8: 7 lines, M4.3: 2, M4.5: 10);
     BUILD requires `setFrozenCarrier` to be compiled (every source of its
     `Make/files`) and present.
   - M2.6f: layer distance is no longer 'not yet implemented'; the closed
     box with `distance 1e-4` (and `areaPerVolume layer`) runs with its 10
     cells, and without a `distance` entry it stops with "Entry 'distance'
     not found" (before: both stopped with 'not yet implemented').
   - `tests/functions`: `makeCase` links `constant/*.dat` (the wall table)
     like `constant/*.csv`; `buildSolver` cleans the utility too;
     `utilitySources`.  `tests/build/Allrun` copies the sources without
     the utility's `Make/<platform>`.
   - M4.7e (new): the restarts of item 2 onto another decomposition.
   - `.gitignore`: `applications/rhoFixedFlowFoam/utilities/*/Make/*/`.

9. **Section 11, amended.**  `Tdep = 680 K` holds for Figs. 8-10; the
   Table 2 runs use Table 1's T_dep^calc (item 4a).  `areaMultiplier 200`
   belongs to the HKS source of T-Flows; the temperature model uses 1
   (section 15, item 3).  The layer of T-Flows (`R - r <= DELTA`) is
   `layer distance 1e-4` (30528 cells, A_I/V_I = 8533.7); the paper's
   text describes one layer (`layer faceCells`, 17681.9); run/014 runs
   both.

M4 status: `tests/Alltest -fresh all` (M0 M1 M2 M3 M4 M5) gives 209 PASS,
0 FAIL, 1 NOTRUN (M0.4, the compile on v2606 under v2412; M4.6 compiles a
copy of the sources against the local v2606 Debug build and passes) in
70 minutes at a machine load of 60-70; 0 warnings; the tracked and
untracked files of the repository (git status and the content of every
file) are unchanged by the run.

M4long status (this round: all 13 runs ended, in about 6 h of wall time
on 13 to 20 cores of the loaded workstation, 01:43 to 07:52;
`tests/Alltest M4long` with `M4LONG_RUNS=ended` evaluates them in 3 min:
12 PASS and 2 FAIL of M4long, plus the gate; exit 1; re-evaluated after
the M4 review, round 1, with the recorded runs revalidated and the gate of
acceptance 5 restated by the researcher: section 26, items 5 and 11):
   - M4.3: n' t_ev on the 123 native bins of 0.25-0.45 m = 3.40681e-4 mol
     s/m at 7 s (t_ev 6 s; max |ratio - 1| = 3.5e-5) and 3.40682e-4 at 16 s
     (t_ev 15 s; 1.1e-6), against n0/U_b,d = 3.40682e-4: PASS (the paper's
     3.38e-4 is 0.8 % lower; n0/U_b with the analytic U_b is 3.4241e-4).
   - M4.4: the solid profiles at 12 s (t_ev 6 s, 0.17 % of n0 still gas)
     and 30 s (t_ev 15 s) differ by 1.7e-3 of the deposit (L1; the largest
     bin by 1.4e-3 of the peak): PASS; T_dep 679.824 K in both.
   - M4.5, both rules, T_dep at the end of the experiment (or its bound):
     T106_6 774.68 / 774.69 K (paper 775: -0.32 / -0.31 K), T106_60 712.08 /
     712.02 K (713: -0.92 / -0.98), T106_300 666.89 / 666.77 K (666: +0.89
     / +0.77), PASS; T5_60 789.25 / 789.22 K (787: +2.25 / +2.22), FAIL by
     0.25 K for both rules; T540_60 670.00 / 670.05 K (670: +0.00 /
     +0.05), PASS.  The rules differ by at most 0.12 K in T_dep, so both
     match where one does: the criterion does not tell the rules apart.
   - M4.5b: T5_60 at dt 10 and 5 ms: 789.25415 K both (difference 2e-7
     K), PASS: dt 10 ms is adequate at 5 mL/min.
   - T5_60 in detail: the onset lies one native bin above the input Tdep
     = 788 K (Table 1).  The last bin without deposit has its wall at
     789.37 K, the first with deposit at 786.91 K and already 22 % of the
     peak (the slow gas deposits within the first bin below Tdep), so the
     1 % crossing, interpolated between the bin centres as the paper
     prescribes, lands at 789.25 K, 1.25 K above Tdep.  The wall falls by
     19 K/cm there (3.1 K per bin), the steepest of the five onsets, so the
     criterion's +-2 K is at the resolution of the binning; the paper's
     787 K is 1 K below its Tdep.  Both rules agree to 0.04 K, and the dt
     test shows no time-step effect.  No parameter was changed: a
     difference of +2.2 K is reported as a FAIL of acceptance 5 for T5_60.
     (Section 26, item 5: under the gate restated by the researcher all
     10 lines pass; the bin of the onset straddles the table node at 0.41
     m and its drop is 2.77 K; the set-up differs from the paper's T5_60
     run in the wall profile, the source position and the state compared,
     and the evidence does not single out one cause.)

---

## 26. Amendments from the M4 review, round 1 (signed off by the researcher, 2026-10-03)

The M4 verification of round 1 failed on one criterion: acceptance 5 for
T5_60 (+2.25 K with layer distance, +2.22 K with faceCells, against +-2 K;
the other 8 lines passed).  The researcher restated that gate on
2026-10-02 (item 5).  The review found three major defects of `layer
distance` (two reviewers found the same defect with several interface
patches; one found the dependence on the decomposition), three major
points of the science (the case is not what Figs. 8-10 show; acceptance 5
is not compared like for like; the SuperBee pin cites code that does not
use it), and fourteen minor points and nits.  Each was assessed first.
All were real; one nit (the round-off fluxes of the pipe's wall faces) is
answered with evidence instead of a change of the carrier (item 3c).
New criteria: M4.1h, M4.1i, M4.1j, M4.7f, M4.13 (2 lines), M4.14 (2
lines) and M4.16 (M4long, 2 lines, not a gate); changed: M4.1d, M4.2,
M4.4 (text), M4.5 (the restated gate), M4.6, M4.7d, M4.11 (the excerpt
of the case as distributed now), M0.3g.  The binary of the round-1
sources fails every new or changed gate that a change of the solver or the
utility addresses (item 15).

1. **`layer distance` searches the exact distance over all interface faces
   of all ranks (major; sections 2 'layer distance', 12 item 10, 25 item
   2).**
   - Defects.  Section 25 built the layer with patchDataWave<
     wallPointData<label>> and turned `cellDistFuncs::useCombinedWallPatch`
     off, because the combined wall patch stores labels of the combined
     patch where patchDataWave expects mesh faces.  The per-patch
     correction that this selects, `correctBoundaryFaceCells`, runs once
     per interface patch: for a cell with faces on two patches, every
     patch overwrites its distance (the last one wins) while
     `Map::insert` keeps the face of the first, so a corner cell was
     dropped from the layer or mapped to the farther face, depending on
     the order of the patches in the boundary file (the reviewers' duct:
     40 or 80 cells, the deposit of 40 corner cells on the farther patch;
     their 2-D box: 9 or 10 cells).  Either correction also sees only the
     wall faces of the cell's own rank: a cell that touches the wall only
     by a point or an edge whose wall faces lie on another rank keeps the
     wave distance (to the nearest face centre), which is larger, so the
     layer depended on the decomposition (a tetrahedral pipe: 1830 cells
     in serial, 1829 on 4 ranks, and the serial state restarted on the 4
     ranks was refused, -4.2e-6 relative).  The consistency check of
     section 25 (a nearest face farther than the distance allows) could
     not see a distance that is too large, or a valid face that is not the
     nearest.  The structured meshes of the M4 tests (one interface patch,
     no point contact) could not show either.
   - Fix (`readDistanceLayer()`): the layer is the definition of section 2
     itself.  d_c is the exact distance of the cell centre from the
     nearest face of the interface patches (`face::nearestPoint`, exact
     for a planar face; the triangles about the face centre otherwise),
     over the faces of all patches and all ranks; f(e) is that face.
     Every rank receives the interface faces whose bounding box overlaps
     the bounding box of its cell centres grown by the search radius 2
     distance (PstreamBuffers: the face's global index and its points in
     their own order), builds an octree over them (indexedOctree of
     treeDataPrimitivePatch) and computes the distances of the faces
     within the radius of each cell centre.  Ties (distances within 1e-10
     relative, e.g. a centre nearest to a point or an edge shared by
     faces) go to the face with the lexicographically smallest centre,
     computed from the face's own points: a boundary face keeps its
     points and their order in every decomposition, so the cells, their
     nearest faces and A_I/V_I depend neither on the order of the patches
     nor on the decomposition (but for the last bits of the centres of
     cells at processor faces: a cell exactly at the distance; the start-
     up margins show how far the rule is from that).  The radius of 2
     distance also gives the nearest cell outside the layer for the
     margin line ('more than <2 distance> m' when there is none).
     patchDataWave, the switch and the consistency check are gone.  Both
     counts of section 12, item 10 are still printed: the layer of the
     exact distance and the centre-distance rule R - r of T-Flows.
   - The pipe is unchanged: 30528 cells, the same distances (28.5-88.6
     um, the next cell at 154.7 um), and the round-1 binary and this one
     write bit-identical mDep_, Y, deposits, regimes and ledger in 2
     steps with a deposit in 9432 cells (only mw_ differs, item 2), so
     the nearest faces of all 30528 cells are those of round 1.  The 13
     long runs are revalidated on that basis (item 11).  The reviewers'
     cases with this binary: the duct 80 cells in both orders, the deposit
     only on the nearer patch (WALLZ 0 on every face); the 2-D box 10
     cells in both orders; the tetrahedral pipe 1830 cells in serial and
     on 4 ranks (21 elements with their nearest face on another rank), the
     serial state restarted on the 4 ranks accepted.
   - Cost: the whole start-up of run/014 in serial takes 1.64 s instead of
     1.49 s (the search about 0.15 s; 260 MB instead of 294 MB).  Every
     rank holds the interface faces within reach of its cells (on the
     pipe with 16 slabs about 1/16 of the 15264 wall faces, each with its
     points).
   - Tests.  New M4.13 (`tests/paperMode/corners`: an L-shaped channel of
     75 cells with the interface patches WALLA, WALLB and WALLC, which
     meet at a concave corner (WALLB, WALLC: a corner cell 0.05 mm from
     WALLB and 0.1 mm from WALLC) and at a re-entrant corner (WALLA,
     WALLB: a cell that touches both only by the corner edge, 0.1118 mm
     from it), the patches in both orders of the boundary file): distance
     7.5e-5 m, 5 cells and A_I/V_I = 20000 1/m in both orders, mDep_ 0 on
     every face of WALLA and WALLC and > 0 on every face of WALLB, the
     same values in both orders; distance 1.25e-4 m, 15 cells and 6666.67
     1/m in both orders, mDep_ the same in both orders and 0 on the
     corner face of WALLC.  New M4.14: (a) the corners at 1.25e-4 m on 3
     ranks, one block each, so that the rank of the cell at the
     re-entrant corner holds no interface face: 15 cells, inventories
     equal to serial to 1e-12, mDep_ (reconstructed) to 1e-12, 1 element
     with its nearest face on another rank; the serial state at 0.001
     restarted on the 3 ranks and the 3-rank state reconstructed and
     restarted in serial, both exact; (b) a tetrahedral pipe made by gmsh
     (R 2.4 mm, L 10 mm, h 0.6 mm; NOTRUN without gmsh), layer 4e-4 m:
     the same cells, A_I and V_I in serial and on 4 ranks (simple (1 2
     2)), the number of cells recomputed by `checks.py exactLayer` from
     constant/polyMesh and the cell centres with the exact point-polygon
     distance, inventories equal to 1e-9, and the serial state restarted
     on the 4 ranks exact.  M4.2 reads the new wording of the comparison
     line ('0 only in the layer').

2. **mw_ holds the deposits of layer distance exactly (minor; section 25
   item 2, interfaceExchange.H).**
   - Defect: `updateDepositField()` wrote the cell value of mw_ as (m''_e
     A_e)/A_e, which is 1 ulp off m''_e for many values.  A restart on the
     same decomposition compares mw_ with the deposits of
     uniform/phaseChange/ at tolerance 0 (a binary field without a
     repeated value) and warned in every restart of the long runs that
     read the exact deposits (T106_60-distance, T106_300-distance,
     T540_60-distance; e.g. 'differs ... by up to 4.33681e-19 kg/m2'),
     and a restart after decomposePar or reconstructPar, which re-reads
     the deposits from mw_, re-read them rounded.  layer faceCells was not
     affected (its patch values are the deposits).
   - Fix: with layer distance the cell value of mw_ is pair.deposit of the
     cell's element, assigned (one element per cell).  mw_ is output and
     the fallback of a re-decomposed restart; the computation is
     unchanged (item 1: everything else bit for bit).
   - Tests: M4.7d (the channel, serial and 4 ranks) also requires no
     warning that mw_ differs from the deposits; new M4.7f restarts the
     serial pipe of M4.7c from 0.02 (deposits in 9432 cells that do not
     all survive the round trip): every file at 0.04 cmp-identical to the
     continuous run, the restart exact (relative 0), no warning.

3. **setFrozenCarrier (minor points and nits; section 25 item 1).**
   a. *-check and the written time (minor).*  `-check` read the phi of
      the start time, and run/014 starts from latestTime: after a run with
      `writeFrozenFields no` the decomposed check stopped ('cannot find
      file .../processorN/0.02/phi'), and there was no -time option.  The
      utility wrote the carrier into the start time as well, i.e. into the
      latest time of a case that had run.  Now `-check` reads the phi that
      the solver reads, the instance that `findInstance` finds from the
      start time (0 for a carrier the solver does not write), and the
      log line names it ('phi read at time 0 (findInstance from 0.02)');
      the carrier is written into the instance of rho found from the
      start time, or into the start time when there is no rho yet.  The
      Time is set to that instance before the writes, because
      `regIOobject::write()` moves an object whose instance is an older
      time to the current time (a first version of the fix wrote into
      0.02 although the fields were created at 0).  New M4.1j: the wedge
      after the 20 steps of M4.1g.  The run/014 readme runs `mpirun -np 16
      setFrozenCarrier -check -parallel` after decomposePar.
   b. *The wall flux (minor).*  The exact flux of the axial profile
      through a wall face is 0 only when the face is parallel to the
      axis.  On a tetrahedral mesh the walls carried 9.5e-4 of the inflow
      (the reviewer's tet pipe), and the solver's paired species, whose
      wall condition must be zeroGradient, was carried through the wall
      (booked as WALL transport); the log's net flux through the patches
      could not show it.  Now every continuity line is followed by a 'Wall
      flux' line (the sum over the wall faces of |phi| relative to the
      inflow, and the largest face), and the utility stops without writing
      when the sum exceeds 1e-12 of the inflow, unless `permeableWalls
      yes` (then a warning).  An exactly conservative construction that
      makes any wall impermeable (face fluxes from an edge vector
      potential, A = U_b (1 - r^2/(2 R^2)) (axis x (p - o)), with the
      edges on the wall given (Q/2 pi) dtheta) would replace the exact
      face integral that deliverable 1 asks for, and the meshes of M4-M8
      (the gmsh pipe, blockMesh wedges; M8 uses physical carriers) have
      walls parallel to the axis; it is recorded here as an option, not
      built.  New M4.1i: a wedge whose wall points are moved radially
      with x (wall faces tilted): the utility stops and writes nothing;
      with permeableWalls yes it writes, warns, continuity <= 1e-15, the
      wall flux reported above 1e-12 of the inflow.
   c. *The round-off fluxes of the pipe's wall (nit; not changed).*  The
      reviewer proposed to set phi = 0 on the 'walls' patches.  The wall
      faces of the gmsh pipe are tilted by |S_x|/|S_f| <= 1.5e-15: the y
      and z of a wall point differ by up to 3.0e-18 m between the two
      axial ends of a face (measured on constant/polyMesh), and the exact
      flux through these faces is +-1e-25 kg/s, 8.4e-16 of the inflow in
      total.  Deliverable 1 asks for the exact integral over every face,
      and the continuity of the cells next to the wall holds with these
      values; setting them to 0 would replace the exact integral by an
      approximation of the same size and change the carrier of run/014
      (and of the 13 long runs) in its last bits.  The solver books the
      flux as transport through WALL (4.1e-25 kg in 12 s).  They are kept,
      and the new wall-flux line reports them (8.38e-16 of the inflow,
      largest face 9.77e-25 kg/s).
   d. *Non-extruded meshes untested (minor).*  The wedge and the pipe are
      extruded along the axis, where even a face-centre flux is
      divergence-free (7e-17 of the inflow on the M4.1c wedge).  New
      M4.1h: the wedge with its inner points moved radially, r' = r (1 +
      0.15 sin(pi r/R) sin(2 pi x/L)), and axially, x' = x + 0.15 mm
      sin(pi x/L) sin(pi r/R) (the axis, the wall, the inlet and the
      outlet unchanged; checkMesh 'Mesh OK'): continuity 1.06e-16 of the
      inflow in serial and 1.07e-16 decomposed (simple (4 1 1)), the wall
      flux 0, while the face-centre flux u(C_f) S_f,x of the same profile
      misses continuity by 4.1e-5 of the inflow (`checks.py
      faceCentreRule`).
   e. *decomposePar and ordinary patches (documentation gap).*  Section
      15 item 7 and the readme named the 10-digit 'uniform' values of
      processor patches only.  decomposePar writes any patch whose values
      are all equal that way, also an ordinary patch with a single face on
      a rank (the verifier's interleaved decomposition put one INLET face
      on rank 3: 2.56e-13 of the inflow).  M4.1d now decomposes the wedge
      with a single INLET face on rank 1: decomposePar writes it as
      'uniform -1.501287666e-09', -check reports 3.97e-13 of the inflow,
      and the carrier written with -parallel on that decomposition is
      exact (6.0e-17); the readmes recommend the check (or -parallel)
      after every decomposition, e.g. scotch on Merlin7.
   f. *80 columns (nit).*  The six new lines over 80 columns
      (setFrozenCarrier.C: 80, 103, 107, 440, 810; interfaceExchange.C:
      1263) are re-wrapped or gone; the only source lines above 80
      columns are again the two inherited lines of rhoFixedFlowFoam.C
      (section 23, item 5).

4. **Build and test plumbing (minor, nit).**
   - M4.6 copied the sources with `--exclude="./Make/*/"`, which matches
     nothing in GNU tar: the copy carried the v2412 objects, and a v2606
     installation with the same WM_OPTIONS (Opt) relinked them without
     compiling anything; M4.6 passed without checking that each source was
     compiled (the reviewer: exit 0, 0 compile lines).  Now the copy
     excludes `./Make/linux*` and `./utilities/*/Make/linux*`, runs wclean
     in the copy, and M4.6 requires every source of both Make/files to be
     compiled, as BUILD does.
   - Allwmake decided whether to clean from the solver's objects and
     stamps only: after a wclean of the solver alone (which removes the
     stamps), objects of the utility from another installation with the
     same WM_OPTIONS were relinked.  Objects of the utility now count:
     without the stamps they are 'unknown' and both are cleaned.  M0.3g
     has the case (stand-in wmake and wclean).

5. **Acceptance 5, as restated by the researcher (decision of
   2026-10-02; plan section 6 M4, item 5; section 25 items 4a, 4e, 4h and
   5).**
   - The gate: |T_dep - T_paper| <= max(2 K, dT_bin), dT_bin the drop of
     the wall temperature across the native bin that holds the onset (the
     Twall column at the bin edges, interpolated linearly between the bin
     centres).  `checks.py tdep` and `tdepDecided` print, per case and
     rule, the raw difference, dT_bin, the verdict of the gate and the
     plain +-2 K comparison for information; M4long re-evaluates the
     recorded runs (no rerun).
   - Results (both rules, distance / faceCells): T106_6 -0.32 / -0.31 K
     (dT_bin 2.64 K); T106_60 -0.92 / -0.98 K (2.44 K); T106_300 +0.89 /
     +0.77 K (2.44 K); T5_60 +2.25 / +2.22 K (2.77 K); T540_60 +0.00 /
     +0.05 K (2.44 K).  All ten lines meet the gate; T5_60 is outside the
     plain +-2 K.  Section 25 quoted '3.1 K per bin' for T5_60: the
     gradient of the table upstream of x = 0.41 m (19 K/cm) times the bin;
     the bin of the onset (0.4085-0.4101 m) straddles the table node at
     0.41 m, beyond which the table falls by 12 K/cm, and its drop is 2.77
     K.
   - The cause of the T5_60 offset, as far as the evidence goes.  The
     offset of +2.25 K from the paper's 787 K has two parts (corrected in
     section 27, item 11): +1.25 K from the binning, relative to the input
     Tdep = 788 K: the last bin without deposit has its wall at 789.37 K,
     the first with deposit holds 22 % of the peak, so the 1 % crossing,
     interpolated between bin centres as the paper prescribes, lies 1.25 K
     above Tdep (no time-step effect, M4.5b); and +1.0 K because the
     paper's own value lies 1 K below its input (787 against 788 K), a cause
     that this case cannot establish.  The set-up
     differs from the paper's T5_60 run in ways of the same order (the
     science review): (i) the wall profile: the paper measured one for
     each operating condition (Section 2.4) and plots it in Fig. 14 iii,
     where it falls by about 15 K/cm at 788 K and lies below this table
     from the inlet on (about 316 K at x = 0, a maximum of about 1050 K
     near 0.13-0.15 m, 793 K at 0.40 m, 755 K at 0.425 m, digitised); a
     different gradient at the onset changes the bin offset; (ii) the
     source position: the arrows of Figs. 13-14 mark 'the initial position
     of the PbI2 source term used in the simulations', about 0.13 m for
     T106 and 0.15 m for T5 and T540, while run/014 releases at 0.04-0.05
     m for every case; (iii) the state compared: the paper's
     temperature-based T5_60 curve integrates to about 0.98 n0, a
     complete deposit (peak 5.1e-5 mol/m at 0.435 m, deposit over
     0.406-0.546 m), while run/014's T5_60 at 60 s holds 33 % of the
     inventory on the wall and 67 % in the gas, spread back to the closed
     inlet (1.1e-6 mol/m at x = 0.001 m), and its deposit reaches 22 % of
     its peak in its first bin.  None of these is established as the
     single cause.
   - Decision needed (sign-off): whether T5_60 is evaluated at the end of
     the experiment (60 s, as now: T_dep of the deposit made by then,
     Table 1's duration) or at completion (all gas deposited, the state
     the paper's curve shows).  (Corrected in section 27, item 10: the
     choice changes the shape of the profile, not T_dep or the verdict.
     The bins upstream of the first bin with a layer cell below Tdep can
     never deposit, so the onset cannot move upstream (T_dep <= 789.37 K),
     and the deposit of 30-60 s grows self-similarly, so at completion the
     onset stays where it is (to about 0.01 K).  The run to completion,
     several hundred seconds at 5 mL/min, 1.5e4 to 3e4 more steps of 10
     ms, about 2 to 5 h on 8 ranks of this machine, is optional; the
     bound of `tdepDecided`, which lets every bin take the remaining gas,
     stays conservative and is not sharpened.)
   - Section 25 item 4e said 'The paper does not say whether Table 2 used
     per-case tables'; Section 2.4 says that each operating condition had
     a dedicated wall-temperature profile, and Figs. 13-14 plot them, so
     the Table 2 runs very likely used them (corrected).  Item 4h said that
     T_dep does not depend on the release position; for the transient
     T106_6 and the incomplete T5_60 it does (corrected).  Item 4a said
     that the 1 % onset 'lies just below Tdep'; it lies within a bin of
     Tdep on either side (666.89 against 666 K, 712.08 against 712 K,
     789.25 against 788 K; the paper's own T106_60 is 713 against 712 K),
     at the resolution dT_bin of 2.4-2.8 K per native bin (corrected).
   - The request for the per-case wall tables (`manuscript/notes-for-
     authors.md`, item F; the file is read-only here) should be extended
     to the source positions of Figs. 13-14 and the state (time) of the
     Table 2 curves (item 14).

6. **Figs. 8-10 are not reproduced (major science; section 25 items 3, 4a
   and 4b).**
   - Section 25 and the readme called run/014 'the reference run of
     Figs. 8-10' without comparing it with them.  Measured now (new
     M4.16, not a gate, on the native profiles of the recorded run
     fig9-tev15, against the paper as digitised by the science review):
     the gas plateau at 16 s 2.2712e-5 against 2.26e-5 mol/m (+0.5 %,
     acceptance 3); the gas at 50 % of the plateau at 0.4987 m against
     0.507 m (-8.3 mm), at 10 % at 0.5212 m against 0.518 m (+3.2 mm);
     the solid onset at 30 s at 0.4761 m (T_dep 679.82 K) against 0.502 m
     (-25.9 mm), the peak 9.11e-4 mol/m at 0.4955 m against 1.14e-3
     mol/m at 0.513 m (Fig. 10; -17.5 mm).  In this wall table 680 K lies
     at 0.476 m and 641 K at 0.50 m: the published onset corresponds to
     about 640 K.  Both T-Flows temperature-model codes use Pb_t_dep =
     638.15 K (`Thesis:User_Mod/Source.fpp:25`, `Project:User_Mod/
     Source.f90:28`) and the layer r > 0.9 R (`Thesis:User_Mod/
     Source.f90:114`, `Project:User_Mod/Source.f90:111`), 61056 cells on
     this mesh, the same as `layer distance; distance 2.4e-4`; `Thesis`
     uses A 150 1/s and k 0.02 1/K (`Thesis:User_Mod/Types.f90:7-8`),
     `Project` A 1.35 and k 2.7 (`Project:User_Mod/Types.f90:6-7`;
     corrected in section 27, item 8: round 1 gave A 150 and k 0.02 to
     both).
   - The T-Flows parameter set as a sensitivity variant (new long run
     `fig9-tev15-TFlows`: Tdep 638.15 K, A 150, k 0.02, distance 2.4e-4,
     otherwise run/014 at t_ev 15 s; 7500 steps of 4 ms (Courant 0.45; a
     4 ms run of run/014 to 8 s, the science review's, gives the gas front
     of the 1 ms run at 8 s and 16 s to 4 digits: 50 % at 0.4987 m, 10 %
     at 0.5212 m), 8 ranks, 44 min (its first step ends at the corrector
     limit at an initial residual of 2.1e-5, as the first step of run/012,
     section 16, item 7); not a gate): the solid onset at 30 s at
     0.5011 m (T_dep 639.14 K) against 0.502 m (-0.9 mm, within a bin),
     its peak 1.62e-3 mol/m at 0.5102 m against 0.513 m (-2.8 mm), all of
     n0 deposited; the gas front at 16 s 50 % at 0.5134 m and 10 % at
     0.5272 m against 0.507 and 0.518 m (+6.4 and +9.2 mm); the plateau
     n' t_ev 3.40682e-4 mol s/m (n0/U_b,d to 3.5e-6, as run/014).  So the
     T-Flows set puts the deposit where Figs. 9-10 show it (onset and peak
     within two bins), while the gas front lies 6-9 mm downstream of Fig.
     9's and the peak is 42 % higher than Fig. 10's (Fig. 9's spike, 3.5e-3
     mol/m, is higher still; notes-for-authors, item D): the published
     figures are consistent with Tdep of about 638 K and the thesis layer,
     not with the 680 K the paper states for them.  Neither set is tuned
     here; the question goes to the authors (item 14).
   - run/014 as distributed is now the configuration of Figs. 8-10 as the
     paper states it, t_ev = 15 s to 30 s (`times (16 30)`), the one Figs.
     8, 9 and 10 share (it was the t_ev = 6 s run to 12 s); the readme
     calls it that, quantifies what it does not reproduce, and states for
     sign-off that the figures appear to use Tdep of about 638 K (and
     possibly the thesis layer).  Section 11's 'Tdep = 680 K holds for
     Figs. 8-10' stays the paper's statement; the runs of acceptance 3-5
     are unchanged (the edits of longRuns set every value they change).
     The out_excerpt of run/014 (M4.11) is regenerated: the release rate
     n0/t_ev and the new start-up line change its numbers, not its
     structure.

7. **Fig. 9's evaporation times (minor; section 25 item 4b).**  Fig. 9
   shows t_ev = 15, 30, 45, 60, 75 and 90 s (gas at t_ev + 1, solid at 2
   t_ev); it has no 6 s curve, and Figs. 8-10 share t_ev = 15 s.  Item 4b
   and the readme are corrected.  M4.4 (t_ev 6 and 15 s) is a check of a
   property of the model (the final deposit does not depend on t_ev,
   which Fig. 9 shows for 15-90 s), and its text says so; the 15 s / 30 s
   pair of Fig. 9 would cost another 60000 steps of 1 ms (about 5 h on 8
   ranks of this machine at its present load) and is not run.

8. **SuperBee (major science; nit).**
   - The run/014 readme and fvSchemes cited `Transport:control:46` and
     `Thesis:control:45` (`ADVECTION_SCHEME_FOR_SCALARS superbee`) for the
     species.  The GPU solver of every LESTO branch never reads that
     keyword (`Gpu:Read_Controls_Mod/Numerical_Schemes.f90:101-105`: 'Line
     reserved for advection scheme'; SuperBee exists only under
     `Sources/Process/Cpu`); it blends central, linear-upwind and upwind
     with `BLENDING_COEFFICIENTS_FOR_SCALARS` (`Gpu:Field_Mod/Core/
     Add_Advection_Term.fpp:41-46`, ramped in over 240 steps), whose
     default depends on the branch (corrected in section 27, item 8):
     0 0 1, pure upwind, on advanced_scalar_model and
     scalar_transport_francesco (`Sources/Shared/Control_Mod/Numerics/
     Blending_Coefficients_For_Scalars.f90:19`; HKS, Thesis, Transport),
     and 1 0 0, pure central differencing, on multiple_scalars_gpu_francesco
     (Project; its `Add_Advection_Term.fpp:95-100` weights the central face
     value with the first coefficient), for the scalars and for the energy
     (`Blending_Coefficients_For_Energy.f90:19`, the same defaults): the
     Thesis, Project and Transport controls set no blending, LESTO_HKS sets
     0.8 0.1 0.1 for the scalars (`HKS:control:45`), so the energy equation
     runs upwind in HKS, Thesis and Transport and central in Project.
     run/014 keeps SuperBee, the scheme the paper states (Section 2.2);
     the source column and the fvSchemes comment are corrected.  The
     plateau and T_dep do not depend on it; the shape of the profiles does
     (upwind's numerical diffusion U dx/2, about 1.4e-4 m2/s on the axis,
     is comparable to D, about 1e-4 m2/s at the onset), so the profile
     comparisons of M6 should be made with an upwind, a blended (0.8 0.1
     0.1) and a central variant as well, and the scheme is a question for
     the authors (item 14).
   - The paper pins SuperBee for the temperature as well (Section 2.2);
     setFrozenCarrier uses linearUpwind (the steady Picard iteration with
     SuperBee does not converge: 200 iterations, a change of 29 K left).
     Upwind, linearUpwind and central differ by at most 0.023 K (mean
     0.003 K) in the layer cells at x = 0.38-0.52 m (the science review),
     so the effect on T_dep is negligible; the readme and
     setFrozenCarrierDict record the deviation.

9. **The geometry of the T-Flows temperature-model cases (minor).**  They
   start the domain at the sample position: `L1 = 0.04` ('sample
   position') to `L2 = 0.69` in `Thesis:lesto_pipe.geo:16-17`,
   `Project:Stainless_Steel.geo:16-17` and `Transport:LESTO_pipe.geo:
   16-17`, with wall tables from 0.04 m and the release into the cells of
   the inflow faces (`Thesis:User_Mod/Source.fpp:61-73`).  run/014 uses
   the HKS domain from x = 0 with the sample at 0.04-0.05 m, so its gas
   can diffuse upstream of the sample (1.1e-6 mol/m at x = 0.001 m in
   T5_60 at 60 s), and the rear flank of its gas at 16 s lies 3-7 mm
   downstream of Fig. 9's (0.102/0.135/0.167 m at 10/50/90 % of the
   plateau against about 0.099/0.129/0.160 m in Fig. 9, whose x axis
   starts at 0.04 m), consistent with the release positions: T-Flows
   injects into the cells of the inflow faces at x = 0.04 m, run/014 over
   the boat, x = 0.04-0.05 m, centre 0.045 m (corrected in section 27,
   item 12: round 1 said the cloud 'starts slightly earlier').  Their mesh
   has 424 axial cells over 0.04-0.69 m (dx 1.533 mm against run/001's
   1.627 mm; `Thesis:lesto_pipe.geo:14-17,40,482`), so their native bins
   quantise T_dep and dT_bin differently (section 27, item 14).  Recorded as a parity deviation in the readme
   (domain row; 'Transport: control file only' corrected) and here; not
   changed (the paper's Table 2 runs release at 0.13-0.15 m, item 5, so
   no single domain serves both).

10. **The pin of the rate per layer cell (nit).**  The readme recorded
    `areaPerVolume layer` for both rules; run/014 sets none, so layer
    faceCells (the one-line switch, and the faceCells runs of M4long) uses
    the default `cell` (A_e = |S_f|), whose rate per cell differs from A
    f(T) by about 1e-6 on this mesh.  The readme records the default of
    each rule, and the dictionary says so in a comment.

11. **The long runs: key, binaries and revalidation (the researcher's
    decision: re-evaluate the recorded runs, no rerun).**
    - The key of a run hashes run/014's inputs without comments now (as
      the sources already were), so a corrected comment no longer
      invalidates a run.
    - The binaries of the runs are kept per key in `long/bin/<key>/`
      (the round-1 layout, one directory, is moved there), and every run
      records the key of the binaries it was made with (`m4bin`; for the
      13 recorded runs it was written before this round's changes, after
      checking that their keys and that of long/bin matched the
      round-1 sources).
    - A run that ended under another key is revalidated instead of run
      again: (1) its spec set up now (run/014 and setFrozenCarrier) must
      give its dictionaries without comments and endTime and its carrier
      bit for bit; (2) the binaries it was made with and the current ones,
      each run for 5 steps from a copy of its state at its last write
      before its end (on its own decomposition), must write the same
      files (cmpFields.sh without uniform/time and mw_PbI2_s, item 2) and
      the same balance and profile files.  Then its key is updated and the
      record kept in `m4revalidated`.  All 13 runs were revalidated, each
      with 5 steps from its last write before its end: fig9-tev6 from 11 s
      (4 ranks) and fig9-tev15 from 29 s (8 ranks), T106_6 (both rules)
      from 5 s in serial, T106_60 and T106_300 from 13 s (4 ranks; T106_60
      faceCells in serial), T5_60 (both rules, and dt 5 ms) from 59 s (4
      ranks), T540_60 from 6.99974 s (distance, 8 ranks) and 6.99952 s
      (faceCells, 4 ranks): 10, 40 or 80 files identical per run, and the
      balance and profile files identical (`m4revalidated` of each run).
      A trial on copies of two runs preceded it.  M4long then re-evaluated
      them (item 5).
    - `M4LONG_REVALIDATE=no` switches it off; a run whose binaries are not
      kept is run again (or prepared, with M4LONG_RUNS=ended).  (Section
      27, item 5: the probe cannot see a change that acts only before the
      run's last write, e.g. in the release window; revalidation is opt-in
      since round 2, and a revalidated run is named so in its criterion
      lines.)

12. **Statements corrected.**
    - Section 25 item 2: T-Flows' 2 pi R L/V_I is 8544.5 1/m (8544.52 on
      this mesh), not 8544.6.
    - Section 25 item 6 said that with every corrector solved to the
      floor the M4.7c iterates 'depend on the decomposition only through
      round-off'.  With SuperBee the limiter amplifies round-off: per cell
      the gas of serial and hierarchical (2 2 1) differs by up to 2.7e-5
      relative and mDep_ per face by 1.9e-7 of its maximum (the
      verifier's uneven split: 4.8e-5 and 4.7e-7), while the inventories
      agree to 1e-10; an upwind control run agrees to 4.7e-14 per cell and
      5.8e-15 per face.  The criterion (inventories to 1e-9) is unchanged;
      the statement is qualified.
    - 'An uneven hierarchical split' (section 6 M4): hierarchical (2 2 1)
      on the pipe has 3, 2, 2 and 3 processor patches per rank, and the
      manual layersApart split is uneven in cells (M4.7c); both run.
    - The solver readme listed 'layer distance (not yet implemented, M4)'
      among the start-up stops and described M2.6f as the 'not yet
      implemented' checks of model HKS, layer distance and mode inventory;
      both are corrected (section 25, item 8 says what M2.6f checks now).

13. **Observation for M8 (corrected and fixed in section 27, item 1).**
    In the T5_60 long runs the closure drifts to 7.6e-13 of the reference
    (distance, 6000 steps of 10 ms; faceCells 7.5e-13; dt 5 ms 5.4e-13 in
    12000 steps; 3.7e-13 in the verifier's 8-rank rerun); the other long
    runs stay at or below 8.2e-14 (measured from their balance files).
    Not 'with the steps' (the 5 ms run has twice the steps and no larger
    drift): it grows linearly with the time that gas or deposit stays in
    place, from double-precision imbalances of the discrete step that the
    ledger did not book (the column sums of the transport matrix, the
    storage products, the rounding of the stored deposits).  They are
    booked since round 2.

14. **Questions for the T-Flows authors (additions; `manuscript/
    notes-for-authors.md` is read-only here).**
    - D: were Figs. 8-10 made with Tdep = 680 K (Section 3.1) or with the
      638.15 K of both temperature-model codes (and their layer r > 0.9 R;
      A 150, k 0.02 in Thesis, A 1.35, k 2.7 in Project's code; section 29,
      item 6: Project's reference.txt lists A 150, k 0.02 for its
      reference case)?  With 680 K and
      this wall table the deposit starts at 0.476 m, not at 0.502 m (item
      6).  Evidence for the authors (section 27, item 14): both codes
      release 3.17e-5 mol, not the printed 3.16e-5; their U_MAX = 0.176
      m/s is the centreline velocity (Q = 95.5 mL/min), whose plateau,
      2.40e-5 mol/m at t_ev 15 s, Fig. 9 (2.26e-5) rules out, while
      FLOW_RATE 1.67e-6 with T-Flows' face fluxes reproduces it (section
      27, item 9); and their mesh has 424 axial cells over 0.04-0.69 m (dx
      1.533 mm, not run/001's 1.627 mm).
    - F: besides the wall tables of 5, 106 and 540 mL/min, the source
      positions of the Table 2 runs (the arrows of Figs. 13-14: about 0.13
      m and 0.15 m), the domain (from x = 0.04 m as in the temperature-
      model cases, or from 0), and the time of the Table 2 curves (the end
      of the experiment, or completion; T5_60's curve holds about the
      whole inventory).
    - E: the advection scheme of the published runs: SuperBee as Section
      2.2 states, or the blending of the GPU solver: by default pure upwind
      on advanced_scalar_model and scalar_transport_francesco and pure
      central differencing on multiple_scalars_gpu_francesco, for the
      scalars and the energy equation; 0.8 0.1 0.1 for the scalars in
      LESTO_HKS (corrected in section 27, item 8).

15. **Tests and negative evidence.**
    - New: M4.1h, M4.1i, M4.1j, M4.7f, M4.13 (2 lines), M4.14 (2 lines),
      M4.16 (M4long, 2 lines, not a gate).  Changed: M4.1d (a single INLET
      face on a rank), M4.2 (the comparison wording), M4.4 (its text),
      M4.5 (the restated gate), M4.6 (the copy, every source compiled),
      M4.7d (no mw_ warning), M4.11 (the regenerated excerpt), M0.3g (the
      utility's objects).  `tests/Alltest` expects them; `tests/paperMode/
      checks.py` has positions, exactLayer, perturbWedge, faceCentreRule,
      wallFlux, mdepPatches, samePatches and patchValue, and tdep with the
      restated gate.
    - Negative evidence: `tests/paperMode/Allrun` with the binaries of the
      round-1 sources (`FOAM_USER_APPBIN` pointed at a copy of them; the
      v2606 compile not set): 21 PASS, 11 FAIL, 1 NOTRUN (M4.6).  The
      failures are those this round addresses: M4.13 (distance 7.5e-5: 4
      cells and A_I/V_I = 25000 in the boundary order A B C, the corner
      cell's deposit on WALLC in the order C B A; 1.25e-4: mDep_ differs
      between the orders), M4.14 (the corners on 3 ranks: 14 cells, 7142.86
      1/m, the inventories 3e-3 from serial, both restarts refused, -2.7e-5
      and -1.4e-3 relative; the tetrahedral pipe: 1829 cells on 4 ranks
      against 1830, the serial state restarted there refused, -3.4e-7 -
      its FatalError on some ranks left mpirun waiting, stopped by hand
      after 14 min), M4.7d and M4.7f (2 and 1 warnings that mw_ differs
      from the deposits), M4.1i (the utility does not stop on the tilted
      wall), M4.1j (-check stops: 'cannot find file .../0.02/phi'; the
      carrier is written into 0.02), M4.1h (no wall-flux line; the exact
      continuity on the moved wedge, 1.06e-16, held in round 1 as well:
      the criterion guards it against a regression, as the face-centre
      flux, 4.1e-5, shows), and by construction M4.2 and M4.11 (the new
      wording of the layer line, and the excerpt of the t_ev = 15 s case).
      M4.1d passes with the round-1 utility (decomposePar's behaviour, and
      the round-1 utility already wrote exact fluxes in parallel).  The old
      Allwmake does not wclean in the case of M0.3g (checked by hand with
      stand-in wmake and wclean: 0 calls of wclean, 2 with the new one).

M4 status after round 1: `tests/Alltest -fresh all M4long` (M0 to M5
with M4, and the evaluation of the long runs with `M4LONG_RUNS=ended`)
gives 233 PASS, 0 FAIL, 1 NOTRUN (M0.4, the compile on v2606 on Merlin7:
no remote access; M4.6 compiles a copy of the sources, setFrozenCarrier
included, against the local v2606 Debug build and requires every source to
be compiled) in 50 minutes at a machine load of 30-40; 0 warnings; the M0
gate identical to main.  `tests/Alltest M4` alone gave 41 PASS and 2 FAIL
before the comparison of the printed A_I/V_I = 6666.67 was set to its
printed precision (the layer was right); then the full run above.  All 14
long runs (the 13 of section 25, revalidated, and the sensitivity variant)
are reused; under the gate restated by the researcher all ten Table 2
lines pass (T5_60 at +2.25 and +2.22 K against dT_bin = 2.77 K).  The
tracked and untracked files of the repository (git status and the content
of all 350 files) are unchanged by the run; the branch is phase-change at
a43cbd2, without commits or stashes.  Pending: the researcher's sign-off
of sections 25 and 26, the decision of item 5 (T5_60 at 60 s or at
completion), and the questions of item 14 to the authors.

---

## 27. Amendments from the M4 review, round 2 (signed off by the researcher, 2026-10-03)

The M4 verification of round 2 passed every criterion (acceptance 5 under
the gate the researcher restated on 2026-10-02; section 26, item 5).  The
review found one major defect (the closure of the ledger drifted over long
runs, against E8), eleven minor points and nits of the code, the tests and
the parity record, and the verifier six notes.  Each was assessed first;
all review findings were real.  One note of the verifier was not (item
15).  A first fix agent of this round was stopped by a usage limit after
part of the work; every one of its edits was reviewed before it was kept
(item 16).  New criteria: M4.1k, M4.7g, M4.10 (a second line), M4.15 and
M4.17 (M4long); changed: M2.2h (closure 1e-14 instead of 1e-12), and the
excerpts of M2.8, M3.7, M5.3b and M4.11 (run/012, 013, 014 rewritten).

1. **The closure drifted over long runs (major; E8, `matrixDefect.H`,
   `recordSolve()` and `realise()` of `interfaceExchange.C`).**
   - Defect.  E8 claims the closure 'at round-off for any linear
     tolerance'; it held for short runs only.  `recordSolve()` booked R'
     with the row-form defect b - A Y of the transport-only matrix, the
     residual of the linear system, and three double-precision imbalances
     of the discrete step were booked nowhere.  Each is a few eps of the
     inventory and repeats alike in every step where the gas or a deposit
     stays in place, so the closure drifted linearly with time:
     a. the column sums: for a conservative operator the off-diagonal
        coefficients of a column add up to the diagonal minus the storage,
        but the assembled diagonal of ddt + div - laplacian is the rounded
        negSumDiag of the convection and the laplacian coefficients plus
        the rounded sums of the three matrices, so delta_c = diag_c - D_c +
        sum_{f: l(f)=c} lower_f + sum_{f: u(f)=c} upper_f is a few ulp of
        the diagonal, and the solved system creates sum_c delta_c Y_c per
        unit time.  The reviewer's diagnostic X = sum_c delta_c Y_c (long
        double, scratch) explains the T5_60 drift: on its state at 40 s, 4
        ranks, 100 steps of 10 ms, the closure changed by -1.758e-20 kg
        (-1.22e-16 of n0 M per step), X predicted -1.908e-20 kg, and
        booking X left no drift; 79 % came from the assembly of ddt + div
        - laplacian, 21 % from the laplacian's own negSumDiag;
     b. the storage products: the Euler scheme computes ((rDeltaT rho)
        V) Y and ((rDeltaT rho) Y^n) V, the ledger (rho Y) V, three
        roundings of the same mass; at a bitwise steady gas they differ in
        every step alike;
     c. the rounding of the stored condensate: `condensate -= dt q/A`
        changes the stored mass by the rounded update, not by dt q; a
        deposit that grows by the same increment step after step rounds
        alike in every step.
     The T5_60 long runs closed to 7.6e-13 at 60 s (6000 steps; the 5 ms
     variant, twice the steps, 5.4e-13: a drift with time, not steps), the
     other long runs to 1.8e-14 to 8.2e-14, and M2.2h to 1.6e-13 (linear
     while its gas stood bitwise steady during the release, 1.5-15 s, flat
     afterwards; under its gate of 1e-12).  Section 26, item 13 ('drifts
     with the steps') and the readme ('a random walk of the per-step
     round-off') misread it; at the T5_60 rate, M8's slow runs could reach
     its gate of 1e-12.  In M2.2h, fixed in three stages, the column sums
     (a) were about 7 % of the drift, the deposit rounding (c) about 85 %
     and the storage products (b) about 8 %; in T5_60 the column sums
     dominate (D dt/h^2 of the wall cells about 300).
   - Fix (a, b; `matrixDefect.H`, `recordSolve()`).  `cellDefect` forms
     the transport-only defect in FLUX form,

         r_c = (b - s)_c - (m_c - m0_c)/dt - sum_{f of c} F_f(Y),

     with F_f = upper_f Y_u - lower_f Y_l computed once per internal face
     and entered into its two cells with opposite signs, ic Y - bc on the
     non-coupled patches (the value of fvMatrix::flux(), as the ledger's
     transport), and on processor (and cyclic) faces the antisymmetric
     mean (F_own - F_other)/2 of the fluxes of the two sides (exchanged
     with syncTools::swapBoundaryFaceList; their coefficients are computed
     separately and are not exactly opposite); s is the source of
     fvm::ddt(rho, Y) (the old-time storage, the same doubles as in the
     solved matrix), m and m0 the masses of the cell as gasInventory()
     computes them, (rho Y) V and (rho Y^n) V (rho is frozen).  Every cell
     is accumulated in long double and the cells in a compensated sum.
     The face fluxes telescope exactly and the storage telescopes with the
     inventory itself, so the imbalances (a) and (b) are booked in
     solverDefect, as the stiff exchange's is (section 21, item 1):
     |solverDefect| grows by a few eps of the inventory per step, far below
     its gates.
   - Fix (c; `realise()`).  The rounding of every condensate update,
     (held - condensate measure) - dt q with both masses as the
     inventories compute them, is booked as `clamped`, the entry that
     already took the round-off residue of a CAPPED element.  `clamped` is
     therefore no longer 0 in runs without CAPPED elements (run/014's
     first 5 steps: 1.6e-57 mol; run/013's sample reservoirs: -7.2e-23
     mol).  The ledger's entries, layout and restart format are unchanged.
   - Only the booking changes.  From a fresh start through the release
     window (T5_60 on 4 ranks, 200 steps of 10 ms), the binaries before
     and after the fix write cmp-identical fields, per-rank state
     (`uniform/phaseChange/`), `mw_`, `mDep_` and profiles (44 files per
     time), and balance files identical but in the columns clamped,
     solverDefect and closure (`checks.py sameBalance`); the same in the
     gas and wall columns of all 136364 steps of M2.2h, in the 20 steps of
     the drift pipe (M4.15), and in the revalidation probes of 11 long
     runs (item 2).  The logs of run/012, run/013 and run/014 differ from
     their excerpts only in the digits of `solver defect`, `clamped`,
     `solverDefect` and `closure` (all 123, 111 and 126 compared lines
     checked token by token), and run/014's in the start-up count of item
     3; the three excerpts are rewritten from the new build (M2.8, M3.7,
     M5.3b, M4.11).
   - Results.  The T5_60 state at 40 s, 100 steps on 4 ranks: closure
     change exactly 0 (was -1.76e-20 kg), fields and profiles
     cmp-identical.  M2.2h: max |closure|/reference 1.1e-15 over its
     136364 steps, 9.5e-16 at the end (was 1.6e-13); its gate is tightened
     from 1e-12 to 1e-14.  New M4.15: the pipe of run/014 with its own
     settings, background Y = 0.1, 20 steps of 0.1 s on 4 ranks (D dt/h^2
     of the wall cells about 3000): 4.2e-16 (the binary before the fix:
     1.16e-14, linear), gate 2e-15.  The three T5_60 long runs, run again
     with the fix (item 2, M4.17): max |closure|/reference 8.4e-16 and
     8.6e-16 over their 6000 steps (distance, faceCells) and 8.0e-16 over
     the 12000 steps of the 5 ms variant (were 7.6e-13, 7.5e-13 and
     5.4e-13), the largest in the first seconds, during the release, none
     growing.
   - Not booked (observation): the rounding of the conservative projection
     of `realise()`, `Y += (qa - q) dt/(rho V)`.  The projections are of
     the size of round-off (none significant in any run), and no drift was
     measured although the fresh T5_60 run made 448632 of them in 200
     steps (closure 6.6e-16 at 2 s, mean increment -1.6e-24 +- 1.9e-24 kg
     per step).
   - Documents: E8 (section 2), section 26 item 13, the readme ('Mass
     ledger', 'Round-off of the stiff exchange'), `matrixDefect.H`,
     `phaseChangeLedger.H`, `interfaceExchange.H`.

2. **The long runs after the fix: three run again, eleven revalidated
   (the handover's question).**
   - The 14 recorded runs were made with the binaries before item 1, so
     their keys changed, and their closure is that of the old booking (max
     |closure|/reference: fig9-tev6 7.7e-14, fig9-tev15 1.8e-14, the
     T-Flows set 4.3e-14, T106_6 2.7e-14 and 2.8e-14 (distance and
     faceCells), T106_60 6.7e-14 and 7.5e-14, T106_300 8.2e-14 and 8.2e-14,
     T540_60 2.1e-14 and 2.4e-14, T5_60 7.6e-13 and 7.5e-13, its 5 ms
     variant 5.4e-13).  Running all 14 again
     costs about 207,000 steps, 6 h on 20 cores of this machine at its
     present load (section 25, item 7).  The cheapest honest treatment: the
     three T5_60 runs, which drifted most, are run again with the fix (4
     ranks each, 6000, 6000 and 12000 steps; 2 h 2 min of wall time on 12
     cores at a machine load of 30-43, the two 10 ms runs 1 h 37 min and 1
     h 43 min, the 5 ms run 2 h 1 min), and the other
     eleven are revalidated for the physics criteria.  The backup of the
     recorded T5_60 runs is kept in the scratch area of this round.
   - Revalidation (section 26, item 11) is extended for a change of the
     ledger's booking: the probe compares the balance files without the
     columns solverDefect, clamped and closure (`checks.py sameBalance`;
     every other column, the inventories, released, transport, removed,
     restart, bit for bit), and leaves out `uniform/phaseChangeBalance`
     (the ledger's state: its other entries are those columns).  This is
     justified for this change only by the evidence of item 1 (it acts
     alike in every step and changes no field); revalidation stays opt-in
     (item 5).  The eleven runs were revalidated with `M4LONG_REVALIDATE=
     yes`, each by 5 steps from its last write before its end, the same
     points as in round 1 (fig9-tev6 from 11 s, fig9-tev15 from 29 s, T106
     from 5 s and 13 s, T540_60 from 6.99974 and 6.99952 s, the T-Flows set
     from 29 s): fields, state and profiles identical, the balance
     files identical but in the booking columns.  Their criterion lines
     name them ('revalidated from key ... by a 5-step probe from ... s, not
     run with the current sources').
   - M4.17 (new, M4long) gates the closure of the runs made with the
     current sources, at most 1e-14 of the reference in every step of
     every balance file of a run (a run restarted in segments has
     several; `checks.py closureRuns`).  A revalidated run stands for the
     physics of the current sources, not for their ledger: M4.17 lists it
     with the closure of its booking and does not gate it; it is NOTRUN
     when no run was made with the current sources, FAIL when a run
     failed.  Result: PASS, the three T5_60 runs 8.4e-16, 8.6e-16 and
     8.0e-16; the eleven revalidated runs listed with 1.8e-14 to 8.2e-14
     (not gated).
   - The T5_60 runs again, against the recorded ones (another history of
     decompositions: 4 ranks throughout instead of serial to 34, 26 and 17
     s, then 4 ranks): T_dep 789.2541540570 K against 789.2541540571 K
     (distance), 789.2150310456 K in both (faceCells) and 789.2541542738 K
     in both (5 ms); the gas and wall inventories at 60 s agree to 1.1e-9
     and 6.6e-10 relative at most, the native wall profiles to 5.2e-9
     (L1): the round-off of the other decomposition history, amplified by
     the lagged SuperBee limiter (section 26, item 12), far below what
     acceptance 5 resolves.  M4.5 for T5_60 (+2.25 and +2.22 K, dT_bin 2.77
     K) and M4.5b (789.2542 K at 10 and at 5 ms) pass on the new runs as
     on the recorded ones.

3. **The order of the sums of the nearest-face map (nit;
   `nearestFaceSums()`).**  The reverse distribute used
   `UPstream::defaultCommsType`, nonBlocking (`etc/controlDict`), whose
   combining branch for a contiguous type adds the received slots in the
   order in which the receives complete (`waitSomeRequests`,
   mapDistributeBaseTemplates.C), so a face with the elements of two or
   more other ranks (and its own) could sum them in a timing-dependent
   order, and `mDep_` differ in its last bit between identical runs.  Now
   `UPstream::commsTypes::scheduled`: its own contribution first, then
   those of the other ranks in the order of the map's schedule, the same
   in every run of a decomposition (the same API in v2606).  The output
   is diagnostic; the ledger was never affected.  (Corrected in section 29,
   item 1: mapDistributeBase::schedule() lists directed pairs, and the
   scheduled branch exchanges both ways for each, so between two ranks
   that map elements to each other's faces every contribution was added
   twice; the push is now explicit, in rank order.)  The start-up line
   counts the interface faces with elements of two or more other ranks
   (the same reverse distribute on flags).  New M4.7g: the channel with
   layer distance 5.5e-4 m (three rows per wall), the three rows of the
   lower wall on ranks 0, 1 and 2 (100 faces with the elements of two
   other ranks): two runs on that decomposition write every file, the
   balance, profile and `mDep_` files bit for bit alike; closure 3.8e-16,
   the inventories equal to serial to 5.7e-14 and `mDep_` to 6.8e-13
   (limits 1e-13 and 1e-9).  The test cannot show the
   old defect deterministically (a timing race; the old build may pass by
   chance), so it guards the fix, not the defect.

4. **setFrozenCarrier checks its input against the mesh (minor; two
   reviewers).**  The sign of the inflow through `inlet` was never
   checked: an axis pointing out of the inlet, or `inlet` naming the
   outflow patch, wrote a reversed carrier with exit status 0 (the
   reviewer's wedge: 'discrete inflow through INLET = -2.3148e-08 m3/s');
   with `rescaleToFacets yes` the rescale factor became negative, flipped
   the profile to carry +Q, and U = max(u, 0) axis was 0 in every cell.  A
   `radius` below the wall radius gave backflow in the exact face fluxes
   beyond r = R while U was clipped to 0 there (the reviewer's 2.3 mm on
   the 2.4 mm wedge: 129 of 750 internal faces with phi < 0, one INLET face
   with outflow), and the rescale hid it.  Now the utility stops without
   writing when the analytic inflow through `inlet` is not positive ('The
   flux of the profile through the patch
   ... is ..., not an inflow: the axis ... must point downstream ...'),
   and when the largest distance of a wall point from the axis exceeds
   `radius` by more than 1e-9 relative ('The radius ... is smaller than
   the largest distance ...'), unless `allowBackflow yes`.  New M4.1k
   (the wedge of M4.1c with `rescaleToFacets yes`): a reversed axis,
   `inlet OUTLET` and `radius 2.3e-3` each stop with their message and
   write no field; with `allowBackflow yes` the last writes.  The pipe,
   the wedges and every carrier of M4 and M4long are unchanged (the long
   runs' carriers bit for bit, item 2).

5. **Revalidation opt-in, named, and the binaries built before they are
   kept (minor, nit; `longRuns`).**  Revalidation was on by default, and
   its 5-step probe from a run's last write cannot see a change that acts
   only before it (the release window, the first steps, a fresh start):
   such a change would have been certified with the old results.  It is
   opt-in now (`M4LONG_REVALIDATE=yes`; otherwise a run whose key changed
   is run again, or prepared and NOTRUN with `M4LONG_RUNS=ended`), and
   every criterion line that uses a revalidated run says so.  The binaries
   of the runs were copied from `$FOAM_USER_APPBIN` under the key of the
   current sources without a build (`tests/paperMode/Allrun long` has no
   BUILD step), so stale executables could be stored under a new key;
   `longRuns` now runs an incremental `./Allwmake` first and records the
   md5 sums of the copies, checked whenever they are used.

6. **The pull of the nearest-face map on several ranks (nit).**  Only
   `interfaceTemperature wall` pulls face values to the elements
   (`nearestFaceValues()`), and it ran only in serial (M4.10).  M4.10 has
   a second line now: the HKS inventory channel on the 4-rank split by
   rows, 200 elements with their nearest wall face on another rank, which
   take its temperature through the map: closure 1.9e-16, gas, wall and
   sample equal to serial to 1.3e-15, 2.1e-15 and 1.9e-16 in every step,
   `areaIntegrate` of `mDep_` = the wall inventory to 5.0e-16 (limits
   1e-13, 1e-9 and 1e-12), the same counts, one removal.

7. **The channel helper's description (nit; tests/paperMode/Allrun).**
   'A converged outer loop (nCorr 20, tolerance 1e-12)' was wrong: every
   step of the channels ends at the corrector limit ('outer iterations
   21, NOT converged'; 200 of 200 steps, the last corrector's initial
   residual 6e-8 in the median step, 2.5e-11 to 1.5e-5), the condition of
   M2.4a-c with the lagged SuperBee limiter.  The comments of the helper
   and of M4.7b say so; the gates (serial = 4 ranks to 1e-9) are
   unchanged and pass (7e-14).

8. **The T-Flows advection defaults and A, k per branch (minor; the
   run/014 readme, its fvSchemes, section 26 items 6, 8 and 14).**  The
   GPU solvers' blending defaults differ by branch:
   `Blending_Coefficients_For_Scalars.f90:19` and `..._For_Energy.f90:19`
   read 0 0 1 (pure upwind; `w_upw = 1 - w_cen - w_lin`) on
   advanced_scalar_model and scalar_transport_francesco (HKS, Thesis,
   Transport), but 1 0 0 on multiple_scalars_gpu_francesco (Project),
   whose `Add_Advection_Term.fpp:95-100` weights the central face value
   with the first coefficient: Project runs central differencing for the
   scalars and the energy.  Its control sets only the ignored
   `ADVECTION_SCHEME_*` keywords (`control:45-47`; no GPU source reads
   them) and no blending.  Round 1 said 'pure upwind by default' for every
   branch and 'the energy equation runs upwind everywhere'.  The M6
   profile comparisons get a central variant besides the upwind and the
   blended (0.8 0.1 0.1) ones, and author question E names the defaults
   per branch.  A = 150, k = 0.02 belong to Thesis only (`Types.f90:7-8`);
   Project has 1.35 and 2.7 (`Project:User_Mod/Types.f90:6-7`).
   (Corrected in section 29, item 6: Project's `reference.txt:4-5` lists A
   = 150, k = 0.02 as the parameters of its reference case, while its
   Types.f90 holds 1.35 and 2.7; they are not Thesis-only.)

9. **The carrier against T-Flows' discrete flow (minor; section 25
   items 1 and 5, the readme's carrier rows).**  T-Flows sets its fluxes
   from the parabola at the vertex-mean centres of its cells and faces
   (`HKS:User_Mod/Initialize_Variables.f90:41-46`, `v_flux = sx 0.5
   (u(c1) + u(c2))`, and `:51-59` at the inflow faces; the centres are
   node averages, `Sources/Shared/Grid_Mod/Calculate/Cell_Centers.f90`,
   `Face_Centers.f90`, the boundary cells at the face centres or the
   projection of the cell centre, `Convert_Mod/Calculate_Geometry.f90:
   249-259`, which on this extruded mesh have the same radius).  On the
   INLET of run/001 (540 faces) that rule gives 1.6764077e-6 m3/s =
   1.003837 Q for FLOW_RATE 1.67e-6 (the midpoint rule at the area
   centroids 1.003236 Q; the exact integral of setFrozenCarrier 0.99997
   Q): U_b,d = 0.0931139 m/s against run/014's 0.0927551 m/s, so at the
   same FLOW_RATE run/014's carrier is 0.39 % slower than T-Flows' (the
   reviewer's computation, repeated here).  With T-Flows' fluxes n0/U_b,d
   = 3.3937e-4 mol s/m, within 0.1-0.4 % of Fig. 9 (the reviewer's
   digitisation: n' t_ev = 3.385e-4 at t_ev 15 s, 3.397e-4 at 30 s,
   3.373e-4 at 60 and 90 s), while run/014 gives 3.4068e-4 (+0.5 to +0.8
   %).  (Corrected in section 29, item 6: against those values -0.1 to
   +0.6 %; within 0.4 % of the paper's 3.38e-4.)  So FLOW_RATE 1.67e-6
   with T-Flows' fluxes explains the plateau of Fig. 9; section 25, item
   5's inference 'U_b = 0.0935 m/s (101.5
   mL/min)' is withdrawn.  Acceptance 3 is unaffected ('n0/U_b of the
   carrier used').  A carrier with T-Flows' total flow on this mesh, for
   the profile comparisons of M6: `flowRate 1.67641e-6` with
   `rescaleToFacets yes` (documented, not applied).

10. **T5_60 at 60 s or at completion (minor; section 26 item 5's
    decision).**  The decision cannot move T_dep upstream and does not
    move it in practice:
    - in both carriers (the same T) the first native bin with a layer
      cell below Tdep = 788 K is bin 252 (x = 0.410094-0.411722 m, its
      coldest layer cell 786.92 K), and the coldest layer cell of bins
      0-251 has 789.378 K: the temperature model never deposits there,
      so the 1 % crossing can never lie upstream of the centre of bin 251,
      T_dep <= 789.37 K, at most +2.37 K from the paper's 787 K, within
      the gate of bin 251 (dT_bin 2.77 K);
    - the deposit grows self-similarly: from 30 to 60 s its peak per
      deposited fraction of n0 is constant to 1 % (1.645e-4 to 1.627e-4
      mol/m with layer distance, 1.347e-4 to 1.328e-4 with faceCells), and
      bin 252 holds 0.225 to 0.218 of the peak (faceCells 0.171 to 0.162),
      so at completion the crossing stays between the centres of bins 251
      and 252, at its 60 s position to about 0.01 K;
    - only a concentration of the deposit into one bin, 22 times (16
      times) the present peak while at most twice the present deposit
      remains to come, could push it past the centre of bin 252: at most
      to 784.96 K (-2.04 K, inside bin 253, whose gate is 2.0 K) with
      layer distance.  Nothing in the runs points that way.
    The choice therefore changes the shape of the profile, not T_dep or
    its verdict; the run to completion stays optional (1.5e4 to 3e4 more
    steps, 2 to 5 h on 8 ranks).  `tdepInterval` (checks.py), which lets
    every bin take all the remaining gas, is not sharpened: its bound is
    conservative and decides nothing that the runs need.

11. **The T5_60 offset decomposed (nit; the readme, section 26 item 5).**
    'Its offset comes from the bin resolution' named one cause, against
    the researcher's decision of 2026-10-02.  The +2.25 K have two parts:
    +1.25 K is the bin placement relative to the input Tdep (789.25 - 788
    K), and +1.0 K is the paper's own value lying 1 K below its input (787
    against 788 K), whose cause this case cannot establish; the set-up
    differences (the wall table, the source position, the axial
    resolution, the state compared) are listed without singling one out.
    The wall of the last bin without deposit is 789.37 K (789.36675 K),
    not 789.35 K.

12. **The rear flank of the gas (nit; the readme's domain row, section 26
    item 9).**  By the numbers given there (Fig. 9 at 0.099/0.129/0.160
    m, run/014 at 0.102/0.135/0.167 m) Fig. 9's flank lies 3-7 mm
    UPSTREAM of run/014's, not downstream, and run/014's cloud does not
    'start slightly earlier'; the offset fits the release positions
    (T-Flows into the cells of the inflow faces at x = 0.04 m, run/014
    over the boat, centre 0.045 m).

13. **Two statements (nits).**  (a) Section 25, item 1: the 66 or 70 rows
    of the wall table differ beyond 0.65 m only, but that is not 'far
    downstream of every deposit': by 8 s T540_60 deposits 1.9 % (layer
    distance) and 7.1 % (faceCells) of its inventory beyond 0.65 m (line
    density up to 8.7e-6 mol/m against a peak of 2.9e-5); its T_dep, at
    0.48 m, is not affected.  (b) `Transport` is a scalar-transport test,
    not a third temperature-model case: no `User_Mod`, so no deposition
    source, the scalar enters through the inlet as a flux (`q_01 1.0`,
    `control:88-91`), 12 time steps; it remains a source for the
    advection keyword, the tolerances and the domain.

14. **Three properties of the T-Flows temperature-model codes (nit; for
    author question D; the readme's n0, flow-rate and domain rows,
    section 25 item 4d, section 26 items 9 and 14).**  (1) Both release n
    = 3.17e-5 mol (`Thesis:User_Mod/Source.f90:27`, `Project:User_Mod/
    Source.f90:26`; 14.6 mg/461 g/mol), not the printed 3.16e-5.  (2)
    Their `U_MAX = 0.176` m/s (`Types.f90:6`, 'Fixed bulk velocity') is
    used as the centreline velocity (`Thesis:User_Mod/
    Initialize_Variables.f90:25`): U_b = 0.088 m/s, Q = 95.5 mL/min, and
    the plateau n/(U_b t_ev) = 2.40e-5 mol/m at t_ev 15 s, which Fig. 9
    (2.26e-5) rules out.  (3) Their mesh has 424 axial cells over
    0.04-0.69 m (`Thesis:lesto_pipe.geo:14-17,40,482`: N_STREAM = L 400 +
    1 = 425 nodes with L = 1.06; confirmed by `Source.f90:28`,
    inlet_volume = 2.7600231e-8 m3 = 1.800384e-5 m2 x 0.65/424), dx =
    1.533 mm against run/001's 1.627 mm: another native-bin quantisation
    of T_dep and dT_bin.  No run changes.

15. **The verifier's notes.**
    - The HKS-only items of the section 11 checklist (sigma,
      kineticScale, the vapour-pressure table, the sample's 2/R_s x 5)
      do not apply to run/014; its readme says so (M6 pins them).
    - Not changed, with evidence: 'checks.py integrals ... would pass
      with zero profile files'.  It cannot: `profileFiles()` raises
      SystemExit('No profile files in <case>'), exit status 1, when no
      file matches, also with a set filter that matches nothing (checked
      on an empty case and with the set 'noSuchSet' on the M4 channel).
    - M4.6 needs `LESTO_OPENFOAM_V2606`, and M0.4 (Merlin7) stays NOTRUN:
      unchanged, documented.
    - The pre-round-2 observations (patchDataWave replaced by the exact
      distance, T5_60 at 60 s) are items 1 and 10 above and section 26.

16. **The handover.**  The first fix agent of this round, stopped by a
    usage limit, had made the code of items 1 (in three stages: the flux
    form with the storage diagonal, then the condensate booking, then the
    storage as the inventory measures it, in long double), 3 and 4, the
    tests M4.1k, M4.7g, M4.10 (second line), M4.15, the tightened M2.2h,
    a first M4.17 that gated every run (the recorded runs failed it), the
    longRuns changes of item 5, and the comments of item 7; nothing of it
    had been built into a suite.  Every edit was reviewed: the flux form
    against the v2412 sources (EulerDdtScheme::fvmDdt, gaussConvectionScheme
    and gaussLaplacianScheme with lduMatrix::negSumDiag, fvMatrix::flux,
    the processor patch values), the condensate booking against the
    ledger's measures (it is the CAPPED residue's formula), the scheduled
    branch of mapDistributeBase::distribute (deterministic; v2606 has the
    same API), and the utility's checks; all were kept.  Added in this
    round: the rework of M4.17 and the revalidation of item 2, the T5_60
    runs, `checks.py sameBalance`, the M4.15 header, the excerpts, items
    8-15 and the documents.

17. **Tests and negative evidence.**  New: M4.1k, M4.7g, M4.10 (second
    line), M4.15, M4.17.  Changed: M2.2h (gate 1e-14), M4.11 and the
    excerpts of M2.8, M3.7 and M5.3b (rewritten), the helper comments of
    item 7.  `tests/Alltest` expects them and says so in its usage text;
    `tests/paperMode/checks.py` has `closureRuns` and `sameBalance`.
    Negative evidence: `tests/paperMode/Allrun` with the binaries of the
    round-2 sources (`FOAM_USER_APPBIN` pointed at a copy of them; the
    v2606 compile not set): 32 PASS, 4 FAIL, 1 NOTRUN (M4.6).  The
    failures are the defects and two lines by construction:
    - M4.15: 1.161e-14 after the 20 steps (limit 2e-15);
    - M4.1k: `inlet OUTLET` writes the reversed carrier ('analytic profile:
      -0.9999980686', rescaled to +Q by a negative factor), and `radius
      2.3e-3` writes the carrier with backflow; the reversed axis stops
      only by accident, at the wall table ('-0.00025 lies outside the
      table [0, 0.02]', the wedge's `mode profile`), as the reviewer
      predicted;
    - M4.7g by construction (the old start-up line has no count of the
      faces with two or more other ranks; its two runs happened to agree,
      the race of item 3 is timing-dependent), M4.11 by construction (the
      excerpt of item 1).
    M4.10's second line passes with the old binary too (it guards the
    pull, which worked).  M2.2h with the old binary (the verifier's run of
    this round): 1.642e-13 over its 136364 steps, above the new 1e-14.

M4 status after round 2: `tests/Alltest -fresh all` (M0 to M5 with M4)
gives 221 PASS, 0 FAIL, 1 NOTRUN (M0.4, the compile on v2606 on Merlin7:
no remote access; M4.6 compiles a copy of the sources, setFrozenCarrier
included, against the local v2606 Debug build, every source compiled, 0
warnings) in 48 minutes at a machine load of 30-43, with three 4-rank long
runs alongside; M4 itself, 37 lines, in about 8 minutes; 0 warnings; the
M0 gate identical to main (the logs of 008, 009 and 010, 256, 398 and 419
lines; the 010 excerpt, 402 lines; the 51 files of 010 at writeInterval
3).  `tests/Alltest M4long` (`M4LONG_RUNS=ended`, every run reused) gives
27 PASS, 0 FAIL, 0 NOTRUN in 3 minutes: M4.3 (2 lines), M4.4, M4.5 (10;
T5_60 +2.25 and +2.22 K under the restated gate), M4.5b, M4.16 (2) and
M4.17 (the three T5_60 runs made with the current sources; the eleven
others revalidated and named so in their lines).  The tracked and
untracked files of the repository (git status and the content of all 364
paths) are unchanged by both runs; the branch is phase-change at a43cbd2,
without commits or stashes.  Pending: the researcher's sign-off of
sections 25 to 27, and the questions of section 26, item 14 to the
authors.

---

## 28. Researcher decisions of 2026-10-02 (signed off)

Sections 20 to 24 (M5) are signed off with the decisions below.  Sections
25 to 27 (M4) are not part of this sign-off.

1. **Round-off imbalance of the stiff exchange (section 21, item 1).**
   Accepted as booked: it stays in solverDefect, and E5 and E6 are
   unchanged.
2. **M5.2j with the loose correctors (section 24, item 1).**  The gate '0
   significant projections' is accepted; the guard is not made to act below
   guardTolerance x inventory.
3. **Observable (section 20, item 13).**  wallPlusGas stays the default for
   comparisons at shutdown.  `observable wall` is used at intermediate
   times and where T_dep is defined on the condensed phase (run/014).
4. **M5 acceptance 3 (section 20, item 12).**  Tested by M5.3a (25 steps of
   4 ms with a background Y = 0.1) together with M5.3c (run/013 at its own
   settings): accepted.
5. **The other points of sections 20 to 24** are acknowledged as written.
6. **Wall resistance (section 9, item 3).**  halfCell becomes the default
   for physical runs (M8).  Evidence from the M7 study (Graetz-Robin, Pe 5,
   Bi 1 and 50): halfCell is second order, with an error of at most 6e-5 in
   the developed decay rate at the resolution of the production first
   layer, against 0.3 % (temperature model) to 1.7 % (HKS) for `none`.
7. **Time-step error of the temperature model while gas is present**
   (M7, acceptance 3): reported, not gated.
8. **Final tolerance of HKS cases.**  `Y_<gas>Final` is tightened to 1e-14.
   The M7 study traced an HKS solverDefect of 2e-8 to 6e-8 (above the M8
   gate of 1e-8) to the Final tolerance of 1e-12; 1e-14 cuts it a hundred
   times.  Applied in the M7 integration run (run/013 and its
   `out_excerpt`).
9. **The final-predictor criterion of M5 acceptance 3 (section 14, item
   1).**  '0 significant regime changes in the final predictor' applies to
   the test.  In a whole HKS run it does not hold while the boat evaporates
   (about 230 steps in the M7 study); the guard leaves no significant miss
   there.
10. **M4.5 (Table 2).**  The gate is |T_dep - T_paper| <= max(2 K, dT_bin),
    with dT_bin the drop of the wall temperature across the native bin
    that contains the onset; the cause of T5_60's offset is documented as
    the evidence supports it (sections 26 and 27).

---

## 29. Amendments from the M4 review, round 3 (signed off by the researcher, 2026-10-03)

The M4 verification of round 3 found one major regression of round 2: the
nearest-face map of `layer distance` added the contributions of other
ranks twice wherever two ranks map elements to each other's faces (item
1).  The reviews found seven minor points and nits (items 2 to 8).  Each
was assessed first; all were real.  Item 8, the sector of a wedge, needed
no code change: the option it asks for exists, and it is now tested.  New
criteria: M4.7h (2 lines), M4.18 and M4.19 (M4), and M4.20 (M4long).
Changed: M4.1c (the exact inflow of a wedge), M4.14 (the second line now
gates mDep_), M4.15 (the drift pipe prints the debug report for M4.18),
and the excerpts of M2.8, M3.7, M5.3b and M4.11 (run/012, run/013 and
run/014 rewritten; item 9).  The build of the sources before this pass
(the snapshot `pre-m4fix`) fails every new or changed gate that a change of
the solver addresses (item 11).

1. **mDep_ with nearest-face maps that run both ways (major; section 27,
   item 3; `nearestFaceSums()`, `readDistanceLayer()`).**
   - Defect.  Round 2 replaced the nonBlocking reverse distribute
     (timing-dependent order) by `mapDistributeBase::distribute` with
     `UPstream::commsTypes::scheduled`, `map.schedule()` and `plusEqOp`.
     `mapDistributeBase::schedule()` lists DIRECTED pairs: a rank inserts
     (me, p) when it sends to p and (p, me) when it receives from p
     (v2412 and v2606 `mapDistributeBase.C:186-262`, the same code).  With
     a two-way map both pairs are in the schedule of each rank, and the
     scheduled branch (`mapDistributeBaseTemplates.C:569-650` of v2412)
     runs a full
     send and receive for each entry, so every contribution of the other
     rank was added twice.  `mDep_` exceeded the wall inventory, and the
     start-up count of faces with elements of two or more other ranks took
     a face with the elements of one such rank twice.  Every test of round
     2 had one-way maps (M4.7b, M4.7c, M4.7g, M4.10, M4.14a); the verifier's
     two-way splits showed it.  `mDep_` is output only: the ledger never
     used it.
   - Measured with the snapshot build (verifier's set-ups, rerun here):
     the channel of M4.7a on 4 ranks with rows 0 and 22 on rank 0 and rows
     1 and 23 on rank 1 (the second row of each wall maps to the other
     rank): `areaIntegrate(mDep_)` 50.5 % above the wall inventory, the
     count 200 (true 0); the pipe of M4.7c on hierarchical (1 1 4):
     7.79e-4 above, count 13 (true 0); the tetrahedral pipe of M4.14b on
     simple (1 2 2): 0.925 % above, count 14 (true 0).  The round-1 build
     (nonBlocking) conserved them to 3e-16.
   - Fix (`pushToFaces()` in `interfaceExchange.C`, used by
     `nearestFaceSums()` and by the start-up count).  The sums of the
     elements per remote slot go in `PstreamBuffers` to the rank of the
     face (to proc p the slots of `constructMap[p]`, received there as the
     faces of `subMap[sender]`, in the same order); the receiving rank adds
     its own slots first and then the other ranks in ascending order, every
     contribution once.  The order of the sums is fixed by the
     decomposition alone, so identical runs give identical bits (M4.7g, two
     runs bit for bit, still passes).  The pull (`nearestFaceValues()`,
     `interfaceTemperature wall`) assigns and needed no change.
   - Results: the channel 5.0e-16, the pipe 8.3e-16, the tetrahedral pipe
     1.4e-16 (limit 1e-12), the counts 200/0, 161/0 and 21/0 (elements
     with their nearest face on another rank / faces with elements of two
     or more other ranks).
   - Not affected (shown): with the snapshot and the new build the
     two-way channel and the pipe (1 1 4) write identical fields, per-rank
     state and profiles (128 and 72 files) but for `mDep_PbI2_s` on the
     ranks with two-way faces, the `functionObjectProperties` of the
     `mDep_` integral and `uniform/phaseChangeBalance` (the booking, item
     2), and balance files identical in every column but `solverDefect`
     and `closure` (item 2) and `evaporatedAboveWarning` (item 4).  So the
     ledger, the inventories, the deposits, the profiles and T_dep were
     never affected.  The long runs are decomposed in slabs along the axis
     (simple (N 1 1)): no element has its nearest face on another rank,
     and their revalidation probes write the same `mDep_` with both
     builds (item 10).
   - Documents: readme ('mDep_ through the nearest face'), the
     `interfaceExchange.H` description of the push, section 25 item 2 and
     section 27 item 3 (annotated).

2. **The face-flux correction of the transport matrix in flux form
   (numerics nit and OpenFOAM nit, the same point; `matrixDefect.H`,
   `recordSolve()`).**
   - Defect.  `cellDefect` kept the explicit part of the transport-only
     matrix, `(b - s)_c`, inside r_c and relied on its sum over the cells
     telescoping to 0.  For the solver's matrix it is the non-orthogonal
     correction of the laplacian (the paired species is `fluxRequired`, so
     OpenFOAM keeps it as `faceFluxCorrectionPtr()` and adds `-V fvc::div`
     of it to the source), and its source form telescopes only to
     round-off: `V (sum_f C_f)/V` per cell, plus the rounding of the
     assembly `s + c`.  The closure changed by `dt sum_c (b - s)_c`, booked
     nowhere: on the drift pipe of M4.15 (non-orthogonality up to about 30
     degrees, 20 steps of 0.1 s, 4 ranks) up to 1.7e-24 kg per step with
     layer distance and 6.8e-25 kg with faceCells (the reviewer's
     diagnostic build, the snapshot plus a print of the sum), against an
     explicit source `sum_c |b - s|` of 1.4e-6 to 1.8e-6 kg/s.
   - Fix.  r_c = E_c - (m_c - m0_c)/dt - sum_f (F_f + C_f), with the
     face values C of the correction entered like the coefficient fluxes
     (each internal face once with opposite signs, a processor face with
     the antisymmetric mean of its two sides, a non-coupled face as in
     `fvMatrix::flux()`), and of the source only the explicit remainder
     `E_c = b_c - fl(s_c + c_c)`, `c = -V fvc::div(C)` recomputed with the
     operations of `gaussLaplacianScheme`.  For ddt + div - laplacian with
     an interpolation scheme without explicit correction (upwind, the
     limited schemes such as SuperBee: every case of the solver) the source
     is assembled as exactly `fl(s + c)` (v2412 and v2606 alike:
     `fvMatrix::operator-=`, `gaussLaplacianSchemes.C`), so E = 0 in every
     cell and the closure changes by exactly 0 for the correction; the
     residue (the source form of C against its face values, and the
     rounding of `s + c`) is part of R, i.e. of solverDefect.  A scheme
     with an explicit correction (linearUpwind, LUST) adds it to the
     source only (`gaussConvectionScheme::fvmDiv`), with no face values in
     the matrix: it stays in E, in source form, and its rounding reaches
     the closure, a few eps of the correction per step.  This is
     documented, not changed (no case uses one; recomputing its face
     values would need the Y of the assembly).  Without a correction (an
     `uncorrected` laplacian) E is b - s as before, bit for bit.
   - Debug report: `balance { reportMatrixResidual yes; }` prints, after
     `gSum(fvMatrix::residual())`, `explicit remainder <sum E> kg/s in <n>
     cells (of an explicit source sum |b - s| <S> kg/s)`.
   - Demonstration on the drift pipe, both rules, 20 steps on 4 ranks,
     against the diagnostic build of the snapshot: the new build prints E
     = 0 in 0 cells in all 20 steps of both; its solverDefect differs from
     the snapshot's by the residue the snapshot left unbooked (after 20
     steps -2.6e-24 kg with distance, +2.0e-24 kg with faceCells; per step
     equal to -sum_c (b - s)_c dt within the 2e-25 kg rounding of that
     diagnostic, which rounds every cell to double); the closure columns of
     both builds are identical, max 4.0e-16 of the reference: the residue
     is 100 times below 1 ulp of the 1.24e-6 kg inventory (2.1e-22 kg) and
     cannot be seen in the closure of so few steps.  Over the 6000 steps
     of the T5_60-distance run of M4long, made again with this build (item
     10), the residue amounts to 2.6e-23 kg (0.12 ulp of the 1.45e-6 kg
     inventory): its solverDefect at 60 s differs from that of the round-2
     run by -2.6e-23 kg, its closure in 372 of the 6000 steps by at most 1
     ulp (2.1e-22 kg), max |closure| 8.47e-22 kg in both (8.4e-16 of the
     reference), and every other column of the balance file, every field,
     deposit and `mDep_` at 60 s and all 120 profile files are identical.
   - Tests: M4.18 (new).  The excerpts of run/012, 013 and 014 change in
     the digits of `solver defect` and `solverDefect` only (with item 4).
   - Documents: `matrixDefect.H`, `recordSolve()`, `phaseChangeLedger.H`,
     readme ('Mass ledger'), E8 (section 2, annotated).

3. **The key of the long runs (OpenFOAM minor; `tests/paperMode/longRuns`).**
   - Defect: `longKey` hashed the solver and the utility, the inputs of
     run/014 (without comments), the run's spec and the OpenFOAM
     installation, but not the mesh (the runs link `tests/work/mesh`) nor
     the shell functions that set a run up, so `tests/Alltest M4long` could
     reuse runs made on another mesh or with an older set-up.
   - Fix: the key adds the checksum of the mesh (`tests/work/mesh/mesh.md5`,
     the md5 of run/001's `pipe.msh` that `makeMesh` converted) and the md5
     of `declare -f setupLong case014 makeCase runCaseDir carrier
     decomposeLong decomposeWith edit setControl`, as the cache keys of the
     regression gate do; the header of `longRuns` lists the parts of the
     key.
   - Tests: M4.20 (new, M4long, no run): the key of the first run, with a
     stand-in mesh directory whose checksum differs, and with `setupLong`
     redefined in a subshell, are three different keys.  With the old
     definition all three were equal (05bd528a0e91 for fig9-tev6; the new:
     2b7f7c165071, 1f7e7d0dfb8d and 4104ae678ce9).
   - The new key changes every run's key: the runs were revalidated or run
     again (item 10); none was run again for this reason alone.

4. **The mole-fraction monitor in paper mode (science minor;
   `realise()`).**
   - Defect: x = (Y/M)/(Y/M + 1/M_carrier) counts 1/M_carrier = 250 mol of
     helium per kg of carrier, i.e. rho/M_carrier = 250 mol/m3 at the
     paper mode's rho = 1 kg/m3, while the carrier holds p/(R T) = 10 to 40
     mol/m3 at 1 bar and 305-1210 K: the printed mole fraction was 6 to 25
     times too small (24.7 to 24.9 times in the first steps of run/014, at
     the boat, 1210 K).
   - Fix: x = c/(c + p/(R T)), c = rho Y/M, with p/(R T) of the frozen p
     and T computed once at start-up (`carrierConcentration_`, R the
     solver's 8.314462618).  For a perfect-gas carrier both forms agree;
     the frozen carrier of run/007 (run/012, run/013) deviates from the
     perfect-gas law by up to 0.3 % (its start-up line), and their printed
     maxima fall by 1.4e-4 to 1.6e-4 relative; the cells above 0.05 and
     the mass entered there are unchanged in their excerpts.
   - Consequence for paper mode (a science finding, recorded, no
     change): the gas of the paper's configurations is not dilute.  At
     the end of the release of the recorded Figs. 8-10 run (fig9-tev15, 15
     s, recomputed from its fields) x exceeds 0.05 in 148530 of the 228960
     cells and 0.1 in 54783, at most 0.121 (its log printed 0.00544): the
     plateau of Fig. 9, 2.26e-5 mol/m over the cross-section of 1.81e-5
     m2, is 1.25 mol/m3 against the 10 to 16 mol/m3 of helium at 1 bar
     and 1210 to 750 K along the plateau.  In the T5_60-distance run of
     this build (5 mL/min) x reaches 0.18 at the end of the release (6 s)
     and exceeds 0.05 in up to 56700 cells in a step, and 92 % of the
     released mass (1.33e-6 of 1.45e-6 kg) entered the gas in cells above
     0.05 (the ledger's `evaporatedAboveWarning` at 60 s; 0 in the round-2
     run).  The passive-carrier (dilute) assumption of the paper mode, and
     of T-Flows, holds only to about 10 to 20 % there.  The diagnostic
     `evaporatedAboveWarning` is not part of the closure and changes no
     field, so the revalidation of the long runs leaves that column out of
     the comparison (`REVALIDATE_EXCLUDE_COLUMNS`, documented in
     `longRuns`; in the 13 probes of item 10 it did not differ anyway).
   - Tests: M4.19 (new): the maximum printed for the last step of the
     serial pipe of M4.7c (0.0674589) equals the maximum recomputed from
     its fields to 7e-7 (limit 1e-5); the old form gives 0.00287, 23.5 times
     less.  The excerpts change in the max mole fraction only (with item
     2).
   - Documents: the code comments, readme ('Every step prints'), E11
     (section 2, annotated), the run readmes.

5. **How T-Flows applies the wall table beyond 0.65 m (science minor;
   record; verified in `origin/advanced_scalar_model`, 35c1f6406e).**
   - The count line of `wall_x_profile_complete.dat` says 66 rows, x <=
     0.65 m, of the 70 in the file, and T-Flows reads 66
     (`Gpu:Read_Controls_Mod/Boundary_Conditions.f90:284`).  For a
     boundary condition given in a file, a face gets a value only when its
     centre lies in one of the intervals of the table (`:538-607`); a wall
     face beyond 0.65 m is never assigned and keeps the boundary value 0 of
     its allocation (`Gpu:Var_Mod/Create_Variable.f90:38`, `phi % b =
     0.0`), which the copy of every boundary value into the field (`:670-
     683`) makes its temperature.  The energy equation treats the whole
     WALL region as a Dirichlet wall (`Gpu:Process_Mod/Insert_Energy_Bc.
     f90:49,68`, `Form_Energy_Matrix.fpp:189-199`): in LESTO_HKS the wall
     beyond x = 0.65 m is at 0 K.  The conductivity lookup
     (`HKS:User_Mod/Beginning_Of_Iteration.fpp:33-36`, `dataIndex =
     int((T - 250)/20) + 1`, clamped above only) extrapolates k_He below
     250 K and reads outside the table below 230 K for the cells that wall
     cools.
   - Neither `rows count` with `outOfRange hold` (the wall beyond 0.65 m at
     T(0.65 m) = 324 K; section 25, item 1 called it what T-Flows reads)
     nor `rows all` (the measured 315 to 305 K at 0.66-0.69 m; run/014)
     reproduces this, and no option of setFrozenCarrier does: the record
     is corrected (section 25 item 1 and section 11 annotated; the
     setFrozenCarrier header and its dictionary comment; the run/014
     readme's wall-table row).  It matters beyond 0.65 m only (T540_60
     deposits 1.9 % and 7.1 % of its inventory there by 8 s with the
     measured wall; a colder wall would deposit more), and it belongs to the
     parity of the HKS code-to-code case (M6), where it is a question for
     the authors.

6. **run/014 numbers and the parity record (science nits).**
   - The deposit of run/014 as distributed begins 2.6 cm (25.9 mm)
     upstream of the published one: the comment of
     `constant/thermochemistryProperties` said 2.5 cm.
   - Courant number 0.11 for the case as distributed (100 mL/min: 2 U_b
     dt/dx = 0.113), 0.12 at the 106 mL/min of Table 2: the readme and the
     controlDict comment said 0.12.
   - The plateau of T-Flows' carrier, n0/U_b,d = 3.3937e-4 mol s/m, lies
     within 0.4 % of the paper's 3.38e-4 and within -0.1 to +0.6 % of the
     plateaus of Fig. 9 as digitised (3.373e-4 to 3.397e-4 for t_ev 15 to
     90 s); the readme said 'within 0.1 %' and section 27, item 9 'within
     0.1-0.4 %' (both corrected).
   - The correctors per step, per run family of M4long (the readme said '2
     correctors in most steps'): Fig. 9 at t_ev 6 s 2 in 89 % of the
     steps, at t_ev 15 s 1, 2 or 3 (37, 38 and 24 %); T106 2 in 67-78 %, 3
     in 9-27 %; T540 (0.22 ms) 1 in 60 %, 2 in 35 %; T5_60 (10 ms) 3 or 4
     (45-56 % and 28-41 %), at 5 ms 2 in 80 %; the T-Flows set (4 ms) 1 in
     34 %, 3 in 30 %, 5 or more in 24 %.  Every step converged except the
     first of the T5_60 runs and of the T-Flows set (11 correctors, the
     limit).
   - A = 150, k = 0.02 are not Thesis-only: Project's `reference.txt:4-5`
     lists them as the "Parameter used in the reference case", while its
     `User_Mod/Types.f90:6-7` holds 1.35 and 2.7 (both of commit
     49806abdd8, "update with different reaction rate models"; the rate
     A[1 - exp(k (T - Tdep))] is the same form, with 1.35 and 2.7 nearly a
     step at Tdep).  Recorded in section 27 item 8 and author question D
     (section 26, item 14; annotated) and in the run/014 readme.
   - The HKS sample (section 11 and the M6 notes, annotated): T-Flows
     zeroes the sample at the start of every iteration of step t_ev/dt
     (`HKS:User_Mod/Beginning_Of_Iteration.fpp:78-94`, `Curr_Dt dt ==
     T_EVAPORATION`; step n spans [(n - 1) dt, n dt]: `Sources/Process/
     Cpu/Time_Mod/Needs_More_Steps.f90`, which the GPU build links, and
     `HKS:Main_Pro_Frozen.f90:149-171`), so the sample
     evaporates over [0, t_ev - dt].  The solver removes an inventory
     sample at the step that starts at t >= removeTime - dt/2: `removeTime`
     t_ev - dt reproduces T-Flows (t_ev itself removes it one step later).

7. **The comment blocks of `interfaceExchange.H` (OpenFOAM nit).**  The
   readme said that its description 'is therefore split into two comment
   blocks'; it has five, each below 16 KB (the largest 11.9 KB): 'split
   into several comment blocks (each below 16 KB)'.

8. **The sector of a wedge (verification nit; my call: no code change, a
   test).**  setFrozenCarrier takes the sector from the normals of the two
   wedge planes, pi - acos(n1 . n2), exact only to their rounding, unless
   `sectorAngle` [deg] is given: that option exists (section 25, item 1;
   the dictionary documents it).  Measured from the written `phi` with
   exact sums: with `rescaleToFacets yes` the inflow of the 5-degree wedge
   is Q 5/360 to 2.5e-13 from the wedge patches and to 5.3e-17 with
   `sectorAngle 5`.  M4.1c now checks both (limits 1e-12 and 1e-15), and
   the header says to give `sectorAngle` where the inflow must be exact to
   round-off.

9. **The excerpts.**  run/012, run/013 and run/014 were run for the 5
   steps of their excerpts with this build, as M2.8, M5.3b and M4.11 do,
   and compared token by token (123, 111 and 126 lines): they differ only
   in the digits of `solver defect` (5, 2 and 1 tokens) and of
   `solverDefect` with its relative value (item 2), and in the max mole
   fraction (5 tokens each; item 4); the closures, the counts of cells
   above 0.05 and the masses entered there are unchanged.  The three
   excerpts are rewritten from those logs (M2.8, M3.7, M5.3b, M4.11).

10. **The long runs: one run again, thirteen revalidated (cost).**
    - Items 1, 2 and 4 change the solver and item 3 the key, so every run's
      key changed.  Item 2 changes the booking only and item 4 a
      diagnostic only; item 1 changes `mDep_`, which the probes compare.
      The thirteen runs other than T5_60-distance were revalidated
      (`M4LONG_REVALIDATE=yes`), each by 5 steps from its last write
      before its end with the builds of round 2 and of this pass, the same
      points as in rounds 1 and 2 (fig9-tev6 from 11 s on 4 ranks,
      fig9-tev15 and the T-Flows set from 29 s on 8 ranks, T106_6 from 5 s
      in serial, T106_60 and T106_300 from 13 s on 4 ranks, T106_60
      faceCells in serial, T5_60-faceCells and the 5 ms variant from 59 s
      on 4 ranks, T540_60 from 6.99974 s on 8 ranks and 6.99952 s on 4
      ranks): 9, 36 or 72 files per probe identical, `mDep_` included, the
      profile files identical, and the balance files identical in every
      column but `solverDefect` (in 11 probes; clamped, closure and
      `evaporatedAboveWarning` identical in all 13).  The 13 probes took 2
      minutes; with the build and the gate of `tests/Alltest M4long`
      before them 6.
    - M4.17 gates the closure of the runs made with the current sources;
      after the revalidation none would be, and item 2 changes the booking
      that M4.17 gates, so T5_60-distance (the run that drifted most
      before round 2) was run again with this build: 6000 steps of 10 ms
      on 4 ranks, 39 minutes of wall time at a machine load of 12-16 (the
      record of round 2 is kept in the scratch area of this pass).  The
      tool that started the run stopped it at its time limit of 30 minutes,
      at 32.37 s; `longRuns` resumed it from its write at 32 s (an exact
      binary restart: the restart column stays 0), so it has two segments
      (24 and 15 minutes, 3237 and 2800 steps, 36 steps done twice).  It is
      the only run made again; the whole set would cost about 6 h on 20
      cores.
    - Results (`tests/Alltest M4long`, every run reused): 18 criteria
      PASS, among them M4.17 with T5_60-distance at 8.4e-16 over its 6036
      steps (largest at 2.09 s), the thirteen revalidated runs listed with
      the closure of their booking (not gated) and named so in their
      lines; M4.5 for T5_60 with layer distance +2.25 K (T_dep 789.2542 K,
      equal to the round-2 run in all 17 printed digits: 789.25415405700539
      K), dT_bin 2.77 K, met; M4.5b 789.2542 K at 10 and at 5 ms; M4.20
      PASS (key 2b7f7c165071, with another mesh checksum 3b83e58c2397, with
      a changed setupLong 4104ae678ce9).

11. **Tests and negative evidence.**
    - New: M4.7h (2 lines), M4.18, M4.19, M4.20 (M4long).  Changed: M4.1c
      (the exact inflow, both sector paths), M4.14 (the tetrahedral pipe:
      mDep_ in serial and on 4 ranks), M4.15 (the debug report on, the
      gate unchanged), M4.11 and the excerpts of M2.8, M3.7 and M5.3b.
      `tests/Alltest` expects them and says so in its usage text;
      `tests/paperMode/checks.py` has `mdepCount`, `explicitRemainder`,
      `moleFraction` and `inflowExact`; `longRuns` has M4.20 and the
      extended key.
    - Negative evidence: `tests/paperMode/Allrun` with the build of the
      snapshot (`FOAM_USER_APPBIN` pointed at it; the v2606 compile not
      set): 34 PASS, 6 FAIL, 1 NOTRUN (M4.6).  The failures are the defects
      and one line by construction: M4.7h, the channel (the count 200
      instead of 0; `mDep_` per face 50 % from serial, its integral 50.5 %
      above the wall inventory) and the pipe on hierarchical (1 1 4) (13
      instead of 0; +7.79e-4); M4.14, the tetrahedral pipe on 4 ranks (14
      instead of 0; +0.925 %; in serial 4.9e-16); M4.18 (no explicit
      remainder printed, in either layer rule or the positive control: the
      report does not exist); M4.19 (0.00287 printed against 0.0675); and
      M4.11 by construction (the excerpt rewritten with the digits of items
      2 and 4).  The criteria that passed include M4.7g and M4.15 (they
      guard round 2's fixes) and M4.1c (item 8 is a test, not a fix).
    - The key of M4.20 without the mesh and the set-up functions gives
      three equal keys (item 3).

M4 status after round 3: `tests/Alltest -fresh all` (M0 to M5 with M4)
gives 225 PASS, 0 FAIL, 1 NOTRUN (M0.4, the compile on v2606 on Merlin7:
no remote access; M4.6 compiles a copy of the sources, setFrozenCarrier
included, against the local v2606 Debug build, every source compiled, 0
warnings) in 30 minutes at a machine load of 15-16, with the 4-rank long
run alongside; M4 itself, 41 lines, in about 8 minutes; 0 warnings; the M0
gate identical to main (the logs of 008, 009 and 010, 256, 398 and 419
lines; the 010 excerpt, 402 lines; the 51 files of 010 at writeInterval
3).  `tests/Alltest M4long` (`M4LONG_RUNS=ended`, every run reused after
the revalidation and the run of item 10) gives 28 PASS, 0 FAIL, 0 NOTRUN
in 2 minutes: M4.3 (2 lines), M4.4, M4.5 (10), M4.5b, M4.16 (2), M4.17 and
M4.20, plus the gate.  The tracked and untracked files of the repository
(git status and the content of all 364 paths) are unchanged by both runs;
the branch is phase-change at a43cbd2, without commits or stashes.
Pending: the researcher's sign-off of sections 25 to 27 and 29, and the
questions to the authors of section 26, item 14, with the additions of
items 4 to 6 above (the 0 K wall beyond 0.65 m of LESTO_HKS, the
parameters of Project's reference case, the dilute assumption).

---

## 30. Researcher decisions of 2026-10-03 (signed off)

Sections 25 to 27 and 29 (M4) are signed off as recommended, with these
points.

1. **`layer distance`** is the exact distance of the cell centre from the
   nearest face of all interface patches on all ranks (section 26, item
   1), not patchDataWave as section 2, section 12 B3 and deliverable 2
   said.
2. **setFrozenCarrier** solves the steady T-Flows energy equation with the
   wall table as a Dirichlet condition (`mode solved`, linearUpwind;
   section 25, item 1).
3. **The carrier's flow.**  No rescale to the faceted inlet by default
   (T-Flows' choice).  T-Flows' own discrete flow is 1.0038 Q on this mesh,
   0.39 % above run/014's at the same FLOW_RATE: a point for M6 (section
   27, item 9).
4. **The wall table beyond 0.65 m.**  run/014 reads all 70 rows; T-Flows
   applies 0 K there, which no option reproduces (section 29, item 5).
5. **Profiles of the distance layer** bin the deposit of a layer cell by
   its own axial extent (section 25, item 2).
6. **A_I and V_I** are compensated sums with `areaPerVolume layer` and
   plain sums with `cell` (section 25, item 2).
7. **run/014** is the configuration of Figs. 8-10 as the paper states it
   (Tdep 680 K, t_ev 15 s).  It places the onset 25.9 mm upstream of the
   published one, while the T-Flows set (Tdep 638.15 K, A 150, k 0.02,
   distance 2.4e-4) reproduces the figures (section 26, item 6): question
   D to the authors.
8. **The Table 2 runs** use Table 1's T_dep^calc per case and the paper's
   flow rates as non-expanding volumetric flows (section 25, item 4).
9. **Acceptance 3-5** as defined in section 25, item 5: A3 with the
   discrete U_b of the carrier; A4 in L1; A5 at the end of the experiment,
   with the irreversibility bound for runs stopped early, 'for both layer
   rules' read as each rule must match, under the gate of section 28, item
   10.
10. **M4.7c** tests the parallel agreement with tight correctors; the
    production SuperBee loop's difference of about 1e-8 on hierarchical
    splits is a property of the outer tolerance (section 25, item 6;
    section 16, item 7).
11. **run/014 keeps SuperBee** (paper Section 2.2).  M6 compares its
    profiles with upwind, blended (0.8 0.1 0.1) and central variants as
    well (section 26, item 8).
12. **The long runs** are reused by opt-in revalidation; M4.17 gates the
    closure only of runs made with the current sources (sections 26 item
    11, 27 item 2, 29 item 10).
13. **E8.**  The transport defect is booked in flux form, the face-flux
    correction too, and the rounding of the stored condensate as `clamped`
    (sections 27 item 1, 29 item 2).  Schemes with an explicit correction
    in source form (linearUpwind, LUST) leave that correction's rounding in
    the closure; no case uses them.
14. **Carried into the M7 integration run:** a criterion that discriminates
    the face-flux correction (M4.18 checks only the debug report), the
    gated start-up count of the tetrahedral pipe (M4.14), and the stale
    comments and wordings of the last review.
15. **For the authors and the paper:** the questions of section 26, item
    14, and the findings of section 29, items 4-6 (the paper-mode gas is
    not dilute: PbI2 mole fractions up to 0.12 in the configuration of
    Figs. 8-10 and 0.18 in T5_60).

---

## 31. Rebase onto Bojan's main (signed off by the researcher, 2026-10-03)

Signed off by the researcher on 2026-10-03: decisions 1 to 6 and the points
of items 7 and 8, as recorded below.

**The new base.**  The main branch moved from a43cbd2, the base of
sections 0 to 30, to 71bf192, nine commits of Bojan Ničeno:

- cosmetic changes (comments, indentation);
- case `run/011-coupled-species-sources` with a conservative gas-to-solid
  mock: `evaluateThermochemistry()` at the beginning of the step returns
  the precipitation PbI2_g -> PbI2_s (proportional to Y_PbI2_g),
  `solveGasSpecies.H` treats the gas sink implicitly and passes the
  realised sink to the solid, a per-species balance
  (`collectSpeciesBalance.H`, `checkSpeciesConservation.H`), integrated
  source lines, and the source fields `source_<species>` written
  (AUTO_WRITE); tagged `0.1` = 7155e6b;
- `applications/xgemsFoamProbe`, a small application that calls xGEMS
  (64c1c51);
- the xGEMS flags in `Make/options`, unconditional, from
  `$(HOME)/Development/GEMS-Related/install` (077ca12), and a call of
  xGEMS on a cement demo system at every start-up of the solver (71bf192).

The researcher's goal: take his commits and rebuild the phase-change work
on them as if it had started from his base; the result is the project as
it was on `phase-change-m0-m5` plus the integration below.  Every file
that only the phase-change work changed is taken as it is; the four files
both sides changed are merged by meaning: `rhoFixedFlowFoam.C` (4
conflicting hunks), `readme.md` (1), `Make/options` and
`solveSolidSpecies.H` (merged by git without conflict; the first then
linked xGEMS unconditionally, the second is his file with the header note
of section 5).

**Decisions.**

1. **Model mock is his coupled mock, bit for bit.**  With
   `constant/thermochemistryProperties` absent or `model mock` the solver
   behaves exactly like the pristine build of the tag 0.1: his loop order
   (the thermochemistry at the beginning of the step; all sources
   suppressed outside the cells next to WALL when solid species exist; the
   gas solved with the sources by his `solveGasSpecies.H`; the solid source
   set to the realised gas sink; the integrated-source lines; the solid
   solve; `checkSpeciesConservation.H`), his AUTO_WRITE source fields and
   his log lines.  This amends section 1 (mock path: "the existing code
   path, no new Info lines, `evaluateThermochemistry()` and
   `solveSolidSpecies.H` unchanged"): the existing code path and its Info
   lines are now those of 0.1, and `evaluateThermochemistry.{H,C}`,
   `solveGasSpecies.H`, `collectSpeciesBalance.H` and
   `checkSpeciesConservation.H` are his files, unchanged (`git diff
   origin/main` of each is empty, so section 17, item 8, "byte-identical to
   main", holds for the new main); `solveSolidSpecies.H` is his file plus
   the header note that it serves model mock only.  In the loop,
   `phaseChange.beginStep()` (nothing in model mock) comes before the
   thermochemistry; the WALL mask stays model mock only (section 17, item
   9).  The one known difference from 0.1 is unchanged since section 15,
   item 7: a decomposed run whose carrier has uniform processor patches
   differs in the last digits, for every model, mock included.
2. **The phase-change models keep their structure and results.**  A paired
   gas is solved by `solvePairedGasSpecies.H` with the exchange; the
   ledger and the profiles are unchanged.  `evaluateThermochemistry()` is
   not called outside model mock, so every source field stays zero; the
   integrated-source lines and the per-species balance
   (`checkSpeciesConservation.H`) are skipped, and the source fields are
   set NO_WRITE after `createPhaseChange.H`, so a phase-change model
   writes and prints exactly what it did before.  An unpaired gas (e.g. a
   tracer) is solved by his `solveGasSpecies.H` with a zero source, which
   leaves it unchanged to the last bit: `fvm::Sp` of a zero field adds
   +-0 to the diagonal, `max(source, 0)` adds 0 to the source,
   `setFluxRequired` stores the face-flux correction without changing a
   value, and `collectSpeciesBalance.H` only fills lists that are not
   reported outside model mock.  Proof: M0.1c (the gas solves of 008 and
   009 reproduce their out_excerpts, recorded before 0.1), M0.1d (every
   Y_tracer line of 010's out_excerpt), M2.9 (run/012 with the tracer of
   run/010: its lines and field equal the pristine build's, the balance
   file of PbI2_g equals run/012's without the tracer), and a one-off
   comparison with the build of `phase-change-m0-m5` on that case (5
   steps: the logs identical, 176 filtered lines; all 13 written files and
   the balance file cmp-identical).
3. **xGEMS is opt-in; his code is kept.**  The environment variable
   `LESTO_XGEMS` (an xGEMS installation prefix) adds his include and
   library flags in `Make/options` (`-I<prefix>/include`,
   `-L<prefix>/lib -lxGEMS -Wl,-rpath,<prefix>/lib`); `Make/xgemsConfig.sh`
   defines
   `LESTO_HAVE_XGEMS` in the generated header `xgemsConfig.H`, as
   `gemsConfig.sh` does for the bridge, so that a plain wmake recompiles
   `rhoFixedFlowFoam.C` after a toggle; `Allwmake` keeps the stamp
   `xgems.state`, runs wclean on a toggle and resolves a relative path.
   His includes and his start-up block in `rhoFixedFlowFoam.C` are
   compiled only with `LESTO_HAVE_XGEMS`, unchanged (his demo path
   included) and with a comment; `-help` reports the state.  Without it
   the solver builds (0 warnings) and runs here.
   `applications/xgemsFoamProbe` stays his separate application, not built
   by our Allwmake (the readme says how to build it).  None of his code or
   files is deleted.  Tests M0.3j-l use a dummy xGEMS installation.
   One-off: with the conda build of xGEMS 2.1.2 (prefix `~/opt/gems-env`,
   built with gcc 16) and his demo path pointed at the local xGEMS
   sources in a scratch copy, the v2412 build (gcc 7.5) compiled, linked
   and ran his call (`xGEMS status: 2, converged: 1, iterations: 157`, the
   rest of the log of case 008 unchanged), with 129 `-Wold-style-cast`
   warnings from the Eigen headers of that prefix.
4. **Case numbering.**  His `run/011-coupled-species-sources` keeps 011;
   ours are renamed (mv): `011-temperature-model` -> `012`,
   `012-hks-inventory` -> `013`, `013-paper-temperature-model` -> `014`.
   Every reference was shifted in the order 013 -> 014, 012 -> 013, 011
   -> 012 (never his 011): the tests (including the names of their work
   directories, e.g. `phaseChange/run012`, and of their variables), the
   case readmes and a dictionary comment, the case paths in the
   out_excerpts of 012-014 (only their `Case` and `Using dynamicCode`
   lines, which regress.sh replaces by `<case>`), the solver readme, the
   source comments, this plan (all sections; its history is otherwise
   unchanged) and the long-run machinery.  The manuscript holds no case
   reference.  Section 5's future cases shift by one: 014/015 paper mode,
   016+ physical carriers.
5. **The regression gate against the tag 0.1.**  The pristine reference is
   the tag 0.1 (`PRISTINE_REF`, pinned; before: the branch main), because
   71bf192 needs his local xGEMS installation, and without xGEMS it equals
   0.1 (`git diff --stat 0.1 origin/main`: `Make/options`, the xGEMS
   start-up block of `rhoFixedFlowFoam.C` and the marker 'NEW' of its first
   comment line, `applications/xgemsFoamProbe`).  His case 011 joins the
   gate: its log is identical to 0.1's and the first block of its
   out_excerpt is reproduced.  Criteria whose expectations encoded the old
   mock, restated:
   - M0.1a: 008-011 (was 008-010).
   - M0.1b: 011's out_excerpt, 519 lines (was 010's, 402 lines: 0.1's
     coupled mock changes the PbI2 lines of 010, whose out_excerpt the
     main branch did not record again).
   - M0.1c: the out_excerpts of 008 and 009 are compared without the
     per-step report of the coupled mock (`regress.sh -noMockReport`); the
     banner line remains the only difference.
   - M0.1d (new): 010's out_excerpt without that report and without the
     lines of Y_PbI2_g and Y_PbI2_s (`-ignore`): 150 lines identical for
     both builds, the start-up block and every Y_tracer line.
   - M0.2: 18 fields per write time (the 3 source fields added), 60 files
     (was 15 and 51).
   - M1.2a: the write directories of `writeFrozenFields no` hold the 3
     source fields too (case 010 is model mock); M1.2b: at most 14 MB and
     a third of `yes` (was 8 MB and a quarter); M1.2c: the source fields
     compared as well, 21 files (was 12).
   - M2.7: model mock against 0.1 for 008-011 and 60 files (was 008-010,
     51 files).
   - "No new Info lines" in model mock (section 1) now means the Info lines
     of 0.1; M0.1a and M2.7 compare the whole logs.
   - New: M2.9 (an unpaired gas in a phase-change model) and M0.3j-l (the
     xGEMS opt-in).  CACHE and the dynamicCode seed are keyed on the
     inputs of 008-011; SELF expects the 5 gate ids.
6. **Documentation.**  The solver readme merges his sections and ours:
   model mock is described as his coupled mock, the limitations of both
   are merged, Building describes the xGEMS opt-in and
   `applications/xgemsFoamProbe`, and Tests the gate against 0.1.  This
   section records the decisions; the case numbers of the plan are
   updated; the rest of its history is not rewritten.
7. **Quality (2026-10-03).**
   - Build: `wclean` and `Allwmake` on v2412 with 0 warnings (BUILD); a
     copy compiles and links against the local v2606 source build with 0
     warnings (M4.6).  M0.4 (v2606 on Merlin7) is NOTRUN, as before.
   - `tests/Alltest -fresh all` against the pristine tag 0.1: 231 passed,
     0 failed, 1 not run (M0.4) in 26 minutes at a load of 4-7 (before:
     225 passed; new: M0.1a of 011, M0.1d, M0.3j-l, M2.9).
   - `tests/Alltest M4long` with `M4LONG_REVALIDATE=yes M4LONG_RUNS=ended`:
     the 14 recorded runs, made with the binaries of `phase-change-m0-m5`
     (key 1172f9ad) and the case path run/013, are revalidated by their
     5-step probes (the set-up made now from run/014 equals each run's
     dictionaries without comments and its carrier bit for bit; both
     binaries write the same files): 29 passed, 0 failed, 1 not run, in 5
     minutes.  M4.17 is NOTRUN: no run was made with the current sources
     (before the rebase T5_60-distance was); the closures of all 14 are
     printed, not gated (8.0e-16 to 8.2e-14).  To make M4.17 a gate again,
     one run must be made with the current sources: the cheapest is
     T5_60-distance, about 0.7 h on 4 ranks (6037 steps; fig9-tev15-TFlows
     about 0.7 h on 8 ranks).  It was not run (no long rerun unless
     unavoidable).  One-off evidence that the ledger is unchanged: probes
     of T5_60-distance (4 ranks), T106_6-faceCells (serial) and
     fig9-tev15-TFlows (8 ranks), 5 steps from their last write with the
     old binaries 1172f9ad and with the new ones, write identical files
     with nothing left out but `uniform/time` (the ledger state
     `uniform/phaseChangeBalance` and `mw_PbI2_s` included), identical
     complete balance files (all columns, `solverDefect`, `clamped` and
     `closure` included), identical profiles and identical logs.
   - The input link `constant/wall_x_profile_complete.dat` of the recorded
     runs (absolute, to run/013-paper-temperature-model, read only by
     setFrozenCarrier at set-up) was repointed to run/014, the same file.
   - The tests leave the repository's files unchanged (work area
     `tests/work/`, git-ignored).

8. **The last points of the review (2026-10-03).**
   - `tests/functions` unsets `LESTO_GEMSBRIDGE` and `LESTO_XGEMS`, so an
     opt-in exported in the user's shell no longer leaks into the build
     under test; the build tests set both explicitly where they test them.
   - Comments and wording: the `writeFrozenFields` note of
     `rhoFixedFlowFoam.C` (model mock also writes `source_*`),
     `createPhaseChange.H` (the logs of 008-011 and the fields of 010 at
     write interval 3 are reproduced bit for bit), and the readme's claim
     on an unpaired gas (the bitwise statement rests on the one-off
     comparison of item 2; M0.1c, M0.1d and M2.9 check the printed lines
     and the equivalence with 0.1).
   - M4.17 is gated again: T5_60-distance was run again with the
     integrated build (6000 steps, 18 min on 8 ranks; its revalidated
     record is kept in the work area as
     `T5_60-distance.revalidated-backup`).  Its closure stays below
     7.9e-16 of the reference in every step (limit 1e-14), and T_dep is
     789.2542 K, as before.  `tests/Alltest M4long` gives 30 PASS, 0 FAIL,
     0 NOTRUN; the other 13 runs are revalidated, with the wall-table links
     of run/014 resolving.
