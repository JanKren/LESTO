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
  (Implemented in M9, section 32: every condensed DC is suppressed; `frozen`
  evaluates p_eq at the nodes of the pair's table in [max(minTemperature,
  T_min), T_max] of the GEMS grid, cold, at the reference composition, and
  interpolates that GEMS table (section 12, item 7); `local` also fails a
  call whose element balance is off by more than 1e-6 and rejects a result
  further than maxLog10Deviation from the table; chi = n_g/c, the fraction
  of the formula units present as the gas species, clipped to [1e-12, 1].)

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
        // as implemented (section 32): defaults are the values above but
        // for system, mode, carrier and (frozen) referenceComposition;
        // new: referencePressure 101325 (frozen [Pa]),
        // writeEquilibriumField no (pEq_<gas>, gemsSource_<gas>, chi_<gas>);
        // the pair's gems block takes an optional excess { I 1e-6; }
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
| `tests/` | `testPhaseChange.C` and `testVapourPressureTable.C` (g++ only); `regress.sh` (filtered full log vs pristine build and vs `out_excerpt`); `cmpFields.sh` (binary cmp of written fields); `phaseChange/{closedBox,channel,pipe,restartAscii,parallel,graetzWedge,layerMesh}` with an `Allrun` that asserts the milestone thresholds; `graetz.py`, `gci.py`, `summarise.py`. (As built: the verification of M7 is `tests/verification/`, section 33.) |
| `utilities/setFrozenCarrier/` (new app) | Paper-mode carrier: rho = 1, parabolic u from Q/A with an optional facet-area rescale, phi exactly divergence-free, T from a wall table. Not solver logic. |
| `thermochemistry/systems/make_pv_table.py` + `data/pv_PbI2_phases.csv`, `data/pv_PbI2_TFlows.csv` | Direct NASA-9 (Gurvich): per-phase cr/l columns every 5 K (270-1250 K) plus the T_m = 683 K node; T-Flows GEMS_VAPOR_PRESSURE at its 95 nodes (283.15 + 10k K, bar). |
| `run/012-...` to `run/017-...` (new) | 012 temperature model (007 carrier, dev); 013 HKS inventory (dev); 014/015 paper-mode temperature/HKS (code-to-code); 016+ physical carriers (M8). Each has a readme and an `out_excerpt` baseline. 008-010 untouched; 011 is the coupled-mock case of the main branch (section 31). (Amended in section 32, item 12: 015 is the GEMS case of M9, `run/015-hks-gems-local`; the paper-mode HKS case of M6 becomes 016, the physical carriers 017+.) |
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
(Implemented and tested by `tests/Alltest M7`, criteria M7.1a-o, part of
`all`, and `tests/Alltest M7long`, criteria M7.2a-f, M7.3a-d, M7.5a-c,
M7.3pa-pc and M7.6a-b, not part of `all`; sections 33 and 35.  The closed
box of the deliverables is verified by the analytic closed boxes of M2.1
and M5.1; section 35, item R.)
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
(Implemented and tested by `tests/Alltest M9`, criteria M9.0 to M9.14;
sections 32 and 35.)
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
     with "not yet implemented" (M9).  (Implemented in M9, section 32:
     with the bridge `equilibrium GEMS` runs, without it the rebuild hint
     stays; `speciation lagged` needs equilibrium GEMS in mode local.)
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

---

## 32. M9 implementation decisions (signed off by the researcher, 2026-10-04)

M9 (the GEMS3K backend of the HKS law) is implemented on top of the
integrated code of section 31 and tested by `tests/Alltest M9`
(`tests/gems/Allrun`, criteria M9.0 to M9.10, and M9.11 to M9.14 since
section 35; part of `all`).  The points
below were decided during the implementation; each sharpens or deviates
from a statement above (E9; sections 1, 4 and 5; section 6 M9; section 9,
items 7-8; section 12, item 7) and needs the researcher's sign-off.
Section 33 is reserved for the M7 integration that follows.

1. **What a call computes (E9).**  One call per element serves every pair
   of the element.  Every condensed species of the system is suppressed
   (metastability bounds 0 and 0; also Pb(cr) and Pb(l), which belong to no
   pair: the HKS law exchanges only the pair's condensate), so GEMS3K
   returns the homogeneous gas of the bulk composition, and the dual
   activity of a suppressed pure condensate is its log10 saturation ratio:
   p_eq = p_g 10^(-max_phi log10 Omega_phi) over the pair's condensates.
   For the ideal PbI2-He gas p_eq depends on T only and equals the table's
   minimum over the phases to the G0 interpolation of GEMS3K on the 10 K
   grid of its DCH file.  Scan of the bridge (`scan3`, cold calls, 500-1250
   K, x = 1e-4 to 0.2): |dlog10 p| <= 6e-5; warm calls in strongly
   supersaturated cells (T < 700 K, x >= 0.01) drift to 1e-3 decades,
   about proportional to log10 Omega (after 1-2 iterations from the warm
   state), well inside maxLog10Deviation; the largest deviation in all M9
   runs is 3.8e-4 decades.

2. **The dictionary (section 4).**  `GEMSCoeffs` as section 4, with the
   values listed there as defaults for every key but `system`, `mode`,
   `carrier` and, in mode frozen, `referenceComposition` (all required);
   the start-up line prints the values used, and a key of the other mode
   is named in a warning (section 35, items G and H: such a key is not
   read either, and an unknown key stops the run).  New keys:
   `referencePressure` (mode frozen, the pressure of its calls, default
   101325 Pa) and `writeEquilibriumField` (item 7).  The pair's
   `gems { gas; condensates; elements; }` block is checked against the
   system: the gas species must be a gas, every
   condensate a condensed species, and all of them must have the formula
   of `elements` (the mass balance of a pair assumes one formula); the
   molar mass of the gas species must equal the solver's `molarMass` to
   1e-6.  `equilibrium GEMS` without `GEMSCoeffs` or without a pair's
   `gems` block, and `speciation lagged` without equilibrium GEMS in mode
   local, stop with a FatalIOError (M9.7, 9 refusals; 16 lines since
   section 35, items H and I).  With model HKS and
   equilibrium table, a `GEMSCoeffs` or `gems` block in the file is named
   in a start-up warning (as the inputs of the other model, section 21,
   item 6).

3. **`excess` (a new key of the gems block).**  At the exact formula PbI2
   (b_I = 2 b_Pb) the interior-point method of GEMS3K occasionally stalls
   (thermochemistry/README.md, pitfall 3).  Scans of the bridge
   (scratchpad of this round, `scan.cpp`, `scan2.cpp`): in a temperature
   sweep (1250 to 270 K, x = 1e-6 to 0.3) 1 of 2364 warm calls failed in the
   guarded range (x >= 1e-6, T >= 500 K) at the exact formula and none with
   a relative iodine excess of 1e-6; outside the guards (x down to 1e-10,
   T down to 270 K) 369 of 3940 failed at the exact formula and 1 with the
   excess.  The prototypes and the frozen `referenceComposition` of
   section 4 (I = 2.000002e-3) carry the same 1e-6.  `excess { I 1e-6; }`
   adds nu_I 1e-6 c iodine to the bulk: p_eq of an ideal gas does not
   change, the speciation by about 1e-6.  An explicit key, default none
   (the exact formula), not a hidden constant; run/015 and the tests set
   it.

4. **mode frozen (section 12, item 7).**  Every rank creates an engine; the
   master evaluates p_eq at the nodes of the pair's table that lie in
   [max(minTemperature, T_min), T_max] of the GEMS3K grid (for the Gurvich
   table 152 nodes in 500-1250 K: every 5 K and the melting point
   683.000002 K), with cold calls at the reference composition and
   pressure; a node whose call fails, or whose p_eq lies further than
   maxLog10Deviation from the table, keeps the table's value (counted in
   the start-up line).  The values form a table of one column, p_eq,
   interpolated as the pair's table (the melting point is a node, so the
   kink stays at a node), broadcast to every rank; the engines are then
   destroyed.  An element below minTemperature or outside that range takes
   the pair's table: the guards and the fallback of mode local.  Results:
   152 nodes from GEMS3K, 0 failed, 0 rejected, max |dlog10 p| 5.7e-5 at
   the nodes (section 35, item E: GEMSB_ERR_FATAL now re-creates the
   engine here too); on the 1150 -> 450 K channel the 212 elements in 500-1100 K
   lie within 5.6e-5 decades of the 5 K table (acceptance 2, limit 1e-3;
   M9.2a), the 14 below 500 K on the table exactly.  Frozen is
   deterministic: a restart is cmp-identical, and on 4 ranks the p_eq of
   every element equals the serial one bit for bit (M9.2b).

5. **mode local.**  One engine per rank for the whole run (RAII; the
   engine of the bridge is created and destroyed with its owner).  At the
   first step of a run and then every updateInterval steps, in
   `beginStep()` before the correctors, every element (wall and sample; a
   sample element for its own pair) is evaluated at its T_e, the pressure
   of its cell and the composition of its cell at the start of the step:
   the carrier p/(R T) (the frozen p and T, as the mole-fraction monitor
   of section 29, item 4) times its formula, plus c_k = rho Y_k/M_k of
   every pair times its formula (with its excess); a step thus uses one
   law throughout its correctors, the guard and `realise()`, and the
   ledger books the realised flows as before.  The guards in this order:
   T_e < minTemperature, T_e or p outside the grid of the system, x_k =
   c_k/(c_k + p/(R T)) < minMoleFraction (no gas: c <= 0), a failed call (a
   status other than OK and OK_RETRIED, or an element balance off by more
   than 1e-6 relative: the README's 'check both the status and the mass
   balance'), a rejected result (|dlog10| > maxLog10Deviation); each gives
   the table and is counted per kind.  Every element keeps its own warm
   state (GEMSB_WARM_CELL); a failed or rejected element's state is
   invalidated.  GEMSB_ERR_FATAL re-creates the engine: destroyed, created
   from the system file, every condensed species suppressed again, every
   warm state of the rank invalidated (section 12, item 7; M9.5 with a
   fault-injecting bridge).  The warm states are not written: a restart
   begins cold (section 8, risk 13); the restart from 0.2 of the ramp
   channel is exact at its start and stays within 5.6e-8 of the
   continuous run's inventories at 0.35 s (M9.9), while frozen restarts are
   exact.  (Section 35, item C: the schedule of the updates follows the
   time index, and the lagged state of the last update is written with
   the per-rank state, so a restart with updateInterval > 1 repeats the
   continuous run up to its next update.)

6. **speciation lagged (E9).**  chi_e = n_g/c_k, the fraction of the pair's
   formula units present as the gas species g in the GEMS3K state of the
   last update (PbI2(g); the rest PbI, I, Pb, I2), clipped to [1e-12, 1]
   (E9 wrote chi = p_i^GEMS/p_i: the same in the dilute limit; GEMS3K holds
   the total pressure at the cell's p, while the solver's p_i = c R T comes
   on top of the carrier's p, so the two differ by the factor 1 + x, and
   the fraction of molecules is the quantity the law needs),
   multiplies the partial pressure of the law, q = A G' (p_eq - chi beta
   Y): g_e is multiplied by chi_e, and with wallResistance halfCell the
   half cell enters with chi beta (G' = G/(1 + G chi beta d/(rhoD)), the
   Robin form of the law with the speciated partial pressure).  An element
   whose p_eq comes from the table has chi = 1.  mode local only (frozen
   has no composition; refused).  The separate run of acceptance 4 (M9.4,
   reported, not gated): on the ramp channel chi = 0.46 to 1 over the GEMS
   elements; at 0.35 s the wall holds +7.6e-5 and the gas +4.0e-6 of the
   inventory more than without speciation, the deposit profile differs by
   1.2e-3 (L1), T_dep (fromInlet, wall) 900.15 against 900.13 K.  Scan of
   chi at 1 atm (cold calls): 1 up to 700 K; at 1100 K 0.97 (x = 0.01),
   0.91 (x = 1e-3), 0.36 (x = 1e-5); at 1200 K down to 0.06 (x = 1e-5):
   the speciation matters only in the hot, dilute gas.

7. **Output.**  The step report of mode local has, after the exchange line,
   one line per pair: the update (calls: warm, retried, cold; iterations
   per call; failures; engine re-creations) or 'no update this step', the
   elements of each source, and the largest |dlog10 p_eq/p_eq,table| of
   the GEMS elements (and with speciation the range of chi).  It holds no
   timing, so excerpts compare (run/015).  `end()` prints the summary with
   the totals and the mean times of a warm and of a cold call (all ranks;
   64-bit counts since section 35, item A).
   Every rank writes `<logDirectory>/gemsEquilibrium.log` below its case
   directory (processorN in parallel; section 35, item D: also a directory
   outside it is per rank): the engine, every update with its
   counts and timing on the rank, every failure (element, T, p, x, status,
   the message of GEMS3K; at most 20 per update), every re-creation;
   GEMS3K's `ipmlog.txt` lies next to it (the log directory and level of
   GEMS3K are set before the first engine).  With `writeEquilibriumField
   yes`: `pEq_<gas>`, `gemsSource_<gas>` (0 GEMS3K, 1 below minTemperature,
   2 below minMoleFraction, 3 outside the GEMS grid or table, 4 rejected, 5
   failed; -1 elsewhere or before the first update) and, with speciation,
   `chi_<gas>`, per element on the interface faces and per cell (the mean;
   the largest source), refreshed in `balance()`, never read.  The fields
   let the tests recompute every guard (M9.3a) and every frozen value
   (M9.2a) independently.

8. **allowNonPhysicalCarrier (section 4).**  The carrier line already gives
   max |T - p W/(R rho)|/T, which is max |rho - p W/(R T)|/rho; above 0.01
   mode local stops with a FatalIOError (paper mode, rho = 1: 0.95 on the
   channel) unless `allowNonPhysicalCarrier yes` (a warning; the composition
   of every cell comes from p/(R T)); mode frozen runs on any carrier
   (M9.3e).

9. **Cost (acceptance 3: <= 70 us per warm call).**  A warm call costs
   about 20-30 us when consecutive calls share their temperature and about
   35 us more when they do not, because GEMS3K re-interpolates G0(T) at
   every new T (the bridge clamps Ttol to 1e-3 K against the Ttol trap of
   the README); every element of the solver lies at its own T, so that is
   the solver's case (`scan4`, load 17: 62 us against 28 us at one T).  On
   the idle workstation the scan gave 42-45 us per call in that pattern;
   the solver measured 65.7 and 60.9 us per warm call on the ramp channel
   (56687 warm calls, 1.05 iterations per call) at a load of 17.4-17.6 (the
   M7 study on 16-20 of the 20 physical cores, hyperthreads shared), and
   73 us in the first 25 steps of run/015.  M9.3b gates the channel's
   value, prints the load, and is NOTRUN, not FAIL, when it is slower on a
   busy machine (a 1-minute load of 4 or more), as acceptance 3 asks for a
   quiet moment.  Mode local adds about 0.05 s to a run/015 step of 2.3 s
   in its first steps (700 calls).  Not done: ordering the calls by T (it
   would create long runs of nearly equal temperatures, the trap the clamp
   guards against), or keeping one engine per temperature.

10. **The Ttol trap and the element order.**  Consecutive calls whose T
    differ by less than 1e-3 K reuse G0 of the previous call (the clamp of
    the bridge), and the reference temperature moves with every call, so a
    long run of slightly increasing temperatures would accumulate an error.
    The elements are called in their order (by cell), where neighbours
    differ by kelvins or by round-off; any drift beyond maxLog10Deviation
    would be rejected and counted.  The largest deviation in all M9 runs is
    3.8e-4 decades (item 1).

11. **The default build stays bridge-free (deliverable 5).**  The M9 tests
    build a copy of the sources with the bridge into
    `tests/work/<platform>/gems/bin` (keyed by the sources, the bridge and
    the OpenFOAM installation; M9.0 checks every source compiled, 0
    warnings, ldd and -help, and that the default solver has no bridge).
    `tests/functions` unsets LESTO_GEMSBRIDGE; the bridge is
    LESTO_GEMSBRIDGE_TEST or the repository's `thermochemistry/gemsbridge`
    (its `libgemsbridge.so` is git-ignored: without one every criterion but
    M9.0 and M9.1 is NOTRUN).  The dummy bridge of the build tests (M0.3)
    is now a copy of the real header with `tests/build/dummyBridge.c`, which
    defines all 30 functions: the solver calls the whole API since M9, so
    the old one-function dummy no longer linked.  `makeCase` links a case's
    `constant/gems` like its tables.

12. **Case numbering (section 5).**  Section 5 had 015 for the paper-mode
    HKS case of M6 and 016+ for the physical carriers; this round's task
    gives 015 to the GEMS case, `run/015-hks-gems-local`, so M6's case
    becomes 016 and the physical carriers 017+ (section 5 annotated).

13. **run/015-hks-gems-local.**  run/013's set-up (the 007 carrier, the
    boat, faceCells, SI, removal at 6 s) with `equilibrium GEMS`, mode
    local, the guards of section 4, `excess { I 1e-6; }`, `constant/gems` a
    link to `thermochemistry/systems/PbI2He`, and the final solve at
    1e-14, the final tolerance of HKS cases of section 28, item 8 (the
    solver defect of 25 steps 2.2e-16 instead of 3e-14; about 20 % more
    time per step); run/013 itself keeps 1e-12 until the M7 integration.
    Its `out_excerpt` holds 5 steps and, after a `...` line that ends the
    compared block, the end-of-run summary, which holds timing (M9.8).
    `**/gemsLog/` is git-ignored.

14. **Tests (`tests/gems/Allrun`, `checks.py`, `faultBridge.c`).**  The
    channel of M5 with the boat (1000 -> 400 K, a perfect-gas carrier) and
    its own loose correctors with the Final solve at 1e-14: equilibrium
    table and GEMS local for 350 steps; local on 4 ranks, with speciation,
    written every step (30 steps), with updateInterval 3, with the
    fault-injecting bridge (faults at the calls 150 and 600), restarted from
    0.2, and with wallResistance halfCell (the default planned for physical
    runs, section 28, item 6; with and without speciation) against the table
    run with halfCell (M9.10: the half cell is stored per element and the
    law recomputed with it at every update; without the half cell the
    deposit differs by 4.6 % (L1), 15 times the gate, so a lost half cell
    fails); the 1150 -> 450 K channel in mode frozen (serial, restarted, 4
    ranks); paper mode; 9 refusals; run/015 for 25 steps.  Results of this
    round (load 17, the M7 study alongside): M9.3a 0 failures in 350 steps,
    the sources of all 240 elements in 30 steps equal to the guards
    recomputed from the fields (0 differences; 34 wall elements below 500 K
    in every step); M9.3c 4 engines, 4 logs, the calls of the rank logs =
    the summary's, inventories within 1e-9 of serial; M9.3d the deposit
    within 6.5e-5 (L1) and 1.7e-4 of the peak (bin) of the table run at
    every profile time (limit 3e-3); closure <= 1.9e-16 and solverDefect <=
    2.5e-13 in every run (corrected in section 35, item B: 2.6e-13 on 4
    ranks, a run whose defect was not gated then; every run gates both
    now).  `faultBridge.c` forwards every call to the real
    bridge (dlopen) and is loaded through LD_LIBRARY_PATH over the RUNPATH
    of the bridge build, so no second build is needed.

15. **Negative evidence.**  `tests/gems/Allrun` with a bridge build of
    the sources with five defects put in on purpose (`LESTO_GEMS_TEST_EXE`
    names the executable to test; the build step is then skipped and the
    first line of M9.0 is NOTRUN): (a) the engine re-created after
    GEMSB_ERR_FATAL without suppressing the condensed species again, (b)
    without invalidating the warm states, (c) the mole-fraction guard
    checked before the temperature guard, (d) mode frozen without the
    fallback below minTemperature, (e) paper mode accepted in mode local.
    Results: M9.5 FAIL ('fault 1: 0 species suppressed on the new engine,
    not 4; fault 1: 81 elements called warm after the re-creation', the
    same for fault 2), M9.3a FAIL twice (1020 element-steps whose source
    differs from the guards recomputed from the fields, while the log and
    the fields agree with each other; 0 elements below minTemperature in
    350 steps), M9.2a FAIL (14 elements below 500 K with the source 'outside
    the GEMS table' instead of 'below minTemperature'), M9.3e FAIL (the run
    on rho = 1 did not stop), M9.8 FAIL by construction (the excerpt's GEMS
    lines count the guards differently); every other line PASS (M9.0 second
    line, M9.1, M9.2b, M9.3b, M9.3c, M9.3d, M9.4, M9.6, the 9 refusals of
    M9.7, the second line of M9.8, M9.9: the defects do not act there;
    M9.10 was added after this run).
    One-off (no criterion): the bridge build compiles and links against
    the local v2606 Debug build with 0 warnings (every source compiled), and
    its binary reproduces the balance file of 20 steps of the ramp channel
    in mode local bit for bit (1836 GEMS3K calls, 0 failures).

16. **Outside the files of this round.**  `thermochemistry/README.md`
    still says that the bridge 'is not yet wired into the solver'; it is
    not among the files implementers edit (section 13).  (Updated by the
    cleanup, section 35, item L.)  M10 (Pb-Bi-I)
    will need several pairs in one call (supported: one call serves every
    pair of an element) and the qualification of the guards at trace
    composition (section 8, risk 7).

M9 status (written by the cleanup of section 35, as section 34, item 2
asked, and brought up to date by its follow-up): M9 is implemented with the
decisions above, signed off in section 34, and with the amendments of
section 35, items A to M and V to Y (pending the sign-off): 64-bit counts;
the closure and the solver defect of every run gated, and a ledger with a
value that is not finite refused; an update schedule that follows the
absolute time index, with the lagged state as restart state; a log
directory per rank, decided alike on every rank and named in the start-up
line; the engine re-created in mode frozen too; unusable results of
converged calls counted as failures and logged with their values; only the
keys of the mode read and unknown keys refused; run/015's excerpt compared
exactly with the bridge that wrote it and, with another bridge, in its
start-up and first step (the rest NOTRUN); a load sampler, and its test,
that leave no process behind.  `tests/Alltest M9` records 44 lines,
criteria M9.0 to M9.14 (two lines each of M9.0, M9.3a, M9.3b, M9.5 and
M9.11, three of M9.6 and M9.8, sixteen of M9.7, one of the others); in the
final `tests/Alltest -fresh all` of section 35 every line passed (M9.3b:
64.2 us per warm call at a 1-minute load of 18.8; in the first round 42.9
us at 2.5).  The default build stays bridge-free (M9.0); the bridge build
compiles with 0 warnings (M9.0; v2606 one-off, item 15).  Open: M10
(Pb-Bi-I: several pairs in one call, the guards at trace composition and
in a gas that is not dilute, section 35, item M), and a bridge for Merlin7
(section 35, item L; with it the first line of M9.8 compares run/015's
start-up and first step only, item V).

---

## 33. M7 integration (signed off by the researcher, 2026-10-04)

M7 (section 6) is integrated on top of M9 (section 32): the verification
suite of the two isolated M7 studies (the first of 2026-10-01/02:
Graetz-Robin, the time-step study on the 007 carrier and the refined pipe
meshes; the second of 2026-10-03: the long runs on the paper-mode carrier)
moved into `tests/verification/`, together with the decisions of section
28 (items 6 and 8), the carry-overs of section 30 (item 14) and the
corrections that the checker of the first study found in its reporting.
`tests/Alltest M7` (part of `all`) runs the Graetz-Robin series, criteria
M7.1a-m (and M7.1n-o since section 35); `tests/Alltest M7long` (not part
of `all`; M7.6a-b since section 35) runs the refined meshes
(M7.2a-b), the time-step study on the 007 carrier (M7.3a-d), the h vs h/2
study (M7.2c-f), the axial x2 study (M7.5a-c) and the time-step study on
the paper-mode carrier (M7.3pa-pc).  The points below were decided in the
integration; each needs the researcher's sign-off.

1. **The suite and its split.**
   - Files (`tests/verification/`): `graetzRobin`, `graetz.py`,
     `graetzChecks.py`, `graetzCase.py`, `graetzWedge/` and the test helper
     `graetzCarrier/` (M7); `pipeMeshes`, `meshChecks.py`, the test helper
     `pipeMeshFacts/`, `longStudies`, `longChecks.py` and `dtChecks.py`
     (M7long).  The staged files were adapted, not applied: the case
     numbers of section 31 (the first study's run/011 and run/012 are
     run/012 and run/013), the Alltest of M9 (merged by hand), and the
     corrections below.  The staged plan texts (numbered 25, 26 and 32 by
     the studies) are this section; their decisions (a)-(e) of the first
     study are those of section 28.
   - `all` holds only the Graetz-Robin series, about 20 min on 4 cores.
     Everything that runs on the pipe is in M7long (the second study's
     decision f): the first study had put its time-step study on the 007
     carrier into M7, which would have made `all` about 3 h.
   - One driver for the long runs.  The first study's `dtStudy` is folded
     into `longStudies` as a second family of runs, the edits of run/012
     and run/013 on the 007 carrier (only `deltaT`, the writes every
     second and a profiles dictionary with the native and the 1 cm bins
     and profile times every 0.5 s differ from the cases).  So every long
     run of M7 is keyed, kept with its binaries, reused, resumed and
     revalidated by M4long's method (sections 26 and 27): the key holds
     the sources and the inputs of the run's own case (run/012, 013 or
     014, the linked carrier and tables included) without comments, the
     md5 of its mesh, the set-up functions, its spec and the OpenFOAM
     installation.  `dtChecks.py` keeps the measures of the time-step
     studies (study, history, plot); its `runs` check is replaced by
     `longChecks.py runs`, which reads segmented runs.
   - Style: `pipeMeshFacts.C` is rewritten in the solver's style (a
     PURPOSE header, 2-space indent, braces on the line of their
     statement); `graetzCarrier.C` already was.
   - The closed box of the deliverables (section 6) is verified before M7,
     by the analytic closed boxes of M2.1 (the temperature model: Y^n =
     Y0/(1 + lambda dt)^n to 1e-14 and closure <= 1e-15, also with halfCell,
     interfaceTemperature wall, areaMultiplier 2, depositionVelocity and a
     relaxation factor) and M5.1 (HKS at 700 K: the bare wall to 1e-14,
     CAPPED with Y_s = 0 exactly, the ample condensate); M7 adds no
     criterion of its own for it (added by section 35, item R).

2. **graetzCarrier stays a test helper** (instead of a wedge option of
   `setFrozenCarrier`).  On the Graetz wedges (vertices at (x, r, -+r
   tan(theta/2))) u must be the parabola of the wedge's radius y: the
   cells then have exactly the metrics of an axisymmetric scheme.
   `setFrozenCarrier` evaluates the parabola at the distance from the
   axis, sqrt(y^2 + z^2) = y/cos(theta/2) on the wedge planes: on these
   wedges that is the axisymmetric problem with the wall at R
   cos(theta/2), a fixed bias of the decay rate of +6.9e-5 (Bi 1) and
   +2.3e-4 (Bi 50) at 2 degrees, larger than the GCI of halfCell (1e-6);
   and its wall points lie at R/cos(1 deg), beyond the radius, which it
   refuses unless `allowBackflow`.  A second velocity definition in the
   utility of the T-Flows carrier was not worth it; the helper (230 lines)
   is built into the work area, as `makeCarrier` of the M2 tests.

3. **Acceptance 1, the Graetz-Robin series (M7.1a-m).**  The set-up of
   the first study, unchanged (`graetzRobin` header): the extended Graetz
   problem with the temperature law as a constant first-order uptake
   (k (Tdep - T) = 100, so f = 1 exactly), R = 2.4 mm, T = 500 K, p = 1
   atm, D = 7.6963e-5 m2/s, Pe 5, Bi 1 and 50, L = 24 R; the reference by
   three methods (series, shooting, Chebyshev; bound 1.7e-12 and 4.5e-14);
   wedges of 2 degrees, Nr 32/64/128/256 with dx = 2 dr; GAMG to 1e-13;
   steady states by implicit Euler at dt = 1e6 s (Co 2e9 to 1.7e10).
   Results of the integration run (25 PASS in 15.5 min on 4 cores, with
   16 cores of M7long alongside); every one of the 44 fitted mu equals the
   first study's bit for bit (its build of the M5 review, round 4, before
   M4): the transport and the wall laws of the temperature model are
   unchanged since.  mu_ref = 107.593854287 (Bi 1) and 273.173651541 1/m
   (Bi 50); zeroGradient outlet:

   | scheme | Bi | law | Nr | mu_fine [1/m] | p | mu_ext/mu_ref - 1 | GCI_fine |
   |---|---|---|---|---|---|---|---|
   | linear | 1 | halfCell | 64/128/256 | 107.594017088 | 2.0010 | +1.07e-9 | 1.9e-6 (0.00019 %) |
   | linear | 1 | none | 64/128/256 | 107.728384521 | 1.0077 | +8.84e-6 | 1.55e-3 (0.155 %) |
   | linear | 50 | halfCell | 64/128/256 | 273.173434995 | 1.9919 | +9.47e-9 | 1.0e-6 (0.0001 %) |
   | linear | 50 | none | 64/128/256 | 273.923038744 | 1.0045 | +1.10e-5 | 3.41e-3 (0.341 %) |
   | linear | 1 | halfCell | 32/64/128 | 107.594505599 | 2.0033 | +1.40e-8 | 7.6e-6 |
   | linear | 1 | none | 32/64/128 | 107.863396655 | 1.0154 | +3.50e-5 | 3.08e-3 (0.308 %) |
   | linear | 50 | halfCell | 32/64/128 | 273.17278248 | 1.9718 | +6.97e-8 | 4.1e-6 |
   | linear | 50 | none | 32/64/128 | 274.674046307 | 1.0081 | +3.88e-5 | 6.78e-3 (0.678 %) |
   | SuperBee | 1 | halfCell | 64/128/256 | 107.594218826 | 2.0024 | +4.68e-9 | 4.2e-6 |
   | SuperBee | 1 | none | 64/128/256 | 107.728586999 | 1.0143 | +1.63e-5 | 1.54e-3 (0.154 %) |
   | SuperBee | 50 | halfCell | 64/128/256 | 273.176516478 | 2.0104 | +6.44e-8 | 1.30e-5 |
   | SuperBee | 50 | none | 64/128/256 | 273.926144881 | 1.0231 | +5.77e-5 | 3.36e-3 (0.336 %) |

   The fixedValue outlet gives the same table to the digits shown (M7.1g:
   mu within 3.5e-10).  M7.1a: the reference bound 1.69e-12 (limit
   1e-10).  M7.1b: continuity of the carrier <= 1e-15 of the inflow on all
   four meshes.  M7.1d: the window's profile changes by 1.1e-7 to 2.2e-7
   between the first two steps and by exactly 0 after them.  M7.1e: every
   mu_fit equals mu_h of the documented scheme to 7.1e-9 (Bi 50; the GAMG
   tolerance at the window's small values), SuperBee to 6.2e-9 (M7.1l).
   - Acceptance 1 holds for both laws, both Bi and both convection
     schemes on Nr 64/128/256: GCI_fine at most 0.341 %; the
     extrapolated mu within its allowance (halfCell 1e-9 to 6e-8 of the
     reference, none 9e-6 to 6e-5); observed orders 2.00, 1.99 (halfCell;
     SuperBee 2.00, 2.01) and 1.008, 1.005 (none; SuperBee 1.014, 1.023).
     On 32/64/128 none at Bi 50 has GCI_fine 0.68 %, above 0.5 %, as the
     discrete predictor said before the first study's runs (none is first
     order with an error of 0.70 dr/R): that mesh is too coarse for the
     first-order law, not a defect (in the table, not a gate).  Negative
     evidence: `GRAETZ_MAIN="32 64 128" graetzRobin check` makes M7.1h
     fail for exactly that line (and pass the 3 others).
   - The corrections of the first study's reporting (its checker), made
     here:
     a. SuperBee's steady state.  The first study wrote 'then exactly 0
        (SuperBee <= 2.2e-13)'.  In fact SuperBee's intermediate steps
        change the window's profile by up to 1.1e-12 (Nr 128, Bi 1, steps 2
        to 4, none and halfCell alike) and is 0 after them there; at Nr 256,
        Bi 1, it changes by 8e-14 to 3.2e-13 in every step, and its last
        step, which M7.1l tests (limit 1e-12), by at most 1.8e-13; at Bi
        50 it is 0 after the first two steps.
     b. The time-accurate control.  It was described as 'written every
        0.05 s up to 1.5 s', but it wrote every 107 steps (0.05005 s), its
        last write was at 1.4514 s and its final state (1.50009 s) was not
        written.  Its steps are now rounded up to whole write intervals:
        3210 steps of dx/u_max = 4.678e-4 s, 30 writes, the last at 1.5015
        s with the final state.
     c. The halfCell rate of M7.1c.  The checker found the agreement 'to
        4.0e-14' only with d from OpenFOAM's centroids of the wall row
        (3.7e-14), and 1.3e-13 with 'the exact analytic centroid'.  That
        1.3e-13 is the centroid formula (2/3)(r2^3 - r1^3)/(r2^2 - r1^2)
        evaluated in metres, whose difference of cubes loses digits: with
        the centroid in exact arithmetic (mpmath) the rate agrees to
        3.4e-14, with the formula on the unit radius times R (graetz.py,
        the study's definition) to 4.0e-14, with OpenFOAM's centroids to
        3.7e-14.  M7.1c states its definition (the unit radius; limit
        1e-12) and prints the agreement with OpenFOAM's centroid too.
     d. The insensitivity to the fit window (M7.1f).  The figures of the
        first study (1.1e-9, 7.5e-9 with the fixedValue outlet) moved one
        end of the window at a time; both ends moved give up to 2.2e-9 and
        1.2e-8 (the checker).  The 8 variants are fitted now (one end or
        both, -+1.5 R each) and the largest change of each kind printed:
        7.5e-9 (one end) and 1.2e-8 (both ends), both with the fixedValue
        outlet at Bi 1, 7.5e-11 (mixing cup) and 8.1e-10 (wall uptake), both
        at Bi 50 (halfCell, Nr 256, with the fixedValue outlet; with
        zeroGradient 7.52e-11 and 8.11e-10); limit 1e-7.  (Corrected in
        section 35, item R: 'all with the fixedValue outlet at Bi 1', and
        item Z: 'zeroGradient outlet'.)
     e. M7.1k compares the time-accurate run with mu_h, the decay rate of
        the discrete mode of its mesh, to 1e-9, not with the
        pseudo-transient run to 1e-8.  The first study had set 1e-9 a
        priori and raised it to 1e-8 after the first run measured 4.5e-9,
        which the checker traced to the pseudo-transient run: the GAMG
        tolerance 1e-13 at the small values of the Bi 50 window puts its
        fit +4.5e-9 above mu_h, while the time-accurate run is within 7e-11
        of it.  Measured now: +7.3e-11 (time-accurate; limit 1e-9) and
        +4.5e-9 (pseudo-transient, printed).

4. **Acceptance 2, h vs h/2 (M7.2c-f).**  The six runs of the second
   study (revalidated, item 7): at 30 s, flux form with none, the onset
   moves by -0.002 native bins (+0.005 K, dT_bin 2.44 K) and the largest
   mDep_ by -0.024 % (limits 1 bin and 5 %); the volumetric form misses by
   its peak (-32 %, onset +0.03 bins), as expected; halfCell moves the
   onset by -0.002 bins and the peak by +0.006 %.  Decisions:
   - 'T_dep shift <= 1 bin' is the shift of the interpolated 1 % onset
     position in native bins (the width of the bin of the onset on h);
     the shift in K is printed against dT_bin (second study, b).
   - M7.2f (halfCell) is a gate with the limits of acceptance 2, as M7.2d:
     halfCell is the wall law of the physical runs (item 9; second study,
     c).
   - M7.2e (the volumetric form) is reported: PASS when computed, its line
     says whether it misses the limits (it does), as section 6 asks
     ('reported, and expected to fail').

5. **The axial x2 study (M7.5a-c): the criterion of the second study is
   a gate** (`M7LONG_AXIAL_GATE` default yes; `=no` reports it): on the
   axial x2 mesh, flux form, none and halfCell, at 30 s, T_dep shifts by
   at most 1 native bin, the 1 cm peak height by at most 5 % and the L1
   change of the 1 cm wall profiles is at most 1 %.  Measured: 0.21 bins,
   7e-6 and 9.8e-4 (none), 0.21 bins, 6e-6 and 9.8e-4 (halfCell): margins
   of 5, 7000 and 10.  Section 6 lists the study as a deliverable without
   a criterion; a gate keeps it from regressing silently, with the limits
   of acceptance 2 (the onset to its resolution, the peak) and 3 (the L1
   norm).

6. **Acceptance 3 (M7.3a-d on the 007 carrier, M7.3pa-pc on the
   paper-mode carrier).**
   - The 007 carrier: run/012 (temperature model) and run/013 (HKS, now
     with the final tolerance 1e-14, item 8) at dt 4, 2 and 1 ms to 16 s,
     made with the current sources (item 7), at 16 s:

     | model | bins | L1(4, 2) | L1(2, 1) | order | Richardson at 1 ms | peak 2 -> 1 | centroid 2 -> 1 | T_dep 2 -> 1 |
     |---|---|---|---|---|---|---|---|---|
     | temperature | native | 5.11e-7 | 2.74e-7 | 0.90 | 3.2e-7 | -2.4e-7, same bin | -4.3e-6 mm | +3.5e-7 K |
     | temperature | 1 cm | 4.41e-7 | 2.22e-7 | 0.99 | 2.2e-7 | -4.9e-8, same bin | -4.4e-6 mm | +1.5e-6 K |
     | HKS | native | 8.89e-4 | 4.51e-4 | 0.98 | 4.6e-4 | -2.2e-5, same bin | -1.5e-2 mm | +8.5e-4 K |
     | HKS | 1 cm | 8.69e-4 | 4.41e-4 | 0.98 | 4.5e-4 | -6.0e-5, same bin | -1.5e-2 mm | +1.7e-3 K |

     Acceptance 3 holds for both models (limit 1 %; M7.3b, M7.3c).  The
     gas of the temperature model has gone by 16 s (2e-14 to 4e-14 of the
     release left); its deposit still depends on dt at first order
     (observed 0.90 and 0.99), but weakly, L1(2 ms, 1 ms) 2.7e-7 (native)
     and 2.2e-7 (1 cm), 1600 times less than that of HKS, which holds 3.6 %
     of its inventory in the gas at 16 s and whose deposit carries the
     first-order error of the transient (4.5e-4, order 0.98).  (Corrected
     in section 35, item R: 'so its deposit no longer depends on dt'.)
     The closure is at most 4.4e-16 of the reference in every step of the
     six runs (gate 1e-14), the solverDefect at the end 5e-13 to 4e-12
     (temperature) and 1.4e-10 to 5.4e-10 (HKS; gate 1e-8, M8's) of the
     reference.

   - The paper-mode carrier: run/014 as distributed (100.2 mL/min, layer
     distance), dt 4, 2 and 1 ms to 30 s, the runs of the second study
     (revalidated): L1(2 ms, 1 ms) 1.29e-7 (native) and 7.14e-8 (1 cm),
     observed order 0.94 and 1.08.  Its 1 ms run reproduced the M4long run
     fig9-tev15 bit for bit (profiles, fields, the balance file but its
     booking columns).
   - Reported, not gated (section 28, item 7): the time-step error while
     gas is present, first order: on the paper-mode carrier 7.1e-3 at 5 s
     (the largest L1(2 ms, 1 ms) with more than 1 % of the release
     deposited), 1.2e-4 at 20 s (the tail of the cloud), the gas profile at
     16 s 3.7e-4; on the 007 carrier the temperature model 7.4e-3 at 4 s
     (the largest; 3.2 % of the release deposited), 1.3e-5 at 6 s (the end
     of the release) and 6.2e-7 at 7 s, the gas profile at 7 s 5.4e-4,
     observed orders 0.98 to 1.07.

7. **Reuse and cost of the long runs.**
   - The 11 runs of the second study (made with the frozen build of
     5852f8d, binaries be159fab...) were imported into the work area (its
     staging helper `importRecordedRuns`, not a repository file: copied
     runs and binaries, the links of the mesh and the wall table
     repointed, the meshes checked md5 for md5) and revalidated
     (`M7LONG_REVALIDATE=yes`): for every run the set-up made now equals
     its dictionaries without comments and its carrier bit for bit, and 5
     steps from 29 s with the recording's binaries and the current ones
     wrote the same 72 files, profiles and balance file, the booking
     columns included (0 differing rows in all 11): the sources of M9
     and of this integration do not change the runs of the temperature
     model.  Cost: 2 min on 2 x 8 ranks.  Saved: about 4.5 h of wall
     clock on 2 x 8 cores.
   - The six runs on the 007 carrier were made again: the first study's
     were made by its scripts with the build of the M5 review, round 4
     (cd53250f, before M4), under no key of this driver, and the HKS
     inputs changed (item 8).  Cost: 1 h 37 min of wall clock on 2 x 8
     bound ranks, the two models side by side at each dt (47, 30 and 20
     min at 1, 2 and 4 ms; 0.18, 0.23 and 0.29 s per step), 26
     core-hours.
   - A point of this machine (2 sockets x 10 cores, hyperthreading: 40
     logical CPUs): two 8-rank runs started unbound (`mpirun --bind-to
     none`, as every test) next to the four serial Graetz runs took about
     1 s per step, with ranks on both hyperthreads of a core; the same runs
     with the driver bound to 16 logical CPUs of distinct physical cores
     (`taskset -c 0-15`) took 0.15 to 0.3 s.  `longStudies` binds itself
     and its runs only on request, `M7LONG_CPUS=0-15` (it re-executes itself
     under taskset), since the CPU list depends on the machine; the readme
     and the usage of Alltest say so.

8. **Section 28, item 8: the final tolerance of HKS cases.**
   `Y_PbI2_gFinal` of run/013 is GAMG 1e-14 (was 1e-12), and its
   `out_excerpt` was regenerated.  Token by token against the excerpt of
   the M4 review, round 3: 61 of its 1656 tokens differ.  Not only the
   solver-defect digits, as expected, but all at or below the residual of
   the old final solve:
   - the final solve's residual and iterations (the tolerance itself: 32
     to 35 iterations to 7e-13, now 44 to 47 to 7e-15);
   - `solver defect` and `solverDefect` (2.6e-16 instead of 3.0e-14 of
     the inventory after 5 steps) and its relative value;
   - the closure by one ulp of the inventory (9.2e-22 mol) in 3 of the 5
     steps, and `clamped`, the rounding of the stored condensate;
   - the INLET and OUTLET transport (1.6e-18 and 1e-60 kg/s, by up to
     1e-3 relative) and the cumulative transport of the balance line;
   - the minimum Y and the negative mass of the overshoot near the boat
     (-6e-7 to -5e-9 in Y, by up to 2e-4 relative, i.e. up to 1e-12 in
     Y).
   Every other token is unchanged: the inventories (gas, wall, sample),
   the max mole fraction and the cells above the warning, the maximum Y,
   the regimes, the outer iterations and the residuals of the loose
   correctors.  Per kind: 10 of the final solve (5 residuals, 5 iteration
   counts), 10 of INLET and OUTLET, 5 of the solver defect, 5 of the
   cumulative transport, 5 of `clamped`, 10 of `solverDefect` with its
   relative value, 6 of the closure (value and relative, 3 steps), 8 of
   the minimum Y (monitor and min/max lines), 2 of the negative mass.  All
   61 come from the tolerance: the build of the cleanup run at the old
   1e-12 reproduces the old excerpt (IDENTICAL, 111 lines; section 35,
   item R, which also completes the readme of run/013).  The M5 criteria
   that run run/013 (M5.2g, M5.2m, M5.3a-d) pass with it (item 12).  On the
   007 carrier, after 16 s, the HKS
   solverDefect is now 1.4e-10, 3.9e-10 and 5.4e-10 of the inventory at 1,
   2 and 4 ms (2.2e-8 to 6.0e-8 with 1e-12,
   above M8's gate of 1e-8).  The other HKS case of ours, run/015, has had
   1e-14 since M9 (section 32, item 13); the test cases of tests/hks and
   tests/gems keep the solvers they test.

9. **Section 28, item 6: halfCell is the wall law of physical runs.**
   Recorded as follows: the code default of `wallResistance` stays `none`,
   the law of the paper and of T-Flows, so the code-to-code cases (run/014,
   M4long, M6) and every case up to run/015 keep it, and a dictionary
   without the key behaves as before; the physical cases of M8 (run/017
   and later, section 32, item 12) set `wallResistance halfCell` in their
   `thermochemistryProperties` (with `layer faceCells`, `areaMultiplier 1`
   and `areaPerVolume cell`, which halfCell needs, section 17, item 12;
   it is refused with layer distance, M4.8).  M7.2f and M7.5c gate its
   mesh independence on the paper-mode carrier.  The readme (dictionary
   and the half-cell paragraph) says so.

10. **Section 30, item 14: the carry-overs of the M4 review.**
    a. **A criterion that discriminates the face-flux correction of the
       transport defect.**  M4.18 checked only the debug report (E = 0 in
       every cell), which a defect formed differently can print as well.
       `matrixDefect.H` now has the cell defect as a template of its
       accumulator (`cellDefectT`; the ledger's `cellDefect` is its long
       double instance, its results unchanged bit for bit: every
       out_excerpt and the revalidation probes of item 7 reproduce), and
       the debug report adds the transport identity

         D = sum_c r_c + sum_c (m_c - m0_c)/dt + sum_b (F_b + C_b)
             - sum_c E_c

       over all ranks (b the faces of the non-coupled patches), with r_c
       of `cellDefectT` accumulated in `exactSum` (a long double with the
       exact rounding error of every addition, Neumaier) and every sum
       exact, the reduction over the ranks included (each long double
       split into two doubles).  In the flux form every internal face
       flux enters its two cells once with opposite signs and a processor
       face its antisymmetric mean, so D = 0 up to the rounding of the
       compensations: on the drift pipe of M4.15 (4 ranks) |D| <= 2.4e-40
       kg/s, at most 1.1e-36 of the magnitude of its terms (2e-4 kg/s), in
       every step of both layer rules; the channel with linearUpwind (E in
       2400 cells) up to 3.7e-41 kg/s, at most 3.8e-38 of its terms
       (corrected in section 35, item R: 9e-42).  (The review suggested D
       in long double.  With the ledger's long double cell defects (a
       one-off build of that variant) D on the same pipe is 4e-27 to
       3.7e-25 kg/s, the rounding of the cell sums: it overlaps the source
       form's 3.4e-25 to 1.8e-23 kg/s over both layer rules and would not
       separate the two in every step.  With the exact accumulation the
       flux form lies 15 orders below.)  Negative control, built by the
       test: a copy of the sources whose `cellDefectT` starts r_c from b_c
       - s_c and leaves out the face values C_f, while the reported E is
       unchanged (`checks.py sourceFormMutant`): it passes the check of E
       (0 in 0 cells, the old M4.18) and violates the identity in every
       step, the double-precision residue of V div(C) that the source form
       left in the closure (section 29, item 2): with layer distance, as
       the test runs it, |D| 1.2e-24 to 1.8e-23 kg/s, 5e-21 to 7.9e-20 of
       its terms; with layer faceCells (measured once by the cleanup)
       3.4e-25 to 7.8e-24 kg/s, 1.3e-21 to 3.3e-20.  (Corrected in section
       35, item R: '3e-25 to 1.8e-23 kg/s, 1e-21 to 8e-20' was the union of
       both rules, attributed to the control, which runs one.)  M4.18
       gates |D| <= 1e-30 of the terms on both drift pipes and the
       channel, and requires the
       control to violate it; the copy is built into the work area when
       the sources change (about 2 min) and runs the drift pipe with layer
       distance (about 1 min).
    b. **M4.14: the start-up count of the tetrahedral pipe** on 4 ranks is
       gated (`checks.py mdepCount` with 21 and 0, as M4.7h), where it
       only had to name some remote elements.
    c. **Stale comments:** the M4.7g block of `tests/paperMode/Allrun`
       (the sums are added in ascending rank order by `pushToFaces()`, not
       in the order of a schedule); the section header of `nearestFaceSums`
       in `interfaceExchange.C` ('push; plusEqOp'); the residue figure of
       `matrixDefect.H` (typically 1e-25 to 8e-25, at most 1.7e-24 kg per
       step, not 'about 1e-25'); the header of `setFrozenCarrier.C` (the
       lines 670-683 belong to `Read_Controls_Mod/Boundary_Conditions.f90`,
       and 'rows all' reaches 0.69 m).
    d. **The readme's two-way statement**: with the build before and after
       the fix of section 29, item 1, the fields, the per-rank state and
       the profiles of those runs are identical but for mDep_ itself; the
       balance files differ only in the booking columns (`solverDefect`,
       `closure`) and the monitor `evaporatedAboveWarning` (it said 'every
       other file of those runs is identical').

11. **Other changes of this round.**
    - `tests/paperMode/longRuns` (M4long) keeps a recorded run that ended
      under another key and is not selected (NOTRUN, with the way to
      revalidate it or to run it again), as M7long does (the second
      study, g).  Before, an evaluation (`M4LONG_RUNS=ended`) after any
      change of the sources without `M4LONG_REVALIDATE=yes` set every
      recorded run up anew, i.e. deleted it: since M9 and this
      integration change the sources, that would have discarded the 14
      runs of M4long.  Table 2 (M4.5) no longer evaluates such a run as a
      partial one.
    - `tests/Alltest` expects M9.10 (recorded by `tests/gems/Allrun` since
      M9, but not in the list of expected ids, so its absence would have
      passed unnoticed).
    - `longStudies` decomposes a prepared run again when the number of
      ranks changed before its first write (found when the runs were
      restarted pinned, item 7).

12. **Quality.**
    - Build: `wclean` and `Allwmake` with 0 warnings (BUILD); a copy
      compiles and links against the local v2606 build with 0 warnings
      (M4.6: OpenFOAM v2606, linux64GccDPInt32Debug, every source of the
      solver and the utility compiled).  M0.4 (v2606 on Merlin7) is NOTRUN.
    - `tests/Alltest -fresh all M7long` (one invocation, against the
      pristine tag 0.1, with `LESTO_OPENFOAM_V2606` set, bound to the 20
      physical cores with `taskset -c 0-19`): 299 passed, 0 failed, 1 not
      run (M0.4) in 29 minutes at a load of 1 to 17.  `all`: 283 (the 258
      of M9 and the 25 of M7.1); M7long: 16 (M7.2a-b and the 14 of the
      long runs, every run reused: the six on the 007 carrier made with
      the current sources, the 11 of the second study revalidated).  The
      Graetz runs (made the same day) were reused as well: their key holds
      the md5 of the solver binary, which the build from scratch reproduced
      (38e51df0...; setFrozenCarrier f7c0df84...).
    - M4long after M9 and this integration: its 14 recorded runs, last
      revalidated after the rebase (section 31), revalidated under the
      current sources (`M4LONG_REVALIDATE=yes M4LONG_RUNS=ended`, 2 min on
      16 cores): every probe wrote the same files, profiles and balance
      file, the booking columns included (0 differing rows); 17 passed,
      M4.17 not run (no M4long run was made with the current sources; the
      closure of runs made with them is gated by M7.3a).
    - The tests leave the repository's files unchanged (the work area
      `tests/work/`, git-ignored; checked by the md5 of every tracked and
      untracked file before and after).
    - Negative evidence: M4.18 (item 10a); the coarse triplet of none at
      Bi 50 (item 3); M7.2e; the evaluation of the ended M7long runs with
      impossible limits (`M7LONG_RUNS=ended`, closure 1e-17, solverDefect
      1e-13, dt L1 1e-8, h vs h/2 0.001 bins and 0.001 %, axial L1 1e-5):
      every gate fails (M7.3a-c, M7.2d, M7.2f, M7.5b-c, M7.3pb), the
      reported and the listed criteria pass (M7.3d, M7.2e, M7.3pc; M7.2c,
      M7.5a, M7.3pa, whose runs are revalidated).  The same-deposit check
      of M7.2c (flux form against A on the base mesh, fixed limits 0.01
      bins and 0.001 %) fails for none against halfCell on the base mesh
      (`longChecks.py pair`: the peak -0.06 %).

13. **Decisions needed.**
    a. The split of M7 (item 1): Graetz-Robin in `all`, every pipe study
       in M7long.
    b. graetzCarrier as a test helper (item 2).
    c. M7.1k against mu_h to 1e-9 (item 3e).
    d. M7.2f a gate (item 4) and the axial criterion a gate (item 5).
    e. The reading of 'T_dep shift <= 1 bin' (item 4).
    f. The record of halfCell for physical runs (item 9).
    g. The transport identity as the criterion of the carry-over (item
       10a): exact sums and a limit of 1e-30 of the terms.
    h. M4long keeps stale recorded runs (item 11).
    i. Section 32 (M9) ends with 'M9 status: see the end of this section'
       without the status: M9's own sign-off.

---

## 34. Researcher decisions of 2026-10-04 (signed off)

Sections 32 (M9) and 33 (M7) are signed off as recommended, in particular:

1. **M9 (section 32).**
   - The `excess` entry of a pair's `gems` block (an iodine excess of
     1e-6) against the GEMS3K stall at exact PbI2 stoichiometry.
   - Lagged speciation as chi = n_g/c, instead of E9's p_i^GEMS/p_i.
   - Mode frozen: the GEMS table at the nodes of the pair's table, computed
     on the master and broadcast.
   - M9.3b (at most 70 us per warm call) is NOTRUN rather than FAIL when
     the machine is busy (1-minute load of 4 or more); on a quiet machine it
     measured 42 us.
   - Case numbering: run/015-hks-gems-local is the GEMS case; the HKS
     code-to-code case of M6 becomes 016 and the physical carriers of M8
     017 and later.
2. **M7 (section 33, item 13 a-h).**  The split of M7 (Graetz-Robin in
   `all`, the pipe studies in M7long); graetzCarrier as a test helper;
   M7.1k against mu_h to 1e-9; M7.2f (halfCell, h vs h/2) and the axial
   x2 criterion as gates; 'T_dep shift <= 1 bin' as the shift of the
   interpolated onset in native bins; halfCell for physical runs through
   the M8 cases, the code default staying `none`; the transport identity
   with exact sums and a limit of 1e-30 of its terms as the criterion of
   the face-flux correction; M4long keeping stale recorded runs.  Item 13i
   (the status of section 32) is completed by the cleanup that follows.
3. **A cleanup pass** for the minor findings of the M9 and M7 reviews
   follows; its amendments are recorded in section 35.

---

## 35. Amendments from the M9 and M7 reviews: the cleanup (signed off by the researcher, 2026-10-04)

The reviews and verifications of M9 (section 32) and M7 (section 33) left
minor findings, items A to S below: A to M of M9 (`gemsEquilibrium`, its
tests and documents), N to R of M7 (`tests/verification` and the
documents), S the gates of the long runs.  Each was assessed first.  All
were real, and M is a record with a decision.  Four came out otherwise in
a detail when measured: A had a second defect behind it (OpenFOAM v2412
writes an int64 as a label); F was more severe (a NaN log10 Omega stopped
the run with a floating-point exception instead of being booked
silently); in L the bridge is linked against glibc 2.34, not 2.38; in N
'cannot be computed' already failed, and `__main__` did not catch
exceptions, which is how a crash read as PASS.  Every fixed defect has a
new or an extended
criterion: M9.3b, M9.5, M9.6, M9.7 and M9.8 (more lines) and M9.11 to
M9.13 (new) in M9; M7.1n and M7.1o in M7; M7.6a and M7.6b in M7long; and
every run of the bridge build gates its closure and its solver defect
(item B).  The sources before this pass are kept as the snapshot
`pre-cleanup9-7` (scratch area of this pass); its bridge build fails
every new or changed criterion of M9 whose defect is in the solver, and
its scripts fail M7.1n, M7.6a and M7.6b (item T).

A. **The counts of mode local are 64-bit integers (`gemsEquilibrium.H`,
   `endUpdate()`, the reports).**
   - Defect: the counts of an update and of the run (calls, warm, cold,
     retried, failures, re-creations, IPM iterations, the sources of each
     pair, the clipped chi) were labels, 32-bit in the Int32 builds.  A
     production run overflows them (signed overflow, undefined behaviour)
     and corrupts the summary: 60 s at 1 ms of the 15264 wall elements of
     the pipe are 9e8 calls and, at about eight per call, 7e9 iterations.
   - Fix: `std::int64_t`; per update one reduction of every integer count
     in 64 bits (one list, one `MPI_Allreduce`), one of the seconds and
     one of the extremes of the pairs (chiMin negated).  The test found a
     second defect on the way out: OpenFOAM v2412 writes an `int64_t` to
     an Ostream as `label(val)` (`int64IO.C`), i.e. truncated to 32 bits
     in an Int32 build (v2606 writes it in full), so a rank log printed
     690588672 for 7.8e10 iterations.  Every count is written through
     `std::to_string` (`num()`), exact in both versions.
   - Tests: M9.13 (new).  The fault bridge reports 1e9 IPM iterations for
     every call (`LESTO_TEST_ITERATIONS`), 3 steps on 4 ranks: every update
     line with calls and the summary print 1e+09 iterations per call, and
     every update line of the 4 rank logs the total calls x 1e9 exactly
     (up to 8e10, 37 times 2^31).  The snapshot: 8.85e+06 and 4.82e+06
     per call, the summary 6.79e+06, and the rank logs 690588672 and
     -1604378624.
   - Documents: `gemsEquilibrium.H` (Logs, the counts), the readme.

B. **The ledger of every GEMS run (`tests/gems/Allrun`).**  The closure and
   the solver defect were gated in M9.2a, M9.3a(ii) and M9.8 only.  Every
   run of the bridge build that ends now gates both at the M5 level, the
   closure at most 1e-13 of the reference in every step and |solverDefect|
   at most 1e-11 at the end (the helper `ledger`): the 30-step sources run
   (M9.3a), 4 ranks (M9.3c), the table run (M9.3d), paper mode allowed and
   frozen (M9.3e), speciation lagged (M9.4), both fault runs (M9.5), the
   updateInterval run and its two restarts (M9.6), the restart from 0.2
   (M9.9), the three halfCell runs (M9.10), the frozen restart and 4 ranks
   (M9.2b), and the 1-step runs of M9.11 to M9.13.  Measured over the 28
   runs of the bridge build in the run of this pass: the closure at most
   1.9e-16 of the reference (the channel), |solverDefect| at most 2.60e-13
   (the 4 ranks; the serial runs at most 2.49e-13, the frozen runs and the
   short ones at most 9.8e-14).  The claim of section 32, item 14
   ('closure <= 1.9e-16 and solverDefect <= 2.5e-13 in every run') held
   for the runs it gated then, not for the 4 ranks, which it did not gate:
   corrected there to 2.6e-13.  (The 1-step run of the bridge build in
   M9.1 was the one run left without the ledger; it is gated since item
   X.  The stand-in run of item V, run015-otherBridge, got its gate in the
   last check of this pass: closure 1.9e-16, solverDefect 3.3e-16.)

C. **The updates of mode local follow the time index, and their lagged
   state is restart state (`gemsEquilibrium::updateDue()`,
   `interfaceExchange::restoreEquilibrium()` and `storeElementState()`).**
   - Defect: `start()` set `stepsSinceUpdate_ = updateInterval - 1` in every
     run, so every run updated at its own first step and counted from
     there, and p_eq, chi and the sources of the last update were not
     written.  With updateInterval > 1 a restarted run applied other,
     lagged values than the continuous one: the snapshot's restart from
     0.004 of the updateInterval-3 channel updated at the steps 5, 8 and 11,
     the continuous run at 1, 4, 7 and 10.
   - Fix: an update is due at the steps whose time index n
     (`Time::timeIndex()`, restored from `uniform/time` on a restart) has
     (n - 1) mod updateInterval = 0, and at the first step of a run whose
     elements hold no lagged state.  That state is written with the per-rank
     state, `uniform/phaseChange/phaseChangeEquilibrium` (mode local only):
     a binary list behind the fingerprint of all elements (as
     `phaseChangeRegimes`), then per pair the counts of the last update for
     its step line (`gemsEquilibrium::lastUpdate()`, 10 values) and per
     pair and element p_eq, chi and the source.  A restart onto the same
     elements restores it, the law of every element recomputed with
     `hksCoefficients()` (the doubles of the update), and says so ('the
     lagged state of the last update ... restored ...; the next update at
     the step of time index 7').  Without it (missing, a collated container
     of another number of ranks, other elements) the first step updates,
     with a warning when the schedule has no update there.  With
     speciation none chi is 1 whatever the file holds.  The warm states
     are still not written (section 8, risk 13): the next update is cold.
   - Tests: M9.6, now 3 lines, on the updateInterval-3 channel written every
     step.  (ii) The restart from 0.004 updates at 7 and 10, not at 5 and
     6; every file of 0.005 and 0.006 (32) is cmp-identical to the
     continuous run's and their step lines equal its lines; the update at 7
     is cold; at 0.012 the inventories lie within 2.6e-8 of the continuous
     run.  (iii) Without the file: the warning, updates at 5, 7 and 10.  The
     snapshot: updates at 5, 8 and 11 in both.
   - Documents: `gemsEquilibrium.H` (mode local), `interfaceExchange.H` (the
     state; the GEMS block), the readme, section 32 item 5 (annotated).

D. **A log directory per rank (`gemsEquilibrium::logDirectory()`).**
   - Defect: a `logDirectory` given absolute, as `<case>/...` or
     `$FOAM_CASE/...`, or relative but leaving processorN with `..`, was
     used as it was on every rank: all ranks appended to one
     `gemsEquilibrium.log`, and the spdlog sink of GEMS3K's `ipmlog.txt` was
     opened by several processes.
   - Fix: the directory is expanded, made absolute against the rank's case
     directory (`Time::path()`) and cleaned; in a parallel run one that does
     not lie below processorN gets a subdirectory `processor<N>`.  The
     start-up line names it ('logs in .../gemsLogCase/processor<N>'; the
     default prints 'processor<N>/gemsLog' as before, so run/015's excerpt
     is unchanged).
   - Tests: M9.12 (new): 4 ranks with `"<case>/gemsLogCase"`,
     `"$FOAM_CASE/gemsLogEnv"` and `"../gemsLogUp"`: every rank has its own
     `gemsEquilibrium.log`, whose only rank line is its first, with
     `ipmlog.txt` next to it, in a directory of its own.  The snapshot: one
     log with the four rank lines in each case.
   - Documents: `gemsEquilibrium.H` (Logs), the readme, section 32 item 7
     (annotated).

E. **Mode frozen re-creates the engine after GEMSB_ERR_FATAL
   (`buildFrozenTables()`).**
   - Defect: the failed node was logged, and the remaining nodes went on
     with the same engine, against the contract of the bridge
     (`GEMSB_ERR_FATAL`, GEMS3K's T_ERROR_GEM: the engine must be
     re-created).
   - Fix: as mode local: the engine destroyed and created again from the
     system file, every condensed species suppressed again
     (`createEngine()`); the calls of mode frozen are cold, so no warm
     state survives.  The start-up line counts the re-creations ('..., 0
     rejected (taken from the table), engine re-created 0 times; ...').
     An unusable result of a converged call (item F) fails the node too.
   - Tests: M9.5, a second line: the hot channel with the fault bridge
     returning GEMSB_ERR_FATAL at the node call 20: the engine destroyed
     and created again, the 4 species suppressed on the new one, every later
     call on it, all 152 calls cold; the start-up line '151 from GEMS3K,
     1 failed, 0 rejected, engine re-created 1 time'.  The snapshot: no
     destroy and create after the fault.  The fault bridge now names the
     engines by their creation number: in a trial run the re-created
     engine had the address of the one destroyed before it.
   - Documents: `gemsEquilibrium.H` (mode frozen), the readme, section 32
     item 4 (annotated).

F. **An unusable result of a converged call is a failure (`evaluate()`).**
   - Defect: a call that converged and balanced but gave the pair p_g <= 0
     or a max log10 Omega that is not finite was booked as source FAILED
     but not counted among the failures and not logged.  Measured with the
     snapshot and the fault bridge with `FOAM_SIGFPE=false`: 4 elements
     'failed' in the step lines, 2 failures in the summary (the 2 fatal
     ones), neither unusable result in the rank log.  With the default
     floating-point trapping the snapshot stopped at the NaN with a
     floating-point exception: Foam::max is (a > b) ? a : b, an ordered
     comparison, which may raise FE_INVALID for a NaN, and OpenFOAM traps
     it (`sigFpe`).
   - Fix: such a call counts once as a failed call and is written to the
     rank log ('converged and balanced (OK, 9 iterations), but PbI2_g p_g
     0 Pa, max log10 Omega -0.035; the table'), its warm state
     invalidated.  log10 Omega, p_g and the sums of the element balance are
     compared only once they are known to be finite (`maxLog10Omega()`,
     `positive()`, `balanced()`), in both modes.
   - Tests: M9.5, first line, extended: besides the fatal faults at the
     calls 150 and 600, the fault bridge reports p_g = 0 after the
     converged call 50 and log10 Omega = NaN after the call 100: 4 failures
     in the step lines and in the summary, 2 re-creations, the rank log
     names 2 ERR_FATAL, 1 'p_g 0 Pa' and 1 'max log10 Omega nan', and the
     next call of such an element is cold.  The snapshot: the run stops
     (SIGFPE).
   - Documents: `gemsEquilibrium.H` (FAILED), the readme.

G. **Only the keys of the mode are read (the constructor).**
   - Defect: every key was read and validated in both modes, so
     referencePressure <= 0 stopped mode local, and minMoleFraction < 0 or
     updateInterval 0 mode frozen, although the warning says that the mode
     does not read them.
   - Fix: the keys of both modes (system, mode, carrier, minTemperature,
     maxLog10Deviation, logDirectory, logLevel, writeEquilibriumField) are
     read and validated in both, the others in their mode only; those of
     the other mode are named in the warning and never read.
   - Tests: M9.11 (new), first line: mode local with referencePressure -1
     and referenceComposition { He -1; }, mode frozen with minMoleFraction
     -1, updateInterval 0 and allowNonPhysicalCarrier maybe run 1 step
     each, with the warning naming them.  The snapshot stops both.

H. **Unknown keys stop the run (`checkKeys()`).**
   - Defect: the header said that the constructor 'refuses unknown or
     inconsistent keys', but a misspelled key (minMolFraction) was ignored
     and its default used.
   - Fix: every key of GEMSCoeffs (13 known) and of a pair's gems block (4)
     must be known; any other stops the run, naming it and the valid keys.
   - Tests: two lines of M9.7: 'GEMSCoeffs: unknown key minMolFraction.
     The valid keys are system mode carrier ...' and 'Pair PbI2_g, gems:
     unknown key excesss.  The valid keys are gas condensates elements
     excess.'  The snapshot runs both.

I. **The untested start-up checks of section 32, item 2.**  Five more
   lines of M9.7: the molar mass of the system against molarMass (0.4610
   for 0.46100894, 1.9e-5 off), a condensed species as the pair's gas
   (PbI2(cr)), a condensate of another formula ('Pb(cr) (Pb 1) has another
   formula than PbI2(g)'), referencePressure 2e5 Pa in mode frozen ('lies
   outside the pressure grid'), an unknown mode (lokal); and the warning of
   equilibrium table about the inputs of equilibrium GEMS ('equilibrium
   table does not read GEMSCoeffs pairs/PbI2_g/gems', M9.11, second line).
   These checks were right (the snapshot passes the six); their negative
   evidence is the input that triggers them.

J. **The excerpt of run/015 without the numbers of GEMS3K (`regress.sh
   -mask`).**
   - Defect: M9.8 compared the iterations per call and the six digits of
     the largest |dlog10 p_eq/p_eq,table| of the GEMS lines, numbers that
     GEMS3K computes itself: a bridge built with another compiler or
     GEMS3K version may print others while the solver is right.
   - Fix: `regress.sh -mask <regex>` replaces every match in both files by
     `<masked>`; M9.8 masks those two (`GEMS_MASK`, 10 tokens of the 116
     lines); they stay printed.  A GEMS3K that gave another p_eq would
     still change the solver's own numbers, which are compared.
   - Tests: M9.8, a third line: the excerpt with those numbers replaced is
     IDENTICAL with the mask and DIFFERENT without; with the gas inventory
     of its first balance line changed as well, DIFFERENT with it.
   - Documents: `regress.sh`, the readme of run/015.
   - Superseded by item V: the mask does not make the comparison pass with
     another bridge (GEMS3K's numbers change only together with its p_eq,
     and then every round-off token of the later steps moves as well); the
     mask and `-mask` were removed again.

K. **The load sampler of M9.3b never outlives its script
   (`startLoadSampler`, `stopLoadSampler` in `tests/functions`).**
   - Defect: an endless background loop, killed only after the timing run
     returned: an interrupted run left it running.
   - Fix: the sampler checks every 2 s that its script lives and ends when
     it does not (a KILL runs no trap); the EXIT, INT and TERM traps of the
     script stop it.
   - Tests: M9.3b, a second line (also without a bridge): a script that
     samples, terminated (TERM) or killed (KILL), leaves no sampler 5 s
     later, while the endless loop it replaced, as the control, outlives a
     killed script (the test then stops it).

L. **Documents of M9.**  The status of section 32 (written below its item
   16, as section 34 asked); the test times of the readme (M9, M7 and all,
   measured in this pass); section 6 (criteria M9.0 to M9.13, M9.14
   since item X);
   `thermochemistry/README.md` (the bridge is wired into the solver; its
   next steps); `libgemsbridge.so` in the readme's Building and in
   `thermochemistry/README.md`: it needs glibc >= 2.34 (it is linked against
   the glibc 2.34 of the conda-forge sysroot, symbol versions up to
   GLIBC_2.34; the 2.38 of the review is that of this workstation), it
   carries an RPATH into the local conda environment, and elsewhere, e.g.
   on Merlin7, it is rebuilt with `build-gems3k-static.sh`; run/015's
   readme (step 1 makes no call, no element having gas at t = 0; 576 to 720
   calls in the steps 2 to 5).

M. **A known limit, recorded: GEMS3K in a gas that is not dilute.**
   Reproduced on the channel from 1300 to 600 K (the boat at the hot end,
   mode local, 10 steps): 37 calls fail with ERR_NOCONV, at x = 0.05 to
   0.73 and 1170-1240 K, after 6 to 7452 IPM iterations (mean 2043); the
   converged warm calls take about 4 ms there (1000-2000 iterations), and a
   step 0.52 s against 0.006 s with the table (a failed call about 50 ms,
   estimated from those times).  Every failed call takes the table and is
   counted and logged, so the results stay correct.  Decision (my call,
   for the sign-off): no guard now.  Such a gas lies outside the dilute,
   passive-carrier assumption of the solver, which the mole-fraction
   monitor already flags (cells above moleFractionWarning, 0.05); a
   `maxMoleFraction` guard would add a source type and a column to every
   GEMS line, and its threshold depends on the system, so it belongs to the
   qualification of the guards with Pb-Bi-I in M10 (section 8, risk 7).
   Recorded in the readme (GEMS section) and `thermochemistry/README.md`.

N. **A crash of `longChecks.py` is not a PASS (`longStudies`,
   `longChecks.py`).**
   - Defect: the reported criteria (`pairCriterion` expectFail and report,
     `axialCriterion` with `M7LONG_AXIAL_GATE=no`) accepted any exit status
     up to 1, and an uncaught exception of python exits with 1, without a
     verdict line.  ('Cannot be computed', exit 2, already failed.  The
     review said that `__main__` catches exceptions; it did not, which is
     how a crash came out as 1.)
   - Fix: `longChecks.py` prints 'CRASHED: <exception>' last and exits with
     3; a reported criterion needs `computed()`: the status 0 or 1 and its
     verdict line ('OK: ...' or 'FAILED: ...') last.
   - Tests: M7.6a (new, M7long, no run): the helpers with a stand-in
     checker (met, not met, cannot, a crash with python's exit 1, a crash
     reported with exit 3) in the 5 modes give the 25 expected verdicts,
     and `longChecks.py` with a number it cannot read 'CRASHED:
     ValueError', exit 3.  The snapshot's helpers give 3 of the 25 wrong
     (a crash PASS in expectFail, report and the reported axial criterion),
     its `longChecks.py` exit 1 without a verdict line.

O. **A fit that is not finite fails M7.1e-g (`graetzChecks.py`).**
   `discrete`, `windows` and `outlet` took the largest deviation with max(),
   and max(0, nan) is 0: a NaN fit passed.  A deviation that is not finite
   now fails the check, naming the run.  M7.1n (new, no run): the three
   checks on stand-in fits, the finite ones met, the five with a NaN (mu,
   a window variant, the wall uptake, the mixing cup) refused; the
   snapshot's `graetzChecks.py` passes all five.

P. **A resumed run of the 007 carrier starts from its latest write
   (`runM7`).**  `setup007` copies run/012 or run/013, whose controlDict
   says `startFrom startTime`, and `runParallel` passes no time option, so
   a stopped run started again from t = 0.  `runM7` now sets `startFrom
   latestTime` when it resumes (the revalidation leaves startFrom out of
   the compared dictionaries, as endTime).  No recorded run was affected:
   all 17 have one segment.  M7.6b (new, M7long): run/012 set up by
   `setupM7` on 8 ranks, stopped at 0.008 s and resumed by `runM7` to 0.016
   s: the second segment starts with 0.012 s, makes 2 steps, books its
   balance from 0.008 s and restarts exactly (relative 0).  With the
   snapshot's `runM7` (a copy of `longStudies` with it) the second segment
   starts again from t = 0: its first step 0.004 s, 4 steps, no balance
   from 0.008 s (M7.6b FAIL).

Q. **The reuse keys of the Graetz runs follow their set-up
   (`graetzRobin`).**  The keys held the executables, the template,
   `graetz.py`, `graetzCase.py` and the options, but not the functions that
   set the runs up: a changed dt, step or write-interval formula of
   `transientCase`, or a changed command of `meshCase`, reused the old
   runs.  The keys now hold the md5 of `declare -f meshCase runCase
   transientCase meshKey runKey transientKey setupKey meshDir runDir
   transientDir runSerial` (as M4.20 and `m7Key`).  M7.1o (new, no run):
   with `meshCase` (the order of its postProcess functions), `runCase` (5
   steps) or `transientCase` (a write every 0.1 s) changed in a subshell,
   the keys of a mesh, a run and the time-accurate run change; the
   snapshot's keys hold no function text and stay.

R. **Documents of M7 (corrected in place, each annotated).**
   - Section 33, item 3d: the largest mixing-cup (7.5e-11) and wall-uptake
     (8.1e-10) changes are at Bi 50 (halfCell, Nr 256; the outlet
     corrected by item Z: fixedValue, not zeroGradient); only the one-end
     and both-ends maxima (7.5e-9, 1.2e-8) are at Bi 1, also with the
     fixedValue outlet.
   - The transport identity (`matrixDefect.H`, section 33 item 10a,
     `tests/paperMode`, the readme): the linearUpwind channel measures up to
     3.7e-41 kg/s (3.8e-38 of its terms), not 9e-42; the negative control
     of M4.18 runs the drift pipe with layer distance (20 steps on 4
     ranks): 1.2e-24 to 1.8e-23 kg/s, 5e-21 to 7.9e-20 of its terms.  The
     figures given for it, 3e-25 to 1.8e-23 kg/s and 1e-21 to 8e-20, were
     the union over both layer rules: the control with layer faceCells,
     measured once in this pass, gives 3.4e-25 to 7.8e-24 kg/s, 1.3e-21 to
     3.3e-20.  The long double variant of the review (4e-27 to 3.7e-25
     kg/s, a one-off of the integration) overlaps that union, as section
     33 said.
   - The readme's row of M7.1c: the gate takes `graetz.half_cell`, the
     centroid formula in double precision on the unit radius times R
     (4.0e-14), not the exact centroid (3.4e-14); OpenFOAM's centroid gives
     3.7e-14.
   - run/013 (its readme, section 33 item 8): the 61 changed tokens per
     kind, the cumulative transport of the balance line included, and all
     from the tolerance: the build of this pass at the old final tolerance
     1e-12 reproduces the old excerpt (IDENTICAL, 111 lines, 5 steps).
   - The temperature model's final deposit (readme, section 33 item 6):
     it still depends on dt at first order, but weakly: L1(2 ms, 1 ms)
     2.7e-7 at 16 s on the 007 carrier (order 0.90) and 1.3e-7 at 30 s in
     paper mode (0.94), 1600 times below the 4.5e-4 of HKS; not 'no longer
     depends on dt'.
   - The closed box of the deliverables of M7: verified by the analytic
     closed boxes of M2.1 and M5.1 (section 33, item 1, and section 6).
   - Stale statements: `longStudies` (M4long keeps such runs too since
     section 33, item 11; M7long always runs pipeMeshes first), section 6
     (M9.0 to M9.13, M9.14 since item X; M7.1a-o, M7.6a-b).

S. **The gates of the long runs, live again (cost).**  This pass changes
   compiled sources (items A to H), so every recorded run of M4long (14)
   and of M7long (17) ended under another key, and M4.17 and M7.3a would
   gate none.
   - The cheapest run of each family was made again with the final
     sources, each on 8 ranks bound to 8 physical cores:
     temperature-dt4ms (M7long; 4000 steps of 4 ms, 17.4 min on cores 8-15,
     alongside the run of `tests/Alltest -fresh all`) and T5_60-distance
     (M4long; 6000 steps of 10 ms in two segments: the first, on cores 0-7,
     ran about eight times slower than its record, since another user's
     job pinned to one of those cores shared it with a rank, and was
     stopped at 4.65 s after 11.7 min; `longRuns` resumed the run from its
     write at 4 s on cores 12-19, 18.9 min); together 48 min of 8 ranks,
     6.4 core-hours (about 5 without the contention).  Both reproduce
     their records bit for bit: temperature-dt4ms every profile file (64),
     every field of the last write on every rank (96 files at 16 s) and
     the balance file, its booking columns included; T5_60-distance every
     profile file (120), the 6000 rows of its balance file merged from the
     two segments in every column, and every field at 60 s but OpenFOAM's
     clock value in uniform/time (8 files: 59.9999999999966747 against
     ...6632, the time accumulated from the restart at 4 s; the ledger
     keeps its own exact clock).  The changes of this pass do not act on
     the temperature model, as expected.
   - The other 29 runs were revalidated (`M4LONG_REVALIDATE=yes
     M4LONG_RUNS=ended M7LONG_REVALIDATE=yes M7LONG_RUNS=ended`, in the
     final invocation; sections 26 and 27): the 13 of M4long and the 16 of
     M7long, HKS included, each by 5 steps from its last write before its
     end with the binaries it was made with and the current ones: every
     probe wrote the same files (9 to 88), the same profiles and a balance
     file identical in every column, the booking columns included (0
     differing rows); none was refused.
   - Results (`tests/Alltest M9 M4long M7long`, below): M4.17 gates
     T5_60-distance, 7.9e-16 of the reference over its 6065 steps (the 65
     of the stopped segment beyond 4 s counted once more), and lists the
     13 revalidated runs; M4.5 for T5_60 with layer distance gives T_dep
     789.2542 K as before; M7.3a gates temperature-dt4ms, closure 2.8e-16
     and solverDefect 5.0e-13 of the reference, and lists the 5 others;
     every criterion of M4long (18 lines) and M7long (18) passed.
   - Made again after the follow-up, item Z: its item Y changed the
     compiled sources once more.

T. **Tests and negative evidence.**
   - New: M9.11 (2 lines), M9.12 and M9.13 (M9); M7.1n and M7.1o (M7);
     M7.6a and M7.6b (M7long).  Extended: M9.3b, M9.5 and M9.8 by a line
     each, M9.6 by two, M9.7 by seven (16 lines); the closure and the
     solver defect gated in every run of the bridge build (item B).
     `tests/Alltest` expects them and says so in its usage text;
     `tests/gems/checks.py` has `frozenFault`, `restartSchedule`,
     `rankLogs`, `iterations` and an extended `fault`;
     `tests/gems/faultBridge.c` injects the unusable results and the
     iterations and names the engines by their creation number;
     `regress.sh` has `-mask` (removed again, item V); `tests/functions`
     the load sampler.
   - Negative evidence, the snapshot `pre-cleanup9-7`: its bridge build
     under the new `tests/gems/Allrun` (`LESTO_GEMS_TEST_EXE`, before the
     line of item I's condensed gas was added) fails 9 of its 42 lines,
     exactly those of its defects: M9.5 (both lines: the local fault run
     stops with SIGFPE at the NaN; no re-creation in mode frozen), M9.6
     (both restarts: updates at 5, 8 and 11), M9.7 (the two unknown keys:
     the runs do not stop), M9.11 (the keys of the other mode stop both
     runs), M9.12 (one shared log per directory) and M9.13 (the 32-bit
     totals);
     the other 32 pass and 1 is NOTRUN (M9.0: the given executable), as
     the changes of items B, I, J and K act in the scripts.  Its scripts:
     `graetzChecks.py` passes the five NaN fits of M7.1n; its criterion
     helpers give 3 of the 25 verdicts of M7.6a wrong, and its
     `longChecks.py` exits 1 without a verdict at a crash; its `runM7`
     restarts the resumed run of M7.6b from 0 (item P); its keys of
     `graetzRobin` hold no set-up function (M7.1o, item Q).  With the
     fault bridge and `FOAM_SIGFPE=false` its fault run counts 2 failures
     for 4 unused results (item F).
   - Build: 0 warnings (BUILD, M0.3, the bridge build of M9.0); a copy with
     the bridge compiles and links against the local v2606 Debug build
     with 0 warnings, every source compiled (one-off; M4.6 compiles the
     default build and passed).  A header comment was reflowed after the
     build of the run of `all`: the executables rebuilt from the final
     sources are bit-identical to those it tested (rhoFixedFlowFoam
     bd7c343f..., setFrozenCarrier f7c0df84..., unchanged).
   - The test runs leave the repository unchanged: the md5 of every
     tracked and untracked file (417) before and after the run of `all`
     differ only in the six documents edited during it (this plan, the
     readme, `gemsEquilibrium.H`, `thermochemistry/README.md` and the
     readmes of run/013 and run/015), with the same set of files; and
     those of the final invocation `tests/Alltest M9 M4long M7long` are
     identical before and after it (417 files, no edit during it).

**The follow-up: the verification and the review of this cleanup.**  The
verification of the cleanup above found two problems (the readme still
said that the final deposit of the temperature model does not depend on
dt; `tests/Alltest` did not expect the new criteria of M7) and four minor
points (the outlet of the largest mixing-cup and wall-uptake changes in
section 33; M9.1's run without the ledger, against item B; a NaN closure
passing the ledger's checker; the start-up line of a log directory that
comes back into processor0), its review three minor defects (the
completeness lists again; the mask of item J does not do what it claims;
the control of M9.3b(ii) is the endless loop that item K removed) and six
nits.  Each was assessed: all are real, and all are fixed, the mask by
the review's second option (the bridge's md5) with an exact comparison
of what no GEMS3K call reaches.  Items U to Z.  The sources of the
cleanup before this follow-up are kept as the snapshot `pre-cleanup9-7b`
(scratch area), the negative evidence of the follow-up.

U. **The completeness lists of `tests/Alltest`.**
   - Defect: the lists still named M7.1a-m and M7.2a-M7.3pc only, so a
     missing or duplicated line of M7.1n, M7.1o, M7.6a or M7.6b would not
     have failed the run, although item T said that Alltest expects them.
   - Fix: M7.1n:1 and M7.1o:1 in the list of M7, M7.6a:1 and M7.6b:1 in
     that of M7long, M9.14:1 in that of M9 (item X); the usage text names
     them.  And the converse check, so that a new criterion cannot be left
     out again: the criterion LIST (every invocation) requires every
     criterion id recorded with PASS or NOTRUN to be in the invocation's
     expected list (`unexpectedResults`), with a stand-in id added to the
     results as its control, which it must name.
   - Negative evidence: `missingResults` with the lists of the first round
     over the results of its own verification without the lines of M7.1n
     and M7.1o (of M7.6a and M7.6b) reports nothing missing; with the new
     lists it reports M7.1n:1 and M7.1o:1 (M7.6a:1 and M7.6b:1), and over
     the complete results nothing.  `unexpectedResults` with the lists of
     the first round over the complete results of its verification names
     M7.1n and M7.1o (`all`) and M7.6a and M7.6b (M7long); with the new
     lists nothing.

V. **run/015's excerpt with another bridge (M9.8; replaces the mask of
   item J).**
   - Defect (review): masking the numbers that GEMS3K computes itself does
     not make M9.8(i) pass with a bridge built elsewhere, as item J
     claimed.  Those numbers change only when GEMS3K's results (or its
     inputs) change, and then every round-off token of the solver in the
     later steps moves with them: the final residuals of GAMG, the solver
     defect, clamped, solverDefect, min Y and the negative mass.
     Measured with the stand-in below (every p_g scaled by 1 + 4.4e-16, 2
     ulp): 22 of the 116 lines of the excerpt differ, all in the steps
     2-5, and 19 with the mask of item J (the first round's M9.8(i): a
     FAIL); the start-up and step 1 are IDENTICAL.
   - Fix: the md5 of the `libgemsbridge.so` with which the excerpt was
     written is recorded with the test (`EXCERPT_BRIDGE_MD5` in
     `tests/gems/Allrun`: fb67103e..., `thermochemistry/gemsbridge` of
     this workstation, built on 2026-09-23).  With that bridge the whole
     excerpt, the start-up and 5 steps, is compared exactly (PASS or FAIL,
     as before item J).  With another bridge only what no GEMS3K call
     reaches is compared exactly, the start-up and step 1 (no element has
     gas at t = 0, so the update of step 1 makes no call): DIFFERENT is a
     FAIL, the solver's own; IDENTICAL leaves the steps 2-5 NOTRUN, naming
     both bridges.  On Merlin7, where the bridge is built anew (item L),
     M9.8(i) is thus NOTRUN unless the start-up or step 1 differ (an
     excerpt written there again needs its bridge's md5).  The
     mask is gone: `regress.sh` is the committed version again, and
     `GEMS_MASK` and the third line of item J were removed.  The md5
     names the bridge, not its arithmetic: libm is a dependency too, so a
     glibc of other numerics under the same bridge would still fail the
     exact comparison (none was seen: the excerpt of 2026-10-03 has been
     reproduced bit for bit since).
   - Tests: M9.8, its first line as above (here: PASS, the whole excerpt
     IDENTICAL, 116 lines), its third line replaced: run/015 for 5 steps
     with the fault bridge scaling every p_g by 1.0000000000000004
     (`LESTO_TEST_PGAS_SCALE`), a stand-in for another bridge: the verdict
     NOTRUN (the start-up and step 1 IDENTICAL, 59 lines), and FAIL with
     the gas inventory of step 1 changed in the excerpt; the whole excerpt
     against its log is printed (DIFFERENT, 22 lines).
   - Documents: the header of `tests/gems/Allrun`, `tests/Alltest` (the
     usage text, its NOTRUN cases), `faultBridge.c`, the readme (M9.8),
     run/015's readme, `thermochemistry/README.md` (Merlin7); decision e.

W. **The sampler test leaves no orphan (M9.3b, second line).**
   - Defect (review): its control was the endless loop of item K, started
     without a check of its parent: a test interrupted while the control
     ran left the loop running and appending to the work area.
     Reproduced with the first round's test: an INT to its process group
     21 s in, during the control, and the loop kept sampling (10 samples,
     14 four seconds later) until stopped by hand.
   - Fix: the control is bounded to 30 s (it still outlives a killed
     script by more than the 5 s of the check, which is what it shows),
     the scripts wait 60 s instead of 300 s, and the test's own EXIT, INT
     and TERM traps stop the scripts, samplers and controls it started
     (`samplerCleanup`), each by its recorded PID and only while its
     command line names the test's script, so that a PID taken over by
     another process is never touched (the bounded sleeps end by
     themselves).
   - Tests: M9.3b(ii) as before (TERM, KILL and the control: 'sampler
     alive 5 s after its script ended: control'); by hand (one-off), an
     INT and a TERM to the test's process group 21 s in leave none of its
     processes 3 s later.

X. **A value that is not finite fails the gates (`tests/phaseChange/
   checks.py`, `graetzChecks.py`; M9.14, M7.1n, M9.1).**
   - Defects: (verification) the ledger of item B gates with
     `tests/phaseChange/checks.py`, whose closure is a max() over the
     steps, which passes over a NaN that is not the first value: a balance
     file with the closure NaN in its second row passed, so did one with
     the gas inf (the reference inf, the ratio 0), and `logClosure` passed
     a NaN closure in a later Balance line.  The steady check of the
     Graetz runs (M7.1d) took its closure the same way; item O had left it
     out.  (Review) M7.1n compared exit statuses only, and a crash of the
     checker is exit status 1 as well: a checker that crashes at a NaN
     passed it.  (Verification, review) M9.1's run of the bridge build was
     the one run without the ledger, against item B's claim.
   - Fix: `readBalance` refuses a balance file with a value that is not
     finite (SystemExit naming the column, the time and the file, exit
     status 1 without a verdict line; it serves every subcommand that
     reads a balance file and the checkers that import it, among them
     `longChecks.py`, whose cmdRuns reports it as 'cannot be computed');
     `logClosure` and `steady` refuse a closure that is not finite and name
     the run ('not finite: <run>').  So every gate of the ledger refuses
     it: the ledger of the M9 runs, the closure and solverDefect gates of
     M2 to M5, M4.17 (`paperMode/checks.py closureRuns`) and M7.3a
     (`longChecks.py runs`).  M7.1n requires the refusal text, and covers
     steady; M9.1's run is gated by the ledger (printed in its line:
     closure 0, solverDefect 1.1e-16).  Not changed (outside the
     findings): the comparisons of fields and profiles in the other
     checkers aggregate with max() as well, e.g. the mDep_ check of
     `paperMode/checks.py` with its own reader of the balance file.
   - Tests: M9.14 (new; no run, also without a bridge): `ledger` on
     stand-in balance files with the header of the channel: finite met;
     the closure NaN in the second of three rows, the solverDefect NaN in
     the last, the gas inf in the first refused ('Not finite: ...');
     `logClosure` on a log with a NaN closure in its second Balance line
     refused, on its finite twin met.  M7.1n extended: 10 calls (steady
     met and refused besides the 8 of item O), and the control,
     `graetzChecks.py` with `math.isfinite` raising at a NaN, told apart
     in each of the 6 calls with a NaN (exit status 1, but no 'not
     finite').
   - Negative evidence: the first round's `tests/phaseChange/checks.py`
     gives 4 of the 6 verdicts of M9.14 wrong (the closure NaN and the gas
     inf pass the ledger, the NaN log passes logClosure, the solverDefect
     NaN is refused without naming it); its `graetzChecks.py` passes the
     NaN closure of steady (M7.1n: 'steady ... exit 0'); its M7.1n passes
     the crashing checker in all 8 calls.  One-off: the gate of the long
     runs (`longChecks.py runs`, M4.17 and M7.3a) on a stand-in run whose
     balance file has the closure NaN in its second row: the first round
     'OK: ... largest |closure|/reference 1.5e-16', now 'Not finite:
     closure nan at t = 0.002 ...' and 'cannot be computed' (exit 2).

Y. **The solver: two log lines and two comments (`gemsEquilibrium`).**
   - Mode frozen (review): the test of an unusable node was `!positive(p_g)
     || !maxLog10Omega(...)`, which short-circuits, so a node with p_g <= 0
     was logged with the start value of its max log10 Omega, -VGREAT
     ('max log10 Omega -1e+300').  It is now computed first, as evaluate()
     does in mode local, and the line names the real value.
   - The log directory (verification): "../processor0/gemsLog" lies below
     processor0 for rank 0 only, so rank 0 wrote into it and the ranks 1-3
     into subdirectories processor<N> of it, and the master's line ('logs
     in processor<N>/gemsLog') named rank 0's directory only.  The decision
     is collective now (`returnReduceOr`): a directory outside processorN
     on any rank gives every rank its subdirectory processor<N>, and the
     line names the directory of every rank ('logs in
     .../processor0/gemsLogBack/processor<N>'); the default and the cases
     of item D are unchanged.
   - Comments (review): `counts::failed` ('failed calls: ...; a rejected
     result, not used either, is counted per pair, in REJECTED'); the
     schedule of mode local is that of the absolute time index: a run that
     starts from a time whose uniform/time holds an index k without a
     lagged state updates at k + 1 and then at the schedule's next, so its
     first interval may be shorter (the header, the readme and two
     comments of `interfaceExchange` said 'the first step of a fresh run
     and then every updateInterval steps'; those comments keep their line
     count, so no __LINE__ of an error message moves, and the executables
     rebuilt from them are bit for bit those tested).
   - Tests: M9.5(ii) extended: besides GEMSB_ERR_FATAL at the node call 20,
     p_g = 0 at 30 and log10 Omega = NaN at 40: 149 nodes from GEMS3K, 3
     failed, 1 re-creation, the rank log naming the real max log10 Omega
     of the first (1.40542) and nan for the second.  M9.12: a fourth
     directory, "../processor0/gemsLogBack", and `rankLogs` checks that
     the start-up line names the directory of every rank.
   - Negative evidence: the first round's sources built with the bridge,
     under the new `tests/gems/Allrun` (`LESTO_GEMS_TEST_EXE`): 40 PASS,
     2 FAIL, 2 NOTRUN (M9.0 by design, and M9.3b at a load of 23), the
     two FAILs exactly the defects: M9.5(ii) ('the p_g 0 nodes with max
     log10 Omega -1e+300') and M9.12 ('rank 1 writes in
     processor0/gemsLogBack/processor1, the line names
     processor1/gemsLogBack', the same for ranks 2 and 3).
   - These change compiled sources (rhoFixedFlowFoam 92589e36... instead
     of bd7c343f...): the gates of the long runs need runs of the final
     sources again (item Z).  Build: 0 warnings (BUILD, M0.3, the bridge
     build of M9.0, M4.6 on v2606); one-off, a copy with the bridge
     compiles and links against the local v2606 Debug build with 0
     warnings, every source compiled (`returnReduceOr` and the frozen path
     on v2606 too).

Z. **Documents of M7, and the long runs again.**
   - The readme's paragraph of acceptance 3 on the 007 carrier said that
     the final deposit of the temperature model does not depend on dt
     (item R had corrected the plan's twin only): now L1(2 ms, 1 ms)
     2.7e-7 at order 0.90, 1600 times below HKS.
   - Section 33, item 3d (and item R): the largest mixing-cup and
     wall-uptake changes, 7.536e-11 and 8.129e-10, are at Bi 50, halfCell,
     Nr 256 with the fixedValue outlet, not zeroGradient (7.524e-11 and
     8.112e-10); the one-end and both-ends maxima at Bi 1 with the
     fixedValue outlet too (recomputed from the fits of the runs).
   - The readme's row of `graetzRobin`: the reuse keys hold the shell
     functions that set the runs up (item Q).
   - The long runs (item S again, cost).  Item Y changed compiled
     sources, so every recorded run of M4long and M7long ended under
     another key again, and the cheapest run of each family was made again
     with the final sources (rhoFixedFlowFoam 92589e36...), each on 8 ranks
     bound to 8 physical cores of one socket (0-2, 4-6, 8 and 9), one after
     the other, on a workstation loaded by other users' jobs (a 1-minute
     load of 26 to 38): T5_60-distance (M4long; 6000 steps of 10 ms, one
     segment) in 38.1 minutes and temperature-dt4ms (M7long; 4000 steps of
     4 ms) in 27.1; together 65 minutes of 8 ranks, 8.7 core-hours, and
     about 0.7 more for two starts of temperature-dt4ms stopped before
     their first write (the first beside T5_60-distance on the
     hyperthreads of its cores, the second on hyperthreads of the other
     socket that other users' jobs shared, 4 steps per minute; each later
     start removed what the stopped one had left, `runM7`'s path 'prepared
     before').  Both reproduce their records bit for bit:
     temperature-dt4ms every file (1408 field files and 128 uniform/time
     files of the 8 ranks at 16 times, 64 profile files, the 4000 rows of
     the balance file); T5_60-distance every field of every rank at its
     60 write times (5280 files), the 120 profile files and the 6000 rows
     of the balance file (the record's two segments merged), and
     OpenFOAM's clock value in uniform/time up to 4 s (from 5 s on the
     record carries the clock of its restart at 4 s, item S:
     59.9999999999966747 against 59.999999999996632 at 60 s).  The changes
     of the follow-up do not act on these runs, as expected.  Results
     (`tests/Alltest M4long M7long` with the revalidation, below): M4.17
     gates T5_60-distance, 7.9e-16 of the reference over its 6000 steps,
     and lists the 13 revalidated runs; M4.5 for T5_60 with layer distance
     gives T_dep 789.2542 K as before; M7.3a gates temperature-dt4ms,
     closure 2.8e-16 and solverDefect 5.0e-13 of the reference, and lists
     the 5 others; the other 29 runs (13 of M4long, 16 of M7long) were
     revalidated by their 5-step probes, every probe the same files,
     profiles and balance rows, the booking columns included; every
     criterion of M4long (18 lines) and M7long (18) passed.

M9 and M7 after the cleanup: `tests/Alltest -fresh all` (M0 to M5, M7 and
M9; with `LESTO_OPENFOAM_V2606`, bound to the 20 physical cores with
`taskset -c 0-19`) gives 300 PASS, 0 FAIL, 1 NOTRUN (M0.4, the compile on
v2606 on Merlin7: no remote access; M4.6 compiles a copy against the
local v2606) in 41 minutes at a load of 2 to 11 (the M7long run of item S
alongside for 18 of them); M9 itself, 42 lines then, in 2 minutes (M9.3b:
42.9 us per warm call at a load of 2.5), M7 in 10 (the Graetz runs made
again: their key holds the solver binary); 0 warnings; the M0 gate
identical to the tag 0.1.  After the line of item I was added (M9.7, 16
lines), `tests/Alltest M9 M4long M7long` (`M4LONG_REVALIDATE=yes
M4LONG_RUNS=ended M7LONG_REVALIDATE=yes M7LONG_RUNS=ended`, one
invocation) gives 91 PASS, 0 FAIL, 0 NOTRUN in 10 minutes: BUILD, the
M0 gate, SELF and CACHE, the 43 lines of M9, the 18 of M4long and the 18
of M7long (M7.2a-b, M7.6a-b and the 14 of the long runs); of the 31 long
runs of both, 2 were made with the final sources and 29 revalidated.  The
branch is phase-change-rebased at 5852f8d, without commits or stashes.

After the follow-up (items U to Z): `tests/Alltest -fresh all` (with
`LESTO_OPENFOAM_V2606`; bound away from the 8 cores of the M7long run of
item Z for the 22 minutes it ran alongside, then unbound; a 1-minute load
of 18 to 30, mostly other users' jobs) gives 303 PASS, 0 FAIL, 1 NOTRUN
(M0.4) in 49 minutes: the 300 lines of the first round, the line of item
I, M9.14 and LIST.  M9 (44 lines) took about 3 minutes (M9.3b: 64.2 us
per warm call at a load of 18.8), M7 (27 lines) 2, with the Graetz runs
reused, made just before with the final sources (20 minutes on 4 cores at
a load of 30: their key holds the solver binary); 0 warnings; the M0 gate
identical to the tag 0.1; LIST names no criterion recorded but not
expected.  A first start of this run, with the Alltest before LIST, was
stopped in M1 and started again with the final one.  Then `tests/Alltest
M4long M7long` (`M4LONG_REVALIDATE=yes M4LONG_RUNS=ended
M7LONG_REVALIDATE=yes M7LONG_RUNS=ended`, bound to the 8 cores of the
long runs and their hyperthreads) gives 49 PASS, 0 FAIL, 0 NOTRUN in 9
minutes: BUILD, the M0 gate, SELF, CACHE and LIST, the 18 lines of M4long
and the 18 of M7long (item Z).  Both invocations leave the repository
unchanged (the md5 of its 414 files; during the first only this plan and
the readme were edited), and the branch is still phase-change-rebased at
5852f8d, without commits or stashes.
Pending: the researcher's sign-off of this section.

Decisions (signed off by the researcher on 2026-10-04, as proposed):
a. The schedule of mode local by the time index, and its lagged state as
   per-rank restart state (item C).
b. Unknown keys of GEMSCoeffs and gems stop the run (item H).
c. A log directory outside processorN gets a subdirectory per rank (item
   D), rather than a refusal; on every rank when it lies outside on any
   (item Y).
d. No `maxMoleFraction` guard now; the limit recorded (item M).
e. The comparison of run/015's excerpt: exact with the bridge that wrote
   it, its start-up and step 1 with another bridge and the rest NOTRUN
   (item V; it replaces the masked comparison of item J).

---

## 36. M10 (Pb-Bi-I) design (pending the researcher's sign-off)

Section 6 gives M10 one paragraph: several pairs and samples, a balance per
element, BiI3 and Bi tables, an optional operator-split re-speciation, and
the qualification of the GEMS3K guards at trace composition, "5+ d".  The
paper plan asks for more: per-cell GEMS3K for multi-species Pb-Bi-I.
Section 35, item M, deferred a `maxMoleFraction` guard to this
qualification.  This section is the design for that extension.  Nothing is
implemented, and the repository (phase-change-rebased at 058e390) was only
read.

The design has three studies and a synthesis:
- A, the data: a Pb-Bi-I-He dataset, the GEMS3K system, the p_v tables and
  their validation, including Liu et al. (2025).
- B, GEMS3K: its robustness and cost with that system.
- C, the solver: the chemistry along the tube, the solver design and a
  prototype re-speciation kernel.
- The synthesis rechecked the key claims of all three with its own code
  (items 36.2.3 and 36.3.8), and cross-checked them against each other
  (item 36.2.5).

Evidence is archived outside the repository, in
`~/opt/gems-data/m10-design-2026-10-04/` (internal: its Bi-I data derive
from HERACLES, i.e. Barin 1995; item 36.1):

| study | directory | contents | rebuild |
|---|---|---|---|
| A | `thermo/` | `data/`, `PbBiIHe/`, `generator/`, `validation/` (logs, figures) | `make_all.sh`: 13 min serial; a rerun gives identical numbers |
| B | `gems/` | `bridge/` (the extended bridge), `driver/`, `maps/`, `logs/`, `cost/`, `figures/` | `run_all.sh` |
| C | `design/` | `notes.md`, `py/` (equilibrium, plug flow, GEMS3K scans), `kernel/` | — |
| synthesis | `synth/` | `check_thermo.py`, `check_gems.py`, `bench_warm.cpp`, `split_bii.py`, with their outputs | — |

The three studies used the same GEMS3K input files (md5 of DCH, IPM and
DBR: 18cdd3b0..., 95b41906..., 3a9b2484...), so their results refer to one
system.

### 36.1 The data and their licences

1. **The dataset: 19 species in two groups.**
   - **Pb-I-He: 10 NASA-9 records (Gurvich 1991).** He, I, I2, Pb, PbI
     and PbI2 (g); PbI2(cr, l); Pb(cr, l).  These are the records of
     `systems/PbI2He`, unchanged.
   - **Bi-I: 9 ThermoFun records** (`thermo/data/bi_i_m10-thermofun.json`).
     - Bi(g), Bi2(g), BiI(g), BiI3(g), Bi(cr), Bi(l) and BiI(cr) come from
       HERACLES-TDB v0.2.  Its numbers equal Barin (1995).  Their data are
       unchanged.
     - BiI3(cr) is HERACLES's record cut at T_fus = 680.85 K.  HERACLES
       continues the crystal past 681.8 K through a zero-enthalpy
       transition to the liquid heat capacity.
     - BiI3(l) is constructed, since no open record exists:
       - the crystal at T_fus;
       - TKV's fusion enthalpy, 31.798 kJ/mol;
       - Cp(l) = 150.624 J/(mol K), Barin's 36.0 cal/(mol K).  TKV's
         vaporisation enthalpies imply 141.6 +- 15.6, which includes it.
   - **Conventions and choices.**
     - Every record uses G = H - T S from its own H298 and S298.
       HERACLES's own G298 values disagree with that by 1.49 J/mol for
       Bi(cr): 56.735 against 56.74 J/(mol K), so Bi melts at 544.58 K by
       G and at 544.51 K by H and S.  With the one convention, the NASA
       and HERACLES records mix exactly.
     - Pb(g) is NASA's.  HERACLES's Pb(g) carries the heat capacity of
       condensed Pb and is 3.9 kJ/mol off at 1000 K.
     - Molar volumes [m3/mol]: Bi 2.1368e-5 and BiI3 1.0206e-4 (CRC
       densities); BiI(cr) 4.826e-5 (an additive estimate).
   - Source, licence and correction of every species are in
     `thermo/data/provenance.csv`.
2. **Licences.**
   - **NASA CEA thermo.inp is Apache-2.0**: redistributable with
     attribution, as the PbI2 table and the PbI2He system already are.
   - **HERACLES** lies in a GPL-3.0 repository, but its numbers are
     Barin's copyrighted compilation.  Everything derived from it stays
     internal until PSI confirms its provenance and licence: the dataset,
     the GEMS3K system PbBiIHe, the Bi, Bi2, BiI and BiI3 tables, and
     BiI3(l).
   - **TKV** is view-only.  Its values are cited individually and never
     copied.  The primary reference of the BiI3 fusion data (31.8 +- 4.2
     kJ/mol) was not retrieved: a gap.
   - **JANAF, Burcat and HSC** served only as internal checks: PbI2
     deposition temperatures and a BiI(g) test variant.
   - **Consequence:** the repository gets code, generator changes and
     tests, never the Bi data (decision 1).  Tests that need those data
     are NOTRUN without them, as M5.4b is without thermo.inp.
3. **The BiI3 data against TKV.**
   - Vapour pressure, log10(p/Pa): -2.987 at 400 K, 0.440 at 500 K, 2.682
     at 600 K and 3.984 at 680.85 K, as the research report gives them.
   - Triple point: log10 p = 3.9845, against TKV's 0.089 atm (3.9551):
     +0.029 decades.
   - **Normal boiling point: 799-801 K, not TKV's 815 +- 1 K.**
     - BiI3 monomer alone: 800.91 K.
     - Whole congruent vapour (96.4 % BiI3): 798.63 K.
     - GEMS3K: 798.76 K.

     The cause is Barin's sublimation enthalpy at 680.85 K, 124.44 kJ/mol
     against TKV's 117.36.  The liquid's Cp does not matter (120-170
     J/(mol K) give 798.1-802.8 K).  Reaching 815 K would need a fusion
     enthalpy of 39.3 kJ/mol, outside TKV's +-4.2.  The liquid's p_v is
     +0.095 decades high at 815 K.  It matters only where BiI3 condenses
     above 681 K, i.e. at more than 9 kPa (decision 10).
4. **BiI(g): the datum that decides the iodine carrier.**
   - Barin and TKV give 74.6 kJ/mol.
   - Liu et al.'s HSC data behave as if BiI(g) were 29-30 kJ/mol less
     stable (item 36.2.4).
     - Burcat's record (a PM3 estimate, dfH 102.80 kJ/mol) comes close.
     - TKV's own D0(BiI) = 238.5 +- 25 kJ/mol covers the difference.
   - The primary measurements were not read: Cubicciotti (1961), J. Phys.
     Chem. 65:521; Oppermann et al. (1991), Z. anorg. allg. Chem. 601:83;
     Cubicciotti and Keneshea (1959).
   - Item 36.2.5 shows that this one value flips the predicted PbI2/BiI3
     split (decision 2).
5. **BiI(cr) is a data problem.**
   - In these data, 3 BiI(cr) -> 2 Bi + BiI3 (all condensed) has dG = +25
     to +32 kJ/mol up to 1000 K, so BiI(cr) never disproportionates.
     Real BiI decomposes near 603 K into BiI3 and a Bi-rich liquid (Liu).
   - With BiI(cr) allowed, a Bi-rich gas deposits 92-95 % of its iodine as
     BiI(cr) (C, plug flow).  Liu found none.
   - **Off by default; an option for sensitivity runs.**
6. **Gaps.**
   - Missing from every open source: Pb2I4(g), Bi4(g), PbBi(g) and
     BiI2(g).
   - KI is not in the system, but took about 7 % of Liu's iodine (a K
     impurity).
   - Left out by choice: Pb2(g) (only HERACLES has it), PbI3(g), PbI4(g)
     and I2(cr, l).
   - The iodine activity in LBE rests on unpublished activity
     coefficients (Liu, citing Aerts 2017).
   - No melt models: neither a Bi-BiI3 liquid nor an LBE phase with
     dissolved iodine.  The ideal Pb-Bi variant
     (`thermo/validation/variants/PbBiIHe-LBE`) is a test system only, and
     it fails in 1-3 of 80 calls at 870-940 C.
   - The diffusivities of BiI3, BiI, PbI, Pb, Bi, Bi2, I and I2 in He are
     estimates; no measurement was found.
7. **The p_v tables** (`thermo/data/`, in the repository's format):

   | table | phases | melting point [K] | largest interpolation error [decades] |
   |---|---|---|---|
   | `pv_BiI3_phases.csv` | cr, l (198 nodes) | 680.850000 | 6.1e-5 |
   | `pv_Bi_phases.csv` | cr, l | 544.512416 | 1.0e-5 |
   | `pv_Bi2_phases.csv` | over 2 Bi(cr, l) | 544.512416 | 3.1e-5 |
   | `pv_BiI_phases.csv` | one phase | — | 9.6e-6 |
   | `pv_Pb_phases.csv` (extra; NASA, redistributable) | cr, l | 600.650043 | 1.2e-5 |

   - The errors are those of interpolation as the solver does it; the
     limit is 2e-3.
   - Reference files sample every 0.25 K.
   - The BiI table is meaningful only below about 600 K (item 5).
8. **GEMS3K's gas constant** is R = 8.31451 (`ms_multi_format.cpp:36`),
   not the solver's 8.314462618.
   - The difference puts GEMS3K up to 1.0e-4 decades (BiI3) and 1.3e-4
     decades (PbI2) off the tables at 300 K.
   - Its 10 K temperature grid itself is exact to 2e-8 decades.
   - Not compensated: scaling G0 would change PbI2He's PbI2 by up to
     1.3e-4 decades, which the "PbI2 unchanged" criterion forbids
     (decision 11).

### 36.2 The GEMS3K system and its validation

1. **PbBiIHe** (`thermo/PbBiIHe/`).
   - Size: 4 elements (Bi, He, I, Pb), 19 species, 10 phases (an ideal gas
     and 9 pure condensates).
   - Grid and solver settings are those of PbI2He:
     - temperature 270-1250 K every 10 K (99 nodes);
     - pressure 0.8, 1.0 and 1.2 bar;
     - pa_DK 1e-5, pa_DHB 1e-10, pa_PSTALL 0, pa_IIM 7000.
   - The 10 shared species are bit-identical to PbI2He: G0, H0, S0, Cp0,
     V0 and molar masses.
   - The extended generator reads ThermoFun records as well as NASA-9.
     - It agrees with the ThermoFun engine to 1.2e-10 J/mol.
     - It still writes the repository's PbI2He system (11 files) and
       `pv_PbI2_phases.csv` byte for byte; `make_all.sh` checks both on
       every run.
     - Diff: `thermo/generator/make_gems3k_input.diff`.
2. **Validation through the solver's bridge** (A, `thermo/validation/logs/`).
   - BiI3 vapour pressure: the report's table to <= 0.0045 decades.
   - log10 Kp of the Bi-I reactions:
     - HERACLES data alone reproduce the report to <= 0.004.
     - With NASA's I and I2 (the M10 mix):
       - BiI3(g) = BiI(g) + I2(g) differs by <= 0.003;
       - 3 BiI(g) = 2 Bi(l) + BiI3(g) does not change;
       - BiI(g) = Bi(g) + I(g) differs by +0.02, because NASA's I(g) is
         0.40 kJ/mol more stable than Barin's.
     - GEMS3K's gas speciation gives the same K to <= 2.5e-4 decades.
   - PbI2 unchanged, against PbI2He:
     - the solver's method: <= 2.9e-5 decades;
     - the dual-based values: identical to 1.8e-15.
   - ThermoFun mode (-o) against the DCH grid (-t): <= 1.2e-4 decades.  A
     +1000 J/mol change in `-fun.json` moves K by exactly exp(1000/RT), so
     -o really evaluates ThermoFun at run time.
3. **Rechecked here with independent code** (`synth/`).  The code
   integrates the ThermoFun cp_ft records itself and reads NASA-9 directly
   from thermo.inp; none of the studies' code is used.

   **Thermodynamics** (`check_thermo.out`; every value equals A's to the
   digits A gives):
   - BiI3(l) meets BiI3(cr) at 680.85 K with dH = 31798.000 J/mol and
     dG = 0.
   - Triple point: 680.8500 K, log10 p = 3.9845.
   - log10 p_v: -2.9868 at 400 K, 0.4402 at 500 K, 2.6824 at 600 K.
   - Boiling point: 800.910 K (monomer); 798.628 K for the congruent
     vapour (BiI3, BiI, Bi, Bi2, I, I2), with 96.44 % BiI3.
   - log10 Kp [bar]:

     | reaction | 700 K | 900 K | 1000 K | 1200 K |
     |---|---|---|---|---|
     | BiI3(g) = BiI(g) + I2(g) | -4.8772 | -2.3933 | -1.5297 | -0.2419 |
     | BiI(g) = Bi(g) + I(g) | | | -7.7099 | -5.5828 |
     | 3 BiI(g) = 2 Bi(l) + BiI3(g) | +2.8286 | -0.5762 | -1.7476 | |

     The last reaction changes sign at 858.10 K.

   **GEMS3K** (`check_gems.py`, `check_gems.out`).  Five cold calls with
   PbBiIHe through a copy of the repository's bridge (md5 fb67103e...),
   each against an exact ideal-gas equilibrium of the same records:

   | call | T [K] | gas | condensed species | status, iterations | element balance | result |
   |---|---|---|---|---|---|---|
   | A | 500 | BiI3 1e-3 mol in 1 mol He | allowed | OK, 54 | 2e-16 | BiI3(cr) 9.727e-4 mol; p(BiI3) 2.75534 Pa, against p_v = 2.75532 Pa (+3.0e-6 decades); the dual activity +3.9e-5 |
   | B | 500 | the same | suppressed | OK, 29 | 1.4e-16 | every species within 1.1e-4 decades of the exact solution with the solver's R, 3.5e-7 with GEMS3K's; M9's p_eq = p_g/Omega = 2.75557 Pa (+4.0e-5) |
   | C | 1000 | the same | suppressed | OK, 26 | 2.2e-16 | every species within 1.5e-4 decades |
   | D | 900 | BiI3 at 1e-7 | suppressed | OK, 22 | 2e-16 | see below |
   | E | 900 | an LBE-I-like gas, x_I 4e-6 | suppressed | OK, 11 | 2.1e-16 | every species within 2.4e-4 decades |

   In call D, BiI (99 % of the Bi) and I are right to 3e-6 decades.  The
   minor species are off by 4.8e-4 (I2, x_g 1e-10) to 1.5e-3 decades (Bi2,
   x_g 2e-13).  The species of the absent Pb, floored at 1e-15 by the
   bridge, are off by 0.038 decades.
   The study maxima of 36.3 are the worst cases over warm chains.  A
   single cold call is usually much better.

4. **Liu et al. (2025)** (A, `thermo/validation/logs/liu2025*.txt`, the
   digitised figures in `liu_figs/`).
   - **Their inputs.**
     - A He molar volume of 24.465 L/mol (25 C, 1 atm) reproduces all 17
       of their Eq. (3) values: max 0.5 K, rms 0.33 K.  So 100 mL/min for
       3 h is 0.7357 mol He, and 45 mL/min is 0.3311 mol.
     - Their "Calc. T_dep (GEM)" is the temperature at which half of the
       substance has evaporated.  Complete evaporation comes 11-15 K
       higher.
   - **Deposition temperatures with the open data** (50 % points, against
     Liu's GEM):
     - BiI3: +0.7 to +1.3 K, so HSC's BiI3 equals Barin's.
     - PbI2: -5.0 to -6.3 K, because Gurvich's PbI2 is more volatile.
       JANAF's (an internal check only) gives +1.6 to +3.0 K.
   - **Fig. 3 (PbI2):** the gas speciation agrees species by species.  For
     example at 900 C (1e-5 mol, HSC/ours): PbI2 0.351/0.352, I
     1.897/1.932, PbI 0.597/0.573, Pb 0.651/0.681.
   - **Figs. 5, 6 and 9 differ completely, and one record explains it:
     BiI(g).**
     - With Barin's BiI(g), the rms of (ours - HSC)/input is 0.270, 0.317
       and 0.368 for the three figures.
     - With BiI(g) raised by 29 kJ/mol it is 0.021, 0.014 and 0.047.  The
       best shift is 30 kJ/mol for Fig. 5 and 29 kJ/mol for Figs. 6 and 9
       (`liu2025_bii.txt`; A's summary said 30 for all three).
     - So with the open data BiI(g) carries the iodine from LBE above
       about 400 C.  With HSC's data PbI2(g) carries it up to 650 C.
   - **Fig. 9 used 3.0e-6 mol I**, twice the experiment's 1.504e-6 mol
     (x_I 9.04e-4 of 346.0 mg; Table 5: 0.191 mg).  The iodine dissolved
     in LBE in that figure cannot be reproduced (unpublished activity
     coefficients).
5. **The chemistry along the tube** (C, `design/notes.md` section 4; the
   plug flow is equilibrium with fractional condensation, 1 K steps, the
   limit of fast wall kinetics, so it gives onsets, not peaks).
   - **Dissociated in the hot zone, recombined where it deposits.**  At
     x_I = 1e-6 and 1 atm:
     - PbI2(g) carries 0.06 % of the Pb at 1173 K, 26 % at 1000 K and
       99.8 % at 800 K;
     - BiI3(g) carries 4.6e-9 of the Bi at 1000 K, 39 % at 600 K and 96 %
       at 500 K.

     At the deposition temperatures the formula units carry >= 99.9 %
     (PbI2) and >= 99.5 % (BiI3) of their metal.  Below 350 K they carry
     >= 99.9 % at every x_I >= 1e-9.
   - **Given the deposited amounts as the source, two independent pairs
     reproduce Liu.**
     - The onsets of four LBE-I runs and two pure-compound runs lie within
       0-7 K of Liu's Eq. 3.  For example, LBE-I_SS_II gives PbI2 293 C
       (Eq. 3: 297) and BiI3 177 C (176).
     - Measured PbI2 lies within its uncertainty.
     - Measured BiI3 lies 23-40 K below the onset in all five runs with
       BiI3.  Liu's T_dep is the temperature of the peak, not of the
       onset.
   - **The release rate moves the deposition temperatures** by ln(10) R
     T^2/dH_sub per decade of partial pressure: about 40 K for PbI2 and
     30 K for BiI3.  Liu's release profile in time is unknown.
   - **New in this synthesis: the BiI(g) value flips the split.**
     C's plug flow was run again (`synth/split_bii.py`, `split_bii.out`),
     with LBE metal vapour at the boat (ideal LBE, 1173 K, 1 %, 10 % and
     100 % of saturation, BiI(cr) excluded).  The share of iodine that
     deposits as PbI2:

     | BiI(g) | 1 % of saturation | 10 % | 100 % |
     |---|---|---|---|
     | Barin | 5.3 % (BiI3 94.7 %) | 5.3 % | 5.3 % |
     | +29 kJ/mol | 98.1 % (BiI3 1.9 %) | 98.1 % | 98.1 % |
     | +30 kJ/mol | 98.9 % | 98.9 % | 98.9 % |
     | Liu, measured | 13-66 % | | |

     With the shift, the PbI2 onset moves from 537 to 585 K.  With the
     measured split as the source (no metal vapour), the onsets do not
     move: 566/450 K (SS_II) and 500/412 K (SiO2_II) with either BiI(g).
     - So the validation of deposition temperatures is robust to the
       BiI(g) question.
     - The prediction of the split is not: equilibrium with metal vapour
       gives all or nothing, and the measured split lies between the two
       data sets.
     - Without metal vapour the split cannot change at all: each formula
       unit keeps its metal-to-iodine ratio.
     - A split between the extremes needs transport separation and finite
       wall kinetics.  These are what the CFD adds, and a science target
       of M10e, not something to validate against.

### 36.3 GEMS3K: robustness, guards and cost

1. **Maps** (B, `gems/driver/gmap.cpp`, `gems/maps/`).
   - **321,480 calls per run:**
     - 8 element mixes in PbBiIHe, plus PbI2 in PbI2He;
     - x (trace formula units per formula units plus He) from 1e-12 to
       0.32, in quarter decades;
     - 95 temperatures, 304-1244 K;
     - per point one cold call, then 3 warm calls, each on a per-cell
       state with the composition changed by 1e-3 and a new T, as in the
       solver;
     - two modes: (a) condensed species suppressed, M9's mode for HKS
       p_eq; (b) condensed species allowed.
   - **Reference:** every call is compared with the exact equilibrium of
     GEMS3K's own G0, checked by an independent KKT certificate.
   - **Deterministic:** two complete runs gave the same status and
     iteration counts.
2. **Failures, mode (a), at T >= 500 K, by decade of x:**

   | mix | 1e-4 | 1e-3 | 1e-2 | 1e-1 |
   |---|---|---|---|---|
   | BiI3 | 0.1 % | 6 % | 25 % | 42 % |
   | PbI2 + BiI3 (Liu's deposit ratio) | | 3.1 % | 22.6 % | 35.7 % |
   | iodine-rich | | | 8.9 % | 30.4 % |
   | metal-rich (LBE vapour with I) | | | | 73-100 % |
   | PbI2 in PbBiIHe | | | 0.3 % | 1 % |
   | PbI2 in PbI2He (from 1060 K; section 35 M) | | | 0.2 % | 0.7 % |

   - **The Bi-iodide region is new.**  It is hot (>= 750 K) and not
     dilute (x >= 1e-3).
   - **The cause:** the Dikin criterion stalls at 7e-5 to 1e-4, above
     pa_DK 1e-5, while G has converged to 1e-12 relative.  The call then
     runs until pa_IIM is spent: 7000 iterations, about 45 ms.
   - **Warm starts are worse than cold ones there.**  BiI3 at x 1e-2,
     904 K:
     - cold: 27 iterations;
     - warm after a 1e-3 change of composition: 7158 iterations (46 ms);
     - warm after a 1e-6 change: 1 iteration.
   - **The failures are speckled.**  I/Bi = 3.00155 fails, while 3.0015
     converges in 22 iterations.
3. **Mode (b): condensed species allowed.**
   - It adds stalls at trace composition along the saturation boundaries:
     29-51 % of the calls below 500 K for Bi:I 1:1 and metal-rich gas at
     x <= 1e-6.
   - The condensed amounts are silently wrong by 30-99 % of their element
     budget at x <= 1e-7.
   - Liu's source compositions (38 temperatures, 4 calls each):

     | composition | failed calls | remark |
     |---|---|---|
     | the inputs of Figs. 3, 5 and 6 | 3-7 each, all at 304-454 K | |
     | LBE source (Fig. 9) | 0 | warm calls up to 3530 iterations, mean 0.45-0.62 ms |
     | boat cell, BiI3 x = 0.5 | 49 of 152, at 854-1229 K | mean 6.3 ms |

4. **The element balance is never wrong.**  None of the 311,976
   converged calls is off by more than 1e-10.  Every error is either a
   reported failure or silent, so M9's balance check cannot catch the
   silent ones.
5. **Silent errors follow the species' own mole fraction x_g, not the
   cell's x.**
   - Mode (a), largest error over the maps: <= 4e-4 decades for x_g >=
     1e-6; 1e-3 at 1e-7; 0.04 at 1e-8; 0.2 at 1e-9; >= 1 below 1e-11.
     Mode (b) needs x_g >= 1e-4 for <= 8e-4 decades.
   - **M9's p_eq** (primal p_g times 10^-max log10 Omega) is off by:
     - 0.33 decades at x 1e-8 and 1.2 at 1e-9 (PbI2He);
     - inside the guards, up to 1.05 decades for the BiI3 pair of a Bi:I
       1:1 gas.  maxLog10Deviation (0.01) rejects such results, so they
       cost calls but not correctness.
   - **The dual p_eq** (the fugacity from the element potentials) is exact
     to <= 1e-9 decades in every converged call.
   - **chi of a minor pair inside the M9 guards** is off by up to 1.9
     decades (LBE source mix).  M9 never checks chi.
   - **At Liu's dilution** (x_I 1e-7 to 4e-6; C,
     `design/py/gems_trace_scan*.out`):
     - the minor species are off by 0.04-0.64 decades;
     - a tighter pa_DK (1e-8 or 1e-10) makes 33-100 % of the calls fail;
     - a rescaled carrier (He reduced, P' >= 20 Pa) reaches 4e-4 decades,
       but needs a DCH pressure grid down to 1 Pa.  One call with P' below
       the grid but inside its Ptol hung.
   - These are consistent: at x_I = 1e-6 the minor species have x_g of
     1e-10 to 1e-8.
6. **Remedies tried** (B; at T >= 504 K and x >= 1e-7 unless stated):
   - **A warm iteration cap of 200, then a cold retry:** BiI3 warm
     iterations fall from 35.8 to 8.6.  The failures stay.
   - **pa_DK 1e-4:** about 2 iterations per call.  But failures rise to
     30-108 per 5100 calls for the PbI2, iodine-rich and metal-rich mixes,
     and speciation errors reach 1.3 decades.  A found 44 balance failures
     per 1000 warm calls of PbI2 with a condensate.  Rejected.
   - **Polish:** Newton from GEMS3K's result.  It makes 311,966 of 311,976
     converged results exact to <= 3.8e-8 decades, in 2.1 iterations on
     average, at +12.8 us.
   - **An IPM-free exact solver with GEMS3K's G0** (`gemsb_ideal_equilibrate`,
     `gems/bridge/idealgem.hpp`):
     - 0 failures in 321,480 warm-chain calls;
     - 0 failures in 54,720 cold calls with a KKT certificate: mass
       balance <= 1.9e-13, gas chemical potentials <= 1.1e-13, absent
       condensates never supersaturated.
   - **C's gas-only kernel** (`design/kernel/gasSpeciation.H`), at 432
     points: 8.4e-11 from the Python reference, element balance 1.6e-14,
     0 failures.  Rebuilt and rerun here: the same accuracy, and 16.2 us
     cold and 8.9 us warm at a load of 27 (C: 13.2 and 7.6 at 3.8).
7. **Guards, if GEMS3K's interior-point method stays in use** (B,
   `gems/logs/guard_table.txt`).
   - **Keep:** minTemperature 500 K and minMoleFraction 1e-6.
   - **New maxMoleFraction** on the total of the trace formula units.  At
     T >= 500 K and x >= 1e-6, for gas with Bi iodides:

     | maxMoleFraction | failed calls | slow warm calls | mean warm iterations |
     |---|---|---|---|
     | 3.2e-4 | 0 | 0.05 % | 3.4 |
     | 1e-3 | 0.05 % | 0.8 % | 21.6 |
     | 1e-2 | 1.4 % | 6.4 % | 174 |
     | none | 8.5 % | 13 % | 534 |

     - PbI2 fails 0 times up to 1e-2.
     - The metal-rich gas fails 0 times up to 1e-3, and 12 % without the
       guard.
     - **Proposal:** 3e-4 (strict) or 1e-3 with Bi iodides, 1e-2 for Pb-I
       only.  This answers section 35, item M: with Bi the guard is
       needed.
   - **New:** use chi only if x_g >= 1e-7 (mode a).
   - **New:** a warm iteration cap of 200, then a cold retry.
   - **New:** p_eq from the duals.
8. **Cost per warm call** (single thread; B's benchmark, 64 cells at
   600-1000 K, a new T every call, x = 1e-4):

   | engine | PbBiIHe [us] | PbI2He [us] |
   |---|---|---|
   | GEMS3K | 90.8 | 57.2 |
   | GEMS3K, H0/S0/Cp0 grids dropped (lean) | 62.6 | |
   | GEMS3K, one T for all calls | 42.3 | |
   | exact solver | 17.0 | |
   | exact solver, the cell's G0 given | 7.5 | 3.6 |

   - **Rechecked here** (`synth/bench_warm.cpp`): an independent driver
     with only the 30 original functions, on the copy of the repository's
     bridge.  On a physical core whose sibling was idle:
     - GEMS3K: 91.9-93.4 us (PbBiIHe) and 51.5-55.3 us (PbI2He), with
       0 failures and element balance <= 9.6e-11;
     - B's benchmark, rerun: exact solver 17.5 and 7.8 us.
   - **On a hyperthread whose sibling ran another job, the same calls
     took 196.8 and 136.1 us: 2.1-2.5 times as long.**  Production runs
     need whole cores (on Merlin7, `--hint=nomultithread`).
   - **A fit of the cost:** about 80 us plus 6.4 us per interior-point
     iteration.
   - **With Bi iodides at x >= 1e-3,** warm calls cost 1.3-3.5 ms (A).
9. **The cost model** (B, `gems/cost/`): the paper pipe, 60 s at 1 ms
   (60,000 steps), with run/007's temperature field.
   - **Calls per step** come from 228,960 cells, 15,264 wall faces (10,908
     of them above 500 K) and 98,280 cells above 700 K.
   - **CPU-hours:**

     | engine | (a) wall elements | (b) walls + cells > 700 K | (c) every cell |
     |---|---|---|---|
     | GEMS3K as M9, PbI2-like gas | 17 | 173 | 363 |
     | GEMS3K, PbI2 + BiI3 (with failures) | 32 | 319 | 669 |
     | GEMS3K, BiI3-rich gas | 94 | 936 | 1960 |
     | exact solver | 3.2 | 32 | 67 |
     | exact solver, G0 cached per cell | 1.4 | 13.7 | 28.6 |

   - **For comparison,** transport costs 40 CPU-h per paired species and
     360 CPU-h for 9 species.
   - **Checked:** for example (c) = 228,960 x 60,000 x 95.2 us = 363 h.
   - **Re-speciation at C's minTemperature 350 K** covers 212,760 cells
     (93 % of the pipe), so it is close to (c): 27-60 CPU-h with the
     exact solver, 340-620 CPU-h with GEMS3K's IPM.
   - **The model's temperature field reaches only 1000 K**; Liu's furnace
     reaches 1173 K.
   - **State per cell:** GEMS3K 76 doubles (139 MB for every cell of the
     pipe); the exact solver 24 doubles (44 MB) plus a G0 cache of 19
     (35 MB).
10. **A surrogate table buys little** (B, `gems/logs/surrogate.txt`).
    - A 2-D chi(T, log x) table for one pair: 0.03 us per lookup, 1.3e-2
      decades of interpolation error.
    - A 4-D Pb-Bi-I table: 109 MB, 0.79 us per lookup, >= 1e-2 decades.
    - At the stoichiometric ridges, at 500-700 K the minor species change
      10-1000 times across I/Pb = 2(1 +- 1e-4).
    - The saving over the cached exact solver is about 26 CPU-h in (c).
    - **Not recommended** for Pb-Bi-I.

### 36.4 The solver design (recommended options)

1. **Species and condensates** (decision 4).
   - **Three sets:**
     - **S-A:** PbI2_g and BiI3_g, with PbI2_s and BiI3_s.
     - **S-B:** S-A plus PbI_g, Pb_g, BiI_g, Bi_g, I_g and I2_g, gas only
       and re-speciated.
     - **S-C:** S-B plus Bi2_g, which holds 60 % of the Bi of LBE vapour
       at 1173 K, and Pb_s and Bi_s.
   - BiI(cr) is an option, off by default.
   - **New inputs per species:**
     - `formula { Bi 1; I 3; }` in `speciesTransportProperties`, checked
       against `molarMass` to 1e-6;
     - a `ChapmanEnskog` diffusivity model with Lennard-Jones parameters
       for every species but PbI2_g, which keeps `PbI2He`.
2. **Pairs become channels** (C, `design/section36-draft.md` 36.3).
   - **A channel** is a reactant gas r, a condensate c (nu_c mol per mol of
     r) and co-products returned to the gas.  The law stays affine in Y_r:
     q_e = A_e G_r (p_r* - beta Y_r).  So the predictor, the guard, the
     defect and the bounds of sections 2-24 apply unchanged.
   - **Three kinds:**
     - **plain:** p_r* = p_eq(T), as now;
     - **dimer:** Bi2_g -> 2 Bi_s;
     - **reactive:** BiI_g -> 2/3 Bi_s + 1/3 BiI3_g, with p_r* =
       (p_BiI3/K(T))^(1/3), where p_BiI3 is that of the cell at the start
       of the step (lagged).
   - **Rejected alternatives for Bi from BiI:**
     - Bi_g channels after re-speciation are 1e8-1e9 times too slow: p_Bi
       is about 1e-9 Pa at 500 K, while p_BiI is 0.1-1 Pa.
     - Local equilibrium at the wall removes the HKS resistance.
   - **Other rules:**
     - per-pair HKS coefficients (`HKS { accommodation; Ce; }` overrides
       `HKSCoeffs`);
     - the channels into a shared condensate are bounded Gauss-Seidel over
       the species order, so their total evaporation never exceeds the
       deposit;
     - co-products are added at the end of the step;
     - reactive channels are irreversible by default.
3. **Ledgers.**
   - **The pair ledger keeps its format,** so every M2-M9 case and restart
     state still reads and closes.
   - **M10a:** the element ledger (mol of Pb, Bi, I) is the nu/M
     combination of the pair ledgers.
   - **From M10b:**
     - every re-speciated species is solved by the paired path, with no
       active element when it has no condensate;
     - per-species ledgers (REACTION, COPRODUCT, EXCHANGED) and
       per-condensate ledgers live in new lists (`uniform/phaseChangeSpecies`,
       `uniform/phaseChangeElements`);
     - a CONVERSION entry per element books the residue of every transfer
       between stores, so the element closure stays at round-off without
       drift, as solverDefect and clamped do (sections 21.1, 27.1).
4. **Re-speciation: an operator split after the species loop**
   (decision 3).
   - Every cell above minTemperature gets the homogeneous ideal-gas
     equilibrium of its element amounts, with the condensed species
     suppressed; the HKS law does the condensation at the walls.
   - It runs once per step after the loop (a Lie split), never inside a
     corrector loop, so predictor, guard, Final solve and defect are
     unchanged.
   - The split error is first order in dt and scales with 1 - chi of the
     depositing species where it deposits: below 0.1 % (PbI2) and 0.5 %
     (BiI3) at Liu's dilution, more for the reactive channel.  A dt study
     gates it.
   - **Engines:**
     - **`kernel`:** the exact solver in the solver itself,
       `gasSpeciation.H`, standalone and g++-testable like
       `phaseChangeKinetics.H`.  It works in the default build, which
       stays bridge-free, with formation constants log10 K_f(T) from a
       table made by the generator from the same sources.
     - **`GEMS` mode `frozen`:** GEMS3K's G0 per cell, read once at
       start-up (T is frozen) through a new getter (`gemsb_species_g0`,
       9.3 us per cell), and the same kernel every step.  This is what
       "per-cell GEMS3K" means in the recommended scope: GEMS3K's
       thermodynamics in every cell, with an exact minimisation.  It
       equals the per-step call (polished) to round-off and costs what
       `kernel` costs.
     - **`GEMS` mode `local`:** a warm GEMS3K call per cell and step,
       guarded (item 36.3.7) and polished by the kernel, with the kernel
       as the fallback.  This is the paper's "direct warm-started
       GEMS3K", kept for the comparison of accuracy, failures and cost; it
       is not the production engine.
     - The engines differ only in the standard-state data: the open-data
       table and GEMS3K's G0 differ by <= 3.1e-4 decades (C,
       `kf_compare.out`).  In cost, `local` is about 10 times `kernel`.
   - **Conservative:**
     - amounts are exponentials of the converged potentials;
     - absent elements are removed exactly, with no 1e-15 floor;
     - cells with an element amount <= 0 are left alone and counted;
     - a slightly negative Y comes out positive by an element-conserving
       re-speciation, never by a clip.
   - **Defaults:** minTemperature 350 K, updateInterval 1.  The kernel is
     deterministic, so restarts repeat the run bit for bit.
   - **The rule for `speciation lagged`:** re-speciation together with
     HKSCoeffs `speciation lagged` stops the run.  The gas species are
     then speciated already, and chi would count the speciation twice.
5. **HKS p_eq of the pairs.**
   - For an ideal gas over a pure condensate of the same formula, p_eq
     depends on T only.  So the table is exact, and `equilibrium GEMS`
     adds nothing to it (section 32, item 1).  M10 cases use `equilibrium
     table`.
   - If `equilibrium GEMS` is used with Bi present, M10c changes its call
     (decision 6):
     - p_eq comes from the duals, `gemsb_pair_peq`;
     - maxMoleFraction joins the guards;
     - chi needs x_g >= 1e-7;
     - warm calls are capped at 200 iterations, then retried cold.
6. **The source at the boat** (decision 5).
   - **Q1:** prescribed release per species, e.g. Liu's deposited amounts
     (Table 5) as PbI2_g and BiI3_g at the mean rate over the 3 h, which
     is the basis of Eq. 3.  This validates the deposition temperatures.
   - **Q3:** inventory samples of pure condensates for Liu's pure-compound
     runs (mode inventory, as in run/013).
   - **Q2:** elements released as monatomic gases, plus metal vapour at a
     chosen fraction of saturation, speciated by the re-speciation of the
     boat cells.  This predicts the split, or infers it.  It needs M10b
     and M10d, and both BiI(g) data sets (item 36.2.5).
   - **Q4:** equilibrium with liquid LBE.  It needs the iodine activity in
     LBE, which is unpublished.  Not in M10.
   - **Multi-species samples:**
     - one sample may release several species: `amounts { PbI2_g ...;
       BiI3_g ...; }` (one window) replaces `pair` and `amount`, whose
       single-pair form stays valid;
     - inventory samples of different pairs may share cells, with one
       exchange cell per cell holding one sample element per pair;
       `cellCoefficients` assigns Su and Sp per cell, so two exchange
       cells on one cell would overwrite each other;
     - samples of the same pair stay disjoint.
7. **The bridge** (decision 7).
   - **Add the functions M10 needs** from B's 22, which are ABI-compatible
     (the 30 original functions bit-identical over 800 calls):
     - `gemsb_species_g0`;
     - `gemsb_element_potentials`;
     - `gemsb_pair_peq` (dual);
     - `gemsb_gas_mole_fractions`;
     - `gemsb_suppress_condensed`;
     - `gemsb_set_ipm_controls` and `gemsb_get_ipm_controls`;
     - `gemsb_set_warm_iteration_limit`;
     - `gemsb_balance_error`;
     - optionally `gemsb_drop_property_grids` (the lean DCH, -28 us per
       call; nothing in the solver reads the H0, S0 or Cp0 grids).
   - **The bridge's own exact solver** (`gemsb_ideal_equilibrate`, which
     handles condensates) stays in the test build as the oracle of the
     kernel and of mode (b) studies.  The production kernel lives in the
     solver, so the default build needs no bridge.
   - **To update with the bridge:**
     - `tests/build/dummyBridge.c`, which defines every function;
     - the readme;
     - `thermochemistry/README.md` ("30 gemsb_*");
     - `make check`.
8. **Tests without the licensed data.**
   - A synthetic Bi-I data set (invented records with the structure and
     magnitudes of the real ones, documented as such) lets every code gate
     run in the repository: pairs, channels, kernel, ledgers.
   - The HERACLES-based data serve only the science checks and the cases,
     which are NOTRUN without them.  The open Pb-I data stay real.
9. **Restart, parallel and outputs.**
   - **Restart:**
     - the layout records formulas, channels, the engine and its species,
       and a restart with another set stops, naming the change (section
       16, item 5: "this matters from M10");
     - shared condensates keep their deposits per condensate;
     - mode `frozen` recomputes its G0 at start-up;
     - mode `local` warm states are not written, as in M9 (decision 13).
   - **Parallel:**
     - the re-speciation is local to the cell, with one reduction of the
       counters per step;
     - cells above 700 K are 43 % of the pipe, all upstream, so an axial
       decomposition unbalances mode `local` by up to 2.3 times; `kernel`
       and `frozen` cost 8-18 % of the transport of 9 species and matter
       little.
   - **Outputs:**
     - element line densities: the iodine profile is the 126I gamma-scan
       observable;
     - T_peak per condensate, by Liu's definition (the temperature of the
       maximum deposit, from a parabola through three bins), besides the
       onset;
     - deposition-rate profiles;
     - a re-speciation line per step: cells, calls, iterations, failures
       and the largest element residue;
     - the E11 monitor over all ledgered species.
10. **Dictionary** (a sketch; values illustrative):

        // constant/speciesTransportProperties (new keys per species)
        BiI3_g { state gas; molarMass 0.58969381; formula { Bi 1; I 3; }
                 diffusivityModel ChapmanEnskog;
                 ChapmanEnskogCoeffs { sigma 6.0e-10; epsilonByK 800; } }
        BiI3_s { state solid; molarMass 0.58969381; formula { Bi 1; I 3; } }

        // constant/thermochemistryProperties
        model HKS;
        pairs {
          PbI2_g { condensed PbI2_s; vapourPressure { ... } }
          BiI3_g { condensed BiI3_s;
                   vapourPressure { file "<constant>/pv_BiI3_phases.csv";
                     phases ("BiI3(cr)" "BiI3(l)"); units log10Pa; }
                   HKS { accommodation 1; Ce 1; } }          // optional
          Bi2_g  { condensed Bi_s; stoichiometry 2; vapourPressure { ... } } // M10d
        }
        wallReactions {                                           // M10d
          BiI_disproportionation {
            stoichiometry { BiI_g -3; Bi_s 2; BiI3_g 1; }
            equilibrium { file "<constant>/logK_BiI_disp.csv"; column logK; }
            reversible no;
          }
        }
        respeciation {                                            // M10b/c
          engine         kernel;          // none (default) | kernel | GEMS
          species        (PbI2_g PbI_g Pb_g BiI3_g BiI_g Bi_g Bi2_g I_g I2_g);
          formationConstants "<constant>/log10Kf_PbBiI.csv";  // engine kernel
          minTemperature 350;             // [K]
          updateInterval 1;
          // engine GEMS: GEMSCoeffs as in M9, plus GEMS { mode frozen; }
          // (G0 per cell at start-up) or { mode local; maxMoleFraction
          // 3e-4; warmIterationLimit 200; } (a polished call per step)
        }
        samples { boat { selection { ... } mode release;
                  amounts { PbI2_g 2.408e-7; BiI3_g 3.069e-7; }  // [mol]
                  startTime 0; duration 10800; } }
        profiles { ... elements (I Pb Bi); peak yes; rate yes; }

### 36.5 Milestones, cases and effort

Every milestone ends with `tests/Alltest <milestone>` and the M0 gate.  The
code gates use the synthetic Bi data (item 36.4.8); criteria marked (data)
are NOTRUN without the HERACLES-based files.

#### M10a: several pairs of different formulas (3-4 d)
**Deliverables:**
- `formula` per species;
- per-pair HKS coefficients;
- the element ledger derived from the pair ledgers;
- multi-species samples, for release and inventory;
- the Chapman-Enskog diffusivity;
- the extended generator and `make_pv_tables.py` in the repository, with
  no Bi data;
- the BiI3, Bi, Bi2 and Pb tables read from outside the repository;
- element profiles and T_peak;
- run/030.

**Acceptance:**
1. Two pairs against two single-pair runs, each pair's fields, balance
   files and profiles cmp-identical to its single-pair run.
   - Cases: (a) PbI2_g plus an identical pair on a second sample; (b)
     PbI2_g plus BiI3_g.
   - Settings: table and temperature model, serial and 4 ranks.
2. Element closure (Pb, Bi, I) <= 1e-14 of its reference in every step, and
   equal to the nu/M combination of the pair closures to 1e-15.
3. Tables (data):
   - the nodes exact;
   - the crossings at 680.85 +- 0.1 K (BiI3), 544.51 (Bi) and 600.65 (Pb);
   - within 2e-3 decades of the records over 300-1250 K (6.1e-5 measured).
4. A channel from 1000 to 350 K with PbI2 and BiI3 released in the ratio of
   LBE-I_SS_II (x_I 3.9e-6) (data):
   - each condensate's onset (1 %, fromInlet) within 10 K of the plug flow
     (566 K and 450 K), the difference printed;
   - released = n0 M of each species to 1e-14.
5. Inventory samples of two pairs on the same cells: the bounded law of
   both recomputed from the fields; closure <= 1e-13.
6. Every earlier criterion unchanged, the excerpts of run/012-015
   reproduced, and the default build bridge-free.

#### M10b: re-speciation with the kernel (4-5 d)
**Deliverables:**
- `gasSpeciation.H` and `tests/testSpeciation.C`;
- the log10 K_f table and its generator;
- ledgered gas-only species;
- the re-speciation step;
- the species and condensate ledgers;
- the element ledger with CONVERSION;
- the refusal of `speciation lagged` with re-speciation;
- run/031.

**Acceptance:**
1. Kernel (g++):
   - the stored reference points to 1e-9 relative (species above 1e-10 of
     the largest element), element balance <= 1e-14, 0 failures, and a
     second call changes nothing above 1e-15;
   - with a bridge (data): equal to the bridge's KKT-certified exact
     solver over the 160,740 gas-mode calls of B's maps (36.3.1) to 1e-9.
2. Closed box at uniform T (450 and 1000 K), with no wall and no
   transport:
   - any initial split reaches the equilibrium of its element amounts in
     one step;
   - element inventories constant to 1e-15 relative per step;
   - the REACTION entries cancel per element to 1e-15.
3. Equivalences, bit for bit: engine none equals M10a; a case entirely
   below minTemperature equals engine none.
4. The channel of M10a.4 from 1173 K, with the source as formula units:
   - (a) element closure <= 1e-14 in every step;
   - (b) with one diffusivity for every species, the deposits within 1 %
     (L1, 1 cm bins) of M10a's run, since recombination precedes
     deposition;
   - (c) the speciation of every cell above minTemperature equals the
     kernel recomputed from the written fields to 1e-12;
   - (d) dt 4, 2 and 1 ms: L1(2 ms, 1 ms) of the deposits <= 1 %.
5. Restart cmp-identical (binary; serial and 4 ranks); 4 ranks equal serial
   to 1e-9 with a converged loop.
6. Kernel cost <= 20 us per cold call, printed.  NOTRUN on a busy machine,
   as M9.3b; measured 13-16 us cold and 7.6-8.9 us warm.

#### M10c: GEMS3K per cell (4-5 d; 2-3 d without mode local)
**Deliverables:**
- the bridge extension of item 36.4.7, with the dummy bridge and the
  documents;
- `respeciation { engine GEMS; }` in mode frozen (G0 per cell, collective
  as M9's frozen tables) and mode local (warm calls, guards, polish, the
  kernel as fallback);
- for pairs with `equilibrium GEMS`: p_eq from the duals, maxMoleFraction,
  the x_g gate of chi, and the warm cap;
- a reduced qualification map as a test;
- run/032.

**Acceptance:**
1. Bridge:
   - the 30 original functions bit-identical to the old .so (800 calls:
     status, iterations, amounts, log10 activities);
   - `make check` lists the new exports and only libc/libm;
   - the dummy bridge defines every function;
   - 0 warnings.
2. Qualification (data; no solver run), on 3 mixes x 24 T x 13 x in both
   modes:
   - failures and the largest |dlog10| per decade of x_g recorded;
   - the polished result equal to the exact equilibrium of the same G0 to
     1e-9 at every converged point (gate);
   - GEMS3K's K_f within 1e-3 decades of the open-data table (3.1e-4
     measured).
3. Element closure <= 1e-14 in both modes.
4. The channel of M10b.4:
   - frozen within 1e-3 (L1 of the deposits) of engine kernel;
   - local (polished) within 1e-9 of frozen;
   - every failure taken by the kernel and counted.
5. Pairs with `equilibrium GEMS` and Bi present:
   - the dual p_eq within 1e-3 decades of the table at every GEMS element;
   - the maxMoleFraction and x_g guards recomputed from the fields (as
     M9.3a).
6. Cost printed: the start-up of frozen, and the warm call of local with
   its polish, <= 150 us at x <= 1e-4 (NOTRUN on a busy machine).
7. Restart: frozen exact; local within 1e-12 of the continuous run's
   inventories.

#### M10d: metal vapour, shared condensates, the reactive channel (4-5 d)
**Deliverables:**
- per-condensate deposits shared by channels, with Gauss-Seidel bounds;
- the dimer and reactive channels (irreversible) and co-products;
- the Q2 source;
- BiI(cr) as an option;
- run/033.

**Acceptance:**
1. Box, supersaturated Bi and Bi2 over a bare wall at 700 K:
   - the implicit-Euler closed forms of both channels to 1e-14;
   - Bi_s gains the Bi of both.
2. Box, BiI over a bare wall at 500 K: Bi_s gains 2/3 and BiI3_g 1/3 of the
   BiI taken up, to 1e-15 per step, with no evaporation through the
   channel.
3. A small Bi deposit evaporating through Bi_g and Bi2_g:
   - the total <= the deposit;
   - CAPPED to exactly 0, with the residue booked;
   - the species order permuted: the end state within 1e-3.
4. The channel with LBE vapour at 1 % of saturation (data):
   - onsets of Pb(l), Bi(l), PbI2, Bi(cr) and BiI3 within 10 K of the plug
     flow, with both BiI(g) data sets;
   - element closure <= 1e-14.
5. M10a-c unchanged.

#### M10e: the Liu (2025) cases (4-6 d + compute)
**Needs:**
- Liu's measured T(x), digitised or from the authors;
- the carriers of M8 (rhoSimpleFoam, mass-flow inlet), at 45 and 100
  mL/min;
- a 2-D axisymmetric wedge of the column: 1.25 m, R 2.4 mm, about 18.5k
  cells.

**Runs:**
- PbI2_SS and BiI3_SiO2_800 (Q3 or Q1);
- LBE-I_SS_II and LBE-I_SiO2_I with the measured split (Q1, S-A);
- the same with the LBE source (Q2, S-C; metal vapour at 1e-3 to 1e-1 of
  saturation) for both BiI(g) data sets;
- each run over a quasi-steady window at the mean rate.

Then the validation table and run/034.

**Acceptance (code gates):**
- every run ends with element closure <= 1e-12 and solverDefect <= 1e-8 of
  the element reference;
- the table is generated: T_peak per condensate, the iodine profile against
  the gamma scan, the split.

**Science targets** (reported, not gated):
- T_peak of PbI2 within +-20 K of Liu, and of BiI3 within +-40 K (their
  uncertainties are +-9 to +-39 K);
- whether HKS and transport move the BiI3 peak 23-40 K below the
  equilibrium onset, as measured;
- the split under Q2 against 13-66 %, for both BiI(g) data sets.

**Cases.**  015 is the GEMS case, 016 M6's, and 017 onwards are M8's
carriers.  M10 takes a block that cannot collide: 030-pbi2-bii3-pairs,
031-pbbii-respeciation, 032-pbbii-gems, 033-lbe-metal-vapour and
034-liu2025.  Each has a readme and an `out_excerpt`; the Bi data are
linked from outside the repository.

**Effort.**

| scope | milestones | days |
|---|---|---|
| pairs only | M10a | 3-4 |
| recommended | M10a-e, M10c with both modes | 19-25 + compute |
| without the per-step GEMS3K engine | M10a-e, M10c in mode frozen only | 17-23 + compute |

Section 6's "5+ d" covered about M10a.

**Compute:**
- **Liu wedge** (9 species, 60 s at 2 ms): about 16 h serial with the
  kernel and 23 h with mode local, or 2-4 h on 8 ranks (C).
- **Paper pipe** (9 species, 60 s at 1 ms): 25-35 h on 16 ranks if the hot
  cells are balanced.  Merlin7, whole cores.

**External dependencies:** PSI's confirmation of the data (decision 1),
Liu's T(x) (decision 15), and M8's carrier tooling.

### 36.6 Risks

1. **The data.**
   - Bi-I comes from HERACLES = Barin only, with its provenance
     unconfirmed.
   - BiI(g) is uncertain by about 30 kJ/mol, which flips the predicted
     split (36.2.5).
   - BiI3(l) misses TKV's boiling point by 14 K.
   - BiI(cr) is doubtful.
   - Pb2I4, Bi4, PbBi, BiI2 and KI are missing.
2. **The source.**
   - The iodine release from LBE (an unpublished activity) and its time
     profile are unknown, and T_dep moves 30-40 K per decade of rate.
   - The metal vapour level and the BiI(g) data together decide the split.
3. **Local gas-phase equilibrium.**  At x_I about 1e-6 three-body
   recombination takes about 1 s, comparable with the residence time in
   the hot zone.  Not checked: no rate data were read.
4. **GEMS3K's interior-point method,** if it is used:
   - failures of 6-42 % in hot (>= 750 K), non-dilute (x >= 1e-3)
     Bi-iodide gas;
   - silent errors at x_g < 1e-7;
   - erratic convergence in the composition: the M9 guards leave 0.13 %
     failures for BiI3 (at up to 119 ms each); maxMoleFraction 3.2e-4
     removed all of them on the maps, but the failures are speckled;
   - another compiler or GEMS3K build, e.g. on Merlin7, moves the failures
     elsewhere.
5. **The exact solvers are prototypes.**
   - B's (about 630 lines; Levenberg-Marquardt, watchdog and simplex
     steps) was tuned on these failures.
   - C's kernel assumes that no species holds both Pb and Bi: PbBi(g)
     would need a 3-D Newton.
   - Both need regression tests before they enter the code.  They solve an
     ideal gas with pure condensates only; non-ideal phases (an LBE melt)
     still need GEMS3K.
6. **The split error** of the operator split with stiff exchange (lambda
   dt 1e3-6e3) is large for the reactive channel.  It is gated by dt
   studies.
7. **Cost and machines.**
   - The timings come from a 40-core workstation at a load of 27-30:
     - absolute microseconds are uncertain by 20-30 %;
     - a shared hyperthread doubles them;
     - the counts of iterations and failures are deterministic.
   - The pipe with 9 species needs Merlin7 and a balanced decomposition.
8. **Liu's measured T(x)** exists only in figures, and the boat is 3-D.
9. **Not modelled:**
   - homogeneous nucleation in the supersaturated metal vapour (Liu's grey
     deposit is particulate);
   - PbI2-BiI3 and Pb-Bi solid solutions;
   - the transport reaction 2 Bi + BiI3 -> 3 BiI.
10. **Non-dilute gas** near a boat of pure compound: the E11 monitor flags
    it, and GEMS3K fails there (BiI3 x = 0.5: 49 of 152 calls).
11. **The diffusivities** of eight new species are estimates.

### 36.7 Decisions needed (recommendation first)

1. **Bi-I data source and licence.**
   - **Recommended:** task A's dataset (HERACLES = Barin 1995 for Bi-I,
     NASA/Gurvich for Pb-I-He) for development, kept outside the
     repository.
   - Ask PSI, which maintains HERACLES (and Liu's co-authors Neuhausen and
     Eichler), to confirm its provenance and licence before any Bi-I file
     enters the repository.
   - Cite TKV values individually, and find the primary reference of the
     BiI3 fusion data.
   - Alternative: wait for a licensed source before M10a's data
     criteria.
2. **BiI(g).**
   - **Recommended:** Barin = TKV (74.6 kJ/mol) as the default, and an
     HSC-like variant (+29-30 kJ/mol) as a required sensitivity in M10d
     and M10e: two systems and tables from the same generator.
   - Ask PSI which BiI(g) HSC's MainDB holds, and read Cubicciotti (1961)
     and Oppermann et al. (1991), before any prediction of the split is
     published.
3. **Scope against the paper plan's per-cell wish.**
   - **Recommended:** independent pairs plus operator-split re-speciation
     with the exact kernel, and GEMS3K's thermodynamics per cell (engine
     GEMS, mode frozen) as what "per-cell GEMS3K" means in the paper.  The
     per-step interior-point call (mode local) is a polished, guarded
     engine for the comparison only.
   - Cost per 60 s of the pipe, every cell: the kernel 29-67 CPU-h, the
     interior-point call 363-669 (up to 1960 in BiI3-rich gas), against
     360 for the transport of 9 species.  The interior-point call also
     fails in 6-42 % of hot, non-dilute Bi-iodide gas and errs silently at
     trace.
   - Alternatives:
     - (a) pairs only (M10a, 3-4 d): this reproduces the deposition
       temperatures with the measured split, but cannot predict the split;
     - (b) the per-step GEMS3K call as the production engine: 5-23 times
       the kernel's cost, the guards, and no exactness without the
       polish.
4. **Species set.**
   - **Recommended:** S-A in M10a, S-B from M10b and S-C in M10d; BiI(cr)
     off (an option).
   - Not modelled: PbI3, PbI4, Pb2, I2(cr, l), and the species without
     open data (Pb2I4, Bi4, PbBi, BiI2, KI).
5. **Source model at the boat.**
   - **Recommended:**
     - Q1 (the measured amounts per species at the mean rate over 3 h)
       for the validation of the deposition temperatures;
     - Q3 for the pure-compound runs;
     - Q2 (elements plus a metal-vapour fraction) as the prediction and
       sensitivity of the split in M10d and M10e.
   - Q4 (equilibrium with liquid LBE) is not in M10, until the iodine
     activity in LBE is available.
6. **GEMS3K guards and the M9 path with Bi.**
   - **Recommended:**
     - keep minTemperature 500 K and minMoleFraction 1e-6;
     - add maxMoleFraction 3e-4 with Bi iodides (1e-2 for Pb-I only);
     - chi only where x_g >= 1e-7;
     - a warm cap of 200 iterations with a cold retry;
     - p_eq from the duals.
   - This closes section 35, item M.
   - Alternative: 1e-3, which costs 0.05 % failures and 0.8 % slow calls.
7. **The bridge.**
   - **Recommended:** add the functions of item 36.4.7 to
     `thermochemistry/gemsbridge`, with the dummy bridge and the documents.
     The exact solver lives in the solver (`gasSpeciation.H`); the
     bridge's `gemsb_ideal_equilibrate` serves as the test oracle only.
   - The lean DCH (`gemsb_drop_property_grids`, -30 %) is optional.
   - Report upstream to GEMS3K:
     - warm starts that need thousands of iterations where a cold start
       needs 20-40;
     - the Dikin stall above pa_DK;
     - IterDone 0 after some failed calls.
8. **Synthetic Bi data for the code tests.**
   - **Recommended:** yes, so that every code gate runs in the repository;
     the licensed data are for the science checks and the cases only.
9. **Order.**
   - **Recommended:** M10a -> M10b -> M10c -> M10d -> M10e.
   - Alternative: M10a, M10b and M10e (Q1, Q3) first, then M10c and M10d.
10. **BiI3(l).**
    - **Recommended:** keep TKV's fusion enthalpy, 31.8 kJ/mol (boiling
      point 801 K), since it matters only above 9 kPa.
    - Alternative: 39.3 kJ/mol (815 K), which lies outside TKV's +-4.2.
11. **GEMS3K's R** (8.31451).
    - **Recommended:** no compensation (1.3e-4 decades, documented), which
      keeps PbI2He unchanged.
12. **Case numbers.**
    - **Recommended:** run/030-034.
    - Alternative: the next numbers after M8's.
13. **Defaults.**
    - **Recommended:**
      - re-speciation minTemperature 350 K;
      - reactive channels irreversible;
      - the warm states of mode local not persisted, as in M9 (restarts
        exact for kernel and frozen).
14. **Diffusivities.**
    - **Recommended:** Chapman-Enskog with Lennard-Jones estimates, the
      source and values recorded per species.  Which estimates to use is
      the researcher's call.
15. **Liu's inputs.**
    - **Recommended:** ask the authors for the measured T(x) (digitising
      is the fallback), a 2-D axisymmetric wedge, and M8's carrier tooling
      at 45 and 100 mL/min.

### 36.8 M10a implementation and qualification (2026-10-05)

Implemented the independent-pair stage described in 36.5: formulas with
molar-mass validation, per-pair `HKS` overrides, multi-species samples,
shared-cell inventory elements, Chapman-Enskog diffusion, element ledgers
and profiles, and bin-based `Tpeak`. Legacy formula-free cases keep their
existing files and defaults. Restart formulas are compared by element
values, independent of dictionary formatting. Inventory samples retain
explicit element maps after grouping shared cells; every cell is divided
by its volume/density once. Pair `wallPlusGas` profiles use that pair's own
sample mask, while elemental profiles use the union.

`tests/multiSpecies/Allrun` records M10a.0–M10a.6. Independence is checked
byte-for-byte against single-pair runs for both PbI2 plus an identical pair
on a second boat and PbI2 plus BiI3 sharing a boat: temperature/HKS,
serial/4 ranks, fields, ledgers and profiles. Element closure and its
stoichiometric relation to pair closure are checked each step. Shared
inventory tests independently recompute each reservoir's bounded HKS law
from the written fields; a tiny BiI3 reservoir exercises depletion beside
a larger PbI2 reservoir, and restarts cross the scheduled removal.

The generators are code-only additions. The default PbI2He generator
reproduces all 11 generated files byte-for-byte. Locally generated external
BiI3, Bi, Bi2 and Pb tables have 198 nodes: the solver returns the exact
file values at nodes, and crossings are 680.850000, 544.512416 and
600.650043 K. The worst error against direct thermodynamic evaluation is
6.074e-5 decades. The PbI2 table remains byte-identical. Restricted Bi
records and derived tables are excluded from the repository.

`run/030-pbi2-bii3-channel` is a synthetic 1000→350 K helium channel with
the prescribed LBE-I_SS_II iodide ratio scaled to x_I = 3.9e-6. The initial
v2412 serial qualification gives PbI2/BiI3 1% onsets 561.164/449.979 K,
versus independent pure-condensate plug-flow crossings 567.057/450.383 K
(differences -5.892/-0.404 K). Release totals equal n0 M and closures are
within 1e-14. These are code/science qualification inputs, not an
experimental prediction: source re-speciation, measured T(x), physical
carrier/boat geometry, evaporation history and validated Bi diffusivities
remain outstanding. The paper's deposit peak is a different observable.

Final v2412 verification (2026-10-06): `tests/Alltest all` with external
M10 data and the local v2606 installation reports **310 PASS, 0 FAIL,
1 NOTRUN**. All seven M10a criteria pass. The NOTRUN is M0.4, which requires
running M0 under v2606; M4.6 separately compiled every solver and utility
source against local v2606 with zero warnings. The original 008–010 gate,
run/012–015 excerpts, M7 Graetz checks and M9 GEMS checks pass unchanged.
The separate Liu benchmark preparation tests also pass (6 tests).

M10b–M10e, M6 HKS parity and M8 physical/experimental validation remain
subsequent work.

Long-study compatibility verification (2026-10-06): the separate M4/M7
scripts revalidated all 31 completed runs using five-step old/new binary
comparisons. Fields, profiles and checked balance columns are identical;
none of the excluded diagnostic columns differ in these probes. The
scripts report **35 PASS, 0 FAIL, 1 NOTRUN**. M4.17 is NOTRUN because no
entire long trajectory was rerun with this source key; these are explicitly
revalidated historical runs, not fresh full-run closure measurements.

### 36.9 M10b implementation and qualification (2026-10-06)

M10b now implements the bridge-free gas kernel, formation-table generator,
gas-only conservative transport, split update, reaction and condensate
accounts, element CONVERSION, element profiles including all participating
gases, restart records and run/031. The default engine remains `none`.
GEMS standard-state engines and reactive wall channels remain M10c/M10d.

The mixed-precision kernel matches the 432 stored points to 7.1e-11
relative for species above the specified 1e-10 cutoff. Element error is
below 6e-16, as is idempotence on the largest-element scale. Cold calls on
that set are about 15–17 us. The full 160740 external gas-map inputs match
the exact GEMS standard-state oracle to 1.8e-10 relative, with no failures.
For this comparison, the oracle is polished with tolerance 2e-15 (maximum
reported KKT residue 1.7e-14). Its original 1e-13 tolerance left a minor I2
difference of 1.8e-9: oracle stopping error amplified by nearly exhausted
iodine, rather than a kernel element-balance error. The broader map has a
higher mean cold cost (~26–30 us); the <=20 us representative timing gate is
the stored 432-point set, as in the design's benchmark.

Closed boxes at 450 and 1000 K, from monatomic and mixed-iodide initial
splits, conserve elements and cancel reactions within 1e-15. Serial and
four-rank synthetic channels close within 5e-16; their converged fields
agree within 1e-9. Continuous versus binary-restarted fields and new/old
accounts compare byte-identical on both decompositions. Engine `none` and
an entirely sub-threshold case reproduce the disabled gas fields exactly.
Restart comparisons read typed metadata rather than dictionary token
formatting; numerical strings and the content hash are stored as strings
to retain every digit across binary/ascii dictionary I/O.

The external-data 1173–350 K study uses the same prescribed formula-unit
source and common D as the independent-pair baseline. Every step satisfies
1e-14 element closure. At 1 ms, 1 cm deposit L1 differences from that
baseline are 8.5e-7 (PbI2) and 1.333e-3 (BiI3). The 2 ms versus 1 ms
changes are 3.31e-4 and 5.76e-4. Written fields re-speciate within 1e-12
on the largest-element scale. Individual near-exhausted minor species can
show larger relative round-off differences; the same scale is used for
idempotence so species approaching zero are not divided by zero.

`tests/Alltest M10b` and `all` include the new gates. External references,
Bi thermodynamics and the extended oracle are explicitly supplied through
environment variables; unavailable checks are NOTRUN. The synthetic code
tests and the default solver need no licensed data or bridge.

The complete `tests/Alltest all` run under OpenFOAM v2412 finished with
319 PASS, 0 FAIL and 1 NOTRUN, including all nine M10b gates with the
external thermodynamic inputs enabled. The full-map species error was
1.71235e-10, maximum oracle KKT residue 1.59397e-14, element error
5.01064e-16 and idempotence 5.79954e-16. The representative cold-call
timing was 16.2136 us. The 4/2/1 ms study reproduced the deposit differences
above and had maximum element closure 6.25909e-16 across those runs.
Clean v2412 and separate v2606 Debug compilations produced zero warnings.
The sole NOTRUN is M0.4, which requires an invocation under v2606;
the separate v2606 compile passed as M4.6. `all` excludes the hours-long
M4long/M7long studies; their earlier revalidation is recorded in 36.8.
Experimental qualification still requires the physical inputs described
in 36.6 and the later M10 stages.


### 36.10 M10c implementation and qualification (2026-10-06)

M10c implements the optional gas GEMS backend in frozen and local modes,
nine bridge additions, dual HKS pressures for Bi systems, upper mole-fraction
and chi guards, input/restart validation, fault fallback and run/032.
The production library exports 39 C functions. It retains the original
30 signatures and default controls, with no prototype exact-oracle or
polish exports. `make check` gates exports and runtime dependencies
(libc, libm and the dynamic loader). The 800-call comparison with the
M10b bridge is bit-identical in statuses, iterations, amounts, pressures
and activities. Bridge and v2412 solver builds have zero warnings.

Frozen mode caches formation constants at unique cell T/P states per rank.
Local mode makes guarded, capped warm IPM calls; valid duals initialize the
kernel refinement against unfloored physical element amounts. A canonical
cold final refinement removes an otherwise observable restart difference
in nearly exhausted iodine. This additional cost is included in timing.
Fatal calls recreate the engine and invalidate warm states; NaN duals and
other failed calls take the exact kernel fallback and are counted.

The reduced map has 3 mixes x 24 temperatures x 13 dilutions. Both cold
and warm kernel refinements match the independent fixed-volume reference
to 6.82e-10 relative (the original 1e-10 species cutoff); element error is
5.03e-16, with zero failures. GEMS/open-data formation constants differ
by at most 3.69e-4 decades, below 1e-3. Raw IPM failures and errors are
recorded by bulk dilution and species mole-fraction decade. A general
70-digit Decimal Newton solve refines the external exact oracle against
the precise fixed-volume inputs. The double oracle's KKT residue of
1.55e-14 alone left a 1.14e-8 error in a near-exhausted minor iodine gas;
this reference error is reported separately and the 1e-9 gate is retained.

Both gas modes close serial/four-rank channels within 1e-14 and agree
within 1e-9. Frozen binary restarts are byte-identical, including the
accounts; local restart fields satisfy 1e-12 and are exact after the
canonical refinement. Injected fatal and NaN-dual failures are counted,
recreate the engine where required and leave fields byte-identical to
frozen kernel fallback. Invalid modes, guards, caps, carriers, unknown keys
and changed backend settings across restart are refused.

Bi HKS dual pressures agree with the vapour tables to 8.37e-5 decades
across 64 values. Start-of-step fields independently reproduce 2000 guard
checks in each of none/lagged mode, including 36 upper-bound refusals and
326 minor-gas chi refusals in lagged mode. Element closure passes 1e-14.

In the 1173–350 K prescribed-source study at 1 ms, frozen/open-kernel deposit
L1 differences in 1 cm bins are 1.0332e-10 (PbI2) and 1.0637e-7 (BiI3),
below the 1e-3 gates; local/frozen differences are zero, below 1e-9.
Maximum element closure is 6.65796e-16. The mean local IPM plus both
refinements is 132.24 us over 5,814,420 calls, below 150 us. Run/032 records
these provisional numerical benchmarks. M10d shared/reactive condensates
and M10e experimental cases remain, with the physical/data uncertainties
of 36.6 and 36.7 unresolved.

The combined `tests/Alltest M0 M9 M10c` qualification finished with
78 PASS, 0 FAIL and 2 NOTRUN, including all nine M10c gates. M10c.4
explicitly rechecked the completed three-mode prescribed-source study,
including every-step closure and the deposit curves. The other M10c
gates used a fresh optional solver build. A narrow fallback-table check
confirms that in-grid GEMS states never evaluate that table. Final default
and optional v2412 builds and a separate v2606 Debug optional build have
zero warnings; the v2606 executable starts with the bridge compiled in.
M0.4 is NOTRUN because this suite was invoked under v2412. M9.8's
historical run/015 excerpt comparison is NOTRUN because the extended bridge
has another binary hash; the original-API bit comparison passes separately.
Run/032's `qualification.json` preserves raw IPM diagnostics by input-bulk
and individual gas mole-fraction decade; `verification.json` records the
benchmark and acceptance summaries.
Additional checks of all eight saved restart cases confirm species/element
closure within 1e-14 in every output-time directory and final inventories
within 1e-12 of the continuous cases; these checks are included in the
channel test for subsequent runs.


### 36.11 M10d implementation and qualification (2026-10-06)

The opt-in `wallChannels` block adds physical per-condensate stores to the
existing HKS solver. Bi_g and Bi2_g share Bi_s, with two Bi atoms deposited
per dimer; Pb_g uses Pb_s. The irreversible, atom-balanced reaction
`3 BiI_g = 2 Bi_s + BiI3_g` shares Bi_s and uses the lagged gas product
pressure with a dimensionless K (gas standard pressure 1 bar). Products
are applied after all gas transport, before homogeneous re-speciation.
An optional plain BiI_g/BiI_s channel can coexist with the reactive channel.
The existing legacy pairs remain available for PbI2 and BiI3, including
their boat samples. Interface geometry, temperature and half-cell diffusion
resistance are reused. Each channel has its own HKS coefficients.

The species loop uses Gauss-Seidel deposit bounds. The implicit predictor,
bounded guard re-solves and final admissible projection enforce the physical
store limit. CAPPED sets the remaining store to exactly zero. Condensate
accounts persist initial/held mass, compensated exchange, clamping, restart
adjustments and cumulative absolute exchange. The last term defines a
reference even after a store has nearly emptied. Species and element
accounts include Q2 release, shared-store exchange and gas co-products;
element profiles count each physical store once.

`gasSources` supplies unpaired gases in mol over a selected volume and
release window, with exact time-step overlap. Q2 supplies monatomic Pb,
Bi and I. It requires enabled re-speciation accounts. Every channel/source
gas and every gas co-product must be selected. Invalid stoichiometry,
kinetics, keys, table units and incompatible ownership are refused, as are
nonfinite accounts or changes to channel/source configuration and table
contents across restart. Exact local binary stores support byte-identical
same-decomposition restarts. Re-decomposition recovers deposits from mw_
and books the inventory rounding adjustment.

Synthetic gates reproduce the 700 K monomer/dimer implicit-Euler solutions
within 1e-14, the 500 K disproportionation ratios within 1e-15, and the
nonzero-product lagged-pressure solution within 1e-14. Shared evaporation
ends at exactly zero, with the species-order permutation within 1e-3.
BiI(cr), half-cell resistance (including zero diffusion), Q2 source totals,
legacy-sample coexistence and profile integrals are tested. Serial/four-rank
fields agree within 1e-9; fields and accounts restart byte-identically.
Changed-decomposition recovery and every-step species/element/condensate
closure pass 1e-14. A separate external Pb-Bi-I smoke check confirms closure
with both GEMS gas modes and byte-identical final fields between them.

Run/033 uses a provisional Q2 LBE source with metal vapour at 1% of ideal
saturation, activities a_Pb=0.447 and a_Bi=0.553, and Liu's LBE-I_SS_II
deposited PbI2/BiI3 amounts as its iodine basis. The numerical carrier is a
0.1 m plane channel, 200 x 2 cells, 1173–350 K, zero molecular diffusion,
with unit HKS coefficients. The independent fractional plug-flow algorithm
uses the same record evaluator and selected nine gases/eight condensates.
The reference is algorithmically independent of the C++ gas kernel and
wall law. BiI(g)'s +29 kJ/mol sensitivity is generated consistently in its
formation, pressure and reaction tables; it is an HSC-like surrogate, not
an HSC record. Licensed records remain external.

At 1 ms, onset temperatures of Pb(l), Bi(l), PbI2(cr), Bi(cr), BiI3(cr) are
919.93, 882.89, 533.12, 541.35, 450.82 K with baseline data, against
921, 885, 537, 544, 451 K from fractional plug flow. With the shifted BiI
data they are 882.89, 899.35, 582.50, 541.35, 409.67 K, against
882, 902, 585, 544, 408 K. All differences are below 3.9 K (gate 10 K).
The onset is the hottest profile bin above 0.1% of the peak. Measurements
are made during the active source window: a subsequent clean-helium flush
evaporates warm metal deposits and is a different physical comparison.
Run/033 itself closes elements to 4.05e-16, species to 2.21e-16 and
condensates to 2.20e-16 relative. Its README and benchmark metadata provide
the inputs and reproduction commands.

M10e remains: measured T(x), physical 45/100 mL/min carriers, release/split
constraints, the axisymmetric column cases and the validation table/run/034.
M10d's plane carrier and ideal-LBE source are numerical qualification;
the scientific targets of 36.5 M10e still need experimental inputs and the
provenance decisions of 36.7.

The final combined `tests/Alltest M0 M10a M10b M10c M10d` invocation records
43 PASS, 0 FAIL and 13 NOTRUN. All four M10d criteria pass, with their
external inputs enabled. Native M10a–c regressions pass; their earlier
external-data studies were not selected, and M0.4 requires a v2606 suite
invocation. Final default/optional v2412 builds and a separate v2606 Debug
build/startup have zero warnings. Run/033's `verification.json` lists the
scope and skips explicitly and fingerprints the final solver sources.

At 0.5 ms every onset still agrees with plug flow within 3.89 K. Baseline
onsets are unchanged; the shifted-data Pb onset moves one 4.115 K bin.
Maximum element closure over both data variants and both time steps is
4.42e-16. The 1 ms/0.5 ms deposit-profile L1 differences span 0.29–6.45%:
these are reported, not gated by M10d's onset criteria. Quantitative M10e
profiles need additional time-step refinement and the measured inputs.

### 36.12 M10e preliminary benchmarks (2026-10-06)

Run/034 and `thermochemistry/benchmarks/liu2025` now digitise the available
paper curves, solve physical 45/100 mL/min laminar helium carriers, prepare
the axisymmetric Q1/Q2 matrix and report conservation, deposition peaks,
iodine profiles and iodide splits. A bounded finite-volume flux projection
removes the low-Mach carrier solve's finite-precision residue; original
defects and correction sizes remain recorded, and final carrier audits use
the 1e-12 inflow-relative gate. Reference flow conditions of 298.15 K and
101325 Pa are an explicit assumption inferred from the paper's Eq. (3).

Only the supplied paper is available. Three temperature profiles and gamma
histograms were digitised with provenance and pixel uncertainties.
`BiI3_SiO2_800` has no plotted temperature profile and remains NOTRUN.
Figure 7C's temperature trace disagrees with its printed peak labels by
approximately 63 K / 108 K for PbI2 / BiI3. Both observations are retained;
the trace is not fitted to the labels. Coordinates, boat/source placement,
release histories and controller reference conditions remain unconfirmed.

The full 384 x 48 wedge carriers and the coarse 128 x 16 carriers pass their
mass-flow audits. Fifteen coarse cases run 60 s at 20 ms and pass numerical
gates (maximum element closure 6.11e-16, solver defect 1.39e-13, zero gas
re-speciation failures). Five selected full-size cases pass 1 s smoke checks.
Two selected 60 s runs at 10 ms leave iodide peaks unchanged; profile L1
changes are 0.0243%, 3.17% and 2.00%. Fourteen coarse cases meet the selected
1% stationarity diagnostic; silica Q1 does not because of its very small
Pb metal deposit. The full-size 60 s / 2 ms matrix is prepared but NOTRUN.

All five available Q1 peak errors fall inside the proposed 20/40 K targets,
but iodine-profile differences span 47.6–185.1 percentage points across the
matrix. Q1 LBE splits are conditioned on measured deposition. Q2 sensitivity
results use an ideal-activity metal source and the provisional +29 kJ/mol
BiI surrogate. Constant diffusivity 1e-4 m2/s, unit accommodation, omitted
boat blockage, unplotted-end temperature extensions and absent shutdown
modelling remain explicit. These are preliminary numerical observations;
M10e experimental validation and its scientific acceptance remain open.

The full silica smoke revealed an FPE when a positive subnormal product
pressure divided by 1 bar became zero before log(). Log-space conversion
now preserves those pressures and rejects overflow. A closed-form test
with FPE traps covers zero, subnormal and ordinary pressures. Final v2412
M0/M10a–d regressions record 43 PASS, 0 FAIL, 13 NOTRUN, and eleven Python
input/failure-path tests pass. No new v2606 qualification is claimed.
Run/034's README, peak/profile artifacts and `verification.json` preserve
the tested scope, source hashes, missing cases and remaining scientific work.

### 36.13 Full-size M10e matrix launch and dilute-bulk precision (2026-10-06)

The fifteen prepared column cases were launched concurrently on 18,432
cells for 60 simulated seconds at 2 ms. The local runner writes live
progress every 30 s and performs conservation reporting and plotting after
the runs finish. Completion and scientific validation remain pending;
`BiI3_SiO2_800` remains NOTRUN for its missing temperature profile.
Run/034's `full_size_launch.json` fingerprints the private solver, core
sources and case manifests and records the exact tested scope.

The first silica Q1 step exposed relative precision loss at a Pb bulk
concentration of 3.0847e-310 mol/m3: double pressure exponentiation rounded
the smaller pressure before conversion back to concentration. Extremely
dilute bulks now bypass the double predictor, and extreme exponents use
extended libm. All species remain in the equilibrium equations, and the
1e-14 element-residue gate is retained. Twenty-four monatomic-limit tests
include the smallest positive double with FPE traps; the original 140-point
reference is unchanged. The corrected ten-step silica smoke closes elements
to 5.32e-16. The first attempt is archived and all cases started fresh with
one corrected binary.

The v2412 regression session recorded 22 PASS, zero FAIL and nine external
studies NOTRUN, then received SIGTERM during the last scientific half-step
run. That remaining M10d.3 check was completed separately, giving 23
successful numerical/build checks. Both variants and both time steps retain
all five onsets within 10 K and ledgers within 1e-14. Default and optional
builds have zero warnings. The SELF/CACHE/LIST framework checks did not run
after the interrupted session; no completed aggregate Alltest or new v2606
qualification is claimed. The eleven benchmark Python tests also pass.
