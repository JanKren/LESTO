# Open questions for the T-Flows authors, and data to request

Collected on 2026-09-23 while comparing Lobresco et al. (2026) with the
OpenFOAM repo and the T-Flows reference code (`/home/jan/src/T-Flows`, branch
`origin/advanced_scalar_model`, `Tests/Laminar/Scalar_Transport_Gpu/LESTO_HKS`,
commit 35c1f6406e, 2026-01-25). That branch may not be the exact version
behind the published figures, so every point here is a question, not a claim.

## A. Carrier-gas kinematics

- `User_Mod/Initialize_Variables.f90:24` sets
  `u_max = 2*FLOW_RATE/(pi*R_PIPE**2)`, with `FLOW_RATE = 1.67e-6 m3/s`
  (100 mL/min at standard conditions). The parabola is uniform along x.
- `User_Mod/Beginning_Of_Iteration.fpp:41` sets `density = 1.0`.
- The paper's Fig. 9 agrees: the gas plateaus are flat over 0.2–0.5 m, and
  n'·t_ev = 3.38e-4 mol·s/m for every t_ev gives U_b ≈ 0.093 m/s.

**Question:** was the velocity deliberately left non-expanding?

**Why it matters:** the physical velocity is larger by T/273.15 (2.6–3.9x in
the deposition zone), and the concentration and partial pressure are smaller
by the same factor. A 1-D equilibrium estimate with physical kinematics gives
T_dep = 807/729/678 K for T106_6/60/300, against the measured 814/723/667 K
(RMS 8 K). The published vapour-pressure model gives 841/760/701 K (RMS 33 K).

## B. Units in the HKS source (`User_Mod/Source.fpp`, `Types.f90`)

- `GEMS_VAPOR_PRESSURE` is in **bar** (`Types.f90:87`).
- `partial_pressure = R_GAS*T*c/1e5` converts to bar.
- The linearisation coefficient is `a_v*k_HKS*(R_GAS*T/1e5)`.
- `k_HKS = 2*SIGMA_C/(2-SIGMA_C)/sqrt(2*PI*R_GAS*MM_PbI2*T)`, with
  `MM_PbI2 = 461.01` in **g/mol**.

Relative to Eq. (6)–(7) in SI units, these conventions scale the rate by
1e-5 × 1/√1000.

- `a_v = a_v_wall*200` in the wall layer and `a_v = 2/R_SAMPLE*5` at the
  sample. Neither multiplier is in the paper. The branch readme reportedly
  says 500.

A 1-D fit to Fig. 11 needs a rate factor of about 1e-5 relative to SI.

**Questions:**
1. Which unit convention and which multipliers produced Figs 11–14 and Table 2?
2. Is C_e the ×200?

**Why it matters:** in consistent SI units the HKS wall conductance
(about 90 m/s) is about 1000x the diffusive supply (k_m ≈ 0.1 m/s). Deposition
is then transport-limited and C_e has almost no effect. The paper's conclusion
that condensation kinetics are "slower than in the experiment" may come from
the units rather than the physics.

## C. Interaction layer

`Initialize_Variables.f90:116` uses a cell-centre test, `R - r <= DELTA`.
On the 228,960-cell mesh that selects **two** layers (30,528 cells,
A/V ≈ 8,545 1/m). The paper's text says control volumes lying *entirely
within* 0.1 mm, which is **one** layer (A_I/V_I = 17,682 1/m).

**Question:** which rule was used for the figures?

## D. Temperature-based model settings

Figs 9 and 10 (the same reference case: A = 220 1/s, k = 8e-3 1/K,
T_dep = 680 K) do not agree with each other:

- Fig. 9 has a solid spike at x ≈ 0.502 m.
- Fig. 10 rises from 0.502 m to its peak at 0.512 m.
- Both start where T ≈ 650 K, but T_dep = 680 K lies at x ≈ 0.48 m.
- Eq. (5) goes to zero at onset, so it cannot produce Fig. 9's spike.

**Questions:**
1. Were different settings or meshes used? `pipe.geo` has L = 1.06 m and a
   "2.5 mm" spacing comment, which hints at a full-length temperature-model
   mesh.
2. The thesis branch (`Lesto_Thesis/User_Mod/Types.f90`) has A = 150,
   k = 0.02 and U_MAX = 0.176. Which values are final?

## E. Other points to confirm

- Which C_e was used for Table 2 and Figs 13–14? The legends say only
  "Vapor-Pressure-Based".
- Does C_e act on sample evaporation? The identical Fig. 11 plateaus suggest
  it does not.
- The sample cell set, its volume, and its A/V.
- The species boundary conditions at x = 0 (back-diffusion in T5_60) and at
  the outlet.
- The GEMS p_v(T) of Fig. 6 is a gas–solid curve with no liquid PbI2 (it
  reaches 22 bar at 1200 K). Is that intended above the ~683 K melting point?
- The Fig. 4 dissociation of PbI2(g) above 1000 K is not in the transport model.
- Brandt correlation: why N = 16 electrons for PbI2 (valence count is 18)?
  Why was the 5-coefficient Ω fit (8–9 % below Neufeld at 666–1100 K) used
  instead of the tabulated values?

## F. Data to request

- [ ] Thermocouple wall-T tables for 5, 106 and 540 mL/min, and for the Liu
      reference case. The repo has only `LESTO_HKS/wall_x_profile_complete.dat`.
- [ ] Measured deposition histograms (1-cm bins) for T106_6/60/300, T5_60,
      T540_60 and Liu PbI2 SS.
- [ ] The GEMS outputs behind Figs 5–6, and the GEM-Selektor project for
      PbI2 + He.
- [ ] Licence status of the HSC-derived MainDB 9.7 for distributing GEMS3K
      input files.
- [ ] The Liu et al. (2025) semi-empirical T_dep formula.
- [ ] Francesco's complete case directories (the repo's `run/001` symlinks to
      `/home/niceno/.../LESTO_Paper` are broken).

## G. Repo fixes found along the way (not paper content)

- `run/007/261/{T,U}` were overwritten by the setExprFields initial fields in
  commit 3ce3791. Restore them with
  `git checkout c187f1e -- run/007-variable-properties-high-temperature/261/{T,U}`
  and use `setExprFields -time 0` in the readmes.
- Cases 007–010 carry 28.3 mL/min (standard), not a paper flow rate, and use
  a linear 1000→300 K wall instead of the measured profiles.
- The mock thermochemistry creates 3.5x the injected PbI2 and has no gas sink.
  No mass-balance monitoring exists.
- `.gitignore` misses the decimal time directories `0.4 … 9.6`. Frozen fields
  are rewritten every output (~36 MB per write).
- The readmes are stale (laminarPipeFlow names; "no solid species"). The wall
  mask looks up the patch by the name `WALL`, whose type is `patch`, not
  `wall`.
