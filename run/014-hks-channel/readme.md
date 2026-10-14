# Case 014 — HKS deposition in a synthetic plane channel

This case is based on case 013.

The HKS pair of case 013 is moved from the pipe and the 007 development
carrier (228,960 cells, 28 mL/min, a carrier whose T and U were
overwritten) to the small synthetic plane channel that every later
development case uses.  The wall temperature falls linearly from 1000 K to
350 K, so the gas released near the inlet condenses where its partial
pressure meets the vapour pressure, and the deposit forms within a few
centimetres.  The case runs in seconds, and it is the first one with axial
profiles and a deposition-onset temperature T_dep.

## Differences from case 013

The carrier.  Instead of the 007 fields linked into `0/`, the channel is
made by `blockMesh` (`system/blockMeshDict`) and filled by the test helper
`makeCarrier` (`system/carrierDict`):

    0.1 m x 4.8 mm x 1 mm, 100 x 24 cells (1 mm x 0.2 mm), WALL both sides
    Tin 1000; Tout 350; L 0.1;            // K, K, m; T linear in x
    p 101325; molarMass 4.0026e-3;        // Pa, kg/mol (helium)
    massFlux 0.01743;                     // kg/(m2 s), plug flow

What follows from it:
- A 2-D channel has no boat cylinder: the sample is a prescribed release
  of 2.7992536908413373e-11 mol PbI2 over 2 s into the box x = 0.01-0.02 m,
  y = 2.0-2.8 mm (40 cells), the PbI2 amount of case 018.
- Constant `D 1e-4` m2/s, since the PbI2-He correlation needs the p and T
  of a pipe; `formula {Pb 1; I 2;}` adds the element ledger.
- A `profiles` block: native 1 mm bins, `onset {fraction .01; searchFrom
  .02;}`, `observable wall` (the condensate, as the paper defines T_dep).
- Upwind advection, PBiCGStab to 1e-15, dt 1 ms, end time 3 s.
- The p_v table is read from the repository (`$LESTO_REPO/...`), as in
  case 018, so that the two cases differ only in the second pair.

## Modified files

    0/T, 0/U, 0/p                    (templates for makeCarrier, not links;
                                      0/rho and 0/phi removed)
    0/Y_PbI2_g, 0/Y_PbI2_s           (patches of the channel)
    constant/speciesTransportProperties, constant/thermochemistryProperties
    system/blockMeshDict, system/carrierDict, system/decomposeParDict  (new)
    system/controlDict, system/fvSchemes, system/fvSolution
    (removed: the links constant/polyMesh and constant/pv_PbI2_phases.csv,
     out_excerpt, laminarPipeFlow.foam)

## Run

    ./Allrun            # blockMesh, makeCarrier and the 3 s (about 10 s),
                        # then the numbers of expect.json
    ./Allclean

Repository data only: anyone can run it.

## What to check

At 3 s (`postProcessing/axialProfiles/3/PbI2_g_native.dat` and the log):

    T_dep of wall at 0.01 of the peak: fromInlet (from x = 0.02 m)
      561.16445524066103 K at x = 0.067513160732205993 m
      (peak 3.3619636112073234e-09 mol/m at x = 0.0695 m)
    Balance PbI2_g [mol]: gas 1.0644342596433258e-13
      wall 2.7863234206888819e-11 ... released 2.7992536908413376e-11 ...
      closure -3.5044594447401219e-27 (relative -1.2519263460136153e-16)

- The 1 % onset, 561.164 K, lies 5.89 K below the ideal plug-flow
  saturation crossing, 567.057 K (the gate of criterion M10a.6 is 10 K);
  the peak bin is at 548.25 K.
- 99.5 % of the release is on the wall at 3 s; 0.08 % has left through
  the outlet.
- Closures in every step: the pair 2.5e-16 of n0 M, the elements Pb and I
  4.05e-16; the released amount equals n0 to 1e-16.
- The PbI2 ledger, profiles and fields are byte-identical to the PbI2 pair
  of case 018 (two pairs equal two single-pair runs, criterion M10a.1);
  `./Allrun` checks this when case 018 has been run.
