# Case 015 — HKS exchange with p_eq from GEMS3K (mode local)

This case is based on case 013 (not 014: GEMS3K is an optional build, and
no later case builds on this one).

Case 013 takes the equilibrium vapour pressure of the HKS law from a table.
Here it comes from a Gibbs-energy minimisation with GEMS3K, for every wall
element and boat cell in every step, at the element's temperature, the
cell's pressure and the cell's composition.  The table stays: it guards
every GEMS3K result and replaces it where a guard applies.  For this ideal
PbI2-He gas the two agree to 1e-5 to 1e-4 decades, so the case shows the
coupling, not a new physical result.

## Differences from case 013

In `constant/thermochemistryProperties`

    HKSCoeffs { ... equilibrium table; ... }

became

    HKSCoeffs { ... equilibrium GEMS; ... }
    GEMSCoeffs { system "<constant>/gems/PbI2He-dat.lst"; mode local;
                 carrier "He(g)"; minMoleFraction 1e-6; minTemperature 500;
                 maxLog10Deviation 0.01; updateInterval 1; ... }

and the pair `PbI2_g` names its species in the GEMS3K system

    gems { gas "PbI2(g)"; condensates ("PbI2(cr)" "PbI2(l)");
           elements { Pb 1; I 2; } excess { I 1e-6; } }

`constant/gems` is a link to the PbI2-He system
`thermochemistry/systems/PbI2He` (Gurvich (1991) via the NASA-9 records).
Below 500 K, below the mole fraction 1e-6, outside the grid of the system,
for a failed call and for a result more than 0.01 decades off the table,
the element takes the table.

## Modified files

    constant/thermochemistryProperties
    constant/gems                       (new, a link into thermochemistry/)
    system/fvSolution                   (comments only)

## Run

The solver must be built with the GEMS3K bridge:

    export LESTO_GEMSBRIDGE=<LESTO>/thermochemistry/gemsbridge
    cd applications/rhoFixedFlowFoam && ./Allwmake

(`libgemsbridge.so` is git-ignored; `thermochemistry/README.md` builds it.)
Then, with `LESTO_SOLVER` pointing to that build:

    ./Allrun            # 5 steps into work/excerpt, against out_excerpt
                        # (about 30 s)
    ./Allrun full 16    # the 16 s case (about 2.3 s per step serial)
    ./Allclean

Without a bridge build `./Allrun` prints NOTRUN.  The whole excerpt is
compared only with the bridge that wrote it (`libgemsbridge.so` md5
fb67103e...); with another bridge only the start-up and step 1 are compared
(no GEMS3K call reaches them) and the rest is NOTRUN, as criterion M9.8
does: the last bits of GEMS3K's p_eq move every round-off token after step
1.

## What to check

One GEMS line per step after the exchange line, with 0 failures:

    Y_PbI2_g GEMS local: update (calls 576: warm 0, retried 0, cold 576;
      ... failures 0, engine re-created 0); elements: GEMS 576, table 14832
      (below minTemperature 4356, below minMoleFraction 10476, ...)

- Step 1 makes no call: 4356 wall faces are below 500 K and no other
  element has gas yet.  Steps 2-5 evaluate the boat's 144 cells and the
  wall cells around it (576 to 720 calls per step).
- The balance line and its closure are those of case 013 to round-off.
- With the bridge of md5 dd8b41fa... (the workstation's current build) the
  start-up and step 1 are identical (59 lines); steps 2-5 are NOTRUN.

Details: notes.md (the GEMS3K system and its guards, the final tolerance,
the output, the excerpt and criterion M9.8).
