# Case 016 — Half-cell diffusion resistance at the wall

This case is based on case 014.

The HKS law of case 014 takes the partial pressure of the wall cell, at
the cell centre, as the pressure that meets the wall.  On a coarse wall
cell the gas must still diffuse across half a cell to reach the wall.
`wallResistance halfCell` puts that diffusion in series with the kinetic
conductance of the law.  It is the wall law decided for every physical run
(plan sections 28.6 and 34.2: cases 017, 023 and 024 use it).

## Differences from case 014

In `constant/thermochemistryProperties`

    interface {patches (WALL);}

became

    interface {patches (WALL); wallResistance halfCell;}

The defaults it needs (`layer faceCells`, `areaPerVolume cell`,
`areaMultiplier 1`) are those of case 014.  Nothing else changes.

On this mesh the half-cell conductance is D/(dy/2) = 1e-4/1e-4 = 1 m/s,
64 to 107 times smaller than the HKS conductance (63.7-107 m/s on the
ramp), so it controls the uptake at the wall.

## Modified files

    constant/thermochemistryProperties

## Run

    ./Allrun            # blockMesh, makeCarrier and the 3 s (about 10 s),
                        # then the numbers of expect.json
    ./Allrun-damkohler  # the Damkoehler sweep of the manuscript: alpha
                        # = 1 ... 1e-6, 13 runs (about 2 min); with 48,
                        # the wall-normal refinement (5 runs).  Records
                        # and figure: manuscript/data/damkohler
    ./Allclean

## What to check

At start-up the interface line names the resistance, and the conductance
of the wall elements drops from 63.7-107 m/s (case 014) to

    areaPerVolume cell, areaMultiplier 1, interfaceTemperature cell,
      wallResistance halfCell
    Pair PbI2_g, HKS on the interface: ... G_e beta_c/rho_c
      0.98454069006370559-0.9907388780457359 m/s (with the half cell);
      lambda = ... 4922.7-4953.7 1/s

(case 014: lambda 3.2e5-5.3e5 1/s).

At 3 s, against case 014 (measured with the same solver binary):

| quantity | case 014 | case 016 |
|---|---|---|
| 1 % onset T_dep (from x = 0.02 m) | 561.1645 K at 67.5132 mm | 561.1593 K at 67.5139 mm |
| peak bin | 548.25 K, 3.362e-9 mol/m at 69.5 mm | 541.75 K, 3.113e-9 mol/m at 70.5 mm |
| wall inventory [mol] | 2.78632e-11 | 2.78499e-11 |
| left through the outlet [mol] | 2.29e-14 | 3.76e-14 |
| pair closure, every step | 2.5e-16 | 2.5e-16 |

- The onset moves by 5 mK, far inside one 1 mm bin (6.5 K per bin on this
  ramp): where the deposit starts, the gas reaches the wall by diffusion
  either way.
- The resistance slows the uptake, so the deposit spreads downstream: the
  peak is 7.4 % lower and one bin further, and the deposit profile differs
  from case 014's by 6.6 % (L1 of the wall profile at 3 s; 6.5 % at 1 and
  2 s).
- The closures stay at round-off; the released amount equals n0.

The analytic checks of the half-cell law are in the tests (M2.1d, M5.1e,
M7.2f).  On the pipe of case 013 the effect would only show after the
deposit has formed, i.e. in the 16 s run (about 9 h serial).
