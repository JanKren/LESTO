# Case 021 — Metal vapour, shared condensates and reactive wall channels

This case is based on case 019 (not 020: the kernel is enough, and the
GEMS3K engine is optional).

An LBE source releases not only iodides but also metal vapour.  Here the
Q2 source releases monatomic Pb, Bi and I into the helium; the gas
re-speciates as in case 019, and the metals condense through wall channels:
Pb_g on its own Pb_s; Bi_g and Bi2_g on one shared Bi_s (the dimer
deposits two atoms); and BiI_g reacts irreversibly at the wall,
`3 BiI_g -> 2 Bi_s + BiI3_g`, the BiI3 product entering the gas after
every gas has been transported.  PbI2 and BiI3 keep their pairs.  One field
stores each metal across its solid/liquid transition; the stable phase is
the lower of the two pressure curves.

## Differences from case 019

In `constant/thermochemistryProperties`

    wallChannels {Pb {...} Bi {...} Bi2 {reactant Bi2_g; condensed Bi_s;
                  stoichiometry 2; ...}
                  BiI {reactant BiI_g; condensed Bi_s; stoichiometry
                  {BiI_g -3; Bi_s 2; BiI3_g 1;} ... reversible no;}}
    gasSources {boat {amounts {...} box ((.01 0 -1) (.02 .0048 1));
                startTime 0; duration 2;}}

The iodine amount follows the deposited PbI2/BiI3 of Liu's LBE-I_SS_II;
the metal vapour is 1 % of the ideal-LBE saturation pressure at 1173 K
(`benchmark.json`).

What follows from it: the condensates Pb_s and Bi_s; the tables of all
pairs and channels from `LESTO_M10D_DATA`; and, for the comparison with
fractional plug flow, the plug-flow-limit carrier (200 x 2 cells, D = 0)
read at 1.2 s, while the source is still active.

## Modified files

    constant/thermochemistryProperties
    constant/speciesTransportProperties   (Pb_s, Bi_s; D 0 for every gas)
    0/Y_Pb_s, 0/Y_Bi_s                    (new)
    0/Y_Pb_g, 0/Y_Bi_g, 0/Y_I_g           (zeroGradient inlet)
    system/blockMeshDict                  (200 x 2 cells)
    system/controlDict                    (1.2 s, binary, outer tolerance
                                           1e-16)
    benchmark.json, verification.json, out_excerpt (the qualification record)

## Run

    python3 ../../thermochemistry/systems/make_wall_tables.py \
      --nasa /external/nasa_thermo.inp --tf /external/bi_i_m10-thermofun.json \
      --out /external/m10d-baseline
    export LESTO_M10D_DATA=/external/m10d-baseline
    ./Allrun            # 1.2 s at 1 ms (about 30-40 s), then expect.json
    ./Allrun halfstep   # the same at 0.5 ms; the deposit L1 against 1 ms
    ./Allclean

`--bii-shift 29000` makes the consistent HSC-like BiI sensitivity tables
(not a licensed record).  `tests/Alltest M10d` runs the full qualification
against an independent constrained-equilibrium reference
(`LESTO_M10D_NASA`, `LESTO_M10D_THERMOFUN`, `LESTO_M10D_REFERENCE`).

## What to check

The onset of each condensate is the hottest 1 mm bin above 0.1 % of its
peak, the metals split at their melting points (`verification.json`,
baseline):

| phase | fractional plug flow | case 021 |
|---|---|---|
| Pb(l) | 921 K | 919.93 K |
| Bi(l) | 885 K | 882.89 K |
| PbI2(cr) | 537 K | 533.12 K |
| Bi(cr) | 544 K | 541.35 K |
| BiI3(cr) | 451 K | 450.82 K |

- Every onset within 3.89 K of plug flow; the half step changes none of
  them.
- Closures of the element, species and channel-condensate ledgers at most
  4.42e-16 over the four runs of the record (3.7-4.6e-16 measured on
  2026-10-14 with the M10e solver; every onset reproduced exactly).
- The deposit L1 of the half step, 0.29-1.44 % (baseline; up to 6.45 %
  with the BiI variant), is reported, not gated.

A clean-helium flush after the source would evaporate the hot deposits and
move their onsets: that is a different physical comparison.  Measured T(x),
carrier flow and release history are needed for experimental validation
(case 024).
