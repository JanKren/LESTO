# Lobresco et al. (2026): five PbI2 thermosublimatography experiments

Source: F. Lobresco, G. Tiebel, G. C. Conigliaro, P. Steinegger, O. Marinich,
N. I. Prasianakis, J. Neuhausen, B. Ničeno, *Assessment of temperature-based
and vapor-pressure-driven models for CFD simulation of volatile species
deposition*, Nucl. Eng. Des. **459**, 115201 (2026),
[doi:10.1016/j.nucengdes.2026.115201](https://doi.org/10.1016/j.nucengdes.2026.115201).
The article is open access under the Creative Commons Attribution licence
(CC BY 4.0, <https://creativecommons.org/licenses/by/4.0/>), as stated on its
first page.  The data here are adapted from its Figs. 13 and 14 and Tables 1
and 2 under that licence; retain this attribution.  The PDF itself is not
included.

## What is here

`digitized/`, produced by `digitize_lobresco.py` from the 300 ppi images
embedded in the PDF (no resampling):

| file (per case T106_6, T106_60, T106_300, T5_60, T540_60) | content |
|---|---|
| `<case>_experiment.csv` | measured deposit (mol/m) at the 1 cm scan points, with the light-blue counting band where it is resolvable |
| `<case>_wall_temperature.csv` | the dashed wall-temperature curve (K) of the same panel |
| `<case>_model_temperature_based.csv` | the paper's temperature-based model curve (reference only) |
| `<case>_model_vapour_pressure_based.csv` | the paper's vapour-pressure-based model curve (reference only) |
| `<case>_overlay.png` | digitised points drawn over the original panel, for checking |
| `provenance.json` | PDF and image hashes, axis calibration per panel, Table 1 and Table 2 verbatim, error estimates, checks |

Digitisation error (from `provenance.json`): about 1 mm in position,
2-3 K in wall temperature (plus about 1.5 K/mm of position error in the
deposition zone), 3e-7 mol/m in deposit.  The Fig. 13 wall curve
(T106_60) is hidden by the legend at 0.105-0.135 m and only bridged there.
The paper normalised each measured profile to the sample amount; the
digitised integrals are 0.994-0.998 of it.

Check: the paper's deposition-temperature rule (first position from
upstream where the deposit exceeds 1 % of its peak, linear interpolation,
wall temperature of the same panel) applied to the digitised profiles gives
723.7, 813.0, 667.4, 847.5 and 699.7 K for T106_60, T106_6, T106_300, T5_60
and T540_60, against 723, 814, 667, 847 and 700 K in Table 2 (largest
difference 1.0 K).

## A different wall table

A wall-temperature table distributed with the T-Flows test case of that
study (branch `advanced_scalar_model`, `LESTO_HKS`) does not reproduce the
published per-flow-rate curves.  It is one profile for all cases, peaking at
1210 K at x = 0.04 m, whereas the published curves peak at 1060-1070 K at
0.13-0.15 m.  Between 0.30 and 0.60 m the published curves minus that table
average -23 K at 106 mL/min (largest difference 68-71 K), +2 K at 5 mL/min
(68 K) and +22 K at 540 mL/min (92 K); no single axial shift reconciles them
(best-shift RMS 53-62 K).  The table is not copied here; the comparison runs
use the published curves (`run/023-lobresco2026`).

## Rerun

    pip install numpy scipy pillow matplotlib
    # put the open-access PDF at <repo>/papers/Lobresco2026.pdf, then, in this directory:
    mkdir -p pdf && pdfimages -j -f 9 -l 10 ../../../papers/Lobresco2026.pdf pdf/img
    python3 digitize_lobresco.py
