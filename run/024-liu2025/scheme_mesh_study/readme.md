# Spatial scheme and mesh study (plan sections 38 and 39)

Two cases of the M10e matrix, run to 60 s at 2 ms with one change at a time:
the species advection scheme and the mesh level.

| key | case | levels |
|---|---|---|
| `pure` | `PbI2_SS_Q1_baseline` | L1-L3 |
| `steel` | `LBE-I_SS_II_Q2_shift29k_sat0.01` (steel, Q2, 1 %, shifted BiI) | L1-L4 |

Levels (axial x radial cells on the 5-degree wedge, ratio 1.5):
L1 256x32 (8,192), L2 384x48 (18,432), L3 576x72 (41,472), L4 864x108 (93,312).
Schemes for all nine gases: `up` = `Gauss upwind` (the original campaign),
`ll` = `Gauss limitedLinear 1` (the production scheme), `vl` = `Gauss vanLeer`.

## Where the inputs came from

Every case is a byte copy of an existing case with only `system/fvSchemes`
changed: L1 and L3 from the follow-up cases `<case>__grid256x32` and
`<case>__grid576x72` (`work/followup-20261007/cases`), L2 from
`work/prepared/cases/<case>`.  L4 takes the L3 inputs with an 864x108 carrier
built and audited exactly as the follow-up carriers (`liuCarrierTool`, the
same 1e-12 continuity audit); its wall-temperature table is identical to L3's.
Solver: the M10e binary, sha256
`16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58`, unchanged.
The upwind L1-L3 profiles are the existing ones in `../full_analysis/profiles`
and `../followup_results/profiles`; L3 limitedLinear is the same run as the
production case of `../production_576x72_limitedLinear` (identical profiles).

## Files

- `profiles/<key>_L<n>_<scheme>/{element_I,PbI2_g,BiI3_g}_native.dat`: final
  (60 s) wall profiles of the new runs.
- `results.json`: all families (pure up/ll/vl L1-L3, steel up/ll L1-L4,
  steel vl L1-L3): per level totals, peaks (1 cm bin, 3 cm smoothed maximum,
  centroid wall temperature), successive differences (relative L1 of the
  1 cm iodine profile, normalised by the finer level), three-grid status of
  scalar functionals, and the scheme-to-scheme differences at equal mesh.
- `levels1-3_limitedLinear_results.json`, `levels1-3_successive_differences.csv`:
  the first stage (L1-L3, upwind and limitedLinear) before L4 and van Leer.
- `analysis.txt`, `vanLeer_vs_limitedLinear.txt`: readable summaries.
- `upwind_L4/`: upwind L4 runs of all four check bases (pure PbI2 and three
  steel Q2 cases) from the submission package; see its README.

`manuscript/prepare_results.py` reads these files for the numerical table and
the refinement figure; `manuscript/check_results.py` recomputes the successive
differences from the profiles.

## Result in one paragraph

The steel iodine profile changes by 15.35, 8.95 and 1.98 % (limitedLinear) and
10.47, 10.58 and 2.13 % (upwind) between successive levels.  The apparent
orders are inconsistent (1.3 and 3.7 for limitedLinear), so no formal GCI is
claimed.  At equal mesh upwind differs from limitedLinear by 20.8, 14.1, 11.3
and 7.7 % (L1-L4); van Leer differs from limitedLinear by at most 2.8 %.
