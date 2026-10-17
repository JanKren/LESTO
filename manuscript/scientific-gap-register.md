# Scientific submission assessment — 12 October 2026, reconciled 14 October

The 14 October pull adds completed production, scheme/mesh and iteration
results and the Lobresco single-species comparison. The observation-width
audit below remains an analysis of the original upwind campaign. The
original next-study specification and runner archive are retained; the
`current_evidence` field in `next-study.json` records which planned checks
now have saved results. No CFD calculation was launched during reconciliation.

**Status: not yet scientifically cleared for submission.** The current
contribution is a conservative method and conditional sensitivity study.
The new work below strengthens that contribution without treating unresolved
validation as complete. No additional CFD simulations were executed.

## Completed in this assessment

- Recomputed iodine inventory comparisons at 1, 2 and 4 cm observation widths,
  including every integer native-bin offset: 28 comparisons across four cases.
  All edge groups and all deposited inventory are retained.
- Computed the discrete, whole-column, normalised transport distance between
  fine and medium profiles. The 0.0115–0.0711 cm range describes the tested
  binned profiles, not experimental accuracy or continuum uncertainty.
- Recomputed the BiI3 iodine fraction directly from saved species inventories
  at coarse, medium, fine and half-step resolutions. The two thermodynamic
  branches differ by 92.555–92.619 percentage points at matched resolutions.
- Added four analytical tests: edge-inventory preservation, signed-change
  cancellation on aggregation, amplitude invariance of shape measures, and
  a known point-mass displacement. These test post-processing, not the CFD solver.
- Added the methods, results and limits to the manuscript, with a new table
  and machine-readable provenance in `generated/resolution-audit.json`.

## Evidence still required

| Gap | Present evidence | Next action and completion condition |
|---|---|---|
| Spatial resolution | Original upwind fine/medium differences 2.06–10.58%; 4 cm aggregation still reaches 10.02%. The pulled shifted-BiI steel L4 runs give 1.98% limited-linear and 2.13% upwind L3/L4 changes, with inconsistent apparent orders | Extend spatial qualification to other branches and competing maxima. The original four-case refinement specification is partly fulfilled; only the shifted-BiI branch has saved L4 evidence here. Do not claim convergence solely from smaller integrals or a chosen bin origin. |
| Full-run iteration sensitivity | Three cold-start 60 s eight-corrector checks now pass; maximum iodine profile change 0.021%. The shifted-BiI check was still queued when the pulled summary was made | Obtain and audit the remaining shifted-BiI comparison. Preserve the matched three-corrector baseline, mesh, step and physical inputs. |
| Thermodynamic qualification | Internal provenance labels and a surrogate perturbation; strong response persists across tested grids | Obtain complete primary equilibrium/thermochemical records and the exact external input used. Reconcile pressure standard, reaction stoichiometry, reference states, entropy and heat capacity over the simulation range before selecting a branch. |
| Experimental coordinates | Two incompatible silica curve-to-label mappings | Obtain original T(x), scan coordinates and their common physical origin. Retain both hypotheses until registration is independently established. |
| Release and cooling history | Assumed three-hour mean rate; 60 s source runs and fixed-temperature flow-stop diagnostics | Obtain time-dependent release or independent source constraints and T(x,t)/flow histories through shutdown. Full predictive validation needs these or defensible bounds on their effect. |
| Controller and source geometry | Nominal controller convention; incomplete boat geometry | Confirm calibration reference, actual volumetric/mass flow, boat cross-section and axial location. Do not infer boat-flow effects from shifting an empty-tube source. |
| Submission scope and declarations | Working two-author draft for the International Journal of Heat and Mass Transfer (`ijhmt/`, checklist there); nothing submitted | Coauthors must accept the methods/sensitivity scope, scientific evidence and author/funding/declaration details for IJHMT, the only target journal since 17 October 2026. |

The observation-width analysis partly characterises the spatial gap; it does
not close it. The tested thermodynamic effect is numerically persistent; it
does not qualify either thermodynamic dataset.

## Execution dependency found here on 12 October

The accessible workspace has the saved results and source code, but no
`run/024-liu2025/work/prepared` directory. The accessible Ubuntu environment
has no discovered OpenFOAM installation or `rhoFixedFlowFoam` executable.
The original run manifests name a private binary with SHA-256
`16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58`.
Its cases, carrier fields, tables and external Bi thermodynamic coefficients
are required to resume matched studies. The user confirms that the original
environment is on their Linux workstation. SSH host/user and repository path
have been requested; no remote connection or new run has been assumed.
The subsequently pulled saved results provide evidence of workstation runs;
they do not install that execution environment in this workspace.

`workstation_preflight.py --repo /path/to/LESTO` checks the initial
prerequisites without running a solver or modifying cases. Run it with
Python after loading the original OpenFOAM environment. It accepts optional
`--prepared` and `--solver` paths if the retained runs are stored elsewhere.
It does not replace the detailed carrier and external-input audit.

## Primary-source follow-up

The accessible [Liu et al. article](https://link.springer.com/article/10.1007/s10967-025-10246-4)
was checked for supplementary material and data availability links. This
check did not locate the missing original coordinate files or release/cooling
time series. This is a limitation of the available evidence, not a claim
that the authors do not retain these data.

The ACS publisher metadata confirms two relevant primary studies:
[Cubicciotti, 1961](https://pubs.acs.org/doi/10.1021/j100821a030), on the
Bi(l)/BiI3(g)/BiI(g) equilibrium, and
[Cubicciotti, 1968](https://pubs.acs.org/doi/10.1021/ic50060a006), on BiI3
formation enthalpy and BiI dissociation energy. Full articles/coefficients
were not retrieved; DOI or title confirmation does not qualify the records.
No coefficient was replaced using an abstract or inferred value.

## Specific records to obtain from collaborators

1. Thermocouple position/temperature arrays, scan bin edges, common coordinate
   origin and detector response for the three represented configurations.
2. Flow-controller reference temperature/pressure and experiment calibration.
3. Source mass loss or release rate/speciation versus time, residual sample
   inventories, and furnace/column temperature and flow during cooling.
4. Boat dimensions and location, plus the omitted BiI3 800 °C temperature profile.
5. Full BiI equilibrium/thermochemistry records, constructed BiI3 liquid-record
   justification, external input provenance and redistribution permissions.

This is a request checklist, not a message sent to collaborators.

## Reproduce the new analysis

From the repository root with Python and NumPy:

```text
python manuscript/test_resolution.py
python manuscript/assess_resolution.py
```

The normal manuscript build regenerates the reassessment. The existing
checker now also recomputes its 28 aggregation comparisons and four matched
thermodynamic pairs and checks its input hashes. The review archive includes
the new script, tests, numerical outputs and next-study specification.
