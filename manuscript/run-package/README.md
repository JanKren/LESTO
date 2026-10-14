# LESTO scientific follow-up run package

Run this **on the Linux workstation holding the original 79-case study**.
The package never edits the original case directories. It prepares eight
new 60 s cases: four 864×108 (93,312-cell) refinements and four original-mesh
cases with eight species correctors instead of three. Each gets a ten-step
smoke test before the full run. All sources, physical coefficients, time
steps and thermodynamic tables remain matched to the original baselines.

## Start

Extract the archive into a new directory. In a terminal there:

```bash
bash launch.sh --repo "$HOME/LESTO" --carrier-tool /full/path/to/liuCarrierTool
```

Replace the two paths with your actual locations. `--carrier-tool` may be
omitted if `liuCarrierTool` is on PATH or in the prepared study's `bin` folder.
The default original cases are under
`$HOME/LESTO/run/024-liu2025/work/prepared`.
The default NEW output is `$HOME/lesto-submission-followup`.
Use `--prepared /path/to/prepared` or `--out /new/output/path` if needed.
The original solver defaults to `prepared/bin/rhoFixedFlowFoam`;
`--solver /path/to/rhoFixedFlowFoam` overrides its location, not its fingerprint.

The launcher loads `$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc` if needed.
If OpenFOAM is elsewhere, source its v2412 environment first or set
`LESTO_FOAM_BASHRC` to its `etc/bashrc` before invoking the launcher.
If `liuCarrierTool` cannot be found, the package includes its original source.
After loading OpenFOAM v2412, build it with `wmake carrierTool` from this
package directory, then rerun the launcher. Its Make configuration installs
the executable in OpenFOAM's user application directory. The new carriers
must still pass their projection and continuity audits.
Python 3 with NumPy and Matplotlib is required. Check with:

```bash
python3 -c 'import numpy, matplotlib; print("Python dependencies OK")'
```

To keep running after closing the terminal, use:

```bash
nohup bash launch.sh --repo "$HOME/LESTO" --carrier-tool /full/path/to/liuCarrierTool > followup-launch.log 2>&1 < /dev/null &
```

Two cases run concurrently by default; use `--jobs 4` if your workstation
has sufficient spare RAM and CPU. Two carrier preparations run sequentially
before the eight transport jobs. Runtime is not estimated without measuring
this workstation; the refined cases have 2.25 times the previous fine-grid
cell count. No deadline or successful scientific convergence is implied.

## Progress and results

```bash
python3 run_study.py status
```

Pass the same `--out` if you used a custom output directory. During carrier
preparation, read `followup-launch.log` and the new `carriers/*/log.*` files.
During transport, `run_state.json` records smoke/running/pass/fail states;
`cases/CASE/log.solver` contains simulation-time progress.

After completion, send back **`lesto-followup-results.tar.gz` from the new
output directory**. It contains the status, final profiles, numerical reports,
full-log diagnostic summaries, material audits, comparisons against the
original medium/fine profiles, and file fingerprints. It excludes native
executables, thermodynamic tables, raw field directories and full solver logs.
Keep those larger files on the workstation for follow-up inspection.

A `COMPLETE` status means all configured numerical gates passed. It does not
mean the mesh has converged, iteration effects are negligible, experimental
validation is complete, or the paper is automatically ready to submit.
We will assess those questions from the returned evidence.

## Interrupted or failed runs

Read the reported reason first. After fixing the environment, rerun the same
launcher with **`--resume`** and the same paths. Completed passing cases are
checked and kept. Incomplete cases are moved into `interrupted/` and restarted
from time zero using the frozen inputs; they are not resumed from a late
physical checkpoint. A live study lock or a still-running solver prevents
conflicting execution. Do not delete locks or run two launchers for one output.

If interruption occurred during initial preparation, choose a NEW `--out`:
partially prepared studies are not automatically reused. All previous files
remain available. Fix physical or solver discrepancies by investigation,
not by editing fingerprints or bypassing checks.

You can collect finished results even when some cases fail:

```bash
python3 run_study.py collect --out "$HOME/lesto-submission-followup"
```

## Prerequisites and reproducibility

The preflight verifies OpenFOAM v2412, the original solver SHA-256, the
original 60 s/2 ms/18,432-cell configuration, source rates, derived-table
hashes, empty initial reactive inventories and the saved baseline profiles.
The retained table paths in the original manifests must still be valid.
No licensed coefficient files are needed: existing fingerprinted derived
tables are reused. If those tables are absent, stop and restore the originals.

Packaged helper code and digitised inputs are frozen and fingerprinted.
New meshes use independently audited frozen carriers and unchanged source
and transport dictionaries. Sources are not recomputed from a different
thermodynamic installation. Every full-run case starts from time zero.

Local package validation covers Python syntax, shell syntax, analytical
profile checks, and six workflow tests using fake solver calls. **Native
OpenFOAM execution has not been tested on this Windows host.** The mandatory
workstation preflight, carrier audits and smoke tests provide the first native
qualification. Failures remain failures; no numerical tolerance is relaxed.
