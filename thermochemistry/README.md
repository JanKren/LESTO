# Thermochemistry: GEMS3K coupling for rhoFixedFlowFoam

This directory holds the pieces that connect the OpenFOAM solver in
`../applications/rhoFixedFlowFoam` to the GEMS3K Gibbs-energy-minimisation
kernel. Status (2026-10-04): the bridge is wired into the solver as the
optional backend `equilibrium GEMS` of its Hertz-Knudsen-Schrage law
(milestone M9 of `applications/rhoFixedFlowFoam/doc/phase-change-plan.md`,
sections 32 and 35; `gemsEquilibrium.{H,C}`, the readme's section "p_eq
from GEMS3K"). A solver built with `LESTO_GEMSBRIDGE` pointing here calls
it; the default build stays without it. Case `run/015-hks-gems-local` uses
it, and `tests/Alltest M9` tests it.

    gemsbridge/   plain-C interface to GEMS3K (libgemsbridge.so)
    systems/      GEMS3K input generator, open thermodynamic data, PbI2He system
    tests/        OpenFOAM application checking the whole chain

## Why a bridge

- GEMS3K 4.6.1 requires C++20. OpenFOAM v2412 on this workstation is built
  with the system gcc 7.5, which cannot compile it.
- Linking `libGEMS3K.so` directly into an OpenFOAM executable fails, because
  two libstdc++ versions meet in one process:
  `undefined reference to __cxa_call_terminate@CXXABI_1.3.15`.
- `libgemsbridge.so` avoids this:
  - GEMS3K, libstdc++ and libgcc are linked in statically and hidden.
  - Only the 39 `gemsb_*` C functions in `gemsbridge.h` are exported.
  - The library depends on libc and libm only.
- Engines are one per MPI rank; there is no MPI inside the bridge.

## Building

```bash
# once: toolchain (conda-forge only; PSI blocks the Anaconda default channel)
conda create -p ~/opt/gems-env --override-channels -c conda-forge \
    python=3.12 cmake ninja gcc_linux-64 gxx_linux-64 sysroot_linux-64 \
    pybind11 nlohmann_json eigen thermofun spdlog numpy scipy
git clone https://github.com/gemshub/GEMS3K ~/src/GEMS3K   # v4.6.1 tested

cd gemsbridge
./build-gems3k-static.sh     # static GEMS3K: no ThermoFun, hidden symbols
make && make check           # expect "39 gemsb_*, 0 other" and libc/libm plus the dynamic loader

cd ../tests/gemsPbI2Check    # with OpenFOAM v2412 sourced
wmake
cd ../../systems/PbI2He && gemsPbI2Check PbI2He-dat.lst
```

The full GEMS3K build with ThermoFun, and the xGEMS Python API, are installed
in `~/opt/gems-env` (`conda activate ~/opt/gems-env; python -c "import xgems"`).

`libgemsbridge.so` is git-ignored and tied to the machine it was built on:

- It is linked against the glibc 2.34 of the conda-forge sysroot (symbol
  versions up to `GLIBC_2.34`), so it needs glibc >= 2.34 where it runs
  (this workstation has 2.38).
- It carries an RPATH into the local conda environment
  (`/home/jan/opt/gems-env/lib`, `readelf -d`). It needs only libc and libm,
  so the RPATH is harmless where that directory does not exist, but it names
  this machine.
- Elsewhere, e.g. on Merlin7, rebuild it there with
  `build-gems3k-static.sh` and `make` (the same recipe, with that machine's
  conda toolchain), rather than copying it.

## The PbI2–He system (`systems/`)

- **Generator.** `make_gems3k_input.py` writes the GEMS3K input set (the DCH,
  IPM and DBR files plus the `.lst` file) without GEM-Selektor.
  - Model: one ideal-gas phase plus pure condensed phases.
  - G0 comes from the NASA-9 records of `thermo.inp`, which you get with
    `fetch_nasa_thermo.sh`.
- **Data.** Gurvich et al. (1991) as distributed with NASA CEA, Apache-2.0,
  so the data can be redistributed.
  - Species: He, I, I2, Pb, PbI, PbI2 (g); PbI2(cr), PbI2(l), Pb(cr), Pb(l).
  - Cross-checked against NIST-JANAF. The internal-only copy is in
    `~/opt/gems-data`.
- **Validation** (`tests/gemsPbI2Check`, through the bridge, in OpenFOAM):
  - p(PbI2) over the stable condensed phase matches the directly evaluated
    Gurvich p_v(T) to |Δlog10 p| < 1e-4 at 500–1100 K. The phase is PbI2(cr)
    below 683 K and PbI2(l) above.
  - CFD-like loop: 1000 wall cells on an 1100→400 K profile × 50 steps, with
    per-cell warm start. Result: **60 µs/call**, 3.9 IPM iterations per call,
    0 failures.
- **Regenerating:**
  `python3 make_gems3k_input.py --nasa nasa_thermo.inp --out PbI2He --name PbI2He --json --thermofun`

## Findings that matter for the model and the paper

1. **p_v of Lobresco et al. (2026), Fig. 6.**
   - Their "GEMS" curve is the gas–*solid* curve extrapolated through the
     melting point (683 K). Their "Binnewies" curve is the supercooled liquid.
   - Over the stable liquid, p_v is lower by a factor of 2.1 at 823 K, about 4
     at 1000 K and about 7 at 1173 K. That range is the evaporation zone.
2. **JANAF PbI4(g) is spurious.** With condensed Pb present it would
   disproportionate PbI2. Keep it out (the generator's `--pbi4` is off by
   default).
3. **Bi–I data.** The only open source is ThermoHub HERACLES, which is GPL-3
   and derived from Barin.
   - It has no BiI3(l), and its Pb(g) record is defective.
   - A provisional subset is in `~/opt/gems-data`. Confirm its provenance with
     PSI before use.
4. **GEMS3K pitfalls found by testing** (the bridge already handles 1 and 2):
   1. **Ttol trap.** G0(T) is re-interpolated only when T changes by ≥ Ttol
      since the *previous call*. Neighbouring cells with small ΔT then never
      trigger an update: p came out 8.5× wrong. The bridge clamps Ttol to 1e-3 K.
   2. **Tiny amounts.** Element amounts below 1e-17 mol are a hard error. The
      bridge normalises to 1 mol and floors at 1e-15.
   3. **Silent success on stalls.** `pa_PSTALL 1` (new default in 4.6.1) can
      report success with a wrong mass balance at exact PbI2 stoichiometry.
      The generated IPM file uses `pa_PSTALL 0`, `pa_DK 1e-5`, `pa_DHB 1e-10`.
      Always check both the status and the mass balance. Worth reporting
      upstream.
   4. **T and P outside the DCH grid are an error.** The grid is 270–1250 K
      and 0.8–1.2 bar.
   5. **Trace Pb at 300 K** (1e-14 to 1e-10 mol per mol He) does not converge.
      Use the p_v table fallback at the cold end.
5. **Cost.** At Δt = 1 ms, 60 s is 9·10⁸ calls in the 15,264 wall cells:
   about 15 CPU-hours, the same as the transport. Direct coupling in
   interfacial cells is affordable. Surrogates are needed only for calls in
   every cell, or for non-ideal multi-component phases.

## Next steps

1. **Done (M9):** the HKS law of `rhoFixedFlowFoam` takes p_eq from the
   bridge (a call with the condensed phases suppressed gives p_i/Ω), with an
   exact gas–solid mass balance; `equilibrium GEMS`, modes `frozen` and
   `local`, plan sections 32 and 35.
2. **Pb–Bi–I system** (M10). Generate it, re-measure the cost, and reproduce
   the Liu et al. (2025) speciation. At high PbI2 mole fractions (x 0.05 to
   0.73 at 1170-1240 K) GEMS3K does not converge (ERR_NOCONV) and its warm
   calls cost about 4 ms instead of about 50 µs; the solver falls back to the
   table there (plan section 35, item M).
3. **Move the source generator and pv table** into the solver's case data;
   decide HSC vs open data with Prasianakis and Marinich.
4. **Merlin7.** Build the bridge there with gcc 14.2 (same recipe): the
   `.so` of this workstation needs glibc ≥ 2.34 and has an RPATH into its
   conda environment (Building, above). For the tests that is another
   bridge: `tests/Alltest M9` compares run/015's excerpt there in its
   start-up and first step only (no GEMS3K call reaches them) and reports
   the rest NOTRUN, since GEMS3K's results then differ in their last bits
   (plan section 35, item V).

## Preliminary M10 experimental benchmarks

`benchmarks/liu2025/` contains the Liu (2025) experimental tables, a
reproducible Eq. (3) baseline, and four preliminary mesh/source manifests.
The silica and steel geometries are separate; the prescribed LBE sources
retain the reported KI exclusions. The package reports missing inputs and
does not claim CFD validation. See
[`benchmarks/liu2025/README.md`](benchmarks/liu2025/README.md) for preparation,
the flow-reference uncertainty, and comparison of future wall-deposit peaks.

## Licences

- GEMS3K: LGPL-3.0.
- NASA CEA `thermo.inp` and the files derived from it here: Apache-2.0.
- HERACLES (ThermoHub): GPL-3.0. Not included here.
- JANAF-derived files: internal use only, kept in `~/opt/gems-data`.
- Data derived from HSC Chemistry (MainDB): not redistributable.


## M10c gas-engine extensions

The bridge adds nine C functions without changing the original 30 signatures
or their default numerical controls: standard-state G0/(RT), element
potentials, gas mole fractions, dual pair equilibrium pressure, condensate
suppression, element-balance error, IPM control setters/getters and an
engine-local warm iteration cap. A cap of zero keeps the IPM file's limit;
failed capped warm attempts retry cold with the full limit. The production
bridge exports no prototype exact-oracle or polish functions. Those remain
external test tools; the production exact gas kernel is in the solver.

`make check` checks every declared/exported symbol and rejects extra runtime
libraries. The original ABI is qualified on 800 cold/warm calls, comparing
statuses, iterations, amounts, pressures and activities bit for bit.
`tests/Alltest M10c` exercises both GEMS gas modes, guarded/faulted calls,
restart and MPI behavior. The default OpenFOAM build remains bridge-free.
See the solver documentation and run/032 for configuration and external-data
requirements. No licensed Bi-I records are included in the repository.
