#!/usr/bin/env python3
"""Physical carriers and cases of Lobresco et al. (2026) in a work directory.

Runs the tooling of the Lobresco study (case 023: make_cases.py, analyse.py)
unchanged; only its working directory is redirected to <work>, so nothing
is written into case 023.  Used by the Allrun of case 017 (T106_6) and of
case 023 (any of the five experiments).

  physicalCarrier.py prepare <lobresco dir> <work> <level>
      case 017: the carrier of T106_6 (blockMesh, liuCarrierTool initialise,
      rhoSimpleFoam, flux projection, audit) in <work>/carriers/He106 and the
      transport case <work>/case; at level coarse its dictionaries must be
      identical to the tracked ones of case 017
  physicalCarrier.py prepare <lobresco dir> <work> <level> <case> ...
      case 023: the carriers of the flows needed and <work>/cases/<case>
  physicalCarrier.py analyse <lobresco dir> <case dir> ...
      the paper's 1 %-of-peak rule etc. (analyse.py) -> <case>/comparison.json
  physicalCarrier.py expect <lobresco dir> <level> <case> <expect.json>
      the checks of a finished case against case 023's tracked records

<work>/bin must hold liuCarrierTool and rhoFixedFlowFoam (the Allrun puts
them there).  LESTO_CPUS is the taskset core list of the carrier steps
(default: all cores).
"""
import filecmp
import importlib.util
import json
import os
import sys
from pathlib import Path

CASE = Path(__file__).resolve().parent
TRACKED = ['0/Y_PbI2_g', '0/Y_PbI2_s',
           'constant/speciesTransportProperties',
           'constant/thermochemistryProperties',
           'constant/thermophysicalProperties',
           'constant/turbulenceProperties',
           'system/blockMeshDict', 'system/controlDict',
           'system/fvSchemes', 'system/fvSolution']


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def redirect(module, work):
    """The module's own folder (G) becomes <work>: inputs/ and bin/ are read there."""
    inputs = work / 'inputs'
    if not inputs.exists():
        inputs.symlink_to(Path(module.G) / 'inputs')
    module.G = work


def carrier(mc, work, flow, nx, nr):
    path = work / 'carriers' / f'He{flow}'
    if not (path / 'carrier-audit.json').exists():
        mc.prepare_carrier(path, flow, nx, nr)
        audit = mc.build_carrier(path)
        print('carrier', flow, 'mL/min:', audit['status'], 'max cell flux defect / inflow',
              audit['max_cell_flux_defect_over_inflow'], flush=True)
    return path


def prepare(lobresco, work, level, names):
    mc = load(lobresco / 'make_cases.py', 'lobresco_make_cases')
    mc.CPUS = os.environ.get('LESTO_CPUS') or f'0-{os.cpu_count() - 1}'
    redirect(mc, work)
    nx, nr = mc.LEVELS[level]
    if names:
        for name in names:
            flow = mc.CASES[name]['flow']
            mc.prepare_case(work / 'cases' / name, name, carrier(mc, work, flow, nx, nr), nx)
            print('prepared', level, name, flush=True)
        return 0
    mc.prepare_case(work / 'case', 'T106_6', carrier(mc, work, 106, nx, nr), nx)
    if level != 'coarse':
        print(f'level {level}: the mesh and bins differ from the dictionaries of case 017 (expected)')
        return 0
    different = [f for f in TRACKED
                 if not filecmp.cmp(CASE / f, work / 'case' / f, shallow=False)]
    if different:
        print('FAIL the dictionaries written by make_cases.py differ from case 017:',
              ' '.join(different))
        return 1
    print('the dictionaries of case 017 = those written by make_cases.py')
    return 0


def analyse(lobresco, cases):
    an = load(lobresco / 'analyse.py', 'lobresco_analyse')
    redirect(an, next(p for p in cases[0].parents if (p / 'bin').is_dir()))
    sys.argv = ['analyse.py'] + [str(c) for c in cases]
    an.main()
    return 0


def expect(lobresco, level, name, out):
    """Checks of a finished case of case 023 against its tracked record."""
    rec = lobresco / 'cases' / level / name
    ref = json.loads((rec / 'comparison.json').read_text())
    end = (rec / 'final_profiles/time').read_text().split()[0]
    slug = lobresco.name.split('-', 1)[1]
    checks = [
        {'id': 'carrier audit', 'type': 'json', 'key': 'carrier_audit.status',
         'file': 'case.json', 'equals': 'PASS'},
        {'id': 'T_dep, solver, from the inlet [K]', 'type': 'json', 'file': 'comparison.json',
         'key': 'solver_onsets.Tdep_fromInlet.T', 'value': ref['solver_onsets']['Tdep_fromInlet']['T'],
         'abs': 0.01},
        {'id': 'wall inventory [mol]', 'type': 'json', 'file': 'comparison.json',
         'key': 'solver_onsets.inventory.wall', 'value': ref['solver_onsets']['inventory']['wall'],
         'rel': 1e-9},
        {'id': 'pair closure, every step', 'type': 'closure',
         'file': 'postProcessing/phaseChangeBalance/0/PbI2_g.dat', 'max': 1e-14},
        {'id': f'final profiles = cases/{level}/{name}/final_profiles', 'type': 'sameFiles',
         'rel': 1e-9, 'pairs': [[f'postProcessing/axialProfiles/{end}/PbI2_g_{s}.dat',
                                 f'case:{slug}/cases/{level}/{name}/final_profiles/PbI2_g_{s}.dat']
                                for s in ['gamma', 'native']]}]
    if 'dT_dep_fromInlet_K' in ref:
        checks.insert(2, {'id': 'T_dep minus the measurement [K]', 'type': 'json',
                          'file': 'comparison.json', 'key': 'dT_dep_fromInlet_K',
                          'value': ref['dT_dep_fromInlet_K'], 'abs': 0.01})
    Path(out).write_text(json.dumps({'description': f'{level} {name} against case 023',
                                     'checks': checks}, indent=1) + '\n')
    return 0


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    mode, lobresco = sys.argv[1], Path(sys.argv[2]).resolve()
    if mode == 'prepare':
        return prepare(lobresco, Path(sys.argv[3]).resolve(), sys.argv[4], sys.argv[5:])
    if mode == 'analyse':
        return analyse(lobresco, [Path(c).resolve() for c in sys.argv[3:]])
    if mode == 'expect':
        return expect(lobresco, *sys.argv[3:6])
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main())
