#!/usr/bin/env python3
"""Read-only prerequisites for the specified scientific follow-up study."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil

EXPECTED = '16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True, type=Path)
    parser.add_argument('--prepared', type=Path)
    parser.add_argument('--solver', type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve()
    prepared = (args.prepared or repo/'run/024-liu2025/work/prepared').resolve()
    solver = (args.solver or prepared/'bin/rhoFixedFlowFoam').resolve()
    checks = []
    def record(name, ok, detail):
        checks.append(dict(check=name,passed=bool(ok),detail=detail))
    record('OpenFOAM environment', bool(os.environ.get('WM_PROJECT_DIR')),
           os.environ.get('WM_PROJECT_VERSION','not loaded'))
    for name in ['blockMesh','checkMesh','rhoSimpleFoam']:
        path = shutil.which(name)
        record(name, bool(path), path)
    sha = hashlib.sha256(solver.read_bytes()).hexdigest() if solver.is_file() else None
    record('Original solver hash',sha==EXPECTED,sha)
    record('Solver executable',solver.is_file() and os.access(solver,os.X_OK),str(solver))
    for base in BASES:
        case = prepared/'cases'/base
        required = [case/'0',case/'constant',case/'system',case/'benchmark.json',
                    case/'postProcessing/axialProfiles/60/element_I_native.dat']
        missing = [str(p) for p in required if not p.exists()]
        record(base+' baseline inputs',not missing,missing)
        manifest = case/'benchmark.json'
        if manifest.is_file():
            data=json.loads(manifest.read_text())
            folder=Path(data['table_directory'])
            bad=[]
            for name,expected in data['table_sha256'].items():
                path=folder/name
                if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
                    bad.append(name)
            record(base+' table fingerprints',not bad,bad)
    # This preflight deliberately neither sources shell files nor executes solvers.
    report={'status':'PREREQUISITES_FOUND' if all(r['passed'] for r in checks) else 'MISSING_PREREQUISITES',
            'scope':'Read-only initial check. Carrier tooling and external thermodynamic inputs still need inspection before preparing finer meshes.',
            'repo':str(repo),'prepared':str(prepared),'checks':checks}
    print(json.dumps(report,indent=2))
    return 0 if report['status']=='PREREQUISITES_FOUND' else 2


if __name__=='__main__':
    raise SystemExit(main())
