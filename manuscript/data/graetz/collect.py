#!/usr/bin/env python3
"""Copy the summary outputs of the M7.1 Graetz-Robin verification into this folder.

The runs themselves live in the git-ignored work area of the solver tests,
applications/rhoFixedFlowFoam/tests/work/<version>-<options>/verification/graetz
(1.2 GB of fields).  This script keeps only small text records:

  reference.json   the reference decay rates (power series, shooting, Chebyshev)
  fits.json        one record per run: set-up, fitted and predicted decay rates
  table.md         the summary table written by graetzChecks.py
  criteria.log     the PASS/FAIL lines of `graetzRobin check` (if given)
  provenance.json  solver identity, OpenFOAM build, run dates, source commit

Usage (from the repository root):
  python3 manuscript/data/graetz/collect.py <work graetz dir> [<check log>]
"""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
KEEP = ['law', 'Bi', 'nr', 'nx', 'scheme', 'outlet', 'dt', 'Co', 'time', 'window', 'nBins', 'mu', 'm',
        'variants', 'muWall', 'muCup', 'muH', 'steady', 'muRef', 'bound']


def main(work, check_log=None):
    work = Path(work).resolve()
    fits = []
    for path in sorted(work.glob('run/**/graetzFit.json')) + sorted(work.glob('transient/**/graetzFit.json')):
        record = json.loads(path.read_text())
        entry = {k: record[k] for k in KEEP if k in record}
        entry['run'] = str(path.parent.relative_to(work))
        if entry['run'].startswith('transient'):
            entry['scheme'] = 'linear (time-accurate)'
        fits.append(entry)
    (HERE / 'fits.json').write_text(json.dumps(fits, indent=1) + '\n')
    for name in ['reference.json', 'table.md']:
        (HERE / name).write_bytes((work / name).read_bytes())
    if check_log:
        lines = [l for l in Path(check_log).read_text().splitlines() if re.match(r'^(PASS|FAIL|NOTRUN)\s', l)]
        lines = list(dict.fromkeys(lines))   # the check prints its last lines again at the end
        (HERE / 'criteria.log').write_text('\n'.join(lines) + '\n')
    # identity of the solver: the first field of every run key is the md5 of the solver executable,
    # the carrier helper, the case template and the two Python scripts (graetzRobin, inputsKey)
    keys = sorted({p.read_text().split()[0] for p in work.glob('run/**/key')} |
                  {p.read_text().split()[0] for p in work.glob('transient/**/key')})
    log = next(work.glob('run/zeroGradient/bi1/halfCell/nr064/log')).read_text()
    header = {k: re.search(rf'^{k}\s*: (.*)$', log, re.MULTILINE).group(1).strip() for k in ['Build', 'Exec', 'Date']}
    exe = Path(header['Exec'])
    tests = REPO / 'applications/rhoFixedFlowFoam/tests/verification'
    inputs = hashlib.md5()
    inputs.update((hashlib.md5(exe.read_bytes()).hexdigest() + '  -\n').encode())
    inputs.update((hashlib.md5((work / 'bin/graetzCarrier').read_bytes()).hexdigest() + '  -\n').encode())
    inputs.update((tests / 'graetzCase.py').read_bytes() + (tests / 'graetz.py').read_bytes())
    for f in sorted((p for p in (tests / 'graetzWedge').rglob('*') if p.is_file()), key=str):
        inputs.update(f.read_bytes())
    head = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    provenance = {
        'study': 'M7.1 Graetz-Robin verification (applications/rhoFixedFlowFoam/doc/phase-change-plan.md, M7 acceptance 1)',
        'openfoam_build': header['Build'], 'run_date': header['Date'],
        'solver_executable': header['Exec'],
        'solver_sha256_at_collection': hashlib.sha256(exe.read_bytes()).hexdigest(),
        'run_input_keys': keys, 'input_key_recomputed_at_collection': inputs.hexdigest(),
        'runs_made_with_this_solver': keys == [inputs.hexdigest()],
        'repository_head_at_collection': head,
        'number_of_runs': len(fits)}
    (HERE / 'provenance.json').write_text(json.dumps(provenance, indent=1) + '\n')
    print(json.dumps(provenance, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
