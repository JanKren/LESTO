#!/usr/bin/env python3
"""Copy the final wall profiles and the end-of-run ledger of the eight production-scheme
sensitivity runs into this folder.

    python3 collect.py <directory holding the cases <base>__<variant>>

Writes profiles/<case>/{element_I,PbI2_g,BiI3_g}_native.dat (60 s) and ledger.json.
"""
import json, re, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = ['element_I_native.dat', 'PbI2_g_native.dat', 'BiI3_g_native.dat']


def ledger(case):
    log = (case/'log.solver').read_text()
    tail = log[-200000:]
    bal = {m.group(1): dict(closure_relative=float(m.group(3)), solver_defect_relative=float(m.group(2)))
           for m in re.finditer(r'^Balance (\S+) \[mol\]: .*?solverDefect \S+ \(relative (\S+)\) closure \S+ \(relative (\S+)\)',
                                tail, re.M)}
    outer = {m.group(1): int(m.group(2)) for m in
             re.finditer(r'^Phase change summary Y_(\S+): (\d+) steps without outer convergence', tail, re.M)}
    resid = [float(v) for v in re.findall(r'max element residue (\S+)', tail)]
    times = re.findall(r'^Time = (\S+)', tail, re.M)
    return dict(end_time=float(times[-1]), ended=log.rstrip().endswith('End'),
                execution_h=float(re.findall(r'ExecutionTime = (\S+)', tail)[-1])/3600,
                final_balance=bal, steps_without_outer_convergence=outer,
                final_max_element_residue=max(resid) if resid else None)


if __name__ == '__main__':
    src = Path(sys.argv[1])
    out = {}
    for case in sorted(p for p in src.iterdir() if '__' in p.name):
        dst = HERE/'profiles'/case.name
        dst.mkdir(parents=True, exist_ok=True)
        for f in FILES:
            shutil.copy(case/'postProcessing/axialProfiles/60'/f, dst/f)
        out[case.name] = ledger(case)
        print(case.name, out[case.name]['end_time'], f"{out[case.name]['execution_h']:.1f} h")
    (HERE/'ledger.json').write_text(json.dumps(out, indent=1) + '\n')
