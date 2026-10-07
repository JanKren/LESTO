#!/usr/bin/env python3
"""Replay the first 10 ms and quantify iodine in transient invalid-bulk cells.

All physical inputs and the private solver are retained. Writing each step
exposes skipped cells without adding a concentration floor or changing gates.
"""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from carriers import REPO
from analyse_matrix import FORMULAS

CASES=['LBE-I_SS_II_Q2_shift29k_sat0.1',
       'LBE-I_SiO2_I_Q2_baseline_sat0.001','LBE-I_SiO2_I_Q2_baseline_sat0.01']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();assert not args.out.exists(),'Use a fresh output directory'
    args.out.mkdir(parents=True)
    spec=importlib.util.spec_from_file_location('checks',REPO/'applications/rhoFixedFlowFoam/tests/phaseChange/checks.py')
    checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(checks)
    def field(path):
        a=np.array(checks.readScalarField(path)[0],dtype=np.longdouble)
        return np.full(18432,a[0],dtype=np.longdouble) if len(a)==1 else a
    for name in CASES:
        source=args.root/name;case=args.out/name;case.mkdir()
        for folder in ['0','constant','system']:shutil.copytree(source/folder,case/folder)
        shutil.copyfile(source/'benchmark.json',case/'benchmark.json')
        control=case/'system/controlDict';text=control.read_text()
        assert 'endTime 60;' in text and 'writeInterval 500;' in text
        control.write_text(text.replace('endTime 60;','endTime .01;').replace('writeInterval 500;','writeInterval 1;').replace('purgeWrite 2;','purgeWrite 0;'))
    def run(name):
        case=args.out/name
        with (case/'log.solver').open('w') as log:
            subprocess.run([str(args.solver.resolve()),'-case',str(case.resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
        print('Finished startup replay',name,flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,CASES))
    output=[]
    for name in CASES:
        case=args.out/name;m=json.loads((case/'benchmark.json').read_text())
        props=(case/'constant/speciesTransportProperties').read_text()
        masses={n:float(re.search(r'\b'+n+r'\s*\{[^}]*?molarMass\s+([^;]+);',props)[1]) for n in FORMULAS}
        rho=field(case/'0/rho');T=field(case/'0/T')
        source_i=sum(FORMULAS[n][2]*rate for n,rate in m['source_wedge_mol_s'].items())
        # Entire uniform axial slice volume is an upper bound on every
        # radial cell volume, avoiding any dependence on cell numbering.
        slice_volume=m['carrier']['radius_m']**2*math.tan(math.radians(2.5))*1.25/384
        events=[]
        for step in range(1,6):
            t=step*.002;bulk=np.zeros((3,18432),dtype=np.longdouble);positive_i=np.zeros(18432,dtype=np.longdouble)
            for gas,nu in FORMULAS.items():
                amount=rho*field(case/f'{t:g}'/('Y_'+gas))/masses[gas]
                for e in range(3):bulk[e]+=nu[e]*amount
                positive_i+=nu[2]*np.maximum(amount,0)
            invalid=(T>=350)&np.any(bulk<0,axis=0)
            if invalid.any():
                events.append({'time_s':t,'invalid_cells':int(invalid.sum()),
                    'positive_I_over_supplied_I_upper_bound':float(positive_i[invalid].sum()*slice_volume/(source_i*t))})
        logged=[int(x) for x in re.findall(r'invalid-element cells (\d+)',(case/'log.solver').read_text())]
        assert max(logged)==max((r['invalid_cells'] for r in events),default=0),name
        output.append({'case':name,'events':events,'maximum_logged_invalid_cells':max(logged)})
    result={'status':'COMPLETED','solver_sha256':hashlib.sha256(args.solver.read_bytes()).hexdigest(),
        'dt_s':.002,'end_time_s':.01,'results':output,
        'scope':'Five-step startup replay; original physical inputs and source rates; positive iodine bound uses entire axial-slice volume per affected radial cell.'}
    (args.out/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
