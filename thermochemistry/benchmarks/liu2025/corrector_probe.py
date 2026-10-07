#!/usr/bin/env python3
"""Short matched restarts to diagnose the full matrix's outer residual floor.

This is an iteration-sensitivity diagnostic at a saved late-time checkpoint,
not a replacement for a complete mesh/time-step convergence study.
"""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from carriers import REPO
from report import rows

CASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.01',
         'LBE-I_SS_II_Q2_shift29k_sat0.01', 'LBE-I_SiO2_I_Q1_baseline']
GASES = ['PbI2_g','BiI3_g','Pb_g','PbI_g','Bi_g','Bi2_g','BiI_g','I_g','I2_g']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--solver',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--reuse-completed',action='store_true',help='Analyse existing completed probes without rerunning them')
    args=parser.parse_args()
    if args.out.exists():
        assert args.reuse_completed, 'Use a fresh output directory'
    else:
        args.out.mkdir(parents=True)
    solver=args.solver.resolve()
    spec=importlib.util.spec_from_file_location('field_checks',REPO/'applications/rhoFixedFlowFoam/tests/phaseChange/checks.py')
    checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(checks)
    configurations=[('control',3,'1e-16'),('more_correctors',8,'1e-16'),('tolerance_1e14',3,'1e-14')]
    prepared=[]
    for name in CASES:
        source=args.root/name
        assert (source/'59').exists(),source
        for label,ncorr,tolerance in configurations:
            case=args.out/name/label
            if case.exists():
                assert args.reuse_completed and (case/'59.04').exists(),case
                assert (case/'log.solver').read_text().rstrip().endswith('End'),case
                continue
            case.mkdir(parents=True)
            for folder in ['0','59','constant','system']:
                shutil.copytree(source/folder,case/folder)
            control=case/'system/controlDict';text=control.read_text()
            assert 'endTime 60;' in text and 'nCorr 3;' in text
            text=text.replace('endTime 60;','endTime 59.04;').replace('writeInterval 500;','writeInterval 20;')
            text=text.replace('nCorr 3;',f'nCorr {ncorr};').replace('tolerance 1e-16;',f'tolerance {tolerance};')
            control.write_text(text)
            prepared.append(case)
    def run(case):
        with (case/'log.solver').open('w') as log:
            subprocess.run([str(solver),'-case',str(case.resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
        assert (case/'59.04').exists(),case
        print('Finished',case.parent.name,case.name,flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,prepared))
    comparisons=[]
    for name in CASES:
        base=args.out/name/'control'
        for label,_,_ in configurations[1:]:
            candidate=args.out/name/label
            fields={}
            for gas in GASES:
                a=np.array(checks.readScalarField(base/'59.04'/('Y_'+gas))[0])
                b=np.array(checks.readScalarField(candidate/'59.04'/('Y_'+gas))[0])
                if len(a)==1:a=np.full(18432,a[0])
                if len(b)==1:b=np.full(18432,b[0])
                assert np.isfinite(a).all() and np.isfinite(b).all()
                scale=float(np.max(np.abs(a)))
                fields[gas]={'max_difference_over_max_control':float(np.max(np.abs(a-b))/scale) if scale else None,
                             'unweighted_cell_L1_relative_to_control':float(np.abs(a-b).sum()/np.abs(a).sum()) if np.abs(a).sum() else None}
            profiles={};increments={}
            for species in ['PbI2_g','BiI3_g','element_I']:
                a=rows(base/'postProcessing/axialProfiles/59.04'/(species+'_native.dat'))
                b=rows(candidate/'postProcessing/axialProfiles/59.04'/(species+'_native.dat'))
                av=np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in a]);bv=np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in b])
                profiles[species]=float(np.abs(av-bv).sum()/av.sum()) if av.sum() else None
                initial=rows(args.root/name/'postProcessing/axialProfiles/59'/(species+'_native.dat'))
                iv=np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in initial])
                change=float(np.abs(av-iv).sum())
                increments[species]=float(np.abs(av-bv).sum()/change) if change else None
            for element in ['Pb','Bi','I']:
                ledger=rows(candidate/'postProcessing/phaseChangeElements/59'/(element+'.dat'))
                assert abs(ledger[-1]['time']-59.04)<1e-7
                assert max(abs(r['closure'])/max(r['reference'],1e-300) for r in ledger)<1e-12
            comparisons.append({'case':name,'variant':label,'gas_fields':fields,'cumulative_wall_profile_L1':profiles,
                                'wall_increment_profile_L1':increments})
    result={'status':'PASS','solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
            'start_time_s':59,'end_time_s':59.04,'steps':20,'dt_s':.002,'cases':CASES,
            'configurations':configurations,'comparisons':comparisons,
            'scope':'Short matched late-time restarts; source, transport discretisation and linear solver tolerances unchanged. No full-run iteration-independence claim.'}
    (args.out/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('PASS: matched corrector/tolerance restart diagnostic',flush=True)


if __name__=='__main__':
    main()
