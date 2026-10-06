#!/usr/bin/env python3
"""M10b prescribed-source channel and Lie-split time-step study.
Requires external Bi tables. Results are code qualification, not predictions
of the unpublished experimental release history or iodine activity in LBE.
"""
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil
import sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('speciation_cases',HERE/'cases.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
base=module.base
sys.path.insert(0,str(HERE.parent/'multiSpecies'))
import scientific as independent
from make_kf_tables import formation,pbbii_system,write

def setup(case,data,nasa,tf,dt=None):
    case,amounts=independent.setup(case,data)
    carrier=case/'system/carrierDict';carrier.write_text(carrier.read_text().replace('Tin       1000;','Tin       1173;'))
    tc=case/'constant/thermochemistryProperties';tc.write_text(tc.read_text().replace('nBins 100;}', 'nBins 100;} cm {min 0; max .1; nBins 10;}'))
    if dt is not None:
        aux=module.prepare(case.parent/(case.name+'_inputs'))
        for n in module.NAMES:
            if n not in ['PbI2_g','BiI3_g']:shutil.copy2(aux/'0'/('Y_'+n),case/'0'/('Y_'+n))
        for name in ['speciesTransportProperties']:shutil.copy2(aux/'constant'/name,case/'constant'/name)
        for name in ['fvSolution','fvSchemes']:shutil.copy2(aux/'system'/name,case/'system'/name)
        tc.write_text(tc.read_text()+'\n'+'respeciation '+(aux/'constant/thermochemistryProperties').read_text().split('respeciation ')[1])
        write(case/'constant/log10Kf.csv',module.NAMES,formation(pbbii_system(nasa,tf)))
        ctrl=case/'system/controlDict';s=ctrl.read_text();s=re.sub(r'deltaT\s+[0-9.]+;',f'deltaT {dt};',s);s=re.sub(r'writeInterval\s+\d+;',f'writeInterval {round(3/dt)};',s);ctrl.write_text(s)
    return case

def deposits(case,name):
    text=(case/f'postProcessing/axialProfiles/3/{name}_cm.dat').read_text()
    keys=re.search(r'# columns:? (.*)',text)
    if not keys:
        # Profiles document the numerical columns in the final comment line.
        header=[l[1:].strip().split() for l in text.splitlines() if l.startswith('#')][-1]
    else:header=keys[1].split()
    rows=[list(map(float,l.split())) for l in text.splitlines() if l and not l.startswith('#')]
    return rows,header

def l1(a,b):return sum(abs(x-y) for x,y in zip(a,b))/max(sum(abs(y) for y in b),1e-300)

if __name__=='__main__':
    work,solver,carrier,data,nasa,tf=sys.argv[1:];work=Path(work);work.mkdir(parents=True,exist_ok=True)
    baseline=setup(work/'baseline',Path(data),nasa,tf);base.run(baseline,solver,carrier)
    import subprocess
    program=work/'testSpeciationFields'
    assert subprocess.run(['g++','-O2','-std=c++17',str(HERE/'testSpeciation.C'),'-o',str(program)]).returncode==0
    report={};curves={}
    for dt in [.004,.002,.001]:
        case=setup(work/('dt'+str(dt)),Path(data),nasa,tf,dt);base.run(case,solver,carrier);module.check(case);module.equilibrium(case,'3',program)
        for n in ['PbI2_g','BiI3_g']:
            a,keys=deposits(case,n);b,_=deposits(baseline,n)
            # Axial profile column densityWall is the deposited mol/m.
            col=keys.index('wall') if 'wall' in keys else keys.index('wall_mol_per_m')
            curve=[r[col] for r in a];reference=[r[col] for r in b]
            curves[(dt,n)]=curve;report[f'{dt}_{n}']={'L1_vs_M10a':l1(curve,reference)}
            assert l1(curve,reference)<=.01,(dt,n,report)
        print('completed dt',dt,flush=True)
    for n in ['PbI2_g','BiI3_g']:
        error=l1(curves[(.002,n)],curves[(.001,n)])
        report[n+'_dt2_vs_dt1']=error;assert error<=.01,(n,error)
    (work/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
