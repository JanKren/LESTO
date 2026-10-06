#!/usr/bin/env python3
"""M10c prescribed-source channel: kernel/frozen/local deposit comparisons."""
import importlib.util
import json
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('speciation_scientific',HERE.parent/'speciation/scientific.py')
science=importlib.util.module_from_spec(spec);spec.loader.exec_module(science)
base=science.base

def setup(case,data,nasa,tf,system,mode):
    science.setup(case,Path(data),nasa,tf,.001)
    if mode!='kernel':
        tc=case/'constant/thermochemistryProperties';s=tc.read_text().replace('engine kernel;', 'engine GEMS;')
        s=s.replace('minTemperature 350;',f'GEMS {{mode {mode}; system "{system}"; carrier "He(g)";}} minTemperature 350;');tc.write_text(s)
    return case

def main():
    checking=len(sys.argv)==3 and sys.argv[1]=='--check'
    if checking:work=Path(sys.argv[2])
    else:
        work,solver,carrier,data,nasa,tf,system=sys.argv[1:];work=Path(work);work.mkdir(parents=True,exist_ok=True)
    curves={};report={};calls=0;seconds=0
    for mode in ['kernel','frozen','local']:
        case=work/mode
        if not checking:setup(case,data,nasa,tf,system,mode);base.run(case,solver,carrier)
        science.module.check(case)
        for n in ['PbI2_g','BiI3_g']:
            rows,keys=science.deposits(case,n);col=keys.index('wall') if 'wall' in keys else keys.index('wall_mol_per_m')
            curves[mode,n]=[r[col] for r in rows]
        if mode=='local':
            log=(case/'log.solver').read_text()
            for m in re.finditer(r'Gas GEMS local: calls (\d+); iterations (\d+); failures (\d+); guards (\d+); mean call us ([0-9.e+-]+)',log):
                c=int(m[1]);calls+=c;seconds+=c*float(m[5])*1e-6
            assert calls>0,'Vacuous local engine: no IPM calls'
            report['local_calls']=calls;report['mean_IPM_plus_polish_us']=1e6*seconds/calls
            # Timing is reported independently from deterministic accuracy gates.
        print('completed mode',mode,flush=True)
    for n in ['PbI2_g','BiI3_g']:
        a=science.l1(curves['frozen',n],curves['kernel',n]);b=science.l1(curves['local',n],curves['frozen',n])
        assert a<=1e-3,(n,a);assert b<=1e-9,(n,b)
        report[n]={'frozen_vs_kernel_L1':a,'local_vs_frozen_L1':b}
    report['rechecked_existing_runs']=checking
    (work/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
