#!/usr/bin/env python3
"""Preliminary independent-iodide channel; prescribed split is not speciation."""
from pathlib import Path
import json
import math
import re
import shutil
import sys
from cases import prepare, run, balance, elements, SPECIES, REPO


def setup(case, data):
    data=Path(data).resolve()
    case=prepare(Path(case),['PbI2','BiI3'])
    shutil.copy2(REPO/'thermochemistry/systems/data/pv_PbI2_phases.csv',case/'constant/pv_PbI2.csv')
    shutil.copy2(data/'pv_BiI3_phases.csv',case/'constant/pv_BiI3.csv')
    # Liu LBE-I_SS_II observed iodide ratio. Scale to a prescribed iodine
    # mole fraction 3.9e-6 in the test carrier; this is a code qualification,
    # not the measured release history or a prediction of its species split.
    nP=2.407762417796063e-7;nB=3.069389519282897e-7
    totalI=2*nP+3*nB
    heFlow=.01743*(.0048*.001)/.0040026
    duration=2.0
    amounts={'PbI2':heFlow*3.9e-6*nP/totalI*duration,
             'BiI3':heFlow*3.9e-6*nB/totalI*duration}
    p=case/'constant/thermochemistryProperties';s=p.read_text()
    s=s.replace('observable wallPlusGas;', 'observable wall;')
    s=s.replace('accommodation .8; Ce .0001; kineticScale 1;','accommodation 1; Ce 1; kineticScale 1;').replace('accommodation .6; Ce .0002; kineticScale .5;','accommodation 1; Ce 1; kineticScale 1;')
    for n,amount in amounts.items():
        s=re.sub(rf'{n}_g [0-9.e+-]+;',f'{n}_g {amount:.17g};',s)
    s=s.replace('startTime .0003; duration .0405;','startTime 0; duration 2;')
    p.write_text(s)
    p=case/'system/controlDict';s=p.read_text().replace('endTime           0.05;','endTime           3;').replace('writeInterval     5;','writeInterval     500;');p.write_text(s)
    (case/'benchmark.json').write_text(json.dumps({'model':'independent PbI2/BiI3 pairs, prescribed source', 'iodine_mole_fraction':3.9e-6,'amounts_mol':amounts,'carrier':'synthetic plane channel, He, uniform axial mass flux','diffusivity':'constant 1e-4 m2/s (provisional)','external_Bi_data':str(data),'observable':'1% onset for code qualification; Tpeak separately reported'},indent=2)+'\n')
    return case,amounts


def check(case,amounts):
    report={}
    for n,amount in amounts.items():
        row=balance(case,n)[-1]
        assert abs(row['released']-amount*SPECIES[n][0]) <= 1e-14*amount*SPECIES[n][0]
        assert abs(row['closure']) <= 1e-14*amount*SPECIES[n][0]
        f=case/f'postProcessing/axialProfiles/3/{n}_g_native.dat'
        text=f.read_text()
        onset=re.search(r'# Tdep.*?status found T ([0-9.e+-]+) x',text)
        assert onset,text[:2000]
        # Independent ideal plug-flow threshold: prescribed partial pressure
        # equals the stable pure-condensate vapour pressure. Bisection with
        # inverse-temperature interpolation on the external table.
        table=case/f'constant/pv_{n}.csv'
        rows=[list(map(float,l.split(','))) for l in table.read_text().splitlines() if l and not l.startswith('#')]
        heFlow=.01743*.0048*.001/.0040026
        partial=101325*amount/(2*heFlow)
        def logp(T):
            for a,b in zip(rows,rows[1:]):
                if a[0] <= T <= b[0]:
                    w=(1/T-1/a[0])/(1/b[0]-1/a[0])
                    return min((1-w)*a[j]+w*b[j] for j in range(1,len(a)))
            raise ValueError(T)
        lo,hi=350.,1000.
        for _ in range(80):
            mid=(lo+hi)/2
            if logp(mid)>math.log10(partial): hi=mid
            else: lo=mid
        crossing=(lo+hi)/2
        report[n]={'onset_K':float(onset[1]),'ideal_plug_flow_saturation_K':crossing,'difference_K':float(onset[1])-crossing,'reference_K':566 if n=='PbI2' else 450}
        assert abs(float(onset[1])-crossing)<=10,(n,report[n])
        assert abs(crossing-report[n]['reference_K'])<=10,(n,report[n])
    (case/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    case,solver,carrier,data=sys.argv[1:]
    case,amounts=setup(Path(case),data)
    run(case,solver,carrier)
    check(case,amounts)
