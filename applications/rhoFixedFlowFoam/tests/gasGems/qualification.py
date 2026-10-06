#!/usr/bin/env python3
"""Reduced map: 3 mixes x 24 temperatures x 13 dilutions, exact GEMS oracle."""
import importlib.util
import json
import math
from pathlib import Path
import sys
from bridge import Engine
from precise import refine
root,library,nasa,tf,out=sys.argv[1:];root=Path(root);out=Path(out)
spec=importlib.util.spec_from_file_location('external_oracle',root/'gems/py/gemsbridge_py.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
oracle=ref.Engine(str(root/'thermo/PbBiIHe/PbBiIHe-dat.lst'),log_dir='/tmp',log_level=6);oracle.suppress()
e=Engine(library,root/'thermo/PbBiIHe/PbBiIHe-dat.lst');e.suppress();e.L.gemsb_set_warm_iteration_limit(e.h,200)
sys.path.insert(0,str(Path(__file__).resolve().parents[4]/'thermochemistry/systems'))
from make_kf_tables import SPECIES,formation,pbbii_system
fun=formation(pbbii_system(nasa,tf));names=list(SPECIES);gas=[n.replace('_g','(g)') for n in names];idx=[e.species.index(n) for n in gas]
mixes={'PbI2':{'Pb':1,'I':2.000002},'BiI3':{'Bi':1,'I':3.000003},'LBEI':{'Pb':.44,'Bi':.56,'I':2.56000256}}
report={};gasErrors={};maxK=0;maxKKT=0;maxOracleDifference=0
with out.open('w') as target:
    target.write('# T,scale,bPb,bBi,bI,9 lnK/reference concentration pairs,3 warm potentials\n')
    for name,mix in mixes.items():
        for ti in range(24):
            T=350+(1173-350)*ti/23;g=e.g0(T,101325)
            ks=[]
            for n,j in zip(names,idx):
                k=-g[j]+sum(v*g[e.species.index(a+'(g)')] for a,v in SPECIES[n].items());ks.append(k)
                maxK=max(maxK,abs(k/math.log(10)-fun(n,T)))
            for xi in range(13):
                x=10**(-12+xi*11/12);nu=x/(1-x);bulk={'He':1,**{a:nu*v for a,v in mix.items()}}
                assert oracle.ideal(T,101325,bulk) in (0,1)
                assert oracle.polish(2e-15) in (0,1)
                maxKKT=max(maxKKT,oracle.L.gemsb_last_kkt_residual(oracle.h))
                amounts=oracle.amounts();total=sum(amounts[s] for s,g in zip(oracle.species,oracle.is_gas) if g)
                volume=total*8.314462618*T/101325;b={a:bulk.get(a,0)/volume for a in ['Pb','Bi','I']}
                trace=sum(amounts[s]/volume for s in gas);carrier=101325/(8.314462618*T);fraction=trace/(trace+carrier)
                key=str(math.floor(math.log10(x)));entry=report.setdefault(name+':'+key,{'points':0,'calls':0,'guards':0,'failures':0,'max_primal_dlog10':0})
                entry['points']+=1
                pi=[0.,0.,0.]
                if fraction<1e-12 or fraction>3.2e-4:entry['guards']+=1
                else:
                    entry['calls']+=1;rc=e.equilibrate(T,101325,{'He':carrier,**b})
                    if rc not in (0,1):entry['failures']+=1;e.valid.value=0
                    else:
                        u=e.values('element_potentials',len(e.elements));pi=[u[e.elements.index(a)]-g[e.species.index(a+'(g)')] if b[a]>0 else 0 for a in ['Pb','Bi','I']]
                        values=e.values('species_amounts');pressures=e.values('gas_partial_pressures')
                        # Raw IPM differences are diagnostics; the polished values are gated by the C++ kernel against this oracle.
                        for s,j in zip(gas,idx):
                            p=amounts[s]/total*101325
                            if p>0 and pressures[j]>0:
                                difference=abs(math.log10(pressures[j]/p));bucket=name+':'+s+':'+str(math.floor(math.log10(p/101325)))
                                gasErrors[bucket]=max(gasErrors.get(bucket,0),difference)
                                if p/101325>1e-10:entry['max_primal_dlog10']=max(entry['max_primal_dlog10'],difference)
                scale=1e5/(8.314462618*T)
                seed=[math.log(amounts[a+'(g)']/volume/scale) if b[a]>0 else 0 for a in ['Pb','Bi','I']]
                exact=refine(ks,scale,[b[a] for a in ['Pb','Bi','I']],seed)
                for s,c in zip(gas,exact):
                    if c>1e-10*max(b.values()):maxOracleDifference=max(maxOracleDifference,abs(amounts[s]/volume/c-1))
                row=[T,scale,b['Pb'],b['Bi'],b['I']]
                for k,c in zip(ks,exact):row.extend([k,c])
                row.extend(pi);target.write(','.join(format(v,'.17g') for v in row)+'\n')
report['max_raw_error_per_species_mole_fraction_decade']=gasErrors;report['maximum_double_oracle_vs_Decimal_species_error']=maxOracleDifference;report['maximum_Kf_difference_decades']=maxK;report['maximum_oracle_KKT_residue']=maxKKT
assert maxK<=1e-3,maxK
out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print('936 reduced-map points; max |GEMS - open log10Kf|',maxK,'max oracle KKT',maxKKT,'double oracle vs Decimal',maxOracleDifference,flush=True)
e.close()
