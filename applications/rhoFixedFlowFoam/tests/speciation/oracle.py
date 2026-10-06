#!/usr/bin/env python3
"""Export the external exact GEMS oracle on every stored gas-map input.
Nothing from the external dataset is copied into the repository. The bridge
must have the design's ideal-equilibrium/G0 getters (M10c production support
is a separate milestone). The solver's default build never links it.
"""
import csv
import importlib.util
import math
from pathlib import Path
import sys
root,out=sys.argv[1:];root=Path(root)
spec=importlib.util.spec_from_file_location('m10_bridge',root/'gems/py/gemsbridge_py.py')
gb=importlib.util.module_from_spec(spec);spec.loader.exec_module(gb)
eng=gb.Engine(str(root/'thermo/PbBiIHe/PbBiIHe-dat.lst'),log_dir='/tmp',log_level=0)
eng.suppress()
sp=['Pb(g)','Pb2(g)','PbI(g)','PbI2(g)','PbI3(g)','PbI4(g)','I(g)','I2(g)','Bi(g)','Bi2(g)','BiI(g)','BiI3(g)']
# The extended system intentionally has no Pb2/PbI3/PbI4 data: suppress any
# absent species in the exported kernel by assigning negligible Kf.
formulas=[{'Pb':1},{'Pb':2},{'Pb':1,'I':1},{'Pb':1,'I':2},{'Pb':1,'I':3},{'Pb':1,'I':4},{'I':1},{'I':2},{'Bi':1},{'Bi':2},{'Bi':1,'I':1},{'Bi':1,'I':3}]
ex=1+1e-6
mixes={'PbI2':{'Pb':1,'I':2*ex},'BiI3':{'Bi':1,'I':3*ex},'LBEI':{'Pb':.44,'Bi':.56,'I':(2*.44+3*.56)*ex},'LBEV':{'Pb':.2224,'Bi':.7727,'I':.0049},'LBE9':{'Pb':7.43e-4/1.665e-3,'Bi':9.19e-4/1.665e-3,'I':3e-6/1.665e-3},'BiI':{'Bi':1,'I':1},'IRICH':{'Pb':1/12,'Bi':1/12,'I':10/12},'METAL':{'Pb':.45,'Bi':.55}}
count=0;worst=0;cache={}
with open(out,'w') as target:
    target.write('# T,scale,bPb,bBi,bI,'+','.join('lnKf_'+n+',c_'+n for n in sp)+'\n')
    for path in sorted((root/'gems/maps').glob('map_*_gas_*.csv')):
        # PbI2He and PbBiIHe use the same open Pb-I records; the full-system
        # oracle removes exactly absent Bi. Use the unfloored physical bulk.
        for r in csv.DictReader(path.open()):
            T=float(r['T']);x=float(r['x']);step=int(r['step'])
            # Recover the design's quarter-decade grid from the rounded CSV.
            exponent=round(4*math.log10(x))/4
            assert abs(math.log10(x)-exponent)<1e-5
            x=10**exponent;nu=x/(1-x)
            original_elements=['He','I','Pb'] if path.name.startswith('map_PbI2He_') else eng.elements
            b={'He':1.0}
            b.update({e:nu*v*(1+1e-3*math.sin(1.3*original_elements.index(e)+.7*step)) for e,v in mixes[r['mix']].items()})
            assert eng.ideal(T,101325,b) in (0,1),(path,r,eng.error())
            assert eng.polish(2e-15) in (0,1),(path,r,eng.error())
            kkt=eng.L.gemsb_last_kkt_residual(eng.h);worst=max(worst,kkt)
            assert kkt<=1e-11,(path,r,kkt)
            n=eng.amounts();p=eng.partial_pressures();total=sum(n[s] for s,g in zip(eng.species,eng.is_gas) if g)
            scale=1e5/(8.314462618*T)
            # concentration scale is set by the exact gas volume, including
            # association/dissociation. This equals the frozen-carrier model
            # when the element concentrations, rather than carrier moles, are inputs.
            volume=total*8.314462618*T/101325
            if T not in cache:
                g0=eng.g0(T,101325)
                cache[T]=[(-g0[s]+sum(v*g0[e+'(g)'] for e,v in f.items())) if s in g0 else -700.0 for s,f in zip(sp,formulas)]
            row=[T,scale,b.get('Pb',0)/volume,b.get('Bi',0)/volume,b.get('I',0)/volume]
            for s,k in zip(sp,cache[T]):row.extend([k,n.get(s,0)/volume])
            target.write(','.join(format(v,'.17g') for v in row)+'\n');count+=1
        print(path.name,count,flush=True)
assert count==160740,count
print('oracle points',count,'maximum KKT residual',worst)
