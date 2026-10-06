#!/usr/bin/env python3
"""Compare the original ABI on 800 cold and warm calls, without new controls."""
import math
from pathlib import Path
import struct
import sys
from bridge import Engine
old,new,system=sys.argv[1:];a=Engine(old,system);b=Engine(new,system)
assert a.elements==b.elements and a.species==b.species and a.is_gas==b.is_gas
for fn in ['num_elements','num_species','num_phases','state_size']:
    assert getattr(a.L,'gemsb_'+fn)(a.h)==getattr(b.L,'gemsb_'+fn)(b.h),fn
for i,e in enumerate(a.elements):assert a.L.gemsb_element_index(a.h,e.encode())==b.L.gemsb_element_index(b.h,e.encode())==i
for j,s in enumerate(a.species):
    assert a.L.gemsb_species_index(a.h,s.encode())==b.L.gemsb_species_index(b.h,s.encode())==j
    for fn in ['species_phase','species_is_gas','species_molar_mass']:
        assert getattr(a.L,'gemsb_'+fn)(a.h,j)==getattr(b.L,'gemsb_'+fn)(b.h,j)
    for i in range(len(a.elements)):assert a.L.gemsb_stoich(a.h,j,i)==b.L.gemsb_stoich(b.h,j,i)
a.suppress();b.suppress()
for call in range(800):
    T=500+10*(call%70);x=10**(-9+(call//70)%6)
    bulk={'He':1,'Pb':x,'I':2.000002*x*(1+1e-4*math.sin(call))}
    mode=0 if call%4==0 else 2
    rc=a.equilibrate(T,101325,bulk,mode);assert rc==b.equilibrate(T,101325,bulk,mode)
    assert a.L.gemsb_last_iterations(a.h)==b.L.gemsb_last_iterations(b.h)
    if rc in (0,1,2):
        for fn in ['species_amounts','gas_partial_pressures','log10_activities','phase_amounts','phase_log10_saturation']:
            n=a.L.gemsb_num_phases(a.h) if fn.startswith('phase_') else len(a.species)
            va=a.values(fn,n);vb=b.values(fn,n)
            assert struct.pack('='+str(n)+'d',*va)==struct.pack('='+str(n)+'d',*vb),(call,fn)
    for e in [a,b]:assert e.L.gemsb_last_seconds(e.h)>=0 and e.L.gemsb_last_error(e.h) is not None
a.close();b.close();print('800 original-ABI calls: statuses, iterations, amounts, pressures and activities bit-identical')
