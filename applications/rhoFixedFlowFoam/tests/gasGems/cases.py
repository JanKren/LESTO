#!/usr/bin/env python3
"""Frozen/local gas GEMS conservation, MPI agreement, restart and input gates."""
import importlib.util
from pathlib import Path
import re
import shutil
import sys
from bridge import Engine
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('gas_spec_cases',HERE.parent/'speciation/cases.py')
sc=importlib.util.module_from_spec(spec);spec.loader.exec_module(sc)
base=sc.base

def prepare(case,system,library,mode):
    e=Engine(library,system)
    paired=['PbI2','BiI3'] if 'Bi' in e.elements else ['PbI2']
    base.prepare(case,paired)
    names=[];props=case/'constant/speciesTransportProperties';s=props.read_text()
    template=(case/'0/Y_PbI2_g').read_text();extra=[]
    for j,gas in enumerate(e.species):
        if not e.is_gas[j] or gas=='He(g)':continue
        name=gas.replace('(g)','_g');names.append(name)
        if name[:-2] in paired:continue
        M=e.L.gemsb_species_molar_mass(e.h,j)
        atoms={a:e.L.gemsb_stoich(e.h,j,i) for i,a in enumerate(e.elements) if e.L.gemsb_stoich(e.h,j,i)}
        extra.append(name)
        s+=f'\n{name} {{state gas; molarMass {M:.17g}; diffusivityModel constant; D 1e-4; formula {{'+''.join(f'{a} {n:g};' for a,n in atoms.items())+'}}\n'
        (case/'0'/('Y_'+name)).write_text(template.replace('PbI2_g',name))
    s=re.sub(r'species \(([^;]+)\);',lambda m:'species ('+m[1]+' '+' '.join(extra)+');',s)
    props.write_text(s)
    sol=case/'system/fvSolution';s=sol.read_text();position=s.rfind('}')
    s=s[:position]+''.join(f'Y_{n} {{solver PBiCGStab; preconditioner DILU; tolerance 1e-18; relTol 0;}}\n' for n in extra)+s[position:]
    sol.write_text(s.replace('tolerance 1e-15;','tolerance 1e-18;'))
    schemes=case/'system/fvSchemes';schemes.write_text(schemes.read_text().replace('default none;','default none; '+''.join(f'div(phi,Y_{n}) Gauss upwind; ' for n in extra)))
    tc=case/'constant/thermochemistryProperties';s=tc.read_text()
    s=re.sub(r'(PbI2_g|BiI3_g) ([0-9.e+-]+);',lambda m:m[1]+' '+format(float(m[2])*1e-4,'.17g')+';',s)
    s+='\nrespeciation {engine GEMS; species ('+' '.join(names)+f'); minTemperature 350; updateInterval 1; GEMS {{mode {mode}; system "{Path(system).resolve()}"; carrier "He(g)";}} }}\n'
    tc.write_text(s)
    ctrl=case/'system/controlDict';s=ctrl.read_text().replace('endTime           0.05;','endTime           0.01;').replace('writeInterval     5;','writeInterval     1;').replace('writeFormat       ascii;','writeFormat       binary;')
    ctrl.write_text(s.replace('tolerance          1e-12;','tolerance          1e-16;').replace('nCorr              3;','nCorr              10;'))
    e.close();return names

def compare(a,b,names,exact=False,tol=1e-9):
    for n in names:
        pa=a/('Y_'+n);pb=b/('Y_'+n)
        if exact:assert pa.read_bytes()==pb.read_bytes(),(pa,pb)
        else:
            va=sc.internal(pa);vb=sc.internal(pb);scale=max(max(va),max(vb),1e-50)
            assert max(abs(x-y) for x,y in zip(va,vb))<=tol*scale,(n,scale,max(abs(x-y) for x,y in zip(va,vb)))

def checkRestart(case,split):
    for group in ['phaseChangeElements','phaseChangeSpecies']:
        files=list((split/'postProcessing'/group).glob('*/*.dat'))
        assert files,(split,group)
        for p in files:
            for row in sc.rows(p):
                assert abs(row['closure'])<=1e-14*row['reference'],(p,row)
        for p in (case/'postProcessing'/group).glob('0/*.dat'):
            continuous=sc.rows(p)[-1]
            resumed=max((row for q in (split/'postProcessing'/group).glob('*/'+p.name)
                         for row in sc.rows(q)),key=lambda row:row['time'])
            scale=max(continuous['reference'],resumed['reference'],1e-300)
            inventories=['gas','wall','sample','held','supplied'] if group=='phaseChangeElements' else ['held','supplied']
            for key in inventories:
                assert abs(continuous[key]-resumed[key])<=1e-12*scale,(p,key,continuous,resumed)

def channel(work,solver,carrier,system,library):
    for mode in ['frozen','local']:
        for ranks in [1,4]:
            case=work/f'{mode}{ranks}';names=prepare(case,system,library,mode)
            base.run(case,solver,carrier,ranks);sc.check(case)
            if ranks>1:base.call(case,['reconstructPar','-latestTime'])
            split=work/f'{mode}restart{ranks}';prepare(split,system,library,mode)
            ctrl=split/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.005;'))
            base.run(split,solver,carrier,ranks)
            ctrl.write_text(ctrl.read_text().replace('endTime           0.005;','endTime           0.01;'))
            command=[solver] if ranks==1 else ['mpirun','-np',str(ranks),solver,'-parallel']
            base.call(split,command,'log.restart')
            checkRestart(case,split)
            prefixes=[Path()] if ranks==1 else [Path(f'processor{i}') for i in range(ranks)]
            for pre in prefixes:
                compare(case/pre/'0.01',split/pre/'0.01',names,exact=mode=='frozen',tol=1e-12)
                if mode=='frozen':
                    for f in ['phaseChangeSpecies','phaseChangeElements','phaseChangeCondensates','phaseChange','phaseChangeLayout']:
                        pa=case/pre/'0.01/uniform'/f;pb=split/pre/'0.01/uniform'/f
                        if pa.is_file():assert pa.read_bytes()==pb.read_bytes(),f
                        else:
                            for path in pa.rglob('*'):
                                if path.is_file():assert path.read_bytes()==(pb/path.relative_to(pa)).read_bytes(),path
        compare(work/f'{mode}1/0.01',work/f'{mode}4/0.01',names)
    compare(work/'frozen1/0.01',work/'local1/0.01',names)
    print('Frozen/local element closure <=1e-14; serial/4-rank fields <=1e-9; frozen binary restart exact; local restart <=1e-12')


def prime(case,carrier):
    base.call(case,['blockMesh']);base.call(case,[carrier])

def inputs(work,solver,carrier,system,library):
    mutations=[('mode frozen;','mode typo;','Invalid gas GEMS'),
      ('mode frozen;','mode frozen; maxMoleFraction -1;','Invalid gas GEMS'),
      ('mode frozen;','mode frozen; warmIterationLimit 0;','Invalid gas GEMS'),
      ('carrier "He(g)";','carrier "I(g)";','Carrier contains'),
      ('mode frozen;','mode frozen; maxMoleFracton 1;','Unknown gas GEMS')]
    for i,(a,b,error) in enumerate(mutations):
        c=work/f'input{i}';prepare(c,system,library,'frozen');prime(c,carrier)
        tc=c/'constant/thermochemistryProperties';tc.write_text(tc.read_text().replace(a,b))
        base.call(c,[solver],'log.refusal',fail=True);assert error in (c/'log.refusal').read_text(),c
    c=work/'changed';prepare(c,system,library,'frozen');base.run(c,solver,carrier)
    tc=c/'constant/thermochemistryProperties';tc.write_text(tc.read_text().replace('mode frozen;', 'mode local;'))
    base.call(c,[solver],'log.refusal',fail=True)
    assert 'configuration changed across restart' in (c/'log.refusal').read_text()
    # GEMS owns in-grid states even if a fallback CSV covers only another
    # temperature range. The fallback must not be evaluated eagerly.
    fallback=work/'narrowFallback';names=prepare(fallback,system,library,'frozen')
    (fallback/'constant/narrow.csv').write_text('# T,'+','.join(names)+'\n'+
        '\n'.join(str(t)+','+','.join('0' for _ in names) for t in [350,351])+'\n')
    tc=fallback/'constant/thermochemistryProperties'
    tc.write_text(tc.read_text().replace('engine GEMS;',
        'engine GEMS; formationConstants "<constant>/narrow.csv";'))
    base.run(fallback,solver,carrier)
    compare(c/'0.01',fallback/'0.01',names,exact=True)
    print('Invalid mode/guards/cap/carrier/unknown key and changed backend restart refused; in-grid GEMS bypasses narrow fallback CSV')

def hks(work,solver,carrier,system,library,data):
    import csv,math,ctypes as C
    tables=work/'tables';tables.mkdir(exist_ok=True)
    shutil.copy2(base.REPO/'thermochemistry/systems/data/pv_PbI2_phases.csv',tables/'pv_PbI2_phases.csv')
    shutil.copy2(Path(data)/'pv_BiI3_phases.csv',tables/'pv_BiI3_phases.csv');data=tables
    e=Engine(library,system);e.suppress();e.L.gemsb_set_warm_iteration_limit(e.h,200)
    def pressure(name,T):
        rows=[]
        with open(Path(data)/('pv_'+name+'_phases.csv')) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):rows.append(list(map(float,line.split(','))))
        for a,b in zip(rows,rows[1:]):
            if a[0]<=T<=b[0]:
                w=(1/T-1/a[0])/(1/b[0]-1/a[0]);return 10**min(x+w*(y-x) for x,y in zip(a[1:],b[1:]))
        raise AssertionError(T)
    worst=0;count=0
    for T in [400+25*i for i in range(32)]:
        rc=e.equilibrate(T,101325,{'He':1,'Pb':1e-4,'Bi':1e-4,'I':5.000005e-4})
        if rc not in (0,1):continue
        for n in ['PbI2','BiI3']:
            j=e.species.index(n+'(g)');cond=(C.c_int*2)(*(e.species.index(n+'('+p+')') for p in ['cr','l']))
            peq=C.c_double();best=C.c_int()
            assert e.L.gemsb_pair_peq(e.h,j,cond,2,1,C.byref(peq),C.byref(best))==0
            worst=max(worst,abs(math.log10(peq.value/pressure(n,T))));count+=1
    assert count>0 and worst<=1e-3,(count,worst);e.close()
    for mode in ['none','lagged']:
        c=work/f'hks_{mode}';base.prepare(c,['PbI2','BiI3'],external=Path(data))
        tc=c/'constant/thermochemistryProperties';s=tc.read_text()
        for n,formula in [('PbI2','Pb 1; I 2;'),('BiI3','Bi 1; I 3;')]:
            s=s.replace(f'{n}_g {{condensed {n}_s;', f'{n}_g {{condensed {n}_s; gems {{gas "{n}(g)"; condensates ("{n}(cr)" "{n}(l)"); elements {{{formula}}} excess {{I 1e-6;}}}}')
        s=s.replace('HKSCoeffs {','HKSCoeffs {equilibrium GEMS; speciation '+mode+';')
        s+='\nGEMSCoeffs {mode local; system "'+str(system)+'"; carrier "He(g)"; minTemperature 350; minMoleFraction 1e-12; maxMoleFraction 3.2e-4; allowNonPhysicalCarrier yes; writeEquilibriumField yes;}\n'
        tc.write_text(s)
        ctrl=c/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.05;','endTime           0.005;').replace('writeInterval     5;','writeInterval     1;'))
        base.run(c,solver,carrier)
        for path in (c/'postProcessing/phaseChangeElements/0').glob('*.dat'):
            for r in sc.rows(path):assert abs(r['closure'])<=1e-14*r['reference'],(path,r)
        # Recompute both the upper formula-unit guard and the chi gas
        # cutoff from the start-of-step fields and an independent cold call.
        e=Engine(library,system);e.suppress();checked=0;upper=0;minor=0
        def field(time,name,count):
            values=sc.internal(c/time/name);return values*count if len(values)==1 else values
        T=sc.internal(c/'0/T');nCells=len(T);P=field('0','p',nCells);rho=field('0','rho',nCells)
        previous='0'
        for time in ['0.001','0.002','0.003','0.004','0.005']:
            Y={n:field(previous,'Y_'+n+'_g',nCells) for n in ['PbI2','BiI3']}
            source={n:field(time,'gemsSource_'+n+'_g',nCells) for n in ['PbI2','BiI3']}
            for cell in range(nCells):
                if all(source[n][cell]<0 for n in source):continue
                cs={n:rho[cell]*Y[n][cell]/base.SPECIES[n][0] for n in Y};carrierAmount=P[cell]/(8.314462618*T[cell])
                bulk={'He':carrierAmount,'Pb':max(cs['PbI2'],0),'Bi':max(cs['BiI3'],0),
                      'I':(2*max(cs['PbI2'],0)+3*max(cs['BiI3'],0))*1.000001}
                raw=None
                for n in source:
                    code=source[n][cell]
                    if code<0:continue
                    x=cs[n]/(cs[n]+carrierAmount)
                    if x>3.2e-4:assert code==4,(mode,time,cell,n,x,code);upper+=1
                    elif x<1e-12 or cs[n]<=0:assert code==2,(mode,time,cell,n,x,code)
                    elif mode=='lagged':
                        if raw is None:
                            rc=e.equilibrate(T[cell],P[cell],bulk,0)
                            assert rc in (0,1),(mode,time,cell,rc)
                            raw=e.values('gas_partial_pressures')
                        xg=raw[e.species.index(n+'(g)')]/P[cell]
                        if xg<1e-7:assert code==4,(mode,time,cell,n,xg,code);minor+=1
                        else:assert code==0,(mode,time,cell,n,xg,code)
                    else:assert code==0,(mode,time,cell,n,x,code)
                    checked+=1
            previous=time
        assert checked>0 and upper>0 and (mode!='lagged' or minor>0),(mode,checked,upper,minor)
        e.close();print(mode,'guards checked',checked,'upper',upper,'minor chi',minor)
        log=(c/'log.solver').read_text();assert 'dual p_eq, warm iteration cap 200' in log
        assert 'maxMoleFraction 0.00032' in log
    print('Dual p_eq:',count,'values, max table deviation',worst,'decades; Bi HKS none/lagged solver paths and element closure passed')


def fault(work,solver,carrier,system,library,faultdir):
    import os
    reference=work/'reference';names=prepare(reference,system,library,'frozen');base.run(reference,solver,carrier)
    c=work/'fatal';prepare(c,system,library,'local');prime(c,carrier)
    values={'LD_LIBRARY_PATH':str(faultdir)+':'+os.environ.get('LD_LIBRARY_PATH',''),
      'LESTO_TEST_REAL_BRIDGE':str(Path(library).resolve()),'LESTO_TEST_FATAL_CALLS':'1 20',
      'LESTO_TEST_NAN_U_CALLS':'3','LESTO_TEST_BRIDGE_LOG':str(c/'fault.log')}
    old={k:os.environ.get(k) for k in values};os.environ.update(values)
    try:base.call(c,[solver],'log.solver')
    finally:
        for k,v in old.items():
            if v is None:os.environ.pop(k,None)
            else:os.environ[k]=v
    sc.check(c);compare(reference/'0.01',c/'0.01',names,exact=True)
    log=(c/'log.solver').read_text();failures=sum(map(int,re.findall(r'Gas GEMS local: calls \d+; iterations \d+; failures (\d+)',log)))
    assert failures>=3,(failures,log[-2000:])
    bridge=(c/'fault.log').read_text();assert bridge.count('fatal injected')==2 and bridge.count('nan potentials injected')==1
    print('Two fatal calls recreate the engine; NaN duals refused; all faults counted; fields byte-identical to frozen fallback')

if __name__=='__main__':
    action,work,solver,carrier,system,library,*extra=sys.argv[1:];work=Path(work);work.mkdir(parents=True,exist_ok=True)
    if action=='channel':channel(work,solver,carrier,system,library)
    elif action=='inputs':inputs(work,solver,carrier,system,library)
    elif action=='hks':hks(work,solver,carrier,system,library,*extra)
    elif action=='fault':fault(work,solver,carrier,system,library,*extra)
