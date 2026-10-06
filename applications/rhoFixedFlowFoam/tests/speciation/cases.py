#!/usr/bin/env python3
"""M10b synthetic integration gates; all thermodynamic constants are invented."""
import importlib.util
import math
from pathlib import Path
import re
import shutil
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'multiSpecies'))
import cases as base
sys.path.insert(0,str(base.REPO/'thermochemistry/systems'))
from make_kf_tables import SPECIES, write, synthetic
from checks import readScalarField
def internal(path): return readScalarField(path)[0]
MASSES={'Pb':.2072,'Bi':.2089804,'I':.12690447}
NAMES=list(SPECIES)

def prepare(case, engine='kernel', temperature=None, binary=False):
    base.prepare(case,['PbI2','BiI3'])
    props=case/'constant/speciesTransportProperties'
    s=props.read_text(); gases=[n for n in NAMES if n not in ['PbI2_g','BiI3_g']]
    s=re.sub(r'species \([^;]+;', 'species ('+' '.join(['PbI2_g','PbI2_s','BiI3_g','BiI3_s']+gases)+');',s)
    template=(case/'0/Y_PbI2_g').read_text()
    for n in gases:
        f=SPECIES[n];M=sum(MASSES[e]*v for e,v in f.items())
        s+=f'\n{n} {{state gas; molarMass {M:.17g}; diffusivityModel constant; D 1e-4; formula {{'+''.join(f'{e} {v};' for e,v in f.items())+'}}\n'
        (case/'0'/('Y_'+n)).write_text(template.replace('PbI2_g',n))
    props.write_text(s)
    write(case/'constant/log10Kf.csv',NAMES,synthetic)
    tc=case/'constant/thermochemistryProperties';s=tc.read_text()
    s+=f'\nrespeciation {{engine {engine}; species ('+' '.join(NAMES)+'); formationConstants "<constant>/log10Kf.csv"; minTemperature 350; updateInterval 1;}\n'
    s=re.sub(r'(PbI2_g|BiI3_g) ([0-9.e+-]+);',lambda m:m[1]+' '+format(float(m[2])*1e-4,'.17g')+';',s)
    tc.write_text(s)
    sol=case/'system/fvSolution';s=sol.read_text();position=s.rfind('}')
    s=s[:position]+''.join(f'Y_{n} {{solver PBiCGStab; preconditioner DILU; tolerance 1e-18; relTol 0;}}\n' for n in gases)+s[position:];sol.write_text(s.replace('tolerance 1e-15;', 'tolerance 1e-18;'))
    schemes=case/'system/fvSchemes';s=schemes.read_text().replace('default none;','default none; '+''.join(f'div(phi,Y_{n}) Gauss upwind; ' for n in gases));schemes.write_text(s)
    ctrl=case/'system/controlDict';s=ctrl.read_text().replace('endTime           0.05;','endTime           0.01;').replace('writeInterval     5;','writeInterval     1;')
    s=s.replace('tolerance          1e-12;','tolerance          1e-16;').replace('nCorr              3;','nCorr              10;')
    if binary:s=s.replace('writeFormat       ascii;','writeFormat       binary;')
    ctrl.write_text(s)
    if temperature is not None:
        tc.write_text(re.sub(r'^samples [^\n]*$', 'samples {}',tc.read_text(),flags=re.M).replace('Ce .0001;','Ce 0;').replace('Ce .0002;','Ce 0;'))
        carrier=case/'system/carrierDict';s=carrier.read_text();s=re.sub(r'Tin\s+\d+;',f'Tin {temperature};',s);s=re.sub(r'Tout\s+\d+;',f'Tout {temperature};',s);s=re.sub(r'massFlux\s+[0-9.]+;','massFlux 0;',s);carrier.write_text(s)
        props.write_text(props.read_text().replace('D 1e-4;','D 0;'))
        for n in ['Pb_g','Bi_g','I_g']:
            field=case/'0'/('Y_'+n);field.write_text(field.read_text().replace('internalField uniform 0;',f'internalField uniform {1e-7 if n!="I_g" else 3e-7};').replace('type fixedValue; value uniform 0;','type zeroGradient;'))
    return case

def rows(path):
    s=path.read_text();keys=re.search(r'# columns (.*)',s)[1].split()
    return [dict(zip(keys,map(float,l.split()))) for l in s.splitlines() if l and not l.startswith('#')]

def check(case):
    worst=0
    for p in (case/'postProcessing/phaseChangeElements/0').glob('*.dat'):
        for r in rows(p):
            rel=abs(r['closure'])/max(r['reference'],1e-300);worst=max(worst,rel)
            assert rel<=1e-14,(p,r,rel)
            assert abs(r['reaction'])<=1e-14*max(r['reference'],1e-300),(p,r)
    for p in (case/'postProcessing/phaseChangeSpecies/0').glob('*.dat'):
        for r in rows(p):assert abs(r['closure'])<=1e-14*r['reference'],(p,r)
    print(case.name,'max element closure',worst)

def closed(work,solver,carrier):
    for T,split in [(450,'atoms'),(450,'iodides'),(1000,'atoms'),(1000,'iodides')]:
        c=prepare(work/f'box{T}_{split}',temperature=T)
        if split=='iodides':
            amounts={'Pb_g':0,'Bi_g':0,'PbI2_g':1e-7/MASSES['Pb']*.46100894,
                     'BiI3_g':1e-7/MASSES['Bi']*.58969381,
                     'I_g':3e-7-(2e-7/MASSES['Pb']+3e-7/MASSES['Bi'])*MASSES['I']}
            # Use only part of each metal in iodides so iodine remains positive.
            for n in ['PbI2_g','BiI3_g']:amounts[n]*=.25
            amounts['Pb_g']=.75e-7;amounts['Bi_g']=.75e-7
            amounts['I_g']=3e-7-.25*(2e-7/MASSES['Pb']+3e-7/MASSES['Bi'])*MASSES['I']
            for n,Y in amounts.items():
                p=c/'0'/('Y_'+n);p.write_text(re.sub(r'internalField uniform [^;]+;',f'internalField uniform {Y:.17g};',p.read_text()).replace('type fixedValue; value uniform 0;','type zeroGradient;'))
        base.run(c,solver,carrier);check(c)
        for p in (c/'postProcessing/phaseChangeElements/0').glob('*.dat'):
            for r in rows(p):
                assert abs(r['held']/r['initial']-1)<=1e-15,(p,r)
                assert abs(r['reaction'])<=1e-15*r['initial'],(p,r)
        first={n:internal(c/'0.001'/('Y_'+n)) for n in NAMES}
        for n in NAMES:
            assert min(first[n])>=0
            last=internal(c/'0.01'/('Y_'+n))
            scale=max(max(v) for v in first.values())
            assert max(abs(a-b) for a,b in zip(first[n],last))/scale<1e-13,(T,n)
        for atom in MASSES:
            initial=sum((1e-7 if n!='I_g' else 3e-7)/MASSES[atom] for n in [atom+'_g'])
            # molar concentrations share the same rho before and after.
            actual=sum(SPECIES[n].get(atom,0)*first[n][0]/sum(MASSES[e]*v for e,v in SPECIES[n].items()) for n in NAMES)
            assert abs(actual/initial-1)<1e-15,(T,atom,actual,initial)

def channel(work,solver,carrier):
    for ranks in [1,4]:
        c=prepare(work/f'channel{ranks}',binary=True)
        base.run(c,solver,carrier,ranks);check(c)
        if ranks>1:base.call(c,['reconstructPar','-latestTime'])
    import subprocess
    program=work/'testSpeciationFields'
    assert subprocess.run(['g++','-O2','-std=c++17',str(HERE/'testSpeciation.C'),'-o',str(program)]).returncode==0
    equilibrium(work/'channel1','0.01',program)
    equilibrium(work/'channel4','0.01',program)
    # Compare serial/parallel converged fields numerically (OpenFOAM utility
    # writes binary: readInternal handles both formats).
    for n in NAMES:
        a=internal(work/'channel1/0.01'/('Y_'+n));b=internal(work/'channel4/0.01'/('Y_'+n))
        scale=max(max(a),max(b),1e-30)
        assert max(abs(x-y) for x,y in zip(a,b))/scale<1e-9,n
    for ranks in [1,4]:
        full=work/f'channel{ranks}';split=work/f'restart{ranks}';prepare(split,binary=True)
        ctrl=split/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.005;'))
        base.run(split,solver,carrier,ranks)
        ctrl.write_text(ctrl.read_text().replace('endTime           0.005;','endTime           0.01;'))
        command=[solver] if ranks==1 else ['mpirun','-np',str(ranks),solver,'-parallel']
        base.call(split,command,'log.restart')
        prefixes=[Path()] if ranks==1 else [Path(f'processor{i}') for i in range(4)]
        for pre in prefixes:
            for path in (full/pre/'0.01').rglob('*'):
                if path.is_file() and (path.name.startswith(('Y_','mw_','mDep_','phaseChange'))):
                    other=split/pre/'0.01'/path.relative_to(full/pre/'0.01')
                    assert path.read_bytes()==other.read_bytes(),(ranks,path,other)
    print('serial/4-rank fields agree; same-decomposition binary restart exact')

def disabled(work,solver,carrier):
    for mode,T in [('none',None),('kernel',300)]:
        c=prepare(work/f'disabled_{mode}',engine=mode)
        if T is not None:
            carrierDict=c/'system/carrierDict';content=carrierDict.read_text()
            content=re.sub(r'Tin\s+\d+;',f'Tin {T};',content);content=re.sub(r'Tout\s+\d+;',f'Tout {T};',content);carrierDict.write_text(content)
            iodine=c/'0/Y_I_g';iodine.write_text(iodine.read_text().replace('type fixedValue; value uniform 0;','type fixedValue; value uniform 1e-8;'))
        control=work/f'control_{mode}'
        if control.exists():shutil.rmtree(control)
        shutil.copytree(c,control)
        tc=control/'constant/thermochemistryProperties'
        content=tc.read_text()
        tc.write_text(content.split('\nrespeciation ')[0]+'\n' if mode=='none' else content.replace('engine kernel;','engine none;'))
        base.run(c,solver,carrier);base.run(control,solver,carrier)
        for p in (c/'0.01').glob('Y_*'):assert p.read_bytes()==(control/'0.01'/p.name).read_bytes(),p
    print('engine none and below-minTemperature fields exact')

def equilibrium(case,time,program):
    import subprocess
    table=[[float(v) for v in l.split(',')] for l in (case/'constant/log10Kf.csv').read_text().splitlines() if l and not l.startswith('#')]
    names=[l[1:].strip().split(',')[1:] for l in (case/'constant/log10Kf.csv').read_text().splitlines() if l.startswith('#')][-1]
    T=internal(case/'0/T');ncells=len(T)
    concentrations={n:internal(case/time/('c_'+n)) for n in NAMES}
    concentrations={n:v*ncells if len(v)==1 else v for n,v in concentrations.items()}
    target=case/'equilibrium.csv'
    with target.open('w') as out:
        out.write('# written-field equilibrium; T,scale,bPb,bBi,bI,(lnKf,c)\n')
        for cell,temperature in enumerate(T):
            if temperature<350:continue
            b=[math.fsum(SPECIES[n].get(e,0)*concentrations[n][cell] for n in NAMES) for e in ['Pb','Bi','I']]
            for lo,hi in zip(table,table[1:]):
                if lo[0]<=temperature<=hi[0]:break
            w=(1/temperature-1/lo[0])/(1/hi[0]-1/lo[0])
            row=[temperature,1e5/(8.314462618*temperature),*b]
            for n in NAMES:
                column=names.index(n)+1
                k=math.log(math.pow(10,(1-w)*lo[column]+w*hi[column]))
                row.extend([k,concentrations[n][cell]])
            out.write(','.join(format(v,'.17g') for v in row)+'\n')
    assert subprocess.run([str(program),str(target),'--fields']).returncode==0,target

def inputs(work,solver,carrier):
    mutations = {
      'engine': lambda s:s.replace('engine kernel;', 'engine typo;'),
      'duplicate':lambda s:s.replace('species (Pb_g', 'species (Pb_g Pb_g'),
      'noMonatomic':lambda s:s.replace('species (Pb_g ', 'species ('),
      'interval':lambda s:s.replace('updateInterval 1;', 'updateInterval 0;'),
      'lagged':lambda s:s+'\nHKSCoeffs {equilibrium table; speciation lagged;}\n',
    }
    for name,change in mutations.items():
        c=prepare(work/('invalid_'+name));p=c/'constant/thermochemistryProperties';p.write_text(change(p.read_text()))
        base.call(c,['blockMesh']);base.call(c,[carrier]);base.call(c,[solver],'log.refusal',fail=True)
    for name,change in {
        'changedEngine':lambda s:s.replace('engine kernel;','engine none;'),
        'changedInterval':lambda s:s.replace('updateInterval 1;','updateInterval 2;'),
        'changedSpecies':lambda s:s.replace('Pb_g PbI_g','PbI_g Pb_g'),
    }.items():
        c=prepare(work/('restart_input_'+name));base.run(c,solver,carrier)
        p=c/'constant/thermochemistryProperties';p.write_text(change(p.read_text()))
        ctrl=c/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.011;'))
        base.call(c,[solver],'log.refusal',fail=True)
        assert 'across restart' in (c/'log.refusal').read_text(),name
    c=prepare(work/'restart_mass');base.run(c,solver,carrier)
    props=c/'constant/speciesTransportProperties';props.write_text(props.read_text().replace('molarMass 0.2072;', 'molarMass 0.2072000001;'))
    ctrl=c/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.011;'))
    base.call(c,[solver],'log.refusal',fail=True)
    assert 'across restart' in (c/'log.refusal').read_text()
    for name in ['phaseChangeSpecies','phaseChangeElements','phaseChangeCondensates','phaseChangeBalance']:
        c=prepare(work/('missing_'+name));base.run(c,solver,carrier)
        (c/'0.01/uniform'/name).unlink()
        ctrl=c/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.011;'))
        base.call(c,[solver],'log.refusal',fail=True)
        if name!='phaseChangeBalance':assert 'Missing or invalid' in (c/'log.refusal').read_text(),name
    print('invalid engine/species/schedule and restart changes refused')

if __name__=='__main__':
    mode,work,solver,carrier=sys.argv[1:];work=Path(work);work.mkdir(parents=True,exist_ok=True)
    globals()[mode](work,solver,carrier)
