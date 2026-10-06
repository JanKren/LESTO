#!/usr/bin/env python3
"""M10d conservative shared-store and reactive-wall code gates (synthetic data)."""
import importlib.util,math,re,sys,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('speciation_cases',HERE.parent/'speciation/cases.py')
sc=importlib.util.module_from_spec(spec);spec.loader.exec_module(sc)
base=sc.base
MASSES={'Pb':.2072,'Bi':.2089804,'I':.12690447}

def solid(case,name,atoms):
    mass=sum(MASSES[a]*n for a,n in atoms.items())
    props=case/'constant/speciesTransportProperties';s=props.read_text()
    s=re.sub(r'species \(([^;]+)\);',lambda m:'species ('+m[1]+' '+name+');',s)
    s+=f'\n{name} {{state solid; molarMass {mass:.17g}; formula {{'+''.join(f'{a} {n};' for a,n in atoms.items())+'}}\n'
    props.write_text(s)
    template=(case/'0/Y_BiI3_s').read_text();(case/'0'/('Y_'+name)).write_text(template.replace('BiI3_s',name))
    return mass

def prepare(case,T=700,reactive=False,evaporate=False,active=False,optional=False):
    sc.prepare(case,temperature=T,binary=True)
    solid(case,'Bi_s',{'Bi':1})
    p=case/'system/blockMeshDict';p.write_text(re.sub(r'\(100 24 1\)','(1 1 1)',p.read_text()))
    table=case/'constant/pv_metal.csv';logp=5 if evaporate else -30
    table.write_text(f'# SYNTHETIC code-test table\n# T,metal\n270,{logp}\n1250,{logp}\n')
    (case/'constant/logK.csv').write_text('# SYNTHETIC dimensionless reaction K; p0=1 bar\n# T,logK\n270,30\n1250,30\n')
    channels='''Bi {reactant Bi_g; condensed Bi_s; vapourPressure {file "<constant>/pv_metal.csv"; phases (metal); units log10Pa;}}
Bi2 {reactant Bi2_g; condensed Bi_s; stoichiometry 2; vapourPressure {file "<constant>/pv_metal.csv"; phases (metal); units log10Pa;}}
'''
    if reactive:channels='''BiI {reactant BiI_g; condensed Bi_s; stoichiometry {BiI_g -3; Bi_s 2; BiI3_g 1;} equilibrium {file "<constant>/logK.csv"; column logK;} reversible no;}'''
    if optional:
        solid(case,'BiI_s',{'Bi':1,'I':1})
        channels+='''BiI_plain {reactant BiI_g; condensed BiI_s; vapourPressure {file "<constant>/pv_metal.csv"; phases (metal); units log10Pa;}}'''
    p=case/'constant/thermochemistryProperties';s=p.read_text()
    # Legacy iodide kinetics are disabled; their existing accounting is exercised.
    s=s.replace('kineticScale 1;','kineticScale 0;').replace('minTemperature 350;','minTemperature '+('350' if active else '2000')+';')
    s+='\nwallChannels {\n'+channels+'}\n';p.write_text(s)
    for name in sc.NAMES:
        Y=0 if evaporate else (1e-6 if name in (['BiI_g'] if reactive else ['Bi_g','Bi2_g']) else 0)
        f=case/'0'/('Y_'+name);f.write_text(re.sub(r'internalField uniform [^;]+;',f'internalField uniform {Y:.17g};',f.read_text()))
    ctrl=case/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.001;'))
    if evaporate:
        mw=(case/'0/Y_BiI3_s').read_text().replace('Y_BiI3_s','mw_Bi_s').replace('[0 0 0 0 0 0 0]','[1 -2 0 0 0 0 0]')
        mw=re.sub(r'internalField uniform [^;]+;','internalField uniform 1e-12;',mw)
        mw=mw.replace('type zeroGradient;','type calculated; value uniform 1e-12;').replace('value uniform 0;','value uniform 1e-12;')
        (case/'0/mw_Bi_s').write_text(mw)
    return case

def closure(case):
    for group in ['phaseChangeElements','phaseChangeSpecies','phaseChangeChannelCondensates']:
        paths=list((case/'postProcessing'/group).glob('*/*.dat'));assert paths,(case,group)
        for p in paths:
            for r in sc.rows(p):assert abs(r['closure'])<=1e-14*max(r['reference'],1e-300),(p,r)

def value(case,name,time='0.001'):
    return sc.internal(case/time/('Y_'+name))[0]

def boxes(work,solver,carrier):
    c=prepare(work/'metals');base.run(c,solver,carrier);closure(c)
    T=700;dt=.001;av=2/.0048
    lost=0
    for name,mass in [('Bi_g',MASSES['Bi']),('Bi2_g',2*MASSES['Bi'])]:
        lam=av*2*math.sqrt(8.314462618*T/(2*math.pi*mass))
        exact=1e-6/(1+dt*lam);actual=value(c,name)
        assert abs(actual/exact-1)<=1e-14,(name,actual,exact)
        lost+=1e-6-actual
    assert abs(value(c,'Bi_s')/lost-1)<=1e-14
    c=prepare(work/'reactive',500,reactive=True);base.run(c,solver,carrier);closure(c)
    lost=1e-6-value(c,'BiI_g');M=MASSES['Bi']+MASSES['I']
    assert abs(value(c,'Bi_s')/MASSES['Bi']/(lost/M)-2/3)<=1e-15
    assert abs(value(c,'BiI3_g')/.58969381/(lost/M)-1/3)<=1e-15
    # The product does not traverse its own wall channel until the next step.
    assert value(c,'BiI3_s')==0
    # Nonzero lagged product pressure exercises the dimensionless K and 1 bar
    # conversion, independently of the zero-product supersaturation box.
    c=prepare(work/'laggedPressure',500,reactive=True)
    p=c/'0/Y_BiI3_g';p.write_text(re.sub(r'internalField uniform [^;]+;','internalField uniform 1e-6;',p.read_text()))
    base.run(c,solver,carrier);closure(c)
    M=MASSES['Bi']+MASSES['I'];p3=101325*.0040026*1e-6/.58969381
    pstar=1e5*(p3/1e5/1e30)**(1/3);Yeq=pstar*M/(101325*.0040026)
    lam=av*2*math.sqrt(8.314462618*500/(2*math.pi*M))
    assert abs(value(c,'BiI_g')/((1e-6+dt*lam*Yeq)/(1+dt*lam))-1)<=1e-14
    # Product-rich gas cannot reverse the irreversible channel into evaporation.
    c=prepare(work/'irreversible',500,reactive=True)
    p=c/'0/Y_BiI_g';p.write_text(p.read_text().replace('internalField uniform 9.9999999999999995e-07;','internalField uniform 0;').replace('internalField uniform 1e-06;','internalField uniform 0;'))
    p=c/'0/Y_BiI3_g';p.write_text(re.sub(r'internalField uniform [^;]+;','internalField uniform 1e-6;',p.read_text()))
    base.run(c,solver,carrier);closure(c);assert value(c,'Bi_s')==0 and value(c,'BiI_g')==0
    for order in [0,1]:
        c=prepare(work/f'capped{order}',evaporate=True,active=True)
        if order:
            p=c/'constant/speciesTransportProperties';s=p.read_text();s=s.replace('Bi_g Bi2_g','Bi2_g Bi_g');p.write_text(s)
        base.run(c,solver,carrier);closure(c)
        assert value(c,'Bi_s')==0,(c,value(c,'Bi_s'))
        assert 'CAPPED 2' in (c/'log.solver').read_text()
    for name in sc.NAMES:
        a=value(work/'capped0',name);b=value(work/'capped1',name)
        assert abs(a-b)<=1e-3*max(a,b,1e-100),(name,a,b)
    c=prepare(work/'optional',500,reactive=True,optional=True);base.run(c,solver,carrier);closure(c)
    assert value(c,'Bi_s')>0 and value(c,'BiI_s')>0
    c=prepare(work/'evaporationWithoutRespeciation',evaporate=True)
    base.run(c,solver,carrier);closure(c);assert value(c,'Bi_s')==0
    for diffusion in [0,1e-4]:
        c=prepare(work/('halfCell'+str(diffusion)))
        p=c/'constant/thermochemistryProperties';p.write_text(p.read_text().replace('patches (WALL);','patches (WALL); wallResistance halfCell;'))
        p=c/'constant/speciesTransportProperties';p.write_text(p.read_text().replace('D 0;',f'D {diffusion};'))
        for name in sc.NAMES:
            p=c/'0'/('Y_'+name);p.write_text(p.read_text().replace('type fixedValue; value uniform 0;','type zeroGradient;'))
        base.run(c,solver,carrier);closure(c)
        for name,mass in [('Bi_g',MASSES['Bi']),('Bi2_g',2*MASSES['Bi'])]:
            v=2*math.sqrt(8.314462618*700/(2*math.pi*mass))
            lam=0 if diffusion==0 else av*v/(1+v*.0024/diffusion)
            assert abs(value(c,name)/(1e-6/(1+dt*lam))-1)<=1e-14
    print('Closed forms <=1e-14; reactive 2/3+1/3 <=1e-15; irreversible; shared CAPPED exact zero and order gate; optional BiI(cr); closure <=1e-14')

def inputs(work,solver,carrier):
    mutations={
      'dimer':lambda s:s.replace('stoichiometry 2;','stoichiometry 1;'),
      'unknown':lambda s:s.replace('reactant Bi_g;','reactant Bi_g; typo 1;'),
      'kinetics':lambda s:s.replace('reactant Bi_g;','reactant Bi_g; HKS {accommodation 2;}'),
      'empty':lambda s:s[:s.index('\nwallChannels')]+ '\nwallChannels {}\n',
      'engine':lambda s:s.replace('engine kernel;','engine none;'),
      'productNotSelected':lambda s:s.replace('BiI3_g I_g','I_g'),
      'unbalanced':lambda s:s.replace('Bi_s 2;','Bi_s 3;'),
      'reversible':lambda s:s.replace('reversible no;','reversible yes;'),
      'wrongUnits':lambda s:s.replace('column logK;','column logK; units Pa;'),
    }
    for name,change in mutations.items():
        c=prepare(work/('invalid_'+name),reactive=name in ['productNotSelected','unbalanced','reversible','wrongUnits'])
        p=c/'constant/thermochemistryProperties';p.write_text(change(p.read_text()))
        base.call(c,['blockMesh']);base.call(c,[carrier]);base.call(c,[solver],'log.refusal',fail=True)
        assert 'FOAM FATAL' in (c/'log.refusal').read_text(),name
    for name in ['table','kinetics','removed','nonfiniteAccount']:
        c=prepare(work/('restartChange_'+name));base.run(c,solver,carrier)
        p=c/'constant/thermochemistryProperties'
        if name=='table':
            p=c/'constant/pv_metal.csv';p.write_text(p.read_text().replace('1250,-30','1250,-29'))
        elif name=='kinetics':p.write_text(p.read_text().replace('reactant Bi_g;','reactant Bi_g; HKS {Ce .5;}'))
        elif name=='removed':p.write_text(p.read_text().split('\nwallChannels')[0])
        else:
            import struct
            p=c/'0.001/uniform/phaseChangeChannelCondensates';b=bytearray(p.read_bytes())
            offset=re.search(rb'\n10\n\(',b).end();b[offset:offset+8]=struct.pack('d',float('nan'));p.write_bytes(b)
        base.call(c,[solver],'log.refusal',fail=True)
        assert 'FOAM FATAL' in (c/'log.refusal').read_text(),name
    print('Invalid channel formulas/keys/kinetics/products/units and changed restart configuration/data refused')

def flowPrepare(case):
    prepare(case,active=True)
    solid(case,'Pb_s',{'Pb':1})
    p=case/'system/blockMeshDict';p.write_text(p.read_text().replace('(1 1 1)','(20 4 1)',1))
    p=case/'system/carrierDict';s=p.read_text();s=re.sub(r'Tin\s+[^;]+;','Tin 1000;',s);s=re.sub(r'Tout\s+[^;]+;','Tout 400;',s);s=re.sub(r'massFlux\s+[^;]+;','massFlux .01743;',s);p.write_text(s)
    p=case/'constant/speciesTransportProperties';p.write_text(p.read_text().replace('D 0;','D 1e-4;'))
    p=case/'constant/thermochemistryProperties';s=p.read_text().replace('wallChannels {','wallChannels { Pb {reactant Pb_g; condensed Pb_s; vapourPressure {file "<constant>/pv_metal.csv"; phases (metal); units log10Pa;}} BiI {reactant BiI_g; condensed Bi_s; stoichiometry {BiI_g -3; Bi_s 2; BiI3_g 1;} equilibrium {file "<constant>/logK.csv"; column logK;} reversible no;}')
    s+='\ngasSources {boat {amounts {Pb_g 2e-12; Bi_g 3e-12; I_g 4e-12;} box ((.01 0 -1) (.03 .0048 1)); startTime 0; duration .006;}}\n';p.write_text(s)
    p=case/'system/controlDict';p.write_text(p.read_text().replace('endTime           0.001;','endTime           0.01;'))
    for name in sc.NAMES:
        p=case/'0'/('Y_'+name);p.write_text(re.sub(r'internalField uniform [^;]+;','internalField uniform 0;',p.read_text()))

def channel(work,solver,carrier):
    c=work/'legacySamples';flowPrepare(c)
    p=c/'constant/thermochemistryProperties';p.write_text(p.read_text().replace('samples {}','samples {boat {amounts {PbI2_g 1e-12; BiI3_g 1e-12;} selection {boat {action use; source box; box (.01 .0012 -1) (.03 .0036 1);}} mode release; startTime 0; duration .006;}}'))
    base.run(c,solver,carrier);closure(c)
    for ranks in [1,4]:
        c=work/f'flow{ranks}';flowPrepare(c);base.run(c,solver,carrier,ranks);closure(c)
        split=work/f'restart{ranks}';flowPrepare(split);ctrl=split/'system/controlDict';ctrl.write_text(ctrl.read_text().replace('endTime           0.01;','endTime           0.005;'));base.run(split,solver,carrier,ranks)
        ctrl.write_text(ctrl.read_text().replace('endTime           0.005;','endTime           0.01;'));base.call(split,[solver] if ranks==1 else ['mpirun','-np',str(ranks),solver,'-parallel'],'log.restart');closure(split)
        prefixes=[Path()] if ranks==1 else [Path(f'processor{i}') for i in range(ranks)]
        for pre in prefixes:
            for name in sc.NAMES+['Pb_s','Bi_s']:
                a=c/pre/'0.01'/('Y_'+name);b=split/pre/'0.01'/('Y_'+name);assert a.read_bytes()==b.read_bytes(),(a,b)
            for file in ['phaseChangeSpecies','phaseChangeElements','phaseChangeChannelCondensates','phaseChangeWallChannels']:
                assert (c/pre/'0.01/uniform'/file).read_bytes()==(split/pre/'0.01/uniform'/file).read_bytes(),file
        if ranks>1:base.call(c,['reconstructPar','-latestTime'])
    for name in sc.NAMES+['Pb_s','Bi_s']:
        a=sc.internal(work/'flow1/0.01'/('Y_'+name));b=sc.internal(work/'flow4/0.01'/('Y_'+name));scale=max(max(a),max(b),1e-100);assert max(abs(x-y) for x,y in zip(a,b))<=1e-9*scale,name
    for atom,nu in [('Pb',2e-12),('Bi',3e-12),('I',4e-12)]:
        r=sc.rows(work/f'flow1/postProcessing/phaseChangeElements/0/{atom}.dat')[-1];assert abs(r['released']/nu-1)<=1e-14,(atom,r)
    # Decompose the serial mid-run state: the exact local geometry no longer
    # matches, so the shared deposits must be recovered from the mw_ fields.
    c=work/'changedDecomposition'
    if c.exists():shutil.rmtree(c)
    c.mkdir()
    for name in ['0','0.005','constant','system']:shutil.copytree(work/'restart1'/name,c/name,dirs_exist_ok=True)
    base.call(c,['decomposePar','-force','-time','0,0.005'])
    base.call(c,['mpirun','-np','4',solver,'-parallel'],'log.restart');closure(c)
    base.call(c,['reconstructPar','-latestTime'])
    for name in sc.NAMES+['Pb_s','Bi_s']:
        a=sc.internal(work/'flow1/0.01'/('Y_'+name));b=sc.internal(c/'0.01'/('Y_'+name));scale=max(max(a),max(b),1e-100)
        assert max(abs(x-y) for x,y in zip(a,b))<=1e-9*scale,(name,'changed decomposition')
    for case in [work/'flow1',work/'flow4',c]:
        for solidName,mass in [('Pb_s',MASSES['Pb']),('Bi_s',MASSES['Bi'])]:
            s=(case/f'postProcessing/axialProfiles/0.01/condensate_{solidName}_native.dat').read_text()
            rows=[list(map(float,l.split())) for l in s.splitlines() if l and not l.startswith('#')]
            integrated=math.fsum((r[1]-r[0])*r[4] for r in rows)
            held=sc.rows(next((case/'postProcessing/phaseChangeChannelCondensates').glob(f'*/{solidName}.dat')))[-1]['held']/mass
            assert abs(integrated-held)<=1e-14*max(held,1e-100),(case,solidName,integrated,held)
    print('Q2 monatomic source amounts exact; all ledgers <=1e-14; serial/4-rank <=1e-9; fields and accounts restart byte-identical')

if __name__=='__main__':
    work=Path(sys.argv[2]);work.mkdir(parents=True,exist_ok=True)
    globals()[sys.argv[1]](work,*sys.argv[3:])
