#!/usr/bin/env python3
"""M10a case preparation and checks. Stdlib; data remain outside the repository."""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
TESTS = HERE.parent
REPO = HERE.parents[3]
sys.path.insert(0, str(TESTS/'phaseChange'))
from checks import readBalance, readInternal

SPECIES = {
    'PbI2': (.46100894, {'Pb': 1, 'I': 2}),
    'PbI2copy': (.46100894, {'Pb': 1, 'I': 2}),
    'BiI3': (.58969381, {'Bi': 1, 'I': 3}),
}
HEADER = 'FoamFile {version 2.0; format ascii; class dictionary; object %s;}\n'


def call(case, command, log=None, fail=False):
    log = log or 'log.'+Path(command[0]).name
    with open(case/log, 'w') as f:
        rc = subprocess.run(command, cwd=case, stdout=f, stderr=subprocess.STDOUT).returncode
    if fail:
        assert rc != 0, f'Expected refusal: {case}'
    else:
        assert rc == 0, f'{command} failed; see {case/log}'


def prepare(case, names, model='HKS', mode='release', external=None, identical=False):
    if case.exists(): shutil.rmtree(case)
    shutil.copytree(TESTS/'phaseChange/channel', case)
    pairs, props, amounts, solvers, divs = [], [], [], [], []
    for name in names:
        M, atoms = SPECIES[name]
        formula = 'formula {'+' '.join(f'{e} {n};' for e,n in atoms.items())+'}'
        diffusion = 'diffusivityModel constant; D 1e-4;'
        props.append(f'{name}_g {{state gas; {diffusion} molarMass {M}; {formula}}}')
        props.append(f'{name}_s {{state solid; molarMass {M}; {formula}}}')
        base = 'PbI2' if name == 'PbI2copy' else name
        for state in ['g','s']:
            template = (case/'0'/f'Y_PbI2_{state}').read_text()
            (case/'0'/f'Y_{name}_{state}').write_text(template.replace(f'PbI2_{state}', f'{name}_{state}'))
        phases = f'"{base}(cr)" "{base}(l)"'
        table = f'pv_{name}.csv'
        if external:
            shutil.copy2(external/f'pv_{base}_phases.csv', case/'constant'/table)
        else:
            (case/'constant'/table).write_text(
                f'# SYNTHETIC M10a code-test data; no experimental claim\n# T,{base}(cr),{base}(l)\n'
                '270,-16,-12\n350,-12,-9\n550,-5,-4\n680.85,-1,-1\n800,1,0.5\n1000,3,2\n1250,5,4\n')
        kinetics = 'accommodation .8; Ce .0001; kineticScale 1;'
        if base == 'BiI3' and not identical: kinetics = 'accommodation .6; Ce .0002; kineticScale .5;'
        hks = f'HKS {{{kinetics}}} vapourPressure {{file "<constant>/{table}"; phases ({phases}); units log10Pa; outOfRange fatal;}}' if model == 'HKS' else ''
        pairs.append(f'{name}_g {{condensed {name}_s; {hks}}}')
        amount = 4.84e-8 if base == 'PbI2' else (1e-12 if mode == 'inventory' else 6.17e-8)
        amounts.append(f'{name}_g {amount};')
        solvers.append(f'Y_{name}_g {{solver PBiCGStab; preconditioner DILU; tolerance 1e-15; relTol 0;}}')
        solvers.append(f'Y_{name}_s {{solver diagonal;}}')
        divs.append(f'div(phi,Y_{name}_g) Gauss upwind;')
    for old in ['PbI2_g','PbI2_s']:
        if old.split('_')[0] not in names: (case/'0'/f'Y_{old}').unlink()
    (case/'constant/speciesTransportProperties').write_text(
        HEADER%'speciesTransportProperties'+'species ('+' '.join(f'{n}_{s}' for n in names for s in ['g','s'])+');\n'+'\n'.join(props)+'\n')
    schedule = 'startTime .0003; duration .0405;' if mode == 'release' else 'areaPerVolume 4000; removeTime .045;'
    sampleGroups = f"boat {{amounts {{{' '.join(amounts)}}} selection {{boat {{action use; source box; box (.01 .002 -1) (.02 .0028 1);}}}} mode {mode}; {schedule}}}"
    if 'PbI2copy' in names:
        # Identical physics in a second physical boat. wallPlusGas must
        # exclude only each pair's own sample cells, in serial and MPI.
        groups=[]
        for name,amount in zip(names,amounts):
            position='(.03 .002 -1) (.04 .0028 1)' if name=='PbI2copy' else '(.01 .002 -1) (.02 .0028 1)'
            value=amount.split()[1]
            groups.append(f"boat__{name}_g {{pair {name}_g; amount {value} selection {{boat {{action use; source box; box {position};}}}} mode {mode}; {schedule}}}")
        sampleGroups=' '.join(groups)
    tc = HEADER%'thermochemistryProperties'+f'''
model {model};
pairs {{ {' '.join(pairs)} }}
interface {{patches (WALL);}}
HKSCoeffs {{accommodation 1; Ce 1; kineticScale 1;}}
temperatureCoeffs {{A 220; k 8e-3; Tdep 680;}}
samples {{{sampleGroups}}}
balance {{nGuard 5; guardTolerance 1e-13; writeExchangeField yes;}}
profiles {{axis (1 0 0); origin (0 0 0); observable wallPlusGas;
 bins {{native {{min 0; max .1; nBins 100;}}}}
 onset {{fraction .01; searchFrom .02;}}
}}
'''
    (case/'constant/thermochemistryProperties').write_text(tc)
    ctrl=(case/'system/controlDict').read_text().replace('endTime           0.2;', 'endTime           0.05;').replace('writeInterval     50;', 'writeInterval     5;').replace('writeFormat       binary;', 'writeFormat       ascii;').replace('tolerance          1e-7;', 'tolerance          1e-12;')
    if mode == 'inventory': ctrl = ctrl.replace('writeInterval     5;', 'writeInterval     1;')
    ctrl += '\nwritePrecision 17;\n'
    (case/'system/controlDict').write_text(ctrl)
    (case/'system/carrierDict').write_text((case/'system/carrierDict').read_text().replace('Tout      400;', 'Tout      350;'))
    (case/'system/fvSolution').write_text(HEADER%'fvSolution'+'solvers {'+' '.join(solvers)+'}\n')
    (case/'system/fvSchemes').write_text(HEADER%'fvSchemes'+'''ddtSchemes {default Euler;}
gradSchemes {default Gauss linear;}
divSchemes {default none; '''+' '.join(divs)+'''}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;}
snGradSchemes {default corrected;}
''')
    (case/'system/decomposeParDict').write_text(HEADER%'decomposeParDict'+'numberOfSubdomains 4; method simple; simpleCoeffs {n (4 1 1); delta .001;}\n')
    return case


def run(case, solver, carrier, ranks=1):
    call(case,['blockMesh'])
    call(case,[carrier])
    if 'mode inventory;' in (case/'constant/thermochemistryProperties').read_text():
        call(case,['postProcess','-func','writeCellCentres','-time','0'])
    if ranks > 1:
        call(case,['decomposePar','-force'])
        call(case,['mpirun','-np',str(ranks),solver,'-parallel'], 'log.solver')
    else: call(case,[solver], 'log.solver')


def balance(case, name):
    path=next((case/'postProcessing/phaseChangeBalance/0').glob(name+'_g.dat'))
    keys,rows=readBalance(path)
    return [dict(zip(keys,map(float,r))) for r in rows]


def elements(case, names):
    ledgers={n:balance(case,n) for n in names}
    for path in (case/'postProcessing/phaseChangeElements/0').glob('*.dat'):
        text=path.read_text(); keys=re.search(r'# columns (.*)',text)[1].split()
        rows=[dict(zip(keys,map(float,l.split()))) for l in text.splitlines() if l and not l.startswith('#')]
        for k,row in enumerate(rows):
            pairs=[(SPECIES[n][1].get(path.stem,0)/SPECIES[n][0],ledgers[n][k]) for n in names]
            for key in ['gas','wall','sample','initial','released','removed','clamped','restart','solverDefect','closure']:
                if key == 'closure': continue # subtraction round-off assessed below
                expected=math.fsum(factor*r[key] for factor,r in pairs)
                assert abs(row[key]-expected) <= 2e-15*max(row['reference'],1e-300),(path,key,row[key],expected)
            expected=math.fsum(factor*r['closure'] for factor,r in pairs)
            assert abs(row['pairClosure']-expected) <= 1e-15*max(row['reference'],1e-300)
            assert abs(row['closure']-expected) <= 1e-15*max(row['reference'],1e-300)
            assert abs(row['closure']) <= 1e-14*max(row['reference'],1e-300),(path,row)
    assert len(list((case/'postProcessing/phaseChangeElements/0').glob('*.dat'))) == len(set(e for n in names for e in SPECIES[n][1]))
    profiles=list((case/'postProcessing/axialProfiles/0.05').glob('element_*.dat'))
    assert profiles
    for path in profiles:
        assert '# Tpeak ' in path.read_text()
        element=path.stem.removeprefix('element_').removesuffix('_native')
        rows=[list(map(float,l.split())) for l in path.read_text().splitlines() if l and not l.startswith('#')]
        for b,row in enumerate(rows):
            for c in [3,4,5]:
                expected=0
                for name in names:
                    f=case/f'postProcessing/axialProfiles/0.05/{name}_g_native.dat'
                    r=[list(map(float,l.split())) for l in f.read_text().splitlines() if l and not l.startswith('#')][b]
                    expected+=SPECIES[name][1].get(element,0)*r[c]
                assert abs(row[c]-expected) <= 1e-13*max(abs(expected),1e-300)


def independence(multi, singles, names, ranks):
    dirs=[''] if ranks == 1 else [f'processor{i}' for i in range(ranks)]
    for name,single in zip(names,singles):
        for proc in dirs:
            for field in [f'Y_{name}_g',f'Y_{name}_s',f'mw_{name}_s',f'Cs_{name}_s']:
                a=multi/proc/'0.05'/field;b=single/proc/'0.05'/field
                if a.exists(): assert a.read_bytes() == b.read_bytes(),(a,b)
        a=multi/f'postProcessing/phaseChangeBalance/0/{name}_g.dat'
        b=single/f'postProcessing/phaseChangeBalance/0/{name}_g.dat'
        assert a.read_bytes() == b.read_bytes(),(a,b)
        a=multi/f'postProcessing/axialProfiles/0.05/{name}_g_native.dat'
        b=single/f'postProcessing/axialProfiles/0.05/{name}_g_native.dat'
        assert a.read_bytes() == b.read_bytes(),(a,b)


def bounded(case, names, ranks):
    # Recompute the HKS bounded law from fields, rather than solver counters.
    # A tiny BiI3 reservoir exercises CAPPED and NONE alongside PbI2 IMPLICIT.
    table=[list(map(float,l.split(','))) for l in (case/'constant/pv_BiI3.csv').read_text().splitlines() if l and not l.startswith('#')]
    def pressure(T):
        for a,b in zip(table, table[1:]):
            if a[0] <= T <= b[0]:
                w=(1/T-1/a[0])/(1/b[0]-1/a[0])
                return min(10**((1-w)*a[j]+w*b[j]) for j in [1,2])
        raise ValueError(T)
    regimes=set()
    for proc in [''] if ranks == 1 else [f'processor{i}' for i in range(ranks)]:
        root=case/proc
        x=readInternal(root/'0/Cx');y=readInternal(root/'0/Cy');N=len(x)
        selected=[i for i in range(N) if .01 <= x[i] <= .02 and .002 <= y[i] <= .0028]
        T=readInternal(root/'0/T',N);rho=readInternal(root/'0/rho',N)
        for name in names:
            M=SPECIES[name][0]
            sigma,Ce,scale=(.8,.0001,1) if name=='PbI2' else (.6,.0002,.5)
            previous=[(4.84e-8 if name=='PbI2' else 1e-12)*M/8e-9]*N
            for step in range(1,51):
                time=f'{step/1000:g}'
                Y=readInternal(root/time/f'Y_{name}_g',N)
                Cs=readInternal(root/time/f'Cs_{name}_s',N)
                exch=readInternal(root/time/f'exch_{name}_g',N)
                for i in selected:
                    G=scale*Ce*2*sigma/(2-sigma)*math.sqrt(M/(2*math.pi*8.314462618*T[i]))
                    raw=4000*G*(pressure(T[i])-rho[i]*8.314462618*T[i]/M*Y[i])
                    available=previous[i] if step <= 45 else 0
                    upper=available/.001
                    expected=min(max(raw,0),upper)
                    tolerance=3e-12*max(abs(raw),upper,1e-300)
                    assert abs(exch[i]-expected) <= tolerance,(name,step,i,exch[i],expected)
                    assert Cs[i]>=0
                    assert abs(Cs[i]-(available-.001*exch[i])) <= 3e-12*max(available,abs(raw)*.001,1e-300)
                    if step < 45:
                        regimes.add('CAPPED' if raw >= upper and upper>0 else 'NONE' if expected==0 else 'IMPLICIT')
                previous=Cs
    assert {'CAPPED','NONE','IMPLICIT'} <= regimes,regimes


def main():
    task,work,solver,carrier=sys.argv[1:5]
    work=Path(work)
    if task == 'independence':
        for model in ['temperature','HKS']:
            for names in [('PbI2','PbI2copy'),('PbI2','BiI3')]:
                for ranks in [1,4]:
                    prefix=f'{model}_{names[1]}_{ranks}'
                    multi=prepare(work/(prefix+'_multi'),names,model)
                    singles=[prepare(work/(prefix+'_'+n),[n],model) for n in names]
                    for c in [multi]+singles: run(c,solver,carrier,ranks)
                    independence(multi,singles,names,ranks);elements(multi,names)
        print('Both species combinations, temperature/HKS, serial/4 ranks: byte-identical fields, balances, profiles; element closure <=1e-14')
    elif task == 'inventory':
        names=('PbI2','BiI3')
        for ranks in [1,4]:
            multi=prepare(work/f'inventory_{ranks}_multi',names,mode='inventory')
            singles=[prepare(work/f'inventory_{ranks}_{n}',[n],mode='inventory') for n in names]
            for c in [multi]+singles: run(c,solver,carrier,ranks)
            independence(multi,singles,names,ranks);elements(multi,names)
            bounded(multi, names, ranks)
            for n in names:
                for row in balance(multi,n):
                    assert row['sample'] >= 0
                    assert abs(row['closure']) < 1e-13*max(row['initial'],1e-300)
            # restart before removal; final fields and ledgers checked numerically
            restart=work/f'inventory_{ranks}_restart'
            if restart.exists(): shutil.rmtree(restart)
            shutil.copytree(multi,restart)
            for parent in [restart]+([restart/f'processor{i}' for i in range(ranks)] if ranks > 1 else []):
                for p in parent.iterdir():
                    if p.is_dir() and re.fullmatch(r'[0-9.]+',p.name) and float(p.name) > .02: shutil.rmtree(p)
            if ranks == 1: call(restart,[solver],'log.restart')
            else: call(restart,['mpirun','-np','4',solver,'-parallel'],'log.restart')
            for proc in [''] if ranks == 1 else [f'processor{i}' for i in range(4)]:
                for n in names:
                    for field in [f'Y_{n}_g',f'Y_{n}_s',f'Cs_{n}_s']:
                        a=multi/proc/'0.05'/field;b=restart/proc/'0.05'/field
                        assert a.read_bytes() == b.read_bytes(),(a,b)
        print('Shared-cell inventories: independent pairs, bounded reservoirs, closure <=1e-13; serial/MPI restart before removal identical')
    elif task == 'inputs':
        c=prepare(work/'diffusion', ['BiI3'])
        p=c/'constant/speciesTransportProperties';s=p.read_text();s=s.replace('diffusivityModel constant; D 1e-4;','diffusivityModel ChapmanEnskog; ChapmanEnskogCoeffs {carrierMolarMass .004002602; sigma 5; carrierSigma 2.576; epsilonOverK 300; carrierEpsilonOverK 10.22;}');p.write_text(s)
        run(c,solver,carrier)
        assert 'ChapmanEnskog' in (c/'log.solver').read_text()
        tests=[('mass',lambda t:t.replace('molarMass 0.58969381','molarMass .5'),'formula molarMass'),('formula',lambda t:t.replace('I 3;','I -3;'),'positive integer'),('lj',lambda t:t.replace('sigma 5;','sigma 0;'),'finite positive')]
        for label,edit,message in tests:
            p.write_text(edit(s));call(c,[solver],f'log.bad_{label}',True);assert message in (c/f'log.bad_{label}').read_text()
        c=prepare(work/'overlap',['PbI2','BiI3'],mode='inventory')
        call(c,['blockMesh']);call(c,[carrier])
        p=c/'constant/thermochemistryProperties';s=p.read_text();p.write_text(s.replace('samples {boat {', 'samples { other {pair PbI2_g; amount 1e-8; mode inventory; areaPerVolume 1; selection {boat {action use; source box; box (.01 .002 -1) (.02 .0028 1);}}} boat {'))
        call(c,[solver],'log.bad_overlap',True);assert 'same pair' in (c/'log.bad_overlap').read_text()
        c=prepare(work/'restart_formula',['PbI2','BiI3'])
        run(c,solver,carrier)
        p=c/'constant/speciesTransportProperties'
        p.write_text(re.sub(r'formula \{[^}]*\}', '', p.read_text()))
        call(c,[solver],'log.bad_formula_restart',True)
        assert 'formulas changed across restart' in (c/'log.bad_formula_restart').read_text()
        c=prepare(work/'bad_hks',['PbI2'])
        call(c,['blockMesh']);call(c,[carrier])
        p=c/'constant/thermochemistryProperties';text=p.read_text()
        for bad in ['accommodation 1.2;', 'accommodation -.1;']:
            p.write_text(text.replace('accommodation .8;',bad))
            call(c,[solver],'log.bad_hks',True)
            assert '0 < accommodation <= 1' in (c/'log.bad_hks').read_text()
        # Legacy pair/amount syntax still works with explicit formulas.
        p.write_text(text.replace('amounts {PbI2_g 4.84e-08;}', 'pair PbI2_g; amount 4.84e-08;'))
        call(c,[solver],'log.legacy')
        print('Chapman-Enskog solver path; invalid formula/mass/LJ and same-pair overlap refused')
    else: raise SystemExit(task)

if __name__ == '__main__': main()
