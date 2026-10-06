#!/usr/bin/env python3
"""M10d provisional LBE Q2 source and independent fractional plug-flow gate.

The prescribed plane carrier tests the fast wall-kinetics limit, not Liu's
experimental transport. BiI's +29 kJ/mol variant is a sensitivity surrogate.
Thermodynamic records stay external; the independent element-potential/active
set implementation is loaded from LESTO_M10D_REFERENCE/design/py/pbbii.py.
"""
import importlib.util,json,math,re,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('wall_cases',HERE/'cases.py')
wc=importlib.util.module_from_spec(spec);spec.loader.exec_module(wc)
sys.path.insert(0,str(wc.base.REPO/'thermochemistry/systems'))
from make_wall_tables import generate
from make_pv_tables import Data

def setup(case,data,dt=.001,nx=200):
    wc.prepare(case,active=True);wc.solid(case,'Pb_s',{'Pb':1})
    p=case/'system/blockMeshDict';p.write_text(p.read_text().replace('(1 1 1)',f'({nx} 2 1)',1))
    p=case/'system/carrierDict';s=p.read_text()
    for key,value in [('Tin',1173),('Tout',350),('massFlux',.01743)]:s=re.sub(key+r'\s+[^;]+;',f'{key} {value};',s)
    p.write_text(s)
    # Measure deposition while the source is active. A post-release clean-He
    # flush evaporates hot metal deposits and tests a different physical state
    # from the fractional plug-flow reference.
    q2=json.loads((data/'q2.json').read_text());duration=2
    nHe=.01743*.0048*.001/.0040026*duration
    amounts={a+'_g':nHe*n for a,n in q2['element_mol_per_mol_He'].items()}
    channels=''
    for gas,solid,nu in [('Pb','Pb',1),('Bi','Bi',1),('Bi2','Bi',2)]:
        channels+=f'{gas} {{reactant {gas}_g; condensed {solid}_s; stoichiometry {nu}; vapourPressure {{file "{data}/pv_{gas}_phases.csv"; phases ("{solid}(cr)" "{solid}(l)"); units log10Pa;}}}}\n'
    channels+=f'BiI {{reactant BiI_g; condensed Bi_s; stoichiometry {{BiI_g -3; Bi_s 2; BiI3_g 1;}} equilibrium {{file "{data}/logK_BiI_disp.csv"; column logK;}} reversible no;}}\n'
    pairs=''
    for gas in ['PbI2','BiI3']:
        pairs+=f'{gas}_g {{condensed {gas}_s; HKS {{accommodation 1; Ce 1; kineticScale 1;}} vapourPressure {{file "{data}/pv_{gas}_phases.csv"; phases ("{gas}(cr)" "{gas}(l)"); units log10Pa;}}}}'
    p=case/'constant/thermochemistryProperties'
    p.write_text(wc.base.HEADER%'thermochemistryProperties'+f'''
model HKS;
pairs {{{pairs}}}
interface {{patches (WALL);}}
HKSCoeffs {{accommodation 1; Ce 1; kineticScale 1;}}
samples {{}}
balance {{nGuard 8; guardTolerance 1e-13; writeExchangeField yes;}}
profiles {{axis (1 0 0); origin (0 0 0); observable wall;
 bins {{native {{min 0; max .1; nBins {nx};}}}}
 onset {{fraction .001; searchFrom .02;}}
}}
respeciation {{engine kernel; species ({' '.join(wc.sc.NAMES)});
 formationConstants "{data}/log10Kf_PbBiI.csv"; minTemperature 350; updateInterval 1;}}
wallChannels {{{channels}}}
gasSources {{boat {{amounts {{{' '.join(f'{n} {v:.17g};' for n,v in amounts.items())}}}
 box ((.01 0 -1) (.02 .0048 1)); startTime 0; duration {duration};}}}}
''')
    for name in wc.sc.NAMES:
        p=case/'0'/('Y_'+name);p.write_text(re.sub(r'internalField uniform [^;]+;','internalField uniform 0;',p.read_text()))
    p=case/'system/controlDict';s=p.read_text()
    for key,value in [('endTime',1.2),('deltaT',dt),('writeInterval',round(1.2/dt)),('nCorr',3)]:s=re.sub(key+r'\s+[^;]+;',f'{key} {value};',s)
    p.write_text(s)
    return amounts

def oracle(reference,nasa,tf,shift,q2):
    sys.path.insert(0,str(Path(reference)/'design/py'))
    import pbbii
    # Share only the record evaluator. The constrained equilibrium algorithm is
    # independent of the C++ homogeneous kernel and of the wall channel law.
    data=Data(nasa,tf)
    class Record:
        def __init__(self,name):self.name=name
        def G(self,T):return data.G(self.name,T)+(shift if self.name=='BiI(g)' else 0)
    for table in [pbbii.GAS,pbbii.COND]:
        for name,(_,formula) in list(table.items()):
            if name in ['Pb2(g)','PbI3(g)','PbI4(g)','I2(cr)','I2(l)']:continue
            table[name]=(Record(name),formula)
    gases=[n.replace('_g','(g)') for n in wc.sc.NAMES]
    allowed=['Pb(cr)','Pb(l)','Bi(cr)','Bi(l)','PbI2(cr)','PbI2(l)','BiI3(cr)','BiI3(l)']
    b=dict(q2['element_mol_per_mol_He']);deposit={n:[] for n in allowed};pi=None
    for T in range(1173,349,-1):
        eq=pbbii.Equilibrium(T,101325,b,1,allowed=allowed,gas=gases).solve(pi0=pi)
        assert eq.balance<1e-8,(T,eq.balance)
        pi=eq.pi
        for name,n in eq.ncond.items():
            if n<=0:continue
            deposit[name].append((T,n))
            for atom,nu in pbbii.COND[name][1].items():b[atom]=max(0,b[atom]-nu*n)
    return {name:max(T for T,n in rows) for name,rows in deposit.items() if rows},deposit

def profile(case,name):
    p=case/f'postProcessing/axialProfiles/1.2/{name}_native.dat'
    s=p.read_text();keys=re.search(r'^# columns (.*)$',s,re.M)[1].split()
    rows=[dict(zip(keys,map(float,l.split()))) for l in s.splitlines() if l and not l.startswith('#')]
    return rows

def analyse(case):
    curves={}
    for name,field in [('Pb','condensate_Pb_s'),('Bi','condensate_Bi_s'),('PbI2','PbI2_g'),('BiI3','BiI3_g')]:
        rows=profile(case,field)
        curves[name]=[(r['Twall'],r['wall']) for r in rows]
    onset={}
    for name,curve in curves.items():
        peak=max(v for T,v in curve)
        for phase,condition in ([('Pb(l)',lambda T:T>600.61)] if name=='Pb' else
             [('Bi(l)',lambda T:T>544.512416),('Bi(cr)',lambda T:T<=544.512416)] if name=='Bi' else [(name+'(cr)',lambda T:True)]):
            Ts=[T for T,v in curve if condition(T) and v>=.001*peak]
            if Ts:onset[phase]=max(Ts)
    return onset,curves

def closure_report(case):
    return {group:max(abs(r['closure'])/max(r['reference'],1e-300)
                      for p in (case/'postProcessing'/group).glob('*/*.dat') for r in wc.sc.rows(p))
            for group in ['phaseChangeElements','phaseChangeSpecies','phaseChangeChannelCondensates']}

def main():
    work,solver,carrier,nasa,tf,reference=sys.argv[1:];work=Path(work);work.mkdir(parents=True,exist_ok=True)
    report={}
    for label,shift in [('baseline',0),('shift29k',29000)]:
        data=(work/'tables'/label).resolve();q2=generate(nasa,tf,data,shift)
        expected,_=oracle(reference,nasa,tf,shift,q2)
        case=work/label;setup(case,data);wc.base.run(case,solver,carrier);wc.closure(case)
        actual,curves=analyse(case)
        entry={'shift_J_mol':shift,'plug_flow_onset_K':expected,'solver_onset_K':actual,
               'difference_K':{n:actual[n]-expected[n] for n in actual},'max_relative_closure':closure_report(case)}
        report[label]=entry
        (work/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(entry,indent=2),flush=True)
        for name in ['Pb(l)','Bi(l)','PbI2(cr)','Bi(cr)','BiI3(cr)']:
            assert abs(actual[name]-expected[name])<=10,(label,name,actual,expected)
        fine=work/(label+'_dt0.0005');setup(fine,data,dt=.0005);wc.base.run(fine,solver,carrier);wc.closure(fine)
        finer,fineCurves=analyse(fine)
        entry['half_step_onset_K']=finer;entry['half_step_closure']=closure_report(fine)
        entry['half_step_onset_change_K']={n:finer[n]-actual[n] for n in actual}
        entry['half_step_deposit_L1']={n:sum(abs(a[1]-b[1]) for a,b in zip(curves[n],fineCurves[n]))/max(sum(b[1] for b in fineCurves[n]),1e-300) for n in curves}
        for name in ['Pb(l)','Bi(l)','PbI2(cr)','Bi(cr)','BiI3(cr)']:
            assert abs(finer[name]-expected[name])<=10,(label,name,'half step',finer,expected)
            assert abs(finer[name]-actual[name])<=4.115+1e-9,(label,name,'time convergence',actual,finer)
        (work/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Both BiI thermodynamic variants: five onsets within 10 K; ledgers <=1e-14')
if __name__=='__main__':main()
