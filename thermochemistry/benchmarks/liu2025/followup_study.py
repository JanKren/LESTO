#!/usr/bin/env python3
"""Prepare and execute isolated refinement and explicitly provisional sensitivities.

The completed matrix remains immutable. New jobs use its private solver and
tables. No experimental temperature correction or physical coefficient is fitted.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading

from carriers import REPO, prepare_carrier, call, field, header
from cases import prepare_case
from report import analyse

FULL = REPO/'run/034-liu2025/work/prepared'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
SENSITIVITY = BASES[:3]+['LBE-I_SiO2_I_Q2_baseline_sat0.001']


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def change_control(case, dt=.002, end=60):
    path=case/'system/controlDict';text=path.read_text()
    for name,value in [('deltaT',format(dt,'.17g')),('endTime',format(end,'.17g')),
                       ('writeInterval',str(round(1/dt))),('purgeWrite','2')]:
        text,count=re.subn(r'\b'+name+r'\s+[^;]+;',name+' '+value+';',text)
        assert count==1,(path,name,count)
    path.write_text(text)
    manifest=json.loads((case/'benchmark.json').read_text())
    manifest.update(dt_s=dt,end_time_s=end)
    manifest['execution']={'field_write_period_s':1,'purge_write':2}
    write_json(case/'benchmark.json',manifest)


def clone_case(base, destination):
    assert not destination.exists(),destination
    destination.mkdir(parents=True)
    for folder in ['0','constant','system']:
        shutil.copytree(base/folder,destination/folder)
    shutil.copyfile(base/'benchmark.json',destination/'benchmark.json')
    change_control(destination)


def geometry_text(text, x0, x1):
    # Preserve radial coordinates and source amounts, changing only the boat box.
    # Match numerical coordinates from the actual dictionary rather than spelling.
    pattern=r'(box\s+\(?\s*\()\s*([^\s()]+)\s+([^()]+)\)\s*\(\s*([^\s()]+)\s+([^()]+)\)'
    def replacement(match):
        assert abs(float(match[2])-.09)<1e-10 and abs(float(match[4])-.105)<1e-10
        return f'{match[1]}{x0:.17g} {match[3]}) ({x1:.17g} {match[5]})'
    changed,count=re.subn(pattern,replacement,text)
    assert count==1,count
    return changed


def append_shutdown(out, entries):
    for base in SENSITIVITY:
        name=base+'__flowStop60to80'
        if any(entry['case']==name for entry in entries):continue
        source=FULL/'cases'/base;case=out/'cases'/name;clone_case(source,case)
        changes={'start_time_s':60,'end_time_s':80,'release_end_time_s':60,'carrier_flow':'zero',
                 'ends':'zero normal species flux','thermal_scope':'Frozen measured temperature field; isothermal relaxation diagnostic, not experimental cooling.'}
        manifest=json.loads((case/'benchmark.json').read_text())
        manifest['followup']={'base_case':base,'category':'shutdown','changes':changes,
                              'scientific_scope':'Isothermal relaxation only; experiment cooling history unavailable.'}
        write_json(case/'benchmark.json',manifest)
        shutil.copytree(source/'60',case/'60');change_control(case,end=80)
        field(case,'U','0 1 -1 0 0 0 0','(0 0 0)',
            'type fixedValue; value uniform (0 0 0);','type fixedValue; value uniform (0 0 0);','type noSlip;',vector=True)
        phi=header('phi','surfaceScalarField')+'dimensions [1 0 -1 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n'
        for patch in ['INLET','OUTLET','WALL']:phi+=patch+' {type calculated; value uniform 0;}\n'
        phi+='BACK {type wedge; value uniform 0;} FRONT {type wedge; value uniform 0;} AXIS {type empty; value uniform 0;}\n}\n'
        (case/'0/phi').write_text(phi)
        for path in (case/'60').glob('Y_*'):
            text,count=re.subn(rb'\bINLET\s*\{[^}]*\}',b'INLET {type zeroGradient;}',path.read_bytes())
            assert count==1,(path,count);path.write_bytes(text)
        entries.append({'case':name,'base_case':base,'category':'shutdown','changes':changes,
                        'status':'PREPARED_PROVISIONAL','experiment':manifest['experiment'],
                        'source_mode':manifest['source_mode'],'data_variant':Path(manifest['table_directory']).name,
                        'saturation':manifest['metal_saturation_fraction'],'reason':''})
        print('Prepared',name,flush=True)


def finalise(args):
    out=args.out.resolve();assert not (out/'run_state.json').exists()
    plan=json.loads((out/'study.json').read_text());entries=plan['cases']
    append_shutdown(out,entries)
    for entry in entries:
        case=out/'cases'/entry['case'];metadata=json.loads((case/'benchmark.json').read_text())
        if entry['category']=='kinetics':
            alpha=entry['changes'].get('mass_accommodation',entry['changes'].get('iodide_mass_accommodation'))
            assert alpha is not None
            entry['changes']={'iodide_mass_accommodation':alpha,'metal_channel_mass_accommodation':1,
                              'scope':'PbI2/BiI3 pair accommodation changed; wall-channel defaults retained.'}
            metadata['followup']['changes']=entry['changes'];write_json(case/'benchmark.json',metadata)
        if entry['changes'].get('PbI2_diffusivity_model')=='PbI2He':
            metadata['gas_diffusivity_m2_s']=None
            metadata['species_diffusivity_models']={'PbI2_g':'PbI2He(T,p)','other_gases':'constant 1e-4 m2/s'}
            write_json(case/'benchmark.json',metadata)
        base=json.loads((FULL/'cases'/entry['base_case']/'benchmark.json').read_text())
        if entry['category'] not in ['source_loading','flow_reference']:
            assert metadata['source_wedge_mol_s']==base['source_wedge_mol_s'],entry['case']
        if entry['category']=='numerical':
            assert metadata['table_sha256']==base['table_sha256']
            assert metadata['gas_diffusivity_m2_s']==base['gas_diffusivity_m2_s']
            assert metadata['end_time_s']==base['end_time_s']==60
        entry['manifest_sha256']=digest(case/'benchmark.json')
        entry['input_sha256']={str(p.relative_to(case)):digest(p) for directory in ['0','constant','system']
                               for p in sorted((case/directory).rglob('*')) if p.is_file()}
    plan.update(cases=entries,status='PREPARED_AND_AUDITED',finalised_utc=utc(),
                finalisation_script_sha256=digest(Path(__file__)))
    write_json(out/'study.json',plan);write_json(out/'cases/matrix.json',entries)
    print('Audited',len(entries),'configurations; numerical sources, tables and end times unchanged.',flush=True)


def prepare(args):
    out=args.out.resolve();assert not out.exists(),'Use a fresh study directory'
    out.mkdir(parents=True);(out/'cases').mkdir();(out/'bin').mkdir()
    for name,source in [('rhoFixedFlowFoam',FULL/'bin/rhoFixedFlowFoam'),
                        ('liuCarrierTool',args.carrier_tool)]:
        shutil.copy2(source,out/'bin'/name)
    assert digest(out/'bin/rhoFixedFlowFoam')==json.loads((FULL/'run_state.json').read_text())['solver_sha256']
    sys.path.insert(0,str(REPO/'thermochemistry/systems'))
    from make_pv_tables import Data
    data=Data(str(args.nasa),str(args.tf))
    entries=[]
    def add(base,label,category,changes,carrier=None,saturation=None):
        source=FULL/'cases'/base;metadata=json.loads((source/'benchmark.json').read_text())
        name=base+'__'+label;case=out/'cases'/name
        if carrier:
            prepare_case(case,carrier,Path(metadata['table_directory']),metadata['source_mode'],
                         saturation if saturation is not None else metadata['metal_saturation_fraction'],
                         .002,60,1e-4,data)
            change_control(case)
        else:clone_case(source,case)
        if saturation is not None and not carrier:
            raise ValueError('A new source saturation needs an audited carrier')
        manifest=json.loads((case/'benchmark.json').read_text())
        manifest['followup']={'base_case':base,'category':category,'changes':changes,
                              'scientific_scope':'Sensitivity only; no experimental fit or dataset qualification.' if category!='numerical' else 'Isolated numerical refinement.'}
        write_json(case/'benchmark.json',manifest)
        entries.append({'case':name,'base_case':base,'category':category,'changes':changes,
                        'status':'PREPARED_PROVISIONAL','experiment':manifest['experiment'],
                        'source_mode':manifest['source_mode'],'data_variant':Path(manifest['table_directory']).name,
                        'saturation':manifest['metal_saturation_fraction'],'reason':''})
        print('Prepared',name,flush=True)
        return case
    def carrier(name,experiment,nx=384,nr=48,flow_T=298.15,temperature_root=None):
        directory=out/'carriers'/name/experiment
        if not directory.exists():
            prepare_carrier(directory,experiment,nx,nr,flow_T,101325,temperature_root)
            for command,label in [(['blockMesh'],'blockMesh'),(['checkMesh','-constant'],'checkMesh'),
                    ([str(out/'bin/liuCarrierTool')],'initialise'),(['rhoSimpleFoam'],'carrier')]:
                call(directory,command,label)
            times=sorted((p for p in directory.iterdir() if p.is_dir() and p.name.replace('.','',1).isdigit()),key=lambda p:float(p.name))
            for command,label in [(['-project'],'project'),(['-audit'],'audit')]:
                call(directory,[str(out/'bin/liuCarrierTool')]+command+['-time',times[-1].name],label)
            assert json.loads((directory/'carrier-audit.json').read_text())['status']=='PASS'
        return directory
    # Four representative full-duration time checks and a uniform-ratio grid sequence.
    for base in BASES:
        case=add(base,'dt1ms','numerical',{'dt_s':.001,'cells':18432})
        change_control(case,.001)
    for nx,nr in [(256,32),(576,72)]:
        for base in BASES:
            experiment=json.loads((FULL/'cases'/base/'benchmark.json').read_text())['experiment']
            add(base,f'grid{nx}x{nr}','numerical',{'nx':nx,'nr':nr,'cells':nx*nr,'dt_s':.002},
                carrier(f'grid{nx}x{nr}',experiment,nx,nr))
    for base in SENSITIVITY:
        for diffusivity in [5e-5,2e-4]:
            case=add(base,f'D{diffusivity:g}','transport',{'constant_D_m2_s':diffusivity})
            path=case/'constant/speciesTransportProperties'
            text,count=re.subn(r'\bD\s+0\.0001;',f'D {diffusivity:.17g};',path.read_text())
            # Original D is written in full precision by prepare_case.
            if not count:
                text,count=re.subn(r'\bD\s+[^;]+;',f'D {diffusivity:.17g};',path.read_text())
            assert count==9,count;path.write_text(text)
            manifest=json.loads((case/'benchmark.json').read_text());manifest['gas_diffusivity_m2_s']=diffusivity
            write_json(case/'benchmark.json',manifest)
        for alpha in [.1,.01]:
            case=add(base,f'alpha{alpha:g}','kinetics',{'mass_accommodation':alpha})
            path=case/'constant/thermochemistryProperties'
            text,count=re.subn(r'\baccommodation\s+1;',f'accommodation {alpha:.17g};',path.read_text())
            assert count==3,count;path.write_text(text)
        for dx in [-.015,.015]:
            case=add(base,'boat'+('upstream' if dx<0 else 'downstream'),'source_geometry',{'source_axial_shift_m':dx})
            path=case/'constant/thermochemistryProperties'
            path.write_text(geometry_text(path.read_text(),.09+dx,.105+dx))
            manifest=json.loads((case/'benchmark.json').read_text())
            manifest['source_box_m'][0][0]+=(dx);manifest['source_box_m'][1][0]+=(dx)
            write_json(case/'benchmark.json',manifest)
        case=add(base,'PbI2He','transport',{'PbI2_diffusivity_model':'PbI2He','other_species_D_m2_s':1e-4})
        path=case/'constant/speciesTransportProperties'
        text,count=re.subn(r'(\bPbI2_g\s*\{[^}]*?)\bdiffusivityModel\s+constant;\s*D\s+[^;]+;',r'\1diffusivityModel PbI2He;',path.read_text())
        assert count==1,count
        path.write_text(text+'\npressureUnit atm;\n')
    # Refine source loading between existing 0.1% and 1% branches, without fitting.
    for experiment in ['LBE-I_SS_II','LBE-I_SiO2_I']:
        original_carrier=Path(args.full_carriers)/experiment
        for variant in ['baseline','shift29k']:
            for saturation in [.002,.0031622776601683794,.005623413251903491]:
                base=experiment+'_Q2_'+variant+'_sat0.01'
                add(base,f'sat{saturation:.9g}','source_loading',{'metal_saturation_fraction':saturation},original_carrier,saturation)
    # Common reference conventions are alternatives, not a calibration correction.
    for base in SENSITIVITY:
        experiment=json.loads((FULL/'cases'/base/'benchmark.json').read_text())['experiment']
        add(base,'flow0C','flow_reference',{'flow_reference_temperature_K':273.15,'flow_reference_pressure_Pa':101325},
            carrier('flow0C',experiment,flow_T=273.15))
    # Species-derived alignment hypotheses bracket the internal Figure 7C conflict.
    diagnostics=json.loads((REPO/'run/034-liu2025/full_analysis/analysis.json').read_text())['silica_figure_diagnostics']
    for diagnostic in diagnostics:
        shift=diagnostic['diagnostic_offset_cm'][0]/100
        temperature_root=out/'temperature_hypotheses'/diagnostic['species'];temperature_root.mkdir(parents=True)
        import csv
        source=REPO/'thermochemistry/benchmarks/liu2025/digitized/LBE-I_SiO2_I_temperature.csv'
        with source.open() as stream:points=list(csv.DictReader(stream))
        with (temperature_root/source.name).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=['x_m','T_K'],lineterminator='\n');writer.writeheader()
            writer.writerows({'x_m':float(p['x_m'])-shift,'T_K':p['T_K']} for p in points if 0<=float(p['x_m'])-shift<=1.25)
        new_carrier=carrier('silica_'+diagnostic['species'],'LBE-I_SiO2_I',temperature_root=temperature_root)
        for base in ['LBE-I_SiO2_I_Q1_baseline','LBE-I_SiO2_I_Q2_baseline_sat0.001']:
            add(base,'curveShift'+diagnostic['species'],'silica_alignment',{'temperature_curve_shift_m':-shift,
                'scope':'Explicit hypothesis derived from curve/label disagreement; not a corrected measurement.'},new_carrier)
    # Closed isothermal relaxation preserves rho/T and the restart inventory.
    append_shutdown(out,entries)
    write_json(out/'cases/matrix.json',entries)
    write_json(out/'study.json',{'status':'PREPARED','created_utc':utc(),'cases':entries,
        'solver_sha256':digest(out/'bin/rhoFixedFlowFoam'),'solver_source_commit':'fe214c6',
        'base_root':str(FULL/'cases'),'source_files_sha256':{p.name:digest(p) for p in [Path(__file__),REPO/'thermochemistry/benchmarks/liu2025/cases.py',REPO/'thermochemistry/benchmarks/liu2025/carriers.py']},
        'scope':'60 s cold-start numerical studies and physical sensitivity scenarios. Cooling history remains unavailable.'})
    print('Study prepared:',len(entries),'cases',flush=True)


def run(args):
    out=args.out.resolve();plan=json.loads((out/'study.json').read_text())
    assert not (out/'run_state.json').exists(),'Existing execution state; use explicit recovery rather than overwriting'
    assert digest(out/'bin/rhoFixedFlowFoam')==plan['solver_sha256']
    assert plan['status']=='PREPARED_AND_AUDITED'
    smoke_record=json.loads((out/'smoke.json').read_text())
    assert smoke_record['status']=='PASS'
    assert smoke_record['study_manifest_sha256']==digest(out/'study.json')
    assert len(smoke_record['cases'])==len(plan['cases'])
    for entry in plan['cases']:
        case=out/'cases'/entry['case']
        assert digest(case/'benchmark.json')==entry['manifest_sha256'],entry['case']
        for path,expected in entry['input_sha256'].items():assert digest(case/path)==expected,(entry['case'],path)
    state={'status':'RUNNING','started_utc':utc(),'pid':os.getpid(),'jobs':args.jobs,
           'solver_sha256':plan['solver_sha256'],'cases':{r['case']:{'status':'QUEUED','category':r['category']} for r in plan['cases']}}
    lock=threading.RLock();stop=threading.Event()
    def save():
        with lock:
            state['updated_utc']=utc()
            for name,record in state['cases'].items():
                log=out/'cases'/name/'log.solver'
                if record['status']=='RUNNING' and log.exists():
                    with log.open('rb') as stream:
                        stream.seek(max(0,log.stat().st_size-65536));tail=stream.read().decode(errors='replace')
                    times=re.findall(r'^Time = ([0-9.eE+\-]+)$',tail,re.M)
                    if times:record['latest_time_s']=float(times[-1])
            write_json(out/'run_state.json',state)
    def job(entry):
        name=entry['case'];case=out/'cases'/name
        assert not (case/'log.solver').exists(),case
        with (case/'log.solver').open('w') as log:
            process=subprocess.Popen([str(out/'bin/rhoFixedFlowFoam'),'-case',str(case)],stdout=log,stderr=subprocess.STDOUT)
            with lock:state['cases'][name].update(status='RUNNING',pid=process.pid,started_utc=utc())
            save();code=process.wait()
        assert code==0,(name,code)
        with (case/'log.solver').open('rb') as log:
            log.seek(max(0,log.seek(0,2)-4096));assert log.read().rstrip().endswith(b'End'),name
        result,_=analyse(case)
        with lock:state['cases'][name].update(status=result['status'],finished_utc=utc())
        save();print(name,result['status'],flush=True)
    def monitor():
        while not stop.wait(30):save()
    save();threading.Thread(target=monitor,daemon=True).start()
    try:
        if args.report_out:
            initial_report=subprocess.run([sys.executable,str(REPO/'thermochemistry/benchmarks/liu2025/report_followup.py'),
                                           '--root',str(out),'--out',str(args.report_out)])
            assert initial_report.returncode==0,'Initial pending report failed'
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending={pool.submit(job,entry):entry['case'] for entry in plan['cases']}
            for future in concurrent.futures.as_completed(pending):
                try:future.result()
                except Exception as error:
                    with lock:state['cases'][pending[future]].update(status='FAIL',reason=str(error),finished_utc=utc())
                    save();print('FAIL',pending[future],error,flush=True)
        result=subprocess.run([sys.executable,str(REPO/'thermochemistry/benchmarks/liu2025/report.py'),str(out/'cases')])
        state.update(status='COMPLETE' if result.returncode==0 and all(r['status']=='PASS' for r in state['cases'].values()) else 'FAILED',finished_utc=utc())
        save()
        analysis=subprocess.run([sys.executable,str(REPO/'thermochemistry/benchmarks/liu2025/report_followup.py'),
                                 '--root',str(out),'--out',str(args.report_out or out/'report')])
        state['analysis_returncode']=analysis.returncode
        if analysis.returncode:state['status']='FAILED_ANALYSIS'
    finally:stop.set();save()


def smoke(args):
    out=args.out.resolve();plan=json.loads((out/'study.json').read_text());root=out/'smoke'
    if root.exists():assert args.reuse_completed,'Use a fresh smoke directory or explicit reuse'
    else:root.mkdir()
    def job(entry):
        source=out/'cases'/entry['case'];case=root/entry['case']
        if case.exists() and args.reuse_completed:
            for directory in ['0','constant']:
                for path in (source/directory).rglob('*'):
                    if path.is_file():assert digest(path)==digest(case/path.relative_to(source)),(entry['case'],path)
            result,_=analyse(case)
            with (case/'log.solver').open('rb') as stream:
                stream.seek(max(0,stream.seek(0,2)-4096));assert stream.read().rstrip().endswith(b'End'),case
            assert result['status']=='PASS',case
            metadata=json.loads((case/'benchmark.json').read_text())
            return {'case':entry['case'],'status':'PASS','steps':10,
                    'start_time_s':60 if entry['category']=='shutdown' else 0,'end_time_s':metadata['end_time_s'],
                    'element_gates':result['element_gates'],'max_respeciation_failures':result['max_respeciation_failures']}
        clone_case(source,case)
        metadata=json.loads((source/'benchmark.json').read_text());dt=metadata['dt_s']
        start=60 if entry['category']=='shutdown' else 0
        if start:shutil.copytree(source/'60',case/'60')
        change_control(case,dt,start+10*dt)
        control=case/'system/controlDict'
        control.write_text(re.sub(r'\bwriteInterval\s+[^;]+;','writeInterval 1;',control.read_text()))
        call(case,[str(out/'bin/rhoFixedFlowFoam')],'solver')
        result,_=analyse(case)
        with (case/'log.solver').open('rb') as stream:
            stream.seek(max(0,stream.seek(0,2)-4096));assert stream.read().rstrip().endswith(b'End'),case
        assert result['status']=='PASS',(case,result)
        print('Smoke PASS',entry['case'],flush=True)
        return {'case':entry['case'],'status':'PASS','steps':10,'start_time_s':start,'end_time_s':start+10*dt,
                'element_gates':result['element_gates'],'max_respeciation_failures':result['max_respeciation_failures']}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:results=list(pool.map(job,plan['cases']))
    write_json(out/'smoke.json',{'status':'PASS','cases':results,'solver_sha256':plan['solver_sha256'],
                               'study_manifest_sha256':digest(out/'study.json'),
                               'scope':'Ten native steps per prepared configuration; physical inputs and source rates retained.'})


def launch(args):
    out=args.out.resolve();assert out.exists()
    if (out/'pipeline.json').exists():
        assert (args.retry_failed_launch or args.resume_after_smoke) and not (out/'run_state.json').exists(),'Existing execution; inspect before recovery'
        previous=json.loads((out/'pipeline.json').read_text())
        process_file=Path('/proc')/str(previous['pid'])/'cmdline'
        assert not process_file.exists() or b'launch_followup.sh' not in process_file.read_bytes(),'Previous pipeline is still active'
        if args.resume_after_smoke:
            smoke_record=json.loads((out/'smoke.json').read_text())
            assert smoke_record['status']=='PASS' and smoke_record['study_manifest_sha256']==digest(out/'study.json')
            assert not any((out/'cases'/entry['case']/'log.solver').exists() for entry in json.loads((out/'study.json').read_text())['cases'])
            label='smoke'
        else:
            assert 'unbound variable' in (out/'log.pipeline').read_text(),'Recovery only supports failed environment initialisation'
            label='environment'
        assert not (out/f'pipeline-failed-{label}.json').exists()
        (out/'pipeline.json').rename(out/f'pipeline-failed-{label}.json')
        (out/'log.pipeline').rename(out/f'log.pipeline-failed-{label}')
    with (out/'log.pipeline').open('w') as log:
        process=subprocess.Popen(['bash',str(Path(__file__).with_name('launch_followup.sh')),
            str(out),str(args.jobs),str(args.report_out.resolve())],stdin=subprocess.DEVNULL,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=REPO)
    write_json(out/'pipeline.json',{'status':'LAUNCHED','pid':process.pid,'launched_utc':utc(),
        'jobs':args.jobs,'report_out':str(args.report_out.resolve()),
        'scope':'Wait for preparation; audit all inputs; smoke every case; run full study; generate report. No automatic git writes.'})
    print('Launched durable study pipeline PID',process.pid,flush=True)


def status(args):
    out=args.out.resolve()
    payload={'prepared_case_directories':sum((case/'benchmark.json').exists() for case in (out/'cases').iterdir()),
             'audited_carriers':len(list((out/'carriers').glob('*/*/carrier-audit.json'))),
             'pipeline':json.loads((out/'pipeline.json').read_text()) if (out/'pipeline.json').exists() else None}
    if (out/'run_state.json').exists():
        state=json.loads((out/'run_state.json').read_text());counts={}
        for record in state['cases'].values():counts[record['status']]=counts.get(record['status'],0)+1
        payload.update(status=state['status'],updated_utc=state['updated_utc'],case_status_counts=counts,
            active_cases={name:record.get('latest_time_s') for name,record in state['cases'].items() if record['status']=='RUNNING'})
    else:payload['status']='PREPARING_OR_SMOKE_PENDING'
    print(json.dumps(payload,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare');p.add_argument('--out',type=Path,required=True)
    p.add_argument('--nasa',type=Path,required=True);p.add_argument('--tf',type=Path,required=True)
    p.add_argument('--carrier-tool',type=Path,required=True);p.add_argument('--full-carriers',type=Path,required=True)
    p=sub.add_parser('run');p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,default=12)
    p.add_argument('--report-out',type=Path,help='Generate a report automatically after the jobs finish')
    p=sub.add_parser('smoke');p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--reuse-completed',action='store_true')
    p=sub.add_parser('finalise');p.add_argument('--out',type=Path,required=True)
    p=sub.add_parser('launch');p.add_argument('--out',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=16);p.add_argument('--report-out',type=Path,required=True)
    p.add_argument('--retry-failed-launch',action='store_true')
    p.add_argument('--resume-after-smoke',action='store_true')
    p=sub.add_parser('status');p.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='prepare':prepare(args)
    elif args.command=='finalise':finalise(args)
    elif args.command=='launch':
        assert args.jobs>0
        launch(args)
    elif args.command=='status':status(args)
    elif args.command=='smoke':
        assert args.jobs>0
        smoke(args)
    else:
        assert args.jobs>0
        run(args)


if __name__=='__main__':main()
