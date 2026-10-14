#!/usr/bin/env python3
"""Isolated, fingerprinted follow-up runs using original workstation inputs."""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import tarfile

PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE/'payload/thermochemistry/benchmarks/liu2025'))
import numpy as np
from carriers import prepare_carrier, call
from report import analyse, rows
from analyse_matrix import ledger_audit, log_audit

EXPECTED = '16fd804e593f39ed7d25f42ea3459d9d8de196c3456b7d3f5702c54a83219a58'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
LOCK = threading.RLock()


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text())
def save(p, value):
    p=Path(p); temporary=p.with_suffix(p.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n'); temporary.replace(p)
def require(ok, message):
    if not ok: raise RuntimeError(message)
def tree(p): return {str(f.relative_to(p)):sha(f) for f in sorted(p.rglob('*')) if f.is_file()}
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()


def solver_call(case, solver):
    with (case/'log.solver').open('w') as log:
        process=subprocess.Popen([str(solver),'-case',str(case)],stdout=log,stderr=subprocess.STDOUT)
        save(case/'solver-process.json',{'pid':process.pid,'case':str(case),'started_utc':utc()})
        code=process.wait()
    require(code==0,f'Solver exited with code {code}: {case}')


def refuse_live_solver(case):
    record=case/'solver-process.json'
    if record.exists():
        process=read(record); proc=Path('/proc')/str(process['pid'])/'cmdline'
        if proc.exists():
            command=proc.read_bytes().split(b'\0')
            require(str(case).encode() not in command,'A solver is still running in '+str(case))


def control(text, correctors=3, end=60., smoke=False):
    for key,value in [('startFrom','startTime'),('startTime','0'),('endTime',str(end)),
                      ('deltaT','0.002'),('writeInterval','1' if smoke else '500'),
                      ('purgeWrite','2'),('nCorr',str(correctors))]:
        text,n=re.subn(r'\b'+key+r'\s+[^;]+;',key+' '+value+';',text)
        require(n==1,f'Expected exactly one {key} entry, got {n}')
    require(bool(re.search(r'\btolerance\s+1e-16\s*;',text)), 'Original outer tolerance not found')
    return text


def locations(args):
    prepared=(args.prepared or args.repo/'run/024-liu2025/work/prepared').resolve()
    solver=(args.solver or prepared/'bin/rhoFixedFlowFoam').resolve()
    tool=args.carrier_tool
    if tool is None:
        found=shutil.which('liuCarrierTool')
        tool=Path(found) if found else prepared/'bin/liuCarrierTool'
    return prepared,solver,tool.resolve()


def check(args):
    # Pin the packaged helpers and reference files before executing native tools.
    for name,digest in read(PACKAGE/'package-manifest.json').items():
        require(sha(PACKAGE/name)==digest,'Package file changed: '+name)
    require(os.environ.get('WM_PROJECT_VERSION')=='v2412','Load the original OpenFOAM v2412 environment first')
    for command in ['blockMesh','checkMesh','rhoSimpleFoam']:
        require(shutil.which(command),command+' not found in PATH')
    prepared,solver,tool=locations(args)
    require(solver.is_file() and sha(solver)==EXPECTED,'Original solver missing or fingerprint mismatch: '+str(solver))
    for binary in [solver,tool]:
        require(binary.is_file() and os.access(binary,os.X_OK),'Executable missing: '+str(binary))
    for base in BASES:
        case=prepared/'cases'/base
        meta=read(case/'benchmark.json')
        require(meta['dt_s']==.002 and meta['end_time_s']==60,'Baseline resolution differs: '+base)
        require(meta['carrier']['cells']==18432,'Baseline mesh differs: '+base)
        for folder in ['0','constant','system']:
            require((case/folder).is_dir(),'Missing input: '+str(case/folder))
        require((case/'constant/polyMesh/points').is_file(),'Missing baseline mesh: '+base)
        control((case/'system/controlDict').read_text())
        require(re.search(r'\bnCorr\s+3\s*;',(case/'system/controlDict').read_text()),'Original baseline does not use three correctors')
        table=Path(meta['table_directory'])
        for name,digest in meta['table_sha256'].items():
            require((table/name).is_file() and sha(table/name)==digest,'Table fingerprint mismatch: '+str(table/name))
        expected=read(PACKAGE/'references'/base/'benchmark.json')
        for key in ['table_sha256','source_wedge_mol_s','gas_diffusivity_m2_s','source_box_m']:
            require(meta[key]==expected[key],f'Baseline {key} differs: {base}')
        for species in ['element_I','PbI2_g','BiI3_g']:
            path=case/'postProcessing/axialProfiles/60'/f'{species}_native.dat'
            require(path.is_file() and sha(path)==sha(PACKAGE/'references'/base/path.name),
                    'Baseline final profile differs: '+str(path))
        # Remeshing only supports initially empty, uniform trace inventories.
        for path in (case/'0').glob('Y_*'):
            require(re.search(rb'internalField\s+uniform\s+0\s*;',path.read_bytes()),
                    'Cannot remesh nonuniform/nonzero initial species: '+str(path))
    print('Prerequisites PASS: original solver, inputs, tables and reference profiles.',flush=True)
    return prepared,solver,tool


def prepare(args):
    prepared,solver,tool=check(args)
    out=args.out.resolve()
    require(not out.exists(),'Output exists; choose a fresh --out or use --resume after an interrupted run')
    out.mkdir(parents=True); (out/'bin').mkdir(); (out/'inputs').mkdir(); (out/'cases').mkdir()
    for source in [solver,tool]: shutil.copy2(source,out/'bin'/source.name)
    # Original tables are used on this workstation; never added to the return archive.
    entries=[]; carriers={}
    for base in BASES:
        source=prepared/'cases'/base; original=read(source/'benchmark.json')
        experiment=original['experiment']
        if experiment not in carriers:
            carrier=out/'carriers'/experiment
            cm=prepare_carrier(carrier,experiment,864,108,
                original['carrier']['flow_reference_temperature_K'],original['carrier']['flow_reference_pressure_Pa'])
            for key in ['temperature_csv_sha256','radius_m','wedge_mass_flow_kg_s','wedge_fraction']:
                require(cm[key]==original['carrier'][key],f'Carrier physical input changed: {key}')
            require(sha(carrier/'constant/thermophysicalProperties')==sha(source/'constant/thermophysicalProperties'),
                    'Carrier thermophysical properties differ from original baseline')
            for cmd,label in [(['blockMesh'],'blockMesh'),(['checkMesh','-constant'],'checkMesh'),
                              ([str(out/'bin'/tool.name)],'initialise'),(['rhoSimpleFoam'],'carrier')]:
                call(carrier,cmd,label)
            times=[p for p in carrier.iterdir() if p.is_dir() and re.fullmatch(r'[0-9]+(?:\.[0-9]+)?',p.name)]
            latest=max(times,key=lambda p:float(p.name))
            for flag in ['-project','-audit']:
                call(carrier,[str(out/'bin'/tool.name),flag,'-time',latest.name],flag[1:])
            require(read(carrier/'carrier-audit.json')['status']=='PASS','Carrier audit failed')
            carriers[experiment]=(carrier,latest)
        for kind,ncorr in [('grid864x108',3),('correctors8_full60s',8)]:
            name=base+'__'+kind; dest=out/'inputs'/name; dest.mkdir()
            for folder in ['0','constant','system']:
                ignore=(lambda path,names:['polyMesh'] if 'polyMesh' in names else []) if kind=='grid864x108' and folder=='constant' else None
                shutil.copytree(source/folder,dest/folder,ignore=ignore)
            meta=json.loads(json.dumps(original))
            if kind=='grid864x108':
                carrier,latest=carriers[experiment]
                shutil.copytree(carrier/'constant/polyMesh',dest/'constant/polyMesh')
                shutil.copyfile(carrier/'system/blockMeshDict',dest/'system/blockMeshDict')
                for field in ['T','p','rho','U','phi']:
                    shutil.copyfile(latest/field,dest/'0'/field)
                meta['carrier']=read(carrier/'benchmark.json')
                meta['carrier_audit']=read(carrier/'carrier-audit.json')
                meta['carrier_projection']=read(carrier/'carrier-projection.json')
            cp=dest/'system/controlDict'; cp.write_text(control(cp.read_text(),ncorr))
            # Changes must not touch source/thermochemistry/transport dictionaries.
            for name0 in ['thermochemistryProperties','speciesTransportProperties','thermophysicalProperties']:
                require(sha(dest/'constant'/name0)==sha(source/'constant'/name0),'Physical input changed: '+name0)
            meta['followup']={'base_case':base,'category':'numerical','changes':{'study':kind,'correctors':ncorr}}
            save(dest/'benchmark.json',meta)
            entries.append({'case':name,'base_case':base,'kind':kind,'correctors':ncorr,'input_sha256':tree(dest)})
    save(out/'study.json',{'created_utc':utc(),'solver_sha256':EXPECTED,'prepared':str(prepared),
                         'cases':entries,'package_sha256':sha(PACKAGE/'package-manifest.json')})
    return out


def run_case(out,entry,state):
    name=entry['case']; inputs=out/'inputs'/name; case=out/'cases'/name
    require(tree(inputs)==entry['input_sha256'],'Prepared inputs changed: '+name)
    if state['cases'].get(name,{}).get('status')=='PASS':
        require(sha(case/'benchmark.json')==sha(inputs/'benchmark.json'),'Completed case manifest changed: '+name)
        require(tree(case/'0')==tree(inputs/'0') and tree(case/'constant')==tree(inputs/'constant') and
                tree(case/'system')==tree(inputs/'system'),'Completed case inputs changed: '+name)
        result,_=analyse(case); require(result['status']=='PASS','Previously completed case no longer passes')
        return
    def record(status,**values):
        with LOCK:
            state['cases'][name]={'status':status,**values};state['updated_utc']=utc();save(out/'run_state.json',state)
        print(name,status,flush=True)
    try:
        refuse_live_solver(case)
        refuse_live_solver(out/'smoke'/name)
        if case.exists():
            archive=out/'interrupted';archive.mkdir(exist_ok=True)
            case.rename(archive/(name+'-'+datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')))
        smoke=out/'smoke'/name
        if smoke.exists():
            archive=out/'interrupted';archive.mkdir(exist_ok=True)
            smoke.rename(archive/(name+'-smoke-'+datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')))
        shutil.copytree(inputs,smoke)
        cp=smoke/'system/controlDict'; cp.write_text(control(cp.read_text(),entry['correctors'],.02,True))
        meta=read(smoke/'benchmark.json');meta['end_time_s']=.02;save(smoke/'benchmark.json',meta)
        record('SMOKE_RUNNING')
        solver_call(smoke,out/'bin/rhoFixedFlowFoam')
        result,_=analyse(smoke);require(result['status']=='PASS','Ten-step smoke failed')
        shutil.copytree(inputs,case);record('RUNNING')
        solver_call(case,out/'bin/rhoFixedFlowFoam')
        with (case/'log.solver').open('rb') as f:
            f.seek(max(0,f.seek(0,2)-4096));require(f.read().rstrip().endswith(b'End'),'Solver log lacks normal End')
        result,_=analyse(case);require(result['status']=='PASS','Full-run numerical gates failed')
        save(case/'material-audit.json',ledger_audit(case))
        save(case/'log-audit.json',log_audit(case/'log.solver',read(case/'benchmark.json')))
        record('PASS',finished_utc=utc())
    except Exception as exc:
        record('FAIL',reason=str(exc));raise


def collect(out):
    study=read(out/'study.json');state=read(out/'run_state.json')
    report=out/('return-results-'+datetime.datetime.now().strftime('%Y%m%d%H%M%S%f'))
    report.mkdir();comparisons=[]
    for entry in study['cases']:
        if state['cases'].get(entry['case'],{}).get('status')!='PASS': continue
        case=out/'cases'/entry['case'];dest=report/entry['case'];dest.mkdir(exist_ok=True)
        for filename in ['benchmark.json','numerical-report.json','material-audit.json','log-audit.json','iodine-comparison.csv']:
            shutil.copyfile(case/filename,dest/filename)
        for path in (case/'postProcessing/axialProfiles/60').glob('*.dat'):shutil.copyfile(path,dest/path.name)
        for species in ['element_I','PbI2_g','BiI3_g']:
            filename=species+'_native.dat';a=rows(dest/filename)
            current=np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in a])
            references=[('medium',PACKAGE/'references'/entry['base_case']/filename)]
            if entry['kind']=='grid864x108':references.append(('fine',PACKAGE/'fine-references'/entry['base_case']/filename))
            for level,path in references:
                b=rows(path);reference=np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in b])
                require(len(a)==len(b) and np.allclose([r['x'] for r in a],[r['x'] for r in b]),'Profile coordinates differ')
                comparisons.append({'case':entry['case'],'reference':level,'species':species,
                    'inventory_L1_pct':float(100*np.abs(current-reference).sum()/reference.sum()) if reference.sum() else None,
                    'total_change_pct':float(100*(current.sum()-reference.sum())/reference.sum()) if reference.sum() else None})
    save(report/'comparisons.json',comparisons)
    save(report/'run_state.json',state);save(report/'study.json',study)
    save(report/'sha256.json',{str(p.relative_to(report)):sha(p) for p in report.rglob('*') if p.is_file() and p.name!='sha256.json'})
    archive=out/'lesto-followup-results.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:tar.add(report,arcname='lesto-followup-results')
    print('Return archive:',archive,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['check','run','status','collect'])
    p.add_argument('--repo',type=Path,default=Path.home()/'LESTO')
    p.add_argument('--prepared',type=Path);p.add_argument('--solver',type=Path);p.add_argument('--carrier-tool',type=Path)
    p.add_argument('--out',type=Path,default=Path.home()/'lesto-submission-followup')
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--resume',action='store_true')
    args=p.parse_args();args.out=args.out.resolve()
    require(1<=args.jobs<=16,'--jobs must be between 1 and 16')
    if args.command=='check':check(args);return
    if args.command=='status':
        path=args.out/'run_state.json';print(json.dumps(read(path) if path.exists() else {'status':'PREPARING_OR_NOT_STARTED'},indent=2));return
    if args.command=='collect':collect(args.out);return
    # Lock held by this process; automatically released on exit. No overlapping resume.
    import fcntl
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with (args.out.parent/(args.out.name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.resume:
            check(args);require((args.out/'study.json').is_file(),'Preparation incomplete; choose a fresh output directory')
            require(read(args.out/'study.json')['package_sha256']==sha(PACKAGE/'package-manifest.json'),'Package changed since preparation')
        else:prepare(args)
        require(sha(args.out/'bin/rhoFixedFlowFoam')==EXPECTED,'Private solver fingerprint changed')
        study=read(args.out/'study.json')
        state=read(args.out/'run_state.json') if args.resume and (args.out/'run_state.json').exists() else {'cases':{}}
        for entry in study['cases']:state['cases'].setdefault(entry['case'],{'status':'QUEUED'})
        state.update(status='RUNNING',updated_utc=utc());save(args.out/'run_state.json',state)
        failed=False
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures=[pool.submit(run_case,args.out,e,state) for e in study['cases']]
            for future in concurrent.futures.as_completed(futures):
                try:future.result()
                except Exception as exc:failed=True;print('FAILED:',exc,flush=True)
        state['status']='FAILED' if failed else 'COMPLETE';save(args.out/'run_state.json',state)
        collect(args.out)
        if failed:raise SystemExit(1)


if __name__=='__main__':
    try: main()
    except Exception as exc:
        print('ERROR:',exc,file=sys.stderr);raise SystemExit(1)
