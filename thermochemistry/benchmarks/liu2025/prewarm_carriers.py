#!/usr/bin/env python3
"""Build independent study carriers concurrently and publish only audited outputs.

Staging directories prevent the study preparer from seeing partial carriers.
If preparation already started the same destination, keep staging as an archive.
"""
import argparse
import concurrent.futures
import csv
import json
from pathlib import Path

from carriers import REPO, call, prepare_carrier


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=6);args=p.parse_args();out=args.out.resolve()
    targets=[('grid576x72','LBE-I_SS_II',576,72,298.15,None)]
    targets += [('flow0C',experiment,384,48,273.15,None) for experiment in ['PbI2_SS','LBE-I_SS_II','LBE-I_SiO2_I']]
    diagnostics=json.loads((REPO/'run/034-liu2025/full_analysis/analysis.json').read_text())['silica_figure_diagnostics']
    for d in diagnostics:
        source=REPO/'thermochemistry/benchmarks/liu2025/digitized/LBE-I_SiO2_I_temperature.csv'
        root=out/'carrier_staging/temperature_hypotheses'/d['species'];root.mkdir(parents=True,exist_ok=True)
        shift=d['diagnostic_offset_cm'][0]/100
        with source.open() as stream:points=list(csv.DictReader(stream))
        with (root/source.name).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=['x_m','T_K'],lineterminator='\n');writer.writeheader()
            writer.writerows({'x_m':float(row['x_m'])-shift,'T_K':row['T_K']} for row in points if 0<=float(row['x_m'])-shift<=1.25)
        targets.append(('silica_'+d['species'],'LBE-I_SiO2_I',384,48,298.15,root))
    def job(target):
        name,experiment,nx,nr,flow_T,temperature_root=target
        destination=out/'carriers'/name/experiment
        if destination.exists():return {'carrier':str(destination),'status':'EXISTING_SKIP'}
        staged=out/'carrier_staging'/name/experiment
        prepare_carrier(staged,experiment,nx,nr,flow_T,101325,temperature_root)
        tool=str(out/'bin/liuCarrierTool')
        for command,label in [(['blockMesh'],'blockMesh'),(['checkMesh','-constant'],'checkMesh'),
                              ([tool],'initialise'),(['rhoSimpleFoam'],'carrier')]:call(staged,command,label)
        times=sorted((path for path in staged.iterdir() if path.is_dir() and path.name.replace('.','',1).isdigit()),key=lambda path:float(path.name))
        for flag,label in [('-project','project'),('-audit','audit')]:call(staged,[tool,flag,'-time',times[-1].name],label)
        assert json.loads((staged/'carrier-audit.json').read_text())['status']=='PASS'
        destination.parent.mkdir(parents=True,exist_ok=True)
        try:
            if destination.exists():status='AUDITED_STAGING_SUPERSEDED'
            else:staged.rename(destination);status='AUDITED_PUBLISHED'
        except OSError:
            if not destination.exists():raise
            status='AUDITED_STAGING_SUPERSEDED'
        print(experiment,name,status,flush=True)
        return {'carrier':str(destination),'status':status}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:results=list(pool.map(job,targets))
    (out/'carrier_prewarm.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
