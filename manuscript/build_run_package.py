#!/usr/bin/env python3
"""Package a frozen scientific runner; never include external coefficient tables."""
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

HERE=Path(__file__).resolve().parent
REPO=HERE.parent


def main():
    stage=HERE/'build'/('run-package-'+datetime.datetime.now().strftime('%Y%m%d%H%M%S%f'))/'lesto-submission-run'
    stage.mkdir(parents=True)
    def copy(source,dest,text=False):
        dest=stage/dest;dest.parent.mkdir(parents=True,exist_ok=True)
        if text:dest.write_bytes(source.read_bytes().replace(b'\r\n',b'\n'))
        else:shutil.copyfile(source,dest)
    for p in (HERE/'run-package').iterdir():
        if p.is_file():copy(p,Path(p.name),True)
    copy(HERE/'next-study.json',Path('next-study.json'))
    bench=Path('thermochemistry/benchmarks/liu2025')
    for name in ['carriers.py','prepare.py','report.py','analyse_matrix.py']:
        copy(REPO/bench/name,Path('payload')/bench/name,True)
    for p in (REPO/bench/'carrierTool').rglob('*'):
        if p.is_file():copy(p,Path('carrierTool')/p.relative_to(REPO/bench/'carrierTool'),True)
    for p in (REPO/bench).glob('*.csv'):copy(p,Path('payload')/bench/p.name)
    for p in (REPO/bench/'digitized').iterdir():
        if p.suffix in ['.csv','.json']:copy(p,Path('payload')/bench/'digitized'/p.name)
    prop=Path('run/007-variable-properties-high-temperature/constant/thermophysicalProperties')
    copy(REPO/prop,Path('payload')/prop)
    bases=['PbI2_SS_Q1_baseline','LBE-I_SS_II_Q2_baseline_sat0.001',
           'LBE-I_SS_II_Q2_baseline_sat0.01','LBE-I_SS_II_Q2_shift29k_sat0.01']
    for base in bases:
        reference=REPO/'run/024-liu2025/followup_interim_20261008/reference_profiles'/base
        copy(reference/'benchmark.json',Path('references')/base/'benchmark.json')
        for species in ['element_I','PbI2_g','BiI3_g']:
            name=species+'_native.dat'
            copy(REPO/'run/024-liu2025/full_analysis/profiles'/base/name,Path('references')/base/name)
            copy(REPO/'run/024-liu2025/followup_results/profiles'/(base+'__grid576x72')/name,
                 Path('fine-references')/base/name)
    hashes={p.relative_to(stage).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(stage.rglob('*')) if p.is_file()}
    (stage/'package-manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    output=HERE/'lesto-submission-run.tar.gz'
    with tarfile.open(output,'w:gz') as tar:
        for p in sorted(stage.rglob('*')):
            if p.is_file():tar.add(p,arcname='lesto-submission-run/'+p.relative_to(stage).as_posix())
    print(stage)
    print(output)
    print('SHA256',hashlib.sha256(output.read_bytes()).hexdigest())


if __name__=='__main__':main()
