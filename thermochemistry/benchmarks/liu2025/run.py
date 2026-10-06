#!/usr/bin/env python3
"""Run prepared Liu cases, recording a log and report for each completed case."""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from carriers import call
from report import analyse


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--solver',default='rhoFixedFlowFoam')
    parser.add_argument('--case',action='append',help='Case names; omit to run all prepared cases')
    parser.add_argument('--jobs',type=int,default=1,help='Independent serial cases to run concurrently')
    args=parser.parse_args()
    matrix=json.loads((args.root/'matrix.json').read_text())
    names={r['case'] for r in matrix if r['status']=='PREPARED_PROVISIONAL'}
    selected=set(args.case) if args.case else names
    if not selected<=names:
        parser.error('unknown or unprepared case: '+str(sorted(selected-names)))
    if args.jobs<1:
        parser.error('--jobs must be positive')
    def run(name):
        case=(args.root/name).resolve()
        call(case,[args.solver],'solver')
        result,_=analyse(case)
        print(case.name,result['status'],flush=True)
        return result['status']
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        statuses=list(pool.map(run,[r['case'] for r in matrix if r['case'] in selected]))
    if any(status!='PASS' for status in statuses):
        raise SystemExit('Element closure or solver-defect gate failed')


if __name__=='__main__':
    main()
