#!/usr/bin/env python3
"""Regenerate default and extended systems; compare every default output."""
from pathlib import Path
import subprocess
import sys
work,nasa,tf=sys.argv[1:]
repo=Path(__file__).resolve().parents[4]
systems=repo/'thermochemistry/systems'
out=Path(work)/'generatorPbI2'
cmd=[sys.executable,str(systems/'make_gems3k_input.py'),'--nasa',nasa,'--out',str(out),'--json','--thermofun']
subprocess.run(cmd,check=True)
files=list(out.iterdir())
assert len(files)>=10
for p in files:
    original=systems/'PbI2He'/p.name
    assert original.exists() and p.read_bytes()==original.read_bytes(),p
subprocess.run(cmd[:cmd.index('--out')]+['--out',str(Path(work)/'generatorPbBiI'),'--system','PbBiIHe','--tf',tf,'--json'],check=True)
import json
record=json.loads((Path(work)/'generatorPbBiI/PbBiIHe-dch.json').read_text())[0]['dch']
assert record['nIC']==4
assert len(record['A']) == record['nIC']*record['nDC']
assert 'Bi' in record['ICNL']
print(f'{len(files)} PbI2He files byte-identical; PbBiIHe has 4 elements and validated element-matrix dimensions')
