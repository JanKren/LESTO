#!/usr/bin/env python3
"""Generate gas formation constants from monatomic gases, standard state 1 bar.

Uses the same NASA/ThermoFun evaluation as make_gems3k_input.py. Bi data
remain external. --synthetic writes invented constants for code tests only.
CSV columns are solver species names; values are log10 Kf, interpolated in
1/T. No condensed phase contributes to this homogeneous equilibrium.
"""
import argparse
import math
from pathlib import Path
from make_gems3k_input import R, species_props, pbi2_system, pbbii_system

SPECIES = {'Pb_g': {'Pb': 1}, 'PbI_g': {'Pb': 1, 'I': 1},
           'PbI2_g': {'Pb': 1, 'I': 2}, 'Bi_g': {'Bi': 1},
           'Bi2_g': {'Bi': 2}, 'BiI_g': {'Bi': 1, 'I': 1},
           'BiI3_g': {'Bi': 1, 'I': 3}, 'I_g': {'I': 1}, 'I2_g': {'I': 2}}
SYNTHETIC = {'Pb_g':0, 'PbI_g':12, 'PbI2_g':26, 'Bi_g':0,
             'Bi2_g':15, 'BiI_g':16, 'BiI3_g':32, 'I_g':0, 'I2_g':8}

def write(path, names, function, step=5):
    grid = [270.0+i*step for i in range(math.floor((1250-270)/step)+1)]
    if grid[-1] != 1250: grid.append(1250.0)
    with Path(path).open('w') as f:
        f.write('# Gas formation constants; standard state p0 = 1 bar\n')
        f.write('# '+('SYNTHETIC invented constants; code tests only' if function is synthetic else
                    'NASA Pb-I; external ThermoFun Bi-I; same evaluator as the GEMS3K generator')+'\n')
        f.write('# T_K,'+','.join(names)+'\n')
        for T in grid:
            f.write(','.join(format(x,'.17g') for x in [T]+[function(n,T) for n in names])+'\n')

def synthetic(name,T):
    return SYNTHETIC[name]*700/T

def formation(system):
    species = {s['name']:s for ph in system['phases'] for s in ph['species']}
    def evaluate(name,T):
        s=species[name.removesuffix('_g')+'(g)']
        mono=sum(n*species_props(species[e+'(g)'],T)[3] for e,n in s['formula'].items())
        return -(species_props(s,T)[3]-mono)/(R*T*math.log(10))
    return evaluate

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--nasa'); ap.add_argument('--tf'); ap.add_argument('--synthetic',action='store_true')
    ap.add_argument('--out',required=True);ap.add_argument('--step',type=float,default=5)
    a=ap.parse_args()
    if not math.isfinite(a.step) or a.step<=0:ap.error('--step must be finite and positive')
    if a.synthetic:names=list(SPECIES);fun=synthetic
    else:
        if not a.nasa:ap.error('--nasa is required without --synthetic')
        system=pbbii_system(a.nasa,a.tf) if a.tf else pbi2_system(a.nasa)
        names=[n for n in SPECIES if a.tf or 'Bi' not in n];fun=formation(system)
    write(a.out,names,fun,a.step)
if __name__=='__main__':main()
