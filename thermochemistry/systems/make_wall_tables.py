#!/usr/bin/env python3
"""M10d channel tables, dimensionless disproportionation K and Q2 source.

External NASA/ThermoFun inputs use the existing common evaluator. No records
are bundled. --bii-shift changes G(BiI_g) by a constant J/mol and therefore
changes both formation and wall-reaction constants consistently.
"""
import argparse,json,math
from pathlib import Path
from make_pv_tables import Data,TABLES,write_table,crossing,R,P0
from make_kf_tables import SPECIES,formation,write as write_kf
from make_gems3k_input import pbbii_system

def generate(nasa,tf,out,shift=0,saturation=.01):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    data=Data(nasa,tf);original=data.G
    data.G=lambda name,T:original(name,T)+(shift if name=='BiI(g)' else 0)
    for filename,gas,cols,nu,bracket,provenance in TABLES:
        cross=crossing(data,gas,cols,nu,*bracket) if bracket else None
        write_table(data,out/filename,gas,cols,nu,cross,provenance+[f'BiI(g) constant Gibbs shift {shift:g} J/mol'])
    cols=[('PbI2(cr)','PbI2(cr)'),('PbI2(l)','PbI2(l)')]
    write_table(data,out/'pv_PbI2_phases.csv','PbI2(g)',cols,1,
                crossing(data,'PbI2(g)',cols,1,600,750),['NASA/Gurvich PbI2 records'])
    system=pbbii_system(nasa,tf);fun=formation(system)
    def constants(name,T):return fun(name,T)-(shift/(R*T*math.log(10)) if name=='BiI_g' else 0)
    write_kf(out/'log10Kf_PbBiI.csv',list(SPECIES),constants)
    with (out/'logK_BiI_disp.csv').open('w') as f:
        f.write('# 3 BiI(g) = 2 Bi(stable) + BiI3(g); dimensionless K, gas p0=1 bar\n# T_K,logK\n')
        grid=sorted(set([270+5*i for i in range(197)]+[544.51]))
        for T in grid:
            dg=data.G('BiI3(g)',T)+2*min(data.G('Bi(cr)',T),data.G('Bi(l)',T))-3*data.G('BiI(g)',T)
            f.write(f'{T:.17g},{-dg/(R*T*math.log(10)):.17g}\n')
    T=1173;P=101325;aPb=.447;aBi=.553
    pPb=aPb*P0*math.exp(-(data.G('Pb(g)',T)-min(data.G('Pb(cr)',T),data.G('Pb(l)',T)))/(R*T))
    metalBi=min(data.G('Bi(cr)',T),data.G('Bi(l)',T))
    pBi=aBi*P0*math.exp(-(data.G('Bi(g)',T)-metalBi)/(R*T))
    pBi2=aBi*aBi*P0*math.exp(-(data.G('Bi2(g)',T)-2*metalBi)/(R*T))
    nHe=45e-6*180*101325/(R*273.15)
    Pb=.111e-6/.46100894;Bi=.181e-6/.58969381
    bulk={'Pb':Pb/nHe+saturation*pPb/P,'Bi':Bi/nHe+saturation*(pBi+2*pBi2)/P,'I':(2*Pb+3*Bi)/nHe}
    result={'status':'provisional numerical source, ideal LBE activities','Tin_K':T,'P_Pa':P,'saturation_fraction':saturation,'BiI_G_shift_J_mol':shift,
            'metal_saturation_Pa':{'Pb':pPb,'Bi':pBi,'Bi2':pBi2},'element_mol_per_mol_He':bulk,
            'iodine_basis':'Liu LBE-I_SS_II deposited amounts, 45 mL/min STP over 180 min'}
    (out/'q2.json').write_text(json.dumps(result,indent=2)+'\n');return result

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--nasa',required=True);ap.add_argument('--tf',required=True);ap.add_argument('--out',required=True);ap.add_argument('--bii-shift',type=float,default=0);ap.add_argument('--saturation',type=float,default=.01)
    a=ap.parse_args()
    if not math.isfinite(a.bii_shift) or not math.isfinite(a.saturation) or not 0<=a.saturation<=1:ap.error('Require finite shift and saturation in [0,1]')
    print(json.dumps(generate(a.nasa,a.tf,a.out,a.bii_shift,a.saturation),indent=2))
if __name__=='__main__':main()
