#!/usr/bin/env python3
"""Independent Decimal reference: prescribe potentials, then sum element totals.
Invented formation constants only; no external thermodynamic records.
"""
from decimal import Decimal as D, localcontext
import math
import sys
FORMULAS=[(1,0,0),(1,0,1),(1,0,2),(0,1,0),(0,2,0),(0,1,1),(0,1,3),(0,0,1),(0,0,2)]
K=[0,12,26,0,15,16,32,0,8]
with open(sys.argv[1],'w') as out,localcontext() as ctx:
    ctx.prec=70
    out.write('# SYNTHETIC Decimal reference; T,scale,bPb,bBi,bI,(lnKf,c) per species\n')
    for T in [350,400,450,700,1000,1173,1250]:
        scale=float(1e5/(8.314462618*T))
        lnK=[float(k*700/T*math.log(10)) for k in K]
        for level in [-30,-25,-20,-15]:
            for absent in [None,0,1,2,'all']:
                iodine=D(level)
                lead=D(level)-max(D.from_float(k)+nu[2]*iodine for nu,k in zip(FORMULAS,lnK) if nu[0])
                bismuth=min((D(level)-D.from_float(k)-nu[2]*iodine)/nu[1]
                            for nu,k in zip(FORMULAS,lnK) if nu[1])
                pi=[lead,bismuth,iodine]
                c=[]
                for nu,k in zip(FORMULAS,lnK):
                    if absent=='all' or (isinstance(absent,int) and nu[absent]):c.append(D(0))
                    else:c.append(D.from_float(scale)*(D.from_float(k)+sum(n*x for n,x in zip(nu,pi))).exp())
                b=[float(sum(nu[e]*x for nu,x in zip(FORMULAS,c))) for e in range(3)]
                row=[T,scale,*b]
                for k,x in zip(lnK,c):row.extend([k,float(x)])
                out.write(','.join(format(x,'.17g') for x in row)+'\n')
