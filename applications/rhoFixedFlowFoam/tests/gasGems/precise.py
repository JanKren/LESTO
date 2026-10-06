"""Independent high-precision full Newton solve of gas element balances.
The bridge's double-precision KKT residual can be amplified in nearly
exhausted iodine. Refine its potentials against the exact fixed-volume
input in Decimal, using a general Jacobian rather than the kernel's nested
one-dimensional solvers. No thermodynamic data are shipped here.
"""
from decimal import Decimal as D, localcontext
import math
FORMULAS=[(1,0,0),(1,0,1),(1,0,2),(0,1,0),(0,2,0),(0,1,1),(0,1,3),(0,0,1),(0,0,2)]
def refine(ks,scale,b,start):
    with localcontext() as ctx:
        ctx.prec=70;active=[a for a in range(3) if b[a]>0]
        target=[D.from_float(v) for v in b];constants=list(map(D.from_float,ks));scale=D.from_float(scale)
        pi=list(map(D.from_float,start))
        for step in range(100):
            c=[scale*(k+sum(n*pi[a] for a,n in enumerate(nu))).exp()
               if all(not nu[a] for a in range(3) if a not in active) else D(0)
               for k,nu in zip(constants,FORMULAS)]
            totals=[sum(nu[a]*v for nu,v in zip(FORMULAS,c)) for a in active]
            if max(abs(totals[i]/target[a]-1) for i,a in enumerate(active))<D('1e-50'):return list(map(float,c))
            matrix=[[sum(nu[a]*nu[z]*v for nu,v in zip(FORMULAS,c)) for z in active]+[target[a]-totals[i]] for i,a in enumerate(active)]
            n=len(active)
            for col in range(n):
                pivot=max(range(col,n),key=lambda r:abs(matrix[r][col]));matrix[col],matrix[pivot]=matrix[pivot],matrix[col]
                value=matrix[col][col];assert value!=0
                matrix[col]=[v/value for v in matrix[col]]
                for row in range(n):
                    if row==col:continue
                    factor=matrix[row][col];matrix[row]=[v-factor*w for v,w in zip(matrix[row],matrix[col])]
            for i,a in enumerate(active):pi[a]+=matrix[i][-1]
        raise AssertionError('Decimal reference did not converge')
