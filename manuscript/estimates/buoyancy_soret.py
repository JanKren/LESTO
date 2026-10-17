"""Buoyancy and thermal-diffusion bounds for the manuscript discussion.

Inputs (case inputs, App. A / Sec. 4.1): helium perfect gas, M=4.0026e-3,
Cp=5193, Sutherland As=1.4683e-6, Ts=83.40 (OpenFOAM modified Eucken
conductivity -> Pr = Cp/(Cv(1.32+1.77R/Cv)) = 0.666); absolute pressure
101.33 kPa (case 0/p); mass flows 2.727e-7 (100 mL/min) and 1.227e-7 kg/s
(45 mL/min); d = 4.8 / 5.0 mm; D = 1e-4 m2/s; digitised wall curves;
largest baseline trace mole fraction 9.19e-5 (generated/results.json).
Thermal diffusion factor alpha_T taken as 1 (generous upper bound).
"""
import csv, numpy as np
Ru=8.314462618; M=4.0026e-3; Rs=Ru/M; Cp=5193.; Cv=Cp-Rs
As,Ts=1.4683e-6,83.40; p=101330.; g=9.81; D=1e-4
mu=lambda T: As*np.sqrt(T)/(1+Ts/T)
kap=lambda T: mu(T)*Cv*(1.32+1.77*Rs/Cv)
rho=lambda T: p/(Rs*T)
print('Pr =', Cp*mu(700)/kap(700))
import pathlib
B=str(pathlib.Path(__file__).resolve().parents[2]/'thermochemistry/benchmarks/liu2025/digitized')+'/'
cases={'PbI2_SS':(4.8e-3,2.727e-7),'LBE-I_SS_II':(4.8e-3,1.227e-7),'LBE-I_SiO2_I':(5.0e-3,1.227e-7)}
out={}
for name,(d,m) in cases.items():
    r=list(csv.DictReader(open(B+name+'_temperature.csv')))
    x=np.array([float(v['x_m']) for v in r]); T=np.array([float(v['T_K']) for v in r])
    dTdx=np.gradient(T,x)
    # smooth over 1 cm to suppress digitisation noise
    k=max(1,int(round(0.01/np.median(np.diff(x))))); dTdx=np.convolve(dTdx,np.ones(k)/k,'same')
    R=d/2; A=np.pi*R**2
    U=m/(rho(T)*A); nu=mu(T)/rho(T); alpha=kap(T)/(rho(T)*Cp); Re=4*m/(np.pi*d*mu(T))
    # developed laminar flow, linear wall T: T_w - T_bulk = (11/192) U d^2 |dT/dx| / alpha
    dTr=11/192*U*d**2*np.abs(dTdx)/alpha
    Gr_r=g*dTr/T*d**3/nu**2               # radial thermal Grashof
    Gr_ax=g*np.abs(dTdx)*d/T*d**3/nu**2   # axial temperature change over one diameter
    xmax=9.19e-5; Mh=0.58969381           # BiI3, heaviest transported iodide
    Gr_c=g*xmax*(Mh/M-1)*d**3/nu**2        # compositional, mole-fraction bound
    soret_rad=dTr/T                        # alpha_T*dT_r/T vs Fick flux x D/R
    soret_ax=D*np.abs(dTdx)/T/U            # Soret drift / carrier velocity
    sel=(T>350)
    out[name]=dict(maxdTdx=np.abs(dTdx).max(),dTr=dTr.max(),Ri_r=(Gr_r/Re**2)[sel].max(),
        Ri_ax=(Gr_ax/Re**2)[sel].max(),Gr_ax=Gr_ax[sel].max(),Ri_c=(Gr_c/Re**2)[sel].max(),
        soret_rad=soret_rad.max(),soret_ax=soret_ax.max())
for k,v in out.items(): print(k,{a:float(f'{b:.3g}') for a,b in v.items()})

# where, and in the deposition zone (wall T 400-750 K)
for name,(d,m) in cases.items():
    r=list(csv.DictReader(open(B+name+'_temperature.csv')))
    x=np.array([float(v['x_m']) for v in r]); T=np.array([float(v['T_K']) for v in r])
    dTdx=np.gradient(T,x); k=max(1,int(round(0.01/np.median(np.diff(x))))); dTdx=np.convolve(dTdx,np.ones(k)/k,'same')
    R=d/2; A=np.pi*R**2; U=m/(rho(T)*A); nu=mu(T)/rho(T); Re=4*m/(np.pi*d*mu(T))
    Gr_ax=g*np.abs(dTdx)*d/T*d**3/nu**2; Ri=Gr_ax/Re**2
    i=np.argmax(np.where(T>350,Ri,0))
    z=(T>400)&(T<750)
    print(name,'max Ri_ax at x=%.3f m T=%.0f K dTdx=%.0f'%(x[i],T[i],dTdx[i]),
          '| deposition zone: Gr_ax<=%.2f Ri_ax<=%.2f Re %.2f-%.2f'%(Gr_ax[z].max(),Ri[z].max(),Re[z].min(),Re[z].max()))
