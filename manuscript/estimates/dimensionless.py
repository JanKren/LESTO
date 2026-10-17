"""Dimensionless numbers for the LESTO manuscript framing paragraph.

Inputs (all stated in the manuscript):
  diameters 4.8 mm (steel), 5.0 mm (silica)                Sec. 2, Tab. 1
  helium mass flow 2.727e-7 kg/s (100 mL/min, steel PbI2),
                   1.227e-7 kg/s (45 mL/min, LBE)          Sec. 4.1
  Sutherland helium: A_s=1.4683e-6, T_s=83.40 K             App. A
  M_He = 0.0040026 kg/mol                                   App. A
  D_i = 1e-4 m2/s (baseline, all gases)                     Sec. 4.2
  HKS conductance, alpha in {1, 0.1, 0.01}, C_e=s_i=1       Sec. 3.3, 4.2
  M_PbI2 = 0.46100894, M_BiI3 = 0.58969381 kg/mol           Tab. A1
  temperature range: 400 K (127 C, lowest printed peak label) to
                     1173.15 K (T_h, source furnace)        Sec. 3.5, 7.2
  axial cell size 1.25 m / 384 (medium mesh)                Sec. 4.1
Not stated in the manuscript (used only for the optional Pe/Sc lines):
  column absolute pressure 101325 Pa (case 0/p internalField).
"""
import numpy as np

R = 8.314462618
As, Ts = 1.4683e-6, 83.40
M_He = 0.0040026
D = 1e-4
mu = lambda T: As*np.sqrt(T)/(1+Ts/T)
T = np.array([400.0, 1173.15])
cases = {"PbI2 steel 100 mL/min": (4.8e-3, 2.727e-7),
         "LBE steel 45 mL/min": (4.8e-3, 1.227e-7),
         "LBE silica 45 mL/min": (5.0e-3, 1.227e-7)}
print("mu(T) [Pa s]:", mu(T))
print("\nReynolds number Re = 4 mdot/(pi d mu)  (pressure-independent)")
for k, (d, m) in cases.items():
    print(f"  {k:24s} Re(400K)={4*m/(np.pi*d*mu(T[0])):.2f}  Re(1173K)={4*m/(np.pi*d*mu(T[1])):.2f}")
print("\nRadial diffusion time R^2/D (pressure-independent)")
for d in (4.8e-3, 5.0e-3):
    print(f"  d={d*1e3} mm: {(d/2)**2/D*1e3:.1f} ms")
print("\nWall Damkoehler Da = G chi R/(rho D) = f(alpha) sqrt(R_u T/(2 pi M)) R/D")
print("  (ratio of radial-diffusion resistance R/(rho D) to interfacial resistance 1/(G chi);")
print("   chi = rho R_u T/M so rho cancels: pressure-independent)")
for name, M in (("PbI2", 0.46100894), ("BiI3", 0.58969381)):
    for alpha in (1.0, 0.1, 0.01):
        f = 2*alpha/(2-alpha)
        vals = [f*np.sqrt(R*t/(2*np.pi*M))*(d/2)/D for t in T for d in (4.8e-3, 5.0e-3)]
        print(f"  {name} alpha={alpha:<5}: Da = {min(vals):.3g} .. {max(vals):.3g}")
print("\nNumerical diffusion threshold: |u_x| dx/2 = D  ->  u_x =", 2*D/(1.25/384), "m/s")

print("\n--- optional, needs column pressure p = 101325 Pa (NOT stated in paper) ---")
p = 101325.0
rho = lambda T: p*M_He/(R*T)
for k, (d, m) in cases.items():
    A = np.pi*d**2/4
    for t in T:
        U = m/(rho(t)*A); Sc = mu(t)/(rho(t)*D); Pe = U*d/D
        Re = 4*m/(np.pi*d*mu(t))
        print(f"  {k:24s} T={t:7.2f} K  U={U:.3f} m/s  Sc={Sc:.3f}  Pe_d={Pe:.2f}  Re*Sc={Re*Sc:.2f}"
              f"  L_D=U R^2/D={U*(d/2)**2/D*100:.2f} cm")
