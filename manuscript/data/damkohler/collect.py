#!/usr/bin/env python3
"""Damkoehler sweep of the half-cell HKS channel (case 016): per-run metrics.

Reads the runs of run/016-hks-channel-halfcell/Allrun-damkohler and writes
the tracked records of this folder:

    python3 manuscript/data/damkohler/collect.py \
        run/016-hks-channel-halfcell/work/damkohler/ny24 \
        run/016-hks-channel-halfcell/work/damkohler/ny48

sweep_summary.json  per run: alpha, Da, ledger fractions, onset/peak/centroid,
                    wall profile at 3 s; at t = 1.5 s the per-column wall flux
                    J, mean, wall and equilibrium concentrations, J^, eta and
                    local Da, and three wall-normal profiles (72.5, 84.5,
                    94.5 mm)
damkohler_sweep.csv the table of the 100 x 24 sweep
"""
import csv, glob, json, os, re, sys
import numpy as np
import foamio

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '../../..'))
PV_TABLE = os.path.join(REPO, 'thermochemistry/systems/data/pv_PbI2_phases.csv')

Ru = 8.314462618       # J/(mol K)
M = 0.46100894         # PbI2 molar mass [kg/mol]
D = 1e-4               # PbI2 diffusivity of case 016 [m2/s]
h = 2.4e-3             # channel half-height [m] (wall-to-centreline diffusion length)
dx = 1e-3; NX = 100    # axial cells
N0 = 2.7992536908413373e-11   # released PbI2 [mol] (case 016 expect.json)
T_ONSET = 561.16       # 1 % onset of the alpha = 1 run [K], the same for every alpha
STATIONS = (72, 84, 94)       # axial columns of the wall-normal profiles (x = 72.5, 84.5, 94.5 mm)

rows = [l for l in open(PV_TABLE) if not l.startswith('#')]
TAB = np.array([[float(v) for v in r.split(',')] for r in rows[1:]])


def p_eq(T):
    """Equilibrium pressure [Pa]: minimum over the condensed phases of the table."""
    return 10**np.interp(T, TAB[:, 0], np.min(TAB[:, 1:], axis=1))


def Gk(alpha, T):
    """Kinetic conductance G_e beta/rho [m/s] of the HKS law with C_e = s = 1."""
    return 2*alpha/(2 - alpha)*np.sqrt(Ru*T/(2*np.pi*M))


def Da(alpha, T):
    return Gk(alpha, T)*h/D


def header(path):
    d = {}
    for l in open(path):
        if l.startswith('# inventory'):
            v = l.split(); d['gas'] = float(v[3]); d['wall'] = float(v[5])
        if l.startswith('# Tdep fromInlet'):
            m = re.search(r'status (\S+) T (\S+) x (\S+) peak (\S+) xPeak (\S+)', l)
            d['onset_status'] = m.group(1)
            if m.group(1) == 'found':
                d.update(T_onset=float(m.group(2)), x_onset=float(m.group(3)),
                         peak=float(m.group(4)), x_peak=float(m.group(5)))
        if l.startswith('# Tpeak'):
            d['T_peak'] = float(l.split()[2])
    a = np.loadtxt(path, comments='#'); d['x'] = a[:, 2]; d['wallprof'] = a[:, 4]
    return d


def run(case):
    ny = int(re.search(r'\(100 (\d+) 1\)', open(case + '/system/blockMeshDict').read()).group(1))
    dy = 4.8e-3/ny
    alpha = float(re.search(r'accommodation ([^;]+);',
                            open(case + '/constant/thermochemistryProperties').read()).group(1))
    log = open(case + '/log.solver').read()
    bal = re.findall(r'^Balance PbI2_g \[mol\]: gas (\S+) wall (\S+) .* released (\S+) transport (\S+)', log, re.M)[-1]
    clos = max(abs(float(v)) for v in re.findall(r'closure \S+ \(relative (\S+)\)', log))
    et = float(re.findall(r'ExecutionTime = (\S+)', log)[-1])
    hd = header(case + '/postProcessing/axialProfiles/3/PbI2_g_native.dat')
    x, w = hd['x'], hd['wallprof']; tot = w.sum()
    cen = (x*w).sum()/tot; sig = np.sqrt(((x - cen)**2*w).sum()/tot)
    cw = np.cumsum(w)/tot; x10, x50, x90 = (np.interp(q, cw, x) for q in (.1, .5, .9))
    T0 = foamio.read(case + '/0/T', ny=ny)['internal'][0]     # wall-cell temperature
    out = dict(case=os.path.basename(case), ny=ny, alpha=alpha,
               Da_at_561K=float(Da(alpha, T_ONSET)), Da_at_541K=float(Da(alpha, 541.75)),
               Da_400_1000=[float(Da(alpha, 400.)), float(Da(alpha, 1000.))],
               deposited_fraction=float(bal[1])/N0, gas_fraction=float(bal[0])/N0,
               outlet_fraction=float(bal[3])/N0, released_over_n0=float(bal[2])/N0,
               max_relative_closure=clos, execution_s=et, onset_status=hd['onset_status'],
               T_onset=hd.get('T_onset'), x_onset_mm=1e3*hd['x_onset'] if 'x_onset' in hd else None,
               T_peak=hd.get('T_peak'), x_peak_mm=1e3*hd['x_peak'] if 'x_peak' in hd else None,
               peak_mol_per_m=hd.get('peak'), centroid_mm=1e3*cen, sigma_mm=1e3*sig,
               width_10_90_mm=1e3*(x90 - x10), x50_mm=1e3*x50, wall_profile=w.tolist())
    t = '1.5'                                            # during the release (0-2 s)
    c = foamio.read(f'{case}/{t}/c_PbI2_g', ny=ny)['internal']
    ex = foamio.read(f'{case}/{t}/exch_PbI2_g', ny=ny)['internal']
    J = -0.5*(ex[0] + ex[-1])*dy/M                       # mol/(m2 s), both walls (symmetric)
    cc = 0.5*(c[0] + c[-1]); cb = c.mean(0); ceq = p_eq(T0)/(Ru*T0)
    cwall = cc - J*dy/(2*D)                              # half-cell split
    ceq_law = cwall - J/Gk(alpha, T0)                    # equilibrium implied by the solver's own law
    resid = np.where(J > 0, (ceq_law - ceq)/np.maximum(ceq, 1e-300), 0)   # table interpolation check
    ceq = np.where(J > 0, ceq_law, ceq)
    drive = cb - ceq
    ok = (J > 0) & (drive > 1e-6*cb.max())
    Jhat = np.where(ok, J*h/(D*np.where(ok, drive, 1)), np.nan)
    eta = np.where(ok, (cwall - ceq)/np.where(ok, drive, 1), np.nan)
    wgt = np.where(ok, J, 0)
    out['t' + t] = dict(
        J=J.tolist(), c_bulk=cb.tolist(), c_wall=cwall.tolist(), c_eq=ceq.tolist(),
        Jhat=Jhat.tolist(), eta=eta.tolist(), Da_loc=Da(alpha, T0).tolist(),
        zone_Jhat=float(np.sum(J[ok])*h/(D*np.sum(drive[ok]))) if ok.any() else None,
        zone_eta=float(np.sum(wgt*np.nan_to_num(eta))/wgt.sum()) if ok.any() else None,
        zone_Da=float(np.sum(wgt*Da(alpha, T0))/wgt.sum()) if ok.any() else None,
        ceq_table_vs_law_max_rel_diff=float(np.max(np.abs(resid[ok]))) if ok.any() else None,
        profiles={str(i): (c[:, i]/cb[i]).tolist() for i in STATIONS if cb[i] > 0})
    return out


def sweep(directory):
    cases = sorted(glob.glob(directory + '/a*'), key=lambda p: float(p.split('/a')[-1]))
    res = sorted((run(c) for c in cases), key=lambda r: r['alpha'])
    ny = res[0]['ny']
    return dict(T_wall_cell=foamio.read(cases[0] + '/0/T', ny=ny)['internal'][0].tolist(),
                x_mm=(1e3*(np.arange(NX) + .5)*dx).tolist(),
                y_mm=(1e3*(np.arange(ny) + .5)*4.8e-3/ny).tolist(), runs=res)


def g(v):
    return '' if v is None else f'{v:.5g}'


CSV = [('alpha', lambda r: r['alpha']), ('Da_561K', lambda r: r['Da_at_561K']),
       ('Da_400K', lambda r: r['Da_400_1000'][0]), ('Da_1000K', lambda r: r['Da_400_1000'][1]),
       ('deposited_fraction', lambda r: r['deposited_fraction']),
       ('outlet_fraction', lambda r: r['outlet_fraction']), ('gas_fraction_3s', lambda r: r['gas_fraction']),
       ('T_onset_K', lambda r: r['T_onset']), ('x_onset_mm', lambda r: r['x_onset_mm']),
       ('T_peak_bin_K', lambda r: r['T_peak']), ('x_peak_mm', lambda r: r['x_peak_mm']),
       ('centroid_mm', lambda r: r['centroid_mm']), ('sigma_mm', lambda r: r['sigma_mm']),
       ('width_10_90_mm', lambda r: r['width_10_90_mm']),
       ('wall_profile_L1_vs_alpha1', None),
       ('Jhat_72p5', lambda r: r['t1.5']['Jhat'][72]), ('Jhat_84p5', lambda r: r['t1.5']['Jhat'][84]),
       ('Jhat_94p5', lambda r: r['t1.5']['Jhat'][94]), ('eta_84p5', lambda r: r['t1.5']['eta'][84]),
       ('Jhat_zone', lambda r: r['t1.5']['zone_Jhat']), ('eta_zone', lambda r: r['t1.5']['zone_eta']),
       ('max_rel_closure', lambda r: r['max_relative_closure']), ('exec_s', lambda r: r['execution_s'])]


def write_csv(runs, path):
    ref = np.array(next(r for r in runs if r['alpha'] == 1.0)['wall_profile'])
    with open(path, 'w', newline='') as f:
        wr = csv.writer(f, lineterminator='\n')
        wr.writerow([k for k, _ in CSV])
        for r in runs:
            l1 = np.abs(np.array(r['wall_profile']) - ref).sum()/ref.sum()   # unnormalised, relative to alpha = 1
            wr.writerow([g(l1) if fn is None else g(fn(r)) for _, fn in CSV])


if __name__ == '__main__':
    base, refined = sys.argv[1], sys.argv[2]
    out = dict(sweep=sweep(base), refined_100x48=sweep(refined))
    json.dump(out, open(os.path.join(HERE, 'sweep_summary.json'), 'w'), indent=0)
    write_csv(out['sweep']['runs'], os.path.join(HERE, 'damkohler_sweep.csv'))
    print(f"{'alpha':>9} {'Da561':>9} {'f_dep':>7} {'T_on':>7} {'x_pk':>6} {'cen':>6} {'w1090':>6} {'J84':>7}")
    for r in out['sweep']['runs']:
        print(f"{r['alpha']:9.2e} {r['Da_at_561K']:9.3g} {r['deposited_fraction']:7.4f} {r['T_onset'] or 0:7.2f}"
              f" {r['x_peak_mm'] or 0:6.1f} {r['centroid_mm']:6.2f} {r['width_10_90_mm']:6.2f} {r['t1.5']['Jhat'][84]:7.3g}")
