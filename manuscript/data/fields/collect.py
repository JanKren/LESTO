#!/usr/bin/env python3
"""Gas-phase PbI2 field of the pure PbI2 steel column at 60 s: the tracked records of Fig. deposition_field.

Two runs of PbI2_SS_Q1_baseline, each checked against the tracked wall profile of the same run:

  production   576 x 72 (41,472 cells), Gauss limitedLinear 1; the figure and the manuscript numbers.
               Profiles: run/024-liu2025/production_576x72_limitedLinear/profiles/PbI2_SS_Q1_baseline
  upwind_L4    864 x 108 (93,312 cells), Gauss upwind; a scheme and mesh cross-check of the numbers.
               Profiles: run/024-liu2025/scheme_mesh_study/upwind_L4/profiles/PbI2_SS_Q1_baseline

    python3 manuscript/data/fields/collect.py production <case directory of the production run>
    python3 manuscript/data/fields/collect.py upwind_L4 \
        ~/lesto-submission-followup/cases/PbI2_SS_Q1_baseline__grid864x108

pbi2_field_<record>.npz   window x = 30-50 cm: PbI2(g) mole fraction at full radial resolution; whole
                          column: wall temperature, PbI2(g) mole fraction on the axis, in the wall row and
                          flow-weighted over the section, the equilibrium mole fraction at the wall
                          temperature, and the PbI2 wall inventory per unit length of the 5 degree wedge
provenance.json           per record: case, solver, mesh, scheme and the sha256 of every file read
"""
import hashlib, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '../../..'))
RECORDS = {
    'production': {'profile': 'run/024-liu2025/production_576x72_limitedLinear/profiles/PbI2_SS_Q1_baseline',
                   'solver_record': 'run/024-liu2025/production_576x72_limitedLinear/ledger.json',
                   'scheme': 'Gauss limitedLinear 1'},
    'upwind_L4': {'profile': 'run/024-liu2025/scheme_mesh_study/upwind_L4/profiles/PbI2_SS_Q1_baseline',
                  'solver_record': 'run/024-liu2025/scheme_mesh_study/upwind_L4/summary.json',
                  'scheme': 'Gauss upwind'},
}
L, R = 1.25, 2.4e-3                 # column length and inner radius [m]
THETA = np.deg2rad(5.0)             # wedge angle [rad]
WALL_WIDTH = 2 * R * np.tan(THETA / 2)   # wall face width of the blockMesh wedge [m] (chord, not R theta)
M_PBI2 = 0.46100894                 # [kg/mol]
Ru = 8.314462618                    # [J/(mol K)]
WINDOW_CM = (30.0, 50.0)


def internal(path, n, components=1):
    """internalField of an OpenFOAM v2412 volField, ascii or binary (LSB, 64-bit scalars)."""
    b = open(path, 'rb').read()
    i = b.index(b'internalField')
    m = re.compile(rb'(nonuniform List<(?:scalar|vector)>|uniform)\s*').search(b, i)
    if m.group(1) == b'uniform':
        return np.full(n, float(re.match(rb'([^;]+);', b[m.end():]).group(1)))
    k = int(re.match(rb'\s*(\d+)', b[m.end():]).group(1)); assert k == n, (path, k)
    p = b.index(b'(', m.end())
    if re.search(rb'format\s+binary', b[:i]):
        a = np.frombuffer(b[p + 1:p + 1 + 8 * n * components], '<f8').copy()
    else:
        a = np.array(b[p + 1:b.index(b'\n)', p)].replace(b'(', b' ').replace(b')', b' ').split(), float)
    return a.reshape(n, components) if components > 1 else a


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def main(name, case):
    spec = RECORDS[name]
    case = os.path.abspath(os.path.expanduser(case))
    block = open(os.path.join(case, 'system/blockMeshDict')).read()
    NX, NY = map(int, re.search(r'hex \([^)]*\) \((\d+) (\d+) 1\)', block).groups())
    assert '(1.25 0.0024' in block, 'unexpected geometry'
    schemes = open(os.path.join(case, 'system/fvSchemes')).read()
    assert f"div(phi,Y_PbI2_g) {spec['scheme']};" in schemes, 'unexpected scheme'
    read = []
    def field(fname, time='60', components=1):
        path = os.path.join(case, time, fname); read.append(path)
        a = internal(path, NX * NY, components)
        return a.reshape(NY, NX, components) if components > 1 else a.reshape(NY, NX)   # [j (radial), i (axial)]

    x = (np.arange(NX) + 0.5) * L / NX                       # cell centres [m]
    r = (np.arange(NY) + 0.5) * R / NY
    T, p = field('T', '0'), field('p', '0')                  # frozen carrier
    U = field('U', '0', 3)[..., 0]
    c = field('c_PbI2_g')                                    # [mol/m3]
    deposit = field('mDep_PbI2_s')[-1]                       # wall row [kg/m2]
    mole = c * Ru * T / p                                    # PbI2(g) mole fraction
    flow = U * r[:, None]                                    # u_x r: axisymmetric flow weight
    cup = (flow * mole).sum(0) / flow.sum(0)

    # equilibrium at the wall temperature from the case's own table (the minimum over the phases)
    props = open(os.path.join(case, 'constant/thermochemistryProperties')).read()
    table = re.search(r'PbI2_g \{.*?file "([^"]+)"', props).group(1)
    if not os.path.exists(table):                            # case 034 was renamed 024 (plan section 40)
        table = table.replace('run/034-liu2025', 'run/024-liu2025')
    read.append(table)
    pv = np.loadtxt(table, delimiter=',', comments='#')
    wall_T = T[-1]                                           # gas temperature = wall temperature per section
    x_eq = 10 ** np.interp(wall_T, pv[:, 0], pv[:, 1:].min(1)) / p[-1]

    w = (x >= WINDOW_CM[0] / 100) & (x <= WINDOW_CM[1] / 100)
    out = os.path.join(HERE, f'pbi2_field_{name}.npz')
    np.savez_compressed(out,
        x_cm=100 * x, r_mm=1e3 * r, window=np.flatnonzero(w),
        mole_window=mole[:, w],
        wall_T_K=wall_T, mole_axis=mole[0], mole_wall_row=mole[-1], mole_cup=cup, mole_eq_wall=x_eq,
        wall_mol_per_m=deposit / M_PBI2 * WALL_WIDTH)

    # the inventory against the tracked element-I profile of the same run
    inventory = float((deposit / M_PBI2 * WALL_WIDTH * L / NX).sum())
    tracked = os.path.join(REPO, spec['profile'], 'element_I_native.dat')
    wall_I = float(re.search(r'# inventory gas \S+ wall (\S+)', open(tracked).read()).group(1))
    solver = json.load(open(os.path.join(REPO, spec['solver_record'])))['solver_sha256']
    provenance_path = os.path.join(HERE, 'provenance.json')
    provenance = json.load(open(provenance_path)) if os.path.exists(provenance_path) else {}
    provenance = {k: v for k, v in provenance.items() if k in RECORDS}   # drop the single-record layout
    provenance[name] = {
        'case': os.path.basename(case),
        'case_directory_at_collection': case,
        'tracked_profile': spec['profile'],
        'solver_sha256': solver,
        'solver_sha256_from': spec['solver_record'],
        'mesh': f'{NX} x {NY} = {NX * NY} cells, 5 degree wedge, R = 2.4 mm, L = 1.25 m',
        'scheme': f"{spec['scheme']} (species advection); 2 ms, 60 s",
        'window_cm': list(WINDOW_CM),
        'mole_fraction': 'x = c R T / p with the frozen-carrier T and p of 0/',
        'cup_mean_weight': 'u_x r (axisymmetric flow-weighted mean over the section)',
        'equilibrium': 'p_eq = 10**min(phases) of the case table at the wall-row temperature, x_eq = p_eq / p',
        'wall_inventory_PbI2_mol': inventory,
        'tracked_profile_wall_I_mol': wall_I,
        'inventory_ratio_2PbI2_over_I': 2 * inventory / wall_I,
        'read_sha256': {os.path.relpath(f, case) if f.startswith(case) else os.path.relpath(f, REPO): sha(f)
                        for f in read},
        'npz_sha256': sha(out),
    }
    json.dump(dict(sorted(provenance.items())), open(provenance_path, 'w'), indent=1)
    print(name, json.dumps({k: provenance[name][k] for k in ['mesh', 'inventory_ratio_2PbI2_over_I']}))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
