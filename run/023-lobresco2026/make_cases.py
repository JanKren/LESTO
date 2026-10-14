#!/usr/bin/env python3
"""Lobresco et al. (2026) five PbI2 cases: physical He carriers + HKS transport cases.

Only existing solver features are used (rhoSimpleFoam carrier as run/024; HKS pair,
inventory sample, halfCell wall law and profiles of rhoFixedFlowFoam).  Nothing is
written into the repository; repository files are copied, never linked.
"""
import argparse, csv, hashlib, json, math, os, shutil, subprocess
from pathlib import Path

G = Path(__file__).resolve().parent
REPO = G.parents[1]
CPUS = os.environ.get('LESTO_CPUS', '0-3')   # taskset core list for carrier builds
R = 8.31446261815324
M_HE = 0.0040026
M_PBI2 = 0.46100894
RADIUS = 2.4e-3          # m, 4.8 mm ID (paper Sec. 2.1)
LENGTH = 0.69            # m, x range of the published wall-T curves and T-Flows mesh
ANGLE = 5.0              # degrees, wedge as run/024
FRACTION = ANGLE / 360.0
P_REF, T_REF = 101325.0, 273.15   # paper Sec. 2.1: standard flow at 273.15 K, 1.01325 bar
LEVELS = {'prod': (318, 72), 'coarse': (212, 48)}   # same dx, dr as run/024 576x72 / 384x48 on 1.25 m
# Table 1 (verbatim) + arrow positions digitised from Figs. 13-14
CASES = {
  'T106_6':   dict(flow=106, duration=6,   mass_mg=2.23, moles=4.84e-6, arrow=0.130, dt=0.002, write=1.0),
  'T106_60':  dict(flow=106, duration=60,  mass_mg=2.04, moles=4.43e-6, arrow=0.130, dt=0.002, write=10.0),
  'T106_300': dict(flow=106, duration=300, mass_mg=1.37, moles=2.97e-6, arrow=0.130, dt=0.002, write=30.0),
  'T5_60':    dict(flow=5,   duration=60,  mass_mg=1.45, moles=3.15e-6, arrow=0.150, dt=0.002, write=10.0),
  'T540_60':  dict(flow=540, duration=60,  mass_mg=1.81, moles=3.93e-6, arrow=0.150, dt=0.0005, write=10.0),
}
BOAT_LENGTH = 0.01       # m, as the T-Flows sample (0.04-0.05 m); assumption
BOAT_RADIUS = 4.8e-4     # m, as run/013 (T-Flows R_SAMPLE 0.5 mm, between cell rings)
BOAT_AV = 4000.0         # 1/m, 2/R_sample as run/013


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def header(name, cls='dictionary'):
    return f'FoamFile {{version 2.0; format ascii; class {cls}; object {name};}}\n'


def field(case, name, dims, value, inlet, outlet, wall, vector=False, where='0'):
    cls = 'volVectorField' if vector else 'volScalarField'
    text = header(name, cls) + f'dimensions [{dims}];\ninternalField uniform {value};\nboundaryField {{\n'
    for patch, bc in [('INLET', inlet), ('OUTLET', outlet), ('WALL', wall),
                      ('BACK', 'type wedge;'), ('FRONT', 'type wedge;'), ('AXIS', 'type empty;')]:
        text += f'{patch} {{{bc}}}\n'
    (case / where / name).write_text(text + '}\n')


def mesh_dictionary(nx, nr):
    rt = RADIUS * math.tan(math.radians(ANGLE / 2))
    return header('blockMeshDict') + f"""mergeType points;
scale 1;
vertices (
  (0 0 0) ({LENGTH} 0 0) ({LENGTH} {RADIUS} {-rt}) (0 {RADIUS} {-rt})
  ({LENGTH} {RADIUS} {rt}) (0 {RADIUS} {rt})
);
blocks (hex (0 1 2 3 0 1 4 5) ({nx} {nr} 1) simpleGrading (1 1 1));
boundary (
  INLET {{ type patch; faces ((0 0 5 3)); }}
  OUTLET {{ type patch; faces ((1 2 4 1)); }}
  WALL {{ type wall; faces ((3 5 4 2)); }}
  BACK {{ type wedge; faces ((0 3 2 1)); }}
  FRONT {{ type wedge; faces ((0 1 4 5)); }}
  AXIS {{ type empty; faces ((0 1 1 0)); }}
);
"""


def prepare_carrier(case, flow, nx, nr):
    if case.exists():
        raise ValueError(f'{case} exists')
    src = G / 'inputs' / f'wallT_{flow}mLmin.csv'
    pts = [(float(r['x_m']), float(r['T_K'])) for r in csv.DictReader(open(src))]
    q = flow * 1e-6 / 60
    mass_flow = q * P_REF * M_HE / (R * T_REF) * FRACTION
    rho_in = 101325 * M_HE / (R * pts[0][1])
    umean = mass_flow / (rho_in * RADIUS**2 * math.tan(math.radians(ANGLE / 2)))
    for d in ['0', 'constant', 'system']:
        (case / d).mkdir(parents=True)
    (case / 'constant/wallTemperature').write_text('(\n' + '\n'.join(f'({x:.17g} {t:.17g})' for x, t in pts) + '\n)\n')
    (case / 'system/blockMeshDict').write_text(mesh_dictionary(nx, nr))
    (case / 'system/liuCarrierDict').write_text(header('liuCarrierDict') +
        f'massFlowRate {mass_flow:.17g}; radius {RADIUS}; pressure 101325;\n')
    (case / 'system/controlDict').write_text(header('controlDict') + '''
application rhoSimpleFoam; startFrom startTime; startTime 0; stopAt endTime;
endTime 3000; deltaT 1; writeControl timeStep; writeInterval 500;
purgeWrite 1; writeFormat ascii; writePrecision 17; timePrecision 10;
runTimeModifiable yes;
''')
    # identical to thermochemistry/benchmarks/liu2025/carriers.py (run/024)
    (case / 'system/fvSchemes').write_text(header('fvSchemes') + '''
ddtSchemes {default steadyState;} gradSchemes {default Gauss linear;}
divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,e) bounded Gauss upwind; div(phi,Ekp) bounded Gauss upwind;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;} snGradSchemes {default corrected;}
''')
    (case / 'system/fvSolution').write_text(header('fvSolution') + '''
solvers {
 p {solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.0001; maxIter 2000;}
 "(U|e)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0.01;}
}
SIMPLE {nNonOrthogonalCorrectors 0;
 residualControl {p 1e-10; U 1e-10; e 1e-10;}}
relaxationFactors {fields {p 0.3;} equations {U 0.5; e 0.7;}}
''')
    shutil.copyfile(REPO / 'run/007-variable-properties-high-temperature/constant/thermophysicalProperties',
                    case / 'constant/thermophysicalProperties')
    (case / 'constant/turbulenceProperties').write_text(header('turbulenceProperties') + 'simulationType laminar;\n')
    field(case, 'p', '1 -1 -2 0 0 0 0', '101325', 'type zeroGradient;',
          'type fixedValue; value uniform 101325;', 'type zeroGradient;')
    field(case, 'T', '0 0 0 1 0 0 0', str(pts[0][1]),
          f'type fixedValue; value uniform {pts[0][1]};', 'type zeroGradient;',
          f'type fixedValue; value uniform {pts[0][1]};')
    field(case, 'U', '0 1 -1 0 0 0 0', f'({umean} 0 0)',
          f'type flowRateInletVelocity; massFlowRate {mass_flow:.17g}; rho rho; rhoInlet {rho_in:.17g}; '
          f'extrapolateProfile yes; value uniform ({umean} 0 0);',
          'type zeroGradient;', 'type noSlip;', vector=True)
    man = dict(flow_mL_min_standard=flow, flow_reference_K=T_REF, flow_reference_Pa=P_REF,
               flow_reference_source='Lobresco et al. (2026) Sec. 2.1: standard volumetric flow at 273.15 K and 1.01325 bar',
               full_column_mass_flow_kg_s=mass_flow / FRACTION, wedge_mass_flow_kg_s=mass_flow,
               wedge_angle_deg=ANGLE, radius_m=RADIUS, length_m=LENGTH, cells=nx * nr, nx=nx, nr=nr,
               wall_T_csv=str(src.relative_to(G)), wall_T_sha256=sha(src),
               thermophysical='He, hConst Cp 5193, Sutherland (run/007), perfect gas at 101325 Pa',
               assumptions=['Figure x = 0 is the inlet (clean He, wall T of the curve at x=0).',
                            'Empty tube: the boat obstruction is omitted.',
                            'He properties of run/007, not the NIST tables of the paper.'])
    (case / 'carrier.json').write_text(json.dumps(man, indent=1) + '\n')
    return man


def run(case, cmd, label, env_bin):
    with (case / ('log.' + label)).open('w') as s:
        subprocess.run(['taskset', '-c', CPUS] + cmd + ['-case', str(case)], stdout=s,
                       stderr=subprocess.STDOUT, check=True)


def build_carrier(case):
    tool = str(G / 'bin/liuCarrierTool')
    for cmd, label in [(['blockMesh'], 'blockMesh'), (['checkMesh', '-constant'], 'checkMesh'),
                       ([tool], 'initialise'), (['rhoSimpleFoam'], 'carrier')]:
        run(case, cmd, label, None)
    times = sorted((p for p in case.iterdir() if p.is_dir() and p.name.replace('.', '', 1).isdigit() and p.name != '0'),
                   key=lambda p: float(p.name))
    latest = times[-1].name
    run(case, [tool, '-project', '-time', latest], 'project', None)
    run(case, [tool, '-audit', '-time', latest], 'audit', None)
    return json.loads((case / 'carrier-audit.json').read_text())


INLET_DANCKWERTS = '''type codedMixed; name zeroTotalFluxInlet;
 refValue uniform 0; refGradient uniform 0; valueFraction uniform 1; value uniform 0;
 code
 #{
   const fvsPatchScalarField& phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
   const fvPatchScalarField& rhoDp = patch().lookupPatchField<volScalarField, scalar>("rhoD_PbI2_g");
   const scalarField Pe(max(-phip, scalar(0))/(rhoDp*patch().magSf()*patch().deltaCoeffs()));
   valueFraction() = Pe/(1 + Pe);
 #};'''


def prepare_case(case, name, carrier, nx, end_time=None, inlet='danckwerts'):
    c = CASES[name]
    if case.exists():
        raise ValueError(f'{case} exists')
    audit = json.loads((carrier / 'carrier-audit.json').read_text())
    if audit['status'] != 'PASS':
        raise ValueError('carrier audit failed')
    case.mkdir(parents=True)
    for d in ['constant', 'system']:
        shutil.copytree(carrier / d, case / d)
    (case / '0').mkdir()
    t = format(audit['time'], 'g')
    for f in ['T', 'p', 'rho', 'U', 'phi']:
        shutil.copyfile(carrier / t / f, case / '0' / f)
    inlet_bc = INLET_DANCKWERTS if inlet == 'danckwerts' else 'type fixedValue; value uniform 0;'
    field(case, 'Y_PbI2_g', '0 0 0 0 0 0 0', '0', inlet_bc, 'type zeroGradient;', 'type zeroGradient;')
    field(case, 'Y_PbI2_s', '0 0 0 0 0 0 0', '0', 'type zeroGradient;', 'type zeroGradient;', 'type zeroGradient;')
    run013 = REPO / 'run/013-hks-inventory'
    shutil.copyfile(run013 / 'constant/speciesTransportProperties', case / 'constant/speciesTransportProperties')
    pv = REPO / 'thermochemistry/systems/data/pv_PbI2_phases.csv'
    shutil.copyfile(pv, case / 'constant/pv_PbI2_phases.csv')
    for f in ['fvSolution']:
        shutil.copyfile(run013 / 'system' / f, case / 'system' / f)
    (case / 'system/fvSchemes').write_text(header('fvSchemes') + '''
ddtSchemes {default Euler;}
gradSchemes {default Gauss linear;}
divSchemes {default none; div(phi,Y_PbI2_g) Gauss limitedLinear 1;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;}
snGradSchemes {default corrected;}
''')
    for f in ['controlDict', 'fvSchemes', 'fvSolution', 'blockMeshDict', 'liuCarrierDict']:
        pass
    (case / 'system/liuCarrierDict').unlink()
    end = c['duration'] if end_time is None else end_time
    dt = c['dt']
    write_steps = max(1, round(min(c['write'], end) / dt))
    (case / 'system/controlDict').write_text(header('controlDict') + f'''
application rhoFixedFlowFoam; startFrom latestTime; startTime 0;
stopAt endTime; endTime {end:.17g}; deltaT {dt:.17g};
writeControl timeStep; writeInterval {write_steps}; purgeWrite 0;
writeFormat binary; writeCompression off; writePrecision 17; timePrecision 10;
runTimeModifiable false;
speciesTransport {{tolerance 1e-5; nCorr 10; writeFrozenFields no;}}
''')
    n_wedge = c['moles'] * FRACTION
    xa = c['arrow']
    times = sorted({round(k * c['write'], 9) for k in range(1, int(round(c['duration'] / c['write'])) + 1)} | {c['duration']})
    (case / 'constant/thermochemistryProperties').write_text(header('thermochemistryProperties') + f'''
model HKS;
pairs {{
  PbI2_g {{
    condensed PbI2_s;
    vapourPressure {{
      file "<constant>/pv_PbI2_phases.csv"; phases ("PbI2(cr)" "PbI2(l)");
      units log10Pa; interpolation logInverseT; outOfRange extrapolate;
    }}
  }}
}}
interface {{
  patches (WALL); layer faceCells; areaPerVolume cell; areaMultiplier 1;
  interfaceTemperature wall; wallResistance halfCell;
}}
HKSCoeffs {{ accommodation 1; Ce 1; kineticScale 1; equilibrium table; speciation none; }}
samples {{
  boat {{
    pair PbI2_g;
    selection {{ boat {{ action use; source cylinder; point1 ({xa - BOAT_LENGTH/2:.6g} 0 0);
                       point2 ({xa + BOAT_LENGTH/2:.6g} 0 0); radius {BOAT_RADIUS}; }} }}
    mode inventory;
    amount {n_wedge:.17g};   /* mol in the 5-degree wedge = Table 1 moles x 5/360 */
    areaPerVolume {BOAT_AV};
    kineticScale 1; Ce 1;
    removeTime {c['duration']:.17g};
  }}
}}
balance {{ nGuard 3; guardTolerance 1e-12; moleFractionWarning 0.05; writeFile yes; }}
profiles {{
  bins {{
    native {{ min 0; max {LENGTH}; nBins {nx}; }}
    gamma  {{ min -0.005; max 0.695; nBins 70; }}   /* 1 cm bins centred on the scan points */
  }}
  times ({' '.join(f'{x:g}' for x in times)});
  onset {{ fraction 0.01; search (fromInlet fromPeak); }}
  observable wallPlusGas;   /* shutdown comparison (plan section 28, item 3) */
}}
''')
    man = dict(case=name, table1=dict(flow_mL_min=c['flow'], duration_s=c['duration'], mass_mg=c['mass_mg'],
               moles=c['moles']), wedge_moles=n_wedge, arrow_x_m=xa,
               boat=dict(x_m=[xa - BOAT_LENGTH/2, xa + BOAT_LENGTH/2], radius_m=BOAT_RADIUS, areaPerVolume=BOAT_AV),
               dt_s=dt, end_time_s=end, inlet_bc=inlet, scheme='Gauss limitedLinear 1',
               carrier=json.loads((carrier / 'carrier.json').read_text()), carrier_audit=audit,
               pv_table_sha256=sha(pv), solver_sha256=sha(G / 'bin/rhoFixedFlowFoam'))
    (case / 'case.json').write_text(json.dumps(man, indent=1) + '\n')
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=G / 'cases')
    ap.add_argument('--levels', nargs='+', default=['prod', 'coarse'])
    ap.add_argument('--carriers-only', action='store_true')
    ap.add_argument('--flows', nargs='+', type=int, default=[5, 106, 540])
    ap.add_argument('--cases', nargs='+', default=list(CASES))
    args = ap.parse_args()
    for lev in args.levels:
        nx, nr = LEVELS[lev]
        for flow in args.flows:
            car = args.out / lev / 'carriers' / f'He{flow}'
            if not (car / 'carrier-audit.json').exists():
                prepare_carrier(car, flow, nx, nr)
                a = build_carrier(car)
                print(lev, flow, a['status'], a['max_cell_flux_defect_over_inflow'], flush=True)
        if args.carriers_only:
            continue
        for name in args.cases:
            c = CASES[name]
            prepare_case(args.out / lev / name, name, args.out / lev / 'carriers' / f'He{c["flow"]}', nx)
            print('prepared', lev, name, flush=True)


if __name__ == '__main__':
    main()
