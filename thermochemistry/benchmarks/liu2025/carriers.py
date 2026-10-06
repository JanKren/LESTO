#!/usr/bin/env python3
"""Prepare/run physical helium carriers using explicitly assumed flow conditions."""
import argparse
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

from prepare import HERE, R, HE_MOLAR_MASS, load_references, mesh_dictionary, read_csv

REPO = HERE.parents[2]


def header(name, cls='dictionary'):
    return f'FoamFile {{version 2.0; format ascii; class {cls}; object {name};}}\n'


def field(case, name, dims, value, inlet, outlet, wall, vector=False):
    cls = 'volVectorField' if vector else 'volScalarField'
    text = header(name, cls) + f'dimensions [{dims}];\ninternalField uniform {value};\nboundaryField {{\n'
    for patch, bc in [('INLET', inlet), ('OUTLET', outlet), ('WALL', wall),
                      ('BACK', 'type wedge;'), ('FRONT', 'type wedge;'), ('AXIS', 'type empty;')]:
        text += f'{patch} {{{bc}}}\n'
    (case / '0' / name).write_text(text + '}\n')


def prepare_carrier(case, name, nx, nr, flow_T, flow_p, temperature_root=None):
    if case.exists():
        raise ValueError(f'{case} exists; use a fresh work directory')
    experiments, _, _ = load_references()
    experiment = experiments[name]
    source = (temperature_root or HERE / 'digitized') / (name + '_temperature.csv')
    if not source.exists():
        raise ValueError('No published or supplied temperature profile for ' + name)
    rows = read_csv(source)
    points = [(float(r['x_m']), float(r['T_K'])) for r in rows]
    length, radius, angle = float(experiment['column_length_m']), float(experiment['inner_diameter_m'])/2, 5
    if len(points)<2 or any(not math.isfinite(x) or not math.isfinite(t) or t<=0
                            or not 0<=x<=length for x,t in points):
        raise ValueError('temperature profile needs at least two finite, positive temperatures within the column')
    if any(b[0]<=a[0] for a,b in zip(points,points[1:])):
        raise ValueError('temperature coordinates must increase strictly')
    # Explicit extrapolation assumption; original digitised coverage is retained.
    extended = ([(0., points[0][1])] if points[0][0]>0 else []) + points
    if points[-1][0]<length:
        extended += [(length, points[-1][1])]
    flow = float(experiment['flow_ml_min'])*1e-6/60
    mass_flow = flow*flow_p*HE_MOLAR_MASS/(R*flow_T)*angle/360
    rho_inlet = 101325*HE_MOLAR_MASS/(R*extended[0][1])
    umean = mass_flow / (rho_inlet*radius**2*math.tan(math.radians(angle/2)))
    for directory in ['0', 'constant', 'system']:
        (case / directory).mkdir(parents=True)
    (case/'constant/wallTemperature').write_text('(\n' + '\n'.join(f'({x:.17g} {t:.17g})' for x,t in extended) + '\n)\n')
    (case/'system/blockMeshDict').write_text(mesh_dictionary(length, radius, nx, nr, angle))
    (case/'system/liuCarrierDict').write_text(header('liuCarrierDict') +
        f'massFlowRate {mass_flow:.17g}; radius {radius}; pressure 101325;\n')
    (case/'system/controlDict').write_text(header('controlDict') + '''
application rhoSimpleFoam; startFrom startTime; startTime 0; stopAt endTime;
endTime 3000; deltaT 1; writeControl timeStep; writeInterval 500;
purgeWrite 1; writeFormat ascii; writePrecision 17; timePrecision 10;
runTimeModifiable yes;
''')
    (case/'system/fvSchemes').write_text(header('fvSchemes') + '''
ddtSchemes {default steadyState;} gradSchemes {default Gauss linear;}
divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,e) bounded Gauss upwind; div(phi,Ekp) bounded Gauss upwind;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;} snGradSchemes {default corrected;}
''')
    (case/'system/fvSolution').write_text(header('fvSolution') + '''
solvers {
 p {solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.0001; maxIter 2000;}
 "(U|e)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0.01;}
}
SIMPLE {nNonOrthogonalCorrectors 0;
 residualControl {p 1e-10; U 1e-10; e 1e-10;}}
relaxationFactors {fields {p 0.3;} equations {U 0.5; e 0.7;}}
''')
    shutil.copyfile(REPO/'run/007-variable-properties-high-temperature/constant/thermophysicalProperties',
                    case/'constant/thermophysicalProperties')
    (case/'constant/turbulenceProperties').write_text(header('turbulenceProperties') + 'simulationType laminar;\n')
    field(case, 'p', '1 -1 -2 0 0 0 0', '101325', 'type zeroGradient;',
          'type fixedValue; value uniform 101325;', 'type zeroGradient;')
    field(case, 'T', '0 0 0 1 0 0 0', str(extended[0][1]),
          f'type fixedValue; value uniform {extended[0][1]};', 'type zeroGradient;',
          f'type fixedValue; value uniform {extended[0][1]};')
    field(case, 'U', '0 1 -1 0 0 0 0', f'({umean} 0 0)',
          f'type flowRateInletVelocity; massFlowRate {mass_flow:.17g}; rho rho; rhoInlet {rho_inlet:.17g}; '
          f'extrapolateProfile yes; value uniform ({umean} 0 0);',
          'type zeroGradient;', 'type noSlip;', vector=True)
    manifest = dict(experiment=name, status='PROVISIONAL_INPUTS', cells=nx*nr,
        flow_reference_temperature_K=flow_T, flow_reference_pressure_Pa=flow_p,
        flow_reference_status='ASSUMED, not confirmed by authors',
        full_column_mass_flow_kg_s=mass_flow*360/angle, wedge_mass_flow_kg_s=mass_flow,
        wedge_angle_degrees=angle, wedge_fraction=angle/360, radius_m=radius,
        temperature_csv_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        temperature_source=str(source.resolve()),
        temperature_coverage_m=[points[0][0], points[-1][0]],
        assumptions=['Plot x=0 is mesh x=0; coordinate origin unconfirmed.',
                     'Constant temperature extension at both unplotted ends.',
                     'Helium hConst/Sutherland properties from case 007; ideal gas at 101325 Pa.',
                     'Empty tube: boat obstruction omitted.'],
        scientific_validation_status='NOTRUN')
    (case/'benchmark.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def call(case, command, label):
    with (case / ('log.' + label)).open('w') as stream:
        subprocess.run(command + ['-case', str(case)], stdout=stream, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--nx', type=int, default=384)
    parser.add_argument('--nr', type=int, default=48)
    parser.add_argument('--flow-reference-temperature', type=float, required=True)
    parser.add_argument('--flow-reference-pressure', type=float, required=True)
    parser.add_argument('--experiment', choices=list(load_references()[0]), action='append')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--tool', default='liuCarrierTool')
    parser.add_argument('--temperature-root',type=Path,help='Override figure CSVs with supplied profiles; same filenames and x_m,T_K columns')
    args = parser.parse_args()
    for value in [args.nx, args.nr, args.flow_reference_temperature, args.flow_reference_pressure]:
        if not math.isfinite(value) or value <= 0:
            parser.error('mesh sizes and flow reference conditions must be positive and finite')
    for name in args.experiment or ['PbI2_SS', 'LBE-I_SS_II', 'LBE-I_SiO2_I']:
        case = (args.out / name).resolve()
        prepare_carrier(case, name, args.nx, args.nr, args.flow_reference_temperature, args.flow_reference_pressure,args.temperature_root)
        if args.run:
            for command, label in [(['blockMesh'], 'blockMesh'), (['checkMesh', '-constant'], 'checkMesh'),
                                   ([args.tool], 'initialise'), (['rhoSimpleFoam'], 'carrier')]:
                call(case, command, label)
            times = sorted((p for p in case.iterdir() if p.is_dir() and p.name.replace('.', '', 1).isdigit()),
                           key=lambda p:float(p.name))
            latest = times[-1].name
            call(case, [args.tool, '-project', '-time', latest], 'project')
            call(case, [args.tool, '-audit', '-time', latest], 'audit')
            print(name, (case/'carrier-audit.json').read_text(), flush=True)


if __name__ == '__main__':
    main()
