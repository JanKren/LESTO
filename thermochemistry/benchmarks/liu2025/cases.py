#!/usr/bin/env python3
"""M10e provisional Q1/Q2 cases on audited physical helium carriers.

Thermodynamic tables stay outside git. Sources retain the experimental mean
rate; shortening a CFD window never concentrates the three-hour inventory.
"""
import argparse
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

from carriers import REPO, field, header
from prepare import R, HE_MOLAR_MASS, load_references

GASES = ['Pb_g', 'PbI_g', 'PbI2_g', 'Bi_g', 'Bi2_g', 'BiI_g', 'BiI3_g', 'I_g', 'I2_g']
PRIMARY = ['PbI2_SS', 'BiI3_SiO2_800', 'LBE-I_SS_II', 'LBE-I_SiO2_I']


def positive(value, name):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + ' must be positive and finite')
    return value


def source_rates(experiment, targets, species, fraction):
    return {r['species']+'_g': float(r['source_species_mass_mg'])*1e-6
            /float(species[r['species']]['molar_mass_kg_mol'])/float(experiment['duration_s'])*fraction
            for r in targets if r['experiment'] == experiment['experiment'] and r['scope'] == 'M10'}


def q2_rates(q1, helium_rate, temperature, saturation, data):
    # Ideal LBE activity is an explicit sensitivity assumption, not a measured
    # release law. Atomic Pb/Bi/I sources re-speciate in the homogeneous kernel.
    pressure = 101325
    pPb = .447*1e5*math.exp(-(data.G('Pb(g)',temperature)-min(data.G('Pb(cr)',temperature),data.G('Pb(l)',temperature)))/(R*temperature))
    bi = min(data.G('Bi(cr)',temperature),data.G('Bi(l)',temperature))
    pBi = .553*1e5*math.exp(-(data.G('Bi(g)',temperature)-bi)/(R*temperature))
    pBi2 = .553**2*1e5*math.exp(-(data.G('Bi2(g)',temperature)-2*bi)/(R*temperature))
    return {'Pb_g': q1.get('PbI2_g',0)+helium_rate*saturation*pPb/pressure,
            'Bi_g': q1.get('BiI3_g',0)+helium_rate*saturation*(pBi+2*pBi2)/pressure,
            'I_g': 2*q1.get('PbI2_g',0)+3*q1.get('BiI3_g',0)}


def prepare_case(case, carrier, tables, mode, saturation, dt, end_time, diffusivity, data=None):
    if case.exists():
        raise ValueError(f'{case} exists; use a fresh output directory')
    metadata = json.loads((carrier/'benchmark.json').read_text())
    audit = json.loads((carrier/'carrier-audit.json').read_text())
    if audit['status'] != 'PASS':
        raise ValueError('Refusing an unaudited or failed carrier: '+str(carrier))
    experiment_name = metadata['experiment']
    experiments, targets, species = load_references()
    experiment = experiments[experiment_name]
    q1 = source_rates(experiment, targets, species, metadata['wedge_fraction'])
    rates = q1
    temperature = None
    if mode == 'Q2':
        if data is None:
            raise ValueError('Q2 needs external thermodynamic evaluator')
        temperature = float(experiment['hot_zone_temperature_C'])+273.15
        rates = q2_rates(q1, metadata['wedge_mass_flow_kg_s']/HE_MOLAR_MASS,
                         temperature, saturation, data)
    case.mkdir(parents=True)
    for folder in ['constant','system']:
        shutil.copytree(carrier/folder, case/folder)
    (case/'0').mkdir()
    time = format(audit['time'], 'g')
    for name in ['T','p','rho','U','phi']:
        shutil.copyfile(carrier/time/name, case/'0'/name)
    # The transport solver freezes the carrier; standard inlet/outlet/wall
    # species BCs do not inject an additional source or impose a wall sink.
    for name in GASES+['PbI2_s','BiI3_s','Pb_s','Bi_s']:
        field(case, 'Y_'+name, '0 0 0 0 0 0 0', '0',
              'type fixedValue; value uniform 0;', 'type zeroGradient;', 'type zeroGradient;')
    template = REPO/'run/033-lbe-metal-vapour'
    props = (template/'constant/speciesTransportProperties').read_text()
    props = props.replace('D 0;', f'D {diffusivity:.17g};')
    (case/'constant/speciesTransportProperties').write_text(props)
    shutil.copyfile(template/'system/fvSchemes',case/'system/fvSchemes')
    shutil.copyfile(template/'system/fvSolution',case/'system/fvSolution')
    steps = math.ceil(end_time/dt)
    write_steps = max(1, round(10/dt)) if end_time >= 20 else steps
    (case/'system/controlDict').write_text(header('controlDict')+f'''
application rhoFixedFlowFoam; startFrom latestTime; startTime 0;
stopAt endTime; endTime {end_time:.17g}; deltaT {dt:.17g};
writeControl timeStep; writeInterval {write_steps}; purgeWrite 0;
writeFormat binary; writeCompression off; writePrecision 17; timePrecision 10;
runTimeModifiable false;
speciesTransport {{tolerance 1e-16; nCorr 3; writeFrozenFields no;}}
''')
    thermo = (template/'constant/thermochemistryProperties').read_text()
    thermo = thermo.replace('$LESTO_M10D_DATA',str(tables.resolve()))
    thermo = thermo.replace('interface {patches (WALL);}',
                            'interface {patches (WALL); interfaceTemperature wall; wallResistance halfCell;}')
    thermo = thermo.replace('min 0; max .1; nBins 200;', 'min -.005; max 1.255; nBins 126;')
    thermo = thermo[:thermo.index('\ngasSources')]
    radius = metadata['radius_m']
    box = f'(.09 {radius*.25:.17g} -1) (.105 {radius*.75:.17g} 1)'
    amounts = ' '.join(f'{name} {rate*end_time:.17g};' for name,rate in rates.items())
    if mode == 'Q1':
        thermo = thermo.replace('samples {}',f'''samples {{boat {{amounts {{{amounts}}}
 selection {{boat {{action use; source box; box {box};}}}}
 mode release; startTime 0; duration {end_time:.17g};}}}}''')
    else:
        thermo += f'\ngasSources {{boat {{amounts {{{amounts}}} box ({box}); startTime 0; duration {end_time:.17g};}}}}\n'
    (case/'constant/thermochemistryProperties').write_text(thermo)
    included = [r for r in targets if r['experiment']==experiment_name and r['scope']=='M10']
    excluded = sum(float(r['iodine_yield_pct']) for r in targets
                   if r['experiment']==experiment_name and r['scope']!='M10' and r['iodine_yield_pct'])
    manifest = dict(experiment=experiment_name, status='PREPARED_PROVISIONAL', source_mode=mode,
        table_directory=str(tables.resolve()), table_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in tables.glob('*.csv')}, carrier=metadata, carrier_audit=audit,
        carrier_projection=json.loads((carrier/'carrier-projection.json').read_text()),
        source_wedge_mol_s=rates, source_full_column_mol_s={n:r/metadata['wedge_fraction'] for n,r in rates.items()},
        source_temperature_K=temperature, metal_saturation_fraction=saturation if mode=='Q2' else None,
        source_box_m=[[.09,radius*.25,-1],[.105,radius*.75,1]], excluded_KI_iodine_yield_pct=excluded,
        dt_s=dt, end_time_s=end_time, gas_diffusivity_m2_s=diffusivity,
        targets=included, scientific_validation_status='NOTRUN',
        assumptions=['Source at 9–10.5 cm in interior radial cells; exact boat position unconfirmed.',
                     'Mean three-hour release rate, no terminal clean-helium flush.',
                     'Constant gas diffusivity is provisional; transport properties need qualification.',
                     'Q1 LBE sources are conditioned on deposited masses and cannot validate split prediction.',
                     'Q2 uses ideal LBE activities; iodine source conditioned on deposits, release history unknown.',
                     'BiI +29 kJ/mol is a sensitivity surrogate, not an independently qualified dataset.'])
    (case/'benchmark.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--carrier-root',type=Path,required=True)
    parser.add_argument('--tables-root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--nasa',type=Path,required=True)
    parser.add_argument('--tf',type=Path,required=True)
    parser.add_argument('--dt',type=float,default=.002)
    parser.add_argument('--end-time',type=float,default=60)
    parser.add_argument('--diffusivity',type=float,default=1e-4)
    args = parser.parse_args()
    for name,value in [('dt',args.dt),('end time',args.end_time),('diffusivity',args.diffusivity)]:
        positive(value,name)
    sys.path.insert(0,str(REPO/'thermochemistry/systems'))
    from make_wall_tables import generate
    from make_pv_tables import Data
    data=Data(str(args.nasa),str(args.tf))
    for label,shift in [('baseline',0),('shift29k',29000)]:
        generate(str(args.nasa),str(args.tf),args.tables_root/label,shift)
    matrix=[]
    for experiment in PRIMARY:
        variants=[('Q1','baseline',None)]
        if experiment.startswith('LBE'):
            variants += [('Q2',label,sat) for label in ['baseline','shift29k'] for sat in [.001,.01,.1]]
        for mode,label,saturation in variants:
            name=experiment+'_'+mode+'_'+label+(f'_sat{saturation:g}' if saturation else '')
            carrier=args.carrier_root/experiment
            entry=dict(case=name, experiment=experiment, source_mode=mode, data_variant=label,
                       saturation=saturation, status='NOTRUN', reason='')
            if not (carrier/'carrier-audit.json').exists():
                entry['reason']='Missing measured temperature profile or audited physical carrier'
            else:
                prepare_case(args.out/name,carrier,args.tables_root/label,mode,saturation,args.dt,
                             args.end_time,args.diffusivity,data)
                entry['status']='PREPARED_PROVISIONAL'
            matrix.append(entry)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'matrix.json').write_text(json.dumps(matrix,indent=2)+'\n')
    print(f'Prepared {sum(r["status"]=="PREPARED_PROVISIONAL" for r in matrix)} of {len(matrix)} cases; missing inputs remain NOTRUN')


if __name__ == '__main__':
    main()
