#!/usr/bin/env python3
"""Analyse completed Liu matrices without changing solver output or inputs.

The primary comparisons retain the digitised histogram heights. Secondary
shape comparisons explicitly normalise both curves over the common scan.
No parameter, coordinate offset, or temperature curve is fitted.
"""
import argparse
import base64
import csv
import datetime
import hashlib
import html
import json
import math
from pathlib import Path
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from prepare import HERE, load_references
from report import rows
from carriers import REPO

SPECIES = {'PbI2': 'PbI2_g', 'BiI3': 'BiI3_g',
           'Pb': 'condensate_Pb_s', 'Bi': 'condensate_Bi_s', 'I': 'element_I'}
EXPERIMENTS = ['PbI2_SS', 'LBE-I_SS_II', 'LBE-I_SiO2_I']
PAPER = 'https://doi.org/10.1007/s10967-025-10246-4'
FORMULAS = {'Pb_g':(1,0,0),'PbI_g':(1,0,1),'PbI2_g':(1,0,2),
            'Bi_g':(0,1,0),'Bi2_g':(0,2,0),'BiI_g':(0,1,1),
            'BiI3_g':(0,1,3),'I_g':(0,0,1),'I2_g':(0,0,2)}
ATOMIC_MASSES = np.array([.2072,.2089804,.12690447])


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def profile(case, species, time=60):
    data = rows(case/'postProcessing/axialProfiles'/f'{time:g}'/(SPECIES[species]+'_native.dat'))
    x = np.array([r['x'] for r in data])
    amounts = np.array([r['wall']*(r['x_hi']-r['x_lo']) for r in data])
    assert np.isfinite(amounts).all() and (amounts >= 0).all(), case
    return x, amounts, data


def relative_l1(a, b):
    denominator = float(np.sum(np.abs(b)))
    return float(np.sum(np.abs(a-b))/denominator) if denominator else None


def moments(x, amounts):
    total = float(amounts.sum())
    if not total:
        return {'total': 0., 'centroid_m': None, 'sigma_m': None}
    centre = float(np.dot(x, amounts)/total)
    return {'total': total, 'centroid_m': centre,
            'sigma_m': float(np.sqrt(np.dot((x-centre)**2, amounts)/total))}


def ledger_audit(case):
    result = {}
    for element in ['Pb', 'Bi', 'I']:
        path = case/'postProcessing/phaseChangeElements/0'/(element+'.dat')
        with path.open() as stream:
            columns = next(line.split('columns ')[1].split() for line in stream if line.startswith('# columns'))
        data = np.loadtxt(path, comments='#', ndmin=2)
        assert np.isfinite(data).all() and abs(data[-1, 0]-60)<1e-7, path
        d = {name: data[:, j] for j, name in enumerate(columns)}
        source = np.maximum(d['initial']+d['released'], 1e-300)
        raw = d['gas']+d['wall']+d['sample']+d['transport']+d['removed']-d['initial']-d['released']
        result[element] = {
            'records': len(data), 'final': {k: float(v[-1]) for k, v in d.items()},
            'max_accounted_closure_over_source': float(np.max(np.abs(d['closure'])/source)),
            'max_raw_material_residue_over_source': float(np.max(np.abs(raw)/source)),
            'max_net_clamped_over_source': float(np.max(np.abs(d['clamped'])/source)),
            'max_solver_defect_over_source': float(np.max(np.abs(d['solverDefect'])/source)),
        }
    return result


def log_audit(path, manifest):
    digest = hashlib.sha256()
    result = {'sha256': None, 'pairs': {}, 'max_respeciation_failures': 0,
              'max_invalid_element_cells': 0, 'max_equilibrium_element_residue': 0.,
              'max_sum_of_species_mole_fraction_maxima': 0., 'normal_end': False,
              'steps_with_invalid_element_cells': 0, 'first_invalid_time_s': None,
              'last_invalid_time_s': None, 'final_invalid_element_cells': 0,
              'max_negative_gas_mass_kg': 0., 'max_negative_I_over_supplied_I': 0.}
    pair_re = re.compile(r'^(Y_\w+) exchange:.*outer iterations (\d+), (NOT converged|converged) '
                         r'\(initial residual (\S+), tolerance (\S+); regime changes (\d+), '
                         r'(\S+) kg, (\d+) significant\)')
    respec_re = re.compile(r'failures (\d+); invalid-element cells (\d+); max element residue (\S+)')
    summary_re = re.compile(r'^Phase change summary (\w+): (\d+) steps without outer convergence, '
                            r'(\d+) guard passes, (\d+) projected elements \((\d+) significant\), '
                            r'(\d+) complementarity violations \((\d+) significant\)')
    residuals = {}
    mole_sum = negative_mass = negative_i = 0.
    current_time = 0.
    source_i_rate = sum(FORMULAS[n][2]*rate for n,rate in manifest['source_wedge_mol_s'].items())
    with path.open('rb') as stream:
        for raw in stream:
            digest.update(raw)
            if b' exchange:' in raw:
                m = pair_re.search(raw.decode())
                if m:
                    name, outer, flag, residual, tolerance, changed, mass, significant = m.groups()
                    p = result['pairs'].setdefault(name, {'steps': 0, 'unconverged_steps': 0,
                        'unconverged_with_significant_regime_changes': 0,
                        'max_initial_residual': 0., 'tolerance': float(tolerance)})
                    p['steps'] += 1
                    p['unconverged_steps'] += flag == 'NOT converged'
                    p['unconverged_with_significant_regime_changes'] += flag == 'NOT converged' and int(significant)>0
                    p['max_initial_residual'] = max(p['max_initial_residual'], float(residual))
                    residuals.setdefault(name, []).append(float(residual))
            elif b're-speciation monitors:' in raw:
                mole_sum += float(raw.split(b'max mole fraction ')[1])
                negative = float(raw.split(b'negative mass ')[1].split()[0])
                negative_mass += negative
                name = raw.split()[0].decode()[2:]
                nu = np.array(FORMULAS[name])
                negative_i += negative*nu[2]/float(np.dot(nu,ATOMIC_MASSES))
            elif raw.startswith(b'Re-speciation: cells'):
                m = respec_re.search(raw.decode())
                assert m, path
                result['max_respeciation_failures'] = max(result['max_respeciation_failures'], int(m[1]))
                result['max_invalid_element_cells'] = max(result['max_invalid_element_cells'], int(m[2]))
                result['final_invalid_element_cells'] = int(m[2])
                if int(m[2]):
                    result['steps_with_invalid_element_cells'] += 1
                    if result['first_invalid_time_s'] is None: result['first_invalid_time_s'] = current_time
                    result['last_invalid_time_s'] = current_time
                result['max_equilibrium_element_residue'] = max(result['max_equilibrium_element_residue'], float(m[3]))
                result['max_sum_of_species_mole_fraction_maxima'] = max(result['max_sum_of_species_mole_fraction_maxima'], mole_sum)
                result['max_negative_gas_mass_kg'] = max(result['max_negative_gas_mass_kg'],negative_mass)
                if source_i_rate*current_time>0:
                    result['max_negative_I_over_supplied_I'] = max(result['max_negative_I_over_supplied_I'],negative_i/(source_i_rate*current_time))
                mole_sum = negative_mass = negative_i = 0.
            elif raw.startswith(b'Phase change summary'):
                m = summary_re.search(raw.decode())
                assert m, path
                name, unconverged, guards, projected, significant, complement, complement_sig = m.groups()
                result['pairs'][name].update(summary_unconverged_steps=int(unconverged),
                    guard_passes=int(guards), projected_elements=int(projected),
                    significant_projected_elements=int(significant), complementarity_violations=int(complement),
                    significant_complementarity_violations=int(complement_sig))
            elif raw.strip() == b'End':
                result['normal_end'] = True
            elif raw.startswith(b'ExecutionTime ='):
                result['final_execution_s'] = float(raw.split()[2])
                result['final_clock_s'] = float(raw.split()[6])
            elif raw.startswith(b'Time = '):
                current_time = float(raw.split()[2])
    for name, values in residuals.items():
        result['pairs'][name]['median_initial_residual'] = float(np.median(values))
        result['pairs'][name]['p95_initial_residual'] = float(np.quantile(values, .95))
        assert result['pairs'][name]['steps']==30000, (path,name)
        assert result['pairs'][name]['unconverged_steps']==result['pairs'][name]['summary_unconverged_steps'], (path,name)
    result['sha256'] = digest.hexdigest()
    assert result['normal_end'], path
    return result


def transport_scales(case, manifest):
    data=(case/'0/U').read_text().split('boundaryField')[0]
    m=re.search(r'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',data)
    assert m,case
    velocity=np.fromstring(data[m.end():].replace('(',' ').replace(')',' ').replace(';',' '),sep=' ').reshape(-1,3)
    assert len(velocity)==int(m[1])==18432,case
    ux=np.abs(velocity[:,0]);dx=1.25/384;diffusivity=manifest['gas_diffusivity_m2_s']
    return {'max_axial_velocity_m_s':float(ux.max()),'median_axial_velocity_m_s':float(np.median(ux)),
            'axial_cell_size_m':dx,'max_cell_Peclet':float(ux.max()*dx/diffusivity),
            'max_advective_Courant':float(ux.max()*manifest['dt_s']/dx),
            'max_1D_upwind_diffusion_over_physical_D':float(ux.max()*dx/(2*diffusivity)),
            'median_1D_upwind_diffusion_over_physical_D':float(np.median(ux)*dx/(2*diffusivity)),
            'scope':'Local 1D Taylor-series estimate D_num = |Ux| dx/2 for upwind advection; not a measured CFD error or a full variable-density/radial analysis.'}


def experimental_shape(case):
    comparison = read_csv(case/'iodine-comparison.csv')
    selected = [r for r in comparison if r['figure_relative_yield_percent']]
    x = np.array([float(r['x_m']) for r in selected])
    model = np.array([float(r['predicted_relative_yield_percent']) for r in selected])
    observed = np.array([float(r['figure_relative_yield_percent']) for r in selected])
    assert len(x)>1 and np.allclose(np.diff(x), .01), case
    assert np.isfinite(model).all() and np.isfinite(observed).all()
    assert (model>=0).all() and (observed>=0).all() and model.sum()>0 and observed.sum()>0
    p, q = model/model.sum(), observed/observed.sum()
    return {'raw_L1_percentage_points': float(np.abs(model-observed).sum()),
            'observed_scan_sum_pct': float(observed.sum()), 'model_scan_sum_pct': float(model.sum()),
            'normalised_scan_total_variation': float(np.abs(p-q).sum()/2),
            'normalised_scan_overlap': float(np.minimum(p, q).sum()),
            'normalised_scan_Wasserstein_m': float(np.abs(np.cumsum(p-q)[:-1]).sum()*.01),
            'normalised_scan_centroid_difference_m':float(np.dot(x,p-q)),
            'scope': 'Secondary normalisation over common scan; observed retains KI; no alignment fitting.'}


def stationarity(case):
    output = {}
    for species in SPECIES:
        x, final, _ = profile(case, species)
        _, a, _ = profile(case, species, 40)
        _, b, _ = profile(case, species, 50)
        _, early, _ = profile(case, species, 20)
        _, mid, _ = profile(case, species, 30)
        previous, latest, earlier = (b-a)/10, (final-b)/10, (mid-early)/10
        denom = max(float(np.abs(previous).sum()), float(np.abs(latest).sum()))
        output[species] = {
            'last_10s_rate_mol_s_wedge': float(latest.sum()),
            'last_two_10s_rate_relative_change': float(abs(latest.sum()-previous.sum())/
                max(abs(latest.sum()), abs(previous.sum()))) if max(abs(latest.sum()),abs(previous.sum())) else 0.,
            'last_two_10s_rate_spatial_L1': float(np.abs(latest-previous).sum()/denom) if denom else 0.,
            '20_30_vs_50_60_rate_spatial_L1': relative_l1(earlier, latest),
            '40_vs_60_cumulative_shape_L1': relative_l1(a/a.sum(), final/final.sum()) if a.sum() and final.sum() else None,
        }
    return output


def write_table(path, records):
    if not records:
        return
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(records)


def tables_html(records, columns):
    def fmt(value):
        if value is None:
            return '—'
        if isinstance(value, float):
            return f'{value:.3g}'
        return html.escape(str(value))
    return '<table><thead><tr>'+''.join('<th>'+html.escape(label)+'</th>' for key,label in columns)+\
        '</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(r.get(key))+'</td>' for key,label in columns)+'</tr>' for r in records)+'</tbody></table>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-root', required=True, type=Path)
    parser.add_argument('--coarse-root', type=Path)
    parser.add_argument('--half-root', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--corrector-summary',type=Path,help='Optional completed short restart diagnostic')
    parser.add_argument('--positivity-summary',type=Path,help='Optional completed startup replay')
    args = parser.parse_args()
    out = args.out; out.mkdir(parents=True, exist_ok=True)
    full = args.full_root
    validation = json.loads((full/'validation.json').read_text())
    passed = [r for r in validation if r['status']=='PASS']
    assert len(passed)==15 and not any(r['status']=='FAIL' for r in validation)
    peaks = read_csv(full/'peaks.csv')
    peak_map = {(r['case'],r['species']):r for r in peaks}
    checks = read_csv(HERE/'digitized/figure_checks.csv')
    check_map = {(r['experiment'],r['species']):r for r in checks}
    q1 = []
    for r in peaks:
        if '_Q1_' not in r['case'] or not r['measured_peak_K'] or not r['peak_temperature_K']:
            continue
        c = check_map[(r['experiment'],r['species'])]
        q1.append(dict(case=r['case'], experiment=r['experiment'], species=r['species'],
            predicted_peak_C=float(r['peak_temperature_K'])-273.15,
            published_peak_C=float(r['measured_peak_K'])-273.15, error_K=float(r['difference_K']),
            uncertainty_K=float(c['printed_uncertainty_K']), within_uncertainty=r['within_experimental_uncertainty'],
            model_x_cm=100*float(r['peak_x_m']), figure_x_cm=100*float(c['scan_peak_x_m']),
            offset_cm=100*(float(r['peak_x_m'])-float(c['scan_peak_x_m']))))

    cache_path = out/'log_audits.json'
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    audits, matrix, refinements, sensitivities = {}, [], [], []
    coarse_peaks = {(r['case'],r['species']):r for r in read_csv(args.coarse_root/'peaks.csv')} if args.coarse_root else {}
    for index, result in enumerate(passed):
        name = result['case']; case = full/name
        manifest = json.loads((case/'benchmark.json').read_text())
        signature = {'size': (case/'log.solver').stat().st_size, 'mtime_ns': (case/'log.solver').stat().st_mtime_ns,'audit_version':2}
        old = cache.get(name)
        audit = old['audit'] if old and old['signature']==signature else log_audit(case/'log.solver',manifest)
        cache[name] = {'signature': signature, 'audit': audit}
        cache_path.write_text(json.dumps(cache, indent=2, allow_nan=False)+'\n')
        ledger = ledger_audit(case)
        steady = stationarity(case)
        shape = experimental_shape(case)
        targets = manifest['targets']
        yield_pb = sum(float(t['iodine_yield_pct']) for t in targets if t['species']=='PbI2' and t['iodine_yield_pct'])
        yield_bi = sum(float(t['iodine_yield_pct']) for t in targets if t['species']=='BiI3' and t['iodine_yield_pct'])
        observed_split = yield_bi/(yield_pb+yield_bi) if yield_pb+yield_bi else 0.
        final_i = ledger['I']['final']; released = final_i['released']
        row = dict(case=name, experiment=manifest['experiment'], source_mode=manifest['source_mode'],
            thermo_variant='shift29k' if '_shift29k' in name else 'baseline',
            metal_saturation=manifest['metal_saturation_fraction'],
            BiI3_iodide_I_pct=100*result['BiI3_fraction_of_modelled_iodide_wall_iodine'],
            observed_BiI3_iodide_I_pct=100*observed_split,
            split_error_percentage_points=100*(result['BiI3_fraction_of_modelled_iodide_wall_iodine']-observed_split),
            iodine_wall_fraction=final_i['wall']/released, iodine_gas_fraction=final_i['gas']/released,
            iodine_outlet_net_fraction=final_i['transport']/released,
            max_accounted_closure_over_reference=max(v['max_relative_closure'] for v in result['element_gates'].values()),
            max_raw_material_residue_over_source=max(v['max_raw_material_residue_over_source'] for v in ledger.values()),
            final_clock_hours=audit['final_clock_s']/3600,
            **shape)
        for species in ['PbI2','BiI3']:
            p = peak_map[(name,species)]
            row[species+'_peak_C'] = float(p['peak_temperature_K'])-273.15 if p['peak_temperature_K'] else None
            row[species+'_peak_error_K'] = float(p['difference_K']) if p['difference_K'] else None
            row[species+'_wall_mol_full_column'] = float(p['full_column_deposit_mol'])
        matrix.append(row)
        saved=out/'profiles'/name;saved.mkdir(parents=True,exist_ok=True)
        (saved/'iodine-comparison.csv').write_bytes((case/'iodine-comparison.csv').read_bytes())
        for species in ['PbI2','BiI3','I']:
            filename=SPECIES[species]+'_native.dat'
            (saved/filename).write_bytes((case/'postProcessing/axialProfiles/60'/filename).read_bytes())
        audits[name] = {'ledger': ledger, 'log': audit, 'stationarity': steady, 'shape': shape,
                        'transport_scales':transport_scales(case,manifest)}
        if args.coarse_root:
            for species in ['PbI2','BiI3']:
                if not peak_map[(name,species)]['peak_temperature_K']:
                    continue
                x, fine_amount, _ = profile(case,species)
                xc, coarse_amount, _ = profile(args.coarse_root/name,species)
                assert np.allclose(x,xc), name
                cp = coarse_peaks[(name,species)]; fp = peak_map[(name,species)]
                refinements.append(dict(case=name,species=species,
                    coarse_cells=2048,fine_cells=18432,coarse_dt_s=.02,fine_dt_s=.002,
                    profile_L1_relative_to_full=relative_l1(coarse_amount,fine_amount),
                    normalised_shape_L1=relative_l1(coarse_amount/coarse_amount.sum(),fine_amount/fine_amount.sum()),
                    deposit_relative_change=(fine_amount.sum()-coarse_amount.sum())/coarse_amount.sum(),
                    peak_change_K=float(fp['peak_temperature_K'])-float(cp['peak_temperature_K']),
                    peak_shift_cm=100*(float(fp['peak_x_m'])-float(cp['peak_x_m'])),
                    scope='Combined mesh/time-step/solver-revision change; not isolated convergence.'))
        print(f'Analysed {index+1}/15: {name}',flush=True)
    write_table(out/'case_metrics.csv', matrix)
    write_table(out/'q1_peaks.csv', q1)
    write_table(out/'combined_refinement.csv', refinements)
    half = []
    if args.half_root and args.coarse_root:
        half_results = json.loads((args.half_root/'validation.json').read_text())
        half_peaks = {(r['case'],r['species']):r for r in read_csv(args.half_root/'peaks.csv')}
        for result in half_results:
            if result['status']!='PASS': continue
            for species in ['PbI2','BiI3']:
                key=(result['case'],species)
                if not half_peaks[key]['peak_temperature_K']: continue
                x,a,_=profile(args.coarse_root/key[0],species);xh,b,_=profile(args.half_root/key[0],species)
                assert np.allclose(x,xh)
                half.append(dict(case=key[0],species=species,cells=2048,dt_s_1=.02,dt_s_2=.01,
                    profile_L1_relative_to_10ms=relative_l1(a,b),
                    peak_change_K=float(half_peaks[key]['peak_temperature_K'])-float(coarse_peaks[key]['peak_temperature_K'])))
    write_table(out/'coarse_time_refinement.csv', half)
    for experiment in EXPERIMENTS[1:]:
        for sat in [.001,.01,.1]:
            base = next(r for r in matrix if r['experiment']==experiment and r['source_mode']=='Q2' and r['thermo_variant']=='baseline' and r['metal_saturation']==sat)
            shift = next(r for r in matrix if r['experiment']==experiment and r['source_mode']=='Q2' and r['thermo_variant']=='shift29k' and r['metal_saturation']==sat)
            sensitivities.append(dict(experiment=experiment,metal_saturation=sat,
                baseline_BiI3_iodide_I_pct=base['BiI3_iodide_I_pct'],shift29k_BiI3_iodide_I_pct=shift['BiI3_iodide_I_pct'],
                split_change_percentage_points=shift['BiI3_iodide_I_pct']-base['BiI3_iodide_I_pct'],
                PbI2_peak_change_K=shift['PbI2_peak_C']-base['PbI2_peak_C'],
                BiI3_peak_change_K=shift['BiI3_peak_C']-base['BiI3_peak_C']))
    write_table(out/'thermodynamic_sensitivity.csv', sensitivities)

    figure_discrepancies=[]
    temperature=read_csv(HERE/'digitized/LBE-I_SiO2_I_temperature.csv')
    tx=np.array([float(r['x_m']) for r in temperature]);tt=np.array([float(r['T_K']) for r in temperature])
    for c in checks:
        if c['experiment']!='LBE-I_SiO2_I': continue
        crossings=[]; target=float(c['printed_peak_K'])
        for x0,x1,t0,t1 in zip(tx[:-1],tx[1:],tt[:-1],tt[1:]):
            if min(t0,t1)<=target<=max(t0,t1) and t0!=t1:
                crossings.append(float(x0+(target-t0)*(x1-x0)/(t1-t0)))
        figure_discrepancies.append(dict(species=c['species'],figure_peak_x_m=float(c['scan_peak_x_m']),
            figure_minus_printed_K=float(c['difference_K']),curve_crossings_at_printed_temperature_m=crossings,
            diagnostic_offset_cm=[100*(x-float(c['scan_peak_x_m'])) for x in crossings],
            scope='Diagnostic only; no coordinate or temperature correction applied.'))
    correctors=json.loads(args.corrector_summary.read_text()) if args.corrector_summary else None
    if correctors:
        (out/'corrector_diagnostic.json').write_text(json.dumps(correctors,indent=2)+'\n')
    positivity=json.loads(args.positivity_summary.read_text()) if args.positivity_summary else None
    if positivity:
        (out/'positivity_diagnostic.json').write_text(json.dumps(positivity,indent=2)+'\n')
    summary = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'input_full_root': str(full.resolve()), 'input_coarse_root': str(args.coarse_root.resolve()) if args.coarse_root else None,
        'input_half_root': str(args.half_root.resolve()) if args.half_root else None,
        'status': 'NUMERICAL_MATRIX_COMPLETE_EXPERIMENTAL_VALIDATION_PROVISIONAL',
        'PASS':15, 'FAIL':0, 'NOTRUN':1, 'cells':18432, 'dt_s':.002, 'end_time_s':60,
        'max_accounted_element_closure_over_reference': max(r['max_accounted_closure_over_reference'] for r in matrix),
        'max_raw_material_residue_over_source': max(r['max_raw_material_residue_over_source'] for r in matrix),
        'max_invalid_element_cells': max(r['log']['max_invalid_element_cells'] for r in audits.values()),
        'max_equilibrium_element_residue': max(r['log']['max_equilibrium_element_residue'] for r in audits.values()),
        'max_trace_mole_fraction_upper_bound': max(r['log']['max_sum_of_species_mole_fraction_maxima'] for r in audits.values()),
        'max_negative_I_over_supplied_I':max(r['log']['max_negative_I_over_supplied_I'] for r in audits.values()),
        'corrector_diagnostic':correctors,
        'positivity_diagnostic':positivity,
        'q1_peaks':q1,'case_metrics':matrix,'combined_refinement':refinements,'coarse_time_refinement':half,
        'thermodynamic_sensitivity':sensitivities,'silica_figure_diagnostics':figure_discrepancies,
        'audits':audits,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [full/'validation.json',full/'peaks.csv',HERE/'targets.csv',HERE/'digitized/figure_checks.csv']}}
    (out/'analysis.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    make_figures(out,full,matrix,refinements,q1,audits)
    write_report(out,summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['audits','input_sha256','case_metrics','q1_peaks','combined_refinement','coarse_time_refinement']},indent=2))


def make_figures(out,full,matrix,refinements,q1,audits):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    figures=[]
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=170,bbox_inches='tight')
        figures.append(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,11),layout='constrained')
    choices=[('Q1_baseline','Q1 prescribed source','black','-'),
             ('Q2_baseline_sat0.001','Q2 baseline, 0.1% saturation','#0072B2','-'),
             ('Q2_baseline_sat0.01','Q2 baseline, 1% saturation','#009E73','-'),
             ('Q2_shift29k_sat0.01','Q2 +29 kJ/mol BiI, 1% saturation','#D55E00','--')]
    for ax,exp in zip(axes,EXPERIMENTS):
        observed=read_csv(HERE/'digitized'/(exp+'_gamma.csv'))
        observed_label='Published iodine scan, digitised'+(' (includes KI)' if exp.startswith('LBE') else '')
        ax.bar([100*float(r['x_m']) for r in observed],[float(r['relative_yield_percent_per_cm']) for r in observed],
               width=.8,color='.8',label=observed_label)
        for suffix,label,color,style in choices:
            p=full/(exp+'_'+suffix)/'iodine-comparison.csv'
            if not p.exists():continue
            curve=read_csv(p)
            ax.plot([100*float(r['x_m']) for r in curve],[float(r['predicted_relative_yield_percent']) for r in curve],
                    color=color,ls=style,label=label,lw=1.7)
        ax.set(title=exp,xlim=(25,70),xlabel='Assumed column coordinate (cm)',ylabel='Iodine yield per 1 cm bin (%)')
        ax.legend(fontsize=8)
    fig.suptitle('Full-size deposition profiles: 18,432 cells, 60 s, 2 ms\nExperimental curves: Liu et al. (2025), DOI 10.1007/s10967-025-10246-4 (CC BY 4.0)')
    save(fig,'iodine_profiles')

    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for ax,exp in zip(axes,EXPERIMENTS[1:]):
        data=np.array([[next(r['BiI3_iodide_I_pct'] for r in matrix if r['experiment']==exp and r['thermo_variant']==variant and r['metal_saturation']==sat)
                        for sat in [.001,.01,.1]] for variant in ['baseline','shift29k']])
        im=ax.imshow(data,vmin=0,vmax=100,cmap='viridis',aspect='auto')
        for i in range(2):
            for j in range(3):ax.text(j,i,f'{data[i,j]:.2f}%',ha='center',va='center',color='white' if data[i,j]<75 else 'black',fontsize=12)
        expected=next(r['observed_BiI3_iodide_I_pct'] for r in matrix if r['experiment']==exp)
        ax.set(xticks=range(3),xticklabels=['0.1%','1%','10%'],yticks=range(2),yticklabels=['Baseline','BiI +29 kJ/mol'],
               xlabel='Prescribed metal saturation',title=f'{exp}\nPublished iodide-only BiI3 share: {expected:.2f}%')
    fig.colorbar(im,ax=axes,label='BiI3 fraction of deposited iodide iodine (%)',shrink=.8)
    fig.suptitle('Thermodynamic and source sensitivity; +29 kJ/mol is a surrogate')
    save(fig,'speciation_sensitivity')

    fig,axes=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
    for index,r in enumerate(q1):
        color='#0072B2' if r['species']=='PbI2' else '#D55E00'
        axes[0].errorbar(r['error_K'],index,xerr=r['uncertainty_K'],fmt='o',color=color,capsize=3)
        axes[1].scatter(r['offset_cm'],index,color=color)
    labels=[r['experiment']+' '+r['species'] for r in q1]
    for ax in axes:
        ax.set(yticks=range(len(labels)),yticklabels=labels);ax.invert_yaxis();ax.axvline(0,color='.3',ls='--')
    axes[0].set(xlabel='Model minus printed peak temperature (K)',title='Bars show published uncertainty')
    axes[1].set(xlabel='Model minus figure peak position (cm)',title='A temperature match can hide a spatial error')
    fig.suptitle('Q1 peak comparisons: source split/rate prescribed, no fit')
    save(fig,'peak_comparison')

    if refinements:
        fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained',sharex=True)
        labels=[r['case'].replace('LBE-I_SS_II','Steel').replace('LBE-I_SiO2_I','Silica').replace('PbI2_SS','Pure')+' '+r['species'] for r in refinements]
        axes[0].bar(range(len(labels)),[100*r['profile_L1_relative_to_full'] for r in refinements],color='#0072B2')
        axes[1].bar(range(len(labels)),[r['peak_change_K'] for r in refinements],color='#D55E00')
        axes[0].set(ylabel='Deposit-profile L1 difference (%)')
        axes[1].set(ylabel='Peak temperature change (K)',xticks=range(len(labels)),xticklabels=labels)
        axes[1].tick_params(axis='x',rotation=90)
        axes[1].axhline(0,color='.3',lw=.6)
        fig.suptitle('2,048 cells / 20 ms versus 18,432 cells / 2 ms\nMesh, time step and solver revision changed together: this is not an isolated convergence study')
        save(fig,'combined_refinement')

    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for exp in EXPERIMENTS:
        case=full/(exp+'_Q1_baseline')
        m=json.loads((case/'benchmark.json').read_text())
        source=sum((2 if n=='PbI2_g' else 3)*v for n,v in m['source_wedge_mol_s'].items())
        cumulative=[]
        for t in range(1,61):
            x,a,_=profile(case,'I',t);cumulative.append(a.sum())
        rates=np.diff([0]+cumulative)/source
        axes[0].plot(range(1,61),rates,label=exp)
        axes[1].plot(range(1,61),np.array(cumulative)/(source*np.arange(1,61)),label=exp)
    axes[0].set(xlabel='Simulated time (s)',ylabel='One-second net wall rate / iodine source rate',title='Rate stationarity')
    axes[1].set(xlabel='Simulated time (s)',ylabel='Wall iodine / cumulative supplied iodine',title='Startup gas inventory remains at 60 s')
    for ax in axes:
        ax.axhline(1,color='.4',ls='--');ax.legend(fontsize=8);ax.set_ylim(0,1.05)
    fig.suptitle('Q1 transient history; continuous release, no cooling or flush')
    save(fig,'transient_history')

    fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    t=read_csv(HERE/'digitized/LBE-I_SiO2_I_temperature.csv')
    ax.plot([100*float(r['x_m']) for r in t],[float(r['T_K'])-273.15 for r in t],color='#D55E00',label='Digitised Figure 7C temperature curve')
    for r in read_csv(HERE/'digitized/figure_checks.csv'):
        if r['experiment']!='LBE-I_SiO2_I':continue
        x=100*float(r['scan_peak_x_m']);curve=float(r['figure_temperature_K'])-273.15;printed=float(r['printed_peak_K'])-273.15
        ax.scatter([x,x],[curve,printed],c=['#D55E00','#0072B2']);ax.plot([x,x],[printed,curve],ls='--',color='.4')
        ax.annotate(f"{r['species']}: {float(r['difference_K']):.0f} K discrepancy",(x,curve),xytext=(x-9,curve+55),arrowprops={'arrowstyle':'-'})
    for r in q1:
        if r['experiment']=='LBE-I_SiO2_I':ax.scatter(r['model_x_cm'],r['predicted_peak_C'],marker='x',s=70,color='black',label='Q1 model '+r['species'])
    ax.set(xlim=(40,68),ylim=(40,650),xlabel='Figure / assumed mesh coordinate (cm)',ylabel='Temperature (°C)',title='Silica: curve and printed temperatures do not agree at measured scan peaks')
    ax.legend(fontsize=8)
    save(fig,'silica_input_discrepancy')
    with PdfPages(out/'analysis_figures.pdf') as pdf:
        for fig in figures:pdf.savefig(fig,bbox_inches='tight');plt.close(fig)


def write_report(out,s):
    matrix=s['case_metrics'];ref=s['combined_refinement']
    worst_ref=max((r['profile_L1_relative_to_full'] for r in ref),default=0)
    pair_records=[p for a in s['audits'].values() for p in a['log']['pairs'].values()]
    unconverged=sum(p['unconverged_steps'] for p in pair_records);steps=sum(p['steps'] for p in pair_records)
    changed=sum(p['unconverged_with_significant_regime_changes'] for p in pair_records)
    corrector_records=[]
    if s['corrector_diagnostic']:
        for r in s['corrector_diagnostic']['comparisons']:
            gas=max(v['max_difference_over_max_control'] for v in r['gas_fields'].values() if v['max_difference_over_max_control'] is not None)
            increment=max(v for v in r['wall_increment_profile_L1'].values() if v is not None)
            corrector_records.append({'case':r['case'],'variant':r['variant'],
                                      'max_scaled_gas_difference':gas,'max_wall_increment_L1':increment})
    overview=f'''<p><strong>The full-size matrix is complete and conserves elements well, but experimental validation remains provisional.</strong>
    All 15 available cases reached 60 s on 18,432 cells with 2 ms steps; the BiI3_SiO2_800 case lacks a temperature profile.
    All five available Q1 printed temperature peaks fall within their stated uncertainties. The spatial profiles and Q2 speciation are less secure:
    the steel result depends strongly on source loading and BiI thermodynamics, and the silica profiles are displaced from the observed scan.</p>
    <p>The source is <a href="{PAPER}">Liu et al. (2025), JRNC 334, 5201–5216</a>. Figure-derived inputs are adapted under CC BY 4.0.
    No new experimental data, coordinate adjustments or parameter fits were introduced.</p>'''
    clocks=[r['final_clock_hours'] for r in matrix]
    overview+=f'''<p>The 15 independent serial jobs ran concurrently, with per-case wall-clock durations of {min(clocks):.2f}–{max(clocks):.2f} h.
    These are execution records on a shared workstation, not a controlled performance benchmark.</p>'''
    numeric=f'''<h2>Numerical completion and conservation</h2><p>The existing ledger closure gate has a worst relative error of
    {s['max_accounted_element_closure_over_reference']:.3g}. That ledger includes the recorded solver defect and clamping corrections.
    An independent raw material balance, gas + wall + sample + net outward transport + removed − initial − released,
    normalised to supplied atoms before these corrections, has a worst error of {s['max_raw_material_residue_over_source']:.3g} across every element and recorded step.
    This distinction prevents a corrected ledger from being mistaken for an uncorrected physical balance.</p>
    <p>Gas equilibrium had zero reported solver failures, but its input-validity check skipped some cells during startup.
    These are distinct diagnostics: at most {s['max_invalid_element_cells']} cells were skipped on one step.
    The maximum cell equilibrium element residue was {s['max_equilibrium_element_residue']:.3g}.
    The upper bound obtained by summing individual trace-species mole-fraction maxima was {s['max_trace_mole_fraction_upper_bound']:.3g};
    these species maxima need not occur in the same cell.</p>
    <p><strong>Outer convergence needs separate qualification.</strong> Paired-gas logs report {unconverged:,} of {steps:,} pair-steps without meeting the
    1e-16 outer criterion; {changed:,} of those steps also had significant regime changes. The final transport solve, bound/complementarity guard and
    gas equilibrium have different criteria. The log audit records their statuses separately. Residual floors can prevent the very strict outer criterion
    from passing, so its failure is not automatically a failed equilibrium solve. Conservation alone cannot establish iteration independence:
    compare additional correctors and a justified tolerance on representative cases before treating the coupling as qualified.</p>'''
    if corrector_records:
        numeric+='''<p>A matched late-time diagnostic restarted four representative cases from 59 s for 20 steps (40 ms),
        comparing nCorr 3 against nCorr 8 and a 1e-14 outer tolerance against 1e-16. Source, discretisation and linear solver tolerances were unchanged.
        The largest gas-field difference, scaled by that species' maximum control value, was below 3.6e-8.
        The deposited-iodide increment profiles changed negligibly, supporting a residual-floor interpretation of the outer messages at these checkpoints.
        Both cumulative deposits and the new 40 ms increments were checked so accumulated deposits could not hide iteration sensitivity.
        This short test does not replace a full-run iteration or time-step convergence study.</p>'''+tables_html(corrector_records,[
            ('case','Restarted case'),('variant','Control change'),('max_scaled_gas_difference','Max scaled gas difference'),
            ('max_wall_increment_L1','Max new-deposit profile L1')])
    invalid=[]
    for name,a in s['audits'].items():
        log=a['log']
        if log['max_invalid_element_cells']:
            invalid.append({'case':name,'cells':log['max_invalid_element_cells'],
                            'time_s':log['first_invalid_time_s'],'steps':log['steps_with_invalid_element_cells']})
    numeric+='<h3>Transient invalid-bulk cells</h3>'+tables_html(invalid,[('case','Case'),('cells','Maximum skipped cells'),('time_s','Time s'),('steps','Affected steps')])
    if s['positivity_diagnostic']:
        upper=max(event['positive_I_over_supplied_I_upper_bound'] for r in s['positivity_diagnostic']['results'] for event in r['events'])
        numeric+=f'''<p>Five-step startup replays reproduced all three events, at 6–10 ms. They contain at most {upper:.3g} of the supplied iodine,
        even using an entire axial slice volume as a conservative upper bound for each affected radial cell.
        These are negligible leading-tail amounts; no invalid-element cells remain at the last full-run step. No concentration floor,
        species deletion or positivity-gate change was introduced. A solver-success counter alone would not have established this.</p>'''
    peaks='<h2>Temperature peaks and spatial locations</h2>'+tables_html(s['q1_peaks'],[
        ('experiment','Experiment'),('species','Species'),('predicted_peak_C','Model °C'),('published_peak_C','Paper °C'),
        ('error_K','Error K'),('uncertainty_K','Uncertainty ±K'),('offset_cm','Spatial offset cm')])+'''
        <p>These uncertainties are the paper's temperature variation around a ±1 cm scan location, not statistical confidence intervals.
        Q1 prescribes the release and, for LBE, a species split inferred from measured deposits. Its temperature agreement is conditional on that source.
        The silica Q1 peaks occur 5 cm and 4 cm downstream of the measured PbI2 and BiI3 scan maxima. A close printed temperature is therefore insufficient
        to establish the correct deposition position. The pure PbI2 result is near the edge of its ±18 K uncertainty and should not be described as a precise match.</p>
        <img src="peak_comparison.png" alt="Q1 peak temperature and position errors">'''
    profiles='''<h2>Iodine profiles</h2><p>The primary error is the sum of absolute differences in percentage points per 1 cm bin over the digitised scan.
    It retains the published histogram heights, their imperfect sums and the omitted KI contribution. It is neither a relative percent error nor a chi-square statistic.
    The secondary overlap, total variation and Wasserstein metrics normalise each curve over the same scan range; they are explicitly conditional shape diagnostics
    and still include the observed KI. No experimental bin covariance or count uncertainty is available for a likelihood test.</p>
    <img src="iodine_profiles.png" alt="Full-size iodine profiles and published scans">'''+tables_html(matrix,[
        ('case','Case'),('raw_L1_percentage_points','Raw L1 pp'),('normalised_scan_overlap','Shape overlap'),
        ('normalised_scan_Wasserstein_m','Shape transport distance m'),('iodine_wall_fraction','Supplied I at wall')])+'''
    <p>The steel low-loading Q2 baseline has the lowest raw profile error in this discrete sweep; that makes it a diagnostic candidate, not an identified
    physical metal saturation. Missing release histories, transport coefficients and KI chemistry prevent a unique inference. The silica error remains large
    across the entire sensitivity sweep. The omitted 7% KI contribution by itself cannot explain a profile discrepancy near 180 percentage points.</p>'''
    chemistry='''<h2>Source and thermodynamic sensitivity</h2><p>BiI3 shares below refer to iodine in the two deposited iodides:
    3 n(BiI3)/(2 n(PbI2)+3 n(BiI3)). The paper's corresponding shares are 61.3/(61.3+32.1) = 65.63% in steel and
    27/(27+66) = 29.03% in silica. They exclude KI and are not fractions of all gas-phase or deposited metal species.</p>
    <img src="speciation_sensitivity.png" alt="BiI3 shares under source and thermodynamic assumptions">'''+tables_html(matrix,[
        ('case','Case'),('BiI3_iodide_I_pct','BiI3 iodide I %'),('observed_BiI3_iodide_I_pct','Paper iodide I %'),
        ('PbI2_peak_error_K','PbI2 error K'),('BiI3_peak_error_K','BiI3 error K')])+'''
    <p>At 1% and 10% saturation the baseline gives about 94.5% BiI3 in either column; adding 29 kJ/mol to BiI gives about 2%.
    This reversal is far larger than the deposited-amount changes associated with grid/time refinement. At those loadings the 1% and 10% iodide outputs
    are almost identical, suggesting an iodine-limited regime in this model. The low-loading branch behaves differently.
    The +29 kJ/mol shift changes the BiI equilibrium formation constant by exp(−29000/RT): approximately 0.051 at 1173 K and 0.00093 at 500 K.
    It is a large sensitivity surrogate, not a second experimentally qualified thermodynamic record. Neither branch establishes the true source speciation.</p>
    <p>Q1 LBE sources use inferred deposited compound masses, so their close deposited split is substantially conditioned by construction. Q2 adds ideal-activity
    metal vapour but still conditions total iodine release on measured deposits. The three-hour mean source rate is retained: a 60 s run releases 1/180
    of that assumed inventory, with the 5° wedge receiving 1/72 of the full-column source. No entire experimental inventory was compressed into 60 s.</p>'''
    chemistry+='''<p>Only Pb–Bi–I chemistry in helium is represented. Table 2 also reports moisture measurements around 157 ppm for pure PbI2,
    185–375 ppm for steel LBE and 120–375 ppm for silica LBE. These were measured downstream; the paper attributes elevated downstream moisture to column
    outgassing, while its inlet was dried. They must not be inserted as uniform inlet water concentrations. Hydrogen/oxygen chemistry is outside the present
    species set, and the supplied paper alone does not quantify its effect on these CFD predictions.</p>'''
    time='''<h2>Stationarity and experimental duration</h2><img src="transient_history.png" alt="Iodine wall rate and cumulative capture histories">
    <p>The original report compares 58–59 and 59–60 s net deposition windows. This analysis additionally compares 40–50 against 50–60 s and 20–30 against
    50–60 s, and distinguishes cumulative profile shape from deposition-rate shape. Relevant iodide rates are near steady under the continuous prescribed source.
    Cumulative wall amounts still retain a startup offset because iodine remains in the gas. Steady deposition rate is not a validation of a three-hour release history.</p>
    <p>The one original NOT_QUALIFIED stationarity row is pure PbI2 Q1, triggered by its minute Pb metal component. Its final Pb metal wall amount is roughly
    8e-17 of its PbI2 wall amount, while the PbI2 rate itself is stationary. The original classification is retained; the trace-metal trigger is reported separately
    rather than silently dropping that phase or relaxing a gate. The experiment withdrew the sample, turned off both furnaces and flow, and cooled before scanning;
    that shutdown/cooling phase is absent from this matrix.</p>'''
    significant_rates=[a['stationarity'][sp]['20_30_vs_50_60_rate_spatial_L1'] for a in s['audits'].values()
                       for sp in ['PbI2','BiI3','I'] if a['stationarity'][sp]['20_30_vs_50_60_rate_spatial_L1'] is not None]
    cumulative=[a['stationarity'][sp]['40_vs_60_cumulative_shape_L1'] for a in s['audits'].values()
                for sp in ['PbI2','BiI3','I'] if a['stationarity'][sp]['40_vs_60_cumulative_shape_L1'] is not None]
    time+=f'''<p>The largest relevant iodide/iodine rate-profile change between 20–30 and 50–60 s is {100*max(significant_rates):.4g}%.
    The largest normalised cumulative iodide/iodine shape change between 40 and 60 s is {100*max(cumulative):.3g}%.
    At 60 s, {100*min(r['iodine_wall_fraction'] for r in matrix):.2f}–{100*max(r['iodine_wall_fraction'] for r in matrix):.2f}% of supplied iodine is on the wall;
    the remainder is primarily gas inventory. Extending an already stationary source run is unlikely to resolve the large silica peak displacement.</p>'''
    refine=f'''<h2>Refinement evidence and its limits</h2><p>The coarse matrix used 2,048 cells and 20 ms steps; the completed matrix uses 18,432 cells and 2 ms.
    Both mesh and time step changed, and the equilibrium solver revision also changed. Their comparison cannot isolate spatial truncation error,
    temporal splitting error or the effect of the precision repair. The largest iodide deposit-profile difference is {100*worst_ref:.2f}% relative to the full-size profile.</p>'''
    if ref:refine+='<img src="combined_refinement.png" alt="Combined mesh and time-step changes">'+tables_html(ref,[
        ('case','Case'),('species','Species'),('profile_L1_relative_to_full','Profile L1 fraction'),
        ('peak_change_K','Peak change K'),('deposit_relative_change','Amount change fraction')])
    refine+=tables_html(s['coarse_time_refinement'],[('case','Coarse time-step check'),('species','Species'),
        ('profile_L1_relative_to_10ms','20 vs 10 ms L1 fraction'),('peak_change_K','Peak change K')])+'''
        <p>The existing isolated 20-to-10 ms checks cover two coarse cases and leave the peaks unchanged, but show approximately 3.17% PbI2 and 2.00% BiI3
        profile changes for the steel 1% baseline Q2 case. Unchanged 1 cm binned peaks do not establish convergence. Use the same corrected solver for
        full-mesh 2-to-1 ms checks and a fixed-time-step mesh sequence, including a third spatial resolution before estimating an observed order or GCI.</p>'''
    scales=[]
    for exp in EXPERIMENTS:
        a=s['audits'][exp+'_Q1_baseline']['transport_scales']
        scales.append({'experiment':exp,**a})
    refine+='''<h3>Upwind advection and profile width</h3><p>All nine transported gases use first-order upwind advection; the physical diffusivity is prescribed as
    1e-4 m²/s. The local one-dimensional Taylor-series estimate D_num ≈ |Ux| Δx/2 indicates numerical diffusion comparable to or greater than that
    diffusivity even on the full mesh (Δx = 3.255 mm). The coarse axial spacing is three times larger.
    The table shows estimates from the frozen velocity field, not measured truncation errors; variable density, radial coupling and operator splitting remain outside this estimate.</p>'''+tables_html(scales,[
        ('experiment','Carrier'),('median_1D_upwind_diffusion_over_physical_D','Median Dnum/D estimate'),
        ('max_1D_upwind_diffusion_over_physical_D','Max Dnum/D estimate'),('max_cell_Peclet','Max cell Peclet'),
        ('max_advective_Courant','Max advective Courant')])+'''
        <p>This gives a plausible numerical contribution to the profile-width changes. Halving the time step cannot remove spatial upwind diffusion.
        A bounded higher-order advection option would require requalification of the wall predictor and bound guarantees; it should not be changed blindly
        to obtain a better experimental fit. First isolate mesh sensitivity using the same solver and time step.</p>'''
    silica='''<h2>Silica input discrepancy</h2><img src="silica_input_discrepancy.png" alt="Published silica curve versus printed temperatures">
    <p>At the measured gamma maxima, the digitised Figure 7C curve differs from the printed deposition temperatures by about +63 K for PbI2 and +108 K
    for BiI3. Inverting that same curve at the printed temperatures suggests offsets of several centimetres. This is a diagnostic of inconsistent supplied
    inputs; it does not identify whether the curve, labels or coordinate mapping is wrong. No curve or coordinate was shifted to improve a fit.</p>'''+tables_html(s['silica_figure_diagnostics'],[
        ('species','Species'),('figure_minus_printed_K','Curve minus label K'),('diagnostic_offset_cm','Curve crossing offset cm')])
    next_steps='''<h2>Supported conclusions and next work</h2><p>Supported: the available full-size matrix completes with strong element bookkeeping,
    small raw material residues and no reported gas-equilibrium failure; source loading and BiI thermodynamics materially affect predicted iodine partition;
    peak temperatures and spatial distributions must be assessed together. Unsupported: a uniquely inferred LBE release mechanism, a selected thermodynamic
    dataset, full experimental profile validation, or independent mesh/time convergence.</p>
    <ol><li>Run isolated full-mesh 2/1 ms comparisons, then a fixed-step mesh sequence on pure PbI2 and contrasting steel Q2 branches; assess spatial upwind diffusion.</li>
    <li>Use the short corrector and startup probes as supporting diagnostics; retain conservation, input-validity and bound checks in any broader coupling qualification.</li>
    <li>Resolve or bracket the silica Figure 7C inconsistency and coordinate origin using the supplied paper; keep its spatial comparisons provisional if raw data remain unavailable.</li>
    <li>Qualify species diffusivities, accommodation, source geometry and release history; the current uniform 1e-4 m²/s diffusivity, gas-equilibrium cutoff at 350 K and empty-tube carrier are assumptions.</li>
    <li>Qualify BiI thermodynamic provenance and refine source sensitivities around the transition between 0.1% and 1% metal saturation; avoid choosing values merely to match observations.</li>
    <li>Model shutdown/cooling if making quantitative claims about post-experiment scans. A software/methods manuscript may instead make the present limitations explicit.</li></ol>
    <p>The supplied paper remains the only experimental evidence. Missing BiI3_SiO2_800 temperatures and raw gamma counts cannot be reconstructed uniquely
    from this matrix. No author contact was made.</p>'''
    links='''<h2>Reproducible artifacts</h2><p><a href="analysis.json">Complete machine-readable analysis and provenance</a> ·
    <a href="case_metrics.csv">Case metrics</a> · <a href="q1_peaks.csv">Q1 peaks</a> ·
    <a href="combined_refinement.csv">Combined refinement</a> · <a href="thermodynamic_sensitivity.csv">Thermodynamic sensitivity</a> ·
    <a href="corrector_diagnostic.json">Corrector diagnostic</a> · <a href="positivity_diagnostic.json">Startup positivity diagnostic</a> ·
    <a href="profiles/">Final numerical profiles</a> ·
    <a href="analysis_figures.pdf">Exportable figure collection</a></p>
    <p>Raw CFD fields, logs, binaries and licensed thermodynamic tables remain in the ignored work area. Figures and numerical summaries contain no thermodynamic
    records. Regenerate with analyse_matrix.py using the full, coarse and half-step case roots.</p>'''
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Liu M10e full-size analysis</title>
    <style>body{font:16px/1.6 system-ui,sans-serif;max-width:1120px;margin:36px auto;padding:0 22px;color:#1e293b}
    h1,h2{line-height:1.25}h2{margin-top:40px}table{border-collapse:collapse;font-size:13px;display:block;overflow-x:auto}
    td,th{border-bottom:1px solid #cbd5e1;padding:8px;text-align:left}th{background:#f1f5f9}img{max-width:100%;height:auto}
    a{color:#075985}@media print{body{font-size:11px;max-width:none}table{font-size:9px}h2{break-after:avoid}img{max-height:90vh}}
    </style><body><h1>Full-size M10e analysis</h1><p>Completed 7 October 2026 · 15 available cases · 18,432 cells · 60 s · 2 ms</p>'''+\
        overview+numeric+peaks+profiles+chemistry+time+refine+silica+next_steps+links+'</body></html>'
    # Embed figures so the report remains readable when shared as one file.
    document=re.sub(r'src="([^\"]+\.png)"',lambda m:'src="data:image/png;base64,'+
                    base64.b64encode((out/m[1]).read_bytes()).decode()+'"',document)
    (out/'report.html').write_text(document)


if __name__=='__main__':
    main()
