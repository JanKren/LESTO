#!/usr/bin/env python3
"""Audit physical evidence and make explicitly conditional duration projections.

Restricted thermodynamic coefficients stay external. This audit reports
provenance, coverage and missing evidence, not a fitted or selected dataset.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from analyse_matrix import profile, moments, write_table, tables_html
from carriers import REPO
from prepare import load_references

BROOKS='https://sales.brooksinstrument.com/hubfs/Installation-Manual-5850E-EN.pdf?hsLang=en'
PAPER='https://doi.org/10.1007/s10967-025-10246-4'


def pbI2_D(T,p=101325):
    reduced=T/71.11785702515257
    omega=.5313-1.2946*math.exp(-.9922*reduced)+1.5591/reduced-.1918/reduced**2
    return 1.86e-7*T**1.5*math.sqrt(1/4.0026+1/461)/(p/101325*4.2584682449728355**2*omega)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--tf',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    dataset=json.loads(args.tf.read_text());coverage=[]
    for record in dataset['substances']:
        limits=[method['limitsTP'] for method in record.get('TPMethods',[]) if 'limitsTP' in method]
        coverage.append({'species':record['name'],'lower_K':min(limit['lowerT'] for limit in limits),
                         'upper_K':max(limit['upperT'] for limit in limits),
                         'source_status':'Constructed phase record' if record['name']=='BiI3(l)' else 'HERACLES/Barin compilation; independent confirmation pending'})
    diffusion=[{'T_K':T,'p_Pa':101325,'PbI2He_D_m2_s':pbI2_D(T),'ratio_to_constant_1e4':pbI2_D(T)/1e-4}
               for T in [293.15,400.15,427.15,561.15,579.15,641.15,1173.15]]
    experiments,_,_=load_references();projections=[];profiles=args.out/'duration_profiles';profiles.mkdir(exist_ok=True)
    matrix=json.loads((args.root/'matrix.json').read_text())
    for entry in matrix:
        if not (args.root/entry['case']/'numerical-report.json').exists():continue
        case=args.root/entry['case'];metadata=json.loads((case/'benchmark.json').read_text())
        assert metadata['end_time_s']==60
        duration=float(experiments[metadata['experiment']]['duration_s'])
        for species in ['PbI2','BiI3','I']:
            x,a,data=profile(case,species,60);_,previous,_=profile(case,species,50)
            rate=(a-previous)/10
            projected=a+(duration-60)*rate
            total=float(projected.sum());initial_total=float(a.sum())
            if not initial_total:continue
            # Report a failed extrapolation rather than applying a positivity floor.
            nonnegative=bool((projected>=0).all())
            target=profiles/(entry['case']+'_'+species+'.csv')
            write_table(target,[{'x_m':float(xx),'wall_60s_wedge_mol':float(aa),
                'late_rate_wedge_mol_s':float(rr),'projected_wall_wedge_mol':float(pp)}
                for xx,aa,rr,pp in zip(x,a,rate,projected)])
            mean=moments(x,projected) if nonnegative else None
            source_i=sum((2 if n=='PbI2_g' else 3 if n=='BiI3_g' else 1 if n=='I_g' else 2 if n=='I2_g' else 0)*v
                         for n,v in metadata['source_wedge_mol_s'].items())
            # Q2 sources use Pb, Bi and I atoms; no unrepresented iodine species.
            ledger=args.root/entry['case']/'postProcessing/phaseChangeElements/0/I.dat'
            from report import rows
            end=rows(ledger)[-1]
            projections.append({'case':entry['case'],'species':species,'projection_duration_s':duration,
                'status':'CONDITIONAL_NONNEGATIVE_PROJECTION' if nonnegative else 'INVALID_NEGATIVE_PROJECTED_BIN',
                'projected_full_column_wall_mol':total/metadata['carrier']['wedge_fraction'],
                'peak_x_cm':100*float(x[np.argmax(projected)]) if nonnegative else None,
                'centroid_cm':100*mean['centroid_m'] if mean and total else None,
                'residual_gas_I_over_3h_supply_if_inventory_steady':end['gas']/(source_i*duration),
                'scope':'Analytical continuation of 50–60 s net wall rates; not a 3 h CFD run or cooling simulation.'})
    write_table(args.out/'diffusion_screening.csv',diffusion)
    write_table(args.out/'duration_projection.csv',projections)
    findings=[
        {'input':'Flow controller reference','status':'NOMINAL_PRODUCT_CONVENTION_IDENTIFIED_EXPERIMENT_UNCONFIRMED',
         'evidence':'Brooks 5850E manual printed p.7: standard flow uses 0 C and 101.3 kPa (760 Torr); laboratory zero adjustment at 21.1 C is a different quantity.',
         'source':BROOKS,'action':'Run 273.15 K / 101325 Pa alternative; retain original 298.15 K assumption as reference.',
         'helium_molar_flow_ratio_0C_to_25C':298.15/273.15},
        {'input':'Silica temperature coordinates','status':'INTERNALLY_INCONSISTENT_PAPER_INPUTS',
         'evidence':'Figure curve and printed temperatures disagree at deposition scan maxima; two species imply different shifts.',
         'source':PAPER,'action':'Run two explicitly labelled alignment hypotheses; no unique correction selected.'},
        {'input':'Species diffusion','status':'PbI2_CORRELATION_AVAILABLE_OTHER_SPECIES_PROVISIONAL',
         'evidence':'Repository utilities/diffusivity.py derives PbI2-He collision parameters from estimated molecular properties. The C++ curve reproduces that model; this is not an independent gas-diffusion measurement.',
         'source':'utilities/diffusivity.py; utilities/interpolation.py; applications/rhoFixedFlowFoam/PbI2HeDiffusivity.H',
         'action':'Compare PbI2He against constant D, and bracket common D by factors 0.5 and 2. No BiI3-specific coefficient fabricated.'},
        {'input':'Mass accommodation','status':'SURFACE_SPECIFIC_COEFFICIENTS_UNCONFIRMED',
         'evidence':'A primary BiI3 effusion study derives condensation coefficients, but accessible abstract does not provide transferable values for these steel/silica walls.',
         'source':'https://doi.org/10.1016/0021-9614(90)90072-X',
         'action':'Pair accommodation 0.1 and 0.01 are sensitivity assumptions. Do not substitute thermal accommodation coefficients for mass accommodation.'},
        {'input':'BiI thermodynamics','status':'PROVENANCE_AUDITED_INDEPENDENT_COEFFICIENT_CHECK_PENDING',
         'evidence':'Internal gas record has explicit temperature coverage, an absolute Gibbs convention and HERACLES/Barin source labels. The 29 kJ/mol shift is not a second measured record.',
         'source':'https://doi.org/10.1021/j100821a030; https://doi.org/10.1021/ic50060a006; external ThermoFun metadata',
         'action':'Retain both branches as sensitivity scenarios. Identified primary papers require full coefficient/standard-state reconciliation before replacing the baseline.'},
        {'input':'Boat geometry','status':'PARTIAL_PUBLISHED_GEOMETRY',
         'evidence':'Paper specifies about 3 mm boat width and 15 mm length; height, detailed cross-section and exact axial position are not supplied.',
         'source':PAPER,'action':'Shift source by one boat length in either direction; empty-tube carrier still omits boat obstruction.'},
        {'input':'Release history','status':'UNIDENTIFIABLE_FROM_FINAL_SCAN',
         'evidence':'The paper reports a three-hour evaporation window and final deposits; it does not provide a measured species-release time series.',
         'source':PAPER,'action':'Retain three-hour mean release rates; refine saturation 0.2%, 0.316% and 0.562% between existing branches without fitting.'},
        {'input':'Shutdown and cooling','status':'PROTOCOL_KNOWN_THERMAL_HISTORY_UNAVAILABLE',
         'evidence':'Sample withdrawal, furnace shutdown and carrier-flow shutdown occur together before cooldown and scanning. No temperature-versus-time trajectory or complete thermal apparatus is given.',
         'source':PAPER,'action':'Run closed isothermal flow-stop relaxation from 60 to 80 s and conditional three-hour rate projections. Neither reconstructs experimental cooling.'}]
    result={'status':'PHYSICAL_EVIDENCE_AUDITED_WITH_EXPLICIT_REMAINING_GAPS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'findings':findings,
        'thermofun_file_sha256':hashlib.sha256(args.tf.read_bytes()).hexdigest(),'thermodynamic_temperature_coverage':coverage,
        'diffusion_screening':diffusion,'duration_projections':projections,
        'restricted_data_policy':'No external thermodynamic coefficients or evaluated thermodynamic tables included.',
        'experimental_evidence':'Supplied Liu paper only; manufacturer documentation and other primary publications inform model qualification, not new experimental inputs.'}
    (args.out/'qualification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><title>M10e physical input audit</title>
    <style>body{font:16px/1.6 system-ui;max-width:1100px;margin:32px auto;padding:0 20px}table{display:block;overflow-x:auto;border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ccc}</style><body>
    <h1>Physical input qualification and duration projections</h1><p><strong>Evidence was audited; unresolved coefficients and histories remain explicitly provisional.</strong>
    The supplied paper does not determine unique source release, silica coordinate alignment or cooling history. No missing experiment was reconstructed as measured data.</p>'''
    for finding in findings:
        document+='<h2>'+finding['input']+'</h2><p><strong>'+finding['status'].replace('_',' ')+'</strong>. '+finding['evidence']+'</p><p>'+finding['action']+'</p>'
        sources=finding['source'].split('; ')
        document+='<p>Sources: '+', '.join('<a href="'+s+'">'+s+'</a>' if s.startswith('https:') else s for s in sources)+'</p>'
    document+='<h2>PbI2 diffusivity screening</h2>'+tables_html(diffusion,[(k,k.replace('_',' ')) for k in diffusion[0]])
    document+='''<p>The estimated diffusivity changes with temperature and pressure. The comparison runs change only PbI2 when selecting PbI2He.
    Agreement between Python and C++ establishes implementation fidelity, not physical accuracy of collision parameters.</p>
    <h2>Conditional three-hour continuation</h2><p>These profiles continue the 50–60 s net wall rates analytically to the paper's 10800 s exposure.
    They are not native three-hour simulations. Source rate, flow, wall geometry and thermodynamics are assumed constant; no measured release history or cooling is inserted.
    The steady residual gas fraction can bound additional deposition from that gas under these assumptions, but cannot bound redistribution of existing wall material during cooling.</p>'''
    document+=tables_html(projections,[(k,k.replace('_',' ')) for k in
        ['case','species','projected_full_column_wall_mol','peak_x_cm','residual_gas_I_over_3h_supply_if_inventory_steady','status']])
    document+='''<p><a href="qualification.json">Audit and provenance JSON</a> · <a href="duration_projection.csv">Duration projection CSV</a> ·
    <a href="diffusion_screening.csv">Diffusion screening CSV</a> · <a href="duration_profiles/">Projected profiles</a></p></body></html>'''
    (args.out/'report.html').write_text(document)
    artifacts={str(path.relative_to(args.out)):hashlib.sha256(path.read_bytes()).hexdigest()
               for path in sorted(args.out.rglob('*')) if path.is_file() and path.name!='artifact_sha256.json'}
    (args.out/'artifact_sha256.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    print('Audited eight physical input topics and generated',len(projections),'conditional duration profiles.')


if __name__=='__main__':main()
