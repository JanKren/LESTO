#!/usr/bin/env python3
"""Report completed follow-up cases; incomplete results never imply qualification."""
import argparse
import base64
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyse_matrix import profile, moments, relative_l1, tables_html, write_table, ledger_audit
from followup_study import BASES, FULL
from report import analyse


def richardson(coarse, medium, fine, ratio=1.5):
    """Conditional scalar estimate, refusing oscillation and non-decreasing errors."""
    values=[coarse,medium,fine]
    if not all(math.isfinite(v) for v in values):raise ValueError('Nonfinite functional')
    if ratio<=1:raise ValueError('Refinement ratio must exceed one')
    older=coarse-medium;newer=medium-fine
    scale=max(abs(v) for v in values)
    floor=256*np.finfo(float).eps*max(scale,1e-300)
    if min(abs(older),abs(newer))<=floor:
        return {'status':'UNRESOLVED_AT_ROUNDOFF','observed_order':None,'fine_GCI_fraction':None}
    if math.copysign(1,older)!=math.copysign(1,newer):
        return {'status':'OSCILLATORY_NO_EXTRAPOLATION','observed_order':None,'fine_GCI_fraction':None}
    p=math.log(abs(older/newer))/math.log(ratio)
    if p<=0:
        return {'status':'NONDECREASING_NO_EXTRAPOLATION','observed_order':p,'fine_GCI_fraction':None}
    denominator=ratio**p-1
    return {'status':'CONDITIONAL_ESTIMATE_NOT_ASYMPTOTICALLY_CONFIRMED','observed_order':p,
            'extrapolated_value':fine+(fine-medium)/denominator,
            'fine_GCI_fraction':1.25*abs((fine-medium)/fine)/denominator if fine else None,
            'safety_factor':1.25,'refinement_ratio':ratio}


def case_profiles(case,time):
    curves={}
    for species in ['PbI2','BiI3','I']:
        x,amounts,_=profile(case,species,time)
        curves[species]=(x,amounts)
    return curves


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    root=args.root.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    plan=json.loads((root/'study.json').read_text())
    state=json.loads((root/'run_state.json').read_text()) if (root/'run_state.json').exists() else {'status':'NOTRUN','cases':{}}
    completed={};metrics=[];comparisons=[];missing=[];audits={};gci=[]
    baseline={}
    for entry in plan['cases']:
        name=entry['case'];case=root/'cases'/name
        if state['cases'].get(name,{}).get('status')!='PASS':
            missing.append({'case':name,'status':state['cases'].get(name,{}).get('status','NOTRUN')});continue
        metadata=json.loads((case/'benchmark.json').read_text());end=metadata['end_time_s']
        report,peaks=analyse(case);assert report['status']=='PASS'
        destination=out/'profiles'/name;destination.mkdir(parents=True,exist_ok=True)
        for filename in ['PbI2_g_native.dat','BiI3_g_native.dat','element_I_native.dat']:
            shutil.copyfile(case/'postProcessing/axialProfiles'/f'{end:g}'/filename,destination/filename)
        shutil.copyfile(case/'iodine-comparison.csv',destination/'iodine-comparison.csv')
        curves=case_profiles(case,end);completed[name]=curves
        base=entry['base_case']
        if base not in baseline:baseline[base]=case_profiles(FULL/'cases'/base,60)
        peak_by_species={r['species']:r for r in peaks}
        metrics.append({'case':name,'category':entry['category'],'base_case':base,
            'cells':metadata['carrier']['cells'],'dt_s':metadata['dt_s'],'end_time_s':end,
            'BiI3_iodide_I_pct':100*report['BiI3_fraction_of_modelled_iodide_wall_iodine'],
            'raw_profile_L1_percentage_points':report['iodine_profile_L1_percent_over_digitized_scan'],
            'PbI2_peak_C':peak_by_species['PbI2']['peak_temperature_K']-273.15 if peak_by_species['PbI2']['peak_temperature_K'] else None,
            'BiI3_peak_C':peak_by_species['BiI3']['peak_temperature_K']-273.15 if peak_by_species['BiI3']['peak_temperature_K'] else None})
        if end==60:audits[name]=ledger_audit(case)
        else:audits[name]={'scope':'Restarted shutdown; existing original element ledger values retained. See numerical-report.json.'}
        for species,(x,amounts) in curves.items():
            xb,reference=baseline[base][species];assert np.allclose(x,xb)
            total=amounts.sum();reference_total=reference.sum()
            comparisons.append({'case':name,'category':entry['category'],'species':species,
                'profile_L1_relative_to_base':relative_l1(amounts,reference),
                'normalised_shape_L1':relative_l1(amounts/total,reference/reference_total) if total and reference_total else None,
                'deposit_relative_change':(total-reference_total)/reference_total if reference_total else None,
                'centroid_shift_cm':100*(moments(x,amounts)['centroid_m']-moments(x,reference)['centroid_m']) if total and reference_total else None,
                'comparison_scope':'80 s isothermal shutdown vs original 60 s inventory' if end!=60 else '60 s comparison'})
    for base in BASES:
        coarse=base+'__grid256x32';fine=base+'__grid576x72'
        if coarse not in completed or fine not in completed:continue
        for species in ['PbI2','BiI3','I']:
            triplet=[completed[coarse][species],baseline[base][species],completed[fine][species]]
            for functional in ['total','centroid_m','sigma_m']:
                values=[moments(x,a)[functional] for x,a in triplet]
                if any(v is None for v in values):continue
                gci.append({'base_case':base,'species':species,'functional':functional,'coarse':values[0],
                            'medium':values[1],'fine':values[2],**richardson(*values)})
    write_table(out/'case_metrics.csv',metrics);write_table(out/'comparisons.csv',comparisons)
    if not metrics:(out/'case_metrics.csv').write_text('case,category,status\n')
    if not comparisons:(out/'comparisons.csv').write_text('case,category,species,status\n')
    write_json=lambda path,value:path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    write_json(out/'grid_estimates.json',gci)
    summary={'execution_status':state['status'],'completed_cases':len(metrics),'planned_cases':len(plan['cases']),
             'missing_or_failed':missing,'case_metrics':metrics,'comparisons':comparisons,'conditional_grid_estimates':gci,
             'material_audits':audits,'solver_sha256':plan['solver_sha256'],
             'study_manifest_sha256':hashlib.sha256((root/'study.json').read_bytes()).hexdigest(),
             'scope':'Numerical and sensitivity evidence; no parameter fitting. Thermal cooling history and independent BiI coefficients remain unqualified.'}
    write_json(out/'analysis.json',summary)
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><title>M10e follow-up studies</title>
    <style>body{font:16px/1.6 system-ui;max-width:1150px;margin:32px auto;padding:0 20px}table{display:block;overflow-x:auto;border-collapse:collapse;font-size:12px}td,th{padding:6px;border-bottom:1px solid #ccc}img{max-width:100%}</style><body>
    <h1>M10e refinement and physical sensitivity studies</h1>'''
    document+=f'<p><strong>{len(metrics)} of {len(plan["cases"])} follow-up cases complete.</strong> Execution status: {html.escape(state["status"])}.</p>'
    document+='''<p>Numerical cases use the same corrected private solver, source rates, thermodynamic tables and 60 s window.
    Time-step checks compare 2 versus 1 ms on 18,432 cells. The three-grid sequence is 256×32, 384×48 and 576×72 cells, at 2 ms,
    with a constant axial/radial refinement ratio of 1.5. The 384×48 results are the completed original matrix.
    Physical scenarios vary individual assumptions; they are not experimentally qualified alternatives.</p>'''
    document+='<h2>Case metrics</h2>'+tables_html(metrics,[(k,k.replace('_',' ')) for k in
        ['case','category','cells','dt_s','BiI3_iodide_I_pct','raw_profile_L1_percentage_points','PbI2_peak_C','BiI3_peak_C']])
    document+='<h2>Changes from the original case</h2>'+tables_html(comparisons,[(k,k.replace('_',' ')) for k in
        ['case','species','profile_L1_relative_to_base','normalised_shape_L1','deposit_relative_change','centroid_shift_cm']])
    if comparisons:
        fig,ax=plt.subplots(figsize=(12,max(5,len(metrics)*.24)),layout='constrained')
        iodine=[r for r in comparisons if r['species']=='I' and r['profile_L1_relative_to_base'] is not None]
        ax.barh(range(len(iodine)),[100*r['profile_L1_relative_to_base'] for r in iodine])
        ax.set(yticks=range(len(iodine)),yticklabels=[r['case'] for r in iodine],xlabel='Iodine deposit profile L1 change from original (%)')
        ax.tick_params(axis='y',labelsize=6);ax.invert_yaxis()
        fig.savefig(out/'profile_changes.png',dpi=150,bbox_inches='tight');plt.close(fig)
        document+='<img src="data:image/png;base64,'+base64.b64encode((out/'profile_changes.png').read_bytes()).decode()+'" alt="Follow-up profile changes">'
    numerical=any(r['category']=='numerical' for r in metrics)
    if numerical:
        fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
        for ax,base in zip(axes.flat,BASES):
            if base not in baseline:continue
            for label,curves in [('Original 384×48, 2 ms',baseline[base]),
                ('256×32, 2 ms',completed.get(base+'__grid256x32')),
                ('576×72, 2 ms',completed.get(base+'__grid576x72')),
                ('384×48, 1 ms',completed.get(base+'__dt1ms'))]:
                if curves is None:continue
                x,amounts=curves['I'];ax.plot(100*x,100*amounts/amounts.sum(),label=label)
            ax.set(title=base.replace('LBE-I_SS_II','Steel'),xlim=(28,67),xlabel='Column coordinate (cm)',ylabel='Normalised wall iodine per 1 cm bin (%)')
            ax.legend(fontsize=8)
        fig.suptitle('Isolated time-step and systematic mesh comparisons; same solver and source rates')
        fig.savefig(out/'refinement_profiles.png',dpi=160,bbox_inches='tight');plt.close(fig)
        document+='<img src="data:image/png;base64,'+base64.b64encode((out/'refinement_profiles.png').read_bytes()).decode()+'" alt="Refinement profile overlays">'
    document+='''<h2>Conditional grid estimates</h2><p>Scalar Richardson/GCI calculations follow the
    <a href="https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html">NASA grid-convergence procedure</a>.
    Oscillatory, non-decreasing and round-off-scale sequences receive no extrapolation.
    Three grids alone do not establish an asymptotic regime for the complete distribution. A stable total amount cannot establish profile convergence.
    Time-discretisation, source selection and non-linear wall/chemical coupling remain relevant.</p>'''
    document+=tables_html(gci,[(k,k.replace('_',' ')) for k in ['base_case','species','functional','status','observed_order','fine_GCI_fraction']])
    document+='''<h2>Physical scope</h2><p>Silica shifts are explicit hypotheses derived separately from the PbI2 and BiI3 curve/label discrepancy.
    They must not replace the published figure as a corrected measurement. The 0 °C flow alternative tests the nominal convention in the Brooks 5850E manual;
    experiment-specific calibration is still unknown. Iodide accommodation changes leave metal-wall-channel accommodation at unity.
    Boat-position changes omit the boat obstruction, whose complete geometry is unavailable. PbI2He changes only PbI2 transport.
    All other trace diffusivities remain provisional.</p>
    <p>The four shutdown cases restart the original 60 s state with release expired, zero carrier flux and zero normal species flux at both ends,
    then relax isothermally until 80 s. Gas inventory and deposited state are preserved. This tests post-flow redistribution at fixed temperature,
    not the experimental cooling trajectory or a three-hour exposure. Cooling cannot be uniquely reconstructed from the paper.</p>
    <p><a href="analysis.json">Machine-readable analysis</a> · <a href="case_metrics.csv">Case metrics CSV</a> ·
    <a href="comparisons.csv">Comparison CSV</a> · <a href="grid_estimates.json">Conditional grid estimates</a></p></body></html>'''
    (out/'report.html').write_text(document)
    artifact_hashes={str(path.relative_to(out)):hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in sorted(out.rglob('*')) if path.is_file() and path.name!='artifact_sha256.json'}
    write_json(out/'artifact_sha256.json',artifact_hashes)
    print(f'Reported {len(metrics)}/{len(plan["cases"])} completed cases',flush=True)


if __name__=='__main__':main()
