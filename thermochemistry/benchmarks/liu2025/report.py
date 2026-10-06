#!/usr/bin/env python3
"""Report numerical gates and provisional peaks, iodine profiles and splits."""
import argparse
import json
import math
import re
from pathlib import Path

from prepare import HERE, read_csv, write_csv, load_references


def rows(path):
    text = path.read_text()
    columns = re.search(r'^# columns (.+)$',text,re.M)
    if not columns:
        raise ValueError('Missing column names: '+str(path))
    keys = columns[1].split()
    result = []
    for line in text.splitlines():
        if line and not line.startswith('#'):
            values = list(map(float,line.split()))
            if len(values)!=len(keys) or not all(math.isfinite(x) for x in values):
                raise ValueError('Malformed/nonfinite result: '+str(path))
            result.append(dict(zip(keys,values)))
    if not result:
        raise ValueError('Empty result: '+str(path))
    return result


def numerical_gates(case, end_time):
    report = {}
    for element in ['Pb','Bi','I']:
        paths = sorted((case/'postProcessing/phaseChangeElements').glob('*/'+element+'.dat'))
        if not paths:
            raise ValueError('Missing element ledger '+element)
        ledger = sorted((r for p in paths for r in rows(p)),key=lambda r:r['time'])
        if any(r['reference']<0 for r in ledger):
            raise ValueError('Negative element reference: '+str(case))
        if abs(ledger[-1]['time']-end_time)>1e-7*max(1,end_time):
            raise ValueError('Incomplete run: '+str(case))
        closure = max(abs(r['closure'])/max(r['reference'],1e-300) for r in ledger)
        defect = max(abs(r['solverDefect'])/max(r['reference'],1e-300) for r in ledger)
        report[element] = dict(max_relative_closure=closure, max_relative_solver_defect=defect,
                               status='PASS' if closure<=1e-12 and defect<=1e-8 else 'FAIL')
    return report


def analyse(case):
    manifest=json.loads((case/'benchmark.json').read_text())
    gates=numerical_gates(case,manifest['end_time_s'])
    times=sorted((p for p in (case/'postProcessing/axialProfiles').iterdir() if p.is_dir()),key=lambda p:float(p.name))
    if not times or abs(float(times[-1].name)-manifest['end_time_s'])>1e-7*max(1,manifest['end_time_s']):
        raise ValueError('Missing final profiles: '+str(case))
    end=times[-1]
    curves={name:rows(end/(name+'_native.dat')) for name in
            ['PbI2_g','BiI3_g','condensate_Pb_s','condensate_Bi_s','element_I']}
    older=None
    interval=None
    if len(times)>1:
        older=times[-2]
        interval=float(end.name)-float(older.name)
    comparison=[]
    rate_changes={}
    spatial_changes={}
    for species,name in [('PbI2','PbI2_g'),('BiI3','BiI3_g'),('Pb','condensate_Pb_s'),('Bi','condensate_Bi_s')]:
        curve=curves[name]
        target=next((r for r in manifest['targets'] if r['species']==species),None)
        deposited=math.fsum(r['wall']*(r['x_hi']-r['x_lo']) for r in curve)
        peak=max(curve,key=lambda r:r['wall'])
        value=peak['Twall'] if deposited>0 else None
        printed=float(target['measured_peak_C'])+273.15 if target else None
        historical_equilibrium=float(target['published_eq3_C'])+273.15 if target else None
        change=None
        if older:
            prev=rows(older/(name+'_native.dat'))
            previous=math.fsum(r['wall']*(r['x_hi']-r['x_lo']) for r in prev)
            change=(deposited-previous)/interval/manifest['carrier']['wedge_fraction']
            if len(times)>=3:
                earliest=rows(times[-3]/(name+'_native.dat'))
                first=math.fsum(r['wall']*(r['x_hi']-r['x_lo']) for r in earliest)
                previous_rate=(previous-first)/(float(older.name)-float(times[-3].name))/manifest['carrier']['wedge_fraction']
                rate_changes[species]=abs(change-previous_rate)/max(abs(change),abs(previous_rate),1e-300)
                first_interval=float(older.name)-float(times[-3].name)
                rates_now=[(a['wall']-b['wall'])/interval for a,b in zip(curve,prev)]
                rates_before=[(a['wall']-b['wall'])/first_interval for a,b in zip(prev,earliest)]
                weights=[r['x_hi']-r['x_lo'] for r in curve]
                spatial_changes[species]=math.fsum(abs(a-b)*w for a,b,w in zip(rates_now,rates_before,weights))/max(
                    math.fsum(abs(a)*w for a,w in zip(rates_now,weights)),
                    math.fsum(abs(a)*w for a,w in zip(rates_before,weights)),1e-300)
        difference=value-printed if printed and value else None
        science_target=20 if species=='PbI2' else 40 if species=='BiI3' else None
        comparison.append(dict(case=case.name,experiment=manifest['experiment'],species=species,
            status='PROVISIONAL_NUMERICS' if deposited>0 else 'NO_DEPOSIT',
            peak_temperature_K=value,peak_x_m=peak['x'] if deposited>0 else None,
            measured_peak_K=printed,difference_K=difference,
            published_constant_HS_eq3_K=historical_equilibrium,
            peak_minus_published_eq3_K=value-historical_equilibrium if value and historical_equilibrium else None,
            within_experimental_uncertainty=abs(difference)<=float(target['uncertainty_K']) if difference is not None else None,
            science_target_K=science_target,
            within_science_target=abs(difference)<=science_target if difference is not None and science_target else None,
            full_column_deposit_mol=deposited/manifest['carrier']['wedge_fraction'],
            last_window_net_deposition_mol_s=change))
    pb=next(r['full_column_deposit_mol'] for r in comparison if r['species']=='PbI2')
    bi=next(r['full_column_deposit_mol'] for r in comparison if r['species']=='BiI3')
    split=3*bi/(2*pb+3*bi) if 2*pb+3*bi>0 else None
    iodine=curves['element_I']
    total=math.fsum(r['wall']*(r['x_hi']-r['x_lo']) for r in iodine)
    gamma_path=HERE/'digitized'/(manifest['experiment']+'_gamma.csv')
    measured=read_csv(gamma_path) if gamma_path.exists() else []
    profile=[]
    for r in iodine:
        observed=next((float(v['relative_yield_percent_per_cm']) for v in measured
                       if abs(float(v['x_m'])-r['x'])<1e-6),None)
        predicted=(100-manifest['excluded_KI_iodine_yield_pct'])*r['wall']*(r['x_hi']-r['x_lo'])/total if total>0 else None
        profile.append(dict(x_m=r['x'],full_column_iodine_wall_mol_m=r['wall']/manifest['carrier']['wedge_fraction'],
                            predicted_relative_yield_percent=predicted,figure_relative_yield_percent=observed))
    write_csv(case/'iodine-comparison.csv',profile,list(profile[0]))
    log=(case/'log.solver').read_text()
    failures=[int(v) for v in re.findall(r'Re-speciation:.*?failures (\d+)',log)]
    if not failures:
        raise ValueError('Missing re-speciation diagnostics: '+str(case))
    result=dict(case=case.name,status='PASS' if all(r['status']=='PASS' for r in gates.values()) and max(failures)==0 else 'FAIL',
        scientific_validation_status='PROVISIONAL_ONLY',element_gates=gates,
        max_respeciation_failures=max(failures),
        BiI3_fraction_of_modelled_iodide_wall_iodine=split,
        excluded_KI_yield_pct=manifest['excluded_KI_iodine_yield_pct'],
        Q1_split_is_conditioned=manifest['source_mode']=='Q1',
        iodine_profile_L1_percent_over_digitized_scan=sum(abs(r['predicted_relative_yield_percent']-r['figure_relative_yield_percent'])
            for r in profile if r['predicted_relative_yield_percent'] is not None and r['figure_relative_yield_percent'] is not None),
        iodine_profile_comparison_scope='Digitized scan range; observed includes KI, model excludes its stated yield; no shape fitting.',
        last_deposition_rate_window_s=interval,
        relative_change_in_last_two_net_deposition_rates=rate_changes,
        spatial_L1_change_in_last_two_deposition_rate_profiles=spatial_changes,
        integrated_rate_stationarity='REPORTED' if rate_changes else 'NOTRUN',
        quasi_steady_qualification=('PRELIMINARY_WITHIN_1_PERCENT' if rate_changes and
             max(rate_changes.values())<=.01 and max(spatial_changes.values())<=.01
             else 'NOT_QUALIFIED'),
        note='End-of-window wall inventory; no cooling/shutdown or clean-He flush simulated.')
    (case/'numerical-report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result,comparison


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--plot',action='store_true',help='Write a standalone preliminary iodine-profile plot (matplotlib)')
    args=parser.parse_args()
    matrix=json.loads((args.root/'matrix.json').read_text())
    results=[]; comparisons=[]
    _,targets,_=load_references()
    def missing_peaks(entry,status):
        values=[]
        for species in ['PbI2','BiI3','Pb','Bi']:
            target=next((r for r in targets if r['experiment']==entry['experiment'] and r['species']==species),None)
            values.append(dict(case=entry['case'],experiment=entry['experiment'],species=species,status=status,
                peak_temperature_K=None,peak_x_m=None,
                measured_peak_K=float(target['measured_peak_C'])+273.15 if target else None,difference_K=None,
                published_constant_HS_eq3_K=float(target['published_eq3_C'])+273.15 if target else None,
                peak_minus_published_eq3_K=None,within_experimental_uncertainty=None,
                science_target_K=20 if species=='PbI2' else 40 if species=='BiI3' else None,
                within_science_target=None,full_column_deposit_mol=None,last_window_net_deposition_mol_s=None))
        return values
    for entry in matrix:
        case=args.root/entry['case']
        if not (case/'postProcessing').exists():
            results.append(dict(case=entry['case'],status='NOTRUN',reason=entry['reason'] or 'No CFD output'))
            comparisons+=missing_peaks(entry,'NOTRUN')
        else:
            try:
                result, peaks=analyse(case)
                results.append(result);comparisons+=peaks
            except (ValueError,FileNotFoundError) as error:
                results.append(dict(case=entry['case'],status='FAIL',reason=str(error)))
                comparisons+=missing_peaks(entry,'FAIL')
    (args.root/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
    if comparisons:
        write_csv(args.root/'peaks.csv',comparisons,list(comparisons[0]))
    if args.plot:
        plot(args.root,matrix)
    print(json.dumps({s:sum(r['status']==s for r in results) for s in ['PASS','FAIL','NOTRUN']}))
    if any(r['status']=='FAIL' for r in results):
        raise SystemExit(1)


def plot(root,matrix):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=['PbI2_SS','LBE-I_SS_II','LBE-I_SiO2_I']
    figure,axes=plt.subplots(3,1,figsize=(9,10),constrained_layout=True)
    colors={.001:'tab:blue',.01:'tab:orange',.1:'tab:green'}
    for ax,name in zip(axes,names):
        observed=read_csv(HERE/'digitized'/(name+'_gamma.csv'))
        ax.bar([float(r['x_m']) for r in observed],
               [float(r['relative_yield_percent_per_cm']) for r in observed],
               width=.01,color='0.85',edgecolor='0.5',label='Published iodine scan (digitised)')
        for entry in matrix:
            path=root/entry['case']/'iodine-comparison.csv'
            if entry['experiment']!=name or not path.exists():
                continue
            curve=read_csv(path)
            xs=[float(r['x_m']) for r in curve]
            ys=[float(r['predicted_relative_yield_percent']) if r['predicted_relative_yield_percent'] else 0 for r in curve]
            q1=entry['source_mode']=='Q1'
            label='Q1 prescribed split/rate' if q1 else f"Q2 {entry['saturation']:g}, {entry['data_variant']}"
            ax.plot(xs,ys,label=label,color='black' if q1 else colors[entry['saturation']],
                    linestyle='--' if entry['data_variant']=='shift29k' else '-')
        ax.set(xlim=(.25,.7),xlabel='Assumed shared column coordinate (m)',
               ylabel='Iodine yield per 1 cm bin (%)',title=name)
        ax.legend(fontsize=7,ncol=2)
    manifest=next((json.loads((root/r['case']/'benchmark.json').read_text()) for r in matrix
                   if (root/r['case']/'benchmark.json').exists()),None)
    scope=(f"{manifest['carrier']['cells']} cells, {manifest['end_time_s']:g} s, {manifest['dt_s']*1000:g} ms"
           if manifest else 'No CFD cases')
    figure.suptitle('Preliminary Liu (2025) comparison: '+scope+'\n'
                   'Assumed sources/transport; silica Figure 7C temperature labels disagree with its curve\n'
                   'Experimental data: DOI 10.1007/s10967-025-10246-4 (CC BY 4.0)',fontsize=10)
    figure.savefig(root/'preliminary_profiles.png',dpi=150)
    plt.close(figure)


if __name__=='__main__':
    main()
