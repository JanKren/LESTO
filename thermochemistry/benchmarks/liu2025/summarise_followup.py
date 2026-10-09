#!/usr/bin/env python3
"""Interpret a saved completed-case report without writing to CFD case folders."""
import argparse
import base64
from collections import Counter
import datetime
import hashlib
import html
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from analyse_matrix import experimental_shape, tables_html, write_table
from followup_study import BASES, FULL
from report import rows


LABELS = {
    'PbI2_SS_Q1_baseline': 'Pure PbI2, steel',
    'LBE-I_SS_II_Q2_baseline_sat0.001': 'LBE, steel, 0.1% source saturation',
    'LBE-I_SS_II_Q2_baseline_sat0.01': 'LBE, steel, 1% source saturation',
    'LBE-I_SS_II_Q2_shift29k_sat0.01': 'LBE, steel, 1%, BiI +29 kJ/mol',
    'LBE-I_SiO2_I_Q2_baseline_sat0.001': 'LBE, silica, 0.1% source saturation',
}
NATIVE = {'I': 'element_I_native.dat', 'PbI2': 'PbI2_g_native.dat',
          'BiI3': 'BiI3_g_native.dat'}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def peak_records(folder, name, base, category, targets, positions):
    output = []
    for species in ['PbI2', 'BiI3']:
        data = rows(folder / NATIVE[species])
        total = sum(r['wall'] * (r['x_hi'] - r['x_lo']) for r in data)
        largest = sorted(data, key=lambda r: r['wall'], reverse=True)
        peak = largest[0] if total > 0 else None
        target = next((r for r in targets if r['species'] == species), None)
        measured = float(target['measured_peak_C']) if target else None
        uncertainty = float(target['uncertainty_K']) if target else None
        temperature = peak['Twall'] - 273.15 if peak else None
        x = 100 * peak['x'] if peak else None
        published_x = positions.get((target['experiment'], species)) if target else None
        output.append({
            'case': name, 'base_case': base, 'category': category, 'species': species,
            'model_peak_x_cm': x, 'published_peak_x_cm': published_x,
            'position_error_cm': x - published_x if x is not None and published_x is not None else None,
            'model_peak_C': temperature, 'published_peak_C': measured,
            'temperature_error_K': temperature - measured if temperature is not None and measured is not None else None,
            'published_uncertainty_K': uncertainty,
            'second_largest_bin_over_peak': largest[1]['wall'] / peak['wall'] if peak else None,
            'second_largest_bin_x_cm': 100 * largest[1]['x'] if peak else None,
            'scope': 'Maximum wall-density bin; native scan spacing is 1 cm. No sub-bin peak fit.',
        })
    return output


def scan_curve(folder):
    data = rows_from_comparison(folder / 'iodine-comparison.csv')
    x = np.array([float(r['x_m']) for r in data]) * 100
    model = np.array([float(r['predicted_relative_yield_percent']) for r in data])
    observed = np.array([float(r['figure_relative_yield_percent']) for r in data])
    return x, 100 * model / model.sum(), 100 * observed / observed.sum()


def rows_from_comparison(path):
    import csv
    with path.open() as stream:
        return [r for r in csv.DictReader(stream) if r['figure_relative_yield_percent']]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--baseline-root', type=Path, default=FULL)
    args = parser.parse_args()
    root, out, full = args.root.resolve(), args.report.resolve(), args.baseline_root.resolve()
    analysis = json.loads((out / 'analysis.json').read_text())
    plan = json.loads((root / 'study.json').read_text())
    entries = {r['case']: r for r in plan['cases']}
    metrics = analysis['case_metrics']
    assert len(metrics) == analysis['completed_cases']
    positions = {(r['experiment'], r['species']): r['figure_x_cm'] for r in
                 json.loads((full.parent.parent / 'full_analysis/analysis.json').read_text())['q1_peaks']}
    bases = list(dict.fromkeys(r['base_case'] for r in metrics))
    agreement, peaks, reports, inputs = [], [], {}, {}
    folders = {}
    for base in bases:
        source = full / 'cases' / base
        destination = out / 'reference_profiles' / base
        destination.mkdir(parents=True, exist_ok=True)
        for filename in NATIVE.values():
            shutil.copyfile(source / 'postProcessing/axialProfiles/60' / filename, destination / filename)
        for filename in ['iodine-comparison.csv', 'numerical-report.json', 'benchmark.json']:
            shutil.copyfile(source / filename, destination / filename)
        folders[base] = destination
    for name, base, category, folder, source in [
        (base, base, 'original_baseline', folders[base], full / 'cases' / base) for base in bases
    ] + [(r['case'], r['base_case'], r['category'], out / 'profiles' / r['case'],
          root / 'cases' / r['case']) for r in metrics]:
        metadata = json.loads((source / 'benchmark.json').read_text())
        numerical = json.loads((source / 'numerical-report.json').read_text())
        assert numerical['status'] == 'PASS'
        reports[name] = numerical
        inputs[name] = {'manifest_sha256': hashlib.sha256((source / 'benchmark.json').read_bytes()).hexdigest(),
                        'source_mode': metadata['source_mode'],
                        'source_split_conditioned': numerical['Q1_split_is_conditioned']}
        for filename in ['numerical-report.json', 'benchmark.json']:
            if source != folder:
                shutil.copyfile(source / filename, folder / filename)
        shape = experimental_shape(folder)
        iodide_yields = {r['species']: float(r['iodine_yield_pct']) for r in metadata['targets']
                         if r['species'] in ['PbI2', 'BiI3'] and r['iodine_yield_pct']}
        measured_split = (100 * iodide_yields.get('BiI3', 0) / sum(iodide_yields.values())
                          if iodide_yields else None)
        model_split = 100 * numerical['BiI3_fraction_of_modelled_iodide_wall_iodine']
        agreement.append({
            'case': name, 'base_case': base, 'category': category,
            'raw_L1_percentage_points': shape['raw_L1_percentage_points'],
            'normalised_common_scan_overlap_pct': 100 * shape['normalised_scan_overlap'],
            'normalised_common_scan_Wasserstein_cm': 100 * shape['normalised_scan_Wasserstein_m'],
            'normalised_common_scan_centroid_error_cm': 100 * shape['normalised_scan_centroid_difference_m'],
            'BiI3_fraction_of_iodide_I_pct': model_split,
            'published_BiI3_fraction_of_iodide_I_pct': measured_split,
            'BiI3_split_error_percentage_points': model_split - measured_split if measured_split is not None else None,
            'Q1_source_split_conditioned': numerical['Q1_split_is_conditioned'],
            'observed_scan_sum_pct': shape['observed_scan_sum_pct'],
            'model_scan_sum_pct': shape['model_scan_sum_pct'],
        })
        peaks += peak_records(folder, name, base, category, metadata['targets'], positions)
    iodine = [r for r in analysis['comparisons'] if r['species'] == 'I']
    by_name = {r['case']: r for r in iodine}
    numerical_rows = []
    for base in BASES:
        selected = {suffix: by_name.get(base + '__' + suffix) for suffix in
                    ['dt1ms', 'grid256x32', 'grid576x72']}
        if not all(selected.values()):
            continue
        numerical_rows.append({
            'base_case': base,
            'dt_2_to_1_ms_I_profile_L1_pct': 100 * selected['dt1ms']['profile_L1_relative_to_base'],
            'coarse_vs_medium_I_profile_L1_pct': 100 * selected['grid256x32']['profile_L1_relative_to_base'],
            'fine_vs_medium_I_profile_L1_pct': 100 * selected['grid576x72']['profile_L1_relative_to_base'],
            'fine_vs_medium_total_I_change_pct': 100 * selected['grid576x72']['deposit_relative_change'],
            'fine_vs_medium_I_centroid_shift_cm': selected['grid576x72']['centroid_shift_cm'],
        })
    sensitivity = [{
        'case': r['case'], 'base_case': entries[r['case']]['base_case'], 'category': r['category'],
        'I_profile_L1_change_pct': 100 * r['profile_L1_relative_to_base'],
        'I_normalised_shape_L1_change_pct': 100 * r['normalised_shape_L1'] if r['normalised_shape_L1'] is not None else None,
        'I_total_change_pct': 100 * r['deposit_relative_change'] if r['deposit_relative_change'] is not None else None,
        'I_centroid_shift_cm': r['centroid_shift_cm'],
    } for r in iodine if r['category'] != 'numerical']
    audit_fields = ['max_accounted_closure_over_source', 'max_raw_material_residue_over_source',
                    'max_net_clamped_over_source', 'max_solver_defect_over_source']
    conservation = {k: max(r[k] for a in analysis['material_audits'].values() for r in a.values() if k in r)
                    for k in audit_fields}
    time_changes = [r['dt_2_to_1_ms_I_profile_L1_pct'] for r in numerical_rows]
    mesh_changes = [r['fine_vs_medium_I_profile_L1_pct'] for r in numerical_rows]
    category_max = {category: max(r['I_profile_L1_change_pct'] for r in sensitivity if r['category'] == category)
                    for category in {r['category'] for r in sensitivity}}
    findings = []
    if time_changes:
        findings.append(f'Halving the time step from 2 to 1 ms changes the total iodine deposit profile by '
                        f'{min(time_changes):.4f}–{max(time_changes):.3f}% across the four numerical cases. '
                        'This is a two-step comparison, not a temporal extrapolation.')
    if mesh_changes:
        findings.append(f'Refining 18,432 to 41,472 cells changes the iodine deposit profile by '
                        f'{min(mesh_changes):.2f}–{max(mesh_changes):.2f}%. '
                        'The complete distribution is not established as mesh independent.')
    if category_max:
        findings.append('Largest completed iodine-profile changes from the original case: ' +
                        '; '.join(f'{key}: {value:.2f}%' for key, value in sorted(category_max.items())) +
                        '. These are scenario effects, not measured parameter uncertainties.')
    pure_base = next((r for r in peaks if r['case'] == BASES[0] and r['species'] == 'PbI2'), None)
    pure_fine = next((r for r in peaks if r['case'] == BASES[0] + '__grid576x72' and r['species'] == 'PbI2'), None)
    if pure_base and pure_fine:
        findings.append('The 1 cm scan bins can change which bin is identified as the peak. '
                        f'For pure PbI2, the original {pure_base["second_largest_bin_x_cm"]:.0f} and '
                        f'{pure_base["model_peak_x_cm"]:.0f} cm bins differ in height by only '
                        f'{100 * (1 - pure_base["second_largest_bin_over_peak"]):.2f}%; '
                        f'the fine mesh selects {pure_fine["model_peak_x_cm"]:.0f} cm, changing the '
                        f'reported peak temperature from about {pure_base["model_peak_C"]:.0f} to '
                        f'{pure_fine["model_peak_C"]:.0f} °C. Use full profiles alongside peak temperatures.')
    fine_split = {r['base_case']: r['BiI3_fraction_of_iodide_I_pct'] for r in agreement
                  if r['case'].endswith('__grid576x72')}
    if BASES[2] in fine_split and BASES[3] in fine_split:
        findings.append('The large chemistry sensitivity survives refinement: at 1% source saturation, '
                        f'fine-mesh BiI3 carries {fine_split[BASES[2]]:.2f}% of deposited iodide iodine '
                        f'with baseline thermodynamics versus {fine_split[BASES[3]]:.2f}% with the '
                        'BiI +29 kJ/mol surrogate. This surrogate is not an independently measured thermodynamic dataset.')
    findings.append('Grid estimates are conditional scalar estimates. Oscillatory and non-decreasing sequences '
                    'receive no extrapolation; stable total deposition does not establish profile convergence.')
    pending_categories = Counter(entries[r['case']]['category'] for r in analysis['missing_or_failed'])
    findings.append('Results still pending in this snapshot: ' +
                    ', '.join(f'{k}: {v}' for k, v in sorted(pending_categories.items())) + '.')
    summary = {
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'Interim completed-case snapshot; no fitting, no writes to CFD inputs or running case outputs.',
        'completed_cases': analysis['completed_cases'], 'planned_cases': analysis['planned_cases'],
        'completed_by_category': dict(Counter(r['category'] for r in metrics)),
        'pending_by_category': dict(pending_categories), 'findings': findings,
        'conservation_maxima': conservation,
        'max_respeciation_failures': max(reports[r['case']]['max_respeciation_failures'] for r in metrics),
        'conditional_grid_estimate_status_counts': dict(Counter(r['status'] for r in analysis['conditional_grid_estimates'])),
        'numerical_comparison': numerical_rows, 'completed_sensitivities': sensitivity,
        'experimental_comparison': agreement, 'peaks': peaks, 'case_inputs': inputs,
        'source_study_manifest_sha256': analysis['study_manifest_sha256'],
        'postprocessing_source_sha256': {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__), Path(__file__).with_name('report_followup.py'),
                      Path(__file__).with_name('analyse_matrix.py'), Path(__file__).with_name('report.py')]},
        'raw_profile_metric': '100 * sum(abs(I_variant - I_baseline)) / sum(I_baseline), in identical 1 cm bins.',
        'experimental_metric': 'Raw published yields retained for primary L1; overlap normalises each curve separately over the common scan. Observed scan retains KI.',
    }
    write_json(out / 'interpretation.json', summary)
    write_table(out / 'numerical_comparison.csv', numerical_rows)
    write_table(out / 'experimental_comparison.csv', agreement)
    write_table(out / 'sensitivity_summary.csv', sensitivity)
    write_table(out / 'peak_positions.csv', peaks)
    figure_dir = out / 'figures'
    figure_dir.mkdir(exist_ok=True)
    with PdfPages(out / 'postprocessing_figures.pdf') as pdf:
        if numerical_rows:
            fig, ax = plt.subplots(figsize=(11, 5), layout='constrained')
            x = np.arange(len(numerical_rows))
            for offset, key, label in [(-.25, 'dt_2_to_1_ms_I_profile_L1_pct', '2 to 1 ms'),
                                       (0, 'coarse_vs_medium_I_profile_L1_pct', '8,192 vs 18,432 cells'),
                                       (.25, 'fine_vs_medium_I_profile_L1_pct', '41,472 vs 18,432 cells')]:
                ax.bar(x + offset, [r[key] for r in numerical_rows], width=.24, label=label)
            ax.set(yscale='log', ylabel='Iodine deposit profile L1 change (%)',
                   xticks=x, xticklabels=[LABELS[r['base_case']].replace(', ', '\n', 1) for r in numerical_rows])
            ax.set_title('Time-step effects are small; spatial profile sensitivity remains')
            ax.set_ylim(top=4 * max(r[k] for r in numerical_rows for k in
                                   ['coarse_vs_medium_I_profile_L1_pct', 'fine_vs_medium_I_profile_L1_pct']))
            ax.grid(axis='y', alpha=.2); ax.legend(fontsize=9, loc='upper center', ncol=3)
            fig.savefig(figure_dir / 'numerical_comparison.png', dpi=150)
            fig.savefig(figure_dir / 'numerical_comparison.pdf'); pdf.savefig(fig); plt.close(fig)
        fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
        for ax, base in zip(axes.flat, BASES):
            if base not in folders:
                ax.set_visible(False); continue
            for suffix, label in [('', '18,432 cells, 2 ms'), ('grid256x32', '8,192 cells, 2 ms'),
                                  ('grid576x72', '41,472 cells, 2 ms'), ('dt1ms', '18,432 cells, 1 ms')]:
                name = base + '__' + suffix if suffix else base
                folder = out / 'profiles' / name if suffix else folders[base]
                if not folder.exists(): continue
                x, model, observed = scan_curve(folder)
                ax.plot(x, model, label=label)
            ax.fill_between(x, observed, color='.8', step='mid', label='Published gamma scan')
            ax.set(title=LABELS[base], xlabel='Column position (cm)',
                   ylabel='Iodine per 1 cm bin (% of common scan)', xlim=(28, 68))
            ax.legend(fontsize=7)
        fig.suptitle('Refinement comparisons; each curve normalised over the common published scan')
        fig.savefig(figure_dir / 'refinement_common_scan.png', dpi=150)
        fig.savefig(figure_dir / 'refinement_common_scan.pdf'); pdf.savefig(fig); plt.close(fig)
        fig, axes = plt.subplots(3, 3, figsize=(13, 11), layout='constrained')
        variations = [(['D5e-05', 'D0.0002', 'PbI2He'], ['D = 5e-5', 'D = 2e-4', 'PbI2He(T,p)'], 'Transport'),
                      (['alpha0.1', 'alpha0.01'], ['Accommodation = 0.1', 'Accommodation = 0.01'], 'Iodide accommodation'),
                      (['boatupstream', 'boatdownstream'], ['Source -1.5 cm', 'Source +1.5 cm'], 'Source position')]
        for column, base in enumerate(BASES[:3]):
            for row, (suffixes, labels, title) in enumerate(variations):
                ax = axes[row, column]
                if base not in folders: ax.set_visible(False); continue
                x, model, observed = scan_curve(folders[base])
                ax.fill_between(x, observed, color='.85', step='mid', label='Published scan')
                ax.plot(x, model, color='black', linewidth=1.5, label='Original baseline')
                for suffix, label in zip(suffixes, labels):
                    folder = out / 'profiles' / (base + '__' + suffix)
                    if folder.exists():
                        x, model, _ = scan_curve(folder); ax.plot(x, model, label=label)
                ax.set(xlim=(28, 68), xlabel='Position (cm)', ylabel='% of common scan',
                       title=(LABELS[base] + '\n' if row == 0 else '') + title)
                ax.legend(fontsize=7)
        fig.suptitle('Completed physical scenarios; separately normalised curves, no parameter fitting')
        fig.savefig(figure_dir / 'physical_sensitivity_profiles.png', dpi=150)
        fig.savefig(figure_dir / 'physical_sensitivity_profiles.pdf'); pdf.savefig(fig); plt.close(fig)
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>M10e interim postprocessing</title>
    <style>body{font:16px/1.6 system-ui;max-width:1150px;margin:32px auto;padding:0 20px}
    table{display:block;overflow-x:auto;border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ccc}
    img{max-width:100%}li{margin:10px 0}</style><body><h1>M10e interim postprocessing</h1>'''
    document += f'<p><strong>{analysis["completed_cases"]}/{analysis["planned_cases"]} follow-up cases analysed.</strong> '
    document += f'Generated {html.escape(summary["created_utc"])}. This is a frozen snapshot; remaining jobs continue separately.</p>'
    document += '<ul>' + ''.join('<li>' + html.escape(t) + '</li>' for t in findings) + '</ul>'
    document += '<h2>Numerical comparison</h2>' + tables_html(numerical_rows, [(k, k.replace('_', ' ')) for k in
                    ['base_case', 'dt_2_to_1_ms_I_profile_L1_pct', 'coarse_vs_medium_I_profile_L1_pct',
                     'fine_vs_medium_I_profile_L1_pct', 'fine_vs_medium_total_I_change_pct']])
    document += '<p>L1 above measures the summed absolute changes in deposited iodine across all native bins, '
    document += 'divided by the original total. It is not the error against the experiment or a confidence interval.</p>'
    for filename in ['numerical_comparison.png', 'refinement_common_scan.png', 'physical_sensitivity_profiles.png']:
        if (figure_dir / filename).exists():
            encoded = base64.b64encode((figure_dir / filename).read_bytes()).decode()
            document += f'<img src="data:image/png;base64,{encoded}" alt="{html.escape(filename.replace("_", " "))}">'
    document += '<h2>Agreement with the published scan</h2><p>The raw L1 comparison retains the published yields. '
    document += 'Overlap is a secondary shape metric: both curves are separately normalised over the common scan. '
    document += 'The measured scan includes KI, which the model excludes using the stated yield. '
    document += 'Q1 source splits are conditioned on the paper; agreement is not independent source-speciation validation. '
    document += 'Q2 source loading and the +29 kJ/mol BiI thermodynamic surrogate remain unqualified hypotheses.</p>'
    document += tables_html(agreement, [(k, k.replace('_', ' ')) for k in ['case', 'raw_L1_percentage_points',
                    'normalised_common_scan_overlap_pct', 'normalised_common_scan_Wasserstein_cm',
                    'BiI3_fraction_of_iodide_I_pct', 'published_BiI3_fraction_of_iodide_I_pct',
                    'Q1_source_split_conditioned']])
    document += '<h2>Numerical integrity</h2><p>All included cases passed element and re-speciation checks. '
    document += f'Maximum raw material residue relative to supplied material: {conservation["max_raw_material_residue_over_source"]:.3g}. '
    document += f'Maximum accounted closure: {conservation["max_accounted_closure_over_source"]:.3g}. '
    document += 'This supports numerical integrity; the experimental model remains provisional.</p>'
    document += '<h2>Files and reproducibility</h2><p>' + ' · '.join(
        f'<a href="{filename}">{label}</a>' for filename, label in [
            ('report.html', 'Detailed case report and conditional grid estimates'),
            ('postprocessing_figures.pdf', 'Figure PDF'),
            ('numerical_comparison.csv', 'Numerical comparison CSV'),
            ('experimental_comparison.csv', 'Experimental comparison CSV'),
            ('sensitivity_summary.csv', 'Sensitivity CSV'), ('peak_positions.csv', 'Peak positions CSV'),
            ('interpretation.json', 'Machine-readable interpretation')]) + '</p>'
    document += '<p>Reproduce by running report_followup.py for the selected study, then summarise_followup.py '
    document += 'with --root pointing to that study and --report pointing to the saved report. '
    document += 'Copied profiles, baseline references, metadata and checksums are included. '
    document += 'The original experimental temperature curves are retained without fitted shifts.</p></body></html>'
    (out / 'overview.html').write_text(document)
    write_json(out / 'artifact_sha256.json', {
        str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'artifact_sha256.json'})
    print(json.dumps({'completed_cases': analysis['completed_cases'], 'overview': str(out / 'overview.html'),
                      'findings': findings}, indent=2))


if __name__ == '__main__':
    main()
