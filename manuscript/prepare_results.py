#!/usr/bin/env python3
"""Build manuscript tables and figures from saved, completed CFD reports."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
BENCH = REPO / 'thermochemistry/benchmarks/liu2025'
sys.path.insert(0, str(BENCH))
from analyse_matrix import experimental_shape
from report import rows

FULL = REPO / 'run/034-liu2025/full_analysis'
FOLLOW = REPO / 'run/034-liu2025/followup_results'
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
SHORT = ['Pure PbI2', 'Steel Q2, 0.1%', 'Steel Q2, 1%', 'Steel Q2, 1%, shifted BiI']


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def curve(folder, species='element_I'):
    data = rows(folder / (species + '_native.dat'))
    return np.array([r['x'] for r in data]) * 100, np.array([
        r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])


def scan(folder, normalise=False):
    data = [r for r in read_csv(folder / 'iodine-comparison.csv') if r['figure_relative_yield_percent']]
    x = np.array([float(r['x_m']) for r in data]) * 100
    model = np.array([float(r['predicted_relative_yield_percent']) for r in data])
    observed = np.array([float(r['figure_relative_yield_percent']) for r in data])
    if normalise:
        model *= 100 / model.sum(); observed *= 100 / observed.sum()
    return x, model, observed


def peak(folder, species):
    data = rows(folder / (species + '_g_native.dat'))
    if not any(r['wall'] > 0 for r in data):
        return None, None
    r = max(data, key=lambda r: r['wall'])
    return 100 * r['x'], r['Twall'] - 273.15


def table(path, columns, header, records):
    text = '\\begin{tabular}{' + columns + '}\n\\toprule\n' + ' & '.join(header) + ' \\\\\n\\midrule\n'
    text += ''.join(' & '.join(map(str, r)) + ' \\\\\n' for r in records)
    path.write_text(text + '\\bottomrule\n\\end{tabular}\n')


def main():
    generated = HERE / 'generated'; generated.mkdir(exist_ok=True)
    figures = HERE / 'figures'; figures.mkdir(exist_ok=True)
    original = json.loads((FULL / 'analysis.json').read_text())
    follow = json.loads((FOLLOW / 'analysis.json').read_text())
    assert original['PASS'] == 15 and follow['completed_cases'] == 64
    assert follow['execution_status'] == 'COMPLETE' and not follow['missing_or_failed']
    metrics = follow['case_metrics']
    compares = [r for r in follow['comparisons'] if r['species'] == 'I']
    by_name = {r['case']: r for r in compares}
    original_metrics = {r['case']: r for r in original['case_metrics']}
    numeric = []
    for base, label in zip(BASES, SHORT):
        selected = [by_name[base + '__' + suffix] for suffix in ['dt1ms', 'grid256x32', 'grid576x72']]
        numeric.append({'case': base, 'label': label, 'time_L1_pct': 100 * selected[0]['profile_L1_relative_to_base'],
                        'coarse_L1_pct': 100 * selected[1]['profile_L1_relative_to_base'],
                        'fine_L1_pct': 100 * selected[2]['profile_L1_relative_to_base'],
                        'fine_total_change_pct': 100 * selected[2]['deposit_relative_change'],
                        'fine_centroid_shift_cm': selected[2]['centroid_shift_cm']})
    table(generated / 'numerical-table.tex', 'lrrrr',
          ['Case', '2 to 1 ms', 'Coarse/medium', 'Fine/medium', r'$\Delta N_I$'],
          [[r['label'].replace('%', r'\%'), f"{r['time_L1_pct']:.4f}", f"{r['coarse_L1_pct']:.2f}",
            f"{r['fine_L1_pct']:.2f}", f"{r['fine_total_change_pct']:.4f}"] for r in numeric])
    categories = ['numerical', 'transport', 'kinetics', 'source_geometry', 'source_loading',
                  'flow_reference', 'silica_alignment', 'shutdown']
    labels = ['Mesh and time step', 'Gas diffusion', 'Iodide accommodation', 'Source position',
              'Q2 metal loading', 'Flow reference', 'Silica temperature coordinates', 'Flow-stop relaxation']
    variation = {'numerical': 'Three grids; 2 and 1 ms', 'transport': r'$D/2$, $2D$, PbI2--He correlation',
                 'kinetics': r'$\alpha=0.1$, 0.01', 'source_geometry': r'$\Delta x_s=\pm1.5$ cm',
                 'source_loading': r'$\beta=0.002$, 0.003162, 0.005623',
                 'flow_reference': '273.15 versus 298.15 K', 'silica_alignment': r'$\Delta x_T=-4.10$, $-4.65$ cm',
                 'shutdown': '60--80 s; fixed temperature'}
    counts = {c: sum(r['category'] == c for r in metrics) for c in categories}
    table(generated / 'study-table.tex', 'lrl', ['Study', 'Cases', 'Variation'],
          [[label, counts[c], variation[c]] for c, label in zip(categories, labels)])
    sensitivity = []
    for category, label in zip(categories[1:], labels[1:]):
        values = [100 * r['profile_L1_relative_to_base'] for r in compares if r['category'] == category]
        sensitivity.append({'category': category, 'label': label, 'min_L1_pct': min(values), 'max_L1_pct': max(values)})
    table(generated / 'sensitivity-table.tex', 'lrr', ['Assumption', 'Minimum', 'Maximum'],
          [[r['label'], f"{r['min_L1_pct']:.2f}", f"{r['max_L1_pct']:.2f}"] for r in sensitivity])
    q1_labels = {'PbI2_SS': 'Pure PbI2, steel', 'LBE-I_SS_II': 'LBE, steel', 'LBE-I_SiO2_I': 'LBE, silica'}
    table(generated / 'peaks-table.tex', 'llrrrr',
          ['Case', 'Species', r'$x_{\rm exp}$', r'$x_{\rm CFD}$', r'$T_{\rm exp}$', r'$T_{\rm CFD}$'],
          [[q1_labels[r['experiment']], r'\ce{' + r['species'] + '}', f"{r['figure_x_cm']:.0f}",
            f"{r['model_x_cm']:.0f}", f"{r['published_peak_C']:.0f}" + r'$\pm$' + f"{r['uncertainty_K']:.0f}",
            f"{r['predicted_peak_C']:.1f}"] for r in original['q1_peaks']])
    alignment = []
    for mode, base in [('Q1', 'LBE-I_SiO2_I_Q1_baseline'), ('Q2', 'LBE-I_SiO2_I_Q2_baseline_sat0.001')]:
        for suffix, shift in [('', 0), ('curveShiftPbI2', -4.096367986080896), ('curveShiftBiI3', -4.6531462799495715)]:
            name = base + '__' + suffix if suffix else base
            folder = FOLLOW / 'profiles' / name if suffix else FULL / 'profiles' / name
            shape = experimental_shape(folder)
            xp, tp = peak(folder, 'PbI2'); xb, tb = peak(folder, 'BiI3')
            alignment.append({'case': name, 'mode': mode, 'shift_cm': shift, 'PbI2_x_cm': xp, 'BiI3_x_cm': xb,
                              'PbI2_C': tp, 'BiI3_C': tb, 'raw_L1_pp': shape['raw_L1_percentage_points'],
                              'overlap_pct': 100 * shape['normalised_scan_overlap']})
    table(generated / 'alignment-table.tex', 'lrrrrr',
          ['Source', r'$\Delta x_T$', r'$x_{\rm PbI_2}$', r'$x_{\rm BiI_3}$', 'Raw L1', 'Overlap'],
          [[r['mode'], f"{r['shift_cm']:.2f}", f"{r['PbI2_x_cm']:.0f}", f"{r['BiI3_x_cm']:.0f}",
            f"{r['raw_L1_pp']:.2f}", f"{r['overlap_pct']:.2f}"] for r in alignment])
    audit_keys = ['max_accounted_closure_over_source', 'max_raw_material_residue_over_source',
                  'max_net_clamped_over_source', 'max_solver_defect_over_source']
    integrity = {k: max(r[k] for a in follow['material_audits'].values() for r in a.values()
                        if isinstance(r, dict) and k in r) for k in audit_keys}
    macros = {'BaselineCases': 15, 'FollowupCases': 64, 'TimeMin': f"{min(r['time_L1_pct'] for r in numeric):.4f}",
              'TimeMax': f"{max(r['time_L1_pct'] for r in numeric):.3f}",
              'MeshMin': f"{min(r['fine_L1_pct'] for r in numeric):.2f}",
              'MeshMax': f"{max(r['fine_L1_pct'] for r in numeric):.2f}"}
    (generated / 'numbers.tex').write_text('% Generated from saved completed reports.\n' + ''.join(
        '\\newcommand{\\' + k + '}{' + str(v) + '}\n' for k, v in macros.items()))
    data = {'baseline_cases': original['PASS'], 'followup_cases': follow['completed_cases'], 'numerical': numeric,
            'sensitivity': sensitivity, 'alignment': alignment, 'integrity_60s_followup': integrity,
            'baseline_trace_mole_fraction_upper_bound': original['max_trace_mole_fraction_upper_bound'],
            'baseline_max_negative_I_over_supplied_I': original['max_negative_I_over_supplied_I'],
            'solver_sha256': follow['solver_sha256'], 'study_manifest_sha256': follow['study_manifest_sha256'],
            'scope': 'Completed numerical studies and conditional input sensitivities; no fitting or experimental qualification.'}
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    def save(fig, name):
        fig.savefig(figures / (name + '.pdf'), bbox_inches='tight',
                    metadata={'CreationDate': None, 'ModDate': None})
        fig.savefig(figures / (name + '.png'), dpi=150, bbox_inches='tight'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 2.0), layout='constrained')
    ax.add_patch(Rectangle((0, .45), 1.25, .16, facecolor='.93', edgecolor='black'))
    ax.add_patch(Rectangle((.09, .46), .015, .14, facecolor='#d89026', edgecolor='black'))
    ax.annotate('Helium', (.02, .53), xytext=(-.12, .53), arrowprops={'arrowstyle': '->'}, va='center')
    ax.annotate('Prescribed release region\nboat obstruction omitted', (.0975, .59), xytext=(.3, .87),
                arrowprops={'arrowstyle': '->'}, ha='center')
    ax.annotate('1.25 m computational column', (0, .17), xytext=(.625, .17), ha='center')
    ax.annotate('', (0, .32), (1.25, .32), arrowprops={'arrowstyle': '<->'})
    ax.text(.7, -.12, 'Steel diameter: 4.8 mm; silica diameter: 5.0 mm\n5 degree wedge; geometry sketch not to scale', ha='center')
    ax.set(xlim=(-.18, 1.32), ylim=(-.32, 1.1)); ax.axis('off'); save(fig, 'experiment_schematic')
    fig, ax = plt.subplots(figsize=(6.5, 3.1), layout='constrained')
    for name, label in [('PbI2_SS', 'Pure PbI2, steel'), ('LBE-I_SS_II', 'LBE, steel'), ('LBE-I_SiO2_I', 'LBE, silica')]:
        r = read_csv(BENCH / 'digitized' / (name + '_temperature.csv'))
        ax.plot([100 * float(v['x_m']) for v in r], [float(v['T_K']) - 273.15 for v in r], label=label)
    ax.axvspan(9, 10.5, alpha=.15, color='#d89026', label='Source region')
    ax.set(xlabel='Published column coordinate (cm)', ylabel='Temperature (degree C)', xlim=(0, 80)); ax.legend()
    save(fig, 'temperature_profiles')
    fig, axes = plt.subplots(3, 1, figsize=(6.5, 6.8), layout='constrained')
    for ax, base, label, limits in zip(axes, ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q1_baseline',
            'LBE-I_SiO2_I_Q1_baseline'], ['Pure PbI2, steel', 'LBE, steel', 'LBE, silica'], [(30, 43), (39, 65), (40, 68)]):
        folder = FULL / 'profiles' / base
        x, model, observed = scan(folder)
        ax.fill_between(x, observed, color='.8', step='mid', label='Published scan')
        ax.plot(x, model, color='black', label='Modelled total I')
        _, total = curve(folder)
        for species, atoms, color in [('PbI2_g', 2, '#007b9e'), ('BiI3_g', 3, '#974b9f')]:
            xp, amount = curve(folder, species)
            y = (100 - (0 if base.startswith('PbI2') else 6.6 if '_SS_' in base else 7)) * atoms * amount / total.sum()
            ax.plot(xp, y, color=color, linestyle='--', label=species.replace('_g', ''))
        ax.set(title=label, xlim=limits, xlabel='Position (cm)', ylabel='I yield per 1 cm bin (%)'); ax.legend(fontsize=8)
    save(fig, 'baseline_profiles')
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.4), layout='constrained')
    for ax, base, label in zip(axes.flat, BASES, SHORT):
        for suffix, tag in [('', '18,432 cells, 2 ms'), ('grid256x32', '8,192 cells, 2 ms'),
                            ('grid576x72', '41,472 cells, 2 ms'), ('dt1ms', '18,432 cells, 1 ms')]:
            folder = FOLLOW / 'profiles' / (base + '__' + suffix) if suffix else FULL / 'profiles' / base
            x, y, obs = scan(folder, normalise=True); ax.plot(x, y, label=tag)
        ax.fill_between(x, obs, color='.85', step='mid', label='Published scan')
        ax.set(title=label, xlim=(30, 43) if base.startswith('PbI2') else (39, 65),
               xlabel='Position (cm)', ylabel='I per bin (% of scan)')
    fig.legend(*axes.flat[0].get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.5, 1.09), ncol=2, fontsize=8)
    save(fig, 'refinement')
    chemistry = []
    for name, r in original_metrics.items():
        if r['source_mode'] == 'Q2':
            chemistry.append({'case': name, 'experiment': r['experiment'], 'variant': r['thermo_variant'],
                              'beta': r['metal_saturation'], 'BiI3_pct': r['BiI3_iodide_I_pct']})
    for r in metrics:
        if r['category'] == 'source_loading':
            chemistry.append({'case': r['case'], 'experiment': 'LBE-I_SiO2_I' if 'SiO2' in r['case'] else 'LBE-I_SS_II',
                              'variant': 'shift29k' if 'shift29k' in r['case'] else 'baseline',
                              'beta': float(r['case'].split('__sat')[1]), 'BiI3_pct': r['BiI3_iodide_I_pct']})
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.3), layout='constrained')
    for ax, experiment, label in zip(axes, ['LBE-I_SS_II', 'LBE-I_SiO2_I'], ['Steel', 'Silica']):
        for variant, tag in [('baseline', 'Baseline thermodynamics'), ('shift29k', 'BiI +29 kJ/mol surrogate')]:
            selected = sorted([r for r in chemistry if r['experiment'] == experiment and r['variant'] == variant], key=lambda r: r['beta'])
            ax.plot([100 * r['beta'] for r in selected], [r['BiI3_pct'] for r in selected], 'o-', label=tag)
        measured = next(r['observed_BiI3_iodide_I_pct'] for r in original_metrics.values() if r['experiment'] == experiment)
        ax.axhline(measured, color='.4', linestyle='--', label='Published iodide-I split')
        ax.set(xscale='log', xlabel='Metal saturation (%)', ylabel='BiI3 share of iodide I (%)',
               ylim=(0, 105), title=label)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.5, 1.16), ncol=2, fontsize=8)
    save(fig, 'source_thermodynamics')
    fig = plt.figure(figsize=(6.5, 5.9), layout='constrained')
    grid = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    r = read_csv(BENCH / 'digitized/LBE-I_SiO2_I_temperature.csv')
    for shift, label in [(0, 'Published curve'), (-4.096367986, 'Curve -4.10 cm'), (-4.65314628, 'Curve -4.65 cm')]:
        axes[0].plot([100 * float(v['x_m']) + shift for v in r], [float(v['T_K']) - 273.15 for v in r], label=label)
    axes[0].scatter([51, 59], [306, 127], marker='x', s=55, color='black', label='Printed peak labels')
    axes[0].set(title='Temperature-coordinate hypotheses', xlim=(45, 68), ylim=(50, 500),
                xlabel='Position (cm)', ylabel='Temperature (degree C)'); axes[0].legend(fontsize=8)
    for ax, mode, base in zip(axes[1:], ['Q1', 'Q2'], ['LBE-I_SiO2_I_Q1_baseline', 'LBE-I_SiO2_I_Q2_baseline_sat0.001']):
        for suffix, label in [('', 'Original inputs'), ('curveShiftPbI2', 'Curve -4.10 cm'), ('curveShiftBiI3', 'Curve -4.65 cm')]:
            folder = FOLLOW / 'profiles' / (base + '__' + suffix) if suffix else FULL / 'profiles' / base
            x, y, obs = scan(folder, normalise=True); ax.plot(x, y, label=label)
        ax.fill_between(x, obs, color='.85', step='mid', label='Published scan')
        ax.set(title=mode + ('; 0.1% saturation' if mode == 'Q2' else ''), xlim=(40, 68),
               xlabel='Position (cm)', ylabel='I per bin (% of scan)')
    fig.legend(*axes[1].get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.5, -.02), ncol=2, fontsize=8)
    save(fig, 'silica_alignment')
    shown = [r for r in sensitivity if r['category'] in ['transport', 'kinetics', 'source_geometry', 'flow_reference', 'shutdown']]
    fig, ax = plt.subplots(figsize=(6.5, 3.1), layout='constrained')
    ax.barh([r['label'] for r in shown], [r['max_L1_pct'] for r in shown], color='#2878a6')
    ax.invert_yaxis(); ax.set(xlabel='Maximum iodine profile L1 change from baseline (%)')
    save(fig, 'physical_sensitivity')
    data['chemistry'] = chemistry
    input_files = [FULL / 'analysis.json', FOLLOW / 'analysis.json',
                   BENCH / 'experiments.csv', BENCH / 'targets.csv', BENCH / 'digitized/provenance.json']
    input_files += list((BENCH / 'digitized').glob('*.csv'))
    input_files += [Path(__file__).resolve(), BENCH / 'analyse_matrix.py', BENCH / 'report.py',
                    BENCH / 'prepare.py', BENCH / 'carriers.py']
    input_files += list((FULL / 'profiles').rglob('*.dat')) + list((FULL / 'profiles').rglob('*.csv'))
    input_files += list((FOLLOW / 'profiles').rglob('*.dat')) + list((FOLLOW / 'profiles').rglob('*.csv'))
    data['input_sha256'] = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(input_files)}
    dump(generated / 'results.json', data)
    dump(generated / 'artifact_sha256.json', {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in [generated, figures] for p in sorted(directory.iterdir())
        if p.is_file() and p.name != 'artifact_sha256.json'})
    print('Generated 7 figures, 5 tables and numerical provenance for 15 baseline + 64 follow-up cases.')


if __name__ == '__main__':
    main()
