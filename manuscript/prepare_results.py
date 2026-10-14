#!/usr/bin/env python3
"""Build manuscript tables and figures from saved, completed CFD reports."""
import csv
import hashlib
import json
import math
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

FULL = REPO / 'run/024-liu2025/full_analysis'
FOLLOW = REPO / 'run/024-liu2025/followup_results'
# Production baselines: 576 x 72 (41,472 cells), Gauss limitedLinear 1 for the species advection, otherwise the
# inputs of the 18,432-cell cases.  The directory holds analysis.json (q1_peaks, case_metrics, ...) and
# profiles/<case>/ for the 15 baseline and the 4 silica alignment cases; FULL keeps the 18,432-cell upwind
# originals, which remain the reference of the 64 follow-up sensitivities.
PROD = REPO / 'run/024-liu2025/production_576x72_limitedLinear'
# Scheme and mesh study (plan sections 38-39): two cases, three schemes, up to four levels.
STUDY = REPO / 'run/024-liu2025/scheme_mesh_study'
# Full-run corrector check: eight versus three species correctors over 60 s.
ITER = REPO / 'run/024-liu2025/iteration_check'
# Head-to-head on Lobresco et al. (2026): digitised published data and our ten runs.
LOB_BENCH = REPO / 'thermochemistry/benchmarks/lobresco2026'
LOB_RUN = REPO / 'run/023-lobresco2026'
LOB_CASES = ['T106_6', 'T106_60', 'T106_300', 'T5_60', 'T540_60']
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.001',
         'LBE-I_SS_II_Q2_baseline_sat0.01', 'LBE-I_SS_II_Q2_shift29k_sat0.01']
SHORT = ['Pure PbI2', 'Steel Q2, 0.1%', 'Steel Q2, 1%', 'Steel Q2, 1%, shifted BiI']
STUDY_BASES = {'pure': 'PbI2_SS_Q1_baseline', 'steel': 'LBE-I_SS_II_Q2_shift29k_sat0.01'}
CELLS = {1: 8192, 2: 18432, 3: 41472, 4: 93312}
RATIO = 1.5
# Final figure width in inches.  The PoF reprint text width is 396 pt = 5.48 in; every figure is drawn at
# its printed size and included without scaling, so the font sizes below are the printed sizes.
WIDTH = 5.4
SPECIES_TEX = {'PbI2': r'PbI$_2$', 'BiI3': r'BiI$_3$'}


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


def scan_from_profile(folder, experiment):
    """Common-scan-normalised model and experimental iodine yields from a bare profile folder."""
    x, amount = curve(folder)
    gamma = read_csv(BENCH / 'digitized' / (experiment + '_gamma.csv'))
    xs = np.array([100 * float(r['x_m']) for r in gamma])
    index = [int(np.argmin(abs(x - v))) for v in xs]
    assert np.allclose(x[index], xs)
    model = amount[index]; observed = np.array([float(r['relative_yield_percent_per_cm']) for r in gamma])
    return xs, 100 * model / model.sum(), 100 * observed / observed.sum()


def smoothed_peak(folder, species):
    """3 cm (three-bin) smoothed maximum and the wall temperature at the deposit centroid."""
    data = rows(folder / (species + '_g_native.dat'))
    amount = np.array([r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])
    x = np.array([r['x'] for r in data]); wall_t = np.array([r['Twall'] for r in data])
    j = int(np.argmax(np.convolve(amount, np.ones(3) / 3, 'same')))
    centroid = float(x @ amount / amount.sum())
    return 100 * x[j], wall_t[j] - 273.15, float(np.interp(centroid, x, wall_t)) - 273.15


def peak(folder, species):
    data = rows(folder / (species + '_g_native.dat'))
    if not any(r['wall'] > 0 for r in data):
        return None, None
    r = max(data, key=lambda r: r['wall'])
    return 100 * r['x'], r['Twall'] - 273.15


def table(path, columns, header, records, midrules=()):
    text = '\\begin{tabular}{' + columns + '}\n\\toprule\n' + ' & '.join(header) + ' \\\\\n\\midrule\n'
    for i, r in enumerate(records):
        if i in midrules:
            text += '\\midrule\n'
        text += ' & '.join(map(str, r)) + ' \\\\\n'
    path.write_text(text + '\\bottomrule\n\\end{tabular}\n')


def study_folder(key, level, scheme):
    """Final 60 s profile folder of one member of the scheme/mesh study."""
    base = STUDY_BASES[key]
    if scheme == 'up' and level < 4:
        return {1: FOLLOW / 'profiles' / (base + '__grid256x32'), 2: FULL / 'profiles' / base,
                3: FOLLOW / 'profiles' / (base + '__grid576x72')}[level]
    return STUDY / 'profiles' / f'{key}_L{level}_{scheme}'


def relative_l1(a, b):
    return 100 * float(np.abs(a - b).sum() / np.abs(b).sum())


def onset(x, y, wall_t, fraction=0.01):
    """Lobresco et al. (2026) rule: first 1 cm bin from upstream above 1 % of the peak, linear interpolation."""
    threshold = fraction * y.max()
    k = int(np.flatnonzero(y >= threshold)[0])
    assert k > 0
    xo = x[k - 1] + (threshold - y[k - 1]) / (y[k] - y[k - 1]) * (x[k] - x[k - 1])
    return float(xo), float(np.interp(xo, *wall_t)), float(x[int(np.argmax(y))])


def lobresco_case(level, name):
    folder = LOB_RUN / 'cases' / level / name
    info = json.loads((folder / 'case.json').read_text())
    wall = read_csv(LOB_RUN / 'inputs' / f"wallT_{info['table1']['flow_mL_min']}mLmin.csv")
    wall_t = (np.array([float(r['x_m']) for r in wall]), np.array([float(r['T_K']) for r in wall]))
    final = folder / 'final_profiles'   # postProcessing/axialProfiles/<end time>/ of the run, flattened
    assert abs(float((final / 'time').read_text()) - info['end_time_s']) < 1e-6
    data = rows(final / 'PbI2_g_gamma.dat')
    scale = 360 / 5   # wedge to full column
    x = np.array([r['x'] for r in data])
    deposit = np.array([r['wall'] for r in data]) * scale           # mol/m, full column
    with_gas = np.array([r['observable'] for r in data]) * scale    # deposit plus gas left at the end time
    gas = np.array([r['gas'] for r in data]) * scale
    boat_end = info['boat']['x_m'][1]
    main = x > boat_end   # the main (downstream) deposit; T5_60 also has a back-diffused upstream deposit
    exp = read_csv(LOB_BENCH / 'digitized' / f'{name}_experiment.csv')
    xe = np.array([float(r['x_bin_m']) for r in exp]); ye = np.array([float(r['n_mol_per_m']) for r in exp])
    ye_full = np.array([ye[np.argmin(abs(xe - v))] if np.min(abs(xe - v)) < 1e-6 else 0.0 for v in x])
    n0 = info['table1']['moles']
    x_main, t_main, peak_main = onset(x[main], deposit[main], wall_t)
    x_inlet, t_inlet, _ = onset(x, deposit, wall_t)
    _, t_gas, _ = onset(x[main], with_gas[main], wall_t)
    x_exp, t_exp, peak_exp = onset(xe, ye, wall_t)
    return {'case': name, 'level': level, 'cells': info['carrier']['cells'], 'flow_mL_min': info['table1']['flow_mL_min'],
            'duration_s': info['table1']['duration_s'], 'n0_mol': n0,
            'T_onset_deposit_K': t_main, 'x_onset_deposit_m': x_main, 'T_onset_from_inlet_K': t_inlet,
            'x_onset_from_inlet_m': x_inlet, 'T_onset_deposit_plus_gas_K': t_gas,
            'x_peak_m': peak_main, 'x_peak_experiment_m': peak_exp, 'T_onset_digitised_experiment_K': t_exp,
            'deposit_over_n0': float(deposit.sum() * 0.01 / n0), 'gas_fraction': float(gas.sum() / (gas.sum() + deposit.sum())),
            'L1_deposit_minus_experiment_over_n0': float(np.abs(deposit - ye_full).sum() * 0.01 / n0),
            'x_m': x, 'deposit_mol_per_m': deposit, 'wall_T': wall_t, 'experiment': (xe, ye)}


def main():
    generated = HERE / 'generated'; generated.mkdir(exist_ok=True)
    figures = HERE / 'figures'; figures.mkdir(exist_ok=True)
    legacy = json.loads((FULL / 'analysis.json').read_text())   # 18,432 cells, upwind
    original = json.loads((PROD / 'analysis.json').read_text())  # production baselines
    assert original['PASS'] == 15 and original['cells'] == 41472
    assert original['scheme'] == 'Gauss limitedLinear 1' and len(original['silica_alignment_cases']) == 4
    follow = json.loads((FOLLOW / 'analysis.json').read_text())
    assert follow['completed_cases'] == 64
    assert follow['execution_status'] == 'COMPLETE' and not follow['missing_or_failed']
    assert original['solver_sha256'] == follow['solver_sha256']
    metrics = follow['case_metrics']
    compares = [r for r in follow['comparisons'] if r['species'] == 'I']
    by_name = {r['case']: r for r in compares}
    original_metrics = {r['case']: r for r in original['case_metrics']}
    ledger = json.loads((PROD / 'ledger.json').read_text())
    assert ledger['n_pass'] == 19 and set(ledger['status'].values()) == {'PASS'}
    assert ledger['solver_sha256'] == original['solver_sha256']
    production_residue = max(r['max_raw_residue_new'] for r in ledger['ledger'])

    # ---- numerical table: time step (follow-up) and the scheme/mesh study --------------------------------
    numeric = []
    for base, label in zip(BASES, SHORT):
        selected = [by_name[base + '__' + suffix] for suffix in ['dt1ms', 'grid256x32', 'grid576x72']]
        numeric.append({'case': base, 'label': label, 'time_L1_pct': 100 * selected[0]['profile_L1_relative_to_base'],
                        'coarse_L1_pct': 100 * selected[1]['profile_L1_relative_to_base'],
                        'fine_L1_pct': 100 * selected[2]['profile_L1_relative_to_base'],
                        'fine_total_change_pct': 100 * selected[2]['deposit_relative_change'],
                        'fine_centroid_shift_cm': selected[2]['centroid_shift_cm']})
    time_step = {r['case']: r['time_L1_pct'] for r in numeric}
    study = json.loads((STUDY / 'results.json').read_text())
    families = {}
    for key in STUDY_BASES:
        for scheme in ['up', 'll', 'vl']:
            levels = [L for L in [1, 2, 3, 4] if study_folder(key, L, scheme).exists()]
            amounts = {L: curve(study_folder(key, L, scheme))[1] for L in levels}
            steps = [relative_l1(amounts[a], amounts[b]) for a, b in zip(levels[:-1], levels[1:])]
            reported = study[f'{key}.{scheme}']['successive']
            assert np.allclose(steps, [reported[f'L{a}->L{b}']['I_profile_L1_pct']
                                       for a, b in zip(levels[:-1], levels[1:])], rtol=1e-9), (key, scheme)
            families[f'{key}.{scheme}'] = {'levels': levels, 'successive_L1_pct': steps}
    for base in BASES[1:3]:   # the other two steel cases: upwind follow-up levels only
        amounts = [curve(FOLLOW / 'profiles' / (base + '__grid256x32'))[1], curve(FULL / 'profiles' / base)[1],
                   curve(FOLLOW / 'profiles' / (base + '__grid576x72'))[1]]
        families[base + '.up'] = {'levels': [1, 2, 3], 'successive_L1_pct': [
            relative_l1(amounts[0], amounts[1]), relative_l1(amounts[1], amounts[2])]}
    scheme_names = {'up': 'upwind', 'll': 'limitedLinear', 'vl': 'van Leer'}
    records = []
    for key, family, label, case in [('pure', 'pure', 'Pure PbI$_2$', BASES[0]),
                                     ('steel', 'steel', r'Steel, shifted BiI', BASES[3]),
                                     (BASES[1], None, r'Steel Q2, 0.1\%', BASES[1]),
                                     (BASES[2], None, r'Steel Q2, 1\%', BASES[2])]:
        for i, scheme in enumerate(['up', 'll', 'vl'] if family else ['up']):
            steps = families[f'{key}.{scheme}']['successive_L1_pct']
            cells = [f'{v:.2f}' for v in steps] + ['--'] * (3 - len(steps))
            records.append([label if i == 0 else '', scheme_names[scheme]] + cells +
                           [f'{time_step[case]:.4f}' if scheme == 'up' else '--'])
    table(generated / 'numerical-table.tex', 'llrrrr',
          ['Case', 'Scheme', r'8,192$\to$18,432', r'18,432$\to$41,472', r'41,472$\to$93,312', r'2$\to$1 ms'],
          records, midrules=(3, 6, 7))
    cross = study['cross_scheme_I_profile_L1_pct']
    vl_ll = {}
    for key, levels in [('pure', [1, 2, 3]), ('steel', [1, 2, 3])]:
        for L in levels:
            vl_ll[f'{key}.L{L}'] = relative_l1(curve(study_folder(key, L, 'vl'))[1], curve(study_folder(key, L, 'll'))[1])
    steel_ll = study['steel.ll']
    d23, d34 = families['steel.ll']['successive_L1_pct'][1:]
    order_last = math.log(d23 / d34) / math.log(RATIO)
    error_l3_apparent = d34 * RATIO ** order_last / (RATIO ** order_last - 1)
    error_l3_second = d34 * RATIO ** 2 / (RATIO ** 2 - 1)
    centroids = [steel_ll['functionals']['I.centroid_m']['values'][str(L)] for L in [1, 2, 3, 4]]
    centroid_t = [steel_ll['per_level'][str(L)]['peaks']['PbI2']['T_at_centroid_C'] for L in [1, 2, 3, 4]]
    totals = [steel_ll['per_level'][str(L)]['total_I_wedge_mol'] for L in [1, 2, 3, 4]]
    scheme_study = {'families': families, 'vanLeer_vs_limitedLinear_L1_pct': vl_ll,
                    'upwind_vs_limitedLinear_L1_pct': {k: v for k, v in cross.items() if 'll_vs_up' in k},
                    'steel_ll_apparent_orders': study['steel.ll']['apparent_order_from_differences'],
                    'steel_ll_L3_error_estimate_apparent_order_pct': error_l3_apparent,
                    'steel_ll_L3_error_estimate_second_order_pct': error_l3_second,
                    'steel_ll_centroid_range_mm': 1000 * (max(centroids) - min(centroids)),
                    'steel_ll_centroid_wall_T_range_K': max(centroid_t) - min(centroid_t),
                    'steel_ll_total_range_pct': 100 * (max(totals) - min(totals)) / totals[-1],
                    'steel_ll_L3_L4_3cm_L1_pct': steel_ll['successive']['L3->L4']['I_profile_L1_3cm_pct'],
                    'steel_ll_L3_L4_total_change_pct': steel_ll['successive']['L3->L4']['total_I_change_pct']}

    # ---- iteration check ------------------------------------------------------------------------------
    iteration = json.loads((ITER / 'summary.json').read_text())
    iteration_l1 = {r['base_case']: r['eight_minus_three']['element_I']['profile_L1_pct'] for r in iteration['cases']}
    assert set(iteration_l1) == set(BASES[:3])

    # ---- study, sensitivity, peaks and alignment tables -------------------------------------------------
    categories = ['numerical', 'transport', 'kinetics', 'source_geometry', 'source_loading',
                  'flow_reference', 'silica_alignment', 'shutdown']
    labels = ['Mesh and time step', 'Gas diffusion', 'Iodide accommodation', 'Source position',
              'Q2 metal loading', 'Flow reference', 'Silica temperature coordinates', 'Flow-stop relaxation']
    variation = {'numerical': 'Three grids; 2 and 1 ms', 'transport': r'$D/2$, $2D$, PbI$_2$--He correlation',
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
        if category == 'silica_alignment':   # production pair: shifted versus unshifted, both 41,472 cells, limitedLinear
            values = []
            for base in ['LBE-I_SiO2_I_Q1_baseline', 'LBE-I_SiO2_I_Q2_baseline_sat0.001']:
                _, reference = curve(PROD / 'profiles' / base)
                for suffix in ['curveShiftPbI2', 'curveShiftBiI3']:
                    _, shifted = curve(PROD / 'profiles' / (base + '__' + suffix))
                    values.append(100 * np.abs(shifted - reference).sum() / reference.sum())
        sensitivity.append({'category': category, 'label': label, 'min_L1_pct': min(values), 'max_L1_pct': max(values)})
    table(generated / 'sensitivity-table.tex', 'lrr', ['Assumption', 'Minimum', 'Maximum'],
          [[r['label'], f"{r['min_L1_pct']:.2f}", f"{r['max_L1_pct']:.2f}"] for r in sensitivity])
    q1_labels = {'PbI2_SS': r'Pure PbI$_2$, steel', 'LBE-I_SS_II': 'LBE, steel', 'LBE-I_SiO2_I': 'LBE, silica'}
    table(generated / 'peaks-table.tex', 'llrrrrrrr',
          ['Case', 'Species', r'$x_{\rm exp}$', r'$x_{\rm bin}$', r'$x_{3}$', r'$T_{\rm exp}$',
           r'$T_{\rm bin}$', r'$T_{3}$', r'$T_{\rm c}$'],
          [[q1_labels[r['experiment']], r'\ce{' + r['species'] + '}', f"{r['figure_x_cm']:.0f}",
            f"{r['model_x_cm']:.0f}", f"{r['smoothed_peak_x_cm']:.0f}",
            f"{r['published_peak_C']:.0f}" + r'$\pm$' + f"{r['uncertainty_K']:.0f}",
            f"{r['predicted_peak_C']:.1f}", f"{r['smoothed_peak_C']:.1f}", f"{r['centroid_T_C']:.1f}"]
           for r in original['q1_peaks']])
    alignment = []
    for mode, base in [('Q1', 'LBE-I_SiO2_I_Q1_baseline'), ('Q2', 'LBE-I_SiO2_I_Q2_baseline_sat0.001')]:
        for suffix, shift in [('', 0), ('curveShiftPbI2', -4.096367986080896), ('curveShiftBiI3', -4.6531462799495715)]:
            name = base + '__' + suffix if suffix else base
            folder = PROD / 'profiles' / name
            shape = experimental_shape(folder)
            xp, tp = peak(folder, 'PbI2'); xb, tb = peak(folder, 'BiI3')
            xps, tps, _ = smoothed_peak(folder, 'PbI2'); xbs, tbs, _ = smoothed_peak(folder, 'BiI3')
            alignment.append({'case': name, 'mode': mode, 'shift_cm': shift, 'PbI2_x_cm': xp, 'BiI3_x_cm': xb,
                              'PbI2_C': tp, 'BiI3_C': tb, 'PbI2_smoothed_x_cm': xps, 'BiI3_smoothed_x_cm': xbs,
                              'PbI2_smoothed_C': tps, 'BiI3_smoothed_C': tbs,
                              'raw_L1_pp': shape['raw_L1_percentage_points'],
                              'overlap_pct': 100 * shape['normalised_scan_overlap']})
    table(generated / 'alignment-table.tex', 'lrrrrr',
          ['Source', r'$\Delta x_T$', r'$x_{\rm PbI_2}$', r'$x_{\rm BiI_3}$', 'Raw L1', 'Overlap'],
          [[r['mode'], f"{r['shift_cm']:.2f}", f"{r['PbI2_x_cm']:.0f}", f"{r['BiI3_x_cm']:.0f}",
            f"{r['raw_L1_pp']:.2f}", f"{r['overlap_pct']:.2f}"] for r in alignment])

    # ---- Lobresco et al. (2026) comparison ----------------------------------------------------------------
    provenance = json.loads((LOB_BENCH / 'digitized/provenance.json').read_text())
    table2 = provenance['table2_verbatim']
    lob = {level: [lobresco_case(level, name) for name in LOB_CASES] for level in ['prod', 'coarse']}
    lob_rows = []
    for r, c in zip(lob['prod'], lob['coarse']):
        t = table2[r['case']]
        lob_rows.append({'case': r['case'], 'experiment_K': t['experiment_K'],
                         'ours_K': r['T_onset_deposit_K'], 'ours_minus_exp_K': r['T_onset_deposit_K'] - t['experiment_K'],
                         'ours_coarse_minus_exp_K': c['T_onset_deposit_K'] - t['experiment_K'],
                         'ours_with_gas_minus_exp_K': r['T_onset_deposit_plus_gas_K'] - t['experiment_K'],
                         'ours_from_inlet_K': r['T_onset_from_inlet_K'],
                         'temperature_based_minus_exp_K': t['temperature_based_K'] - t['experiment_K'],
                         'vapour_pressure_based_minus_exp_K': t['vapour_pressure_based_K'] - t['experiment_K'],
                         'peak_shift_cm': 100 * (r['x_peak_m'] - r['x_peak_experiment_m']),
                         'L1_over_n0': r['L1_deposit_minus_experiment_over_n0'],
                         'deposit_over_n0': r['deposit_over_n0'], 'gas_fraction': r['gas_fraction'],
                         'cells_prod': r['cells'], 'cells_coarse': c['cells']})
    rms = lambda key: math.sqrt(sum(r[key] ** 2 for r in lob_rows) / len(lob_rows))
    lob_summary = {'rms_ours_K': rms('ours_minus_exp_K'), 'rms_ours_coarse_K': rms('ours_coarse_minus_exp_K'),
                   'rms_ours_with_gas_K': rms('ours_with_gas_minus_exp_K'),
                   'rms_temperature_based_K': rms('temperature_based_minus_exp_K'),
                   'rms_vapour_pressure_based_K': rms('vapour_pressure_based_minus_exp_K'),
                   'max_coarse_vs_prod_K': max(abs(r['ours_coarse_minus_exp_K'] - r['ours_minus_exp_K']) for r in lob_rows),
                   'max_digitised_vs_table2_K': max(abs(c['difference_K']) for c in provenance['tdep_checks'])}
    flows = {r['case']: (r['flow_mL_min'], r['duration_s']) for r in lob['prod']}
    sign = lambda v: f'{v:+.0f}'.replace('-', '$-$')
    table(generated / 'lobresco-table.tex', 'lrrrrrrrr',
          ['Case', r'$\dot V$', '$t$', r'$T_{\rm exp}$', 'This work', 'Temp.-based', 'Vap.-press.',
           r'$\Delta x_{\rm pk}$', r'$L_1/n_0$'],
          [[r['case'].replace('_', r'\_'), flows[r['case']][0], flows[r['case']][1], r['experiment_K'],
            sign(r['ours_minus_exp_K']), sign(r['temperature_based_minus_exp_K']),
            sign(r['vapour_pressure_based_minus_exp_K']), f"{r['peak_shift_cm']:.0f}".replace('-', '$-$'),
            f"{r['L1_over_n0']:.2f}"] for r in lob_rows] +
          [['RMS', '', '', '', f"{lob_summary['rms_ours_K']:.0f}", f"{lob_summary['rms_temperature_based_K']:.0f}",
            f"{lob_summary['rms_vapour_pressure_based_K']:.0f}", '', '']], midrules=(5,))

    # ---- macros --------------------------------------------------------------------------------------------
    audit_keys = ['max_accounted_closure_over_source', 'max_raw_material_residue_over_source',
                  'max_net_clamped_over_source', 'max_solver_defect_over_source']
    integrity = {k: max(r[k] for a in follow['material_audits'].values() for r in a.values()
                        if isinstance(r, dict) and k in r) for k in audit_keys}
    sci = lambda v: (lambda m, e: f'{m}\\times10^{{{int(e)}}}')(*f'{v:.2e}'.split('e'))
    ll_steps = families['steel.ll']['successive_L1_pct']; up_steps = families['steel.up']['successive_L1_pct']
    pure_ll = families['pure.ll']['successive_L1_pct']; pure_up = families['pure.up']['successive_L1_pct']
    up_ll = [cross[f'steel.L{L}.ll_vs_up'] for L in [1, 2, 3, 4]]
    t_lob = {r['case']: r for r in lob_rows}
    lob_prod = {r['case']: r for r in lob['prod']}
    macros = {'BaselineCases': 15, 'FollowupCases': 64, 'ProductionCases': 19,
              'TimeMin': f"{min(r['time_L1_pct'] for r in numeric):.4f}",
              'TimeMax': f"{max(r['time_L1_pct'] for r in numeric):.3f}",
              'MeshMin': f"{min(r['fine_L1_pct'] for r in numeric):.2f}",
              'MeshMax': f"{max(r['fine_L1_pct'] for r in numeric):.2f}",
              'ProductionResidue': sci(production_residue),
              'SteelLLOne': f'{ll_steps[0]:.2f}', 'SteelLLTwo': f'{ll_steps[1]:.2f}', 'SteelLLThree': f'{ll_steps[2]:.2f}',
              'SteelUpOne': f'{up_steps[0]:.2f}', 'SteelUpTwo': f'{up_steps[1]:.2f}', 'SteelUpThree': f'{up_steps[2]:.2f}',
              'PureLLOne': f'{pure_ll[0]:.2f}', 'PureLLTwo': f'{pure_ll[1]:.2f}',
              'PureUpOne': f'{pure_up[0]:.2f}', 'PureUpTwo': f'{pure_up[1]:.2f}',
              'SteelOrderOne': f"{scheme_study['steel_ll_apparent_orders'][0]:.1f}",
              'SteelOrderTwo': f"{scheme_study['steel_ll_apparent_orders'][1]:.1f}",
              'SteelErrApparent': f'{error_l3_apparent:.1f}', 'SteelErrSecond': f'{error_l3_second:.1f}',
              'SteelSmoothThree': f"{scheme_study['steel_ll_L3_L4_3cm_L1_pct']:.2f}",
              'SteelTotalThree': f"{abs(scheme_study['steel_ll_L3_L4_total_change_pct']):.3f}",
              'SteelCentroidRange': f"{scheme_study['steel_ll_centroid_range_mm']:.1f}",
              'SteelCentroidTRange': f"{scheme_study['steel_ll_centroid_wall_T_range_K']:.1f}",
              'UpLLOne': f'{up_ll[0]:.1f}', 'UpLLTwo': f'{up_ll[1]:.1f}', 'UpLLThree': f'{up_ll[2]:.1f}',
              'UpLLFour': f'{up_ll[3]:.1f}',
              'PureUpLLMin': f"{min(cross[f'pure.L{L}.ll_vs_up'] for L in [1, 2, 3]):.1f}",
              'PureUpLLMax': f"{max(cross[f'pure.L{L}.ll_vs_up'] for L in [1, 2, 3]):.1f}",
              'VLSteelMax': f"{max(vl_ll[f'steel.L{L}'] for L in [1, 2, 3]):.1f}",
              'VLPureMax': f"{max(vl_ll[f'pure.L{L}'] for L in [1, 2, 3]):.1f}",
              'VLSteelThree': f"{vl_ll['steel.L3']:.2f}", 'VLPureThree': f"{vl_ll['pure.L3']:.2f}",
              'IterPure': sci(iteration_l1[BASES[0]]), 'IterSteelLow': sci(iteration_l1[BASES[1]]),
              'IterSteelHigh': f'{iteration_l1[BASES[2]]:.3f}',
              'LobRMS': f"{lob_summary['rms_ours_K']:.0f}", 'LobRMSTemp': f"{lob_summary['rms_temperature_based_K']:.0f}",
              'LobRMSVap': f"{lob_summary['rms_vapour_pressure_based_K']:.0f}",
              'LobRMSGas': f"{lob_summary['rms_ours_with_gas_K']:.0f}",
              'LobCoarse': f"{lob_summary['max_coarse_vs_prod_K']:.1f}",
              'LobDigitised': f"{lob_summary['max_digitised_vs_table2_K']:.1f}",
              'LobFiveInlet': f"{t_lob['T5_60']['ours_from_inlet_K']:.0f}",
              'LobFiveGas': sign(t_lob['T5_60']['ours_with_gas_minus_exp_K']),
              'LobFiveGasFraction': f"{100 * t_lob['T5_60']['gas_fraction']:.1f}",
              'LobLOneMin': f"{min(r['L1_over_n0'] for r in lob_rows):.1f}",
              'LobLOneMax': f"{max(r['L1_over_n0'] for r in lob_rows):.1f}",
              'LobShiftMin': f"{min(-r['peak_shift_cm'] for r in lob_rows):.0f}",
              'LobShiftMax': f"{max(-r['peak_shift_cm'] for r in lob_rows):.0f}",
              'LobCellsProd': f"{lob_prod['T106_60']['cells']:,}".replace(',', '{,}'),
              'LobCellsCoarse': f"{lob['coarse'][0]['cells']:,}".replace(',', '{,}')}
    for r, tag in zip(lob_rows, 'ABCDE'):   # model minus experiment onset, in LOB_CASES order
        macros['LobDiff' + tag] = sign(r['ours_minus_exp_K'])
    (generated / 'numbers.tex').write_text('% Generated from saved completed reports.\n' + ''.join(
        '\\newcommand{\\' + k + '}{' + str(v) + '}\n' for k, v in macros.items()))
    data = {'baseline_cases': original['PASS'], 'baseline_cells': original['cells'],
            'baseline_scheme': original['scheme'], 'q1_peaks': original['q1_peaks'],
            'followup_cases': follow['completed_cases'], 'numerical': numeric, 'scheme_mesh_study': scheme_study,
            'iteration_check_iodine_L1_pct': iteration_l1, 'production_max_raw_residue': production_residue,
            'sensitivity': sensitivity, 'alignment': alignment, 'integrity_60s_followup': integrity,
            'lobresco2026': {'rows': lob_rows, 'summary': lob_summary},
            'baseline_trace_mole_fraction_upper_bound': original['max_trace_mole_fraction_upper_bound'],
            'legacy_18432_upwind_q1_peaks': legacy['q1_peaks'],
            'baseline_max_negative_I_over_supplied_I': original['max_negative_I_over_supplied_I'],
            'solver_sha256': follow['solver_sha256'], 'study_manifest_sha256': follow['study_manifest_sha256'],
            'macros': macros,
            'scope': 'Completed numerical studies and conditional input sensitivities; no fitting or experimental qualification.'}

    # ---- figures: drawn at printed size, no text below 8 pt ------------------------------------------------
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'legend.fontsize': 8, 'legend.title_fontsize': 8,
                         'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.2,
                         'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'ps.fonttype': 42})
    def save(fig, name):
        fig.savefig(figures / (name + '.pdf'), metadata={'CreationDate': None, 'ModDate': None})
        fig.savefig(figures / (name + '.png'), dpi=200); plt.close(fig)
    fig, ax = plt.subplots(figsize=(WIDTH, 1.9), layout='constrained')
    ax.add_patch(Rectangle((0, .45), 1.25, .16, facecolor='.93', edgecolor='black'))
    ax.add_patch(Rectangle((.09, .46), .015, .14, facecolor='#d89026', edgecolor='black'))
    ax.annotate('Helium', (.02, .53), xytext=(-.14, .53), arrowprops={'arrowstyle': '->'}, va='center')
    ax.annotate('Prescribed release region;\nboat obstruction omitted', (.0975, .59), xytext=(.36, .9),
                arrowprops={'arrowstyle': '->'}, ha='center')
    ax.annotate('1.25 m computational column', (0, .17), xytext=(.625, .17), ha='center')
    ax.annotate('', (0, .32), (1.25, .32), arrowprops={'arrowstyle': '<->'})
    ax.text(.625, -.17, 'Inner diameter: 4.8 mm (steel), 5.0 mm (silica); 5° wedge; not to scale', ha='center')
    ax.set(xlim=(-.18, 1.3), ylim=(-.3, 1.1)); ax.axis('off'); save(fig, 'experiment_schematic')
    fig, ax = plt.subplots(figsize=(4.2, 2.7), layout='constrained')
    for name, label in [('PbI2_SS', r'Pure PbI$_2$, steel'), ('LBE-I_SS_II', 'LBE, steel'), ('LBE-I_SiO2_I', 'LBE, silica')]:
        r = read_csv(BENCH / 'digitized' / (name + '_temperature.csv'))
        ax.plot([100 * float(v['x_m']) for v in r], [float(v['T_K']) - 273.15 for v in r], label=label)
    ax.axvspan(9, 10.5, alpha=.15, color='#d89026', label='Source region')
    ax.set(xlabel='Published column coordinate (cm)', ylabel='Temperature (°C)', xlim=(0, 80)); ax.legend()
    save(fig, 'temperature_profiles')
    fig, axes = plt.subplots(3, 1, figsize=(WIDTH, 6.6), layout='constrained')
    for ax, base, label, limits in zip(axes, ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q1_baseline',
            'LBE-I_SiO2_I_Q1_baseline'], [r'Pure PbI$_2$, steel', 'LBE, steel', 'LBE, silica'], [(30, 43), (39, 65), (40, 68)]):
        folder = PROD / 'profiles' / base
        x, model, observed = scan(folder)
        ax.fill_between(x, observed, color='.8', step='mid', label='Published scan')
        ax.plot(x, model, color='black', label='Modelled total I')
        _, total = curve(folder)
        for species, atoms, color in [('PbI2_g', 2, '#007b9e'), ('BiI3_g', 3, '#974b9f')]:
            xp, amount = curve(folder, species)
            y = (100 - (0 if base.startswith('PbI2') else 6.6 if '_SS_' in base else 7)) * atoms * amount / total.sum()
            ax.plot(xp, y, color=color, linestyle='--', label=SPECIES_TEX[species.replace('_g', '')])
        ax.set(title=label, xlim=limits, xlabel='Position (cm)', ylabel='I yield per 1 cm bin (%)'); ax.legend()
    save(fig, 'baseline_profiles')
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 4.9), layout='constrained', sharey='row')
    colors = {1: '#9ecae1', 2: '#4292c6', 3: '#08519c', 4: '#08306b'}
    styles = {1: ':', 2: '--', 3: '-', 4: '-'}
    for row, (key, experiment, title, limits) in enumerate([
            ('pure', 'PbI2_SS', r'Pure PbI$_2$', (30, 43)),
            ('steel', 'LBE-I_SS_II', r'Steel Q2, 1%, shifted BiI', (39, 65))]):
        for col, scheme in enumerate(['up', 'll']):
            ax = axes[row, col]
            for level in families[f'{key}.{scheme}']['levels']:
                x, y, obs = scan_from_profile(study_folder(key, level, scheme), experiment)
                ax.plot(x, y, color=colors[level], linestyle=styles[level], lw=1.5 if level == 4 else 1.2,
                        label=f'{CELLS[level]:,} cells')
            ax.fill_between(x, obs, color='.85', step='mid', label='Published scan', zorder=0)
            ax.set(title=f"{title}; {scheme_names[scheme]}", xlim=limits, xlabel='Position (cm)')
            if col == 0:
                ax.set_ylabel('I per bin (% of scan)')
    handles, names = axes[1, 1].get_legend_handles_labels()
    fig.legend(handles, names, loc='outside upper center', ncol=3)
    save(fig, 'refinement')
    chemistry = []
    for name, r in original_metrics.items():
        if r['source_mode'] == 'Q2':
            chemistry.append({'case': name, 'experiment': r['experiment'], 'variant': r['thermo_variant'],
                              'beta': r['metal_saturation'], 'BiI3_pct': r['BiI3_iodide_I_pct']})
    # intermediate saturations 0.2--0.56 % are follow-up cases (18,432 cells, upwind); the share differs by
    # at most 0.03 percentage points between schemes for the three production saturations.
    for r in metrics:
        if r['category'] == 'source_loading':
            chemistry.append({'case': r['case'], 'experiment': 'LBE-I_SiO2_I' if 'SiO2' in r['case'] else 'LBE-I_SS_II',
                              'variant': 'shift29k' if 'shift29k' in r['case'] else 'baseline',
                              'beta': float(r['case'].split('__sat')[1]), 'BiI3_pct': r['BiI3_iodide_I_pct']})
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.9), layout='constrained', sharey=True)
    for ax, experiment, label in zip(axes, ['LBE-I_SS_II', 'LBE-I_SiO2_I'], ['Steel', 'Silica']):
        for variant, tag in [('baseline', 'Baseline thermodynamics'), ('shift29k', 'BiI +29 kJ/mol surrogate')]:
            selected = sorted([r for r in chemistry if r['experiment'] == experiment and r['variant'] == variant], key=lambda r: r['beta'])
            ax.plot([100 * r['beta'] for r in selected], [r['BiI3_pct'] for r in selected], 'o-', ms=4, label=tag)
        measured = next(r['observed_BiI3_iodide_I_pct'] for r in original_metrics.values() if r['experiment'] == experiment)
        ax.axhline(measured, color='.4', linestyle='--', label='Published iodide-I split')
        ax.set(xscale='log', xlabel='Metal saturation (%)', ylim=(0, 105), title=label)
    axes[0].set_ylabel(r'BiI$_3$ share of iodide I (%)')
    fig.legend(*axes[0].get_legend_handles_labels(), loc='outside upper center', ncol=2)
    save(fig, 'source_thermodynamics')
    fig = plt.figure(figsize=(WIDTH, 5.3), layout='constrained')
    grid = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    r = read_csv(BENCH / 'digitized/LBE-I_SiO2_I_temperature.csv')
    for shift, label in [(0, 'Published curve'), (-4.096367986, r'Curve $-$4.10 cm'), (-4.65314628, r'Curve $-$4.65 cm')]:
        axes[0].plot([100 * float(v['x_m']) + shift for v in r], [float(v['T_K']) - 273.15 for v in r], label=label)
    axes[0].scatter([51, 59], [306, 127], marker='x', s=45, color='black', label='Printed peak labels')
    axes[0].set(title='Temperature-coordinate hypotheses', xlim=(45, 68), ylim=(50, 500),
                xlabel='Position (cm)', ylabel='Temperature (°C)'); axes[0].legend()
    for ax, mode, base in zip(axes[1:], ['Q1', 'Q2'], ['LBE-I_SiO2_I_Q1_baseline', 'LBE-I_SiO2_I_Q2_baseline_sat0.001']):
        for suffix, label in [('', 'Original inputs'), ('curveShiftPbI2', r'Curve $-$4.10 cm'), ('curveShiftBiI3', r'Curve $-$4.65 cm')]:
            folder = PROD / 'profiles' / (base + '__' + suffix if suffix else base)
            x, y, obs = scan(folder, normalise=True); ax.plot(x, y, label=label)
        ax.fill_between(x, obs, color='.85', step='mid', label='Published scan', zorder=0)
        ax.set(title=mode + ('; 0.1% saturation' if mode == 'Q2' else ''), xlim=(40, 68),
               xlabel='Position (cm)', ylabel='I per bin (% of scan)')
    fig.legend(*axes[1].get_legend_handles_labels(), loc='outside lower center', ncol=2)
    save(fig, 'silica_alignment')
    shown = [r for r in sensitivity if r['category'] in ['transport', 'kinetics', 'source_geometry', 'flow_reference', 'shutdown']]
    fig, ax = plt.subplots(figsize=(4.6, 2.3), layout='constrained')
    ax.barh([r['label'] for r in shown], [r['max_L1_pct'] for r in shown], color='#2878a6')
    ax.invert_yaxis(); ax.set(xlabel='Maximum iodine profile change $E_{L1}$ (%)')
    save(fig, 'physical_sensitivity')
    fig, axes = plt.subplots(3, 2, figsize=(WIDTH, 6.9), layout='constrained')
    titles = {r['case']: f"{r['case']}: {r['flow_mL_min']} mL/min, {r['duration_s']} s" for r in lob['prod']}
    for ax, r in zip(axes.flat, lob['prod']):
        name = r['case']
        xe, ye = r['experiment']
        x, deposit = r['x_m'], r['deposit_mol_per_m']
        scale = 1e4   # mol/m to umol/cm
        ax.plot(100 * xe, scale * ye, 'o-', color='#1f5fa8', ms=3, lw=.8, label='Measured (digitised)')
        ax.step(100 * x, scale * deposit, where='mid', color='black', lw=1.3, label='This work (deposit)')
        for kind, color, label in [('temperature_based', '#2a9d3c', 'Temperature-based model'),
                                   ('vapour_pressure_based', '#e07b00', 'Vapour-pressure model')]:
            m = read_csv(LOB_BENCH / 'digitized' / f'{name}_model_{kind}.csv')
            ax.plot([100 * float(v['x_m']) for v in m], [scale * float(v['n_mol_per_m']) for v in m],
                    color=color, lw=1.0, label=label)
        support = np.concatenate([x[deposit > 0.01 * deposit.max()], xe[ye > 0.01 * ye.max()]])
        if name == 'T5_60':
            low = 0.0
        else:
            low = max(0.0, support.min() - 0.03)
        ax.set(xlim=(100 * low, 100 * min(0.69, support.max() + 0.03)), title=titles[name], xlabel='Position (cm)')
        ax.set_ylabel(r'PbI$_2$ ($\mu$mol/cm)')
        twin = ax.twinx(); twin.spines['right'].set_visible(True)
        twin.plot(100 * r['wall_T'][0], r['wall_T'][1], color='#c0392b', ls='--', lw=.9, label='Wall temperature')
        twin.set_ylabel('Wall T (K)', color='#c0392b'); twin.tick_params(axis='y', colors='#c0392b')
        twin.set_xlim(ax.get_xlim())
        lims = [i for i, v in enumerate(r['wall_T'][0]) if ax.get_xlim()[0] <= 100 * v <= ax.get_xlim()[1]]
        twin.set_ylim(r['wall_T'][1][lims].min() - 20, r['wall_T'][1][lims].max() + 20)
    axes.flat[-1].axis('off')
    handles = axes.flat[0].get_legend_handles_labels()[0] + [plt.Line2D([], [], color='#c0392b', ls='--', lw=.9)]
    names = axes.flat[0].get_legend_handles_labels()[1] + ['Wall temperature (right axis)']
    axes.flat[-1].legend(handles, names, loc='center', frameon=False)
    save(fig, 'lobresco_comparison')
    data['chemistry'] = chemistry
    input_files = [FULL / 'analysis.json', FOLLOW / 'analysis.json', PROD / 'analysis.json', PROD / 'ledger.json',
                   BENCH / 'experiments.csv', BENCH / 'targets.csv', BENCH / 'digitized/provenance.json',
                   STUDY / 'results.json', ITER / 'summary.json', LOB_BENCH / 'digitized/provenance.json']
    input_files += list((BENCH / 'digitized').glob('*.csv')) + list((LOB_BENCH / 'digitized').glob('*.csv'))
    input_files += [Path(__file__).resolve(), BENCH / 'analyse_matrix.py', BENCH / 'report.py',
                    BENCH / 'prepare.py', BENCH / 'carriers.py']
    for root in [FULL / 'profiles', FOLLOW / 'profiles', PROD / 'profiles', STUDY / 'profiles', ITER / 'profiles']:
        input_files += list(root.rglob('*.dat')) + list(root.rglob('*.csv'))
    input_files += list((LOB_RUN / 'cases').rglob('*.dat')) + list((LOB_RUN / 'cases').rglob('case.json'))
    input_files += list((LOB_RUN / 'inputs').glob('*.csv'))
    data['input_sha256'] = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(input_files)}
    dump(generated / 'results.json', data)
    dump(generated / 'artifact_sha256.json', {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in [generated, figures] for p in sorted(directory.iterdir())
        if p.is_file() and p.name != 'artifact_sha256.json'})
    print('Generated 8 figures, 6 tables and numerical provenance for 15 baseline + 64 follow-up cases, '
          '19 production cases, the scheme/mesh study, the iteration check and the Lobresco comparison.')


if __name__ == '__main__':
    main()
