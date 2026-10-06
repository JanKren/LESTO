#!/usr/bin/env python3
"""Digitise Liu's published curves; requires numpy/Pillow and pdfimages.

Calibrations refer to the embedded images, not a screenshot or resized page.
The CSVs are figure-derived observations, not the authors' numerical records.
"""
import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PANELS = {
    'PbI2_SS': dict(page=6, figure='2B', size=[969, 872],
                   frame=[113, 289, 850, 761], x_ticks=[204, 808.5],
                   T_ticks=[761, 341.5], gamma_ticks=[761, 341.5], gamma_max=40,
                   legend=[420, 610], deposits=[33, 42]),
    'LBE-I_SS_II': dict(page=11, figure='7B', size=[2032, 1451],
                      frame=[1093, 23, 1885, 531], x_ticks=[1189.5, 1839.5],
                      T_ticks=[531, 80], gamma_ticks=[531, 94.5], gamma_max=30,
                      legend=[1450, 340], deposits=[40, 63]),
    'LBE-I_SiO2_I': dict(page=11, figure='7C', size=[2032, 1451],
                       frame=[135, 625, 932, 1137], x_ticks=[226, 831],
                       T_ticks=[1137, 681.5], gamma_ticks=[1137, 696], gamma_max=30,
                       legend=[500, 950], deposits=[41, 65]),
}


def extract(image, calibration):
    a = np.asarray(image).astype(float)
    if list(image.size) != calibration['size']:
        raise ValueError('unexpected figure resolution; recalibrate before digitising')
    left, top, right, bottom = calibration['frame']
    red = (a[:, :, 0] > 155) & (a[:, :, 0] - a[:, :, 1] > 65) & (a[:, :, 0] - a[:, :, 2] > 65)
    lx, ly = calibration['legend']
    red[ly:, :lx] = False
    xpix, ypix = [], []
    for x in range(left + 5, right - 4):
        ys = np.flatnonzero(red[top + 4:bottom - 3, x]) + top + 4
        if len(ys):
            xpix.append(x)
            ypix.append(float(np.median(ys)))
    xpix, ypix = np.array(xpix), np.array(ypix)
    if len(xpix) < 200 or np.diff(xpix).max() > 20:
        raise ValueError('red temperature trace missing or interrupted')
    x10, x70 = calibration['x_ticks']
    cm = lambda xp: 10 + (xp - x10) * 60 / (x70 - x10)
    pix = lambda xc: x10 + (xc - 10) * (x70 - x10) / 60
    y0, y800 = calibration['T_ticks']
    # One-millimetre sampling is an interpolation of pixels, not measurement precision.
    xs = np.arange(np.ceil(cm(xpix[0]) * 10) / 10, cm(xpix[-1]), .1)
    temperatures = 273.15 + (y0 - np.interp(pix(xs), xpix, ypix)) * 800 / (y0 - y800)
    profile = [{'x_m': float(x / 100), 'T_K': float(t)} for x, t in zip(xs, temperatures)]
    low, high = calibration['deposits']
    gray = (a.max(2) - a.min(2) < 8) & (a.mean(2) >= 220) & (a.mean(2) <= 245)
    gy0, gymax = calibration['gamma_ticks']
    scan = []
    for xc in np.arange(low, high + .05, .1):
        xp = int(round(pix(xc)))
        tops = []
        for x in range(xp - 1, xp + 2):
            ys = np.flatnonzero(gray[top + 5:bottom - 3, x]) + top + 5
            # A bar must contain a vertical run, not isolated JPEG/text pixels.
            runs = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
            runs = [r for r in runs if len(r) >= 5]
            if runs:
                tops.append(float(runs[0][0]) - 1.5)
        height = 0 if not tops else (gy0 - np.median(tops)) * calibration['gamma_max'] / (gy0 - gymax)
        scan.append({'x_m': float(xc / 100), 'relative_yield_percent_per_cm': float(height)})
    # Black vertical outlines leave short holes in the gray-fill extraction.
    # Bridge only holes <= 0.4 cm; distinct deposits remain separated.
    heights = np.array([r['relative_yield_percent_per_cm'] for r in scan])
    zeros = np.flatnonzero(heights == 0)
    for gap in np.split(zeros, np.flatnonzero(np.diff(zeros) > 1) + 1):
        if len(gap) and len(gap) <= 4 and gap[0] > 0 and gap[-1] < len(heights) - 1:
            heights[gap] = np.interp(gap, [gap[0] - 1, gap[-1] + 1],
                                    heights[[gap[0] - 1, gap[-1] + 1]])
    # Read one height per printed 1 cm bin, using its interior to avoid outlines.
    # Do not integrate the missing gray pixels occupied by black outlines.
    binned = []
    for xc in range(low, high + 1):
        near = [h for r, h in zip(scan, heights)
                if abs(r['x_m'] * 100 - xc) < .001]
        binned.append(dict(x_m=xc / 100,
                           relative_yield_percent_per_cm=float(np.median(near)) if near else 0))
    scan = binned
    integral = float(sum(r['relative_yield_percent_per_cm'] for r in scan))
    if not 90 < integral < 110:
        raise ValueError(f'gamma histogram area {integral}%; inspect calibration')
    return profile, scan, integral


def main():
    from prepare import write_csv, load_references
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf', type=Path)
    parser.add_argument('--out', type=Path, default=HERE / 'digitized')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    provenance = dict(doi='10.1007/s10967-025-10246-4', licence='CC BY 4.0',
                      pdf_sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
                      method='red-pixel median and linear dash interpolation; gray histogram fill',
                      temperature_pixel_uncertainty_K=4,
                      position_pixel_uncertainty_m=.001,
                      gamma_height_pixel_uncertainty_percent=.4,
                      notes=['Axes calibrated manually using printed major ticks.',
                             'Gamma CSV reads the interiors of printed 1 cm bins, not raw scan data.',
                             'Histogram areas are reported without normalization.',
                             'Pixel height uncertainty excludes bin alignment, black outlines and arrow-obscuration error.',
                             'BiI3_SiO2_800 has no plotted T(x) in this paper.',
                             'Plot origin and unplotted column ends require author confirmation.'],
                      experiments={})
    observations = {}
    with tempfile.TemporaryDirectory(prefix='liu-figures-') as tmp:
        for name, c in PANELS.items():
            prefix = Path(tmp) / name
            subprocess.run(['pdfimages', '-f', str(c['page']), '-l', str(c['page']),
                            '-png', str(args.pdf), str(prefix)], check=True)
            path = next(Path(tmp).glob(name + '-*.png'))
            profile, scan, area = extract(Image.open(path).convert('RGB'), c)
            observations[name] = (profile, scan)
            write_csv(args.out / (name + '_temperature.csv'), profile, ['x_m', 'T_K'])
            write_csv(args.out / (name + '_gamma.csv'), scan,
                      ['x_m', 'relative_yield_percent_per_cm'])
            provenance['experiments'][name] = dict(calibration=c, gamma_area_percent=area,
                image_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                covered_x_m=[profile[0]['x_m'], profile[-1]['x_m']])
    _, targets, _ = load_references()
    zones = {'PbI2_SS': {'PbI2': [.34, .40]},
             'LBE-I_SS_II': {'PbI2': [.49, .54], 'BiI3': [.56, .61]},
             'LBE-I_SiO2_I': {'PbI2': [.49, .55], 'BiI3': [.57, .63]}}
    checks = []
    for name, (profile, scan) in observations.items():
        for species, (lo, hi) in zones[name].items():
            target = next(r for r in targets if r['experiment'] == name and r['species'] == species)
            peak = max((r for r in scan if lo <= r['x_m'] <= hi),
                       key=lambda r:r['relative_yield_percent_per_cm'])
            x = peak['x_m']
            temperature = float(np.interp(x, [r['x_m'] for r in profile], [r['T_K'] for r in profile]))
            stated = float(target['measured_peak_C']) + 273.15
            difference = temperature - stated
            uncertainty = float(target['uncertainty_K']) + 4
            checks.append(dict(experiment=name, species=species, scan_peak_x_m=x,
                figure_temperature_K=temperature, printed_peak_K=stated,
                difference_K=difference, printed_uncertainty_K=float(target['uncertainty_K']),
                status='CONSISTENT_WITHIN_UNCERTAINTY' if abs(difference) <= uncertainty
                       else 'FIGURE_LABEL_DISCREPANCY'))
    write_csv(args.out/'figure_checks.csv', checks, list(checks[0]))
    provenance['figure_checks'] = checks
    provenance['notes'].append('Figure 7C trace and printed peak labels disagree; do not fit the trace to the labels.')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), constrained_layout=True)
    for ax, (name, (profile, scan)) in zip(axes, observations.items()):
        ax.plot([r['x_m'] for r in profile], [r['T_K']-273.15 for r in profile], 'r--')
        ax.set(xlabel='Published column coordinate (m)', ylabel='Wall temperature (°C)', title=name)
        second = ax.twinx()
        second.bar([r['x_m'] for r in scan], [r['relative_yield_percent_per_cm'] for r in scan],
                   width=.01, color='0.85', edgecolor='0.4', alpha=.7)
        second.set_ylabel('Iodine yield per 1 cm bin (%)')
        for check in checks:
            if check['experiment'] == name:
                x = check['scan_peak_x_m']
                ax.plot(x, check['printed_peak_K']-273.15, 'bo')
                ax.annotate(f"{check['species']}: printed {check['printed_peak_K']-273.15:g} °C",
                            (x, check['printed_peak_K']-273.15), xytext=(-145, -20), textcoords='offset points', fontsize=8)
    fig.suptitle('Liu et al. 2025, DOI 10.1007/s10967-025-10246-4 (CC BY 4.0)\nFigure digitisation; blue markers = printed peak labels')
    fig.savefig(args.out/'figure_checks.png', dpi=140)
    plt.close(fig)
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({n: v['gamma_area_percent'] for n, v in provenance['experiments'].items()}, indent=2))


if __name__ == '__main__':
    main()
