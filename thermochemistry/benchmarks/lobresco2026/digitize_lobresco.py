#!/usr/bin/env python3
"""Digitise Lobresco et al. (2026) NED 459:115201, Figs. 13 and 14(i)-(iv).

Calibration from the embedded 300 ppi JPEGs (pdfimages -j), axis ticks/gridlines
detected in the image and fitted linearly.  Figure-derived observations only.

The paper PDF and the extracted page images are not part of the repository.
To rerun: put the open-access PDF (doi:10.1016/j.nucengdes.2026.115201) at
papers/Lobresco2026.pdf in the repository root, then
    mkdir -p pdf && pdfimages -j -f 9 -l 10 ../../../papers/Lobresco2026.pdf pdf/img
and run this script from this directory; it rewrites digitized/.
"""
import hashlib, json, csv, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage

G = Path(__file__).resolve().parent
PDF = G.parents[2] / 'papers' / 'Lobresco2026.pdf'
OUT = G / 'digitized'
IMG = G / 'pdf'

# Gridline / tick pixel positions were detected automatically (see detect())
# and are checked again below; values: (pixel, data value)
PANELS = {
 'T106_60': dict(image='img-001.jpg', figure='13', frame=(98, 36, 807, 525),
    legend=(112, 58, 182, 182), markers_cm=[(40, 60)]),
 'T106_6': dict(image='img-002.jpg', figure='14(i)', frame=(110, 184, 845, 621), legend=None, markers_cm=[(35, 55)]),
 'T106_300': dict(image='img-002.jpg', figure='14(ii)', frame=(1050, 184, 1785, 621), legend=None, markers_cm=[(45, 58)]),
 'T5_60': dict(image='img-002.jpg', figure='14(iii)', frame=(110, 807, 845, 1244), legend=None, markers_cm=[(2, 5), (15, 16), (28, 48)]),
 'T540_60': dict(image='img-002.jpg', figure='14(iv)', frame=(1050, 807, 1785, 1244), legend=None, markers_cm=[(45, 65)]),
}
TABLE2 = {'T106_6': (814, 775, 841), 'T106_60': (723, 713, 760), 'T106_300': (667, 666, 701),
          'T5_60': (847, 787, 797), 'T540_60': (700, 670, 685)}
TABLE1 = {'T106_6': (106, 6, 2.23, 4.84, 775), 'T106_60': (106, 60, 2.04, 4.43, 712),
          'T106_300': (106, 300, 1.37, 2.97, 666), 'T5_60': (5, 60, 1.45, 3.15, 788),
          'T540_60': (540, 60, 1.81, 3.93, 671)}


def runs(idx):
    out = []
    for i in idx:
        if out and i - out[-1][-1] <= 1:
            out[-1].append(i)
        else:
            out.append([i])
    return [float(np.mean(r)) for r in out]


def detect(a, frame):
    l, t, r, b = frame
    gray = (a.max(2) - a.min(2) < 15) & (a.mean(2) > 150) & (a.mean(2) < 225)
    v = runs(np.flatnonzero(gray[t + 5:b - 5, :].sum(0) > 0.45 * (b - t)))
    v = [x for x in v if l + 5 < x < r - 5]
    h = runs(np.flatnonzero(gray[:, l + 5:r - 5].sum(1) > 0.4 * (r - l)))
    h = [y for y in h if t + 3 < y < b - 3]
    pink = (a[:, :, 0] > 180) & (a[:, :, 0] - a[:, :, 2] > 30) & (a[:, :, 0] - a[:, :, 1] > 30)
    ticks = runs(np.flatnonzero(pink[:, r + 2:r + 10].sum(1) >= 6))
    ticks = [y for y in ticks if t - 5 < y < b + 5]
    return v, h, ticks


def fit(pix, val):
    A = np.vstack([pix, np.ones_like(pix)]).T
    coef, res, *_ = np.linalg.lstsq(A, val, rcond=None)
    resid = val - A @ coef
    return coef, float(np.abs(resid).max())


def trace(mask, frame, colmin=1):
    l, t, r, b = frame
    xs, ys = [], []
    for x in range(l + 3, r - 2):
        yy = np.flatnonzero(mask[t + 2:b - 1, x]) + t + 2
        if len(yy) >= colmin:
            xs.append(x); ys.append(float(np.median(yy)))
    return np.array(xs, float), np.array(ys, float)


def main():
    OUT.mkdir(exist_ok=True)
    prov = dict(source='Lobresco et al. (2026), Nucl. Eng. Des. 459, 115201',
                pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(),
                extraction='pdfimages -j -f 9 -l 10 (embedded 300 ppi JPEG, no resampling)',
                method=('axis calibration: least-squares line through detected gridlines (x, left y) '
                        'and right-spine major ticks (T); wall T = median of salmon dashed-curve pixels per '
                        'image column, gaps between dashes linearly interpolated; experiment = centroids of blue '
                        'circular markers after binary opening (removes the dash-dot line); band = vertical extent '
                        'of light-blue fill at each marker column (marker pixels bridged); model curves = median of '
                        'green/orange pixels per column; arrow = column of the black arrow shaft.'),
                table1_verbatim={k: dict(flow_mL_min=v[0], duration_s=v[1], mass_mg=v[2],
                                         moles_1e_6_mol=v[3], T_dep_calc_K=v[4]) for k, v in TABLE1.items()},
                table2_verbatim={k: dict(experiment_K=v[0], temperature_based_K=v[1],
                                         vapour_pressure_based_K=v[2]) for k, v in TABLE2.items()},
                stated_experimental_uncertainty='paper Sec. 2.4: +-15 K temperature, +-5 % He flow, '
                    '+-0.05 mg sample mass, +-1 cm axial position',
                panels={},
                estimated_digitisation_error=dict(
                    x_mm='+-1 (calibration residual <= 0.44 mm, 1.03-1.07 mm/px; nonzero markers sit 0.1-1.7 mm right of the nominal cm, systematic +0.6 mm)',
                    wall_T_K='+-2 to 3 (1.68-1.90 K/px, tick-fit residual <= 1.44 K) plus slope x position error: about +-1.5 K/mm in the deposition zone; T106_60 (Fig. 13) curve hidden by the legend at 0.105-0.135 m and linearly bridged there (do not use it in that range)',
                    deposit_mol_per_m='+-3e-7 (one marker radius is about 6 px = 7.5e-7; 1.0-1.25e-7 per px)',
                    band='light-blue fill at the marker column only; where the band is not wider than the marker the counting uncertainty is unresolved (band_note); clipped at the frame top for T106_6 0.40 m and T106_300 0.49-0.50 m'),
                notes=['Zero-valued points are assigned at the nominal 1 cm positions of the plotted marker ranges (markers_cm); only the nonzero markers are located individually.',
                       'T5_60 scan points exist at 0.02-0.05, 0.15-0.16 and 0.28-0.48 m only (all plotted markers); there is no plotted point at 0.06-0.14 m.',
                       'Model curves (green, orange) are digitised for reference only; where curves overlap (orange/green/blue) the per-column median can follow the other curve.',
                       'The integral of each digitised experimental profile is 0.994-0.998 n0 (Table 1 moles): the paper normalised the scan to n0.'])
    checks = []
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for name, c in PANELS.items():
        im = Image.open(IMG / c['image']).convert('RGB')
        a = np.asarray(im).astype(int)
        l, t, r, b = c['frame']
        v, h, ticks = detect(a, c['frame'])
        def regular(lines, first_value, step, first_px=None):
            lines = sorted(lines)
            d = np.diff(lines); sp = np.median(d[d > 0.6 * np.median(d)])
            ref = lines[0] if first_px is None else first_px
            idx = np.round((np.array(lines) - ref) / sp)
            good = np.abs((np.array(lines) - ref) / sp - idx) < 0.06
            return list(np.array(lines)[good]), list(first_value + step * idx[good])
        v, xvals = regular(v, 0.0, 0.1)
        assert len(v) >= 7 and abs(xvals[0]) < 1e-12, (name, v, xvals)
        xcoef, xres = fit(np.array(v), np.array(xvals))
        top_val = 1.6e-4 if name == 'T106_60' else 1.0e-4
        h, vals = regular(h, top_val, -0.2e-4)
        ycoef, yres = fit(np.array(h), np.array(vals))
        top_T = 1100 if name == 'T106_60' else 1000
        dT = -100 if name == 'T106_60' else -200
        ticks, tvals = regular(ticks, top_T, dT)
        tcoef, tres = fit(np.array(ticks), np.array(tvals, float))
        X = lambda px: xcoef[0] * px + xcoef[1]
        Y = lambda py: ycoef[0] * py + ycoef[1]
        TT = lambda py: tcoef[0] * py + tcoef[1]
        Xinv = lambda xm: (xm - xcoef[1]) / xcoef[0]
        R, Gc, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]
        legend = np.zeros(a.shape[:2], bool)
        if c['legend']:
            L0, T0, L1, T1 = c['legend']; legend[T0:T1, L0:L1] = True
        # --- wall temperature: salmon dashed
        salmon = (R > 215) & (R - Gc > 50) & (np.abs(Gc - B) < 25) & (Gc > 95) & (Gc < 205) & ~legend
        # track the dashed curve column by column from its left end; a column joins the
        # trace only within 4 px of the running extrapolation (crossing curves rejected)
        xs, ys = [], []
        for x in range(l + 3, r - 2):
            yy = np.flatnonzero(salmon[t + 2:b - 1, x]) + t + 2
            if not len(yy):
                continue
            if len(ys) < 3:
                if len(ys) == 0 and x > l + 60:
                    continue
                xs.append(x); ys.append(float(np.median(yy))); continue
            # linear extrapolation from the last accepted points
            k0 = max(0, len(xs) - 6)
            slope = (ys[-1] - ys[k0]) / max(1, xs[-1] - xs[k0])
            pred = ys[-1] + slope * (x - xs[-1])
            tol = 4 + 0.15 * (x - xs[-1])
            near = yy[np.abs(yy - pred) <= tol]
            if len(near):
                xs.append(x); ys.append(float(np.median(near)))
            elif x - xs[-1] > 25:
                xs.append(x); ys.append(float(np.median(yy)))
        xs, ys = np.array(xs, float), np.array(ys, float)
        xm_all = np.arange(np.ceil(X(xs[0]) * 1000) / 1000, X(xs[-1]) + 1e-9, 0.001)
        T_all = TT(np.interp(Xinv(xm_all), xs, ys))
        gaps = np.diff(xs)
        # --- experiment markers
        blue = (B > 180) & (R < 90) & (Gc < 90) & ~legend
        inner = np.zeros_like(blue); inner[t + 1:b + 8, l + 1:r - 1] = True
        blue &= inner
        green_m = (Gc > 90) & (R < 80) & (B < 80) & (Gc - R > 50)
        orange_m = (R > 220) & (Gc > 90) & (Gc < 170) & (B < 70)
        disk = (np.hypot(*np.mgrid[-5:6, -5:6]) <= 5.0).astype(float)
        area = disk.sum()
        sb = ndimage.convolve(blue.astype(float), disk, mode='constant') / area
        colored = (a.max(2) - a.min(2) > 90)
        sa = ndimage.convolve(colored.astype(float), disk, mode='constant') / area
        y0px = (0.0 - ycoef[1]) / ycoef[0]
        # Markers lie on the 1 cm grid.  Every 1 cm position of the plotted marker range
        # gets the best disk match above y = 0, or 0 when no marker is found there.
        pts = []
        for k in [k for lo_, hi_ in c['markers_cm'] for k in range(lo_, hi_ + 1)]:
            xm0 = k / 100.0
            px = int(round(Xinv(xm0)))
            top_ = t + 2; bot_ = int(round(y0px)) - 8
            win_sa = sa[top_:bot_, px - 3:px + 4]
            win_sb = sb[top_:bot_, px - 3:px + 4]
            score = np.where(win_sb >= 0.2, win_sa + 0.5 * win_sb, -1)
            iy, ix = np.unravel_index(np.argmax(score), score.shape)
            if win_sa[iy, ix] >= 0.85 and win_sb[iy, ix] >= 0.3:
                cy0, cx0 = iy + top_, ix + px - 3
                yy, xx = np.mgrid[cy0 - 6:cy0 + 7, cx0 - 6:cx0 + 7]
                inside = (np.hypot(yy - cy0, xx - cx0) <= 6) & colored[cy0 - 6:cy0 + 7, cx0 - 6:cx0 + 7]
                cy, cx = yy[inside].mean(), xx[inside].mean()
                if False: print('   cand', name, k, 'sa %.2f sb %.2f' % (win_sa[iy, ix], win_sb[iy, ix]), 'y %.2e off %.2f' % (Y(cy), (X(cx)-xm0)*1000))
                pts.append((X(cx), Y(cy), cx, cy, int(inside.sum())))
            else:
                pts.append((xm0, 0.0, Xinv(xm0), float(y0px), 0))
        # band (light blue fill)
        band = (B > 225) & (R > 150) & (R < 235) & (Gc > 150) & (Gc < 235) & (B - R > 15) & ~legend
        exp_rows = []
        for xm, ym, cx, cy, s in pts:
            col = int(round(cx))
            cols = [col]
            lo = hi = None
            mask = band[t - 2:b + 2, :] | (blue[t - 2:b + 2, :])
            yy = np.flatnonzero(mask[:, cols].any(1)) + t - 2
            # contiguous run containing marker centre
            if len(yy):
                rr = np.split(yy, np.flatnonzero(np.diff(yy) > 2) + 1)
                for run_ in rr:
                    if run_[0] - 1 <= cy <= run_[-1] + 1:
                        hi, lo = Y(run_[0]), Y(run_[-1])
            bandw = None if hi is None else (hi - lo)
            exp_rows.append(dict(x_m=round(xm, 4), x_bin_m=round(round(xm * 100) / 100, 2),
                                 x_offset_mm=round((xm - round(xm * 100) / 100) * 1000, 2),
                                 n_mol_per_m=max(ym, 0.0), band_low_mol_per_m=lo, band_high_mol_per_m=hi,
                                 band_clipped_at_frame_top=bool(hi is not None and hi >= Y(t + 1) - 1e-12),
                                 marker_px=int(s)))
        # band only meaningful if wider than marker diameter (~11 px)
        mpx = abs(ycoef[0]) * 12
        for row in exp_rows:
            if row['band_high_mol_per_m'] is not None and row['band_high_mol_per_m'] - row['band_low_mol_per_m'] <= mpx:
                row['band_note'] = 'band not wider than marker; uncertainty below marker size'
        # --- models
        green = (Gc > 90) & (R < 80) & (B < 80) & (Gc - R > 50) & ~legend
        orange = (R > 220) & (Gc > 90) & (Gc < 170) & (B < 70) & ~legend
        models = {}
        for key, m in [('vapour_pressure_based', green), ('temperature_based', orange)]:
            mx, my = trace(m, c['frame'])
            models[key] = [dict(x_m=round(X(x), 5), n_mol_per_m=max(Y(y), 0.0)) for x, y in zip(mx, my)]
        # --- arrow: black pixels inside frame, away from spines
        black = (a.max(2) < 70)
        sub = black[b - 130:b - 3, l + 6:r - 6]
        colsum = sub.sum(0)
        acol = np.flatnonzero(colsum >= 25)
        arrow_x = X(np.mean(acol) + l + 6) if len(acol) else None
        arrow_rows = np.flatnonzero(sub[:, acol].any(1)) + b - 130 if len(acol) else []
        # --- write
        def wcsv(fn, rows, cols):
            with open(OUT / fn, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n', extrasaction='ignore')
                w.writeheader(); w.writerows(rows)
        wcsv(f'{name}_wall_temperature.csv', [dict(x_m=round(x, 4), T_K=round(T, 2)) for x, T in zip(xm_all, T_all)], ['x_m', 'T_K'])
        wcsv(f'{name}_experiment.csv', exp_rows, ['x_bin_m', 'x_m', 'x_offset_mm', 'n_mol_per_m',
             'band_low_mol_per_m', 'band_high_mol_per_m', 'band_clipped_at_frame_top', 'band_note'])
        for key, rows in models.items():
            wcsv(f'{name}_model_{key}.csv', rows, ['x_m', 'n_mol_per_m'])
        # --- T_dep check (paper rule: first position from inlet where deposit > 1% of peak, lin. interp.)
        e = sorted(exp_rows, key=lambda q: q['x_bin_m'])
        xe = np.array([q['x_bin_m'] for q in e]); ye = np.array([q['n_mol_per_m'] for q in e])
        peak = ye.max(); thr = 0.01 * peak
        k = int(np.flatnonzero(ye > thr)[0])
        if k > 0:
            x_on = xe[k - 1] + (thr - ye[k - 1]) / (ye[k] - ye[k - 1]) * (xe[k] - xe[k - 1])
        else:
            x_on = xe[0]
        T_on = float(np.interp(x_on, xm_all, T_all))
        T_first = float(np.interp(xe[k], xm_all, T_all))
        T_prev = float(np.interp(xe[k - 1], xm_all, T_all)) if k > 0 else None
        integral = float(np.trapz(ye, xe))  # mol over the 1 cm points
        sum_bins = float(ye.sum() * 0.01)
        # model T_dep from digitised model curves (reference only)
        mt = {}
        for key, rows in models.items():
            xx = np.array([q['x_m'] for q in rows]); yv = np.array([q['n_mol_per_m'] for q in rows])
            # restrict to x>0.2 except T5_60 green has an upstream spike
            pk = yv.max(); kk = np.flatnonzero(yv > 0.01 * pk)
            mt[key] = dict(x_onset_m=float(xx[kk[0]]), T_K=float(np.interp(xx[kk[0]], xm_all, T_all)),
                           peak_x_m=float(xx[np.argmax(yv)]), peak_mol_per_m=float(pk))
        chk = dict(case=name, peak_mol_per_m=float(peak), x_peak_m=float(xe[np.argmax(ye)]),
                   x_onset_m=float(x_on), T_dep_digitised_K=T_on,
                   T_at_first_bin_above_K=T_first, T_at_last_bin_below_K=T_prev,
                   T_dep_table2_K=TABLE2[name][0], difference_K=T_on - TABLE2[name][0],
                   exp_integral_trapz_mol=integral, exp_sum_bins_mol=sum_bins,
                   table1_moles=TABLE1[name][3] * 1e-6,
                   integral_over_n0=sum_bins / (TABLE1[name][3] * 1e-6),
                   model_onsets_digitised_reference_only=mt,
                   model_T_dep_table2_K=dict(temperature_based=TABLE2[name][1], vapour_pressure_based=TABLE2[name][2]))
        checks.append(chk)
        prov['panels'][name] = dict(figure=c['figure'], image=c['image'],
            image_sha256=hashlib.sha256((IMG / c['image']).read_bytes()).hexdigest(),
            image_size=list(im.size), frame_px=c['frame'],
            x_gridlines_px=v, x_m_per_px=xcoef[0], x_fit_max_residual_m=xres,
            y_gridlines_px=h, y_values=vals, y_per_px=ycoef[0], y_fit_max_residual=yres,
            T_ticks_px=ticks, T_values=tvals, T_K_per_px=tcoef[0], T_fit_max_residual_K=tres,
            wall_T_coverage_m=[float(xm_all[0]), float(xm_all[-1])],
            wall_T_max_dash_gap_px=float(gaps.max()), wall_T_max_dash_gap_mm=float(gaps.max() * xcoef[0] * 1000),
            arrow_x_m=arrow_x, arrow_tip_px=int(max(arrow_rows)) if len(arrow_rows) else None,
            n_markers=len(exp_rows), max_marker_x_offset_mm=float(max(abs(q['x_offset_mm']) for q in exp_rows)))
        # overlay
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.imshow(a.astype(np.uint8))
        ax.plot(Xinv(xm_all), (T_all - tcoef[1]) / tcoef[0], 'c-', lw=0.8, label='digitised wall T')
        ax.plot([Xinv(q['x_m']) for q in exp_rows], [(q['n_mol_per_m'] - ycoef[1]) / ycoef[0] for q in exp_rows],
                'kx', ms=6, label='digitised markers')
        for q in exp_rows:
            if q['band_high_mol_per_m'] is not None:
                ax.plot([Xinv(q['x_m'])] * 2, [(q['band_low_mol_per_m'] - ycoef[1]) / ycoef[0],
                         (q['band_high_mol_per_m'] - ycoef[1]) / ycoef[0]], 'm-', lw=0.8)
        for key, col in [('vapour_pressure_based', 'k:'), ('temperature_based', 'k--')]:
            rows = models[key]
            ax.plot([Xinv(q['x_m']) for q in rows], [(q['n_mol_per_m'] - ycoef[1]) / ycoef[0] for q in rows], col, lw=0.6)
        if arrow_x is not None:
            ax.axvline(Xinv(arrow_x), color='y', lw=0.6)
        ax.axvline(Xinv(x_on), color='r', lw=0.6)
        ax.set_xlim(l - 40, r + 60); ax.set_ylim(b + 40, t - 40)
        ax.set_title(f'{name}: digitised (overlay) vs embedded image; red line = 1 %-of-peak onset')
        ax.legend(fontsize=7, loc='upper right'); ax.axis('off')
        fig.savefig(OUT / f'{name}_overlay.png', dpi=160, bbox_inches='tight'); plt.close(fig)
    prov['tdep_checks'] = checks
    (OUT / 'provenance.json').write_text(json.dumps(prov, indent=2, default=float) + '\n')
    for q in checks:
        print(f"{q['case']:9s} onset x={q['x_onset_m']:.4f} T={q['T_dep_digitised_K']:.1f} "
              f"(table {q['T_dep_table2_K']}, diff {q['difference_K']:+.1f}; bins {q['T_at_last_bin_below_K']} / {q['T_at_first_bin_above_K']:.1f}) "
              f"peak {q['peak_mol_per_m']:.3e} at {q['x_peak_m']:.2f}  sum/n0 {q['integral_over_n0']:.3f}")
    for k, p in prov['panels'].items():
        print(k, 'xres', p['x_fit_max_residual_m'], 'yres', p['y_fit_max_residual'], 'Tres', p['T_fit_max_residual_K'],
              'arrow', p['arrow_x_m'], 'nmark', p['n_markers'], 'maxoff', p['max_marker_x_offset_mm'],
              'Tcov', p['wall_T_coverage_m'], 'gap', p['wall_T_max_dash_gap_mm'], 'K/px', p['T_K_per_px'], 'mm/px', p['x_m_per_px']*1000)


if __name__ == '__main__':
    main()
