#!/usr/bin/env python3
"""Figures of the partner-update deck, drawn at their size on the slide.

Every figure is drawn from a tracked repository record (manuscript/data,
manuscript/generated/results.json, run/024 production records).  The
column figures are 3.05 x 1.85 in, the width of a 0.55 slide column, so
that they are placed unscaled and their text stays readable;
liu_profiles spans the text width.  figure_sources.json lists the files
each figure reads with their sha256.  No result is typed in here.

    MPLCONFIGDIR=/tmp/lesto-matplotlib python3 make_figures.py
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FIG = HERE / "figures"
MS = REPO / "manuscript"
PROD = REPO / "run" / "024-liu2025" / "production_576x72_limitedLinear" / "profiles"

# Same categorical order as the case-ladder deck
BLUE, ORANGE, TEAL, PURPLE = "#2a3cf0", "#e08a00", "#00a3a3", "#6e14dc"
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#d9d9d9"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.3, "legend.frameon": False, "savefig.dpi": 300,
    "figure.dpi": 150, "pdf.fonttype": 42,
})

COL = (3.05, 1.72)          # inches: one 0.55 column of the 16:9 slide
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def used(name, *paths, note=""):
    entry = SOURCES.setdefault(name, {"inputs": {}, "note": note})
    for p in paths:
        entry["inputs"][str(Path(p).relative_to(REPO))] = sha(p)


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", metadata={"CreationDate": None})
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan
    return {k: np.array([num(r[k]) for r in rows]) for k in rows[0]}


def results():
    return json.loads((MS / "generated" / "results.json").read_text())


def graetz():
    """Relative error of the Graetz-Robin decay rate against radial cells."""
    root = MS / "data" / "graetz"
    fits = json.loads((root / "fits.json").read_text())
    ref = {r["Bi"]: r["m"] / 0.0024 for r in json.loads((root / "reference.json").read_text())}
    fig, ax = plt.subplots(figsize=COL, layout="constrained")
    for law, color, label in (("halfCell", BLUE, "half-cell wall law"),
                              ("none", ORANGE, "law at the cell centre")):
        for bi, marker in ((1.0, "o"), (50.0, "s")):
            pts = sorted((f["nr"], abs(f["mu"] / ref[bi] - 1)) for f in fits
                         if f["scheme"] == "linear" and f["outlet"] == "zeroGradient"
                         and f["Bi"] == bi and f["law"] == law)
            n, e = np.array(pts).T
            ax.loglog(n, e, marker=marker, color=color, ms=3.5,
                      label=label if bi == 1.0 else None)
    n = np.array([32, 256])
    ax.loglog(n, 3e-3 * (n / 32) ** -1.0, color=MUTED, ls=":", lw=1)
    ax.loglog(n, 1.5e-6 * (n / 32) ** -2.0, color=MUTED, ls=":", lw=1)
    ax.text(120, 1.3e-4, "order 1", color=MUTED)
    ax.text(150, 1.2e-8, "order 2", color=MUTED)
    ax.set_ylim(4e-10, 6e-2)
    ax.set_xlabel("radial cells")
    ax.set_ylabel("error of the decay rate")
    ax.set_xticks([32, 64, 128, 256], ["32", "64", "128", "256"])
    ax.minorticks_off()
    ax.legend(loc="lower left", fontsize=7)
    save(fig, "graetz")
    used("graetz", root / "fits.json", root / "reference.json",
         note="Graetz-Robin records (manuscript Sec. 5.1): linear advection, zeroGradient outlet; circles Bi = 1, squares Bi = 50.")


def damkohler():
    """Normalised wall flux against Da, with both limits and the column bands."""
    f = MS / "data" / "damkohler" / "damkohler_sweep.csv"
    d = read_csv(f)
    da = d["Da_561K"]
    plateau = np.pi ** 2 / 4
    g = np.logspace(-3.2, 3.6, 200)
    fig, ax = plt.subplots(figsize=COL, layout="constrained")
    for lo, hi, lab in ((7, 15, "0.01"), (76, 150, "0.1"), (1.4e3, 2.9e3, "1")):
        ax.axvspan(lo, hi, color=GRID, alpha=0.8, lw=0)
        ax.text(np.sqrt(lo * hi), 1.2e-3, r"$\alpha$=" + lab, rotation=90, ha="center",
                va="bottom", color=MUTED, fontsize=7)
    ax.loglog(g, g, color=MUTED, ls=":", lw=1, label="kinetic limit")
    ax.axhline(plateau, color=MUTED, ls="--", lw=1, label=r"diffusion limit $\pi^2/4$")
    ax.loglog(g, 1 / (1 / g + 1 / plateau), color=ORANGE, lw=1, label="series resistances")
    ax.loglog(da, d["Jhat_84p5"], "o", color=BLUE, ms=3.5, label="simulations")
    ax.set_ylim(5e-4, 8)
    ax.set_xlabel("wall Damköhler number Da")
    ax.set_ylabel("normalised wall flux")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=2,
               fontsize=6.8, handlelength=1.6, columnspacing=1.0)
    save(fig, "damkohler")
    used("damkohler", f, note="Damkoehler sweep of run/016 (manuscript Sec. 4.2): flux at x = 84.5 mm against Da at the onset temperature; grey: tube Da of the column cases for alpha = 0.01, 0.1, 1 (manuscript Sec. 4.2 text).")


def sensitivity():
    """Iodine profile change for diffusivity and accommodation variations."""
    f = REPO / "run" / "024-liu2025" / "production_sensitivities" / "comparison.csv"
    with open(f) as fh:
        rows = list(csv.DictReader(fh))
    order = ["D5e-05", "D0.0002", "alpha0.1", "alpha0.01"]
    names = {"D5e-05": "D / 2", "D0.0002": "2 D", "alpha0.1": "α = 0.1", "alpha0.01": "α = 0.01"}
    fig, axes = plt.subplots(1, 2, figsize=COL, sharey=True, layout="constrained")
    for ax, base, title in ((axes[0], "PbI2_SS_Q1_baseline", "pure PbI$_2$"),
                            (axes[1], "LBE-I_SS_II_Q2_baseline_sat0.01", "LBE, steel")):
        r = {x["variant"]: x for x in rows if x["base"] == base}
        y = np.arange(len(order))
        ax.barh(y, [float(r[v]["I_L1_pct"]) for v in order], 0.6,
                color=[BLUE, BLUE, ORANGE, ORANGE])
        ax.set_yticks(y, [names[v] for v in order])
        ax.set_title(title, color=INK, pad=2)
        ax.set_xlabel("profile change [%]")
        ax.grid(axis="y", visible=False)
    axes[0].invert_yaxis()
    save(fig, "sensitivity")
    used("sensitivity", f, note="Production-scheme sensitivities of run/024 (41,472 cells, limitedLinear), each against its own baseline; blue: diffusivity, orange: accommodation.")


def lobresco():
    """Onset-temperature deviation from the measurement, five experiments."""
    f = MS / "generated" / "results.json"
    res = json.loads(f.read_text())["lobresco2026"]
    rows, summ = res["rows"], res["summary"]
    series = (("this work", "ours_minus_exp_K", BLUE, summ["rms_ours_K"]),
              ("paper, T-based", "temperature_based_minus_exp_K", ORANGE, summ["rms_temperature_based_K"]),
              (r"paper, $p_v$-based", "vapour_pressure_based_minus_exp_K", TEAL, summ["rms_vapour_pressure_based_K"]))
    fig, ax = plt.subplots(figsize=COL, layout="constrained")
    x = np.arange(len(rows))
    for i, (label, key, color, rms) in enumerate(series):
        ax.bar(x + (i - 1) * 0.27, [r[key] for r in rows], 0.27, color=color,
               label=f"{label} (RMS {rms:.0f} K)")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xticks(x, [r["case"] for r in rows], fontsize=6.8, rotation=20)
    ax.set_ylabel("model $-$ measured [K]")
    ax.set_ylim(-65, 45)
    ax.grid(axis="x", visible=False)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=2,
               fontsize=6.8, columnspacing=0.8, handlelength=1.2)
    save(fig, "lobresco_onset")
    used("lobresco_onset", f, note="Manuscript Table V (results.json, lobresco2026): onset by the 1 % rule on the main deposit; published model values from Lobresco et al. (2026), Table 2.")


def liu():
    """Iodine wall profiles of the two steel configurations of Liu et al. (2025)."""
    panels = (("PbI2_SS_Q1_baseline", "pure PbI$_2$, steel", (30, 45)),
              ("LBE-I_SS_II_Q1_baseline", "iodine from LBE, steel", (38, 64)))
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 1.72), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1.4]})
    files = []
    for ax, (case, title, xlim) in zip(axes, panels):
        f = PROD / case / "iodine-comparison.csv"
        d = read_csv(f)
        x = d["x_m"] * 100
        ax.bar(x, np.nan_to_num(d["figure_relative_yield_percent"]), width=0.92,
               color=GRID, edgecolor=MUTED, lw=0.3, label="measured")
        ax.plot(x, d["predicted_relative_yield_percent"], color=BLUE,
                drawstyle="steps-mid", label="this work")
        ax.set_xlim(*xlim)
        ax.set_title(title, color=INK, pad=2)
        ax.set_xlabel("position [cm]")
        files.append(f)
    axes[0].set_ylabel("iodine [% per cm]")
    axes[1].legend(loc="upper left")
    save(fig, "liu_profiles")
    used("liu_profiles", *files,
         note="Production runs of run/024 (576 x 72, limitedLinear, steel Q1 baselines); measured scan digitised from Liu et al. (2025) as in the manuscript.")


def silica():
    """Silica Q1 iodine profile with the published and the two shifted temperature curves."""
    rows = {r["case"]: r for r in results()["alignment"] if r["mode"] == "Q1"}
    fig, ax = plt.subplots(figsize=COL, layout="constrained")
    files = []
    first = True
    for case, color, ls in (("LBE-I_SiO2_I_Q1_baseline", BLUE, "-"),
                            ("LBE-I_SiO2_I_Q1_baseline__curveShiftPbI2", ORANGE, "-"),
                            ("LBE-I_SiO2_I_Q1_baseline__curveShiftBiI3", TEAL, "--")):
        f = PROD / case / "iodine-comparison.csv"
        d = read_csv(f)
        x = d["x_m"] * 100
        if first:
            ax.bar(x, np.nan_to_num(d["figure_relative_yield_percent"]), width=0.92,
                   color=GRID, edgecolor=MUTED, lw=0.3, label="measured")
            first = False
        shift = rows[case]["shift_cm"]
        lab = "published curve" if shift == 0 else f"curve shifted {shift:.2f} cm"
        ax.plot(x, d["predicted_relative_yield_percent"], color=color, ls=ls,
                drawstyle="steps-mid", label=f"{lab}: {rows[case]['overlap_pct']:.0f} % overlap")
        files.append(f)
    ax.set_xlim(40, 68)
    ax.set_ylim(0, 58)
    ax.set_xlabel("position [cm]")
    ax.set_ylabel("iodine [% per cm]")
    ax.legend(loc="upper left", fontsize=6.8)
    save(fig, "silica")
    used("silica", *files, MS / "generated" / "results.json",
         note="Production runs of run/024, silica Q1, published and shifted temperature curves; overlap from results.json (alignment); measured scan digitised from Liu et al. (2025).")


def thermodynamics():
    """BiI3 share of iodide iodine against metal loading, two thermodynamic branches."""
    f = MS / "generated" / "results.json"
    fa = PROD.parent / "analysis.json"
    chem = results()["chemistry"]
    observed = {r["experiment"]: r["observed_BiI3_iodide_I_pct"]
                for r in json.loads(fa.read_text())["case_metrics"] if r["case"].endswith("Q1_baseline")}
    fig, axes = plt.subplots(1, 2, figsize=COL, sharey=True, layout="constrained")
    for ax, exp, title in ((axes[0], "LBE-I_SS_II", "steel"), (axes[1], "LBE-I_SiO2_I", "silica")):
        for variant, color, lab in (("baseline", BLUE, "baseline data"),
                                    ("shift29k", ORANGE, "BiI(g) +29 kJ/mol")):
            pts = sorted((100 * r["beta"], r["BiI3_pct"]) for r in chem
                         if r["experiment"] == exp and r["variant"] == variant)
            b, v = np.array(pts).T
            ax.semilogx(b, v, "o-", color=color, ms=3, label=lab)
        ax.axhline(observed[exp], color=MUTED, ls="--", lw=1, label="published split")
        ax.set_title(title, color=INK, pad=2)
        ax.set_xlabel("metal saturation [%]")
        ax.set_ylim(0, 105)
    axes[0].set_ylabel("BiI$_3$ share of iodine [%]")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=3, fontsize=6.5,
               columnspacing=0.8, handlelength=1.4)
    save(fig, "thermodynamics")
    used("thermodynamics", f, fa, note="results.json (chemistry): Q2 cases; dashed: published iodide-I split (analysis.json of the production runs); the +29 kJ/mol branch is a surrogate perturbation, not a qualified dataset (manuscript Fig. thermo).")


def temperatures():
    """Digitised wall-temperature curves of the three Liu et al. (2025) configurations."""
    bench = REPO / "thermochemistry" / "benchmarks" / "liu2025" / "digitized"
    fig, ax = plt.subplots(figsize=COL, layout="constrained")
    files = []
    for name, label, color in (("PbI2_SS", "pure PbI$_2$, steel", BLUE),
                               ("LBE-I_SS_II", "LBE, steel", ORANGE),
                               ("LBE-I_SiO2_I", "LBE, silica", TEAL)):
        f = bench / f"{name}_temperature.csv"
        d = read_csv(f)
        ax.plot(100 * d["x_m"], d["T_K"] - 273.15, color=color, label=label)
        files.append(f)
    # source region x = 0.09-0.105 m (manuscript Sec. 4.1)
    ax.axvspan(9, 10.5, color=ORANGE, alpha=0.18, lw=0, label="source region")
    ax.set_xlim(0, 80)
    ax.set_xlabel("position along the column [cm]")
    ax.set_ylabel("wall temperature [°C]")
    ax.legend(loc="upper right", fontsize=7)
    save(fig, "temperature_profiles")
    used("temperature_profiles", *files, note="Wall temperatures digitised from Liu et al. (2025), as in the manuscript; source region of the computations (manuscript Sec. 4.1).")


def main():
    FIG.mkdir(exist_ok=True)
    for f in (temperatures, graetz, damkohler, sensitivity, lobresco, liu, silica, thermodynamics):
        f()
    (HERE / "figure_sources.json").write_text(json.dumps(
        {"generated_by": "make_figures.py", "figures": SOURCES}, indent=1) + "\n")


if __name__ == "__main__":
    main()
