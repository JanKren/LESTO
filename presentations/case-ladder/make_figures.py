#!/usr/bin/env python3
"""Draw the figures of the case-ladder deck into figures/.

Every figure is drawn from a tracked repository file or from the small
outputs of the runs made for this deck (runs/, written by extract_runs.py
from a scratch mirror of run/; commands, CPUs and wall times in
runs/runs.json).  figure_sources.json lists, for each figure, the files it
reads and their sha256.  No number is typed in here except the
dictionary values the cases themselves state (inlet profile, wall
temperature ramps), each marked where it is used.

    MPLCONFIGDIR=/tmp/lesto-matplotlib python3 make_figures.py
"""
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RUN = REPO / "run"
RUNS = HERE / "runs"
FIG = HERE / "figures"

# PSI-derived categorical order (validated: CVD and normal-vision separation)
BLUE, ORANGE, MAGENTA, TEAL, PURPLE = "#2a3cf0", "#e08a00", "#dc005a", "#00a3a3", "#6e14dc"
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#d9d9d9"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.4, "legend.frameon": False, "savefig.dpi": 300,
    "figure.dpi": 150, "pdf.fonttype": 42,
})

SIZE = (3.0, 2.36)          # inches; the right column is about 6.8 cm wide
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def used(fig_name, *paths, note=""):
    entry = SOURCES.setdefault(fig_name, {"inputs": {}, "note": note})
    for p in paths:
        p = Path(p)
        entry["inputs"][str(p.relative_to(REPO))] = sha(p)
    if note:
        entry["note"] = note


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight", pad_inches=0.02, dpi=200)
    plt.close(fig)


def residuals(log):
    """First initial residual of every field in every SIMPLE iteration."""
    hist = {}
    blocks = re.split(r"^Time = ", Path(log).read_text(), flags=re.M)[1:]
    for b in blocks:
        seen = set()
        for f, r in re.findall(r"Solving for (\w+), Initial residual = (\S+),", b):
            if f not in seen:
                hist.setdefault(f, []).append(float(r))
                seen.add(f)
    return {k: np.array(v) for k, v in hist.items()}


def converged_in(log):
    m = re.search(r"converged in (\d+) iterations", Path(log).read_text())
    return int(m.group(1)) if m else None


def profile(path):
    """An axialProfiles .dat: header dict and the columns."""
    head, rows = {}, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# columns"):
            head["columns"] = line.split()[2:]
        elif line.startswith("# Tdep fromInlet"):
            t = line.split()
            head["Tdep"] = float(t[t.index("T") + 1]) if "found" in t else None
            head["xdep"] = float(t[t.index("x") + 1]) if "found" in t else None
        elif line.startswith("# Tpeak"):
            head["Tpeak"] = float(line.split()[2])
        elif line.startswith("# inventory"):
            t = line.split()
            head["inventory"] = {t[i]: float(t[i + 1]) for i in range(2, len(t) - 1, 2)}
        elif line and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    return head, {c: a[:, i] for i, c in enumerate(head["columns"])}


def read_csv(path):
    with open(path) as f:
        r = list(csv.DictReader(f))
    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan
    return {k: np.array([num(row[k]) for row in r]) for k in r[0]}


def label_end(ax, x, y, text, color, dx=3, dy=0):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                color=INK, fontsize=7.5, va="center")


def balance(path):
    names = None
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# time"):
            names = line[2:].split()
        elif line and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    return {n: a[:, i] for i, n in enumerate(names)}


# --------------------------------------------------------------------------
# 001-007: carrier flow (Bojan Niceno)
# --------------------------------------------------------------------------
def case(n):
    return next(RUN.glob(f"{n}-*"))


def fig001():
    e = np.load(RUNS / "001" / "inlet_edges.npz")
    p, lines = e["points"] * 1e3, e["lines"]
    log = case("001") / "out"
    h = residuals(log)
    fig, (a, b) = plt.subplots(1, 2, figsize=SIZE, gridspec_kw={"width_ratios": [1, 1.25]})
    seg = np.stack([p[lines[:, 0]], p[lines[:, 1]]], axis=1)
    from matplotlib.collections import LineCollection
    a.add_collection(LineCollection(seg, colors=INK, linewidths=0.25))
    a.set_aspect("equal"); a.set_xlim(-2.5, 2.5); a.set_ylim(-2.5, 2.5)
    a.set_xlabel("y [mm]"); a.set_ylabel("z [mm]"); a.grid(False)
    a.set_title(f"INLET patch, {int(e['n_faces'])} faces", color=INK)
    for f, c in (("Ux", BLUE), ("p", MAGENTA)):
        b.semilogy(np.arange(1, len(h[f]) + 1), h[f], color=c, label=f)
    b.set_xlabel("SIMPLE iteration"); b.set_ylabel("initial residual")
    b.set_title(f"converged in {converged_in(log)}", color=INK)
    b.legend(loc="upper right")
    fig.tight_layout(w_pad=0.6)
    save(fig, "c001")
    used("c001", RUNS / "001" / "inlet_edges.npz", log,
         note="Cross-section: INLET patch of the gmshToFoam conversion of run/001/pipe.msh (foamToVTK); residuals: run/001/out.")


def fig002():
    l1, l2 = case("001") / "out", case("002") / "out"
    h1, h2 = residuals(l1), residuals(l2)
    R, Umax, Ub = 2.4, 0.19, 0.095            # mm, m/s: 0/U of cases 001 and 002
    r = np.linspace(-R, R, 101)
    fig, (a, b) = plt.subplots(1, 2, figsize=SIZE, gridspec_kw={"width_ratios": [1, 1.25]})
    a.plot(np.full_like(r, Ub), r, color=MUTED, ls="--", label="001 uniform")
    a.plot(Umax * (1 - (r / R) ** 2), r, color=BLUE, label="002 parabolic")
    a.set_xlabel("inlet $U_x$ [m/s]"); a.set_ylabel("r [mm]")
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1)
    for h, lg, c, n in ((h1, l1, MUTED, "001"), (h2, l2, BLUE, "002")):
        b.semilogy(np.arange(1, len(h["p"]) + 1), h["p"], color=c, label=f"{n}: {converged_in(lg)} it.")
    b.set_xlabel("SIMPLE iteration"); b.set_ylabel("p initial residual")
    b.legend(loc="upper right")
    fig.tight_layout(w_pad=0.6)
    save(fig, "c002")
    used("c002", l1, l2, case("002") / "0" / "U", note="Inlet profiles from 0/U (uniform 0.095 m/s; 0.19(1-r^2/R^2)); residuals from the out logs.")


def fig003():
    l2, l3 = case("002") / "out", case("003") / "out"
    h2, h3 = residuals(l2), residuals(l3)
    fig, ax = plt.subplots(figsize=SIZE)
    for h, lg, c, n, ls in ((h2, l2, MUTED, "002 uniform initial field", "--"), (h3, l3, BLUE, "003 parabolic initial field", "-")):
        ax.semilogy(np.arange(1, len(h["Ux"]) + 1), h["Ux"], color=c, ls=ls,
                    label=f"{n}: {converged_in(lg)} it.")
    ax.set_xlabel("SIMPLE iteration"); ax.set_ylabel("$U_x$ initial residual")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "c003")
    used("c003", l2, l3, note="Ux residual histories of run/002/out and run/003/out.")


def fig004():
    lg = case("004") / "out"
    h = residuals(lg)
    fig, ax = plt.subplots(figsize=SIZE)
    fvs = (case("004") / "system" / "fvSolution").read_text()
    rc = dict(re.findall(r"^\s*(p|U|e)\s+(\S+);", fvs.split("residualControl")[1].split("}")[0], flags=re.M))
    for f, c, key in (("Ux", BLUE, "U"), ("p", MAGENTA, "p"), ("e", ORANGE, "e")):
        it = np.arange(1, len(h[f]) + 1)
        tol = float(rc[key])
        ax.semilogy(it, h[f], color=c, label=f"{f} (criterion {tol:g})")
        ax.axhline(tol, color=c, lw=0.7, ls=":")
    ax.set_xlabel("SIMPLE iteration"); ax.set_ylabel("initial residual")
    ax.set_title(f"rhoSimpleFoam, converged in {converged_in(lg)}", color=INK)
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "c004")
    used("c004", lg, case("004") / "system" / "fvSolution", note="Residuals of run/004/out; criteria from system/fvSolution (residualControl).")


def fig005():
    l4, l5 = case("004") / "out", case("005") / "out"
    h4, h5 = residuals(l4), residuals(l5)
    x = np.linspace(0, 0.69, 50)
    fig, (a, b) = plt.subplots(1, 2, figsize=SIZE)
    a.plot(x, np.full_like(x, 300.0), color=MUTED, ls="--", label="004: 300 K")
    a.plot(x, 400 - 100 * x / 0.69, color=MAGENTA, label="005: linear")   # 0.orig/T WALL code
    a.set_xlabel("x [m]"); a.set_ylabel("wall T [K]"); a.set_ylim(280, 420)
    a.legend(loc="upper right")
    for h, c, n, ls in ((h4, MUTED, "004", "--"), (h5, MAGENTA, "005", "-")):
        b.semilogy(np.arange(1, len(h["e"]) + 1), h["e"], color=c, ls=ls, label=n)
    b.set_xlabel("SIMPLE iteration"); b.set_ylabel("e initial residual")
    b.legend(loc="upper right")
    fig.tight_layout(w_pad=0.6)
    save(fig, "c005")
    used("c005", l4, l5, case("005") / "0.orig" / "T", note="Wall ramp 400-100x/0.69 K from 0.orig/T (codedFixedValue); e residuals of the out logs.")


def fig006():
    l5, l6 = case("005") / "out", case("006") / "out"
    h5, h6 = residuals(l5), residuals(l6)
    fig, ax = plt.subplots(figsize=SIZE)
    for f, c in (("Ux", BLUE), ("p", MAGENTA), ("e", ORANGE)):
        ax.semilogy(np.arange(1, len(h6[f]) + 1), h6[f], color=c, label=f"006 {f}")
    ax.semilogy(np.arange(1, len(h5["p"]) + 1), h5["p"], color=MUTED, ls="--", label="005 p (rhoConst)")
    ax.set_xlabel("SIMPLE iteration"); ax.set_ylabel("initial residual")
    ax.set_title(f"perfectGas: {converged_in(l6)} it. (005: {converged_in(l5)})", color=INK)
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout()
    save(fig, "c006")
    used("c006", l5, l6, note="Residual histories of run/005/out and run/006/out.")


def fig007():
    s = read_csv(RUNS / "007" / "carrier261_centre_slice.csv")
    ax_ = read_csv(RUNS / "007" / "carrier261_axis.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=SIZE, sharex=True, gridspec_kw={"height_ratios": [1, 1.3]})
    # bin the slice cells on a regular grid (2.5 mm x 0.12 mm) for display
    xe = np.linspace(0, 0.69, 277); ye = np.linspace(-2.4e-3, 2.4e-3, 13)
    num, _, _ = np.histogram2d(s["x_m"], s["r_signed_m"], bins=(xe, ye), weights=s["T_K"])
    cnt, _, _ = np.histogram2d(s["x_m"], s["r_signed_m"], bins=(xe, ye))
    grid = np.where(cnt > 0, num / np.maximum(cnt, 1), np.nan)
    pm = a.pcolormesh(xe, ye * 1e3, grid.T, cmap="inferno", vmin=300, vmax=1000, shading="flat", rasterized=True)
    a.set_ylabel("r [mm]"); a.grid(False)
    cax = a.inset_axes([0.62, 1.04, 0.36, 0.08])
    cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
    cb.set_ticks([300, 650, 1000]); cb.ax.tick_params(labelsize=6.5)
    a.set_title("T [K], centre plane (radius stretched)", color=INK, loc="left", fontsize=7.5)
    b.plot(ax_["x_m"], ax_["T_axis_K"], color=MAGENTA, label="axis (r < 0.3 mm)")
    b.plot(ax_["x_m"], 1000 - 700 * ax_["x_m"] / 0.69, color=MUTED, ls="--", label="wall, 0.orig/T")
    b.set_xlim(0, 0.69); b.set_xlabel("x [m]"); b.set_ylabel("T [K]")
    b.legend(loc="upper right")
    fig.tight_layout(h_pad=0.4)
    save(fig, "c007")
    used("c007", RUNS / "007" / "carrier261_centre_slice.csv", RUNS / "007" / "carrier261_axis.csv",
         case("007") / "261" / "T", case("007") / "0.orig" / "T",
         note="foamToVTK of run/007/261 on the converted pipe mesh; wedge |z|<0.27|y|+0.12 mm at signed radius, binned 2.5 mm x 0.12 mm; axis cells r<0.3 mm per axial cell layer.")


# --------------------------------------------------------------------------
# 008-011: species transport (Bojan Niceno), runs to t = 1 s made for the deck
# --------------------------------------------------------------------------
def p0811(n):
    return read_csv(RUNS / n / "profiles_t1s.csv")


def fig008():
    p = p0811("008")
    x = p["x_m"]
    g, t = p["Y_PbI2_g_kg_per_m"] * 1e9, p["Y_tracer_kg_per_m"] * 1e9
    fig, ax = plt.subplots(figsize=SIZE)
    ax.plot(x, g, color=BLUE, lw=2.4, label="PbI2_g")
    ax.plot(x, t, color=ORANGE, lw=1.0, ls="--", label="tracer")
    d = np.max(np.abs(g - t)) / np.max(g)
    ax.set_xlabel("x [m]"); ax.set_ylabel(r"$\int \rho Y\,dA$ [$\mu$g/m]")
    ax.set_title("t = 1 s: identical" if d == 0 else f"t = 1 s: max |PbI2_g - tracer| / max = {d:.1e}", color=INK)
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "c008")
    used("c008", RUNS / "008" / "profiles_t1s.csv", note="Run made for the deck: run/008 to t = 1 s on 3 ranks; cross-section integrals per axial cell layer.")
    return d


def fig009():
    p = p0811("009")
    x = p["x_m"]
    fig, (a, b) = plt.subplots(2, 1, figsize=SIZE, sharex=True, gridspec_kw={"height_ratios": [1.5, 1]})
    a.plot(x, p["Y_PbI2_g_kg_per_m"] * 1e9, color=BLUE, label="PbI2_g, D(T, p)")
    a.plot(x, p["Y_tracer_kg_per_m"] * 1e9, color=ORANGE, ls="--", label="tracer, D = 1e-5")
    a.set_ylabel(r"[$\mu$g/m]"); a.legend(loc="upper right")
    a.set_title("t = 1 s", color=INK)
    b.plot(x, p["D_PbI2_g_axis_m2_s"] * 1e5, color=BLUE)
    b.axhline(1.0, color=ORANGE, ls="--", lw=1.0)
    b.set_ylabel(r"D [$10^{-5}$ m$^2$/s]"); b.set_xlabel("x [m]")
    fig.tight_layout(h_pad=0.3)
    save(fig, "c009")
    used("c009", RUNS / "009" / "profiles_t1s.csv", note="Run made for the deck: run/009 to t = 1 s on 3 ranks; D_PbI2_g on the axis cells.")


def fig010():
    # The current solver carries the coupled mock of case 011 (commit 0be8d2b
    # replaced the source of 010), so a run made today is not case 010: the
    # figure uses the tracked out_excerpt of Bojan Nicenos 16 s run instead.
    f = case("010") / "out_excerpt"
    t, y = [], []
    tnow = None
    for line in f.read_text().splitlines():
        m = re.match(r"^Time = (\S+)$", line)
        if m:
            tnow = float(m.group(1))
        m = re.match(r"^Y_PbI2_s min/max = \S+ (\S+)", line)
        if m and tnow is not None:
            t.append(tnow); y.append(float(m.group(1)))
    t, y = np.array(t), np.array(y)
    fig, (a, b) = plt.subplots(2, 1, figsize=SIZE, gridspec_kw={"height_ratios": [1.3, 1]})
    a.plot([0, t[-1]], [0, y[-1] * 1e6], color=MUTED, ls="--", lw=0.8)
    a.plot(t, y * 1e6, "o", color=MAGENTA, ms=3, label="max $Y$ of PbI2_s")
    a.set_xlabel("t [s]"); a.set_ylabel(r"max $Y_{PbI2,s}$ [$10^{-6}$]")
    a.set_title(f"16 s run: grows linearly to {y[-1]:.3e}", color=INK)
    xx = np.linspace(0, 0.69, 100)
    Tw = 1000 - 700 * xx / 0.69                       # wall ramp of case 007
    b.plot(xx, (1000 - Tw) / (1000 - 300), color=MAGENTA)
    b.set_ylabel(r"$S/S_0$ at the wall"); b.set_xlabel("x [m]")
    fig.tight_layout(h_pad=0.4)
    save(fig, "c010")
    used("c010", f, case("007") / "0.orig" / "T",
         note="Tracked out_excerpt of run/010 (first 10 and last steps of the 16 s run, v2606); S/S0 = (Thot - T)/(Thot - Tcold) of the readme on the wall ramp of 007. Not a deck run: the current solver implements the coupled mock of 011.")


def fig011():
    p = p0811("011")
    x = p["x_m"]
    fig, (a, b) = plt.subplots(2, 1, figsize=SIZE, sharex=True)
    a.plot(x, p["Y_PbI2_g_kg_per_m"] * 1e9, color=BLUE)
    a.set_ylabel(r"gas [$\mu$g/m]"); a.set_title("t = 1 s", color=INK)
    b.plot(x, p["Y_PbI2_s_kg_per_m"] * 1e12, color=MAGENTA)
    b.set_ylabel("solid [ng/m]"); b.set_xlabel("x [m]")
    fig.tight_layout(h_pad=0.3)
    save(fig, "c011")
    used("c011", RUNS / "011" / "profiles_t1s.csv", note="Run made for the deck: run/011 to t = 1 s on 3 ranks.")


# --------------------------------------------------------------------------
# 012-021: gas/condensate exchange (Jan Kren)
# --------------------------------------------------------------------------
MPbI2 = 0.46100894          # kg/mol, constant/speciesTransportProperties


def fig012():
    f = RUNS / "012" / "steps100" / "PbI2_g_balance.dat"
    b = balance(f)
    t = b["time"] * 1e3
    fig, (a, c) = plt.subplots(2, 1, figsize=SIZE, sharex=True, gridspec_kw={"height_ratios": [1.4, 1]})
    a.plot(t, b["released"] / MPbI2 * 1e9, color=MUTED, lw=2.4, label="released")
    a.plot(t, b["gas"] / MPbI2 * 1e9, color=BLUE, ls="--", label="gas")
    a.set_ylabel("[nmol]"); a.legend(loc="upper left")
    a.set_title("first 100 steps of the 6 s release", color=INK)
    rel = np.abs(b["closure"]) / np.maximum(b["released"], 1e-300)
    nz = rel > 0
    c.semilogy(t[nz], rel[nz], ls="none", color=MAGENTA, marker="o", ms=3)
    c.set_ylabel("|closure|/released"); c.set_xlabel("t [ms]")
    c.set_ylim(1e-18, 1e-14)
    c.text(0.98, 0.9, f"{int(nz.sum())} of {len(rel)} steps non-zero, max {rel.max():.1e}",
           transform=c.transAxes, ha="right", va="top", fontsize=7.5, color=INK)
    fig.tight_layout(h_pad=0.3)
    save(fig, "c012")
    used("c012", f, note="Run made for the deck: ./Allrun steps 100 of run/012 (serial). Zero closures drawn at 1e-20.")
    return float(np.max(rel))


def fig013():
    tab = REPO / "thermochemistry" / "systems" / "data" / "pv_PbI2_phases.csv"
    lines = tab.read_text().splitlines()
    hdr = [l for l in lines if l.startswith("# T_K")][0][2:].split(",")
    a_ = np.array([[float(v) for v in l.split(",")] for l in lines if l and not l.startswith("#")])
    f = RUNS / "013" / "steps100" / "PbI2_g_balance.dat"
    b = balance(f)
    fig, (a, c) = plt.subplots(1, 2, figsize=SIZE)
    T = a_[:, 0]
    for i, (name, col) in enumerate(zip(hdr[1:3], (BLUE, MAGENTA))):
        a.plot(T, a_[:, 1 + i], color=col, label=name)
    a.set_xlim(400, 1000); a.set_ylim(-6, 5)
    a.axvline(683, color=MUTED, lw=0.8, ls="--")      # T_m of the table header
    a.text(695, -5.3, "683 K", fontsize=7.5, color=INK, ha="left")
    a.legend(loc="upper left", fontsize=7, handlelength=1.2)
    a.set_xlabel("T [K]"); a.set_ylabel(r"log$_{10}$ $p_{eq}$ [Pa]")
    t = b["time"] * 1e3
    c.plot(t, b["sample"] / MPbI2 * 1e6, color=ORANGE, label="boat")
    c.plot(t, b["gas"] / MPbI2 * 1e6, color=BLUE, label="gas")
    c.set_xlabel("t [ms]"); c.set_ylabel(r"[$\mu$mol]")
    c.legend(loc="center right")
    fig.tight_layout(w_pad=0.6)
    save(fig, "c013")
    used("c013", tab, f, note="p_eq table of the case (link constant/pv_PbI2_phases.csv); run made for the deck: ./Allrun steps 100 of run/013.")
    return hdr


def deposit_vs_T(ax, p, color, label, scale=1e12):
    ax.plot(p["Twall"], p["wall"] * scale, color=color, label=label, drawstyle="steps-mid")


def fig014():
    f = RUNS / "014" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    h, p = profile(f)
    fig, ax = plt.subplots(figsize=SIZE)
    deposit_vs_T(ax, p, BLUE, "wall deposit at 3 s")
    ax.axvline(h["Tdep"], color=MAGENTA, lw=1.0)
    ax.axvline(567.057, color=MUTED, lw=1.0, ls="--")       # readme: plug-flow saturation crossing
    ymax = ax.get_ylim()[1]
    ax.annotate(f"1 % onset\n{h['Tdep']:.2f} K", (h["Tdep"], ymax * 0.45), xytext=(605, ymax * 0.45),
                color=INK, fontsize=7.5, va="center", arrowprops=dict(arrowstyle="->", color=INK, lw=0.7))
    ax.text(567.057 + 2, ymax * 0.97, "plug flow\n567.06 K", color=MUTED, fontsize=7.5, ha="right", va="top")
    ax.set_xlim(620, 480)
    ax.set_xlabel(r"wall temperature [K], flow $\rightarrow$")
    ax.set_ylabel("PbI2 deposit [pmol/m]")
    fig.tight_layout()
    save(fig, "c014")
    used("c014", f, note="Run made for the deck: ./Allrun of run/014 (PASS 9/9); plug-flow crossing 567.057 K from the readme.")


def fig016():
    f4 = RUNS / "014" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    f6 = RUNS / "016" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    h4, p4 = profile(f4)
    h6, p6 = profile(f6)
    fig, ax = plt.subplots(figsize=SIZE)
    ax.plot(p4["x"] * 1e3, p4["wall"] * 1e12, color=MUTED, ls="--", drawstyle="steps-mid", label="014 cell-centre p")
    ax.plot(p6["x"] * 1e3, p6["wall"] * 1e12, color=BLUE, drawstyle="steps-mid", label="016 halfCell")
    l1 = np.sum(np.abs(p6["wall"] - p4["wall"])) / np.sum(p4["wall"])
    ax.set_xlim(55, 95)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("PbI2 deposit at 3 s [pmol/m]")
    ax.set_title(f"deposit L1(016 vs 014) = {100 * l1:.1f} %", color=INK)
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "c016")
    used("c016", f4, f6, note="Runs made for the deck: ./Allrun of run/014 and run/016 (both PASS).")
    return l1


def fig018():
    fP = RUNS / "018" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    fB = RUNS / "018" / "axialProfiles" / "3" / "BiI3_g_native.dat"
    hP, pP = profile(fP)
    hB, pB = profile(fB)
    fig, ax = plt.subplots(figsize=SIZE)
    deposit_vs_T(ax, pP, BLUE, "PbI2")
    deposit_vs_T(ax, pB, ORANGE, "BiI3")
    for h, crossing in ((hP, 567.057), (hB, 450.383)):        # readme: plug-flow crossings
        ax.axvline(h["Tdep"], color=INK, lw=0.8)
        ax.axvline(crossing, color=MUTED, lw=0.8, ls="--")
    ymax = ax.get_ylim()[1]
    ax.text(hP["Tdep"] - 3, ymax * 0.95, f"PbI2 onset\n{hP['Tdep']:.1f} K", fontsize=7.5, va="top", ha="left", color=INK)
    ax.text(hB["Tdep"] + 3, ymax * 0.95, f"BiI3 onset\n{hB['Tdep']:.1f} K", fontsize=7.5, va="top", ha="right", color=INK)
    ax.set_xlim(620, 400)
    ax.set_xlabel("wall temperature [K]"); ax.set_ylabel("deposit at 3 s [pmol/m]")
    ax.legend(loc="center", bbox_to_anchor=(0.62, 0.55))
    fig.tight_layout()
    save(fig, "c018")
    used("c018", fP, fB, note="Run made for the deck: ./Allrun of run/018 with the m10 tables (pv_BiI3_phases.csv md5 0429d3fc...), PASS 8/8; dashed: plug-flow crossings of the readme.")


def fig019():
    fP = RUNS / "019" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    fB = RUNS / "019" / "axialProfiles" / "3" / "BiI3_g_native.dat"
    f8P = RUNS / "018" / "axialProfiles" / "3" / "PbI2_g_native.dat"
    hP, pP = profile(fP)
    hB, pB = profile(fB)
    fig, ax = plt.subplots(figsize=SIZE)
    deposit_vs_T(ax, pP, BLUE, "PbI2")
    deposit_vs_T(ax, pB, ORANGE, "BiI3")
    for h in (hP, hB):
        ax.axvline(h["Tdep"], color=INK, lw=0.8)
    ymax = ax.get_ylim()[1]
    ax.text(hP["Tdep"] - 3, ymax * 0.95, f"PbI2 onset\n{hP['Tdep']:.1f} K", fontsize=7.5, va="top", color=INK)
    ax.text(hB["Tdep"] + 3, ymax * 0.95, f"BiI3 onset\n{hB['Tdep']:.1f} K", fontsize=7.5, va="top", ha="right", color=INK)
    ax.set_xlim(620, 400)
    ax.set_xlabel("wall temperature [K]"); ax.set_ylabel("deposit at 3 s [pmol/m]")
    ax.set_title("nine gases in equilibrium, inlet 1173 K", color=INK)
    ax.legend(loc="center", bbox_to_anchor=(0.62, 0.55))
    fig.tight_layout()
    save(fig, "c019")
    used("c019", fP, fB, note="Run made for the deck: ./Allrun of run/019 with the m10 tables (kernel, 1 ms).")


def fig015():
    f = case("015") / "out_excerpt"
    rows = re.findall(r"elements: GEMS (\d+), table (\d+) \(below minTemperature (\d+), below minMoleFraction (\d+)", f.read_text())
    a = np.array(rows, dtype=float)
    dl = [float(v) for v in re.findall(r"max \|dlog10 p_eq/p_eq,table\| (\S+)", f.read_text())]
    dl = dl[:len(a)]
    steps = np.arange(1, len(a) + 1)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=SIZE, gridspec_kw={"width_ratios": [1.4, 1]})
    ax.bar(steps, a[:, 0], color=BLUE, label="GEMS3K call")
    ax.set_xlabel("step"); ax.set_ylabel("elements per step")
    for s, v in zip(steps, a[:, 0]):
        ax.text(s, v + 10, f"{int(v)}", ha="center", fontsize=7, color=INK)
    ax.set_ylim(0, 820)
    ax.set_title(f"table: {int(a[0, 2])} below 500 K\n+ the cells without gas", color=INK)
    bx.semilogy(steps[1:], dl[1:], color=MAGENTA, marker="o", ms=3)
    bx.set_xlabel("step"); bx.set_ylabel(r"max |$\Delta$log$_{10}p_{eq}$|")
    bx.set_ylim(1e-6, 1e-2)
    bx.axhline(0.01, color=MUTED, ls="--", lw=0.8)
    bx.set_title("GEMS vs table", color=INK)
    fig.tight_layout(w_pad=0.6)
    save(fig, "c015")
    used("c015", f, note="Per-step GEMS lines of the tracked out_excerpt (5 steps; written with the bridge md5 fb67103e...). Dashed: the guard maxLog10Deviation 0.01.")


def fig020():
    f = case("020") / "qualification.json"
    q = json.loads(f.read_text())
    fig, ax = plt.subplots(figsize=SIZE)
    for sysname, col, mk in (("PbI2", BLUE, "o"), ("BiI3", ORANGE, "s"), ("LBEI", TEAL, "^")):
        dec, val, guard = [], [], []
        for k, v in q.items():
            if k.startswith(sysname + ":"):
                d = int(k.split(":")[1])
                if v.get("calls", 0) > 0:
                    dec.append(d); val.append(max(v["max_primal_dlog10"], 1e-6))
                else:
                    guard.append(d)
        o = np.argsort(dec)
        ax.semilogy(np.array(dec)[o], np.array(val)[o], color=col, marker=mk, ms=3.5, label=sysname)
    ax.axvspan(-3.5, -0.5, color=GRID, alpha=0.6, lw=0)
    ax.text(-2, 2e-6, "guarded:\nno call", ha="center", fontsize=7.5, color=INK)
    ax.set_xlabel(r"log$_{10}$ mole fraction of the iodide (decade)")
    ax.set_ylabel(r"max |$\Delta$log$_{10}$| of the raw IPM")
    ax.set_ylim(1e-6, 60)
    ax.legend(loc="upper left", ncol=3)
    fig.tight_layout()
    save(fig, "c020")
    used("c020", f, note="qualification.json of run/020: largest primal error of the raw GEMS3K result by mole-fraction decade (zero drawn at 1e-6); the exact kernel refines every local result.")


def fig021():
    base = RUNS / "021" / "axialProfiles" / "1.2"
    ver = json.loads((case("021") / "verification.json").read_text())["baseline"]
    files = (("condensate_Pb_s_native.dat", "Pb", BLUE), ("condensate_Bi_s_native.dat", "Bi", PURPLE),
             ("PbI2_g_native.dat", "PbI2", ORANGE), ("BiI3_g_native.dat", "BiI3", TEAL))
    fig, (a, b) = plt.subplots(2, 1, figsize=SIZE, sharex=True, gridspec_kw={"height_ratios": [1.6, 1]})
    T = None
    for name, lab, col in files:
        h, p = profile(base / name)
        a.plot(p["x"] * 1e3, p["wall"] / p["wall"].max(), color=col, label=lab, drawstyle="steps-mid")
        T = p
    a.set_ylabel("deposit / peak"); a.legend(loc="upper right", ncol=2)
    a.set_title("t = 1.2 s, plug-flow-limit carrier", color=INK)
    b.plot(T["x"] * 1e3, T["Twall"], color=MUTED)
    for k, v in ver["solver_onset_K"].items():
        xk = np.interp(-v, -T["Twall"], T["x"]) * 1e3
        b.plot(xk, v, "o", color=INK, ms=3)
        off = {"Pb(l)": (-22, -11), "Bi(l)": (4, 4), "Bi(cr)": (-28, -11), "PbI2(cr)": (4, 5), "BiI3(cr)": (4, 4)}.get(k, (3, 3))
        b.annotate(k, (xk, v), xytext=off, textcoords="offset points", fontsize=7, color=INK)
    b.set_ylabel("wall T [K]"); b.set_xlabel("x [mm]")
    fig.tight_layout(h_pad=0.3)
    save(fig, "c021")
    used("c021", *(base / n for n, _, _ in files), case("021") / "verification.json",
         note="Run made for the deck: ./Allrun of run/021 with the m10d baseline tables, PASS 8/8; onsets (dots) from verification.json.")


# --------------------------------------------------------------------------
# 022-024: benchmarks (Jan Kren)
# --------------------------------------------------------------------------
def fig022():
    f = case("022") / "constant" / "wall_x_profile_complete.dat"
    rows = [l.split() for l in f.read_text().splitlines() if l and not l.startswith("#")][1:]
    a = np.array([[float(r[0]), float(r[4])] for r in rows if len(r) >= 5])
    x = a[:, 0]          # the header says cm, the values are metres (680 K at 0.476 m, notes.md)
    T = a[:, 1]
    fig, ax = plt.subplots(figsize=SIZE)
    ax.plot(x, T, color=MAGENTA, marker=".", ms=3, lw=1.0, label="wall table (T-Flows)")
    # readme / notes: onsets of this case, of the paper and of the T-Flows parameter set
    for xv, lab, col, ls in ((0.4761, "022: 0.476 m", BLUE, "-"), (0.502, "paper: 0.502 m", INK, "--"),
                             (0.5011, "T-Flows set: 0.501 m", ORANGE, ":")):
        ax.axvline(xv, color=col, ls=ls, lw=1.1, label=lab)
    ax.axhline(680, color=MUTED, lw=0.7)
    ax.text(0.37, 683, "Tdep 680 K", fontsize=7.5, color=INK)
    ax.set_xlim(0.35, 0.70); ax.set_ylim(450, 950)
    ax.set_xlabel("x [m]"); ax.set_ylabel("wall T [K]")
    ax.set_title("1 % solid onset at 30 s", color=INK)
    ax.legend(loc="upper right", fontsize=7, handlelength=1.5)
    fig.tight_layout()
    save(fig, "c022")
    used("c022", f, case("022") / "readme.md", case("022") / "notes.md",
         note="Wall table constant/wall_x_profile_complete.dat (x in cm); onsets from the readme and notes (M4long runs fig9-tev15 and fig9-tev15-TFlows).")


LOB = REPO / "thermochemistry" / "benchmarks" / "lobresco2026" / "digitized"


def lobresco_panel(ax, level, name, legend=False):
    cdir = case("023") / "cases" / level / name
    h, p = profile(cdir / "final_profiles" / "PbI2_g_gamma.dat")
    comp = json.loads((cdir / "comparison.json").read_text())
    scale = comp["model_wallPlusGas_fromInlet"]["peak_mol_per_m"] / p["observable"].max()
    e = read_csv(LOB / f"{name}_experiment.csv")
    ax.plot(e["x_m"], e["n_mol_per_m"] * 1e6, color=INK, marker="o", ms=2.2, lw=0.8, label="measured")
    ax.plot(p["x"], p["observable"] * scale * 1e6, color=BLUE, drawstyle="steps-mid", label="this work")
    for kind, col, ls in (("temperature_based", ORANGE, "--"), ("vapour_pressure_based", TEAL, ":")):
        m = read_csv(LOB / f"{name}_model_{kind}.csv")
        ax.plot(m["x_m"], m["n_mol_per_m"] * 1e6, color=col, ls=ls, lw=1.0,
                label="paper, " + ("T-based" if kind.startswith("temp") else "$p_v$-based"))
    return cdir, scale, e


def fig017():
    fig, ax = plt.subplots(figsize=SIZE)
    cdir, scale, _ = lobresco_panel(ax, "coarse", "T106_6")
    comp = json.loads((cdir / "comparison.json").read_text())
    ax.set_xlim(0.25, 0.6)
    ax.set_xlabel("x [m]"); ax.set_ylabel(r"PbI2 [$\mu$mol/m]")
    ax.set_title(f"T106_6 at 6 s: $\\Delta T_{{dep}}$ = {comp['dT_dep_fromInlet_K']:+.1f} K", color=INK)
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "c017")
    used("c017", cdir / "final_profiles" / "PbI2_g_gamma.dat", cdir / "comparison.json",
         LOB / "T106_6_experiment.csv", LOB / "T106_6_model_temperature_based.csv",
         LOB / "T106_6_model_vapour_pressure_based.csv",
         note="Tracked record of case 023 coarse T106_6 (identical to the output of case 017's ./Allrun); wedge scaled to the full column by the factor of comparison.json.")
    return scale


def fig023():
    names = ("T106_6", "T106_60", "T106_300", "T5_60", "T540_60")
    fig, axes = plt.subplots(2, 3, figsize=(SIZE[0], SIZE[1] * 1.05))
    files = []
    for ax, n in zip(axes.flat, names):
        cdir, _, e = lobresco_panel(ax, "prod", n)
        ax.set_title(n, color=INK, fontsize=7.5, pad=2)
        ax.tick_params(labelsize=6.5)
        ax.set_xlim(0.0, 0.69)
        ax.set_yticklabels([])
        files += [cdir / "final_profiles" / "PbI2_g_gamma.dat", cdir / "comparison.json", LOB / f"{n}_experiment.csv"]
    leg = axes.flat[-1]
    leg.axis("off")
    hs, ls = axes.flat[0].get_legend_handles_labels()
    leg.legend(hs, ls, loc="center", fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("x [m]", fontsize=7.5)
    fig.tight_layout(pad=0.2, h_pad=0.5, w_pad=0.3)
    save(fig, "c023")
    used("c023", *files, note="Tracked records cases/prod/<case> of run/023 (production wedge, 318 x 72) against the digitised measurements and both published models; each panel on its own scale.")


def fig024():
    base = case("024") / "production_576x72_limitedLinear" / "profiles"
    fig, axes = plt.subplots(2, 1, figsize=SIZE, sharex=True)
    files = []
    for ax, n, title in ((axes[0], "PbI2_SS_Q1_baseline", "pure PbI2, steel"),
                         (axes[1], "LBE-I_SS_II_Q1_baseline", "LBE-I, steel (Q1)")):
        f = base / n / "iodine-comparison.csv"
        d = read_csv(f)
        ax.bar(d["x_m"] * 100, np.nan_to_num(d["figure_relative_yield_percent"]), width=0.9, color=GRID, edgecolor=MUTED, lw=0.3, label="Liu et al. (2025)")
        ax.plot(d["x_m"] * 100, d["predicted_relative_yield_percent"], color=BLUE, drawstyle="steps-mid", label="this work")
        ax.set_title(title, color=INK, fontsize=7.5, pad=2)
        ax.set_ylabel("I [% per cm]")
        files.append(f)
    axes[0].legend(loc="upper right")
    axes[1].set_xlim(25, 75); axes[1].set_xlabel("x [cm]")
    fig.tight_layout(h_pad=0.3)
    save(fig, "c024")
    used("c024", *files, note="Production runs (576 x 72, limitedLinear) of run/024: iodine wall profile against the digitised scan, as written by the analysis (iodine-comparison.csv).")


def main():
    FIG.mkdir(exist_ok=True)
    out = {}
    for f in (fig001, fig002, fig003, fig004, fig005, fig006, fig007,
              fig008, fig009, fig010, fig011, fig012, fig013, fig014, fig015,
              fig016, fig017, fig018, fig019, fig020, fig021, fig022, fig023, fig024):
        try:
            r = f()
        except FileNotFoundError as e:
            print(f"{f.__name__}: missing input {e.filename}")
            continue
        if r is not None:
            out[f.__name__] = r if isinstance(r, (int, float, str)) else str(r)
        print(f"{f.__name__}: ok")
    vals = []
    if "fig008" in out:
        d = out["fig008"]
        vals.append(r"\newcommand{\VAReightdiff}{max $|$PbI$_2-$tracer$|$/max = %s}" % (("%.1e" % d) if d > 0 else "0"))
    else:
        vals.append(r"\newcommand{\VAReightdiff}{(run pending)}")
    vals.append(r"\newcommand{\VARtwelveclosure}{$|$closure$|$/released %.1e}" % out["fig012"])
    runs = {r["label"]: r for r in json.loads((RUNS / "runs.json").read_text())}
    for n, w in (("008", "eight"), ("009", "nine"), ("010", "ten"), ("011", "eleven")):
        r = runs.get(f"{n}-1s")
        vals.append(r"\newcommand{\VAR%sT}{%s}" % (w, r["wall_s"] if r else "?"))
    (FIG / "values.tex").write_text("% written by make_figures.py\n" + "\n".join(vals) + "\n")
    (HERE / "figure_sources.json").write_text(json.dumps(
        {"generated_by": "make_figures.py", "runs": "runs/runs.json",
         "derived_values": out, "figures": SOURCES}, indent=1) + "\n")


if __name__ == "__main__":
    main()
