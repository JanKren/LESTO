#!/usr/bin/env python3
"""Create slide-sized figures and a dated status snapshot from saved evidence.

Only published digitisation and completed, tracked baseline analysis are used
for result plots. The live follow-up state supplies counts, never unfinished
scientific results. Run from any working directory.
"""
import collections
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ANALYSIS = REPO / "run/034-liu2025/full_analysis"
DIGITISED = REPO / "thermochemistry/benchmarks/liu2025/digitized"
WORK = REPO / "run/034-liu2025/work/followup-20261007"
FIGURES = HERE / "figures"
BLUE, RED, GREEN, GRAY = "#0014E6", "#DC005A", "#008F68", "#4B4B4B"
INPUTS = set()


def read_json(path):
    INPUTS.add(path)
    return json.loads(path.read_text())


def read_csv(path):
    INPUTS.add(path)
    with path.open() as stream:
        return list(csv.DictReader(stream))


def save(fig, name):
    fig.savefig(FIGURES / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(FIGURES / (name + ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)


def csv_xy(path, x="x_m", y="predicted_relative_yield_percent"):
    rows = read_csv(path)
    return np.array([float(r[x]) for r in rows]), np.array([float(r[y]) for r in rows])


def gamma(ax, experiment):
    x, y = csv_xy(DIGITISED / (experiment + "_gamma.csv"), y="relative_yield_percent_per_cm")
    ax.bar(x * 100, y, width=.8, color="#D4D4D4", label="Published scan")


def prediction(ax, experiment, suffix, label, colour, style="-"):
    x, y = csv_xy(ANALYSIS / "profiles" / (experiment + "_" + suffix) / "iodine-comparison.csv")
    ax.plot(x * 100, y, color=colour, lw=2.2, ls=style, label=label)


def main():
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 13,
        "axes.titlesize": 14, "axes.labelsize": 12,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": GRAY, "text.color": GRAY,
        "axes.labelcolor": GRAY, "xtick.color": GRAY, "ytick.color": GRAY,
        "legend.frameon": False, "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    analysis = read_json(ANALYSIS / "analysis.json")
    assert analysis["PASS"] == 15 and analysis["FAIL"] == 0

    # Original schematic, with geometry from the experimental summary.
    fig, ax = plt.subplots(figsize=(8, 2.8), layout="constrained")
    ax.set(xlim=(-.13, 1.36), ylim=(-.7, .7)); ax.axis("off")
    gradient = np.linspace(1, 0, 600).reshape(1, -1)
    ax.imshow(gradient, extent=(0, 1.25, -.12, .12), aspect="auto", cmap="coolwarm", alpha=.5)
    ax.add_patch(Rectangle((0, -.12), 1.25, .24, fill=False, ec=GRAY, lw=1.5))
    ax.add_patch(Rectangle((.09, -.12), .08, .075, color=GRAY))
    ax.text(.1, -.29, "Sample boat\nposition schematic", ha="center", va="top", fontsize=11)
    ax.text(.21, .28, "Hot source", color=RED, ha="center", weight="bold")
    ax.text(.87, .28, "Deposition along a\ntemperature gradient", color=BLUE, ha="center", weight="bold")
    for x in [.32, .61, .9]:
        ax.add_patch(FancyArrowPatch((x, 0), (x+.14, 0), arrowstyle="-|>", mutation_scale=18, color=GRAY))
    ax.add_patch(FancyArrowPatch((-.12, 0), (0, 0), arrowstyle="-|>", mutation_scale=18, color=BLUE))
    ax.text(-.12, .12, "He", color=BLUE, ha="left", weight="bold")
    ax.annotate("", xy=(0, -.54), xytext=(1.25, -.54), arrowprops={"arrowstyle":"<->", "color":GRAY})
    ax.text(.625, -.58, "Column length 1.25 m   |   inner diameter 4.8 / 5.0 mm", ha="center", va="top", fontsize=11)
    save(fig, "experiment")

    fig, ax = plt.subplots(figsize=(8, 3), layout="constrained")
    for name, label, colour in [("PbI2_SS", "Pure PbI$_2$ / steel", BLUE),
                                ("LBE-I_SS_II", "LBE--I / steel", GREEN),
                                ("LBE-I_SiO2_I", "LBE--I / silica", RED)]:
        x, y = csv_xy(DIGITISED/(name+"_temperature.csv"), y="T_K")
        ax.plot(x*100, y-273.15, lw=2.2, color=colour, label=label)
    ax.set(xlim=(0,125), ylim=(0,1000), xlabel="Published column coordinate (cm)", ylabel="Temperature ($^\\circ$C)")
    ax.legend(loc="upper right", fontsize=11); ax.grid(alpha=.15)
    save(fig, "temperatures")

    fig, ax = plt.subplots(figsize=(8, 3.15), layout="constrained")
    labels=[]
    for i, r in enumerate(analysis["q1_peaks"]):
        prefix={"PbI2_SS":"Pure PbI$_2$", "LBE-I_SS_II":"LBE / steel", "LBE-I_SiO2_I":"LBE / silica"}[r["experiment"]]
        labels.append(prefix+" / "+("PbI$_2$" if r["species"]=="PbI2" else "BiI$_3$"))
        u=float(r["uncertainty_K"]); error=float(r["error_K"])
        ax.barh(i, 2*u, left=-u, height=.56, color="#0014E6", alpha=.09)
        ax.plot(error, i, "o", ms=8, color=BLUE if r["species"]=="PbI2" else RED)
        ax.text(43, i, f"{error:+.1f} K", va="center", fontsize=12)
    ax.axvline(0, color=GRAY, lw=1)
    ax.set(yticks=range(5), yticklabels=labels, xlim=(-45,61), xlabel="Model minus published peak temperature (K)")
    ax.invert_yaxis(); ax.set_xticks([-40,-20,0,20,40]); ax.grid(axis="x", alpha=.15)
    save(fig, "peak_errors")

    fig, axes=plt.subplots(1,2,figsize=(8.4,3.4),layout="constrained")
    for ax, name, title, limits in zip(axes,["PbI2_SS","LBE-I_SS_II"],["Pure PbI$_2$ / steel","LBE--I / steel"],[(29,45),(43,65)]):
        gamma(ax,name)
        prediction(ax,name,"Q1_baseline","Q1 prescribed iodides",BLUE)
        if name.startswith("LBE"):
            prediction(ax,name,"Q2_baseline_sat0.001","Q2 / 0.1% saturation",GREEN)
        ax.set(xlim=limits, title=title, xlabel="Column coordinate (cm)")
        ax.grid(axis="y",alpha=.15)
    axes[0].set_ylabel("Iodine yield per 1 cm bin (%)")
    axes[1].legend(fontsize=9,loc="upper right")
    save(fig,"steel_profiles")

    fig,axes=plt.subplots(1,2,figsize=(8.4,3.4),layout="constrained")
    gamma(axes[0],"LBE-I_SiO2_I")
    prediction(axes[0],"LBE-I_SiO2_I","Q1_baseline","Q1 model",BLUE)
    axes[0].set(xlim=(45,70),title="Silica iodine profile",xlabel="Column coordinate (cm)",ylabel="Yield per 1 cm bin (%)")
    axes[0].legend(fontsize=10);axes[0].grid(axis="y",alpha=.15)
    x,y=csv_xy(DIGITISED/"LBE-I_SiO2_I_temperature.csv",y="T_K")
    axes[1].plot(100*x,y-273.15,color=RED,lw=2.2,label="Published curve")
    for r in analysis["q1_peaks"]:
        if r["experiment"]!="LBE-I_SiO2_I":continue
        px=float(r["figure_x_cm"]);py=float(r["published_peak_C"])
        axes[1].plot(px,py,"o",color=BLUE,ms=7)
        axes[1].annotate(("PbI$_2$" if r["species"]=="PbI2" else "BiI$_3$")+f"\n{py:.0f}$^\\circ$C",(px,py),xytext=(px-5,py-72),fontsize=10,
                         arrowprops={"arrowstyle":"-", "color":GRAY})
    axes[1].set(xlim=(45,70),ylim=(40,600),xlabel="Column coordinate (cm)",ylabel="Temperature ($^\\circ$C)",title="Printed peaks and curve conflict")
    axes[1].legend(fontsize=10); axes[1].grid(alpha=.15)
    save(fig,"silica_discrepancy")

    fig,axes=plt.subplots(1,2,figsize=(8.4,3.4),layout="constrained")
    for ax,exp,title in zip(axes,["LBE-I_SS_II","LBE-I_SiO2_I"],["LBE--I / steel","LBE--I / silica"]):
        rows=[r for r in analysis["thermodynamic_sensitivity"] if r["experiment"]==exp]
        x=np.arange(3);w=.34
        ax.bar(x-w/2,[r["baseline_BiI3_iodide_I_pct"] for r in rows],w,color=BLUE,label="Baseline BiI")
        ax.bar(x+w/2,[r["shift29k_BiI3_iodide_I_pct"] for r in rows],w,color=RED,label="BiI +29 kJ/mol")
        expected=next(r["observed_BiI3_iodide_I_pct"] for r in analysis["case_metrics"] if r["experiment"]==exp)
        ax.axhline(expected,color=GRAY,ls="--",lw=1.2)
        ax.text(2.45,expected+3,f"Paper: {expected:.1f}%",ha="right",fontsize=10,
                bbox={"facecolor":"white", "edgecolor":"none", "alpha":.9, "pad":1})
        ax.set(xticks=x,xticklabels=["0.1%","1%","10%"],ylim=(0,105),xlabel="Prescribed metal saturation",title=title)
        ax.grid(axis="y",alpha=.15)
    axes[0].set_ylabel("BiI$_3$ share of iodide iodine (%)")
    fig.legend(*axes[0].get_legend_handles_labels(),loc="outside lower center",ncol=2,fontsize=11)
    save(fig,"thermodynamic_sensitivity")

    examples=[("PbI2_SS_Q1_baseline","Pure PbI$_2$ / Q1"),
              ("LBE-I_SS_II_Q1_baseline","Steel / Q1"),
              ("LBE-I_SiO2_I_Q1_baseline","Silica / Q1"),
              ("LBE-I_SiO2_I_Q2_baseline_sat0.01","Silica / Q2 / 1%")]
    fig,ax=plt.subplots(figsize=(8,2.8),layout="constrained")
    y=np.arange(4);w=.33
    for index,species,colour in [(0,"PbI2",BLUE),(1,"BiI3",RED)]:
        values=[]
        for name,_ in examples:
            r=next((r for r in analysis["combined_refinement"] if r["case"]==name and r["species"]==species),None)
            values.append(100*r["profile_L1_relative_to_full"] if r else np.nan)
        ax.barh(y+(index-.5)*w,values,height=w,color=colour,label="PbI$_2$" if index==0 else "BiI$_3$")
    ax.set(yticks=y,yticklabels=[label for _,label in examples],xlim=(0,75),xlabel="Profile L$_1$ difference relative to full-size result (%)")
    ax.invert_yaxis();ax.legend(loc="lower right",fontsize=11);ax.grid(axis="x",alpha=.15)
    save(fig,"mixed_refinement")

    # Read one atomic live state; freeze it for a reproducible presentation.
    state_path=WORK/"run_state.json"
    if state_path.exists():
        state=read_json(state_path)
    else:
        state=read_json(REPO/"run/034-liu2025/followup_launch.json")
    counts=collections.Counter(r["status"] for r in state.get("cases",{}).values()) if isinstance(state.get("cases"),dict) else state["case_status_counts"]
    stamp=state.get("updated_utc",state.get("captured_utc"))
    local=dt.datetime.fromisoformat(stamp).astimezone(ZoneInfo("Europe/Zurich"))
    smoke=read_json(WORK/"smoke.json") if (WORK/"smoke.json").exists() else None
    snapshot={"captured_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"state_updated_utc":stamp,
              "local_snapshot_time":local.strftime("%d %B %Y, %H:%M %Z"),"counts":dict(counts),
              "case_count":64,"scope":"Status snapshot; no unfinished follow-up scientific result plotted.",
              "baseline_solver_source_commit":"fe214c6","study_commit":"361f02d",
              "baseline_summary":{"PASS":analysis["PASS"],"FAIL":analysis["FAIL"],"NOTRUN":analysis["NOTRUN"]}}
    if smoke:snapshot["native_smoke_PASS"]=sum(r["status"]=="PASS" for r in smoke["cases"])
    (HERE/"status_snapshot.json").write_text(json.dumps(snapshot,indent=2)+"\n")
    macros={"SnapshotTime":snapshot["local_snapshot_time"],"RunningCases":counts.get("RUNNING",0),
            "QueuedCases":counts.get("QUEUED",0),"CompletedCases":counts.get("PASS",0),"FailedCases":counts.get("FAIL",0)}
    (HERE/"status_snapshot.tex").write_text("% Generated by prepare_figures.py; frozen until explicitly regenerated.\n"+"".join(
        "\\newcommand{\\"+k+"}{"+str(v)+"}\n" for k,v in macros.items()))
    manifest={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUTS)}
    (HERE/"figure_sources.json").write_text(json.dumps({"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256":manifest,"note":"Digitised Liu curves retain published coordinates; no fitting or alignment applied."},indent=2)+"\n")
    print("Generated 7 vector figures, 7 previews and a frozen status snapshot:",snapshot["local_snapshot_time"],dict(counts))


if __name__=="__main__":
    main()
