"""Figure 4(a): case-023 wall-only inventory against published deposits."""
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "ladder_figures", HERE.parent / "case-ladder" / "make_figures.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
plt = source.plt
plt.rcParams.update({"font.size": 11, "axes.labelsize": 11,
                     "axes.titlesize": 12, "xtick.labelsize": 10,
                     "ytick.labelsize": 10, "legend.fontsize": 11})
fig, axes = plt.subplots(3, 2, figsize=(7.0, 8.2))
for ax, name in zip(axes.flat, ("T106_6", "T106_60", "T106_300", "T5_60", "T540_60")):
    cdir, _, _ = source.lobresco_panel(ax, "prod", name)
    _, profile = source.profile(cdir / "final_profiles" / "PbI2_g_gamma.dat")
    # Convert the 5-degree wedge wall inventory to the full column and
    # mol/m to micromol/m. Residual gas is deliberately excluded.
    ax.lines[1].set_ydata(profile["wall"] * (360 / 5) * 1e6)
    ax.relim()
    ax.autoscale_view()
    flow, duration = name[1:].split("_")
    ax.set_title(f"{flow} std mL/min; {duration} s\n({name})", pad=9, fontsize=11)
    ax.set_xlim(0, 0.69)
    ax.set_xticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("Tube position, x [m]")
    ax.set_ylabel(r"Wall deposit [$\mu$mol m$^{-1}$]")
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(source.matplotlib.ticker.MaxNLocator(nbins=5))
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    for line in ax.lines:
        line.set_linewidth(1.6)
        if line.get_marker() == "o":
            line.set_markersize(3.5)
legend = axes.flat[-1]
legend.axis("off")
handles, _ = axes.flat[0].get_legend_handles_labels()
legend.legend(handles, ["Measured deposit", "This work: wall only", "Lobresco: temperature model",
                        "Lobresco: vapour-pressure model"],
              loc="center left", fontsize=10, handlelength=2.5, labelspacing=1.4)
fig.subplots_adjust(left=.12, right=.99, bottom=.07, top=.95, hspace=.55, wspace=.42)
out = HERE / "figures"
out.mkdir(exist_ok=True)
fig.savefig(out / "c023-report.pdf", bbox_inches="tight", pad_inches=.08)
fig.savefig(out / "c023-report.png", bbox_inches="tight", pad_inches=.08, dpi=200)
