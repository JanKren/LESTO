#!/usr/bin/env python3
"""Copy the small outputs of the runs made for this deck into runs/.

The runs were made on 2026-10-15 in a scratch mirror of run/ (copies of the
cases, links to applications/ and thermochemistry/), never inside run/; see
runs/runs.json for the commands, CPUs and wall times.  This script reads
that mirror and writes only below presentations/case-ladder/runs/, so that
make_figures.py can redraw every figure without the raw CFD output.

    python3 extract_runs.py <mirror>/run

Needs numpy and pyvista (VTK of cases 001, 007 and 008-011, written by
foamToVTK in the mirror).
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs"


def copy(src, dst):
    dst = OUT / dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)


def case(mirror, number):
    matches = sorted(mirror.glob(f"{number}-*"))
    assert len(matches) == 1, (number, matches)
    return matches[0]


def channel_outputs(mirror):
    """Profiles and ledgers of the channel cases 014, 016, 018, 019, 021."""
    for number, time in (("014", "3"), ("016", "3"), ("018", "3"),
                         ("019", "3"), ("021", "1.2")):
        work = case(mirror, number) / "work" / "case"
        for f in sorted((work / "postProcessing" / "axialProfiles" / time).glob("*.dat")):
            copy(f, f"{number}/axialProfiles/{time}/{f.name}")
        lines = [l for l in (work / "log.solver").read_text().splitlines()
                 if l.startswith(("Build", "Exec", "Date", "Pair ", "Interface",
                                  "areaPerVolume", "Phase change summary"))
                 or "G_e beta_c/rho_c" in l]
        (OUT / number / "log.solver.startup").write_text("\n".join(lines) + "\n")


def pipe_ledgers(mirror):
    """Ledgers of the 100-step runs of cases 012 and 013 (pipe, 007 carrier)."""
    for number in ("012", "013"):
        work = case(mirror, number) / "work" / "steps100"
        copy(work / "postProcessing" / "phaseChangeBalance" / "0" / "PbI2_g.dat",
             f"{number}/steps100/PbI2_g_balance.dat")
        lines = [l for l in (work / "log.solver").read_text().splitlines()
                 if l.startswith(("Build", "Exec", "Date", "Balance", "Sample",
                                  "Interface", "Temperature model", "Pair"))]
        (OUT / number / "steps100" / "log.summary").write_text("\n".join(lines) + "\n")


def axial_bins(mesh, fields, dx=0.0025, length=0.69):
    """Axial integrals (per metre) of rho*Y over cross-sections, and the
    wall temperature of the outermost cells, in bins of dx."""
    import pyvista as pv  # noqa: F401
    sized = mesh.compute_cell_sizes(length=False, area=False, volume=True)
    vol = np.abs(sized.cell_data["Volume"])
    cen = mesh.cell_centers().points
    x = cen[:, 0]
    r = np.hypot(cen[:, 1], cen[:, 2])
    # one bin per axial cell layer (the mesh is structured along x), so that
    # the per-metre integrals do not alias against the cell spacing
    layers = np.unique(np.round(x, 7))
    edges = np.concatenate(([0.0], 0.5 * (layers[1:] + layers[:-1]), [length]))
    nb = len(layers)
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, nb - 1)
    widths = np.diff(edges)
    out = {"x_m": layers}
    rho = mesh.cell_data["rho"]
    for name in fields:
        if name not in mesh.cell_data:
            continue
        m = rho * mesh.cell_data[name] * vol
        out[f"{name}_kg_per_m"] = np.bincount(idx, weights=m, minlength=nb) / widths
    # centreline temperature: cells within 0.3 mm of the axis
    core = r < 3e-4
    w = np.bincount(idx[core], minlength=nb)
    out["T_axis_K"] = np.bincount(idx[core], weights=mesh.cell_data["T"][core], minlength=nb) / np.maximum(w, 1)
    return out


def write_csv(path, cols, fmt="%.10g"):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(cols)
    data = np.column_stack([cols[k] for k in keys])
    np.savetxt(path, data, delimiter=",", header=",".join(keys), comments="", fmt=fmt)


def carrier_007(mirror):
    """A centre slice of the frozen carrier at 261 and
    the radial profile of U at three stations, from foamToVTK of zz-c007."""
    import pyvista as pv
    m = pv.read(str(mirror / "zz-c007" / "VTK" / "zz-c007_0" / "internal.vtu"))
    cen = m.cell_centers().points
    out = {"x_m": cen[:, 0], "T_K": m.cell_data["T"]}
    # a thin wedge around the y axis (|z| < 0.27|y| + 0.12 mm), drawn at the
    # signed radius, so that the O-grid cells near the wall are included
    sel = np.abs(cen[:, 2]) < 0.27 * np.abs(cen[:, 1]) + 1.2e-4
    out["r_signed_m"] = np.sign(cen[:, 1]) * np.hypot(cen[:, 1], cen[:, 2])
    write_csv(OUT / "007" / "carrier261_centre_slice.csv", {k: v[sel] for k, v in out.items()}, fmt="%.6g")
    bins = axial_bins(m, [])
    write_csv(OUT / "007" / "carrier261_axis.csv", bins)


def mesh_001(mirror):
    """Edges of the INLET patch of the pipe mesh (cross-section of run/001)."""
    import pyvista as pv
    b = pv.read(str(mirror / "zz-c007" / "VTK" / "zz-c007_0" / "boundary" / "INLET.vtp"))
    edges = b.extract_all_edges()
    pts = edges.points[:, 1:3]
    lines = edges.lines.reshape(-1, 3)[:, 1:]
    (OUT / "001").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "001" / "inlet_edges.npz", points=pts, lines=lines,
                        n_faces=b.n_cells)


def species_0811(mirror):
    """Axial profiles of the 1 s runs of cases 008-011 (foamToVTK of the
    reconstructed time 1)."""
    import pyvista as pv
    for number in ("008", "009", "010", "011"):
        d = case(mirror, number)
        vtm = sorted((d / "VTK").glob("*_*/internal.vtu"))
        if not vtm:
            print(f"{number}: no VTK yet")
            continue
        m = pv.read(str(vtm[-1]))
        carrier = pv.read(str(mirror / "zz-c007" / "VTK" / "zz-c007_0" / "internal.vtu"))
        for k in ("rho", "T"):
            if k not in m.cell_data:  # frozen fields are not written again
                m.cell_data[k] = carrier.cell_data[k]
        fields = [k for k in m.cell_data.keys() if k.startswith("Y_")]
        prof = axial_bins(m, fields)
        if "D_PbI2_g" in m.cell_data:
            # D_PbI2_g on the axis cells, by bin
            cen = m.cell_centers().points
            core = np.hypot(cen[:, 1], cen[:, 2]) < 3e-4
            edges = np.concatenate(([0.0], 0.5 * (prof["x_m"][1:] + prof["x_m"][:-1]), [0.69]))
            n = len(prof["x_m"])
            idx = np.clip(np.searchsorted(edges, cen[:, 0], side="right") - 1, 0, n - 1)
            w = np.bincount(idx[core], minlength=n)
            prof["D_PbI2_g_axis_m2_s"] = np.bincount(idx[core], weights=m.cell_data["D_PbI2_g"][core], minlength=n) / np.maximum(w, 1)
        write_csv(OUT / number / "profiles_t1s.csv", prof)
        log = mirror.parent.parent / f"log.{number}-1s"
        if log.exists():
            keep = [l for l in log.read_text().splitlines()
                    if l.startswith(("Build", "Exec", "Date", "D_PbI2_g", "Integrated source",
                                     "Net thermochemistry", "ExecutionTime", "Cells adjacent"))]
            (OUT / number / "log.summary").write_text("\n".join(keep[:12] + keep[-6:]) + "\n")


def main():
    mirror = Path(sys.argv[1]).resolve()
    OUT.mkdir(exist_ok=True)
    channel_outputs(mirror)
    pipe_ledgers(mirror)
    carrier_007(mirror)
    mesh_001(mirror)
    species_0811(mirror)
    runs = [json.loads(l) for l in (mirror.parent.parent / "runs.jsonl").read_text().splitlines()]
    for r in runs:
        r["dir"] = "<mirror>/run/" + r["dir"].split("/mirror/run/")[-1]
    (OUT / "runs.json").write_text(json.dumps(runs, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
