#!/usr/bin/env python3
"""Prepare Liu (2025) M10 references and mesh/input staging; standard library only.

No CFD solution or measured temperature profile is fabricated. All supplied
peak predictions are user inputs. Eq. (3) is a historical zero-dimensional
baseline with constant sublimation H/S, not the phase-resolved HKS model.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = 8.31446261815324
HE_MOLAR_MASS = 0.0040026
P_REFERENCE = 101325.0
REPRODUCTION_VOLUME_L_MOL = 24.465  # inferred; approximately 25 C, 1 atm
DOI = "10.1007/s10967-025-10246-4"
FORMULAS = {"PbI2": {"Pb": 1, "I": 2}, "BiI3": {"Bi": 1, "I": 3},
            "BiI": {"Bi": 1, "I": 1}, "KI": {"K": 1, "I": 1}}


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, columns):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def positive(value, name):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(name + " must be finite and positive")
    return number


def eq3(species, mass_kg, helium_mol):
    """K; p0 = 1 bar and p_species/p0 approximated by n_species/n_He."""
    ratio = mass_kg / float(species["molar_mass_kg_mol"]) / helium_mol
    denominator = float(species["delta_S_J_mol_K"]) - R * math.log(ratio)
    return float(species["delta_H_J_mol"]) / positive(denominator, "Eq. 3 denominator")


def load_references():
    experiments = read_csv(HERE / "experiments.csv")
    targets = read_csv(HERE / "targets.csv")
    thermo = read_csv(HERE / "sublimation.csv")
    by_experiment = {row["experiment"]: row for row in experiments}
    by_species = {row["species"]: row for row in thermo}
    if len(by_experiment) != len(experiments) or len(by_species) != len(thermo):
        raise ValueError("duplicate experiment or thermodynamic species")
    keys = [(row["experiment"], row["species"]) for row in targets]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate target")
    for row in experiments:
        for field in ("sample_mass_mg", "inner_diameter_m", "column_length_m",
                      "flow_ml_min", "duration_s"):
            positive(row[field], field)
    for row in thermo:
        for field in ("molar_mass_kg_mol", "iodine_atoms", "delta_H_J_mol",
                      "delta_S_J_mol_K"):
            positive(row[field], field)
    for row in targets:
        if row["experiment"] not in by_experiment or row["species"] not in by_species:
            raise ValueError("unknown experiment/species in target")
        for field in ("source_species_mass_mg", "uncertainty_K"):
            positive(row[field], field)
        positive(float(row["measured_peak_C"]) + 273.15, "measured peak K")
    for experiment in experiments:
        rows = [row for row in targets if row["experiment"] == experiment["experiment"]]
        if experiment["total_deposited_iodine_mg"]:
            if abs(sum(float(row["iodine_yield_pct"]) for row in rows) - 100) > 1e-9:
                raise ValueError("published iodine yields must include all species")
            for row in rows:
                species = by_species[row["species"]]
                derived = (float(row["source_species_mass_mg"])
                           / float(species["molar_mass_kg_mol"])
                           * int(species["iodine_atoms"]) * 0.12690447)
                stated = (float(experiment["total_deposited_iodine_mg"])
                          * float(row["iodine_yield_pct"]) / 100)
                if abs(derived / stated - 1) > 0.02:
                    raise ValueError("species mass and iodine yield disagree beyond rounding")
    return by_experiment, targets, by_species


def baseline(experiments, targets, species):
    rows = []
    for target in targets:
        experiment = experiments[target["experiment"]]
        volume_l = float(experiment["flow_ml_min"]) * float(experiment["duration_s"]) / 60000
        mass_kg = float(target["source_species_mass_mg"]) * 1e-6
        temperature = eq3(species[target["species"]], mass_kg,
                          volume_l / REPRODUCTION_VOLUME_L_MOL) - 273.15
        cold_volume = R * 273.15 / P_REFERENCE * 1000
        cold_temperature = eq3(species[target["species"]], mass_kg,
                               volume_l / cold_volume) - 273.15
        difference = temperature - float(target["published_eq3_C"])
        rows.append(dict(experiment=target["experiment"], species=target["species"],
                         eq3_reproduction_C=temperature,
                         published_eq3_C=float(target["published_eq3_C"]),
                         reproduction_difference_K=difference,
                         eq3_273p15K_flow_reference_C=cold_temperature,
                         measured_peak_C=float(target["measured_peak_C"]),
                         uncertainty_K=float(target["uncertainty_K"]),
                         eq3_minus_measured_peak_K=temperature - float(target["measured_peak_C"])))
    return rows


def mesh_dictionary(length, radius, nx, nr, angle):
    rt = radius * math.tan(math.radians(angle / 2))
    return f"""// Preliminary axisymmetric mesh; not a mesh-convergence result.
FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
mergeType points;
scale 1;
vertices (
  (0 0 0) ({length} 0 0) ({length} {radius} {-rt}) (0 {radius} {-rt})
  ({length} {radius} {rt}) (0 {radius} {rt})
);
blocks (hex (0 1 2 3 0 1 4 5) ({nx} {nr} 1) simpleGrading (1 1 1));
boundary (
  INLET {{ type patch; faces ((0 0 5 3)); }}
  OUTLET {{ type patch; faces ((1 2 4 1)); }}
  WALL {{ type wall; faces ((3 5 4 2)); }}
  BACK {{ type wedge; faces ((0 3 2 1)); }}
  FRONT {{ type wedge; faces ((0 1 4 5)); }}
  AXIS {{ type empty; faces ((0 1 1 0)); }}
);
"""


def compare(path, experiments, targets):
    predictions = read_csv(path) if path else []
    expected = {(row["experiment"], row["species"]): row for row in targets
                if experiments[row["experiment"]]["primary_m10"] == "yes"
                and row["scope"] == "M10"}
    supplied = {}
    for row in predictions:
        key = (row["experiment"], row["species"])
        if key not in expected or key in supplied:
            raise ValueError("unknown or duplicate primary prediction: " + str(key))
        supplied[key] = positive(row["peak_temperature_K"], "predicted peak K")
    rows = []
    for key, target in expected.items():
        row = dict(experiment=key[0], species=key[1], status="NOTRUN",
                   predicted_peak_K="", measured_peak_K=float(target["measured_peak_C"]) + 273.15,
                   difference_K="", uncertainty_K=float(target["uncertainty_K"]),
                   within_experimental_uncertainty="", science_target_K=20 if key[1] == "PbI2" else 40,
                   within_science_target="")
        if key in supplied:
            difference = supplied[key] - row["measured_peak_K"]
            row.update(status="REPORTED", predicted_peak_K=supplied[key], difference_K=difference,
                       within_experimental_uncertainty=abs(difference) <= row["uncertainty_K"],
                       within_science_target=abs(difference) <= row["science_target_K"])
        rows.append(row)
    return rows


def prepare(out, experiments, targets, species, args):
    comparisons = compare(args.predictions, experiments, targets)
    out.mkdir(parents=True, exist_ok=True)
    checks = baseline(experiments, targets, species)
    write_csv(out / "eq3_baseline.csv", checks, list(checks[0]))
    selected = [row for row in experiments.values() if row["primary_m10"] == "yes"]
    fraction = args.wedge_angle / 360
    for experiment in selected:
        name = experiment["experiment"]
        directory = out / name
        (directory / "system").mkdir(parents=True, exist_ok=True)
        (directory / "system" / "blockMeshDict").write_text(mesh_dictionary(
            float(experiment["column_length_m"]), float(experiment["inner_diameter_m"]) / 2,
            args.nx, args.nr, args.wedge_angle))
        (directory / "system" / "controlDict").write_text(
            "// Meshing only; benchmark.json lists missing CFD inputs.\n"
            "FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\n"
            "application blockMesh; startFrom startTime; startTime 0;\n"
            "stopAt endTime; endTime 0; deltaT 1; writeControl timeStep; writeInterval 1;\n")
        (directory / "system" / "fvSchemes").write_text(
            "// Required by checkMesh; CFD transport schemes are not configured.\n"
            "FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }\n"
            "ddtSchemes { default Euler; } gradSchemes { default Gauss linear; }\n"
            "divSchemes { default none; } laplacianSchemes { default Gauss linear corrected; }\n"
            "interpolationSchemes { default linear; } snGradSchemes { default corrected; }\n")
        (directory / "system" / "fvSolution").write_text(
            "// Required by checkMesh; no CFD solvers are configured.\n"
            "FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }\n"
            "solvers {}\n")
        write_csv(directory / "wall_temperature.csv", [], ["x_m", "T_K"])
        sources = []
        excluded_yield = 0
        element_amounts = dict(Pb=0.0, Bi=0.0, I=0.0)
        for target in targets:
            if target["experiment"] != name:
                continue
            thermo = species[target["species"]]
            mass = float(target["source_species_mass_mg"]) * 1e-6
            mol = mass / float(thermo["molar_mass_kg_mol"])
            sources.append(dict(species=target["species"], scope=target["scope"],
                                formula=FORMULAS[target["species"]],
                                included_in_pbbii_model=target["scope"] == "M10",
                                mass_basis=target["mass_basis"], full_column_mass_kg=mass,
                                full_column_amount_mol=mol, mean_release_mol_s=mol / float(experiment["duration_s"]),
                                wedge_amount_mol=mol * fraction, wedge_mean_release_mol_s=mol * fraction / float(experiment["duration_s"]),
                                measured_peak_K=float(target["measured_peak_C"]) + 273.15,
                                uncertainty_K=float(target["uncertainty_K"]),
                                published_iodine_yield_pct=float(target["iodine_yield_pct"]) if target["iodine_yield_pct"] else None))
            if target["scope"] == "diagnostic":
                excluded_yield += float(target["iodine_yield_pct"])
            else:
                for element, count in FORMULAS[target["species"]].items():
                    element_amounts[element] += count * mol
        flow = float(experiment["flow_ml_min"]) * 1e-6 / 60
        temperature = args.flow_reference_temperature
        pressure = args.flow_reference_pressure
        manifest = dict(schema_version=1, source_doi=DOI, experiment=experiment,
                        ready_for_cfd=False, observable="wall deposit peak temperature; not 1% onset",
                        preliminary_source_model="Q1: mean prescribed release over 10800 s; measured evaporation history unknown",
                        physical_source_model="Q3 inventory for pure compounds; requires sample geometry and kinetic qualification" if experiment["sample_material"] != "LBE-I" else "Q2 element/metal-vapour source needed to predict species split",
                        flow_reference_temperature_K=temperature,
                        flow_reference_pressure_Pa=pressure,
                        full_column_volume_flow_at_reference_m3_s=flow,
                        full_column_helium_mass_flow_kg_s=flow * pressure * HE_MOLAR_MASS / (R * temperature) if temperature else None,
                        wedge_helium_mass_flow_kg_s=flow * pressure * HE_MOLAR_MASS / (R * temperature) * fraction if temperature else None,
                        mesh=dict(axial_cells=args.nx, radial_cells=args.nr, cells=args.nx * args.nr,
                                  wedge_angle_deg=args.wedge_angle, full_column_fraction=fraction),
                        source_terms=sources, excluded_KI_iodine_yield_pct=excluded_yield,
                        modeled_full_column_element_amounts_mol=element_amounts,
                        modeled_wedge_element_amounts_mol={key: amount * fraction for key, amount in element_amounts.items()},
                        missing_inputs=["measured wall T(x) and its coordinate origin", "boat axial position and source selection geometry",
                                        "physical frozen helium carrier", "raw deposition/gamma profiles", "evaporation history",
                                        "Bi-I data provenance confirmation and external tables (Bi cases)", "M10 implementation/qualification (Bi cases)"]
                        + ([] if temperature else ["confirm mass-flow-controller reference temperature and pressure"]))
        (directory / "benchmark.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    write_csv(out / "predictions-template.csv", [], ["experiment", "species", "peak_temperature_K"])
    write_csv(out / "comparison.csv", comparisons, list(comparisons[0]))
    hashes = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
              for name in ("experiments.csv", "targets.csv", "sublimation.csv", "prepare.py")}
    (out / "provenance.json").write_text(json.dumps(dict(source_doi=DOI, source_license="CC-BY-4.0",
        input_sha256=hashes, eq3_reproduction_molar_volume_L_mol=REPRODUCTION_VOLUME_L_MOL,
        eq3_flow_reference_note="24.465 L/mol inferred from rounded published Eq. 3 values; not a verified flow calibration",
        args=dict(nx=args.nx, nr=args.nr, wedge_angle=args.wedge_angle,
                  flow_reference_temperature_K=args.flow_reference_temperature,
                  flow_reference_pressure_Pa=args.flow_reference_pressure)), indent=2) + "\n")
    print(f"Prepared {len(selected)} M10 benchmark inputs; CFD status NOTRUN (missing inputs).")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="validate references and Eq. 3 reproduction")
    parser.add_argument("--out", type=Path, help="generate staging files into this directory")
    parser.add_argument("--predictions", type=Path, help="compare primary CFD peak predictions (CSV)")
    parser.add_argument("--nx", type=int, default=384)
    parser.add_argument("--nr", type=int, default=48)
    parser.add_argument("--wedge-angle", type=float, default=5)
    parser.add_argument("--flow-reference-temperature", type=float, help="explicit assumed/confirmed reference temperature in K")
    parser.add_argument("--flow-reference-pressure", type=float, help="explicit assumed/confirmed reference pressure in Pa")
    args = parser.parse_args()
    try:
        if not args.check and not args.out:
            raise ValueError("choose --check or --out")
        if args.predictions and not args.out:
            raise ValueError("--predictions requires --out")
        if args.nx < 1 or args.nr < 1 or not 0 < args.wedge_angle <= 10:
            raise ValueError("positive mesh sizes and a wedge angle in (0, 10] are required")
        if args.flow_reference_temperature is not None:
            positive(args.flow_reference_temperature, "flow reference temperature")
        if args.flow_reference_pressure is not None:
            positive(args.flow_reference_pressure, "flow reference pressure")
        if (args.flow_reference_temperature is None) != (args.flow_reference_pressure is None):
            raise ValueError("specify both flow reference temperature and pressure, or neither")
        experiments, targets, species = load_references()
        checks = baseline(experiments, targets, species)
        error = max(abs(row["reproduction_difference_K"]) for row in checks)
        if error > 0.55:
            raise ValueError(f"Eq. 3 reproduction error {error:.6g} K exceeds 0.55 K")
        print(f"References: {len(experiments)} experiments, {len(targets)} peaks; Eq. 3 max difference {error:.4f} K.")
        if args.out:
            prepare(args.out, experiments, targets, species, args)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()
