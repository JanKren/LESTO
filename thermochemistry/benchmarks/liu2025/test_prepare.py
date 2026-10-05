"""Scientific-unit, source-budget and reporting checks for the preparation."""

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("liu_prepare", HERE / "prepare.py")
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


class Benchmarks(unittest.TestCase):
    def setUp(self):
        self.experiments, self.targets, self.species = PREPARE.load_references()

    def test_published_equation_and_iodine_budget(self):
        rows = PREPARE.baseline(self.experiments, self.targets, self.species)
        self.assertEqual(len(rows), 17)
        self.assertLess(max(abs(row["reproduction_difference_K"]) for row in rows), 0.55)
        # Published 100 mL/min for 3 h passes 18 L; 45 mL/min passes 8.1 L.
        # PbI2_SS at the inferred 24.465 L/mol yields about 641.86 K.
        pure = next(row for row in rows if row["experiment"] == "PbI2_SS")
        self.assertAlmostEqual(pure["eq3_reproduction_C"] + 273.15, 641.86, delta=0.02)
        rows = [row for row in self.targets if row["experiment"] == "LBE-I_SS_II"]
        self.assertEqual(sum(float(row["iodine_yield_pct"]) for row in rows), 100)
        self.assertAlmostEqual(sum(float(row["iodine_yield_pct"]) for row in rows
                                   if row["scope"] == "M10"), 93.4)

    def test_geometry_units_wedge_and_missing_inputs(self):
        with tempfile.TemporaryDirectory() as work:
            subprocess.run([sys.executable, str(HERE / "prepare.py"), "--out", work], check=True,
                           capture_output=True, text=True)
            steel = json.loads((Path(work) / "PbI2_SS/benchmark.json").read_text())
            silica = json.loads((Path(work) / "BiI3_SiO2_800/benchmark.json").read_text())
            self.assertEqual(float(steel["experiment"]["inner_diameter_m"]), 0.0048)
            self.assertEqual(float(silica["experiment"]["inner_diameter_m"]), 0.005)
            self.assertFalse(steel["ready_for_cfd"])
            self.assertIsNone(steel["full_column_helium_mass_flow_kg_s"])
            source = steel["source_terms"][0]
            self.assertAlmostEqual(source["full_column_amount_mol"], 3.167e-5, delta=1e-8)
            self.assertAlmostEqual(source["wedge_amount_mol"] * 72, source["full_column_amount_mol"])
            self.assertAlmostEqual(source["mean_release_mol_s"] * 10800, source["full_column_amount_mol"])
            self.assertEqual(source["measured_peak_K"], 641.15)
            lbe = json.loads((Path(work) / "LBE-I_SS_II/benchmark.json").read_text())
            self.assertEqual(lbe["excluded_KI_iodine_yield_pct"], 6.6)
            atoms = lbe["modeled_full_column_element_amounts_mol"]
            self.assertAlmostEqual(atoms["I"], 2 * atoms["Pb"] + 3 * atoms["Bi"])
            self.assertFalse(next(row for row in lbe["source_terms"] if row["species"] == "KI")["included_in_pbbii_model"])
            self.assertEqual((Path(work) / "PbI2_SS/wall_temperature.csv").read_text().strip(), "x_m,T_K")
            comparisons = PREPARE.read_csv(Path(work) / "comparison.csv")
            self.assertEqual(len(comparisons), 6)
            self.assertTrue(all(row["status"] == "NOTRUN" for row in comparisons))

    def prediction_file(self, path, rows):
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["experiment", "species", "peak_temperature_K"])
            writer.writerows(rows)

    def test_science_target_and_uncertainty_are_separate(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "peaks.csv"
            self.prediction_file(path, [("BiI3_SiO2_800", "BiI3", 509.15)])
            rows = PREPARE.compare(path, self.experiments, self.targets)
            row = next(row for row in rows if row["status"] == "REPORTED")
            self.assertAlmostEqual(row["difference_K"], 30)
            self.assertFalse(row["within_experimental_uncertainty"])
            self.assertTrue(row["within_science_target"])
            self.assertEqual(sum(row["status"] == "NOTRUN" for row in rows), 5)

    def test_invalid_predictions_refused(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "peaks.csv"
            for value in ("nan", "inf", "0", "-1"):
                self.prediction_file(path, [("PbI2_SS", "PbI2", value)])
                with self.assertRaises(ValueError):
                    PREPARE.compare(path, self.experiments, self.targets)
            for rows in ([('unknown', 'PbI2', 641.15)],
                         [('PbI2_SS', 'PbI2', 641.15), ('PbI2_SS', 'PbI2', 641.15)]):
                self.prediction_file(path, rows)
                with self.assertRaises(ValueError):
                    PREPARE.compare(path, self.experiments, self.targets)

    def test_flow_reference_requires_both_parameters(self):
        result = subprocess.run([sys.executable, str(HERE / "prepare.py"), "--check",
                                 "--flow-reference-temperature", "273.15"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("specify both", result.stderr)

    def test_explicit_mass_flow_and_preserved_predictions(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "predictions.csv"
            self.prediction_file(path, [("PbI2_SS", "PbI2", 641.15)])
            original = path.read_bytes()
            subprocess.run([sys.executable, str(HERE / "prepare.py"), "--out", work,
                            "--predictions", str(path), "--flow-reference-temperature", "273.15",
                            "--flow-reference-pressure", "101325"], check=True, capture_output=True)
            self.assertEqual(path.read_bytes(), original)
            manifest = json.loads((Path(work) / "PbI2_SS/benchmark.json").read_text())
            mass_flow = manifest["full_column_helium_mass_flow_kg_s"]
            self.assertAlmostEqual(mass_flow, 2.976e-7, delta=1e-10)
            self.assertAlmostEqual(manifest["wedge_helium_mass_flow_kg_s"] * 72, mass_flow)
            rows = PREPARE.read_csv(Path(work) / "comparison.csv")
            self.assertEqual(sum(row["status"] == "REPORTED" for row in rows), 1)


if __name__ == "__main__":
    unittest.main()
