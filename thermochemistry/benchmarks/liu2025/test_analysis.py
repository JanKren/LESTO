"""Independent small examples for the scientific comparison conventions."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from analyse_matrix import experimental_shape, ledger_audit, moments, relative_l1


class AnalysisConventions(unittest.TestCase):
    def test_histogram_transport_and_overlap(self):
        # Half the probability moves exactly two 1 cm bins: W1 = 1 cm.
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            (case/'iodine-comparison.csv').write_text(
                'x_m,predicted_relative_yield_percent,figure_relative_yield_percent\n'
                '0,50,0\n.01,50,50\n.02,0,50\n')
            r=experimental_shape(case)
            self.assertAlmostEqual(r['raw_L1_percentage_points'],100)
            self.assertAlmostEqual(r['normalised_scan_overlap'],.5)
            self.assertAlmostEqual(r['normalised_scan_total_variation'],.5)
            self.assertAlmostEqual(r['normalised_scan_Wasserstein_m'],.01)
            self.assertAlmostEqual(r['normalised_scan_centroid_difference_m'],-.01)

    def test_shape_normalisation_is_explicitly_conditional(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            (case/'iodine-comparison.csv').write_text(
                'x_m,predicted_relative_yield_percent,figure_relative_yield_percent\n'
                '0,20,40\n.01,30,60\n')
            r=experimental_shape(case)
            self.assertAlmostEqual(r['raw_L1_percentage_points'],50)
            self.assertAlmostEqual(r['normalised_scan_total_variation'],0)
            self.assertAlmostEqual(r['model_scan_sum_pct'],50)
            self.assertAlmostEqual(r['observed_scan_sum_pct'],100)

    def test_corrected_ledger_does_not_hide_raw_material_residue(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            for element in ['Pb','Bi','I']:
                path=case/'postProcessing/phaseChangeElements/0'/(element+'.dat')
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text('# columns time initial released transport removed gas wall sample closure clamped solverDefect\n'
                                '60 0 1 0 0 .1 .895 0 0 0 .005\n')
            r=ledger_audit(case)
            for element in r.values():
                self.assertEqual(element['max_accounted_closure_over_source'],0)
                self.assertAlmostEqual(element['max_raw_material_residue_over_source'],.005)

    def test_equal_shape_does_not_imply_equal_deposit(self):
        a=np.array([1.,1.]);b=np.array([2.,2.])
        self.assertAlmostEqual(relative_l1(a,b),.5)
        self.assertAlmostEqual(relative_l1(a/a.sum(),b/b.sum()),0)
        r=moments(np.array([0.,.02]),b)
        self.assertAlmostEqual(r['centroid_m'],.01)
        self.assertAlmostEqual(r['sigma_m'],.01)


if __name__=='__main__':unittest.main()
