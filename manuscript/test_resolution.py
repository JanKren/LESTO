"""Analytical checks for the saved-profile reassessment."""
import unittest
import numpy as np
from assess_resolution import aggregate, relative_l1, shape_metrics


class ResolutionTests(unittest.TestCase):
    def test_aggregation_keeps_edge_inventory(self):
        a=np.arange(1.,8.)
        for width in (1,2,4):
            for offset in range(width):
                self.assertEqual(aggregate(a,width,offset).sum(),28.)

    def test_opposite_bin_changes_can_cancel(self):
        a,b=np.array([2.,0.,0.,2.]),np.ones(4)
        self.assertEqual(relative_l1(a,b),100.)
        self.assertEqual(relative_l1(aggregate(a,2,0),aggregate(b,2,0)),0.)

    def test_shape_metric_ignores_amplitude(self):
        a=np.array([1.,2.,1.]); x=np.array([0.,.01,.02])
        self.assertEqual(shape_metrics(3*a,a,x),{'total_variation_pct':0.,'binned_W1_cm':0.})
        self.assertEqual(relative_l1(3*a,a),200.)

    def test_displaced_point_mass(self):
        metrics=shape_metrics(np.array([0.,0.,1.]),np.array([1.,0.,0.]),np.array([0.,.01,.02]))
        self.assertEqual(metrics['total_variation_pct'],100.)
        self.assertAlmostEqual(metrics['binned_W1_cm'],2.)


if __name__=='__main__': unittest.main()
