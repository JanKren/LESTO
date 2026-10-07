"""Independent convergence examples and preservation of experimental source rates."""
import json
from pathlib import Path
import tempfile
import unittest

from followup_study import change_control, geometry_text
from report_followup import richardson


class FollowupEvidence(unittest.TestCase):
    def test_known_first_order_sequence_at_different_scales(self):
        # f(h)=2+0.3h, h=(2.25,1.5,1), known continuum limit=2.
        for scale in [1.,1e-200]:
            r=richardson(2.675*scale,2.45*scale,2.3*scale)
            self.assertAlmostEqual(r['observed_order'],1)
            self.assertAlmostEqual(r['extrapolated_value']/scale,2)
            self.assertAlmostEqual(r['fine_GCI_fraction'],1.25*.3/2.3)

    def test_unqualified_sequences_are_not_extrapolated(self):
        for values,status in [((1,1.2,1.1),'OSCILLATORY_NO_EXTRAPOLATION'),
                              ((1,2,4),'NONDECREASING_NO_EXTRAPOLATION'),
                              ((1,1,1),'UNRESOLVED_AT_ROUNDOFF')]:
            r=richardson(*values)
            self.assertEqual(r['status'],status)
            self.assertIsNone(r['fine_GCI_fraction'])

    def test_short_window_does_not_compress_experimental_inventory(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory);(case/'system').mkdir();(case/'constant').mkdir()
            (case/'system/controlDict').write_text('deltaT .002; endTime 60; writeInterval 500; purgeWrite 2;')
            source=b'samples {boat {amounts {PbI2_g 6e-12;} startTime 0; duration 60;}}'
            (case/'constant/thermochemistryProperties').write_bytes(source)
            (case/'benchmark.json').write_text(json.dumps({'dt_s':.002,'end_time_s':60,'source_wedge_mol_s':{'PbI2_g':1e-13}}))
            change_control(case,.001,.01)
            self.assertEqual((case/'constant/thermochemistryProperties').read_bytes(),source)
            self.assertEqual(json.loads((case/'benchmark.json').read_text())['source_wedge_mol_s'],{'PbI2_g':1e-13})

    def test_boat_shift_preserves_radial_geometry_and_amount(self):
        for text in ['box (.09 .0006 -1) (.105 .0018 1); amount 7e-12;',
                     'box ((.09 .0006 -1) (.105 .0018 1)); amount 7e-12;']:
            result=geometry_text(text,.075,.09)
            self.assertIn('.0006 -1',result);self.assertIn('.0018 1',result)
            self.assertIn('amount 7e-12;',result)
        with self.assertRaises(AssertionError):geometry_text('box (.1 0 -1) (.2 1 1);',.075,.09)


if __name__=='__main__':unittest.main()
