"""Failure-path checks for scientific input provenance and numerical gates."""
import tempfile
import unittest
from pathlib import Path

from carriers import prepare_carrier
from cases import source_rates
from prepare import load_references
from report import numerical_gates


class Validation(unittest.TestCase):
    def ledger(self, case, time=60, closure=0, defect=0):
        for atom in ['Pb','Bi','I']:
            path=case/'postProcessing/phaseChangeElements/0'/(atom+'.dat')
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(f'# columns time reference closure solverDefect\n{time} 1e-9 {closure} {defect}\n')

    def test_truncated_run_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            self.ledger(case,time=1)
            with self.assertRaisesRegex(ValueError,'Incomplete run'):
                numerical_gates(case,60)

    def test_conservation_and_defect_fail_independently(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            self.ledger(case)
            self.assertTrue(all(v['status']=='PASS' for v in numerical_gates(case,60).values()))
            self.ledger(case,closure=2e-21)
            self.assertTrue(all(v['status']=='FAIL' for v in numerical_gates(case,60).values()))
            self.ledger(case,defect=2e-17)
            self.assertTrue(all(v['status']=='FAIL' for v in numerical_gates(case,60).values()))

    def test_nonfinite_output_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory);self.ledger(case,closure='nan')
            with self.assertRaisesRegex(ValueError,'nonfinite'):
                numerical_gates(case,60)

    def test_mean_release_uses_three_hours_and_sector(self):
        experiments,targets,species=load_references()
        rates=source_rates(experiments['PbI2_SS'],targets,species,1/72)
        self.assertAlmostEqual(rates['PbI2_g']*72*10800*.46100894,14.6e-6,places=14)

    def test_profile_assumptions_and_missing_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)/'physical'
            result=prepare_carrier(case,'PbI2_SS',384,48,298.15,101325)
            self.assertEqual(result['scientific_validation_status'],'NOTRUN')
            self.assertIn('ASSUMED',result['flow_reference_status'])
            self.assertGreater(result['temperature_coverage_m'][0],0)
            self.assertLess(result['temperature_coverage_m'][1],1.25)
            with self.assertRaisesRegex(ValueError,'exists'):
                prepare_carrier(case,'PbI2_SS',384,48,298.15,101325)
            with self.assertRaisesRegex(ValueError,'No published'):
                prepare_carrier(Path(directory)/'missing','BiI3_SiO2_800',384,48,298.15,101325)


if __name__=='__main__':
    unittest.main()
