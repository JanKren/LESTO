"""Workflow tests with fake solver outputs; not native CFD qualification."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import run_study as r

CONTROL='startFrom latestTime; startTime 0; endTime 60; deltaT .002; writeInterval 500; purgeWrite 2; speciesTransport {tolerance 1e-16; nCorr 3;}'


class RunnerTests(unittest.TestCase):
    def test_control_only_changes_requested_settings(self):
        text=r.control(CONTROL,8)
        self.assertIn('nCorr 8;',text);self.assertIn('startFrom startTime;',text)
        self.assertIn('tolerance 1e-16;',text);self.assertIn('endTime 60.0;',text)

    def test_duplicate_and_missing_settings_fail(self):
        for text in [CONTROL+' nCorr 3;',CONTROL.replace('deltaT .002;','')]:
            with self.assertRaises(RuntimeError):r.control(text)

    def fixture(self,out):
        inputs=out/'inputs/case'
        for folder in ['0','constant','system']:(inputs/folder).mkdir(parents=True)
        (inputs/'0/Y_test').write_text('internalField uniform 0;')
        (inputs/'constant/physics').write_text('unchanged source')
        (inputs/'system/controlDict').write_text(CONTROL)
        r.save(inputs/'benchmark.json',{'end_time_s':60})
        return inputs,{'case':'case','correctors':8,'input_sha256':r.tree(inputs)}

    def test_smoke_then_full_and_completed_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);inputs,entry=self.fixture(out);before=r.tree(inputs);state={'cases':{}}
            def fake_solver(case,solver): (case/'log.solver').write_text('End\n')
            with patch.object(r,'solver_call',side_effect=fake_solver) as calls,patch.object(r,'analyse',return_value=({'status':'PASS'},[])),patch.object(r,'ledger_audit',return_value={}),patch.object(r,'log_audit',return_value={}):
                r.run_case(out,entry,state);self.assertEqual(calls.call_count,2)
                r.run_case(out,entry,state);self.assertEqual(calls.call_count,2)
            self.assertEqual(state['cases']['case']['status'],'PASS');self.assertEqual(r.tree(inputs),before)
            self.assertEqual(r.read(out/'smoke/case/benchmark.json')['end_time_s'],.02)
            self.assertEqual(r.read(out/'cases/case/benchmark.json')['end_time_s'],60)

    def test_failed_smoke_never_launches_full_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);inputs,entry=self.fixture(out);state={'cases':{}}
            with patch.object(r,'solver_call'),patch.object(r,'analyse',return_value=({'status':'FAIL'},[])):
                with self.assertRaises(RuntimeError):r.run_case(out,entry,state)
            self.assertFalse((out/'cases/case').exists());self.assertEqual(state['cases']['case']['status'],'FAIL')

    def test_tampered_input_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);inputs,entry=self.fixture(out)
            (inputs/'constant/physics').write_text('changed')
            with self.assertRaises(RuntimeError):r.run_case(out,entry,{'cases':{}})

    def test_interrupted_case_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);inputs,entry=self.fixture(out);old=out/'cases/case';old.mkdir(parents=True)
            (old/'evidence.txt').write_text('keep me')
            with patch.object(r,'solver_call',side_effect=RuntimeError('test stop')):
                with self.assertRaises(RuntimeError):r.run_case(out,entry,{'cases':{}})
            self.assertEqual(len(list((out/'interrupted').glob('*/evidence.txt'))),1)


if __name__=='__main__':unittest.main()
