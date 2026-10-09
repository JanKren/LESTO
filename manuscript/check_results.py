#!/usr/bin/env python3
"""Check saved-result arithmetic, provenance and the compiled review manuscript."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / 'thermochemistry/benchmarks/liu2025'))
from report import rows


def profile(root, name, species):
    filename = {'I': 'element_I', 'PbI2': 'PbI2_g', 'BiI3': 'BiI3_g'}[species]
    data = rows(root / 'profiles' / name / (filename + '_native.dat'))
    return np.array([r['x'] for r in data]), np.array([
        r['wall'] * (r['x_hi'] - r['x_lo']) for r in data])


def main():
    full = REPO / 'run/034-liu2025/full_analysis'
    follow = REPO / 'run/034-liu2025/followup_results'
    original = json.loads((full / 'analysis.json').read_text())
    result = json.loads((follow / 'analysis.json').read_text())
    assert (original['PASS'], original['FAIL'], original['NOTRUN']) == (15, 0, 1)
    assert result['completed_cases'] == 64 and result['execution_status'] == 'COMPLETE'
    assert not result['missing_or_failed']
    by_name = {r['case']: r for r in result['case_metrics']}
    nonzero = 0
    zero = 0
    for record in result['comparisons']:
        name = record['case']
        x, current = profile(follow, name, record['species'])
        xb, baseline = profile(full, by_name[name]['base_case'], record['species'])
        assert np.allclose(x, xb), name
        if baseline.sum():
            value = np.abs(current - baseline).sum() / baseline.sum()
            assert np.isclose(value, record['profile_L1_relative_to_base'], rtol=2e-12, atol=1e-15), name
            nonzero += 1
        else:
            assert record['profile_L1_relative_to_base'] is None, name
            zero += 1
    assert (nonzero, zero) == (180, 12)
    stoichiometric = 0
    for root, metrics in [(full, original['case_metrics']), (follow, result['case_metrics'])]:
        for case in metrics:
            _, iodine = profile(root, case['case'], 'I')
            _, lead = profile(root, case['case'], 'PbI2')
            _, bismuth = profile(root, case['case'], 'BiI3')
            defect = np.max(np.abs(iodine - 2 * lead - 3 * bismuth)) / max(iodine.max(), 1e-300)
            assert defect < 1e-12, (case['case'], defect)
            stoichiometric += 1
    hashes = 0
    for root, filename in [(follow, 'artifact_sha256.json'), (HERE, 'generated/artifact_sha256.json')]:
        for name, expected in json.loads((root / filename).read_text()).items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
            hashes += 1
    generated = json.loads((HERE / 'generated/results.json').read_text())
    for name, expected in generated['input_sha256'].items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected, name
        hashes += 1
    source = '\n'.join(p.read_text() for p in [HERE / 'main.tex', *sorted((HERE / 'sections').glob('*.tex'))])
    assert not re.search(r'\\(?:todo|decide)\b', source)
    log = (HERE / 'build/main.log').read_text()
    assert not re.search(r'undefined|Overfull|LaTeX Error|^!', log, re.MULTILINE), 'LaTeX issue'
    bibliography = (HERE / 'build/main.blg').read_text()
    assert 'Warning' not in bibliography and 'error' not in bibliography.lower(), 'Bibliography issue'
    pdf = HERE / 'lesto-iodine-methods.pdf'
    pdf_info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)', pdf_info, re.MULTILINE).group(1))
    text = subprocess.check_output(['pdftotext', str(pdf), '-'], text=True)
    assert '[?]' not in text and '[TODO' not in text and '[DECIDE' not in text
    checks = {
        'status': 'PASS', 'baseline_cases': 15, 'followup_cases': 64,
        'independently_recomputed_profile_comparisons': nonzero,
        'undefined_zero_reference_comparisons_confirmed': zero,
        'iodine_stoichiometry_profiles_checked': stoichiometric,
        'input_and_artifact_hashes_checked': hashes, 'pdf_pages': pages,
        'undefined_latex_references': 0, 'overfull_latex_boxes': 0,
        'bibliography_warnings': 0, 'draft_placeholders': 0,
        'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Saved-results and manuscript consistency checks; not new CFD runs or experimental validation.'}
    (HERE / 'review-checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
