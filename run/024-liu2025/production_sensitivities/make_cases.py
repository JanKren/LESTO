#!/usr/bin/env python3
"""Write the eight production-scheme sensitivity cases from the two production base cases.

    python3 make_cases.py <production cases> <out> [<tables>]

<production cases> holds PbI2_SS_Q1_baseline and LBE-I_SS_II_Q2_baseline_sat0.01 as run for
../production_576x72_limitedLinear (576 x 72, Gauss limitedLinear 1); <tables> is the baseline
table directory (default ../work/prepared/tables/baseline).  The changes are those of
thermochemistry/benchmarks/liu2025/followup_study.py for the upwind follow-up cases:

- __D5e-05, __D0.0002: all nine `D 0.0001;` of constant/speciesTransportProperties;
- __alpha0.1, __alpha0.01: all three `accommodation 1;` of constant/thermochemistryProperties
  (the PbI2_g and BiI3_g pairs and HKSCoeffs).

The table paths of the base are repointed to <tables>; nothing else changes.
"""
import json, re, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASES = ['PbI2_SS_Q1_baseline', 'LBE-I_SS_II_Q2_baseline_sat0.01']
VARIANTS = ['D5e-05', 'D0.0002', 'alpha0.1', 'alpha0.01']
TABLE_PATH = re.compile(r'"[^"]*/prepared/tables/baseline')

src, out = Path(sys.argv[1]), Path(sys.argv[2])
tables = Path(sys.argv[3] if len(sys.argv) > 3 else HERE.parent/'work/prepared/tables/baseline').resolve()
assert tables.is_dir(), tables
for base in BASES:
    for label in VARIANTS:
        c = out/f'{base}__{label}'
        assert not c.exists(), c
        c.mkdir(parents=True)
        for d in ['0', 'constant', 'system']:
            shutil.copytree(src/base/d, c/d, symlinks=False)
        shutil.copy(src/base/'benchmark.json', c/'benchmark.json')
        tp = c/'constant/thermochemistryProperties'
        t, n = TABLE_PATH.subn(f'"{tables}', tp.read_text()); assert n > 0; tp.write_text(t)
        m = json.loads((c/'benchmark.json').read_text()); m['table_directory'] = str(tables)
        if label.startswith('D'):
            D = float(label[1:]); p = c/'constant/speciesTransportProperties'
            t, n = re.subn(r'\bD\s+0\.0001;', f'D {D:.17g};', p.read_text()); assert n == 9, n; p.write_text(t)
            ch = {'constant_D_m2_s': D}
        else:
            a = float(label[5:]); p = c/'constant/thermochemistryProperties'
            t, n = re.subn(r'\baccommodation\s+1;', f'accommodation {a:.17g};', p.read_text()); assert n == 3, n
            p.write_text(t)
            ch = {'iodide_mass_accommodation': a,
                  'replaced': 'all three accommodation entries (PbI2_g, BiI3_g pairs and HKSCoeffs), as in followup_study.py'}
        m['followup'] = {'base_case': base, 'changes': ch, 'scheme': 'Gauss limitedLinear 1', 'cells': 576*72,
                         'origin': 'production-scheme sensitivity (review round, plan section 39)'}
        (c/'benchmark.json').write_text(json.dumps(m, indent=1))
        print('made', c.name)
