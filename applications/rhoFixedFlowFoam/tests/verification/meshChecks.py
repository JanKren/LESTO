#!/usr/bin/env python3
"""
tests/verification/meshChecks.py

PURPOSE

Checks of the refined pipe meshes of milestone M7 (tests/verification/
pipeMeshes): the h/2 first-layer mesh of acceptance 2 and the axial x2
mesh, each compared with the base mesh of the tests (gmshToFoam of
run/001.../pipe.msh).  The facts are the lines of pipeMeshFacts
(tests/verification/pipeMeshFacts) and of checkMesh.  Every subcommand
prints one summary line and exits with 0 (criterion met) or 1.

  wallH2 <base case> <refined case>
      refineWallLayer '(WALL)' 0.5: 'Mesh OK'; cells = base + wall faces
      (one more layer); the same wall faces and wall area (1e-12); the
      wall-normal thickness of the first layer half of the base's (1e-9)
      in every element (min and max); A_I/V_I of layer faceCells about
      twice the base's (1.9 to 2.1; printed); the cells within 1e-4 m of
      the wall by the centre distance and by the patchWave distance one
      layer more than the base (3 instead of 2) and their volume unchanged
      (1e-12; the distance layer covers the same region)
  axial2 <base case> <refined case>
      refineMesh directions (tan1) of every cell: 'Mesh OK'; cells and
      wall faces twice the base's; twice the axial positions with half the
      axial extent (1e-9); first layer, A_I/V_I and the distance layer
      (cells twice, volume equal) unchanged (1e-9)
  table <case> ...
      the facts of the meshes as a markdown table (always 0)
"""

import os
import re
import sys


def done(ok, text):
    print(('OK: ' if ok else 'FAILED: ') + text)
    sys.exit(0 if ok else 1)


def _numbers(words):
    out = {}
    for key, value in zip(words[0::2], words[1::2]):
        try:
            out[key] = float(value)
        except ValueError:
            out[key] = value
    return out


def facts(case):
    """the facts of log.pipeMeshFacts and log.checkMesh of a case"""
    f = {'case': case}
    with open(os.path.join(case, 'log.pipeMeshFacts')) as log:
        for line in log:
            w = line.split()
            if not w:
                continue
            if w[0] == 'cells':
                f['cells'] = int(w[1])
            elif w[0] == 'patch' and len(w) >= 6:
                f['wallFaces'] = int(w[3])
                f['wallArea'] = float(w[5])
            elif w[0] == 'faceCells' and w[1] == 'elements':
                f.update({'faceCells.' + k: v
                          for k, v in _numbers(w[1:]).items()})
            elif w[0] == 'firstLayer' or w[0] == 'axialExtent':
                name = w[0] if w[0] == 'axialExtent' else w[1]
                f[name] = _numbers(w[2:] if w[0] == 'firstLayer' else w[1:])
            elif w[0] == 'axialPositions':
                f['axialPositions'] = int(w[1])
            elif w[0] == 'layer':
                f.setdefault('layers', []).append(_numbers(w[2:]))
            elif w[0] == 'distanceLayer' and w[1].startswith('centre('):
                f['distance.centre'] = _numbers(w[2:])
            elif w[0] == 'distanceLayer' and w[1].startswith('patchWave('):
                f['distance.patchWave'] = _numbers(w[2:])
            elif w[0] == 'distanceLayer':
                # 'distanceLayer centre cells per layer index (1..6, beyond)'
                f['distance.perLayer'] = [int(v) for v in w[8:]]
    with open(os.path.join(case, 'log.checkMesh')) as log:
        text = log.read()
    f['meshOK'] = bool(re.search(r'^Mesh OK\.', text, re.M))
    m = re.search(r'Max aspect ratio = ([0-9.eE+-]+)', text)
    f['aspectRatio'] = float(m.group(1)) if m else float('nan')
    m = re.search(r'non-orthogonality Max: ([0-9.eE+-]+) average: '
                  r'([0-9.eE+-]+)', text)
    f['nonOrthoMax'] = float(m.group(1)) if m else float('nan')
    f['nonOrthoAverage'] = float(m.group(2)) if m else float('nan')
    m = re.search(r'Max skewness = ([0-9.eE+-]+)', text)
    f['skewness'] = float(m.group(1)) if m else float('nan')
    return f


def rel(a, b):
    return abs(a - b)/abs(b) if b else abs(a - b)


def thickness(f, key='min'):
    return f['wallNormalThickness'][key]


def describe(f):
    return (f'{f["cells"]} cells; non-orthogonality max {f["nonOrthoMax"]:.4g}'
            f' (average {f["nonOrthoAverage"]:.4g}), skewness '
            f'{f["skewness"]:.4g}, aspect ratio {f["aspectRatio"]:.4g}; first'
            f' layer {thickness(f)*1e6:.4f} um (half cell '
            f'{f["halfCellDistance(1/deltaCoeffs)"]["min"]*1e6:.4f} um);'
            f' A_I/V_I {f["faceCells.A_I/V_I"]:.6g} 1/m; within 1e-4 m: '
            f'{int(f["distance.centre"]["cells"])} cells by centre distance'
            f' ({f["distance.centre"]["layers"]:g} layers), '
            f'{int(f["distance.patchWave"]["cells"])} by patchWave distance,'
            f' A_I/V_layer {f["distance.centre"]["A_I/V_layer"]:.6g} 1/m')


def cmdWallH2(baseCase, case):
    b, r = facts(baseCase), facts(case)
    checks = {
        'Mesh OK': r['meshOK'],
        'cells = base + wall faces': r['cells'] == b['cells'] + b['wallFaces'],
        'wall faces and area unchanged': r['wallFaces'] == b['wallFaces']
        and rel(r['wallArea'], b['wallArea']) <= 1e-12,
        'first layer = base/2': all(
            rel(thickness(r, k), 0.5*thickness(b, k)) <= 1e-9
            for k in ('min', 'max')),
        'A_I/V_I about 2x': 1.9 <= r['faceCells.A_I/V_I']
        / b['faceCells.A_I/V_I'] <= 2.1,
        'one more layer within 1e-4 m': all(
            r['distance.' + k]['cells'] == b['distance.' + k]['cells']
            + b['wallFaces'] for k in ('centre', 'patchWave')),
        'distance-layer volume unchanged': all(
            rel(r['distance.' + k]['volume'], b['distance.' + k]['volume'])
            <= 1e-12 for k in ('centre', 'patchWave')),
    }
    failed = [k for k, v in checks.items() if not v]
    done(not failed,
         f'h/2 mesh: {describe(r)}; A_I/V_I ratio '
         f'{r["faceCells.A_I/V_I"]/b["faceCells.A_I/V_I"]:.6f} to the base'
         f' ({describe(b)})'
         + (f'; failed: {", ".join(failed)}' if failed else ''))


def cmdAxial2(baseCase, case):
    b, r = facts(baseCase), facts(case)
    checks = {
        'Mesh OK': r['meshOK'],
        'cells 2x': r['cells'] == 2*b['cells'],
        'wall faces 2x, area unchanged': r['wallFaces'] == 2*b['wallFaces']
        and rel(r['wallArea'], b['wallArea']) <= 1e-12,
        'axial positions 2x': r['axialPositions'] == 2*b['axialPositions'],
        'axial extent /2': all(
            rel(r['axialExtent'][k], 0.5*b['axialExtent'][k]) <= 1e-9
            for k in ('min', 'max')),
        'first layer unchanged': all(
            rel(thickness(r, k), thickness(b, k)) <= 1e-9
            for k in ('min', 'max')),
        'A_I/V_I unchanged': rel(r['faceCells.A_I/V_I'],
                                 b['faceCells.A_I/V_I']) <= 1e-9,
        'distance layer: 2x cells, same volume': all(
            r['distance.' + k]['cells'] == 2*b['distance.' + k]['cells']
            and rel(r['distance.' + k]['volume'],
                    b['distance.' + k]['volume']) <= 1e-12
            for k in ('centre', 'patchWave')),
    }
    failed = [k for k, v in checks.items() if not v]
    done(not failed,
         f'axial x2 mesh: {describe(r)}; {r["axialPositions"]} axial cells of'
         f' {r["axialExtent"]["min"]*1e3:.6f} mm (base {b["axialPositions"]}'
         f' of {b["axialExtent"]["min"]*1e3:.6f} mm)'
         + (f'; failed: {", ".join(failed)}' if failed else ''))


def cmdTable(*cases):
    print('| mesh | cells | non-orth. max (avg) | skewness | aspect ratio |'
          ' first layer [um] | 1/deltaCoeffs [um] | A_I/V_I faceCells [1/m] |'
          ' axial cells | cells (layers) with R - r_c <= 1e-4 m |'
          ' cells with patchWave d <= 1e-4 m | A_I/V_layer [1/m] |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for c in cases:
        f = facts(c)
        print(f'| {os.path.basename(c.rstrip("/"))} | {f["cells"]} | '
              f'{f["nonOrthoMax"]:.2f} ({f["nonOrthoAverage"]:.2f}) | '
              f'{f["skewness"]:.3f} | {f["aspectRatio"]:.2f} | '
              f'{thickness(f)*1e6:.3f} | '
              f'{f["halfCellDistance(1/deltaCoeffs)"]["min"]*1e6:.3f} | '
              f'{f["faceCells.A_I/V_I"]:.2f} | {f["axialPositions"]} | '
              f'{int(f["distance.centre"]["cells"])} '
              f'({f["distance.centre"]["layers"]:g}) | '
              f'{int(f["distance.patchWave"]["cells"])} | '
              f'{f["distance.centre"]["A_I/V_layer"]:.2f} |')
    sys.exit(0)


COMMANDS = {'wallH2': cmdWallH2, 'axial2': cmdAxial2, 'table': cmdTable}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    COMMANDS[sys.argv[1]](*sys.argv[2:])
