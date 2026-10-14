#!/usr/bin/env python3
"""Compare a finished run of a case with the numbers of its expect.json.

Usage: lestoCheck.py <case directory> <expect.json> [<label>]

expect.json holds a list "checks"; every check names a file relative to the
case directory (or to the repository with the prefix "repo:", or to another
case of the ladder with "case:<slug>/") and what to compare:

  header    "file", "pattern": a regular expression with one group, matched
            against the '#' lines of a profile or ledger file
  log       "file", "pattern": the last match in a log
  json      "file", "key": a dotted path into a JSON file
     each of these with "value" and "abs" or "rel" (a number), or with
     "equals" (a string), or with "max" (an upper bound of the magnitude)
  closure   "file": a pair ledger postProcessing/phaseChangeBalance/0/<pair>.dat;
            max |closure| over every step relative to the largest
            initial + released, at most "max"
  elementClosure  "file": postProcessing/phaseChangeElements/0/<element>.dat;
            max |closure|/reference over every step, at most "max"
  sameFiles "pairs": [[file, reference], ...]; "rel": 0 (default) requires
            identical bytes, otherwise every number may differ by rel
  onsetBin  "file": an axial profile; the wall temperature of the hottest
            bin whose "column" (default wall) is at least "fraction" of its
            peak, among the bins with "Tmin" < Twall <= "Tmax" (optional)
  l1        "file", "reference": two axial profiles of the same bins;
            sum |a - b|/sum |b| of "column" (default wall)

A check with "optional": true whose files are missing is NOTRUN, not FAIL.
Prints one line per check and a summary line; exit status 1 on any FAIL.
"""
import json
import math
import re
import sys
from pathlib import Path

CASE_DIR = None
RUN_DIR = Path(__file__).resolve().parent
REPO = RUN_DIR.parent


def resolve(name):
    if name.startswith('repo:'):
        return REPO / name[5:]
    if name.startswith('case:'):
        slug, _, rest = name[5:].partition('/')
        matches = sorted(RUN_DIR.glob('[0-9][0-9][0-9]-' + slug))
        if len(matches) != 1:
            raise FileNotFoundError(f'case {slug}: {len(matches)} directories')
        return matches[0] / rest
    return CASE_DIR / name


def number(text):
    return float(text)


def compare(check, value):
    """(ok, text) for a found value."""
    if 'equals' in check:
        return str(value) == str(check['equals']), f"{value} (expected {check['equals']})"
    v = number(value)
    if not math.isfinite(v):
        return False, f'{value} (not finite)'
    if 'max' in check:
        return abs(v) <= check['max'], f'{v:.6g} (at most {check["max"]:.3g})'
    ref = check['value']
    if 'abs' in check:
        tol = check['abs']
        return abs(v - ref) <= tol, f'{v:.12g} (expected {ref:.12g} +- {tol:.3g})'
    tol = check.get('rel', 0.0)
    diff = abs(v - ref) / abs(ref) if ref else abs(v)
    return diff <= tol, f'{v:.12g} (expected {ref:.12g}, relative difference {diff:.2g}, at most {tol:.3g})'


def regex_value(path, pattern, header_only):
    found = None
    for line in path.read_text(errors='replace').splitlines():
        if header_only and not line.startswith('#'):
            continue
        m = re.search(pattern, line)
        if m:
            found = m.group(1)
            if header_only:
                return found
    if found is None:
        raise ValueError(f'no match of {pattern!r} in {path}')
    return found


def json_value(path, key):
    value = json.loads(path.read_text())
    for part in key.split('.'):
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value


def table(path):
    columns, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            words = line.lstrip('#').split()
            if words and words[0] in ('time', 'columns'):
                columns = words[1:] if words[0] == 'columns' else words
            continue
        if line.strip():
            rows.append([float(x) for x in line.split()])
    return columns, rows


def closure(path):
    columns, rows = table(path)
    c, i, r = columns.index('closure'), columns.index('initial'), columns.index('released')
    scale = max(abs(row[i] + row[r]) for row in rows)
    return max(abs(row[c]) for row in rows) / scale


def element_closure(path):
    columns, rows = table(path)
    c, ref = columns.index('closure'), columns.index('reference')
    return max(abs(row[c]) / row[ref] for row in rows if row[ref] > 0)


def same_files(pairs, rel):
    worst = 0.0
    for a, b in pairs:
        pa, pb = resolve(a), resolve(b)
        if rel == 0:
            if pa.read_bytes() != pb.read_bytes():
                return False, f'{a} differs from {b}'
            continue
        ta = [t for t in re.split(r'\s+', pa.read_text()) if t]
        tb = [t for t in re.split(r'\s+', pb.read_text()) if t]
        if len(ta) != len(tb):
            return False, f'{a}: {len(ta)} tokens, {b}: {len(tb)}'
        for x, y in zip(ta, tb):
            if x == y:
                continue
            try:
                fx, fy = float(x), float(y)
            except ValueError:
                return False, f'{a}: token {x!r} against {y!r}'
            d = abs(fx - fy) / max(abs(fx), abs(fy))
            worst = max(worst, d)
    if rel == 0:
        return True, f'{len(pairs)} file(s) identical'
    return worst <= rel, f'{len(pairs)} file(s), largest relative difference {worst:.2g} (at most {rel:.3g})'


def profile_columns(path):
    columns, rows = table(path)
    return {c: [row[i] for row in rows] for i, c in enumerate(columns)}


def onset_bin(path, check):
    p = profile_columns(path)
    values = p[check.get('column', 'wall')]
    peak = max(values)
    lo, hi = check.get('Tmin', -math.inf), check.get('Tmax', math.inf)
    hits = [T for T, v in zip(p['Twall'], values)
            if lo < T <= hi and v >= check['fraction'] * peak]
    if not hits:
        raise ValueError(f'no bin above {check["fraction"]} of the peak in {path}')
    return max(hits)


def l1(a, b, column):
    va, vb = profile_columns(a)[column], profile_columns(b)[column]
    if len(va) != len(vb):
        raise ValueError(f'{a}: {len(va)} bins, {b}: {len(vb)}')
    return sum(abs(x - y) for x, y in zip(va, vb)) / sum(abs(y) for y in vb)


def run_check(check):
    kind = check['type']
    if kind == 'sameFiles':
        return same_files(check['pairs'], check.get('rel', 0))
    path = resolve(check['file'])
    if kind == 'header':
        return compare(check, regex_value(path, check['pattern'], True))
    if kind == 'log':
        return compare(check, regex_value(path, check['pattern'], False))
    if kind == 'json':
        return compare(check, json_value(path, check['key']))
    if kind == 'closure':
        return compare(check, closure(path))
    if kind == 'elementClosure':
        return compare(check, element_closure(path))
    if kind == 'onsetBin':
        return compare(check, onset_bin(path, check))
    if kind == 'l1':
        return compare(check, l1(path, resolve(check['reference']), check.get('column', 'wall')))
    raise ValueError(f'unknown check type {kind}')


def main():
    global CASE_DIR
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    CASE_DIR = Path(sys.argv[1]).resolve()
    spec = json.loads(Path(sys.argv[2]).read_text())
    label = sys.argv[3] if len(sys.argv) > 3 else CASE_DIR.name
    counts = {'PASS': 0, 'FAIL': 0, 'NOTRUN': 0}
    for check in spec['checks']:
        try:
            ok, text = run_check(check)
            status = 'PASS' if ok else 'FAIL'
        except (OSError, FileNotFoundError) as error:
            status = 'NOTRUN' if check.get('optional') else 'FAIL'
            text = f'missing: {error}'
        except (ValueError, KeyError, IndexError) as error:
            status, text = 'FAIL', f'{type(error).__name__}: {error}'
        counts[status] += 1
        print(f'  {status:6s} {check["id"]}: {text}')
    status = 'FAIL' if counts['FAIL'] else 'PASS'
    print(f'{status} {label}: {counts["PASS"]} PASS, {counts["FAIL"]} FAIL, {counts["NOTRUN"]} NOTRUN')
    return 1 if counts['FAIL'] else 0


if __name__ == '__main__':
    sys.exit(main())
