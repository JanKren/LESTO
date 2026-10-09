#!/usr/bin/env python3
"""Package the working manuscript and saved results, excluding raw CFD inputs."""
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def main():
    archive = HERE / 'lesto-iodine-review.zip'
    manifest = HERE / 'review-package-manifest.json'
    files = {HERE / name for name in [
        'main.tex', 'arxiv.sty', 'references.bib', 'README.md', 'review-notes.md',
        'prepare_results.py', 'check_results.py', 'build.sh', 'package_review.py', 'lesto-iodine-methods.pdf',
        'notes-for-authors.md', 'literature-review.md',
        'review-checks.json']}
    for folder in ['sections', 'figures', 'generated']:
        files.update(p for p in (HERE / folder).iterdir() if p.is_file())
    for result in ['full_analysis', 'followup_results']:
        root = REPO / 'run/034-liu2025' / result
        files.update(p for p in root.glob('*') if p.is_file() and not p.name.startswith('.'))
        files.update(p for p in (root / 'profiles').rglob('*') if p.is_file())
    files.update(REPO / 'run/034-liu2025' / name for name in [
        'full_size_launch.json', 'full_size_completion.json', 'followup_launch.json'])
    bench = REPO / 'thermochemistry/benchmarks/liu2025'
    files.update(bench / name for name in [
        'analyse_matrix.py', 'report.py', 'prepare.py', 'carriers.py', 'README.md',
        'experiments.csv', 'targets.csv', 'sublimation.csv'])
    files.update(p for p in (bench / 'digitized').iterdir() if p.suffix in {'.csv', '.json'})
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(files)}
    manifest.write_text(json.dumps({
        'scope': 'Working manuscript and saved results; not a complete CFD input distribution.',
        'excluded': ['Source-paper PDFs', 'Raw CFD time directories', 'External thermodynamic coefficients'],
        'sha256': hashes}, indent=2) + '\n')
    files.add(manifest)
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
        for path in sorted(files):
            output.write(path, path.relative_to(REPO))
    print(f'Packaged {len(files)} files in {archive.name} ({archive.stat().st_size:,} bytes).')


if __name__ == '__main__':
    main()
