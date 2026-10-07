#!/usr/bin/env python3
"""Package the prepared Beamer source and extract speaker notes."""
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent


def main():
    source = (HERE / "main.tex").read_text()
    frames = re.findall(r"\\begin\{frame\}(.*?)\\end\{frame\}", source, re.S)
    assert len(frames) == 20, len(frames)
    notes = []
    for i, frame in enumerate(frames, 1):
        title = re.match(r"(?:\[[^]]*\])?\{([^}]+)\}", frame)
        title = title.group(1) if title else "Iodine deposition in lead–bismuth systems"
        # Notes are plain prose with no nested TeX groups in this deck.
        note = re.search(r"\\note\{([^{}]*)\}", frame, re.S)
        assert note, i
        plain_note = note.group(1).replace('--', '–').replace(r'\%', '%')
        notes.append(f"{i:02d}. {title}\n\n{plain_note}\n")
    (HERE / "speaker-notes.txt").write_text("PSI Beamer — presenter notes\n7 October 2026\n\n" + "\n".join(notes))
    pdf = HERE / "lesto-m10-review.pdf"
    metadata = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)$", metadata, re.M)[1])
    assert pages == 20, pages
    files = [p for p in HERE.iterdir() if p.is_file() and
             (p.suffix in (".tex", ".sty", ".py", ".sh", ".md", ".txt", ".json") or p.name == pdf.name)
             and not p.name.startswith(".") and p.name != "deck_manifest.json"]
    files += [p for folder in ("figures", "psi-assets") for p in (HERE / folder).rglob("*") if p.is_file()]
    bundle = HERE / "lesto-m10-review-source.zip"
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(files):
            archive.write(p, "m10-project-review/" + str(p.relative_to(HERE)))
    manifest = {"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "slides": 20, "speaker_notes": len(notes), "pdf_pages": pages,
                "template_origin": "/home/jan/Desktop/PSI-Beamer.zip",
                "template_archive_sha256": "628b010d5537d200095c1483db89ad0872ad0fcf6cdf19d763a53aadf6a9f304",
                "build_engine": "Tectonic 0.17.0 / XeTeX, offline cached resources",
                "files_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(files + [bundle])}}
    (HERE / "deck_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Packaged {len(files)} files, {pages} slides and {len(notes)} speaker notes.")


if __name__ == "__main__":
    main()
