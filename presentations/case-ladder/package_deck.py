#!/usr/bin/env python3
"""Extract the speaker notes of main.tex and package the deck source.

Writes speaker-notes.txt (one entry per content frame, with its PDF page;
the four section dividers are automatic pages without notes),
deck_manifest.json (sha256 of every file) and lesto-case-ladder-source.zip.
"""
import datetime
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PDF = HERE / "lesto-case-ladder.pdf"
FRAMES, PAGES = 30, 34


def main():
    source = (HERE / "main.tex").read_text()
    body = source.split(r"\begin{document}", 1)[1]
    # walk frames and section dividers in order to know each frame's page
    tokens = re.finditer(r"\\section\{|\\begin\{frame\}(.*?)\\end\{frame\}", body, re.S)
    notes, page = [], 0
    for t in tokens:
        page += 1
        if t.group(0).startswith(r"\section"):
            continue
        frame = t.group(1)
        title = re.match(r"(?:\[[^]]*\])?\{([^}]+)\}", frame)
        title = title.group(1) if title else "The LESTO case ladder (title)"
        note = re.search(r"\\note\{([^{}]*)\}", frame, re.S)
        assert note, (page, title)
        text = " ".join(note.group(1).split()).replace("--", "–").replace(r"\%", "%")
        notes.append(f"Page {page:02d}. {title}\n\n{text}\n")
    assert len(notes) == FRAMES and page == PAGES, (len(notes), page)
    (HERE / "speaker-notes.txt").write_text(
        "The LESTO case ladder — presenter notes\nJan Kren, 15 October 2026\n"
        "(pages 4, 12, 17 and 29 are the section dividers)\n\n" + "\n".join(notes))
    info = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)$", info, re.M)[1])
    assert pages == PAGES, pages
    files = [p for p in HERE.iterdir() if p.is_file() and not p.name.startswith(".")
             and (p.suffix in (".tex", ".sty", ".py", ".sh", ".md", ".txt", ".json") or p == PDF)
             and p.name != "deck_manifest.json"]
    for folder in ("figures", "psi-assets", "runs"):
        files += [p for p in (HERE / folder).rglob("*") if p.is_file()]
    bundle = HERE / "lesto-case-ladder-source.zip"
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(files):
            z.write(p, "case-ladder/" + str(p.relative_to(HERE)))
    manifest = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "pdf_pages": pages, "content_frames": FRAMES, "speaker_notes": len(notes),
        "theme": "copied unchanged from presentations/m10-project-review",
        "build_engine": "Tectonic 0.17.0 / XeTeX",
        "files_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted(files + [bundle])},
    }
    (HERE / "deck_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"{pages} pages, {len(notes)} speaker notes, {len(files)} files packaged")


if __name__ == "__main__":
    main()
