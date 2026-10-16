#!/usr/bin/env python3
"""Write speaker-notes.txt from the \\note{...} of every frame of main.tex.

Each entry carries the PDF page; the backup divider is an automatic page
without a note.  The minute budget at the start of each main-part note is
summed and printed (target: 30 minutes).
"""
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
PDF = HERE / "lesto-partner-update.pdf"


def main():
    body = (HERE / "main.tex").read_text().split(r"\begin{document}", 1)[1]
    tokens = re.finditer(r"\\section\{|\\begin\{frame\}(.*?)\\end\{frame\}", body, re.S)
    notes, page, minutes = [], 0, 0.0
    for t in tokens:
        page += 1
        if t.group(0).startswith(r"\section"):
            continue
        frame = t.group(1)
        title = re.match(r"(?:\[[^]]*\])?\{([^}]+)\}", frame)
        title = title.group(1) if title else "CFD-GEMS3K modelling of iodine transport and deposition (title)"
        note = re.search(r"\\note\{([^{}]*)\}", frame, re.S)
        assert note, (page, title)
        text = " ".join(note.group(1).split()).replace("--", "–").replace(r"\%", "%")
        m = re.match(r"\((\d+(?:\.\d+)?) min\)", text)
        if m:
            minutes += float(m.group(1))
        title = title.replace(r"\PbI{}", "PbI2").replace(r"\BiI{}", "BiI3")
        notes.append(f"Page {page:02d}. {title}\n\n{text}\n")
    info = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)$", info, re.M)[1])
    assert pages == page, (pages, page)
    (HERE / "speaker-notes.txt").write_text(
        "CFD-GEMS3K modelling of iodine transport and deposition - presenter notes\n"
        "Jan Kren, LESTO Task 3.4 partner update, 15 October 2026\n"
        f"Main part: about {minutes:g} minutes (minute budget in brackets); "
        "pages after the backup divider are not presented by default.\n\n" + "\n".join(notes))
    print(f"{pages} pages, {len(notes)} notes, main part {minutes:g} min")


if __name__ == "__main__":
    main()
