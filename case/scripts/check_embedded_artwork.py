#!/usr/bin/env python3
"""Check the figures EMBEDDED in a manuscript, not just the files beside it.

check_journal_artwork validates the separate artwork files a journal asks for.
Nothing was checking the copies inside the .docx, and Word downscales an image to
whatever it was first inserted at — so a manuscript can sit beside a set of
spec-perfect files while displaying 151 dpi previews of them. A reviewer reads the
manuscript, not the upload folder.

What matters for an embedded figure is its EFFECTIVE resolution: the pixels it
carries divided by the width it is displayed at. That is what a reader sees and
what a printer gets.

    python3 check_embedded_artwork.py [file.docx ...]

Exit status is 1 if any figure falls below the minimum.
"""
import io
import os
import sys

from PIL import Image
from docx import Document

MIN_DPI = 300.0
DEFAULTS = ["/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5.docx",
            "/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5_typeset.docx"]


def check(path):
    d = Document(path)
    sec = d.sections[0]
    text_w = (sec.page_width - sec.left_margin - sec.right_margin) / 914400
    rows, widths = [], set()
    for i, sh in enumerate(d.inline_shapes, 1):
        try:
            rid = sh._inline.graphic.graphicData.pic.blipFill.blip.embed
            blob = d.part.related_parts[rid].blob
            with Image.open(io.BytesIO(blob)) as im:
                px = im.size[0]
        except Exception as exc:
            rows.append((i, None, None, f"unreadable: {exc}"))
            continue
        disp = sh.width / 914400
        widths.add(round(disp, 2))
        eff = px / disp if disp else 0.0
        bad = eff < MIN_DPI - 1
        rows.append((i, px, eff, "below minimum" if bad else ""))
    bad = [r for r in rows if r[3]]
    print(f"=== {os.path.basename(path)} — {len(rows)} figure(s), "
          f"text width {text_w:.2f} in ===")
    #  a uniform display width is the house style; report it either way
    if len(widths) == 1:
        print(f"  [ok  ] every figure displayed at {next(iter(widths))} in "
              f"({next(iter(widths))*25.4:.0f} mm)")
    else:
        print(f"  [warn] {len(widths)} different display widths: "
              f"{sorted(widths)} — figures should share a standard width")
    for i, px, eff, why in bad:
        print(f"  [FAIL] Fig {i}: {px} px at display size = {eff:.0f} dpi "
              f"(minimum {MIN_DPI:.0f})")
    if not bad:
        lo = min(r[2] for r in rows if r[2])
        print(f"  [ok  ] every embedded figure is at least {MIN_DPI:.0f} dpi "
              f"at its display size (lowest {lo:.0f})")
    return len(bad) + (0 if len(widths) <= 1 else 1)


def main(argv):
    files = argv or DEFAULTS
    bad = 0
    for f in files:
        if not os.path.exists(f):
            print(f"  [skip] {f}: not found")
            continue
        bad += check(f)
        print()
    print("Embedded artwork meets the specification." if not bad
          else f"{bad} embedded-artwork problem(s).")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
