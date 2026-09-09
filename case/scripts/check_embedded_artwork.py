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

from docx import Document
from PIL import Image

MIN_DPI = 300.0
DEFAULTS = ["/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5.docx",
            "/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5_typeset.docx"]


def _emu(v, default=0):
    """A length in EMU, or `default` when python-docx reports None.

    Section page size and margins are optional in the OOXML: a section that inherits
    them returns None, and `sec.page_width - sec.left_margin` then raises TypeError
    before this script can report anything at all. The same is true of an inline
    shape's width. Fall back rather than crash, and say which default was assumed.
    """
    return int(default) if v is None else int(v)


LETTER_W_EMU = 7772400          # 8.5 in — the assumed page width when a section omits one


def check(path):
    d = Document(path)
    sec = d.sections[0]
    if sec.page_width is None:
        print(f"  [note] {os.path.basename(path)}: section sets no page width; "
              f"assuming US Letter (8.5 in) for the text-width figure")
    text_w = (_emu(sec.page_width, LETTER_W_EMU)
              - _emu(sec.left_margin) - _emu(sec.right_margin)) / 914400
    rows: list[tuple[int, "int | None", "float | None", str]] = []
    widths: set[float] = set()
    for i, sh in enumerate(d.inline_shapes, 1):
        try:
            rid = sh._inline.graphic.graphicData.pic.blipFill.blip.embed
            blob = d.part.related_parts[rid].blob
            with Image.open(io.BytesIO(blob)) as im:
                px = im.size[0]
        except Exception as exc:
            rows.append((i, None, None, f"unreadable: {exc}"))
            continue
        disp = _emu(sh.width) / 914400
        widths.add(round(disp, 2))
        eff = px / disp if disp else 0.0
        bad = eff < MIN_DPI - 1
        rows.append((i, px, eff, "below minimum" if bad else ""))
    bad = [r for r in rows if r[3]]
    print(f"=== {os.path.basename(path)} — {len(rows)} figure(s), "
          f"text width {text_w:.2f} in ===")
    #  a uniform display width is the house style; report it either way
    if not widths:
        pass                          # nothing embedded: the width check has no subject
    elif len(widths) == 1:
        print(f"  [ok  ] every figure displayed at {next(iter(widths))} in "
              f"({next(iter(widths))*25.4:.0f} mm)")
    else:
        print(f"  [warn] {len(widths)} different display widths: "
              f"{sorted(widths)} — figures should share a standard width")
    for i, px, eff, why in bad:
        #  an image that could not be opened is recorded with px = eff = None, and
        #  formatting None with "%.0f" raised TypeError -- the checker crashed on the
        #  one case it exists to report. Print the recorded reason instead.
        if px is None or eff is None:
            print(f"  [FAIL] Fig {i}: {why}")
        else:
            print(f"  [FAIL] Fig {i}: {px} px at display size = {eff:.0f} dpi "
                  f"(minimum {MIN_DPI:.0f}) — {why}")
    if not bad:
        #  min() over an empty generator raises; a document with no inline shapes is
        #  a legitimate (if useless) input, so say so rather than fail.
        eff_vals = [r[2] for r in rows if r[2]]
        if eff_vals:
            print(f"  [ok  ] every embedded figure is at least {MIN_DPI:.0f} dpi "
                  f"at its display size (lowest {min(eff_vals):.0f})")
        else:
            print("  [ok  ] no embedded figures to check")
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
