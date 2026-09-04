#!/usr/bin/env python3
# =============================================================================
#  check_slides.py — audit a deck for the faults that are obvious on the screen
#  and invisible in the file.
# -----------------------------------------------------------------------------
#  Three things go wrong when figures are placed programmatically, and none of
#  them raises an error:
#
#    EMPTY FRAME    a panel shape that was drawn to hold a picture, with no
#                   picture in it — an empty box on the slide, often still
#                   carrying the caption of the figure that used to be there.
#    TINY FIGURE    a detailed figure scaled into a frame built for a thumbnail,
#                   unreadable at presentation size.
#    OVERLAP        a text box sitting on top of a picture, or two text boxes on
#                   top of each other.
#
#  Everything here is measured from the shape geometry, so the report says which
#  slide and by how much rather than "looks wrong".
#
#      python3 check_slides.py [deck.pptx]
#
#  Exit status is 1 if anything is flagged.
# =============================================================================
import os
import sys

from pptx import Presentation
from pptx.util import Emu

#  a figure smaller than this is not readable from the back of a room
MIN_AREA_SQIN = 3.0
#  a container this big that holds nothing reads as an empty box
MIN_EMPTY_PANEL_SQIN = 2.0
#  ignore incidental touching; flag a real covering
#  two text boxes that merely abut share a sliver of each other's padding; only a
#  real covering is a fault
MIN_OVERLAP_FRAC = 0.25


def box(sh):
    return (Emu(sh.left).inches, Emu(sh.top).inches,
            Emu(sh.width).inches, Emu(sh.height).inches)


def overlap_frac(a, b):
    """Fraction of the SMALLER box that the two boxes share."""
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    x0, x1 = max(ax, bx), min(ax + aw, bx + bw)
    y0, y1 = max(ay, by), min(ay + ah, by + bh)
    if x1 <= x0 or y1 <= y0:
        return 0.0
    inter = (x1 - x0) * (y1 - y0)
    return inter / max(min(aw * ah, bw * bh), 1e-9)


def audit(path):
    prs = Presentation(path)
    print(f"=== {os.path.basename(path)} — {len(prs.slides._sldIdLst)} slides ===")
    n_fail = 0

    for i, slide in enumerate(prs.slides, 1):
        pics = [sh for sh in slide.shapes if sh.__class__.__name__ == "Picture"]
        texts = [sh for sh in slide.shapes
                 if sh.has_text_frame and sh.text_frame.text.strip()]
        #  container shapes: no text, no picture, but a real area — these are the
        #  frames the design draws behind a figure
        panels = [sh for sh in slide.shapes
                  if sh.__class__.__name__ != "Picture"
                  and (not sh.has_text_frame or not sh.text_frame.text.strip())]

        issues = []

        for p in pics:
            _l, _t, w, h = box(p)
            if w * h < MIN_AREA_SQIN:
                issues.append(f"figure only {w:.2f}x{h:.2f} in "
                              f"({w*h:.1f} sq in) — too small to read")

        for pan in panels:
            l, t, w, h = box(pan)
            if w * h < MIN_EMPTY_PANEL_SQIN:
                continue
            if w > 12.0 or h < 0.35:            # rules, dividers, backgrounds
                continue
            #  A panel is only EMPTY if nothing sits in it. Most of these shapes are
            #  the design's styled card behind a block of commentary, so a text box
            #  on top of one means it is doing its job, not that a figure is missing.
            holds_pic = any(overlap_frac(box(pan), box(p)) > 0.30 for p in pics)
            holds_txt = any(overlap_frac(box(pan), box(t)) > 0.25 for t in texts)
            if not (holds_pic or holds_txt):
                issues.append(f"empty {w:.2f}x{h:.2f} in panel at "
                              f"({l:.2f},{t:.2f}) — a frame with nothing in it")

        for t1 in texts:
            for p in pics:
                f = overlap_frac(box(t1), box(p))
                if f >= MIN_OVERLAP_FRAC:
                    issues.append(f"text {t1.text_frame.text.strip()[:34]!r} covers "
                                  f"{f*100:.0f} % of a figure")
        for a in range(len(texts)):
            for b in range(a + 1, len(texts)):
                f = overlap_frac(box(texts[a]), box(texts[b]))
                if f >= MIN_OVERLAP_FRAC:
                    issues.append(
                        f"text {texts[a].text_frame.text.strip()[:26]!r} overlaps "
                        f"{texts[b].text_frame.text.strip()[:26]!r} ({f*100:.0f} %)")

        if issues:
            print(f"  slide {i:2d}  ({len(pics)} figure(s))")
            for msg in issues:
                print(f"        - {msg}")
                n_fail += 1

    total_pics = sum(1 for s in prs.slides for sh in s.shapes
                     if sh.__class__.__name__ == "Picture")
    print()
    print(f"{total_pics} figure(s) across the deck, {n_fail} issue(s)")
    return 1 if n_fail else 0


def main(argv):
    path = argv[0] if argv else "/mnt/c/Users/user/Desktop/slides3.pptx"
    if not os.path.exists(path):
        print(f"not found: {path}")
        return 2
    return audit(path)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
