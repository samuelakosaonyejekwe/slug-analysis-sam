#!/usr/bin/env python3
"""Replace each figure embedded in the manuscript with the current exported file.

python-docx embeds the file it is given and only sets the DISPLAY width, so an
embedded figure is a snapshot of whatever existed when it was inserted. Regenerating
the outputs updates Figure_NN.png on disk and leaves the manuscript showing the old
rendering at whatever resolution it had — which is how this manuscript ended up
displaying figures at 151 dpi that correlate poorly with their own source files.

Figure N in the manuscript is Figure_N.png: the captions match the export map in
export_paper_figures.py entry for entry, which is the authority here. Content
matching is NOT used — two coupling-number maps of different scenarios look alike
enough to fool it, and it was proposing to swap them.

Images are downscaled if that is needed to keep the file manageable, but never
below the resolution the journal requires at the display width.

    python3 reembed_figures.py [file.docx ...]
"""
import io
import os
import sys

from docx import Document
from PIL import Image

#  derived from this file, not hardcoded to one machine's home directory: every
#  other script in this folder locates the repository through _paths.py, and an
#  absolute /home/... path here breaks on any other checkout.
FIGS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "figures_paper")
MIN_DPI = 320.0             # comfortably over the journal's 300
MAX_PX = 3600               # beyond this the file grows for no visible gain
DEFAULTS = ["/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5.docx",
            "/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/paper5_typeset.docx"]


def main(argv):
    files = argv or DEFAULTS
    for path in files:
        if not os.path.exists(path):
            print(f"  [skip] {path}")
            continue
        d = Document(path)
        n_repl = n_same = 0
        for i, sh in enumerate(d.inline_shapes, 1):
            src = os.path.join(FIGS, f"Figure_{i}.png")
            if not os.path.exists(src):
                print(f"  Fig {i}: Figure_{i}.png missing — left as is")
                continue
            disp = sh.width / 914400
            need = int(round(MIN_DPI * disp))
            with Image.open(src) as im:
                w, h = im.size
                target = min(max(need, 0), MAX_PX)
                if w > target:
                    im2 = im.resize((target, int(round(h * target / w))),
                                    Image.LANCZOS)
                else:
                    im2 = im.copy()
                buf = io.BytesIO()
                im2.save(buf, format="PNG", dpi=(MIN_DPI, MIN_DPI), optimize=True)
            rid = sh._inline.graphic.graphicData.pic.blipFill.blip.embed
            part = d.part.related_parts[rid]
            before = part.blob
            new = buf.getvalue()
            if new == before:
                n_same += 1
                continue                     # identical bytes: nothing to write
            part._blob = new
            n_repl += 1
        if not n_repl:
            #  Nothing CHANGED, so do not rewrite the document: saving anyway bumps its
            #  mtime and the PDF-freshness check then calls the export stale against a
            #  manuscript that did not move.
            #
            #  This guard used to fire only when nothing MATCHED. Re-running the script on
            #  an already-current manuscript matched all 37 figures, re-encoded every one to
            #  the same bytes, and saved regardless — so a no-op run still bumped the mtime
            #  and still made paper5.pdf look stale. Byte-equality is the test that answers
            #  the question the guard was written for.
            print(f"  {os.path.basename(path)}: all {n_same} figure(s) already current — "
                  f"left untouched")
            continue
        d.save(path)
        mb = os.path.getsize(path) / 1e6
        print(f"  {os.path.basename(path)}: {n_repl} figure(s) re-embedded "
              f"({n_same} already identical, left alone), file now {mb:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
