#!/usr/bin/env python3
# =============================================================================
#  check_journal_artwork.py — verify the manuscript figure set against the
#  journal's artwork specification before anything is uploaded.
# -----------------------------------------------------------------------------
#  The target is the International Journal of Multiphase Flow (Elsevier). Its
#  Guide for Authors sets artwork requirements that a figure can fail silently:
#  a chart written at 150 dpi looks identical on screen and is rejected at
#  submission. This checks every numbered figure against them.
#
#      python3 check_journal_artwork.py [figures-dir]
#
#  Checks, from the Elsevier artwork guide:
#    * resolution   — >= 300 dpi for halftone and combination artwork
#    * width        — >= 1063 px, the pixel width of a 90 mm single column at
#                     300 dpi; a narrower figure will be upscaled and blur
#    * colour mode  — RGB (CMYK conversion is done by the publisher)
#    * format       — TIFF, EPS, PDF or a high-resolution PNG/JPEG
#    * naming       — Figure_1 ... Figure_N, numbered in citation order with no
#                     gaps, so the artwork items map onto the captions
#    * file size    — flagged over 10 MB, which the submission system throttles
#    * FRESHNESS    — every figure must be newer than what it was made from
#
#  The freshness check exists because the others cannot catch the worst failure.
#  A figure is spec-perfect at 320 dpi and correctly sized while showing results
#  from a solver revision that no longer exists, and this script passed exactly
#  that set for hours: outputs_paper_*/ had been rebuilt, but the numbered
#  Figure_N set assembled FROM it had not, so the artwork was current in
#  appearance and stale in content. Nothing about a PNG says which run produced
#  it. Modification times do, along the chain that actually produced the file:
#
#      solver.py -> outputs_<scenario>/ -> outputs_paper_<scenario>/ -> Figure_N
#
#  so a figure older than any of its upstreams is reported as stale. It is a
#  proxy -- touching a file defeats it -- but it catches the real mistake, which
#  is forgetting a step in a multi-stage rebuild rather than trying to cheat.
#
#  Exit status is 1 if any figure fails, so it can gate a submission pack build.
# =============================================================================
import os
import sys

try:
    from PIL import Image
except Exception:                                          # pragma: no cover
    Image = None

#  Where this file sits, so the upstream chain can be located when the checker is
#  pointed at a copy of the figures somewhere else (e.g. the submission folder).
_HERE = os.path.dirname(os.path.abspath(__file__))

MIN_DPI = 300
MIN_WIDTH_PX = 1063          # 90 mm single column at 300 dpi
MAX_MB = 10.0
OK_FORMATS = {"PNG", "TIFF", "JPEG", "EPS", "PDF"}


def check_one(path):
    """Return (list of failures, list of notes) for one figure."""
    fails, notes = [], []
    mb = os.path.getsize(path) / 1e6
    if mb > MAX_MB:
        fails.append(f"{mb:.1f} MB exceeds the {MAX_MB:.0f} MB upload limit")
    if Image is None:
        return fails, ["Pillow unavailable — only the file size was checked"]
    try:
        im = Image.open(path)
    except Exception as exc:
        return [f"unreadable: {exc}"], notes

    if im.format not in OK_FORMATS:
        fails.append(f"format {im.format} is not an accepted artwork type")

    dpi = im.info.get("dpi", (0, 0))
    dpi_x = float(dpi[0]) if dpi and dpi[0] else 0.0
    if dpi_x < MIN_DPI - 1:                # -1 absorbs 319.99 style rounding
        fails.append(f"{dpi_x:.0f} dpi is below the {MIN_DPI} dpi minimum")

    w, h = im.size
    if w < MIN_WIDTH_PX:
        fails.append(f"{w} px wide is below the {MIN_WIDTH_PX} px single-column "
                     f"minimum")

    if im.mode not in ("RGB", "RGBA", "L"):
        fails.append(f"colour mode {im.mode} — submit RGB")

    #  a figure far wider than a double column will be scaled down, wasting the
    #  resolution it was rendered at; worth knowing, not a failure
    if w > 4 * MIN_WIDTH_PX:
        notes.append(f"{w} px is very wide; it will be scaled down in typesetting")
    notes.append(f"{w}x{h} px, {dpi_x:.0f} dpi, {im.mode}, {mb:.1f} MB")
    return fails, notes


def _mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def freshness(nums, figdir):
    """Report figures older than the artefacts they were generated from.

    Returns (failures, skipped_reason). A missing upstream is NOT a failure --
    the checker must still be usable against a bare directory of figures -- it
    simply means freshness cannot be judged, and that is said out loud rather
    than passing silently.
    """
    sys.path.insert(0, _HERE)
    try:
        from export_paper_figures import FIGURE_MAP
        from _paths import CASE
    except Exception as exc:
        return [], f"upstream sources not importable ({type(exc).__name__})"

    solver_py = os.path.join(os.path.dirname(CASE), "solver.py")
    t_solver = _mtime(solver_py)
    fails = []

    for num, primary, fallback in FIGURE_MAP:
        if num not in nums:
            continue
        fig = os.path.join(figdir, f"Figure_{num}.png")
        t_fig = _mtime(fig)
        if t_fig is None:
            continue
        src = next((os.path.join(CASE, rel) for rel in (primary, fallback)
                    if rel and os.path.exists(os.path.join(CASE, rel))), None)
        if src is None:
            continue
        #  the figure must be at least as new as the intermediate it was copied from
        t_src = _mtime(src)
        if t_src is not None and t_fig < t_src - 1.0:
            fails.append((num, f"older than its source {os.path.relpath(src, CASE)} "
                               f"by {(t_src - t_fig) / 60.0:.0f} min — rerun "
                               f"export_paper_figures.py"))
            continue
        #  and the intermediate must be at least as new as the scenario data it
        #  was drawn from, or the whole branch is stale however fresh the copy is
        scen = os.path.basename(os.path.dirname(src)).replace("outputs_paper_", "outputs_")
        summary = os.path.join(CASE, scen, "summary.json")
        t_sum = _mtime(summary)
        if t_sum is not None and t_src < t_sum - 1.0:
            fails.append((num, f"source {os.path.relpath(src, CASE)} predates "
                               f"{scen}/summary.json by {(t_sum - t_src) / 60.0:.0f} min "
                               f"— rerun make_paper_figures.py"))
            continue
        if t_solver is not None and t_fig < t_solver - 1.0:
            fails.append((num, f"older than solver.py by "
                               f"{(t_solver - t_fig) / 60.0:.0f} min — the code changed "
                               f"after this figure was made"))
    return fails, None


def main(argv):
    figdir = argv[0] if argv else \
        "/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/figures"
    if not os.path.isdir(figdir):
        print(f"figures directory not found: {figdir}")
        return 2

    files = [f for f in os.listdir(figdir)
             if f.startswith("Figure_") and not f.startswith(".")]
    nums = []
    for f in files:
        try:
            nums.append(int(os.path.splitext(f)[0].split("_")[1]))
        except (IndexError, ValueError):
            print(f"  [FAIL] {f}: not named Figure_<n>")
    nums.sort()

    print(f"=== {figdir} ===")
    print(f"{len(files)} artwork file(s), Figure_1 .. Figure_{max(nums) if nums else 0}")
    gaps = [n for n in range(1, (max(nums) if nums else 0) + 1) if n not in nums]
    if gaps:
        print(f"  [FAIL] numbering gaps: {gaps} — artwork must be numbered in "
              f"citation order with no gaps")
    print()

    n_fail = 0
    for n in nums:
        path = os.path.join(figdir, f"Figure_{n}.png")
        if not os.path.exists(path):
            cand = [f for f in files if f.startswith(f"Figure_{n}.")]
            if not cand:
                continue
            path = os.path.join(figdir, cand[0])
        fails, notes = check_one(path)
        tag = "FAIL" if fails else "ok  "
        print(f"  [{tag}] Figure_{n:<2d}  {notes[-1] if notes else ''}")
        for f in fails:
            print(f"           - {f}")
            n_fail += 1

    stale, skipped = freshness(nums, figdir)
    print()
    if skipped:
        print(f"  [NOTE] freshness not checked — {skipped}")
    elif stale:
        print(f"  [FAIL] {len(stale)} figure(s) are STALE — they meet the artwork "
              f"specification but do not show the current results:")
        for num, why in stale:
            print(f"           Figure_{num}: {why}")
        n_fail += len(stale)
    else:
        print("  [ok  ] every figure is newer than the artefacts it was generated from")

    print()
    if n_fail or gaps:
        print(f"{n_fail} artwork problem(s) — do not upload until these are fixed")
        return 1
    print("Every figure meets the journal artwork specification.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
