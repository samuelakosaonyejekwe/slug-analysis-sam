#!/usr/bin/env python3
"""Regenerate the SLIDE-LEGIBLE twin of every case-study figure.

A figure drawn for a journal column is read at 3.5 in from arm's length; the same
figure projected on a slide is read from the back of a room, where its 8 pt tick
labels land nearer 4 pt. The deck therefore carries its own rendering of each
figure — same run, same numbers, larger type — and the deck builder prefers that
twin over the print figure whenever one exists.

Those twins were previously produced by hand, which meant that after a solver
change the deck could be refreshed from figures nobody could reproduce. This
script is that recipe, committed: same scenarios and same resolution as the case
study (12 realisations, 70 cells, seed 13), rendered at SHCT_FIG_FONTSCALE=1.8
and 200 dpi.

    python3 make_slide_figures.py            # all three scenarios
    python3 make_slide_figures.py steady     # just one

Writes case/outputs_slides{,_shutin,_mitigated}/ — the three names audit_deck.py
and check_slides.py already know as the scaled set.
"""
import os
import sys

#  Both must precede the solver/style import.
#
#  SIZE, not font size, is the lever. Type on a slide renders at
#      effective_pt = base_pt x (displayed_width / natural_width)
#  so a figure drawn 10 in wide and shown in a 3.3 in frame puts its 18 pt labels
#  on the wall at 5.9 pt however large the fonts were set. Shrinking the FIGURE
#  brings natural_width toward the frame width, and the type then arrives at
#  nearly the size it was set in.
#
#  One scale cannot serve the whole deck: a global factor simply divides every
#  figure's effective size by the same number, and the deck's frames run from
#  2.4 in to 8.0 in. So this script renders a SET at one scale, and is run once
#  per scale; fit_deck_figures.py then picks, per frame, the set that lands
#  closest to the target type size.
#
#      python3 make_slide_figures.py --scale 0.45
#
#  writes case/outputs_slides45{,_shutin,_mitigated}.
_SCALE = "1.0"
if "--scale" in sys.argv:
    _SCALE = sys.argv[sys.argv.index("--scale") + 1]
    del sys.argv[sys.argv.index("--scale"):sys.argv.index("--scale") + 2]
_TAG = "" if abs(float(_SCALE) - 1.0) < 1e-9 else f"{round(float(_SCALE) * 100):02d}"

os.environ["SHCT_FIG_FONTSCALE"] = "1.8"
os.environ["SHCT_FIG_SIZESCALE"] = _SCALE
#  a shrunk figure shown at or above its natural size is upscaled on the wall, so
#  the pixel budget has to rise as the scale falls to keep it from softening
os.environ.setdefault("SHCT_FIG_DPI", str(int(round(200 / max(float(_SCALE), 0.25)))))

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib                             # noqa: E402
matplotlib.use("Agg")

from concurrent.futures import ProcessPoolExecutor   # noqa: E402

from _paths import CASE                              # noqa: E402
import run_case_study10 as R                         # noqa: E402
import solver                                        # noqa: E402
import shct_spacetime                                # noqa: E402

#  (output directory, build_case variant, t_end_h) — the directory names are the
#  ones already wired into audit_deck.SCALED_DIRS and check_slides.
SCENARIOS = [(f"outputs_slides{_TAG}", "asoperated", 48.0),
             (f"outputs_slides{_TAG}_shutin", "shutin", 24.0),
             (f"outputs_slides{_TAG}_mitigated", "mitigated", 48.0)]


def one(job):
    dirname, variant, t_end = job
    outdir = os.path.join(CASE, dirname)
    os.makedirs(outdir, exist_ok=True)
    case = R.build_case(f"slide figures {variant}", variant, t_end,
                        n_ensemble=12, n_cells=70)
    sv = solver.TransientSHCT(case)
    sv.run()
    eng = sv.engineering()
    solver.make_charts(sv, eng, outdir)
    shct_spacetime.spacetime_outputs(sv, eng, outdir)
    if variant == "asoperated":
        for fn, label in ((R.slug_chart, "slug"), (R.riser_chart, "riser")):
            try:
                fn(sv, outdir)
            except Exception as exc:
                print(f"  [{dirname}] {label} chart skipped: {exc}", flush=True)
        #  The deck also carries the cross-section, compositional, 3-D and
        #  closure-validation figures. They come from other modules, so an earlier
        #  version of this script left them out — and a figure with no small
        #  variant simply cannot be made legible in a small frame, whatever the
        #  deck builder does. Sixteen of the deck's unreadable figures were of
        #  exactly that kind, so generate every family the deck uses.
        for label, call in (
            ("cross-section",
             lambda: __import__("shct_crosssection").crosssection_outputs(sv, outdir)),
            ("compositional PVT",
             lambda: __import__("shct_compositional").compositional_report(sv, outdir)),
            ("compositional transport",
             lambda: __import__("shct_compositional_sim").simulate_composition(sv, outdir)),
            ("3-D field",
             lambda: __import__("shct_threed").threed_outputs(sv, outdir)),
            ("closure validation",
             lambda: solver.validate_closures(
                 outdir=outdir,
                 datadir=os.path.join(os.path.dirname(CASE), "validation", "data"))),
            #  only the one verification figure the deck actually shows — running
            #  the whole suite per size variant would cost more than it returns
            ("water-faucet verification",
             lambda: __import__("shct_verification").check_water_faucet(outdir)),
        ):
            try:
                call()
            except Exception as exc:
                print(f"  [{dirname}] {label} skipped: {exc}", flush=True)
    n = len([f for f in os.listdir(outdir) if f.endswith(".png")])
    print(f"  done {dirname}: {n} figures", flush=True)
    return outdir


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    todo = [s for s in SCENARIOS
            if not argv or s[0] in argv or s[1] in argv
            or s[0].replace("outputs_slides", "steady").strip("_") in argv]
    if not todo:
        print(f"no scenario matched {argv}; known: "
              f"{', '.join(s[0] for s in SCENARIOS)}", flush=True)
        return 2
    print(f"[slide-figs] regenerating {len(todo)} scenario(s) at fontscale 1.8, "
          f"sizescale {_SCALE}, {os.environ['SHCT_FIG_DPI']} dpi", flush=True)
    with ProcessPoolExecutor(max_workers=min(3, len(todo))) as ex:
        list(ex.map(one, todo))
    #  the cross-scenario comparison chart, from metrics already on disk
    try:
        import json
        base = json.load(open(os.path.join(CASE, "outputs_steady", "summary.json")))
        mit = json.load(open(os.path.join(CASE, "outputs_mitigated", "summary.json")))
        R.mitigation_chart(base, mit, os.path.join(CASE, f"outputs_slides{_TAG}"))
        print("  mitigation comparison rebuilt", flush=True)
    except Exception as exc:
        print(f"  mitigation chart skipped: {exc}", flush=True)
    print("[slide-figs] COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
