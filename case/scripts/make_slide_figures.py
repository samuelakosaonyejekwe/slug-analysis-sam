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

#  both must precede the solver/style import
os.environ["SHCT_FIG_FONTSCALE"] = "1.8"
os.environ.setdefault("SHCT_FIG_DPI", "200")

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
SCENARIOS = [("outputs_slides", "asoperated", 48.0),
             ("outputs_slides_shutin", "shutin", 24.0),
             ("outputs_slides_mitigated", "mitigated", 48.0)]


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
        try:
            import shct_crosssection
            shct_crosssection.crosssection_outputs(sv, outdir)
        except Exception as exc:
            print(f"  [{dirname}] cross-section skipped: {exc}", flush=True)
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
    print(f"[slide-figs] regenerating {len(todo)} scenario(s) at "
          f"fontscale 1.8 / {os.environ['SHCT_FIG_DPI']} dpi", flush=True)
    with ProcessPoolExecutor(max_workers=min(3, len(todo))) as ex:
        list(ex.map(one, todo))
    #  the cross-scenario comparison chart, from metrics already on disk
    try:
        import json
        base = json.load(open(os.path.join(CASE, "outputs_steady", "summary.json")))
        mit = json.load(open(os.path.join(CASE, "outputs_mitigated", "summary.json")))
        R.mitigation_chart(base, mit, os.path.join(CASE, "outputs_slides"))
        print("  mitigation comparison rebuilt", flush=True)
    except Exception as exc:
        print(f"  mitigation chart skipped: {exc}", flush=True)
    print("[slide-figs] COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
