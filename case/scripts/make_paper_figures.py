#!/usr/bin/env python3
"""Regenerate the case-study figures WITHOUT the software-flavoured chart titles,
for use in a journal manuscript where the caption carries the description.

Writes to case/outputs_paper_<scenario>/ so the titled figures the case-study
report uses are left untouched. Run with SHCT_FIG_TITLES=0 (this script sets it
before importing the solver).
"""
import os
import sys

os.environ["SHCT_FIG_TITLES"] = "0"          # must precede the solver import
os.environ.setdefault("SHCT_FIG_DPI", "320")  # IJMF: >=300 dpi, >=1063 px single column

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib                             # noqa: E402
matplotlib.use("Agg")

from concurrent.futures import ProcessPoolExecutor   # noqa: E402

from _paths import CASE                              # noqa: E402
import run_case_study10 as R                         # noqa: E402
import solver                                        # noqa: E402

SCENARIOS = [("steady", "asoperated", 48.0),
             ("shutin", "shutin", 24.0),
             ("mitigated", "mitigated", 48.0)]


def one(job):
    name, variant, t_end = job
    outdir = os.path.join(CASE, f"outputs_paper_{name}")
    os.makedirs(outdir, exist_ok=True)
    case = R.build_case(f"paper figures {name}", variant, t_end,
                        n_ensemble=12, n_cells=70)
    sv = solver.TransientSHCT(case)
    sv.run()
    eng = sv.engineering()
    solver.make_charts(sv, eng, outdir)
    #  the driver's bespoke charts (slug prediction, riser screen, mitigation)
    if name == "steady":
        try:
            R.slug_chart(sv, outdir)
            R.riser_chart(sv, outdir)
        except Exception as exc:
            print(f"  [{name}] driver charts skipped: {exc}", flush=True)
    #  the steady case also carries the cross-section and compositional figures
    if name == "steady":
        try:
            import shct_crosssection
            shct_crosssection.crosssection_outputs(sv, outdir)
        except Exception as exc:
            print(f"  [{name}] cross-section skipped: {exc}", flush=True)
        try:
            import shct_compositional
            shct_compositional.compositional_report(sv, outdir)
        except Exception as exc:
            print(f"  [{name}] compositional skipped: {exc}", flush=True)
    print(f"  done {name}: {len([f for f in os.listdir(outdir) if f.endswith('.png')])} "
          f"figures -> {outdir}", flush=True)
    return outdir


def main(argv=None):
    #  Re-running every scenario costs ~13 min each. Accept a filter so that a run
    #  which only needs, say, the steady figures does not redo the other two:
    #      python3 make_paper_figures.py steady
    argv = argv if argv is not None else sys.argv[1:]
    todo = [s for s in SCENARIOS if not argv or s[0] in argv]
    if not todo:
        print(f"no scenario matched {argv}; known: "
              f"{', '.join(s[0] for s in SCENARIOS)}", flush=True)
        return 2
    print(f"[paper-figs] regenerating {len(todo)} scenario(s) without chart titles: "
          f"{', '.join(s[0] for s in todo)}", flush=True)
    with ProcessPoolExecutor(max_workers=min(3, len(todo))) as ex:
        list(ex.map(one, todo))
    #  mitigation comparison needs the as-operated and mitigated metrics together;
    #  they are already on disk, so no extra transient run is required
    out = os.path.join(CASE, "outputs_paper_steady")
    try:
        import json
        eng_base = json.load(open(os.path.join(CASE, "outputs_steady", "summary.json")))
        eng_mit = json.load(open(os.path.join(CASE, "outputs_mitigated", "summary.json")))
        R.mitigation_chart(eng_base, eng_mit, out)
        print("  mitigation comparison rebuilt", flush=True)
    except Exception as exc:
        print(f"  mitigation chart skipped: {exc}", flush=True)

    #  validation charts (fast, no transient run needed)
    out = os.path.join(CASE, "outputs_paper_steady")
    try:
        solver.validate_hydrate_curve(
            os.path.join(os.path.dirname(os.path.abspath(solver.__file__)),
                         "validation", "data", "hydrate_equilibrium_published.json"),
            outdir=out)
    except Exception as exc:
        print(f"  hydrate validation skipped: {exc}", flush=True)
    try:
        solver.validate_closures(outdir=out)
    except Exception as exc:
        print(f"  closure validation skipped: {exc}", flush=True)
    print("[paper-figs] complete", flush=True)


if __name__ == "__main__":
    sys.exit(main() or 0)
