#!/usr/bin/env python3
"""Run the as-operated case at the LOWER BOUND of the published subcooling exponent.

The base case uses growth_exp_n = 1.0, the leading term of a power series in subcooling.
Bhattacharjee et al. (2021), Chem. Eng. Sci. 234:116417, report that n = 1.5-2.5 fits
experimental film-growth data. n = 1.0 is therefore defensible as a FORM and below the
range as a VALUE, and a verdict resting on it is a verdict resting on the most favourable
end of an unfitted parameter.

This runs the same duty at n = 1.5 -- the lower bound of that range -- so the paper can
report the verdict across the interval rather than at a point. It writes to its own
folder and does not touch the base outputs.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import solver  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "outputs_bounding_n15")


def main():
    import run_case_study10 as cs
    #  the SAME builder and the SAME arguments the base case uses, so the only difference
    #  between this run and outputs_steady is the exponent
    c = cs.build_case("bounding — as-operated at growth_exp_n = 1.5", "asoperated", 48.0)
    c.kinetics.growth_exp_n = 1.5
    os.makedirs(OUT, exist_ok=True)
    print("[bounding] as-operated duty at growth_exp_n = 1.5 "
          "(lower bound of Bhattacharjee et al. 2021)", flush=True)
    sv = solver.TransientSHCT(c)
    sv.run(verbose=True)
    eng = sv.engineering()
    solver.make_charts(sv, eng, OUT)
    with open(os.path.join(OUT, "key_metrics.json"), "w") as fh:
        solver.dump_json(eng, fh)
    print("[bounding] Phi_SH %.3g  P_plug %.3g  deposit %.3g mm  MEG %.3g wt%%"
          % (eng.get("max_Phi_SH", float("nan")), eng.get("P_plug", float("nan")),
             eng.get("peak_deposit_mm", float("nan")), eng.get("MEG_wt_pct", float("nan"))),
          flush=True)
    print("[bounding] ->", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
