#!/usr/bin/env python3
"""Run the as-operated case at the BOUNDS of the published subcooling exponent.

The base case now ships growth_exp_n = 1.5 -- Mori (2001) dT^(3/2), inside the 1.5-2.5
range Bhattacharjee et al. (2021), Chem. Eng. Sci. 234:116417 report for experimental
film-growth data. It used to ship 1.0, the leading Taylor term, which is defensible as a
FORM and below the range as a VALUE; the headline then rested on the most favourable end
of an unfitted parameter, so the default was moved.

That move makes THIS script's old job obsolete: it ran n = 1.5 as the bound, and 1.5 is
now the base. The bounds either side of the shipped value are what a reader needs instead:

    n = 1.0   the old default and the leading-order form -- the OPTIMISTIC bound
    n = 2.0   inside the published range, above the shipped floor -- the PESSIMISTIC bound

Each writes its own folder and none touches the base outputs.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import solver  # noqa: E402

BOUNDS = ((1.0, "outputs_bounding_n10", "optimistic — the old default, leading-order form"),
          (2.0, "outputs_bounding_n20", "pessimistic — inside the published range"))


def run_one(n, folder, why):
    import run_case_study10 as cs
    #  the SAME builder and the SAME arguments the base case uses, so the only difference
    #  between this run and outputs_steady is the exponent
    out = os.path.join(os.path.dirname(HERE), folder)
    c = cs.build_case("bounding — as-operated at growth_exp_n = %.1f" % n, "asoperated", 48.0)
    c.kinetics.growth_exp_n = n
    os.makedirs(out, exist_ok=True)
    print("[bounding] as-operated duty at growth_exp_n = %.1f  (%s)" % (n, why), flush=True)
    sv = solver.TransientSHCT(c)
    sv.run(verbose=True)
    eng = sv.engineering()
    solver.make_charts(sv, eng, out)
    with open(os.path.join(out, "key_metrics.json"), "w") as fh:
        solver.dump_json(eng, fh)
    print("[bounding] n=%.1f  Phi_SH %.3g  P_plug %.3g  deposit %.3g mm  MEG %.3g wt%%"
          % (n, eng.get("max_Phi_SH", float("nan")), eng.get("P_plug", float("nan")),
             eng.get("peak_deposit_mm", float("nan")), eng.get("MEG_wt_pct", float("nan"))),
          flush=True)
    return eng


def main():
    for n, folder, why in BOUNDS:
        run_one(n, folder, why)
        print("[bounding] ->", os.path.join(os.path.dirname(HERE), folder), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
