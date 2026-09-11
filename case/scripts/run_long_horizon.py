#!/usr/bin/env python3
"""Run the as-operated duty long enough to find out whether "does not plug" is a
PHYSICAL result or a WINDOW result.

The case study reports 48 h. At growth_exp_n = 1.5 that run ends with Phi_SH = 1.66,
which is super-critical against Phi_crit = 1.08, and P_plug = 0 -- read naively, "the
coupling number says critical and the line does not plug anyway". That reading is wrong,
and the arithmetic says why:

    Phi_SH = delta_eq / delta_ref, delta_ref = 21.2 mm  ->  delta_eq = 35.19 mm
    consolidation binds at restr > consol_restriction = 0.18, i.e. delta > 22.90 mm
    deposit actually reached in 48 h                    =  8.81 mm

So the EQUILIBRIUM deposit is 54 % above the consolidation threshold, and 48 h reaches
only a quarter of that equilibrium. The line is not sitting below the threshold; it is
still climbing towards a value above it. "Does not plug" is a statement about the
simulated window, not about the duty -- and a real tie-back runs for months.

This script runs the same duty out to a horizon where the answer is not window-limited,
and reports WHEN the deposit crosses the consolidation restriction and whether the line
then goes on to block. If it never crosses, that is the physical result and it is worth
having; if it does, the 48 h headline needs its window stated beside it.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import solver  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "outputs_long_horizon")


def main():
    import run_case_study10 as cs
    hours = float(sys.argv[1]) if len(sys.argv) > 1 else 480.0
    c = cs.build_case("long horizon — as-operated", "asoperated", hours)
    os.makedirs(OUT, exist_ok=True)
    k = c.kinetics
    print("[long] as-operated duty to %.0f h at growth_exp_n = %.2f" % (hours, k.growth_exp_n),
          flush=True)
    sv = solver.TransientSHCT(c)
    sv.run(verbose=True)
    eng = sv.engineering()
    r = sv.results

    D = c.pipeline.diameter_m
    #  snap_phi is hydrate volume fraction; the deposit itself is r["delta"] (final state),
    #  so take the restriction history from the monitor series, which carries delta.
    ts_t = np.asarray(r.get("ts_t", np.empty(0)), float)
    ts_d = np.asarray(r.get("ts", {}).get("delta", []), float)
    restr = 2.0 * ts_d / D if ts_d.size else np.empty(0)
    thr = k.consol_restriction
    cross = None
    if restr.size and np.any(restr > thr):
        cross = float(ts_t[int(np.argmax(restr > thr))])

    rec = {"hours": hours, "growth_exp_n": k.growth_exp_n,
           "consol_restriction": thr, "plug_restriction_trip": c.numerics.plug_restriction_trip,
           "max_Phi_SH": eng.get("max_Phi_SH"), "P_plug": eng.get("P_plug"),
           "peak_deposit_mm": eng.get("peak_deposit_mm"),
           "final_monitor_restriction": float(restr[-1]) if restr.size else None,
           "max_monitor_restriction": float(restr.max()) if restr.size else None,
           "hours_to_cross_consol_restriction": cross,
           "time_to_plug_P50_h": eng.get("time_to_plug_P50_h")}
    with open(os.path.join(OUT, "long_horizon.json"), "w") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    with open(os.path.join(OUT, "key_metrics.json"), "w") as fh:
        solver.dump_json(eng, fh)

    print("\n" + "=" * 72)
    print(" DOES 'DOES NOT PLUG' SURVIVE A LONGER WINDOW?")
    print("=" * 72)
    print("  horizon                    %8.0f h" % hours)
    print("  max Phi_SH                 %8.4g   (Phi_crit = 1.08)" % (eng.get("max_Phi_SH") or float("nan")))
    print("  peak wall deposit          %8.2f mm" % (eng.get("peak_deposit_mm") or float("nan")))
    if restr.size:
        print("  monitor restriction max    %8.3f   (consolidates above %.2f, blocks at %.2f)"
              % (restr.max(), thr, c.numerics.plug_restriction_trip))
    print("  crosses consolidation at   %8s h" % ("%.0f" % cross if cross else "never"))
    print("  P_plug                     %8.3f" % (eng.get("P_plug") or 0.0))
    ttp = eng.get("time_to_plug_P50_h")
    print("  time to plug (P50)         %8s h" % ("%.0f" % ttp if ttp and ttp == ttp else "n/a"))
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
