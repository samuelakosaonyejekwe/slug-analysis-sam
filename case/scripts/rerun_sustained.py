#!/usr/bin/env python3
"""Re-run the three case-study scenarios and record the SUSTAINED coupling number.

Adds nothing to the physics: the solver is deterministic on a fixed seed, so the fields
and every existing figure are unchanged. Only the new sustained-Phi_SH diagnostics are
computed, and they are merged into each scenario's key_metrics.json / summary.json.
Figures are deliberately NOT regenerated, so nothing the manuscript already cites moves.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
import run_case_study10 as R

import solver

NEW = ["sustained_Phi_SH", "sustained_Phi_SH_hotspot_km", "sustained_supercritical_km",
       "final_Phi_SH", "Phi_SH_supercritical_time_frac", "Phi_SH_peak_time_h", "max_Phi_SH",
       "coupled_hotspot_km"]

SCEN = [("outputs_steady",    "asoperated", 48.0, " — normal production (as-operated)"),
        ("outputs_shutin",    "shutin",     24.0, " — unplanned shut-in"),
        ("outputs_mitigated", "mitigated",  48.0, " — engineered fix (insulation + MEG)")]
BASE = "Deepwater Medium-Crude-Oil Subsea Tie-back (32 km, 10.75-in)"
OUTROOT = os.path.join(os.path.dirname(HERE))


def main():
    """Re-run the three scenarios and merge the sustained-Phi_SH diagnostics.

    THIS RUNS THE SOLVER AND OVERWRITES TRACKED OUTPUTS, so it lives behind a
    __main__ guard. It used to sit at module level, which meant that merely
    IMPORTING this file -- an editor's auto-import, a documentation tool, a
    `python -c "import rerun_sustained"` smoke test -- silently launched three
    full transient runs and rewrote key_metrics.json, summary.json and
    sustained_phiSH_profile.csv in all three output folders. A module that
    rewrites the repository's published results as a side effect of being
    imported is a trap regardless of whether the values happen to come back the
    same.
    """
    for folder, variant, tend, suffix in SCEN:
        t0 = time.time()
        print(f"\n=== {folder} ({variant}) ===", flush=True)
        case = R.build_case(BASE + suffix, variant, tend)
        sv = solver.TransientSHCT(case)
        sv.run(verbose=False)
        eng = sv.engineering()

        out = os.path.join(OUTROOT, folder)
        for name in ("key_metrics.json", "summary.json"):
            path = os.path.join(out, name)
            with open(path) as fh:
                j = json.load(fh)
            #  Sanity: the re-run must reproduce what is already on disk. This was an
            #  `assert`, which `python -O` strips -- and this script OVERWRITES tracked
            #  output files, so under -O it would have rewritten them from a run that no
            #  longer reproduced them, with the guard silently absent. A guard on a
            #  destructive step has to be a real check.
            for k in ("max_subcooling_C", "peak_deposit_mm", "time_to_plug_P50_h", "MEG_wt_pct"):
                old, new = j.get(k), eng.get(k)
                if isinstance(old, (int, float)) and isinstance(new, (int, float)):
                    if abs(old - new) > 1e-6 * max(1.0, abs(old)):
                        raise SystemExit(
                            f"{folder}/{name}: {k} moved {old} -> {new}; the re-run is NOT "
                            f"deterministic, so these outputs are not being rewritten")
            for k in NEW:
                v = eng.get(k)
                j[k] = float(v) if isinstance(v, (int, float, np.floating)) else v
            with open(path, "w") as fh:
                #  through the solver's strict writer: sustained_Phi_SH is NaN on a
                #  scenario that never forms hydrate (the mitigated line), and json.dump
                #  writes NaN as a bare literal the JSON standard does not permit --
                #  key_metrics.json and summary.json were shipping unparseable by any
                #  strict reader. json_safe turns an undefined value into null.
                solver.dump_json(j, fh)

        # the sustained field itself, for the manuscript
        snapP = np.asarray(sv.results["snap_PhiSH"], float)
        x_km = sv.x / 1000.0
        sust = np.nanmedian(snapP, axis=0)
        forming = np.nanmax(sv.results["max_Tsub"], axis=1) > 0.0
        sust = np.where(forming, sust, np.nan)
        solver._save_csv(os.path.join(out, "sustained_phiSH_profile.csv"),
                         ["x_km", "hydrate_forming", "Phi_SH_sustained", "Phi_SH_running_max",
                          "Phi_SH_final"],
                         [[f"{a:.4f}", int(m), f"{b:.6g}", f"{c:.6g}", f"{d:.6g}"]
                          for a, m, b, c, d in zip(x_km, forming, sust,
                                                   np.nanmedian(sv.results["max_PhiSH"], 1), snapP[-1])],
                         all_str=True)
        print(f"  running max {eng['max_Phi_SH']:.6f} (peaks at t={eng['Phi_SH_peak_time_h']:.3f} h)")
        print(f"  SUSTAINED   {eng['sustained_Phi_SH']:.4g} at {eng['sustained_Phi_SH_hotspot_km']:.2f} km"
              f"  | final {eng['final_Phi_SH']:.4g}"
              f"  | supercritical {100*eng['Phi_SH_supercritical_time_frac']:.0f}% of the window")
        print(f"  ({time.time()-t0:.0f} s)", flush=True)
    print("\nDONE")


if __name__ == "__main__":
    main()

