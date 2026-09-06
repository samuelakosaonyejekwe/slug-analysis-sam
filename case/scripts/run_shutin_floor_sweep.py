#!/usr/bin/env python3
# =============================================================================
#  run_shutin_floor_sweep.py — does the slug-frequency floor actually bind?
# -----------------------------------------------------------------------------
#  run_sensitivity.py sweeps f_slug_floor_Hz on the AS-OPERATED line and finds
#  that nothing moves at all: every column is bit-identical across two decades.
#  That is not evidence that the floor is harmless. It is evidence that the floor
#  is never reached, because the as-operated line slugs from the first step at
#  frequencies of order 0.1 Hz, so the denominator of Phi_SH is always the
#  correlated slug frequency and never the guard beneath it.
#
#  The floor binds where the line is NOT slugging. This script runs the same
#  two-decade sweep on the SHUT-IN scenario, at the same reduced fidelity as
#  run_sensitivity.py so the two are comparable, and it is what Section 6.5 of
#  the manuscript reports: the peak Phi_SH and Psi go exactly as 1/f0 while the
#  time-to-plug, the inhibitor dose and the deposit do not move at all.
#
#      python3 run_shutin_floor_sweep.py
#
#  -> case/outputs_shutin/sensitivity_phiSH_floor.csv
# =============================================================================
import csv
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import OUT                       # noqa: E402
import run_case_study10 as R                 # noqa: E402
import solver                                # noqa: E402

N_ENSEMBLE, N_CELLS, T_END_H = 6, 70, 24.0   # as run_sensitivity.py
FLOORS = (1e-5, 3e-5, 1e-4, 3e-4, 1e-3)      # two decades about the 1e-4 default

METRICS = ["max_Phi_SH_uncapped", "max_Psi_kinetic_ratio", "sustained_Phi_SH",
           "sustained_supercritical_km", "Phi_SH_above_critical_frac",
           "time_to_plug_P50_h", "MEG_wt_pct", "peak_deposit_mm"]


def main():
    rows = []
    for ff in FLOORS:
        t0 = time.time()
        case = R.build_case(f"shutin floor {ff:g}", "shutin", T_END_H,
                            n_ensemble=N_ENSEMBLE, n_cells=N_CELLS)
        case.kinetics.f_slug_floor_Hz = ff
        sv = solver.TransientSHCT(case)
        sv.run(verbose=False)
        e = sv.engineering()
        row = {"f_slug_floor_Hz": ff, "runtime_s": round(time.time() - t0, 1)}
        row.update({m: e.get(m) for m in METRICS})
        rows.append(row)
        print(f"  f0={ff:.0e}  Phi_max={row['max_Phi_SH_uncapped']:12.2f}  "
              f"Psi={row['max_Psi_kinetic_ratio']:10.4f}  "
              f"P50={row['time_to_plug_P50_h']:.2f} h", flush=True)

    base = [r for r in rows if r["f_slug_floor_Hz"] == 1e-4][0]
    for r in rows:
        r["Phi_max_rel"] = round(r["max_Phi_SH_uncapped"] / base["max_Phi_SH_uncapped"], 4)
        r["f0_rel_inverse"] = round(1e-4 / r["f_slug_floor_Hz"], 4)

    out = os.path.join(OUT["shutin"], "sensitivity_phiSH_floor.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", out)
    print("\n  Phi_max/base against 1e-4/f0 — equal columns mean exact 1/f0 scaling:")
    for r in rows:
        print("   f0=%-8.0e  %8.4f   %8.4f" % (r["f_slug_floor_Hz"],
                                               r["Phi_max_rel"], r["f0_rel_inverse"]))


if __name__ == "__main__":
    main()
