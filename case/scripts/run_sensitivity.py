#!/usr/bin/env python3
# =============================================================================
#  run_sensitivity.py — one-at-a-time sensitivity of the coupled predictions to
#  the three assumed kinetic/coupling constants.
#
#  Phi_SH = C_phi * k_g,wall * a_i * dT_sub,wall**n / f_slug, and k_g,wall is set
#  by kg0.  None of C_phi, n or kg0 is fitted to data in this study, so the
#  absolute magnitudes of Phi_SH, of the time-to-plug and of the required MEG
#  dose inherit whatever uncertainty those three constants carry.  This script
#  quantifies that inheritance instead of leaving it to the reader.
#
#  Ranges:
#    kg0   x0.2 .. x5.0   — the solver's own documented calibration bounds
#                           (solver.py CALIB_BOUNDS: "kg0": (0.2, 5.0))
#    n     1.0 .. 2.0      — unity for heat/mass-transfer-controlled growth,
#                           quadratic also reported in the hydrate literature
#    C_phi 500 .. 4500     — +/-3x about the assumed 1500 (no measured value exists)
#
#  Reduced but SELF-CONSISTENT fidelity: every run here, the baseline included,
#  uses the same ensemble size and simulated window, so the movements are
#  comparable with one another.  They are NOT directly comparable with the
#  headline case-study numbers, which use n_ensemble=12 / t_end=48 h.
# =============================================================================
import os, sys, json, time, csv
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import OUT
import run_case_study10 as R
import solver

N_ENSEMBLE = 6
N_CELLS    = 70
T_END_H    = 24.0

METRICS = ["max_Phi_SH", "max_Phi_SH_uncapped", "max_Psi_kinetic_ratio",
           "Phi_SH_gate_saturated_frac", "P_plug", "time_to_plug_P50_h", "time_to_plug_P10_h",
           "time_to_plug_P90_h", "MEG_wt_pct", "under_inhibited_km",
           "peak_deposit_mm", "max_subcooling_C", "cooldown_to_hydrate_h"]

BASE = {"kg0_mult": 1.0, "growth_exp_n": 1.0, "C_phi": 1500.0}


def make_runs():
    runs = [("baseline", dict(BASE))]
    for m in (0.2, 0.5, 2.0, 5.0):
        d = dict(BASE); d["kg0_mult"] = m
        runs.append((f"kg0_x{m:g}", d))
    for n in (1.25, 1.5, 1.75, 2.0):
        d = dict(BASE); d["growth_exp_n"] = n
        runs.append((f"n_{n:g}", d))
    for c in (500.0, 1000.0, 3000.0, 4500.0):
        d = dict(BASE); d["C_phi"] = c
        runs.append((f"C_{c:g}", d))
    return runs


def one(job):
    label, p = job
    t0 = time.time()
    case = R.build_case(f"sensitivity {label}", "asoperated",
                        T_END_H, n_ensemble=N_ENSEMBLE, n_cells=N_CELLS)
    k = case.kinetics
    k.kg0 = 6.0e-7 * p["kg0_mult"]
    k.growth_exp_n = p["growth_exp_n"]
    k.C_phi = p["C_phi"]
    sv = solver.TransientSHCT(case)
    sv.run()
    eng = sv.engineering()
    row = {"label": label, "kg0_mult": p["kg0_mult"],
           "growth_exp_n": p["growth_exp_n"], "C_phi": p["C_phi"],
           "runtime_s": round(time.time() - t0, 1)}
    row.update({m: eng.get(m) for m in METRICS})
    print(f"  done {label:12s} ({row['runtime_s']:.0f}s)  "
          f"Phi_max={row['max_Phi_SH']}  P50={row['time_to_plug_P50_h']}  "
          f"MEG={row['MEG_wt_pct']}", flush=True)
    return row


def main():
    outdir = OUT["steady"]
    os.makedirs(outdir, exist_ok=True)
    jobs = make_runs()
    print(f"[sensitivity] {len(jobs)} runs, n_ensemble={N_ENSEMBLE}, "
          f"t_end={T_END_H} h", flush=True)
    workers = int(os.environ.get("SENS_WORKERS", "2"))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(one, jobs))

    base = next(r for r in rows if r["label"] == "baseline")
    for r in rows:
        for m in METRICS:
            b, v = base.get(m), r.get(m)
            r[m + "_rel"] = (round(v / b, 4)
                             if (isinstance(b, (int, float)) and
                                 isinstance(v, (int, float)) and b) else None)

    cols = (["label", "kg0_mult", "growth_exp_n", "C_phi", "runtime_s"] +
            METRICS + [m + "_rel" for m in METRICS])
    csv_path = os.path.join(outdir, "sensitivity_phiSH.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})
    with open(os.path.join(outdir, "sensitivity_phiSH.json"), "w") as fh:
        json.dump({"settings": {"n_ensemble": N_ENSEMBLE, "n_cells": N_CELLS,
                                "t_end_h": T_END_H, "scenario": "asoperated"},
                   "baseline": BASE, "rows": rows}, fh, indent=2, default=str)
    print(f"[sensitivity] -> {csv_path}", flush=True)
    return rows


if __name__ == "__main__":
    main()
