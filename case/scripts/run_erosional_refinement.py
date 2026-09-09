#!/usr/bin/env python3
# =============================================================================
#  run_erosional_refinement.py — measure how the erosional statistics move under
#  grid refinement, and write the result as an artefact.
# -----------------------------------------------------------------------------
#  The 70/105/140-cell table these numbers form was typed into README.md AND into a
#  comment in solver.py, in two places, with no script behind either. Nothing in the
#  repository produced it, so nobody cloning the repository could check it, reproduce
#  it, or notice when it went stale — which is the exact failure mode the table itself
#  exists to warn about. It is computed here instead, from the same definitions the
#  solver uses for `Vm_peak_mps` and `erosional_exceedance_km`:
#
#      max        max of the ensemble-median volumetric flux j along the route
#      p99/p95    percentiles of |j| weighted by CELL LENGTH, so a refined grid does
#                 not change the weighting simply by having more cells
#      km over    route length where |j| exceeds the API RP 14E limit
#      1 km mean  peak of a 1 km running mean of |j|
#
#  Run:  python3 case/scripts/run_erosional_refinement.py
#  Writes case/outputs_steady/erosional_refinement.json.
# =============================================================================
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _paths import CASE  # noqa: E402

import solver  # noqa: E402

GRIDS = (70, 105, 140)
OUT = os.path.join(CASE, "outputs_steady", "erosional_refinement.json")


def _case(n_cells):
    """The as-operated case study at a given resolution — built by the case study's own
    `build_case`, so this measures the shipped case and not a second definition of it."""
    import run_case_study10 as RCS
    return RCS.build_case("erosional refinement — as-operated", "asoperated", 48.0,
                          n_cells=n_cells)


def _stats(sv):
    """The five erosional statistics, on the solver's own definitions."""
    r, f = sv.results, sv.case.fluids
    rho_m = np.nanmedian(r["alpha_l"] * sv.rho_l + (1 - r["alpha_l"])
                         * solver.gas_density(r["p"], r["T"], f), 1)
    j_med = np.nanmedian(r["j"], 1)
    jm = np.abs(j_med)
    eros = f.api14e_C_factor / math.sqrt(max(np.nanmean(rho_m), 1.0))
    cell_m = np.gradient(sv.x)
    #  LENGTH-weighted percentiles: an unweighted percentile over cells would change
    #  meaning with the grid, which is the very thing being measured here.
    order = np.argsort(jm)
    w = cell_m[order] / np.sum(cell_m)
    cw = np.cumsum(w)

    def pct(q):
        return float(np.interp(q / 100.0, cw, jm[order]))

    over = jm > eros
    km_over = float(np.sum(cell_m[over]) / 1000.0)
    #  peak of a 1 km running mean, in cells (the window is a length, not a cell count)
    dx = float(np.mean(cell_m))
    win = max(1, int(round(1000.0 / dx)))
    if win <= jm.size:
        run = np.convolve(jm, np.ones(win) / win, mode="valid")
        peak_1km = float(np.max(run))
    else:
        peak_1km = float(np.mean(jm))
    return {"n_cells": int(sv.case.pipeline.n_cells), "dx_m": dx,
            "erosional_limit_mps": float(eros),
            "max_mps": float(np.nanmax(j_med)), "p99_mps": pct(99.0), "p95_mps": pct(95.0),
            "km_over_limit": km_over, "peak_1km_mean_mps": peak_1km,
            "max_over_limit_ratio": float(np.nanmax(j_med) / eros),
            "argmax_km": float(sv.x[int(np.nanargmax(j_med))] / 1000.0)}


def main():
    rows = []
    for n in GRIDS:
        t0 = time.time()
        sv = solver.TransientSHCT(_case(n))
        sv.run(verbose=False)
        row = _stats(sv)
        row["runtime_s"] = round(time.time() - t0, 1)
        rows.append(row)
        print(f"  n_cells={n:4d}  dx={row['dx_m']:6.1f} m  max={row['max_mps']:6.3f}  "
              f"p99={row['p99_mps']:6.3f}  p95={row['p95_mps']:6.3f}  "
              f"km_over={row['km_over_limit']:6.3f}  1km={row['peak_1km_mean_mps']:6.3f}  "
              f"({row['runtime_s']:.0f} s)", flush=True)

    def spread(key):
        v = [r[key] for r in rows]
        return float(max(v) / max(min(v), 1e-12))

    spreads = {k: round(spread(k), 3) for k in
               ("max_mps", "p99_mps", "p95_mps", "km_over_limit", "peak_1km_mean_mps")}
    best = min(spreads, key=lambda k: spreads[k])
    rep: dict = {
        "what": "grid sensitivity of the erosional statistics on the as-operated case",
        "why": ("Vm_peak_mps is a point maximum taken next to a flow reversal. This measures "
                "how far it, and four alternatives to it, move under refinement, so the "
                "README quotes a measurement rather than a recollection."),
        "definitions": {
            "max_mps": "max of the ensemble-median flux j (the solver's Vm_peak_mps)",
            "p99_mps": "99th percentile of |j|, weighted by cell length",
            "p95_mps": "95th percentile of |j|, weighted by cell length",
            "km_over_limit": "route length with |j| above the API RP 14E limit "
                             "(the solver's erosional_exceedance_km)",
            "peak_1km_mean_mps": "peak of a 1 km running mean of |j|",
        },
        "grids": list(GRIDS),
        "rows": rows,
        "spread_max_over_min": spreads,
        "most_stable_statistic": best,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(rep, fh, indent=1)
        fh.write("\n")
    print("\n  spreads (max/min):", spreads)
    print(f"  most stable under refinement: {best}")
    print("  wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
