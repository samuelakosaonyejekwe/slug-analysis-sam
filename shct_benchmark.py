#!/usr/bin/env python3
# =============================================================================
#  shct_benchmark.py — compare SHCT against a reference simulator (OLGA,
#  LedaFlow, or any transient multiphase code) on an identical case.
# -----------------------------------------------------------------------------
#  WHAT THIS IS, AND WHAT IT IS NOT
#  --------------------------------
#  This module runs the comparison. It does NOT ship reference results: no OLGA
#  licence is available in this environment, and inventing benchmark numbers
#  would be worse than having none — a fabricated agreement is indistinguishable
#  from a real one until someone tries to reproduce it, and then the whole paper
#  is in question. So the reference data is supplied by whoever has the licence,
#  in the documented schema below, and this module does the rest: it builds the
#  matching SHCT case, runs it, interpolates both onto a common grid, computes
#  the agreement metrics a reviewer will ask for, and draws the comparison.
#
#  HOW TO PRODUCE THE REFERENCE FILE
#  ---------------------------------
#  Set the reference tool up with the SAME geometry, fluid and boundary
#  conditions as the case (case/outputs_*/input_data_deck.csv and
#  feed_composition.csv give every number), run it to the same end time, and
#  export the along-line profiles. Then write them as JSON:
#
#      {
#        "tool":      "OLGA 2023.1",              # name and version
#        "case":      "deepwater medium-crude tie-back, as-operated, 48 h",
#        "operator":  "who ran it, and when",
#        "notes":     "closure choices, grid, anything that would change the answer",
#        "t_end_h":   48.0,
#        "x_km":      [0.0, 0.46, ...],           # profile stations
#        "holdup":    [0.35, 0.36, ...],          # liquid holdup [-]
#        "P_bar":     [150.0, 148.2, ...],        # pressure [bar]
#        "T_C":       [58.0, 55.1, ...],          # temperature [degC]
#        "vm_mps":    [4.9, 4.8, ...]             # mixture velocity [m/s]  (optional)
#      }
#
#  Any of holdup / P_bar / T_C / vm_mps may be omitted; whatever is present is
#  compared. Put the file in validation/data/ and run:
#
#      python3 shct_benchmark.py validation/data/olga_asoperated.json
#
#  WHAT THE METRICS MEAN
#  ---------------------
#  For each field: mean absolute error, RMSE, the largest single deviation and
#  where it occurs, and the normalised RMSE (RMSE over the reference range) so
#  fields with different units can be compared. Two independent transient codes
#  agreeing to a few per cent on holdup along a 32 km line is a strong result;
#  the point of the comparison is to state the number rather than to claim it.
#
#  Author: Akosa Samuel Onyejekwe.
# =============================================================================
from __future__ import annotations

import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import shct_style as S

S.apply_style()

_DPI = int(os.environ.get("SHCT_FIG_DPI", "320"))

FIELDS = [
    ("holdup", "liquid holdup  α$_l$  [-]", "alpha_l"),
    ("P_bar", "pressure  [bar]", "p"),
    ("T_C", "temperature  [°C]", "T"),
    ("vm_mps", "mixture velocity  [m s$^{-1}$]", "j"),
]


def load_reference(path):
    """Read a reference-tool export and check it carries what it claims to."""
    with open(path) as fh:
        ref = json.load(fh)
    for required in ("tool", "x_km"):
        if required not in ref:
            raise ValueError(f"{path}: reference file has no '{required}' entry — "
                             f"a benchmark must say which tool produced it and at "
                             f"which stations")
    x = np.asarray(ref["x_km"], float)
    if x.ndim != 1 or x.size < 3:
        raise ValueError(f"{path}: 'x_km' must be a profile of at least 3 stations")
    present = [k for k, _lab, _s in FIELDS if k in ref]
    if not present:
        raise ValueError(f"{path}: none of {[k for k,_,_ in FIELDS]} present — "
                         f"there is nothing to compare")
    for k in present:
        if len(ref[k]) != x.size:
            raise ValueError(f"{path}: '{k}' has {len(ref[k])} values against "
                             f"{x.size} stations")
    return ref, present


def compare(sv, ref, present):
    """Interpolate both onto the reference stations and score the agreement."""
    x_ref = np.asarray(ref["x_km"], float)
    x_shct = np.asarray(sv.x, float) / 1000.0
    med = lambda A: np.nanmedian(np.asarray(A, float), 1)

    out = {}
    for key, _label, res_key in FIELDS:
        if key not in present:
            continue
        r = np.asarray(ref[key], float)
        try:
            mine_full = med(sv.results[res_key])
        except Exception:
            continue
        mine = np.interp(x_ref, x_shct, mine_full)
        d = mine - r
        fin = np.isfinite(d) & np.isfinite(r)
        if not fin.any():
            continue
        d, r_f = d[fin], r[fin]
        rng = float(np.nanmax(r_f) - np.nanmin(r_f))
        i_worst = int(np.argmax(np.abs(d)))
        out[key] = {
            "mae": float(np.mean(np.abs(d))),
            "rmse": float(np.sqrt(np.mean(d ** 2))),
            "max_abs_dev": float(np.max(np.abs(d))),
            "max_abs_dev_at_km": float(x_ref[fin][i_worst]),
            "nrmse_pct": float(100.0 * np.sqrt(np.mean(d ** 2)) / rng) if rng > 0 else float("nan"),
            "reference": r.tolist(),
            "shct": mine.tolist(),
            "x_km": x_ref.tolist(),
        }
    return out


def figure(ref, scores, outdir, tag="benchmark"):
    """Profiles overlaid, and the deviation beneath, for every compared field."""
    keys = [k for k, _l, _s in FIELDS if k in scores]
    if not keys:
        return None
    labels = {k: lab for k, lab, _s in FIELDS}
    fig, axes = plt.subplots(2, len(keys), figsize=(4.6 * len(keys), 6.4),
                             squeeze=False,
                             gridspec_kw=dict(height_ratios=[2.0, 1.0]))
    for j, k in enumerate(keys):
        sc = scores[k]
        x = np.asarray(sc["x_km"], float)
        a, b = axes[0][j], axes[1][j]
        a.plot(x, sc["reference"], color=S.RED, lw=2.0, marker="o", ms=3.4,
               markerfacecolor="white", label=ref.get("tool", "reference"))
        a.plot(x, sc["shct"], color=S.BLUE, lw=1.8, ls="--", label="SHCT")
        a.set_ylabel(labels[k], fontsize=9)
        a.set_xlim(x.min(), x.max())
        for sp in a.spines.values():
            sp.set_color(S.INK)
        a.grid(True, color=S.GRIDC, lw=0.6, ls=":")
        a.set_axisbelow(True)
        a.tick_params(labelsize=8)
        a.set_title(f"NRMSE {sc['nrmse_pct']:.1f} %", fontsize=9,
                    color=S.TITLE, fontweight="bold", pad=5)
        if j == 0:
            a.legend(loc="upper left", bbox_to_anchor=(0.0, -0.02), fontsize=7.5,
                     framealpha=1.0, facecolor="white", edgecolor=S.INK, ncol=2)

        dev = np.asarray(sc["shct"], float) - np.asarray(sc["reference"], float)
        b.axhline(0.0, color=S.INK, lw=0.9)
        b.plot(x, dev, color=S.TEAL, lw=1.5)
        b.fill_between(x, 0, dev, color=S.TEAL, alpha=0.18)
        b.set_xlabel("distance from wellhead  [km]", fontsize=9)
        b.set_ylabel("SHCT − reference", fontsize=8.5)
        b.set_xlim(x.min(), x.max())
        for sp in b.spines.values():
            sp.set_color(S.INK)
        b.grid(True, color=S.GRIDC, lw=0.6, ls=":")
        b.set_axisbelow(True)
        b.tick_params(labelsize=8)

    fig.suptitle(f"SHCT against {ref.get('tool', 'a reference simulator')} — "
                 f"{ref.get('case', 'identical case')}",
                 color=S.TITLE, fontweight="bold", fontsize=10.5, y=0.995)
    fig.text(0.5, 0.005,
             f"Reference produced by {ref.get('operator', 'an independent run')}; "
             f"identical geometry, fluid and boundary conditions. "
             f"{ref.get('notes', '')}".strip(),
             ha="center", fontsize=6.8, color=S.INK, style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.955))
    p = os.path.join(outdir, f"{tag}_vs_reference.png")
    fig.savefig(p, dpi=_DPI)
    plt.close(fig)
    return p


def run(ref_path, outdir=None, case_builder=None, t_end_h=None):
    """Load a reference export, run the matching SHCT case, score and plot."""
    import solver

    ref, present = load_reference(ref_path)
    outdir = outdir or os.path.dirname(os.path.abspath(ref_path))
    os.makedirs(outdir, exist_ok=True)

    if case_builder is None:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "case", "scripts"))
        import run_case_study10 as R
        case_builder = lambda t: R.build_case("benchmark", "asoperated", t)

    case = case_builder(float(t_end_h or ref.get("t_end_h", 48.0)))
    sv = solver.TransientSHCT(case)
    sv.run(verbose=False)

    scores = compare(sv, ref, present)
    fig_path = figure(ref, scores, outdir)

    report = {
        "reference_tool": ref.get("tool"),
        "reference_case": ref.get("case"),
        "reference_operator": ref.get("operator"),
        "reference_notes": ref.get("notes"),
        "fields_compared": list(scores.keys()),
        "metrics": {k: {m: v for m, v in sc.items()
                        if m not in ("reference", "shct", "x_km")}
                    for k, sc in scores.items()},
        "figure": os.path.basename(fig_path) if fig_path else None,
    }
    with open(os.path.join(outdir, "benchmark_report.json"), "w") as fh:
        solver.dump_json(report, fh)

    print(f"SHCT vs {ref.get('tool')} — {ref.get('case', '')}")
    for k, m in report["metrics"].items():
        print(f"  {k:9s} MAE {m['mae']:.4g}   RMSE {m['rmse']:.4g}   "
              f"NRMSE {m['nrmse_pct']:.1f} %   worst {m['max_abs_dev']:.4g} "
              f"at {m['max_abs_dev_at_km']:.1f} km")
    if fig_path:
        print(f"  -> {fig_path}")
    return report


def main(argv):
    if not argv:
        print(__doc__ or "")
        print("usage: python3 shct_benchmark.py <reference.json> [outdir]")
        print("\nNo reference file is shipped: this environment has no licence for "
              "a reference simulator, and fabricated benchmark numbers would be "
              "worse than none. Export your own run in the schema documented at "
              "the top of this file.")
        return 2
    return 0 if run(argv[0], argv[1] if len(argv) > 1 else None) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
