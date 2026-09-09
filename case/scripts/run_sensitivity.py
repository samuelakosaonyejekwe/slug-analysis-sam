#!/usr/bin/env python3
# =============================================================================
#  run_sensitivity.py — one-at-a-time sensitivity of the coupled predictions to
#  the four assumed kinetic/coupling constants.
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
#    f_slug_floor_Hz 1e-5 .. 1e-3 — two decades either side of the 1e-4 default.
#                           Phi_SH goes as 1/f_slug, so wherever the line is not slugging
#                           (startup, and the whole of a shut-in) this numerical floor —
#                           not the physics — sets the magnitude of the coupling number.
#
#  Reduced but SELF-CONSISTENT fidelity: every run here, the baseline included,
#  uses the same ensemble size and simulated window, so the movements are
#  comparable with one another.  They are NOT directly comparable with the
#  headline case-study numbers, which use n_ensemble=12 / t_end=48 h.
# =============================================================================
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_case_study10 as R
from _paths import OUT

import solver

N_ENSEMBLE = 6
N_CELLS    = 70
T_END_H    = 24.0

METRICS = ["max_Phi_SH", "sustained_Phi_SH", "sustained_Phi_SH_hotspot_km",
           "sustained_supercritical_km", "final_Phi_SH",
           "Phi_SH_supercritical_time_frac", "Phi_SH_peak_time_h",
           "max_Phi_SH_uncapped", "max_Psi_kinetic_ratio",
           "Phi_SH_above_critical_frac", "Phi_SH_critical", "deposit_ref_mm",
           "P_plug", "time_to_plug_P50_h", "time_to_plug_P10_h",
           "time_to_plug_P90_h", "MEG_wt_pct", "under_inhibited_km",
           "peak_deposit_mm", "max_subcooling_C", "cooldown_to_hydrate_h"]

BASE = {"kg0_mult": 1.0, "growth_exp_n": 1.0, "C_phi": 1500.0,
        "f_slug_floor_Hz": 1.0e-4}


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
    for ff in (1e-5, 3e-5, 3e-4, 1e-3):
        d = dict(BASE); d["f_slug_floor_Hz"] = ff
        runs.append((f"ffloor_{ff:g}", d))
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
    k.f_slug_floor_Hz = p["f_slug_floor_Hz"]
    sv = solver.TransientSHCT(case)
    sv.run()
    eng = sv.engineering()
    row = {"label": label, "kg0_mult": p["kg0_mult"],
           "growth_exp_n": p["growth_exp_n"], "C_phi": p["C_phi"],
           "f_slug_floor_Hz": p["f_slug_floor_Hz"],
           "runtime_s": round(time.time() - t0, 1)}
    row.update({m: eng.get(m) for m in METRICS})
    print(f"  done {label:12s} ({row['runtime_s']:.0f}s)  "
          f"Phi_max={row['max_Phi_SH']}  P50={row['time_to_plug_P50_h']}  "
          f"MEG={row['MEG_wt_pct']}", flush=True)
    return row


def plot(rows, outdir):
    """Fig. 18 of the IJMF manuscript: how far the predictions move with the four
    unfitted constants. All four P50 panels share one y-scale so that a flat
    response (C) cannot be mistaken for a strong one by axis choice alone."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    #  one hue for the left-axis series (see the loop below) and one for the MEG axis
    BLUE, RED = "#2E5BBF", "#E0463C"

    #  NOTHING NEED PLUG. With the gas-gravity defect fixed the as-operated line is
    #  sub-critical, so time_to_plug_P50_h is null on most rows of this sweep and can be
    #  null on ALL of them -- at which point `max()` over the filtered generator raises
    #  ValueError, main() catches it, and Fig. 18 silently stops being produced while the
    #  CSV beside it looks fine. Fall back to the coupling number, which is what varies
    #  when the duty does not plug, and say on the figure that that is what happened.
    def _p50(r):
        v = r.get("time_to_plug_P50_h")
        return v if isinstance(v, (int, float)) and v == v else None

    _p50s = [_p50(r) for r in rows if _p50(r) is not None]
    plugs = bool(_p50s)
    if plugs:
        ykey, ylab = "time_to_plug_P50_h", r"$t_{\rm plug,P50}$  (h)"
        ymax = max(_p50s) * 1.08
    else:
        ykey, ylab = "max_Phi_SH_uncapped", r"peak $\Phi_{SH}$  (-)"
        _phis = [r.get(ykey) for r in rows
                 if isinstance(r.get(ykey), (int, float)) and r.get(ykey) == r.get(ykey)]
        ymax = (max(_phis) * 1.08) if _phis else 1.0

    def grp(pfx, xk):
        rs = [r for r in rows if r["label"].startswith(pfx) or r["label"] == "baseline"]
        rs.sort(key=lambda r: r[xk])
        return [r[xk] for r in rs], rs

    _meg = [float(r["MEG_wt_pct"]) for r in rows
            if r.get("MEG_wt_pct") not in (None, "")]
    if _meg:
        _lo, _hi = min(_meg), max(_meg)
        _pad = max(0.08 * (_hi - _lo), 1.0)
        meg_lo, meg_hi = max(0.0, _lo - _pad), min(100.0, _hi + _pad)
    else:
        meg_lo, meg_hi = 0.0, 100.0

    fig, ax = plt.subplots(1, 4, figsize=(14.4, 3.6), sharey=True)
    specs = [("kg0", "kg0_mult",     r"$k_{g0}$ multiplier",
              "(a) growth-rate prefactor $k_{g0}$"),
             ("n_",  "growth_exp_n", r"subcooling exponent $n$",
              "(b) subcooling exponent $n$"),
             ("C_",  "C_phi",        r"coupling coefficient $C$",
              "(c) coupling coefficient $C$"),
             ("ffloor", "f_slug_floor_Hz", r"slug-frequency floor  (Hz)",
              "(d) slug-frequency floor")]
    for i, (a, (pfx, xk, xlab, ttl)) in enumerate(zip(ax, specs)):
        x, rs = grp(pfx, xk)
        yv = [r.get(ykey) for r in rs]
        #  ONE COLOUR FOR THE LEFT-AXIS SERIES, because there is one legend and it is built
        #  from panel (a). Giving each panel its own hue put an orange curve in (b), a green
        #  one in (c) and a purple one in (d) that the legend never named, while claiming the
        #  quantity was blue. The panels are already told apart by their titles and x-axes.
        a.plot(x, yv, "o-", color=BLUE, lw=2, ms=5,
               label=(r"$t_{\rm plug,P50}$" if plugs else r"peak $\Phi_{SH}$"))
        #  A PANEL WITH NO FINITE POINT DRAWS NOTHING, and an empty panel beside three full
        #  ones reads as lost data rather than as "this constant changes nothing". On the
        #  as-operated sweep C_phi and the slug-frequency floor never produce a plug at all,
        #  so both panels were blank. Say which it is.
        if plugs and not any(isinstance(v, (int, float)) and v == v for v in yv):
            #  ABOVE the axes, not in the middle of them. Centred at (0.5, 0.5) this note
            #  landed squarely on the MEG-dose curve it shares the panel with -- the one
            #  series these two panels DO have -- which is the thing the project's figures
            #  are not allowed to do. There is nothing above the axes but the title.
            a.text(0.5, 1.02, "no realisation plugs anywhere on this sweep",
                   transform=a.transAxes, ha="center", va="bottom", fontsize=7.5,
                   style="italic", color="#3A5BA8",
                   bbox={"boxstyle": "round,pad=0.25", "fc": "white",
                         "ec": "#D2DCF2", "lw": 0.8})
        a.set_xlabel(xlab); a.set_ylim(0, ymax); a.grid(alpha=0.3)
        if i == 0:
            a.set_ylabel(ylab)
        if pfx in ("kg0", "ffloor"):
            a.set_xscale("log")
        b = a.twinx()
        b.plot(x, [r["MEG_wt_pct"] for r in rs], "s--", color=RED, lw=1.5, ms=4, alpha=0.85)
        #  the MEG axis must follow the SWEPT doses: a fixed 50-70 wt% window puts
        #  the curve off the plot as soon as a sweep moves the requirement outside
        #  it, and the panels share this axis so it is computed once over all rows.
        b.set_ylim(meg_lo, meg_hi)
        if i == len(specs) - 1:
            b.set_ylabel("required MEG dose (wt%)", color=RED)
            b.tick_params(axis="y", labelcolor=RED)
        else:
            b.set_yticklabels([])
        #  extra pad on the panels that carry the "no realisation plugs" note, which now
        #  sits between the axes and the title
        _no_plug = plugs and not any(isinstance(v, (int, float)) and v == v for v in yv)
        a.set_title(ttl, fontsize=10, pad=20 if _no_plug else 6)
    if not plugs:
        fig.text(0.5, 1.005, "no realisation plugs anywhere in this sweep — the left axis "
                             "shows the peak coupling number instead of a time-to-plug",
                 ha="center", va="bottom", fontsize=8.5, style="italic", color="#3A5BA8")
    h1, l1 = ax[0].get_legend_handles_labels()
    fig.legend(h1 + [plt.Line2D([], [], color=RED, ls="--", marker="s", ms=4)],
               l1 + ["required MEG dose"], loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.09), frameon=False, fontsize=9)
    fig.tight_layout()
    path = os.path.join(outdir, "13_sensitivity.png")
    #  IJMF artwork: >=300 dpi and >=1063 px for single column. 200 dpi cleared the
    #  pixel floor but failed the dpi metadata check, which Elsevier reads. The value
    #  comes from shct_style (default 320, SHCT_FIG_DPI overrides) so this figure moves
    #  with the rest of the set instead of carrying its own literal.
    import shct_style as _S
    fig.savefig(path, dpi=_S.FIG_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[sensitivity] -> {path}", flush=True)
    return path


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

    #  DROP THE COLUMNS THAT ARE EMPTY IN EVERY ROW, AND SAY WHICH. A ratio to a baseline
    #  of zero is undefined, so on a sub-critical duty eight of the _rel columns come out
    #  blank from top to bottom -- a header promising a sensitivity, over nothing. An empty
    #  column is not a measurement and it invites the reader to think the run failed. They
    #  are omitted from the CSV and listed in the JSON with the reason, so the omission is
    #  a statement rather than a gap.
    _rel_cols = [m + "_rel" for m in METRICS]
    _empty = [c for c in _rel_cols
              if all(r.get(c) is None or r.get(c) != r.get(c) for r in rows)]
    _kept = [c for c in _rel_cols if c not in _empty]
    _why = {c: (f"baseline {c[:-4]} = {base.get(c[:-4])!r}; a ratio to it is undefined"
                if not base.get(c[:-4]) else "undefined on every run")
            for c in _empty}
    cols = (["label", "kg0_mult", "growth_exp_n", "C_phi", "f_slug_floor_Hz", "runtime_s"] +
            METRICS + _kept)
    csv_path = os.path.join(outdir, "sensitivity_phiSH.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})
    if _empty:
        print(f"[sensitivity] {len(_empty)} relative column(s) omitted as undefined on "
              f"every run: {', '.join(_empty)}", flush=True)
    with open(os.path.join(outdir, "sensitivity_phiSH.json"), "w") as fh:
        solver.dump_json({"settings": {"n_ensemble": N_ENSEMBLE, "n_cells": N_CELLS,
                                "t_end_h": T_END_H, "scenario": "asoperated"},
                   "baseline": BASE, "rows": rows,
                   "relative_columns_omitted": _why}, fh, indent=2, default=str)
    print(f"[sensitivity] -> {csv_path}", flush=True)
    try:
        plot(rows, outdir)
    except Exception as exc:                       # plotting must never lose the data
        print(f"[sensitivity] plot skipped: {exc}", flush=True)
    return rows


if __name__ == "__main__":
    main()
