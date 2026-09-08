#!/usr/bin/env python3
# =============================================================================
#  shct_compositional.py — DEEP COMPOSITIONAL / PVT TRACKING & REPORT
# -----------------------------------------------------------------------------
#  Walks the solved pressure/temperature trajectory of a run and, at each axial
#  station, performs a Peng-Robinson vapour-liquid flash (shct_eos) to report the
#  full compositional/PVT state ALONG THE LINE:
#     * vapour mole fraction V (gas/liquid split) — where the fluid is two-phase,
#       where it is single-phase liquid (retrograde/bubble behaviour),
#     * per-component K-values K_i = y_i / x_i (which components partition to gas),
#     * phase densities & viscosities (PR + Peneloux shift + Lee/LBC),
#     * gas specific gravity and Z-factor.
#  Output: csv_compositional.csv (one row per station) + charts. Pure
#  post-processing; the core solver is unchanged. Composition defaults to the
#  EOS DEFAULT_COMPOSITION when the case carries none.
# =============================================================================
from __future__ import annotations

import os

import numpy as np

import shct_eos

#  DPI follows SHCT_FIG_DPI (default 320) so every generated figure meets the
#  journal artwork minimum of 300 dpi; a hard-coded 150/155 silently fell short.
import shct_style as _S

_FIG_DPI = _S.FIG_DPI

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:                                       # pragma: no cover
    _HAVE_MPL = False

NAVY = "#2E5BBF"; ACCENT = "#1F8AC0"; RED = "#E0463C"; ORANGE = "#E8842B"; TEAL = "#1AA0A0"
GREEN = "#3FA65A"; PURPLE = "#8E5CC8"
_KEY = ["C1", "C2", "C3", "CO2", "N2", "nC4", "nC5", "C7+"]   # components to highlight if present


def compositional_report(sv, outdir, n_stations=40):
    """Flash the PR EOS along the solved (P,T) line and write the compositional report."""
    os.makedirs(outdir, exist_ok=True)
    comp = getattr(sv.case.fluids, "composition", None) or shct_eos.DEFAULT_COMPOSITION
    names, _ = shct_eos._normalise(comp)
    r = sv.results
    def med(A):
        return np.nanmedian(A, 1)
    x_km = sv.x / 1000.0
    P = med(r["p"]); T = med(r["T"])
    # sample stations evenly
    idx = np.unique(np.linspace(0, len(x_km) - 1, min(n_stations, len(x_km))).astype(int))

    recs = []
    for i in idx:
        fl = shct_eos.flash(float(P[i]), float(T[i]), comp)
        props = shct_eos.eos_properties(float(P[i]), float(T[i]), comp)
        with np.errstate(divide="ignore", invalid="ignore"):
            K = np.where(fl["x"] > 1e-12, fl["y"] / np.maximum(fl["x"], 1e-12), np.nan)
        recs.append({"i": i, "x_km": float(x_km[i]), "P": float(P[i]), "T": float(T[i]),
                         "V": float(fl["V"]), "rho_g": props["rho_gas"], "rho_l": props["rho_oil"],
                         "mu_g": props["mu_gas"], "mu_l": props["mu_oil"], "Z": props["Z_gas"],
                         "sg": props["gas_sg"], "K": {n: float(K[j]) for j, n in enumerate(names)}})

    # --- CSV ---
    kcols = [f"K_{n}" for n in names]
    cols = ["x_km", "P_bar", "T_C", "vapour_frac_V", "rho_gas_kgm3", "rho_liq_kgm3",
            "mu_gas_Pas", "mu_liq_Pas", "Z_gas", "gas_sg"] + kcols
    with open(os.path.join(outdir, "csv_compositional.csv"), "w") as fh:
        fh.write(",".join(cols) + "\n")
        for d in recs:
            base = [d["x_km"], d["P"], d["T"], d["V"], d["rho_g"], d["rho_l"],
                    d["mu_g"], d["mu_l"], d["Z"], d["sg"]]
            kk = [d["K"][n] for n in names]
            fh.write(",".join(f"{v:.5g}" for v in base + kk) + "\n")

    if not _HAVE_MPL:
        return os.path.join(outdir, "csv_compositional.csv")

    _plot_pvt(recs, names, outdir)
    return os.path.join(outdir, "csv_compositional.csv")


def _phase_mask(recs, want):
    """NaN out the stations at which the requested phase does not exist.

    Below the bubble point the flash returns V = 0 and every K = 1, i.e. no split
    at all, and `eos_properties` then hands back the SINGLE-PHASE root for both
    `rho_gas` and `rho_oil`. Plotting that as "gas" drew a 511-559 kg/m3 gas
    density and a gas viscosity that was really the liquid's over 30 of the 32 km
    -- and, because the two curves coincide exactly, the liquid curve was hidden
    underneath the gas one for the whole trunk. Neither phase reading is wrong in
    the CSV; drawing the absent one is what misleads.
    """
    out = []
    for d in recs:
        two_phase = 1e-9 < d["V"] < 1.0 - 1e-9
        out.append(d[want] if two_phase else float("nan"))
    return np.array(out, float)


def _plot_pvt(recs, names, outdir):
    #  Journal artwork carries no chart titles -- the caption does that work. Every
    #  other figure generator routes its titles through solver._ttl; this one drew
    #  them unconditionally, so Fig. 1 was the only figure in the manuscript still
    #  showing them.
    from solver import _ttl
    xs = np.array([d["x_km"] for d in recs])
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.0))
    #  where the fluid is single-phase, say so on the figure rather than drawing
    #  a phantom second phase
    Vv = np.array([d["V"] for d in recs], float)
    sp = Vv <= 1e-9
    x_bub = float(xs[~sp][0]) if (~sp).any() and sp.any() else None
    # (a) vapour fraction
    ax[0, 0].plot(xs, [d["V"] for d in recs], color=NAVY, lw=1.8)
    ax[0, 0].set_ylabel("vapour mole fraction V"); ax[0, 0].set_ylim(-0.02, 1.02)
    ax[0, 0].set_title(_ttl("Gas/liquid split V(x) (PR flash)"), color=NAVY, fontweight="bold", fontsize=9.5)
    ax[0, 0].grid(alpha=.25)
    # (b) K-values of the key components (log)
    palette = [RED, ORANGE, GREEN, TEAL, PURPLE, NAVY, ACCENT, "#9AA8C7"]
    shown = [n for n in _KEY if n in names] or names[:6]
    for ci, n in enumerate(shown):
        Kn = np.array([d["K"].get(n, np.nan) for d in recs], float)
        ax[0, 1].plot(xs, np.where(sp, np.nan, Kn), lw=1.5,
                      color=palette[ci % len(palette)], label=n)
    ax[0, 1].set_yscale("log"); ax[0, 1].axhline(1.0, color="#3A5BA8", ls=":", lw=0.8)
    ax[0, 1].set_ylabel("K-value (y/x)"); ax[0, 1].legend(fontsize=7, ncol=2, framealpha=.85)
    ax[0, 1].set_title(_ttl("Component K-values along line"), color=NAVY, fontweight="bold", fontsize=9.5)
    ax[0, 1].grid(alpha=.25, which="both")
    # (c) phase densities
    ax[1, 0].plot(xs, [d["rho_l"] for d in recs], color=ACCENT, lw=1.8, label="liquid ρ_l")
    ax[1, 0].plot(xs, _phase_mask(recs, "rho_g"), color=RED, lw=1.8,
                  label="gas ρ_g (two-phase only)")
    ax[1, 0].set_ylabel("density (kg/m³)"); ax[1, 0].set_xlabel("distance (km)")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.25)
    ax[1, 0].set_title(_ttl("Phase densities (PR + Peneloux)"), color=NAVY, fontweight="bold", fontsize=9.5)
    # (d) phase viscosities
    ax[1, 1].plot(xs, [d["mu_l"] * 1000 for d in recs], color=ACCENT, lw=1.8, label="liquid μ_l (cP)")
    ax[1, 1].plot(xs, _phase_mask(recs, "mu_g") * 1e6, color=RED, lw=1.8,
                  label="gas μ_g (µPa·s, two-phase only)")
    ax[1, 1].set_xlabel("distance (km)"); ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.25)
    ax[1, 1].set_title(_ttl("Phase viscosities (Lee / LBC)"), color=NAVY, fontweight="bold", fontsize=9.5)
    #  every panel keeps the FULL route on the x axis. Masking the single-phase
    #  reach leaves the K-value panel with data only over the last kilometre, and
    #  matplotlib then autoscales it to that sliver -- four panels of the same
    #  figure ended up on two different distance scales.
    for a in ax.ravel():
        a.set_xlim(float(xs.min()) - 0.5, float(xs.max()) + 0.5)
    if x_bub is not None:
        #  Shade the two-phase reach. Without it the K-value and gas panels read as
        #  broken plots -- the curves occupy the last kilometre of a 32 km axis
        #  because that is the only place a vapour phase exists. The band says so.
        for a in ax.ravel():
            a.axvspan(x_bub, float(xs.max()) + 0.5, color="#DCE4F2", alpha=0.55,
                      zorder=0, lw=0)
            a.axvline(x_bub, color="#6B7A99", ls="--", lw=1.0, zorder=1)
        ax[0, 1].text(0.985, 0.04, "two-phase", transform=ax[0, 1].transAxes,
                      ha="right", va="bottom", fontsize=7.5, color="#3A4A6B")
        ax[0, 0].annotate(f"bubble point ≈ {x_bub:.1f} km\nsingle-phase liquid upstream",
                          xy=(x_bub, 0.12), xytext=(0.06, 0.55),
                          textcoords="axes fraction", fontsize=7.5, color="#3A4A6B",
                          arrowprops={"arrowstyle": "->", "color": "#6B7A99", "lw": 0.9})
    fig.suptitle(_ttl("Compositional / PVT tracking along the line (Peng-Robinson EOS)"),
                 color=NAVY, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(outdir, "compo_pvt.png"), dpi=_FIG_DPI); plt.close(fig)
    return os.path.join(outdir, "compo_pvt.png")


def replot_from_csv(csv_path, outdir=None):
    """Redraw compo_pvt.png from an existing csv_compositional.csv.

    The figure and the CSV come out of the same run, so a corrected DRAWING of an
    unchanged dataset does not need the solver re-run. Same plotting routine as
    the live path, so the two cannot drift.
    """
    import csv as _csv
    outdir = outdir or os.path.dirname(os.path.abspath(csv_path))
    with open(csv_path, newline="") as _fh:
        rows = list(_csv.DictReader(_fh))
    names = [k[2:] for k in rows[0] if k.startswith("K_")]
    recs = [{"x_km": float(r["x_km"]), "P": float(r["P_bar"]), "T": float(r["T_C"]),
                 "V": float(r["vapour_frac_V"]), "rho_g": float(r["rho_gas_kgm3"]),
                 "rho_l": float(r["rho_liq_kgm3"]), "mu_g": float(r["mu_gas_Pas"]),
                 "mu_l": float(r["mu_liq_Pas"]), "Z": float(r["Z_gas"]), "sg": float(r["gas_sg"]),
                 "K": {n: float(r[f"K_{n}"]) for n in names}} for r in rows]
    return _plot_pvt(recs, names, outdir)
