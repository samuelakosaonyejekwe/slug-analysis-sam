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

import json
import os
import textwrap

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
    #  A COLUMN HEADED "gas" MUST NOT CARRY THE LIQUID. Below the bubble point the flash
    #  returns V = 0 and y = z, so eos_properties hands back the FEED at its vapour root
    #  for every "gas" property: on the as-operated case 39 of 40 stations were written
    #  with rho_gas_kgm3 = rho_liq_kgm3 = 513 kg/m3, Z_gas = 0.54 and gas_sg = 1.7645 --
    #  the last of which is the very number this project's release notes identify as the
    #  defect ("the flash returns y = z when a mixture does not split ... on this live
    #  crude that is the liquid"). The figure beside this file already masks the absent
    #  phase; the CSV did not, and it is the CSV the report embeds, so 1.7645 appeared
    #  twelve times in report.docx under a gas-gravity heading. Every property of a phase
    #  that does not exist is now written EMPTY, and the state is named in its own column.
    #  K = y/x is masked for the same reason: with no split it is identically 1, which is
    #  an identity, not a measurement.
    def _state(V):
        if V <= 1e-9:
            return "single-phase liquid"
        if V >= 1.0 - 1e-9:
            return "single-phase gas"
        return "two-phase"

    def _fmt(v, defined=True):
        return f"{v:.5g}" if (defined and v is not None and v == v) else ""

    kcols = [f"K_{n}" for n in names]
    cols = ["x_km", "P_bar", "T_C", "phase_state", "vapour_frac_V", "rho_gas_kgm3",
            "rho_liq_kgm3", "mu_gas_Pas", "mu_liq_Pas", "Z_gas", "gas_sg"] + kcols
    with open(os.path.join(outdir, "csv_compositional.csv"), "w") as fh:
        fh.write(",".join(cols) + "\n")
        for d in recs:
            st = _state(d["V"])
            has_gas = st != "single-phase liquid"
            has_liq = st != "single-phase gas"
            two = st == "two-phase"
            row = [_fmt(d["x_km"]), _fmt(d["P"]), _fmt(d["T"]), st, _fmt(d["V"]),
                   _fmt(d["rho_g"], has_gas), _fmt(d["rho_l"], has_liq),
                   _fmt(d["mu_g"], has_gas), _fmt(d["mu_l"], has_liq),
                   _fmt(d["Z"], has_gas), _fmt(d["sg"], has_gas)]
            row += [_fmt(d["K"][n], two) for n in names]
            fh.write(",".join(row) + "\n")

    if not _HAVE_MPL:
        return os.path.join(outdir, "csv_compositional.csv")

    #  the flow model's own view of the same fluid, so the figure can state the difference.
    #  Persisted beside the CSV because replot_from_csv() has to reproduce this caption
    #  exactly; it used to call _plot_pvt() with no flow at all, which silently dropped
    #  the comparison and gave the redraw a different note from the live figure -- the one
    #  thing that routine's docstring promises cannot happen.
    _fl = sv.case.fluids
    _flow = {
        "gas_pct": float(np.nanmean(1.0 - med(r["alpha_l"]))) * 100.0,
        "rho_l": _fl.rho_oil * (1.0 - _fl.water_cut) + _fl.rho_water * _fl.water_cut,
        "mu_l_cP": _fl.mu_liquid * 1000.0,
        #  the OIL alone, which is what the EOS liquid is comparable to. The composite
        #  above is heavier only because it carries the water cut, and quoting it against
        #  the EOS hydrocarbon made a correct water cut read as a fluid mismatch.
        "rho_oil": float(_fl.rho_oil), "water_cut": float(_fl.water_cut),
        "mu_oil_cP": float(getattr(_fl, "mu_oil", _fl.mu_liquid)) * 1000.0}
    try:
        with open(os.path.join(outdir, "compo_flow_basis.json"), "w") as _bh:
            json.dump(_flow, _bh, indent=2)
    except OSError:
        pass
    _plot_pvt(recs, names, outdir, flow=_flow)
    return os.path.join(outdir, "csv_compositional.csv")


def _phase_mask(recs, want):
    """NaN out the stations at which the requested phase does not exist.

    Below the bubble point the flash returns V = 0 and every K = 1, i.e. no split
    at all, and `eos_properties` then hands back the SINGLE-PHASE root for both
    `rho_gas` and `rho_oil`. Plotting that as "gas" drew a 511-559 kg/m3 gas
    density and a gas viscosity that was really the liquid's over 30 of the 32 km
    -- and, because the two curves coincide exactly, the liquid curve was hidden
    underneath the gas one for the whole trunk.

    This docstring used to add "neither phase reading is wrong in the CSV". It was: a
    column headed rho_gas_kgm3 carrying the liquid's 513 kg/m3, and one headed gas_sg
    carrying the feed's 1.7645, are wrong wherever no gas exists, and the report embeds
    that CSV. The writer above now masks those columns on the same test this uses.
    """
    out = []
    for d in recs:
        two_phase = 1e-9 < d["V"] < 1.0 - 1e-9
        out.append(d[want] if two_phase else float("nan"))
    return np.array(out, float)


def _plot_pvt(recs, names, outdir, flow=None):
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
    sp = (Vv <= 1e-9) | (Vv >= 1.0 - 1e-9)
    x_bub = float(xs[~sp][0]) if (~sp).any() and sp.any() else None
    # (a) vapour fraction
    ax[0, 0].plot(xs, [d["V"] for d in recs], color=NAVY, lw=1.8)
    ax[0, 0].set_ylabel("vapour mole fraction V"); ax[0, 0].set_ylim(-0.02, 1.02)
    ax[0, 0].set_title(_ttl("Gas/liquid split V(x) (PR flash)"), color=NAVY, fontweight="bold", fontsize=9.5)
    ax[0, 0].grid(alpha=.25)
    # (b) K-values of the key components (log)
    palette = [RED, ORANGE, GREEN, TEAL, PURPLE, NAVY, ACCENT, "#9AA8C7"]
    shown = [n for n in _KEY if n in names] or names[:6]
    #  MARKERS, not lines alone. On this crude the fluid is two-phase at ONE of the 40
    #  stations, and a one-point line renders as nothing at all: panel (b) went out with a
    #  full eight-entry legend over a completely empty axis, and panels (c) and (d) drew a
    #  legend entry for a gas curve that was not visible either. A marker makes a
    #  single-station series a single visible point.
    _n_two = int(np.count_nonzero(~sp))
    #  a single-station gas series renders as nothing without a marker (panels c and d)
    _mk = {"marker": "o", "ms": 4.0} if _n_two <= 3 else {}
    #  A K-VALUE PROFILE NEEDS A PROFILE. Plotted against distance, a fluid that splits at
    #  one station out of forty gives an empty axis with eight dots stacked against its
    #  right-hand edge and an eight-entry legend for curves that do not exist -- the honest
    #  output of a degenerate case, but not a graph anyone can read. Where there is no
    #  profile, show the thing that does carry the information: the K-values of that one
    #  state, against COMPONENT, ordered heavy-to-light, which is how a split is actually
    #  read. Above a few two-phase stations the distance profile is meaningful and is kept.
    if _n_two and _n_two <= 3:
        _j = int(np.argmax(~sp))                       # the (first) station that splits
        _rec = recs[_j]
        _ks = [(n, float(_rec["K"].get(n, np.nan))) for n in names]
        _ks = [(n, k) for n, k in _ks if np.isfinite(k) and k > 0]
        _ks.sort(key=lambda t: t[1])                   # heaviest (smallest K) first
        _xi = np.arange(len(_ks))
        ax[0, 1].plot(_xi, [k for _, k in _ks], color=NAVY, lw=1.4, marker="o", ms=5,
                      zorder=3)
        ax[0, 1].set_xticks(_xi)
        ax[0, 1].set_xticklabels([n for n, _ in _ks], rotation=45, fontsize=7)
        ax[0, 1].set_xlabel(_S.label("component (heavy → light)", "component"), fontsize=8)
        ax[0, 1].set_title(_ttl(_S.label(f"K-values at the one station that splits "
                                         f"(x = {_rec['x_km']:.1f} km)",
                                         f"K at x = {_rec['x_km']:.1f} km")),
                           color=NAVY, fontweight="bold", fontsize=9.5)
        ax[0, 1].text(0.02, 0.965,
                      _S.label(f"the fluid is single-phase at {len(recs) - _n_two} of "
                               f"{len(recs)} stations,\nso there is no profile to draw",
                               f"single-phase at {len(recs) - _n_two}/{len(recs)}"),
                      transform=ax[0, 1].transAxes, ha="left", va="top", fontsize=7,
                      style="italic", color="#8A4B2A", linespacing=1.3)
    else:
        for ci, n in enumerate(shown):
            Kn = np.array([d["K"].get(n, np.nan) for d in recs], float)
            ax[0, 1].plot(xs, np.where(sp, np.nan, Kn), lw=1.5,
                          color=palette[ci % len(palette)], label=n)
        ax[0, 1].legend(fontsize=7, ncol=2, framealpha=.85, loc="upper left",
                        bbox_to_anchor=(1.012, 1.0), borderaxespad=0.0)
    ax[0, 1].set_yscale("log")
    ax[0, 1].axhline(1.0, color="#3A5BA8", ls=":", lw=0.8)
    ax[0, 1].set_ylabel("K-value (y/x)")
    if not (_n_two and _n_two <= 3):
        ax[0, 1].set_title(_ttl("Component K-values along line"),
                           color=NAVY, fontweight="bold", fontsize=9.5)
    ax[0, 1].grid(alpha=.25, which="both")
    # (c) phase densities
    ax[1, 0].plot(xs, [d["rho_l"] for d in recs], color=ACCENT, lw=1.8, label="liquid ρ_l")
    #  A TWO-POINT "PROFILE" IS NOT A PROFILE. Where the fluid splits at only a station
    #  or two, the gas series is two markers over 2.5 % of a 32 km axis, under a legend
    #  entry that promises a curve -- the same fault as the K-value panel. State the value
    #  instead, and leave the axis to the quantity that does vary along the line.
    _rg = _phase_mask(recs, "rho_g")
    if _n_two <= 3:
        _v = np.asarray(_rg, float); _v = _v[np.isfinite(_v)]
        if _v.size:
            ax[1, 0].text(0.975, 0.055,
                          _S.label(f"gas ρ_g exists at {_n_two} station(s) only: "
                                   f"{float(np.nanmean(_v)):.0f} kg/m³",
                                   f"ρ_g: {float(np.nanmean(_v)):.0f} kg/m³ "
                                   f"({_n_two} stn)"),
                          transform=ax[1, 0].transAxes, ha="right", va="bottom",
                          fontsize=7, style="italic", color="#8A4B2A")
    else:
        ax[1, 0].plot(xs, _rg, color=RED, lw=1.8, label="gas ρ_g (two-phase only)")
    ax[1, 0].set_ylabel(_S.label("density (kg/m³)", "ρ (kg/m³)"))
    ax[1, 0].set_xlabel(_S.label("distance (km)", "x (km)"))
    _S.legend_outside(ax[1, 0], fontsize=8); ax[1, 0].grid(alpha=.25)
    ax[1, 0].set_title(_ttl("Phase densities (PR + Peneloux)"), color=NAVY, fontweight="bold", fontsize=9.5)
    # (d) phase viscosities — BOTH IN cP, ON A LOG AXIS. The liquid was drawn in cP and
    #     the gas in µPa·s on the SAME linear axis, so a 0.16 cP liquid shared a scale
    #     with a 12.9 "µPa·s" gas: the axis ran 0-13 and the liquid curve lay flat on
    #     zero, unreadable. One unit, and a log axis because the two phases differ by an
    #     order of magnitude.
    ax[1, 1].plot(xs, [d["mu_l"] * 1000 for d in recs], color=ACCENT, lw=1.8,
                  label="liquid μ_l")
    _mg = _phase_mask(recs, "mu_g") * 1000.0
    if _n_two <= 3:
        _v = np.asarray(_mg, float); _v = _v[np.isfinite(_v)]
        if _v.size:
            ax[1, 1].text(0.975, 0.055,
                          _S.label(f"gas μ_g exists at {_n_two} station(s) only: "
                                   f"{float(np.nanmean(_v)):.4f} cP",
                                   f"μ_g: {float(np.nanmean(_v)):.4f} cP "
                                   f"({_n_two} stn)"),
                          transform=ax[1, 1].transAxes, ha="right", va="bottom",
                          fontsize=7, style="italic", color="#8A4B2A")
    else:
        ax[1, 1].plot(xs, _mg, color=RED, lw=1.8, label="gas μ_g (two-phase only)")
    #  LOG ONLY WHERE THERE IS A DECADE TO SPAN. The log axis is here because gas and
    #  liquid viscosity differ by an order of magnitude -- but when the gas is not a curve
    #  (it exists at one station and is stated in words instead) the liquid alone spans
    #  0.09-0.16 cP, a factor of 1.8, and a log axis over less than one decade produces
    #  ticks like "1.1 x 10^-1, 1.2 x 10^-1" that are harder to read than the numbers.
    #  Scale by what is actually plotted.
    _mu_vals = np.asarray([d["mu_l"] * 1000 for d in recs], float)
    _mu_vals = _mu_vals[np.isfinite(_mu_vals) & (_mu_vals > 0)]
    _decades = (np.log10(_mu_vals.max() / _mu_vals.min())
                if _mu_vals.size and _mu_vals.min() > 0 else 0.0)
    if _n_two > 3 or _decades >= 1.0:
        ax[1, 1].set_yscale("log")
    ax[1, 1].set_ylabel(_S.label("viscosity (cP)", "μ (cP)"))
    ax[1, 1].set_xlabel(_S.label("distance (km)", "x (km)"))
    _S.legend_outside(ax[1, 1], fontsize=8)
    ax[1, 1].grid(alpha=.25, which="both")
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
        #  ...but NOT on the K panel when it has been switched to a component axis: a
        #  bubble-point at "x = 30.8" would there mark component number 30.8, which is
        #  nothing at all.
        _component_axis = bool(_n_two and _n_two <= 3)
        for a in ax.ravel():
            if _component_axis and a is ax[0, 1]:
                continue
            a.axvspan(x_bub, float(xs.max()) + 0.5, color="#DCE4F2", alpha=0.55,
                      zorder=0, lw=0)
            a.axvline(x_bub, color="#6B7A99", ls="--", lw=1.0, zorder=1)
        #  the component-axis panel already carries its own, fuller note
        if not _component_axis:
            ax[0, 1].text(0.985, 0.04, "two-phase",
                          transform=ax[0, 1].transAxes,
                          ha="right", va="bottom", fontsize=7.5, color="#3A4A6B")
        _xlo, _xhi = float(np.nanmin(xs)), float(np.nanmax(xs))
        _fb = (x_bub - _xlo) / max(_xhi - _xlo, 1e-9)          # bubble point, axes fraction
        _left = _fb > 0.45                                     # room on the left?
        ax[0, 0].annotate(f"bubble point ≈ {x_bub:.1f} km\nsingle-phase liquid upstream",
                          xy=(x_bub, 0.12),
                          xytext=(_fb - 0.05 if _left else _fb + 0.05, 0.42),
                          textcoords="axes fraction", fontsize=7.5, color="#3A4A6B",
                          ha="right" if _left else "left",
                          arrowprops={"arrowstyle": "->", "color": "#6B7A99", "lw": 0.9,
                                      "shrinkB": 2})
    fig.suptitle(_ttl("Compositional / PVT tracking along the line (Peng-Robinson EOS)"),
                 color=NAVY, fontweight="bold")

    #  SAY WHAT THIS FIGURE DOES AND DOES NOT SHARE WITH THE REST OF THE SET. This note
    #  used to read "these are not the same fluid" on every run, which was true when the
    #  C7+ pseudo-component was still n-heptane's constants and the flash gave a 549 kg/m3
    #  liquid against the flow model's 858. With C7+ characterised (Kesler-Lee on a
    #  Riazi-Daubert T_b) the EOS liquid now lands on the flow model's OIL, and the only
    #  reason the flow model's LIQUID is heavier is the water cut it carries -- correct
    #  physics that the old wording reported as a defect. What genuinely still differs is
    #  the PHASE SPLIT, so the note now says that and nothing wider. The wording is derived
    #  from the numbers, not asserted, so it cannot go stale behind a later change.
    _rl = float(np.nanmedian([d["rho_l"] for d in recs]))
    _ml = float(np.nanmedian([d["mu_l"] for d in recs])) * 1000.0
    _note = (f"EOS on the hydrocarbon composition alone: two-phase at {_n_two} of "
             f"{len(recs)} stations, liquid {_rl:.0f} kg/m\u00b3 and {_ml:.2f} cP.")
    if flow:
        #  like against like: the EOS hydrocarbon liquid vs the flow model's OIL
        _ro = float(flow.get("rho_oil", flow["rho_l"]))
        _gap = 100.0 * (_ro - _rl) / max(_rl, 1e-9)
        _agree = abs(_gap) <= 5.0
        _wc = 100.0 * float(flow.get("water_cut", 0.0))
        if _S.compact():
            _note = (f"EOS on the hydrocarbon alone ({_n_two}/{len(recs)} two-phase); "
                     + (f"\u03c1_l within {abs(_gap):.1f} % of the flow oil \u2014 the phase "
                        f"split differs, not the oil." if _agree else
                        f"the FLOW model runs {flow['gas_pct']:.0f} % gas \u2014 not the "
                        f"same fluid."))
        elif _agree:
            #  WHETHER THE PHASE SPLIT STILL DISAGREES IS A FACT ABOUT THIS RUN, not a
            #  standing caveat. This branch used to assert "the EOS holds the hydrocarbon
            #  single-phase to the bubble point" whenever the DENSITIES agreed, which was
            #  true only of the composition that was cut against density alone. Re-cutting
            #  it against phase behaviour as well put the line two-phase at every station,
            #  and the sentence then described a defect that no longer existed -- on the
            #  very figure whose job is to show it had been fixed.
            _two = _n_two >= 0.8 * max(len(recs), 1)
            _note += (f"  That is within {abs(_gap):.1f} % of the flow model's oil "
                      f"({_ro:.0f} kg/m\u00b3); its liquid reads {flow['rho_l']:.0f} kg/m\u00b3 "
                      f"only because it carries the {_wc:.0f} % water cut.")
            if _two:
                _note += (f"  The PHASE SPLIT agrees in character too: the EOS flashes "
                          f"two-phase at {_n_two} of {len(recs)} stations against the flow "
                          f"model's {flow['gas_pct']:.0f} % gas, so this panel and the flow "
                          f"figures describe one fluid.")
            else:
                _note += (f"  What differs is the PHASE SPLIT: the EOS holds the hydrocarbon "
                          f"single-phase over {len(recs) - _n_two} of {len(recs)} stations, "
                          f"while the flow model runs {flow['gas_pct']:.0f} % gas over the "
                          f"whole route from a fixed inlet gas rate.  Read the gas content "
                          f"here as EOS phase behaviour, not as the pipeline's.")
        else:
            _note += (f"  The FLOW model runs {flow['gas_pct']:.0f} % gas over the whole route "
                      f"on an oil of {_ro:.0f} kg/m\u00b3 ({_gap:+.0f} % against the EOS "
                      f"liquid).  These are not the same fluid \u2014 read this panel as EOS "
                      f"phase behaviour, not as the pipeline's gas content.")
    #  WRAP IT, AND RESERVE FOR WHAT IT ACTUALLY BECOMES. This note is built from the
    #  numbers, so its length is not fixed: on a case where the EOS and the flow oil agree
    #  it grows to explain the water cut and the phase split. Left as one fig.text line
    #  against a hardcoded 4.5 % bottom rect, a longer note is exactly what crushed a
    #  panel earlier in this audit -- the reserve has to follow the line count, not a
    #  number chosen when the note was short.
    _fw_in, _fh_in = (float(v) for v in fig.get_size_inches())
    _fs = 7.0 * (0.85 if _S.compact() else 1.0)
    _cpl = max(40, int(_fw_in * 72.0 / (_fs * 0.52)))          # chars that fit on a line
    _lines = textwrap.wrap(_note, width=_cpl) or [_note]
    #  GROW THE CANVAS, do not eat the axes. Reserving the caption's height out of a fixed
    #  figure is what put the legend on top of the x-label at the 0.30 slide scale: the
    #  note is the same size in points at every scale, so on the smallest canvas it claims
    #  a far bigger share of it. Adding the strip instead keeps every panel the size the
    #  scale asked for, and the bottom rect is then just that strip.
    _need_in = len(_lines) * _fs * 1.45 / 72.0 + 0.10
    fig.set_size_inches(_fw_in, _fh_in + _need_in, forward=True)
    _band = _need_in / (_fh_in + _need_in)
    fig.text(0.5, 0.012, "\n".join(_lines), ha="center", va="bottom", fontsize=_fs,
             color="#8A4B2A", style="italic", linespacing=1.35)
    fig.tight_layout(rect=(0, _band, 1, 0.97))
    _p = os.path.join(outdir, "compo_pvt.png")
    __import__("shct_style").screen(fig, _p)
    fig.savefig(_p, dpi=_FIG_DPI); plt.close(fig)
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

    def _f(v):
        """A blank cell means the phase does not exist at that station, not zero.

        The writer's _fmt() emits "" for any quantity a station has no phase for, so on
        this case 39 of 40 rows carry blank gas columns. This reader used to call float()
        on them straight, which made replot_from_csv() raise ValueError on the project's
        own compositional CSV -- a documented redraw path that could not run on the data
        it ships with. NaN is what the live path holds for the same stations, so the two
        now build identical records.
        """
        v = (v or "").strip()
        if not v:
            return float("nan")
        try:
            return float(v)
        except ValueError:
            return float("nan")

    recs = [{"x_km": _f(r["x_km"]), "P": _f(r["P_bar"]), "T": _f(r["T_C"]),
                 "V": _f(r["vapour_frac_V"]), "rho_g": _f(r["rho_gas_kgm3"]),
                 "rho_l": _f(r["rho_liq_kgm3"]), "mu_g": _f(r["mu_gas_Pas"]),
                 "mu_l": _f(r["mu_liq_Pas"]), "Z": _f(r["Z_gas"]), "sg": _f(r["gas_sg"]),
                 "K": {n: _f(r[f"K_{n}"]) for n in names}} for r in rows]
    #  the flow basis the live path wrote next to this CSV; without it the redraw would
    #  quietly lose the flow-vs-EOS comparison the caption is built on
    _flow = None
    _bp = os.path.join(os.path.dirname(os.path.abspath(csv_path)), "compo_flow_basis.json")
    if os.path.exists(_bp):
        try:
            with open(_bp) as _bh:
                _flow = json.load(_bh)
        except (OSError, ValueError):
            _flow = None
    return _plot_pvt(recs, names, outdir, flow=_flow)
