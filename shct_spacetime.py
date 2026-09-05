#!/usr/bin/env python3
# =============================================================================
#  shct_spacetime.py — the PUBLISHED-SCHEME space-time / multi-time figure set.
# -----------------------------------------------------------------------------
#  The flow-assurance and transient-multiphase literature presents a transient
#  pipeline calculation through a small, stable family of figures.  This module
#  renders the SHCT case study in exactly that family, so the case study can be
#  read side-by-side with the published work it should be compared against:
#
#    14  liquid holdup vs distance at successive times (early / late panels)
#    15  slug growth & propagation — stacked snapshots with tracked fronts
#    16  slug-train space-time waterfall + celerity scan + moveout-corrected stack
#    17  in-pipe volume fractions along the line + phase rates into the host
#    18  shut-in / late-time P, T, T_eq, water-holdup profile + deposit growth
#    19  TRUE space-time fields — holdup, pressure, velocity, subcooling (2x2)
#    20  liquid holdup along the line after successive shut-in / production times
#    21  riser depth-time waterfall with slug boundaries and slug lengths
#    22  pipeline cloud maps — phase distribution + temperature at successive times
#    23  DTS-style thermal waterfall — T(x,t) + pressure trace + stage annotations
#    24  temperature-gradient waterfall — dT/dx(x,t), which localises the front
#    25  DAS-style flow-noise waterfall — holdup fluctuation rate |dα/dt|(x,t)
#    26  along-pipeline parameter diagram — P, T, holdup, velocity at successive times
#    27  well-posedness / Kelvin-Helmholtz map with the case's operating states
#
#  WHAT IS SOLVER OUTPUT AND WHAT IS RECONSTRUCTION
#  ------------------------------------------------
#  Figures 14, 17, 18, 19, 20 and 22 are drawn DIRECTLY from the solver's
#  space-time snapshot history (alpha_l, p, T, j, delta, Tsub, regime, f_slug).
#
#  Figures 15, 16 and 21 show INDIVIDUAL slug units.  The transport grid is
#  dx ~ 460 m while a slug unit is ~10-40 m, so individual slugs are a SUB-GRID
#  quantity that the solver carries statistically (slug frequency f_slug, slug
#  unit length L_u = V_t/f_slug, slug-body holdup alpha_ls).  These three figures
#  therefore render a KINEMATIC RECONSTRUCTION built entirely from those solver
#  outputs: at every station the reconstructed square wave has the solver's local
#  slug frequency, the solver's local translational celerity V_t, and a body /
#  film split chosen so that the unit-averaged holdup reproduces the solver's
#  cell-average alpha_l EXACTLY (mass-consistent by construction).  Nothing is
#  invented: period, celerity, length and holdup all come from the run.  Every
#  such figure says so on its face.
#
#  A slug train only exists while the line is actually flowing intermittently, so
#  those three figures are built from the latest snapshot at which it is (for a
#  shut-in, that is a state before the line stops, and the time is printed on the
#  figure), and are skipped with a stated reason when no such state exists.
#
#  RE-RENDERING WITHOUT RE-RUNNING
#  -------------------------------
#  spacetime_outputs() also archives the snapshot history it reads as
#  spacetime_state.npz next to the figures.  rerender(folder) rebuilds the whole
#  set from that archive, so a figure can be restyled, rescaled or re-cropped
#  without repeating a multi-hour transient.
#
#  STYLE
#  -----
#  The rendering follows the published scheme (full box frame, light dotted grid,
#  plain descriptive titles, per-panel colourbars on the space-time fields,
#  filled contours with thin overlaid contour lines) while keeping this project's
#  two standing style rules: NO black and NO dark colours anywhere, and legends
#  placed OUTSIDE the axes so that no label ever covers a curve or a band.
#
#  Colour alone never carries a series: a line that has reached a steady profile
#  draws every time level on top of the others, so each also gets its own dash
#  pattern and marker, the markers are staggered along x, and the panel says when
#  the curves coincide rather than leaving one visible colour to explain itself.
#
#  Author: Akosa Samuel Onyejekwe.
# =============================================================================
from __future__ import annotations

import math
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, ListedColormap, BoundaryNorm

import shct_style as S

S.apply_style()

_DPI = int(os.environ.get("SHCT_FIG_DPI", "320"))
_TITLES = os.environ.get("SHCT_FIG_TITLES", "1") not in ("0", "false", "False")

#  the ordered, saturated, never-dark line palette used for the multi-time
#  overlays (the published figures use one strong hue per time level)
TIME_COLORS = [S.BLUE, S.RED, S.GREEN, S.PURPLE, S.ORANGE, S.TEAL,
               S.MAGENTA, S.AMBER, S.SKY, S.BROWN]

#  thin contour lines drawn over the filled space-time fields.  The published
#  scheme draws these in black; a medium periwinkle keeps the same read without
#  breaking the no-black / no-dark rule.
CONTOUR_LINE = "#7E93D6"

G = 9.80665


# =============================================================================
#  small shared helpers
# =============================================================================
def _frame(ax, grid=True, minor=False):
    """The published look: a full four-sided box, thin outward ticks, light grid."""
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
        sp.set_color(S.INK)
        sp.set_visible(True)
    ax.tick_params(direction="out", length=3.4, width=0.8, labelsize=8)
    if minor:
        ax.minorticks_on()
        ax.tick_params(which="minor", direction="out", length=1.8, width=0.6)
    if grid:
        ax.grid(True, color=S.GRIDC, lw=0.6, ls=":", alpha=0.95)
        ax.set_axisbelow(True)
    return ax


def _title(ax, text, size=9.5):
    if _TITLES:
        ax.set_title(text, color=S.TITLE, fontweight="bold", fontsize=size, pad=6)
    return ax


def _legend(ax, ncol=1, size=7.5, anchor=(1.012, 1.0), handles=None, title=None):
    """Legend OUTSIDE the axes (standing project rule: never over the data)."""
    kw = dict(loc="upper left", bbox_to_anchor=anchor, borderaxespad=0.0,
              fontsize=size, ncol=ncol, framealpha=1.0, facecolor="white",
              edgecolor=S.INK, fancybox=True)
    leg = ax.legend(handles=handles, title=title, **kw) if handles is not None \
        else ax.legend(title=title, **kw)
    if leg and leg.get_title() is not None:
        leg.get_title().set_fontsize(size)
        leg.get_title().set_color(S.INK)
    if leg:
        leg.get_frame().set_linewidth(0.8)
    return leg


#  Series in a multi-time overlay routinely coincide — a line that has reached a
#  steady profile draws six curves on top of one another, and only the last colour
#  is then visible.  Colour alone therefore cannot carry the series identity: each
#  time level also gets its own dash pattern and its own marker, the markers are
#  staggered along x so overlapping curves still show separately, and the line
#  weight decreases with time so an earlier, thicker curve is visible under a
#  later, thinner one.
_DASHES = [(None, None), (7, 2), (2, 2), (9, 2, 2, 2), (4, 1, 1, 1), (1, 1.4),
           (12, 3), (5, 2, 1, 2), (3, 3), (10, 2, 1, 2)]
_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]


def _series_style(i, n, nx):
    """(linewidth, dashes, marker, markevery) for series i of n over nx points."""
    lw = 2.5 - 1.5 * (i / max(n - 1, 1))
    dash = _DASHES[i % len(_DASHES)]
    step = max(4, nx // 9)
    return (lw, dash, _MARKERS[i % len(_MARKERS)],
            (int(round(i * step / max(n, 1))) % step, step))


def _plot_series(ax, x, y, i, n, color, label):
    lw, dash, marker, mev = _series_style(i, n, len(x))
    ln, = ax.plot(x, y, color=color, lw=lw, label=label, marker=marker,
                  markevery=mev, markersize=4.6, markeredgewidth=0.8,
                  markerfacecolor="white", markeredgecolor=color,
                  solid_capstyle="round")
    if dash[0] is not None:
        ln.set_dashes(list(dash))
    return ln


def _coincidence_note(ax, curves, what="profiles"):
    """If the series lie on top of one another, say so instead of leaving the
    reader to wonder why only one curve is visible."""
    A = np.asarray(curves, float)
    if A.ndim != 2 or A.shape[0] < 2:
        return
    spread = float(np.nanmax(np.nanmax(A, 0) - np.nanmin(A, 0)))
    if spread < 0.02 * max(float(np.nanmax(np.abs(A))), 1e-9) or spread < 0.01:
        ax.annotate(f"the {what} coincide to within {spread:.3g} — the line holds "
                    f"a steady profile over this window",
                    xy=(0.5, -0.16), xycoords="axes fraction", ha="center",
                    va="top", fontsize=7.2, style="italic", color=S.INK,
                    annotation_clip=False,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=S.GRIDC,
                              lw=0.7))


def _margin_note(ax, y_data, text, side="right", color=None, pad=0.02,
                 size=8.5, leader=True):
    """Label a feature from OUTSIDE the axes, with a leader line back to it.

    On a colour-mapped field every pixel is a value, so a label placed inside the
    axes hides part of the result. This puts the text in the margin, aligned with
    the feature's own coordinate, and draws a short leader to it.
    """
    color = color or S.INK
    if side == "right":
        x_text, ha, xy_ax = 1.0 + pad, "left", 1.0
    else:
        x_text, ha, xy_ax = -pad, "right", 0.0
    trans = ax.get_yaxis_transform()          # x in axes fraction, y in data
    ax.annotate(text, xy=(xy_ax, y_data), xycoords=trans,
                xytext=(x_text, y_data), textcoords=trans,
                ha=ha, va="center", fontsize=size, fontweight="bold",
                color=color, annotation_clip=False,
                arrowprops=(dict(arrowstyle="-", color=color, lw=0.9,
                                 shrinkA=0, shrinkB=0) if leader else None),
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=S.GRIDC,
                          lw=0.8))


def _stage_header(ax, stages, size=9.0):
    """Name the operating stages ABOVE the axes, as a timeline header, so the
    field itself is never written over.

    A run with a SINGLE stage needs no header: the figure title already names the
    scenario, and a lone centred label would only collide with it. Returns the
    extra title padding the caller must leave when a header was drawn.
    """
    if len(stages) < 2:
        return 0.0
    trans = ax.get_xaxis_transform()          # x in data, y in axes fraction
    for i, (ta, tb, lab) in enumerate(stages):
        ax.annotate(lab, xy=(0.5 * (ta + tb), 1.012), xycoords=trans,
                    ha="center", va="bottom", fontsize=size, fontweight="bold",
                    color=S.INK, annotation_clip=False,
                    bbox=dict(boxstyle="round,pad=0.28", fc="white",
                              ec=S.GRIDC, lw=0.8))
        if i:
            ax.axvline(ta, color="white", lw=1.6, ls="--")
    return 20.0                                # points of title pad to leave


def _save(fig, path, check=True):
    """Save a figure, having first checked that no text overlaps any other text."""
    if check:
        S.report_text_overlaps(fig, os.path.basename(path))
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    return path


def _tlab(t_h):
    """Time label in the published style: minutes while short, hours once long."""
    t_h = float(t_h)
    if t_h < 1.0:
        return f"{t_h * 60.0:.0f} min"
    if t_h < 10.0:
        return f"{t_h:.1f} h"
    return f"{t_h:.0f} h"


def _pick(snap_t, n, frac_lo=0.06, frac_hi=1.0):
    """n snapshot indices spread evenly over [frac_lo, frac_hi] of the run."""
    snap_t = np.asarray(snap_t, float)
    m = snap_t.size
    if m == 0:
        return np.zeros(0, int)
    lo = int(round(frac_lo * (m - 1)))
    hi = int(round(frac_hi * (m - 1)))
    lo = max(0, min(lo, m - 1))
    hi = max(lo, min(hi, m - 1))
    if n >= (hi - lo + 1):
        return np.arange(lo, hi + 1)
    return np.unique(np.linspace(lo, hi, n).round().astype(int))


def _scenario_label(sv):
    kind = getattr(sv.case.scenario, "kind", "steady")
    meg = float(getattr(sv.case.operating, "MEG_wt_inlet", 0.0) or 0.0)
    if kind == "shutin":
        return "unplanned shut-in"
    return "engineered fix (insulation + MEG)" if meg > 0 else "as-operated"


def _theta(sv):
    """Local inclination [rad] of every cell from the elevation profile."""
    z = np.asarray(sv.z, float)
    x = np.asarray(sv.x, float)
    dz = np.gradient(z, x)
    return np.arctan(dz)


def _event_time_h(sv):
    """Time of the scenario event (shut-in), or None for a plain production run."""
    sc = sv.case.scenario
    if getattr(sc, "kind", "steady") == "shutin":
        return float(getattr(sc, "event_time_h", 0.0))
    return None


def _slug_snapshot(sv):
    """The snapshot at which to build the resolved-slug figures, or None.

    A slug train only exists while the line is actually flowing intermittently.
    In a shut-in the flow stops, the slug frequency falls to its floor and the
    "slug unit length" grows past the length of the pipe, so a reconstruction at
    the final state would be meaningless.  This scores every snapshot by how much
    intermittent, flowing, liquid-carrying pipe it has and returns the LATEST
    snapshot that is still within 20 % of the best one — a developed, actively
    slugging state — or None if the line never slugs.
    """
    r = sv.results
    reg = np.asarray(r.get("snap_regime", np.empty(0)), float)
    fs = np.asarray(r.get("snap_fslug", np.empty(0)), float)
    al = np.asarray(r.get("snap_holdup", np.empty(0)), float)
    j = np.asarray(r.get("snap_j", np.empty(0)), float)
    if any(a.ndim != 2 or a.size == 0 for a in (reg, fs, al, j)):
        return None
    x = np.asarray(sv.x, float)
    flowline = (x < 0.90 * x[-1])[None, :]
    inter = np.isin(np.round(reg), [2, 5]) & (j > 0.05) & flowline
    score = np.where(inter, fs * al, 0.0).sum(axis=1)
    if not np.isfinite(score).any() or float(np.nanmax(score)) <= 0.0:
        return None
    good = np.where(score >= 0.8 * float(np.nanmax(score)))[0]
    return int(good[-1]) if good.size else int(np.nanargmax(score))


def _riser_slug_snapshot(sv, i0):
    """The snapshot at which the RISER is most actively slugging, or None.

    _slug_snapshot() scores `x < 0.90 * x[-1]` — the flowline, deliberately
    excluding the riser — which is right for the flowline figures and wrong for
    this one. Choosing the instant by flowline activity and then asking whether
    the riser has a slug train at that same instant conflates two different
    questions, and the riser figure was being rejected for a state it had never
    been asked about. Score the riser on its own terms instead.
    """
    r = sv.results
    reg = np.asarray(r.get("snap_regime", np.empty(0)), float)
    fs = np.asarray(r.get("snap_fslug", np.empty(0)), float)
    al = np.asarray(r.get("snap_holdup", np.empty(0)), float)
    j = np.asarray(r.get("snap_j", np.empty(0)), float)
    if any(a.ndim != 2 or a.size == 0 for a in (reg, fs, al, j)):
        return None
    riser = np.zeros(np.asarray(sv.x, float).size, bool)
    riser[i0:] = True
    inter = np.isin(np.round(reg), [2, 5]) & (j > 0.05) & riser[None, :]
    score = np.where(inter, fs * al, 0.0).sum(axis=1)
    if not np.isfinite(score).any() or float(np.nanmax(score)) <= 0.0:
        return None
    good = np.where(score >= 0.8 * float(np.nanmax(score)))[0]
    return int(good[-1]) if good.size else int(np.nanargmax(score))


def _slug_train_ok(sv, F, ic):
    """Is a resolved slug train physically meaningful at this cell?"""
    route = float(np.asarray(sv.x, float)[-1])
    Lu = float(F["Lu"][ic])
    if not np.isfinite(Lu) or Lu <= 0.0:
        return False, "slug unit length is not finite"
    D_pipe = float(sv.case.pipeline.diameter_m)
    if Lu > 0.02 * route or Lu > 400.0 * D_pipe:
        return False, (f"slug unit length {Lu:.0f} m is not a slug scale (over 2 % of "
                       f"the route, or over 400 pipe diameters) — the line is not "
                       f"slugging at this state")
    if float(F["Vt"][ic]) < 0.2:
        return False, "translational celerity is below 0.2 m/s — no travelling train"
    return True, ""


def _slugging_cell(sv, F, flowline_only=True):
    """Index of the most actively slugging cell: intermittent regime, high f_slug
    and high holdup.  This is the reach the resolved-slug figures zoom into.  The
    riser is excluded by default because it has its own depth-time figure, and
    because f_slug varies steeply there."""
    r = sv.results
    x = np.asarray(sv.x, float)
    reg = np.asarray(r.get("snap_regime", np.empty(0)), float)
    reg = reg[-1] if (reg.ndim == 2 and reg.size) else np.nanmedian(
        np.asarray(r["regime"], float), 1)
    inter = np.isin(np.round(reg), [2, 5])
    score = np.where(inter, F["fslug"] * F["alpha"], 0.0)
    if flowline_only:
        score = np.where(x < 0.90 * x[-1], score, 0.0)
    if not np.isfinite(score).any() or float(np.nanmax(score)) <= 0.0:
        return int(x.size // 2)
    return int(np.nanargmax(score))


# =============================================================================
#  the mass-consistent sub-grid slug-unit reconstruction (figs 15, 16, 21)
# =============================================================================
def slug_unit_fields(sv, k_snap=-1):
    """Per-cell slug-unit descriptors at snapshot k_snap, all from solver output.

    Returns dict with (per cell):
        alpha   solver cell-average liquid holdup            [-]
        fslug   solver slug frequency                        [Hz]
        Vt      slug translational celerity 1.2*Vm + drift   [m/s]
        Lu      slug unit length  Vt / fslug                 [m]
        als     slug-body holdup (Gregory et al. 1978)       [-]
        alf     residual film holdup                         [-]
        beta    slug-body length fraction L_s / L_u          [-]
        Ls      slug body length                             [m]

    beta is solved from the unit-average mass balance
        alpha = beta*als + (1-beta)*alf
    so the reconstructed square wave integrates back to the solver's alpha
    exactly.  alf is the residual film holdup left behind a passing slug.
    """
    import solver as _solver

    r = sv.results
    D = float(sv.case.pipeline.diameter_m)

    def _snap(key, fallback):
        A = np.asarray(r.get(key, np.empty(0)), float)
        if A.ndim == 2 and A.shape[0] > 0:
            return A[k_snap].copy()
        return np.nanmedian(np.asarray(r[fallback], float), 1)

    alpha = np.clip(_snap("snap_holdup", "alpha_l"), 1e-3, 0.999)
    fslug = np.maximum(_snap("snap_fslug", "fslug"), 1e-4)
    j = np.maximum(_snap("snap_j", "j"), 1e-3)

    Vt = 1.2 * j + 0.35 * math.sqrt(G * D)
    Lu = np.clip(Vt / fslug, D, 5000.0)
    als = _solver.slug_body_holdup(j)

    #  residual film left behind the slug: a thin drained film, capped so that a
    #  body fraction always exists.  The published snapshots (holdup ~1 in the
    #  body over a ~0.03-0.10 film) sit in exactly this range.
    alf = np.minimum(0.05, 0.5 * alpha)
    als = np.maximum(als, alpha + 1e-3)                 # keep beta <= 1
    beta = np.clip((alpha - alf) / np.maximum(als - alf, 1e-6), 0.02, 0.98)
    Ls = beta * Lu
    return dict(alpha=alpha, fslug=fslug, Vt=Vt, Lu=Lu, als=als, alf=alf,
                beta=beta, Ls=Ls)


def reconstruct_slug_field(sv, xq_m, tq_s, k_snap=-1, t0_s=0.0):
    """Reconstructed holdup field alpha_l(x, t) resolving individual slug units.

    The phase of the slug train at (x, t) is

        Theta(x, t) = f_slug(x) * ( t - tau(x) ),      tau(x) = int_0^x dx'/V_t(x')

    so that at any station the passage period is the solver's 1/f_slug, the
    front trajectories in the (x, t) plane have the solver's celerity V_t, and
    the spatial wavelength is the solver's slug unit length L_u = V_t/f_slug.
    A cell is inside a slug body while frac(Theta) < beta.

    Returns (field, meta) with field of shape (len(tq_s), len(xq_m)).
    """
    F = slug_unit_fields(sv, k_snap=k_snap)
    x = np.asarray(sv.x, float)

    #  cumulative transit time tau(x) on the solver grid, then interpolated
    tau = np.concatenate([[0.0], np.cumsum(np.diff(x) / np.maximum(
        0.5 * (F["Vt"][:-1] + F["Vt"][1:]), 1e-3))])

    xq = np.asarray(xq_m, float)
    itp = lambda A: np.interp(xq, x, A)
    #  reference the transit time to the UPSTREAM END OF THE WINDOW.  Referencing
    #  it to the inlet makes tau ~ 1e4 s, so a 1 % along-line variation of f_slug
    #  would swing the phase by ~100 cycles and alias the whole window; relative
    #  to the window the phase is well conditioned and the train reads correctly.
    tauq = np.interp(xq, x, tau)
    tauq = tauq - tauq[0]
    fq, bq, alsq, alfq = (itp(F["fslug"]), itp(F["beta"]),
                          itp(F["als"]), itp(F["alf"]))
    Vtq, Luq = itp(F["Vt"]), itp(F["Lu"])

    T = np.asarray(tq_s, float)[:, None] + float(t0_s)
    theta = fq[None, :] * (T - tauq[None, :])
    phase = theta - np.floor(theta)
    body = phase < bq[None, :]
    field = np.where(body, alsq[None, :], alfq[None, :])
    meta = dict(Vt=Vtq, Lu=Luq, fslug=fq, beta=bq, tau=tauq, x=xq)
    return field, meta


# =============================================================================
#  14 — liquid holdup vs distance at successive times   (scheme: A / C / J)
# =============================================================================
def fig_holdup_multitime(sv, outdir):
    r = sv.results
    H, ts = np.asarray(r["snap_holdup"], float), np.asarray(r["snap_t"], float)
    if H.size == 0:
        return None
    x = sv.x / 1000.0
    lab = _scenario_label(sv)

    fig, ax = plt.subplots(2, 1, figsize=(8.4, 7.0), sharex=True)
    halves = [_pick(ts, 6, 0.06, 0.50), _pick(ts, 6, 0.54, 1.00)]
    names = ["early transient", "late transient / quasi-developed"]
    for a, idx, nm in zip(ax, halves, names):
        for i, (c, k) in enumerate(zip(TIME_COLORS, idx)):
            _plot_series(a, x, H[k], i, len(idx), c, _tlab(ts[k]))
        _coincidence_note(a, [H[k] for k in idx], "holdup profiles")
        a.set_ylabel("liquid holdup  α_l  [-]", fontsize=9)
        a.set_ylim(0, 1.02)
        a.set_xlim(x.min(), x.max())
        _frame(a, minor=True)
        _legend(a, size=7.5, title=nm)
    ax[1].set_xlabel("distance from wellhead  [km]", fontsize=9)
    _title(ax[0], f"Liquid-holdup evolution along the 32 km tie-back — {lab}")

    fig.tight_layout()
    p = os.path.join(outdir, "14_holdup_multitime.png")
    return _save(fig, p)


# =============================================================================
#  15 — slug growth & propagation, stacked snapshots       (scheme: D)
# =============================================================================
def fig_slug_growth(sv, outdir):
    """Three successive snapshots of the resolved slug train over a short reach,
    with ONE front tracked across the panels so its arrival time and position are
    annotated exactly as the published 'slug growth & propagation' figure does."""
    k = _slug_snapshot(sv)
    if k is None:
        print("    [space-time] 15_slug_growth_propagation skipped: the line never "
              "reaches an intermittent, flowing state in this scenario", flush=True)
        return None
    F = slug_unit_fields(sv, k_snap=k)
    x = np.asarray(sv.x, float)
    ic = _slugging_cell(sv, F)
    ok, why = _slug_train_ok(sv, F, ic)
    if not ok:
        print(f"    [space-time] 15_slug_growth_propagation skipped: {why}", flush=True)
        return None
    t_snap = float(np.asarray(sv.results["snap_t"], float)[k])

    Vt_c = float(F["Vt"][ic])
    f_c = float(F["fslug"][ic])
    Lu_c = float(F["Lu"][ic])

    #  a window ~8 slug units long: long enough to show growth and coalescence,
    #  short enough that the individual units are resolved on the page
    window_m = float(np.clip(8.0 * Lu_c, 40.0, 0.05 * x[-1]))
    x0 = float(np.clip(x[ic] - 0.5 * window_m, x[0], x[-1] - window_m))
    xq = np.linspace(x0, x0 + window_m, 2400)

    #  advance by 0.45 of a slug period between panels so the train visibly moves
    dt = 0.45 / max(f_c, 1e-4)
    times = np.array([0.0, dt, 2.0 * dt])

    fig, ax = plt.subplots(3, 1, figsize=(8.6, 7.4))
    x_track = None
    for i, (a, tt) in enumerate(zip(ax, times)):
        fld, meta = reconstruct_slug_field(sv, xq, np.array([0.0]), k_snap=k,
                                           t0_s=tt)
        y = fld[0]
        a.plot(xq - x0, y, color=S.RED, lw=1.5)
        a.set_ylim(0, 1.12)
        a.set_xlim(0, window_m)
        a.set_ylabel("liquid holdup", fontsize=8.5)
        _frame(a, minor=True)
        a.set_title(f"({_scenario_label(sv)}), time = {tt:.1f} s   "
                    f"[reach {x0/1000:.2f}–{(x0+window_m)/1000:.2f} km]",
                    fontsize=8.5, color=S.INK, pad=4)

        #  follow ONE front: pick the front nearest the window centre in panel 1,
        #  then advance it at the solver's own celerity
        rises = np.where(np.diff(y) > 0.25)[0]
        if i == 0 and rises.size:
            x_track = float(xq[rises[np.argmin(np.abs(xq[rises] - (x0 + 0.35 * window_m)))]])
        if x_track is not None:
            xb = x_track + Vt_c * tt - x0
            if 0.0 <= xb <= window_m:
                a.axvline(xb, color=S.BLUE, lw=1.1, ls="--", alpha=0.9)
                a.plot([xb], [1.04], marker="v", ms=6, color=S.BLUE)
            a.annotate(f"T$_{{b{i+1}}}$ ~ {tt:.1f} s\nX$_{{b{i+1}}}$ ~ {xb:.1f} m",
                       xy=(1.012, 0.5), xycoords="axes fraction",
                       ha="left", va="center", fontsize=8.5, color=S.BLUE,
                       fontweight="bold", annotation_clip=False,
                       bbox=dict(boxstyle="round,pad=0.28", fc="white",
                                 ec=S.GRIDC, lw=0.7))
    ax[-1].set_xlabel("distance within the reach  [m]", fontsize=9)
    if _TITLES:
        fig.suptitle(f"Slug propagation and front tracking — resolved slug units "
                     f"at t = {t_snap:.1f} h (V$_t$ = {Vt_c:.2f} m/s, "
                     f"f$_{{slug}}$ = {f_c:.2f} Hz, L$_u$ = {Lu_c:.1f} m)",
                     color=S.TITLE, fontweight="bold", fontsize=10, y=0.995)
    fig.text(0.5, 0.005,
             "Mass-consistent sub-grid reconstruction: period, celerity, unit length and "
             "unit-average holdup are the solver's own fields.",
             ha="center", fontsize=6.8, color=S.INK, style="italic")
    fig.tight_layout(rect=(0, 0.022, 1, 0.975))
    p = os.path.join(outdir, "15_slug_growth_propagation.png")
    return _save(fig, p)


# =============================================================================
#  16 — slug-train space-time waterfall + celerity scan     (scheme: E)
# =============================================================================
def _moveout(fluc, tq, shift_s):
    """Vectorised linear-moveout: sample each station's trace at t + shift[station]."""
    nt, nx = fluc.shape
    dt = float(tq[1] - tq[0])
    base = np.arange(nt)[:, None] + (np.asarray(shift_s, float)[None, :] / dt)
    i0 = np.floor(base).astype(int)
    w = base - i0
    ok = (i0 >= 0) & (i0 + 1 < nt)
    i0c = np.clip(i0, 0, nt - 2)
    cols = np.arange(nx)[None, :]
    out = fluc[i0c, cols] * (1.0 - w) + fluc[i0c + 1, cols] * w
    out[~ok] = np.nan
    return out


def _single_slug_field(sv, xq, tq, k_snap=-1):
    """The space-time signature of ONE slug unit crossing the reach.

    Same kinematics as the train (phase Theta = f_slug*(t - tau(x))) but only the
    k = 0 unit is filled, so the figure carries a single coherent event — the
    condition under which a moveout / semblance scan has one unambiguous peak.
    """
    F = slug_unit_fields(sv, k_snap=k_snap)
    x = np.asarray(sv.x, float)
    tau = np.concatenate([[0.0], np.cumsum(np.diff(x) / np.maximum(
        0.5 * (F["Vt"][:-1] + F["Vt"][1:]), 1e-3))])
    itp = lambda A: np.interp(xq, x, A)
    fq, tauq = itp(F["fslug"]), np.interp(xq, x, tau)
    bq, alsq, alfq = itp(F["beta"]), itp(F["als"]), itp(F["alf"])
    theta = fq[None, :] * (np.asarray(tq, float)[:, None] - (tauq - tauq[0])[None, :])
    body = (theta >= 0.0) & (theta < bq[None, :])
    return np.where(body, alsq[None, :], alfq[None, :]), itp(F["Vt"]), itp(F["Lu"])


def fig_slug_waterfall(sv, outdir):
    """Track one slug across a reach in the space-time plane: waterfall, celerity
    scan, moveout-corrected waterfall and the distance-stacked trace whose width
    gives the slug body length — the published slug-tracking quality-check set."""
    k = _slug_snapshot(sv)
    if k is None:
        print("    [space-time] 16_slug_train_waterfall skipped: the line never "
              "reaches an intermittent, flowing state in this scenario", flush=True)
        return None
    F = slug_unit_fields(sv, k_snap=k)
    x = np.asarray(sv.x, float)
    ic = _slugging_cell(sv, F)
    ok, why = _slug_train_ok(sv, F, ic)
    if not ok:
        print(f"    [space-time] 16_slug_train_waterfall skipped: {why}", flush=True)
        return None
    t_snap = float(np.asarray(sv.results["snap_t"], float)[k])
    Vt_c, Lu_c, f_c = float(F["Vt"][ic]), float(F["Lu"][ic]), float(F["fslug"][ic])
    Ls_c = float(F["beta"][ic]) * Lu_c

    reach = float(np.clip(12.0 * Lu_c, 50.0, 0.05 * x[-1]))
    x0 = float(np.clip(x[ic] - 0.5 * reach, x[0], x[-1] - reach))
    xq = np.linspace(x0, x0 + reach, 420)
    t_transit = reach / max(Vt_c, 1e-3)
    window_s = 1.7 * t_transit
    tq = np.linspace(0.0, window_s, 700)

    fld, Vtq, Luq = _single_slug_field(sv, xq, tq, k_snap=k)
    fluc = fld - fld.mean(axis=0, keepdims=True)
    dxs = xq - xq[0]

    #  (b) semblance over trial celerities — one event, so one peak
    vs = np.linspace(0.45 * Vt_c, 2.0 * Vt_c, 130)
    sem = np.zeros_like(vs)
    for i, v in enumerate(vs):
        stk = _moveout(fluc, tq, dxs / v)
        ok = np.isfinite(stk).all(axis=1)
        if ok.sum() < 10:
            continue
        srow = stk[ok]
        den = srow.shape[1] * np.sum(srow ** 2)
        sem[i] = float(np.sum(np.sum(srow, axis=1) ** 2) / den) if den > 0 else 0.0
    v_best = float(vs[int(np.argmax(sem))])

    corr = _moveout(fluc, tq, dxs / v_best)
    ok = np.isfinite(corr).all(axis=1)
    stacked = np.full(tq.size, np.nan)
    stacked[ok] = np.nanmean(corr[ok], axis=1)

    t_valid = max(window_s - reach / max(v_best, 1e-3), 0.15 * window_s)
    fig, ax = plt.subplots(1, 4, figsize=(13.2, 4.4),
                           gridspec_kw=dict(width_ratios=[1.18, 1.0, 1.18, 1.0]))
    vlim = float(np.nanmax(np.abs(fluc))) or 1.0
    for a, Fld, ttl in ((ax[0], fluc, "(a) slug waterfall"),
                        (ax[2], corr, f"(c) after moveout at {v_best:.2f} m/s")):
        pcm = a.pcolormesh(tq, dxs, np.ma.masked_invalid(Fld).T, cmap="shct_seq",
                           shading="gouraud", vmin=-vlim, vmax=vlim)
        a.set_xlabel("time  [s]", fontsize=8.5)
        a.set_ylabel("distance along reach  [m]", fontsize=8.5)
        a.set_xlim(tq.min(), tq.max() if a is ax[0] else t_valid)
        a.set_ylim(0, reach)
        a.set_title(ttl, fontsize=8.5, color=S.INK, pad=4)
        _frame(a, grid=False)
        cb = fig.colorbar(pcm, ax=a, pad=0.02, fraction=0.05)
        cb.set_label("α'$_l$  [-]", fontsize=7.5)
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_edgecolor(S.INK)

    ax[1].plot(vs, sem, color=S.TEAL, lw=1.5, marker="x", ms=3.0, mew=0.8)
    ax[1].axvline(v_best, color=S.MAGENTA, lw=1.2, ls="--")
    ax[1].axvline(Vt_c, color=S.ORANGE, lw=1.0, ls=":")
    ax[1].plot([], [], color=S.MAGENTA, ls="--", lw=1.2,
               label=f"recovered {v_best:.2f} m/s")
    ax[1].plot([], [], color=S.ORANGE, ls=":", lw=1.0,
               label=f"solver V$_t$ {Vt_c:.2f} m/s")
    ax[1].set_xlabel("trial celerity  [m/s]", fontsize=8.5)
    ax[1].set_ylabel("semblance  [-]", fontsize=8.5)
    ax[1].set_title("(b) semblance vs celerity", fontsize=8.5, color=S.INK, pad=4)
    _frame(ax[1], minor=True)
    _legend(ax[1], size=7.0)

    ax[3].plot(tq, stacked, color=S.BLUE, lw=1.3)
    ax[3].set_xlabel("time  [s]", fontsize=8.5)
    ax[3].set_ylabel("stacked α'$_l$  [-]", fontsize=8.5)
    ax[3].set_xlim(0.0, t_valid)
    ax[3].set_title("(d) distance-stacked trace", fontsize=8.5, color=S.INK, pad=4)
    _frame(ax[3], minor=True)
    fin = np.isfinite(stacked)
    if fin.any():
        thr = 0.5 * float(np.nanmax(stacked))
        wide = tq[fin][stacked[fin] > thr]
        if wide.size > 1:
            ax[3].annotate("", xy=(wide[-1], thr), xytext=(wide[0], thr),
                           arrowprops=dict(arrowstyle="<->", color=S.RED, lw=1.3))
            #  panel (d) already carries a title at the top, so this goes in the
            #  free right margin — never stacked on the title.
            ax[3].annotate(f"L$_s$ ≈\n{Ls_c:.1f} m", xy=(1.03, 0.5),
                           xycoords="axes fraction", ha="left", va="center",
                           fontsize=8, color=S.RED, fontweight="bold",
                           annotation_clip=False,
                           bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                     ec=S.GRIDC, lw=0.7))

    if _TITLES:
        fig.suptitle(f"Slug tracking in the space-time plane — {_scenario_label(sv)}, "
                     f"reach {x0/1000:.2f}–{(x0+reach)/1000:.2f} km at "
                     f"t = {t_snap:.1f} h (L$_u$ = {Lu_c:.1f} m, "
                     f"f$_{{slug}}$ = {f_c:.2f} Hz)",
                     color=S.TITLE, fontweight="bold", fontsize=10, y=0.995)
    fig.text(0.5, 0.005,
             "One slug unit of the mass-consistent sub-grid reconstruction; the celerity "
             "recovered by the moveout scan (b) returns the solver's own V$_t$.",
             ha="center", fontsize=6.8, color=S.INK, style="italic")
    fig.tight_layout(rect=(0, 0.035, 1, 0.955))
    p = os.path.join(outdir, "16_slug_train_waterfall.png")
    return _save(fig, p)


# =============================================================================
#  17 — in-pipe volume fractions + phase rates into the host   (scheme: F)
# =============================================================================
def fig_hydrate_distribution(sv, eng, outdir):
    import solver as _solver

    r, c = sv.results, sv.case
    ts = np.asarray(r["snap_t"], float)
    if ts.size == 0:
        return None
    x = sv.x / 1000.0
    D = float(c.pipeline.diameter_m)
    WC = float(c.fluids.water_cut)

    k = -1
    t_show = float(ts[k])
    alpha = np.asarray(r["snap_holdup"], float)[k]
    delta = np.asarray(r["snap_delta"], float)[k]
    phi = np.asarray(r["snap_phi"], float)[k]

    #  annular wall deposit of thickness delta -> in-pipe volume fraction
    ratio = np.clip(2.0 * delta / D, 0.0, 1.0)
    phi_dep = 100.0 * (1.0 - (1.0 - ratio) ** 2)
    phi_hyd_liq = 100.0 * np.clip(phi, 0.0, 1.0) * alpha
    phi_water = 100.0 * alpha * WC * np.clip(1.0 - phi, 0.0, 1.0)

    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.2))

    a = ax[0]
    a.fill_between(x, 0, phi_dep, color=S.SKY, alpha=0.55, lw=0)
    a.plot(x, phi_dep, color=S.SKY, lw=1.4,
           label="hydrate deposit at the wall")
    a.plot(x, phi_water, color=S.TEAL, lw=1.6, label="unconverted water")
    a.plot(x, phi_hyd_liq, color=S.RED, lw=1.6, label="hydrate in the liquids")
    a.set_xlabel("distance from wellhead  [km]", fontsize=9)
    a.set_ylabel("volume fraction in pipe, φ  [vol %]", fontsize=9)
    a.set_xlim(x.min(), x.max())
    a.set_ylim(bottom=0)
    _frame(a, minor=True)
    a.set_title(f"Simulation time = {t_show:.1f} hours", fontsize=9, color=S.INK, pad=4)
    _legend(a, size=7.5)

    #  (b) phase mass rates delivered into the host separator vs time.  The phase
    #      velocities come straight from the solver's momentum/pressure solve, so
    #      the delivered rates are the run's own outlet fluxes, not a re-closure.
    A = math.pi * D ** 2 / 4.0
    al_o = np.clip(np.asarray(r["snap_holdup"], float)[:, -1], 1e-3, 0.999)
    P_o = np.asarray(r["snap_P"], float)[:, -1]
    T_o = np.asarray(r["snap_T"], float)[:, -1]
    vl_a = np.asarray(r.get("snap_vl", np.empty(0)), float)
    vg_a = np.asarray(r.get("snap_vg", np.empty(0)), float)
    if vl_a.ndim == 2 and vl_a.size:
        v_l = np.clip(vl_a[:, -1], 0.0, None)
        v_g = np.clip(vg_a[:, -1], 0.0, None)
    else:                                   # older result set: drift-flux fallback
        th = _theta(sv)[-1]
        C0 = 1.05 + 0.15 * math.sin(abs(th))
        vd = (0.35 * math.sqrt(G * D) * math.sin(th)
              + 0.20 * math.sqrt(G * D) * math.cos(th))
        j_o = np.maximum(np.asarray(r["snap_j"], float)[:, -1], 0.0)
        v_g = np.clip(C0 * j_o + vd, 0.0, None)
        v_l = np.clip((j_o - (1.0 - al_o) * v_g) / al_o, 0.0, None)
    try:
        rho_g = _solver.gas_density(P_o, T_o, c.fluids)
    except Exception:                       # a restored state without a full Fluids
        rho_g = (P_o * 1e5) / (8.3145 / 0.019 * (T_o + 273.15))
    m_gas = rho_g * (1.0 - al_o) * v_g * A
    m_oil = float(c.fluids.rho_oil) * (1.0 - WC) * al_o * v_l * A
    m_wat = float(c.fluids.rho_water) * WC * al_o * v_l * A

    b = ax[1]
    b.plot(ts, m_gas, color=S.BLUE, lw=1.6, label="gas")
    b.plot(ts, m_oil, color=S.ORANGE, lw=1.6, ls="--", label="oil (live crude)")
    b.plot(ts, m_wat, color=S.TEAL, lw=1.6, ls=":", label="water")
    b.set_xlabel("time  [h]", fontsize=9)
    b.set_ylabel("mass flow rate into separator, ṁ  [kg s$^{-1}$]", fontsize=9)
    b.set_xlim(ts.min(), ts.max())
    _frame(b, minor=True)
    b.set_title("delivery into the host separator", fontsize=9, color=S.INK, pad=4)
    _legend(b, size=7.5)

    if _TITLES:
        fig.suptitle(f"In-pipe hydrate/water distribution and host delivery — "
                     f"{_scenario_label(sv)}",
                     color=S.TITLE, fontweight="bold", fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    p = os.path.join(outdir, "17_hydrate_distribution.png")
    return _save(fig, p)


# =============================================================================
#  18 — late-time / shut-in P, T, T_eq, water holdup + deposit growth (scheme: G)
# =============================================================================
def fig_shutin_profile(sv, outdir):
    r, c = sv.results, sv.case
    ts = np.asarray(r["snap_t"], float)
    if ts.size == 0:
        return None
    x = sv.x / 1000.0
    D = float(c.pipeline.diameter_m)
    WC = float(c.fluids.water_cut)

    k = -1
    P = np.asarray(r["snap_P"], float)[k]
    T = np.asarray(r["snap_T"], float)[k]
    Tsub = np.asarray(r["snap_Tsub"], float)[k]
    Teq = T + Tsub
    alpha = np.asarray(r["snap_holdup"], float)[k]
    phi_w = 100.0 * alpha * WC

    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.2))

    a = ax[0]
    lp, = a.plot(x, P, color=S.BLUE, lw=1.6, label="pressure P")
    lt, = a.plot(x, T, color=S.RED, lw=1.6, ls="--", label="temperature T")
    le, = a.plot(x, Teq, color=S.MAGENTA, lw=1.3, ls=":",
                 label="hydrate equilibrium T$_{eq}$")
    a.set_xlabel("distance from wellhead  [km]", fontsize=9)
    a.set_ylabel("pressure P [bar]  or  temperature T [°C]", fontsize=9)
    a.set_xlim(x.min(), x.max())
    _frame(a, minor=True)
    a2 = a.twinx()
    lw, = a2.plot(x, phi_w, color=S.TEAL, lw=1.3, label="water holdup φ$_w$")
    a2.set_ylabel("water volume fraction, φ$_w$  [vol %]", fontsize=9, color=S.TEAL)
    a2.tick_params(axis="y", labelsize=8, colors=S.TEAL)
    a2.set_ylim(0, max(5.0, float(np.nanmax(phi_w)) * 1.25))
    for sp in a2.spines.values():
        sp.set_color(S.INK)
    ev = _event_time_h(sv)
    when = (f"{ts[k] - ev:.0f} h after shut-in" if ev is not None
            else f"after {ts[k]:.0f} h of production")
    a.set_title(f"pipeline profile {when}", fontsize=9, color=S.INK, pad=4)
    _legend(a, handles=[lp, lt, le, lw], size=7.5, anchor=(1.17, 1.0))

    #  (b) deposit volume fraction at successive times
    b = ax[1]
    idx = _pick(ts, 4, 0.25, 1.0)
    dep = np.asarray(r["snap_delta"], float)
    #  four translucent fills stacked on top of one another muddy into a single
    #  brown mass and hide the ordering. Draw the times as distinct curves and
    #  shade only the LAST one, so the final extent still reads as an area.
    cols = [S.BLUE, S.TEAL, S.ORANGE, S.RED]
    prof = []
    for i, (col, kk) in enumerate(zip(cols, idx)):
        ratio = np.clip(2.0 * dep[kk] / D, 0.0, 1.0)
        ph = 100.0 * (1.0 - (1.0 - ratio) ** 2)
        prof.append(ph)
        if i == len(idx) - 1:
            b.fill_between(x, 0, ph, color=col, alpha=0.16, lw=0)
        _plot_series(b, x, ph, i, len(idx), col, _tlab(ts[kk]))
    if np.nanmax(np.asarray(prof, float)) <= 1e-9:
        b.set_ylim(0, 1.0)
        b.text(0.5, 0.5, "no wall deposit forms anywhere on the line",
               transform=b.transAxes, ha="center", va="center", fontsize=9.5,
               fontweight="bold", color=S.INK,
               bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=S.GRIDC, lw=0.9))
    else:
        _coincidence_note(b, prof, "deposit profiles")
    b.set_xlabel("distance from wellhead  [km]", fontsize=9)
    b.set_ylabel("hydrate deposit fraction, φ$_h$  [vol %]", fontsize=9)
    b.set_xlim(x.min(), x.max())
    b.set_ylim(bottom=0)
    _frame(b, minor=True)
    b.set_title("wall-deposit growth along the line", fontsize=9, color=S.INK, pad=4)
    _legend(b, size=7.5, title="elapsed")

    if _TITLES:
        fig.suptitle(f"Late-time pipeline state and deposit growth — {_scenario_label(sv)}",
                     color=S.TITLE, fontweight="bold", fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    p = os.path.join(outdir, "18_shutin_profile_deposit.png")
    return _save(fig, p)


# =============================================================================
#  19 — TRUE space-time fields, 2x2 filled contours          (scheme: H / I)
# =============================================================================
def fig_spacetime_fields(sv, outdir):
    """The TRUE space-time solution: every transported field as a filled-contour
    map over (distance, time), in the layout the two-fluid literature uses."""
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    if ts.size < 3:
        return None
    s_km = sv.x / 1000.0
    D = float(sv.case.pipeline.diameter_m)

    dep = np.asarray(r["snap_delta"], float)
    ratio = np.clip(2.0 * dep / D, 0.0, 1.0)
    dep_pct = 100.0 * (1.0 - (1.0 - ratio) ** 2)

    def _get(key):
        A = np.asarray(r.get(key, np.empty(0)), float)
        return A if (A.ndim == 2 and A.size) else None

    ug = _get("snap_vg")
    ul = _get("snap_vl")
    if ug is None or ul is None:                    # older result set
        ug = ul = np.asarray(r["snap_j"], float)
        ug_lab, ul_lab = "mixture velocity j [m/s]", "mixture velocity j [m/s]"
    else:
        ug_lab, ul_lab = "u$_g$ [m/s]", "u$_l$ [m/s]"

    panels = [
        ("hold-up fraction [-]", np.asarray(r["snap_holdup"], float), "shct_seq"),
        ("pressure [bar]", np.asarray(r["snap_P"], float), "shct_seq"),
        (ug_lab, ug, "shct_seq"),
        (ul_lab, ul, "shct_seq"),
        ("subcooling ΔT$_{sub}$ [°C]", np.asarray(r["snap_Tsub"], float), "shct_temp"),
        ("wall deposit φ$_h$ [vol %]", dep_pct, "shct_heat"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(10.6, 11.0))
    for a, (ttl, Fld, cm) in zip(axes.ravel(), panels):
        finite = Fld[np.isfinite(Fld)]
        if finite.size == 0:
            a.axis("off")
            continue
        lo = float(np.nanpercentile(finite, 0.5))
        hi = float(np.nanpercentile(finite, 99.5))
        #  a field that never varies (no deposit forms at all once the line is
        #  insulated and inhibited, say) has no contours to draw: state that on the
        #  panel rather than printing a colourbar of six identical ticks.
        flat = (hi - lo) <= 1e-9 * max(1.0, abs(lo))
        if flat:
            hi = lo + 1.0
        lv = np.linspace(lo, hi, 26)
        Z, _sf, _tf = S.smooth_field(np.clip(Fld, lo, hi), s_km, ts)
        cf = a.contourf(_sf, _tf, Z, levels=lv, cmap=cm,
                        extend=("neither" if flat else "both"), antialiased=True)
        for coll in cf.collections:
            coll.set_edgecolor("face")           # smooth, seam-free filled bands
        if not flat:
            a.contour(_sf, _tf, Z, levels=lv[::5], colors=CONTOUR_LINE,
                      linewidths=0.45)
        a.set_xlabel("distance from wellhead  [km]", fontsize=9)
        a.set_ylabel("t  [h]", fontsize=9)
        a.set_xlim(s_km.min(), s_km.max())
        a.set_ylim(ts.min(), ts.max())
        a.set_title(ttl, fontsize=10, fontweight="bold", color=S.TITLE, pad=5)
        _frame(a, grid=False)
        ticks = [lo] if flat else list(np.linspace(lo, hi, 6))
        cb = fig.colorbar(cf, ax=a, pad=0.025, fraction=0.055, ticks=ticks)
        cb.ax.tick_params(labelsize=7.5)
        span = hi - lo
        fmt = "%.0f" if (flat or span >= 20) else ("%.1f" if span >= 2 else "%.2f")
        cb.ax.set_yticklabels([fmt % v for v in ticks])
        cb.outline.set_edgecolor(S.INK)
        cb.outline.set_linewidth(0.8)
        if flat:
            a.text(0.5, 0.5, f"uniform at {lo:g} for the whole run",
                   transform=a.transAxes, ha="center", va="center",
                   fontsize=10, fontweight="bold", color=S.INK,
                   bbox=dict(boxstyle="round,pad=0.4", fc="white",
                             ec=S.GRIDC, lw=0.9))

    if _TITLES:
        fig.suptitle(f"Space-time solution of the 32 km tie-back — {_scenario_label(sv)} "
                     f"(N = {sv.x.size} cells, {ts.size} snapshots)",
                     color=S.TITLE, fontweight="bold", fontsize=10.5, y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.972))
    p = os.path.join(outdir, "19_spacetime_fields.png")
    return _save(fig, p)


# =============================================================================
#  20 — holdup along the line after successive durations     (scheme: J)
# =============================================================================
def fig_holdup_durations(sv, outdir):
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    H = np.asarray(r["snap_holdup"], float)
    if H.size == 0:
        return None
    x = sv.x / 1000.0
    ev = _event_time_h(sv)

    if ev is not None:
        keep = np.where(ts >= ev)[0]
        if keep.size < 3:
            keep = np.arange(ts.size)
        sel = keep[_pick(ts[keep], 5, 0.08, 1.0)]
        legend_title = "shut-in duration"
        durations = ts[sel] - ev
    else:
        sel = _pick(ts, 5, 0.15, 1.0)
        legend_title = "production time"
        durations = ts[sel]

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    for i, (col, k, d) in enumerate(zip(TIME_COLORS, sel, durations)):
        _plot_series(ax, x, H[k], i, len(sel), col, _tlab(max(d, 0.0)))
    _coincidence_note(ax, [H[k] for k in sel], "holdup profiles")
    ax.set_xlabel("distance from wellhead  [km]", fontsize=10)
    ax.set_ylabel("liquid holdup  [-]", fontsize=10)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1.02)
    _frame(ax, minor=True)
    _legend(ax, size=8.5, title=legend_title)
    _title(ax, f"Distribution of liquid holdup along the pipeline after different "
               f"{'shut-in' if ev is not None else 'production'} durations — "
               f"{_scenario_label(sv)}", size=9.5)

    fig.tight_layout()
    p = os.path.join(outdir, "20_holdup_durations.png")
    return _save(fig, p)


# =============================================================================
#  21 — riser depth-time waterfall with slug boundaries      (scheme: K)
# =============================================================================
def fig_riser_depth_time(sv, outdir):
    """Depth-time waterfall over the steel-catenary riser, with the slug-boundary
    trajectories and the slug length annotated — the published waterfall scheme."""
    x = np.asarray(sv.x, float)
    z = np.asarray(sv.z, float)

    climb = (np.gradient(z, x) > 0.02) & (x > 0.8 * x[-1])
    i0 = int(np.argmax(climb)) if climb.any() else int(0.94 * (x.size - 1))
    i0 = int(np.clip(i0, 1, x.size - 4))

    #  Choose the instant by what the RISER is doing, not the flowline — see
    #  _riser_slug_snapshot. Fall back to the flowline choice so a scenario whose
    #  riser never slugs still reports through the guard below rather than here.
    k = _riser_slug_snapshot(sv, i0)
    if k is None:
        k = _slug_snapshot(sv)
    if k is None:
        print("    [space-time] 21_riser_depth_time skipped: the line never reaches "
              "an intermittent, flowing state in this scenario", flush=True)
        return None
    F = slug_unit_fields(sv, k_snap=k)
    t_snap = float(np.asarray(sv.results["snap_t"], float)[k])

    Lu_c = float(np.mean(F["Lu"][i0:]))
    Vt_c = float(np.mean(F["Vt"][i0:]))
    Ls_c = float(np.mean(F["beta"][i0:])) * Lu_c
    riser_len = abs(x[-1] - x[i0])
    #  a slug unit is a metre-scale structure: tens of pipe diameters, not
    #  hundreds of metres. Reject anything that cannot be a slug train before it
    #  is drawn, so an aliased moire pattern can never reach the page.
    D_pipe = float(sv.case.pipeline.diameter_m)
    if (not np.isfinite(Lu_c) or Lu_c <= 0.0 or Lu_c > 0.10 * riser_len
            or Lu_c > 400.0 * D_pipe or Vt_c < 0.2):
        print(f"    [space-time] 21_riser_depth_time skipped: no travelling slug "
              f"train in the riser at this state (L_u = {Lu_c:.0f} m over a "
              f"{riser_len:.0f} m riser, V_t = {Vt_c:.2f} m/s)", flush=True)
        plt.close("all")
        return None

    #  a depth window ~6 slug units deep and ~9 slug periods long: enough units to
    #  read the train, coarse enough that every unit is resolved on the page
    depth_top = float(-z[-1])
    depth_span = float(np.clip(4.0 * Lu_c, 15.0, abs(z[i0] - z[-1])))
    depth_lo, depth_hi = depth_top, depth_top + depth_span
    zq = -np.linspace(depth_hi, depth_lo, 460)          # deep -> shallow
    xq = np.interp(zq, z[i0:], x[i0:])
    depth = -zq

    t_transit = depth_span / max(Vt_c, 1e-3)
    period = Lu_c / max(Vt_c, 1e-3)
    window_s = 1.4 * t_transit
    tq = np.linspace(0.0, window_s, 620)

    fld, meta = reconstruct_slug_field(sv, xq, tq, k_snap=k)

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    pcm = ax.pcolormesh(tq, depth, fld.T, cmap="shct_seq", shading="gouraud",
                        vmin=0.0, vmax=1.0)
    ax.set_xlim(tq.min(), tq.max())
    ax.set_ylim(depth_hi, depth_lo)                     # depth increases downward
    ax.set_xlabel("time  [s]", fontsize=10)
    ax.set_ylabel("depth from host  [m]", fontsize=10)
    _frame(ax, grid=False)
    cb = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("liquid holdup  α$_l$  [-]", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(S.INK)

    #  slug-boundary trajectories: each front leaves depth_hi at t_a and reaches
    #  depth_lo one transit later; clip every line to the plotted window
    #  annotate only every third front, so the trajectories read as guides
    t_a = -t_transit
    n = 0
    while t_a < window_s and n < 40:
        if n % 3 == 0:
            ta, tb = max(t_a, 0.0), min(t_a + t_transit, window_s)
            if tb > ta:
                ax.plot([ta, tb],
                        [depth_hi - Vt_c * (ta - t_a), depth_hi - Vt_c * (tb - t_a)],
                        color=S.MAGENTA, lw=1.2, ls="--", alpha=0.95)
        t_a += period
        n += 1
    ax.plot([], [], color=S.MAGENTA, lw=1.0, ls="--",
            label=f"slug-boundary trajectories (V$_t$ = {Vt_c:.2f} m/s)")

    #  annotate the slug body length the way the published waterfall does
    #  the front-to-front separation IS the slug unit length; mark it, and give
    #  the body length inside it (the published waterfall marks the same span)
    #  the y-axis is DEPTH while L_u is a length along the (steeply inclined)
    #  riser, so project the unit onto the depth axis before marking it
    sin_theta = depth_span / max(abs(xq[0] - xq[-1]), 1e-6)
    Lu_d, Ls_d = Lu_c * sin_theta, Ls_c * sin_theta
    t_mid = 0.62 * window_s
    d_mid = depth_lo + 0.22 * depth_span
    #  the measurement ARROW must stay in data coordinates -- a scale bar drawn
    #  anywhere else measures nothing -- but its VALUE is stated in the margin.
    ax.annotate("", xy=(t_mid, d_mid + Lu_d), xytext=(t_mid, d_mid),
                arrowprops=dict(arrowstyle="<->", color=S.RED, lw=1.8))
    #  the right margin carries the colourbar, so the value goes in the LEFT one
    _margin_note(ax, d_mid + 0.5 * Lu_d,
                 f"L$_u$ = {Lu_c:.1f} m along riser\n"
                 f"({Lu_d:.1f} m of depth, body {Ls_c:.1f} m)",
                 side="left", color=S.RED, pad=0.055, size=9.0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=8,
              framealpha=1.0, facecolor="white", edgecolor=S.INK,
              borderaxespad=0.0)
    _title(ax, f"Riser depth–time waterfall — slug boundaries during upward motion "
               f"({_scenario_label(sv)}, t = {t_snap:.1f} h)", size=9.5)
    fig.text(0.5, 0.005,
             f"Mass-consistent sub-grid reconstruction over the steel-catenary riser "
             f"(L$_u$ = {Lu_c:.1f} m, slug body {Ls_c:.1f} m, period {period:.1f} s); "
             f"celerity, period and length are solver outputs.",
             ha="center", fontsize=6.8, color=S.INK, style="italic")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    p = os.path.join(outdir, "21_riser_depth_time.png")
    return _save(fig, p)


# =============================================================================
#  22 — pipeline cloud maps at successive times              (scheme: L)
# =============================================================================
def fig_cloud_maps(sv, outdir, n_times=3, ny=90):
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    if ts.size < 3:
        return None
    x = sv.x / 1000.0
    H = np.asarray(r["snap_holdup"], float)
    Tf = np.asarray(r["snap_T"], float)

    idx = _pick(ts, n_times, 0.35, 1.0)[:n_times]
    yy = np.linspace(0.0, 1.0, ny)

    #  two-colour phase map: liquid below the interface, gas above (neither dark)
    phase_cmap = ListedColormap([S.BLUE, S.ORANGE])
    phase_norm = BoundaryNorm([0, 0.5, 1], phase_cmap.N)

    tmin = float(np.nanpercentile(Tf[idx], 1))
    tmax = float(np.nanpercentile(Tf[idx], 99))
    if tmax <= tmin:
        tmax = tmin + 1.0

    fig = plt.figure(figsize=(10.6, 1.70 * len(idx) + 2.6))
    gs = fig.add_gridspec(2 * len(idx), 1,
                          height_ratios=[1.0, 0.62] * len(idx),
                          hspace=0.46, top=0.86, bottom=0.20,
                          left=0.055, right=0.985)

    for row, k in enumerate(idx):
        alpha = np.clip(H[k], 0.0, 1.0)
        #  (upper strip) oil-water / gas phase distribution inside the bore
        ap = fig.add_subplot(gs[2 * row])
        liquid = (yy[:, None] <= alpha[None, :]).astype(float)
        ap.pcolormesh(x, yy, 1.0 - liquid, cmap=phase_cmap, norm=phase_norm,
                      shading="auto")
        ap.plot(x, alpha, color="white", lw=1.0)
        ap.set_yticks([])
        ap.set_xticks([])
        ap.set_xlim(x.min(), x.max())
        ap.set_ylim(0, 1)
        for sp in ap.spines.values():
            sp.set_color(S.INK)
            sp.set_linewidth(0.9)
        ap.text(-0.018, 0.5, f"({chr(97 + row)})", transform=ap.transAxes,
                ha="right", va="center", fontsize=10, fontweight="bold",
                color=S.INK)
        ap.text(0.5, 1.06, f"t = {ts[k]:.1f} h", transform=ap.transAxes,
                ha="center", va="bottom", fontsize=8.5, color=S.RED,
                fontweight="bold")

        #  (lower strip) bulk-temperature contour along the same reach
        at = fig.add_subplot(gs[2 * row + 1])
        Tstrip = np.repeat(Tf[k][None, :], ny, axis=0)
        pcm = at.pcolormesh(x, yy, Tstrip, cmap="shct_temp", shading="gouraud",
                            vmin=tmin, vmax=tmax)
        at.set_yticks([])
        at.set_xlim(x.min(), x.max())
        at.set_ylim(0, 1)
        for sp in at.spines.values():
            sp.set_color(S.INK)
            sp.set_linewidth(0.9)
        if row == len(idx) - 1:
            at.set_xlabel("distance from wellhead  [km]", fontsize=9)
            at.tick_params(labelsize=8)
        else:
            at.set_xticks([])

    cax = fig.add_axes([0.16, 0.075, 0.70, 0.022])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cb.set_label("bulk temperature  [°C]", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_edgecolor(S.INK)

    handles = [plt.Line2D([], [], color=S.ORANGE, lw=6, label="gas phase"),
               plt.Line2D([], [], color=S.BLUE, lw=6, label="liquid phase (oil + water)")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=2, fontsize=8, framealpha=1.0, facecolor="white",
               edgecolor=S.INK)
    if _TITLES:
        fig.suptitle(f"Pipeline cloud maps at successive times — phase distribution "
                     f"(upper) and temperature (lower), {_scenario_label(sv)}",
                     color=S.TITLE, fontweight="bold", fontsize=10, y=0.988)
    p = os.path.join(outdir, "22_cloud_maps.png")
    return _save(fig, p)


# =============================================================================
#  23 — DTS-style thermal waterfall with the pressure trace and the stages
# =============================================================================
def _stages(sv):
    """The operating stages of the run, as (t_start, t_end, label) in hours."""
    ts = np.asarray(sv.results["snap_t"], float)
    if ts.size == 0:
        return []
    t0, t1 = float(ts[0]), float(ts[-1])
    ev = _event_time_h(sv)
    if ev is None or not (t0 < ev < t1):
        lab = ("Engineered fix — insulation + MEG"
               if float(getattr(sv.case.operating, "MEG_wt_inlet", 0.0) or 0.0) > 0
               else "As-operated production")
        return [(t0, t1, lab)]
    return [(t0, ev, "Production"), (ev, t1, "Shut-in cooldown")]


def fig_dts_waterfall(sv, outdir):
    """Temperature over (distance, time) rendered as a distributed-temperature
    waterfall: the operating stages are marked along the top, and the monitored
    pressure is overlaid so the thermal response can be read against it."""
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    if ts.size < 3:
        return None
    x = sv.x / 1000.0
    T = np.asarray(r["snap_T"], float)
    P = np.asarray(r["snap_P"], float)

    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    lo = float(np.nanpercentile(T, 0.5))
    hi = float(np.nanpercentile(T, 99.5))
    if hi <= lo:
        hi = lo + 1.0
    #  time on x and distance on y, increasing downward, is the waterfall
    #  convention: a front then reads as a sloping edge whose slope is its speed
    _Ts, _tf, _xf = S.smooth_field(np.clip(T, lo, hi).T, ts, x)
    pcm = ax.pcolormesh(_tf, _xf, _Ts, cmap="shct_dts",
                        shading="gouraud", vmin=lo, vmax=hi)
    ax.set_ylim(x.max(), x.min())
    ax.set_xlim(ts.min(), ts.max())
    ax.set_xlabel("time  [h]", fontsize=10)
    ax.set_ylabel("distance from wellhead  [km]", fontsize=10)
    _frame(ax, grid=False)
    #  the pressure trace needs the right-hand spine, so the colourbar is pushed
    #  clear of it rather than sharing the same margin
    cb = fig.colorbar(pcm, ax=ax, pad=0.115, fraction=0.042)
    cb.set_label("temperature  [°C]", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(S.INK)

    #  the monitored pressure, overlaid on its own axis (the DTS convention)
    _num = getattr(sv.case, "numerics", None)
    mon = int(float(getattr(_num, "monitor_frac", 0.8) if _num is not None else 0.8)
              * (x.size - 1))
    pax = ax.twinx()
    pax.plot(ts, P[:, mon], color=S.INK, lw=1.6, solid_capstyle="round")
    pax.set_ylabel(f"pressure at the monitor ({x[mon]:.0f} km)  [bar]",
                   fontsize=9, color=S.INK)
    pax.tick_params(axis="y", labelsize=8, colors=S.INK)
    pax.set_ylim(0, float(np.nanmax(P)) * 1.35)
    for sp in pax.spines.values():
        sp.set_color(S.INK)
    pax.plot([], [], color=S.INK, lw=1.6, label="monitored pressure")
    pax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.20), fontsize=8,
               framealpha=1.0, facecolor="white", edgecolor=S.INK,
               borderaxespad=0.0)

    #  stage names ABOVE the axes; the hydrate-onset distance labelled in the
    #  LEFT margin against its own coordinate. Only the dotted feature line is
    #  drawn on the field itself.
    _title_pad = _stage_header(ax, _stages(sv))
    Tsub = np.asarray(r.get("snap_Tsub", np.empty(0)), float)
    if Tsub.ndim == 2 and Tsub.size:
        sub = Tsub[-1] > 0.0
        frac = float(sub.mean())
        if frac >= 0.90:
            #  A cooled-down line is subcooled END TO END, so "the first subcooled
            #  cell" is cell 0 and quoting it as an onset distance is meaningless.
            #  State the real result instead.
            _margin_note(ax, float(np.median(x)),
                         f"whole line inside the\nhydrate region ({frac*100:.0f} %)",
                         side="left", pad=0.06, leader=False)
        elif sub.any():
            x_on = float(x[int(np.argmax(sub))])
            ax.axhline(x_on, color="white", lw=1.4, ls=":")
            _margin_note(ax, x_on, f"hydrate onset\n≈ {x_on:.1f} km", side="left",
                         pad=0.06)

    if _TITLES:
        ax.set_title(f"Distributed-temperature waterfall T(x, t) — "
                     f"{_scenario_label(sv)}", color=S.TITLE, fontweight="bold",
                     fontsize=10, pad=6 + _title_pad)
    fig.tight_layout()
    p = os.path.join(outdir, "23_dts_thermal_waterfall.png")
    return _save(fig, p)


# =============================================================================
#  24 — temperature-gradient waterfall (the front detector)
# =============================================================================
def fig_gradient_waterfall(sv, outdir):
    """dT/dx over (distance, time). A travelling thermal front is a narrow band
    of steep gradient, so it stands out here even where the temperature map
    itself looks smooth — the same reason a distributed-sensing record is
    differentiated before it is read."""
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    if ts.size < 3:
        return None
    x_km = sv.x / 1000.0
    T = np.asarray(r["snap_T"], float)
    G = np.gradient(T, sv.x, axis=1) * 1000.0            # °C per km

    lim = float(np.nanpercentile(np.abs(G), 99.0))
    if not np.isfinite(lim) or lim <= 0:
        lim = max(float(np.nanmax(np.abs(G))), 1e-6)

    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    _Gs, _xf, _tf = S.smooth_field(np.clip(G, -lim, lim), x_km, ts)
    pcm = ax.pcolormesh(_xf, _tf, _Gs, cmap="shct_grad",
                        shading="gouraud", vmin=-lim, vmax=lim)
    ax.set_xlabel("distance from wellhead  [km]", fontsize=10)
    ax.set_ylabel("time  [h]", fontsize=10)
    ax.set_xlim(x_km.min(), x_km.max())
    ax.set_ylim(ts.min(), ts.max())
    _frame(ax, grid=False)
    cb = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("temperature gradient  ∂T/∂x  [°C km$^{-1}$]", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(S.INK)

    #  the steepest-gradient location at each time IS the front; track it
    front = x_km[np.nanargmin(G, axis=1)]
    ax.plot(front, ts, color=S.INK, lw=1.5, ls="--")
    ax.plot([], [], color=S.INK, lw=1.5, ls="--",
            label="steepest cooling front")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.26), fontsize=8,
              framealpha=1.0, facecolor="white", edgecolor=S.INK,
              borderaxespad=0.0)

    _title(ax, f"Temperature-gradient waterfall ∂T/∂x (x, t) — "
               f"{_scenario_label(sv)}", size=10)
    fig.text(0.5, 0.005,
             "A travelling thermal front is a narrow band of steep gradient; the "
             "dashed line tracks the steepest cooling at each instant.",
             ha="center", fontsize=6.8, color=S.INK, style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    p = os.path.join(outdir, "24_temperature_gradient.png")
    return _save(fig, p)


# =============================================================================
#  25 — DAS-style flow-noise waterfall
# =============================================================================
def fig_das_waterfall(sv, outdir):
    """|dα_l/dt| over (distance, time): where the holdup is changing fastest is
    where the flow is most unsteady, which is what a distributed-acoustic record
    of a producing line actually shows."""
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    H = np.asarray(r["snap_holdup"], float)
    if ts.size < 4 or H.size == 0:
        return None
    x_km = sv.x / 1000.0
    dt = np.gradient(ts)
    A = np.abs(np.gradient(H, axis=0) / np.maximum(dt[:, None], 1e-9))   # 1/h

    hi = float(np.nanpercentile(A, 99.0))
    if not np.isfinite(hi) or hi <= 0:
        hi = max(float(np.nanmax(A)), 1e-9)

    fig, ax = plt.subplots(figsize=(10.4, 5.2))
    _As, _tf, _xf = S.smooth_field(np.clip(A, 0.0, hi).T, ts, x_km)
    pcm = ax.pcolormesh(_tf, _xf, _As, cmap="shct_dts",
                        shading="gouraud", vmin=0.0, vmax=hi)
    ax.set_ylim(x_km.max(), x_km.min())
    ax.set_xlim(ts.min(), ts.max())
    ax.set_xlabel("time  [h]", fontsize=10)
    ax.set_ylabel("distance from wellhead  [km]", fontsize=10)
    _frame(ax, grid=False)
    cb = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("holdup fluctuation rate  |∂α$_l$/∂t|  [h$^{-1}$]", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(S.INK)

    #  mark the intermittent reach and the riser, the two noise sources.
    #  Take the regime from a snapshot where the line is actually FLOWING: after a
    #  shut-in nothing is intermittent, so the final state would collapse the
    #  bracket onto a single cell and label a reach that does not exist.
    reg = np.asarray(r.get("snap_regime", np.empty(0)), float)
    if reg.ndim == 2 and reg.size:
        _k = _slug_snapshot(sv)
        inter = np.isin(np.round(reg[_k if _k is not None else -1]), [2, 5])
        #  only bracket a reach that is actually a reach
        if inter.mean() < 0.05:
            inter = np.zeros_like(inter)
        if inter.any():
            i0, i1 = int(np.argmax(inter)), int(len(inter) - 1 - np.argmax(inter[::-1]))
            #  bracket the reach on the axis spine, label it in the margin
            trans = ax.get_yaxis_transform()
            ax.annotate("", xy=(-0.012, x_km[i0]), xycoords=trans,
                        xytext=(-0.012, x_km[i1]), textcoords=trans,
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="<->", color=S.INK, lw=2.0))
            _margin_note(ax, 0.5 * (x_km[i0] + x_km[i1]),
                         "intermittent\n(slug / churn)", side="left", pad=0.055,
                         leader=False)
    z = np.asarray(sv.z, float)
    climb = (np.gradient(z, np.asarray(sv.x, float)) > 0.02) & (x_km > 0.8 * x_km.max())
    if climb.any():
        _xb = float(x_km[int(np.argmax(climb))])
        ax.axhline(_xb, color="white", lw=1.4, ls=":")
        #  the RIGHT margin is occupied by the colourbar, so this goes left
        _margin_note(ax, _xb, "riser base", side="left", pad=0.055)

    _title(ax, f"Flow-noise waterfall |∂α$_l$/∂t| (x, t) — {_scenario_label(sv)}",
           size=10)
    fig.tight_layout()
    p = os.path.join(outdir, "25_das_flow_noise.png")
    return _save(fig, p)


# =============================================================================
#  26 — the along-pipeline parameter diagram at successive times
# =============================================================================
def fig_parameter_panels(sv, outdir, n_times=5):
    """Pressure, temperature, liquid holdup and mixture velocity along the route,
    each at the same successive times — the standard way a transient pipeline
    study reports the state of the line as an event develops."""
    r = sv.results
    ts = np.asarray(r["snap_t"], float)
    if ts.size < 3:
        return None
    x = sv.x / 1000.0
    ev = _event_time_h(sv)

    if ev is not None:
        keep = np.where(ts >= ev)[0]
        idx = keep[_pick(ts[keep], n_times, 0.02, 1.0)] if keep.size >= 3 \
            else _pick(ts, n_times, 0.1, 1.0)
        lab = lambda k: _tlab(max(ts[k] - ev, 0.0))
        legend_title = "elapsed since shut-in"
    else:
        idx = _pick(ts, n_times, 0.12, 1.0)
        lab = lambda k: _tlab(ts[k])
        legend_title = "production time"

    panels = [
        ("(a) pressure", "pressure  P  [bar]", np.asarray(r["snap_P"], float), None),
        ("(b) temperature", "temperature  T  [°C]", np.asarray(r["snap_T"], float),
         np.asarray(r["snap_T"], float) + np.asarray(r["snap_Tsub"], float)),
        ("(c) liquid holdup", "holdup  α$_l$  [-]",
         np.asarray(r["snap_holdup"], float), None),
        ("(d) mixture velocity", "velocity  j  [m s$^{-1}$]",
         np.asarray(r["snap_j"], float), None),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.4))
    for a, (ttl, ylab, F, ref) in zip(axes.ravel(), panels):
        for i, k in enumerate(idx):
            _plot_series(a, x, F[k], i, len(idx), TIME_COLORS[i % len(TIME_COLORS)],
                         lab(k))
        if ref is not None:
            #  the hydrate-equilibrium temperature, so the crossing is visible
            a.plot(x, ref[idx[-1]], color=S.MAGENTA, lw=1.4, ls=":",
                   label="hydrate T$_{eq}$")
        a.set_xlabel("distance from wellhead  [km]", fontsize=9)
        a.set_ylabel(ylab, fontsize=9)
        a.set_xlim(x.min(), x.max())
        if "holdup" in ylab:
            a.set_ylim(0, 1.02)
        _frame(a, minor=True)
        a.set_title(ttl, fontsize=9.5, color=S.TITLE, fontweight="bold", pad=5)
        _coincidence_note(a, [F[k] for k in idx])
    _legend(axes[0, 1], size=7.5, title=legend_title)

    if _TITLES:
        fig.suptitle(f"Parameters along the pipeline at successive times — "
                     f"{_scenario_label(sv)}",
                     color=S.TITLE, fontweight="bold", fontsize=10.5, y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    p = os.path.join(outdir, "26_parameter_panels.png")
    return _save(fig, p)


# =============================================================================
#  27 — well-posedness / Kelvin-Helmholtz map
# =============================================================================
def _ikh_slip_limit(alpha_l, D, rho_l, rho_g, theta):
    """The inviscid Kelvin-Helmholtz slip limit for stratified flow in a pipe.

    The 1-D two-fluid model stays hyperbolic (and the initial-value problem well
    posed) while the phase slip satisfies

        (u_g - u_l)^2  <  (rho_l - rho_g) g cos(theta) (A_g/rho_g + A_l/rho_l)
                          / (dA_l/dh)

    with h the liquid level and A_l(h) the liquid area. Beyond it the
    characteristics become complex: the growth rate then rises without bound as
    the grid is refined, so the "instability" is a property of the equations, not
    of the flow. Returns the limiting |u_g - u_l| in m/s.
    """
    alpha_l = np.clip(np.asarray(alpha_l, float), 1e-4, 1.0 - 1e-4)
    A = math.pi * D ** 2 / 4.0
    A_l, A_g = alpha_l * A, (1.0 - alpha_l) * A
    #  liquid level from the area fraction, and the interface width dA_l/dh = w
    gamma = 2.0 * np.arccos(np.clip(1.0 - 2.0 * alpha_l, -1.0, 1.0))   # wetted angle
    w = np.maximum(D * np.sin(gamma / 2.0), 1e-6)                      # chord width
    num = (rho_l - rho_g) * G * np.cos(np.asarray(theta, float)) * \
        (A_g / max(rho_g, 1e-6) + A_l / max(rho_l, 1e-6))
    return np.sqrt(np.maximum(num / w, 0.0))


def fig_wellposedness(sv, outdir):
    """Where the case sits relative to the two-fluid well-posedness boundary.

    The left panel maps the boundary over the superficial-velocity plane at the
    line's mean inclination; the right panel reports the margin cell by cell
    along the route, so any reach that has crossed it is identified.
    """
    import solver as _solver

    r, c = sv.results, sv.case
    x = sv.x / 1000.0
    D = float(c.pipeline.diameter_m)
    rho_l = float(c.fluids.rho_oil) * (1.0 - float(c.fluids.water_cut)) + \
        float(c.fluids.rho_water) * float(c.fluids.water_cut)

    #  prefer the snapshot history (it is what a restored state carries); fall
    #  back to the ensemble-median final state of a live run
    def _last(key_snap, key_final):
        A = np.asarray(r.get(key_snap, np.empty(0)), float)
        if A.ndim == 2 and A.size:
            return A[-1]
        return np.nanmedian(np.asarray(r[key_final], float), 1)

    alpha = np.clip(_last("snap_holdup", "alpha_l"), 1e-3, 0.999)
    P = _last("snap_P", "p")
    T = _last("snap_T", "T")
    try:
        rho_g = _solver.gas_density(P, T, c.fluids)
    except Exception:
        rho_g = (P * 1e5) / (8.3145 / 0.019 * (T + 273.15))
    rho_g = np.maximum(np.asarray(rho_g, float), 1e-3)
    th = _theta(sv)

    #  the actual slip along the line, from the solver's own phase velocities
    vl = np.asarray(r.get("snap_vl", np.empty(0)), float)
    vg = np.asarray(r.get("snap_vg", np.empty(0)), float)
    if vl.ndim == 2 and vl.size:
        slip = np.abs(vg[-1] - vl[-1])
        v_sg = (1.0 - alpha) * vg[-1]
        v_sl = alpha * vl[-1]
    else:
        j = _last("snap_j", "j")
        slip = 0.25 * j
        v_sg, v_sl = (1.0 - alpha) * j, alpha * j
    limit = _ikh_slip_limit(alpha, D, rho_l, float(np.mean(rho_g)), th)
    margin = slip / np.maximum(limit, 1e-9)

    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.6),
                           gridspec_kw=dict(width_ratios=[1.0, 1.15]))

    #  ---- (a) the boundary over the superficial-velocity plane ----------------
    a = ax[0]
    vsg = np.logspace(-2, 1.2, 220)
    vsl = np.logspace(-2, 1.0, 200)
    Vg, Vl = np.meshgrid(vsg, vsl)
    #  a homogeneous-slip estimate of the holdup on the map, closed by the same
    #  drift-flux the solver uses, so the map and the run share one closure
    C0, vd = _solver.drift_params(float(np.mean(th)), D)
    Vm = Vg + Vl
    al_map = np.clip(1.0 - Vg / np.maximum(C0 * Vm + vd, 1e-6), 1e-3, 0.999)
    lim_map = _ikh_slip_limit(al_map, D, rho_l, float(np.mean(rho_g)),
                              float(np.mean(th)))
    slip_map = np.abs(Vg / np.maximum(1.0 - al_map, 1e-3)
                      - Vl / np.maximum(al_map, 1e-3))
    ratio = slip_map / np.maximum(lim_map, 1e-9)

    lv = np.linspace(0.0, 2.0, 21)
    cf = a.contourf(vsg, vsl, np.clip(ratio, 0, 2.0), levels=lv, cmap="shct_seq",
                    extend="max")
    for coll in cf.collections:
        coll.set_edgecolor("face")
    a.contour(vsg, vsl, ratio, levels=[1.0], colors=[S.MAGENTA], linewidths=2.0)
    a.plot([], [], color=S.MAGENTA, lw=2.0, label="well-posedness boundary")
    a.scatter(np.maximum(v_sg, 1e-2), np.maximum(v_sl, 1e-2), s=16,
              facecolor="white", edgecolor=S.INK, linewidth=0.8, zorder=5,
              label="case states along the route")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("superficial gas velocity  V$_{sg}$  [m s$^{-1}$]", fontsize=9)
    a.set_ylabel("superficial liquid velocity  V$_{sl}$  [m s$^{-1}$]", fontsize=9)
    a.set_title("(a) two-fluid well-posedness map", fontsize=9.5, color=S.TITLE,
                fontweight="bold", pad=5)
    _frame(a, grid=False)
    cb = fig.colorbar(cf, ax=a, pad=0.02, fraction=0.05)
    cb.set_label("slip / Kelvin-Helmholtz limit", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_edgecolor(S.INK)
    #  Below the axes at -0.22 the legend sits exactly where the x-axis label is,
    #  which clears at print size and collides once the type is scaled up for a
    #  slide — a 100 % overlap of the legend on "superficial gas velocity". The
    #  project rule for every other panel is legend OUTSIDE to the right, where
    #  nothing it can collide with lives, so this one follows it too.
    _legend(a, size=7.0, anchor=(1.012, 1.0))

    #  ---- (b) the margin along the route -------------------------------------
    b = ax[1]
    b.plot(x, margin, color=S.BLUE, lw=1.8, label="slip / KH limit")
    b.axhline(1.0, color=S.MAGENTA, lw=1.6, ls="--", label="well-posedness limit")
    b.fill_between(x, 0, margin, where=margin >= 1.0, color=S.HYDFILL, alpha=0.8,
                   label="ill-posed reach")
    b.set_xlabel("distance from wellhead  [km]", fontsize=9)
    b.set_ylabel("slip / Kelvin-Helmholtz limit  [-]", fontsize=9)
    b.set_xlim(x.min(), x.max())
    b.set_ylim(0, max(1.35, float(np.nanmax(margin)) * 1.15))
    _frame(b, minor=True)
    b.set_title("(b) margin along the route", fontsize=9.5, color=S.TITLE,
                fontweight="bold", pad=5)
    _legend(b, size=7.5)

    frac = float(np.mean(margin >= 1.0)) * 100.0
    verdict = ("the model stays hyperbolic over the whole route, so the predicted "
               "slug activity is a property of the flow, not of the discretisation"
               if frac < 0.5 else
               f"{frac:.0f} % of the route is past the limit — over that reach the "
               f"two-fluid initial-value problem is ill-posed and the growth rate "
               f"is grid-dependent")
    fig.text(0.5, 0.005, verdict, ha="center", fontsize=7.2, style="italic",
             color=S.INK)

    if _TITLES:
        fig.suptitle(f"Well-posedness of the two-fluid description — "
                     f"{_scenario_label(sv)}",
                     color=S.TITLE, fontweight="bold", fontsize=10.5, y=0.997)
    fig.tight_layout(rect=(0, 0.035, 1, 0.955))
    p = os.path.join(outdir, "27_wellposedness_map.png")
    return _save(fig, p)


# =============================================================================
#  re-rendering without re-running: the snapshot state is saved alongside the
#  figures, so any of them can be rebuilt (restyled, rescaled, re-cropped) from
#  the archived run instead of repeating a multi-hour transient.
# =============================================================================
STATE_FILE = "spacetime_state.npz"

_STATE_KEYS = ("snap_t", "snap_holdup", "snap_P", "snap_T", "snap_phi",
               "snap_PhiSH", "snap_delta", "snap_Tsub", "snap_j", "snap_vl",
               "snap_vg", "snap_regime", "snap_fslug")


def save_state(sv, outdir):
    """Write everything the figures in this module read into one .npz."""
    r = sv.results
    data = {"x": np.asarray(sv.x, float), "z": np.asarray(sv.z, float)}
    for k in _STATE_KEYS:
        A = np.asarray(r.get(k, np.empty(0)), float)
        if A.size:
            data[k] = A
    c = sv.case
    data["_case"] = np.array([
        float(c.pipeline.diameter_m), float(c.fluids.water_cut),
        float(c.fluids.rho_oil), float(c.fluids.rho_water),
        float(getattr(c.operating, "MEG_wt_inlet", 0.0) or 0.0),
        float(getattr(c.scenario, "event_time_h", 0.0) or 0.0),
        1.0 if getattr(c.scenario, "kind", "steady") == "shutin" else 0.0,
    ], float)
    path = os.path.join(outdir, STATE_FILE)
    np.savez_compressed(path, **data)
    return path


class _State:
    """A minimal stand-in for a solved TransientSHCT, restored from the .npz.

    It exposes exactly the attributes the figure functions use: x, z, results
    and a case whose pipeline / fluids / operating / scenario groups carry the
    handful of scalars the figures read.  A real solver Case is rebuilt from the
    run's own case_config.json when that sits next to the state file, so the
    compositional helpers (gas density, in particular) see the true fluid.
    """

    class _Grp:
        def __init__(self, **kw):
            self.__dict__.update(kw)

    def __init__(self, npz_path):
        d = np.load(npz_path, allow_pickle=False)
        self.x = d["x"]
        self.z = d["z"]
        self.results = {k: d[k] for k in d.files if k.startswith("snap_")}
        D, wc, ro, rw, meg, ev, is_shutin = d["_case"]
        self.case = _State._Grp(
            pipeline=_State._Grp(diameter_m=float(D), length_m=float(self.x[-1]),
                                 n_cells=int(self.x.size)),
            fluids=_State._Grp(water_cut=float(wc), rho_oil=float(ro),
                               rho_water=float(rw)),
            operating=_State._Grp(MEG_wt_inlet=float(meg)),
            numerics=_State._Grp(monitor_frac=0.8),
            scenario=_State._Grp(kind="shutin" if is_shutin else "steady",
                                 event_time_h=float(ev)))
        cfg = os.path.join(os.path.dirname(npz_path), "case_config.json")
        if os.path.exists(cfg):
            try:
                import json
                import solver as _solver
                raw = json.load(open(cfg))
                real = _solver.Case()
                for grp, vals in raw.items():
                    tgt = getattr(real, grp, None)
                    if tgt is None or not isinstance(vals, dict):
                        continue
                    for kk, vv in vals.items():
                        if hasattr(tgt, kk):
                            setattr(tgt, kk, vv)
                self.case = real
            except Exception:
                pass


def load_state(path):
    """Restore a run's snapshot state; `path` may be the .npz or its folder."""
    if os.path.isdir(path):
        path = os.path.join(path, STATE_FILE)
    return _State(path)


def rerender(outdir, verbose=True):
    """Rebuild the whole figure set for an output folder from its saved state."""
    sv = load_state(outdir)
    return spacetime_outputs(sv, None, outdir, verbose=verbose, save=False)


# =============================================================================
#  driver
# =============================================================================
FIGURES = [
    ("14_holdup_multitime.png", fig_holdup_multitime, "sv"),
    ("15_slug_growth_propagation.png", fig_slug_growth, "sv"),
    ("16_slug_train_waterfall.png", fig_slug_waterfall, "sv"),
    ("17_hydrate_distribution.png", fig_hydrate_distribution, "sv_eng"),
    ("18_shutin_profile_deposit.png", fig_shutin_profile, "sv"),
    ("19_spacetime_fields.png", fig_spacetime_fields, "sv"),
    ("20_holdup_durations.png", fig_holdup_durations, "sv"),
    ("21_riser_depth_time.png", fig_riser_depth_time, "sv"),
    ("22_cloud_maps.png", fig_cloud_maps, "sv"),
    ("23_dts_thermal_waterfall.png", fig_dts_waterfall, "sv"),
    ("24_temperature_gradient.png", fig_gradient_waterfall, "sv"),
    ("25_das_flow_noise.png", fig_das_waterfall, "sv"),
    ("26_parameter_panels.png", fig_parameter_panels, "sv"),
    ("27_wellposedness_map.png", fig_wellposedness, "sv"),
]


def spacetime_outputs(sv, eng, outdir, verbose=True, save=True):
    """Render the whole space-time figure set into outdir.

    Also archives the snapshot state next to the figures (save=False when
    re-rendering from that archive) so any figure can be rebuilt later without
    repeating the transient.
    """
    os.makedirs(outdir, exist_ok=True)
    if save:
        try:
            save_state(sv, outdir)
        except Exception as exc:
            print(f"    [space-time] state not saved: {type(exc).__name__}: {exc}",
                  flush=True)
    made = []
    for name, fn, sig in FIGURES:
        try:
            p = fn(sv, eng, outdir) if sig == "sv_eng" else fn(sv, outdir)
            if p:
                made.append(p)
                if verbose:
                    print(f"    [space-time] {os.path.basename(p)}", flush=True)
        except Exception as exc:                       # never break a run for a figure
            print(f"    [space-time] {name} skipped: {type(exc).__name__}: {exc}",
                  flush=True)
    return made
