#!/usr/bin/env python3
# =============================================================================
#  shct_style.py  —  shared plotting style + palette for the SHCT case study.
#
#  DESIGN RULE (per request): NO BLACK and NO DARK colours anywhere in any
#  generated figure, chart, graph, curve, contour or map.  This module sets the
#  global matplotlib rcParams so that EVERY foreground element (text, axis
#  labels, tick labels, spines, tick marks, legend frame, patch edges, the
#  default line-colour cycle) is drawn in a medium, clearly-coloured hue rather
#  than the matplotlib default black.  Import it (or call apply_style()) at the
#  top of every plotting / report-building script BEFORE any figure is made.
#
#      import shct_style as S
#      S.apply_style()
#
#  All named colours below are medium-value (L* ~ 0.45-0.7) — saturated, legible
#  on white, and deliberately none of them is black, near-black or a dark
#  grey/brown.
#  Author: Akosa Samuel Onyejekwe.
# =============================================================================
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

# --- the medium, non-black, non-dark palette ---------------------------------
BLUE    = "#2E5BBF"   # primary royal blue
TEAL    = "#1AA0A0"   # teal / cyan-green
ORANGE  = "#E8842B"   # warm orange
RED     = "#E0463C"   # clear red (NOT maroon / dark)
GREEN   = "#3FA65A"   # medium green
PURPLE  = "#8E5CC8"   # medium violet
AMBER   = "#E2B13C"   # golden amber
MAGENTA = "#D24A8E"   # rose magenta
BROWN   = "#B07A33"   # medium ochre/brown (for seabed / terrain lines)
SKY     = "#4FA8E0"   # light sky blue (secondary)

# foreground "ink" used for ALL text, axes, ticks, spines — a medium blue, never
# black and never a dark grey.
INK     = "#3A5BA8"
INK_HEX = INK
TITLE   = "#2E5BBF"   # figure titles
GRIDC   = "#D2DCF2"   # very light blue grid
TAN     = "#E7D7B6"   # light tan terrain fill (light, not dark)
TAN_EDGE = BROWN      # terrain outline
HYDFILL = "#F6D6D2"   # light rose hydrate-stability fill
SLUGFILL = "#FBE2DD"  # light rose intermittent-flow band
CRIT    = "#D24A8E"   # magenta — Phi_SH = 1 critical contour (replaces black)

# ordered cycle used for multi-series line plots
PALETTE = [BLUE, ORANGE, TEAL, GREEN, PURPLE, RED, AMBER, MAGENTA, SKY, BROWN]

# --- no-black / no-dark COLORMAPS for heatmaps, contours and 3-D surfaces -----
#  Every colour stop below is light or medium (L* ~ 0.5-0.95): NONE is black,
#  near-black or a dark hue.  These replace cividis / inferno / viridis /
#  coolwarm / RdYlBu_r (all of which run into black or dark ends) everywhere a
#  field is colour-mapped, so the zero/low end of every map is LIGHT, not black.
#  darkest stop used anywhere is the saturated royal blue at the low end and the
#  saturated red at the high end (relative luminance ~0.29-0.30): both are clearly
#  coloured, neither is a dark or near-black hue.
# STRONG, high-contrast, multi-hue map (selected): blue -> cyan -> green -> amber ->
# orange -> red. Saturated and bold so values that cluster in the mid-range get
# DISTINCT strong colours (the soft single-hue maps washed those out); every stop is
# still medium/saturated — none is black, near-black or dark. Used with smooth
# (gouraud) shading for EVERY colour-mapped field: holdup, Φ_SH, deposit, temperature,
# velocity. Low = strong blue, high = strong red (intuitive: red = hot / critical / full).
#  Built as a HUE SWEEP rather than as a handful of RGB waypoints. Interpolating
#  between waypoints in RGB desaturates every midpoint — the blue/cyan and the
#  green/yellow crossings lost up to a fifth of their brightness and a tenth of
#  their saturation, which is what made the greens read as muddy olive. Sweeping
#  the HUE at full saturation and full value instead keeps every intermediate as
#  vivid as the anchors: pure blue -> cyan -> green -> yellow -> orange -> red.
def _hue_sweep(h0=232.0, h1=0.0, n=256, sat=1.0, val=1.0):
    import colorsys
    return [colorsys.hsv_to_rgb(((h0 + (h1 - h0) * i / (n - 1)) % 360.0) / 360.0,
                                sat, val) for i in range(n)]


_STRONG = _hue_sweep()          # MAX-vivid rainbow, fully saturated throughout
#  ONE scheme for EVERY colour-mapped field in this project (per request): the
#  distributed-sensing rainbow — saturated blue at the low end, through cyan and
#  green, to yellow/orange/red at the high end. Deep blue reads as "cold / empty /
#  quiet" and red as "hot / full / critical" without any further explanation, the
#  gradient is smooth and continuous, and no stop is black, near-black or dark.
#  Holdup, temperature, deposit, velocity, Phi_SH, the waterfalls and the gradient
#  maps all share it, so a colour means the same thing across the whole figure set.
_DTS = list(_STRONG)
_GRAD = list(_STRONG)
CMAP_SEQ  = LinearSegmentedColormap.from_list("shct_seq",  _STRONG)   # holdup, velocity
CMAP_HEAT = LinearSegmentedColormap.from_list("shct_heat", _STRONG)   # deposit thickness
CMAP_TEMP = LinearSegmentedColormap.from_list("shct_temp", _STRONG)   # temperature
CMAP_DIV  = LinearSegmentedColormap.from_list("shct_div",  _STRONG)   # Φ_SH (critical contour drawn on top)
CMAP_DTS  = LinearSegmentedColormap.from_list("shct_dts",  _DTS)      # DTS/DAS-style waterfalls
CMAP_GRAD = LinearSegmentedColormap.from_list("shct_grad", _GRAD)     # signed gradient maps

for _cm in (CMAP_SEQ, CMAP_HEAT, CMAP_TEMP, CMAP_DIV, CMAP_DTS, CMAP_GRAD):
    try:
        mpl.colormaps.register(_cm, force=True)       # register by name (shct_seq, ...)
    except Exception:                                 # pragma: no cover (old mpl)
        try:
            mpl.cm.register_cmap(name=_cm.name, cmap=_cm)
        except Exception:
            pass

# --- smooth rendering of a coarse field ---------------------------------------
def smooth_field(F, x=None, y=None, target=560, order=3):
    """Resample a 2-D field onto a fine grid with a smooth cubic interpolant.

    A 70-cell transport grid drawn directly shows its cells: hard vertical
    banding and stepped edges, rather than the continuous field the numbers
    describe. This upsamples for RENDERING only -- the data is untouched, in the
    same sense that a contour plot draws smooth contours through coarse samples.

    Returns (F_fine, x_fine, y_fine); x and y may be None if only the field is
    wanted. Falls back to the original arrays when SciPy is unavailable, so a
    figure never fails for want of smoothing.
    """
    import numpy as _np
    F = _np.asarray(F, float)
    if F.ndim != 2 or F.size == 0:
        return F, x, y
    ny, nx = F.shape
    zy = max(1.0, float(target) / max(ny, 1))
    zx = max(1.0, float(target) / max(nx, 1))
    if zy <= 1.0 and zx <= 1.0:
        return F, x, y
    try:
        from scipy.ndimage import zoom as _zoom
        #  fill non-finite cells before interpolating, else they smear
        Ff = _np.array(F, dtype=float, copy=True)
        bad = ~_np.isfinite(Ff)
        if bad.any():
            Ff[bad] = _np.nanmedian(Ff[~bad]) if (~bad).any() else 0.0
        out = _zoom(Ff, (zy, zx), order=order, mode="nearest", grid_mode=False)
        #  the interpolant can overshoot at sharp fronts; hold it to the data range
        out = _np.clip(out, _np.nanmin(F), _np.nanmax(F))
    except Exception:
        return F, x, y
    fy, fx = out.shape
    xf = _np.linspace(float(_np.min(x)), float(_np.max(x)), fx) if x is not None else None
    yf = _np.linspace(float(_np.min(y)), float(_np.max(y)), fy) if y is not None else None
    return out, xf, yf


# --- overlapping-text detector -------------------------------------------------
def find_text_overlaps(fig, min_overlap_frac=0.18, ignore_empty=True):
    """Return the pairs of text artists whose drawn boxes overlap.

    Text that lands on top of other text is invisible in the source and glaring on
    the page. Every drawn text artist has a bounding box in display coordinates,
    so the collisions can simply be measured. Returns a list of
    (text_a, text_b, overlap_fraction) with the fraction relative to the SMALLER
    box, worst first; an empty list means nothing collides.

    min_overlap_frac ignores the incidental one- or two-pixel touches that tight
    layouts produce and that no reader would notice.
    """
    fig.canvas.draw()                      # boxes only exist once drawn
    renderer = fig.canvas.get_renderer()
    items = []
    for ax in fig.get_axes():
        cand = list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        leg = ax.get_legend()
        if leg is not None:
            cand += list(leg.texts)
        items += cand
    items += list(fig.texts)

    boxes = []
    for t in items:
        if t is None:
            continue
        try:
            if ignore_empty and not str(t.get_text()).strip():
                continue
            if not t.get_visible():
                continue
            boxes.append((t, t.get_window_extent(renderer=renderer)))
        except Exception:
            continue

    hits = []
    for i in range(len(boxes)):
        ta, ba = boxes[i]
        for j in range(i + 1, len(boxes)):
            tb, bb = boxes[j]
            x0 = max(ba.x0, bb.x0); x1 = min(ba.x1, bb.x1)
            y0 = max(ba.y0, bb.y0); y1 = min(ba.y1, bb.y1)
            if x1 <= x0 or y1 <= y0:
                continue
            inter = (x1 - x0) * (y1 - y0)
            small = max(min(ba.width * ba.height, bb.width * bb.height), 1e-9)
            frac = inter / small
            if frac >= min_overlap_frac:
                hits.append((ta, tb, frac))
    hits.sort(key=lambda h: -h[2])
    return hits


def report_text_overlaps(fig, name="figure", raise_on=None):
    """Print any overlapping text in `fig`. Returns the number of collisions."""
    try:
        hits = find_text_overlaps(fig)
    except Exception:
        return 0
    for ta, tb, frac in hits:
        a = " ".join(str(ta.get_text()).split())[:44]
        b = " ".join(str(tb.get_text()).split())[:44]
        print(f"    [overlap] {name}: {frac*100:.0f} % — {a!r} over {b!r}",
              flush=True)
    if hits and raise_on:
        raise RuntimeError(f"{name}: {len(hits)} overlapping text item(s)")
    return len(hits)


# RGB 0-255 tuples (handy for python-docx RGBColor or other consumers)
def _rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))

RGB = {name: _rgb(val) for name, val in
       dict(BLUE=BLUE, TEAL=TEAL, ORANGE=ORANGE, RED=RED, GREEN=GREEN,
            PURPLE=PURPLE, AMBER=AMBER, MAGENTA=MAGENTA, BROWN=BROWN,
            SKY=SKY, INK=INK, TITLE=TITLE).items()}


def apply_style():
    """Install the no-black / no-dark style into matplotlib's global rcParams."""
    mpl.rcParams.update({
        # --- foreground: every default-black element recoloured to medium ink ---
        "text.color":        INK,
        "axes.labelcolor":   INK,
        "axes.edgecolor":    INK,
        "axes.titlecolor":   TITLE,
        "xtick.color":       INK,
        "ytick.color":       INK,
        "xtick.labelcolor":  INK,
        "ytick.labelcolor":  INK,
        "patch.edgecolor":   INK,
        "hatch.color":       INK,
        "legend.edgecolor":  INK,
        # --- backgrounds stay white (never a dark theme) ---
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
        "savefig.facecolor": "white",
        # --- grid ---
        "axes.grid":         False,
        "grid.color":        GRIDC,
        "grid.alpha":        0.5,
        # --- the line / marker colour cycle (no black) ---
        "axes.prop_cycle":   mpl.cycler(color=PALETTE),
        "lines.color":       BLUE,
        # --- legends: opaque so a legend NEVER lets a bar/curve show through, and
        #     savefig in 'tight' mode so a legend placed OUTSIDE the axes (the
        #     default placement used throughout, so text never overlaps data) is
        #     never clipped. ---
        "legend.framealpha": 1.0,
        "legend.facecolor":  "white",
        "legend.fancybox":   True,
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.06,
    })
    return mpl.rcParams


# apply on import so a bare `import shct_style` is enough.
apply_style()
