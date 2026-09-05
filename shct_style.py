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
import os
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
    #  TYPE SIZE FOR THE MEDIUM THE FIGURE IS READ IN.
    #  A figure drawn for a journal column is read at 3.5 in wide from arm's length;
    #  the same figure projected on a slide is read at 2-4 in from the back of a
    #  room, and its 8 pt tick labels land nearer 4 pt on the wall. Multi-panel
    #  figures are worst, because each panel takes a fraction of the frame. Setting
    #  SHCT_FIG_FONTSCALE scales every text element so a figure can be regenerated
    #  for the slide without redrawing it: the deck build uses ~1.8, print uses 1.0.
    try:
        _fs = float(os.environ.get("SHCT_FIG_FONTSCALE", "1.0"))
    except ValueError:
        _fs = 1.0
    if abs(_fs - 1.0) > 1e-9:
        #  Only font.size is pre-scaled here. Matplotlib builds axis labels, tick
        #  labels, titles and legends by passing the OTHER rcParams through as
        #  explicit sizes, which the Text wrapper below already scales — scaling
        #  both would land them at 1.8 x 1.8 = 3.24 times their intended size.
        mpl.rcParams["font.size"] = 10.0 * _fs
        for _k, _base in (("axes.titlesize", 11.0), ("axes.labelsize", 10.0),
                          ("xtick.labelsize", 8.5), ("ytick.labelsize", 8.5)):
            mpl.rcParams[_k] = _base
        #  legend text is built from FontProperties rather than through
        #  Text.set_fontsize, so it alone still has to be pre-scaled; suptitle does
        #  route through the wrapper and would otherwise land at 1.8 x 1.8
        mpl.rcParams["legend.fontsize"] = 8.5 * _fs
        mpl.rcParams["figure.titlesize"] = 12.0
        #  rcParams only govern text that does NOT carry an explicit size, and this
        #  project sets one on 134 call sites — 101 in shct_spacetime alone. Those
        #  ignored the scale entirely, so a "slide" rendering of the six-panel
        #  space-time figure kept 9 pt axis labels and 7.5 pt colorbar ticks while
        #  audit_deck credited the whole figure with an 18 pt base: the deck was
        #  about half as legible as it was being measured. Scale the explicit sizes
        #  at their single choke point instead of at every call site.
        #
        #  Text.__init__ takes its default size from FontProperties (already scaled
        #  through rcParams above) without routing it here, so this multiplies only
        #  sizes a caller passed deliberately — no double scaling.
        from matplotlib.text import Text as _Text
        if not getattr(_Text, "_shct_fontsize_wrapped", False):
            _orig_set_fontsize = _Text.set_fontsize

            def _scaled_set_fontsize(self, size):
                if isinstance(size, (int, float)) and not isinstance(size, bool):
                    size = size * _fs
                return _orig_set_fontsize(self, size)

            _Text.set_fontsize = _scaled_set_fontsize
            _Text.set_size = _scaled_set_fontsize
            _Text._shct_fontsize_wrapped = True

        #  thicker strokes too, or the lines vanish before the labels do
        for _k, _base in (("lines.linewidth", 1.5), ("axes.linewidth", 0.9),
                          ("xtick.major.width", 0.9), ("ytick.major.width", 0.9),
                          ("grid.linewidth", 0.8)):
            mpl.rcParams[_k] = _base * min(_fs, 1.6)
    #  FIGURE SIZE FOR THE MEDIUM, which is the lever that actually decides legibility.
    #  Type size on a slide is base_pt x (displayed_width / natural_width), so a figure
    #  drawn 13 in wide and shown in a 2.4 in frame renders its 10 pt labels at 1.8 pt
    #  however large the fonts were set: you would need a 65 pt base to recover 12 pt,
    #  which would obliterate the plot. Scaling the FIGURE down instead brings the
    #  natural width toward the frame width, so the figure is displayed near 1:1 and
    #  the type arrives at very nearly the size it was set in.
    #
    #  Most call sites pass figsize explicitly to subplots()/figure(), which overrides
    #  any rcParam, so the size is applied by wrapping those two calls rather than by
    #  setting figure.figsize. The wrapper is installed once and is idempotent.
    try:
        _sz = float(os.environ.get("SHCT_FIG_SIZESCALE", "1.0"))
    except ValueError:
        _sz = 1.0
    if abs(_sz - 1.0) > 1e-9:
        import matplotlib.pyplot as _plt
        if not getattr(_plt, "_shct_size_wrapped", False):
            def _scale(kw):
                fs = kw.get("figsize")
                if fs and len(fs) == 2:
                    kw["figsize"] = (fs[0] * _sz, fs[1] * _sz)
                return kw
            #  A smaller figure with the same number of ticks is how labels collide:
            #  the axis keeps eight tick labels while the axis itself has shrunk to
            #  40 % of its width, so they run into one another. Thinning the ticks in
            #  proportion is what keeps a shrunk figure legible rather than crowded,
            #  and it is the difference between "small" and "small and simple".
            from matplotlib.ticker import MaxNLocator
            _nb = max(3, int(round(6 * _sz)) + 1)      # ~4 ticks at 0.42, 7 at full size

            def _thin(ax):
                #  a 3-D axes carries a third axis, and pruning its ends throws on
                #  some versions, so each axis is thinned independently
                for name in ("xaxis", "yaxis", "zaxis"):
                    axis = getattr(ax, name, None)
                    if axis is None:
                        continue
                    try:
                        axis.set_major_locator(MaxNLocator(nbins=_nb, prune="both"))
                    except Exception:
                        try:
                            axis.set_major_locator(MaxNLocator(nbins=_nb))
                        except Exception:
                            pass

            def _thin_all(res):
                axes = res[1] if isinstance(res, tuple) and len(res) == 2 else None
                if axes is None:
                    return res
                try:
                    for ax in (axes.ravel() if hasattr(axes, "ravel") else [axes]):
                        _thin(ax)
                except Exception:
                    pass
                return res

            #  Thinning only through plt.subplots misses every figure built as
            #  plt.figure() + add_subplot — which is how the 3-D fields and the
            #  closure-validation charts are drawn. Their tick labels are what sets
            #  the tight-bbox floor, so those figures did not shrink AT ALL under a
            #  size scale (9.29 in at both 1.0 and 0.45) and stayed unreadable on a
            #  slide however they were placed. Wrap the Figure method as well.
            from matplotlib.figure import Figure as _Fig
            if not getattr(_Fig, "_shct_axes_wrapped", False):
                _add_sub = _Fig.add_subplot

                def _add_sub_thin(self, *a, **k):
                    ax = _add_sub(self, *a, **k)
                    _thin(ax)
                    return ax

                _Fig.add_subplot = _add_sub_thin
                _Fig._shct_axes_wrapped = True

            _f, _s = _plt.figure, _plt.subplots
            _plt.figure = lambda *a, **k: _f(*a, **_scale(k))
            #  subplots() delegates to figure(), so scaling in BOTH applies the
            #  factor twice and the figure lands at 0.42^2 = 18 % of its intended
            #  size. Only figure() carries the size; subplots() thins ticks only.
            _plt.subplots = lambda *a, **k: _thin_all(_s(*a, **k))
            _plt.rcParams["figure.figsize"] = [v * _sz for v in _plt.rcParams["figure.figsize"]]
            _plt._shct_size_wrapped = True
    return mpl.rcParams



def compact():
    """True when a figure is being rendered small, for a slide rather than a page.

    At a slide's size the LABELS, not the plot, set the figure's width: savefig
    uses a tight bounding box, so shrinking the canvas leaves the 18 pt axis and
    colorbar labels sticking out and the box simply grows back around them. A
    3-D tube view measured 9.29 in wide at every size scale for exactly this
    reason. The only way to make such a figure narrower is to give it less text,
    so drawing code asks this and uses short labels when it is true.
    """
    try:
        return abs(float(os.environ.get("SHCT_FIG_SIZESCALE", "1.0")) - 1.0) > 1e-9
    except ValueError:
        return False


def label(long_form, short_form):
    """The long label for print, the short one for a slide-sized rendering."""
    return short_form if compact() else long_form


# apply on import so a bare `import shct_style` is enough.
apply_style()
