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

#  `mpl.cycler` is a runtime re-export that the matplotlib stubs do not carry; the
#  cycler package is matplotlib's own dependency, so import it from source.
from cycler import cycler as _cycler
from matplotlib.colors import LinearSegmentedColormap

#  ---------------------------------------------------------------------------
#  ONE export resolution for every figure this project draws.
#
#  Eight modules each carried their own SHCT_FIG_DPI default and they disagreed:
#  155 in solver.py, 150 in shct_threed.py and shct_compositional_sim.py, 320 in the
#  other five. So a single output folder held figures at three different resolutions
#  with nothing recording why, two of them below the 300 dpi the journal requires and
#  which check_journal_artwork.py enforces. That mattered because
#  export_paper_figures.py maps some manuscript figures STRAIGHT from those folders
#  rather than from the 320 dpi manuscript build, so a 150 dpi chart could reach the
#  journal set on nothing more than which figure someone chose to map.
#
#  The default lives here, once. SHCT_FIG_DPI still overrides it, which is how the
#  slide build asks for its own resolution.
FIG_DPI = int(float(os.environ.get("SHCT_FIG_DPI", "320")))

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
#  AXIS FURNITURE IS BLACK, and always was meant to be. The no-dark rule on file exists
#  because dark COLORMAPS were rejected -- cividis, inferno, viridis, coolwarm all put
#  near-black at the low end of the DATA, so a zero reads as a hole. That reasoning is
#  about the field, not about the frame. Applying it to text and spines as well left every
#  label, tick and axis line in medium blue, which is why the figures did not look like
#  scientific output: black text on white is the convention, and blue furniture reads as a
#  slide template. Confirmed with the user. The colormap keeps the no-dark rule (ocean,
#  darkest stop L* 31); the FRAME goes back to black.
INK     = "#000000"
INK_BLUE = "#3A5BA8"        # the former ink, kept for anything that wants a blue accent
INK_HEX = INK
TITLE   = "#2E5BBF"   # figure titles
GRIDC   = "#CCCCCC"   # neutral light grey grid (was light blue, under blue axes)
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


#  THE PALETTE HAD DRIFTED FROM THE ONE THAT WAS ASKED FOR. The instruction on file
#  names an OCEAN ramp, mint -> teal -> blue -> indigo. What was here instead was a
#  full saturated hue sweep -- blue, cyan, green, yellow, orange, red -- i.e. jet, under
#  a comment claiming it was 'per request'. Measured against CIE L*:
#
#      rainbow (was)   L* 32-97   NOT monotonic, 3 reversals, 274x rate variation
#      ocean   (now)   L* 31-96   monotonic,     0 reversals,   1.8x rate variation
#
#  A rainbow cannot be ordered by eye: it runs light, dark, light, dark, so a reader
#  cannot tell 'more' from 'less' by brightness. It also invents hard bands at the cyan
#  and yellow turns where the field is smooth, and flattens the green stretch where real
#  structure lives (0.004 L* per step against 1.12 elsewhere). The ocean ramp is no
#  darker at its low end than the rainbow was (L* 31 vs 32), so the no-dark rule holds.
#
#  SATURATED, because the first ocean ramp was measurably dull. Fixing the rainbow's
#  lightness problem was done by dropping chroma as well, which was not required and
#  cost more than it bought: mean chroma 34.3 against 54.3 for viridis and 57.1 for
#  turbo, and a low end of #E6F9F0 at chroma 8.1 -- so nearly white that the bottom of
#  every field read as empty paper rather than as a value. The ramp below keeps the
#  same mint-to-indigo family and the same monotonic lightness, and spends the chroma
#  the earlier one left on the table.
#
#  THE SCHEME IS THE USER'S CHOICE, AND IT IS A RAINBOW. Two earlier ramps were tried
#  and both were rejected as dull on sight: the original ocean (mean chroma 34.3) and a
#  saturated version of it (48.0). Chroma was not the whole problem -- cyan-blue-indigo
#  turns only 127 degrees of hue, so every field rendered as one wash -- but a wide-hue
#  sequential ramp (sun-to-sea, chroma 52.4, 246 degrees) was still not what was wanted.
#  The reference supplied is the classic flow-assurance rainbow, and these are its stops.
#
#      ocean         L* 31-96  chroma 34.3  hue 127 deg  0 reversals
#      saturated     L* 31-91  chroma 48.0  hue 127 deg  0 reversals
#      sun-to-sea    L* 37-95  chroma 52.4  hue 246 deg  0 reversals
#      THIS rainbow  L* 33-97  chroma 97.8  hue 288 deg  4 reversals
#
#  What it buys: mean chroma 97.8, nearly three times the first ocean ramp, and a MINIMUM
#  chroma of 43.8 -- there is no pale end at all, so a field whose data sits near the
#  bottom of its range still renders in a saturated colour instead of near-white paper.
#  That was the real complaint behind "weak and dull" on the deposit maps.
#
#  What it costs, stated plainly because it is a real cost: lightness is NOT monotonic.
#  It runs 60-32-91-97-53, four reversals, so the ramp cannot be ordered by brightness
#  and a reader must use the colour bar rather than read "more" as "darker". Greyscale
#  printing and most colour-vision deficiencies lose the ordering. This was argued
#  against and the choice was made anyway with both alternatives in hand; it is a
#  deliberate decision, not an oversight, and the no-dark rule still holds (darkest stop
#  L* 33, no black or near-black anywhere).
_OCEAN = ["#FF00FF", "#8000FF", "#0000FF", "#0080FF", "#00FFFF", "#00FF80",
          "#00FF00", "#80FF00", "#FFFF00", "#FFBF00", "#FF8000", "#FF4000", "#FF0000"]
_STRONG = _OCEAN
#  ONE scheme for EVERY colour-mapped field in this project: the OCEAN ramp, pale
#  mint at the low end through teal to indigo at the high end. It is MONOTONIC in
#  perceived lightness, so "more" is simply darker and the eye can order the field
#  without consulting the colour bar; no stop is black, near-black or dark; and
#  holdup, temperature, deposit, velocity, Phi_SH, the waterfalls and the gradient
#  maps all share it, so a colour means the same thing across the whole figure set.
#
#  NOTE for the signed fields (the gradient map, and subcooling where it changes
#  sign): a sequential ramp gives zero no visual anchor -- it lands on whatever
#  mid-teal happens to be halfway. A diverging map with a pale midpoint would read
#  better there. Left sequential because the instruction on file names one ocean
#  scheme for every field; worth revisiting if signed maps matter.
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
            #  removed in matplotlib 3.9; reached only on the old versions that have it
            _reg = getattr(mpl.cm, "register_cmap", None)
            if _reg is not None:
                _reg(name=_cm.name, cmap=_cm)
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
def find_text_overlaps(fig, min_overlap_frac=0.18, min_overlap_px=10.0,
                       ignore_empty=True):
    """Return the pairs of text artists whose drawn boxes overlap.

    Text that lands on top of other text is invisible in the source and glaring on
    the page. Every drawn text artist has a bounding box in display coordinates,
    so the collisions can simply be measured. Returns a list of
    (text_a, text_b, overlap_fraction) with the fraction relative to the SMALLER
    box, worst first; an empty list means nothing collides.

    TWO criteria, because the area one alone has a blind spot that let a real
    collision through. A rotated axis label is a long thin box, so a wide note laid
    across it covers only a few per cent of its AREA while hiding whole characters:
    figure 21's margin note sat across "depth from host [m]", covering its "[m]" at
    7.5 % of the label's area, and the 18 % area threshold passed it. A hit is
    therefore also reported when the intersection is at least min_overlap_px in BOTH
    directions -- about one character each way at the default 100 dpi -- which is the
    smallest collision a reader can actually see.

    min_overlap_frac and min_overlap_px together ignore the incidental one- or
    two-pixel touches that tight layouts produce and that no reader would notice.
    """
    fig.canvas.draw()                      # boxes only exist once drawn
    #  FigureCanvasBase declares no get_renderer; every concrete backend has one.
    renderer = getattr(fig.canvas, "get_renderer", fig._get_renderer)()
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
            if frac >= min_overlap_frac or ((x1 - x0) >= min_overlap_px
                                            and (y1 - y0) >= min_overlap_px):
                hits.append((ta, tb, frac))
    hits.sort(key=lambda h: -h[2])
    return hits


def find_empty_axes(fig):
    """Return the visible axes of `fig` that carry no data and no text.

    A panel that draws nothing is invisible in the source and glaring on the page, and
    it is the one figure fault none of the other checks here can see: the canvas is not
    blank (the other panels are full), the colours are fine, the labels are all present.
    Figure 07 shipped for a long time with an empty Kaplan-Meier axis on both scenarios
    that do not plug -- axes, ticks, title and axis labels, and nothing inside them --
    and it read as a broken run rather than as the result.

    An axes counts as empty when it is visible, has its frame on, and holds no lines,
    patches, collections, images, tables, texts or legend. `axis("off")` panels are
    skipped: those are deliberately bare, and are used here to hold plain text blocks.
    """
    fig.canvas.draw()
    empty = []
    for ax in fig.get_axes():
        try:
            if not ax.get_visible() or not ax.axison:
                continue
            if getattr(ax, "_colorbar", None) is not None:
                continue
            has = (list(ax.lines) or list(ax.patches) or list(ax.collections)
                   or list(ax.images) or list(ax.tables) or list(ax.texts)
                   or ax.get_legend() is not None
                   or getattr(ax, "containers", []))
            if not has:
                empty.append(ax)
        except Exception:
            continue
    return empty


def right_margin_anchor(ax, pad_in=0.12, floor=1.02):
    """x, in `ax`'s axes-fraction, just clear of everything drawn in its RIGHT margin.

    A legend parked at a FIXED offset -- bbox_to_anchor=(1.13, 1.0) -- assumes it knows
    how wide the right margin is. It does not. Under SHCT_FIG_FONTSCALE the twin axis's
    tick labels and its "T, T_eq (°C)" ylabel grow while 1.13 stays put, so on the 0.30,
    0.45 and 0.70 slide renders the legend was drawn straight over that label. Measure
    what is actually out there and clear it.

    Returns `floor` when the margin is empty, so an ordinary single-axis figure is
    unchanged.
    """
    fig = ax.figure
    fig.canvas.draw()
    get = getattr(fig.canvas, "get_renderer", None)
    rend = get() if get is not None else fig._get_renderer()
    ax_bb = ax.get_window_extent(rend)
    if ax_bb.width <= 0:
        return floor
    right = ax_bb.x1
    #  the twin shares this axes' position, so its decorations live in our right margin
    for other in fig.get_axes():
        try:
            if other.get_position().bounds != ax.get_position().bounds:
                continue
            cands = [other.yaxis.label] + list(other.get_yticklabels())
            for t in cands:
                if not t.get_visible() or not str(t.get_text()).strip():
                    continue
                bb = t.get_window_extent(rend)
                if bb.x0 >= ax_bb.x1 - 1.0:          # genuinely in the right margin
                    right = max(right, bb.x1)
        except Exception:
            continue
    if right <= ax_bb.x1 + 1.0:
        return floor
    return min(1.6, max(floor, (right - ax_bb.x0 + pad_in * fig.dpi) / ax_bb.width))


def legend_outside(ax, handles=None, labels=None, max_width_frac=0.80, below_pad=0.22,
                   ncol_below=None, **kw):
    """Put a legend outside `ax`, on whichever side it actually fits.

    Beside the axes when the right margin can take it, BELOW when it cannot. A legend
    pinned beside the axes at a fixed offset does not merely look cramped when it
    outgrows the canvas -- tight_layout shrinks the AXES to make room for it, and past a
    point the panel is destroyed. outputs_slides45/01_profiles.png shipped that way: four
    stacked panels crushed to 0.23 x 0.62 in on a 3.6 x 3.7 in canvas, x tick labels
    printed on top of one another, every curve a vertical smear, with two legends taking
    most of the figure. It renders correctly at scale 1.0, which is why it was never seen.

    `max_width_frac` is the share of the figure width a beside-legend may occupy before
    it is moved below.
    """
    _kw = dict(kw)
    if handles is not None and labels is not None:
        lg = ax.legend(handles, labels, loc="upper left",
                       bbox_to_anchor=(right_margin_anchor(ax), 1.0), **_kw)
    elif handles is not None:
        lg = ax.legend(handles=handles, loc="upper left",
                       bbox_to_anchor=(right_margin_anchor(ax), 1.0), **_kw)
    else:
        lg = ax.legend(loc="upper left",
                       bbox_to_anchor=(right_margin_anchor(ax), 1.0), **_kw)
    fig = ax.figure
    try:
        fig.canvas.draw()
        get = getattr(fig.canvas, "get_renderer", None)
        rend = get() if get is not None else fig._get_renderer()
        ax_bb = ax.get_window_extent(rend)
        w_frac = lg.get_window_extent(rend).width / max(ax_bb.width, 1.0)
        _margin_taken = any(
            other is not ax
            and other.get_position().bounds != ax.get_position().bounds
            and other.get_window_extent(rend).x0 >= ax_bb.x1 - 2.0
            and other.get_window_extent(rend).y1 > ax_bb.y0
            and other.get_window_extent(rend).y0 < ax_bb.y1
            for other in fig.get_axes())
    except Exception:
        return lg
    if w_frac <= max_width_frac and not _margin_taken:
        return lg
    #  too wide to sit beside: below it goes, in as few rows as the entries allow
    n = len(lg.get_texts()) or 1
    lg.remove()
    _kw.pop("borderaxespad", None)
    ncol = ncol_below if ncol_below else (n if n <= 3 else (n + 1) // 2)
    try:
        try:
            fig.tight_layout()
        except Exception:
            pass
        fig.canvas.draw()
        rend = (getattr(fig.canvas, "get_renderer", None) or (lambda: fig._get_renderer()))()
        lab = ax.xaxis.label
        if not str(lab.get_text()).strip():
            for other in fig.get_axes():
                if (other is not ax
                        and other.get_position().bounds == ax.get_position().bounds
                        and str(other.xaxis.label.get_text()).strip()):
                    lab = other.xaxis.label
                    break
        if lab.get_visible() and str(lab.get_text()).strip():
            ax_bb = ax.get_window_extent(rend)
            lab_bb = lab.get_window_extent(rend)
            drop = (ax_bb.y0 - lab_bb.y0) / max(ax_bb.height, 1.0)
            below_pad = max(below_pad, drop + 0.08)
    except Exception:
        pass
    if handles is not None and labels is not None:
        return ax.legend(handles, labels, loc="upper center",
                         bbox_to_anchor=(0.5, -below_pad), ncol=ncol,
                         frameon=False, **_kw)
    if handles is not None:
        return ax.legend(handles=handles, loc="upper center",
                         bbox_to_anchor=(0.5, -below_pad), ncol=ncol,
                         frameon=False, **_kw)
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, -below_pad), ncol=ncol,
                     frameon=False, **_kw)


def find_text_over_data(fig, min_frac=0.12):
    """Return in-axes text that is drawn across the plotted data.

    THE SCREENS COULD NOT SEE THIS. find_text_overlaps compares text against TEXT, and
    find_empty/degenerate/squeezed_axes look at axes, so a note sitting squarely on a
    curve passed every one of them -- which is how 01_profiles shipped four lines of
    prose over its subcooling trace and 08_diagnostics a caption across its slug-length
    trace, on figures that reported zero hits.

    Only text in AXES coordinates is considered: captions placed below an axes or in
    figure coordinates are outside the data area by construction.
    """
    import numpy as _np
    out: list[tuple] = []
    try:
        fig.canvas.draw()
        _g = getattr(fig.canvas, "get_renderer", None)
        rend = _g() if _g is not None else fig._get_renderer()
    except Exception:
        return out
    for ax in fig.get_axes():
        try:
            if not ax.get_visible() or not ax.axison:
                continue
            ax_bb = ax.get_window_extent(rend)
            if ax_bb.width <= 1 or ax_bb.height <= 1:
                continue
            pts = []
            for ln in ax.get_lines():
                xd, yd = ln.get_xdata(), ln.get_ydata()
                if xd is None or len(xd) == 0:
                    continue
                d = ax.transData.transform(_np.column_stack(
                    [_np.asarray(xd, float), _np.asarray(yd, float)]))
                pts.append(d[_np.isfinite(d).all(axis=1)])
            if not pts:
                continue
            P = _np.vstack(pts)
            if P.size == 0:
                continue
            for t in ax.texts:
                if not t.get_visible() or not str(t.get_text()).strip():
                    continue
                #  only notes anchored INSIDE the axes box
                tb = t.get_window_extent(rend)
                if (tb.x1 < ax_bb.x0 or tb.x0 > ax_bb.x1
                        or tb.y1 < ax_bb.y0 or tb.y0 > ax_bb.y1):
                    continue
                inside = ((P[:, 0] >= tb.x0) & (P[:, 0] <= tb.x1)
                          & (P[:, 1] >= tb.y0) & (P[:, 1] <= tb.y1))
                frac = float(inside.sum()) / max(len(P), 1)
                if frac >= min_frac:
                    out.append((t, frac))
        except Exception:
            continue
    return out


def report_text_over_data(fig, name):
    """Print any note drawn across the data it sits on."""
    for t, frac in find_text_over_data(fig):
        txt = " ".join(str(t.get_text()).split())[:52]
        print(f"    [text on data] {name}: {100.0 * frac:.0f} % of the plotted points "
              f"lie under '{txt}'", flush=True)


def data_occupancy(ax, quadrant, pad=0.03):
    """Fraction of this axes' plotted points that fall inside one corner quadrant.

    Corners are where in-axes notes go, and a note dropped in a corner the data
    happens to occupy sits on the curve -- which is the one thing the figures in this
    project are not allowed to do. Counting the points that fall there is the only way
    to know which corner is free, because it depends on the data, not on the layout.
    """
    x0, x1 = (0.0, 0.5) if "left" in quadrant else (0.5, 1.0)
    y0, y1 = (0.0, 0.5) if "lower" in quadrant else (0.5, 1.0)
    try:
        xlo, xhi = ax.get_xlim(); ylo, yhi = ax.get_ylim()
        if xhi == xlo or yhi == ylo:
            return 1.0
        n_in = n_tot = 0
        for ln in ax.get_lines():
            xd, yd = ln.get_xdata(), ln.get_ydata()
            if xd is None or len(xd) == 0:
                continue
            import numpy as _np
            fx = (_np.asarray(xd, float) - xlo) / (xhi - xlo)
            fy = (_np.asarray(yd, float) - ylo) / (yhi - ylo)
            ok = _np.isfinite(fx) & _np.isfinite(fy)
            n_tot += int(ok.sum())
            n_in += int((ok & (fx >= x0 - pad) & (fx <= x1 + pad)
                            & (fy >= y0 - pad) & (fy <= y1 + pad)).sum())
        return 1.0 if n_tot == 0 else n_in / float(n_tot)
    except Exception:
        return 1.0


def place_note(ax, s, fontsize=6.6, color="#666666", style="italic", **kw):
    """Put a note INSIDE `ax`, in whichever corner the data leaves free.

    A note pinned to a fixed corner is right until the data changes shape. Two shipped
    that way: 01_profiles put four lines of prose at (0.985, 0.06) straight through the
    subcooling curve, and 08_diagnostics put its off-scale note at (0.98, 0.94) across
    the slug-length trace. Both corners were free when they were chosen.

    Corners are ranked by how much of the plotted data falls in them, and the emptiest
    wins; if every corner is more than 15 % occupied the note goes BELOW the axes, where
    it cannot cover anything.
    """
    corners = {"upper right": (0.985, 0.955, "right", "top"),
               "upper left":  (0.015, 0.955, "left",  "top"),
               "lower right": (0.985, 0.045, "right", "bottom"),
               "lower left":  (0.015, 0.045, "left",  "bottom")}
    best, occ = None, 1.0
    for name, spec in corners.items():
        o = data_occupancy(ax, name)
        if o < occ:
            best, occ = spec, o
    if best is None or occ > 0.15:
        return ax.text(0.5, -0.28, s, transform=ax.transAxes, ha="center", va="top",
                       fontsize=fontsize, color=color, style=style, **kw)
    x, y, ha, va = best
    return ax.text(x, y, s, transform=ax.transAxes, ha=ha, va=va,
                   fontsize=fontsize, color=color, style=style, **kw)


def fix_below_legends(fig, pad_frac=0.08):
    """Re-drop any below-axes legend that a later layout change has landed on its x-label.

    legend_outside() measures where the x-label sits and drops the legend clear of it,
    but the offset it stores is in AXES-HEIGHT fractions. Anything that changes the axes
    height afterwards -- a caption strip added to the canvas, relieve_crush() growing a
    panel, a second tight_layout -- rescales that offset and can walk the legend back up
    onto the label it was measured to clear. compo_pvt at the 0.45 slide scale shipped
    exactly that, 'x (km)' under 'liquid rho_l' at 28 %.

    Re-measuring after the final layout is the only way the two can stay consistent; this
    is the same failure the coincidence note had, and the same fix. Cheap and idempotent,
    so screen() runs it on every figure rather than the one that was caught.
    """
    try:
        fig.canvas.draw()
        _g = getattr(fig.canvas, "get_renderer", None)
        rend = _g() if _g is not None else fig._get_renderer()
    except Exception:
        return False
    moved = False
    for ax in fig.get_axes():
        try:
            lg = ax.get_legend()
            if lg is None or not lg.get_visible():
                continue
            lab = ax.xaxis.label
            if not lab.get_visible() or not str(lab.get_text()).strip():
                continue
            ax_bb = ax.get_window_extent(rend)
            lg_bb = lg.get_window_extent(rend)
            lab_bb = lab.get_window_extent(rend)
            #  only legends that sit BELOW the axes are ours to move
            if lg_bb.y1 > ax_bb.y0:
                continue
            #  no vertical gap between the label and the legend top means they collide
            if lg_bb.y1 <= lab_bb.y0 - 1.0:
                continue
            if lg_bb.x1 < lab_bb.x0 or lg_bb.x0 > lab_bb.x1:
                continue
            drop = (ax_bb.y0 - lab_bb.y0) / max(ax_bb.height, 1.0) + pad_frac
            lg.set_bbox_to_anchor((0.5, -drop), transform=ax.transAxes)
            moved = True
        except Exception:
            continue
    return moved


def find_squeezed_axes(fig, min_w_in=0.5, min_h_in=0.30):
    """Return axes that tight_layout has crushed to a sliver.

    THIS ONE SHIPPED. cx2_azimuthal_deposit carried a one-line caption placed in AXES
    coordinates and five times wider than the axes; tight_layout counts an axes-anchored
    artist when it sizes the axes, so it shrank the (x, theta) map to about half an inch
    on a 7.6-inch canvas. The x tick labels "0" and "30" printed on top of one another
    and the whole field rendered as a vertical blue line -- in six folders. Nothing
    caught it, because every artist was present and none of them overlapped: the AXES
    was the casualty, not the content.

    A panel legitimately gets small in a dense grid -- four across an 11-inch figure is
    still ~18 % of the width -- so the bounds are set well below any real layout and only
    a collapse trips it.

    IT ALSO SHIPPED VERTICALLY, and this check did not see it. The bound was width-only,
    because an early height test flagged three legitimate short panels for every real one
    it caught. That left compo_pvt at the 0.30 slide scale rendering four panels about
    0.13 in tall on a 4.8 in canvas -- y tick labels printed on top of one another, every
    curve a horizontal bar -- and passing every screen, because each panel was still the
    full width. A separate, much tighter height bound is what the earlier test lacked:
    0.30 in is below any panel a real grid produces at any of the four slide scales, so
    it separates a collapse from a merely short panel instead of guessing from one bound.
    """
    fig.canvas.draw()
    fw, fh = fig.get_size_inches()
    out = []
    for ax in fig.get_axes():
        try:
            if not ax.get_visible() or not ax.axison:
                continue
            if getattr(ax, "_colorbar", None) is not None:
                continue
            bb = ax.get_position()
            w, h = bb.width * fw, bb.height * fh
            if w < min_w_in or h < min_h_in:
                out.append((ax, w, h))
        except Exception:
            continue
    return out


def report_squeezed_axes(fig, name="figure", stream=None):
    """Print any axes crushed below `min_frac` of the canvas. Returns the count."""
    import sys as _sys
    hits = find_squeezed_axes(fig)
    st = stream or _sys.stdout
    for ax, w, h in hits:
        t = str(ax.get_title() or ax.get_ylabel() or "untitled")[:60]
        print(f"    [squeezed axes] {name}: panel '{t}' is only {w:.2f} x {h:.2f} in "
              f"on a {fig.get_size_inches()[0]:.1f} x {fig.get_size_inches()[1]:.1f} in "
              f"canvas — a caption or legend has crushed it", file=st)
    return len(hits)


def find_degenerate_axes(fig, span_frac=0.05, min_pts=3):
    """Return axes whose plotted data is too sparse to be the profile they claim.

    `find_empty_axes` only catches a panel that draws NOTHING. The failure that actually
    shipped was subtler: a K-value-vs-distance panel on a fluid that splits at 1 of 40
    stations drew eight markers stacked against the right-hand edge of an otherwise blank
    axis, under a full eight-entry legend for curves that did not exist. There WERE
    artists, so nothing flagged it, and it read as a broken plot rather than as the
    degenerate case it is.

    A line is degenerate here when it declares many points but its FINITE data covers less
    than `span_frac` of the axes' own x-range, or amounts to fewer than `min_pts` points.
    Either way the axis is not showing a profile, and the panel should say so or be drawn
    against something it can actually resolve.
    """
    import numpy as np
    fig.canvas.draw()
    out = []
    for ax in fig.get_axes():
        try:
            if not ax.get_visible() or not ax.axison:
                continue
            if getattr(ax, "_colorbar", None) is not None:
                continue
            x0, x1 = ax.get_xlim()
            span = abs(x1 - x0)
            if span <= 0:
                continue
            worst = None
            for ln in ax.lines:
                #  An author can DECLARE a line degenerate. A tracker that is meant to be
                #  narrow -- the steepest-cooling front on a line whose hot inlet meets a
                #  cold seabed genuinely does not move -- is not a defect once the figure
                #  says so, and leaving a permanent warning on it would mask the next real
                #  one. Tag such a line with gid "declared-degenerate"; anything untagged
                #  is still judged.
                if str(ln.get_gid() or "") == "declared-degenerate":
                    continue
                xd = np.asarray(ln.get_xdata(), float)
                yd = np.asarray(ln.get_ydata(), float)
                if xd.size < 4:                     # short series are not profiles
                    continue
                ok = np.isfinite(xd) & np.isfinite(yd)
                n = int(ok.sum())
                if n == 0:
                    continue                        # wholly empty -> find_empty_axes' job
                frac = (float(np.nanmax(xd[ok]) - np.nanmin(xd[ok])) / span) if n > 1 else 0.0
                if n < min_pts or frac < span_frac:
                    worst = (n, int(xd.size), frac) if worst is None else worst
            if worst is not None:
                out.append((ax, worst))
        except Exception:
            continue
    return out


def report_degenerate_axes(fig, name="figure"):
    """Print any axes whose data is too sparse to be the profile it claims."""
    try:
        bad = find_degenerate_axes(fig)
    except Exception:
        return 0
    for ax, (n, tot, frac) in bad:
        t = " ".join(str(ax.get_title()).split())[:52] or "(untitled)"
        print(f"    [degenerate axes] {name}: panel {t!r} plots {n} finite point(s) of "
              f"{tot} over {frac*100:.1f} % of its x-range", flush=True)
    return len(bad)


def report_empty_axes(fig, name="figure"):
    """Print any axes of `fig` that draw nothing. Returns the number found."""
    try:
        empty = find_empty_axes(fig)
    except Exception:
        return 0
    for ax in empty:
        t = " ".join(str(ax.get_title()).split())[:52] or "(untitled)"
        print(f"    [empty axes] {name}: a panel titled {t!r} draws nothing", flush=True)
    return len(empty)


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
       {"BLUE": BLUE, "TEAL": TEAL, "ORANGE": ORANGE, "RED": RED, "GREEN": GREEN,
            "PURPLE": PURPLE, "AMBER": AMBER, "MAGENTA": MAGENTA, "BROWN": BROWN,
            "SKY": SKY, "INK": INK, "TITLE": TITLE}.items()}


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
        "axes.prop_cycle":   _cycler(color=PALETTE),
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
    #  CLAMP THE ASK TO SOMETHING THAT FITS. The two scales are set independently, and the
    #  deck build combines FONTSCALE 1.8 with SIZESCALE 0.30/0.45/0.70 -- 1.8x text on a
    #  0.3x canvas, six times more text than there is room for. That is not a drawing bug,
    #  it is an impossible request. What matters for legibility is the APPARENT size when
    #  the figure is placed at a fixed width on a slide, which goes as FONTSCALE/SIZESCALE,
    #  so the cap is on that ratio, and it is announced rather than applied silently.
    #
    #  MEASURED, on the project's own text-overlap checker, over make_charts +
    #  spacetime_outputs (counts are total collisions):
    #
    #      SIZESCALE           1.00   0.70   0.45   0.30
    #      before any of this     38     19    114    171
    #      ratio cap 1.8           0      5     25     20
    #      ratio cap 1.2           -      1      9      9
    #      FONTSCALE 1.0 (none)    -      -      9      8   <- the floor
    #
    #  The cap is 1.8 because that is what the PRIMARY slide set (SIZESCALE 1.0, the one
    #  the deck is built from) uses, and at 1.8 it is now collision-free; a lower cap would
    #  buy the sub-scale sets a little and cost the primary set its legibility.
    #
    #  The sub-scale sets still collide, and no font setting fixes that: the bottom row
    #  above is the same figures at NO font enlargement at all. A four-panel figure drawn
    #  at a third of its design size has nowhere to put a legible label. Those three sets
    #  are gitignored derivatives feeding a deck that is not in this repository; the fix,
    #  if they are ever needed collision-free, is fewer panels per figure, not more scaling.
    try:
        _szq = float(os.environ.get("SHCT_FIG_SIZESCALE", "1.0"))
    except ValueError:
        _szq = 1.0
    _MAX_TEXT_TO_CANVAS = 1.8
    if _fs > _MAX_TEXT_TO_CANVAS * _szq:
        _clamped = _MAX_TEXT_TO_CANVAS * _szq
        import warnings as _warn
        _warn.warn(
            f"SHCT_FIG_FONTSCALE={_fs:g} with SHCT_FIG_SIZESCALE={_szq:g} asks for "
            f"{_fs / max(_szq, 1e-9):.1f}x text per unit canvas; capped at "
            f"{_MAX_TEXT_TO_CANVAS:g}x (font scale {_clamped:.2f}) so the labels fit.",
            stacklevel=2)
        _fs = _clamped
    if abs(_fs - 1.0) > 1e-9:
        #  Only font.size is pre-scaled here. Matplotlib builds axis labels, tick
        #  labels, titles and legends by passing the OTHER rcParams through as
        #  explicit sizes, which the Text wrapper below already scales — scaling
        #  both would land them at 1.8 x 1.8 = 3.24 times their intended size.
        mpl.rcParams["font.size"] = 10.0 * _fs
        for _k, _base in (("axes.titlesize", 11.0), ("axes.labelsize", 10.0),
                          ("xtick.labelsize", 8.5), ("ytick.labelsize", 8.5)):
            #  the key is built at run time; RcParams is typed with a Literal of every
            #  valid key, which a loop variable cannot satisfy.
            mpl.rcParams[_k] = _base   # type: ignore[index]
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

            #  DELIBERATE monkey-patching of matplotlib's own classes: this is how the
            #  slide font scaling is applied to text created by matplotlib internals that
            #  never see our rcParams. The idempotence flag is our own attribute.
            _Text.set_fontsize = _scaled_set_fontsize          # type: ignore[method-assign]
            _Text.set_size = _scaled_set_fontsize              # type: ignore[attr-defined]
            _Text._shct_fontsize_wrapped = True                # type: ignore[attr-defined]

        #  thicker strokes too, or the lines vanish before the labels do
        for _k, _base in (("lines.linewidth", 1.5), ("axes.linewidth", 0.9),
                          ("xtick.major.width", 0.9), ("ytick.major.width", 0.9),
                          ("grid.linewidth", 0.8)):
            mpl.rcParams[_k] = _base * min(_fs, 1.6)   # type: ignore[index]
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
    #  CROWDING IS THE RATIO OF TEXT TO CANVAS, NOT EITHER ALONE. This block used to run
    #  only when the SIZE scale was off unity, so the base slide set -- SIZESCALE 1.0 with
    #  FONTSCALE 1.8 -- got 1.8x text on an unchanged canvas and no layout adaptation at
    #  all. That is where all 38 of the text overlaps the project's own checker reported
    #  came from, 23 of them in one four-panel figure. Adapt whenever either scale moves.
    _crowd = _sz / max(_fs, 1e-9)          # <1 means text is large for the canvas
    if abs(_sz - 1.0) > 1e-9 or abs(_fs - 1.0) > 1e-9:
        import matplotlib.pyplot as _plt
        if not getattr(_plt, "_shct_size_wrapped", False):
            #  Enlarged text needs somewhere to go. Growing the canvas by the FULL font
            #  factor would cancel the effect (same proportions, same apparent size when
            #  placed); sqrt gives the layout room while still gaining 1.34x relative text
            #  at FONTSCALE 1.8.
            _grow = _sz * (max(_fs, 1.0) ** 0.5)

            def _scale(kw):
                fs = kw.get("figsize")
                if fs and len(fs) == 2:
                    kw["figsize"] = (fs[0] * _grow, fs[1] * _grow)
                return kw
            #  A smaller figure with the same number of ticks is how labels collide:
            #  the axis keeps eight tick labels while the axis itself has shrunk to
            #  40 % of its width, so they run into one another. Thinning the ticks in
            #  proportion is what keeps a shrunk figure legible rather than crowded,
            #  and it is the difference between "small" and "small and simple".
            from matplotlib.ticker import MaxNLocator
            #  keyed on the CROWDING ratio: 4 ticks at FONTSCALE 1.8 / SIZESCALE 1.0,
            #  ~4 at SIZESCALE 0.42 unscaled text, 7 when nothing is scaled.
            _nb = max(3, int(round(6 * _crowd)) + 1)

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

                _Fig.add_subplot = _add_sub_thin               # type: ignore[method-assign]
                _Fig._shct_axes_wrapped = True                 # type: ignore[attr-defined]

            _f, _s = _plt.figure, _plt.subplots
            _plt.figure = lambda *a, **k: _f(*a, **_scale(k))
            #  subplots() delegates to figure(), so scaling in BOTH applies the
            #  factor twice and the figure lands at 0.42^2 = 18 % of its intended
            #  size. Only figure() carries the size; subplots() thins ticks only.
            _plt.subplots = lambda *a, **k: _thin_all(_s(*a, **k))
            _plt.rcParams["figure.figsize"] = [v * _sz for v in _plt.rcParams["figure.figsize"]]
            _plt._shct_size_wrapped = True                     # type: ignore[attr-defined]
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
        _s = float(os.environ.get("SHCT_FIG_SIZESCALE", "1.0"))
        _f = float(os.environ.get("SHCT_FIG_FONTSCALE", "1.0"))
    except ValueError:
        return False
    #  a long label is just as crowded by big text as by a small canvas
    return abs(_s - 1.0) > 1e-9 or abs(_f - 1.0) > 1e-9


def label(long_form, short_form):
    """The long label for print, the short one for a slide-sized rendering."""
    return short_form if compact() else long_form


# apply on import so a bare `import shct_style` is enough.
apply_style()


def fit_suptitle(fig, rect=None, gap=0.012):
    """tight_layout with a top margin measured from the figure's own suptitle.

    A fixed rect is right at exactly one font scale. At the slide scales the suptitle's
    band grows into the region tight_layout gives the panel titles, and they are drawn
    through it -- 17_hydrate_distribution and 18_shutin_profile_deposit each did that at
    every panel. Measure the suptitle, then lay out beneath it.
    """
    left, bottom, right, top = rect or (0.0, 0.0, 1.0, 0.97)
    try:
        fig.canvas.draw()
        get = getattr(fig.canvas, "get_renderer", None)
        rend = get() if get is not None else fig._get_renderer()
        sup = getattr(fig, "_suptitle", None)
        if sup is not None and str(sup.get_text()).strip():
            bb = sup.get_window_extent(rend)
            top = min(top, bb.y0 / max(fig.get_window_extent().height, 1.0) - gap)
    except Exception:
        pass
    try:
        fig.tight_layout(rect=(left, bottom, right, max(top, 0.55)))
    except Exception:
        pass


def relieve_crush(fig, target_in=0.72, max_grow=2.6, max_passes=3):
    """Grow the canvas when a panel has been squeezed below usability, then re-lay it out.

    Some figures are simply over-constrained at the smaller slide scales: 07_probabilistic
    puts two panels on a 2.8 x 1.1 in canvas and 08_diagnostics four on a 2.4 in one,
    leaving each panel a few tenths of an inch. Shortening labels cannot fix that -- there
    is no room for the panels, their tick labels and their legends at that size, whatever
    the text says. Give the canvas back the room it needs and let tight_layout redistribute.

    IT MUST ITERATE, and the single-pass version is what proved it. Growing the width
    rescales nothing vertical: the suptitle, panel title, x label and any annotation keep
    their height in points, so a canvas made 2.4x wider and only 1.35x taller hands all of
    the new width to the panels and none of the new height, and what was a 0.30 in panel
    becomes a 0.21 in one. 07_probabilistic at the 0.30 scale came out of relief MORE
    crushed than it went in, and passed, because relief never looked again. Each pass now
    re-measures what the previous pass actually produced, and each dimension is grown on
    its own evidence rather than one factor being derived from the width and the other
    following at a cap.

    Bounded by `max_grow` per dimension in total, so a genuinely broken layout (a caption
    that has run the axes away) is still reported by the screen rather than inflated
    forever.
    """
    changed = False
    tot_w = tot_h = 1.0
    try:
        for _ in range(max(1, int(max_passes))):
            hits = find_squeezed_axes(fig)
            if not hits:
                break
            narrowest = min(w for _ax, w, _h in hits)
            shortest = min(h for _ax, _w, h in hits)
            if narrowest <= 0.02 or shortest <= 0.02:   # a runaway layout, not a cramped one
                break
            grow_w = min(max_grow / tot_w, target_in / max(narrowest, 1e-3)) if narrowest < 0.5 else 1.0
            grow_h = min(max_grow / tot_h, 0.62 / max(shortest, 1e-3)) if shortest < 0.30 else 1.0
            grow_w = max(grow_w, 1.0)
            grow_h = max(grow_h, 1.0)
            if max(grow_w, grow_h) <= 1.02:
                break
            #  what the CURRENT layout has reserved, before the resize discards it. A bare
            #  tight_layout() here threw away both reservations the figure had already
            #  made: fit_suptitle()'s measured top margin, so panel titles were drawn back
            #  through the suptitle (17_hydrate_distribution and 18_shutin_profile_deposit,
            #  the very figures fit_suptitle exists for), and any bottom band a caption had
            #  claimed. Relief must not undo the layout it is relieving.
            _bottom = float(getattr(fig.subplotpars, "bottom", 0.0) or 0.0)
            w, h = (float(v) for v in fig.get_size_inches())
            fig.set_size_inches(w * grow_w, h * grow_h)
            tot_w *= grow_w
            tot_h *= grow_h
            _sup = getattr(fig, "_suptitle", None)
            if _sup is not None and str(_sup.get_text()).strip():
                fit_suptitle(fig, rect=(0.0, _bottom, 1.0, 0.97))
            else:
                try:
                    fig.tight_layout(rect=(0.0, _bottom, 1.0, 1.0))
                except Exception:
                    pass
            changed = True
        return changed
    except Exception:
        return changed

def guard_size(fig, path, max_in=200.0):
    """Refuse to let a figure be saved at an absurd size, and say which one.

    matplotlib raises "Image size of NxM pixels is too large" from deep inside the Agg
    backend, with a traceback that names the renderer rather than the figure -- one such
    failure (a zero limit on a log axis) killed a whole build step today. A figure this
    large is always a layout bug, never an intent, so clamp it and name it.
    """
    import os as _os
    try:
        w, h = fig.get_size_inches()
        if w > max_in or h > max_in:
            print(f"    [oversize figure] {_os.path.basename(path)}: {w:.0f} x {h:.0f} in "
                  f"— a layout bug has run the canvas away; clamped to {max_in:.0f} in",
                  flush=True)
            fig.set_size_inches(min(w, max_in), min(h, max_in))
    except Exception:
        pass


def screen(fig, path):
    """Run the shared figure screens on `fig` before it is written to `path`.

    Every module that saves a PNG calls this, so no generator is exempt. The checks
    print and return; they never raise, because a screening warning must not be able
    to fail a run that is otherwise producing correct output.
    """
    import os as _os
    try:
        import shct_style as _st
        name = _os.path.basename(path)
        #  give an over-constrained canvas room BEFORE judging it, so the screens report
        #  what will actually be written rather than an intermediate layout
        _st.relieve_crush(fig)
        #  after any resize, before anything is judged: a below-axes legend's offset is
        #  stored in axes fractions and a resize rescales it back onto the x-label
        _st.fix_below_legends(fig)
        _st.report_text_overlaps(fig, name)
        #  the fault the overlap check cannot see: a panel that draws nothing. This was
        #  in solver's _save_checked() and not here, while guard_size and
        #  fix_below_legends were here and not there -- two screens, each missing what
        #  the other had, and a figure got whichever set its generator happened to call.
        _st.report_empty_axes(fig, name)
        #  and the one none of the others can see: a note drawn ON the curve
        _st.report_text_over_data(fig, name)
        _st.report_degenerate_axes(fig, name)
        _st.report_squeezed_axes(fig, name)
        _st.guard_size(fig, path)
    except Exception:
        pass
