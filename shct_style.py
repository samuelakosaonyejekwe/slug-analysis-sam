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
#  darkest stop used anywhere is the medium royal blue / medium magenta anchors.
# STRONG, high-contrast, multi-hue map (selected): blue -> cyan -> green -> amber ->
# orange -> red. Saturated and bold so values that cluster in the mid-range get
# DISTINCT strong colours (the soft single-hue maps washed those out); every stop is
# still medium/saturated — none is black, near-black or dark. Used with smooth
# (gouraud) shading for EVERY colour-mapped field: holdup, Φ_SH, deposit, temperature,
# velocity. Low = strong blue, high = strong red (intuitive: red = hot / critical / full).
_STRONG = ["#2348FF", "#00B6F5", "#00D63C", "#FFE400", "#FF8A00", "#FF1A1A"]   # MAX-vivid: pure saturated blue/cyan/green/yellow/orange/red, no dark
CMAP_SEQ  = LinearSegmentedColormap.from_list("shct_seq",  _STRONG)   # holdup, velocity
CMAP_HEAT = LinearSegmentedColormap.from_list("shct_heat", _STRONG)   # deposit thickness
CMAP_TEMP = LinearSegmentedColormap.from_list("shct_temp", _STRONG)   # temperature
CMAP_DIV  = LinearSegmentedColormap.from_list("shct_div",  _STRONG)   # Φ_SH (critical contour drawn on top)

for _cm in (CMAP_SEQ, CMAP_HEAT, CMAP_TEMP, CMAP_DIV):
    try:
        mpl.colormaps.register(_cm, force=True)       # register by name (shct_seq, ...)
    except Exception:                                 # pragma: no cover (old mpl)
        try:
            mpl.cm.register_cmap(name=_cm.name, cmap=_cm)
        except Exception:
            pass

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
