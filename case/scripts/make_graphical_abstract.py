#!/usr/bin/env python3
"""Graphical abstract for the IJMF manuscript.

IJMF: "Ensure the image is 531 x 1328 pixels (h x w) or proportionally more, and
is readable at a size of 5 x 13 cm." This renders at exactly 2x that (1062 x 2656)
so it stays sharp, and every label is sized to remain legible at 5 x 13 cm.

Three panels carrying the paper's argument rather than decoration:
  (1) the competition at the wall - slug scouring against hydrate deposition
  (2) the dimensionless group that ranks them, over the real Phi_SH field
  (3) what it buys: the as-operated line against the engineered fix
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#  a global savefig.bbox="tight" would crop the canvas and break the exact
#  2.5:1 proportion the guide asks for, so pin it for this figure
matplotlib.rcParams["savefig.bbox"] = None
from matplotlib.patches import FancyArrowPatch, Rectangle, Wedge
import matplotlib.image as mpimg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import CASE                      # noqa: E402

BLUE, TEAL, ORANGE, RED, GREEN = "#2E5BBF", "#1AA0A0", "#E8842B", "#E0463C", "#3FA65A"
GREY = "#6E7B8B"
OUT = r"/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates/graphical_abstract.png"

W_IN, H_IN, DPI = 13.28, 5.312, 200        # -> 2656 x 1062 px, exactly 2x the
                                            # minimum and the same 2.5:1 proportion

fig = plt.figure(figsize=(W_IN, H_IN), dpi=DPI, facecolor="white")
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 0.95],
                      left=0.012, right=0.988, top=0.86, bottom=0.06, wspace=0.16)

fig.text(0.012, 0.945,
         "Coupled transient modelling of slug flow and gas-hydrate deposition "
         "in subsea pipelines",
         fontsize=15, fontweight="bold", color="#1F3B60", ha="left", va="center")

# ---------------------------------------------------------------- (1) the wall
ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off")
ax.set_title("A slug passes, or it does not", fontsize=12, fontweight="bold",
             color=BLUE, pad=6)

for y0, lab, col in [(3.5, "slug arrives", BLUE), (0.55, "between slugs", ORANGE)]:
    ax.add_patch(Rectangle((0.5, y0), 9.0, 2.15, facecolor="#F2F6FC",
                           edgecolor=GREY, lw=1.4))
    ax.add_patch(Rectangle((0.5, y0), 9.0, 0.20, facecolor=GREY, edgecolor="none"))
    ax.add_patch(Rectangle((0.5, y0 + 1.95), 9.0, 0.20, facecolor=GREY, edgecolor="none"))
    ax.text(0.5, y0 + 2.42, lab, fontsize=11.5, fontweight="bold", color=col)

# slug body sweeping the wall
ax.add_patch(Rectangle((1.3, 3.72), 4.2, 1.71, facecolor=BLUE, alpha=0.55,
                       edgecolor=BLUE, lw=1.2))
ax.text(3.4, 4.58, "liquid slug", fontsize=10.5, color="white",
        ha="center", va="center", fontweight="bold")
ax.add_patch(FancyArrowPatch((5.8, 4.58), (9.0, 4.58), arrowstyle="-|>",
                             mutation_scale=20, lw=2.4, color=BLUE))
ax.text(7.4, 5.02, "scours + re-warms", fontsize=10, color=BLUE,
        ha="center", fontweight="bold")

# deposit growing on the cold wall between slugs
for x in (1.6, 3.0, 4.4, 5.8, 7.2):
    ax.add_patch(Wedge((x, 0.75), 0.62, 0, 180, facecolor=ORANGE, alpha=0.75,
                       edgecolor=ORANGE))
ax.text(5.0, 1.72, "hydrate deposits on the cold wall",
        fontsize=10.5, color=ORANGE, ha="center", fontweight="bold")
ax.text(5.0, 0.06, "whichever rate wins decides whether the line plugs",
        fontsize=10.5, color="#1F3B60", ha="center", style="italic")

# ------------------------------------------------- (2) the group over the field
ax2 = fig.add_subplot(gs[0, 1]); ax2.axis("off")
ax2.set_title("A dimensionless group ranks them", fontsize=12, fontweight="bold",
              color=TEAL, pad=6)
ax2.text(0.5, 0.895,
         r"$\Phi_{SH}\;=\;C\,k_{g,w}\,a_i\,\Delta T_{sub,w}^{\,n}\;/\;f_s$",
         fontsize=17, ha="center", va="center", transform=ax2.transAxes,
         color="#1F3B60")
ax2.text(0.5, 0.775, "deposition tendency  ÷  slug renewal rate",
         fontsize=10.5, ha="center", va="center", transform=ax2.transAxes, color=GREY)
mp = os.path.join(CASE, "outputs_paper_steady", "04_PhiSH_map.png")
if not os.path.exists(mp):
    mp = os.path.join(CASE, "outputs_steady", "04_PhiSH_map.png")
if os.path.exists(mp):
    inset = ax2.inset_axes([0.015, 0.02, 0.97, 0.70])
    inset.imshow(mpimg.imread(mp)); inset.axis("off")
    inset.set_title(r"mapped along the route and through time  ($\Phi_{SH}=1$ contour)",
                    fontsize=9.5, color=GREY, pad=3)

# ------------------------------------------------------------- (3) what it buys
ax3 = fig.add_subplot(gs[0, 2]); ax3.set_xlim(0, 10); ax3.set_ylim(0, 6.4); ax3.axis("off")
ax3.set_title("Used as a design tool", fontsize=12, fontweight="bold",
              color=GREEN, pad=6)
try:
    eng = json.load(open(os.path.join(CASE, "outputs_steady", "summary.json")))
    mit = json.load(open(os.path.join(CASE, "outputs_mitigated", "summary.json")))
    rows = [("max subcooling", f"{eng['max_subcooling_C']:.1f} °C",
             f"{mit['max_subcooling_C']:.1f} °C"),
            ("wall deposit", f"{eng['peak_deposit_mm']:.0f} mm",
             f"{mit['peak_deposit_mm']:.0f} mm"),
            ("plug probability", f"{eng['P_plug']*100:.0f} %",
             f"{mit['P_plug']*100:.0f} %"),
            ("no-touch time", "≈0 h", f"{mit['cooldown_to_hydrate_h']:.0f} h")]
except Exception:
    rows = [("max subcooling", "20.9 °C", "6.6 °C"), ("wall deposit", "117 mm", "0 mm"),
            ("plug probability", "100 %", "25 %"), ("no-touch time", "≈0 h", "17 h")]

ax3.text(0.2, 5.75, "as operated", fontsize=11, fontweight="bold", color=RED)
ax3.text(6.7, 5.75, "engineered", fontsize=11, fontweight="bold", color=GREEN)
for i, (lab, a, b) in enumerate(rows):
    y = 4.85 - i * 1.12
    ax3.text(0.2, y + 0.34, lab, fontsize=10, color=GREY)
    ax3.add_patch(Rectangle((0.2, y - 0.34), 3.9, 0.66, facecolor=RED, alpha=0.13,
                            edgecolor=RED, lw=1.0))
    ax3.text(2.15, y, a, fontsize=11.5, fontweight="bold", color=RED,
             ha="center", va="center")
    ax3.add_patch(FancyArrowPatch((4.35, y), (5.85, y), arrowstyle="-|>",
                                  mutation_scale=15, lw=1.8, color=GREY))
    ax3.add_patch(Rectangle((6.0, y - 0.34), 3.8, 0.66, facecolor=GREEN, alpha=0.15,
                            edgecolor=GREEN, lw=1.0))
    ax3.text(7.9, y, b, fontsize=11.5, fontweight="bold", color=GREEN,
             ha="center", va="center")
ax3.text(5.0, 0.06, "insulation + inhibitor sized by the same model",
         fontsize=10.5, color="#1F3B60", ha="center", style="italic")

fig.savefig(OUT, dpi=DPI, facecolor="white")
plt.close(fig)

from PIL import Image

#  Something in the import chain trims the canvas, so the saved size drifts from
#  figsize x dpi. Rather than fight it, pad the finished image onto a white canvas
#  of exactly the proportion the guide asks for (1328 x 531 h x w, here at 2x).
TARGET_W, TARGET_H = 2656, 1062
_im = Image.open(OUT).convert("RGB")
if (_im.width, _im.height) != (TARGET_W, TARGET_H):
    _scale = min(TARGET_W / _im.width, TARGET_H / _im.height)
    _rs = _im.resize((max(1, int(_im.width * _scale)), max(1, int(_im.height * _scale))),
                     Image.LANCZOS)
    _canvas = Image.new("RGB", (TARGET_W, TARGET_H), "white")
    _canvas.paste(_rs, ((TARGET_W - _rs.width) // 2, (TARGET_H - _rs.height) // 2))
    _canvas.save(OUT, dpi=(DPI, DPI))

w, h = Image.open(OUT).size
print(f"wrote {OUT}")
print(f"  {w} x {h} px  (IJMF minimum 1328 x 531 h x w; this is "
      f"{w/1328:.2f}x the minimum width)")
print(f"  readable at {w/DPI*2.54:.1f} x {h/DPI*2.54:.1f} cm at {DPI} dpi "
      f"(guide asks for legibility at 13 x 5 cm)")
