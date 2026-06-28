#!/usr/bin/env python3
# =============================================================================
#  shct_crosssection.py — REDUCED-ORDER CROSS-SECTION (quasi-3-D) RECONSTRUCTION
# -----------------------------------------------------------------------------
#  The SHCT core solves 1-D area-averaged transport along the pipe axis x. This
#  module RECONSTRUCTS the cross-sectional (y-z plane) structure of the flow from
#  that 1-D solution — per axial cell — to resolve effects an area-averaged model
#  alone cannot show:
#     * the stratified gas/liquid interface level h/D (from the holdup, via the
#       exact circular-segment geometry — Taitel-Dukler),
#     * the wetted-perimeter fraction and interface width,
#     * a vertical/radial velocity profile across the section (turbulent 1/7 law),
#     * a radial temperature profile (cold wall -> warm core),
#     * the AZIMUTHAL distribution of the hydrate wall deposit (thicker at the
#       cold, liquid-/water-wetted bottom of line — the real bottom-of-line
#       corrosion/hydrate location).
#  Combined with the axial coordinate this is a (y, z, x) = quasi-3-D field.
#
#  HONEST SCOPE (read before citing): this is a fast REDUCED-ORDER reconstruction
#  consistent with the 1-D conservation laws — NOT a 3-D Navier-Stokes solve. It
#  does NOT resolve turbulence, secondary flows or eddies the way a full 3-D CFD
#  code (ANSYS Fluent, OpenFOAM) does, and is NOT a substitute for or "better
#  than" CFD. It is COMPLEMENTARY: seconds vs. hours, whole-line vs. a single
#  short section, and it gives the cross-sectional picture flow-assurance design
#  needs (interface, bottom-of-line deposit, velocity skew) at screening cost.
#  Pure post-processing: it reads the solved fields and changes nothing in the core.
# =============================================================================
from __future__ import annotations
import os, math
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:                                       # pragma: no cover
    _HAVE_MPL = False

NAVY = "#2E5BBF"; ACCENT = "#1F8AC0"; RED = "#E0463C"; ORANGE = "#E8842B"


# ---------------------------------------------------------------------------
#  Circular-segment geometry: invert the holdup -> liquid level h/D
# ---------------------------------------------------------------------------
def _area_fraction(h):
    """Liquid area fraction of a horizontal circular pipe filled to height fraction
    h = h_liquid / D, h in [0,1].  A_L/A = (acos(1-2h) - (1-2h)*sqrt(1-(1-2h)^2))/pi."""
    s = 1.0 - 2.0 * np.clip(h, 0.0, 1.0)
    return (np.arccos(s) - s * np.sqrt(np.maximum(1.0 - s * s, 0.0))) / math.pi


def liquid_level(alpha_l):
    """Invert A_L/A = alpha_l for the dimensionless liquid level h/D (vectorised bisection)."""
    alpha_l = np.clip(np.asarray(alpha_l, float), 1e-4, 1.0 - 1e-4)
    lo = np.zeros_like(alpha_l); hi = np.ones_like(alpha_l)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = _area_fraction(mid) - alpha_l
        hi = np.where(f > 0, mid, hi)
        lo = np.where(f > 0, lo, mid)
    return 0.5 * (lo + hi)


def section_geometry(alpha_l, D):
    """Return per-cell stratified geometry from holdup & diameter:
    h_over_D, wetted-perimeter fraction (liquid), interface width (m)."""
    h = liquid_level(alpha_l)
    s = 1.0 - 2.0 * h
    gamma = np.arccos(s)                       # half-angle subtended by the wetted wall
    wetted_frac = gamma / math.pi             # S_liquid / (pi D)
    interface_w = D * np.sqrt(np.maximum(1.0 - s * s, 0.0))
    return h, wetted_frac, interface_w


# ---------------------------------------------------------------------------
#  Azimuthal deposit distribution (bottom-of-line weighting)
# ---------------------------------------------------------------------------
def azimuthal_deposit(delta_mean, h_over_D, skew=1.6):
    """Distribute the area-mean deposit thickness delta around the circumference,
    weighted toward the cold liquid-wetted BOTTOM of the line. theta_az = 0 at the
    bottom, pi at the top. Returns (theta_az[m], delta_profile[m,ncell]) and the
    bottom/top thicknesses. Conserves the azimuthal mean = delta_mean."""
    theta = np.linspace(0.0, math.pi, 37)                 # 0=bottom .. pi=top
    # weight: high at bottom (liquid + water settling + coldest), decays to top;
    # the liquid covers up to angle gamma from the bottom -> stronger weight there.
    w = (1.0 + skew * np.cos(theta)) / 2.0                # 0..(1+skew)/2, cos: 1 at bottom
    w = np.clip(w, 0.05, None)
    w = w / w.mean()                                      # normalise so azimuthal mean = 1
    prof = np.outer(w, np.asarray(delta_mean, float))     # (ntheta, ncell)
    return theta, prof, prof[0], prof[-1]


# ---------------------------------------------------------------------------
#  2-D cross-section field reconstruction at one axial station
# ---------------------------------------------------------------------------
def reconstruct_section(D, alpha_l, u_mix, T_bulk, T_wall, delta_mean, n=120,
                        velocity_exp=1.0 / 7.0, gas_factor=1.25, liq_factor=0.85, skew=1.6):
    """Build (n x n) cross-section fields at one station: phase (1=liquid,0=gas),
    velocity (turbulent profile, liquid slower than gas), temperature (cold wall ->
    warm core), and a wall-deposit mask (azimuthally bottom-weighted). The closure
    constants (velocity_exp, gas/liq factors, azimuthal skew) are calibratable (item 9)."""
    R = 0.5 * D
    yy = np.linspace(-R, R, n); zz = np.linspace(-R, R, n)
    Z, Y = np.meshgrid(zz, yy)
    r = np.sqrt(Y * Y + Z * Z)
    inside = r <= R
    h = float(liquid_level(np.array([alpha_l]))[0])
    y_int = R * (2.0 * h - 1.0)                           # interface height (y)
    liquid = inside & (Y <= y_int)
    # velocity: turbulent power-law on distance from wall; gas faster than liquid
    dist_wall = np.maximum(R - r, 0.0)
    prof = (dist_wall / R) ** velocity_exp
    vel = np.where(liquid, liq_factor * u_mix, gas_factor * u_mix) * prof
    vel = np.where(inside, vel, np.nan)
    # temperature: cold wall -> warm centre (radial), reduced-order
    temp = T_wall + (T_bulk - T_wall) * (1.0 - (r / R) ** 2)
    temp = np.where(inside, temp, np.nan)
    # azimuthal deposit ring at the wall (bottom-weighted)
    theta_pt = np.arctan2(-(Y), Z)                        # 0 at +z; we want 0 at bottom (-y)
    theta_bottom = np.arctan2(R, 0.0)                     # unused; compute bottom-referenced angle
    ang = np.arctan2(Z, -Y)                               # 0 at bottom (-y), +-pi at top
    w_az = (1.0 + skew * np.cos(ang)) / 2.0
    w_az = np.clip(w_az / np.mean((1.0 + skew * np.cos(np.linspace(-math.pi, math.pi, 200))) / 2.0),
                   0.05, None)
    local_delta = delta_mean * w_az
    deposit = inside & (r >= (R - local_delta))
    return dict(Y=Y, Z=Z, inside=inside, liquid=liquid, vel=vel, temp=temp,
                deposit=deposit, y_int=y_int, R=R, h=h)


# ---------------------------------------------------------------------------
#  Top-level: write the cross-section CSV + charts from a solved run
# ---------------------------------------------------------------------------
def crosssection_outputs(sv, outdir, stations_km=None):
    """Post-process a solved TransientSHCT `sv`: write csv_crosssection.csv and the
    cross-section charts (longitudinal geometry, azimuthal deposit map, and 2-D
    section reconstructions at representative stations) into `outdir`."""
    os.makedirs(outdir, exist_ok=True)
    r = sv.results
    med = lambda A: np.nanmedian(A, 1)
    x_km = sv.x / 1000.0
    D = med(r["D"]) if "D" in r else np.full_like(x_km, sv.case.pipeline.diameter_m)
    alpha_l = med(r["alpha_l"])
    u_mix = med(r["j"])
    T_bulk = med(r["T"])
    delta = med(r["delta"])
    T_wall = float(sv.case.operating.T_seabed_C)
    #  item 9: reduced-order closure constants are read from the case (calibratable), with the
    #  original literature defaults preserved when unset.
    nm = sv.case.numerics
    skew = float(getattr(nm, "cx_deposit_skew", 1.6))
    vexp = float(getattr(nm, "cx_velocity_exp", 1.0 / 7.0))
    gfac = float(getattr(nm, "cx_gas_vel_factor", 1.25))
    lfac = float(getattr(nm, "cx_liq_vel_factor", 0.85))

    h, wetted_frac, interface_w = section_geometry(alpha_l, D)
    theta, depo_prof, depo_bot, depo_top = azimuthal_deposit(delta, h, skew=skew)

    # --- CSV ---
    cols = ["x_km", "diameter_mm", "holdup_alpha_l", "liquid_level_h_over_D", "liquid_level_mm",
            "wetted_perim_frac", "interface_width_mm", "deposit_mean_mm",
            "deposit_bottom_mm", "deposit_top_mm", "u_mix_mps"]
    rows = np.column_stack([x_km, D * 1000.0, alpha_l, h, h * D * 1000.0, wetted_frac,
                            interface_w * 1000.0, delta * 1000.0, depo_bot * 1000.0,
                            depo_top * 1000.0, u_mix])
    with open(os.path.join(outdir, "csv_crosssection.csv"), "w") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(f"{v:.5g}" for v in row) + "\n")

    if not _HAVE_MPL:
        return os.path.join(outdir, "csv_crosssection.csv")

    # --- chart 1: longitudinal cross-section geometry curves ---
    fig, ax = plt.subplots(3, 1, figsize=(8, 6.2), sharex=True)
    ax[0].fill_between(x_km, sv.z, sv.z.min() - 20, color="#C9B79B", alpha=.55)
    ax[0].plot(x_km, sv.z, color="#B07A33"); ax[0].set_ylabel("elev (m)")
    ax[0].set_title("Cross-section reconstruction along the line (quasi-3-D)",
                    color=NAVY, fontweight="bold")
    ax[1].plot(x_km, h, color=ACCENT, lw=1.8, label="liquid level h/D")
    ax[1].plot(x_km, wetted_frac, color=ORANGE, lw=1.4, label="wetted-perimeter fraction")
    ax[1].set_ylabel("h/D , wetted frac"); ax[1].set_ylim(0, 1)
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    ax[2].plot(x_km, depo_bot * 1000.0, color=RED, lw=1.8, label="bottom-of-line deposit (mm)")
    ax[2].plot(x_km, depo_top * 1000.0, color=NAVY, lw=1.2, ls="--", label="top-of-line deposit (mm)")
    ax[2].set_ylabel("deposit (mm)"); ax[2].set_xlabel("distance (km)")
    ax[2].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "cx1_geometry.png"), dpi=150); plt.close(fig)

    # --- chart 2: azimuthal deposit "unrolled" map (x vs azimuth) ---
    fig, axm = plt.subplots(figsize=(7.6, 4.2))
    pcm = axm.pcolormesh(x_km, np.degrees(theta), depo_prof * 1000.0, cmap="shct_heat", shading="auto")
    axm.set_xlabel("distance (km)"); axm.set_ylabel("azimuth (deg: 0=bottom, 180=top)")
    axm.set_title("Azimuthal hydrate-deposit distribution δ(x, θ) — bottom-of-line accumulation",
                  color=NAVY, fontweight="bold")
    fig.colorbar(pcm, ax=axm, label="deposit thickness (mm)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "cx2_azimuthal_deposit.png"), dpi=150); plt.close(fig)

    # --- chart 3: 2-D section reconstructions at representative stations ---
    if stations_km is None:
        # inlet, mid-line, and the coldest/most-critical cell (max deposit)
        i_hot = int(np.argmax(delta)) if np.nanmax(delta) > 0 else int(0.5 * len(x_km))
        idxs = sorted(set([2, len(x_km) // 2, i_hot, len(x_km) - 3]))
    else:
        idxs = [int(np.argmin(np.abs(x_km - s))) for s in stations_km]
    idxs = idxs[:4]
    fig, axes = plt.subplots(1, len(idxs), figsize=(3.4 * len(idxs), 3.6))
    axes = np.atleast_1d(axes)
    for ax, i in zip(axes, idxs):
        sec = reconstruct_section(float(D[i]), float(alpha_l[i]), float(max(u_mix[i], 0.05)),
                                  float(T_bulk[i]), T_wall, float(delta[i]),
                                  velocity_exp=vexp, gas_factor=gfac, liq_factor=lfac, skew=skew)
        v = np.ma.masked_invalid(sec["vel"])
        pcm = ax.pcolormesh(sec["Z"], sec["Y"], v, cmap="shct_seq", shading="auto")
        # phase interface line
        ax.axhline(sec["y_int"], color="white", lw=1.2, ls="--")
        # deposit ring
        ax.contourf(sec["Z"], sec["Y"], sec["deposit"].astype(float), levels=[0.5, 1.5], colors=[RED], alpha=.8)
        ax.add_patch(plt.Circle((0, 0), sec["R"], fill=False, color="#3A5BA8", lw=1.0))
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"x={x_km[i]:.1f} km\nα_l={alpha_l[i]:.2f}, δ={delta[i]*1000:.0f}mm",
                     fontsize=8.5, color=NAVY)
    fig.suptitle("2-D cross-section reconstruction — velocity field, gas/liquid interface "
                 "(dashed), wall deposit (red)", color=NAVY, fontweight="bold", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(outdir, "cx3_sections.png"), dpi=150); plt.close(fig)

    return os.path.join(outdir, "csv_crosssection.csv")
