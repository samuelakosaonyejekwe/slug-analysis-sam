#!/usr/bin/env python3
# =============================================================================
#  shct_threed.py — 3-D FIELD RECONSTRUCTION, EXPORT & VISUALISATION
# -----------------------------------------------------------------------------
#  Assembles a full THREE-DIMENSIONAL field over the pipe volume from the solved
#  1-D axial solution and the reduced-order cross-section model (shct_crosssection):
#  for every axial cell it builds the (azimuth a, radial r) cross-section and maps
#  it to world coordinates (X axial, Y/Z cross-section + terrain elevation), giving
#  3-D fields of velocity, temperature, phase (gas/liquid), holdup, subcooling,
#  Phi_SH and wall-hydrate deposit. It then:
#     * exports a STRUCTURED-GRID VTK file (pipe_3d.vtk) openable in ParaView for
#       true interactive 3-D inspection / slicing / iso-surfaces, and
#     * renders 3-D tube views (deposit, temperature) as PNGs.
#
#  HONEST SCOPE (cite this verbatim): the 3-D fields are RECONSTRUCTED from the
#  validated 1-D conservation laws (mass/momentum/energy/hydrate) plus the
#  cross-section closures (Taitel-Dukler interface geometry, turbulent velocity
#  profile, radial thermal profile, bottom-of-line azimuthal deposit). The
#  GOVERNING PHYSICS remains 1-D area-averaged: this is a reduced-order / quasi-3-D
#  reconstruction, NOT a 3-D Navier-Stokes (CFD) solve on a 3-D mesh. It does not
#  resolve turbulence, secondary flow or eddies the way ANSYS Fluent / OpenFOAM do,
#  is not a substitute for or "better than" CFD, and is intended as a fast,
#  whole-line, screening-grade 3-D picture (seconds, 22 km) complementary to CFD
#  (hours, a short section). Pure post-processing; the core solver is unchanged.
# =============================================================================
from __future__ import annotations
import os, math
import numpy as np
import shct_crosssection as cx

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    _HAVE_MPL = True
except Exception:                                       # pragma: no cover
    _HAVE_MPL = False

NAVY = "#2E5BBF"


# ---------------------------------------------------------------------------
#  Assemble the 3-D field on a cylindrical (axial, azimuth, radial) grid
# ---------------------------------------------------------------------------
def build_3d_field(sv, n_axial=60, n_theta=24, n_r=6):
    """Return a structured 3-D field. Coordinates: X = axial distance (m); the
    cross-section (azimuth a measured from the BOTTOM of line, radial fraction rf)
    is mapped to (Yc, Zc) and offset by the pipe-centreline elevation. Arrays are
    shaped (n_axial, n_theta, n_r). Returns a dict of coordinate & field arrays."""
    r = sv.results
    med = lambda A: np.nanmedian(A, 1)
    xall = sv.x
    D = med(r["D"]) if "D" in r else np.full_like(xall, sv.case.pipeline.diameter_m)
    alpha = med(r["alpha_l"]); umix = med(r["j"]); T = med(r["T"])
    delta = med(r["delta"]); sub = med(r["Tsub"]); phish = med(r["max_PhiSH"])
    elev = sv.z
    Twall = float(sv.case.operating.T_seabed_C)

    ia = np.unique(np.linspace(0, len(xall) - 1, min(n_axial, len(xall))).astype(int))
    a = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)        # 0 = bottom of line
    rf = np.linspace(0.0, 1.0, n_r)                                     # 0 = centre, 1 = wall

    shape = (len(ia), n_theta, n_r)
    X = np.zeros(shape); Y = np.zeros(shape); Z = np.zeros(shape)
    VEL = np.zeros(shape); TEMP = np.zeros(shape); PHASE = np.zeros(shape)
    DEP = np.zeros(shape); HOLD = np.zeros(shape); SUB = np.zeros(shape); PHI = np.zeros(shape)

    # azimuthal deposit weight (bottom-of-line), normalised to mean 1
    wnorm = np.mean((1.0 + 1.6 * np.cos(np.linspace(0, 2 * math.pi, 400, endpoint=False))) / 2.0)

    for ki, i in enumerate(ia):
        R = 0.5 * float(D[i])
        h = float(cx.liquid_level(np.array([alpha[i]]))[0])
        y_int = R * (2.0 * h - 1.0)
        for j, aj in enumerate(a):
            w_az = max((1.0 + 1.6 * math.cos(aj)) / 2.0 / wnorm, 0.05)
            local_delta = float(delta[i]) * w_az
            for kk, rr in enumerate(rf):
                yc = -R * rr * math.cos(aj)             # bottom (a=0) -> negative y
                zc = R * rr * math.sin(aj)
                X[ki, j, kk] = xall[i]
                Y[ki, j, kk] = elev[i] + yc
                Z[ki, j, kk] = zc
                liquid = yc <= y_int
                wall_dist_frac = max(1.0 - rr, 1e-3)
                VEL[ki, j, kk] = (0.85 if liquid else 1.25) * max(umix[i], 0.0) * wall_dist_frac ** (1.0 / 7.0)
                TEMP[ki, j, kk] = Twall + (T[i] - Twall) * (1.0 - rr ** 2)
                PHASE[ki, j, kk] = 1.0 if liquid else 0.0
                DEP[ki, j, kk] = 1.0 if (1.0 - rr) * R < local_delta else 0.0
                HOLD[ki, j, kk] = alpha[i]
                SUB[ki, j, kk] = sub[i]
                PHI[ki, j, kk] = phish[i]
    return dict(dims=(n_r, n_theta, len(ia)), X=X, Y=Y, Z=Z,
                fields={"velocity_mps": VEL, "temperature_C": TEMP, "phase_liquid": PHASE,
                        "deposit": DEP, "holdup": HOLD, "subcooling_C": SUB, "Phi_SH": PHI},
                ia=ia, a=a, rf=rf)


# ---------------------------------------------------------------------------
#  VTK STRUCTURED_GRID export (ASCII; opens in ParaView / VisIt)
# ---------------------------------------------------------------------------
def write_vtk(field, path):
    """Write an ASCII VTK STRUCTURED_GRID. Point ordering: r fastest, then theta,
    then axial (i.e. dims = (n_r, n_theta, n_axial))."""
    n_r, n_theta, n_ax = field["dims"]
    X, Y, Z = field["X"], field["Y"], field["Z"]
    npts = n_r * n_theta * n_ax
    # build the flat point list in (axial outer, theta, r inner) order so r varies fastest
    lines = ["# vtk DataFile Version 3.0", "SHCT quasi-3-D reconstructed field",
             "ASCII", "DATASET STRUCTURED_GRID",
             f"DIMENSIONS {n_r} {n_theta} {n_ax}", f"POINTS {npts} float"]
    for k in range(n_ax):
        for j in range(n_theta):
            for i in range(n_r):
                lines.append(f"{X[k, j, i]:.4f} {Y[k, j, i]:.5f} {Z[k, j, i]:.5f}")
    lines.append(f"POINT_DATA {npts}")
    for name, arr in field["fields"].items():
        lines.append(f"SCALARS {name} float 1")
        lines.append("LOOKUP_TABLE default")
        for k in range(n_ax):
            for j in range(n_theta):
                for i in range(n_r):
                    lines.append(f"{arr[k, j, i]:.5g}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
#  3-D tube renders
# ---------------------------------------------------------------------------
def _tube_surface(sv, wall_value, title, cbar_label, cmap, out, r_vis=18.0):
    """Render the pipe as a 3-D tube (axial X following terrain elevation as the
    vertical axis), the wall coloured by `wall_value[a, x]`. The tube radius is
    exaggerated (r_vis, in metres) for visibility against the km-scale length."""
    x_km = sv.x / 1000.0
    elev = sv.z
    n_ax = len(x_km)
    a = np.linspace(0.0, 2.0 * math.pi, wall_value.shape[0])
    U, A = np.meshgrid(x_km, a)                       # (n_a, n_ax)
    Xs = U
    Zs = elev[None, :] + r_vis * np.cos(A)            # vertical (elevation + tube)
    Ys = r_vis * np.sin(A)                            # transverse
    norm = plt.Normalize(np.nanmin(wall_value), max(np.nanmax(wall_value), np.nanmin(wall_value) + 1e-9))
    colors = cm.get_cmap(cmap)(norm(wall_value))
    fig = plt.figure(figsize=(11, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(Xs, Ys, Zs, facecolors=colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False)
    from matplotlib.ticker import MaxNLocator
    ax.set_xlabel("axial distance (km)", labelpad=12)
    ax.set_ylabel("transverse (m, exaggerated)", labelpad=16)
    ax.set_zlabel("elevation (m)", labelpad=10)
    ax.set_title(title, color=NAVY, fontweight="bold")
    # keep the 3-D axes uncluttered: few, well-spaced ticks so labels never jam
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(3))     # transverse axis is the worst offender
    ax.zaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(labelsize=7, pad=1.5)
    m = cm.ScalarMappable(norm=norm, cmap=cmap); m.set_array(wall_value)
    fig.colorbar(m, ax=ax, shrink=0.6, pad=0.10, label=cbar_label)
    ax.view_init(elev=22, azim=-60)
    try:
        ax.set_box_aspect((4, 1, 1.4))
    except Exception:
        pass
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


def threed_outputs(sv, outdir, n_axial=60, n_theta=24, n_r=6):
    """Build the 3-D field, export the VTK volume, and render the 3-D tube views."""
    os.makedirs(outdir, exist_ok=True)
    field = build_3d_field(sv, n_axial=n_axial, n_theta=n_theta, n_r=n_r)
    vtk_path = write_vtk(field, os.path.join(outdir, "pipe_3d.vtk"))
    if not _HAVE_MPL:
        return vtk_path

    # wall-value arrays (azimuth x axial) for the tube renders
    med = lambda A: np.nanmedian(A, 1)
    r = sv.results
    delta = med(r["delta"]); T = med(r["T"])
    a = np.linspace(0.0, 2.0 * math.pi, n_theta)
    wnorm = np.mean((1.0 + 1.6 * np.cos(np.linspace(0, 2 * math.pi, 400, endpoint=False))) / 2.0)
    w_az = np.clip((1.0 + 1.6 * np.cos(a)) / 2.0 / wnorm, 0.05, None)
    depo_wall = np.outer(w_az, delta) * 1000.0                    # (n_theta, n_ax) mm
    Twall = float(sv.case.operating.T_seabed_C)
    # wall temperature ~ between seabed (bottom, water-wetted, coldest) and a bit warmer at top
    temp_wall = np.outer(1.0 - 0.15 * np.cos(a), np.ones_like(T)) * 0 + \
        (Twall + 0.25 * (T[None, :] - Twall) * (0.5 + 0.5 * np.cos(a)[:, None]))

    p1 = _tube_surface(sv, depo_wall, "3-D reconstructed pipe — hydrate wall-deposit distribution",
                       "deposit (mm)", "shct_heat", os.path.join(outdir, "threed_deposit.png"))
    p2 = _tube_surface(sv, temp_wall, "3-D reconstructed pipe — wall temperature distribution",
                       "wall T (°C)", "shct_temp", os.path.join(outdir, "threed_temperature.png"))
    return vtk_path
