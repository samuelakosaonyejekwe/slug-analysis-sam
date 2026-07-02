#!/usr/bin/env python3
# =============================================================================
#  make_animations.py  —  the canonical flow-assurance GIF animations for the
#  deepwater medium-crude slug + hydrate case study.
#
#  Three non-redundant animations per scenario (each shows something the static
#  outputs do NOT):
#    1. anim_flow_line.gif    — flow IN the line: a terrain-following pipe ribbon
#                               coloured by liquid holdup α_l, with slugs
#                               physically travelling along the route in time.
#    2. anim_crosssection.gif — the 2-D pipe bore at the monitor station: the
#                               stratified liquid level and the hydrate deposit
#                               ring closing the bore toward a plug, in time.
#    3. anim_PT_cooldown.gif  — the monitor P–T operating point marching across
#                               the hydrate-stability envelope (into the hydrate
#                               region) as the line cools / responds in time.
#
#  All three are rendered from ONE solver run per scenario.  Style obeys the
#  project rule: NO black / NO dark colours (shct_style palette + shct_* light
#  colormaps); every legend / annotation sits OUTSIDE the data.
#
#  Author: Akosa Samuel Onyejekwe.
# =============================================================================
import os, sys, time, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors

from _paths import CASE                      # shared layout + no-black style hook
import shct_style as S
S.apply_style()
import solver
import shct_crosssection as CX
import run_case_study10 as R                 # reuse the exact case configuration

FPS = 12
DPI = 105

#  scenario -> (build_case variant, t_end_h, output folder). t_end_h matches the
#  committed case_config.json of each scenario.
SCENARIOS = [
    ("asoperated", 48.0, os.path.join(CASE, "outputs_steady")),
    ("shutin",     24.0, os.path.join(CASE, "outputs_shutin")),
    ("mitigated",  48.0, os.path.join(CASE, "outputs_mitigated")),
]

TAN, BROWN = "#E7D7B6", S.BROWN
NAVY, RED, ORANGE, GREEN, TEAL, SKY = S.BLUE, S.RED, S.ORANGE, S.GREEN, S.TEAL, S.SKY
CMAP = "shct_seq"


def _suptitle(fig, text):
    """Centred figure title — kept SHORT so it never overruns the fixed animation
    canvas (bbox='tight' is not applied to animation frames, so an over-wide title
    or any artist outside the canvas would be clipped)."""
    fig.suptitle(text, color=NAVY, fontweight="bold", fontsize=11.5, y=0.965)


def _save(anim, fig, path, fps=FPS):
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=DPI)
    plt.close(fig)
    print(f"    wrote {os.path.basename(path)}  ({os.path.getsize(path)//1024} kB)")


def _subsample(t, k=120):
    """Indices that thin a dense history down to ~k evenly-spaced frames."""
    if t.size <= k:
        return np.arange(t.size)
    return np.linspace(0, t.size - 1, k).astype(int)


# -----------------------------------------------------------------------------
#  1. Flow IN the line — terrain-following pipe ribbon coloured by holdup
# -----------------------------------------------------------------------------
def anim_flow_line(sv, outdir, title):
    r = sv.results
    H, t = r["snap_holdup"], r["snap_t"]
    if not H.size:
        return
    x = sv.x / 1000.0
    z = np.asarray(sv.z, float)
    zr = max(z.max() - z.min(), 1.0)
    half = 0.045 * zr                                   # display half-thickness (exaggerated)
    Xg = np.vstack([x, x])                              # (2, nx)
    Yg = np.vstack([z + half, z - half])               # pipe walls, terrain-following
    norm = mcolors.PowerNorm(gamma=0.5, vmin=float(np.nanmin(H)),
                             vmax=float(np.nanmax(H)))

    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    fig.subplots_adjust(left=0.085, right=0.88, top=0.88, bottom=0.12)
    _suptitle(fig, f"{title} — flow in the line: liquid holdup α_l(x,t)")
    ax.fill_between(x, z - half, z.min() - 0.15 * zr, color=TAN, alpha=.65, zorder=0)
    ax.plot(x, z, color=BROWN, lw=0.8, alpha=.5, zorder=1)
    C0 = np.vstack([H[0], H[0]])
    pcm = ax.pcolormesh(Xg, Yg, C0, cmap=CMAP, norm=norm, shading="gouraud", zorder=2)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(z.min() - 0.15 * zr, z.max() + 0.12 * zr)
    ax.set_xlabel("distance along flowline (km)")
    ax.set_ylabel("seabed elevation (m)")
    cb = fig.colorbar(pcm, ax=ax, pad=.012, fraction=.043)
    cb.set_label("liquid holdup α_l")
    #  time stamp placed in the empty white space ABOVE the pipe (never over the
    #  title or the ribbon), in a light box so it stays legible over the ribbon.
    tag = ax.text(0.015, 0.95, "", transform=ax.transAxes, ha="left", va="top",
                  fontsize=10, color=NAVY, fontweight="bold",
                  bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=S.INK, alpha=.85))
    #  caption kept INSIDE the axes (over the light seabed) — animation frames are
    #  saved at a fixed canvas, so anything outside the axes would be clipped.
    ax.text(0.5, 0.045, "pipe drawn terrain-following; vertical thickness exaggerated",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.5,
            color=S.INK, style="italic")

    def update(k):
        pcm.set_array(np.vstack([H[k], H[k]]).ravel())
        tag.set_text(f"t = {t[k]:5.2f} h")
        return pcm, tag

    anim = animation.FuncAnimation(fig, update, frames=len(t), blit=False)
    _save(anim, fig, os.path.join(outdir, "anim_flow_line.gif"))


# -----------------------------------------------------------------------------
#  2. Pipe cross-section at the monitor — liquid level + hydrate deposit ring
# -----------------------------------------------------------------------------
def anim_crosssection(sv, outdir, title):
    r = sv.results
    ts, tt = r["ts"], r["ts_t"]
    if not tt.size:
        return
    idx = _subsample(tt, 120)
    D = float(sv.case.pipeline.diameter_m)
    R = 0.5 * D
    Tw = float(sv.case.operating.T_seabed_C)
    al = np.clip(np.asarray(ts["alpha_l"], float)[idx], 1e-3, 1 - 1e-3)
    um = np.asarray(ts["j"], float)[idx]
    Tb = np.asarray(ts["T"], float)[idx]
    dl = np.asarray(ts["delta"], float)[idx]                # deposit thickness (m)
    tim = tt[idx]
    vmax = max(1e-3, float(np.nanmax(np.abs(um))) * 1.25)
    seq = plt.get_cmap(CMAP)

    def frame_rgba(k):
        f = CX.reconstruct_section(D, al[k], um[k], Tb[k], Tw, dl[k], n=140)
        n = f["inside"].shape[0]
        img = np.zeros((n, n, 4))                            # transparent background
        gas = f["inside"] & ~f["liquid"] & ~f["deposit"]
        liq = f["liquid"] & ~f["deposit"]
        img[gas] = mcolors.to_rgba(S.SKY, 0.30)             # light gas core
        vnorm = np.clip(np.abs(f["vel"]) / vmax, 0, 1)
        img[liq] = seq(vnorm[liq])                          # liquid coloured by velocity
        img[f["deposit"]] = mcolors.to_rgba(RED, 0.95)      # hydrate deposit ring
        return img, f

    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    fig.subplots_adjust(left=0.08, right=0.72, top=0.88, bottom=0.09)
    #  equal aspect anchored WEST: the square bore sits on the left, leaving real
    #  empty canvas on the right for the legend + read-out (never clipped).
    ax.set_aspect("equal"); ax.set_anchor("W")
    ax.set_xlim(-R * 1.08, R * 1.08); ax.set_ylim(-R * 1.08, R * 1.08)
    ax.set_xlabel("z (m)"); ax.set_ylabel("y (m)")
    _suptitle(fig, f"{title} — pipe bore + hydrate deposit at monitor")
    wall = plt.Circle((0, 0), R, fill=False, color=BROWN, lw=2.0, zorder=5)
    ax.add_patch(wall)
    img0, f0 = frame_rgba(0)
    im = ax.imshow(img0, extent=[-R, R, -R, R], origin="lower", zorder=3,
                   interpolation="bilinear")
    (interface,) = ax.plot([-R, R], [f0["y_int"], f0["y_int"]], color=NAVY, lw=1.2,
                           ls="--", alpha=.7, zorder=4)
    # legend (in the reserved right-hand canvas): proxy handles
    import matplotlib.patches as mpatches
    handles = [mpatches.Patch(color=S.SKY, alpha=.5, label="gas"),
               mpatches.Patch(color=TEAL, label="liquid (colour=velocity)"),
               mpatches.Patch(color=RED, label="hydrate deposit"),
               plt.Line2D([], [], color=BROWN, lw=2, label="pipe wall")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 1.0),
              fontsize=8, borderaxespad=0.0)
    tag = ax.text(1.04, 0.42, "", transform=ax.transAxes, ha="left", va="top",
                  fontsize=8.5, color=NAVY, fontweight="bold")

    def update(k):
        img, f = frame_rgba(k)
        im.set_data(img)
        interface.set_ydata([f["y_int"], f["y_int"]])
        bore_left = D * 1000.0 - 2 * dl[k] * 1000.0
        tag.set_text(f"t = {tim[k]:5.2f} h\nα_l = {al[k]:.2f}\n"
                     f"deposit = {dl[k]*1000:5.1f} mm\nbore left = {max(bore_left,0):5.1f} mm")
        return im, interface, tag

    anim = animation.FuncAnimation(fig, update, frames=len(idx), blit=False)
    _save(anim, fig, os.path.join(outdir, "anim_crosssection.gif"), fps=14)


# -----------------------------------------------------------------------------
#  3. P–T operating point crossing the hydrate-stability envelope, in time
# -----------------------------------------------------------------------------
def anim_PT_cooldown(sv, outdir, title):
    r, c = sv.results, sv.case
    ts, tt = r["ts"], r["ts_t"]
    if not tt.size:
        return
    idx = _subsample(tt, 130)
    P = np.asarray(ts["P"], float)[idx]
    T = np.asarray(ts["T"], float)[idx]
    sub = np.asarray(ts["Tsub"], float)[idx]                # >0 => inside hydrate region
    tim = tt[idx]
    Pc = np.linspace(max(1.0, np.nanmin(P) * 0.6), np.nanmax(P) * 1.12, 200)
    Teq = solver.hydrate_equilibrium_T(Pc, gas_sg=c.fluids.gas_sg,
                                       salinity_wt=c.fluids.salinity_wt,
                                       table=c.fluids.hyd_Teq_table)

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    fig.subplots_adjust(left=0.085, right=0.70, top=0.88, bottom=0.11)
    ax.plot(Teq, Pc, color=RED, lw=2.2)
    ax.fill_betweenx(Pc, Teq.min() - 8, Teq, color=S.HYDFILL, alpha=.55, zorder=0)
    ax.set_xlim(min(np.nanmin(T), Teq.min()) - 2, max(np.nanmax(T), Teq.max()) + 3)
    ax.set_ylim(Pc.min(), Pc.max())
    ax.set_xlabel("temperature T (°C)"); ax.set_ylabel("pressure P (bar)")
    _suptitle(fig, f"{title} — monitor P–T vs hydrate envelope")
    ax.text(0.03, 0.06, "shaded = hydrate-stable\n(left of the red curve)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color=S.INK)
    (trail,) = ax.plot([], [], color=NAVY, lw=1.6, alpha=.85)
    #  the marker colour is a LIVE state indicator (green = safe, red = inside the
    #  hydrate zone), so the legend shows BOTH states rather than a single colour.
    (dot,) = ax.plot([], [], "o", ms=10, color=GREEN, mec=NAVY)
    tag = ax.text(1.04, 0.98, "", transform=ax.transAxes, ha="left", va="top",
                  fontsize=8.5, color=NAVY, fontweight="bold")
    import matplotlib.lines as mlines
    handles = [
        mlines.Line2D([], [], color=RED, lw=2.2, label="hydrate equilibrium T_eq"),
        mlines.Line2D([], [], color=NAVY, lw=1.6, label="P–T trajectory"),
        mlines.Line2D([], [], marker="o", ls="none", ms=9, mfc=GREEN, mec=NAVY,
                      label="now — safe"),
        mlines.Line2D([], [], marker="o", ls="none", ms=9, mfc=RED, mec=NAVY,
                      label="now — in hydrate zone"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 0.62),
              fontsize=8, borderaxespad=0.0)

    def update(k):
        trail.set_data(T[: k + 1], P[: k + 1])
        inside = sub[k] > 0
        dot.set_data([T[k]], [P[k]])
        dot.set_color(RED if inside else GREEN)
        state = "IN HYDRATE ZONE" if inside else "safe"
        tag.set_text(f"t = {tim[k]:5.2f} h\nT = {T[k]:4.1f} °C\nP = {P[k]:5.1f} bar\n"
                     f"ΔT_sub = {sub[k]:+4.1f} °C\n({state})")
        return trail, dot, tag

    anim = animation.FuncAnimation(fig, update, frames=len(idx), blit=False)
    _save(anim, fig, os.path.join(outdir, "anim_PT_cooldown.gif"), fps=14)


# -----------------------------------------------------------------------------
#  4. Riser-region flow behaviour — α_l–P trajectory at the monitor
#     (repeating loops = intermittent / slug flow;  a settled point = stable flow)
# -----------------------------------------------------------------------------
def anim_riser_cycle(sv, outdir, title):
    ts, tt = sv.results["ts"], sv.results["ts_t"]
    if not tt.size:
        return
    #  trim the START-UP PRESSURISATION ramp (P climbing from the initial state to
    #  the operating band) so the phase portrait shows the DEVELOPED flow, not the
    #  fill transient — otherwise every case is dominated by one long diagonal.
    Pall = np.asarray(ts["P"], float)
    Pmed = float(np.nanmedian(Pall))
    i0 = int(np.argmax(Pall >= 0.85 * Pmed)) if np.any(Pall >= 0.85 * Pmed) else 0
    sel = np.arange(i0, tt.size)
    idx = sel[_subsample(tt[sel], 150)]
    hold = np.asarray(ts["alpha_l"], float)[idx]
    P = Pall[idx]
    tim = tt[idx]
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    fig.subplots_adjust(left=0.10, right=0.72, top=0.88, bottom=0.12)
    _suptitle(fig, f"{title} — riser-region monitor: α_l–P trajectory")
    hpad = max(0.02, 0.05 * (np.nanmax(hold) - np.nanmin(hold)))
    ax.set_xlim(max(0.0, np.nanmin(hold) - hpad), min(1.0, np.nanmax(hold) + hpad))
    ax.set_ylim(np.nanmin(P) - 1, np.nanmax(P) + 1)
    ax.set_xlabel("liquid holdup  α_l  (monitor)")
    ax.set_ylabel("pressure  P (bar)")
    (trail,) = ax.plot([], [], color=S.BLUE, lw=1.2, alpha=.8)
    (dot,) = ax.plot([], [], "o", ms=9, color=ORANGE, mec=NAVY)
    tag = ax.text(1.03, 0.98, "", transform=ax.transAxes, ha="left", va="top",
                  fontsize=8.5, color=NAVY, fontweight="bold")
    ax.text(1.03, 0.58, "repeating loops =\nintermittent (slug)\nflow;\na settled point =\nstable flow",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=S.INK)

    def update(k):
        trail.set_data(hold[: k + 1], P[: k + 1])
        dot.set_data([hold[k]], [P[k]])
        tag.set_text(f"t = {tim[k]:5.2f} h\nα_l = {hold[k]:.2f}\nP = {P[k]:5.1f} bar")
        return trail, dot, tag

    anim = animation.FuncAnimation(fig, update, frames=len(idx), blit=False)
    _save(anim, fig, os.path.join(outdir, "anim_riser_cycle.gif"), fps=14)


# -----------------------------------------------------------------------------
#  5. P(x,t) / T(x,t) profile wave marching along the line
# -----------------------------------------------------------------------------
def anim_profile_wave(sv, outdir, title):
    r = sv.results
    SP, ST, st = r.get("snap_P"), r.get("snap_T"), r["snap_t"]
    if SP is None or ST is None or not np.asarray(SP).size:
        return
    SP, ST = np.asarray(SP, float), np.asarray(ST, float)
    x = np.asarray(sv.x, float) / 1000.0
    Tsea = float(sv.case.operating.T_seabed_C)
    fig, ax = plt.subplots(2, 1, figsize=(9.0, 5.8), sharex=True)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.9, bottom=0.1, hspace=0.12)
    _suptitle(fig, f"{title} — P(x,t) & T(x,t) profile wave along the line")
    ax[0].set_ylim(min(Tsea - 2, float(np.nanmin(ST))) - 1, float(np.nanmax(ST)) + 3)
    ax[0].set_ylabel("temperature T (°C)")
    ax[0].axhline(Tsea, color=S.TEAL, ls=":", lw=1.1, label=f"seabed {Tsea:.0f} °C")
    ax[0].legend(loc="upper right", fontsize=7.5)
    ax[1].set_ylim(float(np.nanmin(SP)) - 2, float(np.nanmax(SP)) + 2)
    ax[1].set_ylabel("pressure P (bar)"); ax[1].set_xlabel("distance along flowline (km)")
    ax[0].set_xlim(x.min(), x.max())
    (lT,) = ax[0].plot(x, ST[0], color=RED, lw=1.8)
    (lP,) = ax[1].plot(x, SP[0], color=NAVY, lw=1.8)
    tag = ax[0].text(0.015, 0.05, "", transform=ax[0].transAxes, ha="left", va="bottom",
                     fontsize=9, color=NAVY, fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=S.INK, alpha=.85))

    def update(k):
        lT.set_ydata(ST[k]); lP.set_ydata(SP[k])
        tag.set_text(f"t = {st[k]:5.2f} h")
        return lT, lP, tag

    anim = animation.FuncAnimation(fig, update, frames=len(st), blit=False)
    _save(anim, fig, os.path.join(outdir, "anim_profile_wave.gif"))


# -----------------------------------------------------------------------------
#  Render-data cache — the solver run is the only expensive step, so the minimal
#  arrays the three animations need are cached to a gitignored .npz.  With the
#  cache present, `--render-only` rebuilds every GIF in seconds (no re-solve),
#  which is how cosmetic tweaks are iterated.
# -----------------------------------------------------------------------------
import types
CACHE = os.path.join(CASE, "scripts", "_anim_cache")
PRETTY = {"asoperated": "As-operated (steady)",
          "shutin": "Shut-in cooldown",
          "mitigated": "Mitigated (insulation + MEG)"}
#  SHORT labels for the on-figure titles (the long PRETTY names overran the canvas)
SHORT = {"asoperated": "As-operated", "shutin": "Shut-in", "mitigated": "Mitigated"}


def _cache_path(variant):
    return os.path.join(CACHE, f"{variant}.npz")


def _save_cache(variant, sv):
    os.makedirs(CACHE, exist_ok=True)
    r, ts = sv.results, sv.results["ts"]
    f = sv.case.fluids
    tbl = f.hyd_Teq_table
    np.savez(_cache_path(variant),
             snap_holdup=r["snap_holdup"], snap_t=r["snap_t"], x=sv.x, z=sv.z,
             snap_P=r.get("snap_P", np.zeros((0, 0))), snap_T=r.get("snap_T", np.zeros((0, 0))),
             ts_t=r["ts_t"], ts_alpha_l=ts["alpha_l"], ts_j=ts["j"], ts_T=ts["T"],
             ts_delta=ts["delta"], ts_P=ts["P"], ts_Tsub=ts["Tsub"],
             diameter_m=sv.case.pipeline.diameter_m,
             T_seabed_C=sv.case.operating.T_seabed_C,
             gas_sg=f.gas_sg, salinity_wt=f.salinity_wt,
             has_table=(tbl is not None),
             table=(np.asarray(tbl, float) if tbl is not None else np.zeros((0, 2))))


def _load_shim(variant):
    """Rebuild a lightweight object exposing exactly what the render fns read."""
    d = np.load(_cache_path(variant))
    results = dict(snap_holdup=d["snap_holdup"], snap_t=d["snap_t"],
                   snap_P=d["snap_P"], snap_T=d["snap_T"],
                   ts_t=d["ts_t"],
                   ts=dict(alpha_l=d["ts_alpha_l"], j=d["ts_j"], T=d["ts_T"],
                           delta=d["ts_delta"], P=d["ts_P"], Tsub=d["ts_Tsub"]))
    fluids = types.SimpleNamespace(gas_sg=float(d["gas_sg"]),
                                   salinity_wt=float(d["salinity_wt"]),
                                   hyd_Teq_table=(d["table"] if bool(d["has_table"]) else None))
    case = types.SimpleNamespace(
        pipeline=types.SimpleNamespace(diameter_m=float(d["diameter_m"])),
        operating=types.SimpleNamespace(T_seabed_C=float(d["T_seabed_C"])),
        fluids=fluids)
    return types.SimpleNamespace(results=results, x=d["x"], z=d["z"], case=case)


def _render_all(sv, outdir, pretty):
    os.makedirs(outdir, exist_ok=True)
    anim_flow_line(sv, outdir, pretty)
    anim_crosssection(sv, outdir, pretty)
    anim_PT_cooldown(sv, outdir, pretty)
    anim_riser_cycle(sv, outdir, pretty)
    anim_profile_wave(sv, outdir, pretty)


def run_scenario(variant, t_end_h, outdir, render_only=False):
    short = SHORT[variant]
    if render_only:
        if not os.path.exists(_cache_path(variant)):
            print(f"[{variant}] no cache — run a full solve first; skipping.", flush=True)
            return
        print(f"[{variant}] render-only from cache -> {os.path.relpath(outdir, CASE)}", flush=True)
        _render_all(_load_shim(variant), outdir, short)
        return
    print(f"[{variant}] running solver -> {os.path.relpath(outdir, CASE)}", flush=True)
    t0 = time.time()
    case = R.build_case(f"anim-{variant}", variant, t_end_h=t_end_h,
                        n_ensemble=12, n_cells=70)
    sv = solver.TransientSHCT(case)
    sv.run(verbose=False)
    try:
        km = sv.engineering()
        print(f"    solved in {time.time()-t0:5.1f}s  |  mass_err="
              f"{km.get('mass_conservation_err', float('nan')):.2e}"
              f"  gas_mass_err={km.get('gas_mass_conservation_err', float('nan')):.2e}"
              f"  fallbacks={int(km.get('fallbacks', 0))}", flush=True)
    except Exception as e:                                   # pragma: no cover
        print(f"    solved in {time.time()-t0:5.1f}s  |  (metrics check skipped: {e})", flush=True)
    try:
        _save_cache(variant, sv)
    except Exception as e:                                   # pragma: no cover
        print(f"    (cache save skipped: {e})", flush=True)
    _render_all(sv, outdir, short)


def main():
    args = sys.argv[1:]
    render_only = "--render-only" in args
    want = [a for a in args if not a.startswith("-")]
    for variant, t_end_h, outdir in SCENARIOS:
        if want and variant not in want:
            continue
        run_scenario(variant, t_end_h, outdir, render_only=render_only)
    print("done.", flush=True)


if __name__ == "__main__":
    main()
