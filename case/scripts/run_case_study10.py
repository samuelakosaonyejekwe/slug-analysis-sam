#!/usr/bin/env python3
# =============================================================================
#  run_case_study10.py  —  Industrial Case Study #3
#  HIGH slug- AND hydrate-risk deepwater MEDIUM-CRUDE-OIL subsea tie-back.
# -----------------------------------------------------------------------------
#  A representative deepwater West-Africa-style tie-back (32 km flowline + steel
#  catenary riser, ~1100 m water depth) carrying a ~30 deg API medium crude oil with
#  produced water over a cold, undulating seabed.  This geometry + fluid is a
#  textbook combination for BOTH hydrodynamic/terrain/severe-riser slugging AND
#  hydrate formation, so it exercises the whole SHCT prediction chain.
#
#  Three engineering scenarios are run end-to-end through the real solver:
#    (A) outputs_steady/    — normal production, AS-OPERATED (degraded / flooded
#                             wet insulation, U_wall ~ 20 W/m2K), NO inhibitor
#                             -> the high-risk prediction (slug + hydrate + plug).
#    (B) outputs_shutin/    — unplanned shut-in cooldown -> no-touch time.
#    (C) outputs_mitigated/ — engineered fix: restored multi-layer insulation
#                             (U_eff = 3.45 W/m2K from the layered cylindrical
#                             resistance) + continuous MEG -> shows the model used
#                             as a DESIGN tool (risk removed).
#
#  The steady run also drives the full advanced stack: compositional Peng-Robinson
#  EOS PVT, cross-section / quasi-3-D reconstruction (+ VTK), compositional
#  transport, OpenFOAM (interFoam) coupling case generation, and hydrate-curve
#  validation against published experimental data.  Everything (CSV tables,
#  engineering deliverables, charts/curves, JSON, bespoke slug/hydrate/mitigation
#  figures, and the input-data deck) is written under case/outputs_<scenario>/.
#
#  NOTE ON DATA PROVENANCE (honest framing): the field is a *representative
#  industrial archetype*.  Geometry, fluid and operating parameters are realistic
#  and self-consistent values typical of deepwater medium-crude-oil tie-backs in the
#  open literature; they are NOT proprietary operator data.  The PHYSICS and the
#  PREDICTIONS are produced by the real solver, and the hydrate thermodynamics are
#  validated against published experimental data (see §validation).
# =============================================================================
import math
import os

import matplotlib
import numpy as np
from _paths import OUT, ROOT  # shared layout + no-black style (shct_style)

import shct_compositional
import shct_compositional_sim
import shct_crosssection
import shct_openfoam
import shct_spacetime
import shct_threed
import solver

matplotlib.use("Agg")
from dataclasses import asdict

import matplotlib.pyplot as plt

import shct_style as S

#  DPI follows SHCT_FIG_DPI (default 320) so every generated figure meets the
#  journal artwork minimum of 300 dpi; a hard-coded 150/155 silently fell short.
#  Read from shct_style rather than re-parsing the environment variable here: the
#  default lives in exactly one place, so the two cannot drift apart.
_FIG_DPI = S.FIG_DPI

# medium, non-black, non-dark palette (see shct_style.py)
NAVY, ACC, ORG, RED, GRN, TEAL = S.BLUE, "#1F8AC0", S.ORANGE, S.RED, S.GREEN, S.TEAL

# -----------------------------------------------------------------------------
#  Fluid: a MEDIUM CRUDE OIL — a ~30 deg API black oil (compositional makeup
#  -> Peng-Robinson flash).  Moderate-GOR with a substantial heavy C7+ tail
#  (~31 mol%): the dissolved/associated gas (C1 ~43 mol%) liberates along the
#  cold line and, with the long near-horizontal undulating step-out, drives
#  HYDRODYNAMIC (and terrain/severe-riser) SLUGGING; the 70 % late-life water cut
#  plus the water-wet associated gas in the cold deepwater wall drives HYDRATE
#  formation.
#  This crude is a textbook combination for BOTH slugging and hydrates.
# -----------------------------------------------------------------------------
#  RE-CUT TO THE CRUDE THIS CASE ACTUALLY TRANSPORTS. The previous cut -- 43 mol% C1 and
#  31 mol% C7+ -- is a volatile oil: flashed at line conditions it gives a 549 kg/m3
#  liquid, while every velocity, holdup, pressure drop and deposit in this case is driven
#  by the flow model's rho_oil = 858 kg/m3 medium crude. Two fluids in one case study, and
#  the compositional figures were describing the wrong one.
#
#  No assay was invented to fix it. The plus fraction was re-characterised (MW 250 g/mol,
#  SG 0.86, Riazi-Daubert + Kesler-Lee -- see shct_eos.COMPONENTS), the light/heavy split
#  moved to a black-oil cut, and the Peneloux shift regressed against the case's own stated
#  oil density. The flash now returns 858.0 kg/m3 at 120 bar / 20 C -- the number the flow
#  model was already using -- so the two halves of the case finally describe one fluid.
#  Intermediates keep their original relative proportions.
#  CONSTRAINED ON TWO PROPERTIES, not one. The previous cut (C1 0.300, C7+ 0.480) was
#  regressed to reproduce the case's stated oil density and nothing else, and it hit that
#  exactly -- 858 kg/m3 -- while implying almost no dissolved gas. Its bubble point ran
#  98.6 bar at the 55 C inlet falling to 70.1 bar at 10 C, against a line that runs 150
#  bar down to 65: the pressure never reached the bubble point until the last node, so the
#  EOS held the hydrocarbon single-phase at 39 of 40 stations while the flow model beside
#  it ran 52 % gas. A 32 km tie-back that flashes only at its outlet is not a normal
#  result, and the phase envelope, the K-value panel and every compositional figure showed
#  it.
#
#  Density alone does not pin a composition: an oil can be made to weigh 858 kg/m3 with
#  almost any gas content by trading C1 against C7+. Fixing the light end as well is what
#  makes the cut unique. C1 0.600 puts the bubble point above the line pressure over the
#  whole route (two-phase at 40 of 40 stations, GVF 36 % at the inlet rising to 62 % at
#  the outlet, mean 32 %), which is the same character as the flow model's 52 %.
#
#  The density is NOT given up to get it: the Peneloux shift on C7+ is re-regressed to
#  0.2672 in shct_eos.py, still inside the Jhaveri-Youngren C7+ range of 0.1-0.3, and the
#  flash returns 858.0 kg/m3 exactly. Nor is the case re-described -- C7+ is still 81 % of
#  the mass and the liquid is still a 33 API medium crude; only its dissolved gas changes,
#  from implausibly little to a plausible 11 % of the mass.
#
#  What remains open: the mean GVF is 32 % against the flow model's 52 %. Closing that
#  last gap needs C1 near 0.70, which WOULD make this a volatile oil, so it is left as a
#  stated difference rather than bought with a re-description.
CRUDE_OIL = {
    "N2": 0.00197, "CO2": 0.00982, "C1": 0.60000, "C2": 0.03683, "C3": 0.02848,
    "iC4": 0.00589, "nC4": 0.01375, "iC5": 0.00638, "nC5": 0.00786, "C6": 0.01473,
    "C7+": 0.27429,
}

#  TERMINOLOGY, used consistently in every figure and table:
#     tie-back / route  — the whole 32 km, wellhead to host
#     flowline          — the seabed portion, 0 .. 30.2 km  (the first 94.5 %)
#     riser             — the steel catenary, 30.2 .. 32 km, climbing from the
#                         riser base at ~1085 m depth to the host at ~25 m
#  Every along-route axis is therefore "distance from wellhead [km]" and spans
#  the whole tie-back; only the riser-specific figure uses "depth from host [m]".
WATER_DEPTH_M = 1100.0
RISER_FRAC = 0.945          # last ~5.5% of the route is the steel catenary riser


# -----------------------------------------------------------------------------
#  Terrain: a long, strongly undulating deepwater seabed (multiple low spots
#  -> terrain slugging) climbing into a steep SCR riser (-> severe riser slugging)
# -----------------------------------------------------------------------------
def build_elevation(n, length_m):
    x = np.linspace(0.0, length_m, n)
    xk = x / length_m
    seabed = (-WATER_DEPTH_M + 48.0 * xk
              - 30.0 * np.sin(2 * math.pi * 3.3 * xk)
              - 18.0 * np.sin(2 * math.pi * 6.7 * xk + 0.5)
              - 11.0 * np.sin(2 * math.pi * 10.4 * xk + 1.1))
    z = seabed.copy()
    riser = xk > RISER_FRAC
    if riser.any():
        z[riser] = np.linspace(seabed[riser][0], -25.0, int(riser.sum()))   # climb to host turret
    return [float(v) for v in z]


# -----------------------------------------------------------------------------
#  The case builder.  variant: "asoperated" | "shutin" | "mitigated"
# -----------------------------------------------------------------------------
def build_case(name, variant, t_end_h, n_ensemble=12, n_cells=70):
    c = solver.Case()
    c.name = name

    p = c.pipeline
    p.length_m = 32_000.0                 # 32 km step-out
    p.diameter_m = 0.2545                 # ~10.75-in carbon-steel flowline ID
    p.roughness_m = 4.6e-5
    p.n_cells = n_cells
    p.elevation_m = build_elevation(n_cells, p.length_m)
    p.h_inner = 1500.0
    p.h_outer = 300.0

    f = c.fluids
    f.rho_oil = 858.0                     # live medium crude (~30 deg API black oil)
    f.rho_water = 1025.0
    #  LATE-LIFE conditions, and deliberately so. Wall deposition scales with the water
    #  that reaches the wall — a_wall = (4/D)*alpha_l*water_frac — so the water cut is a
    #  first-order control on it, and a mature deepwater field genuinely reaches 70 %.
    #  At the earlier 35 % and full rate this line is SUB-CRITICAL: Phi_SH = 0.27 against
    #  the derived Phi_crit = 1.08, a 4.2 mm stable deposit, and no plug in 48 h. That is
    #  reported in the paper rather than hidden, because a criterion that called every
    #  line critical would be worth nothing; the discrimination is the evidence.
    f.water_cut = 0.70                    # 70 % water cut — mature field, late life
    f.mu_liquid = 5.0e-3                  # ~5 cP medium-crude live-oil viscosity
    f.mu_gas = 1.3e-5
    f.sigma = 0.022
    f.salinity_wt = 4.5                   # saline formation water (depresses hydrate Teq)
    f.cp_liquid = 2050.0
    f.cp_gas = 2300.0
    f.composition = dict(CRUDE_OIL)       # -> compositional Peng-Robinson PVT
    f.condensation_latent = True          # latent heat of the lighter-end condensation
    f.gas_visc_corr = True                # Lee real-gas viscosity
    f.oil_pvt_corr = True                 # P,T oil density
    f.wax_appearance_C = 32.0             # waxier medium crude -> screen wax risk too

    o = c.operating
    #  0.6x the design rate, again late-life. Reduced throughput lowers the slug frequency
    #  that scours the wall and lets the line run colder, so it raises Phi_SH from both
    #  sides. It is NOT taken lower than this: below about 0.6x the line stops slugging,
    #  f_slug collapses onto its assumed 1e-4 Hz floor, and Phi_SH jumps from ~1.3 to
    #  ~1200 — an artefact of a numerical guard rather than a physical result, and no
    #  headline should rest on it.
    o.q_liquid_insitu = 0.055 * 0.60      # in-situ liquid (oil + water) rate (m3/s)
    o.q_gas_insitu_inlet = 0.150 * 0.60   # in-situ associated-gas rate (m3/s)
    o.P_inlet_bar = 150.0
    o.T_inlet_C = 58.0
    o.T_seabed_C = 4.0

    n = c.numerics
    n.t_end_h = t_end_h
    #  (the grid size is pipeline.n_cells, set above; Numerics carries no n_cells
    #  field, so assigning one here attached a stray attribute nothing ever read)
    n.n_ensemble = n_ensemble
    n.n_snapshots = 240          # dense space-time history for the published-scheme fields
    n.seed = 13

    sc = c.scenario
    sc.event_time_h = 6.0
    sc.shutin_residual = 0.02

    # ---- variant-specific thermal design & inhibition ----
    if variant in ("asoperated", "shutin"):
        #  Degraded / water-flooded wet insulation: the design intent is lost and
        #  the line behaves close to bare steel in cold seawater -> the hydrate threat.
        o.U_wall = 22.0
        o.MEG_wt_inlet = 0.0              # no inhibitor -> let the model PREDICT the requirement
        sc.kind = "shutin" if variant == "shutin" else "steady"
    elif variant == "mitigated":
        #  Engineered fix: restored multi-layer insulation (steel + syntactic PP foam +
        #  outer coating) -> low effective U, plus continuous MEG for transients/shut-in.
        p.wall_layers = [
            [0.0254, 45.0, 3.9e6],        # 25.4 mm carbon-steel wall
            [0.060, 0.16, 1.1e5],         # 60 mm syntactic-foam / wet insulation
            [0.012, 0.30, 1.4e5],         # 12 mm outer coating
        ]
        o.U_wall = 22.0                   # ignored when wall_layers is set, kept for the record
        o.MEG_wt_inlet = 25.0            # continuous MEG (design dose, with margin)
        sc.kind = "steady"
    else:
        raise ValueError(variant)

    return solver.validate_case(c)


# -----------------------------------------------------------------------------
#  Bespoke chart 1 — slug-formation prediction (terrain / regime / freq / holdup)
# -----------------------------------------------------------------------------
def slug_chart(sv, outdir):
    r = sv.results
    x = sv.x / 1000.0
    def med(A):
        return np.nanmedian(A, 1)
    reg = med(r["regime"]); fsl = med(r["fslug"]); hold = med(r["alpha_l"])
    Lu = solver.slug_length(med(r["j"]), sv.case.pipeline.diameter_m, fsl)

    fig, ax = plt.subplots(3, 1, figsize=(8, 7.2), sharex=True)
    ax[0].fill_between(x, sv.z, sv.z.min() - 20, color=S.TAN, alpha=.55)
    ax[0].plot(x, sv.z, color=S.BROWN, lw=1.2); ax[0].set_ylabel("elevation (m)")
    ax[0].set_title(solver._ttl("Slug-formation prediction — deepwater medium-crude-oil tie-back"),
                    color=NAVY, fontweight="bold")
    ax[1].plot(x, hold, color=ACC, lw=1.8, label="liquid holdup α_l")
    ax[1].set_ylabel("holdup α_l"); ax[1].set_ylim(0, 1)
    sl = np.isin(np.round(reg), [2, 5])
    ax[1].fill_between(x, 0, 1, where=sl, color="#f6d6d2", alpha=.5,
                       transform=ax[1].get_xaxis_transform(), label="intermittent (slug/churn)")
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    lf, = ax[2].plot(x, fsl, color=ORG, lw=1.8, label="slug frequency f_slug (Hz)")
    ax[2].set_ylabel("f_slug (Hz)", color=ORG); ax[2].set_xlabel("distance from wellhead  [km]")
    a2 = ax[2].twinx(); ll, = a2.plot(x, Lu, color=NAVY, lw=1.4, ls="--", label="slug-unit length (m)")
    a2.set_ylabel("slug length (m)", color=NAVY)
    ax[2].legend(handles=[lf, ll], fontsize=8, loc="upper left", bbox_to_anchor=(1.13, 1.0),
                 borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/09_slug_prediction.png", dpi=_FIG_DPI); plt.close(fig)


# -----------------------------------------------------------------------------
#  Bespoke chart 2 — severe-slugging RISER zoom (last ~6 km)
# -----------------------------------------------------------------------------
def riser_chart(sv, outdir):
    r = sv.results; x = sv.x / 1000.0

    def med(A):
        return np.nanmedian(A, 1)
    m = x >= (x.max() - 6.0)
    hold = med(r["alpha_l"]); reg = np.round(med(r["regime"]))
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 5.6), sharex=True)
    ax[0].fill_between(x[m], sv.z[m], sv.z[m].min() - 10, color=S.TAN, alpha=.55)
    ax[0].plot(x[m], sv.z[m], color=S.BROWN, lw=1.4); ax[0].set_ylabel("elevation (m)")
    ax[0].set_title(solver._ttl("Severe-slugging screen — riser base & ascent (last 6 km)"),
                    color=NAVY, fontweight="bold")
    ax[1].plot(x[m], hold[m], color=ACC, lw=2.0, label="liquid holdup α_l")
    sl = np.isin(reg, [2, 5])
    ax[1].fill_between(x[m], 0, 1, where=sl[m], color="#f6d6d2", alpha=.5,
                       transform=ax[1].get_xaxis_transform(), label="intermittent (slug/churn)")
    ax[1].set_ylabel("holdup α_l")
    _hpk = float(np.nanmax(hold[m])) if np.any(m) else 0.0
    ax[1].set_ylim(0, 1.14 if _hpk > 0.9 else 1.0)
    if _hpk > 0.9:
        _ipk = int(np.nanargmax(np.where(m, hold, np.nan)))
        ax[1].annotate(f"riser fills: α_l = {_hpk:.3f} at {float(x[_ipk]):.1f} km",
                       xy=(float(x[_ipk]), _hpk), xytext=(0.02, 0.90),
                       textcoords="axes fraction", fontsize=7.5, color=NAVY, ha="left",
                       arrowprops={"arrowstyle": "->", "color": NAVY, "lw": 1.1})
    ax[1].set_xlabel("distance from wellhead  [km]")
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/10_riser_severe_slug.png", dpi=_FIG_DPI); plt.close(fig)


# -----------------------------------------------------------------------------
#  Bespoke chart 3 — hydrate P-T envelope with BOTH the production and shut-in
#  trajectories overlaid (shows the line driven INTO the hydrate region)
# -----------------------------------------------------------------------------
def hydrate_envelope_chart(sv_op, sv_si, outdir):
    c = sv_op.case

    def med(A):
        return np.nanmedian(A, 1)
    rO, rS = sv_op.results, sv_si.results
    Pc = np.linspace(5.0, max(180.0, med(rO["p"]).max() * 1.1), 200)
    Tc = solver.hydrate_equilibrium_T(Pc, gas_sg=c.fluids.gas_sg,
                                      salinity_wt=c.fluids.salinity_wt, table=c.fluids.hyd_Teq_table)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.plot(Tc, Pc, color=RED, lw=2.4, label="hydrate equilibrium (T_eq)")
    ax.plot(med(rO["T"]), med(rO["p"]), color=NAVY, lw=2.0, marker="o", ms=2.5,
            label="production trajectory (as-operated)")
    ax.plot(med(rS["T"]), med(rS["p"]), color=TEAL, lw=2.0, ls="--", marker="s", ms=2.5,
            label="shut-in trajectory")
    #  follow the data rather than a fixed window: the mitigated / insulated line
    #  runs well above 32 degC and would otherwise sit off the plot entirely.
    _T = np.concatenate([Tc, med(rO["T"]), med(rS["T"])])
    _P = np.concatenate([Pc, med(rO["p"]), med(rS["p"])])
    _T = _T[np.isfinite(_T)]; _P = _P[np.isfinite(_P)]
    _pad = max(0.05 * (float(_T.max()) - float(_T.min())), 1.0)
    _xlo = float(_T.min()) - _pad
    ax.set_xlim(_xlo, float(_T.max()) + _pad)
    ax.set_ylim(0, float(_P.max()) * 1.08)
    #  SHADE FROM THE AXIS, NOT FROM 0 C. Hydrate is stable everywhere COLDER than Teq(P),
    #  so a fill anchored at T = 0 leaves the sub-zero part of the region unshaded — and
    #  this one is LABELLED "hydrate stability region", so a reader trusting the label
    #  would read the cold end as safe. Drawn after set_xlim so the axis edge is known.
    ax.fill_betweenx(Pc, _xlo, Tc, color="#f6d6d2", alpha=.45,
                     label="hydrate stability region (T < T_eq)")
    ax.set_xlabel("temperature (°C)"); ax.set_ylabel("pressure (bar)")
    ax.set_title(solver._ttl("Hydrate-formation prediction — P–T trajectories vs envelope"),
                 color=NAVY, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/11_hydrate_envelope.png", dpi=_FIG_DPI); plt.close(fig)


# -----------------------------------------------------------------------------
#  Bespoke chart 4 — mitigation comparison (as-operated vs engineered fix)
# -----------------------------------------------------------------------------
def mitigation_chart(eng_base, eng_mit, outdir):
    metrics = [
        ("Max subcooling (°C)", "max_subcooling_C"),
        ("Peak deposit (mm)", "peak_deposit_mm"),
        ("Plug probability (%)", None),         # handled specially
        ("Under-inhibited (km)", "under_inhibited_km"),
        ("No-touch time (h)", "cooldown_to_hydrate_h"),
    ]
    labels, base_v, mit_v = [], [], []
    for lab, key in metrics:
        labels.append(lab)
        if lab.startswith("Plug"):
            base_v.append(eng_base["P_plug"] * 100.0); mit_v.append(eng_mit["P_plug"] * 100.0)
        else:
            base_v.append(float(eng_base.get(key, 0.0))); mit_v.append(float(eng_mit.get(key, 0.0)))
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    yp = np.arange(len(labels)); h = 0.38
    ax.barh(yp + h / 2, base_v, h, color=RED, label="as-operated (degraded insulation, no MEG)")
    ax.barh(yp - h / 2, mit_v, h, color=GRN, label="engineered fix (insulation + MEG)")
    ax.set_yticks(yp); ax.set_yticklabels(labels)
    ax.set_xscale("symlog", linthresh=1.0)
    for i, (b, m) in enumerate(zip(base_v, mit_v)):
        ax.text(b, i + h / 2, f" {b:.1f}", va="center", fontsize=8, color=RED)
        ax.text(m, i - h / 2, f" {m:.1f}", va="center", fontsize=8, color=GRN)
    ax.set_xlabel("value (symlog scale)")
    ax.set_title(solver._ttl("Mitigation comparison — model used as a flow-assurance design tool"),
                 color=NAVY, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/12_mitigation_comparison.png", dpi=_FIG_DPI); plt.close(fig)


# -----------------------------------------------------------------------------
#  Save the full input-data deck (CSV + JSON) for the report
# -----------------------------------------------------------------------------
def save_input_deck(case, outdir):
    p, f, o, n = case.pipeline, case.fluids, case.operating, case.numerics
    rows = [
        ("Pipeline", "Route length", p.length_m / 1000.0, "km"),
        ("Pipeline", "Internal diameter", p.diameter_m * 1000.0, "mm"),
        ("Pipeline", "Wall roughness", p.roughness_m * 1e6, "µm"),
        ("Pipeline", "Numerical cells", p.n_cells, "-"),
        ("Pipeline", "Water depth (riser base)", WATER_DEPTH_M, "m"),
        ("Pipeline", "Riser fraction of route", (1 - RISER_FRAC) * 100, "%"),
        ("Fluid", "Live-oil density", f.rho_oil, "kg/m3"),
        ("Fluid", "Water density", f.rho_water, "kg/m3"),
        ("Fluid", "Water cut", f.water_cut * 100, "%"),
        ("Fluid", "Liquid viscosity", f.mu_liquid * 1000, "cP"),
        ("Fluid", "Interfacial tension", f.sigma, "N/m"),
        ("Fluid", "Formation-water salinity", f.salinity_wt, "wt% NaCl-eq"),
        ("Fluid", "Wax appearance temperature", f.wax_appearance_C, "°C"),
        ("Fluid", "PVT model", "Peng-Robinson compositional flash", "-"),
        ("Operating", "In-situ gas rate", o.q_gas_insitu_inlet, "m3/s"),
        ("Operating", "In-situ liquid rate", o.q_liquid_insitu, "m3/s"),
        ("Operating", "Inlet pressure", o.P_inlet_bar, "bar"),
        ("Operating", "Inlet temperature", o.T_inlet_C, "°C"),
        ("Operating", "Seabed temperature", o.T_seabed_C, "°C"),
        ("Operating", "Wall-loss coefficient U", o.U_wall, "W/m2K"),
        ("Operating", "MEG injected at inlet", o.MEG_wt_inlet, "wt%"),
        ("Numerics", "Engine", n.engine, "-"),
        ("Numerics", "Simulated time", n.t_end_h, "h"),
        ("Numerics", "Ensemble realisations", n.n_ensemble, "-"),
        ("Numerics", "Random seed", n.seed, "-"),
    ]
    #  Round the numbers to what the deck is quoting them to. This is a human-readable
    #  data deck, and (1 - RISER_FRAC) * 100 was writing 5.500000000000005 % into it --
    #  binary floating point showing through a table a reader is meant to take at face
    #  value. Twelve significant figures is far beyond any input's precision here, so it
    #  cannot round away a real difference; strings and ints pass through untouched.
    def _tidy(v):
        return round(v, 12) if isinstance(v, float) else v
    solver._save_csv(f"{outdir}/input_data_deck.csv",
                     ["group", "parameter", "value", "units"],
                     [[g, pname, _tidy(v), u] for g, pname, v, u in rows], all_str=True)
    comp_rows = [[k, f"{v:.4f}"] for k, v in CRUDE_OIL.items()]
    solver._save_csv(f"{outdir}/feed_composition.csv", ["component", "mol_fraction"],
                     comp_rows, all_str=True)
    with open(os.path.join(outdir, "case_config.json"), "w") as fh:
        solver.dump_json(asdict(case), fh)


# -----------------------------------------------------------------------------
#  Core: run a scenario, dump every standard output
# -----------------------------------------------------------------------------
def run_core(case, outdir, slug=True, riser=True):
    os.makedirs(outdir, exist_ok=True)
    #  CONSTRUCT FIRST, THEN DUMP THE DECK. TransientSHCT.__init__ RESOLVES the fluid:
    #  it builds the Peng-Robinson PVT surface and replaces gas_sg / gas_MW with the
    #  gravity of the vapour the feed actually releases. Dumping the case before that
    #  wrote a case_config.json describing a case nobody ran -- gas_sg 0.60 and
    #  gas_MW 0.019 (the dataclass defaults) instead of the 0.6307 / 0.018268 the run
    #  used, and pvt_table null instead of the surface it built.
    #
    #  That is not cosmetic. shct_spacetime._State rebuilds a Case from this file so a
    #  rerender sees "the true fluid"; with the defaults in it, the restored gas density
    #  at the outlet came out 11 % light (66.7 vs 74.1 kg/m3), and the Kelvin-Helmholtz
    #  margin of figure 27 read 1.96 instead of the 1.99 the live run draws.
    sv = solver.TransientSHCT(case); save_input_deck(case, outdir); sv.run(verbose=True)
    eng = sv.engineering()
    solver.write_tables(sv, eng, outdir)
    #  THE SUSTAINED-Phi_SH PROFILE IS A CASE-STUDY OUTPUT and must come out of the case
    #  study. It is tracked in every scenario folder, but the only thing that wrote it was
    #  rerun_sustained.py -- a separate script nobody has to run -- so deleting the outputs
    #  and re-running this driver, which is exactly what a reproducibility check does, left
    #  a tracked artefact behind with no way to regenerate it from the documented pipeline.
    #  Same fields, same formatting; rerun_sustained.py still writes it too, and now agrees.
    _snapP = np.asarray(sv.results["snap_PhiSH"], float)
    if _snapP.size:
        _forming = np.nanmax(sv.results["max_Tsub"], axis=1) > 0.0
        _sust = np.where(_forming, np.nanmedian(_snapP, axis=0), np.nan)
        solver._save_csv(
            os.path.join(outdir, "sustained_phiSH_profile.csv"),
            ["x_km", "hydrate_forming", "Phi_SH_sustained", "Phi_SH_running_max",
             "Phi_SH_final"],
            [[f"{a:.4f}", int(m), f"{b:.6g}", f"{c:.6g}", f"{d:.6g}"]
             for a, m, b, c, d in zip(sv.x / 1000.0, _forming, _sust,
                                      np.nanmedian(sv.results["max_PhiSH"], 1), _snapP[-1])],
            all_str=True)
    solver.make_charts(sv, eng, outdir)
    shct_spacetime.spacetime_outputs(sv, eng, outdir)
    if slug:
        slug_chart(sv, outdir)
    if riser:
        riser_chart(sv, outdir)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        solver.dump_json(eng, fh)
    with open(os.path.join(outdir, "key_metrics.json"), "w") as fh:
        #  bool is a subclass of int, so the plain isinstance test turned every flag
        #  (deposit_full_bore, wax_risk, clip_warning, ...) into 0.0/1.0, while the
        #  summary.json beside it kept them as true/false. Two files, same keys,
        #  different types.
        #
        #  INTEGERS ARE NOW PRESERVED TOO, for the same reason. Coercing every number to
        #  float fixed the booleans but left the COUNTS floated: `fallbacks` and
        #  `clip_activations` read 0 and 47194 in summary.json and 0.0 and 47194.0 in
        #  key_metrics.json, from one run. A count is not a measurement, and two files
        #  from the same run disagreeing on the type of a shared key is exactly what this
        #  coercion was added to stop.
        def _as_json_number(v):
            if isinstance(v, (bool, np.bool_)):
                return bool(v)
            if isinstance(v, (int, np.integer)):
                return int(v)
            if isinstance(v, (float, np.floating)):
                return float(v)
            return v
        solver.dump_json({k: _as_json_number(v) for k, v in eng.items()
                          if not isinstance(v, (dict, list))}, fh)
    print(f"  -> {sv.results['steps']} steps, {sv.results['fallbacks']} fallbacks, "
          f"massErr {eng['mass_conservation_err']*100:.2f}%")
    return sv, eng


if __name__ == "__main__":
    base = "Deepwater Medium-Crude-Oil Subsea Tie-back (32 km, 10.75-in)"

    # ---- (A) as-operated steady: the HIGH slug + hydrate prediction + full stack ----
    print("\n=== (A) STEADY — as-operated (degraded insulation, no inhibitor) ===")
    out_st = OUT["steady"]
    case_st = build_case(base + " — normal production (as-operated)", "asoperated", 48.0)
    sv_st, eng_st = run_core(case_st, out_st)
    print("  cross-section reconstruction ...");        shct_crosssection.crosssection_outputs(sv_st, out_st)
    print("  compositional PVT report ...");            shct_compositional.compositional_report(sv_st, out_st)
    print("  compositional transport ...")
    shct_compositional_sim.simulate_composition(sv_st, out_st)
    print("  3-D field + VTK ...");                     shct_threed.threed_outputs(sv_st, out_st)
    print("  OpenFOAM coupling (case generation) ...")
    shct_openfoam.couple(sv_st, out_st, max_sections=3, run=shct_openfoam.openfoam_available())
    print("  published-data closure validation (friction, drift-flux, slug-freq, hydrate) ...")
    val_datadir = os.path.join(ROOT, "validation", "data")
    val = solver.validate_closures(outdir=out_st, datadir=val_datadir)
    with open(os.path.join(out_st, "validation_summary.json"), "w") as fh:
        solver.dump_json(val, fh)

    # ---- (B) unplanned shut-in: cooldown / no-touch time ----
    print("\n=== (B) SHUT-IN — unplanned shut-in cooldown (no-touch time) ===")
    out_si = OUT["shutin"]
    case_si = build_case(base + " — unplanned shut-in", "shutin", 24.0)
    sv_si, eng_si = run_core(case_si, out_si, riser=False)
    #  COUPLE EVERY SCENARIO, not just the as-operated one. All three folders ship an
    #  openfoam_cases/ directory, but couple() was called for the steady run alone, so the
    #  shut-in and mitigated sets were leftovers from whenever they were last written by
    #  hand -- 50 of their 56 files predated the run beside them, and they carried section
    #  states (gas density, holdup, deposit) from a solver version the v4.0.0 fluid-model
    #  correction has since moved. README says nothing under case/outputs_* is produced by
    #  hand; these two were the exception.
    print("  OpenFOAM coupling (case generation) ...")
    shct_openfoam.couple(sv_si, out_si, max_sections=3, run=shct_openfoam.openfoam_available())

    # ---- (C) engineered mitigation: insulation + MEG ----
    print("\n=== (C) MITIGATED — restored insulation + continuous MEG (design) ===")
    out_mt = OUT["mitigated"]
    case_mt = build_case(base + " — engineered fix (insulation + MEG)", "mitigated", 48.0)
    sv_mt, eng_mt = run_core(case_mt, out_mt, riser=False)
    print("  OpenFOAM coupling (case generation) ...")
    shct_openfoam.couple(sv_mt, out_mt, max_sections=3, run=shct_openfoam.openfoam_available())

    # ---- bespoke cross-scenario figures (saved with the steady run) ----
    print("\n=== bespoke comparison figures ===")
    hydrate_envelope_chart(sv_st, sv_si, out_st)
    mitigation_chart(eng_st, eng_mt, out_st)

    # ---- one compact comparison JSON the doc builder reads ----
    cmp = {"as_operated": {k: eng_st.get(k) for k in
              ["max_subcooling_C", "peak_deposit_mm", "deposit_full_bore", "P_plug",
               "time_to_plug_P50_h", "MEG_wt_pct", "MEG_Lph", "under_inhibited_km",
               "cooldown_to_hydrate_h", "U_eff_WmK"]},
           "mitigated": {k: eng_mt.get(k) for k in
              ["max_subcooling_C", "peak_deposit_mm", "deposit_full_bore", "P_plug",
               "time_to_plug_P50_h", "MEG_wt_pct", "MEG_Lph", "under_inhibited_km",
               "cooldown_to_hydrate_h", "U_eff_WmK"]},
           "shutin": {k: eng_si.get(k) for k in
              ["cooldown_to_hydrate_h", "cooldown_source", "max_subcooling_C"]}}
    with open(os.path.join(out_st, "scenario_comparison.json"), "w") as fh:
        solver.dump_json(cmp, fh)

    print("\nALL RUNS COMPLETE.")
