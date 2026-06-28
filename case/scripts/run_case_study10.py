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
#                             (U_eff ~ 2.4 W/m2K) + continuous MEG -> shows the
#                             model used as a DESIGN tool (risk removed).
#
#  The steady run also drives the full advanced stack: compositional Peng-Robinson
#  EOS PVT, cross-section / quasi-3-D reconstruction (+ VTK), compositional
#  transport, OpenFOAM (interFoam) coupling case generation, and hydrate-curve
#  validation against published experimental data.  Everything (CSV tables,
#  engineering deliverables, charts/curves, JSON, bespoke slug/hydrate/mitigation
#  figures, and the input-data deck) is written under 10/.
#
#  NOTE ON DATA PROVENANCE (honest framing): the field is a *representative
#  industrial archetype*.  Geometry, fluid and operating parameters are realistic
#  and self-consistent values typical of deepwater medium-crude-oil tie-backs in the
#  open literature; they are NOT proprietary operator data.  The PHYSICS and the
#  PREDICTIONS are produced by the real solver, and the hydrate thermodynamics are
#  validated against published experimental data (see §validation).
# =============================================================================
import os, sys, json, math
import numpy as np

from _paths import HERE, CASE, ROOT, OUT    # shared layout + no-black style (shct_style)
import solver
import shct_crosssection, shct_compositional, shct_compositional_sim, shct_threed, shct_openfoam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import asdict
import shct_style as S

# medium, non-black, non-dark palette (see shct_style.py)
NAVY, ACC, ORG, RED, GRN, TEAL = S.BLUE, "#1F8AC0", S.ORANGE, S.RED, S.GREEN, S.TEAL

# -----------------------------------------------------------------------------
#  Fluid: a MEDIUM CRUDE OIL — a ~30 deg API black oil (compositional makeup
#  -> Peng-Robinson flash).  Moderate-GOR with a substantial heavy C7+ tail
#  (~31 mol%): the dissolved/associated gas (C1 ~43 mol%) liberates along the
#  cold line and, with the long near-horizontal undulating step-out, drives
#  HYDRODYNAMIC (and terrain/severe-riser) SLUGGING; the 35 % water cut plus the
#  water-wet associated gas in the cold deepwater wall drives HYDRATE formation.
#  This crude is a textbook combination for BOTH slugging and hydrates.
# -----------------------------------------------------------------------------
CRUDE_OIL = {
    "N2": 0.004, "CO2": 0.020, "C1": 0.430, "C2": 0.075, "C3": 0.058,
    "iC4": 0.012, "nC4": 0.028, "iC5": 0.013, "nC5": 0.016, "C6": 0.030, "C7+": 0.314,
}

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
    f.water_cut = 0.35                    # 35 % water cut -> ample free water for hydrates
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
    o.q_liquid_insitu = 0.055             # in-situ liquid (oil + water) rate (m3/s)
    o.q_gas_insitu_inlet = 0.150          # in-situ associated-gas rate (m3/s) -> gassy, intermittent
    o.P_inlet_bar = 150.0
    o.T_inlet_C = 58.0
    o.T_seabed_C = 4.0

    n = c.numerics
    n.t_end_h = t_end_h
    n.n_cells = n_cells
    n.n_ensemble = n_ensemble
    n.n_snapshots = 80
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
    med = lambda A: np.nanmedian(A, 1)
    reg = med(r["regime"]); fsl = med(r["fslug"]); hold = med(r["alpha_l"])
    Lu = solver.slug_length(med(r["j"]), sv.case.pipeline.diameter_m, fsl)

    fig, ax = plt.subplots(3, 1, figsize=(8, 7.2), sharex=True)
    ax[0].fill_between(x, sv.z, sv.z.min() - 20, color=S.TAN, alpha=.55)
    ax[0].plot(x, sv.z, color=S.BROWN, lw=1.2); ax[0].set_ylabel("elevation (m)")
    ax[0].set_title("Slug-formation prediction — deepwater medium-crude-oil tie-back",
                    color=NAVY, fontweight="bold")
    ax[1].plot(x, hold, color=ACC, lw=1.8, label="liquid holdup α_l")
    ax[1].set_ylabel("holdup α_l"); ax[1].set_ylim(0, 1)
    sl = np.isin(np.round(reg), [2, 5])
    ax[1].fill_between(x, 0, 1, where=sl, color="#f6d6d2", alpha=.5,
                       transform=ax[1].get_xaxis_transform(), label="intermittent (slug/churn)")
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    lf, = ax[2].plot(x, fsl, color=ORG, lw=1.8, label="slug frequency f_slug (Hz)")
    ax[2].set_ylabel("f_slug (Hz)", color=ORG); ax[2].set_xlabel("distance along route (km)")
    a2 = ax[2].twinx(); ll, = a2.plot(x, Lu, color=NAVY, lw=1.4, ls="--", label="slug-unit length (m)")
    a2.set_ylabel("slug length (m)", color=NAVY)
    ax[2].legend(handles=[lf, ll], fontsize=8, loc="upper left", bbox_to_anchor=(1.13, 1.0),
                 borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/09_slug_prediction.png", dpi=155); plt.close(fig)


# -----------------------------------------------------------------------------
#  Bespoke chart 2 — severe-slugging RISER zoom (last ~6 km)
# -----------------------------------------------------------------------------
def riser_chart(sv, outdir):
    r = sv.results; x = sv.x / 1000.0; med = lambda A: np.nanmedian(A, 1)
    m = x >= (x.max() - 6.0)
    hold = med(r["alpha_l"]); reg = np.round(med(r["regime"]))
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 5.6), sharex=True)
    ax[0].fill_between(x[m], sv.z[m], sv.z[m].min() - 10, color=S.TAN, alpha=.55)
    ax[0].plot(x[m], sv.z[m], color=S.BROWN, lw=1.4); ax[0].set_ylabel("elevation (m)")
    ax[0].set_title("Severe-slugging screen — riser base & ascent (last 6 km)",
                    color=NAVY, fontweight="bold")
    ax[1].plot(x[m], hold[m], color=ACC, lw=2.0, label="liquid holdup α_l")
    sl = np.isin(reg, [2, 5])
    ax[1].fill_between(x[m], 0, 1, where=sl[m], color="#f6d6d2", alpha=.5,
                       transform=ax[1].get_xaxis_transform(), label="intermittent (slug/churn)")
    ax[1].set_ylabel("holdup α_l"); ax[1].set_ylim(0, 1); ax[1].set_xlabel("distance (km)")
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/10_riser_severe_slug.png", dpi=155); plt.close(fig)


# -----------------------------------------------------------------------------
#  Bespoke chart 3 — hydrate P-T envelope with BOTH the production and shut-in
#  trajectories overlaid (shows the line driven INTO the hydrate region)
# -----------------------------------------------------------------------------
def hydrate_envelope_chart(sv_op, sv_si, outdir):
    c = sv_op.case; med = lambda A: np.nanmedian(A, 1)
    rO, rS = sv_op.results, sv_si.results
    Pc = np.linspace(5.0, max(180.0, med(rO["p"]).max() * 1.1), 200)
    Tc = solver.hydrate_equilibrium_T(Pc, gas_sg=c.fluids.gas_sg,
                                      salinity_wt=c.fluids.salinity_wt, table=c.fluids.hyd_Teq_table)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.plot(Tc, Pc, color=RED, lw=2.4, label="hydrate equilibrium (T_eq)")
    ax.fill_betweenx(Pc, 0, Tc, color="#f6d6d2", alpha=.45, label="hydrate stability region")
    ax.plot(med(rO["T"]), med(rO["p"]), color=NAVY, lw=2.0, marker="o", ms=2.5,
            label="production trajectory (as-operated)")
    ax.plot(med(rS["T"]), med(rS["p"]), color=TEAL, lw=2.0, ls="--", marker="s", ms=2.5,
            label="shut-in trajectory")
    ax.set_xlim(2, 32); ax.set_ylim(0, max(180, med(rO["p"]).max() * 1.1))
    ax.set_xlabel("temperature (°C)"); ax.set_ylabel("pressure (bar)")
    ax.set_title("Hydrate-formation prediction — P–T trajectories vs envelope",
                 color=NAVY, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/11_hydrate_envelope.png", dpi=155); plt.close(fig)


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
    ax.set_title("Mitigation comparison — model used as a flow-assurance design tool",
                 color=NAVY, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/12_mitigation_comparison.png", dpi=155); plt.close(fig)


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
    solver._save_csv(f"{outdir}/input_data_deck.csv",
                     ["group", "parameter", "value", "units"],
                     [[g, pname, v, u] for g, pname, v, u in rows], all_str=True)
    comp_rows = [[k, f"{v:.4f}"] for k, v in CRUDE_OIL.items()]
    solver._save_csv(f"{outdir}/feed_composition.csv", ["component", "mol_fraction"],
                     comp_rows, all_str=True)
    with open(os.path.join(outdir, "case_config.json"), "w") as fh:
        json.dump(asdict(case), fh, indent=2, default=str)


# -----------------------------------------------------------------------------
#  Core: run a scenario, dump every standard output
# -----------------------------------------------------------------------------
def run_core(case, outdir, slug=True, riser=True):
    os.makedirs(outdir, exist_ok=True)
    save_input_deck(case, outdir)
    sv = solver.TransientSHCT(case); sv.run(verbose=True)
    eng = sv.engineering()
    solver.write_tables(sv, eng, outdir)
    solver.make_charts(sv, eng, outdir)
    if slug:
        slug_chart(sv, outdir)
    if riser:
        riser_chart(sv, outdir)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(eng, fh, indent=2, default=str)
    with open(os.path.join(outdir, "key_metrics.json"), "w") as fh:
        json.dump({k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v)
                   for k, v in eng.items() if not isinstance(v, (dict, list))}, fh,
                  indent=2, default=str)
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
    print("  compositional transport ...");             shct_compositional_sim.simulate_composition(sv_st, out_st)
    print("  3-D field + VTK ...");                     shct_threed.threed_outputs(sv_st, out_st)
    print("  OpenFOAM coupling (case generation) ...")
    shct_openfoam.couple(sv_st, out_st, max_sections=3, run=shct_openfoam.openfoam_available())
    print("  published-data closure validation (friction, drift-flux, slug-freq, hydrate) ...")
    val_datadir = os.path.join(ROOT, "validation", "data")
    val = solver.validate_closures(outdir=out_st, datadir=val_datadir)
    with open(os.path.join(out_st, "validation_summary.json"), "w") as fh:
        json.dump(val, fh, indent=2, default=str)

    # ---- (B) unplanned shut-in: cooldown / no-touch time ----
    print("\n=== (B) SHUT-IN — unplanned shut-in cooldown (no-touch time) ===")
    out_si = OUT["shutin"]
    case_si = build_case(base + " — unplanned shut-in", "shutin", 24.0)
    sv_si, eng_si = run_core(case_si, out_si, riser=False)

    # ---- (C) engineered mitigation: insulation + MEG ----
    print("\n=== (C) MITIGATED — restored insulation + continuous MEG (design) ===")
    out_mt = OUT["mitigated"]
    case_mt = build_case(base + " — engineered fix (insulation + MEG)", "mitigated", 48.0)
    sv_mt, eng_mt = run_core(case_mt, out_mt, riser=False)

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
        json.dump(cmp, fh, indent=2, default=str)

    print("\nALL RUNS COMPLETE.")
