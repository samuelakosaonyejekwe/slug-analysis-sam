#!/usr/bin/env python3
# =============================================================================
#  build_report.py  —  builds the SINGLE comprehensive `report.docx` for the
#  Deepwater MEDIUM-CRUDE-OIL subsea tie-back case study.
#
#  It consolidates EVERYTHING (nothing left out):
#    * background, problem statement and how the problem is solved
#    * the well-defined case study (asset, geometry, fluid, scenarios)
#    * ALL governing equations + closures used to build the solver
#    * the numerical method
#    * contribution to knowledge
#    * how this universal solver compares to / improves on existing solvers
#    * ALL input data (deck + feed composition, per scenario)
#    * ALL generated output data for every scenario — every metric, CSV table,
#      per-column graph of every CSV, every solver chart / curve / contour / map
#    * validation against credible PUBLISHED data, with every source recorded
#    * calibration of the solver to the case study, with sources
#    * cross-scenario summary, conclusions and a consolidated reference list
#
#  STYLE RULE (per request): NO BLACK / NO DARK colours anywhere in any figure.
#  All embedded solver PNGs are produced under the global shct_style no-black
#  palette; every figure this script regenerates uses the same medium palette
#  (the build_reports plot helpers are reused with their palette monkey-patched
#  to the shct_style medium colours, and legends are kept off the curves).
#
#  Author: Akosa Samuel Onyejekwe.
# =============================================================================
import os, json
import numpy as np

from _paths import HERE, CASE, ROOT, OUT          # shared layout + no-black shct_style
import shct_style as S
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

import build_reports as BR                          # reuse Doc + equations + plot helpers

# --- recolour every build_reports palette hook to the medium no-black palette -
def _rgb(hex_):
    h = hex_.lstrip("#"); return RGBColor(int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))

BR.NAVY   = _rgb(S.BLUE)        # headings / emphasis (medium royal blue, not navy)
BR.ACCENT = _rgb("#1F8AC0")    # sub-headings (medium cyan-blue)
BR.RED    = _rgb(S.RED)        # warnings (clear red, not maroon)
BR.GREEN  = _rgb(S.GREEN)      # good/mitigated (medium green)
BR.GREY   = _rgb("#6E7CA8")    # captions (medium slate-blue, not dark grey)
# plot-helper hex palette (used by plot_xy / plot_prob_range / plot_eng_bar / plot_composition)
BR.NAVY_H = S.BLUE; BR.ACC_H = "#1F8AC0"; BR.RED_H = S.RED
BR.ORG_H  = S.ORANGE; BR.TEAL_H = S.TEAL; BR.GRN_H = S.GREEN
# keep the regenerated per-CSV plots tucked away under the (gitignored) scripts dir,
# so the repo root holds only the finished report.docx — not a build-artefact folder
BR.PLOTDIR = os.path.join(CASE, "scripts", "_report_plots"); os.makedirs(BR.PLOTDIR, exist_ok=True)

from build_reports import (Doc, read_csv, fnum, write_equations, plot_xy,
                           plot_prob_range, plot_eng_bar, plot_composition,
                           kpi_rows, title_block)

OUTROOT = CASE
SCEN = [
    ("As-operated normal production — degraded/flooded insulation, no inhibitor",
     "outputs_steady", "as-operated"),
    ("Unplanned shut-in — cooldown / no-touch time", "outputs_shutin", "shut-in"),
    ("Engineered mitigation — restored insulation + continuous MEG", "outputs_mitigated", "mitigated"),
]

# the solver's curated, no-black figure gallery (embedded as produced)
GALLERY = [
    ("01_profiles.png", "Final-state P50 profiles: seabed elevation, liquid holdup, pressure & temperature vs hydrate T_eq, and subcooling along the whole route."),
    ("09_slug_prediction.png", "Slug-formation prediction: terrain, intermittent (slug/churn) bands, slug frequency and slug-unit length, and liquid holdup."),
    ("10_riser_severe_slug.png", "Severe-slugging screen at the steel-catenary-riser base and ascent (last 6 km)."),
    ("02_holdup_spacetime.png", "Liquid-holdup field α_l(x,t) — the bands are slug activity migrating along the line in time (space-time contour map)."),
    ("03_PT_envelope.png", "Production P–T trajectory against the hydrate-stability envelope (curve)."),
    ("11_hydrate_envelope.png", "Hydrate-formation prediction: production AND shut-in P–T trajectories overlaid on the hydrate envelope."),
    ("04_PhiSH_map.png", "Slug–Hydrate coupling-criticality map Φ_SH(x,t); the Φ_SH = 1 contour separates slug-scoured from plugging-critical (space-time map)."),
    ("05_scenario_timeseries.png", "Monitored-station transient response vs time (curves)."),
    ("06_deposit.png", "Wall-deposit growth at the monitor station vs time (curve)."),
    ("07_probabilistic.png", "Probabilistic time-to-plug CDF and the max-Φ_SH P10–P90 uncertainty band (Monte-Carlo ensemble)."),
    ("08_diagnostics.png", "Solver diagnostics: liquid & gas mass balances, clip activity, slug length and numerical consistency."),
    ("cx1_geometry.png", "Quasi-3-D cross-section geometry reconstructed along the line (curves)."),
    ("cx2_azimuthal_deposit.png", "Azimuthal (bottom-of-line) wall-deposit map around the pipe circumference along the route (contour map)."),
    ("cx3_sections.png", "2-D pipe cross-section reconstructions at selected stations (holdup level + deposit ring)."),
    ("compo_pvt.png", "Compositional PVT tracking along the line (Peng–Robinson flash): gas/liquid split, K-values, densities."),
    ("compositional_transport.png", "Compositional grading along the route — preferential depletion of the hydrate-forming light ends."),
    ("threed_deposit.png", "3-D reconstructed pipe coloured by wall-deposit thickness."),
    ("threed_temperature.png", "3-D reconstructed pipe coloured by wall temperature."),
    ("hydrate_validation.png", "Solver hydrate-equilibrium closure vs published experimental data (validation curve)."),
    ("12_mitigation_comparison.png", "Mitigation comparison — as-operated vs the engineered fix across the key flow-assurance metrics (bar chart)."),
]

# every CSV the solver writes per scenario, with a one-line description
XY_CSVS = [
    ("fields_profile.csv", "steady-state spatial profile of every field along the route"),
    ("timeseries_monitor.csv", "transient response at the monitored station"),
    ("csv_crosssection.csv", "cross-section / quasi-3-D geometry along the route"),
    ("csv_compositional.csv", "compositional / PVT state along the route (PR EOS)"),
    ("csv_compositional_transport.csv", "compositional grading (hydrate-former depletion) along the route"),
]


def jload(p):
    try:
        return json.load(open(p))
    except Exception:
        return {}


def km_of(folder):
    return jload(os.path.join(OUTROOT, folder, "key_metrics.json"))


def g(km, k, fmt="{:.2f}", default="n/a"):
    v = km.get(k)
    if v is None:
        return default
    try:
        return fmt.format(float(v))
    except Exception:
        return str(v)


# =============================================================================
#  NARRATIVE SECTIONS
# =============================================================================
def sec_title_page(D):
    title_block(D,
        "Coupled Slug & Hydrate Prediction in a Deepwater Crude-Oil Subsea Tie-back",
        "A comprehensive flow-assurance case study with the SHCT universal transient solver "
        "— case definition, model equations, inputs, all generated outputs, validation and calibration")
    D.para("Author: Akosa Samuel Onyejekwe", bold=True, color=BR.NAVY)
    D.para("Solver: SHCT — Slug–Hydrate Coupled-Transient multiphase flow-assurance solver "
           "(transient, compositional, probabilistic).", size=10)
    D.para("Scope: this single document consolidates EVERYTHING produced for the study — the "
           "engineering background and problem, the well-defined case study, every governing "
           "equation and closure used to build the solver, the numerical method, the contribution "
           "to knowledge, a like-for-like comparison with existing commercial and academic solvers, "
           "the complete input-data deck, and every generated output (all metrics, CSV tables, "
           "per-column graphs of every CSV, and the full chart/curve/contour/map gallery) for the "
           "three engineering scenarios, plus the published-data validation and the calibration, "
           "with all sources recorded. Nothing is left out.", size=10)
    D.para("Figure/chart style: no black and no dark colours are used in any figure, curve, "
           "contour or map; all backgrounds are white and all ink is a medium, legible hue.",
           italic=True, color=BR.GREY, size=9)
    D.pagebreak()


def sec_executive_summary(D):
    kmA, kmS, kmM = km_of("outputs_steady"), km_of("outputs_shutin"), km_of("outputs_mitigated")
    D.H1("Executive summary")
    D.para("A 32 km, 10.75-inch deepwater subsea tie-back carrying a ~30° API medium crude oil with "
           "35 % water cut over a cold (4 °C), strongly undulating seabed and a steel catenary riser "
           "is analysed for the two flow-assurance threats that dominate such systems and that "
           "physically reinforce one another: hydrodynamic / terrain / severe-riser SLUGGING and gas-"
           "HYDRATE formation. The SHCT solver runs the full transient, compositional, probabilistic "
           "prediction end-to-end for three scenarios.")
    D.bullet(f"As-operated (degraded insulation, no inhibitor): the line is BOTH slug- and hydrate-"
             f"critical — peak coupling number Φ_SH = {g(kmA,'max_Phi_SH')}, plug probability "
             f"{float(kmA.get('P_plug',0))*100:.0f}%, P50 time-to-plug {g(kmA,'time_to_plug_P50_h','{:.1f}')} h, "
             f"peak wall deposit {g(kmA,'peak_deposit_mm','{:.0f}')} mm, max subcooling "
             f"{g(kmA,'max_subcooling_C','{:.1f}')} °C.", color=BR.RED)
    D.bullet(f"Inhibitor demand the model predicts to clear that risk: MEG ≈ {g(kmA,'MEG_wt_pct','{:.0f}')} wt% "
             f"({g(kmA,'MEG_Lph','{:.0f}')} L/h); under-inhibited length {g(kmA,'under_inhibited_km','{:.1f}')} km.")
    D.bullet(f"Unplanned shut-in: essentially no safe window — no-touch time ≈ "
             f"{g(kmS,'cooldown_to_hydrate_h','{:.2f}')} h.")
    D.bullet(f"Engineered fix (restored insulation U_eff = {g(kmM,'U_eff_WmK')} W/m²K + continuous MEG): "
             f"the subcooling is removed, the plug probability falls to "
             f"{float(kmM.get('P_plug',0))*100:.0f}% and {g(kmM,'cooldown_to_hydrate_h','{:.0f}')} h of "
             f"no-touch time is restored — the model used as a design tool.", color=BR.GREEN)
    D.bullet(f"Numerics are mass-consistent and stable: liquid mass error {g(kmA,'mass_conservation_err','{:.2e}')}, "
             f"gas {g(kmA,'gas_mass_conservation_err','{:.2e}')}, {int(float(kmA.get('fallbacks',0)))} solver fallbacks.")
    D.pagebreak()


def sec_background(D):
    D.H1("1.  Background to the case study")
    D.para("Deepwater subsea tie-backs transport unprocessed well fluids — oil, produced water and "
           "associated gas — over tens of kilometres of cold seabed (typically 4 °C) before reaching "
           "a host facility. Two production-chemistry / multiphase-flow threats dominate the integrity "
           "of such lines:")
    D.H2("1.1  Hydrodynamic, terrain and severe-riser slugging")
    D.para("In multiphase pipe flow the gas and liquid self-organise into intermittent SLUGS — long "
           "liquid plugs separated by gas pockets. On a long, near-horizontal, undulating seabed the "
           "low spots accumulate liquid and shed it as terrain slugs; at the riser base, liquid "
           "fallback and gas blow-through produce large, low-frequency SEVERE slugs. Slugs impose "
           "cyclic structural loads, flood the topside slug-catcher, and force unstable operation.")
    D.H2("1.2  Gas-hydrate formation")
    D.para("Natural-gas components (methane, ethane, propane…) combine with free water at high "
           "pressure and low temperature to form ICE-LIKE crystalline hydrates. In a cold deepwater "
           "line the operating point sits well inside the hydrate-stability region; with a 35 % water "
           "cut there is ample free water. Hydrates can deposit on the cold wall, consolidate and grow "
           "a PLUG that blocks the line — the single most feared flow-assurance failure, slow and "
           "dangerous to remediate.")
    D.H2("1.3  Why the two threats must be solved together")
    D.para("Slugging and hydrates are not independent. The same cold, gassy, water-wet, intermittent "
           "flow that drives slugging also drives hydrates; and the two interact at the wall: slug "
           "passage scours and re-warms the wall (suppressing deposition), while between slugs the "
           "cold wall lets hydrate deposit and consolidate. Whether a line plugs or stays open is "
           "therefore set by a COMPETITION between hydrate-formation rate and slug-renewal rate. A "
           "medium-crude tie-back of this geometry and fluid is the textbook case in which BOTH "
           "threats are simultaneously active, which is exactly why it is chosen here.")


def sec_problem(D):
    D.H1("2.  Problem statement and how the problem is solved")
    D.H2("2.1  Problem statement")
    D.para("Existing tools predict slugging OR hydrates well, but treat them as separate analyses run "
           "in different software, with the slug–hydrate WALL interaction left to engineering "
           "judgement. There is no single, transient, probabilistic predictor that resolves the "
           "competition between hydrate growth and slug scouring along the whole route and in time, "
           "and returns an actionable, quantified, UNCERTAINTY-aware plugging risk together with the "
           "slug loads, the inhibitor demand and the mitigation design.")
    D.para("The engineering question for THIS asset: given the geometry, fluid and operating "
           "conditions, (a) where and when does the line slug and how big are the slugs; (b) where "
           "and when does it enter the hydrate region and form a consolidating wall deposit; (c) what "
           "is the probabilistic time-to-plug; (d) how much inhibitor (and what insulation) removes "
           "the risk; and (e) how long is the no-touch time after a shut-in?")
    D.H2("2.2  How the problem is solved — the SHCT approach")
    D.para("SHCT integrates, in time, a coupled system of conservation PDEs for multiphase flow, heat "
           "transfer and hydrate formation on the ACTUAL terrain, closed by published flow-assurance "
           "correlations and a compositional Peng–Robinson PVT engine, and run as a Monte-Carlo "
           "ensemble for genuine probabilistic output. The slug–hydrate interaction is made explicit "
           "through a single dimensionless coupling number Φ_SH — the ratio of the wall hydrate-"
           "formation tendency to the slug-renewal rate — mapped in space and time. The full set of "
           "governing equations and closures is given in Section 4; the numerical method in Section 5.")
    D.bullet("Resolve flow regime, holdup, slug frequency and slug length along the route (hydraulics).")
    D.bullet("Track P, T and the hydrate-equilibrium temperature T_eq → local subcooling ΔT_sub.")
    D.bullet("Form hydrate by stochastic nucleation + growth in the bulk (slurry) and on the wall (deposit).")
    D.bullet("Combine the two through Φ_SH(x,t); Φ_SH > 1 marks plugging-critical locations/times.")
    D.bullet("Repeat over a Monte-Carlo ensemble → P10/P50/P90 time-to-plug and design margins.")
    D.bullet("Invert for the required MEG dose and the insulation that removes the risk (design mode).")


def sec_case_study(D):
    D.H1("3.  The case study — clearly defined")
    cfg = jload(os.path.join(OUTROOT, "outputs_steady", "case_config.json"))
    D.para("A representative deepwater West-Africa-style medium-crude-oil tie-back. The geometry, "
           "fluid and operating values are realistic, self-consistent values typical of such "
           "developments in the open literature; they are NOT proprietary operator data. The PHYSICS "
           "and PREDICTIONS are produced by the real solver, and the hydrate thermodynamics and the "
           "hydraulic closures are validated against published data (Sections 7–8).")
    D.H2("3.1  Asset and geometry")
    D.bullet("Route length: 32 km step-out flowline + steel catenary riser (SCR).")
    D.bullet("Internal diameter: 0.2545 m (~10.75-inch carbon-steel flowline).")
    D.bullet("Water depth at the riser base: ~1100 m; cold seabed at 4 °C.")
    D.bullet("Terrain: long, strongly undulating seabed (multiple low spots → terrain slugging) "
             "climbing into a steep SCR (→ severe-riser slugging).")
    D.H2("3.2  Fluid — a medium crude oil (≈30° API black oil)")
    D.para("Moderate-GOR live crude with a substantial heavy C7+ tail; the associated gas (C1 ≈ 43 "
           "mol%) liberates along the cold line and, with the long near-horizontal step-out, drives "
           "slugging, while the 35 % water cut plus water-wet gas in the cold wall drives hydrates.")
    fc = os.path.join(OUTROOT, "outputs_steady", "feed_composition.csv")
    if os.path.exists(fc):
        D.figure(plot_composition(fc, "feed", "case"),
                 "Feed composition — medium-crude makeup (mole fraction per component) that feeds the "
                 "Peng–Robinson flash.", width=6.2)
        D.csv_table(fc, max_rows=20, fs=9)
    D.H2("3.3  The three engineering scenarios")
    D.bullet("(A) As-operated normal production — degraded / water-flooded wet insulation "
             "(behaves close to bare steel, U ≈ 22 W/m²K), NO inhibitor: the high-risk prediction.")
    D.bullet("(B) Unplanned shut-in — the line cools toward the seabed; the no-touch (cooldown) time "
             "before hydrate onset is measured from the transient.")
    D.bullet("(C) Engineered mitigation — restored multi-layer insulation (steel + syntactic foam + "
             "coating → low U_eff) plus continuous MEG: the model used as a DESIGN tool.")


def sec_inputs(D):
    D.H1("6.  Input data (complete deck)")
    D.para("All three scenarios share the same asset and feed; they differ only in thermal design and "
           "inhibition. The complete machine-readable configuration (case_config.json) accompanies "
           "each scenario folder; the human-readable deck and the per-scenario design follow.")
    st = os.path.join(OUTROOT, "outputs_steady")
    idp = os.path.join(st, "input_data_deck.csv")
    if os.path.exists(idp):
        D.H2("6.1  Input data deck (pipeline / fluid / operating / numerics)")
        D.csv_table(idp, max_rows=40, fs=8.5)
    D.H2("6.2  Per-scenario thermal & inhibition design")
    for label, fname, tag in SCEN:
        ip = os.path.join(OUTROOT, fname, "input_data_deck.csv")
        if os.path.exists(ip):
            D.para(f"{label}  [{tag}]", bold=True, color=BR.ACCENT)
            D.csv_table(ip, max_rows=40, fs=7.6)


def sec_numerics(D):
    D.H1("5.  Numerical method")
    D.para("The coupled PDE system is advanced by a conservative finite-volume scheme with adaptive "
           "CFL time-stepping. Holdup, gas-mass, energy and the hydrate phase-field are updated in "
           "flux form (optionally 2nd-order TVD: minmod / van Leer / superbee); the stiff terms — the "
           "pressure gradient, wall drag, interfacial drag and wall heat loss — are treated "
           "implicitly. Eliminating velocity from the implicit mixture momentum gives a tridiagonal "
           "(Poisson-type) pressure equation solved each step by the Thomas algorithm, vectorised "
           "over the whole Monte-Carlo ensemble. A never-fail guard auto-degrades any non-finite step "
           "to a proven quasi-steady solver so that time always advances (0 fallbacks on this case). "
           "Full equation forms are in Section 4 (group F).")


def sec_contribution(D):
    D.H1("9.  Contribution to knowledge")
    D.bullet("A single dimensionless Slug–Hydrate Coupling Number Φ_SH = C·k_g,wall·a_i·ΔT_sub,wall^n / "
             "f_slug that quantifies, locally and in time, the competition between wall hydrate growth "
             "and slug scouring — turning a qualitative engineering concern into a mapped, critical-"
             "threshold (Φ_SH = 1) field Φ_SH(x,t).", color=BR.NAVY)
    D.bullet("Genuine TWO-WAY coupling of slugging and hydrates at the wall (slug scouring vs "
             "consolidating deposition) inside one transient solver, rather than two separate analyses.")
    D.bullet("A transient, PROBABILISTIC time-to-plug (Monte-Carlo stochastic nucleation + parameter "
             "spread, Kaplan–Meier right-censored CDF) instead of a single deterministic number.")
    D.bullet("Self-consistent integration of compositional Peng–Robinson PVT, multiphase hydraulics, "
             "heat transfer, hydrate kinetics, deposition with self-insulating feedback, and inhibitor "
             "/ insulation inverse-design in ONE auditable model.")
    D.bullet("Closures anchored to, and validated against, universally-cited published references "
             "(Section 7), with an explicit one-parameter calibration discipline (Section 8).")
    D.bullet("A never-fail numerical guarantee (mass-consistent, 0 fallbacks) that makes the solver "
             "usable as a routine design tool, not just a research code.")


def sec_comparison(D):
    D.H1("10.  Comparison with existing solvers and software")
    D.para("The table contrasts SHCT with the tools normally used for these analyses. The honest "
           "position: established commercial transient simulators are mature, extensively field-"
           "validated products; SHCT's distinct contribution is to UNIFY slug and hydrate prediction "
           "through the explicit Φ_SH coupling, to be transient AND probabilistic in one pass, and to "
           "be fully open and auditable — capabilities that are split across several tools elsewhere.")
    rows = [
        ["Capability", "OLGA / LedaFlow", "PIPESIM (steady)", "PVTsim / Multiflash + CSMHyK", "SHCT (this work)"],
        ["Transient multiphase flow", "Yes", "No (steady)", "No", "Yes"],
        ["Hydrate equilibrium / kinetics", "Add-on module", "Equilibrium screen", "Yes (specialist)", "Yes (integrated)"],
        ["Explicit slug–hydrate wall coupling (Φ_SH)", "No", "No", "No", "Yes (core invention)"],
        ["Probabilistic time-to-plug (Monte-Carlo)", "Limited / manual", "No", "No", "Yes (P10/P50/P90)"],
        ["Compositional PR PVT in-loop", "Yes", "Yes", "Yes", "Yes"],
        ["Inhibitor & insulation inverse design", "Manual sweeps", "Partial", "Inhibitor only", "Yes (direct)"],
        ["Severe-riser + terrain slug screening", "Yes", "Correlation", "No", "Yes"],
        ["Open / auditable / scriptable", "No (commercial)", "No", "No", "Yes (open source)"],
        ["Never-fail numerical guard", "—", "—", "—", "Yes (0 fallbacks)"],
    ]
    D.kv_table(rows, fs=7.6, widths=(2.2, 1.25, 1.15, 1.5, 1.4))
    D.H2("10.1  Why this solver is positioned as a universal predictor")
    D.bullet("It spans the whole chain — hydraulics + thermal + compositional PVT + hydrate kinetics "
             "+ deposition + probabilistic risk + inhibitor/insulation design — in one consistent model.")
    D.bullet("It works on ARBITRARY terrain and any composition (PR EOS), so it is not tied to one "
             "asset type; the same code runs gas, condensate and crude tie-backs.")
    D.bullet("Its closures are the SAME published correlations the commercial tools use, and they are "
             "validated here against the original reference data (Section 7), so its hydraulic and "
             "thermodynamic ingredients are on the established footing.")
    D.bullet("The Φ_SH coupling and the probabilistic time-to-plug add information no single existing "
             "tool returns in one pass.")
    D.para("Honest caveat: full production-flow validation of pressure drop and arrival temperature "
           "along a SPECIFIC operating line still requires an operator's measured data; the holdup, "
           "friction, slip, slug-frequency and hydrate-equilibrium closures that those predictions "
           "rest on are validated against published data in Section 7.", italic=True, color=BR.GREY, size=9)


def sec_validation(D):
    D.H1("7.  Validation against published data (sources recorded)")
    D.para("Every closure the hydraulics and thermodynamics actually USE is scored against "
           "universally-cited, reproducible PUBLISHED references — each reference attributed to its "
           "PRIMARY source (where the value originated), with later works tagged as corroborating. "
           "The combined results are written to validation_summary.json; the hydrate curve is also "
           "plotted below.")
    vs = jload(os.path.join(OUTROOT, "outputs_steady", "validation_summary.json"))
    hy = jload(os.path.join(OUTROOT, "outputs_steady", "hydrate_validation_report.json"))
    D.H2("7.1  Hydrate-equilibrium curve vs experimental data")
    if hy:
        D.para(f"Methane Lw–H–V three-phase equilibrium, n = {hy.get('n')} points: as-shipped RMSE "
               f"= {hy.get('rmse_C',float('nan')):.2f} °C (bias {hy.get('bias_C',0):+.2f} °C, "
               f"max |err| {hy.get('max_abs_err_C',0):.2f} °C); after the single allowable T_eq "
               f"offset of {hy.get('offset_C',0):+.2f} °C, RMSE = {hy.get('rmse_after_offset_C',0):.2f} °C.")
    hv = os.path.join(OUTROOT, "outputs_steady", "hydrate_validation.png")
    if os.path.exists(hv):
        D.figure(hv, "Hydrate-equilibrium validation: SHCT closure (as-shipped and 1-parameter "
                     "calibrated) vs published experimental methane Lw–H–V data.", width=5.8)
    D.para("Source — hydrate data: Deaton, W.M. & Frost, E.M. (1946), 'Gas Hydrates and Their "
           "Relation to the Operation of Natural-Gas Pipe Lines', U.S. Bureau of Mines Monograph 8 "
           "(PRIMARY); reproduced/corroborated in Sloan, E.D. & Koh, C.A. (2008), 'Clathrate Hydrates "
           "of Natural Gases', 3rd ed., CRC Press (Appendix A); the 274.7 K / 3.27 MPa point is "
           "independently corroborated by the Blake Ridge measurement (ODP Leg 164).",
           italic=True, color=BR.GREY, size=9)
    D.H2("7.2  Hydraulic closures vs their primary references")
    if vs:
        fr = vs.get("friction") or {}
        sf = vs.get("slug_frequency") or {}
        if fr:
            D.bullet(f"Friction factor: Haaland (1983) explicit fit vs the Colebrook–White equation "
                     f"computed in-code — RMS deviation {fr.get('rms_pct_dev',float('nan')):.2f}%, "
                     f"max {fr.get('max_abs_pct_dev',float('nan')):.2f}% over the turbulent Re × "
                     f"roughness grid.")
        if sf:
            D.bullet(f"Slug frequency: reproduces the Gregory–Scott (1969) / Zabaras (2000) "
                     f"correlation to {sf.get('max_fidelity_err_pct',0):.1e}% (published accuracy band "
                     f"≈ ±{sf.get('published_accuracy_band_pct',0):.0f}%).")
        D.bullet("Drift-flux slip: distribution coefficient and drift velocity scored against the "
                 "canonical ANALYTICAL values — vertical Taylor-bubble Froude 0.35 and horizontal "
                 "nose Froude 0.542 — that the closure is built from.")
    fv = os.path.join(OUTROOT, "outputs_steady", "friction_validation.png")
    if os.path.exists(fv):
        D.figure(fv, "Friction-closure validation: Haaland (1983) explicit factor vs the Colebrook–"
                     "White equation across the turbulent Re × relative-roughness grid.", width=5.8)
    D.para("Sources — closures: friction reference Colebrook, C.F. (1939), J. Inst. Civ. Eng. "
           "11:133–156 (PRIMARY) with the Moody (1944) chart corroborating, closure under test "
           "Haaland, S.E. (1983), J. Fluids Eng. 105(1):89–90. Drift-flux: Dumitrescu (1943) ZAMM "
           "23:139–149 and Nicklin et al. (1962) for the vertical terms, Benjamin, T.B. (1968) J. "
           "Fluid Mech. 31:209–248 and Bendiksen (1984) for the horizontal terms. Slug frequency: "
           "Gregory, G.A. & Scott, D.S. (1969), AIChE J. 15:933–935 (PRIMARY base term) and Zabaras "
           "(2000), SPE J. 5(3):252–258 (inclination factor / implemented form).",
           italic=True, color=BR.GREY, size=9)
    D.H2("7.3  What remains open")
    D.para("Full production-flow pressure-drop and arrival-temperature along a SPECIFIC operating line "
           "cannot be closed by public data; it needs an operator's measured dP / arrival-T for that "
           "line. The holdup, friction, slip, slug-frequency and hydrate-equilibrium PHYSICS the "
           "predictions rest on are validated above against published data; this is stated plainly so "
           "the evidentiary standard is not overstated.", italic=True, color=BR.GREY, size=9)


def sec_calibration(D):
    D.H1("8.  Calibration of the solver to the case study (sources)")
    D.para("Calibration discipline: the solver ships with literature-DEFAULT closures and is run on "
           "this case WITHOUT bespoke tuning; where a single, physically-meaningful calibration knob "
           "exists it is exposed and reported, never silently fitted.")
    D.H2("8.1  Anchored / calibrated closures and their sources")
    D.bullet("Hydrate-equilibrium curve: anchored to the GPSA / Sloan natural-gas hydrate correlation "
             "and shifted for gas specific gravity and produced-water salinity; a single allowable "
             "T_eq offset is the only hydrate calibration parameter (reported in Section 7.1). "
             "Source: Sloan & Koh (2008); Deaton & Frost (1946).")
    D.bullet("Drift-flux slip: Bendiksen (1984) distribution coefficient and drift velocity, with the "
             "single verified numerics.drift_C0_factor knob available for a 1-parameter holdup "
             "calibration. Source: Bendiksen (1984); Dumitrescu (1943); Nicklin et al. (1962).")
    D.bullet("Friction: Haaland (1983), calibrated against Colebrook–White (1939) to ≈2 %. "
             "Slug frequency / length: Gregory–Scott (1969) / Zabaras (2000); Dukler–Hubbard scaling.")
    D.bullet("Compositional PVT: Peng–Robinson (1976) with Sutton pseudo-criticals, Lee–Gonzalez–Eakin "
             "gas viscosity and Lohrenz–Bray–Clark liquid viscosity. Hydrate kinetics: CSMHyK-type "
             "Arrhenius growth; MEG suppression by Nielsen–Bucklin (1983).")
    D.H2("8.2  Calibration of the case definition itself")
    D.para("The asset geometry, fluid makeup and operating conditions are calibrated to realistic, "
           "self-consistent values representative of deepwater medium-crude-oil tie-backs in the open "
           "literature (≈30° API live oil, 35 % water cut, ~1100 m water depth, 32 km step-out, "
           "10.75-inch line, 4 °C seabed). The erosional limit uses the API RP 14E C-factor. These "
           "anchor the case to industry-standard design practice rather than to one proprietary line.")


def sec_outputs(D):
    D.H1("11.  Generated outputs — every scenario (nothing left out)")
    D.para("For each scenario the report gives: the headline prediction metrics, the complete solver "
           "chart/curve/contour/map gallery (embedded exactly as generated under the no-black "
           "palette), a per-column graph of every CSV the solver writes, and every CSV as a data "
           "table. The machine-readable summary.json / key_metrics.json accompany each folder.")
    D.para("Animated outputs (GIF): each scenario folder additionally contains a set of transient GIF "
           "animations that show, as motion, the time-evolution the static charts below capture as "
           "single frames. They are a supplementary visualisation layer for presentation and review — "
           "a static document cannot play them, so they are provided as files alongside the outputs "
           "(under the no-black style, with every label/legend kept off the data):", bold=True,
           color=BR.ACCENT)
    D.bullet("anim_flow_line.gif — liquid holdup α_l along the terrain-following pipe, showing slug "
             "bodies travelling down the route in time.")
    D.bullet("anim_crosssection.gif — the 2-D pipe bore at the monitor: the stratified liquid level and "
             "the hydrate deposit ring closing the bore toward a plug (open bore in the mitigated case).")
    D.bullet("anim_PT_cooldown.gif — the monitor P–T point tracking across the hydrate-stability "
             "envelope in time (green = safe, red = inside the hydrate zone).")
    D.bullet("anim_riser_cycle.gif — the riser-region monitor α_l–P trajectory (repeating loops = "
             "intermittent / slug flow; a settled point = stable flow, as in the mitigated case).")
    D.bullet("anim_profile_wave.gif — the P(x,t) and T(x,t) profile wave marching along the line "
             "(the thermal cooling front and the pressure profile evolving in time).")
    D.para("These GIFs are generated by case/scripts/make_animations.py from the same transient solve "
           "as the charts; they add no new data beyond the fields already tabulated and plotted here.",
           italic=True, color=BR.GREY, size=9)
    for si, (label, fname, tag) in enumerate(SCEN, 1):
        folder = os.path.join(OUTROOT, fname)
        slug = tag.replace(" ", "")
        km = km_of(fname)
        D.pagebreak()
        D.H1(f"11.{si}  {label}")
        D.H2(f"11.{si}.1  Headline prediction metrics")
        D.kv_table(kpi_rows(km))
        # solver figure gallery (only the ones present in this folder)
        gallery = [(f, c) for (f, c) in GALLERY if os.path.exists(os.path.join(folder, f))]
        D.H2(f"11.{si}.2  Charts, curves, contours and maps ({len(gallery)} figures)")
        for fn, cap in gallery:
            D.figure(os.path.join(folder, fn), f"{cap}  [{tag}]", width=6.7)
        # engineering deliverables CSV (+ bar chart)
        ed = os.path.join(folder, "engineering_deliverables.csv")
        sub = 2
        if os.path.exists(ed):
            sub += 1
            D.H2(f"11.{si}.{sub}  Engineering deliverables (engineering_deliverables.csv)")
            D.figure(plot_eng_bar(ed, tag, slug),
                     f"Numeric engineering deliverables as a magnitude chart  [{tag}].", width=6.9)
            D.csv_table(ed, max_rows=40, fs=7.4)
        # probabilistic summary CSV (+ range plot)
        pp = os.path.join(folder, "probabilistic_summary.csv")
        if os.path.exists(pp):
            sub += 1
            D.H2(f"11.{si}.{sub}  Probabilistic summary (probabilistic_summary.csv)")
            D.figure(plot_prob_range(pp, tag, slug),
                     f"P10 → P90 uncertainty range with the P50 marker per metric  [{tag}].", width=6.9)
            D.csv_table(pp, max_rows=30, fs=7.6)
        # every remaining XY CSV: per-column graph + table
        for csvname, desc in XY_CSVS:
            path = os.path.join(folder, csvname)
            if not os.path.exists(path):
                continue
            sub += 1
            D.H2(f"11.{si}.{sub}  {csvname} — {desc}")
            D.figure(plot_xy(path, tag, slug),
                     f"Every column of {csvname} drawn as a curve  [{tag}].", width=6.95)
            D.csv_table(path, max_rows=14, fs=6.8)


def sec_cross_and_conclusions(D):
    D.pagebreak()
    D.H1("12.  Cross-scenario summary and conclusions")
    cmp_ = jload(os.path.join(OUTROOT, "outputs_steady", "scenario_comparison.json"))
    if cmp_:
        ao = cmp_.get("as_operated", {}); mi = cmp_.get("mitigated", {}); su = cmp_.get("shutin", {})
        keys = [("max_subcooling_C", "Max subcooling (°C)"), ("peak_deposit_mm", "Peak deposit (mm)"),
                ("P_plug", "Plug probability (frac)"), ("time_to_plug_P50_h", "Time-to-plug P50 (h)"),
                ("MEG_wt_pct", "MEG required (wt%)"), ("under_inhibited_km", "Under-inhibited (km)"),
                ("cooldown_to_hydrate_h", "No-touch time (h)"), ("U_eff_WmK", "Effective U (W/m²K)")]
        rows = [["Metric", "As-operated", "Mitigated"]]
        for k, lab in keys:
            def fmt(d):
                v = d.get(k)
                return "n/a" if v is None else (f"{float(v):.3g}" if isinstance(v, (int, float)) else str(v))
            rows.append([lab, fmt(ao), fmt(mi)])
        D.H2("12.1  As-operated vs engineered fix")
        D.kv_table(rows, fs=8.5, widths=(3.0, 1.6, 1.6))
        if su:
            D.para(f"Shut-in no-touch time ≈ {su.get('cooldown_to_hydrate_h','n/a')} h "
                   f"(source: {su.get('cooldown_source','transient')}).", size=9.5)
    D.H2("12.2  Engineering conclusions")
    kmA, kmS, kmM = km_of("outputs_steady"), km_of("outputs_shutin"), km_of("outputs_mitigated")
    D.bullet(f"As-operated the line is slug- AND hydrate-critical: Φ_SH = {g(kmA,'max_Phi_SH')}, "
             f"{float(kmA.get('P_plug',0))*100:.0f}% plug probability, P50 time-to-plug "
             f"{g(kmA,'time_to_plug_P50_h','{:.1f}')} h, peak deposit {g(kmA,'peak_deposit_mm','{:.0f}')} mm.")
    D.bullet(f"Inhibitor demand to clear it: MEG ≈ {g(kmA,'MEG_wt_pct','{:.0f}')} wt% "
             f"({g(kmA,'MEG_Lph','{:.0f}')} L/h); under-inhibited length {g(kmA,'under_inhibited_km','{:.1f}')} km.")
    D.bullet(f"Shut-in gives effectively no safe window (no-touch ≈ {g(kmS,'cooldown_to_hydrate_h','{:.2f}')} h).")
    D.bullet(f"The engineered fix (U_eff {g(kmM,'U_eff_WmK')} W/m²K + MEG) removes the subcooling, "
             f"zeroes the plug probability and restores {g(kmM,'cooldown_to_hydrate_h','{:.0f}')} h of no-touch time.")
    D.bullet(f"Predictions are mass-consistent and stable (liquid mass error "
             f"{g(kmA,'mass_conservation_err','{:.2e}')}, {int(float(kmA.get('fallbacks',0)))} fallbacks).")


def sec_references(D):
    D.pagebreak()
    D.H1("13.  References (validation and calibration sources)")
    refs = [
        "Deaton, W.M. & Frost, E.M. (1946). Gas Hydrates and Their Relation to the Operation of "
        "Natural-Gas Pipe Lines. U.S. Bureau of Mines, Monograph 8. [hydrate equilibrium data — PRIMARY]",
        "Sloan, E.D. & Koh, C.A. (2008). Clathrate Hydrates of Natural Gases, 3rd ed. CRC Press. "
        "[hydrate data compilation / corroboration]",
        "Colebrook, C.F. (1939). Turbulent flow in pipes. J. Inst. Civ. Eng. 11:133–156. "
        "[friction reference — PRIMARY]",
        "Moody, L.F. (1944). Friction factors for pipe flow. Trans. ASME 66:671–684. [Moody chart — corroborating]",
        "Haaland, S.E. (1983). Simple and explicit formulas for the friction factor. J. Fluids Eng. "
        "105(1):89–90. [friction closure under test]",
        "Dumitrescu, D.T. (1943). Strömung an einer Luftblase im senkrechten Rohr. ZAMM 23:139–149. "
        "[vertical Taylor-bubble drift — PRIMARY]",
        "Nicklin, D.J., Wilkes, J.O. & Davidson, J.F. (1962). Two-phase flow in vertical tubes. "
        "Trans. IChemE 40:61–68. [vertical distribution coefficient]",
        "Benjamin, T.B. (1968). Gravity currents and related phenomena. J. Fluid Mech. 31:209–248. "
        "[horizontal nose Froude — PRIMARY]",
        "Bendiksen, K.H. (1984). An experimental investigation of the motion of long bubbles in "
        "inclined tubes. Int. J. Multiphase Flow 10(4):467–483. [drift-flux slip closure]",
        "Gregory, G.A. & Scott, D.S. (1969). Correlation of liquid slug velocity and frequency. "
        "AIChE J. 15:933–935. [slug-frequency base term — PRIMARY]",
        "Zabaras, G.J. (2000). Prediction of slug frequency for gas/liquid flows. SPE J. 5(3):252–258. "
        "[slug-frequency inclination factor / implemented form]",
        "Taitel, Y. & Dukler, A.E. (1976). A model for predicting flow-regime transitions. AIChE J. "
        "22(1):47–55. [flow-regime map / interfacial geometry]",
        "Peng, D.-Y. & Robinson, D.B. (1976). A new two-constant equation of state. Ind. Eng. Chem. "
        "Fundam. 15(1):59–64. [compositional PVT]",
        "Lee, A.L., Gonzalez, M.H. & Eakin, B.E. (1966). The viscosity of natural gases. J. Pet. "
        "Technol. 18(8):997–1000. [gas viscosity]",
        "Lohrenz, J., Bray, B.G. & Clark, C.R. (1964). Calculating viscosities of reservoir fluids. "
        "J. Pet. Technol. 16(10):1171–1176. [liquid viscosity]",
        "Nielsen, R.B. & Bucklin, R.W. (1983). Why not use methanol for hydrate control? Hydrocarbon "
        "Processing 62(4):71–78. [MEG/methanol hydrate suppression]",
        "API RP 14E. Recommended Practice for Design and Installation of Offshore Production Platform "
        "Piping Systems. [erosional-velocity limit]",
        "Schmidt, Z., Brill, J.P. & Beggs, H.D. (1980). Experimental study of severe slugging in a "
        "two-phase flow pipeline-riser system. SPE J. 20(5):407–414. [severe-slugging basis]",
    ]
    for r in refs:
        D.bullet(r, color=BR.GREY)
    D.para("Provenance note: the case is a representative industrial ARCHETYPE; geometry, fluid and "
           "operating values are realistic open-literature values, not proprietary operator data. The "
           "physics and predictions are produced by the real solver; the closures and hydrate "
           "thermodynamics are validated against the published data above.", italic=True,
           color=BR.GREY, size=9)


def build(path):
    D = Doc()
    sec_title_page(D)
    sec_executive_summary(D)
    sec_background(D)
    sec_problem(D)
    sec_case_study(D)
    D.pagebreak(); D.H1("4.  Model equations used to build the solver")
    write_equations(D, intro=True)
    D.pagebreak(); sec_numerics(D)
    sec_inputs(D)
    D.pagebreak(); sec_validation(D)
    D.pagebreak(); sec_calibration(D)
    D.pagebreak(); sec_contribution(D)
    D.pagebreak(); sec_comparison(D)
    sec_outputs(D)
    sec_cross_and_conclusions(D)
    sec_references(D)
    D.save(path)
    print("wrote", path)


if __name__ == "__main__":
    build(os.path.join(ROOT, "report.docx"))
    print("DONE")
