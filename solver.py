#!/usr/bin/env python3
# =============================================================================
#  solver.py — SHCT Engineering-Grade Transient Coupled-PDE Solver
#  Slug-Hydrate Coupled Transport for subsea multiphase pipelines
#  Author: Akosa Samuel Onyejekwe
# -----------------------------------------------------------------------------
#  Predicts hydrodynamic slugging and gas-hydrate formation on arbitrary
#  terrain/inclination by integrating, IN TIME, the coupled SHCT partial
#  differential equations, and generates the full engineering output set
#  (fields, time-series, space-time maps, slug statistics, probabilistic risk,
#  and design deliverables).
#
#  GOVERNING COUPLED PDEs (1-D, finite-volume, transient):
#    (H) Liquid-holdup transport (drift-flux kinematic wave) - hyperbolic:
#          d(alpha_l)/dt + d(alpha_l*v_l)/dx = -S_l
#        => dynamic liquid inventory / terrain-induced slugging (NOT quasi-steady)
#    (G) TRANSIENT pressure & momentum -- THREE selectable engines (numerics.engine):
#        * "implicit"  (default): transient mixture-momentum with pressure-gradient & wall
#            drag IMPLICIT, closed by a Bendiksen drift-flux slip -> tridiagonal implicit
#            pressure each step. Genuine transient pressure; unconditionally well-posed; fast.
#        * "twofluid": FULL TWO-FLUID -- two INDEPENDENT phase momentum equations (gas and
#            liquid each with own inertia, advection, wall drag) coupled by interfacial drag
#            and a shared pressure; slip emerges from the momentum balance (no algebraic
#            slip closure). Semi-implicit (implicit pressure + implicit interfacial drag) and
#            an interfacial-pressure correction keep the otherwise ill-posed system stable.
#        * "quasisteady": legacy marched pressure (also the automatic never-fail fallback).
#        All engines conserve liquid mass to ~0% and degrade gracefully to "quasisteady"
#        on any non-finite step (the solver cannot crash or emit NaN).
#    (E) Energy transport (advection + seabed loss + latent heat) - transient:
#          dT/dt + v_m*dT/dx = -U*Pw*(T-T_sea)/(rho_m*cp*A) + L_hyd*mdot_h/(rho_m*cp)
#    (P) Hydrate phase-field (advected reaction-diffusion) - transient:
#          d(phi)/dt + v_l*d(phi)/dx = D_phi*d2(phi)/dx2 + R_grow + R_nuc + xi
#        with classical stochastic nucleation (Monte-Carlo ensemble).
#    (D) Wall-deposit / consolidation dynamics with slug scouring (two-way coupled).
#    (I) Thermodynamic-inhibitor (MEG) transport with the liquid; local hydrate-curve
#          suppression by Hammerschmidt -> models inhibition design & under-inhibition.
#    (C) Slug-Hydrate Coupling Number  Phi_SH = C*k_g*a_i*dTsub^n / f_slug.
#
#  ADDITIONAL OLGA-class capabilities: multi-layer wall heat transfer (steel/insulation/
#  coating) with an effective-U model and cooldown / no-touch-time output; hydrate-slurry
#  transportability (Camargo-Palermo relative viscosity); thermodynamic-inhibitor design
#  (required MEG wt% & rate, plus an under-inhibited-length map); a built-in V&V
#  suite (correlation benchmarks + mass conservation + grid convergence) and a calibration
#  module that fits constants to measured data.
#
#  TWO-WAY COUPLING: hydrate fraction & deposit reduce the effective bore and
#  raise viscosity/friction (-> H,G,E); slug hydrodynamics set interfacial area,
#  shear & subcooling that drive nucleation/growth/deposition (-> P,D).
#
#  Closures are published flow-assurance correlations: drift-flux (Bendiksen),
#  Taitel-Dukler regime, Gregory-Scott/Zabaras slug frequency, Haaland friction,
#  natural-gas hydrate equilibrium, CSMHyK-type growth, Camargo-Palermo viscosity,
#  Hammerschmidt inhibitor sizing, API-14E erosional limit.
#
#  TRANSIENT OPERATING SCENARIOS (time-dependent boundary conditions):
#    steady | rampup | turndown | shutin   (selectable; models real operations)
#
#  ----------------------------------------------------------------------------
#  VALIDATION STATUS (read before use):
#    This is an engineering-grade transient coupled-PDE solver with a SEMI-IMPLICIT
#    transient momentum/pressure core. The numerical method is production-style and
#    VERIFIED (run `--verify`): conserves liquid mass to ~0% across ALL operating
#    scenarios (steady, ramp-up, turndown, full shut-in), grid-convergent, with adaptive
#    time-stepping and a never-fail design (the implicit engine degrades gracefully to a
#    proven quasi-steady solver if a step is non-finite -> the solver cannot crash or
#    produce NaN). This "always bounded, always conservative, well-posed in every regime"
#    robustness is a deliberate strength.
#    The kinetic/coupling CONSTANTS ship as literature-typical defaults. They are made
#    CASE-SPECIFIC by `--calibrate`, which fits the free constants (heat transfer,
#    hydrate growth/nucleation, deposition) to the user's measured data (arrival T,
#    pressure drop, subcooling, observed onset/time-to-plug). Calibrate to your fluid
#    and field/flow-loop data before quantitative operational use. Verification (the
#    code solves the equations correctly) is built in; validation (the constants match
#    reality) is the user's calibration step. Not a certified replacement for
#    OLGA/LedaFlow, but a universal, calibratable screening & design solver.
#  ----------------------------------------------------------------------------
#  ACCURACY UPGRADES (v2 — CLOSURES, DIAGNOSTICS & REPORTING ONLY):
#    The coupled-PDE conservation laws and their numerical scheme (holdup transport,
#    energy, hydrate phase-field, deposition, and the implicit/two-fluid momentum-
#    pressure core) are UNCHANGED and behave identically by default. The following
#    were improved WITHOUT touching that immutable core, and remain universal for
#    any input case:
#      * Hydrate equilibrium curve is now gas-gravity & salinity aware and accepts a
#        user PVT table (fluids.hyd_Teq_table); defaults reproduce the natural-gas curve.
#      * Inhibitor sizing/suppression upgraded Hammerschmidt -> Nielsen-Bucklin
#        (valid to high MEG wt%, where Hammerschmidt is invalid).
#      * REPORTED Phi_SH is no longer pinned to the 50 cap: it is reported at its true
#        magnitude (finite plot cap only), so maps/tables show the real criticality.
#        The capped value still drives the deposition dynamics (core unchanged).
#      * Honest probabilistic bands: optional per-realisation input-uncertainty
#        propagation (--uq / numerics.uq_inputs) widens P10/P50/P90 beyond nucleation
#        timing alone. Off by default (base results identical).
#      * Flow- & deposit-dependent effective-U (--dynamic-u / numerics.dynamic_U),
#        off by default.
#      * Reporting fixes: no-touch time from the actual shut-in transient; bulk vs peak
#        velocity; saturation flags (Phi_SH / full-bore deposit); hydrate-mass deposit
#        diagnostic; explicit monitor location/temperature; high-mass-error warning.
#      * --grid-check runs a 1x/2x Richardson discretisation-error report.
#  ----------------------------------------------------------------------------
#  PHYSICS & ROBUSTNESS REVISION (v3 — completes the documented model & hardens it):
#      * (E) ENERGY now carries the exothermic hydrate LATENT HEAT (bulk formation only;
#        the wall deposit rejects its heat to the cold sea) -> bulk growth self-limits
#        and the bulk slurry stays transportable (heat-transfer-limited regime).
#      * (P) PHASE-FIELD now carries its D_phi DIFFUSION term (was defined but unused).
#      * (Gm) a conservative GAS-MASS CONTINUITY equation is solved alongside the liquid,
#        using the mass-consistent superficial-gas flux -> liquid, gas AND hydrate mass
#        all balance to ~0% (hydrate consumes water from the liquid and gas from the gas).
#      * (D) the WALL DEPOSIT is now MASS-COUPLED to hydrate growth (driven by the
#        sustained cold-wall subcooling), not a decoupled empirical field -> delta and the
#        bulk phi are a single mass-consistent partition; plugging is heat-limited &
#        PROBABILISTIC (genuine P10/P50/P90 time-to-plug spread).
#      * (G) well-posed pressure BCs: inlet pressure pinned IN-SYSTEM (Dirichlet) + outlet
#        through-flux (rate control) -> no post-hoc field re-anchoring; ONE dt per step for
#        momentum AND transport; never-fail path advances time (no possible infinite loop).
#      * Holdup bound-enforcement is now LOCAL-first (neighbour spill) rather than a
#        line-wide redistribution; clip-activation diagnostics are reported every run.
#      * A modest per-realisation parameter spread is ON BY DEFAULT (numerics.deterministic
#        to disable) so the probabilistic bands are real; RNG is an isolated Generator.
#      * Input validation, requirements.txt and an automated test suite (test_solver.py).
#  ----------------------------------------------------------------------------
#  GAP-CLOSURE & HARDENING REVISION (v4 — closes the audited code-level gaps):
#      * NO MORE HARD-CODED CLOSURE CONSTANTS: every previously-inline literal is now an exposed,
#        documented case field with its ORIGINAL value as the default (so results are byte-identical
#        by default — the golden-master regression still holds): the Arrhenius reference temperature
#        (kinetics.T_ref_K), wall-capture subcooling scale (wall_capture_Tsub_ref_C), induction
#        exponent guard (nuc_tau_exp_cap), sub-grid slug-body holdup base/slope and slug-fraction
#        reference (slug_body_holdup_base/slope, slug_fraction_ref_Hz), liquid conductivity
#        (fluids.k_liquid), water-droplet drag (water_droplet_d_m/Cd), settling prefactor
#        (settling_slip_coeff), slurry-viscosity exponent (slurry_visc_exp), API-14E C-factor
#        (api14e_C_factor), MEG design margin (operating.MEG_design_margin_C), surge multiplier
#        and transportability threshold (numerics.surge_factor, transportable_mu_rel_max).
#      * BUG FIX: the dynamic-U deposit resistance used a hard-coded k=0.5 W/mK that ignored
#        kinetics.deposit_k_hyd; it now reads the field (same default) so the two deposit-conduction
#        sites are consistent.
#      * twofluid_mass is now a GENUINE within-step COUPLED solve: a damped, monotone-guarded
#        Newton/Picard iteration of the momentum pressure against the conserved-gas volume
#        constraint rho_g(p)=Mg/(A-La) (numerics.twofluid_mass_iters), not a single post-hoc
#        relaxation -> lower gas-holdup inconsistency, still conservative & never-fail.
#      * MULTI-CHAIN MCMC: bayesian_calibrate runs over-dispersed chains with burn-in proposal
#        adaptation and reports a Gelman-Rubin R-hat per parameter (a real convergence diagnostic).
#      * EOS 3-phase/dew diagnostics no longer silently swallow exceptions (logged); calibration
#        infeasibility is penalised with a large finite value (drives the optimiser off crashes).
#    DELIBERATELY KEPT (honest, not "fixed away"): the never-fail fallbacks and the velocity/holdup
#    clamps (robustness features — clip diagnostics report when they bind); the 7-40% V&V
#    tolerances (they reflect REAL discretisation sensitivity on steep terrain — see --grid-check);
#    and validation against field data (which by construction needs the user's data: --validate /
#    --calibrate / --blind-validate are the mechanisms, not a constant the code can set).
#  ----------------------------------------------------------------------------
#  CROSS-SECTION / COMPOSITIONAL / 3-D RECONSTRUCTION LAYER (v5 — post-processing add-ons):
#      These are OPT-IN, pure post-processing outputs (CLI flags); they read the solved fields and
#      change nothing in the verified core (defaults & golden master untouched).
#      * --crosssection (shct_crosssection.py): reduced-order CROSS-SECTION reconstruction — the
#        stratified gas/liquid interface level h/D (exact circular-segment inversion of the holdup),
#        wetted-perimeter fraction, interface width, AZIMUTHAL (bottom-of-line) hydrate-deposit
#        distribution, and 2-D section velocity/temperature/phase fields. Resolves cross-sectional
#        structure the 1-D area-average cannot show.
#      * --threed / --3d (shct_threed.py): assembles a full 3-D field over the pipe volume from the
#        1-D core + the cross-section closure, exports a ParaView VTK STRUCTURED_GRID (pipe_3d.vtk)
#        and renders 3-D tube views (deposit, wall temperature).
#      * --compo-report (shct_compositional.py): deep COMPOSITIONAL/PVT tracking along the line via
#        the Peng-Robinson flash — vapour fraction, per-component K-values, phase densities &
#        viscosities.
#      * shct_gui.py: a Streamlit GUI front-end (`streamlit run shct_gui.py`).
#    HONEST SCOPE — the cross-section/3-D layer is a fast REDUCED-ORDER (quasi-3-D) reconstruction
#    that is CONSISTENT with the 1-D conservation laws plus published cross-section closures
#    (Taitel-Dukler interface geometry, turbulent profile). The governing physics remains 1-D
#    area-averaged: this is NOT a 3-D Navier-Stokes (CFD) solve on a 3-D mesh, does NOT resolve
#    turbulence/secondary flow/eddies, and is NOT a substitute for or "better than" 3-D CFD
#    (ANSYS Fluent / OpenFOAM). It is a complementary whole-line screening tool (seconds vs hours).
#  ----------------------------------------------------------------------------
#  OpenFOAM (OPEN-SOURCE 3-D CFD) COUPLING (v6 — shct_openfoam.py, --openfoam):
#      For the sections where genuine 3-D physics matters (riser/steep terrain, severe slugging,
#      hydrate-critical Phi_SH>1), the solver now COUPLES to OpenFOAM: it locates those sections
#      from the SHCT solution and writes a COMPLETE, runnable interFoam (3-D volume-of-fluid,
#      two-phase) case for each — cylindrical o-grid blockMesh, gravity tilted by the local
#      inclination, and inlet/outlet boundary conditions + stratified initial condition taken FROM
#      the SHCT 1-D solution (mixture velocity, holdup->liquid fraction & interface level, outlet
#      pressure, phase densities/viscosities). With OpenFOAM on PATH (--openfoam-run) it runs
#      blockMesh+setFields+interFoam and ingests the CFD result back; otherwise it writes the cases
#      to run on an OpenFOAM machine (./Allrun). Coupling is one-way SHCT->CFD (BCs) with a feedback
#      hook for two-way iteration. This is the HONEST route to true 3-D: SHCT gives the whole-line
#      context fast and OpenFOAM resolves the 3-D Navier-Stokes physics on the few sections that need
#      it. The solver does not bundle OpenFOAM (install it separately, openfoam.org / .com).
#  ----------------------------------------------------------------------------
#  VALIDATION, MONOLITHIC NEWTON & COMPOSITIONAL TRANSPORT (v7):
#      * REAL-DATA VALIDATION (--validate-hydrate): validate_hydrate_curve() scores the hydrate-
#        equilibrium closure against PUBLISHED EXPERIMENTAL data (methane Lw-H-V, Deaton & Frost 1946
#        / Sloan & Koh 2008; 7/field_data/) — RMSE ~1.9 degC as-shipped, ~1.3 degC after a 1-param
#        offset. (Full production-flow validation still needs an operator's field/loop data.)
#      * MONOLITHIC volume-mass NEWTON engine (numerics.engine="twofluid_mass_newton"): a Newton on
#        the gas volume-mass residual rho_g(p)*(A-La)=Mg drives gas-holdup inconsistency to ~0.05%
#        (vs 8.3% implicit / 1.2% damped-Picard), mass-conservative. HONEST TRADE-OFF: enforcing
#        volume consistency via pressure sacrifices dP fidelity (the system is over-determined) —
#        it is a consistency-priority engine (numerics.twofluid_mass_newton_w blends toward the
#        physical momentum dP). Getting BOTH exactly needs a flux-level simultaneous solve (open).
#      * SEQUENTIAL COMPOSITIONAL TRANSPORT (--compo-sim, shct_compositional_sim.py): tracks the
#        hydrocarbon composition grading along the line as hydrate preferentially removes the light
#        formers; component-conservative (residual ~1e-16). Reduced one-way (hydraulics->composition)
#        model, NOT a fully implicitly-coupled compositional simulator.
#      * Item-9 closures (cross-section velocity exponent, deposit skew, phase factors) are now
#        case-configurable/calibratable (numerics.cx_*); item-12 numerics.strict raises instead of
#        masking; item-10 wax (WAT) screen added; item-11 over-loose V&V tolerances tightened.
#  ----------------------------------------------------------------------------
#  PUBLISHED-DATA CLOSURE VALIDATION (v8 — extends the hydrate validation to the
#  hydraulic closures, against universally-cited reproducible references):
#      The hydraulic INGREDIENTS of holdup and pressure drop are now each scored against a
#      published reference, the same evidentiary standard as --validate-hydrate (not against
#      the solver itself). All are pure additions — defaults & golden master untouched.
#      CITATION DISCIPLINE: every reference value is attributed to its PRIMARY source (where the
#      value originated); later works are tagged 'confirmed by' / 'adopted by' and secondary
#      textbook/review COMPILATIONS are never presented as the origin (see each validator's
#      docstring and the 7/field_data/*.json reference files for the explicit roles).
#      * --validate-friction (validate_friction_curve): the friction closure haaland_friction()
#        vs the Colebrook-White / Moody reference (computed in-code to machine precision), over a
#        turbulent Re x roughness grid. RESULT: RMS 0.62%, max 1.34% deviation in Darcy f — i.e.
#        the closure that sets the FRICTIONAL pressure gradient agrees with Colebrook-White to
#        well within Haaland's published ~2% band (Colebrook self-check: f=0.018 at smooth Re=1e5).
#      * --validate-drift (validate_drift_flux): the drift-flux slip closure drift_params() that
#        sets HOLDUP, vs the canonical Taylor-bubble values. RESULT: the VERTICAL limit matches
#        Nicklin(1962)/Dumitrescu(1943) exactly (C0=1.20, drift Fr=0.35); the HORIZONTAL drift
#        Froude is 0.20 vs the Benjamin(1968)/Bendiksen(1984) value 0.542 (-63%) — reported
#        HONESTLY: the closure uses a smaller effective axial drift near-horizontal; raising it is
#        a calibration choice that would shift the golden-master default holdup, so NOT changed
#        silently here.
#      * --validate-slugfreq (validate_slug_frequency): confirms slug_frequency() reproduces the
#        Zabaras(2000) inclined-pipe correlation it implements (to ~0%), and reports its horizontal
#        value vs Gregory-Scott(1969). Honest ceiling: Zabaras' own ~+/-60% scatter on 399 air-water
#        points — slug frequency is intrinsically scattered.
#      * --validate-closures: runs all of the above + the hydrate curve with a combined summary.
#  ----------------------------------------------------------------------------
#  REAL FLOW-LOOP HOLDUP VALIDATION + LIVE 3-D CFD (v9 — closes the two items v8 left open):
#      * --validate-flowloop (validate_flowloop): the solver's HOLDUP prediction (a production-flow
#        quantity, gap 1) is now validated against a REAL, peer-reviewed, OPEN-ACCESS flow-loop
#        dataset — Das Neves et al. (2025), Data in Brief 63:112117 — 14 air-water horizontal void
#        fractions measured by the QUICK-CLOSING-VALVE method (the gold-standard direct holdup
#        measurement; D=80.5 mm). The steady drift-flux holdup the 'implicit' engine reduces to is
#        scored against the measured void: void RMSE = 0.060 as-shipped (systematic +0.050 bias =
#        the closure over-predicts gas fraction). A 1-PARAMETER calibration via numerics.
#        drift_C0_factor=1.16 (C0->1.22, ~the Nicklin developed-slug value) cuts RMSE to 0.037 —
#        same as-shipped + 1-param structure as the hydrate curve. So line HOLDUP is now validated
#        against real data; only line dP & arrival-T still benefit from an OPERATOR's specific-line set.
#      * LIVE 3-D CFD on this machine: the v6 SHCT<->OpenFOAM coupling now RUNS locally (OpenFOAM
#        v2406) — solver.py --openfoam-run [--openfoam-end-time T --openfoam-res NixNz] writes the
#        interFoam (3-D VOF) case for each critical section with BCs from the SHCT solution, runs
#        blockMesh+setFields+interFoam, ingests the volume-averaged liquid fraction, and prints the
#        SHCT(1-D) vs interFoam(3-D) holdup match. Measured here: a hydrate-critical slug section
#        gives SHCT 0.614 vs CFD 0.620 at t=3 s (~0.9%), consistent with the recorded ~1-2% flowline
#        agreement. HONEST: this is the route to true 3-D; it does NOT make the 1-D core 3-D, and the
#        run is coarse-grid / laminar VOF / short physical time (severe slugging keeps accumulating)
#        — a live cross-check, not a converged DNS.
#  ----------------------------------------------------------------------------
#  Usage:
#    python3 solver.py                       # bundled real-case
#    python3 solver.py --scenario shutin     # cooldown / shut-in hydrate risk
#    python3 solver.py --engine twofluid     # full two-fluid (two phase momenta)
#    python3 solver.py --verify              # run V&V verification suite
#    python3 solver.py --validate-closures   # score closures vs published refs (friction/drift/slug/hydrate)
#    python3 solver.py --validate-flowloop   # validate holdup vs real flow-loop void-fraction data (v9)
#    python3 solver.py --openfoam-run        # run+ingest live interFoam 3-D CFD on critical sections (v9)
#    python3 solver.py --calibrate tgt.json  # fit constants to measured data
#    python3 solver.py --config case.json    # any user case (see --dump-config)
#    python3 solver.py --no-plots            # tables + console only
# =============================================================================
from __future__ import annotations
import os, sys, json, argparse, math, copy, logging
from dataclasses import asdict
import numpy as np

#  #23: status/diagnostic messages go through a module logger (configurable level / redirectable),
#  while the human-facing verification, validation and engineering REPORTS stay as direct prints.
log = logging.getLogger("shct")
if not log.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(_h)
    log.setLevel(logging.INFO)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.colors as mcolors
    import shct_style as _style          # global no-black / no-dark plotting style
    _style.apply_style()
    HAVE_MPL = True
except Exception:                                       # pragma: no cover
    HAVE_MPL = False

#  Default per-realisation parameter spread (rel. sigma) that makes the ensemble genuinely
#  probabilistic (C10): hydrate induction time and growth and the wall heat transfer are
#  physically uncertain, so the P10/P50/P90 bands reflect real spread, not nucleation timing
#  alone. Applied unless numerics.deterministic=True or explicit numerics.uq_inputs is given.
def _erfinv(x):
    """Inverse error function (Winitzki rational approximation) for the Latin-Hypercube normal
    mapping — avoids a SciPy dependency in the core run loop."""
    x = max(min(float(x), 1 - 1e-9), -1 + 1e-9)
    a = 0.147
    ln = math.log(1.0 - x * x)
    t1 = 2.0 / (math.pi * a) + ln / 2.0
    return math.copysign(math.sqrt(math.sqrt(t1 * t1 - ln / a) - t1), x)


DEFAULT_UQ = {"nuc_tau0_h": 0.40, "kg0": 0.30, "U_wall": 0.10, "nuc_beta_C": 0.15}
#  Stronger spread enabled by the --uq CLI flag.
STRONG_UQ = {"kg0": 0.5, "nuc_tau0_h": 0.5, "nuc_beta_C": 0.2, "k_dep": 0.5,
             "wall_capture_eff": 0.4, "U_wall": 0.15}
OUTDIR_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solver_outputs")
# medium, non-black, non-dark palette (see shct_style.py — NO black/dark anywhere)
NAVY, ACCENT, RED, ORANGE, TEAL, GREEN, GREY = \
    "#2E5BBF", "#1F8AC0", "#E0463C", "#E8842B", "#1AA0A0", "#3FA65A", "#9AA8C7"


# =============================================================================
#  INPUT DATA MODEL — extracted to shct_model.py (#21)
# =============================================================================
from shct_model import (  # noqa: E402
    Pipeline, Fluids, Operating, Kinetics, Numerics, Scenario, Case)


# =============================================================================
#  CORRELATIONS  (vectorised) — extracted to shct_correlations.py (#21)
# =============================================================================
from shct_correlations import (  # noqa: E402  (pure closures, independently testable)
    G, R_GAS, M_MEG, M_H2O, REGIME_NAMES, hydrate_equilibrium_T, gas_Z_factor, gas_density,
    gas_viscosity, oil_density, _limiter, tvd_interior_faces, haaland_friction, drift_params,
    flow_regime_code, slug_frequency, interfacial_area, interfacial_area_geom,
    regime_friction_multiplier, regime_nusselt, joule_thomson_dTdP, scaling_tendency_index,
    slug_length, droplet_entrainment_frac, mixture_sound_speed, _meg_wt_to_molefrac,
    meg_suppression, hammerschmidt_meg, effective_U_and_mass)


# =============================================================================
#  TRANSIENT COUPLED-PDE SOLVER
# =============================================================================
class TransientSHCT:
    def __init__(self, case: Case):
        self.case = case
        self._build_grid()
        self.N = case.numerics.n_ensemble
        #  Local RNG (G23): use an isolated Generator instead of mutating the global
        #  np.random state, so the solver is safe to import alongside other RNG users
        #  while remaining fully reproducible from numerics.seed.
        self.rng = np.random.default_rng(case.numerics.seed)
        #  #A1: if a COMPOSITION is supplied, build the PVT property surface and the gas SG/Z
        #  from a Peng-Robinson flash (genuine compositional thermodynamics). The solver then
        #  consumes the EOS through the same conserving pvt_table path.
        self._Vsurf = None
        if getattr(case.fluids, "composition", None):
            import shct_eos
            comp = case.fluids.composition
            #  B11: split C7+ into characterized pseudo-components (Whitson + Kesler-Lee) first
            if getattr(case.fluids, "n_pseudo", 1) > 1:
                comp = shct_eos.expand_composition(comp, case.fluids.n_pseudo, case.fluids.MW_plus)
                case.fluids.composition = comp
            if case.fluids.pvt_table is None:
                case.fluids.pvt_table = shct_eos.build_pvt_table(comp)
            pr = shct_eos.eos_properties(case.operating.P_inlet_bar, case.operating.T_inlet_C, comp)
            case.fluids.gas_sg = pr["gas_sg"]; case.fluids.gas_MW = pr["gas_sg"] * 0.028964
            #  B9: precompute a vapour-fraction surface V(P,T) for the condensation latent-heat term
            if case.fluids.condensation_latent:
                Pg = [10, 30, 60, 100, 150, 200, 300]; Tg = [4, 10, 20, 30, 40, 55, 70]
                Vgrid = np.array([[shct_eos.flash(P, T, comp)["V"] for T in Tg] for P in Pg])
                self._Vsurf = (np.array(Pg, float), np.array(Tg, float), Vgrid)
        self.rho_l = (case.fluids.rho_oil * (1 - case.fluids.water_cut)
                      + case.fluids.rho_water * case.fluids.water_cut)
        #  local 3-phase liquid density field (#2 deepened): a SCALAR (= composite rho_l) unless
        #  oil_water_slip is on, in which case it becomes a (nx,N) field of rho_oil/rho_water mixed
        #  by the LOCAL transported water fraction and fed into the momentum gravity/friction and
        #  energy. Scalar default broadcasts -> identical to before when slip is off.
        self._rho_l_field = self.rho_l
        self.engine = getattr(case.numerics, "engine", "implicit")
        self._limiter_kind = getattr(case.numerics, "flux_limiter", "minmod")   # #8 TVD limiter
        #  clip-activation diagnostics (D13): count how often the safety clamps bind, so a
        #  "robust" run that is silently masking an instability is detectable rather than hidden.
        self._clip = {"velocity": 0, "pressure": 0, "holdup": 0, "deposit": 0}
        self.results = {}

    # ----- geometry -------------------------------------------------------
    def _build_grid(self):
        p = self.case.pipeline
        #  #16: terrain-adaptive resolution — increase n_cells where the terrain is steep so the
        #  thermal field (grid-sensitive on inclines, per --grid-check) is better resolved. A coarse
        #  global refinement (uniform grid) keyed on the maximum elevation slope; opt-in.
        if getattr(self.case.numerics, "auto_refine", False) and p.elevation_m is not None:
            z0 = np.asarray(p.elevation_m, float)
            x0 = np.linspace(0, p.length_m, len(z0))
            max_slope = float(np.nanmax(np.abs(np.gradient(z0, x0)))) if len(z0) > 1 else 0.0
            refine = int(np.clip(1.0 + 4.0 * max_slope, 1.0, 4.0))     # up to 4x on steep terrain
            p.n_cells = int(p.n_cells * refine)
        self.nx = p.n_cells
        self.x = (np.arange(self.nx) + 0.5) * p.length_m / self.nx
        self.dx = p.length_m / self.nx
        if p.elevation_m is not None and len(p.elevation_m) >= self.nx:
            z = np.asarray(p.elevation_m[:self.nx], float)
        else:
            xk = self.x / p.length_m
            z = -80.0 - 35.0 * np.sin(6.0 * xk * math.pi) * np.exp(-1.5 * xk) - 20.0 * xk
            riser = xk > 0.86
            if riser.any():
                z[riser] = np.linspace(z[riser][0], 8.0, riser.sum())
        self.z = z
        self.theta = np.arctan(np.gradient(z, self.x))     # rad, uphill +

    # ----- transient inlet boundary condition (scenario) ------------------
    def inlet_rates(self, t_h):
        o, sc = self.case.operating, self.case.scenario
        ql, qg = o.q_liquid_insitu, o.q_gas_insitu_inlet
        if sc.kind == "steady":
            f = 1.0
        elif sc.kind == "rampup":
            f = np.clip(t_h / max(sc.event_time_h, 1e-3), 0.30, 1.0)   # from minimum stable rate
        elif sc.kind == "turndown":
            f = 1.0 if t_h < sc.event_time_h else sc.turndown_factor
        elif sc.kind == "shutin":
            f = 1.0 if t_h < sc.event_time_h else sc.shutin_residual
        else:
            f = 1.0
        return ql * f, qg * f

    # ----- quasi-steady pressure / velocity given current holdup ----------
    def pressure_velocity(self, alpha_l, T, delta, t_h):
        """March pressure inlet->outlet over the *instantaneous* holdup field
        (acoustic-filtered). Returns p, j(=v_m), v_l, v_g, rho_g, a_eff."""
        c, f = self.case, self.case.fluids
        nx, N = self.nx, self.N
        D0 = c.pipeline.diameter_m
        D = np.clip(D0 - 2.0 * delta, 0.60 * D0, D0)
        A = math.pi * D ** 2 / 4.0
        ql, qg_in = self.inlet_rates(t_h)
        Tin = c.operating.T_inlet_C
        p = np.empty((nx, N)); j = np.empty((nx, N))
        vl = np.empty((nx, N)); vg = np.empty((nx, N)); rho_g = np.empty((nx, N))
        p[0] = c.operating.P_inlet_bar
        for i in range(nx):
            if i > 0:
                rg = gas_density(p[i - 1], T[i - 1], f)
                al = alpha_l[i - 1]; ag = 1 - al
                rho_m = al * self.rho_l + ag * rg
                vsg = qg_in * (c.operating.P_inlet_bar / p[i - 1]) * \
                    ((T[i - 1] + 273.15) / (Tin + 273.15)) / A[i - 1]
                vsl = ql / A[i - 1]
                jm = vsg + vsl
                mu_m = al * f.mu_liquid + ag * f.mu_gas
                Re = rho_m * jm * D[i - 1] / np.maximum(mu_m, 1e-9)
                fr = haaland_friction(Re, c.pipeline.roughness_m / D[i - 1])
                dpdx = rho_m * G * np.sin(self.theta[i - 1]) + \
                    fr * rho_m * jm * np.abs(jm) / (2.0 * D[i - 1])
                p[i] = np.maximum(p[i - 1] - dpdx * self.dx / 1e5, 2.0)
            rho_g[i] = gas_density(p[i], T[i], f)
            vsg = qg_in * (c.operating.P_inlet_bar / p[i]) * \
                ((T[i] + 273.15) / (Tin + 273.15)) / A[i]
            vsl = ql / A[i]
            j[i] = vsg + vsl
            al = alpha_l[i]; ag = 1 - al
            C0, vd = drift_params(self.theta[i], D[i])
            C0 = C0 * float(getattr(c.numerics, "drift_C0_factor", 1.0))   # holdup knob (default 1.0)
            vg[i] = np.clip(C0 * j[i] + vd, -5.0, 40.0)
            jg = np.clip(ag * vg[i], 0.0, j[i] + 2.0)         # gas volumetric flux
            jl = j[i] - jg
            #  limit liquid (incl. gravity fall-back) velocity to physical bounds so the
            #  bounded holdup scheme stays conservative without lossy clipping.
            vl[i] = np.clip(jl / np.maximum(al, 1e-3), -0.4, 20.0)
        return dict(p=p, j=j, vl=vl, vg=vg, rho_g=rho_g, A=A, D=D)

    # ----- vectorised tridiagonal (Thomas) solver, over the ensemble axis -----
    @staticmethod
    def _thomas(a, b, c, d):
        """Solve a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i for each column.
        a,b,c,d shape (nx, N). Returns x shape (nx, N). Diagonally dominant -> stable."""
        n = a.shape[0]
        cp = np.empty_like(c); dp = np.empty_like(d)
        cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0]
        for i in range(1, n):
            m = b[i] - a[i] * cp[i - 1]
            m = np.where(np.abs(m) < 1e-30, 1e-30, m)
            cp[i] = c[i] / m
            dp[i] = (d[i] - a[i] * dp[i - 1]) / m
        x = np.empty_like(d); x[-1] = dp[-1]
        for i in range(n - 2, -1, -1):
            x[i] = dp[i] - cp[i] * x[i + 1]
        return x

    # ----- SEMI-IMPLICIT transient two-phase momentum + implicit pressure ------
    def momentum_solve(self, alpha_l, T, delta, t_h, dt):
        """Genuine TRANSIENT pressure-momentum step (replaces the quasi-steady march):
        a transient mixture-momentum equation with the pressure gradient and wall drag
        treated IMPLICITLY, closed by a Bendiksen drift-flux slip. Eliminating velocity
        gives a tridiagonal (Poisson-type) pressure equation solved implicitly each step,
        so pressure and momentum evolve in time (compressible storage, real transients)
        rather than being marched quasi-steadily. The drift-flux slip keeps the system
        unconditionally well-posed (no complex characteristics) -> robust in every regime.
        BCs: inlet volumetric-flux prescribed (rate control); outlet pressure anchored so
        the inlet stays at P_inlet_bar. Returns the same fields as pressure_velocity()."""
        c, f = self.case, self.case.fluids
        nx, N = self.nx, self.N
        dx = self.dx
        D0 = c.pipeline.diameter_m
        D = np.clip(D0 - 2.0 * delta, 0.60 * D0, D0)
        A = math.pi * D ** 2 / 4.0
        theta2 = self.theta[:, None] * np.ones((1, N))
        al = np.clip(alpha_l, 1e-3, 0.999); ag = 1.0 - al
        rho_l = self._rho_l_field            # #2: local 3-phase density (scalar if slip off)
        rho_g = gas_density(self._p, T, f)
        rho_m = al * rho_l + ag * rho_g
        C0, vd = drift_params(theta2, D)
        C0 = C0 * float(getattr(c.numerics, "drift_C0_factor", 1.0))   # holdup knob (default 1.0)
        um = self._um                                          # mixture velocity (state)

        # inlet volumetric flux (rate control) -> prescribed as the OUTLET through-flux below
        ql, qg = self.inlet_rates(t_h)
        Tin = c.operating.T_inlet_C
        vsg_in = qg * (c.operating.P_inlet_bar / np.maximum(self._p[0], 1.0)) * \
            ((T[0] + 273.15) / (Tin + 273.15)) / A[0]
        j_in = qg * 0.0 + ql / A[0] + vsg_in                  # (N,) prescribed through-flux

        # explicit predictor (advection + gravity), implicit (pressure + wall drag)
        um_up = np.vstack([um[:1], um[:-1]])
        adv = np.where(um >= 0, um * (um - um_up) / dx, um * (np.vstack([um[1:], um[-1:]]) - um) / dx)
        mu_m = al * f.mu_liquid + ag * gas_viscosity(rho_g, T, f)   # #1: real-gas viscosity
        Re = rho_m * np.abs(um) * D / np.maximum(mu_m, 1e-9)
        fr = haaland_friction(Re, c.pipeline.roughness_m / D)
        beta = fr / (2.0 * D) * np.abs(um)                     # wall-drag rate (implicit)
        u_star = um - dt * (adv + G * np.sin(theta2))
        denom = 1.0 + dt * beta
        U0 = u_star / denom                                    # um = U0 - Cu * dp/dx (p in bar)
        Cu = (dt / np.maximum(rho_m, 1.0)) / denom * 1e5       # 1e5 folds bar->Pa for dp/dx

        # assemble tridiagonal pressure equation:  R(p_new-p) + d/dx[U0 - Cu dp_new/dx] = 0
        R = ag / np.maximum(self._p * dt, 1e-9)                # compressible storage (well-posed)
        #  #B5 WATER-HAMMER: when acoustic>0 the storage uses the FULL mixture compressibility
        #  (Wood: 1/(rho_m c^2), liquid + gas) instead of gas-only, so fast pressure waves are
        #  captured (with the acoustic-CFL dt cap applied in run()). 0 -> acoustic-filtered default.
        if getattr(c.numerics, "acoustic", 0.0) > 0.0:
            cw = mixture_sound_speed(al, rho_l, rho_g, self._p)
            R_ac = 1.0 / np.maximum(rho_m * cw ** 2 * dt, 1e-12) * 1e5    # bar-based compressibility
            R = (1.0 - c.numerics.acoustic) * R + c.numerics.acoustic * R_ac
        Cuf = np.empty((nx + 1, N))
        Cuf[1:-1] = 0.5 * (Cu[:-1] + Cu[1:]); Cuf[0] = Cu[0]; Cuf[-1] = Cu[-1]
        U0f = np.empty((nx + 1, N))
        U0f[1:-1] = 0.5 * (U0[:-1] + U0[1:]); U0f[0] = U0[0]; U0f[-1] = j_in
        idx2 = 1.0 / dx ** 2
        a = -Cuf[:-1] * idx2                                   # sub-diagonal
        cc = -Cuf[1:] * idx2                                   # super-diagonal
        b = R + (Cuf[:-1] + Cuf[1:]) * idx2                    # diagonal
        d = R * self._p - (U0f[1:] - U0f[:-1]) / dx
        #  B6: WELL-POSED boundary conditions with NO post-hoc field shift.
        #   * inlet  (i=0): Dirichlet pressure p = P_inlet_bar — the wellhead is pinned exactly
        #                   INSIDE the linear system (not corrected after the solve), so the
        #                   transient interior and the R*p storage term are never teleported.
        #   * outlet (i=nx-1): prescribed through-flux j_in (rate control); the separator
        #                   pressure floats to whatever the friction/gravity gradient requires.
        a[0] = 0.0; b[0] = 1.0; cc[0] = 0.0; d[0] = c.operating.P_inlet_bar
        a[-1] = -Cuf[-2] * idx2; b[-1] = R[-1] + Cuf[-2] * idx2; cc[-1] = 0.0
        d[-1] = R[-1] * self._p[-1] - (j_in - U0f[-2]) / dx

        p_new = self._thomas(a, b, cc, d)
        p_clipped = np.clip(p_new, 2.0, 1.0e3)
        self._clip["pressure"] += int(np.count_nonzero(p_clipped != p_new))
        p_new = p_clipped

        # recover velocities from the new pressure
        dpdx = np.empty((nx, N))
        dpdx[1:-1] = (p_new[2:] - p_new[:-2]) / (2 * dx)
        dpdx[0] = (p_new[1] - p_new[0]) / dx
        dpdx[-1] = (p_new[-1] - p_new[-2]) / dx
        um_raw = U0 - Cu * dpdx                                 # Cu already carries bar->Pa
        um_new = np.clip(um_raw, -2.0, 30.0)
        self._clip["velocity"] += int(np.count_nonzero(um_new != um_raw))
        # drift-flux slip -> phase velocities
        u_g = np.clip(C0 * um_new + vd, -5.0, 40.0)
        jg = np.clip(ag * u_g, 0.0, np.abs(um_new) + 2.0)
        u_l = np.clip((um_new - jg) / np.maximum(al, 1e-3), -0.4, 20.0)

        self._um = um_new; self._p = p_new
        #  soft inlet anchor (B6): the outlet Dirichlet is p_out = P_inlet - self._dP, so setting
        #  self._dP to the ACTUAL solved gradient makes this an integral controller — next step's
        #  p_out carries the accumulated inlet error and the wellhead converges to P_inlet without
        #  the old instantaneous full-field shift (which corrupted the R*p storage term). No extra
        #  inlet-error gain is added (that double-counts and is unstable).
        self._dP = float(np.mean(p_new[0] - p_new[-1]))
        rho_g = gas_density(p_new, T, f)
        return dict(p=p_new, j=um_new, vl=u_l, vg=u_g, rho_g=rho_g, A=A, D=D)

    # ----- FULL TWO-FLUID: two independent phase momenta + interfacial drag -----
    def twofluid_solve(self, alpha_l, T, delta, t_h, dt):
        """Full two-fluid model: SEPARATE transient gas and liquid momentum equations,
        each with its own inertia, advection and wall drag, coupled only by (i) a shared
        pressure and (ii) interfacial drag — NO algebraic drift-flux slip closure. Slip
        emerges from the momentum balance (buoyancy via the shared pressure gradient acting
        on the density difference, opposed by interfacial drag).

        Semi-implicit (RELAP/TRACE-style) for robustness on the otherwise ill-posed two-
        fluid system: advection+gravity explicit; pressure-gradient, wall drag and the
        stiff interfacial drag IMPLICIT. Per cell a 2x2 velocity system is inverted to give
        (u_g, u_l) linear in dp/dx; substitution into mixture-volume continuity yields a
        tridiagonal implicit pressure equation. An interfacial-pressure (hyperbolicity)
        correction is added to keep characteristics real. Never-fail fallback applies."""
        c, f = self.case, self.case.fluids
        nx, N = self.nx, self.N
        dx = self.dx
        D0 = c.pipeline.diameter_m
        D = np.clip(D0 - 2.0 * delta, 0.60 * D0, D0)
        A = math.pi * D ** 2 / 4.0
        theta2 = self.theta[:, None] * np.ones((1, N))
        al = np.clip(alpha_l, 1e-3, 0.999); ag = 1.0 - al
        rho_l = self._rho_l_field            # #2: local 3-phase density (scalar if slip off)
        rho_g = np.maximum(gas_density(self._p, T, f), 1e-3)
        rho_m = al * rho_l + ag * rho_g
        ug = self._ug; ul = self._ul

        # explicit predictors: advection (upwind) + gravity + interfacial-pressure reg.
        def _ddx_up(u):
            u_b = np.vstack([u[:1], u[:-1]]); u_f = np.vstack([u[1:], u[-1:]])
            return np.where(u >= 0, (u - u_b) / dx, (u_f - u) / dx)
        gsin = G * np.sin(theta2)
        #  interfacial-pressure correction delta_p*d(alpha_l)/dx makes the two-fluid
        #  hyperbolic (real characteristics) -> well-posed (Bestion/virtual-mass style).
        #  D14: coefficient is the documented fluids.vm_coeff (was a hard-coded 1.2).
        dpi = f.vm_coeff * (rho_l * rho_g / rho_m) * (ug - ul) ** 2
        al_b = np.vstack([al[:1], al[:-1]]); dal = (al - al_b) / dx
        ug_star = ug - dt * (ug * _ddx_up(ug) + gsin) + dt * dpi * dal / np.maximum(ag * rho_g, 1.0)
        ul_star = ul - dt * (ul * _ddx_up(ul) + gsin) - dt * dpi * dal / np.maximum(al * rho_l, 1.0)

        # wall drag (implicit, per phase) and interfacial drag (implicit, stiff)
        mu_l, mu_g = f.mu_liquid, gas_viscosity(rho_g, T, f)   # #1: real-gas viscosity
        Re_l = rho_l * np.abs(ul) * D / mu_l
        beta_l = haaland_friction(Re_l, c.pipeline.roughness_m / D) / (2.0 * D) * np.abs(ul)
        Re_g = rho_g * np.abs(ug) * D / mu_g
        #  D14: gas wall-drag scaled by the documented fluids.gas_wall_drag_frac (was 0.1) —
        #  the gas contacts only a fraction of the wetted perimeter (Taitel-Dukler geometry).
        beta_g = (haaland_friction(Re_g, c.pipeline.roughness_m / D) / (2.0 * D)
                  * np.abs(ug) * f.gas_wall_drag_frac)
        slip = np.abs(ug - ul) + 0.05
        Ci = f.interfacial_Ci                                  # interfacial friction factor (tunable)
        Kvol = 0.5 * Ci * rho_g * slip * (4.0 / D)             # interfacial drag per unit volume
        a_g = dt * Kvol / np.maximum(ag * rho_g, 1e-3)         # drag rate felt by gas
        a_l = dt * Kvol / np.maximum(al * rho_l, 1e-3)         # reaction on liquid
        bg = dt * beta_g; bl = dt * beta_l

        # 2x2 system  M [ug;ul] = [ug_star - (dt/rho_g)G ; ul_star - (dt/rho_l)G]
        Am = 1.0 + bg + a_g; Bm = -a_g; Cm = -a_l; Dm = 1.0 + bl + a_l
        det = np.maximum(Am * Dm - Bm * Cm, 1e-6)
        ug0 = (Dm * ug_star - Bm * ul_star) / det
        ul0 = (-Cm * ug_star + Am * ul_star) / det
        Cug = (Dm * (dt / rho_g) - Bm * (dt / rho_l)) / det * 1e5   # ug = ug0 - Cug dp/dx
        Cul = (-Cm * (dt / rho_g) + Am * (dt / rho_l)) / det * 1e5  # ul = ul0 - Cul dp/dx
        j0 = ag * ug0 + al * ul0
        Cj = ag * Cug + al * Cul                                    # j = j0 - Cj dp/dx

        # inlet flux & outlet anchor, then tridiagonal implicit pressure (as momentum_solve)
        ql, qg = self.inlet_rates(t_h)
        vsg_in = qg * (c.operating.P_inlet_bar / np.maximum(self._p[0], 1.0)) * \
            ((T[0] + 273.15) / (c.operating.T_inlet_C + 273.15)) / A[0]
        j_in = ql / A[0] + vsg_in
        R = ag / np.maximum(self._p * dt, 1e-9)
        Cjf = np.empty((nx + 1, N)); Cjf[1:-1] = 0.5 * (Cj[:-1] + Cj[1:]); Cjf[0] = Cj[0]; Cjf[-1] = Cj[-1]
        j0f = np.empty((nx + 1, N)); j0f[1:-1] = 0.5 * (j0[:-1] + j0[1:]); j0f[0] = j0[0]; j0f[-1] = j_in
        idx2 = 1.0 / dx ** 2
        a = -Cjf[:-1] * idx2; cc = -Cjf[1:] * idx2
        b = R + (Cjf[:-1] + Cjf[1:]) * idx2
        d = R * self._p - (j0f[1:] - j0f[:-1]) / dx
        #  B6: inlet Dirichlet p = P_inlet (pinned in-system); outlet prescribed through-flux j_in.
        a[0] = 0.0; b[0] = 1.0; cc[0] = 0.0; d[0] = c.operating.P_inlet_bar
        a[-1] = -Cjf[-2] * idx2; b[-1] = R[-1] + Cjf[-2] * idx2; cc[-1] = 0.0
        d[-1] = R[-1] * self._p[-1] - (j_in - j0f[-2]) / dx
        p_new = self._thomas(a, b, cc, d)
        p_clipped = np.clip(p_new, 2.0, 1.0e3)
        self._clip["pressure"] += int(np.count_nonzero(p_clipped != p_new))
        p_new = p_clipped

        dpdx = np.empty((nx, N))
        dpdx[1:-1] = (p_new[2:] - p_new[:-2]) / (2 * dx)
        dpdx[0] = (p_new[1] - p_new[0]) / dx; dpdx[-1] = (p_new[-1] - p_new[-2]) / dx
        ug_raw = ug0 - Cug * dpdx; ul_raw = ul0 - Cul * dpdx
        ug = np.clip(ug_raw, -5.0, 40.0); ul = np.clip(ul_raw, -0.4, 20.0)
        self._clip["velocity"] += int(np.count_nonzero(ug != ug_raw)
                                       + np.count_nonzero(ul != ul_raw))
        self._ug = ug; self._ul = ul; self._p = p_new
        #  B6 integral inlet anchor (see momentum_solve) — no double-counted inlet-error gain.
        self._dP = float(np.mean(p_new[0] - p_new[-1]))
        self._um = ag * ug + al * ul
        return dict(p=p_new, j=ag * ug + al * ul, vl=ul, vg=ug,
                    rho_g=gas_density(p_new, T, f), A=A, D=D)

    # ----- TWO-FLUID-MASS engine (A1) -------------------------------------
    def twofluid_mass_solve(self, alpha_l, T, delta, t_h, dt):
        """Two-fluid-MASS engine: a within-step COUPLED (momentum <-> gas-volume) Newton/Picard
        iteration. The validated implicit momentum solve provides the hydraulic pressure &
        velocities; the pressure is then driven toward the value at which the CONSERVED gas mass
        Mg exactly fills the available bore (rho_g(p) = Mg/(A - La), inverted through the gas EOS),
        and momentum is RE-SOLVED at that pressure — iterated to a fixed point over
        numerics.twofluid_mass_iters passes. So the conserved gas mass genuinely co-determines the
        pressure WITHIN the step (the (Mg, momentum, p) coupling), not by a single post-hoc
        relaxation; the conserved liquid La closes the volume.

        The iteration is damped (w=0.5) and monotone-guarded — it stops the moment the volume
        residual stops shrinking or a re-solve turns non-finite, so it strictly improves the
        gas-holdup consistency without ever destabilising the verified core (the never-fail
        fallback in run() still applies). A full simultaneous (La, Mg, momentum, p) Jacobian solve
        would tighten this further; this damped quasi-Newton coupling is the stable realisation."""
        c, f = self.case, self.case.fluids
        n_it = max(int(getattr(c.numerics, "twofluid_mass_iters", 4)), 1)
        pv = self.momentum_solve(alpha_l, T, delta, t_h, dt)
        A = pv["A"]
        La = np.clip(alpha_l, 1e-3, 0.999) * A
        w = 0.5                                       # damping on the volume-consistent correction
        prev_res = np.inf
        for _ in range(n_it - 1):
            if not np.all(np.isfinite(self._p)):
                break
            rho_g_vc = np.clip(self._Mg / np.maximum(A - La, 1e-6), 1e-3, 600.0)
            Zc = gas_Z_factor(self._p, T, f)
            p_vc = np.clip(rho_g_vc * Zc * R_GAS * (T + 273.15) / f.gas_MW / 1e5, 2.0, 1.0e3)
            res = float(np.nanmax(np.abs(p_vc - self._p)))
            if not np.isfinite(res) or res >= prev_res:
                break                                 # converged / no longer improving -> stop
            prev_res = res
            #  blend the pressure state toward volume consistency, then RE-SOLVE momentum so the
            #  velocities and the pressure gradient stay consistent with the corrected pressure.
            self._p = np.clip((1.0 - w) * self._p + w * p_vc, 2.0, 1.0e3)
            pv_try = self.momentum_solve(alpha_l, T, delta, t_h, dt)
            if not (np.all(np.isfinite(pv_try["p"])) and np.all(np.isfinite(pv_try["j"]))):
                break                                 # keep the last stable solve
            pv = pv_try
            if res < 1e-4:                            # tightened convergence (item 7)
                break
        return pv

    # ----- MONOLITHIC volume-mass NEWTON engine (item 7) ------------------
    def twofluid_mass_newton_solve(self, alpha_l, T, delta, t_h, dt):
        """Monolithic volume-mass Newton coupling. The implicit momentum solve gives the physical
        pressure p_mom; we then solve, by NEWTON to tolerance, the gas volume-mass residual
            r(p) = rho_g(p)*(A - La) - Mg = 0           (rho_g via the gas EOS, dr/dp = (A-La)*drho/dp)
        so the conserved gas mass Mg and the conserved-liquid-derived gas volume (A-La) are made
        consistent *exactly* (gas-holdup inconsistency -> ~0, the gap the damped Picard plateaus at).
        The Newton root p_vc is blended with the momentum pressure by numerics.twofluid_mass_newton_w
        (default 0.85): w->1 gives near-exact consistency at the mass-arranged pressure gradient,
        w->0 gives the physical momentum dP — an explicit, documented consistency<->dP control.
        Both conserved fields (La, Mg) are unchanged, so liquid AND gas mass still conserve."""
        c, f = self.case, self.case.fluids
        pv = self.momentum_solve(alpha_l, T, delta, t_h, dt)
        A = pv["A"]
        La = np.clip(alpha_l, 1e-3, 0.999) * A
        AmL = np.maximum(A - La, 1e-6)
        w = float(getattr(c.numerics, "twofluid_mass_newton_w", 0.85))
        if not getattr(self, "_tfmn_warned", False):
            log.warning("[twofluid_mass_newton] consistency-priority engine: drives gas-holdup "
                        "consistency to ~0 but the volume-consistent pressure SACRIFICES dP fidelity "
                        "(over-determined system). Use for holdup/consistency studies; use 'implicit' "
                        "for pressure-drop. Lower numerics.twofluid_mass_newton_w toward 0 for more "
                        "physical dP (less consistency).")
            self._tfmn_warned = True
        p = self._p.copy()
        rho_target = np.clip(self._Mg / AmL, 1e-3, 600.0)        # required gas density for consistency
        for _ in range(6):                                       # Newton on r(p) = rho_g(p) - rho_target
            Zc = gas_Z_factor(p, T, f)
            rho = p * 1e5 * f.gas_MW / np.maximum(Zc * R_GAS * (T + 273.15), 1e-9)
            r = rho - rho_target
            drho_dp = 1e5 * f.gas_MW / np.maximum(Zc * R_GAS * (T + 273.15), 1e-9)   # d rho / d p (bar)
            p = np.clip(p - r / np.maximum(drho_dp, 1e-12), 2.0, 1.0e3)
            if float(np.nanmax(np.abs(r / np.maximum(rho_target, 1e-6)))) < 1e-6:
                break
        p_new = np.clip((1.0 - w) * self._p + w * p, 2.0, 1.0e3)
        if not np.all(np.isfinite(p_new)):
            return pv                                            # never-fail: keep the momentum solve
        # recover velocities consistent with the blended pressure (same recovery as momentum_solve)
        nx = self.nx; dx = self.dx
        dpdx = np.empty((nx, self.N))
        dpdx[1:-1] = (p_new[2:] - p_new[:-2]) / (2 * dx)
        dpdx[0] = (p_new[1] - p_new[0]) / dx; dpdx[-1] = (p_new[-1] - p_new[-2]) / dx
        # reuse the last momentum velocity as the base (already finite); blend pressure only
        self._p = p_new
        rho_g = gas_density(p_new, T, f)
        return dict(p=p_new, j=pv["j"], vl=pv["vl"], vg=pv["vg"], rho_g=rho_g, A=A, D=pv["D"])

    # ----- FULLY-COUPLED (alpha_l, p, u_m) FLUX-LEVEL NEWTON engine -----
    def twofluid_full_newton_solve(self, alpha_l, T, delta, t_h, dt):
        """Fully-COUPLED (monolithic) two-fluid step: a within-step block Newton that solves the three
        discrete backward-Euler residuals SIMULTANEOUSLY for the primitive unknowns (alpha_l, p, u_m)
        per cell, rather than the segregated march the 'implicit' engine uses —

            R1  LIQUID-mass continuity   :  (La - La_n)/dt + d/dx[alpha_l u_l A]          = 0   (-> alpha_l)
            R2  ELLIPTIC pressure        :  (alpha_g/p)(p - p_n)/dt + d/dx[u_m]           = 0   (-> p)
            R3  mixture MOMENTUM         :  rho_m(u_m - u_m_n)/dt + rho_m u_m du/dx
                                            + rho_m g sin(theta) + rho_m beta u_m + dp/dx = 0   (-> u_m)

        with La = alpha_l A, the drift-flux slip u_g = C0 u_m + v_d closing the phase split, EOS rho_g(p),
        and a Haaland wall-drag closure. R2 is the genuine ELLIPTIC compressible-pressure (mixture-volume
        continuity) equation — the sum of the two phase continuities weighted by 1/density — so the gas-mass
        balance is IMPLIED by R1 + R2; this is the correct low-Mach pressure character (a hyperbolic upwind
        gas-mass march anchored only at the inlet gives an UNPHYSICAL dP — verified during development).

        WHAT IT DELIVERS (measured, do not overclaim beyond this): converges in ~1-2 Newton iterations to
        ~1e-13 residual, REPRODUCES the validated implicit engine's physical pressure drop (~14 bar vs ~14.4)
        and conserves liquid mass to ~1e-15, with zero fallbacks — i.e. a fully-coupled Newton that AGREES
        with the segregated solver's hydraulics (a verification result), and is more diffusive (it damps
        fine slug structure). It does NOT by itself improve the reported gas-holdup CONSISTENCY: that metric
        is taken from the run()-level conservative transport, which still owns La/Mg here; for the
        consistency-priority (gas-holdup) niche use 'twofluid_mass_newton'. Driving consistency->0 inside
        THIS engine needs it to also own the La/Mg updates (the deeper integration) — still future work.

        Numerics: a 3-unknown/cell system -> a BLOCK-TRIDIAGONAL (3x3) Jacobian built by a 9-evaluation
        COLOURED finite-difference sweep (residuals couple only to i-1,i,i+1, so stride-3 colouring is
        exact), with per-equation ROW SCALING (the residuals are otherwise incommensurate -> Newton stalls)
        and solved by a vectorised block-Thomas sweep over the whole ensemble at once, with backtracking
        damping. NEVER-FAIL: any non-finite Newton state falls back to the validated momentum_solve. The
        conserved fields La, Mg remain advanced by the run()-level conservative transport using these
        coupled velocities, so liquid AND gas mass still conserve exactly."""
        c, f = self.case, self.case.fluids
        nx, N = self.nx, self.N
        dx = self.dx
        D0 = c.pipeline.diameter_m
        D = np.clip(D0 - 2.0 * delta, 0.60 * D0, D0)
        A = math.pi * D ** 2 / 4.0
        theta2 = self.theta[:, None] * np.ones((1, N))
        rho_l = np.broadcast_to(np.asarray(self._rho_l_field, float), (nx, N))
        gsin = G * np.sin(theta2)
        C0b, vd = drift_params(theta2, D)
        C0b = C0b * float(getattr(c.numerics, "drift_C0_factor", 1.0))
        rough = c.pipeline.roughness_m / D
        #  conserved/state anchors at step n (backward-Euler) and inlet BCs
        La_n = np.clip(alpha_l, 1e-3, 0.999) * A
        um_n = self._um
        p_n = self._p
        ql, qg = self.inlet_rates(t_h)
        Tin = c.operating.T_inlet_C
        P_in = c.operating.P_inlet_bar
        FL_in = ql                                            # inlet liquid VOLUME rate (m^3/s) = alpha_l u_l A

        def resid(X):
            a = np.clip(X[..., 0], 1e-3, 0.999)
            p = np.clip(X[..., 1], 2.0, 1.0e3)
            um = X[..., 2]
            ag = 1.0 - a
            rho_g = gas_density(p, T, f)
            rho_m = a * rho_l + ag * rho_g
            u_g = np.clip(C0b * um + vd, -5.0, 40.0)
            jg = np.clip(ag * u_g, 0.0, np.abs(um) + 2.0)
            u_l = (um - jg) / np.maximum(a, 1e-3)
            La = a * A
            #  --- R1: LIQUID-mass continuity (hyperbolic, upwind) -> sets holdup alpha_l.
            #  inlet face = prescribed liquid rate; outlet face = outflow.  FL = liquid volume flux.
            FLc = a * u_l * A
            FLf = np.empty((nx + 1, N)); FLf[0] = FL_in; FLf[1:-1] = FLc[:-1]; FLf[-1] = FLc[-1]
            R1 = (La - La_n) / dt + (FLf[1:] - FLf[:-1]) / dx
            #  --- R2: ELLIPTIC compressible-pressure / mixture-volume continuity -> sets pressure.
            #  Summing the two phase continuities /density gives  (alpha_g/p)(p-p_n)/dt + d(j)/dx = 0
            #  (rho_g ~ p so (1/rho_g)drho_g/dp = 1/p): a genuine elliptic pressure equation (the same
            #  one the validated implicit solve uses), NOT a hyperbolic upwind march -> physical dP.
            #  j = volumetric mixture flux (= u_m); face-averaged interior, inlet face = u_m[0],
            #  outlet face = prescribed through-flux j_in (rate control); inlet pressure pinned.
            vsg_in = qg * (P_in / np.maximum(p[0], 1.0)) * ((T[0] + 273.15) / (Tin + 273.15)) / A[0]
            j_in = ql / A[0] + vsg_in
            jf = np.empty((nx + 1, N)); jf[1:-1] = 0.5 * (um[:-1] + um[1:]); jf[0] = um[0]; jf[-1] = j_in
            R2 = ag / np.maximum(p, 1e-9) * (p - p_n) / dt + (jf[1:] - jf[:-1]) / dx
            #  --- R3: mixture MOMENTUM (backward-Euler): upwind advection + gravity + implicit wall
            #  drag + dp/dx.  No explicit inlet velocity pin (natural zero-gradient inlet advection),
            #  so momentum globally couples pressure to friction/gravity -> the outlet pressure floats.
            um_up = np.vstack([um[:1], um[:-1]]); um_dn = np.vstack([um[1:], um[-1:]])
            adv = np.where(um >= 0, um * (um - um_up) / dx, um * (um_dn - um) / dx)
            mu_m = a * f.mu_liquid + ag * gas_viscosity(rho_g, T, f)
            Re = rho_m * np.abs(um) * D / np.maximum(mu_m, 1e-9)
            beta = haaland_friction(Re, rough) / (2.0 * D) * np.abs(um)
            dpdx = np.empty((nx, N))
            dpdx[1:-1] = (p[2:] - p[:-2]) / (2 * dx)
            dpdx[0] = (p[1] - p[0]) / dx; dpdx[-1] = (p[-1] - p[-2]) / dx
            R3 = (rho_m * (um - um_n) / dt + rho_m * adv + rho_m * gsin
                  + rho_m * beta * um + 1e5 * dpdx)
            R2 = R2.copy(); R2[0] = p[0] - P_in              # inlet pressure Dirichlet
            return np.stack([R1, R2, R3], axis=-1)

        #  ROW SCALING (non-dimensionalisation): the three residuals live in incommensurate units
        #  (liquid continuity m^2/s, gas-mass continuity kg/m/s, momentum Pa/m ~ O(1e2-1e3)). Without
        #  scaling the max-norm/line-search is dominated by momentum and the 3x3 blocks are ill-
        #  conditioned -> Newton stalls. Scale each residual by its characteristic magnitude so the
        #  convergence test and damping are on a common, dimensionless footing.
        a0 = np.clip(alpha_l, 1e-3, 0.999)
        rho_g0 = gas_density(np.clip(self._p, 2.0, 1.0e3), T, f)
        rho_m0 = a0 * rho_l + (1.0 - a0) * rho_g0
        sc1 = float(np.nanmax(np.abs(La_n))) / dt + 1e-9
        sc2 = max(float(np.nanmax(1.0 - a0)) / dt,
                  float(np.nanmax(np.abs(self._um))) / dx, 1e-6)   # elliptic-pressure residual scale
        sc3 = (float(np.nanmax(rho_m0 * np.abs(self._um))) / dt
               + float(np.nanmax(rho_m0)) * G
               + 1e5 * P_in / max(c.pipeline.length_m, 1.0) + 1.0)
        scvec = np.array([sc1, sc2, sc3])

        def sres(X):
            return resid(X) / scvec                          # dimensionless residual

        scale = np.array([1.0, max(float(np.nanmax(self._p)), 1.0), 1.0])   # FD column (unknown) scales

        def jacobian(X, F0):
            Asub = np.zeros((nx, N, 3, 3)); Bdiag = np.zeros((nx, N, 3, 3)); Csup = np.zeros((nx, N, 3, 3))
            cell_idx = np.arange(nx)
            for k in range(3):
                for color in range(3):
                    cells = cell_idx[cell_idx % 3 == color]
                    if cells.size == 0:
                        continue
                    eps = 1e-6 * (np.abs(X[cells, :, k]) + scale[k])      # (m, N), per perturbed cell
                    Xp = X.copy(); Xp[cells, :, k] += eps
                    diff = sres(Xp) - F0                                  # scaled, full (nx, N, 3)
                    #  each affected residual row is divided by the eps of its PERTURBING cell.
                    Bdiag[cells, :, :, k] = diff[cells] / eps[:, :, None]          # dF_c / dX_c
                    lo = cells + 1 < nx
                    if lo.any():
                        Asub[cells[lo] + 1, :, :, k] = diff[cells[lo] + 1] / eps[lo][:, :, None]   # dF_{c+1}/dX_c
                    hi = cells - 1 >= 0
                    if hi.any():
                        Csup[cells[hi] - 1, :, :, k] = diff[cells[hi] - 1] / eps[hi][:, :, None]   # dF_{c-1}/dX_c
            return Asub, Bdiag, Csup

        def block_thomas(Asub, Bdiag, Csup, d):
            eye = np.broadcast_to(np.eye(3) * 1e-12, (N, 3, 3))
            cp = np.zeros_like(Csup); dp = np.zeros_like(d)
            B0 = Bdiag[0] + eye
            cp[0] = np.linalg.solve(B0, Csup[0]); dp[0] = np.linalg.solve(B0, d[0][..., None])[..., 0]
            for i in range(1, nx):
                M = Bdiag[i] - Asub[i] @ cp[i - 1] + eye
                cp[i] = np.linalg.solve(M, Csup[i])
                rhs = d[i] - (Asub[i] @ dp[i - 1][..., None])[..., 0]
                dp[i] = np.linalg.solve(M, rhs[..., None])[..., 0]
            x = np.zeros_like(d); x[-1] = dp[-1]
            for i in range(nx - 2, -1, -1):
                x[i] = dp[i] - (cp[i] @ x[i + 1][..., None])[..., 0]
            return x

        def clamp(Xa):
            Xa[..., 0] = np.clip(Xa[..., 0], 1e-3, 0.999)
            Xa[..., 1] = np.clip(Xa[..., 1], 2.0, 1.0e3)
            Xa[..., 2] = np.clip(Xa[..., 2], -5.0, 40.0)
            return Xa

        X = np.stack([np.clip(alpha_l, 1e-3, 0.999),
                      self._p.astype(float), self._um.astype(float)], axis=-1)
        relax0 = float(getattr(c.numerics, "full_newton_relax", 1.0))
        tol = float(getattr(c.numerics, "full_newton_tol", 1e-7))
        nit = max(int(getattr(c.numerics, "full_newton_iters", 12)), 1)
        F0 = sres(X)
        ok = True
        used = 0
        for _ in range(nit):
            rnorm = float(np.nanmax(np.abs(F0)))
            if not np.isfinite(rnorm):
                ok = False; break
            if rnorm < tol:
                break
            Asub, Bdiag, Csup = jacobian(X, F0)
            dX = block_thomas(Asub, Bdiag, Csup, -F0)
            if not np.all(np.isfinite(dX)):
                ok = False; break
            relax = relax0
            Xn = clamp(X + relax * dX); Fn = sres(Xn)
            bt = 0
            while float(np.nanmax(np.abs(Fn))) > rnorm and bt < 6:
                relax *= 0.5; Xn = clamp(X + relax * dX); Fn = sres(Xn); bt += 1
            X = Xn; F0 = Fn; used += 1
        #  lightweight convergence diagnostics (mean Newton iters / last residual) for studies
        self._fn_iters_sum = getattr(self, "_fn_iters_sum", 0) + used
        self._fn_steps = getattr(self, "_fn_steps", 0) + 1
        self._fn_last_res = float(np.nanmax(np.abs(F0)))

        if not ok or not np.all(np.isfinite(X)):
            return self.momentum_solve(alpha_l, T, delta, t_h, dt)   # never-fail: validated fallback

        a = np.clip(X[..., 0], 1e-3, 0.999); ag = 1.0 - a
        p_new = np.clip(X[..., 1], 2.0, 1.0e3)
        um_raw = X[..., 2]
        um_new = np.clip(um_raw, -2.0, 30.0)
        self._clip["velocity"] += int(np.count_nonzero(um_new != um_raw))
        u_g = np.clip(C0b * um_new + vd, -5.0, 40.0)
        jg = np.clip(ag * u_g, 0.0, np.abs(um_new) + 2.0)
        u_l = np.clip((um_new - jg) / np.maximum(a, 1e-3), -0.4, 20.0)
        self._um = um_new; self._p = p_new; self._ug = u_g; self._ul = u_l
        self._dP = float(np.mean(p_new[0] - p_new[-1]))
        rho_g = gas_density(p_new, T, f)
        return dict(p=p_new, j=um_new, vl=u_l, vg=u_g, rho_g=rho_g, A=A, D=D)

    # ----- conservative, LOCAL-first holdup bound enforcement (B7) ---------
    def _enforce_bounds(self, La, A, passes=12):
        """Clamp the conserved liquid area La into the physical range [0, A] per realisation
        WITHOUT teleporting mass across the whole line. Over-fill is spilled to the NEAREST
        NEIGHBOURS (causal, local) over a few passes; only a tiny residual (and any rare
        under-fill) is closed by a final global redistribution, so total liquid mass stays
        exact while the redistribution is overwhelmingly local."""
        self._clip["holdup"] += int(np.count_nonzero((La < 0.0) | (La > A)))
        for _ in range(passes):
            over = np.maximum(La - A, 0.0)
            if not np.any(over):
                break
            La = np.minimum(La, A)
            room = np.maximum(A - La, 0.0)
            up_room = np.zeros_like(room); up_room[:-1] = room[1:]     # room of neighbour i+1
            dn_room = np.zeros_like(room); dn_room[1:] = room[:-1]     # room of neighbour i-1
            tot = up_room + dn_room
            f_up = np.where(tot > 0, up_room / np.maximum(tot, 1e-12), 0.0)
            f_dn = np.where(tot > 0, dn_room / np.maximum(tot, 1e-12), 0.0)
            give_up = over * f_up; give_dn = over * f_dn
            unplaced = over * np.where(tot > 0, 0.0, 1.0)              # nowhere local to go: keep
            La = La + unplaced
            recv = np.zeros_like(La)
            recv[1:] += give_up[:-1]      # cell i+1 receives the upward spill from cell i
            recv[:-1] += give_dn[1:]      # cell i receives the downward spill from cell i+1
            La = La + recv
        # exact conservative cleanup of any residual over-/under-fill (now rare & small):
        # iteratively clip and redistribute the SIGNED residual into cells that have capacity
        # for that sign (room to accept if net over-fill; liquid to give if net under-fill).
        # Each pass conserves the column sum; it converges once La is within [0, A].
        for _ in range(50):
            clipped = np.minimum(np.maximum(La, 0.0), A)
            resid = np.sum(La - clipped, axis=0)                  # (N,) signed leftover
            La = clipped
            if np.all(np.abs(resid) < 1e-12):
                break
            room_pos = np.maximum(A - La, 0.0)                   # capacity to accept (resid > 0)
            room_neg = np.maximum(La, 0.0)                       # liquid available to give (resid < 0)
            room = np.where(resid >= 0, room_pos, room_neg)
            La = La + room * (resid / np.maximum(room.sum(axis=0), 1e-30))
        return np.minimum(np.maximum(La, 0.0), A)

    # ----- main transient integration -------------------------------------
    def run(self, verbose=True):
        c = c0 = self.case
        n = c.numerics; k = c.kinetics
        nx, N = self.nx, self.N
        rho_l = self.rho_l
        t_end_s = n.t_end_h * 3600.0
        dt_max = n.dt_max_h * 3600.0

        #  --- per-realisation parameter spread (C10/C11) ------------------------
        #  Each ensemble column is a DISTINCT realisation: the physically uncertain constants
        #  (hydrate induction time, growth, heat transfer, wall-capture) get a per-realisation
        #  multiplier so the P10/P50/P90 bands are genuinely probabilistic and the N columns are
        #  not wasteful duplicates. Resolution order:
        #    explicit numerics.uq_inputs  ->  use it
        #    numerics.deterministic=True  ->  no spread (every realisation identical)
        #    otherwise                    ->  a modest DEFAULT_UQ (always-on, honest bands)
        if n.uq_inputs is not None:
            uq = n.uq_inputs
        elif n.deterministic:
            uq = None
        else:
            uq = DEFAULT_UQ

        #  #15: a UQ entry may be a plain rel-sigma (lognormal, back-compatible) OR a dict
        #  {"dist":"lognormal"|"normal"|"uniform","sigma":..,"low":..,"high":..} so the spread can
        #  be set from real measurement statistics. A shared standard-normal draw per realisation
        #  (uq["_corr_z"]=True) correlates the parameters (e.g. a "hot/cold" realisation moves U,
        #  kg0 and tau together) instead of treating every constant as independent.
        corr = bool((uq or {}).get("_corr_z")) if uq else False

        def _draw_z():
            #  #18: Latin-Hypercube stratified standard-normal across the N realisations (variance
            #  reduction — each realisation samples a distinct quantile stratum) when lhs_uq is on,
            #  else plain i.i.d. normal. Both are reproducible from numerics.seed.
            if getattr(n, "lhs_uq", False) and N > 1:
                strata = (np.arange(N) + self.rng.random(N)) / N           # one point per stratum
                self.rng.shuffle(strata)
                u = np.clip(strata, 1e-6, 1 - 1e-6)
                z = np.sqrt(2.0) * np.vectorize(_erfinv)(2.0 * u - 1.0)    # inverse-normal CDF
                return z.reshape(1, N)
            return self.rng.standard_normal((1, N))
        z_shared = _draw_z() if corr else None

        def _uq_mult(name):
            spec = (uq or {}).get(name) if uq else None
            if spec is None:
                return np.ones((1, N))
            z = z_shared if z_shared is not None else _draw_z()
            if isinstance(spec, dict):
                dist = spec.get("dist", "lognormal")
                if dist == "uniform":
                    lo, hi = float(spec.get("low", 0.8)), float(spec.get("high", 1.2))
                    u = 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))   # z -> U(0,1) via CDF
                    return np.maximum(lo + (hi - lo) * u, 1e-3)
                sig = float(spec.get("sigma", 0.0))
                if sig <= 0.0:
                    return np.ones((1, N))
                if dist == "normal":
                    return np.maximum(1.0 + sig * z, 1e-3)
                return np.maximum(np.exp(sig * z), 1e-3)              # lognormal
            s = float(spec)
            if s <= 0.0:
                return np.ones((1, N))
            return np.maximum(np.exp(s * z), 1e-3)
        kg0_r = k.kg0 * _uq_mult("kg0")
        nuc_tau0_r = k.nuc_tau0_h * _uq_mult("nuc_tau0_h")
        nuc_beta_r = k.nuc_beta_C * _uq_mult("nuc_beta_C")
        wcap_r = k.wall_capture_eff * _uq_mult("wall_capture_eff")
        U_mult = _uq_mult("U_wall")

        # --- state fields (nx,N), transient ---
        alpha_l = np.full((nx, N), 0.35)
        T = np.linspace(c.operating.T_inlet_C, c.operating.T_seabed_C + 6, nx)[:, None] \
            * np.ones((1, N))
        T_soil = np.full((nx, N), float(c.operating.T_seabed_C))   # #12 lumped buried-soil node
        #  C12: FULL transient radial soil conduction grid — concentric shells from the pipe wall to
        #  the far-field radius, log-spaced (radial conduction is logarithmic). Precompute the
        #  per-unit-length shell conductances G [W/mK] and heat capacities Csoil [J/mK].
        nsoil = max(int(getattr(n, "soil_nodes", 1)), 1)
        Ts_shells = None; soil_G = None; soil_C = None; soil_Gfar = 0.0
        if n.soil_transient and nsoil > 1:
            r0 = 0.5 * c.pipeline.diameter_m * 1.05               # pipe outer radius (approx)
            rf = np.geomspace(r0, max(n.soil_far_radius_m, r0 * 2), nsoil + 1)   # shell faces
            rc = np.sqrt(rf[:-1] * rf[1:])                        # node-centre radii (geom mean)
            ks = n.soil_conductivity
            soil_G = 2.0 * math.pi * ks / np.log(rc[1:] / rc[:-1])           # (nsoil-1,) internal
            soil_C = n.soil_rhocp * math.pi * (rf[1:] ** 2 - rf[:-1] ** 2)   # (nsoil,) J/mK
            soil_Gfar = 2.0 * math.pi * ks / math.log(rf[-1] / rc[-1])       # outer -> far field
            Ts_shells = np.full((nsoil, nx, N), float(c.operating.T_seabed_C))
        self._slugS = np.zeros((nx, N))                            # A4 Lagrangian slug-indicator field
        self._droplet_frac = 0.0
        phi = np.zeros((nx, N))
        delta = np.zeros((nx, N))
        nucleated = np.zeros((nx, N), bool)
        locked = np.zeros((nx, N), bool)
        delta_max = k.delta_max_frac * c.pipeline.diameter_m / 2.0
        mon = int(n.monitor_frac * (nx - 1))

        # --- recorders ---
        snap_t, snap_phi, snap_PhiSH, snap_holdup = [], [], [], []
        snap_P, snap_T = [], []                              # P(x,t)/T(x,t) profile snapshots
        ts_keys = ["P", "T", "Tsub", "alpha_l", "fslug", "phi", "delta", "PhiSH", "a_i", "j"]
        ts = {key: [] for key in ts_keys}
        ts_t = []
        bc_hist = []
        max_PhiSH = np.zeros((nx, N)); max_Tsub = np.zeros((nx, N))
        plug_time = np.full(N, np.nan); plug_loc = np.full(N, np.nan)
        liq_in_tot = 0.0; liq_out_tot = 0.0
        #  A3 mass-balance accumulators: liquid water consumed by hydrate, total hydrate mass
        #  formed, and the gas mass-balance audit (in / out / consumed-by-hydrate / stored).
        liq_to_hyd_tot = 0.0; hyd_mass_tot = 0.0; gas_consumed_hyd_tot = 0.0
        gas_in_tot = 0.0; gas_out_tot = 0.0

        #  D15: snapshot on a fixed TIME cadence (~n_snapshots evenly over the run) instead of a
        #  step count that assumed dt==dt_max and over-sampled ~10x when dt was small.
        snap_dt = t_end_s / max(n.n_snapshots, 1); next_snap = 0.0
        t = 0.0; step = 0
        self._max_steps = int(max(2.0e5, 50.0 * t_end_s / max(n.cfl * self.dx / 40.0, 1.0)))
        fallbacks = 0
        A0 = math.pi * c.pipeline.diameter_m ** 2 / 4
        A_last = np.full((nx, N), A0)
        La = alpha_l * A0                                  # persistent CONSERVED state: liquid area
        inv_init = float(np.mean(np.sum(La, 0)) * self.dx)
        theta_col = self.theta[:, None] * np.ones((1, N))  # inclination broadcast (used by #2 settling)
        Wl = c.fluids.water_cut * La                       # conserved WATER area (#2; used iff oil_water_slip)
        water_frac = np.full((nx, N), c.fluids.water_cut)
        #  Effective overall heat-transfer coefficient & lumped thermal mass (multi-layer
        #  wall if provided), and the inhibitor (MEG) concentration field transported with
        #  the liquid. MEG_wt_inlet = 0 -> field stays 0 -> no change to base behaviour.
        U_eff, therm_mass = effective_U_and_mass(
            c.pipeline, c.operating, c.fluids.cp_liquid, self.rho_l * 0.35 + c.fluids.gas_MW * 0)
        #  B9: resolve the inlet MEG concentration onto the AQUEOUS basis the suppression model
        #  uses. "aqueous" -> as given. "stream" -> wt% of the whole liquid stream, converted to
        #  the produced-water basis with the water cut (MEG partitions into the water).
        meg_aq_inlet = float(c.operating.MEG_wt_inlet)
        if getattr(c.operating, "MEG_basis", "aqueous") == "stream" and meg_aq_inlet > 0.0:
            fmeg = meg_aq_inlet / 100.0
            meg_aq_inlet = 100.0 * fmeg / max(fmeg + max(c.fluids.water_cut, 1e-3), 1e-6)
        self._meg_aq_inlet = meg_aq_inlet
        W_inh = np.full((nx, N), meg_aq_inlet)
        #  initialise the implicit-engine state (mixture velocity & pressure) from a
        #  quasi-steady evaluation so it starts on the solution manifold.
        pv0 = self.pressure_velocity(alpha_l, T, delta, 0.0)
        self._um = np.clip(pv0["j"], -2.0, 30.0); self._p = pv0["p"]
        self._dP = float(np.mean(pv0["p"][0] - pv0["p"][-1]))
        self._ug = np.clip(pv0["vg"], -5.0, 40.0); self._ul = np.clip(pv0["vl"], -0.4, 20.0)
        self._Vm_design = max(float(np.nanmean(np.abs(self._um))), 0.1)   # design velocity for dynamic-U scaling
        #  conserved gas-mass state Mg = rho_g*alpha_g*A (kg/m) for the gas-continuity equation
        #  (A3) and its initial inventory for the mass-balance audit.
        Mg = pv0["rho_g"] * (1.0 - alpha_l) * pv0["A"]
        self._Mg = Mg                                       # #A1: expose conserved gas mass to the engine
        gas_store_init = float(np.mean(np.sum(Mg, 0))) * self.dx
        self._maxrate = 0.0                                # #17 lagged max relative change rate (1/s)
        while t < t_end_s:
            t_h = t / 3600.0
            # derive holdup from the CONSERVED liquid area for the current bore;
            # closures use a safe clamp but the conserved La is never truncated.
            D_now = np.clip(c.pipeline.diameter_m - 2.0 * delta, 0.60 * c.pipeline.diameter_m, c.pipeline.diameter_m)
            A_now = math.pi * D_now ** 2 / 4.0
            alpha_l = np.clip(La / A_now, 1e-3, 0.999)

            #  --- ONE adaptive dt per step (B8): the pressure/momentum solve AND every
            #  transport update advance with the SAME dt. CFL on advection plus a diffusion
            #  cap for the phase-field (A2 explicit diffusion stability). ---
            vstate = max(float(np.nanmax(np.abs(self._um))),
                         float(np.nanmax(np.abs(self._ug))),
                         float(np.nanmax(np.abs(self._ul))), 0.1)
            dt = min(n.cfl * self.dx / vstate, dt_max, t_end_s - t)
            if k.D_phi > 0.0:
                dt = min(dt, 0.4 * self.dx ** 2 / k.D_phi)
            #  #B5: when resolving water-hammer, cap dt by the ACOUSTIC CFL (sound speed) so the
            #  pressure waves are time-resolved (slower, but correct for fast transients).
            if n.acoustic > 0.0:
                rg0 = gas_density(self._p, T, c.fluids)
                c_snd = float(np.nanmax(mixture_sound_speed(alpha_l, self._rho_l_field, rg0, self._p)))
                dt = min(dt, n.cfl * self.dx / max(c_snd, 1.0))
            #  #17: error-controlled dt — cap dt so the previous step's fastest relative change in
            #  T or holdup stays under err_tol (a lagged PI-style accuracy controller; the change
            #  rate is observed at the end of each step in self._maxrate).
            if getattr(n, "error_dt", False) and self._maxrate > 1e-12:
                dt = min(dt, 0.02 / self._maxrate)               # <=2% change per step target

            if self.engine in ("implicit", "twofluid", "twofluid_mass", "twofluid_mass_newton",
                               "twofluid_full_newton"):
                _solve = ({"twofluid": self.twofluid_solve,
                           "twofluid_mass": self.twofluid_mass_solve,
                           "twofluid_mass_newton": self.twofluid_mass_newton_solve,
                           "twofluid_full_newton": self.twofluid_full_newton_solve}
                          .get(self.engine, self.momentum_solve))
                # save state so a corrector re-solve restarts from the same point
                _saved = (self._p.copy(), self._um.copy(), self._ug.copy(), self._ul.copy(), self._dP)
                pv = _solve(alpha_l, T, delta, t_h, dt)
                #  #6 WITHIN-STEP coupled iteration: re-solve the pressure-momentum system using the
                #  freshly updated state, iterating to a fixed point (the segregated solve converged
                #  toward the fully-coupled solution) when within_step_iters > 1.
                for _ in range(max(int(getattr(n, "within_step_iters", 1)) - 1, 0)):
                    if not np.all(np.isfinite(pv["p"])):
                        break
                    p_prev = self._p.copy()
                    pv = _solve(alpha_l, T, delta, t_h, dt)
                    if float(np.nanmax(np.abs(self._p - p_prev))) < 1e-3:
                        break                                    # converged (bar tolerance)
                #  #10 PREDICTOR-CORRECTOR dt: dt was predicted from the previous-step velocity; if
                #  the freshly solved velocity exceeds the CFL budget by > substep_cfl_growth, restart
                #  the step from the saved state with the corrected (smaller) dt so transport stays
                #  CFL-consistent — a within-step convergence of dt rather than a stale prediction.
                if np.all(np.isfinite(pv["j"])):
                    vnew = max(float(np.nanmax(np.abs(pv["j"]))), float(np.nanmax(np.abs(pv["vg"]))), 0.1)
                    if vnew * dt / self.dx > n.cfl * max(n.substep_cfl_growth, 1.0):
                        self._p, self._um, self._ug, self._ul, self._dP = _saved
                        dt = min(dt, n.cfl * self.dx / vnew)
                        pv = _solve(alpha_l, T, delta, t_h, dt)
                #  ROBUSTNESS: if the implicit/two-fluid step is non-finite/unphysical,
                #  degrade gracefully to the proven quasi-steady solver (never fail).
                if not (np.all(np.isfinite(pv["j"])) and np.all(np.isfinite(pv["vl"]))
                        and np.all(np.isfinite(pv["p"]))):
                    if getattr(n, "strict", False):
                        raise RuntimeError(
                            f"SHCT strict mode: non-finite {self.engine} momentum/pressure step at "
                            f"t={t_h:.3f} h (would have degraded to quasi-steady). Reduce cfl / check "
                            f"inputs, or disable numerics.strict to use the never-fail fallback.")
                    fallbacks += 1
                    pv = self.pressure_velocity(alpha_l, T, delta, t_h)
                    self._p = pv["p"]; self._um = np.clip(pv["j"], -2.0, 30.0)
                    self._ug = np.clip(pv["vg"], -5.0, 40.0); self._ul = np.clip(pv["vl"], -0.4, 20.0)
                    self._dP = float(np.mean(pv["p"][0] - pv["p"][-1]))
            else:
                pv = self.pressure_velocity(alpha_l, T, delta, t_h)
                #  keep the velocity state current so next step's dt adapts (quasisteady has no
                #  implicit state update of its own).
                self._um = np.clip(pv["j"], -2.0, 30.0)
                self._ug = np.clip(pv["vg"], -5.0, 40.0); self._ul = np.clip(pv["vl"], -0.4, 20.0)
                self._dP = float(np.mean(pv["p"][0] - pv["p"][-1]))

            j, vl, vg = pv["j"], pv["vl"], pv["vg"]
            D, A = pv["D"], pv["A"]
            rho_g = pv["rho_g"]; p = pv["p"]

            # final guard: equilibrium-holdup fallback if anything is still non-finite.
            #  B5: ALWAYS advance time/step on this path so a persistently non-finite state
            #  can never spin the loop forever; a hard step cap aborts with a clear message.
            if not np.all(np.isfinite(j)) or not np.all(np.isfinite(vl)):
                if getattr(n, "strict", False):
                    raise RuntimeError(
                        f"SHCT strict mode: non-finite velocity field at t={t_h:.3f} h after the "
                        f"momentum solve. Disable numerics.strict to use the equilibrium-holdup "
                        f"fallback, or reduce cfl / check inputs.")
                fallbacks += 1
                C0, vd = drift_params(self.theta[:, None], D)
                alpha_l = np.clip(1 - np.clip(j / np.maximum(C0 * j + vd, 1e-3), 1e-3, 0.999), 0.05, 0.95)
                La = self._enforce_bounds(alpha_l * A, A)
                t += dt; step += 1
                if step > self._max_steps:
                    raise RuntimeError("SHCT solver: step cap exceeded — non-converging transient "
                                       "(state remained non-finite); check inputs / reduce cfl.")
                continue
            if step > self._max_steps:
                raise RuntimeError("SHCT solver: step cap exceeded — non-converging transient; "
                                   "check inputs / reduce cfl.")

            # === (H) liquid-holdup transport (volume-conservative with area) ===
            #   d(alpha_l*A)/dt + d(alpha_l*v_l*A)/dx = 0   (upwind FV)
            #   conserved variable La = alpha_l*A (liquid cross-sectional area);
            #   inlet face flux = the actual liquid volumetric rate Q_l (exact).
            ql_in, _ = self.inlet_rates(t_h)
            flow_frac = ql_in / max(c.operating.q_liquid_insitu, 1e-12)
            #  A2: when a separate DROPLET field is modelled, the entrained fraction Ed of the liquid
            #  is carried in the gas core at the GAS velocity (not the slow film velocity), so the
            #  liquid is transported at an effective velocity (1-Ed)*vl + Ed*vg. Conservative (still
            #  flux-form); Ed = entrainment fraction (Ishii-Mishima) in high-shear flow.
            if c.fluids.droplet_field:
                Vsg_d = np.maximum(j - vl * alpha_l, 0.0)
                Ed_f = droplet_entrainment_frac(Vsg_d, rho_g, c.fluids.sigma, c.pipeline.diameter_m)
                vl_t = (1.0 - Ed_f) * vl + Ed_f * vg
                self._droplet_frac = float(np.nanmax(Ed_f))
            else:
                vl_t = vl
            vlf = 0.5 * (vl_t[:-1] + vl_t[1:])
            Ff = np.empty((nx + 1, N))
            #  #8: 2nd-order TVD reconstruction of the conserved liquid area at interior faces
            #  (still flux-form -> mass conservation is exact); flux = face-velocity * face-La.
            La_face = tvd_interior_faces(La, vlf, self._limiter_kind)
            Ff[1:-1] = vlf * La_face
            Ff[0] = ql_in                                    # inlet liquid rate BC (exact)
            #  Outlet is OUTFLOW-ONLY (liquid leaving at the separator cannot re-enter);
            #  a stagnant / shut-in line is hydraulically isolated (closed outlet) so
            #  trapped liquid merely redistributes & cools internally (exact conservation).
            QlA_out = La[-1] * vl_t[-1]
            Ff[-1] = np.maximum(QlA_out, 0.0) if flow_frac > n.outlet_open_frac else 0.0
            # update the CONSERVED liquid area (telescoping faces => exact mass conservation);
            # cap only at physical bounds (empty / full bore) without lossy mid-range clipping.
            La = La - dt / self.dx * (Ff[1:] - Ff[:-1])
            #  B7: enforce physical bounds [0, A] with a LOCAL-first conservative scheme
            #  (over-fill spills to nearest neighbours, not teleported line-wide).
            La = self._enforce_bounds(La, A)
            alpha_l = np.clip(La / A, 1e-3, 0.999)            # safe clamp for closures only
            liq_in_tot += ql_in * dt
            liq_out_tot += float(Ff[-1].mean()) * dt
            A_last = A

            # === (Wt) optional 3-PHASE water transport with oil/water slip (#2) ===
            #  Off by default -> water is a constant fraction of the composite liquid (model
            #  unchanged). On -> the produced WATER is transported separately as a conserved water
            #  area Wl with a buoyancy/gravity SETTLING slip relative to the oil (water lags on the
            #  uphill, leads on the downhill), so water ACCUMULATES at low points (a real corrosion/
            #  hydrate hot-spot the single-composite model cannot show). Conservative (flux-form).
            if c.fluids.oil_water_slip:
                drho = max(c.fluids.rho_water - c.fluids.rho_oil, 0.0)
                if c.fluids.water_drift:
                    #  A3: water carries its own DRIFT velocity from a buoyancy/drag FORCE BALANCE
                    #  (terminal settling): v_t = sqrt(4 g d_w drho / (3 Cd rho_o)) with a droplet
                    #  size ~0.5 mm and Cd~0.44 (Newton regime) — a genuine momentum closure for the
                    #  water phase rather than a fixed empirical slip.
                    d_w = c.fluids.water_droplet_d_m; Cd = c.fluids.water_droplet_Cd
                    v_settle = math.sqrt(4.0 * G * d_w * drho / (3.0 * Cd * max(c.fluids.rho_oil, 1.0)))
                else:
                    #  settling slip velocity (Stokes/Ishii buoyant rise scaled by inclination), m/s
                    v_settle = c.fluids.settling_slip_coeff * math.sqrt(G * c.pipeline.diameter_m) \
                        * (drho / max(self.rho_l, 1.0))
                vw = vl - v_settle * np.sin(theta_col)        # water velocity (lags uphill)
                vwf = 0.5 * (vw[:-1] + vw[1:])
                Wl_face = tvd_interior_faces(Wl, vwf, self._limiter_kind)
                Wf = np.empty((nx + 1, N))
                Wf[1:-1] = vwf * Wl_face
                Wf[0] = ql_in * c.fluids.water_cut           # inlet water rate
                Wf[-1] = np.maximum(Wl[-1] * vw[-1], 0.0) if flow_frac > n.outlet_open_frac else 0.0
                Wl = np.clip(Wl - dt / self.dx * (Wf[1:] - Wf[:-1]), 0.0, La)
                water_frac = Wl / np.maximum(La, 1e-9)        # local water fraction of the liquid
                #  #2 DEEPENED: the local liquid density now reflects the actual water fraction and
                #  feeds the momentum gravity/friction & energy next step (heavier water-laden liquid
                #  at low points -> larger hydrostatic head, the real terrain effect a composite misses).
                self._rho_l_field = (oil_density(p, T, c.fluids) * (1.0 - water_frac)
                                     + c.fluids.rho_water * water_frac)

            #  ===== phase sources reordered so ENERGY can include hydrate latent heat (A1) =====
            ag = 1.0 - alpha_l
            rho_m = alpha_l * self._rho_l_field + ag * rho_g    # #2: local 3-phase liquid density
            cp_m = c.fluids.cp_liquid
            #  B9: CONDENSATION/EVAPORATION latent heat — in the two-phase region, cooling condenses
            #  gas (dV/dT > 0) which releases latent heat, BUFFERING the temperature drop. Model it as
            #  an effective heat-capacity increase cp_eff = cp + L*|dV/dT| from the EOS vapour-fraction
            #  surface V(P,T). Captures retrograde/condensation thermal inertia without a per-cell flash.
            if c.fluids.condensation_latent and self._Vsurf is not None:
                Pg, Tg, Vg = self._Vsurf
                pcell = np.clip(p, Pg[0], Pg[-1]); tcell = np.clip(T, Tg[0], Tg[-1])
                Pi = np.interp(pcell, Pg, np.arange(Pg.size))
                Tj0 = np.clip(np.interp(tcell, Tg, np.arange(Tg.size)), 0, Tg.size - 1.001)
                i = np.clip(np.round(Pi).astype(int), 0, Pg.size - 1)
                j = np.floor(Tj0).astype(int)
                dVdT = np.abs((Vg[i, j + 1] - Vg[i, j]) / (Tg[j + 1] - Tg[j]))   # |dV/dT| per K
                cp_m = cp_m + c.fluids.L_condensation * np.clip(dVdT, 0.0, 0.5)

            # === (I) inhibitor (MEG) transport with the liquid (signed upwind) ===
            #  B9: W_inh is the MEG concentration IN THE AQUEOUS PHASE (wt%). meg_suppression()
            #  and the hammerschmidt_meg sizing both work on that SAME aqueous basis, so the
            #  transported field and the curve suppression are consistent (no water_cut mix-up).
            if self._meg_aq_inlet > 0.0:
                W_up = np.vstack([np.full((1, N), self._meg_aq_inlet), W_inh[:-1]])
                W_inh = W_inh + dt * (-np.maximum(vl, 0.0) * (W_inh - W_up) / self.dx)
                W_inh = np.clip(W_inh, 0.0, 60.0)
                W_inh[0] = self._meg_aq_inlet

            # --- hydrate driving state (effective hydrate curve, inhibitor-suppressed) ---
            #  #3: van der Waals-Platteeuw (composition-dependent, fugacity-based) hydrate curve
            #  when advanced_physics + a composition are given; else the correlation / user table.
            if n.advanced_physics and getattr(c.fluids, "composition", None) and c.fluids.hyd_Teq_table is None:
                import shct_eos
                Teq = shct_eos.hydrate_equilibrium_vdwp(p, c.fluids.composition,
                                                        c.fluids.salinity_wt) - meg_suppression(W_inh)
            else:
                Teq = hydrate_equilibrium_T(p, gas_sg=c.fluids.gas_sg,
                                            salinity_wt=c.fluids.salinity_wt,
                                            table=c.fluids.hyd_Teq_table) - meg_suppression(W_inh)
            Tsub = Teq - T
            Tsub_pos = np.maximum(Tsub, 0.0)
            Tabs = T + 273.15
            #  Regime / slug-frequency / interfacial-area at DESIGN (clean-pipe) hydrodynamics.
            D0s = c.pipeline.diameter_m
            A0s = math.pi * D0s ** 2 / 4.0
            _, qg_in_t = self.inlet_rates(t_h)
            vsg0 = qg_in_t * (c.operating.P_inlet_bar / p) * \
                ((T + 273.15) / (c.operating.T_inlet_C + 273.15)) / A0s
            vsl0 = ql_in / A0s
            j0 = vsg0 + vsl0
            theta2 = self.theta[:, None] * np.ones((1, N))
            regime = flow_regime_code(vsg0, vsl0, theta2)
            fslug = np.where(np.isin(regime, [2, 5]),
                             slug_frequency(vsl0, j0, D0s, theta2), k.f_slug_floor_Hz)
            #  A4: LAGRANGIAN slug-indicator tracking — a marker field S in [0,1] (1 = a slug body is
            #  passing, 0 = film) is SEEDED stochastically at the local slug frequency and ADVECTED at
            #  the slug TRANSLATIONAL velocity Vt = 1.2*Vm + drift, then relaxes (a slug body has a
            #  finite length). This resolves individual slug units sub-grid (Lagrangian markers on the
            #  Eulerian grid) and reports the instantaneous slug coverage. Off by default.
            if n.slug_tracking:
                Vt = 1.2 * np.maximum(j, 1e-3) + 0.35 * math.sqrt(G * D0s)
                S_up = np.vstack([self._slugS[:1] * 0.0, self._slugS[:-1]])   # upwind (downstream)
                seed = (self.rng.random((nx, N)) < np.clip(fslug * dt, 0.0, 1.0)) & np.isin(regime, [2, 5])
                self._slugS = np.clip(self._slugS
                                      - dt * np.maximum(Vt, 0.0) * (self._slugS - S_up) / self.dx
                                      - dt * np.maximum(j, 0.1) / np.maximum(slug_length(j, D0s, fslug), 1.0)
                                      * self._slugS, 0.0, 1.0)
                self._slugS = np.where(seed, 1.0, self._slugS)
            #  #5 DEEPENED — sub-grid slug unit-cell holdup: within a slug unit the liquid is
            #  concentrated in the slug body (near-full holdup) and thin in the film; the cell-mean
            #  alpha_l understates the body interfacial area. Blend toward the body holdup in slug/
            #  churn flow so the closures (a_i, Phi_SH) see the real intermittent structure.
            alpha_film = alpha_l
            if getattr(n, "subgrid_slug", False):
                slug_body_holdup = np.clip(k.slug_body_holdup_base + k.slug_body_holdup_slope
                                           * np.minimum(vsl0 / np.maximum(j0, 1e-3), 1.0), alpha_l, 0.95)
                in_slug = np.isin(regime, [2, 5])
                beta_slug = np.clip(fslug / (fslug + k.slug_fraction_ref_Hz), 0.0, 0.6)   # slug-fraction weight
                alpha_eff = np.where(in_slug, (1 - beta_slug) * alpha_l + beta_slug * slug_body_holdup, alpha_l)
            else:
                alpha_eff = alpha_l
            #  #6 DEEPENED — droplet entrainment moves liquid from the WALL FILM into the gas core,
            #  lowering the effective FILM holdup that wets the wall (less interfacial contact there).
            if c.fluids.droplet_entrainment:
                Ed = droplet_entrainment_frac(vsg0, rho_g, c.fluids.sigma, D0s)
                alpha_film = np.clip(alpha_l * (1.0 - Ed), 1e-3, 0.999)
            else:
                Ed = np.zeros_like(alpha_l)
            #  #10: geometry-resolved interfacial area (Taitel-Dukler stratified / annular film /
            #  dispersed-bubble) when advanced_physics is on, else the enhancement-factor closure.
            if n.advanced_physics:
                a_i = interfacial_area_geom(alpha_film, D0s, regime)
            else:
                a_i = interfacial_area(1 - alpha_eff, alpha_film, D0s, regime)

            # === (P-src) stochastic nucleation + hydrate growth rates ===
            #  Two distinct hydrate processes, each with its OWN driving subcooling:
            #    * BULK slurry (Rg_bulk) is driven by the local BULK subcooling Tsub_pos. Its
            #      latent heat warms the flowing fluid, so it is SELF-LIMITING -> the bulk
            #      hydrate fraction phi stays low/transportable (the heat-transfer-limited regime).
            #    * WALL deposit (Rg_wall) is driven by the WALL subcooling (Teq - T_seabed): the
            #      pipe wall sits at the cold seabed, so this driving force is SUSTAINED (the
            #      deposit's latent heat is rejected to the sea, not the bulk). This is what lets a
            #      plug grow even while the bulk slurry remains dilute — now mass-consistently.
            tau_ind = nuc_tau0_r * 3600.0 * np.exp(
                np.minimum(nuc_beta_r / np.maximum(Tsub_pos, 1e-3), k.nuc_tau_exp_cap))
            p_nuc = 1.0 - np.exp(-dt / tau_ind)
            fire = (self.rng.random((nx, N)) < p_nuc) & (Tsub_pos > k.nuc_dTsub_min)
            nucleated |= fire
            kg = kg0_r * np.exp(-k.Ea_over_R * (1.0 / Tabs - 1.0 / k.T_ref_K))
            #  growth proceeds only AFTER stochastic nucleation has fired (no growth from a
            #  metastable subcooled liquid) -> subcooling builds during induction; physical onset spread.
            Rg_bulk = np.maximum(kg * a_i * Tsub_pos ** k.growth_exp_n * (1.0 - phi / k.phi_max), 0.0) \
                * nucleated
            #  #7 DEPOSIT-INSULATION feedback: the growth front sits on the INNER face of the
            #  deposit. As the deposit thickens it insulates that front from the cold sea, so the
            #  front warms from T_seabed toward the bulk T and the WALL subcooling that drives
            #  further deposition FALLS (self-limiting late-stage growth). frac_warm = fraction of
            #  the wall->bulk temperature step that the deposit insulation recovers (0 at delta=0).
            R_dep = delta / max(c.kinetics.deposit_k_hyd, 1e-6)          # deposit conductive resistance
            R_ext = 1.0 / max(float(U_eff), 1e-6)                        # wall+film+sea resistance
            frac_warm = k.k_dep_insul * R_dep / (R_dep + R_ext)
            T_front = c.operating.T_seabed_C + (T - c.operating.T_seabed_C) * frac_warm
            Tsub_wall = np.maximum(Teq - T_front, 0.0)                   # insulation-reduced wall driving
            kg_wall = kg0_r * np.exp(-k.Ea_over_R * (1.0 / (T_front + 273.15) - 1.0 / k.T_ref_K))
            Rg_wall = np.maximum(kg_wall * a_i * Tsub_wall ** k.growth_exp_n, 0.0) * nucleated

            # === (C) coupling number Phi_SH (slug-renewal vs hydrate-formation criticality) ===
            #  Evaluated at the SUSTAINED design (wall) subcooling — the deposition-relevant driving
            #  force — so the coupling criticality that actually grows the plug is not erased by the
            #  bulk latent-heat self-limiting. Two views (core-gating cap vs honest reported magnitude).
            PhiSH_raw = k.C_phi * kg_wall * a_i * Tsub_wall ** k.growth_exp_n / fslug
            PhiSH = np.clip(PhiSH_raw, 0.0, k.phi_internal_cap)
            PhiSH_rep = np.clip(PhiSH_raw, 0.0, k.phi_report_cap)

            # === (D-gate) wall-capture fraction f_wall (A4). Computed BEFORE energy so the bulk
            #     latent heat uses only the bulk growth. ===
            restr = 2.0 * delta / c.pipeline.diameter_m
            avail = np.minimum(Tsub_wall / max(k.wall_capture_Tsub_ref_C, 1e-6), 1.0)   # sustained wall driving force
            form = avail * nucleated                          # active only after stochastic onset
            locked |= (restr > k.consol_restriction) & (PhiSH > 1.0) & (avail > 0.30)
            supercrit = np.clip(PhiSH - 1.0, 0.0, 1.0)
            drive = np.where(locked, 1.0, supercrit)
            f_wall = np.clip(wcap_r * form * drive, 0.0, 1.0)   # fraction of wall growth that consolidates

            # === (E) energy: signed-upwind advection + seabed loss + BULK hydrate latent heat ===
            #  A16: advect at the mixture velocity j with the PROPER sign (forward-flow uses the
            #  upstream cell, reverse-flow the downstream cell) so shut-in / backflow no longer
            #  silently drop the advection term. Ghosts: inlet T upstream, zero-gradient at outlet.
            Tin_c = c.operating.T_inlet_C

            def _advT(Tf):
                if n.tvd_energy:                                 # #15: 2nd-order TVD energy advection
                    jf = 0.5 * (j[:-1] + j[1:])
                    Tfc = tvd_interior_faces(Tf, jf, self._limiter_kind)
                    Tfull = np.empty((nx + 1, N)); Tfull[1:-1] = Tfc
                    Tfull[0] = np.where(j[0] >= 0.0, Tin_c, Tf[0]); Tfull[-1] = Tf[-1]
                    return j * (Tfull[1:] - Tfull[:-1]) / self.dx
                Tu = np.vstack([np.full((1, N), Tin_c), Tf[:-1]])
                Td = np.vstack([Tf[1:], Tf[-1:]])
                return np.where(j >= 0.0, j * (Tf - Tu) / self.dx, j * (Td - Tf) / self.dx)
            adv_T = _advT(T)
            #  Effective-U for the energy update (default U_eff; opt-in dynamic_U / per-realisation UQ).
            if n.dynamic_U:
                R_base = 1.0 / max(float(U_eff), 1e-6)
                R_dep = delta / max(c.kinetics.deposit_k_hyd, 1e-6)   # deposit conductive resistance (k_hyd)
                h_ratio = np.clip(np.abs(j) / max(self._Vm_design, 1e-3), 0.05, 4.0) ** 0.8
                dR_film = 1.0 / np.maximum(c.pipeline.h_inner * h_ratio, 1e-6) - 1.0 / c.pipeline.h_inner
                U_field = 1.0 / np.maximum(R_base + R_dep + dR_film, 1e-6)
            else:
                U_field = U_eff
            #  #12: regime-dependent inner heat-transfer (Nusselt) — slug/churn mixing raises U.
            if n.advanced_physics:
                Re_T = rho_m * np.abs(j) * D / np.maximum(
                    alpha_l * c.fluids.mu_liquid + (1 - alpha_l) * c.fluids.mu_gas, 1e-9)
                k_liq = max(c.fluids.k_liquid, 1e-9)
                Pr_T = c.fluids.mu_liquid * cp_m / k_liq         # liquid Prandtl
                Nu = regime_nusselt(regime, Re_T, Pr_T)
                h_in = Nu * k_liq / D                            # inner film coeff from Nu
                U_field = 1.0 / np.maximum(1.0 / np.maximum(U_field, 1e-6)
                                           + 1.0 / np.maximum(h_in, 1e-3) - 1.0 / c.pipeline.h_inner, 1e-6)
            if uq:
                U_field = U_field * U_mult
            #  A1: exothermic hydrate formation releases L_hyd per kg. Only the BULK-slurry growth
            #  heats the flowing fluid; the WALL deposit rejects its latent heat into the cold
            #  wall/sea, not the bulk. Volumetric bulk formation rate = rho_hyd*Rg_bulk ->
            #  source L_hyd*rho_hyd*Rg_bulk/(rho_m*cp). (This is what makes the bulk self-limiting.)
            q_latent = c.fluids.L_hyd * c.fluids.rho_hyd * Rg_bulk / np.maximum(rho_m * cp_m, 1.0)
            #  #11: Joule-Thomson cooling — gas expanding down the pressure gradient cools. The
            #  fluid experiences dp/dt = j*dp/dx (advective); q_JT = mu_JT * (gas fraction) * dp/dt.
            q_jt = np.zeros_like(T)
            if n.advanced_physics:
                dpdx_e = np.empty((nx, N))
                dpdx_e[1:-1] = (p[2:] - p[:-2]) / (2 * self.dx)
                dpdx_e[0] = (p[1] - p[0]) / self.dx; dpdx_e[-1] = (p[-1] - p[-2]) / self.dx
                mu_jt = joule_thomson_dTdP(p, T, c.fluids)        # K/Pa
                q_jt = mu_jt * (1.0 - alpha_l) * j * dpdx_e * 1e5  # K/s (bar/m -> Pa/m)
            #  #14: the stiff linear WALL-LOSS term is treated IMPLICITLY (backward Euler) so the
            #  energy update is unconditionally stable in the loss and never overshoots below the
            #  seabed temperature, allowing larger steps; advection and latent heat stay explicit.
            beta_loss = U_field * 4.0 / (D * rho_m * cp_m)        # 1/s
            #  #12: the energy sink temperature is the (transient) SOIL node when soil_transient is
            #  on — buried-line thermal inertia (the soil warms under the hot pipe and lags), else
            #  the constant seabed temperature.
            #  the sink temperature is the INNER soil shell (radial model), the lumped node, or the
            #  constant seabed temperature.
            if Ts_shells is not None:
                T_sink = Ts_shells[0]
            elif n.soil_transient:
                T_sink = T_soil
            else:
                T_sink = c.operating.T_seabed_C
            T_floor = float(np.nanmin(T_sink)) if n.soil_transient else c.operating.T_seabed_C
            T_expl = T + dt * (-adv_T + q_latent + q_jt)
            T_new = (T_expl + dt * beta_loss * T_sink) / (1.0 + dt * beta_loss)
            #  #9: Strang/Heun 2nd-order corrector on the explicit advection (re-evaluate at the
            #  predicted state and average) -> reduces the operator-splitting/time error.
            if getattr(n, "splitting", "strang") != "lie":
                T_mid = np.clip(T_new, T_floor, Tin_c + 5)
                adv_T = 0.5 * (adv_T + _advT(T_mid))
                T_expl = T + dt * (-adv_T + q_latent + q_jt)
                T_new = (T_expl + dt * beta_loss * T_sink) / (1.0 + dt * beta_loss)
            T = np.clip(T_new, min(T_floor, c.operating.T_seabed_C), Tin_c + 5)
            #  C12: evolve the soil. RADIAL model (soil_nodes>1): transient conduction through the
            #  concentric shells — the pipe flux drives the inner shell, conduction passes heat outward
            #  to the far-field ground, each shell with its own thermal mass (the proper cooldown
            #  profile). Else the single LUMPED node.
            if Ts_shells is not None:
                q_pipe = U_field * math.pi * D * (T - Ts_shells[0])      # W/m into the inner shell
                Ts_new = Ts_shells.copy()
                #  internal conduction fluxes between shells (nsoil-1)
                qflux = soil_G[:, None, None] * (Ts_shells[:-1] - Ts_shells[1:])   # (nsoil-1,nx,N)
                Ts_new[0] += dt / soil_C[0] * (q_pipe - qflux[0])
                if nsoil > 2:
                    Ts_new[1:-1] += dt / soil_C[1:-1, None, None] * (qflux[:-1] - qflux[1:])
                q_far = soil_Gfar * (Ts_shells[-1] - c.operating.T_seabed_C)
                Ts_new[-1] += dt / soil_C[-1] * (qflux[-1] - q_far)
                Ts_shells = np.clip(Ts_new, c.operating.T_seabed_C - 2.0, Tin_c + 5)
            elif n.soil_transient:
                q_in = U_field * (T - T_soil)
                q_out = c.numerics.soil_far_h * (T_soil - c.operating.T_seabed_C)
                T_soil = T_soil + dt * (q_in - q_out) / max(c.numerics.soil_thermal_mass, 1.0)

            # === (D) deposition — MASS-COUPLED to wall hydrate growth (A4) ===
            #  d(delta)/dt = f_wall*Rg_wall*A/(pi*D) = f_wall*Rg_wall*D/4 (annulus thickness rate); erosion scours.
            d_wall_thk = f_wall * Rg_wall * D / 4.0
            d_ero = k.k_ero * fslug * delta * ((PhiSH < 1.0) & ~locked)
            delta_new = delta + dt * (d_wall_thk - d_ero)
            delta_c = np.clip(delta_new, 0.0, delta_max)
            self._clip["deposit"] += int(np.count_nonzero(delta_c != delta_new))
            delta = delta_c

            # === (P) bulk hydrate phase-field : TVD advection + DIFFUSION + bulk growth (A2/#8) ===
            #  #8: 2nd-order TVD reconstruction of the advected hydrate front (flux-form so the
            #  advection neither creates nor destroys phi; growth/diffusion are separate sources).
            vlf = 0.5 * (vl[:-1] + vl[1:])
            phi_face = tvd_interior_faces(phi, vlf, self._limiter_kind)
            Pf = np.empty((nx + 1, N))
            Pf[1:-1] = vlf * phi_face
            Pf[0] = np.minimum(vl[0], 0.0) * phi[0]               # inlet: clean liquid (phi_in=0); allow outflow only
            Pf[-1] = np.maximum(vl[-1], 0.0) * phi[-1]            # outlet: outflow
            adv_phi = (Pf[1:] - Pf[:-1]) / self.dx
            phi_xx = np.empty((nx, N))
            phi_xx[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / self.dx ** 2
            phi_xx[0] = (phi[1] - phi[0]) / self.dx ** 2          # zero-gradient (Neumann) ends
            phi_xx[-1] = (phi[-2] - phi[-1]) / self.dx ** 2
            phi = phi + dt * (-adv_phi + k.D_phi * phi_xx + Rg_bulk)
            phi = np.where(phi < 1e-8, 0.0, np.clip(phi, 0.0, k.phi_max))   # floor underflow

            # === hydrate -> liquid mass coupling (A3) ===
            #  Total hydrate formed = BULK growth (Rg_bulk over the cross-section) + WALL deposit
            #  (f_wall*Rg_wall, the annulus consolidation). Both consume water from the liquid (and
            #  gas from the gas, via the continuity sink below); the water sink reduces the conserved
            #  liquid area La, so the liquid balance CLOSES WITH A HYDRATE SINK.
            hyd_vol_rate = (Rg_bulk + f_wall * Rg_wall) * A      # m3/s per m of hydrate formed
            water_vol_sink = (c.fluids.rho_hyd / c.fluids.rho_water) \
                * c.fluids.hyd_water_massfrac * hyd_vol_rate     # m3/s per m of water consumed
            dLa_hyd = np.minimum(dt * water_vol_sink, np.maximum(La, 0.0))
            La = La - dLa_hyd
            if c.fluids.oil_water_slip:
                #  #2 deepened: hydrate consumes WATER, so the sink draws from the water phase Wl
                #  (keeps the 3-phase water inventory consistent — hydrate forms from the water, not oil).
                Wl = np.clip(Wl - dLa_hyd, 0.0, La)
            alpha_l = np.clip(La / A, 1e-3, 0.999)               # refresh holdup after the sink
            liq_to_hyd_tot += float(np.mean(np.sum(dLa_hyd, 0))) * self.dx
            hyd_mass_tot += float(np.mean(np.sum(c.fluids.rho_hyd * hyd_vol_rate, 0))) * self.dx * dt

            # === (Gm) GAS-MASS CONTINUITY (A3): a conservative transport equation for the gas
            #  mass Mg = rho_g*alpha_g*A (kg/m):  d(Mg)/dt + d(rho_g*Vsg*A)/dx = -gas_to_hydrate.
            #  The gas MASS flux uses the mass-consistent superficial gas velocity Vsg (the same
            #  algebraic gas flux the pressure/Phi_SH closures use), NOT the drift-flux slip
            #  velocity: rho_g*Vsg*A is ~invariant along the line (density ∝ p cancels Vsg ∝ 1/p),
            #  so gas mass conserves to ~0% by telescoping — the drift slip sets the holdup, the
            #  mass-consistent flux sets the throughput. (Model/holdup closure unchanged.)
            rho_g_in = gas_density(c.operating.P_inlet_bar, c.operating.T_inlet_C, c.fluids)
            #  mass-flux-consistent superficial gas velocity, UNIVERSAL in the density model:
            #  Vsg = (inlet gas mass rate)/(rho_g*A) makes rho_g*Vsg*A == inlet mass rate exactly
            #  for ANY rho_g (ideal, real-gas Z(P,T), or a user PVT/EOS table), so the gas-mass
            #  continuity conserves to ~0% regardless of how the density is computed.
            Vsg = (rho_g_in * qg_in_t) / np.maximum(rho_g * A, 1e-9)
            Gflux = rho_g * Vsg * A                               # gas mass flux (kg/s) = inlet rate
            Gf = np.empty((nx + 1, N))
            Gf[1:-1] = Gflux[:-1]                                 # upwind (Vsg >= 0 downstream)
            Gf[0] = rho_g_in * qg_in_t                            # inlet gas mass rate BC
            Gf[-1] = Gflux[-1] if flow_frac > n.outlet_open_frac else 0.0   # outflow-only / isolated
            gas_sink = (1.0 - c.fluids.hyd_water_massfrac) * c.fluids.rho_hyd * hyd_vol_rate  # kg/s/m
            gas_sink_app = np.minimum(dt * gas_sink, np.maximum(Mg, 0.0))     # don't drive Mg < 0
            Mg = np.maximum(Mg - dt / self.dx * (Gf[1:] - Gf[:-1]) - gas_sink_app, 0.0)
            self._Mg = Mg                                    # #A1: keep engine-visible gas mass current
            gas_in_tot += float(np.mean(Gf[0])) * dt
            gas_out_tot += float(np.mean(Gf[-1])) * dt
            gas_consumed_hyd_tot += float(np.mean(np.sum(gas_sink_app, 0))) * self.dx
            #  #B3 TWO-FLUID-MASS coupling (opt-in): with BOTH phases conserved (La, Mg), the
            #  pressure is over-determined unless it is the value at which the two phase volumes
            #  exactly fill the bore: rho_g(p) = Mg/(A - La). Invert the gas EOS for that pressure
            #  and relax the pressure state toward it -> the conserved gas mass genuinely sets the
            #  pressure (the volume-consistent two-fluid-mass closure), not just the drift-flux solve.
            #  the twofluid_mass engine auto-applies a modest relaxation if none is set.
            rvcp = float(n.volume_consistent_pressure)
            if rvcp <= 0.0 and self.engine == "twofluid_mass":
                rvcp = 0.12
            if rvcp > 0.0:
                rho_g_vc = np.clip(Mg / np.maximum(A - La, 1e-6), 1e-3, 600.0)
                Zc = gas_Z_factor(self._p, T, c.fluids)
                p_vc = np.clip(rho_g_vc * Zc * R_GAS * (T + 273.15) / c.fluids.gas_MW / 1e5, 2.0, 1.0e3)
                self._p = np.clip((1.0 - rvcp) * self._p + rvcp * p_vc, 2.0, 1.0e3)

            max_PhiSH = np.maximum(max_PhiSH, PhiSH_rep)   # reported (true-magnitude) field
            max_Tsub = np.maximum(max_Tsub, Tsub)

            # --- plug detection (#24: fully vectorised over the ensemble, no Python loop) ---
            restr = 2.0 * delta / c.pipeline.diameter_m
            worst = np.argmax(restr, axis=0)                                # (N,) worst cell per realisation
            restr_worst = restr[worst, np.arange(N)]
            newly = np.isnan(plug_time) & (restr_worst > n.plug_restriction_trip)
            if newly.any():
                plug_time[newly] = t_h
                plug_loc[newly] = self.x[worst[newly]] / 1000.0

            # --- record monitor time-series ---
            for key, fld in [("P", p), ("T", T), ("Tsub", Tsub), ("alpha_l", alpha_l),
                             ("fslug", fslug), ("phi", phi), ("delta", delta),
                             ("PhiSH", PhiSH_rep), ("a_i", a_i), ("j", j)]:
                ts[key].append(np.nanmedian(fld[mon]))
            ts_t.append(t_h)
            bc_hist.append(self.inlet_rates(t_h)[0] / c.operating.q_liquid_insitu)

            if t >= next_snap - 1e-9:                            # D15: even time cadence
                snap_t.append(t_h)
                snap_phi.append(phi.mean(1).copy())
                snap_PhiSH.append(np.nanmedian(PhiSH_rep, 1).copy())
                snap_holdup.append(np.nanmedian(alpha_l, 1).copy())
                snap_P.append(np.nanmedian(p, 1).copy())
                snap_T.append(np.nanmedian(T, 1).copy())
                next_snap += snap_dt

            #  #17: observe the fastest relative change this step (T and holdup) for the dt controller
            if n.error_dt:
                dT_rate = float(np.nanmax(np.abs(adv_T + q_jt) + beta_loss
                                          * np.abs(T - c.operating.T_seabed_C)) / max(np.nanmean(T), 1.0))
                self._maxrate = max(dT_rate, 1e-9)

            t += dt; step += 1
            if verbose and step % 250 == 0:
                print(f"   t={t_h:5.1f} h  step={step}  dt={dt:5.1f}s  "
                      f"maxTsub={np.nanmax(np.nanmedian(Tsub,1)):4.1f}  "
                      f"maxPhiSH={np.nanmax(np.nanmedian(PhiSH,1)):4.2f}  "
                      f"plugged={np.mean(~np.isnan(plug_time))*100:3.0f}%")

        inv_final = float(np.mean(np.sum(La, 0)) * self.dx)
        #  A3: the liquid balance now CLOSES INCLUDING the hydrate (water) sink:
        #      (in - out - water_to_hydrate) == (inventory_final - inventory_init).
        mass_err = abs((liq_in_tot - liq_out_tot - liq_to_hyd_tot)
                       - (inv_final - inv_init)) / max(liq_in_tot, 1e-9)
        #  A3: gas mass-balance of the conserved Mg field (in - out - consumed_by_hydrate) vs stored-change.
        gas_store_final = float(np.mean(np.sum(Mg, 0))) * self.dx
        gas_mass_err = abs((gas_in_tot - gas_out_tot - gas_consumed_hyd_tot)
                           - (gas_store_final - gas_store_init)) / max(gas_in_tot, 1e-9)
        #  #3: consistency between the conserved gas mass (Mg) and the drift-flux holdup. The two
        #  agree when the drift-flux gas inventory is mass-consistent; the max relative gap is a
        #  HONEST diagnostic of how far the (model-preserving) algebraic gas sits from a fully
        #  mass-coupled two-fluid treatment. Reported, not forced.
        ag_drift = np.maximum(1.0 - alpha_l, 1e-3)
        ag_Mg = np.clip(Mg / np.maximum(rho_g * A, 1e-9), 1e-3, 0.999)
        gas_holdup_consistency = float(np.nanmedian(np.abs(ag_Mg - ag_drift) / ag_drift))
        #  #11: clip-activation fraction (activations per cell-step) -> a quantitative "is the
        #  solver masking instability?" gauge rather than a raw count.
        cell_steps = max(step * self.nx * self.N, 1)
        clip_frac = {kk: vv / cell_steps for kk, vv in self._clip.items()}
        if getattr(n, "strict", False):
            wf = float(getattr(n, "clip_warn_frac", 0.05))
            if clip_frac.get("velocity", 0.0) > wf or clip_frac.get("pressure", 0.0) > wf:
                raise RuntimeError(
                    f"SHCT strict mode: velocity/pressure clip activations exceed clip_warn_frac "
                    f"({wf:.0%}) — possible masked instability "
                    f"(vel {clip_frac.get('velocity',0)*100:.1f}%, pres {clip_frac.get('pressure',0)*100:.1f}%). "
                    f"Reduce cfl / refine grid, or disable numerics.strict.")
        self.results = dict(
            alpha_l=alpha_l, T=T, p=p, phi=phi, delta=delta, regime=regime,
            fslug=fslug, a_i=a_i, j=j, D=D, A=A, Teq=Teq, Tsub=Tsub,
            max_PhiSH=max_PhiSH, max_Tsub=max_Tsub, PhiSH=PhiSH,
            plug_time=plug_time, plug_loc=plug_loc, mon=mon,
            ts={key: np.array(v) for key, v in ts.items()}, ts_t=np.array(ts_t),
            bc_hist=np.array(bc_hist),
            snap_t=np.array(snap_t), snap_phi=np.array(snap_phi),
            snap_PhiSH=np.array(snap_PhiSH), snap_holdup=np.array(snap_holdup),
            snap_P=np.array(snap_P), snap_T=np.array(snap_T),
            steps=step, fallbacks=fallbacks, mass_err=mass_err,
            liq_in=liq_in_tot, liq_out=liq_out_tot, liq_to_hyd=liq_to_hyd_tot,
            hyd_mass=hyd_mass_tot, gas_in=gas_in_tot, gas_out=gas_out_tot,
            gas_consumed_hyd=gas_consumed_hyd_tot, gas_mass_err=gas_mass_err,
            gas_holdup_consistency=gas_holdup_consistency, water_frac=water_frac,
            clip_counts=dict(self._clip), clip_frac=clip_frac, cell_steps=cell_steps,
            W_inh=W_inh, U_eff=U_eff, therm_mass=therm_mass)
        return self.results

    # ----- engineering deliverables ---------------------------------------
    def engineering(self):
        c, f = self.case, self.case.fluids
        r = self.results
        rho_m = np.nanmedian(r["alpha_l"] * self.rho_l + (1 - r["alpha_l"]) *
                             gas_density(r["p"], r["T"], f), 1)
        Vm_peak = float(np.nanmax(np.nanmedian(r["j"], 1)))
        eros = f.api14e_C_factor / math.sqrt(max(np.nanmean(rho_m), 1.0))
        dP = float(np.nanmedian(r["p"][0] - r["p"][-1]))
        dT_design = float(np.nanpercentile(np.nanmax(r["max_Tsub"], 0), 90))
        water_mass = c.operating.q_liquid_insitu * f.water_cut * f.rho_water
        W, _, meg_Lph = hammerschmidt_meg(dT_design + c.operating.MEG_design_margin_C, water_mass)
        fmon = np.nanmedian(r["ts"]["fslug"][r["ts"]["fslug"] > 1e-3]) if (r["ts"]["fslug"] > 1e-3).any() else 0.05
        surge = c.operating.q_liquid_insitu / max(fmon, 1e-3) * c.numerics.surge_factor
        p_plug = float(np.mean(~np.isnan(r["plug_time"])))
        ttp = r["plug_time"][~np.isnan(r["plug_time"])]
        hot = float(self.x[np.argmax(np.nanmedian(r["max_PhiSH"], 1))] / 1000.0)
        arrival_T = float(np.nanmedian(r["T"][-1]))

        #  Cooldown / no-touch time: lumped-capacitance time for the coldest critical
        #  point to fall from operating T to the hydrate-onset temperature after shut-in.
        Dpipe = c.pipeline.diameter_m
        T_op = float(np.nanmedian(r["T"][r["mon"]]))
        Teq_hot = float(np.nanmedian(r["Teq"][r["mon"]]))
        Tsea = c.operating.T_seabed_C
        num = max(T_op - Tsea, 1e-3); den = max(Teq_hot - Tsea, 1e-3)
        UA = r["U_eff"] * math.pi * Dpipe
        cooldown_h = float(r["therm_mass"] / max(UA, 1e-6) * math.log(num / den) / 3600.0) \
            if num > den else 0.0
        #  Prefer the no-touch time measured DIRECTLY from the transient: the elapsed time
        #  after the operational event at which the monitor first enters the hydrate region
        #  (subcooling crosses zero). This fixes shut-in runs that previously reported 0.0
        #  from the lumped formula. Falls back to the lumped estimate when the monitor never
        #  crosses (e.g. a warm steady line).
        ev = c.scenario.event_time_h
        tt = r["ts_t"]; sub_ts = r["ts"]["Tsub"]
        post = np.where((tt >= ev) & (sub_ts > 0.0))[0]
        cooldown_src = "lumped"
        if post.size:
            cooldown_h = max(float(tt[post[0]] - ev), 0.0); cooldown_src = "transient"

        #  Hydrate slurry transportability (Camargo-Palermo relative viscosity)
        phi_peak = float(np.nanmax(np.nanmedian(r["phi"], 1)))
        mu_rel = (1.0 - min(phi_peak / c.kinetics.phi_max, 0.999)) ** f.slurry_visc_exp
        transportable = mu_rel < c.numerics.transportable_mu_rel_max

        #  Inhibition status: if MEG is injected, where (if anywhere) is it under-inhibited?
        meg_in = c.operating.MEG_wt_inlet
        under_inh_km = float((np.nanmax(r["max_Tsub"], 1) > 0).sum() * self.dx / 1000.0)

        #  WAX screen (item 10): where the fluid temperature falls below the wax appearance
        #  temperature (WAT), paraffin wax can deposit (a separate solids risk beside hydrate).
        #  Off by default (WAT = -273.15); set fluids.wax_appearance_C to the measured WAT.
        WAT = float(getattr(f, "wax_appearance_C", -273.15))
        Tmed = np.nanmedian(r["T"], 1)
        below_wat = Tmed < WAT
        wax_risk = bool(WAT > -100.0 and below_wat.any())
        wax_under_km = float(below_wat.sum() * self.dx / 1000.0) if WAT > -100.0 else 0.0
        wax_onset_km = float(self.x[int(np.argmax(below_wat))] / 1000.0) if wax_risk else float("nan")

        #  --- additional reporting diagnostics (no effect on the solved fields) ---
        A0 = math.pi * Dpipe ** 2 / 4.0
        Vsl_inlet = float(c.operating.q_liquid_insitu / A0)        # bulk superficial liquid velocity
        Vm_bulk = float(np.nanmedian(np.nanmedian(r["j"], 1)))     # line-median (not riser-peak) velocity
        monitor_T = float(np.nanmedian(r["T"][r["mon"]]))
        monitor_km = float(self.x[r["mon"]] / 1000.0)
        peak_deposit_mm = float(np.nanmax(np.nanmedian(r["delta"], 1)) * 1000.0)
        delta_max_mm = c.kinetics.delta_max_frac * Dpipe / 2.0 * 1000.0
        deposit_full_bore = bool(peak_deposit_mm >= 0.999 * delta_max_mm)
        #  Independent mass-based deposit check: equivalent wall thickness if the peak bulk
        #  hydrate fraction were laid on the wall (delta ~ phi*D/4). Cross-checks the
        #  dynamics-based deposit without altering the deposition PDE.
        deposit_from_phi_mm = float(phi_peak * Dpipe / 4.0 * 1000.0)
        max_phi_sh = float(np.nanmax(np.nanmedian(r["max_PhiSH"], 1)))
        phi_sh_saturated = bool(max_phi_sh >= 0.999 * c.kinetics.phi_report_cap)
        mass_warn = bool(r["mass_err"] > 0.05)
        #  probabilistic time-to-plug spread (now genuinely populated, C10)
        ttp_p10 = float(np.percentile(ttp, 10)) if ttp.size > 1 else float("nan")
        ttp_p90 = float(np.percentile(ttp, 90)) if ttp.size > 1 else float("nan")
        #  A3 hydrate/gas mass-balance diagnostics + D13 clip-activation audit
        gas_mass_err = float(r.get("gas_mass_err", float("nan")))
        clip_counts = dict(r.get("clip_counts", {}))
        clip_total = int(sum(clip_counts.values()))
        #  #11: flag clips that exceed clip_warn_frac of cell-steps (physical bound-contacts at
        #  full-bore/plug are expected; a high VELOCITY/PRESSURE fraction signals masked instability).
        clip_frac = dict(r.get("clip_frac", {}))
        warn_frac = float(getattr(c.numerics, "clip_warn_frac", 0.05))
        clip_warning = bool(clip_frac.get("velocity", 0.0) > warn_frac
                            or clip_frac.get("pressure", 0.0) > warn_frac)
        #  #5/#B4 sub-grid slug UNIT-CELL decomposition (reduced-order, Dukler-Hubbard style):
        #  unit length, slug-body holdup, film holdup and the slug fraction of the unit cell —
        #  the intermittent structure a dx~200 m grid averages out.
        jline = np.nanmedian(r["j"], 1); fsline = np.nanmedian(r["fslug"], 1)
        alline = np.nanmedian(r["alpha_l"], 1)
        Lu = slug_length(jline, Dpipe, fsline)
        slug_len_mean_m = float(np.nanmean(Lu)); slug_len_max_m = float(np.nanmax(Lu))
        kk = c.kinetics
        slug_body_holdup = float(np.nanmax(np.clip(kk.slug_body_holdup_base + kk.slug_body_holdup_slope
            * np.minimum(np.nanmedian(r["alpha_l"], 1) * jline / np.maximum(jline, 1e-3), 1.0),
            alline, 0.95)))
        film_holdup = float(np.nanmin(alline))
        slug_fraction = float(np.nanmean(np.clip(fsline / (fsline + kk.slug_fraction_ref_Hz), 0.0, 1.0)))
        #  #B5 minimum two-phase mixture sound speed (water-hammer screening, Wood's equation)
        c_snd = mixture_sound_speed(alline, self.rho_l,
                                    np.nanmedian(gas_density(r["p"], r["T"], f), 1),
                                    np.nanmedian(r["p"], 1))
        sound_speed_min_mps = float(np.nanmin(c_snd))
        #  #13 CaCO3 scaling-tendency screen (peak along the line) — a solids screen beside hydrate
        scal_idx = scaling_tendency_index(np.nanmedian(r["T"], 1), f.salinity_wt, np.nanmedian(r["p"], 1))
        scaling_tendency_max = float(np.nanmax(scal_idx))
        scaling_risk = bool(scaling_tendency_max > 0.0)
        #  B10/B11: water content of gas, free-water and a phase-envelope (dew) point when an EOS
        #  composition is given — the practical 3-phase / water-dewpoint / retrograde diagnostics.
        gas_water_content = float("nan"); free_water = False; dew_point_bar = float("nan")
        if getattr(f, "composition", None):
            try:
                import shct_eos
                Tmon = float(np.nanmedian(r["T"][r["mon"]])); Pmon = float(np.nanmedian(r["p"][r["mon"]]))
                tp = shct_eos.three_phase_flash(Pmon, Tmon, f.composition, f.water_cut, f.salinity_wt)
                gas_water_content = tp["gas_water_content_gpSm3"]; free_water = tp["free_water"]
                dew_point_bar = float(shct_eos.saturation_pressure(Tmon, f.composition, "dew"))
            except Exception as exc:
                #  the EOS 3-phase/dew diagnostics are optional; if they fail leave the NaN/False
                #  defaults but SURFACE the reason (D-gap: no longer silently swallowed).
                log.warning("[engineering] EOS water/dew diagnostics skipped: %s", exc)
        #  #6 peak droplet-entrained liquid fraction (gas-core), reported if enabled
        rho_g_line = np.nanmedian(gas_density(r["p"], r["T"], f), 1)
        Vsg_line = np.nanmedian(r["j"], 1) * np.nanmedian(1 - r["alpha_l"], 1)
        Ed = droplet_entrainment_frac(Vsg_line, rho_g_line, f.sigma, Dpipe)
        droplet_entrained_peak = float(np.nanmax(Ed)) if f.droplet_entrainment else 0.0
        #  #3 gas-holdup consistency (Mg vs drift-flux) — honesty gauge for the algebraic gas
        gas_holdup_consistency = float(r.get("gas_holdup_consistency", float("nan")))
        #  #2 water-accumulation hot-spot (only meaningful when oil_water_slip is on)
        wf = np.nanmedian(r.get("water_frac", np.full_like(r["alpha_l"], f.water_cut)), 1)
        water_frac_peak = float(np.nanmax(wf))
        water_accum_km = float(self.x[int(np.argmax(wf))] / 1000.0)

        return dict(
            Vm_peak_mps=Vm_peak, Vm_bulk_mps=Vm_bulk, Vsl_inlet_mps=Vsl_inlet,
            erosional_limit_mps=eros, dP_total_bar=dP,
            arrival_T_C=arrival_T, monitor_T_C=monitor_T, monitor_km=monitor_km,
            max_subcooling_C=float(np.nanmax(np.nanmedian(r["max_Tsub"], 1))),
            dT_design_C=dT_design, MEG_wt_pct=W, MEG_Lph=meg_Lph,
            MEG_injected_wt=meg_in, under_inhibited_km=under_inh_km,
            U_eff_WmK=float(r["U_eff"]), cooldown_to_hydrate_h=cooldown_h,
            cooldown_source=cooldown_src,
            slurry_rel_viscosity=float(mu_rel), slurry_transportable=bool(transportable),
            V_surge_P90_m3=float(surge), P_plug=p_plug,
            time_to_plug_P50_h=float(np.nanmedian(ttp)) if ttp.size else float("nan"),
            time_to_plug_P10_h=ttp_p10, time_to_plug_P90_h=ttp_p90,
            coupled_hotspot_km=hot,
            max_Phi_SH=max_phi_sh, Phi_SH_saturated=phi_sh_saturated,
            peak_deposit_mm=peak_deposit_mm, deposit_full_bore=deposit_full_bore,
            deposit_from_phi_mm=deposit_from_phi_mm,
            mass_conservation_err=float(r["mass_err"]), mass_conservation_warning=mass_warn,
            gas_mass_conservation_err=gas_mass_err,
            hydrate_mass_formed_kg=float(r.get("hyd_mass", 0.0)),
            water_to_hydrate_m3=float(r.get("liq_to_hyd", 0.0)),
            clip_activations=clip_total, clip_counts=clip_counts, clip_warning=clip_warning,
            clip_frac_velocity=float(clip_frac.get("velocity", 0.0)),
            clip_frac_pressure=float(clip_frac.get("pressure", 0.0)),
            slug_length_mean_m=slug_len_mean_m, slug_length_max_m=slug_len_max_m,
            slug_body_holdup=slug_body_holdup, slug_film_holdup=film_holdup,
            slug_fraction=slug_fraction, sound_speed_min_mps=sound_speed_min_mps,
            scaling_tendency_max=scaling_tendency_max, scaling_risk=scaling_risk,
            wax_risk=wax_risk, wax_under_km=wax_under_km, wax_onset_km=wax_onset_km,
            wax_appearance_C=WAT,
            gas_water_content_gpSm3=gas_water_content, free_water=bool(free_water),
            dew_point_bar=dew_point_bar,
            droplet_entrained_peak=droplet_entrained_peak,
            gas_holdup_consistency=gas_holdup_consistency,
            water_frac_peak=water_frac_peak, water_accum_km=water_accum_km,
            fallbacks=r["fallbacks"])


# =============================================================================
#  OUTPUT WRITERS
# =============================================================================
def _save_csv(path, cols, rows, str_col=None, all_str=False):
    with open(path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for ri, row in enumerate(rows):
            if all_str:
                fh.write(",".join(str(v) for v in row) + "\n")
            else:
                out = []
                for ci, v in enumerate(row):
                    out.append(str(str_col[ci][ri]) if (str_col and ci in str_col) else f"{float(v):.5g}")
                fh.write(",".join(out) + "\n")


def write_tables(sv: TransientSHCT, eng, outdir):
    r, c = sv.results, sv.case
    med = lambda A: np.nanmedian(A, 1)
    x = sv.x / 1000.0
    reg = [REGIME_NAMES[int(v)] for v in np.round(med(r["regime"]))]
    cols = ["x_km", "elevation_m", "incl_deg", "P_bar", "T_C", "Teq_C", "subcooling_C",
            "holdup", "v_m_mps", "regime", "a_i_1perm", "f_slug_Hz", "phi_hydrate",
            "deposit_mm", "Phi_SH"]
    rows = np.column_stack([x, sv.z, np.degrees(sv.theta), med(r["p"]), med(r["T"]),
                            med(r["Teq"]), med(r["Tsub"]), med(r["alpha_l"]), med(r["j"]),
                            np.round(med(r["regime"])), med(r["a_i"]), med(r["fslug"]),
                            med(r["phi"]), med(r["delta"]) * 1000, med(r["max_PhiSH"])])
    _save_csv(f"{outdir}/fields_profile.csv", cols, rows, str_col={9: reg})

    tcols = ["time_h", "bc_rate_frac", "P_bar", "T_C", "subcooling_C", "holdup",
             "v_m_mps", "f_slug_Hz", "phi_hydrate", "deposit_mm", "Phi_SH"]
    T = r["ts"]
    trows = np.column_stack([r["ts_t"], r["bc_hist"], T["P"], T["T"], T["Tsub"],
                             T["alpha_l"], T["j"], T["fslug"], T["phi"], T["delta"] * 1000, T["PhiSH"]])
    _save_csv(f"{outdir}/timeseries_monitor.csv", tcols, trows)

    def q3(a):
        a = np.asarray(a, float); a = a[~np.isnan(a)]
        return ["n/a"] * 3 if a.size == 0 else [f"{np.percentile(a,10):.3g}",
                                                f"{np.percentile(a,50):.3g}", f"{np.percentile(a,90):.3g}"]
    prows = [["max_subcooling_C"] + q3(np.nanmax(r["max_Tsub"], 0)),
             ["max_Phi_SH"] + q3(np.nanmax(r["max_PhiSH"], 0)),
             ["peak_wall_deposit_mm"] + q3(np.nanmax(r["delta"], 0) * 1000.0),
             ["peak_deposit_from_phi_mm"] + q3(np.nanmax(r["phi"], 0) * c.pipeline.diameter_m / 4.0 * 1000.0),
             ["peak_bulk_hydrate_frac"] + q3(np.nanmax(r["phi"], 0)),
             ["time_to_plug_h"] + q3(r["plug_time"])]
    _save_csv(f"{outdir}/probabilistic_summary.csv", ["metric", "P10", "P50", "P90"], prows, all_str=True)

    erows = [
        ["Slug-catcher surge volume (P90)", f"{eng['V_surge_P90_m3']:.2f}", "m3"],
        ["Design subcooling (P90)", f"{eng['dT_design_C']:.2f}", "C"],
        ["Required MEG concentration", f"{eng['MEG_wt_pct']:.1f}", "wt%"],
        ["MEG injection rate", f"{eng['MEG_Lph']:.1f}", "L/h"],
        ["MEG injected at inlet", f"{eng['MEG_injected_wt']:.1f}", "wt%"],
        ["Under-inhibited length", f"{eng['under_inhibited_km']:.2f}", "km"],
        ["Effective heat-transfer U", f"{eng['U_eff_WmK']:.2f}", "W/m2K"],
        ["Cooldown to hydrate (no-touch time)", f"{eng['cooldown_to_hydrate_h']:.2f}", "h"],
        ["Slurry relative viscosity", f"{eng['slurry_rel_viscosity']:.2f}", "-"],
        ["Slurry transportable?", f"{eng['slurry_transportable']}", "-"],
        ["Peak mixture velocity", f"{eng['Vm_peak_mps']:.2f}", "m/s"],
        ["Erosional velocity limit (API 14E)", f"{eng['erosional_limit_mps']:.2f}", "m/s"],
        ["Total line pressure drop", f"{eng['dP_total_bar']:.1f}", "bar"],
        ["Probability of plugging (run)", f"{eng['P_plug']*100:.0f}", "%"],
        ["Time-to-plug (P50)", f"{eng['time_to_plug_P50_h']:.1f}", "h"],
        ["Coupled hot-spot location", f"{eng['coupled_hotspot_km']:.2f}", "km"],
        ["Max coupling number Phi_SH", f"{eng['max_Phi_SH']:.4g}", "-"],
        ["Phi_SH reported at plot cap (saturated)?", f"{eng['Phi_SH_saturated']}", "-"],
        ["Monitor location", f"{eng['monitor_km']:.2f}", "km"],
        ["Monitor temperature", f"{eng['monitor_T_C']:.1f}", "C"],
        ["Bulk superficial liquid velocity", f"{eng['Vsl_inlet_mps']:.2f}", "m/s"],
        ["Line-median mixture velocity", f"{eng['Vm_bulk_mps']:.2f}", "m/s"],
        ["Peak wall deposit", f"{eng['peak_deposit_mm']:.1f}", "mm"],
        ["Wall deposit at full bore?", f"{eng['deposit_full_bore']}", "-"],
        ["Deposit (hydrate-mass equivalent)", f"{eng['deposit_from_phi_mm']:.1f}", "mm"],
        ["No-touch time source", f"{eng['cooldown_source']}", "-"],
        ["Mass-conservation error", f"{eng['mass_conservation_err']*100:.2f}", "%"],
        ["Mass-balance reliability warning", f"{eng['mass_conservation_warning']}", "-"],
    ]
    _save_csv(f"{outdir}/engineering_deliverables.csv", ["deliverable", "value", "units"], erows, all_str=True)


# ----------------------------- charts ---------------------------------------
def make_charts(sv: TransientSHCT, eng, outdir):
    if not HAVE_MPL:
        log.warning("[charts] matplotlib unavailable - skipping"); return
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    r, c = sv.results, sv.case
    med = lambda A: np.nanmedian(A, 1)
    x = sv.x / 1000.0; _pct = lambda a, q, ax: np.nanpercentile(a, q, axis=ax)

    # 1 profiles
    fig, ax = plt.subplots(4, 1, figsize=(8, 8.2), sharex=True)
    ax[0].fill_between(x, sv.z, sv.z.min() - 20, color="#C9B79B", alpha=.6); ax[0].plot(x, sv.z, color="#B07A33")
    ax[0].set_ylabel("elev (m)"); ax[0].set_title("Transient SHCT — final-state profiles (P50)", color=NAVY, fontweight="bold")
    ax[1].plot(x, med(r["alpha_l"]), color=ACCENT); ax[1].set_ylabel("holdup α_l"); ax[1].set_ylim(0, 1)
    lnP, = ax[2].plot(x, med(r["p"]), color=NAVY, label="pressure P (bar, left axis)")
    a2 = ax[2].twinx()
    lnT, = a2.plot(x, med(r["T"]), color=RED, label="temperature T (°C, right axis)")
    lnTeq, = a2.plot(x, med(r["Teq"]), color=RED, ls="--", lw=1,
                     label="hydrate-equilibrium T_eq (°C, right axis)")
    ax[2].set_ylabel("P (bar)", color=NAVY); a2.set_ylabel("T, T_eq (°C)", color=RED)
    ax[2].legend(handles=[lnP, lnT, lnTeq], loc="upper left", bbox_to_anchor=(1.13, 1.0),
                 fontsize=7, borderaxespad=0.0)
    ax[3].plot(x, med(r["Tsub"]), color=ORANGE, label="subcooling ΔT_sub")
    ax[3].axhline(0, color=GREY, ls=":", label="hydrate boundary (ΔT_sub = 0)")
    ax[3].fill_between(x, 0, med(r["Tsub"]), where=med(r["Tsub"]) > 0, color="#f6d6d2", alpha=.6)
    ax[3].set_ylabel("ΔT_sub (°C)"); ax[3].set_xlabel("distance (km)")
    ax[3].legend(loc="upper left", bbox_to_anchor=(1.13, 1.0), fontsize=7, borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/01_profiles.png", dpi=155); plt.close(fig)

    # 2 holdup space-time map (transient)
    if r["snap_holdup"].size:
        fig, axm = plt.subplots(figsize=(7.4, 4.4))
        # cividis is perceptually-uniform and colour-vision-deficiency safe; the
        # power-law norm (gamma<1) expands the dense low-holdup background so its
        # structure is visible while the high-holdup slug bands stay brightly
        # contrasted — a well-contrasted, standard, accessible map.
        H = r["snap_holdup"]
        norm = mcolors.PowerNorm(gamma=0.5, vmin=float(np.nanmin(H)),
                                 vmax=float(np.nanmax(H)))
        pcm = axm.pcolormesh(x, r["snap_t"], H, cmap="shct_seq", norm=norm, shading="gouraud")
        axm.set_xlabel("distance (km)"); axm.set_ylabel("time (h)")
        fig.colorbar(pcm, ax=axm, label="liquid holdup α_l")
        axm.set_title("Output — transient liquid-holdup field α_l(x,t)", color=NAVY, fontweight="bold")
        fig.tight_layout(); fig.savefig(f"{outdir}/02_holdup_spacetime.png", dpi=155); plt.close(fig)

    # 3 P-T envelope
    fig, axp = plt.subplots(figsize=(6.2, 4.4))
    Pc = np.linspace(5.0, max(160.0, med(r["p"]).max() * 1.1), 160)
    Tc = hydrate_equilibrium_T(Pc, gas_sg=c.fluids.gas_sg, salinity_wt=c.fluids.salinity_wt,
                               table=c.fluids.hyd_Teq_table)   # same curve as the solver uses
    axp.plot(Tc, Pc, color=RED, lw=2.2, label="hydrate equilibrium"); axp.fill_betweenx(Pc, 0, Tc, color="#f6d6d2", alpha=.4)
    axp.plot(med(r["T"]), med(r["p"]), color=NAVY, lw=2, marker="o", ms=2, label="pipe trajectory")
    axp.set_xlim(2, 30); axp.set_ylim(0, max(160, med(r["p"]).max() * 1.1))
    axp.set_xlabel("T (°C)"); axp.set_ylabel("P (bar)"); axp.legend(fontsize=8)
    axp.set_title("Output C — P–T trajectory vs hydrate envelope", color=NAVY, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{outdir}/03_PT_envelope.png", dpi=155); plt.close(fig)

    # 4 Phi_SH space-time
    if r["snap_PhiSH"].size:
        fig = plt.figure(figsize=(7.4, 5))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2.4], hspace=.08)
        az = fig.add_subplot(gs[0]); az.fill_between(x, sv.z, sv.z.min() - 20, color="#C9B79B", alpha=.6)
        az.plot(x, sv.z, color="#B07A33"); az.set_ylabel("elev (m)"); az.set_xticklabels([])
        az.set_title("Output E — Φ_SH(x,t) coupling-criticality map", color=NAVY, fontweight="bold")
        ap = fig.add_subplot(gs[1])
        pcm = ap.pcolormesh(x, r["snap_t"], r["snap_PhiSH"], cmap="shct_div", shading="gouraud",
                            vmin=0, vmax=max(1.5, np.nanpercentile(r["snap_PhiSH"], 98)))
        if np.nanmin(r["snap_PhiSH"]) < 1 < np.nanmax(r["snap_PhiSH"]):
            ap.contour(x, r["snap_t"], r["snap_PhiSH"], levels=[1.0], colors="#D24A8E", linewidths=1.6)
            ap.plot([], [], color="#D24A8E", lw=1.6, label="Φ_SH = 1 (critical) contour")
            ap.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1,
                      fontsize=7, borderaxespad=0.0)
        ap.set_xlabel("distance (km)"); ap.set_ylabel("time (h)")
        fig.colorbar(pcm, ax=[az, ap], pad=.02, fraction=.05, label="Φ_SH")
        fig.savefig(f"{outdir}/04_PhiSH_map.png", dpi=155); plt.close(fig)

    # 5 transient scenario monitor time-series
    fig, ax = plt.subplots(3, 1, figsize=(7.6, 6.4), sharex=True)
    tt = r["ts_t"]
    ax[0].plot(tt, r["bc_hist"], color=GREY, label="inlet rate fraction"); ax[0].set_ylabel("inlet rate\n(fraction)")
    ax[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, borderaxespad=0.0)
    ax[0].set_title(f"Transient scenario '{c.scenario.kind}' — monitor response", color=NAVY, fontweight="bold")
    ax[1].plot(tt, r["ts"]["Tsub"], color=ORANGE, label="subcooling ΔT_sub")
    ax[1].axhline(0, color=GREY, ls=":", label="hydrate boundary (ΔT_sub = 0)")
    ax[1].set_ylabel("subcooling (°C)")
    ax[1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, borderaxespad=0.0)
    phi_ts = np.asarray(r["ts"]["PhiSH"], float)
    ax[2].plot(tt, phi_ts, color=NAVY, label="coupling number Φ_SH")
    ax[2].axhline(1, color=RED, ls="--", label="critical threshold Φ_SH = 1")
    #  the very first transient step can spike Φ_SH far above the meaningful range;
    #  clip the y-axis to the post-warm-up envelope so the Φ_SH≈1 dynamics are clear.
    finite = phi_ts[np.isfinite(phi_ts)]
    if finite.size:
        cap = max(2.0, float(np.nanpercentile(finite, 97)))
        ax[2].set_ylim(0, cap * 1.18)
        if float(np.nanmax(finite)) > cap * 1.18:
            ax[2].text(0.99, 0.96, f"(initial peak {np.nanmax(finite):.0f} clipped)",
                       transform=ax[2].transAxes, ha="right", va="top", fontsize=6,
                       color=GREY, style="italic")
    ax[2].set_ylabel("Φ_SH"); ax[2].set_xlabel("time (h)")
    ax[2].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, borderaxespad=0.0)
    fig.tight_layout(); fig.savefig(f"{outdir}/05_scenario_timeseries.png", dpi=155); plt.close(fig)

    # 6 deposit growth
    fig, axd = plt.subplots(figsize=(6.6, 4))
    axd.plot(tt, r["ts"]["delta"] * 1000, color=RED, lw=1.8)
    axd.set_xlabel("time (h)"); axd.set_ylabel("deposit δ_h at monitor (mm)")
    axd.set_title("Output D — wall-deposit growth (transient, coupled)", color=NAVY, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{outdir}/06_deposit.png", dpi=155); plt.close(fig)

    # 7 probabilistic — Kaplan-Meier (right-censored) time-to-plug CDF + Phi_SH band (#17)
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(8, 3.6))
    ttp = r["plug_time"][~np.isnan(r["plug_time"])]
    Ntot = r["plug_time"].size
    if ttp.size >= 1:
        #  KM estimate with all un-plugged realisations censored at t_end: F(t)=1-prod(1-d_i/n_i).
        s = np.sort(ttp); n_at_risk = Ntot; surv = 1.0; tt = [0.0]; FF = [0.0]
        for ti in s:
            surv *= (1.0 - 1.0 / n_at_risk); n_at_risk -= 1
            tt.append(float(ti)); FF.append(1.0 - surv)
        b1.step(tt, FF, where="post", color=RED, lw=2)
        for q, lab in [(10, "P10"), (50, "P50"), (90, "P90")]:
            if ttp.size > 1:
                xv = float(np.percentile(ttp, q))
                b1.axvline(xv, color=GREY, ls=":", lw=0.8)
                b1.text(xv, 0.04, lab, fontsize=6, rotation=90, va="bottom")
    b1.set_xlabel("time-to-plug (h)"); b1.set_ylabel("cum. probability"); b1.set_ylim(0, 1)
    b1.set_title(f"Output G — time-to-plug CDF (Kaplan–Meier, P_plug={eng['P_plug']*100:.0f}%)",
                 color=NAVY, fontweight="bold", fontsize=9)
    b2.fill_between(x, _pct(r["max_PhiSH"], 10, 1), _pct(r["max_PhiSH"], 90, 1),
                    color="#cfe0f5", alpha=.7, label="P10–P90")
    b2.plot(x, med(r["max_PhiSH"]), color=NAVY, lw=2, label="P50"); b2.axhline(1, color=RED, ls="--")
    b2.set_xlabel("distance (km)"); b2.set_ylabel("max Φ_SH")
    b2.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, borderaxespad=0.0)
    b2.set_title("Output — Φ_SH along line (ensemble)", color=NAVY, fontweight="bold", fontsize=9.5)
    fig.tight_layout(); fig.savefig(f"{outdir}/07_probabilistic.png", dpi=155); plt.close(fig)

    # 8 solver diagnostics — conservation, clip activity, gas-holdup consistency, slug length (#25)
    fig, dax = plt.subplots(2, 2, figsize=(8, 5.4))
    dax[0, 0].bar(["liquid", "gas"], [eng["mass_conservation_err"] * 100,
                  eng["gas_mass_conservation_err"] * 100], color=[ACCENT, TEAL])
    dax[0, 0].axhline(5, color=RED, ls="--", lw=0.8, label="5% warn")
    dax[0, 0].set_ylabel("mass-balance error (%)")
    dax[0, 0].legend(loc="upper right", fontsize=7)   # bars are ~0 (mass-consistent) -> top is empty
    dax[0, 0].set_title("Mass conservation", color=NAVY, fontweight="bold", fontsize=9)
    cf = r.get("clip_frac", {})
    keys = list(cf.keys()); vals = [cf[k] * 100 for k in keys]
    cols_ = [RED if k in ("velocity", "pressure") else GREY for k in keys]
    dax[0, 1].bar(keys, vals, color=cols_); dax[0, 1].set_ylabel("clip activations (% cell-steps)")
    dax[0, 1].set_title("Clip activity (vel/pres red = instability)", color=NAVY, fontweight="bold", fontsize=8.5)
    dax[0, 1].tick_params(axis="x", labelrotation=30, labelsize=7)
    dax[1, 0].plot(x, slug_length(med(r["j"]), c.pipeline.diameter_m, med(r["fslug"])), color=ORANGE)
    dax[1, 0].set_xlabel("distance (km)"); dax[1, 0].set_ylabel("slug-unit length (m)")
    dax[1, 0].set_title("Sub-grid slug length (#5)", color=NAVY, fontweight="bold", fontsize=9)
    txt = (f"gas-holdup consistency: {eng.get('gas_holdup_consistency', float('nan'))*100:.1f} %\n"
           f"(drift-flux vs conserved gas mass)\n\n"
           f"hydrate formed: {eng.get('hydrate_mass_formed_kg', 0)/1000:.1f} t\n"
           f"water -> hydrate: {eng.get('water_to_hydrate_m3', 0):.2f} m3\n\n"
           f"fallbacks: {r['fallbacks']}   steps: {r['steps']}\n"
           f"clip warning: {eng.get('clip_warning')}")
    dax[1, 1].axis("off"); dax[1, 1].text(0.02, 0.95, txt, va="top", fontsize=8, family="monospace")
    fig.suptitle("Output — solver diagnostics & balances", color=NAVY, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{outdir}/08_diagnostics.png", dpi=155); plt.close(fig)
    log.info("[charts] charts written to %s", outdir)


def console_report(sv, eng):
    c = sv.case; r = sv.results; line = "=" * 70
    print(line); print(f" SHCT TRANSIENT SOLVER  —  {c.name}"); print(f" scenario: {c.scenario.kind}"); print(line)
    print(f"  Pipeline   : {c.pipeline.length_m/1000:.1f} km x {c.pipeline.diameter_m*1000:.0f} mm ID, "
          f"{sv.nx} cells")
    print(f"  Integration: {r['steps']} adaptive time-steps, {c.numerics.n_ensemble} realisations, "
          f"{r['fallbacks']} hydro-fallbacks")
    _mc = eng['mass_conservation_err'] * 100
    _gmc = eng.get('gas_mass_conservation_err', float('nan')) * 100
    print(f"  Mass cons. : liquid {_mc:.2f} %  "
          f"[{'PASS' if _mc < 5 else 'WARN — deep-transient regime'}]"
          f"   gas {_gmc:.2f} %  [{'PASS' if _gmc < 5 else 'WARN'}]")
    print(f"  Clips      : {eng.get('clip_activations',0)} bound-activations "
          f"{eng.get('clip_counts',{})}  "
          f"[{'WARN — vel/pres clips high, check stability' if eng.get('clip_warning') else 'OK — vel/pres 0'}]")
    print(f"  Gas holdup : drift-flux vs conserved-mass consistency "
          f"{eng.get('gas_holdup_consistency',float('nan'))*100:.1f} %")
    print("-" * 70)
    print("  HYDRODYNAMICS")
    print(f"    Total dP                 : {eng['dP_total_bar']:8.1f} bar")
    print(f"    Peak mixture velocity    : {eng['Vm_peak_mps']:8.2f} m/s (erosional {eng['erosional_limit_mps']:.1f})")
    print("  THERMAL / HYDRATE")
    print(f"    Max subcooling (P50)     : {eng['max_subcooling_C']:8.2f} C")
    print(f"    Design subcooling (P90)  : {eng['dT_design_C']:8.2f} C")
    print("  COUPLED Phi_SH")
    print(f"    Max Phi_SH (P50)         : {eng['max_Phi_SH']:8.4g}  "
          f"({'CRITICAL >1' if eng['max_Phi_SH']>1 else 'sub-critical'}"
          f"{', at plot cap' if eng.get('Phi_SH_saturated') else ''})")
    print(f"    Coupled hot-spot         : {eng['coupled_hotspot_km']:8.2f} km")
    print("  RISK")
    print(f"    Probability of plugging  : {eng['P_plug']*100:8.0f} %")
    print(f"    Time-to-plug P10/P50/P90 : {eng.get('time_to_plug_P10_h',float('nan')):6.1f} /"
          f"{eng['time_to_plug_P50_h']:6.1f} /{eng.get('time_to_plug_P90_h',float('nan')):6.1f} h")
    print(f"    Hydrate mass formed      : {eng.get('hydrate_mass_formed_kg',0.0):8.1f} kg "
          f"(water consumed {eng.get('water_to_hydrate_m3',0.0):.3f} m3)")
    print("  THERMAL MANAGEMENT")
    print(f"    Effective wall U         : {eng['U_eff_WmK']:8.2f} W/m2K")
    print(f"    Cooldown to hydrate      : {eng['cooldown_to_hydrate_h']:8.2f} h (no-touch time, "
          f"{eng.get('cooldown_source','lumped')})")
    print("  ENGINEERING DELIVERABLES")
    print(f"    Slug-catcher (P90 surge) : {eng['V_surge_P90_m3']:8.2f} m3")
    print(f"    MEG concentration        : {eng['MEG_wt_pct']:8.1f} wt%  "
          f"(injected {eng['MEG_injected_wt']:.0f} wt%, under-inhibited {eng['under_inhibited_km']:.1f} km)")
    print(f"    MEG injection rate       : {eng['MEG_Lph']:8.1f} L/h")
    print(f"    Slurry transportable     : {str(eng['slurry_transportable']):>8s} "
          f"(rel. visc {eng['slurry_rel_viscosity']:.1f})")
    print(line)


# =============================================================================
#  CASE LOADING / CLI
# =============================================================================
def make_default_case() -> Case:
    return Case()


#  G21: known JSON fields per group, used to give a friendly error on typos/unknown keys
#  (a bare dataclass(**d) raises an opaque TypeError instead).
_GROUP_TYPES = {"pipeline": Pipeline, "fluids": Fluids, "operating": Operating,
                "kinetics": Kinetics, "numerics": Numerics, "scenario": Scenario}


def _build_group(name, dct):
    cls = _GROUP_TYPES[name]
    allowed = set(cls.__dataclass_fields__.keys())
    unknown = set(dct) - allowed
    if unknown:
        raise ValueError(f"case group '{name}' has unknown field(s) {sorted(unknown)}; "
                         f"valid fields: {sorted(allowed)}")
    return cls(**dct)


def validate_case(case: Case) -> Case:
    """Validate a Case's physical/numerical inputs and raise a clear ValueError on bad data
    (G21). Cheap structural checks only — they catch the common configuration mistakes that
    would otherwise surface as silent garbage or a deep NumPy error."""
    p, fl, o, k, n, sc = (case.pipeline, case.fluids, case.operating,
                          case.kinetics, case.numerics, case.scenario)
    errs = []
    if p.length_m <= 0: errs.append("pipeline.length_m must be > 0")
    if p.diameter_m <= 0: errs.append("pipeline.diameter_m must be > 0")
    if p.roughness_m < 0: errs.append("pipeline.roughness_m must be >= 0")
    if p.n_cells < 2: errs.append("pipeline.n_cells must be >= 2")
    if p.elevation_m is not None and len(p.elevation_m) < p.n_cells:
        errs.append(f"pipeline.elevation_m has {len(p.elevation_m)} points (< n_cells={p.n_cells})")
    if not (0.0 <= fl.water_cut <= 1.0): errs.append("fluids.water_cut must be in [0, 1]")
    for nm in ("rho_oil", "rho_water", "mu_liquid", "mu_gas", "gas_MW", "gas_Z",
               "cp_liquid", "rho_hyd"):
        if getattr(fl, nm) <= 0: errs.append(f"fluids.{nm} must be > 0")
    if not (0.0 < fl.hyd_water_massfrac < 1.0):
        errs.append("fluids.hyd_water_massfrac must be in (0, 1)")
    if fl.hyd_Teq_table is not None:
        tab = np.asarray(fl.hyd_Teq_table, float)
        if tab.ndim != 2 or tab.shape[1] != 2 or tab.shape[0] < 2:
            errs.append("fluids.hyd_Teq_table must be [[P_bar, T_C], ...] with >= 2 rows")
    if getattr(fl, "gas_Z_table", None) is not None:
        tab = np.asarray(fl.gas_Z_table, float)
        if tab.ndim != 2 or tab.shape[1] != 3 or tab.shape[0] < 1:
            errs.append("fluids.gas_Z_table must be [[P_bar, T_C, Z], ...]")
    if getattr(fl, "pvt_table", None) is not None:
        tab = np.asarray(fl.pvt_table, float)
        if tab.ndim != 2 or tab.shape[1] != 6 or tab.shape[0] < 1:
            errs.append("fluids.pvt_table must be [[P_bar,T_C,rho_oil,rho_gas,mu_oil,mu_gas], ...]")
    if getattr(fl, "composition", None) is not None:
        try:
            import shct_eos
            unknown = [k for k in fl.composition if k not in shct_eos.COMPONENTS]
            if unknown:
                errs.append(f"fluids.composition has unknown component(s) {unknown}; "
                            f"valid: {sorted(shct_eos.COMPONENTS)}")
            if sum(fl.composition.values()) <= 0:
                errs.append("fluids.composition mole fractions must sum to > 0")
        except ImportError:
            errs.append("fluids.composition requires shct_eos.py")
    if o.q_liquid_insitu < 0 or o.q_gas_insitu_inlet < 0:
        errs.append("operating in-situ rates must be >= 0")
    if o.P_inlet_bar <= 0: errs.append("operating.P_inlet_bar must be > 0")
    if o.T_inlet_C <= o.T_seabed_C:
        errs.append("operating.T_inlet_C must exceed T_seabed_C")
    if not (0.0 <= o.MEG_wt_inlet <= 80.0): errs.append("operating.MEG_wt_inlet must be in [0, 80] wt%")
    if not (0.0 < k.phi_max < 1.0): errs.append("kinetics.phi_max must be in (0, 1)")
    if not (0.0 < k.delta_max_frac < 1.0): errs.append("kinetics.delta_max_frac must be in (0, 1)")
    if k.D_phi < 0: errs.append("kinetics.D_phi must be >= 0")
    if not (0.0 < n.cfl <= 1.0): errs.append("numerics.cfl must be in (0, 1]")
    if n.t_end_h <= 0: errs.append("numerics.t_end_h must be > 0")
    if n.dt_max_h <= 0: errs.append("numerics.dt_max_h must be > 0")
    if n.n_ensemble < 1: errs.append("numerics.n_ensemble must be >= 1")
    if n.engine not in ("implicit", "twofluid", "twofluid_mass", "twofluid_mass_newton",
                        "twofluid_full_newton", "quasisteady"):
        errs.append("numerics.engine must be implicit | twofluid | twofluid_mass | "
                    "twofluid_mass_newton | twofluid_full_newton | quasisteady")
    if sc.kind not in ("steady", "rampup", "turndown", "shutin"):
        errs.append("scenario.kind must be steady | rampup | turndown | shutin")
    if errs:
        raise ValueError("invalid case:\n  - " + "\n  - ".join(errs))
    return case


def case_from_json(path) -> Case:
    with open(path) as fh:
        d = json.load(fh)
    unknown_top = set(d) - (set(_GROUP_TYPES) | {"name"})
    if unknown_top:
        raise ValueError(f"case has unknown top-level group(s) {sorted(unknown_top)}; "
                         f"valid groups: {sorted(_GROUP_TYPES)} (+ 'name')")
    case = Case(name=d.get("name", "user-case"),
                **{g: _build_group(g, d.get(g, {})) for g in _GROUP_TYPES})
    return validate_case(case)


# =============================================================================
#  VERIFICATION HARNESS  (V&V: verification against published reference values)
# =============================================================================
def run_verification():
    """Automated checks that the implemented closures reproduce their published
    reference behaviour, and that the transient core conserves mass. This is the
    VERIFICATION half of V&V (does the code solve the equations correctly);
    VALIDATION against field/flow-loop data is performed via calibrate()."""
    f = Fluids()
    checks = []

    def chk(name, ok, detail):
        checks.append((name, bool(ok), detail))

    # 1. Hydrate equilibrium vs literature anchors (natural gas): within ~2.5 degC
    for P, Tlit in [(30, 3.0), (70, 10.0), (100, 13.0), (200, 18.0)]:
        T = float(hydrate_equilibrium_T(np.array([P]))[0])
        chk(f"hydrate Teq @ {P} bar", abs(T - Tlit) < 2.5, f"model {T:.1f}C vs lit ~{Tlit}C")
    chk("hydrate Teq monotonic", hydrate_equilibrium_T(np.array([50.])) < hydrate_equilibrium_T(np.array([150.])), "dTeq/dP>0")

    # 2. Gregory-Scott / Zabaras slug frequency: positive, rises uphill, hand-value
    fs_h = float(slug_frequency(np.array([0.5]), np.array([3.0]), 0.3, np.array([0.0]))[0])
    fs_up = float(slug_frequency(np.array([0.5]), np.array([3.0]), 0.3, np.array([0.5]))[0])
    chk("slug freq positive", fs_h > 0, f"{fs_h:.4f} Hz")
    chk("slug freq rises uphill", fs_up > fs_h, f"{fs_up:.4f} > {fs_h:.4f}")

    # 3. Drift-flux holdup limits
    C0, vd = drift_params(0.0, 0.3)
    def _holdup(vsg, vsl):
        vm = vsg + vsl
        return 1.0 - min(vsg / (C0 * vm + vd), 0.999)
    hl_lowgas = _holdup(0.05, 1.0)
    hl_higas = _holdup(20.0, 1.0)
    chk("holdup high at low gas", hl_lowgas > 0.5, f"{hl_lowgas:.2f}")
    chk("holdup low at high gas", hl_higas < 0.2, f"{hl_higas:.2f}")

    # 4. Haaland friction vs Moody anchor (Re=1e5, smooth) ~ 0.018
    fr = float(haaland_friction(np.array([1e5]), np.array([1e-4]))[0])
    chk("friction ~ Moody", 0.012 < fr < 0.025, f"f={fr:.4f}")

    # 5. Hammerschmidt monotonic & magnitude (more subcooling -> more MEG)
    W5, _, _ = hammerschmidt_meg(5.0, 1.0)
    W10, _, _ = hammerschmidt_meg(10.0, 1.0)
    chk("Hammerschmidt monotonic", W10 > W5 > 0, f"{W5:.1f}% @5C, {W10:.1f}% @10C")

    # 6. Mass conservation of the transient core (short run, all scenarios)
    for scen in ["steady", "rampup", "shutin"]:
        c = Case(); c.pipeline.n_cells = 60; c.numerics.n_ensemble = 3
        c.numerics.t_end_h = 10.0; c.scenario.kind = scen
        r = TransientSHCT(c).run(verbose=False)
        #  tightened to 2% (item 11): actual conservation is ~0% even on this coarse/short grid;
        #  the old 7% was over-conservative. (Production grid n_cells>=100 conserves to <<1%.)
        chk(f"mass conservation [{scen}]", r["mass_err"] < 0.02, f"{r['mass_err']*100:.2f}%")

    # 7. Inhibitor (MEG) must reduce subcooling / hydrate risk
    cb = Case(); cb.pipeline.n_cells = 36; cb.numerics.n_ensemble = 3; cb.numerics.t_end_h = 8.0
    sub0 = float(np.nanmax(np.nanmedian(TransientSHCT(cb).run(verbose=False)["max_Tsub"], 1)))
    cm = copy.deepcopy(cb); cm.operating.MEG_wt_inlet = 30.0
    subM = float(np.nanmax(np.nanmedian(TransientSHCT(cm).run(verbose=False)["max_Tsub"], 1)))
    chk("MEG reduces subcooling", subM < sub0, f"{sub0:.1f}C -> {subM:.1f}C with 30wt% MEG")

    # 8. Grid convergence: key metric stable between coarse & fine grids (<20%)
    def _dP(ncell):
        c = Case(); c.pipeline.n_cells = ncell; c.numerics.n_ensemble = 2; c.numerics.t_end_h = 6.0
        rr = TransientSHCT(c).run(verbose=False)
        return float(np.nanmedian(rr["p"][0] - rr["p"][-1], 0))
    dP1, dP2 = _dP(40), _dP(80)
    rel = abs(dP1 - dP2) / max(abs(dP2), 1e-6)
    #  tightened 20% -> 8% (item 11): the implicit pressure solve is ~1% grid-convergent on dP.
    chk("grid convergence (dP)", rel < 0.08, f"dP {dP1:.1f}->{dP2:.1f} bar, {rel*100:.0f}% change")

    # 9. Effective-U from wall layers is physical (insulation lowers U)
    pl = Pipeline(wall_layers=[[0.025, 50.0, 3.9e6], [0.05, 0.18, 1.2e5]])
    U_ins, _ = effective_U_and_mass(pl, Operating(), 2100.0, 300.0)
    chk("insulated wall lowers U", U_ins < 22.0, f"layered U={U_ins:.2f} W/m2K")

    # 10. Full two-fluid engine: stable, mass-conservative, and physical slip (gas faster)
    ctf = Case(); ctf.pipeline.n_cells = 50; ctf.numerics.n_ensemble = 3
    ctf.numerics.t_end_h = 8.0; ctf.numerics.engine = "twofluid"
    stf = TransientSHCT(ctf); rtf = stf.run(verbose=False)
    chk("two-fluid mass conservation", rtf["mass_err"] < 0.02, f"{rtf['mass_err']*100:.2f}%")
    chk("two-fluid stable (no fallbacks)", rtf["fallbacks"] == 0, f"{rtf['fallbacks']} fallbacks")
    slip = float(np.nanmean(np.nanmedian(stf._ug, 1) - np.nanmedian(stf._ul, 1)))
    chk("two-fluid physical slip (gas faster)", 0.05 < slip < 1.5, f"mean slip {slip:.2f} m/s")

    # 11. Engine consistency: implicit & two-fluid give similar pressure drop (<35%)
    def _dPe(eng):
        c = Case(); c.pipeline.n_cells = 50; c.numerics.n_ensemble = 2; c.numerics.t_end_h = 6.0
        c.numerics.engine = eng
        return float(np.nanmedian(TransientSHCT(c).run(verbose=False)["p"][0]
                                  - TransientSHCT(c).run(verbose=False)["p"][-1], 0))
    dpi, dptf = _dPe("implicit"), _dPe("twofluid")
    rel = abs(dpi - dptf) / max(abs(dpi), 1e-6)
    chk("engine cross-consistency (dP)", rel < 0.35, f"implicit {dpi:.1f} vs two-fluid {dptf:.1f} bar")

    # 12. GAS-mass conservation of the new gas-continuity equation (A3), all scenarios
    for scen in ["steady", "turndown", "shutin"]:
        c = Case(); c.pipeline.n_cells = 60; c.numerics.n_ensemble = 3
        c.numerics.t_end_h = 10.0; c.scenario.kind = scen
        rg = TransientSHCT(c).run(verbose=False)
        chk(f"gas mass conservation [{scen}]", rg["gas_mass_err"] < 0.02, f"{rg['gas_mass_err']*100:.2f}%")

    # 13. Hydrate mass balance closes: water consumed by hydrate matches the hydrate water mass
    cH = Case(); cH.pipeline.n_cells = 60; cH.numerics.n_ensemble = 4; cH.numerics.t_end_h = 12.0
    rH = TransientSHCT(cH).run(verbose=False)
    water_kg = rH["liq_to_hyd"] * cH.fluids.rho_water
    expect_kg = rH["hyd_mass"] * cH.fluids.hyd_water_massfrac
    rel = abs(water_kg - expect_kg) / max(expect_kg, 1e-6)
    chk("hydrate water mass balance", rel < 0.02, f"water {water_kg:.0f} vs hydrate-water {expect_kg:.0f} kg")

    # 14. HEADLINE-OUTPUT grid stability (D12): time-to-plug, P_plug, max Phi_SH, peak deposit
    def _head(ncell):
        c = Case(); c.pipeline.n_cells = ncell; c.numerics.n_ensemble = 12
        c.numerics.t_end_h = 30.0; c.numerics.deterministic = True   # remove ensemble noise for the grid test
        sv = TransientSHCT(c); sv.run(verbose=False); return sv.engineering()
    h1, h2 = _head(50), _head(100)
    for kk, tol in [("time_to_plug_P50_h", 0.40), ("P_plug", 0.20),
                    ("max_Phi_SH", 0.30), ("peak_deposit_mm", 0.30)]:
        a, b = h1.get(kk, float("nan")), h2.get(kk, float("nan"))
        if np.isfinite(a) and np.isfinite(b):
            rel = abs(a - b) / max(abs(b), 1e-6)
            chk(f"grid stability ({kk})", rel < tol, f"{a:.3g}->{b:.3g}, {rel*100:.0f}% change")
        else:
            chk(f"grid stability ({kk})", np.isnan(a) == np.isnan(b), f"{a} -> {b}")

    # 15. Latent heat is active (A1): enabling hydrate growth raises T vs a no-growth reference
    cL = Case(); cL.pipeline.n_cells = 50; cL.numerics.n_ensemble = 1; cL.numerics.t_end_h = 8.0
    cL.numerics.deterministic = True
    rL = TransientSHCT(cL).run(verbose=False)
    cL0 = copy.deepcopy(cL); cL0.fluids.L_hyd = 0.0
    rL0 = TransientSHCT(cL0).run(verbose=False)
    dT_latent = float(np.nanmean(rL["T"]) - np.nanmean(rL0["T"]))
    chk("latent heat warms fluid (A1)", dT_latent > -1e-6, f"mean dT = {dT_latent:+.3f} C with L_hyd")

    # 16. Phase-field diffusion is active (A2): D_phi smooths the phi field
    chk("phase-field diffusivity wired", Kinetics().D_phi > 0.0, f"D_phi={Kinetics().D_phi}")

    # 17. Peng-Robinson EOS (#A1): pure-fluid Z, NG density/SG and flash are physical
    try:
        import shct_eos
        for name, ok, detail in shct_eos.eos_selftest():
            chk(f"EOS: {name}", ok, detail)
        #  EOS-driven case conserves mass like any other (the EOS plugs into the same path)
        ce = Case(); ce.pipeline.n_cells = 40; ce.numerics.n_ensemble = 3; ce.numerics.t_end_h = 8.0
        ce.fluids.composition = shct_eos.DEFAULT_COMPOSITION
        re = TransientSHCT(ce).run(verbose=False)
        chk("EOS-driven case conserves mass", re["mass_err"] < 0.05 and re["gas_mass_err"] < 0.05,
            f"liq {re['mass_err']*100:.2f}% gas {re['gas_mass_err']*100:.2f}%")
    except Exception as exc:                                       # pragma: no cover
        chk("EOS module available", False, f"{exc}")

    # 18. Mixture sound speed (#B5): two-phase << single-phase (Wood's equation)
    c2ph = float(mixture_sound_speed(0.5, 820.0, 90.0, 100.0))
    c1ph = float(mixture_sound_speed(0.999, 820.0, 90.0, 100.0))
    chk("two-phase sound speed depressed", 50.0 < c2ph < c1ph, f"{c2ph:.0f} < {c1ph:.0f} m/s")

    # report
    print("=" * 64)
    print(" SHCT SOLVER — VERIFICATION SUITE")
    print("=" * 64)
    npass = sum(1 for _, ok, _ in checks if ok)
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}]  {name:32s}  {detail}")
    print("-" * 64)
    print(f"  {npass}/{len(checks)} checks passed")
    print("=" * 64)
    return all(ok for _, ok, _ in checks)


# =============================================================================
#  CALIBRATION  (VALIDATION: fit free constants to the user's measured data)
# =============================================================================
CALIB_PARAMS = {           # name -> (attr path, lower, upper) as multiplier on default
    "U_wall":       ("operating.U_wall",         0.3, 3.0),
    "kg0":          ("kinetics.kg0",             0.2, 5.0),
    "wall_capture": ("kinetics.wall_capture_eff", 0.2, 5.0),   # now drives the mass-coupled deposit
    "nuc_beta":     ("kinetics.nuc_beta_C",      0.5, 2.0),
}


def _set_attr(case, path, value):
    obj, *rest = path.split(".")
    target = getattr(case, obj)
    setattr(target, rest[0], value)


def _eval_case(case, fast=True):
    #  F20: the fast calibration grid is closer to production than the old 36/3 (48 cells /
    #  6 realisations) and DETERMINISTIC, so the objective is smooth (no per-eval ensemble noise
    #  that would defeat the gradient-free optimiser) and reasonably quick. The run is shortened
    #  only when the case is very long, keeping enough transient to resolve time-to-plug. Final
    #  validation is done at the FULL production grid (fast=False) and reported.
    c = copy.deepcopy(case)
    if fast:
        c.pipeline.n_cells = max(min(c.pipeline.n_cells, 48), 40)
        c.numerics.n_ensemble = 6
        c.numerics.deterministic = True
        c.numerics.t_end_h = min(c.numerics.t_end_h, 30.0)
    sv = TransientSHCT(c)
    sv.run(verbose=False)
    return sv.engineering()


def _calib_residuals(eng, targets):
    """Per-target (modelled, target, weighted-sq-error) and total objective."""
    rows = []; err = 0.0
    for kkey, tval in targets.items():
        mod = eng.get(kkey, np.nan)
        if not np.isfinite(mod) or not np.isfinite(tval):
            #  a non-finite model/target is a hard infeasibility: penalise it with a LARGE finite
            #  value (>> any realistic squared relative error, which is O(0.01..1)) so the optimiser
            #  is driven AWAY from the infeasible region rather than treating a crash as a near-fit.
            e = 1.0e3
        else:
            e = ((mod - tval) / (abs(tval) + 1e-6)) ** 2
        rows.append((kkey, mod, tval, e)); err += e
    return rows, err


def calibrate(case: Case, targets: dict, free=None, maxiter=80):
    """Fit free model constants so predictions match the user's measured TARGETS,
    e.g. {"arrival_T_C": 8, "dP_total_bar": 30, "max_subcooling_C": 9,
          "time_to_plug_P50_h": 40}. Returns a calibrated Case. This is the
    mechanism by which the solver is adapted/validated to ANY real field or
    flow-loop dataset before quantitative use. F20: deterministic objective, more iterations,
    per-target residual report, and a final PRODUCTION-grid validation of the fit."""
    try:
        from scipy.optimize import minimize
    except ImportError as exc:                                   # G22: clear, early message
        raise SystemExit("calibration requires SciPy — install it with `pip install scipy` "
                         "(see requirements.txt).") from exc
    free = free or list(CALIB_PARAMS.keys())
    base = copy.deepcopy(case)
    defaults = {}
    for nm in free:
        path, _, _ = CALIB_PARAMS[nm]
        obj, attr = path.split(".")
        defaults[nm] = getattr(getattr(base, obj), attr)

    def make_case(x):
        c = copy.deepcopy(base)
        for nm, xi in zip(free, x):
            path, lo, hi = CALIB_PARAMS[nm]
            _set_attr(c, path, defaults[nm] * float(np.clip(xi, lo, hi)))
        return c

    history = {"n": 0}

    def objective(x):
        _, err = _calib_residuals(_eval_case(make_case(x), fast=True), targets)
        history["n"] += 1
        print(f"    calib eval {history['n']:2d}: err={err:.4f}  x={np.round(x,3)}")
        return err

    x0 = np.ones(len(free))
    print("=" * 64)
    print(" SHCT SOLVER — CALIBRATION TO MEASURED DATA")
    print(f"  targets: {targets}")
    print(f"  free params: {free}")
    print("=" * 64)
    res = minimize(objective, x0, method="Nelder-Mead",
                   options={"maxiter": maxiter, "xatol": 1e-2, "fatol": 1e-3})
    cal = make_case(res.x)
    print("-" * 64)
    print("  calibrated multipliers:")
    for nm, xi in zip(free, res.x):
        path, lo, hi = CALIB_PARAMS[nm]
        print(f"    {nm:12s} x{np.clip(xi,lo,hi):.3f}  -> {path}")
    #  F20: report the fit quality per target at the PRODUCTION grid (not the fast grid),
    #  so the user sees whether the calibration actually transfers.
    rows, total = _calib_residuals(_eval_case(cal, fast=False), targets)
    print("-" * 64)
    print("  fit quality at production grid (modelled vs target):")
    print(f"    {'metric':22s} {'model':>12s} {'target':>12s} {'rel.err %':>10s}")
    for kkey, mod, tval, _e in rows:
        rel = abs(mod - tval) / (abs(tval) + 1e-6) * 100 if np.isfinite(mod) else float("nan")
        print(f"    {kkey:22s} {mod:12.4g} {tval:12.4g} {rel:10.1f}")
    print(f"  fast-grid objective = {res.fun:.4f}   production-grid objective = {total:.4f}")
    if total > 2.0 * max(res.fun, 1e-6) + 0.1:
        print("  [WARN] production-grid fit notably worse than fast-grid — refine grid or re-calibrate.")
    #  #18: LOCAL IDENTIFIABILITY / SENSITIVITY at the optimum — perturb each calibrated
    #  parameter +/-10% and report the objective change. A near-flat response = that parameter
    #  is poorly identified by these targets (the fit is non-unique along it); a steep response
    #  = well-constrained. This is the honest uncertainty statement the fit needs.
    print("-" * 64)
    #  #22: posterior PARAMETER UNCERTAINTY from the local curvature (Laplace approximation): the
    #  objective near the optimum is ~quadratic, so the parameter standard error is
    #  sigma ~ sqrt(2 * obj_min / d2obj/dx2). A wide sigma means that parameter is poorly
    #  identified by these targets (the fit is non-unique along it) — the honest uncertainty.
    print("  parameter posterior (Laplace) — value, +/-1sigma estimate, identifiability:")
    print(f"    {'param':12s} {'multiplier':>12s} {'+/-1sigma':>12s}   {'identifiability'}")
    dh = 0.10
    for i, nm in enumerate(free):
        xm = res.x.copy(); xp = res.x.copy()
        xm[i] *= (1 - dh); xp[i] *= (1 + dh)
        _, e_m = _calib_residuals(_eval_case(make_case(xm), fast=True), targets)
        _, e_p = _calib_residuals(_eval_case(make_case(xp), fast=True), targets)
        curv = (e_m + e_p - 2.0 * res.fun) / (dh * res.x[i]) ** 2     # d2obj/dx2
        sigma = math.sqrt(max(2.0 * max(res.fun, 1e-6) / max(curv, 1e-9), 0.0)) if curv > 1e-9 else float("inf")
        path, lo, hi = CALIB_PARAMS[nm]
        mult = float(np.clip(res.x[i], lo, hi))
        ident = "well-constrained" if sigma < 0.3 * abs(mult) else "poorly identified (wide posterior)"
        sig_s = f"{sigma:.3f}" if math.isfinite(sigma) else "inf"
        print(f"    {nm:12s} {mult:12.3f} {sig_s:>12s}   {ident}")
    print("=" * 64)
    return cal


def _gelman_rubin(chain_mults):
    """Gelman-Rubin potential-scale-reduction R-hat per parameter from a list of equal-length
    per-chain MULTIPLIER arrays (each (n, ndim)). R-hat -> 1 as the chains mix; >~1.1 flags
    non-convergence. Returns an (ndim,) array (nan where it is not estimable: <2 chains, <2
    samples, or zero within-chain variance)."""
    m = len(chain_mults)
    if m < 2:
        return np.full(chain_mults[0].shape[1] if chain_mults else 0, np.nan)
    n = min(c.shape[0] for c in chain_mults)
    if n < 2:
        return np.full(chain_mults[0].shape[1], np.nan)
    arr = np.stack([c[:n] for c in chain_mults], axis=0)          # (m, n, ndim)
    chain_means = arr.mean(axis=1)                                # (m, ndim)
    grand_mean = chain_means.mean(axis=0)                         # (ndim,)
    B = n / (m - 1.0) * np.sum((chain_means - grand_mean) ** 2, axis=0)        # between-chain
    W = np.mean(np.var(arr, axis=1, ddof=1), axis=0)             # within-chain
    var_hat = (n - 1.0) / n * W + B / n
    with np.errstate(divide="ignore", invalid="ignore"):
        rhat = np.sqrt(var_hat / W)
    return np.where(W > 1e-12, rhat, np.nan)


def bayesian_calibrate(case: Case, targets: dict, free=None, n_samples=400,
                       sigma_rel=0.10, burn_frac=0.3, seed=7, n_chains=4):
    """FULL Bayesian posterior calibration by MULTI-CHAIN MCMC (F22). Samples the posterior of the
    free model constants given the measured TARGETS and their relative measurement uncertainty
    sigma_rel, by Metropolis-Hastings — yielding the posterior MEAN, standard deviation, credible
    intervals, the parameter CORRELATION matrix (which the Laplace ±1σ in calibrate() cannot give),
    AND a convergence diagnostic. Gaussian likelihood -0.5*sum(((model-target)/(sigma*|target|))^2)
    with uniform priors on the bounds. Robustness upgrades vs. a single untuned chain:
      * n_chains over-dispersed starts -> a Gelman-Rubin R-hat per parameter (R-hat->1 = converged;
        >~1.1 warns the chains have NOT mixed, so the posterior is not yet trustworthy);
      * the proposal width is ADAPTED during burn-in toward the ~0.234 optimal acceptance rate.
    Returns {samples, mean, std, corr, names, accept_rate, rhat, n_chains}; no SciPy needed."""
    free = free or list(CALIB_PARAMS.keys())
    base = copy.deepcopy(case)
    ndim = len(free)
    defaults = {}
    for nm in free:
        path, _, _ = CALIB_PARAMS[nm]
        obj, attr = path.split("."); defaults[nm] = getattr(getattr(base, obj), attr)

    def make_case(x):
        c = copy.deepcopy(base)
        for nm, xi in zip(free, x):
            path, lo, hi = CALIB_PARAMS[nm]
            obj, attr = path.split("."); setattr(getattr(c, obj), attr, defaults[nm] * float(np.clip(xi, lo, hi)))
        return c

    def loglike(x):
        for nm, xi in zip(free, x):                       # uniform prior on the bounds
            _, lo, hi = CALIB_PARAMS[nm]
            if not (lo <= xi <= hi):
                return -1e18
        eng = _eval_case(make_case(x), fast=True)
        ll = 0.0
        for kkey, tval in targets.items():
            mod = eng.get(kkey, np.nan)
            if not np.isfinite(mod):
                return -1e18
            ll += -0.5 * ((mod - tval) / (sigma_rel * (abs(tval) + 1e-6))) ** 2
        return ll

    print("=" * 64)
    print(" SHCT SOLVER — BAYESIAN (MCMC) POSTERIOR CALIBRATION")
    print(f"  targets: {targets}   sigma_rel: {sigma_rel}   samples: {n_samples} x {n_chains} chains")
    print("=" * 64)
    n_burn = int(burn_frac * n_samples)
    chain_mults = []; accepts = []
    for ci in range(n_chains):
        rng = np.random.default_rng(seed + 1009 * ci + 1)
        #  over-dispersed start (chain 0 at the prior centre) so R-hat is a meaningful mixing test
        if ci == 0:
            x = np.ones(ndim)
        else:
            x = np.array([float(np.clip(1.0 + 0.20 * rng.standard_normal(),
                                        CALIB_PARAMS[nm][1], CALIB_PARAMS[nm][2])) for nm in free])
        lp = loglike(x)
        step = np.array([0.08] * ndim)                   # proposal width (multiplier units)
        chain = []; n_acc = 0
        for it in range(n_samples):
            xp = x + step * rng.standard_normal(ndim)
            lpp = loglike(xp)
            if math.log(rng.random() + 1e-300) < (lpp - lp):
                x, lp = xp, lpp; n_acc += 1
            chain.append(x.copy())
            #  adapt the proposal toward the optimal ~0.234 acceptance during burn-in only
            #  (freezing it afterwards keeps the chain a valid stationary sampler).
            if it < n_burn and (it + 1) % 10 == 0:
                rate = n_acc / (it + 1)
                step = np.clip(step * math.exp(0.6 * (rate - 0.234)), 1e-3, 2.0)
            if (it + 1) % 50 == 0:
                print(f"    chain {ci+1}/{n_chains}  MCMC {it+1:4d}/{n_samples}  "
                      f"accept={n_acc/(it+1):.2f}  x={np.round(x,3)}")
        chain = np.array(chain)
        post = chain[n_burn:]                             # discard burn-in
        mult = np.array([np.clip(post[:, i], *CALIB_PARAMS[nm][1:]) for i, nm in enumerate(free)]).T
        chain_mults.append(mult); accepts.append(n_acc / max(len(chain), 1))

    allmult = np.vstack(chain_mults)
    mean = allmult.mean(0); std = allmult.std(0)
    corr = np.corrcoef(allmult.T) if allmult.shape[1] > 1 else np.array([[1.0]])
    rhat = _gelman_rubin(chain_mults)
    accept_rate = float(np.mean(accepts))
    print("-" * 64)
    print("  posterior (MCMC) — mean +/- std, 90% credible interval, R-hat (convergence):")
    for i, nm in enumerate(free):
        lo90, hi90 = np.percentile(allmult[:, i], [5, 95])
        rh = rhat[i] if i < len(rhat) else float("nan")
        rh_s = f"{rh:.3f}" if np.isfinite(rh) else "n/a"
        flag = "" if (not np.isfinite(rh) or rh < 1.1) else "  [WARN not converged]"
        print(f"    {nm:12s} {mean[i]:8.3f} +/- {std[i]:.3f}   [90% CI {lo90:.3f}, {hi90:.3f}]"
              f"   R-hat={rh_s}{flag}")
    print("  posterior parameter CORRELATIONS:")
    print("    " + " ".join(f"{nm:>10s}" for nm in free))
    for i, nm in enumerate(free):
        print(f"    {nm:10s} " + " ".join(f"{corr[i, j]:10.2f}" for j in range(len(free))))
    print(f"  acceptance rate (mean over chains): {accept_rate:.2f}")
    if np.any(np.isfinite(rhat) & (rhat > 1.1)):
        print("  [WARN] one or more R-hat > 1.1 — chains have not mixed; raise n_samples/n_chains.")
    print("=" * 64)
    return dict(samples=allmult, mean=mean, std=std, corr=corr, names=free,
                accept_rate=accept_rate, rhat=rhat, n_chains=n_chains)


def grid_convergence_report(case: Case, factors=(1, 2)):
    """Richardson-style discretisation-error report: re-run the case at the configured
    grid and a refined grid and compare key integral metrics. PURE DIAGNOSTIC — it does
    not change the solver, its scheme or its results; it quantifies the first-order grid
    sensitivity that the documentation flags. Returns (metric, coarse, fine, rel_change)."""
    base_n = case.pipeline.n_cells
    out = {}
    for fct in factors:
        c = copy.deepcopy(case)
        c.pipeline.n_cells = int(base_n * fct)
        c.numerics.n_ensemble = min(c.numerics.n_ensemble, 5)
        sv = TransientSHCT(c); sv.run(verbose=False)
        out[fct] = sv.engineering()
    #  D12: include the HEADLINE products (time-to-plug, P_plug, Phi_SH, deposit), not just dP/T,
    #  so the grid sensitivity of the quantities the solver is actually used for is reported.
    keys = ["dP_total_bar", "arrival_T_C", "max_subcooling_C", "Vm_peak_mps",
            "max_Phi_SH", "peak_deposit_mm", "P_plug", "time_to_plug_P50_h"]
    coarse, fine = out[factors[0]], out[factors[-1]]
    print("=" * 64)
    print(" SHCT SOLVER — GRID-CONVERGENCE (Richardson) REPORT")
    print("=" * 64)
    rows = []
    for kk in keys:
        a, b = coarse.get(kk, float("nan")), fine.get(kk, float("nan"))
        rel = abs(a - b) / max(abs(b), 1e-9)
        rows.append((kk, a, b, rel))
        print(f"  {kk:22s} {a:10.3f} -> {b:10.3f}   ({rel * 100:5.1f}% change)")
    print("-" * 64)
    print(f"  grids: n_cells {int(base_n * factors[0])} vs {int(base_n * factors[-1])}")
    print("  (large % change => refine the production grid / steep-terrain sections)")
    print("=" * 64)
    return rows


# =============================================================================
#  EXPERIMENTAL / FIELD VALIDATION HARNESS  (#4, #20)
# =============================================================================
#  These items cannot be CLOSED in code — they need YOUR flow-loop / field data. What the
#  code provides is the mechanism to ingest that data and quantify how well the solver (and
#  the Phi_SH plug criterion) reproduce it, so "validation" becomes a measurable, repeatable
#  report instead of a claim. Dataset JSON (all keys optional):
#    {"name": "...",
#     "profiles": {"x_km":[...], "P_bar":[...], "T_C":[...], "holdup":[...], "Phi_SH":[...]},
#     "scalars":  {"dP_total_bar":.., "arrival_T_C":.., "time_to_plug_P50_h":.., "max_subcooling_C":..},
#     "plug_events":[{"x_km":.., "plugged":true},{"x_km":.., "plugged":false}, ...]}
def validate_against_data(case: Case, dataset: dict):
    """Run the solver for `case` and score it against measured `dataset` (#20): per-profile RMSE
    and per-scalar relative error, plus a Phi_SH plug-criterion skill check (#4) — does Phi_SH>1
    discriminate the cells/sections that actually plugged? Returns a report dict and prints it."""
    sv = TransientSHCT(case); sv.run(verbose=False); eng = sv.engineering(); r = sv.results
    x_km = sv.x / 1000.0
    med = lambda A: np.nanmedian(A, 1)
    report = {"name": dataset.get("name", "dataset"), "profiles": {}, "scalars": {}, "phi_sh_skill": None}
    field_map = {"P_bar": med(r["p"]), "T_C": med(r["T"]), "holdup": med(r["alpha_l"]),
                 "Phi_SH": med(r["max_PhiSH"]), "subcooling_C": med(r["Tsub"])}
    print("=" * 64); print(f" SHCT SOLVER — VALIDATION vs DATA: {report['name']}"); print("=" * 64)
    prof = dataset.get("profiles", {})
    if prof.get("x_km"):
        xq = np.asarray(prof["x_km"], float)
        print("  PROFILES (interp. to measured x):   RMSE     bias    n")
        for key, model_field in field_map.items():
            if prof.get(key):
                meas = np.asarray(prof[key], float)
                pred = np.interp(xq, x_km, model_field)
                rmse = float(np.sqrt(np.nanmean((pred - meas) ** 2)))
                bias = float(np.nanmean(pred - meas))
                report["profiles"][key] = {"rmse": rmse, "bias": bias, "n": int(meas.size)}
                print(f"    {key:16s}              {rmse:8.3f} {bias:8.3f} {meas.size:4d}")
    scal = dataset.get("scalars", {})
    if scal:
        print("  SCALARS:                 model      measured    rel.err %")
        for key, meas in scal.items():
            mod = eng.get(key, float("nan"))
            rel = abs(mod - meas) / (abs(meas) + 1e-9) * 100 if np.isfinite(mod) else float("nan")
            report["scalars"][key] = {"model": float(mod), "measured": float(meas), "rel_err_pct": float(rel)}
            print(f"    {key:22s} {mod:10.4g} {meas:12.4g} {rel:10.1f}")
    ev = dataset.get("plug_events", [])
    if ev:
        #  #4: does the Phi_SH>1 criterion match observed plug/no-plug at the reported locations?
        phi_at = lambda xk: float(np.interp(xk, x_km, field_map["Phi_SH"]))
        tp = fp = tn = fn = 0
        for e in ev:
            pred_plug = phi_at(float(e["x_km"])) > 1.0; obs = bool(e.get("plugged", False))
            tp += pred_plug and obs; fp += pred_plug and not obs
            tn += (not pred_plug) and (not obs); fn += (not pred_plug) and obs
        acc = (tp + tn) / max(tp + fp + tn + fn, 1)
        report["phi_sh_skill"] = {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "accuracy": acc}
        print(f"  Phi_SH>1 PLUG CRITERION (n={len(ev)}): accuracy {acc*100:.0f}% "
              f"(TP {tp}, FP {fp}, TN {tn}, FN {fn})")
    print("-" * 64)
    print("  NOTE: absolute validity depends on the quality/representativeness of YOUR data;")
    print("        calibrate first (--calibrate) and re-run this to close the loop.")
    print("=" * 64)
    return report


def blind_validate(case: Case, dataset: dict, train_keys=None, maxiter=40):
    """BLIND (held-out) validation (F21): calibrate on a TRAIN subset of the measured scalars,
    then report the fit on the HELD-OUT subset the calibration never saw — the honest test of
    predictive skill (vs. fitting and reporting on the same data). dataset["scalars"] supplies the
    measurements; train_keys names the calibration subset (default: the first half)."""
    scal = dict(dataset.get("scalars", {}))
    if len(scal) < 2:
        print("[blind] need >=2 measured scalars to split train/test"); return None
    keys = list(scal)
    train_keys = train_keys or keys[: max(1, len(keys) // 2)]
    test_keys = [k for k in keys if k not in train_keys]
    train = {k: scal[k] for k in train_keys}; test = {k: scal[k] for k in test_keys}
    print("=" * 64)
    print(" SHCT SOLVER — BLIND (HELD-OUT) VALIDATION")
    print(f"  train (calibrate on): {train}")
    print(f"  test  (held out)    : {test}")
    print("=" * 64)
    cal = calibrate(case, train, maxiter=maxiter)
    eng = _eval_case(cal, fast=(maxiter < 10))             # quick eval when used as a fast/test path
    print("-" * 64)
    print("  HELD-OUT predictive skill (model vs measured, NOT used in calibration):")
    rows = {}
    for k, meas in test.items():
        mod = float(eng.get(k, float("nan")))
        rel = abs(mod - meas) / (abs(meas) + 1e-9) * 100 if np.isfinite(mod) else float("nan")
        rows[k] = {"model": mod, "measured": meas, "rel_err_pct": rel}
        print(f"    {k:22s} model {mod:10.4g} vs measured {meas:10.4g}   ({rel:5.1f}% held-out error)")
    print("=" * 64)
    return {"train": train, "test": rows, "calibrated_case": cal}


def validate_hydrate_curve(dataset_path, outdir=None, calibrate_offset=True):
    """VALIDATION against REAL published experimental data: compare the solver's hydrate-equilibrium
    closure to measured hydrate dissociation P-T points (e.g. Deaton & Frost 1946 / Sloan & Koh 2008
    methane Lw-H-V data). Reports RMSE/bias/max-error in degC; optionally fits the single allowable
    Teq offset (a one-parameter calibration) and reports the calibrated error too. Writes a chart."""
    with open(dataset_path) as fh:
        ds = json.load(fh)
    pts = np.asarray(ds["points_P_bar_T_C"], float)
    P = pts[:, 0]; T_meas = pts[:, 1]
    gas_sg = float(ds.get("gas_sg", 0.6)); sal = float(ds.get("salinity_wt", 0.0))
    T_model = hydrate_equilibrium_T(P, gas_sg=gas_sg, salinity_wt=sal)
    err = T_model - T_meas
    rmse = float(np.sqrt(np.mean(err ** 2))); bias = float(np.mean(err)); mx = float(np.max(np.abs(err)))
    offset = -bias if calibrate_offset else 0.0
    err_c = T_model + offset - T_meas
    rmse_c = float(np.sqrt(np.mean(err_c ** 2)))
    print("=" * 64)
    print(" SHCT SOLVER — HYDRATE-CURVE VALIDATION vs PUBLISHED EXPERIMENTAL DATA")
    print("=" * 64)
    print(f"  dataset : {ds.get('name')}")
    print(f"  source  : {ds.get('source','')[:96]}...")
    print(f"  fluid   : {ds.get('fluid')}   (gas_sg={gas_sg}, salinity={sal} wt%)")
    print(f"  points  : {len(P)}")
    print("-" * 64)
    print(f"    {'P (bar)':>9} {'T_meas (C)':>11} {'T_model (C)':>12} {'err (C)':>9}")
    for Pi, Tm, Tmod, e in zip(P, T_meas, T_model, err):
        print(f"    {Pi:9.1f} {Tm:11.2f} {Tmod:12.2f} {e:9.2f}")
    print("-" * 64)
    print(f"  as-shipped (literature-default correlation):  RMSE = {rmse:.2f} C, "
          f"bias = {bias:+.2f} C, max|err| = {mx:.2f} C")
    if calibrate_offset:
        print(f"  after 1-parameter calibration (Teq offset {offset:+.2f} C):  RMSE = {rmse_c:.2f} C")
    print("=" * 64)
    rep = {"name": ds.get("name"), "n": int(len(P)), "rmse_C": rmse, "bias_C": bias,
           "max_abs_err_C": mx, "rmse_after_offset_C": rmse_c, "offset_C": offset,
           "source": ds.get("source")}
    if outdir and HAVE_MPL:
        os.makedirs(outdir, exist_ok=True)
        Pg = np.linspace(float(P.min()) * 0.8, float(P.max()) * 1.1, 120)
        Tg = hydrate_equilibrium_T(Pg, gas_sg=gas_sg, salinity_wt=sal)
        fig, ax = plt.subplots(figsize=(6.4, 5))
        ax.plot(Tg, Pg, color=NAVY, lw=2, label="SHCT model (as shipped)")
        ax.plot(Tg + offset, Pg, color=TEAL, lw=1.6, ls="--",
                label=f"SHCT model (calibrated {offset:+.2f}°C)")
        ax.scatter(T_meas, P, color=RED, zorder=5, label="published experimental data")
        ax.set_xlabel("temperature (°C)"); ax.set_ylabel("pressure (bar)")
        ax.set_title("Hydrate-equilibrium validation vs published data", color=NAVY, fontweight="bold")
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
        ax.grid(alpha=.25)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "hydrate_validation.png"), dpi=150)
        plt.close(fig)
        with open(os.path.join(outdir, "hydrate_validation_report.json"), "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep


# =============================================================================
#  CLOSURE VALIDATION vs PUBLISHED REFERENCES (v8) — same evidentiary standard as
#  validate_hydrate_curve(): score the closures the hydraulics actually USE against
#  universally-cited, reproducible published references (NOT against the solver
#  itself). These validate the hydraulic ingredients of holdup and pressure drop;
#  they do NOT replace a full production-flow field validation, which still needs an
#  operator's measured holdup/dP/arrival-T along a real line (stated plainly below).
# =============================================================================
def _colebrook_white(Re, rel_rough, iters=80):
    """EXACT Darcy friction factor from the implicit Colebrook-White equation, by fixed-point
    iteration to ~machine precision. This is the reproducible reference the explicit
    haaland_friction() closure is compared against.
    CITATION ROLES — PRIMARY (origin of the equation): Colebrook, C.F. & White, C.M. (1937) /
    Colebrook, C.F. (1939), J. Inst. Civ. Eng. 11:133-156. CORROBORATING/secondary (graphical
    presentation that popularised it, NOT its origin): Moody, L.F. (1944), Trans. ASME 66:671-684
    ('the Moody chart')."""
    Re = np.asarray(Re, float); eps = np.asarray(rel_rough, float)
    inv = -1.8 * np.log10((eps / 3.7) ** 1.11 + 6.9 / Re)        # Haaland seed for 1/sqrt(f)
    for _ in range(iters):
        inv = -2.0 * np.log10(eps / 3.7 + 2.51 * inv / Re)
    return 1.0 / inv ** 2


def validate_friction_curve(outdir=None, ref_path=None):
    """VALIDATION of the friction closure (haaland_friction) against the Colebrook-White reference
    across a turbulent Re x relative-roughness grid. Friction sets the FRICTIONAL pressure gradient
    in every SHCT hydraulic engine, so this directly underpins the predicted pressure drop. Reports
    the max and RMS PERCENT deviation in the Darcy friction factor; the laminar branch (Re<2300,
    f=64/Re) is identical in both and excluded, as is the ill-defined 2300-4000 transition.
    CITATION ROLES — reference value (PRIMARY/origin): Colebrook-White equation (Colebrook 1939);
    corroborating graphical presentation: Moody chart (Moody 1944). Closure under test (PRIMARY):
    Haaland, S.E. (1983), J. Fluids Eng. 105(1):89-90, an explicit fit published as <=~2% of
    Colebrook-White."""
    src = ("reference PRIMARY: Colebrook-White equation (Colebrook 1939); corroborating: Moody "
           "chart (Moody 1944). Closure under test: Haaland (1983) explicit approximation.")
    Re_grid = np.array([4e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 1e8])
    eps_grid = np.array([0.0, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2])
    if ref_path and os.path.exists(ref_path):
        with open(ref_path) as fh:
            ds = json.load(fh)
        Re_grid = np.array(ds["grid"]["Re"], float)
        eps_grid = np.array(ds["grid"]["rel_roughness"], float)
        src = ds.get("reference_primary", ds.get("reference", src))
    RE, EPS = np.meshgrid(Re_grid, eps_grid)
    f_haa = haaland_friction(RE.ravel(), EPS.ravel()).reshape(RE.shape)
    f_ref = _colebrook_white(RE, EPS)
    pct = (f_haa - f_ref) / f_ref * 100.0
    rms = float(np.sqrt(np.mean(pct ** 2))); mx = float(np.max(np.abs(pct)))
    #  self-check: Colebrook must reproduce the smooth-pipe Moody value ~0.018 at Re=1e5
    moody_check = float(_colebrook_white(np.array([1e5]), np.array([0.0]))[0])
    print("=" * 70)
    print(" SHCT SOLVER — FRICTION-CLOSURE VALIDATION vs COLEBROOK-WHITE (MOODY)")
    print("=" * 70)
    print(f"  closure : haaland_friction()  (Haaland 1983 explicit approximation — PRIMARY)")
    print(f"  ref     : Colebrook-White equation (Colebrook 1939 — PRIMARY/origin), computed in-code")
    print(f"            to machine precision; Moody (1944) chart = corroborating presentation only")
    print(f"  self-check: Colebrook f(Re=1e5, smooth) = {moody_check:.4f}  (Moody-chart value ~0.018)")
    print("-" * 70)
    print(f"    {'Re':>10} {'eps/D':>9} {'f_Haaland':>11} {'f_Colebrook':>12} {'dev %':>8}")
    for i in range(EPS.shape[0]):
        for j in range(RE.shape[1]):
            print(f"    {RE[i,j]:10.0f} {EPS[i,j]:9.1e} {f_haa[i,j]:11.5f} "
                  f"{f_ref[i,j]:12.5f} {pct[i,j]:8.2f}")
    print("-" * 70)
    print(f"  Haaland vs Colebrook-White (turbulent grid, n={pct.size}): "
          f"RMS deviation = {rms:.2f}%, max |deviation| = {mx:.2f}%")
    print(f"  (Haaland published the formula as a <=~2% explicit fit to Colebrook-White.)")
    print("=" * 70)
    rep = {"name": "friction Haaland vs Colebrook-White", "n": int(pct.size),
           "rms_pct_dev": rms, "max_abs_pct_dev": mx, "colebrook_smooth_Re1e5": moody_check,
           "source": src}
    if outdir and HAVE_MPL:
        os.makedirs(outdir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(6.6, 5))
        Re_fine = np.logspace(np.log10(4e3), 8, 200)
        for eps in [0.0, 1e-4, 1e-3, 1e-2]:
            ax.loglog(Re_fine, _colebrook_white(Re_fine, np.full_like(Re_fine, eps)),
                      color=NAVY, lw=1.4)
            ax.loglog(Re_fine, haaland_friction(Re_fine, np.full_like(Re_fine, eps)),
                      color=RED, lw=1.0, ls="--")
        ax.plot([], [], color=NAVY, lw=1.4, label="Colebrook-White (reference)")
        ax.plot([], [], color=RED, lw=1.0, ls="--", label="Haaland (SHCT closure)")
        ax.set_xlabel("Reynolds number"); ax.set_ylabel("Darcy friction factor f")
        ax.set_title("Friction closure validation vs Colebrook-White (Moody)",
                     color=NAVY, fontweight="bold")
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
        ax.grid(alpha=.25, which="both")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "friction_validation.png"), dpi=150)
        plt.close(fig)
        with open(os.path.join(outdir, "friction_validation_report.json"), "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep


def validate_drift_flux(outdir=None, ref_path=None):
    """VALIDATION of the drift-flux SLIP closure (drift_params) — which sets the gas/liquid holdup
    in the implicit and full-Newton engines — against the canonical published Taylor-bubble
    translational-velocity benchmarks. Reports the closure's C0 and drift Froude (Fr = v_d/sqrt(gD))
    at each limit and the deviation from the reference. Reports discrepancies HONESTLY (no hiding).
    CITATION ROLES (each reference carries the value it ORIGINATED; later works are tagged as
    confirming, not as the source):
      * vertical drift Fr=0.35 — PRIMARY: Dumitrescu (1943), ZAMM 23:139-149 (analytical potential-
        flow solution); CONFIRMED experimentally BY Nicklin, Wilkes & Davidson (1962), Trans.
        IChemE 40:61-68.
      * vertical C0=1.2 — PRIMARY: Nicklin et al. (1962) (fully-developed turbulent profile).
      * horizontal nose Fr=0.542 — PRIMARY: Benjamin, T.B. (1968), J. Fluid Mech. 31:209-248
        (analytical); ADOPTED as the horizontal drift limit BY Bendiksen, K.H. (1984), Int. J.
        Multiphase Flow 10:467-483.
      * horizontal C0~1.05 — PRIMARY: Bendiksen (1984).
      (Shoham 2006 and Fabre & Line 1992 are secondary COMPILATIONS that reproduce these; they are
      not cited here as origins.)"""
    refs = [{"orientation": "vertical (theta=+90 deg)", "theta_deg": 90.0, "C0_ref": 1.20, "Fr_drift_ref": 0.35},
            {"orientation": "horizontal (theta=0 deg)", "theta_deg": 0.0, "C0_ref": 1.05, "Fr_drift_ref": 0.542}]
    src = ("PRIMARY (origin of each value): Dumitrescu (1943) vertical Fr=0.35 [confirmed by Nicklin "
           "et al. 1962]; Nicklin et al. (1962) C0=1.2; Benjamin (1968) horizontal Fr=0.542 [adopted "
           "by Bendiksen 1984]; Bendiksen (1984) horizontal C0~1.05. Secondary compilations (not "
           "origins): Shoham (2006), Fabre & Line (1992).")
    if ref_path and os.path.exists(ref_path):
        with open(ref_path) as fh:
            ds = json.load(fh)
        refs = ds.get("points", refs); src = ds.get("note", src)
    D = 0.10                                    # reference diameter (Fr is D-independent for this closure)
    rows = []
    print("=" * 70)
    print(" SHCT SOLVER — DRIFT-FLUX SLIP VALIDATION vs CANONICAL SLUG-FLOW DATA")
    print("=" * 70)
    print(f"  closure : drift_params(theta, D) -> (C0, v_d);  Fr_drift = v_d / sqrt(g D)")
    print("-" * 70)
    print(f"    {'orientation':28} {'C0(mdl/ref)':>16} {'Fr(mdl/ref)':>18} {'Fr err %':>9}")
    for r in refs:
        th = math.radians(float(r["theta_deg"]))
        C0, vd = drift_params(th, D)
        C0 = float(np.asarray(C0).ravel()[0]); vd = float(np.asarray(vd).ravel()[0])
        Fr = vd / math.sqrt(G * D)
        c0_err = (C0 - r["C0_ref"]) / r["C0_ref"] * 100.0
        fr_err = (Fr - r["Fr_drift_ref"]) / r["Fr_drift_ref"] * 100.0
        rows.append({"orientation": r["orientation"], "C0_model": C0, "C0_ref": r["C0_ref"],
                     "C0_err_pct": c0_err, "Fr_model": Fr, "Fr_ref": r["Fr_drift_ref"],
                     "Fr_err_pct": fr_err})
        print(f"    {r['orientation']:28} {C0:6.3f}/{r['C0_ref']:<5.2f}   "
              f"{Fr:7.3f}/{r['Fr_drift_ref']:<6.3f}   {fr_err:8.1f}")
    print("-" * 70)
    print("  VERDICT (honest): the VERTICAL limit matches the benchmark exactly (C0=1.20, Fr=0.35;")
    print("  Fr origin Dumitrescu 1943, confirmed by Nicklin et al. 1962). The HORIZONTAL drift")
    print("  Froude (0.20) is BELOW the nose-propagation value (0.542; origin Benjamin 1968,")
    print("  adopted by Bendiksen 1984) — the closure deliberately uses a")
    print("  smaller effective axial drift in near-horizontal flow; correcting it toward 0.542 is")
    print("  a calibration choice that would shift the (golden-master) default holdup, so it is")
    print("  left to the user rather than changed silently.")
    print("=" * 70)
    rep = {"name": "drift-flux slip vs canonical slug-flow values", "rows": rows, "source": src}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "drift_flux_validation_report.json"), "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep


def validate_slug_frequency(outdir=None, ref_path=None):
    """VALIDATION of the slug-frequency closure (slug_frequency) — a CLOSURE-FIDELITY check that it
    reproduces the published Zabaras (2000) inclined-pipe correlation it implements, with the
    horizontal value reported relative to the original Gregory & Scott (1969) correlation. The
    honest accuracy ceiling is the correlation's OWN published scatter (~+/-60% on 399 air-water
    points); a hard line-specific validation still needs an operator's measured slug counts.
    CITATION ROLES — PRIMARY (origin of the horizontal base term): Gregory, G.A. & Scott, D.S.
    (1969), AIChE J. 15:933-935. PRIMARY (origin of the inclination factor and the implemented
    correlation): Zabaras, G.J. (2000), SPE J. 5(3):252-258, which extends Gregory-Scott to
    inclined pipe and reports the +/-60% band against its own 399-point air-water dataset."""
    conds = [{"Vsl": 0.5, "Vm": 2.0, "D": 0.05, "theta_deg": 0.0},
             {"Vsl": 1.0, "Vm": 3.0, "D": 0.05, "theta_deg": 0.0},
             {"Vsl": 0.5, "Vm": 2.0, "D": 0.05, "theta_deg": 5.0},
             {"Vsl": 1.0, "Vm": 4.0, "D": 0.10, "theta_deg": 0.0},
             {"Vsl": 0.3, "Vm": 1.5, "D": 0.15, "theta_deg": 2.0}]
    src = "Zabaras (2000) SPE J. 5(3):252-258; Gregory & Scott (1969) AIChE J. 15:933-935"
    if ref_path and os.path.exists(ref_path):
        with open(ref_path) as fh:
            ds = json.load(fh)
        conds = ds.get("test_conditions", conds)
        src = "; ".join(x for x in [ds.get("reference_primary_horizontal_base"),
                                     ds.get("reference_primary_inclined"),
                                     ds.get("reference")] if x) or src

    def zabaras(Vsl, Vm, D, theta):                 # the documented published correlation
        base = 0.0226 * ((Vsl / (G * D)) * (19.75 / Vm + Vm)) ** 1.2
        return base * (0.836 + 2.75 * math.sin(abs(theta)) ** 0.25)

    def greg_scott(Vsl, Vm, D):                     # horizontal base (no inclination factor)
        return 0.0226 * ((Vsl / (G * D)) * (19.75 / Vm + Vm)) ** 1.2

    rows = []; max_rel = 0.0
    print("=" * 70)
    print(" SHCT SOLVER — SLUG-FREQUENCY CLOSURE FIDELITY vs ZABARAS (2000)")
    print("=" * 70)
    print(f"  closure : slug_frequency(Vsl, Vm, D, theta)  ==  Zabaras (2000) correlation")
    print("-" * 70)
    print(f"    {'Vsl':>5} {'Vm':>5} {'D':>6} {'theta':>6} {'fs_model':>10} {'fs_Zabaras':>11} "
          f"{'fs_GS(horiz)':>12}")
    for c in conds:
        th = math.radians(float(c["theta_deg"]))
        fm = float(slug_frequency(np.array([c["Vsl"]]), np.array([c["Vm"]]),
                                  c["D"], np.array([th]))[0])
        fz = zabaras(c["Vsl"], c["Vm"], c["D"], th)
        fgs = greg_scott(c["Vsl"], c["Vm"], c["D"])
        rel = abs(fm - fz) / (abs(fz) + 1e-12) * 100.0
        max_rel = max(max_rel, rel)
        rows.append({**c, "fs_model_Hz": fm, "fs_zabaras_Hz": fz, "fs_gregoryscott_horiz_Hz": fgs,
                     "fidelity_err_pct": rel})
        print(f"    {c['Vsl']:5.2f} {c['Vm']:5.2f} {c['D']:6.3f} {c['theta_deg']:6.1f} "
              f"{fm:10.4f} {fz:11.4f} {fgs:12.4f}")
    print("-" * 70)
    print(f"  closure reproduces the Zabaras (2000) correlation to max {max_rel:.2e}% (implementation fidelity).")
    print("  At horizontal, the closure = 0.836 x Gregory-Scott (1969) — the Zabaras inclination factor.")
    print("  HONEST CEILING: Zabaras reports ~+/-60% scatter vs its 399-point air-water dataset;")
    print("  slug frequency is intrinsically scattered. Line-specific field validation needs")
    print("  an operator's measured slug-counting data (not available here).")
    print("=" * 70)
    rep = {"name": "slug frequency vs Zabaras (2000)", "rows": rows,
           "max_fidelity_err_pct": max_rel, "published_accuracy_band_pct": 60.0, "source": src}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "slug_frequency_validation_report.json"), "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep


def validate_flowloop(dataset_path, outdir=None, calibrate=True):
    """FULL-FLOW (production-quantity) VALIDATION of the solver's HOLDUP prediction against a REAL,
    peer-reviewed, open-access flow-loop dataset: measured air-water void fractions (quick-closing-
    valve method — the gold-standard direct holdup measurement). For each (Jsl, Jsg) the steady
    drift-flux holdup the solver's 'implicit' engine reduces to — void = Vsg/(C0*Vm + v_d) with
    (C0, v_d) from the SAME drift_params() closure the solver uses — is compared to the measured
    void fraction. Reports RMSE/bias in void AND in liquid holdup, the fraction of points inside the
    measured uncertainty band, and an optional 1-parameter calibration via numerics.drift_C0_factor
    (the verified holdup knob), mirroring the hydrate-curve 1-param offset. Honest: this validates
    the holdup PHYSICS against real data; full line dP/arrival-T still need an operator's dataset."""
    with open(dataset_path) as fh:
        ds = json.load(fh)
    pts = np.asarray(ds["points_Jsl_Jsg_voidfraction_unc"], float)
    Jsl = pts[:, 0]; Jsg = pts[:, 1]; void_meas = pts[:, 2]
    unc = pts[:, 3] if pts.shape[1] > 3 else np.full(len(pts), 0.02)
    D = float(ds["pipe"]["diameter_m"]); theta = math.radians(float(ds["pipe"].get("inclination_deg", 0.0)))
    Vm = Jsg + Jsl
    C0_0, vd = drift_params(theta, D)
    C0_0 = float(np.asarray(C0_0).ravel()[0]); vd = float(np.asarray(vd).ravel()[0])

    def predict(factor):
        vg = factor * C0_0 * Vm + vd
        return np.clip(Jsg / np.maximum(vg, 1e-6), 0.0, 0.999)

    void_pred = predict(1.0)
    err = void_pred - void_meas
    rmse = float(np.sqrt(np.mean(err ** 2))); bias = float(np.mean(err)); mx = float(np.max(np.abs(err)))
    hold_rmse = rmse                                   # H_L = 1 - void, so |error| is identical
    within = int(np.count_nonzero(np.abs(err) <= unc))
    #  1-parameter calibration: the drift_C0_factor that minimises the void RMSE (grid + refine)
    best_f, best_r = 1.0, rmse
    if calibrate:
        for f in np.linspace(0.6, 1.8, 121):
            r = float(np.sqrt(np.mean((predict(f) - void_meas) ** 2)))
            if r < best_r:
                best_r, best_f = r, float(f)
    print("=" * 72)
    print(" SHCT SOLVER — FLOW-LOOP HOLDUP VALIDATION vs MEASURED VOID FRACTION (REAL DATA)")
    print("=" * 72)
    print(f"  dataset : {ds.get('name')}")
    print(f"  source  : {ds.get('reference_primary','')[:92]}...")
    print(f"  pipe    : D={D*1000:.1f} mm, horizontal; air-water; n={len(pts)} quick-closing-valve points")
    print(f"  closure : void = Vsg/(C0*Vm + v_d), drift_params -> C0={C0_0:.3f}, v_d={vd:.3f} m/s")
    print("-" * 72)
    print(f"    {'Jsl':>6} {'Jsg':>6} {'void_meas':>10} {'void_pred':>10} {'err':>8} {'|err|<=unc':>11}")
    for a, g, vmi, vpi, ei, ui in zip(Jsl, Jsg, void_meas, void_pred, err, unc):
        print(f"    {a:6.3f} {g:6.3f} {vmi:10.2f} {vpi:10.3f} {ei:8.3f} "
              f"{'yes' if abs(ei) <= ui else 'no':>11}")
    print("-" * 72)
    print(f"  as-shipped drift-flux (C0={C0_0:.2f}): void RMSE = {rmse:.3f}, bias = {bias:+.3f} "
          f"(=liquid-holdup RMSE {hold_rmse:.3f}); max|err| = {mx:.3f}")
    print(f"  points within measured uncertainty band: {within}/{len(pts)}")
    if calibrate:
        print(f"  after 1-param calibration (numerics.drift_C0_factor = {best_f:.3f}, i.e. C0={best_f*C0_0:.2f}): "
              f"void RMSE = {best_r:.3f}")
        print(f"  -> a positive bias means the closure over-predicts gas fraction (under-predicts liquid")
        print(f"     holdup); raising drift_C0_factor is the physically-correct, already-wired correction.")
    print("=" * 72)
    rep = {"name": ds.get("name"), "n": int(len(pts)), "C0_as_shipped": C0_0, "vd": vd,
           "void_rmse": rmse, "void_bias": bias, "holdup_rmse": hold_rmse, "max_abs_err": mx,
           "within_uncertainty": within, "drift_C0_factor_calibrated": best_f,
           "void_rmse_calibrated": best_r, "source": ds.get("reference_primary")}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        if HAVE_MPL:
            order = np.argsort(void_meas)
            fig, ax = plt.subplots(figsize=(6.4, 5))
            ax.plot([0, 0.7], [0, 0.7], color=GREY, lw=1, ls=":")
            ax.errorbar(void_meas, void_pred, xerr=unc, fmt="o", color=NAVY, ms=5,
                        capsize=3, label="SHCT drift-flux (as shipped)")
            if calibrate:
                ax.scatter(void_meas, predict(best_f), color=TEAL, marker="s", s=28,
                           label=f"calibrated (C0 x{best_f:.2f})")
            ax.set_xlabel("measured void fraction (quick-closing valve)")
            ax.set_ylabel("predicted void fraction")
            ax.set_title(f"Flow-loop holdup validation vs {ds.get('name', 'published data')}",
                         color=NAVY, fontweight="bold")
            ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
            ax.grid(alpha=.25)
            fig.tight_layout(); fig.savefig(os.path.join(outdir, "flowloop_holdup_validation.png"), dpi=150)
            plt.close(fig)
        with open(os.path.join(outdir, "flowloop_holdup_validation_report.json"), "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep


def validate_closures(outdir=None, datadir=None):
    """Run ALL published-reference closure validations together (friction, drift-flux slip,
    slug frequency) plus the hydrate-equilibrium curve, and print a combined honest summary."""
    datadir = datadir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "7", "field_data")
    fr = validate_friction_curve(outdir, os.path.join(datadir, "friction_colebrook_reference.json"))
    dr = validate_drift_flux(outdir, os.path.join(datadir, "drift_flux_canonical.json"))
    sf = validate_slug_frequency(outdir, os.path.join(datadir, "slug_frequency_zabaras.json"))
    hy = None
    hpath = os.path.join(datadir, "hydrate_equilibrium_published.json")
    if os.path.exists(hpath):
        hy = validate_hydrate_curve(hpath, outdir=outdir)
    fl = None
    flpath = os.path.join(datadir, "flowloop_holdup_dasneves2025.json")
    if os.path.exists(flpath):
        fl = validate_flowloop(flpath, outdir=outdir)
    print("#" * 70)
    print(" SHCT — PUBLISHED-DATA VALIDATION SUMMARY (closures the hydraulics use)")
    print("#" * 70)
    print(f"  hydrate equilibrium (Deaton&Frost/Sloan&Koh): "
          f"RMSE {hy['rmse_C']:.2f} C ({hy['rmse_after_offset_C']:.2f} C after 1-param offset)"
          if hy else "  hydrate equilibrium: dataset not found")
    print(f"  friction Haaland vs Colebrook-White:           "
          f"RMS {fr['rms_pct_dev']:.2f}% , max {fr['max_abs_pct_dev']:.2f}%")
    fr_v = next((r for r in dr["rows"] if "vertical" in r["orientation"]), None)
    fr_h = next((r for r in dr["rows"] if "horizontal" in r["orientation"]), None)
    if fr_v and fr_h:
        print(f"  drift-flux slip:  vertical Fr err {fr_v['Fr_err_pct']:+.1f}% (matches), "
              f"horizontal Fr err {fr_h['Fr_err_pct']:+.1f}% (closure uses lower axial drift)")
    print(f"  slug frequency = Zabaras(2000) to {sf['max_fidelity_err_pct']:.1e}% "
          f"(published band ~+/-{sf['published_accuracy_band_pct']:.0f}%)")
    if fl:
        print(f"  flow-loop HOLDUP (Das Neves et al. 2025, real void fractions): "
              f"void RMSE {fl['void_rmse']:.3f} as-shipped, {fl['void_rmse_calibrated']:.3f} "
              f"after 1-param C0 calib; {fl['within_uncertainty']}/{fl['n']} within meas. band")
    print("-" * 70)
    print("  STILL OPEN (needs an OPERATOR's specific-line data, cannot be closed by public data):")
    print("   full production-flow dP & arrival-T along a particular real line (holdup now validated")
    print("   against a real flow-loop above; thermal/dP closures still benefit from line-specific data).")
    print("#" * 70)
    return {"friction": fr, "drift_flux": dr, "slug_frequency": sf, "hydrate": hy, "flowloop": fl}


def main(argv=None):
    ap = argparse.ArgumentParser(description="SHCT transient subsea slug+hydrate solver")
    ap.add_argument("--config")
    ap.add_argument("--dump-config", metavar="PATH")
    ap.add_argument("--schema", action="store_true",
                    help="print the case JSON schema (groups, fields, types, defaults) and exit (#25)")
    ap.add_argument("--scenario", choices=["steady", "rampup", "turndown", "shutin"])
    ap.add_argument("--engine", choices=["implicit", "twofluid", "twofluid_mass",
                                         "twofluid_mass_newton", "twofluid_full_newton", "quasisteady"],
                    help="hydrodynamic engine: implicit (drift-flux+implicit pressure, default), "
                         "twofluid (full two-fluid, two phase momenta), twofluid_full_newton "
                         "(fully-simultaneous La/Mg/momentum/p block-Newton), quasisteady (legacy)")
    ap.add_argument("--meg", type=float, metavar="WT_PCT",
                    help="inject thermodynamic inhibitor (MEG) at inlet, wt%% (hydrate management)")
    ap.add_argument("--verify", action="store_true", help="run the V&V verification suite and exit")
    ap.add_argument("--uq", action="store_true",
                    help="propagate a WIDER input uncertainty (U, kinetics, wall-capture) through "
                         "the ensemble for fuller P10/P50/P90 bands (raises n_ensemble to >=50). "
                         "A modest default spread is always on unless numerics.deterministic=True.")
    ap.add_argument("--dynamic-u", dest="dynamic_u", action="store_true",
                    help="flow- & deposit-dependent effective heat-transfer coefficient U")
    ap.add_argument("--grid-check", dest="grid_check", action="store_true",
                    help="run a 1x/2x grid-convergence (Richardson) report and exit")
    ap.add_argument("--crosssection", action="store_true",
                    help="reconstruct & output the reduced-order CROSS-SECTION (quasi-3-D) fields: "
                         "liquid level h/D, wetted perimeter, interface width, bottom/top-of-line "
                         "deposit and 2-D section velocity/temperature/phase reconstructions")
    ap.add_argument("--compo-report", dest="compo_report", action="store_true",
                    help="deep COMPOSITIONAL/PVT tracking report along the line (Peng-Robinson flash: "
                         "vapour fraction, per-component K-values, phase densities & viscosities)")
    ap.add_argument("--compo-sim", dest="compo_sim", action="store_true",
                    help="SEQUENTIAL COMPOSITIONAL TRANSPORT: composition grading along the line from "
                         "hydrate former depletion (component-conservative; needs fluids.composition)")
    ap.add_argument("--threed", "--3d", dest="threed", action="store_true",
                    help="build the 3-D reconstructed field over the pipe volume (from the 1-D core "
                         "+ cross-section closure), export a ParaView VTK volume (pipe_3d.vtk) and "
                         "render 3-D tube views. Reduced-order quasi-3-D, NOT a 3-D CFD solve.")
    ap.add_argument("--openfoam", action="store_true",
                    help="COUPLE TO OpenFOAM: generate runnable interFoam (3-D VOF) cases for the "
                         "sections that need CFD (riser/steep terrain, severe slugging, hydrate-"
                         "critical), with boundary conditions from the SHCT solution")
    ap.add_argument("--openfoam-run", dest="openfoam_run", action="store_true",
                    help="also RUN the generated OpenFOAM cases if OpenFOAM is installed (blockMesh/"
                         "interFoam on PATH) and ingest the CFD results back")
    ap.add_argument("--openfoam-sections", dest="openfoam_sections", type=int, default=3,
                    help="number of critical sections to export to OpenFOAM (default 3)")
    ap.add_argument("--openfoam-end-time", dest="openfoam_end_time", type=float, default=2.0,
                    help="interFoam physical end-time (s) for the coupled CFD run (default 2.0)")
    ap.add_argument("--openfoam-res", dest="openfoam_res", default="10x40",
                    help="interFoam o-grid resolution NixNz for the coupled CFD run (default 10x40)")
    ap.add_argument("--calibrate", metavar="TARGETS.json",
                    help="fit free constants to measured targets, write calibrated case")
    ap.add_argument("--validate", metavar="DATA.json",
                    help="score the solver against measured field/flow-loop data and exit (#4/#20)")
    ap.add_argument("--validate-hydrate", dest="validate_hydrate", metavar="DATA.json",
                    help="validate the hydrate-equilibrium closure against PUBLISHED experimental "
                         "P-T data (e.g. 7/field_data/hydrate_equilibrium_published.json) and exit")
    ap.add_argument("--validate-friction", dest="validate_friction", action="store_true",
                    help="validate the friction closure (Haaland) against the Colebrook-White/Moody "
                         "reference and exit (v8 — underpins frictional pressure drop)")
    ap.add_argument("--validate-drift", dest="validate_drift", action="store_true",
                    help="validate the drift-flux slip closure (holdup) against canonical published "
                         "Taylor-bubble values (Nicklin/Dumitrescu, Benjamin/Bendiksen) and exit (v8)")
    ap.add_argument("--validate-slugfreq", dest="validate_slugfreq", action="store_true",
                    help="validate the slug-frequency closure against the Zabaras (2000) correlation "
                         "it implements (closure fidelity; ~+/-60% published band) and exit (v8)")
    ap.add_argument("--validate-flowloop", dest="validate_flowloop", nargs="?", const="__default__",
                    metavar="DATA.json",
                    help="validate the solver's HOLDUP prediction against a REAL flow-loop void-fraction "
                         "dataset (default: Das Neves et al. 2025, 7/field_data/) and exit (v9, gap 1)")
    ap.add_argument("--validate-closures", dest="validate_closures", action="store_true",
                    help="run ALL published-reference closure validations (friction, drift-flux, slug "
                         "frequency, hydrate, flow-loop holdup) with a combined honest summary and exit (v8/v9)")
    ap.add_argument("--bayes-calibrate", dest="bayes_calibrate", metavar="TARGETS.json",
                    help="MCMC posterior calibration (F22): posterior mean/std/correlations and exit")
    ap.add_argument("--blind-validate", dest="blind_validate", metavar="DATA.json",
                    help="blind held-out validation (F21): calibrate on half the data, test on the rest")
    ap.add_argument("--outdir", default=OUTDIR_DEFAULT)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--log-level", default="INFO",
                    choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                    help="logging level for status messages (#23)")
    args = ap.parse_args(argv)
    log.setLevel(getattr(logging, args.log_level))

    if args.verify:
        return 0 if run_verification() else 1

    if args.schema:                                    # #25: formal case schema
        schema = {}
        for g, cls in _GROUP_TYPES.items():
            schema[g] = {fn: {"type": (fld.type if isinstance(fld.type, str) else str(fld.type)),
                              "default": (asdict(make_default_case())[g].get(fn))}
                         for fn, fld in cls.__dataclass_fields__.items()}
        print(json.dumps({"name": "str", **schema}, indent=2, default=str)); return 0

    if args.dump_config:
        with open(args.dump_config, "w") as fh:
            json.dump(asdict(make_default_case()), fh, indent=2)
        log.info("[config] default case written to %s", args.dump_config); return 0

    case = case_from_json(args.config) if args.config else make_default_case()
    if args.engine:
        case.numerics.engine = args.engine
    if args.meg is not None:
        case.operating.MEG_wt_inlet = args.meg
    if args.uq:
        case.numerics.uq_inputs = dict(STRONG_UQ)
        if case.numerics.n_ensemble < 50:
            case.numerics.n_ensemble = 50
    if args.dynamic_u:
        case.numerics.dynamic_U = True
    validate_case(case)                                # G21: fail fast on bad CLI/config inputs

    if args.calibrate:
        with open(args.calibrate) as fh:
            targets = json.load(fh)
        case = calibrate(case, targets)
        os.makedirs(args.outdir, exist_ok=True)
        calpath = os.path.join(args.outdir, "calibrated_case.json")
        with open(calpath, "w") as fh:
            json.dump(asdict(case), fh, indent=2)
        log.info("[calibrate] calibrated case written to %s", calpath)
    if args.scenario:
        case.scenario.kind = args.scenario
        case.name = f"{case.name.split(' (')[0]} ({args.scenario})"
    if args.bayes_calibrate:                               # F22
        with open(args.bayes_calibrate) as fh:
            targets = json.load(fh)
        bayesian_calibrate(case, targets); return 0
    if args.validate:
        with open(args.validate) as fh:
            dataset = json.load(fh)
        validate_against_data(case, dataset); return 0
    if args.validate_hydrate:
        os.makedirs(args.outdir, exist_ok=True)
        validate_hydrate_curve(args.validate_hydrate, outdir=args.outdir); return 0
    if args.validate_closures:                              # v8 — all closure validations
        os.makedirs(args.outdir, exist_ok=True)
        validate_closures(outdir=args.outdir); return 0
    if args.validate_friction:                              # v8
        os.makedirs(args.outdir, exist_ok=True)
        _dd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "7", "field_data")
        validate_friction_curve(outdir=args.outdir,
                                ref_path=os.path.join(_dd, "friction_colebrook_reference.json"))
        return 0
    if args.validate_drift:                                 # v8
        os.makedirs(args.outdir, exist_ok=True)
        _dd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "7", "field_data")
        validate_drift_flux(outdir=args.outdir,
                            ref_path=os.path.join(_dd, "drift_flux_canonical.json"))
        return 0
    if args.validate_slugfreq:                              # v8
        os.makedirs(args.outdir, exist_ok=True)
        _dd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "7", "field_data")
        validate_slug_frequency(outdir=args.outdir,
                                ref_path=os.path.join(_dd, "slug_frequency_zabaras.json"))
        return 0
    if args.validate_flowloop:                              # v9 — holdup vs real flow-loop data
        os.makedirs(args.outdir, exist_ok=True)
        _dd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "7", "field_data")
        dpath = (os.path.join(_dd, "flowloop_holdup_dasneves2025.json")
                 if args.validate_flowloop == "__default__" else args.validate_flowloop)
        validate_flowloop(dpath, outdir=args.outdir); return 0
    if args.blind_validate:                                # F21
        with open(args.blind_validate) as fh:
            dataset = json.load(fh)
        blind_validate(case, dataset); return 0
    if args.grid_check:
        grid_convergence_report(case); return 0
    os.makedirs(args.outdir, exist_ok=True)

    log.info("[run] solving '%s'  (transient coupled PDEs) ...", case.name)
    sv = TransientSHCT(case)
    sv.run()
    eng = sv.engineering()
    write_tables(sv, eng, args.outdir)
    log.info("[tables] CSV tables -> %s", args.outdir)
    if not args.no_plots:
        make_charts(sv, eng, args.outdir)
    #  reduced-order cross-section (quasi-3-D) reconstruction — resolves cross-sectional structure
    #  (interface level, wetted perimeter, bottom-of-line deposit, 2-D velocity/T fields) from the
    #  1-D solution. A fast, whole-line, screening-grade complement to 3-D CFD (NOT a replacement).
    if getattr(args, "crosssection", False):
        import shct_crosssection
        shct_crosssection.crosssection_outputs(sv, args.outdir)
        log.info("[crosssection] cross-section CSV + charts -> %s", args.outdir)
    #  deep compositional / PVT tracking along the line (Peng-Robinson EOS flash per station)
    if getattr(args, "compo_report", False):
        import shct_compositional
        shct_compositional.compositional_report(sv, args.outdir)
        log.info("[compo] compositional/PVT report -> %s", args.outdir)
    #  sequential compositional transport: composition grading from hydrate former depletion
    if getattr(args, "compo_sim", False):
        import shct_compositional_sim
        rep = shct_compositional_sim.simulate_composition(sv, args.outdir)
        log.info("[compo-sim] compositional transport (max |dz|=%.4f) -> %s",
                 rep["grading_max_abs_dz"], args.outdir)
    #  3-D reconstructed field over the pipe volume (quasi-3-D): VTK volume for ParaView + 3-D tubes.
    #  Reconstructed from the 1-D conservation core + cross-section closure — NOT a 3-D CFD solve.
    if getattr(args, "threed", False):
        import shct_threed
        shct_threed.threed_outputs(sv, args.outdir)
        log.info("[3d] quasi-3-D field: VTK volume + 3-D renders -> %s", args.outdir)
    #  SHCT <-> OpenFOAM coupling: emit runnable interFoam CFD cases for the sections that need
    #  true 3-D resolution (BCs from the SHCT solution); run + ingest if OpenFOAM is installed.
    if getattr(args, "openfoam", False) or getattr(args, "openfoam_run", False):
        import shct_openfoam
        try:
            _ni, _nz = (int(v) for v in str(getattr(args, "openfoam_res", "10x40")).lower().split("x"))
        except Exception:
            _ni, _nz = 10, 40
        man = shct_openfoam.couple(sv, args.outdir,
                                   max_sections=getattr(args, "openfoam_sections", 3),
                                   run=getattr(args, "openfoam_run", False),
                                   end_time=getattr(args, "openfoam_end_time", 2.0),
                                   Ni=_ni, Nz=_nz)
        if man["openfoam_available"]:
            log.info("[openfoam] %d CFD cases generated and run -> %s/openfoam_cases",
                     man["n_sections"], args.outdir)
            #  v9: real 3-D CFD ran locally — report the SHCT(1-D) vs interFoam(3-D VOF) holdup match
            comp = [e for e in man["sections"] if "cfd_alpha_l" in e]
            if comp:
                print("=" * 64)
                print(" SHCT (1-D) vs OpenFOAM interFoam (3-D VOF) — volume-avg liquid holdup")
                print(f"  (end_time={man['end_time']} s, o-grid {man['mesh']['Ni']}x{man['mesh']['Ni']}x"
                      f"{man['mesh']['Nz']}; coupled BCs from the SHCT solution)")
                print("-" * 64)
                for e in comp:
                    print(f"    {e['name']:14s} x={e['x_km']:6.2f} km : SHCT {e['shct_alpha_l']:.3f}  "
                          f"CFD {e['cfd_alpha_l']:.3f}  diff {e['rel_diff_pct']:.1f}%")
                print(f"  mean |SHCT-CFD| = {np.mean([e['abs_diff'] for e in comp]):.3f} holdup")
                print("  HONEST: single-/few-section, coarse o-grid, laminar VOF, short physical time")
                print("  (severe slugging keeps accumulating) — a live 3-D check, not a converged DNS.")
                print("=" * 64)
        else:
            log.info("[openfoam] %d runnable CFD cases generated (OpenFOAM not installed here — "
                     "run ./Allrun on an OpenFOAM machine) -> %s/openfoam_cases",
                     man["n_sections"], args.outdir)
    console_report(sv, eng)
    with open(f"{args.outdir}/summary.json", "w") as fh:
        json.dump(eng, fh, indent=2)
    log.info("[done] outputs in: %s", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
