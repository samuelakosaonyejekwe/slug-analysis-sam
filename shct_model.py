#!/usr/bin/env python3
# =============================================================================
#  shct_model.py — SHCT input data model (dataclasses) (#21)
# -----------------------------------------------------------------------------
#  Extracted from solver.py so the case schema is importable and validatable on
#  its own. Pure dataclasses; no numpy/solver dependency.
# =============================================================================
from __future__ import annotations
from dataclasses import dataclass, field


# =============================================================================
#  INPUT DATA MODEL
# =============================================================================
@dataclass
class Pipeline:
    length_m: float = 20_000.0
    diameter_m: float = 0.3048
    roughness_m: float = 4.5e-5
    n_cells: int = 100
    elevation_m: list | None = None       # length >= n_cells; default profile if None
    #  Multi-layer wall (steel + insulation + coating ...) for proper heat transfer
    #  & cooldown. Each layer: [thickness_m, conductivity_WmK, volumetric_heatcap_Jm3K].
    #  If None, the constant Operating.U_wall and a default steel+insulation mass are used.
    wall_layers: list | None = None
    h_inner: float = 1500.0               # inner film coeff (W/m2K)
    h_outer: float = 250.0                # outer (seawater) film coeff (W/m2K)


@dataclass
class Fluids:
    rho_oil: float = 820.0
    rho_water: float = 1025.0
    water_cut: float = 0.30
    mu_liquid: float = 2.0e-3
    mu_gas: float = 1.4e-5
    sigma: float = 0.020
    gas_MW: float = 0.019
    gas_Z: float = 0.90
    cp_liquid: float = 2100.0
    cp_gas: float = 2300.0
    L_hyd: float = 4.5e5                   # J/kg hydrate latent heat (exothermic formation)
    rho_hyd: float = 920.0                 # hydrate (sI/sII) bulk density (kg/m3)
    hyd_water_massfrac: float = 0.86       # mass fraction of a hydrate that is water (rest = gas)
    interfacial_Ci: float = 0.012         # interfacial friction factor (two-fluid slip)
    #  --- two-fluid closure constants (were hard-coded; exposed & documented) ---
    gas_wall_drag_frac: float = 0.10       # gas wall-shear vs liquid: gas wets a small fraction of
    #                                        the perimeter in stratified/slug flow (Taitel-Dukler geometry)
    vm_coeff: float = 1.2                  # interfacial-/virtual-mass pressure coeff (Bestion/RELAP)
    #                                        that renders the two-fluid system hyperbolic (real chars.)
    #  --- hydrate-equilibrium accuracy (closure inputs; defaults = natural-gas curve) ---
    gas_sg: float = 0.60                   # gas specific gravity (air=1); shifts hydrate Teq
    salinity_wt: float = 0.0               # produced-water salinity (wt% NaCl-eq); depresses Teq
    hyd_Teq_table: list | None = None      # optional [[P_bar,T_C],...] measured/PVT hydrate curve
    #                                        (overrides the correlation when provided -> exact, universal)
    #  --- gas compressibility (#1): Z either constant (gas_Z) or a Standing-Katz-style
    #  correlation Z(P,T) (gas_Z_corr=True); or a user [[P_bar,T_C,Z],...] table (exact). ---
    gas_Z_corr: bool = False               # True -> real-gas Z(P,T) instead of the constant gas_Z
    gas_Z_table: list | None = None        # optional [[P_bar, T_C, Z], ...] measured/EOS Z surface
    #  --- 3-phase oil/water slip (#2): off -> single composite liquid (back-compatible);
    #  on -> the water settles relative to oil (gravity/inclination), transported separately. ---
    oil_water_slip: bool = False
    mu_oil: float = 2.0e-3                  # used only when oil_water_slip; else mu_liquid is used
    mu_water: float = 0.8e-3
    #  --- droplet entrainment (#6): optional gas-core liquid-droplet fraction in high-shear
    #  (annular/churn) flow; lowers the effective wall-film holdup. Off by default. ---
    droplet_entrainment: bool = False
    droplet_field: bool = False            # A2: transport a SEPARATE conserved liquid-droplet field
    #                                        (gas-core droplets advected at the gas velocity, exchanging
    #                                        with the wall film via entrainment/deposition)
    water_drift: bool = False              # A3: water phase carries its own buoyancy/drag DRIFT velocity
    #                                        (force balance) rather than a fixed settling slip
    #  --- deeper PVT (#1): real-gas viscosity (Lee) and P,T-dependent oil density; or a full
    #  user PVT property table that overrides the correlations -> the EOS/black-oil enabler. ---
    gas_visc_corr: bool = False            # True -> Lee gas viscosity mu_g(rho_g,T,MW)
    oil_pvt_corr: bool = False             # True -> thermal-expansion/compressibility oil density
    pvt_table: list | None = None          # optional [[P_bar,T_C, rho_oil, rho_gas, mu_oil, mu_gas],...]
    #                                        property surface (nearest-point), for a measured/EOS PVT
    #  --- full compositional EOS (#A1): a {component: mole_fraction} makeup. When given, a
    #  Peng-Robinson flash (shct_eos) builds the PVT surface above AND sets the gas SG/Z from
    #  the EOS -> genuine compositional thermodynamics, not just a correlation. ---
    composition: dict | None = None
    n_pseudo: int = 1                      # B11: split C7+ into this many characterized pseudo-components
    #                                        (Whitson gamma + Kesler-Lee) before the flash; 1 = lumped
    MW_plus: float = 140.0                 # B11: average molar mass (g/mol) of the C7+ plus-fraction
    L_condensation: float = 3.0e5          # B9: hydrocarbon condensation/vaporisation latent heat (J/kg)
    condensation_latent: bool = False      # B9: include the condensation/evaporation latent heat in
    #                                        the energy balance (two-phase region buffers cooling)
    #  --- closure constants previously hard-coded inline (now exposed; defaults unchanged) ---
    k_liquid: float = 0.6                   # liquid thermal conductivity (W/mK) — sets the inner-film
    #                                        Prandtl/Nusselt in the advanced/dynamic-U heat-transfer path
    water_droplet_d_m: float = 5.0e-4       # water-droplet diameter for the Newton-regime settling drag
    water_droplet_Cd: float = 0.44         # water-droplet drag coefficient (Newton regime)
    settling_slip_coeff: float = 0.05      # prefactor on the inclination-scaled water settling slip
    slurry_visc_exp: float = -2.5          # Camargo-Palermo/Krieger-Dougherty relative-viscosity exponent
    api14e_C_factor: float = 122.0         # API-14E erosional-velocity C-factor (service-dependent)
    wax_appearance_C: float = -273.15      # wax appearance temperature (WAT, degC). Default off
    #                                        (absolute zero); set to the fluid's WAT to screen wax risk.


@dataclass
class Operating:
    q_liquid_insitu: float = 0.035
    q_gas_insitu_inlet: float = 0.140
    P_inlet_bar: float = 150.0
    T_inlet_C: float = 55.0
    T_seabed_C: float = 4.0
    U_wall: float = 22.0                  # used if pipeline.wall_layers is None
    MEG_wt_inlet: float = 0.0            # thermodynamic inhibitor injected at inlet (wt%)
    #  B9: how MEG_wt_inlet is interpreted. "aqueous" (default) = wt% IN the produced water
    #  (the basis Nielsen-Bucklin/Hammerschmidt use). "stream" = wt% of the whole liquid
    #  stream -> converted to the aqueous basis using the water cut before transport.
    MEG_basis: str = "aqueous"
    MEG_design_margin_C: float = 2.0     # subcooling design margin added before MEG sizing


@dataclass
class Kinetics:
    kg0: float = 6.0e-7
    growth_exp_n: float = 1.0
    phi_max: float = 0.55
    Ea_over_R: float = 1200.0
    D_phi: float = 5.0e-3                  # phase-field diffusivity (m2/s, numerical/physical)
    nuc_tau0_h: float = 0.6
    nuc_beta_C: float = 14.0
    nuc_dTsub_min: float = 1.0
    k_dep: float = 1.2e-4                   # legacy deposition gain (retained for back-compat;
    #                                        the deposit is now mass-coupled to hydrate growth, see wall_capture_eff)
    wall_capture_eff: float = 1.0          # efficiency with which wall-adjacent hydrate growth
    #                                        consolidates into the wall deposit (couples delta <-> phi mass)
    k_ero: float = 2.0e-3
    delta_max_frac: float = 0.92
    consol_restriction: float = 0.18
    C_phi: float = 1500.0
    #  Interface-renewal floor sets the SCALE of Phi_SH in stratified/stagnant flow
    #  (where slugs no longer renew the interface). Default 1e-4 preserves the original
    #  dynamics exactly; expose it so a case can use a physically-motivated minimum.
    f_slug_floor_Hz: float = 1.0e-4
    phi_internal_cap: float = 50.0         # cap on the Phi_SH that DRIVES the deposition core (unchanged)
    phi_report_cap: float = 1.0e4          # finite cap for the REPORTED Phi_SH (plot/IO only)
    #  --- deposit-insulation feedback (#7): as the wall deposit thickens it insulates the
    #  pipe wall, so the wall warms toward the bulk fluid and the wall subcooling that drives
    #  further deposition FALLS. k_dep_insul in [0,1] sets how strongly a full-bore deposit
    #  warms the wall (0 = old behaviour: wall pinned at T_seabed; 1 = wall reaches bulk T). ---
    k_dep_insul: float = 0.6
    deposit_k_hyd: float = 0.5             # hydrate-deposit thermal conductivity (W/mK), insulation
    #  --- kinetics/closure constants previously hard-coded inline (now exposed; defaults unchanged) ---
    T_ref_K: float = 283.15               # Arrhenius reference temperature for the growth-rate kg(T)
    nuc_tau_exp_cap: float = 50.0         # overflow guard on the induction-time exponent (numerical)
    wall_capture_Tsub_ref_C: float = 8.0  # wall subcooling that saturates the wall-capture availability
    slug_body_holdup_base: float = 0.70   # sub-grid slug-body holdup intercept (Dukler-Hubbard)
    slug_body_holdup_slope: float = 0.25  # sub-grid slug-body holdup slope on vsl/j
    slug_fraction_ref_Hz: float = 0.05    # reference frequency in the slug-fraction weight fslug/(fslug+ref)


@dataclass
class Numerics:
    t_end_h: float = 40.0
    cfl: float = 0.40
    dt_max_h: float = 0.05
    n_ensemble: int = 20
    seed: int = 7
    monitor_frac: float = 0.92
    plug_restriction_trip: float = 0.85
    n_snapshots: int = 80
    engine: str = "implicit"             # "implicit" = semi-implicit transient momentum
    #                                      "quasisteady" = legacy marched pressure (fallback)
    #  --- accuracy options (opt-in; default None/False reproduce the base results exactly) ---
    uq_inputs: dict | None = None        # per-realisation input uncertainty {param: rel_sigma},
    #                                      e.g. {"kg0":0.5,"nuc_tau0_h":0.5,"U_wall":0.15} -> honest bands
    dynamic_U: bool = False              # recompute effective-U from local flow & deposit each step
    #  --- ensemble realism (C10/C11): a MODEST per-realisation parameter spread is applied by
    #  default so the P10/P50/P90 bands are genuinely probabilistic (induction time, heat
    #  transfer and growth are physically uncertain) and the N ensemble columns are distinct
    #  realisations rather than identical duplicates. Set deterministic=True (or supply explicit
    #  uq_inputs) to override. --uq widens it further.
    deterministic: bool = False          # True -> no default spread (every realisation identical)
    #  --- numerical-scheme options (default values reproduce the verified base scheme) ---
    flux_limiter: str = "minmod"         # #8 TVD limiter for holdup/energy/phi advection:
    #                                      "upwind" (1st order) | "minmod" | "vanleer" | "superbee"
    splitting: str = "strang"            # #9 source/transport splitting: "lie" (1st) | "strang" (2nd)
    substep_cfl_growth: float = 2.0      # #10 if post-solve velocity exceeds this x the predicted
    #                                      CFL budget, the transport sub-steps to stay CFL-stable
    clip_warn_frac: float = 0.05         # #11 warn if clip activations exceed this fraction of cell-steps
    outlet_open_frac: float = 0.10       # outlet through-flow opens above this inlet-rate fraction
    #                                      (was a hard-coded 0.10 magic number)
    subgrid_slug: bool = False           # #5 deepened (OPT-IN): blend cell-mean holdup toward the slug
    #                                      BODY holdup in slug/churn flow so the closures see the unit-cell
    #                                      structure. Off by default — the interfacial-area closure already
    #                                      carries a slug/churn enhancement; enable to study the sensitivity.
    volume_consistent_pressure: float = 0.0  # #B3 two-fluid-mass coupling (OPT-IN, 0..1): relax pressure
    #                                      toward the value at which the CONSERVED liquid (La) and gas (Mg)
    #                                      volumes exactly fill the bore (rho_g = Mg/(A-La), inverted via the
    #                                      EOS). 0 = drift-flux pressure only; 0.2-0.4 typical.
    acoustic: float = 0.0                # #B5 water-hammer (OPT-IN, 0..1): weight of a compressible acoustic
    #                                      pressure-wave term (finite sound speed) in the pressure update.
    #  --- Tier-2/3/4 precision options (OPT-IN; defaults reproduce the verified results) ---
    advanced_physics: bool = False       # #9/#10/#11/#12/#13: geometry-resolved interfacial area,
    #                                      regime-dependent friction & Nusselt, Joule-Thomson cooling,
    #                                      and a CaCO3 scaling-tendency screen.
    tvd_energy: bool = False             # #15: 2nd-order TVD reconstruction of the energy advection too
    #                                      (the energy field is smooth + already Heun-corrected; opt-in)
    error_dt: bool = False               # #17: embedded error-controlled adaptive time step
    lhs_uq: bool = False                 # #18: Latin-Hypercube ensemble sampling (variance reduction)
    auto_refine: bool = False            # #16: auto-increase n_cells for steep-terrain thermal resolution
    soil_transient: bool = False         # #12: transient buried-soil thermal inertia
    soil_thermal_mass: float = 2.0e6     # J/m2K lumped-node heat capacity (used when soil_nodes==1)
    soil_far_h: float = 5.0              # W/m2K lumped-node conductance to the far-field seabed
    #  C12 FULL transient radial soil conduction (the axisymmetric heat equation around a buried
    #  pipe): soil_nodes>1 discretises the ground into concentric shells from the wall to a far-field
    #  radius, each with thermal mass, solving the radial conduction transient (a proper 2-D-soil
    #  cooldown profile, not a single lumped node).
    soil_nodes: int = 1                  # number of radial soil shells (1 = lumped node)
    soil_far_radius_m: float = 2.0       # far-field (undisturbed-ground) radius
    soil_conductivity: float = 2.0       # soil thermal conductivity (W/mK, wet soil)
    soil_rhocp: float = 2.5e6            # soil volumetric heat capacity (J/m3K)
    within_step_iters: int = 1           # #6: pressure-momentum iterations to convergence per step (>1)
    slug_tracking: bool = False          # A4: advance a Lagrangian slug-indicator field that seeds at
    #                                      the slug frequency and advects at the translational velocity
    surge_factor: float = 1.6            # slug-catcher surge-volume design multiplier (reporting)
    transportable_mu_rel_max: float = 10.0  # slurry relative-viscosity transportability threshold (reporting)
    twofluid_mass_iters: int = 8         # within-step coupled (momentum<->gas-volume) Picard/Newton
    #                                      iterations for the twofluid_mass engine (volume-consistent p)
    twofluid_mass_newton_w: float = 0.85  # blend weight for the monolithic volume-mass Newton engine
    #                                      (twofluid_mass_newton): 1=exact gas consistency, 0=momentum dP
    #  --- fully-coupled (alpha_l, p, u_m) flux-level Newton engine (twofluid_full_newton) ---
    #  a genuine within-step block Newton on the primitives (alpha_l, p, u_m) per cell: LIQUID-mass
    #  continuity + ELLIPTIC compressible-pressure (mixture-volume) continuity + mixture momentum solved
    #  SIMULTANEOUSLY (block-tridiagonal 3x3 Newton). Converges in ~1-2 iters; reproduces the validated
    #  implicit engine's physical dP and conserves mass exactly (a fully-coupled-vs-segregated verification
    #  result). It does NOT target gas-holdup consistency -> for that niche use twofluid_mass_newton.
    full_newton_iters: int = 12          # max Newton iterations per step (typically converges in ~1-2)
    full_newton_tol: float = 1e-7        # convergence tolerance on max |residual|
    full_newton_relax: float = 1.0       # Newton under-relaxation (1.0 = full step; lowered auto on backtrack)
    strict: bool = False                 # STRICT mode: raise on any hydro fallback or excessive clip
    #                                      activation instead of silently degrading (no masked instability)
    #  --- reduced-order cross-section/3-D reconstruction closures (now CALIBRATABLE, item 9) ---
    cx_velocity_exp: float = 0.142857     # turbulent velocity-profile exponent (1/7 power law)
    cx_deposit_skew: float = 1.6          # azimuthal bottom-of-line deposit weighting skew
    cx_gas_vel_factor: float = 1.25       # gas-region velocity factor in the section reconstruction
    cx_liq_vel_factor: float = 0.85       # liquid-region velocity factor in the section reconstruction
    drift_C0_factor: float = 1.0          # multiplier on the drift-flux distribution parameter C0
    #                                       (default 1.0 = unchanged). >1 raises liquid holdup; used as
    #                                       the strong, physically-grounded knob the SHCT<->CFD two-way
    #                                       coupling tunes to match a CFD-resolved holdup.


@dataclass
class Scenario:
    """Time-dependent inlet boundary conditions (real operations)."""
    kind: str = "steady"                   # steady | rampup | turndown | shutin
    event_time_h: float = 12.0             # when the operational change occurs
    turndown_factor: float = 0.5           # rate multiplier for turndown
    shutin_residual: float = 0.03          # residual flow fraction at shut-in


@dataclass
class Case:
    name: str = "North-Sea-style 20 km tie-back"
    pipeline: Pipeline = field(default_factory=Pipeline)
    fluids: Fluids = field(default_factory=Fluids)
    operating: Operating = field(default_factory=Operating)
    kinetics: Kinetics = field(default_factory=Kinetics)
    numerics: Numerics = field(default_factory=Numerics)
    scenario: Scenario = field(default_factory=Scenario)

