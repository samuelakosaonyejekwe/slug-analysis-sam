#!/usr/bin/env python3
# =============================================================================
#  shct_correlations.py — pure, side-effect-free flow-assurance closures (#21)
# -----------------------------------------------------------------------------
#  Extracted from solver.py so the published correlations are independently
#  importable and unit-testable in isolation (no solver state, no I/O). The
#  dataclass type annotations (f: Fluids, pipe: Pipeline, ...) are lazy strings
#  via `from __future__ import annotations`, so this module needs no import of
#  the data model — it only reads attributes by duck-typing at call time.
# =============================================================================
from __future__ import annotations
import math
import numpy as np

G = 9.81
R_GAS = 8.314
M_MEG, M_H2O = 62.07, 18.015                       # MEG / water molar masses (g/mol)

REGIME_NAMES = {0: "stratified-smooth", 1: "stratified-wavy", 2: "slug",
                3: "annular", 4: "dispersed-bubble", 5: "churn"}


def hydrate_equilibrium_T(P_bar, gas_sg=0.60, salinity_wt=0.0, table=None):
    """Hydrate equilibrium temperature [degC] vs pressure.

    Base natural-gas correlation (fitted to typical sI/sII data: ~3 C @ 30 bar,
    ~10 C @ 70 bar, ~13 C @ 100 bar, ~18 C @ 200 bar), made CASE-SPECIFIC and
    universal without changing the default behaviour:
      * gas_sg : gas specific gravity (air=1). Heavier gas -> hydrates stable at higher
        T. The shift is zero at the reference sg=0.60, so the default curve (and the
        verification anchors) are reproduced EXACTLY.
      * salinity_wt : produced-water salinity (wt% NaCl-eq) depresses Teq (thermodynamic
        inhibition by dissolved salt); zero by default.
      * table : optional [[P_bar, T_C], ...] measured/PVT hydrate curve. When supplied it
        OVERRIDES the correlation (monotonic-interpolated) -> exact, fluid-specific."""
    P = np.maximum(np.asarray(P_bar, float), 1.0)
    if table:
        tab = np.asarray(table, float)
        order = np.argsort(tab[:, 0])
        return np.interp(P, tab[order, 0], tab[order, 1])
    base = 7.7 * np.log(P) - 23.2
    sg_shift = 18.0 * (gas_sg - 0.60)                 # ~+1.8 C per 0.10 sg; 0 at reference
    salt_shift = -0.74 * max(float(salinity_wt), 0.0)  # NaCl depression (~Hammerschmidt-like)
    return base + sg_shift + salt_shift


def gas_Z_factor(P_bar, T_C, f: "Fluids"):
    """Gas compressibility factor Z (#1). Three levels, all back-compatible:
      * default (gas_Z_corr=False, no table): the constant f.gas_Z (identical to before);
      * gas_Z_table given: nearest-point interpolation of a measured/EOS [[P,T,Z],...] surface;
      * gas_Z_corr=True: a Standing-style pseudo-reduced correlation from the gas SG, giving a
        smooth Z(P,T) that captures the dip below 1 at high pressure (real-gas, not ideal)."""
    if getattr(f, "gas_Z_table", None):
        tab = np.asarray(f.gas_Z_table, float)
        P = np.asarray(P_bar, float); T = np.asarray(T_C, float)
        from_pts = tab[:, :2]; zvals = tab[:, 2]
        Pf = np.atleast_1d(P).ravel()[:, None]; Tf = np.atleast_1d(T).ravel()[:, None]
        d2 = (Pf - from_pts[:, 0]) ** 2 + ((Tf - from_pts[:, 1]) * 10.0) ** 2
        z = zvals[np.argmin(d2, axis=1)].reshape(np.shape(P) if np.ndim(P) else (1,))
        return z if np.ndim(P) else float(z[0])
    if not getattr(f, "gas_Z_corr", False):
        return f.gas_Z
    #  Sutton pseudo-criticals from SG, then the standard Beggs-Brill explicit Z(Tr,Pr)
    #  correlation (reproduces Standing-Katz to a few % over 1.2<=Tr<=2.4, Pr<=10).
    sg = max(float(f.gas_sg), 0.55)
    Tpc = 169.2 + 349.5 * sg - 74.0 * sg ** 2          # R (Sutton)
    Ppc = 756.8 - 131.0 * sg - 3.6 * sg ** 2           # psia
    Tr = np.maximum(((np.asarray(T_C, float) + 273.15) * 1.8) / Tpc, 1.02)
    Pr = np.clip((np.asarray(P_bar, float) * 14.5038) / Ppc, 0.0, 12.0)
    A = 1.39 * np.sqrt(np.maximum(Tr - 0.92, 1e-6)) - 0.36 * Tr - 0.101
    B = ((0.62 - 0.23 * Tr) * Pr
         + (0.066 / np.maximum(Tr - 0.86, 1e-3) - 0.037) * Pr ** 2
         + 0.32 * Pr ** 6 / (10.0 ** (9.0 * (Tr - 1.0))))
    C = 0.132 - 0.32 * np.log10(Tr)
    Dd = 10.0 ** (0.3106 - 0.49 * Tr + 0.1824 * Tr ** 2)
    Z = A + (1.0 - A) / np.exp(np.minimum(B, 50.0)) + C * Pr ** Dd
    return np.clip(Z, 0.25, 1.15)


def _pvt_lookup(P_bar, T_C, table, col):
    """Interpolate column `col` from a user PVT property surface
    table = [[P_bar, T_C, rho_oil, rho_gas, mu_oil, mu_gas], ...]  (cols 2..5). Uses smooth
    BILINEAR interpolation (#24) when the table lies on a regular (P,T) grid — no property
    discontinuities — and falls back to nearest-point for a scattered cloud."""
    tab = np.asarray(table, float)
    Ps = np.unique(tab[:, 0]); Ts = np.unique(tab[:, 1])
    P = np.asarray(P_bar, float); T = np.asarray(T_C, float)
    if Ps.size * Ts.size == tab.shape[0] and Ps.size >= 2 and Ts.size >= 2:
        #  regular grid -> reshape to (nP, nT) and bilinear-interpolate
        order = np.lexsort((tab[:, 1], tab[:, 0]))
        grid = tab[order, col].reshape(Ps.size, Ts.size)
        Pi = np.clip(np.interp(P, Ps, np.arange(Ps.size)), 0, Ps.size - 1.0001)
        Ti = np.clip(np.interp(T, Ts, np.arange(Ts.size)), 0, Ts.size - 1.0001)
        i0 = np.floor(Pi).astype(int); j0 = np.floor(Ti).astype(int)
        fp = Pi - i0; ft = Ti - j0
        v = ((1 - fp) * (1 - ft) * grid[i0, j0] + fp * (1 - ft) * grid[i0 + 1, j0]
             + (1 - fp) * ft * grid[i0, j0 + 1] + fp * ft * grid[i0 + 1, j0 + 1])
        return v if np.ndim(P) else float(np.atleast_1d(v)[0])
    Pf = np.atleast_1d(P).ravel()[:, None]; Tf = np.atleast_1d(T).ravel()[:, None]
    d2 = (Pf - tab[:, 0]) ** 2 + ((Tf - tab[:, 1]) * 10.0) ** 2
    v = tab[np.argmin(d2, axis=1), col]
    return v.reshape(np.shape(P)) if np.ndim(P) else float(v[0])


def gas_density(P_bar, T_C, f: "Fluids"):
    if getattr(f, "pvt_table", None):                     # #1: user PVT/EOS property surface
        return _pvt_lookup(P_bar, T_C, f.pvt_table, 3)
    Z = gas_Z_factor(P_bar, T_C, f)
    return (P_bar * 1e5) * f.gas_MW / (Z * R_GAS * (T_C + 273.15))


def gas_viscosity(rho_g, T_C, f: "Fluids"):
    """Gas viscosity (Pa·s) (#1 depth). Default: the constant f.mu_gas. With f.gas_visc_corr,
    the Lee–Gonzalez–Eakin correlation μg(ρg, T, MW) — viscosity rises with pressure (density)
    and temperature, which the constant value misses; matters for friction at high P."""
    if not getattr(f, "gas_visc_corr", False):
        return f.mu_gas
    M = f.gas_MW * 1000.0                                  # kg/kmol -> g/mol units for Lee
    Tr = (np.asarray(T_C, float) + 273.15) * 1.8           # Rankine
    rho_gcc = np.asarray(rho_g, float) / 1000.0            # kg/m3 -> g/cc
    K = (9.4 + 0.02 * M) * Tr ** 1.5 / (209.0 + 19.0 * M + Tr)
    X = 3.5 + 986.0 / Tr + 0.01 * M
    Y = 2.4 - 0.2 * X
    mu_cp = 1.0e-4 * K * np.exp(X * rho_gcc ** Y)          # centipoise
    return np.clip(mu_cp * 1.0e-3, 1e-6, 1e-3)             # cP -> Pa·s


def oil_density(P_bar, T_C, f: "Fluids"):
    """Oil density (kg/m3) (#1 depth). Default constant f.rho_oil; with f.oil_pvt_corr a simple
    thermal-expansion + isothermal-compressibility correction about the reference state, so the
    liquid density varies with P and T (a black-oil-style surface; replace by a user table/EOS)."""
    if getattr(f, "pvt_table", None):                     # #1: user PVT/EOS property surface
        return _pvt_lookup(P_bar, T_C, f.pvt_table, 2)
    if not getattr(f, "oil_pvt_corr", False):
        return f.rho_oil
    beta = 7.0e-4      # 1/K thermal expansion
    cP = 1.0e-4        # 1/bar isothermal compressibility (oil)
    return f.rho_oil * (1.0 + cP * (np.asarray(P_bar, float) - 1.0)
                        - beta * (np.asarray(T_C, float) - 15.0))


#  ------- TVD flux limiters (#8) -------------------------------------------
def _limiter(r, kind="minmod"):
    """Slope-ratio flux limiter phi(r) for 2nd-order TVD advection. r is the ratio of
    consecutive gradients. 'upwind' -> 0 (recovers 1st-order upwind)."""
    if kind == "upwind":
        return np.zeros_like(r)
    r = np.where(np.isfinite(r), r, 0.0)
    if kind == "minmod":
        return np.maximum(0.0, np.minimum(1.0, r))
    if kind == "vanleer":
        return (r + np.abs(r)) / (1.0 + np.abs(r))
    if kind == "superbee":
        return np.maximum.reduce([np.zeros_like(r), np.minimum(2.0 * r, 1.0), np.minimum(r, 2.0)])
    return np.maximum(0.0, np.minimum(1.0, r))                    # default minmod


def tvd_interior_faces(q, vface, kind="minmod"):
    """Limited (2nd-order TVD) face VALUES of a cell quantity q at the nx-1 INTERIOR faces,
    given the interior face velocities vface (sign sets the upwind side). The caller forms
    flux = vface * value and applies the BC faces, so the scheme stays exactly flux-conservative
    (the limiter only sharpens the reconstruction). kind='upwind' reproduces 1st order exactly."""
    qm = q[:-1]; qp = q[1:]
    dq = qp - qm
    dq_safe = np.where(np.abs(dq) < 1e-30, 1e-30, dq)
    r_pos = np.zeros_like(dq); r_neg = np.zeros_like(dq)
    r_pos[1:] = (qm[1:] - qm[:-1]) / dq_safe[1:]                  # vel>0: upstream gradient / local
    r_neg[:-1] = (qp[1:] - qp[:-1]) / dq_safe[:-1]               # vel<0
    qf_pos = qm + 0.5 * _limiter(r_pos, kind) * dq
    qf_neg = qp - 0.5 * _limiter(r_neg, kind) * dq
    return np.where(vface >= 0.0, qf_pos, qf_neg)


def mixture_sound_speed(alpha_l, rho_l, rho_g, p_bar, gamma=1.3, c_liq=1200.0):
    """Two-phase (Wood's-equation) mixture sound speed [m/s] (#B5). The mixture compressibility
    is the holdup-weighted sum of phase compressibilities; even a little gas drops the sound
    speed to O(10-100 m/s) — the key quantity for water-hammer / fast-transient screening.
       1/(rho_m c^2) = alpha_l/(rho_l c_l^2) + alpha_g/(gamma p)
    (gas compressibility ~ gamma*p; liquid c_l ~ 1200 m/s)."""
    ag = np.clip(1.0 - alpha_l, 1e-4, 1.0)
    al = np.clip(alpha_l, 0.0, 1.0)
    rho_m = al * rho_l + ag * rho_g
    p = np.maximum(p_bar, 1.0) * 1e5
    comp = al / np.maximum(rho_l * c_liq ** 2, 1.0) + ag / np.maximum(gamma * p, 1.0)
    return np.sqrt(1.0 / np.maximum(rho_m * comp, 1e-12))


def haaland_friction(Re, rel_rough):
    Re = np.maximum(Re, 1.0)
    inv = -1.8 * np.log10((rel_rough / 3.7) ** 1.11 + 6.9 / Re)
    turb = (1.0 / inv) ** 2
    return np.where(Re < 2300.0, 64.0 / Re, turb)


def drift_params(theta, D):
    """Bendiksen-type distribution coefficient and drift velocity."""
    C0 = 1.05 + 0.15 * np.sin(np.abs(theta))
    vd = 0.35 * np.sqrt(G * D) * np.sin(theta) + 0.20 * np.sqrt(G * D) * np.cos(theta)
    return C0, vd


def flow_regime_code(Vsg, Vsl, theta):
    cond = [(Vsg > 12.0),
            (Vsl > 2.5) & (Vsg < 1.5),
            (Vsg > 5.0) & (Vsg <= 12.0),
            (Vsg < 0.10) & (Vsl < 0.4) & (np.abs(theta) < 0.09),
            (Vsg < 1.0) & (Vsl < 0.3) & (theta < 0.04)]
    return np.select(cond, [3, 4, 5, 0, 1], default=2)


def slug_frequency(Vsl, Vm, D, theta):
    Vm = np.maximum(Vm, 1e-3)
    base = 0.0226 * ((np.maximum(Vsl, 1e-4) / (G * D)) * (19.75 / Vm + Vm)) ** 1.2
    incl = 0.836 + 2.75 * np.sin(np.abs(theta)) ** 0.25
    return np.maximum(base * incl, 1e-4)


def interfacial_area(alpha_g, alpha_l, D, regime):
    enh = np.where(np.isin(regime, [2, 5]), 7.0, np.where(regime == 3, 10.0, 3.0))
    return enh * alpha_g * alpha_l * 4.0 / D


def interfacial_area_geom(alpha_l, D, regime):
    """GEOMETRY-resolved interfacial area per unit volume [1/m] (#10). Instead of a constant
    enhancement factor, use the actual gas-liquid interface geometry of each regime:
      * stratified: flat interface of width from the wetted-angle (Taitel-Dukler);
      * annular: cylindrical film interface ~ 4*sqrt(alpha_l)/D;
      * slug/churn: dispersed bubbles -> high area ~ 6*alpha_g/d_b with d_b ~ 0.1 D;
      * bubble: dispersed bubbles in liquid."""
    al = np.clip(alpha_l, 1e-3, 0.999); ag = 1.0 - al
    #  stratified wetted half-angle from holdup (Biberg approximation of Taitel-Dukler)
    gamma = np.pi * al + (1.5 * np.pi) ** (1.0 / 3.0) * (1.0 - 2.0 * al + al ** (1.0 / 3.0) - ag ** (1.0 / 3.0))
    a_strat = np.sin(np.clip(gamma, 0.0, np.pi)) / D                   # interface width / area
    a_ann = 4.0 * np.sqrt(al) / D                                      # annular film interface
    db = 0.10 * D                                                     # dispersed bubble diameter
    a_disp = 6.0 * ag / db                                            # slug/churn/bubble dispersed
    return np.where(np.isin(regime, [0, 1]), a_strat,
                    np.where(regime == 3, a_ann, a_disp))


def regime_friction_multiplier(regime):
    """Two-phase wall-friction multiplier by flow regime (#9): slug/churn flow has a much higher
    effective wall shear than smooth stratified flow (intermittent contact, mixing). A coarse,
    regime-dependent multiplier on the mixture friction (vs a single regime-independent value)."""
    #  0 strat-smooth, 1 strat-wavy, 2 slug, 3 annular, 4 disp-bubble, 5 churn
    return np.select([regime == 0, regime == 1, regime == 2, regime == 3, regime == 4, regime == 5],
                     [0.9, 1.0, 1.6, 1.2, 1.1, 1.7], default=1.0)


def regime_nusselt(regime, Re, Pr):
    """Inner-wall Nusselt number by flow regime (#12): turbulent Dittus-Boelter baseline with a
    regime enhancement (slug/churn mixing boosts the inner heat-transfer coefficient)."""
    Nu0 = 0.023 * np.maximum(Re, 1.0) ** 0.8 * np.maximum(Pr, 0.1) ** 0.4
    enh = np.select([np.isin(regime, [2, 5]), regime == 3], [1.5, 1.3], default=1.0)
    return np.maximum(Nu0 * enh, 3.66)                                # laminar floor


def joule_thomson_dTdP(P_bar, T_C, f: "Fluids"):
    """Joule-Thomson coefficient mu_JT = (dT/dP) [K/Pa] for the gas (#11): real-gas cooling on
    expansion. From the EOS-style relation mu_JT ~ (RT^2/Pcp)(dZ/dT)/Z; positive for gas at
    pipeline conditions (cools as it depressurises). Reduced form via the Z-factor slope."""
    T = np.asarray(T_C, float) + 273.15
    dT = 1.0
    Z1 = gas_Z_factor(P_bar, T_C, f); Z2 = gas_Z_factor(P_bar, np.asarray(T_C, float) + dT, f)
    dZdT = (np.asarray(Z2, float) - np.asarray(Z1, float)) / dT
    cp = getattr(f, "cp_gas", 2300.0)
    rho_mol = (np.asarray(P_bar, float) * 1e5) / (np.maximum(Z1, 0.1) * R_GAS * T)   # mol/m3 (P/ZRT)
    return (T / (cp * f.gas_MW * np.maximum(rho_mol, 1e-6) * np.maximum(Z1, 0.1))) * dZdT


def scaling_tendency_index(T_C, salinity_wt, P_bar):
    """Calcium-carbonate scaling-tendency index (#13, reduced Stiff-Davis style): a screening
    saturation index that rises with temperature and salinity and falls with pressure (CO2).
    >0 indicates a scaling tendency (CaCO3 deposition risk) — a flow-assurance solids screen
    alongside hydrate. Reduced model: replace with full PHREEQC-style speciation for rigour."""
    return (0.02 * (np.asarray(T_C, float) - 20.0)
            + 0.03 * np.maximum(salinity_wt, 0.0)
            - 0.004 * (np.asarray(P_bar, float) - 50.0))


def slug_length(Vm, D, fslug):
    """Mean slug unit length L_u [m] (#5): a reduced-order, sub-grid descriptor of the
    individual metre-scale slugs a dx~200 m grid cannot resolve. L_u = V_t / f_slug with the
    slug translational velocity V_t ~ 1.2*Vm + drift (Dukler-Hubbard / Gregory-Scott scaling)."""
    Vt = 1.2 * np.maximum(Vm, 1e-3) + 0.35 * np.sqrt(G * D)
    return np.clip(Vt / np.maximum(fslug, 1e-4), D, 5000.0)


def droplet_entrainment_frac(Vsg, rho_g, sigma, D):
    """Liquid fraction entrained as droplets in the gas core (#6), Ishii-Mishima-style:
    rises with gas inertia (Weber) toward a plateau; ~0 in stratified/slug, significant only
    in high-shear annular/churn flow."""
    We = rho_g * np.maximum(Vsg, 0.0) ** 2 * D / max(float(sigma), 1e-4)
    return np.clip(1.0 - np.exp(-We / 1.0e5), 0.0, 0.95)


#  MEG inhibition: Nielsen-Bucklin (mole-fraction) model, valid to high wt% where the
#  classical Hammerschmidt linear form breaks down. ΔT = -72 ln(1 - x_MEG_in_water).
def _meg_wt_to_molefrac(W_wt):
    W = np.clip(W_wt, 0.0, 80.0)
    n_meg = W / M_MEG
    n_w = np.maximum(100.0 - W, 0.0) / M_H2O
    return n_meg / np.maximum(n_meg + n_w, 1e-12)


def meg_suppression(W_wt):
    """Hydrate-equilibrium temperature depression (degC) for a local inhibitor concentration
    W (wt% in the AQUEOUS phase) — Nielsen-Bucklin (valid to ~50 wt%)."""
    x = _meg_wt_to_molefrac(W_wt)
    return -72.0 * np.log(np.maximum(1.0 - x, 1e-6))


def hammerschmidt_meg(dT_required_C, water_massrate_kgps):
    """Required MEG wt% in the aqueous phase, mass rate and volume rate to achieve a hydrate-
    curve depression of dT_required_C. Inverts the Nielsen-Bucklin model (kept under the
    historical name for API compatibility)."""
    dT = max(float(dT_required_C), 0.0)
    x = 1.0 - math.exp(-dT / 72.0)                       # required MEG mole fraction in water
    x = min(x, 0.60)
    W = 100.0 * x * M_MEG / ((1.0 - x) * M_H2O + x * M_MEG)   # back to wt%
    W = min(W, 60.0)
    meg_kgps = water_massrate_kgps * (W / max(100.0 - W, 1e-3))
    return W, meg_kgps, meg_kgps / 1113.0 * 1000.0 * 3600.0


def effective_U_and_mass(pipe: "Pipeline", op: "Operating", cp_fluid, rho_fluid_eff):
    """Effective overall heat-transfer coefficient (W/m2K, referred to inner area) and lumped
    thermal mass per metre (J/mK) of fluid + wall. Uses pipe.wall_layers if provided, else
    Operating.U_wall with a default steel+insulation wall mass."""
    D = pipe.diameter_m
    A_in = math.pi * D ** 2 / 4.0
    fluid_mass = rho_fluid_eff * cp_fluid * A_in           # J/mK (per unit length)
    if pipe.wall_layers:
        Rinv = 1.0 / pipe.h_inner
        wall_mass = 0.0
        for layer in pipe.wall_layers:
            th, k = float(layer[0]), float(layer[1])
            rhoCp = float(layer[2]) if len(layer) > 2 else 3.5e6   # default steel-ish
            Rinv += th / max(k, 1e-6)
            wall_mass += rhoCp * (math.pi * (D + th) * th)        # ring mass approx
        Rinv += 1.0 / pipe.h_outer
        U = 1.0 / Rinv
    else:
        U = op.U_wall
        wall_mass = (3.5e6 * math.pi * D * 0.025) + (1.2e5 * math.pi * D * 0.050)
    return U, fluid_mass + wall_mass
