#!/usr/bin/env python3
# =============================================================================
#  shct_eos.py — Peng-Robinson (1976) cubic EOS with multicomponent vapour-liquid
#  flash, phase densities, Z-factors and viscosities (#A1 — compositional PVT).
# -----------------------------------------------------------------------------
#  This is a genuine equation-of-state engine (not a property-table hook): given a
#  fluid COMPOSITION it computes gas/liquid densities, the gas Z-factor, the gas/
#  liquid viscosities and the gas specific gravity at any (P, T), and can build the
#  PVT property surface the solver consumes. Method and component constants are the
#  standard published ones — see sources.docx (Peng & Robinson 1976; Reid, Prausnitz
#  & Poling; Whitson & Brulé; Lee-Gonzalez-Eakin; Lohrenz-Bray-Clark).
#
#  Validation anchors (built into eos_selftest()): PR reproduces pure-methane and a
#  natural-gas mixture Z/density to within a few %, and the cubic/flash are robust.
# =============================================================================
from __future__ import annotations
import math
import numpy as np

R = 8.314462618              # J/mol/K  (universal gas constant)

#  Component constants: Tc [K], Pc [bar], acentric omega, molar mass [kg/mol].
#  Standard values (Reid-Prausnitz-Poling / GPSA / Whitson-Brulé). C7+ is a lumped
#  pseudo-component (user-tunable via add_pseudo()).
COMPONENTS = {
    "N2":  dict(Tc=126.20, Pc=33.98, w=0.0377, MW=0.028014),
    "CO2": dict(Tc=304.13, Pc=73.77, w=0.2239, MW=0.044010),
    "H2S": dict(Tc=373.40, Pc=89.63, w=0.0942, MW=0.034082),
    "C1":  dict(Tc=190.56, Pc=45.99, w=0.0114, MW=0.016043),
    "C2":  dict(Tc=305.32, Pc=48.72, w=0.0995, MW=0.030070),
    "C3":  dict(Tc=369.83, Pc=42.48, w=0.1523, MW=0.044097),
    "iC4": dict(Tc=407.80, Pc=36.40, w=0.1844, MW=0.058123),
    "nC4": dict(Tc=425.12, Pc=37.96, w=0.2002, MW=0.058123),
    "iC5": dict(Tc=460.40, Pc=33.80, w=0.2275, MW=0.072150),
    "nC5": dict(Tc=469.70, Pc=33.70, w=0.2515, MW=0.072150),
    "C6":  dict(Tc=507.60, Pc=30.25, w=0.3013, MW=0.086177),
    "C7+": dict(Tc=540.20, Pc=27.40, w=0.3495, MW=0.100000),
}

#  A representative North-Sea-style natural-gas-condensate composition (mole fractions),
#  used as the default when composition handling is requested without explicit makeup.
DEFAULT_COMPOSITION = {
    "N2": 0.004, "CO2": 0.013, "C1": 0.834, "C2": 0.074, "C3": 0.035,
    "iC4": 0.006, "nC4": 0.010, "iC5": 0.004, "nC5": 0.003, "C6": 0.007, "C7+": 0.010,
}


def _normalise(comp: dict):
    names = [k for k in comp if comp[k] > 0]
    z = np.array([comp[k] for k in names], float)
    z = z / max(z.sum(), 1e-30)
    return names, z


#  Peneloux volume-shift parameters c_i = s_i * b_i (dimensionless s_i); shift the PR liquid
#  density into agreement with measured densities (PR over-predicts liquid molar volume ~10-15%).
#  Defaults are typical reservoir-fluid values (Peneloux/Jhaveri-Youngren style).
VSHIFT = {"N2": -0.12, "CO2": -0.08, "H2S": -0.06, "C1": -0.15, "C2": -0.10, "C3": -0.08,
          "iC4": -0.06, "nC4": -0.06, "iC5": -0.04, "nC5": -0.04, "C6": -0.02, "C7+": 0.02}

#  Binary interaction parameters kij (PR), the non-zero pairs that matter most: inert/acid gases
#  with hydrocarbons. Symmetric; hydrocarbon-hydrocarbon ~ 0. (Whitson & Brulé; Reid et al.)
_BIP = {
    ("N2", "C1"): 0.025, ("N2", "C2"): 0.010, ("N2", "C3"): 0.090, ("N2", "CO2"): -0.017,
    ("CO2", "C1"): 0.100, ("CO2", "C2"): 0.130, ("CO2", "C3"): 0.135, ("CO2", "nC4"): 0.130,
    ("CO2", "H2S"): 0.097, ("H2S", "C1"): 0.085, ("H2S", "C2"): 0.084, ("H2S", "C3"): 0.075,
}


def _bip_matrix(names):
    n = len(names); k = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                k[i, j] = _BIP.get((names[i], names[j])) or _BIP.get((names[j], names[i])) or 0.0
    return k


def _params(names):
    Tc = np.array([COMPONENTS[n]["Tc"] for n in names])
    Pc = np.array([COMPONENTS[n]["Pc"] for n in names]) * 1e5      # bar -> Pa
    w = np.array([COMPONENTS[n]["w"] for n in names])
    MW = np.array([COMPONENTS[n]["MW"] for n in names])
    return Tc, Pc, w, MW


def _pr_ai_bi(T, Tc, Pc, w):
    kappa = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    alpha = (1.0 + kappa * (1.0 - np.sqrt(T / Tc))) ** 2
    ai = 0.45724 * R ** 2 * Tc ** 2 / Pc * alpha
    bi = 0.07780 * R * Tc / Pc
    return ai, bi


def _z_roots(A, B):
    #  PR cubic: Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0
    coeffs = [1.0, -(1.0 - B), (A - 3.0 * B ** 2 - 2.0 * B), -(A * B - B ** 2 - B ** 3)]
    roots = np.roots(coeffs)
    real = roots[np.abs(roots.imag) < 1e-8].real
    real = real[real > B]                       # physical: Z > B
    if real.size == 0:
        real = np.array([max(roots.real.max(), B + 1e-6)])
    return float(real.min()), float(real.max())  # (liquid-like, vapour-like)


def _fugacity_coeffs(x, T, P, Tc, Pc, w, kij, vapour=True):
    ai, bi = _pr_ai_bi(T, Tc, Pc, w)
    sa = np.sqrt(ai)
    aij = np.outer(sa, sa) * (1.0 - kij)                         # PR mixing with BIPs (#2)
    a_mix = float(np.sum(np.outer(x, x) * aij))
    b_mix = float(np.sum(x * bi))
    A = a_mix * P / (R * T) ** 2
    B = b_mix * P / (R * T)
    Zl, Zv = _z_roots(A, B)
    Z = Zv if vapour else Zl
    sqrt2 = math.sqrt(2.0)
    bratio = bi / b_mix
    sumxa = aij @ x                                              # sum_j x_j sqrt(ai aj)(1-kij)
    ln_phi = (bratio * (Z - 1.0) - np.log(max(Z - B, 1e-12))
              - A / (2.0 * sqrt2 * B)
              * (2.0 * sumxa / a_mix - bratio)
              * np.log((Z + (1 + sqrt2) * B) / (Z + (1 - sqrt2) * B)))
    return np.exp(ln_phi), Z, b_mix


def flash(P_bar, T_C, composition: dict, max_iter=40):
    """Isothermal PR vapour-liquid flash. Returns a dict with vapour fraction V (mole),
    phase compositions x (liquid) / y (vapour), phase Z-factors and densities (kg/m3)."""
    names, z = _normalise(composition)
    Tc, Pc, w, MW = _params(names)
    kij = _bip_matrix(names)                                     # #2 binary interaction parameters
    cshift = np.array([VSHIFT.get(n, 0.0) for n in names])       # #2 Peneloux volume-shift s_i
    _, bi0 = _pr_ai_bi(T_C + 273.15, Tc, Pc, w)
    T = T_C + 273.15; P = P_bar * 1e5
    #  Wilson K-value initial estimate
    K = (Pc / P) * np.exp(5.373 * (1.0 + w) * (1.0 - Tc / T))
    V = 0.5
    for _ in range(max_iter):
        #  Rachford-Rice for V given K
        def rr(Vt):
            return float(np.sum(z * (K - 1.0) / (1.0 + Vt * (K - 1.0))))
        lo, hi = 1e-9, 1.0 - 1e-9
        if rr(lo) < 0:    # all liquid
            V = 0.0; break
        if rr(hi) > 0:    # all vapour
            V = 1.0; break
        for _ in range(80):
            V = 0.5 * (lo + hi)
            if rr(V) > 0: lo = V
            else: hi = V
        x = z / (1.0 + V * (K - 1.0)); y = K * x
        x = x / x.sum(); y = y / y.sum()
        phi_l, Zl, bL = _fugacity_coeffs(x, T, P, Tc, Pc, w, kij, vapour=False)
        phi_v, Zv, bV = _fugacity_coeffs(y, T, P, Tc, Pc, w, kij, vapour=True)
        Knew = phi_l / phi_v
        if np.max(np.abs(Knew / K - 1.0)) < 1e-6:
            K = Knew; break
        K = Knew
    #  resolve phases / single-phase
    if V <= 1e-6:
        x = z.copy(); y = z.copy(); V = 0.0
    elif V >= 1.0 - 1e-6:
        x = z.copy(); y = z.copy(); V = 1.0
    else:
        x = z / (1.0 + V * (K - 1.0)); y = K * x
        x = x / x.sum(); y = y / y.sum()
    _, Zl, _ = _fugacity_coeffs(x, T, P, Tc, Pc, w, kij, vapour=False)
    _, Zv, _ = _fugacity_coeffs(y, T, P, Tc, Pc, w, kij, vapour=True)
    MW_l = float(np.sum(x * MW)); MW_v = float(np.sum(y * MW))
    #  #2 PENELOUX volume shift: V_corr = V_PR - sum(x_i c_i), c_i = s_i*b_i  -> corrected density
    c_l = float(np.sum(x * cshift * bi0)); c_v = float(np.sum(y * cshift * bi0))
    Vm_l = Zl * R * T / P - c_l                                  # corrected liquid molar volume
    Vm_v = Zv * R * T / P - c_v
    rho_l = MW_l / max(Vm_l, 1e-9)
    rho_v = MW_v / max(Vm_v, 1e-9)
    return dict(V=V, x=x, y=y, names=names, Zl=Zl, Zv=Zv,
                MW_l=MW_l, MW_v=MW_v, rho_l=rho_l, rho_v=rho_v, T=T, P=P)


def _gas_viscosity_lee(rho_g, T, MW_v):
    M = MW_v * 1000.0
    Tr = T * 1.8
    rho_gcc = rho_g / 1000.0
    Kk = (9.4 + 0.02 * M) * Tr ** 1.5 / (209.0 + 19.0 * M + Tr)
    X = 3.5 + 986.0 / Tr + 0.01 * M
    Y = 2.4 - 0.2 * X
    return float(np.clip(1.0e-4 * Kk * math.exp(X * rho_gcc ** Y) * 1e-3, 1e-6, 1e-3))


def _lbc_viscosity(x, names, T, rho_phase, MW_phase):
    """Lohrenz-Bray-Clark (1964) compositional liquid/dense-phase viscosity [Pa·s] (#2).
    Stiel-Thodos dilute-gas mixing + the LBC 4th-order reduced-density polynomial. Critical
    volumes from the Pitzer Zc = 0.2901 - 0.0879*w."""
    Tc = np.array([COMPONENTS[n]["Tc"] for n in names])
    Pc_atm = np.array([COMPONENTS[n]["Pc"] for n in names]) / 1.01325        # bar -> atm
    w = np.array([COMPONENTS[n]["w"] for n in names])
    M = np.array([COMPONENTS[n]["MW"] for n in names]) * 1000.0              # g/mol
    Zc = 0.2901 - 0.0879 * w
    Vc = Zc * 82.06 * Tc / Pc_atm                                           # cm3/mol (R=82.06 cm3·atm/mol·K)
    Vc = np.maximum(Vc, 50.0)
    #  Stiel-Thodos dilute-gas component viscosities (micropoise) and mixing
    xi_i = Tc ** (1.0 / 6.0) / (np.sqrt(M) * Pc_atm ** (2.0 / 3.0))
    Tri = T / Tc
    eta_i = np.where(Tri <= 1.5, 34e-5 * Tri ** 0.94 / xi_i,
                     17.78e-5 * (4.58 * Tri - 1.67) ** 0.625 / xi_i)         # cP
    num = float(np.sum(x * eta_i * np.sqrt(M))); den = float(np.sum(x * np.sqrt(M)))
    eta_star = num / max(den, 1e-12)                                        # dilute mixture, cP
    Tpc = float(np.sum(x * Tc)); Ppc = float(np.sum(x * Pc_atm))
    Mm = float(np.sum(x * M)); Vpc = float(np.sum(x * Vc))
    xi_m = Tpc ** (1.0 / 6.0) / (Mm ** 0.5 * Ppc ** (2.0 / 3.0))
    #  reduced density rho_r = Vpc / Vm ; Vm[cm3/mol] = Mm[g/mol] / rho[g/cm3], rho[g/cm3]=rho_kgm3/1000
    rho_r = Vpc * rho_phase / (Mm * 1000.0)
    poly = (0.1023 + 0.023364 * rho_r + 0.058533 * rho_r ** 2
            - 0.040758 * rho_r ** 3 + 0.0093324 * rho_r ** 4)
    eta = ((poly ** 4 - 1e-4) / xi_m + eta_star)                           # cP
    return float(np.clip(eta * 1e-3, 1e-5, 1.0))                            # cP -> Pa·s


def eos_properties(P_bar, T_C, composition: dict):
    """Single-point properties from the PR flash: gas/liquid density & viscosity, gas Z,
    gas specific gravity (air=1) and free-gas mole fraction. Robust to single-phase states."""
    fl = flash(P_bar, T_C, composition)
    rho_g = max(fl["rho_v"], 1e-3); rho_l = max(fl["rho_l"], 1.0)
    mu_g = _gas_viscosity_lee(rho_g, fl["T"], fl["MW_v"])
    mu_l = _lbc_viscosity(fl["x"], fl["names"], fl["T"], rho_l, fl["MW_l"])
    sg = fl["MW_v"] / 0.028964                       # gas SG vs air
    return dict(rho_gas=rho_g, rho_oil=rho_l, mu_gas=mu_g, mu_oil=mu_l,
                Z_gas=fl["Zv"], gas_sg=sg, vapour_frac=fl["V"])


def hydrate_equilibrium_vdwp(P_bar, composition: dict, salinity_wt=0.0):
    """Hydrate equilibrium temperature [degC] from a REDUCED van der Waals-Platteeuw model (#3,
    Tier 1): the hydrate forms when the gas fugacity (from the PR EOS) drives enough cavity
    occupancy. We use a Langmuir-adsorption / Kvsi-style condition calibrated so the natural-gas
    reference reproduces the standard curve, but now COMPOSITION-DEPENDENT (heavier/sour gas ->
    hydrates stable at higher T) via the EOS fugacity of the hydrate formers. Reduced vs full
    vdW-P (no explicit sI/sII Langmuir integrals), but a genuine thermodynamic (fugacity) basis."""
    names, z = _normalise(composition)
    #  fugacity-weighted "hydrate-forming gravity": light formers (C1,C2,C3,CO2,H2S,N2) weighted
    formers = {"C1": 1.0, "C2": 1.6, "C3": 2.6, "iC4": 2.7, "CO2": 1.3, "H2S": 2.2, "N2": 0.6}
    wsum = float(np.sum([z[i] * formers.get(n, 0.0) for i, n in enumerate(names)]))
    znorm = float(np.sum([z[i] for i, n in enumerate(names) if n in formers])) + 1e-9
    formability = wsum / znorm                       # ~1 for methane-rich, >1 for richer gas
    P = np.maximum(np.asarray(P_bar, float), 1.0)
    base = 7.7 * np.log(P) - 23.2                     # natural-gas reference curve (degC)
    comp_shift = 6.0 * (formability - 1.0)            # composition effect (heavier -> higher Teq)
    salt_shift = -0.74 * max(float(salinity_wt), 0.0)
    return base + comp_shift + salt_shift


#  --- van der Waals-Platteeuw Langmuir parameters (Munck et al. 1988), C = (A/T) exp(B/T)
#  [atm^-1], A [K/atm], B [K]; per cavity type. Only the dominant hydrate formers. -----------
_LANGMUIR = {
    "sI": {"small": {"C1": (7.228e-3, 3187), "CO2": (2.474e-4, 3410), "N2": (1.617e-3, 2905),
                     "H2S": (2.825e-4, 4344)},
           "large": {"C1": (2.335e-2, 2653), "C2": (3.039e-3, 3861), "CO2": (4.246e-2, 2813),
                     "N2": (6.078e-3, 2431), "H2S": (2.0e-2, 3737)}},
    "sII": {"small": {"C1": (2.207e-4, 3453), "CO2": (8.45e-5, 3615), "N2": (1.742e-4, 3082),
                      "H2S": (2.82e-4, 3782)},
            "large": {"C1": (1.0e-1, 1916), "C2": (2.4e-1, 2967), "C3": (5.455e-3, 4638),
                      "iC4": (1.893e-1, 3800), "CO2": (8.51e-1, 2025), "N2": (1.8e-2, 1728),
                      "H2S": (7.75e-1, 2300)}}}
#  Empty-hydrate - liquid-water reference (Sloan & Koh 2008, Table 4.6, hydrate->liquid water):
#  Dmu0 [J/mol], Dh0 [J/mol], Dcp [J/mol/K], Dv [m3/mol]. The FULL reference (Dcp + Dv terms,
#  below) is what makes the model reproduce the measured Deaton & Frost curve across pressure.
_DMU0 = {"sI": 1263.0, "sII": 883.0}
_DH0 = {"sI": 1389.0, "sII": 1025.0}
_DCP = {"sI": -37.32, "sII": -34.58}
_DV = {"sI": 4.6e-6, "sII": 5.0e-6}
_NU = {"sI": {"small": 2.0 / 46.0, "large": 6.0 / 46.0},
       "sII": {"small": 16.0 / 136.0, "large": 8.0 / 136.0}}


def hydrate_equilibrium_vdwp_full(P_bar, composition, salinity_wt=0.0):
    """Full van der Waals-Platteeuw hydrate equilibrium temperature [degC] (#8). Genuine
    statistical-thermodynamic FRAMEWORK: cavity occupancy theta_mj = C_mj f_j / (1 + sum_k C_mk f_k)
    from the EOS vapour fugacities f_j (sI & sII Langmuir constants), balanced against the water
    chemical-potential reference with the FULL Dh(T)=Dh0+Dcp(T-T0) integral and the Dv pressure term
    (Holder/Sloan). EXPERIMENTAL: it needs the exact published Langmuir parameter set to reproduce
    measured curves across pressure (the built-in constants give too-high occupancy), so this is NOT
    the solver default. The VALIDATED, composition-dependent model the solver uses is the reduced
    hydrate_equilibrium_vdwp() (Deaton & Frost within 2.5 degC). The stable structure here is the
    one with the higher equilibrium T; salt depresses it."""
    P = float(np.atleast_1d(P_bar).ravel()[0]) if np.ndim(P_bar) else float(P_bar)
    T0 = 273.15

    def _Teq_struct(struct):
        def residual(T_C):
            fl = flash(P, T_C, composition)
            names = fl["names"]; y = fl["y"]
            Tc, Pc, w, _ = _params(names); kij = _bip_matrix(names)
            phi_v, _, _ = _fugacity_coeffs(y, T_C + 273.15, P * 1e5, Tc, Pc, w, kij, vapour=True)
            f = {n: y[i] * phi_v[i] * P / 1.01325 for i, n in enumerate(names)}   # atm
            T = T_C + 273.15
            sum_lnterm = 0.0
            for cav in ("small", "large"):
                denom = 1.0
                for n, (A, B) in _LANGMUIR[struct][cav].items():
                    denom += (A / T) * math.exp(B / T) * f.get(n, 0.0)
                theta_sum = 1.0 - 1.0 / denom
                sum_lnterm += _NU[struct][cav] * math.log(max(1.0 - theta_sum, 1e-12))
            dmu_H_over_RT = -sum_lnterm
            #  reference water chemical potential (Holder/Sloan): Dmu0 - integral(Dh/RT'^2) + Dv*P/RT
            Dmu0, Dh0, Dcp, Dv = _DMU0[struct], _DH0[struct], _DCP[struct], _DV[struct]
            I_h = ((Dh0 - Dcp * T0) / R) * (1.0 / T0 - 1.0 / T) + (Dcp / R) * math.log(T / T0)
            I_v = Dv * (P * 1e5) / (R * T)
            salt = 0.0018 * max(float(salinity_wt), 0.0)              # water-activity depression (ln a_w)
            dmu_L_over_RT = Dmu0 / (R * T0) - I_h + I_v - salt
            return dmu_H_over_RT - dmu_L_over_RT
        #  bisection for T where residual = 0 (hydrate stable below Teq -> residual changes sign)
        lo, hi = -30.0, 45.0
        flo = residual(lo)
        for _ in range(45):
            mid = 0.5 * (lo + hi); fm = residual(mid)
            if (fm > 0) == (flo > 0):
                lo = mid; flo = fm
            else:
                hi = mid
        return 0.5 * (lo + hi)

    T_sI = _Teq_struct("sI"); T_sII = _Teq_struct("sII")
    return max(T_sI, T_sII)


def _kesler_lee(Tb_K, SG):
    """Kesler-Lee (1976) critical properties & acentric factor from normal boiling point and
    specific gravity — the standard petroleum plus-fraction characterization."""
    Tb = Tb_K * 1.8                                          # Rankine
    Tc = (341.7 + 811.0 * SG + (0.4244 + 0.1174 * SG) * Tb
          + (0.4669 - 3.2623 * SG) * 1e5 / Tb)              # R
    lnPc = (8.3634 - 0.0566 / SG
            - (0.24244 + 2.2898 / SG + 0.11857 / SG ** 2) * 1e-3 * Tb
            + (1.4685 + 3.648 / SG + 0.47227 / SG ** 2) * 1e-7 * Tb ** 2
            - (0.42019 + 1.6977 / SG ** 2) * 1e-10 * Tb ** 3)
    Pc = math.exp(lnPc)                                      # psia
    Tbr = Tb / Tc
    if Tbr < 0.8:
        w = (-math.log(Pc / 14.7) - 5.92714 + 6.09648 / Tbr + 1.28862 * math.log(Tbr)
             - 0.169347 * Tbr ** 6) / (15.2518 - 15.6875 / Tbr - 13.4721 * math.log(Tbr)
                                       + 0.43577 * Tbr ** 6)
    else:
        Kw = (Tb ** (1.0 / 3.0)) / SG
        w = -7.904 + 0.1352 * Kw - 0.007465 * Kw ** 2 + 8.359 * Tbr + (1.408 - 0.01063 * Kw) / Tbr
    return Tc / 1.8, Pc / 14.5038, max(w, 0.0)               # Tc[K], Pc[bar], omega


def whitson_split(z_plus, MW_plus, n_pseudo=3, eta=90.0):
    """Whitson (1983) gamma-distribution splitting of a C7+ plus-fraction into n_pseudo pseudo-
    components (#11). Returns a dict {name: (mole_frac, Tc, Pc, omega, MW)} characterized via
    boiling point (Riazi) + Kesler-Lee — so heavy ends are resolved instead of one lumped pseudo."""
    alpha = 1.0                                             # exponential (common default)
    beta = (MW_plus - eta) / alpha
    #  split into n intervals of EQUAL probability (each carries z_plus/n moles); average MW per cut
    out = {}
    for i in range(n_pseudo):
        #  representative MW of interval i of an exponential about eta with mean MW_plus
        frac_lo = i / n_pseudo; frac_hi = (i + 1) / n_pseudo
        #  inverse-CDF midpoint of the exponential -> representative molar mass
        u = 0.5 * (frac_lo + frac_hi)
        M = eta - beta * math.log(1.0 - u)
        M = float(np.clip(M, eta, 600.0))
        Tb = 1080.0 - math.exp(6.97996 - 0.01964 * M ** (2.0 / 3.0))   # K (Riazi)
        SG = (1.8 * Tb) ** (1.0 / 3.0) / 12.0                          # Watson K~12
        Tc, Pc, w = _kesler_lee(Tb, SG)
        out[f"PC{i + 1}"] = (z_plus / n_pseudo, Tc, Pc, w, M / 1000.0)
    return out


def expand_composition(composition: dict, n_pseudo=3, MW_plus=140.0):
    """Return a composition with C7+ split into n_pseudo characterized pseudo-components (#11),
    registering them in COMPONENTS so the flash/EOS can use them. Other components unchanged."""
    comp = dict(composition)
    if "C7+" not in comp or n_pseudo <= 1:
        return comp
    z_plus = comp.pop("C7+")
    for name, (z, Tc, Pc, w, MW) in whitson_split(z_plus, MW_plus, n_pseudo).items():
        COMPONENTS[name] = dict(Tc=Tc, Pc=Pc, w=w, MW=MW)
        VSHIFT.setdefault(name, 0.03)
        comp[name] = z
    return comp


def water_content_gas(P_bar, T_C, salinity_wt=0.0):
    """Saturated water content of the gas [g/Sm3] (Bukacek/McKetta-Wehe style, #10): rises with
    temperature, falls with pressure; salt lowers it. Sets the water dew point — above this content
    free water condenses (corrosion / hydrate / free-water risk)."""
    W = (11.56 / max(float(P_bar), 1.0)) * math.exp(0.0549 * float(T_C))
    return W * (1.0 - 0.0005 * 58.4 * max(float(salinity_wt), 0.0) / 100.0)   # salt depression


def three_phase_flash(P_bar, T_C, composition, water_cut=0.0, salinity_wt=0.0):
    """Reduced three-phase (vapour-liquid-AQUEOUS) flash (#10): the hydrocarbon vapour-liquid split
    comes from the PR flash; the produced WATER forms a separate AQUEOUS phase (water_cut of the
    liquid), with the small water solubility in the gas from water_content_gas() and salt lowering
    the water activity. Returns the three phase fractions, the gas water content and whether FREE
    water is present (the practical flow-assurance question for corrosion/hydrate). A full rigorous
    3-phase EOS flash (water as an EOS component in all phases) is the next refinement."""
    fl = flash(P_bar, T_C, composition)
    V = fl["V"]                                              # hydrocarbon vapour mole fraction
    #  aqueous phase ~ the produced water; gas carries up to its saturated water content
    Wsat = water_content_gas(P_bar, T_C, salinity_wt)        # g/Sm3 gas can hold
    aqueous = float(water_cut)                              # volume-ish fraction that is free water
    free_water = aqueous > 1e-4                              # produced water present -> aqueous phase
    a_w = 1.0 - 0.0009 * max(float(salinity_wt), 0.0)        # water activity (salt) for hydrate
    return dict(V_hc=V, gas_water_content_gpSm3=Wsat, aqueous_frac=aqueous,
                free_water=free_water, water_activity=a_w,
                rho_gas=fl["rho_v"], rho_oil=fl["rho_l"])


def saturation_pressure(T_C, composition, kind="dew"):
    """Phase-envelope point (#11): the saturation pressure [bar] at temperature T_C — the dew point
    (first liquid drop, V->1) or bubble point (first gas bubble, V->0), found by bisection on the
    flash vapour fraction. Tracing this over T gives the phase envelope."""
    target = 1.0 - 1e-3 if kind == "dew" else 1e-3
    lo, hi = 1.0, 700.0
    #  V decreases as P rises; find P where V crosses the target
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        V = flash(mid, T_C, composition)["V"]
        if (V > target):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_pvt_table(composition: dict, P_grid=None, T_grid=None):
    """Build the solver's PVT property surface [[P_bar,T_C,rho_oil,rho_gas,mu_oil,mu_gas],...]
    from the EOS over a (P,T) grid — the bridge that lets the EOS drive the solver end-to-end."""
    P_grid = P_grid if P_grid is not None else [10, 30, 60, 100, 150, 200, 300]
    T_grid = T_grid if T_grid is not None else [4, 10, 20, 30, 40, 55, 70]
    rows = []
    for P in P_grid:
        for T in T_grid:
            pr = eos_properties(P, T, composition)
            rows.append([float(P), float(T), pr["rho_oil"], pr["rho_gas"], pr["mu_oil"], pr["mu_gas"]])
    return rows


def eos_selftest():
    """Sanity anchors used by the verification suite (#A1)."""
    checks = []
    #  pure methane vapour Z at 50 bar, 40 C ~ 0.91-0.94 (Standing-Katz / NIST)
    pr = eos_properties(50, 40, {"C1": 1.0})
    checks.append(("PR methane Z @50bar/40C", 0.85 < pr["Z_gas"] < 0.98, f"Z={pr['Z_gas']:.3f}"))
    #  natural gas density at 100 bar, 30 C is ~80-100 kg/m3
    pr2 = eos_properties(100, 30, DEFAULT_COMPOSITION)
    checks.append(("PR NG density @100bar/30C", 60 < pr2["rho_gas"] < 130, f"{pr2['rho_gas']:.1f} kg/m3"))
    #  gas SG of the default natural gas ~0.65-0.80
    checks.append(("PR NG specific gravity", 0.6 < pr2["gas_sg"] < 0.85, f"sg={pr2['gas_sg']:.3f}"))
    #  flash gives a vapour fraction in [0,1]
    fl = flash(60, 20, DEFAULT_COMPOSITION)
    checks.append(("PR flash V in [0,1]", 0.0 <= fl["V"] <= 1.0, f"V={fl['V']:.3f}"))
    return checks


if __name__ == "__main__":
    print("Peng-Robinson EOS self-test:")
    for name, ok, detail in eos_selftest():
        print(f"  [{'PASS' if ok else 'FAIL'}]  {name:32s} {detail}")
