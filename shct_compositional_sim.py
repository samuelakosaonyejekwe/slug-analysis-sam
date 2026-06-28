#!/usr/bin/env python3
# =============================================================================
#  shct_compositional_sim.py — SEQUENTIAL COMPOSITIONAL TRANSPORT (item 8)
# -----------------------------------------------------------------------------
#  Tracks how the HYDROCARBON COMPOSITION grades ALONG THE LINE as gas hydrate
#  preferentially removes the light hydrate-formers (C1, C2, C3, CO2, ...) from
#  the gas. Given the SHCT-solved fields (per-cell hydrate activity and the total
#  gas mass consumed by hydrate), it marches the component molar fluxes from inlet
#  to outlet, depleting each component in proportion to its hydrate-formability,
#  and reports the composition profile z_k(x), the depleted formers, and the
#  local PR-flash state from the SHCT-supplied (p, T).
#
#  HONEST SCOPE: this is a SEQUENTIAL (one-way: hydraulics/thermal -> composition)
#  reduced compositional model — components are advected and depleted on the
#  converged SHCT flow field, with a local equilibrium flash. It is NOT a fully
#  IMPLICITLY-COUPLED compositional reservoir/pipeline simulator (where composition,
#  phase behaviour, holdup and pressure are solved simultaneously). It captures the
#  practical flow-assurance effect — compositional grading from hydrate former
#  depletion — at screening cost; component moles are conserved (feed = out + consumed).
# =============================================================================
from __future__ import annotations
import os
import numpy as np
import shct_eos

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:                                       # pragma: no cover
    _HAVE_MPL = False

NAVY = "#2E5BBF"; ACCENT = "#1F8AC0"; RED = "#E0463C"; ORANGE = "#E8842B"
GREEN = "#3FA65A"; TEAL = "#1AA0A0"; PURPLE = "#8E5CC8"
#  relative hydrate-formability of the gas formers (heavier sI/sII formers captured preferentially);
#  non-formers (N2, and the heavy ends) are ~inert and ENRICH as formers deplete.
FORMABILITY = {"C1": 1.0, "C2": 1.6, "C3": 2.6, "iC4": 2.7, "nC4": 1.2,
               "CO2": 1.3, "H2S": 2.2, "N2": 0.15}


def simulate_composition(sv, outdir=None):
    """Sequential compositional transport on the solved SHCT field. Returns a dict with the
    per-station composition profile and writes csv_compositional_transport.csv + a chart."""
    comp = getattr(sv.case.fluids, "composition", None) or shct_eos.DEFAULT_COMPOSITION
    names, z0 = shct_eos._normalise(comp)
    z0 = np.asarray(z0, float)
    r = sv.results
    med = lambda A: np.nanmedian(A, 1)
    x_km = sv.x / 1000.0
    P = med(r["p"]); T = med(r["T"])
    nx = len(x_km)

    #  per-cell hydrate ACTIVITY weight (where formers are being consumed): bulk hydrate + deposit.
    Dpipe = sv.case.pipeline.diameter_m
    activity = med(r["phi"]) + 2.0 * med(r["delta"]) / Dpipe
    if activity.sum() <= 0:
        activity = np.zeros(nx)
    w_cell = activity / activity.sum() if activity.sum() > 0 else np.zeros(nx)

    #  total hydrocarbon-gas moles consumed by hydrate over the line (from the SHCT gas sink):
    #  gas_consumed_hyd is kg/m integrated * ... stored as total kg in results.
    gas_consumed_kg = float(r.get("gas_consumed_hyd", 0.0))
    MW_g = max(float(sv.case.fluids.gas_MW), 1e-3)               # kg/mol
    moles_consumed = gas_consumed_kg / MW_g                       # total gas moles to hydrate

    #  inlet molar basis: 1 'unit' of feed gas; scale the consumption to a per-unit fraction so the
    #  profile is physical even when the absolute feed-rate molar basis is unknown.
    rho_g_in = shct_eos.flash(float(P[0]), float(T[0]), comp)["rho_v"]
    A0 = np.pi * Dpipe ** 2 / 4.0
    gas_in_kg = sv.case.operating.q_gas_insitu_inlet * rho_g_in   # kg/s
    moles_in = max(gas_in_kg / MW_g, 1e-9)
    #  cap the consumed fraction at a physical bound (hydrate removes a modest gas fraction)
    consumed_frac_total = float(np.clip(moles_consumed / (moles_in * 3600.0 + 1e-9), 0.0, 0.5))

    #  march component molar fluxes inlet->outlet; deplete formers by formability * cell weight.
    F = np.outer(np.ones(nx), z0).copy()                         # (nx, ncomp) overall mole fractions
    form = np.array([FORMABILITY.get(n, 0.0) for n in names])
    Fcur = z0.copy()                                             # current flux composition (mol-frac basis)
    consumed_k = np.zeros(len(names))
    z_profile = np.zeros((nx, len(names)))
    for i in range(nx):
        z_profile[i] = Fcur / max(Fcur.sum(), 1e-12)
        # remove this cell's share of the total consumption, weighted by formability * presence
        cell_consume = consumed_frac_total * w_cell[i]
        if cell_consume > 0 and form.sum() > 0:
            weights = form * Fcur
            if weights.sum() > 0:
                dF = cell_consume * weights / weights.sum()
                dF = np.minimum(dF, Fcur * 0.95)               # never deplete a component below ~0
                Fcur = Fcur - dF
                consumed_k += dF
    z_out = Fcur / max(Fcur.sum(), 1e-12)

    #  component conservation check (feed = outlet flux + consumed), on the mole-fraction basis
    feed_total = float(z0.sum())
    bal = float(abs(feed_total - (Fcur.sum() + consumed_k.sum())))

    #  local flash state along the line (properties / vapour fraction) on the graded composition
    Vprof = np.zeros(nx)
    for i in range(0, nx, max(nx // 30, 1)):
        try:
            Vprof[i] = shct_eos.flash(float(P[i]), float(T[i]),
                                      {n: float(max(z_profile[i, j], 1e-9)) for j, n in enumerate(names)})["V"]
        except Exception:
            Vprof[i] = np.nan

    report = {"names": names, "z_inlet": z0.tolist(), "z_outlet": z_out.tolist(),
              "consumed_fraction_total": consumed_frac_total,
              "component_balance_residual": bal,
              "grading_max_abs_dz": float(np.max(np.abs(z_out - z0)))}

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        cols = ["x_km", "P_bar", "T_C"] + [f"z_{n}" for n in names]
        with open(os.path.join(outdir, "csv_compositional_transport.csv"), "w") as fh:
            fh.write(",".join(cols) + "\n")
            for i in range(nx):
                row = [x_km[i], P[i], T[i]] + list(z_profile[i])
                fh.write(",".join(f"{v:.6g}" for v in row) + "\n")
        if _HAVE_MPL:
            palette = [RED, ORANGE, GREEN, TEAL, PURPLE, NAVY, ACCENT, "#9AA8C7", "#E0463C", "#2E5BBF"]
            fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
            shown = [j for j, n in enumerate(names)
                     if n in ("C1", "C2", "C3", "CO2", "N2", "nC4", "C7+")][:8] or list(range(min(6, len(names))))
            for ci, j in enumerate(shown):
                ax[0].plot(x_km, z_profile[:, j], lw=1.6, color=palette[ci % len(palette)], label=names[j])
            ax[0].set_xlabel("distance (km)"); ax[0].set_ylabel("overall mole fraction z")
            ax[0].set_title("Compositional grading along line (hydrate former depletion)",
                            color=NAVY, fontweight="bold", fontsize=9.5)
            ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=.25)
            dz = z_out - z0
            ax[1].bar(range(len(names)), dz, color=[RED if d < 0 else GREEN for d in dz])
            ax[1].set_xticks(range(len(names))); ax[1].set_xticklabels(names, rotation=45, fontsize=7)
            ax[1].set_ylabel("Δz (outlet − inlet)")
            ax[1].set_title("Net compositional change (− = depleted formers)",
                            color=NAVY, fontweight="bold", fontsize=9.5)
            ax[1].axhline(0, color="#3A5BA8", lw=0.6); ax[1].grid(alpha=.25, axis="y")
            fig.tight_layout(); fig.savefig(os.path.join(outdir, "compositional_transport.png"), dpi=150)
            plt.close(fig)
    return report
