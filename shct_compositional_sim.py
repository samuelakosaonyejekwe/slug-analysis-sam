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

#  Journal artwork carries no chart titles and needs the journal dpi; both come
#  from the environment so one switch covers every figure generator.
import os as _os


def _ttl(t):
    return t if _os.environ.get('SHCT_FIG_TITLES', '1') != '0' else ''
import shct_style as _S

_S.apply_style()
#  one export resolution for the whole project; this file defaulted to 150
_FIGDPI = _S.FIG_DPI

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
    def med(A):
        return np.nanmedian(A, 1)
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

    #  The consumed fraction is (gas to hydrate) / (gas fed), both over the SAME window.
    #  gas_consumed_hyd is a total over the whole run, and it used to be divided by ONE
    #  HOUR of inlet feed (moles_in * 3600), so the depletion a 48 h run reported was 48x
    #  the real one and the answer moved with t_end_h for a fixed physical case. The run
    #  already records the matching total, gas_in, so use it; fall back to the inlet rate
    #  times the run length for a result set that predates it.
    gas_in_kg_total = float(r.get("gas_in", 0.0))
    if gas_in_kg_total <= 0.0:
        rho_g_in = shct_eos.flash(float(P[0]), float(T[0]), comp)["rho_v"]
        gas_in_kg_total = (sv.case.operating.q_gas_insitu_inlet * rho_g_in
                           * float(sv.case.numerics.t_end_h) * 3600.0)
    moles_in = max(gas_in_kg_total / MW_g, 1e-9)
    #  cap the consumed fraction at a physical bound (hydrate removes a modest gas fraction)
    consumed_frac_total = float(np.clip(moles_consumed / moles_in, 0.0, 0.5))

    #  march component molar fluxes inlet->outlet; deplete formers by formability * cell weight.
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

    #  local flash state on the GRADED composition, at a sample of stations (a flash per
    #  cell is the expensive part). Sampled points are interpolated onto every station so
    #  the profile is usable; it was previously computed, left zero in between, and then
    #  dropped on the floor without ever reaching the report.
    stride = max(nx // 30, 1)
    isamp = list(range(0, nx, stride))
    if isamp[-1] != nx - 1:
        isamp.append(nx - 1)
    Vsamp = []
    for i in isamp:
        try:
            Vsamp.append(shct_eos.flash(float(P[i]), float(T[i]),
                                        {n: float(max(z_profile[i, j], 1e-9))
                                         for j, n in enumerate(names)})["V"])
        except Exception:
            Vsamp.append(np.nan)
    Vprof = np.interp(np.arange(nx), np.asarray(isamp, float), np.asarray(Vsamp, float))

    report = {"names": names, "z_inlet": z0.tolist(), "z_outlet": z_out.tolist(),
              "consumed_fraction_total": consumed_frac_total,
              "component_balance_residual": bal,
              "vapour_fraction_profile": [float(v) for v in Vprof],
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
                     if n in ("C1", "C2", "C3", "CO2", "N2", "nC4", "C7+")][:8] \
                or list(range(min(6, len(names))))
            #  PLOT THE GRADING, NOT THE LEVEL. On absolute z this panel drew seven flat
            #  horizontal lines: the grading is ~0.3 % of the mole fraction while the axis
            #  has to span 0 to 0.43, so the depletion -- the entire subject of the figure --
            #  was invisible and the panel said only "the composition is roughly constant".
            #  Normalising to each component's own inlet value keeps this panel's unique
            #  content, WHERE along the route the change happens, which the bar chart beside
            #  it (inlet vs outlet only) cannot show.
            for ci, j in enumerate(shown):
                _z0 = float(z_profile[0, j])
                _rel = 100.0 * (z_profile[:, j] / _z0 - 1.0) if abs(_z0) > 1e-12 \
                    else np.zeros(len(x_km))
                ax[0].plot(x_km, _rel, lw=1.6, color=palette[ci % len(palette)],
                           label=f"{names[j]}  (z$_0$ = {_z0:.4g})")
            ax[0].axhline(0.0, color="#3A5BA8", lw=0.6, ls=":")
            ax[0].set_xlabel("distance from wellhead  [km]")
            ax[0].set_ylabel(_S.label("change in z from inlet  [% of inlet value]",
                                      "Δz from inlet [%]"))
            ax[0].set_title(_ttl(_S.label(
                                "Compositional grading along line (hydrate former depletion)",
                                "Grading along line")),
                            color=NAVY, fontweight="bold", fontsize=9.5)
            #  OUTSIDE the axes (the project rule), and say that this panel draws a
            #  SUBSET: seven of the eleven components are plotted while the bar chart
            #  beside it shows all eleven, and nothing said so.
            ax[0].legend(fontsize=7, ncol=1, loc="upper left", bbox_to_anchor=(1.012, 1.0),
                         borderaxespad=0.0, framealpha=1.0, facecolor="white",
                         title=f"{len(shown)} of {len(names)} shown\n(omitted: "
                               + ", ".join(nm for k, nm in enumerate(names) if k not in shown)
                               + ")", title_fontsize=6.5)
            ax[0].grid(alpha=.25)
            dz = z_out - z0
            ax[1].bar(range(len(names)), dz, color=[RED if d < 0 else GREEN for d in dz])
            ax[1].set_xticks(range(len(names))); ax[1].set_xticklabels(names, rotation=45, fontsize=7)
            ax[1].set_ylabel("Δz (outlet − inlet)")
            ax[1].set_title(_ttl(_S.label(
                                "Net compositional change (− = depleted formers)",
                                "Net Δz (− = depleted)")),
                            color=NAVY, fontweight="bold", fontsize=9.5)
            ax[1].axhline(0, color="#3A5BA8", lw=0.6); ax[1].grid(alpha=.25, axis="y")
            #  "- = depleted formers" is only half true, and on a gas-rich composition it is
            #  actively misleading: z is a mole FRACTION, so consuming C2/C3 raises C1's
            #  fraction even while C1 moles are being consumed too. Say what the sign means.
            ax[1].text(0.5, -0.30, "Δz is a change in mole FRACTION: a component can rise here "
                       "while still being consumed,\nbecause removing the others renormalises "
                       "what is left. Total formers consumed: "
                       f"{100.0 * consumed_frac_total:.2f} % of feed moles.",
                       transform=ax[1].transAxes, ha="center", va="top",
                       fontsize=6.6, color="#555", linespacing=1.3)
            fig.tight_layout()
            _p = os.path.join(outdir, "compositional_transport.png")
            __import__("shct_style").screen(fig, _p)
            fig.savefig(_p, dpi=_FIGDPI)
            plt.close(fig)
    return report
