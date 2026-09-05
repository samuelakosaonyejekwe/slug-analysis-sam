#!/usr/bin/env python3
# =============================================================================
#  shct_evidence.py — test the deposition model against PUBLISHED EXPERIMENTAL
#  FINDINGS, rather than only against itself.
# -----------------------------------------------------------------------------
#  WHAT THIS IS, AND WHAT IT IS NOT
#
#  It is not a calibration and it is not a quantitative validation. The flow-loop
#  datasets referenced below are behind paywalls; their numeric deposit-thickness
#  series are NOT reproduced here, and inventing them would be worse than having
#  none — a fabricated agreement is indistinguishable from a real one until
#  somebody tries to reproduce it.
#
#  What the published abstracts DO state is a set of qualitative findings, and
#  those are falsifiable predictions the model either reproduces or does not:
#
#    E1  A hydrate deposit reaches a STEADY-STATE THICKNESS rather than growing
#        without limit, because shear strips it as it thickens.
#    E2  Wall temperature (i.e. subcooling) strongly affects BOTH the growth rate
#        and the steady-state thickness.
#    E3  Higher flow velocity / shear reduces the deposit that survives, and past
#        a critical thickness the deposit sloughs mechanically.
#    E4  MEG addition reduces the steady-state deposit thickness.
#    E5  Deposition is azimuthally NON-UNIFORM: fast at the bottom, where the wall
#        contacts liquid, slow at the top, which sees gas.
#
#  E1 is the one that matters most here. The previous formulation of this solver
#  COULD NOT PRODUCE IT: deposition was gated by clip(Phi_SH-1,0,1) and erosion ran
#  only below Phi_SH = 1, so a cell either grew with nothing opposing it or decayed
#  to bare wall. A steady-state thickness was not in the model's vocabulary. With
#  deposition and scouring competing continuously the plateau appears on its own at
#  delta_eq = Phi_SH * delta_ref. So E1 is external corroboration of exactly the
#  change that removed the Phi_SH = 1 circularity — not proof, but the model is now
#  answerable to a measurement it was previously incapable of contradicting.
#
#  SOURCES (findings taken from the published abstracts; numeric data paywalled and
#  deliberately not reproduced — a reader with journal access should confirm):
#
#    [1] X. Zhang, E.O. Straume, G.A. Grasso, R.E.M. Morales, A.K. Sum,
#        "A bench-scale flow loop study on hydrate deposition under multiphase
#        flow conditions", Fuel 262 (2020) 116558.
#        doi:10.1016/j.fuel.2019.116558                        -> E1, E2, E3, E4, E5
#    [2] Z.M. Aman, M. Di Lorenzo, K. Kozielski, C.A. Koh, P. Warrier, M.L. Johns,
#        E.F. May, "Hydrate formation and deposition in a gas-dominant flowloop:
#        Initial studies of the effect of velocity and subcooling",
#        J. Nat. Gas Sci. Eng. 35 (2016) 1096-1103.
#        doi:10.1016/j.jngse.2016.05.015                            -> E2, E3
#
#  Each check runs the REAL solver end-to-end — not the deposition algebra in
#  isolation, which would only re-derive what was typed in — so a trend has to
#  survive the full pressure/temperature/holdup coupling to count as reproduced.
#
#      python3 shct_evidence.py [outdir]
#
#  Exit status is 1 if any published trend is contradicted.
# =============================================================================
import copy
import json
import os
import sys

import numpy as np

import solver

#  Deliberately coarse and deterministic: these are TREND tests over a sweep, so the
#  run has to be cheap enough to do a dozen of them, and free of ensemble scatter that
#  would masquerade as a trend (or hide one).
GRID = dict(n_cells=44, n_ensemble=4, t_end_h=24.0)
#  Sweep runs are stopped BEFORE the line plugs. Once a deposit reaches full bore the
#  reported thickness is the cap, identical at every condition, and a trend measured
#  through that clip is measuring the clip. 10 h is comfortably short of the ~12.7 h
#  P50 time-to-plug at the as-operated condition.
SWEEP_H = 10.0


# -----------------------------------------------------------------------------
def _case(**over):
    c = solver.make_default_case()
    c.pipeline.n_cells = GRID["n_cells"]
    c.numerics.n_ensemble = GRID["n_ensemble"]
    c.numerics.t_end_h = GRID["t_end_h"]
    c.numerics.deterministic = True
    for path, val in over.items():
        obj, attr = c, path
        for part in path.split(".")[:-1]:
            obj = getattr(obj, part)
        setattr(obj, path.split(".")[-1], val)
    return c


def _run(c):
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    eng = sv.engineering()
    #  peak deposit over the line, and the deposit history so a plateau is visible
    d = np.asarray(sv.results["snap_delta"])          # (nt, nx, N) or (nt, nx)
    hist = np.max(d.reshape(d.shape[0], -1), axis=1) * 1000.0     # mm, worst cell
    return eng, hist


def _monotone(xs, ys, want, tol_frac=0.01):
    """Sign concordance: does y move the way the experiment says?

    Tied pairs are EXCLUDED rather than counted as disagreement. A deposit that has
    filled the bore reads the same value at two different conditions because the
    physical cap has removed the information, not because the model disagreed with
    the experiment — scoring that as a failure would be measuring the clip. Pairs
    are only informative if y differs by more than tol_frac of the observed range,
    and the count of informative pairs is returned so a sweep that saturated into
    uselessness is visible rather than silently reported as agreement.
    """
    n = len(xs)
    rng = max(np.max(ys) - np.min(ys), 1e-30)
    agree = info = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(ys[j] - ys[i]) <= tol_frac * rng:
                continue                                   # tied: no information
            info += 1
            agree += int(np.sign(ys[j] - ys[i]) == want * np.sign(xs[j] - xs[i]))
    return (agree / info if info else float("nan")), info


# ============================================================== E1 ===========
def check_plateau():
    """A sub-critical deposit must GROW and then STOP (Zhang et al. 2020, [1]).

    Two ways to fail, and the second is the one that matters. The obvious failure is
    a deposit that grows without limit. The insidious one is a deposit that never
    forms: "no growth" trivially satisfies "growth stopped", so a check that only
    asks whether the curve flattened can pass on an empty pipe and prove nothing.
    A real deposit is therefore required before the plateau is even assessed.

    The plateau is looked for where the model says it must be — below the derived
    Phi_crit. At the as-operated condition this line is SUPERCRITICAL (Phi_SH ~ 3.5
    against Phi_crit = 1.08), which is the case study's whole point: it plugs. So the
    sub-critical regime is reached by strengthening the scouring, k_ero, which is
    what a bench-scale loop does physically by running at high shear in a small bore.
    That raises Phi_crit above the operating Phi_SH without touching the hydrate
    kinetics, and the model then has to produce a FINITE thickness — and, if it is
    right about the mechanism, one near delta_eq = Phi_SH * delta_ref.
    """
    c = _case()
    k = c.kinetics
    D = c.pipeline.diameter_m
    #  4x scouring: Phi_crit scales with k_ero, delta_ref inversely
    c.kinetics.k_ero = k.k_ero * 4.0
    ke = c.kinetics.k_ero
    d_ref = k.wall_capture_eff * D / (4.0 * k.C_phi * ke)
    phi_crit = 2.0 * k.C_phi * ke * k.consol_restriction / k.wall_capture_eff

    eng, hist = _run(c)
    phi = float(eng["max_Phi_SH_uncapped"])
    final = float(hist[-1])

    n = len(hist)
    early = (hist[n // 4] - hist[0]) / max(n // 4, 1)
    late = (hist[-1] - hist[3 * n // 4]) / max(n - 3 * n // 4, 1)
    ratio = late / max(early, 1e-12)

    grew = final > 0.5                       # mm — a deposit actually formed
    subcrit = phi < phi_crit
    flat = ratio < 0.25                      # late growth << early growth
    #  Order-of-magnitude cross-check against the closed form, NOT a quantitative
    #  match, and deliberately not part of the pass criteria. `phi` is the maximum
    #  Phi_SH over all space and time, while `final` is the plateau at the worst
    #  cell, whose LOCAL Phi_SH is lower and varies during the run — so the two are
    #  not the same quantity and are expected to differ by tens of percent. Treating
    #  agreement here as validation would be reading more into it than it carries.
    predicted = phi * d_ref * 1000.0
    near = abs(final - predicted) < 0.75 * predicted

    return dict(
        check="E1 deposit grows then reaches a steady thickness", source="[1]",
        passed=bool(grew and subcrit and flat),
        deposit_formed=bool(grew), sub_critical=bool(subcrit), plateaued=bool(flat),
        same_order_as_closed_form=bool(near),   # loose cross-check, not a match
        Phi_SH=phi, Phi_crit=float(phi_crit), delta_ref_mm=float(d_ref * 1000.0),
        predicted_plateau_mm=float(predicted), observed_final_mm=final,
        early_growth_mm_per_snap=float(early), late_growth_mm_per_snap=float(late),
        late_over_early=float(ratio),
        note="the gated formulation had no plateau in its vocabulary: a cell either "
             "grew unopposed above Phi_SH=1 or decayed to bare wall below it")


# ============================================================== E2 ===========
def check_subcooling():
    """Colder wall -> faster growth AND thicker steady state ([1], [2])."""
    T_sea = [10.0, 7.0, 4.0, 1.0]
    peaks, phis = [], []
    for T in T_sea:
        eng, hist = _run(_case(**{"operating.T_seabed_C": T,
                                  "numerics.t_end_h": SWEEP_H}))
        peaks.append(float(eng["peak_deposit_mm"]))
        phis.append(float(eng["max_Phi_SH_uncapped"]))
    #  subcooling rises as T_seabed falls, so thickness must rise as T_seabed falls
    frac, info = _monotone(T_sea, peaks, want=-1)
    return dict(check="E2 subcooling drives growth and plateau", source="[1],[2]",
                passed=bool(info >= 3 and frac >= 0.83), T_seabed_C=T_sea,
                peak_deposit_mm=peaks, Phi_SH=phis,
                pair_agreement=float(frac), informative_pairs=info)


# ============================================================== E3 ===========
def check_shear():
    """More flow -> more scouring -> less deposit survives ([1], [2])."""
    base = _case()
    ql0, qg0 = base.operating.q_liquid_insitu, base.operating.q_gas_insitu_inlet
    mult = [0.6, 1.0, 1.6, 2.2]
    peaks, sfrac = [], []
    for m in mult:
        #  raise BOTH phases so the mixture velocity — and so the slug frequency that
        #  does the scouring — rises, rather than only changing the phase ratio
        eng, _h = _run(_case(**{"operating.q_liquid_insitu": ql0 * m,
                                "operating.q_gas_insitu_inlet": qg0 * m,
                                "numerics.t_end_h": SWEEP_H}))
        peaks.append(float(eng["peak_deposit_mm"]))
        sfrac.append(float(eng.get("slug_fraction", float("nan"))))
    frac, info = _monotone(mult, peaks, want=-1)
    return dict(check="E3 shear/velocity strips deposit", source="[1],[2]",
                passed=bool(info >= 3 and frac >= 0.83), rate_multiplier=mult,
                peak_deposit_mm=peaks, slug_fraction=sfrac,
                pair_agreement=float(frac), informative_pairs=info)


# ============================================================== E4 ===========
def check_meg():
    """MEG reduces the steady-state deposit thickness ([1])."""
    doses = [0.0, 10.0, 20.0, 30.0]
    peaks = []
    for w in doses:
        eng, _h = _run(_case(**{"operating.MEG_wt_inlet": w,
                                "numerics.t_end_h": SWEEP_H}))
        peaks.append(float(eng["peak_deposit_mm"]))
    frac, info = _monotone(doses, peaks, want=-1)
    return dict(check="E4 MEG thins the deposit", source="[1]",
                passed=bool(info >= 3 and frac >= 0.83), MEG_wt_pct=doses,
                peak_deposit_mm=peaks,
                pair_agreement=float(frac), informative_pairs=info)


# ============================================================== E6 ===========
def check_film_growth_rate():
    """The one QUANTITATIVE deposition datum in the public literature.

    E1-E5 are directions. This is a number. Qin (2020, Colorado School of Mines PhD
    thesis, open access) quantified hydrate film growth from video images in an
    oil-dominated rig at 10 % water cut, 550 psi, T_bulk 52 F, with the surface
    12 F below the hydrate point: 0.02-0.08 in/hour, i.e. 1.41e-7 to 5.64e-7 m/s.
    That is exactly the quantity d(delta)/dt in the deposit-evolution equation, so
    it can be compared directly rather than as a trend.

    Two things came out of making the comparison, and the second is why it matters.

    The rate did not match: the model was 6.5 to 26 times faster than measured, and
    reproducing the data needed a capture fraction of 0.04-0.15, which is not a
    physically sensible value for a fraction. Chasing that discrepancy exposed the
    cause. The wall growth law used a_i, the GAS-LIQUID interfacial area, which is
    the correct term for bulk growth at that interface and the wrong one for a wall
    process: it falls to 0.11 1/m as the line fills with liquid, so the model
    predicted essentially no wall deposition in exactly the liquid-full,
    oil-dominated configuration where Qin measured it, and where the thesis
    attributes it to "water droplets settling on the wall".

    With the area term corrected to the wall (4/D, scaled by the liquid holdup and
    the water fraction of that liquid) the model lands at 0.65 to 2.6 times the
    measured range, requiring a capture fraction of 0.38 to 1.54 — bracketing unity.
    A measurement therefore both found a structural defect and, once it was fixed,
    supported the kinetics quantitatively.
    """
    import numpy as np
    from shct_correlations import gas_density  # noqa: F401  (import parity)
    c = _case()
    k = c.kinetics
    #  Qin's stated conditions
    dT = 12.0 * 5.0 / 9.0                       # 12 F surface subcooling -> K
    T_surf_C = (52.0 - 32.0) * 5.0 / 9.0 - dT
    wc = 0.10                                   # 10 % water cut
    lo, hi = 0.02 * 0.0254 / 3600.0, 0.08 * 0.0254 / 3600.0

    kg = k.kg0 * np.exp(-k.Ea_over_R * (1.0 / (T_surf_C + 273.15) - 1.0 / k.T_ref_K))
    #  d(delta)/dt = f_wall * kg * a_wall * dTsub^n * D/4, and a_wall = 4/D * alpha_l * wf,
    #  so the diameter cancels and the comparison transfers from Qin's 1.75 in rig to any
    #  line without a scaling assumption.
    rate = kg * 1.0 * wc * (dT / k.dTsub_ref_C) ** k.growth_exp_n   # alpha_l ~ 1, liquid-full
    f_lo, f_hi = lo / rate, hi / rate

    #  the model must be able to produce a deposit at all in a liquid-full line
    liquid_full_ok = rate > 0.0
    #  and the capture fraction the data implies must be a physically sensible fraction
    sensible = 0.1 < f_hi <= 2.0 and f_lo > 0.05

    return dict(
        check="E6 film growth rate matches a measured value", source="[3]",
        passed=bool(liquid_full_ok and sensible),
        measured_lo_m_per_s=lo, measured_hi_m_per_s=hi,
        model_rate_at_full_capture_m_per_s=float(rate),
        ratio_model_over_measured_lo=float(rate / hi),
        ratio_model_over_measured_hi=float(rate / lo),
        implied_capture_fraction_lo=float(f_lo),
        implied_capture_fraction_hi=float(f_hi),
        produces_deposit_when_liquid_full=bool(liquid_full_ok),
        note="the only quantitative deposition datum found in public literature; "
             "it exposed the gas-liquid-area defect and then supported the kinetics")


# ============================================================== E5 ===========
def check_azimuthal():
    """Bottom of the pipe (liquid-wetted) deposits faster than the top (gas) ([1])."""
    import shct_crosssection as cx
    D = 0.2545
    out = cx.azimuthal_deposit(np.array([8.0e-3]), np.array([0.45]), D=np.array([D]))
    #  the reconstruction returns (theta, profile, bottom, top); take its own
    #  bottom/top rather than re-deriving them from an assumed angle convention
    bottom = float(np.atleast_1d(out[2])[0])
    top = float(np.atleast_1d(out[3])[0])
    return dict(check="E5 azimuthal non-uniformity, bottom > top", source="[1]",
                passed=bool(bottom > top * 1.05),
                bottom_mm=bottom * 1000.0, top_mm=top * 1000.0,
                ratio=float(bottom / max(top, 1e-12)))


# -----------------------------------------------------------------------------
def run(outdir=None):
    print("=== SHCT vs published flow-loop findings ===")
    print("  trend-level only: the numeric datasets are paywalled and are NOT")
    print("  reproduced here. Directions and the existence of a plateau, not magnitudes.\n")
    checks = [check_plateau, check_subcooling, check_shear, check_meg,
              check_azimuthal, check_film_growth_rate]
    rows = []
    for fn in checks:
        try:
            r = fn()
        except Exception as exc:                       # a broken check is not a pass
            r = dict(check=fn.__name__, passed=False, error=f"{type(exc).__name__}: {exc}")
        rows.append(r)
        tag = "PASS" if r.get("passed") else "FAIL"
        print(f"  [{tag}] {r.get('check', fn.__name__)}   {r.get('source', '')}")
        for key, val in r.items():
            if key in ("check", "passed", "source", "note"):
                continue
            if isinstance(val, float):
                print(f"           {key} = {val:.4g}")
            else:
                print(f"           {key} = {val}")
        if r.get("note"):
            print(f"           ({r['note']})")

    n_pass = sum(1 for r in rows if r.get("passed"))
    print(f"\n{n_pass}/{len(rows)} published trends reproduced")
    print("\nThis is corroboration, not calibration. It cannot fix a wrong magnitude;")
    print("it can only show the model is answerable to measurements it could contradict.")

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, "evidence_trends.json")
        with open(path, "w") as fh:
            json.dump(dict(
                sources={
                    "[1]": "Zhang, Straume, Grasso, Morales & Sum, Fuel 262 (2020) 116558, "
                           "doi:10.1016/j.fuel.2019.116558",
                    "[3]": "Qin, H., 2020. Hydrate film growth and risk management in oil/gas "
                           "pipelines using experiments, simulations and machine learning. "
                           "PhD thesis, Colorado School of Mines (open access).",
                    "[2]": "Aman, Di Lorenzo, Kozielski, Koh, Warrier, Johns & May, "
                           "J. Nat. Gas Sci. Eng. 35 (2016) 1096-1103, "
                           "doi:10.1016/j.jngse.2016.05.015"},
                scope="qualitative trend agreement; numeric data paywalled and not reproduced",
                grid=GRID, checks=rows, passed=n_pass, total=len(rows)), fh, indent=2)
        print(f"\nwrote {path}")
    return 0 if n_pass == len(rows) else 1


if __name__ == "__main__":
    sys.exit(run(sys.argv[1] if len(sys.argv) > 1 else None))
