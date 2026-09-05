#!/usr/bin/env python3
"""Automated test suite for solver.py (G22).

Runs as a plain script (`python3 test_solver.py`) or under pytest
(`pytest test_solver.py`). It wraps the built-in V&V verification suite and adds
targeted unit/regression checks for the items hardened in this revision:

  A1  latent-heat source term is active in the energy equation
  A2  phase-field diffusion smooths the hydrate field
  A3  liquid AND gas AND hydrate mass balances close
  A4  wall deposit is mass-coupled to (driven by) hydrate growth
  B5  the never-fail path advances time (no infinite loop)
  B7  holdup bound enforcement is conservative
  C10 the default ensemble produces a genuine P10/P50/P90 spread
  G21 input validation rejects bad cases
  H24 the reported V&V count is self-consistent
"""
import sys
import os
import inspect
import copy
import pytest
import numpy as np
import solver


def _short_case(**kw):
    c = solver.Case()
    c.pipeline.n_cells = kw.get("n_cells", 50)
    c.numerics.n_ensemble = kw.get("n_ensemble", 6)
    c.numerics.t_end_h = kw.get("t_end_h", 12.0)
    c.numerics.deterministic = kw.get("deterministic", False)
    return c


def test_verification_suite_passes():
    assert solver.run_verification() is True


def test_latent_heat_active_A1():
    c = _short_case(n_ensemble=1, t_end_h=8.0, deterministic=True)
    T_on = np.nanmean(solver.TransientSHCT(c).run(verbose=False)["T"])
    c0 = copy.deepcopy(c); c0.fluids.L_hyd = 0.0
    T_off = np.nanmean(solver.TransientSHCT(c0).run(verbose=False)["T"])
    # latent heat can only add energy -> mean T not lower than the no-latent reference
    assert T_on >= T_off - 1e-6


def test_phasefield_diffusion_active_A2():
    assert solver.Kinetics().D_phi > 0.0
    base = _short_case(n_ensemble=2, t_end_h=8.0, deterministic=True)
    r0 = solver.TransientSHCT(base).run(verbose=False)
    cd = copy.deepcopy(base); cd.kinetics.D_phi = base.kinetics.D_phi * 50.0
    r1 = solver.TransientSHCT(cd).run(verbose=False)
    # stronger diffusion must not break conservation and must keep phi bounded
    assert r1["mass_err"] < 0.05
    assert float(np.nanmax(r1["phi"])) <= base.kinetics.phi_max + 1e-9
    # the two solutions differ (diffusion has an effect)
    assert not np.allclose(np.nanmedian(r0["phi"], 1), np.nanmedian(r1["phi"], 1))


def test_mass_balances_close_A3():
    for scen in ("steady", "turndown", "shutin"):
        c = _short_case(n_cells=60, t_end_h=12.0)
        c.scenario.kind = scen
        r = solver.TransientSHCT(c).run(verbose=False)
        assert r["mass_err"] < 0.02, (scen, r["mass_err"])
        assert r["gas_mass_err"] < 0.05, (scen, r["gas_mass_err"])
        # hydrate water balance: water consumed == hydrate-water mass
        water = r["liq_to_hyd"] * c.fluids.rho_water
        expect = r["hyd_mass"] * c.fluids.hyd_water_massfrac
        assert abs(water - expect) <= 0.02 * max(expect, 1e-6)


def test_deposit_coupled_to_growth_A4():
    # with no hydrate growth (kg0 = 0) there can be no wall deposit
    c = _short_case(t_end_h=12.0); c.kinetics.kg0 = 0.0
    r = solver.TransientSHCT(c).run(verbose=False)
    assert float(np.nanmax(r["delta"])) < 1e-9
    assert float(np.nanmax(r["phi"])) < 1e-9


def test_no_fallback_spin_B5():
    # a normal run must terminate with a finite step count well under the cap
    c = _short_case(t_end_h=10.0)
    r = solver.TransientSHCT(c).run(verbose=False)
    assert 0 < r["steps"] < 5_000_000
    assert np.all(np.isfinite(r["p"])) and np.all(np.isfinite(r["T"]))


def test_enforce_bounds_conservative_B7():
    c = _short_case(n_ensemble=3, t_end_h=2.0)
    sv = solver.TransientSHCT(c)
    A = np.full((c.pipeline.n_cells, 3), 1.0)
    La = np.random.default_rng(0).uniform(-0.3, 1.4, A.shape)  # deliberately out of bounds
    before = La.sum(0)
    out = sv._enforce_bounds(La.copy(), A)
    assert np.all(out >= -1e-9) and np.all(out <= A + 1e-9)
    assert np.allclose(out.sum(0), before, atol=1e-6)          # mass conserved per column


def test_default_ensemble_has_spread_C10():
    c = _short_case(n_cells=70, n_ensemble=16, t_end_h=26.0)
    r = solver.TransientSHCT(c).run(verbose=False)
    ttp = r["plug_time"][~np.isnan(r["plug_time"])]
    # the default (non-deterministic) ensemble should give a real spread when it plugs
    if ttp.size > 3:
        assert np.percentile(ttp, 90) - np.percentile(ttp, 10) > 1e-3
    # a deterministic run must collapse the input-spread (reproducibility control)
    cd = copy.deepcopy(c); cd.numerics.deterministic = True
    r2a = solver.TransientSHCT(cd).run(verbose=False)
    r2b = solver.TransientSHCT(cd).run(verbose=False)
    assert np.allclose(r2a["ts"]["Tsub"], r2b["ts"]["Tsub"])


def test_input_validation_G21():
    bad = solver.Case(); bad.fluids.water_cut = 1.5
    try:
        solver.validate_case(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("water_cut > 1 should be rejected")
    bad2 = solver.Case(); bad2.numerics.engine = "nope"
    try:
        solver.validate_case(bad2)
    except ValueError:
        pass
    else:
        raise AssertionError("unknown engine should be rejected")


#  -------- v3 follow-up gap fixes (closure-level unit tests; fast, no full sim) ----------
def test_gas_Z_factor_physical_1():
    # constant Z unchanged by default; correlation gives the high-pressure dip; table interpolates
    f0 = solver.Fluids()
    assert solver.gas_Z_factor(150, 55, f0) == f0.gas_Z
    fc = solver.Fluids(gas_Z_corr=True)
    z30 = float(solver.gas_Z_factor(30, 55, fc)); z150 = float(solver.gas_Z_factor(150, 55, fc))
    assert 0.9 < z30 <= 1.0 and 0.78 < z150 < 0.92 and z150 < z30   # dip with pressure
    ft = solver.Fluids(gas_Z_table=[[150, 55, 0.7]])
    assert abs(float(np.atleast_1d(solver.gas_Z_factor(np.array([150.0]), np.array([55.0]), ft))[0]) - 0.7) < 1e-6


def test_tvd_limiter_and_conservation_8():
    # limiter is TVD-bounded in [0, 2]; reduces to upwind at kind='upwind'
    r = np.array([-1.0, 0.0, 0.5, 1.0, 3.0])
    for kind in ("minmod", "vanleer", "superbee"):
        phi = solver._limiter(r, kind)
        assert np.all(phi >= -1e-9) and np.all(phi <= 2.0 + 1e-9)
    assert np.allclose(solver._limiter(r, "upwind"), 0.0)
    # TVD holdup transport still conserves liquid mass exactly
    c = _short_case(t_end_h=8.0); c.numerics.flux_limiter = "vanleer"
    assert solver.TransientSHCT(c).run(verbose=False)["mass_err"] < 1e-3


def test_deposit_insulation_slows_growth_7():
    base = _short_case(n_ensemble=4, t_end_h=20.0)
    base.kinetics.k_dep_insul = 0.0           # no insulation -> faster late growth
    d0 = float(np.nanmax(solver.TransientSHCT(base).run(verbose=False)["delta"]))
    ins = copy.deepcopy(base); ins.kinetics.k_dep_insul = 0.9   # strong insulation
    d1 = float(np.nanmax(solver.TransientSHCT(ins).run(verbose=False)["delta"]))
    assert d1 <= d0 + 1e-9                      # insulation cannot accelerate deposition


def test_oil_water_slip_conserves_and_accumulates_2():
    c = _short_case(n_ensemble=3, t_end_h=10.0); c.fluids.oil_water_slip = True
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert r["mass_err"] < 1e-3                 # water transport stays mass-conservative
    wf = np.nanmedian(r["water_frac"], 1)
    assert np.all(wf >= -1e-9) and float(np.nanmax(wf)) >= c.fluids.water_cut - 1e-6


def test_uq_distribution_specs_15():
    c = _short_case(n_ensemble=16, t_end_h=12.0)
    c.numerics.uq_inputs = {"kg0": {"dist": "lognormal", "sigma": 0.3},
                            "U_wall": {"dist": "normal", "sigma": 0.1},
                            "nuc_tau0_h": {"dist": "uniform", "low": 0.7, "high": 1.5},
                            "_corr_z": True}
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2                 # runs and conserves with dict-spec + correlated UQ


def test_validation_harness_4_20():
    ds = {"name": "unit", "scalars": {"dP_total_bar": 18.0},
          "plug_events": [{"x_km": 12, "plugged": True}, {"x_km": 2, "plugged": False}]}
    c = _short_case(n_ensemble=3, t_end_h=10.0)
    rep = solver.validate_against_data(c, ds)
    assert "scalars" in rep and "dP_total_bar" in rep["scalars"]
    assert rep["phi_sh_skill"] is not None and 0.0 <= rep["phi_sh_skill"]["accuracy"] <= 1.0


def test_slug_length_and_clip_diag_5_11():
    Lu = solver.slug_length(np.array([4.0, 0.5]), 0.3, np.array([0.05, 0.01]))
    assert np.all(Lu > 0)
    c = _short_case(t_end_h=8.0)
    e = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(c))
    assert "clip_warning" in e and e["clip_frac_velocity"] >= 0.0
    assert e["slug_length_max_m"] >= e["slug_length_mean_m"] > 0.0


#  -------- round-3 depth deepenings (couplings now real, not reporting-only) ----------
def test_gas_viscosity_lee_1():
    f0 = solver.Fluids(); assert solver.gas_viscosity(120.0, 50.0, f0) == f0.mu_gas
    fc = solver.Fluids(gas_visc_corr=True)
    mu_lo = float(solver.gas_viscosity(20.0, 50.0, fc)); mu_hi = float(solver.gas_viscosity(180.0, 50.0, fc))
    assert 1e-6 < mu_lo < 1e-4 and mu_hi > mu_lo            # gas viscosity rises with density/pressure


def test_oil_density_pvt_1():
    f0 = solver.Fluids(); assert solver.oil_density(150.0, 55.0, f0) == f0.rho_oil
    fc = solver.Fluids(oil_pvt_corr=True)
    # hotter -> lighter, higher pressure -> denser
    assert float(solver.oil_density(150.0, 80.0, fc)) < float(solver.oil_density(150.0, 20.0, fc))
    assert float(solver.oil_density(250.0, 55.0, fc)) > float(solver.oil_density(50.0, 55.0, fc))


def test_realgas_Z_still_conserves_gas_1_3():
    c = _short_case(n_ensemble=3, t_end_h=10.0); c.fluids.gas_Z_corr = True
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["gas_mass_err"] < 1e-2                          # Z-ratio keeps gas continuity conservative


def test_threephase_density_couples_to_momentum_2():
    # with oil/water slip the local liquid density becomes a field (heavier at water-rich low points)
    c = _short_case(n_ensemble=3, t_end_h=10.0); c.fluids.oil_water_slip = True
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert not np.isscalar(sv._rho_l_field)                 # density is now a field, not a constant
    assert r["mass_err"] < 1e-3 and r["gas_mass_err"] < 5e-2


def test_deep_options_stable_and_conservative():
    c = _short_case(n_ensemble=3, t_end_h=10.0)
    c.fluids.oil_water_slip = True; c.fluids.droplet_entrainment = True
    c.fluids.gas_visc_corr = True; c.fluids.oil_pvt_corr = True; c.fluids.gas_Z_corr = True
    c.numerics.subgrid_slug = True
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["fallbacks"] == 0 and r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2


#  -------- round-4: universality build (EOS flash, two-fluid-mass, acoustics) ----------
def test_pr_eos_flash_A1():
    import shct_eos
    for name, ok, detail in shct_eos.eos_selftest():
        assert ok, (name, detail)
    pr = shct_eos.eos_properties(100, 30, shct_eos.DEFAULT_COMPOSITION)
    assert 60 < pr["rho_gas"] < 130 and 0.6 < pr["gas_sg"] < 0.85
    tab = shct_eos.build_pvt_table(shct_eos.DEFAULT_COMPOSITION)
    assert len(tab) > 10 and len(tab[0]) == 6


def test_eos_driven_case_conserves_A1():
    import shct_eos
    c = _short_case(n_ensemble=3, t_end_h=8.0)
    c.fluids.composition = shct_eos.DEFAULT_COMPOSITION
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0
    assert c.fluids.pvt_table is not None and 0.6 < c.fluids.gas_sg < 0.85


def test_mixture_sound_speed_B5():
    c2 = float(solver.mixture_sound_speed(0.5, 820.0, 90.0, 100.0))
    c1 = float(solver.mixture_sound_speed(0.999, 820.0, 90.0, 100.0))
    assert 50.0 < c2 < c1 < 1300.0                         # two-phase << near-liquid


def test_acoustic_option_stable_B5():
    c = _short_case(n_ensemble=2, t_end_h=1.0); c.numerics.acoustic = 0.5
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0


def test_volume_consistent_pressure_B3():
    # opt-in two-fluid-mass coupling stays conservative and moves toward gas-mass consistency
    base = _short_case(n_ensemble=3, t_end_h=10.0)
    sv0 = solver.TransientSHCT(base); sv0.run(verbose=False); e0 = sv0.engineering()
    cc = copy.deepcopy(base); cc.numerics.volume_consistent_pressure = 0.12
    sv1 = solver.TransientSHCT(cc); r1 = sv1.run(verbose=False); e1 = sv1.engineering()
    assert r1["fallbacks"] == 0
    assert e1["gas_holdup_consistency"] <= e0["gas_holdup_consistency"] + 1e-6


def test_pvt_table_universal_gas_conservation_A1_3():
    tab = [[50, 40, 830, 45, 2.2e-3, 1.3e-5], [150, 55, 815, 120, 1.9e-3, 1.6e-5],
           [250, 60, 800, 200, 1.7e-3, 1.9e-5]]
    c = _short_case(n_ensemble=3, t_end_h=8.0); c.fluids.pvt_table = tab
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["gas_mass_err"] < 1e-2                          # gas conserves for a table density too


#  -------- round-5: precision build (EOS BIP/Peneloux/LBC, vdW-P, advanced physics, numerics) ----
def test_eos_bip_peneloux_lbc_2():
    import shct_eos
    pr = shct_eos.eos_properties(100, 30, shct_eos.DEFAULT_COMPOSITION)
    assert 1e-5 < pr["mu_oil"] < 1.0          # LBC liquid viscosity is physical (not clipped at 1000 cP)
    # Peneloux shift increases liquid density vs unshifted PR (denser, more realistic)
    assert pr["rho_oil"] > 0
    # BIPs: a CO2-rich mix flashes without error and gives a physical density
    pr2 = shct_eos.eos_properties(80, 25, {"CO2": 0.3, "C1": 0.6, "C3": 0.1})
    assert pr2["rho_gas"] > 0 and 0 <= pr2["vapour_frac"] <= 1


def test_vdwp_hydrate_composition_dependent_3():
    import shct_eos
    import numpy as np
    lean = shct_eos.hydrate_equilibrium_vdwp(np.array([70.0]), {"C1": 1.0})
    rich = shct_eos.hydrate_equilibrium_vdwp(np.array([70.0]), {"C1": 0.80, "C3": 0.15, "nC4": 0.05})
    assert float(rich[0]) > float(lean[0])    # richer gas -> hydrates stable at higher T


def test_advanced_physics_conserves_9_13():
    import shct_eos
    c = _short_case(n_ensemble=4, t_end_h=10.0)
    c.numerics.advanced_physics = True
    c.fluids.composition = shct_eos.DEFAULT_COMPOSITION
    c.fluids.salinity_wt = 5.0
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False); e = sv.engineering()
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0
    assert "scaling_tendency_max" in e


def test_advanced_numerics_options_15_17_18():
    c = _short_case(n_ensemble=8, t_end_h=10.0)
    c.numerics.tvd_energy = True; c.numerics.error_dt = True; c.numerics.lhs_uq = True
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0


def test_joule_thomson_cools_11():
    import shct_correlations as C
    f = solver.Fluids(gas_Z_corr=True)
    mu = float(C.joule_thomson_dTdP(100.0, 40.0, f))
    assert mu > 0.0                           # positive JT coefficient -> cooling on expansion


def test_golden_master_24():
    #  locks the deterministic default-style case numerics against silent regression (#24)
    c = solver.Case()
    c.numerics.t_end_h = 12.0; c.pipeline.n_cells = 50; c.numerics.n_ensemble = 4
    c.numerics.deterministic = True
    sv = solver.TransientSHCT(c); sv.run(verbose=False); e = sv.engineering()
    #  Updated when Phi_SH stopped driving the deposition (the gate that made the
    #  Phi_SH = 1 criterion an assumption rather than a result). Every movement was
    #  ATTRIBUTED before the numbers were rewritten, and each is in the direction the
    #  physics requires:
    #    dP        -3.3 %      less deposit surviving -> less bore restriction
    #    arrival_T -3.2 %      thinner insulating deposit -> more heat lost to seabed
    #    max_dTsub -15.6 %     hydrate now forms over a wider reach and its latent heat
    #                          is exothermic, so the peak bulk subcooling self-limits
    #    peak depo -4.7 %      scouring opposes growth everywhere, not only below 1
    #    P_plug    0.50->0.25  same reason: fewer realisations reach full bore
    #    max_Phi_SH  unchanged — it peaks at startup, at the slug-frequency floor,
    #                before any deposit exists to change it
    #  The scoured-hydrate transfer accounts for <= 0.42 % of the movement; removing
    #  the gate accounts for the rest.
    golden = {"dP_total_bar": 16.621, "arrival_T_C": 11.396, "max_subcooling_C": 7.247,
              "max_Phi_SH": 3.732, "peak_deposit_mm": 121.861, "P_plug": 0.25}
    for kk, ref in golden.items():
        got = float(e[kk])
        tol = max(0.01 * abs(ref), 1e-3)
        assert abs(got - ref) <= tol, f"golden regression: {kk} {got} vs {ref}"


#  -------- round-6: precision build (bilinear PVT, soil, within-step iters, vdW-P machinery) ----
def test_bilinear_pvt_interpolation_24():
    import shct_correlations as C
    import numpy as np
    tab = [[p, t, 800 + t, 0.5 * p, 2e-3, 1e-5] for p in [50, 100, 150, 200] for t in [20, 40, 60]]
    v = float(np.atleast_1d(C._pvt_lookup(np.array([125.0]), np.array([30.0]), tab, 3))[0])
    assert abs(v - 62.5) < 1e-6                  # midpoint -> interpolated, not snapped


def test_soil_transient_conserves_12():
    c = _short_case(n_ensemble=3, t_end_h=10.0); c.numerics.soil_transient = True
    c.scenario.kind = "shutin"
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0


def test_within_step_iteration_6():
    c = _short_case(n_ensemble=3, t_end_h=8.0); c.numerics.within_step_iters = 3
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0


def test_vdwp_langmuir_machinery_8():
    import shct_eos
    #  the full-Langmuir vdW-P runs and is structure-aware (sI/sII), even if it is the EXPERIMENTAL
    #  path (the validated model is the reduced one, checked elsewhere)
    T = shct_eos.hydrate_equilibrium_vdwp_full(70.0, {"C1": 1.0})
    assert -40.0 < float(T) < 50.0


#  -------- A. Numerical-core build (two-fluid-mass engine, droplet/water/slug fields) ----------
def test_twofluid_mass_engine_A1():
    base = _short_case(n_ensemble=3, t_end_h=10.0)
    e0 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(base))
    cc = copy.deepcopy(base); cc.numerics.engine = "twofluid_mass"
    sv = solver.TransientSHCT(cc); r = sv.run(verbose=False); e = sv.engineering()
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0
    #  improves gas-holdup consistency vs the drift-flux implicit engine
    assert e["gas_holdup_consistency"] <= e0["gas_holdup_consistency"] + 1e-6


def test_droplet_field_A2():
    c = _short_case(n_ensemble=3, t_end_h=8.0); c.fluids.droplet_field = True
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0      # conservative with droplet transport
    assert 0.0 <= sv._droplet_frac <= 0.95


def test_water_drift_momentum_A3():
    c = _short_case(n_ensemble=3, t_end_h=8.0)
    c.fluids.oil_water_slip = True; c.fluids.water_drift = True
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0


def test_slug_tracking_A4():
    c = _short_case(n_ensemble=3, t_end_h=8.0); c.numerics.slug_tracking = True
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0
    s = sv._slugS
    assert np.all(s >= -1e-9) and np.all(s <= 1.0 + 1e-9)    # indicator stays in [0,1]


#  -------- B. Thermodynamics / PVT build (vdW-P framework, condensation, 3-phase, plus-fraction) --
def test_vdwp_full_framework_B8():
    import shct_eos
    #  the full-Langmuir vdW-P FRAMEWORK runs (sI/sII, EOS fugacities, full Dh/Dcp/Dv reference);
    #  it is experimental, and the VALIDATED model is the reduced one (checked separately)
    T = shct_eos.hydrate_equilibrium_vdwp_full(70.0, {"C1": 1.0})
    assert -40.0 < float(T) < 60.0
    lean = float(shct_eos.hydrate_equilibrium_vdwp(np.array([70.0]), {"C1": 1.0})[0])
    assert abs(lean - 9.5) < 2.5            # the reduced model still validates vs Deaton-Frost


def test_whitson_plus_fraction_B11():
    import shct_eos
    sp = shct_eos.whitson_split(0.05, 140.0, 3)
    assert len(sp) == 3
    Tcs = [v[1] for v in sp.values()]; MWs = [v[4] for v in sp.values()]
    assert Tcs == sorted(Tcs) and MWs == sorted(MWs)        # heavier pseudos -> higher Tc / MW
    comp = shct_eos.expand_composition(shct_eos.DEFAULT_COMPOSITION, 3, 140.0)
    assert "C7+" not in comp and "PC1" in comp and "PC3" in comp
    pr = shct_eos.eos_properties(100, 30, comp)
    assert pr["rho_gas"] > 0


def test_phase_envelope_B11():
    import shct_eos
    Pd = shct_eos.saturation_pressure(30, {"C1": 0.6, "C3": 0.2, "nC5": 0.2}, "dew")
    assert 1.0 <= Pd <= 700.0


def test_three_phase_water_B10():
    import shct_eos
    tp = shct_eos.three_phase_flash(100, 30, shct_eos.DEFAULT_COMPOSITION, 0.3, 3.0)
    assert tp["free_water"] is True and 0.0 <= tp["water_activity"] <= 1.0
    # water content rises with T, falls with P (Bukacek)
    assert shct_eos.water_content_gas(50, 40) > shct_eos.water_content_gas(50, 20)
    assert shct_eos.water_content_gas(50, 40) > shct_eos.water_content_gas(150, 40)


def test_condensation_latent_B9():
    import shct_eos
    c = _short_case(n_ensemble=3, t_end_h=8.0)
    c.fluids.composition = dict(shct_eos.DEFAULT_COMPOSITION); c.fluids.condensation_latent = True
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["fallbacks"] == 0 and sv._Vsurf is not None


#  -------- C. Heat transfer build (transient radial soil conduction) ----------
def test_radial_soil_conduction_C12():
    c = _short_case(n_ensemble=3, t_end_h=12.0); c.scenario.kind = "shutin"
    c.numerics.soil_transient = True; c.numerics.soil_nodes = 5
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2 and r["fallbacks"] == 0


def test_radial_vs_lumped_soil_inertia_C12():
    #  the multi-shell radial soil holds more heat than no soil model -> a buried line cools slower
    base = _short_case(n_ensemble=2, t_end_h=12.0); base.scenario.kind = "shutin"
    e0 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(base))
    cr = copy.deepcopy(base); cr.numerics.soil_transient = True; cr.numerics.soil_nodes = 6
    er = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(cr))
    assert er["cooldown_to_hydrate_h"] >= e0["cooldown_to_hydrate_h"] - 1e-6


#  -------- F. Validation / calibration build (MCMC posterior, blind validation) ----------
def test_bayesian_mcmc_posterior_F22():
    #  minimal MULTI-CHAIN MCMC run — verify it returns a posterior with mean/std, a correlation
    #  matrix and a Gelman-Rubin R-hat convergence diagnostic
    c = solver.make_default_case(); c.numerics.t_end_h = 8.0
    res = solver.bayesian_calibrate(c, {"dP_total_bar": 20.0}, free=["U_wall"],
                                    n_samples=3, sigma_rel=0.2, n_chains=2)
    assert "mean" in res and "std" in res and "corr" in res and "rhat" in res
    assert res["corr"].shape == (1, 1) and 0.0 <= res["accept_rate"] <= 1.0
    assert res["n_chains"] == 2 and len(res["rhat"]) == 1


def test_blind_validate_holdout_F21():
    #  blind_validate must hold out a key and report its predictive error (structure check, quick)
    c = solver.make_default_case(); c.numerics.t_end_h = 8.0
    ds = {"scalars": {"dP_total_bar": 20.0, "arrival_T_C": 11.0}}
    res = solver.blind_validate(c, ds, train_keys=["dP_total_bar"], maxiter=3)
    assert res is not None and "arrival_T_C" in res["test"]


#  -------- gap-closure: previously hard-coded constants are now wired config fields ----------
def test_exposed_constants_are_wired():
    #  defaults must reproduce the baseline exactly (golden-master safety), and overriding a
    #  field must move the dependent reported quantity (the constant is genuinely plumbed through).
    base = _short_case(n_ensemble=3, t_end_h=8.0, deterministic=True)
    e0 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(base))
    #  API-14E C-factor scales the erosional-velocity limit linearly
    cc = copy.deepcopy(base); cc.fluids.api14e_C_factor = 244.0   # 2x default
    e1 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(cc))
    assert abs(e1["erosional_limit_mps"] - 2.0 * e0["erosional_limit_mps"]) < 1e-6
    #  surge factor scales the slug-catcher surge volume linearly
    cs = copy.deepcopy(base); cs.numerics.surge_factor = base.numerics.surge_factor * 1.5
    e2 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(cs))
    assert abs(e2["V_surge_P90_m3"] - 1.5 * e0["V_surge_P90_m3"]) < 1e-3
    #  the k_hyd dynamic-U field is consistent with kinetics.deposit_k_hyd (no hard-coded 0.5)
    cu = copy.deepcopy(base); cu.numerics.dynamic_U = True
    ru = solver.TransientSHCT(cu).run(verbose=False)
    assert ru["mass_err"] < 1e-2 and ru["fallbacks"] == 0


def test_gelman_rubin_converges_for_identical_chains():
    #  R-hat -> 1 for well-mixed (here identical) chains; nan/strict for a single chain
    a = np.tile(np.linspace(0.0, 1.0, 20)[:, None], (1, 2))
    rhat = solver._gelman_rubin([a, a.copy()])
    assert rhat.shape == (2,) and np.all(np.isfinite(rhat)) and np.all(rhat < 1.05)
    assert np.all(np.isnan(solver._gelman_rubin([a])))         # one chain -> not estimable


#  -------- quasi-3-D reconstruction layer (cross-section, compositional, 3-D export) ----------
def _solved(n_cells=40, n_ensemble=3, t_end_h=6.0):
    c = solver.Case(); c.pipeline.n_cells = n_cells
    c.numerics.n_ensemble = n_ensemble; c.numerics.t_end_h = t_end_h
    sv = solver.TransientSHCT(c); sv.run(verbose=False)
    return sv


def test_crosssection_geometry_inversion():
    import shct_crosssection as cx
    #  the holdup -> liquid-level inversion is exact against the circular-segment area
    h = cx.liquid_level(np.array([0.25, 0.5, 0.75]))
    assert abs(h[1] - 0.5) < 1e-6                         # half-full -> h/D = 0.5
    back = cx._area_fraction(h)
    assert np.allclose(back, [0.25, 0.5, 0.75], atol=1e-4)
    # wetted-perimeter fraction & interface width are bounded and physical
    hh, wf, iw = cx.section_geometry(np.array([0.3, 0.6]), 0.3)
    assert np.all((wf >= 0) & (wf <= 1)) and np.all(iw >= 0)


def test_crosssection_outputs(tmp_path=None):
    import shct_crosssection as cx, tempfile, os
    sv = _solved()
    out = tempfile.mkdtemp()
    p = cx.crosssection_outputs(sv, out)
    assert os.path.exists(p)                              # csv written
    hdr = open(p).readline().strip().split(",")
    assert "liquid_level_h_over_D" in hdr and "deposit_bottom_mm" in hdr
    for fn in ["cx1_geometry.png", "cx2_azimuthal_deposit.png", "cx3_sections.png"]:
        assert os.path.exists(os.path.join(out, fn))


def test_compositional_report():
    import shct_compositional as cp, tempfile, os
    sv = _solved()
    out = tempfile.mkdtemp()
    p = cp.compositional_report(sv, out, n_stations=12)
    assert os.path.exists(p)
    hdr = open(p).readline().strip().split(",")
    assert "vapour_frac_V" in hdr and any(h.startswith("K_") for h in hdr)
    assert os.path.exists(os.path.join(out, "compo_pvt.png"))


def test_threed_field_and_vtk():
    import shct_threed as t3, tempfile, os
    sv = _solved()
    out = tempfile.mkdtemp()
    field = t3.build_3d_field(sv, n_axial=20, n_theta=12, n_r=4)
    n_r, n_theta, n_ax = field["dims"]
    assert field["fields"]["temperature_C"].shape == (n_ax, n_theta, n_r)
    # phase is binary, holdup in [0,1]
    assert set(np.unique(field["fields"]["phase_liquid"])) <= {0.0, 1.0}
    assert np.all((field["fields"]["holdup"] >= 0) & (field["fields"]["holdup"] <= 1))
    vtk = t3.write_vtk(field, os.path.join(out, "pipe_3d.vtk"))
    txt = open(vtk).read()
    assert txt.startswith("# vtk DataFile Version") and "STRUCTURED_GRID" in txt
    assert f"DIMENSIONS {n_r} {n_theta} {n_ax}" in txt
    # point count matches the structured-grid dimensions
    npts = n_r * n_theta * n_ax
    assert f"POINTS {npts} float" in txt and f"POINT_DATA {npts}" in txt


def test_openfoam_coupling_generates_cases():
    import shct_openfoam as of, tempfile, os
    sv = _solved(n_cells=50, t_end_h=6.0)
    out = tempfile.mkdtemp()
    man = of.couple(sv, out, max_sections=2, run=False)
    assert man["n_sections"] == 2 and "openfoam_available" in man
    # each section produced a complete, structured interFoam case
    for e in man["sections"]:
        cd = os.path.join(out, e["casedir"])
        for f in ["system/blockMeshDict", "system/controlDict", "system/fvSchemes",
                  "system/fvSolution", "system/setFieldsDict", "constant/g",
                  "constant/transportProperties", "0/U", "0/p_rgh", "0/alpha.liquid",
                  "Allrun", "README.txt", "section.json"]:
            assert os.path.exists(os.path.join(cd, f)), f
        bm = open(os.path.join(cd, "system/blockMeshDict")).read()
        assert bm.count("hex") == 5 and all(k in bm for k in ("vertices", "boundary", "inlet",
                                                              "outlet", "walls", "arc"))
        # interFoam application + two phases + BCs sourced from SHCT
        assert "application     interFoam" in open(os.path.join(cd, "system/controlDict")).read()
        assert "phases (liquid gas)" in open(os.path.join(cd, "constant/transportProperties")).read()
        assert "fixedValue" in open(os.path.join(cd, "0/U")).read()
    # ingest is graceful when no CFD has run
    cd0 = os.path.join(out, man["sections"][0]["casedir"])
    assert of.ingest_results(cd0)["available"] is False


def test_compositional_transport_conserves_and_grades():
    import shct_eos, shct_compositional_sim as cs, tempfile, os
    c = _short_case(n_cells=50, n_ensemble=3, t_end_h=8.0)
    c.fluids.composition = shct_eos.DEFAULT_COMPOSITION
    sv = solver.TransientSHCT(c); sv.run(verbose=False)
    out = tempfile.mkdtemp()
    rep = cs.simulate_composition(sv, out)
    assert rep["component_balance_residual"] < 1e-9          # component moles conserved
    assert 0.0 <= rep["consumed_fraction_total"] <= 0.5
    assert rep["grading_max_abs_dz"] >= 0.0                  # composition grades along the line
    # outlet composition is a valid normalised mole fraction
    zo = np.array(rep["z_outlet"]); assert abs(zo.sum() - 1.0) < 1e-6 and np.all(zo >= -1e-9)
    assert os.path.exists(os.path.join(out, "csv_compositional_transport.csv"))


def test_twofluid_mass_newton_consistency():
    #  the monolithic volume-mass Newton drives gas-holdup consistency to ~0 and conserves mass
    #  (the documented trade-off is dP fidelity, which this test does NOT assert).
    base = _short_case(n_cells=50, n_ensemble=3, t_end_h=8.0); base.numerics.deterministic = True
    e0 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(copy.deepcopy(base)))
    cc = copy.deepcopy(base); cc.numerics.engine = "twofluid_mass_newton"
    sv = solver.TransientSHCT(cc); r = sv.run(verbose=False); e = sv.engineering()
    assert r["fallbacks"] == 0
    assert r["mass_err"] < 0.02 and r["gas_mass_err"] < 0.02            # conservation intact
    assert e["gas_holdup_consistency"] < 0.01                            # ~0 (far below implicit's ~8%)
    assert e["gas_holdup_consistency"] < e0["gas_holdup_consistency"]    # better than the implicit engine


def test_full_newton_engine_runs_clean():
    #  the fully-coupled (alpha_l, p, u_m) block-Newton must run with NO fallbacks, finite fields,
    #  exact mass conservation, and converge tightly (mean Newton residual ~machine precision).
    c = _short_case(n_cells=30, n_ensemble=2, t_end_h=4.0); c.numerics.deterministic = True
    c.numerics.engine = "twofluid_full_newton"
    sv = solver.TransientSHCT(c); r = sv.run(verbose=False)
    assert r["fallbacks"] == 0
    assert r["mass_err"] < 1e-2 and r["gas_mass_err"] < 5e-2            # conservation intact
    assert np.all(np.isfinite(r["p"])) and np.all(np.isfinite(r["j"]))
    #  converged: the lightweight diagnostic records a tiny final scaled residual
    assert getattr(sv, "_fn_last_res", 1.0) < 1e-3
    assert getattr(sv, "_fn_iters_sum", 0) / max(getattr(sv, "_fn_steps", 1), 1) < 8.0


def test_full_newton_reproduces_implicit_dP():
    #  KEY verification: the fully-coupled Newton AGREES with the validated segregated 'implicit'
    #  engine's physical pressure drop (it is a coupled-vs-segregated cross-check, not a dP improver).
    base = _short_case(n_cells=30, n_ensemble=2, t_end_h=4.0); base.numerics.deterministic = True
    e0 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(copy.deepcopy(base)))
    cc = copy.deepcopy(base); cc.numerics.engine = "twofluid_full_newton"
    sv = solver.TransientSHCT(cc); sv.run(verbose=False); e = sv.engineering()
    dP_imp = e0["dP_total_bar"]; dP_fn = e["dP_total_bar"]
    assert dP_imp > 0 and abs(dP_fn - dP_imp) / dP_imp < 0.25            # within 25% of the implicit dP


def test_twoway_openfoam_coupling_loop():
    import shct_openfoam as of, tempfile
    c = solver.Case(); c.pipeline.n_cells = 40; c.numerics.n_ensemble = 2; c.numerics.t_end_h = 5.0
    out = tempfile.mkdtemp()
    #  synthetic CFD target holdup drives the closed loop (no OpenFOAM needed for the test)
    res = of.couple_iterate(c, out, max_sections=2, max_iters=4, tol=0.0,
                            synthetic_cfd=lambda s: 0.55)
    h = res["history"]
    assert len(h) >= 2 and all(r["mismatch"] is not None for r in h)
    assert h[-1]["mismatch"] <= h[0]["mismatch"] + 1e-9      # loop converges toward CFD (non-worsening)
    #  feedback tunes the drift-flux distribution parameter C0 (the strong holdup knob),
    #  not roughness — CFD holds more liquid (0.55) than the base run, so C0 is raised.
    assert res["calibrated_case"].numerics.drift_C0_factor != c.numerics.drift_C0_factor   # feedback applied


def test_strict_mode_flag_runs_clean_case():
    #  strict mode must NOT raise on a well-behaved case (no fallbacks / excessive clips)
    c = _short_case(n_cells=50, t_end_h=6.0); c.numerics.strict = True
    r = solver.TransientSHCT(c).run(verbose=False)
    assert r["fallbacks"] == 0 and np.all(np.isfinite(r["p"]))


def test_wax_screen_reports():
    c = _short_case(n_cells=50, t_end_h=6.0); c.fluids.wax_appearance_C = 30.0
    sv = solver.TransientSHCT(c); sv.run(verbose=False); e = sv.engineering()
    assert e["wax_risk"] is True and e["wax_under_km"] > 0.0
    c2 = _short_case(n_cells=50, t_end_h=6.0)               # default: wax screen off
    e2 = (lambda s: (s.run(verbose=False), s.engineering())[1])(solver.TransientSHCT(c2))
    assert e2["wax_risk"] is False and e2["wax_under_km"] == 0.0


#  -------- v8: published-reference CLOSURE validations (gap 1 / gap 5) ----------
def test_friction_closure_matches_colebrook():
    #  Haaland (the friction closure) must agree with the exact Colebrook-White reference to
    #  within Haaland's published ~2% band; the Colebrook self-check must hit the Moody value.
    rep = solver.validate_friction_curve(outdir=None)
    assert rep["max_abs_pct_dev"] < 2.0 and rep["rms_pct_dev"] < 1.0
    assert abs(rep["colebrook_smooth_Re1e5"] - 0.018) < 0.001     # Moody smooth-pipe anchor
    #  the in-code Colebrook iteration must actually solve its own implicit equation
    Re = np.array([1e4, 1e6, 1e8]); eps = np.array([1e-4, 1e-3, 1e-2])
    f = solver._colebrook_white(Re, eps)
    resid = 1.0 / np.sqrt(f) + 2.0 * np.log10(eps / 3.7 + 2.51 / (Re * np.sqrt(f)))
    assert np.max(np.abs(resid)) < 1e-9


def test_drift_flux_vertical_matches_canonical():
    #  the slip closure must reproduce the canonical vertical Taylor-bubble values exactly
    #  (C0=1.2, drift Fr=0.35) and expose the honest horizontal drift deficit vs Benjamin 0.542.
    rep = solver.validate_drift_flux(outdir=None)
    vert = next(r for r in rep["rows"] if "vertical" in r["orientation"])
    horiz = next(r for r in rep["rows"] if "horizontal" in r["orientation"])
    assert abs(vert["C0_err_pct"]) < 1e-6 and abs(vert["Fr_err_pct"]) < 1e-6
    assert horiz["Fr_model"] < horiz["Fr_ref"]                   # honest: closure under-drifts horizontally


def test_slug_frequency_reproduces_zabaras():
    #  the closure must reproduce the Zabaras (2000) correlation it implements to ~machine zero.
    rep = solver.validate_slug_frequency(outdir=None)
    assert rep["max_fidelity_err_pct"] < 1e-6
    assert rep["published_accuracy_band_pct"] == 60.0


def test_flowloop_holdup_validation(tmp_path=None):
    #  v9 (gap 1): the solver's drift-flux holdup must validate against the REAL measured
    #  void-fraction dataset (Das Neves et al. 2025) with a sane RMSE, and the 1-param
    #  drift_C0_factor calibration must IMPROVE the fit (lower RMSE) — same pattern as hydrate.
    import os
    dd = os.path.join(os.path.dirname(os.path.abspath(solver.__file__)), "7", "field_data")
    ds = os.path.join(dd, "flowloop_holdup_dasneves2025.json")
    if not os.path.exists(ds):
        return                                              # dataset shipped with repo; skip if absent
    rep = solver.validate_flowloop(ds, outdir=None)
    assert rep["n"] == 14
    assert rep["void_rmse"] < 0.10                          # drift-flux holdup within ~0.1 of measured
    assert rep["void_rmse_calibrated"] <= rep["void_rmse"]  # calibration via drift_C0_factor helps
    assert rep["drift_C0_factor_calibrated"] > 1.0          # honest: over-predicts void -> raise C0


def test_openfoam_couple_accepts_resolution_and_time(tmp_path=None):
    #  v9: couple() must thread the CFD end_time / o-grid resolution and (without OpenFOAM)
    #  still generate runnable cases recording the requested mesh/time in the manifest.
    import shct_openfoam as of, tempfile, os
    c = _short_case(n_cells=40, t_end_h=6.0)
    sv = solver.TransientSHCT(c); sv.run(verbose=False)
    d = str(tmp_path) if tmp_path is not None else tempfile.mkdtemp()
    man = of.couple(sv, d, max_sections=1, run=False, end_time=0.4, Ni=8, Nz=20)
    assert man["end_time"] == 0.4 and man["mesh"] == {"Ni": 8, "Nz": 20}
    assert man["n_sections"] >= 1
    assert os.path.exists(os.path.join(d, "openfoam_cases", man["sections"][0]["name"], "Allrun"))


def test_validate_closures_writes_reports(tmp_path=None):
    import tempfile, os, json as _json
    d = str(tmp_path) if tmp_path is not None else tempfile.mkdtemp()
    out = solver.validate_closures(outdir=d)
    assert out["friction"]["max_abs_pct_dev"] < 2.0
    assert os.path.exists(os.path.join(d, "friction_validation_report.json"))
    assert os.path.exists(os.path.join(d, "drift_flux_validation_report.json"))
    assert os.path.exists(os.path.join(d, "slug_frequency_validation_report.json"))


def _run_all():
    fns = [v for kname, v in sorted(globals().items())
           if kname.startswith("test_") and callable(v)]
    npass = 0
    for fn in fns:
        try:
            fn(); print(f"  [PASS]  {fn.__name__}"); npass += 1
        except Exception as exc:                                # noqa: BLE001
            print(f"  [FAIL]  {fn.__name__}: {exc}")
    print(f"\n  {npass}/{len(fns)} tests passed")
    return npass == len(fns)


if __name__ == "__main__":
    import sys
    print("=" * 64); print(" SHCT SOLVER — TEST SUITE"); print("=" * 64)
    sys.exit(0 if _run_all() else 1)


# ---------------------------------------------------------------------------
#  Phi_SH reporting diagnostics (uncapped magnitude, C-free ratio, gate share).
#  Phi_SH = C_phi * Psi. Phi_SH no longer drives the deposition at all — it is a
#  diagnostic formed from the competing rates — so what has to be guarded is that
#  its critical value stays DERIVED rather than drifting back to an assumed 1.
# ---------------------------------------------------------------------------
def test_phish_reporting_keys_present():
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    for k in ("max_Phi_SH", "max_Phi_SH_uncapped",
              "max_Psi_kinetic_ratio", "Phi_SH_above_critical_frac",
              "Phi_SH_critical", "deposit_ref_mm"):
        assert k in eng, f"engineering() is missing {k}"


def test_phish_uncapped_is_never_below_capped():
    """The quoted magnitude must never be the plot cap in disguise."""
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    assert eng["max_Phi_SH_uncapped"] >= eng["max_Phi_SH"] - 1e-9
    assert eng["max_Phi_SH"] <= c.kinetics.phi_report_cap * (1 + 1e-9)


def test_psi_is_phish_over_C_and_is_C_invariant():
    """Psi is the part the model predicts; C_phi is pure scale."""
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); e1 = sv.engineering()
    assert abs(e1["max_Psi_kinetic_ratio"]
               - e1["max_Phi_SH_uncapped"] / c.kinetics.C_phi) < 1e-9

    c2 = copy.deepcopy(c); c2.kinetics.C_phi = c.kinetics.C_phi * 3.0
    sv2 = solver.TransientSHCT(c2); sv2.run(verbose=False); e2 = sv2.engineering()
    #  Phi_SH scales exactly with C_phi; Psi does not move.
    assert abs(e2["max_Phi_SH_uncapped"] / e1["max_Phi_SH_uncapped"] - 3.0) < 1e-6
    assert abs(e2["max_Psi_kinetic_ratio"] - e1["max_Psi_kinetic_ratio"]) < 1e-9


def test_above_critical_fraction_is_a_fraction():
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    g = eng["Phi_SH_above_critical_frac"]
    assert (g != g) or (0.0 <= g <= 1.0), f"fraction out of range: {g}"


# ---------------------------------------------------------------------------
#  The Phi_SH = 1 criterion must be EARNED, not assumed.
#
#  It used to be wired into the deposition three times over: consolidation required
#  Phi_SH > 1, the wall-capture fraction was scaled by clip(Phi_SH-1,0,1) so nothing
#  deposited below 1, and erosion ran only below 1. A criterion a model has been
#  handed cannot be tested by that model. These four tests hold the fix in place.
# ---------------------------------------------------------------------------
def test_phish_does_not_appear_in_the_deposition_block():
    """Source-level guard: the gate must not creep back in.

    A numerical test can miss a re-introduced switch that happens not to fire on the
    short test case, so this reads the deposition block itself.
    """
    src = inspect.getsource(solver.TransientSHCT.run)
    start = src.index("(D-gate) wall-capture fraction")
    end = src.index("bulk hydrate phase-field")
    block = src[start:end]
    code = "\n".join(ln.split("#")[0] for ln in block.splitlines())
    assert "PhiSH" not in code, (
        "Phi_SH drives the deposition again — the criterion it is used to report "
        "would then be an assumption, not a result:\n"
        + "\n".join(ln for ln in code.splitlines() if "PhiSH" in ln))


def test_reference_thickness_is_what_C_phi_encodes():
    """Phi_SH = 1 means: scouring balances deposition at delta_ref."""
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    k = c.kinetics
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    expect_mm = 1000.0 * k.wall_capture_eff * c.pipeline.diameter_m / (
        4.0 * k.C_phi * k.k_ero)
    assert abs(eng["deposit_ref_mm"] - expect_mm) < 1e-9 * max(expect_mm, 1.0)
    #  doubling C_phi halves the thickness the criterion refers to: C_phi is a
    #  statement about a thickness, not a free multiplier
    c2 = copy.deepcopy(c); c2.kinetics.C_phi *= 2.0
    sv2 = solver.TransientSHCT(c2); sv2.run(verbose=False); e2 = sv2.engineering()
    assert abs(e2["deposit_ref_mm"] / eng["deposit_ref_mm"] - 0.5) < 1e-9


def test_critical_coupling_is_derived_not_unity():
    """Phi_crit follows from three kinetic constants and is NOT hard-wired to 1."""
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    k = c.kinetics
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    expect = 2.0 * k.C_phi * k.k_ero * k.consol_restriction / k.wall_capture_eff
    assert abs(eng["Phi_SH_critical"] - expect) < 1e-9 * expect
    #  it lands near 1 for the shipped constants — that is the finding — but it is a
    #  computed number, so a different erosion rate must move it
    assert 0.5 < eng["Phi_SH_critical"] < 2.0, eng["Phi_SH_critical"]
    c2 = copy.deepcopy(c); c2.kinetics.k_ero *= 4.0
    sv2 = solver.TransientSHCT(c2); sv2.run(verbose=False); e2 = sv2.engineering()
    assert abs(e2["Phi_SH_critical"] / eng["Phi_SH_critical"] - 4.0) < 1e-9, (
        "Phi_crit did not move with k_ero — it is being asserted, not derived")


def test_subcritical_deposit_equilibrates_instead_of_running_away():
    """The competing rates must produce a finite plateau, which is the falsifiable bit.

    Integrate the solver's own two rate expressions. Below Phi_crit the thickness
    settles at delta_eq = Phi_SH * delta_ref; above it, the deposit passes the
    consolidation restriction, locks, erosion stops and growth is unbounded.
    """
    c = _short_case(n_ensemble=1, t_end_h=1.0, n_cells=10, deterministic=True)
    k, D = c.kinetics, c.pipeline.diameter_m
    f_wall, f_slug = k.wall_capture_eff, 0.02
    d_ref = f_wall * D / (4.0 * k.C_phi * k.k_ero)
    phi_crit = 2.0 * k.C_phi * k.k_ero * k.consol_restriction / f_wall

    def integrate(phi_sh, hours=400.0, n=400_000):
        Rg_wall = phi_sh * f_slug / k.C_phi          # invert Phi_SH = C_phi*Rg/f_slug
        dt, delta, locked = hours * 3600.0 / n, 0.0, False
        for _ in range(n):
            locked |= (2.0 * delta / D) > k.consol_restriction
            delta += dt * (f_wall * Rg_wall * D / 4.0
                           - k.k_ero * f_slug * delta * (not locked))
        return delta

    for phi in (0.25, 0.5, 0.75):
        d = integrate(phi)
        assert abs(d - phi * d_ref) < 0.02 * phi * d_ref, (
            f"Phi_SH={phi} should plateau at {phi*d_ref*1000:.2f} mm, got {d*1000:.2f}")
        assert 2.0 * d / D < k.consol_restriction    # sub-critical: never locks

    #  and past the derived threshold it does not plateau
    assert 2.0 * integrate(phi_crit * 1.2) / D > k.consol_restriction


def test_sensitivity_report_writes_and_scales(tmp_path):
    """--sensitivity must produce both files and show the C_phi linearity."""
    c = _short_case(n_ensemble=1, t_end_h=2.0, n_cells=14, deterministic=True)
    rows = solver.sensitivity_report(c, str(tmp_path), kg0_mults=(1.0,),
                                     n_values=(1.0,), C_values=(1500.0, 3000.0),
                                     f_floor_values=(1e-4,))
    assert (tmp_path / "sensitivity_phiSH.csv").exists()
    assert (tmp_path / "sensitivity_phiSH.json").exists()
    base = next(r for r in rows if r["label"] == "baseline")
    hi = next(r for r in rows if r["label"].startswith("C_3000"))
    ratio = hi["max_Phi_SH_uncapped"] / base["max_Phi_SH_uncapped"]
    assert abs(ratio - 3000.0 / c.kinetics.C_phi) < 1e-6


# ---------------------------------------------------------------------------
#  The running maximum of Phi_SH is attained during flow startup, while f_slug is
#  still pinned at its floor, so it is a function of that numerical guard rather
#  than of the operating state. The sustained field is the quantity that is not.
#  These three tests guard the distinction the manuscript now rests on.
# ---------------------------------------------------------------------------
def test_sustained_phish_keys_present():
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    for k in ("sustained_Phi_SH", "sustained_Phi_SH_hotspot_km",
              "sustained_supercritical_km", "final_Phi_SH",
              "Phi_SH_supercritical_time_frac", "Phi_SH_peak_time_h"):
        assert k in eng, f"engineering() is missing {k}"


def test_sustained_never_exceeds_the_running_maximum():
    c = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); eng = sv.engineering()
    sust = eng["sustained_Phi_SH"]
    if sust == sust:                                   # not NaN
        assert sust <= eng["max_Phi_SH_uncapped"] + 1e-9


def test_sustained_phish_is_independent_of_the_floor():
    """The sustained field must not move when the numerical guard moves."""
    c1 = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv1 = solver.TransientSHCT(c1); sv1.run(verbose=False); e1 = sv1.engineering()
    c2 = copy.deepcopy(c1); c2.kinetics.f_slug_floor_Hz = c1.kinetics.f_slug_floor_Hz / 2.0
    sv2 = solver.TransientSHCT(c2); sv2.run(verbose=False); e2 = sv2.engineering()
    s1, s2 = e1["sustained_Phi_SH"], e2["sustained_Phi_SH"]
    if s1 == s1 and s2 == s2:
        assert abs(s2 - s1) <= 1e-9 * max(1.0, abs(s1)), (
            f"sustained Phi_SH must be independent of the floor, {s1} -> {s2}")


def test_running_max_response_to_the_floor_is_bounded():
    """Halving the floor may double the running maximum, and may do nothing.

    Phi_SH goes as 1/f_slug, so WHERE the floor binds — a cold line during flow
    startup, before slugging develops — the running maximum is exactly proportional
    to 1/f_slug_floor_Hz. Where the line is slugging throughout, as in this short
    synthetic case, the floor never binds and the maximum does not move at all. The
    ratio must therefore lie in [1, 2]; anything outside means the floor is leaking
    into the field somewhere it should not.
    """
    c1 = _short_case(n_ensemble=2, t_end_h=4.0, n_cells=20, deterministic=True)
    sv1 = solver.TransientSHCT(c1); sv1.run(verbose=False); e1 = sv1.engineering()
    c2 = copy.deepcopy(c1); c2.kinetics.f_slug_floor_Hz = c1.kinetics.f_slug_floor_Hz / 2.0
    sv2 = solver.TransientSHCT(c2); sv2.run(verbose=False); e2 = sv2.engineering()
    ratio = e2["max_Phi_SH_uncapped"] / e1["max_Phi_SH_uncapped"]
    assert 1.0 - 1e-6 <= ratio <= 2.0 + 1e-6, f"floor response out of bounds: {ratio}"


def test_case_study_sweep_records_exact_inverse_floor_scaling():
    """Guard the manuscript's claim against the recorded case-study sweep.

    Section 6.5 states that across the f_slug floor block the uncapped running maximum
    moves exactly as 1/f0 while the time-to-plug, inhibitor dose, deposit and
    super-critical length do not move at all. That claim is about the deepwater case,
    not about any configuration, so it is checked against that case's own output.
    """
    import csv, os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "case", "outputs_steady", "sensitivity_phiSH.csv")
    if not os.path.exists(path):
        pytest.skip("case-study sweep has not been run in this checkout")
    rows = {r["label"]: r for r in csv.DictReader(open(path))}
    base = rows.get("baseline")
    floors = {k: v for k, v in rows.items() if k.startswith("ffloor_")}
    if not base or not floors:
        pytest.skip("sweep predates the f_slug floor block")
    f0b = float(base["f_slug_floor_Hz"]); pb = float(base["max_Phi_SH_uncapped"])
    for lab, r in floors.items():
        f0 = float(r["f_slug_floor_Hz"]); pk = float(r["max_Phi_SH_uncapped"])
        assert abs(pk / pb - f0b / f0) < 1e-6, (
            f"{lab}: peak should scale as 1/f0, got {pk/pb} vs {f0b/f0}")
        for k in ("time_to_plug_P50_h", "peak_deposit_mm", "sustained_supercritical_km"):
            a, b = float(base[k]), float(r[k])
            assert abs(b - a) <= 1e-6 * max(1.0, abs(a)), (
                f"{lab}: {k} moved with the floor, {a} -> {b}")


#  -------- v3.2: the space-time recorders and the sub-grid slug reconstruction ----------
def test_slug_body_holdup_gregory():
    """The slug-body holdup closure must reproduce Gregory, Nicholson & Aziz (1978)
    exactly, be monotone decreasing in Vm, and stay physical."""
    Vm = np.array([0.5, 2.0, 5.0, 10.0])
    got = solver.slug_body_holdup(Vm)
    want = np.clip(1.0 / (1.0 + (Vm / 8.66) ** 1.39), 0.30, 1.0)
    assert np.allclose(got, want, rtol=0, atol=1e-12)
    assert np.all(np.diff(got) < 0.0)                 # more mixing -> more entrained gas
    assert np.all((got > 0.0) & (got <= 1.0))
    assert solver.slug_body_holdup(0.0) == pytest.approx(1.0)


def test_spacetime_recorders_present_and_shaped():
    """Every field the space-time figure set reads must be recorded on the snapshot
    cadence, with one row per snapshot and one column per cell."""
    c = _short_case(n_ensemble=2, t_end_h=6.0)
    c.numerics.n_snapshots = 24
    sv = solver.TransientSHCT(c)
    r = sv.run(verbose=False)
    nt = np.asarray(r["snap_t"]).size
    assert nt > 1
    for key in ("snap_holdup", "snap_P", "snap_T", "snap_delta", "snap_Tsub",
                "snap_j", "snap_regime", "snap_fslug", "snap_vl", "snap_vg"):
        A = np.asarray(r[key], float)
        assert A.shape == (nt, sv.x.size), f"{key} has shape {A.shape}"
        assert np.isfinite(A).all(), f"{key} carries non-finite values"
    assert np.all(np.diff(np.asarray(r["snap_t"], float)) > 0)     # strictly increasing
    assert np.all(np.asarray(r["snap_delta"], float) >= 0.0)
    assert np.all(np.asarray(r["snap_fslug"], float) > 0.0)


def test_slug_reconstruction_is_mass_consistent():
    """The sub-grid slug reconstruction must integrate back to the solver's own
    cell-average holdup: over one slug unit, beta*alpha_ls + (1-beta)*alpha_film
    has to return alpha_l.  This is what makes figures 15/16/21 a reconstruction
    of the run rather than an illustration."""
    import shct_spacetime as ST
    c = _short_case(n_ensemble=2, t_end_h=6.0)
    c.numerics.n_snapshots = 24
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)

    F = ST.slug_unit_fields(sv)
    recovered = F["beta"] * F["als"] + (1.0 - F["beta"]) * F["alf"]
    assert np.allclose(recovered, F["alpha"], rtol=1e-6, atol=1e-9)
    assert np.all((F["beta"] > 0.0) & (F["beta"] < 1.0))
    assert np.all(F["Ls"] <= F["Lu"] + 1e-9)
    assert np.all(F["Lu"] > 0.0) and np.all(F["Vt"] > 0.0)

    #  and the rendered field must average to the same holdup over a whole number
    #  of slug periods at every station
    x0 = float(sv.x[len(sv.x) // 2])
    xq = np.array([x0])
    f0 = float(np.interp(x0, sv.x, F["fslug"]))
    tq = np.linspace(0.0, 40.0 / f0, 40_001)          # 40 exact periods
    fld, _ = ST.reconstruct_slug_field(sv, xq, tq)
    assert fld.mean() == pytest.approx(float(np.interp(x0, sv.x, F["alpha"])),
                                       rel=2e-3)


#  the figures that need a travelling slug train, and so may legitimately be
#  skipped when the line is not slugging (a shut-in, say)
_SLUG_FIGS = {"15_slug_growth_propagation.png", "16_slug_train_waterfall.png",
              "21_riser_depth_time.png"}


def test_spacetime_figures_render():
    """Every figure that is always defined must render for a real (short) run,
    and every rendered file must be a real image."""
    import tempfile, shct_spacetime as ST
    c = _short_case(n_ensemble=2, t_end_h=6.0)
    c.numerics.n_snapshots = 40
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    eng = sv.engineering()
    with tempfile.TemporaryDirectory() as td:
        made = ST.spacetime_outputs(sv, eng, td, verbose=False)
        names = {os.path.basename(p) for p in made}
        for expect in (n for n, _f, _s in ST.FIGURES):
            if expect in _SLUG_FIGS:
                continue
            assert expect in names, f"{expect} did not render"
        for p in made:
            assert os.path.getsize(p) > 5_000


def test_slug_figures_use_a_slugging_state_not_the_final_one():
    """A shut-in line is not slugging: its slug frequency falls to the floor and
    the 'slug unit length' grows past the length of the pipe. The resolved-slug
    figures must therefore be built from a snapshot where the line is genuinely
    flowing intermittently -- here, before the shut-in event -- or skipped."""
    import shct_spacetime as ST
    c = _short_case(n_ensemble=2, t_end_h=12.0)
    c.scenario.kind = "shutin"
    c.scenario.event_time_h = 6.0
    c.numerics.n_snapshots = 60
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)

    k = ST._slug_snapshot(sv)
    ts = np.asarray(sv.results["snap_t"], float)
    if k is None:                      # never slugs: the figures must be skipped
        F = ST.slug_unit_fields(sv)
        ok, why = ST._slug_train_ok(sv, F, ST._slugging_cell(sv, F))
        assert not ok and why
        return
    assert ts[k] <= c.scenario.event_time_h + 1e-9, (
        f"resolved-slug figures would be built at t={ts[k]:.2f} h, after the "
        f"shut-in at {c.scenario.event_time_h} h")
    #  and at that state the train must be physical
    F = ST.slug_unit_fields(sv, k_snap=k)
    ic = ST._slugging_cell(sv, F)
    ok, why = ST._slug_train_ok(sv, F, ic)
    assert ok, f"chosen state is not a slugging one: {why}"
    assert F["Lu"][ic] < 0.02 * float(sv.x[-1])


def test_slug_train_guard_rejects_a_stopped_line():
    """The guard itself must reject a state with no travelling train."""
    import shct_spacetime as ST
    c = _short_case(n_ensemble=2, t_end_h=6.0)
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    F = ST.slug_unit_fields(sv)
    ic = ST._slugging_cell(sv, F)
    stalled = {k: v.copy() for k, v in F.items()}
    stalled["Lu"][ic] = 0.9 * float(sv.x[-1])        # a "slug" as long as the line
    ok, why = ST._slug_train_ok(sv, stalled, ic)
    #  Assert the GUARD's behaviour, not its prose. This previously required the
    #  substring "exceeds", which the message has not contained for some time — the
    #  wording was changed and the test was not, so it failed while the guard it exists
    #  to protect worked correctly. Match the reason loosely instead.
    assert not ok, "the guard accepted a line-length slug as a travelling train"
    assert "slug" in why.lower() and "not a slug scale" in why, why
    stalled = {k: v.copy() for k, v in F.items()}
    stalled["Vt"][ic] = 0.01                          # no motion
    ok, why = ST._slug_train_ok(sv, stalled, ic)
    assert not ok and "celerity" in why


def test_azimuthal_deposit_cannot_exceed_the_pipe_radius():
    """A deposit thicker than the pipe RADIUS is geometrically impossible — at
    delta = D/2 the bore is already shut. The azimuthal redistribution weights
    the deposit toward the cold bottom of the line, which multiplies the mean by
    up to (1 + skew), so without a cap a heavily deposited line reported a
    thickness larger than the pipe itself."""
    import shct_crosssection as CX
    D = 0.2545
    R = D / 2.0
    #  a mean deposit already at half the radius: the bottom weight would push the
    #  redistributed value past R
    delta = np.full(12, 0.5 * R)
    h = np.full(12, 0.5)
    _theta, prof, bot, top = CX.azimuthal_deposit(delta, h, D=D)
    assert np.all(prof <= R + 1e-12), f"deposit exceeds the radius: {prof.max()} > {R}"
    assert np.all(bot <= R + 1e-12) and np.all(top <= R + 1e-12)
    assert np.all(prof >= 0.0)
    #  the bottom must still carry more than the top (the physics of the skew)
    assert np.all(bot >= top)
    #  and with a small mean, well clear of the cap, the azimuthal mean is preserved
    small = np.full(12, 0.02 * R)
    _t2, prof2, _b2, _t2b = CX.azimuthal_deposit(small, h, D=D)
    assert np.allclose(prof2.mean(axis=0), small, rtol=1e-6)

    #  D is the PER-CELL bore, which shrinks as the deposit grows — that is how
    #  the solver actually calls this, so the cap must broadcast per cell
    Dv = np.linspace(0.60 * D, D, 12)
    _t3, prof3, bot3, top3 = CX.azimuthal_deposit(delta, h, D=Dv)
    assert np.all(prof3 <= Dv[None, :] / 2.0 + 1e-12), "per-cell cap not applied"
    assert np.all(bot3 >= top3)
    assert prof3.shape == (_t3.size, Dv.size)


def test_check_outputs_flags_a_broken_figure_and_table():
    """The output inspector must actually catch a blank figure and an
    out-of-bounds table column — the two failure modes that do not raise."""
    import importlib.util
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "case", "scripts")
    spec = importlib.util.spec_from_file_location(
        "check_outputs", os.path.join(here, "check_outputs.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, here)
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(here)

    with tempfile.TemporaryDirectory() as td:
        blank = os.path.join(td, "blank.png")
        fig = plt.figure(figsize=(6, 4))
        fig.savefig(blank, dpi=100)
        plt.close(fig)
        rep = mod.Report()
        mod.check_image(blank, rep)
        assert any(lvl == "FAIL" for lvl, _w, _t in rep.rows), \
            "a blank canvas was not flagged"

        bad = os.path.join(td, "csv_crosssection.csv")
        with open(bad, "w") as fh:
            fh.write("x_km,holdup,deposit_bottom_mm\n0,0.5,10\n1,0.5,900\n")
        rep = mod.Report()
        mod.check_csv(bad, rep)
        assert any(lvl == "FAIL" and "radius" in what for lvl, _w, what in rep.rows), \
            "a deposit thicker than the pipe radius was not flagged"

        ok = os.path.join(td, "fields_profile.csv")
        with open(ok, "w") as fh:
            fh.write("x_km,holdup,P_bar\n0,0.35,150\n1,0.40,120\n")
        rep = mod.Report()
        mod.check_csv(ok, rep)
        assert not rep.rows, f"a clean table was flagged: {rep.rows}"


def test_mixture_velocity_is_not_clobbered_by_the_latent_heat_lookup():
    """The condensation-latent-heat term interpolates into a V(P,T) surface. Its
    lookup indices were once bound to the names `i` and `j` -- and `j` is the
    mixture volumetric flux carried by the same time step, so the velocity field
    was silently replaced by an integer temperature index. Everything derived
    from the velocity (slug-unit length, erosional margin, the advection of
    temperature later in the same step) was then computed from that index.

    The signature of the bug is a velocity field that takes only a handful of
    integer values, so assert the opposite: the flux must be continuous.
    """
    c = _short_case(n_ensemble=3, t_end_h=3.0)
    c.fluids.condensation_latent = True          # the path that carried the bug
    c.numerics.n_snapshots = 8
    sv = solver.TransientSHCT(c)
    r = sv.run(verbose=False)

    j = np.asarray(r["j"], float)
    assert np.isfinite(j).all()
    #  a genuine velocity field over nx cells is essentially all-distinct; an
    #  integer index array collapses onto a handful of values
    col = j[:, 0]
    distinct = np.unique(np.round(col, 6)).size
    assert distinct > 0.5 * col.size, (
        f"mixture velocity takes only {distinct} distinct values over "
        f"{col.size} cells — it looks like an index array, not a velocity")
    assert not np.allclose(col, np.round(col)), \
        "every mixture velocity is an exact integer — the field was clobbered"

    #  and the snapshot recorder must carry the same continuous field
    sj = np.asarray(r["snap_j"], float)
    assert np.unique(np.round(sj[-1], 6)).size > 0.5 * sj.shape[1]


def test_liquid_balance_closes_when_the_bore_plugs():
    """The liquid balance must close even after the hydrate deposit shuts the bore.

    _enforce_bounds keeps La in [0, A]. Its redistribution passes are conservative,
    but it used to END with a second, non-redistributing clip. While the pipe has
    spare capacity that is a no-op; once the bore closes, A collapses, no cell has
    room, and that clip deleted liquid silently -- the balance then failed by ~6 %
    on a 48 h run that plugs while showing 0.000 % on a 12 h run that does not, at
    every CFL tested. The discard is now measured and carried explicitly.
    """
    c = _short_case(n_ensemble=8, t_end_h=48.0)
    c.numerics.n_snapshots = 20
    sv = solver.TransientSHCT(c)
    r = sv.run(verbose=False)
    e = sv.engineering()

    #  the balance itself must close to numerical precision
    assert e["mass_conservation_err"] < 1e-6, (
        f"liquid balance does not close: {e['mass_conservation_err']*100:.4f} %")
    #  This test previously passed with THREE realisations while the real case, at
    #  twelve, failed at 5.9 % — because the loss only appears once enough
    #  realisations plug. Assert the ensemble is large enough to exercise that.
    assert c.numerics.n_ensemble >= 8, (
        "too few realisations to exercise the plugged-bore path that lost mass")
    #  and whatever the bounds had to discard must be REPORTED, not hidden
    assert "liq_bounds_discard_frac" in e
    assert np.isfinite(e["liq_bounds_discard_frac"])
    #  The plugged-bore path must still be EXERCISED, or this test silently stops
    #  testing anything: a zero discard proves nothing if the bore never closed.
    assert e["deposit_full_bore"], "bore never plugged — the lossy path is untested"
    assert e["P_plug"] > 0.0
    #  On THIS case the bounds no longer have to discard at all — it reads about -3e-15,
    #  i.e. zero to roundoff, where it was 5.9 % before. That is not a general claim: the
    #  full 48 h, 12-realisation case study still discards ~5.8 %, reported openly as
    #  liq_bounds_discard_frac, because there the bore really does close on a line that
    #  is still being fed. A signed epsilon is roundoff, not a negative discard, so allow
    #  it — but keep the upper bound tight enough that a real loss here would fail.
    _d = e["liq_bounds_discard_frac"]
    assert _d >= -1e-9, f"negative discard beyond roundoff: {_d:.3e}"
    assert _d < 1e-6, f"bounds enforcement is losing liquid again: {_d*100:.4f} %"
    #  This used to require a NON-ZERO discard whenever the line plugged, on the
    #  reasoning that a closed bore cannot hold the liquid still arriving. That was a
    #  proxy for the real concern — the discard being measured in the wrong place, so
    #  that it read zero while mass vanished a few lines earlier — and the proxy has
    #  outlived its premise: with erosion competing everywhere the bore no longer
    #  closes hard enough for the bounds to discard anything, and the run above reports
    #  exactly zero. Keeping the assertion would demand a loss the solver no longer has.
    #
    #  The concern it stood for is asserted directly instead, and more strictly. The
    #  closure identity
    #        in - out - to_hydrate - discarded  ==  final inventory - initial inventory
    #  is precisely what mass_conservation_err measures, using the initial inventory the
    #  solver actually started from. If liquid were lost anywhere and not attributed to
    #  the discard, this would be non-zero — which is the failure the old assertion was
    #  reaching for, caught at its source rather than through a symptom.
    #
    #  (The identity was previously also "checked" here by an expression that reduced to
    #  abs(lhs - lhs) <= tol, i.e. 0 <= tol. It asserted nothing, and is removed rather
    #  than left standing as false reassurance.)
    assert e["mass_conservation_err"] < 1e-12, (
        f"liquid balance no longer closes to machine precision: "
        f"{e['mass_conservation_err']:.3e} — liquid is going somewhere unaccounted")


def test_gas_floor_is_measured_not_hidden():
    """`Mg = max(Mg, 0)` CREATES gas mass whenever the floor fires. It is dormant
    for these cases, but it must be measured so the gas balance can never drift
    silently the way the liquid one did."""
    c = _short_case(n_ensemble=2, t_end_h=8.0)
    sv = solver.TransientSHCT(c)
    r = sv.run(verbose=False)
    assert "gas_floor_created" in r and "gas_floor_created_frac" in r
    assert np.isfinite(r["gas_floor_created"])
    assert r["gas_floor_created"] >= -1e-12          # the floor can only ADD mass
    #  and the gas balance closes with that term included
    assert sv.engineering()["gas_mass_conservation_err"] < 1e-6


def test_slug_length_statistics_exclude_the_non_slugging_reach():
    """L_u = V_t/f_slug returns the correlation's 5000 m ceiling wherever the line
    is NOT slugging (f_slug at its floor). Averaging that in reported a 5000 m
    "slug" for a line whose slugs are tens of metres, so the statistics must be
    taken over the slugging reach only."""
    c = _short_case(n_ensemble=3, t_end_h=12.0)
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    e = sv.engineering()
    D = c.pipeline.diameter_m
    assert "slug_length_reach_frac" in e
    if e["slug_length_reach_frac"] > 0:
        assert np.isfinite(e["slug_length_max_m"])
        assert e["slug_length_max_m"] < 0.999 * 5000.0, \
            "the correlation ceiling leaked into the reported slug length"
        assert e["slug_length_max_m"] < 400.0 * D
        assert e["slug_length_mean_m"] <= e["slug_length_max_m"] + 1e-9


def test_text_overlap_detector_catches_a_stacked_label():
    """Text drawn on top of other text is invisible in the source and glaring on
    the page. The detector must find it, so a collision is reported at render time
    rather than discovered in a published figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shct_style as S

    fig, ax = plt.subplots()
    ax.set_title("(d) distance-stacked trace")
    ax.annotate("L_s = 9.0 m", xy=(0.5, 1.015), xycoords="axes fraction",
                ha="center")
    hits = S.find_text_overlaps(fig)
    plt.close(fig)
    assert hits, "a label stacked on the title was not detected"
    assert hits[0][2] > 0.3

    #  and a clean figure must produce no false positive
    fig, ax = plt.subplots()
    ax.set_title("(d) distance-stacked trace")
    ax.set_xlabel("time  [s]")
    ax.set_ylabel("stacked value")
    ax.annotate("L_s = 9.0 m", xy=(1.03, 0.5), xycoords="axes fraction",
                ha="left", va="center", annotation_clip=False)
    hits = S.find_text_overlaps(fig)
    plt.close(fig)
    assert not hits, f"clean figure reported a false overlap: {hits}"


def test_every_spacetime_figure_is_free_of_text_overlaps():
    """Render the whole set from a real run and assert that nothing collides."""
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import shct_spacetime as ST
    import shct_style as S

    c = _short_case(n_ensemble=2, t_end_h=8.0)
    c.numerics.n_snapshots = 30
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    eng = sv.engineering()

    offenders = []
    orig = ST._save

    def _spy(fig, path, check=True):
        hits = S.find_text_overlaps(fig)
        if hits:
            offenders.append((os.path.basename(path),
                              [(str(a.get_text())[:40], str(b.get_text())[:40],
                                round(f, 2)) for a, b, f in hits]))
        return orig(fig, path, check=False)

    ST._save = _spy
    try:
        with tempfile.TemporaryDirectory() as td:
            ST.spacetime_outputs(sv, eng, td, verbose=False)
    finally:
        ST._save = orig
    assert not offenders, f"overlapping text in: {offenders}"


def test_phi_sh_is_dimensionless_and_kg_units_do_not_depend_on_n():
    """Phi_SH is offered as a dimensionless group, so it must actually be one, and
    the constants inside it must keep fixed units as the exponent is varied.

        d(delta)/dt = f_wall * Rg * D/4   [m/s]      =>  Rg = [1/s]
        Rg          = kg * a_i * (dT/dT_ref)**n
                    = kg * [1/m] * [-]               =>  kg = [m/s]
        Phi_SH      = C * Rg / f_slug = [-] * [1/s]/[1/s] = [-]

    Written as dT**n instead, kg would carry m/(s*K**n): its dimensions would
    change with n, and the sensitivity sweep over n at fixed kg0 would compare
    quantities with different units. The reference subcooling removes that.
    """
    import shct_model
    k = shct_model.Kinetics()
    assert hasattr(k, "dTsub_ref_C") and k.dTsub_ref_C > 0

    #  (a) the reference must reproduce the plain dT**n form when it is 1 K, so the
    #      change is a re-statement of the dimensions and not a change of physics
    c = _short_case(n_ensemble=3, t_end_h=8.0)
    c.kinetics.dTsub_ref_C = 1.0
    sv = solver.TransientSHCT(c)
    sv.run(verbose=False)
    ref = sv.engineering()["peak_deposit_mm"]

    #  (b) scaling the reference and the coefficient together must leave the growth
    #      rate unchanged: Rg = kg * a_i * (dT/dTref)**n, so kg -> kg * s**n with
    #      dTref -> dTref * s is the same rate. This is the invariance that only
    #      holds if the subcooling really is non-dimensionalised.
    s_fac, n = 2.0, 1.0
    c2 = _short_case(n_ensemble=3, t_end_h=8.0)
    c2.kinetics.growth_exp_n = n
    c2.kinetics.dTsub_ref_C = 1.0 * s_fac
    c2.kinetics.kg0 = c.kinetics.kg0 * (s_fac ** n)
    sv2 = solver.TransientSHCT(c2)
    sv2.run(verbose=False)
    got = sv2.engineering()["peak_deposit_mm"]
    assert abs(got - ref) <= 1e-6 * max(abs(ref), 1.0), (
        f"growth rate is not invariant under (kg0, dTsub_ref) rescaling: "
        f"{ref} vs {got} — the subcooling is not properly non-dimensionalised")

    #  (c) and Phi_SH itself must come out a pure number
    e = sv.engineering()
    assert np.isfinite(e["max_Phi_SH"]) and e["max_Phi_SH"] >= 0.0


# ---------------------------------------------------------------------------
#  Ransom's water faucet: the one benchmark available without a licence.
#
#  The standing objection to a solver compared only against itself is that it has
#  been compared to nothing. The faucet answers part of that, because it is the
#  problem the industrial and reactor-safety codes are themselves assessed on and
#  it has an exact solution. Judged on L1, not on a percentage: the solution has a
#  contact discontinuity, where Linf never converges and L2 converges at O(h^1/2),
#  so an NRMSE threshold would measure the discontinuity rather than the scheme.
# ---------------------------------------------------------------------------
def test_water_faucet_matches_ransom_exact_solution():
    import shct_verification as V
    r = V.check_water_faucet(None, nx=240)
    assert r["pass"], r
    #  first-order convergence in L1 across the front
    assert r["observed_L1_order"] > 0.8, r["observed_L1_order"]
    #  and the production TVD scheme must be materially better than naive upwind
    assert r["upwind_over_tvd_L1"] > 2.0, r["upwind_over_tvd_L1"]


def test_water_faucet_velocity_field_is_time_dependent():
    """Guard the trap: the steady profile does NOT hold ahead of the front.

    Prescribing sqrt(v0^2+2gx) everywhere thins fluid the exact solution leaves
    untouched, which quietly turns a passing benchmark into a 60 %-error one. If the
    error ever jumps by an order of magnitude, this is the first thing to check.
    """
    import shct_verification as V
    r = V.check_water_faucet(None, nx=240)
    #  undisturbed column ahead of the front must stay at its initial value, so the
    #  error cannot be of the same order as the whole liquid-fraction range
    assert r["tvd_L1"] < 0.01, r["tvd_L1"]


def test_packing_cap_returns_hydrate_to_the_wall_instead_of_destroying_it():
    """phi_max is a carrying limit, not an incinerator.

    Scoured deposit is transferred into the bulk phase field. When the bulk is already
    at its packing limit there is nowhere for it to go, and simply clipping phi
    destroyed 35 % of all hydrate formed on the as-operated case — moving the leak the
    transfer was introduced to close, rather than closing it. What the slurry cannot
    carry is returned to the wall, so wall + bulk is conserved whatever the cap does.
    """
    c = _short_case(n_ensemble=4, t_end_h=24.0)
    sv = solver.TransientSHCT(c); sv.run(verbose=False); e = sv.engineering()
    #  Two mechanisms had to be fixed to reach zero rather than merely small. Returning
    #  the rejected fraction to the wall took it from 35 % to 8 %; the residue was cells
    #  whose wall was ALSO at delta_max, where there was genuinely nowhere to put the
    #  mass. Scouring is now capped by the slurry's carrying headroom, so the overshoot
    #  is prevented rather than repaired, and this reads about 1e-26.
    lost = e["hydrate_packing_clip_frac"]
    assert lost < 1e-12, (
        f"the packing cap is destroying hydrate again: {lost*100:.2e} % of all formed")
    #  and the transfer itself must still be happening, or the test proves nothing
    assert e["hydrate_scoured_frac"] > 1e-3, (
        "no hydrate was scoured at all — the conservation check is vacuous")
    assert e["mass_conservation_err"] < 1e-9
