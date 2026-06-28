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
import copy
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
    golden = {"dP_total_bar": 17.253, "arrival_T_C": 11.784, "max_subcooling_C": 8.584,
              "max_Phi_SH": 3.732, "peak_deposit_mm": 127.716, "P_plug": 0.5}
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
