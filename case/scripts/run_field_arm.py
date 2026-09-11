#!/usr/bin/env python3
"""Run one ARM of the Roberts real-field spool case and record it as JSON.

Arm "asfound"    -- case 1, the rigid spool as it actually was: no inhibitor. It PLUGGED.
Arm "inhibited"  -- the same spool with the methanol treatment that the operating
                    procedure called for and that was not performed.
Arm "case2"      -- case 2, a DIFFERENT line on the same development: the 2.2 km 10 in
                    flexible flowline DC1 -> FPSO, shut in for 11 h, five hours beyond
                    its no-touch time. It did NOT block. A real confirmed negative, and
                    the harder of the two tests -- a model that predicts plugging
                    everywhere passes case 1 and fails here.

The two arms exist because a single positive proves very little. A model that
predicts plugging everywhere passes the plugged case and is worthless. The
inhibited arm is the control that catches exactly that failure: it is the same
geometry, the same fluid and the same shut-in, changed only where the real
difference between blocking and not blocking was said to lie.

It is a CONTROL, not a second field observation -- no one measured this spool
with methanol in it. Do not report it as a validated negative.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _ROOT)
import solver  # noqa: E402

REF = os.path.join(_ROOT, "validation", "data", "field_olga_csmhyk_roberts.json")
OUT = os.path.join(_ROOT, "validation", "data")


def field_assay():
    """Lump the published field assay onto this solver's component set."""
    with open(REF) as fh:
        raw = json.load(fh)["field_crude_assay"]["composition_mol_pct"]
    keys = {"N2": ["N2"], "CO2": ["CO2"], "C1": ["C1"], "C2": ["C2"], "C3": ["C3"],
            "iC4": ["iC4"], "nC4": ["nC4"], "iC5": ["iC5", "neoC5"], "nC5": ["nC5"],
            "C6": ["C6", "methylcyclopentane", "benzene", "cyclohexane"]}
    comp = {k: sum(raw.get(s, 0.0) for s in src) for k, src in keys.items()}
    named = {s for src in keys.values() for s in src}
    comp["C7+"] = sum(v for k, v in raw.items() if k not in named and k != "H2S")
    tot = sum(comp.values())
    return {k: v / tot for k, v in comp.items()}


#  Table 6, the 10 in FLEXIBLE RISER wall. Note what it does NOT contain: any insulation
#  layer at all. The flowline (Table 7) carries 11 mm of Insulation MO1 at k = 0.0605; the
#  riser carries none. So the last ~180 m of segment 1a cools far faster than the rest, and
#  a single-wall-stack solver run with the flowline wall is OPTIMISTIC exactly at the riser
#  base, which is one of the two places the reference model puts hydrate.
#
#  This stack is the other bound: run the whole line as if it were riser. The truth is
#  between the two.
#
#  WHAT THE BOUND ACTUALLY TESTS, stated precisely, because the loose version of this claim
#  is wrong. Running both bounds and finding P_plug = 0 either way shows the BIT is robust.
#  It does NOT show the wall simplification is harmless, and measured across the two bounds
#  it plainly is not:
#
#      max_subcooling_C     7.434 -> 7.465    x1.00   robust
#      P_plug               0     -> 0        --      robust
#      max_Phi_SH           385.7 -> 486.9    x1.26
#      peak_deposit_mm      17.89 -> 43.56    x2.43   NOT robust
#      restriction          0.141 -> 0.343            (trip 0.85)
#
#  The deposit MAGNITUDE moves by a factor of 2.4, and magnitude is one of the three things
#  the real-line comparison scores -- against a reference bound of <1 vol% hydrate at the
#  riser base, which is exactly where the two walls differ. The bit survives only because
#  BOTH bounds sit far below the trip; on a line closer to it the wall would decide the
#  answer outright. That is why the segmented wall is the model and these two are bounds.
_RISER_WALL = [[0.0060, 1161.65, 4800.0 * 502.32],
               [0.0096, 0.1849, 1700.0 * 1500.0],
               [0.0080, 0.9300, 7850.0 * 502.32],
               [0.0075, 0.9300, 7850.0 * 502.32],
               [0.0015, 0.3350, 1040.0 * 2428.0],
               [0.0050, 0.9300, 7850.0 * 502.32],
               [0.0015, 0.3350, 1040.0 * 2428.0],
               [0.0050, 0.9300, 7850.0 * 502.32],
               [0.00155, 1.1628, 571.0 * 2511.60],
               [0.0093, 0.2000, 1010.0 * 2220.0],
               [0.0097, 0.4070, 950.0 * 2500.0]]


_FLOWLINE_WALL = [[0.0060, 1161.65, 4800.0 * 502.32],
                  [0.0070, 0.1849, 1700.0 * 1500.0],
                  [0.0100, 0.1849, 1700.0 * 1500.0],
                  [0.0050, 0.9300, 7850.0 * 502.32],
                  [0.0050, 0.9300, 7850.0 * 502.32],
                  [0.00145, 1.1628, 571.0 * 2511.60],
                  [0.0087, 0.4070, 955.0 * 2302.0],
                  [0.0110, 0.0605, 730.0 * 1214.0],
                  [0.00105, 1.1628, 765.0 * 2511.60],
                  [0.0114, 0.4070, 950.0 * 2500.0]]

#  Where the flowline ends and the riser begins, from the Figure 19 flowpath: the seabed
#  run holds to about 2020 m and the lazy-wave riser carries the last ~180 m.
_RISER_START_M = 2020.0


def build_case2(comp, n_ensemble, t_end_h, seed=7, wall="segmented"):
    """Case 2: flowline segment 1a, DC1 to the FPSO, at the conditions prevailing when the
    plant tripped. Every number here is from the thesis; see field_olga_csmhyk_roberts.json
    field_cases[1] for where each one is read from."""
    c = solver.Case()
    c.name = "Roberts case 2 - flowline 1a, DC1 to FPSO (real, did not block)"
    #  --- geometry: 2.2 km of 10 in flexible, seabed run then a lazy-wave riser
    c.pipeline.length_m = 2200.0
    c.pipeline.diameter_m = 0.254
    c.pipeline.n_cells = 60
    #  Figure 19 flowpath, digitised: a dip over the first ~60 m, seabed at about +8 m to
    #  ~2020 m, then the riser -- up to ~70, down to ~32 (the lazy-wave sag), then ~150 at
    #  the FPSO. Interpolated onto the cell centres.
    _x = [0.0, 30.0, 60.0, 90.0, 2020.0, 2100.0, 2130.0, 2180.0, 2200.0]
    _z = [0.0, -4.0, -2.0, 8.0, 8.0, 70.0, 32.0, 140.0, 150.0]
    _xc = [(i + 0.5) * c.pipeline.length_m / c.pipeline.n_cells
           for i in range(c.pipeline.n_cells)]
    c.pipeline.elevation_m = list(np.interp(_xc, _x, _z))
    #  --- wall: Table 7, 10 in original flexible flowline. [thickness, k, rho*cp] per layer.
    #  The last ~180 m is the RISER, Table 6, which carries no insulation layer at all; this
    #  solver takes one wall stack for the whole pipe, so the riser is modelled with the
    #  flowline's insulation and this run is OPTIMISTIC about the riser base -- one of the
    #  two places the thesis puts hydrate. Stated, not corrected.
    #  THE REAL WALL, resolved along the line. Segment 1a is insulated flexible flowline
    #  (Table 7, 11 mm of MO1 at k = 0.0605) for ~2020 m and then flexible RISER (Table 6)
    #  which carries no insulation layer at all. Running one stack for both was a real
    #  limitation of this comparison, not a presentational one: the riser's U comes out
    #  3.6x the flowline's, and the riser base is one of the two places the reference model
    #  puts hydrate. The solver now takes a segmented wall, so the pipe is described as it
    #  is. The two single-stack options remain as BOUNDS.
    if wall == "segmented":
        c.pipeline.wall_layers = [
            {"from_m": 0.0, "to_m": _RISER_START_M, "layers": _FLOWLINE_WALL},
            {"from_m": _RISER_START_M, "to_m": c.pipeline.length_m, "layers": _RISER_WALL},
        ]
    elif wall == "riser":
        c.pipeline.wall_layers = _RISER_WALL
    else:
        c.pipeline.wall_layers = _FLOWLINE_WALL
    #  --- ambient: the same seawater as case 1
    c.operating.T_seabed_C = -1.8
    c.pipeline.h_outer = 500.0
    #  --- rates prior to the trip, Table 11 row "DC1 - 1a": oil 2.93, water 2.50,
    #  gas 1.60 + 0.13 gas lift kg/s at 62 C. Converted to in-situ volumetric here.
    _rho_o, _rho_w = 883.0, 1025.0            # stock-tank oil density from the field assay
    _q_o, _q_w = 2.93 / _rho_o, 2.50 / _rho_w
    c.operating.q_liquid_insitu = _q_o + _q_w
    c.fluids.water_cut = _q_w / (_q_o + _q_w)
    c.fluids.rho_oil, c.fluids.rho_water = _rho_o, _rho_w
    c.operating.P_inlet_bar = 46.0            # DC1a manifold, 4600 kPaa, field = OLGA
    c.operating.T_inlet_C = 62.0
    #  gas: 1.73 kg/s at inlet P,T through the solver's own gas density, not a hand figure
    _rho_g = float(np.ravel(solver.gas_density(np.array([c.operating.P_inlet_bar]),
                                               np.array([c.operating.T_inlet_C]),
                                               c.fluids))[0])
    c.operating.q_gas_insitu_inlet = 1.73 / _rho_g
    c.operating.MEG_wt_inlet = 0.0
    #  the thesis built its hydrate curves on lift gas with PURE water, no salt. Match that
    #  rather than quietly giving this run a salinity suppression the reference did not have.
    c.fluids.salinity_wt = 0.0
    c.fluids.composition = comp
    #  --- the event: plant trip, then stagnant. Real duration 11 h before blowdown; the
    #  thesis extended its own run to 20 h as a sensitivity and saw growth accelerate from
    #  ~17 h, so run to 20 h and read BOTH the real 11 h point and that acceleration.
    c.scenario.kind = "shutin"
    c.scenario.event_time_h = 0.5
    c.numerics.t_end_h = float(t_end_h)
    c.numerics.n_ensemble = int(n_ensemble)
    c.numerics.seed = int(seed)
    return c


def build_case(comp, meg_wt, n_ensemble, t_end_h, seed=7):
    c = solver.Case()
    c.name = "Roberts case 1 - rigid production spool (real, blocked)"
    #  --- geometry: 6 in rigid production spool, ~35 m, tree to manifold
    c.pipeline.length_m = 35.0
    c.pipeline.diameter_m = 0.1241
    c.pipeline.n_cells = 30
    c.pipeline.elevation_m = None
    #  --- wall: duplex 21.95 mm then SPU 55 mm. volumetric heat capacity = rho * cp
    c.pipeline.wall_layers = [[0.02195, 14.4, 7600.0 * 480.0],
                              [0.0550, 0.135, 700.0 * 1700.0]]
    #  --- ambient: seawater -1.8 C, swept by a 0.94 m/s current
    c.operating.T_seabed_C = -1.8
    c.pipeline.h_outer = 500.0
    #  --- rates at shut-in: 600 Sm3/d liquid, 7000 Sm3/d gas
    c.operating.q_liquid_insitu = 600.0 / 86400.0
    c.operating.q_gas_insitu_inlet = 7000.0 / 86400.0
    c.operating.P_inlet_bar = 53.0          # spool open to flowline pressure ~5300 kPa
    c.operating.T_inlet_C = 60.0
    c.operating.MEG_wt_inlet = float(meg_wt)
    c.fluids.composition = comp
    #  --- the event: the well is shut in and the spool cools to seabed
    c.scenario.kind = "shutin"
    c.scenario.event_time_h = 0.5
    #  A 35 m spool is three orders shorter than the 32 km case study, so its CFL
    #  timestep is correspondingly tiny. The physics of interest -- cool down,
    #  subcool, deposit -- is complete well inside 12 h.
    c.numerics.t_end_h = float(t_end_h)
    c.numerics.n_ensemble = int(n_ensemble)
    c.numerics.seed = int(seed)
    return c


def _clean(v):
    """JSON-safe: bools stay bools, NaN becomes null, everything else a float."""
    if isinstance(v, bool):
        return v
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f == f else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=["asfound", "inhibited", "case2",
                                     "case2_riserwall", "case2_flowlinewall"],
                    required=True)
    ap.add_argument("--meg", type=float, default=30.0,
                    help="inhibitor wt%% on the aqueous phase for the inhibited arm")
    ap.add_argument("--ensemble", type=int, default=8)
    ap.add_argument("--cells", type=int, default=None,
                    help="override n_cells (the spool's default 30 gives dx = 1.17 m, and the "
                         "timestep follows dx, so this is the main cost knob)")
    ap.add_argument("--quiet", action="store_true", help="suppress the per-step progress line")
    ap.add_argument("--n", type=float, default=None, dest="growth_n",
                    help="override kinetics.growth_exp_n. The published 1.5-2.5 range is for "
                         "hydrate FILM growth at a gas-liquid interface; this model grows a WALL "
                         "deposit, and no wall exponent has ever been measured. A real line with a "
                         "known outcome can constrain one -- that is what this flag is for.")
    ap.add_argument("--hours", type=float, default=None,
                    help="simulated duration; defaults to 12 h for the spool, 20 h for case 2")
    a = ap.parse_args()

    comp = field_assay()
    hours = a.hours if a.hours is not None else (20.0 if a.arm.startswith("case2") else 12.0)
    meg = a.meg if a.arm == "inhibited" else 0.0
    t0 = time.time()
    _wall = {"case2": "segmented", "case2_riserwall": "riser",
             "case2_flowlinewall": "flowline"}.get(a.arm, "segmented")
    case = (build_case2(comp, a.ensemble, hours, wall=_wall)
            if a.arm.startswith("case2") else build_case(comp, meg, a.ensemble, hours))
    if a.cells:
        case.pipeline.n_cells = int(a.cells)
    if a.growth_n is not None:
        case.kinetics.growth_exp_n = float(a.growth_n)
    #  VERBOSE. These arms were first run with verbose=False and the 35 m spool then sat for
    #  84 minutes with no output at all -- long past the 58 minutes its own CFL scaling
    #  predicts (steps ~ t_end/dx, cost/step ~ cells x realisations, calibrated against the
    #  2.2 km arm's 366 s). With no progress line there was no way to tell a slow run from a
    #  stuck one, or to see sub-stepping and fallbacks piling up. A run that can take an hour
    #  must say where it is.
    sv = solver.TransientSHCT(case)
    sv.run(verbose=not a.quiet)
    e = sv.engineering()
    wall_s = time.time() - t0

    _line = {
        "asfound": "Roberts case 1 -- rigid production spool, 35 m x 124.1 mm ID, shut-in",
        "inhibited": "Roberts case 1 spool, CONTROL arm with the skipped methanol treatment",
        "case2": "Roberts case 2 -- flowline 1a, DC1 to FPSO, 2200 m x 254 mm 10 in flexible, "
                 "shut in 11 h (run to 20 h), five hours beyond the 6 h no-touch time",
        "case2_flowlinewall": "Roberts case 2, BOUNDING arm -- the same line run entirely with "
                              "the INSULATED FLOWLINE wall (Table 7). The warm bound, and what "
                              "this comparison used before the solver could vary the wall",
        "case2_riserwall": "Roberts case 2, COLD BOUND -- the same line run entirely with the "
                           "UNINSULATED riser wall (Table 6). Not the real line. It bounds the "
                           "wall simplification; it does not excuse it. Against the warm bound "
                           "the plug/no-plug bit and the subcooling are unchanged, but the peak "
                           "deposit moves x2.43 (17.89 -> 43.56 mm), so the wall DOES decide the "
                           "magnitude the real-line comparison scores",
    }[a.arm]
    _outcome = {
        "asfound": "PLUGGED -- hydrate blockage, pressure built ~5300 -> 8000 kPa",
        "inhibited": "NOT OBSERVED -- control arm, the treatment was never performed",
        "case2": "DID NOT PLUG -- restarted without incident. Reference model gave <2 % hydrate "
                 "at the drill-centre manifold and <1 % at the riser base at 11 h",
        "case2_flowlinewall": "DID NOT PLUG (the real line). This arm is a WARM BOUND on it, not "
                              "a separate observation",
        "case2_riserwall": "DID NOT PLUG (the real line). This arm is a COLD BOUND on it, not a "
                           "separate observation",
    }[a.arm]
    rec = {
        "arm": a.arm,
        "line": _line,
        "is_real_observation": a.arm in ("asfound", "case2"),
        "observed_plugged": {"asfound": True, "case2": False}.get(a.arm),
        "is_bounding_arm": a.arm in ("case2_riserwall", "case2_flowlinewall"),
        "wall_model": _wall,
        "observed_outcome": _outcome,
        "MEG_wt_inlet": meg,
        "n_ensemble": a.ensemble,
        "n_cells": case.pipeline.n_cells,
        "growth_exp_n": case.kinetics.growth_exp_n,
        "t_end_h": hours,
        "wall_clock_s": round(wall_s, 1),
        "predicted": {k: _clean(e[k])
                      for k in ("P_plug", "max_Phi_SH", "max_subcooling_C",
                                "peak_deposit_mm", "time_to_plug_P50_h",
                                "Phi_SH_supercritical_time_frac",
                                "sustained_supercritical_km", "deposit_full_bore",
                                "MEG_wt_pct", "under_inhibited_km")
                      if k in e},
    }
    #  the tag carries the DURATION too. The real line was stagnant 11 h before
    #  blowdown, so an outcome comparison must be read at 11 h; the 20 h default
    #  exists to see the ~17 h acceleration the reference reports, and scoring
    #  "did not block" against a 20 h deposit compares the model against nine
    #  hours that never happened.
    _tag = a.arm if a.growth_n is None else "%s_n%g_%gh" % (a.arm, a.growth_n, hours)
    path = os.path.join(OUT, "field_roberts_arm_%s.json" % _tag)
    with open(path, "w") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    print(json.dumps(rec, indent=2))
    print("written %s in %.0f s" % (path, wall_s))
    return 0


if __name__ == "__main__":
    sys.exit(main())
