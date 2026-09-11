#!/usr/bin/env python3
"""CALIBRATE the plugging prediction against real, measured outcomes -- and state
precisely how far that calibration reaches.

Until now every constant in the plugging chain was either a literature-typical default or
an assumption. `numerics.plug_restriction_trip` -- the bore fraction at which the solver
declares a plug -- shipped at 0.85 with nothing behind it at all. This script replaces
that with a number bounded by things that actually happened.

HOW A PLUG/NO-PLUG OUTCOME CALIBRATES A TRIP
--------------------------------------------
The solver computes a bore restriction R = 2*delta/D in every cell of every realisation,
and declares a plug when R exceeds the trip. So a real outcome constrains the trip from
ONE side, and which side depends on what the line did:

    a line that DID plug     -> the model must call it plugged  -> trip <= R+
    a line that did NOT plug -> the model must call it clear    -> trip >  R-

One observation gives one inequality and an unbounded interval. TWO observations in
OPPOSITE directions give a bracket: the trip must lie in (R-, R+]. That is the whole
reason case 2 was worth the work -- a second positive would have added nothing.

WHAT THIS IS NOT
----------------
It is not a fit of C_phi, k_ero, k_g0 or growth_exp_n. Two outcomes cannot determine four
constants, and pretending otherwise would be worse than leaving them declared. It is one
parameter, bounded by two real lines, reported with the width of its bracket.
"""
import json
import math
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(_ROOT, "validation", "data")
SHIPPED_TRIP = 0.85


def _load(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def main():
    rule = "=" * 78
    #  USE EACH ARM AT ITS REAL DURATION. case 2's line was stagnant 11 h before blowdown,
    #  so the 20 h default record scores the model against nine hours that never happened
    #  and inflates its restriction -- f = 0.308 at 20 h against 0.002 at 11 h. The bracket
    #  is built from restrictions, so the wrong window moves the bound it sets. Prefer the
    #  11 h record where one exists.
    arms = {a: _load("field_roberts_arm_%s.json" % a)
            for a in ("asfound", "case2", "inhibited")}
    _c2_real = (_load("field_roberts_arm_case2_n1.5_11h.json")
                or _load("field_roberts_arm_case2_n1_11h.json"))
    if _c2_real is not None:
        arms["case2"] = _c2_real
    bores = {"asfound": 0.1241, "case2": 0.254, "inhibited": 0.1241}

    print(rule)
    print(" CALIBRATING THE PLUG TRIP ON REAL OUTCOMES")
    print(rule)

    obs = []
    for arm in ("asfound", "case2"):
        rec = arms.get(arm)
        if rec is None:
            print("  [missing] %s -- run: python case/scripts/run_field_arm.py --arm %s"
                  % (arm, arm))
            continue
        dep = rec["predicted"].get("peak_deposit_mm")
        R = 2.0 * (dep / 1000.0) / bores[arm] if dep else float("nan")
        obs.append({"arm": arm, "plugged": bool(rec["observed_plugged"]),
                    "peak_deposit_mm": dep, "restriction": R,
                    "P_plug_at_shipped_trip": rec["predicted"].get("P_plug")})
        print("  %-9s observed %-8s  model peak deposit %6.2f mm  ->  restriction %.3f"
              "   [%.0f h, %s realisations]"
              % (arm, "PLUG" if rec["observed_plugged"] else "no plug", dep, R,
                 rec.get("t_end_h") or 0, rec.get("n_ensemble") or "?"))

    pos = [o["restriction"] for o in obs if o["plugged"]]
    neg = [o["restriction"] for o in obs if not o["plugged"]]
    out = {"shipped_trip": SHIPPED_TRIP, "observations": obs}

    print()
    print("  Durations: case 2 is scored at its REAL 11 h stagnant period. The spool arm's")
    print("  12 h is a MODELLING CHOICE -- the thesis records that the spool blocked after")
    print("  shut-in but not how long it took, so the upper bound below inherits that choice")
    print("  and is not purely observational.")
    print()
    if pos and neg:
        lo, hi = max(neg), min(pos)
        out["bracket"] = {"must_exceed": lo, "must_not_exceed": hi,
                          "is_consistent": bool(lo < hi)}
        if lo < hi:
            #  a point estimate inside the bracket. GEOMETRIC mean, because the trip is a
            #  ratio and the bracket spans most of an order of magnitude; the arithmetic
            #  mean of 0.14 and 0.92 would sit at 0.53 purely because the upper bound is
            #  clipped by delta_max_frac, which is a numerical cap and not a measurement.
            pt = math.sqrt(lo * hi)
            out["bracket"]["point_estimate_geometric"] = pt
            out["bracket"]["shipped_inside"] = bool(lo < SHIPPED_TRIP <= hi)
            print("  BRACKET FROM REAL DATA:  %.3f  <  trip  <=  %.3f" % (lo, hi))
            print("    below %.3f the line that did NOT block would be called plugged;" % lo)
            print("    above %.3f the line that DID block would be called clear." % hi)
            print("    width %.3f of the unit interval; geometric midpoint %.3f."
                  % (hi - lo, pt))
            print("    SHIPPED 0.85 -> %s"
                  % ("INSIDE the bracket. It was an assumption; it is now an assumption "
                     "that real data does not contradict."
                     if out["bracket"]["shipped_inside"] else
                     "OUTSIDE the bracket. Real data contradicts the shipped value."))
        else:
            print("  NO CONSISTENT TRIP EXISTS.  The line that did not block reaches a")
            print("  restriction of %.3f, at or above the %.3f of the line that did." % (lo, hi))
            print("  No threshold on bore restriction separates these two real outcomes, so")
            print("  the trip is not the right decision variable. Reported, not tuned away.")
    else:
        print("  Cannot bracket: need one real PLUG and one real NO-PLUG. Have %d and %d."
              % (len(pos), len(neg)))

    # ------------------------------- an INDEPENDENT measured bracket on the same variable
    #  The solver computes restr = 2*delta/D. With delta = f*R = f*D/2 that is exactly f,
    #  the deposit thickness as a fraction of the bore RADIUS -- so a measured pipe-volume
    #  occupancy converts into this model's own trip variable with no fudge.
    print()
    print(" A SECOND BRACKET, FROM A LAB LOOP, ON THE SAME VARIABLE")
    LAB = {"plugged_f": 0.117, "no_plug_f": 0.041, "bore_m": 0.0254,
           "source": "DOE DE-FE0031578 (OSTI 1986259): 22 % of the cooled pipe volume "
                     "occupied by the water-hydrate mass on the uncoated line, which trended "
                     "to plugging; 8 % on the coated line, which could never reach plugging "
                     "conditions. An annular deposit at radius-fraction f occupies "
                     "1-(1-f)^2, so 22 % -> f = 0.117 and 8 % -> f = 0.041."}
    print("   measured: PLUGGED at f = %.3f,  did NOT plug at f = %.3f  (bore %.4f m)"
          % (LAB["plugged_f"], LAB["no_plug_f"], LAB["bore_m"]))
    print("   shipped : plug_restriction_trip = %.2f  ->  %.1fx the f at which that loop plugged."
          % (SHIPPED_TRIP, SHIPPED_TRIP / LAB["plugged_f"]))
    print("   Read carefully: the loop's 'plugging conditions' is a PRESSURE-DROP TREND, while")
    print("   this trip declares a bore blockage. A line becomes inoperable long before 85 % of")
    print("   its radius is deposit, so the two are not the same event and 7.3x is not an error")
    print("   factor. What it does say is that this trip is an OPERABILITY-LATE criterion.")
    out["lab_loop_bracket"] = LAB

    #  and the two real sources DISAGREE, which is the point
    _c2 = arms.get("case2")
    if _c2 and _c2["predicted"].get("peak_deposit_mm"):
        f_c2 = 2.0 * (_c2["predicted"]["peak_deposit_mm"] / 1000.0) / bores["case2"]
        #  COMPUTE the relation; do not assert it. This block used to print "ABOVE the lab
        #  loop's plugging fraction" as fixed text, and the number it printed beside those
        #  words was BELOW it. The claim was inherited from an earlier run at f = 0.141 --
        #  the arm that put the INSULATED flowline wall on the bare riser too. Once the wall
        #  was modelled segment by segment the field figure fell to 0.075 and the conflict
        #  went with it, but the sentence stayed. Hard-coded conclusions survive the
        #  evidence that produced them; computed ones cannot.
        _consistent = f_c2 < LAB["plugged_f"]
        _marginal = LAB["no_plug_f"] <= f_c2 < LAB["plugged_f"]
        print()
        print("   DO THE TWO SCALES AGREE?  Computed, not assumed:")
        print("     lab loop (0.0254 m bore)   no plug at f = %.3f, PLUGGED at f = %.3f"
              % (LAB["no_plug_f"], LAB["plugged_f"]))
        print("     Roberts flowline (0.254 m) did NOT plug; this solver puts it at f = %.4f"
              % f_c2)
        if _consistent:
            print("     %.4f < %.3f, so the field observation sits BELOW the lab plugging"
                  % (f_c2, LAB["plugged_f"]))
            print("     fraction and the two are CONSISTENT with one threshold band%s."
                  % (" -- inside its marginal range" if _marginal else ""))
            print("     A factor of ten in bore does not break the band. That is not proof it")
            print("     transfers -- two points cannot establish a scaling -- but the")
            print("     transferability objection is no longer supported by the data.")
        else:
            print("     %.4f >= %.3f, so a single threshold cannot satisfy both: calibrating"
                  % (f_c2, LAB["plugged_f"]))
            print("     the trip to the lab number would call a real, un-plugged field line")
            print("     plugged. The band stays a band and is not fitted into a trip.")
        out["scale_transfer"] = {"lab_plugged_f": LAB["plugged_f"],
                                 "lab_no_plug_f": LAB["no_plug_f"],
                                 "field_no_plug_f": f_c2,
                                 "consistent_with_one_band": bool(_consistent),
                                 "field_inside_marginal_band": bool(_marginal),
                                 "note": "the earlier INCONSISTENT reading (f = 0.141) came "
                                         "from running the bare riser with the insulated "
                                         "flowline wall; the segmented wall removes it"}

    # ---------------------------------------------------------------- magnitude vs OLGA
    rl = None
    p_rl = os.path.join(_ROOT, "case", "outputs_realline", "real_line_score.json")
    if os.path.exists(p_rl):
        with open(p_rl) as fh:
            rl = json.load(fh)
    print()
    print(" MAGNITUDE AGAINST A LICENSED TRANSIENT CODE")
    if rl and rl.get("magnitude_at_11h"):
        m = rl["magnitude_at_11h"]
        print("  Roberts case 2 at %.1f h, OLGA + CSMHyK reference bounds:" % m["snapshot_time_h"])
        for place, key, bound in (("drill-centre manifold", "manifold_phi_volfrac", 0.02),
                                  ("riser base", "riser_base_phi_volfrac", 0.01)):
            v = m.get(key)
            if v is None:
                continue
            ratio = v / bound if bound else float("nan")
            print("    %-22s this solver %.4f vol  vs reference < %.2f  -> %.2fx the bound%s"
                  % (place, v, bound, ratio, "" if ratio <= 1.0 else "  ABOVE"))
        out["magnitude_vs_olga"] = m
        print("  This is a MAGNITUDE check, not a bit: two independent codes on the same")
        print("  real line at the same hour. It constrains the kinetics in a way that")
        print("  plug/no-plug cannot.")
    else:
        print("  [pending] run: python case/scripts/run_real_line_case.py")

    # ------------------------------------- the calibrated prediction this data DOES support
    #  The lab bracket and the field outcomes disagree about a single blockage threshold, so
    #  the trip cannot be fitted. But they disagree about ONE event and agree about another,
    #  and separating the two is what turns this into a calibrated prediction rather than a
    #  bracketed assumption.
    #
    #      plug_restriction_trip = 0.85   BORE BLOCKAGE. Assumed. Bracketed by field
    #                                     outcomes, never measured, operability-LATE.
    #      f_operability = 0.041-0.117    OPERABILITY LOSS. MEASURED, on a real loop, in
    #                                     this model's own variable: 8 % occupancy still
    #                                     transported, 22 % trended to plugging.
    #
    #  A line that reaches f = 0.117 is not blocked; it is losing the pressure-drop battle.
    #  That is the event an operator acts on, it is the one with a measured threshold, and
    #  the model already computes the variable it is expressed in.
    print()
    print(rule)
    print(" THE CALIBRATED PREDICTION THIS DATA SUPPORTS")
    print(rule)
    print(" Two different events, and only one of them has a measured threshold:")
    print("   bore blockage   trip = %.2f of bore radius   ASSUMED, bracketed, operability-late"
          % SHIPPED_TRIP)
    print("   operability     f    = %.3f-%.3f            MEASURED on a real loop"
          % (LAB["no_plug_f"], LAB["plugged_f"]))
    print()
    scen = {}
    for name in ("steady", "shutin", "mitigated"):
        km = os.path.join(_ROOT, "case", "outputs_%s" % name, "key_metrics.json")
        if not os.path.exists(km):
            continue
        with open(km) as fh:
            d = json.load(fh)
        dep = d.get("peak_deposit_mm")
        dia = d.get("diameter_m") or 0.2545
        if not isinstance(dep, (int, float)):
            continue
        f = 2.0 * (dep / 1000.0) / dia
        verdict = ("BLOCKED" if f >= SHIPPED_TRIP else
                   "operability-limited" if f >= LAB["plugged_f"] else
                   "marginal" if f >= LAB["no_plug_f"] else "transports")
        scen[name] = {"peak_deposit_mm": dep, "f": f, "verdict": verdict,
                      "P_plug_at_assumed_trip": d.get("P_plug")}
        print("   %-11s peak deposit %7.3f mm -> f = %.4f   P_plug(assumed) %.2f   -> %s"
              % (name, dep, f, d.get("P_plug", float("nan")), verdict))
    out["calibrated_operability"] = {
        "threshold_measured_f": [LAB["no_plug_f"], LAB["plugged_f"]],
        "assumed_blockage_trip": SHIPPED_TRIP,
        "scenarios": scen,
        "meaning": "f is the deposit as a fraction of bore RADIUS, which is exactly the "
                   "solver's restr = 2*delta/D. The operability band is measured; the "
                   "blockage trip is not."}
    print()
    print(" This is the calibrated part: whether the line stays OPERABLE is now judged")
    print(" against a threshold measured on a real loop, in the model's own variable,")
    print(" instead of against a number nobody measured. Whether it fully BLOCKS is still")
    print(" judged against an assumption, and that assumption is late.")
    print(" The band is from a 0.0254 m loop; transfer to a 0.2545 m line is NOT established")
    print(" -- see the disagreement above -- so it is reported as a band, not a trip.")

    # ---------------------------------------------------------------- the honest ledger
    print()
    print(rule)
    print(" WHAT IS NOW CALIBRATED ON REAL VALUES, AND WHAT IS NOT")
    print(rule)
    print(" CALIBRATED ON REAL MEASUREMENTS")
    print("   hydrate equilibrium curve   Deaton & Frost (1946), 11 points, RMSE 0.165 C")
    #  do not print a calibration that did not happen. This block asserted the trip was
    #  bracketed regardless of whether the runs that bracket it had been made -- the same
    #  class of unconditional claim this audit exists to remove.
    _br = out.get("bracket")
    if _br and _br.get("is_consistent"):
        print("   plug_restriction_trip       bracketed by two real field outcomes: "
              "%.3f < trip <= %.3f" % (_br["must_exceed"], _br["must_not_exceed"]))
    elif _br:
        print("   plug_restriction_trip       NOT calibrated -- the two real outcomes do "
              "not bracket it (see above)")
    else:
        print("   plug_restriction_trip       NOT calibrated -- both bracketing runs have "
              "not been made")
    print(" SCORED ON REAL MEASUREMENTS, RESIDUAL REPORTED, NOT FITTED")
    print("   line holdup                 das Neves (2025), 12 QCV points, void RMSE 0.056,")
    print("                               residual +0.82 correlated with mixture Froude")
    print(" ANCHORED TO PUBLISHED RANGES, NOT FITTED TO THIS LINE")
    print("   k_g0                        Englezos (1987), reconciled to order of magnitude")
    print("   growth_exp_n                inside the published 1.0-2.0 range; n = 1.5 bounded")
    print("   tau_deposit                 Di Lorenzo (2018), 100-200 Pa measured in situ")
    print(" NOT CALIBRATED AT ALL")
    print("   C_phi, k_ero                no measured counterpart exists in open literature")
    print("   Phi_crit = 1.08             a decomposition of the two constants above, so it")
    print("                               inherits their status: derived, falsifiable, unfitted")
    print(rule)

    p = os.path.join(DATA, "plugging_calibration.json")
    with open(p, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print("written %s" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
