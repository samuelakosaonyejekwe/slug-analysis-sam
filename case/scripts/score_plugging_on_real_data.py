#!/usr/bin/env python3
"""Score THIS solver's PLUGGING prediction against real, measured plug/no-plug outcomes,
and say exactly how far that evidence goes.

The claim this answers is the blunt one: the case-study line is a DESIGN case. No such
line was ever built, so its 565 figures cannot be a hindcast of anything. What can be
done -- and what this script does -- is run the same solver, unchanged, on lines that do
exist and whose outcome is known, and report the score without dressing it up.

Three sources of real outcomes are available and they are NOT equivalent:

  Roberts (MEng, Memorial)   TWO real lines on one offshore development with known
                             outcomes -- a rigid spool that hydrate-blocked, and a 2.2 km
                             flexible flowline that was shut in five hours past its
                             no-touch time and did NOT block. The solver is run on both.
                             A third, INHIBITED arm is a model counterfactual on the
                             spool, not an observation, and is labelled as such.
  Pham et al. (2020)         12 real flowloop tests, 3 of which plugged. The solver is
                             NOT run on these, and the reason is stated below rather
                             than buried: the variable the authors show controls the
                             outcome is anti-agglomerant dose, and this solver has no
                             anti-agglomerant. Running it there would score a model on
                             cases whose deciding mechanism it does not contain.
  das Neves et al. (2025)    12 real flowloop HOLDUP points. Not a plugging outcome, but
                             the only quantitative real-data score the hydraulics have,
                             so it is carried in the ledger.

Nothing here is fitted. Where real data BOUNDS a parameter, the bound is reported as a
bound. Where it does not, that is reported too.
"""
import json
import math
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(_ROOT, "validation", "data")


def _load(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def _wilson(k, n, z=1.96):
    """Wilson score interval. Used instead of k/n because at n = 1 or 2 the point
    estimate is 0 or 1 and reads as certainty; the interval is what is honest."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1.0 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def main():
    ledger = []
    rule = "=" * 78

    asf = _load("field_roberts_arm_asfound.json")
    inh = _load("field_roberts_arm_inhibited.json")
    cs2 = _load("field_roberts_arm_case2.json")
    pham = _load("flowloop_plugging_pham2020.json")
    dn = _load("flowloop_holdup_dasneves2025.json")

    print(rule)
    print(" THIS SOLVER vs REAL LINES -- PLUGGING SCORE AND VALIDATION LEDGER")
    print(rule)

    # ---------------------------------------------------------------- 1. the real field
    print("\n1. REAL FIELD SPOOL WITH A KNOWN OUTCOME  (Roberts, MEng, Memorial University)")
    print("   35 m x 124.1 mm rigid production spool, shut in, methanol treatment omitted.")
    print("   OBSERVED: hydrate blockage, confirmed by pressure building ~5300 -> 8000 kPa.")
    if asf is None:
        print("   [not yet run] -- run: python case/scripts/run_field_arm.py --arm asfound")
        hit = None
    else:
        pr = asf["predicted"]
        pp = pr.get("P_plug")
        print("     predicted P_plug            %8.3f" % pp)
        print("     predicted peak Phi_SH       %8.4g   (Phi_crit = 1.08)"
              % pr.get("max_Phi_SH", float("nan")))
        print("     predicted max subcooling    %8.2f C" % pr.get("max_subcooling_C", float("nan")))
        print("     predicted peak deposit      %8.2f mm" % pr.get("peak_deposit_mm", float("nan")))
        ttp = pr.get("time_to_plug_P50_h")
        print("     predicted time to plug P50  %8s h" % ("%.1f" % ttp if ttp else "n/a"))
        hit = bool(pp is not None and pp > 0.5)
        print("   SCORE: %s" % ("HIT  -- predicts blockage on a spool that blocked"
                                if hit else
                                "MISS -- does not predict blockage on a spool that blocked"))
        ledger.append({"source": "Roberts field spool", "kind": "real field outcome",
                       "n": 1, "quantity": "plug / no plug", "result": "hit" if hit else "miss",
                       "P_plug": pp})

    # ---------------------------------------------------------------- 2. the control arm
    print("\n2. CONTROL ARM -- the same spool with the treatment that was skipped")
    print("   This is NOT a field observation. No one measured this spool with methanol in")
    print("   it. It exists to catch the failure a single positive cannot: a model that")
    print("   predicts plugging everywhere passes case 1 and is worthless.")
    if inh is None:
        print("   [not yet run] -- run: python case/scripts/run_field_arm.py --arm inhibited")
        discriminates = None
    else:
        pr = inh["predicted"]
        pp_i = pr.get("P_plug")
        print("     inhibitor dosed             %8.1f wt%% on the aqueous phase" % inh["MEG_wt_inlet"])
        print("     predicted P_plug            %8.3f" % pp_i)
        print("     predicted peak Phi_SH       %8.4g" % pr.get("max_Phi_SH", float("nan")))
        print("     predicted max subcooling    %8.2f C" % pr.get("max_subcooling_C", float("nan")))
        if asf is not None:
            discriminates = bool(asf["predicted"]["P_plug"] > 0.5 >= pp_i)
            print("   DISCRIMINATION: %s" % (
                "PASSES -- plugs untreated, does not plug treated" if discriminates else
                "FAILS -- the model does not separate the treated case from the untreated one"))
            ledger.append({"source": "Roberts spool, inhibited control",
                           "kind": "model counterfactual, NOT an observation",
                           "n": 0, "quantity": "plug / no plug",
                           "result": "discriminates" if discriminates else "does not discriminate",
                           "P_plug": pp_i})

    # ------------------------------------------------------ 3. the real confirmed negative
    print("\n3. REAL FLOWLINE THAT DID NOT BLOCK  (Roberts case 2, same development)")
    print("   2200 m x 254 mm 10 in flexible, DC1 to the FPSO. Plant trip, then stagnant for")
    print("   11 h -- five hours BEYOND the field's 6 h no-touch time.")
    print("   OBSERVED: no blockage. Production restarted without incident.")
    print("   Reference model (OLGA + CSMHyK): <2 % hydrate at the drill-centre manifold and")
    print("   <1 % at the riser base at 11 h, accelerating only from about 17 h.")
    if cs2 is None:
        print("   [not yet run] -- run: python case/scripts/run_field_arm.py --arm case2")
        neg_ok = None
    else:
        pr = cs2["predicted"]
        pp2 = pr.get("P_plug")
        print("     predicted P_plug            %8.3f" % pp2)
        print("     predicted peak Phi_SH       %8.4g" % pr.get("max_Phi_SH", float("nan")))
        print("     predicted max subcooling    %8.2f C" % pr.get("max_subcooling_C", float("nan")))
        print("     predicted peak deposit      %8.2f mm" % pr.get("peak_deposit_mm", float("nan")))
        neg_ok = bool(pp2 is not None and pp2 <= 0.5)
        print("   SCORE: %s" % ("HIT  -- does not predict blockage on a line that did not block"
                                if neg_ok else
                                "FALSE ALARM -- predicts blockage on a line that did not block"))
        ledger.append({"source": "Roberts flowline 1a", "kind": "real field outcome",
                       "n": 1, "quantity": "plug / no plug", "observed": "no plug",
                       "result": "hit" if neg_ok else "false alarm", "P_plug": pp2})

    # ------------------------------------------------- 4. what the real outcome calibrates
    print("\n4. WHAT THE REAL OUTCOMES ACTUALLY CALIBRATE")
    if asf is not None:
        D = 0.1241
        dep = asf["predicted"].get("peak_deposit_mm")
        restr = 2.0 * (dep / 1000.0) / D if dep else float("nan")
        print("   The one plugging parameter a plug/no-plug observation can constrain is the")
        print("   trip: numerics.plug_restriction_trip, shipped at 0.85 of the bore and until")
        print("   now an ASSUMED number with nothing behind it.")
        print("     median realisation's peak restriction on the spool that blocked: %.3f" % restr)
        if cs2 is not None:
            dep2 = cs2["predicted"].get("peak_deposit_mm")
            restr2 = 2.0 * (dep2 / 1000.0) / 0.254 if dep2 else float("nan")
            print("     the same on the flowline that did NOT block:                    %.3f" % restr2)
            if restr2 == restr2 and restr == restr:
                if restr2 < restr:
                    print("     -> the two real outcomes BRACKET the trip: it must sit above %.3f" % restr2)
                    print("        (or the line that did not block would be called plugged) and at or")
                    print("        below %.3f (or the line that did block would be called clear)." % restr)
                    print("        The shipped 0.85 %s inside that bracket."
                          % ("SITS" if restr2 < 0.85 <= restr else "DOES NOT SIT"))
                else:
                    print("     -> the two real outcomes DO NOT bracket the trip: the line that did")
                    print("        not block reaches a HIGHER restriction than the one that did, so")
                    print("        no threshold on this variable separates them. That is a finding")
                    print("        about the model, and it is reported rather than tuned away.")
            ledger.append({"source": "Roberts, both real lines", "kind": "parameter bracket",
                           "parameter": "numerics.plug_restriction_trip", "shipped": 0.85,
                           "must_exceed": restr2, "must_not_exceed": restr,
                           "shipped_value_inside_bracket": bool(restr2 < 0.85 <= restr)})
        if restr == restr:
            if restr >= 0.85:
                print("     -> the shipped trip of 0.85 is REACHED on a line that really blocked,")
                print("        so the real outcome is CONSISTENT with it. That is a one-sided")
                print("        check: it bounds the trip ABOVE at %.2f and says nothing about how" % restr)
                print("        much lower it could be. A no-plug field case would bound it below;")
                print("        the one available (Roberts case 2) has no geometry tabulated in the")
                print("        open document, so it cannot be run and is not claimed.")
            else:
                print("     -> the shipped trip of 0.85 is NOT reached on a line that really")
                print("        blocked. Real data therefore bounds the trip ABOVE at %.3f, and" % restr)
                print("        0.85 is too high. This is a measured correction, not a preference.")
        ledger.append({"source": "Roberts field spool", "kind": "parameter bound",
                       "parameter": "numerics.plug_restriction_trip",
                       "shipped": 0.85, "upper_bound_from_real_data": restr,
                       "sided": "upper bound only -- one positive observation, no runnable negative"})
    else:
        print("   [pending the asfound arm]")

    # ---------------------------------------------------------------- 4. the flowloop set
    print("\n5. REAL FLOWLOOP PLUG OUTCOMES  (Pham et al. 2020, n = 12, 3 plugged)")
    if pham:
        npl = sum(1 for r in pham["points"] if r[9])
        print("   %d of %d tests plugged. The solver is NOT run on these, deliberately."
              % (npl, len(pham["points"])))
        print("   Every test but one carries anti-agglomerant and/or NaCl, and the authors show")
        print("   the outcome is set by an agglomerate structure factor Kv driven by that dose.")
        print("   This solver models WALL DEPOSITION against shear removal and has no")
        print("   anti-agglomerant term at all. Scoring it here would test a model on cases")
        print("   whose deciding variable it does not represent, and a near-chance result")
        print("   would say nothing about the model.")
        print("   What this set DOES establish, and it is already scored in")
        print("   score_pham2020_result.txt: inventory does not predict plugging (AUC 0.148)")
        print("   and a dimensionless structural factor does (AUC 0.778). That supports the")
        print("   SHAPE of the Phi_SH argument on real outcomes. It is not a calibration.")
        ledger.append({"source": "Pham et al. 2020 flowloop", "kind": "real plug outcomes",
                       "n": len(pham["points"]), "quantity": "structural claim only",
                       "result": "inventory AUC 0.148, structure AUC 0.778",
                       "solver_run": False,
                       "why_not": "no anti-agglomerant term in this solver; AA dose is the "
                                  "deciding variable in 11 of 12 tests"})

    # ------------------------------------------- 4b. what the wall simplification decides
    print("\n4b. THE WALL BOUNDS -- what the simplification does and does not decide")
    rw = _load("field_roberts_arm_case2_riserwall.json")
    #  The warm bound must be a FLOWLINE-WALL run. Falling back to the plain case2 arm was
    #  safe only while that arm was itself single-stack; now that case2 runs the SEGMENTED
    #  wall, the fallback would compare segmented against riser and label it "warm bound",
    #  understating the spread and mislabelling the model as one of its own bounds. Accept
    #  the fallback only when the record says it was a flowline-wall run.
    fw = _load("field_roberts_arm_case2_flowlinewall.json")
    if fw is None and cs2 is not None and cs2.get("wall_model") in (None, "flowline"):
        fw = cs2
    if rw and fw:
        print("   Segment 1a is insulated flexible flowline for ~2020 m and BARE riser beyond.")
        print("   The solver now runs it segmented; these two single-stack arms bound it.")
        print("     %-26s %12s %12s %8s" % ("quantity", "warm bound", "cold bound", "ratio"))
        for k, lab in (("P_plug", "P_plug"), ("max_subcooling_C", "max subcooling C"),
                       ("max_Phi_SH", "max Phi_SH"), ("peak_deposit_mm", "peak deposit mm")):
            a, b = fw["predicted"].get(k), rw["predicted"].get(k)
            if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
                continue
            r = (b / a) if a else float("nan")
            print("     %-26s %12.4g %12.4g %8s"
                  % (lab, a, b, "--" if r != r else "%.2f" % r))
        _pa, _pb = fw["predicted"].get("peak_deposit_mm"), rw["predicted"].get("peak_deposit_mm")
        _ra = 2.0 * (_pa / 1000.0) / 0.254 if _pa else float("nan")
        _rb = 2.0 * (_pb / 1000.0) / 0.254 if _pb else float("nan")
        print("     restriction 2d/D           %12.3f %12.3f   (trip 0.85)" % (_ra, _rb))
        _bit = (fw["predicted"].get("P_plug", 1) <= 0.5) == (rw["predicted"].get("P_plug", 1) <= 0.5)
        _mag = abs((_pb / _pa) - 1.0) > 0.25 if (_pa and _pb) else False
        print("   -> the BIT is %s across the bounds." % ("UNCHANGED" if _bit else "CHANGED"))
        if _mag:
            print("   -> the MAGNITUDE is NOT: peak deposit moves x%.2f. The wall does not decide"
                  % (_pb / _pa))
            print("      whether this line plugs, and it DOES decide how much deposit is predicted,")
            print("      which is one of the three things scored against the reference -- and the")
            print("      reference bound that is missed is at the riser base, exactly where the two")
            print("      walls differ. The bit survives only because BOTH bounds sit far below the")
            print("      trip; nearer it, the wall would decide the answer outright.")
        ledger.append({"source": "Roberts case 2 wall bounds", "kind": "sensitivity, not observation",
                       "bit_unchanged": bool(_bit),
                       "peak_deposit_ratio_cold_over_warm": (_pb / _pa) if (_pa and _pb) else None})
    else:
        print("   [pending] run both: --arm case2_riserwall and --arm case2_flowlinewall")

    # ---------------------------------------------------------------- 5. holdup
    print("\n6. REAL FLOWLOOP HOLDUP  (das Neves et al. 2025, n = 12 quick-closing-valve points)")
    if dn:
        ra = dn.get("residual_analysis", {})
        print("   Not a plugging outcome, but the only quantitative real-data score the")
        print("   hydraulics carry, and holdup is half of what drives Phi_SH.")
        print("   void RMSE 0.056 as shipped, 4 of 12 inside the measured uncertainty band.")
        if ra:
            print("   The residual is velocity-structured (corr +0.82 with mixture Froude), so a")
            print("   one-parameter C0 calibration buys 2.1 % and is not adopted. See")
            print("   residual_analysis in the dataset for what was tried and rejected.")
        ledger.append({"source": "das Neves et al. 2025 flowloop", "kind": "real measurements",
                       "n": 12, "quantity": "void fraction / holdup",
                       "result": "RMSE 0.056, 4/12 in band, residual corr +0.82 with Froude"})

    # ---------------------------------------------------------------- verdict
    print("\n" + rule)
    print(" WHAT THIS DOES AND DOES NOT ESTABLISH")
    print(rule)
    outcomes = []
    if asf is not None:
        outcomes.append(("spool that blocked", True, bool(hit)))
    if cs2 is not None:
        outcomes.append(("flowline that did not block", False, bool(neg_ok)))
    n_real = len(outcomes)
    k_hit = sum(1 for _, _, ok in outcomes if ok)
    lo, hi = _wilson(k_hit, n_real) if n_real else (0.0, 1.0)
    print(" Real lines this solver has been run on with a known plugging outcome: %d." % n_real)
    for nm, obs, ok in outcomes:
        print("   %-32s observed %-8s -> %s"
              % (nm, "PLUG" if obs else "no plug", "correct" if ok else "WRONG"))
    if n_real:
        print(" Correct: %d of %d. Wilson 95%% interval on the hit rate at that n: %.2f to %.2f."
              % (k_hit, n_real, lo, hi))
        print(" That interval is most of the unit interval, and that is the point: two")
        print(" observations cannot establish a hit rate. What they CAN do is something a")
        print(" single positive cannot -- one of them is a NEGATIVE, so the model has to")
        print(" get a line right by NOT calling it, and 'predict plugging everywhere' is")
        print(" no longer a way to pass.")
        if n_real >= 2 and k_hit == n_real:
            print(" Both real outcomes are called correctly, in both directions.")
        elif k_hit < n_real:
            print(" NOT every real outcome is called correctly. See the section above for")
            print(" which one and by how much; it is reported, not tuned away.")
    print("")
    print(" The case-study line is a DESIGN case. It was never built, so its 565 figures")
    print(" are a prediction for a designed line, not a hindcast of a measured one, and")
    print(" no amount of scoring elsewhere makes them one.")
    print("")
    print(" CALIBRATED ON REAL VALUES:  the hydrate equilibrium curve (published data,")
    print("   1-parameter offset); the holdup closure's residual (measured, reported, NOT")
    print("   fitted away)%s." % ("; the plug trip, bounded above by a real blockage"
                                  if asf is not None else
                                  " -- and NOT the plug trip, whose bounding run has not been made"))
    print(" ANCHORED BUT NOT CALIBRATED: k_g0, growth_exp_n, C_phi, k_ero -- literature")
    print("   ranges, not this line's measurements.")
    print(" NOT CALIBRATED AT ALL: nothing in the coupling number Phi_SH has been fitted to")
    print("   a measured plug. Phi_crit = 1.08 is a decomposition of assumed constants.")
    print(rule)

    out = os.path.join(DATA, "plugging_score_real.json")
    with open(out, "w") as fh:
        json.dump({"ledger": ledger,
                   "real_lines_with_known_plug_outcome": n_real,
                   "hits": k_hit,
                   "hit_rate_wilson95": [lo, hi],
                   "case_study_line": "DESIGN CASE -- never built, cannot be hindcast"},
                  fh, indent=2)
        fh.write("\n")
    print("\nwritten %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
