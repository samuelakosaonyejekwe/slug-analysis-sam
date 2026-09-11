#!/usr/bin/env python3
"""Run the ENTIRE case-study output pipeline on a line that actually exists.

The standing objection to this project is exact and fair: the 32 km tie-back the 565
figures describe is a DESIGN case. It was never built. So those figures are a prediction
for a line nobody can check, and no amount of closure validation changes that.

This script is the answer that is actually available. It takes a REAL line -- Roberts
case 2, a 2.2 km 10 in flexible flowline on a producing offshore development, with a
published geometry, a published wall build-up, published rates measured by multiphase
meters the morning of the event, published pressures, and a KNOWN OUTCOME -- and runs the
same solver, the same settings and the same figure pipeline that produced the case study.

Same code, same charts, real line, known answer. Whatever the figures say here, they said
it about something checkable.

What the line did: an emergency shutdown left it stagnant for 11 h, five hours beyond its
6 h no-touch time. It did NOT block. Production restarted without incident. The reference
model (OLGA + CSMHyK, a licensed transient code with the Colorado School of Mines hydrate
kinetics) put hydrate below 2 vol% at the drill-centre manifold and below 1 vol% at the
riser base at 11 h, accelerating only from about 17 h.

So there are three things to check, not one:
  * the BIT      -- does this solver say "no plug" on a line that did not plug?
  * the MAGNITUDE-- does its hydrate fraction at 11 h land near OLGA's <2 % and <1 %?
  * the TIMING   -- does it also show growth accelerating from ~17 h and not before?

A model can pass the first and fail the other two, which is why all three are reported.
"""
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import run_case_study10 as cs  # noqa: E402  (the case-study pipeline itself)
import run_field_arm as rfa  # noqa: E402  (the real-line case builder)

import shct_crosssection  # noqa: E402

OUT = os.path.join(_ROOT, "case", "outputs_realline")
REF = os.path.join(_ROOT, "validation", "data", "field_olga_csmhyk_roberts.json")


def main():
    os.makedirs(OUT, exist_ok=True)
    comp = rfa.field_assay()
    #  20 h, so the 11 h real point AND the ~17 h acceleration the reference reports are
    #  both inside the window. 20 realisations, not the 8 the quick arm used, because
    #  P_plug is quantised at 1/n and a figure set should not carry a 12.5 % quantum.
    case = rfa.build_case2(comp, n_ensemble=20, t_end_h=20.0)

    print("=" * 78)
    print(" THE FULL CASE-STUDY PIPELINE, RUN ON A REAL LINE")
    print("=" * 78)
    print("  line     Roberts case 2 -- flowline 1a, DC1 to FPSO, 2200 m x 254 mm")
    print("  event    plant trip, stagnant 11 h (5 h beyond the 6 h no-touch time)")
    print("  OUTCOME  DID NOT BLOCK. Restarted without incident.")
    print("  ref      OLGA + CSMHyK: <2 vol%% hydrate at the manifold, <1 %% at the riser")
    print("           base at 11 h; acceleration only from ~17 h")
    print("=" * 78)

    sv, eng = cs.run_core(case, OUT, slug=True, riser=True)
    shct_crosssection.crosssection_outputs(sv, OUT)

    # ---------------------------------------------------------------- the three checks
    r = sv.results
    st = np.asarray(r.get("snap_t", np.empty(0)), float)
    sp = np.asarray(r.get("snap_phi", np.empty(0)), float)      # (nsnap, nx) mean over ens
    x_km = sv.x / 1000.0
    out = {"line": "Roberts case 2 -- flowline 1a, DC1 to FPSO (REAL, did not block)",
           "observed_plugged": False}

    #  --- 1. the bit
    P_plug = float(eng.get("P_plug", float("nan")))
    out["bit"] = {"P_plug": P_plug, "observed": "no plug",
                  "correct": bool(P_plug <= 0.5)}

    #  --- 2. the magnitude, at the two places the reference names and the time it names
    if st.size and sp.size:
        i11 = int(np.argmin(np.abs(st - 11.0)))
        #  The drill-centre manifold is the inlet end. The RISER BASE is where the line
        #  leaves the seabed for the FPSO -- and the obvious test for that, "first cell
        #  where the elevation rises", is wrong on this geometry: flowline 1a dips below
        #  datum over its first 60 m and climbs back to the seabed plateau by ~90 m, so
        #  the first rise is that recovery and the riser base came out at x = 0.09 km
        #  instead of ~2.02 km. Scoring then compared this solver at 90 m against a
        #  reference value quoted 2 km away.
        #
        #  Anchor on the plateau instead: the seabed run is the modal elevation over the
        #  first 80 % of the line, and the riser base is the LAST cell still sitting on it.
        z = np.asarray(case.pipeline.elevation_m, float)
        if z.size:
            plateau = float(np.median(z[: max(1, int(0.8 * z.size))]))
            on_bed = np.flatnonzero(np.abs(z - plateau) <= 2.0)
            i_riser = int(on_bed[-1]) if on_bed.size else z.size - 1
        else:
            i_riser = 0
        out["magnitude_at_11h"] = {
            "snapshot_time_h": float(st[i11]),
            "manifold_x_km": float(x_km[0]),
            "manifold_phi_volfrac": float(sp[i11, 0]),
            "reference_manifold_volfrac": "< 0.02 (OLGA + CSMHyK)",
            "riser_base_x_km": float(x_km[i_riser]),
            "riser_base_phi_volfrac": float(sp[i11, i_riser]),
            "reference_riser_base_volfrac": "< 0.01 (OLGA + CSMHyK)",
            "line_max_phi_volfrac": float(np.nanmax(sp[i11])),
            "line_max_x_km": float(x_km[int(np.nanargmax(sp[i11]))]),
            "reference_locations": "OLGA + CSMHyK put the hydrate at the drill-centre "
                                   "manifold and the base of the riser",
        }
        out["magnitude_at_11h"]["manifold_within_reference"] = bool(sp[i11, 0] < 0.02)
        out["magnitude_at_11h"]["riser_base_within_reference"] = bool(sp[i11, i_riser] < 0.01)

        #  --- 3. the timing: when does line-max hydrate start climbing fast?
        mx = np.nanmax(sp, axis=1)
        if mx.size > 3 and np.nanmax(mx) > 1e-9:
            #  "acceleration" = first snapshot where the growth rate exceeds twice its
            #  median over the record. A definition, stated, not tuned to the answer.
            d = np.gradient(mx, st)
            med = float(np.nanmedian(d[d > 0])) if np.any(d > 0) else 0.0
            fast = np.flatnonzero(d > 2.0 * med) if med > 0 else np.array([], int)
            out["timing"] = {
                "phi_max_final": float(mx[-1]),
                "acceleration_onset_h": float(st[fast[0]]) if fast.size else None,
                "reference_acceleration_onset_h": 17.0,
                "definition": "first snapshot whose growth rate exceeds 2x the median "
                              "positive growth rate over the record",
            }

    with open(os.path.join(OUT, "real_line_score.json"), "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    print("\n" + "=" * 78)
    print(" SCORE ON THE REAL LINE")
    print("=" * 78)
    print("  1. THE BIT        P_plug = %.3f   observed: no plug   -> %s"
          % (P_plug, "correct" if out["bit"]["correct"] else "WRONG"))
    m = out.get("magnitude_at_11h")
    if m:
        print("  2. THE MAGNITUDE  at %.1f h" % m["snapshot_time_h"])
        print("       manifold   (x = %.2f km)  phi = %.4f vol   ref < 0.02   -> %s"
              % (m["manifold_x_km"], m["manifold_phi_volfrac"],
                 "agrees" if m["manifold_within_reference"] else "ABOVE the reference"))
        print("       riser base (x = %.2f km)  phi = %.4f vol   ref < 0.01   -> %s"
              % (m["riser_base_x_km"], m["riser_base_phi_volfrac"],
                 "agrees" if m["riser_base_within_reference"] else "ABOVE the reference"))
        print("       line maximum (x = %.2f km)  phi = %.4f vol"
              % (m.get("line_max_x_km", float("nan")), m["line_max_phi_volfrac"]))
        #  WHERE the peak sits is a result in its own right, and it is the one that can
        #  disagree with the reference while both named points agree.
        _near_ref = (m.get("line_max_x_km", 1e9) < 0.15
                     or abs(m.get("line_max_x_km", 1e9) - m["riser_base_x_km"]) < 0.15)
        print("       -> peak location %s the reference, which puts hydrate at the"
              % ("AGREES with" if _near_ref else "DISAGREES with"))
        print("          manifold and the riser base.")
    t = out.get("timing")
    if t:
        on = t["acceleration_onset_h"]
        print("  3. THE TIMING     growth accelerates at %s h; reference says ~17 h"
              % ("%.1f" % on if on is not None else "never in window"))
    print("=" * 78)
    print("  Figures, tables and space-time maps for this REAL line are in")
    print("  case/outputs_realline/ -- the same set, from the same code, as the")
    print("  design case study. That is the comparison this project did not have.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
