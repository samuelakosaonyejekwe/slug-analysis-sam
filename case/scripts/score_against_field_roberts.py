#!/usr/bin/env python3
"""Run THIS solver on a REAL field case with a known outcome, and score it.

Roberts (MEng, Memorial University) documents a rigid production spool on an offshore
development that PLUGGED with hydrate after a well shut-in when the planned methanol
treatment was not performed. Everything needed to set the case up is tabulated: spool
bore and length, the full wall build-up with thermal properties, seawater ambient, the
rates at shut-in, and a laboratory-validated crude assay.

This is the first time anything in this repository is run against a real line with a
known answer. The score is a single bit -- does the solver predict hydrate blockage in a
spool that blocked -- plus where and how fast it says it happens.

The fluid here is the FIELD assay, not this project's own re-cut composition, so the
solver is being asked about someone else's oil on someone else's geometry.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import solver  # noqa: E402

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                   "validation", "data", "field_olga_csmhyk_roberts.json")


def build_case(assay):
    c = solver.Case()
    #  --- geometry: 6 in rigid production spool, ~35 m, tree to manifold
    c.pipeline.length_m = 35.0
    c.pipeline.diameter_m = 0.1241
    c.pipeline.n_cells = 30
    c.pipeline.elevation_m = None
    #  --- wall: duplex 21.95 mm then SPU 55 mm. volumetric heat capacity = rho * cp
    c.pipeline.wall_layers = [[0.02195, 14.4, 7600.0 * 480.0],
                              [0.0550, 0.135, 700.0 * 1700.0]]
    #  --- ambient: seawater -1.8 C, 0.94 m/s current (a swept outer film)
    c.operating.T_seabed_C = -1.8
    c.pipeline.h_outer = 500.0
    #  --- rates at shut-in: 600 Sm3/d liquid, 7000 Sm3/d gas
    c.operating.q_liquid_insitu = 600.0 / 86400.0
    c.operating.q_gas_insitu_inlet = 7000.0 / 86400.0
    c.operating.P_inlet_bar = 53.0          # spool open to flowline pressure ~5300 kPa
    c.operating.T_inlet_C = 60.0
    c.operating.MEG_wt_inlet = 0.0          # the methanol treatment was NOT performed
    c.fluids.composition = assay
    #  --- the event: the well is shut in and the spool cools to seabed
    c.scenario.kind = "shutin"
    c.scenario.event_time_h = 0.5
    #  12 h, not 48. A 35 m spool is three orders of magnitude shorter than the case
    #  study's 32 km line, so its CFL timestep is correspondingly tiny and 48 h of
    #  simulated time is an enormous number of steps for a body that reaches seabed
    #  temperature within a couple of hours. The physics of interest -- cool down,
    #  subcool, deposit -- is complete well inside 12 h. The first attempt ran 1 h 40 m
    #  of wall clock without finishing, which is a modelling error, not a slow machine.
    c.numerics.t_end_h = 12.0
    c.numerics.n_ensemble = 8
    c.numerics.seed = 7
    return c


def main():
    with open(REF) as fh:
        ref = json.load(fh)
    raw = ref["field_crude_assay"]["composition_mol_pct"]
    #  lump the assay onto this solver's component set; everything C7 and heavier -> C7+
    keys = {"N2": ["N2"], "CO2": ["CO2"], "C1": ["C1"], "C2": ["C2"], "C3": ["C3"],
            "iC4": ["iC4"], "nC4": ["nC4"], "iC5": ["iC5", "neoC5"], "nC5": ["nC5"],
            "C6": ["C6", "methylcyclopentane", "benzene", "cyclohexane"]}
    comp = {k: sum(raw.get(s, 0.0) for s in src) for k, src in keys.items()}
    heavy = sum(v for k, v in raw.items()
                if k not in {s for src in keys.values() for s in src} and k != "H2S")
    comp["C7+"] = heavy
    tot = sum(comp.values())
    comp = {k: v / tot for k, v in comp.items()}

    print("=" * 74)
    print(" SCORING THIS SOLVER ON A REAL FIELD CASE THAT PLUGGED")
    print(" Roberts, MEng, Memorial University -- rigid production spool, shut-in")
    print("=" * 74)
    print("  spool      35 m x 124.1 mm ID, duplex 21.95 mm + SPU 55 mm")
    print("  ambient    seawater -1.8 C")
    print("  at shut-in 600 Sm3/d liquid, 7000 Sm3/d gas, NO methanol")
    print("  assay      field crude, C1 %.2f mol%% -> lumped C7+ %.3f" % (raw["C1"], comp["C7+"]))
    print("  KNOWN OUTCOME: hydrate blockage. Pressure built ~5300 -> 8000 kPa.")
    print()
    sv = solver.TransientSHCT(build_case(comp))
    sv.run(verbose=False)
    e = sv.engineering()
    P_plug = float(e.get("P_plug", float("nan")))
    print("  SOLVER SAYS")
    print("    max subcooling        %8.2f C" % float(e.get("max_subcooling_C", float("nan"))))
    print("    peak Phi_SH           %8.3g   (Phi_crit = 1.08)" % float(e.get("max_Phi_SH", float("nan"))))
    print("    peak wall deposit     %8.2f mm" % float(e.get("peak_deposit_mm", float("nan"))))
    print("    plug probability      %8.2f" % P_plug)
    ttp = e.get("time_to_plug_P50_h")
    print("    time to plug (P50)    %8s h" % ("%.1f" % ttp if ttp and ttp == ttp else "n/a"))
    print()
    hit = P_plug > 0.5
    print("  SCORE: %s" % ("HIT -- predicts blockage on a spool that blocked"
                           if hit else
                           "MISS -- does not predict blockage on a spool that blocked"))
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
