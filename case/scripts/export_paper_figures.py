#!/usr/bin/env python3
# =============================================================================
#  export_paper_figures.py — assemble the manuscript figure set from the case
#  study outputs.
# -----------------------------------------------------------------------------
#  The manuscript refers to its figures as Figure_1 ... Figure_N.  This script is
#  the single, version-controlled definition of WHICH generated output becomes
#  WHICH manuscript figure, so the figure folder can always be rebuilt from a
#  fresh run instead of being curated by hand.
#
#      python3 export_paper_figures.py                 # -> case/figures_paper/
#      python3 export_paper_figures.py <target-dir>    # -> anywhere else
#
#  The untitled ("paper") variants are preferred where they exist, because the
#  journal caption carries the description; the titled outputs are the fallback.
# =============================================================================
import os
import shutil
import sys

from _paths import CASE      # noqa: E402  (also installs the no-black style)

#  (manuscript figure number, preferred source, fallback source)
#  Figures 1-13 are the as-operated case, 14-16 the shut-in and mitigated
#  scenarios, 17-18 the cross-scenario comparison and the sensitivity sweep.
#  19-25 are the space-time / multi-time set added in v3.2.
FIGURE_MAP = [
    (1,  "outputs_paper_steady/compo_pvt.png",              "outputs_steady/compo_pvt.png"),
    (2,  "outputs_paper_steady/08_diagnostics.png",         "outputs_steady/08_diagnostics.png"),
    (3,  "outputs_paper_steady/hydrate_validation.png",     "outputs_steady/hydrate_validation.png"),
    (4,  "outputs_paper_steady/friction_validation.png",    "outputs_steady/friction_validation.png"),
    (5,  "outputs_paper_steady/01_profiles.png",            "outputs_steady/01_profiles.png"),
    (6,  "outputs_paper_steady/02_holdup_spacetime.png",    "outputs_steady/02_holdup_spacetime.png"),
    (7,  "outputs_paper_steady/03_PT_envelope.png",         "outputs_steady/03_PT_envelope.png"),
    (8,  "outputs_paper_steady/09_slug_prediction.png",     "outputs_steady/09_slug_prediction.png"),
    (9,  "outputs_paper_steady/10_riser_severe_slug.png",   "outputs_steady/10_riser_severe_slug.png"),
    (10, "outputs_paper_steady/04_PhiSH_map.png",           "outputs_steady/04_PhiSH_map.png"),
    (11, "outputs_paper_steady/06_deposit.png",             "outputs_steady/06_deposit.png"),
    (12, "outputs_paper_steady/cx2_azimuthal_deposit.png",  "outputs_steady/cx2_azimuthal_deposit.png"),
    (13, "outputs_paper_steady/07_probabilistic.png",       "outputs_steady/07_probabilistic.png"),
    (14, "outputs_paper_shutin/01_profiles.png",            "outputs_shutin/01_profiles.png"),
    (15, "outputs_paper_shutin/04_PhiSH_map.png",           "outputs_shutin/04_PhiSH_map.png"),
    (16, "outputs_paper_mitigated/04_PhiSH_map.png",        "outputs_mitigated/04_PhiSH_map.png"),
    (17, "outputs_paper_steady/12_mitigation_comparison.png", "outputs_steady/12_mitigation_comparison.png"),
    (18, "outputs_steady/13_sensitivity.png",               "outputs_steady/13_sensitivity.png"),
    #  ---- v3.2: the space-time / multi-time set --------------------------------
    (19, "outputs_paper_steady/19_spacetime_fields.png",        "outputs_steady/19_spacetime_fields.png"),
    (20, "outputs_paper_steady/14_holdup_multitime.png",        "outputs_steady/14_holdup_multitime.png"),
    (21, "outputs_paper_steady/15_slug_growth_propagation.png", "outputs_steady/15_slug_growth_propagation.png"),
    (22, "outputs_paper_steady/16_slug_train_waterfall.png",    "outputs_steady/16_slug_train_waterfall.png"),
    (23, "outputs_paper_steady/21_riser_depth_time.png",        "outputs_steady/21_riser_depth_time.png"),
    (24, "outputs_paper_steady/17_hydrate_distribution.png",    "outputs_steady/17_hydrate_distribution.png"),
    (25, "outputs_paper_steady/22_cloud_maps.png",              "outputs_steady/22_cloud_maps.png"),
    (26, "outputs_paper_shutin/20_holdup_durations.png",        "outputs_shutin/20_holdup_durations.png"),
    (27, "outputs_paper_shutin/18_shutin_profile_deposit.png",  "outputs_shutin/18_shutin_profile_deposit.png"),
    (28, "outputs_paper_steady/23_dts_thermal_waterfall.png",   "outputs_steady/23_dts_thermal_waterfall.png"),
    (29, "outputs_paper_steady/24_temperature_gradient.png",    "outputs_steady/24_temperature_gradient.png"),
    (30, "outputs_paper_steady/25_das_flow_noise.png",          "outputs_steady/25_das_flow_noise.png"),
    (31, "outputs_paper_shutin/26_parameter_panels.png",        "outputs_shutin/26_parameter_panels.png"),
    (32, "outputs_paper_steady/27_wellposedness_map.png",       "outputs_steady/27_wellposedness_map.png"),
]

#  the caption of every figure, so the manuscript and the deck stay in step with
#  what the run actually produced
CAPTIONS = {
    19: ("The space-time solution of the tie-back, as-operated: liquid holdup, pressure, "
         "gas and liquid velocity, wall subcooling and the wall-deposit volume fraction, "
         "each as a filled-contour field over distance and time. The deposit panel shows "
         "the bore closing from roughly 10 km onward once the subcooling establishes."),
    20: ("Liquid holdup along the whole route at successive times: the early transient "
         "(upper) and the late, quasi-developed state (lower). Terrain-locked accumulation "
         "in the 10-20 km band and drainage toward the riser base are both visible."),
    21: ("Slug propagation and front tracking over a short reach at three successive times. "
         "One front is followed across the panels; its arrival time and position give the "
         "translational celerity directly. Sub-grid reconstruction (see text)."),
    22: ("Slug tracking in the space-time plane: (a) the distance-time waterfall of a single "
         "slug unit, (b) semblance against trial celerity, (c) the waterfall after linear "
         "moveout at the recovered celerity, and (d) the distance-stacked trace. The "
         "recovered celerity returns the solver's own translational velocity. Sub-grid "
         "reconstruction (see text)."),
    23: ("Depth-time waterfall over the steel-catenary riser: slug boundaries during upward "
         "motion, their trajectories, and the slug unit length projected onto the depth "
         "axis. Sub-grid reconstruction (see text)."),
    24: ("(a) In-pipe volume fractions along the route: unconverted water, hydrate carried "
         "in the liquids, and the hydrate deposit standing on the wall. (b) The gas, oil and "
         "water mass rates delivered into the host separator against time."),
    25: ("Pipeline cloud maps at successive times: the gas-liquid phase distribution inside "
         "the bore (upper strip of each pair) above the bulk-temperature field along the "
         "same reach (lower strip), on a shared temperature scale."),
    26: ("Distribution of liquid holdup along the pipeline after different shut-in "
         "durations. Liquid drains from the crests into the low spots as the line cools."),
    27: ("(a) The pipeline profile late in the shut-in: pressure, temperature against the "
         "hydrate-equilibrium temperature, and the water volume fraction. (b) The "
         "wall-deposit volume fraction along the line at successive elapsed times."),
    28: ("Distributed-temperature waterfall of the as-operated line: the thermal field "
         "over distance and time, with the monitored pressure overlaid and the "
         "hydrate-onset distance annotated. The cold section develops from about 10 km "
         "outward and holds for the rest of the run."),
    29: ("Temperature-gradient waterfall. A travelling thermal front appears as a narrow "
         "band of steep gradient, so it is localised here even where the temperature "
         "field itself is smooth; the dashed line tracks the steepest cooling at each "
         "instant."),
    30: ("Flow-noise waterfall: the rate of change of liquid holdup over distance and "
         "time. The unsteadiness is concentrated in the intermittent reach and at the "
         "riser, and decays late in the run as the deposit closes the bore."),
    31: ("Pressure, temperature, liquid holdup and mixture velocity along the route at "
         "successive times after the shut-in."),
    32: ("Well-posedness of the two-fluid description. (a) The inviscid Kelvin-Helmholtz "
         "boundary over the superficial-velocity plane, with the states the case "
         "actually occupies. (b) The slip against that limit along the route. The margin "
         "stays below unity everywhere, so the predicted slug activity is a property of "
         "the flow and not a grid-dependent artefact of an ill-posed problem."),
}


def export(target=None, verbose=True):
    target = target or os.path.join(CASE, "figures_paper")
    os.makedirs(target, exist_ok=True)
    written, missing = [], []
    for num, primary, fallback in FIGURE_MAP:
        for rel in (primary, fallback):
            src = os.path.join(CASE, rel)
            if os.path.exists(src):
                dst = os.path.join(target, f"Figure_{num}.png")
                shutil.copyfile(src, dst)
                written.append((num, rel))
                if verbose:
                    print(f"  Figure_{num:<2d} <- {rel}", flush=True)
                break
        else:
            missing.append((num, primary))
            if verbose:
                print(f"  Figure_{num:<2d} MISSING ({primary})", flush=True)
    if verbose:
        print(f"[export] {len(written)} figures -> {target}"
              f"{f'  ({len(missing)} missing)' if missing else ''}", flush=True)
    return written, missing


if __name__ == "__main__":
    export(sys.argv[1] if len(sys.argv) > 1 else None)
