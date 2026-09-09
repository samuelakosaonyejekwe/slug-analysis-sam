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

from _paths import CASE  # noqa: E402  (also installs the no-black style)

#  (manuscript figure number, preferred source, fallback source)
#  1-8 are the PVT, diagnostics and verification set; 9-17 the as-operated case;
#  18-20 the cross-section reconstruction; 21-23 the quasi-3-D and compositional
#  set; 24-27 the shut-in, mitigated and cross-scenario comparison; 28 the
#  sensitivity sweep; 29-42 the space-time / multi-time set.
#
#  DUPLICATE, LEFT AS IT IS AND SAID OUT LOUD: 16 and 19 both export
#  cx2_azimuthal_deposit.png, so Figure_16.png and Figure_19.png are byte-identical.
#  Which of the two the manuscript means is a question about the manuscript, not about
#  this map, so it is not guessed here -- export() reports every duplicate it writes.
#
#  THE MANUSCRIPT ENDS AT 27.  Entries beyond N_MANUSCRIPT are generated and
#  captioned but are cited by no document; they are exported only on request,
#  because the journal figure folder must contain exactly the figures the paper
#  cites and nothing else -- 42 files against 37 citations is a submission defect
#  an editor will bounce. Use --all (or export(extras=True)) to get them.
FIGURE_MAP = [
    ( 1, "outputs_paper_steady/compo_pvt.png",            "outputs_steady/compo_pvt.png"),
    ( 2, "outputs_paper_steady/08_diagnostics.png",       "outputs_steady/08_diagnostics.png"),
    ( 3, "outputs_paper_steady/hydrate_validation.png",   "outputs_steady/hydrate_validation.png"),
    ( 4, "outputs_paper_steady/friction_validation.png",  "outputs_steady/friction_validation.png"),
    ( 5, "outputs_paper_steady/verif_thermal_exact.png",  "outputs_steady/verif_thermal_exact.png"),
    ( 6, "outputs_paper_steady/verif_grid_convergence.png", "outputs_steady/verif_grid_convergence.png"),
    ( 7, "outputs_paper_steady/verif_cross_engine.png",   "outputs_steady/verif_cross_engine.png"),
    ( 8, "outputs_paper_steady/verif_water_faucet.png",   "outputs_steady/verif_water_faucet.png"),
    ( 9, "outputs_paper_steady/01_profiles.png",          "outputs_steady/01_profiles.png"),
    (10, "outputs_paper_steady/02_holdup_spacetime.png",  "outputs_steady/02_holdup_spacetime.png"),
    (11, "outputs_paper_steady/03_PT_envelope.png",       "outputs_steady/03_PT_envelope.png"),
    (12, "outputs_paper_steady/09_slug_prediction.png",   "outputs_steady/09_slug_prediction.png"),
    (13, "outputs_paper_steady/10_riser_severe_slug.png", "outputs_steady/10_riser_severe_slug.png"),
    (14, "outputs_paper_steady/04_PhiSH_map.png",         "outputs_steady/04_PhiSH_map.png"),
    (15, "outputs_paper_steady/06_deposit.png",           "outputs_steady/06_deposit.png"),
    (16, "outputs_paper_steady/cx2_azimuthal_deposit.png", "outputs_steady/cx2_azimuthal_deposit.png"),
    (17, "outputs_paper_steady/07_probabilistic.png",     "outputs_steady/07_probabilistic.png"),
    (18, "outputs_paper_steady/cx1_geometry.png",         "outputs_steady/cx1_geometry.png"),
    (19, "outputs_paper_steady/cx2_azimuthal_deposit.png", "outputs_steady/cx2_azimuthal_deposit.png"),
    (20, "outputs_paper_steady/cx3_sections.png",         "outputs_steady/cx3_sections.png"),
    (21, "outputs_paper_steady/threed_deposit.png",       "outputs_steady/threed_deposit.png"),
    (22, "outputs_paper_steady/threed_temperature.png",   "outputs_steady/threed_temperature.png"),
    (23, "outputs_paper_steady/compositional_transport.png", "outputs_steady/compositional_transport.png"),
    (24, "outputs_paper_shutin/01_profiles.png",          "outputs_shutin/01_profiles.png"),
    (25, "outputs_paper_shutin/04_PhiSH_map.png",         "outputs_shutin/04_PhiSH_map.png"),
    (26, "outputs_paper_mitigated/04_PhiSH_map.png",      "outputs_mitigated/04_PhiSH_map.png"),
    (27, "outputs_paper_steady/12_mitigation_comparison.png", "outputs_steady/12_mitigation_comparison.png"),
    (28, "outputs_steady/13_sensitivity.png",             "outputs_steady/13_sensitivity.png"),
    (29, "outputs_paper_steady/19_spacetime_fields.png",  "outputs_steady/19_spacetime_fields.png"),
    (30, "outputs_paper_steady/14_holdup_multitime.png",  "outputs_steady/14_holdup_multitime.png"),
    (31, "outputs_paper_steady/15_slug_growth_propagation.png",
         "outputs_steady/15_slug_growth_propagation.png"),
    (32, "outputs_paper_steady/16_slug_train_waterfall.png", "outputs_steady/16_slug_train_waterfall.png"),
    (33, "outputs_paper_steady/21_riser_depth_time.png",  "outputs_steady/21_riser_depth_time.png"),
    (34, "outputs_paper_steady/17_hydrate_distribution.png", "outputs_steady/17_hydrate_distribution.png"),
    (35, "outputs_paper_steady/22_cloud_maps.png",        "outputs_steady/22_cloud_maps.png"),
    (36, "outputs_paper_shutin/20_holdup_durations.png",  "outputs_shutin/20_holdup_durations.png"),
    (37, "outputs_paper_shutin/18_shutin_profile_deposit.png",
         "outputs_shutin/18_shutin_profile_deposit.png"),
    (38, "outputs_paper_steady/23_dts_thermal_waterfall.png", "outputs_steady/23_dts_thermal_waterfall.png"),
    (39, "outputs_paper_steady/24_temperature_gradient.png", "outputs_steady/24_temperature_gradient.png"),
    (40, "outputs_paper_steady/25_das_flow_noise.png",    "outputs_steady/25_das_flow_noise.png"),
    (41, "outputs_paper_shutin/26_parameter_panels.png",  "outputs_shutin/26_parameter_panels.png"),
    (42, "outputs_paper_steady/27_wellposedness_map.png", "outputs_steady/27_wellposedness_map.png"),
]

#  The last figure the manuscript actually cites. Everything above this number is
#  an extra: real output, real caption, no citation anywhere in paper5.docx.
N_MANUSCRIPT = 37

def check_provenance(verbose=True):
    """Warn for every manuscript figure taken from the FALLBACK source.

    The map carries a journal source (outputs_paper_*, 320 dpi, chart titles
    suppressed) and a fallback (outputs_*, the report set: lower dpi, titles ON).
    The fallback is silent, so a figure missing from the journal set ships with a
    chart title no other figure has, duplicating its own caption -- which is how
    Figs 5-8 (the verification set) nearly went out. Sizes and counts cannot see
    it, because the fallback files are large enough to pass the pixel rule.
    """
    import os
    bad = []
    for n, primary, fallback in FIGURE_MAP:
        if n > N_MANUSCRIPT:
            continue
        if not os.path.exists(os.path.join(CASE, primary)):
            bad.append((n, primary, fallback))
    if verbose:
        if bad:
            print(f"[provenance] {len(bad)} manuscript figure(s) fell back to the report set:")
            for n, pr, fb in bad:
                print(f"    Figure_{n}: missing {pr} -> using {fb}")
        else:
            print(f"[provenance] all {N_MANUSCRIPT} manuscript figures come from the journal set")
    return bad


#  the caption of every figure, so the manuscript and the deck stay in step with
#  what the run actually produced
#  KEYED BY THE MANUSCRIPT FIGURE NUMBER IN FIGURE_MAP ABOVE, and they were not: every
#  one of these fourteen captions was keyed ten low, so "the space-time solution of the
#  tie-back" sat on 19 (the azimuthal-deposit section) instead of 29, and so on down the
#  list. Nothing in this repository reads the table, which is why it drifted unnoticed --
#  but reembed_figures.py names it "the authority" for which output is which figure.
CAPTIONS = {
    29: ("The space-time solution of the tie-back, as-operated: liquid holdup, pressure, "
         "gas and liquid velocity, wall subcooling and the wall-deposit volume fraction, "
         "each as a filled-contour field over distance and time. The deposit panel shows the "
         "wall deposit establishing from roughly 15 km onward once the subcooling does; "
         "it stabilises at 4.1 mm and the bore does not close."),
    30: ("Liquid holdup along the whole route at successive times: the early transient "
         "(upper) and the late, quasi-developed state (lower). Terrain-locked accumulation "
         "in the 10-20 km band and drainage toward the riser base are both visible."),
    31: ("Slug propagation and front tracking over a short reach at three successive times. "
         "One front is followed across the panels; its arrival time and position give the "
         "translational celerity directly. Sub-grid reconstruction (see text)."),
    32: ("Slug tracking in the space-time plane: (a) the distance-time waterfall of a single "
         "slug unit, (b) semblance against trial celerity, (c) the waterfall after linear "
         "moveout at the recovered celerity, and (d) the distance-stacked trace. The "
         "recovered celerity returns the solver's own translational velocity. Sub-grid "
         "reconstruction (see text)."),
    33: ("Depth-time waterfall over the steel-catenary riser: slug boundaries during upward "
         "motion, their trajectories, and the slug unit length projected onto the depth "
         "axis. Sub-grid reconstruction (see text)."),
    34: ("(a) In-pipe volume fractions along the route: unconverted water, hydrate carried "
         "in the liquids, and the hydrate deposit standing on the wall. (b) The gas, oil and "
         "water mass rates delivered into the host separator against time."),
    35: ("Pipeline cloud maps at successive times: the gas-liquid phase distribution inside "
         "the bore (upper strip of each pair) above the bulk-temperature field along the "
         "same reach (lower strip), on a shared temperature scale."),
    36: ("Distribution of liquid holdup along the pipeline after different shut-in "
         "durations. Liquid drains from the crests into the low spots as the line cools."),
    37: ("(a) The pipeline profile late in the shut-in: pressure, temperature against the "
         "hydrate-equilibrium temperature, and the water volume fraction. (b) The "
         "wall-deposit volume fraction along the line at successive elapsed times."),
    38: ("Distributed-temperature waterfall of the as-operated line: the thermal field "
         "over distance and time, with the monitored pressure overlaid and the "
         "hydrate-onset distance annotated. The cold section develops from about 10 km "
         "outward and holds for the rest of the run."),
    39: ("Temperature-gradient waterfall. A travelling thermal front appears as a narrow "
         "band of steep gradient, so it is localised here even where the temperature "
         "field itself is smooth; the dashed line tracks the steepest cooling at each "
         "instant."),
    40: ("Flow-noise waterfall: the rate of change of liquid holdup over distance and "
         "time. The unsteadiness is concentrated in the intermittent reach and at the "
         "riser, and settles as the line reaches its quasi-developed state."),
    41: ("Pressure, temperature, liquid holdup and mixture velocity along the route at "
         "successive times after the shut-in."),
    42: ("Well-posedness of the two-fluid description. (a) The inviscid Kelvin-Helmholtz "
         "boundary over the superficial-velocity plane, with the states the case "
         "actually occupies. (b) The slip against that limit along the route. The margin "
         "peaks at 1.99 and exceeds unity over 4.3 % of the route -- a short reach at "
         "the riser base where the film is thin and fast. Over that reach the "
         "initial-value problem is ill-posed and the growth rate is grid-dependent, so "
         "slug activity localised there should not be read as a property of the flow; "
         "over the remaining 96 % it can be."),
}


def export(target=None, verbose=True, extras=False):
    """Write the manuscript figure set. `extras=True` also writes 28+."""
    target = target or os.path.join(CASE, "figures_paper")
    os.makedirs(target, exist_ok=True)
    written, missing, skipped = [], [], []
    for num, primary, fallback in FIGURE_MAP:
        if num > N_MANUSCRIPT and not extras:
            skipped.append(num)
            continue
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
    #  Two manuscript numbers can be mapped to the SAME source, which writes two
    #  byte-identical Figure_N.png files. That is a submission defect an editor will
    #  bounce, and nothing here reported it: 16 and 19 have both pointed at
    #  cx2_azimuthal_deposit.png. Report it; do not guess which one is wrong.
    dup: dict[str, list[int]] = {}
    for num, rel in written:
        dup.setdefault(rel, []).append(num)
    dups = {rel: ns for rel, ns in dup.items() if len(ns) > 1}
    if verbose:
        check_provenance()
        for rel, ns in sorted(dups.items()):
            print("[export] DUPLICATE: figures "
                  + ", ".join(str(n) for n in ns)
                  + f" are all {rel} — the manuscript must not carry the same image twice",
                  flush=True)
        print(f"[export] {len(written)} figures -> {target}"
              f"{f'  ({len(missing)} missing)' if missing else ''}", flush=True)
        if skipped:
            print(f"[export] {len(skipped)} beyond the manuscript not written "
                  f"({skipped[0]}-{skipped[-1]}); pass --all to include them", flush=True)
    return written, missing


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--all"]
    export(args[0] if args else None, extras="--all" in sys.argv)
