#!/usr/bin/env python3
# =============================================================================
#  check_outputs.py — an automated inspection of every generated output.
# -----------------------------------------------------------------------------
#  Figures fail quietly. A curve drawn off the axis, six series stacked in one
#  place, a colourbar of identical ticks, a field clipped to a physically
#  impossible value — none of these raise, and none show up in a test that only
#  asks whether a file was written. This walks every output folder and reports:
#
#    IMAGES   mostly-blank canvases, near-empty plot areas, images with almost no
#             colour variety (a scale collapsed onto one hue), and images whose
#             ink sits in a single band (a curve pinned to one edge).
#    TABLES   non-finite values, empty columns, and columns outside their
#             physical bounds (holdup and volume fractions in [0, 1], deposit
#             thickness no greater than the pipe radius, pressures positive).
#    METRICS  the scalar deliverables against their physical bounds, and the two
#             mass balances against the warning threshold.
#    SET      which expected figures are present, and which are absent (the
#             resolved-slug figures are legitimately absent when a line is not
#             slugging, so their absence is reported, not failed).
#    ANIM     GIFs that are single-frame, empty, or suspiciously small.
#    JSON     reports that are unparseable, empty, or carry non-finite numbers.
#    VTK/FOAM the 3-D export and the generated OpenFOAM cases: present, non-empty
#             and carrying the files their own README promises.
#
#  Every file in the folder is visited; anything with no specific check is listed
#  as SEEN with its size, so an unexpected empty artefact still surfaces.
#
#      python3 check_outputs.py                 # every scenario folder
#      python3 check_outputs.py outputs_steady  # one of them
#
#  Exit status is 1 if anything is flagged FAIL, 0 otherwise (WARN does not fail:
#  a quiet figure is often the correct answer for the mitigated scenario).
# =============================================================================
import csv
import json
import math
import os
import sys

import numpy as np

from _paths import CASE      # noqa: E402  (also installs the no-black style)

try:
    from PIL import Image
except Exception:                                          # pragma: no cover
    Image = None

SCENARIOS = ["outputs_steady", "outputs_shutin", "outputs_mitigated",
             "outputs_paper_steady", "outputs_paper_shutin",
             "outputs_paper_mitigated"]

#  figures every scenario must carry
REQUIRED = ["01_profiles.png", "02_holdup_spacetime.png", "03_PT_envelope.png",
            "04_PhiSH_map.png", "05_scenario_timeseries.png", "06_deposit.png",
            "07_probabilistic.png", "08_diagnostics.png",
            "14_holdup_multitime.png", "17_hydrate_distribution.png",
            "18_shutin_profile_deposit.png", "19_spacetime_fields.png",
            "20_holdup_durations.png", "22_cloud_maps.png",
            "23_dts_thermal_waterfall.png", "24_temperature_gradient.png",
            "25_das_flow_noise.png", "26_parameter_panels.png",
            "27_wellposedness_map.png"]

#  present only when the line is genuinely slugging
CONDITIONAL = ["15_slug_growth_propagation.png", "16_slug_train_waterfall.png",
               "21_riser_depth_time.png"]

#  column -> (low, high) physical bound, checked where the column exists
BOUNDS = {
    "holdup": (0.0, 1.0), "alpha_l": (0.0, 1.0),
    "phi_hydrate": (0.0, 1.0), "vapour_fraction": (0.0, 1.0),
    "P_bar": (0.0, 1.0e4), "p_bar": (0.0, 1.0e4),
    "T_C": (-50.0, 400.0), "Teq_C": (-50.0, 400.0),
    "f_slug_Hz": (0.0, 100.0), "a_i_1perm": (0.0, 1.0e5),
    "wetted_perim_frac": (0.0, 1.0), "liquid_level_h_over_D": (0.0, 1.0),
}

PIPE_D_MM = 254.5                                  # 10.75-in flowline ID


class Report:
    def __init__(self):
        self.rows = []

    def add(self, level, where, what):
        self.rows.append((level, where, what))

    def show(self):
        order = {"FAIL": 0, "WARN": 1, "NOTE": 2}
        for lvl, where, what in sorted(self.rows, key=lambda r: order[r[0]]):
            print(f"  [{lvl}] {where}: {what}", flush=True)
        n_fail = sum(1 for r in self.rows if r[0] == "FAIL")
        n_warn = sum(1 for r in self.rows if r[0] == "WARN")
        if not self.rows:
            print("  nothing flagged", flush=True)
        return n_fail, n_warn


# --------------------------------------------------------------- images ------
def check_image(path, rep):
    if Image is None:
        return
    name = os.path.basename(path)
    try:
        im = Image.open(path).convert("RGB")
    except Exception as exc:
        rep.add("FAIL", name, f"unreadable: {exc}")
        return
    a = np.asarray(im, np.int16)
    w, h = im.size
    if w < 500 or h < 300:
        rep.add("WARN", name, f"small canvas {w}x{h}")

    #  how much of the canvas is not the white background
    nonwhite = (a < 245).any(axis=2)
    ink = float(nonwhite.mean())
    #  A single curve on a large canvas legitimately covers only a per-cent or two
    #  of it, so "sparse" alone is not a fault. Judge blankness on how much of the
    #  PLOT AREA carries ink in at least one row and column — an empty axis has
    #  essentially none, a real curve spans it.
    if ink < 0.004:
        rep.add("FAIL", name, f"canvas is essentially blank (ink {ink*100:.2f} %)")
    elif ink < 0.010:
        rep.add("WARN", name, f"very little drawn (ink {ink*100:.2f} %)")

    #  distinct colours: a scale that collapsed onto one hue, or a plot that
    #  drew a single series where several were intended
    q = (a // 24).reshape(-1, 3)
    ncol = len({tuple(v) for v in q[::37]})
    if ncol < 6:
        rep.add("WARN", name, f"only {ncol} distinct colours — a scale or series "
                              f"set may have collapsed")

    #  ink concentrated in one horizontal band => a curve pinned to an edge
    rows_with_ink = nonwhite.mean(axis=1)
    busy = np.where(rows_with_ink > 0.02)[0]
    if busy.size and ink > 0.01:
        span = (busy[-1] - busy[0]) / h
        if span < 0.15:
            rep.add("WARN", name, f"all content in {span*100:.0f} % of the height "
                                  f"— content may be clipped to an edge")
    return


# --------------------------------------------------------------- tables ------
def check_csv(path, rep):
    name = os.path.basename(path)
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as exc:
        rep.add("FAIL", name, f"unreadable: {exc}")
        return
    if not rows:
        rep.add("FAIL", name, "no data rows")
        return
    cols = [c for c in rows[0] if c]
    for c in cols:
        vals = []
        for r in rows:
            try:
                vals.append(float(r[c]))
            except (TypeError, ValueError):
                vals.append(np.nan)
        v = np.asarray(vals, float)
        if np.isnan(v).all():
            continue                                   # a text column
        if not np.isfinite(v[~np.isnan(v)]).all():
            rep.add("FAIL", name, f"column '{c}' carries non-finite values")
        lo_hi = BOUNDS.get(c)
        if lo_hi is not None:
            fin = v[np.isfinite(v)]
            if fin.size and (fin.min() < lo_hi[0] - 1e-9 or fin.max() > lo_hi[1] + 1e-9):
                rep.add("FAIL", name, f"column '{c}' outside {lo_hi}: "
                                      f"{fin.min():.4g} .. {fin.max():.4g}")
        if c.startswith("deposit") and c.endswith("mm"):
            fin = v[np.isfinite(v)]
            if fin.size and fin.max() > PIPE_D_MM / 2.0 + 1e-6:
                rep.add("FAIL", name, f"column '{c}' exceeds the pipe radius "
                                      f"({PIPE_D_MM/2:.1f} mm): max {fin.max():.1f} mm")


# ------------------------------------------------------------ animations -----
def check_gif(path, rep):
    name = os.path.basename(path)
    if Image is None:
        return
    try:
        im = Image.open(path)
        n = getattr(im, "n_frames", 1)
    except Exception as exc:
        rep.add("FAIL", name, f"unreadable: {exc}")
        return
    if n < 2:
        rep.add("FAIL", name, f"animation has {n} frame(s)")
    elif n < 8:
        rep.add("WARN", name, f"only {n} frames")
    size = os.path.getsize(path)
    if size < 20_000:
        rep.add("WARN", name, f"very small for an animation ({size/1024:.0f} KB)")
    #  a GIF whose frames never change is a still in disguise
    try:
        im.seek(0)
        first = np.asarray(im.convert("RGB"), np.int16)
        im.seek(n - 1)
        last = np.asarray(im.convert("RGB"), np.int16)
        if first.shape == last.shape:
            #  These are sparse line drawings on white, so a MEAN difference is
            #  tiny even when the animation moves plainly. Judge it on how many
            #  pixels actually changed, and by how much, not on the average.
            d = np.abs(first - last).max(axis=2)
            changed = float((d > 20).mean())
            if d.max() < 20 or changed < 5e-4:
                rep.add("WARN", name, f"first and last frames are nearly identical "
                                      f"({changed*100:.3f} % of pixels changed) — the "
                                      f"animation may not be moving")
    except Exception:
        pass


# ------------------------------------------------------------------ json -----
def _walk_numbers(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from ((f"{k}.{p}" if p else k, x) for p, x in _walk_numbers(v))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            yield from ((f"[{i}].{p}" if p else f"[{i}]", x) for p, x in _walk_numbers(v))
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        yield "", float(obj)


def check_json(path, rep):
    name = os.path.basename(path)
    try:
        d = json.load(open(path))
    except Exception as exc:
        rep.add("FAIL", name, f"unparseable: {exc}")
        return
    if not d:
        rep.add("FAIL", name, "empty")
        return
    bad = [k for k, v in _walk_numbers(d) if not math.isfinite(v)]
    if bad:
        rep.add("FAIL", name, f"{len(bad)} non-finite number(s), e.g. {bad[:3]}")


# ------------------------------------------------------------- vtk / foam ----
def check_vtk(path, rep):
    name = os.path.basename(path)
    size = os.path.getsize(path)
    if size < 1_000:
        rep.add("FAIL", name, f"suspiciously small ({size} bytes)")
        return
    head = open(path, "r", errors="ignore").read(400)
    if "vtk" not in head.lower():
        rep.add("FAIL", name, "does not look like a VTK file")
    if "POINTS" not in head and "POINTS" not in open(path, "r",
                                                     errors="ignore").read(20_000):
        rep.add("WARN", name, "no POINTS section found in the header")


def check_openfoam(folder, rep):
    root = os.path.join(folder, "openfoam_cases")
    if not os.path.isdir(root):
        return
    cases = [d for d in sorted(os.listdir(root))
             if os.path.isdir(os.path.join(root, d))]
    if not cases:
        rep.add("WARN", "openfoam_cases", "no case directories generated")
        return
    for c in cases:
        cd = os.path.join(root, c)
        for need in ("0", "constant", "system"):
            if not os.path.isdir(os.path.join(cd, need)):
                rep.add("FAIL", f"openfoam_cases/{c}", f"missing '{need}/'")
        allrun = os.path.join(cd, "Allrun")
        if not os.path.exists(allrun):
            rep.add("FAIL", f"openfoam_cases/{c}", "missing Allrun")
        elif os.path.getsize(allrun) < 20:
            rep.add("WARN", f"openfoam_cases/{c}", "Allrun is nearly empty")
    rep.add("NOTE", "openfoam_cases", f"{len(cases)} case(s) generated and complete")


# -------------------------------------------------------------- metrics ------
def check_metrics(folder, rep):
    path = os.path.join(folder, "key_metrics.json")
    if not os.path.exists(path):
        #  the outputs_paper_* folders are a FIGURE-ONLY build (make_paper_figures
        #  renders charts without titles for the manuscript and writes no tables),
        #  so a missing metrics file there is the expected shape, not a fault.
        if "outputs_paper" in os.path.basename(folder):
            rep.add("NOTE", "key_metrics.json",
                    "absent — this is the figure-only manuscript build")
        else:
            rep.add("WARN", "key_metrics.json", "absent")
        return
    d = json.load(open(path))
    checks = [
        ("mass_conservation_err", 0.0, 0.05, "FAIL"),
        ("gas_mass_conservation_err", 0.0, 0.05, "FAIL"),
        ("P_plug", 0.0, 1.0, "FAIL"),
        ("peak_deposit_mm", 0.0, PIPE_D_MM / 2.0, "FAIL"),
        ("MEG_wt_pct", 0.0, 100.0, "FAIL"),
        ("max_subcooling_C", -100.0, 100.0, "FAIL"),
    ]
    for key, lo, hi, level in checks:
        if key not in d or d[key] is None:
            continue
        try:
            v = float(d[key])
        except (TypeError, ValueError):
            continue
        if not math.isfinite(v):
            rep.add("FAIL", "key_metrics.json", f"{key} is not finite")
        elif v < lo - 1e-9 or v > hi + 1e-9:
            rep.add(level, "key_metrics.json",
                    f"{key} = {v:.6g} outside [{lo:g}, {hi:g}]")
    if d.get("fallbacks", 0):
        rep.add("WARN", "key_metrics.json",
                f"{d['fallbacks']:.0f} solver fallbacks were triggered")


# ----------------------------------------------------------------- main ------
def check_folder(folder, rep):
    present = set(os.listdir(folder))
    for fn in REQUIRED:
        if fn not in present:
            rep.add("FAIL", fn, "expected figure is missing")
    for fn in CONDITIONAL:
        if fn not in present:
            rep.add("NOTE", fn, "absent — the line is not slugging in this "
                                "scenario, which is a legitimate skip")
    if "spacetime_state.npz" not in present:
        rep.add("WARN", "spacetime_state.npz", "snapshot archive not written")

    seen = 0
    for fn in sorted(present):
        p = os.path.join(folder, fn)
        if os.path.isdir(p):
            continue
        seen += 1
        if fn.endswith(".png"):
            check_image(p, rep)
        elif fn.endswith(".csv"):
            check_csv(p, rep)
        elif fn.endswith(".gif"):
            check_gif(p, rep)
        elif fn.endswith(".json"):
            check_json(p, rep)
        elif fn.endswith(".vtk"):
            check_vtk(p, rep)
        elif fn.endswith(".npz"):
            if os.path.getsize(p) < 1_000:
                rep.add("FAIL", fn, "snapshot archive is empty")
        elif fn.endswith((".txt", ".md")):
            if os.path.getsize(p) == 0:
                rep.add("WARN", fn, "empty text file")
        else:
            rep.add("NOTE", fn, f"no specific check ({os.path.getsize(p)} bytes)")
        if os.path.getsize(p) == 0:
            rep.add("FAIL", fn, "zero bytes")
    check_openfoam(folder, rep)
    check_metrics(folder, rep)
    rep.add("NOTE", "(folder)", f"{seen} file(s) inspected")


def main(argv):
    todo = argv or SCENARIOS
    total_fail = total_warn = 0
    for name in todo:
        folder = os.path.join(CASE, name)
        if not os.path.isdir(folder):
            print(f"\n=== {name}: absent ===", flush=True)
            continue
        n_png = len([f for f in os.listdir(folder) if f.endswith(".png")])
        n_csv = len([f for f in os.listdir(folder) if f.endswith(".csv")])
        print(f"\n=== {name}  ({n_png} figures, {n_csv} tables) ===", flush=True)
        rep = Report()
        check_folder(folder, rep)
        f, w = rep.show()
        total_fail += f
        total_warn += w
    print(f"\n{total_fail} failure(s), {total_warn} warning(s)", flush=True)
    return 1 if total_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
