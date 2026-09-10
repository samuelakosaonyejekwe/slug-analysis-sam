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
from _paths import CASE  # noqa: E402  (also installs the no-black style)

try:
    from PIL import Image
except Exception:                                          # pragma: no cover
    Image = None   # type: ignore[assignment]  # optional dependency; guarded at every use

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

def _pipe_d_mm(default=254.5):
    """Bore in mm, from whichever run wrote a case_config.json.

    It was a literal 254.5, which silently mis-bounds the deposit checks the moment
    the case geometry moves — the same drift the report's prose had. The default is
    kept for a folder that carries no config.
    """
    for scen in SCENARIOS:
        cfg = os.path.join(CASE, scen, "case_config.json")
        try:
            with open(cfg) as fh:
                return float(json.load(fh)["pipeline"]["diameter_m"]) * 1000.0
        except Exception:
            continue
    return default


PIPE_D_MM = _pipe_d_mm()                           # flowline ID, from the run's own config


def _line_km(default=32.0):
    """Route length in km, from whichever run wrote a case_config.json."""
    for scen in SCENARIOS:
        cfg = os.path.join(CASE, scen, "case_config.json")
        try:
            with open(cfg) as fh:
                return float(json.load(fh)["pipeline"]["length_m"]) / 1000.0
        except Exception:
            continue
    return default


def _line_liquid_m3(default=1.0e4):
    """Total volume of the line, in m3 — the ceiling on any surge volume.

    No slug, and no riser blowdown, can deliver more liquid than the pipe holds. Taken
    from the run's own config so it tracks a geometry change instead of drifting.
    """
    for scen in SCENARIOS:
        cfg = os.path.join(CASE, scen, "case_config.json")
        try:
            with open(cfg) as fh:
                c = json.load(fh)["pipeline"]
            return math.pi / 4.0 * float(c["diameter_m"]) ** 2 * float(c["length_m"])
        except Exception:
            continue
    return default


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
        #  convert() returns a NEW image, so the file handle can be released here
        #  rather than left to the garbage collector — this runs over every figure
        #  in every output folder, and the descriptors accumulate.
        with Image.open(path) as _src:
            im = _src.convert("RGB")
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
        numeric_text = 0            # cells that PARSE as a number, nan/inf included
        for r in rows:
            cell = (r[c] or "").strip()
            if cell.lower() in ("nan", "inf", "-inf", "+inf", "infinity", "-infinity", ""):
                numeric_text += 1
                vals.append(np.nan)
                continue
            try:
                vals.append(float(cell))
                numeric_text += 1
            except (TypeError, ValueError):
                vals.append(np.nan)
        v = np.asarray(vals, float)
        if np.isnan(v).all():
            #  An all-NaN column was treated as a text column and skipped outright. It is
            #  only a text column if its cells do not parse as numbers; if they DO -- a
            #  column of "nan" -- then it is a numeric column that is undefined
            #  everywhere, and skipping it silently means the gate never looks at it. On
            #  the mitigated scenario Phi_SH_sustained is 70/70 nan, which is correct (no
            #  hydrate forms, so a sustained coupling number does not exist there) but is
            #  worth saying rather than passing over.
            if numeric_text == len(rows) and rows:
                rep.add("NOTE", name, f"column '{c}' is undefined in every row "
                                      f"({len(rows)} of {len(rows)})")
            continue                                   # text, or wholly undefined
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
    #  The image is opened ONCE and held for the whole check. It used to be opened in a
    #  `with` block that read only the frame count, which closed it -- and every line
    #  below then ran against a closed file, so `im.seek(0)` raised "'NoneType' object has
    #  no attribute 'read'" on every animation. The surrounding handler swallowed that,
    #  so the folder reported clean and the frames-never-change check had in fact never
    #  run on a single GIF. Closed in the `finally` at the end instead.
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
    except Exception as exc:
        #  A gate that swallows its own errors passes everything it could not read.
        #  This one covered the whole GIF inspection, so an unreadable or malformed
        #  animation produced no row at all and the folder came back clean.
        rep.add("WARN", name, f"could not be inspected ({type(exc).__name__}: {exc}) "
                              f"— this check did not run on it")
    finally:
        try:
            im.close()
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
        with open(path) as fh:
            d = json.load(fh)
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
    with open(path, "r", errors="ignore") as fh:
        head = fh.read(400)
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
    with open(path) as fh:
        d = json.load(fh)
    checks = [
        ("mass_conservation_err", 0.0, 0.05, "FAIL"),
        ("gas_mass_conservation_err", 0.0, 0.05, "FAIL"),
        ("P_plug", 0.0, 1.0, "FAIL"),
        ("peak_deposit_mm", 0.0, PIPE_D_MM / 2.0, "FAIL"),
        ("MEG_wt_pct", 0.0, 100.0, "FAIL"),
        ("max_subcooling_C", -100.0, 100.0, "FAIL"),
        #  the erosional exceedance is a fraction of route and the length it implies;
        #  bounded here so a sign error or a unit slip surfaces rather than being
        #  read as a design number, which is exactly what the raw peak velocity was
        ("erosional_exceedance_frac", 0.0, 1.0, "FAIL"),
        ("erosional_exceedance_km", 0.0, _line_km(), "FAIL"),
        #  Bounds that would have caught real defects and did not exist to.
        #
        #  dew_point_bar was 1.0000000000000004 on all three scenarios -- the LOWER END of
        #  the EOS bisection bracket, returned because the root was never bracketed and the
        #  loop simply converged on where it started. Three different monitor temperatures
        #  giving one identical value is the signature. The saturation search now returns
        #  NaN (written as null, and skipped below) when the root is not in range, so a
        #  value AT the bracket edge is the thing to refuse.
        ("dew_point_bar", 1.0 + 1e-6, 700.0 - 1e-6, "FAIL"),
        #  a fraction, and a length that cannot exceed the route
        ("Phi_SH_above_critical_frac", 0.0, 1.0, "FAIL"),
        ("Phi_SH_supercritical_time_frac", 0.0, 1.0, "FAIL"),
        ("hydrate_scoured_frac", 0.0, 1.0, "FAIL"),
        ("hydrate_packing_clip_frac", 0.0, 1.0, "FAIL"),
        ("hydrate_outflow_frac", 0.0, 1.0, "FAIL"),
        ("sustained_supercritical_km", 0.0, _line_km(), "FAIL"),
        ("under_inhibited_km", 0.0, _line_km(), "FAIL"),
        #  a no-touch time is an elapsed time; it cannot be negative, and it cannot exceed
        #  the window it was measured in
        ("cooldown_to_hydrate_h", 0.0, 1.0e5, "FAIL"),
        #  a hydrate slurry is never LESS viscous than its carrier
        ("slurry_rel_viscosity", 1.0 - 1e-9, math.inf, "FAIL"),
        #  The slug-catcher duty is bounded BELOW by nothing useful but bounded ABOVE by the
        #  liquid the line can physically hold: no cycle can deliver more liquid than exists
        #  in the pipe. The riser-inventory basis is a sum over cells, so a unit slip or a
        #  mask that ran away would show up here and nowhere else.
        ("V_surge_P90_m3", 0.0, _line_liquid_m3(), "FAIL"),
        ("V_riser_liquid_m3", 0.0, _line_liquid_m3(), "FAIL"),
        ("V_surge_hydrodynamic_m3", 0.0, _line_liquid_m3(), "FAIL"),
        #  an inclination in degrees, of an ascent that exists
        ("riser_incline_deg", 0.0, 90.0, "FAIL"),
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
    #  The Camargo-Palermo relative viscosity DIVERGES at the packing limit, so once
    #  phi_peak reaches phi_max the number is set by the 0.999 clip and by nothing else
    #  (3.16e7 = 0.001**-2.5 on the shut-in case). It is the honest output of the
    #  correlation at that packing, but it is not a resolved magnitude, and it should
    #  never be quoted as one.
    #  Phi_SH divides by shear removal, so at zero flow it diverges: the shut-in reports
    #  1467 against Phi_crit 1.08. Honest output of the definition, meaningless magnitude.
    if d.get("Phi_SH_shear_limited"):
        rep.add("WARN", "key_metrics.json",
                f"max_Phi_SH = {float(d.get('max_Phi_SH', float('nan'))):.4g} is "
                f"SHEAR-LIMITED — the line is at rest, so the removal term in Phi_SH goes "
                f"to zero and the magnitude is set by that, not by resolved physics; quote "
                f"the super-critical EXTENT, not this number")
    if d.get("slurry_visc_saturated"):
        rep.add("WARN", "key_metrics.json",
                f"slurry_rel_viscosity = {float(d.get('slurry_rel_viscosity', float('nan'))):.3g} "
                f"is SATURATED at the packing limit — set by the 0.999 clip, not resolved")


def check_fluid_identity(folder, rep):
    """The EOS composition and the flow model must describe the SAME fluid.

    They now do, and this check is what proves it every run rather than a claim in a
    README. The gap used to be 56 %: the composition flashed to a ~549 kg/m3 liquid at
    line conditions -- a light volatile oil -- while the flow model ran rho_oil = 858
    kg/m3, a medium crude, and every velocity, holdup and pressure drop came from the
    latter. The reason was not the volume shift but the pseudo-component: C7+ carried
    n-heptane's own constants (Tc 540.20 K, Pc 27.40 bar, w 0.3495, MW 100 g/mol), so no
    Peneloux shift inside its physical range could reach 858 -- it would have needed
    s_C7+ ~ 1.33 against about +-0.2, a correction larger than the co-volume.

    Characterising C7+ properly (Riazi-Daubert T_b, then Kesler-Lee for Tc/Pc/w, MW 250
    g/mol) moves the flash onto a medium crude on its own, and s_C7+ = 0.2253 -- inside
    the Jhaveri-Youngren C7+ range of 0.1-0.3 -- then reproduces the case's own stated
    858 kg/m3. The check stays, at WARN, because the two descriptions are still
    independently specified and a later edit to either could separate them again.

    NOTE what this does NOT check. It compares the EOS liquid against rho_oil, the OIL
    alone. The flow model's LIQUID is heavier than both because it carries the water cut,
    which is correct and not a mismatch. The phase split is also still two descriptions:
    the EOS holds the hydrocarbon single-phase to its bubble point while the flow model
    runs a fixed inlet gas rate. compo_pvt.png states that on its face.
    """
    cfg_p = os.path.join(folder, "case_config.json")
    if not os.path.exists(cfg_p):
        return
    try:
        with open(cfg_p) as fh:
            cfg = json.load(fh)
        comp = cfg.get("fluids", {}).get("composition")
        rho_flow = float(cfg["fluids"]["rho_oil"])
        if not comp:
            return
        sys.path.insert(0, os.path.dirname(CASE))
        import shct_eos
        rho_eos = float(shct_eos.eos_properties(120.0, 20.0, comp)["rho_oil"])
        gap = 100.0 * (rho_flow - rho_eos) / max(rho_eos, 1e-9)
        if abs(gap) > 10.0:
            rep.add("WARN", "case_config.json",
                    f"FLUID IDENTITY: the EOS composition gives a {rho_eos:.0f} kg/m3 liquid "
                    f"at 120 bar/20 C while the flow model runs rho_oil = {rho_flow:.0f} "
                    f"kg/m3 — {gap:+.0f} %. These are different fluids; the compositional "
                    f"figures describe EOS phase behaviour, NOT the pipeline's medium crude")
    except Exception:
        return


def check_folder_freshness(folder, rep, tol_h=6.0):
    """Flag files a regeneration LEFT BEHIND.

    A generated folder is rebuilt as a set, so its files should share a timestamp. When a
    generator silently stops producing one of them, that file simply persists — and nothing
    notices, because every check that exists looks at the files that ARE written. That is
    how outputs_paper_steady/verif_*.png and threed_*.png sat two days stale while
    export_paper_figures.py copied them into the manuscript set with a fresh mtime, so the
    PDF-freshness check passed on stale content. Compare each file with the NEWEST file
    beside it instead: anything left far behind was not rewritten by the last run.
    """
    try:
        files = [os.path.join(folder, f) for f in os.listdir(folder)]
        files = [f for f in files if os.path.isfile(f)]
    except OSError:
        return
    if len(files) < 3:
        return
    newest = max(os.path.getmtime(f) for f in files)
    stale = [f for f in files if (newest - os.path.getmtime(f)) > tol_h * 3600.0]
    if stale:
        names = ", ".join(sorted(os.path.basename(f) for f in stale)[:6])
        more = "" if len(stale) <= 6 else f" (+{len(stale) - 6} more)"
        rep.add("WARN", "(folder)",
                f"{len(stale)} file(s) more than {tol_h:g} h older than the newest file "
                f"here — the last regeneration did not rewrite them: {names}{more}")


# ----------------------------------------------------------------- main ------
def check_folder(folder, rep):
    check_folder_freshness(folder, rep)
    check_fluid_identity(folder, rep)
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
