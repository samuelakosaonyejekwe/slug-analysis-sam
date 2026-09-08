#!/usr/bin/env python3
# =============================================================================
#  shct_openfoam.py — SHCT  <->  OpenFOAM (open-source CFD) COUPLING
# -----------------------------------------------------------------------------
#  Couples the fast 1-D SHCT solver to full 3-D CFD (OpenFOAM) for ONLY the
#  sections that need it. Workflow:
#     1. identify_critical_sections(sv): from the solved SHCT fields, pick the
#        sections where 3-D resolution actually matters — riser base / steep
#        terrain, severe (intermittent) slugging, and the hydrate-critical
#        (Phi_SH>1, high subcooling, max deposit) cells.
#     2. write_case(section): emit a COMPLETE, runnable OpenFOAM `interFoam`
#        (volume-of-fluid, two-phase) case for that pipe segment, with the mesh
#        (cylindrical o-grid blockMesh, gravity tilted by the section inclination)
#        and the inlet/outlet boundary conditions taken FROM the SHCT solution
#        (mixture velocity, liquid fraction = holdup, outlet pressure, stratified
#        initial condition from the interface level h/D).
#     3. run_case(): if an OpenFOAM installation is detected (blockMesh on PATH),
#        run blockMesh + setFields + interFoam; otherwise the case is written for
#        the user to run on an OpenFOAM machine (./Allrun).
#     4. ingest_results(): read the CFD result back (volume-averaged liquid
#        fraction, etc.) to refine / validate the SHCT prediction for that section.
#
#  COUPLING TYPE: one-way SHCT -> CFD (boundary conditions) by default, with a
#  feedback hook (CFD -> SHCT effective closure) for two-way iteration.
#
#  HONEST SCOPE: this module GENERATES correct, standard interFoam cases and runs
#  them when OpenFOAM is available. It does not bundle OpenFOAM. The 3-D physics
#  is then OpenFOAM's (a genuine 3-D Navier-Stokes VOF solve on the section);
#  SHCT supplies the whole-line context and the per-section boundary conditions.
#
#  WHAT THE COMPARISON CAN AND CANNOT SHOW. Until 2026-09-08 this module had never
#  been run against a real OpenFOAM install -- the tests wrote cases and asserted on
#  the files, which is why several defects survived in code that looked correct. Run
#  for real (OpenFOAM v2406, as-operated case, 20 000-cell o-grid, 2 s), it establishes:
#
#    * In the default inlet_mode="holdup" the SHCT liquid fraction is imposed at the
#      inlet AND used as the initial condition, so the CFD returns it. Measured:
#      imposing 0.20 / 0.3637 / 0.60 on the same section gives 0.2015 / 0.3637 /
#      0.5988. The manifest's rel_diff_pct is therefore a consistency check on the
#      case writer, not a validation of the holdup closure, and couple_iterate driven
#      by a real run converges at iteration 0 having changed nothing. inlet_mode=
#      "noslip" injects the volumetric split lambda_l = Vsl/(Vsl+Vsg) instead and
#      leaves the holdup for the CFD to predict.
#    * The generated case IS a valid interFoam case: blockMesh reports max non-
#      orthogonality 31.5, max skewness 0.97, "Mesh OK", and setFields selects
#      7120/20000 cells against a target liquid fraction of 0.3637 -- 0.356, which is
#      the cell-centre discretisation of it.
#    * Three things that made the selected sections unfit to test the slip closure have
#      since been corrected. The mixture reads single-phase at the states it picks (EOS
#      flash gives V = 0), and the PVT table's gas column returned the LIQUID root there,
#      so cases were written with a density ratio of 1.7-2.0; build_pvt_table now fills
#      those nodes from the vapour the fluid actually releases and the ratio is ~8. The
#      riser-base velocity was pinned on the 1-D solver's own momentum clip; the clip has
#      been widened to a bound that does not bind. And the case was written laminar at
#      Re = 1e5-2.6e5; it now carries k-omega SST above RE_TURBULENT.
#    * What remains open is DEVELOPMENT LENGTH, and it is a property of the segment, not
#      a defect. Injecting lambda_l = 0.2537 where the 1-D closure predicts alpha_l =
#      0.3637 (C0 ~ 1.17), interFoam returns a holdup of 0.2546 laminar and 0.2554 with
#      k-omega SST at 12 diameters, and 0.2561 at 50 -- C0 ~ 1.01. Twelve diameters fed
#      from a uniform inlet profile cannot build the velocity/void correlation that C0
#      IS, whatever closure is switched on. Reading these runs as a refutation of the
#      drift-flux closure would be reading the segment length.
# =============================================================================
from __future__ import annotations

import copy
import json
import math
import os
import shutil
import subprocess
from typing import Any, Callable, Optional

import numpy as np

gas_density: Optional[Callable[..., Any]]
try:
    from shct_correlations import gas_density as _gas_density
    gas_density = _gas_density
except Exception:                                       # pragma: no cover
    gas_density = None
import shct_crosssection as cx

G = 9.81


# ---------------------------------------------------------------------------
#  FoamFile dictionary header
# ---------------------------------------------------------------------------
def _foam(cls, obj, location):
    return ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
            "| SHCT->OpenFOAM coupling — auto-generated case                              |\n"
            "\\*---------------------------------------------------------------------------*/\n"
            "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            f"    class       {cls};\n    location    \"{location}\";\n    object      {obj};\n"
            "}\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")


# ---------------------------------------------------------------------------
#  Cylindrical o-grid blockMeshDict (pipe segment, axis = z)
# ---------------------------------------------------------------------------
def _blockmeshdict(R, L, Ni=10, Nz=40, box=0.5):
    oc = R / math.sqrt(2.0)
    a = box * oc
    V = []
    for z in (0.0, L):
        V += [(-a, -a, z), (a, -a, z), (a, a, z), (-a, a, z),
              (-oc, -oc, z), (oc, -oc, z), (oc, oc, z), (-oc, oc, z)]
    verts = "\n".join(f"    ({x:.6g} {y:.6g} {z:.6g})" for (x, y, z) in V)
    blocks = [
        (0, 1, 2, 3, 8, 9, 10, 11),       # central
        (4, 5, 1, 0, 12, 13, 9, 8),       # bottom
        (5, 6, 2, 1, 13, 14, 10, 9),      # right
        (6, 7, 3, 2, 14, 15, 11, 10),     # top
        (7, 4, 0, 3, 15, 12, 8, 11),      # left
    ]
    blk = "\n".join(f"    hex ({' '.join(map(str, b))}) ({Ni} {Ni} {Nz}) simpleGrading (1 1 1)"
                    for b in blocks)
    arcs = [
        (4, 5, (0.0, -R, 0.0)), (5, 6, (R, 0.0, 0.0)), (6, 7, (0.0, R, 0.0)), (7, 4, (-R, 0.0, 0.0)),
        (12, 13, (0.0, -R, L)), (13, 14, (R, 0.0, L)), (14, 15, (0.0, R, L)), (15, 12, (-R, 0.0, L)),
    ]
    edg = "\n".join(f"    arc {i} {j} ({p[0]:.6g} {p[1]:.6g} {p[2]:.6g})" for (i, j, p) in arcs)
    inlet = ["(0 1 2 3)", "(4 5 1 0)", "(5 6 2 1)", "(6 7 3 2)", "(7 4 0 3)"]
    outlet = ["(8 9 10 11)", "(12 13 9 8)", "(13 14 10 9)", "(14 15 11 10)", "(15 12 8 11)"]
    walls = ["(4 5 13 12)", "(5 6 14 13)", "(6 7 15 14)", "(7 4 12 15)"]
    def _patch(name, typ, faces):
        f = "\n".join(f"            {x}" for x in faces)
        return (f"    {name}\n    {{\n        type {typ};\n"
                f"        faces\n        (\n{f}\n        );\n    }}\n")
    bnd = (_patch("inlet", "patch", inlet) + _patch("outlet", "patch", outlet)
           + _patch("walls", "wall", walls))
    return (_foam("dictionary", "blockMeshDict", "system")
            + "scale 1;\n\nvertices\n(\n" + verts + "\n);\n\n"
            + "blocks\n(\n" + blk + "\n);\n\n"
            + "edges\n(\n" + edg + "\n);\n\n"
            + "boundary\n(\n" + bnd + ");\n\n"
            + "mergePatchPairs ();\n")


# ---------------------------------------------------------------------------
#  Identify the sections that need CFD
# ---------------------------------------------------------------------------
def identify_critical_sections(sv, max_sections=3, seg_len_factor=12.0):
    """Score each cell for 'needs 3-D CFD' and return up to max_sections distinct
    sections. Score rewards: hydrate criticality (Phi_SH>1), steep terrain / riser,
    intermittent (slug/churn) regime, high subcooling and deposit."""
    r = sv.results
    def med(A):
        return np.nanmedian(A, 1)
    x = sv.x
    D = med(r["D"]) if "D" in r else np.full_like(x, sv.case.pipeline.diameter_m)
    phish = med(r["max_PhiSH"]); sub = med(r["Tsub"]); delta = med(r["delta"])
    alpha = med(r["alpha_l"]); regime = med(r["regime"])
    theta = np.abs(sv.theta)
    # normalised scores
    def nz(a):
        a = np.asarray(a, float); rng = np.nanmax(a) - np.nanmin(a)
        return (a - np.nanmin(a)) / rng if rng > 1e-12 else np.zeros_like(a)
    intermittent = np.isin(np.round(regime), [2, 5]).astype(float)
    score = (1.2 * nz(np.clip(phish, 0, None)) + 1.0 * nz(theta) + 0.8 * intermittent
             + 0.7 * nz(np.clip(sub, 0, None)) + 0.9 * nz(delta))
    # greedy pick of separated peaks
    picked: list[int] = []
    order = np.argsort(score)[::-1]
    min_sep = max(int(0.04 * len(x)), 3)
    for i in order:
        if all(abs(i - j) > min_sep for j in picked):
            picked.append(int(i))
        if len(picked) >= max_sections:
            break
    picked = sorted(picked)
    rho_l = float(sv.rho_l)
    #  the momentum clip the 1-D solver applies, read from solver rather than repeated
    try:
        from solver import UM_CLIP_HI, UM_CLIP_LO
    except Exception:                                    # pragma: no cover
        UM_CLIP_LO, UM_CLIP_HI = -2.0, 30.0
    j_all = np.asarray(r["j"], float)
    #  NO-SLIP (input) liquid fraction lambda_l = Vsl / (Vsl + Vsg), from the phase
    #  velocities. Unlike alpha_l this is the volumetric split actually being INJECTED;
    #  the difference alpha_l - lambda_l is the slip the drift-flux closure predicts,
    #  and it is the only part of the holdup a CFD run can independently confirm.
    lam = None
    try:
        h_s = np.asarray(r["snap_holdup"], float)
        vl_s = np.asarray(r["snap_vl"], float)
        vg_s = np.asarray(r["snap_vg"], float)
        if h_s.shape == vl_s.shape == vg_s.shape and h_s.ndim == 2 and h_s.shape[1] == x.size:
            Vsl = np.nanmedian(h_s * vl_s, axis=0)
            Vsg = np.nanmedian((1.0 - h_s) * vg_s, axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                lam = np.clip(Vsl / (Vsl + Vsg), 0.0, 1.0)
    except (KeyError, ValueError, TypeError):               # pragma: no cover
        lam = None
    secs = []
    for rank, i in enumerate(picked):
        p_bar = float(med(r["p"])[i]); T_C = float(med(r["T"])[i])
        rho_g = (float(gas_density(p_bar, T_C, sv.case.fluids))
                 if gas_density is not None else 30.0)
        Di = float(D[i])
        #  SIGNED mixture velocity. This was `max(med(r["j"])[i], 0.05)`, which mapped a
        #  section flowing BACKWARDS onto a trickle forwards. It is not a corner case: at
        #  the riser base of the as-operated run the median j is -2.00 m/s, so the CFD
        #  case for the single section this module scores highest -- steep terrain, the
        #  reason the coupling exists -- was written with an inlet of +0.05 m/s, reversed
        #  in direction and 40x too small. Reversing the flow is instead expressed by
        #  reversing the SEGMENT: running -|j| along +z up a slope +theta is the same
        #  physics as running +|j| along +z down a slope -theta, and interFoam's inlet /
        #  outlet patches then sit at the ends the fluid actually enters and leaves.
        j_i = float(med(r["j"])[i])
        reversed_flow = j_i < 0.0
        theta_i = float(sv.theta[i])
        Vmix = max(abs(j_i), 0.05)
        #  Flag a boundary condition taken from a cell pinned on the solver's velocity
        #  clip. The state there is a numerical guard, not a converged physical velocity,
        #  so a CFD-vs-SHCT holdup difference for such a section measures the clip.
        a_i_l = float(np.clip(alpha[i], 1e-3, 0.999))
        rho_m = a_i_l * rho_l + (1.0 - a_i_l) * rho_g
        mu_m = (a_i_l * float(sv.case.fluids.mu_liquid)
                + (1.0 - a_i_l) * float(sv.case.fluids.mu_gas))
        Re_i = float(rho_m * Vmix * Di / max(mu_m, 1e-9))
        col = j_all[i] if j_all.ndim > 1 else np.array([j_i])
        clip_frac = float(np.mean(np.isclose(col, UM_CLIP_LO)
                                  | np.isclose(col, UM_CLIP_HI)))
        secs.append({
            "name": f"section_{rank+1}_x{ x[i]/1000:.1f}km".replace(".", "p"),
            "index": i, "x_km": float(x[i] / 1000.0), "length_m": float(seg_len_factor * Di),
            "D": Di, "theta_rad": -theta_i if reversed_flow else theta_i,
            "Vmix": Vmix, "alpha_l": float(np.clip(alpha[i], 1e-3, 0.999)),
            "lambda_l": (float(lam[i]) if lam is not None and np.isfinite(lam[i])
                         else float(np.clip(alpha[i], 1e-3, 0.999))),
            "Vmix_shct_signed": j_i, "flow_reversed": bool(reversed_flow),
            "theta_rad_shct": theta_i,
            "velocity_clip_fraction": clip_frac,
            "bc_from_clipped_state": bool(clip_frac > 0.5),
            "p_bar": p_bar, "T_C": T_C, "Phi_SH": float(phish[i]), "subcooling_C": float(sub[i]),
            "deposit_mm": float(delta[i] * 1000.0), "regime": int(round(regime[i])),
            "rho_l": rho_l, "rho_g": rho_g,
            "rho_ratio": float(rho_l / max(rho_g, 1e-6)),
            #  Mixture Reynolds number of the segment. The case is written laminar (see
            #  constant/turbulenceProperties), which for these sections means solving a
            #  Re ~ 1.5e5 flow with no turbulence model at all. That is not a detail: the
            #  drift-flux distribution parameter C0, which is where nearly all of SHCT's
            #  predicted slip comes from, IS the correlation between the turbulent
            #  velocity profile and the void profile. Measured, injecting lambda_l =
            #  0.2537 where SHCT predicts alpha_l = 0.3637 (C0 ~ 1.17): interFoam returns
            #  0.2546 at 12 diameters and 0.2553 at 50 diameters -- C0 ~ 1.00, i.e. no
            #  profile slip forms, at either length and at either density ratio. A
            #  laminar solve with a flat inlet profile cannot produce C0, so these cases
            #  cannot yet adjudicate that closure however long they are run.
            "Re_mixture": Re_i,
            #  interFoam is a TWO-phase VOF solver; it needs two phases that differ. When
            #  the mixture is single-phase at the section state the PVT table has no
            #  vapour to report and its gas column returns the liquid/dense root, so the
            #  case is written with a "gas" barely lighter than the liquid. Measured on
            #  the as-operated run: at all three selected sections (115.2/106.2/95.2 bar,
            #  9.5/8.0/7.4 C) the EOS flash gives vapour fraction V = 0.0000 and the table
            #  returns rho_g = 560/560/499 kg/m3 against rho_l = 975 -- density ratios of
            #  1.7-2.0 where a real gas-liquid system at these pressures is nearer 10.
            #  The case still RUNS and its numbers still look plausible, which is why this
            #  survived until the coupling was exercised against real OpenFOAM. Flagged
            #  rather than silently corrected: which density is right is a question about
            #  the fluid model, not about this module.
            "phases_distinct": bool(rho_l / max(rho_g, 1e-6) > 4.0),
            "mu_l": float(sv.case.fluids.mu_liquid), "mu_g": float(sv.case.fluids.mu_gas),
            "sigma": float(sv.case.fluids.sigma),
            "reason": _why(phish[i], theta[i], intermittent[i], sub[i], delta[i])})
    return secs


def _why(phish, theta, inter, sub, delta):
    bits = []
    if phish > 1:
        bits.append("Phi_SH>1 (hydrate-critical)")
    if theta > 0.15:
        bits.append("steep terrain / riser")
    if inter > 0.5:
        bits.append("intermittent (slug/churn)")
    if sub > 0:
        bits.append("subcooled")
    if delta > 1e-4:
        bits.append("wall deposit")
    return ", ".join(bits) or "elevated coupled score"


# ---------------------------------------------------------------------------
#  Write a complete interFoam case for one section
# ---------------------------------------------------------------------------
#  What the inlet imposes, and therefore what the run can independently tell you.
#
#  MEASURED, on real interFoam v2406 runs of the as-operated section at x = 23.1 km
#  (20 000 cells, 2 s, three runs identical but for the imposed liquid fraction):
#
#      imposed alpha_l   0.2000   0.3637   0.6000
#      CFD volume mean   0.2015   0.3637   0.5988
#
#  The CFD returns whatever it is given, to within its own +-0.3 % oscillation, over a
#  three-fold range. That is not agreement — with alpha_l fixed at the inlet AND used as
#  the initial condition, a 3 m segment flushed 1.7 times in 2 s can do nothing else. So
#  the manifest's SHCT-vs-CFD difference in "holdup" mode is a consistency check on the
#  case writer, NOT a validation of the holdup closure, and couple_iterate driven by a
#  real run converges at iteration 0 having changed nothing.
#
#  "noslip" mode instead injects the volumetric split the line is actually delivering,
#  lambda_l = Vsl/(Vsl+Vsg), and lets gravity and drag decide how much liquid the segment
#  holds. The CFD holdup is then a prediction, and comparing it to SHCT's alpha_l tests
#  the drift-flux slip closure — the thing the coupling claims to calibrate.
#
#  "noslip" is now the DEFAULT. "holdup" cannot be a default for a comparison whose
#  whole purpose is to check the holdup: it imposes the answer. The tracked cases were
#  regenerated in "noslip" and their manifest numbers moved accordingly.
INLET_MODES = ("holdup", "noslip")

#  Above this mixture Reynolds number the case is written with a k-omega SST closure rather
#  than laminar. 4000 is the usual upper end of pipe transition; the sections this module
#  picks run 1e5-2.6e5, so in practice they are all turbulent and only a contrived case is
#  laminar.
RE_TURBULENT = 4000.0


def write_case(section, casedir, end_time=2.0, Ni=10, Nz=40, inlet_mode="noslip",
               n_procs=1):
    """Write a complete interFoam case. Ni (cells across the cross-section o-grid) and Nz (axial
    cells) control the mesh resolution — raise them for higher-fidelity CFD.

    inlet_mode: "holdup" imposes the SHCT liquid fraction at the inlet (default);
    "noslip" imposes the injected volumetric split lambda_l instead, so the segment's
    holdup is the CFD's own answer. See INLET_MODES above for what each can show."""
    if inlet_mode not in INLET_MODES:
        raise ValueError(f"inlet_mode must be one of {INLET_MODES}, got {inlet_mode!r}")
    n_procs = max(1, int(n_procs))
    s = section
    for d in ("0", "constant", "system"):
        os.makedirs(os.path.join(casedir, d), exist_ok=True)
    R = 0.5 * s["D"]; L = s["length_m"]
    th = s["theta_rad"]
    gy = -G * math.cos(th); gz = -G * math.sin(th)          # gravity tilted by inclination
    nu_l = s["mu_l"] / max(s["rho_l"], 1.0)
    nu_g = s["mu_g"] / max(s["rho_g"], 1e-3)
    Vm = s["Vmix"]
    #  the fraction the INLET imposes; the initial stratified level matches it so the run
    #  does not start by having to expel a slug of fluid it was never fed
    alpha = float(s["alpha_l"] if inlet_mode == "holdup"
                  else s.get("lambda_l", s["alpha_l"]))
    h = float(cx.liquid_level(np.array([alpha]))[0]); y_int = R * (2.0 * h - 1.0)

    def w(rel, txt):
        with open(os.path.join(casedir, rel), "w") as fh:
            fh.write(txt)

    # system/blockMeshDict — cylindrical o-grid pipe segment (resolution Ni x Ni x Nz)
    w("system/blockMeshDict", _blockmeshdict(R, L, Ni=Ni, Nz=Nz))

    # constant/g, transportProperties
    w("constant/g", _foam("uniformDimensionedVectorField", "g", "constant")
      + f"dimensions [0 1 -2 0 0 0 0];\nvalue ( 0 {gy:.4f} {gz:.4f} );\n")
    w("constant/transportProperties", _foam("dictionary", "transportProperties", "constant")
      + "phases (liquid gas);\n\n"
      + f"liquid\n{{\n    transportModel  Newtonian;\n"
        f"    nu      {nu_l:.6g};\n    rho     {s['rho_l']:.6g};\n}}\n\n"
      + f"gas\n{{\n    transportModel  Newtonian;\n"
        f"    nu      {nu_g:.6g};\n    rho     {s['rho_g']:.6g};\n}}\n\n"
      + f"sigma   {s['sigma']:.6g};\n")
    #  TURBULENCE. This was unconditionally `laminar`, described as the screening default.
    #  For these sections that means solving a Re = 1e5-2.6e5 flow with no turbulence model,
    #  and it is why the coupling could not adjudicate the slip closure: the drift-flux
    #  distribution parameter C0 -- where nearly all of SHCT's predicted slip lives -- IS the
    #  correlation between the turbulent velocity profile and the void profile. Measured with
    #  laminar: injecting lambda_l = 0.2537 where the 1-D closure predicts alpha_l = 0.3637
    #  returned 0.2546 at 12 diameters and 0.2561 at 50, i.e. C0 = 1.00, and neither the
    #  density ratio (1.7 -> 9.4) nor the length (12 D -> 50 D) moved it.
    #  k-omega SST is the standard RAS closure for wall-bounded pipe flow and is what
    #  interFoam tutorials use for stratified/slug cases. Laminar is still selected below the
    #  transition Reynolds number, where a model would be wrong.
    turbulent = bool(float(s.get("Re_mixture", 0.0)) > RE_TURBULENT)
    if turbulent:
        w("constant/turbulenceProperties", _foam("dictionary", "turbulenceProperties", "constant")
          + "simulationType  RAS;\n\nRAS\n{\n    RASModel        kOmegaSST;\n"
          + "    turbulence      on;\n    printCoeffs     on;\n}\n")
    else:
        w("constant/turbulenceProperties", _foam("dictionary", "turbulenceProperties", "constant")
          + "simulationType  laminar;\n")

    if turbulent:
        #  Inlet k and omega from the standard pipe estimates: I = 0.16 Re^-1/8 (fully
        #  developed duct flow), l = 0.07 D, omega = k^0.5 / (Cmu^0.25 l).
        Re = float(s["Re_mixture"])
        I_turb = 0.16 * Re ** (-0.125)
        k_in = max(1.5 * (Vm * I_turb) ** 2, 1e-8)
        l_turb = 0.07 * (2.0 * R)
        omega_in = max(k_in ** 0.5 / (0.09 ** 0.25 * l_turb), 1e-8)
        nut_in = max(k_in / omega_in, 1e-12)
        w("0/k", _foam("volScalarField", "k", "0")
          + "dimensions [0 2 -2 0 0 0 0];\n"
          + f"internalField uniform {k_in:.6g};\n\nboundaryField\n{{\n"
          + f"    inlet   {{ type fixedValue; value uniform {k_in:.6g}; }}\n"
          + f"    outlet  {{ type inletOutlet; inletValue uniform {k_in:.6g}; "
            f"value uniform {k_in:.6g}; }}\n"
          + f"    walls   {{ type kqRWallFunction; value uniform {k_in:.6g}; }}\n}}\n")
        w("0/omega", _foam("volScalarField", "omega", "0")
          + "dimensions [0 0 -1 0 0 0 0];\n"
          + f"internalField uniform {omega_in:.6g};\n\nboundaryField\n{{\n"
          + f"    inlet   {{ type fixedValue; value uniform {omega_in:.6g}; }}\n"
          + f"    outlet  {{ type inletOutlet; inletValue uniform {omega_in:.6g}; "
            f"value uniform {omega_in:.6g}; }}\n"
          + f"    walls   {{ type omegaWallFunction; value uniform {omega_in:.6g}; }}\n}}\n")
        w("0/nut", _foam("volScalarField", "nut", "0")
          + "dimensions [0 2 -1 0 0 0 0];\n"
          + f"internalField uniform {nut_in:.6g};\n\nboundaryField\n{{\n"
          + "    inlet   { type calculated; value uniform 0; }\n"
          + "    outlet  { type calculated; value uniform 0; }\n"
          + "    walls   { type nutkWallFunction; value uniform 0; }\n}\n")

    # 0/U
    w("0/U", _foam("volVectorField", "U", "0")
      + "dimensions [0 1 -1 0 0 0 0];\n"
      + f"internalField uniform (0 0 {Vm:.5g});\n\nboundaryField\n{{\n"
      + f"    inlet   {{ type fixedValue; value uniform (0 0 {Vm:.5g}); }}\n"
      + "    outlet  { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
      + "    walls   { type noSlip; }\n}\n")
    # 0/p_rgh
    w("0/p_rgh", _foam("volScalarField", "p_rgh", "0")
      + "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform 0;\n\nboundaryField\n{\n"
      + "    inlet   { type fixedFluxPressure; value uniform 0; }\n"
      + "    outlet  { type prghPressure; p uniform 0; value uniform 0; }\n"
      + "    walls   { type fixedFluxPressure; value uniform 0; }\n}\n")
    # 0/alpha.liquid
    w("0/alpha.liquid", _foam("volScalarField", "alpha.liquid", "0")
      + "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n\nboundaryField\n{\n"
      + f"    inlet   {{ type fixedValue; value uniform {alpha:.5g}; }}\n"
      + "    outlet  { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
      + "    walls   { type zeroGradient; }\n}\n")

    # system/controlDict
    w("system/controlDict", _foam("dictionary", "controlDict", "system")
      + "application     interFoam;\nstartFrom       startTime;\nstartTime       0;\n"
      + f"stopAt          endTime;\nendTime         {end_time};\ndeltaT          1e-4;\n"
      #  writeInterval was fixed at 0.1 s regardless of endTime. Any run shorter than
      #  that reached no write time at all, so the volFieldValue function object never
      #  fired: the case ran, Allrun exited 0, and ingest_results found "result file
      #  empty". couple_iterate's own default end_time=0.5 gave just five samples, of
      #  which the settled window kept three. Twenty samples per run, whatever the
      #  duration, and never longer than the run itself.
      + f"writeControl    adjustableRunTime;\nwriteInterval   {max(end_time / 20.0, 1e-6):.6g};\n"
      + "purgeWrite      0;\n"
      + "writeFormat     ascii;\nwritePrecision  6;\nwriteCompression off;\n"
      + "timeFormat      general;\ntimePrecision   6;\nrunTimeModifiable yes;\n"
      + "adjustTimeStep  yes;\nmaxCo           1;\nmaxAlphaCo      1;\nmaxDeltaT       0.01;\n\n"
      + "functions\n{\n    liquidVolAvg\n    {\n        type            volFieldValue;\n"
      + "        libs            (fieldFunctionObjects);\n        writeControl    writeTime;\n"
      + "        fields          (alpha.liquid);\n        operation       volAverage;\n"
      + "        regionType      all;\n        writeFields     false;\n    }\n}\n")
    # system/fvSchemes
    w("system/fvSchemes", _foam("dictionary", "fvSchemes", "system")
      + "ddtSchemes      { default Euler; }\n"
      + "gradSchemes     { default Gauss linear; }\n"
      + "divSchemes\n{\n    div(rhoPhi,U)        Gauss linearUpwind grad(U);\n"
      + "    div(phi,alpha)       Gauss vanLeer;\n    div(phirb,alpha)     Gauss linear;\n"
      + "    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;\n"
      + ("    div(phi,k)           Gauss linearUpwind grad(k);\n"
         "    div(phi,omega)       Gauss linearUpwind grad(omega);\n" if turbulent else "")
      + "}\n"
      + "laplacianSchemes { default Gauss linear corrected; }\n"
      + "interpolationSchemes { default linear; }\n"
      + "snGradSchemes   { default corrected; }\n"
      #  k-omega SST needs y+, so it needs a wall-distance method; without this entry
      #  interFoam aborts with "Entry 'method' not found in dictionary .../wallDist"
      #  after it has already selected the model, which reads as a turbulence failure.
      + ("wallDist        { method meshWave; }\n" if turbulent else ""))
    # system/fvSolution
    w("system/fvSolution", _foam("dictionary", "fvSolution", "system")
      + "solvers\n{\n"
      + '    "alpha.liquid.*"\n    {\n        nAlphaCorr      2;\n        nAlphaSubCycles 1;\n'
      + "        cAlpha          1;\n        MULESCorr       yes;\n        nLimiterIter    5;\n"
      + "        solver          smoothSolver;\n        smoother        symGaussSeidel;\n"
      + "        tolerance       1e-8;\n        relTol          0;\n    }\n"
      + '    "pcorr.*"        { solver PCG; preconditioner DIC; tolerance 1e-5; relTol 0; }\n'
      + "    p_rgh           { solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0.05; }\n"
      + "    p_rghFinal      { $p_rgh; relTol 0; }\n"
      + "    U               { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }\n"
      #  the trailing .* matters: PIMPLE asks for kFinal/omegaFinal on the last corrector,
      #  and a bare "(k|omega)" leaves interFoam aborting on 'omegaFinal' not found AFTER
      #  it has already taken a time step
      + ('    "(k|omega).*"   { solver smoothSolver; smoother symGaussSeidel; '
         'tolerance 1e-8; relTol 0; }\n' if turbulent else "")
      + "}\n"
      + "PIMPLE\n{\n    momentumPredictor no;\n    nOuterCorrectors 1;\n    nCorrectors 3;\n"
      + "    nNonOrthogonalCorrectors 0;\n}\n")
    # system/setFieldsDict (stratified init: liquid below the interface level)
    w("system/setFieldsDict", _foam("dictionary", "setFieldsDict", "system")
      + "defaultFieldValues ( volScalarFieldValue alpha.liquid 0 );\n\nregions\n(\n"
      + f"    boxToCell\n    {{\n        box ({-R:.5g} {-R:.5g} {-0.01:.5g}) "
        f"({R:.5g} {y_int:.5g} {L+0.01:.5g});\n"
      + "        fieldValues ( volScalarFieldValue alpha.liquid 1 );\n    }\n);\n")

    #  system/decomposeParDict — only meaningful when Allrun is asked to run in parallel,
    #  but written unconditionally so `decomposePar` works on any generated case without
    #  the user having to author one. Simple (hierarchical) decomposition along the pipe
    #  axis: the mesh is a long o-grid, so splitting it in z gives the smallest interfaces.
    w("system/decomposeParDict", _foam("dictionary", "decomposeParDict", "system")
      + f"numberOfSubdomains {n_procs};\n\nmethod          simple;\n\n"
      + f"coeffs\n{{\n    n   (1 1 {n_procs});\n}}\n")

    # Allrun
    #  Rewritten after the first REAL OpenFOAM run of these cases (v2406). Three defects
    #  in the previous four-line version, each measured:
    #
    #   1. `. ${WM_PROJECT_DIR:?}/...` — the `:?` makes an unset variable a FATAL
    #      expansion error for the `.` special builtin, so with WM_PROJECT_DIR unset the
    #      script died on line 3 with exit 2 and never reached blockMesh. `2>/dev/null ||
    #      true` does not catch it. That is exactly what run_case() produces when the
    #      caller has the solvers on PATH but has not sourced the OpenFOAM bashrc — the
    #      condition openfoam_available() reports as available.
    #   2. The three stages ran unconditionally and the script ended with `echo done`, so
    #      a case in which blockMesh, setFields AND interFoam had all failed still exited
    #      0 (measured on an empty PATH). run_case(check=True) reads that as success and
    #      records "ran": true for a case where nothing ran.
    #   3. `cp -r 0 0.orig` is not idempotent: setFields overwrites 0/alpha.liquid, so on
    #      a second run the already-initialised field was copied INTO the existing backup
    #      as 0.orig/0 and the pristine state was lost.
    allrun = ('#!/bin/sh\ncd "${0%/*}" || exit 1\n'
              'if [ -n "$WM_PROJECT_DIR" ] && [ -f "$WM_PROJECT_DIR/bin/tools/RunFunctions" ]; then\n'
              '    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"\n'
              'fi\n'
              'if [ -d 0.orig ]; then rm -rf 0 && cp -r 0.orig 0; else cp -r 0 0.orig; fi\n'
              'run() { "$1" > "log.$1" 2>&1 || { echo "FAILED: $1 (see log.$1)" >&2; exit 1; }; }\n'
              'run blockMesh\n'
              'run setFields\n'
              + (f'NP={n_procs}\n'
                 'rm -rf processor*\n'
                 'run decomposePar\n'
                 '#  reconstructPar is required: ingest_results reads postProcessing, which a\n'
                 '#  parallel run writes from the master anyway, but the time directories it\n'
                 '#  checks for as evidence live under processor*/ until reconstruction.\n'
                 'mpirun -np $NP interFoam -parallel > log.interFoam 2>&1 \\\n'
                 '    || { echo "FAILED: interFoam (see log.interFoam)" >&2; exit 1; }\n'
                 'run reconstructPar\n'
                 if n_procs > 1 else 'run interFoam\n')
              + 'echo done\n')
    w("Allrun", allrun)
    os.chmod(os.path.join(casedir, "Allrun"), 0o755)

    # README + section metadata
    w("README.txt",
      f"OpenFOAM interFoam case auto-generated by the SHCT->OpenFOAM coupling.\n"
      f"Section: {s['name']}  (x = {s['x_km']:.2f} km)\n"
      f"Why this section needs 3-D CFD: {s['reason']}\n\n"
      f"Geometry: pipe segment D = {s['D']*1000:.0f} mm, L = {L:.2f} m, "
      f"inclination = {math.degrees(th):.1f} deg (gravity tilted accordingly).\n\n"
      f"Boundary conditions FROM the SHCT 1-D solution at this section:\n"
      f"  inlet mixture velocity Vm = {Vm:.3f} m/s\n"
      f"  inlet liquid fraction = {alpha:.3f}  (stratified init to h/D = {h:.3f})\n"
      f"  inlet_mode = {inlet_mode} — "
      + ("the SHCT holdup is IMPOSED here, so the run cannot independently confirm it\n"
         if inlet_mode == "holdup" else
         "the injected volumetric split is imposed; the holdup is the CFD's own answer\n")
      + f"  SHCT alpha_l = {s['alpha_l']:.3f}, injected lambda_l = {s.get('lambda_l', float('nan')):.3f}\n"
      f"  pressure (context) = {s['p_bar']:.1f} bar, temperature = {s['T_C']:.1f} C\n"
      f"  liquid rho/mu = {s['rho_l']:.0f}/{s['mu_l']:.2e}, gas rho/mu = {s['rho_g']:.1f}/{s['mu_g']:.2e}\n\n"
      + (f"\nFLOW REVERSED: the SHCT mixture velocity here is {s['Vmix_shct_signed']:+.3f} m/s.\n"
         f"  The segment is therefore written in the flow's own direction: inclination\n"
         f"  {math.degrees(th):+.1f} deg is the SHCT inclination "
         f"{math.degrees(s['theta_rad_shct']):+.1f} deg\n"
         f"  with the axis reversed, and the inlet patch is the downstream end of the line.\n"
         if s.get("flow_reversed") else "")
      + (f"\nTURBULENCE: k-omega SST (RAS), Re = {s['Re_mixture']:,.0f}.\n"
         f"  This case was previously written simulationType laminar at this Reynolds\n"
         f"  number, and that is why it could not check the slip closure: the drift-flux\n"
         f"  distribution parameter C0 -- where nearly all of the 1-D model's predicted\n"
         f"  slip lives -- IS the turbulent velocity/void profile correlation. Measured\n"
         f"  laminar, injecting lambda_l = 0.254 where the closure predicts alpha_l =\n"
         f"  0.364 returned 0.255 at 12 diameters and 0.256 at 50, and neither the\n"
         f"  density ratio nor the length moved it.\n"
         if turbulent else
         f"\nTURBULENCE: laminar, Re = {s.get('Re_mixture', 0.0):,.0f} "
         f"(below the {RE_TURBULENT:,.0f} transition).\n")
      + (f"\nWARNING - THE TWO PHASES ARE NOT DISTINCT: rho_gas = {s['rho_g']:.1f} against\n"
         f"  rho_liquid = {s['rho_l']:.0f} kg/m3, a density ratio of {s['rho_ratio']:.1f}. The mixture is\n"
         f"  single-phase at {s['p_bar']:.1f} bar / {s['T_C']:.1f} C, so the PVT table has no vapour to\n"
         f"  report and its gas column returns the dense/liquid root. interFoam will solve\n"
         f"  this case happily, but it is not the gas-liquid system the 1-D model assumed\n"
         f"  when it predicted a liquid fraction of {s['alpha_l']:.3f} here.\n"
         if not s.get("phases_distinct", True) else "")
      + (f"\nWARNING - BOUNDARY CONDITION FROM A CLIPPED STATE: {s['velocity_clip_fraction']:.0%} of the\n"
         f"  SHCT realisations at this cell sit exactly on the 1-D solver's mixture-velocity\n"
         f"  clip, so the inlet velocity above is a numerical guard rather than a converged\n"
         f"  physical state. A CFD-vs-SHCT holdup difference for this section measures the\n"
         f"  clip, not the closure. Treat it as an exploratory run.\n"
         if s.get("bc_from_clipped_state") else "")
      + "\nRun on an OpenFOAM machine:   ./Allrun     (needs blockMesh, setFields, interFoam)\n"
      "Then ingest with shct_openfoam.ingest_results('<this dir>').\n")
    w("section.json", json.dumps(dict(s, inlet_mode=inlet_mode,
                                      inlet_alpha=alpha), indent=2))
    return casedir


# ---------------------------------------------------------------------------
#  Run a case if OpenFOAM is available
# ---------------------------------------------------------------------------
def openfoam_available():
    #  setFields is checked too: Allrun runs it between blockMesh and interFoam, so a
    #  partial installation missing it would be reported available and then fail mid-run.
    return all(shutil.which(t) is not None
               for t in ("blockMesh", "setFields", "interFoam"))


def _log_tail(casedir, name, n=12):
    """Last n non-blank lines of a stage log, for the failure reason."""
    try:
        with open(os.path.join(casedir, f"log.{name}"), errors="replace") as fh:
            lines = [ln.rstrip() for ln in fh if ln.strip()]
        return "\n".join(lines[-n:])
    except OSError:
        return ""


def run_case(casedir, timeout=3600):
    """Run ./Allrun if OpenFOAM is installed; return a status dict.

    A zero exit status is NOT taken as proof the case ran. Allrun is a shell script a
    user may have edited, and the version this module shipped before the coupling was
    first exercised against real OpenFOAM returned 0 even when every stage had failed.
    The evidence that interFoam ran is its own output: a time directory beyond 0 and a
    log ending in "End". Both are checked here, and the tail of the log that actually
    failed is returned as the reason so the caller does not have to go looking."""
    if not openfoam_available():
        return {"ran": False, "reason": "OpenFOAM not found on PATH "
                "(blockMesh/setFields/interFoam). Case written; run ./Allrun on an "
                "OpenFOAM machine."}
    try:
        proc = subprocess.run(["./Allrun"], cwd=casedir, timeout=timeout,
                              capture_output=True, text=True, errors="replace")
    except Exception as exc:
        return {"ran": False, "reason": f"Allrun could not be executed: {exc}"}
    if proc.returncode != 0:
        stage = "interFoam"
        for name in ("blockMesh", "setFields", "interFoam"):
            if "FAILED: " + name in (proc.stderr or ""):
                stage = name
                break
        return {"ran": False, "stage": stage, "returncode": proc.returncode,
                "reason": f"{stage} failed:\n{_log_tail(casedir, stage)}"}
    #  exit 0 is necessary but not sufficient — verify interFoam's own evidence
    times = [t for t in os.listdir(casedir)
             if _as_time(t) is not None and _as_time(t) > 0.0]
    log = _log_tail(casedir, "interFoam", n=4)
    if not times or not log.endswith("End"):
        return {"ran": False, "stage": "interFoam", "returncode": 0,
                "reason": "Allrun exited 0 but interFoam produced no time directory "
                          f"beyond 0 and/or did not reach End:\n{log}"}
    return {"ran": True, "casedir": casedir, "latest_time": max(map(_as_time, times))}


def _as_time(name):
    """OpenFOAM time-directory name -> float, or None if it is not one."""
    try:
        return float(name)
    except (TypeError, ValueError):
        return None


#  Fraction of the record discarded as start-up before the holdup is time-averaged.
#  The case is initialised stratified at the SHCT holdup and then has to establish its
#  own interface shape and slip; in the first real interFoam run of these cases the
#  volume-averaged liquid fraction of section 1 moved 0.3598 -> 0.3657 over the first
#  half-second and only then began to oscillate about its mean.
SETTLE_FRAC = 0.5


def ingest_results(casedir, settle_frac=SETTLE_FRAC):
    """Read the CFD result back: the volume-averaged liquid fraction from the function
    object (if the case has run).

    Two things the first real interFoam run showed, which the synthetic tests could not:

      * The reported holdup OSCILLATES. Section 1 passed through 0.3598, 0.3657 and
        0.3639 at 0.1 s, 1.2 s and 1.6 s — a +-0.3 % swing about its mean. Reporting the
        LAST sample, as this function used to, therefore returned wherever in that swing
        the run happened to stop, and `abs_diff` against SHCT inherited the phase of the
        oscillation rather than measuring the closure. The comparison value is now the
        time MEAN over the settled window, with the spread returned alongside so a
        difference smaller than the oscillation cannot be read as agreement.
      * Time directories sort LEXICOGRAPHICALLY. `sorted(os.listdir(base))` puts "10"
        before "2" and "0.5" before "0.9" is fine but "0.5" after "0.15"; a restarted
        run (OpenFOAM writes one postProcessing directory per start time) was therefore
        concatenated out of order and `vals[-1]` was not the latest time at all. The
        samples are now sorted by their numeric time.
    """
    base = os.path.join(casedir, "postProcessing", "liquidVolAvg")
    if not os.path.isdir(base):
        return {"available": False, "reason": "no postProcessing output (case not run yet)"}
    vals = []
    for tdir in sorted(os.listdir(base), key=lambda d: (_as_time(d) is None, _as_time(d) or 0.0)):
        d = os.path.join(base, tdir)
        if not os.path.isdir(d):
            continue
        #  Re-running a case in place does NOT overwrite volFieldValue.dat. OpenFOAM keeps
        #  the existing file and writes the new run to volFieldValue_0.dat (then _1, ...),
        #  so reading only the canonical name returns the STALE run -- or, when the previous
        #  attempt died before its first write, an empty header and "result file empty" for
        #  a case that in fact ran to End. Measured: a k-omega SST case re-run after a failed
        #  attempt left a 4-line volFieldValue.dat beside a complete 24-line
        #  volFieldValue_0.dat. Read them all, oldest file first, so later writes win.
        files = sorted((os.path.join(d, fn) for fn in os.listdir(d)
                        if fn.startswith("volFieldValue") and fn.endswith(".dat")),
                       key=lambda f: os.path.getmtime(f))
        for f in files:
            with open(f, errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    try:
                        vals.append((float(parts[0]), float(parts[-1])))
                    except (ValueError, IndexError):
                        pass
    if not vals:
        return {"available": False, "reason": "result file empty"}
    #  a restart re-reports times already present; keep the last value written for each
    by_t: dict[float, float] = {}
    for t_s, a_s in vals:
        by_t[t_s] = a_s
    ts = sorted(by_t)
    t_end = ts[-1]
    t0 = float(settle_frac) * t_end
    window = [t for t in ts if t >= t0] or [t_end]
    a_win = np.array([by_t[t] for t in window], float)
    return {"available": True,
            "time": t_end,
            "cfd_mean_liquid_fraction": float(a_win.mean()),   # settled-window time mean
            "cfd_final_liquid_fraction": float(by_t[t_end]),   # last instantaneous sample
            "cfd_min_liquid_fraction": float(a_win.min()),
            "cfd_max_liquid_fraction": float(a_win.max()),
            "cfd_swing": float(a_win.max() - a_win.min()),
            "window_start_s": float(window[0]),
            "n_samples": int(a_win.size)}


# ---------------------------------------------------------------------------
#  Orchestrator
# ---------------------------------------------------------------------------
def couple(sv, outdir, max_sections=3, run=False, end_time=2.0, Ni=10, Nz=40,
           inlet_mode="noslip", seg_len_factor=12.0, n_procs=1):
    """Identify critical sections, write an OpenFOAM case for each, optionally run
    them (if OpenFOAM is installed), and write a manifest. Returns the manifest dict.
    When run=True and OpenFOAM is on PATH the interFoam result is ingested and the
    SHCT-vs-CFD volume-averaged liquid-holdup difference recorded per section."""
    root = os.path.join(outdir, "openfoam_cases")
    os.makedirs(root, exist_ok=True)
    #  seg_len_factor was a parameter of identify_critical_sections that no caller could
    #  reach: couple() did not accept it and the CLI did not expose it, so every case ever
    #  generated was 12 diameters long. That is the binding limit on what the comparison
    #  can show. Measured on interFoam v2406, section 1 at x = 23.1 km, inlet_mode=noslip
    #  (injecting lambda_l = 0.2537 where SHCT predicts alpha_l = 0.3637, a slip of
    #  +0.110): the CFD returns 0.2546 -- a slip of +0.0009. Raising the density ratio
    #  from 1.7 to 9.4 does not change that (0.2502, slip -0.0036), so it is not buoyancy;
    #  over 12 diameters the injected split simply convects through before wall friction
    #  can redistribute it. The knob is now reachable so a segment long enough to develop
    #  can be asked for.
    secs = identify_critical_sections(sv, max_sections=max_sections,
                                      seg_len_factor=seg_len_factor)
    #  DROP THE SECTIONS THIS RUN DID NOT PRODUCE. A section directory is named after the
    #  station it was cut at, so when the physics moves the critical sections move and the
    #  new run writes new names beside the old ones -- nothing ever removed them. The
    #  as-operated folder had accumulated fourteen section directories from six historical
    #  runs while the manifest named three, so eleven runnable interFoam cases were sitting
    #  there carrying boundary conditions from superseded solutions and looking current.
    #  Only directories matching the name this function itself generates are removed.
    _keep = {s["name"] for s in secs}
    for _d in sorted(os.listdir(root)):
        _p = os.path.join(root, _d)
        if os.path.isdir(_p) and _d.startswith("section_") and _d not in _keep:
            shutil.rmtree(_p, ignore_errors=True)
    manifest: dict = {"n_sections": len(secs), "openfoam_available": openfoam_available(),
                      "end_time": end_time, "mesh": {"Ni": Ni, "Nz": Nz},
                      "inlet_mode": inlet_mode, "seg_len_factor": seg_len_factor,
                      "n_procs": int(max(1, n_procs)),
                      #  say what a difference in this manifest can and cannot show
                      "comparison_meaning": (
                          "inlet_mode=holdup imposes alpha_l at the inlet, so rel_diff_pct is a "
                          "consistency check on the case writer, not a validation of the holdup "
                          "closure (measured: imposing 0.20/0.36/0.60 returns 0.2015/0.3637/0.5988)"
                          if inlet_mode == "holdup" else
                          "inlet_mode=noslip injects lambda_l and lets the segment find its own "
                          "holdup, so rel_diff_pct tests the drift-flux slip closure"),
                      "sections": []}
    for s in secs:
        casedir = os.path.join(root, s["name"])
        write_case(s, casedir, end_time=end_time, Ni=Ni, Nz=Nz, inlet_mode=inlet_mode,
                   n_procs=n_procs)
        entry = {"name": s["name"], "x_km": s["x_km"], "reason": s["reason"],
                 "Vmix": s["Vmix"], "alpha_l": s["alpha_l"], "Phi_SH": s["Phi_SH"],
                 "Vmix_shct_signed": s["Vmix_shct_signed"],
                 "flow_reversed": s["flow_reversed"],
                 "velocity_clip_fraction": s["velocity_clip_fraction"],
                 "bc_from_clipped_state": s["bc_from_clipped_state"],
                 "rho_ratio": s["rho_ratio"], "phases_distinct": s["phases_distinct"],
                 "Re_mixture": s["Re_mixture"],
                 "lambda_l": s["lambda_l"],
                 "casedir": os.path.relpath(casedir, outdir)}
        if run:
            status = run_case(casedir)
            entry["run"] = status
            if status.get("ran"):
                cfd = ingest_results(casedir)
                entry["cfd"] = cfd
                if cfd.get("available"):
                    entry["shct_alpha_l"] = float(s["alpha_l"])
                    entry["cfd_alpha_l"] = float(cfd["cfd_mean_liquid_fraction"])
                    entry["abs_diff"] = abs(entry["shct_alpha_l"] - entry["cfd_alpha_l"])
                    entry["rel_diff_pct"] = entry["abs_diff"] / max(entry["shct_alpha_l"], 1e-6) * 100.0
                    #  A difference smaller than the run's own oscillation is not
                    #  agreement, it is the resolution floor of the comparison. Say so
                    #  in the manifest rather than leaving abs_diff to be read as a
                    #  validation error.
                    entry["cfd_swing"] = float(cfd.get("cfd_swing", 0.0))
                    entry["resolved"] = bool(entry["abs_diff"] > entry["cfd_swing"])
        manifest["sections"].append(entry)
    with open(os.path.join(root, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    # a top-level index/README
    with open(os.path.join(root, "README.txt"), "w") as fh:
        fh.write("SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.\n"
                 f"OpenFOAM detected: {manifest['openfoam_available']}.\n"
                 "Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the\n"
                 "SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call\n"
                 "shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.\n\n")
        for e in manifest["sections"]:
            fh.write(f"  - {e['name']}: x={e['x_km']:.2f} km — {e['reason']}\n")
            if not e.get("phases_distinct", True):
                fh.write(f"      WARNING: density ratio {e['rho_ratio']:.1f} — the mixture is "
                         "single-phase here, so this is not a gas-liquid VOF case.\n")
            if e.get("bc_from_clipped_state"):
                fh.write(f"      WARNING: {e['velocity_clip_fraction']:.0%} of realisations here "
                         "sit on the 1-D velocity clip; its BC is a numerical guard.\n")
            elif e.get("flow_reversed"):
                fh.write(f"      note: flow is reversed ({e['Vmix_shct_signed']:+.3f} m/s); "
                         "the segment is written in the flow's own direction.\n")
    return manifest


# ---------------------------------------------------------------------------
#  TWO-WAY (closed-loop) coupling: CFD result -> SHCT closure -> re-run (item 6)
# ---------------------------------------------------------------------------
def couple_iterate(case, outdir, max_sections=3, max_iters=4, tol=0.02, gain=0.8,
                   run=False, synthetic_cfd=None, end_time=0.5, Ni=10, Nz=40,
                   inlet_mode="noslip"):
    """Closed-loop SHCT<->OpenFOAM coupling. Each iteration:
      1. run SHCT for `case`;
      2. locate the critical sections and write their OpenFOAM cases;
      3. obtain the CFD liquid holdup per section — from a real interFoam run when
         OpenFOAM is installed (run=True), else from `synthetic_cfd(section)` if given
         (for testing / when CFD is unavailable);
      4. compute the SHCT-vs-CFD holdup mismatch and, while it exceeds `tol`, apply a
         damped feedback to the SHCT drift-flux distribution parameter C0 via
         numerics.drift_C0_factor (higher C0 -> faster gas -> lower gas holdup ->
         higher liquid holdup) and re-run.
    Returns {history, calibrated_case}. This is a REDUCED two-way coupling — CFD informs a
    global SHCT closure knob — not a full domain-decomposition co-simulation; it converges
    the whole-line model toward the CFD-resolved holdup of the critical sections.

    inlet_mode defaults to "noslip" HERE, unlike couple(). With "holdup" the inlet imposes
    the very quantity the loop is trying to calibrate, so against a real interFoam run the
    mismatch is ~0 at iteration 0 (measured: 0.3637 imposed returns 0.3637) and the loop
    exits having changed nothing. Only "noslip" leaves the holdup for the CFD to predict,
    which is what makes the feedback meaningful. The synthetic_cfd path is unaffected."""
    import solver  # lazy (avoid import cycle)
    cur = copy.deepcopy(case)
    history = []
    root = os.path.join(outdir, "openfoam_coupling_iter")
    os.makedirs(root, exist_ok=True)
    for it in range(max_iters):
        sv = solver.TransientSHCT(cur); sv.run(verbose=False)
        secs = identify_critical_sections(sv, max_sections=max_sections)
        pairs: list[tuple[float, float]] = []
        failures: list[dict] = []
        itdir = os.path.join(root, f"iter{it}")
        for s in secs:
            cd = os.path.join(itdir, s["name"])
            write_case(s, cd, end_time=end_time, Ni=Ni, Nz=Nz, inlet_mode=inlet_mode)
            shct_h = s["alpha_l"]
            if synthetic_cfd is not None:
                cfd_h = float(synthetic_cfd(s))
            elif run and openfoam_available():
                #  The run status was previously discarded, so a section whose interFoam
                #  had failed simply dropped out of `pairs` and the iteration carried on
                #  averaging the ones that survived -- silently coupling to a subset.
                st = run_case(cd)
                if not st.get("ran"):
                    failures.append({"section": s["name"], "reason": st.get("reason", "")})
                    cfd_h = None
                else:
                    got = ingest_results(cd)
                    cfd_h = got.get("cfd_mean_liquid_fraction")
                    if cfd_h is None:
                        failures.append({"section": s["name"],
                                         "reason": got.get("reason", "no CFD result")})
            else:
                cfd_h = None
            if cfd_h is not None:
                pairs.append((shct_h, float(cfd_h)))
        if not pairs:
            note = ("OpenFOAM not available — CFD cases generated only "
                    "(pass synthetic_cfd or run on an OpenFOAM machine)")
            if failures:
                note = "every CFD section failed to produce a holdup"
            history.append({"iter": it, "mismatch": None,
                            "roughness_m": cur.pipeline.roughness_m,
                            "cfd_failures": failures, "note": note})
            break
        mism = float(np.mean([abs(a - b) for a, b in pairs]))
        signed = float(np.mean([b - a for a, b in pairs]))     # cfd - shct
        mean_h = float(np.mean([a for a, _ in pairs]))
        history.append({"iter": it, "mismatch": mism, "signed": signed,
                        "n_sections_coupled": len(pairs), "cfd_failures": failures,
                        "drift_C0_factor": float(getattr(cur.numerics, "drift_C0_factor", 1.0)),
                        "roughness_m": cur.pipeline.roughness_m})
        if mism < tol:
            break
        # FEEDBACK: tune the drift-flux distribution parameter C0 — the strong, physically-grounded
        # holdup knob (higher C0 -> higher liquid holdup). signed>0 (CFD holds more liquid) -> raise C0.
        factor = float(np.clip(1.0 + gain * signed / max(mean_h, 1e-3), 0.5, 2.0))
        cur = copy.deepcopy(cur)
        cur.numerics.drift_C0_factor = float(np.clip(
            getattr(cur.numerics, "drift_C0_factor", 1.0) * factor, 0.5, 2.0))
    with open(os.path.join(root, "coupling_history.json"), "w") as fh:
        json.dump(history, fh, indent=2)
    return {"history": history, "calibrated_case": cur}
