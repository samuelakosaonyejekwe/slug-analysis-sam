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
# =============================================================================
from __future__ import annotations
import os, math, json, shutil, subprocess, copy
import numpy as np

try:
    from shct_correlations import gas_density
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
        return (f"    {name}\n    {{\n        type {typ};\n        faces\n        (\n{f}\n        );\n    }}\n")
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
    med = lambda A: np.nanmedian(A, 1)
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
    picked = []
    order = np.argsort(score)[::-1]
    min_sep = max(int(0.04 * len(x)), 3)
    for i in order:
        if all(abs(i - j) > min_sep for j in picked):
            picked.append(int(i))
        if len(picked) >= max_sections:
            break
    picked = sorted(picked)
    rho_l = float(sv.rho_l)
    secs = []
    for rank, i in enumerate(picked):
        p_bar = float(med(r["p"])[i]); T_C = float(med(r["T"])[i])
        rho_g = float(gas_density(p_bar, T_C, sv.case.fluids)) if gas_density else 30.0
        Di = float(D[i]); A = math.pi * Di ** 2 / 4.0
        secs.append(dict(
            name=f"section_{rank+1}_x{ x[i]/1000:.1f}km".replace(".", "p"),
            index=i, x_km=float(x[i] / 1000.0), length_m=float(seg_len_factor * Di),
            D=Di, theta_rad=float(sv.theta[i]),
            Vmix=float(max(med(r["j"])[i], 0.05)), alpha_l=float(np.clip(alpha[i], 1e-3, 0.999)),
            p_bar=p_bar, T_C=T_C, Phi_SH=float(phish[i]), subcooling_C=float(sub[i]),
            deposit_mm=float(delta[i] * 1000.0), regime=int(round(regime[i])),
            rho_l=rho_l, rho_g=rho_g,
            mu_l=float(sv.case.fluids.mu_liquid), mu_g=float(sv.case.fluids.mu_gas),
            sigma=float(sv.case.fluids.sigma),
            reason=_why(phish[i], theta[i], intermittent[i], sub[i], delta[i])))
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
def write_case(section, casedir, end_time=2.0, Ni=10, Nz=40):
    """Write a complete interFoam case. Ni (cells across the cross-section o-grid) and Nz (axial
    cells) control the mesh resolution — raise them for higher-fidelity CFD."""
    s = section
    for d in ("0", "constant", "system"):
        os.makedirs(os.path.join(casedir, d), exist_ok=True)
    R = 0.5 * s["D"]; L = s["length_m"]
    th = s["theta_rad"]
    gy = -G * math.cos(th); gz = -G * math.sin(th)          # gravity tilted by inclination
    nu_l = s["mu_l"] / max(s["rho_l"], 1.0)
    nu_g = s["mu_g"] / max(s["rho_g"], 1e-3)
    Vm = s["Vmix"]; alpha = s["alpha_l"]
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
      + f"liquid\n{{\n    transportModel  Newtonian;\n    nu      {nu_l:.6g};\n    rho     {s['rho_l']:.6g};\n}}\n\n"
      + f"gas\n{{\n    transportModel  Newtonian;\n    nu      {nu_g:.6g};\n    rho     {s['rho_g']:.6g};\n}}\n\n"
      + f"sigma   {s['sigma']:.6g};\n")
    #  required by interFoam (turbulence model selector); laminar is the screening default
    w("constant/turbulenceProperties", _foam("dictionary", "turbulenceProperties", "constant")
      + "simulationType  laminar;\n")

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
      + "writeControl    adjustableRunTime;\nwriteInterval   0.1;\npurgeWrite      0;\n"
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
      + "    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;\n}\n"
      + "laplacianSchemes { default Gauss linear corrected; }\n"
      + "interpolationSchemes { default linear; }\n"
      + "snGradSchemes   { default corrected; }\n")
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
      + "    U               { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }\n}\n"
      + "PIMPLE\n{\n    momentumPredictor no;\n    nOuterCorrectors 1;\n    nCorrectors 3;\n"
      + "    nNonOrthogonalCorrectors 0;\n}\n")
    # system/setFieldsDict (stratified init: liquid below the interface level)
    w("system/setFieldsDict", _foam("dictionary", "setFieldsDict", "system")
      + "defaultFieldValues ( volScalarFieldValue alpha.liquid 0 );\n\nregions\n(\n"
      + f"    boxToCell\n    {{\n        box ({-R:.5g} {-R:.5g} {-0.01:.5g}) ({R:.5g} {y_int:.5g} {L+0.01:.5g});\n"
      + "        fieldValues ( volScalarFieldValue alpha.liquid 1 );\n    }\n);\n")

    # Allrun
    allrun = ("#!/bin/sh\ncd \"${0%/*}\" || exit 1\n"
              ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions 2>/dev/null || true\n"
              "blockMesh > log.blockMesh 2>&1\n"
              "cp -r 0 0.orig 2>/dev/null || true\n"
              "setFields > log.setFields 2>&1\n"
              "interFoam > log.interFoam 2>&1\n"
              "echo done\n")
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
      f"  inlet liquid fraction (holdup) alpha_l = {alpha:.3f}  (stratified init to h/D = {h:.3f})\n"
      f"  pressure (context) = {s['p_bar']:.1f} bar, temperature = {s['T_C']:.1f} C\n"
      f"  liquid rho/mu = {s['rho_l']:.0f}/{s['mu_l']:.2e}, gas rho/mu = {s['rho_g']:.1f}/{s['mu_g']:.2e}\n\n"
      f"Run on an OpenFOAM machine:   ./Allrun     (needs blockMesh, setFields, interFoam)\n"
      f"Then ingest with shct_openfoam.ingest_results('<this dir>').\n")
    w("section.json", json.dumps(s, indent=2))
    return casedir


# ---------------------------------------------------------------------------
#  Run a case if OpenFOAM is available
# ---------------------------------------------------------------------------
def openfoam_available():
    return shutil.which("blockMesh") is not None and shutil.which("interFoam") is not None


def run_case(casedir, timeout=3600):
    """Run ./Allrun if OpenFOAM is installed; return a status dict."""
    if not openfoam_available():
        return {"ran": False, "reason": "OpenFOAM not found on PATH (blockMesh/interFoam). "
                "Case written; run ./Allrun on an OpenFOAM machine."}
    try:
        subprocess.run(["./Allrun"], cwd=casedir, check=True, timeout=timeout,
                       capture_output=True)
        return {"ran": True, "casedir": casedir}
    except Exception as exc:                             # pragma: no cover
        return {"ran": False, "reason": f"run failed: {exc}"}


def ingest_results(casedir):
    """Read the CFD result back: the volume-averaged liquid fraction from the
    function object (if the case has run). Returns a dict; empty if no results."""
    base = os.path.join(casedir, "postProcessing", "liquidVolAvg")
    if not os.path.isdir(base):
        return {"available": False, "reason": "no postProcessing output (case not run yet)"}
    # find the latest time dir with a volFieldValue.dat
    vals = []
    for t in sorted(os.listdir(base)):
        f = os.path.join(base, t, "volFieldValue.dat")
        if os.path.isfile(f):
            for line in open(f):
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    try:
                        vals.append((float(parts[0]), float(parts[-1])))
                    except Exception:
                        pass
    if not vals:
        return {"available": False, "reason": "result file empty"}
    t_last, alpha_cfd = vals[-1]
    return {"available": True, "time": t_last, "cfd_mean_liquid_fraction": alpha_cfd}


# ---------------------------------------------------------------------------
#  Orchestrator
# ---------------------------------------------------------------------------
def couple(sv, outdir, max_sections=3, run=False, end_time=2.0, Ni=10, Nz=40):
    """Identify critical sections, write an OpenFOAM case for each, optionally run
    them (if OpenFOAM is installed), and write a manifest. Returns the manifest dict.
    When run=True and OpenFOAM is on PATH the interFoam result is ingested and the
    SHCT-vs-CFD volume-averaged liquid-holdup difference recorded per section."""
    root = os.path.join(outdir, "openfoam_cases")
    os.makedirs(root, exist_ok=True)
    secs = identify_critical_sections(sv, max_sections=max_sections)
    manifest = {"n_sections": len(secs), "openfoam_available": openfoam_available(),
                "end_time": end_time, "mesh": {"Ni": Ni, "Nz": Nz}, "sections": []}
    for s in secs:
        casedir = os.path.join(root, s["name"])
        write_case(s, casedir, end_time=end_time, Ni=Ni, Nz=Nz)
        entry = {"name": s["name"], "x_km": s["x_km"], "reason": s["reason"],
                 "Vmix": s["Vmix"], "alpha_l": s["alpha_l"], "Phi_SH": s["Phi_SH"],
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
        manifest["sections"].append(entry)
    manifest["sections_detail"] = secs
    with open(os.path.join(root, "manifest.json"), "w") as fh:
        json.dump({k: v for k, v in manifest.items() if k != "sections_detail"}, fh, indent=2)
    # a top-level index/README
    manifest.pop("sections_detail", None)
    with open(os.path.join(root, "README.txt"), "w") as fh:
        fh.write("SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.\n"
                 f"OpenFOAM detected: {manifest['openfoam_available']}.\n"
                 "Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the\n"
                 "SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call\n"
                 "shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.\n\n")
        for e in manifest["sections"]:
            fh.write(f"  - {e['name']}: x={e['x_km']:.2f} km — {e['reason']}\n")
    return manifest


# ---------------------------------------------------------------------------
#  TWO-WAY (closed-loop) coupling: CFD result -> SHCT closure -> re-run (item 6)
# ---------------------------------------------------------------------------
def couple_iterate(case, outdir, max_sections=3, max_iters=4, tol=0.02, gain=0.8,
                   run=False, synthetic_cfd=None, end_time=0.5, Ni=10, Nz=40):
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
    the whole-line model toward the CFD-resolved holdup of the critical sections."""
    import solver                                          # lazy (avoid import cycle)
    cur = copy.deepcopy(case)
    history = []
    root = os.path.join(outdir, "openfoam_coupling_iter")
    os.makedirs(root, exist_ok=True)
    for it in range(max_iters):
        sv = solver.TransientSHCT(cur); sv.run(verbose=False)
        secs = identify_critical_sections(sv, max_sections=max_sections)
        pairs = []
        itdir = os.path.join(root, f"iter{it}")
        for s in secs:
            cd = os.path.join(itdir, s["name"]); write_case(s, cd, end_time=end_time, Ni=Ni, Nz=Nz)
            shct_h = s["alpha_l"]
            if synthetic_cfd is not None:
                cfd_h = float(synthetic_cfd(s))
            elif run and openfoam_available():
                run_case(cd); cfd_h = ingest_results(cd).get("cfd_mean_liquid_fraction")
            else:
                cfd_h = None
            if cfd_h is not None:
                pairs.append((shct_h, cfd_h))
        if not pairs:
            history.append({"iter": it, "mismatch": None,
                            "roughness_m": cur.pipeline.roughness_m,
                            "note": "OpenFOAM not available — CFD cases generated only "
                                    "(pass synthetic_cfd or run on an OpenFOAM machine)"})
            break
        mism = float(np.mean([abs(a - b) for a, b in pairs]))
        signed = float(np.mean([b - a for a, b in pairs]))     # cfd - shct
        mean_h = float(np.mean([a for a, _ in pairs]))
        history.append({"iter": it, "mismatch": mism, "signed": signed,
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
