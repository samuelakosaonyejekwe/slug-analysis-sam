#!/usr/bin/env python3
"""Test delta_eq = Phi_SH * delta_ref POINTWISE, which is how it is stated.

The identity is LOCAL and per-realisation. It was previously "checked" by comparing
max_Phi_SH against peak_deposit_mm -- a running maximum that the code's own comment says
lands on the startup spike, against a median-over-ensemble maximum over cells. Those are
different places at different times, and the factor-of-two mismatch that produced was an
artefact of the comparison, not a property of the model. Twice.

So: same cell, same realisation, same instant, at steady state.

    growth   d(delta)/dt = f_wall * Rg_wall * D/4
    erosion  d(delta)/dt = k_ero * f_slug * delta          (where not locked)
    balance  delta_eq    = f_wall*Rg_wall*D / (4*k_ero*f_slug)
    and      Phi_SH      = C_phi * Rg_wall / f_slug
    so       delta_eq    = Phi_SH * f_wall*D/(4*C_phi*k_ero) = Phi_SH * delta_ref_local

Every term on the right is now exported per cell, so the check is arithmetic.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import solver  # noqa: E402


def main():
    import run_case_study10 as cs
    hours = float(sys.argv[1]) if len(sys.argv) > 1 else 480.0
    c = cs.build_case("Phi_SH identity test", "asoperated", hours)
    c.numerics.n_ensemble = 4
    sv = solver.TransientSHCT(c); sv.run(verbose=False)
    r, k = sv.results, c.kinetics
    D = c.pipeline.diameter_m

    fw = np.asarray(r["f_wall_final"], float)
    rg = np.asarray(r["Rg_wall_final"], float)
    fs = np.asarray(r["fslug_final"], float)
    dl = np.asarray(r["delta"], float)
    ph = np.asarray(r["PhiSH_final_field"], float)
    lk = np.asarray(r["locked_final"], bool)

    #  only cells that are actually depositing and NOT locked -- the identity assumes an
    #  uncapped growth/erosion balance, and a locked cell has erosion switched off.
    live = (fw > 1e-9) & (rg > 1e-12) & (fs > 1e-12) & (~lk) & (dl > 1e-6)
    print("=" * 74)
    print(" delta_eq = Phi_SH * delta_ref, TESTED POINTWISE AT %.0f h" % hours)
    print("=" * 74)
    print("  cells depositing, unlocked: %d of %d" % (int(live.sum()), dl.size))
    if not live.any():
        print("  nothing to test.")
        return 1

    dref = fw * D / (4.0 * k.C_phi * k.k_ero)          # local delta_ref
    d_bal = fw * rg * D / (4.0 * k.k_ero * fs)          # from the balance directly
    d_idn = ph * dref                                   # via the identity
    a, b, cc = dl[live], d_bal[live], d_idn[live]

    def stat(name, pred):
        ratio = pred / np.maximum(a, 1e-12)
        print("  %-28s median %7.3f   p10 %7.3f   p90 %7.3f" %
              (name, np.median(ratio), np.percentile(ratio, 10), np.percentile(ratio, 90)))
    print("  ratios of PREDICTED to ACTUAL deposit:")
    stat("balance  f*Rg*D/(4*kero*fs)", b)
    stat("identity Phi_SH * delta_ref", cc)
    same = np.allclose(b, cc, rtol=1e-9)
    print("  balance and identity agree algebraically: %s" % same)
    print()
    med = float(np.median(cc / np.maximum(a, 1e-12)))
    if 0.8 <= med <= 1.25:
        print("  IDENTITY HOLDS pointwise (median ratio %.3f). The earlier factor-of-two" % med)
        print("  mismatch was an artefact of comparing non-co-located maxima.")
    else:
        print("  IDENTITY DOES NOT HOLD pointwise: median ratio %.3f." % med)
        print("  The cells are still approaching equilibrium, or a term in the balance is")
        print("  missing from the identity as documented. Reported, not explained away.")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
