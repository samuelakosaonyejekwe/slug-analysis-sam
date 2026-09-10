#!/usr/bin/env python3
"""Score THIS project's structural claim against real flowloop plug outcomes.

The claim Phi_SH encodes is that whether a line plugs is set by a dimensionless
COMPETITION, not by how much hydrate exists. Pham et al. (2020) measured 12 high-water-
cut flowloop tests and recorded, for each, the hydrate volume fraction at the end and
whether it plugged. That is exactly the discriminating test: if plugging were an
inventory threshold, hydrate volume fraction would separate the outcomes.

This does NOT score Phi_SH itself. Phi_SH needs a wall subcooling and a slug frequency,
and Pham tabulates neither; and their mechanism is agglomeration in a water-continuous
suspension, whereas Phi_SH is wall deposition against shear removal. What is scored is
the SHAPE of the argument, on real data, against the obvious alternative.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DS = os.path.join(os.path.dirname(HERE), "..", "validation", "data",
                  "flowloop_plugging_pham2020.json")


def auc(scores, labels):
    """Probability a randomly chosen plugged test outranks a randomly chosen intact one."""
    pos = [s for s, y in zip(scores, labels) if y]
    neg = [s for s, y in zip(scores, labels) if not y]
    if not pos or not neg:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def main():
    with open(os.path.normpath(DS)) as fh:
        ds = json.load(fh)
    P = ds["points"]
    hyd = np.array([r[10] for r in P], float)      # hydrate volume fraction at end
    kv = np.array([r[7] for r in P], float)        # agglomerate structure factor
    rate = np.array([r[6] for r in P], float)      # crystallisation rate
    plug = np.array([bool(r[9]) for r in P])

    print("=" * 72)
    print(" SCORING THE STRUCTURAL CLAIM AGAINST REAL PLUG OUTCOMES")
    print(" Pham et al. (2020) Energies 13(3):686 -- Archimede flowloop, n =", len(P))
    print("=" * 72)
    print(f"  plugged: {int(plug.sum())} of {len(P)}")
    print()
    print("  candidate predictor        AUC     separation (plug vs intact)")
    print("  " + "-" * 62)
    for name, v in (("hydrate volume fraction", hyd),
                    ("agglomerate factor Kv", kv),
                    ("crystallisation rate", rate)):
        a = auc(v, plug)
        print("  %-24s  %5.3f   plug %.3f +/- %.3f | intact %.3f +/- %.3f"
              % (name, a, v[plug].mean(), v[plug].std(),
                 v[~plug].mean(), v[~plug].std()))
    print()
    a_inv, a_str = auc(hyd, plug), auc(kv, plug)
    print("  INVENTORY (hydrate volume fraction) AUC = %.3f" % a_inv)
    print("  STRUCTURE (agglomerate factor Kv)   AUC = %.3f" % a_str)
    print()
    if a_str > a_inv:
        print("  -> Structure outranks inventory on real data. An inventory threshold does")
        print("     NOT separate these outcomes; a dimensionless structural factor does.")
        print("     That is the shape Phi_SH asserts, corroborated on a different mechanism.")
    else:
        print("  -> Inventory outranks structure here. The structural argument is NOT")
        print("     supported by this dataset and the claim should be weakened.")
    print()
    print("  AUC 0.5 = no discrimination, 1.0 = perfect. Inventory BELOW 0.5 means larger")
    print("  hydrate volume is associated with NOT plugging -- the opposite of a threshold.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
