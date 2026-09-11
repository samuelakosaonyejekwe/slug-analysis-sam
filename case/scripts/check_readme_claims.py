"""Every number the README asserts, checked against the run that is actually on disk.

The physics corrections move T_eq up ~1 C, which moves subcooling, deposit and possibly
the plug verdict with it. A README that still says "sub-critical and does not plug" when
the corrected solver says otherwise would be worse than the gap it fixed. This re-reads
each claim from key_metrics.json and reports any that no longer hold.
"""
import json
import os
import re

R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if not os.path.exists(os.path.join(R, "README.md")):
    R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
km = {s: json.load(open(f"{R}/case/outputs_{s}/key_metrics.json"))
      for s in ("steady", "shutin", "mitigated")}
txt = open(f"{R}/README.md").read()

#  (label, scenario, metric key, regex that captures the number the README states, tol)
CLAIMS = [
 ("max subcooling (as-operated)", "steady", "max_subcooling_C",
  r"drives the fluid ([\d.]+) °C into the hydrate", 0.05),
 ("peak wall deposit (as-operated)", "steady", "peak_deposit_mm",
  r"a ([\d.]+) mm peak wall\s*\n?deposit", 0.05),
 ("MEG dose", "steady", "MEG_wt_pct",
  r"inhibition at ([\d.]+) wt% MEG", 0.02),
 ("under-inhibited length", "steady", "under_inhibited_km",
  r"over a ([\d.]+) km under-inhibited length", 0.05),
 ("peak Phi_SH (as-operated)", "steady", "max_Phi_SH",
  r"peak Φ_SH ([\d.]+)", 0.05),
 ("plug probability (as-operated)", "steady", "P_plug",
  r"(\d+) % plug probability", 0.01),
 ("no-touch time (mitigated)", "mitigated", "cooldown_to_hydrate_h",
  r"buys a ([\d.]+) h no-touch time", 0.05),
]
bad = 0
print(f"{'claim':34s} {'README':>10s} {'actual':>10s}   status")
for lbl, sc, key, rx, tol in CLAIMS:
    m = re.search(rx, txt)
    if not m:
        #  A CLAIM THAT VANISHED IS A FAILURE, NOT A PASS. This printed a note and moved
        #  on, so rewording the README -- or deleting a headline number outright -- left
        #  the checker reporting "all README numbers hold". Three deliverables (peak wall
        #  deposit, MEG dose, under-inhibited length) were dropped in a rewrite and this
        #  is what let it through: the check that exists to stop the README drifting from
        #  the outputs was silent about the README no longer making the claim at all.
        print(f"{lbl:34s} {'—':>10s} {'':>10s}   CLAIM MISSING — the README no longer "
              f"states this; restore it or delete the check")
        bad += 1
        continue
    said = float(m.group(1))
    got = km[sc].get(key)
    if got is None:
        print(f"{lbl:34s} {said:10.4g} {'—':>10s}   metric absent"); continue
    got = float(got) * (100.0 if key == "P_plug" else 1.0)
    ok = abs(got - said) <= tol * max(abs(said), 1e-9) + (0.5 if key == "P_plug" else 0.0)
    print(f"{lbl:34s} {said:10.4g} {got:10.4g}   {'ok' if ok else 'STALE — README must change'}")
    bad += 0 if ok else 1
#  ---- the growth_exp_n table, which is now the HEADLINE result -------------------
#  It states the central engineering finding -- whether this duty plugs is decided by the
#  subcooling exponent -- and until now nothing checked it against the sweep it is read
#  from. It also used to MIX SOURCES: the n = 1.00 row was taken from key_metrics.json
#  (the full 48 h headline run, 0.531) while the other four came from the reduced 24 h /
#  6-realisation sweep (baseline 0.576), so the row anchoring the table was not comparable
#  with the rows below it. Every row must now come from the sweep, and match it.
_sens = os.path.join(R, "case", "outputs_steady", "sensitivity_phiSH.csv")
if os.path.exists(_sens):
    import csv as _csv
    with open(_sens) as _fh:
        _rows = list(_csv.DictReader(_fh))
    _by_n = {}
    for _r in _rows:
        if _r["label"] == "baseline" or _r["label"].startswith("n_"):
            _by_n[round(float(_r["growth_exp_n"]), 2)] = _r
    #  pull "| 1.75 | 3.18 | **0.50** |"-shaped rows out of the README
    _pat = re.compile(r"^\|[^|]*?(\d\.\d\d)\b[^|]*\|\s*\**([0-9.]+)\**\s*\|"
                      r"\s*\**([0-9.]+)\**\s*\|", re.M)
    _seen = 0
    for _m in _pat.finditer(txt):
        _n = round(float(_m.group(1)), 2)
        if _n not in _by_n:
            continue
        _seen += 1
        _phi_said, _pp_said = float(_m.group(2)), float(_m.group(3))
        _phi_got = float(_by_n[_n]["max_Phi_SH"])
        _pp_got = float(_by_n[_n]["P_plug"])
        _ok_phi = abs(_phi_got - _phi_said) <= 0.01 * max(_phi_said, 1e-9) + 0.006
        _ok_pp = abs(_pp_got - _pp_said) <= 0.01
        print(f"growth_exp_n={_n:<5} Phi_SH {_phi_said:>6.2f} vs {_phi_got:>6.3f}"
              f"   P_plug {_pp_said:>4.2f} vs {_pp_got:>5.3f}   "
              f"{'ok' if (_ok_phi and _ok_pp) else 'STALE — README must change'}")
        bad += 0 if (_ok_phi and _ok_pp) else 1
    if _seen < 4:
        print(f"  *** the growth_exp_n table matched only {_seen} sweep rows — "
              f"the headline table is not being checked ***")
        bad += 1
else:
    print("  *** sensitivity_phiSH.csv absent — the headline n-table cannot be checked ***")
    bad += 1

#  the verdict sentence itself
verdict_says_no_plug = "sub-critical and does not plug" in txt
actually_plugs = float(km["steady"].get("P_plug", 0.0)) > 0.0
print(f"\n  README verdict 'does not plug' : {verdict_says_no_plug}")
print(f"  solver P_plug (as-operated)    : {km['steady'].get('P_plug')}")
if verdict_says_no_plug and actually_plugs:
    print("  *** VERDICT REVERSED — the README's central conclusion is now wrong ***"); bad += 1
print(f"\n{'all README numbers hold' if bad == 0 else f'{bad} README claim(s) need updating'}")
raise SystemExit(1 if bad else 0)
