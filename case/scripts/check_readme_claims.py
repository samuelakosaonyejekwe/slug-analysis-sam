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
        print(f"{lbl:34s} {'—':>10s} {'':>10s}   claim text not found (reworded?)")
        continue
    said = float(m.group(1))
    got = km[sc].get(key)
    if got is None:
        print(f"{lbl:34s} {said:10.4g} {'—':>10s}   metric absent"); continue
    got = float(got) * (100.0 if key == "P_plug" else 1.0)
    ok = abs(got - said) <= tol * max(abs(said), 1e-9) + (0.5 if key == "P_plug" else 0.0)
    print(f"{lbl:34s} {said:10.4g} {got:10.4g}   {'ok' if ok else 'STALE — README must change'}")
    bad += 0 if ok else 1
#  the verdict sentence itself
verdict_says_no_plug = "sub-critical and does not plug" in txt
actually_plugs = float(km["steady"].get("P_plug", 0.0)) > 0.0
print(f"\n  README verdict 'does not plug' : {verdict_says_no_plug}")
print(f"  solver P_plug (as-operated)    : {km['steady'].get('P_plug')}")
if verdict_says_no_plug and actually_plugs:
    print("  *** VERDICT REVERSED — the README's central conclusion is now wrong ***"); bad += 1
print(f"\n{'all README numbers hold' if bad == 0 else f'{bad} README claim(s) need updating'}")
raise SystemExit(1 if bad else 0)
