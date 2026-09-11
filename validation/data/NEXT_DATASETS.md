# Datasets that would close the remaining validation gaps

> **Status, read this first.** Items 1 and 4 are RETRIEVED and in
> `validation/data/`; item 4 is now RUN and SCORED. Items 2 and 3 remain open
> and are blocked by paywalls, not by effort. The per-item notes below are the
> original scouting notes, kept because they record what was tried.

Each entry is open-access or transcribable, directly relevant, and blocked only by
automated-access restrictions — not by availability. Each takes minutes in a browser.

## 1. RETRIEVED — Oil-dominated / high-water-cut flowloop plug outcomes

**Pham, T.-K., Cameirão, A., Melchuna, A., Herri, J.-M., Glénat, P. (2020),**
*"Relative Pressure Drop Model for Hydrate Formation and Transportability in Flowlines
in High Water Cut Systems"*, **Energies 13(3):686**, doi:10.3390/en13030686.
**CC-BY, gold open access** (confirmed via Unpaywall).

Why this one: high water cut is this case's duty (70 %), and the measured quantity is
**relative pressure drop during hydrate formation** — the coupled slug/hydrate response
this solver predicts and currently has no measurement for. Werner-Bolley, by contrast,
is a gas-condensate well (4 MMSCFD gas, 100 BPD condensate, 10 BPD water) and would not
validate a medium-crude line at 70 % water cut whatever its data showed.

What to extract: water cut, liquid loading, mixture velocity, temperature, pressure,
oil and gas used, pipe diameter and length, and the relative-pressure-drop series
against hydrate volume fraction or time.

Mirrors that refused automated access: mdpi.com (403), hal.science / hal-emse (denied),
doaj.org. **The attachment host `res.mdpi.com` did not**, and Table 1 was read from the
publisher PDF there — 12 tests with plug/no-plug outcomes, now in
`flowloop_plugging_pham2020.json` and scored in `score_pham2020_result.txt`.
What the paper does NOT tabulate is the relative-pressure-drop SERIES against time; the
reduced Kv values are what is published, so the coupled dP response is still unmeasured.

## 2. ExxonMobil Friendswood flowloop, Conroe crude — OIL-DOMINATED

Reported conditions: 70 % liquid loading, 75 mol% methane / 25 mol% ethane, 1400 rpm,
**37 vol% water cut**, cooled to 4 °C; a companion series in a 4-inch loop covers
1–2.5 m/s mixture velocity and 50–90 vol% liquid loading at 100 % water cut. Used as the
validation set for CSMHyK-OLGA. Instrumented for pressure, temperature, density and
differential pressure. Mostly SPE / Elsevier; check institutional access.

## 3. Englezos et al. (1987) — the SPECIFIC AREA, to finish the kg0 reconciliation

`kinetics_englezos1987_comparison.json` reconciles this project's kg0 with the measured
intrinsic constant to within an order of magnitude, and the residual factor is entirely
the AREA BASIS: K* multiplies hydrate particle surface in a stirred slurry, kg0
multiplies a_wall = 4/D (14.6 m2/m3 at 10.75 in). Closing it needs the particle size
distribution or specific surface area those experiments were reduced against — Englezos
used a population-balance treatment, so the number exists in the primary paper
(Chem. Eng. Sci. 42(11):2647-2658) but is not carried by any citing source found here.

## 4. RETRIEVED AND SCORED — see field_olga_csmhyk_roberts.json

Roberts, T., *Oil and Gas Field Application of Hydrate Kinetics Modeling*, MEng thesis,
Memorial University of Newfoundland (open repository). Couples **OLGA** with **CSMHyK**
and applies it to two real field cases: one confirmed hydrate blockage (spool pressure
built from ~5,300 to 8,000 kPa; model predicted the location correctly) and one confirmed
NO-blockage on restart. Also publishes a full field crude assay, C1 48.83 mol% through
C30+ 5.10, with lab PVT: bubble point 29,150 kPag and stock-tank density 883 kg/m3.

## 5. What is STILL open after all of the above

- **A point-by-point profile from a real line.** Roberts gives outcomes — a pressure build
  and a blockage location — not a dP / arrival-T / holdup series this solver could be
  scored against, term by term. That specific artefact is still missing.
- **This solver has not been run against ANY of the above.** Every dataset here is
  recorded as evidence and as a route. Nothing in this repository is scored against the
  field cases, the OLGA coupling, or the Pham flowloop set.
- **Phi_crit and C_phi have no measured counterpart at all.** kg0 is now reconciled to
  order of magnitude and growth_exp_n is bounded by published data; these two are not.
- **The Englezos specific area**, which would turn the kg0 reconciliation into a number.
