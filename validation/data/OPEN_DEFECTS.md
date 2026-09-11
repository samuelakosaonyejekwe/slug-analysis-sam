# Open defects and unfinished work

> Closed entries are struck through and kept for one revision, so the record of what
> was wrong and how it was settled survives the fix.

Everything flagged during the audit that is **not** closed. Each entry says what is wrong,
how it was found, and what would close it. Anything closed is removed from this file, not
ticked — the record of what was fixed lives in the artefact it was fixed in.

## ~~1. Φ_SH's equilibrium identity~~ — CLOSED

Tested POINTWISE (same cell, same realisation, same instant, 480 h, on the 21 of 280 cells
depositing and not locked): **median ratio of predicted to actual deposit 1.050**, p10
0.824, p90 1.228. The identity holds; the balance and identity forms agree to machine
precision, as the algebra requires.

The reported factor-of-two was an artefact of comparing `max_Phi_SH` — a running maximum
the code's own comment says lands on the startup spike — against `peak_deposit_mm`, a
median-over-ensemble maximum over cells. A local identity cannot be tested with non-local
aggregates, and two further wrong explanations were built on that artefact before it was
caught. `case/scripts/test_phish_identity.py`.

## 2. `consol_restriction = 0.18` — bracketed, not measured

Has a measured counterpart now (DOE occupancy: plugged at f = 0.117, transported at
f = 0.041) and the shipped value sits **above** that bracket. What is still unmeasured is
the quantity the parameter literally encodes: the thickness at which a deposit stops being
**removable**. The DOE work reports sloughing as events in time, never as a critical
thickness. **To close:** a flow-loop measurement of deposit thickness at the onset of
irreversible retention.

## 3. dP / arrival-T profile from a real line — OPEN

Roberts gives field OUTCOMES, no series. das Neves' shipped pressure channels are
piezoelectric and AC-coupled — fluctuation, not level. **Route identified and partly
walked:** the paper's Table 2 lists per-point `Pressure gradient [Pa/m]` in the metadata
`.pkl` files on the Unicamp share (doi:10.25824/redu/ISMWP4). The share serves files but
refuses `PROPFIND`, so paths must be guessed; `/Dataset_25_6kHz/Point_8` is confirmed to
exist, the metadata filename is not yet known. `case/scripts/fetch_flowloop_pressure_gradient.py`
is written and waiting on that filename. **To close:** read one folder listing.

## ~~4. `sustained_supercritical_km` uses a third threshold~~ — CLOSED

It counts cells where Φ_SH > **1.0** — not Φ_crit nominal (1.08) and not Φ_crit operating
(1.456). Φ_SH > 1 means δ_eq > δ_ref, a different statement from "consolidates", while the
NAME implies the second. Left unchanged because shipped outputs and the README quote it;
`sustained_above_phicrit_km` now ships alongside, against the operating Φ_crit.

## 5. Wall-exponent constraint has no discriminating power

The real-line negative (Roberts case 2) is reproduced by every exponent in 1.0–1.5 at the
real 11 h duration, because nucleation onset (~4 h at 7.4 °C subcooling) leaves almost
nothing deposited by then. So the one real line available does **not** constrain the wall
exponent. **To close:** a real line whose outcome is sensitive to the exponent — i.e. one
that ran long enough, or cold enough, to deposit measurably.

## 6. Anti-agglomerant is not modelled at all

Why Pham's 12 real flowloop outcomes cannot score this solver: AA dose is the deciding
variable in 11 of 12, and there is no AA term. A scope limit, not a missing dataset.

## 7. Real-line magnitude and timing disagree with the reference code

Against OLGA + CSMHyK on Roberts case 2: riser-base hydrate **2.66× the reference bound**,
peak location mid-flowline where the reference puts it at the manifold and riser base, and
growth accelerating at 4 h against the reference's ~17 h. The bit (no plug) is right.
**To close:** unknown — this is a real disagreement with a licensed code, not a bug found.
