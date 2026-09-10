# SHCT Solver — Transient Coupled-PDE Prediction of Slugging & Hydrates in Subsea Pipelines

_Author: **Akosa Samuel Onyejekwe**_

`solver.py` is an engineering-grade **transient, coupled partial-differential-equation**
solver that predicts hydrodynamic slugging and gas-hydrate formation in subsea
multiphase pipelines on arbitrary terrain/inclination, and emits the full output
set (fields, time-series, space-time maps, slug statistics, probabilistic risk and
engineering design deliverables). It is the computational realisation of the SHCT
invention and of its output catalogue. (Earlier revisions pointed here at `1/work.docx`
and `2/output.docx`; neither directory is part of this repository — the equation set and
the full output catalogue are in `report.pdf`, Sections 4 and 11.)

---

## 1. What it solves — the coupled PDE system

Integrated **in time** by a conservative finite-volume method with adaptive
CFL-controlled stepping and a stochastic Monte-Carlo ensemble:

| Eq. | Physics | Form |
|----|---------|------|
| **H** | Liquid-holdup transport (drift-flux kinematic wave) | `∂(αl·A)/∂t + ∂(αl·vl·A)/∂x = 0` (volume-conservative, upwind) |
| **G** | **Transient momentum & pressure — 3 engines** | `implicit` (default): mixture momentum + implicit pressure + drift-flux slip. `twofluid`: **full two-fluid** — two independent phase momenta + interfacial drag + shared pressure (slip emerges from momentum balance). `quasisteady`: legacy + auto-fallback. All → tridiagonal **implicit pressure** each step |
| **E** | Energy transport | `∂T/∂t + vm·∂T/∂x = −U·Pw·(T−Tsea)/(ρm·cp·A) + L·ṁh/(ρm·cp)` |
| **P** | Hydrate phase-field (advected reaction–diffusion) | `∂φ/∂t + vl·∂φ/∂x = Dφ·∂²φ/∂x² + Rgrow + Rnuc + ξ` (stochastic nucleation) |
| **D** | Wall-deposit / consolidation with slug scouring | two-way coupled to H,G,E |
| **I** | Thermodynamic-inhibitor (MEG) transport + local hydrate-curve suppression | `Teq_eff = Teq(p) − ΔT_NielsenBucklin(W)` |
| **Gm** | Gas-mass continuity (conservative, mass-consistent flux) | `∂(ρg·αg·A)/∂t + ∂(ρg·Vsg·A)/∂x = −ṁ_gas→hydrate` |
| **C** | Slug–Hydrate Coupling Number | `Φ_SH = C·kg,wall·a_wall·ΔTsub,wall^n / f_slug`, `a_wall = (4/D)·αl·f_water` |

**OLGA-class capabilities added:** multi-layer **wall heat transfer** (steel/insulation/
coating → effective U) and **cooldown / no-touch-time**; **inhibitor design** (required
MEG wt% & rate + an *under-inhibited-length* map; inject with `--meg 30`); hydrate-slurry
**transportability** (Camargo–Palermo relative viscosity); an automated **V&V suite**
(correlation benchmarks + liquid/gas/hydrate mass conservation + headline-output grid
convergence) and **calibration**.

**Conservation & coupling (this revision):** the energy equation carries the **exothermic
hydrate latent heat** (bulk formation self-limits → the bulk slurry stays transportable);
the hydrate phase-field carries its **`D_phi` diffusion** term; a conservative **gas-mass
continuity** equation is solved alongside the liquid so **liquid, gas and hydrate mass all
balance to ~0 %**; the **wall deposit is mass-coupled to hydrate growth** (driven by the
sustained cold-wall subcooling, with the water/gas it consumes removed from the liquid/gas
inventories). Plugging is therefore **heat-transfer-limited and probabilistic** (a genuine
P10/P50/P90 time-to-plug spread), not a single deterministic time.

**Two-way coupling:** hydrate fraction & deposit reduce the effective bore and raise
friction/viscosity (→ H,G,E); slug hydrodynamics set interfacial area, shear and
subcooling that drive nucleation/growth/deposition (→ P,D).

**Published closures:** Bendiksen drift-flux, Taitel–Dukler regime, Gregory–Scott /
Zabaras slug frequency, Haaland friction, natural-gas hydrate equilibrium,
CSMHyK-type growth, Camargo–Palermo viscosity, **Nielsen–Bucklin** inhibitor
suppression/sizing (valid to high MEG wt%, where the classical Hammerschmidt linear
form breaks down — the API is kept under the historical `hammerschmidt_meg` name),
API-RP-14E erosional limit.

---

## 2. Outputs generated (written to `solver_outputs/`)

**Tables (CSV)**
- `fields_profile.csv` — full along-line P50 profile: P, T, Teq, subcooling, holdup,
  velocity, regime, interfacial area, slug frequency, hydrate fraction, deposit, Φ_SH.
- `timeseries_monitor.csv` — transient monitor-point history of all key variables + BC.
- `probabilistic_summary.csv` — P10/P50/P90 of subcooling, Φ_SH, peak hydrate fraction, time-to-plug.
- `engineering_deliverables.csv` — slug-catcher surge volume, MEG concentration & rate,
  peak vs erosional velocity, total ΔP, plug probability, time-to-plug, hot-spot, Φ_SH, mass-error.
  The surge volume takes whichever of two bases governs — a hydrodynamic slug period
  `(q_l/f_slug)·surge_factor`, or the liquid inventory `A·Σ(α_l·dx)` held in the contiguous
  ascent that reaches the outlet (P90 across the ensemble), the latter only where that ascent
  is steeper than 10°. Both are reported alongside the governing value, with the basis named
  (`V_surge_hydrodynamic_m3`, `V_riser_liquid_m3`, `riser_incline_deg`, `V_surge_basis`).
- `summary.json` — machine-readable scalar metrics.

**Charts (PNG)**
1. `01_profiles.png` — elevation / holdup / P–T–Teq / subcooling profiles.
2. `02_holdup_spacetime.png` — transient liquid-holdup field α_l(x,t).
3. `03_PT_envelope.png` — P–T trajectory vs hydrate envelope.
4. `04_PhiSH_map.png` — Φ_SH(x,t) coupling-criticality map with the Φ_SH = 1 contour.
5. `05_scenario_timeseries.png` — BC, subcooling and Φ_SH response to the operating scenario.
6. `06_deposit.png` — transient wall-deposit growth.
7. `07_probabilistic.png` — time-to-plug CDF + Φ_SH ensemble band along the line.
8. `08_diagnostics.png` — liquid/gas mass-balance error, clip activity, sub-grid slug
   length, and the hydrate/fallback/step counters.

**Space-time / multi-time set (PNG)** — written by the `shct_spacetime` extension
module (`shct_spacetime.spacetime_outputs(sv, eng, outdir)`), which the case-study
driver calls after `make_charts`. It renders the run in the figure family the
transient-multiphase literature uses:

9.  `14_holdup_multitime.png` — holdup along the route at successive times (early / late).
10. `15_slug_growth_propagation.png` — resolved slug units at three times, one front tracked.
11. `16_slug_train_waterfall.png` — distance–time slug tracking, celerity scan, moveout, stack.
12. `17_hydrate_distribution.png` — in-pipe volume fractions + phase rates into the host.
13. `18_shutin_profile_deposit.png` — late-time P/T/T_eq/water profile + deposit growth.
14. `19_spacetime_fields.png` — the true space-time solution: α_l, p, u_g, u_l, ΔT_sub, deposit.
15. `20_holdup_durations.png` — holdup along the line after successive elapsed durations.
16. `21_riser_depth_time.png` — riser depth–time waterfall with slug boundaries.
17. `22_cloud_maps.png` — bore phase distribution + temperature strips at successive times.
18. `23_dts_thermal_waterfall.png` — distributed-temperature waterfall T(x,t) with the
    monitored pressure overlaid and the operating stages marked.
19. `24_temperature_gradient.png` — temperature-gradient waterfall ∂T/∂x(x,t), which
    localises a travelling front where the temperature map itself looks smooth.
20. `25_das_flow_noise.png` — flow-noise waterfall |∂α_l/∂t|(x,t), with the intermittent
    reach and the riser base marked.
21. `26_parameter_panels.png` — P, T, holdup and mixture velocity along the route at the
    same successive times.
22. `27_wellposedness_map.png` — the two-fluid (Kelvin–Helmholtz) well-posedness boundary
    over the superficial-velocity plane, and the margin along the route.

Figures 15, 16 and 21 show individual slug units, which are **sub-grid** at dx ≈ 460 m.
They are a kinematic reconstruction built from the solver's own slug statistics
(`f_slug`, `L_u = V_t/f_slug`, `α_ls`), with the body/film split solved so the
unit-averaged holdup returns the cell-average `α_l` exactly — mass-consistent by
construction, and labelled as such on every figure. See `shct_spacetime.py`.

**Snapshot history in `results`** — on the `n_snapshots` cadence the solver records
`snap_t`, `snap_holdup`, `snap_P`, `snap_T`, `snap_phi`, `snap_PhiSH`, and (from v3.2)
`snap_delta` (wall-deposit thickness), `snap_Tsub` (subcooling), `snap_j` (mixture
velocity), `snap_vl` / `snap_vg` (phase velocities from the momentum/pressure solve),
`snap_regime` and `snap_fslug`. These are what the space-time figures read.

**Console report** — hydrodynamics, thermal/hydrate, coupling, risk and engineering
deliverables, plus solver diagnostics (step count, hydro fallbacks, mass-conservation error).

---

## 3. Operating scenarios (time-dependent boundary conditions)

Selectable with `--scenario`: `steady` (default), `rampup`, `turndown`, `shutin`.
These exercise the genuinely transient core (e.g. a shut-in cools the line and
collapses slug renewal → hydrate plugging), demonstrating real operational cases.
**All scenarios conserve liquid AND gas mass** (steady/turndown/shut-in: ~0.00 %;
ramp-up: <1 % at production resolution), including full shut-in, via an
outflow-only / isolated-line outlet and a **local-first** (nearest-neighbour spill)
conservative bound-enforcement step — plus the hydrate water/gas sinks, so the
balances close *with* the hydrate consumption, not merely in = out.

## 3a. Verification & Calibration (V&V)

The honest path from "tuned defaults" to "trusted on a real asset" is built in:

```bash
python3 solver.py --verify                 # VERIFICATION: automated checks that the
                                           # closures reproduce published reference values
                                           # (hydrate Teq, Gregory-Scott, Haaland, Hammerschmidt,
                                           # drift-flux limits) and that the core conserves mass.

python3 solver.py --calibrate targets.json # VALIDATION: fit the free constants (heat transfer,
                                           # hydrate growth/nucleation, deposition) to YOUR
                                           # measured data, then writes calibrated_case.json.
```

`targets.json` holds whatever you measured, e.g.
`{"arrival_T_C": 8.0, "dP_total_bar": 30.0, "max_subcooling_C": 9.0, "time_to_plug_P50_h": 40.0}`.

**Check the target-reachability table the run prints before believing a residual.** The four
free constants are thermal and kinetic (`U_wall`, `kg0`, `wall_capture_eff`, `nuc_beta_C`), and
they cannot move every metric. Measured on the bundled case, perturbing each by ±10 % around
the optimum moves `max_subcooling_C` by 28 %, `arrival_T_C` by 1.0 % and `dP_total_bar` by
0.02 % — so two of the three targets in the example line above are **unreachable with the
default free set**, and an 80 % residual against `dP_total_bar` is structural, not a bad fit.
(Arrival temperature is insensitive because a 32 km line has already equilibrated to the
seabed; pressure drop is insensitive because none of the four constants touches friction.)
The run now prints that table and warns by name. Either drop such a target or widen `free` to
a parameter that controls it.
The optimiser (Nelder–Mead) adjusts the constants to match, so the same solver is
adapted to a specific fluid / field / flow-loop dataset before quantitative use. This
is what allows the same solver to be applied beyond the single bundled case — after
calibration to data from the system in question, not before.

---

## 4. Usage

```bash
python3 solver.py                      # bundled real case (20 km tie-back, steady)
python3 solver.py --scenario shutin    # shut-in cooldown / hydrate-risk transient
python3 solver.py --engine twofluid    # full two-fluid (two independent phase momenta)
python3 solver.py --meg 30             # inject 30 wt% MEG inhibitor
python3 solver.py --config case.json   # any user case
python3 solver.py --dump-config c.json # write an editable default case template
python3 solver.py --no-plots           # tables + console only
python3 solver.py --outdir results/    # choose output directory
```

**Engines** (`--engine` or `numerics.engine`) — six, not three; the CLI accepts all of them:
`implicit` (default — drift-flux + implicit pressure, fast & robust), `twofluid` (full
two-fluid: two independent phase momenta + interfacial drag; slip emerges from the momentum
balance), `twofluid_mass` (damped within-step momentum/gas-volume coupling),
`twofluid_mass_newton` (monolithic volume-mass Newton — a consistency-priority engine that
trades dP fidelity for gas-holdup consistency, and warns at run time that it does),
`twofluid_full_newton` (block-tridiagonal Newton on the primitives simultaneously), and
`quasisteady` (legacy + auto-fallback). All conserve mass to ~0 %; `implicit` and `twofluid`
are the two the `--verify` suite exercises directly.

**Figure self-checks.** Every figure is saved through a helper that first runs three
checks and prints what it finds, so a bad figure is reported at the moment it is written
rather than found by eye later:

* `find_text_overlaps` — any text drawn over other text. This is what caught the 38
  collisions in the slide build and the annotation sitting on the x-axis label in
  `25_das_flow_noise.png`.
* `find_empty_axes` — a panel that draws nothing at all.
* `find_degenerate_axes` — **a panel whose data is too sparse to be the profile it
  claims**: a line declaring many points whose finite data covers under 5 % of the axis,
  or fewer than three points. The first two checks could not see this one. On the
  case-study crude the EOS splits at 1 of 40 stations, so `compo_pvt.png` drew a
  K-value-*vs-distance* panel as eight markers stacked against the right-hand edge of a
  blank axis — under a full eight-entry legend for curves that did not exist — and two
  gas-property panels as two markers over 2.5 % of a 32 km axis. There were artists, so
  nothing flagged them. Those panels now plot against COMPONENT (heavy → light), which is
  what a single split state can actually be read on, or state the value in words; and the
  check runs over every figure this project draws, so the next one is caught.

**Slide renderings.** `SHCT_FIG_FONTSCALE` enlarges every text element so a figure survives
being shrunk onto a slide, and `SHCT_FIG_SIZESCALE` shrinks the canvas. They were applied
independently, so the deck build's 1.8x text on a 1.0x canvas got no layout adaptation at all
and produced 38 text collisions, and its 1.8x-on-0.30x variant produced 171. Both scales now
drive the same crowding ratio — tick thinning, short labels and canvas headroom key on
FONTSCALE/SIZESCALE — and the ratio is capped at 1.8 with a warning when the ask exceeds it.
The primary slide set is now collision-free. The 0.30/0.45/0.70 sub-scale sets are not, and no
font setting will fix them: with **no** font enlargement at all those same figures still
collide 8-9 times at 0.30-0.45, because a four-panel figure drawn at a third of its design
size has nowhere to put a legible label. They are gitignored derivatives; the fix, if they are
ever needed clean, is fewer panels per figure.

**OpenFOAM.** `openfoam_available()` checked PATH only, and OpenFOAM is not on PATH until its
`etc/bashrc` is sourced — so on a machine with v2406 installed it returned False, the coupling
wrote cases without running them, and `test_real_interfoam_run_end_to_end` skipped with "not on
PATH". Sourced, that test passes. An installed-but-unsourced OpenFOAM is now detected and
reported, so the skip says which of the two it is.

`key_metrics.gas_holdup_consistency` measures what the default engine gives up for that
speed: the median relative gap between the CONSERVED gas mass `Mg` and the gas holdup the
algebraic drift-flux closure returns. It is a gap, not a conservation error — the gas mass
balance itself closes to ~1e-14 on every scenario — but it is not small, and it is not
uniform: on the case study it reads 0.055 as-operated, 0.091 mitigated and **0.59 on the
shut-in**, where the flow has stopped and an algebraic slip closure has least to work with.
The shut-in is the scenario the study's headline hazard comes from, so read its holdup field
with that number beside it. `twofluid_mass_newton` drives the same gap below 0.01 and says at
run time what it costs in pressure-drop fidelity.

**Engine-consistent, mass-coupled plug prediction.** The wall deposit is now driven by the
hydrate growth at the **cold-wall subcooling** (Teq − Tsea), a robust, sustained quantity
that is nearly engine-invariant, and the water/gas it consumes is removed from the conserved
liquid/gas inventories. So all three engines give consistent plugging at ~0% liquid/gas mass
error with 0 fallbacks. Because the **bulk** hydrate φ self-limits via its own latent heat,
it stays low (transportable slurry) while the **wall** deposit builds to a plug — the same
field picture as before, but now a single **mass-consistent partition** of one hydrate source.
Plugging is **heat-transfer-limited and probabilistic**: the stochastic nucleation induction
time plus the default parameter spread give a genuine **P10/P50/P90 time-to-plug** band rather
than one deterministic time.

A case is fully described by JSON groups `pipeline`, `fluids`, `operating`,
`kinetics`, `numerics`, `scenario` — including an arbitrary `elevation_m[]` terrain
profile. Run `--dump-config` to see every tunable field.

---

## 5. Numerical quality / engineering-software practices (VERIFIED)

- Conservative finite-volume holdup **and gas-mass** transport → liquid & gas mass conserved
  to **~0.00 %** (steady/turndown/shut-in) and **<1 %** (ramp-up) at production resolution;
  the hydrate water/gas sinks are included so the balances close *with* the consumption.
- Outflow-only / isolated-line outlet + **local-first** (nearest-neighbour) conservative
  bound-enforcement → exact conservation even in full **shut-in**, without acausal line-wide
  teleporting of liquid.
- **Well-posed pressure BCs** (default implicit engine): inlet pressure pinned **in-system**
  (Dirichlet) with an outlet through-flux (rate control) — genuine transient pressure with no
  post-hoc field re-anchoring; a **single dt per step** for the momentum and transport updates.
- **Bounded fallback**: on any non-finite step the implicit engine degrades that step to a
  quasi-steady update (`numerics.engine="quasisteady"` to force it), and that path always advances
  time. 0 fallbacks were triggered on the test cases.
- Liquid **and gas** mass-conservation error + **clip-activation diagnostics** reported every
  run (`[PASS]`/`[WARN]`) as a self-check (so "robust" never silently masks an instability).
- Stochastic-nucleation Monte-Carlo ensemble with an always-on modest parameter spread →
  genuine P10/P50/P90 bands (`numerics.deterministic=True` to disable; `--uq` to widen).
- Built-in **verification suite** (`--verify`), **input validation**, an automated
  **test suite** (`test_solver.py`), a pinned **requirements.txt**, and **calibration**
  (`--calibrate`).

---

## 6. STATUS — verification vs validation (read before use)

A genuine transient coupled-PDE solver with production-style, **verified** numerics:

- **Verification (the code solves the equations correctly): built in & passing.**
  `--verify` confirms the closures reproduce published reference values (hydrate
  equilibrium, Gregory–Scott slug frequency, Haaland friction, Hammerschmidt, drift-flux
  limits) and that the transient core conserves liquid, gas and hydrate mass. All checks pass.
- **Against a community benchmark, not only itself.** `shct_verification.py` runs
  **Ransom's water-faucet problem** — the standard assessment case for RELAP5-3D,
  MARS and TRACE, with a closed-form solution. The holdup transport reproduces it at
  an observed **L1 order of 1.04**, **6.1×** better than first-order upwind. Scope is
  stated, not glossed: this is drift-flux, so it cannot close the faucet's *momentum*
  problem (gas at rest, liquid in free fall); what is tested is the conservative TVD
  scheme that carries the holdup, which is the part that transports slugs.
- **Against published experiment.** `shct_evidence.py` checks the deposition model
  against six findings reported in flow-loop studies. Five are directional (steady-state
  thickness, subcooling, shear stripping, MEG, azimuthal skew) — the numeric series are
  paywalled and are not reproduced. The steady-state thickness is the notable one: the
  previous gated formulation could not have produced a plateau at all, so the
  most-reported observation in that literature corroborates precisely the change that
  made Φ_SH = 1 an emergent balance instead of a switch. The sixth is quantitative: the
  film growth rate against Qin (2020), 0.02–0.08 in/hr, which the wall-area growth law
  meets at 0.65–2.6× — the check that exposed the earlier interfacial-area form.
- **Validation (the constants match a specific reality): your calibration step.**
  The kinetic/coupling constants ship as literature-typical defaults; `--calibrate`
  fits them to **your** measured data so the solver is adapted to any fluid/field.
  Run it against your dataset before relying on absolute numbers.
- The Slug–Hydrate Coupling Number Φ_SH and the consolidation/plug mechanism are the
  central contribution — physically reasoned and mass-consistent, but their
  quantitative law still warrants experimental confirmation (flow-loop).
- **How to read a Φ_SH magnitude.** Φ_SH = C·Ψ, where Ψ = kg,wall·a_wall·ΔTsub,wall^n/f_slug
  is formed on the WALL area the deposit actually grows on, not the gas–liquid
  interfacial area a_i (which closes the BULK growth); the two were the same symbol
  until the wall-area correction, and Φ_SH has to be built from the rate the deposit
  follows or the balance below is not true of the model.

  Φ_SH **drives nothing in the solver.** It used to: consolidation required Φ_SH > 1, the
  wall-capture fraction was scaled by `clip(Φ_SH − 1, 0, 1)` so nothing deposited below 1,
  and erosion ran only below 1. That made Φ_SH = 1 a switch between "no deposit" and
  "unopposed deposit" — the criterion this work reports was an *input*, and no result the
  solver produced could have contradicted it. Deposition and slug scouring now simply
  compete, and Φ_SH is a diagnostic formed from the two rates.

  What follows is the part worth quoting. Setting d(δ)/dt = 0 gives a finite equilibrium
  thickness, and substituting Φ_SH = C·kg,wall·a_wall·ΔTsub,wall^n/f_slug collapses it to

  > δ_eq = Φ_SH · δ_ref,  δ_ref = f_wall·D / (4·C·k_ero)

  so **Φ_SH is the equilibrium deposit thickness in units of δ_ref** — 21.2 mm on the
  10.75-in case study — which is what `C` has always encoded without saying so. Growth runs
  away when δ_eq exceeds the consolidation restriction, because the deposit then locks and
  erosion stops. That happens above

  > Φ_crit = 2·C·k_ero·`consol_restriction` / f_wall  = **1.08** at the shipped constants,

  *derived* from three kinetic constants rather than asserted. It lands near 1, which is the
  finding; it is not 1 by construction, which is the point. The solver reports:
  - `max_Phi_SH_uncapped` — the true peak, so a quoted value is never the plot cap
    (`phi_report_cap`, default 1e4) in disguise;
  - `max_Psi_kinetic_ratio` — the C-free ratio Ψ, the part the model actually predicts and
    the quantity a future calibration should fit `C` against;
  - `Phi_SH_critical` and `deposit_ref_mm` — Φ_crit and δ_ref, so a reader can see what the
    threshold means and check it against a measurement;
  - `Phi_SH_above_critical_frac` — the fraction of WALL-SUBCOOLED cell-timesteps past Φ_crit,
    where the denominator is every cell-step with a positive wall subcooling (not only those
    that have nucleated; un-nucleated cells have f_wall = 0, so their local Φ_crit is
    infinite and they can never enter the numerator),
    i.e. past the point of no return. (This replaces `Phi_SH_gate_saturated_frac`, which
    counted cells above Φ_SH = 2 because the old gate saturated there. There is no gate now.)
  - `cooldown_to_hydrate_h` with `cooldown_source` — the no-touch time and how it was
    obtained. `transient` is a measured crossing into the hydrate region after the event;
    `lumped` is the lumped-capacitance estimate when no crossing occurs in the window;
    `already-subcooled` means the monitor was ALREADY inside the hydrate region when the
    event happened, so the no-touch time is zero. That third case is the as-operated and
    shut-in answer here, and it is a result rather than a missing value: a line that runs
    inside the hydrate envelope has no touch-free window to spend.
  - `dew_point_bar` — the water/hydrocarbon dew point at the monitor, or **null** where no
    dew point exists between 1 and 700 bar. On this live crude it is null at every monitor
    temperature: the flash vapour fraction is 0.63–0.71 at 1 bar and falls with pressure, so
    it never reaches the V → 1 target. The bubble point is defined and is reported by the
    EOS module. A number here that sits exactly on 1 bar or 700 bar is a bracket endpoint,
    not a dew point, and `check_outputs.py` now refuses it.
  - `slurry_rel_viscosity` with `slurry_visc_saturated` — the Camargo–Palermo relative
    viscosity, and whether it is pinned. The correlation diverges at the packing limit, so
    once φ reaches φ_max the magnitude is set by an internal 0.999 clip and by nothing else
    (3.16e7 = 0.001^−2.5 on the shut-in). When the flag is true the number is a bound, not a
    measurement, and should not be quoted as one.
  - `shear_margin_vs_deposit_strength` and `..._hi` — the sustained wall shear over the
    100 Pa and 200 Pa ends of Di Lorenzo's measured consolidated-deposit strength, so the
    margin is reported as the band the measurement actually is (0.66 and 0.33 as-operated).
  - `hydrate_scoured_frac` — hydrate moved from the wall deposit into the bulk phase field by
    slug scouring. Erosion once discarded this mass; it is now transferred and audited rather
    than lost. The fraction is duty-dependent, not a constant of the model: 25 % on the
    bundled default case, 50 % on the as-operated case study (which slugs over its whole
    length), 3.5 % on the shut-in and 0 % on the mitigated line, which forms no hydrate at
    all. This used to be quoted as a flat “~9 %”, which was true of none of them.
  - `hydrate_outflow_frac` — hydrate that a cell already at the packing limit handed
    downstream and that then left the pipe at the outlet. It closes the audit
    (formed = stored + scoured + clipped + out) and it is ~1e-9 on every scenario here: the
    cascade almost always finds headroom before the outlet. The accumulator existed but was
    never reported, so the balance could not be checked from the outputs.

  Run `python3 solver.py --sensitivity` to see how far each headline number moves when the
  four unfitted constants — `kg0`, `growth_exp_n`, `C_phi` and `f_slug_floor_Hz` — are
  swept across their plausible ranges.
- Individual metre-scale slugs are **sub-grid** on any production mesh (dx ≈ 200 m on the
  bundled 20 km default case, ≈ 460 m on the 32 km case study), so they are carried as slug
  statistics from correlations; terrain/void-wave dynamics and all transients are resolved.

**Appropriate use:** screening, design, scenario ranking, sensitivity and risk analysis
on a realistic case once calibrated to its data.

### How it now compares to OLGA / LedaFlow — and the honest remaining gap

**Closed:** transient coupled PDEs; terrain/inclination; multi-scenario operations
(steady/ramp/turndown/shut-in); multi-layer wall heat + cooldown time; thermodynamic-
inhibitor transport & design; slug statistics & loads; hydrate kinetics, deposition &
transportability; the coupled Φ_SH risk; probabilistic ensembles; exact mass
conservation; built-in V&V and data calibration.

**Where it is arguably MORE robust than OLGA/LedaFlow (robustness ≠ accuracy):**
- **Bounded fallback by design** — on any non-finite step the implicit engine degrades
  that step to a quasi-steady update so that time continues to advance; 0 fallbacks were
  triggered across the reported runs.
- **Conservation in every regime** — liquid, gas **and hydrate** mass to ~0% in steady,
  ramp-up, turndown and **full shut-in** (gas via a mass-consistent continuity equation; the
  hydrate water/gas sinks close the balance), with **clip-activation diagnostics** so "robust"
  never silently masks an instability.
- **Unconditionally well-posed** — the drift-flux slip closure avoids the ill-posed
  complex-characteristic instability of the bare two-fluid model; inlet/outlet BCs are
  well-posed (in-system inlet datum + outlet through-flux), with predictor–corrector `dt`.
- **Self-checking** — built-in V&V (correlation benchmarks, liquid/gas/hydrate conservation,
  **headline-output grid convergence**) reported every run.

**Where a certified OLGA/LedaFlow still leads (be clear with stakeholders):**
1. **Compositional PVT/EOS** — OLGA carries full multi-component thermodynamics; this solver
   uses black-oil-style correlations. *Partially closed:* a real-gas **Z(P,T)** correlation
   (`fluids.gas_Z_corr`) and a user **Z/PVT table** (`gas_Z_table`, `hyd_Teq_table`) are now
   accepted — but it is still not a multi-component flash.
2. **Droplet/film fields & fast acoustics** — the `twofluid` engine solves two INDEPENDENT
   phase momenta and an optional **droplet-entrainment** fraction is reported
   (`fluids.droplet_entrainment`), but OLGA additionally carries fully-coupled droplet and
   liquid-film fields and resolves fast acoustics; this solver acoustic-filters.
3. **3-phase oil/water slip** — single composite liquid by default; an **opt-in separate
   water transport with gravity settling** (`fluids.oil_water_slip`) now resolves water
   accumulation at low points, but it is a reduced-order model, not a full 3-field treatment.
4. **Decades of experimental validation** — OLGA is tuned against thousands of lab/field
   cases. This solver is *verified* and *calibratable*, and ships a **validation harness**
   (`--validate data.json`) that scores predictions and the Φ_SH plug criterion against YOUR
   measured data — but absolute accuracy still depends on that data, which only you can supply.

So: a verified, calibratable, robust OLGA-*style* engineering
solver for slug & hydrate flow assurance — arguably more robust (bounded fallback, always
conservative, well-posed) than a bare two-fluid code, while **not** matching OLGA's
compositional PVT, multi-field momentum, or breadth of validation. The path to full
parity is the programme set out in `report.pdf`.

---

## 7. v3 follow-up hardening (numerics, statistics, packaging)

Beyond the coupled-PDE completion above, this revision adds (defaults keep the verified
behaviour; most are opt-in):

- **Numerics:** 2nd-order **TVD** holdup/hydrate advection (`numerics.flux_limiter`),
  **implicit** wall-loss + **Heun** energy corrector (`numerics.splitting`),
  **predictor–corrector dt** (`numerics.substep_cfl_growth`), and **clip-activation
  diagnostics** with a warning threshold (`numerics.clip_warn_frac`).
- **Physics options:** real-gas **Z(P,T)** / Z-table (#1), opt-in **oil/water slip** with
  water-accumulation reporting (#2), a **conserved gas-mass continuity** equation with a
  drift-flux-vs-mass **consistency** diagnostic (#3), **deposit-insulation** self-limiting of
  late growth (`kinetics.k_dep_insul`, #7), reduced-order **slug-length** (#5) and
  **droplet-entrainment** (#6) statistics, and a **MEG basis** option (aqueous/stream).
- **Statistics:** UQ entries accept **distribution specs** (`{"dist":"normal|lognormal|
  uniform",...}`) and a shared-z **correlation** flag (#15); a **Kaplan–Meier** (right-
  censored) time-to-plug CDF (#17).
- **Calibration:** deterministic objective, production-grid validation, and a **parameter
  identifiability/sensitivity** report at the optimum (#18).
- **Packaging:** the code is split into **`shct_model.py`** (schema), **`shct_correlations.py`**
  (pure closures, independently testable), and `solver.py`; with **`requirements.txt`**,
  **`pyproject.toml`** (ruff + mypy + pytest), **`test_solver.py`** (closure + regression
  tests) and **`.github/workflows/ci.yml`**.

**Deeper PVT/3-phase (round 3):** real-gas **Lee viscosity** (`gas_visc_corr`), **oil density**
ρo(P,T) (`oil_pvt_corr`), and a drop-in **PVT property table** (`fluids.pvt_table` =
`[[P,T,ρ_oil,ρ_gas,μ_oil,μ_gas],…]`) that overrides the correlations, and the gas-mass continuity
stays **0%** for *any* density model. With `oil_water_slip` the local 3-phase liquid density now
feeds the **momentum gravity/friction & energy**, and hydrate draws its water from the **water
phase**; droplet entrainment lowers the **wall-film holdup**.

**Compositional EOS build (round 4):** a genuine **Peng–Robinson** equation of
state with multicomponent **vapour-liquid flash** ships in **`shct_eos.py`** — set
`fluids.composition = {"C1":0.83, "C2":0.07, …}` and the solver runs on EOS-computed densities,
Z-factors and viscosities (validated vs Standing-Katz/NIST). A first-generation **two-fluid-mass**
coupling (`numerics.volume_consistent_pressure`), a resolved **water-hammer / acoustic** option
(`numerics.acoustic`, with the Wood mixture sound speed reported) complete the build.
The published-reference scoring now lives in the solver itself — `--validate-closures` runs
the friction, drift-flux, slug-frequency and hydrate-curve checks against the data in
`validation/data/` and prints a combined summary. It also runs a flow-loop **holdup** check
when a void-fraction dataset is present, and **none ships**: the summary says so rather than
claiming a score it did not produce, and `--validate-flowloop` names the file it expects. (Earlier revisions pointed at
`6/validate_against_literature.py`, `6/sources.docx` and `UNIVERSALITY.md`; none of those is
part of this repository. The sources are cited in place, in each validator's docstring and in
`report.pdf`.)

A handful of items remain genuinely out of code scope (full compositional flash; fully-coupled
multi-field two-fluid mass with a volume-consistent pressure; resolving metre-scale slugs;
water-hammer acoustics; and **experimental validation**, which needs your lab/field data — feed it
via `--validate`). Section 6 above and the status table in `README.md` set out what is verified,
what is validated and what is neither, each with its code-side entry point.
