# slug-analysis-sam

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22259744.svg)](https://doi.org/10.5281/zenodo.22259744)

Archived on Zenodo. The badge above is the *concept* DOI and always resolves to the
newest release, which is
[10.5281/zenodo.22666728](https://doi.org/10.5281/zenodo.22666728) (SHCT v4.0.0) — the
archive this tree corresponds to, and the one to cite.

**Do not cite [10.5281/zenodo.22348213](https://doi.org/10.5281/zenodo.22348213) (v3.4.0)
for any number in this repository.** Two corrections since that deposit move every reported
quantity: the PVT surface returned the *liquid* density wherever the mixture was
single-phase, so the gas at the critical sections was about eight times too heavy, and the
mixture-velocity clip was setting the riser velocity rather than bounding it. v3.4.0 is
internally consistent — those files really do produce those numbers — but the numbers are
superseded. Releases before v3.4.0 archive superseded physics of a different kind again:
they deposit hydrate on the gas–liquid interfacial area rather than the pipe wall.

`.zenodo.json` holds the deposit metadata, and a published record can be edited by
hand on zenodo.org, so the two can drift. `python3 case/scripts/sync_zenodo_metadata.py
--check` reports any drift and `sync_zenodo_metadata.py` pulls the record back into the
file. **Run the check before tagging a release**, since the file — not the record — is
what a new version deposits. An Action
(`.github/workflows/zenodo-metadata-sync.yml`) also polls the record every five minutes —
GitHub's floor for a scheduled workflow — and commits the file when it moves, so a hand
edit on zenodo.org lands here within the hour.

**Transient, coupled-PDE prediction of hydrodynamic slugging and gas-hydrate
formation in subsea multiphase pipelines — the SHCT solver, plus a full deepwater
flow-assurance case study.**

_Author: **Akosa Samuel Onyejekwe**_

This repository contains an engineering-grade simulator (`solver.py` and the
`shct_*` extension modules) and a complete, reproducible case study of a deepwater
medium-crude-oil subsea tie-back, with every generated output (fields, time-series,
space-time maps, slug statistics, probabilistic risk, cross-section / quasi-3-D
reconstructions, compositional PVT) and the engineering reports that document them.

---

## 1. What the solver does

The **SHCT** (Slug–Hydrate Coupled-Transient) solver integrates, in time, a system
of coupled partial-differential equations for multiphase flow, heat transfer and
hydrate formation on arbitrary terrain, closed by published flow-assurance
correlations and a compositional Peng–Robinson PVT engine. It predicts, end-to-end:

- **Slugging** — hydrodynamic, terrain and severe-riser slugging; slug frequency,
  length, holdup, surge volume and loads.
- **Hydrates** — formation, wall deposition, consolidation and plugging, with a
  genuinely probabilistic (P10/P50/P90) time-to-plug.
- **The coupling between them** — quantified by the **Slug–Hydrate Coupling Number
  Φ_SH**, the central risk metric (Φ_SH > 1 ⇒ hydrate formation outruns slug
  scouring ⇒ plugging criticality).
- **Thermal & inhibitor design** — multi-layer wall heat transfer, effective U,
  cooldown / no-touch time, and required MEG dose (Nielsen–Bucklin).
- **Compositional PVT** — multicomponent Peng–Robinson vapour-liquid flash, real-gas
  Z-factor, Lee gas viscosity, Lohrenz–Bray–Clark liquid viscosity.

### Governing equations (summary)

| Eq. | Physics | Form |
|----|---------|------|
| **H** | Liquid-holdup transport (drift-flux kinematic wave) | `∂(α_l·A)/∂t + ∂(α_l·v_l·A)/∂x = −S_l` |
| **Gm** | Gas-mass continuity (conservative) | `∂(ρ_g·α_g·A)/∂t + ∂(ρ_g·V_sg·A)/∂x = −ṁ_(gas→hyd)` |
| **G** | Mixture momentum → implicit pressure | `u_m = U₀ − C_u·∂p/∂x`, tridiagonal Poisson pressure |
| **E** | Energy transport | `∂T/∂t + j·∂T/∂x = −U·(4/D)(T−T_sink)/(ρ_m c_p) + q_latent + q_JT` |
| **P** | Hydrate phase-field (advected reaction–diffusion) | `∂φ/∂t + v_l·∂φ/∂x = D_φ·∂²φ/∂x² + R_grow + R_nuc + ξ` |
| **C** | Slug–Hydrate Coupling Number | `Φ_SH = C·k_g·a_i·ΔT_sub^n / f_slug` |

The complete, transcribed equation set (governing PDEs **and** every closure) is in
**`report.pdf`** (Section 4 — Model equations).

**Published closures:** Bendiksen drift-flux, Taitel–Dukler regime map,
Gregory–Scott / Zabaras slug frequency, Haaland friction, natural-gas hydrate
equilibrium, CSMHyK-type growth, Camargo–Palermo slurry viscosity,
Nielsen–Bucklin inhibitor suppression/sizing, API RP 14E erosional limit,
Peng–Robinson EOS.

**Numerics:** conservative finite-volume with adaptive CFL stepping, optional
2nd-order TVD advection, implicit (tridiagonal/Thomas) pressure, implicit wall-loss
with a Heun corrector, a stochastic Monte-Carlo ensemble for P10/P50/P90 bands, and
a bounded fallback (any non-finite step degrades to a quasi-steady update so that time
continues to advance; 0 fallbacks were triggered across the reported runs). Liquid, gas
and hydrate mass conserve to ~0 %.

---

## 2. Repository layout

```
.
├── solver.py                     # the SHCT transient coupled-PDE solver
├── shct_model.py                 # input data model (dataclasses / case schema)
├── shct_correlations.py          # pure, unit-testable flow-assurance closures
├── shct_eos.py                   # Peng–Robinson EOS + multicomponent VLE flash
├── shct_crosssection.py          # reduced-order cross-section / quasi-3-D reconstruction
├── shct_compositional.py         # compositional / PVT tracking along the line
├── shct_compositional_sim.py     # compositional-transport (hydrate-former depletion)
├── shct_spacetime.py             # space-time / multi-time figure set (holdup vs distance at
│                                 #   successive times, true space-time fields, resolved slug
│                                 #   propagation & tracking, riser waterfall, cloud maps)
├── shct_threed.py                # 3-D field reconstruction + VTK export
├── shct_openfoam.py              # OpenFOAM (interFoam) coupling: case generation, running and
│                                 #   the drift-flux distribution parameter read off the CFD field
├── shct_verification.py          # verification against exact solutions (thermal relaxation,
│                                 #   order of accuracy, cross-engine, Ransom water faucet, MMS)
├── shct_evidence.py              # the deposition model against published flow-loop findings
├── shct_benchmark.py             # harness for comparing against a reference simulator
│                                 #   (OLGA / LedaFlow); ships no reference data
├── shct_style.py                 # the shared palette, the no-black/no-dark rcParams, and the
│                                 #   single figure export resolution every module draws at
├── test_solver.py                # closure + regression test suite
├── README_solver.md              # in-depth solver documentation
├── pyproject.toml / requirements.txt
├── report.pdf                     # THE comprehensive report (new medium-crude case): background,
│                                  #   problem, all equations, inputs, every output, validation, calibration
├── validation/data/               # credible published validation datasets (+ recorded sources)
└── case/                          # the deepwater medium-crude-oil case study
    ├── scripts/                   # active pipeline
    │   ├── run_case_study10.py    #   runs the 3 scenarios + advanced stack + validation
    │   ├── build_report.py        #   assembles the report (report.docx → report.pdf) from the outputs
    │   ├── build_reports.py       #   shared docx helpers + the full equation catalogue
    │   ├── make_animations.py     #   renders the transient GIF animations (per scenario)
    │   ├── run_sensitivity.py     #   parallel Phi_SH sweep over kg0 / growth_exp_n / C_phi
    │   ├── export_paper_figures.py#   assembles the numbered manuscript figure set from the outputs
    │   ├── check_outputs.py       #   inspects EVERY generated file (blank/collapsed figures,
    │   │                          #     non-finite or out-of-bounds table columns, metric bounds)
    │   ├── docx2pdf_safe.py       #   docx->pdf via Word COM; refreshes SEQ/TOC fields first
    │   └── _paths.py              #   shared layout + no-black/no-dark style hook
    ├── outputs_steady/            # (A) as-operated normal production
    ├── outputs_shutin/            # (B) unplanned shut-in cooldown
    └── outputs_mitigated/         # (C) engineered mitigation (insulation + MEG)
```

---

## 3. The case study

A representative deepwater medium-crude-oil subsea tie-back — **32 km, 10.75-in
carbon-steel flowline + steel catenary riser, ~1100 m water depth** — carrying a
~30° API medium crude oil (C1 ≈ 43 mol%, ~31 mol% C7+ tail) over a cold (4 °C),
undulating seabed, at the late-life duty of 70 % water cut and 0.6× design rate
(see the v3.4.0 note below for why). This geometry and fluid is a textbook combination for
**both** slugging and hydrates, so it exercises the whole prediction chain.

Three scenarios are run end-to-end through the real solver:

| Scenario | Folder | Description |
|----------|--------|-------------|
| **A — as-operated** | `case/outputs_steady/` | normal production, degraded (water-flooded) insulation, no inhibitor → the high-risk prediction |
| **B — shut-in** | `case/outputs_shutin/` | unplanned shut-in cooldown → no-touch time |
| **C — mitigated** | `case/outputs_mitigated/` | restored multi-layer insulation + continuous MEG → risk removed: 0 % plug, no deposit, no under-inhibited length (design tool) |

**Headline result (as-operated):** intermittent flow over the whole line with slugs
up to ~78 m; the cold under-insulated wall drives the fluid 6.4 °C into the hydrate
region — peak Φ_SH 0.34, sustained 0.29 against the derived Φ_crit = 1.08 — so the
line is **sub-critical and does not plug**: 0 % plug probability, a 4.1 mm peak wall
deposit, and no reach above Φ_SH = 1. The model sizes the inhibition at 30.2 wt% MEG
over a 20.1 km under-inhibited length.

**The hazard is the unplanned shut-in, not production.** Once the flow stops the
interface stops being renewed, and the same line plugs in 11 of 12 realisations at a
P50 of 16.1 h, with 34 % of the route above Φ_crit over 23.3 km and the bore closed to
the 117 mm full-bore cap. The engineered insulation + MEG fix removes the subcooling
entirely (peak deposit 0.0 mm, 0 % plug probability) and buys a 65.9 h no-touch time
at an effective U of 3.45 W/m²K.

These are the numbers from the outputs in this tree, regenerated against the corrected
solver. The much larger hydrate numbers that earlier versions of this README quoted are
artefacts of the two defects described immediately below.

> [!WARNING]
> **The tracked `case/outputs_*` have been regenerated against the corrected solver. They no
> longer match the v3.4.0 archive; the numbers in this tree are the corrected ones, are
> archived as v4.0.0, and the numbers on the v3.4.0 record are artefacts of the defects
> described below and of the fluid-model defect described at the top of this file.**
>
> **The hydrate curve in those outputs is ~20 °C too high.** The gas specific gravity was
> read from an EOS flash at inlet conditions, where this live oil is undersaturated and the
> flash returns no vapour, so the value used was the whole feed's (1.7645) rather than a gas
> gravity (0.6307). `hydrate_equilibrium_T` applies `18*(sg-0.60)`, so every cell carried
> about +20.4 °C. Fixed in `solver.py`, and the case study re-run. Corrected, the as-operated
> line does not plug: subcooling 24.4 → 6.44 °C, P_plug 1.00 → 0.00, peak deposit 117 mm →
> 4.1 mm, MEG 60 → 30.2 wt%. **The hazard moves to the unplanned shut-in, which plugs 11 of 12
> realisations at a P50 of 16.1 h with 34 % of the route above Φ_crit over 23.3 km.**
>
> **Φ_crit does not bound the deposition.** Bulk hydrate formed above the slurry packing limit
> is returned to the wall as deposit — a channel outside the wall-growth-versus-scouring
> competition that Φ_SH measures, so Φ_crit does not limit it. Across 18 duties with the
> gravity corrected, Φ_SH never reaches 1.08 (highest 1.038) yet 9 of 18 close the bore, all
> below threshold, in every realisation. This is the unattributed mechanism behind
> `hydrate_packing_clip_frac`. **Fixed after v3.4.0 was archived**: the excess is now handed
> downstream, cascading while there is headroom and leaving the pipe at the outlet, so mass is
> conserved by transport instead of by deposition (`reject_mode`, default `advect`; `plate`
> reproduces the archived behaviour). The decisive test was switching wall capture off
> entirely — with `f_wall = 0` and therefore zero wall growth, the archived code still closed
> the bore to `delta_max` on the same schedule, so the plug it reported contained no wall
> growth at all. With the fix the same duty settles at 15.2 mm with the bore open, and at
> `f_wall = 0` the deposit is zero.
>
> A residual remains and is stated rather than smoothed over: on the reference duty, three of
> eight realisations close above Φ_crit as the theory requires and three correctly stay open
> below it, but two close marginally early, at Φ_SH 0.90 and 1.04 against 1.08. Closure
> locks irreversibly once the deposit passes the consolidation restriction, so a transient
> excursion is enough. The gross artefact is gone; a near-threshold band is not.
>
> The numerical core is unaffected — balances close to 2.5e-15 / 3.3e-18, the five
> exact-solution checks pass, six published trends are reproduced, 129/129 tests pass.

> **v3.4.0 — hydrate deposits on the WALL area, and the case study moves to late life.**
> The wall growth law used `a_i`, the gas–liquid interfacial area. That is the right term for
> growth at that interface and the wrong one for a wall process: `a_i` falls toward zero as a
> line fills with liquid, so the model predicted almost no wall deposit in exactly the
> liquid-full, oil-dominated configuration where Qin (2020) *measures* film growth at
> 0.02–0.08 in/hr. It now uses the wall area, `4/D` scaled by the liquid holdup and the water
> fraction of that liquid, which puts the model at 0.65–2.6× the measured rate instead of
> 6.5–26×. Φ_SH is computed from the same area, so the coupling number and the process it
> describes finally refer to the same surface.
>
> With that correction 35 % water cut at full rate is **sub-critical** (Φ_SH = 0.27 against
> Φ_crit = 1.08) and does not plug, so the case study moved to late-life conditions — 70 %
> water cut at 0.6× design rate — which is the duty it still runs.
>
> *What v3.4.0 reported at that duty* — peak Φ_SH 1.95, sustained 1.06, 1.37 km above the
> Φ_SH = 1 contour, P50 3.72 h, peak deposit 117 mm, max subcooling 24.4 °C, engineered fix
> P_plug 0.33 — **is superseded**: those figures carry the gas-gravity defect described in the
> warning above. The current numbers for the same duty are in §3. Verification passes 5/5,
> six published trends are reproduced, and 129/129 tests pass.

> **Solver corrections in v3.2.0 — read the numbers from this release.** Three defects
> in the previous release moved every velocity-derived quantity. (i) The condensation
> latent-heat term bound its lookup indices to the names `i` and `j`, and `j` is the
> mixture volumetric flux, so the velocity field was overwritten by an integer
> temperature index on every step — corrupting slug length, the erosional margin, the
> interfacial area and the advection of temperature. (ii) The liquid-bounds enforcement
> destroyed liquid once the hydrate deposit closed the bore, which surfaced as an
> unexplained 5.9 % conservation error on a 48 h run that plugs — 0.000 % on a 12 h run
> that does not, at every timestep tested, so not a resolution problem. The loss is now
> measured and reported as `liq_bounds_discard_frac`, and the balance closes to ~1e-11 %.
> **Read that term as a limit of the model, not a solver defect:** once the deposit shuts
> the bore the pipe cannot hold the liquid arriving, and a one-dimensional model carries
> no representation of the pressure that would build behind a closing plug. For the
> as-operated case 5.92 % of the injected liquid (≈ 563 m³) had nowhere to go and was
> dropped at the bounds. It was previously invisible. **That discard is gone as of
> v3.4.0** — the as-operated case reports `liq_bounds_discard_frac` = 1.2e-15 and the
> liquid balance closes to 2.5e-15. (The bore no longer shuts on this case at all: with
> the gas gravity corrected the line is sub-critical, so the discard has nothing left to
> discard.)
> (iii) The slug-length statistics averaged in the correlation's 5000 m "not slugging"
> ceiling. Slug lengths, velocities and the erosional check should be taken from this
> release rather than the last; the hydrate and coupling results are unchanged in
> character.

> **The two-fluid description is well posed over 96 % of the route, and is not everywhere.**
> `27_wellposedness_map.png` reports the slip against the inviscid Kelvin–Helmholtz limit at
> which the one-dimensional two-fluid model loses hyperbolicity. The margin peaks at **1.96**
> and exceeds 1 over **4.3 %** of the route at the final state — a short reach near the riser,
> where the film is thin and fast. Over that reach the initial-value problem is ill-posed and
> the growth rate is grid-dependent, so slug activity localised *there* should not be read as
> a property of the flow; everywhere else it can be.
>
> This claim previously read "the margin peaks at ~0.86 and never reaches 1". That was an
> artefact of the limit itself: `dA_l/dh` is the chord width at the liquid LEVEL, and the code
> was passing the area fraction in its place, which understates the chord on a thin film
> (0.44 D against a true 0.64 D at α_l = 0.05) and so overstates the limit. Corrected — the
> exact circular-segment inversion the cross-section module already used — the margin is what
> it is.

> **Read these magnitudes with care.** This block used to warn that a 60 wt% MEG requirement,
> a 3.72 h P50 and a ~0 h no-touch time sat beyond reported field experience. With the
> gas-gravity defect fixed the as-operated figures are unremarkable — 30.2 wt% MEG, inside the
> 20–50 wt% band of normal continuous dosing, and no plug at all — so that warning no longer
> describes this case and has been withdrawn rather than left to lend the old numbers weight.
>
> The caution that does still apply is about provenance, not size: every kinetic and coupling
> constant (C_φ, n, k_g0, the slug-frequency floor) is a literature-typical value fitted to no
> dataset, so the ABSOLUTE magnitude of Φ_SH, of the time-to-plug and of the required dose
> inherits whatever uncertainty those four carry. `--sensitivity` (see §5) measures that
> inheritance directly. Treat them as model outputs, not as calibrated predictions.

> **Data provenance (honest framing):** the field is a representative *industrial
> archetype*. Geometry, fluid and operating parameters are realistic,
> self-consistent, literature-typical values for deepwater medium-crude-oil tie-backs —
> **not** proprietary operator data. The physics and predictions are produced by the
> real solver; the hydrate thermodynamics are anchored to published data.

---

## 4. Outputs

Each scenario folder contains the full output set:

- **Tables (CSV):** `fields_profile.csv` (along-line profile), `timeseries_monitor.csv`
  (transient history), `probabilistic_summary.csv` (P10/P50/P90),
  `engineering_deliverables.csv`, `feed_composition.csv`, `input_data_deck.csv`.
- **Metrics (JSON):** `summary.json`, `key_metrics.json`, `case_config.json`.
- **Charts (PNG):** profiles, the transient liquid-holdup field α_l(x,t), P–T vs the
  hydrate envelope, the Φ_SH(x,t) coupling-criticality map, slug prediction, deposit
  growth, probabilistic time-to-plug, diagnostics, cross-section / quasi-3-D
  reconstructions, compositional PVT, and the mitigation comparison.
- **Space-time / multi-time set (PNG)** — the figure family the transient-multiphase
  and flow-assurance literature uses to present a transient pipeline calculation, so
  the case study can be read directly against published work. Figures 23–25 adopt the
  distributed-sensing (DTS/DAS) waterfall convention — distance against time —
  because the fields they show (temperature, its gradient, and flow unsteadiness)
  are exactly what a fibre installed on such a line measures:

  | File | What it shows |
  |------|---------------|
  | `14_holdup_multitime.png` | liquid holdup along the whole route at six successive times, early transient and late quasi-developed state |
  | `15_slug_growth_propagation.png` | resolved slug units over a short reach at three successive times, one front tracked across the panels (T_b, X_b) |
  | `16_slug_train_waterfall.png` | slug tracking in the space-time plane: waterfall, semblance vs trial celerity, moveout-corrected waterfall, distance-stacked trace |
  | `17_hydrate_distribution.png` | in-pipe volume fractions along the line (unconverted water, hydrate in the liquids, wall deposit) + gas/oil/water rates into the host |
  | `18_shutin_profile_deposit.png` | late-time P, T vs T_eq and water holdup along the line + deposit volume fraction at successive times |
  | `19_spacetime_fields.png` | the **true space-time solution**: α_l, p, u_g, u_l, ΔT_sub and the wall deposit, each as a filled-contour field over (distance, time) |
  | `20_holdup_durations.png` | holdup along the pipeline after successive shut-in (or production) durations |
  | `21_riser_depth_time.png` | riser depth–time waterfall — slug boundaries during upward motion, their trajectories and the slug unit length |
  | `22_cloud_maps.png` | pipeline cloud maps at successive times: bore phase distribution above the bulk-temperature field, shared scale |
  | `23_dts_thermal_waterfall.png` | distributed-temperature waterfall T(x,t) with the monitored pressure overlaid, the operating stages marked and the hydrate-onset distance annotated |
  | `24_temperature_gradient.png` | temperature-gradient waterfall ∂T/∂x(x,t) — a travelling thermal front is a narrow band of steep gradient, so this localises it where the temperature map itself looks smooth |
  | `25_das_flow_noise.png` | flow-noise waterfall \|∂α_l/∂t\|(x,t) — where the holdup changes fastest is where the flow is most unsteady, with the intermittent reach and the riser base marked |
  | `26_parameter_panels.png` | pressure, temperature, holdup and mixture velocity along the route, each at the same successive times |
  | `27_wellposedness_map.png` | the two-fluid well-posedness (Kelvin–Helmholtz) boundary over the (V_sg, V_sl) plane with the case's own states, and the margin along the route |

  > **Resolved-slug figures — what is computed and what is reconstructed.** The transport
  > grid is `dx ≈ 460 m` while a slug unit is ~10–40 m, so individual slugs are a
  > *sub-grid* quantity the solver carries statistically (slug frequency `f_slug`, unit
  > length `L_u = V_t/f_slug`, slug-body holdup `α_ls`). Figures **15, 16 and 21**
  > therefore render a **kinematic reconstruction built entirely from those solver
  > outputs**: at every station the reconstructed square wave has the solver's local
  > slug frequency and translational celerity, and the body/film split is solved so
  > that the unit-averaged holdup reproduces the solver's cell-average `α_l` *exactly*
  > (mass-consistent by construction). Period, celerity, length and holdup are all run
  > outputs — nothing is assumed, and each figure states this on its face. A slug train
  > only exists while the line is actually flowing intermittently, so these three are
  > built from the latest snapshot at which it is (for the shut-in scenario, a state
  > before the line stops — the time is printed on the figure) and are skipped, with a
  > stated reason, when no such state exists. Every other figure in the set is plotted
  > directly from the solver's space-time history.
- **Snapshot archive:** `spacetime_state.npz` — the space-time history the figures above
  read (holdup, pressure, temperature, subcooling, deposit, phase velocities, regime and
  slug frequency, on the snapshot cadence). `shct_spacetime.rerender(folder)` rebuilds the
  whole figure set from it, so a figure can be restyled or rescaled without repeating the
  transient.
- **Animations (GIF):** `anim_flow_line.gif` (slugs travelling along the terrain-following
  pipe — liquid holdup α_l), `anim_crosssection.gif` (the pipe bore filling and the hydrate
  deposit ring closing toward a plug at the monitor), `anim_PT_cooldown.gif` (the monitor
  P–T point crossing the hydrate-stability envelope in time), `anim_riser_cycle.gif`
  (the riser-region monitor α_l–P trajectory — repeating loops = intermittent/slug flow,
  a settled point = stable flow), and
  `anim_profile_wave.gif` (the P(x,t) & T(x,t) cooling/pressure wave marching along the line).
  These are a supplementary visualisation layer — the transient story the static charts above
  capture as single frames.

**`report.pdf`** (repo root) assembles all of these into a single comprehensive
report — background, problem statement, the case study, every model equation, the
full input deck, every generated output (all metrics, CSV tables, per-CSV graphs and
the complete chart/curve/contour/map gallery), the published-data validation and the
calibration, with all sources recorded. Nothing is left out.

---

## 5. Usage

```bash
pip install -r requirements.txt

python3 solver.py                      # bundled real case
python3 solver.py --scenario shutin    # shut-in cooldown / hydrate-risk transient
python3 solver.py --engine twofluid    # full two-fluid (two independent phase momenta)
python3 solver.py --meg 30             # inject 30 wt% MEG inhibitor
python3 solver.py --config case.json   # any user case
python3 solver.py --verify             # verification: closures vs published values + mass conservation
python3 solver.py --sensitivity        # one-at-a-time sensitivity of Phi_SH, time-to-plug, MEG dose
                                       #   and deposit to the four ASSUMED constants kg0, n,
                                       #   C_phi and f_slug_floor_Hz
python3 solver.py --calibrate t.json   # validation: fit free constants to measured data

pytest test_solver.py                  # run the test suite

python3 case/scripts/check_outputs.py  # inspect every generated output: blank or
                                       #   colour-collapsed figures, non-finite or
                                       #   out-of-bounds table columns, metrics
                                       #   outside their physical bounds
```

A case is fully described by the JSON groups `pipeline`, `fluids`, `operating`,
`kinetics`, `numerics`, `scenario` (run `--dump-config` for an editable template).

---

## 6. Status — verification vs validation

- **Verification (the code solves the equations correctly):** built in and passing.
  `--verify` confirms the closures reproduce published reference values and that the
  transient core conserves liquid, gas and hydrate mass.
- **Validation (the constants match a specific reality):** the kinetic/coupling
  constants ship as literature-typical defaults; `--calibrate` fits them to *your*
  measured data. Run it against your dataset before relying on absolute numbers.

### What is and is not validated — read this before citing a number

| Element | Status | Against what |
|---|---|---|
| Haaland friction closure | **verified** | Colebrook–White (1939); 0.62 % RMS deviation |
| Slug-frequency closure | **verified** | reproduces Zabaras (2000) to machine zero |
| Drift-flux parameters | **verified** | Dumitrescu (1943), Bendiksen (1984) source values |
| Drift-flux slip vs 3-D CFD | **validated** | OpenFOAM v2406 interFoam, streamwise-periodic pipe, k–ω SST; distribution parameter measured 1.146 against the closure's 1.172 — 2.3 % |
| Hydrate equilibrium curve | **validated** | Deaton & Frost (1946) measurements; 1.72 °C RMSE |
| Mass conservation (liquid, gas) | **verified** | liquid 2.5e-15, gas 3.3e-18; bounds discard 1.2e-15 |
| Hydrate mass conservation | **partial — measured** | zero loss unless the bore plugs; 1.1 % unplaceable in plugged cells (shut-in), reported as `hydrate_packing_clip_frac` |
| Two-fluid well-posedness | **partial — measured** | inviscid Kelvin–Helmholtz limit; margin peaks at 1.96, above 1 over 4.3 % of the route |
| Holdup transport vs Ransom water faucet | **verified** | exact solution; observed L1 order 1.04, 6.1× better than upwind |
| Lumped thermal relaxation | **verified** | analytical decay; 0.0797 % NRMSE |
| Order of accuracy | **verified** | three-level refinement; observed order 0.995 on outlet T |
| Cross-engine agreement | **verified** | drift-flux vs two-fluid; holdup to 3.3e-3 |
| Smooth manufactured solution | **verified** | limited transport at order 1.17; 8.8× upwind, 3.9× Euler |
| Deposition trends vs published flow-loop findings | **corroborated — 6/6** | plateau, subcooling, shear, MEG, azimuthal skew (qualitative); film growth rate against Qin (2020), 0.65–2.6× measured (quantitative) |
| Φ_SH dimensional consistency | **verified** | dimensionless for any *n*; invariance test in the suite |
| Φ_SH criterion is *derived*, not imposed | **verified** | Φ_SH drives no term; Φ_crit = 1.08 follows from `C`, `k_ero`, `consol_restriction`; four tests |
| **Φ_SH magnitude and the value of Φ_crit** | **NOT validated** | *no dataset* |
| **C, n, k_g0** | **NOT fitted** | literature-typical values only |
| **Whole-system prediction vs a reference simulator** | **NOT yet run** | see below |

The central proposal of this work — that the competition between hydrate
deposition and slug renewal is captured by a single dimensionless group, and that
a threshold in that group separates scoured from plugging-critical — rests on
physical reasoning, dimensional consistency and internal consistency. **It has not
been tested against experiment.** Treat it as a hypothesis with a solver behind it,
not a calibrated predictor, and read §"Read these magnitudes with care" above
alongside it.

One thing did change, and it changes what kind of claim this is. Φ_SH was
previously wired into the deposition in three places — consolidation required
Φ_SH > 1, wall capture was scaled by `clip(Φ_SH − 1, 0, 1)`, and erosion ran only
below 1 — so the threshold was an **assumption of the model wearing the clothes of
a result**, and no output the solver produced could have contradicted it. Those
three are gone. Deposition and slug scouring now compete continuously, Φ_SH drives
no term, and the balance d(δ)/dt = 0 yields

> δ_eq = Φ_SH · δ_ref  with  δ_ref = f_wall·D / (4·C·k_ero) = **21.2 mm** here,

i.e. Φ_SH is the equilibrium deposit thickness in units of δ_ref — the physical
content `C` carried all along. Runaway begins where that equilibrium passes the
consolidation restriction and the deposit locks, at

> Φ_crit = 2·C·k_ero·`consol_restriction` / f_wall = **1.08**,

computed from three kinetic constants. That it lands within 8 % of 1 is a result rather
than a definition. The figures that used to sit here — plug probability 100 %, peak
deposit 117 mm, P50 moving 2.8 → 3.2 h — recorded what *that* rewrite did to the
then-current case, and all three are superseded by the gas-gravity correction: the
as-operated line no longer plugs at all (§3), so there is no P50 to quote for it, and the
plugging case is now the shut-in at a P50 of 16.1 h. The threshold itself is unchanged,
because it is computed from constants the correction did not touch.

It is still **not validated** — a derived threshold is falsifiable, which is not the same
as confirmed. The experiment below is what would settle it.

### Benchmarking against a reference simulator

`shct_benchmark.py` runs SHCT against OLGA, LedaFlow or any transient multiphase
code on an identical case and reports MAE, RMSE, normalised RMSE and the worst
deviation for holdup, pressure, temperature and mixture velocity, with a
comparison figure.

**No reference dataset ships with this repository.** There is no licence for such
a tool in the development environment, and fabricated benchmark numbers would be
worse than none — a made-up agreement is indistinguishable from a real one until
somebody tries to reproduce it.

What *can* be done without a licence, and is done: **Ransom's water-faucet problem**
(`shct_verification.py`). It is the standard assessment case that RELAP5-3D, MARS
and TRACE are themselves published against, and it has a closed-form solution, so
agreeing with it puts this solver's transport on the same yardstick as those codes
rather than one of its own devising. The holdup transport scheme reproduces the
exact liquid-fraction profile at an observed L1 convergence order of **1.04**, and
**6.1×** more accurately than first-order upwind at the same resolution. The scope is
stated rather than glossed: this is drift-flux, carrying one mixture momentum
equation, so it cannot close the faucet's *momentum* problem (which requires the gas
to stay at rest while the liquid falls freely). What is tested is the conservative
TVD scheme that carries the holdup — the production code path, not a
reimplementation written to pass — which is the part that transports slugs.

For a whole-system comparison, export your own reference run in the schema
documented at the top of `shct_benchmark.py` (the geometry, fluid and boundary
conditions are all in `case/outputs_*/input_data_deck.csv` and
`feed_composition.csv`), then:

```bash
python3 shct_benchmark.py validation/data/olga_asoperated.json
```

The loader refuses a file that does not name the tool that produced it, so a
benchmark in this repository always carries its provenance.

### A hydrate-mass loss that is measured rather than hidden

Scoured wall deposit is transferred into the bulk phase field rather than discarded, and
what a cell already at the packing limit cannot hold is handed **downstream** — cascading
while there is headroom and leaving the pipe at the outlet (`reject_mode`, default
`advect`). Only what the domain can neither carry nor pass on is genuinely lost, and that
is reported rather than absorbed.

On the two scenarios that do **not** plug it is exactly zero: as-operated and mitigated
both report `hydrate_packing_clip_frac` = 0.0. On the **shut-in**, which plugs 11 of 12
realisations, it is **1.1 %** — confined to cells where the wall has reached `delta_max`
*and* the bulk has reached `phi_max`, both full, in cells the model has already declared
solid. It does not touch the engineering answer: the line blocks at a P50 of 16.1 h,
long before those cells saturate.

Earlier versions of this section quoted 26.3 % and attributed it to the as-operated case.
That was the `plate` behaviour, which returned the excess to the wall as deposit; the
figure and the scenario both moved when the excess was handed downstream instead.

The liquid balance is unaffected and still closes to roundoff — 5e-15 on this case.

### Corroboration against published flow-loop findings

`shct_evidence.py` tests the deposition model against findings reported in the
published flow-loop literature, rather than only against itself:

| # | published finding | source |
|---|---|---|
| E1 | a deposit reaches a **steady-state thickness** rather than growing without limit | [1] |
| E2 | wall temperature strongly affects **both** growth rate and steady-state thickness | [1], [2] |
| E3 | higher velocity/shear reduces the surviving deposit; past a critical thickness it sloughs | [1], [2] |
| E4 | MEG reduces the steady-state deposit thickness | [1] |
| E5 | deposition is azimuthally non-uniform — fast at the liquid-wetted invert, slow at the gas-swept crown | [1] |
| E6 | **measured film growth rate**, 0.02–0.08 in/hr in a liquid-full, oil-dominated loop — the model returns 0.65–2.6× that band | [3] |

[1] X. Zhang, E. O. Straume, G. A. Grasso, R. E. M. Morales, A. K. Sum, *Fuel* **262**
(2020) 116558, [doi:10.1016/j.fuel.2019.116558](https://doi.org/10.1016/j.fuel.2019.116558).
[2] Z. M. Aman *et al.*, *J. Nat. Gas Sci. Eng.* **35** (2016) 1096–1103,
[doi:10.1016/j.jngse.2016.05.015](https://doi.org/10.1016/j.jngse.2016.05.015).
[3] H. Qin, *Hydrate film growth and risk management in oil/gas pipelines using
experiments, simulations and machine learning*, PhD thesis, Colorado School of Mines
(2020), open access.

**Read this for exactly what it is.** E1–E5 are *directions*, not magnitudes. The
numeric deposit-thickness series in those papers are paywalled; they are **not
reproduced here**, and inventing them would be worse than having none. Each check
runs the real solver end-to-end and asks whether the trend survives the full
coupling — it is corroboration, and it cannot rescue a wrong magnitude. **E6 is the
exception**: it is the one quantitative deposition rate found in the open literature,
and it is what exposed the interfacial-area form of the wall growth term (the model
sat at 6.5–26× the measured rate) and then supported the wall-area form that replaced
it (0.65–2.6×). It anchors the growth magnitude; it does not fit C, n or k_g0.

E1 is the one worth pausing on. **The previous formulation could not have produced
it.** With deposition gated by `clip(Φ_SH − 1, 0, 1)` and erosion running only below
Φ_SH = 1, a cell either grew with nothing opposing it or decayed to bare wall; a
steady-state thickness was not in the model's vocabulary. It appears on its own once
the two rates compete. So the single most-reported observation in the flow-loop
literature is reproduced by exactly the change that removed the Φ_SH = 1
circularity — the model is now answerable to a measurement it was previously
incapable of contradicting.

```bash
python3 shct_evidence.py case/outputs_steady
```

### One side of the coupling is now anchored to a measurement

Φ_crit is derived from three constants, and until now none of them had been measured.
One now has been, and it constrains the model from outside.

Di Lorenzo, Aman, Kozielski, Norris, Johns & May (*J. Chem. Thermodyn.* **117** (2018)
81–90, [doi:10.1016/j.jct.2017.08.038](https://doi.org/10.1016/j.jct.2017.08.038))
determined the effective shear strength of a **consolidated** hydrate deposit in situ,
from sloughing events in a gas-dominant flow loop: **100–200 Pa**. The solver now
computes the wall shear stress the line actually raises and reports it against that
number:

| | |
|---|---|
| wall shear stress, as-operated (mean / sustained / startup peak) | **7.2 / 66 / 75 Pa** |
| measured consolidated-deposit shear strength | **100–200 Pa** |
| margin on the sustained figure (`shear_margin_vs_deposit_strength`) | **0.66** |

**This table used to read 4.4 / 14.2 Pa and a margin of 0.14 — "seven-fold short".** Those
figures predate the move to 70 % water cut: τ goes as ρ_m·j², and the liquid density rises
from 858 to 975 kg/m³ with the water. The margin is 0.66, not 0.14 — and it fell from ~0.8
when the PVT gas density was corrected, because a gas eight times lighter carries less of
the mixture momentum.

**The bound still holds for an operable line, and this line is not operable at the riser.**
Checked directly against the friction closure at the current fluid: at the API RP 14E
erosional limit for this case — **5.69 m/s** — the wall shear reaches only **63 Pa**, below
the measured 100 Pa lower bound, and **7.25 m/s** would be needed to reach 100 Pa.
The model's own peak mixture velocity is **8.13 m/s**, above that limit — but see the
grid-convergence note below before reading a ratio off it. So the correct statement is narrower than the one this
section used to make: *within* the erosional envelope, flow cannot strip a consolidated
deposit; the riser of this case study is predicted to run outside that envelope, and there
the shear does reach the measured strength.

**That exceedance is one to two cells wide at the riser base, and it is not grid-converged.**
`Vm_peak_mps` is a point maximum taken next to a flow reversal, which is the least
convergent statistic the model produces. Refined from 70 to 210 cells on the as-operated
case it reads 8.05, 9.59, 7.23, 6.81 m/s, its location moves 2 km, and the ratio to the
limit swings between 1.19 and 1.68. Extent does not rescue it — that was tried, on the
expectation that a length would be mesh-independent, and measured across 70/105/140 cells
the route-length over the limit is the *worst*-behaved of the candidates:

| statistic | 70 | 105 | 140 | spread |
|---|---|---|---|---|
| peak velocity (m/s) | 8.05 | 9.59 | 7.23 | 1.33× |
| 99th percentile of route (m/s) | 6.34 | 7.22 | 4.69 | 1.54× |
| **95th percentile of route (m/s)** | **2.88** | **2.53** | **2.50** | **1.15×** |
| route over the limit (km) | 0.46 | 0.61 | 0.23 | 2.67× |
| peak of a 1 km running mean (m/s) | 6.83 | 4.13 | 4.45 | 1.65× |

So the honest statement is narrower than a design finding. **The line as a whole sits well
inside its erosional envelope** — the 95th percentile of route length is 2.5–2.9 m/s against
a 5.69 m/s limit, and that is the one quantity here stable under refinement. **A short reach
at the riser base is predicted above the limit on every grid tried**, which is a real flag,
but its magnitude and its extent are both properties of the mesh at this resolution. It
warrants a locally refined study, not a number. `Vm_peak_mps`, `erosional_exceedance_km` and
`erosional_exceedance_frac` are reported in every summary; read them with this table.

Two consequences, and the second corrects something this project previously implied.

1. **`locked` is released where the shear reaches the measured strength — and it turns out
   not to matter.** A cell is now freed when its own wall shear reaches `tau_deposit_Pa`,
   which is the condition Di Lorenzo measured sloughing at, so the flag follows the
   measurement instead of an assumption. But it gates `d_ero ∝ f_slug·δ·(~locked)`, and
   consolidation needs δ > 27.4 mm: those conditions do not overlap on any duty this model
   reaches. A flowing line peaks at 4–15 mm and never consolidates; a shut-in consolidates
   to 85 mm but has stopped flowing, so `f_slug` sits on its 1e-4 Hz floor. Measured, a
   100× change in the release strength moves the peak deposit by **0.037 %**. `locked` was
   documented for a long time as a live modelling assumption; at that size it is not one,
   and the suite pins the inertness so a future duty that changes it gets surfaced.
2. **The erosion term is not mechanical stripping of consolidated deposit.** Over the
   flowline the shear is an order of magnitude short of the measured strength (mean
   7.2 Pa against 100–200 Pa), so what the term removes is *nascent, weakly-adhered*
   deposit, and adhesion prevented before consolidation. Text throughout this project
   described slugs "scouring" and "shearing away" the deposit, which reads as the
   stronger claim; over the flowline the measurement rules that reading out.

This does **not** measure Φ_crit. It measures one constant on one side of the balance,
and it makes the model answerable to a number it did not choose. The threshold itself
still requires the experiment below.

### The measurement that would settle it

Φ_SH is falsifiable, and cheaply — and now it is falsifiable in two independent
ways rather than one, because removing the gate turned a switch into a prediction
about a *thickness*.

1. **The plateau.** At fixed subcooling and sub-critical Φ_SH, the model predicts
   the wall deposit to stop growing at δ_eq = Φ_SH · δ_ref rather than continue —
   a number, at a measurable place, that a flow loop either sees or does not.
   Sweeping f_slug should trace δ_eq ∝ 1/f_slug through that plateau.
2. **The threshold.** Growth should stop plateauing and run away as Φ_SH crosses
   Φ_crit ≈ 1.08, at a deposit restriction of `consol_restriction` ≈ 18 % of bore.

The first test is the more searching of the two: it checks the *form* of the
balance, and it does not require reaching plugging conditions. Both need only wall
deposit thickness against time, at a few slug frequencies and one subcooling. That
experiment, not more simulation, is what would turn this from a proposal into a
result.

See **`README_solver.md`** for the in-depth solver documentation.

---

## 7. Author & license

Created, authored and solely maintained by **Akosa Samuel Onyejekwe**.

Released under the [MIT License](LICENSE) — © 2026 Akosa Samuel Onyejekwe.
