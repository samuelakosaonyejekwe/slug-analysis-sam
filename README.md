# slug-analysis-sam

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22259744.svg)](https://doi.org/10.5281/zenodo.22259744)

Archived on Zenodo. The badge above is the *concept* DOI and always resolves to the
newest release; the version archived for the journal submission is
[10.5281/zenodo.22311139](https://doi.org/10.5281/zenodo.22311139) (SHCT v3.2.0).

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
├── shct_openfoam.py              # OpenFOAM (interFoam) coupling case generation
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
~30° API medium crude oil (C1 ≈ 43 mol%, ~31 mol% C7+ tail) at 35 % water cut over a
cold (4 °C), undulating seabed. This geometry and fluid is a textbook combination for
**both** slugging and hydrates, so it exercises the whole prediction chain.

Three scenarios are run end-to-end through the real solver:

| Scenario | Folder | Description |
|----------|--------|-------------|
| **A — as-operated** | `case/outputs_steady/` | normal production, degraded (water-flooded) insulation, no inhibitor → the high-risk prediction |
| **B — shut-in** | `case/outputs_shutin/` | unplanned shut-in cooldown → no-touch time |
| **C — mitigated** | `case/outputs_mitigated/` | restored multi-layer insulation + continuous MEG → risk removed (design tool) |

**Headline result (as-operated):** intermittent flow over the whole line with slugs
up to ~37 m; the cold under-insulated wall drives the fluid ~21 °C into the hydrate
region (Φ_SH ≫ 1 over the cold section), giving a 100 % plug probability with a P50
time-to-plug of only ~2.8 h and a peak wall deposit of ~117 mm. The model sizes the
remedy at ~60 wt% MEG over a ~24 km under-inhibited length, and the engineered insulation +
MEG fix removes the subcooling and zeroes the plug probability.

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
> as-operated case 5.92 % of the injected liquid (≈ 563 m³) has nowhere to go and is
> dropped at the bounds. It was previously invisible.
> (iii) The slug-length statistics averaged in the correlation's 5000 m "not slugging"
> ceiling. Slug lengths, velocities and the erosional check should be taken from this
> release rather than the last; the hydrate and coupling results are unchanged in
> character.

> **The two-fluid description stays well posed.** `27_wellposedness_map.png` reports the
> slip against the inviscid Kelvin–Helmholtz limit at which the one-dimensional two-fluid
> model loses hyperbolicity. For this case the margin peaks at ~0.86 and never reaches 1,
> so the predicted slug activity is a property of the flow rather than a grid-dependent
> artefact of an ill-posed initial-value problem.

> **Read these magnitudes with care.** A ~60 wt% MEG requirement sits well outside normal
> field practice (typical continuous doses are 20–50 wt%), and the peak Φ_SH, the 2.8 h P50
> and the ~0 h no-touch time are all at or beyond the edge of reported field experience.
> They are the model's answer for a deliberately extreme, water-flooded-insulation,
> uninhibited scenario, produced with literature-typical kinetic constants that have been
> fitted to no dataset. `--sensitivity` (see §5) reports how far each of them moves when
> C, n and k_g0 are swept across their plausible ranges. Treat them as model outputs, not
> as calibrated predictions.

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
                                       #   and deposit to the ASSUMED constants kg0, n and C_phi
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
| Hydrate equilibrium curve | **validated** | Deaton & Frost (1946) measurements; 1.72 °C RMSE |
| Mass conservation | **verified** | liquid ~1e-11 %, gas ~1e-14 % |
| Two-fluid well-posedness | **verified** | inviscid Kelvin–Helmholtz limit; margin ≤ 0.86 |
| Φ_SH dimensional consistency | **verified** | dimensionless for any *n*; invariance test in the suite |
| **Φ_SH magnitude and the Φ_SH = 1 criterion** | **NOT validated** | *no dataset* |
| **C, n, k_g0** | **NOT fitted** | literature-typical values only |
| **Whole-system prediction vs a reference simulator** | **NOT yet run** | see below |

The central proposal of this work — that the competition between hydrate
deposition and slug renewal is captured by a single dimensionless group, and that
Φ_SH = 1 separates scoured from plugging-critical — rests on physical reasoning,
dimensional consistency and internal consistency. **It has not been tested against
experiment.** Treat it as a hypothesis with a solver behind it, not a calibrated
predictor, and read §"Read these magnitudes with care" above alongside it.

### Benchmarking against a reference simulator

`shct_benchmark.py` runs SHCT against OLGA, LedaFlow or any transient multiphase
code on an identical case and reports MAE, RMSE, normalised RMSE and the worst
deviation for holdup, pressure, temperature and mixture velocity, with a
comparison figure.

**No reference dataset ships with this repository.** There is no licence for such
a tool in the development environment, and fabricated benchmark numbers would be
worse than none — a made-up agreement is indistinguishable from a real one until
somebody tries to reproduce it. Export your own reference run in the schema
documented at the top of `shct_benchmark.py` (the geometry, fluid and boundary
conditions are all in `case/outputs_*/input_data_deck.csv` and
`feed_composition.csv`), then:

```bash
python3 shct_benchmark.py validation/data/olga_asoperated.json
```

The loader refuses a file that does not name the tool that produced it, so a
benchmark in this repository always carries its provenance.

### The measurement that would settle it

Φ_SH is falsifiable, and cheaply. In a flow loop at fixed subcooling, the model
predicts the deposit growth rate to fall as the slug frequency rises, crossing
from net accumulation to net removal at Φ_SH = 1. Sweeping f_slug at constant
ΔT_sub while measuring wall deposit over time would confirm or refute both the
form of the group and the location of the threshold. That experiment, not more
simulation, is what would turn this from a proposal into a result.

See **`README_solver.md`** for the in-depth solver documentation.

---

## 7. Author & license

Created, authored and solely maintained by **Akosa Samuel Onyejekwe**.

Released under the [MIT License](LICENSE) — © 2026 Akosa Samuel Onyejekwe.
