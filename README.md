# slug-analysis-sam

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
**`report.docx`** (Section 4 — Model equations).

**Published closures:** Bendiksen drift-flux, Taitel–Dukler regime map,
Gregory–Scott / Zabaras slug frequency, Haaland friction, natural-gas hydrate
equilibrium, CSMHyK-type growth, Camargo–Palermo slurry viscosity,
Nielsen–Bucklin inhibitor suppression/sizing, API RP 14E erosional limit,
Peng–Robinson EOS.

**Numerics:** conservative finite-volume with adaptive CFL stepping, optional
2nd-order TVD advection, implicit (tridiagonal/Thomas) pressure, implicit wall-loss
with a Heun corrector, a stochastic Monte-Carlo ensemble for P10/P50/P90 bands, and
a never-fail fallback (auto-degrade to a proven quasi-steady solver on any
non-finite step). Liquid, gas and hydrate mass conserve to ~0 %.

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
├── shct_threed.py                # 3-D field reconstruction + VTK export
├── shct_openfoam.py              # OpenFOAM (interFoam) coupling case generation
├── test_solver.py                # closure + regression test suite
├── README_solver.md              # in-depth solver documentation
├── pyproject.toml / requirements.txt
├── report.docx                    # THE comprehensive report (new medium-crude case): background,
│                                  #   problem, all equations, inputs, every output, validation, calibration
├── validation/data/               # credible published validation datasets (+ recorded sources)
└── case/                          # the deepwater medium-crude-oil case study
    ├── scripts/                   # active pipeline
    │   ├── run_case_study10.py    #   runs the 3 scenarios + advanced stack + validation
    │   ├── build_report.py        #   assembles report.docx from the generated outputs
    │   ├── build_reports.py       #   shared docx helpers + the full equation catalogue
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
cure at ~60 wt% MEG over a ~24 km under-inhibited length; the engineered insulation +
MEG fix removes the subcooling and zeroes the plug probability — the model quantifies
both the threat and the cure.

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

**`report.docx`** (repo root) assembles all of these into a single comprehensive
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
python3 solver.py --calibrate t.json   # validation: fit free constants to measured data

pytest test_solver.py                  # run the test suite
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
- The Φ_SH coupling law and the consolidation/plug mechanism are physically reasoned
  and mass-consistent, but their quantitative form still warrants experimental
  (flow-loop) confirmation.

See **`README_solver.md`** for the in-depth solver documentation.

---

## 7. Author & license

Created, authored and solely maintained by **Akosa Samuel Onyejekwe**.

Released under the [MIT License](LICENSE) — © 2026 Akosa Samuel Onyejekwe.
