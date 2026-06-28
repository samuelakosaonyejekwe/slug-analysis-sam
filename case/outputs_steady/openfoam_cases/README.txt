SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.
OpenFOAM detected: False.
Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the
SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call
shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.

  - section_1_x27p2km: x=27.20 km — Phi_SH>1 (hydrate-critical), intermittent (slug/churn), subcooled, wall deposit
  - section_2_x29p0km: x=29.03 km — Phi_SH>1 (hydrate-critical), intermittent (slug/churn), subcooled, wall deposit
  - section_3_x31p3km: x=31.31 km — Phi_SH>1 (hydrate-critical), steep terrain / riser, intermittent (slug/churn), subcooled, wall deposit
