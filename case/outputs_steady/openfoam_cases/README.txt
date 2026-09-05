SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.
OpenFOAM detected: False.
Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the
SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call
shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.

  - section_1_x16p2km: x=16.23 km — Phi_SH>1 (hydrate-critical), intermittent (slug/churn), subcooled, wall deposit
  - section_2_x27p7km: x=27.66 km — Phi_SH>1 (hydrate-critical), intermittent (slug/churn), subcooled, wall deposit
  - section_3_x30p4km: x=30.40 km — Phi_SH>1 (hydrate-critical), steep terrain / riser, intermittent (slug/churn), subcooled, wall deposit
