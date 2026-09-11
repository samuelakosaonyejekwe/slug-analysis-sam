SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.
OpenFOAM detected: False.
Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the
SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call
shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.

  - section_1_x1p1km: x=1.14 km — Phi_SH>1 (hydrate-critical), subcooled, wall deposit
      note: flow is reversed (-0.918 m/s); the segment is written in the flow's own direction.
  - section_2_x12p6km: x=12.57 km — Phi_SH>1 (hydrate-critical), subcooled, wall deposit
      note: flow is reversed (-0.042 m/s); the segment is written in the flow's own direction.
  - section_3_x30p4km: x=30.40 km — Phi_SH>1 (hydrate-critical), steep terrain / riser, intermittent (slug/churn), subcooled, wall deposit
      note: flow is reversed (-1.174 m/s); the segment is written in the flow's own direction.
