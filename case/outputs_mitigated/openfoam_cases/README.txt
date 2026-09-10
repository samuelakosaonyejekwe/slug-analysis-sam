SHCT -> OpenFOAM coupling: CFD cases for the sections that need 3-D resolution.
OpenFOAM detected: False.
Each subfolder is a runnable interFoam (VOF two-phase) case with BCs from the
SHCT 1-D solution. Run each with ./Allrun on an OpenFOAM machine, then call
shct_openfoam.ingest_results(<casedir>) to feed the CFD result back to SHCT.

  - section_1_x1p1km: x=1.14 km — intermittent (slug/churn)
  - section_2_x3p0km: x=2.97 km — intermittent (slug/churn)
  - section_3_x30p4km: x=30.40 km — steep terrain / riser, intermittent (slug/churn), subcooled
      note: flow is reversed (-4.773 m/s); the segment is written in the flow's own direction.
