#!/usr/bin/env python3
# =============================================================================
#  _paths.py  —  one shared, layout-independent path map for every build/run
#  script of the crude-oil case study, plus installation of the global
#  no-black / no-dark plotting style (shct_style).
#
#      from _paths import ROOT, CASE, OUTROOT, OUT, FIGURES, SIMU_PLOTS, \
#                         REPORT_PLOTS, REPORTS
#
#  Layout:
#      <ROOT>/                  repo root (solver.py, shct_*.py, shct_style.py)
#      <ROOT>/case/             CASE
#      <ROOT>/case/scripts/     this folder (build/run scripts)  = HERE
#      <ROOT>/case/outputs_*/   solver outputs (OUT)
#      <ROOT>/case/figures/     dissertation figures (FIGURES)
#      <ROOT>/case/simu_plots/  simulation-output curves (SIMU_PLOTS)
#      <ROOT>/case/reports/     .docx report sources (REPORTS)
# =============================================================================
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))     # case/scripts
CASE = os.path.dirname(HERE)                           # case
ROOT = os.path.dirname(CASE)                           # repo root

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)                           # import solver / shct_* / shct_style

OUTROOT = CASE
OUT = {
    "steady":    os.path.join(CASE, "outputs_steady"),
    "shutin":    os.path.join(CASE, "outputs_shutin"),
    "mitigated": os.path.join(CASE, "outputs_mitigated"),
}
# only the three scenario output directories the ACTIVE new-case pipeline writes are
# ensured here; build_report.py keeps its own intermediate plots under case/scripts/.
for _d in OUT.values():
    os.makedirs(_d, exist_ok=True)

import shct_style          # noqa: E402  — installs the no-black / no-dark rcParams
shct_style.apply_style()
