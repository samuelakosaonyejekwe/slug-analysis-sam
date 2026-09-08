#!/usr/bin/env python3
# =============================================================================
#  _paths.py  —  one shared, layout-independent path map for every build/run
#  script of the crude-oil case study, plus installation of the global
#  no-black / no-dark plotting style (shct_style).
#
#      from _paths import ROOT, CASE, HERE, OUTROOT, OUT
#
#  Those five are what this module defines. It used to advertise FIGURES,
#  SIMU_PLOTS, REPORT_PLOTS and REPORTS as well, naming three directories that do
#  not exist in this repository -- the documented import line raised ImportError,
#  and only the fact that no caller ever used those names kept it from being noticed.
#
#  Layout:
#      <ROOT>/                  repo root (solver.py, shct_*.py, shct_style.py)
#      <ROOT>/case/             CASE
#      <ROOT>/case/scripts/     this folder (build/run scripts)  = HERE
#      <ROOT>/case/outputs_*/   solver outputs (OUT: steady, shutin, mitigated)
# =============================================================================
import os
import sys

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

import shct_style  # noqa: E402  — installs the no-black / no-dark rcParams

shct_style.apply_style()
