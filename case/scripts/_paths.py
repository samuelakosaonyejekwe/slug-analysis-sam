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

#  ---------------------------------------------------------------------------
#  EXTERNAL DOCUMENTS (manuscript, thesis, slide deck) — NOT part of this repository.
#
#  Fifteen scripts in this folder carried an absolute path into one machine's Windows
#  Desktop as their DEFAULT target. On any other checkout those defaults name nothing,
#  so the scripts either silently did nothing or -- for check_docs.py and
#  check_journal_artwork.py -- reported failures about documents this repository does
#  not contain and cannot fix, drowning the one document it does build.
#
#  The documents are opt-in now. Point SHCT_DOCS_DIR at a folder to work on them:
#
#      SHCT_DOCS_DIR=/path/to/manuscripts python3 case/scripts/check_docs.py
#
#  or pass the files as arguments. With neither, each script confines itself to what is
#  in this repository and says so.
#  ---------------------------------------------------------------------------
def docs_dir():
    """The external manuscript/deck folder, or None when none is configured."""
    d = os.environ.get("SHCT_DOCS_DIR", "").strip()
    return d if d and os.path.isdir(d) else None


#  When SHCT_DOCS_DIR is unset, docs_path() returns a path under this sentinel rather
#  than None: every caller already tests it with os.path.exists / os.path.isdir, which
#  is False for a path that cannot exist, so each script takes its own "not found"
#  branch and prints a message that names the variable to set. Returning None instead
#  made all twelve of them die in os.path.exists with a TypeError.
DOCS_UNSET = "«set SHCT_DOCS_DIR to the folder holding these documents»"


def docs_path(*parts):
    """A path inside the external document folder.

    When no folder is configured this is a path that deliberately does not exist, so a
    caller's existence check fails cleanly and its own not-found message is printed.
    """
    return os.path.join(docs_dir() or DOCS_UNSET, *parts)


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
