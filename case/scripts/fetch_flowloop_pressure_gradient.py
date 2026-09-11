#!/usr/bin/env python3
"""Retrieve the MEASURED PRESSURE GRADIENT for the das Neves flow-loop points.

Why this matters. The one gap this project has never closed is a dP profile from a real
facility: Roberts gives field OUTCOMES but no series, and the das Neves pressure channels
that ship in the time-series files are PIEZOELECTRIC and AC-coupled -- they measure
fluctuation, not level, so no gradient can be formed from them. That was recorded as
"needs an operator's dataset".

It may not. The paper's Table 2 lists "Pressure gradient [Pa/m]" as a per-point METADATA
variable, and states gradients exist for points 8-19, 21 and 22 -- exactly the points this
project already scores holdup against. The values are not in the paper; they are in the
per-point metadata .pkl files in the dataset repository. If they can be read, the solver's
pressure-gradient closure can be scored against real measurements on the same points, at
the same conditions, as the holdup score.

Access route, established by probing: the Dataverse record (doi:10.25824/redu/ISMWP4)
carries only documentation. The data sits on a Nextcloud public share whose DAV endpoint
refuses PROPFIND (401) but whose /download endpoint SERVES files when carrying a session
cookie from the share page -- it 303-redirects to the same DAV path and succeeds. So
listing is impossible and retrieval is not: paths must be guessed and confirmed by a
GET that returns content rather than an error page.
"""
import io
import os
import pickle
import sys
import urllib.parse

import requests

TOKEN = "LSGjo8bnfnT6cEX"
BASE = "https://unidrive.unicamp.br"
SHARE = f"{BASE}/index.php/s/{TOKEN}"
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/126.0"}
POINTS = [8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 21, 22]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                   "validation", "data", "flowloop_pressure_gradient_dasneves2025.json")


def session():
    s = requests.Session()
    s.headers.update(UA)
    s.get(SHARE, timeout=60)          # sets the share cookie
    return s


def fetch(s, path, name, timeout=90):
    """GET one file through the share's download endpoint. Returns bytes or None."""
    q = urllib.parse.urlencode({"path": path, "files": name})
    try:
        r = s.get(f"{SHARE}/download?{q}", timeout=timeout, allow_redirects=True)
    except Exception:
        return None
    if r.status_code != 200 or not r.content:
        return None
    #  an HTML error page is not a file
    head = r.content[:200].lstrip().lower()
    if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
        return None
    return r.content


def try_layouts(s, pt):
    """The folder/file naming is not documented, so try the plausible conventions."""
    roots = ["/Dataset_800Hz", "/Dataset_25_6kHz"]
    folders = ["Point_%d" % pt, "point_%d" % pt, "P%d" % pt, "%d" % pt,
               "Ponto_%d" % pt, "Point %d" % pt]
    files = ["metadata.pkl", "meta.pkl", "metadata_%d.pkl" % pt, "Metadata.pkl",
             "metadados.pkl", "metadata_point_%d.pkl" % pt]
    for root in roots:
        for fol in folders:
            for fn in files:
                blob = fetch(s, "%s/%s" % (root, fol), fn)
                if blob:
                    return "%s/%s/%s" % (root, fol, fn), blob
    return None, None


def main():
    s = session()
    got, missing = {}, []
    for pt in POINTS:
        where, blob = try_layouts(s, pt)
        if blob is None:
            missing.append(pt)
            print("  point %-3d not reached" % pt, flush=True)
            continue
        try:
            meta = pickle.load(io.BytesIO(blob))
        except Exception as e:
            print("  point %-3d fetched %s but could not unpickle (%s)"
                  % (pt, where, type(e).__name__), flush=True)
            missing.append(pt)
            continue
        got[pt] = {"source_path": where,
                   "keys": sorted(map(str, meta.keys())) if hasattr(meta, "keys") else None,
                   "metadata": {str(k): (v if isinstance(v, (int, float, str, bool)) else str(v))
                                for k, v in meta.items()} if hasattr(meta, "items") else str(meta)}
        print("  point %-3d OK  <- %s" % (pt, where), flush=True)

    print()
    if not got:
        print("NOTHING RETRIEVED. The share serves files but cannot be listed, so this")
        print("depends on guessing the folder layout and none of the tried conventions hit.")
        print("The gap stays open and stays honestly described: no dP profile from a real")
        print("facility is in this repository.")
        return 1
    import json
    with open(OUT, "w") as fh:
        json.dump({"source": "das Neves et al. (2025) Data in Brief, dataset "
                             "doi:10.25824/redu/ISMWP4, Unicamp research data repository",
                   "what": "per-point metadata including mean absolute pressure and pressure "
                           "gradient [Pa/m] -- the quantity the shipped time-series files cannot "
                           "give, because their pressure channels are AC-coupled piezoelectric",
                   "points_retrieved": sorted(got), "points_missing": missing,
                   "records": got}, fh, indent=2)
        fh.write("\n")
    print("written %s  (%d of %d points)" % (OUT, len(got), len(POINTS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
