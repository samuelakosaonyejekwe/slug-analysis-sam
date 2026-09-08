#!/usr/bin/env python3
"""Mirror the published Zenodo record's metadata into .zenodo.json.

The record is the thing people read, and it can be edited by hand on zenodo.org
(Edit -> change -> Publish, which does NOT make a new version). This pulls that
record back down so the repo's deposit metadata never drifts from it.

    python3 case/scripts/sync_zenodo_metadata.py --check   # report drift only
    python3 case/scripts/sync_zenodo_metadata.py           # write .zenodo.json

Editor cruft is stripped from the description — Word's class="MsoNormal" and the
bare <span> wrappers a paste leaves behind — while <p>, <strong>, <em>, <sub> and
<sup> are kept, so the text and its formatting mirror the record exactly.
"""
import argparse
import difflib
import html
import json
import os
import re
import urllib.request

CONCEPT = "22259744"                      # concept DOI 10.5281/zenodo.22259744
RECORD = "22666728"                       # any version of it; the latest is resolved from here
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LOCAL = os.path.join(ROOT, ".zenodo.json")
API = "https://zenodo.org/api/records/"


class Unreachable(Exception):
    """Zenodo could not be read. Not drift, and not a failure of this repository."""


def get(url):
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return json.loads(r.read().decode())
    except Exception as exc:                       # HTTPError, URLError, timeout, bad JSON
        raise Unreachable(f"{type(exc).__name__}: {exc}") from exc


def fetch(record):
    """Return the NEWEST version of the record's family, not the pinned one.

    A record id names one version. Following links.latest means that when a new
    version is published the sync tracks it instead of quietly mirroring an old
    release forever.
    """
    rec = get(API + record)
    latest = (rec.get("links") or {}).get("latest")
    if latest:
        try:
            rec = get(latest)
        except Unreachable as exc:                    # network hiccup: use what we have
            print(f"note: could not resolve latest version ({exc}); using {record}")
    return rec


def _ver(v):
    """A dotted version as a tuple of ints; None when it is not one."""
    try:
        parts = [int(x) for x in str(v).strip().split(".")]
    except (TypeError, ValueError):
        return None
    return tuple(parts) if parts else None


def _ahead(local_v, live_v):
    """True when the local version is strictly newer than the published one."""
    a, b = _ver(local_v), _ver(live_v)
    return bool(a and b and a > b)


def clean(desc):
    desc = re.sub(r'\s+class="[^"]*"', "", desc)      # MsoNormal and friends
    desc = re.sub(r"</?span[^>]*>", "", desc)          # bare paste wrappers
    desc = re.sub(r"<p>\s*</p>", "", desc)
    return re.sub(r"\s+", " ", desc).replace("> <", "><").strip()


def plain(s):
    return re.sub(r"\s+", " ", html.unescape(re.sub("<[^>]+>", " ", s))).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="report drift, write nothing")
    ap.add_argument("--record", default=RECORD)
    a = ap.parse_args()

    #  A scheduled job that polls a third party must not go red when that third party
    #  is having an outage. Zenodo returned 504 for about two hours on 2026-09-08 and
    #  this script crashed with an unhandled HTTPError on every five-minute run: 4 of the
    #  last 29 runs failed that way. Unreachable is not drift -- there is simply nothing
    #  to compare against -- so it reports and exits 0, and a real difference still exits
    #  1 under --check.
    try:
        rec = fetch(a.record)
    except Unreachable as exc:
        print(f"zenodo is not reachable right now ({exc}); nothing to compare, "
              f"leaving {os.path.basename(LOCAL)} untouched")
        return 0
    md = rec.get("metadata", {})
    with open(LOCAL, encoding="utf-8") as fh:
        local = json.load(fh)
    live = {
        "description": clean(md.get("description", "")),
        "title": md.get("title", local.get("title")),
        "version": md.get("version", local.get("version")),
        "keywords": md.get("keywords", local.get("keywords")),
    }

    #  A VERSION AHEAD OF THE RECORD IS A PENDING RELEASE, NOT DRIFT. This script is run
    #  every five minutes by a workflow that commits what it writes, so editing
    #  .zenodo.json for a release you have not deposited yet gets your edit reverted
    #  within five minutes. That happened on 2026-09-08: the file was bumped to 4.0.0,
    #  the workflow pulled the still-published 3.4.0 back over it, and the bump had to be
    #  restored by hand. The deposit itself was unaffected -- Zenodo archives the tag, not
    #  the tip of main -- but the working tree fought the author.
    #
    #  The record stays authoritative in the direction that matters: anything edited on
    #  zenodo.org is still pulled down. Only a local version STRICTLY AHEAD of the
    #  published one is left alone, and it is announced rather than passed over.
    if _ahead(local.get("version"), live.get("version")):
        print(f"local version {local.get('version')} is ahead of the published "
              f"{live.get('version')}: treating this as a release not yet deposited and "
              f"leaving {os.path.basename(LOCAL)} alone. Re-run after the deposit.")
        return 0

    drift = {k: v for k, v in live.items() if local.get(k) != v}
    text_drift = plain(live["description"]) != plain(local.get("description", ""))
    print(f"record {rec.get('id', a.record)} (latest of concept {CONCEPT}): "
          f"version {live['version']}, doi {rec.get('doi')}")
    print("fields differing:", ", ".join(drift) or "none")
    print("description TEXT differs:", text_drift)
    if text_drift:
        d = difflib.unified_diff(plain(local.get("description", "")).split(". "),
                                 plain(live["description"]).split(". "),
                                 "repo", "zenodo", lineterm="", n=0)
        print("\n".join(list(d)[:40]))
    if not drift:
        print("in sync — nothing to write")
        return 0
    if a.check:
        print("--check: not writing")
        return 1
    local.update(drift)
    with open(LOCAL, "w", encoding="utf-8") as fh:
        json.dump(local, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    print("wrote", LOCAL)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
