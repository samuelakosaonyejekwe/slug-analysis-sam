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
import argparse, difflib, html, json, os, re, sys, urllib.request

CONCEPT = "22259744"                      # concept DOI 10.5281/zenodo.22259744
RECORD = "22348213"                       # any version of it; the latest is resolved from here
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LOCAL = os.path.join(ROOT, ".zenodo.json")
API = "https://zenodo.org/api/records/"


def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())


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
        except Exception as exc:                      # network hiccup: use what we have
            print(f"note: could not resolve latest version ({exc}); using {record}")
    return rec


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

    rec = fetch(a.record)
    md = rec.get("metadata", {})
    local = json.load(open(LOCAL, encoding="utf-8"))
    live = {
        "description": clean(md.get("description", "")),
        "title": md.get("title", local.get("title")),
        "version": md.get("version", local.get("version")),
        "keywords": md.get("keywords", local.get("keywords")),
    }

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
