#!/usr/bin/env python3
"""Prepend the SUPERSEDED note to each older Zenodo record's description.

A published Zenodo record cannot simply be edited: the deposition must be reopened
(actions/edit), its metadata replaced, and the record published again. The DOI, the
version and the files are untouched — only the description gains a leading notice —
so every existing citation keeps working, which is the whole reason for marking
these superseded rather than withdrawing them.

The token is read from a FILE, never from the command line, so it does not reach
the shell history or a process listing:

    printf %s 'YOUR_TOKEN' > ~/.zenodo_token && chmod 600 ~/.zenodo_token
    python3 apply_zenodo_notes.py            # add --dry-run to preview

Create the token at https://zenodo.org/account/settings/applications/tokens/new
with the scopes deposit:write and deposit:actions.
"""
import json
import os
import sys
import urllib.error
import urllib.request

NOTES = "/mnt/c/Users/user/Desktop/paperinfo-slugs_hydrates"
TOKEN_FILE = os.path.expanduser("~/.zenodo_token")
API = "https://zenodo.org/api/deposit/depositions"

#  record id -> the note file written for it
RECORDS = [
    ("22259745", "zenodo_note_v3.1.0_22259745.html", "v3.1.0"),
    ("22311139", "zenodo_note_v3.2.0_22311139.html", "v3.2.0"),
    ("22311939", "zenodo_note_v3.2.1_22311939.html", "v3.2.1"),
    ("22323776", "zenodo_note_v3.3.0_22323776.html", "v3.3.0"),
]
MARKER = "SUPERSEDED"


def call(method, url, token, payload=None):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    if data:
        req.add_header("Content-Type", "application/json")
    with urllib.request.urlopen(req, timeout=60) as r:
        body = r.read().decode()
        return json.loads(body) if body.strip() else {}


def main(argv):
    dry = "--dry-run" in argv
    if not os.path.exists(TOKEN_FILE):
        print(f"  no token at {TOKEN_FILE} — see the module docstring")
        return 2
    token = open(TOKEN_FILE).read().strip()

    for rec, note_file, ver in RECORDS:
        path = os.path.join(NOTES, note_file)
        if not os.path.exists(path):
            print(f"  {ver}: note file missing ({note_file})")
            continue
        note = open(path, encoding="utf-8").read().strip()
        try:
            dep = call("GET", f"{API}/{rec}", token)
        except urllib.error.HTTPError as e:
            print(f"  {ver}: cannot read record {rec} — HTTP {e.code}")
            continue
        md = dep.get("metadata", {})
        desc = md.get("description", "")
        if MARKER in desc[:4000]:
            print(f"  {ver}: already marked superseded, left alone")
            continue
        md["description"] = (f'<p><strong>{MARKER}.</strong></p>{note}'
                             f'<hr/>{desc}')
        if dry:
            print(f"  {ver}: would prepend {len(note)} chars to record {rec}")
            continue
        try:
            edit = call("POST", f"{API}/{rec}/actions/edit", token)
            call("PUT", f"{API}/{edit.get('id', rec)}", token, {"metadata": md})
            call("POST", f"{API}/{edit.get('id', rec)}/actions/publish", token)
            print(f"  {ver}: superseded note applied to record {rec}")
        except urllib.error.HTTPError as e:
            print(f"  {ver}: FAILED on record {rec} — HTTP {e.code} {e.reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
