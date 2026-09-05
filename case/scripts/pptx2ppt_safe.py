#!/usr/bin/env python3
"""Convert .pptx -> legacy .ppt via PowerPoint COM from WSL, WITHOUT quitting it.

The companion of docx2pdf_safe.py, and it exists for the same reason: COM
attaches to whatever PowerPoint instance is already running, so a naive
$p.Quit() would close the user's live session and take unsaved work with it.
This opens the deck invisibly, saves a copy in the 97-2003 format, closes only
that presentation, and leaves the application as it found it.

    python3 pptx2ppt_safe.py <deck.pptx> [<deck2.pptx> ...]

Each output is written next to its source as <name>.ppt. The .ppt twin exists so
the deck opens on machines without a modern PowerPoint; it must be regenerated
whenever the .pptx changes, or the two drift apart silently — which is exactly
what happened once already, leaving a .ppt stamped three hours before its .pptx.
"""
import os
import shutil
import subprocess
import sys

PPSAVEASPRESENTATION = 1          # PpSaveAsFileType: 97-2003 .ppt


def _win_temp():
    wtemp = subprocess.check_output(["cmd.exe", "/c", "echo %TEMP%"],
                                    stderr=subprocess.DEVNULL).decode().strip()
    ltemp = subprocess.check_output(["wslpath", "-u", wtemp]).decode().strip()
    if not os.path.isdir(ltemp):
        raise RuntimeError(f"Windows TEMP not reachable: {ltemp}")
    return ltemp


def convert(pptx_path, ppt_path=None, timeout=1800):
    pptx_path = os.path.abspath(pptx_path)
    if ppt_path is None:
        ppt_path = os.path.splitext(pptx_path)[0] + ".ppt"
    if not shutil.which("powershell.exe"):
        raise RuntimeError("powershell.exe not on PATH — not running under WSL?")

    ltemp = _win_temp()
    stem = f"_conv_{os.getpid()}_" + os.path.splitext(os.path.basename(pptx_path))[0]
    tmp_pptx = os.path.join(ltemp, stem + ".pptx")
    tmp_ppt = os.path.join(ltemp, stem + ".ppt")
    for f in (tmp_pptx, tmp_ppt):
        if os.path.exists(f):
            os.remove(f)
    shutil.copyfile(pptx_path, tmp_pptx)

    win_pptx = subprocess.check_output(["wslpath", "-w", tmp_pptx]).decode().strip()
    win_ppt = subprocess.check_output(["wslpath", "-w", tmp_ppt]).decode().strip()

    #  PowerPoint has no invisible Open (WithWindow:=0 is unsupported on some
    #  builds), so the window is opened minimised rather than hidden, and the
    #  application is only quit if it was not already running.
    ps = (
        "$ErrorActionPreference='Stop';"
        "try { $p = [Runtime.InteropServices.Marshal]::GetActiveObject('PowerPoint.Application'); "
        "$preexisting = $true } "
        "catch { $p = New-Object -ComObject PowerPoint.Application; $preexisting = $false }"
        f"$d = $p.Presentations.Open('{win_pptx}', $true, $false, $false);"
        f"$d.SaveAs('{win_ppt}', {PPSAVEASPRESENTATION});"
        "$slides = $d.Slides.Count;"
        "$d.Close();"
        "if (-not $preexisting) { $p.Quit() }"
        "Write-Output ('SLIDES=' + $slides)"
    )
    res = subprocess.run(["powershell.exe", "-NoProfile", "-Command", ps],
                         capture_output=True, text=True, timeout=timeout)
    slides = None
    for line in (res.stdout or "").splitlines():
        if line.startswith("SLIDES="):
            slides = line.split("=", 1)[1].strip()
    if not os.path.exists(tmp_ppt):
        raise RuntimeError(f"conversion produced no .ppt.\nstdout: {res.stdout}\n"
                           f"stderr: {res.stderr}")
    shutil.copyfile(tmp_ppt, ppt_path)
    for f in (tmp_pptx, tmp_ppt):
        try:
            os.remove(f)
        except OSError:
            pass
    return ppt_path, slides


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    rc = 0
    for src in argv:
        try:
            out, slides = convert(src)
            print(f"OK   {os.path.basename(src)} -> {out}  "
                  f"({os.path.getsize(out)/1e6:.2f} MB"
                  f"{', ' + slides + ' slides' if slides else ''})")
        except Exception as e:
            print(f"FAIL {src}: {e}")
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
