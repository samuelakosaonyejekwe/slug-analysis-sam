#!/usr/bin/env python3
"""Export .docx -> .pdf through a NEW Word instance, leaving the user's alone."""
import os, shutil, subprocess, sys

def wintemp():
    w = subprocess.check_output(["cmd.exe","/c","echo %TEMP%"],stderr=subprocess.DEVNULL).decode().strip()
    return subprocess.check_output(["wslpath","-u",w]).decode().strip()

def export(src):
    src = os.path.abspath(src)
    dst = os.path.splitext(src)[0] + ".pdf"
    lt = wintemp()
    stem = f"_fx_{os.getpid()}_{os.path.basename(src)}"
    tdocx = os.path.join(lt, stem)
    tpdf  = os.path.splitext(tdocx)[0] + ".pdf"
    for p in (tdocx, tpdf):
        if os.path.exists(p): os.remove(p)
    shutil.copyfile(src, tdocx)
    W = lambda p: subprocess.check_output(["wslpath","-w",p]).decode().strip()
    ps1 = os.path.join(lt, f"_fx_{os.getpid()}.ps1")
    open(ps1,"w",encoding="utf-8").write(f"""
$ErrorActionPreference='Stop'
$w = New-Object -ComObject Word.Application
$w.Visible = $false
$w.DisplayAlerts = 0
try {{
  $d = $w.Documents.Open('{W(tdocx)}', $false, $false)
  foreach ($sr in $d.StoryRanges) {{ $null = $sr.Fields.Update() }}
  $d.ExportAsFixedFormat('{W(tpdf)}', 17)
  $n = $d.ComputeStatistics(2)
  $d.Close(0)
  Write-Output ("PAGES=" + $n)
}} finally {{ $w.Quit() }}
""")
    r = subprocess.run(["powershell.exe","-NoProfile","-ExecutionPolicy","Bypass","-File", W(ps1)],
                       capture_output=True, text=True, timeout=1500)
    print(r.stdout.strip()[:200], r.stderr.strip()[:300])
    if not os.path.exists(tpdf):
        print("FAIL", src); return False
    shutil.copyfile(tpdf, dst)
    print(f"OK  {os.path.basename(src)} -> {dst}  ({os.path.getsize(dst)/1e6:.2f} MB)")
    for p in (tdocx, tpdf, ps1):
        try: os.remove(p)
        except OSError: pass
    return True

for a in sys.argv[1:]:
    export(a)
