#!/usr/bin/env python3
"""Convert .docx -> .pdf via Word COM from WSL, WITHOUT quitting Word.

The converter in build_report.py calls $word.Quit(), which will close a Word
session the user already has open (COM attaches to the running instance) and can
take unsaved work with it. This variant opens the document read-only and
invisible, exports it, closes only that document, and leaves the application
running. Safe to use while Word is open.

    python3 docx2pdf_safe.py <file.docx> [<file2.docx> ...]

Each output is written next to its source as <name>.pdf.
"""
import os
import subprocess
import sys
import shutil


WRITE_BACK = False


def _win_temp():
    wtemp = subprocess.check_output(["cmd.exe", "/c", "echo %TEMP%"],
                                    stderr=subprocess.DEVNULL).decode().strip()
    ltemp = subprocess.check_output(["wslpath", "-u", wtemp]).decode().strip()
    if not os.path.isdir(ltemp):
        raise RuntimeError(f"Windows TEMP not reachable: {ltemp}")
    return ltemp


def convert(docx_path, pdf_path=None, timeout=900, write_back=False):
    docx_path = os.path.abspath(docx_path)
    if pdf_path is None:
        pdf_path = os.path.splitext(docx_path)[0] + ".pdf"
    if not shutil.which("powershell.exe"):
        raise RuntimeError("powershell.exe not on PATH — not running under WSL?")

    ltemp = _win_temp()
    stem = "_conv_" + os.path.splitext(os.path.basename(docx_path))[0]
    tmp_docx = os.path.join(ltemp, stem + ".docx")
    tmp_pdf = os.path.join(ltemp, stem + ".pdf")
    for f in (tmp_docx, tmp_pdf):
        if os.path.exists(f):
            os.remove(f)
    shutil.copyfile(docx_path, tmp_docx)

    win_docx = subprocess.check_output(["wslpath", "-w", tmp_docx]).decode().strip()
    win_pdf = subprocess.check_output(["wslpath", "-w", tmp_pdf]).decode().strip()

    #  NOTE: no $w.Quit(). If Word was already running we are attached to the
    #  user's live instance; quitting it would discard their unsaved work.
    #
    #  Fields are updated before export. Caption numbers are SEQ fields with no
    #  cached result, so without this every caption exports as "Figure ." with the
    #  number missing; the same applies to the TOC and to cross-references. The
    #  document is opened writable (not read-only) so the fields can be refreshed,
    #  and saved back only when write_back is set.
    save_clause = ("$d.Save();" if write_back else "")
    readonly = "$false" if write_back else "$true"
    ps = (
        "$ErrorActionPreference='Stop';"
        "try { $w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application'); "
        "$preexisting = $true } "
        "catch { $w = New-Object -ComObject Word.Application; $preexisting = $false }"
        "$w.DisplayAlerts = 0;"
        f"$d = $w.Documents.Open([ref]'{win_docx}', [ref]$false, [ref]{readonly}, [ref]$false,"
        " [ref]'', [ref]'', [ref]$true, [ref]'', [ref]'', [ref]0, [ref]0, [ref]$false);"
        #  refresh every story (body, headers, footers, textboxes) then the TOCs
        "foreach ($sr in $d.StoryRanges) { $null = $sr.Fields.Update() };"
        "foreach ($t in $d.TablesOfContents) { $t.Update() };"
        "foreach ($t in $d.TablesOfFigures) { $t.Update() };"
        f"{save_clause}"
        f"$d.ExportAsFixedFormat([ref]'{win_pdf}', [ref]17);"
        "$pages = $d.ComputeStatistics(2);"
        "$d.Close([ref]$false);"
        "if (-not $preexisting) { $w.Quit() }"
        "Write-Output ('PAGES=' + $pages)"
    )
    res = subprocess.run(["powershell.exe", "-NoProfile", "-Command", ps],
                         capture_output=True, text=True, timeout=timeout)
    pages = None
    for line in (res.stdout or "").splitlines():
        if line.startswith("PAGES="):
            pages = line.split("=", 1)[1].strip()
    if not os.path.exists(tmp_pdf):
        raise RuntimeError(f"conversion produced no PDF.\nstdout: {res.stdout}\n"
                           f"stderr: {res.stderr}")
    shutil.copyfile(tmp_pdf, pdf_path)
    if write_back:
        shutil.copyfile(tmp_docx, docx_path)   # keeps the refreshed field results
    for f in (tmp_docx, tmp_pdf):
        try:
            os.remove(f)
        except OSError:
            pass
    return pdf_path, pages


def main(argv):
    global WRITE_BACK
    if argv and argv[0] == "--write-back":
        WRITE_BACK = True
        argv = argv[1:]
    if not argv:
        print(__doc__)
        return 2
    rc = 0
    for src in argv:
        try:
            out, pages = convert(src, write_back=WRITE_BACK)
            size = os.path.getsize(out)
            print(f"OK   {os.path.basename(src)} -> {out}  "
                  f"({size/1e6:.2f} MB{', ' + pages + ' pages' if pages else ''})")
        except Exception as e:
            print(f"FAIL {src}: {e}")
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
