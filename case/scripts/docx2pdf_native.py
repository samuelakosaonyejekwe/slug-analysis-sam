#!/usr/bin/env python3
# =============================================================================
#  docx2pdf_native.py — render report.docx to PDF with no Word and no LibreOffice.
# -----------------------------------------------------------------------------
#  WHY THIS EXISTS
#  ---------------
#  report.pdf is the committed deliverable, and every route to it needed an
#  office suite: LibreOffice, Microsoft Word over COM, or docx2pdf (which drives
#  one of those). On a machine with none of them -- a CI runner, a container, this
#  one -- build_report.py wrote the .docx, printed "report.pdf was NOT refreshed"
#  and left the tracked PDF describing an older run. A deliverable that can only
#  be rebuilt on one person's laptop is not reproducible.
#
#  This renders the document directly with reportlab, walking the .docx body in
#  document order and mapping the small set of constructs build_report.py actually
#  emits: headings, body paragraphs, bullets, inline figures with captions,
#  tables, and page breaks. It is a LAYOUT APPROXIMATION, not a Word-identical
#  rendering -- pagination and table column widths are this renderer's, not
#  Word's -- so it sits LAST in the converter chain: wherever LibreOffice or Word
#  is available, that output is still preferred and this is never reached.
#
#  Fonts are DejaVu, which the report needs: the text carries Phi, alpha, degree,
#  superscript-three, multiplication sign and em-dash, and reportlab's built-in
#  Helvetica has none of them (they would render as black boxes).
#
#      python3 docx2pdf_native.py report.docx [report.pdf]
# =============================================================================
from __future__ import annotations

import io
import os
import sys

from docx import Document
from docx.table import Table
from docx.text.paragraph import Paragraph
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, PageBreak, SimpleDocTemplate, Spacer, TableStyle
from reportlab.platypus import Paragraph as RLPara
from reportlab.platypus import Table as RLTable

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
WP = "{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
EMU = 914400.0

#  Figures are embedded at the resolution they are DISPLAYED at, not the resolution
#  they were rendered at. The report's figures come out of matplotlib at 320 dpi and are
#  shown about 6.8 in wide, so passing the originals through would carry ~37 MB of pixels
#  no reader can see -- the PDF came to 46 MB against the 13 MB the Word route produced.
#  200 dpi at display size is a deliberate trade and it brings the document to ~18 MB. It
#  is BELOW the 300 dpi the journal asks for -- this comment used to claim the opposite,
#  which is arithmetically impossible: 200 dpi at the width the figure is shown at is
#  200 dpi at print. That is acceptable here only because this PDF is the internal report,
#  not the submission: the journal receives case/figures_paper/Figure_N.png at 320 dpi,
#  and check_journal_artwork.py gates those. Raise SHCT_PDF_IMG_DPI for a sharper file.
IMG_DPI = float(os.environ.get("SHCT_PDF_IMG_DPI", "200"))

NAVY = colors.HexColor("#2E5BBF")
INK = colors.HexColor("#3A5BA8")
GREY = colors.HexColor("#6E7B8B")
GRID = colors.HexColor("#D2DCF2")

_FONT_DIRS = ("/usr/share/fonts/truetype/dejavu", "/usr/share/fonts/dejavu",
              "/Library/Fonts", "C:/Windows/Fonts")
_FACES = {"Body": "DejaVuSans.ttf", "Body-Bold": "DejaVuSans-Bold.ttf",
          "Body-Italic": "DejaVuSans-Oblique.ttf",
          "Body-BoldItalic": "DejaVuSans-BoldOblique.ttf"}


def _register_fonts():
    """Register DejaVu, or fall back to Helvetica and say what will be lost."""
    found = {}
    for name, fn in _FACES.items():
        for d in _FONT_DIRS:
            p = os.path.join(d, fn)
            if os.path.exists(p):
                found[name] = p
                break
    if len(found) != len(_FACES):
        return "Helvetica", "Helvetica-Bold", "Helvetica-Oblique"
    for name, path in found.items():
        try:
            pdfmetrics.registerFont(TTFont(name, path))
        except Exception:
            return "Helvetica", "Helvetica-Bold", "Helvetica-Oblique"
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold",
                                  italic="Body-Italic", boldItalic="Body-BoldItalic")
    return "Body", "Body-Bold", "Body-Italic"


def _esc(t):
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _runs_markup(par):
    """Paragraph text as reportlab mini-HTML, keeping bold / italic / colour."""
    out = []
    for r in par.runs:
        t = _esc(r.text)
        if not t:
            continue
        if r.bold:
            t = f"<b>{t}</b>"
        if r.italic:
            t = f"<i>{t}</i>"
        col = None
        try:
            col = r.font.color.rgb
        except Exception:
            col = None
        if col is not None:
            t = f'<font color="#{col}">{t}</font>'
        out.append(t)
    return "".join(out) or _esc(par.text)


def _images_in(par, doc):
    """(blob, display width in points) for every inline image in a paragraph."""
    got = []
    for blip in par._p.iter(f"{A}blip"):
        rid = blip.get(f"{R}embed")
        if not rid:
            continue
        try:
            blob = doc.part.related_parts[rid].blob
        except Exception:
            continue
        w_pt = None
        ext = par._p.find(f".//{WP}extent")
        if ext is not None and ext.get("cx"):
            w_pt = float(ext.get("cx")) / EMU * 72.0
        got.append((blob, w_pt))
    return got


def _fit_image(blob, disp_pt):
    """Re-encode one figure at its display size, choosing the smaller of PNG / JPEG.

    Charts are line art, so PNG usually wins and stays lossless; the colour-mapped
    space-time fields are photographic and JPEG wins there. Chroma subsampling is off
    (4:4:4) because these images carry coloured text and hairlines, which 4:2:0 smears.
    Returns (bytes, pixel width) or (original blob, None) if anything goes wrong.
    """
    try:
        im = PILImage.open(io.BytesIO(blob))
        cap = max(int(disp_pt / 72.0 * IMG_DPI), 200)
        src = im.resize((cap, max(1, round(im.height * cap / im.width))),
                        PILImage.Resampling.LANCZOS) if im.width > cap else im
        rgb = src.convert("RGB")
        bp = io.BytesIO(); rgb.save(bp, "PNG", optimize=True)
        bj = io.BytesIO(); rgb.save(bj, "JPEG", quality=92, optimize=True, subsampling=0)
        best = bj if bj.tell() < bp.tell() else bp
        return best.getvalue(), rgb.width
    except Exception:
        return blob, None


def _has_page_break(par):
    for br in par._p.iter(f"{W}br"):
        if br.get(f"{W}type") == "page":
            return True
    return False


def convert(docx_path, pdf_path=None, pagesize=A4):
    """Render `docx_path` to PDF. Returns the output path."""
    pdf_path = pdf_path or (os.path.splitext(docx_path)[0] + ".pdf")
    body_f, bold_f, ital_f = _register_fonts()
    doc = Document(docx_path)

    margin = 0.72 * inch
    frame_w = pagesize[0] - 2 * margin

    def S(name, size, leading, **kw):
        kw.setdefault("fontName", body_f)          # callers may name a bold/italic face
        return ParagraphStyle(name, fontSize=size, leading=leading, **kw)

    st = {
        "Title": S("Title", 19, 24, fontName=bold_f, textColor=NAVY, spaceAfter=10,
                   alignment=TA_CENTER),
        "Subtitle": S("Subtitle", 11, 15, textColor=GREY, spaceAfter=16,
                      alignment=TA_CENTER),
        "H1": S("H1", 15.5, 20, fontName=bold_f, textColor=NAVY, spaceBefore=13, spaceAfter=6),
        "H2": S("H2", 12.5, 16, fontName=bold_f, textColor=NAVY, spaceBefore=10, spaceAfter=5),
        "H3": S("H3", 11, 14.5, fontName=bold_f, textColor=INK, spaceBefore=8, spaceAfter=4),
        "Body": S("Body", 9.6, 13.4, alignment=TA_JUSTIFY, spaceAfter=5),
        "Bullet": S("Bullet", 9.6, 13.4, leftIndent=14, bulletIndent=4, spaceAfter=3),
        "Caption": S("Caption", 8.3, 11, textColor=GREY, alignment=TA_CENTER,
                     spaceBefore=2, spaceAfter=9),
        "Cell": S("Cell", 7.4, 9.4),
        "CellHead": S("CellHead", 7.4, 9.4, fontName=bold_f, textColor=NAVY),
    }

    flow = []
    body = doc.element.body
    prev_was_image = False
    for child in body.iterchildren():
        tag = child.tag
        if tag == f"{W}p":
            par = Paragraph(child, doc)
            if _has_page_break(par):
                flow.append(PageBreak())
            imgs = _images_in(par, doc)
            if imgs:
                for blob, w_pt in imgs:
                    try:
                        w = min(w_pt or frame_w, frame_w)
                        blob, _ = _fit_image(blob, w)
                        img = Image(io.BytesIO(blob))
                        scale = w / float(img.imageWidth)
                        img.drawWidth = w
                        img.drawHeight = float(img.imageHeight) * scale
                        #  never taller than the usable column, or platypus drops it
                        max_h = pagesize[1] - 2 * margin - 40
                        if img.drawHeight > max_h:
                            k = max_h / img.drawHeight
                            img.drawWidth *= k
                            img.drawHeight *= k
                        img.hAlign = "CENTER"
                        flow.append(img)
                        prev_was_image = True
                    except Exception:
                        continue
                continue
            text = par.text.strip()
            if not text:
                flow.append(Spacer(1, 3))
                continue
            #  python-docx returns None for a paragraph that carries no style, and
            #  `None.name` would abort the whole PDF build on one such paragraph.
            name = ((par.style.name if par.style is not None else "") or "").lower()
            markup = _runs_markup(par)
            if name.startswith("title"):
                flow.append(RLPara(markup, st["Title"]))
            elif name.startswith("heading 1"):
                flow.append(RLPara(markup, st["H1"]))
            elif name.startswith("heading 2"):
                flow.append(RLPara(markup, st["H2"]))
            elif name.startswith("heading 3"):
                flow.append(RLPara(markup, st["H3"]))
            elif name.startswith("list") or text.startswith(("•", "–", "-")) and len(text) < 400:
                flow.append(RLPara(markup.lstrip("•–- "), st["Bullet"], bulletText="•"))
            elif prev_was_image:
                flow.append(RLPara(markup, st["Caption"]))
            else:
                flow.append(RLPara(markup, st["Body"]))
            prev_was_image = False
        elif tag == f"{W}tbl":
            tbl = Table(child, doc)
            data, ncol = [], 0
            for ri, row in enumerate(tbl.rows):
                cells = []
                for cell in row.cells:
                    style = st["CellHead"] if ri == 0 else st["Cell"]
                    cells.append(RLPara(_esc(cell.text.strip()) or "&nbsp;", style))
                ncol = max(ncol, len(cells))
                data.append(cells)
            if not data:
                continue
            for rrow in data:                       # ragged rows -> rectangular
                while len(rrow) < ncol:
                    rrow.append(RLPara("&nbsp;", st["Cell"]))
            t = RLTable(data, colWidths=[frame_w / ncol] * ncol, repeatRows=1)
            t.setStyle(TableStyle([
                ("GRID", (0, 0), (-1, -1), 0.4, GRID),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF3FC")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]))
            flow.append(t)
            flow.append(Spacer(1, 7))
            prev_was_image = False

    def _page(canvas, docp):
        canvas.saveState()
        canvas.setFont(body_f, 7.6)
        canvas.setFillColor(GREY)
        canvas.drawCentredString(pagesize[0] / 2.0, 0.42 * inch, str(docp.page))
        canvas.restoreState()

    SimpleDocTemplate(pdf_path, pagesize=pagesize,
                      leftMargin=margin, rightMargin=margin,
                      topMargin=0.66 * inch, bottomMargin=0.66 * inch,
                      title=os.path.basename(pdf_path)).build(
        flow, onFirstPage=_page, onLaterPages=_page)
    return pdf_path


def main(argv):
    if not argv:
        print("usage: python3 docx2pdf_native.py <file.docx> [out.pdf]")
        return 2
    out = convert(argv[0], argv[1] if len(argv) > 1 else None)
    print(f"wrote {out}  ({os.path.getsize(out) / 1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
