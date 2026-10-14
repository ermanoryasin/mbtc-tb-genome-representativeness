#!/usr/bin/env python3
"""Packages the Epidemiology and Infection submission files with figures numbered by first citation, and builds
(1) the journal manuscript (figures/tables at the end) and (2) a single-file review copy with figures and tables embedded."""
import json, re, shutil, os
from pathlib import Path
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
SUB = HERE / "submission_EI"; UP = SUB / "upload"; UP.mkdir(exist_ok=True)
FIGMAP = {int(k): v for k, v in json.load(open(HERE / "figure_map.json")).items()}
SLUG = {1: "concentration", 2: "coverage", 3: "temporal", 4: "regression", 5: "turkey", 6: "study_dependence"}
STEM = {1: "Fig1_concentration", 2: "Fig2_coverage", 3: "Fig3_temporal", 4: "Fig4_regression", 5: "Fig5_turkey", 6: "Fig6_study_dependence"}

# ---- 1. figure files with new numbers ----
for f in UP.glob("03_Figure_*"): f.unlink()
newpng = {}
for old, new in FIGMAP.items():
    for ext in ("tiff", "png"):
        src = HERE / "figures" / f"{STEM[old]}.{ext}"
        dst = (UP if ext == "tiff" else SUB / "figures_by_manuscript_number") / f"03_Figure_{new}_{SLUG[old]}.{ext}"
        dst.parent.mkdir(exist_ok=True); shutil.copy(src, dst)
        if ext == "png": newpng[new] = dst

# ---- 2. alt text renumbered ----
alt = open(SUB / "figure_alt_text.md", encoding="utf-8").read().split("\n")
head = [l for l in alt if l.startswith("# ")]
paras = [l for l in alt if l.startswith("**Figure")]
if all(re.match(r"\*\*Figure (\d)\.\*\*", p) for p in paras):
    # alt text is stored by OLD numbers; renumber once and mark file
    mapped = [re.sub(r"^\*\*Figure (\d)\.\*\*", lambda m: f"**Figure {FIGMAP[int(m.group(1))]}.**", p) for p in paras]
    if not any(l.startswith("<!--renumbered") for l in alt):
        mapped.sort(key=lambda p: int(re.match(r"\*\*Figure (\d)", p).group(1)))
        open(SUB / "figure_alt_text.md", "w", encoding="utf-8").write("\n".join(head + [""] + ["<!--renumbered by first citation-->"] + [""] + [x + "\n" for x in mapped]))
alt_text = open(SUB / "figure_alt_text.md", encoding="utf-8").read().split("\n")

def base_doc(font=12, spacing=2.0, line_numbers=True, header=None):
    d = Document(); s = d.sections[0]
    for a in ("left_margin", "right_margin", "top_margin", "bottom_margin"): setattr(s, a, Inches(1))
    st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(font)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    st.paragraph_format.line_spacing = spacing; st.paragraph_format.space_after = Pt(0 if spacing >= 2 else 6)
    if line_numbers:
        ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); s._sectPr.append(ln)
    if header:
        hp = s.header.paragraphs[0]; hp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = hp.add_run(header); r.italic = True; r.font.size = Pt(9)
    p = s.footer.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.RIGHT; r = p.add_run()
    for t, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if t: e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), t)
        else: e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt
        r._r.append(e)
    return d

def runs(par, text, bold=False, size=None):
    for ch in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
        if not ch: continue
        if ch.startswith("**"): r = par.add_run(ch[2:-2]); r.bold = True
        elif ch.startswith("*") and len(ch) > 2: r = par.add_run(ch[1:-1]); r.italic = True; r.bold = bold
        else: r = par.add_run(ch); r.bold = bold
        if size: r.font.size = Pt(size)

def add_par(d, text, bold=False, size=None, align=WD_ALIGN_PARAGRAPH.JUSTIFY, single=False, italic_all=False):
    p = d.add_paragraph(); runs(p, text, bold, size); p.alignment = align
    if single: p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    if italic_all:
        for r in p.runs: r.italic = True
    return p

def add_table(d, rows, size=9, widths=None):
    t = d.add_table(rows=len(rows), cols=len(rows[0])); t.style = "Table Grid" if False else t.style
    for i, row in enumerate(rows):
        for j, c in enumerate(row):
            cell = t.cell(i, j); cell.text = ""; p = cell.paragraphs[0]
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE; p.paragraph_format.space_after = Pt(2)
            runs(p, str(c), bold=(i == 0), size=size)
    return t

md = open(HERE / "manuscript_EI.md", encoding="utf-8").read()
lines = md.split("\n")
# split md sections
def block(a, b=None):
    i = md.index(a); j = md.index(b, i) if b else len(md); return md[i:j]
front = block("# ", "## Introduction")
main = block("## Introduction", "**Acknowledgements.**")
tail = block("**Acknowledgements.**", "## References")
refs = block("## References", "## Tables")
tables_md = block("## Tables", "## Figure legends")
legends_md = block("## Figure legends")
legends = {int(re.match(r"\*\*Figure (\d)", p).group(1)): p for p in legends_md.split("\n\n") if p.startswith("**Figure")}
tbl_blocks = re.split(r"(?=\*\*Table \d\.)", tables_md)[1:]
def parse_table(blk):
    ls = [l for l in blk.strip().split("\n")]
    title = ls[0]; rows = []; foot = []
    for l in ls[1:]:
        if l.startswith("|"):
            cells = [c.strip() for c in l.strip().strip("|").split("|")]
            if not all(set(c) <= set("- ") for c in cells): rows.append(cells)
        elif l.strip(): foot.append(l.strip())
    return title, rows, " ".join(foot)
TABLES = {int(re.match(r"\*\*Table (\d)", b).group(1)): parse_table(b) for b in tbl_blocks}

def render_md_block(d, text, embed=None, spacing_single=False):
    """Render markdown-ish text to doc. embed(paragraph_text)-> callable invoked after a paragraph."""
    shown_fig, shown_tab = set(), set()
    for blk in text.strip().split("\n\n"):
        blk = blk.strip()
        if not blk: continue
        if blk.startswith("### "): p = add_par(d, blk[4:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT, single=spacing_single); [setattr(r, "italic", True) for r in p.runs]
        elif blk.startswith("## "): add_par(d, blk[3:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT, single=spacing_single)
        elif blk.startswith("# "): add_par(d, blk[2:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT, single=spacing_single)
        elif blk.startswith("- "): [add_par(d, "• " + l[2:], single=spacing_single) for l in blk.split("\n")]
        else:
            add_par(d, blk.replace("\n", " "), single=spacing_single)
            if embed:
                for kind, n in embed(blk, shown_fig, shown_tab): yield_ = None
    return

# =============== A. journal manuscript (as before: tables and legends at the end) ===============
# (kept by build_docx_ei.py; this script produces only the review copy below)

# =============== B. review copy with embedded figures and tables ===============
d = base_doc(font=11, spacing=1.5, line_numbers=True, header="Manuscript for external review. Confidential; please do not distribute.")
def figure_block(n):
    p = d.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(newpng[n]), width=Inches(6.2)); p.paragraph_format.keep_with_next = True
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    add_par(d, legends[n], size=9.5, single=True)
def table_block(n):
    title, rows, foot = TABLES[n]
    add_par(d, title, size=10, single=True, align=WD_ALIGN_PARAGRAPH.LEFT).paragraph_format.keep_with_next = True
    add_table(d, rows, size=8 if n == 2 else 9)
    add_par(d, foot, size=9, single=True)
    d.add_paragraph().paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE

# front matter
for blk in front.strip().split("\n\n"):
    blk = blk.strip()
    if blk.startswith("# "): add_par(d, blk[2:], bold=True, size=14, align=WD_ALIGN_PARAGRAPH.LEFT)
    elif blk.startswith("## "): add_par(d, blk[3:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT)
    elif blk.startswith("- "): [add_par(d, "• " + l[2:]) for l in blk.split("\n")]
    else: add_par(d, blk.replace("\n", " "), align=WD_ALIGN_PARAGRAPH.LEFT if len(blk) < 220 else WD_ALIGN_PARAGRAPH.JUSTIFY)
# main text with embedded items placed after the first citing paragraph
done_fig, done_tab = set(), set()
for blk in main.strip().split("\n\n"):
    blk = blk.strip()
    if not blk: continue
    if blk.startswith("## "): add_par(d, blk[3:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT)
    elif blk.startswith("### "): p = add_par(d, blk[4:], bold=True, align=WD_ALIGN_PARAGRAPH.LEFT); [setattr(r, "italic", True) for r in p.runs]
    else:
        add_par(d, blk.replace("\n", " "))
        for n in sorted({int(x) for x in re.findall(r"Figure[s]? (\d)", blk)}):
            if n not in done_fig: done_fig.add(n); figure_block(n)
        for n in sorted({int(x) for x in re.findall(r"Table (\d)(?![\d-])", re.sub(r"Supplementary Table S\d", "", blk))}):
            if n in TABLES and n not in done_tab: done_tab.add(n); table_block(n)
for blk in tail.strip().split("\n\n"):
    add_par(d, blk.strip().replace("\n", " "))
# references
add_par(d, "References", bold=True, align=WD_ALIGN_PARAGRAPH.LEFT)
for l in refs.split("\n"):
    if re.match(r"^\d+\. ", l): add_par(d, l, size=10, align=WD_ALIGN_PARAGRAPH.LEFT, single=True)
# supplementary tables
d.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
add_par(d, "Supplementary tables (selected)", bold=True, align=WD_ALIGN_PARAGRAPH.LEFT)
add_par(d, "Full set (S1 to S9, including regression output S5, negative binomial models S7 and country-level expected counts S9) is in the Supplementary_Tables.xlsx workbook and the repository (doi:10.5281/zenodo.23376282).", size=10, single=True)
wb = load_workbook(SUB / "Supplementary_Tables.xlsx")
for name in ["S1 study cap", "S2 region x release", "S3 model specification", "S4 study bootstrap", "S6 Turkey provenance", "S8 Turkey NB expectations"]:
    ws = wb[name]; rows = [[("" if v is None else v) for v in r] for r in ws.iter_rows(values_only=True)]
    rows = [r for r in rows if any(str(c).strip() for c in r)]
    add_par(d, f"Supplementary Table {name.split()[0]}: {' '.join(name.split()[1:])}", bold=True, size=10, single=True, align=WD_ALIGN_PARAGRAPH.LEFT).paragraph_format.keep_with_next = True
    add_table(d, [[(f"{c:g}" if isinstance(c, float) else str(c))[:58] for c in r] for r in rows], size=6.5 if len(rows[0]) > 7 else 8)
    d.add_paragraph().paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
d.save(SUB / "Manuscript_for_external_review.docx")
print("review copy ok; figures embedded:", sorted(done_fig), "tables:", sorted(done_tab))
