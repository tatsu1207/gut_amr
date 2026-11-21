#!/usr/bin/env python3
"""Rebuild the submission manuscript.docx from manuscript.txt.

Uses an existing manuscript.docx as a style template (page setup, Normal /
Heading styles, footer page numbers), clears its body, and re-fills it from the
plain-text source. Headings are recognised from explicit lists, not guessed from
line length, so body paragraphs can never be mis-styled as headings.

Usage: python3 build_manuscript_docx.py manuscript.txt template.docx out.docx
"""
import sys
import copy
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

src, template, out = sys.argv[1:4]

H1 = {"Abstract", "Background", "Results", "Discussion", "Conclusions", "Methods",
      "List of abbreviations", "Declarations", "References", "Figure legends",
      "Supplementary information"}
H2_SECTIONS = {"Results", "Methods", "Declarations"}
DECL_H2 = {"Ethics approval and consent to participate", "Consent for publication",
           "Competing interests", "Funding", "Authors' contributions",
           "Acknowledgements", "Availability of data and materials"}
ABSTRACT_LABELS = ("Background:", "Results:", "Conclusions:")

lines = [l.rstrip("\n") for l in open(src, encoding="utf-8")]

doc = Document(template)
body = doc.element.body
sectPr = body.find('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}sectPr')
for el in list(body):
    if el is not sectPr:
        body.remove(el)


def para(text, style="Normal", align=None, size=None, bold=None, label=None):
    p = doc.add_paragraph(style=style)
    if align is not None:
        p.alignment = align
    if label:  # bold run label followed by plain text
        r = p.add_run(label)
        r.bold = True
        if size: r.font.size = Pt(size)
        text = text[len(label):]
    r = p.add_run(text)
    if bold: r.bold = True
    if size: r.font.size = Pt(size)
    return p


def is_h2(text, section):
    if section == "Declarations":
        return text in DECL_H2
    if section not in H2_SECTIONS:
        return False
    return (len(text) < 100 and not text.endswith(".") and not text.startswith("(")
            and ":" not in text)


# --- Title block: everything before "Abstract" ---
first_h1 = next(i for i, l in enumerate(lines) if l.strip() == "Abstract")
title_block = [l.strip() for l in lines[:first_h1] if l.strip()]
para(title_block[0], align=WD_ALIGN_PARAGRAPH.CENTER, size=16, bold=True)
for l in title_block[1:]:
    size = 11 if l.startswith("Jungman") else 10
    para(l, align=WD_ALIGN_PARAGRAPH.CENTER, size=size)

section = None
counts = {"Heading 1": 0, "Heading 2": 0, "Normal": 0}
for raw in lines[first_h1:]:
    t = raw.strip()
    if not t:
        continue
    if t in H1:
        section = t
        para(t, style="Heading 1"); counts["Heading 1"] += 1
    elif is_h2(t, section):
        para(t, style="Heading 2"); counts["Heading 2"] += 1
    elif section == "Abstract" and t.startswith(ABSTRACT_LABELS):
        para(t, label=t.split(":")[0] + ":"); counts["Normal"] += 1
    elif t.startswith("Keywords:"):
        para(t, label="Keywords:"); counts["Normal"] += 1
    else:
        para(t); counts["Normal"] += 1

doc.save(out)
print(out, counts)
