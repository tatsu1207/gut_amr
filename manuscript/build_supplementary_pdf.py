#!/usr/bin/env python3
"""Build the single merged Supplementary Information PDF required by npj Biofilms and Microbiomes.

Large data tables stay in supplementary_tables.xlsx as Supplementary Data 1-14 (titles and
legends listed here); figures and Supplementary Table 1 are rendered inline.
"""
import os
import sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tplabels
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                PageBreak, Table, TableStyle)

XL = "/home/unnot/vault/KJM_info/manuscript/supplementary_tables.xlsx"
FIGM = "/home/unnot/vault/kjm/miseq/dada2_output/figures"
OUT = "/home/unnot/vault/KJM_info/manuscript/Supplementary_Information.pdf"

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontName="Times-Bold", fontSize=13, spaceAfter=10)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontName="Times-Bold", fontSize=11, spaceBefore=12, spaceAfter=6)
BODY = ParagraphStyle("B", parent=ss["BodyText"], fontName="Times-Roman", fontSize=10, leading=14)
LEG = ParagraphStyle("L", parent=BODY, fontSize=9, leading=12, spaceBefore=6)
SMALL = ParagraphStyle("S", parent=BODY, fontSize=8, leading=10)

doc = SimpleDocTemplate(OUT, pagesize=A4, topMargin=0.8*inch, bottomMargin=0.8*inch,
                        leftMargin=0.8*inch, rightMargin=0.8*inch,
                        title="Supplementary Information")
S = []

S.append(Paragraph("Supplementary Information", H1))
S.append(Paragraph(
    "Intrinsic, chromosomally encoded resistance dominates the gut resistome, and high-fat "
    "diet is associated with failed recovery of butyrate-synthesis capacity after amoxicillin in mice", BODY))
S.append(Paragraph("Jungman Kim, Tatsuya Unno", BODY))
S.append(Spacer(1, 16))

S.append(Paragraph("Contents", H2))
contents = [
    ["Supplementary Fig. 1", "Cecal short-chain fatty acid concentrations by group"],
    ["Supplementary Fig. 2", "Rarefaction curves and library size distribution"],
    ["Supplementary Fig. 3", "Genus composition of each group at each timepoint"],
    ["Supplementary Fig. 4", "The 25 individual ARGs that changed most"],
    ["Supplementary Table 1", "PERMANOVA results across timepoints"],
    ["Supplementary Data 1-14", "Titles and legends; the data are provided as a separate Excel workbook "
                                "(supplementary_tables.xlsx), one or more sheets per item, sheet names "
                                "prefixed SD1-SD14"],
]
t = Table([[Paragraph(f"<b>{a}</b>", SMALL), Paragraph(b, SMALL)] for a, b in contents],
          colWidths=[1.6*inch, 4.6*inch])
t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                       ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
S.append(t)
S.append(PageBreak())


def fig(path, w_in, title, legend):
    img = Image(path)
    r = (w_in * inch) / img.imageWidth
    img._restrictSize(img.imageWidth * r, img.imageHeight * r)
    S.append(Paragraph(title, H2))
    S.append(img)
    S.append(Paragraph(legend, LEG))
    S.append(PageBreak())


fig(f"{FIGM}/figS2_scfa.png", 6.6, "Supplementary Fig. 1",
    "<b>Supplementary Fig. 1. Cecal short-chain fatty acid concentrations.</b> "
    "(a) Acetate, (b) propionate, and (c) butyrate measured by gas chromatography at the end of the "
    "experiment. Each point is one mouse; horizontal bars are group means. Blue, normal diet; red, "
    "high-fat diet. P values are from a Kruskal-Wallis test across all eight groups. This figure shows "
    "the per-animal distribution underlying the summary box plots of Fig. 1d-f.")

fig(f"{FIGM}/figS1_rarefaction.png", 6.3, "Supplementary Fig. 2",
    "<b>Supplementary Fig. 2. Sequencing depth and rarefaction.</b> "
    "(a) Rarefaction curves for a random subset of 150 of the 451 retained samples, showing expected "
    "ASV richness as a function of sequencing depth; the dashed line marks the rarefaction depth of "
    "5,159 reads used for alpha diversity. (b) Distribution of library sizes across all retained "
    "samples. Samples with fewer than 5,000 reads were excluded before rarefaction (Supplementary Data 1).")

fig(os.path.expanduser("~/vault/kjm/miseq/dada2_output/figures/taxonomy_stacked_bars.png"), 6.6, "Supplementary Fig. 3",
    "<b>Supplementary Fig. 3. Genus-level composition by group and timepoint.</b> "
    "Stacked bars show the mean relative abundance of the 15 most abundant genera (remaining genera as 'Other') "
    "for each diet-antibiotic group at each 16S rRNA gene sampling timepoint. Relative abundances are "
    "percentages of the reads assigned at genus level in each sample. CTL and HFD are the untreated "
    "normal-diet and high-fat-diet groups. R7 is shown for completeness only: its samples were sequenced alone "
    "on a separate MiSeq run that performed poorly (33% of reads passing quality filtering, compared with "
    "52-79% at the other timepoints), and untreated mice showed abrupt, transient changes at R7 (for example, "
    "Akkermansia falling to 0%); R7 is therefore not interpreted (Methods).")

fig(os.path.expanduser("~/vault/kjm/hiseq/figures/fig4d_top_args_heatmap.png"), 6.6, "Supplementary Fig. 4",
    "<b>Supplementary Fig. 4. Individual ARGs with the largest changes.</b> "
    "Heatmap of mean read-level abundance (ARG copies per 16S rRNA gene copy; CARD protein homolog models "
    "only) for the 25 ARGs with the largest change in group mean between T0 and T3 in any "
    "antibiotic-treated group, shown for every group and timepoint (n = 3 mice each). Descriptive; "
    "no statistical test was applied.")


def table(df, title, legend, maxrows=None, fontsize=7, landscape_page=False, colw=None):
    S.append(Paragraph(title, H2))
    d = df if maxrows is None else df.head(maxrows)
    data = [[Paragraph(f"<b>{c}</b>", ParagraphStyle('h', parent=SMALL, fontSize=fontsize))
             for c in d.columns]]
    for _, r in d.iterrows():
        row = []
        for v in r:
            if isinstance(v, float):
                v = f"{v:.3g}"
            v = str(v)
            v = tplabels.label(v) if not v.startswith("ARG_") else v
            if len(v) > 110:                       # keep free-text cells from blowing up row height
                v = v[:107].rsplit(";", 1)[0] + " ..."
            row.append(Paragraph(v, ParagraphStyle('c', parent=SMALL, fontSize=fontsize)))
        data.append(row)
    tb = Table(data, colWidths=colw, repeatRows=1)
    tb.setStyle(TableStyle([
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#999999")),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    S.append(tb)
    S.append(Paragraph(legend, LEG))
    S.append(PageBreak())


perm = pd.read_excel(XL, sheet_name="SupplTable1_PERMANOVA")
table(perm, "Supplementary Table 1",
      "<b>Supplementary Table 1.</b> PERMANOVA of gut microbiome community composition among the eight "
      "groups at each timepoint (Bray-Curtis dissimilarity, 999 permutations).")

# Supplementary Data: title and legend for each item (npj: large tables are labelled
# "Supplementary Data N" and their titles/legends go in this PDF). Legends are read from
# the Supplementary information list in manuscript.txt so the two cannot drift apart.
import re
from openpyxl import load_workbook
MS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "manuscript.txt")
items = [l.strip() for l in open(MS) if re.match(r"Supplementary Data \d+\. ", l)]
sheets = load_workbook(XL, read_only=True).sheetnames
S.append(Paragraph("Supplementary Data", H2))
S.append(Paragraph("Provided as a separate Excel workbook (supplementary_tables.xlsx). Each item occupies "
                   "the sheet(s) whose name begins with its number (SD1, SD2a, ...).", BODY))
for it in items:
    n = re.match(r"Supplementary Data (\d+)\.", it).group(1)
    title, legend = it.split(". ", 1)
    sh = [x for x in sheets if re.match(rf"SD{n}(?:[a-z])?_", x)]
    assert sh, n
    S.append(Paragraph(f"<b>{title}.</b> {legend} <i>Sheets: {', '.join(sh)}.</i>", LEG))

doc.build(S)
print("wrote", OUT)
