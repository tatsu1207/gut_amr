#!/usr/bin/env python3
"""Build the single merged Supplementary Information PDF required by npj Biofilms and Microbiomes.

Large data tables stay in supplementary_tables.xlsx and are cited here as
Supplementary Data; only figures and small tables are rendered inline.
"""
import os
import pandas as pd
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
    "diet impairs recovery of butyrate-synthesis capacity after antibiotics in mice", BODY))
S.append(Paragraph("Jungman Kim, Tatsuya Unno", BODY))
S.append(Spacer(1, 16))

S.append(Paragraph("Contents", H2))
contents = [
    ["Supplementary Figure S1", "Rarefaction curves and library size distribution"],
    ["Supplementary Figure S2", "Cecal short-chain fatty acid concentrations by group"],
    ["Supplementary Figure S3", "Genus composition of each group at each timepoint"],
    ["Supplementary Figure S4", "The 25 individual ARGs that changed most"],
    ["Supplementary Table 1", "PERMANOVA results across timepoints"],
    ["Supplementary Table S4", "ARG counts by bacterial genus (assembly-based)"],
    ["Supplementary Table S6c", "SD14 butyrate-gene abundance vs measured cecal SCFAs"],
    ["Supplementary Table S6d", "Butyrate gene / SCFA correlation statistics"],
    ["Supplementary Table S7a", "ARG association with mobile genetic elements by drug class"],
    ["Supplementary Data", "Tables S1, S2 (CARD protein homolog models; rRNA-variant-excluded and unfiltered versions alongside), "
                           "S3, S5, S6a, S6b, S6e, S7b, S7c, S8-S14 "
                           "(provided as a separate Excel workbook; too large to render here)"],
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


fig(f"{FIGM}/figS1_rarefaction.png", 6.3, "Supplementary Figure S1",
    "<b>Supplementary Fig. S1. Sequencing depth and rarefaction.</b> "
    "(a) Rarefaction curves for a random subset of 150 of the 498 retained samples, showing expected "
    "ASV richness as a function of sequencing depth; the dashed line marks the rarefaction depth of "
    "5,159 reads used for alpha diversity. (b) Distribution of library sizes across all retained "
    "samples. Samples with fewer than 5,000 reads were excluded before rarefaction.")

fig(f"{FIGM}/figS2_scfa.png", 6.6, "Supplementary Figure S2",
    "<b>Supplementary Fig. S2. Cecal short-chain fatty acid concentrations.</b> "
    "(a) Acetate, (b) propionate, and (c) butyrate measured by gas chromatography at the end of the "
    "experiment. Each point is one mouse; horizontal bars are group means. Blue, normal diet; red, "
    "high-fat diet. P values are from a Kruskal-Wallis test across all eight groups. This figure shows "
    "the per-animal distribution underlying the summary box plots of Fig. 1d-f.")

fig(os.path.expanduser("~/vault/kjm/miseq/dada2_output/figures/taxonomy_stacked_bars.png"), 6.6, "Supplementary Figure S3",
    "<b>Supplementary Fig. S3. Genus-level composition by group and timepoint.</b> "
    "Stacked bars show the mean relative abundance of the 15 most abundant genera (remaining genera as 'Other') "
    "for each diet-antibiotic group at each 16S rRNA gene sampling timepoint. Relative abundances are "
    "percentages of the reads assigned at genus level in each sample. CTL and HFD are the untreated "
    "normal-diet and high-fat-diet groups.")

fig(os.path.expanduser("~/vault/kjm/hiseq/figures/fig4d_top_args_heatmap.png"), 6.6, "Supplementary Figure S4",
    "<b>Supplementary Fig. S4. Individual ARGs with the largest changes.</b> "
    "Heatmap of mean read-level abundance (ARG copies per 16S rRNA gene copy; CARD protein homolog models "
    "only) for the 25 ARGs with the largest change in group mean between T0h and T72h in any "
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


perm = pd.read_excel(XL, sheet_name="Table1_PERMANOVA")
table(perm, "Supplementary Table 1",
      "<b>Supplementary Table 1.</b> PERMANOVA results for gut microbiome community composition "
      "across timepoints (Bray-Curtis dissimilarity, 999 permutations).")

genus = pd.read_excel(XL, sheet_name="TableS4_ARG_by_genus")
table(genus, "Supplementary Table S4",
      "<b>Supplementary Table S4.</b> Assembly-based ARG counts by bacterial genus, from GTDB-Tk "
      "classification of medium-quality or better MAGs. The 40 genera with the most ARGs are shown; "
      f"the complete table of {len(genus)} genera is provided as Supplementary Data.",
      maxrows=40, colw=[1.05*inch, 0.6*inch, 0.6*inch, 0.5*inch, 3.55*inch])

s7c = pd.read_excel(XL, sheet_name="TableS6c_SD14_vs_cecal_SCFA")
table(s7c, "Supplementary Table S6c",
      "<b>Supplementary Table S6c.</b> Mean butyrate-synthesis gene abundance (copies per genome) at "
      "SD14 for each experimental group, alongside the mean cecal acetate, propionate, and butyrate "
      "concentrations (micromolar) measured by gas chromatography.")

s7d = pd.read_excel(XL, sheet_name="TableS6d_SCFA_correlation")
table(s7d, "Supplementary Table S6d",
      "<b>Supplementary Table S6d.</b> Correlation between butyrate-synthesis gene abundance at SD14 "
      "and measured cecal SCFA concentrations, across the eight group means.")

s8a = pd.read_excel(XL, sheet_name="TableS7a_mobilome_by_class")
table(s8a, "Supplementary Table S7a",
      "<b>Supplementary Table S7a.</b> Association of assembly-based ARGs with mobile genetic elements "
      "by antibiotic class. An ARG was counted as mobile if its contig was predicted to be a plasmid or "
      "phage, or carried a transposase, integrase, recombinase, relaxase, or conjugation gene. The 25 "
      "most abundant drug-class assignments are shown; the full table is provided as Supplementary Data.",
      maxrows=25, fontsize=6.5, colw=[2.9*inch, 0.6*inch, 0.5*inch, 0.6*inch, 0.6*inch, 0.75*inch, 0.45*inch])

doc.build(S)
print("wrote", OUT)
