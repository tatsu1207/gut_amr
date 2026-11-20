#!/usr/bin/env python3
"""Combine the separate panel PNGs of Figs 3, 4 and 5 into one labelled image per figure.

npj guidelines (2026-11): each figure legend must sit on the same page as its figure in
the combined review PDF, and at acceptance each figure is one file containing all panels.
Panels are scaled so that every image in a row shares a height, rows are scaled to a
common width, and a bold panel letter matching the legend is drawn at each panel's
top-left corner. Output: <name>_combined.png (300 dpi) next to the panel files.
"""
import os
from PIL import Image, ImageDraw, ImageFont

H = os.path.expanduser
MF = H("~/vault/kjm/miseq/dada2_output/figures/")
HF = H("~/vault/kjm/hiseq/figures/")
FONT = "/usr/share/fonts/truetype/msttcorefonts/arialbd.ttf"
WIDTH = 4200          # px; 14 in at 300 dpi, reduced to page width in the PDF
GAP = 40              # px between panels and rows

FIGS = {
    MF + "fig3_combined.png": [
        # former panel A (stacked bars by group x timepoint) is now Supplementary Fig. S3
        [("A", MF + "genus_trajectories.png"), ("B", MF + "genus_heatmap.png")],
    ],
    HF + "fig4_combined.png": [
        [("A", HF + "fig4a_total_arg.png")],
        [("B", HF + "fig4b_arg_class_composition.png")],
        [("C", HF + "fig4c_arg_foldchange_heatmap.png")],
        # former panels: D (top-25 gene heatmap) is now Supplementary Fig. S3; E (diet bars) was
        # merged into A, which now plots ND and HFD side by side.
    ],
    HF + "fig5_combined.png": [
        [("A", HF + "fig5a_assembly_arg_classes.png"), ("B", HF + "fig5b_top_arg_genera.png"),
         ("D", HF + "fig5d_genomic_context.png")],
        [("C", HF + "fig5c_lactococcus_vs_mls.png"), ("E", HF + "fig5e_assembly_args_by_diet.png")],
    ],
}


def build_row(panels):
    ims = [Image.open(p).convert("RGB") for _, p in panels]
    h = max(im.height for im in ims)
    ims = [im.resize((round(im.width * h / im.height), h), Image.LANCZOS) for im in ims]
    w = sum(im.width for im in ims) + GAP * (len(ims) - 1)
    scale = WIDTH / w
    h2 = round(h * scale)
    TOP = 100  # white strip above the panels so letters never overlap panel titles
    row = Image.new("RGB", (WIDTH, h2 + TOP), "white")
    draw = ImageDraw.Draw(row)
    font = ImageFont.truetype(FONT, 90)
    x = 0
    for (label, _), im in zip(panels, ims):
        im = im.resize((round(im.width * scale), h2), Image.LANCZOS)
        row.paste(im, (x, TOP))
        draw.text((x + 8, 0), label, fill="black", font=font)
        x += im.width + round(GAP * scale)
    return row


for out, rows in FIGS.items():
    built = [build_row(r) for r in rows]
    total_h = sum(r.height for r in built) + GAP * 2 * (len(built) - 1)
    canvas = Image.new("RGB", (WIDTH, total_h), "white")
    y = 0
    for r in built:
        canvas.paste(r, (0, y))
        y += r.height + GAP * 2
    canvas.save(out, dpi=(300, 300))
    print(f"{out}: {canvas.size[0]}x{canvas.size[1]} px")
