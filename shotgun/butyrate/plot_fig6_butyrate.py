#!/usr/bin/env python3
"""Fig. 6 - Shotgun-metagenomic butyrate-synthesis capacity vs measured cecal butyrate."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.expanduser("~/vault/kjm/miseq"))
import tplabels; tplabels.install()  # show T0/T1/T3/R3...R14 instead of the data codes
from scipy import stats
from PIL import Image

B = "/home/unnot/vault/kjm/hiseq/butyrate_analysis"
T = f"{B}/tables"
FIG = "/home/unnot/vault/kjm/hiseq/figures"
SUB = "/home/unnot/vault/KJM_info/manuscript"

ND, HFD = "#3070B8", "#CC6677"          # validated pair (chroma, CVD, contrast)
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})

ab = pd.read_csv(f"{T}/butyrate_gene_abundance_per_sample.tsv", sep="\t", keep_default_na=False)
link = pd.read_csv(f"{T}/sd14_gene_vs_cecal_scfa_group_means.tsv", sep="\t", index_col=0)
gen = pd.read_csv(f"{T}/but_buk_carrier_genera.tsv", sep="\t")

TP = ["T0h", "T72h", "SD14"]
ABX = [("None", "No antibiotic"), ("AMX", "Amoxicillin"), ("CIP", "Ciprofloxacin"), ("ERY", "Erythromycin")]

fig = plt.figure(figsize=(7.2, 5.6))
outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.45)
gs_top = outer[0].subgridspec(1, 4, wspace=0.35)
gs_bot = outer[1].subgridspec(1, 2, wspace=0.55)

# ---- A: but over time, small multiples by antibiotic ----
ymax = ab["but"].max() * 1.1
for i, (abx, title) in enumerate(ABX):
    ax = fig.add_subplot(gs_top[0, i])
    for diet, col in [("ND", ND), ("HFD", HFD)]:
        d = ab[(ab.antibiotic == abx) & (ab.diet == diet)]
        m = d.groupby("timepoint")["but"].mean().reindex(TP)
        x = np.arange(3) + (-0.06 if diet == "ND" else 0.06)
        for _, mouse in d.groupby("mouse"):
            mm = mouse.set_index("timepoint")["but"].reindex(TP)
            ax.plot(x, mm.values, color=col, lw=0.6, alpha=0.35, zorder=1)
        ax.plot(x, m.values, color=col, lw=2, zorder=3, solid_capstyle="round")
        ax.scatter(x, m.values, s=22, color=col, edgecolor="white", linewidth=1, zorder=4)
        if abx == "AMX":  # direct labels where ND/HFD are well separated (T0h)
            ax.text(x[0] + 0.12, m.values[0] + (0.012 if diet == "HFD" else -0.018), diet,
                    color=INK, fontsize=7, va="center")
    ax.set_xticks(range(3)); ax.set_xticklabels(TP)
    ax.set_xlim(-0.3, 2.45); ax.set_ylim(0, ymax)
    ax.set_title(title, fontsize=8, color=INK, pad=4)
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    if i == 0:
        ax.set_ylabel("but gene copies per genome", fontsize=8)
    else:
        ax.set_yticklabels([])
fig.text(0.07, 0.94, "A", fontsize=11, fontweight="bold")

# ---- B: SD14 gene abundance vs measured cecal butyrate (group means) ----
ax = fig.add_subplot(gs_bot[0, 0])
for g, r in link.iterrows():
    col = HFD if g.endswith("HFD") else ND
    ax.scatter(r["but"], r["butyrate"], s=46, color=col, edgecolor="white", linewidth=1.5, zorder=3)
    dx, ha = (0.006, "left")
    if g in ("CIP_CTL", "ERY_HFD"):
        dx, ha = (-0.006, "right")
    ax.text(r["but"] + dx, r["butyrate"], g, fontsize=6.5, color=INK, va="center", ha=ha)
rho, p = stats.spearmanr(link["but"], link["butyrate"])
ax.text(0.03, 0.96, f"Spearman ρ = {rho:.2f}, p = {p:.3f}\n(n = 8 group means, SD14)",
        transform=ax.transAxes, fontsize=7, color=MUTED, va="top")
ax.set_xlabel("but copies per genome (shotgun, SD14)")
ax.set_ylabel("Cecal butyrate (μM, GC)")
ax.set_xlim(-0.01, link["but"].max() * 1.25); ax.set_ylim(-3, link["butyrate"].max() * 1.15)
ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
h = [plt.Line2D([], [], marker="o", ls="", color=c, label=l, markersize=6) for l, c in [("ND", ND), ("HFD", HFD)]]
ax.legend(handles=h, frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0.5, 0.02))
fig.text(0.07, 0.47, "B", fontsize=11, fontweight="bold")

# ---- C: genera carrying but (share of total but abundance) ----
ax = fig.add_subplot(gs_bot[0, 1])
gb = gen[gen.gene == "but"].sort_values("share", ascending=False)
named = gb[~gb.genus.isin(["Unbinned", "Unclassified MAG"])]
top = named.head(5)
other = named.iloc[5:]["share"].sum()
rows = list(zip(top.genus, top.share)) + [("Other MAGs", other),
        ("Unbinned contigs", gb.loc[gb.genus == "Unbinned", "share"].sum())]
labels, vals = zip(*rows)
y = np.arange(len(vals))[::-1]
colors = [ND] * len(top) + ["#a9a9a9", "#cfcfcf"]
ax.barh(y, np.array(vals) * 100, color=colors, height=0.62, edgecolor="white", linewidth=1)
for yi, v in zip(y, vals):
    ax.text(v * 100 + 1, yi, f"{v*100:.1f}%", va="center", fontsize=6.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("Share of total but abundance (%)")
ax.set_xlim(0, max(vals) * 100 * 1.2)
ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
for t in ax.get_yticklabels()[:len(top)]:
    t.set_fontstyle("italic")
fig.text(0.55, 0.47, "C", fontsize=11, fontweight="bold")

for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/fig6_butyrate_capacity.{ext}", dpi=300, bbox_inches="tight")
fig.savefig(f"{B}/fig6_tmp.png", dpi=300, bbox_inches="tight")
Image.open(f"{B}/fig6_tmp.png").convert("RGB").save(f"{SUB}/fig6_butyrate_capacity.tiff", dpi=(300, 300), compression="tiff_lzw")
print("saved")
