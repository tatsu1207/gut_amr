#!/usr/bin/env python3
"""Supplementary Fig. 3: genus composition stacked bars by group and timepoint (Fig. 3 is fig3_taxonomy.py)."""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.expanduser("~/vault/kjm/miseq"))
import tplabels; tplabels.install()  # show T0/T1/T3/R3...R14 instead of the data codes
import seaborn as sns
from pathlib import Path

OUT = Path("/home/unnot/vault/kjm/miseq/dada2_output/figures")
OUT.mkdir(exist_ok=True)

# Load data
asv = pd.read_csv("/home/unnot/vault/kjm/miseq/dada2_output/asv_table.tsv", sep="\t", index_col=0)
tax = pd.read_csv("/home/unnot/vault/kjm/miseq/dada2_output/taxonomy.tsv", sep="\t", index_col=0)
meta = pd.read_csv("/home/unnot/vault/kjm/miseq/dada2_output/metadata.tsv", sep="\t", index_col=0)

# Drop sequence column from ASV table
asv = asv.drop(columns=["sequence"], errors="ignore")

# Filter to ARG experiment samples only (have proper metadata)
arg_samples = meta[meta["experiment"] == "ARG"].index  # main experiment only, as in the text (ARG_D0 day-0 batch excluded)
arg_samples = arg_samples.intersection(asv.columns)
# same read threshold as every other 16S analysis (>= 5,000 reads; added 2026-11-21 -- the 33 low-read
# libraries, five of them empty, had been included in Fig. 3 and Supplementary Fig. 3)
arg_samples = arg_samples[asv[arg_samples].sum() >= 5000]
# T6h (untreated groups only, own MiSeq run) and R7 (failed MiSeq run) left out (2026-11-21), as in Figs. 2 and 3
arg_samples = arg_samples[~meta.loc[arg_samples, "timepoint"].isin(["T6h", "SD7"])]
meta_arg = meta.loc[arg_samples]
asv_arg = asv[arg_samples]

# Merge taxonomy at genus level
asv_arg_tax = asv_arg.copy()
asv_arg_tax["Genus"] = tax.loc[asv_arg_tax.index, "Genus"]
genus_counts = asv_arg_tax.groupby("Genus").sum()

# Convert to relative abundance
rel_abund = genus_counts.div(genus_counts.sum(axis=0), axis=1) * 100

# Identify top genera (top 15 by mean relative abundance)
top_genera = rel_abund.mean(axis=1).sort_values(ascending=False).head(15).index.tolist()
# Collapse rest into "Other"
rel_plot = rel_abund.loc[top_genera].T.copy()
rel_plot["Other"] = 100 - rel_plot.sum(axis=1)
rel_plot = rel_plot.clip(lower=0)

# Add metadata
rel_plot = rel_plot.join(meta_arg[["group", "timepoint", "day", "diet", "antibiotic", "phase"]])

# Define timepoint order
tp_order = ["T0h", "T24h", "T72h", "S72h", "SD10", "SD14"]
rel_plot["timepoint"] = pd.Categorical(rel_plot["timepoint"], categories=tp_order, ordered=True)

# Color palette for genera
n_colors = len(top_genera) + 1
colors = sns.color_palette("tab20", n_colors)

# ============================================================
# Figure 1: Stacked bar plot - mean composition by group x timepoint
# ============================================================
groups_to_plot = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD",
                  "ERY_CTL", "ERY_HFD"]  # the AMX_Early/CIP_Early cohorts are not part of this study
groups_to_plot = [g for g in groups_to_plot if g in rel_plot["group"].unique()]

fig, axes = plt.subplots(2, len(groups_to_plot) // 2 + len(groups_to_plot) % 2,
                         figsize=(24, 10), sharey=True)
axes = axes.flatten()

genera_cols = top_genera + ["Other"]

for idx, grp in enumerate(groups_to_plot):
    ax = axes[idx]
    grp_data = rel_plot[rel_plot["group"] == grp].copy()
    grp_data = grp_data.dropna(subset=["timepoint"])
    if grp_data.empty:
        ax.set_title(grp)
        continue
    mean_comp = grp_data.groupby("timepoint", observed=False)[genera_cols].mean()
    mean_comp = mean_comp.loc[mean_comp.index.isin(grp_data["timepoint"].unique())]

    bottom = np.zeros(len(mean_comp))
    for i, genus in enumerate(genera_cols):
        vals = mean_comp[genus].values
        ax.bar(range(len(mean_comp)), vals, bottom=bottom,
               color=colors[i], label=genus if idx == 0 else None, width=0.8)
        bottom += vals

    ax.set_xticks(range(len(mean_comp)))
    ax.set_xticklabels(mean_comp.index, rotation=45, ha="right", fontsize=8)
    ax.set_title(grp, fontsize=10, fontweight="bold")
    ax.set_ylim(0, 100)
    if idx == 0:
        ax.set_ylabel("Relative abundance (%)")

# Hide unused axes
for idx in range(len(groups_to_plot), len(axes)):
    axes[idx].set_visible(False)

fig.legend(genera_cols, loc="center left", bbox_to_anchor=(1.0, 0.5),
           fontsize=12, title="Genus", title_fontsize=13)
fig.suptitle("Gut Microbiota Composition by Treatment Group Over Time", fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(OUT / "taxonomy_stacked_bars.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "taxonomy_stacked_bars.pdf", bbox_inches="tight")
print(f"Saved: {OUT / 'taxonomy_stacked_bars.png'}")
plt.close()

# Fig. 3A (genus trajectories) and Fig. 3B (genus heatmap by timepoint) are drawn by fig3_taxonomy.py
# since 2026-11-21; this script now makes only the stacked bars (Supplementary Fig. 3).
print("All taxonomy plots done.")
