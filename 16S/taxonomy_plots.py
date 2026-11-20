#!/usr/bin/env python3
"""Taxonomy composition plots for KJM project.
- Stacked bar plots by group and timepoint (genus level)
- Lactococcus and Akkermansia relative abundance trajectories over time
"""

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
tp_order = ["T0h", "T6h", "T24h", "T72h", "S72h", "SD7", "SD10", "SD14"]
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

# ============================================================
# Figure 2: Lactococcus and Akkermansia trajectories
# ============================================================
target_genera = ["Lactococcus", "Akkermansia"]
available_targets = [g for g in target_genera if g in rel_abund.index]

if available_targets:
    fig, axes = plt.subplots(1, len(available_targets), figsize=(7 * len(available_targets), 5))
    if len(available_targets) == 1:
        axes = [axes]

    # Group by antibiotic + diet for cleaner trajectory lines
    rel_plot["abx_diet"] = rel_plot["antibiotic"].fillna("None") + "_" + rel_plot["diet"].fillna("None")

    # Define trajectory groups: antibiotic-based
    traj_groups = {
        "AMX_ND": ["AMX_CTL"],
        "AMX_HFD": ["AMX_HFD"],
        "CIP_ND": ["CIP_CTL"],
        "CIP_HFD": ["CIP_HFD"],
        "ERY_ND": ["ERY_CTL"],
        "ERY_HFD": ["ERY_HFD"],
        "Control_ND": ["CTL"],
        "Control_HFD": ["HFD"],
    }

    traj_colors = sns.color_palette("Set2", len(traj_groups))

    for ax_idx, genus in enumerate(available_targets):
        ax = axes[ax_idx]
        for t_idx, (traj_name, grp_list) in enumerate(traj_groups.items()):
            valid_grps = [g for g in grp_list if g in rel_plot["group"].unique()]
            if not valid_grps:
                continue
            subset = rel_plot[rel_plot["group"].isin(valid_grps)].copy()
            subset = subset.dropna(subset=["timepoint"])
            if subset.empty:
                continue

            # Mean + SEM per timepoint
            stats = subset.groupby("timepoint", observed=True)[genus].agg(["mean", "sem"]).reset_index()
            stats["day"] = stats["timepoint"].map(
                dict(zip(tp_order, [0, 0.25, 1, 3, 6, 10, 13, 17])))
            stats = stats.dropna(subset=["day"]).sort_values("day")

            ax.errorbar(stats["day"], stats["mean"], yerr=stats["sem"],
                       marker="o", markersize=4, capsize=3, linewidth=1.5,
                       color=traj_colors[t_idx], label=traj_name)

        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title(f"{genus} dynamics", fontsize=12, fontweight="bold")
        ax.legend(fontsize=7, loc="best")
        ax.axvspan(0, 3, alpha=0.1, color="red", label="Treatment" if ax_idx == 0 else None)
        ax.set_xticks([0, 1, 3, 6, 10, 13, 17])
        ax.set_xticklabels(tp_order[:1] + tp_order[2:])  # skip T6h for readability

    fig.suptitle("Key Genus Dynamics During Antibiotic Treatment and Recovery", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / "genus_trajectories.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "genus_trajectories.pdf", bbox_inches="tight")
    print(f"Saved: {OUT / 'genus_trajectories.png'}")
    plt.close()

# ============================================================
# Figure 3: Overall composition heatmap - top genera across groups
# ============================================================
mean_by_group = rel_plot.groupby("group")[top_genera].mean()
mean_by_group = mean_by_group.reindex([g for g in groups_to_plot if g in mean_by_group.index])

fig, ax = plt.subplots(figsize=(12, 6))
sns.heatmap(mean_by_group.T, cmap="YlOrRd", annot=True, fmt=".1f",
            ax=ax, cbar_kws={"label": "Mean relative abundance (%)"})
ax.set_title("Mean genus relative abundance by group (all timepoints pooled)", fontsize=12)
ax.set_ylabel("Genus")
ax.set_xlabel("Treatment Group")
plt.tight_layout()
fig.savefig(OUT / "genus_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "genus_heatmap.pdf", bbox_inches="tight")
print(f"Saved: {OUT / 'genus_heatmap.png'}")
plt.close()

print("All taxonomy plots done.")
