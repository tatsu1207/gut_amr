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
# same read threshold as every other 16S analysis (>= 5,000 reads; added 2026-11-21 -- the 33 low-read
# libraries, five of them empty, had been included in Fig. 3 and Supplementary Fig. 3)
arg_samples = arg_samples[asv[arg_samples].sum() >= 5000]
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
    fig, axes = plt.subplots(1, len(available_targets), figsize=(7 * len(available_targets), 5.8))
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

    # colour = antibiotic (Okabe-Ito, colour-blind safe), style = diet: HFD solid line + filled marker,
    # ND dashed line + open marker; ND/HFD nudged apart slightly on x so their error bars do not overlap.
    ABX_COL = {"Control": "#555555", "AMX": "#E69F00", "CIP": "#0072B2", "ERY": "#009E73"}
    DIET_STYLE = {"ND": dict(linestyle="--", mfc="white", dx=-0.12), "HFD": dict(linestyle="-", mfc=None, dx=0.12)}

    # Pre-diet baseline (diet start, before HFD feeding; ARG_D0_* samples with CTL/HFD = diet assignment;
    # 16 per diet; same genus % of genus-assigned reads). Plotted at day XPRE, left of an axis break.
    XPRE = -2.5
    pre_ids = meta[(meta["experiment"] == "ARG_D0") & meta["group"].isin(["CTL", "HFD"])].index.intersection(asv.columns)
    pre = asv[pre_ids].copy(); pre["Genus"] = tax.loc[pre.index, "Genus"]
    pre = pre.groupby("Genus").sum(); pre = pre.div(pre.sum(axis=0), axis=1) * 100
    pre_diet = meta.loc[pre_ids, "diet"]
    PRE_COL = {"ND": "#4477AA", "HFD": "#CC6677"}

    for ax_idx, genus in enumerate(available_targets):
        ax = axes[ax_idx]
        for t_idx, (traj_name, grp_list) in enumerate(traj_groups.items()):
            valid_grps = [g for g in grp_list if g in rel_plot["group"].unique()]
            if not valid_grps:
                continue
            subset = rel_plot[rel_plot["group"].isin(valid_grps)].copy()
            subset = subset.dropna(subset=["timepoint"])
            # R7 (SD7) omitted (2026-11-21): sequenced alone on a failed MiSeq run (KDFCK); see Methods
            subset = subset[subset["timepoint"] != "SD7"]
            if subset.empty:
                continue

            # Mean + SEM per timepoint
            stats = subset.groupby("timepoint", observed=True)[genus].agg(["mean", "sem"]).reset_index()
            stats["day"] = stats["timepoint"].map(
                dict(zip(tp_order, [0, 0.25, 1, 3, 6, 10, 13, 17])))
            stats = stats.dropna(subset=["day"]).sort_values("day")

            abx, diet = traj_name.split("_")
            col, sty = ABX_COL[abx], DIET_STYLE[diet]
            ax.errorbar(stats["day"].astype(float) + sty["dx"], stats["mean"], yerr=stats["sem"],
                       marker="o", markersize=5, capsize=2.5, linewidth=1.6, linestyle=sty["linestyle"],
                       color=col, mec=col, mfc=sty["mfc"] or col, mew=1.3, elinewidth=0.9)

        for diet, off in (("ND", -0.25), ("HFD", 0.25)):
            v = pre.loc[genus, pre_diet[pre_diet == diet].index]
            c = ABX_COL["Control"]
            ax.errorbar(XPRE + off, v.mean(), yerr=v.sem(), marker="o", markersize=6, capsize=3,
                        mfc=DIET_STYLE[diet]["mfc"] or c, mec=c, mew=1.4, color=c, linestyle="none")
        ax.axvline(-1.2, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.text(-1.2, 1.0, "13-wk diet", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=7, color="gray")
        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title(f"{genus} dynamics", fontsize=12, fontweight="bold")
        from matplotlib.lines import Line2D
        h = [Line2D([], [], color=c, lw=2.5, label=a if a != "Control" else "No antibiotic")
             for a, c in ABX_COL.items()]
        h += [Line2D([], [], color="black", lw=1.6, ls="-", marker="o", mfc="black", label="HFD"),
              Line2D([], [], color="black", lw=1.6, ls="--", marker="o", mfc="white", label="ND")]
        legend_handles = h   # one shared legend below the panels (an in-panel legend hid data points)
        ax.axvspan(0, 3, alpha=0.1, color="red", label="Treatment" if ax_idx == 0 else None)
        ax.set_xticks([XPRE, 0, 1, 3, 6, 13, 17])
        ax.set_xticklabels(["Pre-diet", "T0h", "T24h", "T72h", "S72h", "SD10", "SD14"])  # T6h and R7 not shown

    fig.suptitle("Key Genus Dynamics During Antibiotic Treatment and Recovery", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=6, fontsize=10,
               frameon=False, handlelength=3)
    fig.savefig(OUT / "genus_trajectories.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "genus_trajectories.pdf", bbox_inches="tight")
    print(f"Saved: {OUT / 'genus_trajectories.png'}")
    plt.close()

# ============================================================
# Figure 3: Overall composition heatmap - top genera across groups
# ============================================================
# R7 (SD7) excluded (2026-11-21): sequenced alone on a failed MiSeq run (KDFCK); top 15 genera ranked
# and averaged over the remaining timepoints.
_ok = meta_arg.loc[rel_abund.columns, "timepoint"] != "SD7"
_rel7 = rel_abund.loc[:, _ok.values]
top_hm = _rel7.mean(axis=1).sort_values(ascending=False).head(15).index.tolist()
print("heatmap genera changed vs all-timepoint top 15:", sorted(set(top_hm) ^ set(top_genera)))
mean_by_group = _rel7.loc[top_hm].T.join(meta_arg[["group"]]).groupby("group")[top_hm].mean()
mean_by_group = mean_by_group.reindex([g for g in groups_to_plot if g in mean_by_group.index])

fig, ax = plt.subplots(figsize=(12, 7.5))
sns.heatmap(mean_by_group.T, cmap="YlOrRd", annot=True, fmt=".1f",
            ax=ax, cbar_kws={"label": "Mean relative abundance (%)"})
ax.set_title("Mean genus relative abundance by group (all timepoints except R7 pooled)", fontsize=12)
ax.set_ylabel("Genus")
ax.set_xlabel("Treatment Group")
plt.tight_layout()
fig.savefig(OUT / "genus_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "genus_heatmap.pdf", bbox_inches="tight")
print(f"Saved: {OUT / 'genus_heatmap.png'}")
plt.close()

print("All taxonomy plots done.")
