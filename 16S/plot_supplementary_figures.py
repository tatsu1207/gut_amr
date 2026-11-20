#!/usr/bin/env python3
"""Supplementary figures S1 (rarefaction curves) and S2 (cecal SCFA by group).

S1 justifies the rarefaction depth of 5,159 reads used for alpha diversity.
S2 shows per-mouse SCFA values with group-level tests, complementing the
box-plot summary in Fig. 1D-F.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

OUT = "/home/unnot/vault/kjm/miseq/dada2_output/figures"
ASV = "/home/unnot/vault/kjm/miseq/dada2_output/asv_table_arg.tsv"
RAW = "/home/unnot/vault/kjm/Raw_data_bio-prameter.xlsx"
DEPTH = 5159
GROUPS = ["CTL", "HFD", "CIP_CTL", "CIP_HFD", "AMX_CTL", "AMX_HFD", "ERY_CTL", "ERY_HFD"]
ND = {"CTL", "CIP_CTL", "AMX_CTL", "ERY_CTL"}
BLUE, RED = "#3070B8", "#CC6677"
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})

# ---------------- S1: rarefaction curves ----------------
d = pd.read_csv(ASV, sep="\t", index_col=0)
if d.shape[0] > d.shape[1]:
    d = d.T                      # samples as rows
counts = d.values.astype(int)
depths = d.sum(axis=1).values
keep = depths >= 5000
counts, depths = counts[keep], depths[keep]

def rarefy_richness(row, n, reps=1):
    """Expected observed ASVs at depth n (analytic Hurlbert estimator)."""
    row = row[row > 0]
    N = row.sum()
    if N < n:
        return np.nan
    # E[S] = sum_i 1 - C(N-xi, n)/C(N, n), computed in log space
    from scipy.special import gammaln
    lg = gammaln
    term = np.exp(lg(N - row + 1) - lg(N - row - n + 1) - (lg(N + 1) - lg(N - n + 1)))
    return float(np.sum(1 - term))

xs = [100, 250, 500, 1000, 2000, 3000, 4000, DEPTH, 7500, 10000, 15000, 20000, 30000]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
ax = axes[0]
rng = np.random.default_rng(42)
idx = rng.choice(len(counts), size=min(150, len(counts)), replace=False)  # 150 curves for legibility
for i in idx:
    ys = [rarefy_richness(counts[i], x) for x in xs if x <= depths[i]]
    ax.plot(xs[:len(ys)], ys, color="#9fb6cc", lw=0.5, alpha=0.5)
ax.axvline(DEPTH, color=RED, lw=1.5, ls="--")
ax.text(DEPTH * 1.15, ax.get_ylim()[1] * 0.08, f"rarefaction depth\n{DEPTH:,} reads",
        color=RED, fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("Sequencing depth (reads)")
ax.set_ylabel("Observed ASVs"); ax.set_title("a  Rarefaction curves", loc="left", fontweight="bold")

ax = axes[1]
ax.hist(depths, bins=50, color="#9fb6cc", edgecolor="white", linewidth=0.4)
ax.axvline(DEPTH, color=RED, lw=1.5, ls="--")
ax.set_xscale("log"); ax.set_xlabel("Library size (reads)"); ax.set_ylabel("Samples")
ax.set_title("b  Library size distribution", loc="left", fontweight="bold")
ax.text(0.97, 0.95, f"n = {len(depths)} samples\nmin = {depths.min():,}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/figS1_rarefaction.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("S1 done:", len(depths), "samples, min depth", depths.min())

# ---------------- S2: cecal SCFA with individual mice ----------------
sc = pd.read_excel(RAW, sheet_name="SCFA")
panels = [("Acetate", "group", "value(μM)"), ("Propionate", "group.1", "value(μM).1"),
          ("Butyrate", "group.2", "value(μM).2")]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
for ax, (name, gcol, vcol) in zip(axes, panels):
    df = sc[[gcol, vcol]].dropna()
    df["g"] = df[gcol].str.replace(r"^\d+_", "", regex=True)
    groups = [df.loc[df.g == g, vcol].values for g in GROUPS]
    H, p = stats.kruskal(*groups)
    for i, g in enumerate(GROUPS):
        v = df.loc[df.g == g, vcol].values
        col = BLUE if g in ND else RED
        ax.scatter(np.full(len(v), i) + rng.uniform(-.12, .12, len(v)), v,
                   s=16, color=col, alpha=.75, edgecolor="white", linewidth=.4, zorder=3)
        ax.hlines(v.mean(), i - .28, i + .28, color="#333", lw=1.6, zorder=4)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GROUPS, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Concentration (μM)" if name == "Acetate" else "")
    ax.set_title(f"{'abc'[panels.index((name,gcol,vcol))]}  {name}", loc="left", fontweight="bold")
    ax.text(0.97, 0.95, f"Kruskal-Wallis\np = {p:.1e}" if p < 0.001 else f"Kruskal-Wallis\np = {p:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#555")
    ax.grid(axis="y", color="#e6e6e6", lw=.5); ax.set_axisbelow(True)
h = [plt.Line2D([], [], marker="o", ls="", color=c, label=l, markersize=6) for l, c in [("ND", BLUE), ("HFD", RED)]]
fig.legend(handles=h, frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 1.04), ncol=2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/figS2_scfa.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("S2 done")
