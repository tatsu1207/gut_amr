#!/usr/bin/env python3
"""Fig. 3 (2026-11-21). A: five genera over time (group mean +/- SEM), evenly spaced timepoints as in Fig. 2.
B: top 15 genera in each group at T0, T3 and R14 (log colour scale).
Samples: experiment ARG, eight groups, >= 5,000 reads; genus % of genus-assigned reads. T6h (untreated groups only,
own MiSeq run) and R7 (failed run) are not shown. Pre-diet: ARG_D0 samples, 16 per diet assignment."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import seaborn as sns
D = "/home/unnot/vault/kjm/miseq/dada2_output/"; OUT = D + "figures/"
asv = pd.read_csv(D + "asv_table.tsv", sep="\t", index_col=0).drop(columns="sequence")
tax = pd.read_csv(D + "taxonomy.tsv", sep="\t", index_col=0)
meta = pd.read_csv(D + "metadata.tsv", sep="\t", index_col=0)
G = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD"]
def pct(ids):
    g = asv[ids].groupby(tax.loc[asv.index, "Genus"]).sum()
    return g.div(g.sum()) * 100
m = meta[(meta.experiment == "ARG") & meta.group.isin(G)]
m = m.loc[m.index.intersection(asv.columns)]; m = m[asv[m.index].sum() >= 5000]
P = pct(m.index)
pre = meta[(meta.experiment == "ARG_D0") & meta.group.isin(["CTL", "HFD"])]
pre = pre.loc[pre.index.intersection(asv.columns)]; PP = pct(pre.index)
TP = ["T0h", "T24h", "T72h", "S72h", "SD10", "SD14"]; LAB = ["T0", "T1", "T3", "R3", "R10", "R14"]
ABX_COL = {"None": "#555555", "AMX": "#E69F00", "CIP": "#0072B2", "ERY": "#009E73"}
def gabx(g): return g.split("_")[0] if "_" in g else "None"
def hfd(g): return g.endswith("HFD")
XPRE = -1.4
plt.rcParams.update({"font.size": 9})

# ---------------- A ----------------
GEN = ["Lactococcus", "Akkermansia", "Faecalibaculum", "Bacteroides", "Escherichia-Shigella"]
fig, axs = plt.subplots(3, 2, figsize=(11, 11)); axs = axs.flatten()
for k, gen in enumerate(GEN):
    ax = axs[k]
    for g in G:
        c = ABX_COL[gabx(g)]; h = hfd(g)
        st = P.loc[gen, m.index[m.group == g]].to_frame("v").join(m.timepoint).groupby("timepoint").v.agg(["mean", "sem"]).reindex(TP)
        ax.errorbar(np.arange(len(TP)) + (0.07 if h else -0.07), st["mean"], yerr=st["sem"], color=c, mec=c,
                    mfc=(c if h else "white"), marker=("^" if h else "o"), ls=("-" if h else "--"), lw=1.5, ms=5,
                    mew=1.2, capsize=2.5, elinewidth=0.9)
    for diet, off in (("ND", -0.12), ("HFD", 0.12)):
        v = PP.loc[gen, pre.index[pre.diet == diet]] if gen in PP.index else pd.Series([0.0])
        ax.errorbar(XPRE + off, v.mean(), yerr=v.sem(), marker=("^" if diet == "HFD" else "o"), ms=6, capsize=3,
                    mfc=("#555555" if diet == "HFD" else "white"), mec="#555555", mew=1.4, color="#555555", ls="none")
    ax.axvline(-0.7, color="gray", ls=":", lw=0.8, alpha=0.7)
    ax.axvline(2.5, color="gray", ls="--", lw=0.8, alpha=0.5)
    ax.text(-0.7, 1.0, "13-wk diet", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color="gray")
    ax.text(1.0, 1.0, "Treatment", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.5, color="gray")
    ax.text(4.0, 1.0, "Recovery", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.5, color="gray")
    ax.set_xticks([XPRE] + list(range(len(TP)))); ax.set_xticklabels(["Pre-diet"] + LAB)
    ax.set_ylim(bottom=0); ax.set_ylabel("% of genus-assigned reads")
    ax.set_title(gen, fontsize=11, fontweight="bold", fontstyle="italic", pad=14)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
axs[5].axis("off")
h = [Line2D([], [], color=c, lw=3, label=("No antibiotic" if a == "None" else a)) for a, c in ABX_COL.items()]
h += [Line2D([], [], color="black", lw=1.5, ls="--", marker="o", mfc="white", label="ND"),
      Line2D([], [], color="black", lw=1.5, ls="-", marker="^", mfc="black", label="HFD"),
      Line2D([], [], color="#555555", ls="none", marker="o", mfc="white", label="Pre-diet (diet assignment)")]
axs[5].legend(handles=h, loc="center", fontsize=10.5, frameon=False, handlelength=3)
plt.tight_layout(h_pad=2.2)
fig.savefig(OUT + "genus_trajectories.png", dpi=300, bbox_inches="tight"); fig.savefig(OUT + "genus_trajectories.pdf", bbox_inches="tight")
plt.close()

# ---------------- B ----------------
show = m[m.timepoint.isin(["T0h", "T72h", "SD14"])]
top = P[show.index].mean(axis=1).sort_values(ascending=False).head(15).index.tolist()
cols, labs = [], []
for tp, tl in (("T0h", "T0"), ("T72h", "T3"), ("SD14", "R14")):
    for g in G:
        cols.append(P.loc[top, show.index[(show.timepoint == tp) & (show.group == g)]].mean(axis=1)); labs.append((tl, g))
H = pd.concat(cols, axis=1); H.columns = pd.MultiIndex.from_tuples(labs)
H.to_csv(D + "diversity/fig3b_genus_by_group_timepoint.tsv", sep="\t", float_format="%.2f")
ann = H.map(lambda v: "0" if v < 0.05 else (f"{v:.0f}" if v >= 10 else f"{v:.1f}"))
fig, ax = plt.subplots(figsize=(17, 6.6))
sns.heatmap(H.clip(lower=0.1), norm=LogNorm(vmin=0.1, vmax=100), cmap="YlOrRd", annot=ann.values, fmt="",
            annot_kws={"fontsize": 8}, linewidths=0.4, linecolor="white", ax=ax,
            cbar_kws={"label": "% of genus-assigned reads (log scale)", "ticks": [0.1, 1, 10, 100]})
ax.collections[0].colorbar.set_ticklabels(["≤0.1", "1", "10", "100"])
ax.set_xticks(np.arange(24) + 0.5); ax.set_xticklabels([g for _, g in labs], rotation=90, fontsize=9)
for t in ax.get_yticklabels(): t.set_fontstyle("italic")
for b in (8, 16): ax.axvline(b, color="black", lw=2)
for j, tl in enumerate(["T0 (baseline)", "T3 (end of treatment)", "R14 (end of recovery)"]):
    ax.text(j * 8 + 4, -0.4, tl, ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout()
fig.savefig(OUT + "genus_heatmap.png", dpi=300, bbox_inches="tight"); fig.savefig(OUT + "genus_heatmap.pdf", bbox_inches="tight")
plt.close()
print("Fig 3 done; heatmap genera:", top)
