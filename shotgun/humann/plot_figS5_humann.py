#!/usr/bin/env python3
"""Supplementary Fig. 5 — community-wide pathway profiling (HUMAnN 3.9) and why butyrate is
reported from the assemblies instead.

(A) Number of MetaCyc pathways detected per sample at T0, T3 and R14, by group.
(B) Read-level EC 2.8.3.8 (CPM) against the assembly-based terminal butyrate measure
    (copies per genome equivalent, Fig. 6A), one point per sample.

Colour scheme is the one shared by Figs. 1, 2, 3, 4 and 6: colour = antibiotic (Okabe-Ito),
ND = dashed line / open circle, HFD = solid line / filled triangle.

Input: tables/humann_butyrate_reconcile.tsv (humann_butyrate_reconcile.py).
Output: ../figures/figS5_humann_pathways.{png,pdf} and the submission TIFF.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image
from scipy.stats import spearmanr

H = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(H, "tables")
FIG = os.path.abspath(os.path.join(H, "../figures"))
SUB = os.path.expanduser("~/vault/KJM_info/manuscript/figures")
os.makedirs(FIG, exist_ok=True)

ABX_COL = {"None": "#555555", "AMX": "#E69F00", "CIP": "#0072B2", "ERY": "#009E73"}
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})


def gabx(g):
    return g.split("_")[0] if "_" in g else "None"


def is_hfd(g):
    return g == "HFD" or g.endswith("_HFD")


m = pd.read_csv(os.path.join(T, "humann_butyrate_reconcile.tsv"), sep="\t", index_col=0)
TP = ["T0", "T3", "R14"]
GROUPS = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD"]

fig = plt.figure(figsize=(7.2, 9.0))
gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.95, 0.95], hspace=0.46)
axA = fig.add_subplot(gs[0])
axD = fig.add_subplot(gs[1])   # (B) depth vs richness at T3
axB = fig.add_subplot(gs[2])   # (C) read-level vs assembly butyrate

# ---- (A) pathway richness over the three shotgun timepoints -------------------------------
x = np.arange(len(TP))
rng = np.random.default_rng(42)
for g in GROUPS:
    c = ABX_COL[gabx(g)]
    hfd = is_hfd(g)
    sub = m[m["group"] == g]
    mu = [sub[sub.timepoint == t]["n_pathways"].mean() for t in TP]
    sd = [sub[sub.timepoint == t]["n_pathways"].std() for t in TP]
    off = 0.055 * (1 if hfd else -1)
    axA.errorbar(x + off, mu, yerr=sd, color=c, lw=1.3, capsize=2.2, elinewidth=0.7,
                 ls=("-" if hfd else "--"), marker=("^" if hfd else "o"), ms=5,
                 mfc=(c if hfd else "white"), mec=c, zorder=4)
    for k, t in enumerate(TP):
        v = sub[sub.timepoint == t]["n_pathways"].values
        axA.scatter(np.full(len(v), x[k] + off) + rng.uniform(-0.03, 0.03, len(v)), v,
                    s=7, marker=("^" if hfd else "o"), facecolor=(c if hfd else "white"),
                    edgecolor=c, linewidth=0.5, alpha=0.75, zorder=3)
axA.set_xticks(x)
axA.set_xticklabels(TP)
axA.set_xlim(-0.35, len(TP) - 0.65)
axA.set_ylabel("MetaCyc pathways detected per sample")
axA.set_xlabel("Timepoint")
axA.yaxis.grid(True, color=GRID, lw=0.6)
axA.set_axisbelow(True)
axA.set_title("A", loc="left", fontweight="bold", fontsize=11, pad=8)
axA.annotate("T3: antibiotic $p$ = 0.001, diet $p$ = 0.004\n(two-way ANOVA, interaction $p$ = 0.09)",
             xy=(0.985, 0.04), xycoords="axes fraction", ha="right", va="bottom",
             fontsize=7, color=MUTED)
handles = [Line2D([], [], color=c, lw=1.3, ls="-", marker="s", ms=0, label=lab)
           for lab, c in [("No antibiotic", ABX_COL["None"]), ("Amoxicillin", ABX_COL["AMX"]),
                          ("Ciprofloxacin", ABX_COL["CIP"]), ("Erythromycin", ABX_COL["ERY"])]]
handles += [Line2D([], [], color=MUTED, lw=1.1, ls="--", marker="o", ms=5, mfc="white",
                   mec=MUTED, label="ND"),
            Line2D([], [], color=MUTED, lw=1.1, ls="-", marker="^", ms=5, mfc=MUTED,
                   mec=MUTED, label="HFD")]
axA.legend(handles=handles, frameon=False, fontsize=7, ncol=3, loc="upper center",
           bbox_to_anchor=(0.5, 1.30), handlelength=2.0, columnspacing=1.4)

# ---- (B) sequencing depth against richness at T3 ------------------------------------------
t3 = m[m["timepoint"] == "T3"]
rho3, pv3 = spearmanr(t3["genome_equivalents_depth"], t3["n_pathways"])
for g in GROUPS:
    c = ABX_COL[gabx(g)]
    hfd = is_hfd(g)
    sub = t3[t3["group"] == g]
    axD.scatter(sub["genome_equivalents_depth"], sub["n_pathways"], s=34,
                marker=("^" if hfd else "o"), facecolor=(c if hfd else "white"),
                edgecolor=c, linewidth=0.9, alpha=0.95, zorder=3)
amx = t3[t3["group"].isin(["AMX_CTL", "AMX_HFD"])]
if len(amx):
    for grp, ls in (("AMX_CTL", "--"), ("AMX_HFD", "-")):
        mu = t3[t3["group"] == grp]["genome_equivalents_depth"].mean()
        axD.axvline(mu, color=ABX_COL["AMX"], ls=ls, lw=1.0, alpha=0.85, zorder=1)
    axD.annotate("amoxicillin group means (lines):\nthe two diets are matched for depth\n"
                 "(1,062 vs 1,001), so 283 vs 201\npathways is not a depth artefact",
                 xy=(0.015, 0.95), xycoords="axes fraction", ha="left", va="top", fontsize=7,
                 color=INK)
axD.set_xlabel("Effective sequencing depth at T3 (genome equivalents)")
axD.set_ylabel("MetaCyc pathways detected")
axD.grid(True, color=GRID, lw=0.6)
axD.set_axisbelow(True)
axD.set_title("B", loc="left", fontweight="bold", fontsize=11, pad=8)
axD.annotate(f"T3: Spearman $\\rho$ = {rho3:.2f}, $p$ = {pv3:.1g} ($n$ = {len(t3)})\n"
             f"antibiotic term after adjusting for depth: $p$ = 0.15; diet term: $p$ = 0.003",
             xy=(0.985, 0.06), xycoords="axes fraction", ha="right", va="bottom", fontsize=7,
             color=MUTED)

# ---- (C) read-level EC 2.8.3.8 against the assembly measure -------------------------------
rho, pv = spearmanr(m["humann_but_cpm"], m["terminal_but_buk"])
for g in GROUPS:
    c = ABX_COL[gabx(g)]
    hfd = is_hfd(g)
    sub = m[m["group"] == g]
    axB.scatter(sub["terminal_but_buk"], sub["humann_but_cpm"], s=26,
                marker=("^" if hfd else "o"), facecolor=(c if hfd else "white"),
                edgecolor=c, linewidth=0.8, alpha=0.9, zorder=3)
ac7 = "ARG_AC7_T72h"
if ac7 in m.index:
    r = m.loc[ac7]
    axB.annotate("AMX_CTL mouse 7, T3\nno $but$ gene in the assembly,\n203 CPM of EC 2.8.3.8",
                 xy=(r["terminal_but_buk"], r["humann_but_cpm"]), xytext=(0.30, 0.93),
                 textcoords="axes fraction", fontsize=7, color=INK, va="top",
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
axB.set_yscale("symlog", linthresh=1, linscale=0.5)
axB.set_yticks([0, 1, 10, 100])
axB.set_yticklabels(["0", "1", "10", "100"])
axB.set_ylim(-0.25, 400)
axB.set_xlabel("Assembly: terminal butyrate genes ($but$ + $buk$, copies per genome equivalent)")
axB.set_ylabel("HUMAnN: EC 2.8.3.8 (CPM, symlog)")
axB.grid(True, color=GRID, lw=0.6)
axB.set_axisbelow(True)
axB.set_title("C", loc="left", fontweight="bold", fontsize=11, pad=8)
axB.annotate(f"Spearman $\\rho$ = {rho:.2f}, $p$ = {pv:.2f} ($n$ = {len(m)})",
             xy=(0.985, 0.94), xycoords="axes fraction", ha="right", va="top", fontsize=7.5,
             color=INK)

for ext in ("png", "pdf"):
    fig.savefig(f"{FIG}/figS5_humann_pathways.{ext}", dpi=300, bbox_inches="tight")
tmp = f"{T}/figS5_tmp.png"
fig.savefig(tmp, dpi=300, bbox_inches="tight")
if os.path.isdir(SUB):
    Image.open(tmp).convert("RGB").save(f"{SUB}/figS5_humann_pathways.tiff",
                                        dpi=(300, 300), compression="tiff_lzw")
os.remove(tmp)
print(f"wrote {FIG}/figS5_humann_pathways.png  (rho={rho:.2f}, p={pv:.2f})")
