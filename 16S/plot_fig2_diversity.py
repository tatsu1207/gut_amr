#!/usr/bin/env python3
"""Generate Fig 2: Alpha and beta diversity panels."""

import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.expanduser("~/vault/kjm/miseq"))
import tplabels; tplabels.install()  # show T0/T1/T3/R3...R14 instead of the data codes
import matplotlib as mpl
import numpy as np

mpl.rcParams.update({
    'font.size': 9, 'font.family': 'sans-serif',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 300, 'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

DIVDIR = "/home/unnot/vault/kjm/miseq/dada2_output/diversity"
OUTDIR = "/home/unnot/vault/kjm/miseq/dada2_output/figures"

alpha = pd.read_csv(f"{DIVDIR}/alpha_diversity.csv")
# R7 (SD7) is left out of every panel (2026-11-21): it was sequenced alone on a failed MiSeq run (KDFCK),
# and untreated mice lost two thirds of their ASVs at R7 only. PCoA recomputed without R7 by run_diversity.R.
pcoa = pd.read_csv(f"{DIVDIR}/pcoa_bray_curtis_noR7.csv")
eigen = pd.read_csv(f"{DIVDIR}/pcoa_eigenvalues_noR7.csv")
bc_base = pd.read_csv(f"{DIVDIR}/bc_from_baseline.csv")
# Pre-diet baseline (diet start, before HFD feeding; ARG_D0_* samples, 16 per diet assignment),
# rarefied to the same 5,159 reads by miseq/d0_baseline.R. Shown left of T0 in panels A and B.
_d0 = pd.read_csv(f"{DIVDIR}/d0_baseline_alpha.tsv", sep="\t")
prediet = _d0[_d0["stage"] == "D0_prediet"]
XPRE = -1.4   # x position of the pre-diet points (T0 = 0)

def add_prediet(ax, col):
    """Pre-diet ND/HFD means +/- SEM, open markers, with an axis-break mark before T0."""
    for diet, c, off in (("ND", "#555555", -0.12), ("HFD", "#555555", 0.12)):
        v = prediet[prediet["diet"] == diet][col]
        ax.errorbar(XPRE + off, v.mean(), yerr=v.sem(), marker=("o" if diet == "ND" else "^"), mfc=("white" if diet == "ND" else c), mec=c, mew=1.4,
                    color=c, capsize=3, markersize=6, linestyle="none", label=f"Pre-diet, {diet}")
    ax.axvline(x=-0.7, color="gray", linestyle=":", alpha=0.6, linewidth=0.8)
    ax.text(-0.7, 1.0, "13-wk diet", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=7, color="gray")

# Filter to main 8 groups
GROUPS = ['CTL', 'HFD', 'AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD']
TP_ORDER = ['T0h', 'T24h', 'T72h', 'S72h', 'SD10', 'SD14']   # R7 omitted (failed run); T6h omitted (untreated groups only, own run)

alpha = alpha[alpha['group'].isin(GROUPS)].copy()
pcoa = pcoa[pcoa['group'].isin(GROUPS)].copy()
bc_base = bc_base[bc_base['group'].isin(GROUPS)].copy()

# Style (same scheme as Fig. 3A): colour = antibiotic (Okabe-Ito); A, B, D: HFD solid line + filled
# marker, ND dashed line + open marker; C: circles = ND, triangles = HFD.
ABX_COL = {"None": "#555555", "AMX": "#E69F00", "CIP": "#0072B2", "ERY": "#009E73"}
def gabx(g): return g.split("_")[0] if "_" in g else "None"
def gdiet(g): return "HFD" if g.endswith("HFD") else "ND"
def line_style(g):
    d = gdiet(g); c = ABX_COL[gabx(g)]
    return dict(color=c, mec=c, mfc=(c if d == "HFD" else "white"), marker=("^" if d == "HFD" else "o"), mew=1.2,
                linestyle=("-" if d == "HFD" else "--"), linewidth=1.5, capsize=2.5, markersize=5)
DX = {"ND": -0.07, "HFD": 0.07}   # nudge ND/HFD apart so error bars do not overlap

# ═══════════════════════════════════════════════════════
# Fig 2: 4-panel figure
# ═══════════════════════════════════════════════════════
fig = plt.figure(figsize=(11, 17))  # portrait: A B / C (T0, T3, R14) / D
gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 0.95, 1.0], hspace=0.38)
g1 = gs[0].subgridspec(1, 2, wspace=0.22); g2 = gs[1].subgridspec(1, 3, wspace=0.28)
axes = np.empty((2, 2), dtype=object)
axes[0, 0] = fig.add_subplot(g1[0, 0]); axes[0, 1] = fig.add_subplot(g1[0, 1])
axes_c = [fig.add_subplot(g2[0, k]) for k in range(3)]
axes[1, 1] = fig.add_subplot(gs[2])

# ── Panel A: Shannon diversity over time ──
ax = axes[0, 0]
alpha['timepoint'] = pd.Categorical(alpha['timepoint'], categories=TP_ORDER, ordered=True)

for group in GROUPS:
    sub = alpha[alpha['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['shannon'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) + DX[gdiet(group)] for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'], label=group, **line_style(group))

add_prediet(ax, "shannon")
ax.set_xticks([XPRE] + list(range(len(TP_ORDER))))
ax.set_xticklabels(["Pre-diet"] + TP_ORDER, rotation=45, ha='right')
ax.set_ylabel("Shannon index")
ax.set_title("A) Shannon diversity", fontweight='bold')
ax.axvline(x=2.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.text(1.0, ax.get_ylim()[1]*0.98, 'Treatment', ha='center', fontsize=8, color='gray')
ax.text(4.0, ax.get_ylim()[1]*0.98, 'Recovery', ha='center', fontsize=8, color='gray')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel B: Observed ASVs (richness) over time ──
ax = axes[0, 1]
for group in GROUPS:
    sub = alpha[alpha['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['observed_ASVs'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) + DX[gdiet(group)] for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'], label=group, **line_style(group))

add_prediet(ax, "richness")
ax.set_xticks([XPRE] + list(range(len(TP_ORDER))))
ax.set_xticklabels(["Pre-diet"] + TP_ORDER, rotation=45, ha='right')
ax.set_ylabel("Observed ASVs")
ax.set_title("B) Richness (observed ASVs)", fontweight='bold')
ax.axvline(x=2.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel C: PCoA per timepoint (each computed on that timepoint only; permanova_factorial.R) ──
from matplotlib.patches import Ellipse
pt = pd.read_csv(f"{DIVDIR}/pcoa_per_timepoint.tsv", sep="\t")
fac = pd.read_csv(f"{DIVDIR}/permanova_factorial_by_timepoint.tsv", sep="\t").set_index("timepoint")
def star(p): return "***" if p <= 0.001 else ("**" if p <= 0.01 else ("*" if p <= 0.05 else " n.s."))
for k, tp in enumerate(["T0h", "T72h", "SD14"]):
    ax = axes_c[k]; d = pt[pt.timepoint == tp]
    for group in GROUPS:
        sub = d[d.group == group]; c = ABX_COL[gabx(group)]; h = gdiet(group) == "HFD"
        ax.scatter(sub.PCoA1, sub.PCoA2, s=16, marker=("^" if h else "o"), facecolor=(c if h else "white"),
                   edgecolor=c, linewidth=0.8, alpha=0.6, zorder=2)
        if len(sub) >= 3:   # 1-SD ellipse of the group's points
            cov = np.cov(sub.PCoA1, sub.PCoA2); vals, vecs = np.linalg.eigh(cov)
            ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
            ax.add_patch(Ellipse((sub.PCoA1.mean(), sub.PCoA2.mean()), 2 * np.sqrt(vals[1]), 2 * np.sqrt(vals[0]),
                                 angle=ang, fill=False, edgecolor=c, lw=1.1, ls=("-" if h else "--"), zorder=1))
        ax.scatter(sub.PCoA1.mean(), sub.PCoA2.mean(), s=70, marker=("^" if h else "o"),
                   facecolor=(c if h else "white"), edgecolor=c, linewidth=1.8, zorder=4)
    r = fac.loc[tp]
    ax.text(0.02, 0.02, f"Antibiotic R² {r.R2_antibiotic:.2f}{star(r.p_antibiotic)}\nDiet R² {r.R2_diet:.2f}{star(r.p_diet)}\n"
            f"Abx × diet R² {r.R2_interaction:.2f}{star(r.p_interaction)}", transform=ax.transAxes, fontsize=7.5,
            va="bottom", ha="left", color="#333333", bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
    ax.set_xlabel(f"PCoA1 ({d.pct1.iloc[0]:.1f}%)"); ax.set_ylabel(f"PCoA2 ({d.pct2.iloc[0]:.1f}%)")
    ax.set_title({"T0h": "T0", "T72h": "T3", "SD14": "R14"}[tp], fontsize=10, fontweight="bold")
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
fig.text(axes_c[0].get_position().x0 - 0.05, axes_c[0].get_position().y1 + 0.022,
         "C) PCoA (Bray-Curtis) at T0, T3 and R14", fontsize=11, fontweight="bold")

# ── Panel D: Bray-Curtis dissimilarity from baseline ──
ax = axes[1, 1]
bc_base['timepoint'] = pd.Categorical(bc_base['timepoint'], categories=TP_ORDER, ordered=True)

for group in GROUPS:
    sub = bc_base[bc_base['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['bc_from_baseline'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) + DX[gdiet(group)] for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'], label=group, **line_style(group))

ax.set_xticks(range(1, len(TP_ORDER)))
ax.set_xticklabels(TP_ORDER[1:], rotation=45, ha='right')
ax.set_xlim(0.6, len(TP_ORDER) - 0.6)
ax.set_ylabel("Bray-Curtis dissimilarity from baseline (T0)")
ax.set_title("D) Community displacement from baseline", fontweight='bold')
ax.axvline(x=2.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

from matplotlib.lines import Line2D
_h = [Line2D([], [], color=c, lw=3, label=("No antibiotic" if a == "None" else a)) for a, c in ABX_COL.items()]
_h += [Line2D([], [], color="black", lw=1.5, ls="--", marker="o", mfc="white", label="ND"),
       Line2D([], [], color="black", lw=1.5, ls="-", marker="^", mfc="black", label="HFD")]
fig.legend(handles=_h, loc="upper center", bbox_to_anchor=(0.5, 0.075), ncol=6, fontsize=9, frameon=False,
           handlelength=2.6, columnspacing=1.2)
plt.savefig(f"{OUTDIR}/fig2_diversity.pdf")
plt.savefig(f"{OUTDIR}/fig2_diversity.png")
plt.close()
print("Fig 2 done")
