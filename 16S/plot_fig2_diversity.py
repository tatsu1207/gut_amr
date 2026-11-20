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
pcoa = pd.read_csv(f"{DIVDIR}/pcoa_bray_curtis.csv")
eigen = pd.read_csv(f"{DIVDIR}/pcoa_eigenvalues.csv")
bc_base = pd.read_csv(f"{DIVDIR}/bc_from_baseline.csv")

# Filter to main 8 groups
GROUPS = ['CTL', 'HFD', 'AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD']
TP_ORDER = ['T0h', 'T6h', 'T24h', 'T72h', 'S72h', 'SD7', 'SD10', 'SD14']

alpha = alpha[alpha['group'].isin(GROUPS)].copy()
pcoa = pcoa[pcoa['group'].isin(GROUPS)].copy()
bc_base = bc_base[bc_base['group'].isin(GROUPS)].copy()

# Color scheme: ND=blue shades, HFD=red shades
GROUP_COLORS = {
    'CTL': '#3498db', 'HFD': '#e74c3c',
    'AMX_CTL': '#2980b9', 'AMX_HFD': '#c0392b',
    'CIP_CTL': '#1abc9c', 'CIP_HFD': '#e67e22',
    'ERY_CTL': '#9b59b6', 'ERY_HFD': '#f39c12',
}

GROUP_MARKERS = {
    'CTL': 'o', 'HFD': 'o',
    'AMX_CTL': 's', 'AMX_HFD': 's',
    'CIP_CTL': '^', 'CIP_HFD': '^',
    'ERY_CTL': 'D', 'ERY_HFD': 'D',
}

# ═══════════════════════════════════════════════════════
# Fig 2: 4-panel figure
# ═══════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# ── Panel A: Shannon diversity over time ──
ax = axes[0, 0]
alpha['timepoint'] = pd.Categorical(alpha['timepoint'], categories=TP_ORDER, ordered=True)

for group in GROUPS:
    sub = alpha[alpha['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['shannon'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'],
                marker=GROUP_MARKERS[group], color=GROUP_COLORS[group],
                label=group, linewidth=1.5, capsize=3, markersize=5)

ax.set_xticks(range(len(TP_ORDER)))
ax.set_xticklabels(TP_ORDER, rotation=45, ha='right')
ax.set_ylabel("Shannon index")
ax.set_title("A) Shannon diversity", fontweight='bold')
ax.axvline(x=3.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.text(1.5, ax.get_ylim()[1]*0.98, 'Treatment', ha='center', fontsize=8, color='gray')
ax.text(5.5, ax.get_ylim()[1]*0.98, 'Recovery', ha='center', fontsize=8, color='gray')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=6, ncol=2, loc='lower left', framealpha=0.9)

# ── Panel B: Observed ASVs (richness) over time ──
ax = axes[0, 1]
for group in GROUPS:
    sub = alpha[alpha['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['observed_ASVs'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'],
                marker=GROUP_MARKERS[group], color=GROUP_COLORS[group],
                label=group, linewidth=1.5, capsize=3, markersize=5)

ax.set_xticks(range(len(TP_ORDER)))
ax.set_xticklabels(TP_ORDER, rotation=45, ha='right')
ax.set_ylabel("Observed ASVs")
ax.set_title("B) Richness (observed ASVs)", fontweight='bold')
ax.axvline(x=3.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=6, ncol=2, loc='lower left', framealpha=0.9)

# ── Panel C: PCoA (Bray-Curtis) at T72h ──
ax = axes[1, 0]

# Get variance explained (the variance_pct column is already a percentage)
var_exp = eigen['variance_pct'].values
pct1, pct2 = var_exp[0], var_exp[1]

for group in GROUPS:
    sub = pcoa[pcoa['group'] == group]
    ax.scatter(sub['PCoA1'], sub['PCoA2'],
               c=GROUP_COLORS[group], marker=GROUP_MARKERS[group],
               label=group, s=30, alpha=0.7, edgecolors='white', linewidth=0.3)

ax.set_xlabel(f"PCoA1 ({pct1:.1f}%)")
ax.set_ylabel(f"PCoA2 ({pct2:.1f}%)")
ax.set_title("C) PCoA (Bray-Curtis, all timepoints)", fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=6, ncol=2, loc='best', framealpha=0.9)

# ── Panel D: Bray-Curtis dissimilarity from baseline ──
ax = axes[1, 1]
bc_base['timepoint'] = pd.Categorical(bc_base['timepoint'], categories=TP_ORDER, ordered=True)

for group in GROUPS:
    sub = bc_base[bc_base['group'] == group]
    grouped = sub.groupby('timepoint', observed=True)['bc_from_baseline'].agg(['mean', 'sem']).reindex(TP_ORDER)
    valid = grouped.dropna()
    if len(valid) == 0:
        continue
    x = [TP_ORDER.index(t) for t in valid.index]
    ax.errorbar(x, valid['mean'], yerr=valid['sem'],
                marker=GROUP_MARKERS[group], color=GROUP_COLORS[group],
                label=group, linewidth=1.5, capsize=3, markersize=5)

ax.set_xticks(range(len(TP_ORDER)))
ax.set_xticklabels(TP_ORDER, rotation=45, ha='right')
ax.set_ylabel("Bray-Curtis dissimilarity from baseline")
ax.set_title("D) Community displacement from baseline", fontweight='bold')
ax.axvline(x=3.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=6, ncol=2, loc='lower right', framealpha=0.9)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig2_diversity.pdf")
plt.savefig(f"{OUTDIR}/fig2_diversity.png")
plt.close()
print("Fig 2 done")
