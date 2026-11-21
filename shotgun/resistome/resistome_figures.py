#!/usr/bin/env python3
"""Resistome figures (Fig. 4A-E) for the KJM manuscript.

Usage: python3 resistome_figures_rrna_filtered.py [homolog_only|no_rrna|all]
Default homolog_only (the main analysis since rev. 10.6). Panels are written to
hiseq/figures/, where compose_multipanel_figures.py reads them.
"""
import sys

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

LEVEL = sys.argv[1] if len(sys.argv) > 1 else "homolog_only"
OUT = Path("/home/unnot/vault/kjm/hiseq/figures") if LEVEL == "homolog_only" else Path(f"/home/unnot/vault/kjm/hiseq/figures/{LEVEL}")
OUT.mkdir(parents=True, exist_ok=True)

# ── Load data ──
simple = pd.read_csv(f"/home/unnot/vault/kjm/hiseq/arg_results/model_filtered/arg_abundance_simple_class_{LEVEL}.tsv",
                      sep="\t", index_col=0)
pergene = pd.read_csv(f"/home/unnot/vault/kjm/hiseq/arg_results/model_filtered/arg_abundance_per_gene_{LEVEL}.tsv",
                       sep="\t")

# ── Sample metadata from sample names ──
group_map = {
    'AC': ('AMX_CTL', 'Amoxicillin', 'ND'),
    'AH': ('AMX_HFD', 'Amoxicillin', 'HFD'),
    'CC': ('CIP_CTL', 'Ciprofloxacin', 'ND'),
    'CH': ('CIP_HFD', 'Ciprofloxacin', 'HFD'),
    'EC': ('ERY_CTL', 'Erythromycin', 'ND'),
    'EH': ('ERY_HFD', 'Erythromycin', 'HFD'),
    'CTL': ('CTL', 'None', 'ND'),
    'HFD': ('HFD', 'None', 'HFD'),
    'CSL': ('CTL', 'None', 'ND'),  # CSL = mislabelled CTL at recovery timepoints (same mice)
}

def parse_sample(s):
    s2 = s.replace('ARG_', '')
    for prefix in sorted(group_map.keys(), key=len, reverse=True):
        rest = s2[len(prefix):]
        if s2.startswith(prefix) and (rest == '' or rest[0] in '0123456789_'):
            tp = s.rsplit('_', 1)[1]
            grp, abx, diet = group_map[prefix]
            return grp, abx, diet, tp
    return None, None, None, None

meta = pd.DataFrame([parse_sample(s) for s in simple.index],
                    columns=['group', 'antibiotic', 'diet', 'timepoint'],
                    index=simple.index)

# Merge
df = simple.join(meta)
# Total ARG = sum over GENES (per-gene table). Do NOT sum the simple_class columns:
# a gene whose CARD drug-class string lists k classes appears in all k columns, so a
# row sum over classes counts it k times (fixed 2026-11-20; see review point R9).
_pg_num = pergene.drop(columns=['gene', 'gene_name', 'aro', 'drug_class',
                                'resistance_mechanism', 'gene_family'], errors='ignore')
_pg_num = _pg_num.select_dtypes('number')
df['total_ARG'] = _pg_num.sum(axis=0).reindex(df.index)
assert df['total_ARG'].notna().all()

# Fractional class allocation for the composition plot (Fig 4B): each gene's abundance
# is split equally among the k drug classes in its CARD string, so stacked class
# segments sum exactly to the per-gene total.
_k = pergene['drug_class'].str.count(';') + 1
frac = {}
for _i, _r in pergene.iterrows():
    for _c in _r['drug_class'].split(';'):
        frac[_c] = frac.get(_c, 0) + _pg_num.loc[_i] / _k[_i]
frac = pd.DataFrame(frac).reindex(df.index)
assert np.allclose(frac.sum(axis=1), df['total_ARG'])

# Timepoint order
tp_order = ['T0h', 'T72h', 'SD14']
df['timepoint'] = pd.Categorical(df['timepoint'], categories=tp_order, ordered=True)

# Group order
grp_order = ['CTL', 'HFD', 'AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD']
df['group'] = pd.Categorical(df['group'], categories=grp_order, ordered=True)

# Color scheme
grp_colors = {
    'CTL': '#4477AA', 'HFD': '#CC6677',
    'AMX_CTL': '#228833', 'AMX_HFD': '#EE6677',
    'CIP_CTL': '#66CCEE', 'CIP_HFD': '#AA3377',
    'ERY_CTL': '#CCBB44', 'ERY_HFD': '#EE8866',
}

# ============================================================
# Fig 4A: Total ARG abundance, ND and HFD side by side for each antibiotic and timepoint
# (replaces the former Fig 4E diet panel, which replotted the same values)
# ============================================================
DIET_COL = {'ND': '#4477AA', 'HFD': '#CC6677'}
PAIRS = [('CTL', 'HFD', 'None'), ('AMX_CTL', 'AMX_HFD', 'AMX'),
         ('CIP_CTL', 'CIP_HFD', 'CIP'), ('ERY_CTL', 'ERY_HFD', 'ERY')]
fig, ax = plt.subplots(figsize=(11, 4.8))
rng = np.random.default_rng(1)
pos, xt, xl, tp_mid = 0.0, [], [], []
for tp in tp_order:
    start = pos
    for nd_g, hfd_g, abx in PAIRS:
        for off, g, diet in ((-0.2, nd_g, 'ND'), (0.2, hfd_g, 'HFD')):
            vals = df[(df['timepoint'] == tp) & (df['group'] == g)]['total_ARG']
            ax.bar(pos + off, vals.mean(), width=0.38, color=DIET_COL[diet], edgecolor='black', linewidth=0.5)
            ax.errorbar(pos + off, vals.mean(), yerr=vals.std(), fmt='none', color='black', capsize=2.5, linewidth=0.9)
            ax.scatter(pos + off + rng.uniform(-0.06, 0.06, len(vals)), vals, color='black', s=10, zorder=5, alpha=0.7)
        xt.append(pos); xl.append(abx); pos += 1
    tp_mid.append((start + pos - 1) / 2); pos += 0.7
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=9)
for m, tp in zip(tp_mid, tp_order):
    ax.text(m, -0.12, tp, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=11, fontweight='bold')
ax.set_ylabel("Total ARG abundance\n(copies / 16S rRNA gene copy)", fontsize=10)
ax.set_ylim(bottom=0)  # SD bars of skewed n = 3 groups would otherwise run below zero
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=DIET_COL['ND'], edgecolor='black', label='Normal diet (ND)'),
                   Patch(facecolor=DIET_COL['HFD'], edgecolor='black', label='High-fat diet (HFD)')],
          fontsize=9, loc='upper left', frameon=False)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_title("Total ARG abundance by diet, antibiotic and timepoint", fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(OUT / "fig4a_total_arg.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4a_total_arg.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4a_total_arg.png'}")
plt.close()

# ============================================================
# Fig 4B: ARG composition, one bar per mouse (2026-11-21). Group means had hidden that one AMX_CTL
# animal at T3 made up 85% of its group's bar, and splitting multidrug efflux genes across their
# classes made them look like class-specific resistance; see class_categories.py.
# ============================================================
import sys as _s; _s.path.insert(0, "/home/unnot/vault/kjm/hiseq")
from class_categories import categorise, ORDER as CAT_ORDER, MULTI
cat = categorise(pergene, list(df.index))
assert np.allclose(cat.sum(axis=1), df['total_ARG'])
CAT_COL = {"tetracycline": "#88CCEE", "lincosamide": "#AA4499", "macrolide": "#CC6677",
           "streptogramin": "#882255", "pleuromutilin": "#DDCC77", "glycopeptide": "#117733",
           "fluoroquinolone": "#332288", "beta-lactam": "#E69F00", "aminoglycoside": "#44AA99",
           "phenicol": "#999933", "peptide": "#DDDDDD", "other": "#AAAAAA", MULTI: "#444444"}
cats = [c for c in CAT_ORDER if cat[c].max() > 0]

fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
for ax_idx, tp in enumerate(tp_order):
    ax = axes[ax_idx]
    pos, xt, xl, x = 0.0, [], [], []
    order = []
    for g in grp_order:
        ids = sorted(df.index[(df['timepoint'] == tp) & (df['group'] == g)])
        xs = [pos + i for i in range(len(ids))]
        order += ids; x += xs
        if ids:
            xt.append(np.mean(xs)); xl.append(g)
        pos += len(ids) + 0.8
    bottom = np.zeros(len(order))
    for c in cats:
        v = cat.loc[order, c].values
        ax.bar(x, v, bottom=bottom, width=0.85, color=CAT_COL[c], edgecolor='white', linewidth=0.3,
               label=c if ax_idx == 0 else None)
        bottom += v
    YMAX = 0.6   # one T3 bar (AMX_CTL mouse 7, 1.58) is truncated and labelled with its value
    for xi, tot in zip(x, bottom):
        if tot > YMAX:
            ax.text(xi, YMAX * 1.005, f"\u2191 {tot:.2f}", ha='center', va='bottom', fontsize=9,
                    color='black', clip_on=False)
    ax.set_ylim(0, YMAX)
    ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=45, ha='right', fontsize=8)
    ax.set_title(tp, fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    if ax_idx == 0:
        ax.set_ylabel("ARG abundance\n(copies / 16S rRNA gene copy)", fontsize=10)
fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12, title='Drug class', title_fontsize=13)
fig.suptitle("ARG composition by drug class (one bar per mouse)", fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(OUT / "fig4b_arg_class_composition.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4b_arg_class_composition.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4b_arg_class_composition.png'}")
plt.close()

# ============================================================
# Fig 4C: Fold change from baseline (2026-11-21). Same categories as Fig. 4B (class_categories.py;
# multidrug genes no longer counted in full under every class they list); each mouse compared with
# its own T0 (the same three mice were sampled at T0, T3 and R14), then the median of the three, so a
# single animal cannot set a group's value; one pseudocount (0.001 copies/16S) for every cell, so
# near-absent classes do not produce large spurious ratios.
# ============================================================
PSEUDO = 1e-3
FC_CATS = [c for c in ["tetracycline", "lincosamide", "macrolide", "streptogramin", "glycopeptide",
                       "fluoroquinolone", "beta-lactam", "aminoglycoside", MULTI]]   # phenicol: no class-specific genes once multidrug is separate
mouse = pd.Series({sid: sid.rsplit('_', 1)[0].replace('CSL', 'CTL') for sid in df.index})   # CSL = untreated CTL mice, relabelled at recovery
fc_rows = []
for grp in grp_order:
    for tp in ['T72h', 'SD14']:
        ids = [i for i in df.index[(df['group'] == grp) & (df['timepoint'] == tp)]]
        for c in FC_CATS:
            v = []
            for i in ids:
                b0 = [j for j in df.index[(df['timepoint'] == 'T0h')] if mouse[j] == mouse[i]]
                if b0:
                    v.append(np.log2((cat.loc[i, c] + PSEUDO) / (cat.loc[b0[0], c] + PSEUDO)))
            if v:
                fc_rows.append({'group': grp, 'timepoint': tp, 'drug_class': c, 'log2FC': float(np.median(v))})
fc_df = pd.DataFrame(fc_rows)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax_idx, tp in enumerate(['T72h', 'SD14']):
    ax = axes[ax_idx]
    pivot = fc_df[fc_df['timepoint'] == tp].pivot(index='drug_class', columns='group', values='log2FC')
    pivot = pivot[[g for g in grp_order if g in pivot.columns]].reindex(FC_CATS)

    sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.1f',
                ax=ax, cbar_kws={'label': 'log2 fold change vs own T0\n(median of 3 mice)'},
                linewidths=0.5, vmin=-4, vmax=4)
    # Outline each antibiotic's own target class(es): AMX -> beta-lactam, CIP -> fluoroquinolone,
    # ERY -> macrolide, lincosamide, streptogramin (MLS).
    from matplotlib.patches import Rectangle
    OWN = {"AMX": ["beta-lactam"], "CIP": ["fluoroquinolone"], "ERY": ["lincosamide", "macrolide", "streptogramin"]}
    for ci, g in enumerate(pivot.columns):
        for cls in OWN.get(g.split("_")[0], []):
            ri = list(pivot.index).index(cls)
            ax.add_patch(Rectangle((ci + 0.04, ri + 0.04), 0.92, 0.92, fill=False, edgecolor="black", lw=2.2, zorder=5))
    ax.set_title(f"Fold change at {tp} vs baseline", fontsize=11, fontweight='bold')
    ax.set_ylabel("Drug class" if ax_idx == 0 else "")
    ax.set_xlabel("")

plt.tight_layout()
fig.text(0.5, -0.04, "Black outline: the drug class targeted by that group's antibiotic (AMX, beta-lactam; CIP, fluoroquinolone; ERY, MLS)",
         ha="center", fontsize=10)
fig.savefig(OUT / "fig4c_arg_foldchange_heatmap.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4c_arg_foldchange_heatmap.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4c_arg_foldchange_heatmap.png'}")
plt.close()

# ============================================================
# Fig 4D: Top differentially abundant ARGs (per-gene level)
# ============================================================
gene_meta = pergene[['gene_name', 'drug_class', 'resistance_mechanism', 'gene_family']].copy()
gene_abund = pergene.drop(columns=['gene', 'gene_name', 'aro', 'drug_class',
                                    'resistance_mechanism', 'gene_family'])
gene_abund.index = pergene['gene_name']

# Calculate mean per group per timepoint
gene_grp = {}
for sample in gene_abund.columns:
    grp, abx, diet, tp = parse_sample(sample)
    if grp:
        key = f"{grp}_{tp}"
        gene_grp.setdefault(key, []).append(sample)

gene_means = pd.DataFrame({k: gene_abund[v].mean(axis=1) for k, v in gene_grp.items()})

# Find genes with biggest changes: T72h vs T0h across antibiotic groups
change_scores = {}
for grp in ['AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD']:
    t0_col = f"{grp}_T0h"
    t72_col = f"{grp}_T72h"
    if t0_col in gene_means.columns and t72_col in gene_means.columns:
        diff = gene_means[t72_col] - gene_means[t0_col]
        for gene, val in diff.items():
            change_scores.setdefault(gene, []).append(abs(val))

# Top 25 genes by max absolute change
top_genes_score = {g: max(vals) for g, vals in change_scores.items()}
top_genes = sorted(top_genes_score, key=top_genes_score.get, reverse=True)[:25]

# Build heatmap: genes × (group_timepoint) for antibiotic-treated groups
cols_for_heatmap = []
for grp in ['AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD', 'CTL', 'HFD']:
    for tp in tp_order:
        col = f"{grp}_{tp}"
        if col in gene_means.columns:
            cols_for_heatmap.append(col)

heatmap_data = gene_means.loc[top_genes, cols_for_heatmap]

fig, ax = plt.subplots(figsize=(16, 8))
sns.heatmap(heatmap_data, cmap='YlOrRd', ax=ax, linewidths=0.3,
            cbar_kws={'label': 'ARG abundance (copies / 16S copy)'},
            xticklabels=True, yticklabels=True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
ax.set_title("25 ARGs with the largest change in group mean, T0h to T72h",
             fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(OUT / "fig4d_top_args_heatmap.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4d_top_args_heatmap.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4d_top_args_heatmap.png'}")
plt.close()

print("\nAll resistome figures done.")
