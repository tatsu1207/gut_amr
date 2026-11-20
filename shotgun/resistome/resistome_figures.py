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
# Fig 4A: Total ARG abundance by group across timepoints
# ============================================================
fig, ax = plt.subplots(figsize=(10, 5))

x_positions = []
x_labels = []
pos = 0
for tp in tp_order:
    for grp in grp_order:
        subset = df[(df['timepoint'] == tp) & (df['group'] == grp)]
        if subset.empty:
            continue
        vals = subset['total_ARG']
        bar = ax.bar(pos, vals.mean(), width=0.7, color=grp_colors[grp],
                     edgecolor='black', linewidth=0.5)
        ax.errorbar(pos, vals.mean(), yerr=vals.std(), fmt='none',
                    color='black', capsize=3, linewidth=1)
        # Individual points
        ax.scatter([pos] * len(vals), vals, color='black', s=12, zorder=5, alpha=0.6)
        x_positions.append(pos)
        x_labels.append(f"{grp}\n{tp}")
        pos += 1
    pos += 0.8  # gap between timepoints

ax.set_xticks(x_positions)
ax.set_xticklabels(x_labels, fontsize=6, rotation=45, ha='right')
ax.set_ylabel("Total ARG abundance\n(copies / 16S rRNA gene copy)", fontsize=10)
ax.set_title("Total ARG Abundance Across Groups and Timepoints", fontsize=12, fontweight='bold')

# Add timepoint labels
tp_starts = [0, 0, 0]
idx = 0
for i, tp in enumerate(tp_order):
    start = idx
    for grp in grp_order:
        if not df[(df['timepoint'] == tp) & (df['group'] == grp)].empty:
            idx += 1
    tp_starts[i] = (start + idx - 1) / 2
    idx += 0  # account for gap differently

ax.set_ylim(bottom=0)  # SD bars of skewed n = 3 groups would otherwise run below zero
plt.tight_layout()
fig.savefig(OUT / "fig4a_total_arg.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4a_total_arg.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4a_total_arg.png'}")
plt.close()

# ============================================================
# Fig 4B: ARG class composition stacked bars
# ============================================================
# Top drug classes by mean abundance
top_classes = ['tetracycline antibiotic', 'macrolide antibiotic',
               'lincosamide antibiotic', 'streptogramin antibiotic',
               'pleuromutilin antibiotic', 'phenicol antibiotic',
               'glycopeptide antibiotic', 'oxazolidinone antibiotic',
               'penicillin beta-lactam', 'fluoroquinolone antibiotic',
               'cephalosporin', 'aminoglycoside antibiotic']

# Verify which exist
top_classes = [c for c in top_classes if c in frac.columns]

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
class_colors = sns.color_palette("Set3", len(top_classes)) + [(.4,.4,.4)]

for ax_idx, tp in enumerate(tp_order):
    ax = axes[ax_idx]
    tp_data = df[df['timepoint'] == tp].copy()
    tp_data = tp_data.sort_values('group')

    groups_present = [g for g in grp_order if g in tp_data['group'].values]
    means = frac.loc[tp_data.index].groupby(tp_data['group'], observed=True)[top_classes].mean()
    means = means.loc[groups_present]
    means['Other'] = tp_data.groupby('group', observed=True)['total_ARG'].mean().loc[groups_present] - means.sum(axis=1)
    means['Other'] = means['Other'].clip(lower=0)

    bottom = np.zeros(len(means))
    for i, cls in enumerate(top_classes + ['Other']):
        vals = means[cls].values
        ax.bar(range(len(means)), vals, bottom=bottom,
               color=class_colors[i], edgecolor='white', linewidth=0.3,
               label=cls.replace(' antibiotic', '') if ax_idx == 0 else None,
               width=0.75)
        bottom += vals

    ax.set_xticks(range(len(means)))
    ax.set_xticklabels(means.index, rotation=45, ha='right', fontsize=8)
    ax.set_title(tp, fontsize=12, fontweight='bold')
    if ax_idx == 0:
        ax.set_ylabel("ARG abundance\n(copies / 16S copy;\nmulti-class genes split equally)", fontsize=10)

fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=8,
           title='Drug class', title_fontsize=9)
fig.suptitle("ARG Composition by Drug Class", fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(OUT / "fig4b_arg_class_composition.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4b_arg_class_composition.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4b_arg_class_composition.png'}")
plt.close()

# ============================================================
# Fig 4C: Fold change from baseline heatmap (by group)
# ============================================================
# Calculate fold change T72h/T0h and SD14/T0h for key classes
key_classes = ['tetracycline antibiotic', 'macrolide antibiotic',
               'penicillin beta-lactam', 'fluoroquinolone antibiotic',
               'lincosamide antibiotic', 'phenicol antibiotic',
               'aminoglycoside antibiotic', 'cephalosporin']
key_classes = [c for c in key_classes if c in simple.columns]

fc_rows = []
for grp in grp_order:
    t0 = df[(df['group'] == grp) & (df['timepoint'] == 'T0h')]
    for tp in ['T72h', 'SD14']:
        tx = df[(df['group'] == grp) & (df['timepoint'] == tp)]
        if t0.empty or tx.empty:
            continue
        for cls in key_classes:
            baseline = t0[cls].mean()
            treated = tx[cls].mean()
            if baseline > 0.0001:
                fc = np.log2(treated / baseline)
            else:
                fc = np.log2((treated + 0.001) / (baseline + 0.001))
            fc_rows.append({'group': grp, 'timepoint': tp,
                           'drug_class': cls.replace(' antibiotic', ''), 'log2FC': fc})

fc_df = pd.DataFrame(fc_rows)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax_idx, tp in enumerate(['T72h', 'SD14']):
    ax = axes[ax_idx]
    pivot = fc_df[fc_df['timepoint'] == tp].pivot(index='drug_class', columns='group', values='log2FC')
    pivot = pivot[[g for g in grp_order if g in pivot.columns]]

    sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.1f',
                ax=ax, cbar_kws={'label': 'log2(fold change vs T0h)'},
                linewidths=0.5, vmin=-4, vmax=4)
    ax.set_title(f"Fold Change at {tp} vs Baseline", fontsize=11, fontweight='bold')
    ax.set_ylabel("Drug Class" if ax_idx == 0 else "")
    ax.set_xlabel("")

plt.tight_layout()
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
ax.set_title("25 ARGs with the Largest Change in Group Mean, T0h to T72h",
             fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(OUT / "fig4d_top_args_heatmap.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4d_top_args_heatmap.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4d_top_args_heatmap.png'}")
plt.close()

# ============================================================
# Fig 4E: HFD vs ND paired comparison (total ARG, T72h and SD14)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

pairs = [('AMX_CTL', 'AMX_HFD', 'Amoxicillin'),
         ('CIP_CTL', 'CIP_HFD', 'Ciprofloxacin'),
         ('ERY_CTL', 'ERY_HFD', 'Erythromycin'),
         ('CTL', 'HFD', 'No antibiotic')]

for ax_idx, tp in enumerate(['T72h', 'SD14']):
    ax = axes[ax_idx]
    x = np.arange(len(pairs))
    width = 0.35

    nd_vals = []
    hfd_vals = []
    nd_err = []
    hfd_err = []

    for nd_grp, hfd_grp, label in pairs:
        nd_data = df[(df['group'] == nd_grp) & (df['timepoint'] == tp)]['total_ARG']
        hfd_data = df[(df['group'] == hfd_grp) & (df['timepoint'] == tp)]['total_ARG']
        nd_vals.append(nd_data.mean() if len(nd_data) > 0 else 0)
        hfd_vals.append(hfd_data.mean() if len(hfd_data) > 0 else 0)
        nd_err.append(nd_data.std() if len(nd_data) > 0 else 0)
        hfd_err.append(hfd_data.std() if len(hfd_data) > 0 else 0)

    bars1 = ax.bar(x - width/2, nd_vals, width, yerr=nd_err, capsize=4,
                   color='#4477AA', edgecolor='black', linewidth=0.5, label='ND')
    bars2 = ax.bar(x + width/2, hfd_vals, width, yerr=hfd_err, capsize=4,
                   color='#CC6677', edgecolor='black', linewidth=0.5, label='HFD')
    # npj guideline: bar charts must show individual data points (n = 3 mice per bar)
    _rng = np.random.default_rng(1)
    for k, (nd_grp, hfd_grp, _) in enumerate(pairs):
        for off, g in ((-width/2, nd_grp), (width/2, hfd_grp)):
            v = df[(df['group'] == g) & (df['timepoint'] == tp)]['total_ARG'].values
            ax.scatter(k + off + _rng.uniform(-0.06, 0.06, len(v)), v, s=14,
                       color='black', alpha=0.7, zorder=5)
    ax.set_ylim(bottom=0)  # SD bar of the skewed AMX_CTL T72h group would run below zero

    ax.set_xticks(x)
    ax.set_xticklabels([p[2] for p in pairs], fontsize=9)
    ax.set_ylabel("Total ARG abundance\n(copies / 16S copy)" if ax_idx == 0 else "")
    ax.set_title(f"{tp}", fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)

fig.suptitle("Diet Effect on Total ARG Abundance", fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(OUT / "fig4e_hfd_vs_nd_arg.png", dpi=300, bbox_inches='tight')
fig.savefig(OUT / "fig4e_hfd_vs_nd_arg.pdf", bbox_inches='tight')
print(f"Saved: {OUT / 'fig4e_hfd_vs_nd_arg.png'}")
plt.close()

print("\nAll resistome figures done.")
