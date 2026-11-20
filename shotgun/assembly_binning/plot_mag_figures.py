#!/usr/bin/env python3
"""Generate figures for MAG/assembly results."""

import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.expanduser("~/vault/kjm/miseq"))
import tplabels; tplabels.install()  # show T0/T1/T3/R3...R14 instead of the data codes
import matplotlib as mpl
import numpy as np

mpl.rcParams.update({
    'font.size': 10, 'font.family': 'sans-serif',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 300, 'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

OUTDIR = "/home/unnot/vault/kjm/hiseq/figures"
MAG_DIR = "/home/unnot/vault/kjm/hiseq/mag_analysis"

# ── Load data ──
full = pd.read_csv(f"{MAG_DIR}/table_arg_mag_taxonomy.tsv", sep="\t")
genus_summary = pd.read_csv(f"{MAG_DIR}/table_genus_arg_summary.tsv", sep="\t")
mls_by_genus = pd.read_csv(f"{MAG_DIR}/table_mls_arg_by_genus.tsv", sep="\t")
lacto = pd.read_csv(f"{MAG_DIR}/table_lactococcus_args.tsv", sep="\t")
per_sample = pd.read_csv(f"{MAG_DIR}/table_per_sample_summary.tsv", sep="\t")


# ═══════════════════════════════════════════════════════
# Fig 5A: Assembly-based ARG drug class composition (pie)
# ═══════════════════════════════════════════════════════
def simplify_drug_class(dc):
    dc = str(dc).lower()
    if "glycopeptide" in dc:
        return "Glycopeptide"
    elif "nitroimidazole" in dc:
        return "Nitroimidazole"
    elif "tetracycline" in dc:
        return "Tetracycline"
    elif "rifamycin" in dc:
        return "Rifamycin"
    elif "macrolide" in dc or "lincosamide" in dc or "streptogramin" in dc:
        return "MLS"
    elif "fluoroquinolone" in dc:
        return "Fluoroquinolone"
    elif "beta-lactam" in dc or "cephalosporin" in dc or "carbapenem" in dc or "penicillin" in dc or "monobactam" in dc:
        return "Beta-lactam"
    elif "aminoglycoside" in dc:
        return "Aminoglycoside"
    elif "phenicol" in dc:
        return "Phenicol"
    elif "disinfect" in dc:
        return "Disinfectants"
    elif "peptide antibiotic" in dc:
        return "Peptide"
    elif "diaminopyrimidine" in dc:
        return "Diaminopyrimidine"
    elif "oxazolidinone" in dc:
        return "Oxazolidinone"
    else:
        return "Other"

# Same rule as the Results text: a gene counts toward a class only if CARD assigns it that
# single class; genes with several classes in their CARD string form their own slice.
def single_class(dc):
    dc = str(dc)
    if ";" in dc:
        return "Multi-class"
    return dc.replace(" antibiotic", "").replace("disinfecting agents and antiseptics", "disinfectants").capitalize()
full["simple_class"] = full["drug_class"].apply(single_class)
class_counts = full["simple_class"].value_counts()

# Merge small categories
threshold = 200
major = class_counts[class_counts >= threshold]
other_sum = class_counts[class_counts < threshold].sum()
if other_sum > 0:
    major["Other"] = other_sum

colors_pie = plt.cm.Set3(np.linspace(0, 1, len(major)))

fig, ax = plt.subplots(figsize=(6, 6))
# Small slices are labelled in a legend, not on the pie, so their labels cannot overlap.
wedges, _ = ax.pie(
    major.values, colors=colors_pie, startangle=90,
    wedgeprops=dict(linewidth=0.5, edgecolor='white')
)
ax.legend(wedges, [f"{k} ({v:,}; {v / major.sum() * 100:.1f}%)" for k, v in major.items()],
          loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=2, fontsize=9, frameon=False)
ax.set_title("Assembly-based ARG drug class distribution\n(n = 17,371)", fontsize=11, fontweight='bold')
plt.savefig(f"{OUTDIR}/fig5a_assembly_arg_classes.pdf")
plt.savefig(f"{OUTDIR}/fig5a_assembly_arg_classes.png")
plt.close()
print("Fig 5A done")


# ═══════════════════════════════════════════════════════
# Fig 5B: Top 20 ARG-carrying genera (horizontal bar)
# ═══════════════════════════════════════════════════════
top20 = genus_summary.head(20).copy()
top20 = top20.iloc[::-1]  # reverse for horizontal bar

fig, ax = plt.subplots(figsize=(7, 6))
colors_bar = ['#e74c3c' if g == 'Lactococcus' else '#3498db' for g in top20['genus']]
bars = ax.barh(top20['genus'], top20['total_ARGs'], color=colors_bar, edgecolor='white', linewidth=0.5)

# Add MAG counts
for i, (_, row) in enumerate(top20.iterrows()):
    ax.text(row['total_ARGs'] + 5, i, f"({row['n_MAGs']} MAGs)", va='center', fontsize=7, color='gray')

ax.set_xlabel("Total ARGs", fontsize=10)
ax.set_title("Top 20 ARG-carrying genera\n(assembly-based)", fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add Lactococcus annotation if in top 20
if 'Lactococcus' in top20['genus'].values:
    ax.annotate('Lactococcus\n(rank 39)', xy=(0, 0), fontsize=8, color='#e74c3c')

plt.savefig(f"{OUTDIR}/fig5b_top_arg_genera.pdf")
plt.savefig(f"{OUTDIR}/fig5b_top_arg_genera.png")
plt.close()
print("Fig 5B done")


# ═══════════════════════════════════════════════════════
# Fig 5C: Lactococcus ARGs vs MLS ARGs by genus
# ═══════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios': [1, 1.2]})

# Left: Lactococcus ARG types
lacto_args = lacto.groupby('best_hit_ARG')['best_hit_ARG'].count().sort_values(ascending=False)
lacto_classes = lacto.groupby('drug_class')['drug_class'].count().sort_values(ascending=True)

colors_lacto = {'glycopeptide antibiotic': '#2ecc71', 'rifamycin antibiotic': '#9b59b6', 'aminoglycoside antibiotic': '#e67e22'}
bar_colors = [colors_lacto.get(c, '#95a5a6') for c in lacto_classes.index]

# lmrD (lincosamide efflux, CARD reference = L. lactis) is missed by RGI Strict under the
# DIAMOND aligner (bitscore 1248 < cutoff 1320); found by a direct DIAMOND search
# (mag_analysis/lmrD/, added 2026-11-20). Shown hatched, separately from RGI calls.
lmrd = pd.read_csv(f"{MAG_DIR}/lmrD/table_lmrD_lactococcus_MAGs_prodigal.tsv", sep="\t")
n_lmrd = int((lmrd['status'] != 'absent').sum())
n_full = int((lmrd['status'] == 'full-length').sum()); n_part = int((lmrd['status'] == 'partial').sum())
lacto_classes = pd.concat([pd.Series({'lincosamide (lmrD)*': n_lmrd}), lacto_classes.sort_values()]).sort_values()
bar_colors = ['#c0392b' if c.startswith('lincosamide') else colors_lacto.get(c, '#95a5a6') for c in lacto_classes.index]
bars = ax1.barh(lacto_classes.index, lacto_classes.values, color=bar_colors, edgecolor='white')
for b, c in zip(bars, lacto_classes.index):
    if c.startswith('lincosamide'):
        b.set_hatch('///'); b.set_edgecolor('black')
ax1.set_xlabel("ARG ORF count", fontsize=10)
ax1.set_title(f"Lactococcus MAG ARGs\n({len(lmrd)} MAGs)", fontsize=11, fontweight='bold')
ax1.text(0, -0.22, f"*lmrD: direct DIAMOND search; RGI cannot call it with the DIAMOND aligner,\n so absent from RGI calls ({n_full} full-length, {n_part} partial)",
         transform=ax1.transAxes, fontsize=7, va='top')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for i, v in enumerate(lacto_classes.values):
    ax1.text(v + 0.3, i, str(v), va='center', fontsize=9)

# Right: MLS ARGs by genus
mls_plot = mls_by_genus.head(10).copy()
mls_plot = mls_plot.iloc[::-1]

# Color by type: Gram-negative efflux vs Firmicutes mef(H) vs other
def mls_color(row):
    names = str(row['ARG_names'])
    if 'CRP' in names or 'emrE' in names or 'mdtE' in names:
        return '#e74c3c'  # Gram-negative efflux
    elif 'mef' in names or 'Mef' in names:
        return '#f39c12'  # mef(H) efflux
    elif 'Erm' in names or 'efr' in names:
        return '#3498db'  # other MLS
    return '#95a5a6'

bar_colors2 = [mls_color(row) for _, row in mls_plot.iterrows()]
ax2.barh(mls_plot['genus'], mls_plot['n_MLS_ARGs'], color=bar_colors2, edgecolor='white')
ax2.set_xlabel("MLS ARG count", fontsize=10)
ax2.set_title("MLS ARGs by genus\n(RGI Strict calls; lmrD not included)", fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#e74c3c', label='Gram-negative efflux'),
    Patch(facecolor='#f39c12', label='mef(H) efflux'),
    Patch(facecolor='#3498db', label='Other MLS'),
]
ax2.legend(handles=legend_elements, loc='lower right', fontsize=8, framealpha=0.9)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig5c_lactococcus_vs_mls.pdf")
plt.savefig(f"{OUTDIR}/fig5c_lactococcus_vs_mls.png")
plt.close()
print("Fig 5C done")


# ═══════════════════════════════════════════════════════
# Fig 5D: Genomic context (chromosome/plasmid/phage)
# ═══════════════════════════════════════════════════════
context_counts = full["genomic_location"].value_counts()
context_labels = []
context_vals = []
context_colors = []
color_map = {'chromosome': '#3498db', 'plasmid': '#e74c3c', 'phage': '#2ecc71'}

for loc in ['chromosome', 'plasmid', 'phage']:
    if loc in context_counts.index:
        context_labels.append(f"{loc.capitalize()}\n({context_counts[loc]:,}; {context_counts[loc]/len(full)*100:.1f}%)")
        context_vals.append(context_counts[loc])
        context_colors.append(color_map[loc])

fig, ax = plt.subplots(figsize=(5, 5))
wedges, _ = ax.pie(
    context_vals, colors=context_colors,
    startangle=90, wedgeprops=dict(linewidth=0.5, edgecolor='white')
)
ax.legend(wedges, [l.replace("\n", " ") for l in context_labels],
          loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=1, fontsize=10, frameon=False)
ax.set_title("Genomic context of assembly-based ARGs\n(n = 17,371)", fontsize=11, fontweight='bold')
plt.savefig(f"{OUTDIR}/fig5d_genomic_context.pdf")
plt.savefig(f"{OUTDIR}/fig5d_genomic_context.png")
plt.close()
print("Fig 5D done")


# ═══════════════════════════════════════════════════════
# Fig 5E: ARG-carrying ORFs per sample, ND and HFD side by side for each antibiotic and timepoint
# (same layout as Fig 4A; the former version pooled all four treatments into one line per diet)
# ═══════════════════════════════════════════════════════
from matplotlib.patches import Patch
tp_order = ['T0h', 'T72h', 'SD14']
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
            vals = per_sample[(per_sample['timepoint'] == tp) & (per_sample['group'] == g)]['total_ARGs']
            assert len(vals) == 3, (tp, g, len(vals))
            ax.bar(pos + off, vals.mean(), width=0.38, color=DIET_COL[diet], edgecolor='black', linewidth=0.5)
            ax.errorbar(pos + off, vals.mean(), yerr=vals.std(), fmt='none', color='black', capsize=2.5, linewidth=0.9)
            ax.scatter(pos + off + rng.uniform(-0.06, 0.06, len(vals)), vals, color='black', s=10, zorder=5, alpha=0.7)
        xt.append(pos); xl.append(abx); pos += 1
    tp_mid.append((start + pos - 1) / 2); pos += 0.7
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=9)
for m, tp in zip(tp_mid, tp_order):
    ax.text(m, -0.12, tp, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=11, fontweight='bold')
ax.set_ylabel("ARG-carrying ORFs per sample\n(assembly-based)", fontsize=10)
ax.set_ylim(0, 560)
ax.legend(handles=[Patch(facecolor=DIET_COL['ND'], edgecolor='black', label='Normal diet (ND)'),
                   Patch(facecolor=DIET_COL['HFD'], edgecolor='black', label='High-fat diet (HFD)')],
          fontsize=9, loc='upper left', frameon=False, ncol=2)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_title("Assembly-based ARG-carrying ORFs by diet, antibiotic and timepoint", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig5e_assembly_args_by_diet.pdf")
plt.savefig(f"{OUTDIR}/fig5e_assembly_args_by_diet.png")
plt.close()
print("Fig 5E done")


print("\nAll figures saved to", OUTDIR)
