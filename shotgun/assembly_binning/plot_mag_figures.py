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
# Fig 5 style (2026-11-21): neutral greys for panels A-D (colour is reserved for the antibiotics, as in
# Figs 1-4 and 6); E uses the shared scheme (colour = antibiotic, ND light fill, HFD solid fill).
# ═══════════════════════════════════════════════════════
_sys.path.insert(0, _os.path.expanduser("~/vault/kjm/hiseq"))
from class_categories import is_multidrug
from matplotlib.patches import Patch
DARK, MID, LIGHT, PALE = "#3a3a3a", "#7a7a7a", "#b5b5b5", "#dedede"
def clean(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

full["multidrug"] = full["drug_class"].map(is_multidrug)
LIGASE = ("vanG", "vanD", "vanC", "vanA", "vanB", "vanE", "vanL", "vanM", "vanN")
def is_ligase(name):
    return name in LIGASE or any(name.startswith(l + " ") for l in ("vanC",))

# ── Fig 5A: drug-class composition, with the van clusters split into ligases and accessory genes ──
def a_cat(r):
    dc = str(r["drug_class"])
    if dc == "glycopeptide antibiotic":
        return "Glycopeptide: van ligase" if is_ligase(r["best_hit_ARG"]) else "Glycopeptide: van accessory/regulatory genes"
    if r["multidrug"]:
        return "Multidrug (>=3 class families)"
    if ";" in dc:
        return "Two classes"
    if dc == "rifamycin antibiotic":
        return "Rifamycin (mostly rpoC variants)"
    return dc.replace(" antibiotic", "").replace("disinfecting agents and antiseptics", "disinfectants").capitalize()
full["a_cat"] = full.apply(a_cat, axis=1)
cc = full["a_cat"].value_counts()
small = cc[cc < 200].sum()
cc = cc[cc >= 200]
order = ["Glycopeptide: van accessory/regulatory genes", "Glycopeptide: van ligase"] + \
        [k for k in cc.index if not k.startswith("Glycopeptide")]
cc = cc.reindex(order); cc["Other"] = small
fig, ax = plt.subplots(figsize=(7, 6))
y = np.arange(len(cc))[::-1]
cols = [PALE if k.startswith("Glycopeptide: van acc") else (MID if k.startswith("Rifamycin") else DARK) for k in cc.index]
bars = ax.barh(y, cc.values, color=cols, edgecolor=DARK, linewidth=0.6, height=0.7)
bars[0].set_hatch("///")
for yi, v in zip(y, cc.values):
    ax.text(v + 150, yi, f"{v:,} ({v / len(full) * 100:.1f}%)", va="center", fontsize=9)
ax.set_yticks(y); ax.set_yticklabels(cc.index, fontsize=9.5)
ax.set_xlabel("ARG-carrying ORFs (n = 17,371)")
ax.set_xlim(0, cc.max() * 1.32)
ax.set_title("Drug class of assembly-based ARGs", fontsize=11, fontweight="bold")
ax.text(0.98, 0.02, "Hatched: CARD glycopeptide-class genes that do not\nconfer resistance on their own (vanW, vanT, vanY, ...)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")
clean(ax)
plt.savefig(f"{OUTDIR}/fig5a_assembly_arg_classes.pdf"); plt.savefig(f"{OUTDIR}/fig5a_assembly_arg_classes.png")
plt.close(); print("Fig 5A done", dict(cc))

# ── Fig 5B: top 20 genera by total ARGs, plotted as ARGs per MAG (bins of the same organism recur) ──
top20 = genus_summary.head(20).copy()
top20["per_mag"] = top20["total_ARGs"] / top20["n_MAGs"]
top20 = top20.sort_values("per_mag")
fig, ax = plt.subplots(figsize=(7, 6))
ax.barh(top20["genus"], top20["per_mag"], color=MID, edgecolor="white", height=0.72)
for i, (_, r) in enumerate(top20.iterrows()):
    ax.text(r["per_mag"] + 0.4, i, f"{r['total_ARGs']} ARGs / {r['n_MAGs']} MAGs", va="center", fontsize=7.5, color="#555555")
for t in ax.get_yticklabels():
    t.set_fontstyle("italic")
ax.set_xlabel("ARG-carrying ORFs per MAG")
ax.set_xlim(0, top20["per_mag"].max() * 1.45)
ax.set_title("Top 20 ARG-carrying genera (by total)", fontsize=11, fontweight="bold")
clean(ax)
plt.savefig(f"{OUTDIR}/fig5b_top_arg_genera.pdf"); plt.savefig(f"{OUTDIR}/fig5b_top_arg_genera.png")
plt.close(); print("Fig 5B done")

# ── Fig 5C: genomic context (geNomad) ──
ctx = full["genomic_location"].value_counts().reindex(["chromosome", "plasmid", "phage"])
fig, ax = plt.subplots(figsize=(5, 6))
yy = np.arange(3)[::-1]
ax.barh(yy, ctx.values, color=[DARK, MID, LIGHT], height=0.6)
for yi, v in zip(yy, ctx.values):
    ax.text(v + 250, yi, f"{v:,}\n({v / len(full) * 100:.1f}%)", va="center", fontsize=9)
ax.set_yticks(yy); ax.set_yticklabels(["Chromosome", "Plasmid", "Phage"], fontsize=10)
ax.set_xlim(0, ctx.max() * 1.4)
ax.set_xlabel("ARG-carrying ORFs")
ax.set_title("Genomic context", fontsize=11, fontweight="bold")
clean(ax)
plt.savefig(f"{OUTDIR}/fig5c_genomic_context.pdf"); plt.savefig(f"{OUTDIR}/fig5c_genomic_context.png")
plt.close(); print("Fig 5C done")

# ── Fig 5D: Lactococcus MAG ARGs (left) and MLS-class ARGs by genus (right) ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios': [1, 1.2]})
# lmrD (lincosamide efflux, CARD reference = L. lactis) is missed by RGI Strict under the
# DIAMOND aligner (bitscore 1248 < cutoff 1320); found by a direct DIAMOND search (mag_analysis/lmrD/).
lmrd = pd.read_csv(f"{MAG_DIR}/lmrD/table_lmrD_lactococcus_MAGs_prodigal.tsv", sep="\t")
n_lmrd = int((lmrd['status'] != 'absent').sum())
n_full = int((lmrd['status'] == 'full-length').sum()); n_part = int((lmrd['status'] == 'partial').sum())
GENE = {'glycopeptide antibiotic': 'glycopeptide (vanY)', 'rifamycin antibiotic': 'rifamycin (rpoC variant)',
        'aminoglycoside antibiotic': 'aminoglycoside (aadS)'}
lc = lacto.groupby('drug_class').size().rename(index=GENE)
lc = pd.concat([pd.Series({'lincosamide (lmrD)*': n_lmrd}), lc]).sort_values()
bars = ax1.barh(lc.index, lc.values, color=[DARK if k.startswith('lincos') else (LIGHT if k.startswith('rifa') else MID) for k in lc.index],
                edgecolor='white')
for bb, k in zip(bars, lc.index):
    if k.startswith('lincos'):
        bb.set_hatch('///'); bb.set_edgecolor('white')
for i, v in enumerate(lc.values):
    ax1.text(v + 0.3, i, str(v), va='center', fontsize=9)
ax1.set_xlabel("ARG ORF count", fontsize=10)
ax1.set_title(f"Lactococcus MAG ARGs ({len(lmrd)} MAGs)", fontsize=11, fontweight='bold')
ax1.text(0, -0.2, f"*lmrD: direct DIAMOND search; RGI cannot call it with the DIAMOND aligner\n({n_full} full-length, {n_part} partial)",
         transform=ax1.transAxes, fontsize=7.5, va='top')
clean(ax1)

# MLS-class = CARD drug class containing a macrolide, lincosamide, streptogramin or pleuromutilin term;
# binned ORFs only. Multidrug efflux genes (>=3 class families; same rule as Fig. 4B) shown separately.
mls = full[full['bin'].notna() & full['genus'].notna() &
           full['drug_class'].str.contains('macrolide|lincosamide|streptogramin|pleuromutilin', na=False)]
mg = mls.groupby('genus').agg(n=('orf_id', 'size'), md=('multidrug', 'sum'))
mg['spec'] = mg['n'] - mg['md']
mg = mg[(mg['spec'] > 0) | (mg['n'] >= 3)].sort_values(['n', 'spec'], ascending=False).iloc[::-1]   # every genus with an MLS-specific gene
ax2.barh(mg.index, mg['spec'], color=DARK, edgecolor='white', label='MLS-specific (mef, Mef(En2), ErmF, emrE)')
ax2.barh(mg.index, mg['md'], left=mg['spec'], color=LIGHT, edgecolor='white',
         label='Multidrug efflux / regulators (>=3 classes)')
for t in ax2.get_yticklabels():
    t.set_fontstyle('italic')
ax2.set_xlabel("MLS-class ARG ORFs in MAGs", fontsize=10)
ax2.set_title("MLS-class ARGs by genus (RGI; lmrD not included)", fontsize=11, fontweight='bold')
ax2.legend(loc='lower right', fontsize=8, frameon=False)
clean(ax2)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig5d_lactococcus_vs_mls.pdf"); plt.savefig(f"{OUTDIR}/fig5d_lactococcus_vs_mls.png")
plt.close(); print("Fig 5D done"); print(mg)


# ═══════════════════════════════════════════════════════
# Fig 5E: ARG-carrying ORFs per sample, ND and HFD side by side for each antibiotic and timepoint
# (same layout as Fig 4A; the former version pooled all four treatments into one line per diet)
# ═══════════════════════════════════════════════════════
from matplotlib.patches import Patch
tp_order = ['T0h', 'T72h', 'SD14']
ABX_COL = {'None': '#555555', 'AMX': '#E69F00', 'CIP': '#0072B2', 'ERY': '#009E73'}
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
            c = ABX_COL[abx]
            ax.bar(pos + off, vals.mean(), width=0.38, color=c, alpha=(0.9 if diet == 'HFD' else 0.3),
                   edgecolor=c, linewidth=1.2)
            ax.errorbar(pos + off, vals.mean(), yerr=vals.std(), fmt='none', color='black', capsize=2.5, linewidth=0.9)
            ax.scatter(pos + off + rng.uniform(-0.06, 0.06, len(vals)), vals, s=16, zorder=5,
                       marker=('^' if diet == 'HFD' else 'o'), facecolor=(c if diet == 'HFD' else 'white'),
                       edgecolor='black', linewidth=0.6)
        xt.append(pos); xl.append('CTL' if abx == 'None' else abx); pos += 1
    tp_mid.append((start + pos - 1) / 2); pos += 0.7
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=9)
for m, tp in zip(tp_mid, tp_order):
    ax.text(m, -0.12, tp, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=11, fontweight='bold')
ax.set_ylabel("ARG-carrying ORFs per sample\n(assembly-based)", fontsize=10)
ax.set_ylim(0, 560)
from matplotlib.lines import Line2D
ax.legend(handles=[Patch(facecolor=c, edgecolor=c, label=('No antibiotic' if a == 'None' else a)) for a, c in ABX_COL.items()]
          + [Patch(facecolor='#bbbbbb', alpha=0.4, edgecolor='#555555', label='ND (light fill, open circles)'),
             Patch(facecolor='#555555', edgecolor='#555555', label='HFD (solid fill, triangles)')],
          fontsize=8.5, loc='upper left', frameon=False, ncol=6, bbox_to_anchor=(0, 1.02), columnspacing=1.0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_title("Assembly-based ARG-carrying ORFs per sample", fontsize=12, fontweight='bold', pad=14)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig5e_assembly_args_by_diet.pdf")
plt.savefig(f"{OUTDIR}/fig5e_assembly_args_by_diet.png")
plt.close()
print("Fig 5E done")


print("\nAll figures saved to", OUTDIR)
