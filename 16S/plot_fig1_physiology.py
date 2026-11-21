#!/usr/bin/env python3
"""Generate Fig 1: Physiology data from raw experimental measurements.
Data source: ~/vault/kjm/Raw_data_bio-prameter.xlsx"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.rcParams.update({
    'font.size': 9, 'font.family': 'sans-serif',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 300, 'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

OUTDIR = "/home/unnot/vault/kjm/miseq/dada2_output/figures"
DATAFILE = "/home/unnot/vault/kjm/Raw_data_bio-prameter.xlsx"

GROUPS = ['CTL', 'HFD', 'AMX_CTL', 'AMX_HFD', 'CIP_CTL', 'CIP_HFD', 'ERY_CTL', 'ERY_HFD']  # same order as Figs 2-6
ND_GROUPS = ['CTL', 'CIP_CTL', 'AMX_CTL', 'ERY_CTL']
HFD_GROUPS = ['HFD', 'CIP_HFD', 'AMX_HFD', 'ERY_HFD']
colors = {g: '#3498db' if g in ND_GROUPS else '#e74c3c' for g in GROUPS}

# ══════════════════════════════════════════════════
# Load raw data
# ══════════════════════════════════════════════════
wb = pd.ExcelFile(DATAFILE)

# === WEIGHT ===
wt_raw = pd.read_excel(wb, 'Weight', header=None)
# Row 0: dates, Row 1: week labels, Row 2: column headers, Row 3+: data
# Body weight columns: 1,3,5,...,31 (odd indices after col 0)
# Weight gain columns: 2,4,6,...,32

weeks_labels = ['-13W', '-12W', '-11W', '-10W', '-9W', '-8W', '-7W',
                '-6W', '-5W', '-4W', '-3W', '-2W', '-1W', '0d', '7d', '14d']
bw_cols = list(range(1, 33, 2))  # body weight columns (odd indices: 1,3,5,...,31)

# Parse weight data by group
def get_group_from_sample(sample):
    s = str(sample).strip()
    if s.startswith('CTL'): return 'CTL'
    if s.startswith('HFD'): return 'HFD'
    if s.startswith('CIP_CTL'): return 'CIP_CTL'
    if s.startswith('CIP_HFD'): return 'CIP_HFD'
    if s.startswith('AMX_CTL'): return 'AMX_CTL'
    if s.startswith('AMX_HFD'): return 'AMX_HFD'
    if s.startswith('ERY_CTL'): return 'ERY_CTL'
    if s.startswith('ERY_HFD'): return 'ERY_HFD'
    return None

weight_by_group = {g: [] for g in GROUPS}
for i in range(3, wt_raw.shape[0]):
    sample = wt_raw.iloc[i, 0]
    if sample is None or str(sample).strip() == '':
        continue
    grp = get_group_from_sample(sample)
    if grp is None:
        continue
    # Read body weights and compute gain relative to baseline (first timepoint)
    bws = []
    for c in bw_cols:
        v = wt_raw.iloc[i, c]
        if v is None or v == '-' or (isinstance(v, str) and v.strip() == '-'):
            bws.append(np.nan)
        else:
            bws.append(float(v))
    baseline = bws[0]
    gains = [bw - baseline if not np.isnan(bw) else np.nan for bw in bws]
    weight_by_group[grp].append(gains)

# Compute mean and SEM per group per timepoint
weight_mean = {}
weight_sem = {}
for g in GROUPS:
    arr = np.array(weight_by_group[g])
    weight_mean[g] = np.nanmean(arr, axis=0)
    weight_sem[g] = np.nanstd(arr, axis=0, ddof=1) / np.sqrt(np.sum(~np.isnan(arr), axis=0))

# === CHOLESTEROL ===
chol = pd.read_excel(wb, 'Cholesterol', header=1)
chol = chol.dropna(subset=['group'])
# Filter to main 8 groups
chol = chol[chol['group'].isin(GROUPS)]

ldl_mean = {}; ldl_sem = {}
hdl_mean = {}; hdl_sem = {}
for g in GROUPS:
    sub = chol[chol['group'] == g]
    ldl_mean[g] = sub['LDL'].mean()
    ldl_sem[g] = sub['LDL'].std(ddof=1) / np.sqrt(len(sub))
    hdl_mean[g] = sub['HDL'].mean()
    hdl_sem[g] = sub['HDL'].std(ddof=1) / np.sqrt(len(sub))

# === SCFA ===
scfa_raw = pd.read_excel(wb, 'SCFA', header=0)
acetate = scfa_raw[['group', 'value(μM)']].rename(columns={'value(μM)': 'value'})
acetate['scfa'] = 'Acetate'
propionate = scfa_raw[['group.1', 'value(μM).1']].rename(columns={'group.1': 'group', 'value(μM).1': 'value'})
propionate['scfa'] = 'Propionate'
butyrate = scfa_raw[['group.2', 'value(μM).2']].rename(columns={'group.2': 'group', 'value(μM).2': 'value'})
butyrate['scfa'] = 'Butyrate'

all_scfa = pd.concat([acetate, propionate, butyrate], ignore_index=True)
all_scfa = all_scfa.dropna(subset=['group', 'value'])
all_scfa['group_clean'] = all_scfa['group'].str.replace(r'^\d+_', '', regex=True)
all_scfa = all_scfa[all_scfa['group_clean'].isin(GROUPS)]

# Get raw values per group for real boxplots
scfa_data = {}
for s in ['Acetate', 'Propionate', 'Butyrate']:
    scfa_data[s] = {}
    for g in GROUPS:
        vals = all_scfa[(all_scfa['scfa'] == s) & (all_scfa['group_clean'] == g)]['value'].values
        scfa_data[s][g] = vals

# ══════════════════════════════════════════════════
# Plot Fig 1: 6 panels (3x2). Same scheme as Figs 2, 3 and 6: colour = antibiotic (Okabe-Ito);
# ND = dashed line / open circle / light fill, HFD = solid line / filled triangle / solid fill.
# ══════════════════════════════════════════════════
ABX_COL = {"None": "#555555", "AMX": "#E69F00", "CIP": "#0072B2", "ERY": "#009E73"}
def gabx(g): return g.split("_")[0] if "_" in g else "None"
def is_hfd(g): return g.endswith("HFD")
def gcol(g): return ABX_COL[gabx(g)]
TICKLAB = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD"]

fig, axes = plt.subplots(3, 2, figsize=(10, 12.5))  # portrait: A B / C D / E F
x = np.arange(len(GROUPS))

def style_ax(ax, title, ylabel, xt=True):
    if xt:
        ax.set_xticks(x); ax.set_xticklabels(TICKLAB, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel(ylabel); ax.set_title(title, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

def points(ax, k, g, v, rng):
    """Individual mice: ND open circles, HFD filled triangles, in the group colour (darkened edge)."""
    c = gcol(g)
    ax.scatter(k + rng.uniform(-0.14, 0.14, len(v)), v, s=14, marker=("^" if is_hfd(g) else "o"),
               facecolor=(c if is_hfd(g) else "white"), edgecolor="black", linewidth=0.6, zorder=5)

# ── Panel A: Body weight gain ──
ax = axes.flat[0]
x_weeks = np.arange(len(weeks_labels))
for g in GROUPS:
    c = gcol(g); h = is_hfd(g)
    ax.errorbar(x_weeks + (0.08 if h else -0.08), weight_mean[g], yerr=weight_sem[g],
                color=c, mec=c, mfc=(c if h else "white"), marker=("^" if h else "o"), mew=1.0,
                ls=("-" if h else "--"), linewidth=1.3, capsize=2, markersize=4.5)
TK = [0, 3, 6, 9, 13, 14, 15]
ax.set_xticks(TK); ax.set_xticklabels([weeks_labels[i] for i in TK], rotation=45, ha='right', fontsize=8)
ax.axvline(x=13, color='gray', linestyle=':', alpha=0.7, linewidth=0.9)
ax.axvspan(13, 15.5, color="#eeeeee", zorder=0, lw=0)
ax.set_xlim(-0.6, 15.6)
ax.text(6.5, 1.0, "Diet (13 weeks)", transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8, color='gray')
ax.text(14.25, 1.0, "Abx +\nrecovery", transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7.5, color='gray')
style_ax(ax, "A) Body weight gain", "Weight gain (g)", xt=False)
ax.set_title("A) Body weight gain", fontweight='bold', pad=22)

# ── Panels B, C: LDL / HDL cholesterol ──
for i, (col, title, ylab) in enumerate([("LDL", "B) LDL cholesterol", "Serum LDL-C (mg/dL)"),
                                        ("HDL", "C) HDL cholesterol", "Serum HDL-C (mg/dL)")], start=1):
    ax = axes.flat[i]; rng = np.random.default_rng(1)
    m = [chol.loc[chol.group == g, col].mean() for g in GROUPS]
    se = [chol.loc[chol.group == g, col].std(ddof=1) / np.sqrt(chol.loc[chol.group == g, col].notna().sum()) for g in GROUPS]
    for k, g in enumerate(GROUPS):
        ax.bar(k, m[k], color=gcol(g), alpha=(0.9 if is_hfd(g) else 0.3), edgecolor=gcol(g), linewidth=1.2, width=0.75)
        ax.errorbar(k, m[k], yerr=se[k], color="black", capsize=3, lw=1, zorder=4)
        points(ax, k, g, chol.loc[chol.group == g, col].dropna().values, rng)
    style_ax(ax, title, ylab)
axes.flat[2].set_ylim(0, max(60, chol['HDL'].max() * 1.05))

# ── Panels D-F: cecal SCFAs ──
for i, (scfa, title) in enumerate([("Acetate", "D) Acetate"), ("Propionate", "E) Propionate"), ("Butyrate", "F) Butyrate")], start=3):
    ax = axes.flat[i]; rng = np.random.default_rng(2)
    bp = ax.boxplot([scfa_data[scfa][g] for g in GROUPS], positions=x, widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color='black', linewidth=1.5))
    for k, (patch, g) in enumerate(zip(bp['boxes'], GROUPS)):
        patch.set_facecolor(gcol(g)); patch.set_alpha(0.85 if is_hfd(g) else 0.3); patch.set_edgecolor(gcol(g))
        points(ax, k, g, scfa_data[scfa][g], rng)
    style_ax(ax, title, "Concentration (μM)")

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
_h = [Patch(facecolor=c, edgecolor=c, label=("No antibiotic" if a == "None" else a)) for a, c in ABX_COL.items()]
_h += [Line2D([], [], color="black", lw=1.3, ls="--", marker="o", mfc="white", label="ND (light fill)"),
       Line2D([], [], color="black", lw=1.3, ls="-", marker="^", mfc="black", label="HFD (solid fill)")]
plt.tight_layout(h_pad=2.0)
fig.legend(handles=_h, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=6, fontsize=9, frameon=False,
           handlelength=2.6, columnspacing=1.2)
plt.savefig(f"{OUTDIR}/fig1_physiology.pdf")
plt.savefig(f"{OUTDIR}/fig1_physiology.png")
plt.close()
print("Fig 1 done")

# Print summary for manuscript verification
print("\n=== Summary for manuscript ===")
for g in GROUPS:
    print(f"  {g}: final weight gain = {weight_mean[g][-1]:.1f} ± {weight_sem[g][-1]:.1f} g (n={len(weight_by_group[g])})")
print()
for g in GROUPS:
    print(f"  {g}: LDL = {ldl_mean[g]:.1f} ± {ldl_sem[g]:.1f}, HDL = {hdl_mean[g]:.1f} ± {hdl_sem[g]:.1f} mg/dL")
print()
for s in ['Acetate', 'Propionate', 'Butyrate']:
    print(f"  {s}:")
    for g in GROUPS:
        vals = scfa_data[s][g]
        print(f"    {g}: median={np.median(vals):.1f}, mean={np.mean(vals):.1f} ± {np.std(vals,ddof=1):.1f} μM (n={len(vals)})")
