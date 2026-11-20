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

GROUPS = ['CTL', 'HFD', 'CIP_CTL', 'CIP_HFD', 'AMX_CTL', 'AMX_HFD', 'ERY_CTL', 'ERY_HFD']
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
# Plot Fig 1: 6 panels (2x3)
# ══════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(14, 9))

# ── Panel A: Body weight gain ──
ax = axes[0, 0]
markers = {'CTL': 'o', 'HFD': 'o', 'CIP_CTL': 'v', 'CIP_HFD': 'v',
           'AMX_CTL': 's', 'AMX_HFD': 's', 'ERY_CTL': 'D', 'ERY_HFD': 'D'}
fillstyle = {g: 'full' if g in ND_GROUPS else 'none' for g in GROUPS}
x_weeks = list(range(len(weeks_labels)))

for g in GROUPS:
    ax.errorbar(x_weeks, weight_mean[g], yerr=weight_sem[g],
                marker=markers[g], color=colors[g], label=g,
                linewidth=1.2, capsize=2, markersize=4,
                fillstyle=fillstyle[g], markeredgewidth=0.8)

ax.set_xticks(x_weeks[::2])
ax.set_xticklabels([weeks_labels[i] for i in range(0, len(weeks_labels), 2)], rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Weight gain (g)")
ax.set_title("A) Body weight gain", fontweight='bold')
ax.axvline(x=13, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.text(6, max(weight_mean['AMX_HFD']) + 2, 'Diet period', ha='center', fontsize=7, color='gray')
ax.text(14.5, max(weight_mean['AMX_HFD']) + 2, 'Abx', ha='center', fontsize=7, color='gray')
ax.legend(fontsize=5.5, ncol=2, loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel B: LDL cholesterol ──
ax = axes[0, 1]
x = np.arange(len(GROUPS))
bar_colors = [colors[g] for g in GROUPS]
bars = ax.bar(x, [ldl_mean[g] for g in GROUPS], yerr=[ldl_sem[g] for g in GROUPS],
              color=bar_colors, edgecolor='white', capsize=3, alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Serum LDL-C (mg/dL)")
ax.set_title("B) LDL cholesterol", fontweight='bold')

# npj guideline: bar charts must overlay individual data points
_rng = np.random.default_rng(1)
for _k, _g in enumerate(GROUPS):
    _v = chol.loc[chol['group'] == _g, 'LDL'].dropna().values
    ax.scatter(_k + _rng.uniform(-0.15, 0.15, len(_v)), _v, s=8, color='black', alpha=0.7, zorder=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel C: HDL cholesterol ──
ax = axes[0, 2]
bars = ax.bar(x, [hdl_mean[g] for g in GROUPS], yerr=[hdl_sem[g] for g in GROUPS],
              color=bar_colors, edgecolor='white', capsize=3, alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Serum HDL-C (mg/dL)")
ax.set_title("C) HDL cholesterol", fontweight='bold')

# npj guideline: bar charts must overlay individual data points
_rng = np.random.default_rng(1)
for _k, _g in enumerate(GROUPS):
    _v = chol.loc[chol['group'] == _g, 'HDL'].dropna().values
    ax.scatter(_k + _rng.uniform(-0.15, 0.15, len(_v)), _v, s=8, color='black', alpha=0.7, zorder=5)
ax.set_ylim(0, max(60, chol['HDL'].max() * 1.05))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel D: Acetate ──
ax = axes[1, 0]
bp = ax.boxplot([scfa_data['Acetate'][g] for g in GROUPS],
                positions=x, widths=0.6, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.5))
for patch, g in zip(bp['boxes'], GROUPS):
    patch.set_facecolor(colors[g])
    patch.set_alpha(0.7)
ax.set_xticks(x)
ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Concentration (μM)")
ax.set_title("D) Acetate", fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel E: Propionate ──
ax = axes[1, 1]
bp = ax.boxplot([scfa_data['Propionate'][g] for g in GROUPS],
                positions=x, widths=0.6, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.5))
for patch, g in zip(bp['boxes'], GROUPS):
    patch.set_facecolor(colors[g])
    patch.set_alpha(0.7)
ax.set_xticks(x)
ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Concentration (μM)")
ax.set_title("E) Propionate", fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel F: Butyrate ──
ax = axes[1, 2]
bp = ax.boxplot([scfa_data['Butyrate'][g] for g in GROUPS],
                positions=x, widths=0.6, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.5))
for patch, g in zip(bp['boxes'], GROUPS):
    patch.set_facecolor(colors[g])
    patch.set_alpha(0.7)
ax.set_xticks(x)
ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=7)
ax.set_ylabel("Concentration (μM)")
ax.set_title("F) Butyrate", fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
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
