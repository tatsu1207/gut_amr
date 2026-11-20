#!/usr/bin/env python3
"""
Mobilome analysis: Do ARGs on mobile elements (plasmids/phages) persist
differently in HFD vs ND samples across timepoints?

Reads all *_plasmid_summary.tsv files (72 samples)
Outputs summary tables and statistical tests to ~/vault/kjm/hiseq/mobilome_analysis/
"""

import os
import glob
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = os.path.expanduser("~/vault/kjm/hiseq/results")
OUT_DIR = os.path.expanduser("~/vault/kjm/hiseq/mobilome_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# --- 1. Load all plasmid summary files ---
print("=== Loading plasmid summary files ===")
all_rows = []
for f in sorted(glob.glob(f"{RESULTS_DIR}/ARG_*/*_plasmid_summary.tsv")):
    df = pd.read_csv(f, sep='\t')
    if len(df) > 0:
        all_rows.append(df)

df_all = pd.concat(all_rows, ignore_index=True)
print(f"Total ARG entries: {len(df_all)}")
print(f"Unique samples: {df_all['sample'].nunique()}")
print(f"Genomic locations: {df_all['genomic_location'].value_counts().to_dict()}")

# --- 2. Parse sample metadata from sample name ---
def parse_sample(s):
    """ARG_AC1_T0h -> group=AMX_CTL, diet=ND, antibiotic=AMX, timepoint=T0h"""
    s = s.replace("ARG_", "")
    parts = s.rsplit("_", 1)
    prefix, timepoint = parts[0], parts[1]

    group_map = {
        'AC': ('AMX_CTL', 'ND', 'Amoxicillin'),
        'AH': ('AMX_HFD', 'HFD', 'Amoxicillin'),
        'CC': ('CIP_CTL', 'ND', 'Ciprofloxacin'),
        'CH': ('CIP_HFD', 'HFD', 'Ciprofloxacin'),
        'EC': ('ERY_CTL', 'ND', 'Erythromycin'),
        'EH': ('ERY_HFD', 'HFD', 'Erythromycin'),
        'CTL': ('CTL', 'ND', 'None'),
        'HFD': ('HFD', 'HFD', 'None'),
    }

    # Extract group prefix (letters before the digit)
    gp = ''
    for c in prefix:
        if c.isdigit():
            break
        gp += c

    group, diet, abx = group_map.get(gp, ('Unknown', 'Unknown', 'Unknown'))
    return pd.Series({'group': group, 'diet': diet, 'antibiotic': abx, 'timepoint': timepoint})

meta = df_all['sample'].apply(parse_sample)
df_all = pd.concat([df_all, meta], axis=1)

# --- 3. Simplify drug classes ---
def simplify_drug_class(dc):
    """Map detailed drug classes to simple categories"""
    dc_lower = str(dc).lower()
    if 'macrolide' in dc_lower and 'lincosamide' not in dc_lower:
        return 'macrolide'
    elif 'lincosamide' in dc_lower:
        return 'MLS'
    elif 'macrolide' in dc_lower and 'streptogramin' in dc_lower:
        return 'MLS'
    elif 'tetracycline' in dc_lower and 'fluoroquinolone' not in dc_lower:
        return 'tetracycline'
    elif 'beta-lactam' in dc_lower or 'penam' in dc_lower or 'cephalosporin' in dc_lower:
        return 'beta-lactam'
    elif 'glycopeptide' in dc_lower:
        return 'glycopeptide'
    elif 'fluoroquinolone' in dc_lower:
        return 'fluoroquinolone'
    elif 'aminoglycoside' in dc_lower:
        return 'aminoglycoside'
    elif 'nitroimidazole' in dc_lower:
        return 'nitroimidazole'
    elif 'phenicol' in dc_lower:
        return 'phenicol'
    elif 'rifamycin' in dc_lower:
        return 'rifamycin'
    elif 'peptide' in dc_lower:
        return 'peptide'
    elif 'diaminopyrimidine' in dc_lower:
        return 'diaminopyrimidine'
    elif 'sulfonamide' in dc_lower:
        return 'sulfonamide'
    else:
        # Multi-drug or complex: check for tet + fluoro combo
        if 'tetracycline' in dc_lower and 'fluoroquinolone' in dc_lower:
            return 'multidrug_efflux'
        return 'other'

df_all['simple_class'] = df_all['drug_class'].apply(simplify_drug_class)

# --- 4. Summary stats ---
print("\n=== ARGs by genomic location ===")
loc_counts = df_all.groupby('genomic_location').size()
print(loc_counts)

print("\n=== Plasmid-borne ARGs by drug class ===")
plasmid_df = df_all[df_all['genomic_location'] == 'plasmid']
print(plasmid_df['simple_class'].value_counts())

print("\n=== Phage-borne ARGs by drug class ===")
phage_df = df_all[df_all['genomic_location'] == 'phage']
if len(phage_df) > 0:
    print(phage_df['simple_class'].value_counts())
else:
    print("None found")

# --- 5. Per-sample counts by location ---
print("\n=== Per-sample ARG counts by genomic location ===")

# Count ARGs per sample per location
sample_loc = df_all.groupby(['sample', 'genomic_location']).size().unstack(fill_value=0)
sample_loc = sample_loc.reset_index()

# Add metadata
for col in ['group', 'diet', 'antibiotic', 'timepoint']:
    sample_loc[col] = sample_loc['sample'].map(
        df_all.drop_duplicates('sample').set_index('sample')[col]
    )

# Calculate plasmid fraction
if 'plasmid' in sample_loc.columns:
    total_cols = [c for c in sample_loc.columns if c in ['chromosome', 'plasmid', 'phage']]
    sample_loc['total_ARGs'] = sample_loc[total_cols].sum(axis=1)
    sample_loc['plasmid_fraction'] = sample_loc.get('plasmid', 0) / sample_loc['total_ARGs']
    sample_loc['plasmid_fraction'] = sample_loc['plasmid_fraction'].fillna(0)

sample_loc.to_csv(f"{OUT_DIR}/per_sample_arg_location.tsv", sep='\t', index=False)
print(f"Saved: {OUT_DIR}/per_sample_arg_location.tsv")

# --- 6. HFD vs ND comparisons at each timepoint ---
print("\n" + "="*60)
print("=== HFD vs ND: Plasmid-borne ARG counts ===")
print("="*60)

# Only antibiotic-treated samples (exclude CTL/HFD-only groups)
abx_samples = sample_loc[sample_loc['antibiotic'].isin(['Amoxicillin', 'Ciprofloxacin', 'Erythromycin'])]

results_rows = []
for tp in ['T0h', 'T72h', 'SD14']:
    tp_data = abx_samples[abx_samples['timepoint'] == tp]
    if len(tp_data) == 0:
        continue

    hfd = tp_data[tp_data['diet'] == 'HFD']
    nd = tp_data[tp_data['diet'] == 'ND']

    if 'plasmid' in tp_data.columns:
        hfd_plasmid = hfd['plasmid'].values if 'plasmid' in hfd.columns else np.array([0])
        nd_plasmid = nd['plasmid'].values if 'plasmid' in nd.columns else np.array([0])
    else:
        hfd_plasmid = np.array([0])
        nd_plasmid = np.array([0])

    if len(hfd_plasmid) >= 2 and len(nd_plasmid) >= 2:
        stat, pval = stats.mannwhitneyu(hfd_plasmid, nd_plasmid, alternative='two-sided')
    else:
        stat, pval = np.nan, np.nan

    print(f"\n{tp}:")
    print(f"  HFD (n={len(hfd)}): plasmid ARGs = {hfd_plasmid.mean():.1f} +/- {hfd_plasmid.std():.1f}")
    print(f"  ND  (n={len(nd)}):  plasmid ARGs = {nd_plasmid.mean():.1f} +/- {nd_plasmid.std():.1f}")
    print(f"  Mann-Whitney p = {pval:.4f}" if not np.isnan(pval) else "  Not enough samples")

    results_rows.append({
        'timepoint': tp,
        'HFD_mean': hfd_plasmid.mean(), 'HFD_sd': hfd_plasmid.std(), 'HFD_n': len(hfd),
        'ND_mean': nd_plasmid.mean(), 'ND_sd': nd_plasmid.std(), 'ND_n': len(nd),
        'p_value': pval
    })

    # Also by antibiotic
    for abx in ['Amoxicillin', 'Ciprofloxacin', 'Erythromycin']:
        abx_tp = tp_data[tp_data['antibiotic'] == abx]
        h = abx_tp[abx_tp['diet'] == 'HFD']['plasmid'].values if 'plasmid' in abx_tp.columns else np.array([0])
        n = abx_tp[abx_tp['diet'] == 'ND']['plasmid'].values if 'plasmid' in abx_tp.columns else np.array([0])
        if len(h) >= 2 and len(n) >= 2:
            _, p = stats.mannwhitneyu(h, n, alternative='two-sided')
        else:
            p = np.nan
        print(f"  {abx}: HFD={h.mean():.1f}+/-{h.std():.1f} vs ND={n.mean():.1f}+/-{n.std():.1f}, p={p:.3f}" if not np.isnan(p) else f"  {abx}: HFD={h.mean():.1f} vs ND={n.mean():.1f} (n too small)")

pd.DataFrame(results_rows).to_csv(f"{OUT_DIR}/hfd_vs_nd_plasmid_args.tsv", sep='\t', index=False)

# --- 7. Which specific ARGs are on plasmids? ---
print("\n" + "="*60)
print("=== Top plasmid-borne ARGs ===")
print("="*60)

plasmid_args = df_all[df_all['genomic_location'] == 'plasmid']
top_plasmid = plasmid_args.groupby('best_hit_ARG').agg(
    count=('sample', 'size'),
    n_samples=('sample', 'nunique'),
    drug_classes=('simple_class', 'first')
).sort_values('count', ascending=False)
print(top_plasmid.head(30).to_string())
top_plasmid.to_csv(f"{OUT_DIR}/top_plasmid_args.tsv", sep='\t')

# --- 8. Plasmid ARGs by drug class and diet across timepoints ---
print("\n" + "="*60)
print("=== Plasmid ARGs by drug class: HFD vs ND at each timepoint ===")
print("="*60)

# Count plasmid ARGs per sample per drug class
plasmid_sample_class = plasmid_args.groupby(['sample', 'simple_class']).size().unstack(fill_value=0)
plasmid_sample_class = plasmid_sample_class.reset_index()

# Add metadata
for col in ['group', 'diet', 'antibiotic', 'timepoint']:
    plasmid_sample_class[col] = plasmid_sample_class['sample'].map(
        df_all.drop_duplicates('sample').set_index('sample')[col]
    )

drug_classes = [c for c in plasmid_sample_class.columns if c not in ['sample', 'group', 'diet', 'antibiotic', 'timepoint']]

class_results = []
for tp in ['T0h', 'T72h', 'SD14']:
    print(f"\n--- {tp} ---")
    tp_data = plasmid_sample_class[
        (plasmid_sample_class['timepoint'] == tp) &
        (plasmid_sample_class['antibiotic'].isin(['Amoxicillin', 'Ciprofloxacin', 'Erythromycin']))
    ]
    for dc in sorted(drug_classes):
        if dc in tp_data.columns:
            hfd = tp_data[tp_data['diet'] == 'HFD'][dc].values
            nd = tp_data[tp_data['diet'] == 'ND'][dc].values
            if hfd.sum() + nd.sum() > 0:
                if len(hfd) >= 2 and len(nd) >= 2:
                    _, p = stats.mannwhitneyu(hfd, nd, alternative='two-sided')
                else:
                    p = np.nan
                print(f"  {dc:20s}: HFD={hfd.mean():.2f}+/-{hfd.std():.2f} vs ND={nd.mean():.2f}+/-{nd.std():.2f}  p={p:.4f}" if not np.isnan(p) else f"  {dc:20s}: HFD={hfd.mean():.2f} vs ND={nd.mean():.2f}")
                class_results.append({
                    'timepoint': tp, 'drug_class': dc,
                    'HFD_mean': hfd.mean(), 'HFD_sd': hfd.std(),
                    'ND_mean': nd.mean(), 'ND_sd': nd.std(),
                    'p_value': p
                })

pd.DataFrame(class_results).to_csv(f"{OUT_DIR}/plasmid_args_by_class_diet.tsv", sep='\t', index=False)

# --- 9. Specific ARGs on plasmids: HFD vs ND at SD14 ---
print("\n" + "="*60)
print("=== SD14: Which plasmid ARGs persist more in HFD? ===")
print("="*60)

sd14_plasmid = plasmid_args[
    (plasmid_args['timepoint'] == 'SD14') &
    (plasmid_args['antibiotic'].isin(['Amoxicillin', 'Ciprofloxacin', 'Erythromycin']))
]

# Count per sample per ARG
sd14_arg_sample = sd14_plasmid.groupby(['sample', 'best_hit_ARG', 'diet']).size().reset_index(name='count')

# For each ARG found on plasmids at SD14, compare HFD vs ND
sd14_args = sd14_plasmid['best_hit_ARG'].value_counts()
print(f"\nTotal unique plasmid ARGs at SD14: {len(sd14_args)}")
print(f"\nARGs found on plasmids at SD14 (sorted by frequency):")

for arg in sd14_args.index:
    arg_data = sd14_plasmid[sd14_plasmid['best_hit_ARG'] == arg]
    hfd_samples = arg_data[arg_data['diet'] == 'HFD']['sample'].nunique()
    nd_samples = arg_data[arg_data['diet'] == 'ND']['sample'].nunique()
    dc = arg_data['simple_class'].iloc[0]
    print(f"  {arg:45s} [{dc:15s}]: HFD={hfd_samples} samples, ND={nd_samples} samples")

# --- 10. Mobilome dynamics: total plasmid contigs from geNomad ---
print("\n" + "="*60)
print("=== geNomad plasmid contig counts per sample ===")
print("="*60)

genomad_counts = []
for d in sorted(glob.glob(f"{RESULTS_DIR}/ARG_*")):
    sample = os.path.basename(d)
    plasmid_fna = f"{d}/genomad_output/final.contigs_summary/final.contigs_plasmid_summary.tsv"
    if os.path.exists(plasmid_fna):
        try:
            gdf = pd.read_csv(plasmid_fna, sep='\t')
            n_plasmids = len(gdf)
            conj_genes = gdf['conjugation_genes'].notna().sum() if 'conjugation_genes' in gdf.columns else 0
            has_conj = (gdf['conjugation_genes'].str.len() > 0).sum() if 'conjugation_genes' in gdf.columns else 0
        except:
            n_plasmids = 0
            has_conj = 0
        meta = parse_sample(sample)
        genomad_counts.append({
            'sample': sample, 'n_plasmid_contigs': n_plasmids,
            'n_conjugative': has_conj,
            **meta.to_dict()
        })

genomad_df = pd.DataFrame(genomad_counts)
genomad_df.to_csv(f"{OUT_DIR}/genomad_plasmid_counts.tsv", sep='\t', index=False)

# HFD vs ND plasmid contig counts
for tp in ['T0h', 'T72h', 'SD14']:
    tp_data = genomad_df[
        (genomad_df['timepoint'] == tp) &
        (genomad_df['antibiotic'].isin(['Amoxicillin', 'Ciprofloxacin', 'Erythromycin']))
    ]
    hfd = tp_data[tp_data['diet'] == 'HFD']['n_plasmid_contigs'].values
    nd = tp_data[tp_data['diet'] == 'ND']['n_plasmid_contigs'].values
    if len(hfd) >= 2 and len(nd) >= 2:
        _, p = stats.mannwhitneyu(hfd, nd, alternative='two-sided')
    else:
        p = np.nan
    print(f"\n{tp}: Total plasmid contigs")
    print(f"  HFD: {hfd.mean():.0f} +/- {hfd.std():.0f} (n={len(hfd)})")
    print(f"  ND:  {nd.mean():.0f} +/- {nd.std():.0f} (n={len(nd)})")
    print(f"  p = {p:.4f}" if not np.isnan(p) else "  n too small")

# --- 11. Conjugative plasmid analysis ---
print("\n" + "="*60)
print("=== Conjugative plasmid analysis ===")
print("="*60)

for tp in ['T0h', 'T72h', 'SD14']:
    tp_data = genomad_df[
        (genomad_df['timepoint'] == tp) &
        (genomad_df['antibiotic'].isin(['Amoxicillin', 'Ciprofloxacin', 'Erythromycin']))
    ]
    hfd = tp_data[tp_data['diet'] == 'HFD']['n_conjugative'].values
    nd = tp_data[tp_data['diet'] == 'ND']['n_conjugative'].values
    if len(hfd) >= 2 and len(nd) >= 2:
        _, p = stats.mannwhitneyu(hfd, nd, alternative='two-sided')
    else:
        p = np.nan
    print(f"\n{tp}: Conjugative plasmids")
    print(f"  HFD: {hfd.mean():.1f} +/- {hfd.std():.1f}")
    print(f"  ND:  {nd.mean():.1f} +/- {nd.std():.1f}")
    print(f"  p = {p:.4f}" if not np.isnan(p) else "  n too small")

# --- 12. Summary report ---
print("\n" + "="*60)
print("=== SUMMARY ===")
print("="*60)

total = len(df_all)
n_plasmid = len(df_all[df_all['genomic_location'] == 'plasmid'])
n_phage = len(df_all[df_all['genomic_location'] == 'phage'])
n_chrom = len(df_all[df_all['genomic_location'] == 'chromosome'])

print(f"Total assembly-based ARGs: {total}")
print(f"  Chromosome: {n_chrom} ({100*n_chrom/total:.1f}%)")
print(f"  Plasmid: {n_plasmid} ({100*n_plasmid/total:.1f}%)")
print(f"  Phage: {n_phage} ({100*n_phage/total:.1f}%)")

print(f"\nOutput files saved to: {OUT_DIR}/")
print("  per_sample_arg_location.tsv")
print("  hfd_vs_nd_plasmid_args.tsv")
print("  top_plasmid_args.tsv")
print("  plasmid_args_by_class_diet.tsv")
print("  genomad_plasmid_counts.tsv")

print("\n=== DONE ===")
