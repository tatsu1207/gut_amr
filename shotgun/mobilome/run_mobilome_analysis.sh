#!/bin/bash
# Mobilome-ARG association analysis
# Links ARGs to plasmids (geNomad), conjugation machinery, transposons/IS elements, and phages
# Safe to run in background with nohup
set -euo pipefail

RESDIR=~/vault/kjm/hiseq/results
OUTDIR=~/vault/kjm/hiseq/mobilome_analysis
mkdir -p "$OUTDIR"
exec > >(tee -a "$OUTDIR/mobilome.log") 2>&1

timestamp() { date "+[%Y-%m-%d %H:%M:%S]"; }
echo "$(timestamp) === Starting mobilome-ARG analysis ==="

python3 <<'PYEOF'
import os, csv, glob, json, sys
from collections import defaultdict

RESDIR = os.path.expanduser("~/vault/kjm/hiseq/results")
OUTDIR = os.path.expanduser("~/vault/kjm/hiseq/mobilome_analysis")

def parse_sample(s):
    """Parse sample name -> group, diet, timepoint"""
    parts = s.replace('ARG_', '').split('_')
    prefix = parts[0][:2]
    tp = parts[1]
    group_map = {'AC':'AMX_CTL','AH':'AMX_HFD','CC':'CIP_CTL','CH':'CIP_HFD',
                 'EC':'ERY_CTL','EH':'ERY_HFD','CT':'CTL','HF':'HFD'}
    group = group_map.get(prefix, prefix)
    diet = 'HFD' if 'HFD' in group or group == 'HFD' else 'ND'
    abx = group.split('_')[0] if '_' in group else 'none'
    return group, diet, abx, tp

samples = sorted([d for d in os.listdir(RESDIR) if d.startswith('ARG_')])
print(f"Processing {len(samples)} samples\n")

###########################################################################
# 1. Load all ARGs from RGI main results
###########################################################################
print("Step 1: Loading ARGs from RGI results...")
all_args = []  # list of dicts
for sample in samples:
    rgi_file = os.path.join(RESDIR, sample, f"{sample}_rgi_main.txt")
    if not os.path.isfile(rgi_file):
        print(f"  WARNING: no RGI for {sample}")
        continue
    with open(rgi_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Parse contig from ORF_ID (format: k141_12345_6 # ...)
            orf_full = row['ORF_ID']
            orf_id = orf_full.split(' # ')[0].strip()
            # Contig = everything up to the last underscore
            contig = '_'.join(orf_id.split('_')[:-1])
            all_args.append({
                'sample': sample,
                'contig': contig,
                'orf_id': orf_id,
                'best_hit_ARG': row['Best_Hit_ARO'],
                'drug_class': row['Drug Class'],
                'resistance_mechanism': row['Resistance Mechanism'],
                'amr_gene_family': row['AMR Gene Family'],
                'cut_off': row['Cut_Off'],
            })
print(f"  Loaded {len(all_args)} ARGs from {len(samples)} samples")

###########################################################################
# 2. Load geNomad plasmid contigs
###########################################################################
print("\nStep 2: Loading geNomad plasmid classifications...")
plasmid_contigs = {}  # (sample, contig) -> {score, conjugation_genes, ...}
for sample in samples:
    pf = os.path.join(RESDIR, sample, "genomad_output",
                      "final.contigs_summary", "final.contigs_plasmid_summary.tsv")
    if not os.path.isfile(pf):
        continue
    with open(pf) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            plasmid_contigs[(sample, row['seq_name'])] = {
                'plasmid_score': float(row['plasmid_score']),
                'topology': row['topology'],
                'conjugation_genes': row['conjugation_genes'],
                'length': int(row['length']),
            }
print(f"  Total plasmid contigs: {len(plasmid_contigs)}")

###########################################################################
# 3. Load geNomad virus/phage contigs
###########################################################################
print("\nStep 3: Loading geNomad virus/phage classifications...")
virus_contigs = {}  # (sample, contig) -> {score, taxonomy}
for sample in samples:
    vf = os.path.join(RESDIR, sample, "genomad_output",
                      "final.contigs_summary", "final.contigs_virus_summary.tsv")
    if not os.path.isfile(vf):
        continue
    with open(vf) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Virus contigs may have provirus coordinates; use seq_name
            virus_contigs[(sample, row['seq_name'])] = {
                'virus_score': float(row['virus_score']),
                'taxonomy': row['taxonomy'],
                'length': int(row['length']),
            }
print(f"  Total virus/phage contigs: {len(virus_contigs)}")

###########################################################################
# 4. Load geNomad gene annotations for transposase/IS/integrase detection
###########################################################################
print("\nStep 4: Scanning geNomad gene annotations for transposases/IS/integrases...")
# Build a set of contigs that carry mobile element genes
# Also track which ORFs are mobile element genes for co-location with ARGs
mobile_genes = {}  # (sample, contig) -> list of {gene, description, type}
mobile_gene_set = set()  # (sample, gene_id) for quick lookup

import re
mobile_patterns = [
    (re.compile(r'transpos', re.I), 'transposase'),
    (re.compile(r'IS\d+|IS[A-Z][a-z]\d+|insertion.sequence', re.I), 'IS_element'),
    (re.compile(r'integrase|intI|int[A-Z]\b', re.I), 'integrase'),
    (re.compile(r'recombinase|invertase|resolvase', re.I), 'recombinase'),
    (re.compile(r'conjugat|tra[A-Z]\b|trb[A-Z]\b|mob[A-Z]\b|T4SS|type.IV.secretion', re.I), 'conjugation'),
    (re.compile(r'relaxase', re.I), 'relaxase'),
]

for sample in samples:
    gf = os.path.join(RESDIR, sample, "genomad_output",
                      "final.contigs_annotate", "final.contigs_genes.tsv")
    if not os.path.isfile(gf):
        continue
    with open(gf) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            desc = row.get('annotation_description', '') or ''
            conjscan = row.get('annotation_conjscan', '') or ''
            gene_id = row['gene']
            contig = '_'.join(gene_id.split('_')[:-1])

            mge_types = []
            for pattern, mge_type in mobile_patterns:
                if pattern.search(desc) or pattern.search(conjscan):
                    mge_types.append(mge_type)

            # Also flag conjugation from conjscan column
            if conjscan and conjscan != 'NA':
                if 'conjugation' not in mge_types:
                    mge_types.append('conjugation')

            if mge_types:
                key = (sample, contig)
                if key not in mobile_genes:
                    mobile_genes[key] = []
                mobile_genes[key].append({
                    'gene': gene_id,
                    'description': desc[:100],
                    'types': mge_types,
                })
                mobile_gene_set.add((sample, gene_id))

n_contigs_with_mge = len(mobile_genes)
print(f"  Contigs with mobile element genes: {n_contigs_with_mge}")

###########################################################################
# 5. Load plasmid_summary.tsv (pre-computed plasmid/chromosome classification)
###########################################################################
print("\nStep 5: Loading plasmid_summary classifications...")
plasmid_summary = {}  # (sample, contig) -> genomic_location
for sample in samples:
    psf = os.path.join(RESDIR, sample, f"{sample}_plasmid_summary.tsv")
    if not os.path.isfile(psf):
        continue
    with open(psf) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            plasmid_summary[(sample, row['contig'])] = row['genomic_location']

###########################################################################
# 6. Classify each ARG by mobilome association
###########################################################################
print("\nStep 6: Classifying ARG mobilome associations...")

results = []
for a in all_args:
    s, c = a['sample'], a['contig']
    group, diet, abx, tp = parse_sample(s)

    # Mobilome flags
    on_plasmid = (s, c) in plasmid_contigs
    on_phage = (s, c) in virus_contigs
    near_mge = (s, c) in mobile_genes
    ps_class = plasmid_summary.get((s, c), 'unknown')

    # Determine mobilome category
    categories = []
    if on_plasmid:
        pinfo = plasmid_contigs[(s, c)]
        categories.append('plasmid')
        if pinfo['conjugation_genes'] and pinfo['conjugation_genes'] != 'NA':
            categories.append('conjugative_plasmid')
    if on_phage:
        categories.append('phage')
    if near_mge:
        mge_info = mobile_genes[(s, c)]
        mge_types = set()
        for g in mge_info:
            mge_types.update(g['types'])
        categories.extend(sorted(mge_types))
    if ps_class == 'plasmid' and 'plasmid' not in categories:
        categories.append('plasmid_mobtyper')

    if not categories:
        mobilome = 'chromosomal'
    else:
        mobilome = ';'.join(sorted(set(categories)))

    # Conjugation info
    conj_genes = ''
    if on_plasmid:
        cg = plasmid_contigs[(s, c)].get('conjugation_genes', 'NA')
        if cg and cg != 'NA':
            conj_genes = cg

    # MGE gene details
    mge_detail = ''
    if near_mge:
        mge_info = mobile_genes[(s, c)]
        details = []
        for g in mge_info:
            details.append(f"{g['gene']}({','.join(g['types'])})")
        mge_detail = '; '.join(details[:5])  # limit to 5

    results.append({
        **a,
        'group': group,
        'diet': diet,
        'antibiotic': abx,
        'timepoint': tp,
        'mobilome': mobilome,
        'on_plasmid_genomad': on_plasmid,
        'on_phage': on_phage,
        'near_mge_gene': near_mge,
        'plasmid_summary_class': ps_class,
        'conjugation_genes': conj_genes,
        'mge_genes_on_contig': mge_detail,
    })

###########################################################################
# 7. Write full results table
###########################################################################
print("\nStep 7: Writing output files...")

outf = os.path.join(OUTDIR, 'arg_mobilome_full.tsv')
cols = ['sample','group','diet','antibiotic','timepoint','contig','orf_id',
        'best_hit_ARG','drug_class','resistance_mechanism','amr_gene_family',
        'cut_off','mobilome','on_plasmid_genomad','on_phage','near_mge_gene',
        'plasmid_summary_class','conjugation_genes','mge_genes_on_contig']
with open(outf, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in results:
        f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
print(f"  Full table: {outf} ({len(results)} rows)")

###########################################################################
# 8. Summary statistics
###########################################################################
print("\n" + "="*70)
print("=== MOBILOME-ARG SUMMARY ===")
print("="*70)

# Overall mobilome categories
cat_counts = defaultdict(int)
for r in results:
    cat_counts[r['mobilome']] += 1
print("\n--- ARGs by mobilome category ---")
for cat, count in sorted(cat_counts.items(), key=lambda x: -x[1]):
    pct = count * 100 / len(results)
    print(f"  {cat}: {count} ({pct:.1f}%)")

# Plasmid-borne ARGs
plasmid_args = [r for r in results if 'plasmid' in r['mobilome']]
print(f"\n--- Plasmid-associated ARGs: {len(plasmid_args)} ---")
if plasmid_args:
    drug_counts = defaultdict(int)
    for r in plasmid_args:
        drug_counts[r['drug_class']] += 1
    print("  By drug class:")
    for drug, count in sorted(drug_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"    {drug}: {count}")

    diet_counts = defaultdict(int)
    for r in plasmid_args:
        diet_counts[r['diet']] += 1
    print(f"  By diet: ND={diet_counts.get('ND',0)}, HFD={diet_counts.get('HFD',0)}")

# Conjugative plasmid ARGs
conj_args = [r for r in results if 'conjugative_plasmid' in r['mobilome']]
print(f"\n--- Conjugative plasmid ARGs: {len(conj_args)} ---")
if conj_args:
    drug_counts = defaultdict(int)
    for r in conj_args:
        drug_counts[r['drug_class']] += 1
    print("  By drug class:")
    for drug, count in sorted(drug_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    {drug}: {count}")

# Transposase/IS-associated ARGs
is_args = [r for r in results if any(x in r['mobilome'] for x in ['transposase','IS_element','integrase','recombinase'])]
print(f"\n--- Transposase/IS/integrase-associated ARGs: {len(is_args)} ---")
if is_args:
    drug_counts = defaultdict(int)
    for r in is_args:
        drug_counts[r['drug_class']] += 1
    print("  By drug class:")
    for drug, count in sorted(drug_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    {drug}: {count}")

# Phage-associated ARGs
phage_args = [r for r in results if 'phage' in r['mobilome']]
print(f"\n--- Phage-associated ARGs: {len(phage_args)} ---")
if phage_args:
    drug_counts = defaultdict(int)
    for r in phage_args:
        drug_counts[r['drug_class']] += 1
    print("  By drug class:")
    for drug, count in sorted(drug_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    {drug}: {count}")

###########################################################################
# 9. Diet comparison: mobile ARGs in HFD vs ND
###########################################################################
print("\n" + "="*70)
print("=== DIET COMPARISON: Mobile ARGs ===")
print("="*70)

# By diet x timepoint: fraction of ARGs on mobile elements
for tp in ['T0h', 'T72h', 'SD14']:
    print(f"\n--- Timepoint: {tp} ---")
    for diet in ['ND', 'HFD']:
        tp_diet_args = [r for r in results if r['timepoint'] == tp and r['diet'] == diet]
        if not tp_diet_args:
            continue
        mobile = sum(1 for r in tp_diet_args if r['mobilome'] != 'chromosomal')
        total = len(tp_diet_args)
        pct = mobile * 100 / total if total else 0
        plasmid_n = sum(1 for r in tp_diet_args if 'plasmid' in r['mobilome'])
        is_n = sum(1 for r in tp_diet_args if any(x in r['mobilome'] for x in ['transposase','IS_element','integrase']))
        phage_n = sum(1 for r in tp_diet_args if 'phage' in r['mobilome'])
        print(f"  {diet}: {mobile}/{total} mobile ({pct:.1f}%) — plasmid:{plasmid_n}, IS/transposase:{is_n}, phage:{phage_n}")

###########################################################################
# 10. Drug class x mobilome summary table
###########################################################################
print("\n" + "="*70)
print("=== DRUG CLASS × MOBILOME ===")
print("="*70)

drug_mob = defaultdict(lambda: defaultdict(int))
drug_total = defaultdict(int)
for r in results:
    drug = r['drug_class']
    drug_total[drug] += 1
    if r['mobilome'] != 'chromosomal':
        drug_mob[drug]['mobile'] += 1
    if 'plasmid' in r['mobilome']:
        drug_mob[drug]['plasmid'] += 1
    if any(x in r['mobilome'] for x in ['transposase','IS_element','integrase','recombinase']):
        drug_mob[drug]['IS_transposase'] += 1
    if 'phage' in r['mobilome']:
        drug_mob[drug]['phage'] += 1

outf2 = os.path.join(OUTDIR, 'drug_class_mobilome_summary.tsv')
with open(outf2, 'w') as f:
    f.write('drug_class\ttotal_ARGs\tmobile\tmobile_pct\tplasmid\tIS_transposase\tphage\n')
    for drug, total in sorted(drug_total.items(), key=lambda x: -x[1]):
        mob = drug_mob[drug]
        mobile = mob.get('mobile', 0)
        pct = mobile * 100 / total if total else 0
        f.write(f"{drug}\t{total}\t{mobile}\t{pct:.1f}\t{mob.get('plasmid',0)}\t{mob.get('IS_transposase',0)}\t{mob.get('phage',0)}\n")
        if total >= 50:  # only print major classes
            print(f"  {drug}: {mobile}/{total} mobile ({pct:.1f}%) — plasmid:{mob.get('plasmid',0)}, IS:{mob.get('IS_transposase',0)}, phage:{mob.get('phage',0)}")
print(f"\n  Full table: {outf2}")

###########################################################################
# 11. Top mobile ARGs (specific genes)
###########################################################################
print("\n" + "="*70)
print("=== TOP 30 MOBILE ARGs ===")
print("="*70)

mobile_arg_counts = defaultdict(lambda: {'total': 0, 'plasmid': 0, 'IS': 0, 'phage': 0, 'drugs': set()})
for r in results:
    if r['mobilome'] == 'chromosomal':
        continue
    gene = r['best_hit_ARG']
    mobile_arg_counts[gene]['total'] += 1
    mobile_arg_counts[gene]['drugs'].add(r['drug_class'])
    if 'plasmid' in r['mobilome']:
        mobile_arg_counts[gene]['plasmid'] += 1
    if any(x in r['mobilome'] for x in ['transposase','IS_element','integrase']):
        mobile_arg_counts[gene]['IS'] += 1
    if 'phage' in r['mobilome']:
        mobile_arg_counts[gene]['phage'] += 1

for gene, info in sorted(mobile_arg_counts.items(), key=lambda x: -x[1]['total'])[:30]:
    drugs = '; '.join(sorted(info['drugs']))[:60]
    print(f"  {gene}: {info['total']} (plasmid:{info['plasmid']}, IS:{info['IS']}, phage:{info['phage']}) — {drugs}")

###########################################################################
# 12. MLS ARGs on mobile elements (key for manuscript)
###########################################################################
print("\n" + "="*70)
print("=== MLS ARGs ON MOBILE ELEMENTS (key finding) ===")
print("="*70)

mls_keywords = ['macrolide', 'lincosamide', 'streptogramin', 'MLS']
mls_args = [r for r in results if any(k.lower() in r['drug_class'].lower() for k in mls_keywords)]
mls_mobile = [r for r in mls_args if r['mobilome'] != 'chromosomal']
print(f"Total MLS ARGs: {len(mls_args)}")
print(f"Mobile MLS ARGs: {len(mls_mobile)} ({len(mls_mobile)*100/len(mls_args):.1f}%)" if mls_args else "")

if mls_mobile:
    print("\nMobile MLS ARGs by gene:")
    gene_counts = defaultdict(int)
    for r in mls_mobile:
        gene_counts[r['best_hit_ARG']] += 1
    for gene, count in sorted(gene_counts.items(), key=lambda x: -x[1]):
        print(f"  {gene}: {count}")

    print("\nMobile MLS ARGs by diet:")
    for diet in ['ND', 'HFD']:
        d_args = [r for r in mls_mobile if r['diet'] == diet]
        print(f"  {diet}: {len(d_args)}")

    print("\nMobile MLS ARGs by timepoint × diet:")
    for tp in ['T0h', 'T72h', 'SD14']:
        nd = sum(1 for r in mls_mobile if r['timepoint'] == tp and r['diet'] == 'ND')
        hfd = sum(1 for r in mls_mobile if r['timepoint'] == tp and r['diet'] == 'HFD')
        print(f"  {tp}: ND={nd}, HFD={hfd}")

# Write MLS mobile detail
if mls_mobile:
    outf3 = os.path.join(OUTDIR, 'mls_mobile_args.tsv')
    with open(outf3, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in mls_mobile:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    print(f"\n  Detail: {outf3}")

###########################################################################
# 13. Conjugation machinery detail
###########################################################################
print("\n" + "="*70)
print("=== CONJUGATION MACHINERY ON PLASMIDS ===")
print("="*70)

# Summarize conjugation gene types across all plasmid contigs
conj_type_counts = defaultdict(int)
conj_type_args = defaultdict(list)
conj_plasmid_info = []

for (s, c), pinfo in plasmid_contigs.items():
    cg = pinfo.get('conjugation_genes', 'NA')
    if not cg or cg == 'NA':
        continue
    conj_type_counts[cg] += 1
    # Find ARGs on this conjugative plasmid contig
    for r in results:
        if r['sample'] == s and r['contig'] == c:
            conj_type_args[cg].append(r)
            conj_plasmid_info.append({
                'sample': s,
                'contig': c,
                'conjugation_type': cg,
                'plasmid_score': pinfo['plasmid_score'],
                'plasmid_length': pinfo['length'],
                'topology': pinfo['topology'],
                'arg': r['best_hit_ARG'],
                'drug_class': r['drug_class'],
                'diet': r['diet'],
                'timepoint': r['timepoint'],
                'group': r['group'],
            })

print("\nConjugation gene types on plasmid contigs:")
for cg, count in sorted(conj_type_counts.items(), key=lambda x: -x[1]):
    n_args = len(conj_type_args[cg])
    print(f"  {cg}: {count} plasmid contigs, {n_args} ARGs")

if conj_plasmid_info:
    # Write conjugative plasmid ARG detail
    outf_conj = os.path.join(OUTDIR, 'conjugative_plasmid_args.tsv')
    with open(outf_conj, 'w') as f:
        f.write('sample\tgroup\tdiet\ttimepoint\tcontig\tconjugation_type\tplasmid_score\tplasmid_length\ttopology\targ\tdrug_class\n')
        for ci in conj_plasmid_info:
            f.write(f"{ci['sample']}\t{ci['group']}\t{ci['diet']}\t{ci['timepoint']}\t{ci['contig']}\t{ci['conjugation_type']}\t{ci['plasmid_score']}\t{ci['plasmid_length']}\t{ci['topology']}\t{ci['arg']}\t{ci['drug_class']}\n")
    print(f"\n  Conjugative plasmid ARG detail: {outf_conj}")

    print("\nConjugative plasmid ARGs by diet:")
    for diet in ['ND', 'HFD']:
        d_args = [ci for ci in conj_plasmid_info if ci['diet'] == diet]
        print(f"  {diet}: {len(d_args)} ARGs")

    print("\nConjugative plasmid ARGs by diet × timepoint:")
    for tp in ['T0h', 'T72h', 'SD14']:
        nd = sum(1 for ci in conj_plasmid_info if ci['timepoint'] == tp and ci['diet'] == 'ND')
        hfd = sum(1 for ci in conj_plasmid_info if ci['timepoint'] == tp and ci['diet'] == 'HFD')
        print(f"  {tp}: ND={nd}, HFD={hfd}")

    print("\nTop ARGs on conjugative plasmids:")
    cpa_counts = defaultdict(int)
    for ci in conj_plasmid_info:
        cpa_counts[ci['arg']] += 1
    for arg, count in sorted(cpa_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {arg}: {count}")

###########################################################################
# 14. Per-sample mobile ARG counts (for statistical testing)
###########################################################################
print("\n" + "="*70)
print("=== PER-SAMPLE MOBILE ARG COUNTS ===")
print("="*70)

outf4 = os.path.join(OUTDIR, 'per_sample_mobile_counts.tsv')
with open(outf4, 'w') as f:
    f.write('sample\tgroup\tdiet\tantibiotic\ttimepoint\ttotal_ARGs\tmobile_ARGs\tplasmid_ARGs\tIS_ARGs\tphage_ARGs\tchromosomal_ARGs\tmobile_fraction\n')
    sample_data = defaultdict(lambda: {'total': 0, 'mobile': 0, 'plasmid': 0, 'IS': 0, 'phage': 0, 'chromosomal': 0})
    for r in results:
        key = r['sample']
        sample_data[key]['total'] += 1
        sample_data[key]['group'] = r['group']
        sample_data[key]['diet'] = r['diet']
        sample_data[key]['abx'] = r['antibiotic']
        sample_data[key]['tp'] = r['timepoint']
        if r['mobilome'] == 'chromosomal':
            sample_data[key]['chromosomal'] += 1
        else:
            sample_data[key]['mobile'] += 1
        if 'plasmid' in r['mobilome']:
            sample_data[key]['plasmid'] += 1
        if any(x in r['mobilome'] for x in ['transposase','IS_element','integrase','recombinase']):
            sample_data[key]['IS'] += 1
        if 'phage' in r['mobilome']:
            sample_data[key]['phage'] += 1

    for s in sorted(sample_data):
        d = sample_data[s]
        frac = d['mobile'] / d['total'] if d['total'] else 0
        f.write(f"{s}\t{d['group']}\t{d['diet']}\t{d['abx']}\t{d['tp']}\t{d['total']}\t{d['mobile']}\t{d['plasmid']}\t{d['IS']}\t{d['phage']}\t{d['chromosomal']}\t{frac:.4f}\n")

print(f"  Per-sample counts: {outf4}")

print("\n" + "="*70)
print("=== ANALYSIS COMPLETE ===")
print("="*70)
print(f"\nOutput files in {OUTDIR}:")
for fn in sorted(os.listdir(OUTDIR)):
    if fn.endswith('.tsv'):
        fp = os.path.join(OUTDIR, fn)
        size = os.path.getsize(fp)
        print(f"  {fn} ({size/1024:.0f} KB)")

PYEOF

echo "$(timestamp) === Mobilome analysis finished ==="
