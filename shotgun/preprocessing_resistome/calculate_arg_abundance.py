#!/usr/bin/env python3
"""
Calculate ARG abundance normalized to 16S rRNA gene copies.
Logic: ARG copies per 16S copy = (ARG_reads / ARG_gene_length) / (16S_reads / 16S_gene_length)

Inputs:
- CARD idxstats files (*_card.txt) for ARG read counts per gene
- 16S read count files (*_16s_count.txt)
- CARD gene length file
- CARD metadata (aro_categories.tsv) for antibiotic class annotation

Outputs:
- arg_abundance_per_gene.tsv: per-gene ARG copies/16S copy for all samples
- arg_abundance_per_class.tsv: aggregated by antibiotic class
"""

import os
import json
import pandas as pd
import numpy as np

RESDIR = "/home/unnot/vault/kjm/hiseq/results"
CARD_DIR = "/home/unnot/CDC/refdb/card"
OUTDIR = "/home/unnot/vault/kjm/hiseq/arg_results"
os.makedirs(OUTDIR, exist_ok=True)

AVG_16S_LENGTH = 1500  # average 16S rRNA gene length in bp

# Load CARD gene lengths
gene_lengths = {}
with open("/home/unnot/vault/kjm/hiseq/card_bt2/card_gene_lengths.tsv") as f:
    for line in f:
        name, length = line.strip().split("\t")
        gene_lengths[name] = int(length)

# Load CARD ARO categories for antibiotic class mapping
# Parse gene name to ARO mapping from the FASTA headers
# idxstats uses full header up to first space as gene ID
# e.g., gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1
gene_to_info = {}
with open("/home/unnot/vault/kjm/hiseq/card_bt2/card_all.fasta") as f:
    for line in f:
        if line.startswith(">"):
            header = line.strip()[1:]
            gene_id = header.split()[0]  # everything before first space = idxstats ID
            parts = gene_id.split("|")
            aro = None
            gene_name = parts[-1] if parts else gene_id
            for p in parts:
                if p.startswith("ARO:"):
                    aro = p
                    break
            gene_to_info[gene_id] = {"aro": aro, "name": gene_name, "full_header": header}

# Load ARO index for drug class mapping (aro_index.tsv has Drug Class, Resistance Mechanism)
aro_categories = {}
aro_idx_file = os.path.join(CARD_DIR, "aro_index.tsv")
if os.path.exists(aro_idx_file):
    df_idx = pd.read_csv(aro_idx_file, sep="\t", low_memory=False)
    for _, row in df_idx.iterrows():
        aro = str(row.get("ARO Accession", ""))
        if aro.startswith("ARO:") and pd.notna(row.get("Drug Class")):
            aro_categories[aro] = {
                "drug_class": str(row["Drug Class"]),
                "resistance_mechanism": str(row.get("Resistance Mechanism", "Unknown")),
                "gene_family": str(row.get("AMR Gene Family", "Unknown")),
            }
    print(f"Loaded {len(aro_categories)} ARO-to-drug-class mappings")

# Collect data from all samples
samples = []
arg_counts = {}  # {sample: {gene: read_count}}
s16_counts = {}  # {sample: 16S_read_count}

for sdir in sorted(os.listdir(RESDIR)):
    if not sdir.startswith("ARG_"):
        continue
    spath = os.path.join(RESDIR, sdir)

    # CARD read counts
    card_file = os.path.join(spath, f"{sdir}_card.txt")
    s16_file = os.path.join(spath, f"{sdir}_16s_count.txt")

    if not os.path.exists(card_file) or not os.path.exists(s16_file):
        continue

    # Read CARD idxstats
    counts = {}
    with open(card_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                gene = parts[0]
                mapped = int(parts[2])
                if mapped > 0:
                    counts[gene] = mapped
    arg_counts[sdir] = counts

    # Read 16S count
    with open(s16_file) as f:
        s16 = int(f.read().strip())
    s16_counts[sdir] = s16

    samples.append(sdir)

print(f"Samples with both CARD and 16S data: {len(samples)}")

# Get all ARG genes detected
all_genes = set()
for counts in arg_counts.values():
    all_genes.update(counts.keys())
print(f"Total ARG genes detected: {len(all_genes)}")

# Calculate ARG copies per 16S copy
# Formula: (ARG_reads / ARG_length_kb) / (16S_reads / 16S_length_kb)
results = {}
for sample in samples:
    s16 = s16_counts[sample]
    if s16 == 0:
        print(f"WARNING: {sample} has 0 16S reads, skipping")
        continue
    s16_norm = s16 / (AVG_16S_LENGTH / 1000)  # 16S reads per kb

    sample_results = {}
    for gene in all_genes:
        reads = arg_counts[sample].get(gene, 0)
        if reads == 0:
            sample_results[gene] = 0
            continue
        gene_len = gene_lengths.get(gene, 1000)  # default 1kb if unknown
        arg_norm = reads / (gene_len / 1000)  # ARG reads per kb
        sample_results[gene] = arg_norm / s16_norm  # copies per 16S copy

    results[sample] = sample_results

# Build per-gene abundance table
df_gene = pd.DataFrame(results).fillna(0)
df_gene.index.name = "gene"

# Add gene annotations
annotations = []
for gene in df_gene.index:
    info = gene_to_info.get(gene, {})
    aro = info.get("aro", "")
    cat = aro_categories.get(aro, {})
    annotations.append({
        "gene": gene,
        "gene_name": info.get("name", gene),
        "aro": aro,
        "drug_class": cat.get("drug_class", "Unknown"),
        "resistance_mechanism": cat.get("resistance_mechanism", "Unknown"),
        "gene_family": cat.get("gene_family", "Unknown"),
    })
df_annot = pd.DataFrame(annotations).set_index("gene")

# Save per-gene table
df_out = pd.concat([df_annot, df_gene], axis=1)
df_out.to_csv(os.path.join(OUTDIR, "arg_abundance_per_gene.tsv"), sep="\t")
print(f"Per-gene table: {df_out.shape}")

# Aggregate by drug class
drug_class_results = {}
for sample in results:
    class_sums = {}
    for gene, abundance in results[sample].items():
        info = gene_to_info.get(gene, {})
        aro = info.get("aro", "")
        cat = aro_categories.get(aro, {})
        drug_class = cat.get("drug_class", "Unknown")
        # Drug class can have multiple classes separated by "; "
        for dc in drug_class.split("; "):
            dc = dc.strip()
            class_sums[dc] = class_sums.get(dc, 0) + abundance
    drug_class_results[sample] = class_sums

df_class = pd.DataFrame(drug_class_results).fillna(0).T
df_class.index.name = "sample"
df_class.to_csv(os.path.join(OUTDIR, "arg_abundance_per_class.tsv"), sep="\t")
print(f"Per-class table: {df_class.shape}")

# Summary
print(f"\n=== Top 15 drug classes by mean abundance ===")
means = df_class.mean().sort_values(ascending=False)
for dc, val in means.head(15).items():
    print(f"  {dc}: {val:.6f}")

# Save 16S counts for reference
with open(os.path.join(OUTDIR, "16s_read_counts.tsv"), "w") as f:
    f.write("sample\t16s_read_pairs\n")
    for s in sorted(s16_counts):
        f.write(f"{s}\t{s16_counts[s]}\n")

print(f"\nOutput saved to {OUTDIR}")
