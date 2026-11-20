#!/bin/bash
# GTDB-Tk classification of medium-quality MAGs (>=50% completeness, <10% contamination)
# Collects qualifying bins from all 72 samples, runs GTDB-Tk classify_wf,
# then merges taxonomy with ARG-bin linkage data.

set -euo pipefail

RESULTS=~/vault/kjm/hiseq/results
OUTDIR=~/vault/kjm/hiseq/gtdbtk_output
MAGDIR=$OUTDIR/all_mags
GTDBTK_DB=~/DB/gtdbtk
THREADS=64
LOG=~/vault/kjm/hiseq/gtdbtk.log

export GTDBTK_DATA_PATH=$GTDBTK_DB

echo "[$(date)] Starting GTDB-Tk pipeline" | tee $LOG

# Step 1: Collect medium-quality MAGs with unique names
mkdir -p $MAGDIR

echo "[$(date)] Collecting medium-quality MAGs..." | tee -a $LOG
COUNT=0
for SAMPLE_DIR in $RESULTS/ARG_*; do
    SAMPLE=$(basename $SAMPLE_DIR)
    QR=$SAMPLE_DIR/checkm2_output/quality_report.tsv
    BINDIR=$SAMPLE_DIR/dastool_output/das_DASTool_bins

    if [[ ! -f "$QR" ]] || [[ ! -d "$BINDIR" ]]; then
        echo "  Skipping $SAMPLE (missing checkm2 or bins)" | tee -a $LOG
        continue
    fi

    # Get medium+ quality bins (completeness>=50, contamination<10)
    while IFS=$'\t' read -r binname comp contam rest; do
        BINFILE="$BINDIR/${binname}.fa"
        if [[ -f "$BINFILE" ]]; then
            # Prefix with sample name to avoid collisions
            cp "$BINFILE" "$MAGDIR/${SAMPLE}__${binname}.fa"
            COUNT=$((COUNT + 1))
        fi
    done < <(tail -n+2 "$QR" | awk -F'\t' '$2>=50 && $3<10')
done

echo "[$(date)] Collected $COUNT medium-quality MAGs" | tee -a $LOG

# Step 2: Run GTDB-Tk classify_wf
echo "[$(date)] Running GTDB-Tk classify_wf..." | tee -a $LOG
conda run -n gtdbtk gtdbtk classify_wf \
    --genome_dir $MAGDIR \
    --out_dir $OUTDIR/classify \
    --extension fa \
    --cpus $THREADS \
    --pplacer_cpus 16 \
    --skip_ani_screen \
    2>&1 | tee -a $LOG

echo "[$(date)] GTDB-Tk classification complete" | tee -a $LOG

# Step 3: Merge taxonomy with ARG-bin linkage
echo "[$(date)] Merging taxonomy with ARG data..." | tee -a $LOG

python3 << 'PYEOF'
import pandas as pd
import os, glob

outdir = os.path.expanduser("~/vault/kjm/hiseq/gtdbtk_output")
results = os.path.expanduser("~/vault/kjm/hiseq/results")

# Load GTDB-Tk results (bacteria and archaea)
tax_dfs = []
for f in ["classify/gtdbtk.bac120.summary.tsv", "classify/gtdbtk.ar53.summary.tsv"]:
    fp = os.path.join(outdir, f)
    if os.path.exists(fp):
        df = pd.read_csv(fp, sep='\t')
        tax_dfs.append(df)

if not tax_dfs:
    print("ERROR: No GTDB-Tk output found")
    exit(1)

tax = pd.concat(tax_dfs, ignore_index=True)
print(f"GTDB-Tk classified {len(tax)} MAGs")

# Parse sample and bin from user_genome (format: ARG_XX_timepoint__binname)
tax['sample'] = tax['user_genome'].apply(lambda x: '__'.join(x.split('__')[:-1]))
tax['bin'] = tax['user_genome'].apply(lambda x: x.split('__')[-1])

# Load all ARG-bin linkage files
arg_dfs = []
for f in glob.glob(os.path.join(results, "ARG_*/*_arg_bins.tsv")):
    df = pd.read_csv(f, sep='\t')
    arg_dfs.append(df)

args = pd.concat(arg_dfs, ignore_index=True)
print(f"Total ARG-bin linkages: {len(args)}")

# Merge: add taxonomy to ARGs that are binned
# In arg_bins.tsv, 'bin' column has values like 'bin.2' or 'unbinned'
# In tax, 'bin' has values like 'bin.002_sub'
# Need to match sample+bin

# Create merge key
args['merge_sample'] = args['sample']
args['merge_bin'] = args['bin']

tax_lookup = tax[['sample', 'bin', 'classification']].copy()
tax_lookup.columns = ['merge_sample', 'merge_bin', 'classification']

merged = args.merge(tax_lookup, on=['merge_sample', 'merge_bin'], how='left')
merged.drop(columns=['merge_sample', 'merge_bin'], inplace=True)

# Save full merged table
out_file = os.path.join(outdir, "arg_mag_taxonomy.tsv")
merged.to_csv(out_file, sep='\t', index=False)
print(f"Saved merged ARG-MAG-taxonomy table: {out_file}")

# Summary: ARGs by taxonomy
binned = merged[merged['bin'] != 'unbinned'].copy()
binned_with_tax = binned.dropna(subset=['classification'])

print(f"\nBinned ARGs: {len(binned)}")
print(f"Binned ARGs with taxonomy: {len(binned_with_tax)}")
print(f"Unbinned ARGs: {len(merged[merged['bin'] == 'unbinned'])}")

# Parse genus from classification
binned_with_tax['genus'] = binned_with_tax['classification'].apply(
    lambda x: [t for t in x.split(';') if t.startswith('g__')][-1] if 'g__' in x else 'unclassified'
)

# Summary by drug class and genus
summary = binned_with_tax.groupby(['drug_class', 'genus']).size().reset_index(name='count')
summary = summary.sort_values(['drug_class', 'count'], ascending=[True, False])
summary_file = os.path.join(outdir, "arg_taxonomy_summary.tsv")
summary.to_csv(summary_file, sep='\t', index=False)
print(f"Saved summary: {summary_file}")

# Key question: Are MLS ARGs linked to Lactococcus MAGs?
mls_classes = ['macrolide antibiotic', 'lincosamide antibiotic', 'streptogramin antibiotic',
               'macrolide antibiotic;lincosamide antibiotic;streptogramin antibiotic']
mls = binned_with_tax[binned_with_tax['drug_class'].str.contains('macrolide|lincosamide|streptogramin', case=False, na=False)]
print(f"\n=== MLS ARGs in MAGs ===")
print(f"Total MLS ARGs in MAGs: {len(mls)}")
if len(mls) > 0:
    print(mls.groupby('genus')['best_hit_ARG'].count().sort_values(ascending=False).head(20))

# Lactococcus specifically
lacto = binned_with_tax[binned_with_tax['classification'].str.contains('Lactococcus', case=False, na=False)]
print(f"\n=== All ARGs in Lactococcus MAGs ===")
print(f"Total ARGs in Lactococcus MAGs: {len(lacto)}")
if len(lacto) > 0:
    print(lacto.groupby('drug_class')['best_hit_ARG'].count().sort_values(ascending=False).head(20))
PYEOF

echo "[$(date)] All done" | tee -a $LOG
