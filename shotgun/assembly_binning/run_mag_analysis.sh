#!/bin/bash
# MAG analysis: GTDB-Tk taxonomy + ARG-MAG linkage + plasmid integration
# Goal: confirm/refute whether Lactococcus MAGs carry MLS ARGs
#
# Steps:
#   1. Collect all MQ+ bins (completeness>=50, contamination<10) across 72 samples
#   2. Run GTDB-Tk classify_wf on all MQ+ bins
#   3. Join taxonomy → bins → ARGs → plasmid data
#   4. Generate summary tables for manuscript
#
# Usage: nohup bash run_mag_analysis.sh &> mag_analysis.log &

set -euo pipefail

BASEDIR="$HOME/vault/kjm/hiseq"
RESULTS="$BASEDIR/results"
OUTDIR="$BASEDIR/mag_analysis"
BINDIR="$OUTDIR/all_mq_bins"
GTDBTK_ENV="gtdbtk"
export GTDBTK_DATA_PATH="$HOME/DB/gtdbtk"
CPUS=32
PPLACER_CPUS=4  # each pplacer thread uses ~10GB RAM

echo "=== MAG Analysis Pipeline ==="
echo "Started: $(date)"
echo ""

# ─── Step 1: Collect all MQ+ bins ───────────────────────────────────────────
echo "[Step 1] Collecting medium+ quality bins..."
mkdir -p "$BINDIR"

# Create master ARG-bin linkage table
ARG_BIN_TABLE="$OUTDIR/all_arg_bins.tsv"
echo -e "sample\tcontig\torf_id\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tbin\tbin_global_id" > "$ARG_BIN_TABLE"

# Create master CheckM2 quality table
QUALITY_TABLE="$OUTDIR/all_checkm2_quality.tsv"
echo -e "bin_global_id\tsample\toriginal_bin\tcompleteness\tcontamination\tgenome_size\tgc_content\tn50\ttotal_contigs\tquality" > "$QUALITY_TABLE"

# Create master plasmid table
PLASMID_TABLE="$OUTDIR/all_plasmid_summary.tsv"
echo -e "sample\tcontig\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tgenomic_location" > "$PLASMID_TABLE"

n_samples=0
n_bins_total=0
n_bins_mq=0

for sample_dir in "$RESULTS"/ARG_*; do
    SAMPLE=$(basename "$sample_dir")
    checkm="$sample_dir/checkm2_output/quality_report.tsv"
    arg_bins="$sample_dir/${SAMPLE}_arg_bins.tsv"
    das_bins="$sample_dir/dastool_output/das_DASTool_bins"
    plasmid="$sample_dir/${SAMPLE}_plasmid_summary.tsv"

    if [[ ! -f "$checkm" ]]; then
        echo "  WARNING: no CheckM2 for $SAMPLE, skipping"
        continue
    fi

    n_samples=$((n_samples + 1))

    # Process each bin in CheckM2 output
    while IFS=$'\t' read -r bin_name completeness contamination _ _ _ n50 _ genome_size gc _ total_contigs _ _; do
        n_bins_total=$((n_bins_total + 1))
        global_id="${SAMPLE}__${bin_name}"

        if (( $(echo "$completeness >= 90" | bc -l) )) && (( $(echo "$contamination < 5" | bc -l) )); then
            quality="HQ"
        elif (( $(echo "$completeness >= 50" | bc -l) )) && (( $(echo "$contamination < 10" | bc -l) )); then
            quality="MQ"
        else
            quality="LQ"
        fi

        echo -e "${global_id}\t${SAMPLE}\t${bin_name}\t${completeness}\t${contamination}\t${genome_size}\t${gc}\t${n50}\t${total_contigs}\t${quality}" >> "$QUALITY_TABLE"

        if [[ "$quality" == "HQ" || "$quality" == "MQ" ]]; then
            n_bins_mq=$((n_bins_mq + 1))
            bin_fa="$das_bins/${bin_name}.fa"
            if [[ -f "$bin_fa" ]]; then
                cp "$bin_fa" "$BINDIR/${global_id}.fa"
            else
                echo "  WARNING: bin file not found: $bin_fa"
            fi
        fi
    done < <(tail -n +2 "$checkm")

    # Append ARG-bin linkages with global bin IDs
    if [[ -f "$arg_bins" ]]; then
        tail -n +2 "$arg_bins" | while IFS=$'\t' read -r s contig orf arg drug mech bin; do
            if [[ "$bin" == "unbinned" ]]; then
                global_id="unbinned"
            else
                global_id="${SAMPLE}__${bin}"
            fi
            echo -e "${s}\t${contig}\t${orf}\t${arg}\t${drug}\t${mech}\t${bin}\t${global_id}"
        done >> "$ARG_BIN_TABLE"
    fi

    # Append plasmid summary
    if [[ -f "$plasmid" ]]; then
        tail -n +2 "$plasmid" >> "$PLASMID_TABLE"
    fi
done

echo "  Samples processed: $n_samples"
echo "  Total bins: $n_bins_total"
echo "  MQ+ bins copied: $n_bins_mq"
echo "  Bin files in $BINDIR: $(ls "$BINDIR"/*.fa 2>/dev/null | wc -l)"
echo ""

# ─── Step 2: Run GTDB-Tk ────────────────────────────────────────────────────
echo "[Step 2] Running GTDB-Tk classify_wf on $n_bins_mq MQ+ bins..."
echo "  This is the slow step (pplacer). Expect 1-3 hours."
GTDBTK_OUT="$OUTDIR/gtdbtk_output"
mkdir -p "$GTDBTK_OUT"

conda run -n "$GTDBTK_ENV" \
    gtdbtk classify_wf \
    --genome_dir "$BINDIR" \
    --out_dir "$GTDBTK_OUT" \
    --extension fa \
    --cpus "$CPUS" \
    --pplacer_cpus "$PPLACER_CPUS" \
    --prefix mag \
    --skip_ani_screen

echo "  GTDB-Tk done: $(date)"
echo ""

# ─── Step 3: Parse GTDB-Tk results and join with ARG data ───────────────────
echo "[Step 3] Joining taxonomy with ARG-bin linkages..."

TAXONOMY_TABLE="$OUTDIR/all_taxonomy.tsv"
echo -e "bin_global_id\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tclassification" > "$TAXONOMY_TABLE"

for summary in "$GTDBTK_OUT"/mag.*.summary.tsv; do
    [[ -f "$summary" ]] || continue
    tail -n +2 "$summary" | while IFS=$'\t' read -r bin_id classification rest; do
        domain=$(echo "$classification" | sed 's/;.*//' | sed 's/d__//')
        phylum=$(echo "$classification" | grep -o 'p__[^;]*' | sed 's/p__//')
        class=$(echo "$classification" | grep -o 'c__[^;]*' | sed 's/c__//')
        order=$(echo "$classification" | grep -o 'o__[^;]*' | sed 's/o__//')
        family=$(echo "$classification" | grep -o 'f__[^;]*' | sed 's/f__//')
        genus=$(echo "$classification" | grep -o 'g__[^;]*' | sed 's/g__//')
        species=$(echo "$classification" | grep -o 's__[^;]*' | sed 's/s__//')
        echo -e "${bin_id}\t${domain}\t${phylum}\t${class}\t${order}\t${family}\t${genus}\t${species}\t${classification}"
    done >> "$TAXONOMY_TABLE"
done

echo "  Classified MAGs: $(( $(wc -l < "$TAXONOMY_TABLE") - 1 ))"
echo ""

# ─── Step 4: Generate summary tables ────────────────────────────────────────
echo "[Step 4] Generating summary tables..."

python3 << 'PYEOF'
import pandas as pd
import os, re, sys

outdir = os.path.expanduser("~/vault/kjm/hiseq/mag_analysis")

# Load tables
arg_bins = pd.read_csv(f"{outdir}/all_arg_bins.tsv", sep="\t")
taxonomy = pd.read_csv(f"{outdir}/all_taxonomy.tsv", sep="\t")
quality  = pd.read_csv(f"{outdir}/all_checkm2_quality.tsv", sep="\t")
plasmid  = pd.read_csv(f"{outdir}/all_plasmid_summary.tsv", sep="\t")

print(f"  ARG-bin rows: {len(arg_bins)}")
print(f"  Taxonomy rows: {len(taxonomy)}")
print(f"  Quality rows: {len(quality)}")
print(f"  Plasmid rows: {len(plasmid)}")

# ── Parse sample metadata ──
prefix_map = {
    'AC': ('AMX', 'CTL'), 'AH': ('AMX', 'HFD'),
    'CC': ('CIP', 'CTL'), 'CH': ('CIP', 'HFD'),
    'EC': ('ERY', 'CTL'), 'EH': ('ERY', 'HFD'),
    'CTL': ('NONE', 'CTL'), 'HFD': ('NONE', 'HFD'),
    'CSL': ('NONE', 'CTL'),  # CSL = mislabelled CTL at recovery timepoints (same mice)
}

def parse_sample(sample_id):
    parts = sample_id.replace("ARG_", "").split("_")
    timepoint = parts[-1]
    prefix_part = "_".join(parts[:-1])
    m = re.match(r'([A-Za-z]+)(\d+)', prefix_part)
    if m:
        prefix, mouse = m.group(1), m.group(2)
    else:
        prefix, mouse = prefix_part, "?"
    abx, diet = prefix_map.get(prefix, ('UNK', 'UNK'))
    group = f"{abx}_{diet}" if abx != "NONE" else diet
    return abx, diet, group, timepoint, mouse

meta = arg_bins['sample'].apply(parse_sample)
arg_bins[['antibiotic', 'diet', 'group', 'timepoint', 'mouse_id']] = pd.DataFrame(
    meta.tolist(), index=arg_bins.index)

# ── Merge ARGs with taxonomy + quality ──
arg_tax = arg_bins.merge(taxonomy, on="bin_global_id", how="left")
arg_tax = arg_tax.merge(
    quality[['bin_global_id', 'completeness', 'contamination', 'quality']],
    on="bin_global_id", how="left")

# ── Merge with plasmid data (contig-level) ──
plasmid_loc = plasmid[['sample', 'contig', 'genomic_location']].drop_duplicates()
arg_tax = arg_tax.merge(plasmid_loc, on=['sample', 'contig'], how='left')

arg_tax.to_csv(f"{outdir}/table_arg_mag_taxonomy.tsv", sep="\t", index=False)

# ── Binning rate summary ──
total_args = len(arg_tax)
binned_args = len(arg_tax[arg_tax['bin_global_id'] != 'unbinned'])
mq_binned = arg_tax[(arg_tax['bin_global_id'] != 'unbinned') &
                     (arg_tax['quality'].isin(['HQ', 'MQ']))]
print(f"\n  === Binning Summary ===")
print(f"  Total assembly-based ARGs: {total_args}")
print(f"  Binned (any quality): {binned_args} ({binned_args/total_args*100:.1f}%)")
print(f"  In MQ+ bins: {len(mq_binned)} ({len(mq_binned)/total_args*100:.1f}%)")

# ── Plasmid summary ──
print(f"\n  === Plasmid Summary ===")
for loc in ['chromosome', 'plasmid']:
    n = len(arg_tax[arg_tax['genomic_location'] == loc])
    print(f"  {loc}: {n}")
plasmid_args = arg_tax[arg_tax['genomic_location'] == 'plasmid']
if len(plasmid_args) > 0:
    print(f"  Plasmid ARG drug classes:")
    for dc, cnt in plasmid_args['drug_class'].value_counts().head(10).items():
        print(f"    {dc}: {cnt}")

# ── Lactococcus MAGs ──
lacto = mq_binned[mq_binned['genus'].str.contains('Lactococcus', na=False, case=False)]
print(f"\n  === Lactococcus Analysis ===")
print(f"  ARGs in Lactococcus MQ+ MAGs: {len(lacto)}")
if len(lacto) > 0:
    print(f"  Unique ARGs: {lacto['best_hit_ARG'].nunique()}")
    print(f"  Drug classes:")
    for dc, cnt in lacto['drug_class'].value_counts().items():
        print(f"    {dc}: {cnt}")
    print(f"  By timepoint and diet:")
    for tp in ['T0h', 'T72h', 'SD14']:
        for diet in ['CTL', 'HFD']:
            n = len(lacto[(lacto['timepoint'] == tp) & (lacto['diet'] == diet)])
            if n > 0:
                print(f"    {tp} {diet}: {n} ARGs")
    lacto.to_csv(f"{outdir}/table_lactococcus_args.tsv", sep="\t", index=False)

# ── Top genera by ARG count (MQ+) ──
genus_summary = mq_binned.groupby('genus').agg(
    total_ARGs=('best_hit_ARG', 'count'),
    unique_ARGs=('best_hit_ARG', 'nunique'),
    n_MAGs=('bin_global_id', 'nunique'),
    drug_classes=('drug_class', lambda x: '; '.join(sorted(set(x.dropna())))),
).sort_values('total_ARGs', ascending=False)
genus_summary.to_csv(f"{outdir}/table_genus_arg_summary.tsv", sep="\t")
print(f"\n  === Top 15 Genera by ARG Count (MQ+ bins) ===")
print(genus_summary[['total_ARGs', 'unique_ARGs', 'n_MAGs']].head(15).to_string())

# ── MLS ARGs by genus ──
mls_keywords = ['macrolide', 'lincosamide', 'streptogramin']
mls_mask = mq_binned['drug_class'].str.lower().apply(
    lambda x: any(k in str(x) for k in mls_keywords) if pd.notna(x) else False
)
mls_args = mq_binned[mls_mask]
print(f"\n  === MLS ARGs in MQ+ MAGs ===")
print(f"  Total MLS ARGs: {len(mls_args)}")
if len(mls_args) > 0:
    mls_genus = mls_args.groupby('genus').agg(
        n_MLS_ARGs=('best_hit_ARG', 'count'),
        unique_MLS_ARGs=('best_hit_ARG', 'nunique'),
        n_MAGs=('bin_global_id', 'nunique'),
        ARG_names=('best_hit_ARG', lambda x: '; '.join(sorted(set(x)))),
    ).sort_values('n_MLS_ARGs', ascending=False)
    mls_genus.to_csv(f"{outdir}/table_mls_arg_by_genus.tsv", sep="\t")
    print(mls_genus.to_string())

# ── HFD vs ND comparison by timepoint ──
print(f"\n  === ARGs in MQ+ MAGs by Diet and Timepoint ===")
diet_tp = mq_binned.groupby(['timepoint', 'diet']).agg(
    total_ARGs=('best_hit_ARG', 'count'),
    n_MAGs=('bin_global_id', 'nunique'),
    n_genera=('genus', 'nunique'),
).reset_index()
for tp in ['T0h', 'T72h', 'SD14']:
    print(f"  {tp}:")
    sub = diet_tp[diet_tp['timepoint'] == tp]
    for _, row in sub.iterrows():
        print(f"    {row['diet']}: {row['total_ARGs']} ARGs in {row['n_MAGs']} MAGs ({row['n_genera']} genera)")

# ── Plasmid ARGs on MQ+ MAGs ──
plasmid_mq = mq_binned[mq_binned['genomic_location'] == 'plasmid']
print(f"\n  === Plasmid-borne ARGs in MQ+ MAGs ===")
print(f"  Total: {len(plasmid_mq)}")
if len(plasmid_mq) > 0:
    print(f"  By genus:")
    for genus, cnt in plasmid_mq['genus'].value_counts().head(10).items():
        print(f"    {genus}: {cnt}")

# ── Per-sample summary ──
sample_summary = arg_tax.groupby('sample').agg(
    total_ARGs=('best_hit_ARG', 'count'),
    binned_ARGs=('bin_global_id', lambda x: (x != 'unbinned').sum()),
    mq_binned_ARGs=('quality', lambda x: x.isin(['HQ', 'MQ']).sum()),
    n_plasmid_ARGs=('genomic_location', lambda x: (x == 'plasmid').sum()),
    n_genera=('genus', lambda x: x.dropna().nunique()),
).reset_index()
meta2 = sample_summary['sample'].apply(parse_sample)
sample_summary[['antibiotic', 'diet', 'group', 'timepoint', 'mouse_id']] = pd.DataFrame(
    meta2.tolist(), index=sample_summary.index)
sample_summary.to_csv(f"{outdir}/table_per_sample_summary.tsv", sep="\t", index=False)

# ── MAG quality summary ──
print(f"\n  === MAG Quality Summary ===")
q_counts = quality['quality'].value_counts()
for q in ['HQ', 'MQ', 'LQ']:
    print(f"  {q}: {q_counts.get(q, 0)}")

print(f"\n  All tables saved to: {outdir}")
print("  Done!")
PYEOF

echo ""
echo "=== Pipeline complete ==="
echo "Finished: $(date)"
echo ""
echo "Key output files:"
echo "  $OUTDIR/all_arg_bins.tsv             - All ARG-bin linkages (72 samples)"
echo "  $OUTDIR/all_taxonomy.tsv             - GTDB-Tk taxonomy for all MQ+ MAGs"
echo "  $OUTDIR/all_checkm2_quality.tsv      - CheckM2 quality for all bins"
echo "  $OUTDIR/all_plasmid_summary.tsv      - Plasmid/chromosome classification"
echo "  $OUTDIR/table_arg_mag_taxonomy.tsv   - Full joined table (ARG+MAG+taxonomy+plasmid)"
echo "  $OUTDIR/table_lactococcus_args.tsv   - ARGs in Lactococcus MAGs"
echo "  $OUTDIR/table_genus_arg_summary.tsv  - ARG counts by genus"
echo "  $OUTDIR/table_mls_arg_by_genus.tsv   - MLS ARGs by genus"
echo "  $OUTDIR/table_per_sample_summary.tsv - Per-sample summary"
