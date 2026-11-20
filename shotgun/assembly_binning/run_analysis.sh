#!/bin/bash
# MAG analysis: resume GTDB-Tk (identify done, align+classify needed) + ARG-MAG linkage
# Safe to run in background with nohup
set -euo pipefail

RESDIR=~/vault/kjm/hiseq/results
OUTDIR=~/vault/kjm/hiseq/mag_analysis
GTDB_DATA=/home/unnot/CDC/refdb/gtdb/data
THREADS=8

# Previous GTDB-Tk run (identify already completed for 1535 MAGs)
PREV_GTDBTK=~/vault/kjm/hiseq/gtdbtk_output/classify

mkdir -p "$OUTDIR"
exec > >(tee -a "$OUTDIR/analysis.log") 2>&1

timestamp() { date "+[%Y-%m-%d %H:%M:%S]"; }
echo "$(timestamp) === Starting MAG analysis (resuming) ==="

###############################################################################
# STEP 1: Collect quality MAGs + ARG-bin data (quick)
###############################################################################
echo "$(timestamp) Step 1: Collecting quality MAGs from CheckM2 results"

QUALITY_MAGS="$OUTDIR/quality_mags.tsv"
echo -e "sample\tbin\tcompleteness\tcontamination\tgenome_size\tN50\tGC\tcontigs" > "$QUALITY_MAGS"

for SDIR in "$RESDIR"/ARG_*; do
    SAMPLE=$(basename "$SDIR")
    QR="$SDIR/checkm2_output/quality_report.tsv"
    [[ -f "$QR" ]] || continue
    awk -F'\t' -v s="$SAMPLE" 'NR>1 && $2>=50 && $3<10 {
        print s"\t"$1"\t"$2"\t"$3"\t"$9"\t"$7"\t"$10"\t"$12
    }' "$QR" >> "$QUALITY_MAGS"
done
echo "$(timestamp)   Found $(tail -n+2 "$QUALITY_MAGS" | wc -l) quality MAGs"

echo "$(timestamp) Step 2: Consolidating ARG-bin assignments"
ARG_BINS="$OUTDIR/all_arg_bins.tsv"
echo -e "sample\tcontig\torf_id\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tbin" > "$ARG_BINS"
for SDIR in "$RESDIR"/ARG_*; do
    SAMPLE=$(basename "$SDIR")
    ABF="$SDIR/${SAMPLE}_arg_bins.tsv"
    [[ -f "$ABF" ]] && tail -n+2 "$ABF" >> "$ARG_BINS"
done
TOTAL_ARGS=$(tail -n+2 "$ARG_BINS" | wc -l)
BINNED_ARGS=$(tail -n+2 "$ARG_BINS" | awk -F'\t' '$7!="unbinned"' | wc -l)
echo "$(timestamp)   Total ARGs: $TOTAL_ARGS, Binned: $BINNED_ARGS ($(echo "scale=1; $BINNED_ARGS*100/$TOTAL_ARGS" | bc)%)"

###############################################################################
# STEP 3: Resume GTDB-Tk — reuse identify, run classify_wf with correct DB
###############################################################################
echo "$(timestamp) Step 3: Running GTDB-Tk classify_wf (identify already done, correct DB path)"

export GTDBTK_DATA_PATH="$GTDB_DATA"

# Need genome_dir with the bins — check if previous run had one
# The previous run used a genome_dir; let's find it
PREV_GENOME_DIR=""
if grep -q "genome_dir" "$PREV_GTDBTK/gtdbtk.json" 2>/dev/null; then
    PREV_GENOME_DIR=$(python3 -c "import json; d=json.load(open('$PREV_GTDBTK/gtdbtk.json')); print(d.get('genome_dir',''))" 2>/dev/null)
fi

# If we can't find it, prepare bins ourselves
GTDBTK_INDIR="$OUTDIR/gtdbtk_bins"
if [[ -n "$PREV_GENOME_DIR" && -d "$PREV_GENOME_DIR" ]]; then
    GTDBTK_INDIR="$PREV_GENOME_DIR"
    echo "$(timestamp)   Reusing genome dir: $GTDBTK_INDIR"
else
    echo "$(timestamp)   Preparing bin files..."
    mkdir -p "$GTDBTK_INDIR"
    while IFS=$'\t' read -r sample bin comp cont gsize n50 gc ncontigs; do
        [[ "$sample" == "sample" ]] && continue
        BINFA="$RESDIR/$sample/dastool_output/das_DASTool_bins/${bin}.fa"
        DEST="$GTDBTK_INDIR/${sample}__${bin}.fa"
        [[ -f "$BINFA" && ! -f "$DEST" ]] && cp "$BINFA" "$DEST"
    done < "$QUALITY_MAGS"
    echo "$(timestamp)   Prepared $(ls "$GTDBTK_INDIR"/*.fa | wc -l) bin files"
fi

GTDBTK_OUT="$OUTDIR/gtdbtk_output"
mkdir -p "$GTDBTK_OUT"

conda run -n gtdbtk \
    gtdbtk classify_wf \
    --genome_dir "$GTDBTK_INDIR" \
    --out_dir "$GTDBTK_OUT" \
    --extension fa \
    --cpus "$THREADS" \
    --skip_ani_screen \
    --pplacer_cpus 4

echo "$(timestamp)   GTDB-Tk classification complete"

###############################################################################
# STEP 4: Merge taxonomy with ARG-bin data
###############################################################################
echo "$(timestamp) Step 4: Merging taxonomy with ARG-bin linkage"

python3 - "$ARG_BINS" "$GTDBTK_OUT" "$QUALITY_MAGS" "$OUTDIR" <<'PYEOF'
import sys, csv, os, glob
from collections import defaultdict

arg_bins_file = sys.argv[1]
gtdbtk_dir = sys.argv[2]
quality_file = sys.argv[3]
outdir = sys.argv[4]

# Load GTDB-Tk taxonomy from summary files
tax = {}  # (sample, bin) -> classification
for summary_f in glob.glob(os.path.join(gtdbtk_dir, "classify", "gtdbtk.*.summary.tsv")):
    with open(summary_f) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mag_id = row['user_genome']
            # Format: ARG_XX#_timepoint__bin.###
            idx = mag_id.index('__')
            sample = mag_id[:idx]
            binname = mag_id[idx+2:]
            tax[(sample, binname)] = row['classification']

print(f"Loaded taxonomy for {len(tax)} MAGs")

# Load quality info
quality = {}
with open(quality_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        quality[(row['sample'], row['bin'])] = row

# Load ARG-bin assignments
args = []
with open(arg_bins_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        args.append(row)

# Write full ARG-MAG-taxonomy table
outf = os.path.join(outdir, 'arg_mag_taxonomy.tsv')
with open(outf, 'w') as f:
    f.write('sample\tcontig\torf_id\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tbin\ttaxonomy\tcompleteness\tcontamination\n')
    for a in args:
        s, b = a['sample'], a['bin']
        t = tax.get((s, b), 'unclassified')
        q = quality.get((s, b), {})
        comp = q.get('completeness', 'NA')
        cont = q.get('contamination', 'NA')
        if b == 'unbinned':
            t = 'unbinned'
            comp = cont = 'NA'
        f.write(f"{s}\t{a['contig']}\t{a['orf_id']}\t{a['best_hit_ARG']}\t{a['drug_class']}\t{a['resistance_mechanism']}\t{b}\t{t}\t{comp}\t{cont}\n")
print(f"Full table: {outf}")

# Parse group/diet from sample name
def parse_sample(s):
    parts = s.replace('ARG_', '').split('_')
    prefix = parts[0][:2]
    tp = parts[1]
    group_map = {'AC':'AMX_CTL','AH':'AMX_HFD','CC':'CIP_CTL','CH':'CIP_HFD',
                 'EC':'ERY_CTL','EH':'ERY_HFD','CT':'CTL','HF':'HFD'}
    group = group_map.get(prefix, prefix)
    diet = 'HFD' if 'HFD' in group or group == 'HFD' else 'ND'
    return group, diet, tp

# Genus-drug summary
genus_drug = defaultdict(int)
genus_diet_drug = defaultdict(lambda: defaultdict(int))
for a in args:
    if a['bin'] == 'unbinned':
        continue
    s, b = a['sample'], a['bin']
    t = tax.get((s, b), 'unclassified')
    genus = 'unknown'
    for level in t.split(';'):
        if level.startswith('g__') and len(level) > 3:
            genus = level[3:]
    group, diet, tp = parse_sample(s)
    drug = a['drug_class']
    genus_drug[(genus, drug)] += 1
    genus_diet_drug[(genus, drug)][diet] += 1

outf2 = os.path.join(outdir, 'genus_drug_summary.tsv')
with open(outf2, 'w') as f:
    f.write('genus\tdrug_class\ttotal_ARGs\tND_count\tHFD_count\n')
    for (genus, drug), total in sorted(genus_drug.items(), key=lambda x: -x[1]):
        nd = genus_diet_drug[(genus, drug)].get('ND', 0)
        hfd = genus_diet_drug[(genus, drug)].get('HFD', 0)
        f.write(f"{genus}\t{drug}\t{total}\t{nd}\t{hfd}\n")
print(f"Genus-drug summary: {outf2}")

# === Lactococcus-specific analysis ===
print("\n" + "="*60)
print("=== Lactococcus ARG analysis ===")
print("="*60)
lacto_args = []
for a in args:
    if a['bin'] == 'unbinned':
        continue
    s, b = a['sample'], a['bin']
    t = tax.get((s, b), 'unclassified')
    if 'Lactococcus' in t:
        group, diet, tp = parse_sample(s)
        lacto_args.append({**a, 'taxonomy': t, 'group': group, 'diet': diet, 'timepoint': tp})

if lacto_args:
    print(f"Lactococcus-linked ARGs: {len(lacto_args)}")
    drug_counts = defaultdict(int)
    for a in lacto_args:
        drug_counts[a['drug_class']] += 1
    print("Drug classes:")
    for drug, count in sorted(drug_counts.items(), key=lambda x: -x[1]):
        print(f"  {drug}: {count}")
    diet_counts = defaultdict(int)
    for a in lacto_args:
        diet_counts[a['diet']] += 1
    print(f"By diet: ND={diet_counts.get('ND',0)}, HFD={diet_counts.get('HFD',0)}")

    # By timepoint
    tp_diet = defaultdict(lambda: defaultdict(int))
    for a in lacto_args:
        tp_diet[a['timepoint']][a['diet']] += 1
    print("By timepoint × diet:")
    for tp in ['T0h', 'T72h', 'SD14']:
        print(f"  {tp}: ND={tp_diet[tp].get('ND',0)}, HFD={tp_diet[tp].get('HFD',0)}")

    outf3 = os.path.join(outdir, 'lactococcus_args.tsv')
    with open(outf3, 'w') as f:
        f.write('sample\tgroup\tdiet\ttimepoint\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tbin\ttaxonomy\n')
        for a in lacto_args:
            f.write(f"{a['sample']}\t{a['group']}\t{a['diet']}\t{a['timepoint']}\t{a['best_hit_ARG']}\t{a['drug_class']}\t{a['resistance_mechanism']}\t{a['bin']}\t{a['taxonomy']}\n")
    print(f"Detail: {outf3}")
else:
    print("No Lactococcus-linked ARGs found in binned data")

# === Overall stats ===
print(f"\n{'='*60}")
print("=== Overall MAG stats ===")
print(f"{'='*60}")
total_binned = sum(1 for a in args if a['bin'] != 'unbinned')
total_unbinned = sum(1 for a in args if a['bin'] == 'unbinned')
print(f"Total ARGs: {len(args)}")
print(f"Binned: {total_binned} ({total_binned*100//len(args)}%)")
print(f"Unbinned: {total_unbinned}")

genus_total = defaultdict(int)
for a in args:
    if a['bin'] == 'unbinned':
        continue
    s, b = a['sample'], a['bin']
    t = tax.get((s, b), 'unclassified')
    genus = 'unknown'
    for level in t.split(';'):
        if level.startswith('g__') and len(level) > 3:
            genus = level[3:]
    genus_total[genus] += 1
print("\nTop 20 ARG-carrying genera:")
for genus, count in sorted(genus_total.items(), key=lambda x: -x[1])[:20]:
    print(f"  {genus}: {count}")

PYEOF

echo "$(timestamp) === All analysis complete ==="
echo ""
echo "Output files in $OUTDIR:"
ls -lh "$OUTDIR"/*.tsv 2>/dev/null
