#!/bin/bash
# Re-analyze ARGs using standard tools:
#   1. ARGs-OAP (read-based) — using host-removed reads
#   2. AMRFinderPlus (contig-based) — using assembled contigs + Prodigal proteins/GFF
#
# Run: nohup bash ~/vault/kjm/hiseq/run_arg_reanalysis.sh &> ~/vault/kjm/hiseq/arg_reanalysis.log &

set -euo pipefail

HISEQ_DIR="$HOME/vault/kjm/hiseq"
RESULTS_DIR="$HISEQ_DIR/results"
OUTDIR="$HISEQ_DIR/arg_reanalysis"
THREADS=16

mkdir -p "$OUTDIR"

echo "=== ARG Re-analysis Pipeline ==="
echo "Started: $(date)"
echo ""

# ═══════════════════════════════════════════════════════
# PART 1: ARGs-OAP (read-based)
# ═══════════════════════════════════════════════════════
echo "[PART 1] ARGs-OAP read-based analysis"

ARGSOAP_IN="$OUTDIR/argsoap_input"
ARGSOAP_OUT="$OUTDIR/argsoap_output"
ARGSOAP_META="$OUTDIR/argsoap_metadata.txt"

mkdir -p "$ARGSOAP_IN"

# Prepare input: symlink host-removed reads with correct naming
# ARGs-OAP expects paired files named *_1.fq and *_2.fq in input dir
echo "  Preparing input reads..."
for sample_dir in "$RESULTS_DIR"/ARG_*/; do
    SAMPLE=$(basename "$sample_dir")
    R1="${sample_dir}/${SAMPLE}_nhost_1.fastq.gz"
    R2="${sample_dir}/${SAMPLE}_nhost_2.fastq.gz"
    if [ -f "$R1" ] && [ -f "$R2" ]; then
        ln -sf "$R1" "$ARGSOAP_IN/${SAMPLE}_1.fq.gz"
        ln -sf "$R2" "$ARGSOAP_IN/${SAMPLE}_2.fq.gz"
    fi
done

NSAMP=$(ls "$ARGSOAP_IN"/*_1.fq.gz 2>/dev/null | wc -l)
echo "  Samples with reads: $NSAMP"

# Create metadata file (SampleID\tCategory)
echo -e "SampleID\tCategory" > "$ARGSOAP_META"
for f in "$ARGSOAP_IN"/*_1.fq.gz; do
    SAMPLE=$(basename "$f" | sed 's/_1\.fq\.gz//')
    echo -e "${SAMPLE}\t${SAMPLE}" >> "$ARGSOAP_META"
done

# Run Stage 1
echo "  Running ARGs-OAP stage_one..."
ARGs_OAP stage_one \
    -i "$ARGSOAP_IN" \
    -o "$ARGSOAP_OUT" \
    -m "$ARGSOAP_META" \
    -n "$THREADS" \
    -f fq \
    -z 2>&1 | tee "$OUTDIR/argsoap_stage1.log"

echo "  Running ARGs-OAP stage_two..."
ARGs_OAP stage_two \
    -i "$ARGSOAP_OUT" \
    -o "$ARGSOAP_OUT/stage_two" \
    -m "$ARGSOAP_OUT/metadata.txt" \
    -n "$THREADS" 2>&1 | tee "$OUTDIR/argsoap_stage2.log"

echo "  ARGs-OAP complete: $(date)"
echo ""

# ═══════════════════════════════════════════════════════
# PART 2: AMRFinderPlus (contig-based)
# ═══════════════════════════════════════════════════════
echo "[PART 2] AMRFinderPlus contig-based analysis"

AMRFINDER_OUT="$OUTDIR/amrfinderplus"
mkdir -p "$AMRFINDER_OUT"

# Process each sample
DONE=0
TOTAL=$(ls -d "$RESULTS_DIR"/ARG_*/ | wc -l)

for sample_dir in "$RESULTS_DIR"/ARG_*/; do
    SAMPLE=$(basename "$sample_dir")
    CONTIGS="${sample_dir}/megahit/final.contigs.fa"
    PROTEINS="${sample_dir}/${SAMPLE}_proteins.faa"
    GFF="${sample_dir}/${SAMPLE}_genes.gff"
    OUTFILE="$AMRFINDER_OUT/${SAMPLE}_amrfinder.tsv"

    if [ -f "$OUTFILE" ]; then
        DONE=$((DONE + 1))
        echo "  [$DONE/$TOTAL] $SAMPLE — already done, skipping"
        continue
    fi

    if [ ! -f "$CONTIGS" ]; then
        echo "  [$DONE/$TOTAL] $SAMPLE — no contigs, skipping"
        DONE=$((DONE + 1))
        continue
    fi

    DONE=$((DONE + 1))
    echo "  [$DONE/$TOTAL] $SAMPLE — running AMRFinderPlus..."

    # Run AMRFinderPlus with both protein and nucleotide input
    # --plus enables stress/virulence/biocide detection
    # --annotation_format prodigal for Prodigal GFF
    if [ -f "$PROTEINS" ] && [ -f "$GFF" ]; then
        conda run -n amr amrfinder \
            --protein "$PROTEINS" \
            --nucleotide "$CONTIGS" \
            --gff "$GFF" \
            --annotation_format prodigal \
            --plus \
            --threads 4 \
            --output "$OUTFILE" 2>/dev/null
    else
        # Nucleotide-only mode
        conda run -n amr amrfinder \
            --nucleotide "$CONTIGS" \
            --plus \
            --threads 4 \
            --output "$OUTFILE" 2>/dev/null
    fi
done

echo "  AMRFinderPlus complete: $(date)"
echo ""

# ═══════════════════════════════════════════════════════
# PART 3: Summarize AMRFinderPlus results
# ═══════════════════════════════════════════════════════
echo "[PART 3] Summarizing AMRFinderPlus results"

# Merge all AMRFinderPlus outputs
MERGED="$AMRFINDER_OUT/all_amrfinder_results.tsv"
first=true
for f in "$AMRFINDER_OUT"/*_amrfinder.tsv; do
    SAMPLE=$(basename "$f" | sed 's/_amrfinder\.tsv//')
    if $first; then
        echo -e "sample\t$(head -1 "$f")" > "$MERGED"
        first=false
    fi
    tail -n+2 "$f" | while IFS= read -r line; do
        echo -e "${SAMPLE}\t${line}"
    done >> "$MERGED"
done

TOTAL_HITS=$(tail -n+2 "$MERGED" | wc -l)
UNIQUE_GENES=$(tail -n+2 "$MERGED" | cut -f8 | sort -u | wc -l)
echo "  Total AMRFinderPlus hits: $TOTAL_HITS"
echo "  Unique gene symbols: $UNIQUE_GENES"

# Summary by element type and subclass
echo "  Generating summary tables..."
python3 << 'PYEOF'
import pandas as pd
import os

outdir = os.path.expanduser("~/vault/kjm/hiseq/arg_reanalysis/amrfinderplus")
df = pd.read_csv(f"{outdir}/all_amrfinder_results.tsv", sep="\t")

# Summary by sample and element type
summary = df.groupby(["sample", "Element type"]).size().unstack(fill_value=0)
summary.to_csv(f"{outdir}/summary_by_element_type.tsv", sep="\t")

# Summary by sample and subclass
subclass = df.groupby(["sample", "Element subclass"]).size().unstack(fill_value=0)
subclass.to_csv(f"{outdir}/summary_by_subclass.tsv", sep="\t")

# Point mutations specifically
if "Method" in df.columns:
    mutations = df[df["Method"].str.contains("POINT", na=False)]
    if len(mutations) > 0:
        mutations.to_csv(f"{outdir}/point_mutations.tsv", sep="\t", index=False)
        print(f"  Point mutations found: {len(mutations)}")
    else:
        print("  No point mutations detected")

# AMR genes only (exclude stress/virulence)
amr_only = df[df["Element type"] == "AMR"]
amr_summary = amr_only.groupby(["sample", "Class"]).size().unstack(fill_value=0)
amr_summary.to_csv(f"{outdir}/amr_by_class.tsv", sep="\t")
print(f"  AMR-only hits: {len(amr_only)}")
print(f"  Unique AMR classes: {amr_only['Class'].nunique()}")
PYEOF

echo ""
echo "=== Pipeline Complete ==="
echo "Finished: $(date)"
echo ""
echo "Output files:"
echo "  ARGs-OAP:      $ARGSOAP_OUT/"
echo "  AMRFinderPlus:  $AMRFINDER_OUT/"
echo "  Merged results: $AMRFINDER_OUT/all_amrfinder_results.tsv"
echo "  Summaries:      $AMRFINDER_OUT/summary_by_*.tsv"
