#!/bin/bash
# HiSeq Resistome Pipeline for KJM study
# Usage: bash run_resistome.sh <sample_name>
# e.g., bash run_resistome.sh ARG_AC6_T0h

set -e

SAMPLE=$1
INDIR="/home/unnot/vault/KJM_info/HiSeq"
OUTDIR="/home/unnot/vault/kjm/hiseq/results"
THREADS=8

# Mouse genome index
MOUSE_IDX="/home/unnot/vault/kjm/hiseq/mouse_genome/mm39"

# CARD database
CARD_JSON="/home/unnot/CDC/refdb/card/card.json"

R1="${INDIR}/${SAMPLE}_1.fastq.gz"
R2="${INDIR}/${SAMPLE}_2.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: Missing files for $SAMPLE"
    exit 1
fi

mkdir -p ${OUTDIR}/${SAMPLE}
cd ${OUTDIR}/${SAMPLE}

echo "=== Processing $SAMPLE ==="

# Step 1: QC with fastp
if [ ! -f "${SAMPLE}_qc_1.fastq.gz" ]; then
    echo "[$(date)] Step 1: Quality trimming with fastp..."
    fastp \
        -i $R1 -I $R2 \
        -o ${SAMPLE}_qc_1.fastq.gz -O ${SAMPLE}_qc_2.fastq.gz \
        -h ${SAMPLE}_fastp.html -j ${SAMPLE}_fastp.json \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread $THREADS \
        --detect_adapter_for_pe \
        2>&1 | tail -5
else
    echo "Step 1: fastp already done, skipping."
fi

# Step 2: Host removal with bowtie2
if [ ! -f "${SAMPLE}_nhost_1.fastq.gz" ]; then
    echo "[$(date)] Step 2: Removing host reads..."
    bowtie2 -x $MOUSE_IDX \
        -1 ${SAMPLE}_qc_1.fastq.gz -2 ${SAMPLE}_qc_2.fastq.gz \
        --un-conc-gz ${SAMPLE}_nhost_%.fastq.gz \
        -S /dev/null \
        --threads $THREADS \
        --very-sensitive \
        2>&1 | tail -3
    # Rename bowtie2 output (it uses _1 and _2 suffix)
    mv ${SAMPLE}_nhost_1.fastq.gz ${SAMPLE}_nhost_1.fastq.gz 2>/dev/null || true
    mv ${SAMPLE}_nhost_2.fastq.gz ${SAMPLE}_nhost_2.fastq.gz 2>/dev/null || true
else
    echo "Step 2: Host removal already done, skipping."
fi

# Step 3: RGI for ARG detection
if [ ! -f "${SAMPLE}_rgi.txt" ]; then
    echo "[$(date)] Step 3: Running RGI bwt..."
    rgi bwt \
        --read_one ${SAMPLE}_nhost_1.fastq.gz \
        --read_two ${SAMPLE}_nhost_2.fastq.gz \
        --output_file ${SAMPLE}_rgi \
        --threads $THREADS \
        --local \
        --include_wildcard \
        --clean \
        2>&1 | tail -5
else
    echo "Step 3: RGI already done, skipping."
fi

echo "[$(date)] === Done: $SAMPLE ==="
