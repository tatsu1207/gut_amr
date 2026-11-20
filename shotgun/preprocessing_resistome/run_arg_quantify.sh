#!/bin/bash
# ARG-OAP-like quantification pipeline
# Step 1: Map reads to 16S (SILVA) for normalization
# Step 2: Use existing CARD mapping results
# Step 3: Calculate ARG copies per 16S copy
# Run 8 samples in parallel

RESDIR="/home/unnot/vault/kjm/hiseq/results"
SILVA_IDX="/home/unnot/vault/kjm/hiseq/silva_16s/silva_16s"
THREADS=8

# Step 1: Map to 16S
for SDIR in ${RESDIR}/ARG_*; do
  SAMPLE=$(basename "$SDIR")
  [ ! -f "${SDIR}/${SAMPLE}_nhost_1.fastq.gz" ] && continue
  [ -f "${SDIR}/${SAMPLE}_16s_count.txt" ] && continue

  (
    echo "[$(date +%H:%M)] 16S mapping: $SAMPLE"
    bowtie2 -x $SILVA_IDX \
      -1 ${SDIR}/${SAMPLE}_nhost_1.fastq.gz \
      -2 ${SDIR}/${SAMPLE}_nhost_2.fastq.gz \
      -S /dev/null \
      --threads $THREADS \
      --no-unal \
      --very-sensitive 2>&1 | tee ${SDIR}/${SAMPLE}_16s_bt2.log | grep "overall alignment rate"

    # Extract 16S read count from bowtie2 log
    ALIGNED=$(grep "aligned concordantly exactly 1 time\|aligned concordantly >1 times" ${SDIR}/${SAMPLE}_16s_bt2.log | awk '{sum+=$1}END{print sum}')
    echo "$ALIGNED" > ${SDIR}/${SAMPLE}_16s_count.txt
    echo "  $SAMPLE: $ALIGNED 16S read pairs"
  ) &
  [ $(jobs -r | wc -l) -ge 8 ] && wait -n
done
wait
echo "=== All 16S mapping done ==="
