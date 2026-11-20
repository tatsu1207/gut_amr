#!/bin/bash
# MEGAHIT assembly + Prodigal + RGI for all HiSeq samples
# Run 4 assemblies in parallel

RESDIR="/home/unnot/vault/kjm/hiseq/results"
THREADS=16

for SDIR in ${RESDIR}/ARG_*; do
  SAMPLE=$(basename "$SDIR")
  [ ! -f "${SDIR}/${SAMPLE}_nhost_1.fastq.gz" ] && continue
  [ -d "${SDIR}/megahit" ] && continue

  (
    echo "[$(date +%H:%M)] MEGAHIT: $SAMPLE"
    megahit \
      -1 ${SDIR}/${SAMPLE}_nhost_1.fastq.gz \
      -2 ${SDIR}/${SAMPLE}_nhost_2.fastq.gz \
      -o ${SDIR}/megahit \
      -t $THREADS \
      --min-contig-len 500 \
      --presets meta-sensitive 2>/dev/null

    # Prodigal gene prediction
    echo "[$(date +%H:%M)] Prodigal: $SAMPLE"
    prodigal \
      -i ${SDIR}/megahit/final.contigs.fa \
      -o ${SDIR}/${SAMPLE}_genes.gff \
      -a ${SDIR}/${SAMPLE}_proteins.faa \
      -d ${SDIR}/${SAMPLE}_genes.fna \
      -p meta -f gff 2>/dev/null

    # RGI on contigs (main mode, not bwt)
    echo "[$(date +%H:%M)] RGI main: $SAMPLE"
    rgi main \
      -i ${SDIR}/${SAMPLE}_proteins.faa \
      -o ${SDIR}/${SAMPLE}_rgi_main \
      -t protein \
      -a DIAMOND \
      -n $THREADS \
      --local \
      --clean 2>/dev/null

    # Count results
    if [ -f "${SDIR}/${SAMPLE}_rgi_main.txt" ]; then
      ARGS=$(tail -n+2 ${SDIR}/${SAMPLE}_rgi_main.txt | wc -l)
      CONTIGS=$(grep -c "^>" ${SDIR}/megahit/final.contigs.fa)
      echo "  $SAMPLE: $CONTIGS contigs, $ARGS ARGs"
    fi
  ) &
  [ $(jobs -r | wc -l) -ge 4 ] && wait -n
done
wait
echo "=== All assembly + RGI done ==="
