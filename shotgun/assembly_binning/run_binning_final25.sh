#!/bin/bash
# Bin the 25 remaining samples that were assembled after the previous binning run.

RESDIR="/home/unnot/vault/kjm/hiseq/results"
THREADS=16

echo "=== Binning final 25 samples ==="
for SDIR in ${RESDIR}/ARG_*; do
  SAMPLE=$(basename "$SDIR")
  [ ! -f "${SDIR}/megahit/final.contigs.fa" ] && continue
  [ -f "${SDIR}/bins_done" ] && continue

  (
    cd "$SDIR"
    CONTIGS="megahit/final.contigs.fa"

    # Step 1: Map reads
    if [ ! -f "${SAMPLE}_sorted.bam" ]; then
      echo "[$(date +%H:%M)] Mapping: $SAMPLE"
      bowtie2-build --threads $THREADS $CONTIGS ${SAMPLE}_bt2idx 2>/dev/null &&
      bowtie2 -x ${SAMPLE}_bt2idx \
        -1 ${SAMPLE}_nhost_1.fastq.gz \
        -2 ${SAMPLE}_nhost_2.fastq.gz \
        --threads $THREADS --no-unal 2>/dev/null | \
      samtools sort -@ 4 -o ${SAMPLE}_sorted.bam - &&
      samtools index ${SAMPLE}_sorted.bam
      rm -f ${SAMPLE}_bt2idx.*.bt2
    fi

    [ ! -f "${SAMPLE}_sorted.bam" ] && echo "MAPPING FAILED: $SAMPLE" && exit 1

    # Step 2: MetaBAT2
    if [ ! -d "metabat2_bins" ]; then
      echo "[$(date +%H:%M)] MetaBAT2: $SAMPLE"
      conda run -n metabat2 bash -c "
        jgi_summarize_bam_contig_depths --outputDepth ${SAMPLE}_depth.txt ${SAMPLE}_sorted.bam 2>/dev/null &&
        metabat2 -i $CONTIGS -a ${SAMPLE}_depth.txt -o metabat2_bins/bin \
          -t $THREADS --minContig 1500 2>/dev/null
      "
      [ ! -d "metabat2_bins" ] && mkdir -p metabat2_bins
    fi

    # Step 3: MaxBin2
    if [ ! -d "maxbin2_bins" ]; then
      echo "[$(date +%H:%M)] MaxBin2: $SAMPLE"
      mkdir -p maxbin2_bins
      conda run -n maxbin2 bash -c "
        run_MaxBin.pl -contig $CONTIGS -out maxbin2_bins/bin \
          -abund ${SAMPLE}_depth.txt \
          -thread $THREADS 2>/dev/null
      "
    fi

    # Step 4: DAS Tool
    if [ ! -d "dastool_output" ]; then
      echo "[$(date +%H:%M)] DAS Tool: $SAMPLE"

      for f in metabat2_bins/bin.*.fa; do
        [ ! -f "$f" ] && continue
        bin=$(basename "$f" .fa)
        grep '^>' "$f" | sed "s/^>//" | awk -v b="$bin" '{print $1"\t"b}'
      done > metabat2_contigs2bin.tsv

      for f in maxbin2_bins/bin.*.fasta; do
        [ ! -f "$f" ] && continue
        bin=$(basename "$f" .fasta)
        grep '^>' "$f" | sed "s/^>//" | awk -v b="$bin" '{print $1"\t"b}'
      done > maxbin2_contigs2bin.tsv

      if [ -s "metabat2_contigs2bin.tsv" ] || [ -s "maxbin2_contigs2bin.tsv" ]; then
        BINS_INPUT=""
        LABELS=""
        if [ -s "metabat2_contigs2bin.tsv" ]; then
          BINS_INPUT="metabat2_contigs2bin.tsv"
          LABELS="metabat2"
        fi
        if [ -s "maxbin2_contigs2bin.tsv" ]; then
          [ -n "$BINS_INPUT" ] && BINS_INPUT="${BINS_INPUT},maxbin2_contigs2bin.tsv" && LABELS="${LABELS},maxbin2"
          [ -z "$BINS_INPUT" ] && BINS_INPUT="maxbin2_contigs2bin.tsv" && LABELS="maxbin2"
        fi

        conda run -n dastool bash -c "
          DAS_Tool -i $BINS_INPUT -l $LABELS \
            -c $CONTIGS -o dastool_output/das \
            --threads $THREADS --write_bins \
            --score_threshold 0.3 2>/dev/null
        "
      else
        echo "  No bins produced for $SAMPLE"
        mkdir -p dastool_output
      fi
    fi

    # Step 5: CheckM2 (with -x fa fix)
    if [ ! -d "checkm2_output" ] && [ -d "dastool_output/das_DASTool_bins" ]; then
      n_bins=$(ls dastool_output/das_DASTool_bins/*.fa 2>/dev/null | wc -l)
      if [ "$n_bins" -gt 0 ]; then
        echo "[$(date +%H:%M)] CheckM2: $SAMPLE ($n_bins bins)"
        conda run -n checkm2 bash -c "
          checkm2 predict --input dastool_output/das_DASTool_bins/ \
            --output-directory checkm2_output \
            -x fa --threads $THREADS 2>/dev/null
        "
      else
        echo "  No DAS Tool bins for $SAMPLE"
        mkdir -p checkm2_output
      fi
    fi

    # Step 6: Link ARGs to bins
    if [ -f "${SAMPLE}_rgi_main.txt" ] && [ -d "dastool_output/das_DASTool_bins" ]; then
      echo "[$(date +%H:%M)] Linking ARGs to bins: $SAMPLE"
      export SAMPLE
      python3 - <<'PYEOF'
import os, glob, collections

sample = os.environ["SAMPLE"]

contig2bin = {}
for fa in glob.glob("dastool_output/das_DASTool_bins/*.fa"):
    bin_name = os.path.basename(fa).replace(".fa", "")
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                contig = line.strip()[1:].split()[0]
                contig2bin[contig] = bin_name

with open(f"{sample}_rgi_main.txt") as f:
    header = f.readline()
    with open(f"{sample}_arg_bins.tsv", "w") as out:
        out.write("sample\tcontig\torf_id\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tbin\n")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            orf_id = parts[0].split(" # ")[0]
            contig = orf_id.rsplit("_", 1)[0]
            best_hit = parts[8] if len(parts) > 8 else "NA"
            drug_class = parts[14] if len(parts) > 14 else "NA"
            mechanism = parts[15] if len(parts) > 15 else "NA"
            bin_name = contig2bin.get(contig, "unbinned")
            out.write(f"{sample}\t{contig}\t{orf_id}\t{best_hit}\t{drug_class}\t{mechanism}\t{bin_name}\n")

bins = collections.Counter()
with open(f"{sample}_arg_bins.tsv") as f:
    next(f)
    for line in f:
        b = line.strip().split("\t")[-1]
        bins[b] += 1

n_binned = sum(v for k, v in bins.items() if k != "unbinned")
n_unbinned = bins.get("unbinned", 0)
n_bins_with_args = len([k for k in bins if k != "unbinned"])
print(f"  {sample}: {n_binned + n_unbinned} ARGs, {n_binned} binned ({n_bins_with_args} bins), {n_unbinned} unbinned")
PYEOF
    fi

    touch bins_done
    echo "[$(date +%H:%M)] Binning complete: $SAMPLE"
  ) &
  [ $(jobs -r | wc -l) -ge 2 ] && wait -n
done
wait

echo "=== All done ==="
echo "Total binned samples: $(ls ${RESDIR}/ARG_*/bins_done 2>/dev/null | wc -l)/72"
echo "Samples with CheckM2: $(ls ${RESDIR}/ARG_*/checkm2_output/quality_report.tsv 2>/dev/null | wc -l)"
