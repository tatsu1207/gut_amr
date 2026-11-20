#!/bin/bash
# Full assembly pipeline: MEGAHIT → Prodigal → RGI → geNomad → MOB-typer
# Run 4 samples in parallel

RESDIR="/home/unnot/vault/kjm/hiseq/results"
GENOMAD_DB="/home/unnot/DB/genomad_db"
THREADS=16

for SDIR in ${RESDIR}/ARG_*; do
  SAMPLE=$(basename "$SDIR")
  [ ! -f "${SDIR}/${SAMPLE}_nhost_1.fastq.gz" ] && continue
  # Skip if already fully done
  [ -f "${SDIR}/${SAMPLE}_plasmid_summary.tsv" ] && continue

  (
    cd "$SDIR"

    # Step 1: MEGAHIT assembly
    if [ ! -d "megahit" ]; then
      echo "[$(date +%H:%M)] MEGAHIT: $SAMPLE"
      megahit \
        -1 ${SAMPLE}_nhost_1.fastq.gz \
        -2 ${SAMPLE}_nhost_2.fastq.gz \
        -o megahit \
        -t $THREADS \
        --min-contig-len 500 \
        --presets meta-sensitive 2>/dev/null
    fi

    [ ! -f "megahit/final.contigs.fa" ] && echo "MEGAHIT FAILED: $SAMPLE" && exit 1

    # Step 2: Prodigal gene prediction
    if [ ! -f "${SAMPLE}_proteins.faa" ]; then
      echo "[$(date +%H:%M)] Prodigal: $SAMPLE"
      prodigal \
        -i megahit/final.contigs.fa \
        -o ${SAMPLE}_genes.gff \
        -a ${SAMPLE}_proteins.faa \
        -d ${SAMPLE}_genes.fna \
        -p meta -f gff 2>/dev/null
    fi

    # Step 3: RGI main on predicted proteins
    if [ ! -f "${SAMPLE}_rgi_main.txt" ]; then
      echo "[$(date +%H:%M)] RGI main: $SAMPLE"
      cd /home/unnot/vault/kjm/hiseq/card_db
      rgi main \
        -i ${SDIR}/${SAMPLE}_proteins.faa \
        -o ${SDIR}/${SAMPLE}_rgi_main \
        -t protein \
        -a DIAMOND \
        -n $THREADS \
        --local \
        --clean 2>/dev/null
      cd "$SDIR"
    fi

    # Step 4: geNomad - identify plasmids and viruses
    if [ ! -d "genomad_output" ]; then
      echo "[$(date +%H:%M)] geNomad: $SAMPLE"
      genomad end-to-end \
        megahit/final.contigs.fa \
        genomad_output \
        $GENOMAD_DB \
        --threads $THREADS \
        --cleanup 2>/dev/null
    fi

    # Step 5: MOB-typer on plasmid contigs identified by geNomad
    if [ ! -f "${SAMPLE}_mob_typer.tsv" ]; then
      # Extract plasmid contigs from geNomad results
      PLASMID_FASTA="genomad_output/final.contigs_summary/final.contigs_plasmid.fna"
      if [ -f "$PLASMID_FASTA" ] && [ -s "$PLASMID_FASTA" ]; then
        echo "[$(date +%H:%M)] MOB-typer: $SAMPLE ($(grep -c '^>' $PLASMID_FASTA) plasmid contigs)"
        # MOB-typer works on individual sequences, run on each plasmid contig
        > ${SAMPLE}_mob_typer.tsv
        mob_typer \
          --infile "$PLASMID_FASTA" \
          --out_file ${SAMPLE}_mob_typer.tsv \
          --num_threads $THREADS 2>/dev/null
      else
        echo "  $SAMPLE: No plasmid contigs found by geNomad"
        echo "no_plasmids" > ${SAMPLE}_mob_typer.tsv
      fi
    fi

    # Step 6: Create summary - link ARGs to plasmids
    echo "[$(date +%H:%M)] Summary: $SAMPLE"
    export SAMPLE SDIR
    python3 - <<'PYEOF'
import os, sys

sample = os.environ.get("SAMPLE", "unknown")
sdir = os.environ.get("SDIR", ".")

# Load plasmid contigs from geNomad
plasmid_contigs = set()
plasmid_file = f"genomad_output/final.contigs_summary/final.contigs_plasmid.fna"
if os.path.exists(plasmid_file):
    with open(plasmid_file) as f:
        for line in f:
            if line.startswith(">"):
                plasmid_contigs.add(line.strip()[1:].split()[0].replace("|", "_"))

# Load virus contigs
virus_contigs = set()
virus_file = f"genomad_output/final.contigs_summary/final.contigs_virus.fna"
if os.path.exists(virus_file):
    with open(virus_file) as f:
        for line in f:
            if line.startswith(">"):
                virus_contigs.add(line.strip()[1:].split()[0].replace("|", "_"))

# Load RGI results
rgi_file = f"{sample}_rgi_main.txt"
arg_results = []
if os.path.exists(rgi_file):
    with open(rgi_file) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            # ORF_ID format: k141_225762_1 # 202 # 714 # ... (contig_orf # start # stop ...)
            # Contig column (parts[1]) is empty in protein mode, so parse from ORF_ID
            orf_id = parts[0].split(" # ")[0]  # "k141_225762_1"
            contig = orf_id.rsplit("_", 1)[0]  # "k141_225762"
            best_hit = parts[8] if len(parts) > 8 else "NA"
            drug_class = parts[14] if len(parts) > 14 else "NA"
            mechanism = parts[15] if len(parts) > 15 else "NA"

            location = "chromosome"
            if contig in plasmid_contigs:
                location = "plasmid"
            elif contig in virus_contigs:
                location = "phage"

            arg_results.append({
                "contig": contig,
                "best_hit": best_hit,
                "drug_class": drug_class,
                "mechanism": mechanism,
                "location": location,
            })

# Write summary
with open(f"{sample}_plasmid_summary.tsv", "w") as f:
    f.write("sample\tcontig\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tgenomic_location\n")
    for r in arg_results:
        f.write(f"{sample}\t{r['contig']}\t{r['best_hit']}\t{r['drug_class']}\t{r['mechanism']}\t{r['location']}\n")

n_total = len(arg_results)
n_plasmid = sum(1 for r in arg_results if r["location"] == "plasmid")
n_phage = sum(1 for r in arg_results if r["location"] == "phage")
print(f"  {sample}: {n_total} ARGs ({n_plasmid} plasmid, {n_phage} phage, {n_total-n_plasmid-n_phage} chromosome), {len(plasmid_contigs)} plasmid contigs, {len(virus_contigs)} virus contigs")
PYEOF

  ) &
  [ $(jobs -r | wc -l) -ge 4 ] && wait -n
done
wait
echo "=== All assembly + annotation done ==="
echo "Assemblies: $(ls -d ${RESDIR}/*/megahit 2>/dev/null | wc -l)/72"
echo "RGI main: $(ls ${RESDIR}/*/*_rgi_main.txt 2>/dev/null | wc -l)/72"
echo "geNomad: $(ls -d ${RESDIR}/*/genomad_output 2>/dev/null | wc -l)/72"
echo "Summaries: $(ls ${RESDIR}/*/*_plasmid_summary.tsv 2>/dev/null | wc -l)/72"
