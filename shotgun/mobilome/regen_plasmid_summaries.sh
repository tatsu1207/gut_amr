#!/bin/bash
# Regenerate all plasmid summaries with fixed contig-matching logic
RESDIR="/home/unnot/vault/kjm/hiseq/results"

for SDIR in ${RESDIR}/ARG_*; do
  SAMPLE=$(basename "$SDIR")
  # Only redo samples that have RGI + geNomad done
  [ ! -f "${SDIR}/${SAMPLE}_rgi_main.txt" ] && continue
  [ ! -d "${SDIR}/genomad_output" ] && continue

  # Remove old summary
  rm -f "${SDIR}/${SAMPLE}_plasmid_summary.tsv"

  cd "$SDIR"
  export SAMPLE SDIR
  python3 - <<'PYEOF'
import os

sample = os.environ["SAMPLE"]

# Load plasmid contigs from geNomad
plasmid_contigs = set()
plasmid_file = "genomad_output/final.contigs_summary/final.contigs_plasmid.fna"
if os.path.exists(plasmid_file):
    with open(plasmid_file) as f:
        for line in f:
            if line.startswith(">"):
                plasmid_contigs.add(line.strip()[1:].split()[0])

# Load virus contigs
virus_contigs = set()
virus_file = "genomad_output/final.contigs_summary/final.contigs_virus.fna"
if os.path.exists(virus_file):
    with open(virus_file) as f:
        for line in f:
            if line.startswith(">"):
                virus_contigs.add(line.strip()[1:].split()[0])

# Load RGI results
rgi_file = f"{sample}_rgi_main.txt"
arg_results = []
if os.path.exists(rgi_file):
    with open(rgi_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            orf_id = parts[0].split(" # ")[0]
            contig = orf_id.rsplit("_", 1)[0]
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

with open(f"{sample}_plasmid_summary.tsv", "w") as f:
    f.write("sample\tcontig\tbest_hit_ARG\tdrug_class\tresistance_mechanism\tgenomic_location\n")
    for r in arg_results:
        f.write(f"{sample}\t{r['contig']}\t{r['best_hit']}\t{r['drug_class']}\t{r['mechanism']}\t{r['location']}\n")

n_total = len(arg_results)
n_plasmid = sum(1 for r in arg_results if r["location"] == "plasmid")
n_phage = sum(1 for r in arg_results if r["location"] == "phage")
n_chr = n_total - n_plasmid - n_phage
print(f"  {sample}: {n_total} ARGs ({n_plasmid} plasmid, {n_phage} phage, {n_chr} chromosome)")
PYEOF
done
echo "Done."
