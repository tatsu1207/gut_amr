#!/bin/bash
# Re-generate plasmid summaries for samples that have unknown_plasmid_summary.tsv
RESDIR="/home/unnot/vault/kjm/hiseq/results"

for bad_file in ${RESDIR}/*/unknown_plasmid_summary.tsv; do
  [ ! -f "$bad_file" ] && continue
  SDIR=$(dirname "$bad_file")
  SAMPLE=$(basename "$SDIR")

  echo "Fixing: $SAMPLE"
  rm "$bad_file"

  cd "$SDIR"
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
            contig = parts[1].rsplit("_", 1)[0] if "_" in parts[1] else parts[1]
            best_hit = parts[8] if len(parts) > 8 else "NA"
            drug_class = parts[14] if len(parts) > 14 else "NA"
            mechanism = parts[15] if len(parts) > 15 else "NA"

            location = "chromosome"
            if any(contig.startswith(pc.split("|")[0]) or pc.startswith(contig) for pc in plasmid_contigs):
                location = "plasmid"
            elif any(contig.startswith(vc.split("|")[0]) or vc.startswith(contig) for vc in virus_contigs):
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
done
echo "Done fixing plasmid summaries"
