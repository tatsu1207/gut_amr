#!/usr/bin/env python3
"""Re-implementation of the Methods' ARG mobility rule (the original script is not in the
project): mobile if the contig is a geNomad plasmid or virus, or any gene on it has a
geNomad annotation matching a mobile-element keyword or a ConjScan hit. Validated against
mobilome_analysis/arg_mobilome_full.tsv before use (2026-11-20)."""
import re, sys, functools
import pandas as pd

R = "/home/unnot/vault/kjm/hiseq/results"
KW = re.compile(r"transposase|\bIS\d|insertion sequence|integrase|recombinase|invertase|resolvase|"
                r"relaxase|conjugat|\btra[A-Z]\b|\btrb[A-Z]\b|\bmob[A-Z]\b|type IV secretion", re.I)

@functools.lru_cache(maxsize=None)
def load(sample):
    g = f"{R}/{sample}/genomad_output"
    genes = pd.read_csv(f"{g}/final.contigs_annotate/final.contigs_genes.tsv", sep="\t",
                        usecols=["gene", "annotation_conjscan", "annotation_description"])
    genes["contig"] = genes.gene.str.rsplit("_", n=1).str[0]
    hit = genes.annotation_description.fillna("").str.contains(KW) | genes.annotation_conjscan.notna()
    marker = set(genes.loc[hit, "contig"])
    pl = set(pd.read_csv(f"{g}/final.contigs_summary/final.contigs_plasmid_summary.tsv", sep="\t").seq_name)
    vi = set(pd.read_csv(f"{g}/final.contigs_summary/final.contigs_virus_summary.tsv", sep="\t")
             .seq_name.str.split("|").str[0])
    return marker, pl, vi

def mobile(sample, contig):
    marker, pl, vi = load(sample)
    return contig in marker or contig in pl or contig in vi

if __name__ == "__main__":
    full = pd.read_csv("/home/unnot/vault/kjm/hiseq/mobilome_analysis/arg_mobilome_full.tsv", sep="\t")
    test = full[full["sample"].isin(sorted(full["sample"].unique())[:: int(sys.argv[1]) if len(sys.argv) > 1 else 6])]
    pred = [mobile(s, c) for s, c in zip(test["sample"], test.contig)]
    truth = (test.mobilome != "chromosomal").tolist()
    agree = sum(p == t for p, t in zip(pred, truth))
    print(f"validation: {agree}/{len(truth)} agree; truth mobile {sum(truth)}, predicted mobile {sum(pred)}")
