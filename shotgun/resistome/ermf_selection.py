#!/usr/bin/env python3
"""ErmF in the erythromycin groups (Results, Fig. 4C; 2026-11-21): per-mouse read-level abundance, share of total
ARG abundance, assembly context, and 16S Bacteroides in the same mice."""
import pandas as pd
H = "/home/unnot/vault/kjm/hiseq/"; D = "/home/unnot/vault/kjm/miseq/dada2_output/"
pg = pd.read_csv(H + "arg_results/model_filtered/arg_abundance_per_gene_homolog_only.tsv", sep="\t")
S = [c for c in pg.columns if c.startswith("ARG_")]
e = pg[pg.gene_name == "ErmF"][S].sum(); tot = pg[S].sum()
asv = pd.read_csv(D + "asv_table.tsv", sep="\t", index_col=0).drop(columns="sequence")
tax = pd.read_csv(D + "taxonomy.tsv", sep="\t", index_col=0)
g = asv.groupby(tax.loc[asv.index, "Genus"]).sum(); pct = g.div(g.sum()) * 100
for m in ["EH3", "EH4", "EH7", "EC3", "EC6", "EC8"]:
    v = [e[f"ARG_{m}_{t}"] for t in ("T0h", "T72h", "SD14")]
    b = [pct.at["Bacteroides", f"ARG_{m}_{t}"] for t in ("T0h", "T72h", "SD14")]
    fold = f"{v[1] / v[0]:.0f}x" if v[0] > 0 else "-"
    print(f"{m} ErmF {v[0]:.5f} {v[1]:.4f} {v[2]:.5f} ({fold}; {v[1] / tot[f'ARG_{m}_T72h'] * 100:.1f}% of ARG at T3) | "
          f"Bacteroides % {b[0]:.1f} {b[1]:.1f} {b[2]:.1f} ({b[1] / b[0]:.1f}x)")
other = e[[c for c in S if not c.startswith("ARG_EH")]]
print("max ErmF in any other sample:", round(other.max(), 5), other.idxmax())
f = pd.read_csv(H + "mag_analysis/table_arg_mag_taxonomy.tsv", sep="\t")
print(f[f.best_hit_ARG == "ErmF"][["sample", "genus", "species", "genomic_location"]].to_string(index=False))
