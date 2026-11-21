#!/usr/bin/env python3
"""but (butyryl-CoA:acetate CoA-transferase) copies per genome split by carrier genus, per sample (2026-11-21).
Each sample's but abundance (butyrate_gene_abundance_per_sample.tsv) is divided among carrier genera in proportion
to the summed contig depth of its but hits (butyrate_gene_hits_with_taxonomy.tsv; genera from GTDB-Tk of the MAG,
"Unbinned" = contig not in a MAG). Writes tables/but_carriers_per_sample.tsv (Supplementary Data 12f)."""
import pandas as pd
T = "/home/unnot/vault/kjm/hiseq/butyrate_analysis/tables/"
hits = pd.read_csv(T + "butyrate_gene_hits_with_taxonomy.tsv", sep="\t")
hits = hits[hits.gene == "but"]
ab = pd.read_csv(T + "butyrate_gene_abundance_per_sample.tsv", sep="\t", keep_default_na=False).set_index("sample")
rows = []
for s, x in hits.groupby("sample"):
    sh = x.groupby("genus").depth.sum(); sh = sh / sh.sum()
    meta = x.iloc[0]
    for g, v in sh.sort_values(ascending=False).items():
        rows.append(dict(sample=s, group=meta.group, mouse=meta.mouse, timepoint=meta.timepoint, carrier_genus=g,
                         but_copies_per_genome=round(v * ab.at[s, "but"], 4), share_of_sample_but=round(v, 3),
                         sample_total_but=round(ab.at[s, "but"], 4)))
d = pd.DataFrame(rows)
d.to_csv(T + "but_carriers_per_sample.tsv", sep="\t", index=False)
for m in ["AMX_CTL_6", "AMX_CTL_7", "AMX_CTL_8", "AMX_HFD_3", "AMX_HFD_4", "AMX_HFD_5"]:
    for tp in ("T0h", "SD14"):
        x = d[(d.mouse == m) & (d.timepoint == tp)]
        print(m, tp, x.sample_total_but.iloc[0] if len(x) else 0,
              "; ".join(f"{g} {v} ({p:.0%})" for g, v, p in zip(x.carrier_genus, x.but_copies_per_genome, x.share_of_sample_but)))
