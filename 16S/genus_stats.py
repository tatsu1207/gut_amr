#!/usr/bin/env python3
"""Genus-level numbers quoted in the Results (group mean +/- SD of % of genus-assigned reads).

Samples: experiment == "ARG", the eight study groups, >= 5,000 reads (the same set as every other 16S
analysis; 451 samples). Writes dada2_output/diversity/genus_group_stats.tsv (all groups x timepoints x
top genera) and prints the values the text quotes. Run with --nofilter to reproduce the pre-fix values
(which wrongly kept the 33 low-read libraries, five of them empty).
"""
import sys
import pandas as pd
from scipy.stats import mannwhitneyu

D = "/home/unnot/vault/kjm/miseq/dada2_output/"
G = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD"]
asv = pd.read_csv(D + "asv_table.tsv", sep="\t", index_col=0).drop(columns="sequence")
tax = pd.read_csv(D + "taxonomy.tsv", sep="\t", index_col=0)
meta = pd.read_csv(D + "metadata.tsv", sep="\t", index_col=0)
meta = meta[(meta.experiment == "ARG") & meta.group.isin(G)]
meta = meta.loc[meta.index.intersection(asv.columns)]
if "--nofilter" not in sys.argv:
    meta = meta[asv[meta.index].sum() >= 5000]
g = asv[meta.index].groupby(tax.loc[asv.index, "Genus"]).sum()          # NaN genus dropped
pct = g.div(g.sum(), axis=1) * 100
GEN = ["Lactococcus", "Faecalibaculum", "Bacteroides", "Ligilactobacillus", "Akkermansia"]
long = pct.loc[GEN].T.join(meta[["group", "timepoint"]])
st = long.groupby(["group", "timepoint"])[GEN].agg(["mean", "std", "count"])
if "--nofilter" not in sys.argv:
    st.round(2).to_csv(D + "diversity/genus_group_stats.tsv", sep="\t")


def q(gen, grp, tp):
    v = long[(long.group == grp) & (long.timepoint == tp)][gen]
    return f"{v.mean():.1f} +/- {v.std():.1f} (n={len(v)})"


print(f"samples: {len(meta)}")
for gen, tp, grps in [("Lactococcus", "T72h", G), ("Faecalibaculum", "T72h", ["AMX_CTL", "AMX_HFD"]),
                      ("Bacteroides", "T72h", ["ERY_CTL", "ERY_HFD"]),
                      ("Lactococcus", "T0h", G), ("Lactococcus", "SD14", G),
                      ("Akkermansia", "SD14", G), ("Bacteroides", "SD14", ["AMX_HFD"])]:
    for grp in grps:
        print(f"{gen:15s} {tp:5s} {grp:8s} {q(gen, grp, tp)}")
for gen in GEN:
    m = long[long.timepoint == "T0h"].groupby("group")[gen].mean()
    print(f"T0 range of group means {gen}: ND {m[[x for x in G if 'HFD' not in x]].min():.1f}-{m[[x for x in G if 'HFD' not in x]].max():.1f}"
          f"  HFD {m[[x for x in G if 'HFD' in x]].min():.1f}-{m[[x for x in G if 'HFD' in x]].max():.1f}  all {m.min():.1f}-{m.max():.1f}")
for tp in ("T0h", "SD14"):
    a = long[(long.group == "CTL") & (long.timepoint == tp)]["Lactococcus"]
    b = long[(long.group == "HFD") & (long.timepoint == tp)]["Lactococcus"]
    print(f"untreated Lactococcus {tp}: Mann-Whitney exact p = {mannwhitneyu(a, b, method='exact').pvalue:.2g}")
