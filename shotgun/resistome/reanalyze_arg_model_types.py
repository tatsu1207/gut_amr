#!/usr/bin/env python3
"""Recompute read-level ARG abundance after removing CARD variant models.

Bowtie2 + samtools idxstats counts every read mapping to a CARD reference as an ARG.
For "rRNA gene variant model" and "protein variant model" entries, resistance is
conferred by a point mutation, so wild-type reads are counted as resistance genes.
This inflates totals and can create spurious signal (notably Spne_23S_* in the MLS class).

Outputs, in arg_results/model_filtered/:
  arg_abundance_per_gene_<level>.tsv      filtered per-gene table
  arg_abundance_per_class_<level>.tsv     drug-class totals
  total_arg_by_sample.tsv                 totals per sample at each filter level
  comparison_stats.tsv                    HFD vs ND per timepoint at each filter level
where <level> is all / no_rrna / homolog_only.
"""
import json, os, re
import pandas as pd
from scipy import stats

H = "/home/unnot/vault/kjm/hiseq"
OUT = f"{H}/arg_results/model_filtered"
os.makedirs(OUT, exist_ok=True)

# ---- CARD model type per ARO ----
card = json.load(open(f"{H}/card_db/localDB/card.json"))
aro2model = {}
for k, v in card.items():
    if isinstance(v, dict) and "model_type" in v and "ARO_accession" in v:
        aro2model[str(v["ARO_accession"])] = v["model_type"]

d = pd.read_csv(f"{H}/arg_results/arg_abundance_per_gene.tsv", sep="\t")
META = ["gene", "gene_name", "aro", "drug_class", "resistance_mechanism", "gene_family"]
SAMP = [c for c in d.columns if c not in META]

d["aro_str"] = d["aro"].astype(str).str.replace(r"^ARO:", "", regex=True).str.strip()
d["model_type"] = d["aro_str"].map(aro2model)
# fall back to the gene_family text where the ARO lookup misses
fam_rrna = d["gene_family"].astype(str).str.contains("rRNA with mutation", case=False)
d.loc[d["model_type"].isna() & fam_rrna, "model_type"] = "rRNA gene variant model"
print("model_type coverage:", d["model_type"].notna().mean().round(3))
print(d["model_type"].value_counts(dropna=False).to_string())

LEVELS = {
    "all":          lambda df: pd.Series(True, index=df.index),
    "no_rrna":      lambda df: df["model_type"] != "rRNA gene variant model",
    "homolog_only": lambda df: ~df["model_type"].isin(
                        ["rRNA gene variant model", "protein variant model",
                         "protein overexpression model", "protein knockout model"]),
}

GM = {"AC": ("AMX_CTL", "ND", "AMX"), "AH": ("AMX_HFD", "HFD", "AMX"),
      "CC": ("CIP_CTL", "ND", "CIP"), "CH": ("CIP_HFD", "HFD", "CIP"),
      "EC": ("ERY_CTL", "ND", "ERY"), "EH": ("ERY_HFD", "HFD", "ERY"),
      "CTL": ("CTL", "ND", "None"), "HFD": ("HFD", "HFD", "None"),
      "CSL": ("CTL", "ND", "None")}

def meta_of(s):
    m = re.match(r"ARG_([A-Za-z]+?)(\d+)_(\w+)$", s)
    pre, mouse, tp = m.groups()
    g, diet, abx = GM[pre]
    return dict(sample=s, group=g, diet=diet, antibiotic=abx, timepoint=tp, mouse=f"{g}_{mouse}")

totals, stats_rows = {}, []
for level, fn in LEVELS.items():
    keep = fn(d)
    sub = d[keep]
    sub.drop(columns=["aro_str"]).to_csv(f"{OUT}/arg_abundance_per_gene_{level}.tsv", sep="\t", index=False)
    cls = sub.groupby(sub["drug_class"].astype(str))[SAMP].sum()
    cls.to_csv(f"{OUT}/arg_abundance_per_class_{level}.tsv", sep="\t")
    tot = sub[SAMP].sum()
    totals[level] = tot
    print(f"\n[{level}] genes kept {keep.sum()}/{len(d)} | total abundance {tot.sum():.1f}")

    md = pd.DataFrame([meta_of(s) for s in SAMP]).set_index("sample")
    md["total"] = tot
    for tp in ["T0h", "T72h", "SD14"]:
        x = md[md.timepoint == tp]
        a, b = x[x.diet == "HFD"].total, x[x.diet == "ND"].total
        u = stats.mannwhitneyu(a, b)
        stats_rows.append(dict(level=level, timepoint=tp, n_HFD=len(a), n_ND=len(b),
                               mean_HFD=round(a.mean(), 3), sd_HFD=round(a.std(), 3),
                               mean_ND=round(b.mean(), 3), sd_ND=round(b.std(), 3),
                               ratio=round(a.mean() / b.mean(), 2), p=u.pvalue))
        print(f"   {tp}: HFD {a.mean():.3f} vs ND {b.mean():.3f}  ratio {a.mean()/b.mean():.2f}  p={u.pvalue:.4f}")

md = pd.DataFrame([meta_of(s) for s in SAMP]).set_index("sample")
for level, tot in totals.items():
    md[f"total_{level}"] = tot
md.to_csv(f"{OUT}/total_arg_by_sample.tsv", sep="\t")
st = pd.DataFrame(stats_rows)
st.to_csv(f"{OUT}/comparison_stats.tsv", sep="\t", index=False)

print("\n=== group means by level ===")
for level in LEVELS:
    piv = md.pivot_table(index="group", columns="timepoint", values=f"total_{level}")
    print(f"\n[{level}]"); print(piv[["T0h", "T72h", "SD14"]].round(3).to_string())
print("\nwrote", OUT)
