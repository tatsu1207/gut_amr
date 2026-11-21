#!/usr/bin/env python3
"""Assembly ARG-carrying ORFs by CARD model type (numbers in the Methods, 2026-11-21)."""
import glob
import pandas as pd
rows = []
for f in glob.glob("/home/unnot/vault/kjm/hiseq/results/*/*_rgi_main.txt"):
    d = pd.read_csv(f, sep="\t", usecols=["Model_type", "Best_Hit_ARO", "SNPs_in_Best_Hit_ARO"])
    s = f.split("/")[-2]; d["sample"] = s; d["tp"] = s.rsplit("_", 1)[1]; d["grp"] = s.split("_")[1].rstrip("0123456789")
    rows.append(d)
r = pd.concat(rows)
print("total", len(r)); print((r.Model_type.value_counts() / len(r) * 100).round(1).to_string())
v = r[r.Model_type == "protein variant model"]
print("variant share by timepoint (%):", (v.groupby("tp").size() / r.groupby("tp").size() * 100).round(1).to_dict())
print("rpoC (M. tuberculosis model) hits:", v.Best_Hit_ARO.str.contains("rpoC").sum(), "of", len(v))
q = r[r.Best_Hit_ARO.str.contains(r"gyrA|parC", case=False)]
print("gyrA/parC hits (any model):", len(q), q.groupby("grp").size().to_dict())
x = v[v.Best_Hit_ARO.str.contains("gyrB|PBP3|EF-Tu|UhpT")]
print("drug-relevant variants by group/timepoint:", x.groupby(["grp", "tp"]).size().to_dict())
