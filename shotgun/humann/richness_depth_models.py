#!/usr/bin/env python3
"""Pathway richness: the models quoted in the Results, and the sequencing-depth sensitivity analysis.

Pathway richness is a detection-limited statistic, so it rises with effective sequencing depth. This
script quantifies that dependence, refits the per-timepoint models with log10 depth as a covariate, and
reports the depth balance of each diet contrast, so that the claims in the text are restricted to the
contrasts that depth cannot explain.

Factor coding: the antibiotic factor has FOUR levels - AMX, CIP, ERY and None, the last being the
untreated CTL and HFD groups. Dropping the untreated groups (a three-level factor on 18 samples) gives
materially different p-values and answers a different question; both are reported below so the
distinction is explicit.

Input:  ../../../KJM_info/manuscript/supplementary_tables.xlsx  sheet SD15_humann_pathways
Output: tables/richness_depth_models.tsv  (sheet SD15b of the supplementary workbook)
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

H = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(H, "../../../KJM_info/manuscript/supplementary_tables.xlsx")
OUT = os.path.join(H, "tables")
os.makedirs(OUT, exist_ok=True)

UNTREATED = "NONE"   # the label the published sheet uses for the CTL and HFD groups

d = pd.read_excel(XLSX, "SD15_humann_pathways", index_col=0)
# The published sheet labels the untreated rows "NONE". That spelling matters: the literal string
# "None" is in pandas' default na_values, so a sheet holding "None" (or a blank cell) is read back as
# NaN and statsmodels then DROPS those rows, silently fitting a three-level antibiotic factor on 18
# samples. fillna is kept so the script also runs against older copies of the workbook.
d["abx"] = d["antibiotic"].fillna(UNTREATED).replace({"None": UNTREATED, "none": UNTREATED})
assert d["abx"].nunique() == 4, f"expected 4 antibiotic levels, got {sorted(d['abx'].unique())}"
d["logp"] = np.log10(d["n_pathways"])
d["logd"] = np.log10(d["genome_equivalents_depth"])

rows = []


def fmt_p(p):
    """Keep three significant digits. round(p, 4) turned 2.55e-5 into 0, which reads as p = 0."""
    return float(f"{float(p):.3g}")


def add(timepoint, model, term, p, note=""):
    rows.append(dict(timepoint=timepoint, model=model, term=term, p_value=fmt_p(p), note=note))


# --- how strongly does richness depend on depth? ------------------------------------------------
print("Spearman correlation between genome equivalents and pathway richness")
for tp, sub in [("all", d)] + [(t, d[d.timepoint == t]) for t in ("T0", "T3", "R14")]:
    rho, p = spearmanr(sub["genome_equivalents_depth"], sub["n_pathways"])
    print(f"  {tp:4s} n={len(sub):2d}  rho={rho:5.2f}  p={p:.2g}")
    add(tp, "Spearman depth vs richness", "rho", p, f"rho = {rho:.2f}, n = {len(sub)}")

# --- the models, with and without the depth covariate -------------------------------------------
MODELS = [("two-way ANOVA (reported)", "n_pathways ~ C(abx)*C(diet)", None),
          ("two-way ANOVA + log10 depth", "n_pathways ~ logd + C(abx)*C(diet)", None),
          ("two-way ANOVA, untreated dropped", "n_pathways ~ C(abx)*C(diet)", "drop")]
print("\nTwo-way ANOVA (type II) per timepoint")
for tp in ("T0", "T3", "R14"):
    for label, formula, how in MODELS:
        sub = d[d.timepoint == tp]
        if how == "drop":
            sub = sub[sub.abx != UNTREATED]
        a = sm.stats.anova_lm(smf.ols(formula, data=sub).fit(), typ=2)
        ps = {t: a.loc[t, "PR(>F)"] for t in a.index if t != "Residual"}
        print(f"  {tp:4s} {label:34s} n={len(sub):2d}  " +
              "  ".join(f"{k.replace('C(','').replace(')',''):14s}{v:.4f}" for k, v in ps.items()))
        for k, v in ps.items():
            add(tp, label, k.replace("C(", "").replace(")", ""), v, f"n = {len(sub)}")

# --- is each diet contrast balanced for depth? --------------------------------------------------
print("\nDepth balance of the ND vs HFD contrast within each antibiotic, at T3")
t3 = d[d.timepoint == "T3"]
for ab in ("AMX", "CIP", "ERY", UNTREATED):
    s = t3[t3.abx == ab]
    nd, hf = s[s.diet == "ND"], s[s.diet == "HFD"]
    dn, dh = nd["genome_equivalents_depth"].mean(), hf["genome_equivalents_depth"].mean()
    ratio = max(dn, dh) / min(dn, dh)
    verdict = "matched" if ratio < 1.15 else "NOT matched"
    print(f"  {ab:5s} depth {dn:7.0f} vs {dh:7.0f} ({ratio:.2f}x, {verdict})   "
          f"richness {nd['n_pathways'].mean():5.0f} vs {hf['n_pathways'].mean():5.0f}")
    rows.append(dict(timepoint="T3", model="depth balance of ND vs HFD", term=ab, p_value=np.nan,
                     note=f"depth ND {dn:.0f} vs HFD {dh:.0f} ({ratio:.2f}x, {verdict}); "
                          f"richness {nd['n_pathways'].mean():.0f} vs {hf['n_pathways'].mean():.0f}"))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(OUT, "richness_depth_models.tsv"), sep="\t", index=False)
print(f"\nwrote {OUT}/richness_depth_models.tsv ({len(out)} rows)")
