#!/usr/bin/env python3
"""Reconcile the HUMAnN 3.9 read-level butyrate signal with the assembly-based measure of Fig. 6A.

Background: Fig. 6A quantifies butyrate synthesis from the assemblies, as copies per genome equivalent
of a 20-gene pathway, and reports the terminal step (but, EC 2.8.3.8; buk, EC 2.7.2.7) as
terminal_but_buk (Supplementary Data 12a). HUMAnN gives a read-level alternative. The two agree at T3
and disagree at R14. This script shows why, so that the choice of measure in the manuscript is
documented rather than assumed.

Inputs (built by run_humann.sh then humann_join_tables / humann_renorm_table / humann_regroup_table):
    joined/pathabundance_cpm.tsv   MetaCyc pathways, CPM, stratified
    joined/ec_cpm.tsv              level-4 EC, CPM, stratified
    ../../miseq/dada2_output/metadata.tsv
    ../../../KJM_info/manuscript/supplementary_tables.xlsx  sheet SD12a_butyrate_abundance

Writes tables/humann_butyrate_reconcile.tsv and prints the numbers quoted in STATUS.md rev. 10.29.
"""
import os
import pandas as pd
from scipy.stats import spearmanr

H = os.path.dirname(os.path.abspath(__file__))
META = os.path.join(H, "../../miseq/dada2_output/metadata.tsv")
XLSX = os.path.join(H, "../../../KJM_info/manuscript/supplementary_tables.xlsx")

# MetaCyc pathways that end in butanoate
BUT_PWY = ["CENTFERM-PWY: pyruvate fermentation to butanoate",
           "PWY-5676: acetyl-CoA fermentation to butanoate II",
           "P163-PWY: L-lysine fermentation to acetate and butanoate",
           "PWY-5677: succinate fermentation to butanoate"]
# terminal enzymes of the butyrate pathway, as level-4 EC
EC = {"2.8.3.8": "but", "2.7.2.7": "buk", "2.3.1.19": "ptb", "1.3.8.1": "bcd"}


def load(path):
    d = pd.read_csv(path, sep="\t", index_col=0)
    d.columns = [c.replace("_Abundance", "").replace("-CPM", "") for c in d.columns]
    return d


def unstratified(d):
    return d[~d.index.str.contains(r"\|")]


p = load(os.path.join(H, "joined/pathabundance_cpm.tsv"))
e = load(os.path.join(H, "joined/ec_cpm.tsv"))
up, ue = unstratified(p), unstratified(e)
meta = pd.read_csv(META, sep="\t", index_col=0)


def ec_sum(key):
    rows = [i for i in ue.index if i.split(":")[0] == key]
    return ue.loc[rows].sum() if rows else pd.Series(0.0, index=ue.columns)


h = pd.DataFrame({f"humann_{v}_cpm": ec_sum(k) for k, v in EC.items()})
h["humann_terminal_cpm"] = h["humann_but_cpm"] + h["humann_buk_cpm"]
h["humann_but_pathways_cpm"] = up.loc[[x for x in BUT_PWY if x in up.index]].sum()
richness = (up.drop(index=["UNMAPPED", "UNINTEGRATED"], errors="ignore") > 0).sum()
h["n_pathways"] = richness

asm = pd.read_excel(XLSX, "SD12a_butyrate_abundance", index_col=0)
m = h.join(asm[["group", "diet", "antibiotic", "timepoint", "mouse",
                "terminal_but_buk", "genome_equivalents_depth"]], how="inner")
os.makedirs(os.path.join(H, "tables"), exist_ok=True)
m.round(4).to_csv(os.path.join(H, "tables/humann_butyrate_reconcile.tsv"), sep="\t")

print(f"n = {len(m)} samples\n")
print("### pathway richness (mean unstratified MetaCyc pathways with CPM > 0)")
print(m.pivot_table(index="group", columns="timepoint", values="n_pathways", aggfunc="mean").round(0).to_string())

print("\n### agreement with the assembly measure, sample by sample (Spearman, n = 72)")
for col in ("humann_but_cpm", "humann_buk_cpm", "humann_terminal_cpm", "humann_but_pathways_cpm"):
    rho, pv = spearmanr(m[col], m["terminal_but_buk"])
    print(f"  {col:26s} vs terminal_but_buk: rho = {rho:5.2f}  p = {pv:.1e}")

print("\n### the amoxicillin mice: assembly copies/genome vs HUMAnN CPM")
amx = m[m["antibiotic"].astype(str).eq("AMX")].sort_values(["group", "mouse", "timepoint"])
print(amx[["group", "mouse", "timepoint", "terminal_but_buk",
           "humann_but_cpm", "humann_buk_cpm", "humann_but_pathways_cpm"]].round(3).to_string(index=False))

print("\n### the EC-homolog artefact")
ac7 = "ARG_AC7_T72h"
if ac7 in m.index:
    r = m.loc[ac7]
    print(f"  {ac7} (the Enterobacterales bloom): assembly terminal_but_buk = {r.terminal_but_buk:.3f}, "
          f"HUMAnN EC 2.8.3.8 = {r.humann_but_cpm:.1f} CPM.")
    print("  E. coli atoDA (acetoacetyl-CoA transferase) carries EC 2.8.3.8, so the read-level route scores")
    print("  a butyrate terminal enzyme in a sample whose assembly contains none.")
strat = p[p.index.str.contains(r"\|")]
for s in ("ARG_AH3_SD14", "ARG_AH4_SD14", "ARG_AH5_SD14"):
    if s not in strat.columns:
        continue
    sub = strat[strat.index.str.split("|").str[0].isin(BUT_PWY)][s]
    print(f"  {s}: butyrate pathways {m.loc[s, 'humann_but_pathways_cpm']:.1f} CPM total, "
          f"of which {sub[sub > 0].sum():.1f} CPM is attributed to any species.")
print("  An unattributed community-level call is what remains when no single genome carries the pathway.")
