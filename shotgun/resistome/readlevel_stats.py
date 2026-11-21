#!/usr/bin/env python3
"""Every read-level resistome statistic quoted in the KJM manuscript, from one place.

Usage: python3 readlevel_stats.py [homolog_only|no_rrna|all]   (default homolog_only)

Input: arg_results/model_filtered/arg_abundance_per_gene_<level>.tsv (ARG copies per
16S rRNA gene copy, per CARD model, per sample) and the 16S DADA2 genus table.
Totals are sums over genes (never over simple_class columns, which double-count
multi-class genes). Class abundances are the sum over genes whose CARD drug-class
string contains the class. Tests: exact two-sided Mann-Whitney U; three-way ANOVA
(type II) on log10 total abundance with diet x antibiotic x timepoint.
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

LEVEL = sys.argv[1] if len(sys.argv) > 1 else "homolog_only"
# second argument "noAC7": sensitivity analysis without ARG_AC7_T72h, the AMX_CTL T3 sample dominated by
# Enterobacterales (15.7% of pooled read-level abundance); output goes to readlevel_stats_<level>_noAC7.txt
NOAC7 = len(sys.argv) > 2 and sys.argv[2] == "noAC7"
H = "/home/unnot/vault/kjm/hiseq/arg_results/model_filtered/"
D = "/home/unnot/vault/kjm/miseq/dada2_output/"

pg = pd.read_csv(f"{H}arg_abundance_per_gene_{LEVEL}.tsv", sep="\t")
if NOAC7:
    pg = pg.drop(columns=["ARG_AC7_T72h"])
S = [c for c in pg.columns if c.startswith("ARG_")]
GROUP = {"AC": "AMX_CTL", "AH": "AMX_HFD", "CC": "CIP_CTL", "CH": "CIP_HFD", "EC": "ERY_CTL",
         "EH": "ERY_HFD", "CTL": "CTL", "CSL": "CTL", "HFD": "HFD"}


def parse(s):
    core = s.replace("ARG_", "")
    for p in sorted(GROUP, key=len, reverse=True):
        if core.startswith(p) and core[len(p)].isdigit():
            return GROUP[p], s.rsplit("_", 1)[1]
    raise ValueError(s)


meta = pd.DataFrame([parse(s) for s in S], index=S, columns=["group", "tp"])
meta["diet"] = np.where(meta.group.str.endswith("HFD"), "HFD", "ND")
meta["abx"] = meta.group.str.split("_").str[0].where(meta.group.str.contains("_"), "none")
tot = pg[S].sum()


def has(cls):
    return pg.drop_duplicates("gene").drop_duplicates(subset=None).loc[
        pg.drop_duplicates("gene").drop_duplicates(subset=None).drug_class.str.contains(cls)].index


def class_ab(pattern):
    return pg.loc[pg.drug_class.str.contains(pattern), S].sum()


def gm(x, g, tp):
    return x[meta.index[(meta.group == g) & (meta.tp == tp)]]


def mw(a, b):
    return mannwhitneyu(a, b, method="exact").pvalue


out = []
P = out.append
P(f"LEVEL {LEVEL}: {len(pg)} models")
# --- baseline
b = tot[meta.tp == "T0h"]
means = b.groupby(meta.group).mean()
P(f"T0h group-mean range {means.min():.3f}-{means.max():.3f}")
h, n = b[meta.diet == "HFD"], b[meta.diet == "ND"]
P(f"T0h HFD {h.mean():.3f}+/-{h.std():.3f} vs ND {n.mean():.3f}+/-{n.std():.3f}; {h.mean()/n.mean():.2f}x p={mw(h, n):.3f}")
for g in ["AMX", "CIP", "ERY", ""]:
    gh, gn = (f"{g}_HFD", f"{g}_CTL") if g else ("HFD", "CTL")
    P(f"  T0h {gh}/{gn} {gm(tot, gh, 'T0h').mean()/gm(tot, gn, 'T0h').mean():.2f}x")
# --- T72h
ac = gm(tot, "AMX_CTL", "T72h"); base = gm(tot, "AMX_CTL", "T0h").mean()
P(f"AMX_CTL T72h values {sorted(ac.round(2))}; baseline {base:.3f}; mean {ac.mean():.2f} ({ac.mean()/base:.1f}x) median {ac.median():.2f} ({ac.median()/base:.1f}x)")
for g in ["AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD", "CTL", "HFD"]:
    P(f"  T72h {g} {gm(tot, g, 'T72h').mean():.3f} {gm(tot, g, 'T72h').mean()/gm(tot, g, 'T0h').mean():.1f}x")
t = tot[meta.tp == "T72h"]
h, n = t[meta.diet == "HFD"], t[meta.diet == "ND"]
P(f"T72h ND {n.mean():.3f}+/-{n.std():.3f} vs HFD {h.mean():.3f}+/-{h.std():.3f} p={mw(h, n):.2f}")
# --- classes
beta = class_ab("penicillin beta-lactam|cephalosporin|carbapenem|monobactam")
mls = class_ab("macrolide|lincosamide|streptogramin")
pen = class_ab("penicillin beta-lactam")
P(f"AMX_CTL beta-lactam {gm(beta,'AMX_CTL','T0h').mean():.3f}->{gm(beta,'AMX_CTL','T72h').mean():.3f}; "
  f"penicillin log2FC {np.log2(gm(pen,'AMX_CTL','T72h').mean()/gm(pen,'AMX_CTL','T0h').mean()):.1f}; "
  f"MLS {gm(mls,'AMX_CTL','T0h').mean():.3f}->{gm(mls,'AMX_CTL','T72h').mean():.3f}")
for g in ["CIP_HFD", "ERY_HFD"]:
    P(f"  {g} MLS {gm(mls,g,'T0h').mean():.3f}->{gm(mls,g,'T72h').mean():.3f}; T72h MLS/total {gm(mls,g,'T72h').mean()/gm(tot,g,'T72h').mean()*100:.0f}%")
lmrd = pg.loc[pg.gene_name == "lmrD", S].sum()
for g in ["CIP_HFD", "ERY_HFD", "AMX_CTL"]:
    P(f"  lmrD share of MLS T72h {g} {gm(lmrd,g,'T72h').mean()/gm(mls,g,'T72h').mean()*100:.0f}%")
P(f"lmrD share of all MLS (pooled) {lmrd.sum()/mls.sum()*100:.1f}%; of total {lmrd.sum()/tot.sum()*100:.1f}%")
# --- SD14
for g in ["AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD", "CTL", "HFD"]:
    P(f"  SD14 {g} {gm(tot,g,'SD14').mean():.3f} {gm(tot,g,'SD14').mean()/gm(tot,g,'T0h').mean():.2f}x")
s = tot[meta.tp == "SD14"]
a, o = s[meta.abx == "AMX"], s[meta.abx != "AMX"]
P(f"SD14 AMX {a.mean():.3f} vs others {o.mean():.3f} ({a.mean()/o.mean():.1f}x) p={mw(a, o):.2g}; min AMX > max other: {a.min() > o.max()}")
h, n = s[meta.diet == "HFD"], s[meta.diet == "ND"]
P(f"SD14 HFD {h.mean():.3f}+/-{h.std():.3f} vs ND {n.mean():.3f}+/-{n.std():.3f} {h.mean()/n.mean():.2f}x p={mw(h, n):.2f}")
for g in ["AMX", "CIP", "ERY"]:
    x, y = gm(tot, f"{g}_HFD", "SD14"), gm(tot, f"{g}_CTL", "SD14")
    P(f"  SD14 {g} {x.mean()/y.mean():.2f}x p={mw(x, y):.2f}")
uh0, uh14 = gm(tot, "HFD", "T0h").mean(), gm(tot, "HFD", "SD14").mean()
un0, un14 = gm(tot, "CTL", "T0h").mean(), gm(tot, "CTL", "SD14").mean()
P(f"untreated HFD {uh0:.3f}->{uh14:.3f} {uh14/uh0:.2f}x; ND {un0:.3f}->{un14:.3f} {un14/un0:.2f}x; SD14 apart {uh14/un14:.2f}x")
for g in ["AMX", "CIP", "ERY"]:
    hx = gm(tot, f"{g}_HFD", "SD14").mean() - uh14
    nx = gm(tot, f"{g}_CTL", "SD14").mean() - un14
    P(f"  SD14 excess over untreated {g}: HFD {hx:.4f} ND {nx:.4f} adjusted ratio {hx/nx:.2f}")
d = meta.assign(y=np.log10(tot))
an = anova_lm(smf.ols("y ~ C(diet) * C(abx) * C(tp)", d).fit(), typ=2)
for k in an.index[:-1]:
    P(f"  ANOVA {k}: F={an.loc[k,'F']:.2f} p={an.loc[k,'PR(>F)']:.2g}")
an.to_csv(f"{H}threeway_anova_{LEVEL}.tsv", sep="\t")
# --- fractional class shares (Fig 4B)
k = pg.drug_class.str.count(";") + 1
frac = {}
for i, r in pg.iterrows():
    for c in r.drug_class.split(";"):
        frac[c] = frac.get(c, 0) + pg.loc[i, S] / k[i]
frac = pd.DataFrame(frac)
share = frac.sum() / frac.sum().sum() * 100
P("class shares: " + ", ".join(f"{c} {v:.1f}%" for c, v in share.sort_values(ascending=False).head(5).items()))
gmeans = frac.groupby([meta.group, meta.tp]).mean()
lead = gmeans.idxmax(axis=1)
P(f"leading class counts over 24 group-timepoints: {lead.value_counts().to_dict()}")
P(f"  lincosamide-led cells: {list(lead[lead == 'lincosamide antibiotic'].index)}")
# --- Fig. 4B categories (class_categories.py: genes spanning >=3 class families = "multidrug") ---
sys.path.insert(0, "/home/unnot/vault/kjm/hiseq")
from class_categories import categorise
cat = categorise(pg, S)
cshare = cat.sum() / cat.sum().sum() * 100
P("Fig4B category shares: " + ", ".join(f"{c} {v:.1f}%" for c, v in cshare.sort_values(ascending=False).head(7).items()))
clead = cat.groupby([meta.group, meta.tp]).mean().idxmax(axis=1)
P(f"Fig4B leading category counts over 24 group-timepoints: {clead.value_counts().to_dict()}")
P(f"  non-tetracycline-led cells: {[(k, v) for k, v in clead.items() if v != 'tetracycline']}")
# AMX classes rising T0h->SD14 (Fig 4C classes)
simple = pd.read_csv(f"{H}arg_abundance_simple_class_{LEVEL}.tsv", sep="\t", index_col=0)
rise = {}
KEY = ["tetracycline antibiotic", "macrolide antibiotic", "penicillin beta-lactam", "fluoroquinolone antibiotic",
       "lincosamide antibiotic", "phenicol antibiotic", "aminoglycoside antibiotic", "cephalosporin"]
for c in [c for c in KEY if c in simple.columns]:  # same formula as Fig. 4C
    v = []
    for g in ["AMX_CTL", "AMX_HFD"]:
        b0 = simple.loc[gm(tot, g, "T0h").index, c].mean(); b1 = simple.loc[gm(tot, g, "SD14").index, c].mean()
        v.append(np.log2(b1 / b0) if b0 > 1e-4 else np.log2((b1 + 1e-3) / (b0 + 1e-3)))
    rise[c] = np.mean(v)
P("AMX classes rising most T0h->SD14: " + ", ".join(f"{c} {v:.1f}" for c, v in sorted(rise.items(), key=lambda x: -x[1])[:6]))
# --- correlations with 16S genera (genus-assigned denominator, as in the text)
a16 = pd.read_csv(f"{D}asv_table.tsv", sep="\t", index_col=0).drop(columns="sequence")
tx = pd.read_csv(f"{D}taxonomy.tsv", sep="\t", index_col=0)
g16 = tx.loc[a16.index, "Genus"]
a16 = a16[g16.notna()]
gen = (a16 / a16.sum()).groupby(g16[g16.notna()]).sum()
common = [x for x in S if x in gen.columns]
P(f"matched 16S samples {len(common)}")
for g in ["Lactococcus", "Faecalibaculum", "Bifidobacterium", "Bacteroides"]:
    r, p = spearmanr(gen.loc[g, common], tot[common])
    P(f"  rho total ~ {g} {r:+.2f} p={p:.1g}")
tet = class_ab("tetracycline")
for nm, v in [("MLS", mls), ("tetracycline", tet), ("beta-lactam", beta)]:
    r, p = spearmanr(gen.loc["Lactococcus", common], v[common])
    P(f"  rho Lactococcus ~ {nm} {r:+.3f} p={p:.2g}")
other = mls - lmrd
r, p = spearmanr(gen.loc["Lactococcus", common], lmrd[common]); P(f"  rho Lactococcus ~ lmrD {r:+.2f} p={p:.2g}")
r, p = spearmanr(gen.loc["Lactococcus", common], other[common]); P(f"  rho Lactococcus ~ other MLS {r:+.2f} p={p:.2g}")
# --- amoxicillin: per-mouse class totals and what the penicillin beta-lactam signal consists of (added 2026-11-21)
AMXM = [m for m in ["AC6", "AC7", "AC8", "AH3", "AH4", "AH5"] if f"ARG_{m}_T72h" in pg.columns]
def cls_sum(pat, m, tp): return pg[pg.drug_class.str.contains(pat, na=False)][f"ARG_{m}_{tp}"].sum()
for nm, pat in [("MLS", "macrolide|lincosamide|streptogramin"), ("beta-lactam", "penicillin beta-lactam|cephalosporin|carbapenem|monobactam")]:
    P(f"AMX per mouse {nm} T0h -> T72h: " + ", ".join(f"{m} {cls_sum(pat, m, 'T0h'):.4f}->{cls_sum(pat, m, 'T72h'):.4f}" for m in AMXM))
pb = pg[pg.drug_class.str.contains("penicillin beta-lactam", na=False)]
sd = pb[[f"ARG_{m}_SD14" for m in AMXM]].mean(axis=1)
mech = sd.groupby(pb.resistance_mechanism).sum(); mech = mech / mech.sum() * 100
P("AMX SD14 penicillin beta-lactam signal by mechanism (%): " + ", ".join(f"{k} {v:.1f}" for k, v in mech.items()))
if not NOAC7:
  ac7 = pg.set_index("gene_name")["ARG_AC7_T72h"]
  P(f"AC7_T72h share of pooled abundance {100 * ac7.sum() / pg[S].sum().sum():.1f}%; of pooled qnrD "
    f"{100 * (ac7.get('QnrD1', 0) + ac7.get('QnrD2', 0)) / pg[pg.gene_name.isin(['QnrD1', 'QnrD2'])][S].sum().sum():.0f}%")
  P(f"AC7_T72h total {ac7.sum():.3f}; qnrD1+qnrD2 share {100 * (ac7.get('QnrD1', 0) + ac7.get('QnrD2', 0)) / ac7.sum():.1f}%")
print("\n".join(out))
open(f"{H}readlevel_stats_{LEVEL}{'_noAC7' if NOAC7 else ''}.txt", "w").write("\n".join(out) + "\n")
