#!/usr/bin/env python3
"""Shotgun-metagenome butyrate-synthesis gene analysis (KJM).

Inputs (per sample, from run_search.sh and the existing assembly pipeline):
  hmm_out/<s>.tbl          hmmsearch tblout vs 36 KOfam profiles
  diamond_out/<s>.tsv      DIAMOND best hit vs But positive/negative reference
  results/<s>/<s>_depth.txt  contig mean depth (jgi_summarize_bam_contig_depths)
Gene abundance = summed depth of contigs carrying the gene, divided by the
median summed depth of 12 single-copy marker genes (= copies per genome).
"""
import glob, os, re
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

B = "/home/unnot/vault/kjm/hiseq/butyrate_analysis"
R = "/home/unnot/vault/kjm/hiseq/results"
MAG = "/home/unnot/vault/kjm/hiseq/mag_analysis"
RAW = "/home/unnot/vault/kjm/Raw_data_bio-prameter.xlsx"
OUT = f"{B}/tables"; os.makedirs(OUT, exist_ok=True)

import sys
BUT_MIN_ID = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
BUT_MIN_COV = 0.80
if len(sys.argv) > 1 and BUT_MIN_ID != 60.0:
    OUT = f"{B}/tables_sensitivity_id{int(BUT_MIN_ID)}"; os.makedirs(OUT, exist_ok=True)
# curated positive But references (documented butyrate producers); all others act as negatives
CUR = pd.read_csv(f"{B}/ref/but_reference_curation.tsv", sep="\t")
POS_ACC = set(CUR.loc[CUR.curated_positive, "Entry"])

MARKERS = ["K02886", "K02906", "K02926", "K02931", "K02933", "K02874", "K02878",
           "K02982", "K02988", "K02946", "K02994", "K02996"]
GENES = {  # gene label -> (pathway, KOs)
    "thl":  ("Acetyl-CoA", ["K00626"]),
    "hbd":  ("Acetyl-CoA", ["K00074"]),
    "crt":  ("Acetyl-CoA", ["K01715", "K17865"]),
    "bcd":  ("Acetyl-CoA", ["K00248"]),
    "etfA": ("Acetyl-CoA", ["K03522"]),
    "etfB": ("Acetyl-CoA", ["K03521"]),
    "kamA": ("Lysine", ["K01843"]),
    "kamE": ("Lysine", ["K18011"]),
    "kdd":  ("Lysine", ["K18012"]),
    "kce":  ("Lysine", ["K18013"]),
    "kal":  ("Lysine", ["K18014"]),
    "4hbt": ("4-Aminobutyrate", ["K18122"]),
    "abfD": ("4-Aminobutyrate", ["K14534"]),
    "4hbd": ("4-Aminobutyrate", ["K00043", "K18120"]),
    "gctA": ("Glutarate", ["K01039"]),
    "gctB": ("Glutarate", ["K01040"]),
    "hgdA": ("Glutarate", ["K20903"]),
    "hgdB": ("Glutarate", ["K20904"]),
    "gcdA": ("Glutarate", ["K01615"]),
    "ptb":  ("Terminal", ["K00634"]),
    "buk":  ("Terminal", ["K00929"]),
    "but":  ("Terminal", []),  # from DIAMOND
    "buk_orphan": ("Other", []),  # buk-like kinase without adjacent ptb (not counted)
}
KO2GENE = {k: g for g, (_, kos) in GENES.items() for k in kos}

# ---------------- metadata ----------------
GROUP = {"AC": "AMX_CTL", "AH": "AMX_HFD", "CC": "CIP_CTL", "CH": "CIP_HFD",
         "EC": "ERY_CTL", "EH": "ERY_HFD", "CTL": "CTL", "HFD": "HFD",
         "CSL": "CTL"}  # CSL = CTL mice sampled at recovery timepoints (16S metadata)
def parse(s):
    m = re.match(r"ARG_([A-Z]+)(\d+)_(\w+)$", s)
    pre, num, tp = m.groups()
    g = GROUP[pre]
    return dict(sample=s, group=g, mouse=f"{g}_{num}", timepoint=tp,
                diet="HFD" if g.endswith("HFD") else "ND",
                antibiotic=g.split("_")[0] if "_" in g else "None")

samples = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{B}/hmm_out/*.tbl"))
meta = pd.DataFrame([parse(s) for s in samples]).set_index("sample")

# ---------------- thresholds ----------------
ko = pd.read_csv(f"{B}/ref/ko_thresholds.tsv", sep="\t", header=None)
THR = {r[0]: (float(r[1]), r[2]) for r in ko.itertuples(index=False)}

def contig_of(prot):
    return prot.rsplit("_", 1)[0]

# ---------------- MAG taxonomy (contig -> genus) ----------------
tax = pd.read_csv(f"{MAG}/all_taxonomy.tsv", sep="\t").set_index("bin_global_id")
contig2mag = {}
for fa in glob.glob(f"{MAG}/all_mq_bins/*.fa"):
    mag = os.path.basename(fa)[:-3]
    smp = mag.split("__")[0]
    with open(fa) as fh:
        for l in fh:
            if l.startswith(">"):
                contig2mag[(smp, l[1:].split()[0])] = mag

rows, gene_hits = [], []
for s in samples:
    depth = pd.read_csv(f"{R}/{s}/{s}_depth.txt", sep="\t", usecols=[0, 2])
    depth = dict(zip(depth.iloc[:, 0], depth.iloc[:, 1]))
    # --- KOfam hits (best KO per protein, KOfam adaptive thresholds) ---
    best = {}
    for l in open(f"{B}/hmm_out/{s}.tbl"):
        if l.startswith("#"):
            continue
        f = l.split()
        prot, k, full, dom = f[0], f[2], float(f[5]), float(f[8])
        thr, typ = THR[k]
        sc = dom if typ == "domain" else full
        if sc < thr:
            continue
        ratio = sc / thr
        if prot not in best or ratio > best[prot][1]:
            best[prot] = (k, ratio)
    assign = {p: (KO2GENE.get(k, k)) for p, (k, _) in best.items()}
    # --- but from DIAMOND (best hit must be a POS reference) ---
    if os.path.getsize(f"{B}/diamond_out/{s}.tsv") > 0:
        dm = pd.read_csv(f"{B}/diamond_out/{s}.tsv", sep="\t", header=None,
                         names="q s pid len qlen slen ev bits".split())
        dm = dm.sort_values("bits", ascending=False).drop_duplicates("q")
        dm["acc"] = dm.s.str.split("|").str[2]
        ok = dm[dm.s.str.startswith("POS|") & dm.acc.isin(POS_ACC) & (dm.pid >= BUT_MIN_ID)]
        ok = ok[(ok.len / (ok.qlen) >= BUT_MIN_COV) | (ok.len / ok.slen >= BUT_MIN_COV)]
        for q in ok.q:
            assign[q] = "but"  # overrides 4hbt/other CoA-transferase KO calls
    # --- buk counts only when ptb is on the same contig (ptb-buk operon) ---
    ptb_contigs = {contig_of(p) for p, g in assign.items() if g == "ptb"}
    for p, g in list(assign.items()):
        if g == "buk" and contig_of(p) not in ptb_contigs:
            assign[p] = "buk_orphan"
    # --- abundance ---
    sums = defaultdict(float)
    for p, g in assign.items():
        d = depth.get(contig_of(p), 0.0)
        sums[g] += d
        if g in GENES:
            mag = contig2mag.get((s, contig_of(p)))
            genus = tax.loc[mag, "genus"] if mag in tax.index else None
            gene_hits.append(dict(sample=s, protein=p, gene=g, depth=d, mag=mag,
                                  genus=genus if isinstance(genus, str) and genus else
                                  ("Unclassified MAG" if mag else "Unbinned")))
    marker = np.median([sums.get(k, 0.0) for k in MARKERS])
    row = dict(sample=s, genome_equivalents_depth=marker)
    for g in GENES:
        row[g] = sums.get(g, 0.0) / marker if marker > 0 else np.nan
    rows.append(row)

ab = pd.DataFrame(rows).set_index("sample").join(meta)
ab["terminal_but_buk"] = ab["but"] + ab["buk"]
ab.to_csv(f"{OUT}/butyrate_gene_abundance_per_sample.tsv", sep="\t")
hits = pd.DataFrame(gene_hits).join(meta, on="sample")
hits.to_csv(f"{OUT}/butyrate_gene_hits_with_taxonomy.tsv", sep="\t", index=False)

# ---------------- statistics ----------------
TP = ["T0h", "T72h", "SD14"]
GRP = ["CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD"]
KEY = ["but", "buk", "terminal_but_buk", "bcd", "ptb"]
eps = 1e-4
summ = (ab.groupby(["timepoint", "group"])[KEY]
          .agg(["mean", "std"]).round(4))
summ.to_csv(f"{OUT}/butyrate_gene_summary_by_group.tsv", sep="\t")

lines = []
for gene in KEY:
    for tp in TP:
        d = ab[ab.timepoint == tp].copy()
        d["y"] = np.log10(d[gene] + eps)
        m = smf.ols("y ~ C(diet) * C(antibiotic)", data=d).fit()
        a = anova_lm(m, typ=2)
        for term in a.index[:-1]:
            lines.append(dict(gene=gene, timepoint=tp, term=term,
                              F=round(a.loc[term, "F"], 2), p=a.loc[term, "PR(>F)"]))
anova = pd.DataFrame(lines)
anova.to_csv(f"{OUT}/butyrate_twoway_anova.tsv", sep="\t", index=False)

# Within-mouse change T0h -> T72h -> SD14 (paired, all mice and per antibiotic)
paired = []
for gene in KEY:
    w = ab.pivot_table(index="mouse", columns="timepoint", values=gene)
    w = w.join(ab.groupby("mouse")[["group", "antibiotic", "diet"]].first())
    for sub, dd in [("all", w)] + [(a, w[w.antibiotic == a]) for a in ["AMX", "CIP", "ERY", "None"]]:
        for tp in ["T72h", "SD14"]:
            x = dd[["T0h", tp]].dropna()
            if len(x) >= 3:
                st = stats.wilcoxon(x["T0h"], x[tp]) if len(x) >= 6 else None
                paired.append(dict(gene=gene, subset=sub, comparison=f"T0h vs {tp}", n=len(x),
                                   median_T0h=x["T0h"].median(), median_tp=x[tp].median(),
                                   fold_change=(x[tp].mean() + eps) / (x["T0h"].mean() + eps),
                                   wilcoxon_p=st.pvalue if st else np.nan))
pd.DataFrame(paired).to_csv(f"{OUT}/butyrate_paired_change.tsv", sep="\t", index=False)

# ---------------- link to measured cecal SCFA (group means; end of experiment = SD14) ----------------
sc = pd.read_excel(RAW, sheet_name="SCFA")
scfa = {}
for gcol, vcol, name in [("group", "value(μM)", "acetate"), ("group.1", "value(μM).1", "propionate"),
                         ("group.2", "value(μM).2", "butyrate")]:
    d = sc[[gcol, vcol]].dropna()
    d["g"] = d[gcol].str.replace(r"^\d+_", "", regex=True)
    scfa[name] = d.groupby("g")[vcol].mean()
scfa = pd.DataFrame(scfa)
sd14 = ab[ab.timepoint == "SD14"].groupby("group")[KEY].mean()
link = sd14.join(scfa)
link.to_csv(f"{OUT}/sd14_gene_vs_cecal_scfa_group_means.tsv", sep="\t")
cor = []
for gene in KEY:
    for acid in ["butyrate", "acetate", "propionate"]:
        r, p = stats.spearmanr(link[gene], link[acid])
        pr, pp = stats.pearsonr(link[gene], link[acid])
        cor.append(dict(gene=gene, scfa=acid, n_groups=len(link), spearman_rho=round(r, 3),
                        spearman_p=round(p, 4), pearson_r=round(pr, 3), pearson_p=round(pp, 4)))
pd.DataFrame(cor).to_csv(f"{OUT}/sd14_gene_vs_scfa_correlation.tsv", sep="\t", index=False)

# ---------------- taxonomy of but/buk carriers ----------------
tx = hits[hits.gene.isin(["but", "buk"])].copy()
tx = tx.merge(ab[["genome_equivalents_depth"]], left_on="sample", right_index=True)
tx["copies_per_genome"] = tx["depth"] / tx["genome_equivalents_depth"]
gen = (tx.groupby(["gene", "genus"])["copies_per_genome"].sum() /
       tx.groupby("gene")["copies_per_genome"].sum()).rename("share").reset_index()
gen.sort_values(["gene", "share"], ascending=[True, False]).to_csv(
    f"{OUT}/but_buk_carrier_genera.tsv", sep="\t", index=False)
tg = (tx[tx.gene == "but"].groupby(["timepoint", "group", "genus"])["copies_per_genome"].sum()
      / 3).rename("mean_copies_per_genome").reset_index()
tg.to_csv(f"{OUT}/but_carrier_genera_by_group.tsv", sep="\t", index=False)

# ---------------- console report ----------------
pd.set_option("display.width", 200)
print("Samples:", len(ab), "| groups x timepoints:", ab.groupby(["group", "timepoint"]).size().unstack().to_dict())
print("\nMarker depth (genome equivalents) range:", ab.genome_equivalents_depth.min().round(1), "-", ab.genome_equivalents_depth.max().round(1))
print("\nMean copies/genome by group x timepoint:")
for gene in ["but", "buk", "bcd"]:
    print(f"\n[{gene}]"); print(ab.pivot_table(index="group", columns="timepoint", values=gene).reindex(GRP)[TP].round(3))
print("\nTwo-way ANOVA (log10):"); print(anova[anova.p < 0.1].to_string())
print("\nSD14 group means vs cecal SCFA:"); print(link.round(3).to_string())
print(pd.DataFrame(cor).to_string())
print("\nTop but carriers:"); print(gen[gen.gene == "but"].sort_values("share", ascending=False).head(12).to_string())
print("\nTop buk carriers:"); print(gen[gen.gene == "buk"].sort_values("share", ascending=False).head(8).to_string())
print("\nPaired change:"); print(pd.DataFrame(paired).query("gene in ['but','terminal_but_buk']").round(4).to_string())
