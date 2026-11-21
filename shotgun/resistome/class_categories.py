"""Drug-class categories for the read-level composition (Fig. 4B and the class shares in the text).

Rule (2026-11-21): CARD drug classes are grouped into families (MLS = macrolide, lincosamide,
streptogramin, pleuromutilin; beta-lactam = all beta-lactam subclasses). A gene whose CARD drug-class
string spans three or more families -- in these data mostly the chromosomal multidrug efflux systems of
Escherichia (mdtEF, acrEF, acrB, tolC and their regulators) -- is shown as one "multidrug (>=3 classes)"
category instead of being split across classes, which made efflux pumps look like class-specific
resistance. Every other gene is split equally among its display classes, as before, so the categories
of a sample still sum exactly to its total. Display classes are deduplicated first, so a gene is divided
among the distinct classes it is shown under, not among CARD's subclass terms.
"""
import pandas as pd

FAMILY = {"macrolide antibiotic": "MLS", "lincosamide antibiotic": "MLS",
          "streptogramin antibiotic": "MLS", "streptogramin A antibiotic": "MLS",
          "streptogramin B antibiotic": "MLS", "pleuromutilin antibiotic": "MLS"}
BETA = ("penicillin beta-lactam", "cephalosporin", "carbapenem", "monobactam", "cephamycin", "penem")
DISPLAY = {"tetracycline antibiotic": "tetracycline", "lincosamide antibiotic": "lincosamide",
           "macrolide antibiotic": "macrolide", "streptogramin antibiotic": "streptogramin",
           "streptogramin A antibiotic": "streptogramin", "streptogramin B antibiotic": "streptogramin",
           "pleuromutilin antibiotic": "pleuromutilin", "glycopeptide antibiotic": "glycopeptide",
           "fluoroquinolone antibiotic": "fluoroquinolone", "aminoglycoside antibiotic": "aminoglycoside",
           "phenicol antibiotic": "phenicol", "peptide antibiotic": "peptide"}
ORDER = ["tetracycline", "lincosamide", "macrolide", "streptogramin", "pleuromutilin", "glycopeptide",
         "fluoroquinolone", "beta-lactam", "aminoglycoside", "phenicol", "peptide", "other",
         "multidrug (>=3 classes)"]
MULTI = "multidrug (>=3 classes)"


def family(c):
    if c in FAMILY:
        return FAMILY[c]
    if c in BETA:
        return "beta-lactam"
    return c


def display(c):
    if c in BETA:
        return "beta-lactam"
    return DISPLAY.get(c, "other")


def categorise(pergene, samples):
    """Return a samples x categories DataFrame (columns in ORDER) that sums to each sample's total."""
    out = {k: 0 for k in ORDER}
    out = {k: pd.Series(0.0, index=samples) for k in ORDER}
    for _, r in pergene.iterrows():
        cls = [c.strip() for c in str(r["drug_class"]).split(";") if c.strip()]
        v = r[samples].astype(float)
        if len({family(c) for c in cls}) >= 3:
            out[MULTI] = out[MULTI] + v
            continue
        # Deduplicate display classes before splitting (2026-11-21): CARD lists several terms that map to
        # one display class (streptogramin / streptogramin A / streptogramin B; the beta-lactam subclasses),
        # and counting each term separately gave that class extra weight -- ErmF, a single MLS gene, was
        # split 3/5 streptogramin, 1/5 macrolide, 1/5 lincosamide instead of one third each.
        disp = sorted({display(c) for c in cls})
        for d in disp:
            out[d] = out[d] + v / len(disp)
    return pd.DataFrame(out)[ORDER]


def is_multidrug(drug_class):
    cls = [c.strip() for c in str(drug_class).split(";") if c.strip()]
    return len({family(c) for c in cls}) >= 3
