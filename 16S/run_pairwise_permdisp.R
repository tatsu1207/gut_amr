#!/usr/bin/env Rscript
# Review point R8 (2026-11-20): pairwise ND vs HFD PERMANOVA *with* pairwise PERMDISP,
# within each antibiotic, at the recovery timepoints SD7, SD10, SD14. The 8-group
# PERMDISP in run_sseries_stats.R does not test the ND-vs-HFD contrast the paper claims.
# Same filtering, rarefaction (seed 42) and Bray-Curtis as run_sseries_stats.R.
suppressPackageStartupMessages(library(vegan))
BASE <- "/home/unnot/vault/kjm/miseq/dada2_output"; OUT <- file.path(BASE, "diversity")
asv <- read.table(file.path(BASE, "asv_table.tsv"), sep="\t", header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="", quote="")
asv$sequence <- NULL; asv_t <- t(as.matrix(asv))
meta <- read.table(file.path(BASE, "metadata.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
meta <- meta[meta$experiment == "ARG", ]
keep <- intersect(rownames(asv_t)[rowSums(asv_t) >= 5000], meta$sample_id)
asv_t <- asv_t[keep, ]; meta <- meta[match(keep, meta$sample_id), ]
GROUPS8 <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
sel <- meta$group %in% GROUPS8; asv_t <- asv_t[sel, ]; meta <- meta[sel, ]
set.seed(42); asv_rare <- rrarefy(asv_t, min(rowSums(asv_t)))
PAIRS <- list(Amoxicillin=c("AMX_CTL","AMX_HFD"), Ciprofloxacin=c("CIP_CTL","CIP_HFD"),
              Erythromycin=c("ERY_CTL","ERY_HFD"), None=c("CTL","HFD"))
res <- data.frame()
for (tp in c("SD7","SD10","SD14")) for (ab in names(PAIRS)) {
  idx <- which(meta$timepoint == tp & meta$group %in% PAIRS[[ab]])
  g <- factor(meta$group[idx], levels=PAIRS[[ab]])
  if (min(table(g)) < 2) next
  sub <- asv_rare[idx, ]; sub <- sub[, colSums(sub) > 0]
  d <- vegdist(sub, method="bray")
  set.seed(42); pm <- adonis2(d ~ g, permutations=999)
  bd <- betadisper(d, g); set.seed(42); pd <- permutest(bd, permutations=999)
  disp <- tapply(bd$distances, g, mean)
  res <- rbind(res, data.frame(timepoint=tp, antibiotic=ab,
    n_ND=sum(g==PAIRS[[ab]][1]), n_HFD=sum(g==PAIRS[[ab]][2]),
    PERMANOVA_R2=round(pm$R2[1],3), PERMANOVA_p=pm$`Pr(>F)`[1],
    disp_ND=round(disp[1],3), disp_HFD=round(disp[2],3),
    PERMDISP_F=round(pd$tab$F[1],3), PERMDISP_p=pd$tab$`Pr(>F)`[1]))
}
print(res, row.names=FALSE)
write.table(res, file.path(OUT,"pairwise_diet_permanova_permdisp.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
