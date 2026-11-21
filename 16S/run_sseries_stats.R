#!/usr/bin/env Rscript
# Additional statistics requested in review (S1 repeated measures, S2 PERMDISP).
#   S1: linear mixed model on alpha diversity with mouse as a random effect
#   S2: PERMDISP (betadisper) alongside PERMANOVA, to show whether a significant
#       R2 reflects a location shift or merely heterogeneous dispersion
suppressPackageStartupMessages({library(vegan); library(lme4); library(lmerTest)})

BASE <- "/home/unnot/vault/kjm/miseq/dada2_output"
OUT  <- file.path(BASE, "diversity")

asv <- read.table(file.path(BASE, "asv_table.tsv"), sep="\t", header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="", quote="")
asv$sequence <- NULL
asv_t <- t(as.matrix(asv))
meta <- read.table(file.path(BASE, "metadata.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
meta <- meta[meta$experiment == "ARG", ]
# CSL = the same untreated control animals, relabelled at recovery timepoints
meta$prefix <- sub("^ARG_([A-Za-z]+)[0-9]+_.*$", "\\1", meta$sample_id)
meta$mouse_uid <- paste0(meta$group, "_", meta$mouse_id)

keep <- rownames(asv_t)[rowSums(asv_t) >= 5000]
keep <- intersect(keep, meta$sample_id)
asv_t <- asv_t[keep, ]
meta <- meta[match(keep, meta$sample_id), ]
GROUPS8 <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
sel <- meta$group %in% GROUPS8
asv_t <- asv_t[sel, ]; meta <- meta[sel, ]
cat("samples:", nrow(asv_t), " mice:", length(unique(meta$mouse_uid)),
    " groups:", length(unique(meta$group)), "\n")

set.seed(42)
asv_rare <- rrarefy(asv_t, min(rowSums(asv_t)))
meta$shannon  <- diversity(asv_rare, index="shannon")
meta$observed <- specnumber(asv_rare)

TP <- c("T0h","T6h","T24h","T72h","S72h","SD7","SD10","SD14")
meta$timepoint <- factor(meta$timepoint, levels=TP)
meta$group <- factor(meta$group)
meta$diet  <- factor(meta$diet)

# ---- S1: mixed model, mouse as random effect ----
cat("\n=== S1: linear mixed model, Shannon ~ group * timepoint + (1|mouse), R7 excluded ===\n")
# R7 (SD7) excluded (2026-11-21): sequenced alone on run KDFCK, which failed (33% of reads passed
# filtering; untreated mice lost two thirds of their ASVs at R7 only). Rarefaction above still uses the
# full sample set, so every other result is unchanged.
mm <- droplevels(meta[meta$timepoint != "SD7", ])
cat("mixed model samples (R7 excluded):", nrow(mm), "\n")
m <- lmer(shannon ~ group * timepoint + (1|mouse_uid), data=mm, REML=TRUE)
an <- anova(m)
print(an)
write.table(data.frame(term=rownames(an), an), file.path(OUT,"mixed_model_shannon_anova.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
vc <- as.data.frame(VarCorr(m))
cat("\nvariance components:\n"); print(vc[, c("grp","vcov","sdcor")])
icc <- vc$vcov[vc$grp=="mouse_uid"] / sum(vc$vcov)
cat(sprintf("intraclass correlation (mouse): %.3f\n", icc))
write.table(data.frame(component=vc$grp, variance=vc$vcov, sd=vc$sdcor, ICC=c(icc, NA)[seq_len(nrow(vc))]),
            file.path(OUT,"mixed_model_shannon_varcomp.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ---- S2: PERMANOVA + PERMDISP per timepoint ----
cat("\n=== S2: PERMANOVA with PERMDISP per timepoint ===\n")
res <- data.frame()
for (tp in TP) {
  idx <- which(meta$timepoint == tp)
  if (length(idx) < 10) next
  sub <- asv_rare[idx, ]; sub <- sub[, colSums(sub) > 0]
  g <- droplevels(meta$group[idx])
  if (nlevels(g) < 2) next
  d <- vegdist(sub, method="bray")
  set.seed(42)
  pm <- adonis2(d ~ g, permutations=999)
  bd <- betadisper(d, g)
  set.seed(42)
  pd <- permutest(bd, permutations=999)
  res <- rbind(res, data.frame(
    timepoint = tp, n = length(idx), n_groups = nlevels(g),
    PERMANOVA_R2 = round(pm$R2[1], 4), PERMANOVA_F = round(pm$F[1], 3),
    PERMANOVA_p  = pm$`Pr(>F)`[1],
    PERMDISP_F   = round(pd$tab$F[1], 3), PERMDISP_p = pd$tab$`Pr(>F)`[1]))
}
print(res)
write.table(res, file.path(OUT,"permanova_with_permdisp.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\nwritten to", OUT, "\n")
