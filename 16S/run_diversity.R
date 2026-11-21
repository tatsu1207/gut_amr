#!/usr/bin/env Rscript
# Alpha and Beta diversity analysis for KJM obesity+antibiotics study
# Uses only ARG experiment samples (excludes other experiments)

library(vegan)

cat("Loading ASV table...\n")
asv <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/asv_table.tsv", row.names=1, check.names=FALSE)
seqs <- asv[, 1, drop=FALSE]  # sequence column
asv <- asv[, -1]  # remove sequence column, keep only counts

cat("Loading metadata...\n")
meta <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/metadata.tsv", stringsAsFactors=FALSE)

# Filter to ARG experiment only (main 8 groups, exclude Early and D0)
main_groups <- c("CTL", "HFD", "AMX_CTL", "AMX_HFD", "CIP_CTL", "CIP_HFD", "ERY_CTL", "ERY_HFD")
meta_arg <- meta[meta$experiment == "ARG" & meta$group %in% main_groups, ]
cat("ARG main group samples:", nrow(meta_arg), "\n")

# Subset ASV table to these samples
samples_keep <- intersect(meta_arg$sample_id, colnames(asv))
cat("Samples in ASV table:", length(samples_keep), "\n")
asv_arg <- asv[, samples_keep]

# Remove empty ASVs
asv_arg <- asv_arg[rowSums(asv_arg) > 0, ]
cat("ASVs after filtering:", nrow(asv_arg), "\n")

# Transpose for vegan (samples as rows)
asv_t <- t(asv_arg)

# Per-sample read counts
read_counts <- rowSums(asv_t)
cat("\nRead count summary:\n")
print(summary(read_counts))
cat("Samples < 5000 reads:", sum(read_counts < 5000), "\n")
cat("Samples < 10000 reads:", sum(read_counts < 10000), "\n")

# Filter samples with < 5000 reads
min_reads <- 5000
keep <- read_counts >= min_reads
cat("\nFiltering samples with <", min_reads, "reads:", sum(!keep), "removed\n")
asv_t <- asv_t[keep, ]
meta_filt <- meta_arg[meta_arg$sample_id %in% rownames(asv_t), ]
meta_filt <- meta_filt[match(rownames(asv_t), meta_filt$sample_id), ]

cat("Final samples:", nrow(asv_t), "\n")
cat("Final ASVs:", ncol(asv_t), "\n")

# ── Alpha Diversity ──────────────────────────────────────────────────────────
cat("\n=== Alpha Diversity ===\n")

# Rarefy to minimum depth for alpha diversity
raredepth <- min(rowSums(asv_t))
cat("Rarefaction depth:", raredepth, "\n")

set.seed(42)
asv_rare <- rrarefy(asv_t, raredepth)

# Calculate indices
chao1 <- estimateR(asv_rare)["S.chao1", ]
shannon <- diversity(asv_rare, index="shannon")
observed <- specnumber(asv_rare)
simpson <- diversity(asv_rare, index="simpson")

alpha_df <- data.frame(
  sample_id = rownames(asv_rare),
  observed_ASVs = observed,
  chao1 = chao1,
  shannon = shannon,
  simpson = simpson,
  stringsAsFactors = FALSE
)
alpha_df <- merge(alpha_df, meta_filt, by="sample_id")

outdir <- "/home/unnot/vault/kjm/miseq/dada2_output/diversity"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

write.csv(alpha_df, file.path(outdir, "alpha_diversity.csv"), row.names=FALSE)
cat("Alpha diversity saved.\n")

# Summary by group and timepoint
cat("\nAlpha diversity summary (Shannon) by group × timepoint:\n")
agg <- aggregate(shannon ~ group + timepoint, data=alpha_df, FUN=function(x) round(mean(x), 2))
agg <- agg[order(agg$group, agg$timepoint), ]
print(agg)

# ── Beta Diversity ───────────────────────────────────────────────────────────
cat("\n=== Beta Diversity ===\n")

# Bray-Curtis on rarefied data
bc_dist <- vegdist(asv_rare, method="bray")

# PCoA
pcoa_res <- cmdscale(bc_dist, k=3, eig=TRUE)
eig_pct <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  PCoA3 = pcoa_res$points[, 3],
  stringsAsFactors = FALSE
)
pcoa_df <- merge(pcoa_df, meta_filt, by="sample_id")

write.csv(pcoa_df, file.path(outdir, "pcoa_bray_curtis.csv"), row.names=FALSE)
cat("PCoA axes: PC1=", eig_pct[1], "%, PC2=", eig_pct[2], "%, PC3=", eig_pct[3], "%\n")

# Save eigenvalues
write.csv(data.frame(axis=1:length(eig_pct), variance_pct=eig_pct),
          file.path(outdir, "pcoa_eigenvalues.csv"), row.names=FALSE)

# PCoA without R7 (SD7) for Fig. 2C (2026-11-21): R7 was sequenced alone on a failed run (KDFCK);
# same rarefied table and Bray-Curtis distances, R7 samples dropped before ordination.
keep7 <- meta_filt$sample_id[meta_filt$timepoint != "SD7"]
keep7 <- intersect(labels(bc_dist), keep7)
bc7 <- as.dist(as.matrix(bc_dist)[keep7, keep7])
p7 <- cmdscale(bc7, k=3, eig=TRUE)
e7 <- round(p7$eig / sum(p7$eig[p7$eig > 0]) * 100, 1)
p7df <- merge(data.frame(sample_id=rownames(p7$points), PCoA1=p7$points[,1], PCoA2=p7$points[,2],
                         PCoA3=p7$points[,3], stringsAsFactors=FALSE), meta_filt, by="sample_id")
write.csv(p7df, file.path(outdir, "pcoa_bray_curtis_noR7.csv"), row.names=FALSE)
write.csv(data.frame(axis=1:length(e7), variance_pct=e7), file.path(outdir, "pcoa_eigenvalues_noR7.csv"), row.names=FALSE)
cat("PCoA without R7:", length(keep7), "samples; PC1=", e7[1], "%, PC2=", e7[2], "%\n")

# PERMANOVA by group at each timepoint
cat("\nPERMANOVA by group at each timepoint:\n")
timepoints <- unique(meta_filt$timepoint)
permanova_results <- list()
for (tp in timepoints) {
  idx <- meta_filt$timepoint == tp
  if (sum(idx) < 10) next
  tp_dist <- vegdist(asv_rare[idx, ], method="bray")
  tp_meta <- meta_filt[idx, ]
  perm <- adonis2(tp_dist ~ group, data=tp_meta, permutations=999)
  r2 <- round(perm$R2[1], 3)
  pval <- perm$`Pr(>F)`[1]
  cat(sprintf("  %s: R2=%.3f, p=%.4f (n=%d)\n", tp, r2, pval, sum(idx)))
  permanova_results[[tp]] <- data.frame(timepoint=tp, R2=r2, p_value=pval, n=sum(idx))
}
permanova_df <- do.call(rbind, permanova_results)
write.csv(permanova_df, file.path(outdir, "permanova_by_timepoint.csv"), row.names=FALSE)

# PERMANOVA: diet effect within antibiotic groups
cat("\nPERMANOVA: diet effect (ND vs HFD) within each antibiotic at SD14:\n")
antibiotics <- c("Amoxicillin", "Ciprofloxacin", "Erythromycin")
sd14_idx <- meta_filt$timepoint == "SD14"
for (ab in antibiotics) {
  idx <- sd14_idx & meta_filt$antibiotic == ab
  if (sum(idx) < 4) next
  ab_dist <- vegdist(asv_rare[idx, ], method="bray")
  ab_meta <- meta_filt[idx, ]
  perm <- adonis2(ab_dist ~ diet, data=ab_meta, permutations=999)
  cat(sprintf("  %s at SD14: R2=%.3f, p=%.4f (n=%d)\n", ab, perm$R2[1], perm$`Pr(>F)`[1], sum(idx)))
}

# Bray-Curtis dissimilarity from baseline (T0h) per mouse
cat("\nCalculating dissimilarity from baseline...\n")
bc_matrix <- as.matrix(bc_dist)
baseline_dissim <- data.frame()
for (i in 1:nrow(meta_filt)) {
  s <- meta_filt$sample_id[i]
  grp <- meta_filt$group[i]
  mid <- meta_filt$mouse_id[i]
  tp <- meta_filt$timepoint[i]
  if (tp == "T0h") next

  # Find this mouse's T0h sample
  t0_idx <- meta_filt$group == grp & meta_filt$mouse_id == mid & meta_filt$timepoint == "T0h"
  if (sum(t0_idx) != 1) next
  t0_sample <- meta_filt$sample_id[t0_idx]

  if (s %in% rownames(bc_matrix) & t0_sample %in% rownames(bc_matrix)) {
    dissim <- bc_matrix[s, t0_sample]
    baseline_dissim <- rbind(baseline_dissim, data.frame(
      sample_id=s, group=grp, diet=meta_filt$diet[i],
      antibiotic=meta_filt$antibiotic[i], timepoint=tp,
      day=meta_filt$day[i], mouse_id=mid, bc_from_baseline=dissim
    ))
  }
}
write.csv(baseline_dissim, file.path(outdir, "bc_from_baseline.csv"), row.names=FALSE)
cat("Baseline dissimilarity calculated for", nrow(baseline_dissim), "sample pairs.\n")

cat("\n=== Done ===\n")
cat("Output directory:", outdir, "\n")
