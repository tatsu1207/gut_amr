#!/usr/bin/env Rscript
# Factorial PERMANOVA (antibiotic x diet) at each timepoint, and per-timepoint PCoA for Fig. 2C (2026-11-21).
# Same samples and rarefaction as run_diversity.R: experiment ARG, 8 groups, >= 5,000 reads, rrarefy to the
# minimum depth (5,159) with seed 42, Bray-Curtis, 999 permutations. T6h (untreated groups only) is skipped.
suppressMessages(library(vegan))
D <- "/home/unnot/vault/kjm/miseq/dada2_output/"
asv <- read.delim(paste0(D, "asv_table.tsv"), row.names = 1, check.names = FALSE)[, -1]
meta <- read.delim(paste0(D, "metadata.tsv"), stringsAsFactors = FALSE)
G <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
m <- meta[meta$experiment == "ARG" & meta$group %in% G, ]
x <- t(asv[, intersect(m$sample_id, colnames(asv))]); x <- x[, colSums(x) > 0]
x <- x[rowSums(x) >= 5000, ]
m <- m[match(rownames(x), m$sample_id), ]
set.seed(42); xr <- rrarefy(x, min(rowSums(x)))
cat("samples", nrow(xr), "depth", min(rowSums(x)), "\n")
m$abx <- ifelse(grepl("_", m$group), sub("_.*", "", m$group), "None")
m$abx <- factor(m$abx, levels = c("None", "AMX", "CIP", "ERY"))
m$diet2 <- factor(ifelse(grepl("HFD$", m$group), "HFD", "ND"), levels = c("ND", "HFD"))
res <- list(); pts <- list()
for (tp in c("T0h","T24h","T72h","S72h","SD7","SD10","SD14")) {
  i <- m$timepoint == tp; d <- vegdist(xr[i, ], "bray"); mm <- m[i, ]
  set.seed(42); a <- adonis2(d ~ abx * diet2, data = mm, permutations = 999, by = "terms")
  set.seed(42); bda <- permutest(betadisper(d, mm$abx), permutations = 999)$tab$`Pr(>F)`[1]
  set.seed(42); bdd <- permutest(betadisper(d, mm$diet2), permutations = 999)$tab$`Pr(>F)`[1]
  set.seed(42); g8 <- adonis2(d ~ group, data = mm, permutations = 999)
  set.seed(42); bdg <- permutest(betadisper(d, mm$group), permutations = 999)$tab$`Pr(>F)`[1]
  res[[tp]] <- data.frame(timepoint = tp, n = sum(i),
    R2_antibiotic = a$R2[1], p_antibiotic = a$`Pr(>F)`[1], R2_diet = a$R2[2], p_diet = a$`Pr(>F)`[2],
    R2_interaction = a$R2[3], p_interaction = a$`Pr(>F)`[3], R2_residual = a$R2[4],
    PERMDISP_p_antibiotic = bda, PERMDISP_p_diet = bdd, R2_group8 = g8$R2[1], p_group8 = g8$`Pr(>F)`[1], PERMDISP_p_group8 = bdg)
  if (tp %in% c("T0h", "T72h", "SD14")) {
    pc <- cmdscale(d, k = 2, eig = TRUE); ev <- pc$eig / sum(pc$eig[pc$eig > 0]) * 100
    pts[[tp]] <- data.frame(sample_id = rownames(pc$points), timepoint = tp, PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2],
                            group = mm$group, pct1 = ev[1], pct2 = ev[2])
  }
}
r <- do.call(rbind, res); print(format(r, digits = 3))
write.table(r, paste0(D, "diversity/permanova_factorial_by_timepoint.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, pts), paste0(D, "diversity/pcoa_per_timepoint.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
