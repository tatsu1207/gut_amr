#!/usr/bin/env Rscript
# Pre-diet baseline (D0 = diet start, before HFD feeding): were ND- and HFD-assigned mice alike before
# the diet, and how far did 13 weeks of diet move them (D0 vs T0h = antibiotic baseline)?
# D0 samples: experiment == "ARG_D0", groups CTL/HFD = diet assignment (Early-cohort "EAR" samples excluded).
# Same processing as run_diversity.R: rrarefy to 5,159 reads (seed 42), Bray-Curtis, adonis2 999 perm.
suppressMessages(library(vegan))
set.seed(42)
D <- "/home/unnot/vault/kjm/miseq/dada2_output/"
asv  <- read.delim(paste0(D, "asv_table.tsv"), row.names = 1, check.names = FALSE)[, -1]
meta <- read.delim(paste0(D, "metadata.tsv"), stringsAsFactors = FALSE)
G <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
d0 <- meta[meta$experiment == "ARG_D0" & meta$group %in% c("CTL","HFD"), ]
t0 <- meta[meta$experiment == "ARG" & meta$group %in% G & meta$timepoint == "T0h", ]
d0$stage <- "D0_prediet"; t0$stage <- "T0_postdiet"
m <- rbind(d0, t0)
x <- t(asv[, m$sample_id]); m <- m[rowSums(x) >= 5000, ]; x <- x[m$sample_id, ]
x <- rrarefy(x, 5159); x <- x[, colSums(x) > 0]
m$shannon <- diversity(x, "shannon"); m$richness <- rowSums(x > 0)
bc <- vegdist(x, "bray")
out <- file(paste0(D, "diversity/d0_baseline_stats.txt"), "w")
p <- function(...) { cat(..., "\n"); cat(..., "\n", file = out) }
p("n by stage x diet:"); tb <- table(m$stage, m$diet); for (r in rownames(tb)) p(r, paste(colnames(tb), tb[r, ], collapse = "  "))
for (st in c("D0_prediet", "T0_postdiet")) {
  i <- m$stage == st
  a <- adonis2(as.dist(as.matrix(bc)[i, i]) ~ diet, data = m[i, ], permutations = 999)
  bd <- permutest(betadisper(as.dist(as.matrix(bc)[i, i]), m$diet[i]), permutations = 999)
  w1 <- wilcox.test(shannon ~ diet, data = m[i, ], exact = FALSE)
  w2 <- wilcox.test(richness ~ diet, data = m[i, ], exact = FALSE)
  p(sprintf("%s  PERMANOVA diet R2 = %.3f p = %.3f | PERMDISP p = %.3f", st, a$R2[1], a$`Pr(>F)`[1], bd$tab$`Pr(>F)`[1]))
  for (dt in c("ND", "HFD")) p(sprintf("   %s Shannon %.2f +/- %.2f  richness %.0f +/- %.0f", dt,
        mean(m$shannon[i & m$diet == dt]), sd(m$shannon[i & m$diet == dt]),
        mean(m$richness[i & m$diet == dt]), sd(m$richness[i & m$diet == dt])))
  p(sprintf("   Shannon ND vs HFD p = %.3f ; richness p = %.3f", w1$p.value, w2$p.value))
}
a <- adonis2(bc ~ stage * diet, data = m, permutations = 999, by = "terms")
p("Combined D0 + T0, terms:"); for (k in 1:3) p(sprintf("   %s R2 = %.3f p = %.3f", rownames(a)[k], a$R2[k], a$`Pr(>F)`[k]))
# genus-level Lactococcus at D0 vs T0 (% of genus-assigned reads, as elsewhere)
tax <- read.delim(paste0(D, "taxonomy.tsv"), row.names = 1, stringsAsFactors = FALSE)
gcol <- grep("^genus$", colnames(tax), ignore.case = TRUE, value = TRUE)
raw <- t(asv[, m$sample_id]); gen <- tax[colnames(raw), gcol]; ok <- !is.na(gen) & gen != ""
gsum <- rowsum(t(raw[, ok]), gen[ok]); pct <- sweep(gsum, 2, colSums(gsum), "/") * 100
for (g in c("Lactococcus", "Faecalibaculum", "Akkermansia")) for (st in c("D0_prediet", "T0_postdiet")) for (dt in c("ND", "HFD")) {
  v <- pct[g, m$stage == st & m$diet == dt]; p(sprintf("%-15s %-12s %-3s %5.1f +/- %4.1f %%", g, st, dt, mean(v), sd(v))) }
L <- pct["Lactococcus", ]
for (st in c("D0_prediet", "T0_postdiet")) { i <- m$stage == st
  p(sprintf("Lactococcus %s ND vs HFD Wilcoxon p = %.2g", st, wilcox.test(L[i & m$diet == "ND"], L[i & m$diet == "HFD"], exact = FALSE)$p.value)) }
for (dt in c("ND", "HFD"))
  p(sprintf("Lactococcus %s pre-diet vs T0 Wilcoxon p = %.2g", dt, wilcox.test(L[m$stage == "D0_prediet" & m$diet == dt], L[m$stage == "T0_postdiet" & m$diet == dt], exact = FALSE)$p.value))
close(out)
write.table(m[, c("sample_id","stage","diet","group","shannon","richness")], paste0(D, "diversity/d0_baseline_alpha.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
