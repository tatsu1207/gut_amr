#!/usr/bin/env Rscript
# ALDEx2 diet comparisons: HFD vs ND within each antibiotic at key timepoints

library(ALDEx2)

outdir <- "/home/unnot/vault/kjm/miseq/dada2_output/differential_abundance"
dir.create(outdir, showWarnings = FALSE)

asv <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/asv_table.tsv",
                  row.names = 1, check.names = FALSE)
asv$sequence <- NULL

tax <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/taxonomy.tsv",
                  row.names = 1, check.names = FALSE)

meta <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/metadata.tsv",
                   row.names = 1, check.names = FALSE)

main_groups <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
# main experiment, the eight study groups only (Early cohorts, ARG_D0 batch and excluded duplicates out; 2026-11-20)
meta_arg <- meta[meta$experiment == "ARG" & meta$group %in% main_groups, ]
shared <- intersect(colnames(asv), rownames(meta_arg))
asv_arg <- asv[, shared]
meta_arg <- meta_arg[shared, ]

# Filter ASVs
min_samples <- ceiling(ncol(asv_arg) * 0.05)
keep <- rowSums(asv_arg > 5) >= min_samples
asv_filt <- asv_arg[keep, ]
cat(sprintf("Retained %d ASVs\n", sum(keep)))

# HFD vs ND group pairs
diet_pairs <- list(
  Amoxicillin = c(nd = "AMX_CTL", hfd = "AMX_HFD"),
  Ciprofloxacin = c(nd = "CIP_CTL", hfd = "CIP_HFD"),
  Erythromycin = c(nd = "ERY_CTL", hfd = "ERY_HFD"),
  None = c(nd = "CTL", hfd = "HFD")
)

timepoints <- c("T0h", "T72h", "S72h", "SD7", "SD14")

for (abx_name in names(diet_pairs)) {
  nd_grp <- diet_pairs[[abx_name]]["nd"]
  hfd_grp <- diet_pairs[[abx_name]]["hfd"]

  for (tp in timepoints) {
    nd_mask <- meta_arg$group == nd_grp & meta_arg$timepoint == tp
    hfd_mask <- meta_arg$group == hfd_grp & meta_arg$timepoint == tp
    nd_mask[is.na(nd_mask)] <- FALSE
    hfd_mask[is.na(hfd_mask)] <- FALSE

    if (sum(nd_mask) < 3 || sum(hfd_mask) < 3) {
      cat(sprintf("Skipping %s @ %s: ND=%d, HFD=%d\n", abx_name, tp, sum(nd_mask), sum(hfd_mask)))
      next
    }

    samples_use <- rownames(meta_arg)[nd_mask | hfd_mask]
    samples_use <- intersect(samples_use, colnames(asv_filt))
    conds <- ifelse(meta_arg[samples_use, "diet"] == "HFD", "HFD", "ND")

    cat(sprintf("\nALDEx2: %s @ %s — ND vs HFD (n=%d vs %d)\n",
                abx_name, tp, sum(conds == "ND"), sum(conds == "HFD")))

    counts <- asv_filt[, samples_use]
    counts <- counts[rowSums(counts) > 0, ]

    tryCatch({
      aldex_out <- aldex(counts, conds, mc.samples = 128, test = "t",
                         effect = TRUE, verbose = FALSE)
      aldex_out$ASV_ID <- rownames(aldex_out)
      aldex_out <- merge(aldex_out, tax, by.x = "ASV_ID", by.y = "row.names")
      aldex_out$antibiotic <- abx_name
      aldex_out$comparison <- paste0("ND_vs_HFD_", tp)

      fname <- sprintf("%s/aldex2_%s_diet_%s.tsv", outdir, abx_name, tp)
      write.table(aldex_out, fname, sep = "\t", row.names = FALSE, quote = FALSE)

      sig <- aldex_out[aldex_out$wi.eBH < 0.05, ]
      cat(sprintf("  Significant ASVs (BH < 0.05): %d\n", nrow(sig)))
      if (nrow(sig) > 0) {
        sig_sorted <- sig[order(sig$wi.eBH), ]
        top_n <- min(5, nrow(sig_sorted))
        for (i in 1:top_n) {
          cat(sprintf("    %s (%s): effect=%.2f, BH=%.4f\n",
                      sig_sorted$Genus[i], sig_sorted$ASV_ID[i],
                      sig_sorted$effect[i], sig_sorted$wi.eBH[i]))
        }
      }
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
    })
  }
}

cat("\nDiet comparisons complete.\n")
