#!/usr/bin/env Rscript
# Differential abundance analysis using ALDEx2 on 16S ASV data
# Compares treatment vs baseline (T0h) for each antibiotic group
# Also compares HFD vs ND within antibiotic treatments

library(ALDEx2)

outdir <- "/home/unnot/vault/kjm/miseq/dada2_output/differential_abundance"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load data
asv <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/asv_table.tsv",
                  row.names = 1, check.names = FALSE)
asv$sequence <- NULL  # drop sequence column

tax <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/taxonomy.tsv",
                  row.names = 1, check.names = FALSE)

meta <- read.delim("/home/unnot/vault/kjm/miseq/dada2_output/metadata.tsv",
                   row.names = 1, check.names = FALSE)

# Filter to ARG experiment
main_groups <- c("CTL","HFD","AMX_CTL","AMX_HFD","CIP_CTL","CIP_HFD","ERY_CTL","ERY_HFD")
# main experiment, the eight study groups only (Early cohorts, ARG_D0 batch and excluded duplicates out; 2026-11-20)
meta_arg <- meta[meta$experiment == "ARG" & meta$group %in% main_groups, ]
shared_samples <- intersect(colnames(asv), rownames(meta_arg))
asv_arg <- asv[, shared_samples]
meta_arg <- meta_arg[shared_samples, ]

# Filter low-abundance ASVs (present in at least 5% of samples with >5 reads)
min_samples <- ceiling(ncol(asv_arg) * 0.05)
keep <- rowSums(asv_arg > 5) >= min_samples
asv_filt <- asv_arg[keep, ]
cat(sprintf("Retained %d / %d ASVs after filtering\n", sum(keep), nrow(asv_arg)))

# ============================================================
# Analysis 1: Treatment timepoints vs Baseline (T0h) per antibiotic
# ============================================================
antibiotics <- unique(meta_arg$antibiotic[!is.na(meta_arg$antibiotic)])
timepoints_vs_baseline <- c("T24h", "T72h", "S72h", "SD7", "SD14")

all_results <- list()

for (abx in antibiotics) {
  for (tp in timepoints_vs_baseline) {
    # Get samples for this comparison
    baseline <- meta_arg$timepoint == "T0h" & meta_arg$antibiotic == abx
    treatment <- meta_arg$timepoint == tp & meta_arg$antibiotic == abx

    if (sum(baseline, na.rm=TRUE) < 3 || sum(treatment, na.rm=TRUE) < 3) next

    samples_use <- rownames(meta_arg)[baseline | treatment]
    samples_use <- intersect(samples_use, colnames(asv_filt))
    if (length(samples_use) < 6) next

    conds <- ifelse(meta_arg[samples_use, "timepoint"] == "T0h", "Baseline", tp)

    cat(sprintf("\nALDEx2: %s — Baseline vs %s (n=%d vs %d)\n",
                abx, tp, sum(conds == "Baseline"), sum(conds == tp)))

    counts <- asv_filt[, samples_use]
    # Remove zero-sum rows
    counts <- counts[rowSums(counts) > 0, ]

    tryCatch({
      aldex_out <- aldex(counts, conds, mc.samples = 128, test = "t",
                         effect = TRUE, verbose = FALSE)
      aldex_out$ASV_ID <- rownames(aldex_out)
      aldex_out <- merge(aldex_out, tax, by.x = "ASV_ID", by.y = "row.names")
      aldex_out$antibiotic <- abx
      aldex_out$comparison <- paste0("T0h_vs_", tp)

      # Save individual result
      fname <- sprintf("%s/aldex2_%s_T0h_vs_%s.tsv", outdir, abx, tp)
      write.table(aldex_out, fname, sep = "\t", row.names = FALSE, quote = FALSE)

      # Summarize significant hits
      sig <- aldex_out[aldex_out$wi.eBH < 0.05, ]
      cat(sprintf("  Significant ASVs (BH < 0.05): %d\n", nrow(sig)))
      if (nrow(sig) > 0) {
        cat("  Top hits:\n")
        sig_sorted <- sig[order(sig$wi.eBH), ]
        top_n <- min(5, nrow(sig_sorted))
        for (i in 1:top_n) {
          cat(sprintf("    %s (%s): effect=%.2f, BH=%.4f\n",
                      sig_sorted$Genus[i], sig_sorted$ASV_ID[i],
                      sig_sorted$effect[i], sig_sorted$wi.eBH[i]))
        }
      }

      all_results[[paste0(abx, "_", tp)]] <- aldex_out
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
    })
  }
}

# ============================================================
# Analysis 2: HFD vs ND within each antibiotic at each timepoint
# ============================================================
cat("\n\n=== Diet comparisons (HFD vs ND) ===\n")

for (abx in antibiotics) {
  abx_hfd <- paste0(abx, "_HFD") # e.g., AMX_HFD
  # Find matching ND group
  abx_nd <- paste0(abx, "_CTL")

  for (tp in c("T0h", "T72h", "S72h", "SD14")) {
    hfd_mask <- meta_arg$group == abx_hfd & meta_arg$timepoint == tp
    nd_mask <- meta_arg$group == abx_nd & meta_arg$timepoint == tp

    hfd_mask[is.na(hfd_mask)] <- FALSE
    nd_mask[is.na(nd_mask)] <- FALSE

    if (sum(hfd_mask) < 3 || sum(nd_mask) < 3) next

    samples_use <- rownames(meta_arg)[hfd_mask | nd_mask]
    samples_use <- intersect(samples_use, colnames(asv_filt))
    if (length(samples_use) < 6) next

    conds <- ifelse(meta_arg[samples_use, "diet"] == "HFD", "HFD", "ND")

    cat(sprintf("\nALDEx2: %s @ %s — ND vs HFD (n=%d vs %d)\n",
                abx, tp, sum(conds == "ND"), sum(conds == "HFD")))

    counts <- asv_filt[, samples_use]
    counts <- counts[rowSums(counts) > 0, ]

    tryCatch({
      aldex_out <- aldex(counts, conds, mc.samples = 128, test = "t",
                         effect = TRUE, verbose = FALSE)
      aldex_out$ASV_ID <- rownames(aldex_out)
      aldex_out <- merge(aldex_out, tax, by.x = "ASV_ID", by.y = "row.names")
      aldex_out$antibiotic <- abx
      aldex_out$comparison <- paste0("ND_vs_HFD_", tp)

      fname <- sprintf("%s/aldex2_%s_diet_%s.tsv", outdir, abx, tp)
      write.table(aldex_out, fname, sep = "\t", row.names = FALSE, quote = FALSE)

      sig <- aldex_out[aldex_out$wi.eBH < 0.05, ]
      cat(sprintf("  Significant ASVs (BH < 0.05): %d\n", nrow(sig)))
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
    })
  }
}

# ============================================================
# Summary table of all significant results
# ============================================================
if (length(all_results) > 0) {
  common <- Reduce(intersect, lapply(all_results, names))  # per-comparison column names differ (rab.win.<group>)
  combined <- do.call(rbind, lapply(all_results, function(x) x[, common]))
  sig_all <- combined[combined$wi.eBH < 0.05, ]
  if (nrow(sig_all) > 0) {
    sig_summary <- sig_all[, c("ASV_ID", "antibiotic", "comparison",
                                "effect", "wi.eBH", "Phylum", "Genus")]
    sig_summary <- sig_summary[order(sig_summary$wi.eBH), ]
    write.table(sig_summary,
                file.path(outdir, "aldex2_all_significant.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("\n\nTotal significant results across all comparisons: %d\n", nrow(sig_summary)))
    cat(sprintf("Saved to: %s/aldex2_all_significant.tsv\n", outdir))
  } else {
    cat("\nNo significant results found at BH < 0.05\n")
  }
}

cat("\nDifferential abundance analysis complete.\n")
