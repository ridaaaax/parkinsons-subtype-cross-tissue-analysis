# ============================================================
# GSE165082 Blood Differential Expression
# limma-voom: Parkinson's disease vs Control
# ============================================================

library(limma)

INPUT_FILE <- "data/GSE165082_PD-CC.counts.txt.gz"
RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")

dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Read raw count matrix
# ------------------------------------------------------------

cat("Reading blood count matrix...\n")

counts <- read.delim(
  INPUT_FILE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

cat("Genes:", nrow(counts), "\n")
cat("Samples:", ncol(counts), "\n")


# ------------------------------------------------------------
# 2. Define sample groups
# ------------------------------------------------------------

sample_names <- colnames(counts)

group <- ifelse(
  grepl("_PD$", sample_names),
  "PD",
  ifelse(grepl("_CC$", sample_names), "Control", NA)
)

group <- factor(
  group,
  levels = c("Control", "PD")
)

cat("\nSample groups:\n")
print(table(group))


# ------------------------------------------------------------
# 3. Calculate CPM
# ------------------------------------------------------------

library_sizes <- colSums(counts)

cpm_matrix <- sweep(
  counts,
  2,
  library_sizes,
  "/"
) * 1000000


# ------------------------------------------------------------
# 4. Low-expression filtering
# ------------------------------------------------------------

# Keep genes with CPM > 1 in at least 6 samples

keep <- rowSums(cpm_matrix > 1) >= 6

filtered_counts <- counts[keep, ]

cat("\nGenes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(filtered_counts), "\n")
cat("Genes removed:", sum(!keep), "\n")


# ------------------------------------------------------------
# 5. Design matrix
# ------------------------------------------------------------

design <- model.matrix(~ group)

colnames(design) <- c(
  "Intercept",
  "PD_vs_Control"
)

cat("\nDesign matrix:\n")
print(design)


# ------------------------------------------------------------
# 6. limma-voom transformation
# ------------------------------------------------------------

cat("\nRunning limma-voom...\n")

v <- voom(
  filtered_counts,
  design = design,
  plot = TRUE
)

# Save the voom mean-variance plot
dev.copy(
  png,
  file.path(
    FIGURES_DIR,
    "GSE165082_blood_voom_mean_variance.png"
  ),
  width = 1000,
  height = 800,
  res = 150
)

dev.off()


# ------------------------------------------------------------
# 7. Fit linear model
# ------------------------------------------------------------

fit <- lmFit(
  v,
  design
)


# ------------------------------------------------------------
# 8. Empirical Bayes moderation
# ------------------------------------------------------------

fit <- eBayes(
  fit,
  robust = TRUE
)


# ------------------------------------------------------------
# 9. Extract PD vs Control results
# ------------------------------------------------------------

results <- topTable(
  fit,
  coef = "PD_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$Gene <- rownames(results)

results <- results[
  ,
  c(
    "Gene",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


# ------------------------------------------------------------
# 10. Save complete DE results
# ------------------------------------------------------------

full_output <- file.path(
  RESULTS_DIR,
  "GSE165082_blood_DE_limma_voom.csv"
)

write.csv(
  results,
  full_output,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Define significant genes
# ------------------------------------------------------------

significant <- subset(
  results,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

sig_output <- file.path(
  RESULTS_DIR,
  "GSE165082_blood_DE_limma_voom_significant.csv"
)

write.csv(
  significant,
  sig_output,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

up <- sum(
  significant$logFC > 0
)

down <- sum(
  significant$logFC < 0
)

cat("\nDifferential-expression summary:\n")
cat("Genes tested:", nrow(results), "\n")
cat("Significant genes:", nrow(significant), "\n")
cat("Up in PD:", up, "\n")
cat("Down in PD:", down, "\n")


if (nrow(results) > 0) {

  cat(
    "Smallest raw P value:",
    min(results$P.Value),
    "\n"
  )

  cat(
    "Smallest adjusted P value:",
    min(results$adj.P.Val),
    "\n"
  )
}


# ------------------------------------------------------------
# 13. Save summary table
# ------------------------------------------------------------

summary <- data.frame(
  metric = c(
    "genes_before_filtering",
    "genes_after_filtering",
    "genes_tested",
    "PD_samples",
    "Control_samples",
    "significant_genes_FDR_0.05_abs_logFC_1",
    "up_in_PD",
    "down_in_PD"
  ),
  value = c(
    nrow(counts),
    nrow(filtered_counts),
    nrow(results),
    sum(group == "PD"),
    sum(group == "Control"),
    nrow(significant),
    up,
    down
  )
)

write.csv(
  summary,
  file.path(
    RESULTS_DIR,
    "GSE165082_blood_DE_limma_voom_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Volcano plot
# ------------------------------------------------------------

results$significance <- "Not significant"

results$significance[
  results$adj.P.Val < 0.05 &
    results$logFC >= 1
] <- "Up"

results$significance[
  results$adj.P.Val < 0.05 &
    results$logFC <= -1
] <- "Down"

plot(
  results$logFC,
  -log10(results$P.Value),
  pch = 20,
  cex = 0.6,
  xlab = "log2 fold change (PD vs Control)",
  ylab = "-log10(P value)",
  main = "GSE165082 blood differential expression"
)

abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)

dev.copy(
  png,
  file.path(
    FIGURES_DIR,
    "GSE165082_blood_volcano.png"
  ),
  width = 1000,
  height = 800,
  res = 150
)

dev.off()


cat("\nAnalysis complete.\n")

cat("Created:\n")
cat(full_output, "\n")
cat(sig_output, "\n")
cat(
  file.path(
    RESULTS_DIR,
    "GSE165082_blood_DE_limma_voom_summary.csv"
  ),
  "\n"
)