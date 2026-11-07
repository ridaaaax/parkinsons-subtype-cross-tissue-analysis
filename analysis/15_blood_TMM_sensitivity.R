# ============================================================
# GSE165082 Blood TMM Sensitivity Analysis
# edgeR TMM + limma-voom
# Original CPM filtering and PD vs Control design preserved
# ============================================================

.libPaths(c("~/R/library", .libPaths()))

library(edgeR)
library(limma)

INPUT_FILE <- "data/GSE165082_PD-CC.counts.txt.gz"
ORIGINAL_FILE <- "results/GSE165082_blood_DE_limma_voom.csv"
RESULTS_DIR <- "results"

# ------------------------------------------------------------
# 1. Read and validate counts
# ------------------------------------------------------------

counts <- read.delim(
  INPUT_FILE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

stopifnot(
  !anyNA(counts),
  all(counts >= 0),
  all(counts == floor(counts)),
  all(colSums(counts) > 0)
)

cat("Genes:", nrow(counts), "\n")
cat("Samples:", ncol(counts), "\n")

# ------------------------------------------------------------
# 2. Define groups
# ------------------------------------------------------------

sample_names <- colnames(counts)

group <- ifelse(
  grepl("_PD$", sample_names), "PD",
  ifelse(grepl("_CC$", sample_names), "Control", NA)
)

stopifnot(!anyNA(group))

group <- factor(group, levels = c("Control", "PD"))

cat("\nSample groups:\n")
print(table(group))

# ------------------------------------------------------------
# 3. Reproduce original CPM filtering
# ------------------------------------------------------------

library_sizes <- colSums(counts)

cpm_matrix <- sweep(
  counts,
  2,
  library_sizes,
  "/"
) * 1000000

keep <- rowSums(cpm_matrix > 1) >= 6

filtered_counts <- counts[keep, , drop = FALSE]

cat("\nGenes retained:", nrow(filtered_counts), "\n")

# ------------------------------------------------------------
# 4. TMM normalisation
# ------------------------------------------------------------

dge <- DGEList(
  counts = filtered_counts,
  group = group
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

cat("\nTMM normalisation factors:\n")
print(dge$samples[, c("lib.size", "norm.factors")])

# ------------------------------------------------------------
# 5. Design matrix
# ------------------------------------------------------------

design <- model.matrix(~ group)

colnames(design) <- c("Intercept", "PD_vs_Control")

# ------------------------------------------------------------
# 6. limma-voom
# ------------------------------------------------------------

v <- voom(
  dge,
  design,
  plot = FALSE
)

fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# ------------------------------------------------------------
# 7. Extract complete results
# ------------------------------------------------------------

tmm_results <- topTable(
  fit,
  coef = "PD_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

tmm_results$Gene <- rownames(tmm_results)

tmm_results <- tmm_results[
  ,
  c("Gene", "logFC", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B")
]

write.csv(
  tmm_results,
  file.path(RESULTS_DIR, "GSE165082_blood_DE_TMM_limma_voom.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Compare against original analysis
# ------------------------------------------------------------

original <- read.csv(
  ORIGINAL_FILE,
  stringsAsFactors = FALSE
)

comparison <- merge(
  original[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val")],
  tmm_results[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val")],
  by = "Gene",
  suffixes = c("_original", "_TMM")
)

write.csv(
  comparison,
  file.path(RESULTS_DIR, "GSE165082_blood_TMM_comparison.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------

summary <- data.frame(
  metric = c(
    "genes_tested",
    "samples",
    "PD_samples",
    "control_samples",
    "positive_t_original",
    "negative_t_original",
    "positive_t_TMM",
    "negative_t_TMM",
    "median_t_original",
    "median_t_TMM",
    "correlation_original_vs_TMM_t",
    "significant_genes_original",
    "significant_genes_TMM"
  ),
  value = c(
    nrow(tmm_results),
    ncol(counts),
    sum(group == "PD"),
    sum(group == "Control"),
    sum(original$t > 0, na.rm = TRUE),
    sum(original$t < 0, na.rm = TRUE),
    sum(tmm_results$t > 0, na.rm = TRUE),
    sum(tmm_results$t < 0, na.rm = TRUE),
    median(original$t, na.rm = TRUE),
    median(tmm_results$t, na.rm = TRUE),
    cor(comparison$t_original, comparison$t_TMM),
    sum(original$adj.P.Val < 0.05, na.rm = TRUE),
    sum(tmm_results$adj.P.Val < 0.05, na.rm = TRUE)
  )
)

write.csv(
  summary,
  file.path(RESULTS_DIR, "GSE165082_blood_TMM_sensitivity_summary.csv"),
  row.names = FALSE
)

cat("\nTMM sensitivity analysis completed.\n")
print(summary)
