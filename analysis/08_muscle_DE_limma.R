# 08_muscle_DE_limma.R
# Proper limma analysis of GSE128177 skeletal muscle
# PD vs age/sex-matched older controls

library(limma)

# -----------------------------
# 1. Load expression data
# -----------------------------

data_file <- "data/GSE128177_featureCounts.txt.gz"

expr_df <- read.delim(
  data_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

cat("Expression matrix loaded.\n")
cat("Genes:", nrow(expr_df), "\n")
cat("Samples:", ncol(expr_df), "\n")

# -----------------------------
# 2. Define primary samples
# -----------------------------

groups <- c(
  "M-M-01" = "PD",
  "M-M-02" = "Old_control",
  "M-M-03" = "Young_control",
  "M-M-04" = "PD",
  "M-M-06" = "Old_control",
  "M-M-07" = "Young_control",
  "M-M-08" = "PD",
  "M-M-10" = "Old_control",
  "M-M-11" = "Young_control",
  "M-M-12" = "PD",
  "M-M-14" = "Old_control",
  "M-M-15" = "Young_control",
  "M-M-16" = "PD",
  "M-M-18" = "Old_control",
  "M-M-19" = "Young_control",
  "M-M-20" = "PD",
  "M-M-22" = "Old_control",
  "M-M-23" = "Young_control",
  "M-M-24" = "PD",
  "M-M-25" = "Old_control",
  "M-M-26" = "Young_control",
  "M-M-27" = "PD",
  "M-M-28" = "Old_control",
  "M-M-29" = "Young_control",
  "M-M-30" = "PD",
  "M-M-31" = "Old_control",
  "M-M-32" = "Young_control",
  "M-M-33" = "PD",
  "M-M-34" = "Old_control",
  "M-M-35" = "Young_control",
  "M-M-36" = "PD",
  "M-M-37" = "Old_control",
  "M-M-38" = "Young_control",
  "M-M-39" = "PD",
  "M-M-40" = "Old_control",
  "M-M-41" = "Young_control"
)

# Use only PD and age/sex-matched older controls.
primary_samples <- names(groups)[groups %in% c("PD", "Old_control")]

expr <- expr_df[, primary_samples]

group <- factor(
  groups[primary_samples],
  levels = c("Old_control", "PD")
)

cat("\nPrimary analysis samples:", ncol(expr), "\n")
cat("PD:", sum(group == "PD"), "\n")
cat("Old controls:", sum(group == "Old_control"), "\n")

# -----------------------------
# 3. CPM normalization
# -----------------------------

library_sizes <- colSums(expr)

cpm <- sweep(
  expr,
  2,
  library_sizes,
  "/"
) * 1e6

# -----------------------------
# 4. Expression filtering
# -----------------------------

keep <- rowSums(cpm > 1) >= 6

cpm_filtered <- cpm[keep, ]

cat("\nGenes before filtering:", nrow(expr), "\n")
cat("Genes after filtering:", nrow(cpm_filtered), "\n")

# -----------------------------
# 5. Log2 CPM transformation
# -----------------------------

log_cpm <- log2(cpm_filtered + 1)

# -----------------------------
# 6. Design matrix
# -----------------------------

design <- model.matrix(~ group)

colnames(design) <- c(
  "Intercept",
  "PD_vs_Old"
)

# -----------------------------
# 7. limma model
# -----------------------------

fit <- lmFit(log_cpm, design)

fit <- eBayes(
  fit,
  trend = TRUE
)

# -----------------------------
# 8. Extract PD vs old-control results
# -----------------------------

results <- topTable(
  fit,
  coef = "PD_vs_Old",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$Geneid <- rownames(results)

# Reorder columns
results <- results[, c(
  "Geneid",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)]

# -----------------------------
# 9. Significant genes
# -----------------------------

significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) >= 1,
]

upregulated <- significant[
  significant$logFC > 0,
]

downregulated <- significant[
  significant$logFC < 0,
]

cat("\n--- LIMMA RESULTS ---\n")
cat("Genes tested:", nrow(results), "\n")
cat("Significant genes:", nrow(significant), "\n")
cat("Upregulated in PD:", nrow(upregulated), "\n")
cat("Downregulated in PD:", nrow(downregulated), "\n")

# -----------------------------
# 10. Save results
# -----------------------------

dir.create("results", showWarnings = FALSE)

write.csv(
  results,
  "results/GSE128177_muscle_DE_limma.csv",
  row.names = FALSE
)

write.csv(
  significant,
  "results/GSE128177_muscle_DE_limma_significant.csv",
  row.names = FALSE
)

cat("\nResults saved.\n")

# -----------------------------
# 11. Show top genes
# -----------------------------

cat("\n--- TOP 20 GENES ---\n")
print(
  head(
    results,
    20
  )
)

cat("\nAnalysis complete.\n")