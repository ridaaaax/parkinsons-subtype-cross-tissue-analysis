.libPaths(c("/home/codespace/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

# 1. Load the verified brain-subtype gene lists
brain_sets <- read.csv(
  "data/brain_signature/brain_gene_sets.csv",
  stringsAsFactors = FALSE
)

gene_sets <- split(brain_sets$gene, brain_sets$signature)
gene_sets <- lapply(gene_sets, unique)

cat("\\nBRAIN SIGNATURE SIZES\\n")
print(sapply(gene_sets, length))

# 2. Prepare ranked peripheral-tissue statistics
prepare_stats <- function(file) {
  x <- read.csv(file, check.names = FALSE)

  gene_col <- intersect(
    c("Geneid", "Gene", "gene_id", "gene"),
    names(x)
  )[1]

  stat_col <- intersect(
    c("t", "statistic", "t_stat"),
    names(x)
  )[1]

  if (is.na(gene_col) || is.na(stat_col)) {
    stop("Could not identify gene ID or t-statistic columns in: ", file)
  }

  d <- data.frame(
    gene = sub("\\.[0-9]+$", "", as.character(x[[gene_col]])),
    stat = as.numeric(x[[stat_col]])
  )

  d <- d[!is.na(d$gene) & d$gene != "" & is.finite(d$stat), ]
  d <- d[order(abs(d$stat), decreasing = TRUE), ]
  d <- d[!duplicated(d$gene), ]

  stats <- setNames(d$stat, d$gene)
  sort(stats, decreasing = TRUE)
}

blood_stats <- prepare_stats(
  "results/GSE165082_blood_DE_TMM_limma_voom.csv"
)

muscle_stats <- prepare_stats(
  "results/GSE128177_muscle_DE_limma.csv"
)

tissues <- list(
  Blood = blood_stats,
  Muscle = muscle_stats
)

# 3. Run the six signature-enrichment tests
all_results <- list()
leading_results <- list()
k <- 1

for (tissue_name in names(tissues)) {

  stats <- tissues[[tissue_name]]

  matched_sets <- lapply(
    gene_sets,
    function(genes) intersect(genes, names(stats))
  )

  cat("\n", tissue_name, "MATCHING\n", sep = "")
  print(sapply(matched_sets, length))

  usable_sets <- matched_sets[
    lengths(matched_sets) >= 15
  ]

  if (length(usable_sets) == 0) {
    stop("No gene sets had enough matched genes in ", tissue_name)
  }

  fg <- fgseaMultilevel(
    pathways = usable_sets,
    stats = stats,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

  for (i in seq_len(nrow(fg))) {

    signature <- fg$pathway[i]
    edge <- fg$leadingEdge[[i]]

    all_results[[k]] <- data.frame(
      tissue = tissue_name,
      signature = signature,
      n_signature = length(gene_sets[[signature]]),
      n_matched = length(matched_sets[[signature]]),
      size = fg$size[i],
      ES = fg$ES[i],
      NES = fg$NES[i],
      pval = fg$pval[i],
      leading_edge = paste(edge, collapse = ";"),
      stringsAsFactors = FALSE
    )

    if (length(edge) > 0) {
      leading_results[[k]] <- data.frame(
        tissue = tissue_name,
        signature = signature,
        gene = edge
      )
    }

    k <- k + 1
  }
}

# 4. Combine results and adjust p-values across all six tests
results <- do.call(rbind, all_results)
results$padj_global <- p.adjust(results$pval, method = "BH")

results <- results[order(results$padj_global), ]

dir.create("results", showWarnings = FALSE)

write.csv(
  results,
  "results/brain_signature_cross_tissue_enrichment.csv",
  row.names = FALSE
)

leading <- do.call(rbind, leading_results)

write.csv(
  leading,
  "results/brain_signature_leading_edge_genes.csv",
  row.names = FALSE
)

# 5. Generate the comparison figure
plot_data <- results
plot_data$signature <- factor(
  plot_data$signature,
  levels = c("PDC1_unique", "PDC2_unique", "Shared_brain_DEGs")
)

p <- ggplot(
  plot_data,
  aes(x = signature, y = NES, fill = tissue)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Enrichment of Brain-Derived Parkinson's Signatures",
    subtitle = "Ranked peripheral-tissue differential-expression statistics",
    x = "Brain-derived gene signature",
    y = "Normalized Enrichment Score (NES)",
    fill = "Tissue"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "results/brain_signature_cross_tissue_enrichment.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

# 6. Display the results
cat("\nFINAL ENRICHMENT RESULTS\n")

print(
  results[, c(
    "tissue", "signature", "n_signature",
    "n_matched", "NES", "pval", "padj_global"
  )],
  row.names = FALSE
)

cat("\nFILES CREATED:\n")
cat("results/brain_signature_cross_tissue_enrichment.csv\n")
cat("results/brain_signature_leading_edge_genes.csv\n")
cat("results/brain_signature_cross_tissue_enrichment.png\n")
