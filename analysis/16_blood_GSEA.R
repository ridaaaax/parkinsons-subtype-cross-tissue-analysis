# ============================================================
# 16. Blood transcriptome pathway analysis using fgsea
# Parkinson's disease cross-tissue project
#
# Primary model: TMM-normalised limma-voom
# Sensitivity model: non-TMM limma-voom
# Ranking: signed moderated t-statistic
# Pathways: Hallmark, Reactome, GO Biological Process
# ============================================================

# 1. Package library
.libPaths(c("/home/codespace/R/library", .libPaths()))

required_packages <- c("fgsea", "msigdbr", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)

cat("\nPackages loaded successfully.\n")

# 2. File paths
tmm_file <- "results/GSE165082_blood_DE_TMM_limma_voom.csv"
original_file <- "results/GSE165082_blood_DE_limma_voom.csv"

output_dir <- "results/blood_GSEA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tmm_file) || !file.exists(original_file)) {
  stop("One or both blood differential-expression files are missing.")
}

# 3. Retrieve MSigDB pathway collections
cat("\nRetrieving pathway databases...\n")

hallmark <- msigdbr::msigdbr(
  species = "human",
  collection = "H"
)

reactome <- msigdbr::msigdbr(
  species = "human",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

gobp <- msigdbr::msigdbr(
  species = "human",
  collection = "C5",
  subcollection = "GO:BP"
)

make_pathways <- function(database) {

  database <- database[
    !is.na(database$ensembl_gene) &
    database$ensembl_gene != "" &
    !is.na(database$gs_name),
  ]

  pathways <- split(
    database$ensembl_gene,
    database$gs_name
  )

  lapply(pathways, unique)
}

pathway_collections <- list(
  Hallmark = make_pathways(hallmark),
  Reactome = make_pathways(reactome),
  GO_BP = make_pathways(gobp)
)

cat("\nPathway collection sizes:\n")
for (nm in names(pathway_collections)) {
  cat(nm, ":", length(pathway_collections[[nm]]), "pathways\n")
}

# 4. Prepare ranked gene lists
prepare_ranks <- function(file) {

  dat <- read.csv(file, stringsAsFactors = FALSE)

  required_columns <- c("Gene", "t")

  if (!all(required_columns %in% names(dat))) {
    stop("Required columns Gene and t are missing in ", file)
  }

  dat <- dat[
    !is.na(dat$Gene) &
    dat$Gene != "" &
    is.finite(dat$t),
  ]

  # Remove Ensembl version suffixes if present
  dat$Gene <- sub("\\.[0-9]+$", "", dat$Gene)

  # Resolve duplicate Ensembl identifiers by retaining
  # the entry with the largest absolute moderated t-statistic
  dat <- dat[order(abs(dat$t), decreasing = TRUE), ]
  dat <- dat[!duplicated(dat$Gene), ]

  stats <- dat$t
  names(stats) <- dat$Gene

  stats <- sort(stats, decreasing = TRUE)

  if (length(stats) < 100) {
    stop("Unexpectedly few genes available for GSEA in ", file)
  }

  cat("Genes ranked:", length(stats), "\n")

  stats
}

# 5. Run GSEA for one model and one collection
run_gsea <- function(stats, pathways, model_name, collection_name) {

  cat("\nRunning:", model_name, "-", collection_name, "\n")

result <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000
)

  result <- as.data.frame(result)

  result <- result[order(result$padj, result$pval), ]

  result$model <- model_name
  result$collection <- collection_name

  # Make leading-edge genes suitable for CSV export
  result$leadingEdge <- vapply(
    result$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )

  result
}

# 6. Save results and plots
save_gsea_results <- function(result, model_name, collection_name) {

  prefix <- paste(model_name, collection_name, sep = "_")

  output_file <- file.path(
    output_dir,
    paste0(prefix, "_GSEA.csv")
  )

  write.csv(result, output_file, row.names = FALSE)

  significant <- result[
    !is.na(result$padj) & result$padj < 0.05,
  ]

  write.csv(
    significant,
    file.path(output_dir, paste0(prefix, "_significant.csv")),
    row.names = FALSE
  )

  # Plot top 20 pathways by adjusted p-value
  plot_data <- result[
    !is.na(result$padj) & is.finite(result$NES),
  ]

  plot_data <- head(plot_data, 20)

  if (nrow(plot_data) > 0) {

    plot_data$pathway <- factor(
      plot_data$pathway,
      levels = rev(plot_data$pathway)
    )

    p <- ggplot(
      plot_data,
      aes(x = NES, y = pathway)
    ) +
      geom_col() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(
        title = paste(model_name, collection_name, "GSEA"),
        subtitle = paste(
          "Top 20 pathways by adjusted p-value;",
          sum(result$padj < 0.05, na.rm = TRUE),
          "pathways with FDR < 0.05"
        ),
        x = "Normalized Enrichment Score (NES)",
        y = "Pathway"
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold")
      )

    ggsave(
      filename = file.path(
        output_dir,
        paste0(prefix, "_top_pathways.png")
      ),
      plot = p,
      width = 11,
      height = 8,
      dpi = 300
    )
  }

  cat(
    "Completed:", prefix,
    "| Significant pathways:", nrow(significant), "\n"
  )
}

# 7. Run both blood models
models <- list(
  TMM = tmm_file,
  nonTMM = original_file
)

all_results <- list()

for (model_name in names(models)) {

  cat("\n====================================\n")
  cat("MODEL:", model_name, "\n")
  cat("====================================\n")

  stats <- prepare_ranks(models[[model_name]])

  for (collection_name in names(pathway_collections)) {

    result <- run_gsea(
      stats = stats,
      pathways = pathway_collections[[collection_name]],
      model_name = model_name,
      collection_name = collection_name
    )

    save_gsea_results(
      result,
      model_name,
      collection_name
    )

    all_results[[paste(model_name, collection_name, sep = "_")]] <-
      result
  }
}

# 8. Combined summary
combined <- do.call(rbind, all_results)

write.csv(
  combined,
  file.path(output_dir, "ALL_blood_GSEA_results.csv"),
  row.names = FALSE
)

summary_table <- aggregate(
  padj ~ model + collection,
  data = combined,
  FUN = function(x) sum(x < 0.05, na.rm = TRUE)
)

names(summary_table)[3] <- "significant_pathways_FDR_0.05"

write.csv(
  summary_table,
  file.path(output_dir, "GSEA_summary.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("GSEA ANALYSIS COMPLETED\n")
cat("====================================\n")

print(summary_table)

cat("\nResults saved in:", output_dir, "\n")
