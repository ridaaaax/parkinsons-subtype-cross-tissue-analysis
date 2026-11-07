# ============================================================
# 17. Skeletal muscle pathway analysis using fgsea
# Parkinson's disease cross-tissue project
#
# Contrast: PD vs older controls
# Ranking: signed moderated t-statistic
# Pathways: Hallmark, Reactome, GO Biological Process
# ============================================================

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

# 1. Input and output paths
input_file <- "results/GSE128177_muscle_DE_limma.csv"
output_dir <- "results/muscle_GSEA"

if (!file.exists(input_file)) {
  stop("Muscle differential-expression file not found.")
}

if (dir.exists(output_dir) && length(list.files(output_dir)) > 0) {
  stop("Output directory already contains files. Refusing to overwrite.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load muscle differential-expression results
dat <- read.csv(input_file, stringsAsFactors = FALSE)

required_columns <- c("Geneid", "t", "logFC")

if (!all(required_columns %in% names(dat))) {
  stop("Required columns missing: ",
       paste(setdiff(required_columns, names(dat)), collapse = ", "))
}

dat$t <- as.numeric(dat$t)
dat$logFC <- as.numeric(dat$logFC)

dat <- dat[
  !is.na(dat$Geneid) &
  dat$Geneid != "" &
  is.finite(dat$t) &
  is.finite(dat$logFC),
]

# Remove Ensembl version suffixes
dat$Geneid <- sub("\\.[0-9]+$", "", dat$Geneid)

# Resolve duplicate identifiers by retaining the largest
# absolute moderated t-statistic
dat <- dat[order(abs(dat$t), decreasing = TRUE), ]
dat <- dat[!duplicated(dat$Geneid), ]

stats <- dat$t
names(stats) <- dat$Geneid
stats <- sort(stats, decreasing = TRUE)

if (length(stats) < 100) {
  stop("Unexpectedly few genes available for GSEA.")
}

cat("Muscle genes ranked:", length(stats), "\n")
cat("Positive t-statistics:", sum(stats > 0), "\n")
cat("Negative t-statistics:", sum(stats < 0), "\n")

# 3. Retrieve pathway collections
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

# 4. Run GSEA
run_gsea <- function(stats, pathways, collection_name) {

  cat("\nRunning muscle", collection_name, "GSEA...\n")

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

  result$leadingEdge <- vapply(
    result$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )

  result
}

# 5. Save results and plots
summary_rows <- list()

for (collection_name in names(pathway_collections)) {

  result <- run_gsea(
    stats,
    pathway_collections[[collection_name]],
    collection_name
  )

  write.csv(
    result,
    file.path(
      output_dir,
      paste0("muscle_", collection_name, "_GSEA.csv")
    ),
    row.names = FALSE
  )

  significant <- result[
    !is.na(result$padj) & result$padj < 0.05,
  ]

  write.csv(
    significant,
    file.path(
      output_dir,
      paste0("muscle_", collection_name, "_significant.csv")
    ),
    row.names = FALSE
  )

  summary_rows[[collection_name]] <- data.frame(
    collection = collection_name,
    pathways_tested = nrow(result),
    significant_FDR_0.05 = nrow(significant),
    positive_NES = sum(result$NES > 0, na.rm = TRUE),
    negative_NES = sum(result$NES < 0, na.rm = TRUE),
    uncalculated_p_values = sum(is.na(result$pval)),
    uncalculated_log2err = sum(is.na(result$log2err))
  )

  plot_data <- result[
    !is.na(result$padj) &
    is.finite(result$NES),
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
        title = paste("Muscle", collection_name, "GSEA"),
        subtitle = paste(
          nrow(significant),
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
      file.path(
        output_dir,
        paste0("muscle_", collection_name, "_top_pathways.png")
      ),
      p,
      width = 11,
      height = 8,
      dpi = 300
    )
  }

  cat(
    collection_name,
    "| Tested:", nrow(result),
    "| FDR significant:", nrow(significant), "\n"
  )
}

# 6. Save summary
summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  file.path(output_dir, "muscle_GSEA_summary.csv"),
  row.names = FALSE
)

cat("\nMuscle GSEA complete.\n")
print(summary_table)
cat("\nResults saved in:", output_dir, "\n")
