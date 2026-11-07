
# 18. Corrected mitochondrial leading-edge integration
# Blood: TMM differential-expression model
# Muscle: limma differential-expression model

blood_gsea_file <- "results/blood_GSEA/TMM_GO_BP_GSEA.csv"
muscle_gsea_file <- "results/muscle_GSEA/muscle_GO_BP_GSEA.csv"

blood_de_file <- "results/GSE165082_blood_DE_TMM_limma_voom.csv"
muscle_de_file <- "results/GSE128177_muscle_DE_limma.csv"

output_file <- "results/mitochondrial_leading_edge_TMM_blood_muscle.csv"

pathways <- c(
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
  "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_ELECTRON_TRANSPORT_CHAIN",
  "GOBP_ATP_BIOSYNTHETIC_PROCESS",
  "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE"
)

# Check inputs
files <- c(
  blood_gsea_file, muscle_gsea_file,
  blood_de_file, muscle_de_file
)

if (!all(file.exists(files))) {
  stop("One or more required input files are missing.")
}

if (file.exists(output_file)) {
  stop("Output already exists. Refusing to overwrite: ", output_file)
}

# Read data
blood_gsea <- read.csv(blood_gsea_file, stringsAsFactors = FALSE)
muscle_gsea <- read.csv(muscle_gsea_file, stringsAsFactors = FALSE)

blood_de <- read.csv(blood_de_file, stringsAsFactors = FALSE)
muscle_de <- read.csv(muscle_de_file, stringsAsFactors = FALSE)

# Check required columns
stopifnot(
  all(c("pathway", "NES", "padj", "leadingEdge") %in% names(blood_gsea)),
  all(c("pathway", "NES", "padj", "leadingEdge") %in% names(muscle_gsea)),
  all(c("Gene", "logFC", "t", "P.Value", "adj.P.Val") %in% names(blood_de)),
  all(c("Geneid", "logFC", "t", "P.Value", "adj.P.Val") %in% names(muscle_de))
)

# Select the six pathways
blood_selected <- blood_gsea[blood_gsea$pathway %in% pathways, ]
muscle_selected <- muscle_gsea[muscle_gsea$pathway %in% pathways, ]

if (nrow(blood_selected) != 6 || nrow(muscle_selected) != 6) {
  stop("Expected all six pathways in both GSEA files.")
}

# Parse leading-edge genes and calculate pathway frequency
extract_frequency <- function(dat, tissue) {
  genes <- unlist(
    lapply(dat$leadingEdge, function(x) {
      if (is.na(x) || x == "") character(0)
      else strsplit(x, ";", fixed = TRUE)[[1]]
    }),
    use.names = FALSE
  )

  frequency <- sort(table(genes), decreasing = TRUE)

  data.frame(
    Gene = names(frequency),
    pathway_frequency = as.integer(frequency),
    tissue = tissue,
    stringsAsFactors = FALSE
  )
}

blood_freq <- extract_frequency(blood_selected, "blood")
muscle_freq <- extract_frequency(muscle_selected, "muscle")

# Normalize gene identifiers
blood_de$Gene <- sub("\\.[0-9]+$", "", blood_de$Gene)
muscle_de$Geneid <- sub("\\.[0-9]+$", "", muscle_de$Geneid)

# Match DE statistics
blood_idx <- match(blood_freq$Gene, blood_de$Gene)
muscle_idx <- match(muscle_freq$Gene, muscle_de$Gene)

blood_result <- data.frame(
  Gene = blood_freq$Gene,
  blood_pathway_frequency = blood_freq$pathway_frequency,
  blood_logFC = blood_de$logFC[blood_idx],
  blood_t = blood_de$t[blood_idx],
  blood_P = blood_de$P.Value[blood_idx],
  blood_FDR = blood_de$adj.P.Val[blood_idx]
)

muscle_result <- data.frame(
  Gene = muscle_freq$Gene,
  muscle_pathway_frequency = muscle_freq$pathway_frequency,
  muscle_logFC = muscle_de$logFC[muscle_idx],
  muscle_t = muscle_de$t[muscle_idx],
  muscle_P = muscle_de$P.Value[muscle_idx],
  muscle_FDR = muscle_de$adj.P.Val[muscle_idx]
)

# Merge by Ensembl identifier
integrated <- merge(
  blood_result,
  muscle_result,
  by = "Gene",
  all = TRUE
)

# Compare directions only where both values are available
integrated$logFC_direction <- ifelse(
  is.na(integrated$blood_logFC) | is.na(integrated$muscle_logFC),
  NA,
  ifelse(
    sign(integrated$blood_logFC) == sign(integrated$muscle_logFC),
    "same",
    "opposite"
  )
)

# Save results
write.csv(integrated, output_file, row.names = FALSE)

# Report audit summary
cat("\nCorrected mitochondrial leading-edge integration\n")
cat("Blood pathways:", nrow(blood_selected), "\n")
cat("Muscle pathways:", nrow(muscle_selected), "\n")
cat("Unique blood leading-edge genes:", nrow(blood_freq), "\n")
cat("Unique muscle leading-edge genes:", nrow(muscle_freq), "\n")
cat("Genes in integrated table:", nrow(integrated), "\n")
cat("Matched blood genes:", sum(!is.na(integrated$blood_t)), "\n")
cat("Matched muscle genes:", sum(!is.na(integrated$muscle_t)), "\n")
cat("\nDirection comparison:\n")
print(table(integrated$logFC_direction, useNA = "ifany"))

cat("\nOutput saved to:", output_file, "\n")
