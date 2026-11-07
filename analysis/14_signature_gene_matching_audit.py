
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# 14. Brain signature vs blood gene matching audit
# ============================================================

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

count_file = Path("data/GSE165082_PD-CC.counts.txt.gz")
signature_file = Path("results/brain_signature_42_genes.csv")
limma_file = Path("results/GSE165082_blood_DE_limma_voom.csv")

output_dir = Path("results")
output_dir.mkdir(exist_ok=True)

# ------------------------------------------------------------
# 2. Load datasets
# ------------------------------------------------------------

print("Loading original blood counts...")

counts = pd.read_csv(
    count_file,
    sep="\t",
    compression="gzip",
    index_col=0
)

print("Loading brain signature...")
signature = pd.read_csv(signature_file)

print("Loading blood limma-voom results...")
limma = pd.read_csv(limma_file)

print("All files loaded.")

# ------------------------------------------------------------
# 3. Validate required columns
# ------------------------------------------------------------

if "Gene" not in signature.columns:
    raise ValueError("Brain signature is missing Gene column.")

if "Gene" not in limma.columns:
    raise ValueError("Blood limma results are missing Gene column.")

if signature["Gene"].isna().any():
    raise ValueError("Brain signature contains missing gene IDs.")

# ------------------------------------------------------------
# 4. Prepare gene identifiers
# ------------------------------------------------------------

def clean_gene_ids(series):
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
    )

signature["Gene_base"] = clean_gene_ids(signature["Gene"])

limma["Gene_base"] = clean_gene_ids(limma["Gene"])

counts.index = clean_gene_ids(pd.Series(counts.index)).values

# ------------------------------------------------------------
# 5. Validate and prepare raw counts
# ------------------------------------------------------------

counts = counts.apply(pd.to_numeric, errors="coerce")

if counts.isna().any().any():
    raise ValueError(
        "Non-numeric or missing values found in blood counts."
    )

if (counts < 0).any().any():
    raise ValueError("Negative values found in raw counts.")

if not np.allclose(
    counts.to_numpy(),
    np.round(counts.to_numpy())
):
    raise ValueError("Raw count matrix contains non-integer values.")

library_sizes = counts.sum(axis=0)

if (library_sizes <= 0).any():
    raise ValueError("One or more samples have zero library size.")

# ------------------------------------------------------------
# 6. Reproduce original CPM filtering
# ------------------------------------------------------------

cpm = counts.div(library_sizes, axis=1) * 1_000_000

# Same filtering rule as script 13:
# CPM > 1 in at least 6 samples

passes_filter = cpm.gt(1).sum(axis=1) >= 6

raw_gene_ids = set(counts.index)

filtered_gene_ids = set(
    counts.index[passes_filter]
)

limma_gene_ids = set(limma["Gene_base"])

# ------------------------------------------------------------
# 7. Audit all 42 brain signature genes
# ------------------------------------------------------------

signature = signature.drop_duplicates("Gene_base").copy()

audit = signature.copy()

audit["present_in_raw_counts"] = (
    audit["Gene_base"].isin(raw_gene_ids)
)

audit["passes_CPM_filter"] = (
    audit["Gene_base"].isin(filtered_gene_ids)
)

audit["present_in_limma_results"] = (
    audit["Gene_base"].isin(limma_gene_ids)
)

# Number of samples with CPM > 1 for each matched gene

cpm_sample_counts = cpm.gt(1).sum(axis=1)

cpm_lookup = (
    pd.DataFrame({
        "Gene_base": counts.index,
        "samples_CPM_gt_1": cpm_sample_counts.values,
        "maximum_CPM": cpm.max(axis=1).values
    })
    .groupby("Gene_base")
    .agg({
        "samples_CPM_gt_1": "max",
        "maximum_CPM": "max"
    })
)

audit["samples_CPM_gt_1"] = (
    audit["Gene_base"]
    .map(cpm_lookup["samples_CPM_gt_1"])
)

audit["maximum_CPM"] = (
    audit["Gene_base"]
    .map(cpm_lookup["maximum_CPM"])
)

# ------------------------------------------------------------
# 8. Classify matching status
# ------------------------------------------------------------

def classify_gene(row):

    if not row["present_in_raw_counts"]:
        return "Not found in raw counts by exact Ensembl ID"

    if not row["passes_CPM_filter"]:
        return "Present in raw counts but removed by CPM filter"

    if not row["present_in_limma_results"]:
        return "Passed CPM filter but absent from limma results"

    return "Present in raw counts, passed filter and tested"

audit["matching_status"] = audit.apply(
    classify_gene,
    axis=1
)

# ------------------------------------------------------------
# 9. Print audit summary
# ------------------------------------------------------------

print("\n--- GENE MATCHING AUDIT ---")

print("Brain signature genes:", len(audit))

print(
    "Present in raw blood counts:",
    audit["present_in_raw_counts"].sum()
)

print(
    "Passed CPM filter:",
    audit["passes_CPM_filter"].sum()
)

print(
    "Present in limma results:",
    audit["present_in_limma_results"].sum()
)

print("\n--- MATCHING STATUS ---")

print(
    audit["matching_status"].value_counts().to_string()
)

print("\n--- SAMPLE INFORMATION ---")
print("Blood genes before filtering:", len(raw_gene_ids))
print("Blood genes after filtering:", len(filtered_gene_ids))
print("Genes tested by limma:", len(limma_gene_ids))
print("Blood samples:", counts.shape[1])

# ------------------------------------------------------------
# 10. Save detailed audit
# ------------------------------------------------------------

audit.to_csv(
    output_dir / "42_gene_blood_matching_audit.csv",
    index=False
)

# ------------------------------------------------------------
# 11. Save audit summary
# ------------------------------------------------------------

summary = pd.DataFrame({
    "metric": [
        "brain_signature_genes",
        "blood_genes_before_filtering",
        "blood_genes_after_filtering",
        "blood_genes_tested_by_limma",
        "signature_genes_present_in_raw_counts",
        "signature_genes_passing_CPM_filter",
        "signature_genes_present_in_limma",
        "signature_genes_removed_by_CPM_filter",
        "signature_genes_not_found_in_raw_counts"
    ],
    "value": [
        len(audit),
        len(raw_gene_ids),
        len(filtered_gene_ids),
        len(limma_gene_ids),
        audit["present_in_raw_counts"].sum(),
        audit["passes_CPM_filter"].sum(),
        audit["present_in_limma_results"].sum(),
        (
            audit["present_in_raw_counts"]
            & ~audit["passes_CPM_filter"]
        ).sum(),
        (~audit["present_in_raw_counts"]).sum()
    ]
})

summary.to_csv(
    output_dir / "42_gene_blood_matching_audit_summary.csv",
    index=False
)

print("\nAudit files saved successfully.")
print("Analysis 14 complete.")
