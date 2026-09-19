import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# --------------------------------------------------
# 1. LOAD DATA
# --------------------------------------------------

data_file = Path("data/GSE128177_featureCounts.txt.gz")

print("Loading GSE128177 muscle expression data...")

df = pd.read_csv(
    data_file,
    sep="\t",
    compression="gzip"
)

print("Data loaded successfully.")


# --------------------------------------------------
# 2. SAMPLE GROUPS
# --------------------------------------------------

groups = {
    "M-M-01": "Parkinson_disease",
    "M-M-04": "Parkinson_disease",
    "M-M-08": "Parkinson_disease",
    "M-M-12": "Parkinson_disease",
    "M-M-16": "Parkinson_disease",
    "M-M-20": "Parkinson_disease",
    "M-M-24": "Parkinson_disease",
    "M-M-27": "Parkinson_disease",
    "M-M-30": "Parkinson_disease",
    "M-M-33": "Parkinson_disease",
    "M-M-36": "Parkinson_disease",
    "M-M-39": "Parkinson_disease",

    "M-M-02": "Old_adult",
    "M-M-06": "Old_adult",
    "M-M-10": "Old_adult",
    "M-M-14": "Old_adult",
    "M-M-18": "Old_adult",
    "M-M-22": "Old_adult",
    "M-M-25": "Old_adult",
    "M-M-28": "Old_adult",
    "M-M-31": "Old_adult",
    "M-M-34": "Old_adult",
    "M-M-37": "Old_adult",
    "M-M-40": "Old_adult"
}

primary_samples = [
    sample
    for sample in df.columns[1:]
    if sample in groups
]

print("\nNumber of primary samples:", len(primary_samples))

print(
    "Parkinson's disease samples:",
    sum(groups[s] == "Parkinson_disease" for s in primary_samples)
)

print(
    "Older control samples:",
    sum(groups[s] == "Old_adult" for s in primary_samples)
)


# --------------------------------------------------
# 3. CREATE EXPRESSION MATRIX
# --------------------------------------------------

df["Gene_clean"] = (
    df["Geneid"]
    .astype(str)
    .str.split(".")
    .str[0]
)

expression = df.set_index("Gene_clean")[primary_samples]

# Remove duplicate cleaned gene IDs if any
expression = expression.groupby(expression.index).first()


# --------------------------------------------------
# 4. CONVERT TO CPM
# --------------------------------------------------

library_sizes = expression.sum(axis=0)

cpm = (
    expression
    .div(library_sizes, axis=1)
    * 1_000_000
)


# --------------------------------------------------
# 5. FILTER LOW-EXPRESSION GENES
# --------------------------------------------------

keep_genes = (
    (cpm > 1).sum(axis=1) >= 6
)

filtered_cpm = cpm.loc[keep_genes]

print("\n--- EXPRESSION FILTERING ---")
print("Genes before filtering:", expression.shape[0])
print("Genes after filtering:", filtered_cpm.shape[0])


# --------------------------------------------------
# 6. LOG2 TRANSFORMATION
# --------------------------------------------------

log_cpm = np.log2(filtered_cpm + 1)


# --------------------------------------------------
# 7. CREATE DESIGN MATRIX
# --------------------------------------------------

group_vector = np.array([
    1 if groups[sample] == "Parkinson_disease"
    else 0
    for sample in primary_samples
])

design = sm.add_constant(group_vector)

print("\nDesign matrix:")
print("Rows:", design.shape[0])
print("Columns:", design.shape[1])


# --------------------------------------------------
# 8. GENE-BY-GENE LINEAR MODEL
# --------------------------------------------------

results = []

print("\nRunning linear models...")

for gene in log_cpm.index:

    y = log_cpm.loc[gene].values

    model = sm.OLS(y, design).fit()

    # coefficient for Parkinson's disease
    effect = model.params[1]

    p_value = model.pvalues[1]

    results.append({
        "gene_id": gene,
        "log2FC_PD_vs_old": effect,
        "p_value": p_value
    })


results = pd.DataFrame(results)


# --------------------------------------------------
# 9. MULTIPLE-TEST CORRECTION
# --------------------------------------------------

results["FDR"] = multipletests(
    results["p_value"],
    method="fdr_bh"
)[1]


# --------------------------------------------------
# 10. SIGNIFICANT GENES
# --------------------------------------------------

results["significant"] = (
    (results["FDR"] < 0.05)
    &
    (results["log2FC_PD_vs_old"].abs() >= 1)
)


# --------------------------------------------------
# 11. SORT RESULTS
# --------------------------------------------------

results = results.sort_values(
    by="p_value"
)


# --------------------------------------------------
# 12. PRINT SUMMARY
# --------------------------------------------------

significant = results[
    results["significant"]
]

print("\n--- DIFFERENTIAL EXPRESSION SUMMARY ---")

print(
    "Genes tested:",
    len(results)
)

print(
    "Significant genes:",
    len(significant)
)

print(
    "Upregulated in PD:",
    (
        significant["log2FC_PD_vs_old"] > 0
    ).sum()
)

print(
    "Downregulated in PD:",
    (
        significant["log2FC_PD_vs_old"] < 0
    ).sum()
)


print("\n--- TOP 20 RESULTS ---")

print(
    results.head(20).to_string(index=False)
)


# --------------------------------------------------
# 13. SAVE RESULTS
# --------------------------------------------------

Path("results").mkdir(exist_ok=True)

results.to_csv(
    "results/GSE128177_muscle_DE_final.csv",
    index=False
)

significant.to_csv(
    "results/GSE128177_muscle_DE_final_significant.csv",
    index=False
)


print("\n--- COMPLETE ---")

print(
    "Full results saved to:"
    " results/GSE128177_muscle_DE_final.csv"
)

print(
    "Significant genes saved to:"
    " results/GSE128177_muscle_DE_final_significant.csv"
)