import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from pathlib import Path


# --------------------------------------------------
# 1. Load the expression matrix
# --------------------------------------------------

data_file = Path("data/GSE128177_featureCounts.txt.gz")

print("Loading GSE128177 expression matrix...")

df = pd.read_csv(
    data_file,
    sep="\t",
    compression="gzip"
)

print("Data loaded successfully.")


# --------------------------------------------------
# 2. Define sample groups
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


# --------------------------------------------------
# 3. Select the primary comparison samples
# --------------------------------------------------

primary_samples = [
    sample
    for sample in df.columns[1:]
    if sample in groups
]

expression = df.set_index("Geneid")[primary_samples]

pd_samples = [
    sample
    for sample in primary_samples
    if groups[sample] == "Parkinson_disease"
]

control_samples = [
    sample
    for sample in primary_samples
    if groups[sample] == "Old_adult"
]

print("\n--- SAMPLE COUNTS ---")
print("Parkinson's disease samples:", len(pd_samples))
print("Old-adult control samples:", len(control_samples))


# --------------------------------------------------
# 4. Convert to CPM
# --------------------------------------------------

print("\nCalculating CPM...")

library_sizes = expression.sum(axis=0)

cpm = expression.div(library_sizes, axis=1) * 1_000_000


# --------------------------------------------------
# 5. Filter low-expression genes
# --------------------------------------------------

print("Filtering low-expression genes...")

keep_genes = (cpm > 1).sum(axis=1) >= 6

filtered_cpm = cpm.loc[keep_genes]

print("Genes before filtering:", expression.shape[0])
print("Genes after filtering:", filtered_cpm.shape[0])


# --------------------------------------------------
# 6. Log-transform the expression values
# --------------------------------------------------

log_cpm = np.log2(filtered_cpm + 1)


# --------------------------------------------------
# 7. Perform Welch's t-test for every gene
# --------------------------------------------------

print("\nPerforming gene-by-gene statistical tests...")

results = []

for gene in log_cpm.index:

    pd_values = log_cpm.loc[gene, pd_samples]
    control_values = log_cpm.loc[gene, control_samples]

    mean_pd = pd_values.mean()
    mean_control = control_values.mean()

    log2_fold_change = mean_pd - mean_control

    statistic, p_value = ttest_ind(
        pd_values,
        control_values,
        equal_var=False
    )

    results.append({
        "gene_id": gene,
        "mean_PD": mean_pd,
        "mean_old_control": mean_control,
        "log2_fold_change": log2_fold_change,
        "p_value": p_value
    })


results_df = pd.DataFrame(results)


# --------------------------------------------------
# 8. Benjamini-Hochberg FDR correction
# --------------------------------------------------

print("Applying Benjamini-Hochberg FDR correction...")


def benjamini_hochberg(p_values):

    p_values = np.asarray(p_values)
    number_of_tests = len(p_values)

    order = np.argsort(p_values)
    sorted_p_values = p_values[order]

    adjusted_values = (
        sorted_p_values
        * number_of_tests
        / np.arange(1, number_of_tests + 1)
    )

    adjusted_values = np.minimum.accumulate(
        adjusted_values[::-1]
    )[::-1]

    adjusted_values = np.minimum(
        adjusted_values,
        1.0
    )

    corrected = np.empty(number_of_tests)
    corrected[order] = adjusted_values

    return corrected


results_df["FDR"] = benjamini_hochberg(
    results_df["p_value"].values
)


# --------------------------------------------------
# 9. Add significance labels
# --------------------------------------------------

results_df["significant"] = (
    (results_df["FDR"] < 0.05)
    &
    (results_df["log2_fold_change"].abs() >= 1)
)


results_df = results_df.sort_values(
    by=["FDR", "p_value"]
)


# --------------------------------------------------
# 10. Save the results
# --------------------------------------------------

Path("results").mkdir(exist_ok=True)

results_df.to_csv(
    "results/GSE128177_muscle_differential_expression.csv",
    index=False
)

significant_results = results_df[
    results_df["significant"]
]

significant_results.to_csv(
    "results/GSE128177_muscle_significant_genes.csv",
    index=False
)


# --------------------------------------------------
# 11. Print a summary
# --------------------------------------------------

print("\n--- DIFFERENTIAL-EXPRESSION SUMMARY ---")

print("Total tested genes:", len(results_df))

print(
    "Genes with FDR < 0.05 and absolute log2 fold-change >= 1:",
    len(significant_results)
)

print("\nTop 20 results:")

print(
    results_df.head(20).to_string(index=False)
)


# --------------------------------------------------
# 12. Create a basic volcano plot
# --------------------------------------------------

plot_df = results_df.copy()

plot_df["minus_log10_FDR"] = -np.log10(
    plot_df["FDR"].clip(lower=1e-300)
)

plt.figure(figsize=(9, 7))

plt.scatter(
    plot_df["log2_fold_change"],
    plot_df["minus_log10_FDR"],
    s=8,
    alpha=0.5
)

significant_plot_df = plot_df[
    plot_df["significant"]
]

plt.scatter(
    significant_plot_df["log2_fold_change"],
    significant_plot_df["minus_log10_FDR"],
    s=12,
    label="FDR < 0.05 and |log2FC| >= 1"
)

plt.axvline(
    1,
    linestyle="--"
)

plt.axvline(
    -1,
    linestyle="--"
)

plt.axhline(
    -np.log10(0.05),
    linestyle="--"
)

plt.xlabel("Log2 fold-change: Parkinson's disease vs old controls")

plt.ylabel("-Log10(FDR)")

plt.title(
    "GSE128177 skeletal muscle differential expression"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "results/GSE128177_muscle_volcano_plot.png",
    dpi=300
)

plt.show()


print("\nDifferential-expression analysis complete.")

print(
    "Full results saved to:",
    "results/GSE128177_muscle_differential_expression.csv"
)

print(
    "Significant genes saved to:",
    "results/GSE128177_muscle_significant_genes.csv"
)

print(
    "Volcano plot saved to:",
    "results/GSE128177_muscle_volcano_plot.png"
)