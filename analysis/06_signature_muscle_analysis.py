import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, spearmanr, binomtest
from pathlib import Path
from statsmodels.stats.multitest import multipletests

# --------------------------------------------------
# 1. LOAD BRAIN SIGNATURE
# --------------------------------------------------

signature = pd.read_csv(
    "results/brain_signature_42_genes.csv"
)

signature["Gene_clean"] = (
    signature["Gene"]
    .astype(str)
    .str.split(".")
    .str[0]
)

print("Brain signature loaded.")
print("Number of genes:", len(signature))


# --------------------------------------------------
# 2. LOAD MUSCLE DATA
# --------------------------------------------------

df = pd.read_csv(
    "data/GSE128177_featureCounts.txt.gz",
    sep="\t",
    compression="gzip"
)

df["Gene_clean"] = (
    df["Geneid"]
    .astype(str)
    .str.split(".")
    .str[0]
)


# --------------------------------------------------
# 3. DEFINE PRIMARY MUSCLE SAMPLES
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
    sample for sample in df.columns[1:-1]
    if sample in groups
]

print("\nPrimary samples:", len(primary_samples))


# --------------------------------------------------
# 4. CREATE EXPRESSION MATRIX
# --------------------------------------------------

expression = df.set_index("Gene_clean")[primary_samples]

# In case duplicate cleaned IDs occur
expression = expression.groupby(expression.index).first()

library_sizes = expression.sum(axis=0)

cpm = expression.div(
    library_sizes,
    axis=1
) * 1_000_000

log_cpm = np.log2(cpm + 1)


# --------------------------------------------------
# 5. ANALYSE THE 42 BRAIN SIGNATURE GENES
# --------------------------------------------------

results = []

for _, row in signature.iterrows():

    gene = row["Gene_clean"]

    if gene not in log_cpm.index:
        continue

    values = log_cpm.loc[gene]

    pd_values = values[
        [s for s in primary_samples
         if groups[s] == "Parkinson_disease"]
    ]

    old_values = values[
        [s for s in primary_samples
         if groups[s] == "Old_adult"]
    ]

    # Detectability using the same QC threshold
    detectable_samples = (
        cpm.loc[gene, primary_samples] > 1
    ).sum()

    detectable = detectable_samples >= 6

    # Muscle effect
    muscle_log2fc = (
        pd_values.mean() - old_values.mean()
    )

    t_stat, p_value = ttest_ind(
        pd_values,
        old_values,
        equal_var=False
    )

    results.append({
        "Gene": row["Gene"],
        "HGNC": row["HGNC"],
        "brain_logFC": row["logFC"],
        "muscle_log2FC": muscle_log2fc,
        "muscle_t": t_stat,
        "muscle_p": p_value,
        "muscle_detectable_samples": detectable_samples,
        "muscle_detectable": detectable
    })


results = pd.DataFrame(results)


# --------------------------------------------------
# 6. MULTIPLE-TEST CORRECTION
# --------------------------------------------------

results["muscle_FDR"] = np.nan

valid = results["muscle_p"].notna()

results.loc[valid, "muscle_FDR"] = (
    multipletests(
        results.loc[valid, "muscle_p"],
        method="fdr_bh"
    )[1]
)


# --------------------------------------------------
# 7. DIRECTIONAL CONCORDANCE
# --------------------------------------------------

detectable_results = results[
    results["muscle_detectable"]
].copy()

positive_muscle = (
    detectable_results["muscle_log2FC"] > 0
).sum()

negative_muscle = (
    detectable_results["muscle_log2FC"] < 0
).sum()

n_detectable = len(detectable_results)

if n_detectable > 0:

    sign_test = binomtest(
        positive_muscle,
        n_detectable,
        p=0.5,
        alternative="greater"
    )

    sign_p = sign_test.pvalue

else:
    sign_p = np.nan


# --------------------------------------------------
# 8. CORRELATION BETWEEN BRAIN AND MUSCLE EFFECTS
# --------------------------------------------------

if len(detectable_results) >= 3:

    rho, correlation_p = spearmanr(
        detectable_results["brain_logFC"],
        detectable_results["muscle_log2FC"]
    )

else:
    rho = np.nan
    correlation_p = np.nan


# --------------------------------------------------
# 9. PRINT RESULTS
# --------------------------------------------------

print("\n--- SIGNATURE DETECTABILITY ---")
print(
    "42-gene signature genes detectable in muscle:",
    n_detectable,
    "of",
    len(results)
)

print("\n--- DIRECTION OF MUSCLE EFFECTS ---")
print("Positive muscle log2FC:", positive_muscle)
print("Negative muscle log2FC:", negative_muscle)

print("\n--- SIGN TEST ---")
print("P-value:", sign_p)

print("\n--- BRAIN vs MUSCLE EFFECT CORRELATION ---")
print("Spearman rho:", rho)
print("P-value:", correlation_p)

print("\n--- 42-GENE RESULTS ---")

print(
    results[
        [
            "Gene",
            "HGNC",
            "brain_logFC",
            "muscle_log2FC",
            "muscle_p",
            "muscle_FDR",
            "muscle_detectable_samples",
            "muscle_detectable"
        ]
    ].to_string(index=False)
)


# --------------------------------------------------
# 10. SAVE RESULTS
# --------------------------------------------------

Path("results").mkdir(exist_ok=True)

results.to_csv(
    "results/42_gene_muscle_analysis.csv",
    index=False
)

print("\n--- COMPLETE ---")
print(
    "Results saved to: "
    "results/42_gene_muscle_analysis.csv"
)