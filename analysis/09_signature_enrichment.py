import pandas as pd
import numpy as np
from pathlib import Path

# -----------------------------------
# 1. Load data
# -----------------------------------

limma_file = Path("results/GSE128177_muscle_DE_limma.csv")
signature_file = Path("results/brain_signature_42_genes.csv")

limma = pd.read_csv(limma_file)
signature = pd.read_csv(signature_file)

print("Files loaded successfully.")

# -----------------------------------
# 2. Remove Ensembl version suffixes
# -----------------------------------

limma["Gene_base"] = (
    limma["Geneid"]
    .astype(str)
    .str.replace(r"\.\d+$", "", regex=True)
)

signature["Gene_base"] = (
    signature["Gene"]
    .astype(str)
    .str.replace(r"\.\d+$", "", regex=True)
)

# -----------------------------------
# 3. Match brain signature to muscle
# -----------------------------------

signature_muscle = limma[
    limma["Gene_base"].isin(signature["Gene_base"])
].copy()

print("\n--- SIGNATURE MATCHING ---")
print("Brain signature genes:", len(signature))
print("Matched in muscle limma results:", len(signature_muscle))

# -----------------------------------
# 4. Direction-neutral signal
# -----------------------------------

signature_t = signature_muscle["t"].values
background_t = limma["t"].values

observed_mean_abs_t = np.mean(np.abs(signature_t))

print("\n--- DIRECTION-NEUTRAL SIGNAL ---")
print(
    "Mean absolute muscle t-statistic:",
    observed_mean_abs_t
)
print(
    "Median absolute muscle t-statistic:",
    np.median(np.abs(signature_t))
)

print(
    "Positive t-statistics:",
    np.sum(signature_t > 0)
)

print(
    "Negative t-statistics:",
    np.sum(signature_t < 0)
)

# -----------------------------------
# 5. Permutation test
# -----------------------------------

rng = np.random.default_rng(42)

random_means = []

for _ in range(10000):

    random_genes = rng.choice(
        background_t,
        size=len(signature_t),
        replace=False
    )

    random_means.append(
        np.mean(np.abs(random_genes))
    )

random_means = np.array(random_means)

p_value = (
    np.sum(random_means >= observed_mean_abs_t) + 1
) / (
    len(random_means) + 1
)

print("\n--- PERMUTATION TEST ---")
print(
    "Observed mean absolute t:",
    observed_mean_abs_t
)

print(
    "Permutation p-value:",
    p_value
)

# -----------------------------------
# 6. Save matched gene results
# -----------------------------------

Path("results").mkdir(exist_ok=True)

signature_muscle.to_csv(
    "results/42_gene_limma_results.csv",
    index=False
)

# -----------------------------------
# 7. Save summary
# -----------------------------------

summary = pd.DataFrame({
    "metric": [
        "brain_signature_genes",
        "matched_genes",
        "mean_absolute_t",
        "median_absolute_t",
        "positive_t",
        "negative_t",
        "permutation_p"
    ],
    "value": [
        len(signature),
        len(signature_muscle),
        observed_mean_abs_t,
        np.median(np.abs(signature_t)),
        np.sum(signature_t > 0),
        np.sum(signature_t < 0),
        p_value
    ]
})

summary.to_csv(
    "results/42_gene_signature_enrichment_summary.csv",
    index=False
)

print("\nResults saved.")
print("Analysis complete.")