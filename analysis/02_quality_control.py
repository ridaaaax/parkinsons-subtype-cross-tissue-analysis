import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pathlib import Path


# -----------------------------
# 1. Load the expression data
# -----------------------------

data_file = Path("data/GSE128177_featureCounts.txt.gz")

df = pd.read_csv(
    data_file,
    sep="\t",
    compression="gzip"
)


# -----------------------------
# 2. Define sample groups
# -----------------------------

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
    "M-M-40": "Old_adult",

    "M-M-03": "Young_adult",
    "M-M-07": "Young_adult",
    "M-M-11": "Young_adult",
    "M-M-15": "Young_adult",
    "M-M-19": "Young_adult",
    "M-M-23": "Young_adult",
    "M-M-26": "Young_adult",
    "M-M-29": "Young_adult",
    "M-M-32": "Young_adult",
    "M-M-35": "Young_adult",
    "M-M-38": "Young_adult",
    "M-M-41": "Young_adult"
}


# -----------------------------
# 3. Select PD and older controls
# -----------------------------

primary_samples = [
    sample for sample in df.columns[1:]
    if groups[sample] in ["Parkinson_disease", "Old_adult"]
]

expression = df.set_index("Geneid")[primary_samples]

print("Primary samples:", len(primary_samples))
print("Expression matrix:", expression.shape)


# -----------------------------
# 4. Calculate CPM
# -----------------------------

library_sizes = expression.sum(axis=0)

cpm = expression.div(library_sizes, axis=1) * 1_000_000


# -----------------------------
# 5. Filter low-expression genes
# -----------------------------

keep_genes = (cpm > 1).sum(axis=1) >= 6

filtered_cpm = cpm.loc[keep_genes]

print("Genes before filtering:", expression.shape[0])
print("Genes after filtering:", filtered_cpm.shape[0])


# -----------------------------
# 6. Log2 transformation
# -----------------------------

log_cpm = np.log2(filtered_cpm + 1)


# -----------------------------
# 7. PCA
# -----------------------------

pca = PCA(n_components=2)

pca_coordinates = pca.fit_transform(log_cpm.T)

pca_df = pd.DataFrame(
    pca_coordinates,
    columns=["PC1", "PC2"],
    index=primary_samples
)

pca_df["group"] = [
    groups[sample]
    for sample in primary_samples
]


# -----------------------------
# 8. Print PCA information
# -----------------------------

print("\n--- PCA VARIANCE ---")

print("PC1:", round(pca.explained_variance_ratio_[0] * 100, 2), "%")
print("PC2:", round(pca.explained_variance_ratio_[1] * 100, 2), "%")


# -----------------------------
# 9. Save PCA coordinates
# -----------------------------

Path("results").mkdir(exist_ok=True)

pca_df.to_csv(
    "results/GSE128177_PCA_coordinates.csv"
)


# -----------------------------
# 10. Create PCA plot
# -----------------------------

plt.figure(figsize=(8, 6))

for group in ["Parkinson_disease", "Old_adult"]:

    subset = pca_df[pca_df["group"] == group]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group
    )

    for sample in subset.index:
        plt.annotate(
            sample,
            (subset.loc[sample, "PC1"],
             subset.loc[sample, "PC2"]),
            fontsize=7
        )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
)

plt.title("GSE128177 skeletal muscle PCA")

plt.legend()

plt.tight_layout()

plt.savefig(
    "results/GSE128177_PCA.png",
    dpi=300
)

plt.show()