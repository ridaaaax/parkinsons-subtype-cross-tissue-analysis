import gzip
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


INPUT_FILE = "data/GSE165082_PD-CC.counts.txt.gz"
RESULTS_DIR = "results"
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")


os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


# --------------------------------------------------
# 1. Read count matrix
# --------------------------------------------------

print("Reading blood count matrix...")

with gzip.open(INPUT_FILE, "rt") as f:
    df = pd.read_csv(f, sep="\t")


gene_column = df.columns[0]
genes = df[gene_column].astype(str)
counts = df.drop(columns=[gene_column])

counts.index = genes


print("Genes:", counts.shape[0])
print("Samples:", counts.shape[1])


# --------------------------------------------------
# 2. Identify groups from sample names
# --------------------------------------------------

sample_names = counts.columns.tolist()

groups = []

for sample in sample_names:
    if sample.endswith("_PD"):
        groups.append("PD")
    elif sample.endswith("_CC"):
        groups.append("Control")
    else:
        groups.append("Unknown")


group_series = pd.Series(groups, index=sample_names)

print("\nSample groups:")
print(group_series.value_counts())


# --------------------------------------------------
# 3. Library sizes
# --------------------------------------------------

library_sizes = counts.sum(axis=0)

print("\nLibrary sizes:")
print(library_sizes)


# --------------------------------------------------
# 4. CPM transformation
# --------------------------------------------------

cpm = counts.div(library_sizes, axis=1) * 1_000_000


# --------------------------------------------------
# 5. Low-expression filtering
# --------------------------------------------------

# Keep genes with CPM > 1 in at least 6 samples
keep = (cpm > 1).sum(axis=1) >= 6

filtered_cpm = cpm.loc[keep]

print("\nGenes before filtering:", counts.shape[0])
print("Genes after filtering:", filtered_cpm.shape[0])
print("Genes removed:", counts.shape[0] - filtered_cpm.shape[0])


# --------------------------------------------------
# 6. Log2-CPM transformation
# --------------------------------------------------

log_cpm = np.log2(filtered_cpm + 1)


# --------------------------------------------------
# 7. PCA
# --------------------------------------------------

# PCA requires samples as rows
X = log_cpm.T

pca = PCA(n_components=2)
coordinates = pca.fit_transform(X)

pca_df = pd.DataFrame(
    coordinates,
    columns=["PC1", "PC2"],
    index=sample_names
)

pca_df["Group"] = group_series


print("\nPCA variance explained:")
print("PC1:", round(pca.explained_variance_ratio_[0] * 100, 2), "%")
print("PC2:", round(pca.explained_variance_ratio_[1] * 100, 2), "%")


# --------------------------------------------------
# 8. Save PCA coordinates
# --------------------------------------------------

pca_output = os.path.join(
    RESULTS_DIR,
    "GSE165082_blood_PCA_coordinates.csv"
)

pca_df.to_csv(pca_output)


# --------------------------------------------------
# 9. Save QC summary
# --------------------------------------------------

qc_summary = pd.DataFrame({
    "metric": [
        "genes_before_filtering",
        "genes_after_filtering",
        "genes_removed",
        "number_of_samples",
        "PD_samples",
        "control_samples",
        "PC1_variance_percent",
        "PC2_variance_percent"
    ],
    "value": [
        counts.shape[0],
        filtered_cpm.shape[0],
        counts.shape[0] - filtered_cpm.shape[0],
        counts.shape[1],
        sum(group_series == "PD"),
        sum(group_series == "Control"),
        pca.explained_variance_ratio_[0] * 100,
        pca.explained_variance_ratio_[1] * 100
    ]
})

qc_output = os.path.join(
    RESULTS_DIR,
    "GSE165082_blood_QC_summary.csv"
)

qc_summary.to_csv(qc_output, index=False)


# --------------------------------------------------
# 10. Library-size plot
# --------------------------------------------------

plt.figure(figsize=(12, 6))

plt.bar(
    sample_names,
    library_sizes.values
)

plt.xticks(rotation=90)
plt.ylabel("Total mapped counts")
plt.title("GSE165082 blood sample library sizes")

plt.tight_layout()

plt.savefig(
    os.path.join(FIGURES_DIR, "GSE165082_blood_library_sizes.png"),
    dpi=300
)

plt.close()


# --------------------------------------------------
# 11. PCA plot
# --------------------------------------------------

plt.figure(figsize=(8, 6))

for group in ["PD", "Control"]:
    subset = pca_df[pca_df["Group"] == group]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        s=60
    )

    for sample, row in subset.iterrows():
        plt.annotate(
            sample,
            (row["PC1"], row["PC2"]),
            fontsize=7,
            xytext=(4, 4),
            textcoords="offset points"
        )


plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
)

plt.title("GSE165082 blood PCA")
plt.legend()

plt.tight_layout()

plt.savefig(
    os.path.join(FIGURES_DIR, "GSE165082_blood_PCA.png"),
    dpi=300
)

plt.close()


print("\nQC complete.")

print("Created:")
print(pca_output)
print(qc_output)
print(
    os.path.join(
        FIGURES_DIR,
        "GSE165082_blood_library_sizes.png"
    )
)
print(
    os.path.join(
        FIGURES_DIR,
        "GSE165082_blood_PCA.png"
    )
)